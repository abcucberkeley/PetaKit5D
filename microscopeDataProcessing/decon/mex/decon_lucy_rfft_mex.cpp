// Richardson-Lucy iterations of decon_lucy_omw_function ('omw') and the FFT
// steps of decon_lucy_function ('rlstep') with the same FFTs as MATLAB's own
// fftn and ifftn, so that the results equal MATLAB's bit for bit:
//   fftn of real data: 1D real-to-complex transforms along y, complex ones
//     along x and z (MATLAB's single precision FFTW and its plans), then the
//     Hermitian fill of MATLAB (an element whose mirror has the lower linear
//     index is the conjugate of the mirror, self-conjugate elements are real)
//   ifftn of the (Hermitian) product: complex transforms along y and x, then
//     complex-to-real along z, each divided by its length
//   MATLAB runs the 1D transforms of a stage in groups of four lines and some
//     lines one by one (see mark_alone); for some lengths (e.g., 64 or 128)
//     the two give different last bits, so the mex runs the same lines alone.
// It is faster than MATLAB and needs less memory because only the half spectrum
// (k_z <= nz / 2) that the inverse uses is kept and transformed back, the
// element-wise steps run inside the transforms' copy loops, and no complex
// volume of the full size is allocated. 2D data uses the same steps along y
// and x.
//
// H = decon_lucy_rfft_mex('otf', psf, outSize, nThreads, grain)
//   psf      : single, normalized PSF, at most outSize in each dimension
//   outSize  : data size [ny, nx] or [ny, nx, nz] (ny, nx >= 2)
//   H        : complex single half OTF, the first floor(n / 2) + 1 elements
//              along the last dimension of decon_psf2otf(psf, outSize)
//
// [J, t] = decon_lucy_rfft_mex('omw', I, H_f, H_b, nIter, background, nThreads, grain)
//   I          : uint16 or single raw data; the iterations use max(I - background, 0)
//   H_f, H_b   : half OTFs from 'otf'
//   J          : single result after nIter iterations
//
// [J, t] = decon_lucy_rfft_mex('rlstep', Y, I, H, nThreads, grain)
//   one update of decon_lucy_function: max(Y .* real(ifftn(conj(H) .* fftn(I ./ max(real(ifftn(H .* fftn(Y))), eps)))), 0)
//   Y, I       : single (I background subtracted)
//
// Y = decon_lucy_rfft_mex('fft1', X, nThreads, grain)
//   rows 1 .. n / 2 + 1 of fft(X, [], 1) for single X (n x m, n >= 2), with MATLAB's
//   groups of lines for the given grain (to find out which grain MATLAB uses)
//
//   nThreads   : threads
//   grain      : most lines per piece of MATLAB's split of a stage among its threads: 32 when
//                MATLAB runs its FFTs in parallel, 0 when it doesn't (decon_lucy_rfft finds it out)
//   t          : seconds in [planning, forward transforms, inverse transforms, other]
//
// Build (-R2018a for interleaved complex arrays; no fused multiply-adds, as MATLAB computes products):
// mex -R2018a CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp -ffp-contract=off' LDFLAGS='$LDFLAGS -O3 -fopenmp' decon_lucy_rfft_mex.cpp -L<matlabroot>/bin/glnxa64 -lmwfftw3f
// (compileAllMexFiles builds it on all platforms)

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include "mex.h"
#include "fftw3.h"

typedef ptrdiff_t idx_t;

// transforms per FFTW execute, as in MATLAB's fft
static const int kB = 4;

// byte offsets (modulo 64) of the scratch buffers, as MATLAB's
static const int kInOff = 16;
static const int kOutOff = 32;

struct Timer {
    double t[4];
    double t0;
    Timer() : t0(0) {
        for (int i = 0; i < 4; i++) t[i] = 0;
    }
    void start() { t0 = omp_get_wtime(); }
    void stop(const int k) { t[k] += omp_get_wtime() - t0; }
};

// sizes, plans, the lines MATLAB transforms alone, and per-thread scratch of one call
struct Ctx {
    int ny, nx, nz;      // data size (nz = 1 for 2D)
    bool is3d;
    int nyh;             // ny / 2 + 1, complex elements of the y-half of a column
    int nh;              // kept elements along the last dimension: nz / 2 + 1 (3D), nx / 2 + 1 (2D)
    idx_t ldr;           // floats per column of the real / y-half buffer (2 * nyh)
    idx_t nP;            // complex elements of a half spectrum
    int nt;
    // plans [0]: groups of four, [1]: single transforms
    fftwf_plan fy[2], fx[2], fz[2], iy[2], ix[2], il[2];
    // nonzero for the lines MATLAB transforms alone: along y (column x + nx * z), along x (row y + ny * z),
    // along z (y + ny * x, 3D only)
    uint8_t* aloneY;
    uint8_t* aloneX;
    uint8_t* aloneZ;
    // scratch per thread: in, out (bases for fftwf_free) and the x-line band buffer of the 3D forward
    void** base;
    fftwf_complex** in;
    fftwf_complex** out;
    fftwf_complex** band;
    Ctx() : nt(0), aloneY(NULL), aloneX(NULL), aloneZ(NULL), base(NULL), in(NULL), out(NULL), band(NULL) {
        for (int i = 0; i < 2; i++) fy[i] = fx[i] = fz[i] = iy[i] = ix[i] = il[i] = NULL;
    }
    ~Ctx() {
        release();
    }
    void release() {
        fftwf_plan* ps[6] = {fy, fx, fz, iy, ix, il};
        for (int k = 0; k < 6; k++) {
            for (int i = 0; i < 2; i++) {
                if (ps[k][i] != NULL) fftwf_destroy_plan(ps[k][i]);
                ps[k][i] = NULL;
            }
        }
        uint8_t** as[3] = {&aloneY, &aloneX, &aloneZ};
        for (int k = 0; k < 3; k++) {
            if (*as[k] != NULL) mxFree(*as[k]);
            *as[k] = NULL;
        }
        if (base != NULL) {
            for (int i = 0; i < 3 * nt; i++) {
                if (base[i] != NULL) fftwf_free(base[i]);
            }
            mxFree(base);
            base = NULL;
        }
    }
};

// Marks the lines a .. b - 1 of a stage that MATLAB transforms alone. MATLAB's FFT splits the lines of a
// stage among its threads by halving them (at a + (b - a) / 2) until a piece has at most grain lines (when
// it runs in parallel and the stage is large enough, see stage_grain; grain 0: one piece), and runs a piece
// in groups of four lines from its start, the last (length % 4) lines one by one. Lines along x of 3D data
// (seg = ny) also restart the groups at each z slice.
static void mark_alone(uint8_t* alone, const idx_t a, const idx_t b, const idx_t grain, const idx_t seg) {
    if (grain > 0 && b - a > grain) {
        const idx_t mid = a + (b - a) / 2;
        mark_alone(alone, a, mid, grain, seg);
        mark_alone(alone, mid, b, grain, seg);
        return;
    }
    idx_t s = a;
    while (s < b) {
        idx_t e = b;
        if (seg > 0 && (s / seg + 1) * seg < b) {
            e = (s / seg + 1) * seg;
        }
        for (idx_t L = e - (e - s) % kB; L < e; L++) {
            alone[L] = 1;
        }
        s = e;
    }
}

// MATLAB splits the N lines of length n of a stage only if n * N * (bits of n) > 98304 (its estimate of the
// work of the stage), otherwise one thread transforms them
static idx_t stage_grain(const idx_t grain, const int n, const idx_t N) {
    int bits = 0;
    for (int v = n; v > 0; v >>= 1) bits++;
    return (grain > 0 && (double)n * (double)N * bits > 98304.0) ? grain : 0;
}

static void* alloc_off(const size_t bytes, const int off, void** base) {
    *base = fftwf_malloc(bytes + 128);
    if (*base == NULL) return NULL;
    uintptr_t p = ((uintptr_t)(*base) + 63) & ~(uintptr_t)63;
    return (void*)(p + off);
}

static void setup(Ctx& c, const int* d, const int nThreads, const idx_t grain) {
    c.ny = d[0];
    c.nx = d[1];
    c.nz = d[2];
    c.is3d = d[2] > 1;
    c.nyh = c.ny / 2 + 1;
    c.nh = c.is3d ? c.nz / 2 + 1 : c.nx / 2 + 1;
    c.ldr = 2 * (idx_t)c.nyh;
    c.nP = c.is3d ? (idx_t)c.ny * c.nx * c.nh : (idx_t)c.ny * c.nh;
    c.nt = nThreads;
    // the lines MATLAB transforms alone, per stage (inverse stages as the forward ones of the full spectrum)
    const idx_t nlY = (idx_t)c.nx * c.nz, nlX = (idx_t)c.ny * c.nz;
    c.aloneY = (uint8_t*)mxCalloc(nlY, 1);
    c.aloneX = (uint8_t*)mxCalloc(nlX, 1);
    mark_alone(c.aloneY, 0, nlY, stage_grain(grain, c.ny, nlY), 0);
    mark_alone(c.aloneX, 0, nlX, stage_grain(grain, c.nx, nlX), c.is3d ? c.ny : 0);
    if (c.is3d) {
        const idx_t nlZ = (idx_t)c.ny * c.nx;
        c.aloneZ = (uint8_t*)mxCalloc(nlZ, 1);
        mark_alone(c.aloneZ, 0, nlZ, stage_grain(grain, c.nz, nlZ), 0);
    }
    const int maxn = c.ny > c.nx ? (c.ny > c.nz ? c.ny : c.nz) : (c.nx > c.nz ? c.nx : c.nz);
    const size_t sbytes = sizeof(fftwf_complex) * kB * (size_t)maxn;
    c.base = (void**)mxCalloc(3 * nThreads, sizeof(void*));
    c.in = (fftwf_complex**)mxCalloc(nThreads, sizeof(fftwf_complex*));
    c.out = (fftwf_complex**)mxCalloc(nThreads, sizeof(fftwf_complex*));
    c.band = (fftwf_complex**)mxCalloc(nThreads, sizeof(fftwf_complex*));
    for (int t = 0; t < nThreads; t++) {
        c.in[t] = (fftwf_complex*)alloc_off(sbytes, kInOff, &c.base[3 * t]);
        c.out[t] = (fftwf_complex*)alloc_off(sbytes, kOutOff, &c.base[3 * t + 1]);
        if (c.is3d) {
            c.band[t] = (fftwf_complex*)alloc_off(sizeof(fftwf_complex) * kB * (size_t)c.nx * c.nz, 0, &c.base[3 * t + 2]);
        }
        if (c.in[t] == NULL || c.out[t] == NULL || (c.is3d && c.band[t] == NULL)) {
            c.release();
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:outOfMemory", "Out of memory for the decon buffers.");
        }
    }
    // 1D plans as MATLAB makes them: rank 1, contiguous transforms, FFTW_ESTIMATE
    float* ri = (float*)c.in[0];
    float* ro = (float*)c.out[0];
    fftwf_complex* ci = c.in[0];
    fftwf_complex* co = c.out[0];
    for (int b = 0; b < 2; b++) {
        const ptrdiff_t m = b == 0 ? kB : 1;
        fftwf_iodim64 dy = {c.ny, 1, 1}, hy = {m, c.ny, c.nyh}, hyc = {m, c.ny, c.ny};
        c.fy[b] = fftwf_plan_guru64_dft_r2c(1, &dy, 1, &hy, ri, co, FFTW_ESTIMATE);
        c.iy[b] = fftwf_plan_guru64_dft(1, &dy, 1, &hyc, ci, co, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftwf_iodim64 dx = {c.nx, 1, 1}, hx = {m, c.nx, c.nx};
        c.fx[b] = fftwf_plan_guru64_dft(1, &dx, 1, &hx, ci, co, FFTW_FORWARD, FFTW_ESTIMATE);
        if (c.is3d) {
            c.ix[b] = fftwf_plan_guru64_dft(1, &dx, 1, &hx, ci, co, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftwf_iodim64 dz = {c.nz, 1, 1}, hz = {m, c.nz, c.nz}, hzr = {m, c.nh, c.nz};
            c.fz[b] = fftwf_plan_guru64_dft(1, &dz, 1, &hz, ci, co, FFTW_FORWARD, FFTW_ESTIMATE);
            c.il[b] = fftwf_plan_guru64_dft_c2r(1, &dz, 1, &hzr, ci, ro, FFTW_ESTIMATE);
        } else {
            fftwf_iodim64 hxr = {m, c.nh, c.nx};
            c.il[b] = fftwf_plan_guru64_dft_c2r(1, &dx, 1, &hxr, ci, ro, FFTW_ESTIMATE);
        }
        if (c.fy[b] == NULL || c.iy[b] == NULL || c.fx[b] == NULL || c.il[b] == NULL || (c.is3d && (c.ix[b] == NULL || c.fz[b] == NULL))) {
            c.release();
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:fftwPlan", "FFTW planning failed.");
        }
    }
}

// lanes of a block of m lines (first line j0) that MATLAB transforms in a group: their offsets in gl, the
// count returned
static inline int grouped_lanes(const uint8_t* alone, const idx_t j0, const int m, int* gl) {
    int ng = 0;
    for (int l = 0; l < m; l++) {
        if (!alone[j0 + l]) gl[ng++] = l;
    }
    return ng;
}

// complex element of the full y range of a y-half column (rows above ny / 2 are conjugates)
static inline void yfull(const float* col, const int ny, const int nyh, const int k, float& re, float& im) {
    if (k < nyh) {
        re = col[2 * k];
        im = col[2 * k + 1];
    } else {
        re = col[2 * (ny - k)];
        im = -col[2 * (ny - k) + 1];
    }
}

// line along x (row y, slice z) of the full y range from the y-half columns of A
static inline void x_line(const Ctx& c, const float* A, const int y, const int z, float* li) {
    for (int x = 0; x < c.nx; x++) {
        yfull(A + ((idx_t)x + (idx_t)c.nx * z) * c.ldr, c.ny, c.nyh, y, li[2 * x], li[2 * x + 1]);
    }
}

// forward stage 1: real columns (stride lds floats; nx * nz of them) to y-half columns of A (stride ldr
// floats). A may be src (lds == ldr).
static void fwd_y(Ctx& c, const float* src, const idx_t lds, float* A) {
    const idx_t ncol = (idx_t)c.nx * c.nz;
    const idx_t nb = (ncol + kB - 1) / kB;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const idx_t j0 = b * kB;
        const int m = (int)(ncol - j0 < kB ? ncol - j0 : kB);
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        int gl[kB];
        const int ng = grouped_lanes(c.aloneY, j0, m, gl);
        if (ng > 0) {
            // all inputs are read before A is written (A may be src)
            for (int q = 0; q < kB; q++) {
                memcpy(in + (idx_t)q * c.ny, src + (j0 + gl[q < ng ? q : ng - 1]) * lds, sizeof(float) * c.ny);
            }
            fftwf_execute_dft_r2c(c.fy[0], in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                memcpy(A + (j0 + gl[q]) * c.ldr, out + (idx_t)q * c.ldr, sizeof(float) * c.ldr);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneY[j0 + l]) continue;
            memcpy(in, src + (j0 + l) * lds, sizeof(float) * c.ny);
            fftwf_execute_dft_r2c(c.fy[1], in, (fftwf_complex*)out);
            memcpy(A + (j0 + l) * c.ldr, out, sizeof(float) * c.ldr);
        }
    }
}

// x transform of band row l at slice z into the band buffer [m][nx][nz] (z contiguous)
static inline void to_band(float* band, const float* lo, const int l, const int z, const int nx, const int nz) {
    float* bl = band + 2 * ((idx_t)l * nx) * nz;
    for (int k = 0; k < nx; k++) {
        bl[2 * ((idx_t)k * nz + z)] = lo[2 * k];
        bl[2 * ((idx_t)k * nz + z) + 1] = lo[2 * k + 1];
    }
}

// the kept outputs k_z < nh of the z line of row y, column k to P [ny, nx, nh]
static inline void z_out(const Ctx& c, fftwf_complex* P, const float* lo, const int y, const int k) {
    const idx_t nl = (idx_t)c.ny * c.nx;
    float* dst = (float*)P + 2 * ((idx_t)y + (idx_t)c.ny * k);
    for (int kz = 0; kz < c.nh; kz++) {
        dst[2 * kz * nl] = lo[2 * kz];
        dst[2 * kz * nl + 1] = lo[2 * kz + 1];
    }
}

// a group of np z lines (band rows pl, columns pk), padded with its last line
static void z_group(const Ctx& c, const float* band, const int* pl, const int* pk, const int np, float* in, float* out,
        fftwf_complex* P, const int y0) {
    const int nz = c.nz, nx = c.nx;
    for (int q = 0; q < kB; q++) {
        const int s = q < np ? q : np - 1;
        memcpy(in + 2 * (idx_t)q * nz, band + 2 * ((idx_t)pl[s] * nx + pk[s]) * nz, sizeof(fftwf_complex) * nz);
    }
    fftwf_execute_dft(c.fz[0], (fftwf_complex*)in, (fftwf_complex*)out);
    for (int q = 0; q < np; q++) {
        z_out(c, P, out + 2 * (idx_t)q * nz, y0 + pl[q], pk[q]);
    }
}

// forward stages 2 and 3 (3D): lines along x, then along z; the z outputs k_z < nh go to P [ny, nx, nh].
// Rows are processed in bands of four; lines MATLAB transforms in groups are grouped in any fours (the lanes
// of a group are independent), the others transformed alone.
static void fwd_xz(Ctx& c, const float* A, fftwf_complex* P) {
    const int ny = c.ny, nx = c.nx, nz = c.nz;
    const int nbands = (ny + kB - 1) / kB;
    #pragma omp parallel for num_threads(c.nt) schedule(dynamic, 1)
    for (int b = 0; b < nbands; b++) {
        const int t = omp_get_thread_num();
        const int y0 = b * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        float* band = (float*)c.band[t];
        // lines along x of the band rows, every z slice
        for (int z = 0; z < nz; z++) {
            const idx_t L0 = (idx_t)y0 + (idx_t)ny * z;
            int gl[kB];
            const int ng = grouped_lanes(c.aloneX, L0, m, gl);
            if (ng > 0) {
                for (int q = 0; q < kB; q++) {
                    x_line(c, A, y0 + gl[q < ng ? q : ng - 1], z, in + 2 * (idx_t)q * nx);
                }
                fftwf_execute_dft(c.fx[0], (fftwf_complex*)in, (fftwf_complex*)out);
                for (int q = 0; q < ng; q++) {
                    to_band(band, out + 2 * (idx_t)q * nx, gl[q], z, nx, nz);
                }
            }
            for (int l = 0; l < m; l++) {
                if (!c.aloneX[L0 + l]) continue;
                x_line(c, A, y0 + l, z, in);
                fftwf_execute_dft(c.fx[1], (fftwf_complex*)in, (fftwf_complex*)out);
                to_band(band, out, l, z, nx, nz);
            }
        }
        // lines along z (row y0 + l, column k)
        int pl[kB], pk[kB], np = 0;
        for (int k = 0; k < nx; k++) {
            for (int l = 0; l < m; l++) {
                if (c.aloneZ[(idx_t)(y0 + l) + (idx_t)ny * k]) {
                    memcpy(in, band + 2 * ((idx_t)l * nx + k) * nz, sizeof(fftwf_complex) * nz);
                    fftwf_execute_dft(c.fz[1], (fftwf_complex*)in, (fftwf_complex*)out);
                    z_out(c, P, out, y0 + l, k);
                    continue;
                }
                pl[np] = l;
                pk[np] = k;
                np++;
                if (np == kB) {
                    z_group(c, band, pl, pk, np, in, out, P, y0);
                    np = 0;
                }
            }
        }
        if (np > 0) {
            z_group(c, band, pl, pk, np, in, out, P, y0);
        }
    }
}

// the kept outputs k_x < nh of the x transform of row y (2D) to P [ny, nh]
static inline void x2d_out(const Ctx& c, fftwf_complex* P, const float* lo, const int y) {
    float* dst = (float*)P + 2 * (idx_t)y;
    for (int k = 0; k < c.nh; k++) {
        dst[2 * (idx_t)k * c.ny] = lo[2 * k];
        dst[2 * (idx_t)k * c.ny + 1] = lo[2 * k + 1];
    }
}

// forward stage 2 (2D): lines along x (rows)
static void fwd_x2d(Ctx& c, const float* A, fftwf_complex* P) {
    const int ny = c.ny, nx = c.nx;
    const int nb = (ny + kB - 1) / kB;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (int b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const int y0 = b * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        int gl[kB];
        const int ng = grouped_lanes(c.aloneX, y0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                x_line(c, A, y0 + gl[q < ng ? q : ng - 1], 0, in + 2 * (idx_t)q * nx);
            }
            fftwf_execute_dft(c.fx[0], (fftwf_complex*)in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                x2d_out(c, P, out + 2 * (idx_t)q * nx, y0 + gl[q]);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneX[y0 + l]) continue;
            x_line(c, A, y0 + l, 0, in);
            fftwf_execute_dft(c.fx[1], (fftwf_complex*)in, (fftwf_complex*)out);
            x2d_out(c, P, out, y0 + l);
        }
    }
}

// MATLAB's Hermitian fill on the kept half: in the self-conjugate planes (k_z = 0 and nz / 2 for 3D, k_x = 0
// and nx / 2 for 2D) an element whose mirror has the lower linear index becomes the mirror's conjugate, a
// self-conjugate element its real part
static void hermitian_fill(Ctx& c, fftwf_complex* P) {
    const int ny = c.ny, nx = c.nx;
    float* p = (float*)P;
    if (c.is3d) {
        const int planes[2] = {0, c.nz / 2};
        const int np = (c.nz % 2 == 0) ? 2 : 1;
        for (int q = 0; q < np; q++) {
            float* pl = p + 2 * (idx_t)planes[q] * ny * nx;
            #pragma omp parallel for num_threads(c.nt) schedule(static)
            for (int k2 = 0; k2 < nx; k2++) {
                const int m2 = (nx - k2) % nx;
                for (int k1 = 0; k1 < ny; k1++) {
                    const int m1 = (ny - k1) % ny;
                    const idx_t L = k1 + (idx_t)ny * k2, Lm = m1 + (idx_t)ny * m2;
                    if (Lm < L) {
                        pl[2 * L] = pl[2 * Lm];
                        pl[2 * L + 1] = -pl[2 * Lm + 1];
                    } else if (Lm == L) {
                        pl[2 * L + 1] = 0.0f;
                    }
                }
            }
        }
    } else {
        const int cols[2] = {0, nx / 2};
        const int ncl = (nx % 2 == 0) ? 2 : 1;
        for (int q = 0; q < ncl; q++) {
            float* cl = p + 2 * (idx_t)cols[q] * ny;
            for (int k1 = 0; k1 < ny; k1++) {
                const int m1 = (ny - k1) % ny;
                if (m1 < k1) {
                    cl[2 * k1] = cl[2 * m1];
                    cl[2 * k1 + 1] = -cl[2 * m1 + 1];
                } else if (m1 == k1) {
                    cl[2 * k1 + 1] = 0.0f;
                }
            }
        }
    }
}

// column j of P .* H (or conj(H) .* P) into li (MATLAB's complex product: no fused multiply-add)
static inline void product_column(const fftwf_complex* P, const mxComplexSingle* H, const idx_t j, const int ny,
        const bool conjugate, float* li) {
    const float* pc = (const float*)P + 2 * j * ny;
    const float* hc = (const float*)H + 2 * j * ny;
    for (int k = 0; k < ny; k++) {
        const float a = pc[2 * k], b = pc[2 * k + 1], cr = hc[2 * k], ci = hc[2 * k + 1];
        if (conjugate) {
            li[2 * k] = a * cr + b * ci;
            li[2 * k + 1] = b * cr - a * ci;
        } else {
            li[2 * k] = a * cr - b * ci;
            li[2 * k + 1] = a * ci + b * cr;
        }
    }
}

// inverse stage 1: P .* H (or conj(H) .* P) along y columns, c2c, divided by ny, in place. Only the kept
// columns (the first ones of MATLAB's nx * nz) are transformed.
static void inv_y(Ctx& c, fftwf_complex* P, const mxComplexSingle* H, const bool conjugate) {
    const int ny = c.ny;
    const idx_t ncol = c.is3d ? (idx_t)c.nx * c.nh : (idx_t)c.nh;
    const idx_t nb = (ncol + kB - 1) / kB;
    const float fn = (float)ny;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        const idx_t j0 = b * kB;
        const int m = (int)(ncol - j0 < kB ? ncol - j0 : kB);
        int gl[kB];
        const int ng = grouped_lanes(c.aloneY, j0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                product_column(P, H, j0 + gl[q < ng ? q : ng - 1], ny, conjugate, in + 2 * (idx_t)q * ny);
            }
            fftwf_execute_dft(c.iy[0], (fftwf_complex*)in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                float* pc = (float*)P + 2 * (j0 + gl[q]) * ny;
                const float* lo = out + 2 * (idx_t)q * ny;
                for (int k = 0; k < 2 * ny; k++) pc[k] = lo[k] / fn;
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneY[j0 + l]) continue;
            product_column(P, H, j0 + l, ny, conjugate, in);
            fftwf_execute_dft(c.iy[1], (fftwf_complex*)in, (fftwf_complex*)out);
            float* pc = (float*)P + 2 * (j0 + l) * ny;
            for (int k = 0; k < 2 * ny; k++) pc[k] = out[k] / fn;
        }
    }
}

// inverse stage 2 (3D): lines along x of the kept z planes, c2c, divided by nx
static void inv_x(Ctx& c, fftwf_complex* P) {
    const int ny = c.ny, nx = c.nx;
    const int nb = (ny + kB - 1) / kB;
    const idx_t nl = (idx_t)ny * nx;
    const float fn = (float)nx;
    const idx_t ntask = (idx_t)c.nh * nb;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t task = 0; task < ntask; task++) {
        const int t = omp_get_thread_num();
        const int kz = (int)(task / nb);
        const int y0 = (int)(task % nb) * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        float* pl = (float*)P + 2 * (idx_t)kz * nl;
        const idx_t L0 = (idx_t)y0 + (idx_t)ny * kz;
        int gl[kB];
        const int ng = grouped_lanes(c.aloneX, L0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                const int y = y0 + gl[q < ng ? q : ng - 1];
                float* li = in + 2 * (idx_t)q * nx;
                for (int x = 0; x < nx; x++) {
                    li[2 * x] = pl[2 * ((idx_t)y + (idx_t)ny * x)];
                    li[2 * x + 1] = pl[2 * ((idx_t)y + (idx_t)ny * x) + 1];
                }
            }
            fftwf_execute_dft(c.ix[0], (fftwf_complex*)in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                const int y = y0 + gl[q];
                const float* lo = out + 2 * (idx_t)q * nx;
                for (int x = 0; x < nx; x++) {
                    pl[2 * ((idx_t)y + (idx_t)ny * x)] = lo[2 * x] / fn;
                    pl[2 * ((idx_t)y + (idx_t)ny * x) + 1] = lo[2 * x + 1] / fn;
                }
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneX[L0 + l]) continue;
            const int y = y0 + l;
            for (int x = 0; x < nx; x++) {
                in[2 * x] = pl[2 * ((idx_t)y + (idx_t)ny * x)];
                in[2 * x + 1] = pl[2 * ((idx_t)y + (idx_t)ny * x) + 1];
            }
            fftwf_execute_dft(c.ix[1], (fftwf_complex*)in, (fftwf_complex*)out);
            for (int x = 0; x < nx; x++) {
                pl[2 * ((idx_t)y + (idx_t)ny * x)] = out[2 * x] / fn;
                pl[2 * ((idx_t)y + (idx_t)ny * x) + 1] = out[2 * x + 1] / fn;
            }
        }
    }
}

// element-wise step on the real result v of the inverse at voxel i (column-major index of the data)
template <typename T>
struct Ratio {
    // dst = I ./ max(v, eps), dst with column stride ldd
    float* dst;
    idx_t ldd;
    const T* raw;
    double bg;
    int ny;
    inline void operator()(const idx_t i, const idx_t col, const int y, const float v) const {
        const float epsf = (float)DBL_EPSILON;
        const float cx = v > epsf ? v : epsf;
        float in = (float)((double)raw[i] - bg);
        in = in > 0.0f ? in : 0.0f;
        dst[col * ldd + y] = in / cx;
    }
};

struct Update {
    // out = max(v .* J, 0)
    float* out;
    const float* j;
    inline void operator()(const idx_t i, const idx_t, const int, const float v) const {
        const float w = v * j[i];
        out[i] = w > 0.0f ? w : 0.0f;
    }
};

// the kept half spectrum of line L of the last dimension (k = 0 .. nh - 1, stride nl) into li
static inline void last_line(const fftwf_complex* P, const idx_t L, const idx_t nl, const int nh, float* li) {
    const float* src = (const float*)P + 2 * L;
    for (int k = 0; k < nh; k++) {
        li[2 * k] = src[2 * k * nl];
        li[2 * k + 1] = src[2 * k * nl + 1];
    }
}

// the real line L (lo, nlen values) divided by its length, then the element-wise step
template <typename Op>
static inline void last_out(const Ctx& c, const Op& op, const float* lo, const idx_t L, const int nlen, const idx_t nl,
        const float fn) {
    const int y = (int)(L % c.ny);
    const idx_t col0 = L / c.ny;
    const idx_t colstep = c.is3d ? (idx_t)c.nx : 1;
    for (int k = 0; k < nlen; k++) {
        op(L + (idx_t)k * nl, col0 + (idx_t)k * colstep, y, lo[k] / fn);
    }
}

// inverse last stage: c2r along z (2D: x) of the kept half, divided by its length, then the element-wise step
template <typename Op>
static void inv_last(Ctx& c, const fftwf_complex* P, const Op& op) {
    const int nh = c.nh;
    const int nlen = c.is3d ? c.nz : c.nx;                       // length of the real lines
    const idx_t nl = c.is3d ? (idx_t)c.ny * c.nx : (idx_t)c.ny;  // number of lines
    const uint8_t* alone = c.is3d ? c.aloneZ : c.aloneX;
    const idx_t nb = (nl + kB - 1) / kB;
    const float fn = (float)nlen;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const idx_t L0 = b * kB;
        const int m = (int)(nl - L0 < kB ? nl - L0 : kB);
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        int gl[kB];
        const int ng = grouped_lanes(alone, L0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                last_line(P, L0 + gl[q < ng ? q : ng - 1], nl, nh, in + 2 * (idx_t)q * nh);
            }
            fftwf_execute_dft_c2r(c.il[0], (fftwf_complex*)in, out);
            for (int q = 0; q < ng; q++) {
                last_out(c, op, out + (idx_t)q * nlen, L0 + gl[q], nlen, nl, fn);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!alone[L0 + l]) continue;
            last_line(P, L0 + l, nl, nh, in);
            fftwf_execute_dft_c2r(c.il[1], (fftwf_complex*)in, out);
            last_out(c, op, out, L0 + l, nlen, nl, fn);
        }
    }
}

// fftn of a real volume (columns of stride lds floats) into the kept half P, with MATLAB's Hermitian fill;
// A: y-half buffer (may be src when lds == ldr)
static void fftn_half(Ctx& c, const float* src, const idx_t lds, float* A, fftwf_complex* P, Timer& tm) {
    tm.start();
    fwd_y(c, src, lds, A);
    if (c.is3d) {
        fwd_xz(c, A, P);
    } else {
        fwd_x2d(c, A, P);
    }
    hermitian_fill(c, P);
    tm.stop(1);
}

template <typename Op>
static void ifftn_op(Ctx& c, fftwf_complex* P, const mxComplexSingle* H, const bool conjugate, const Op& op, Timer& tm) {
    tm.start();
    inv_y(c, P, H, conjugate);
    if (c.is3d) {
        inv_x(c, P);
    }
    inv_last(c, P, op);
    tm.stop(2);
}

// buffers of one call (fftwf_malloc), freed at the end and on errors (also MATLAB's)
struct Buffers {
    void* b[4];
    Buffers() {
        for (int i = 0; i < 4; i++) b[i] = NULL;
    }
    ~Buffers() {
        for (int i = 0; i < 4; i++) {
            if (b[i] != NULL) fftwf_free(b[i]);
            b[i] = NULL;
        }
    }
    void* get(const int slot, const size_t bytes) {
        b[slot] = fftwf_malloc(bytes);
        if (b[slot] == NULL) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:outOfMemory", "Out of memory for the decon buffers.");
        }
        return b[slot];
    }
};

template <typename T>
static inline float input_value(const T* raw, const idx_t i, const double bg) {
    const float v = (float)((double)raw[i] - bg);
    return v > 0.0f ? v : 0.0f;
}

// The ratio buffer R (y-half columns, at least n floats) also takes the result: it is free during the last
// inverse transform. It comes from mxMalloc, so that the output array can take it over (result_array).
static float* ratio_buffer(const Ctx& c) {
    return (float*)mxMalloc(sizeof(float) * (size_t)c.ldr * c.nx * c.nz);
}

// single array of the data size with the first n floats of the buffer r (from mxMalloc) as its data
static mxArray* result_array(float* r, const idx_t n, const mwSize nd, const mwSize* dims) {
    r = (float*)mxRealloc(r, sizeof(float) * (size_t)n);
    mxArray* a = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);
    mxSetSingles(a, r);
    mxSetDimensions(a, dims, nd);
    return a;
}

// OMW: J <- max(real(ifftn(fftn(I ./ max(real(ifftn(fftn(J) .* H_f)), eps)) .* H_b)) .* J, 0); returns the
// result in a ratio buffer
template <typename T>
static float* omw_iterations(const T* raw, const double bg, const mxComplexSingle* hf, const mxComplexSingle* hb,
        const int nIter, Ctx& c, Timer& tm) {
    const idx_t n = (idx_t)c.ny * c.nx * c.nz;
    float* r = ratio_buffer(c);
    if (nIter == 0) {
        #pragma omp parallel for num_threads(c.nt) schedule(static)
        for (idx_t i = 0; i < n; i++) r[i] = input_value(raw, i, bg);
        return r;
    }
    Buffers buf;
    float* j = (float*)buf.get(0, sizeof(float) * n);
    fftwf_complex* p = (fftwf_complex*)buf.get(2, sizeof(fftwf_complex) * (size_t)c.nP);
    tm.start();
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t i = 0; i < n; i++) j[i] = input_value(raw, i, bg);
    tm.stop(3);
    Ratio<T> ratio = {r, c.ldr, raw, bg, c.ny};
    for (int k = 1; k <= nIter; k++) {
        fftn_half(c, j, c.ny, r, p, tm);
        ifftn_op(c, p, hf, false, ratio, tm);
        fftn_half(c, r, c.ldr, r, p, tm);
        Update up = {k == nIter ? r : j, j};
        ifftn_op(c, p, hb, false, up, tm);
    }
    return r;
}

// one simplified RL update: max(Y .* real(ifftn(conj(H) .* fftn(I ./ max(real(ifftn(H .* fftn(Y))), eps)))), 0);
// returns the result in a ratio buffer
static float* rl_step(const float* y, const float* im, const mxComplexSingle* h, Ctx& c, Timer& tm) {
    float* r = ratio_buffer(c);
    Buffers buf;
    fftwf_complex* p = (fftwf_complex*)buf.get(2, sizeof(fftwf_complex) * (size_t)c.nP);
    Ratio<float> ratio = {r, c.ldr, im, 0.0, c.ny};
    fftn_half(c, y, c.ny, r, p, tm);
    ifftn_op(c, p, h, false, ratio, tm);
    fftn_half(c, r, c.ldr, r, p, tm);
    Update up = {r, y};
    ifftn_op(c, p, h, true, up, tm);
    return r;
}

// half OTF: fftn (with MATLAB's Hermitian fill) of the PSF zero-padded to the data size and circularly shifted
// by -floor(psfSize / 2), kept half only
static void otf_half(mxComplexSingle* hout, const float* psf, const int* ps, Ctx& c, Timer& tm) {
    Buffers buf;
    float* r = (float*)buf.get(1, sizeof(float) * (size_t)c.ldr * c.nx * c.nz);
    tm.start();
    const idx_t nr = (idx_t)c.ldr * c.nx * c.nz;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t i = 0; i < nr; i++) r[i] = 0.0f;
    const int cy = ps[0] / 2, cx = ps[1] / 2, cz = ps[2] / 2;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (int z = 0; z < ps[2]; z++) {
        const idx_t oz = (idx_t)((z - cz + c.nz) % c.nz) * c.nx;
        for (int x = 0; x < ps[1]; x++) {
            const idx_t col = oz + (x - cx + c.nx) % c.nx;
            const float* pp = psf + ((idx_t)z * ps[1] + x) * ps[0];
            for (int yy = 0; yy < ps[0]; yy++) {
                r[col * c.ldr + (yy - cy + c.ny) % c.ny] = pp[yy];
            }
        }
    }
    tm.stop(3);
    fftn_half(c, r, c.ldr, r, (fftwf_complex*)hout, tm);
}

// data size as 3 int dimensions (2D: nz = 1)
static void get_size(const mxArray* a, int* d) {
    const mwSize nd = mxGetNumberOfDimensions(a);
    const mwSize* dims = mxGetDimensions(a);
    if (nd > 3) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Only 2D and 3D data are supported.");
    }
    for (int i = 0; i < 3; i++) {
        const mwSize v = i < (int)nd ? dims[i] : 1;
        if (v > 2147483647) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Dimensions must be smaller than 2^31.");
        }
        d[i] = (int)v;
    }
}

static void check_shape(const int* d) {
    if (d[0] < 2 || d[1] < 2) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "The data must have at least 2 elements along y and x.");
    }
}

// the half OTF size for data size d
static void half_size(const int* d, mwSize* dh, mwSize& ndh) {
    if (d[2] > 1) {
        dh[0] = d[0];
        dh[1] = d[1];
        dh[2] = d[2] / 2 + 1;
        ndh = 3;
    } else {
        dh[0] = d[0];
        dh[1] = d[1] / 2 + 1;
        dh[2] = 1;
        ndh = 2;
    }
}

static void check_otf(const mxArray* h, const int* d) {
    if (!mxIsSingle(h) || !mxIsComplex(h)) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "The OTFs must be complex single half OTFs.");
    }
    int dh[3];
    get_size(h, dh);
    mwSize e[3], ne;
    half_size(d, e, ne);
    if ((mwSize)dh[0] != e[0] || (mwSize)dh[1] != e[1] || (mwSize)dh[2] != e[2]) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "The OTF size must be that of the half OTF of the data size.");
    }
}

// count from a numeric scalar as a loop 1 : v runs it (0 for v < 1 or NaN)
static int scalar_int(const mxArray* a, const char* name) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "%s must be a real numeric scalar.", name);
    }
    const double v = mxGetScalar(a);
    if (!(v >= 1)) {
        return 0;
    }
    if (v > 2147483647.0) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "%s is too large.", name);
    }
    return (int)v;
}

static int threads_arg(const mxArray* a) {
    const int nt = scalar_int(a, "nThreads");
    return nt > 0 ? nt : 1;
}

// lines per piece of MATLAB's split (0 for no split: below 1, Inf or NaN)
static idx_t grain_arg(const mxArray* a) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "grain must be a real numeric scalar.");
    }
    const double v = mxGetScalar(a);
    if (!(v >= 1) || v > 1e15) {
        return 0;
    }
    return (idx_t)v;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
#ifdef __APPLE__
    // keep the mex file loaded until MATLAB exits: GCC's OpenMP runtime leaves
    // thread exit handlers in it, so unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if (nrhs < 1 || !mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "The first input must be 'otf', 'omw', 'rlstep' or 'fft1'.");
    }
    char method[16];
    mxGetString(prhs[0], method, sizeof(method));
    Timer tm;

    if (strcmp(method, "otf") == 0) {
        if (nrhs != 5) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Usage: H = decon_lucy_rfft_mex('otf', psf, outSize, nThreads, grain)");
        }
        const mxArray* psf = prhs[1];
        if (!mxIsSingle(psf) || mxIsComplex(psf) || mxIsEmpty(psf)) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "The PSF must be a real, non-empty single array.");
        }
        int ps[3], d[3];
        get_size(psf, ps);
        const mxArray* os = prhs[2];
        if (!mxIsDouble(os) || mxIsComplex(os) || mxGetNumberOfElements(os) < 2 || mxGetNumberOfElements(os) > 3) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "outSize must be a double vector of 2 or 3 elements.");
        }
        const double* osv = mxGetDoubles(os);
        for (int i = 0; i < 3; i++) {
            const double v = i < (int)mxGetNumberOfElements(os) ? osv[i] : 1;
            if (!(v >= 1) || v > 2147483647.0 || v != (double)(int)v) {
                mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "outSize must contain positive integers.");
            }
            d[i] = (int)v;
            if (ps[i] > d[i]) {
                mexErrMsgIdAndTxt("decon_lucy_rfft_mex:outSizeIsSmallerThanPsfSize", "outSizeIsSmallerThanPsfSize");
            }
        }
        check_shape(d);
        Ctx c;
        tm.start();
        setup(c, d, threads_arg(prhs[3]), grain_arg(prhs[4]));
        tm.stop(0);
        mwSize dh[3], ndh;
        half_size(d, dh, ndh);
        plhs[0] = mxCreateUninitNumericArray(ndh, dh, mxSINGLE_CLASS, mxCOMPLEX);
        otf_half(mxGetComplexSingles(plhs[0]), mxGetSingles(psf), ps, c, tm);
    } else if (strcmp(method, "omw") == 0) {
        if (nrhs != 8) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Usage: J = decon_lucy_rfft_mex('omw', I, H_f, H_b, nIter, background, nThreads, grain)");
        }
        const mxArray* im = prhs[1];
        if (!(mxIsUint16(im) || mxIsSingle(im)) || mxIsComplex(im) || mxIsEmpty(im)) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "The data must be a real, non-empty uint16 or single array.");
        }
        int d[3];
        get_size(im, d);
        check_shape(d);
        check_otf(prhs[2], d);
        check_otf(prhs[3], d);
        const int nIter = scalar_int(prhs[4], "nIter");
        const mxArray* bga = prhs[5];
        if (!mxIsNumeric(bga) || mxGetNumberOfElements(bga) != 1 || mxIsComplex(bga)) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "background must be a real numeric scalar.");
        }
        const double bg = mxGetScalar(bga);
        Ctx c;
        tm.start();
        setup(c, d, threads_arg(prhs[6]), grain_arg(prhs[7]));
        tm.stop(0);
        const mxComplexSingle* h1 = mxGetComplexSingles(prhs[2]);
        const mxComplexSingle* h2 = mxGetComplexSingles(prhs[3]);
        float* r;
        if (mxIsUint16(im)) {
            r = omw_iterations(mxGetUint16s(im), bg, h1, h2, nIter, c, tm);
        } else {
            r = omw_iterations(mxGetSingles(im), bg, h1, h2, nIter, c, tm);
        }
        plhs[0] = result_array(r, (idx_t)d[0] * d[1] * d[2], mxGetNumberOfDimensions(im), mxGetDimensions(im));
    } else if (strcmp(method, "rlstep") == 0) {
        if (nrhs != 6) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Usage: J = decon_lucy_rfft_mex('rlstep', Y, I, H, nThreads, grain)");
        }
        const mxArray* ya = prhs[1];
        const mxArray* ia = prhs[2];
        if (!mxIsSingle(ya) || mxIsComplex(ya) || mxIsEmpty(ya) || !mxIsSingle(ia) || mxIsComplex(ia)
                || mxGetNumberOfElements(ia) != mxGetNumberOfElements(ya)) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Y and I must be real single arrays of the same size.");
        }
        int d[3], di[3];
        get_size(ya, d);
        get_size(ia, di);
        if (d[0] != di[0] || d[1] != di[1] || d[2] != di[2]) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Y and I must be real single arrays of the same size.");
        }
        check_shape(d);
        check_otf(prhs[3], d);
        Ctx c;
        tm.start();
        setup(c, d, threads_arg(prhs[4]), grain_arg(prhs[5]));
        tm.stop(0);
        float* r = rl_step(mxGetSingles(ya), mxGetSingles(ia), mxGetComplexSingles(prhs[3]), c, tm);
        plhs[0] = result_array(r, (idx_t)d[0] * d[1] * d[2], mxGetNumberOfDimensions(ya), mxGetDimensions(ya));
    } else if (strcmp(method, "fft1") == 0) {
        if (nrhs != 4) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Usage: Y = decon_lucy_rfft_mex('fft1', X, nThreads, grain)");
        }
        const mxArray* xa = prhs[1];
        if (!mxIsSingle(xa) || mxIsComplex(xa) || mxGetNumberOfDimensions(xa) != 2 || mxGetM(xa) < 2 || mxGetN(xa) < 2) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "X must be a real single matrix with at least 2 rows and columns.");
        }
        int d[3];
        get_size(xa, d);
        Ctx c;
        setup(c, d, threads_arg(prhs[2]), grain_arg(prhs[3]));
        const mwSize dh[2] = {(mwSize)c.nyh, (mwSize)d[1]};
        plhs[0] = mxCreateUninitNumericArray(2, (mwSize*)dh, mxSINGLE_CLASS, mxCOMPLEX);
        fwd_y(c, mxGetSingles(xa), d[0], (float*)mxGetComplexSingles(plhs[0]));
    } else {
        mexErrMsgIdAndTxt("decon_lucy_rfft_mex:inputError", "Unknown method %s, use 'otf', 'omw', 'rlstep' or 'fft1'.", method);
    }

    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(1, 4, mxREAL);
        double* t = mxGetDoubles(plhs[1]);
        for (int i = 0; i < 4; i++) {
            t[i] = tm.t[i];
        }
    }
}
