// The parts of the normalized cross-correlation of normxcorr3_fast in a window of its 'full' output, with
// the same FFTs as MATLAB's own fftn and ifftn (as decon_lucy_rfft_mex) and the same window sums as its MATLAB
// code (box_sums), so that they equal normxcorr3_fast's bit for bit:
//   corr = real(ifftn(fftn(A, szOut) .* fftn(T(end : -1 : 1, end : -1 : 1, end : -1 : 1), szOut)))
//   sA   = sums of A over the windows of size szT (zero outside A), in double
//   sA2  = sums of A .* A (in double) over the windows
//   sT   = sums of the elements of T that overlap A at each window position, in double
// in the window, with szOut = szT + szA - 1.
//   fftn of real data: 1D real-to-complex transforms along y, complex ones along x and z (MATLAB's single
//     precision FFTW and its plans), then the Hermitian fill of MATLAB (an element whose mirror has the
//     lower linear index is the conjugate of the mirror, self-conjugate elements are real)
//   ifftn of the (Hermitian) product: complex transforms along y and x, then complex-to-real along z, each
//     divided by its length
//   MATLAB runs the 1D transforms of a stage in groups of four lines and some lines one by one (see
//     mark_alone); for some lengths the two give different last bits, so the mex runs the same lines alone.
// Only what the window needs is computed: the zero padding is not transformed, the forward transforms keep
// the half spectrum (k_z <= nz / 2) that the inverse uses (A's multiplied into T's as it is computed), the
// inverse transforms along x only the rows of the window and along z only its lines, and the window sums
// only the window's lines along x and z. 2D data (nz = 1) uses the same steps along y and x.
//
// [corr, sA, sA2, sT, t] = normxcorr3_window_mex('ncc', T, A, bbox, nThreads, grain[, fast])
//   T, A     : real single arrays (2D or 3D), size(T) <= size(A), szOut with at least 2 elements along y and x
//   bbox     : [s1, s2, s3, t1, t2, t3], the window in the full output (1-based, 1 <= s <= t <= szOut)
//   fast     : true for smaller FFTs instead of MATLAB's (default false): per dimension the smallest 2, 3, 5,
//              7-smooth size n >= size(A) with n >= t and n >= szOut - s + 1, for which the circular
//              correlation equals the linear one in the window; corr then differs from MATLAB's by rounding
//   corr     : single array of the window size
//   sA, sA2, sT : double arrays of the window size (computed only when requested)
//   t        : seconds in [planning, forward transform of T, of A, inverse transform, window sums]
//
// Y = normxcorr3_window_mex('fft1', X, nThreads, grain)
//   rows 1 .. n / 2 + 1 of fft(X, [], 1) for single X (n x m, n >= 2), with MATLAB's groups of lines for
//   the given grain (to find out which grain MATLAB uses)
//
//   nThreads : threads
//   grain    : most lines per piece of MATLAB's split of a stage among its threads: 32 when MATLAB runs
//              its FFTs in parallel, 0 when it doesn't (normxcorr3_fast finds it out)
//
// Build (-R2018a for interleaved complex arrays; no fused multiply-adds, as MATLAB computes products):
// mex -R2018a CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp -ffp-contract=off' LDFLAGS='$LDFLAGS -O3 -fopenmp' normxcorr3_window_mex.cpp -L<matlabroot>/bin/glnxa64 -lmwfftw3f
// (compileAllMexFiles builds it on Linux)

#include <stddef.h>
#include <stdint.h>
#include <string.h>
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

// sizes, plans, the lines MATLAB transforms alone, and per-thread scratch of one call
struct Ctx {
    int ny, nx, nz;      // FFT size (szOut; nz = 1 for 2D)
    bool is3d;
    int nyh;             // ny / 2 + 1, complex elements of the y-half of a column
    int nh;              // kept elements along the last dimension: nz / 2 + 1 (3D), nx / 2 + 1 (2D)
    idx_t ldr;           // floats per y-half column (2 * nyh)
    idx_t nP;            // complex elements of a half spectrum: [nh, ny, nx] (k_z first) for 3D, [ny, nh] for 2D
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
    // the lines MATLAB transforms alone, per stage of the full size (inverse stages as the forward ones)
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
            mexErrMsgIdAndTxt("normxcorr3_window_mex:outOfMemory", "Out of memory for the FFT buffers.");
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
            mexErrMsgIdAndTxt("normxcorr3_window_mex:fftwPlan", "FFTW planning failed.");
        }
    }
}

// a real volume zero-padded to the FFT size: data size v (3 elements), flipped along every dimension
struct Vol {
    const float* p;
    int v[3];
    bool flip;
};

// column j (x + v[1] * z) of the volume's data columns, zero-padded to ny floats
static inline void load_column(const Ctx& c, const Vol& V, const idx_t j, float* li) {
    const idx_t x = j % V.v[1], z = j / V.v[1];
    const int n = V.v[0];
    if (V.flip) {
        const float* src = V.p + ((idx_t)(V.v[1] - 1 - x) + (idx_t)V.v[1] * (V.v[2] - 1 - z)) * n;
        for (int y = 0; y < n; y++) li[y] = src[n - 1 - y];
    } else {
        memcpy(li, V.p + j * n, sizeof(float) * n);
    }
    for (int y = n; y < c.ny; y++) li[y] = 0.0f;
}

// forward stage 1: the y-half columns of the volume's data columns (x < v[1], z < v[2]; the others are zero)
// into Yh (column x + v[1] * z, ldr floats each)
static void fwd_y(Ctx& c, const Vol& V, float* Yh) {
    const idx_t vx = V.v[1];
    const idx_t ncol = vx * V.v[2];
    const idx_t nb = (ncol + kB - 1) / kB;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const idx_t j0 = b * kB;
        const int m = (int)(ncol - j0 < kB ? ncol - j0 : kB);
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        // the lanes of the block in a group (the lines of a group are transformed independently, so any
        // lines MATLAB doesn't transform alone can share one)
        int gl[kB], ng = 0;
        bool al[kB];
        for (int l = 0; l < m; l++) {
            const idx_t j = j0 + l;
            al[l] = c.aloneY[j % vx + (idx_t)c.nx * (j / vx)] != 0;
            if (!al[l]) gl[ng++] = l;
        }
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                load_column(c, V, j0 + gl[q < ng ? q : ng - 1], in + (idx_t)q * c.ny);
            }
            fftwf_execute_dft_r2c(c.fy[0], in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                memcpy(Yh + (j0 + gl[q]) * c.ldr, out + (idx_t)q * c.ldr, sizeof(float) * c.ldr);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            load_column(c, V, j0 + l, in);
            fftwf_execute_dft_r2c(c.fy[1], in, (fftwf_complex*)out);
            memcpy(Yh + (j0 + l) * c.ldr, out, sizeof(float) * c.ldr);
        }
    }
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

// line along x (row y, slice z < v[2]) of the full y range from the y-half columns, zero for x >= v[1]
static inline void x_line(const Ctx& c, const float* Yh, const int* v, const int y, const int z, float* li) {
    for (int x = 0; x < v[1]; x++) {
        yfull(Yh + ((idx_t)x + (idx_t)v[1] * z) * c.ldr, c.ny, c.nyh, y, li[2 * x], li[2 * x + 1]);
    }
    for (int x = v[1]; x < c.nx; x++) {
        li[2 * x] = 0.0f;
        li[2 * x + 1] = 0.0f;
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

// MATLAB's complex product a .* t (no fused multiply-add) into r (r may be t)
static inline void cmul(const float* a, const float* t, float* r) {
    const float re = a[0] * t[0] - a[1] * t[1];
    const float im = a[0] * t[1] + a[1] * t[0];
    r[0] = re;
    r[1] = im;
}

// Where the kept half spectrum of a forward transform goes: into P ([nh, ny, nx] for 3D, [ny, nh] for 2D), or
// for the second transform its product with P (A .* T, into P), except the self-conjugate planes (k_z = 0 and
// nz / 2, 2D: k_x = 0 and nx / 2), which wait in S (S[q], MATLAB's linear order) for MATLAB's Hermitian fill
struct Sink {
    fftwf_complex* P;
    float* S[2];
    bool product;
};

// the self-conjugate plane (0 or 1) of index k along the last dimension, or -1
static inline int side_plane(const Ctx& c, const int k) {
    const int n = c.is3d ? c.nz : c.nx;
    if (k == 0) return 0;
    if (n % 2 == 0 && k == n / 2) return 1;
    return -1;
}

// the kept outputs k_z < nh of the z line of row y, column k
static inline void z_out(const Ctx& c, const Sink& sk, const float* lo, const int y, const int k) {
    const idx_t L = (idx_t)y + (idx_t)c.ny * k;
    float* dst = (float*)sk.P + 2 * (idx_t)c.nh * L;
    if (!sk.product) {
        memcpy(dst, lo, sizeof(fftwf_complex) * c.nh);
        return;
    }
    for (int kz = 0; kz < c.nh; kz++) {
        const int q = side_plane(c, kz);
        if (q >= 0) {
            sk.S[q][2 * L] = lo[2 * kz];
            sk.S[q][2 * L + 1] = lo[2 * kz + 1];
        } else {
            cmul(lo + 2 * kz, dst + 2 * kz, dst + 2 * kz);
        }
    }
}

// a group of np z lines (band rows pl, columns pk), padded with its last line
static void z_group(const Ctx& c, const float* band, const int* pl, const int* pk, const int np, float* in, float* out,
        const Sink& sk, const int y0) {
    const int nz = c.nz, nx = c.nx;
    for (int q = 0; q < kB; q++) {
        const int s = q < np ? q : np - 1;
        memcpy(in + 2 * (idx_t)q * nz, band + 2 * ((idx_t)pl[s] * nx + pk[s]) * nz, sizeof(fftwf_complex) * nz);
    }
    fftwf_execute_dft(c.fz[0], (fftwf_complex*)in, (fftwf_complex*)out);
    for (int q = 0; q < np; q++) {
        z_out(c, sk, out + 2 * (idx_t)q * nz, y0 + pl[q], pk[q]);
    }
}

// forward stages 2 and 3 (3D): lines along x of the slices with data (z < v[2]), then all lines along z; the z
// outputs k_z < nh go to the sink. Rows are processed in bands of four; lines MATLAB transforms in groups are
// grouped in any fours (the lanes of a group are independent), the others transformed alone.
static void fwd_xz(Ctx& c, const float* Yh, const int* v, const Sink& sk) {
    const int ny = c.ny, nx = c.nx, nz = c.nz;
    const int nbands = (ny + kB - 1) / kB;
    // the band entries of the zero slices (z >= v[2]) stay zero
    for (int t = 0; t < c.nt; t++) {
        memset(c.band[t], 0, sizeof(fftwf_complex) * kB * (size_t)nx * nz);
    }
    #pragma omp parallel for num_threads(c.nt) schedule(dynamic, 1)
    for (int b = 0; b < nbands; b++) {
        const int t = omp_get_thread_num();
        const int y0 = b * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        float* band = (float*)c.band[t];
        // lines along x of the band rows, every slice with data
        for (int z = 0; z < v[2]; z++) {
            const idx_t L0 = (idx_t)y0 + (idx_t)ny * z;
            int gl[kB], ng = 0;
            for (int l = 0; l < m; l++) {
                if (!c.aloneX[L0 + l]) gl[ng++] = l;
            }
            if (ng > 0) {
                for (int q = 0; q < kB; q++) {
                    x_line(c, Yh, v, y0 + gl[q < ng ? q : ng - 1], z, in + 2 * (idx_t)q * nx);
                }
                fftwf_execute_dft(c.fx[0], (fftwf_complex*)in, (fftwf_complex*)out);
                for (int q = 0; q < ng; q++) {
                    to_band(band, out + 2 * (idx_t)q * nx, gl[q], z, nx, nz);
                }
            }
            for (int l = 0; l < m; l++) {
                if (!c.aloneX[L0 + l]) continue;
                x_line(c, Yh, v, y0 + l, z, in);
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
                    z_out(c, sk, out, y0 + l, k);
                    continue;
                }
                pl[np] = l;
                pk[np] = k;
                np++;
                if (np == kB) {
                    z_group(c, band, pl, pk, np, in, out, sk, y0);
                    np = 0;
                }
            }
        }
        if (np > 0) {
            z_group(c, band, pl, pk, np, in, out, sk, y0);
        }
    }
}

// the kept outputs k_x < nh of the x transform of row y (2D)
static inline void x2d_out(const Ctx& c, const Sink& sk, const float* lo, const int y) {
    float* dst = (float*)sk.P + 2 * (idx_t)y;
    for (int k = 0; k < c.nh; k++) {
        float* e = dst + 2 * (idx_t)k * c.ny;
        const int q = sk.product ? side_plane(c, k) : -1;
        if (!sk.product) {
            e[0] = lo[2 * k];
            e[1] = lo[2 * k + 1];
        } else if (q >= 0) {
            sk.S[q][2 * y] = lo[2 * k];
            sk.S[q][2 * y + 1] = lo[2 * k + 1];
        } else {
            cmul(lo + 2 * k, e, e);
        }
    }
}

// forward stage 2 (2D): lines along x (rows)
static void fwd_x2d(Ctx& c, const float* Yh, const int* v, const Sink& sk) {
    const int ny = c.ny, nx = c.nx;
    const int nb = (ny + kB - 1) / kB;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (int b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const int y0 = b * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        int gl[kB], ng = 0;
        for (int l = 0; l < m; l++) {
            if (!c.aloneX[y0 + l]) gl[ng++] = l;
        }
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                x_line(c, Yh, v, y0 + gl[q < ng ? q : ng - 1], 0, in + 2 * (idx_t)q * nx);
            }
            fftwf_execute_dft(c.fx[0], (fftwf_complex*)in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                x2d_out(c, sk, out + 2 * (idx_t)q * nx, y0 + gl[q]);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneX[y0 + l]) continue;
            x_line(c, Yh, v, y0 + l, 0, in);
            fftwf_execute_dft(c.fx[1], (fftwf_complex*)in, (fftwf_complex*)out);
            x2d_out(c, sk, out, y0 + l);
        }
    }
}

// MATLAB's Hermitian fill of a self-conjugate plane (3D: ny x nx, 2D: a column of ny) whose element with
// MATLAB's linear index L in the plane is at base + 2 * st * L: an element whose mirror has the lower linear
// index becomes the mirror's conjugate, a self-conjugate element its real part
static void fill_plane(const Ctx& c, float* base, const idx_t st) {
    const int ny = c.ny, nx = c.is3d ? c.nx : 1;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (int k2 = 0; k2 < nx; k2++) {
        const int m2 = (nx - k2) % nx;
        for (int k1 = 0; k1 < ny; k1++) {
            const int m1 = (ny - k1) % ny;
            const idx_t L = k1 + (idx_t)ny * k2, Lm = m1 + (idx_t)ny * m2;
            if (Lm < L) {
                base[2 * st * L] = base[2 * st * Lm];
                base[2 * st * L + 1] = -base[2 * st * Lm + 1];
            } else if (Lm == L) {
                base[2 * st * L + 1] = 0.0f;
            }
        }
    }
}

// fftn of the zero-padded real volume into the sink, with MATLAB's Hermitian fill (for the product, of the
// self-conjugate planes of this transform before they are multiplied); Yh: y-half buffer of at least
// ldr * v[1] * v[2] floats
static void fftn_half(Ctx& c, const Vol& V, float* Yh, const Sink& sk) {
    fwd_y(c, V, Yh);
    if (c.is3d) {
        fwd_xz(c, Yh, V.v, sk);
    } else {
        fwd_x2d(c, Yh, V.v, sk);
    }
    const int n = c.is3d ? c.nz : c.nx;
    const idx_t st = c.is3d ? c.nh : 1;
    const idx_t np = c.is3d ? (idx_t)c.ny * c.nx : (idx_t)c.ny;
    for (int q = 0; q < (n % 2 == 0 ? 2 : 1); q++) {
        const int k = q == 0 ? 0 : n / 2;
        float* pb = (float*)sk.P + 2 * (c.is3d ? (idx_t)k : (idx_t)k * c.ny);
        if (!sk.product) {
            fill_plane(c, pb, st);
            continue;
        }
        fill_plane(c, sk.S[q], 1);
        #pragma omp parallel for num_threads(c.nt) schedule(static)
        for (idx_t L = 0; L < np; L++) {
            cmul(sk.S[q] + 2 * L, pb + 2 * st * L, pb + 2 * st * L);
        }
    }
}

// The y lines of the kept half spectrum, numbered i with k_z fastest (3D: k_z = i % nh, k_x = i / nh; 2D:
// k_x = i): MATLAB's line number (k_x + nx * k_z), the offset of the line's first element in P and the
// distance of its elements
static inline idx_t yline_matlab(const Ctx& c, const idx_t i) {
    return c.is3d ? i / c.nh + (idx_t)c.nx * (i % c.nh) : i;
}
static inline idx_t yline_base(const Ctx& c, const idx_t i) {
    return c.is3d ? i % c.nh + (idx_t)c.nh * c.ny * (i / c.nh) : (idx_t)c.ny * i;
}

// y line i of P into li
static inline void y_line(const Ctx& c, const fftwf_complex* P, const idx_t i, float* li) {
    const idx_t st = c.is3d ? c.nh : 1;
    const float* src = (const float*)P + 2 * yline_base(c, i);
    for (int k = 0; k < c.ny; k++) {
        li[2 * k] = src[2 * k * st];
        li[2 * k + 1] = src[2 * k * st + 1];
    }
}

// the inverse transform of y line i (lo), divided by ny, into P
static inline void y_back(const Ctx& c, fftwf_complex* P, const float* lo, const idx_t i, const float fn) {
    const idx_t st = c.is3d ? c.nh : 1;
    float* pc = (float*)P + 2 * yline_base(c, i);
    for (int k = 0; k < c.ny; k++) {
        pc[2 * k * st] = lo[2 * k] / fn;
        pc[2 * k * st + 1] = lo[2 * k + 1] / fn;
    }
}

// inverse stage 1: P along y, c2c, divided by ny, in place. Only the kept lines (the first ones of MATLAB's
// nx * nz) are transformed.
static void inv_y(Ctx& c, fftwf_complex* P) {
    const idx_t nl = c.is3d ? (idx_t)c.nx * c.nh : (idx_t)c.nh;
    const idx_t nb = (nl + kB - 1) / kB;
    const float fn = (float)c.ny;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        const idx_t i0 = b * kB;
        const int m = (int)(nl - i0 < kB ? nl - i0 : kB);
        int gl[kB], ng = 0;
        bool al[kB];
        for (int l = 0; l < m; l++) {
            al[l] = c.aloneY[yline_matlab(c, i0 + l)] != 0;
            if (!al[l]) gl[ng++] = l;
        }
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                y_line(c, P, i0 + gl[q < ng ? q : ng - 1], in + 2 * (idx_t)q * c.ny);
            }
            fftwf_execute_dft(c.iy[0], (fftwf_complex*)in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                y_back(c, P, out + 2 * (idx_t)q * c.ny, i0 + gl[q], fn);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            y_line(c, P, i0 + l, in);
            fftwf_execute_dft(c.iy[1], (fftwf_complex*)in, (fftwf_complex*)out);
            y_back(c, P, out, i0 + l, fn);
        }
    }
}

// the window of the output: w0 .. w1 (0-based, inclusive) and its size wn per dimension
struct Win {
    int w0[3], w1[3], wn[3];
};

// x line (window row r, plane kz) of P [nh, ny, nx] into li
static inline void x_line_win(const Ctx& c, const fftwf_complex* P, const Win& W, const int r, const int kz, float* li) {
    const idx_t st = (idx_t)c.nh * c.ny;
    const float* src = (const float*)P + 2 * (kz + (idx_t)c.nh * (W.w0[0] + r));
    for (int x = 0; x < c.nx; x++) {
        li[2 * x] = src[2 * x * st];
        li[2 * x + 1] = src[2 * x * st + 1];
    }
}

// the window columns of the inverse x transform (lo) of window row r, plane kz, divided by nx, into Q [nh, wy, wx]
static inline void x_back_win(const Ctx& c, const Win& W, const float* lo, const int r, const int kz, const float fn, float* Q) {
    const idx_t st = (idx_t)c.nh * W.wn[0];
    float* dst = Q + 2 * (kz + (idx_t)c.nh * r);
    for (int x = 0; x < W.wn[1]; x++) {
        dst[2 * x * st] = lo[2 * (W.w0[1] + x)] / fn;
        dst[2 * x * st + 1] = lo[2 * (W.w0[1] + x) + 1] / fn;
    }
}

// inverse stage 2 (3D): lines along x of the window rows in the kept z planes, c2c, divided by nx; the window
// columns go to Q [nh, wy, wx]
static void inv_x_win(Ctx& c, const fftwf_complex* P, const Win& W, float* Q) {
    const int wy = W.wn[0];
    // lines numbered i = kz + nh * r (window row r); MATLAB numbers them y + ny * kz
    const idx_t nl = (idx_t)c.nh * wy;
    const idx_t nb = (nl + kB - 1) / kB;
    const float fn = (float)c.nx;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        const idx_t i0 = b * kB;
        const int m = (int)(nl - i0 < kB ? nl - i0 : kB);
        int gl[kB], ng = 0;
        bool al[kB];
        for (int l = 0; l < m; l++) {
            const idx_t i = i0 + l;
            al[l] = c.aloneX[(idx_t)W.w0[0] + i / c.nh + (idx_t)c.ny * (i % c.nh)] != 0;
            if (!al[l]) gl[ng++] = l;
        }
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                const idx_t i = i0 + gl[q < ng ? q : ng - 1];
                x_line_win(c, P, W, (int)(i / c.nh), (int)(i % c.nh), in + 2 * (idx_t)q * c.nx);
            }
            fftwf_execute_dft(c.ix[0], (fftwf_complex*)in, (fftwf_complex*)out);
            for (int q = 0; q < ng; q++) {
                const idx_t i = i0 + gl[q];
                x_back_win(c, W, out + 2 * (idx_t)q * c.nx, (int)(i / c.nh), (int)(i % c.nh), fn, Q);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            const idx_t i = i0 + l;
            x_line_win(c, P, W, (int)(i / c.nh), (int)(i % c.nh), in);
            fftwf_execute_dft(c.ix[1], (fftwf_complex*)in, (fftwf_complex*)out);
            x_back_win(c, W, out, (int)(i / c.nh), (int)(i % c.nh), fn, Q);
        }
    }
}

// the input of window line L of the inverse last stage (nh complex values): along z from Q [nh, wy, wx] (3D),
// along x from row w0[0] + L of P [ny, nh] (2D)
static inline void last_line(const Ctx& c, const fftwf_complex* P, const float* Q, const Win& W, const idx_t L, float* li) {
    const int nh = c.nh;
    if (c.is3d) {
        memcpy(li, Q + 2 * (idx_t)nh * L, sizeof(fftwf_complex) * nh);
    } else {
        const float* src = (const float*)P + 2 * (W.w0[0] + L);
        for (int k = 0; k < nh; k++) {
            li[2 * k] = src[2 * (idx_t)k * c.ny];
            li[2 * k + 1] = src[2 * (idx_t)k * c.ny + 1];
        }
    }
}

// the window part of the real line L (lo) divided by its length, into R [wy, wx, wz]
static inline void last_out(const Ctx& c, const Win& W, const float* lo, const idx_t L, float* R) {
    const int o0 = c.is3d ? W.w0[2] : W.w0[1];
    const int on = c.is3d ? W.wn[2] : W.wn[1];
    const idx_t ostep = c.is3d ? (idx_t)W.wn[0] * W.wn[1] : (idx_t)W.wn[0];
    const float fn = (float)(c.is3d ? c.nz : c.nx);
    for (int k = 0; k < on; k++) {
        R[L + (idx_t)k * ostep] = lo[o0 + k] / fn;
    }
}

// inverse last stage: c2r along z of the window lines (3D) or along x of the window rows (2D), divided by
// the length; the window part goes to R
static void inv_last_win(Ctx& c, const fftwf_complex* P, const float* Q, const Win& W, float* R) {
    const int nh = c.nh;
    const int nlen = c.is3d ? c.nz : c.nx;
    const int wy = W.wn[0];
    // window lines: (y, x) pairs (3D) or rows (2D); MATLAB numbers them y + ny * x (3D) or y (2D)
    const idx_t nl = c.is3d ? (idx_t)wy * W.wn[1] : (idx_t)wy;
    const uint8_t* alone = c.is3d ? c.aloneZ : c.aloneX;
    const idx_t nb = (nl + kB - 1) / kB;
    #pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        float* in = (float*)c.in[t];
        float* out = (float*)c.out[t];
        const idx_t L0 = b * kB;
        const int m = (int)(nl - L0 < kB ? nl - L0 : kB);
        int gl[kB], ng = 0;
        bool al[kB];
        for (int l = 0; l < m; l++) {
            const idx_t L = L0 + l;
            const idx_t y = W.w0[0] + L % wy;
            al[l] = alone[c.is3d ? y + (idx_t)c.ny * (W.w0[1] + L / wy) : y] != 0;
            if (!al[l]) gl[ng++] = l;
        }
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                last_line(c, P, Q, W, L0 + gl[q < ng ? q : ng - 1], in + 2 * (idx_t)q * nh);
            }
            fftwf_execute_dft_c2r(c.il[0], (fftwf_complex*)in, out);
            for (int q = 0; q < ng; q++) {
                last_out(c, W, out + (idx_t)q * nlen, L0 + gl[q], R);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            last_line(c, P, Q, W, L0 + l, in);
            fftwf_execute_dft_c2r(c.il[1], (fftwf_complex*)in, out);
            last_out(c, W, out, L0 + l, R);
        }
    }
}

// The window sums of normxcorr3_fast's MATLAB code in double (its box_sums): along each dimension d, for the
// window's positions k (k = 0 .. w[d] - 1), the sum of the elements lo[d][k] .. hi[d][k] (1-based, inclusive)
// of the previous stage, as MATLAB's cumsum from the first element (P(1) = x(1), P(i) = P(i - 1) + x(i)) and
// P(hi) - P(lo - 1) (P(hi) for lo = 1), along y, then x, then z, so that the sums are MATLAB's bit for bit.
// The elements are X, or X .* X in double (sq). n: size of X; out: [w0, w1, w2].
static void box_sums_window(const float* X, const int* n, int* const* lo, int* const* hi, const int* w, const bool sq,
                            double* out, const int nt) {
    int H[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
        for (int k = 0; k < w[d]; k++) H[d] = hi[d][k] > H[d] ? hi[d][k] : H[d];
    }
    // stage 1 (y) of every column: S1 [w0, n1, n2]
    const idx_t ncol = (idx_t)n[1] * n[2];
    double* S1 = (double*)mxMalloc(sizeof(double) * (size_t)w[0] * ncol);
    double* pb = (double*)mxMalloc(sizeof(double) * (size_t)H[0] * nt);
    #pragma omp parallel num_threads(nt)
    {
        double* P = pb + (size_t)H[0] * omp_get_thread_num();
        #pragma omp for schedule(static)
        for (idx_t c = 0; c < ncol; c++) {
            const float* f = X + c * n[0];
            for (int i = 0; i < H[0]; i++) {
                const double v = sq ? (double)f[i] * (double)f[i] : (double)f[i];
                P[i] = i == 0 ? v : P[i - 1] + v;
            }
            double* o = S1 + c * w[0];
            for (int k = 0; k < w[0]; k++) {
                o[k] = lo[0][k] > 1 ? P[hi[0][k] - 1] - P[lo[0][k] - 2] : P[hi[0][k] - 1];
            }
        }
    }
    mxFree(pb);
    // stage 2 (x) of the rows of every slice: S2 [w0, w1, n2]; P holds the prefix sums of the w0 rows
    const idx_t w0 = w[0];
    double* S2 = (double*)mxMalloc(sizeof(double) * (size_t)w0 * w[1] * n[2]);
    pb = (double*)mxMalloc(sizeof(double) * (size_t)w0 * H[1] * nt);
    #pragma omp parallel num_threads(nt)
    {
        double* P = pb + (size_t)w0 * H[1] * omp_get_thread_num();
        #pragma omp for schedule(static)
        for (idx_t z = 0; z < n[2]; z++) {
            const double* src = S1 + w0 * n[1] * z;
            for (int i = 0; i < H[1]; i++) {
                const double* s = src + w0 * i;
                double* p = P + w0 * i;
                if (i == 0) {
                    for (idx_t y = 0; y < w0; y++) p[y] = s[y];
                } else {
                    const double* q = p - w0;
                    for (idx_t y = 0; y < w0; y++) p[y] = q[y] + s[y];
                }
            }
            double* o = S2 + w0 * w[1] * z;
            for (int k = 0; k < w[1]; k++) {
                const double* ph = P + w0 * (hi[1][k] - 1);
                double* ok = o + w0 * k;
                if (lo[1][k] > 1) {
                    const double* pl = P + w0 * (lo[1][k] - 2);
                    for (idx_t y = 0; y < w0; y++) ok[y] = ph[y] - pl[y];
                } else {
                    for (idx_t y = 0; y < w0; y++) ok[y] = ph[y];
                }
            }
        }
    }
    mxFree(pb);
    mxFree(S1);
    // stage 3 (z) of the lines of every window column: out [w0, w1, w2]
    const idx_t wyx = w0 * w[1];
    pb = (double*)mxMalloc(sizeof(double) * (size_t)w0 * H[2] * nt);
    #pragma omp parallel num_threads(nt)
    {
        double* P = pb + (size_t)w0 * H[2] * omp_get_thread_num();
        #pragma omp for schedule(static)
        for (idx_t x = 0; x < w[1]; x++) {
            for (int i = 0; i < H[2]; i++) {
                const double* s = S2 + w0 * x + wyx * i;
                double* p = P + w0 * i;
                if (i == 0) {
                    for (idx_t y = 0; y < w0; y++) p[y] = s[y];
                } else {
                    const double* q = p - w0;
                    for (idx_t y = 0; y < w0; y++) p[y] = q[y] + s[y];
                }
            }
            for (int k = 0; k < w[2]; k++) {
                const double* ph = P + w0 * (hi[2][k] - 1);
                double* ok = out + w0 * x + wyx * k;
                if (lo[2][k] > 1) {
                    const double* pl = P + w0 * (lo[2][k] - 2);
                    for (idx_t y = 0; y < w0; y++) ok[y] = ph[y] - pl[y];
                } else {
                    for (idx_t y = 0; y < w0; y++) ok[y] = ph[y];
                }
            }
        }
    }
    mxFree(pb);
    mxFree(S2);
}

// buffers of one call (fftwf_malloc), freed at the end and on errors
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
        b[slot] = fftwf_malloc(bytes > 0 ? bytes : 1);
        if (b[slot] == NULL) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:outOfMemory", "Out of memory for the FFT buffers.");
        }
        return b[slot];
    }
};

// array size as 3 int dimensions
static void get_size(const mxArray* a, int* d, const char* name) {
    const mwSize nd = mxGetNumberOfDimensions(a);
    const mwSize* dims = mxGetDimensions(a);
    if (nd > 3) {
        mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "%s must have at most 3 dimensions.", name);
    }
    for (int i = 0; i < 3; i++) {
        const mwSize v = i < (int)nd ? dims[i] : 1;
        if (v > 1073741823) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "The dimensions of %s must be smaller than 2^30.", name);
        }
        d[i] = (int)v;
    }
}

static void check_real_single(const mxArray* a, const char* name) {
    if (!mxIsSingle(a) || mxIsComplex(a) || mxIsEmpty(a)) {
        mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "%s must be a real, non-empty single array.", name);
    }
}

// count from a numeric scalar (at least 1)
static int threads_arg(const mxArray* a) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "nThreads must be a real numeric scalar.");
    }
    const double v = mxGetScalar(a);
    if (!(v >= 1)) {
        return 1;
    }
    return v > 4096 ? 4096 : (int)v;
}

// smallest n >= m whose prime factors are 2, 3, 5 and 7
static int smooth_size(const int m) {
    for (int n = m;; n++) {
        int r = n;
        const int ps[4] = {2, 3, 5, 7};
        for (int k = 0; k < 4; k++) {
            while (r % ps[k] == 0) r /= ps[k];
        }
        if (r == 1) return n;
    }
}

// lines per piece of MATLAB's split (0 for no split: below 1, Inf or NaN)
static idx_t grain_arg(const mxArray* a) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "grain must be a real numeric scalar.");
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
        mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "The first input must be 'ncc' or 'fft1'.");
    }
    char method[16];
    mxGetString(prhs[0], method, sizeof(method));

    if (strcmp(method, "ncc") == 0) {
        if (nrhs != 6 && nrhs != 7) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "Usage: [corr, sA, sA2, sT, t] = normxcorr3_window_mex('ncc', T, A, bbox, nThreads, grain[, fast])");
        }
        // fast: smaller FFTs (2, 3, 5, 7-smooth sizes that just avoid wrap-around in the window), not MATLAB's
        const bool fast = nrhs == 7 && mxGetNumberOfElements(prhs[6]) == 1 && mxGetScalar(prhs[6]) != 0;
        if (nlhs > 5) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:outputError", "At most 5 outputs.");
        }
        check_real_single(prhs[1], "T");
        check_real_single(prhs[2], "A");
        int ts[3], a[3], d[3];
        get_size(prhs[1], ts, "T");
        get_size(prhs[2], a, "A");
        for (int i = 0; i < 3; i++) {
            if (ts[i] > a[i]) {
                mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "The template must not be larger than the image.");
            }
            d[i] = ts[i] + a[i] - 1;
        }
        if (d[0] < 2 || d[1] < 2) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "The full output must have at least 2 elements along y and x.");
        }
        const mxArray* bb = prhs[3];
        if (!mxIsDouble(bb) || mxIsComplex(bb) || mxGetNumberOfElements(bb) != 6) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "bbox must be a double vector of 6 elements.");
        }
        const double* bv = mxGetDoubles(bb);
        Win W;
        for (int i = 0; i < 3; i++) {
            const double s = bv[i], e = bv[i + 3];
            if (!(s >= 1) || !(e <= d[i]) || !(s <= e) || s != floor(s) || e != floor(e)) {
                mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "bbox must hold integers with 1 <= s <= t <= size of the full output.");
            }
            W.w0[i] = (int)s - 1;
            W.w1[i] = (int)e - 1;
            W.wn[i] = W.w1[i] - W.w0[i] + 1;
        }
        const int nt = threads_arg(prhs[4]);
        idx_t grain = grain_arg(prhs[5]);
        if (fast) {
            // the circular correlation of size n equals the linear one at the window if n >= t + 1 and
            // n >= full size - s (0-based s, t); n >= size(A) keeps the inputs unwrapped
            for (int i = 0; i < 3; i++) {
                int n = d[i] - W.w0[i] > W.w1[i] + 1 ? d[i] - W.w0[i] : W.w1[i] + 1;
                n = n > a[i] ? n : a[i];
                d[i] = d[i] > 1 ? smooth_size(n) : 1;
            }
            grain = 0;
        }
        double tm[5];
        double t0 = omp_get_wtime();
        Ctx c;
        setup(c, d, nt, grain);
        tm[0] = omp_get_wtime() - t0;

        const mwSize od[3] = {(mwSize)W.wn[0], (mwSize)W.wn[1], (mwSize)W.wn[2]};
        plhs[0] = mxCreateNumericArray(3, od, mxSINGLE_CLASS, mxREAL);

        // the y-half columns, the half spectrum (of T, then of the product), the self-conjugate planes of A's
        Buffers buf;
        const idx_t ncolMax = (idx_t)a[1] * a[2] > (idx_t)ts[1] * ts[2] ? (idx_t)a[1] * a[2] : (idx_t)ts[1] * ts[2];
        float* Yh = (float*)buf.get(0, sizeof(float) * (size_t)c.ldr * ncolMax);
        fftwf_complex* P = (fftwf_complex*)buf.get(1, sizeof(fftwf_complex) * (size_t)c.nP);
        const idx_t nplane = c.is3d ? (idx_t)c.ny * c.nx : (idx_t)c.ny;
        float* S = (float*)buf.get(2, sizeof(fftwf_complex) * 2 * (size_t)nplane);
        const Vol VT = {mxGetSingles(prhs[1]), {ts[0], ts[1], ts[2]}, true};
        const Vol VA = {mxGetSingles(prhs[2]), {a[0], a[1], a[2]}, false};
        t0 = omp_get_wtime();
        const Sink skT = {P, {S, S + 2 * nplane}, false};
        fftn_half(c, VT, Yh, skT);
        tm[1] = omp_get_wtime() - t0;
        t0 = omp_get_wtime();
        const Sink skA = {P, {S, S + 2 * nplane}, true};
        fftn_half(c, VA, Yh, skA);
        tm[2] = omp_get_wtime() - t0;
        t0 = omp_get_wtime();
        inv_y(c, P);
        float* R = mxGetSingles(plhs[0]);
        if (c.is3d) {
            float* Q = (float*)buf.get(3, sizeof(fftwf_complex) * (size_t)W.wn[0] * W.wn[1] * c.nh);
            inv_x_win(c, P, W, Q);
            inv_last_win(c, P, Q, W, R);
        } else {
            inv_last_win(c, P, NULL, W, R);
        }
        tm[3] = omp_get_wtime() - t0;
        // the window sums (only the requested ones): at window position k (1-based in the full output) the
        // image's elements max(1, k - ts + 1) .. min(k, a) and the template's max(1, ts - k + 1) ..
        // min(ts, ts - k + a) overlap
        t0 = omp_get_wtime();
        if (nlhs > 1) {
            int* lohi = (int*)mxMalloc(sizeof(int) * 4 * (size_t)(W.wn[0] + W.wn[1] + W.wn[2]));
            int* loA[3];
            int* hiA[3];
            int* loT[3];
            int* hiT[3];
            int* q = lohi;
            for (int i = 0; i < 3; i++) {
                loA[i] = q;
                hiA[i] = q + W.wn[i];
                loT[i] = q + 2 * W.wn[i];
                hiT[i] = q + 3 * W.wn[i];
                q += 4 * W.wn[i];
                for (int j = 0; j < W.wn[i]; j++) {
                    const int k = W.w0[i] + j + 1;
                    loA[i][j] = k - ts[i] + 1 > 1 ? k - ts[i] + 1 : 1;
                    hiA[i][j] = k < a[i] ? k : a[i];
                    loT[i][j] = ts[i] - k + 1 > 1 ? ts[i] - k + 1 : 1;
                    hiT[i][j] = ts[i] < ts[i] - k + a[i] ? ts[i] : ts[i] - k + a[i];
                }
            }
            for (int o = 1; o < nlhs && o < 4; o++) {
                plhs[o] = mxCreateNumericArray(3, od, mxDOUBLE_CLASS, mxREAL);
                if (o < 3) {
                    box_sums_window(mxGetSingles(prhs[2]), a, loA, hiA, W.wn, o == 2, mxGetDoubles(plhs[o]), nt);
                } else {
                    box_sums_window(mxGetSingles(prhs[1]), ts, loT, hiT, W.wn, false, mxGetDoubles(plhs[o]), nt);
                }
            }
            mxFree(lohi);
        }
        tm[4] = omp_get_wtime() - t0;
        if (nlhs > 4) {
            plhs[4] = mxCreateDoubleMatrix(1, 5, mxREAL);
            memcpy(mxGetDoubles(plhs[4]), tm, sizeof(tm));
        }
    } else if (strcmp(method, "fft1") == 0) {
        if (nrhs != 4) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "Usage: Y = normxcorr3_window_mex('fft1', X, nThreads, grain)");
        }
        const mxArray* xa = prhs[1];
        if (!mxIsSingle(xa) || mxIsComplex(xa) || mxGetNumberOfDimensions(xa) != 2 || mxGetM(xa) < 2 || mxGetN(xa) < 2) {
            mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "X must be a real single matrix with at least 2 rows and columns.");
        }
        int d[3];
        get_size(xa, d, "X");
        Ctx c;
        setup(c, d, threads_arg(prhs[2]), grain_arg(prhs[3]));
        const mwSize dh[2] = {(mwSize)c.nyh, (mwSize)d[1]};
        plhs[0] = mxCreateUninitNumericArray(2, (mwSize*)dh, mxSINGLE_CLASS, mxCOMPLEX);
        const Vol V = {mxGetSingles(xa), {d[0], d[1], 1}, false};
        fwd_y(c, V, (float*)mxGetComplexSingles(plhs[0]));
    } else {
        mexErrMsgIdAndTxt("normxcorr3_window_mex:inputError", "Unknown method %s, use 'ncc' or 'fft1'.", method);
    }
}
