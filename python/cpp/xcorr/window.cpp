// normxcorr3_window_mex of PetaKit5D (MATLAB 4a3bef3) on plain arrays in MATLAB's layout: its loops are the
// mex's, unchanged; the gateway is replaced by xcorr.h's functions, MATLAB's allocation by the C library's, its
// errors by exceptions, and the FFTW calls go through the fft group's backend table (backend.h): the exact
// backend (FFTW 3.3.8 with MATLAB's codelets, the optional add-on) gives MATLAB's bits, the permissive one
// (pocketfft) the same transforms to rounding.
//
//   fftn of real data: 1D real-to-complex transforms along y, complex ones along x and z (MATLAB's single
//     precision plans), then the Hermitian fill of MATLAB (an element whose mirror has the lower linear index
//     is the conjugate of the mirror, self-conjugate elements are real)
//   ifftn of the (Hermitian) product: complex transforms along y and x, then complex-to-real along z, each
//     divided by its length
//   MATLAB runs the 1D transforms of a stage in groups of four lines and some lines one by one (see
//     mark_alone); for some lengths the two give different last bits, so the same lines run alone here.
// Only what the window needs is computed: the zero padding is not transformed, the forward transforms keep
// the half spectrum (k_z <= nz / 2) that the inverse uses (A0's multiplied into T0's as it is computed), the
// inverse transforms along x only the rows of the window and along z only its lines, and the window sums only
// the window's lines along x and z. 2D data (nz = 1) uses the same steps along y and x. The results don't
// depend on the thread count.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "backend.h"
#include "xcorr.h"

namespace pk5d::xcorr {
namespace {

using pk5d::fft::FftwApi;

typedef float fcomplex[2];
typedef void* plan_t;

// transforms per FFTW execute, as in MATLAB's fft
constexpr int kB = 4;

// byte offsets (modulo 64) of the scratch buffers, as MATLAB's
constexpr int kInOff = 16;
constexpr int kOutOff = 32;

inline int thread_num() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

inline double wtime() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// sizes, plans, the lines MATLAB transforms alone, and per-thread scratch of one call
struct Ctx {
    int ny, nx, nz;  // FFT size (szOut; nz = 1 for 2D)
    bool is3d;
    int nyh;    // ny / 2 + 1, complex elements of the y-half of a column
    int nh;     // kept elements along the last dimension: nz / 2 + 1 (3D), nx / 2 + 1 (2D)
    idx_t ldr;  // floats per y-half column (2 * nyh)
    idx_t nP;   // complex elements of a half spectrum: [nh, ny, nx] (k_z first) for 3D, [ny, nh] for 2D
    int nt;
    const FftwApi* F;
    // plans [0]: groups of four, [1]: single transforms
    plan_t fy[2], fx[2], fz[2], iy[2], ix[2], il[2];
    // nonzero for the lines MATLAB transforms alone: along y (column x + nx * z), along x (row y + ny * z),
    // along z (y + ny * x, 3D only)
    std::vector<uint8_t> aloneY, aloneX, aloneZ;
    // scratch per thread: in, out (bases for free) and the x-line band buffer of the 3D forward
    std::vector<void*> base;
    std::vector<fcomplex*> in, out, band;
    Ctx() : nt(0), F(nullptr) {
        for (int i = 0; i < 2; i++) fy[i] = fx[i] = fz[i] = iy[i] = ix[i] = il[i] = nullptr;
    }
    ~Ctx() { release(); }
    Ctx(const Ctx&) = delete;
    Ctx& operator=(const Ctx&) = delete;
    // the plans belong to the process's plan cache (pk5d::fft::cached_plan): only the scratch is freed
    void release() {
        for (void*& b : base) {
            std::free(b);
            b = nullptr;
        }
    }
};

// Marks the lines a .. b - 1 of a stage that MATLAB transforms alone. MATLAB's FFT splits the lines of a
// stage among its threads by halving them (at a + (b - a) / 2) until a piece has at most grain lines (when
// it runs in parallel and the stage is large enough, see stage_grain; grain 0: one piece), and runs a piece
// in groups of four lines from its start, the last (length % 4) lines one by one. Lines along x of 3D data
// (seg = ny) also restart the groups at each z slice.
void mark_alone(uint8_t* alone, const idx_t a, const idx_t b, const idx_t grain, const idx_t seg) {
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
idx_t stage_grain(const idx_t grain, const int n, const idx_t N) {
    int bits = 0;
    for (int v = n; v > 0; v >>= 1) bits++;
    return (grain > 0 && (double)n * (double)N * bits > 98304.0) ? grain : 0;
}

void* alloc_off(const size_t bytes, const int off, void** base) {
    *base = std::malloc(bytes + 128);
    if (*base == nullptr) return nullptr;
    uintptr_t p = ((uintptr_t)(*base) + 63) & ~(uintptr_t)63;
    return (void*)(p + off);
}

void setup(Ctx& c, const int* d, const int nThreads, const idx_t grain, const int be) {
    // single precision of the backend
    c.F = &pk5d::fft::backend(be).api[0];
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
    c.aloneY.assign(nlY, 0);
    c.aloneX.assign(nlX, 0);
    mark_alone(c.aloneY.data(), 0, nlY, stage_grain(grain, c.ny, nlY), 0);
    mark_alone(c.aloneX.data(), 0, nlX, stage_grain(grain, c.nx, nlX), c.is3d ? c.ny : 0);
    if (c.is3d) {
        const idx_t nlZ = (idx_t)c.ny * c.nx;
        c.aloneZ.assign(nlZ, 0);
        mark_alone(c.aloneZ.data(), 0, nlZ, stage_grain(grain, c.nz, nlZ), 0);
    }
    const int maxn = c.ny > c.nx ? (c.ny > c.nz ? c.ny : c.nz) : (c.nx > c.nz ? c.nx : c.nz);
    const size_t sbytes = sizeof(fcomplex) * kB * (size_t)maxn;
    c.base.assign(3 * (size_t)nThreads, nullptr);
    c.in.assign(nThreads, nullptr);
    c.out.assign(nThreads, nullptr);
    c.band.assign(nThreads, nullptr);
    for (int t = 0; t < nThreads; t++) {
        c.in[t] = (fcomplex*)alloc_off(sbytes, kInOff, &c.base[3 * t]);
        c.out[t] = (fcomplex*)alloc_off(sbytes, kOutOff, &c.base[3 * t + 1]);
        if (c.is3d) {
            c.band[t] = (fcomplex*)alloc_off(sizeof(fcomplex) * kB * (size_t)c.nx * c.nz, 0, &c.base[3 * t + 2]);
        }
        if (c.in[t] == nullptr || c.out[t] == nullptr || (c.is3d && c.band[t] == nullptr)) {
            throw std::bad_alloc();
        }
    }
    // 1D plans as MATLAB makes them (rank 1, contiguous transforms, FFTW_ESTIMATE, out of place, the scratch
    // alignment of the buffers above): the process's cached plans, made under its one planner lock
    using pk5d::fft::cached_plan;
    for (int b = 0; b < 2; b++) {
        const int m = b == 0 ? kB : 1;
        c.fy[b] = cached_plan(be, 0, pk5d::fft::kR2C, c.ny, m);
        c.iy[b] = cached_plan(be, 0, pk5d::fft::kC2CB, c.ny, m);
        c.fx[b] = cached_plan(be, 0, pk5d::fft::kC2CF, c.nx, m);
        if (c.is3d) {
            c.ix[b] = cached_plan(be, 0, pk5d::fft::kC2CB, c.nx, m);
            c.fz[b] = cached_plan(be, 0, pk5d::fft::kC2CF, c.nz, m);
            c.il[b] = cached_plan(be, 0, pk5d::fft::kC2R, c.nz, m);
        } else {
            c.il[b] = cached_plan(be, 0, pk5d::fft::kC2R, c.nx, m);
        }
        if (c.fy[b] == nullptr || c.iy[b] == nullptr || c.fx[b] == nullptr || c.il[b] == nullptr ||
            (c.is3d && (c.ix[b] == nullptr || c.fz[b] == nullptr))) {
            throw std::runtime_error("FFTW planning failed");
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
inline void load_column(const Ctx& c, const Vol& V, const idx_t j, float* li) {
    const idx_t x = j % V.v[1], z = j / V.v[1];
    const int n = V.v[0];
    if (V.flip) {
        const float* src = V.p + ((idx_t)(V.v[1] - 1 - x) + (idx_t)V.v[1] * (V.v[2] - 1 - z)) * n;
        for (int y = 0; y < n; y++) li[y] = src[n - 1 - y];
    } else {
        std::memcpy(li, V.p + j * n, sizeof(float) * n);
    }
    for (int y = n; y < c.ny; y++) li[y] = 0.0f;
}

// forward stage 1: the y-half columns of the volume's data columns (x < v[1], z < v[2]; the others are zero)
// into Yh (column x + v[1] * z, ldr floats each)
void fwd_y(Ctx& c, const Vol& V, float* Yh) {
    const idx_t vx = V.v[1];
    const idx_t ncol = vx * V.v[2];
    const idx_t nb = (ncol + kB - 1) / kB;
    const FftwApi& F = *c.F;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = thread_num();
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
            F.exec_r2c(c.fy[0], in, out);
            for (int q = 0; q < ng; q++) {
                std::memcpy(Yh + (j0 + gl[q]) * c.ldr, out + (idx_t)q * c.ldr, sizeof(float) * c.ldr);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            load_column(c, V, j0 + l, in);
            F.exec_r2c(c.fy[1], in, out);
            std::memcpy(Yh + (j0 + l) * c.ldr, out, sizeof(float) * c.ldr);
        }
    }
}

// complex element of the full y range of a y-half column (rows above ny / 2 are conjugates)
inline void yfull(const float* col, const int ny, const int nyh, const int k, float& re, float& im) {
    if (k < nyh) {
        re = col[2 * k];
        im = col[2 * k + 1];
    } else {
        re = col[2 * (ny - k)];
        im = -col[2 * (ny - k) + 1];
    }
}

// line along x (row y, slice z < v[2]) of the full y range from the y-half columns, zero for x >= v[1]
inline void x_line(const Ctx& c, const float* Yh, const int* v, const int y, const int z, float* li) {
    for (int x = 0; x < v[1]; x++) {
        yfull(Yh + ((idx_t)x + (idx_t)v[1] * z) * c.ldr, c.ny, c.nyh, y, li[2 * x], li[2 * x + 1]);
    }
    for (int x = v[1]; x < c.nx; x++) {
        li[2 * x] = 0.0f;
        li[2 * x + 1] = 0.0f;
    }
}

// x transform of band row l at slice z into the band buffer [m][nx][nz] (z contiguous)
inline void to_band(float* band, const float* lo, const int l, const int z, const int nx, const int nz) {
    float* bl = band + 2 * ((idx_t)l * nx) * nz;
    for (int k = 0; k < nx; k++) {
        bl[2 * ((idx_t)k * nz + z)] = lo[2 * k];
        bl[2 * ((idx_t)k * nz + z) + 1] = lo[2 * k + 1];
    }
}

// MATLAB's complex product a .* t (no fused multiply-add) into r (r may be t)
inline void cmul(const float* a, const float* t, float* r) {
    const float re = a[0] * t[0] - a[1] * t[1];
    const float im = a[0] * t[1] + a[1] * t[0];
    r[0] = re;
    r[1] = im;
}

// Where the kept half spectrum of a forward transform goes: into P ([nh, ny, nx] for 3D, [ny, nh] for 2D), or
// for the second transform its product with P (A .* T, into P), except the self-conjugate planes (k_z = 0 and
// nz / 2, 2D: k_x = 0 and nx / 2), which wait in S (S[q], MATLAB's linear order) for MATLAB's Hermitian fill
struct Sink {
    fcomplex* P;
    float* S[2];
    bool product;
};

// the self-conjugate plane (0 or 1) of index k along the last dimension, or -1
inline int side_plane(const Ctx& c, const int k) {
    const int n = c.is3d ? c.nz : c.nx;
    if (k == 0) return 0;
    if (n % 2 == 0 && k == n / 2) return 1;
    return -1;
}

// the kept outputs k_z < nh of the z line of row y, column k
inline void z_out(const Ctx& c, const Sink& sk, const float* lo, const int y, const int k) {
    const idx_t L = (idx_t)y + (idx_t)c.ny * k;
    float* dst = (float*)sk.P + 2 * (idx_t)c.nh * L;
    if (!sk.product) {
        std::memcpy(dst, lo, sizeof(fcomplex) * c.nh);
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
void z_group(const Ctx& c, const float* band, const int* pl, const int* pk, const int np, float* in, float* out,
    const Sink& sk, const int y0) {
    const int nz = c.nz, nx = c.nx;
    for (int q = 0; q < kB; q++) {
        const int s = q < np ? q : np - 1;
        std::memcpy(in + 2 * (idx_t)q * nz, band + 2 * ((idx_t)pl[s] * nx + pk[s]) * nz, sizeof(fcomplex) * nz);
    }
    c.F->exec_dft(c.fz[0], in, out);
    for (int q = 0; q < np; q++) {
        z_out(c, sk, out + 2 * (idx_t)q * nz, y0 + pl[q], pk[q]);
    }
}

// forward stages 2 and 3 (3D): lines along x of the slices with data (z < v[2]), then all lines along z; the z
// outputs k_z < nh go to the sink. Rows are processed in bands of four; lines MATLAB transforms in groups are
// grouped in any fours (the lanes of a group are independent), the others transformed alone.
void fwd_xz(Ctx& c, const float* Yh, const int* v, const Sink& sk) {
    const int ny = c.ny, nx = c.nx, nz = c.nz;
    const int nbands = (ny + kB - 1) / kB;
    const FftwApi& F = *c.F;
    // the band entries of the zero slices (z >= v[2]) stay zero
    for (int t = 0; t < c.nt; t++) {
        std::memset(c.band[t], 0, sizeof(fcomplex) * kB * (size_t)nx * nz);
    }
#pragma omp parallel for num_threads(c.nt) schedule(dynamic, 1)
    for (int b = 0; b < nbands; b++) {
        const int t = thread_num();
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
                F.exec_dft(c.fx[0], in, out);
                for (int q = 0; q < ng; q++) {
                    to_band(band, out + 2 * (idx_t)q * nx, gl[q], z, nx, nz);
                }
            }
            for (int l = 0; l < m; l++) {
                if (!c.aloneX[L0 + l]) continue;
                x_line(c, Yh, v, y0 + l, z, in);
                F.exec_dft(c.fx[1], in, out);
                to_band(band, out, l, z, nx, nz);
            }
        }
        // lines along z (row y0 + l, column k)
        int pl[kB], pk[kB], np = 0;
        for (int k = 0; k < nx; k++) {
            for (int l = 0; l < m; l++) {
                if (c.aloneZ[(idx_t)(y0 + l) + (idx_t)ny * k]) {
                    std::memcpy(in, band + 2 * ((idx_t)l * nx + k) * nz, sizeof(fcomplex) * nz);
                    F.exec_dft(c.fz[1], in, out);
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
inline void x2d_out(const Ctx& c, const Sink& sk, const float* lo, const int y) {
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
void fwd_x2d(Ctx& c, const float* Yh, const int* v, const Sink& sk) {
    const int ny = c.ny, nx = c.nx;
    const int nb = (ny + kB - 1) / kB;
    const FftwApi& F = *c.F;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (int b = 0; b < nb; b++) {
        const int t = thread_num();
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
            F.exec_dft(c.fx[0], in, out);
            for (int q = 0; q < ng; q++) {
                x2d_out(c, sk, out + 2 * (idx_t)q * nx, y0 + gl[q]);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneX[y0 + l]) continue;
            x_line(c, Yh, v, y0 + l, 0, in);
            F.exec_dft(c.fx[1], in, out);
            x2d_out(c, sk, out, y0 + l);
        }
    }
}

// MATLAB's Hermitian fill of a self-conjugate plane (3D: ny x nx, 2D: a column of ny) whose element with
// MATLAB's linear index L in the plane is at base + 2 * st * L: an element whose mirror has the lower linear
// index becomes the mirror's conjugate, a self-conjugate element its real part
void fill_plane(const Ctx& c, float* base, const idx_t st) {
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
void fftn_half(Ctx& c, const Vol& V, float* Yh, const Sink& sk) {
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
inline idx_t yline_matlab(const Ctx& c, const idx_t i) { return c.is3d ? i / c.nh + (idx_t)c.nx * (i % c.nh) : i; }
inline idx_t yline_base(const Ctx& c, const idx_t i) {
    return c.is3d ? i % c.nh + (idx_t)c.nh * c.ny * (i / c.nh) : (idx_t)c.ny * i;
}

// y line i of P into li
inline void y_line(const Ctx& c, const fcomplex* P, const idx_t i, float* li) {
    const idx_t st = c.is3d ? c.nh : 1;
    const float* src = (const float*)P + 2 * yline_base(c, i);
    for (int k = 0; k < c.ny; k++) {
        li[2 * k] = src[2 * k * st];
        li[2 * k + 1] = src[2 * k * st + 1];
    }
}

// the inverse transform of y line i (lo), divided by ny, into P
inline void y_back(const Ctx& c, fcomplex* P, const float* lo, const idx_t i, const float fn) {
    const idx_t st = c.is3d ? c.nh : 1;
    float* pc = (float*)P + 2 * yline_base(c, i);
    for (int k = 0; k < c.ny; k++) {
        pc[2 * k * st] = lo[2 * k] / fn;
        pc[2 * k * st + 1] = lo[2 * k + 1] / fn;
    }
}

// inverse stage 1: P along y, c2c, divided by ny, in place. Only the kept lines (the first ones of MATLAB's
// nx * nz) are transformed.
void inv_y(Ctx& c, fcomplex* P) {
    const idx_t nl = c.is3d ? (idx_t)c.nx * c.nh : (idx_t)c.nh;
    const idx_t nb = (nl + kB - 1) / kB;
    const float fn = (float)c.ny;
    const FftwApi& F = *c.F;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = thread_num();
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
            F.exec_dft(c.iy[0], in, out);
            for (int q = 0; q < ng; q++) {
                y_back(c, P, out + 2 * (idx_t)q * c.ny, i0 + gl[q], fn);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            y_line(c, P, i0 + l, in);
            F.exec_dft(c.iy[1], in, out);
            y_back(c, P, out, i0 + l, fn);
        }
    }
}

// x line (window row r, plane kz) of P [nh, ny, nx] into li
inline void x_line_win(const Ctx& c, const fcomplex* P, const Win& W, const int r, const int kz, float* li) {
    const idx_t st = (idx_t)c.nh * c.ny;
    const float* src = (const float*)P + 2 * (kz + (idx_t)c.nh * (W.w0[0] + r));
    for (int x = 0; x < c.nx; x++) {
        li[2 * x] = src[2 * x * st];
        li[2 * x + 1] = src[2 * x * st + 1];
    }
}

// the window columns of the inverse x transform (lo) of window row r, plane kz, divided by nx, into Q [nh, wy, wx]
inline void x_back_win(const Ctx& c, const Win& W, const float* lo, const int r, const int kz, const float fn, float* Q) {
    const idx_t st = (idx_t)c.nh * W.wn[0];
    float* dst = Q + 2 * (kz + (idx_t)c.nh * r);
    for (int x = 0; x < W.wn[1]; x++) {
        dst[2 * x * st] = lo[2 * (W.w0[1] + x)] / fn;
        dst[2 * x * st + 1] = lo[2 * (W.w0[1] + x) + 1] / fn;
    }
}

// inverse stage 2 (3D): lines along x of the window rows in the kept z planes, c2c, divided by nx; the window
// columns go to Q [nh, wy, wx]
void inv_x_win(Ctx& c, const fcomplex* P, const Win& W, float* Q) {
    const int wy = W.wn[0];
    // lines numbered i = kz + nh * r (window row r); MATLAB numbers them y + ny * kz
    const idx_t nl = (idx_t)c.nh * wy;
    const idx_t nb = (nl + kB - 1) / kB;
    const float fn = (float)c.nx;
    const FftwApi& F = *c.F;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = thread_num();
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
            F.exec_dft(c.ix[0], in, out);
            for (int q = 0; q < ng; q++) {
                const idx_t i = i0 + gl[q];
                x_back_win(c, W, out + 2 * (idx_t)q * c.nx, (int)(i / c.nh), (int)(i % c.nh), fn, Q);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            const idx_t i = i0 + l;
            x_line_win(c, P, W, (int)(i / c.nh), (int)(i % c.nh), in);
            F.exec_dft(c.ix[1], in, out);
            x_back_win(c, W, out, (int)(i / c.nh), (int)(i % c.nh), fn, Q);
        }
    }
}

// the input of window line L of the inverse last stage (nh complex values): along z from Q [nh, wy, wx] (3D),
// along x from row w0[0] + L of P [ny, nh] (2D)
inline void last_line(const Ctx& c, const fcomplex* P, const float* Q, const Win& W, const idx_t L, float* li) {
    const int nh = c.nh;
    if (c.is3d) {
        std::memcpy(li, Q + 2 * (idx_t)nh * L, sizeof(fcomplex) * nh);
    } else {
        const float* src = (const float*)P + 2 * (W.w0[0] + L);
        for (int k = 0; k < nh; k++) {
            li[2 * k] = src[2 * (idx_t)k * c.ny];
            li[2 * k + 1] = src[2 * (idx_t)k * c.ny + 1];
        }
    }
}

// the window part of the real line L (lo) divided by its length, into R [wy, wx, wz]
inline void last_out(const Ctx& c, const Win& W, const float* lo, const idx_t L, float* R) {
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
void inv_last_win(Ctx& c, const fcomplex* P, const float* Q, const Win& W, float* R) {
    const int nh = c.nh;
    const int nlen = c.is3d ? c.nz : c.nx;
    const int wy = W.wn[0];
    // window lines: (y, x) pairs (3D) or rows (2D); MATLAB numbers them y + ny * x (3D) or y (2D)
    const idx_t nl = c.is3d ? (idx_t)wy * W.wn[1] : (idx_t)wy;
    const uint8_t* alone = c.is3d ? c.aloneZ.data() : c.aloneX.data();
    const idx_t nb = (nl + kB - 1) / kB;
    const FftwApi& F = *c.F;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = thread_num();
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
            F.exec_c2r(c.il[0], in, out);
            for (int q = 0; q < ng; q++) {
                last_out(c, W, out + (idx_t)q * nlen, L0 + gl[q], R);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!al[l]) continue;
            last_line(c, P, Q, W, L0 + l, in);
            F.exec_c2r(c.il[1], in, out);
            last_out(c, W, out, L0 + l, R);
        }
    }
}

// buffer of n elements, uninitialized, freed at scope exit
template <typename T>
struct Buf {
    T* p = nullptr;
    explicit Buf(const size_t n) {
        p = static_cast<T*>(std::malloc(sizeof(T) * (n > 0 ? n : 1)));
        if (p == nullptr) throw std::bad_alloc();
    }
    ~Buf() { std::free(p); }
    Buf(const Buf&) = delete;
    Buf& operator=(const Buf&) = delete;
};

// smallest n >= m whose prime factors are 2, 3, 5 and 7
int smooth_size(const int m) {
    for (int n = m;; n++) {
        int r = n;
        const int ps[4] = {2, 3, 5, 7};
        for (int k = 0; k < 4; k++) {
            while (r % ps[k] == 0) r /= ps[k];
        }
        if (r == 1) return n;
    }
}


// The window sums of normxcorr3_fast's MATLAB code in double (its box_sums): along each dimension d, for the
// window's positions k (k = 0 .. w[d] - 1), the sum of the elements lo[d][k] .. hi[d][k] (1-based, inclusive)
// of the previous stage, as MATLAB's cumsum from the first element (P(1) = x(1), P(i) = P(i - 1) + x(i)) and
// P(hi) - P(lo - 1) (P(hi) for lo = 1), along y, then x, then z, so that the sums are MATLAB's bit for bit.
// The elements are X, or X .* X in double (sq). n: size of X; out: [w0, w1, w2].
template <typename T>
void box_sums_window(const T* X, const int* n, int* const* lo, int* const* hi, const int* w, const bool sq,
    double* out, const int nt) {
    int H[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
        for (int k = 0; k < w[d]; k++) H[d] = hi[d][k] > H[d] ? hi[d][k] : H[d];
    }
    // stage 1 (y) of every column: S1 [w0, n1, n2]
    const idx_t ncol = (idx_t)n[1] * n[2];
    Buf<double> s1buf((size_t)w[0] * ncol);
    double* S1 = s1buf.p;
    {
        Buf<double> pbuf((size_t)H[0] * nt);
#pragma omp parallel num_threads(nt)
        {
            double* P = pbuf.p + (size_t)H[0] * thread_num();
#pragma omp for schedule(static)
            for (idx_t c = 0; c < ncol; c++) {
                const T* f = X + c * n[0];
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
    }
    // stage 2 (x) of the rows of every slice: S2 [w0, w1, n2]; P holds the prefix sums of the w0 rows
    const idx_t w0 = w[0];
    Buf<double> s2buf((size_t)w0 * w[1] * n[2]);
    double* S2 = s2buf.p;
    {
        Buf<double> pbuf((size_t)w0 * H[1] * nt);
#pragma omp parallel num_threads(nt)
        {
            double* P = pbuf.p + (size_t)w0 * H[1] * thread_num();
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
    }
    // stage 3 (z) of the lines of every window column: out [w0, w1, w2]
    const idx_t wyx = w0 * w[1];
    Buf<double> pbuf((size_t)w0 * H[2] * nt);
#pragma omp parallel num_threads(nt)
    {
        double* P = pbuf.p + (size_t)w0 * H[2] * thread_num();
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
}

}  // namespace

void window_corr(const float* T0, const int* ts, const float* A0, const int* a, const Win& W, const int nt, idx_t grain,
    const bool fast, const int backend, float* corr, double* tm) {
    int d[3];
    for (int i = 0; i < 3; i++) d[i] = ts[i] + a[i] - 1;
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
    double tl[4];
    double t0 = wtime();
    Ctx c;
    setup(c, d, nt, grain, backend);
    tl[0] = wtime() - t0;

    // the y-half columns, the half spectrum (of T0, then of the product), the self-conjugate planes of A0's
    const idx_t ncolMax = (idx_t)a[1] * a[2] > (idx_t)ts[1] * ts[2] ? (idx_t)a[1] * a[2] : (idx_t)ts[1] * ts[2];
    Buf<float> yh((size_t)c.ldr * ncolMax);
    Buf<fcomplex> pbuf((size_t)c.nP);
    const idx_t nplane = c.is3d ? (idx_t)c.ny * c.nx : (idx_t)c.ny;
    Buf<float> sbuf(2 * 2 * (size_t)nplane);
    float* Yh = yh.p;
    fcomplex* P = pbuf.p;
    float* S = sbuf.p;
    const Vol VT = {T0, {ts[0], ts[1], ts[2]}, true};
    const Vol VA = {A0, {a[0], a[1], a[2]}, false};
    t0 = wtime();
    const Sink skT = {P, {S, S + 2 * nplane}, false};
    fftn_half(c, VT, Yh, skT);
    tl[1] = wtime() - t0;
    t0 = wtime();
    const Sink skA = {P, {S, S + 2 * nplane}, true};
    fftn_half(c, VA, Yh, skA);
    tl[2] = wtime() - t0;
    t0 = wtime();
    inv_y(c, P);
    if (c.is3d) {
        Buf<float> qbuf(2 * (size_t)W.wn[0] * W.wn[1] * c.nh);
        inv_x_win(c, P, W, qbuf.p);
        inv_last_win(c, P, qbuf.p, W, corr);
    } else {
        inv_last_win(c, P, nullptr, W, corr);
    }
    tl[3] = wtime() - t0;
    if (tm != nullptr) {
        for (int i = 0; i < 4; i++) tm[i] = tl[i];
    }
}

template <typename T>
void window_sums(const T* T0, const int* ts, const T* A0, const int* a, const Win& W, double* sA, double* sA2,
    double* sT, const int nt) {
    // at window position k (1-based in the full output) the image's elements max(1, k - ts + 1) .. min(k, a)
    // and the template's max(1, ts - k + 1) .. min(ts, ts - k + a) overlap
    std::vector<int> lohi(4 * (size_t)(W.wn[0] + W.wn[1] + W.wn[2]));
    int* loA[3];
    int* hiA[3];
    int* loT[3];
    int* hiT[3];
    int* q = lohi.data();
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
    if (sA != nullptr) box_sums_window(A0, a, loA, hiA, W.wn, false, sA, nt);
    if (sA2 != nullptr) box_sums_window(A0, a, loA, hiA, W.wn, true, sA2, nt);
    if (sT != nullptr) box_sums_window(T0, ts, loT, hiT, W.wn, false, sT, nt);
}

template void window_sums<float>(const float*, const int*, const float*, const int*, const Win&, double*, double*,
    double*, int);
template void window_sums<double>(const double*, const int*, const double*, const int*, const Win&, double*, double*,
    double*, int);

}  // namespace pk5d::xcorr
