// decon_lucy_rfft_mex's Richardson-Lucy steps in MATLAB's layout (the C array (z, x, y)), with the FFTs of
// MATLAB's fftn and ifftn (on the exact backend MATLAB's bits) through a backend's FFT table
// (cpp/fft/backend.h; no FFT library is included here):
//   fftn of real data: 1D real-to-complex transforms along y, complex ones along x and z, then MATLAB's
//     Hermitian fill (an element whose mirror has the lower linear index is the conjugate of the mirror,
//     self-conjugate elements are real)
//   ifftn of the (Hermitian) product: complex transforms along y and x, then complex-to-real along z, each
//     divided by its length
//   MATLAB runs the 1D transforms of a stage in groups of four lines and some lines one by one (mark_alone);
//     for some lengths (e.g., 64 or 128) the two give different last bits, so the same lines run alone here.
// Only the half spectrum (k_z <= nz / 2) the inverse uses is kept and transformed back, the element-wise
// steps run inside the transforms' copy loops, and no complex volume of the full size is allocated. 2D data
// uses the same steps along y and x. The loops are the mex's, statement for statement.
#include <omp.h>

#include <cfloat>
#include <cmath>
#include <cstring>
#include <mutex>
#include <stdexcept>

#include "backend.h"
#include "decon.hpp"

namespace pk5d::decon {
namespace {

using pk5d::fft::FftwApi;
using pk5d::fft::Iodim;

// transforms per execute, as in MATLAB's fft
constexpr int kB = 4;
// byte offsets (modulo 64) of the scratch buffers, as MATLAB's
constexpr int kInOff = 16;
constexpr int kOutOff = 32;
constexpr unsigned kEstimate = 1U << 6;  // the planner flag MATLAB uses
constexpr int kForward = -1;
constexpr int kBackward = 1;

struct Timer {
    Times& t;
    double t0 = 0;
    explicit Timer(Times& tt) : t(tt) {}
    void start() { t0 = omp_get_wtime(); }
    void stop(int k) { t.t[k] += omp_get_wtime() - t0; }
};

// 64-byte aligned storage plus an offset
struct Scratch {
    void* base = nullptr;
    float* p = nullptr;
    Scratch() = default;
    Scratch(const Scratch&) = delete;
    Scratch& operator=(const Scratch&) = delete;
    ~Scratch() { std::free(base); }
    void alloc(size_t bytes, int off) {
        base = std::malloc(bytes + 128);
        if (base == nullptr) throw std::bad_alloc();
        p = reinterpret_cast<float*>(((reinterpret_cast<uintptr_t>(base) + 63) & ~static_cast<uintptr_t>(63)) + off);
    }
};

}  // namespace

// sizes, plans, the lines MATLAB transforms alone, and per-thread scratch of one call
struct RfftPlan {
    int ny = 1, nx = 1, nz = 1;  // data size (nz = 1 for 2D)
    bool is3d = false;
    int nyh = 1;   // ny / 2 + 1, complex elements of the y-half of a column
    int nh = 1;    // kept elements along the last dimension: nz / 2 + 1 (3D), nx / 2 + 1 (2D)
    idx_t ldr = 2; // floats per column of the real / y-half buffer (2 * nyh)
    idx_t nP = 0;  // complex elements of a half spectrum
    int nt = 1;
    const FftwApi* api = nullptr;
    // plans [0]: groups of four, [1]: single transforms
    void* fy[2] = {nullptr, nullptr};
    void* fx[2] = {nullptr, nullptr};
    void* fz[2] = {nullptr, nullptr};
    void* iy[2] = {nullptr, nullptr};
    void* ix[2] = {nullptr, nullptr};
    void* il[2] = {nullptr, nullptr};
    // nonzero for the lines MATLAB transforms alone: along y (column x + nx * z), along x (row y + ny * z),
    // along z (y + ny * x, 3D only)
    std::vector<uint8_t> aloneY, aloneX, aloneZ;
    // scratch per thread: in, out and the x-line band buffer of the 3D forward
    std::vector<Scratch> sin, sout, sband;
    RfftPlan() = default;
    RfftPlan(const RfftPlan&) = delete;
    RfftPlan& operator=(const RfftPlan&) = delete;
    ~RfftPlan() {
        if (api == nullptr) return;
        std::lock_guard<std::mutex> g(pk5d::fft::planner_lock());
        void** ps[6] = {fy, fx, fz, iy, ix, il};
        for (auto* q : ps) {
            for (int i = 0; i < 2; i++) {
                if (q[i] != nullptr) api->destroy(q[i]);
            }
        }
    }
    float* in(int t) const { return sin[t].p; }
    float* out(int t) const { return sout[t].p; }
    float* band(int t) const { return sband[t].p; }
};

namespace {

using Ctx = RfftPlan;

// Marks the lines a .. b - 1 of a stage that MATLAB transforms alone. MATLAB's FFT splits the lines of a
// stage among its threads by halving them (at a + (b - a) / 2) until a piece has at most grain lines (when
// it runs in parallel and the stage is large enough, see stage_grain; grain 0: one piece), and runs a piece
// in groups of four lines from its start, the last (length % 4) lines one by one. Lines along x of 3D data
// (seg = ny) also restart the groups at each z slice.
void mark_alone(uint8_t* alone, idx_t a, idx_t b, idx_t grain, idx_t seg) {
    if (grain > 0 && b - a > grain) {
        const idx_t mid = a + (b - a) / 2;
        mark_alone(alone, a, mid, grain, seg);
        mark_alone(alone, mid, b, grain, seg);
        return;
    }
    idx_t s = a;
    while (s < b) {
        idx_t e = b;
        if (seg > 0 && (s / seg + 1) * seg < b) e = (s / seg + 1) * seg;
        for (idx_t L = e - (e - s) % kB; L < e; L++) alone[L] = 1;
        s = e;
    }
}

// MATLAB splits the N lines of length n of a stage only if n * N * (bits of n) > 98304 (its estimate of the
// work of the stage), otherwise one thread transforms them
idx_t stage_grain(idx_t grain, int n, idx_t N) {
    int bits = 0;
    for (int v = n; v > 0; v >>= 1) bits++;
    return (grain > 0 && static_cast<double>(n) * static_cast<double>(N) * bits > 98304.0) ? grain : 0;
}

void setup(Ctx& c, const Call& call) {
    c.ny = call.d[0];
    c.nx = call.d[1];
    c.nz = call.d[2];
    c.is3d = call.d[2] > 1;
    c.nyh = c.ny / 2 + 1;
    c.nh = c.is3d ? c.nz / 2 + 1 : c.nx / 2 + 1;
    c.ldr = 2 * static_cast<idx_t>(c.nyh);
    c.nP = c.is3d ? static_cast<idx_t>(c.ny) * c.nx * c.nh : static_cast<idx_t>(c.ny) * c.nh;
    c.nt = call.nt > 0 ? call.nt : 1;
    const idx_t grain = call.grain;
    // the lines MATLAB transforms alone, per stage (inverse stages as the forward ones of the full spectrum)
    const idx_t nlY = static_cast<idx_t>(c.nx) * c.nz, nlX = static_cast<idx_t>(c.ny) * c.nz;
    c.aloneY.assign(static_cast<size_t>(nlY), 0);
    c.aloneX.assign(static_cast<size_t>(nlX), 0);
    mark_alone(c.aloneY.data(), 0, nlY, stage_grain(grain, c.ny, nlY), 0);
    mark_alone(c.aloneX.data(), 0, nlX, stage_grain(grain, c.nx, nlX), c.is3d ? c.ny : 0);
    if (c.is3d) {
        const idx_t nlZ = static_cast<idx_t>(c.ny) * c.nx;
        c.aloneZ.assign(static_cast<size_t>(nlZ), 0);
        mark_alone(c.aloneZ.data(), 0, nlZ, stage_grain(grain, c.nz, nlZ), 0);
    }
    const int maxn = c.ny > c.nx ? (c.ny > c.nz ? c.ny : c.nz) : (c.nx > c.nz ? c.nx : c.nz);
    const size_t sbytes = 2 * sizeof(float) * kB * static_cast<size_t>(maxn);
    c.sin = std::vector<Scratch>(static_cast<size_t>(c.nt));
    c.sout = std::vector<Scratch>(static_cast<size_t>(c.nt));
    c.sband = std::vector<Scratch>(static_cast<size_t>(c.nt));
    for (int t = 0; t < c.nt; t++) {
        c.sin[t].alloc(sbytes, kInOff);
        c.sout[t].alloc(sbytes, kOutOff);
        if (c.is3d) c.sband[t].alloc(2 * sizeof(float) * kB * static_cast<size_t>(c.nx) * c.nz, 0);
    }
    // 1D plans as MATLAB makes them: rank 1, contiguous transforms, its planner flag
    const FftwApi& a = pk5d::fft::backend(call.backend).api[0];
    std::lock_guard<std::mutex> g(pk5d::fft::planner_lock());
    c.api = &a;
    float* ri = c.in(0);
    float* ro = c.out(0);
    for (int b = 0; b < 2; b++) {
        const idx_t m = b == 0 ? kB : 1;
        const Iodim dy = {c.ny, 1, 1}, hy = {m, c.ny, c.nyh}, hyc = {m, c.ny, c.ny};
        c.fy[b] = a.plan_r2c(1, &dy, 1, &hy, ri, ro, kEstimate);
        c.iy[b] = a.plan_dft(1, &dy, 1, &hyc, ri, ro, kBackward, kEstimate);
        const Iodim dx = {c.nx, 1, 1}, hx = {m, c.nx, c.nx};
        c.fx[b] = a.plan_dft(1, &dx, 1, &hx, ri, ro, kForward, kEstimate);
        if (c.is3d) {
            c.ix[b] = a.plan_dft(1, &dx, 1, &hx, ri, ro, kBackward, kEstimate);
            const Iodim dz = {c.nz, 1, 1}, hz = {m, c.nz, c.nz}, hzr = {m, c.nh, c.nz};
            c.fz[b] = a.plan_dft(1, &dz, 1, &hz, ri, ro, kForward, kEstimate);
            c.il[b] = a.plan_c2r(1, &dz, 1, &hzr, ri, ro, kEstimate);
        } else {
            const Iodim hxr = {m, c.nh, c.nx};
            c.il[b] = a.plan_c2r(1, &dx, 1, &hxr, ri, ro, kEstimate);
        }
        if (c.fy[b] == nullptr || c.iy[b] == nullptr || c.fx[b] == nullptr || c.il[b] == nullptr ||
            (c.is3d && (c.ix[b] == nullptr || c.fz[b] == nullptr))) {
            throw std::runtime_error("decon: FFT planning failed");
        }
    }
}

inline void exec_r2c(const Ctx& c, void* p, float* in, float* out) { c.api->exec_r2c(p, in, out); }
inline void exec_dft(const Ctx& c, void* p, float* in, float* out) { c.api->exec_dft(p, in, out); }
inline void exec_c2r(const Ctx& c, void* p, float* in, float* out) { c.api->exec_c2r(p, in, out); }

// lanes of a block of m lines (first line j0) that MATLAB transforms in a group: their offsets in gl, the
// count returned
inline int grouped_lanes(const uint8_t* alone, idx_t j0, int m, int* gl) {
    int ng = 0;
    for (int l = 0; l < m; l++) {
        if (!alone[j0 + l]) gl[ng++] = l;
    }
    return ng;
}

// complex element of the full y range of a y-half column (rows above ny / 2 are conjugates)
inline void yfull(const float* col, int ny, int nyh, int k, float& re, float& im) {
    if (k < nyh) {
        re = col[2 * k];
        im = col[2 * k + 1];
    } else {
        re = col[2 * (ny - k)];
        im = -col[2 * (ny - k) + 1];
    }
}

// line along x (row y, slice z) of the full y range from the y-half columns of A
inline void x_line(const Ctx& c, const float* A, int y, int z, float* li) {
    for (int x = 0; x < c.nx; x++) {
        yfull(A + (static_cast<idx_t>(x) + static_cast<idx_t>(c.nx) * z) * c.ldr, c.ny, c.nyh, y, li[2 * x],
              li[2 * x + 1]);
    }
}

// forward stage 1: real columns (stride lds floats; nx * nz of them) to y-half columns of A (stride ldr
// floats). A may be src (lds == ldr).
void fwd_y(Ctx& c, const float* src, idx_t lds, float* A) {
    const idx_t ncol = static_cast<idx_t>(c.nx) * c.nz;
    const idx_t nb = (ncol + kB - 1) / kB;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const idx_t j0 = b * kB;
        const int m = static_cast<int>(ncol - j0 < kB ? ncol - j0 : kB);
        float* in = c.in(t);
        float* out = c.out(t);
        int gl[kB];
        const int ng = grouped_lanes(c.aloneY.data(), j0, m, gl);
        if (ng > 0) {
            // all inputs are read before A is written (A may be src)
            for (int q = 0; q < kB; q++) {
                std::memcpy(in + static_cast<idx_t>(q) * c.ny, src + (j0 + gl[q < ng ? q : ng - 1]) * lds,
                            sizeof(float) * c.ny);
            }
            exec_r2c(c, c.fy[0], in, out);
            for (int q = 0; q < ng; q++) {
                std::memcpy(A + (j0 + gl[q]) * c.ldr, out + static_cast<idx_t>(q) * c.ldr, sizeof(float) * c.ldr);
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneY[j0 + l]) continue;
            std::memcpy(in, src + (j0 + l) * lds, sizeof(float) * c.ny);
            exec_r2c(c, c.fy[1], in, out);
            std::memcpy(A + (j0 + l) * c.ldr, out, sizeof(float) * c.ldr);
        }
    }
}

// x transform of band row l at slice z into the band buffer [m][nx][nz] (z contiguous)
inline void to_band(float* band, const float* lo, int l, int z, int nx, int nz) {
    float* bl = band + 2 * (static_cast<idx_t>(l) * nx) * nz;
    for (int k = 0; k < nx; k++) {
        bl[2 * (static_cast<idx_t>(k) * nz + z)] = lo[2 * k];
        bl[2 * (static_cast<idx_t>(k) * nz + z) + 1] = lo[2 * k + 1];
    }
}

// the kept outputs k_z < nh of the z line of row y, column k to P [ny, nx, nh]
inline void z_out(const Ctx& c, cfloat* P, const float* lo, int y, int k) {
    const idx_t nl = static_cast<idx_t>(c.ny) * c.nx;
    float* dst = reinterpret_cast<float*>(P) + 2 * (static_cast<idx_t>(y) + static_cast<idx_t>(c.ny) * k);
    for (int kz = 0; kz < c.nh; kz++) {
        dst[2 * kz * nl] = lo[2 * kz];
        dst[2 * kz * nl + 1] = lo[2 * kz + 1];
    }
}

// a group of np z lines (band rows pl, columns pk), padded with its last line
void z_group(const Ctx& c, const float* band, const int* pl, const int* pk, int np, float* in, float* out, cfloat* P,
             int y0) {
    const int nz = c.nz, nx = c.nx;
    for (int q = 0; q < kB; q++) {
        const int s = q < np ? q : np - 1;
        std::memcpy(in + 2 * static_cast<idx_t>(q) * nz, band + 2 * (static_cast<idx_t>(pl[s]) * nx + pk[s]) * nz,
                    2 * sizeof(float) * nz);
    }
    exec_dft(c, c.fz[0], in, out);
    for (int q = 0; q < np; q++) z_out(c, P, out + 2 * static_cast<idx_t>(q) * nz, y0 + pl[q], pk[q]);
}

// forward stages 2 and 3 (3D): lines along x, then along z; the z outputs k_z < nh go to P [ny, nx, nh].
// Rows are processed in bands of four; lines MATLAB transforms in groups are grouped in any fours (the lanes
// of a group are independent), the others transformed alone.
void fwd_xz(Ctx& c, const float* A, cfloat* P) {
    const int ny = c.ny, nx = c.nx, nz = c.nz;
    const int nbands = (ny + kB - 1) / kB;
#pragma omp parallel for num_threads(c.nt) schedule(dynamic, 1)
    for (int b = 0; b < nbands; b++) {
        const int t = omp_get_thread_num();
        const int y0 = b * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = c.in(t);
        float* out = c.out(t);
        float* band = c.band(t);
        // lines along x of the band rows, every z slice
        for (int z = 0; z < nz; z++) {
            const idx_t L0 = static_cast<idx_t>(y0) + static_cast<idx_t>(ny) * z;
            int gl[kB];
            const int ng = grouped_lanes(c.aloneX.data(), L0, m, gl);
            if (ng > 0) {
                for (int q = 0; q < kB; q++) {
                    x_line(c, A, y0 + gl[q < ng ? q : ng - 1], z, in + 2 * static_cast<idx_t>(q) * nx);
                }
                exec_dft(c, c.fx[0], in, out);
                for (int q = 0; q < ng; q++) to_band(band, out + 2 * static_cast<idx_t>(q) * nx, gl[q], z, nx, nz);
            }
            for (int l = 0; l < m; l++) {
                if (!c.aloneX[L0 + l]) continue;
                x_line(c, A, y0 + l, z, in);
                exec_dft(c, c.fx[1], in, out);
                to_band(band, out, l, z, nx, nz);
            }
        }
        // lines along z (row y0 + l, column k)
        int pl[kB], pk[kB], np = 0;
        for (int k = 0; k < nx; k++) {
            for (int l = 0; l < m; l++) {
                if (c.aloneZ[static_cast<idx_t>(y0 + l) + static_cast<idx_t>(ny) * k]) {
                    std::memcpy(in, band + 2 * (static_cast<idx_t>(l) * nx + k) * nz, 2 * sizeof(float) * nz);
                    exec_dft(c, c.fz[1], in, out);
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
        if (np > 0) z_group(c, band, pl, pk, np, in, out, P, y0);
    }
}

// the kept outputs k_x < nh of the x transform of row y (2D) to P [ny, nh]
inline void x2d_out(const Ctx& c, cfloat* P, const float* lo, int y) {
    float* dst = reinterpret_cast<float*>(P) + 2 * static_cast<idx_t>(y);
    for (int k = 0; k < c.nh; k++) {
        dst[2 * static_cast<idx_t>(k) * c.ny] = lo[2 * k];
        dst[2 * static_cast<idx_t>(k) * c.ny + 1] = lo[2 * k + 1];
    }
}

// forward stage 2 (2D): lines along x (rows)
void fwd_x2d(Ctx& c, const float* A, cfloat* P) {
    const int ny = c.ny, nx = c.nx;
    const int nb = (ny + kB - 1) / kB;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (int b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const int y0 = b * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = c.in(t);
        float* out = c.out(t);
        int gl[kB];
        const int ng = grouped_lanes(c.aloneX.data(), y0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) x_line(c, A, y0 + gl[q < ng ? q : ng - 1], 0, in + 2 * static_cast<idx_t>(q) * nx);
            exec_dft(c, c.fx[0], in, out);
            for (int q = 0; q < ng; q++) x2d_out(c, P, out + 2 * static_cast<idx_t>(q) * nx, y0 + gl[q]);
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneX[y0 + l]) continue;
            x_line(c, A, y0 + l, 0, in);
            exec_dft(c, c.fx[1], in, out);
            x2d_out(c, P, out, y0 + l);
        }
    }
}

// MATLAB's Hermitian fill on the kept half: in the self-conjugate planes (k_z = 0 and nz / 2 for 3D, k_x = 0
// and nx / 2 for 2D) an element whose mirror has the lower linear index becomes the mirror's conjugate, a
// self-conjugate element its real part
void hermitian_fill(Ctx& c, cfloat* P) {
    const int ny = c.ny, nx = c.nx;
    float* p = reinterpret_cast<float*>(P);
    if (c.is3d) {
        const int planes[2] = {0, c.nz / 2};
        const int np = (c.nz % 2 == 0) ? 2 : 1;
        for (int q = 0; q < np; q++) {
            float* pl = p + 2 * static_cast<idx_t>(planes[q]) * ny * nx;
#pragma omp parallel for num_threads(c.nt) schedule(static)
            for (int k2 = 0; k2 < nx; k2++) {
                const int m2 = (nx - k2) % nx;
                for (int k1 = 0; k1 < ny; k1++) {
                    const int m1 = (ny - k1) % ny;
                    const idx_t L = k1 + static_cast<idx_t>(ny) * k2, Lm = m1 + static_cast<idx_t>(ny) * m2;
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
            float* cl = p + 2 * static_cast<idx_t>(cols[q]) * ny;
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
inline void product_column(const cfloat* P, const cfloat* H, idx_t j, int ny, bool conjugate, float* li) {
    const float* pc = reinterpret_cast<const float*>(P) + 2 * j * ny;
    const float* hc = reinterpret_cast<const float*>(H) + 2 * j * ny;
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
void inv_y(Ctx& c, cfloat* P, const cfloat* H, bool conjugate) {
    const int ny = c.ny;
    const idx_t ncol = c.is3d ? static_cast<idx_t>(c.nx) * c.nh : static_cast<idx_t>(c.nh);
    const idx_t nb = (ncol + kB - 1) / kB;
    const float fn = static_cast<float>(ny);
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        float* in = c.in(t);
        float* out = c.out(t);
        const idx_t j0 = b * kB;
        const int m = static_cast<int>(ncol - j0 < kB ? ncol - j0 : kB);
        int gl[kB];
        const int ng = grouped_lanes(c.aloneY.data(), j0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                product_column(P, H, j0 + gl[q < ng ? q : ng - 1], ny, conjugate, in + 2 * static_cast<idx_t>(q) * ny);
            }
            exec_dft(c, c.iy[0], in, out);
            for (int q = 0; q < ng; q++) {
                float* pc = reinterpret_cast<float*>(P) + 2 * (j0 + gl[q]) * ny;
                const float* lo = out + 2 * static_cast<idx_t>(q) * ny;
                for (int k = 0; k < 2 * ny; k++) pc[k] = lo[k] / fn;
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneY[j0 + l]) continue;
            product_column(P, H, j0 + l, ny, conjugate, in);
            exec_dft(c, c.iy[1], in, out);
            float* pc = reinterpret_cast<float*>(P) + 2 * (j0 + l) * ny;
            for (int k = 0; k < 2 * ny; k++) pc[k] = out[k] / fn;
        }
    }
}

// inverse stage 2 (3D): lines along x of the kept z planes, c2c, divided by nx
void inv_x(Ctx& c, cfloat* P) {
    const int ny = c.ny, nx = c.nx;
    const int nb = (ny + kB - 1) / kB;
    const idx_t nl = static_cast<idx_t>(ny) * nx;
    const float fn = static_cast<float>(nx);
    const idx_t ntask = static_cast<idx_t>(c.nh) * nb;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t task = 0; task < ntask; task++) {
        const int t = omp_get_thread_num();
        const int kz = static_cast<int>(task / nb);
        const int y0 = static_cast<int>(task % nb) * kB;
        const int m = ny - y0 < kB ? ny - y0 : kB;
        float* in = c.in(t);
        float* out = c.out(t);
        float* pl = reinterpret_cast<float*>(P) + 2 * static_cast<idx_t>(kz) * nl;
        const idx_t L0 = static_cast<idx_t>(y0) + static_cast<idx_t>(ny) * kz;
        int gl[kB];
        const int ng = grouped_lanes(c.aloneX.data(), L0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) {
                const int y = y0 + gl[q < ng ? q : ng - 1];
                float* li = in + 2 * static_cast<idx_t>(q) * nx;
                for (int x = 0; x < nx; x++) {
                    li[2 * x] = pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x)];
                    li[2 * x + 1] = pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x) + 1];
                }
            }
            exec_dft(c, c.ix[0], in, out);
            for (int q = 0; q < ng; q++) {
                const int y = y0 + gl[q];
                const float* lo = out + 2 * static_cast<idx_t>(q) * nx;
                for (int x = 0; x < nx; x++) {
                    pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x)] = lo[2 * x] / fn;
                    pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x) + 1] = lo[2 * x + 1] / fn;
                }
            }
        }
        for (int l = 0; l < m; l++) {
            if (!c.aloneX[L0 + l]) continue;
            const int y = y0 + l;
            for (int x = 0; x < nx; x++) {
                in[2 * x] = pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x)];
                in[2 * x + 1] = pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x) + 1];
            }
            exec_dft(c, c.ix[1], in, out);
            for (int x = 0; x < nx; x++) {
                pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x)] = out[2 * x] / fn;
                pl[2 * (static_cast<idx_t>(y) + static_cast<idx_t>(ny) * x) + 1] = out[2 * x + 1] / fn;
            }
        }
    }
}

// element-wise step on the real result v of the inverse at voxel i (MATLAB's linear index of the data)
template <typename T>
struct Ratio {
    // dst = max(I - bg, 0) ./ max(v, eps), dst with column stride ldd
    float* dst;
    idx_t ldd;
    const T* raw;
    double bg;
    inline void operator()(idx_t i, idx_t col, int y, float v) const {
        const float epsf = static_cast<float>(DBL_EPSILON);
        const float cx = v > epsf ? v : epsf;
        float in = static_cast<float>(static_cast<double>(raw[i]) - bg);
        in = in > 0.0f ? in : 0.0f;
        dst[col * ldd + y] = in / cx;
    }
};

struct Update {
    // out = max(v .* J, 0)
    float* out;
    const float* j;
    inline void operator()(idx_t i, idx_t, int, float v) const {
        const float w = v * j[i];
        out[i] = w > 0.0f ? w : 0.0f;
    }
};

// the kept half spectrum of line L of the last dimension (k = 0 .. nh - 1, stride nl) into li
inline void last_line(const cfloat* P, idx_t L, idx_t nl, int nh, float* li) {
    const float* src = reinterpret_cast<const float*>(P) + 2 * L;
    for (int k = 0; k < nh; k++) {
        li[2 * k] = src[2 * k * nl];
        li[2 * k + 1] = src[2 * k * nl + 1];
    }
}

// the real line L (lo, nlen values) divided by its length, then the element-wise step
template <typename Op>
inline void last_out(const Ctx& c, const Op& op, const float* lo, idx_t L, int nlen, idx_t nl, float fn) {
    const int y = static_cast<int>(L % c.ny);
    const idx_t col0 = L / c.ny;
    const idx_t colstep = c.is3d ? static_cast<idx_t>(c.nx) : 1;
    for (int k = 0; k < nlen; k++) op(L + static_cast<idx_t>(k) * nl, col0 + static_cast<idx_t>(k) * colstep, y, lo[k] / fn);
}

// inverse last stage: c2r along z (2D: x) of the kept half, divided by its length, then the element-wise step
template <typename Op>
void inv_last(Ctx& c, const cfloat* P, const Op& op) {
    const int nh = c.nh;
    const int nlen = c.is3d ? c.nz : c.nx;                                       // length of the real lines
    const idx_t nl = c.is3d ? static_cast<idx_t>(c.ny) * c.nx : static_cast<idx_t>(c.ny);  // number of lines
    const uint8_t* alone = c.is3d ? c.aloneZ.data() : c.aloneX.data();
    const idx_t nb = (nl + kB - 1) / kB;
    const float fn = static_cast<float>(nlen);
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const int t = omp_get_thread_num();
        const idx_t L0 = b * kB;
        const int m = static_cast<int>(nl - L0 < kB ? nl - L0 : kB);
        float* in = c.in(t);
        float* out = c.out(t);
        int gl[kB];
        const int ng = grouped_lanes(alone, L0, m, gl);
        if (ng > 0) {
            for (int q = 0; q < kB; q++) last_line(P, L0 + gl[q < ng ? q : ng - 1], nl, nh, in + 2 * static_cast<idx_t>(q) * nh);
            exec_c2r(c, c.il[0], in, out);
            for (int q = 0; q < ng; q++) last_out(c, op, out + static_cast<idx_t>(q) * nlen, L0 + gl[q], nlen, nl, fn);
        }
        for (int l = 0; l < m; l++) {
            if (!alone[L0 + l]) continue;
            last_line(P, L0 + l, nl, nh, in);
            exec_c2r(c, c.il[1], in, out);
            last_out(c, op, out, L0 + l, nlen, nl, fn);
        }
    }
}

// fftn of a real volume (columns of stride lds floats) into the kept half P, with MATLAB's Hermitian fill;
// A: y-half buffer (may be src when lds == ldr)
void fftn_half(Ctx& c, const float* src, idx_t lds, float* A, cfloat* P, Timer& tm) {
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
void ifftn_op(Ctx& c, cfloat* P, const cfloat* H, bool conjugate, const Op& op, Timer& tm) {
    tm.start();
    inv_y(c, P, H, conjugate);
    if (c.is3d) inv_x(c, P);
    inv_last(c, P, op);
    tm.stop(2);
}

template <typename T>
inline float input_value(const T* raw, idx_t i, double bg) {
    const float v = static_cast<float>(static_cast<double>(raw[i]) - bg);
    return v > 0.0f ? v : 0.0f;
}

// the ratio buffer (y-half columns, at least n floats); it also takes the result
float* ratio_buffer(const Ctx& c, Owned<float>& r) { return r.alloc(static_cast<size_t>(c.ldr) * c.nx * c.nz); }

// a copy in MATLAB's layout of a (z, y, x) array, or the array itself
template <typename T>
const T* matlab_input(const Ctx& c, const T* a, bool native_io, Owned<T>& copy) {
    if (!native_io || a == nullptr) return a;
    T* d = copy.alloc(static_cast<size_t>(c.ny) * c.nx * c.nz);
    transpose_planes(a, d, c.nz, c.ny, c.nx, sizeof(T), c.nt);
    return d;
}

}  // namespace

void RfftPlanDeleter::operator()(RfftPlan* p) const { delete p; }

RfftPlanPtr matlab_plan(const Call& call, Times& times) {
    Timer tm(times);
    tm.start();
    RfftPlanPtr c(new RfftPlan);
    setup(*c, call);
    tm.stop(0);
    return c;
}

void matlab_otf(RfftPlan& c, const float* psf, const int* ps, cfloat* h, Times& times) {
    Timer tm(times);
    Owned<float> rb(static_cast<size_t>(c.ldr) * c.nx * c.nz);
    float* r = rb.p;
    tm.start();
    const idx_t nr = c.ldr * c.nx * c.nz;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (idx_t i = 0; i < nr; i++) r[i] = 0.0f;
    // the PSF zero-padded to the data size and circularly shifted by -floor(psfSize / 2); ps = {pz, py, px}
    const int pz = ps[0], py = ps[1], px = ps[2];
    const int cy = py / 2, cx = px / 2, cz = pz / 2;
#pragma omp parallel for num_threads(c.nt) schedule(static)
    for (int z = 0; z < pz; z++) {
        const idx_t oz = static_cast<idx_t>((z - cz + c.nz) % c.nz) * c.nx;
        for (int yy = 0; yy < py; yy++) {
            const float* pp = psf + (static_cast<idx_t>(z) * py + yy) * px;
            const idx_t row = (yy - cy + c.ny) % c.ny;
            for (int x = 0; x < px; x++) r[(oz + (x - cx + c.nx) % c.nx) * c.ldr + row] = pp[x];
        }
    }
    tm.stop(3);
    fftn_half(c, r, c.ldr, r, h, tm);
}

template <typename T>
float* matlab_omw(RfftPlan& c, const T* raw_in, double bg, const cfloat* hf, const cfloat* hb, int n_iter,
                  const float* j0_in, bool native_io, Times& times) {
    Timer tm(times);
    const idx_t n = static_cast<idx_t>(c.ny) * c.nx * c.nz;
    Owned<float> rb;
    float* r = ratio_buffer(c, rb);
    if (n_iter == 0) {
        // element by element: the result in the layout of the inputs
        tm.start();
#pragma omp parallel for num_threads(c.nt) schedule(static)
        for (idx_t i = 0; i < n; i++) r[i] = j0_in != nullptr ? j0_in[i] : input_value(raw_in, i, bg);
        tm.stop(3);
        return rb.release();
    }
    tm.start();
    Owned<T> rawc;
    Owned<float> jb;
    const T* raw = matlab_input(c, raw_in, native_io, rawc);
    float* j = jb.alloc(static_cast<size_t>(n));
    if (j0_in != nullptr && native_io) {
        transpose_planes(j0_in, j, c.nz, c.ny, c.nx, sizeof(float), c.nt);
    } else {
#pragma omp parallel for num_threads(c.nt) schedule(static)
        for (idx_t i = 0; i < n; i++) j[i] = j0_in != nullptr ? j0_in[i] : input_value(raw, i, bg);
    }
    Owned<cfloat> pb(static_cast<size_t>(c.nP));
    cfloat* p = pb.p;
    tm.stop(3);
    const Ratio<T> ratio = {r, c.ldr, raw, bg};
    for (int k = 1; k <= n_iter; k++) {
        fftn_half(c, j, c.ny, r, p, tm);
        ifftn_op(c, p, hf, false, ratio, tm);
        fftn_half(c, r, c.ldr, r, p, tm);
        // the last update into the ratio buffer (free during the last inverse), or in place into J when it is
        // transposed into the result
        const Update up = {k == n_iter && !native_io ? r : j, j};
        ifftn_op(c, p, hb, false, up, tm);
    }
    if (native_io) {
        tm.start();
        transpose_planes(j, r, c.nz, c.nx, c.ny, sizeof(float), c.nt);
        tm.stop(3);
    }
    return rb.release();
}

float* matlab_rlstep(RfftPlan& c, const float* y_in, const float* im_in, const cfloat* h, bool native_io,
                     Times& times) {
    Timer tm(times);
    Owned<float> rb, yc, ic;
    float* r = ratio_buffer(c, rb);
    tm.start();
    const float* y = matlab_input(c, y_in, native_io, yc);
    const float* im = matlab_input(c, im_in, native_io, ic);
    Owned<cfloat> pb(static_cast<size_t>(c.nP));
    cfloat* p = pb.p;
    tm.stop(3);
    const Ratio<float> ratio = {r, c.ldr, im, 0.0};
    fftn_half(c, y, c.ny, r, p, tm);
    ifftn_op(c, p, h, false, ratio, tm);
    fftn_half(c, r, c.ldr, r, p, tm);
    // the update into the ratio buffer, or in place into the copy of Y that is transposed into the result
    const Update up = {native_io ? yc.p : r, y};
    ifftn_op(c, p, h, true, up, tm);
    if (native_io) {
        tm.start();
        transpose_planes(yc.p, r, c.nz, c.nx, c.ny, sizeof(float), c.nt);
        tm.stop(3);
    }
    return rb.release();
}

void matlab_fft1(RfftPlan& c, const float* x, cfloat* y) {
    fwd_y(c, x, c.ny, reinterpret_cast<float*>(y));
}

template float* matlab_omw<uint16_t>(RfftPlan&, const uint16_t*, double, const cfloat*, const cfloat*, int,
                                     const float*, bool, Times&);
template float* matlab_omw<float>(RfftPlan&, const float*, double, const cfloat*, const cfloat*, int,
                                  const float*, bool, Times&);

}  // namespace pk5d::decon
