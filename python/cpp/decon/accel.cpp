// The acceleration step of decon_lucy_function's Richardson-Lucy iterations (decon_lucy_accel_mex) and its
// sums in double, in an order that doesn't depend on the number of threads: MATLAB's at one thread for
// arrays of 89000 or more elements. Blocks of 4096 elements in MATLAB's linear order (the last one shorter),
// each summed in four lanes (element t of a block in lane mod(t, 4), in index order), the lanes added as
// ((l1 + l2) + l3) + l4 and the block sums in index order. Each block is summed by one thread. The
// accumulators start at -0, so a sum of -0s is -0 as in MATLAB; an empty sum is 0.
//
// The arrays are flat in MATLAB's linear order (MATLAB's layout, a vector), or the port's (planes, ny, nx)
// C arrays, read in MATLAB's linear order t = y + ny * (x + nx * plane): a block is a run of columns (y
// lines, strided by nx); consecutive columns share their cache lines, so a column's reads hit the lines
// the previous one loaded.
//
// The element-wise arithmetic is MATLAB's: differences and products in the class of the data (single
// products converted to double for the sums), a double scalar times a single array computed in double and
// rounded to single, max(v, 0) is 0 for NaN and keeps -0 (MATLAB returns the first argument on a tie), and
// min(NaN, 1) is 1. A sum that comes out NaN is summed again with the operands of each addition in the
// mex's order (the running sum first; x86 returns the first NaN operand, quieted), since the fast loops
// leave the operand order to the compiler.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

#include "decon.hpp"

namespace pk5d::decon {
namespace {

constexpr idx_t kBlock = 4096;

inline idx_t n_blocks(idx_t n) { return (n + kBlock - 1) / kBlock; }

// a + b with the NaN of a, else of b, quieted (as x86 adds a NaN a and b)
inline double add_first(double a, double b) {
    if (a != a) return a + a;
    if (b != b) return b + b;
    return a + b;
}

template <bool Exact>
inline double add(double a, double b) {
    if constexpr (Exact) {
        return add_first(a, b);
    } else {
        return a + b;
    }
}

template <bool Exact>
inline double add_lanes(const double* a) {
    return add<Exact>(add<Exact>(add<Exact>(a[0], a[1]), a[2]), a[3]);
}

// the block sums in index order
template <bool Exact>
double add_blocks(const std::vector<double>& b) {
    if (b.empty()) return 0.0;
    double s = -0.0;
    for (double v : b) s = add<Exact>(s, v);
    return s;
}

// the lambda step at element i: d = J_2 - Y, the products d .* J_4 and J_4 .* J_4 (in T) as doubles, and
// (W) J_4 = d written
template <typename T, bool W = true>
struct LambdaElem {
    static constexpr int K = 2;
    const T* j2;
    const T* y;
    const T* j4;
    T* out;
    inline void at(idx_t i, double* v) const {
        const T d = j2[i] - y[i];
        const T w = j4[i];
        const T p1 = d * w;
        const T p2 = w * w;
        v[0] = static_cast<double>(p1);
        v[1] = static_cast<double>(p2);
        if constexpr (W) out[i] = d;
    }
};

template <typename T>
struct SumElem {
    static constexpr int K = 1;
    const T* x;
    inline void at(idx_t i, double* v) const { v[0] = static_cast<double>(x[i]); }
};

// elements o, o + s, ..., o + (m - 1) s, the first in lane l
template <typename E, bool Exact>
inline void run(const E& e, double (*a)[4], idx_t o, idx_t s, idx_t m, int l) {
    constexpr int K = E::K;
    double v[K];
    idx_t k = 0;
    if (s == 1 && l == 0) {
        // the lanes in turn, four elements at a time
        for (; k + 4 <= m; k += 4) {
            for (int q = 0; q < 4; q++) {
                e.at(o + k + q, v);
                for (int c = 0; c < K; c++) a[c][q] = add<Exact>(a[c][q], v[c]);
            }
        }
    }
    for (; k < m; k++) {
        e.at(o + k * s, v);
        for (int c = 0; c < K; c++) a[c][l] = add<Exact>(a[c][l], v[c]);
        l = (l + 1) & 3;
    }
}

// the sums of the blocks [b0, b1) of every stream, flat arrays
template <typename E, bool Exact>
void blocks(const E& e, idx_t n, idx_t b0, idx_t b1, double* const* bs) {
    constexpr int K = E::K;
    for (idx_t b = b0; b < b1; b++) {
        const idx_t t0 = b * kBlock, t1 = std::min(t0 + kBlock, n);
        double a[K][4];
        for (int c = 0; c < K; c++) a[c][0] = a[c][1] = a[c][2] = a[c][3] = -0.0;
        run<E, Exact>(e, a, t0, 1, t1 - t0, 0);
        for (int c = 0; c < K; c++) bs[c][b] = add_lanes<Exact>(a[c]);
    }
}

// columns c0 .. c1 - 1 of a (ny, nx) plane (row stride nx) into dst (column c at (c - c0) * ny), each
// element f(y, x); blocked, so that reads and writes stay in cache lines
template <typename T, typename F>
void columns_of(idx_t ny, idx_t nx, idx_t c0, idx_t c1, T* dst, const F& f) {
    constexpr idx_t kT = 16;
    for (idx_t y0 = 0; y0 < ny; y0 += kT) {
        const idx_t y1 = std::min(y0 + kT, ny);
        for (idx_t x0 = c0; x0 < c1; x0 += kT) {
            const idx_t x1 = std::min(x0 + kT, c1);
            for (idx_t y = y0; y < y1; y++) {
                for (idx_t x = x0; x < x1; x++) dst[(x - c0) * ny + y] = f(y, x);
            }
        }
    }
}

// the block sums of the elements [t_lo, t_hi) (whole blocks) of a (planes, ny, nx) array in MATLAB's
// order: per plane, the columns the range covers are transposed into a scratch per thread (in chunks of
// columns), whose elements are then contiguous in MATLAB's order. K streams: lambda (d = J_2 - Y and J_4,
// their products) or the values of x.
// ceil(a / b) for b > 0
inline idx_t ceil_div(idx_t a, idx_t b) { return a >= 0 ? (a + b - 1) / b : -((-a) / b); }

// d = J_2 - Y of columns c0 .. c1 - 1 of a plane into the scratch (as columns_of), and into out (row by row,
// contiguous) where the element's index in the plane's MATLAB order is in [lo, hi)
template <typename T>
void differences_of(const T* a0, const T* a1, T* out, idx_t ny, idx_t nx, idx_t c0, idx_t c1, idx_t lo, idx_t hi,
                    T* dst) {
    constexpr idx_t kT = 16;
    for (idx_t y0 = 0; y0 < ny; y0 += kT) {
        const idx_t y1 = std::min(y0 + kT, ny);
        for (idx_t x0 = c0; x0 < c1; x0 += kT) {
            const idx_t x1 = std::min(x0 + kT, c1);
            for (idx_t y = y0; y < y1; y++) {
                const T* r0 = a0 + y * nx;
                const T* r1 = a1 + y * nx;
                for (idx_t x = x0; x < x1; x++) dst[(x - c0) * ny + y] = r0[x] - r1[x];
                if (out != nullptr) {
                    const idx_t xa = std::max(x0, ceil_div(lo - y, ny)), xb = std::min(x1, ceil_div(hi - y, ny));
                    T* ro = out + y * nx;
                    for (idx_t x = xa; x < xb; x++) ro[x] = r0[x] - r1[x];
                }
            }
        }
    }
}

// The lanes of m scratch elements from k0, the first in lane l; K streams (lambda: p0 = d, p1 = J_4)
template <typename T, int K, bool Exact>
inline void feed(double (*acc)[4], const T* p0, const T* p1, idx_t m, int l) {
    auto one = [&](idx_t k, int q) {
        if constexpr (K == 2) {
            const T d = p0[k], w = p1[k];
            const T v1 = d * w;
            const T v2 = w * w;
            acc[0][q] = add<Exact>(acc[0][q], static_cast<double>(v1));
            acc[1][q] = add<Exact>(acc[1][q], static_cast<double>(v2));
        } else {
            acc[0][q] = add<Exact>(acc[0][q], static_cast<double>(p0[k]));
        }
    };
    idx_t k = 0;
    for (; k < m && l != 0; k++, l = (l + 1) & 3) one(k, l);
    if constexpr (!Exact) {
        // four at a time in local lanes
        double a[K][4];
        for (int c = 0; c < K; c++)
            for (int q = 0; q < 4; q++) a[c][q] = acc[c][q];
        for (; k + 4 <= m; k += 4) {
            for (int q = 0; q < 4; q++) {
                if constexpr (K == 2) {
                    const T d = p0[k + q], w = p1[k + q];
                    const T v1 = d * w;
                    const T v2 = w * w;
                    a[0][q] = a[0][q] + static_cast<double>(v1);
                    a[1][q] = a[1][q] + static_cast<double>(v2);
                } else {
                    a[0][q] = a[0][q] + static_cast<double>(p0[k + q]);
                }
            }
        }
        for (int c = 0; c < K; c++)
            for (int q = 0; q < 4; q++) acc[c][q] = a[c][q];
    }
    for (; k < m; k++, l = (l + 1) & 3) one(k, l);
}

template <typename T, int K, bool Exact>
void zyx_blocks(const T* a0, const T* a1, const T* a2, T* out, idx_t n, idx_t ny, idx_t nx, idx_t b0, idx_t b1,
                double* const* bs) {
    const idx_t plane = ny * nx;
    const idx_t t_lo = b0 * kBlock, t_hi = std::min(b1 * kBlock, n);
    const idx_t cmax = std::max<idx_t>(1, std::min<idx_t>(nx, (idx_t(1) << 18) / ny));
    std::vector<T> s0(static_cast<size_t>(cmax * ny)), s1(K == 2 ? static_cast<size_t>(cmax * ny) : 0);
    double acc[K][4];
    for (int c = 0; c < K; c++) acc[c][0] = acc[c][1] = acc[c][2] = acc[c][3] = -0.0;
    idx_t b = b0;
    idx_t t = t_lo;
    while (t < t_hi) {
        const idx_t pl = t / plane, r = t - pl * plane;
        const idx_t c0 = r / ny;
        const idx_t c1 = std::min({nx, c0 + cmax, (std::min(t_hi, (pl + 1) * plane) - 1 - pl * plane) / ny + 1});
        const idx_t base = pl * plane;
        // the scratch holds t' in [base + c0 * ny, base + c1 * ny), this thread's part [t, t_end)
        const idx_t s_t0 = base + c0 * ny;
        const idx_t t_end = std::min(t_hi, base + c1 * ny);
        if constexpr (K == 2) {
            // J_4 first: out (J_4 = d, written here for this thread's elements) may be J_4 itself
            columns_of(ny, nx, c0, c1, s1.data(), [&](idx_t y, idx_t x) { return a2[base + y * nx + x]; });
            // whole columns of this thread are written without a check, the partial first and last ones
            // element by element
            differences_of(a0 + base, a1 + base, out != nullptr ? out + base : nullptr, ny, nx, c0, c1, t - base,
                           t_end - base, s0.data());
        } else {
            columns_of(ny, nx, c0, c1, s0.data(), [&](idx_t y, idx_t x) { return a0[base + y * nx + x]; });
        }
        while (t < t_end) {
            const idx_t blk_end = std::min(t_end, (b + 1) * kBlock);
            const idx_t k0 = t - s_t0, m = blk_end - t;
            int l = static_cast<int>(t & 3);
            feed<T, K, Exact>(acc, s0.data() + k0, K == 2 ? s1.data() + k0 : nullptr, m, l);
            t = blk_end;
            if (t == (b + 1) * kBlock || t == n) {
                for (int c = 0; c < K; c++) {
                    bs[c][b] = add_lanes<Exact>(acc[c]);
                    acc[c][0] = acc[c][1] = acc[c][2] = acc[c][3] = -0.0;
                }
                b++;
            }
        }
    }
}

// the block sums of every stream, the blocks split among the threads (contiguous ranges); flat arrays
// through e, (planes, ny, nx) arrays (zyx) through the column scratch
template <typename E, bool Exact, typename T>
void block_sums(const E& e, const T* a0, const T* a1, const T* a2, T* out, idx_t n, bool zyx, idx_t ny, idx_t nx,
                int nt, std::vector<double>* sums) {
    constexpr int K = E::K;
    const idx_t nb = n_blocks(n);
    double* bs[K];
    for (int c = 0; c < K; c++) {
        sums[c].assign(static_cast<size_t>(nb), 0.0);
        bs[c] = sums[c].data();
    }
    const int np = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(nt, nb)));
#pragma omp parallel for num_threads(np) schedule(static, 1)
    for (int q = 0; q < np; q++) {
        const idx_t b0 = nb * q / np, b1 = nb * (q + 1) / np;
        if (b0 >= b1) continue;
        if (zyx) {
            zyx_blocks<T, K, Exact>(a0, a1, a2, out, n, ny, nx, b0, b1, bs);
        } else {
            blocks<E, Exact>(e, n, b0, b1, bs);
        }
    }
}

}  // namespace

double clamped_lambda(const double* s) {
    const double v = s[0] / (s[1] + DBL_EPSILON);
    const double m = std::isnan(v) ? 1.0 : (v < 1.0 ? v : 1.0);
    return m >= 0.0 ? m : 0.0;
}

template <typename T>
void accel_lambda(const T* j2, const T* y, const T* j4, T* out, idx_t n, bool zyx, idx_t ny, idx_t nx, int nt,
                  double* s) {
    // out may be J_4 itself: then the sums run first (a NaN sum is summed again from the same J_4) and
    // J_4 = J_2 - Y is written after them
    const bool inplace = out == j4;
    const LambdaElem<T, true> e{j2, y, j4, out};
    const LambdaElem<T, false> e_sums{j2, y, j4, nullptr};
    std::vector<double> sums[2];
    if (inplace) {
        block_sums<LambdaElem<T, false>, false>(e_sums, j2, y, j4, static_cast<T*>(nullptr), n, zyx, ny, nx, nt,
                                                sums);
    } else {
        block_sums<LambdaElem<T, true>, false>(e, j2, y, j4, out, n, zyx, ny, nx, nt, sums);
    }
    s[0] = add_blocks<false>(sums[0]);
    s[1] = add_blocks<false>(sums[1]);
    if (s[0] != s[0] || s[1] != s[1]) {
        // again with the mex's operand order
        block_sums<LambdaElem<T, false>, true>(e_sums, j2, y, j4, static_cast<T*>(nullptr), n, zyx, ny, nx, nt,
                                               sums);
        s[0] = add_blocks<true>(sums[0]);
        s[1] = add_blocks<true>(sums[1]);
    }
    if (inplace) {
#pragma omp parallel for num_threads(nt) schedule(static)
        for (idx_t i = 0; i < n; i++) out[i] = j2[i] - y[i];
    }
}

template <typename T>
void accel_predict(const T* j2, const T* j3, double lambda, T* y, idx_t n, int nt) {
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        const T d = j2[i] - j3[i];
        const T t = static_cast<T>(lambda * static_cast<double>(d));
        const T u = j2[i] + t;
        // max(u, 0): NaN is 0, -0 is kept
        y[i] = u >= static_cast<T>(0) ? u : static_cast<T>(0);
    }
}

template <typename T>
double accel_sum(const T* x, idx_t n, bool zyx, idx_t ny, idx_t nx, int nt) {
    const SumElem<T> e{x};
    std::vector<double> sums[1];
    block_sums<SumElem<T>, false>(e, x, x, x, static_cast<T*>(nullptr), n, zyx, ny, nx, nt, sums);
    double s = add_blocks<false>(sums[0]);
    if (s != s) {
        block_sums<SumElem<T>, true>(e, x, x, x, static_cast<T*>(nullptr), n, zyx, ny, nx, nt, sums);
        s = add_blocks<true>(sums[0]);
    }
    return s;
}

template void accel_lambda<float>(const float*, const float*, const float*, float*, idx_t, bool, idx_t, idx_t, int,
                                  double*);
template void accel_lambda<double>(const double*, const double*, const double*, double*, idx_t, bool, idx_t, idx_t,
                                   int, double*);
template void accel_predict<float>(const float*, const float*, double, float*, idx_t, int);
template void accel_predict<double>(const double*, const double*, double, double*, idx_t, int);
template double accel_sum<float>(const float*, idx_t, bool, idx_t, idx_t, int);
template double accel_sum<double>(const double*, idx_t, bool, idx_t, idx_t, int);

}  // namespace pk5d::decon
