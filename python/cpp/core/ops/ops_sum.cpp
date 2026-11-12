// sum(x(:), 'double') in MATLAB's order at one thread, at any thread count (K08; the order of
// decon_sum_double and decon_lucy_accel_mex of the MATLAB code).
//
// ops_sum_double(x, num_threads): x 1-D contiguous, summed in its order (MATLAB's linear order
//   of a MATLAB-layout array, or a vector).
// ops_sum_double_yx(x, num_threads): x (planes, ny, nx) C-contiguous, the (z, y, x) volume of
//   the port; summed in MATLAB's linear order t = y + ny * (x + nx * plane) (y fastest): tiles
//   of columns transposed into a scratch per thread (yx_block_sums), strided reads below
//   89000 elements. It costs a transposing read: 1.6 to 1.8 times the contiguous sum at 1
//   and 4 threads (K08's report); a kernel that holds MATLAB's layout sums it contiguously.
//
// The order, by the number n of elements, in MATLAB's linear index t:
//   n <= 43: in index order, ((x1 + x2) + x3) + ...;
//   44 <= n < 89000: four lanes (element t in lane mod(t, 4), each lane in increasing t), added
//     as ((l1 + l2) + l3) + l4;
//   n >= 89000: blocks of 4096 elements (block t / 4096; the last one shorter), each in four
//     lanes as above, the block sums added in index order. Each block is summed by one thread,
//     so the result doesn't depend on the number of threads (MATLAB splits the blocks among its
//     threads, which changes its last bits).
// Accumulators start at -0 (a sum of -0s is -0); an empty sum is 0. Logical, single, double
// and 8 to 32-bit integer data (converted to double exactly).
//
// A NaN sum has the payload MATLAB gives it, which follows the operand order of each addition
// (x86 returns the first NaN operand, quieted). Measured in R2024b: the running sums are the
// first operand, except in the lanes' tails (the last mod(m, 4) elements of the four-lane sum,
// and of the last block of single data) and for the block sums of double data, which are the
// first operand of their additions. A sum that comes out NaN is summed again in that order
// (exact_sum), since the fast loops leave the operand order to the compiler.
#include "ops_common.hpp"

#include <memory>
#include <vector>

#if defined(__SSE2__) || defined(_M_X64)
#include <emmintrin.h>
#define PK5D_OPS_SSE2 1
#endif

namespace pk5d::ops {
namespace {

constexpr idx_t kSequentialUpTo = 43;
constexpr idx_t kBlocksFrom = 89000;
constexpr idx_t kBlock = 4096;

using flat3_in = nb::ndarray<nb::ro, nb::ndim<3>, nb::c_contig, nb::device::cpu>;

// MATLAB's linear order of a contiguous vector: element t is p[t]
template <typename T>
struct linear_source {
    const T *p;
    static constexpr bool direct = true;
    void fill(idx_t t0, idx_t m, double *out) const {
        for (idx_t k = 0; k < m; k++) out[k] = static_cast<double>(p[t0 + k]);
    }
};

// MATLAB's linear order of a (planes, ny, nx) C array: t = y + ny * (x + nx * plane)
template <typename T>
struct yx_source {
    const T *p;
    idx_t ny, nx;
    static constexpr bool direct = false;
    void fill(idx_t t0, idx_t m, double *out) const {
        const idx_t plane = ny * nx;
        idx_t pl = t0 / plane, r = t0 % plane;
        idx_t col = r / ny, y = r % ny;
        idx_t k = 0;
        while (k < m) {
            const idx_t run = ny - y < m - k ? ny - y : m - k;
            const T *c = p + pl * plane + y * nx + col;
            for (idx_t j = 0; j < run; j++) out[k + j] = static_cast<double>(c[j * nx]);
            k += run;
            y = 0;
            if (++col == nx) {
                col = 0;
                ++pl;
            }
        }
    }
    void block_sums(idx_t n, idx_t b0, idx_t b1, double *bs) const;
};

// adds m elements whose first one has linear index pos to the lanes a (element t to lane
// mod(t, 4), each lane in increasing t); lanes 1-2 and 3-4 as vector pairs where possible
template <typename T>
inline void lanes_add(double *a, const T *p, idx_t m, idx_t pos) {
    idx_t j = 0;
    for (; j < m && ((pos + j) & 3) != 0; j++) {
        const idx_t l = (pos + j) & 3;
        a[l] = a[l] + static_cast<double>(p[j]);
    }
#ifdef PK5D_OPS_SSE2
    if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
        __m128d v01 = _mm_loadu_pd(a), v23 = _mm_loadu_pd(a + 2);
        for (; j + 4 <= m; j += 4) {
            if constexpr (std::is_same_v<T, float>) {
                const __m128 f = _mm_loadu_ps(p + j);
                v01 = _mm_add_pd(v01, _mm_cvtps_pd(f));
                v23 = _mm_add_pd(v23, _mm_cvtps_pd(_mm_movehl_ps(f, f)));
            } else {
                v01 = _mm_add_pd(v01, _mm_loadu_pd(p + j));
                v23 = _mm_add_pd(v23, _mm_loadu_pd(p + j + 2));
            }
        }
        _mm_storeu_pd(a, v01);
        _mm_storeu_pd(a + 2, v23);
    }
#endif
    double a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
    for (; j + 4 <= m; j += 4) {
        a0 = a0 + static_cast<double>(p[j]);
        a1 = a1 + static_cast<double>(p[j + 1]);
        a2 = a2 + static_cast<double>(p[j + 2]);
        a3 = a3 + static_cast<double>(p[j + 3]);
    }
    if (j < m) a0 = a0 + static_cast<double>(p[j]);
    if (j + 1 < m) a1 = a1 + static_cast<double>(p[j + 1]);
    if (j + 2 < m) a2 = a2 + static_cast<double>(p[j + 2]);
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    a[3] = a3;
}

inline double add_lanes(const double *a) { return ((a[0] + a[1]) + a[2]) + a[3]; }

// the four lanes of m elements: ((l1 + l2) + l3) + l4
template <typename T>
inline double lanes_sum(const T *p, idx_t m) {
    double a[4] = {-0.0, -0.0, -0.0, -0.0};
    lanes_add(a, p, m, 0);
    return add_lanes(a);
}

// the block sums of a run of elements in MATLAB's order, fed in pieces
struct block_stream {
    double *bs;  // the block sums, from block b
    idx_t b;
    idx_t pos = 0;
    double a[4] = {-0.0, -0.0, -0.0, -0.0};

    template <typename T>
    void feed(const T *p, idx_t m) {
        while (m > 0) {
            const idx_t k = m < kBlock - pos ? m : kBlock - pos;
            lanes_add(a, p, k, pos);
            pos += k;
            p += k;
            m -= k;
            if (pos == kBlock) flush();
        }
    }
    void flush() {
        bs[b++] = add_lanes(a);
        a[0] = a[1] = a[2] = a[3] = -0.0;
        pos = 0;
    }
};

// A tile of a (planes, ny, nx) C array: w columns of one plane transposed into a scratch
// (column c of the tile at s + c * ld; ld padded so that columns don't share cache sets), and
// the run of MATLAB's linear order it holds: t in [tile + c * ny, tile + (c + 1) * ny) for column
// c, of which [t_cur, t_stop) is still to be fed to the lanes.
template <typename T>
struct yx_tile {
    T *s = nullptr;
    idx_t ld = 0, ny = 0, tile = 0, t_cur = 0, t_stop = 0;

    idx_t left() const { return t_stop - t_cur; }

    // feeds up to budget elements, in MATLAB's order
    void feed(block_stream &st, idx_t budget) {
        while (budget > 0 && t_cur < t_stop) {
            const idx_t c = (t_cur - tile) / ny;
            const idx_t col_end = tile + (c + 1) * ny < t_stop ? tile + (c + 1) * ny : t_stop;
            const idx_t take = col_end - t_cur < budget ? col_end - t_cur : budget;
            st.feed(s + c * ld + (t_cur - tile - c * ny), take);
            t_cur += take;
            budget -= take;
        }
    }
};

// rows [y0, y0 + yb) of wc columns from src (row stride nx) into the scratch columns dst (stride
// ld); 4 x 4 transposes in registers for single (a copy: the values don't change)
template <typename T>
inline void transpose_rows(const T *src, idx_t nx, idx_t yb, idx_t wc, T *dst, idx_t ld) {
    idx_t c = 0;
#ifdef PK5D_OPS_SSE2
    if constexpr (std::is_same_v<T, float>) {
        if (yb == 4) {
            for (; c + 4 <= wc; c += 4) {
                __m128 r0 = _mm_loadu_ps(src + c), r1 = _mm_loadu_ps(src + nx + c);
                __m128 r2 = _mm_loadu_ps(src + 2 * nx + c), r3 = _mm_loadu_ps(src + 3 * nx + c);
                _MM_TRANSPOSE4_PS(r0, r1, r2, r3);
                _mm_storeu_ps(dst + c * ld, r0);
                _mm_storeu_ps(dst + (c + 1) * ld, r1);
                _mm_storeu_ps(dst + (c + 2) * ld, r2);
                _mm_storeu_ps(dst + (c + 3) * ld, r3);
            }
        }
    }
#endif
    for (; c < wc; c++)
        for (idx_t k = 0; k < yb; k++) dst[c * ld + k] = src[k * nx + c];
}

// the block sums of blocks [b0, b1) of a (planes, ny, nx) C array in MATLAB's order. Tiles of
// w columns (about 512 KB) are transposed four rows at a time, with the rows ahead prefetched,
// and the previous tile is fed to the lanes in between, so that the reads of one and the
// additions of the other overlap; the tiles are fed in order, each in MATLAB's order.
template <typename T>
void yx_block_sums(const T *p, idx_t ny, idx_t nx, idx_t n, idx_t b0, idx_t b1, double *bs) {
    constexpr idx_t kRows = 4;
    constexpr idx_t kAhead = 8;  // rows prefetched ahead, into L2
    const idx_t plane = ny * nx;
    const idx_t ld = ny + static_cast<idx_t>(64 / sizeof(T));
    idx_t w = static_cast<idx_t>(524288 / sizeof(T)) / ld;
    w = w < 8 ? 8 : w;
    w = w > nx ? nx : w;
    std::unique_ptr<T[]> buf(new T[static_cast<size_t>(2 * w * ld)]);
    block_stream st{bs, b0};
    yx_tile<T> prev;
    const idx_t t_end = b1 * kBlock < n ? b1 * kBlock : n;
    idx_t t = b0 * kBlock;
    for (idx_t k = 0; t < t_end; k++) {
        const idx_t pl = t / plane;
        const idx_t c0 = (t - pl * plane) / ny;
        const idx_t plane_end = (pl + 1) * plane < t_end ? (pl + 1) * plane : t_end;
        const idx_t c_last = (plane_end - 1 - pl * plane) / ny;
        const idx_t c1 = c_last + 1 < c0 + w ? c_last + 1 : c0 + w;
        const idx_t wc = c1 - c0;
        const T *base = p + pl * plane + c0;
        T *s = buf.get() + (k & 1) * w * ld;
        const idx_t groups = (ny + kRows - 1) / kRows;
        const idx_t share = (prev.left() + groups - 1) / groups;
        for (idx_t y0 = 0; y0 < ny; y0 += kRows) {
            const idx_t yb = ny - y0 < kRows ? ny - y0 : kRows;
#ifdef PK5D_OPS_SSE2
            if (y0 + kAhead < ny) {
                const idx_t ya = y0 + kAhead + kRows <= ny ? kRows : ny - y0 - kAhead;
                for (idx_t r = 0; r < ya; r++)
                    for (idx_t c = 0; c < wc; c += static_cast<idx_t>(64 / sizeof(T)))
                        _mm_prefetch(
                            reinterpret_cast<const char *>(base + (y0 + kAhead + r) * nx + c),
                            _MM_HINT_T1);
            }
#endif
            transpose_rows(base + y0 * nx, nx, yb, wc, s + y0, ld);
            prev.feed(st, share);
        }
        prev.feed(st, prev.left());
        const idx_t tile = pl * plane + c0 * ny;
        prev.s = s;
        prev.ld = ld;
        prev.ny = ny;
        prev.tile = tile;
        prev.t_cur = t;
        prev.t_stop = tile + wc * ny < t_end ? tile + wc * ny : t_end;
        t = prev.t_stop;
    }
    prev.feed(st, prev.left());
    if (st.pos > 0) st.flush();
}

template <typename T>
void yx_source<T>::block_sums(idx_t n, idx_t b0, idx_t b1, double *bs) const {
    yx_block_sums(p, ny, nx, n, b0, b1, bs);
}

// a + b with the NaN of a, else of b, quieted (as x86 adds a NaN a and b)
inline double add_first(double a, double b) {
    if (a != a) return a + a;
    if (b != b) return b + b;
    return a + b;
}

// the four lanes of m values with MATLAB's operand order; tail_first: the last mod(m, 4)
// values are the first operand of their additions
double lanes_sum_exact(const double *p, idx_t m, bool tail_first) {
    double a[4] = {-0.0, -0.0, -0.0, -0.0};
    const idx_t main = m - m % 4;
    for (idx_t t = 0; t < main; t++) a[t & 3] = add_first(a[t & 3], p[t]);
    for (idx_t t = main; t < m; t++)
        a[t & 3] = tail_first ? add_first(p[t], a[t & 3]) : add_first(a[t & 3], p[t]);
    return add_first(add_first(add_first(a[0], a[1]), a[2]), a[3]);
}

// below 89000 elements, on the values in MATLAB's order
double small_sum(const double *v, idx_t n) {
    if (n <= kSequentialUpTo) {
        double s = -0.0;
        for (idx_t i = 0; i < n; i++) s = s + v[i];
        return s;
    }
    return lanes_sum(v, n);
}

double small_sum_exact(const double *v, idx_t n) {
    if (n <= kSequentialUpTo) {
        double s = -0.0;
        for (idx_t i = 0; i < n; i++) s = add_first(s, v[i]);
        return s;
    }
    return lanes_sum_exact(v, n, true);
}

// the blocks again on one thread with MATLAB's operand order (the payload of a NaN result)
template <typename Src>
double blocks_exact(const Src &src, idx_t n, bool single) {
    const idx_t nb = (n + kBlock - 1) / kBlock;
    std::vector<double> buf(static_cast<size_t>(kBlock));
    double s = -0.0;
    for (idx_t b = 0; b < nb; b++) {
        const idx_t i0 = b * kBlock;
        const idx_t m = b < nb - 1 ? kBlock : n - i0;
        src.fill(i0, m, buf.data());
        const double bs = lanes_sum_exact(buf.data(), m, single);
        s = single ? add_first(s, bs) : add_first(bs, s);
    }
    return s;
}

template <typename T, typename Src>
double fixed_sum(const Src &src, idx_t n, int num_threads) {
    constexpr bool single = std::is_same_v<T, float>;
    if (n == 0) return 0.0;
    if (n < kBlocksFrom) {
        std::vector<double> v(static_cast<size_t>(n));
        src.fill(0, n, v.data());
        const double s = small_sum(v.data(), n);
        if (is_float_v<T> && s != s) return small_sum_exact(v.data(), n);
        return s;
    }
    const idx_t nb = (n + kBlock - 1) / kBlock;
    std::vector<double> bs(static_cast<size_t>(nb));
    const int nt = team(n, num_threads);
    if constexpr (Src::direct) {
#pragma omp parallel for num_threads(nt) schedule(static)
        for (idx_t b = 0; b < nb; b++) {
            const idx_t i0 = b * kBlock;
            const idx_t m = b < nb - 1 ? kBlock : n - i0;
            bs[static_cast<size_t>(b)] = lanes_sum(src.p + i0, m);
        }
    } else {
        // each thread a range of whole blocks
#pragma omp parallel num_threads(nt)
        {
#ifdef _OPENMP
            const idx_t k = omp_get_thread_num(), nk = omp_get_num_threads();
#else
            const idx_t k = 0, nk = 1;
#endif
            const idx_t b0 = nb * k / nk, b1 = nb * (k + 1) / nk;
            if (b1 > b0) src.block_sums(n, b0, b1, bs.data());
        }
    }
    double s = -0.0;
    for (idx_t b = 0; b < nb; b++) s = s + bs[static_cast<size_t>(b)];
    if (is_float_v<T> && s != s) return blocks_exact(src, n, single);
    return s;
}

template <typename Fn>
void with_sum_class(nb::dlpack::dtype dt, Fn &&fn) {
    with_class(class_of(dt), "ops_sum_double input", [&](auto ts) {
        using T = typename decltype(ts)::type;
        if constexpr (std::is_same_v<T, int64_t> || std::is_same_v<T, uint64_t>) {
            bad_class("ops_sum_double input (64-bit integers)");
        } else {
            fn(ts);
        }
    });
}

double sum_double(const flat_in &x, int num_threads) {
    check_threads(num_threads);
    const idx_t n = static_cast<idx_t>(x.shape(0));
    const void *px = x.data();
    double s = 0.0;
    nb::gil_scoped_release release;
    with_sum_class(x.dtype(), [&](auto ts) {
        using T = typename decltype(ts)::type;
        s = fixed_sum<T>(linear_source<T>{static_cast<const T *>(px)}, n, num_threads);
    });
    return s;
}

double sum_double_yx(const flat3_in &x, int num_threads) {
    check_threads(num_threads);
    const idx_t np_ = static_cast<idx_t>(x.shape(0));
    const idx_t ny = static_cast<idx_t>(x.shape(1));
    const idx_t nx = static_cast<idx_t>(x.shape(2));
    const idx_t n = np_ * ny * nx;
    const void *px = x.data();
    double s = 0.0;
    nb::gil_scoped_release release;
    with_sum_class(x.dtype(), [&](auto ts) {
        using T = typename decltype(ts)::type;
        const T *p = static_cast<const T *>(px);
        if (ny == 1 || nx == 1)
            s = fixed_sum<T>(linear_source<T>{p}, n, num_threads);
        else
            s = fixed_sum<T>(yx_source<T>{p, ny, nx}, n, num_threads);
    });
    return s;
}

}  // namespace
}  // namespace pk5d::ops

PK5D_BINDINGS(bind_ops_sum) {
    namespace nb = nanobind;
    m.def("ops_sum_double", &pk5d::ops::sum_double, nb::arg("x").noconvert(),
          nb::arg("num_threads"),
          "sum(x(:), 'double') of a vector in its order, MATLAB's at one thread, independent of "
          "num_threads.");
    m.def("ops_sum_double_yx", &pk5d::ops::sum_double_yx, nb::arg("x").noconvert(),
          nb::arg("num_threads"),
          "sum(v(:), 'double') of a (planes, ny, nx) C array in MATLAB's linear order (y "
          "fastest), MATLAB's order at one thread, independent of num_threads.");
}
