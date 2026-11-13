// The engine of MATLAB's vector rule (K08): lanes, 4096-blocks, sources in MATLAB's linear
// order, additions in single or double of the values or of abs(x - center).^2. ops_sum.cpp
// documents the order and binds it.
#pragma once

#include "ops_common.hpp"

#include <atomic>
#include <cmath>
#include <memory>
#include <vector>

#if defined(__SSE2__) || defined(_M_X64)
#include <emmintrin.h>
#define PK5D_OPS_SSE2 1
#endif

namespace pk5d::ops {

constexpr idx_t kSequentialUpTo = 43;
constexpr idx_t kBlocksFrom = 89000;
constexpr idx_t kBlock = 4096;

// the value added for an element: x converted to the class of the sum
template <typename T, typename Acc>
struct op_value {
    using src_t = T;
    using acc_t = Acc;
    static constexpr bool sq = false;
    Acc operator()(T v) const { return static_cast<Acc>(v); }
};

// abs(x - q).^2 in the class of x (var's sum; q is the mean in that class)
template <typename T>
struct op_sqdev {
    using src_t = T;
    using acc_t = T;
    static constexpr bool sq = true;
    T q = T(0);
    T operator()(T v) const {
        const T d = std::fabs(v - q);
        return d * d;
    }
};

// adds m elements whose first one has linear index pos to the lanes a (element t to lane
// mod(t, 4), each lane in increasing t); vector registers where the classes allow, with the same
// operations per element
template <typename Op>
inline void lanes_add(typename Op::acc_t *a, const typename Op::src_t *p, idx_t m, idx_t pos,
                      const Op &op) {
    using T = typename Op::src_t;
    using Acc = typename Op::acc_t;
    idx_t j = 0;
    for (; j < m && ((pos + j) & 3) != 0; j++) {
        const idx_t l = (pos + j) & 3;
        a[l] = a[l] + op(p[j]);
    }
#ifdef PK5D_OPS_SSE2
    if constexpr (std::is_same_v<Acc, double> && (std::is_same_v<T, float> ||
                                                 (std::is_same_v<T, double> && !Op::sq))) {
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
    } else if constexpr (std::is_same_v<Acc, double> && std::is_same_v<T, double> && Op::sq) {
        const __m128d q = _mm_set1_pd(op.q), sign = _mm_set1_pd(-0.0);
        __m128d v01 = _mm_loadu_pd(a), v23 = _mm_loadu_pd(a + 2);
        for (; j + 4 <= m; j += 4) {
            __m128d d0 = _mm_andnot_pd(sign, _mm_sub_pd(_mm_loadu_pd(p + j), q));
            __m128d d1 = _mm_andnot_pd(sign, _mm_sub_pd(_mm_loadu_pd(p + j + 2), q));
            v01 = _mm_add_pd(v01, _mm_mul_pd(d0, d0));
            v23 = _mm_add_pd(v23, _mm_mul_pd(d1, d1));
        }
        _mm_storeu_pd(a, v01);
        _mm_storeu_pd(a + 2, v23);
    } else if constexpr (std::is_same_v<Acc, float> && std::is_same_v<T, float>) {
        __m128 v = _mm_loadu_ps(a);
        if constexpr (Op::sq) {
            const __m128 q = _mm_set1_ps(op.q), sign = _mm_set1_ps(-0.0f);
            for (; j + 4 <= m; j += 4) {
                const __m128 d = _mm_andnot_ps(sign, _mm_sub_ps(_mm_loadu_ps(p + j), q));
                v = _mm_add_ps(v, _mm_mul_ps(d, d));
            }
        } else {
            for (; j + 4 <= m; j += 4) v = _mm_add_ps(v, _mm_loadu_ps(p + j));
        }
        _mm_storeu_ps(a, v);
    }
#endif
    Acc a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
    for (; j + 4 <= m; j += 4) {
        a0 = a0 + op(p[j]);
        a1 = a1 + op(p[j + 1]);
        a2 = a2 + op(p[j + 2]);
        a3 = a3 + op(p[j + 3]);
    }
    if (j < m) a0 = a0 + op(p[j]);
    if (j + 1 < m) a1 = a1 + op(p[j + 1]);
    if (j + 2 < m) a2 = a2 + op(p[j + 2]);
    a[0] = a0;
    a[1] = a1;
    a[2] = a2;
    a[3] = a3;
}

template <typename Acc>
inline Acc add_lanes(const Acc *a) {
    return ((a[0] + a[1]) + a[2]) + a[3];
}

// the four lanes of m elements: ((l1 + l2) + l3) + l4
template <typename Op>
inline typename Op::acc_t lanes_sum(const typename Op::src_t *p, idx_t m, const Op &op) {
    using Acc = typename Op::acc_t;
    Acc a[4] = {Acc(-0.0), Acc(-0.0), Acc(-0.0), Acc(-0.0)};
    lanes_add(a, p, m, 0, op);
    return add_lanes(a);
}

// below 89000 elements, on the values (already in the class of the sum) in MATLAB's order
template <typename Acc>
inline Acc small_sum(const Acc *v, idx_t n) {
    if (n <= kSequentialUpTo) {
        Acc s = Acc(-0.0);
        for (idx_t i = 0; i < n; i++) s = s + v[i];
        return s;
    }
    return lanes_sum(v, n, op_value<Acc, Acc>{});
}

// MATLAB's linear order of a contiguous vector: element t is p[t]
template <typename T>
struct linear_source {
    const T *p;
    static constexpr bool direct = true;
    template <typename Op>
    void fill(idx_t t0, idx_t m, typename Op::acc_t *out, const Op &op) const {
        for (idx_t k = 0; k < m; k++) out[k] = op(p[t0 + k]);
    }
};

// the block sums of a run of elements in MATLAB's order, fed in pieces
template <typename Op>
struct block_stream {
    using Acc = typename Op::acc_t;
    const Op &op;
    Acc *bs;  // the block sums, from block b
    idx_t b;
    idx_t pos = 0;
    Acc a[4] = {Acc(-0.0), Acc(-0.0), Acc(-0.0), Acc(-0.0)};

    void feed(const typename Op::src_t *p, idx_t m) {
        while (m > 0) {
            const idx_t k = m < kBlock - pos ? m : kBlock - pos;
            lanes_add(a, p, k, pos, op);
            pos += k;
            p += k;
            m -= k;
            if (pos == kBlock) flush();
        }
    }
    void flush() {
        bs[b++] = add_lanes(a);
        a[0] = a[1] = a[2] = a[3] = Acc(-0.0);
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
    template <typename Stream>
    void feed(Stream &st, idx_t budget) {
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
template <typename Op>
void yx_block_sums(const typename Op::src_t *p, idx_t ny, idx_t nx, idx_t n, idx_t b0, idx_t b1,
                   typename Op::acc_t *bs, const Op &op) {
    using T = typename Op::src_t;
    constexpr idx_t kRows = 4;
    constexpr idx_t kAhead = 8;  // rows prefetched ahead, into L2
    const idx_t plane = ny * nx;
    const idx_t ld = ny + static_cast<idx_t>(64 / sizeof(T));
    idx_t w = static_cast<idx_t>(524288 / sizeof(T)) / ld;
    w = w < 8 ? 8 : w;
    w = w > nx ? nx : w;
    std::unique_ptr<T[]> buf(new T[static_cast<size_t>(2 * w * ld)]);
    block_stream<Op> st{op, bs, b0};
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

// MATLAB's linear order of a (planes, ny, nx) C array: t = y + ny * (x + nx * plane)
template <typename T>
struct yx_source {
    const T *p;
    idx_t ny, nx;
    static constexpr bool direct = false;
    template <typename Op>
    void fill(idx_t t0, idx_t m, typename Op::acc_t *out, const Op &op) const {
        const idx_t plane = ny * nx;
        idx_t pl = t0 / plane, r = t0 % plane;
        idx_t col = r / ny, y = r % ny;
        idx_t k = 0;
        while (k < m) {
            const idx_t run = ny - y < m - k ? ny - y : m - k;
            const T *c = p + pl * plane + y * nx + col;
            for (idx_t j = 0; j < run; j++) out[k + j] = op(c[j * nx]);
            k += run;
            y = 0;
            if (++col == nx) {
                col = 0;
                ++pl;
            }
        }
    }
    template <typename Op>
    void block_sums(idx_t n, idx_t b0, idx_t b1, typename Op::acc_t *bs, const Op &op) const {
        yx_block_sums(p, ny, nx, n, b0, b1, bs, op);
    }
};

// the blocks of 4096 (n >= 89000): the block sums in threads, added in order afterwards. Which
// thread sums a block doesn't change the result: contiguous chunks of blocks are handed out as
// threads ask for them (threads that start late or run slower do less); the transposing source
// takes one range of whole blocks per thread
template <typename Src, typename Op>
typename Op::acc_t fixed_sum(const Src &src, idx_t n, int num_threads, bool bind, const Op &op) {
    using Acc = typename Op::acc_t;
    if (n == 0) return Acc(0);
    const idx_t nb = (n + kBlock - 1) / kBlock;
    std::unique_ptr<Acc[]> bs(new Acc[static_cast<size_t>(nb)]);  // every one is written
    const int nt = team(n, num_threads);
    std::atomic<idx_t> next{0};
    constexpr idx_t kChunk = 16;  // blocks (256 KB of single) per request (16 to 64 measured best)
    parallel(nt, bind, [&](int k, int nk) {
        if constexpr (Src::direct) {
            // one thread alone has too few reads in flight to fill the memory bandwidth:
            // it prefetches the block after next into L2 (about 5-8% at one thread; with
            // more threads the bandwidth is the limit and the prefetches cost, measured)
            const bool ahead = nk == 1;
            constexpr idx_t kAheadBlocks = 2;
            for (;;) {
                const idx_t b0 = next.fetch_add(kChunk, std::memory_order_relaxed);
                if (b0 >= nb) break;
                const idx_t b1 = b0 + kChunk < nb ? b0 + kChunk : nb;
                for (idx_t b = b0; b < b1; b++) {
                    const idx_t i0 = b * kBlock;
                    const idx_t m = b < nb - 1 ? kBlock : n - i0;
#ifdef PK5D_OPS_SSE2
                    if (ahead && b + kAheadBlocks < nb) {
                        const char *q =
                            reinterpret_cast<const char *>(src.p + (b + kAheadBlocks) * kBlock);
                        const idx_t bytes = kBlock * static_cast<idx_t>(sizeof(*src.p));
                        for (idx_t o = 0; o < bytes; o += 64) _mm_prefetch(q + o, _MM_HINT_T1);
                    }
#endif
                    bs[static_cast<size_t>(b)] = lanes_sum(src.p + i0, m, op);
                }
            }
        } else {
            idx_t b0, b1;
            share(nb, k, nk, b0, b1);
            if (b1 > b0) src.block_sums(n, b0, b1, bs.get(), op);
        }
    });
    Acc s = Acc(-0.0);
    for (idx_t b = 0; b < nb; b++) s = s + bs[static_cast<size_t>(b)];
    return s;
}

// the whole vector rule, any n (no payload rule for NaN results)
template <typename Src, typename Op>
typename Op::acc_t ordered_sum(const Src &src, idx_t n, int num_threads, bool bind, const Op &op) {
    using Acc = typename Op::acc_t;
    if (n == 0) return Acc(0);
    if (n < kBlocksFrom) {
        std::vector<Acc> v(static_cast<size_t>(n));
        src.fill(0, n, v.data(), op);
        return small_sum(v.data(), n);
    }
    return fixed_sum(src, n, num_threads, bind, op);
}

}  // namespace pk5d::ops
