// The deskew/rotate warp (volume_deskew_rotate_warp_mex) and the fused skewed space
// interpolation and warp (skewed_space_interp_volume_deskew_rotate_warp_mex) on C-ordered
// (z, y, x) volumes, blocked so that the mex's loop runs on y-fastest copies of the input.
//
// Output slice zo reads input columns col, col + 1 (plan.hpp; interpolated grid slices also the
// raw columns up to xa on either side), its line xo grid slices z0, z0 + 1, and a line is the
// rows y. A unit of work is a block of rows (a multiple of 8, the lanes) and a window of grid
// slices z0. It takes the output slices in the order of their columns; their columns are copied
// for the raw slices the window reads, a chunk of 32 columns at a time (whole cache lines of the
// input, the next slice's rows prefetched) through 8 x 8 transposes, into a ring of columns with
// the rows fastest, as single (the next output slice's columns lie a column or two further, so
// each input value is copied about once per block and window). For each slice the unit computes
// the interpolated grid columns col and col + 1 its lines read, then the lines whose z0 is in the
// window (a run of xo, since z0 moves monotonically with xo), each with the mex's sum over its
// rows, into a y-fastest block that goes to the output's rows through 8 x 8 transposes, whole
// cache lines with streaming stores. The scratch of a unit stays in a core's cache (padded lane
// arrays: no power-of-two strides).
//
// Every voxel has the mex's operands, weights and formula (lanes.hpp), chosen by the voxel's
// place in its output line (and, for interpolated values, its row y), never by the blocking:
// the results do not depend on the block and window sizes or the threads. With `fill` the
// voxels the mex leaves alone (slices and lines outside, rows past the input) are written as 0,
// so that the output need not start zeroed (a buffer reused from file to file). Without AVX2, a
// float value (interpolated, or of a line) that comes out NaN is computed again with x86's NaN
// (nans_x86.hpp).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "lanes.hpp"
#include "nans_x86.hpp"
#include "plan.hpp"
#include "skewed_warp.hpp"
#include "transpose.hpp"

namespace pk5d::deskew {

// The sizes only change the speed (measured on a demo tile, 1 to 24 threads):
// rows per block (by default, and at most) and the scratch a window's raw columns may take
constexpr uint64_t BLOCK_ROWS = 32, MAX_BLOCK_ROWS = 64;
constexpr uint64_t RING_BYTES = 320 * 1024;
// lines computed per pass of a unit's slice (their y-fastest block and its rows)
constexpr uint64_t SEG_LINES = 256;
// lanes added to the lane arrays of the scratch
constexpr uint64_t LANE_PAD = 8;
// columns of a chunk of the ring (a 64-byte line of uint16)
constexpr int64_t CHUNK = 32;
// raw slices ahead whose rows the copies prefetch
constexpr int64_t PREFETCH = 1;

inline void prefetch_line(const void* p) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(p, 0, 3);
#else
    (void)p;
#endif
}

// the partial lines at the ends of a run of output bytes, for writing (the rest is streamed)
inline void prefetch_write(const void* p, const uint64_t nbytes) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(p, 1, 3);
    __builtin_prefetch((const char*)p + nbytes - 1, 1, 3);
#else
    (void)p;
    (void)nbytes;
#endif
}

// the raw slices grid slice q reads: [lo, hi]
inline int64_t raw_lo(const Grid& G, const int64_t q) {
    return G.is_original[q] ? (int64_t)G.interp_index[q] : (int64_t)G.slices[G.interp_index[q]].zorg;
}
inline int64_t raw_hi(const Grid& G, const int64_t q) {
    return G.is_original[q] ? (int64_t)G.interp_index[q]
                            : (int64_t)G.slices[G.interp_index[q]].zorg + 1;
}

// n values to dst, a run of an output row: the whole 64-byte lines inside it with streaming
// stores (the output is written once), the partial lines at its ends with plain ones, so that
// no line takes both (another unit may write the rest of a partial line)
template <typename TO>
inline void write_segment(const TO* src, TO* dst, const uint64_t n) {
#if defined(__AVX2__)
    const uintptr_t a = reinterpret_cast<uintptr_t>(dst);
    const uintptr_t lo = (a + 63) & ~(uintptr_t)63, hi = (a + n * sizeof(TO)) & ~(uintptr_t)63;
    if (hi <= lo) {
        memcpy(dst, src, n * sizeof(TO));
        return;
    }
    const uint64_t i0 = (lo - a) / sizeof(TO), i1 = (hi - a) / sizeof(TO);
    memcpy(dst, src, i0 * sizeof(TO));
    constexpr uint64_t step = 32 / sizeof(TO);
    for (uint64_t i = i0; i < i1; i += step)
        _mm256_stream_si256((__m256i*)(dst + i), _mm256_loadu_si256((const __m256i*)(src + i)));
    memcpy(dst + i1, src + i1, (n - i1) * sizeof(TO));
#else
    memcpy(dst, src, n * sizeof(TO));
#endif
}

// rows [j0, j0 + n) of output plane `plane` (oy x ox) set to 0
template <typename TO>
inline void zero_rows(TO* plane, const uint64_t ox, const uint64_t j0, const uint64_t n) {
    if (n > 0) memset(plane + j0 * ox, 0, n * ox * sizeof(TO));
}

// columns [c0, c0 + 8) of rows [yrow, yrow + nl) of a raw slice (row stride nx) as single lanes:
// dst[(c - c0) * by + j]; columns past nx and lanes nl .. nb8 - 1 are 0
template <typename T>
inline void copy_columns(const T* plane, const uint64_t nx, const uint64_t yrow, const uint64_t nl,
                         const uint64_t nb8, const int64_t c0, float* dst, const uint64_t by) {
    const T* src = plane + yrow * nx + c0;
    const int64_t nc = std::min<int64_t>(8, (int64_t)nx - c0);
    uint64_t j = 0;
    if (nc == 8)
        for (; j + 8 <= nl; j += 8) block_to_lanes(src + j * nx, nx, dst + j, by);
    for (int64_t c = 0; c < 8; ++c) {
        float* d = dst + c * by;
        if (c >= nc) {
            for (uint64_t k = 0; k < nb8; ++k) d[k] = 0.0f;
            continue;
        }
        for (uint64_t k = j; k < nl; ++k) d[k] = (float)src[k * nx + c];
        for (uint64_t k = nl; k < nb8; ++k) d[k] = 0.0f;
    }
}

// ng groups of 8 lanes of an interpolated grid column from raw lanes: two lines (a, b) or four
// (a, b, c, d) with the slice's weights; mode4[g]: the 4-line sum of group g (0 SEQ, 1 PAIR, n + 2
// PAIR in the first n lanes, as pair_mode); round16: uint16 values (rounded)
inline void interp_lanes2(const float* a, const float* b, const V8 wa, const V8 wb, float* o,
                          const uint64_t ng, const bool round16) {
    for (uint64_t g = 0; g < ng; ++g) {
        V8 v = sum2(load8(a + 8 * g), load8(b + 8 * g), wa, wb);
        store8(o + 8 * g, round16 ? round_u16(v) : v);
    }
}

inline void interp_lanes4(const float* a, const float* b, const float* c, const float* d,
                          const V8 (&cw)[4], float* o, const uint64_t ng, const int* mode4,
                          const bool round16) {
    for (uint64_t g = 0; g < ng; ++g) {
        const V8 A = load8(a + 8 * g), B = load8(b + 8 * g), C = load8(c + 8 * g), D = load8(d + 8 * g);
        const int m = mode4[g];
        V8 v;
        if (m == 0)
            v = seq4(A, B, C, D, cw[0], cw[1], cw[2], cw[3]);
        else if (m == 1)
            v = pair4(A, B, C, D, cw[0], cw[1], cw[2], cw[3]);
        else
            v = first_lanes(pair4(A, B, C, D, cw[0], cw[1], cw[2], cw[3]),
                            seq4(A, B, C, D, cw[0], cw[1], cw[2], cw[3]), m - 2);
        store8(o + 8 * g, round16 ? round_u16(v) : v);
    }
}

// column r of interpolated grid slice p (the mex's three row ranges, its rows being columns
// here): raw(z, c) gives the lanes of raw slice z, column c
template <typename Raw>
inline void interp_column_lanes(const Raw& raw, const InterpSlice& p, const uint64_t r,
                                const uint64_t xa, const uint8_t Reverse, const uint64_t nx,
                                float* o, const uint64_t ng, const int* mode4, const bool round16) {
    const uint64_t zorg = p.zorg, s = p.s, t = p.t;
    const V8 cw[4] = {set8(p.c00_w), set8(p.c01_w), set8(p.c10_w), set8(p.c11_w)};
    if (Reverse) {
        if (r < xa - s)
            interp_lanes2(raw(zorg + 1, r + t), raw(zorg + 1, r + t + 1), cw[2], cw[3], o, ng, round16);
        else if (r < nx - t - 1)
            interp_lanes4(raw(zorg, r + s - xa), raw(zorg, r + s - xa + 1), raw(zorg + 1, r + t),
                          raw(zorg + 1, r + t + 1), cw, o, ng, mode4, round16);
        else
            interp_lanes2(raw(zorg, r + s - xa), raw(zorg, r + s - xa + 1), cw[0], cw[1], o, ng, round16);
    } else {
        if (r < xa - s)
            interp_lanes2(raw(zorg, r + s), raw(zorg, r + s + 1), cw[0], cw[1], o, ng, round16);
        else if (r < nx - t - 1)
            interp_lanes4(raw(zorg, r + s), raw(zorg, r + s + 1), raw(zorg + 1, r + t - xa),
                          raw(zorg + 1, r + t - xa + 1), cw, o, ng, mode4, round16);
        else
            interp_lanes2(raw(zorg + 1, r + t - xa), raw(zorg + 1, r + t - xa + 1), cw[2], cw[3], o, ng,
                          round16);
    }
}

#ifdef PK5D_DESKEW_NAN_X86
// the lanes of interp_column_lanes's column r (float data: the sums SEQ) that came out NaN,
// computed again from the same operands
template <typename Raw>
inline void interp_column_nans_x86(const Raw& raw, const InterpSlice& p, const uint64_t r,
                                   const uint64_t xa, const uint8_t Reverse, const uint64_t nx,
                                   float* o, const uint64_t ng) {
    const uint64_t zorg = p.zorg, s = p.s, t = p.t, n = 8 * ng;
    if (!any_nan(o, n)) return;
    if (Reverse) {
        if (r < xa - s)
            nans_x86_2(raw(zorg + 1, r + t), raw(zorg + 1, r + t + 1), o, p.c10_w, p.c11_w, n);
        else if (r < nx - t - 1)
            nans_x86_4(raw(zorg, r + s - xa), raw(zorg, r + s - xa + 1), raw(zorg + 1, r + t),
                       raw(zorg + 1, r + t + 1), o, p.c00_w, p.c01_w, p.c10_w, p.c11_w, n);
        else
            nans_x86_2(raw(zorg, r + s - xa), raw(zorg, r + s - xa + 1), o, p.c00_w, p.c01_w, n);
    } else {
        if (r < xa - s)
            nans_x86_2(raw(zorg, r + s), raw(zorg, r + s + 1), o, p.c00_w, p.c01_w, n);
        else if (r < nx - t - 1)
            nans_x86_4(raw(zorg, r + s), raw(zorg, r + s + 1), raw(zorg + 1, r + t - xa),
                       raw(zorg + 1, r + t - xa + 1), o, p.c00_w, p.c01_w, p.c10_w, p.c11_w, n);
        else
            nans_x86_2(raw(zorg + 1, r + t - xa), raw(zorg + 1, r + t - xa + 1), o, p.c10_w,
                       p.c11_w, n);
    }
}
#endif

// one output line: groups [0, np) with the formula fp, [np, ng) with ft (0 SEQ, 1 PAIR, 2
// SEQ_SWAP); half: + 0.5 for the uint16 output's rounding
template <typename TO>
inline void line_lanes(const float* a, const float* b, const float* c, const float* d, const V8 (&w)[4],
                       TO* o, const uint64_t np, const uint64_t ng, const int fp, const int ft,
                       const bool half) {
    for (uint64_t g = 0; g < ng; ++g) {
        const V8 A = load8(a + 8 * g), B = load8(b + 8 * g), C = load8(c + 8 * g), D = load8(d + 8 * g);
        const int f = g < np ? fp : ft;
        V8 v = f == 1 ? pair4(A, B, C, D, w[0], w[1], w[2], w[3])
                      : (f == 0 ? seq4(A, B, C, D, w[0], w[1], w[2], w[3])
                                : seq4_swap(A, B, C, D, w[0], w[1], w[2], w[3]));
        if (half) v = add_half(v);
        store8_any(o + 8 * g, v);
    }
}

// the first line xo in [x0, x1 + 1) of output slice z whose grid slice z0 is at least v (inc:
// z0 grows with xo) or below v (z0 falls with xo); the lines x0 .. x1 are inside
inline int64_t first_line(const Geometry& GE, const uint64_t z, const float wy, const int64_t x0,
                          const int64_t x1, const bool inc, const int64_t v) {
    int64_t l = x0, h = x1 + 1;
    while (l < h) {
        const int64_t m = l + (h - l) / 2;
        int64_t z0;
        float w[4];
        GE.line(z, (uint64_t)m, wy, z0, w);
        if (inc ? z0 >= v : z0 < v)
            h = m;
        else
            l = m + 1;
    }
    return l;
}

template <typename T, typename TO>
void blocked_warp(const T* const vol, const uint64_t nz, const uint64_t ny, const uint64_t nx,
                  const Grid& G, const bool interp, const uint8_t Reverse,
                  const double (&tmat)[4][4], TO* const out, const uint64_t oz, const uint64_t oy,
                  const uint64_t ox, const float start_y, const float start_x, const float start_z,
                  const uint64_t line_len, const bool fill, const int nthreads,
                  const uint64_t block_rows, const uint64_t window) {
    (void)nz;
    const uint64_t oyx = oy * ox;
    if (line_len == 0 || G.ntotal < 2) {
        if (fill) {
#pragma omp parallel for schedule(static) num_threads(nthreads)
            for (uint64_t z = 0; z < oz; ++z) zero_rows(out + z * oyx, ox, 0, oy);
        }
        return;
    }
    Geometry GE{{}, nx, G.ntotal, start_x, start_z};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) GE.tmat[i][j] = tmat[i][j];

    // per output slice: its column (-1: nothing to compute), wy and the run of lines inside
    std::vector<int64_t> col(oz, -1), xlo(oz, 0), xhi(oz, -1);
    std::vector<float> wys(oz, 0.0f);
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
    for (uint64_t z = 0; z < oz; ++z) {
        int64_t y0;
        float wy;
        if (!GE.slice(z, y0, wy)) continue;
        int64_t lo = -1, hi = -1;
        for (uint64_t xo = 0; xo < ox; ++xo) {
            int64_t z0;
            float w[4];
            if (GE.line(z, xo, wy, z0, w)) {
                if (lo < 0) lo = (int64_t)xo;
                hi = (int64_t)xo;
            }
        }
        if (lo < 0) continue;
        col[z] = y0;
        wys[z] = wy;
        xlo[z] = lo;
        xhi[z] = hi;
    }
    // (z0 = floor(tmat[2][1] * xo + ...): monotonic in xo)
    const bool inc = tmat[2][1] >= 0;
    if (fill) {
#pragma omp parallel for schedule(static) num_threads(nthreads)
        for (uint64_t z = 0; z < oz; ++z) {
            TO* plane = out + z * oyx;
            if (col[z] < 0) {
                zero_rows(plane, ox, 0, oy);
                continue;
            }
            for (uint64_t j = 0; j < line_len; ++j) {
                memset(plane + j * ox, 0, (uint64_t)xlo[z] * sizeof(TO));
                memset(plane + j * ox + xhi[z] + 1, 0, (ox - 1 - (uint64_t)xhi[z]) * sizeof(TO));
            }
            zero_rows(plane, ox, line_len, oy - line_len);
        }
    }
    // the slices that compute lines, in the order of their columns
    std::vector<uint64_t> order;
    for (uint64_t z = 0; z < oz; ++z)
        if (col[z] >= 0) order.push_back(z);
    if (order.empty()) return;
    if (tmat[1][2] < 0) std::reverse(order.begin(), order.end());

    // the ring: a power of two of 8-column groups, at least the groups of one slice
    const int64_t xa = interp ? (int64_t)G.xa : 0;
    const int64_t margin = xa;  // interpolated columns r read raw columns r - xa .. r + xa
    // (the ring holds chunks of CHUNK columns: a chunk is copied whole, so that each cache line
    // of the input is read once)
    int64_t gspan = 1;
    for (const uint64_t z : order) {
        const int64_t lo = std::max<int64_t>(0, col[z] - margin);
        const int64_t hi = std::min<int64_t>((int64_t)nx - 1, col[z] + 1 + margin);
        gspan = std::max(gspan, hi / CHUNK - lo / CHUNK + 1);
    }
    int64_t rg = 1;
    while (rg < gspan) rg *= 2;
    const uint64_t rw = CHUNK * (uint64_t)rg;

    // rows per block; windows of grid slices z0 whose raw columns fit RING_BYTES
    uint64_t by = block_rows > 0 ? std::min<uint64_t>(block_rows / 8 * 8, MAX_BLOCK_ROWS) : BLOCK_ROWS;
    if (by == 0) by = 8;
    by = std::min<uint64_t>(by, (line_len + 7) / 8 * 8);
    const uint64_t nblocks = (line_len + by - 1) / by;
    // lane strides of the scratch (padded: no power-of-two strides, which make the loads of the
    // lines alias earlier stores and the copies share a few cache sets)
    const uint64_t ls = by + LANE_PAD, os = by + 8;
    const uint64_t ngrid = G.ntotal - 1;  // z0 is at most ntotal - 2
    const double raw_per_grid =
        (double)(raw_hi(G, (int64_t)G.ntotal - 1) + 1) / (double)G.ntotal;
    uint64_t win = window;
    if (win == 0) {
        const double per = (double)(rw * (by + LANE_PAD) * sizeof(float)) * raw_per_grid;
        win = (uint64_t)std::max(8.0, std::floor((double)RING_BYTES / per));
        // enough units for the threads
        const uint64_t want = 4 * (uint64_t)nthreads;
        if (nblocks * ((ngrid + win - 1) / win) < want)
            win = std::max<uint64_t>(4, (ngrid * nblocks + want - 1) / want);
    }
    win = std::min<uint64_t>(std::max<uint64_t>(win, 1), ngrid);
    const uint64_t nwin = (ngrid + win - 1) / win;
    const uint64_t nunits = nblocks * nwin;
    // per slice (in order), the line runs of the windows: window w has the lines of z0 in
    // [w win, (w + 1) win), xo in [bnd[w], bnd[w + 1]) (or [bnd[w + 1], bnd[w]) where z0 falls
    // with xo: z0 moves monotonically with xo)
    const uint64_t nb = nwin + 1;
    std::vector<int64_t> bnd(order.size() * nb);
#pragma omp parallel for schedule(dynamic, 4) num_threads(nthreads)
    for (uint64_t k = 0; k < order.size(); ++k) {
        const uint64_t z = order[k];
        for (uint64_t w = 0; w < nb; ++w)
            bnd[k * nb + w] =
                first_line(GE, z, wys[z], xlo[z], xhi[z], inc, (int64_t)std::min(ngrid, w * win));
    }
    // the most raw and interpolated slices a window holds
    uint64_t max_raw = 1, max_int = 1;
    for (uint64_t w = 0; w < nwin; ++w) {
        const int64_t ga = (int64_t)(w * win), gb = (int64_t)std::min(ngrid, (w + 1) * win);
        max_raw = std::max<uint64_t>(max_raw, (uint64_t)(raw_hi(G, gb) - raw_lo(G, ga) + 1));
        uint64_t ni = 0;
        for (int64_t q = ga; q <= gb; ++q) ni += G.is_original[q] ? 0 : 1;
        max_int = std::max(max_int, ni);
    }

    constexpr int kind = sizeof(T) == 4 ? 0 : (sizeof(TO) == 4 ? 1 : 2);
    const uint64_t y0 = (uint64_t)floor(start_y);
    const uint64_t len16 = line_len - (line_len % 16);
    const uint64_t ny16 = ny - (ny % 16);
    const uint64_t nyx = ny * nx;
    const uint64_t* idx = G.interp_index.data();
    const bool* orig = G.is_original.get();
    const uint64_t seg = std::min<uint64_t>(SEG_LINES, ox);

#pragma omp parallel num_threads(nthreads)
    {
        std::vector<float> R(max_raw * rw * ls), I(max_int * 2 * ls);
        std::vector<int64_t> ring_tag(rg);
        std::vector<const float*> P(2 * (win + 1));
        std::vector<TO> O(seg * os), RT(seg * by);

#pragma omp for schedule(dynamic, 1)
        for (uint64_t u = 0; u < nunits; ++u) {
            struct Fence {
                ~Fence() { stream_fence(); }
            } fence;
            const uint64_t b = u % nblocks, wi = u / nblocks;
            const uint64_t i0 = b * by;
            const uint64_t nl = std::min<uint64_t>(by, line_len - i0);
            const uint64_t nb8 = (nl + 7) / 8 * 8, ng = nb8 / 8;
            const uint64_t yrow = y0 + i0;
            const int64_t ga = (int64_t)(wi * win), gb = (int64_t)std::min(ngrid, (wi + 1) * win);
            const int64_t ulo = raw_lo(G, ga), uhi = raw_hi(G, gb);
            // interpolated grid slices of the window: [ilo, ihi] (none: ihi < ilo)
            int64_t ilo = INT64_MAX, ihi = -1;
            for (int64_t q = ga; q <= gb; ++q)
                if (!orig[q]) {
                    ilo = std::min<int64_t>(ilo, (int64_t)idx[q]);
                    ihi = std::max<int64_t>(ihi, (int64_t)idx[q]);
                }
            // the 4-line mode of the interpolated values per group (absolute rows) and the
            // groups of output lines below floor(line_len / 16) * 16 (PAIR for uint16)
            int mode4[MAX_BLOCK_ROWS / 8];
            for (uint64_t g = 0; g < ng; ++g) mode4[g] = kind == 0 ? 0 : pair_mode(yrow + 8 * g, ny16);
            const uint64_t np = i0 >= len16 ? 0 : std::min<uint64_t>(ng, (len16 - i0) / 8);
            std::fill(ring_tag.begin(), ring_tag.end(), -1);
            auto rcol = [&](uint64_t z, uint64_t c) -> const float* {
                return R.data() + ((z - (uint64_t)ulo) * rw + (c & (rw - 1))) * ls;
            };

            for (uint64_t si = 0; si < order.size(); ++si) {
                const uint64_t zo = order[si];
                const int64_t c = col[zo];
                const int64_t a = bnd[si * nb + (inc ? wi : wi + 1)];
                const int64_t e = bnd[si * nb + (inc ? wi + 1 : wi)];
                if (a >= e) continue;
                // the slice's columns in the ring
                const int64_t clo = std::max<int64_t>(0, c - margin);
                const int64_t chi = std::min<int64_t>((int64_t)nx - 1, c + 1 + margin);
                for (int64_t k = clo / CHUNK; k <= chi / CHUNK; ++k) {
                    if (ring_tag[k & (rg - 1)] == k) continue;
                    ring_tag[k & (rg - 1)] = k;
                    const int64_t g1 = std::min<int64_t>((k + 1) * (CHUNK / 8), ((int64_t)nx + 7) / 8);
                    for (int64_t z = ulo; z <= uhi; ++z) {
                        // the rows of the next raw slice on their way (strided rows: no hardware
                        // prefetch)
                        if (PREFETCH && z + PREFETCH <= uhi) {
                            const T* pz = vol + (uint64_t)(z + PREFETCH) * nyx + yrow * nx + k * CHUNK;
                            for (uint64_t j = 0; j < nl; ++j) prefetch_line(pz + j * nx);
                        }
                        float* ring = R.data() + (uint64_t)(z - ulo) * rw * ls;
                        for (int64_t g = k * (CHUNK / 8); g < g1; ++g)
                            copy_columns(vol + (uint64_t)z * nyx, nx, yrow, nl, nb8, 8 * g,
                                         ring + ((uint64_t)(8 * g) & (rw - 1)) * ls, ls);
                    }
                }
                // the grid slices the run's lines read: [qa, qb]
                int64_t qa, qb;
                {
                    int64_t za, zb;
                    float wt[4];
                    GE.line(zo, (uint64_t)a, wys[zo], za, wt);
                    GE.line(zo, (uint64_t)(e - 1), wys[zo], zb, wt);
                    qa = std::min(za, zb);
                    qb = std::max(za, zb) + 1;
                }
                // their interpolated columns c and c + 1
                for (int64_t q = qa; q <= qb; ++q) {
                    if (orig[q]) continue;
                    const int64_t ii = (int64_t)idx[q];
                    for (int k = 0; k < 2; ++k) {
                        interp_column_lanes(rcol, G.slices[ii], (uint64_t)(c + k), G.xa, Reverse, nx,
                                            I.data() + ((uint64_t)(ii - ilo) * 2 + k) * ls, ng, mode4,
                                            kind == 2);
#ifdef PK5D_DESKEW_NAN_X86
                        if constexpr (kind == 0)
                            interp_column_nans_x86(rcol, G.slices[ii], (uint64_t)(c + k), G.xa,
                                                   Reverse, nx,
                                                   I.data() + ((uint64_t)(ii - ilo) * 2 + k) * ls, ng);
#endif
                    }
                }
                // grid slice q, column c + k: raw or interpolated lanes
                for (int64_t q = qa; q <= qb; ++q)
                    for (int k = 0; k < 2; ++k)
                        P[2 * (q - ga) + k] =
                            orig[q] ? rcol(idx[q], (uint64_t)(c + k))
                                    : I.data() + ((uint64_t)((int64_t)idx[q] - ilo) * 2 + k) * ls;
                for (int64_t s0 = a; s0 < e; s0 += (int64_t)seg) {
                    const int64_t s1 = std::min<int64_t>(e, s0 + (int64_t)seg);
                    const uint64_t ns = (uint64_t)(s1 - s0);
                    for (int64_t xo = s0; xo < s1; ++xo) {
                        int64_t z0;
                        float wt[4];
                        GE.line(zo, (uint64_t)xo, wys[zo], z0, wt);
                        const V8 w[4] = {set8(wt[0]), set8(wt[1]), set8(wt[2]), set8(wt[3])};
                        // the line sums the mex picks for the type and, for uint16 to single,
                        // the grid slices' kinds (interpolated values are single there)
                        int fp = 0, ft = 0;
                        if (kind == 2) {
                            fp = 1;
                        } else if (kind == 1) {
                            const bool oa = orig[z0], ob = orig[z0 + 1];
                            if (oa || ob) {
                                fp = 1;
                                ft = oa ? 0 : 2;
                            }
                        }
                        // (z0, col), (z0 + 1, col), (z0, col + 1), (z0 + 1, col + 1)
                        const float* const* p = P.data() + 2 * (z0 - ga);
                        TO* o = O.data() + (uint64_t)(xo - s0) * os;
                        line_lanes(p[0], p[2], p[1], p[3], w, o, np, ng, fp, ft, kind == 2);
#ifdef PK5D_DESKEW_NAN_X86
                        if constexpr (kind == 0)
                            nans_x86_4(p[0], p[2], p[1], p[3], o, wt[0], wt[1], wt[2], wt[3], 8 * ng);
#endif
                    }
                    // the lines to rows i0 .. i0 + nl - 1, columns s0 .. s1 - 1: transposed, streamed
                    uint64_t j8 = 0;
                    for (; j8 + 8 <= nl; j8 += 8) {
                        uint64_t x8 = 0;
                        for (; x8 + 8 <= ns; x8 += 8)
                            transpose_block8(O.data() + x8 * os + j8, os, RT.data() + j8 * ns + x8, ns);
                        for (; x8 < ns; ++x8)
                            for (uint64_t j = j8; j < j8 + 8; ++j) RT[j * ns + x8] = O[x8 * os + j];
                    }
                    for (; j8 < nl; ++j8)
                        for (uint64_t x = 0; x < ns; ++x) RT[j8 * ns + x] = O[x * os + j8];
                    TO* dst = out + zo * oyx + i0 * ox + (uint64_t)s0;
                    for (uint64_t j = 0; j < nl; ++j) {
                        if (PREFETCH && j + 4 < nl) prefetch_write(dst + (j + 4) * ox, ns * sizeof(TO));
                        write_segment(RT.data() + j * ns, dst + j * ox, ns);
                    }
                }
            }
        }
    }
}

// the fused kernel (skewed_space_interp_volume_deskew_rotate_warp_mex) through blocked_warp;
// block_rows and window 0: chosen here
template <typename T, typename TO>
void skewed_space_interp_warp(const T* const vol, const uint64_t nz, const uint64_t ny,
                                      const uint64_t nx, const double xstep, const double stepsize,
                                      const uint8_t Reverse, const double (&tmat)[4][4],
                                      TO* const out, const uint64_t oz, const uint64_t oy,
                                      const uint64_t ox, const float start_y, const float start_x,
                                      const float start_z, const uint64_t line_len, const bool fill,
                                      const int nthreads, const uint64_t block_rows,
                                      const uint64_t window) {
    Grid G;
    make_grid(nz, xstep, stepsize, Reverse, G);
    blocked_warp<T, TO>(vol, nz, ny, nx, G, true, Reverse, tmat, out, oz, oy, ox, start_y,
                        start_x, start_z, line_len, fill, nthreads, block_rows, window);
}

}  // namespace pk5d::deskew
