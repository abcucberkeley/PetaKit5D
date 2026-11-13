// decon_mask_edge_erosion_mex on (z, y, x) masks, with the mex's result (that of MATLAB's imerode): each z
// slice of a mask is eroded by a sequence of 2D neighborhoods as imerode applies a decomposed structuring
// element (e.g., the periodic lines of strel('disk', r)): the slice is padded with true by pad_ul and pad_lr
// (morphop's padding for decompositions), eroded by each neighborhood in turn with pixels outside the padded
// slice counting as true, and cropped back. Then the mask is eroded along z by a line of 2 * rz + 1 voxels
// (voxels outside the volume count as true), as imerode with true(2 * rz + 1, 1) along z. The erosion by a
// neighborhood is B(y, x) = all of A(y + i - cy, x + j - cx) over its nonzero elements (i, j), with the
// center (cy, cx) = floor((size + 1) / 2), as imerode's. With clear_border the first and last rows and
// columns (3D: and slices) are false first.
//
// An AND does not depend on the order it is taken in, so the slices are eroded row by row along x (the
// mex: column by column along y), and the threads split whole slices or the rows of a slice.
#include <omp.h>

#include <algorithm>
#include <cstring>
#include <utility>
#include <vector>

#include "decon.hpp"

namespace pk5d::decon {
namespace {

struct Geom {
    const ErosionPlan* p;
    idx_t h, w, d;          // mask size: rows (y), columns (x), slices
    idx_t gr, gc;           // guard rows / columns around the padded slice (true, never written)
    idx_t hp, wp, hg, wg;   // padded slice and buffer size
    int ns;
};

// buffer element of row r, column c of the padded slice (row-major)
inline idx_t bidx(const Geom& g, idx_t r, idx_t c) { return (r + g.gr) * g.wg + c + g.gc; }

// rows r0 .. r1 - 1 of the padded slice from slice z of A (clear: the whole slice false)
void load_rows(const Geom& g, const uint8_t* src, bool clear, uint8_t* a, idx_t r0, idx_t r1) {
    const ErosionPlan& p = *g.p;
    for (idx_t r = r0; r < r1; r++) {
        uint8_t* row = a + bidx(g, r, 0);
        const idx_t y = r - p.pu[0];
        if (y < 0 || y >= g.h) {
            std::memset(row, 1, static_cast<size_t>(g.wp));
            continue;
        }
        std::memset(row, 1, static_cast<size_t>(p.pu[1]));
        std::memset(row + p.pu[1] + g.w, 1, static_cast<size_t>(p.pl[1]));
        uint8_t* cc = row + p.pu[1];
        if (clear || (p.clear_border && (y == 0 || y == g.h - 1))) {
            std::memset(cc, 0, static_cast<size_t>(g.w));
            continue;
        }
        const uint8_t* s = src + y * g.w;
        for (idx_t x = 0; x < g.w; x++) cc[x] = s[x] != 0;
        if (p.clear_border) {
            cc[0] = 0;
            cc[g.w - 1] = 0;
        }
    }
}

// erosion by neighborhood k of rows r0 .. r1 - 1 of the padded slice, from a to b
void erode_rows(const Geom& g, int k, const uint8_t* a, uint8_t* b, idx_t r0, idx_t r1) {
    const ErosionPlan& p = *g.p;
    const int i0 = p.start[k], n = p.start[k + 1] - i0;
    const idx_t wp = g.wp;
    for (idx_t r = r0; r < r1; r++) {
        uint8_t* __restrict o = b + bidx(g, r, 0);
        std::memcpy(o, a + bidx(g, r + p.dr[i0], p.dc[i0]), static_cast<size_t>(wp));
        for (int i = 1; i < n; i++) {
            const uint8_t* __restrict s = a + bidx(g, r + p.dr[i0 + i], p.dc[i0 + i]);
            for (idx_t c = 0; c < wp; c++) o[c] &= s[c];
        }
    }
}

// rows y0 .. y1 - 1 of the cropped slice from the padded slice a to dst
void store_rows(const Geom& g, const uint8_t* a, uint8_t* dst, idx_t y0, idx_t y1) {
    for (idx_t y = y0; y < y1; y++) {
        std::memcpy(dst + y * g.w, a + bidx(g, g.p->pu[0] + y, g.p->pu[1]), static_cast<size_t>(g.w));
    }
}

bool cleared(const Geom& g, idx_t z) { return g.p->clear_border && g.d > 1 && (z == 0 || z == g.d - 1); }

// the 2D erosions of slice z from A to B with the buffers a and b (guards true); serial
void erode_slice(const Geom& g, const uint8_t* A, uint8_t* B, idx_t z, uint8_t* a, uint8_t* b) {
    const idx_t n2 = g.h * g.w;
    load_rows(g, A + z * n2, cleared(g, z), a, 0, g.hp);
    for (int k = 0; k < g.ns; k++) {
        erode_rows(g, k, a, b, 0, g.hp);
        std::swap(a, b);
    }
    store_rows(g, a, B + z * n2, 0, g.h);
}

// the same with the rows split among the threads (for few slices)
void erode_slice_parallel(const Geom& g, const uint8_t* A, uint8_t* B, idx_t z, uint8_t* a, uint8_t* b, int nt) {
    const idx_t n2 = g.h * g.w;
    const bool clear = cleared(g, z);
#pragma omp parallel num_threads(nt)
    {
        uint8_t* ta = a;
        uint8_t* tb = b;
#pragma omp for schedule(static)
        for (idx_t r = 0; r < g.hp; r++) load_rows(g, A + z * n2, clear, ta, r, r + 1);
        for (int k = 0; k < g.ns; k++) {
#pragma omp for schedule(static)
            for (idx_t r = 0; r < g.hp; r++) erode_rows(g, k, ta, tb, r, r + 1);
            std::swap(ta, tb);
        }
#pragma omp for schedule(static)
        for (idx_t y = 0; y < g.h; y++) store_rows(g, ta, B + z * n2, y, y + 1);
    }
}

// erosion along z by 2 * rz + 1 voxels, in place, with per-thread copies of blocks of voxel columns; C counts
// the false voxels in the window of each voxel column
template <typename C>
void erode_z(const Geom& g, uint8_t* B, idx_t rz, int nt) {
    const idx_t n2 = g.h * g.w, d = g.d;
    const idx_t bs = 1024;
    const idx_t nb = (n2 + bs - 1) / bs;
    std::vector<std::vector<uint8_t>> blk(static_cast<size_t>(nt), std::vector<uint8_t>(static_cast<size_t>(bs * d)));
    std::vector<std::vector<C>> cnt(static_cast<size_t>(nt), std::vector<C>(static_cast<size_t>(bs)));
#pragma omp parallel for num_threads(nt) schedule(dynamic, 4)
    for (idx_t ib = 0; ib < nb; ib++) {
        const int t = omp_get_thread_num();
        const idx_t p0 = ib * bs;
        const idx_t L = n2 - p0 < bs ? n2 - p0 : bs;
        uint8_t* __restrict v = blk[t].data();
        C* __restrict c = cnt[t].data();
        for (idx_t z = 0; z < d; z++) std::memcpy(v + z * L, B + z * n2 + p0, static_cast<size_t>(L));
        for (idx_t i = 0; i < L; i++) c[i] = 0;
        const idx_t zi = rz < d - 1 ? rz : d - 1;
        for (idx_t z = 0; z <= zi; z++) {
            const uint8_t* __restrict s = v + z * L;
            for (idx_t i = 0; i < L; i++) c[i] += static_cast<C>(1 - s[i]);
        }
        for (idx_t z = 0; z < d; z++) {
            uint8_t* __restrict o = B + z * n2 + p0;
            for (idx_t i = 0; i < L; i++) o[i] = c[i] == 0;
            if (z + rz + 1 < d) {
                const uint8_t* __restrict s = v + (z + rz + 1) * L;
                for (idx_t i = 0; i < L; i++) c[i] += static_cast<C>(1 - s[i]);
            }
            if (z - rz >= 0) {
                const uint8_t* __restrict s = v + (z - rz) * L;
                for (idx_t i = 0; i < L; i++) c[i] -= static_cast<C>(1 - s[i]);
            }
        }
    }
}

}  // namespace

void mask_edge_erosion(const uint8_t* A, uint8_t* B, idx_t nz, idx_t ny, idx_t nx, const ErosionPlan& p, int nt) {
    Geom g;
    g.p = &p;
    g.h = ny;
    g.w = nx;
    g.d = nz;
    g.ns = static_cast<int>(p.start.size()) - 1;
    g.gr = 0;
    g.gc = 0;
    for (size_t i = 0; i < p.dr.size(); i++) {
        g.gr = std::max<idx_t>(g.gr, p.dr[i] < 0 ? -p.dr[i] : p.dr[i]);
        g.gc = std::max<idx_t>(g.gc, p.dc[i] < 0 ? -p.dc[i] : p.dc[i]);
    }
    const idx_t n2 = g.h * g.w;
    if (n2 * g.d == 0) return;
    g.hp = g.h + p.pu[0] + p.pl[0];
    g.wp = g.w + p.pu[1] + p.pl[1];
    g.hg = g.hp + 2 * g.gr;
    g.wg = g.wp + 2 * g.gc;
    const size_t bsz = static_cast<size_t>(g.hg * g.wg);
    // slices: per thread with whole slices when there are enough of them, otherwise the rows split
    const bool per_thread = nt == 1 || (g.d >= 2 * static_cast<idx_t>(nt) && 2 * bsz * static_cast<size_t>(nt) <= (static_cast<size_t>(1) << 32));
    if (per_thread) {
        std::vector<std::vector<uint8_t>> bufs(static_cast<size_t>(2 * nt), std::vector<uint8_t>(bsz, 1));
#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
        for (idx_t z = 0; z < g.d; z++) {
            const int t = omp_get_thread_num();
            erode_slice(g, A, B, z, bufs[2 * t].data(), bufs[2 * t + 1].data());
        }
    } else {
        std::vector<uint8_t> a(bsz, 1), b(bsz, 1);
        for (idx_t z = 0; z < g.d; z++) erode_slice_parallel(g, A, B, z, a.data(), b.data(), nt);
    }
    if (g.d > 1 && p.rz > 0) {
        const idx_t wmax = 2 * p.rz + 1 < g.d ? 2 * p.rz + 1 : g.d;
        if (wmax <= 255) {
            erode_z<uint8_t>(g, B, p.rz, nt);
        } else if (wmax <= 65535) {
            erode_z<uint16_t>(g, B, p.rz, nt);
        } else {
            erode_z<uint32_t>(g, B, p.rz, nt);
        }
    }
}

}  // namespace pk5d::decon
