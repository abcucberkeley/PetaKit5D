// Output geometry of the deskew/rotate warps, and the work units of the C-order kernels.
//
// Names follow the mex, whose arrays are MATLAB's (y, x, z): the mex's "x" is MATLAB's y (the
// lines it sums), its "y" MATLAB's x, its "z" z. Here volumes are C-ordered (z, y, x), so the
// mex's lines are y columns of stride nx. The geometry is the mex's: output slice zo reads input
// columns x = col, col + 1 (col from tmat row 2, the nudge at the last column included); output
// line (zo, xo) reads slices z0, z0 + 1 of the (interpolation) grid (tmat row 3, the nudge at
// the last slice included), with the weights w00 = (1 - wy)(1 - wz), w01 = (1 - wy) wz,
// w10 = wy (1 - wz), w11 = wy wz in single precision. tmat keeps MATLAB's convention.
//
// Work: a unit is a block of 16 rows y of the lines (blocks start at multiples of 16 of the
// line, so floor(len / 16) * 16 never splits one). It goes through groups of consecutive output
// slices whose columns lie within a few columns: copies the input values a group needs into a
// scratch with the 16 rows as lanes (8 x 8 transposes; the next group reuses the cache lines),
// computes 16 voxels of a line per step with the mex's arithmetic (the position of a voxel in
// its line decides it, lanes.hpp), and writes the 16 rows of each output slice back through
// 8 x 8 transposes.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "lanes.hpp"

namespace pk5d::deskew {

// The mex's geometry. tmat rows 2 and 3 give, for output slice z (start startZ), the input
// column ty and, for its line y (start startY), the grid slice tz; the nudges keep the last
// column and the last slice inside.
struct Geometry {
    double tmat[4][4];
    uint64_t shapeY, shapeZ;  // input columns, grid slices
    float startY, startZ;

    // column of slice z and its weight wy; false: the slice is outside
    bool slice(const uint64_t z, int64_t& y0, float& wy) const {
        double ty;
        ty = tmat[1][2] * (z + startZ) + tmat[1][3];
        if (ty < 0 || ty + 1 >= shapeY) {
            if (ty + 1 == shapeY && shapeY > 1) {  // on the last row (a single row: outside)
                uint64_t counter = 1;
                while (floor(ty) + 1 == shapeY) {
                    ty -= 1e-9 * counter;
                    counter *= 2;
                }
            } else
                return false;
        }
        y0 = floor(ty);
        wy = ty - y0;
        return true;
    }

    // grid slice of line y of slice z and the weights w00, w01, w10, w11; false: outside
    bool line(const uint64_t z, const uint64_t y, const float wy, int64_t& z0, float (&w)[4]) const {
        double tz;
        tz = tmat[2][1] * (y + startY) + tmat[2][2] * (z + startZ) + tmat[2][3];
        if (tz < 0 || tz + 1 >= shapeZ) {
            // handling of edge cases
            if (tz + 1 == shapeZ && shapeZ > 1) {  // on the last slice (a single slice: outside)
                uint64_t counter = 1;
                while (floor(tz) + 1 == shapeZ) {
                    tz -= 1e-9 * counter;
                    counter *= 2;
                }
            } else
                return false;
        }
        z0 = floor(tz);
        float wz = tz - z0;
        w[0] = (1 - wy) * (1 - wz);
        w[1] = (1 - wy) * wz;
        w[2] = wy * (1 - wz);
        w[3] = wy * wz;
        return true;
    }
};

// per output slice: its column (-1: nothing to compute), wy, and the grid slices its lines
// read ([zlo, zhi]: z0 and z0 + 1)
struct Slices {
    std::vector<int64_t> col, zlo, zhi;
    std::vector<float> wy;
};

inline void compute_slices(const Geometry& G, const uint64_t ox, const uint64_t oz, Slices& S,
                           const int nthreads) {
    S.col.assign(oz, -1);
    S.zlo.assign(oz, 0);
    S.zhi.assign(oz, -1);
    S.wy.assign(oz, 0.0f);
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
    for (uint64_t z = 0; z < oz; ++z) {
        int64_t y0;
        float wy;
        if (!G.slice(z, y0, wy)) continue;
        int64_t lo = INT64_MAX, hi = -1;
        for (uint64_t y = 0; y < ox; ++y) {
            int64_t z0;
            float w[4];
            if (G.line(z, y, wy, z0, w)) {
                lo = std::min(lo, z0);
                hi = std::max(hi, z0 + 1);
            }
        }
        if (hi < 0) continue;
        S.col[z] = y0;
        S.wy[z] = wy;
        S.zlo[z] = lo;
        S.zhi[z] = hi;
    }
}

// consecutive output slices with lines inside, whose columns span at most max_span
struct Group {
    uint64_t zo_begin, zo_end;  // output slices [begin, end)
    int64_t cmin, cmax;         // their columns
    int64_t glo, ghi;           // grid slices their lines read: [glo, ghi]
};

inline std::vector<Group> make_groups(const Slices& S, const int64_t max_span, const uint64_t max_slices) {
    std::vector<Group> groups;
    bool open = false;
    Group g{};
    for (uint64_t z = 0; z < S.col.size(); ++z) {
        const int64_t c = S.col[z];
        if (c < 0) {  // nothing to compute in this slice
            if (open) groups.push_back(g);
            open = false;
            continue;
        }
        if (open && (z - g.zo_begin) < max_slices && std::max(g.cmax, c) - std::min(g.cmin, c) <= max_span) {
            g.zo_end = z + 1;
            g.cmin = std::min(g.cmin, c);
            g.cmax = std::max(g.cmax, c);
            g.glo = std::min(g.glo, S.zlo[z]);
            g.ghi = std::max(g.ghi, S.zhi[z]);
            continue;
        }
        if (open) groups.push_back(g);
        g = Group{z, z + 1, c, c, S.zlo[z], S.zhi[z]};
        open = true;
    }
    if (open) groups.push_back(g);
    return groups;
}

// the lines of a group's slices: grid slice z0 (-1: outside) and the 4 weights of each
// (computed once per unit, used by each of its blocks of rows)
struct GroupLines {
    std::vector<int64_t> z0;
    std::vector<float> w;
};

inline void group_lines(const Geometry& G, const Slices& S, const Group& g, const uint64_t ox,
                        GroupLines& GL) {
    const uint64_t ns = g.zo_end - g.zo_begin;
    GL.z0.assign(ns * ox, -1);
    GL.w.resize(ns * ox * 4);
    for (uint64_t s = 0; s < ns; ++s) {
        const uint64_t zo = g.zo_begin + s;
        if (S.col[zo] < 0) continue;
        for (uint64_t xo = 0; xo < ox; ++xo) {
            int64_t z0;
            float w[4];
            const uint64_t k = s * ox + xo;
            if (!G.line(zo, xo, S.wy[zo], z0, w)) continue;
            GL.z0[k] = z0;
            for (int q = 0; q < 4; ++q) GL.w[4 * k + q] = w[q];
        }
    }
}

// blocks of 16 rows per unit
constexpr uint64_t UNIT_BLOCKS = 32;

// the columns a group's units copy: [cmin, cmax + 1], widened by `extra` on both sides, within
// the input
inline void group_columns(const Group& g, const int64_t extra, const int64_t nx, int64_t& lo, int64_t& hi) {
    lo = std::max<int64_t>(0, g.cmin - extra);
    hi = std::min<int64_t>(nx - 1, g.cmax + 1 + extra);
    // widened to a multiple of 8 columns inside the input, so that the copies are 8 x 8
    // transposes
    const int64_t w8 = (hi - lo + 1 + 7) / 8 * 8;
    if (w8 <= nx) {
        if (lo + w8 - 1 <= nx - 1)
            hi = lo + w8 - 1;
        else
            lo = nx - w8, hi = nx - 1;
    }
}

// lanes per step: 16 rows y of a line (two vectors of 8)
constexpr int LANES = 16;

// copy rows y0 .. y0 + nl - 1 (nl <= 16) of columns [clo, chi] of plane `plane` (row stride
// nx) into dst[(c - clo) * 16 + j] (lanes past nl are 0)
template <typename T>
inline void copy_lanes(const T* plane, const uint64_t nx, const uint64_t y0, const int nl,
                       const int64_t clo, const int64_t chi, float* dst) {
    const int64_t w = chi - clo + 1;
    const T* p = plane + y0 * nx + clo;
    int64_t c = 0;
    if (nl == LANES) {
        for (; c + 8 <= w; c += 8) {
            block_to_lanes(p + c, nx, dst + c * LANES, LANES);
            block_to_lanes(p + 8 * nx + c, nx, dst + c * LANES + 8, LANES);
        }
    }
    for (; c < w; ++c) {
        float* d = dst + c * LANES;
        for (int j = 0; j < LANES; ++j) d[j] = j < nl ? (float)p[j * nx + c] : 0.0f;
    }
}

// write lanes src[xo * 16 + j] (ox lines) to rows j < nl of an output plane (row stride ox):
// transposed into the row tile R (nl x ox) first, then each row streamed to the output
template <typename TO>
inline void write_lanes(const float* src, const uint64_t ox, const int nl, TO* rows, TO* R) {
    uint64_t xo = 0;
    if (nl == LANES) {
        for (; xo + 8 <= ox; xo += 8) {
            lanes_to_block(src + xo * LANES, LANES, R + xo, ox);
            lanes_to_block(src + xo * LANES + 8, LANES, R + 8 * ox + xo, ox);
        }
    }
    for (; xo < ox; ++xo) {
        const float* s = src + xo * LANES;
        for (int j = 0; j < nl; ++j) R[j * ox + xo] = to_out<TO>(s[j]);
    }
    stream_rows(R, rows, (uint64_t)nl * ox);
}

}  // namespace pk5d::deskew
