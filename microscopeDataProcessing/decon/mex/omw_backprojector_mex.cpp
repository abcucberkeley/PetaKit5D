// Steps of the OMW back projector generation (omw_backprojector_generation) with the same results as the
// MATLAB code, bit for bit.
//
// B = omw_backprojector_mex('convhull', A, zinds, yinds, xinds, nThreads)
//   Replaces the slices A(:, :, z) for z in zinds, then squeeze(A(y, :, :)) for y in yinds, then
//   squeeze(A(:, x, :)) for x in xinds (1-based; the slices of each step in parallel) of the 3D logical array A
//   (all dimensions > 1) with their convex hull images, as bwconvhull(slice, 'union') makes them: the convex
//   hull, with all the points on its boundary, of the midpoints of the pixel edges of the perimeter pixels,
//   counterclockwise in (row, column) as convhull returns it, filled within the bounding box of the pixels
//   with the algorithm of poly2mask (vertices on a 5x5 subpixel grid, edge list, parity scan; the portable
//   version in toolbox/images/images/eml/poly2mask.m), the rest of the slice false. Only the edge midpoints
//   that can be on the hull boundary are used (of the first and last pixels of each row, and all of the first
//   and last row), which gives the same boundary points.
//
// [mask, inds] = omw_backprojector_mex('window', psz, E, R, bounds, nThreads)
//   For the voxels of a psz grid with center (psz + 1) / 2, decides where the apodization window of
//   omw_backprojector_generation is 1 or 0 whatever the nearest edge points are: the window of a voxel at
//   radius r is (r / rn < l) + taper (l <= r / rn < u), where rn is a convex combination of the radii R of
//   the 2 edge points nearest to the voxel's direction (E: their unit vectors, [y, x, z] as rows). The unit
//   vectors are binned in cells of 1/64 of [-1, 1]^3; for each cell, the 2 nearest edge points of any
//   direction in the cell are within the distance of the 2nd nearest edge point from the cell center plus
//   twice the cell's half diagonal, so the min and max of R over the edge points within that distance bound
//   rn. mask is 1 for voxels with r < l * min R (with a margin of 1e-9 for rounding), 0 elsewhere; inds are
//   the 1-based linear indices of the voxels in between (r < u * max R), whose window has to be computed.
//   bounds = [l, u] with 0 < l < u.
//
// B = omw_backprojector_mex('components', A, mode, n, nThreads)
//   The 26-connected components of the 3D logical array A that are kept by bwareaopen(A, n) (mode 'minsize':
//   the components with at least n voxels) or that are the n largest ones (mode 'largest'; ties go to the
//   component first in bwconncomp's order, that of their first voxels, as max and maxk choose).
//
// v = omw_backprojector_mex('sortdesc', a, nThreads)
//   sort(double(a(:)), 'descend') for single a without NaN or negative values (radix sort; the sorted values
//   don't depend on the algorithm); empty if a has NaN or negative values.
//
// B = omw_backprojector_mex('medfilt3', A, nThreads)
//   medfilt3(A, [3, 3, 3]) of psf_gen_new for a single 3D array A (symmetric padding, which for one voxel
//   repeats the border voxels): the 14th smallest of the 27 values around each voxel, with a selection
//   network (Batcher's odd-even merge sort of 32 values, 5 of them +Inf, pruned to that output); empty if
//   A has NaN or -0 values, whose order the network doesn't keep as medfilt3 does.
//
// Build: mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp -ffp-contract=off' LDFLAGS='$LDFLAGS -O3 -fopenmp' omw_backprojector_mex.cpp
// (no fused multiply-adds, so that the polygon lines are rasterized as in poly2mask; compileAllMexFiles
// builds it on all platforms)

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "mex.h"

typedef ptrdiff_t idx_t;

// ---------------------------------------------------------------------------------------------------------
// convex hull images of slices

// a hull point in doubled coordinates (2 * row, 2 * column), exact integers
struct Pt {
    int64_t x, y;
    bool operator<(const Pt& o) const { return x < o.x || (x == o.x && y < o.y); }
    bool operator==(const Pt& o) const { return x == o.x && y == o.y; }
};

static inline int64_t cross(const Pt& o, const Pt& a, const Pt& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// a 2D view of a slice: element (r, c) (0-based) at base[r * sr + c * sc]
struct Slice {
    uint8_t* base;
    idx_t h, w, sr, sc;
    inline uint8_t& at(const idx_t r, const idx_t c) const { return base[r * sr + c * sc]; }
};

// the counterclockwise hull of the sorted, distinct points p with all the points on its boundary, closed
// (the first point again at the end)
static void hull_with_boundary_points(const std::vector<Pt>& p, std::vector<Pt>& out) {
    const idx_t n = (idx_t)p.size();
    out.clear();
    if (n == 1) {
        out.push_back(p[0]);
        out.push_back(p[0]);
        return;
    }
    // strict hull (Andrew's monotone chain): lower chain, then upper chain
    std::vector<idx_t> hv(2 * n);
    idx_t k = 0;
    for (idx_t i = 0; i < n; i++) {
        while (k >= 2 && cross(p[hv[k - 2]], p[hv[k - 1]], p[i]) <= 0) {
            k--;
        }
        hv[k++] = i;
    }
    for (idx_t i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross(p[hv[k - 2]], p[hv[k - 1]], p[i]) <= 0) {
            k--;
        }
        hv[k++] = i;
    }
    // hv[0 .. k - 1] is closed (hv[k - 1] == hv[0])
    std::vector<std::pair<int64_t, idx_t> > on;
    for (idx_t e = 0; e + 1 < k; e++) {
        const Pt& a = p[hv[e]];
        const Pt& b = p[hv[e + 1]];
        out.push_back(a);
        // the other points on the edge, ordered from a to b
        const int64_t xl = std::min(a.x, b.x), xh = std::max(a.x, b.x);
        const int64_t yl = std::min(a.y, b.y), yh = std::max(a.y, b.y);
        const int64_t dx = b.x - a.x, dy = b.y - a.y;
        on.clear();
        // the points are sorted by x: scan the x range of the edge
        idx_t i0 = std::lower_bound(p.begin(), p.end(), Pt{xl, INT64_MIN}) - p.begin();
        for (idx_t i = i0; i < n && p[i].x <= xh; i++) {
            const Pt& q = p[i];
            if (q.y < yl || q.y > yh || q == a || q == b) {
                continue;
            }
            if (cross(a, b, q) == 0) {
                on.push_back(std::make_pair((q.x - a.x) * dx + (q.y - a.y) * dy, i));
            }
        }
        std::sort(on.begin(), on.end());
        for (size_t j = 0; j < on.size(); j++) {
            out.push_back(p[on[j].second]);
        }
    }
    out.push_back(out[0]);
}

// poly2mask(x, y, M, N) for a closed polygon (portable algorithm of toolbox/images/images/eml/poly2mask.m)
// into bw (M x N, column-major, zeroed by the caller)
static void poly2mask_portable(const std::vector<double>& xin, const std::vector<double>& yin, const int M,
        const int N, uint8_t* bw, std::vector<double>& xl, std::vector<double>& yl, std::vector<int32_t>& minY,
        std::vector<int32_t>& maxY) {
    const double scale = 5.0;
    const idx_t n = (idx_t)xin.size();
    std::vector<double> x(n), y(n);
    idx_t borderSize = 0;
    for (idx_t i = 0; i < n; i++) {
        x[i] = floor(scale * (xin[i] - 0.5) + 0.5) + 1.0;
        y[i] = floor(scale * (yin[i] - 0.5) + 0.5) + 1.0;
        if (i > 0) {
            const int32_t tdx = (int32_t)fabs(x[i] - x[i - 1]);
            const int32_t tdy = (int32_t)fabs(y[i] - y[i - 1]);
            borderSize += (tdx >= tdy ? tdx : tdy) + 1;
        }
    }
    xl.assign(borderSize, 0.0);
    yl.assign(borderSize, 0.0);
    idx_t pos = 0;
    for (idx_t i = 1; i < n; i++) {
        const double x1 = x[i - 1], y1 = y[i - 1], x2 = x[i], y2 = y[i];
        const int32_t dx = (int32_t)fabs(x2 - x1);
        const int32_t dy = (int32_t)fabs(y2 - y1);
        if (dx == 0 && dy == 0) {
            xl[pos] = x1;
            yl[pos] = y1;
            pos++;
            continue;
        }
        if (dx >= dy) {
            const double m = (y2 - y1) / (x2 - x1);
            if (x2 > x1) {
                for (double xv = x1; xv <= x2; xv += 1.0) {
                    const double t = m * (xv - x1);
                    xl[pos] = xv;
                    yl[pos] = round(y1 + t);
                    pos++;
                }
            } else {
                for (double xv = x1; xv >= x2; xv -= 1.0) {
                    const double t = m * (xv - x1);
                    xl[pos] = xv;
                    yl[pos] = round(y1 + t);
                    pos++;
                }
            }
        } else {
            const double m = (x2 - x1) / (y2 - y1);
            if (y2 > y1) {
                for (double yv = y1; yv <= y2; yv += 1.0) {
                    const double t = m * (yv - y1);
                    xl[pos] = round(x1 + t);
                    yl[pos] = yv;
                    pos++;
                }
            } else {
                for (double yv = y1; yv >= y2; yv -= 1.0) {
                    const double t = m * (yv - y1);
                    xl[pos] = round(x1 + t);
                    yl[pos] = yv;
                    pos++;
                }
            }
        }
    }
    minY.assign(N, 0);
    maxY.assign(N, 0);
    for (idx_t pt = 0; pt + 1 < borderSize; pt++) {
        const double diff = xl[pt + 1] - xl[pt];
        if (fabs(diff) >= 1) {
            double ypt = std::min(yl[pt], yl[pt + 1]);
            double xpt = xl[pt];
            if (diff < 0) {
                xpt = xpt - 1;
            }
            const double scaledDown = (xpt + (scale - 1) / 2) / scale;
            if (fabs(scaledDown - floor(scaledDown)) < 1 / (scale * 50)) {
                xpt = scaledDown;
                ypt = ceil((ypt + (scale - 1.0) / 2.0) / scale);
                const int32_t xUse = (int32_t)lround(xpt) + 1;
                int32_t yUse = (int32_t)lround(ypt) + 1;
                if (!(xUse < 2 || xUse > N + 1)) {
                    int32_t& mn = minY[xUse - 2];
                    int32_t& mx = maxY[xUse - 2];
                    if (yUse > M + 1) {
                        mn = (mn == 0) ? M + 1 : std::min(mn, M + 1);
                        mx = (mx == 0) ? M + 1 : std::max(mx, M + 1);
                    } else {
                        yUse = std::max(2, yUse);
                        uint8_t& o = bw[(idx_t)(xUse - 2) * M + (yUse - 2)];
                        o = !o;
                        mn = (mn == 0) ? yUse - 1 : std::min(mn, yUse - 1);
                        mx = (mx == 0) ? yUse : std::max(mx, yUse);
                    }
                }
            }
        }
    }
    // parity scan of each column
    for (int c = 0; c < N; c++) {
        bool pixel = false;
        uint8_t* col = bw + (idx_t)c * M;
        for (int32_t r = minY[c]; r <= maxY[c] - 1; r++) {
            if (col[r - 1]) {
                pixel = !pixel;
            }
            col[r - 1] = pixel;
        }
    }
}

// replace a slice with bwconvhull(slice, 'union')
static void convhull_slice(const Slice& s) {
    // first and last column of the pixels of each row (columns in memory order)
    std::vector<idx_t> cmin(s.h, -1), cmax(s.h, -1);
    for (idx_t c = 0; c < s.w; c++) {
        for (idx_t r = 0; r < s.h; r++) {
            if (s.at(r, c)) {
                if (cmin[r] < 0) {
                    cmin[r] = c;
                }
                cmax[r] = c;
            }
        }
    }
    idx_t r1 = -1, r2 = -1, c1 = s.w, c2 = -1;
    for (idx_t r = 0; r < s.h; r++) {
        if (cmin[r] >= 0) {
            if (r1 < 0) {
                r1 = r;
            }
            r2 = r;
            c1 = std::min(c1, cmin[r]);
            c2 = std::max(c2, cmax[r]);
        }
    }
    if (r1 < 0) {
        return;
    }
    // edge midpoints that can be on the hull boundary, in doubled 1-based coordinates
    std::vector<Pt> p;
    p.reserve(6 * (r2 - r1 + 1) + 2 * (c2 - c1 + 1));
    for (idx_t r = r1; r <= r2; r++) {
        if (cmin[r] < 0) {
            continue;
        }
        const int64_t R = 2 * (r + 1);
        const idx_t ce[2] = {cmin[r], cmax[r]};
        for (int k = 0; k < 2; k++) {
            const int64_t C = 2 * (ce[k] + 1);
            p.push_back(Pt{R - 1, C});
            p.push_back(Pt{R + 1, C});
        }
        p.push_back(Pt{R, 2 * (cmin[r] + 1) - 1});
        p.push_back(Pt{R, 2 * (cmax[r] + 1) + 1});
        if (r == r1 || r == r2) {
            for (idx_t c = cmin[r]; c <= cmax[r]; c++) {
                if (s.at(r, c)) {
                    p.push_back(Pt{r == r1 ? R - 1 : R + 1, 2 * (c + 1)});
                    if (r1 == r2) {
                        p.push_back(Pt{R + 1, 2 * (c + 1)});
                    }
                }
            }
        }
    }
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    std::vector<Pt> hull;
    hull_with_boundary_points(p, hull);

    // poly2mask of the hull within the bounding box: x = column, y = row
    const int M = (int)(r2 - r1 + 1), N = (int)(c2 - c1 + 1);
    std::vector<double> xs(hull.size()), ys(hull.size());
    for (size_t i = 0; i < hull.size(); i++) {
        xs[i] = (double)hull[i].y / 2 - (double)c1;
        ys[i] = (double)hull[i].x / 2 - (double)r1;
    }
    std::vector<uint8_t> bw((idx_t)M * N, 0);
    std::vector<double> xl, yl;
    std::vector<int32_t> minY, maxY;
    poly2mask_portable(xs, ys, M, N, bw.data(), xl, yl, minY, maxY);

    for (idx_t c = 0; c < s.w; c++) {
        for (idx_t r = 0; r < s.h; r++) {
            uint8_t v = 0;
            if (r >= r1 && r <= r2 && c >= c1 && c <= c2) {
                v = bw[(c - c1) * M + (r - r1)];
            }
            s.at(r, c) = v;
        }
    }
}

static void convhull_slices(uint8_t* B, const idx_t n1, const idx_t n2, const idx_t n3, const int dim,
        const std::vector<idx_t>& inds, const int nThreads) {
    const idx_t ns = (idx_t)inds.size();
    if (dim == 1) {
        // the slices A(y, :, :) are strided by n1: copied to and from contiguous buffers, plane by plane
        const idx_t sn = n2 * n3;
        std::vector<uint8_t> buf(ns * sn);
        #pragma omp parallel for num_threads(nThreads) schedule(static)
        for (idx_t k = 0; k < n3; k++) {
            for (idx_t i = 0; i < ns; i++) {
                const uint8_t* src = B + inds[i] + k * n1 * n2;
                uint8_t* dst = buf.data() + i * sn + k * n2;
                for (idx_t j = 0; j < n2; j++) {
                    dst[j] = src[j * n1];
                }
            }
        }
        #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
        for (idx_t i = 0; i < ns; i++) {
            convhull_slice(Slice{buf.data() + i * sn, n2, n3, 1, n2});
        }
        #pragma omp parallel for num_threads(nThreads) schedule(static)
        for (idx_t k = 0; k < n3; k++) {
            for (idx_t i = 0; i < ns; i++) {
                uint8_t* dst = B + inds[i] + k * n1 * n2;
                const uint8_t* src = buf.data() + i * sn + k * n2;
                for (idx_t j = 0; j < n2; j++) {
                    dst[j * n1] = src[j];
                }
            }
        }
        return;
    }
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (idx_t i = 0; i < ns; i++) {
        const idx_t k = inds[i];
        if (dim == 3) {
            convhull_slice(Slice{B + k * n1 * n2, n1, n2, 1, n1});
        } else {
            convhull_slice(Slice{B + k * n1, n1, n3, 1, n1 * n2});
        }
    }
}

static std::vector<idx_t> slice_indices(const mxArray* a, const idx_t n, const char* name) {
    if (!mxIsDouble(a) || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "%s must be real double indices.", name);
    }
    const double* v = mxGetPr(a);
    const idx_t m = (idx_t)mxGetNumberOfElements(a);
    std::vector<idx_t> out, u;
    out.resize(m);
    u.resize(m);
    for (idx_t i = 0; i < m; i++) {
        if (!(v[i] >= 1 && v[i] <= n && v[i] == floor(v[i]))) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "%s must be indices from 1 to %.0f.", name, (double)n);
        }
        out[i] = (idx_t)v[i] - 1;
        u[i] = out[i];
    }
    // each slice once, so that the parallel slices don't overlap
    std::sort(u.begin(), u.end());
    if (std::adjacent_find(u.begin(), u.end()) != u.end()) {
        mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "%s must not repeat indices.", name);
    }
    return out;
}

// ---------------------------------------------------------------------------------------------------------
// window bounds

static const int G = 128;            // cells per axis of [-1, 1]
static const double CELL = 2.0 / G;  // cell size
static const int MAXSHELL = 8;       // cells searched around a cell before all edge points are

static inline int cell_of(const double v) {
    int c = (int)floor((v + 1.0) * (G / 2.0));
    return c < 0 ? 0 : (c >= G ? G - 1 : c);
}

static inline idx_t cell_index(const int a, const int b, const int c) {
    return ((idx_t)c * G + b) * G + a;
}

// calls f(cell) for the cells at Chebyshev distance sh from cell (a, b, c)
template <class F>
static void visit_shell(const int a, const int b, const int c, const int sh, F f) {
    for (int kc = std::max(0, c - sh); kc <= std::min(G - 1, c + sh); kc++) {
        const bool zFace = abs(kc - c) == sh;
        for (int kb = std::max(0, b - sh); kb <= std::min(G - 1, b + sh); kb++) {
            if (zFace || abs(kb - b) == sh) {
                for (int ka = std::max(0, a - sh); ka <= std::min(G - 1, a + sh); ka++) {
                    f(cell_index(ka, kb, kc));
                }
            } else {
                if (a - sh >= 0) {
                    f(cell_index(a - sh, kb, kc));
                }
                if (a + sh <= G - 1) {
                    f(cell_index(a + sh, kb, kc));
                }
            }
        }
    }
}

static void window_bounds(const idx_t* psz, const double* E, const double* R, const idx_t nE, const double l,
        const double u, double* mask, std::vector<idx_t>& uinds, const int nThreads) {
    const idx_t nCells = (idx_t)G * G * G;
    // edge points by cell (CSR)
    std::vector<int32_t> start(nCells + 1, 0);
    std::vector<idx_t> ecell(nE);
    for (idx_t j = 0; j < nE; j++) {
        ecell[j] = cell_index(cell_of(E[j]), cell_of(E[j + nE]), cell_of(E[j + 2 * nE]));
        start[ecell[j] + 1]++;
    }
    for (idx_t i = 0; i < nCells; i++) {
        start[i + 1] += start[i];
    }
    std::vector<int32_t> fill(start.begin(), start.end() - 1);
    std::vector<int32_t> pts(nE);
    for (idx_t j = 0; j < nE; j++) {
        pts[fill[ecell[j]]++] = (int32_t)j;
    }

    // cells that meet the unit sphere (the directions of the voxels are unit vectors)
    std::vector<int32_t> slot(nCells, -1);
    std::vector<idx_t> sphereCells;
    for (int c = 0; c < G; c++) {
        for (int b = 0; b < G; b++) {
            for (int a = 0; a < G; a++) {
                const int ijk[3] = {a, b, c};
                double dmin = 0, dmax = 0;
                for (int t = 0; t < 3; t++) {
                    const double lo = -1.0 + ijk[t] * CELL, hi = lo + CELL;
                    const double nearv = (lo > 0) ? lo : ((hi < 0) ? hi : 0.0);
                    const double farv = std::max(fabs(lo), fabs(hi));
                    dmin += nearv * nearv;
                    dmax += farv * farv;
                }
                if (dmin <= (1 + 1e-6) * (1 + 1e-6) && dmax >= (1 - 1e-6) * (1 - 1e-6)) {
                    slot[cell_index(a, b, c)] = (int32_t)sphereCells.size();
                    sphereCells.push_back(cell_index(a, b, c));
                }
            }
        }
    }
    const idx_t nS = (idx_t)sphereCells.size();
    std::vector<double> Rlo(nS), Rhi(nS);
    // half diagonal of a cell, with a margin for the directions of the voxels as MATLAB computes them
    const double rho = sqrt(3.0) * (CELL / 2 + 1e-9);

    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 64)
    for (idx_t s = 0; s < nS; s++) {
        const idx_t ci = sphereCells[s];
        const int a = (int)(ci % G), b = (int)((ci / G) % G), c = (int)(ci / ((idx_t)G * G));
        const double ctr[3] = {-1.0 + (a + 0.5) * CELL, -1.0 + (b + 0.5) * CELL, -1.0 + (c + 0.5) * CELL};
        // distance to the 2nd nearest edge point, over shells of cells (points in shells beyond sh are
        // farther than (sh + 0.5) cells), or over all edge points if they are far
        double d1 = INFINITY, d2 = INFINITY;
        auto nearest2 = [&](const idx_t j) {
            const double e0 = E[j] - ctr[0], e1 = E[j + nE] - ctr[1], e2 = E[j + 2 * nE] - ctr[2];
            const double d = sqrt(e0 * e0 + e1 * e1 + e2 * e2);
            if (d < d1) {
                d2 = d1;
                d1 = d;
            } else if (d < d2) {
                d2 = d;
            }
        };
        bool found = false;
        for (int sh = 0; sh <= MAXSHELL && !found; sh++) {
            visit_shell(a, b, c, sh, [&](const idx_t k) {
                for (int32_t q = start[k]; q < start[k + 1]; q++) {
                    nearest2(pts[q]);
                }
            });
            found = d2 <= (sh + 0.5) * CELL;
        }
        if (!found) {
            d1 = d2 = INFINITY;
            for (idx_t j = 0; j < nE; j++) {
                nearest2(j);
            }
        }
        // the 2 nearest edge points of any direction in the cell are within rad of the cell center
        const double rad = d2 + 2 * rho + 1e-12;
        const int S = (int)ceil(rad / CELL + 0.5);
        double lo = INFINITY, hi = -INFINITY;
        auto bound = [&](const idx_t j) {
            const double e0 = E[j] - ctr[0], e1 = E[j + nE] - ctr[1], e2 = E[j + 2 * nE] - ctr[2];
            if (sqrt(e0 * e0 + e1 * e1 + e2 * e2) <= rad) {
                lo = std::min(lo, R[j]);
                hi = std::max(hi, R[j]);
            }
        };
        if (S <= MAXSHELL) {
            for (int sh = 0; sh <= S; sh++) {
                visit_shell(a, b, c, sh, [&](const idx_t k) {
                    for (int32_t q = start[k]; q < start[k + 1]; q++) {
                        bound(pts[q]);
                    }
                });
            }
        } else {
            for (idx_t j = 0; j < nE; j++) {
                bound(j);
            }
        }
        Rlo[s] = lo;
        Rhi[s] = hi;
    }

    // voxels: 1 where the window is 1, 0 where it is 0, 2 where it has to be computed; first with the
    // bounds of all edge radii (on the squared radius, with a margin for the square root)
    const idx_t n1 = psz[0], n2 = psz[1], n3 = psz[2];
    const double cy = (n1 + 1) / 2.0, cx = (n2 + 1) / 2.0, cz = (n3 + 1) / 2.0;
    double Rmin = INFINITY, Rmax = -INFINITY;
    for (idx_t j = 0; j < nE; j++) {
        Rmin = std::min(Rmin, R[j]);
        Rmax = std::max(Rmax, R[j]);
    }
    const double in2 = (l * Rmin * (1 - 1e-9)) * (l * Rmin * (1 - 1e-9)) * (1 - 1e-12);
    const double out2 = (u * Rmax * (1 + 1e-9)) * (u * Rmax * (1 + 1e-9)) * (1 + 1e-12);
    std::vector<uint8_t> cls(n1 * n2 * n3);
    std::vector<idx_t> counts(n3, 0);
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (idx_t k = 0; k < n3; k++) {
        const double Z = (k + 1) - cz;
        idx_t cnt = 0;
        for (idx_t j = 0; j < n2; j++) {
            const double X = (j + 1) - cx;
            const idx_t off = (k * n2 + j) * n1;
            for (idx_t i = 0; i < n1; i++) {
                const double Y = (i + 1) - cy;
                const double r2 = Y * Y + X * X + Z * Z;
                uint8_t v = 2;
                if (r2 > 0 && r2 < in2) {
                    v = 1;
                } else if (r2 >= out2) {
                    v = 0;
                } else if (r2 > 0) {
                    const double r = sqrt(r2);
                    const double inv = 1.0 / r;
                    const int32_t sl = slot[cell_index(cell_of(Y * inv), cell_of(X * inv), cell_of(Z * inv))];
                    if (sl >= 0 && Rlo[sl] <= Rhi[sl]) {
                        if (r < l * Rlo[sl] * (1 - 1e-9)) {
                            v = 1;
                        } else if (r >= u * Rhi[sl] * (1 + 1e-9)) {
                            v = 0;
                        }
                    }
                }
                cls[off + i] = v;
                mask[off + i] = (v == 1) ? 1.0 : 0.0;
                cnt += (v == 2);
            }
        }
        counts[k] = cnt;
    }
    std::vector<idx_t> offs(n3 + 1, 0);
    for (idx_t k = 0; k < n3; k++) {
        offs[k + 1] = offs[k] + counts[k];
    }
    uinds.resize(offs[n3]);
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (idx_t k = 0; k < n3; k++) {
        idx_t p = offs[k];
        const idx_t base = k * n1 * n2;
        for (idx_t i = 0; i < n1 * n2; i++) {
            if (cls[base + i] == 2) {
                uinds[p++] = base + i;
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------------------
// connected components (26-connected)

static inline int32_t find_root(int32_t* parent, int32_t i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

// union by the smaller index, so that the root of a component is its first voxel in linear order (and a
// parent is always before its child)
static inline void unite(int32_t* parent, const int32_t a, const int32_t b) {
    const int32_t ra = find_root(parent, a), rb = find_root(parent, b);
    if (ra < rb) {
        parent[rb] = ra;
    } else if (rb < ra) {
        parent[ra] = rb;
    }
}

// labels of the components of A in lab (1-based, in the order of their first voxels; 0: background) and
// their sizes (sizes[label])
static void label_components(const uint8_t* A, const idx_t n1, const idx_t n2, const idx_t n3, int32_t* lab,
        std::vector<idx_t>& sizes, const int nThreads) {
    const idx_t np = n1 * n2;
    // offsets (x, y, z; x along dimension 1) of the 13 neighbors before a voxel in linear order
    int nb[13][3];
    int m = 0;
    for (int dz = -1; dz <= 0; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dz < 0 || (dz == 0 && dy < 0) || (dz == 0 && dy == 0 && dx < 0)) {
                    nb[m][0] = dx;
                    nb[m][1] = dy;
                    nb[m][2] = dz;
                    m++;
                }
            }
        }
    }
    // components within slabs of planes in parallel, then joined across the slabs
    const int nt = (int)std::max<idx_t>(1, std::min<idx_t>(nThreads, n3));
    std::vector<idx_t> zs(nt + 1);
    for (int t = 0; t <= nt; t++) {
        zs[t] = n3 * t / nt;
    }
    int32_t* parent = lab;
    auto link = [&](const idx_t z, const idx_t y, const idx_t x, const idx_t zmin) {
        const int32_t i = (int32_t)(z * np + y * n1 + x);
        for (int k = 0; k < 13; k++) {
            const idx_t xx = x + nb[k][0], yy = y + nb[k][1], zz = z + nb[k][2];
            if (xx < 0 || xx >= n1 || yy < 0 || yy >= n2 || zz < zmin) {
                continue;
            }
            const int32_t j = (int32_t)(zz * np + yy * n1 + xx);
            if (A[j]) {
                unite(parent, i, j);
            }
        }
    };
    #pragma omp parallel for num_threads(nt) schedule(static, 1)
    for (int t = 0; t < nt; t++) {
        for (idx_t z = zs[t]; z < zs[t + 1]; z++) {
            for (idx_t y = 0; y < n2; y++) {
                for (idx_t x = 0; x < n1; x++) {
                    const idx_t i = z * np + y * n1 + x;
                    parent[i] = A[i] ? (int32_t)i : -1;
                    if (A[i]) {
                        link(z, y, x, zs[t]);
                    }
                }
            }
        }
    }
    for (int t = 1; t < nt; t++) {
        const idx_t z = zs[t];
        for (idx_t y = 0; y < n2; y++) {
            for (idx_t x = 0; x < n1; x++) {
                if (A[z * np + y * n1 + x]) {
                    link(z, y, x, z - 1);
                }
            }
        }
    }
    // in linear order a voxel's parent is final (its root) before the voxel, so one step finds the root;
    // a root gets the next label, and the other voxels the label already in their root's entry
    sizes.assign(1, 0);
    const idx_t n = np * n3;
    int32_t nl = 0;
    for (idx_t i = 0; i < n; i++) {
        const int32_t p = parent[i];
        if (p < 0) {
            lab[i] = 0;
            continue;
        }
        int32_t l;
        if (p == (int32_t)i) {
            l = ++nl;
            sizes.push_back(0);
        } else {
            // the parent was visited: its entry holds -(label) of the root
            l = -lab[p];
        }
        sizes[l]++;
        lab[i] = -l;
    }
    #pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        lab[i] = -lab[i];
    }
}

static void keep_components(const uint8_t* A, const idx_t n1, const idx_t n2, const idx_t n3, const bool largest,
        const idx_t n, uint8_t* B, const int nThreads) {
    const idx_t nv = n1 * n2 * n3;
    std::vector<int32_t> lab;
    lab.resize(nv);
    std::vector<idx_t> sizes;
    label_components(A, n1, n2, n3, lab.data(), sizes, nThreads);
    const idx_t nl = (idx_t)sizes.size() - 1;
    std::vector<uint8_t> keep;
    keep.assign(nl + 1, 0);
    if (largest) {
        // the n largest, ties to the earlier component
        std::vector<idx_t> order;
        order.resize(nl);
        for (idx_t k = 0; k < nl; k++) {
            order[k] = k + 1;
        }
        std::stable_sort(order.begin(), order.end(), [&](const idx_t a, const idx_t b) { return sizes[a] > sizes[b]; });
        for (idx_t k = 0; k < std::min(n, nl); k++) {
            keep[order[k]] = 1;
        }
    } else {
        for (idx_t k = 1; k <= nl; k++) {
            keep[k] = sizes[k] >= n;
        }
    }
    #pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t i = 0; i < nv; i++) {
        B[i] = keep[lab[i]];
    }
}

// ---------------------------------------------------------------------------------------------------------
// descending sort of nonnegative singles

static bool sort_desc(const float* a, const idx_t n, double* out, const int nThreads) {
    // nonnegative floats sort as their bits
    std::vector<uint32_t> k0, k1;
    k0.resize(n);
    k1.resize(n);
    int bad = 0;
    #pragma omp parallel for num_threads(nThreads) schedule(static) reduction(| : bad)
    for (idx_t i = 0; i < n; i++) {
        uint32_t b;
        memcpy(&b, a + i, 4);
        // NaN and negative values (-0 included) are left to MATLAB's sort
        bad |= (b > 0x7F800000u);
        k0[i] = b;
    }
    if (bad) {
        return false;
    }
    // LSD radix sort, 4 passes of 8 bits, ascending
    const int nt = (int)std::max<idx_t>(1, std::min<idx_t>(nThreads, n / 65536 + 1));
    std::vector<idx_t> hist;
    hist.resize((size_t)nt * 256);
    uint32_t* src = k0.data();
    uint32_t* dst = k1.data();
    for (int pass = 0; pass < 4; pass++) {
        const int sh = 8 * pass;
        std::fill(hist.begin(), hist.end(), 0);
        #pragma omp parallel num_threads(nt)
        {
            const int t = omp_get_thread_num();
            const idx_t i0 = n * t / nt, i1 = n * (t + 1) / nt;
            idx_t* h = hist.data() + (size_t)t * 256;
            for (idx_t i = i0; i < i1; i++) {
                h[(src[i] >> sh) & 255]++;
            }
            #pragma omp barrier
            #pragma omp single
            {
                idx_t sum = 0;
                for (int b = 0; b < 256; b++) {
                    for (int tt = 0; tt < nt; tt++) {
                        const idx_t c = hist[(size_t)tt * 256 + b];
                        hist[(size_t)tt * 256 + b] = sum;
                        sum += c;
                    }
                }
            }
            for (idx_t i = i0; i < i1; i++) {
                dst[h[(src[i] >> sh) & 255]++] = src[i];
            }
        }
        std::swap(src, dst);
    }
    #pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        float f;
        memcpy(&f, src + (n - 1 - i), 4);
        out[i] = (double)f;
    }
    return true;
}

// ---------------------------------------------------------------------------------------------------------
// 3x3x3 median filter

static const int ML = 16;  // voxels along dimension 1 per network pass
// ML floats (GCC vector extension: SSE or NEON operations, without branches)
typedef float vfloat __attribute__((vector_size(4 * ML)));

// the 14th smallest of v[0 .. 26] in v[13], for ML voxels at once
static inline void median27(vfloat* v) {
#define CX(a, b) { const vfloat x = v[a], y = v[b]; v[a] = y < x ? y : x; v[b] = y < x ? x : y; }
    CX(0, 1); CX(2, 3); CX(0, 2); CX(1, 3); CX(1, 2); CX(4, 5);
    CX(6, 7); CX(4, 6); CX(5, 7); CX(5, 6); CX(0, 4); CX(2, 6);
    CX(2, 4); CX(1, 5); CX(3, 7); CX(3, 5); CX(1, 2); CX(3, 4);
    CX(5, 6); CX(8, 9); CX(10, 11); CX(8, 10); CX(9, 11); CX(9, 10);
    CX(12, 13); CX(14, 15); CX(12, 14); CX(13, 15); CX(13, 14); CX(8, 12);
    CX(10, 14); CX(10, 12); CX(9, 13); CX(11, 15); CX(11, 13); CX(9, 10);
    CX(11, 12); CX(13, 14); CX(0, 8); CX(4, 12); CX(4, 8); CX(2, 10);
    CX(6, 14); CX(6, 10); CX(2, 4); CX(6, 8); CX(10, 12); CX(1, 9);
    CX(5, 13); CX(5, 9); CX(3, 11); CX(7, 15); CX(7, 11); CX(3, 5);
    CX(7, 9); CX(11, 13); CX(1, 2); CX(3, 4); CX(5, 6); CX(7, 8);
    CX(9, 10); CX(11, 12); CX(13, 14); CX(16, 17); CX(18, 19); CX(16, 18);
    CX(17, 19); CX(17, 18); CX(20, 21); CX(22, 23); CX(20, 22); CX(21, 23);
    CX(21, 22); CX(16, 20); CX(18, 22); CX(18, 20); CX(17, 21); CX(19, 23);
    CX(19, 21); CX(17, 18); CX(19, 20); CX(21, 22); CX(24, 25); CX(24, 26);
    CX(25, 26); CX(25, 26); CX(16, 24); CX(20, 24); CX(18, 26); CX(22, 26);
    CX(18, 20); CX(22, 24); CX(17, 25); CX(21, 25); CX(19, 21); CX(23, 25);
    CX(17, 18); CX(19, 20); CX(21, 22); CX(23, 24); CX(25, 26); CX(0, 16);
    CX(8, 24); CX(8, 16); CX(4, 20); CX(12, 20); CX(12, 16); CX(2, 18);
    CX(10, 26); CX(10, 18); CX(6, 22); CX(14, 22); CX(14, 18); CX(14, 16);
    CX(1, 17); CX(9, 25); CX(9, 17); CX(5, 21); CX(13, 21); CX(13, 17);
    CX(3, 19); CX(11, 19); CX(7, 23); CX(7, 11); CX(11, 13); CX(13, 14);
#undef CX
}

static bool medfilt3(const float* A, const idx_t n1, const idx_t n2, const idx_t n3, float* B, const int nThreads) {
    const idx_t n = n1 * n2 * n3;
    int bad = 0;
    #pragma omp parallel for num_threads(nThreads) schedule(static) reduction(| : bad)
    for (idx_t i = 0; i < n; i++) {
        bad |= (A[i] != A[i]) || (A[i] == 0 && signbit(A[i]));
    }
    if (bad) {
        return false;
    }
    // rows along dimension 1 with the border voxel repeated at both ends (and room for the last pass)
    const idx_t w = n1 + 2;
    std::vector<float> P;
    P.resize(w * n2 * n3 + ML);
    #pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t r = 0; r < n2 * n3; r++) {
        float* d = P.data() + r * w;
        const float* a = A + r * n1;
        d[0] = a[0];
        memcpy(d + 1, a, n1 * sizeof(float));
        d[n1 + 1] = a[n1 - 1];
    }
    #pragma omp parallel for num_threads(nThreads) schedule(dynamic, 4)
    for (idx_t r = 0; r < n2 * n3; r++) {
        const idx_t j = r % n2, k = r / n2;
        const float* rows[9];
        int m = 0;
        for (int dz = -1; dz <= 1; dz++) {
            const idx_t kk = std::min(std::max<idx_t>(k + dz, 0), n3 - 1);
            for (int dy = -1; dy <= 1; dy++) {
                const idx_t jj = std::min(std::max<idx_t>(j + dy, 0), n2 - 1);
                rows[m++] = P.data() + (kk * n2 + jj) * w;
            }
        }
        vfloat v[27];
        float* out = B + r * n1;
        for (idx_t i0 = 0; i0 < n1; i0 += ML) {
            for (int q = 0; q < 9; q++) {
                for (int dx = 0; dx < 3; dx++) {
                    memcpy(&v[q * 3 + dx], rows[q] + i0 + dx, sizeof(vfloat));
                }
            }
            median27(v);
            const idx_t nl = std::min<idx_t>(ML, n1 - i0);
            float med[ML];
            memcpy(med, &v[13], sizeof(vfloat));
            memcpy(out + i0, med, nl * sizeof(float));
        }
    }
    return true;
}

// ---------------------------------------------------------------------------------------------------------

static int get_threads(const mxArray* a) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1) {
        mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "nThreads must be a number.");
    }
    const double t = mxGetScalar(a);
    return t >= 1 ? (int)t : 1;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 1 || !mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "The first input must be 'convhull', 'window', 'components', 'sortdesc' or 'medfilt3'.");
    }
    char mode[16];
    mxGetString(prhs[0], mode, sizeof(mode));

    if (strcmp(mode, "convhull") == 0) {
        if (nrhs != 6) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "Usage: B = omw_backprojector_mex('convhull', A, zinds, yinds, xinds, nThreads)");
        }
        const mxArray* A = prhs[1];
        const mwSize* d = mxGetDimensions(A);
        if (!mxIsLogical(A) || mxGetNumberOfDimensions(A) != 3 || d[0] < 2 || d[1] < 2 || d[2] < 2) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "A must be a 3D logical array with all dimensions > 1.");
        }
        const idx_t n1 = d[0], n2 = d[1], n3 = d[2];
        std::vector<idx_t> zi = slice_indices(prhs[2], n3, "zinds");
        std::vector<idx_t> yi = slice_indices(prhs[3], n1, "yinds");
        std::vector<idx_t> xi = slice_indices(prhs[4], n2, "xinds");
        const int nThreads = get_threads(prhs[5]);
        plhs[0] = mxDuplicateArray(A);
        uint8_t* B = (uint8_t*)mxGetData(plhs[0]);
        convhull_slices(B, n1, n2, n3, 3, zi, nThreads);
        convhull_slices(B, n1, n2, n3, 1, yi, nThreads);
        convhull_slices(B, n1, n2, n3, 2, xi, nThreads);
        return;
    }

    if (strcmp(mode, "window") == 0) {
        if (nrhs != 6) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "Usage: [mask, inds] = omw_backprojector_mex('window', psz, E, R, bounds, nThreads)");
        }
        if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 3) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "psz must be 3 numbers.");
        }
        const double* pv = mxGetPr(prhs[1]);
        idx_t psz[3];
        for (int t = 0; t < 3; t++) {
            if (!(pv[t] >= 1 && pv[t] == floor(pv[t]))) {
                mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "psz must be positive integers.");
            }
            psz[t] = (idx_t)pv[t];
        }
        const mxArray* E = prhs[2];
        const mxArray* R = prhs[3];
        const idx_t nE = (idx_t)mxGetM(E);
        if (!mxIsDouble(E) || mxIsComplex(E) || mxGetN(E) != 3 || nE < 2 || !mxIsDouble(R) || mxIsComplex(R)
                || (idx_t)mxGetNumberOfElements(R) != nE) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "E must be an n x 3 double array (n >= 2) and R n doubles.");
        }
        const double* Ev = mxGetPr(E);
        for (idx_t j = 0; j < 3 * nE; j++) {
            if (!(fabs(Ev[j]) <= 1.5)) {
                mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "E must be unit vectors.");
            }
        }
        if (!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 2) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "bounds must be [l, u].");
        }
        const double l = mxGetPr(prhs[4])[0], u = mxGetPr(prhs[4])[1];
        if (!(l > 0 && l < u && u < INFINITY)) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "bounds must be [l, u] with 0 < l < u < Inf.");
        }
        const int nThreads = get_threads(prhs[5]);
        const mwSize dims[3] = {(mwSize)psz[0], (mwSize)psz[1], (mwSize)psz[2]};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        std::vector<idx_t> uinds;
        window_bounds(psz, Ev, mxGetPr(R), nE, l, u, mxGetPr(plhs[0]), uinds, nThreads);
        plhs[1] = mxCreateDoubleMatrix((mwSize)uinds.size(), 1, mxREAL);
        double* o = mxGetPr(plhs[1]);
        for (size_t i = 0; i < uinds.size(); i++) {
            o[i] = (double)(uinds[i] + 1);
        }
        return;
    }

    if (strcmp(mode, "components") == 0) {
        if (nrhs != 5 || !mxIsChar(prhs[2])) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "Usage: B = omw_backprojector_mex('components', A, 'minsize'|'largest', n, nThreads)");
        }
        const mxArray* A = prhs[1];
        const mwSize* d = mxGetDimensions(A);
        if (!mxIsLogical(A) || mxGetNumberOfDimensions(A) != 3) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "A must be a 3D logical array.");
        }
        if ((double)d[0] * d[1] * d[2] >= 2147483647.0) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "A must have fewer than 2^31 - 1 elements.");
        }
        char cmode[16];
        mxGetString(prhs[2], cmode, sizeof(cmode));
        const bool largest = strcmp(cmode, "largest") == 0;
        if (!largest && strcmp(cmode, "minsize") != 0) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "The mode must be 'minsize' or 'largest'.");
        }
        if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1 || !(mxGetScalar(prhs[3]) >= 0)
                || mxGetScalar(prhs[3]) != floor(mxGetScalar(prhs[3])) || mxGetScalar(prhs[3]) > 9e15) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "n must be a nonnegative integer.");
        }
        const idx_t n = (idx_t)mxGetScalar(prhs[3]);
        const int nThreads = get_threads(prhs[4]);
        plhs[0] = mxCreateLogicalArray(3, d);
        keep_components((const uint8_t*)mxGetData(A), d[0], d[1], d[2], largest, n, (uint8_t*)mxGetData(plhs[0]), nThreads);
        return;
    }

    if (strcmp(mode, "sortdesc") == 0) {
        if (nrhs != 3 || !mxIsSingle(prhs[1]) || mxIsComplex(prhs[1])) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "Usage: v = omw_backprojector_mex('sortdesc', a, nThreads) with real single a");
        }
        const idx_t n = (idx_t)mxGetNumberOfElements(prhs[1]);
        const int nThreads = get_threads(prhs[2]);
        plhs[0] = mxCreateDoubleMatrix((mwSize)n, 1, mxREAL);
        if (!sort_desc((const float*)mxGetData(prhs[1]), n, mxGetPr(plhs[0]), nThreads)) {
            mxDestroyArray(plhs[0]);
            plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
        return;
    }

    if (strcmp(mode, "medfilt3") == 0) {
        if (nrhs != 3 || !mxIsSingle(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) > 3) {
            mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "Usage: B = omw_backprojector_mex('medfilt3', A, nThreads) with real single 3D A");
        }
        const mwSize* d = mxGetDimensions(prhs[1]);
        const idx_t n1 = d[0], n2 = d[1], n3 = mxGetNumberOfDimensions(prhs[1]) == 3 ? d[2] : 1;
        const int nThreads = get_threads(prhs[2]);
        if (n1 * n2 * n3 == 0) {
            plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
            return;
        }
        plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]), d, mxSINGLE_CLASS, mxREAL);
        if (!medfilt3((const float*)mxGetData(prhs[1]), n1, n2, n3, (float*)mxGetData(plhs[0]), nThreads)) {
            mxDestroyArray(plhs[0]);
            plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
        return;
    }

    mexErrMsgIdAndTxt("omw_backprojector_mex:invalidInput", "Unknown mode '%s'.", mode);
}
