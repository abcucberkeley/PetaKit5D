// bwconvhull(BW, 'union') of 2D images, bit for bit: regionprops' ConvexImage of all true pixels,
// placed in their bounding box. The hull is that of the midpoints of the edges of the pixels
// (regionprops' PerimeterCornerPixelList), with all the points on its boundary (convhull keeps
// collinear points), counterclockwise in (row, column) as convhull returns it; it is filled
// with poly2mask's algorithm (vertices on a 5 x 5 subpixel grid, edge list, parity scan; the
// portable version in toolbox/images/images/eml/poly2mask.m). Only the edge midpoints that can
// be on the hull boundary are used (of the first and last pixels of each row, and all of the
// first and last row), which gives the same boundary points. The algorithm and its checks
// against bwconvhull come from omw_backprojector_mex.cpp (its 'convhull' mode).
#include <pk5d/common.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

namespace nb = nanobind;
using idx_t = std::ptrdiff_t;

// a hull point in doubled coordinates (2 * row, 2 * column), exact integers
struct Pt {
    int64_t x, y;
    bool operator<(const Pt& o) const { return x < o.x || (x == o.x && y < o.y); }
    bool operator==(const Pt& o) const { return x == o.x && y == o.y; }
};

inline int64_t cross(const Pt& o, const Pt& a, const Pt& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// a 2D view of a slice: element (r, c) (0-based) at base[r * sr + c * sc]
struct Slice {
    uint8_t* base;
    idx_t h, w, sr, sc;
    inline uint8_t& at(const idx_t r, const idx_t c) const { return base[r * sr + c * sc]; }
};

// the counterclockwise hull of the sorted, distinct points p with all the points on its
// boundary, closed (the first point again at the end)
void hull_with_boundary_points(const std::vector<Pt>& p, std::vector<Pt>& out) {
    const idx_t n = (idx_t)p.size();
    out.clear();
    if (n == 1) {
        out.push_back(p[0]);
        out.push_back(p[0]);
        return;
    }
    // strict hull (Andrew's monotone chain): lower chain, then upper chain
    std::vector<idx_t> hv(2 * (size_t)n);
    idx_t k = 0;
    for (idx_t i = 0; i < n; i++) {
        while (k >= 2 && cross(p[hv[k - 2]], p[hv[k - 1]], p[i]) <= 0) k--;
        hv[k++] = i;
    }
    for (idx_t i = n - 2, t = k + 1; i >= 0; i--) {
        while (k >= t && cross(p[hv[k - 2]], p[hv[k - 1]], p[i]) <= 0) k--;
        hv[k++] = i;
    }
    // hv[0 .. k - 1] is closed (hv[k - 1] == hv[0])
    std::vector<std::pair<int64_t, idx_t>> on;
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
            if (q.y < yl || q.y > yh || q == a || q == b) continue;
            if (cross(a, b, q) == 0) on.push_back(std::make_pair((q.x - a.x) * dx + (q.y - a.y) * dy, i));
        }
        std::sort(on.begin(), on.end());
        for (size_t j = 0; j < on.size(); j++) out.push_back(p[on[j].second]);
    }
    out.push_back(out[0]);
}

// poly2mask(x, y, M, N) for a closed polygon (portable algorithm of
// toolbox/images/images/eml/poly2mask.m) into bw (M x N, column-major, zeroed by the caller)
void poly2mask_portable(const std::vector<double>& xin, const std::vector<double>& yin, const int M, const int N,
                        uint8_t* bw) {
    const double scale = 5.0;
    const idx_t n = (idx_t)xin.size();
    std::vector<double> x(n), y(n);
    idx_t borderSize = 0;
    for (idx_t i = 0; i < n; i++) {
        x[i] = std::floor(scale * (xin[i] - 0.5) + 0.5) + 1.0;
        y[i] = std::floor(scale * (yin[i] - 0.5) + 0.5) + 1.0;
        if (i > 0) {
            const int32_t tdx = (int32_t)std::fabs(x[i] - x[i - 1]);
            const int32_t tdy = (int32_t)std::fabs(y[i] - y[i - 1]);
            borderSize += (tdx >= tdy ? tdx : tdy) + 1;
        }
    }
    std::vector<double> xl(borderSize, 0.0), yl(borderSize, 0.0);
    idx_t pos = 0;
    for (idx_t i = 1; i < n; i++) {
        const double x1 = x[i - 1], y1 = y[i - 1], x2 = x[i], y2 = y[i];
        const int32_t dx = (int32_t)std::fabs(x2 - x1);
        const int32_t dy = (int32_t)std::fabs(y2 - y1);
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
                    yl[pos] = std::round(y1 + t);
                    pos++;
                }
            } else {
                for (double xv = x1; xv >= x2; xv -= 1.0) {
                    const double t = m * (xv - x1);
                    xl[pos] = xv;
                    yl[pos] = std::round(y1 + t);
                    pos++;
                }
            }
        } else {
            const double m = (x2 - x1) / (y2 - y1);
            if (y2 > y1) {
                for (double yv = y1; yv <= y2; yv += 1.0) {
                    const double t = m * (yv - y1);
                    xl[pos] = std::round(x1 + t);
                    yl[pos] = yv;
                    pos++;
                }
            } else {
                for (double yv = y1; yv >= y2; yv -= 1.0) {
                    const double t = m * (yv - y1);
                    xl[pos] = std::round(x1 + t);
                    yl[pos] = yv;
                    pos++;
                }
            }
        }
    }
    std::vector<int32_t> minY(N, 0), maxY(N, 0);
    for (idx_t pt = 0; pt + 1 < borderSize; pt++) {
        const double diff = xl[pt + 1] - xl[pt];
        if (std::fabs(diff) >= 1) {
            double ypt = std::min(yl[pt], yl[pt + 1]);
            double xpt = xl[pt];
            if (diff < 0) xpt = xpt - 1;
            const double scaledDown = (xpt + (scale - 1) / 2) / scale;
            if (std::fabs(scaledDown - std::floor(scaledDown)) < 1 / (scale * 50)) {
                xpt = scaledDown;
                ypt = std::ceil((ypt + (scale - 1.0) / 2.0) / scale);
                const int32_t xUse = (int32_t)std::lround(xpt) + 1;
                int32_t yUse = (int32_t)std::lround(ypt) + 1;
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
            if (col[r - 1]) pixel = !pixel;
            col[r - 1] = pixel;
        }
    }
}

// replace a slice with bwconvhull(slice, 'union')
void convhull_slice(const Slice& s) {
    // first and last column of the pixels of each row (columns in memory order)
    std::vector<idx_t> cmin(s.h, -1), cmax(s.h, -1);
    for (idx_t c = 0; c < s.w; c++) {
        for (idx_t r = 0; r < s.h; r++) {
            if (s.at(r, c)) {
                if (cmin[r] < 0) cmin[r] = c;
                cmax[r] = c;
            }
        }
    }
    idx_t r1 = -1, r2 = -1, c1 = s.w, c2 = -1;
    for (idx_t r = 0; r < s.h; r++) {
        if (cmin[r] >= 0) {
            if (r1 < 0) r1 = r;
            r2 = r;
            c1 = std::min(c1, cmin[r]);
            c2 = std::max(c2, cmax[r]);
        }
    }
    if (r1 < 0) return;
    // edge midpoints that can be on the hull boundary, in doubled 1-based coordinates
    std::vector<Pt> p;
    p.reserve(6 * (size_t)(r2 - r1 + 1) + 2 * (size_t)(c2 - c1 + 1));
    for (idx_t r = r1; r <= r2; r++) {
        if (cmin[r] < 0) continue;
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
                    if (r1 == r2) p.push_back(Pt{R + 1, 2 * (c + 1)});
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
    std::vector<uint8_t> bw((size_t)M * (size_t)N, 0);
    poly2mask_portable(xs, ys, M, N, bw.data());
    for (idx_t c = 0; c < s.w; c++) {
        for (idx_t r = 0; r < s.h; r++) {
            uint8_t v = 0;
            if (r >= r1 && r <= r2 && c >= c1 && c <= c2) v = bw[(size_t)((c - c1) * M + (r - r1))];
            s.at(r, c) = v;
        }
    }
}

// The slices of a C-contiguous (nz, ny, nx) volume B along a NumPy axis at the indices inds,
// replaced by their convex hull images, the slices in parallel. A slice keeps MATLAB's rows and
// columns (poly2mask is not symmetric in them): axis 0 (MATLAB's B(:, :, k)): rows y, columns
// x; axis 1 (squeeze(B(k, :, :))): rows x, columns z; axis 2 (squeeze(B(:, k, :))): rows y,
// columns z. Each slice is copied to a contiguous buffer of its thread and back.
void convhull_slices(uint8_t* B, const idx_t nz, const idx_t ny, const idx_t nx, const int axis,
                     const std::vector<idx_t>& inds, const int nt) {
    const idx_t ns = (idx_t)inds.size();
    const idx_t plane = ny * nx;
    // rows, columns and their strides in B
    const idx_t h = axis == 1 ? nx : ny, w = axis == 0 ? nx : nz;
    const idx_t sr = axis == 1 ? 1 : nx, sc = axis == 0 ? 1 : plane;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
    {
        std::vector<uint8_t> buf((size_t)(h * w));
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (idx_t i = 0; i < ns; i++) {
            const idx_t k = inds[i];
            uint8_t* base = B + (axis == 0 ? k * plane : axis == 1 ? k * nx : k);
            for (idx_t c = 0; c < w; c++) {
                for (idx_t r = 0; r < h; r++) buf[(size_t)(c * h + r)] = base[r * sr + c * sc];
            }
            convhull_slice(Slice{buf.data(), h, w, 1, h});
            for (idx_t c = 0; c < w; c++) {
                for (idx_t r = 0; r < h; r++) base[r * sr + c * sc] = buf[(size_t)(c * h + r)];
            }
        }
    }
}

using bvol = pk5d::cvolume<bool>;
using ivec = nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

void bwconvhull_slices(bvol b, int axis, ivec inds, int num_threads) {
    pk5d::check_threads(num_threads);
    if (axis < 0 || axis > 2) throw std::invalid_argument("bwconvhull_slices: axis must be 0, 1 or 2");
    const idx_t n[3] = {(idx_t)b.shape(0), (idx_t)b.shape(1), (idx_t)b.shape(2)};
    std::vector<idx_t> v((size_t)inds.shape(0)), u;
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = (idx_t)inds((idx_t)i);
        if (v[i] < 0 || v[i] >= n[axis]) throw std::invalid_argument("bwconvhull_slices: index out of range");
    }
    u = v;
    std::sort(u.begin(), u.end());
    if (std::adjacent_find(u.begin(), u.end()) != u.end())
        throw std::invalid_argument("bwconvhull_slices: repeated slice index");
    if (n[0] * n[1] * n[2] == 0) return;
    uint8_t* B = reinterpret_cast<uint8_t*>(b.data());
    nb::gil_scoped_release release;
    convhull_slices(B, n[0], n[1], n[2], axis, v, num_threads);
}

}  // namespace

PK5D_BINDINGS(bind_hull) {
    m.def("bwconvhull_slices", &bwconvhull_slices, nb::arg("b").noconvert(), nb::arg("axis"),
          nb::arg("inds").noconvert(), nb::arg("num_threads"),
          "Replace the slices of a C-contiguous (z, y, x) bool volume along a NumPy axis at the "
          "0-based indices inds by bwconvhull(slice, 'union') of MATLAB's slice (axis 0: rows y, "
          "columns x; axis 1: rows x, columns z; axis 2: rows y, columns z), in place.");
}
