// petakit5d._ext.mip: the loops of max_pooling_3d_mex, range_max_3d_mex and min_bbox_3d_mex
// (microscopeDataProcessing/tools/MIP/mex) on (z, y, x) C-contiguous volumes (x fastest).
// Argument checks and MATLAB's error identifiers are on the Python side
// (petakit5d.kernels.mip); the bindings check class, layout and sizes. Every result is
// independent of the thread count: each output value is computed by one thread in a fixed
// order, or (max pooling, the minimum's value) does not depend on the order at all.
#include <pk5d/common.hpp>

#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace nb = nanobind;
using u64 = uint64_t;

namespace {

template <typename T>
using in_volume = pk5d::cvolume<const T>;
template <typename T>
using out_volume = pk5d::cvolume<T>;

// ---------------------------------------------------------------------------------------
// max pooling (max_pooling_3d_mex): out[k, i, j] is the maximum of the block of pool
// (pz, py, px) at (k, i, j), the edge blocks cut at the volume's end. As the mex, the maximum
// starts at +0 (blocks of negative values give 0: a known issue kept as MATLAB has it) and
// takes a value only if it is greater (NaN never, -0 never over the starting +0). So the
// result cannot depend on the order in which a block's voxels are visited (a value greater
// than +0 has one bit pattern), and the loops visit them in memory order, not MATLAB's.
// A task is a run of output voxels along x of one (z, y) output row: the element-wise maximum
// of the block's input rows goes into a row buffer, then the maximum of each pool of it.
// ---------------------------------------------------------------------------------------

// v if greater, else m: the mex's comparison, without a branch
template <typename T>
inline T take_greater(const T m, const T v) {
    return v > m ? v : m;
}

template <typename T>
void max_pooling_3d(const T *in, T *out, u64 nz, u64 ny, u64 nx, u64 pz, u64 py, u64 px,
                    u64 mz, u64 my, u64 mx, int num_threads) {
    // output voxels along x per task: about 16 K input voxels per row
    const u64 cx = std::max<u64>(1, 16384 / px);
    const u64 nxc = (mx + cx - 1) / cx;
    const int64_t ntask = (int64_t)(nxc * my * mz);
    const u64 nbuf = px == 1 ? 0 : std::min(cx * px, nx);
#pragma omp parallel num_threads(num_threads)
    {
        std::unique_ptr<T[]> buffer(nbuf ? new T[nbuf] : nullptr);
#pragma omp for schedule(static)
        for (int64_t t = 0; t < ntask; t++) {
            const u64 xc = (u64)t % nxc;
            const u64 row = (u64)t / nxc;  // output row index y + my * z
            const u64 yo = row % my;
            const u64 zo = row / my;
            const u64 xo0 = xc * cx;
            const u64 xo1 = std::min(xo0 + cx, mx);
            const u64 x0 = xo0 * px;
            const u64 x1 = std::min(xo1 * px, nx);
            T *o = out + row * mx;
            // the row maximum goes into the output itself for px == 1
            T *r = px == 1 ? o + xo0 : buffer.get();
            const u64 n = x1 - x0;
            for (u64 i = 0; i < n; i++) r[i] = T(0);
            const u64 z1 = std::min((zo + 1) * pz, nz);
            const u64 y1 = std::min((yo + 1) * py, ny);
            for (u64 zi = zo * pz; zi < z1; zi++) {
                for (u64 yi = yo * py; yi < y1; yi++) {
                    const T *p = in + (zi * ny + yi) * nx + x0;
                    for (u64 i = 0; i < n; i++) r[i] = take_greater(r[i], p[i]);
                }
            }
            if (px == 1) continue;
            for (u64 xo = xo0; xo < xo1; xo++) {
                const u64 a = xo * px - x0;
                const u64 b = std::min((xo + 1) * px, nx) - x0;
                T m = T(0);
                for (u64 i = a; i < b; i++) m = take_greater(m, r[i]);
                o[xo] = m;
            }
        }
    }
}

// ---------------------------------------------------------------------------------------
// range maximum (range_max_3d_mex): the max along an axis (0: z, 1: y, 2: x) over
// consecutive ranges [starts[r], ends[r]). As MATLAB's max, NaN are ignored (a range of NaN
// gives its first NaN) and of equal values (0 and -0) the first one is kept, so every output
// value folds its range in increasing index along the axis, as the mex does.
// ---------------------------------------------------------------------------------------

template <typename T>
inline T max_update(const T m, const T v) {
    return v > m ? v : m;
}

// a NaN maximum is replaced by the first number
inline float max_update(const float m, const float v) {
    return (v > m || (m != m && v == v)) ? v : m;
}

inline double max_update(const double m, const double v) {
    return (v > m || (m != m && v == v)) ? v : m;
}

// the fold m = p[0], m = max_update(m, p[i]) of a run along x, in 32 independent lanes (which
// vectorize). For integers equal values have equal bits, so the lanes give the fold's result.
// For floats the lanes give the largest number (NaN never wins); its bits are the fold's
// unless it is a zero (0 or -0: the fold keeps the first) or -Inf / no number at all (all NaN:
// the fold keeps the first NaN), and then the serial fold runs.
template <typename T>
T fold_max(const T *p, u64 n) {
    constexpr int lanes = 32;  // 32: faster than 8 or 16 for uint16 and single (measured)
    auto serial = [&]() {
        T m = p[0];
        for (u64 i = 1; i < n; i++) m = max_update(m, p[i]);
        return m;
    };
    if (n < 2 * lanes) return serial();
    T acc[lanes];
    const T init = std::numeric_limits<T>::is_iec559 ? -std::numeric_limits<T>::infinity() : p[0];
    for (int k = 0; k < lanes; k++) acc[k] = init;
    const u64 nv = n / lanes * lanes;
    for (u64 i = 0; i < nv; i += lanes) {
        for (int k = 0; k < lanes; k++) acc[k] = take_greater(acc[k], p[i + k]);
    }
    for (u64 i = nv; i < n; i++) acc[0] = take_greater(acc[0], p[i]);
    // the lanes merged as a tree (vectorizes; a serial merge cost 18% of a uint16 x MIP)
    for (int w = lanes / 2; w >= 1; w /= 2) {
        for (int k = 0; k < w; k++) acc[k] = take_greater(acc[k], acc[k + w]);
    }
    const T m = acc[0];
    if constexpr (std::numeric_limits<T>::is_iec559) {
        if (m == T(0) || m == init) return serial();
    }
    return m;
}

// length of the parts a row of n values is cut into, given `units` such rows: `most`, smaller
// (down to `least`, multiples of 64) when that leaves fewer than 8 tasks per thread. Only the
// split among threads depends on it, not a result.
inline u64 task_chunk(u64 n, u64 units, int num_threads, u64 least, u64 most) {
    const u64 want = 8 * (u64)num_threads;
    const u64 parts = units >= want ? 1 : (want + units - 1) / units;
    u64 c = (n + parts - 1) / parts;
    c = (c + 63) / 64 * 64;
    return std::min(most, std::max(least, c));
}

template <typename T>
void range_max_3d(const T *im, T *out, u64 nz, u64 ny, u64 nx, int axis,
                  const std::vector<u64> &starts, const std::vector<u64> &ends,
                  int num_threads) {
    const int64_t nr = (int64_t)starts.size();
    if (axis == 2) {
        // along x: each row (z, y) is contiguous
        const int64_t nrow = (int64_t)(nz * ny);
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t c = 0; c < nrow; c++) {
            const T *row = im + c * nx;
            T *o = out + c * nr;
            for (int64_t r = 0; r < nr; r++) o[r] = fold_max(row + starts[r], ends[r] - starts[r]);
        }
    } else if (axis == 1) {
        // along y: per z plane and range, the rows of the range merged into one, in parts
        const u64 chunk = task_chunk(nx, nz * nr, num_threads, 256, 4096);
        const u64 nxc = (nx + chunk - 1) / chunk;
        const int64_t ntask = (int64_t)(nz * nr * nxc);
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t t = 0; t < ntask; t++) {
            const u64 xc = t % nxc;
            const u64 r = (t / nxc) % nr;
            const u64 z = t / (nxc * nr);
            const u64 x0 = xc * chunk;
            const u64 x1 = x0 + chunk < nx ? x0 + chunk : nx;
            T *o = out + (z * nr + r) * nx;
            const T *p = im + (z * ny + starts[r]) * nx;
            for (u64 x = x0; x < x1; x++) o[x] = p[x];
            for (u64 y = starts[r] + 1; y < ends[r]; y++) {
                p = im + (z * ny + y) * nx;
                for (u64 x = x0; x < x1; x++) o[x] = max_update(o[x], p[x]);
            }
        }
    } else {
        // along z: per range, the planes of the range merged into one, in parts (the mex's)
        const u64 nyx = ny * nx;
        const u64 chunk = 16384;
        const u64 nc = (nyx + chunk - 1) / chunk;
        const int64_t ntask = (int64_t)(nr * nc);
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t t = 0; t < ntask; t++) {
            const u64 ci = t % nc;
            const u64 r = t / nc;
            const u64 i0 = ci * chunk;
            const u64 i1 = i0 + chunk < nyx ? i0 + chunk : nyx;
            T *o = out + r * nyx;
            const T *p = im + starts[r] * nyx;
            for (u64 i = i0; i < i1; i++) o[i] = p[i];
            for (u64 z = starts[r] + 1; z < ends[r]; z++) {
                p = im + z * nyx;
                for (u64 i = i0; i < i1; i++) o[i] = max_update(o[i], p[i]);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------
// minimum in a box (min_bbox_3d_mex): starting from the largest finite value of the class, a
// value is taken only if it is smaller, so NaN are skipped and an empty box gives that
// value. The mex runs serially in MATLAB's order. The minimum's value does not depend on the
// order, so here contiguous parts of the box's (z, y) rows run in parallel, each row in 16
// independent lanes (which vectorize), and the parts' minimums are merged. Of equal values
// only 0 and -0 differ in bits: the serial loop keeps the first zero in MATLAB's order (y
// fastest, then x, then z), so a zero minimum is replaced by that zero.
// ---------------------------------------------------------------------------------------

template <typename T>
inline T take_smaller(const T m, const T v) {
    return v < m ? v : m;
}

template <typename T>
T min_box_3d(const T *im, u64 ny, u64 nx, u64 z0, u64 y0, u64 x0, u64 z1, u64 y1, u64 x1,
             int num_threads) {
    constexpr int lanes = 16;
    const T start = std::numeric_limits<T>::max();
    const u64 by = y1 - y0, bz = z1 - z0;
    const u64 nrow = (x1 > x0) ? by * bz : 0;
    if (nrow == 0) return start;
    const u64 nparts = std::min<u64>((u64)num_threads, nrow);
    // (not std::vector: vector<bool> packs bits, which threads cannot write apart)
    std::unique_ptr<T[]> part(new T[nparts]);
    const u64 n = x1 - x0;
    const u64 nv = n / lanes * lanes;
#pragma omp parallel for schedule(static, 1) num_threads((int)nparts)
    for (int64_t q = 0; q < (int64_t)nparts; q++) {
        const u64 c0 = nrow * (u64)q / nparts;
        const u64 c1 = nrow * ((u64)q + 1) / nparts;
        T acc[lanes];
        for (int k = 0; k < lanes; k++) acc[k] = start;
        for (u64 c = c0; c < c1; c++) {
            const u64 y = y0 + c % by;
            const u64 z = z0 + c / by;
            const T *p = im + (z * ny + y) * nx + x0;
            for (u64 i = 0; i < nv; i += lanes) {
                for (int k = 0; k < lanes; k++) acc[k] = take_smaller(acc[k], p[i + k]);
            }
            for (u64 i = nv; i < n; i++) acc[0] = take_smaller(acc[0], p[i]);
        }
        T m = start;
        for (int k = 0; k < lanes; k++) m = take_smaller(m, acc[k]);
        part[q] = m;
    }
    T m = start;
    for (u64 q = 0; q < nparts; q++) m = take_smaller(m, part[q]);
    if constexpr (std::numeric_limits<T>::is_iec559) {
        if (m == T(0)) {
            // the first zero in MATLAB's order (strided; only for a zero minimum)
            for (u64 z = z0; z < z1; z++) {
                for (u64 x = x0; x < x1; x++) {
                    for (u64 y = y0; y < y1; y++) {
                        const T v = im[(z * ny + y) * nx + x];
                        if (v == T(0)) return v;
                    }
                }
            }
        }
    }
    return m;
}

// ---------------------------------------------------------------------------------------
// bindings
// ---------------------------------------------------------------------------------------

template <typename T>
void bind_max_pooling(in_volume<T> orig, out_volume<T> out, u64 pz, u64 py, u64 px,
                      int num_threads) {
    pk5d::check_threads(num_threads);
    if (pz < 1 || py < 1 || px < 1) throw std::invalid_argument("pool sizes must be at least 1");
    const u64 nz = orig.shape(0), ny = orig.shape(1), nx = orig.shape(2);
    const u64 mz = nz / pz + (nz % pz != 0), my = ny / py + (ny % py != 0),
              mx = nx / px + (nx % px != 0);
    if (out.shape(0) != mz || out.shape(1) != my || out.shape(2) != mx)
        throw std::invalid_argument("out has the wrong size");
    if (mz * my * mx == 0) return;
    const T *src = orig.data();
    T *dst = out.data();
    nb::gil_scoped_release release;
    max_pooling_3d<T>(src, dst, nz, ny, nx, pz, py, px, mz, my, mx, num_threads);
}

template <typename T>
void bind_range_max(in_volume<T> im, out_volume<T> out, int axis,
                    const std::vector<u64> &starts, int num_threads) {
    pk5d::check_threads(num_threads);
    if (axis < 0 || axis > 2) throw std::invalid_argument("axis must be 0, 1 or 2");
    const u64 sz[3] = {(u64)im.shape(0), (u64)im.shape(1), (u64)im.shape(2)};
    if (sz[0] * sz[1] * sz[2] == 0) throw std::invalid_argument("im must not be empty");
    const u64 n = sz[axis];
    const size_t nr = starts.size();
    if (nr == 0 || starts[0] != 0) throw std::invalid_argument("starts must begin with 0");
    for (size_t r = 1; r < nr; r++) {
        if (starts[r] <= starts[r - 1]) throw std::invalid_argument("starts must increase");
    }
    if (starts[nr - 1] >= n) throw std::invalid_argument("starts must be below the size along axis");
    u64 osz[3] = {sz[0], sz[1], sz[2]};
    osz[axis] = nr;
    if ((u64)out.shape(0) != osz[0] || (u64)out.shape(1) != osz[1] || (u64)out.shape(2) != osz[2])
        throw std::invalid_argument("out has the wrong size");
    std::vector<u64> ends(nr);
    for (size_t r = 0; r < nr; r++) ends[r] = r + 1 < nr ? starts[r + 1] : n;
    const T *src = im.data();
    T *dst = out.data();
    nb::gil_scoped_release release;
    range_max_3d<T>(src, dst, sz[0], sz[1], sz[2], axis, starts, ends, num_threads);
}

template <typename T>
T bind_min_box(in_volume<T> im, u64 z0, u64 y0, u64 x0, u64 z1, u64 y1, u64 x1,
               int num_threads) {
    pk5d::check_threads(num_threads);
    const u64 nz = im.shape(0), ny = im.shape(1), nx = im.shape(2);
    if (z1 < z0 || y1 < y0 || x1 < x0 || z1 > nz || y1 > ny || x1 > nx)
        throw std::invalid_argument("box is not valid");
    const T *src = im.data();
    nb::gil_scoped_release release;
    return min_box_3d<T>(src, ny, nx, z0, y0, x0, z1, y1, x1, num_threads);
}

template <typename T>
void def_all(nb::module_ &m, bool range_max) {
    using namespace nb::literals;
    m.def("max_pooling_3d", &bind_max_pooling<T>, "orig"_a.noconvert(), "out"_a.noconvert(),
          "pz"_a, "py"_a, "px"_a, "num_threads"_a,
          "Max pooling of orig (nz, ny, nx) into out (ceil(size / pool)), blocks of pool "
          "(pz, py, px); the maximum starts at 0 as max_pooling_3d_mex's.");
    m.def("min_box_3d", &bind_min_box<T>, "im"_a.noconvert(), "z0"_a, "y0"_a, "x0"_a, "z1"_a,
          "y1"_a, "x1"_a, "num_threads"_a,
          "Minimum of im[z0:z1, y0:y1, x0:x1], starting from the largest finite value of the "
          "class (NaN skipped), as min_bbox_3d_mex; of zeros the first in MATLAB's order.");
    if (range_max)
        m.def("range_max_3d", &bind_range_max<T>, "im"_a.noconvert(), "out"_a.noconvert(),
              "axis"_a, "starts"_a, "num_threads"_a,
              "Max of im (nz, ny, nx) along axis over the ranges that begin at starts "
              "(0-based, increasing, from 0), as range_max_3d_mex.");
}

}  // namespace

NB_MODULE(mip, m) {
    m.doc() = "Max pooling, range maximum and box minimum kernels (MIPs), on (z, y, x) volumes";
    def_all<bool>(m, false);
    def_all<uint8_t>(m, true);
    def_all<uint16_t>(m, true);
    def_all<float>(m, true);
    def_all<double>(m, true);
}
