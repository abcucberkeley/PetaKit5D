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

// helpers inlined into each instruction set's kernels (kernels.inc)
#if defined(__GNUC__) || defined(__clang__)
#define PK5D_MIP_INLINE inline __attribute__((always_inline))
#else
#define PK5D_MIP_INLINE __forceinline
#endif

// Run-time dispatch on x86-64 with GCC or Clang: the kernels are also compiled for SSE4.1
// (pmaxuw/pminuw: uint16 max and min in one instruction, which SSE2 lacks) and AVX2, and the
// best set the CPU has is chosen once. Other compilers and processors: the baseline only.
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
#define PK5D_MIP_DISPATCH 1
#endif

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
PK5D_MIP_INLINE T take_greater(const T m, const T v) {
    return v > m ? v : m;
}

// ---------------------------------------------------------------------------------------
// range maximum (range_max_3d_mex): the max along an axis (0: z, 1: y, 2: x) over
// consecutive ranges [starts[r], ends[r]). As MATLAB's max, NaN are ignored (a range of NaN
// gives its first NaN) and of equal values (0 and -0) the first one is kept, so every output
// value folds its range in increasing index along the axis, as the mex does.
// ---------------------------------------------------------------------------------------

template <typename T>
PK5D_MIP_INLINE T max_update(const T m, const T v) {
    return v > m ? v : m;
}

// a NaN maximum is replaced by the first number
PK5D_MIP_INLINE float max_update(const float m, const float v) {
    return (v > m || (m != m && v == v)) ? v : m;
}

PK5D_MIP_INLINE double max_update(const double m, const double v) {
    return (v > m || (m != m && v == v)) ? v : m;
}

// the fold m = p[0], m = max_update(m, p[i]) of a run along x, in 32 independent lanes (which
// vectorize). For integers equal values have equal bits, so the lanes give the fold's result.
// For floats the lanes give the largest number (NaN never wins); its bits are the fold's
// unless it is a zero (0 or -0: the fold keeps the first) or -Inf / no number at all (all NaN:
// the fold keeps the first NaN), and then the serial fold runs.
template <typename T>
PK5D_MIP_INLINE T fold_max(const T *p, u64 n) {
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
PK5D_MIP_INLINE u64 task_chunk(u64 n, u64 units, int num_threads, u64 least, u64 most) {
    const u64 want = 8 * (u64)num_threads;
    const u64 parts = units >= want ? 1 : (want + units - 1) / units;
    u64 c = (n + parts - 1) / parts;
    c = (c + 63) / 64 * 64;
    return std::min(most, std::max(least, c));
}

// ---------------------------------------------------------------------------------------
// minimum in a box (min_bbox_3d_mex): starting from the largest finite value of the class, a
// value is taken only if it is smaller, so NaN are skipped and an empty box gives that
// value. The mex runs serially in MATLAB's order. The minimum's value does not depend on the
// order, so here contiguous parts of the box's (z, y) rows run in parallel, each row in 32
// independent lanes (which vectorize), and the parts' minimums are merged. Of equal values
// only 0 and -0 differ in bits: the serial loop keeps the first zero in MATLAB's order (y
// fastest, then x, then z), so a zero minimum is replaced by that zero.
// ---------------------------------------------------------------------------------------

template <typename T>
PK5D_MIP_INLINE T take_smaller(const T m, const T v) {
    return v < m ? v : m;
}
// ---------------------------------------------------------------------------------------
// the kernels, once per instruction set
// ---------------------------------------------------------------------------------------

#define PK5D_MIP_NS isa_baseline
#define PK5D_MIP_ISA
#include "kernels.inc"
#undef PK5D_MIP_NS
#undef PK5D_MIP_ISA
#ifdef PK5D_MIP_DISPATCH
#define PK5D_MIP_NS isa_sse41
#define PK5D_MIP_ISA __attribute__((target("sse4.1")))
#include "kernels.inc"
#undef PK5D_MIP_NS
#undef PK5D_MIP_ISA
#define PK5D_MIP_NS isa_avx2
#define PK5D_MIP_ISA __attribute__((target("avx2")))
#include "kernels.inc"
#undef PK5D_MIP_NS
#undef PK5D_MIP_ISA
#endif

// instruction sets: 0 baseline, 1 SSE4.1, 2 AVX2
int best_isa() {
#ifdef PK5D_MIP_DISPATCH
    static const int best = [] {
        __builtin_cpu_init();
        if (__builtin_cpu_supports("avx2")) return 2;
        if (__builtin_cpu_supports("sse4.1")) return 1;
        return 0;
    }();
    return best;
#else
    return 0;
#endif
}

// the set a call runs: -1 the best one, else the one asked for (refused if the CPU lacks it)
int pick_isa(int isa) {
    if (isa < 0) return best_isa();
    if (isa > best_isa()) throw std::invalid_argument("this CPU lacks that instruction set");
    return isa;
}

#ifdef PK5D_MIP_DISPATCH
#define PK5D_MIP_CALL(isa, fn, ...)                                      \
    ((isa) == 2   ? isa_avx2::fn(__VA_ARGS__)                          \
     : (isa) == 1 ? isa_sse41::fn(__VA_ARGS__)                         \
                  : isa_baseline::fn(__VA_ARGS__))
#else
#define PK5D_MIP_CALL(isa, fn, ...) isa_baseline::fn(__VA_ARGS__)
#endif

// ---------------------------------------------------------------------------------------
// bindings
// ---------------------------------------------------------------------------------------

template <typename T>
void bind_max_pooling(in_volume<T> orig, out_volume<T> out, u64 pz, u64 py, u64 px,
                      int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
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
    PK5D_MIP_CALL(isa, max_pooling_3d<T>, src, dst, nz, ny, nx, pz, py, px, mz, my, mx,
                  num_threads);
}

template <typename T>
void bind_range_max(in_volume<T> im, out_volume<T> out, int axis,
                    const std::vector<u64> &starts, int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
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
    PK5D_MIP_CALL(isa, range_max_3d<T>, src, dst, sz[0], sz[1], sz[2], axis, starts, ends,
                  num_threads);
}

template <typename T>
T bind_min_box(in_volume<T> im, u64 z0, u64 y0, u64 x0, u64 z1, u64 y1, u64 x1,
               int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    const u64 nz = im.shape(0), ny = im.shape(1), nx = im.shape(2);
    if (z1 < z0 || y1 < y0 || x1 < x0 || z1 > nz || y1 > ny || x1 > nx)
        throw std::invalid_argument("box is not valid");
    const T *src = im.data();
    nb::gil_scoped_release release;
    return PK5D_MIP_CALL(isa, min_box_3d<T>, src, ny, nx, z0, y0, x0, z1, y1, x1, num_threads);
}

template <typename T>
void def_all(nb::module_ &m, bool range_max) {
    using namespace nb::literals;
    m.def("max_pooling_3d", &bind_max_pooling<T>, "orig"_a.noconvert(), "out"_a.noconvert(),
          "pz"_a, "py"_a, "px"_a, "num_threads"_a, "isa"_a = -1,
          "Max pooling of orig (nz, ny, nx) into out (ceil(size / pool)), blocks of pool "
          "(pz, py, px); the maximum starts at 0 as max_pooling_3d_mex's.");
    m.def("min_box_3d", &bind_min_box<T>, "im"_a.noconvert(), "z0"_a, "y0"_a, "x0"_a, "z1"_a,
          "y1"_a, "x1"_a, "num_threads"_a, "isa"_a = -1,
          "Minimum of im[z0:z1, y0:y1, x0:x1], starting from the largest finite value of the "
          "class (NaN skipped), as min_bbox_3d_mex; of zeros the first in MATLAB's order.");
    if (range_max)
        m.def("range_max_3d", &bind_range_max<T>, "im"_a.noconvert(), "out"_a.noconvert(),
              "axis"_a, "starts"_a, "num_threads"_a, "isa"_a = -1,
              "Max of im (nz, ny, nx) along axis over the ranges that begin at starts "
              "(0-based, increasing, from 0), as range_max_3d_mex.");
}

}  // namespace

NB_MODULE(mip, m) {
    m.doc() = "Max pooling, range maximum and box minimum kernels (MIPs), on (z, y, x) volumes";
    m.def("best_isa", &best_isa,
          "Instruction set the kernels run (0 baseline, 1 SSE4.1, 2 AVX2), chosen once from the "
          "CPU; each kernel's isa argument asks for another one (tests).");
    def_all<bool>(m, false);
    def_all<uint8_t>(m, true);
    def_all<uint16_t>(m, true);
    def_all<float>(m, true);
    def_all<double>(m, true);
}
