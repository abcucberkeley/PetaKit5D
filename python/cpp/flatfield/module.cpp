// petakit5d._ext.flatfield: the loop of flat_field_correction_3d_mex
// (microscopeDataProcessing/flatfield_correction/mex) on a (z, y, x) C-contiguous volume with
// (y, x) flat-field and background images. Argument checks and MATLAB's error identifiers are
// on the Python side (petakit5d.kernels.flatfield); the bindings check class, layout and sizes.
//
// Each output voxel is the mex's arithmetic in the mex's order of operations:
//   f = (float(in) - b) / flat; f = f >= 0 ? f : 0; out = round_output(f + offset or b)
// in single precision, where round_output is (U)(f + 0.5) for uint16 (computed in double and
// converted as the mex's GCC build converts: cvttsd2si into 32 bits, then the low 16 bits, so
// values beyond 65535 wrap, negative ones wrap from 65535 down, and NaN and values outside the
// int32 range give 0) and no rounding for single. The value of a voxel depends on nothing
// else, so the loops run in memory order, blocked so that the rows of the flat field and the
// background are read once per block of planes, and results do not depend on the thread count
// or the instruction set.
#include <pk5d/common.hpp>
#include <pk5d/fpexact.hpp>

#include "ops_bind.hpp"

#include <algorithm>
#include <cstdint>
#include <type_traits>

namespace nb = nanobind;
using u64 = uint64_t;

#if defined(__GNUC__) || defined(__clang__)
#define PK5D_FF_INLINE inline __attribute__((always_inline))
#else
#define PK5D_FF_INLINE __forceinline
#endif

// Run-time dispatch on x86-64 with GCC or Clang: the loops are also compiled for AVX2 (8
// voxels per instruction; no FMA, which the policy flags rule out anyway) and the best set the
// CPU has is chosen once. Other compilers and processors: the baseline only.
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
#define PK5D_FF_DISPATCH 1
#endif

namespace {

// the mex's round_output: (U)(value + 0.5) for uint16, the value itself for single
template <typename U>
PK5D_FF_INLINE U round_output(float v) {
    if constexpr (std::is_same_v<U, uint16_t>) {
        // (double)v + 0.5 is exact; cvttsd2si into a 32-bit register gives the truncation for
        // values in (-2^31 - 1, 2^31) and 0x80000000 for the others and NaN; the mex keeps the
        // low 16 bits
        double d = static_cast<double>(v) + 0.5;
        d = (d > -2147483649.0 && d < 2147483648.0) ? d : -2147483648.0;
        return static_cast<uint16_t>(static_cast<uint32_t>(static_cast<int32_t>(d)));
    } else {
        return v;
    }
}

// one call: the arrays, the sizes and the task grid (blocks of `rows` rows of y by chunks of
// `zc` planes of z; nyb blocks per chunk)
template <typename T, typename U>
struct ff_job {
    const T *in;
    const float *flat;
    const float *bg;
    U *out;
    u64 nz, ny, nx;
    bool has_offset;
    float offset;
    u64 rows, nyb, zc, ntask;
};

}  // namespace

// ---------------------------------------------------------------------------------------
// the loops, once per instruction set
// ---------------------------------------------------------------------------------------

namespace {
#define PK5D_FF_NS isa_baseline
#define PK5D_FF_ISA
#include "kernels.inc"
#undef PK5D_FF_NS
#undef PK5D_FF_ISA
#ifdef PK5D_FF_DISPATCH
#define PK5D_FF_NS isa_avx2
#define PK5D_FF_ISA __attribute__((target("avx2")))
#include "kernels.inc"
#undef PK5D_FF_NS
#undef PK5D_FF_ISA
#endif

// instruction sets: 0 baseline, 1 AVX2
int best_isa() {
#ifdef PK5D_FF_DISPATCH
    static const int best = [] {
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2") ? 1 : 0;
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

// rows of y per block: about 64 K voxels (256 KB of flat field and as much of background, in
// a core's L2 cache; runs of 128 KB or more of data per plane, which the prefetchers follow),
// at least one row. 8 K and whole planes measured within 5% of it on a 512 x 1800 plane.
constexpr u64 kBlockVoxels = 65536;

// the task grid: z cut into at least one chunk per thread (more, up to one plane each, when
// that leaves fewer than 4 tasks per thread), tasks ordered chunk by chunk, so that each thread
// gets contiguous planes (as the mex's static schedule: a thread's pages are its own, which
// matters when they are first touched) and reads each block of the flat field once per chunk
template <typename T, typename U>
void plan(ff_job<T, U> &j, int nt) {
    j.rows = std::max<u64>(1, std::min<u64>(j.ny, kBlockVoxels / std::max<u64>(1, j.nx)));
    j.nyb = (j.ny + j.rows - 1) / j.rows;
    const u64 want = 4 * (u64)nt;
    u64 parts = std::max<u64>((u64)nt, (want + j.nyb - 1) / j.nyb);
    parts = std::min<u64>(parts, j.nz);
    j.zc = (j.nz + parts - 1) / parts;
    j.ntask = j.nyb * ((j.nz + j.zc - 1) / j.zc);
}

template <typename T, typename U>
void bind_correct(pk5d::cvolume<const T> im, pk5d::cimage<const float> flat,
                  pk5d::cimage<const float> bg, pk5d::cvolume<U> out, bool has_offset,
                  float offset, int num_threads, int isa, bool bind) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    const u64 nz = im.shape(0), ny = im.shape(1), nx = im.shape(2);
    if (flat.shape(0) != ny || flat.shape(1) != nx)
        throw std::invalid_argument("flat must have the size of a plane of im");
    if (bg.shape(0) != ny || bg.shape(1) != nx)
        throw std::invalid_argument("background must have the size of a plane of im");
    if (out.shape(0) != nz || out.shape(1) != ny || out.shape(2) != nx)
        throw std::invalid_argument("out must have the size of im");
    if (nz * ny * nx == 0) return;
    ff_job<T, U> j{im.data(), flat.data(), bg.data(), out.data(), nz, ny, nx, has_offset,
                   offset, 0, 0, 0, 0};
    plan(j, num_threads);
    const int nt = (int)std::min<u64>((u64)num_threads, j.ntask);
    void (*fn)(const ff_job<T, U> &, u64, u64) = isa_baseline::ff_tasks<T, U>;
#ifdef PK5D_FF_DISPATCH
    if (isa == 1) fn = isa_avx2::ff_tasks<T, U>;
#endif
    nb::gil_scoped_release release;
    pk5d::ops::parallel(nt, bind, [&](int k, int nk) {
        std::ptrdiff_t lo, hi;
        pk5d::ops::share((std::ptrdiff_t)j.ntask, k, nk, lo, hi);
        fn(j, (u64)lo, (u64)hi);
    });
}

template <typename T, typename U>
void def_correct(nb::module_ &m) {
    using namespace nb::literals;
    m.def("flat_field_correct", &bind_correct<T, U>, "im"_a.noconvert(), "flat"_a.noconvert(),
          "background"_a.noconvert(), "out"_a.noconvert(), "has_offset"_a, "offset"_a,
          "num_threads"_a, "isa"_a = -1, "bind"_a = true,
          "out = max((im - background) / flat, 0) + (offset if has_offset else background) "
          "per voxel of im (nz, ny, nx) with the (ny, nx) single images, in single precision; "
          "uint16 out rounds and wraps as flat_field_correction_3d_mex.");
}

}  // namespace

NB_MODULE(flatfield, m) {
    m.doc() = "Flat-field correction kernel (flat_field_correction_3d_mex) on (z, y, x) volumes";
    m.def("best_isa", &best_isa,
          "Instruction set the kernel runs (0 baseline, 1 AVX2), chosen once from the CPU; the "
          "isa argument asks for another one (tests).");
    def_correct<uint16_t, uint16_t>(m);
    def_correct<uint16_t, float>(m);
    def_correct<float, float>(m);
}
