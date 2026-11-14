// petakit5d._ext.unmix: the loop of unmix_channels_sum_mex (microscopeDataProcessing/unmix/mex)
// on 1-D contiguous views of the channels (any shape: the sum is per element, so the layout
// does not matter). Argument checks and MATLAB's error identifiers are on the Python side
// (petakit5d.tools.unmix); the bindings check classes and sizes.
//
// Each element is the mex's arithmetic, in the mex's order:
//   acc = im0 (single) or 0
//   for each channel c: acc = acc + (float)((double)(float)((double)(float)x_c - b_c) * f_c),
//     the background step only when b_c != 0 (x - 0 == x bit for bit), each single-double
//     operation in double and rounded to single once (MATLAB's rule), the sum in single
//   finish: v = acc > 0 ? acc : 0 (max(0, NaN) and max(0, -0) are +0), then
//     (float)((double)v + offset) when offset != 0, then MATLAB's cast: round half away from
//     zero and saturate for uint8 and uint16 (NaN gives 0), the value for single and double
//   partial (no finish): acc as it is, in single
// No fused multiply-adds (the group is built with -ffp-contract=off); a NaN that a step makes
// from numbers is x86's (sign bit set) on every processor. The elements are cut into blocks of
// 8192 (the sum of a block stays in the cache), shared out contiguously between the threads; no
// element depends on another, so results do not depend on the thread count or the instruction
// set (the loops are also compiled for AVX2, chosen at run time on x86-64).
#include <pk5d/common.hpp>
#include <pk5d/fpexact.hpp>

#include "ops_bind.hpp"

#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace nb = nanobind;
using u64 = uint64_t;

#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
#define PK5D_UX_DISPATCH 1
#endif

namespace {

// elements per block of work (as the mex)
constexpr u64 kBlock = 8192;

using flat_in = nb::ndarray<nb::ro, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using flat_out = nb::ndarray<nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using dvec = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using fvec = nb::ndarray<const float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using fout = nb::ndarray<float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

enum class cls { u8, u16, f32, f64 };

cls class_of(nb::dlpack::dtype d, const char *what) {
    if (d == nb::dtype<uint8_t>()) return cls::u8;
    if (d == nb::dtype<uint16_t>()) return cls::u16;
    if (d == nb::dtype<float>()) return cls::f32;
    if (d == nb::dtype<double>()) return cls::f64;
    throw std::invalid_argument(std::string(what) + ": uint8, uint16, single or double only");
}

struct channel {
    const void *data;
    cls c;
    double factor;
    double background;
};

}  // namespace

// ---------------------------------------------------------------------------------------
// the loops, once per instruction set
// ---------------------------------------------------------------------------------------

namespace {
#define PK5D_UX_NS isa_baseline
#define PK5D_UX_ISA
#include "kernels.inc"
#undef PK5D_UX_NS
#undef PK5D_UX_ISA
#ifdef PK5D_UX_DISPATCH
#define PK5D_UX_NS isa_avx2
#define PK5D_UX_ISA __attribute__((target("avx2")))
#include "kernels.inc"
#undef PK5D_UX_NS
#undef PK5D_UX_ISA
#endif

// instruction sets: 0 baseline, 1 AVX2
int best_isa() {
#ifdef PK5D_UX_DISPATCH
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

template <typename O>
void run(const std::vector<channel> &chans, const float *im0, O *out, const u64 n,
         const double offset, const bool finish, const int nt, const bool bind, const int isa) {
    const u64 nblocks = (n + kBlock - 1) / kBlock;
    const int team = (int)std::min<u64>((u64)nt, nblocks);
    void (*fn)(const channel *, u64, const float *, O *, u64, u64, u64, double, bool) =
        isa_baseline::blocks<O>;
#ifdef PK5D_UX_DISPATCH
    if (isa == 1) fn = isa_avx2::blocks<O>;
#else
    (void)isa;
#endif
    pk5d::ops::parallel(team, bind, [&](int k, int nk) {
        std::ptrdiff_t lo, hi;
        pk5d::ops::share((std::ptrdiff_t)nblocks, k, nk, lo, hi);
        fn(chans.data(), chans.size(), im0, out, n, (u64)lo, (u64)hi, offset, finish);
    });
}

// out = finish ? cast(max(0, im0 + sum_c ...) + offset) : im0 + sum_c ...; out may be im0
void unmix_sum(const std::vector<flat_in> &channels, dvec factors, dvec backgrounds,
               double offset, bool finish, std::optional<fvec> im0, flat_out out, int num_threads,
               bool bind, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    const u64 nc = channels.size();
    if (factors.shape(0) != nc || backgrounds.shape(0) != nc)
        throw std::invalid_argument("one factor and one background per channel");
    if (nc == 0 && !im0) throw std::invalid_argument("no channels and no partial sum");
    const u64 n = out.shape(0);
    std::vector<channel> chans(nc);
    for (u64 c = 0; c < nc; c++) {
        if (channels[c].shape(0) != n)
            throw std::invalid_argument("the channels must have the size of out");
        chans[c] = {channels[c].data(), class_of(channels[c].dtype(), "channel"),
                    factors.data()[c], backgrounds.data()[c]};
    }
    const float *p0 = nullptr;
    if (im0) {
        if (im0->shape(0) != n) throw std::invalid_argument("im0 must have the size of out");
        p0 = im0->data();
    }
    const cls oc = class_of(out.dtype(), "out");
    if (!finish && oc != cls::f32) throw std::invalid_argument("a partial sum is single");
    if (n == 0) return;
    void *o = out.data();
    nb::gil_scoped_release release;
    switch (oc) {
        case cls::u8:
            run(chans, p0, (uint8_t *)o, n, offset, finish, num_threads, bind, isa);
            break;
        case cls::u16:
            run(chans, p0, (uint16_t *)o, n, offset, finish, num_threads, bind, isa);
            break;
        case cls::f32:
            run(chans, p0, (float *)o, n, offset, finish, num_threads, bind, isa);
            break;
        default:
            run(chans, p0, (double *)o, n, offset, finish, num_threads, bind, isa);
            break;
    }
}

// out = x - b for single x and a double b (MATLAB's single minus double); out may be x
void subtract(fvec x, double b, fout out, int num_threads, bool bind, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    const u64 n = x.shape(0);
    if (out.shape(0) != n) throw std::invalid_argument("out must have the size of x");
    if (n == 0) return;
    const float *px = x.data();
    float *po = out.data();
    const u64 nblocks = (n + kBlock - 1) / kBlock;
    const int team = (int)std::min<u64>((u64)num_threads, nblocks);
    void (*fn)(const float *, float *, double, u64, u64) = isa_baseline::subtract;
#ifdef PK5D_UX_DISPATCH
    if (isa == 1) fn = isa_avx2::subtract;
#endif
    nb::gil_scoped_release release;
    pk5d::ops::parallel(team, bind, [&](int k, int nk) {
        std::ptrdiff_t lo, hi;
        pk5d::ops::share((std::ptrdiff_t)nblocks, k, nk, lo, hi);
        fn(px, po, b, (u64)lo * kBlock, std::min(n, (u64)hi * kBlock));
    });
}

}  // namespace

NB_MODULE(unmix, m) {
    using namespace nb::literals;
    m.doc() = "The weighted channel sum of unmix_channels_sum_mex on contiguous arrays";
    m.def("best_isa", &best_isa,
          "Instruction set the kernels run (0 baseline, 1 AVX2), chosen once from the CPU; the "
          "isa argument asks for another one (tests).");
    m.def("unmix_sum", &unmix_sum, "channels"_a.noconvert(), "factors"_a.noconvert(),
          "backgrounds"_a.noconvert(), "offset"_a, "finish"_a, "im0"_a.noconvert().none(),
          "out"_a.noconvert(), "num_threads"_a, "bind"_a = true, "isa"_a = -1,
          "out = cast(max(0, im0 + sum_c (single(x_c) - b_c) * f_c) + offset) per element (finish),\n"
          "or the partial sum im0 + sum_c ... in single, in unmix_channels_sum_mex's arithmetic;\n"
          "channels: 1-D contiguous uint8, uint16, single or double arrays of out's size; im0:\n"
          "None (0) or a single array of out's size (it may be out); out: uint8, uint16, single\n"
          "or double (single without finish).");
    m.def("subtract", &subtract, "x"_a.noconvert(), "b"_a, "out"_a.noconvert(), "num_threads"_a,
          "bind"_a = true, "isa"_a = -1,
          "out = x - b for a single x and a double b, in double and rounded to single (MATLAB's\n"
          "single minus double); out may be x.");
}
