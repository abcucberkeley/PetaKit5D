// petakit5d._ext.filters: the arithmetic of MATLAB's linear filters on C-contiguous arrays
// (petakit5d.matlabcompat.filters builds the tap orders, measured on MATLAB R2024b):
// - line_program: a filter along one axis of an array viewed as (outer, n, inner), whose output
//   row i runs a tap program (rowmap[i] selects it): fused or plain multiply-adds in a given tap
//   order, in one or more chains that are then added in order; the input is padded along the
//   axis as MATLAB's padarray pads (constant, replicate, symmetric, circular).
// - taps_nd: an N-D correlation of a padded volume with a list of taps in a given order (one
//   chain per output voxel, fused or plain).
// Every output value is computed by one thread in a fixed order, so the results do not depend
// on the thread count.
#include <pk5d/common.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace nb = nanobind;
using i64 = int64_t;

// Run-time dispatch on x86-64 with GCC or Clang: the kernels are also compiled with FMA (and
// AVX2), where std::fma is one instruction; the baseline calls the C library's fma. Other
// compilers and processors (ARM has the instruction): the baseline only.
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
#define PK5D_FILT_DISPATCH 1
#endif

namespace {

// MATLAB's padarray: the input index of padded index q (any integer) of a line of n values,
// or -1 for the constant
inline i64 pad_index(i64 q, i64 n, int mode) {
    if (q >= 0 && q < n) return q;
    switch (mode) {
        case 1:  // replicate
            return q < 0 ? 0 : n - 1;
        case 2: {  // symmetric: ... 1 0 | 0 1 ... n-1 | n-1 n-2 ...
            i64 r = q % (2 * n);
            if (r < 0) r += 2 * n;
            return r < n ? r : 2 * n - 1 - r;
        }
        case 3: {  // circular
            i64 r = q % n;
            return r < 0 ? r + n : r;
        }
        default:
            return -1;
    }
}

// program codes of the edge taps (python/src/petakit5d/matlabcompat/filters.py)
constexpr i64 kEdgeLow = i64(1) << 21;
constexpr i64 kEdgeHigh = i64(1) << 22;

}  // namespace

#define PK5D_FILT_NS isa_baseline
#define PK5D_FILT_ISA
#include "kernels.inc"
#undef PK5D_FILT_NS
#undef PK5D_FILT_ISA
#ifdef PK5D_FILT_DISPATCH
#define PK5D_FILT_NS isa_fma
#define PK5D_FILT_ISA __attribute__((target("avx2,fma")))
#include "kernels.inc"
#undef PK5D_FILT_NS
#undef PK5D_FILT_ISA
#endif

namespace {

// instruction sets: 0 baseline, 1 AVX2 with FMA
int best_isa() {
#ifdef PK5D_FILT_DISPATCH
    static const int best = [] {
        __builtin_cpu_init();
        return (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) ? 1 : 0;
    }();
    return best;
#else
    return 0;
#endif
}

int pick_isa(int isa) {
    if (isa < 0) return best_isa();
    if (isa > best_isa()) throw std::invalid_argument("this CPU lacks that instruction set");
    return isa;
}

#ifdef PK5D_FILT_DISPATCH
#define PK5D_FILT_CALL(isa, fn, ...) \
    ((isa) == 1 ? isa_fma::fn(__VA_ARGS__) : isa_baseline::fn(__VA_ARGS__))
#else
#define PK5D_FILT_CALL(isa, fn, ...) isa_baseline::fn(__VA_ARGS__)
#endif

template <typename T>
using in_array = nb::ndarray<const T, nb::c_contig, nb::device::cpu>;
template <typename T>
using out_array = nb::ndarray<T, nb::c_contig, nb::device::cpu>;
using i32_array = nb::ndarray<const int32_t, nb::c_contig, nb::device::cpu>;
using i64_array = nb::ndarray<const int64_t, nb::c_contig, nb::device::cpu>;

template <typename T>
void bind_line_program(in_array<T> src, out_array<T> dst, i64 outer, i64 n_in, i64 inner, i64 n_out,
                       in_array<T> w, i32_array prog, i32_array rowmap, i64 pre, i64 span, int mode,
                       double value, int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    if (outer < 0 || n_in < 0 || inner < 0 || n_out < 0)
        throw std::invalid_argument("negative sizes");
    if ((i64)src.size() != outer * n_in * inner) throw std::invalid_argument("src has the wrong size");
    if ((i64)dst.size() != outer * n_out * inner) throw std::invalid_argument("dst has the wrong size");
    if (prog.ndim() != 2) throw std::invalid_argument("prog must be 2-D (programs x length)");
    if ((i64)rowmap.size() != n_out) throw std::invalid_argument("rowmap needs one entry per output row");
    if (mode < 0 || mode > 3) throw std::invalid_argument("mode is 0 to 3");
    if (n_in == 0 && mode != 0 && outer * inner * n_out > 0)
        throw std::invalid_argument("an empty line can only be padded with a constant");
    const i64 nprog = (i64)prog.shape(0), plen = (i64)prog.shape(1);
    const i64 m = (i64)w.size();
    const int32_t *pp = prog.data();
    for (i64 k = 0; k < nprog * plen; k++) {
        const i64 op = pp[k];
        if (op == -1 || op == -2 || op == -3) continue;
        const i64 t = op >= 0 ? op : op <= -kEdgeHigh ? -op - kEdgeHigh : op <= -kEdgeLow ? -op - kEdgeLow : -op - 10;
        if (t < 0 || t >= m) throw std::invalid_argument("a tap out of range");
        if (t >= span && op >= 0 && mode != 0)
            throw std::invalid_argument("taps past the kernel need constant padding");
        if (op <= -kEdgeLow && mode != 0 && mode != 1)
            throw std::invalid_argument("edge taps need replicate or constant padding");
    }
    for (i64 i = 0; i < n_out; i++)
        if (rowmap.data()[i] < 0 || rowmap.data()[i] >= nprog)
            throw std::invalid_argument("rowmap names a program that does not exist");
    {
        nb::gil_scoped_release release;
        PK5D_FILT_CALL(isa, line_program<T>, src.data(), dst.data(), outer, n_in, inner, n_out,
                       w.data(), pp, plen, rowmap.data(), pre, span, mode, (T)value, num_threads);
    }
}

template <typename T>
void bind_taps_nd(pk5d::cvolume<const T> P, pk5d::cvolume<T> out, i64_array offs, in_array<T> w,
                  bool fused, int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    const i64 K = (i64)w.size();
    if ((i64)offs.size() != 3 * K) throw std::invalid_argument("offs must be K x 3");
    const i64 pz = (i64)P.shape(0), py = (i64)P.shape(1), px = (i64)P.shape(2);
    const i64 oz = (i64)out.shape(0), oy = (i64)out.shape(1), ox = (i64)out.shape(2);
    const int64_t *o = offs.data();
    for (i64 k = 0; k < K; k++) {
        if (o[3 * k] < 0 || o[3 * k + 1] < 0 || o[3 * k + 2] < 0 || o[3 * k] + oz > pz ||
            o[3 * k + 1] + oy > py || o[3 * k + 2] + ox > px)
            throw std::invalid_argument("a tap reads outside the padded volume");
    }
    {
        nb::gil_scoped_release release;
        PK5D_FILT_CALL(isa, taps_nd<T>, P.data(), py, px, out.data(), oz, oy, ox, o, w.data(), K,
                       fused, num_threads);
    }
}

}  // namespace

NB_MODULE(PK5D_MODULE_NAME, m) {
    m.doc() = "MATLAB's linear filters: tap programs along an axis and N-D tap correlations";
    m.def("best_isa", &best_isa, "the instruction set the kernels run (0 baseline, 1 AVX2 with FMA)");
    m.def("line_program_f32", &bind_line_program<float>, nb::arg("src").noconvert(),
          nb::arg("dst").noconvert(), nb::arg("outer"), nb::arg("n_in"), nb::arg("inner"), nb::arg("n_out"),
          nb::arg("w").noconvert(), nb::arg("prog").noconvert(), nb::arg("rowmap").noconvert(), nb::arg("pre"),
          nb::arg("span"), nb::arg("mode"), nb::arg("value"), nb::arg("num_threads"), nb::arg("isa") = -1);
    m.def("line_program_f64", &bind_line_program<double>, nb::arg("src").noconvert(),
          nb::arg("dst").noconvert(), nb::arg("outer"), nb::arg("n_in"), nb::arg("inner"), nb::arg("n_out"),
          nb::arg("w").noconvert(), nb::arg("prog").noconvert(), nb::arg("rowmap").noconvert(), nb::arg("pre"),
          nb::arg("span"), nb::arg("mode"), nb::arg("value"), nb::arg("num_threads"), nb::arg("isa") = -1);
    m.def("taps_nd_f32", &bind_taps_nd<float>, nb::arg("P").noconvert(), nb::arg("out").noconvert(),
          nb::arg("offs").noconvert(), nb::arg("w").noconvert(), nb::arg("fused"), nb::arg("num_threads"),
          nb::arg("isa") = -1);
    m.def("taps_nd_f64", &bind_taps_nd<double>, nb::arg("P").noconvert(), nb::arg("out").noconvert(),
          nb::arg("offs").noconvert(), nb::arg("w").noconvert(), nb::arg("fused"), nb::arg("num_threads"),
          nb::arg("isa") = -1);
}
