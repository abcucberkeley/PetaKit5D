// Class conversion and the scale/offset/clamp/mask chain with MATLAB's rounding (K08).
//
// ops_convert(x, out, exact, num_threads)
//   out = cls(x) (exact false: the class function on an array) or cast(x, cls) (exact true),
//   cls the class of out; any class pair of logical and the numeric classes to numeric.
// ops_affine(x, out, has_scale, scale, has_offset, offset, has_hi, hi, has_lo, lo, mask,
//            exact, num_threads)
//   y = x * scale + offset; y = min(y, hi); y = max(y, lo); y = y .* mask; out = cls(y),
//   each step only when its flag is set (mask: None for no step), every step rounded to the
//   class of x as MATLAB rounds it. x is single, double or an 8 to 32-bit integer class; hi
//   and lo are already values of x's class (ops.py converts them); mask is logical and needs
//   single or double data.
#include "ops_common.hpp"

namespace pk5d::ops {
namespace {

template <typename S, typename D, bool Exact>
void convert_loop(const S *x, D *y, idx_t n, int nt, bool bind) {
    parallel(nt, bind, [&](int k, int nk) {
        idx_t lo, hi;
        share(n, k, nk, lo, hi);
        for (idx_t i = lo; i < hi; i++) y[i] = convert_value<D, S, Exact>(x[i]);
    });
}

struct affine_args {
    bool has_scale, has_offset, has_hi, has_lo;
    double scale, offset, hi, lo;
};

// v op scalar s with a NaN v: v itself, quieted (v + v), whatever operand order the compiler
// chose (x86 returns the first NaN operand, and so does MATLAB); Inf * 0 and Inf - Inf give
// MATLAB's NaN on every processor (fpexact.hpp)
template <typename S>
inline S nan_first(S v, S r, double s) {
    if constexpr (is_float_v<S>) {
        return v != v ? v + v : PK5D_X86_NAN(r, v, static_cast<S>(s));
    } else {
        return r;
    }
}

// the chain for one value; the flags are loop-invariant, so the selects vectorize
template <typename S>
inline S affine_value(S v, const affine_args &a, S hi, S lo) {
    const S s1 = nan_first(v, to_work<S>(static_cast<double>(v) * a.scale), a.scale);
    v = a.has_scale ? s1 : v;
    const S s2 = nan_first(v, to_work<S>(static_cast<double>(v) + a.offset), a.offset);
    v = a.has_offset ? s2 : v;
    // min(v, hi): v on a tie, hi for NaN v; max(v, lo) likewise
    const S s3 = PK5D_LE_SELECT(v, hi);
    v = a.has_hi ? s3 : v;
    const S s4 = PK5D_GE_SELECT(v, lo);
    v = a.has_lo ? s4 : v;
    return v;
}

template <typename S, typename D, bool Exact, bool Mask>
void affine_loop(const S *x, D *y, const bool *mask, idx_t n, const affine_args &a, int nt,
                 bool bind) {
    const S hi = static_cast<S>(a.hi);
    const S lo = static_cast<S>(a.lo);
    parallel(nt, bind, [&](int k, int nk) {
        idx_t i0, i1;
        share(n, k, nk, i0, i1);
        for (idx_t i = i0; i < i1; i++) {
            S v = affine_value<S>(x[i], a, hi, lo);
            if constexpr (Mask) v = v * static_cast<S>(mask[i] ? 1 : 0);
            y[i] = convert_value<D, S, Exact>(v);
        }
    });
}

void convert(const flat_in &x, flat_out &out, bool exact, int num_threads, bool bind) {
    check_same_size(x, out, "ops_convert");
    const idx_t n = static_cast<idx_t>(x.shape(0));
    const int nt = team(n, num_threads);
    const mclass cd = class_of(out.dtype());
    if (cd == mclass::logical) bad_class("ops_convert output");
    const void *px = x.data();
    void *py = out.data();
    nb::gil_scoped_release release;
    with_class(class_of(x.dtype()), "ops_convert input", [&](auto ts) {
        using S = typename decltype(ts)::type;
        with_class(cd, "ops_convert output", [&](auto td) {
            using D = typename decltype(td)::type;
            if constexpr (!std::is_same_v<D, bool>) {
                const S *xs = static_cast<const S *>(px);
                D *yd = static_cast<D *>(py);
                if constexpr (has_inexact_v<S, D>) {
                    if (exact)
                        convert_loop<S, D, true>(xs, yd, n, nt, bind);
                    else
                        convert_loop<S, D, false>(xs, yd, n, nt, bind);
                } else {
                    convert_loop<S, D, true>(xs, yd, n, nt, bind);
                }
            }
        });
    });
}

template <typename S, typename D, bool Exact>
void affine_masked(const S *xs, D *yd, const bool *mask, idx_t n, const affine_args &a,
                   int nt, bool bind) {
    if (mask != nullptr)
        affine_loop<S, D, Exact, true>(xs, yd, mask, n, a, nt, bind);
    else
        affine_loop<S, D, Exact, false>(xs, yd, nullptr, n, a, nt, bind);
}

void affine(const flat_in &x, flat_out &out, bool has_scale, double scale, bool has_offset,
            double offset, bool has_hi, double hi, bool has_lo, double lo, const flat_in &mask,
            bool exact, int num_threads, bool bind) {
    check_same_size(x, out, "ops_affine");
    const idx_t n = static_cast<idx_t>(x.shape(0));
    const int nt = team(n, num_threads);
    const mclass cs = class_of(x.dtype());
    const mclass cd = class_of(out.dtype());
    if (cd == mclass::logical) bad_class("ops_affine output");
    const bool *pm = nullptr;
    if (mask.is_valid()) {
        if (class_of(mask.dtype()) != mclass::logical)
            throw std::invalid_argument("ops_affine: mask must be logical");
        if (mask.shape(0) != x.shape(0))
            throw std::invalid_argument("ops_affine: mask and data sizes differ");
        if (cs != mclass::single && cs != mclass::dbl)
            throw std::invalid_argument("ops_affine: a mask needs single or double data");
        pm = static_cast<const bool *>(mask.data());
    }
    if (cs == mclass::logical || cs == mclass::i64 || cs == mclass::u64)
        bad_class("ops_affine input");
    const affine_args a{has_scale, has_offset, has_hi, has_lo, scale, offset, hi, lo};
    const void *px = x.data();
    void *py = out.data();
    nb::gil_scoped_release release;
    with_class(cs, "ops_affine input", [&](auto ts) {
        using S = typename decltype(ts)::type;
        if constexpr (is_float_v<S> || is_small_int_v<S>) {
            with_class(cd, "ops_affine output", [&](auto td) {
                using D = typename decltype(td)::type;
                if constexpr (!std::is_same_v<D, bool>) {
                    const S *xs = static_cast<const S *>(px);
                    D *yd = static_cast<D *>(py);
                    if constexpr (has_inexact_v<S, D>) {
                        if (exact)
                            affine_masked<S, D, true>(xs, yd, pm, n, a, nt, bind);
                        else
                            affine_masked<S, D, false>(xs, yd, pm, n, a, nt, bind);
                    } else {
                        affine_masked<S, D, true>(xs, yd, pm, n, a, nt, bind);
                    }
                }
            });
        }
    });
}

}  // namespace
}  // namespace pk5d::ops

PK5D_BINDINGS(bind_ops_convert) {
    namespace nb = nanobind;
    m.def("ops_convert", &pk5d::ops::convert, nb::arg("x").noconvert(),
          nb::arg("out").noconvert(), nb::arg("exact"), nb::arg("num_threads"),
          nb::arg("bind") = true,
          "out = cls(x) (exact=False) or cast(x, cls) (exact=True) with MATLAB's rounding and "
          "saturation; cls is the class of out (1-D contiguous arrays).");
    m.def("ops_affine", &pk5d::ops::affine, nb::arg("x").noconvert(), nb::arg("out").noconvert(),
          nb::arg("has_scale"), nb::arg("scale"), nb::arg("has_offset"), nb::arg("offset"),
          nb::arg("has_hi"), nb::arg("hi"), nb::arg("has_lo"), nb::arg("lo"),
          nb::arg("mask").noconvert().none(), nb::arg("exact"), nb::arg("num_threads"),
          nb::arg("bind") = true,
          "MATLAB's y = x * scale + offset; min(y, hi); max(y, lo); y .* mask; cls(y), each "
          "step when given, rounded to the class of x at every step.");
}
