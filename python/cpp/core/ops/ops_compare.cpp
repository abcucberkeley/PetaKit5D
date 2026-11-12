// Comparison with a double scalar (K08).
//
// ops_compare(x, out, op, value, num_threads)
//   out = x == value (op 0), ~= (1), < (2), <= (3), > (4), >= (5): every value of x compared
//   exactly with value (NaN: only ~= is true). out is logical.
#include "ops_common.hpp"

namespace pk5d::ops {
namespace {

enum cmp_op { eq = 0, ne, lt, le, gt, ge };

template <typename T, typename P>
void predicate_loop(const T *x, bool *y, idx_t n, P p, int nt, bool bind) {
    parallel(nt, bind, [&](int k, int nk) {
        idx_t lo, hi;
        share(n, k, nk, lo, hi);
        for (idx_t i = lo; i < hi; i++) y[i] = p(x[i]);
    });
}

// x op k, in the class of x
template <typename T, typename K>
void compare_loop(const T *x, bool *y, idx_t n, int op, K k, int nt, bool bind) {
    switch (op) {
        case eq: return predicate_loop(x, y, n, [k](T a) { return static_cast<K>(a) == k; }, nt, bind);
        case ne: return predicate_loop(x, y, n, [k](T a) { return static_cast<K>(a) != k; }, nt, bind);
        case lt: return predicate_loop(x, y, n, [k](T a) { return static_cast<K>(a) < k; }, nt, bind);
        case le: return predicate_loop(x, y, n, [k](T a) { return static_cast<K>(a) <= k; }, nt, bind);
        case gt: return predicate_loop(x, y, n, [k](T a) { return static_cast<K>(a) > k; }, nt, bind);
        case ge: return predicate_loop(x, y, n, [k](T a) { return static_cast<K>(a) >= k; }, nt, bind);
    }
}

// float, double and up to 32-bit integers: compared in double, which is exact
template <typename T>
void compare_double(const T *x, bool *y, idx_t n, int op, double v, int nt, bool bind) {
    compare_loop<T, double>(x, y, n, op, v, nt, bind);
}

void fill(bool *y, idx_t n, bool v, int nt, bool bind) {
    parallel(nt, bind, [&](int k, int nk) {
        idx_t lo, hi;
        share(n, k, nk, lo, hi);
        for (idx_t i = lo; i < hi; i++) y[i] = v;
    });
}

// 64-bit integers against a double: the integer comparison with floor(v) or ceil(v), or a
// constant where those lie outside the class
template <typename T>
void compare_int64(const T *x, bool *y, idx_t n, int op, double v, int nt, bool bind) {
    if (v != v) return fill(y, n, op == ne, nt, bind);
    constexpr double lo = std::is_signed_v<T> ? -9223372036854775808.0 : 0.0;
    constexpr double top = std::is_signed_v<T> ? 9223372036854775808.0 : 18446744073709551616.0;
    const double f = std::floor(v), c = std::ceil(v);
    // where a bound lies: -1 below the class, 1 above it, 0 inside
    auto where = [&](double b) { return b < lo ? -1 : (b >= top ? 1 : 0); };
    switch (op) {
        case eq:
        case ne: {
            if (f != c || where(f) != 0) return fill(y, n, op == ne, nt, bind);
            return compare_loop<T, T>(x, y, n, op, static_cast<T>(f), nt, bind);
        }
        case lt:  // x < ceil(v)
        case ge:  // x >= ceil(v)
            if (where(c) != 0) return fill(y, n, (where(c) > 0) == (op == lt), nt, bind);
            return compare_loop<T, T>(x, y, n, op, static_cast<T>(c), nt, bind);
        case le:  // x <= floor(v)
        case gt:  // x > floor(v)
            if (where(f) != 0) return fill(y, n, (where(f) > 0) == (op == le), nt, bind);
            return compare_loop<T, T>(x, y, n, op, static_cast<T>(f), nt, bind);
    }
}

void compare(const flat_in &x, flat_out &out, int op, double value, int num_threads,
             bool bind) {
    if (x.shape(0) != out.shape(0))
        throw std::invalid_argument("ops_compare: input and output sizes differ");
    if (class_of(out.dtype()) != mclass::logical)
        throw std::invalid_argument("ops_compare: out must be logical");
    if (op < eq || op > ge) throw std::invalid_argument("ops_compare: op must be 0 to 5");
    const idx_t n = static_cast<idx_t>(x.shape(0));
    const int nt = team(n, num_threads);
    const void *px = x.data();
    bool *py = static_cast<bool *>(out.data());
    nb::gil_scoped_release release;
    with_class(class_of(x.dtype()), "ops_compare input", [&](auto ts) {
        using T = typename decltype(ts)::type;
        const T *xs = static_cast<const T *>(px);
        if constexpr (std::is_same_v<T, int64_t> || std::is_same_v<T, uint64_t>) {
            compare_int64<T>(xs, py, n, op, value, nt, bind);
        } else {
            compare_double<T>(xs, py, n, op, value, nt, bind);
        }
    });
}

}  // namespace
}  // namespace pk5d::ops

PK5D_BINDINGS(bind_ops_compare) {
    namespace nb = nanobind;
    m.def("ops_compare", &pk5d::ops::compare, nb::arg("x").noconvert(),
          nb::arg("out").noconvert(), nb::arg("op"), nb::arg("value"), nb::arg("num_threads"),
          nb::arg("bind") = true,
          "out = x op value (0 ==, 1 ~=, 2 <, 3 <=, 4 >, 5 >=), compared exactly; out logical.");
}
