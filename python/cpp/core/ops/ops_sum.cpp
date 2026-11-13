// MATLAB's vector rule of sums in its order at one thread, at any thread count (K08; the order
// of decon_sum_double and decon_lucy_accel_mex of the MATLAB code, and of MATLAB's sum, mean and
// var of single and double data, N05's docs/decisions/0005-sums-and-math-library.md).
//
// ops_sum_double(x, num_threads): x 1-D contiguous, sum(x(:), 'double') in its order (MATLAB's
//   linear order of a MATLAB-layout array, or a vector), with MATLAB's NaN payloads.
// ops_sum_double_yx(x, num_threads): the same of x (planes, ny, nx) C-contiguous, the (z, y, x)
//   volume of the port, in MATLAB's linear order t = y + ny * (x + nx * plane) (y fastest):
//   tiles of columns transposed into a scratch per thread (yx_block_sums), strided reads below
//   89000 elements. It costs a transposing read: 1.6 to 1.8 times the contiguous sum at 1 and 4
//   threads; a kernel that holds MATLAB's layout sums it contiguously.
// ops_sum_ordered(x, single, has_center, center, num_threads) and ops_sum_ordered_yx(...): the
//   same order with the additions in single (single=True, single data) or double, of the values
//   or, with has_center, of abs(x - center).^2 in the class of x (var's sum, center its mean in
//   that class). The payload of a NaN result is left to the caller (N05 canonicalizes it).
//
// The order, by the number n of elements, in MATLAB's linear index t:
//   n <= 43: in index order, ((x1 + x2) + x3) + ...;
//   44 <= n < 89000: four lanes (element t in lane mod(t, 4), each lane in increasing t), added
//     as ((l1 + l2) + l3) + l4;
//   n >= 89000: blocks of 4096 elements (block t / 4096; the last one shorter), each in four
//     lanes as above, the block sums added in index order. Each block is summed by one thread,
//     so the result doesn't depend on the number of threads (MATLAB splits the blocks among its
//     threads, which changes its last bits).
// Accumulators start at -0 (a sum of -0s is -0); an empty sum is 0. Logical, single, double
// and 8 to 32-bit integer data (converted exactly to double; single sums of single data only).
//
// ops_sum_double's NaN sum has the payload MATLAB gives it, which follows the operand order of
// each addition (x86 returns the first NaN operand, quieted). Measured in R2024b: the running
// sums are the first operand, except in the lanes' tails (the last mod(m, 4) elements of the
// four-lane sum, and of the last block of single data) and for the block sums of double data,
// which are the first operand of their additions. A sum that comes out NaN is summed again in
// that order (blocks_exact), since the fast loops leave the operand order to the compiler.
#include "ops_sum.hpp"

#include <cstdint>
#include <cstring>

namespace pk5d::ops {
namespace {

using flat3_in = nb::ndarray<nb::ro, nb::ndim<3>, nb::c_contig, nb::device::cpu>;

// a + b with the NaN of a, else of b, quieted (as x86 adds a NaN a and b); Inf + -Inf is x86's
// default NaN (sign bit set), also on CPUs whose default NaN is another (ARM's)
inline double add_first(double a, double b) {
    if (a != a) return a + a;
    if (b != b) return b + b;
    const double r = a + b;
    if (r == r) return r;
    const uint64_t u = 0xFFF8000000000000ull;
    double nan;
    std::memcpy(&nan, &u, sizeof nan);
    return nan;
}

// the four lanes of m values with MATLAB's operand order; tail_first: the last mod(m, 4)
// values are the first operand of their additions
double lanes_sum_exact(const double *p, idx_t m, bool tail_first) {
    double a[4] = {-0.0, -0.0, -0.0, -0.0};
    const idx_t main = m - m % 4;
    for (idx_t t = 0; t < main; t++) a[t & 3] = add_first(a[t & 3], p[t]);
    for (idx_t t = main; t < m; t++)
        a[t & 3] = tail_first ? add_first(p[t], a[t & 3]) : add_first(a[t & 3], p[t]);
    return add_first(add_first(add_first(a[0], a[1]), a[2]), a[3]);
}

double small_sum_exact(const double *v, idx_t n) {
    if (n <= kSequentialUpTo) {
        double s = -0.0;
        for (idx_t i = 0; i < n; i++) s = add_first(s, v[i]);
        return s;
    }
    return lanes_sum_exact(v, n, true);
}

// the blocks again on one thread with MATLAB's operand order (the payload of a NaN result)
template <typename Src, typename Op>
double blocks_exact(const Src &src, idx_t n, bool single, const Op &op) {
    const idx_t nb = (n + kBlock - 1) / kBlock;
    std::vector<double> buf(static_cast<size_t>(kBlock));
    double s = -0.0;
    for (idx_t b = 0; b < nb; b++) {
        const idx_t i0 = b * kBlock;
        const idx_t m = b < nb - 1 ? kBlock : n - i0;
        src.fill(i0, m, buf.data(), op);
        const double bs = lanes_sum_exact(buf.data(), m, single);
        s = single ? add_first(s, bs) : add_first(bs, s);
    }
    return s;
}

// sum(x(:), 'double') with MATLAB's NaN payloads
template <typename T, typename Src>
double sum_double_exact(const Src &src, idx_t n, int nt, bool bind) {
    constexpr bool single = std::is_same_v<T, float>;
    const op_value<T, double> op{};
    if (n > 0 && n < kBlocksFrom) {
        std::vector<double> v(static_cast<size_t>(n));
        src.fill(0, n, v.data(), op);
        const double s = small_sum(v.data(), n);
        if (is_float_v<T> && s != s) return small_sum_exact(v.data(), n);
        return s;
    }
    const double s = fixed_sum(src, n, nt, bind, op);
    if (is_float_v<T> && s != s) return blocks_exact(src, n, single, op);
    return s;
}

template <typename Fn>
void with_sum_class(nb::dlpack::dtype dt, Fn &&fn) {
    with_class(class_of(dt), "ops_sum input", [&](auto ts) {
        using T = typename decltype(ts)::type;
        if constexpr (std::is_same_v<T, int64_t> || std::is_same_v<T, uint64_t>) {
            bad_class("ops_sum input (64-bit integers)");
        } else {
            fn(ts);
        }
    });
}

double sum_double(const flat_in &x, int num_threads, bool bind) {
    check_threads(num_threads);
    const idx_t n = static_cast<idx_t>(x.shape(0));
    const void *px = x.data();
    double s = 0.0;
    nb::gil_scoped_release release;
    with_sum_class(x.dtype(), [&](auto ts) {
        using T = typename decltype(ts)::type;
        s = sum_double_exact<T>(linear_source<T>{static_cast<const T *>(px)}, n, num_threads,
                                bind);
    });
    return s;
}

double sum_double_yx(const flat3_in &x, int num_threads, bool bind) {
    check_threads(num_threads);
    const idx_t np_ = static_cast<idx_t>(x.shape(0));
    const idx_t ny = static_cast<idx_t>(x.shape(1));
    const idx_t nx = static_cast<idx_t>(x.shape(2));
    const idx_t n = np_ * ny * nx;
    const void *px = x.data();
    double s = 0.0;
    nb::gil_scoped_release release;
    with_sum_class(x.dtype(), [&](auto ts) {
        using T = typename decltype(ts)::type;
        const T *p = static_cast<const T *>(px);
        if (ny == 1 || nx == 1)
            s = sum_double_exact<T>(linear_source<T>{p}, n, num_threads, bind);
        else
            s = sum_double_exact<T>(yx_source<T>{p, ny, nx}, n, num_threads, bind);
    });
    return s;
}

// the sum in single or double of the values or of abs(x - center).^2, from a source
template <typename Fn>
void with_ordered(nb::dlpack::dtype dt, bool single, bool has_center, Fn &&fn) {
    with_sum_class(dt, [&](auto ts) {
        using T = typename decltype(ts)::type;
        if (has_center) {
            if constexpr (is_float_v<T>) {
                if (single != std::is_same_v<T, float>)
                    throw std::invalid_argument(
                        "ops_sum_ordered: abs(x - center).^2 is summed in the class of x");
                fn(op_sqdev<T>{});
            } else {
                throw std::invalid_argument("ops_sum_ordered: a center needs single or double data");
            }
        } else if (single) {
            if constexpr (std::is_same_v<T, float>)
                fn(op_value<float, float>{});
            else
                throw std::invalid_argument("ops_sum_ordered: single sums need single data");
        } else {
            fn(op_value<T, double>{});
        }
    });
}

double sum_ordered(const flat_in &x, bool single, bool has_center, double center, int num_threads,
                   bool bind) {
    check_threads(num_threads);
    const idx_t n = static_cast<idx_t>(x.shape(0));
    const void *px = x.data();
    double s = 0.0;
    nb::gil_scoped_release release;
    with_ordered(x.dtype(), single, has_center, [&](auto op0) {
        using Op = decltype(op0);
        using T = typename Op::src_t;
        Op op = op0;
        if constexpr (Op::sq) op.q = static_cast<T>(center);
        s = static_cast<double>(ordered_sum(linear_source<T>{static_cast<const T *>(px)}, n,
                                            num_threads, bind, op));
    });
    return s;
}

double sum_ordered_yx(const flat3_in &x, bool single, bool has_center, double center,
                      int num_threads, bool bind) {
    check_threads(num_threads);
    const idx_t np_ = static_cast<idx_t>(x.shape(0));
    const idx_t ny = static_cast<idx_t>(x.shape(1));
    const idx_t nx = static_cast<idx_t>(x.shape(2));
    const idx_t n = np_ * ny * nx;
    const void *px = x.data();
    double s = 0.0;
    nb::gil_scoped_release release;
    with_ordered(x.dtype(), single, has_center, [&](auto op0) {
        using Op = decltype(op0);
        using T = typename Op::src_t;
        Op op = op0;
        if constexpr (Op::sq) op.q = static_cast<T>(center);
        const T *p = static_cast<const T *>(px);
        if (ny == 1 || nx == 1)
            s = static_cast<double>(ordered_sum(linear_source<T>{p}, n, num_threads, bind, op));
        else
            s = static_cast<double>(
                ordered_sum(yx_source<T>{p, ny, nx}, n, num_threads, bind, op));
    });
    return s;
}

}  // namespace
}  // namespace pk5d::ops

PK5D_BINDINGS(bind_ops_sum) {
    namespace nb = nanobind;
    m.def("ops_sum_double", &pk5d::ops::sum_double, nb::arg("x").noconvert(),
          nb::arg("num_threads"), nb::arg("bind") = true,
          "sum(x(:), 'double') of a vector in its order, MATLAB's at one thread, independent of "
          "num_threads.");
    m.def("ops_sum_double_yx", &pk5d::ops::sum_double_yx, nb::arg("x").noconvert(),
          nb::arg("num_threads"), nb::arg("bind") = true,
          "sum(v(:), 'double') of a (planes, ny, nx) C array in MATLAB's linear order (y "
          "fastest), MATLAB's order at one thread, independent of num_threads.");
    m.def("ops_sum_ordered", &pk5d::ops::sum_ordered, nb::arg("x").noconvert(), nb::arg("single"),
          nb::arg("has_center"), nb::arg("center"), nb::arg("num_threads"),
          nb::arg("bind") = true,
          "MATLAB's vector-rule sum of a vector (or of abs(x - center).^2) in single or double.");
    m.def("ops_sum_ordered_yx", &pk5d::ops::sum_ordered_yx, nb::arg("x").noconvert(),
          nb::arg("single"), nb::arg("has_center"), nb::arg("center"), nb::arg("num_threads"),
          nb::arg("bind") = true,
          "MATLAB's vector-rule sum of a (planes, ny, nx) C array in MATLAB's linear order (or of "
          "abs(x - center).^2), in single or double.");
}
