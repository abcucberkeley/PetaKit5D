// petakit5d._ext.resize: the arithmetic of MATLAB's imresize and imresize3 on C-contiguous arrays
// (petakit5d.matlabcompat.resize computes the weights; measured on MATLAB R2024b, black box):
// - resize_lines: a weighted sum along one axis of an array viewed as (outer, n, inner), in
//   double from +0 in tap order, rounded to the class as MATLAB's resize rounds (a NaN sum keeps
//   the first NaN in tap order, MATLAB's payload rule);
// - gather_lines: index copies along an axis (MATLAB copies where the weights are one-hot);
// - cubic: MATLAB's cubic kernel, |x|^3 from fdlibm's pow (MATLAB's power).
// Every output value is computed by one thread in a fixed order, so the results do not depend
// on the thread count. The NaN of an invalid operation in pow and the cubic kernel (0 / 0,
// Inf - Inf, Inf * 0) is x86's (MATLAB's) on every processor (PK5D_X86_NAN, the plain operation
// except on ARM).
#include <pk5d/common.hpp>
#include <pk5d/fpexact.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace nb = nanobind;
using i64 = int64_t;

#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
#define PK5D_RS_DISPATCH 1
#endif

namespace {

// MATLAB's NaN (the sign bit set) and a NaN quieted with its payload kept, from the bits (the
// same on every processor)
inline double matlab_nan() {
    const uint64_t b = 0xFFF8000000000000ull;
    double d;
    std::memcpy(&d, &b, 8);
    return d;
}

inline double quiet_nan(double x) {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    b |= 0x0008000000000000ull;
    double d;
    std::memcpy(&d, &b, 8);
    return d;
}

}  // namespace

#define PK5D_RS_NS isa_baseline
#define PK5D_RS_ISA
#include "kernels.inc"
#undef PK5D_RS_NS
#undef PK5D_RS_ISA
#ifdef PK5D_RS_DISPATCH
#define PK5D_RS_NS isa_avx2
#define PK5D_RS_ISA __attribute__((target("avx2")))
#include "kernels.inc"
#undef PK5D_RS_NS
#undef PK5D_RS_ISA
#endif

namespace {

// instruction sets: 0 baseline, 1 AVX2
int best_isa() {
#ifdef PK5D_RS_DISPATCH
    static const int best = [] {
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2") ? 1 : 0;
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

// ---------------------------------------------------------------------------------------------
// fdlibm 5.3's __ieee754_pow (e_pow.c; Sun Microsystems' permissive notice), MATLAB's x .^ y,
// as matlabcompat.elementary writes it out in Python: one IEEE operation per C operation
// ---------------------------------------------------------------------------------------------

inline int32_t hi_word(double x) {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return static_cast<int32_t>(static_cast<uint32_t>(b >> 32));
}
inline uint32_t lo_word(double x) {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return static_cast<uint32_t>(b);
}
inline double make_double(int32_t hi, uint32_t lo) {
    const uint64_t b = (static_cast<uint64_t>(static_cast<uint32_t>(hi)) << 32) | lo;
    double d;
    std::memcpy(&d, &b, 8);
    return d;
}
inline double with_lo0(double v) { return make_double(hi_word(v), 0); }

double fdlibm_pow(double x, double y) {
    static const double bp[2] = {1.0, 1.5};
    static const double dp_h[2] = {0.0, 5.84962487220764160156e-01};
    static const double dp_l[2] = {0.0, 1.35003920212974897128e-08};
    const double two53 = 9007199254740992.0, huge = 1.0e300, tiny = 1.0e-300;
    const double L1 = 5.99999999999994648725e-01, L2 = 4.28571428578550184252e-01,
                 L3 = 3.33333329818377432918e-01, L4 = 2.72728123808534006489e-01,
                 L5 = 2.30660745775561754067e-01, L6 = 2.06975017800338417784e-01;
    const double P1 = 1.66666666666666019037e-01, P2 = -2.77777777770155933842e-03,
                 P3 = 6.61375632143793436117e-05, P4 = -1.65339022054652515390e-06,
                 P5 = 4.13813679705723846039e-08;
    const double lg2 = 6.93147180559945286227e-01, lg2_h = 6.93147182464599609375e-01,
                 lg2_l = -1.90465429995776804525e-09, ovt = 8.0085662595372944372e-17;
    const double cp = 9.61796693925975554329e-01, cp_h = 9.61796700954437255859e-01,
                 cp_l = -7.02846165095275826516e-09;
    const double ivln2 = 1.44269504088896338700e+00, ivln2_h = 1.44269502162933349609e+00,
                 ivln2_l = 1.92596299112661746887e-08;

    double z, ax, z_h, z_l, p_h, p_l, y1, t1, t2, r, s, t, u, v, w;
    int32_t i, j, k, yisint, n;
    int32_t hx = hi_word(x), hy = hi_word(y);
    uint32_t lx = lo_word(x), ly = lo_word(y);
    int32_t ix = hx & 0x7fffffff, iy = hy & 0x7fffffff;

    if ((iy | static_cast<int32_t>(ly)) == 0) return 1.0;
    if (ix > 0x7ff00000 || (ix == 0x7ff00000 && lx != 0) || iy > 0x7ff00000 ||
        (iy == 0x7ff00000 && ly != 0))
        return x + y;
    yisint = 0;
    if (hx < 0) {
        if (iy >= 0x43400000) {
            yisint = 2;
        } else if (iy >= 0x3ff00000) {
            k = (iy >> 20) - 0x3ff;
            if (k > 20) {
                const uint32_t jj = ly >> (52 - k);
                if ((jj << (52 - k)) == ly) yisint = 2 - static_cast<int32_t>(jj & 1);
            } else if (ly == 0) {
                j = iy >> (20 - k);
                if ((j << (20 - k)) == iy) yisint = 2 - (j & 1);
            }
        }
    }
    if (ly == 0) {
        if (iy == 0x7ff00000) {
            if (((ix - 0x3ff00000) | static_cast<int32_t>(lx)) == 0) {
                const double r = y - y;  // (+-1)^(+-Inf): Inf - Inf
                return PK5D_X86_NAN(r, y, y);
            }
            if (ix >= 0x3ff00000) return hy >= 0 ? y : 0.0;
            return hy < 0 ? -y : 0.0;
        }
        if (iy == 0x3ff00000) return hy < 0 ? 1.0 / x : x;
        if (hy == 0x40000000) return x * x;
        if (hy == 0x3fe00000 && hx >= 0) return std::sqrt(x);
    }
    ax = std::fabs(x);
    if (lx == 0 && (ix == 0x7ff00000 || ix == 0 || ix == 0x3ff00000)) {
        z = ax;
        if (hy < 0) z = 1.0 / z;
        if (hx < 0) {
            if (((ix - 0x3ff00000) | yisint) == 0) {
                const double d = z - z, r = d / d;  // (-1)^(not an integer): 0 / 0
                z = PK5D_X86_NAN(r, d, d);
            } else if (yisint == 1) {
                z = -z;
            }
        }
        return z;
    }
    n = (hx >> 31) + 1;
    if ((n | yisint) == 0) {
        const double d = x - x, r = d / d;  // x < 0, y not an integer: 0 / 0
        return PK5D_X86_NAN(r, d, d);
    }
    s = 1.0;
    if ((n | (yisint - 1)) == 0) s = -1.0;
    if (iy > 0x41e00000) {
        if (iy > 0x43f00000) {
            if (ix <= 0x3fefffff) return hy < 0 ? huge * huge : tiny * tiny;
            if (ix >= 0x3ff00000) return hy > 0 ? huge * huge : tiny * tiny;
        }
        if (ix < 0x3fefffff) return hy < 0 ? s * huge * huge : s * tiny * tiny;
        if (ix > 0x3ff00000) return hy > 0 ? s * huge * huge : s * tiny * tiny;
        t = ax - 1.0;
        w = (t * t) * (0.5 - t * (0.3333333333333333333333 - t * 0.25));
        u = ivln2_h * t;
        v = t * ivln2_l - w * ivln2;
        t1 = with_lo0(u + v);
        t2 = v - (t1 - u);
    } else {
        double ss, s2, s_h, s_l, t_h, t_l;
        n = 0;
        if (ix < 0x00100000) {
            ax *= two53;
            n -= 53;
            ix = hi_word(ax);
        }
        n += ((ix) >> 20) - 0x3ff;
        j = ix & 0x000fffff;
        ix = j | 0x3ff00000;
        if (j <= 0x3988E) {
            k = 0;
        } else if (j < 0xBB67A) {
            k = 1;
        } else {
            k = 0;
            n += 1;
            ix -= 0x00100000;
        }
        ax = make_double(ix, lo_word(ax));
        u = ax - bp[k];
        v = 1.0 / (ax + bp[k]);
        ss = u * v;
        s_h = with_lo0(ss);
        t_h = make_double(((ix >> 1) | 0x20000000) + 0x00080000 + (k << 18), 0);
        t_l = ax - (t_h - bp[k]);
        s_l = v * ((u - s_h * t_h) - s_h * t_l);
        s2 = ss * ss;
        r = s2 * s2 * (L1 + s2 * (L2 + s2 * (L3 + s2 * (L4 + s2 * (L5 + s2 * L6)))));
        r += s_l * (s_h + ss);
        s2 = s_h * s_h;
        t_h = with_lo0(3.0 + s2 + r);
        t_l = r - ((t_h - 3.0) - s2);
        u = s_h * t_h;
        v = s_l * t_h + t_l * ss;
        p_h = with_lo0(u + v);
        p_l = v - (p_h - u);
        z_h = cp_h * p_h;
        z_l = cp_l * p_h + p_l * cp + dp_l[k];
        t = static_cast<double>(n);
        t1 = with_lo0(((z_h + z_l) + dp_h[k]) + t);
        t2 = z_l - (((t1 - t) - dp_h[k]) - z_h);
    }
    y1 = with_lo0(y);
    p_l = (y - y1) * t1 + y * t2;
    p_h = y1 * t1;
    z = p_l + p_h;
    j = hi_word(z);
    i = static_cast<int32_t>(lo_word(z));
    if (j >= 0x40900000) {
        if (((j - 0x40900000) | i) != 0) return s * huge * huge;
        if (p_l + ovt > z - p_h) return s * huge * huge;
    } else if ((j & 0x7fffffff) >= 0x4090cc00) {
        if (((j - static_cast<int32_t>(0xc090cc00)) | i) != 0) return s * tiny * tiny;
        if (p_l <= z - p_h) return s * tiny * tiny;
    }
    i = j & 0x7fffffff;
    k = (i >> 20) - 0x3ff;
    n = 0;
    if (i > 0x3fe00000) {
        n = j + (0x00100000 >> (k + 1));
        k = ((n & 0x7fffffff) >> 20) - 0x3ff;
        t = make_double(n & ~(0x000fffff >> k), 0);
        n = ((n & 0x000fffff) | 0x00100000) >> (20 - k);
        if (j < 0) n = -n;
        p_h -= t;
    }
    t = with_lo0(p_l + p_h);
    u = t * lg2_h;
    v = (p_l - (t - p_h)) * lg2 + t * lg2_l;
    z = u + v;
    w = v - (z - u);
    t = z * z;
    t1 = z - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    r = (z * t1) / (t1 - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    j = hi_word(z);
    j += static_cast<int32_t>(static_cast<uint32_t>(n) << 20);
    if ((j >> 20) <= 0)
        z = std::ldexp(z, n);
    else
        z = make_double(j, lo_word(z));
    return s * z;
}

#ifdef PK5D_FP_ARM
// a + b, a - b and a * b with x86's NaN for Inf - Inf and Inf * 0
inline double add_x86(double a, double b) {
    const double r = a + b;
    return PK5D_X86_NAN(r, a, b);
}
inline double sub_x86(double a, double b) {
    const double r = a - b;
    return PK5D_X86_NAN(r, a, b);
}
inline double mul_x86(double a, double b) {
    const double r = a * b;
    return PK5D_X86_NAN(r, a, b);
}
#endif

// MATLAB's cubic kernel (Keys, a = -1/2):
// (1.5|x|^3 - 2.5|x|^2 + 1) (|x| <= 1) + (-0.5|x|^3 + 2.5|x|^2 - 4|x| + 2) (1 < |x| <= 2)
// (NaN where |x|^3 overflows: Inf * 0, and Inf - Inf where |x|^2 does too; on ARM those
// operations written out, so that x86-64 compiles the code it compiled before)
double cubic_value(double x) {
    const double a = std::fabs(x);
    const double a2 = a * a;
    const double a3 = fdlibm_pow(a, 3.0);
#ifdef PK5D_FP_ARM
    const double p1 = mul_x86(sub_x86(1.5 * a3, 2.5 * a2) + 1, a <= 1 ? 1.0 : 0.0);
    const double p2 = mul_x86(add_x86(-0.5 * a3, 2.5 * a2) - 4 * a + 2,
                              (1 < a && a <= 2) ? 1.0 : 0.0);
#else
    const double p1 = (1.5 * a3 - 2.5 * a2 + 1) * (a <= 1 ? 1.0 : 0.0);
    const double p2 = (-0.5 * a3 + 2.5 * a2 - 4 * a + 2) * ((1 < a && a <= 2) ? 1.0 : 0.0);
#endif
    return p1 + p2;
}

using f64_in = nb::ndarray<const double, nb::c_contig, nb::device::cpu>;
using f64_out = nb::ndarray<double, nb::c_contig, nb::device::cpu>;
using i64_in = nb::ndarray<const int64_t, nb::c_contig, nb::device::cpu>;
template <typename T>
using in_array = nb::ndarray<const T, nb::c_contig, nb::device::cpu>;
template <typename T>
using out_array = nb::ndarray<T, nb::c_contig, nb::device::cpu>;

void bind_cubic(f64_in x, f64_out out) {
    if (x.size() != out.size()) throw std::invalid_argument("x and out differ in size");
    const double *px = x.data();
    double *po = out.data();
    const size_t n = x.size();
    nb::gil_scoped_release release;
    for (size_t i = 0; i < n; i++) po[i] = cubic_value(px[i]);
}

void bind_pow(f64_in x, double y, f64_out out) {
    if (x.size() != out.size()) throw std::invalid_argument("x and out differ in size");
    const double *px = x.data();
    double *po = out.data();
    const size_t n = x.size();
    nb::gil_scoped_release release;
    for (size_t i = 0; i < n; i++) po[i] = fdlibm_pow(px[i], y);
}

void check_shape(i64 outer, i64 n_in, i64 inner, i64 n_out, size_t src, size_t dst) {
    if (outer < 0 || n_in < 0 || inner < 0 || n_out < 0) throw std::invalid_argument("negative sizes");
    if (static_cast<i64>(src) != outer * n_in * inner) throw std::invalid_argument("src has the wrong size");
    if (static_cast<i64>(dst) != outer * n_out * inner) throw std::invalid_argument("dst has the wrong size");
}

template <typename T>
void bind_resize(in_array<T> src, out_array<T> dst, i64 outer, i64 n_in, i64 inner, i64 n_out,
                 f64_in w, i64_in ind, int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    check_shape(outer, n_in, inner, n_out, src.size(), dst.size());
    if (w.ndim() != 2 || ind.ndim() != 2 || static_cast<i64>(w.shape(0)) != n_out ||
        w.shape(0) != ind.shape(0) || w.shape(1) != ind.shape(1))
        throw std::invalid_argument("weights and indices must both be n_out x P");
    const i64 P = static_cast<i64>(w.shape(1));
    const int64_t *pi = ind.data();
    for (i64 k = 0; k < n_out * P; k++)
        if (pi[k] < 0 || pi[k] >= n_in) throw std::invalid_argument("an index out of range");
    if (static_cast<const void *>(src.data()) == static_cast<const void *>(dst.data()) && dst.size())
        throw std::invalid_argument("src and dst must not overlap");
    nb::gil_scoped_release release;
#ifdef PK5D_RS_DISPATCH
    if (isa == 1) {
        isa_avx2::resize_lines<T>(src.data(), dst.data(), outer, n_in, inner, n_out, w.data(), pi, P,
                                  num_threads);
        return;
    }
#endif
    isa_baseline::resize_lines<T>(src.data(), dst.data(), outer, n_in, inner, n_out, w.data(), pi, P,
                                  num_threads);
}

template <typename T>
void bind_gather(in_array<T> src, out_array<T> dst, i64 outer, i64 n_in, i64 inner, i64 n_out,
                 i64_in idx, int num_threads) {
    pk5d::check_threads(num_threads);
    check_shape(outer, n_in, inner, n_out, src.size(), dst.size());
    if (static_cast<i64>(idx.size()) != n_out) throw std::invalid_argument("one index per output");
    const int64_t *pi = idx.data();
    for (i64 k = 0; k < n_out; k++)
        if (pi[k] < 0 || pi[k] >= n_in) throw std::invalid_argument("an index out of range");
    nb::gil_scoped_release release;
    isa_baseline::gather_lines<T>(src.data(), dst.data(), outer, n_in, inner, n_out, pi, num_threads);
}

template <typename T>
void def_resize(nb::module_ &m, const char *name) {
    m.def(name, &bind_resize<T>, nb::arg("src").noconvert(), nb::arg("dst").noconvert(), nb::arg("outer"),
          nb::arg("n_in"), nb::arg("inner"), nb::arg("n_out"), nb::arg("w").noconvert(),
          nb::arg("ind").noconvert(), nb::arg("num_threads"), nb::arg("isa") = -1);
}

template <typename T>
void def_gather(nb::module_ &m, const char *name) {
    m.def(name, &bind_gather<T>, nb::arg("src").noconvert(), nb::arg("dst").noconvert(), nb::arg("outer"),
          nb::arg("n_in"), nb::arg("inner"), nb::arg("n_out"), nb::arg("idx").noconvert(),
          nb::arg("num_threads"));
}

}  // namespace

NB_MODULE(PK5D_MODULE_NAME, m) {
    m.doc() = "MATLAB's imresize and imresize3 passes: weighted sums along an axis, index copies, cubic";
    m.def("best_isa", &best_isa, "the instruction set the kernels run (0 baseline, 1 AVX2)");
    m.def("cubic", &bind_cubic, nb::arg("x").noconvert(), nb::arg("out").noconvert(),
          "MATLAB's cubic kernel of x (double) into out");
    m.def("fdlibm_pow", &bind_pow, nb::arg("x").noconvert(), nb::arg("y"), nb::arg("out").noconvert(),
          "fdlibm's pow(x, y) of each x (double) into out");
    def_resize<uint8_t>(m, "resize_uint8");
    def_resize<int8_t>(m, "resize_int8");
    def_resize<uint16_t>(m, "resize_uint16");
    def_resize<int16_t>(m, "resize_int16");
    def_resize<uint32_t>(m, "resize_uint32");
    def_resize<int32_t>(m, "resize_int32");
    def_resize<float>(m, "resize_single");
    def_resize<double>(m, "resize_double");
    def_gather<uint8_t>(m, "gather_1");
    def_gather<uint16_t>(m, "gather_2");
    def_gather<uint32_t>(m, "gather_4");
    def_gather<uint64_t>(m, "gather_8");
}
