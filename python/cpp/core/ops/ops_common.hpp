// Shared parts of the element-wise primitives (K08): array classes, dispatch, thread teams,
// and MATLAB's rounding and saturation as branch-free code the compiler vectorizes.
//
// Arrays arrive as 1-D contiguous views in memory order (petakit5d/kernels/ops.py makes them),
// so every primitive is a loop over n elements; the class of each array is checked here.
#pragma once

#include <pk5d/common.hpp>
#include <pk5d/fpexact.hpp>

#include "ops_bind.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace pk5d::ops {

namespace nb = nanobind;

using idx_t = std::ptrdiff_t;

// a 1-D contiguous array of any class: input (read-only allowed) and output
using flat_in = nb::ndarray<nb::ro, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using flat_out = nb::ndarray<nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// MATLAB's classes (char excluded)
enum class mclass { logical, i8, u8, i16, u16, i32, u32, i64, u64, single, dbl, other };

template <typename T>
struct tag {
    using type = T;
};

inline mclass class_of(nb::dlpack::dtype d) {
    if (d == nb::dtype<bool>()) return mclass::logical;
    if (d == nb::dtype<int8_t>()) return mclass::i8;
    if (d == nb::dtype<uint8_t>()) return mclass::u8;
    if (d == nb::dtype<int16_t>()) return mclass::i16;
    if (d == nb::dtype<uint16_t>()) return mclass::u16;
    if (d == nb::dtype<int32_t>()) return mclass::i32;
    if (d == nb::dtype<uint32_t>()) return mclass::u32;
    if (d == nb::dtype<int64_t>()) return mclass::i64;
    if (d == nb::dtype<uint64_t>()) return mclass::u64;
    if (d == nb::dtype<float>()) return mclass::single;
    if (d == nb::dtype<double>()) return mclass::dbl;
    return mclass::other;
}

[[noreturn]] inline void bad_class(const char *what) {
    throw std::invalid_argument(std::string(what) + ": unsupported array class");
}

// fn(tag<T>{}) for the C++ type of a class
template <typename Fn>
void with_class(mclass c, const char *what, Fn &&fn) {
    switch (c) {
        case mclass::logical: fn(tag<bool>{}); break;
        case mclass::i8: fn(tag<int8_t>{}); break;
        case mclass::u8: fn(tag<uint8_t>{}); break;
        case mclass::i16: fn(tag<int16_t>{}); break;
        case mclass::u16: fn(tag<uint16_t>{}); break;
        case mclass::i32: fn(tag<int32_t>{}); break;
        case mclass::u32: fn(tag<uint32_t>{}); break;
        case mclass::i64: fn(tag<int64_t>{}); break;
        case mclass::u64: fn(tag<uint64_t>{}); break;
        case mclass::single: fn(tag<float>{}); break;
        case mclass::dbl: fn(tag<double>{}); break;
        default: bad_class(what);
    }
}

template <typename T>
inline constexpr bool is_float_v = std::is_floating_point_v<T>;
template <typename T>
inline constexpr bool is_int_v = std::is_integral_v<T> && !std::is_same_v<T, bool>;
// the classes whose double-computed results are rounded and saturated (8 to 32 bits)
template <typename T>
inline constexpr bool is_small_int_v = is_int_v<T> && sizeof(T) <= 4;
// the targets of the floor(|x| + 0.5) rounding of single arrays
template <typename T>
inline constexpr bool is_fast_target_v = is_int_v<T> && sizeof(T) <= 2;

// Threads for n elements: at least kGrain elements per thread (element-wise results never
// depend on it)
constexpr idx_t kGrain = idx_t(1) << 15;

inline int team(idx_t n, int nt) {
    check_threads(nt);
    idx_t t = n / kGrain;
    if (t < 1) t = 1;
    return t < nt ? static_cast<int>(t) : nt;
}

// MATLAB's round: half away from zero, exactly, for any value (NaN and Inf pass). For |x|
// below 2^52 (2^23 in single) (a + big) - big rounds to an integer, corrected to floor(a).
template <typename F>
inline F round_half_away(F x) {
    constexpr F big = std::is_same_v<F, float> ? F(8388608.0f) : F(4503599627370496.0);
    const F a = std::fabs(x);
    F r = (a + big) - big;
    r = r > a ? r - F(1) : r;
    const F frac = a - r;
    r = frac >= F(0.5) ? r + F(1) : r;
    r = a < big ? r : a;
    return std::copysign(r, x);
}

// an integral double (or NaN, +-Inf) to an integer class with saturation, NaN to 0
template <typename D>
inline D saturate_integral(double r) {
    r = r == r ? r : 0.0;
    if constexpr (sizeof(D) <= 4) {
        constexpr double lo = static_cast<double>(std::numeric_limits<D>::min());
        constexpr double hi = static_cast<double>(std::numeric_limits<D>::max());
        r = r > lo ? r : lo;
        r = r < hi ? r : hi;
        if constexpr (std::is_same_v<D, uint32_t>) {
            return static_cast<D>(static_cast<int64_t>(r));
        } else {
            return static_cast<D>(static_cast<int32_t>(r));
        }
    } else if constexpr (std::is_signed_v<D>) {
        // [-2^63, 2^63)
        if (r >= 9223372036854775808.0) return std::numeric_limits<D>::max();
        if (r <= -9223372036854775808.0) return std::numeric_limits<D>::min();
        return static_cast<D>(r);
    } else {
        if (r >= 18446744073709551616.0) return std::numeric_limits<D>::max();
        if (r <= 0.0) return D(0);
        return static_cast<D>(r);
    }
}

// an integral float (or NaN, +-Inf) to an 8 or 16-bit class with saturation, NaN to 0
template <typename D>
inline D saturate_integral_f(float r) {
    constexpr float lo = static_cast<float>(std::numeric_limits<D>::min());
    constexpr float hi = static_cast<float>(std::numeric_limits<D>::max());
    r = r == r ? r : 0.0f;
    r = r > lo ? r : lo;
    r = r < hi ? r : hi;
    return static_cast<D>(static_cast<int32_t>(r));
}

// MATLAB's conversion of one value of class S to class D. Exact: round half away from zero
// (cast, and the class functions except below); otherwise (the class function on an array)
// single to 8 and 16-bit classes is sign(x) * floor(|x| + 0.5) in single, and 64-bit integers
// to single go through double (rounded twice).
template <typename S, typename D>
inline constexpr bool has_inexact_v =
    (std::is_same_v<S, float> && is_fast_target_v<D>) ||
    (is_int_v<S> && sizeof(S) == 8 && std::is_same_v<D, float>);

template <typename D, typename S, bool Exact>
inline D convert_value(S x) {
    if constexpr (std::is_same_v<S, bool>) {
        return static_cast<D>(x ? 1 : 0);
    } else if constexpr (is_float_v<D>) {
        if constexpr (!Exact && std::is_same_v<D, float> && is_int_v<S> && sizeof(S) == 8) {
            // single() of a 64-bit integer array: through double, as MATLAB (cast: directly)
            return static_cast<float>(static_cast<double>(x));
        } else if constexpr (std::is_same_v<D, float> && is_int_v<S> && sizeof(S) == 8) {
            return int64_to_float(x);  // rounded once on every processor (fpexact.hpp)
        } else {
            return static_cast<D>(x);
        }
    } else if constexpr (is_int_v<S>) {
        // integer to integer: saturate
        if constexpr (std::is_same_v<S, uint64_t>) {
            constexpr uint64_t hi = static_cast<uint64_t>(std::numeric_limits<D>::max());
            return static_cast<D>(x < hi ? x : hi);
        } else if constexpr (std::is_same_v<D, uint64_t>) {
            return static_cast<D>(x > 0 ? x : 0);
        } else {
            constexpr int64_t lo = static_cast<int64_t>(std::numeric_limits<D>::min());
            constexpr int64_t hi = static_cast<int64_t>(std::numeric_limits<D>::max());
            int64_t v = static_cast<int64_t>(x);
            v = v > lo ? v : lo;
            v = v < hi ? v : hi;
            return static_cast<D>(v);
        }
    } else if constexpr (std::is_same_v<S, float> && is_fast_target_v<D>) {
        if constexpr (Exact) {
            return saturate_integral_f<D>(round_half_away(x));
        } else {
            // floor(|x| + 0.5) in single: |x| + 0.5 >= 0.5 (NaN: 0), at most 2^16 kept
            float t = std::fabs(x) + 0.5f;
            t = t >= 0.0f ? t : 0.0f;
            t = t < 65536.0f ? t : 65536.0f;
            const int32_t i = static_cast<int32_t>(t);
            int32_t v = x < 0.0f ? -i : i;
            constexpr int32_t lo = std::numeric_limits<D>::min();
            constexpr int32_t hi = std::numeric_limits<D>::max();
            v = v > lo ? v : lo;
            v = v < hi ? v : hi;
            return static_cast<D>(v);
        }
    } else if constexpr (std::is_same_v<S, float>) {
        return saturate_integral<D>(static_cast<double>(round_half_away(x)));
    } else {
        // double to integer
        return saturate_integral<D>(round_half_away(x));
    }
}

// one step of MATLAB arithmetic of data of class W with a double scalar: in double, then
// rounded to single, or rounded half away and saturated for 8 to 32-bit integers
template <typename W>
inline W to_work(double r) {
    if constexpr (std::is_same_v<W, double>) {
        return r;
    } else if constexpr (std::is_same_v<W, float>) {
        return static_cast<float>(r);
    } else {
        return saturate_integral<W>(round_half_away(r));
    }
}

inline void check_same_size(const flat_in &x, const flat_out &y, const char *what) {
    if (x.shape(0) != y.shape(0))
        throw std::invalid_argument(std::string(what) + ": input and output sizes differ");
}

}  // namespace pk5d::ops
