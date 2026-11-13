// Floating-point steps whose bits differ by compiler or processor unless written out
// (docs/platforms.md). On x86-64 with GCC or Clang the macros expand to the plain expressions,
// so those builds compile the code they compiled before.
//
// - PK5D_GE_SELECT(a, b), PK5D_LE_SELECT(a, b): a >= b ? a : b and a <= b ? a : b as written
//   (a on a tie, also -0 against +0; b when either is NaN). MSVC (19.5x, also with /fp:strict)
//   compiles the plain selects to maxss/minss, which give b on a tie, so max(-0, 0) became +0;
//   there the select is made on the bits.
// - PK5D_X86_NAN(r, a, b): r, the result of one IEEE operation on a and b, with x86's default
//   NaN (sign bit set: MATLAB's NaN) where r is NaN and neither operand is (Inf - Inf, 0 * Inf,
//   0 / 0). ARM's default NaN has the sign bit clear; elsewhere r itself.
// - pk5d::int64_to_float(x): x rounded once to the nearest float (ties to even), as x86's
//   cvtsi2ss. AArch64 compilers convert vectors of 64-bit integers through double, which rounds
//   twice for |x| > 2^53; there the value is first rounded to odd at 53 bits, which makes the
//   double exact and its rounding to float the correct one.
// The arguments of the macros are evaluated more than once: pass plain values.
#pragma once

#include <cstdint>
#include <cstring>
#include <type_traits>

#if defined(_MSC_VER) && !defined(__clang__)
#define PK5D_FP_MSVC_SELECT 1
#endif
#if defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__) || defined(_M_ARM)
#define PK5D_FP_ARM 1
#endif

namespace pk5d {

namespace fp_detail {
template <typename T>
using bits_t = std::conditional_t<sizeof(T) == 4, uint32_t, uint64_t>;

// c ? a : b on the bits (MSVC keeps it a select)
template <typename T>
inline T select_bits(bool c, T a, T b) {
    if constexpr (std::is_floating_point_v<T>) {
        using U = bits_t<T>;
        U ua, ub;
        std::memcpy(&ua, &a, sizeof a);
        std::memcpy(&ub, &b, sizeof b);
        const U m = U(0) - static_cast<U>(c);
        const U ur = (ua & m) | (ub & ~m);
        T r;
        std::memcpy(&r, &ur, sizeof r);
        return r;
    } else {
        return c ? a : b;
    }
}

template <typename T>
inline T x86_nan_arm(T r, T a, T b) {
    if (r != r && a == a && b == b) {
        using U = bits_t<T>;
        const U u = sizeof(T) == 4 ? U(0xFFC00000u) : U(0xFFF8000000000000ull);
        std::memcpy(&r, &u, sizeof r);
    }
    return r;
}
}  // namespace fp_detail

template <typename I>
inline float int64_to_float(I x) {
    static_assert(sizeof(I) == 8 && std::is_integral_v<I>, "int64_to_float: 64-bit integers");
#ifdef PK5D_FP_ARM
    // |x| >= 2^53: keep the top 53 bits, the last of them set if any bit below was (round to odd)
    using U = uint64_t;
    const bool neg = std::is_signed_v<I> && x < 0;
    U u = neg ? U(0) - static_cast<U>(x) : static_cast<U>(x);
    if (u >> 53) {
        int s = 0;
        while ((u >> (53 + s)) != 0) s++;
        const U low = u & ((U(1) << s) - 1);
        u = (u >> s << s) | (low ? (U(1) << s) : U(0));
    }
    const float f = static_cast<float>(static_cast<double>(u));  // the double is exact
    return neg ? -f : f;
#else
    return static_cast<float>(x);
#endif
}

}  // namespace pk5d

#ifdef PK5D_FP_MSVC_SELECT
#define PK5D_GE_SELECT(a, b) ::pk5d::fp_detail::select_bits((a) >= (b), (a), (b))
#define PK5D_LE_SELECT(a, b) ::pk5d::fp_detail::select_bits((a) <= (b), (a), (b))
#else
#define PK5D_GE_SELECT(a, b) ((a) >= (b) ? (a) : (b))
#define PK5D_LE_SELECT(a, b) ((a) <= (b) ? (a) : (b))
#endif

#ifdef PK5D_FP_ARM
#define PK5D_X86_NAN(r, a, b) ::pk5d::fp_detail::x86_nan_arm((r), (a), (b))
#else
#define PK5D_X86_NAN(r, a, b) (r)
#endif
