// x86's NaN in the deskew/rotate kernels' sums, for the builds without the mex's AVX2 code (other
// processors, and deskew_scalar on x86-64): a float value that comes out NaN is computed again
// with x86's NaN for each operation (pk5d::x86_nans) of the mex's formula as written (the products
// a * wa, the sums from the left): MATLAB's NaN for Inf * 0 and Inf - Inf, and where NaNs meet the
// first operand's. ARM's default NaN has the sign bit clear. A NaN never becomes a number in these
// sums, so the other values keep their bits. The AVX2 build keeps the mex's code: x86's NaN, and
// where NaNs meet the compiler's choice of operand. uint16 data hold no NaN or Inf.
#pragma once

#if !defined(__AVX2__)
#define PK5D_DESKEW_NAN_X86 1

#include <cstddef>
#include <cstdint>

#include <pk5d/fpexact.hpp>

namespace pk5d::deskew {

inline float mul_x86(float a, float w) { return pk5d::x86_nans(a * w, a, w); }
inline float add_x86(float a, float b) { return pk5d::x86_nans(a + b, a, b); }

// a * wa + b * wb
inline float sum2_x86(float a, float b, float wa, float wb) {
    return add_x86(mul_x86(a, wa), mul_x86(b, wb));
}

// ((a * wa + b * wb) + c * wc) + d * wd (SEQ), the float lines' sum
inline float seq4_x86(float a, float b, float c, float d, float wa, float wb, float wc, float wd) {
    return add_x86(add_x86(add_x86(mul_x86(a, wa), mul_x86(b, wb)), mul_x86(c, wc)), mul_x86(d, wd));
}

inline bool any_nan(const float* r, size_t n) {
    bool any = false;
    for (size_t i = 0; i < n; i++) any |= r[i] != r[i];
    return any;
}

// r[i] = a[i] * wa + b[i] * wb for i < n: the values that came out NaN, computed again
inline void nans_x86_2(const float* a, const float* b, float* r, float wa, float wb, size_t n) {
    if (!any_nan(r, n)) return;
    for (size_t i = 0; i < n; i++)
        if (r[i] != r[i]) r[i] = sum2_x86(a[i], b[i], wa, wb);
}

// r[i] = SEQ of a[i], b[i], c[i], d[i] for i < n: the values that came out NaN, computed again
inline void nans_x86_4(const float* a, const float* b, const float* c, const float* d, float* r,
                       float wa, float wb, float wc, float wd, size_t n) {
    if (!any_nan(r, n)) return;
    for (size_t i = 0; i < n; i++)
        if (r[i] != r[i]) r[i] = seq4_x86(a[i], b[i], c[i], d[i], wa, wb, wc, wd);
}

}  // namespace pk5d::deskew

#endif
