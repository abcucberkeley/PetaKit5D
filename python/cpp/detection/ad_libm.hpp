// The C library's erf, erfc and log of the Anderson-Darling statistics (ad_stat, fit_mex), which
// adtest_mex and the fitGaussian3D mex take from theirs. With MinGW-w64 they are the UCRT's,
// through their imports, as MSVC's build has them (docs/decisions/0011-windows-toolchain.md):
// MinGW-w64's own (libmingwex's x87 code, run at 64 bits by cpp/common/mingw_libm.c) took 62 ns
// for erf and 29 for log where the UCRT's take 18 and 4, most of the time of a fit of A and c.
// All three go together: libmingwex has erf and erfc in one object, and the import library's
// members define erf and erfc too, so a module may not take one of them from each.
// Elsewhere std::erf, std::erfc, std::log.
#pragma once

#include <cmath>

#if defined(__MINGW32__)
extern "C" {
extern double (*__imp_erf)(double);
extern double (*__imp_erfc)(double);
extern double (*__imp_log)(double);
}
#endif

namespace pk5d_ad {

#if defined(__MINGW32__)
inline double erf(double x) { return __imp_erf(x); }
inline double erfc(double x) { return __imp_erfc(x); }
inline double log(double x) { return __imp_log(x); }
#else
inline double erf(double x) { return std::erf(x); }
inline double erfc(double x) { return std::erfc(x); }
inline double log(double x) { return std::log(x); }
#endif

}  // namespace pk5d_ad
