/* MinGW-w64 only: the C math functions that the modules call, made independent of the calling thread
   and right (docs/decisions/0011-windows-toolchain.md). cmake/Pk5dExtension.cmake compiles this file once
   per function, with PK5D_LIBM_<kind> and PK5D_LIBM_NAME defined, into a static library, so that a module
   carries only the ones it calls, and links the modules with --wrap=<name>, which sends their calls of
   <name> to __wrap_<name> (and __real_<name> to MinGW-w64's own).

   - MinGW-w64's own functions (libmingwex) compute on the x87 unit at the precision of the calling thread's
     control word: 53 bits in a Python thread (Windows' default), 64 in OpenMP's threads (winpthreads sets
     it), so that exp, pow, erfc and others gave other bits on other threads (measured: exp in 17% of random
     arguments, pow in 48%). __wrap_<name> runs MinGW-w64's at 64 bits, MinGW's and Linux's precision, where
     its results come closest to glibc's; and an invalid operation gives x86's default NaN (sign bit set,
     as glibc and the UCRT give it: MinGW-w64's log(-1) and sin(Inf) were positive NaNs).
   - fma, fmaf: MinGW-w64's are not fused (they differed from the instruction in 11% and 22% of random
     near-cancelling inputs): the UCRT's, which are.
   - round, lround, llround and their float forms: MinGW-w64 10's round gives 1 for 0.49999999999999994 (the
     largest double below one half): halves are rounded away from zero exactly.
   - nearbyint, ldexp: exact at any precision, but MinGW-w64's run the x87 unit (nearbyint with its control
     word saved and set around frndint: 27 ns a call, against 1.4 here; ldexp 5.8 against 1.4), and the
     detection fits' correctly rounded exp calls them four times a value: the same results in SSE2. */
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

/* ("memory": the call stays between the two) */
static inline unsigned short pk5d_x87_enter(void) {
    unsigned short saved, cw;
    __asm__ __volatile__("fnstcw %0" : "=m"(saved) : : "memory");
    cw = (unsigned short)(saved | 0x0300); /* precision control: 64 bits */
    __asm__ __volatile__("fldcw %0" : : "m"(cw) : "memory");
    return saved;
}

static inline void pk5d_x87_leave(unsigned short saved) {
    __asm__ __volatile__("fldcw %0" : : "m"(saved) : "memory");
}

static inline double pk5d_x86_nan(void) {
    const uint64_t u = 0xFFF8000000000000ull;
    double r;
    memcpy(&r, &u, sizeof r);
    return r;
}

static inline float pk5d_x86_nanf(void) {
    const uint32_t u = 0xFFC00000u;
    float r;
    memcpy(&r, &u, sizeof r);
    return r;
}

#define PK5D_CAT2(a, b) a##b
#define PK5D_CAT(a, b) PK5D_CAT2(a, b)
#define PK5D_REAL PK5D_CAT(__real_, PK5D_LIBM_NAME)
#define PK5D_WRAP PK5D_CAT(__wrap_, PK5D_LIBM_NAME)

#if defined(PK5D_LIBM_D_D)
double PK5D_REAL(double);
double PK5D_WRAP(double x) {
    const unsigned short cw = pk5d_x87_enter();
    double r = PK5D_REAL(x);
    pk5d_x87_leave(cw);
    return (r != r && x == x) ? pk5d_x86_nan() : r;
}
#elif defined(PK5D_LIBM_D_DD)
double PK5D_REAL(double, double);
double PK5D_WRAP(double x, double y) {
    const unsigned short cw = pk5d_x87_enter();
    double r = PK5D_REAL(x, y);
    pk5d_x87_leave(cw);
    return (r != r && x == x && y == y) ? pk5d_x86_nan() : r;
}
#elif defined(PK5D_LIBM_F_F)
float PK5D_REAL(float);
float PK5D_WRAP(float x) {
    const unsigned short cw = pk5d_x87_enter();
    float r = PK5D_REAL(x);
    pk5d_x87_leave(cw);
    return (r != r && x == x) ? pk5d_x86_nanf() : r;
}
#elif defined(PK5D_LIBM_F_FF)
float PK5D_REAL(float, float);
float PK5D_WRAP(float x, float y) {
    const unsigned short cw = pk5d_x87_enter();
    float r = PK5D_REAL(x, y);
    pk5d_x87_leave(cw);
    return (r != r && x == x && y == y) ? pk5d_x86_nanf() : r;
}
#elif defined(PK5D_LIBM_SINCOS)
void PK5D_REAL(double, double *, double *);
void PK5D_WRAP(double x, double *s, double *c) {
    const unsigned short cw = pk5d_x87_enter();
    PK5D_REAL(x, s, c);
    pk5d_x87_leave(cw);
    if (x == x) {
        if (*s != *s) *s = pk5d_x86_nan();
        if (*c != *c) *c = pk5d_x86_nan();
    }
}
#elif defined(PK5D_LIBM_SINCOSF)
void PK5D_REAL(float, float *, float *);
void PK5D_WRAP(float x, float *s, float *c) {
    const unsigned short cw = pk5d_x87_enter();
    PK5D_REAL(x, s, c);
    pk5d_x87_leave(cw);
    if (x == x) {
        if (*s != *s) *s = pk5d_x86_nanf();
        if (*c != *c) *c = pk5d_x86_nanf();
    }
}
#elif defined(PK5D_LIBM_UCRT_FMA)
/* the UCRT's, through its import (MinGW-w64's libucrt.a) */
extern double (*__imp_fma)(double, double, double);
double PK5D_WRAP(double a, double b, double c) { return __imp_fma(a, b, c); }
#elif defined(PK5D_LIBM_UCRT_FMAF)
extern float (*__imp_fmaf)(float, float, float);
float PK5D_WRAP(float a, float b, float c) { return __imp_fmaf(a, b, c); }
#elif defined(PK5D_LIBM_ROUND) || defined(PK5D_LIBM_LROUND) || defined(PK5D_LIBM_LLROUND)
/* x - trunc(x) is exact (|x| < 2^52), so is the comparison with one half */
static double pk5d_round(double x) {
    if (!(fabs(x) < 4503599627370496.0)) return x + 0.0; /* integers, Inf, NaN (quieted) */
    double t = copysign((double)(long long)x, x);
    const double d = x - t;
    if (d >= 0.5)
        t += 1.0;
    else if (d <= -0.5)
        t -= 1.0;
    return t;
}
#if defined(PK5D_LIBM_ROUND)
double PK5D_WRAP(double x) { return pk5d_round(x); }
#elif defined(PK5D_LIBM_LROUND)
long PK5D_WRAP(double x) { return (long)pk5d_round(x); }
#else
long long PK5D_WRAP(double x) { return (long long)pk5d_round(x); }
#endif
#elif defined(PK5D_LIBM_ROUNDF) || defined(PK5D_LIBM_LROUNDF) || defined(PK5D_LIBM_LLROUNDF)
static float pk5d_roundf(float x) {
    if (!(fabsf(x) < 8388608.0f)) return x + 0.0f;
    float t = copysignf((float)(long)x, x);
    const float d = x - t;
    if (d >= 0.5f)
        t += 1.0f;
    else if (d <= -0.5f)
        t -= 1.0f;
    return t;
}
#if defined(PK5D_LIBM_ROUNDF)
float PK5D_WRAP(float x) { return pk5d_roundf(x); }
#elif defined(PK5D_LIBM_LROUNDF)
long PK5D_WRAP(float x) { return (long)pk5d_roundf(x); }
#else
long long PK5D_WRAP(float x) { return (long long)pk5d_roundf(x); }
#endif
#elif defined(PK5D_LIBM_NEARBYINT)
/* in the current rounding mode: from 2^52 on every double is an integer, and below it the sum with 2^52 (of
   x's sign) rounds x to one; NaN quieted, as the x87 unit's frndint returns it */
double PK5D_WRAP(double x) {
    const double big = 4503599627370496.0; /* 2^52 */
    if (x != x) return x + x;
    if (!(fabs(x) < big)) return x;
    const double r = x >= 0 ? (x + big) - big : (x - big) + big;
    return copysign(r, x); /* -0 for -0.5 <= x <= -0 */
}
#elif defined(PK5D_LIBM_LDEXP)
/* x 2^k by one multiplication when 2^k is a normal double (one rounding, as the x87 unit's fscale and store
   give it), MinGW-w64's own otherwise; NaN, infinities and zeros as they are, and errno ERANGE where the
   result overflows or underflows to zero, as MinGW-w64's */
double PK5D_REAL(double, int);
double PK5D_WRAP(double x, int k) {
    if (!isfinite(x) || x == 0.0) return x;
    if (k < -1022 || k > 1023) return PK5D_REAL(x, k);
    const uint64_t u = (uint64_t)(k + 1023) << 52;
    double p;
    memcpy(&p, &u, sizeof p);
    const double r = x * p;
    if (r == 0.0 || !isfinite(r)) errno = ERANGE;
    return r;
}
#else
#error "mingw_libm.c: define PK5D_LIBM_<kind> and PK5D_LIBM_NAME (cmake/Pk5dExtension.cmake)"
#endif
