"""MATLAB's elementary functions, windows and array shifts, bit for bit on every platform (R2024b;
measured, N05 and S03; ``docs/decisions/0005-sums-and-math-library.md``). The port's one home for
them: other modules take their sines, cosines, powers and so on from here.

MATLAB's double ``sin``, ``cos``, ``tan``, ``atan``, ``atan2`` and ``exp`` are fdlibm's and its
``hypot`` msun's (measured on 318,207 values each): glibc 2.28 and NumPy differ from them in the
last bit for 1.5 to 15 % of the values, and other C libraries and NumPy's SIMD paths differ in
other values, so they are written out here in NumPy (one IEEE operation per C operation, no
fused multiply-adds), which gives the same bits everywhere:

- fdlibm's kernels and argument reduction are public: ``kernel_sin``, ``kernel_cos``,
  ``kernel_tan``, ``rem_pio2``; ``fdlibm_sin`` and ``fdlibm_cos`` are fdlibm's functions as they
  are (``x - x`` for Inf and NaN; MATLAB's NaN for Inf on every CPU).
- ``sin``, ``cos``, ``tan``, ``cot`` (``1/tan``), ``atan``, ``atan2`` (+0 for two zeros),
  ``hypot``, ``exp`` (e rounded at 1) and ``power`` (fdlibm's ``e_pow.c``) are MATLAB's;
  ``cos``, ``tan`` and ``cot`` give the positive NaN for Inf and NaN, as MATLAB. fdlibm reduces
  arguments above 2^20 pi/2 by Payne-Hanek, which is not written out: there the C library's
  sine and cosine are used (they may differ from MATLAB's in the last bit).
- ``sind`` and ``cosd`` (S03's rule, measured on 1,979 angles and a dense grid of 59,451): exact
  at multiples of 30 degrees (with sqrt(3)/2 rounded); otherwise ``fmod(|x|, 360)`` reduced to
  ``d`` in [-45, 45] from the nearest multiple of 90, the cosine kernel at +-45, and fdlibm's
  kernel of ``d*(pi/180)`` with the tail ``d*lo(pi/180)`` (pi/180 to 60 digits); ``sind`` is
  odd (``sind(-0)`` is -0), ``cosd`` even (+0 at odd multiples of 90, the positive NaN for Inf).
- ``sqrt`` is IEEE (exact everywhere); ``power(x, p)`` is ``x*x`` for ``p = 2``, ``sqrt(x)`` for
  0.5, ``1/x`` for -1, else fdlibm's pow.
- single data: the double function of the single value, rounded to single. MATLAB's single
  ``sin``, ``exp``, ``cosd``, ... differ from that in the last bit for a few percent of the
  values (as glibc's ``sinf`` does); no Phase 1 code applies them to single data.

``fast_power`` is PetaKit5D's ``fastPower`` (binary powering in the class of the base);
``tukeywin`` MATLAB's Tukey window; ``padarray`` and ``circshift`` move elements only (exact).
"""

from __future__ import annotations

import numpy as np

from petakit5d.core import dtypes

__all__ = [
    "atan",
    "atan2",
    "circshift",
    "cos",
    "cosd",
    "cot",
    "exp",
    "fast_power",
    "fdlibm_cos",
    "fdlibm_sin",
    "hypot",
    "kernel_cos",
    "kernel_sin",
    "kernel_tan",
    "padarray",
    "power",
    "rem_pio2",
    "sin",
    "sind",
    "sqrt",
    "tan",
    "tukeywin",
]

# --- fdlibm and msun in NumPy (E06's, moved here from matlabcompat.neighbors) -----------------
#
# hypot is fdlibm's as FreeBSD's msun has it (the reference platform's C library, glibc 2.28, gives
# the same bits on 24 million argument pairs of every range, and MATLAB on the oracle cases;
# glibc 2.35 and later compute hypot otherwise and differ in the last bit for about one value in
# ten), with MATLAB's NaN for a NaN argument, also next to Inf. atan2, sin and cos are fdlibm's
# (MATLAB's are, bit for bit, measured on 357,997 values; the C library's and NumPy's differ in
# the last bit for a part of them), with MATLAB's +0 for atan2 of two zeros of any sign; sin and
# cos of |x| > 2^20 pi/2 (fdlibm's Payne-Hanek reduction, not ported) are the C library's.
# Element-wise NumPy operations round each operation as the C code does (no fused multiply-adds).
#
# The algorithms and constants of _hypot follow FreeBSD msun's e_hypot.c (fdlibm's, with b below
# 2^-450 scaled, and a and b read again after a subnormal b is scaled), those of _atan, _atan2,
# _rem_pio2, _ksin, _kcos, _sin and _cos fdlibm 5.3 (s_atan.c, e_atan2.c, e_rem_pio2.c, k_sin.c,
# k_cos.c, s_sin.c, s_cos.c); the notice of both is:
#   Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
#   Developed at SunSoft, a Sun Microsystems, Inc. business.
#   Permission to use, copy, modify, and distribute this software is freely granted, provided
#   that this notice is preserved.

_ATANHI = (
    4.63647609000806093515e-01,
    7.85398163397448278999e-01,
    9.82793723247329054082e-01,
    1.57079632679489655800e00,
)
_ATANLO = (
    2.26987774529616870924e-17,
    3.06161699786838301793e-17,
    1.39033110312309984516e-17,
    6.12323399573676603587e-17,
)
_AT = (
    3.33333333333329318027e-01,
    -1.99999999998764832476e-01,
    1.42857142725034663711e-01,
    -1.11111104054623557880e-01,
    9.09088713343650656196e-02,
    -7.69187620504482999495e-02,
    6.66107313738753120669e-02,
    -5.83357013379057348645e-02,
    4.97687799461593236017e-02,
    -3.65315727442169155270e-02,
    1.62858201153657823623e-02,
)
_TINY = 1.0e-300
_PI_O_4 = 7.8539816339744827900e-01
_PI_O_2 = 1.5707963267948965580e00
_PI = 3.1415926535897931160e00
_PI_LO = 1.2246467991473531772e-16
_S = (
    -1.66666666666666324348e-01,
    8.33333333332248946124e-03,
    -1.98412698298579493134e-04,
    2.75573137070700676789e-06,
    -2.50507602534068634195e-08,
    1.58969099521155010221e-10,
)
_C = (
    4.16666666666666019037e-02,
    -1.38888888888741095749e-03,
    2.48015872894767294178e-05,
    -2.75573143513906633035e-07,
    2.08757232129817482790e-09,
    -1.13596475577881948265e-11,
)
_INVPIO2 = 6.36619772367581382433e-01
_PIO2_1 = 1.57079632673412561417e00
_PIO2_1T = 6.07710050650619224932e-11
_PIO2_2 = 6.07710050630396597660e-11
_PIO2_2T = 2.02226624879595063154e-21
_PIO2_3 = 2.02226624871116645580e-21
_PIO2_3T = 8.47842766036889956997e-32
_NPIO2_HW = np.array(
    [
        0x3FF921FB,
        0x400921FB,
        0x4012D97C,
        0x401921FB,
        0x401F6A7A,
        0x4022D97C,
        0x4025FDBB,
        0x402921FB,
        0x402C463A,
        0x402F6A7A,
        0x4031475C,
        0x4032D97C,
        0x40346B9C,
        0x4035FDBB,
        0x40378FDB,
        0x403921FB,
        0x403AB41B,
        0x403C463A,
        0x403DD85A,
        0x403F6A7A,
        0x40407E4C,
        0x4041475C,
        0x4042106C,
        0x4042D97C,
        0x4043A28C,
        0x40446B9C,
        0x404534AC,
        0x4045FDBB,
        0x4046C6CB,
        0x40478FDB,
        0x404858EB,
        0x404921FB,
    ],
    dtype=np.int32,
)


def _hi(x):
    """the signed high 32-bit word of float64 values"""
    return (
        (np.ascontiguousarray(x).view(np.uint64) >> np.uint64(32)).astype(np.uint32).view(np.int32)
    )


def _lo(x):
    return (np.ascontiguousarray(x).view(np.uint64) & np.uint64(0xFFFFFFFF)).astype(np.uint32)


def _with_hi(x, hw):
    """x with its high 32-bit word replaced by hw"""
    lo = np.ascontiguousarray(x).view(np.uint64) & np.uint64(0xFFFFFFFF)
    return ((hw.astype(np.uint32).astype(np.uint64) << np.uint64(32)) | lo).view(np.float64)


def _from_hi(hw):
    """the float64 of high word hw and low word 0"""
    return (hw.astype(np.uint32).astype(np.uint64) << np.uint64(32)).view(np.float64)


_POS_NAN = float(np.array([0x7FF8000000000000], dtype=np.uint64).view(np.float64)[0])
_MATLAB_NAN = float(np.array([0xFFF8000000000000], dtype=np.uint64).view(np.float64)[0])


def _hypot(x, y):
    """MATLAB's hypot: msun's e_hypot.c, +Inf for an infinite argument, and the positive quiet
    NaN (0x7FF8000000000000) when an argument is NaN, whatever its payload and also next to Inf"""
    x, y = np.broadcast_arrays(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
    shape = x.shape
    x = np.ascontiguousarray(x).ravel()
    y = np.ascontiguousarray(y).ravel()
    # a = |the larger|, b = |the smaller| by high word
    ha, hb = _hi(x) & 0x7FFFFFFF, _hi(y) & 0x7FFFFFFF
    sw = hb > ha
    a, b = np.abs(np.where(sw, y, x)), np.abs(np.where(sw, x, y))
    ha, hb = np.where(sw, hb, ha), np.where(sw, ha, hb)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        # a / b > 2^60: a + b
        ratio = (ha - hb) > 0x3C00000
        direct = a + b
        # a > 2^500: both scaled by 2^-600 (k: the exponent of the scale back)
        big = ha > 0x5F300000
        ha = np.where(big, ha - 0x25800000, ha)
        hb = np.where(big, hb - 0x25800000, hb)
        a, b = np.where(big, _with_hi(a, ha), a), np.where(big, _with_hi(b, hb), b)
        # b < 2^-450: a subnormal b with a scaled by 2^1022 and both read again (and swapped
        # if b became the larger), a normal one by 2^600; b == 0: a
        small = hb < 0x23D00000
        zero = small & ((hb | _lo(b).astype(np.int32)) == 0)
        sub = small & ~zero & (hb <= 0x000FFFFF)
        a, b = np.where(sub, a * 2.0**1022, a), np.where(sub, b * 2.0**1022, b)
        ha, hb = np.where(sub, _hi(a), ha), np.where(sub, _hi(b), hb)
        s2 = sub & (hb > ha)
        a, b = np.where(s2, b, a), np.where(s2, a, b)
        ha, hb = np.where(s2, hb, ha), np.where(s2, ha, hb)
        tiny = small & ~zero & ~sub
        ha = np.where(tiny, ha + 0x25800000, ha)
        hb = np.where(tiny, hb + 0x25800000, hb)
        a, b = np.where(tiny, _with_hi(a, ha), a), np.where(tiny, _with_hi(b, hb), b)
        k = np.select([big, sub, tiny], [600, -1022, -600], 0).astype(np.int32)
        # medium size a and b: sqrt(a^2 + b^2) from the products of their 21-bit high parts
        w = a - b
        t1 = _from_hi(ha)
        t2 = a - t1
        far = np.sqrt(t1 * t1 - (b * (-b) - t2 * (a + t1)))
        y1 = _from_hi(hb)
        y2 = b - y1
        t1 = _from_hi(ha + 0x00100000)
        t2 = (a + a) - t1
        near = np.sqrt(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)))
        h = np.where(w > b, far, near)
        h = np.where(k != 0, _from_hi((1023 + k) << 20) * h, h)
        h = np.where(ratio, direct, np.where(zero, a, h))
        h = np.where(np.isinf(x) | np.isinf(y), np.inf, h)
    h = np.where(np.isnan(x) | np.isnan(y), _POS_NAN, h)
    return h.reshape(shape)


def _atan(x):
    """fdlibm atan (s_atan.c)"""
    x = np.asarray(x, dtype=np.float64)
    hx = _hi(x)
    ix = hx & 0x7FFFFFFF
    ax = np.abs(x)
    big = ix >= 0x44100000
    small = ix < 0x3FDC0000
    m0 = ~small & ~big & (ix < 0x3FE60000)
    m1 = ~small & ~big & (ix >= 0x3FE60000) & (ix < 0x3FF30000)
    m2 = ~small & ~big & (ix >= 0x3FF30000) & (ix < 0x40038000)
    m3 = ~small & ~big & (ix >= 0x40038000)
    t = np.where(small, x, 0.0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        t = np.where(m0, (2.0 * ax - 1.0) / (2.0 + ax), t)
        t = np.where(m1, (ax - 1.0) / (ax + 1.0), t)
        t = np.where(m2, (ax - 1.5) / (1.0 + 1.5 * ax), t)
        t = np.where(m3, -1.0 / ax, t)
        z = t * t
        w = z * z
        a = _AT
        s1 = z * (a[0] + w * (a[2] + w * (a[4] + w * (a[6] + w * (a[8] + w * a[10])))))
        s2 = w * (a[1] + w * (a[3] + w * (a[5] + w * (a[7] + w * a[9]))))
        out = t - t * (s1 + s2)
        for k, m in enumerate((m0, m1, m2, m3)):
            zz = _ATANHI[k] - ((t * (s1 + s2) - _ATANLO[k]) - t)
            out = np.where(m, np.where(hx < 0, -zz, zz), out)
    out = np.where(small & (ix < 0x3E200000), x, out)
    out = np.where(big, np.where(hx > 0, _ATANHI[3] + _ATANLO[3], -_ATANHI[3] - _ATANLO[3]), out)
    return np.where(np.isnan(x), x + x, out)


def _atan2(y, x):
    """fdlibm atan2 (e_atan2.c), with MATLAB's +0 for two zeros and MATLAB's NaN for NaN"""
    y, x = np.broadcast_arrays(np.asarray(y, dtype=np.float64), np.asarray(x, dtype=np.float64))
    hx, hy = _hi(x), _hi(y)
    ix, iy = hx & 0x7FFFFFFF, hy & 0x7FFFFFFF
    lx, ly = _lo(x), _lo(y)
    m = ((hy >> 31) & 1) | ((hx >> 30) & 2)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        k = (iy - ix) >> 20
        z = _atan(np.abs(y / x))
        z = np.where(k > 60, _PI_O_2 + 0.5 * _PI_LO, np.where((hx < 0) & (k < -60), 0.0, z))
        out = np.select([m == 0, m == 1, m == 2], [z, -z, _PI - (z - _PI_LO)], (z - _PI_LO) - _PI)
        yinf = iy == 0x7FF00000
        out = np.where(yinf, np.where(hy < 0, -_PI_O_2 - _TINY, _PI_O_2 + _TINY), out)
        xinf = ix == 0x7FF00000
        both = np.choose(
            m, [_PI_O_4 + _TINY, -_PI_O_4 - _TINY, 3.0 * _PI_O_4 + _TINY, -3.0 * _PI_O_4 - _TINY]
        )
        one = np.choose(m, [0.0, -0.0, _PI + _TINY, -_PI - _TINY])
        out = np.where(xinf, np.where(yinf, both, one), out)
        xzero = (ix | lx) == 0
        out = np.where(xzero, np.where(hy < 0, -_PI_O_2 - _TINY, _PI_O_2 + _TINY), out)
        yzero = (iy | ly) == 0
        out = np.where(
            yzero,
            np.choose(m, [y, y, np.full(y.shape, _PI + _TINY), np.full(y.shape, -_PI - _TINY)]),
            out,
        )
        out = np.where(x == 1.0, _atan(y), out)
        out = np.where(np.isnan(x) | np.isnan(y), _MATLAB_NAN, out)
    return np.where(yzero & xzero, 0.0, out)


def _ksin(x, y, iy):
    """fdlibm __kernel_sin (k_sin.c), iy == 0 when y is 0"""
    ix = _hi(x) & 0x7FFFFFFF
    s = _S
    z = x * x
    v = z * x
    r = s[1] + z * (s[2] + z * (s[3] + z * (s[4] + z * s[5])))
    out = x + v * (s[0] + z * r) if iy == 0 else x - ((z * (0.5 * y - v * r) - y) - v * s[0])
    return np.where(ix < 0x3E400000, x, out)


def _kcos(x, y):
    """fdlibm __kernel_cos (k_cos.c)"""
    ix = _hi(x) & 0x7FFFFFFF
    c = _C
    z = x * x
    r = z * (c[0] + z * (c[1] + z * (c[2] + z * (c[3] + z * (c[4] + z * c[5])))))
    small = 1.0 - (0.5 * z - (z * r - x * y))
    qhi = (ix - 0x00200000).astype(np.uint32).astype(np.uint64) << np.uint64(32)
    qx = np.where(ix > 0x3FE90000, 0.28125, qhi.view(np.float64))
    hz = 0.5 * z - qx
    a = 1.0 - qx
    out = np.where(ix < 0x3FD33333, small, a - (hz - (z * r - x * y)))
    return np.where(ix < 0x3E400000, 1.0, out)


def _rem_pio2(x):
    """fdlibm __ieee754_rem_pio2 (e_rem_pio2.c) up to |x| = 2^20 pi/2: (n, y0, y1)"""
    hx = _hi(x)
    ix = hx & 0x7FFFFFFF
    n = np.zeros(x.shape, np.int64)
    y0 = x.copy()
    y1 = np.zeros(x.shape)
    # |x| < 3pi/4: one step
    m = (ix > 0x3FE921FB) & (ix < 0x4002D97C)
    if np.any(m):
        sg = np.where(hx[m] > 0, 1.0, -1.0)
        xm = x[m]
        z = xm - sg * _PIO2_1
        exact = ix[m] == 0x3FF921FB  # the 33+53 bits of pi/2 are not enough near pi/2
        z2 = z - sg * _PIO2_2
        a0 = np.where(exact, z2 - sg * _PIO2_2T, z - sg * _PIO2_1T)
        a1 = np.where(exact, (z2 - a0) - sg * _PIO2_2T, (z - a0) - sg * _PIO2_1T)
        y0[m], y1[m], n[m] = a0, a1, sg.astype(np.int64)
    m = (ix >= 0x4002D97C) & (ix <= 0x413921FB)
    if np.any(m):
        t = np.abs(x[m])
        nn = np.trunc(t * _INVPIO2 + 0.5)
        fn = nn
        r = t - fn * _PIO2_1
        w = fn * _PIO2_1T
        ni = nn.astype(np.int64)
        quick = (ni < 32) & (ix[m] != _NPIO2_HW[np.clip(ni - 1, 0, 31)])
        a0 = r - w
        j = ix[m] >> 20
        i = j - ((_hi(a0) >> 20) & 0x7FF)
        s2 = ~quick & (i > 16)
        t2 = r
        w2 = fn * _PIO2_2
        r2 = t2 - w2
        w2 = fn * _PIO2_2T - ((t2 - r2) - w2)
        b0 = r2 - w2
        i2 = j - ((_hi(b0) >> 20) & 0x7FF)
        s3 = s2 & (i2 > 49)
        t3 = r2
        w3 = fn * _PIO2_3
        r3 = t3 - w3
        w3 = fn * _PIO2_3T - ((t3 - r3) - w3)
        c0 = r3 - w3
        rr = np.where(s3, r3, np.where(s2, r2, r))
        ww = np.where(s3, w3, np.where(s2, w2, w))
        a0 = np.where(s3, c0, np.where(s2, b0, a0))
        a1 = (rr - a0) - ww
        neg = hx[m] < 0
        y0[m] = np.where(neg, -a0, a0)
        y1[m] = np.where(neg, -a1, a1)
        n[m] = np.where(neg, -ni, ni)
    return n, y0, y1


def _sincos(x, which):
    x = np.ascontiguousarray(np.asarray(x, dtype=np.float64))
    ix = _hi(x) & 0x7FFFFFFF
    n, y0, y1 = _rem_pio2(x)
    q = n & 3
    s = _ksin(y0, y1, 1)
    c = _kcos(y0, y1)
    if which == "sin":
        out = np.select([q == 0, q == 1, q == 2], [s, c, -s], -c)
        out = np.where(ix <= 0x3FE921FB, _ksin(x, np.zeros(x.shape), 0), out)
        big = np.sin
    else:
        out = np.select([q == 0, q == 1, q == 2], [c, -s, -c], s)
        out = np.where(ix <= 0x3FE921FB, _kcos(x, np.zeros(x.shape)), out)
        big = np.cos
    large = (ix > 0x413921FB) & (ix < 0x7FF00000)
    if np.any(large):
        out = np.where(large, big(np.where(large, x, 0.0)), out)
    with np.errstate(invalid="ignore"):
        return np.where(ix >= 0x7FF00000, x - x, out)


def _sin(x):
    return _sincos(x, "sin")


def _cos(x):
    return _sincos(x, "cos")


# fdlibm 5.3 (k_tan.c, s_tan.c, e_exp.c, e_pow.c); their notice:
#   Copyright (C) 1993, 2004 by Sun Microsystems, Inc. All rights reserved.
#   Developed at SunSoft, a Sun Microsystems, Inc. business.
#   Permission to use, copy, modify, and distribute this software is freely granted, provided
#   that this notice is preserved.

_T = (
    3.33333333333334091986e-01,
    1.33333333333201242699e-01,
    5.39682539762260521377e-02,
    2.18694882948595424599e-02,
    8.86323982359930005737e-03,
    3.59207910759131235356e-03,
    1.45620945432529025516e-03,
    5.88041240820264096874e-04,
    2.46463134818469906812e-04,
    7.81794442939557092300e-05,
    7.14072491382608190305e-05,
    -1.85586374855275456654e-05,
    2.59073051863633712884e-05,
)
_PIO4 = 7.85398163397448278999e-01
_PIO4LO = 3.06161699786838301793e-17

_LN2HI = 6.93147180369123816490e-01
_LN2LO = 1.90821492927058770002e-10
_INVLN2 = 1.44269504088896338700e00
_P = (
    1.66666666666666019037e-01,
    -2.77777777770155933842e-03,
    6.61375632143793436117e-05,
    -1.65339022054652515390e-06,
    4.13813679705723846039e-08,
)
_O_THRESHOLD = 7.09782712893383973096e02
_U_THRESHOLD = -7.45133219101941108420e02
_TWOM1000 = 9.33263618503218878990e-302


def _f64(x) -> np.ndarray:
    return np.ascontiguousarray(np.asarray(x, dtype=np.float64))


def _with_lo0(x: np.ndarray) -> np.ndarray:
    """x with its low 32-bit word cleared"""
    return (np.ascontiguousarray(x).view(np.uint64) & np.uint64(0xFFFFFFFF00000000)).view(
        np.float64
    )


def _single(fn):
    """the double function fn of single data, rounded to single (module docstring)"""

    def wrapped(*args):
        arrs = [np.asarray(a) for a in args]
        single = any(a.dtype == np.float32 for a in arrs) and not any(
            a.dtype == np.float64 and a.ndim > 0 for a in arrs
        )
        shape = np.broadcast(*arrs).shape
        with np.errstate(all="ignore"):
            r = np.asarray(fn(*[_f64(a) for a in arrs])).reshape(shape)
        r = r.astype(np.float32) if single else r
        return r[()] if r.ndim == 0 else r

    wrapped.__name__ = fn.__name__
    wrapped.__doc__ = fn.__doc__
    return wrapped


# -- fdlibm's kernels and functions as they are ------------------------------------------------


def _vec(*args):
    """float64, C-contiguous, at least 1-D, broadcast; and the shape to return"""
    arrs = np.broadcast_arrays(*(np.asarray(a, dtype=np.float64) for a in args))
    shape = arrs[0].shape
    return [np.ascontiguousarray(np.atleast_1d(a)) for a in arrs], shape


def kernel_sin(x, y, iy: int = 1) -> np.ndarray:
    """fdlibm's ``__kernel_sin(x, y, iy)``: sin(x + y) for |x| <= pi/4, ``y`` the tail of ``x``
    (``iy = 0``: ``y`` is zero and not used)."""
    (x, y), shape = _vec(x, y)
    return _ksin(x, y, iy).reshape(shape)


def kernel_cos(x, y) -> np.ndarray:
    """fdlibm's ``__kernel_cos(x, y)``: cos(x + y) for |x| <= pi/4."""
    (x, y), shape = _vec(x, y)
    return _kcos(x, y).reshape(shape)


def kernel_tan(x, y, iy) -> np.ndarray:
    """fdlibm's ``__kernel_tan(x, y, iy)``: tan(x + y) for ``iy = 1``, -1/tan(x + y) for
    ``iy = -1``, |x| <= pi/4."""
    (x, y, k), shape = _vec(x, y, iy)
    return _ktan(x, y, k.astype(np.int64)).reshape(shape)


def rem_pio2(x):
    """fdlibm's ``__ieee754_rem_pio2`` up to |x| = 2^20 pi/2: ``(n, y0, y1)`` with
    ``x = n*pi/2 + y0 + y1`` (beyond, ``n = 0`` and ``y0 = x``: Payne-Hanek is not written
    out)."""
    (x,), shape = _vec(x)
    n, y0, y1 = _rem_pio2(x)
    return n.reshape(shape), y0.reshape(shape), y1.reshape(shape)


def fdlibm_sin(x) -> np.ndarray:
    """fdlibm's ``sin`` as it is (``x - x`` for NaN; MATLAB's NaN for Inf on every CPU); E06's,
    as ``neighbors.sph2cart`` uses it."""
    x = np.asarray(x, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        return np.where(np.isinf(x), _MATLAB_NAN, _sincos(x, "sin").reshape(x.shape))


def fdlibm_cos(x) -> np.ndarray:
    """fdlibm's ``cos`` as it is (as ``fdlibm_sin``)."""
    x = np.asarray(x, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        return np.where(np.isinf(x), _MATLAB_NAN, _sincos(x, "cos").reshape(x.shape))


# -- sin, cos, tan, cot, atan ------------------------------------------------------------------


@_single
def sin(x):
    """MATLAB's ``sin`` (fdlibm; MATLAB's NaN for Inf, a NaN kept)."""
    return fdlibm_sin(x)


@_single
def cos(x):
    """MATLAB's ``cos`` (fdlibm; the positive NaN for Inf and NaN, as MATLAB)."""
    return _pos_nan(x, _sincos(x, "cos").reshape(np.shape(x)))


def _pos_nan(x, r):
    """r with the positive quiet NaN where x is not finite (MATLAB's cos, tan, cot and cosd)"""
    return np.where(np.isfinite(x), r, _POS_NAN)


def _ktan(x: np.ndarray, y: np.ndarray, iy: np.ndarray) -> np.ndarray:
    """fdlibm __kernel_tan(x, y, iy): tan(x + y) for iy = 1, -1/tan(x + y) for iy = -1"""
    hx = _hi(x)
    ix = hx & 0x7FFFFFFF
    big = ix >= 0x3FE59428  # |x| >= 0.6744
    neg = big & (hx < 0)
    x = np.where(neg, -x, x)
    y = np.where(neg, -y, y)
    xb = np.where(big, (_PIO4 - x) + (_PIO4LO - y), x)
    y = np.where(big, 0.0, y)
    x = xb
    z = x * x
    w = z * z
    t = _T
    r = t[1] + w * (t[3] + w * (t[5] + w * (t[7] + w * (t[9] + w * t[11]))))
    v = z * (t[2] + w * (t[4] + w * (t[6] + w * (t[8] + w * (t[10] + w * t[12])))))
    s = z * x
    r = y + z * (s * (r + v) + y)
    r = r + t[0] * s
    w = x + r
    vb = iy.astype(np.float64)
    sgn = (1 - ((hx >> 30) & 2)).astype(np.float64)
    res_big = sgn * (vb - 2.0 * (x - (w * w / (w + vb) - r)))
    # iy = -1 away from pi/4: -1/(x + r) computed carefully
    z2 = _with_lo0(w)
    v2 = r - (z2 - x)
    a = -1.0 / w
    t2 = _with_lo0(a)
    s2 = 1.0 + t2 * z2
    res_m1 = t2 + a * (s2 + t2 * v2)
    res = np.where(big, res_big, np.where(iy == 1, w, res_m1))
    # |x| < 2^-28
    tiny = ix < 0x3E300000
    if np.any(tiny):
        xt, yt = np.where(tiny, xb, 0.0), np.where(tiny, y, 0.0)
        zero = ((ix | _lo(xt).astype(np.int32)) == 0) & (iy == -1)
        wt = xt + yt
        zt = _with_lo0(wt)
        vt = yt - (zt - xt)
        at = -1.0 / wt
        tt = _with_lo0(at)
        st = 1.0 + tt * zt
        small = np.where(iy == 1, xt, tt + at * (st + tt * vt))
        small = np.where(zero, 1.0 / np.abs(xt), small)
        res = np.where(tiny, small, res)
    return res


def _tan_core(x: np.ndarray, which: str) -> np.ndarray:
    """tan (iy = 1) or -1/tan... of fdlibm's s_tan.c; which 'tan' or 'cot'"""
    x = np.ascontiguousarray(x)
    ix = _hi(x) & 0x7FFFFFFF
    n, y0, y1 = _rem_pio2(x)
    small = ix <= 0x3FE921FB
    y0 = np.where(small, x, y0)
    y1 = np.where(small, 0.0, y1)
    n = np.where(small, 0, n)
    iy = 1 - ((n & 1) << 1)
    with np.errstate(all="ignore"):
        t = _ktan(y0, y1, iy)
        large = (ix > 0x413921FB) & (ix < 0x7FF00000)
        if np.any(large):
            t = np.where(large, np.tan(np.where(large, x, 0.0)), t)
        t = np.where(ix >= 0x7FF00000, x - x, t)
    return t


@_single
def tan(x):
    """MATLAB's ``tan`` (fdlibm; the positive NaN for Inf and NaN)."""
    return _pos_nan(x, _tan_core(x, "tan"))


@_single
def cot(x):
    """MATLAB's ``cot``: ``1 ./ tan(x)`` (the positive NaN for Inf and NaN)."""
    with np.errstate(divide="ignore"):
        return _pos_nan(x, 1.0 / _tan_core(x, "tan"))


@_single
def atan(x):
    """MATLAB's ``atan`` (fdlibm)."""
    return _atan(x)


# -- sind, cosd (S03's rule) -------------------------------------------------------------------

_DEG = float.fromhex("0x1.1df46a2529d39p-6")  # pi/180 rounded
_DEG_LO = float.fromhex("0x1.5c1d8becdd291p-62")  # pi/180 - _DEG (pi/180 to 60 digits)
_R3 = float.fromhex("0x1.bb67ae8584caap-1")  # sqrt(3)/2 rounded
# sind and cosd at 0, 30, ..., 330 degrees
_SIN30 = (0.0, 0.5, _R3, 1.0, _R3, 0.5, 0.0, -0.5, -_R3, -1.0, -_R3, -0.5)
_COS30 = (1.0, _R3, 0.5, 0.0, -0.5, -_R3, -1.0, -_R3, -0.5, 0.0, 0.5, _R3)


def _deg(x, sine: bool):
    x = _f64(x)
    with np.errstate(all="ignore"):
        a = np.fmod(np.abs(x), 360.0)  # exact
        # round(a/90), half away from zero and exact: floor(a/90 + 0.5) rounds the sum up just
        # below 45 (44.99999999999999 would leave [-45, 45]); MATLAB's cosd there is the cosine
        n = np.asarray(dtypes.mround(a / 90.0))
        d = a - 90.0 * n  # exact, in [-45, 45]
        # at +-45 degrees the cosine kernel: an odd quadrant for sind, an even one for cosd
        swap = (np.abs(d) == 45.0) & (np.mod(n, 2) != (1 if sine else 0))
        n = np.where(swap, np.where(d < 0, n - 1, n + 1), n)
        d = np.where(swap, np.where(d < 0, 45.0, -45.0), d)
        t, tl = np.ascontiguousarray(d * _DEG), np.ascontiguousarray(d * _DEG_LO)
        s, c = _ksin(t, tl, 1), _kcos(t, tl)
        q = np.mod(n, 4).astype(np.int64)
        if sine:
            v = np.select([q == 0, q == 1, q == 2], [s, c, -s], -c)
        else:
            v = np.select([q == 0, q == 1, q == 2], [c, -s, -c], s)
        k = np.where((a == np.floor(a)) & (np.mod(a, 30.0) == 0), a / 30.0, -1.0)
        k = np.nan_to_num(k, nan=-1.0).astype(np.int64)
        table = np.asarray(_SIN30 if sine else _COS30)
        v = np.where(k >= 0, table[np.clip(k, 0, 11)], v)
        if sine:
            v = np.copysign(1.0, x) * v  # odd: sind(-0) is -0
        return np.where(np.isfinite(x), v, _MATLAB_NAN if sine else _POS_NAN)


@_single
def sind(x):
    """MATLAB's ``sind`` (module docstring)."""
    return _deg(x, True)


@_single
def cosd(x):
    """MATLAB's ``cosd`` (module docstring)."""
    return _deg(x, False)


# -- exp ---------------------------------------------------------------------------------------


def _exp(x: np.ndarray) -> np.ndarray:
    """fdlibm __ieee754_exp"""
    x = np.ascontiguousarray(x)
    hxs = _hi(x)
    xsb = (hxs >> 31) & 1
    hx = hxs & 0x7FFFFFFF
    with np.errstate(all="ignore"):
        mid = (hx > 0x3FD62E42) & (hx < 0x3FF0A2B2)
        far = hx >= 0x3FF0A2B2
        sgn = np.where(xsb == 1, -1.0, 1.0)
        k = np.zeros(x.shape, np.int64)
        hi = x.copy()
        lo = np.zeros(x.shape)
        # |x| in (0.5 ln2, 1.5 ln2)
        k = np.where(mid, 1 - xsb - xsb, k)
        hi = np.where(mid, x - sgn * _LN2HI, hi)
        lo = np.where(mid, sgn * _LN2LO, lo)
        # |x| >= 1.5 ln2
        kf = np.trunc(_INVLN2 * x + np.where(xsb == 1, -0.5, 0.5))
        k = np.where(far, kf.astype(np.int64), k)
        hi = np.where(far, x - kf * _LN2HI, hi)
        lo = np.where(far, kf * _LN2LO, lo)
        red = mid | far
        xr = np.where(red, hi - lo, x)
        t = xr * xr
        c = xr - t * (_P[0] + t * (_P[1] + t * (_P[2] + t * (_P[3] + t * _P[4]))))
        y0 = 1.0 - ((xr * c) / (c - 2.0) - xr)
        y = 1.0 - ((lo - (xr * c) / (2.0 - c)) - hi)
        kk = np.where(red, k, 0)
        normal = kk >= -1021
        yh = _hi(y).astype(np.int64)
        y_n = _with_hi(y, (yh + (kk << 20)).astype(np.int64))
        y_s = _with_hi(y, (yh + ((kk + 1000) << 20)).astype(np.int64)) * _TWOM1000
        out = np.where(red, np.where(normal, y_n, y_s), y0)
        tiny = hx < 0x3E300000
        out = np.where(tiny, 1.0 + x, out)
        out = np.where(x > _O_THRESHOLD, np.inf, out)
        out = np.where(x < _U_THRESHOLD, 0.0, out)
        out = np.where(np.isinf(x), np.where(x > 0, np.inf, 0.0), out)
        out = np.where(np.isnan(x), x + x, out)
    return out


_E = 2.718281828459045


@_single
def exp(x):
    """MATLAB's ``exp``: fdlibm's, except ``exp(1)``, which is e rounded (fdlibm's is 1 ulp
    above; measured on 107,493 values)."""
    r = _exp(x)
    return np.where(np.asarray(x) == 1.0, _E, r)


# -- sqrt, hypot, atan2 ------------------------------------------------------------------------


def sqrt(x):
    """MATLAB's ``sqrt`` of real data: the IEEE square root (exact on every platform); a negative
    value gives a complex result, as in MATLAB."""
    a = np.asarray(x)
    if a.dtype.kind in "iu" or a.dtype.kind == "b":
        a = a.astype(np.float64)
    if np.any(a < 0):
        return np.sqrt(a.astype(np.complex64 if a.dtype == np.float32 else np.complex128))
    return np.sqrt(a)


@_single
def hypot(x, y):
    """MATLAB's ``hypot`` (msun)."""
    return _hypot(x, y)


@_single
def atan2(y, x):
    """MATLAB's ``atan2`` (fdlibm; +0 for two zeros)."""
    return _atan2(y, x)


# -- power -------------------------------------------------------------------------------------

_BP = (1.0, 1.5)
_DP_H = (0.0, 5.84962487220764160156e-01)
_DP_L = (0.0, 1.35003920212974897128e-08)
_L = (
    5.99999999999994648725e-01,
    4.28571428578550184252e-01,
    3.33333329818377432918e-01,
    2.72728123808534006489e-01,
    2.30660745775561754067e-01,
    2.06975017800338417784e-01,
)
_PP = (
    1.66666666666666019037e-01,
    -2.77777777770155933842e-03,
    6.61375632143793436117e-05,
    -1.65339022054652515390e-06,
    4.13813679705723846039e-08,
)
_LG2 = 6.93147180559945286227e-01
_LG2_H = 6.93147182464599609375e-01
_LG2_L = -1.90465429995776804525e-09
_OVT = 8.0085662595372944372e-17
_CP = 9.61796693925975554329e-01
_CP_H = 9.61796700954437255859e-01
_CP_L = -7.02846165095275826516e-09
_IVLN2 = 1.44269504088896338700e00
_IVLN2_H = 1.44269502162933349609e00
_IVLN2_L = 1.92596299112661746887e-08
_TWO53 = 9007199254740992.0


def _pow_scalar(x: float, y: float) -> float:
    """fdlibm __ieee754_pow for one pair"""
    import struct

    def words(v):
        b = struct.unpack("<Q", struct.pack("<d", v))[0]
        hi = b >> 32
        return (hi - (1 << 32) if hi & 0x80000000 else hi), b & 0xFFFFFFFF

    def make(hi, lo):
        return struct.unpack(
            "<d", struct.pack("<Q", ((hi & 0xFFFFFFFF) << 32) | (lo & 0xFFFFFFFF))
        )[0]

    def set_lo0(v):
        h, _ = words(v)
        return make(h, 0)

    hx, lx = words(x)
    hy, ly = words(y)
    ix, iy = hx & 0x7FFFFFFF, hy & 0x7FFFFFFF
    if (iy | ly) == 0:
        return 1.0
    if (
        ix > 0x7FF00000
        or (ix == 0x7FF00000 and lx != 0)
        or iy > 0x7FF00000
        or (iy == 0x7FF00000 and ly != 0)
    ):
        return x + y
    yisint = 0
    if hx < 0:
        if iy >= 0x43400000:
            yisint = 2
        elif iy >= 0x3FF00000:
            k = (iy >> 20) - 0x3FF
            if k > 20:
                j = ly >> (52 - k)
                if (j << (52 - k)) == ly:
                    yisint = 2 - (j & 1)
            elif ly == 0:
                j = iy >> (20 - k)
                if (j << (20 - k)) == iy:
                    yisint = 2 - (j & 1)
    if ly == 0:
        if iy == 0x7FF00000:
            if ((ix - 0x3FF00000) | lx) == 0:
                return y - y
            if ix >= 0x3FF00000:
                return y if hy >= 0 else 0.0
            return -y if hy < 0 else 0.0
        if iy == 0x3FF00000:
            return 1.0 / x if hy < 0 else x
        if hy == 0x40000000:
            return x * x
        if hy == 0x3FE00000 and hx >= 0:
            return float(np.sqrt(x))
    ax = abs(x)
    if lx == 0 and (ix == 0x7FF00000 or ix == 0 or ix == 0x3FF00000):
        z = ax
        if hy < 0:
            z = 1.0 / z
        if hx < 0:
            if ((ix - 0x3FF00000) | yisint) == 0:
                z = (z - z) / (z - z)
            elif yisint == 1:
                z = -z
        return z
    n = (hx >> 31) + 1
    if (n | yisint) == 0:
        return (x - x) / (x - x)
    s = 1.0
    if (n | (yisint - 1)) == 0:
        s = -1.0
    if iy > 0x41E00000:
        if iy > 0x43F00000:
            if ix <= 0x3FEFFFFF:
                return 1e300 * 1e300 if hy < 0 else 1e-300 * 1e-300
            if ix >= 0x3FF00000:
                return 1e300 * 1e300 if hy > 0 else 1e-300 * 1e-300
        if ix < 0x3FEFFFFF:
            return s * 1e300 * 1e300 if hy < 0 else s * 1e-300 * 1e-300
        if ix > 0x3FF00000:
            return s * 1e300 * 1e300 if hy > 0 else s * 1e-300 * 1e-300
        t = ax - 1.0
        w = (t * t) * (0.5 - t * (0.3333333333333333333333 - t * 0.25))
        u = _IVLN2_H * t
        v = t * _IVLN2_L - w * _IVLN2
        t1 = set_lo0(u + v)
        t2 = v - (t1 - u)
    else:
        n = 0
        if ix < 0x00100000:
            ax *= _TWO53
            n -= 53
            ix, _ = words(ax)
        n += ((ix) >> 20) - 0x3FF
        j = ix & 0x000FFFFF
        ix = j | 0x3FF00000
        if j <= 0x3988E:
            k = 0
        elif j < 0xBB67A:
            k = 1
        else:
            k = 0
            n += 1
            ix -= 0x00100000
        _, axlo = words(ax)
        ax = make(ix, axlo)
        u = ax - _BP[k]
        v = 1.0 / (ax + _BP[k])
        ss = u * v
        s_h = set_lo0(ss)
        t_h = make(((ix >> 1) | 0x20000000) + 0x00080000 + (k << 18), 0)
        t_l = ax - (t_h - _BP[k])
        s_l = v * ((u - s_h * t_h) - s_h * t_l)
        s2 = ss * ss
        r = (
            s2
            * s2
            * (_L[0] + s2 * (_L[1] + s2 * (_L[2] + s2 * (_L[3] + s2 * (_L[4] + s2 * _L[5])))))
        )
        r += s_l * (s_h + ss)
        s2 = s_h * s_h
        t_h = set_lo0(3.0 + s2 + r)
        t_l = r - ((t_h - 3.0) - s2)
        u = s_h * t_h
        v = s_l * t_h + t_l * ss
        p_h = set_lo0(u + v)
        p_l = v - (p_h - u)
        z_h = _CP_H * p_h
        z_l = _CP_L * p_h + p_l * _CP + _DP_L[k]
        t = float(n)
        t1 = set_lo0(((z_h + z_l) + _DP_H[k]) + t)
        t2 = z_l - (((t1 - t) - _DP_H[k]) - z_h)
    y1 = set_lo0(y)
    p_l = (y - y1) * t1 + y * t2
    p_h = y1 * t1
    z = p_l + p_h
    j, i = words(z)
    if j >= 0x40900000:
        if ((j - 0x40900000) | i) != 0:
            return s * 1e300 * 1e300
        if p_l + _OVT > z - p_h:
            return s * 1e300 * 1e300
    elif (j & 0x7FFFFFFF) >= 0x4090CC00:
        if ((j - 0xC090CC00) | i) != 0:
            return s * 1e-300 * 1e-300
        if p_l <= z - p_h:
            return s * 1e-300 * 1e-300
    i = j & 0x7FFFFFFF
    k = (i >> 20) - 0x3FF
    n = 0
    if i > 0x3FE00000:
        n = j + (0x00100000 >> (k + 1))
        k = ((n & 0x7FFFFFFF) >> 20) - 0x3FF
        t = make(n & ~(0x000FFFFF >> k), 0)
        n = ((n & 0x000FFFFF) | 0x00100000) >> (20 - k)
        if j < 0:
            n = -n
        p_h -= t
    t = set_lo0(p_l + p_h)
    u = t * _LG2_H
    v = (p_l - (t - p_h)) * _LG2 + t * _LG2_L
    z = u + v
    w = v - (z - u)
    t = z * z
    t1 = z - t * (_PP[0] + t * (_PP[1] + t * (_PP[2] + t * (_PP[3] + t * _PP[4]))))
    r = (z * t1) / (t1 - 2.0) - (w + z * w)
    z = 1.0 - (r - z)
    j, zl = words(z)
    j += n << 20
    if (j >> 20) <= 0:
        z = float(np.ldexp(z, n))
    else:
        z = make(j, zl)
    return s * z


def _pow(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x, y = np.broadcast_arrays(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
    out = np.empty(x.shape)
    fx, fy, fo = x.reshape(-1), y.reshape(-1), out.reshape(-1)
    with np.errstate(all="ignore"):
        for i in range(fx.size):
            fo[i] = _pow_scalar(float(fx[i]), float(fy[i]))
    return out[()] if out.ndim == 0 else out


def power(x, p):
    """MATLAB's ``x .^ p`` of real data (module docstring): ``p = 2`` is ``x .* x``, 0.5 ``sqrt``,
    -1 ``1 ./ x``, 0 one and 1 ``x``; any other exponent fdlibm's pow (a negative base with a
    non-integer exponent gives NaN here, where MATLAB returns a complex result). Single data:
    the double result rounded to single, except for the exact cases above."""
    a, b = np.asarray(x), np.asarray(p)
    single = (a.dtype == np.float32 or b.dtype == np.float32) and not (
        (a.dtype == np.float64 and a.ndim > 0) or (b.dtype == np.float64 and b.ndim > 0)
    )
    dt = np.float32 if single else np.float64
    if b.ndim == 0:
        e = float(b)
        A = a.astype(dt)
        if e == 2:
            return np.multiply(A, A)
        if e == 1:
            return A.copy()
        if e == 0:
            return np.ones_like(A)
        if e == -1:
            with np.errstate(divide="ignore"):
                return np.divide(dt(1), A)
        if e == 0.5 and not np.any(A < 0):
            return np.sqrt(A)
    r = _pow(a, b)
    return r.astype(dt) if single else r


def fast_power(base, exponent):
    """PetaKit5D's ``fastPower(base, exponent)`` for an integer exponent: binary powering in the
    class of ``base`` (``result .* base`` for each set bit, ``base .* base`` per bit), the base
    inverted first for a negative exponent; exponent 1 returns the (inverted) base and exponent 0
    the double scalar 1, as the MATLAB code does."""
    e = int(exponent)
    if e != exponent:
        raise ValueError("fastPower takes an integer exponent")
    b = np.asarray(base)
    if e < 0:
        b = np.asarray(dtypes.rdivide(1.0, b))
        e = -e
    if e == 1:
        return b
    result = 1.0
    while e > 0:
        if e & 1:
            result = dtypes.times(result, b)
        b = np.asarray(dtypes.times(b, b))
        e >>= 1
    return np.asarray(result) if not isinstance(result, float) else result


# -- windows -----------------------------------------------------------------------------------


def tukeywin(n: int, r: float = 0.5) -> np.ndarray:
    """MATLAB's ``tukeywin(n, r)`` (a column, as a 1-D array of n values): ones for ``r <= 0``;
    for ``0 < r < 1``, with ``t = (0:n-1)/(n-1)`` and ``per = r/2``,
    ``(1 + cos(pi/per*(t - per)))/2`` on the first ``floor(per*(n - 1)) + 1`` points,
    ``(1 + cos(pi/per*(t - 1 + per)))/2`` on as many last points and one between (exact on
    every case measured); for ``r >= 1`` the Hann window, ``0.5 - 0.5*cos(2*pi*t)`` on the first
    half, mirrored (MATLAB's differs from it in the last bit for some ``n`` above 11)."""
    n = int(n)
    if n <= 0:
        return np.zeros(0)
    if n == 1:
        return np.ones(1)
    r = float(r)
    t = np.arange(n, dtype=np.float64) / (n - 1)
    if r <= 0:
        return np.ones(n)
    if r >= 1:
        half = (n + 1) // 2
        w = 0.5 - 0.5 * _cos(2 * np.pi * t[:half])
        return np.concatenate([w, w[: n - half][::-1]])
    per = r / 2
    tl = int(np.floor(per * (n - 1))) + 1
    w = np.ones(n)
    w[:tl] = (1 + _cos(np.pi / per * (t[:tl] - per))) / 2
    w[n - tl :] = (1 + _cos(np.pi / per * (t[n - tl :] - 1 + per))) / 2
    return w


# -- padarray, circshift -----------------------------------------------------------------------

_PAD_MODES = {"replicate": "edge", "symmetric": "symmetric", "circular": "wrap"}


def padarray(a, pad_size, value=0, *, method: str = "constant", direction: str = "both"):
    """MATLAB's ``padarray(A, padsize, padval | method, direction)`` with ``pad_size`` per axis of
    ``a`` (the port's ``(z, y, x)``): ``method`` 'constant' (``value`` in the class of ``a``, as
    MATLAB converts it), 'replicate', 'symmetric' or 'circular'; ``direction`` 'both', 'pre' or
    'post'."""
    x = np.asarray(a)
    ps = [int(p) for p in np.atleast_1d(pad_size)]
    if len(ps) < x.ndim:
        ps = [0] * (x.ndim - len(ps)) + ps
    if len(ps) != x.ndim or any(p < 0 for p in ps):
        raise ValueError("pad_size needs one non-negative size per axis")
    if direction not in ("both", "pre", "post"):
        raise ValueError("direction must be 'both', 'pre' or 'post'")
    widths = [
        (p if direction in ("both", "pre") else 0, p if direction in ("both", "post") else 0)
        for p in ps
    ]
    if method == "constant":
        if x.dtype.kind in "iu":
            v = np.asarray(dtypes.cast(float(value), x.dtype.name))
        else:
            v = np.asarray(value).astype(x.dtype)
        return np.pad(x, widths, mode="constant", constant_values=v)
    if method not in _PAD_MODES:
        raise ValueError(f"unknown method {method!r}")
    if x.size == 0:
        raise ValueError(f"{method} padding of an empty array")
    return np.pad(x, widths, mode=_PAD_MODES[method])


def circshift(a, shift, axis=None):
    """MATLAB's ``circshift(A, K, dim)`` on the given NumPy axes: element ``i`` moves to
    ``i + K`` (mod the size); ``axis=None`` with a scalar ``shift``: MATLAB's first non-singleton
    dimension (``reduce.default_axis``)."""
    x = np.asarray(a)
    if axis is None:
        if np.ndim(shift) == 0:
            from petakit5d.matlabcompat.reduce import default_axis

            axis = default_axis(x)
        else:
            raise ValueError("give the axes of a vector of shifts")
    return np.roll(x, shift, axis=axis)
