"""MATLAB's elementary functions, windows and array shifts, bit for bit on every platform (R2024b;
measured, N05; ``docs/decisions/0005-sums-and-math-library.md``).

MATLAB's double ``sin``, ``cos``, ``tan``, ``atan``, ``atan2`` and ``exp`` are fdlibm's and its
``hypot`` msun's (measured on 318,207 values each): glibc 2.28 and NumPy differ from them in the
last bit for 1.5 to 15 % of the values, and other C libraries and NumPy's SIMD paths differ in
other values, so they are written out here in NumPy (one IEEE operation per C operation, no
fused multiply-adds), which gives the same bits everywhere:

- ``sin``, ``cos``, ``atan2``, ``hypot``: E06's versions (``matlabcompat.neighbors``); ``tan``
  (fdlibm's ``__kernel_tan``), ``cot`` (``1/tan``), ``exp`` and ``power`` (fdlibm's ``e_pow.c``)
  below. fdlibm reduces arguments above 2^20 pi/2 by Payne-Hanek, which is not written out: there
  the C library's sine and cosine are used (they may differ from MATLAB's in the last bit).
- ``sind``, ``cosd``: the angle reduced exactly in degrees (``fmod(x, 360)``, then to
  [-45, 45] around a multiple of 90), the sine or cosine of ``z*(pi/180)``, MATLAB's rounded
  values at multiples of 30 and 45 degrees; exact zeros at multiples of 180 (``sind``, the sign
  of ``x``) and odd multiples of 90 (``cosd``, +0). This is MATLAB's result at every special
  angle measured and at skew angles such as 31.5, 32.45 and 30.1; MATLAB converts to radians
  more accurately than one double product in other places, so 8 % (``sind``) and 4 %
  (``cosd``) of dense arguments differ in the last bit (decision record 0005).
- ``sqrt`` is IEEE (exact everywhere); ``power(x, p)`` is ``x*x`` for ``p = 2``, ``sqrt(x)`` for
  0.5, ``1/x`` for -1, else fdlibm's pow.
- single data: the double function of the single value, rounded to single. MATLAB's single
  ``sin``, ``exp``, ... differ from that in the last bit for a few percent of the values (as
  glibc's ``sinf`` does); no Phase 1 code applies them to single data.

``fast_power`` is PetaKit5D's ``fastPower`` (binary powering in the class of the base);
``tukeywin`` MATLAB's Tukey window; ``padarray`` and ``circshift`` move elements only (exact).
"""

from __future__ import annotations

import numpy as np

from petakit5d.core import dtypes
from petakit5d.matlabcompat import neighbors as _nb

__all__ = [
    "atan2",
    "circshift",
    "cos",
    "cosd",
    "cot",
    "exp",
    "fast_power",
    "hypot",
    "padarray",
    "power",
    "sin",
    "sind",
    "sqrt",
    "tan",
    "tukeywin",
]

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

_hi, _lo, _with_hi, _from_hi = _nb._hi, _nb._lo, _nb._with_hi, _nb._from_hi


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
        with np.errstate(all="ignore"):
            r = fn(*[_f64(a) for a in arrs])
        r = r.astype(np.float32) if single else r
        return r[()] if r.ndim == 0 else r

    wrapped.__name__ = fn.__name__
    wrapped.__doc__ = fn.__doc__
    return wrapped


# -- sin, cos, tan, cot ------------------------------------------------------------------------


@_single
def sin(x):
    """MATLAB's ``sin`` (fdlibm)."""
    return _nb._sin(x)


@_single
def cos(x):
    """MATLAB's ``cos`` (fdlibm; the positive NaN for Inf and NaN, as MATLAB)."""
    return _pos_nan(x, _nb._cos(x))


_POS_NAN = float(np.array([0x7FF8000000000000], dtype=np.uint64).view(np.float64)[0])


def _pos_nan(x, r):
    """r with the positive quiet NaN where x is not finite (MATLAB's cos, tan and cot)"""
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
    n, y0, y1 = _nb._rem_pio2(x)
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


# -- sind, cosd --------------------------------------------------------------------------------


# the values MATLAB returns at multiples of 30 and 45 degrees (rounded exact values)
_DEG_EXACT = ((30.0, 0.5, 0.8660254037844386), (45.0, 0.7071067811865476, 0.7071067811865476))


def _deg(x, sine: bool):
    x = _f64(x)
    with np.errstate(all="ignore"):
        r = np.fmod(x, 360.0)  # exact (IEEE remainder of the division, the sign of x)
        n = np.asarray(dtypes.mround(r / 90.0))
        z = r - 90.0 * n
        rad = z * (np.pi / 180.0)
        s, c = _nb._sin(rad), _nb._cos(rad)
        az = np.abs(z)
        for a, sv, cv in _DEG_EXACT:
            hit = az == a
            s = np.where(hit, np.copysign(sv, z), s)
            c = np.where(hit, cv, c)
        q = np.mod(n, 4).astype(np.int64)
        if sine:
            out = np.select([q == 0, q == 1, q == 2], [s, c, -s], -c)
            out = np.where((z == 0) & (q % 2 == 0), np.copysign(0.0, x), out)
        else:
            out = np.select([q == 0, q == 1, q == 2], [c, -s, -c], s)
            out = np.where((z == 0) & (q % 2 == 1), 0.0, out)
        out = np.where(np.isfinite(x), out, dtypes.NAN if sine else -dtypes.NAN)
    return out


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
    return _nb._hypot(x, y)


@_single
def atan2(y, x):
    """MATLAB's ``atan2`` (fdlibm; +0 for two zeros)."""
    return _nb._atan2(y, x)


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
    x, y = np.broadcast_arrays(_f64(x), _f64(y))
    out = np.empty(x.shape)
    fx, fy, fo = x.reshape(-1), y.reshape(-1), out.reshape(-1)
    with np.errstate(all="ignore"):
        for i in range(fx.size):
            fo[i] = _pow_scalar(float(fx[i]), float(fy[i]))
    return out


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
        w = 0.5 - 0.5 * _nb._cos(2 * np.pi * t[:half])
        return np.concatenate([w, w[: n - half][::-1]])
    per = r / 2
    tl = int(np.floor(per * (n - 1))) + 1
    w = np.ones(n)
    w[:tl] = (1 + _nb._cos(np.pi / per * (t[:tl] - per))) / 2
    w[n - tl :] = (1 + _nb._cos(np.pi / per * (t[n - tl :] - 1 + per))) / 2
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
