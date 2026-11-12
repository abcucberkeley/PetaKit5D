"""MATLAB's ``knnsearch`` (Euclidean, as the OMW back projector and the group shift assignment call
it) and ``cart2sph`` / ``sph2cart``.

``knnsearch(X, Y, k)`` gives, for each row of ``Y``, the ``k`` rows of ``X`` nearest to it and their
distances, as MATLAB R2024b does (its default k-d tree and its exhaustive search agree): the
squared distance is summed over the columns in order, ``((y1 - x1)^2 + (y2 - x2)^2) + ...``,
neighbors are ordered by it and, among equal squared distances, by row of ``X``; the distance is
its square root. (Ordering by the rounded distance instead would differ: two squared distances
can have the same root.) Indices are 0-based. Small problems are searched exhaustively; larger
ones take candidates from a k-d tree (every point within the k-th distance and a margin) and rank
them the same way, so the result is the same.

``cart2sph`` and ``sph2cart`` are MATLAB's formulas (``hypot(hypot(x, y), z)``,
``atan2(z, hypot(x, y))``, ``atan2(y, x)``; ``r .* sin(elev)``, ``r .* cos(elev) .* cos(az)``, ...).
``hypot``, ``atan2``, ``sin`` and ``cos`` are written out here (fdlibm's, as MATLAB computes them),
so that the results are MATLAB's on every C library and CPU: ``hypot`` with MATLAB's NaN for a NaN
argument, also next to Inf, and ``atan2`` of two zeros +0 as in MATLAB (C gives -0 or +-pi for
signed zeros). The C libraries' ``hypot`` (glibc 2.35 and later), ``atan2``, ``sin`` and ``cos``
differ from MATLAB's in the last bit for a part of the inputs (MATLAB's ``atan2`` is also 1 ulp
above +-pi/2 for a negative x tiny next to y); the coordinates of ``sph2cart`` would carry such
differences of the angles, amplified where a cosine or sine is near zero (up to about 5e-16 of the
radius). See ``_hypot``, ``_atan2``, ``_sin``, ``_cos``.
"""

from __future__ import annotations

import numpy as np

__all__ = ["knnsearch", "cart2sph", "sph2cart", "KnnError"]

# exhaustive search up to this many distances per block of queries
_BLOCK = 1 << 22
# above this many distances in all, a k-d tree gives the candidates
_EXHAUSTIVE_MAX = 1 << 26


class KnnError(ValueError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def knnsearch(X, Y, k: int = 1):
    """``[idx, d] = knnsearch(X, Y, 'K', k)``: (idx (ny, min(k, nx)) int64, 0-based; d float64).
    Rows are points; a 1-D array is one point (a row). Rows of ``X`` with a NaN are never
    neighbors; a query with a NaN gets the first rows of ``X`` and NaN distances (MATLAB)."""
    from petakit5d.core.dtypes import NAN

    X = _points(X, "X")
    Y = _points(Y, "Y")
    if X.shape[1] != Y.shape[1]:
        raise KnnError(
            "X and Y must have the same number of columns",
            "stats:KDTreeSearcher:knnsearch:SizeMisMatch",
        )
    if not float(k) == int(k) or int(k) < 1:
        raise KnnError("K must be a positive integer", "stats:KDTreeSearcher:knnsearch:BadK")
    k = int(k)
    valid = ~np.any(np.isnan(X), axis=1)
    rowmap = np.nonzero(valid)[0]
    Xv = X[valid] if not valid.all() else X
    nx, ny = Xv.shape[0], Y.shape[0]
    kk = min(k, nx)
    idx = np.zeros((ny, kk), np.int64)
    dist = np.zeros((ny, kk))
    if ny == 0 or kk == 0:
        return idx, dist
    qok = ~np.any(np.isnan(Y), axis=1)
    Yq = Y[qok]
    if Yq.shape[0]:
        if nx * Yq.shape[0] <= _EXHAUSTIVE_MAX or X.shape[1] > 16:
            i, d = _exhaustive(Xv, Yq, kk)
        else:
            i, d = _tree(Xv, Yq, kk)
        idx[qok] = rowmap[i]
        dist[qok] = d
    idx[~qok] = rowmap[:kk]
    dist[~qok] = NAN
    return idx, dist


def _points(A, name):
    A = np.asarray(A, dtype=np.float64)
    if A.ndim == 1:
        A = A[None, :]
    if A.ndim != 2:
        raise KnnError(f"{name} must be a matrix")
    return A


def _sqdist(Y, X):
    """squared distances (len(Y), len(X)), summed over the columns in order"""
    d2 = np.zeros((Y.shape[0], X.shape[0]))
    for j in range(X.shape[1]):
        df = Y[:, j][:, None] - X[:, j][None, :]
        d2 += df * df
    return d2


def _select(d2, cand, k):
    """the k smallest of each row of d2 by (value, candidate index): positions into the row"""
    # sort by candidate index, then stably by squared distance
    o = np.argsort(cand, axis=1, kind="stable")
    d2o = np.take_along_axis(d2, o, axis=1)
    o2 = np.argsort(d2o, axis=1, kind="stable")[:, :k]
    return np.take_along_axis(o, o2, axis=1)


def _exhaustive(X, Y, k):
    nx = X.shape[0]
    idx = np.empty((Y.shape[0], k), np.int64)
    d2k = np.empty((Y.shape[0], k))
    step = max(1, _BLOCK // max(1, nx))
    for a in range(0, Y.shape[0], step):
        d2 = _sqdist(Y[a : a + step], X)
        rows = np.arange(d2.shape[0])
        if k <= 4 and np.all(np.isfinite(d2)):
            # the first minimum is the lowest index among equal values
            for j in range(k):
                o = np.argmin(d2, axis=1)
                idx[a + rows, j] = o
                d2k[a + rows, j] = d2[rows, o]
                d2[rows, o] = np.inf
        else:
            # columns are in index order: a stable sort keeps index order among equal values
            o = np.argsort(d2, axis=1, kind="stable")[:, :k]
            idx[a : a + step] = o
            d2k[a : a + step] = np.take_along_axis(d2, o, axis=1)
    return idx, np.sqrt(d2k)


def _tree(X, Y, k):
    """candidates from a k-d tree: every point within the k-th distance (and a margin), ranked
    exactly"""
    from scipy.spatial import cKDTree

    tree = cKDTree(X)
    extra = 8
    m = min(X.shape[0], k + extra)
    dist, cand = tree.query(Y, k=m)
    dist = np.asarray(dist).reshape(Y.shape[0], m)
    cand = np.asarray(cand).reshape(Y.shape[0], m)
    radius = dist[:, k - 1] * (1 + 1e-9) + 1e-300
    complete = (m == X.shape[0]) | (dist[:, m - 1] > radius)
    idx = np.empty((Y.shape[0], k), np.int64)
    d2k = np.empty((Y.shape[0], k))
    rows = np.nonzero(complete)[0]
    if rows.size:
        c = cand[rows]
        d2 = np.zeros(c.shape)
        for j in range(X.shape[1]):
            df = Y[rows, j][:, None] - X[c, j]
            d2 += df * df
        pos = _select(d2, c, k)
        idx[rows] = np.take_along_axis(c, pos, axis=1)
        d2k[rows] = np.take_along_axis(d2, pos, axis=1)
    for r in np.nonzero(~complete)[0]:
        c = np.asarray(sorted(tree.query_ball_point(Y[r], radius[r])), np.int64)
        d2 = _sqdist(Y[r : r + 1], X[c])
        pos = _select(d2, c[None, :], k)[0]
        idx[r] = c[pos]
        d2k[r] = d2[0, pos]
    return idx, np.sqrt(d2k)


# --- the elementary functions of cart2sph and sph2cart ------------------------------------------
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


def cart2sph(x, y, z):
    """``[az, elev, r] = cart2sph(x, y, z)`` (float64)"""
    x, y, z = (np.asarray(v, dtype=np.float64) for v in (x, y, z))
    hxy = _hypot(x, y)
    r = _hypot(hxy, z)
    elev = _atan2(z, hxy)
    az = _atan2(y, x)
    return az, elev, r


def sph2cart(az, elev, r):
    """``[x, y, z] = sph2cart(az, elev, r)`` (float64)"""
    az, elev, r = (np.asarray(v, dtype=np.float64) for v in (az, elev, r))
    z = _mul(r, _sin(elev))
    rcoselev = _mul(r, _cos(elev))
    x = _mul(rcoselev, _cos(az))
    y = _mul(rcoselev, _sin(az))
    return x, y, z


def _mul(a, b):
    """a * b with a NaN operand's payload kept as x86 does (the first NaN operand); NumPy may
    swap the operands of a product"""
    with np.errstate(invalid="ignore"):
        p = a * b
    return np.where(np.isnan(a), a, np.where(np.isnan(b), b, p))
