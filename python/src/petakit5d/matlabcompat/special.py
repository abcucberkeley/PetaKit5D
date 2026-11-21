"""MATLAB's ``erf`` and ``erfc`` of double data, bit for bit on every platform (R2024b; measured
on 409,402 values from 1e-300 to the tails, task C05), and ``normcdf`` as the lifetime analysis
calls it.

MATLAB's ``erf`` and ``erfc`` are fdlibm's (``s_erf.c``) with fdlibm's ``exp``; SciPy's differ
from them in the last bit for 7 % of the values (``erf``) and by up to whole subnormal ranges in
the far tail of ``erfc``. They are written out here in NumPy, one IEEE operation per C operation
(no fused multiply-adds), with ``matlabcompat.elementary.exp`` (fdlibm's; ``exp(1)`` never occurs
here). NaN stays NaN (its bits as given), ``erf(+-Inf)`` is +-1, ``erfc(+Inf)`` 0, ``erfc(-Inf)``
2, as in MATLAB. Other classes are converted to double first.

The algorithms and constants follow fdlibm 5.3's ``s_erf.c``, whose notice is:
  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
  Developed at SunSoft, a Sun Microsystems, Inc. business.
  Permission to use, copy, modify, and distribute this software is freely granted, provided
  that this notice is preserved.
"""

from __future__ import annotations

import numpy as np

from petakit5d.matlabcompat import elementary as _el

__all__ = ["erf", "erfc", "normcdf"]

_ERX = 8.45062911510467529297e-01
_EFX = 1.28379167095512586316e-01
_EFX8 = 1.02703333676410069053e00
_PP = (
    1.28379167095512558561e-01,
    -3.25042107247001499370e-01,
    -2.84817495755985104766e-02,
    -5.77027029648944159157e-03,
    -2.37630166566501626084e-05,
)
_QQ = (
    3.97917223959155352819e-01,
    6.50222499887672944485e-02,
    5.08130628187576562776e-03,
    1.32494738004321644526e-04,
    -3.96022827877536812320e-06,
)
_PA = (
    -2.36211856075265944077e-03,
    4.14856118683748331666e-01,
    -3.72207876035701323847e-01,
    3.18346619901161753674e-01,
    -1.10894694282396677476e-01,
    3.54783043256182359371e-02,
    -2.16637559486879084300e-03,
)
_QA = (
    1.06420880400844228286e-01,
    5.40397917702171048937e-01,
    7.18286544141962662868e-02,
    1.26171219808761642112e-01,
    1.36370839120290507362e-02,
    1.19844998467991074170e-02,
)
_RA = (
    -9.86494403484714822705e-03,
    -6.93858572707181764372e-01,
    -1.05586262253232909814e01,
    -6.23753324503260060396e01,
    -1.62396669462573470355e02,
    -1.84605092906711035994e02,
    -8.12874355063065934246e01,
    -9.81432934416914548592e00,
)
_SA = (
    1.96512716674392571292e01,
    1.37657754143519042600e02,
    4.34565877475229228821e02,
    6.45387271733267880336e02,
    4.29008140027567833386e02,
    1.08635005541779435134e02,
    6.57024977031928170135e00,
    -6.04244152148580987438e-02,
)
_RB = (
    -9.86494292470009928597e-03,
    -7.99283237680523006574e-01,
    -1.77579549177547519889e01,
    -1.60636384855821916062e02,
    -6.37566443368389627722e02,
    -1.02509513161107724954e03,
    -4.83519191608651397019e02,
)
_SB = (
    3.03380607434824582924e01,
    3.25792512996573918826e02,
    1.53672958608443695994e03,
    3.19985821950859553908e03,
    2.55305040643316442583e03,
    4.74528541206955367215e02,
    -2.24409524465858183362e01,
)
_TINY = 1e-300


def _words(x: np.ndarray):
    u = x.view(np.uint64)
    hx = (u >> np.uint64(32)).astype(np.uint32).view(np.int32).astype(np.int64)
    return hx, hx & 0x7FFFFFFF, (u >> np.uint64(63)).astype(bool)


def _poly(c, s, first=None):
    """c[0] + s*(c[1] + s*(...)), or 1 + s*(c[0] + ...) with first=1 (as fdlibm writes it)"""
    r = c[-1]
    for k in c[-2::-1]:
        r = k + s * r
    return r if first is None else first + s * r


def _small(x):
    """erf of |x| < 0.84375: x + x*y"""
    z = x * x
    y = _poly(_PP, z) / _poly(_QQ, z, 1.0)
    return y


def _mid(ax):
    """the ratio P/Q of 0.84375 <= |x| < 1.25"""
    s = ax - 1.0
    return _poly(_PA, s) / _poly(_QA, s, 1.0)


def _tail(ax, ix, bound):
    """exp(-z*z - 0.5625)*exp((z - x)*(z + x) + R/S), 1.25 <= |x|; bound: the switch from the ra/sa
    to the rb/sb rational function (erf and erfc differ in its last bit)"""
    s = 1.0 / (ax * ax)
    use_a = ix < bound
    R = np.where(use_a, _poly(_RA, s), _poly(_RB, s))
    S = np.where(use_a, _poly(_SA, s, 1.0), _poly(_SB, s, 1.0))
    z = (ax.view(np.uint64) & np.uint64(0xFFFFFFFF00000000)).view(np.float64)
    return _el.exp(-z * z - 0.5625) * _el.exp((z - ax) * (z + ax) + R / S)


def _as_double(x) -> np.ndarray:
    a = np.asarray(x)
    if a.dtype != np.float64:
        a = a.astype(np.float64)
    return np.ascontiguousarray(a)


def erf(x):
    """MATLAB's ``erf`` (module docstring)."""
    x = _as_double(x)
    shape = x.shape
    x = x.reshape(-1)
    out = np.empty_like(x)
    hx, ix, neg = _words(x)
    with np.errstate(all="ignore"):
        m = ix >= 0x7FF00000
        out[m] = np.where(np.isnan(x[m]), x[m], np.where(neg[m], -1.0, 1.0))
        m1 = ix < 0x3FEB0000
        if m1.any():
            xx = x[m1]
            i1 = ix[m1]
            v = xx + xx * _small(xx)
            v = np.where(
                i1 < 0x3E300000,
                np.where(i1 < 0x00800000, 0.125 * (8.0 * xx + _EFX8 * xx), xx + _EFX * xx),
                v,
            )
            out[m1] = v
        m2 = (ix >= 0x3FEB0000) & (ix < 0x3FF40000)
        if m2.any():
            pq = _mid(np.abs(x[m2]))
            out[m2] = np.where(neg[m2], -_ERX - pq, _ERX + pq)
        m3 = (ix >= 0x40180000) & ~m
        out[m3] = np.where(neg[m3], _TINY - 1.0, 1.0 - _TINY)
        m4 = (ix >= 0x3FF40000) & (ix < 0x40180000)
        if m4.any():
            ax = np.abs(x[m4])
            r = _tail(ax, ix[m4], 0x4006DB6E)
            out[m4] = np.where(neg[m4], r / ax - 1.0, 1.0 - r / ax)
    return out.reshape(shape)


def erfc(x):
    """MATLAB's ``erfc`` (module docstring)."""
    x = _as_double(x)
    shape = x.shape
    x = x.reshape(-1)
    out = np.empty_like(x)
    hx, ix, neg = _words(x)
    with np.errstate(all="ignore"):
        m = ix >= 0x7FF00000
        out[m] = np.where(np.isnan(x[m]), x[m], np.where(neg[m], 2.0, 0.0))
        m1 = ix < 0x3FEB0000
        if m1.any():
            xx = x[m1]
            y = _small(xx)
            r = xx * y
            v = np.where(hx[m1] < 0x3FD00000, 1.0 - (xx + xx * y), 0.5 - (r + (xx - 0.5)))
            out[m1] = np.where(ix[m1] < 0x3C700000, 1.0 - xx, v)
        m2 = (ix >= 0x3FEB0000) & (ix < 0x3FF40000)
        if m2.any():
            pq = _mid(np.abs(x[m2]))
            out[m2] = np.where(neg[m2], 1.0 + (_ERX + pq), (1.0 - _ERX) - pq)
        m3 = (ix >= 0x3FF40000) & (ix < 0x403C0000)
        if m3.any():
            ax = np.abs(x[m3])
            ng = neg[m3]
            r = _tail(ax, ix[m3], 0x4006DB6D)
            v = np.where(ng, 2.0 - r / ax, r / ax)
            out[m3] = np.where(ng & (ix[m3] >= 0x40180000), 2.0 - _TINY, v)
        m4 = (ix >= 0x403C0000) & ~m
        out[m4] = np.where(neg[m4], 2.0 - _TINY, _TINY * _TINY)
    return out.reshape(shape)


_SQRT2 = float(np.sqrt(2.0))


def normcdf(x, mu=0.0, sigma=1.0):
    """MATLAB's ``normcdf(x, mu, sigma)``: ``0.5*erfc(-((x - mu)/sigma)/sqrt(2))``, bit for bit
    (measured on 40,002 values with sigma 1 and a sigma of 1.62; task C05). Arguments broadcast;
    sigma > 0 as the lifetime analysis passes it (other sigma not probed)."""
    x, mu, sigma = (_as_double(v) for v in (x, mu, sigma))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (x - mu) / sigma
    return 0.5 * erfc(-z / _SQRT2)
