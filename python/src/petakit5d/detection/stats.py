"""MATLAB's statistics functions as the point-detection code calls them (R2024b, Statistics and
Machine Learning Toolbox): ``norminv``, ``tcdf``, ``tinv``, ``chi2inv``, ``fcdf``, ``ksdensity``,
``ecdf``, ``mad`` and ``kstest2``.

Written from MATLAB's documentation and black-box probes (task C01); no MathWorks source was
read. MATLAB's special functions (``erfcinv``, ``betainc``, ``gammainc`` and their inverses) are
closed, so the distribution functions take SciPy's (Boost and Cephes, BSD) and differ from
MATLAB's in the last bits. Where they differ more, MATLAB is the less exact one: its ``tcdf`` and
``tinv`` are up to about 1e-12 off the exact values at the non-integer degrees of freedom of a
few hundred that detection's Welch tests have (SciPy's within a few ulp), its ``tcdf`` up to
1e-7 at 1e6 to 1e7 degrees of freedom, and its quantiles lose digits in the far tails (``tinv``
returns -Inf below p = 2^-54 for fewer than about 600 degrees of freedom, ``chi2inv`` NaN below
p = 1e-150 for hundreds of degrees of freedom); the port keeps the exact values there.
``docs/detection_scope.md`` (statistics) lists the measured differences and the proposed
tolerances. Exact where MATLAB's form is known:

- ``tcdf(x, 1)`` is ``atan2(1, -x)/pi`` below 0 and ``1 - atan2(1, x)/pi`` above, and
  ``tinv(p, 1)`` is ``tan(pi*(p - 0.5))``, with MATLAB's fdlibm ``atan2`` and ``tan``
  (``matlabcompat.elementary``): bit for bit on dense grids.
- MATLAB's edge values: ``p`` outside [0, 1], ``sigma <= 0``, ``nu <= 0`` and NaN give NaN;
  ``p`` 0 and 1 give -Inf and Inf (``chi2inv``: 0 and Inf); ``chi2inv(p, 0)`` is 0 inside (0, 1)
  and NaN at 0 and 1, ``chi2inv(p, Inf)`` NaN; ``tinv(p, Inf)`` is ``norminv(p)``; ``tcdf``
  with ``nu`` above 1e7 is the normal cdf ``erfc(-x/sqrt(2))/2``; ``fcdf`` with an infinite
  degree of freedom is the chi-square limit.

The sample statistics (``mad``, ``ecdf``, ``ksdensity``, ``kstest2``) are arithmetic in MATLAB's
order through ``petakit5d.matlabcompat``.

Arguments as MATLAB takes them: arrays of any shape or scalars; the non-scalar ones must have one
shape (MATLAB's ``InputSizeMismatch`` error otherwise); the result has that shape, or is a NumPy
scalar when every argument is one. Single data give single results (computed in double and
rounded; MATLAB computes in single). Data vectors are 1-D arrays.
"""

from __future__ import annotations

import numpy as np

from petakit5d.core import dtypes
from petakit5d.matlabcompat.reduce import MatlabError

__all__ = ["chi2inv", "ecdf", "fcdf", "kstest2", "ksdensity", "mad", "norminv", "tcdf", "tinv"]

NAN = dtypes.NAN
# tcdf: degrees of freedom above this give the normal cdf (MATLAB's own switch lies between 1e7
# and 3e7; the two forms agree within 2e-7 there)
_T_NORMAL_FROM = 1e7


def _prep(fn: str, *args):
    """the arguments as float64 arrays of the common shape, that shape (None: all scalars) and
    whether any was single"""
    arrs = []
    single = False
    for a in args:
        a = np.asarray(a)
        if a.dtype.kind not in "biuf":
            raise TypeError(f"{fn}: real numeric arguments required, got {a.dtype}")
        single = single or a.dtype == np.float32
        arrs.append(a.astype(np.float64))
    shape = None
    for a in arrs:
        if a.size != 1:
            if shape is None:
                shape = a.shape
            elif a.shape != shape:
                raise MatlabError(
                    f"stats:{fn}:InputSizeMismatch", "Non-scalar arguments must match in size."
                )
    tshape = () if shape is None else shape
    out = [np.broadcast_to(a.reshape(()) if a.size == 1 else a, tshape) for a in arrs]
    return out, shape, single


def _finish(r: np.ndarray, shape, single: bool):
    r = np.asarray(r, dtype=np.float64)
    if single:
        r = r.astype(np.float32)
    if shape is None:
        return r.reshape(())[()]
    return np.ascontiguousarray(r.reshape(shape))


def _set(r: np.ndarray, where: np.ndarray, value: float) -> None:
    r[np.broadcast_to(where, r.shape)] = value


# -- the normal, t, chi-square and F distributions --------------------------------------------


def _ndtri(p: np.ndarray) -> np.ndarray:
    from scipy import special

    with np.errstate(all="ignore"):
        return special.ndtri(p)


def norminv(p, mu=0.0, sigma=1.0):
    """MATLAB's ``norminv(p, mu, sigma)``: ``sigma * z + mu`` with ``z`` the standard normal
    quantile of ``p``; NaN for ``p`` outside [0, 1], ``sigma <= 0`` and NaN arguments."""
    (p, mu, sigma), shape, single = _prep("norminv", p, mu, sigma)
    with np.errstate(all="ignore"):
        x = np.array(sigma * _ndtri(p) + mu, dtype=np.float64)
    _set(x, ~((p >= 0) & (p <= 1) & (sigma > 0)) | np.isnan(mu), NAN)
    return _finish(x, shape, single)


def _tcdf_lower(x: np.ndarray, nu: np.ndarray) -> np.ndarray:
    from scipy import special

    from petakit5d.matlabcompat import elementary as E

    with np.errstate(all="ignore"):
        r = special.stdtr(nu, x)
        big = nu > _T_NORMAL_FROM
        if big.any():
            r = np.where(big, 0.5 * special.erfc(-x / np.sqrt(2.0)), r)
        one = nu == 1
        if one.any():
            ones = np.ones(x.shape)
            lo = E.atan2(ones, -x) / np.pi
            hi = 1.0 - E.atan2(ones, x) / np.pi
            r = np.where(one, np.where(x < 0, lo, hi), r)
    return r


def tcdf(x, nu, upper: bool = False):
    """MATLAB's ``tcdf(x, nu)`` (``upper=True``: ``tcdf(x, nu, 'upper')``), Student's t
    cumulative distribution; NaN for ``nu <= 0`` and NaN arguments."""
    (x, nu), shape, single = _prep("tcdf", x, nu)
    r = np.array(_tcdf_lower(-x if upper else x, nu), dtype=np.float64)
    _set(r, ~(nu > 0) | np.isnan(x), NAN)
    return _finish(r, shape, single)


def tinv(p, nu):
    """MATLAB's ``tinv(p, nu)``, the inverse of Student's t cdf; NaN for ``p`` outside [0, 1],
    ``nu <= 0`` and NaN arguments."""
    from scipy import special

    from petakit5d.matlabcompat import elementary as E

    (p, nu), shape, single = _prep("tinv", p, nu)
    with np.errstate(all="ignore"):
        # SciPy's inversion stops at a tolerance (1e-9 relative near p = 1/2): two Newton steps
        # on the t cdf make it exact to the cdf's accuracy
        r = special.stdtrit(nu, p)
        logc = special.gammaln((nu + 1) / 2) - special.gammaln(nu / 2) - 0.5 * np.log(nu * np.pi)
        for _ in range(2):
            pdf = np.exp(logc - (nu + 1) / 2 * np.log1p(r * r / nu))
            step = (special.stdtr(nu, r) - p) / pdf
            ok = np.isfinite(step) & (pdf > 0) & np.isfinite(r)
            r = np.where(ok, r - step, r)
        inf = np.isinf(nu)
        if inf.any():
            r = np.where(inf, _ndtri(p), r)
        one = nu == 1
        if one.any():
            r = np.where(one, E.tan(np.pi * (p - 0.5)), r)
    r = np.array(r, dtype=np.float64)
    _set(r, (p == 0) & (nu > 0), -np.inf)
    _set(r, (p == 1) & (nu > 0), np.inf)
    _set(r, ~((p >= 0) & (p <= 1) & (nu > 0)), NAN)
    return _finish(r, shape, single)


def chi2inv(p, nu):
    """MATLAB's ``chi2inv(p, nu)``, the inverse of the chi-square cdf; NaN for ``p`` outside
    [0, 1], ``nu < 0``, infinite ``nu`` and NaN arguments; ``nu`` 0 gives 0 inside (0, 1)."""
    from scipy import special

    (p, nu), shape, single = _prep("chi2inv", p, nu)
    with np.errstate(all="ignore"):
        r = 2.0 * special.gammaincinv(nu / 2.0, p)
    r = np.array(r, dtype=np.float64)
    _set(r, (p == 0) & (nu > 0), 0.0)
    _set(r, (p == 1) & (nu > 0), np.inf)
    _set(r, (nu == 0) & (p > 0) & (p < 1), 0.0)
    _set(r, ~((p >= 0) & (p <= 1) & (nu >= 0) & np.isfinite(nu)), NAN)
    _set(r, (nu == 0) & ((p == 0) | (p == 1)), NAN)
    return _finish(r, shape, single)


def fcdf(x, nu1, nu2, upper: bool = False):
    """MATLAB's ``fcdf(x, nu1, nu2)`` (``upper=True``: ``'upper'``), the F cumulative
    distribution; NaN for ``nu1 <= 0``, ``nu2 <= 0`` and NaN arguments."""
    from scipy import special

    (x, nu1, nu2), shape, single = _prep("fcdf", x, nu1, nu2)
    with np.errstate(all="ignore"):
        xp = np.maximum(x, 0.0)
        lo = special.fdtr(nu1, nu2, xp)
        hi = special.fdtrc(nu1, nu2, xp)
        i1, i2 = np.isinf(nu1) & np.isfinite(nu2), np.isinf(nu2) & np.isfinite(nu1)
        if i1.any():  # F -> nu2 / chi2(nu2)
            lo = np.where(i1, special.chdtrc(nu2, nu2 / xp), lo)
            hi = np.where(i1, special.chdtr(nu2, nu2 / xp), hi)
        if i2.any():  # F -> chi2(nu1) / nu1
            lo = np.where(i2, special.chdtr(nu1, nu1 * xp), lo)
            hi = np.where(i2, special.chdtrc(nu1, nu1 * xp), hi)
        r = np.array(hi if upper else lo, dtype=np.float64)
    _set(r, x <= 0, 1.0 if upper else 0.0)
    _set(r, np.isposinf(x), 0.0 if upper else 1.0)
    _set(r, ~((nu1 > 0) & (nu2 > 0)) | np.isnan(x), NAN)
    return _finish(r, shape, single)


# -- sample statistics -------------------------------------------------------------------------


def mad(x, flag: int = 0, axis=None, *, matlab_layout: bool = False):
    """MATLAB's ``mad(x, flag)``: the mean absolute deviation ``mean(abs(x - mean(x)))``
    (``flag`` 0) or the median absolute deviation ``median(abs(x - median(x)))`` (``flag`` 1),
    NaN values left out. ``axis`` None: all values of a 1-D array, else MATLAB's default
    dimension; as in ``matlabcompat.reduce`` otherwise."""
    from petakit5d.matlabcompat import reduce as R
    from petakit5d.matlabcompat import stats as S

    a = np.asarray(x)
    if a.dtype.kind != "f":
        a = a.astype(np.float64)
    vector = a.ndim <= 1
    if axis is None:
        axis = -1 if vector else R.default_axis(a, matlab_layout=matlab_layout)
    if flag == 0:
        c = R.mean(a, axis, omitnan=True, keepdims=True, matlab_layout=matlab_layout)
        r = R.mean(np.abs(a - c), axis, omitnan=True, matlab_layout=matlab_layout)
    else:
        c = S.median(a, axis, omitnan=True, keepdims=True, matlab_layout=matlab_layout)
        r = S.median(np.abs(a - c), axis, omitnan=True, matlab_layout=matlab_layout)
    return np.asarray(r).reshape(())[()] if vector else r


def ecdf(y):
    """MATLAB's ``[f, x] = ecdf(y)`` for uncensored data: ``x`` the sorted distinct values of
    ``y`` after its smallest one (which comes twice), ``f`` 0 then the Kaplan-Meier estimate at
    each distinct value, ``1 - cumprod(1 - d./r)`` (``d`` the count of the value, ``r`` the
    count of values not below it; probed: not ``cumsum(d)/n``); NaN values left out. Both are
    1-D (MATLAB's columns)."""
    v = np.asarray(y, dtype=np.float64).reshape(-1, order="F")
    v = v[~np.isnan(v)]
    n = v.size
    s = np.sort(v, kind="stable")
    last = np.append(s[1:] != s[:-1], True) if n else np.zeros(0, bool)
    xs = s[last]
    ends = np.flatnonzero(last) + 1
    d = np.diff(np.concatenate(([0], ends))).astype(np.float64)
    r = (n - np.concatenate(([0], ends[:-1]))).astype(np.float64)
    f = 1.0 - np.cumprod(1.0 - d / r)
    return np.concatenate(([0.0], f)), np.concatenate((xs[:1], xs))


def ksdensity(x, xi=None, *, bandwidth=None, npoints: int = 100):
    """MATLAB's ``[f, xi, bw] = ksdensity(x)`` and ``ksdensity(x, xi)``: the normal-kernel
    density estimate, at ``npoints`` points over the data's range and 3 bandwidths each side, or
    at ``xi``. The bandwidth is MATLAB's default, ``sig * (4/(3n))^(1/5)`` with ``sig`` the
    median absolute deviation over 0.6745 (bit for bit, with MATLAB's ``pow``); the density the
    direct kernel sum, MATLAB's within a few ulp for samples of a few hundred values. For larger
    samples MATLAB evaluates otherwise (1.5% lower in the far tails at 500 values; not probed
    further: the traced code does not call ksdensity)."""
    from petakit5d.matlabcompat import elementary as E
    from petakit5d.matlabcompat import reduce as R
    from petakit5d.matlabcompat import stats as S

    v = np.asarray(x, dtype=np.float64).reshape(-1, order="F")
    v = v[~np.isnan(v)]
    n = v.size
    if bandwidth is None:
        sig = float(S.median(np.abs(v - S.median(v)))) / 0.6745
        bw = sig * float(np.asarray(E.power(np.float64(4.0 / (3.0 * n)), 0.2)).reshape(-1)[0])
    else:
        bw = float(bandwidth)
    if xi is None:
        # the data's range and 3 bandwidths each side, k*(hi - lo)/(n - 1) steps (probed)
        lo, hi = float(np.min(v)) - 3 * bw, float(np.max(v)) + 3 * bw
        xi = lo + np.arange(npoints, dtype=np.float64) * (hi - lo) / (npoints - 1)
    xq = np.asarray(xi, dtype=np.float64)
    z = (xq.reshape(-1, 1) - v.reshape(1, -1)) / bw
    k = E.exp(-0.5 * z * z) / np.sqrt(2.0 * np.pi)
    f = R.sum(k, 1, matlab_layout=True) / (n * bw)
    return f.reshape(xq.shape), xq, bw


def kstest2(x1, x2, *, alpha: float = 0.05, tail: str = "unequal"):
    """MATLAB's ``[h, p, ks2stat] = kstest2(x1, x2)``: the two-sample Kolmogorov-Smirnov test
    with the asymptotic p-value; NaN values left out."""
    from petakit5d.matlabcompat import elementary as E
    from petakit5d.matlabcompat import reduce as R

    a = np.asarray(x1, dtype=np.float64).reshape(-1, order="F")
    b = np.asarray(x2, dtype=np.float64).reshape(-1, order="F")
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    n1, n2 = a.size, b.size
    edges = np.sort(np.concatenate((a, b)), kind="stable")
    c1 = np.searchsorted(np.sort(a), edges, side="right") / n1
    c2 = np.searchsorted(np.sort(b), edges, side="right") / n2
    if tail == "unequal":
        d = np.abs(c1 - c2)
    elif tail == "larger":  # the alternative: x1's cdf is larger than x2's
        d = c1 - c2
    elif tail == "smaller":
        d = c2 - c1
    else:
        raise ValueError(f"kstest2: unknown tail {tail!r}")
    ks = float(np.max(d))
    n = n1 * n2 / (n1 + n2)
    lam = max((np.sqrt(n) + 0.12 + 0.11 / np.sqrt(n)) * ks, 0.0)
    if tail == "unequal":
        j = np.arange(1, 102, dtype=np.float64)
        terms = (-1.0) ** (j - 1) * E.exp(-2.0 * lam * lam * j * j)
        p = 2.0 * float(np.asarray(R.sum(terms)).reshape(-1)[0])
        p = min(max(p, 0.0), 1.0)
    else:
        p = float(E.exp(-2.0 * lam * lam))
    return bool(p < alpha), p, ks
