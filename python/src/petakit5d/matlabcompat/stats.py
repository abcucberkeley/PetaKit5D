"""MATLAB's var, std, median, prctile and mode, exactly (R2024b; measured, N05;
``docs/decisions/0005-sums-and-math-library.md``).

Arrays and axes as in ``matlabcompat.reduce``: the port's ``(z, y, x)`` C arrays, one NumPy axis
or None for all elements in MATLAB's order (``x(:)``), ``matlab_layout=True`` for arrays in
MATLAB's own form.

- ``var(x, w)``: ``sum(abs(x - mean(x)).^2) / (n - 1)`` (``w = 1``: ``/ n``), the mean and the sum
  in ``reduce``'s order, in the class (single data in single; the division in double, rounded);
  ``abs`` clears the sign of a NaN, as in MATLAB. One value has variance 0, none NaN; integer
  data raise ``MATLAB:var:integerClass``. ``std`` is its square root.
- ``median``: the middle element of the stably sorted values (so the -0 or +0 MATLAB picks), for an
  even count the mean of the middle two, ``a + (b - a)/2``, or ``(a + b)/2`` when their signs
  differ or one is infinite (integers: ``(a + b)/2`` in double, rounded half away from zero).
  A NaN gives NaN (``omitnan``: NaN left out); no values give NaN.
- ``prctile`` (MATLAB's default 'exact' method): NaN left out; with the sorted values
  ``x(1..n)``, ``q = p/100*n``, ``k = floor(q + 0.5)``, ``r = q - (k - 0.5)``: ``x(1)`` for
  ``k < 1``, ``x(n)`` for ``k >= n``, ``x(k)`` for ``r = 0`` or ``x(k) = x(k+1)``, else
  ``(1 - r)*x(k) + r*x(k+1)`` (single data: each product rounded to single, the sum in single;
  integers in double, rounded half away from zero). No values give NaN.
- ``mode``: the smallest of the most frequent values (each NaN counted alone); with one or two
  outputs MATLAB returns the last of those equal values in its stable sort (``mode([0 -0])`` is
  -0), with three the first (``first_of_run=True``); ``with_count`` adds the count. No values
  give NaN and count 0.

Along an axis the values are sorted stably (``kind='stable'``, a true stable sort on every
platform); over all elements the k-th values come from ``numpy.partition`` and the -0 or +0 that
MATLAB's stable sort puts there from counting in MATLAB's order, so no result depends on NumPy's
sorting algorithms or SIMD paths.
"""

from __future__ import annotations

import builtins

import numpy as np

from petakit5d.core import dtypes
from petakit5d.matlabcompat import reduce as _r

__all__ = ["median", "mode", "prctile", "std", "var"]


def _nan_of(dt: np.dtype):
    return (np.float64 if dt.kind in "biu" else dt.type)(dtypes.NAN)


# -- var and std -------------------------------------------------------------------------------


def var(
    x,
    w: int = 0,
    axis=None,
    *,
    omitnan: bool = False,
    keepdims: bool = False,
    matlab_layout: bool = False,
    num_threads: int | None = None,
):
    """MATLAB's ``var(x, w, dim, nanflag)`` for ``w`` 0 (``/(n - 1)``) or 1 (``/n``)."""
    a = _r._array(x)
    if a.dtype.kind in "iu":
        raise _r.MatlabError(
            "MATLAB:var:integerClass", "Invalid data type. First argument must be single or double."
        )
    if w not in (0, 1):
        raise ValueError("w must be 0 or 1 (weight vectors are not supported)")
    if a.dtype.kind == "b":
        a = a.astype(np.float64)
    dims = _r._mdims(axis, a.ndim, matlab_layout)
    m = _r._mview(a, matlab_layout)
    nt = _r.resolve_threads(num_threads)
    dt = a.dtype
    q, _ = _r._mean_view(a, m, dims, "default", omitnan, nt, matlab_layout)
    omit = bool(omitnan)
    s2 = None
    if not omit:
        # K08's fused sum of abs(x - mean).^2 (the same values in the same order); a NaN result
        # goes the NumPy way below, which gives MATLAB's payload
        s2 = _r._kernel_sum(a, m, dims, dt, nt, matlab_layout, center=q.astype(dt))
        if s2 is not None and np.isnan(s2).any():
            s2 = None
    if s2 is None:
        d = np.abs(np.subtract(m, q.astype(dt), dtype=dt))
        sq = np.multiply(d, d, dtype=dt)
        s2 = _r._sum_core(sq, sq, dims, dt, omit, nt, True)
    n = np.asarray(_r._count(m, dims, omit), dtype=np.float64)
    den = np.where(n > 1, n - 1, 1.0) if w == 0 else np.where(n > 0, n, 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = (s2.astype(np.float64) / den).astype(dt)
    v = np.where(n == 0, dt.type(dtypes.NAN), v).astype(dt)
    return _r._finish(v, a, dims, matlab_layout, keepdims)


def std(
    x,
    w: int = 0,
    axis=None,
    *,
    omitnan: bool = False,
    keepdims: bool = False,
    matlab_layout: bool = False,
    num_threads: int | None = None,
):
    """MATLAB's ``std(x, w, dim, nanflag)``: ``sqrt(var(...))`` in the class (IEEE square root)."""
    v = var(
        x,
        w,
        axis,
        omitnan=omitnan,
        keepdims=keepdims,
        matlab_layout=matlab_layout,
        num_threads=num_threads,
    )
    return np.sqrt(v)


# -- selection ---------------------------------------------------------------------------------


def _axis(a: np.ndarray, axis) -> int | None:
    if axis is None:
        return None
    if not isinstance(axis, int | np.integer):
        raise TypeError("give one axis, or None for all elements")
    if a.ndim == 0:
        return 0
    return int(axis) % a.ndim


def _values(a: np.ndarray, omitnan: bool):
    """(the values to select from, as a 1-D copy in any order; whether a NaN decides)"""
    v = np.array(a, copy=True).reshape(-1)
    if v.dtype.kind == "f":
        nan = np.isnan(v)
        if nan.any():
            if not omitnan:
                return None, True
            v = v[~nan]
    return v, False


def _kth_whole(a: np.ndarray, v: np.ndarray, ks, matlab_layout: bool) -> list:
    """the elements at 0-based positions ks of MATLAB's stable sort of all of a (no NaN in v,
    v changed)"""
    v.partition(sorted(set(ks)))
    out = []
    for k in ks:
        val = v[k]
        if val == 0 and v.dtype.kind == "f":
            # MATLAB's stable sort: the (k - #negatives)-th zero in MATLAB's order
            j = k - int(np.count_nonzero(a < 0))
            t = _r._nth_true_matlab(a == 0, j, matlab_layout)
            val = a[_r._index_tuple(t, a, matlab_layout)]
        out.append(val)
    return out


def _meanof(a, b, dt):
    """the mean of a <= b as MATLAB's median takes it (arrays)"""
    if dt.kind in "iu":
        return np.asarray(dtypes.cast((a.astype(np.float64) + b) / 2, dt.name))
    T = dt.type
    with np.errstate(invalid="ignore", over="ignore"):
        near = np.add(a, np.subtract(b, a, dtype=dt) / T(2), dtype=dt)
        far = np.add(a, b, dtype=dt) / T(2)
    use_far = (np.sign(a) != np.sign(b)) | np.isinf(a) | np.isinf(b)
    r = np.where(use_far, far, near).astype(dt)
    return np.where(np.isnan(r), T(dtypes.NAN), r).astype(dt)


def _sorted_along(a: np.ndarray, ax: int):
    """(the values sorted stably along ax as the last axis, NaN last; the count of non-NaN)"""
    w = np.moveaxis(a, ax, -1)
    s = np.sort(w, axis=-1, kind="stable")
    if a.dtype.kind == "f":
        nv = np.count_nonzero(~np.isnan(w), axis=-1)
    else:
        nv = np.full(w.shape[:-1], w.shape[-1], dtype=np.int64)
    return s, nv


def _take(s: np.ndarray, k: np.ndarray) -> np.ndarray:
    k = np.clip(k, 0, builtins.max(s.shape[-1] - 1, 0))
    return np.take_along_axis(s, k[..., None].astype(np.intp), axis=-1)[..., 0]


def median(
    x, axis=None, *, omitnan: bool = False, keepdims: bool = False, matlab_layout: bool = False
):
    """MATLAB's ``median(x, dim, nanflag)`` (module docstring); integers keep their class."""
    a = _r._array(x)
    if a.dtype.kind == "b":
        a = a.astype(np.float64)
    dt = a.dtype
    ax = _axis(a, axis)
    if ax is None:
        v, nan = _values(a, omitnan)
        n = 0 if v is None else v.size
        if nan or n == 0:
            r = _nan_of(dt)
        elif n % 2:
            r = _kth_whole(a, v, [n // 2], matlab_layout)[0]
        else:
            lo, hi = _kth_whole(a, v, [n // 2 - 1, n // 2], matlab_layout)
            r = _meanof(np.asarray(lo), np.asarray(hi), dt)[()]
        if keepdims:
            return np.full((1,) * a.ndim, r)
        return r
    if a.ndim == 0:
        a = a.reshape(1)
    s, nv = _sorted_along(a, ax)
    if s.shape[-1] == 0:
        r = np.full(nv.shape, _nan_of(dt))
    else:
        lo, hi = _take(s, (nv - 1) // 2), _take(s, nv // 2)
        r = np.where(nv % 2 == 1, lo, _meanof(lo, hi, dt)).astype(dt)
        if dt.kind == "f":
            bad = (nv == 0) | ((nv < s.shape[-1]) & (not omitnan))
            r = np.where(bad, dt.type(dtypes.NAN), r).astype(dt)
    if keepdims:
        r = np.expand_dims(r, ax)
    return r[()] if r.ndim == 0 else r


# -- prctile -----------------------------------------------------------------------------------


def _prc_plan(n, p: np.ndarray):
    """k (1-based lower neighbor) and r for each p, for n values (arrays broadcast)"""
    q = p / 100 * np.asarray(n, dtype=np.float64)
    k = np.floor(q + 0.5)
    return k.astype(np.int64), q - (k - 0.5)


def _interp(lo: np.ndarray, hi: np.ndarray, r: np.ndarray, dt) -> np.ndarray:
    with np.errstate(over="ignore", invalid="ignore"):
        if dt == np.float32:
            f4 = np.float32
            t1 = ((1 - r) * lo.astype(np.float64)).astype(f4)
            t2 = (r * hi.astype(np.float64)).astype(f4)
            v = np.add(t1, t2, dtype=f4)
        else:
            v = (1 - r) * lo.astype(np.float64) + r * hi.astype(np.float64)
    if dt.kind in "iu":
        v = np.asarray(dtypes.cast(v, dt.name))
    v = np.where((r == 0) | (lo == hi), lo, v)
    return v.astype(dt)


def prctile(x, p, axis=None, *, matlab_layout: bool = False):
    """MATLAB's ``prctile(x, p, dim)`` (module docstring): for ``axis=None`` the values for each
    ``p`` (1-D); along an axis, that axis holds the values for each ``p`` in order. NaN is left
    out; with no values the result is NaN. ``p`` in [0, 100]."""
    a = _r._array(x)
    if a.dtype.kind == "b":
        a = a.astype(np.float64)
    dt = a.dtype
    pv = np.atleast_1d(np.asarray(p, dtype=np.float64)).ravel()
    if np.any((pv < 0) | (pv > 100) | np.isnan(pv)):
        raise ValueError("p must be between 0 and 100")
    ax = _axis(a, axis)
    if ax is None:
        v, _ = _values(a, True)
        n = v.size
        if n == 0:
            return np.full(pv.size, _nan_of(dt))
        k, r = _prc_plan(n, pv)
        klo = np.clip(k, 1, n) - 1
        khi = np.clip(k + 1, 1, n) - 1
        need = sorted({int(i) for i in np.concatenate([klo, khi])})
        kv = dict(zip(need, _kth_whole(a, v, need, matlab_layout), strict=True))
        lo = np.array([kv[int(i)] for i in klo], dtype=dt)
        hi = np.array([kv[int(i)] for i in khi], dtype=dt)
        r = np.where((k < 1) | (k >= n), 0.0, r)
        return _interp(lo, hi, r, dt)
    if a.ndim == 0:
        a = a.reshape(1)
    s, nv = _sorted_along(a, ax)
    k, r = _prc_plan(nv[..., None], pv)  # (..., np)
    nvb = np.broadcast_to(nv[..., None], k.shape)
    klo = np.clip(k, 1, np.maximum(nvb, 1)) - 1
    khi = np.clip(k + 1, 1, np.maximum(nvb, 1)) - 1
    lo = np.take_along_axis(s, klo.astype(np.intp), axis=-1)
    hi = np.take_along_axis(s, khi.astype(np.intp), axis=-1)
    r = np.where((k < 1) | (k >= nvb), 0.0, r)
    out = _interp(lo, hi, r, dt)
    if (nvb == 0).any():
        out = np.where(nvb == 0, _nan_of(dt), out)
    return np.ascontiguousarray(np.moveaxis(out, -1, ax))


# -- mode --------------------------------------------------------------------------------------


def _mode_1d(v: np.ndarray, first_of_run: bool):
    n = v.size
    if n == 0:
        return _nan_of(v.dtype), 0.0
    s = np.sort(v, kind="stable")
    new = np.ones(n, bool)
    new[1:] = s[1:] != s[:-1]
    if s.dtype.kind == "f":
        new[1:] |= np.isnan(s[1:])
    starts = np.flatnonzero(new)
    counts = np.diff(np.append(starts, n))
    i = int(np.argmax(counts))
    j = starts[i] if first_of_run else starts[i] + counts[i] - 1
    return s[j], float(counts[i])


def mode(
    x,
    axis=None,
    *,
    with_count: bool = False,
    first_of_run: bool = False,
    keepdims: bool = False,
    matlab_layout: bool = False,
):
    """MATLAB's ``mode(x, dim)`` (module docstring): ``m``, or ``(m, f)`` with ``with_count``
    (``f`` double). ``first_of_run=True`` gives the value of MATLAB's three-output call."""
    a = _r._array(x)
    ax = _axis(a, axis)
    if ax is None:
        v = _r._fortran_flat(_r._mview(a, matlab_layout), a.dtype) if a.ndim >= 2 else a.ravel()
        val, cnt = _mode_1d(v, first_of_run)
        if keepdims:
            val, cnt = np.full((1,) * a.ndim, val), np.full((1,) * a.ndim, cnt)
        return (val, np.float64(cnt)) if with_count else val
    if a.ndim == 0:
        a = a.reshape(1)
    w = np.moveaxis(a, ax, -1)
    flat = w.reshape(-1, w.shape[-1])
    res = [_mode_1d(np.ascontiguousarray(flat[i]), first_of_run) for i in range(flat.shape[0])]
    out_dt = np.float64 if (w.shape[-1] == 0 and a.dtype.kind in "biu") else a.dtype
    vals = np.array([v for v, _ in res], dtype=out_dt).reshape(w.shape[:-1])
    cnt = np.array([c for _, c in res], dtype=np.float64).reshape(w.shape[:-1])
    if keepdims:
        vals, cnt = np.expand_dims(vals, ax), np.expand_dims(cnt, ax)
    vals = vals[()] if vals.ndim == 0 else vals
    cnt = cnt[()] if cnt.ndim == 0 else cnt
    return (vals, cnt) if with_count else vals
