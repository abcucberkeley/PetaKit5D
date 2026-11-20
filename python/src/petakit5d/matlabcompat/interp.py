"""MATLAB's ``interp1`` with the ``'linear'`` method, as track processing and the lifetime analysis
call it (R2024b; measured by black-box probes, task C04a).

``interp1(x, v, xq)``, ``interp1(x, v, xq, 'linear')``, ``interp1(x, v, xq, 'linear', 'extrap')``
and ``interp1(x, v, xq, 'linear', extrapval)`` for vectors ``x`` and ``v`` of the same length:

- a query in ``[x(k), x(k+1)]`` gives ``(1 - s)*v(k) + s*v(k+1)`` with
  ``s = (xq - x(k)) / (x(k+1) - x(k))``, each operation one IEEE rounding (bit for bit MATLAB's on
  724 probed calls: two-point grids, uneven grids of 3 to 41 points, 999 queries on 910 points);
  a query on a knot uses the interval that starts there (the last one at the end); a term of
  weight 0 whose value is Inf or NaN is left out (a query on a knot next to one gives the
  knot's value, -0 included), a finite one stays (``-0 + 0*v``); equal neighbours give the
  first as it is (MATLAB's -0 between -0 and +0);
- queries outside ``[x(1), x(end)]`` give NaN, the extrapolation value, or with ``'extrap'`` the
  same formula on the first or last interval; a NaN query gives NaN in every form;
- a decreasing or unsorted grid is sorted first (``v`` with it); repeated points, NaN or Inf
  points and a one-point grid raise MATLAB's errors (``identifier``); integer ``x`` or ``v`` raise
  ``MATLAB:interp1:Xnumeric`` and ``MATLAB:interp1:NonFloatValues``;
- the result has the shape of ``xq``; every NaN is MATLAB's (sign bit set). All single: single
  arithmetic; a single grid or single queries with double values: double, then single (both
  measured); single values with a double grid or queries are not ported (the probes did not
  resolve MATLAB's mixed arithmetic there).

Only the linear method is ported (the traced code calls no other). ``x``, ``v`` and ``xq`` are
vectors and matrices of values, not images: their shapes are kept as they are.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from petakit5d.core import dtypes

__all__ = ["InterpError", "interp1"]


class InterpError(ValueError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _float_vector(a: Any, ident: str, what: str) -> np.ndarray:
    v = np.asarray(a)
    if v.dtype.kind not in "f" or v.dtype.itemsize not in (4, 8):
        raise InterpError(f"The {what} must be of type double or single.", ident)
    if v.ndim > 2 or (v.ndim == 2 and min(v.shape) > 1):
        raise NotImplementedError("interp1: only vectors x and v are ported")
    return v.reshape(-1)


def interp1(x: Any, v: Any, xq: Any, method: str = "linear", extrap: Any = None) -> np.ndarray:
    """MATLAB's ``interp1(x, v, xq, method, extrapolation)`` for the linear method (module
    docstring). ``extrap`` is None (NaN outside the grid), ``'extrap'`` or a value."""
    if not isinstance(method, str) or method.lower() != "linear":
        raise NotImplementedError(f"interp1: only the 'linear' method is ported, not {method!r}")
    xa = _float_vector(x, "MATLAB:interp1:Xnumeric", "grid vector x")
    va = _float_vector(v, "MATLAB:interp1:NonFloatValues", "values v")
    q = np.asarray(xq)
    if q.dtype.kind not in "f":
        q = q.astype(np.float64)
    single = any(a.dtype == np.float32 for a in (xa, va, q))
    # all single: single arithmetic; single grid or queries with double values: double, then
    # single (measured); single values with a double grid or queries: not resolved by the probes
    all_single = all(a.dtype == np.float32 for a in (xa, va, q))
    if single and not all_single and va.dtype == np.float32:
        raise NotImplementedError("interp1: single values with a double grid or queries")
    dt = np.float32 if all_single else np.float64
    if xa.size != va.size:
        raise InterpError("The grid vectors and the values do not have the same length.")
    if xa.size < 2:
        raise InterpError(
            "Interpolation requires at least two sample points for each grid dimension.",
            "MATLAB:griddedInterpolant:DegenerateGridErrId",
        )
    if not np.all(np.isfinite(xa)):
        raise InterpError(
            "The grid coordinates must be finite.",
            "MATLAB:griddedInterpolant:NonFiniteInputPtsErrId",
        )
    xs = xa.astype(dt)
    vs = va.astype(dt)
    if not np.all(np.diff(xs) > 0):
        if np.all(np.diff(xs) < 0):
            xs, vs = xs[::-1], vs[::-1]
        else:
            order = np.argsort(xs, kind="stable")
            xs, vs = xs[order], vs[order]
        if np.any(np.diff(xs) == 0):
            raise InterpError(
                "The grid vectors must contain unique points.",
                "MATLAB:griddedInterpolant:NonUniqueCompVecsPtsErrId",
            )
    qq = q.astype(dt)
    flat = qq.reshape(-1)
    n = xs.size
    k = np.clip(np.searchsorted(xs, flat, side="right") - 1, 0, n - 2)
    x0, x1 = xs[k], xs[k + 1]
    v0, v1 = vs[k], vs[k + 1]
    with np.errstate(invalid="ignore", over="ignore"):
        s = (flat - x0) / (x1 - x0)
        t0, t1 = (dt(1) - s) * v0, s * v1
        out = t0 + t1
    # a term of weight 0 whose value is Inf or NaN is left out: a query on a knot next to one
    # gives the knot's value as it is (-0 included); finite terms of weight 0 stay (measured)
    left_out_1 = (s == 0) & ~np.isfinite(v1)
    out[left_out_1] = t0[left_out_1]
    left_out_0 = (s == 1) & ~np.isfinite(v0)
    out[left_out_0] = t1[left_out_0]
    # equal neighbours give the first one as it is (MATLAB's -0 between -0 and +0)
    same = v0 == v1
    out[same] = v0[same]
    nan = np.isnan(flat)
    outside = ~nan & ((flat < xs[0]) | (flat > xs[-1]))
    if isinstance(extrap, str):
        if extrap.lower() != "extrap":
            raise InterpError(f"Invalid extrapolation method {extrap!r}.")
    elif extrap is None:
        out[outside] = dt(dtypes.NAN)
    else:
        out[outside] = dt(np.asarray(extrap).reshape(-1)[0])
    # every NaN result is MATLAB's (sign bit set), as measured
    out[nan | np.isnan(out)] = dt(dtypes.NAN)
    if single and not all_single:
        out = out.astype(np.float32)
    return out.reshape(qq.shape)
