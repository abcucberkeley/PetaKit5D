"""MATLAB's ``corr`` (Statistics and Machine Learning Toolbox) as ``XZoffsetEstimation3D`` calls
it: the Pearson correlation of two columns, ``corr(x, y)`` with the defaults (``'type'``
``'Pearson'``, ``'rows'`` ``'all'``). From the documentation and black-box probes (task C04a; no
MathWorks source read).

MATLAB computes ``(xc' * yc) / norm(xc) / norm(yc)`` on the centered columns ``xc = x -
mean(x)`` (measured: 200 of 200 probes from MATLAB's own pieces), through its BLAS, whose order
of additions depends on the CPU and, from 4,096 values, on MATLAB's thread count (its results
differ there). The port takes ``sum(xc .* yc) / nx / ny`` with ``nx = sqrt(sum(xc .^ 2))``
(scaled by the largest magnitude where the plain sum of squares would overflow or underflow),
the means and sums in MATLAB's order (``matlabcompat.reduce``), so that its result depends on
neither the platform nor the thread count: bit for bit where the sums are exact (the probed edge
values: perfect correlations that MATLAB gives as 1 - eps/2, the overflow-free norm, 0 and NaN),
otherwise within a few ulp (up to 17 where the sum cancels in small integer probes, 16 at 65,536
values: parity row W1 of ``docs/parity/detection.md``). 0/0 gives MATLAB's NaN; a NaN of the
data propagates. Empty columns raise ``stats:corr:TooFewInputs``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from petakit5d.core import dtypes

__all__ = ["CorrError", "corr"]


class CorrError(ValueError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _column(a: Any, what: str) -> np.ndarray:
    v = np.asarray(a)
    if v.ndim > 2 or (v.ndim == 2 and v.shape[1] != 1):
        raise NotImplementedError(f"corr: only column vectors are ported ({what} is {v.shape})")
    if v.dtype.kind != "f" or v.dtype.itemsize not in (4, 8):
        raise NotImplementedError(f"corr: only double and single data are ported ({what})")
    return v.reshape(-1)


def _norm(c: np.ndarray, dt) -> Any:
    from petakit5d.matlabcompat import reduce as mr

    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        ss = mr.sum(c * c).reshape(-1)[0]
        if np.isfinite(ss) and ss >= np.finfo(dt).tiny * 2**52:
            return dt(np.sqrt(ss))
        m = dt(np.max(np.abs(c))) if c.size else dt(0)
        if not np.isfinite(m) or m == 0:
            return dt(np.sqrt(ss))
        t = c / m
        return dt(m * np.sqrt(mr.sum(t * t).reshape(-1)[0]))


def corr(x: Any, y: Any) -> Any:
    """The Pearson correlation of the column vectors ``x`` and ``y`` (module docstring): a
    scalar of double, or single when either is single."""
    from petakit5d.matlabcompat import reduce as mr

    a = _column(x, "x")
    b = _column(y, "y")
    if a.size != b.size:
        raise CorrError(
            "X and Y must have the same number of rows.", "stats:corr:InputSizeMismatch"
        )
    if a.size == 0:
        raise CorrError("Too few inputs.", "stats:corr:TooFewInputs")
    single = a.dtype == np.float32 or b.dtype == np.float32
    dt = np.float32 if single else np.float64
    a = a.astype(dt, copy=False)
    b = b.astype(dt, copy=False)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        ac = a - mr.mean(a).reshape(-1)[0]
        bc = b - mr.mean(b).reshape(-1)[0]
        sab = mr.sum(ac * bc).reshape(-1)[0]
        r = dt(dt(sab / _norm(ac, dt)) / _norm(bc, dt))
    if np.isnan(r) and not (np.isnan(a).any() or np.isnan(b).any()):
        return dt(dtypes.NAN)  # 0/0: MATLAB's NaN; a NaN of the data propagates (measured)
    return r
