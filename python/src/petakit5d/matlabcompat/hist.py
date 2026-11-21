"""MATLAB's ``hist(y, x)`` with bin centers, as measured from MATLAB R2024b (task C05; black-box
probes, no MathWorks source read).

The bins' edges are the midpoints of adjacent centers, ``(x(k) + x(k+1))/2``; a value on an edge
counts in the lower bin (``x(k) < y <= x(k+1)`` between the edges, so bin k holds the values in
``(e(k-1), e(k)]``); the first and last bins are open (``-Inf`` and ``Inf`` count there); NaN is
left out. A vector gives a row of counts (double), a matrix one column of counts per column (an
``n x numel(x)`` result in MATLAB's layout: ``P[i - 1, j - 1] = N(i, j)``). Integer and logical
data count as their values.
"""

from __future__ import annotations

import numpy as np

__all__ = ["hist"]


def _edges(x: np.ndarray) -> np.ndarray:
    return (x[:-1] + x[1:]) / 2


def _counts(y: np.ndarray, edges: np.ndarray, n: int) -> np.ndarray:
    v = y[~np.isnan(y)] if y.dtype.kind == "f" else y
    # bin of each value: the number of edges below it (a value on an edge stays in the lower bin)
    b = np.searchsorted(edges, v, side="left")
    return np.bincount(b, minlength=n).astype(np.float64)[:n]


def hist(y, x, *, matlab_layout: bool = False) -> np.ndarray:
    """MATLAB's ``N = hist(y, x)`` for bin centers ``x`` (a vector, sorted ascending; module
    docstring). ``y``: a vector (1-D, or a MATLAB row or column), giving a 1-D array of counts, or
    with ``matlab_layout=True`` a matrix whose columns are counted separately, giving the
    ``numel(x)``-row matrix of counts."""
    c = np.asarray(x, dtype=np.float64).ravel()
    if c.size == 0:
        raise ValueError("hist needs at least one bin center")
    if np.any(np.diff(c) < 0) or np.isnan(c).any():
        raise ValueError("the bin centers must be sorted ascending")
    e = _edges(c)
    a = np.asarray(y)
    if a.dtype.kind not in "biuf":
        raise TypeError("hist counts numeric or logical data")
    a = a.astype(np.float64)
    if matlab_layout and a.ndim == 2 and min(a.shape) > 1:
        return np.column_stack([_counts(a[:, j], e, c.size) for j in range(a.shape[1])])
    return _counts(a.ravel(), e, c.size)
