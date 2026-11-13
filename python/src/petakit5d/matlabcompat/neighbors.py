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
``hypot``, ``atan2``, ``sin`` and ``cos`` are ``matlabcompat.elementary``'s (fdlibm's and msun's,
as MATLAB computes them, written out in NumPy), so that the results are MATLAB's on every C
library and CPU: ``hypot`` with MATLAB's NaN for a NaN
argument, also next to Inf, and ``atan2`` of two zeros +0 as in MATLAB (C gives -0 or +-pi for
signed zeros). The C libraries' ``hypot`` (glibc 2.35 and later), ``atan2``, ``sin`` and ``cos``
differ from MATLAB's in the last bit for a part of the inputs (MATLAB's ``atan2`` is also 1 ulp
above +-pi/2 for a negative x tiny next to y); the coordinates of ``sph2cart`` would carry such
differences of the angles, amplified where a cosine or sine is near zero (up to about 5e-16 of the
radius). ``sph2cart`` takes fdlibm's sine and cosine as they are (``elementary.fdlibm_sin``,
``fdlibm_cos``).
"""

from __future__ import annotations

import numpy as np

from petakit5d.core import dtypes
from petakit5d.core.threads import resolve_threads
from petakit5d.matlabcompat.elementary import atan2, fdlibm_cos, fdlibm_sin, hypot

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


def knnsearch(X, Y, k: int = 1, *, num_threads: int | None = None):
    """``[idx, d] = knnsearch(X, Y, 'K', k)``: (idx (ny, min(k, nx)) int64, 0-based; d float64).
    Rows are points; a 1-D array is one point (a row). Rows of ``X`` with a NaN are never
    neighbors; a query with a NaN gets the first rows of ``X`` and NaN distances (MATLAB).
    ``num_threads`` (Python only; None: the CPUs available) are the k-d tree's query threads;
    the result doesn't depend on them."""
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
    nt = resolve_threads(num_threads)
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
            i, d = _tree(Xv, Yq, kk, nt)
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


def _tree(X, Y, k, workers=1):
    """candidates from a k-d tree (queries on ``workers`` threads): every point within the k-th
    distance (and a margin), ranked exactly"""
    from scipy.spatial import cKDTree

    tree = cKDTree(X)
    extra = 8
    m = min(X.shape[0], k + extra)
    dist, cand = tree.query(Y, k=m, workers=workers)
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


def cart2sph(x, y, z):
    """``[az, elev, r] = cart2sph(x, y, z)`` (float64)"""
    x, y, z = (np.asarray(v, dtype=np.float64) for v in (x, y, z))
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    elev = atan2(z, hxy)
    az = atan2(y, x)
    return az, elev, r


def sph2cart(az, elev, r):
    """``[x, y, z] = sph2cart(az, elev, r)`` (float64)"""
    az, elev, r = (np.asarray(v, dtype=np.float64) for v in (az, elev, r))
    z = _mul(r, fdlibm_sin(elev))
    rcoselev = _mul(r, fdlibm_cos(elev))
    x = _mul(rcoselev, fdlibm_cos(az))
    y = _mul(rcoselev, fdlibm_sin(az))
    return x, y, z


def _mul(a, b):
    """a * b with a NaN operand's payload kept as x86 does (the first NaN operand); NumPy may
    swap the operands of a product"""
    with np.errstate(invalid="ignore"):
        p = dtypes.x86_nan(a * b, a, b)
    return np.where(np.isnan(a), a, np.where(np.isnan(b), b, p))
