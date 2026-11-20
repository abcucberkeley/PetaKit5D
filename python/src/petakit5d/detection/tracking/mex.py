"""The point-set mex functions of the tracker and the detection code, written from their method and
their measured behaviour (no source was read):

- ``create_distance_matrix(A, B)`` (``createDistanceMatrix``): the Euclidean distances between
  the rows of ``A`` (m, d) and ``B`` (n, d), ``sqrt(((a1 - b1)^2 + (a2 - b2)^2) + ...)`` in that
  order (bit for bit on 8.75 million recorded distances); empty inputs give (m, n) empty.
- ``kd_tree_ball_query(X, Y, r)`` (``KDTreeBallQuery``): for each row of ``Y`` the rows of ``X``
  within its radius (inclusive), nearest first, with their distances (computed as above).
  MATLAB returns 1-based indices as double column vectors and ``[]`` (0 x 0) where none is in
  range; here 0-based ``int64`` arrays. Equal distances: MATLAB's order comes from its k-d tree
  and is not reproduced; here they are in row order (a documented difference that needs two
  neighbours at exactly the same distance).
"""

from __future__ import annotations

import numpy as np

__all__ = ["create_distance_matrix", "kd_tree_ball_query"]


def _points(a, name: str) -> np.ndarray:
    x = np.asarray(a, dtype=np.float64)
    if x.ndim == 1:
        x = x.reshape(1, -1) if x.size else x.reshape(0, 0)
    if x.ndim != 2:
        raise ValueError(f"{name} must be a matrix of points (rows)")
    return x


def _dist(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    d = A.shape[1]
    s = (A[:, None, 0] - B[None, :, 0]) ** 2
    for k in range(1, d):
        s = s + (A[:, None, k] - B[None, :, k]) ** 2
    return np.sqrt(s)


def create_distance_matrix(A, B=None) -> np.ndarray:
    """``createDistanceMatrix(A, B)``: (m, n) Euclidean distances (``B`` defaults to ``A``)"""
    A = _points(A, "A")
    B = A if B is None else _points(B, "B")
    m = A.shape[0]
    n = B.shape[0]
    if m == 0 or n == 0:
        return np.zeros((m, n))
    if A.shape[1] != B.shape[1]:
        raise ValueError("createDistanceMatrix: A and B must have the same number of columns")
    return _dist(A, B)


def kd_tree_ball_query(X, Y, radius) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """``[idx, dist] = KDTreeBallQuery(X, Y, radius)``: per row of ``Y``, the 0-based rows of
    ``X`` within ``radius`` (a scalar or one per row of ``Y``), nearest first, and their
    distances"""
    from scipy.spatial import cKDTree

    X = _points(X, "X")
    Y = _points(Y, "Y")
    r = (
        np.broadcast_to(np.asarray(radius, dtype=np.float64).ravel(), (Y.shape[0],))
        if Y.shape[0]
        else np.zeros(0)
    )
    idx: list[np.ndarray] = []
    dist: list[np.ndarray] = []
    if X.shape[0] == 0 or Y.shape[0] == 0:
        return [np.zeros(0, np.int64) for _ in range(Y.shape[0])], [
            np.zeros(0) for _ in range(Y.shape[0])
        ]
    tree = cKDTree(X)
    # candidates with a margin, then the exact distances decide
    cand = tree.query_ball_point(Y, r * (1 + 1e-12) + 1e-300)
    for q, c in enumerate(cand):
        c = np.asarray(sorted(c), dtype=np.int64)
        if c.size == 0:
            idx.append(c)
            dist.append(np.zeros(0))
            continue
        d = _dist(Y[q : q + 1], X[c])[0]
        keep = d <= r[q]
        c, d = c[keep], d[keep]
        o = np.argsort(d, kind="stable")
        idx.append(c[o])
        dist.append(d[o])
    return idx, dist
