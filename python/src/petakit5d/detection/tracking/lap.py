"""The linear assignment problems of the tracker (MATLAB's ``lap`` with the ``mexLap`` solver).

``mexLap`` (no source) solves the sparse linear assignment problem with the shortest augmenting
path method of Jonker and Volgenant (LAPJVsp). Measured as a black box: ``lap(cc, -5, 0)`` takes
the square cost matrix with ``-5`` where a link is not allowed, keeps the other entries (in
compressed rows, 1-based, a dummy element 0) and returns, for every row, its column. Here the
same problem is solved by SciPy's ``min_weight_full_bipartite_matching`` (BSD, LAPJVsp too).

Tie rule: on equal totals the assignment is the one SciPy's LAPJVsp returns for the rows in
order and the columns of each row in increasing order (the input MATLAB gives ``mexLap``). On
the demo data's 178 problems (tracking, 2,000 to 47,000 rows) SciPy's assignment equals
MATLAB's in 169; in the other 9 the two differ only between equal-cost entries of the
birth/death blocks (no link differs); on 60 small problems with integer costs and many ties
the two are the same. A cost of 0 is a link for MATLAB (a probe); SciPy needs nonzero weights,
so a 0 is given the smallest positive normal number (totals change by less than any cost
difference, so only ties between assignments with different numbers of 0 links could differ).
"""

from __future__ import annotations

import numpy as np

__all__ = ["solve_assignment"]

_TINY = float(np.finfo(np.float64).tiny)


def solve_assignment(rows: np.ndarray, cols: np.ndarray, vals: np.ndarray, n: int) -> np.ndarray:
    """the column assigned to each row of the ``n x n`` problem given by its allowed entries
    (0-based ``rows``, ``cols``, costs ``vals``, no duplicates), as an ``int64`` array"""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import min_weight_full_bipartite_matching

    vals = np.asarray(vals, dtype=np.float64)
    if (vals == 0).any():
        vals = np.where(vals == 0, _TINY, vals)
    m = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    m.sort_indices()
    r, c = min_weight_full_bipartite_matching(m)
    x = np.empty(n, np.int64)
    x[r] = c
    return x
