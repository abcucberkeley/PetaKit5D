"""MATLAB's ``lsqlin`` with its default (interior-point) algorithm, as the stitching shift
assignment calls it: ``x = lsqlin(C, d, A, b)`` or ``lsqlin(C, d, A, b, Aeq, beq)``, the minimum of
``||C x - d||^2`` subject to ``A x <= b`` and ``Aeq x = beq``.

MATLAB's solver is iterative and closed (``ipqpdense``, ``ipqpsparse``), so it is not copied: the
optimum here is exact (an interior-point method, then the equality-constrained problem of its
active set solved again by regularized iterative refinement, so it is reached to rounding). Where
the optimum is not unique, the member returned is the one MATLAB's solvers come closest to:

- **Lineality** (adding a constant to every variable of a group that the rows of ``C``, ``A`` and
  ``Aeq`` connect changes nothing; for shifts: all tiles of a group of overlapping tiles): the mean
  of each group is zero. MATLAB's dense solver starts at mean zero and drifts from it by rounding
  noise that changes with MATLAB's thread count (up to 0.02 on generated grids, 0.065 on the demo
  datasets, more with small weights); its sparse solver ends anywhere, which is why the sparse
  caller subtracts the mean (``x - mean(x)``).
- **Free offsets** (groups of variables that ``C`` and ``Aeq`` don't connect to the rest, held
  only by inequalities that are not active: tiles whose pairs all fall below the correlation
  threshold): the analytic center of the optimal set, the limit of the interior-point path.
  MATLAB stops on its path after a few iterations, a voxel or more from it (up to 8 on the demo
  datasets); these values are flagged (``free``).

``solve`` also flags values near a half integer, where rounding the solution can differ from
MATLAB's (``near_half``; the window is ``LsqlinResult.tol``).

The rows are expected as the stitching code makes them: each row of ``C``, ``A`` and ``Aeq`` is
zero or ``v * (e_j - e_i)`` (a weighted difference of two variables). Other rows are solved the
same way with dense null spaces (for small problems). Results never depend on the thread count
(sparse LU and sparse products only).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "lsqlin",
    "solve",
    "LsqlinResult",
    "LsqlinError",
    "near_half_integer",
    "HALF_TOL",
    "HALF_REL",
    "FREE_TOL",
]

# values within this of a half integer are flagged, by form of MATLAB's solver: the sparse
# solution is moved to mean zero by its caller and agrees to about 3e-5; the dense one drifts from
# mean zero by rounding noise that changes with MATLAB's thread count: up to 0.017 on generated
# grids (shifts up to about 30), 0.065 on the demo datasets (shifts up to about 110), so the window
# is the larger of HALF_TOL and HALF_REL times the largest |x| of the lineality group (callers
# with other scales pass their own: the group method's step, weights 0.01, drifts up to 0.19)
HALF_TOL = {"dense": 0.02, "sparse": 0.005}
HALF_REL = {"dense": 1e-3, "sparse": 0.0}
# how far MATLAB's free values land from the analytic center (measured: up to 1.65 voxels dense,
# 1.51 sparse); the mean-zero rule passes their mean on to the other values of their group
FREE_TOL = 2.0

_TOL = 1e-10  # interior-point stopping tolerance (relative); the active set's solve does the rest
_MAXIT = 200


class LsqlinError(ValueError):
    """An error of the inputs, or a problem without a feasible point; ``identifier`` is MATLAB's
    error id where MATLAB raises one."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


@dataclass
class LsqlinResult:
    """The solution and what is known about it.

    ``x``: the solution. ``free``: values in a group of variables whose offset only inactive
    inequalities hold, other than the largest such group of their lineality group (MATLAB's value
    is where its iterations stopped, ours the analytic center), and, for the sparse form, every
    value when there is more than one lineality group (MATLAB's offsets of the groups are
    arbitrary and its caller subtracts one mean). ``near_half``: values within ``tol`` of a half
    integer, where rounding can differ from MATLAB's; ``tol`` is the window of the value's
    lineality group (the larger of ``half_tol`` and ``HALF_REL[form]`` times the group's largest
    |x|), widened where the group has free values (the mean-zero rule passes their deviation on:
    the window times the group's size over its main part's, plus ``FREE_TOL`` times its free
    values over its main part's). ``active``: rows of ``A`` that hold with equality. ``groups``: the
    lineality group of each variable (0-based, numbered by first variable). ``iterations``:
    interior-point iterations (0 when the unconstrained optimum is feasible)."""

    x: np.ndarray
    near_half: np.ndarray
    free: np.ndarray
    active: np.ndarray
    groups: np.ndarray
    iterations: int = 0
    form: str = "dense"
    tol: np.ndarray = field(default_factory=lambda: np.zeros(0))
    info: dict = field(default_factory=dict)


def lsqlin(C, d, A=None, b=None, Aeq=None, beq=None) -> np.ndarray:
    """``x = lsqlin(C, d, A, b, Aeq, beq)``: see the module. A SciPy sparse ``C`` is MATLAB's
    sparse ``C`` (its sparse solver); the solution is the same member either way."""
    return solve(C, d, A, b, Aeq, beq).x


def solve(C, d, A=None, b=None, Aeq=None, beq=None, *, form=None, half_tol=None) -> LsqlinResult:
    """The solution with its flags. ``form`` is ``'dense'`` or ``'sparse'`` (default: sparse if
    ``C`` is a SciPy sparse matrix, as MATLAB chooses); ``half_tol`` defaults to
    ``HALF_TOL[form]``."""
    import scipy.sparse as sp

    if form is None:
        form = "sparse" if sp.issparse(C) else "dense"
    if form not in HALF_TOL:
        raise ValueError(f"form must be 'dense' or 'sparse', not {form!r}")
    if half_tol is None:
        half_tol = HALF_TOL[form]
    P = _Problem(C, d, A, b, Aeq, beq)
    if P.p == 0 and P.q == 0:
        x, active, it, sup = _unconstrained(P), np.zeros(0, bool), 0, None
    else:
        x, active, it, sup = _solve_constrained(P)
    free, tol = _flags(P, x, sup, form, half_tol)
    return LsqlinResult(
        x,
        near_half_integer(x, tol),
        free,
        active,
        P.lineality,
        it,
        form,
        tol,
        {"structured": P.structured},
    )


def _flags(P, x, sup, form, half_tol):
    """free values (outside the main group of their lineality group, or all for the sparse form
    with several lineality groups) and the rounding tolerance of each value: the window of the
    group (the larger of ``half_tol`` and HALF_REL times its largest |x|) scaled to its main part,
    plus FREE_TOL times its share of free values"""
    n = P.n
    lab = P.lineality
    free = np.zeros(n, bool)
    tol = np.full(n, float(half_tol))
    for g in range(lab.max() + 1):
        m = lab == g
        base = max(float(half_tol), HALF_REL[form] * float(np.abs(x[m]).max()))
        tol[m] = base
        if sup is None:
            continue
        vals, cnt = np.unique(sup[m], return_counts=True)
        if vals.size < 2:
            continue
        main = vals[np.argmax(cnt)]  # the largest; ties to the one of the lowest variable
        nm = int(cnt.max())
        nf = int(m.sum()) - nm
        free[m & (sup != main)] = True
        tol[m] = base * m.sum() / nm + FREE_TOL * nf / nm
    if form == "sparse" and P.ndirections > 1:
        free[:] = True
    return free, tol


def near_half_integer(x, tol: float) -> np.ndarray:
    """True where ``x`` lies within ``tol`` of a half integer (``k + 0.5``)."""
    x = np.asarray(x, dtype=np.float64)
    return np.abs(x - np.floor(x) - 0.5) < tol


# ------------------------------------------------------------------------------------------------
# the problem and its structure
# ------------------------------------------------------------------------------------------------


class _Problem:
    def __init__(self, C, d, A, b, Aeq, beq):
        import scipy.sparse as sp

        C = _csr(C, None, "C")
        if C.shape[0] == 0 or C.shape[1] == 0 or np.size(d) == 0:
            raise LsqlinError(
                "The first two arguments to LSQLIN cannot be empty matrices.",
                "optimlib:lsqlin:FirstTwoArgsEmpty",
            )
        n = C.shape[1]
        d = _vec(d, C.shape[0], "d", "optimlib:lsqlin:InvalidCAndD")
        A = _csr(A, n, "A")
        b = _vec(b, A.shape[0], "b", "optimlib:lsqlin:InvalidAAndB")
        Aeq = _csr(Aeq, n, "Aeq")
        beq = _vec(beq, Aeq.shape[0], "beq", "optimlib:lsqlin:InvalidAeqAndBeq")
        for name, v in (
            ("C", C.data),
            ("d", d),
            ("A", A.data),
            ("b", b),
            ("Aeq", Aeq.data),
            ("beq", beq),
        ):
            if not np.all(np.isfinite(v)):
                raise LsqlinError(f"{name} must be finite")
        self.C, self.d, self.A, self.b, self.Aeq, self.beq = C, d, A, b, Aeq, beq
        self.n, self.p, self.q = n, A.shape[0], Aeq.shape[0]
        cs = _difference_rows(C)
        as_ = _difference_rows(A)
        es = _difference_rows(Aeq)
        self.structured = cs is not None and as_ is not None and es is not None
        if self.structured:
            from scipy.sparse.csgraph import connected_components

            ce, ae, ee = cs[0], as_[0], es[0]
            obj = np.concatenate([ce, ee], axis=0)
            allr = np.concatenate([ce, ae, ee], axis=0)
            self.lineality = _labels(n, allr, connected_components, sp)
            self.ndirections = int(self.lineality.max()) + 1
            self.ccomp = _labels(n, obj, connected_components, sp)
            self.a_ends = ae  # (p, 2) ends of each A row (i, j), -1 for zero rows
            self.a_val = as_[1]  # A row = v * (e_j - e_i)
        else:
            self.lineality = None
            self.ccomp = None
            self._dense_structure()

    def _dense_structure(self):
        """lineality groups and bases for rows that are not differences (small problems)"""
        if self.n > 5000:
            raise LsqlinError(
                "rows that are not weighted differences are supported up to 5000 variables"
            )
        M = np.vstack([self.C.toarray(), self.A.toarray(), self.Aeq.toarray()])
        self.L = _null(M)  # lineality basis, orthonormal
        self.N = _null(np.vstack([self.C.toarray(), self.Aeq.toarray()]))
        # groups: variables that some lineality vector moves (connected through L's support)
        import scipy.sparse as sp
        from scipy.sparse.csgraph import connected_components

        vi, ci = np.nonzero(np.abs(self.L) > 1e-12)
        nn = self.n + self.L.shape[1]
        G = sp.csr_matrix((np.ones(vi.size), (vi, self.n + ci)), shape=(nn, nn))
        _, lab = connected_components(G, directed=False)
        lab = lab[: self.n]
        # the variables no lineality direction moves form one group
        moved = np.zeros(self.n, bool)
        moved[vi] = True
        if not moved.all():
            lab[~moved] = lab[~moved][0]
        _, first = np.unique(lab, return_index=True)
        rank = np.empty(lab.max() + 1, np.int64)
        rank[lab[np.sort(first)]] = np.arange(first.size)
        self.lineality = rank[lab]
        self.ndirections = self.L.shape[1]


def _csr(M, n, name):
    import scipy.sparse as sp

    if M is None:
        return sp.csr_matrix((0, 0 if n is None else n))
    if not sp.issparse(M) and np.size(M) == 0:
        # an empty matrix keeps its columns; MATLAB's [] (0 x 0) takes n
        shape = np.shape(M)
        cols = shape[1] if len(shape) == 2 and shape[1] else (0 if n is None else n)
        rows = shape[0] if len(shape) == 2 else 0
        return sp.csr_matrix((rows, cols))
    if sp.issparse(M):
        M = sp.csr_matrix(M, dtype=np.float64)
    else:
        M = np.asarray(M, dtype=np.float64)
        if M.ndim == 1:
            M = M[None, :]
        if M.ndim != 2:
            raise LsqlinError(f"{name} must be a matrix")
        M = sp.csr_matrix(M)
    if n is not None and M.shape[1] != n:
        ident = {
            "A": "optimlib:lsqlin:CAndA",
            "Aeq": "optimlib:lsqlin:CAndAeq",
        }.get(name, "")
        raise LsqlinError(f"{name} must have {n} columns", ident)
    M.sum_duplicates()
    M.eliminate_zeros()
    return M


def _vec(v, m, name, ident):
    if v is None:
        v = np.zeros(0)
    v = np.asarray(v, dtype=np.float64).reshape(-1)
    if v.size != m:
        raise LsqlinError(f"{name} must have {m} elements", ident)
    return v


def _difference_rows(M):
    """(ends (m, 2) with -1 for zero rows, values v) if every row is zero or v (e_j - e_i), else
    None"""
    m = M.shape[0]
    nnz = np.diff(M.indptr)
    if np.any((nnz != 0) & (nnz != 2)):
        return None
    ends = -np.ones((m, 2), np.int64)
    val = np.zeros(m)
    r = np.nonzero(nnz == 2)[0]
    if r.size:
        st = M.indptr[r]
        c0, c1 = M.indices[st], M.indices[st + 1]
        v0, v1 = M.data[st], M.data[st + 1]
        if np.any(v0 != -v1):
            return None
        # v (e_j - e_i): j has +v
        pos = v1 > 0
        i = np.where(pos, c0, c1)
        j = np.where(pos, c1, c0)
        ends[r, 0], ends[r, 1] = i, j
        val[r] = np.abs(v1)
    return ends, val


def _labels(n, ends, connected_components, sp):
    """connected components of the rows' ends, numbered in order of their first variable"""
    e = ends[ends[:, 0] >= 0]
    G = sp.csr_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
    _, lab = connected_components(G, directed=False)
    # renumber by first occurrence
    first = np.full(lab.max() + 1, n)
    np.minimum.at(first, lab, np.arange(n))
    order = np.argsort(first, kind="stable")
    rank = np.empty_like(order)
    rank[order] = np.arange(order.size)
    return rank[lab].astype(np.int64)


def _null(M, rtol=1e-10):
    if M.shape[0] == 0:
        return np.eye(M.shape[1])
    _, s, vt = np.linalg.svd(M, full_matrices=True)
    r = int(np.sum(s > rtol * max(1.0, s[0] if s.size else 0.0)))
    return vt[r:].T.copy()


# ------------------------------------------------------------------------------------------------
# solving
# ------------------------------------------------------------------------------------------------


def _unconstrained(P):
    """no inequalities and no equalities: MATLAB solves C \\ d (unique for full column rank)"""
    import scipy.sparse.linalg as sla

    H = (P.C.T @ P.C).tocsc()
    try:
        lu = sla.splu(H, permc_spec="MMD_AT_PLUS_A")
        x = lu.solve(P.C.T @ P.d)
    except RuntimeError:
        x = None
    if x is None or not np.all(np.isfinite(x)) or H.nnz == 0:
        raise LsqlinError("without constraints only full-rank C is supported (MATLAB: C \\ d)")
    return x


def _solve_constrained(P):
    """(x, active rows, iterations, super groups: variables whose offsets move together)"""
    if P.structured:
        r = _without_inequalities(P)
        if r is not None:
            x, sup = r
            return _center_groups(P, x), np.zeros(P.p, bool), 0, sup
    T, expand = _reduction(P)
    Cu, Au, Eu = P.C @ T, P.A @ T, P.Aeq @ T
    it = 0
    if P.p == 0:
        u = _eq_ls(Cu, P.d, Eu, P.beq)
        if u is None:
            raise LsqlinError("no feasible point: the equalities are inconsistent")
        active = np.zeros(0, bool)
    else:
        u, s, z, y, it = _ipm(Cu, P.d, Au, P.b, Eu, P.beq)
        active = z > np.maximum(s, 1e-9 * (1.0 + np.abs(z).max(initial=0.0)))
        u = _polish(Cu, P.d, Au, P.b, Eu, P.beq, u, active, z, y)
    x = expand(u)
    x, active, sup = _face_center(P, x, active)
    return _center_groups(P, x), active, it, sup


def _reduction(P):
    """x = T u without the lineality directions: one variable per lineality group is fixed at 0
    (rows of differences), or u spans its orthogonal complement (other rows)"""
    import scipy.sparse as sp

    n = P.n
    if P.structured:
        lab = P.lineality
        root = np.zeros(n, bool)
        first = np.full(lab.max() + 1, n)
        np.minimum.at(first, lab, np.arange(n))
        root[first] = True
        keep = np.nonzero(~root)[0]
        T = sp.csr_matrix((np.ones(keep.size), (keep, np.arange(keep.size))), shape=(n, keep.size))

        def expand(u):
            x = np.zeros(n)
            x[keep] = u
            return x

        return T, expand
    L = P.L
    Q = _null(L.T) if L.shape[1] else np.eye(n)
    T = sp.csr_matrix(Q)
    return T, lambda u: Q @ u


def _splu(K, symmetric=True):
    """sparse LU (SuperLU: one thread, so the bits don't depend on the thread count); symmetric
    mode first for symmetric matrices, partial pivoting if that meets a zero pivot"""
    import scipy.sparse.linalg as sla

    K = K.tocsc()
    if symmetric:
        try:
            return sla.splu(
                K,
                permc_spec="MMD_AT_PLUS_A",
                diag_pivot_thresh=0.0,
                options={"SymmetricMode": True},
            )
        except RuntimeError:
            pass
    return sla.splu(K)


def _without_inequalities(P):
    """The optimum of C and Aeq alone, if some point of it holds every inequality: then it is the
    optimum, returned with its free offsets at the analytic center. None otherwise."""
    import scipy.sparse as sp

    # a particular solution: one variable per group of C and Aeq fixed (their offsets are free)
    n = P.n
    lab = P.ccomp
    first = np.full(lab.max() + 1, n)
    np.minimum.at(first, lab, np.arange(n))
    root = np.zeros(n, bool)
    root[first] = True
    keep = np.nonzero(~root)[0]
    S = sp.csr_matrix((np.ones(keep.size), (keep, np.arange(keep.size))), shape=(n, keep.size))
    v = _eq_ls(P.C @ S, P.d, P.Aeq @ S, P.beq)
    if v is None:
        return None
    x = np.zeros(n)
    x[keep] = v
    r = _center_offsets(P, P.b - P.A @ x, np.zeros(P.p, bool), need_feasible=True)
    if r is None:
        return None
    c, sup = r
    return x + c[P.ccomp], sup[P.ccomp]


def _eq_ls(C, d, E, e, refine=4):
    """min ||C v - d||^2 subject to E v = e, by regularized iterative refinement of the KKT
    system (quasi-definite, so dependent rows are allowed). None if it does not converge."""
    import scipy.sparse as sp

    r = C.shape[1]
    H = (C.T @ C).tocsr()
    g = C.T @ d
    q = E.shape[0]
    scale = max(1.0, abs(H).max() if H.nnz else 1.0)
    rho = 1e-10 * scale
    K = sp.bmat([[H + rho * sp.identity(r), E.T], [E, -rho * sp.identity(q)]], format="csc")
    try:
        lu = _splu(K, symmetric=False)
    except RuntimeError:
        return None
    v = np.zeros(r)
    w = np.zeros(q)
    for _ in range(refine):
        r1 = g - H @ v - E.T @ w
        r2 = e - E @ v
        sol = lu.solve(np.concatenate([r1, r2]))
        v = v + sol[:r]
        w = w + sol[r:]
    if not np.all(np.isfinite(v)):
        return None
    # residual check: stationarity and equalities to rounding
    r1 = g - H @ v - E.T @ w
    r2 = e - E @ v
    tol = 1e-9 * (1.0 + np.abs(g).max(initial=0.0) + np.abs(e).max(initial=0.0))
    if np.abs(r1).max(initial=0.0) > tol or np.abs(r2).max(initial=0.0) > tol:
        return None
    return v


def _ipm(C, d, G, h, E, e):
    """Mehrotra's predictor-corrector for min 1/2 ||C u - d||^2, G u <= h, E u = e.
    Returns (u, s, z, y, iterations); raises LsqlinError if no feasible point is found."""
    import scipy.sparse as sp

    r = C.shape[1]
    m = G.shape[0]
    q = E.shape[0]
    H = (C.T @ C).tocsr()
    f = -(C.T @ d)
    Gt = G.T.tocsr()
    Et = E.T.tocsr()
    scale = 1.0 + max(
        np.abs(f).max(initial=0.0), np.abs(h).max(initial=0.0), np.abs(e).max(initial=0.0)
    )
    rho = 1e-12 * max(1.0, abs(H).max() if H.nnz else 1.0)
    u = np.zeros(r)
    s = np.maximum(np.abs(h - G @ u), 1.0)
    z = np.ones(m)
    y = np.zeros(q)
    I_r = sp.identity(r, format="csr")
    I_q = sp.identity(q, format="csr")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        u, s, z, y, it = _ipm_loop(
            H, f, G, Gt, h, E, Et, e, u, s, z, y, I_r, I_q, rho, scale, m, q, r
        )
    rp = G @ u + s - h
    if (
        np.abs(rp).max(initial=0.0) > 1e-6 * scale
        or np.abs(E @ u - e).max(initial=0.0) > 1e-6 * scale
    ):
        raise LsqlinError("no feasible point: the constraints are inconsistent")
    return u, s, z, y, it


def _ipm_loop(H, f, G, Gt, h, E, Et, e, u, s, z, y, I_r, I_q, rho, scale, m, q, r):
    import scipy.sparse as sp

    best = None
    for it in range(1, _MAXIT + 1):
        rd = H @ u + f + Gt @ z + Et @ y
        rp = G @ u + s - h
        re = E @ u - e
        mu = float(np.sum(s * z)) / m if m else 0.0
        err = (
            max(
                np.abs(rd).max(initial=0.0),
                np.abs(rp).max(initial=0.0),
                np.abs(re).max(initial=0.0),
            )
            / scale
        )
        if err < _TOL and mu < _TOL * scale:
            return u, s, z, y, it - 1
        if best is None or err + mu / scale < best[0]:
            best = (err + mu / scale, u.copy(), s.copy(), z.copy(), y.copy(), it - 1)
        D = z / s
        lu = None
        for reg in (rho, 1e4 * rho, 1e8 * rho):
            M = (H + Gt @ sp.diags(D) @ G + reg * I_r).tocsr()
            try:
                if q:
                    lu = _splu(sp.bmat([[M, Et], [E, -reg * I_q]], format="csc"), symmetric=False)
                else:
                    lu = _splu(M)
                break
            except RuntimeError:
                continue
        if lu is None:
            break

        res = (rd, rp, re, D, s, G, Gt, lu, q, r)
        du, ds, dz, dy = _ipm_step(s * z, *res)
        a = min(_max_step(s, ds), _max_step(z, dz))
        mu_aff = float(np.sum((s + a * ds) * (z + a * dz))) / m
        sigma = (mu_aff / mu) ** 3 if mu > 0 else 0.0
        du, ds, dz, dy = _ipm_step(s * z + ds * dz - sigma * mu, *res)
        tau = max(0.99, 1.0 - mu / scale)
        a = min(1.0, tau * min(_max_step(s, ds), _max_step(z, dz)))
        u = u + a * du
        s = s + a * ds
        z = z + a * dz
        y = y + a * dy
        if not (np.all(np.isfinite(u)) and np.all(s > 0) and np.all(z > 0)):
            break
    return best[1:]


def _ipm_step(rc, rd, rp, re, D, s, G, Gt, lu, q, r):
    """the Newton step for the complementarity target rc (the factorization of the iteration)"""
    rhs = -rd - Gt @ (D * rp - rc / s)
    if q:
        sol = lu.solve(np.concatenate([rhs, -re]))
        du, dy = sol[:r], sol[r:]
    else:
        du, dy = lu.solve(rhs), np.zeros(0)
    dz = D * (G @ du + rp) - rc / s
    ds = -rp - G @ du
    return du, ds, dz, dy


def _max_step(v, dv):
    neg = dv < 0
    if not np.any(neg):
        return np.inf
    return float(np.min(-v[neg] / dv[neg]))


def _polish(C, d, G, h, E, e, u, active, z, y, refine=6, rounds=30):
    """the optimum on the interior-point solution's active set, by proximal iterative
    refinement of the KKT system of its rows (active inequalities as equalities, the
    equalities), then a primal active-set correction: rows with a negative multiplier leave,
    violated rows join. The interior-point solution if that does not end optimal."""
    import scipy.sparse as sp

    r = C.shape[1]
    H = (C.T @ C).tocsr()
    f = -(C.T @ d)
    scale = max(1.0, abs(H).max() if H.nnz else 1.0)
    rho = 1e-9 * scale
    tol = 1e-9 * (1.0 + np.abs(h).max(initial=0.0))
    act = active.copy()
    zw = z.copy()
    for _ in range(rounds):
        B = sp.vstack([G[act], E]).tocsr()
        c = np.concatenate([h[act], e])
        w = np.concatenate([np.maximum(zw[act], 0.0), y])
        k = B.shape[0]
        K = sp.bmat([[H + rho * sp.identity(r), B.T], [B, -rho * sp.identity(k)]], format="csc")
        try:
            lu = _splu(K, symmetric=False)
        except RuntimeError:
            return u
        v = u.copy()
        for _ in range(refine):
            sol = lu.solve(np.concatenate([-(H @ v + f + B.T @ w), c - B @ v]))
            v = v + sol[:r]
            w = w + sol[r:]
        if not np.all(np.isfinite(v)):
            return u
        na = int(np.sum(act))
        lam = w[:na]
        slack = h - G @ v
        neg = lam < -tol
        viol = (slack < -tol) & ~act
        if not neg.any() and not viol.any():
            return v
        idx = np.nonzero(act)[0]
        act[idx[neg]] = False
        zw[idx] = np.where(neg, 0.0, lam)
        if viol.any():
            # the most violated row joins (with a multiplier of 0 to start)
            j = int(np.argmin(np.where(viol, slack, np.inf)))
            act[j] = True
            zw[j] = 0.0
    return u


# ------------------------------------------------------------------------------------------------
# the member of the optimal set: free offsets at the analytic center, groups at mean zero
# ------------------------------------------------------------------------------------------------


def _face_center(P, x, active):
    """x on the optimal set; its free offsets moved to the analytic center. (x, active, free)"""
    if not P.structured:
        return _face_center_dense(P, x, active)
    s0 = P.b - P.A @ x
    r = _center_offsets(P, s0, active, need_feasible=False)
    if r is None:
        return x, active, P.lineality.copy()
    c, sup = r
    return x + c[P.ccomp], active, sup[P.ccomp]


def _center_offsets(P, s0, eq, need_feasible):
    """(c, super groups): offsets c (one per group of C and Aeq) maximizing
    sum log(s0_i - v_i (c_J - c_I)) over the rows between groups that are not tight on the whole
    set (those, and the rows in ``eq``, hold their groups together: the super groups, one label
    per group); one group per lineality group fixed at 0. With ``need_feasible``: None if no c
    holds every row."""
    import scipy.sparse as sp

    lab = P.ccomp
    k = lab.max() + 1
    ends = P.a_ends
    rowok = ends[:, 0] >= 0
    gi = np.where(rowok, lab[np.maximum(ends[:, 0], 0)], -1)
    gj = np.where(rowok, lab[np.maximum(ends[:, 1], 0)], -1)
    cross = rowok & (gi != gj)
    tol = 1e-9 * (1.0 + np.abs(P.b).max(initial=0.0))
    if need_feasible and np.any(s0[~cross] < -tol):
        return None
    eqc = eq & cross
    sup = _join(k, gi[eqc], gj[eqc])
    lin_of_group = np.zeros(k, np.int64)
    lin_of_group[lab] = P.lineality
    rows = np.nonzero(cross & ~eq)[0]
    while True:
        ns = sup.max() + 1
        # fixed per lineality group: the super group of its first group
        fixed = np.zeros(ns, bool)
        seen = set()
        for g in range(k):
            if lin_of_group[g] not in seen:
                seen.add(lin_of_group[g])
                fixed[sup[g]] = True
        var = np.nonzero(~fixed)[0]
        vidx = -np.ones(ns, np.int64)
        vidx[var] = np.arange(var.size)
        a = sup[gi[rows]]
        bb = sup[gj[rows]]
        same = a == bb
        if np.any(same):
            if need_feasible and np.any(s0[rows[same]] < -tol):
                return None
            rows, a, bb = rows[~same], a[~same], bb[~same]
        if var.size == 0:
            if need_feasible and np.any(s0[rows] < -tol):
                return None
            return np.zeros(k), sup
        v = P.a_val[rows]
        # slack_i(c) = s0_i - v_i (c_b - c_a)
        rr = np.arange(rows.size)
        ia, ib = vidx[a], vidx[bb]
        data = np.concatenate([-v[ia >= 0], v[ib >= 0]])
        ri = np.concatenate([rr[ia >= 0], rr[ib >= 0]])
        ci = np.concatenate([ia[ia >= 0], ib[ib >= 0]])
        G = sp.csr_matrix((data, (ri, ci)), shape=(rows.size, var.size))
        st = _strict_point(G, s0[rows])
        if st is None:
            return None
        t, cv, tight = st
        if need_feasible and t < -tol:
            return None
        if t > tol:
            cv = _analytic_center(G, s0[rows], cv)
            c_sup = np.zeros(ns)
            c_sup[var] = cv
            return c_sup[sup], sup
        # rows tight on the whole set: they hold their groups together
        if not np.any(tight):
            tight = s0[rows] - G @ cv <= tol
        tr = rows[tight]
        sup = _join(ns, sup[gi[tr]], sup[gj[tr]])[sup]
        rows = rows[~tight]


def _join(k, a, b):
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    G = sp.csr_matrix((np.ones(a.size), (a, b)), shape=(k, k))
    _, lab = connected_components(G, directed=False)
    return lab


def _strict_point(G, s0):
    """max t such that G c + t <= s0 (t <= 1): (t, c, rows with a positive dual) or None"""
    from scipy.optimize import linprog

    m, k = G.shape
    if m == 0:
        return 1.0, np.zeros(k), np.zeros(0, bool)
    import scipy.sparse as sp

    Aub = sp.hstack([G, sp.csr_matrix(np.ones((m, 1)))]).tocsr()
    cost = np.zeros(k + 1)
    cost[-1] = -1.0
    res = linprog(
        cost, A_ub=Aub, b_ub=s0, bounds=[(None, None)] * k + [(None, 1.0)], method="highs"
    )
    if res.status != 0:
        return None
    t = float(res.x[-1])
    marg = np.asarray(res.ineqlin.marginals)
    return t, res.x[:k], marg < -1e-12


def _analytic_center(G, s0, c, maxit=100):
    """Newton's method for max sum log(s0 - G c) from a strictly feasible c (damped while the
    Newton decrement is large, full steps near the center)"""
    import scipy.sparse as sp

    for _ in range(maxit):
        sl = s0 - G @ c
        g = G.T @ (1.0 / sl)
        Hn = (G.T @ sp.diags(1.0 / sl**2) @ G).tocsc()
        try:
            dc = -_splu(Hn).solve(g)
        except RuntimeError:
            return c
        lam2 = float(-np.sum(g * dc))
        if not np.isfinite(lam2):
            break
        if lam2 < 1e-26:
            return c
        t = 1.0
        fo = -np.sum(np.log(sl))
        while t > 1e-14:
            sn = s0 - G @ (c + t * dc)
            if np.all(sn > 0) and (lam2 < 1e-6 or -np.sum(np.log(sn)) <= fo - 0.25 * t * lam2):
                break
            t *= 0.5
        c = c + t * dc
        if lam2 < 1e-18 and t == 1.0:
            return c
        if not np.all(np.abs(c) < 1e15):
            break
    else:
        return c
    raise LsqlinError("the optimal set is unbounded in a direction the objective leaves free")


def _center_groups(P, x):
    """mean zero in each lineality group (an exactly rounded mean)"""
    if not P.structured:
        return x
    lab = P.lineality
    x = x.copy()
    for g in range(lab.max() + 1):
        m = np.nonzero(lab == g)[0]
        x[m] = x[m] - math.fsum(x[m]) / m.size
    return x


def _face_center_dense(P, x, active):
    """the analytic center over x + N c for rows that are not differences (dense, small)"""
    N = P.N
    L = P.L
    fixed = P.lineality.copy()
    if N.shape[1] == 0:
        return x, active, fixed
    A = P.A.toarray()
    s0 = P.b - A @ x
    G = A @ N
    vary = np.abs(G).max(axis=1) > 1e-12
    E = []
    if L.shape[1]:
        E.append(L.T @ N)
    eqr = active & vary
    if np.any(eqr):
        E.append(G[eqr])
    Z = _null(np.vstack(E)) if E else np.eye(N.shape[1])
    rows = vary & ~active
    if Z.shape[1] == 0 or not np.any(rows):
        return x, active, fixed
    import scipy.sparse as sp

    Gz = sp.csr_matrix(G[rows] @ Z)
    st = _strict_point(Gz, s0[rows])
    if st is None or st[0] <= 1e-9:
        return x, active, fixed
    w = _analytic_center(Gz, s0[rows], st[1])
    x = x + N @ (Z @ w)
    # variables the free directions move are their own group (marked free)
    moved = np.abs(N @ Z).max(axis=1) > 1e-12
    sup = fixed * 2
    sup[moved] = fixed[moved] * 2 + 1
    return x, active, sup
