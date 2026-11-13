"""MATLAB's ``stitch_assignment_lsqlin`` (the solver of the grid, global and group shift
assignment since d67e1b7) and ``lsqlin`` with its default (interior-point) algorithm: the minimum of
``||C x - d||^2`` subject to ``A x <= b`` and ``Aeq x = beq``, whose rows are weighted differences
of the shifts of two tiles.

``assignment(C, d, A, b, Aeq, beq)`` is ``stitch_assignment_lsqlin``, step by step:

1. Groups: of the rows of C and Aeq (objective groups: their offsets are free in the objective)
   and of all rows (lineality groups: a constant added to one changes nothing), numbered by their
   first tile.
2. Without rows in C (every pair below the correlation threshold), the stage positions (x = 0).
   Without rows in Aeq, the least squares of C (the first tile of each objective group fixed),
   the groups placed (3.), kept if every bound holds to 1e-9 (1 + max|b|).
3. Placement: the offsets of the objective groups (the group of each lineality group's first
   tile held) with the least sum of squares of the rows of A between groups, unconstrained, or
   within those rows' bounds where the unconstrained offsets break one. So a free tile gets the
   mean of its neighbors' shifts, a free group the mean over its pairs of the neighbor's shift
   less its own tile's.
4. Otherwise (equalities, a broken bound, or no offsets within the bounds): the exact optimum of
   the whole problem (the first tile of each lineality group fixed), then 3. again. MATLAB solves
   it with lsqlin refined on its active set, here an interior-point method refined the same way:
   the optimum's values within each objective group are unique, and 3. sets the offsets.
5. Mean zero in each lineality group. The callers round with ``round(x + 1e-9)``
   (``ROUND_OFFSET``) and subtract tile 1's value.

Both solve the same linear systems (MATLAB's sparse Cholesky and SuperLU here differ in the last
bits), so the rounded shifts can differ only for a value within about 1e-12 of a rounding
boundary, ``k + 0.5 - 1e-9`` (``near_round_boundary``).

``lsqlin``/``solve`` give MATLAB's ``lsqlin`` itself: the same exact optimum and member (rows of
differences; other rows by dense null spaces, small problems), with flags for where MATLAB's own
lsqlin differs, as it ends where its iterations stop: values near a half integer (its dense solver
drifts from mean zero by rounding noise that changes with the thread count) and values its path
decides (``free``: objective groups other than the largest of a lineality group).

``free_group_offsets`` is the group method's step 3 (``place_free_groups``).

Results never depend on the thread count (sparse LU and sparse products only).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "assignment",
    "free_group_offsets",
    "lsqlin",
    "solve",
    "LsqlinResult",
    "LsqlinError",
    "near_half_integer",
    "near_round_boundary",
    "ROUND_OFFSET",
    "ROUND_TOL",
    "HALF_TOL",
    "HALF_REL",
    "FREE_TOL",
]

# the callers round with round(x + ROUND_OFFSET), so that ties of the problem (exact shifts on
# half integers) round alike; the boundary moves to k + 0.5 - ROUND_OFFSET, and a value within
# ROUND_TOL of it is flagged (MATLAB's solution and this one agree to about 1e-12)
ROUND_OFFSET = 1e-9
ROUND_TOL = 1e-10

# for MATLAB's plain lsqlin: values within this of a half integer are flagged, by form of its
# solver: the sparse solution agrees to about 3e-5 once its mean is removed; the dense one drifts
# from mean zero by rounding noise that changes with MATLAB's thread count (up to 0.017 on
# generated grids, 0.065 on the demo datasets, 0.19 with weights 0.01), so the window is the larger
# of HALF_TOL and HALF_REL times the largest |x| of the lineality group
HALF_TOL = {"dense": 0.02, "sparse": 0.005}
HALF_REL = {"dense": 1e-3, "sparse": 0.0}
# how far MATLAB's plain lsqlin leaves free values from ours (a few voxels); mean zero passes their
# mean on to the other values of their group
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
    """``solve``'s solution and where MATLAB's plain lsqlin can differ from it.

    ``x``: the solution (``assignment``'s for rows of differences). ``free``: values in an
    objective group other than the largest of their lineality group (MATLAB's lsqlin leaves them
    where its iterations stop), and, for the sparse form, every value when there is more than one
    lineality direction (its offsets per group are arbitrary). ``near_half``: values within
    ``tol`` of a half integer; ``tol`` is the window of the value's lineality group (the larger of
    ``half_tol`` and ``HALF_REL[form]`` times the group's largest |x|), widened where the group has
    free values (the window times the group's size over its main part's, plus ``FREE_TOL`` times
    its free values over its main part's). ``active``: rows of ``A`` that hold with equality.
    ``groups``: the lineality group of each variable (0-based, numbered by first variable).
    ``iterations``: interior-point iterations (0 when none was needed)."""

    x: np.ndarray
    near_half: np.ndarray
    free: np.ndarray
    active: np.ndarray
    groups: np.ndarray
    iterations: int = 0
    form: str = "dense"
    tol: np.ndarray = field(default_factory=lambda: np.zeros(0))
    info: dict = field(default_factory=dict)


def assignment(C, d, A=None, b=None, Aeq=None, beq=None, *, info=None) -> np.ndarray:
    """``x = stitch_assignment_lsqlin(C, d, A, b, Aeq, beq)`` (see the module): float64 (n,).
    ``d``, ``b`` and ``beq`` may have columns, one system each, as the callers' three axes (which
    share C, A and Aeq): x has as many, each column bit for bit what a call with that column gives,
    and the work on the matrices is shared. A dict ``info`` receives how it was solved (a list per
    column with columns): ``path`` (``'stage'``: no rows in C; ``'ls'``: the least squares of C;
    ``'fallback'``: the whole problem; ``'empty'``: no variables), ``bounded`` (a bound held the
    last placement of the groups) and ``iterations`` (interior-point iterations)."""
    S = _System(C, d, A, b, Aeq, beq)
    x = np.zeros((0, S.k)) if S.n == 0 else _assign(S)
    if info is not None:
        for key in ("path", "bounded", "iterations"):
            v = getattr(S, key)
            info[key] = v if S.columns else v[0]
    return x if S.columns else x[:, 0]


def near_round_boundary(x, tol: float = ROUND_TOL) -> np.ndarray:
    """True where ``round(x + ROUND_OFFSET)`` is within ``tol`` of changing (``x`` within ``tol``
    of ``k + 0.5 - ROUND_OFFSET``)."""
    y = np.asarray(x, dtype=np.float64) + ROUND_OFFSET
    return np.abs(y - np.floor(y) - 0.5) < tol


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
    it = 0
    sup = None
    extra = np.zeros(P.n, bool)
    if P.p == 0 and P.q == 0:
        x = _unconstrained(P)
    elif P.structured:
        S = _System(P.C, P.d, P.A, P.b, P.Aeq, P.beq)
        x = _assign(S)[:, 0]
        sup, it = S.grp, S.iterations[0]
    else:
        x, extra, it = _solve_generic(P)
    tol0 = 1e-9 * (1.0 + np.abs(P.b).max(initial=0.0))
    active = (P.b - P.A @ x) <= tol0
    free, tol = _flags(P, x, sup, form, half_tol)
    free |= extra
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
    """free values (outside the main objective group of their lineality group, or all for the
    sparse form with several lineality directions) and the rounding window of each value: the
    group's (the larger of ``half_tol`` and HALF_REL times its largest |x|) scaled to its main
    part, plus FREE_TOL times its share of free values"""
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
# stitch_assignment_lsqlin
# ------------------------------------------------------------------------------------------------


class _System:
    """the system as stitch_assignment_lsqlin takes it: n is the most columns of C, A and Aeq (a
    matrix without rows takes n); each row with entries has exactly two. The right-hand sides are
    columns (D, B, BEQ: one system each), with the tolerance of each (``tol``)."""

    def __init__(self, C, d, A, b, Aeq, beq):
        import scipy.sparse as sp
        from scipy.sparse.csgraph import connected_components

        names = ("C", "A", "Aeq")
        mats = [_csr(M, None, name) for M, name in zip((C, A, Aeq), names, strict=True)]
        n = max(M.shape[1] for M in mats)
        for i, name in enumerate(names):
            if mats[i].shape[1] != n:
                if mats[i].shape[0]:
                    raise LsqlinError(f"{name} must have {n} columns")
                mats[i] = sp.csr_matrix((0, n))
        self.C, self.A, self.Aeq = mats
        self.n = n
        k = max([np.shape(v)[1] for v in (d, b, beq) if v is not None and np.ndim(v) == 2] + [1])
        self.k, self.columns = k, k > 1
        self.D = _rhs(d, self.C.shape[0], k, "d")
        self.B = _rhs(b, self.A.shape[0], k, "b")
        self.BEQ = _rhs(beq, self.Aeq.shape[0], k, "beq")
        self.iterations = [0] * k
        self.path = ["empty"] * k
        self.bounded = [False] * k
        if n == 0:
            return
        self.ec, _ = _row_ends(self.C)
        self.ea, self.ra = _row_ends(self.A)
        self.ee, _ = _row_ends(self.Aeq)
        self.grp = _labels(n, np.concatenate([self.ec, self.ee]), connected_components, sp)
        self.lin = _labels(n, np.concatenate([self.ec, self.ea, self.ee]), connected_components, sp)
        self.tol = 1e-9 * (1.0 + np.abs(self.B).max(axis=0, initial=0.0))
        self._place = False

    def placement(self):
        """place_groups' matrices, shared by the columns: (rows of A between objective groups, G,
        the factorization of G' G, G' G, the groups that move, the number of groups), or None if
        there is nothing to place"""
        import scipy.sparse as sp

        if self._place is not False:
            return self._place
        grp = self.grp
        cross = grp[self.ea[:, 0]] != grp[self.ea[:, 1]]
        kg = int(grp.max()) + 1
        mov = np.ones(kg, bool)
        mov[grp[_first_members(self.lin)]] = False
        self._place = None
        if cross.any() and mov.any():
            rows = self.ra[cross]
            Pm = sp.csr_matrix((np.ones(self.n), (np.arange(self.n), grp)), shape=(self.n, kg))
            G = (self.A[rows] @ Pm[:, np.nonzero(mov)[0]]).tocsr()
            H = (G.T @ G).tocsc()
            self._place = (rows, G, _splu(H), H, mov, kg)
        return self._place


def _rhs(v, m, k, name):
    """a right-hand side as (m, k): a vector (or one column) for every column"""
    v = np.zeros((0, 1)) if v is None else np.asarray(v, dtype=np.float64)
    if v.ndim < 2:
        v = v.reshape(-1, 1)
    if v.size == 0 and m == 0:
        return np.zeros((0, k))
    if v.ndim != 2 or v.shape[0] != m or v.shape[1] not in (1, k):
        raise LsqlinError(f"{name} must have {m} rows and 1 or {k} columns")
    return np.repeat(v, k, axis=1) if v.shape[1] != k else v.copy()


def _row_ends(M):
    """the two columns (ascending) of each row with entries, and those rows' indices; a row of one
    or of more than two entries is MATLAB's error"""
    nnz = np.diff(M.indptr)
    if np.any((nnz != 0) & (nnz != 2)):
        raise LsqlinError(
            "Each row must be a difference of the shifts of two tiles.",
            "stitch_assignment_lsqlin:rowError",
        )
    r = np.nonzero(nnz == 2)[0]
    st = M.indptr[r]
    c0, c1 = M.indices[st], M.indices[st + 1]
    ends = np.column_stack([np.minimum(c0, c1), np.maximum(c0, c1)]).astype(np.int64)
    return ends.reshape(-1, 2), r.astype(np.int64)


def _first_members(lab):
    """the first (lowest) index of each group, in order of the groups"""
    first = np.full(int(lab.max()) + 1, lab.size)
    np.minimum.at(first, lab, np.arange(lab.size))
    return first


def _assign(S):
    """stitch_assignment_lsqlin for each column: (n, k)"""
    n, k = S.n, S.k
    x = np.zeros((n, k))
    ok = np.zeros(k, bool)
    if S.ec.shape[0] == 0:
        # no rows in C: the stage positions
        ok[:] = True
        S.path = ["stage"] * k
    elif S.ee.shape[0] == 0:
        # the optimum of C alone (the first tile of each objective group fixed)
        S.path = ["ls"] * k
        keep = np.ones(n, bool)
        keep[_first_members(S.grp)] = False
        x[keep] = _normal_solve(S.C[:, keep], S.D)
        x, ok = _place_groups(S, x, range(k))
        for j in range(k):
            ok[j] = ok[j] and bool(np.all(S.A @ x[:, j] <= S.B[:, j] + S.tol[j]))
    if not ok.all():
        # the whole problem, refined on its active set (the first tile of each lineality group
        # fixed), then the groups placed again
        keep = np.ones(n, bool)
        keep[_first_members(S.lin)] = False
        Ck, Ak, Ek = S.C[:, keep], S.A[:, keep], S.Aeq[:, keep]
        for j in np.nonzero(~ok)[0]:
            S.path[j], S.bounded[j] = "fallback", False
            xj = np.zeros((n, 1))
            xj[keep, 0], it = _optimum(Ck, S.D[:, j], Ak, S.B[:, j], Ek, S.BEQ[:, j])
            S.iterations[j] += it
            x[:, j] = _place_groups(S, xj, [j])[0][:, 0]
    # mean zero in each lineality group (sums in index order, as accumarray's)
    cnt = np.bincount(S.lin)
    for j in range(k):
        x[:, j] -= (np.bincount(S.lin, weights=x[:, j]) / cnt)[S.lin]
    return x


def _normal_solve(C, d):
    """(C' C) \\ (C' d) for C' C positive definite (a sparse LU and a step of refinement); d may
    have columns (one factorization, each column solved alone)"""
    H = (C.T @ C).tocsc()
    lu = _splu(H)
    if d.ndim == 1:
        return _refined(lu, H, C.T @ d)
    return np.column_stack([_refined(lu, H, C.T @ d[:, j]) for j in range(d.shape[1])])


def _refined(lu, H, g):
    v = lu.solve(g)
    return v + lu.solve(g - H @ v)


def _optimum(C, d, G, h, E, e):
    """the exact optimum of min ||C u - d||^2, G u <= h, E u = e (no lineality left): an interior
    point refined on its active set, or the equalities' KKT system alone. (u, iterations)"""
    if G.shape[0] == 0:
        u = _eq_ls(C, d, E, e)
        if u is None:
            raise LsqlinError("no feasible point: the equalities are inconsistent")
        return u, 0
    u, _, z, y, it = _ipm(C, d, G, h, E, e)
    s = h - G @ u
    act = z > np.maximum(s, 1e-9 * (1.0 + np.abs(z).max(initial=0.0)))
    return _polish(C, d, G, h, E, e, u, act, z, y), it


def _place_groups(S, x, cols):
    """stitch_assignment_lsqlin's place_groups for the columns ``cols`` of the system (x: one
    column each): the offsets of the objective groups (the group of each lineality group's first
    tile held) with the least sum of squares of the rows of A between groups, within those rows'
    bounds. (x, ok); ok is False for a column if no offsets hold its bounds (x left as it was)."""
    import scipy.sparse as sp

    ok = np.ones(len(cols), bool)
    P = S.placement()
    if P is None:
        return x, ok
    rows, G, lu, H, mov, kg = P
    Ar = S.A[rows]
    x = x.copy()
    for i, j in enumerate(cols):
        S.bounded[j] = False
        r0 = Ar @ x[:, i]
        s0 = S.B[rows, j] - r0
        c = -_refined(lu, H, G.T @ r0)
        if np.any(G @ c > s0 + S.tol[j]):
            # a bound holds them: the bounded least squares, refined on its active set
            try:
                c, it = _optimum(G, -r0, G, s0, sp.csr_matrix((0, G.shape[1])), np.zeros(0))
            except LsqlinError:
                ok[i] = False
                continue
            S.iterations[j] += it
            S.bounded[j] = True
        cm = np.zeros(kg)
        cm[mov] = c
        x[:, i] = x[:, i] + cm[S.grp]
    return x, ok


# ------------------------------------------------------------------------------------------------
# the group method's step 3
# ------------------------------------------------------------------------------------------------


def free_group_offsets(placed, gid, d_shift, d_shift_pre, overlap) -> np.ndarray:
    """stitch_global_group_assignment's place_free_groups: the offsets (ng, 3) of the groups
    without nodes (``placed`` False), which keep them at their stage positions relative to the
    tiles of the other groups they overlap: per axis, the least squares of the shift differences
    ``d(q) - d(p)`` of the overlapping pairs (p, q) of tiles of two groups, at least one of them
    free (tiles of free groups at ``d_shift_pre`` plus their group's offset, the others at
    ``d_shift``), without bounds. Of each set of free groups linked by pairs and without a pair to
    a placed group, the first keeps offset 0. ``gid``: the 0-based group of each tile;
    ``overlap``: the tiles' overlap matrix."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    placed = np.asarray(placed, bool).reshape(-1)
    gid = np.asarray(gid, np.int64).reshape(-1)
    d_shift = np.asarray(d_shift, np.float64)
    d_shift_pre = np.asarray(d_shift_pre, np.float64)
    ng = placed.size
    offset = np.zeros((ng, 3))
    if placed.all():
        return offset
    qq, pp = np.nonzero(np.asarray(overlap, bool).T)  # find(overlap_matrix): column-major
    use = (gid[pp] != gid[qq]) & ~(placed[gid[pp]] & placed[gid[qq]])
    p, q = pp[use], qq[use]
    free = np.nonzero(~placed)[0]
    k = free.size
    col = -np.ones(ng, np.int64)
    col[free] = np.arange(k)
    cp, cq = col[gid[p]], col[gid[q]]
    # the shift differences d(q) - d(p) = r0 + G * (offsets of the free groups)
    known = d_shift_pre.copy()
    pl = placed[gid]
    known[pl] = d_shift[pl]
    r0 = known[q] - known[p]
    rq, rp = np.nonzero(cq >= 0)[0], np.nonzero(cp >= 0)[0]
    G = sp.csr_matrix(
        (np.r_[np.ones(rq.size), -np.ones(rp.size)], (np.r_[rq, rp], np.r_[cq[rq], cp[rp]])),
        shape=(p.size, k),
    )
    # sets of free groups linked by pairs; those with a pair to a placed group are held by it, the
    # others keep their first group at offset 0
    ff = (cp >= 0) & (cq >= 0)
    comp = _labels(k, np.column_stack([cp[ff], cq[ff]]), connected_components, sp)
    held = np.zeros(int(comp.max()) + 1, bool)
    fp = (cp >= 0) ^ (cq >= 0)
    held[comp[np.where(cp[fp] >= 0, cp[fp], cq[fp])]] = True
    keep = np.ones(k, bool)
    keep[_first_members(comp)[~held]] = False
    c = np.zeros((k, 3))
    if keep.any():
        c[keep] = -_normal_solve(G[:, np.nonzero(keep)[0]], r0)
    offset[free] = c
    return offset


# ------------------------------------------------------------------------------------------------
# lsqlin: the problem
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

            allr = np.concatenate([cs[0], as_[0], es[0]], axis=0)
            self.lineality = _labels(n, allr, connected_components, sp)
            self.ndirections = int(self.lineality.max()) + 1
        else:
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


def _solve_generic(P):
    """rows that are not differences: the exact optimum without the lineality directions (the
    solution orthogonal to them). Directions that only inequalities hold are left where the
    interior point ends; the variables they move are returned as free. (x, free, iterations)"""
    T, expand = _reduction(P)
    u, it = _optimum(P.C @ T, P.d, P.A @ T, P.b, P.Aeq @ T, P.beq)
    x = expand(u)
    N, L = P.N, P.L
    free = np.zeros(P.n, bool)
    if N.shape[1] > L.shape[1]:
        Q = N @ _null(L.T @ N) if L.shape[1] else N
        free = np.abs(Q).max(axis=1) > 1e-12
    return x, free, it


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


def _reduction(P):
    """x = T u without the lineality directions (rows that are not differences): u spans their
    orthogonal complement"""
    import scipy.sparse as sp

    L = P.L
    Q = _null(L.T) if L.shape[1] else np.eye(P.n)
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
