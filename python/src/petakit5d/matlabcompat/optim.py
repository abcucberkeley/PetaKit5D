"""MATLAB's ``lsqnonlin`` with bounds, trust-region-reflective (its default algorithm for bounded
problems), as measured from MATLAB R2024b (task C05: 469 problems with every evaluation recorded,
the lifetime analysis' 69 fits among them; no MathWorks source read).

The method is Coleman and Li's interior trust-region approach for bound-constrained least squares
(SIAM J. Optim. 6, 1996) with the two-dimensional subspace of Branch, Coleman and Li (SIAM J. Sci.
Comput. 21, 1999), for f(x) = |F(x)|^2 / 2, gradient g = J'F, Gauss-Newton Hessian J'J:

- **Jacobian**: forward differences, ``delta_j = sqrt(eps)*max(|x_j|, 1)*sign0(x_j)`` (sign0(0)
  = 1), the opposite sign where ``x_j + delta_j`` leaves the box; ``J(:, j) = (F(x + delta_j e_j)
  - F(x)) / delta_j``. F and J are evaluated at every trial point, in coordinate order after F.
- **Scaling**: Coleman and Li's ``v``: ``x - ub`` (g < 0, finite ub), ``x - lb`` (g >= 0, finite
  lb), -1 or 1 otherwise; ``dv`` 1 where a bound is finite. In the scaled variables
  ``sh = s / |v|^(1/2)``: ``gh = |v|^(1/2) g``, ``Mh = |v|^(1/2) J'J |v|^(1/2) + diag(|g| dv)``.
- **Step**: the Newton step ``-Mh^-1 gh`` when Mh is positive definite and the step is inside the
  radius; otherwise the trust-region solution in span{gh, Newton step} (exact: secular equation to
  rounding), and the Cauchy point along ``-gh``. A step that would leave the box is cut back to
  ``theta = max(0.95, 1 - |v.*g|_inf)`` of the way to the first bound (a Newton step cut back so
  has the Cauchy point as a second candidate); among the candidates the one of least model value
  is taken. Coleman and Li's reflected step is not a candidate: with it, 3 of the 469 recorded
  fits came out further from MATLAB's.
- **Radius**: 10 at the start; ``|sh|/4`` after a trial with ratio <= 0.25 or no decrease; doubled
  after a ratio > 0.75 with ``|sh| >= 0.9 radius``; otherwise kept. The ratio is (actual change of f
  plus ``sh'diag(|g| dv)sh/2``) over the model change.
- **Stops** (exit flags): 1, first-order optimality ``|v.*g|_inf`` below ``TolFun`` (MATLAB's legacy
  ``TolFun`` is also the optimality tolerance), checked at every iterate; 3, an accepted step that
  changed the sum of squares by less than ``TolFun*(1 + ss)`` with ratio > 0.25 and (ratio > 0.75 or
  the step inside 0.9 of the radius); 2, ``|sh| < TolX`` after any trial; 0, the evaluation or
  iteration limit.

**Not reproduced**: MATLAB sometimes takes the dogleg point of the subspace (the rule was not
found), and its subproblem solution is sometimes off the exact one by up to 1e-3 of the step; its
products and solves (BLAS and LAPACK) round otherwise. Where the paths agree the results agree to
about 1e-9 (relative); a non-smooth objective (the lifetime analysis' mode fit) can turn a
different step into a different local minimum (``docs/parity/detection.md``, rows L1 to L4).
Starting values where the objective is undefined (NaN or Inf) raise MATLAB's
``optimlib:commonMsgs:UsrObjUndefAtX0``.
The arithmetic is written out (left-to-right sums, closed-form 2 x 2 eigenvalues, Gaussian
elimination), so the results do not depend on the BLAS, the CPU or the thread count.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["LsqResult", "LsqnonlinError", "lsqnonlin"]

_EPS = float(np.finfo(np.float64).eps)
_SQEPS = float(np.sqrt(_EPS))
_DELTA0 = 10.0


class LsqnonlinError(ValueError):
    """an error MATLAB's lsqnonlin raises, with its identifier"""

    def __init__(self, identifier: str, message: str):
        super().__init__(f"{message} ({identifier})")
        self.identifier = identifier


@dataclass
class LsqResult:
    """lsqnonlin's outputs: x, resnorm (|F|^2), residual, exitflag and MATLAB's output fields"""

    x: np.ndarray
    resnorm: float
    residual: np.ndarray
    exitflag: int
    iterations: int
    func_count: int
    first_order_opt: float
    step_size: float


def _dot(a: np.ndarray, b: np.ndarray) -> float:
    """sum of a.*b, left to right (the same bits on every platform)"""
    p = np.multiply(a, b)
    return float(np.cumsum(p)[-1]) if p.size else 0.0


def _norm(a: np.ndarray) -> float:
    """the 2-norm as sqrt of the left-to-right sum of squares (correctly rounded sqrt)"""
    return float(np.sqrt(_dot(a, a)))


def _mv(M: np.ndarray, x: np.ndarray) -> np.ndarray:
    """M @ x with left-to-right row sums"""
    return np.array([_dot(M[i], x) for i in range(M.shape[0])])


def _mm(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """A @ B with left-to-right sums"""
    out = np.empty((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            out[i, j] = _dot(A[i], B[:, j])
    return out


def _pd(M: np.ndarray) -> bool:
    """positive definite: a Cholesky factorization with positive pivots"""
    n = M.shape[0]
    L = np.zeros((n, n))
    for j in range(n):
        d = M[j, j] - _dot(L[j, :j], L[j, :j])
        if not d > 0:
            return False
        L[j, j] = float(np.sqrt(d))
        for i in range(j + 1, n):
            L[i, j] = (M[i, j] - _dot(L[i, :j], L[j, :j])) / L[j, j]
    return True


def _gram(J: np.ndarray, F: np.ndarray):
    n = J.shape[1]
    H = np.empty((n, n))
    g = np.empty(n)
    for i in range(n):
        g[i] = _dot(J[:, i], F)
        for k in range(i + 1):
            H[i, k] = H[k, i] = _dot(J[:, i], J[:, k])
    return H, g


def _solve(A: np.ndarray, b: np.ndarray):
    """Gaussian elimination with partial pivoting (None for a singular matrix)"""
    n = b.size
    M = np.column_stack([A.astype(np.float64), b.astype(np.float64)])
    for c in range(n):
        p = c + int(np.argmax(np.abs(M[c:, c])))
        if M[p, c] == 0 or not np.isfinite(M[p, c]):
            return None
        if p != c:
            M[[c, p]] = M[[p, c]]
        for r in range(c + 1, n):
            f = M[r, c] / M[c, c]
            M[r, c:] = M[r, c:] - f * M[c, c:]
    x = np.empty(n)
    for r in range(n - 1, -1, -1):
        s = M[r, n]
        for k in range(r + 1, n):
            s = s - M[r, k] * x[k]
        x[r] = s / M[r, r]
    return x if np.all(np.isfinite(x)) else None


def _eig2(M: np.ndarray):
    """eigenvalues (ascending) and eigenvectors of a symmetric 2 x 2 matrix, closed form"""
    a, b, d = M[0, 0], M[0, 1], M[1, 1]
    if b == 0:
        vals = np.array([a, d])
        vecs = np.eye(2)
        if a > d:
            vals, vecs = vals[::-1], vecs[:, ::-1]
        return vals, vecs
    h = 0.5 * (a - d)
    r = float(np.sqrt(h * h + b * b))
    m = 0.5 * (a + d)
    lo, hi = (m - r, m + r)
    # the eigenvector of the larger eigenvalue, from the better-conditioned row
    if h >= 0:
        v1 = np.array([h + r, b])
    else:
        v1 = np.array([b, r - h])
    v1 = v1 / float(np.sqrt(v1[0] * v1[0] + v1[1] * v1[1]))
    v0 = np.array([-v1[1], v1[0]])
    return np.array([lo, hi]), np.column_stack([v0, v1])


def _eigsym(M: np.ndarray):
    if M.shape[0] == 1:
        return np.array([M[0, 0]]), np.ones((1, 1))
    if M.shape[0] == 2:
        return _eig2(M)
    raise ValueError("the subspace problems have dimension 1 or 2")


def _trust(g: np.ndarray, H: np.ndarray, delta: float):
    """min g'w + w'Hw/2 over |w| <= delta (dimension 1 or 2): Newton steps on the secular equation
    1/|w(mu)| = 1/delta, kept in a bracket, to rounding"""
    lam, Q = _eigsym(H)
    a = _mv(Q.T.copy(), g)
    if lam[0] > 0:
        w = -_mv(Q, a / lam)
        if _norm(w) <= delta:
            return w
    lo = max(0.0, -lam[0])
    hi = lo + _norm(g) / delta + 1.0
    for _ in range(200):
        if _norm(a / (lam + hi)) <= delta:
            break
        hi = 2.0 * hi
    mu = hi
    for _ in range(200):
        d = lam + mu
        r = a / d
        n = _norm(r)
        if n == 0:
            break
        phi = 1.0 / n - 1.0 / delta
        if phi > 0:
            hi = mu
        else:
            lo = mu
        dphi = _dot(r, r / d) / n**3
        new = mu - phi / dphi if dphi > 0 else 0.5 * (lo + hi)
        if not lo < new < hi:
            new = 0.5 * (lo + hi)
        if abs(new - mu) <= 4 * _EPS * abs(mu) or hi - lo <= 4 * _EPS * hi:
            mu = new
            break
        mu = new
    return -_mv(Q, a / (lam + mu))


def _definev(x, g, lb, ub):
    v = np.where(
        g < 0, np.where(np.isfinite(ub), x - ub, -1.0), np.where(np.isfinite(lb), x - lb, 1.0)
    )
    dv = np.where(g < 0, np.isfinite(ub), np.isfinite(lb)).astype(np.float64)
    return v, dv


def _to_bound(x, s, lb, ub):
    """the fraction of s to the first bound it meets (inf if none) and that bound's variable"""
    with np.errstate(divide="ignore", invalid="ignore"):
        tl = np.where(s < 0, (lb - x) / s, np.inf)
        tu = np.where(s > 0, (ub - x) / s, np.inf)
    t = np.minimum(tl, tu)
    if not s.size:
        return np.inf, -1
    i = int(np.argmin(t))
    return float(t[i]), i


def _fdjac(fun, x, F, lb, ub, args):
    n = x.size
    J = np.empty((F.size, n))
    for j in range(n):
        d = _SQEPS * max(abs(x[j]), 1.0) * (1.0 if x[j] >= 0 else -1.0)
        if x[j] + d > ub[j] or x[j] + d < lb[j]:
            d = -d
        xp = x.copy()
        xp[j] = x[j] + d
        J[:, j] = (np.asarray(fun(xp, *args), dtype=np.float64).ravel() - F) / d
    return J


def _step(x, g, H, lb, ub, delta):
    """the trial step: returns s, sh, gh, Mh, the C diagonal"""
    v, dv = _definev(x, g, lb, ub)
    Dinv = np.sqrt(np.abs(v))
    C = np.abs(g) * dv
    gh = Dinv * g
    Mh = Dinv[:, None] * H * Dinv[None, :] + np.diag(C)
    theta = max(0.95, 1.0 - float(np.max(np.abs(v * g))))
    sN = _solve(Mh, -gh)
    cands = []
    if sN is not None and _pd(Mh) and _norm(sN) <= delta:
        cands.append(sN)
    else:
        # the subspace of gh and the Newton direction (Gram-Schmidt)
        ng = _norm(gh)
        basis = [gh / ng] if ng > 0 else []
        if sN is not None and basis:
            w2 = sN - basis[0] * _dot(basis[0], sN)
            nw = _norm(w2)
            if nw > 1e-14 * _norm(sN):
                basis.append(w2 / nw)
        if basis:
            Z = np.column_stack(basis)
            Zt = Z.T.copy()
            cands.append(_mv(Z, _trust(_mv(Zt, gh), _mm(Zt, _mm(Mh, Z)), delta)))
            d = -gh / ng
            curv = _dot(d, _mv(Mh, d))
            tc = min(delta, ng / curv) if curv > 0 else delta
            cands.append(tc * d)
    if not cands:
        return np.zeros_like(x), np.zeros_like(x), gh, Mh, C
    out = []
    hit = False
    for sh in cands:
        tb, _ = _to_bound(x, Dinv * sh, lb, ub)
        if tb <= 1.0:
            hit = True
            sh = sh * (theta * tb)
        out.append(sh)
    if hit and len(cands) == 1:
        # the Newton step met a bound: the Cauchy point (cut back the same way) is a candidate too
        ng = _norm(gh)
        if ng > 0:
            d = -gh / ng
            curv = _dot(d, _mv(Mh, d))
            sc = (min(delta, ng / curv) if curv > 0 else delta) * d
            tb, _ = _to_bound(x, Dinv * sc, lb, ub)
            out.append(sc * (theta * tb) if tb <= 1.0 else sc)
    q = [_dot(gh, sh) + 0.5 * _dot(sh, _mv(Mh, sh)) for sh in out]
    sh = out[int(np.argmin(q))]
    return Dinv * sh, sh, gh, Mh, C


def lsqnonlin(
    fun,
    x0,
    lb=None,
    ub=None,
    args=(),
    *,
    tol_fun: float = 1e-6,
    tol_x: float = 1e-6,
    tol_opt: float | None = None,
    max_fun_evals: int | None = None,
    max_iter: int = 400,
) -> LsqResult:
    """MATLAB's ``lsqnonlin(fun, x0, lb, ub, optimset('TolFun', tol_fun, 'TolX', tol_x,
    'MaxFunEvals', max_fun_evals, 'MaxIter', max_iter, 'Jacobian', 'off'), args{:})``
    (module docstring). ``fun(x, *args)`` returns the residual vector. ``tol_opt`` (MATLAB's
    OptimalityTolerance) defaults to ``tol_fun``, as optimset's legacy TolFun sets both;
    ``max_fun_evals`` defaults to MATLAB's 100 times the number of variables."""
    x = np.array(x0, dtype=np.float64).ravel()
    n = x.size
    lb = (
        np.full(n, -np.inf)
        if lb is None
        else np.broadcast_to(np.asarray(lb, dtype=np.float64).ravel(), (n,)).copy()
    )
    ub = (
        np.full(n, np.inf)
        if ub is None
        else np.broadcast_to(np.asarray(ub, dtype=np.float64).ravel(), (n,)).copy()
    )
    tol_opt = tol_fun if tol_opt is None else tol_opt
    max_fun_evals = 100 * n if max_fun_evals is None else int(max_fun_evals)
    F = np.asarray(fun(x, *args), dtype=np.float64).ravel()
    if not np.all(np.isfinite(F)):
        raise LsqnonlinError(
            "optimlib:commonMsgs:UsrObjUndefAtX0",
            "the objective function is undefined (NaN or Inf) at the initial point",
        )
    J = _fdjac(fun, x, F, lb, ub, args)
    nfev = 1 + n
    ss = _dot(F, F)
    delta = _DELTA0
    flag = 0
    it = 0
    nrmsh = 0.0
    while True:
        H, g = _gram(J, F)
        v, _ = _definev(x, g, lb, ub)
        opt = float(np.max(np.abs(v * g))) if n else 0.0
        if opt < tol_opt:
            flag = 1
            break
        if it >= max_iter or nfev >= max_fun_evals:
            break
        s, sh, gh, Mh, C = _step(x, g, H, lb, ub, delta)
        nrmsh = _norm(sh)
        xt = x + s
        Ft = np.asarray(fun(xt, *args), dtype=np.float64).ravel()
        Jt = _fdjac(fun, xt, Ft, lb, ub, args)
        nfev += 1 + n
        it += 1
        sst = _dot(Ft, Ft)
        qp = _dot(gh, sh) + 0.5 * _dot(sh, _mv(Mh, sh))
        ratio = (0.5 * sst - 0.5 * ss + 0.5 * _dot(sh, C * sh)) / qp if qp != 0 else 0.0
        inside = nrmsh < 0.9 * delta
        decrease = sst < ss
        if ratio <= 0.25 or not decrease:
            new_delta = 0.25 * nrmsh
        elif ratio > 0.75 and not inside:
            new_delta = 2.0 * delta
        else:
            new_delta = delta
        if decrease:
            change = ss - sst
            ss_old = ss
            x, F, J, ss = xt, Ft, Jt, sst
            if change < tol_fun * (1.0 + ss_old) and ratio > 0.25 and (ratio > 0.75 or inside):
                flag = 3
                H, g = _gram(J, F)
                v, _ = _definev(x, g, lb, ub)
                opt = float(np.max(np.abs(v * g)))
                break
        delta = new_delta
        if nrmsh < tol_x:
            flag = 2
            H, g = _gram(J, F)
            v, _ = _definev(x, g, lb, ub)
            opt = float(np.max(np.abs(v * g)))
            break
    return LsqResult(
        x=x,
        resnorm=ss,
        residual=F,
        exitflag=flag,
        iterations=it,
        func_count=nfev,
        first_order_opt=opt,
        step_size=nrmsh,
    )
