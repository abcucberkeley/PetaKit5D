"""MATLAB's ``fminsearch``: the Nelder-Mead simplex method as MATLAB documents it.

Written from the published algorithm (J. C. Lagarias, J. A. Reeds, M. H. Wright, P. E. Wright,
"Convergence properties of the Nelder-Mead simplex method in low dimensions", SIAM J. Optim.
9(1), 1998) with its coefficients rho = 1, chi = 2, gamma = 1/2, sigma = 1/2, the initial simplex
MATLAB's documentation gives (each coordinate in turn changed by 5 %, a zero one set to
0.00025), its options (``TolX``, ``TolFun``, ``MaxIter``, ``MaxFunEvals``) and its stopping
test, and checked against MATLAB's iterates (black-box traces through an output function).
MathWorks' code is not read.

Points of the simplex are kept in the order of their function values (stable for ties); a step
evaluates the reflection, then the expansion or a contraction, or shrinks the simplex toward the
best point.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

__all__ = ["FminsearchResult", "fminsearch"]

_RHO, _CHI, _PSI, _SIGMA = 1.0, 2.0, 0.5, 0.5
_USUAL_DELTA, _ZERO_DELTA = 0.05, 0.00025


@dataclass
class FminsearchResult:
    """``[x, fval, exitflag, output]``: 1 converged, 0 stopped at ``max_iter`` or
    ``max_fun_evals``."""

    x: np.ndarray
    fval: float
    exitflag: int
    iterations: int
    func_count: int


def _eps(v: float) -> float:
    return float(np.spacing(abs(np.float64(v))))


def fminsearch(
    fun: Callable[[np.ndarray], float],
    x0: Sequence[float],
    *,
    tol_x: float = 1e-4,
    tol_fun: float = 1e-4,
    max_iter: float | None = None,
    max_fun_evals: float | None = None,
    output_fcn: Callable[[np.ndarray, float, int, int, str, str], None] | None = None,
) -> FminsearchResult:
    """The minimum of ``fun`` near ``x0`` (``fminsearch(fun, x0, optimset(...))``); ``max_iter``
    and ``max_fun_evals`` default to 200 times the number of variables. ``output_fcn(x, fval,
    iteration, funccount, procedure, state)`` is called as MATLAB calls an output function
    ('init', then 'iter' per iteration, 'done')."""
    x0 = np.asarray(x0, dtype=np.float64).ravel()
    n = x0.size
    maxiter = 200 * n if max_iter is None else max_iter
    maxfun = 200 * n if max_fun_evals is None else max_fun_evals

    def call(x: np.ndarray) -> float:
        return float(fun(x.copy()))

    v = np.empty((n, n + 1))
    fv = np.empty(n + 1)
    v[:, 0] = x0
    fv[0] = call(x0)
    func_count = 1
    iterations = 0
    if output_fcn is not None:
        output_fcn(x0.copy(), fv[0], iterations, func_count, "", "init")
        output_fcn(x0.copy(), fv[0], iterations, func_count, "", "iter")
    for j in range(n):
        y = x0.copy()
        y[j] = (1 + _USUAL_DELTA) * y[j] if y[j] != 0 else _ZERO_DELTA
        v[:, j + 1] = y
        fv[j + 1] = call(y)
    func_count = n + 1
    order = np.argsort(fv, kind="stable")
    v, fv = v[:, order], fv[order]
    iterations = 1
    procedure = "initial simplex"
    if output_fcn is not None:
        output_fcn(v[:, 0].copy(), fv[0], iterations, func_count, procedure, "iter")

    exitflag = 1
    while func_count < maxfun and iterations < maxiter:
        if np.max(np.abs(fv[0] - fv[1:])) <= max(tol_fun, 10 * _eps(fv[0])) and np.max(
            np.abs(v[:, 1:] - v[:, [0]])
        ) <= max(tol_x, 10 * _eps(np.max(v[:, 0]))):
            break
        # the centroid of the n best points (summed in index order)
        xbar = v[:, 0].copy()
        for j in range(1, n):
            xbar = xbar + v[:, j]
        xbar = xbar / n
        worst = v[:, n]
        xr = (1 + _RHO) * xbar - _RHO * worst
        fxr = call(xr)
        func_count += 1
        if fxr < fv[0]:
            xe = (1 + _RHO * _CHI) * xbar - _RHO * _CHI * worst
            fxe = call(xe)
            func_count += 1
            if fxe < fxr:
                v[:, n], fv[n], procedure = xe, fxe, "expand"
            else:
                v[:, n], fv[n], procedure = xr, fxr, "reflect"
        elif fxr < fv[n - 1]:
            v[:, n], fv[n], procedure = xr, fxr, "reflect"
        else:
            shrink = False
            if fxr < fv[n]:
                xc = (1 + _PSI * _RHO) * xbar - _PSI * _RHO * worst
                fxc = call(xc)
                func_count += 1
                if fxc <= fxr:
                    v[:, n], fv[n], procedure = xc, fxc, "contract outside"
                else:
                    shrink = True
            else:
                xcc = (1 - _PSI) * xbar + _PSI * worst
                fxcc = call(xcc)
                func_count += 1
                if fxcc < fv[n]:
                    v[:, n], fv[n], procedure = xcc, fxcc, "contract inside"
                else:
                    shrink = True
            if shrink:
                for j in range(1, n + 1):
                    v[:, j] = v[:, 0] + _SIGMA * (v[:, j] - v[:, 0])
                    fv[j] = call(v[:, j])
                func_count += n
                procedure = "shrink"
        order = np.argsort(fv, kind="stable")
        v, fv = v[:, order], fv[order]
        iterations += 1
        if output_fcn is not None:
            output_fcn(v[:, 0].copy(), fv[0], iterations, func_count, procedure, "iter")
    else:
        exitflag = 0
    if output_fcn is not None:
        output_fcn(v[:, 0].copy(), fv[0], iterations, func_count, "", "done")
    return FminsearchResult(v[:, 0].copy(), float(fv[0]), exitflag, iterations, func_count)
