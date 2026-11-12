"""Tests of petakit5d.matlabcompat.lsq: MATLAB's lsqlin (interior-point) as the stitching shift
assignment uses it.

The contract is the rounded result: the oracle cases run MATLAB's own callers
(``stitch_global_grid_assignment``, dense up to 400 tiles and sparse above, and
``stitch_global_assignment``) and ``lsqlin`` itself on generated problems, and compare the rounded
shifts exactly, except the values the solver flags: near a half integer (rounding may differ by
one) and free (held only by inactive bounds; MATLAB's value is where its iterations stopped).
The flags and Python's values at the flagged places are inputs of each case, so a case's key
changes when they do, and MATLAB reports how its flagged values differ from Python's (at most 1
where a value is near a half integer). The callers' matrices are built here as the MATLAB code
builds them (the port of the callers is the stitching task's).

MATLAB tests whose intent is ported: tests/stitch/StitchShiftAssignmentTest.m
largeGridMatchesDenseSolver (test_large_grid_sparse_matches_dense; the class is the stitching
task's).
"""

import math
import os
import subprocess
import sys

import numpy as np
import pytest

from petakit5d.core.dtypes import mround
from petakit5d.matlabcompat import lsq
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.compare import Custom

# --- bit-stable inputs (integer arithmetic and + - * / only) -----------------------------------

_C1 = np.uint64(0xBF58476D1CE4E5B9)
_C2 = np.uint64(0x94D049BB133111EB)
_GOLD = np.uint64(0x9E3779B97F4A7C15)


class Rand:
    """splitmix64 streams: each draw uses its own sub-seed"""

    def __init__(self, seed: int):
        self.seed = seed
        self.k = 0

    def u64(self, n: int) -> np.ndarray:
        self.k += 1
        s = np.uint64((self.seed * 1000003 + self.k * 7919) & 0xFFFFFFFFFFFFFFFF)
        with np.errstate(over="ignore"):
            z = s + np.arange(1, n + 1, dtype=np.uint64) * _GOLD
            z = (z ^ (z >> np.uint64(30))) * _C1
            z = (z ^ (z >> np.uint64(27))) * _C2
            return z ^ (z >> np.uint64(31))

    def rand(self, n: int) -> np.ndarray:
        return (self.u64(n) >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)

    def randi(self, lo: int, hi: int, n: int) -> np.ndarray:
        return lo + np.floor(self.rand(n) * (hi - lo + 1))

    def randn(self, n: int) -> np.ndarray:
        """about normal: the sum of four uniforms, centered and scaled to variance 1"""
        u = self.rand(4 * n).reshape(4, n)
        return (((u[0] + u[1]) + u[2]) + u[3] - 2.0) * 1.7320508075688772


def grid_problem(
    seed,
    gx,
    gy,
    gz=1,
    *,
    filt=0.05,
    noise=0.7,
    off=(15, 15, 5),
    max_off=(100, 100, 30),
    overlap=(100, 100, 30),
    free_tiles=0,
    holes=0,
    perm=False,
    two_groups=False,
    isolated=0,
    thr=0.25,
    aw=(1.0, 0.1, 10.0),
    stitch2d=False,
    corners=False,
    return_true=False,
):
    """stitch_shift_assignment's inputs for a grid of tiles: pairs of neighbors (and of corner
    neighbors with ``corners``), noisy relative shifts, correlations (``filt`` of them below the
    threshold) and the bounds of the 'grid' method (with the true offsets if ``return_true``)"""
    r = Rand(seed)
    idx = [(x, y, z) for z in range(gz) for x in range(gx) for y in range(gy)]
    idx = np.array(idx, dtype=np.float64).reshape(-1, 3)
    if two_groups:
        other = idx.copy()
        other[:, 0] += gx + 2
        idx = np.vstack([idx, other])
    if holes:
        keep = np.ones(len(idx), bool)
        keep[r.randi(1, len(idx) - 1, holes).astype(np.int64)] = False
        idx = idx[keep]
    if isolated:
        far = np.array([[1000.0 + 3 * k, 1000.0, 0.0] for k in range(isolated)])
        idx = np.vstack([idx, far])
    nF = len(idx)
    if perm:
        idx = idx[np.argsort(r.rand(nF), kind="stable")]
    tileIdx = np.hstack([idx, np.ones((nF, 1))])
    xyz = idx * np.array([900.0, 900.0, 300.0])
    true_off = np.column_stack(
        [r.randi(-off[0], off[0], nF), r.randi(-off[1], off[1], nF), r.randi(-off[2], off[2], nF)]
    )
    if stitch2d:
        true_off[:, 2] = 0
    pairs = []
    for i in range(nF):
        dd = np.abs(idx[i + 1 :] - idx[i])
        ok = dd.sum(axis=1) == 1
        if corners:
            ok |= (dd.max(axis=1) == 1) & (dd[:, 2] == 0)
        pairs.extend((i, i + 1 + j) for j in np.nonzero(ok)[0])
    pairs = np.array(pairs, dtype=np.int64).reshape(-1, 2)
    ti, tj = pairs[:, 0], pairs[:, 1]
    nP = len(ti)
    rel = true_off[tj] - true_off[ti] + r.randn(3 * nP).reshape(3, nP).T * noise
    if stitch2d:
        rel[:, 2] = 0.0
    xc = 0.3 + 0.7 * r.rand(nP)
    xc[r.rand(nP) < filt] = 0.1
    if free_tiles:
        for t in r.randi(0, nF - 1, free_tiles).astype(np.int64):
            xc[(ti == t) | (tj == t)] = 0.1
    mo = np.array(max_off, dtype=np.float64)
    l = -np.ones((nP, 1)) * mo
    u = np.minimum(np.array(overlap, dtype=np.float64) - 3, mo) * np.ones((nP, 1))
    if stitch2d:
        u[:, 2] = 0.0
    order = (xyz[tj] >= xyz[ti]).astype(np.float64)
    mas = np.hstack([l * order + u * (order - 1), u * order + l * (order - 1)])
    one = np.column_stack([ti + 1, tj + 1]).astype(np.float64)
    p = {
        "nF": float(nF),
        "mx": np.column_stack([one, xc]),
        "ash": np.column_stack([one, rel]),
        "pairs": one,
        "mas": mas,
        "thr": float(thr),
        "tileIdx": tileIdx,
        "aw": np.array(aw, dtype=np.float64),
    }
    return (p, true_off) if return_true else p


# --- the callers' systems, as the MATLAB code builds them --------------------------------------


def _overlap(p):
    nF = int(p["nF"])
    ov = np.zeros((nF, nF), bool)
    pr = p["pairs"].astype(np.int64) - 1
    ov[pr[:, 0], pr[:, 1]] = True
    return ov


def grid_systems(p):
    """stitch_global_grid_assignment's (C, d, A, b) per axis and its form (sparse above 400)"""
    import scipy.sparse as sp

    nF = int(p["nF"])
    nj, ni = np.nonzero(_overlap(p).T)  # find(): column-major
    neq = ni.size
    mx, tI, aw = p["mx"], p["tileIdx"], p["aw"]
    filt = mx[:, 2] < p["thr"]
    mxf = mx[~filt]
    nP = mxf.shape[0]
    s = mxf[:, 0].astype(np.int64) - 1
    t = mxf[:, 1].astype(np.int64) - 1
    tdiff = np.abs(tI[s, :3] - tI[t, :3])
    w = np.zeros(nP)
    nb = tdiff.sum(axis=1) == 1
    w[nb] = aw[[1, 0, 2]][np.argmax(tdiff[nb], axis=1)]
    sw = np.sqrt(w)
    dw = sw[:, None] * p["ash"][~filt, 2:5]
    rows = np.arange(neq)
    if nF > 400:
        R = sp.csr_matrix(
            (np.r_[-np.ones(neq), np.ones(neq)], (np.r_[rows, rows], np.r_[ni, nj])),
            shape=(neq, nF),
        )
        r2 = np.arange(nP)
        C = sp.diags(sw) @ sp.csr_matrix(
            (np.r_[-np.ones(nP), np.ones(nP)], (np.r_[r2, r2], np.r_[s, t])), shape=(nP, nF)
        )
        A = sp.vstack([R, -R]).tocsr()
        form = "sparse"
    else:
        R = np.zeros((neq, nF))
        R[rows, ni] = -1
        R[rows, nj] = 1
        Rw = np.zeros((nP, nF))
        Rw[np.arange(nP), s] = -1
        Rw[np.arange(nP), t] = 1
        C = sw[:, None] * Rw
        A = np.vstack([R, -R])
        form = "dense"
    mas = p["mas"]
    return [(C, dw[:, i], A, np.r_[mas[:, 3 + i], -mas[:, i]]) for i in range(3)], form


def global_systems(p):
    """stitch_global_assignment's systems (weights: the correlations; always dense)"""
    nF = int(p["nF"])
    nj, ni = np.nonzero(_overlap(p).T)
    neq = ni.size
    R = np.zeros((neq, nF))
    R[np.arange(neq), ni] = -1
    R[np.arange(neq), nj] = 1
    mx = p["mx"]
    filt = mx[:, 2] < p["thr"]
    mxf = mx[~filt]
    nP = mxf.shape[0]
    w = mxf[:, 2]
    Rw = np.zeros((nP, nF))
    Rw[np.arange(nP), mxf[:, 0].astype(np.int64) - 1] = -1
    Rw[np.arange(nP), mxf[:, 1].astype(np.int64) - 1] = 1
    sw = np.sqrt(w)
    C = sw[:, None] * Rw
    dw = sw[:, None] * p["ash"][~filt, 2:5]
    A = np.vstack([R, -R])
    mas = p["mas"]
    return [(C, dw[:, i], A, np.r_[mas[:, 3 + i], -mas[:, i]]) for i in range(3)], "dense"


def assign(p, method):
    """the rounded shifts relative to tile 1, and the flags of each value (near, free), as the
    caller computes them: round, then subtract tile 1's"""
    systems, form = (grid_systems if method == "grid" else global_systems)(p)
    nF = int(p["nF"])
    x = np.zeros((nF, 3))
    near = np.zeros((nF, 3), bool)
    free = np.zeros((nF, 3), bool)
    for i, (C, d, A, b) in enumerate(systems):
        res = lsq.solve(C, d, A, b, form=form)
        xi = res.x
        if form == "sparse":
            xi = xi - math.fsum(xi) / xi.size
        x[:, i] = xi
        near[:, i] = res.near_half
        free[:, i] = res.free
    d_shift = mround(x)
    # + 0: a zero shift is +0 (its sign follows the solver's rounding noise near 0)
    d_shift = d_shift - d_shift[0] + 0.0
    # a value depends on its own rounding and on tile 1's
    fr = free | free[0]
    nh = (near | near[0]) & ~fr
    return d_shift, nh, fr


def flags_input(p, method):
    d, nh, fr = assign(p, method)
    m = nh | fr
    return {"near": nh, "free": fr, "pyv": d.ravel(order="F")[m.ravel(order="F")][:, None]}


def py_assign(p, f, method):
    d, _, _ = assign(p, method)
    m = f["near"] | f["free"]
    d[m] = 0.5
    s = np.array([[np.count_nonzero(f["near"]), 0, 0, np.count_nonzero(f["free"]), 0, 0]], float)
    return [d, s]


_SUMMARY = """
d = d + 0;
m = f.near | f.free;
dv = d(m) - f.pyv;
nm = f.near(m);
fm = f.free(m);
s = [nnz(nm), nnz(dv(nm)), max([0; abs(dv(nm))]), nnz(fm), nnz(dv(fm)), max([0; abs(dv(fm))])];
d(m) = 0.5;
out = {d, s};
"""

_GRID_CODE = (
    """
p = in{1}; f = in{2};
ov = false(p.nF); ov(sub2ind([p.nF, p.nF], p.pairs(:, 1), p.pairs(:, 2))) = true;
d = [];
evalc('d = stitch_global_grid_assignment(p.nF, p.mx, p.ash, ov, p.mas, p.thr, p.tileIdx, p.aw);');
"""
    + _SUMMARY
)

_GLOBAL_CODE = (
    """
p = in{1}; f = in{2};
ov = false(p.nF); ov(sub2ind([p.nF, p.nF], p.pairs(:, 1), p.pairs(:, 2))) = true;
d = [];
evalc('d = stitch_global_assignment(p.nF, p.mx, p.ash, ov, p.mas, p.thr);');
"""
    + _SUMMARY
)


def summary_check(ref, act):
    """out2: MATLAB's differences at the flagged values (counts must agree; near-half values
    differ by at most 1; free values are reported, not limited)"""
    r = np.asarray(ref.data, dtype=np.float64).ravel()
    a = np.asarray(act.data, dtype=np.float64).ravel()
    problems = []
    if r[0] != a[0] or r[3] != a[3]:
        problems.append(f"flag counts differ: MATLAB {r[[0, 3]]}, Python {a[[0, 3]]}")
    if r[2] > 1:
        problems.append(f"a value flagged near a half integer differs by {r[2]:g} (at most 1)")
    return problems


SUMMARY = Custom(summary_check, "flagged values: counts equal, near-half ones within 1")


def _assignment_case(cid, method, spec):
    seed, gx, gy, gz, kw = spec

    def prob():
        return grid_problem(seed, gx, gy, gz, **kw)

    p = prob()
    if not np.any(p["mx"][:, 2] >= p["thr"]):
        # every pair below the threshold: C is empty and lsqlin refuses it
        return Case(
            cid,
            code=_GRID_CODE if method == "grid" else _GLOBAL_CODE,
            args=[In("p"), {}],
            nargout=2,
            python=lambda p, f: py_assign(p, f, method),
            inputs={"p": Value(prob)},
            expect_error="optimlib:lsqlin:FirstTwoArgsEmpty",
            py_error=lsq.LsqlinError,
        )
    return Case(
        cid,
        code=_GRID_CODE if method == "grid" else _GLOBAL_CODE,
        args=[In("p"), In("f")],
        nargout=2,
        python=lambda p, f: py_assign(p, f, method),
        inputs={"p": Value(prob), "f": Value(lambda: flags_input(prob(), method))},
        compare={1: SUMMARY},
    )


def _grid_specs():
    """(id, (seed, gx, gy, gz, options)): sizes 2 to 1,600 tiles, every family"""
    sp = []
    s = 100

    def add(name, gx, gy, gz=1, **kw):
        nonlocal s
        s += 1
        sp.append((f"grid_{name}_{gx}x{gy}x{gz}_s{s}", (s, gx, gy, gz, kw)))

    small = [
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (3, 2),
        (2, 3),
        (3, 3),
        (4, 3),
        (4, 4),
        (5, 4),
        (5, 5),
        (6, 5),
        (6, 6),
        (7, 6),
        (8, 7),
        (9, 8),
        (10, 10),
    ]
    mid = [(12, 11), (14, 14), (17, 16), (20, 20)]
    for gx, gy in small + mid:
        add("conn", gx, gy, filt=0.0)
        add("filt", gx, gy, filt=0.15)
        add("tight", gx, gy, max_off=(6, 6, 3), off=(12, 12, 4))
    for gx, gy in small[4:] + mid:
        add("free", gx, gy, free_tiles=max(1, gx * gy // 12))
        add("bigfilt", gx, gy, filt=0.45)
        add("tightfree", gx, gy, filt=0.25, max_off=(8, 8, 3), free_tiles=2)
        add("noisy", gx, gy, noise=2.5, filt=0.1)
        add("weights", gx, gy, aw=(2.0, 0.5, 3.0), filt=0.1)
    for gx, gy in small[2:]:
        add("stitch2d", gx, gy, stitch2d=True, filt=0.2)
        add("thr", gx, gy, thr=0.6, filt=0.0)
    for gx, gy, gz in [(2, 2, 2), (3, 2, 2), (3, 3, 2), (4, 3, 3), (5, 4, 3), (6, 5, 4), (7, 7, 4)]:
        add("conn3d", gx, gy, gz, filt=0.0)
        add("filt3d", gx, gy, gz, filt=0.2)
        add("tight3d", gx, gy, gz, max_off=(6, 6, 2), off=(12, 12, 4))
        add("free3d", gx, gy, gz, free_tiles=2)
    for gx, gy in [(3, 3), (4, 4), (6, 5), (8, 8), (11, 10)]:
        add("holes", gx, gy, holes=3, perm=True, filt=0.1)
        add("two", gx, gy, two_groups=True, filt=0.1)
        add("isolated", gx, gy, isolated=2, filt=0.1)
        add("activeall", gx, gy, max_off=(2, 2, 1), off=(10, 10, 3))
    # the sparse form, above 400 tiles
    for gx, gy, gz in [(21, 21, 1), (25, 20, 1), (30, 30, 1), (40, 40, 1), (9, 9, 6)]:
        add("conn", gx, gy, gz, filt=0.0)
        add("filt", gx, gy, gz, filt=0.1)
        add("free", gx, gy, gz, free_tiles=6)
        add("tight", gx, gy, gz, max_off=(6, 6, 3), off=(12, 12, 4))
        add("bigfilt", gx, gy, gz, filt=0.4)
    add("two", 21, 11, two_groups=True, filt=0.05)
    add("holes", 25, 20, holes=10, perm=True, filt=0.1)
    # a second draw of the small and middle sizes
    for gx, gy in small[3:] + mid[:3]:
        add("conn", gx, gy, filt=0.0, noise=1.2)
        add("filt", gx, gy, filt=0.3)
        add("tight", gx, gy, max_off=(4, 4, 2), off=(10, 10, 4))
        add("free", gx, gy, free_tiles=max(1, gx * gy // 8), filt=0.05)
        add("bigfilt", gx, gy, filt=0.5)
        add("stitch2dfree", gx, gy, stitch2d=True, free_tiles=1, filt=0.1)
    for gx, gy, gz in [(2, 2, 3), (3, 3, 3), (4, 4, 2), (5, 5, 2), (6, 4, 3)]:
        add("bigfilt3d", gx, gy, gz, filt=0.45)
        add("tightfree3d", gx, gy, gz, filt=0.25, max_off=(8, 8, 3), free_tiles=2)
    # lines of tiles
    for n in (4, 7, 12, 20, 35):
        add("line", 1, n, filt=0.1)
        add("linetight", n, 1, max_off=(5, 5, 2), off=(12, 12, 4))
        add("linefree", 1, n, free_tiles=1)
    # sizes and families in turn
    fams = [
        ("mixconn", {"filt": 0.0}),
        ("mixfilt", {"filt": 0.2}),
        ("mixtight", {"max_off": (7, 7, 3), "off": (12, 12, 4)}),
        ("mixfree", {"free_tiles": 2, "filt": 0.1}),
        ("mixbig", {"filt": 0.4}),
        ("mixnoisy", {"noise": 1.8, "filt": 0.1}),
    ]
    sizes = [
        (2, 5),
        (3, 7),
        (4, 9),
        (5, 11),
        (6, 13),
        (7, 3),
        (8, 5),
        (9, 7),
        (10, 9),
        (11, 4),
        (13, 6),
        (15, 8),
        (16, 3),
        (18, 5),
        (19, 7),
        (6, 17),
        (5, 19),
        (3, 15),
        (12, 14),
        (14, 12),
        (7, 18),
        (17, 7),
        (9, 16),
        (16, 9),
        (4, 20),
        (20, 4),
        (8, 13),
        (13, 8),
        (10, 15),
        (15, 10),
    ]
    for k, (gx, gy) in enumerate(sizes):
        name, kw = fams[k % len(fams)]
        add(name, gx, gy, **kw)
    for gx, gy, gz in [(3, 3, 4), (4, 4, 4), (5, 5, 3), (6, 6, 2), (2, 8, 5)]:
        add("conn3d", gx, gy, gz, filt=0.0, noise=1.2)
        add("tightfree3d", gx, gy, gz, filt=0.2, max_off=(7, 7, 3), free_tiles=1)
    for gx, gy in [(22, 22), (35, 30)]:
        add("conn", gx, gy, filt=0.05)
        add("tightfree", gx, gy, filt=0.2, max_off=(8, 8, 3), free_tiles=4)
    return sp


def _global_specs():
    sp = []
    s = 5000
    for gx, gy in [(2, 2), (3, 2), (3, 3), (4, 4), (5, 5), (6, 5), (7, 7), (9, 9), (12, 10)]:
        for name, kw in (
            ("conn", {"filt": 0.0}),
            ("filt", {"filt": 0.2}),
            ("tight", {"max_off": (6, 6, 3), "off": (12, 12, 4)}),
            ("free", {"free_tiles": 2, "filt": 0.1}),
            ("bigfilt", {"filt": 0.45}),
        ):
            s += 1
            sp.append((f"global_{name}_{gx}x{gy}_s{s}", (s, gx, gy, 1, dict(kw, corners=True))))
    return sp


# --- lsqlin itself: generic rows, equalities, sparse C ------------------------------------------


def _matvec(M, v):
    """M @ v summed over the columns in order (no BLAS: the same bits everywhere)"""
    acc = np.zeros(M.shape[0])
    for j in range(M.shape[1]):
        acc = acc + M[:, j] * v[j]
    return acc


def generic_problem(seed, n, m, p, *, sparse_c=False, eq=0):
    """a full-rank C with rows that are not differences; inequalities that a point xf holds
    (with small slacks: they cut the unconstrained optimum), equalities that xf holds"""
    r = Rand(seed)
    C = r.randn(m * n).reshape(m, n)
    d = r.randn(m) * 10.0
    A = r.randn(p * n).reshape(p, n)
    xf = r.randn(n)
    b = _matvec(A, xf) + r.rand(p) * 0.5
    Aeq = r.randn(eq * n).reshape(eq, n) if eq else np.zeros((0, n))
    beq = _matvec(Aeq, xf) if eq else np.zeros(0)
    return {
        "C": C,
        "d": d[:, None],
        "A": A,
        "b": b[:, None],
        "Aeq": Aeq,
        "beq": beq[:, None],
        "sp": float(sparse_c),
        "ht": 0.0,
    }


def group_problem(seed, ng, *, tight=False):
    """step 2 of stitch_global_group_assignment: pairs of groups as pairs of nodes, the
    distances of their centroids as equalities of the nodes of each group, bounds around true
    group offsets (feasible), distances that pull past them (``tight``: active)"""
    r = Rand(seed)
    pairs = [
        (i, j) for i in range(ng) for j in range(i + 1, ng) if r.rand(1)[0] < 0.7 or j == i + 1
    ]
    nP = len(pairs)
    n = 2 * nP
    w = 0.01
    R = np.zeros((nP, n))
    R[np.arange(nP), np.arange(nP)] = -1
    R[np.arange(nP), nP + np.arange(nP)] = 1
    C = math.sqrt(w) * R
    cent = mround(r.randn(n) * 30.0 * 16.0) / 16.0
    off = mround(r.randn(ng) * 20.0 * 16.0) / 16.0
    gi = np.array([a for a, _ in pairs])
    gj = np.array([bb for _, bb in pairs])
    true = (off[gj] + cent[nP:]) - (off[gi] + cent[:nP])
    lim = 2.0 if tight else 40.0
    lo = true - lim - np.floor(r.rand(nP) * 4)
    hi = true + lim + np.floor(r.rand(nP) * 4)
    dist = true + mround(r.randn(nP) * 6.0 * 64.0) / 64.0
    d = math.sqrt(w) * dist
    A = np.vstack([R, -R])
    b = np.r_[hi, -lo]
    rows = []
    beq = []
    for g in range(ng):
        nodes = sorted(
            {k for k, (a, _) in enumerate(pairs) if a == g}
            | {nP + k for k, (_, bb) in enumerate(pairs) if bb == g}
        )
        for a in range(len(nodes)):
            for bb in range(a + 1, len(nodes)):
                row = np.zeros(n)
                row[nodes[a]] = -1
                row[nodes[bb]] = 1
                rows.append(row)
                beq.append(cent[nodes[bb]] - cent[nodes[a]])
    Aeq = np.array(rows).reshape(-1, n)
    return {
        "C": C,
        "d": d[:, None],
        "A": A,
        "b": b[:, None],
        "Aeq": Aeq,
        "beq": np.array(beq)[:, None],
        "sp": 0.0,
        # weights 0.01: MATLAB's constant drifts up to 0.19 from mean zero (measured), so the
        # window of rounding is wider
        "ht": 0.25,
    }


_LSQLIN_CODE = """
p = in{1}; f = in{2};
C = p.C; if p.sp, C = sparse(C); end
x = [];
evalc('x = lsqlin(C, p.d, p.A, p.b, p.Aeq, p.beq);');
obj = sum((p.C * x - p.d) .^ 2);
d = round(x);
""" + _SUMMARY.replace("out = {d, s};", "out = {d, s, obj};")


def lsqlin_solution(p):
    import scipy.sparse as sp

    C = sp.csr_matrix(p["C"]) if p["sp"] else p["C"]
    res = lsq.solve(C, p["d"], p["A"], p["b"], p["Aeq"], p["beq"], half_tol=p["ht"] or None)
    d = mround(res.x)[:, None] + 0.0
    fr = res.free[:, None]
    nh = res.near_half[:, None] & ~fr
    obj = float(np.sum((p["C"] @ res.x - p["d"][:, 0]) ** 2))
    return d, nh, fr, obj


def lsqlin_flags(p):
    d, nh, fr, _ = lsqlin_solution(p)
    m = nh | fr
    return {"near": nh, "free": fr, "pyv": d.ravel(order="F")[m.ravel(order="F")][:, None]}


def py_lsqlin(p, f):
    d, _, _, obj = lsqlin_solution(p)
    m = f["near"] | f["free"]
    d[m] = 0.5
    s = np.array([[np.count_nonzero(f["near"]), 0, 0, np.count_nonzero(f["free"]), 0, 0]], float)
    return [d, s, obj]


def objective_check(ref, act):
    """Python's objective is MATLAB's or lower, to 1e-9 (MATLAB stops at its tolerances)"""
    r = float(np.asarray(ref.data).ravel()[0])
    a = float(np.asarray(act.data).ravel()[0])
    if a > r + 1e-9 * (1.0 + abs(r)):
        return [f"objective {a!r} above MATLAB's {r!r} by more than 1e-9 relative"]
    return []


def _lsqlin_case(cid, make):
    return Case(
        cid,
        code=_LSQLIN_CODE,
        args=[In("p"), In("f")],
        nargout=3,
        python=py_lsqlin,
        inputs={"p": Value(make), "f": Value(lambda: lsqlin_flags(make()))},
        compare={1: SUMMARY, 2: Custom(objective_check, "objective at most MATLAB's (1e-9)")},
    )


def _lsqlin_cases():
    cs = []
    s = 9000
    for n, m, p, eq in [
        (3, 5, 2, 0),
        (5, 8, 4, 0),
        (8, 12, 10, 0),
        (10, 20, 15, 2),
        (15, 30, 20, 3),
        (20, 40, 30, 0),
        (25, 50, 40, 5),
        (6, 6, 12, 1),
        (12, 15, 30, 0),
        (30, 60, 60, 4),
    ]:
        for sparse_c in (False, True):
            s += 1
            cs.append(
                _lsqlin_case(
                    f"lsqlin_generic_{n}_{m}_{p}_{eq}_{'sp' if sparse_c else 'de'}",
                    lambda s=s, n=n, m=m, p=p, eq=eq, sc=sparse_c: generic_problem(
                        s, n, m, p, sparse_c=sc, eq=eq
                    ),
                )
            )
    for ng in (2, 3, 4, 5, 6, 8, 10):
        for tight in (False, True):
            s += 1
            cs.append(
                _lsqlin_case(
                    f"lsqlin_group_{ng}_{'tight' if tight else 'loose'}",
                    lambda s=s, ng=ng, t=tight: group_problem(s, ng, tight=t),
                )
            )
    # MATLAB's errors for inputs of the wrong size, and for an empty C
    e3, z31, z21 = np.eye(3), np.zeros((3, 1)), np.zeros((2, 1))
    for cid, args, ident in (
        ("err_d_size", [e3, z21, e3, z31], "optimlib:lsqlin:InvalidCAndD"),
        ("err_A_cols", [e3, z31, np.eye(2), z21], "optimlib:lsqlin:CAndA"),
        ("err_b_size", [e3, z31, e3, z21], "optimlib:lsqlin:InvalidAAndB"),
        ("err_Aeq_cols", [e3, z31, e3, np.ones((3, 1)), np.eye(2), z21], "optimlib:lsqlin:CAndAeq"),
        ("err_beq_size", [e3, z31, e3, np.ones((3, 1)), e3, z21], True),
        (
            "err_empty_c",
            [np.zeros((0, 3)), np.zeros((0, 1)), e3, z31],
            "optimlib:lsqlin:FirstTwoArgsEmpty",
        ),
    ):
        cs.append(
            Case(
                f"lsqlin_{cid}",
                matlab="lsqlin",
                args=args,
                python=lsq.lsqlin,
                expect_error=ident,
                py_error=lsq.LsqlinError,
            )
        )
    return cs


CASES = (
    [_assignment_case(cid, "grid", spec) for cid, spec in _grid_specs()]
    + [_assignment_case(cid, "global", spec) for cid, spec in _global_specs()]
    + _lsqlin_cases()
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- unit tests --------------------------------------------------------------------------------


def _kkt(C, d, A, b, res, Aeq=None, beq=None):
    """the largest violation of the optimality conditions, relative: stationarity with
    nonnegative multipliers on the active rows, feasibility"""
    from scipy.optimize import nnls

    C = np.asarray(C.toarray() if hasattr(C, "toarray") else C)
    A = np.asarray(A.toarray() if hasattr(A, "toarray") else A)
    x = res.x
    g = C.T @ (C @ x - d)
    act = A[res.active]
    E = np.zeros((0, x.size)) if Aeq is None else np.asarray(Aeq)
    # g + act' lam + E' mu = 0, lam >= 0 (mu free: as +/- pairs)
    M = np.hstack([act.T, E.T, -E.T])
    if M.shape[1]:
        _, rn = nnls(M, -g, maxiter=50 * M.shape[1])
    else:
        rn = float(np.linalg.norm(g))
    scale = 1.0 + np.abs(C.T @ d).max()
    feas = max(0.0, float(np.max(A @ x - b, initial=0.0)))
    if Aeq is not None and len(beq):
        feas = max(feas, float(np.abs(E @ x - beq).max()))
    return max(rn / scale, feas / (1.0 + np.abs(b).max()))


FAMILIES = [
    ("conn", {"filt": 0.0}),
    ("filt", {"filt": 0.2}),
    ("tight", {"max_off": (5, 5, 2), "off": (12, 12, 4)}),
    ("free", {"free_tiles": 3}),
    ("tightfree", {"filt": 0.3, "max_off": (8, 8, 3), "free_tiles": 2}),
    ("activeall", {"max_off": (2, 2, 1), "off": (10, 10, 3)}),
    ("stitch2d", {"stitch2d": True, "filt": 0.2}),
    ("two", {"two_groups": True}),
]


@pytest.mark.parametrize("name,kw", FAMILIES, ids=[f[0] for f in FAMILIES])
@pytest.mark.parametrize("size", [(3, 3), (6, 5), (12, 12)])
def test_optimum_to_1e9(name, kw, size):
    p = grid_problem(77 + size[0], size[0], size[1], **kw)
    systems, form = grid_systems(p)
    for C, d, A, b in systems:
        res = lsq.solve(C, d, A, b, form=form)
        assert _kkt(C, d, A, b, res) < 1e-9


@pytest.mark.parametrize("ng", [3, 6])
@pytest.mark.parametrize("tight", [False, True])
def test_optimum_with_equalities(ng, tight):
    p = group_problem(5 + ng, ng, tight=tight)
    res = lsq.solve(p["C"], p["d"], p["A"], p["b"], p["Aeq"], p["beq"])
    assert _kkt(p["C"], p["d"][:, 0], p["A"], p["b"][:, 0], res, p["Aeq"], p["beq"][:, 0]) < 1e-9


@pytest.mark.parametrize("seed", range(6))
def test_optimum_generic_rows(seed):
    p = generic_problem(300 + seed, 8 + seed, 16 + 2 * seed, 12 + seed, eq=seed % 3)
    res = lsq.solve(p["C"], p["d"], p["A"], p["b"], p["Aeq"], p["beq"])
    assert not res.info["structured"]
    assert _kkt(p["C"], p["d"][:, 0], p["A"], p["b"][:, 0], res, p["Aeq"], p["beq"][:, 0]) < 1e-9


def test_generic_rows_without_lineality_have_no_free_values():
    """a full-rank C has no lineality direction: the sparse form flags nothing as free"""
    import scipy.sparse as sp

    p = generic_problem(321, 10, 20, 15, eq=1)
    for C in (p["C"], sp.csr_matrix(p["C"])):
        res = lsq.solve(C, p["d"], p["A"], p["b"], p["Aeq"], p["beq"])
        assert not res.free.any()
        assert np.unique(res.groups).size == 1


def test_mean_zero_in_each_group():
    p = grid_problem(11, 5, 4, two_groups=True, isolated=2, filt=0.1)
    C, d, A, b = grid_systems(p)[0][0]
    res = lsq.solve(C, d, A, b)
    for g in np.unique(res.groups):
        assert abs(math.fsum(res.x[res.groups == g])) < 1e-12
    assert np.unique(res.groups).size == 4  # two grids and two isolated tiles


def test_free_offsets_at_the_analytic_center():
    """a tile whose pairs are all below the threshold is set by its bounds alone: its value
    maximizes the sum of the logarithms of the slacks of its rows"""
    p = grid_problem(12, 4, 4, filt=0.0, free_tiles=1)
    (C, d, A, b) = grid_systems(p)[0][0]
    res = lsq.solve(C, d, A, b)
    free = np.nonzero(res.free)[0]
    assert free.size >= 1
    t = free[0]
    rows = A[:, t] != 0
    s = b[rows] - A[rows] @ res.x
    assert np.all(s > 0)
    # d/dx_t sum log s_i = -sum a_it / s_i = 0
    assert abs(np.sum(A[rows, t] / s)) < 1e-10 * np.sum(1 / s)
    # the others are not free; the window of their rounding is widened
    assert np.all(res.tol[~res.free] > lsq.HALF_TOL["dense"])


def test_flags_without_free_values():
    p = grid_problem(13, 6, 6, filt=0.0)
    (C, d, A, b) = grid_systems(p)[0][0]
    res = lsq.solve(C, d, A, b)
    assert not res.free.any()
    w = max(lsq.HALF_TOL["dense"], lsq.HALF_REL["dense"] * np.abs(res.x).max())
    assert np.all(res.tol == w)
    assert np.array_equal(res.near_half, lsq.near_half_integer(res.x, w))
    assert res.near_half.mean() < 0.15


def test_sparse_form_with_two_groups_is_free():
    """MATLAB's sparse solver leaves each group at an arbitrary offset; its caller subtracts one
    mean, so with two groups no value is determined"""
    import scipy.sparse as sp

    p = grid_problem(14, 4, 4, two_groups=True)
    (C, d, A, b) = grid_systems(p)[0][0]
    res = lsq.solve(sp.csr_matrix(C), d, sp.csr_matrix(A), b)
    assert res.form == "sparse" and res.free.all()
    dense = lsq.solve(C, d, A, b)
    assert not dense.free.any()
    assert np.array_equal(dense.x, res.x)


def test_near_half_integer():
    x = np.array([0.5, 1.495, 1.505, -0.5, -2.51, 3.0, 2.49, 2.4799])
    assert lsq.near_half_integer(x, 0.02).tolist() == [
        True,
        True,
        True,
        True,
        True,
        False,
        True,
        False,
    ]


@pytest.mark.parametrize("seed", range(4))
def test_dense_and_sparse_inputs_give_the_same_solution(seed):
    import scipy.sparse as sp

    p = grid_problem(20 + seed, 6, 7, filt=0.3, free_tiles=1, max_off=(9, 9, 3))
    for C, d, A, b in grid_systems(p)[0]:
        a = lsq.solve(C, d, A, b)
        s = lsq.solve(sp.csr_matrix(C), d, sp.csr_matrix(A), b)
        assert np.array_equal(a.x, s.x)


def test_structured_and_generic_paths_agree(monkeypatch):
    """the graph path and the dense null-space path give the same solution"""
    p = grid_problem(31, 4, 3, filt=0.2, free_tiles=1, max_off=(7, 7, 3))
    systems, _ = grid_systems(p)
    ref = [lsq.solve(C, d, A, b).x for C, d, A, b in systems]
    monkeypatch.setattr(lsq, "_difference_rows", lambda M: None)
    for (C, d, A, b), xr in zip(systems, ref, strict=False):
        res = lsq.solve(C, d, A, b)
        assert not res.info["structured"]
        assert np.abs(res.x - xr).max() < 1e-8


def test_large_grid_sparse_matches_dense():
    """StitchShiftAssignmentTest.largeGridMatchesDenseSolver: above 400 tiles the sparse form
    (moved to mean zero) gives the dense form's rounded shifts, and both find the simulated
    positions (noise 0.7; with x pairs weighted 0.1 a position's error has a deviation of about
    1.3 voxels at the far side of the grid)"""
    p, true = grid_problem(1, 21, 21, filt=0.05, return_true=True)
    d_sparse, _, _ = assign(p, "grid")
    systems, form = grid_systems(p)
    assert form == "sparse"
    xd = np.column_stack(
        [
            lsq.solve(np.asarray(C.toarray()), dd, np.asarray(A.toarray()), b, form="dense").x
            for C, dd, A, b in systems
        ]
    )
    d_dense = mround(xd)
    d_dense = d_dense - d_dense[0]
    assert np.array_equal(d_sparse, d_dense)
    err = np.abs(d_sparse - (true - true[0]))
    assert np.mean(err <= 3) >= 0.95 and err.max() <= 7


def test_results_do_not_depend_on_the_thread_count(tmp_path):
    code = (
        f"import sys, numpy as np; sys.path.insert(0, {os.path.dirname(__file__)!r});"
        "import test_lsq as t; from petakit5d.matlabcompat import lsq;"
        "p = t.grid_problem(41, 12, 10, filt=0.3, free_tiles=2, max_off=(8, 8, 3));"
        "out = [lsq.solve(*s, form=f).x for f in ('dense',) for s in t.grid_systems(p)[0]];"
        "q = t.grid_problem(42, 25, 20, filt=0.1);"
        "out += [lsq.solve(*s).x for s in t.grid_systems(q)[0]];"
        "np.save(sys.argv[1], np.concatenate(out))"
    )
    outs = []
    for n in ("1", "4"):
        env = dict(os.environ, OMP_NUM_THREADS=n, OPENBLAS_NUM_THREADS=n, MKL_NUM_THREADS=n)
        f = tmp_path / f"x{n}.npy"
        subprocess.run([sys.executable, "-c", code, str(f)], env=env, check=True)
        outs.append(np.load(f))
    assert outs[0].tobytes() == outs[1].tobytes()


def test_without_constraints_is_least_squares():
    r = Rand(5)
    C = r.randn(40).reshape(10, 4)
    d = r.randn(10)
    x = lsq.lsqlin(C, d)
    assert np.abs(C.T @ (C @ x - d)).max() < 1e-12
    with pytest.raises(lsq.LsqlinError):
        lsq.lsqlin(np.zeros((3, 2)), np.zeros(3))


def test_infeasible_raises():
    C = np.array([[-1.0, 1.0]])
    A = np.array([[-1.0, 1.0], [1.0, -1.0]])
    b = np.array([-2.0, 1.0])  # x2 - x1 <= -2 and x2 - x1 >= -1
    with pytest.raises(lsq.LsqlinError):
        lsq.lsqlin(C, np.array([0.0]), A, b)


def test_sizes_are_checked():
    with pytest.raises(lsq.LsqlinError) as e:
        lsq.lsqlin(np.eye(3), np.zeros(2), np.eye(3), np.zeros(3))
    assert e.value.identifier
    with pytest.raises(lsq.LsqlinError):
        lsq.lsqlin(np.eye(3), np.zeros(3), np.eye(2), np.zeros(2))


def test_assignment_of_1600_tiles_is_fast():
    import time

    p = grid_problem(3, 40, 40, filt=0.05)
    t0 = time.perf_counter()
    assign(p, "grid")
    assert time.perf_counter() - t0 < 5.0
