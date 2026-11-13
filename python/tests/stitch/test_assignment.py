"""Tests of petakit5d.stitch.assignment: the stitching shift assignment.

The oracle cases run MATLAB's own functions on generated problems (2 to 1,600 tiles) and compare
the shifts exactly: stitch_global_grid_assignment (dense up to 400 tiles, sparse above),
stitch_global_assignment, stitch_local_assignment, stitch_global_group_assignment (the groups'
grid problems, pairs of groups, free groups), stitch_global_grid_assignment_wrapper on its .mat
files, and stitch_shift_assignment itself on tiles whose registration results are on disk (so no
registration runs): the pairs of each method, the signs, the bounds, the large-file names, the
group file and the 'test' method. Inputs are in MATLAB's form (1-based tile numbers, [y x z]
axis weights); the Python side converts them by hand.

MATLAB tests whose intent is ported (tests/stitch/StitchShiftAssignmentTest.m, all 12 methods:
the ``test_`` functions below with their names, and their fixed problems as the oracle cases
``mtest_*``), and the assignment part of StitchingTest's xcorrRecoversPerturbedTilePosition,
largeFileRegistrationRecoversPerturbedTilePosition and localShiftMethod (from the registration
results a correct registration gives; the registration is S04's).
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.axes import axis_string_from_matlab, from_yxz
from petakit5d.core.dtypes import mround
from petakit5d.core.text import through_mat2str
from petakit5d.io.matfile import read_mat, write_mat
from petakit5d.matlabcompat import lsq
from petakit5d.stitch import assignment as sa
from petakit5d.stitch import geometry
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value

# --- generated problems (MATLAB's form: tile numbers from 1, axis weights [y x z]) -------------


def _rng(seed, name="p"):
    return synth.CounterRNG(seed, name)


def _neighbors(idx, corners=False):
    """pairs (i < j, by i, then j; 0-based) of tiles one apart (``corners``: also diagonal in
    x and y)"""
    pairs = []
    for i in range(len(idx)):
        dd = np.abs(idx[i + 1 :] - idx[i])
        ok = ((dd[:, 0] + dd[:, 1]) + dd[:, 2]) == 1
        if corners:
            ok |= (dd.max(axis=1, initial=0) == 1) & (dd[:, 2] == 0)
        pairs.extend((i, i + 1 + int(j)) for j in np.nonzero(ok)[0])
    return np.array(pairs, dtype=np.int64).reshape(-1, 2)


def _bounds(lo, hi, order):
    o = order.astype(np.float64)
    return np.hstack([lo * o + hi * (o - 1), hi * o + lo * (o - 1)])


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
    free_bound=None,
    jitter=False,
    return_true=False,
    free_keep_corners=False,
):
    """a caller's inputs for a grid of tiles (x, y, z tile indices; y fastest): pairs of
    neighbors (``corners``: and of diagonal ones), noisy relative shifts of true offsets,
    correlations (a share ``filt`` below the threshold; ``free_tiles`` tiles with all their pairs
    below, or with ``free_keep_corners`` only their pairs with neighbors, so that the grid
    method's weights, 0 for diagonal pairs, leave them free), the bounds of the grid method from
    the overlap size and max offset (x, y, z), turned around where the second tile is not ahead
    (``jitter``: stage positions jittered so that some are not); ``free_bound`` (lo, hi): free
    tiles bounded to lo .. hi from each neighbor"""
    r = _rng(seed, "grid")
    idx = np.array(
        [(x, y, z) for z in range(gz) for x in range(gx) for y in range(gy)], dtype=np.float64
    ).reshape(-1, 3)
    if two_groups:
        other = idx.copy()
        other[:, 0] += gx + 2
        idx = np.vstack([idx, other])
    if holes:
        keep = np.ones(len(idx), bool)
        keep[r.integers(1, len(idx) - 1, holes)] = False
        idx = idx[keep]
    if isolated:
        idx = np.vstack([idx, [[1000.0 + 3 * k, 1000.0, 0.0] for k in range(isolated)]])
    n = len(idx)
    if perm:
        idx = idx[np.argsort(r.uniform(n), kind="stable")]
    tile_idx = np.hstack([idx, np.ones((n, 1))])
    xyz = idx * np.array([900.0, 900.0, 300.0])
    if jitter:
        xyz = xyz + np.floor(r.uniform((n, 3)) * 3) * 400.0
    true = np.column_stack([r.integers(-o, o, n) for o in off]).astype(np.float64)
    if stitch2d:
        true[:, 2] = 0
    pr = _neighbors(idx, corners)
    ti, tj = pr[:, 0], pr[:, 1]
    p = len(ti)
    rel = true[tj] - true[ti] + r.normal((p, 3)) * noise
    if stitch2d:
        rel[:, 2] = 0.0
    xc = 0.3 + 0.7 * r.uniform(p)
    xc[r.uniform(p) < filt] = 0.1
    ft = r.integers(0, n - 1, free_tiles) if free_tiles else np.zeros(0, np.int64)
    one_apart = np.abs(idx[ti] - idx[tj]).sum(axis=1) == 1
    for t in ft:
        xc[((ti == t) | (tj == t)) & (one_apart | (not free_keep_corners))] = 0.1
    mo = np.array(max_off, dtype=np.float64)
    lo = -np.ones((p, 1)) * mo
    hi = np.minimum(np.array(overlap, dtype=np.float64) - 3, mo) * np.ones((p, 1))
    if stitch2d:
        hi[:, 2] = 0.0
    mas = _bounds(lo, hi, xyz[tj] >= xyz[ti])
    if free_bound is not None:
        a, b = free_bound
        for t in ft:
            mas[tj == t] = [a, a, a, b, b, b]
            mas[ti == t] = [-b, -b, -b, -a, -a, -a]
    one = np.column_stack([ti + 1, tj + 1]).astype(np.float64)
    p = {
        "nF": float(n),
        "mx": np.column_stack([one, xc]),
        "ash": np.column_stack([one, rel]),
        "pairs": one,
        "mas": mas,
        "thr": float(thr),
        "tileIdx": tile_idx,
        "aw": np.array(aw, dtype=np.float64),
    }
    return (p, true) if return_true else p


def group_problem(
    seed,
    gx,
    gy,
    *,
    bands=2,
    labels=None,
    filt=0.05,
    noise=0.7,
    free_groups=(),
    perm=False,
    corners=True,
    max_off=(20, 20, 6),
    thr=0.25,
    aw=(1.0, 0.1, 10.0),
    group_off=8,
    islands=0,
    vary_bounds=False,
):
    """stitch_global_group_assignment's inputs: a gx x gy grid in ``bands`` groups along x (group
    ids ``labels``, by band), each group offset by up to ``group_off``; pairs of neighbors (and of
    diagonal ones across groups with ``corners``), bounds within groups from the max offset,
    twice it between groups (as stitch_shift_assignment sets them; ``vary_bounds``: 0 to 3 voxels
    wider, per pair and side, so that the groups' bounds are the widest of their pairs');
    ``free_groups``: bands whose pairs to other groups are all below the threshold; ``perm``:
    tiles numbered at random (pairs whose first tile is in the higher group are not used between
    groups, as MATLAB has it); ``islands``: groups of one tile, overlapping the last band"""
    r = _rng(seed, "group")
    idx = np.array([(x, y, 0) for x in range(gx) for y in range(gy)], dtype=np.float64)
    band = np.minimum((idx[:, 0] * bands) // gx, bands - 1).astype(np.int64)
    if islands:
        isl = np.array([(gx, k, 0) for k in range(islands)], dtype=np.float64)
        idx = np.vstack([idx, isl])
        band = np.r_[band, bands + np.arange(islands)]
    n = len(idx)
    if perm:
        o = np.argsort(r.uniform(n), kind="stable")
        idx, band = idx[o], band[o]
    labels = np.arange(1, bands + islands + 1) if labels is None else np.asarray(labels)
    grp = labels[band].astype(np.float64)
    tile_idx = np.hstack([idx, np.ones((n, 1))])
    xyz = idx * np.array([900.0, 900.0, 300.0])
    goff = np.column_stack([r.integers(-group_off, group_off, bands + islands) for _ in range(3)])
    true = goff[band] + np.column_stack([r.integers(-4, 4, n) for _ in range(3)])
    pr = _neighbors(idx, False)
    if corners:
        pc = _neighbors(idx, True)
        cross = band[pc[:, 0]] != band[pc[:, 1]]
        pr = np.unique(np.vstack([pr, pc[cross]]), axis=0)
    ti, tj = pr[:, 0], pr[:, 1]
    p = len(ti)
    rel = true[tj] - true[ti] + r.normal((p, 3)) * noise
    xc = 0.3 + 0.7 * r.uniform(p)
    xc[r.uniform(p) < filt] = 0.1
    inter = band[ti] != band[tj]
    for b in free_groups:
        xc[inter & ((band[ti] == b) | (band[tj] == b))] = 0.1
    mo = np.array(max_off, dtype=np.float64)
    mas = _bounds(-np.ones((p, 1)) * mo, np.ones((p, 1)) * mo, xyz[tj] >= xyz[ti])
    mas[inter] = np.hstack([-np.ones((inter.sum(), 1)) * mo, np.ones((inter.sum(), 1)) * mo]) * 2
    if vary_bounds:
        wider = r.integers(0, 3, (int(inter.sum()), 6)).astype(np.float64)
        mas[inter] += np.hstack([-wider[:, 0:3], wider[:, 3:6]])
    one = np.column_stack([ti + 1, tj + 1]).astype(np.float64)
    return {
        "nF": float(n),
        "mx": np.column_stack([one, xc]),
        "ash": np.column_stack([one, rel]),
        "pairs": one,
        "mas": mas,
        "thr": float(thr),
        "tileIdx": tile_idx,
        "aw": np.array(aw, dtype=np.float64),
        "grp": grp.reshape(-1, 1),
    }


def py_problem(p):
    """the Python form of a caller's inputs"""
    n = int(p["nF"])
    mx = np.array(p["mx"], dtype=np.float64).reshape(-1, 3)
    ash = np.array(p["ash"], dtype=np.float64).reshape(-1, 5)
    mx[:, 0:2] -= 1
    ash[:, 0:2] -= 1
    om = np.zeros((n, n), bool)
    pr = np.asarray(p["pairs"], dtype=np.int64).reshape(-1, 2) - 1
    om[pr[:, 0], pr[:, 1]] = True
    return n, mx, ash, om, np.asarray(p["mas"]).reshape(-1, 6), p["thr"]


def py_grid(p):
    n, mx, ash, om, mas, thr = py_problem(p)
    return sa.grid_assignment(n, mx, ash, om, mas, thr, p["tileIdx"], from_yxz(p["aw"]))


def py_global(p):
    return sa.global_assignment(*py_problem(p))


def py_local(p):
    n, mx, ash, om, _, thr = py_problem(p)
    return sa.local_assignment(n, mx, ash, om, thr)


def py_group(p):
    n, mx, ash, om, mas, thr = py_problem(p)
    return sa.group_assignment(
        n, mx, ash, om, mas, thr, p["tileIdx"], from_yxz(p["aw"]), np.ravel(p["grp"])
    )


_OV = "ov = false(p.nF); ov(sub2ind([p.nF, p.nF], p.pairs(:, 1), p.pairs(:, 2))) = true;\n"
_CODE = {
    "grid": "d = stitch_global_grid_assignment(p.nF, p.mx, p.ash, ov, p.mas, p.thr, p.tileIdx, "
    "p.aw);",
    "global": "d = stitch_global_assignment(p.nF, p.mx, p.ash, ov, p.mas, p.thr);",
    "local": "d = stitch_local_assignment(p.nF, p.mx, p.ash, ov, p.thr);",
    "group": "d = stitch_global_group_assignment(p.nF, p.mx, p.ash, ov, p.mas, p.thr, p.tileIdx, "
    "p.aw, p.grp, [], '', false, 1, false, '');",
}
_PY = {"grid": py_grid, "global": py_global, "local": py_local, "group": py_group}


def _code(method):
    call = _CODE[method].replace("'", "''")
    return f"p = in{{1}};\n{_OV}d = [];\nevalc('{call}');\nout = {{d}};\n"


def caller_case(cid, method, make, **kw):
    return Case(
        cid,
        code=_code(method),
        args=[In("p")],
        python=_PY[method],
        inputs={"p": Value(make)},
        **kw,
    )


def _grid_specs():
    """(id, (seed, gx, gy, gz, options)): 2 to 1,600 tiles, every family"""
    out = []
    s = 300

    def add(name, gx, gy, gz=1, **kw):
        nonlocal s
        s += 1
        out.append((f"grid.{name}.{gx}x{gy}x{gz}.s{s}", (s, gx, gy, gz, kw)))

    small = [(2, 1), (1, 2), (2, 2), (3, 1), (3, 2), (2, 3), (3, 3), (4, 3), (4, 4), (5, 4)]
    small += [(5, 5), (6, 5), (7, 6), (8, 7), (10, 10)]
    mid = [(12, 11), (14, 14), (17, 16), (20, 20)]
    for gx, gy in small + mid:
        add("conn", gx, gy, filt=0.0)
        add("filt", gx, gy, filt=0.15)
        add("tight", gx, gy, max_off=(6, 6, 3), off=(12, 12, 4))
    for gx, gy in small[4:] + mid:
        add("free", gx, gy, free_tiles=max(1, gx * gy // 12))
        add("bigfilt", gx, gy, filt=0.45)
        add("tightfree", gx, gy, filt=0.25, max_off=(8, 8, 3), free_tiles=2)
        add("jitter", gx, gy, filt=0.1, jitter=True, max_off=(10, 10, 4))
        add("weights", gx, gy, aw=(2.0, 0.5, 3.0), filt=0.1)
    for gx, gy in small[2:]:
        add("stitch2d", gx, gy, stitch2d=True, filt=0.2)
        add("thr", gx, gy, thr=0.6, filt=0.0)
    for gx, gy, gz in [(2, 2, 2), (3, 2, 2), (3, 3, 2), (4, 3, 3), (5, 4, 3), (7, 7, 4)]:
        add("conn3d", gx, gy, gz, filt=0.0)
        add("filt3d", gx, gy, gz, filt=0.2)
        add("tight3d", gx, gy, gz, max_off=(6, 6, 2), off=(12, 12, 4))
        add("free3d", gx, gy, gz, free_tiles=2)
    for gx, gy in [(3, 3), (4, 4), (6, 5), (8, 8), (11, 10)]:
        add("holes", gx, gy, holes=3, perm=True, filt=0.1)
        add("two", gx, gy, two_groups=True, filt=0.1)
        add("isolated", gx, gy, isolated=2, filt=0.1)
        add("activeall", gx, gy, max_off=(2, 2, 1), off=(10, 10, 3))
        add("corners", gx, gy, corners=True, filt=0.1)
    # the sparse form, above 400 tiles
    for gx, gy, gz in [(21, 20, 1), (25, 20, 1), (30, 30, 1), (40, 40, 1), (9, 9, 6)]:
        add("conn", gx, gy, gz, filt=0.0)
        add("filt", gx, gy, gz, filt=0.1)
        add("free", gx, gy, gz, free_tiles=6)
        add("tight", gx, gy, gz, max_off=(6, 6, 3), off=(12, 12, 4))
        add("bigfilt", gx, gy, gz, filt=0.4)
    add("two", 21, 11, two_groups=True, filt=0.05)
    add("holes", 25, 20, holes=10, perm=True, filt=0.1)
    add("jitter", 30, 30, filt=0.1, jitter=True, max_off=(10, 10, 4))
    add("corners", 22, 20, corners=True, filt=0.1)
    # lines
    for k in (4, 7, 12, 20, 35):
        add("line", 1, k, filt=0.1)
        add("linetight", k, 1, max_off=(5, 5, 2), off=(12, 12, 4))
        add("linefree", 1, k, free_tiles=1)
    # free tiles whose stage positions relative to their neighbors are out of their bounds
    for gx, gy in [(2, 2), (3, 3), (4, 3), (5, 5), (8, 6), (12, 10)]:
        add("freebound", gx, gy, filt=0.0, free_tiles=1, free_bound=(2.0, 60.0))
    add("freebound", 22, 20, filt=0.0, free_tiles=4, free_bound=(2.0, 60.0))
    # every pair below the threshold: the stage positions
    add("allbelow", 4, 4, filt=1.0)
    add("allbelow", 21, 20, filt=1.0)
    # free tiles whose diagonal pairs are above the threshold (weight 0: no rows in C)
    for gx, gy in [(3, 3), (4, 4), (6, 5), (21, 20)]:
        add("cornersfree", gx, gy, corners=True, filt=0.0, free_tiles=1, free_keep_corners=True)
    return out


def _global_specs():
    out = []
    s = 6000
    for gx, gy in [(2, 2), (3, 2), (3, 3), (4, 4), (5, 5), (6, 5), (7, 7), (9, 9), (12, 10)]:
        for name, kw in (
            ("conn", {"filt": 0.0}),
            ("filt", {"filt": 0.2}),
            ("tight", {"max_off": (6, 6, 3), "off": (12, 12, 4)}),
            ("free", {"free_tiles": 2, "filt": 0.1}),
            ("jitter", {"filt": 0.1, "jitter": True, "max_off": (10, 10, 4)}),
        ):
            s += 1
            out.append((f"global.{name}.{gx}x{gy}.s{s}", (s, gx, gy, 1, dict(kw, corners=True))))
    for gx, gy, gz in [(3, 3, 2), (21, 20, 1)]:
        s += 1
        out.append((f"global.conn.{gx}x{gy}x{gz}.s{s}", (s, gx, gy, gz, {"filt": 0.05})))
    s += 1
    out.append((f"global.allbelow.4x4.s{s}", (s, 4, 4, 1, {"filt": 1.0})))
    return out


def _local_specs():
    out = []
    s = 7000
    for gx, gy, gz, kw in [
        (2, 1, 1, {}),
        (3, 3, 1, {"filt": 0.2}),
        (4, 4, 1, {"filt": 0.0, "corners": True}),
        (6, 5, 1, {"filt": 0.3}),
        (5, 5, 2, {"filt": 0.1}),
        (8, 8, 1, {"holes": 4, "perm": True, "filt": 0.1}),
        (6, 6, 1, {"two_groups": True}),
        (5, 4, 1, {"isolated": 2}),
        (10, 9, 1, {"corners": True, "filt": 0.25}),
        (21, 20, 1, {"filt": 0.1}),
        (40, 40, 1, {"filt": 0.1}),
        (4, 4, 1, {"filt": 1.0}),
    ]:
        s += 1
        out.append((f"local.{gx}x{gy}x{gz}.s{s}", (s, gx, gy, gz, kw)))
    return out


def _group_specs():
    out = []
    s = 8000
    for name, gx, gy, kw in [
        ("two", 4, 3, {}),
        ("three", 6, 3, {"bands": 3}),
        ("four", 8, 4, {"bands": 4, "filt": 0.1}),
        ("labels", 6, 4, {"bands": 3, "labels": (7, 3, 12)}),
        ("perm", 6, 4, {"bands": 3, "perm": True}),
        ("permlabels", 8, 5, {"bands": 4, "perm": True, "labels": (5, 2, 9, 1)}),
        ("free", 6, 4, {"bands": 3, "free_groups": (1,)}),
        ("freeend", 6, 4, {"bands": 3, "free_groups": (2,)}),
        ("freetwo", 8, 4, {"bands": 4, "free_groups": (1, 2)}),
        ("allfree", 6, 3, {"bands": 3, "free_groups": (0, 1, 2)}),
        ("nocorners", 6, 4, {"bands": 3, "corners": False}),
        ("one", 4, 4, {"bands": 1}),
        ("islands", 6, 3, {"bands": 2, "islands": 2}),
        ("islandsfree", 6, 3, {"bands": 3, "islands": 1, "free_groups": (1,)}),
        ("weights", 6, 4, {"bands": 3, "aw": (2.0, 0.5, 3.0)}),
        ("filt", 9, 6, {"bands": 3, "filt": 0.3}),
        ("noisy", 10, 6, {"bands": 5, "noise": 2.0}),
        ("big", 16, 12, {"bands": 4, "filt": 0.1}),
        ("bigfree", 16, 12, {"bands": 4, "free_groups": (1,), "perm": True}),
        ("sparse", 24, 18, {"bands": 3, "filt": 0.05}),
        # bounds between groups that hold the groups (group offsets up to 8, bounds 4)
        ("tightinter", 6, 4, {"bands": 3, "max_off": (2, 2, 1)}),
        ("tightinterperm", 8, 5, {"bands": 4, "max_off": (3, 3, 1), "perm": True}),
        ("varybounds", 6, 4, {"bands": 3, "max_off": (2, 2, 1), "vary_bounds": True}),
        (
            "varyboundsperm",
            9,
            5,
            {"bands": 3, "max_off": (2, 2, 1), "vary_bounds": True, "perm": True},
        ),
    ]:
        s += 1
        out.append((f"group.{name}.{gx}x{gy}.s{s}", (s, gx, gy, kw)))
    return out


# --- StitchShiftAssignmentTest's fixed problems ------------------------------------------------


def mtest_problem(tile_idx, pairs, xc, rel, mas, aw=(1.0, 1.0, 1.0), thr=0.25, grp=None):
    """a caller's inputs from the MATLAB test's arguments (pairs 1-based)"""
    pr = np.asarray(pairs, dtype=np.float64).reshape(-1, 2)
    p = len(pr)
    t = np.asarray(tile_idx, dtype=np.float64)
    out = {
        "nF": float(len(t)),
        "mx": np.column_stack([pr, np.asarray(xc, dtype=np.float64).reshape(p)]),
        "ash": np.column_stack([pr, np.asarray(rel, dtype=np.float64).reshape(p, 3)]),
        "pairs": pr,
        "mas": np.asarray(mas, dtype=np.float64).reshape(-1, 6) * np.ones((p, 1)),
        "thr": float(thr),
        "tileIdx": t,
        "aw": np.array(aw, dtype=np.float64),
    }
    if grp is not None:
        out["grp"] = np.asarray(grp, dtype=np.float64).reshape(-1, 1)
    return out


def _row(n):
    return np.column_stack([np.arange(n), np.zeros((n, 2)), np.ones(n)])


def _neighbor_pairs(tile_idx):
    return [tuple(int(v) + 1 for v in q) for q in _neighbors(np.asarray(tile_idx)[:, :3])]


def _mtests():
    """id -> (problem, methods, MATLAB test's expected shifts or None, columns)"""
    out = {}
    t2 = np.array([[0, 0, 0, 1], [1, 0, 0, 1]], dtype=np.float64)
    out["one_pair"] = (
        mtest_problem(t2, [(1, 2)], [0.9], [[8, -4, 2]], [-50, -50, -10, 50, 50, 10], (1, 0.1, 10)),
        ("grid", "global"),
        np.array([[0, 0, 0], [8, -4, 2]]),
    )
    t3 = _row(3)
    mas = np.array([[-50, -50, -50, 50, 50, 50], [2, -50, -50, 8, 50, 50]])
    out["free_tile_held_by_its_bounds"] = (
        mtest_problem(t3, _neighbor_pairs(t3), [0.9, 0.1], np.zeros((2, 3)), mas),
        ("grid",),
        np.array([[0, 0, 0], [0, 0, 0], [2, 0, 0]]),
    )
    rel = np.array([[10, 0, 0], [0, 0, 0], [0, 0, 0]])
    mas = np.array(
        [[-50, -50, -50, 4, 50, 50], [-1, -50, -50, 1, 50, 50], [-1, -50, -50, 1, 50, 50]]
    )
    out["bounds_of_a_free_tile_hold_its_neighbors"] = (
        mtest_problem(t3, [(1, 2), (1, 3), (2, 3)], [0.9, 0.1, 0.1], rel, mas),
        ("grid",),
        None,
    )
    t4 = np.array([[0, 0, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1], [1, 1, 0, 1]], dtype=np.float64)
    pos = np.array([[0, 0, 0], [0, 0, 0], [-2, 0, 0], [0, 0, 0]], dtype=np.float64)
    pr = _neighbor_pairs(t4)
    rel = np.array([pos[j - 1] - pos[i - 1] for i, j in pr])
    out["integer_shifts_round_alike"] = (
        mtest_problem(t4, pr, [0.9] * len(pr), rel, [-20, -20, -10, 20, 20, 10], (1, 0.1, 10)),
        ("grid", "global"),
        pos,
    )
    u = 4.5 + 3e-8
    mas = np.array([[-50, -50, -50, u, 50, 50], [-50, -50, -50, 50, 50, 50]])
    out["active_bounds_give_the_exact_optimum"] = (
        mtest_problem(t3, _neighbor_pairs(t3), [0.9, 0.9], [[10, 0, 0], [0, 0, 0]], mas),
        ("grid",),
        None,
    )
    # a correlation equal to the threshold is kept, one just below it is not
    out["threshold_equal"] = (
        mtest_problem(
            t3, _neighbor_pairs(t3), [0.25, 0.2499], [[3, 1, 0], [5, -2, 1]], [-50] * 3 + [50] * 3
        ),
        ("grid", "global", "local"),
        None,
    )
    t6 = _row(6)
    pr6 = _neighbor_pairs(t6)
    sx = np.array([2, 2, -2, 7, 4], dtype=np.float64)
    rel6 = np.column_stack([sx, np.zeros((5, 2))])
    grp6 = [1, 1, 2, 2, 3, 3]
    for name, xc in (
        ("groups_linked", [0.9, 0.9, 0.9, 0.1, 0.9]),
        ("groups_unlinked", [0.9, 0.1, 0.9, 0.1, 0.9]),
        ("groups_all_below", [0.1] * 5),
    ):
        out[name] = (
            mtest_problem(t6, pr6, xc, rel6, [-50, -50, -50, 50, 50, 50], grp=grp6),
            ("group",),
            None,
        )
    return out


_MTESTS = _mtests()


# --- the wrapper: a group's grid problem in .mat files -----------------------------------------

_MAT_HEADER = b"MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Thu Jan  1 00:00:00 2026"


def write_fixed_mat(path, variables):
    """a MAT v6 file whose header has no time or platform (inputs the same everywhere)"""
    write_mat(path, variables, version="6")
    data = bytearray(Path(path).read_bytes())
    data[:116] = _MAT_HEADER.ljust(116, b" ")
    Path(path).write_bytes(bytes(data))


def group_info(p):
    return {
        "nF_g": p["nF"],
        "max_xcorr_mat_g": p["mx"],
        "absolute_shift_mat_g": p["ash"],
        "overlap_matrix_g": _overlap_of(p),
        "max_allow_shift_g": p["mas"],
        "tileIdx_g": p["tileIdx"],
    }


def _overlap_of(p):
    n = int(p["nF"])
    om = np.zeros((n, n), bool)
    pr = np.asarray(p["pairs"], dtype=np.int64).reshape(-1, 2) - 1
    om[pr[:, 0], pr[:, 1]] = True
    return om


def wrapper_tree(make, existing=False):
    def fill(folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        write_fixed_mat(folder / "info.mat", group_info(make()))
        if existing:
            write_fixed_mat(folder / "res.mat", {"d_shift_g": np.zeros((1, 3))})

    return fill


_WRAPPER_CODE = """
p = in{1}; d = in{2};
stitch_global_grid_assignment_wrapper([d 'info.mat'], [d 'res.mat'], p.thr, p.aw, 'u1');
out = {};
"""


def py_wrapper(p, d):
    sa.global_grid_assignment(f"{d}info.mat", f"{d}res.mat", p["thr"], from_yxz(p["aw"]), "u1")
    return []


def wrapper_case(cid, make, existing=False):
    return Case(
        cid,
        code=_WRAPPER_CODE,
        args=[In("p"), InPath("d")],
        nargout=0,
        python=py_wrapper,
        inputs={"p": Value(make), "d": Tree(wrapper_tree(make, existing), "d")},
        files=True,
    )


# --- stitch_shift_assignment on tiles with their registration results on disk ------------------


def ssa_problem(
    seed,
    gx,
    gy,
    gz=1,
    *,
    method="grid",
    data_order="y,x,z",
    size_yxz=(64, 72, 20),
    overlap=(10, 12, 4),
    vs=(0.108, 0.108, 0.25),
    jitter=0.0,
    filt=0.1,
    noise=0.6,
    off=(6, 6, 3),
    max_off=(12, 12, 5),
    thr=0.25,
    aw=(1.0, 0.1, 10.0),
    stitch2d=False,
    groups=None,
    group_file=0,
    large_file=False,
    locs=None,
    even=False,
    zshift_2d=False,
):
    """stitch_shift_assignment's inputs for a grid of tiles (``size_yxz`` voxels, MATLAB's
    [ny nx nz] in ``data_order``'s layout) placed with ``overlap`` voxels (x, y, z) between
    neighbors and stage positions jittered by up to ``jitter`` um (so that tiles change order):
    the overlaps and regions from the geometry, and the registration results of the pairs the
    method registers (relative shifts of true offsets, signed by the tiles' order; in 2D without
    z unless ``zshift_2d``). ``groups``: a group per tile (column x band); ``group_file``: 4 or 5
    columns written to groups.csv (else the groups are the locations ``locs``)."""
    r = _rng(seed, "ssa")
    idx = np.array(
        [(x, y, z) for z in range(gz) for x in range(gx) for y in range(gy)], dtype=np.float64
    ).reshape(-1, 3)
    n = len(idx)
    vs = np.array(vs, dtype=np.float64)
    # sizes along x, y, z of the standard order; the image size in the data order's layout
    yxz = np.array(size_yxz, dtype=np.float64)
    ext_xyz = np.array([yxz[1], yxz[0], yxz[2]])
    dom = np.argsort(list(data_order.replace(",", "")), kind="stable")  # data_order_mat - 1
    im_size = np.zeros(3)
    for k in range(3):
        im_size[dom[k]] = ext_xyz[k]
    step = (ext_xyz - 1 - np.array(overlap, dtype=np.float64)) * vs
    xyz = idx * step
    if jitter:
        xyz = xyz + np.floor(r.uniform((n, 3)) * 64) / 64 * jitter
    xyz = xyz - xyz.min(axis=0)
    cub = np.hstack([xyz, xyz + (ext_xyz - 1) * vs])
    om, regions = geometry.find_overlaps(cub, stitch2d)
    tile_idx = np.hstack([idx, np.ones((n, 1)) if locs is None else np.asarray(locs)[:, None]])
    grp = None
    if groups:
        grp = np.minimum((idx[:, 0] * groups) // gx, groups - 1) + 1
    # the pairs the method registers
    reg = om.copy()
    ti, tj = np.nonzero(om.T)[::-1]
    if method in ("grid", "test"):
        corner = np.abs(tile_idx[ti, :3] - tile_idx[tj, :3]).sum(axis=1) != 1
        reg[ti[corner], tj[corner]] = False
    elif method == "group":
        g = grp if grp is not None else tile_idx[:, 3]
        same = g[ti] == g[tj]
        corner = np.abs(tile_idx[ti] - tile_idx[tj]).sum(axis=1) != 1
        reg[ti[same & corner], tj[same & corner]] = False
    true = np.column_stack([r.integers(-o, o, n) for o in off]).astype(np.float64)
    if even:
        true = 2 * true
    if stitch2d and not zshift_2d:
        true[:, 2] = 0
    ri, rj = np.nonzero(reg)
    pr = len(ri)
    absolute = true[rj] - true[ri] + mround(r.normal((pr, 3)) * noise * 4) / 4
    if stitch2d and not zshift_2d:
        absolute[:, 2] = 0.0
    flag = (0.5 - (xyz[ri] > xyz[rj])) * 2
    xc = 0.3 + 0.7 * r.uniform(pr)
    xc[r.uniform(pr) < filt] = 0.1
    return {
        "nF": float(n),
        "imSizes": np.tile(im_size, (n, 1)),
        "xyz": xyz,
        "vs": vs,
        "dataOrder": data_order,
        "ov": om,
        "regions": regions,
        "mo": np.array(max_off, dtype=np.float64),
        "ds": np.array([1.0, 1.0, 1.0]),
        "thr": float(thr),
        "tileIdx": tile_idx,
        "method": method,
        "s2d": bool(stitch2d),
        "aw": np.array(aw, dtype=np.float64),
        "gf": float(group_file),
        "large": bool(large_file),
        "pool": np.array([2.0, 2.0, 2.0]) if large_file else np.zeros((0, 0)),
        "regpairs": np.column_stack([ri + 1, rj + 1]).astype(np.float64).reshape(-1, 2),
        "rel": (absolute * flag).reshape(-1, 3),
        "xc": xc,
        "grp": np.zeros(0) if grp is None else grp,
    }


def ssa_tree(make):
    """the tiles' folder (tile 1 as a small tiff, for its size), the registration results (MAT
    v6 with a fixed header), the MIP slab folders of large files and the group file"""

    def fill(folder):
        p = make()
        folder = Path(folder)
        (folder / "tiles").mkdir(parents=True)
        (folder / "xcorr").mkdir()
        n = int(p["nF"])
        sz = [int(v) for v in p["imSizes"][0]]
        synth.write_tiff(folder / "tiles" / "tile_1.tif", np.ones((sz[2], sz[0], sz[1]), np.uint16))
        pr = np.asarray(p["regpairs"], dtype=np.int64).reshape(-1, 2)
        if p["large"]:
            mip = folder / "tiles" / "MIP_Slabs_pooling_2_2_2"
            mip.mkdir()
            for k in range(n):
                (mip / f"tile_{k + 1}_z.zarr").mkdir()
            for k in range(len(pr)):
                write_fixed_mat(
                    folder / "xcorr" / f"xcorr_tile_{pr[k, 0]}_{pr[k, 1]}.mat",
                    {
                        "relative_shift_mat": p["rel"][k : k + 1],
                        "max_xcorr_mat": p["xc"][k : k + 1].reshape(1, 1),
                        "pair_indices": pr[k : k + 1].astype(np.float64),
                    },
                )
        else:
            for f in np.unique(pr[:, 0]):
                k = np.nonzero(pr[:, 0] == f)[0]
                # in the task's order: by second tile
                k = k[np.argsort(pr[k, 1], kind="stable")]
                write_fixed_mat(
                    folder / "xcorr" / f"xcorr_tile_{f}.mat",
                    {
                        "relative_shift_mat": p["rel"][k],
                        "max_xcorr_mat": p["xc"][k].reshape(-1, 1),
                        "pair_indices": pr[k].astype(np.float64),
                    },
                )
        if p["gf"]:
            g = np.asarray(p["grp"])
            t = p["tileIdx"]
            rows = ["x,y,z,g"] if p["gf"] == 4 else ["x,y,z,did,g"]
            for k in range(n):
                v = [t[k, 0], t[k, 1], t[k, 2]] + ([] if p["gf"] == 4 else [t[k, 3]]) + [g[k]]
                rows.append(",".join(f"{x:g}" for x in v))
            (folder / "groups.csv").write_text("\n".join(rows) + "\n")

    return fill


_SSA_CODE = """
p = in{1}; d = in{2};
paths = arrayfun(@(k) sprintf('%stiles/tile_%d.tif', d, k), (1 : p.nF)', 'unif', 0);
[~, dom] = sort(strrep(p.dataOrder, ',', ''));
gf = '';
if p.gf > 0
    gf = [d 'groups.csv'];
end
[xs, dxs] = stitch_shift_assignment(paths, [d 'xcorr/'], p.imSizes, p.xyz, p.vs, dom, ...
    p.ov, p.regions, p.mo, p.ds, p.thr, p.tileIdx, p.method, p.s2d, p.aw, gf, p.large, ...
    p.pool, false, 1, false, '', 1);
out = {xs, dxs};
"""


def py_ssa(p, d):
    n = int(p["nF"])
    paths = [f"{d}tiles/tile_{k}.tif" for k in range(1, n + 1)]
    sizes = np.asarray(p["imSizes"])[:, [2, 0, 1]]
    pool = p["pool"]
    return list(
        sa.shift_assignment(
            paths,
            f"{d}xcorr/",
            sizes,
            p["xyz"],
            p["vs"],
            axis_string_from_matlab(p["dataOrder"]),
            p["ov"],
            p["regions"],
            np.asarray(p["mo"])[::-1],
            from_yxz(p["ds"]),
            p["thr"],
            p["tileIdx"],
            assign_method=p["method"],
            stitch_2d=p["s2d"],
            axis_weight=from_yxz(p["aw"]),
            group_file=f"{d}groups.csv" if p["gf"] else "",
            large_file=p["large"],
            pool_size=None if np.size(pool) == 0 else from_yxz(pool),
            node_factor=1,
        )
    )


def ssa_case(cid, make, **kw):
    return Case(
        cid,
        code=_SSA_CODE,
        args=[In("p"), InPath("d")],
        nargout=2,
        python=py_ssa,
        inputs={"p": Value(make), "d": Tree(ssa_tree(make), "d")},
        files=False,
        **kw,
    )


def _ssa_specs():
    out = []
    s = 9000
    for method in ("grid", "global", "local", "group", "test"):
        for name, gx, gy, gz, kw in [
            ("2x1", 2, 1, 1, {}),
            ("3x3", 3, 3, 1, {}),
            ("4x3jit", 4, 3, 1, {"jitter": 1.5}),
            ("3x2x2", 3, 2, 2, {}),
            ("5x4filt", 5, 4, 1, {"filt": 0.3}),
        ]:
            if method == "test":
                continue
            s += 1
            kw = dict(kw)
            if method == "group":
                kw["groups"] = 2 if gx > 1 else 1
                kw["group_file"] = 5
            out.append((f"ssa.{method}.{name}.s{s}", (s, gx, gy, gz, dict(kw, method=method))))
    extra = [
        # 'test' takes the pairs' values from the first columns of absolute_shift_mat (tile
        # numbers, then x) and leaves the constant to plain lsqlin: two tiles along y, x held at
        # 0 by its bounds, even x shifts, so that no value is a half integer
        ("test_1x2", 1, 2, 1, {"method": "test", "max_off": (0, 12, 5), "noise": 0.0,
                               "even": True, "filt": 0.0, "off": (1, 1, 1)}),
        ("grid_xyz_order", 3, 3, 1, {"data_order": "x,y,z"}),
        ("grid_zyx_order", 3, 2, 2, {"data_order": "z,y,x", "size_yxz": (40, 44, 36)}),
        ("grid_2d", 3, 3, 1, {"stitch2d": True, "size_yxz": (64, 72, 1), "overlap": (10, 12, 0)}),
        ("group_2d", 4, 3, 1, {"stitch2d": True, "size_yxz": (64, 72, 1), "overlap": (10, 12, 0),
                               "method": "group", "groups": 2, "group_file": 4}),
        ("grid_tight", 4, 4, 1, {"max_off": (3, 3, 2), "off": (8, 8, 3)}),
        ("grid_small_overlap", 4, 3, 1, {"overlap": (4, 5, 2), "max_off": (20, 20, 8)}),
        ("global_small_overlap", 4, 3, 1, {"overlap": (4, 5, 2), "max_off": (20, 20, 8),
                                           "method": "global"}),
        ("grid_large", 3, 3, 1, {"large_file": True}),
        ("global_large", 3, 2, 1, {"large_file": True, "method": "global"}),
        ("grid_weights", 4, 4, 1, {"aw": (2.0, 0.5, 3.0)}),
        ("grid_thr", 4, 4, 1, {"thr": 0.6}),
        ("group_locations", 4, 3, 1, {"method": "group",
                                      "locs": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]}),
        ("group_one", 3, 3, 1, {"method": "group"}),
        ("group_file4", 6, 3, 1, {"method": "group", "groups": 3, "group_file": 4}),
        ("group_jitter", 6, 3, 1, {"method": "group", "groups": 3, "group_file": 5,
                                   "jitter": 1.5}),
        ("grid_locations", 3, 2, 1, {"locs": [1, 1, 1, 2, 2, 2]}),
        ("local_jitter", 5, 4, 1, {"method": "local", "jitter": 2.0, "filt": 0.3}),
        ("grid_jitter_big", 8, 7, 1, {"jitter": 2.0}),
        # bounds between groups (twice the max offset) that hold
        ("group_tight", 6, 3, 1, {"method": "group", "groups": 3, "group_file": 5,
                                  "off": (10, 10, 4), "max_off": (6, 6, 3)}),
        # z shifts in 2D, held at 0 by the grid method's bounds
        ("grid_2d_zshift", 3, 3, 1, {"stitch2d": True, "size_yxz": (64, 72, 1),
                                     "overlap": (10, 12, 0), "zshift_2d": True}),
        # a group of tiles of two locations: pairs within it whose locations differ are corners
        ("group_mixed_locations", 4, 3, 1, {"method": "group", "groups": 2, "group_file": 5,
                                            "locs": [1, 2] * 6}),
    ]  # fmt: skip
    for name, gx, gy, gz, kw in extra:
        s += 1
        out.append((f"ssa.{name}.s{s}", (s, gx, gy, gz, kw)))
    return out


def _ssa_cases():
    cs = []
    for cid, (seed, gx, gy, gz, kw) in _ssa_specs():
        warns = "ids" if kw.get("locs") is not None or kw.get("method") == "group" else "record"
        cs.append(
            ssa_case(cid, lambda a=(seed, gx, gy, gz), kw=kw: ssa_problem(*a, **kw), warnings=warns)
        )
    return cs


# --- errors ------------------------------------------------------------------------------------


def _half_tile(p):
    """the first tile's index moved by a half along x and y: its pairs with its neighbors are
    one apart in all, along two axes"""
    p["tileIdx"][0, 0:2] += 0.5
    return p


def _error_cases():
    t2 = np.array([[0, 0, 0, 1], [0.5, 0.5, 0, 1]])
    p = mtest_problem(t2, [(1, 2)], [0.9], [[1, 2, 3]], [-50, -50, -10, 50, 50, 10])
    cs = [
        caller_case(
            "error.grid_two_axes_one_apart",
            "grid",
            lambda: p,
            expect_error=True,
            py_error=sa.AssignmentError,
        ),
        # the dense form up to 400 tiles switches on the axis (an error); the sparse form above
        # takes the first largest difference
        caller_case(
            "error.grid_two_axes_one_apart_400_tiles",
            "grid",
            lambda: _half_tile(grid_problem(9903, 20, 20, filt=0.0)),
            expect_error=True,
            py_error=sa.AssignmentError,
        ),
        caller_case(
            "grid.half_tile_index.401_tiles",
            "grid",
            lambda: _half_tile(grid_problem(9904, 20, 20, filt=0.0, isolated=1)),
        ),
        ssa_case(
            "error.ssa_test_pair_index",
            # pairs (1, 2) and (2, 3) above the threshold: the second's pair index, 3, is past
            # the 2 rows of absolute_shift_mat
            lambda: ssa_problem(9901, 3, 1, method="test", filt=0.0),
            expect_error=True,
            py_error=sa.AssignmentError,
        ),
        ssa_case(
            "error.ssa_unknown_method",
            lambda: ssa_problem(9902, 2, 1, method="mst"),
            expect_error=True,
            py_error=sa.AssignmentError,
        ),
    ]
    return cs


CASES = (
    [
        caller_case(cid, "grid", lambda a=spec: grid_problem(*a[:4], **a[4]))
        for cid, spec in _grid_specs()
    ]
    + [
        caller_case(cid, "global", lambda a=spec: grid_problem(*a[:4], **a[4]))
        for cid, spec in _global_specs()
    ]
    + [
        caller_case(cid, "local", lambda a=spec: grid_problem(*a[:4], **a[4]))
        for cid, spec in _local_specs()
    ]
    + [
        caller_case(cid, "group", lambda a=spec: group_problem(*a[:3], **a[3]))
        for cid, spec in _group_specs()
    ]
    + [
        caller_case(f"mtest.{name}.{method}", method, lambda p=p: p)
        for name, (p, methods, _) in _MTESTS.items()
        for method in methods
    ]
    + [
        wrapper_case("wrapper.grid_4x3", lambda: grid_problem(8501, 4, 3, filt=0.1)),
        wrapper_case("wrapper.grid_free", lambda: grid_problem(8502, 5, 4, free_tiles=1)),
        wrapper_case("wrapper.grid_line", lambda: grid_problem(8503, 1, 6, filt=0.2)),
        wrapper_case("wrapper.grid_weights", lambda: grid_problem(8504, 4, 4, aw=(2.0, 0.5, 3.0))),
        wrapper_case("wrapper.grid_allbelow", lambda: grid_problem(8505, 3, 3, filt=1.0)),
        wrapper_case("wrapper.existing", lambda: grid_problem(8506, 3, 2), existing=True),
    ]
    + _ssa_cases()
    + _error_cases()
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- StitchShiftAssignmentTest -----------------------------------------------------------------


def test_allPairsBelowTheThresholdKeepTheStagePositions():
    p = grid_problem(11, 2, 2, filt=1.0)
    n = int(p["nF"])
    assert np.array_equal(py_grid(p), np.zeros((n, 3)))
    assert np.array_equal(py_global(p), np.zeros((n, 3)))
    for v in (py_grid(p), py_global(p)):
        assert not np.signbit(v).any()


def test_onePairOfTiles():
    p, _, want = _MTESTS["one_pair"]
    assert np.array_equal(py_grid(p), want)
    assert np.array_equal(py_global(p), want)


def test_freeTilesKeepTheirStagePositionsRelativeToTheirNeighbors():
    # a 5 x 5 grid with exact shifts; tiles 3, 15 and 17 (1-based; no two neighbors) free; bounds
    # -20 .. 8
    g = 5
    ix, iy = np.meshgrid(np.arange(g), np.arange(g))
    t = np.column_stack([ix.ravel(order="F"), iy.ravel(order="F"), np.zeros(25), np.ones(25)])
    r = _rng(31, "free")
    true = np.column_stack([r.integers(-3, 3, 25), r.integers(-3, 3, 25), r.integers(-2, 2, 25)])
    pr = _neighbor_pairs(t)
    free = [2, 14, 16]
    xc = [0.1 if (i - 1 in free or j - 1 in free) else 0.9 for i, j in pr]
    rel = np.array([true[j - 1] - true[i - 1] for i, j in pr], dtype=np.float64)
    p = mtest_problem(t, pr, xc, rel, [-20, -20, -20, 8, 8, 8], (1, 0.1, 10))
    x = true.astype(np.float64)
    for k in free:
        nb = [j - 1 for i, j in pr if i - 1 == k] + [i - 1 for i, j in pr if j - 1 == k]
        x[k] = x[nb].mean(axis=0)
    x = x - x.mean(axis=0)
    assert np.min(np.abs(np.abs(x - np.floor(x)) - 0.5)) > 1e-3
    want = mround(x) - mround(x[0])
    assert np.array_equal(py_grid(p), want)
    assert np.array_equal(py_global(p), want)


def test_freeTileHeldByItsBounds():
    p, _, want = _MTESTS["free_tile_held_by_its_bounds"]
    assert np.array_equal(py_grid(p), want)


def test_boundsOfAFreeTileHoldItsNeighbors():
    p, _, _ = _MTESTS["bounds_of_a_free_tile_hold_its_neighbors"]
    assert np.array_equal(py_grid(p)[:, 0], [0, 2, 1])


def test_groupsWithoutPairsAboveTheThresholdKeepTheirStagePositions():
    want = {
        "groups_linked": [0, 2, 4, 2, 2, 6],
        "groups_unlinked": [0, 2, 2, 0, 0, 4],
        "groups_all_below": [0, 0, 0, 0, 0, 0],
    }
    for name, x in want.items():
        p, _, _ = _MTESTS[name]
        d = py_group(p)
        assert np.array_equal(d, np.column_stack([x, np.zeros((6, 2))])), name


def test_shiftsRoundTheExactOptimumAtAnyThreadCount(tmp_path):
    # positions 1e-7 above or below half integers; the same shifts at 1, 2 and 4 threads
    g = 10
    ix, iy = np.meshgrid(np.arange(g), np.arange(g))
    n = g * g
    t = np.column_stack([ix.ravel(order="F"), iy.ravel(order="F"), np.zeros(n), np.ones(n)])
    r = _rng(41, "half")
    v = r.integers(-10, 10, (n // 2, 3)) + 0.5 + (2.0 * (r.uniform((n // 2, 3)) > 0.5) - 1) * 1e-7
    x = np.vstack([v, -v])[np.argsort(r.uniform(n), kind="stable")]
    pr = _neighbor_pairs(t)
    rel = np.array([x[j - 1] - x[i - 1] for i, j in pr])
    p = mtest_problem(t, pr, [0.9] * len(pr), rel, [-100, -100, -100, 100, 100, 100], (1, 0.1, 10))
    want = mround(x) - mround(x[0])
    assert np.array_equal(py_grid(p), want)
    out = _in_subprocesses(tmp_path, "grid", p, (1, 2, 4))
    assert all(np.array_equal(o, want) for o in out)


def test_integerShiftsRoundAlike():
    p, _, pos = _MTESTS["integer_shifts_round_alike"]
    assert np.array_equal(py_grid(p), pos)
    assert np.array_equal(py_global(p), pos)


def test_activeBoundsGiveTheExactOptimum():
    p, _, _ = _MTESTS["active_bounds_give_the_exact_optimum"]
    u = 4.5 + 3e-8
    x = np.array([-2 * u, u, u]) / 3
    assert np.array_equal(py_grid(p)[:, 0], mround(x) - mround(x[0]))


def test_boundsHoldForTheirOwnPairs():
    # the first pair whose place in find(overlap_matrix)'s order differs from the list's,
    # measured 10 voxels, allows only 4
    ix, iy = np.meshgrid(np.arange(4.0), np.arange(3.0))
    t = np.column_stack([ix.ravel(order="F"), iy.ravel(order="F"), np.zeros(12), np.ones(12)])
    pr = np.array(_neighbor_pairs(t))
    ov = np.zeros((12, 12), bool)
    ov[pr[:, 0] - 1, pr[:, 1] - 1] = True
    fj, fi = np.nonzero(ov.T)
    k = int(np.nonzero((fi + 1 != pr[:, 0]) | (fj + 1 != pr[:, 1]))[0][0])
    rel = np.zeros((len(pr), 3))
    rel[k, 0] = 10
    mas = np.tile([-50.0, -50, -50, 50, 50, 50], (len(pr), 1))
    mas[k, 3] = 4
    p = mtest_problem(t, pr, [0.9] * len(pr), rel, mas)
    for d in (py_grid(p), py_global(p)):
        shift = d[pr[:, 1] - 1] - d[pr[:, 0] - 1]
        assert shift[k, 0] == 4
        assert np.all((shift >= mas[:, 0:3]) & (shift <= mas[:, 3:6]))


def _sub_problem(p, mask):
    """the problem of the tiles in mask alone (renumbered)"""
    inds = np.nonzero(mask)[0]
    new = np.zeros(int(p["nF"]))
    new[inds] = np.arange(1, inds.size + 1)
    mx, ash, pairs = p["mx"], p["ash"], p["pairs"].astype(np.int64)
    keep = mask[pairs[:, 0] - 1] & mask[pairs[:, 1] - 1]
    q = dict(p)
    q["nF"] = float(inds.size)
    q["mx"] = np.column_stack([new[pairs[keep] - 1], mx[keep, 2]])
    q["ash"] = np.column_stack([new[pairs[keep] - 1], ash[keep, 2:5]])
    q["pairs"] = new[pairs[keep] - 1]
    q["mas"] = p["mas"][keep]
    q["tileIdx"] = p["tileIdx"][inds]
    return q, inds


def test_separateGroupsOfALargeGridHaveMeanZeroEach():
    p = grid_problem(51, 15, 14, two_groups=True, filt=0.05, max_off=(100, 100, 30))
    d = py_grid(p)
    n1 = 15 * 14
    groups = np.r_[np.zeros(n1, bool), np.ones(n1, bool)]
    for mask in (~groups, groups):
        q, inds = _sub_problem(p, mask)
        dg = py_grid(q)
        assert np.array_equal(d[inds] - d[inds[0]], dg - dg[0])


def test_largeGridMatchesDenseSolver(monkeypatch):
    # the sparse form's shifts are the dense form's (the same solver); both find the simulated
    # positions within 3 voxels (noise 0.7)
    p, true = grid_problem(61, 21, 21, filt=0.05, return_true=True)
    sparse = py_grid(p)
    monkeypatch.setattr(sa, "SPARSE_TILE_NUM", 10**6)
    dense = py_grid(p)
    assert np.array_equal(sparse, dense)
    assert np.abs(sparse - (true - true[0])).max() <= 3


# --- StitchingTest: the assignment of registered tiles ----------------------------------------


def _write_results(xcorr_dir, results):
    """registration results per leading tile (1-based pairs) as S04's tasks write them"""
    os.makedirs(xcorr_dir, exist_ok=True)
    by = {}
    for (i, j), (rel, xc) in results.items():
        by.setdefault(i, []).append((j, rel, xc))
    for i, rows in by.items():
        rows.sort()
        write_mat(
            os.path.join(xcorr_dir, f"xcorr_tile_{i}.mat"),
            {
                "relative_shift_mat": np.array([r[1] for r in rows], dtype=np.float64),
                "max_xcorr_mat": np.array([[r[2]] for r in rows], dtype=np.float64),
                "pair_indices": np.array([[i, r[0]] for r in rows], dtype=np.float64),
            },
        )


def _registered(tmp_path, stage_xyz, true_xyz, tile_idx, size_zyx, method, large=False):
    """shift_assignment for tiles at stage positions whose data sit at true positions (um): the
    registration results a correct registration gives (the overlap's shift in voxels, signed by
    the tiles' order)"""
    vs = np.array([0.108, 0.108, 0.25])
    n = len(stage_xyz)
    ext = np.array([size_zyx[2], size_zyx[1], size_zyx[0]], dtype=np.float64)
    cub = np.hstack([stage_xyz, stage_xyz + (ext - 1) * vs])
    om, regions = geometry.find_overlaps(cub)
    paths = [str(tmp_path / f"t{k}.tif") for k in range(n)]
    synth.write_tiff(paths[0], np.ones(size_zyx, np.uint16))
    xd = str(tmp_path / "xcorr")
    res = {}
    for i, j in zip(*np.nonzero(om), strict=True):
        absolute = mround((true_xyz[j] - stage_xyz[j] - (true_xyz[i] - stage_xyz[i])) / vs)
        flag = (0.5 - (stage_xyz[i] > stage_xyz[j])) * 2
        res[int(i) + 1, int(j) + 1] = (absolute * flag, 0.9)
    if large:
        os.makedirs(xd, exist_ok=True)
        for (i, j), (rel, xc) in res.items():
            write_mat(
                os.path.join(xd, f"xcorr_tile_{i}_{j}.mat"),
                {
                    "relative_shift_mat": rel[None, :],
                    "max_xcorr_mat": np.array([[xc]]),
                    "pair_indices": np.array([[i, j]], dtype=np.float64),
                },
            )
        mip = tmp_path / "MIP_Slabs_pooling_2_2_2"
        mip.mkdir()
        for k in range(n):
            (mip / f"t{k}_z.zarr").mkdir()
    else:
        _write_results(xd, res)
    xs, dxs = sa.shift_assignment(
        paths, xd, np.tile(size_zyx, (n, 1)), stage_xyz, vs, "z,y,x", om, regions,
        (10, 20, 20), (1, 1, 1), 0.25, tile_idx, assign_method=method, large_file=large,
        pool_size=(2, 2, 2) if large else None,
    )  # fmt: skip
    return xs, np.asarray(mround(dxs / vs))


def test_xcorrRecoversPerturbedTilePosition(tmp_path):
    # the second tile is listed 4 voxels off in x
    px = 0.108
    stage = np.array([[0, 0, 0], [44 * px, 0, 0]])
    true = np.array([[0, 0, 0], [40 * px, 0, 0]])
    tidx = np.array([[0, 0, 0, 1], [1, 0, 0, 1]], dtype=np.float64)
    for method in ("grid", "global", "local"):
        d = tmp_path / method
        d.mkdir()
        _, shifts = _registered(d, stage, true, tidx, (32, 64, 64), method)
        assert np.array_equal(shifts, [[0, 0, 0], [-4, 0, 0]]), method


def test_largeFileRegistrationRecoversPerturbedTilePosition(tmp_path):
    # the third tile is listed 2 voxels off in x
    px = 0.108
    stage = np.array([[0, 0, 0], [0, 36 * px, 0], [50 * px, 0, 0], [48 * px, 36 * px, 0]])
    true = np.array([[0, 0, 0], [0, 36 * px, 0], [48 * px, 0, 0], [48 * px, 36 * px, 0]])
    tidx = np.array([[0, 0, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1], [1, 1, 0, 1]], dtype=np.float64)
    _, shifts = _registered(tmp_path, stage, true, tidx, (32, 64, 72), "grid", large=True)
    rows = sorted(map(tuple, shifts.tolist()))
    assert rows == [(-2, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]


def test_localShiftMethod(tmp_path):
    px = 0.108
    stage = np.array([[0, 0, 0], [44 * px, 0, 0], [88 * px, 0, 0]])
    true = np.array([[0, 0, 0], [40 * px, 0, 0], [83 * px, 0, 0]])
    tidx = np.array([[0, 0, 0, 1], [1, 0, 0, 1], [2, 0, 0, 1]], dtype=np.float64)
    xs, shifts = _registered(tmp_path, stage, true, tidx, (32, 64, 64), "local")
    # tile 2 is listed 4 voxels off, tile 3 5 voxels
    assert np.array_equal(shifts, [[0, 0, 0], [-4, 0, 0], [-5, 0, 0]])
    assert np.array_equal(xs, stage + shifts * np.array([0.108, 0.108, 0.25]))


# --- the parts of shift_assignment ------------------------------------------------------------


class _Calls:
    def __init__(self):
        self.calls = []

    def __call__(self, inputs, outputs, tasks, **kw):
        self.calls.append((list(inputs), list(outputs), list(tasks), kw))
        return [os.path.exists(o) for o in outputs]


def test_registration_tasks_per_leading_tile(tmp_path, monkeypatch):
    # tasks of tiles with pairs only, pairs in find's order, cuboids as mat2str gives them (15
    # digits), 2, 4 and 8 CPUs, and a missing result is MATLAB's error
    from petakit5d.jobs import framework

    p = ssa_problem(77, 3, 2, jitter=1.37)
    xyz = p["xyz"] + 1.0 / 3.0
    calls = _Calls()
    monkeypatch.setattr(framework, "run_tasks", calls)
    synth.write_tiff(tmp_path / "t0.tif", np.ones((20, 64, 72), np.uint16))
    paths = [str(tmp_path / f"t{k}.tif") for k in range(6)]
    with pytest.raises(sa.AssignmentError, match="Some xcorr files are missing!"):
        sa.shift_assignment(
            paths, str(tmp_path / "x"), np.tile((20, 64, 72), (6, 1)), xyz, p["vs"], "z,y,x",
            p["ov"], p["regions"], (5, 12, 12), (1, 1, 1), 0.25, p["tileIdx"],
        )  # fmt: skip
    assert [c[3]["cpus_per_task"] for c in calls.calls] == [2, 4, 8]
    inputs, outputs, tasks, kw = calls.calls[0]
    om = p["ov"].copy()
    ti, tj = np.nonzero(om.T)[::-1]
    corner = np.abs(p["tileIdx"][ti, :3] - p["tileIdx"][tj, :3]).sum(axis=1) != 1
    om[ti[corner], tj[corner]] = False
    lead = sorted({int(i) for i in np.nonzero(om)[0]})
    assert inputs == [paths[f] for f in lead]
    assert outputs == [os.path.join(str(tmp_path / "x"), f"xcorr_tile_{f + 1}.mat") for f in lead]
    rec = tasks[0].record()
    assert rec["func"] == sa.REGISTRATION_TASK
    f = lead[0]
    js = sorted(int(j) for j in np.nonzero(om[f])[0])
    assert rec["args"][3] == [[f, j] for j in js]
    assert rec["args"][1] == [paths[j] for j in js]
    exact = np.r_[xyz[f], xyz[f] + (np.array([72, 64, 20]) - 1) * p["vs"]]
    assert rec["args"][4] == through_mat2str(exact).tolist()
    assert rec["args"][4] != exact.tolist()  # 1/3 has more than 15 digits
    assert rec["kwargs"]["max_offset"] == [5.0, 12.0, 12.0]
    assert kw["min_task_job_num"] == 1


def test_large_file_tasks_per_pair(tmp_path, monkeypatch):
    from petakit5d.jobs import framework

    p = ssa_problem(78, 3, 1, large_file=True)
    calls = _Calls()
    monkeypatch.setattr(framework, "run_tasks", calls)
    synth.write_tiff(tmp_path / "t0.tif", np.ones((20, 64, 72), np.uint16))
    paths = [str(tmp_path / f"t{k}.tif") for k in range(3)]
    with pytest.raises(sa.AssignmentError):
        sa.shift_assignment(
            paths, str(tmp_path / "x"), np.tile((20, 64, 72), (3, 1)), p["xyz"], p["vs"],
            "z,y,x", p["ov"], p["regions"], (5, 12, 12), (1, 1, 1), 0.25, p["tileIdx"],
            large_file=True, pool_size=(2, 3, 4),
        )  # fmt: skip
    mip_call, reg_call = calls.calls[0], calls.calls[3]
    assert mip_call[1] == [str(tmp_path / f"MIP_Slabs_pooling_3_4_2/t{k}_z.zarr") for k in range(3)]
    assert (tmp_path / "MIP_Slabs_pooling_3_4_2").is_dir()
    assert reg_call[1] == [os.path.join(str(tmp_path / "x"), f"xcorr_tile_{a}_{a + 1}.mat")
                           for a in (1, 2)]  # fmt: skip
    assert reg_call[3]["min_task_job_num"] == 2
    assert reg_call[2][0].record()["kwargs"]["mip_dir_str"] == "MIP_Slabs_pooling_3_4_2"


def test_xcorr_result_paths_in_subfolders_above_10000_tiles(tmp_path):
    paths = sa.xcorr_result_paths(str(tmp_path), 10001)
    assert paths[0] == os.path.join(str(tmp_path), "xcorr_tiles_000001_010000", "xcorr_tile_1.mat")
    assert paths[-1] == os.path.join(
        str(tmp_path), "xcorr_tiles_010001_010001", "xcorr_tile_10001.mat"
    )
    assert sorted(os.listdir(tmp_path)) == [
        "xcorr_tiles_000001_010000",
        "xcorr_tiles_010001_010001",
    ]
    assert sa.xcorr_result_paths("x", 2) == [os.path.join("x", "xcorr_tile_1.mat"),
                                             os.path.join("x", "xcorr_tile_2.mat")]  # fmt: skip


def test_single_tile_has_no_tasks_and_no_shift(tmp_path):
    synth.write_tiff(tmp_path / "t.tif", np.ones((4, 6, 5), np.uint16))
    xs, dxs = sa.shift_assignment(
        [str(tmp_path / "t.tif")], str(tmp_path / "x"), [(4, 6, 5)], [[0.0, 0.0, 0.0]],
        (0.1, 0.1, 0.2), "z,y,x", np.zeros((1, 1), bool), np.zeros((0, 6)), (5, 10, 10),
        (1, 1, 1), 0.25, [[0, 0, 0, 1]],
    )  # fmt: skip
    assert np.array_equal(xs, np.zeros((1, 3))) and np.array_equal(dxs, np.zeros((1, 3)))


def test_group_cluster_tasks_equal_the_local_solve(tmp_path, monkeypatch):
    # more than 100 tiles in groups: the groups' grid problems go through .mat files and
    # global_grid_assignment tasks (the threshold as MATLAB's %d gives it), with the shifts of
    # the local solve
    from petakit5d.jobs import framework

    monkeypatch.chdir(tmp_path)
    p = group_problem(91, 15, 8, bands=4, filt=0.1, thr=0.2512345678)
    local = py_group(p)
    real = framework.run_tasks
    seen = []

    def run(inputs, outputs, tasks, **kw):
        seen.append([t.record() for t in tasks])
        kw["parse_cluster"] = False
        return real(inputs, outputs, tasks, **kw)

    monkeypatch.setattr(framework, "run_tasks", run)
    n, mx, ash, om, mas, thr = py_problem(p)
    d = sa.group_assignment(
        n, mx, ash, om, mas, thr, p["tileIdx"], from_yxz(p["aw"]), np.ravel(p["grp"]),
        xcorr_dir=str(tmp_path), parse_cluster=True,
    )  # fmt: skip
    rec = seen[0][0]
    assert rec["func"] == "petakit5d.stitch.assignment:global_grid_assignment"
    assert rec["args"][2] == 0.2512346
    opt = tmp_path / "xcorr_group_optimization"
    assert sorted(f.name for f in opt.iterdir()) == sorted(
        [f"xcorr_group_{g}_opt_information.mat" for g in range(1, 5)]
        + [f"xcorr_group_{g}_opt_result.mat" for g in range(1, 5)]
    )
    info = read_mat(opt / "xcorr_group_1_opt_information.mat")
    assert info["max_xcorr_mat_g"][:, :2].min() == 1
    # the threshold rounded to 7 digits changes nothing here
    assert np.array_equal(d, local)


def test_group_of_one_tile_in_the_cluster_path_errors_as_matlab(tmp_path):
    p = group_problem(92, 15, 8, bands=4, islands=1)
    n, mx, ash, om, mas, thr = py_problem(p)
    with pytest.raises(sa.AssignmentError, match="group of one tile"):
        sa.group_assignment(
            n, mx, ash, om, mas, thr, p["tileIdx"], from_yxz(p["aw"]), np.ravel(p["grp"]),
            xcorr_dir=str(tmp_path), parse_cluster=True,
        )  # fmt: skip


def test_wrapper_keeps_an_existing_result(tmp_path, capsys):
    p = grid_problem(81, 3, 3)
    write_fixed_mat(tmp_path / "info.mat", group_info(p))
    out = tmp_path / "res.mat"
    sa.global_grid_assignment(str(tmp_path / "info.mat"), str(out), 0.25, (10, 1, 0.1), "u")
    d = read_mat(out)["d_shift_g"]
    assert np.array_equal(d, py_grid(p))
    write_mat(out, {"d_shift_g": np.ones((1, 3))})
    sa.global_grid_assignment(str(tmp_path / "info.mat"), str(out), 0.25, (10, 1, 0.1), "u")
    assert "already exists, skip it!" in capsys.readouterr().out
    assert np.array_equal(read_mat(out)["d_shift_g"], np.ones((1, 3)))
    assert sorted(f.name for f in tmp_path.iterdir()) == ["info.mat", "res.mat"]


def test_wrapper_is_registered_with_matlab_values(tmp_path):
    from petakit5d import matlab

    p = grid_problem(82, 4, 2, aw=(2.0, 0.5, 3.0))
    write_fixed_mat(tmp_path / "info.mat", group_info(p))
    out = tmp_path / "res.mat"
    matlab.stitch_global_grid_assignment_wrapper(
        str(tmp_path / "info.mat"), str(out), 0.25, [2.0, 0.5, 3.0], "u"
    )
    assert np.array_equal(read_mat(out)["d_shift_g"], py_grid(p))


def test_grid_weights_dense_and_sparse_rules():
    t = np.array([[0, 0, 0, 1], [1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [1, 1, 0, 1]], float)
    s, q = np.array([0, 0, 0, 0]), np.array([1, 2, 3, 4])
    aw = (10.0, 1.0, 0.1)  # z, y, x
    for dense in (True, False):
        w = sa._grid_weights(t, s, q, aw, dense)
        assert np.array_equal(w, [0.1, 1.0, 10.0, 0.0])
    t2 = np.array([[0, 0, 0, 1], [0.5, 0.5, 0, 1]])
    assert np.array_equal(sa._grid_weights(t2, [0], [1], aw, False), [0.1])
    with pytest.raises(sa.AssignmentError) as e:
        sa._grid_weights(t2, [0], [1], aw, True)
    assert e.value.identifier.startswith("MATLAB:")


def test_signs_and_bounds_follow_the_tiles_order(tmp_path):
    # a pair whose second tile lies before the first along x: its shift and bounds turn around
    vs = np.array([0.108, 0.108, 0.25])
    xyz = np.array([[5.0, 0.0, 0.0], [0.0, 6.0, 0.0]])
    size = (20, 64, 72)
    ext = np.array([72, 64, 20], float)
    cub = np.hstack([xyz, xyz + (ext - 1) * vs])
    om, regions = geometry.find_overlaps(cub)
    assert om[0, 1]
    synth.write_tiff(tmp_path / "t0.tif", np.ones(size, np.uint16))
    _write_results(str(tmp_path / "x"), {(1, 2): (np.array([3.0, 2.0, 1.0]), 0.9)})
    tidx = np.array([[0, 0, 0, 1], [0, 1, 0, 1]], float)
    _, dxs = sa.shift_assignment(
        [str(tmp_path / "t0.tif"), str(tmp_path / "t1.tif")], str(tmp_path / "x"),
        [size, size], xyz, vs, "z,y,x", om, regions, (5, 50, 50), (1, 1, 1), 0.25, tidx,
        assign_method="global",
    )  # fmt: skip
    d = mround(dxs / vs)
    # x: tile 1 ahead of tile 2, so the shift's sign flips; y and z keep theirs
    assert np.array_equal(d[1] - d[0], [-3, 2, 1])


def test_warnings_of_the_pair_selection(tmp_path):
    p = ssa_problem(79, 2, 1)
    tidx = p["tileIdx"].copy()
    tidx[1, 3] = 2
    with pytest.warns(sa.AssignmentWarning, match="not from the same data folder"):
        _ssa_python(tmp_path / "a", p, tile_idx=tidx)
    with pytest.warns(sa.AssignmentWarning, match="treat all tiles as one group"):
        _ssa_python(tmp_path / "b", dict(p, method="group"))


def _ssa_python(folder, p, tile_idx=None):
    folder = Path(folder)
    ssa_tree(lambda: p)(folder)
    q = dict(p)
    if tile_idx is not None:
        q["tileIdx"] = tile_idx
    return py_ssa(q, f"{folder}/")


def test_results_do_not_depend_on_the_thread_count(tmp_path):
    probs = {
        "grid": grid_problem(101, 12, 11, free_tiles=3, filt=0.1),
        "sparse": grid_problem(102, 25, 20, free_tiles=4, filt=0.1, max_off=(6, 6, 3)),
        "group": group_problem(103, 12, 6, bands=3, free_groups=(1,)),
    }
    for name, p in probs.items():
        method = "group" if name == "group" else "grid"
        out = _in_subprocesses(tmp_path / name, method, p, (1, 4))
        assert np.array_equal(out[0], out[1]), name
        want = py_group(p) if method == "group" else py_grid(p)
        assert np.array_equal(out[0], want), name


def _in_subprocesses(folder, method, p, threads):
    """the shifts of p by method in fresh interpreters at each thread count"""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    np.savez(folder / "p.npz", **{k: np.asarray(v) for k, v in p.items()})
    script = (
        "import sys, numpy as np\n"
        "sys.path.insert(0, sys.argv[3])\n"
        "import test_assignment as t\n"
        "z = np.load(sys.argv[1])\n"
        "p = {k: (z[k].item() if z[k].ndim == 0 else z[k]) for k in z.files}\n"
        "d = t._PY[sys.argv[2]](p)\n"
        "np.save(sys.argv[1] + '.' + sys.argv[4] + '.npy', d)\n"
    )
    out = []
    here = os.path.dirname(os.path.abspath(__file__))
    for n in threads:
        env = dict(os.environ)
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            env[v] = str(n)
        env["PETAKIT5D_NUM_THREADS"] = str(n)
        subprocess.run(
            [sys.executable, "-c", script, str(folder / "p.npz"), method, here, str(n)],
            check=True,
            env=env,
        )
        out.append(np.load(folder / f"p.npz.{n}.npy"))
    return out


def test_xcorr_results_read_as_load_reads_them(tmp_path):
    # results of several pairs, in v7.3 (as MATLAB and the port write them) and in v6
    rel = np.array([[1.5, -2.0, 3.0], [4.0, 5.25, -6.0], [7.0, 8.0, 9.5]])
    v = {
        "relative_shift_mat": rel,
        "max_xcorr_mat": np.array([[0.9], [0.1], [0.5]]),
        "pair_indices": np.array([[2.0, 3.0], [2.0, 5.0], [2.0, 9.0]]),
    }
    for version in ("7.3", "6"):
        path = str(tmp_path / f"x{version}.mat")
        write_mat(path, v, version=version)
        pairs, r, mxv = sa._read_xcorr(path)
        assert np.array_equal(pairs, [[1, 2], [1, 4], [1, 8]])
        assert np.array_equal(r, rel) and np.array_equal(mxv, [0.9, 0.1, 0.5])


def test_test_method_rounds_without_the_offset():
    # 'test' rounds plain lsqlin's solution with plain round: tiles at -0.5 and 0.5 end 2 apart
    tidx = np.array([[0, 0, 0, 1], [1, 0, 0, 1]], dtype=np.float64)
    mx = np.array([[0, 1, 0.9]])
    ash = np.array([[0, 1, 1.0, 0, 0]])
    om = np.array([[False, True], [False, False]])
    mas = np.array([[0.0, -12, -5, 0, 12, 5]])
    d = sa.grid_assignment_test(2, mx, ash, om, mas, 0.25, tidx)
    # x: tile number 1, held at 0 by its bounds; y: tile number 2; z: the x shift, 1
    assert np.array_equal(d, [[0, 0, 0], [0, 2, 2]])


def test_free_group_steps_use_the_modified_overlaps():
    # a free group overlapping two placed groups takes the mean of its pairs' differences to
    # them: pairs within groups are left out of step 3
    p, _, _ = _MTESTS["groups_unlinked"]
    n, mx, ash, om, mas, thr = py_problem(p)
    d = sa.group_assignment(n, mx, ash, om, mas, thr, p["tileIdx"], (1, 1, 1), np.ravel(p["grp"]))
    assert np.array_equal(d[:, 0], [0, 2, 2, 0, 0, 4])


def test_default_axis_weight_is_matlabs():
    assert sa.AXIS_WEIGHT == tuple(from_yxz([1, 0.1, 10]))
    assert sa.METHODS == ("local", "global", "grid", "test", "group")
    assert lsq.ROUND_OFFSET == 1e-9
