"""The stitching shift assignment: the tiles' shifts from the registration of their overlapping
pairs (``stitch_shift_assignment`` and the methods it calls).

==========================  =======================================================================
Python                      MATLAB
==========================  =======================================================================
``shift_assignment``        stitch_shift_assignment: the pairs of each method, the registration
                            tasks (one per leading tile, one per pair for large files), their result
                            files, the signs, the bounds and the assignment
``local_assignment``        stitch_local_assignment ('local': a spanning forest of the strongest
                            pairs, the shifts added up along it from tile 1)
``global_assignment``       stitch_global_assignment ('global': pairs weighted by their correlation)
``grid_assignment``         stitch_global_grid_assignment ('grid', the default: neighbors weighted
                            by their axis; sparse matrices above 400 tiles)
``group_assignment``        stitch_global_group_assignment ('group': the grid method in each group,
                            then the groups by their pairs above the threshold)
``grid_assignment_test``    stitch_shift_assignment's stitch_global_grid_assignment_test ('test')
``global_grid_assignment``  stitch_global_grid_assignment_wrapper: a group's grid assignment as a
                            task, from and to .mat files
==========================  =======================================================================

Conventions (S03's for the stitching geometry):

- Tiles are numbered from 0. The pair tables keep MATLAB's columns: ``max_xcorr_mat`` ``(p, 3)``
  ``[i, j, max xcorr]``, ``absolute_shift_mat`` ``(p, 5)`` ``[i, j, x, y, z]``,
  ``max_allow_shift`` ``(p, 6)`` (the lower bounds x, y, z, then the upper ones), shifts in voxels,
  the rows by first, then second tile; ``overlap_matrix`` is ``(n, n)`` bool, ``[i, j]`` for
  ``i < j``. ``d_shift`` is ``(n, 3)`` x, y, z voxels; ``xyz``, ``xyz_shift`` and ``dxyz_shift``
  ``(n, 3)`` x, y, z um; ``overlap_regions`` the physical overlap cuboids of
  ``geometry.find_overlaps``; ``tile_idx`` ``(n, 4)`` the columns x, y, z, location.
- Per-axis vectors are ``(z, y, x)``: ``axis_weight`` (MATLAB's ``[y x z]``; default
  ``(10, 1, 0.1)``), ``max_offset`` (MATLAB's ``MaxOffset``, ``[x y z]``), ``xcorr_downsample``
  and ``pool_size`` (MATLAB's ``[y x z]``). Tile sizes are the arrays' ``(z, y, x)`` with
  ``data_order`` naming their axes.
- Files keep MATLAB's form: the registration results ``xcorr_tile_<i>.mat`` (``relative_shift_mat``,
  ``max_xcorr_mat``, 1-based ``pair_indices``) and the group method's optimization files.

The grid, global and group methods solve with ``matlabcompat.lsq.assignment`` (MATLAB's
``stitch_assignment_lsqlin``, the three axes in one call), then round with ``round(x + 1e-9)`` and
subtract tile 1's value, as MATLAB does, so the shifts do not depend on the thread count.
"""

from __future__ import annotations

import math
import os
import time
import warnings
from collections.abc import Sequence

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "AXIS_WEIGHT",
    "METHODS",
    "SPARSE_TILE_NUM",
    "XCORR_BATCH_SIZE",
    "AssignmentError",
    "AssignmentWarning",
    "global_assignment",
    "global_grid_assignment",
    "grid_assignment",
    "grid_assignment_test",
    "group_assignment",
    "local_assignment",
    "shift_assignment",
    "xcorr_result_paths",
]

# MATLAB's default axisWeight [1, 0.1, 10] (y, x, z)
AXIS_WEIGHT = (10.0, 1.0, 0.1)
METHODS = ("local", "global", "grid", "test", "group")
# the grid method solves with sparse matrices above this number of tiles
SPARSE_TILE_NUM = 400
# registration results of more tiles go to subfolders of this many tiles
XCORR_BATCH_SIZE = 10000
# the registration task of S04 (stitch_shift_assignment's cross_correlation_registration_wrapper)
REGISTRATION_TASK = "petakit5d.stitch.registration:cross_correlation_registration"


class AssignmentError(ValueError):
    """An error MATLAB raises here; ``identifier`` is its error id (``''`` for none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


class AssignmentWarning(UserWarning):
    """A warning MATLAB gives here; ``identifier`` is its warning id (``''`` for none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


def _warn(message: str, identifier: str = "") -> None:
    warnings.warn(AssignmentWarning(message, identifier), stacklevel=3)


# ------------------------------------------------------------------------------------------------
# tables and systems
# ------------------------------------------------------------------------------------------------


def _table(m, ncol: int, name: str) -> np.ndarray:
    a = np.array(m, dtype=np.float64)
    if a.size == 0:
        return np.zeros((0, ncol))
    a = a.reshape(-1, ncol) if a.ndim < 2 else a
    if a.ndim != 2 or a.shape[1] != ncol:
        raise AssignmentError(f"{name} must have {ncol} columns")
    return a


def _tiles(a: np.ndarray) -> np.ndarray:
    """tile numbers of a table column as indices"""
    return a.astype(np.int64)


def _overlap(overlap_matrix, n: int) -> np.ndarray:
    om = np.array(overlap_matrix, dtype=bool)
    if om.size == 0:
        return np.zeros((n, n), bool)
    if om.shape != (n, n):
        raise AssignmentError(f"overlap_matrix must be {n} x {n}")
    return om


def _find(om: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``[i, j] = find(overlap_matrix)``: the pairs in MATLAB's column-major order (by second,
    then first tile)"""
    j, i = np.nonzero(om.T)
    return i.astype(np.int64), j.astype(np.int64)


def _sorted_pairs(om: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``find(overlap_matrix)`` sorted by first, then second tile (the rows of max_xcorr_mat)"""
    i, j = np.nonzero(om)
    return i.astype(np.int64), j.astype(np.int64)


def _differences(ci, cj, n: int, values=None):
    """the sparse rows ``v (e_j - e_i)``, one per pair (CSR with sorted columns), without zero
    entries (as MATLAB's sparse keeps them: rows of weight 0, and pairs of a tile with itself,
    have none)"""
    import scipy.sparse as sp

    ci = np.asarray(ci, dtype=np.int64).ravel()
    cj = np.asarray(cj, dtype=np.int64).ravel()
    m = ci.size
    v = np.ones(m) if values is None else np.asarray(values, dtype=np.float64).ravel()
    keep = (v != 0) & (ci != cj)
    at_lo = np.where(ci < cj, -v, v)
    indptr = np.zeros(m + 1, dtype=np.int64)
    np.cumsum(np.where(keep, 2, 0), out=indptr[1:])
    indices = np.column_stack([np.minimum(ci, cj), np.maximum(ci, cj)])[keep].ravel()
    data = np.column_stack([at_lo, -at_lo])[keep].ravel()
    return sp.csr_matrix((data, indices, indptr), shape=(m, n))


def _bounds_matrix(ni, nj, n: int):
    """``A = [R; -R]`` of the pairs (ni, nj)"""
    ni = np.asarray(ni, dtype=np.int64)
    nj = np.asarray(nj, dtype=np.int64)
    one = np.ones(ni.size)
    return _differences(np.r_[ni, ni], np.r_[nj, nj], n, np.r_[one, -one])


def _bound_rows(mas: np.ndarray) -> np.ndarray:
    """``b = [u; -l]`` of the three axes, one column each"""
    return np.vstack([mas[:, 3:6], -mas[:, 0:3]])


def _bounded_system(n, ni, nj, C, D, mas):
    """``[R; -R]`` of the pairs (ni, nj) and their bounds, solved along the three axes"""
    from petakit5d.matlabcompat import lsq

    return lsq.assignment(C, D, _bounds_matrix(ni, nj, n), _bound_rows(mas))


def _round_to_first(x: np.ndarray, offset: bool = True) -> np.ndarray:
    """``round(d_shift + 1e-9) - (its first row)`` (``offset`` False: plain ``round``)"""
    from petakit5d.core.dtypes import mround
    from petakit5d.matlabcompat.lsq import ROUND_OFFSET

    x = np.asarray(x, dtype=np.float64)
    d = mround(x + ROUND_OFFSET) if offset else mround(x)
    if d.shape[0] == 0:
        raise AssignmentError(
            "Index in position 1 exceeds array bounds. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    return d - d[0]


def _sum3(a: np.ndarray) -> np.ndarray:
    """``sum(a, 2)`` of rows of three, in index order"""
    return (a[:, 0] + a[:, 1]) + a[:, 2]


# ------------------------------------------------------------------------------------------------
# the methods
# ------------------------------------------------------------------------------------------------


def local_assignment(
    n_tiles: int, max_xcorr_mat, absolute_shift_mat, overlap_matrix, xcorr_thresh: float
) -> np.ndarray:
    """``stitch_local_assignment``: the minimum spanning forest of the pairs weighted by minus
    their correlation (pairs below ``xcorr_thresh`` at 0.001), searched depth first from tile 0;
    each tile's shift is its predecessor's plus their pair's (signed by the pair's order). Tiles
    the search doesn't reach keep 0. ``overlap_matrix`` is not used. Not rounded."""
    from petakit5d.matlabcompat import graph as mg

    del overlap_matrix
    n = int(n_tiles)
    mx = _table(max_xcorr_mat, 3, "max_xcorr_mat")
    ash = _table(absolute_shift_mat, 5, "absolute_shift_mat")
    d = np.zeros((n, 3))
    w = mx[:, 2].copy()
    w[w < xcorr_thresh] = 0.001
    G = mg.graph(_tiles(mx[:, 0]), _tiles(mx[:, 1]), -w)
    T, _ = mg.minspantree(G, type="forest")
    v = mg.dfsearch(T, 0)
    # each tile's predecessor, its only neighbor in the tree visited before it
    order = np.full(T.numnodes, -1, dtype=np.int64)
    order[v] = np.arange(v.size)
    e = T.edges
    reached = (order[e[:, 0]] >= 0) & (order[e[:, 1]] >= 0)
    e = e[reached]
    later = order[e[:, 1]] > order[e[:, 0]]
    child = np.where(later, e[:, 1], e[:, 0])
    pred = np.full(T.numnodes, -1, dtype=np.int64)
    pred[child] = np.where(later, e[:, 0], e[:, 1])
    # the first row of absolute_shift_mat of each pair (s, t), s < t
    si, ti = _tiles(ash[:, 0]), _tiles(ash[:, 1])
    m = int(max(T.numnodes, n, si.max(initial=0) + 1, ti.max(initial=0) + 1))
    keys, first = np.unique(si * m + ti, return_index=True)
    vals = ash[:, 2:5].tolist()
    out = d.tolist()
    for ni in v[1:].tolist():
        npr = int(pred[ni])
        if ni > npr:
            s, t, sign = npr, ni, 1.0
        else:
            s, t, sign = ni, npr, -1.0
        row = vals[int(first[np.searchsorted(keys, s * m + t)])]
        p = out[npr]
        out[ni] = [p[0] + sign * row[0], p[1] + sign * row[1], p[2] + sign * row[2]]
    return np.array(out, dtype=np.float64).reshape(n, 3)


def global_assignment(
    n_tiles: int,
    max_xcorr_mat,
    absolute_shift_mat,
    overlap_matrix,
    max_allow_shift,
    xcorr_thresh: float,
) -> np.ndarray:
    """``stitch_global_assignment``: the least squares of the pairs at or above ``xcorr_thresh``
    weighted by their correlation, every overlapping pair within its bounds (``[R; -R] x <= [u;
    -l]``), solved as ``stitch_assignment_lsqlin``, rounded and relative to tile 0."""
    n = int(n_tiles)
    mx = _table(max_xcorr_mat, 3, "max_xcorr_mat")
    ash = _table(absolute_shift_mat, 5, "absolute_shift_mat")
    mas = _table(max_allow_shift, 6, "max_allow_shift")
    ni, nj = _sorted_pairs(_overlap(overlap_matrix, n))
    filt = mx[:, 2] < xcorr_thresh
    mxf = mx[~filt]
    sw = np.sqrt(mxf[:, 2])
    C = _differences(_tiles(mxf[:, 0]), _tiles(mxf[:, 1]), n, sw)
    D = sw[:, None] * ash[~filt, 2:5]
    return _round_to_first(_bounded_system(n, ni, nj, C, D, mas))


def _grid_weights(tile_idx: np.ndarray, s, t, axis_weight, dense: bool) -> np.ndarray:
    """the weights of the pairs (s, t): neighbors along tile index x, y or z (one apart in all)
    get that axis's weight, other pairs 0. MATLAB's dense form switches on the axis that differs
    (an error if more than one does); the sparse form takes the first largest difference."""
    aw = np.asarray(axis_weight, dtype=np.float64).ravel()
    if aw.size != 3:
        raise AssignmentError("axis_weight must have 3 values (z, y, x)")
    w_xyz = np.array([aw[2], aw[1], aw[0]])
    tdiff = np.abs(tile_idx[s, 0:3] - tile_idx[t, 0:3])
    nb = _sum3(tdiff) == 1
    w = np.zeros(len(s))
    if not nb.any():
        return w
    if dense:
        nz = tdiff[nb] != 0
        if np.any(nz.sum(axis=1) != 1):
            raise AssignmentError(
                "SWITCH expression must be a scalar or a character vector.",
                "MATLAB:badSwitchExpression",
            )
        w[nb] = w_xyz[np.argmax(nz, axis=1)]
    else:
        w[nb] = w_xyz[np.argmax(tdiff[nb], axis=1)]
    return w


def grid_assignment(
    n_tiles: int,
    max_xcorr_mat,
    absolute_shift_mat,
    overlap_matrix,
    max_allow_shift,
    xcorr_thresh: float,
    tile_idx,
    axis_weight=AXIS_WEIGHT,
) -> np.ndarray:
    """``stitch_global_grid_assignment``: as ``global_assignment``, with the pairs at or above
    ``xcorr_thresh`` weighted by the axis along which they are neighbors in ``tile_idx`` (x, y,
    z columns; ``axis_weight`` (z, y, x)), 0 for other pairs. MATLAB's dense form up to 400 tiles
    and its sparse form above give the same shifts (both solve as ``stitch_assignment_lsqlin``);
    they differ only in the rule for non-integer tile indices (``_grid_weights``)."""
    n = int(n_tiles)
    mx = _table(max_xcorr_mat, 3, "max_xcorr_mat")
    ash = _table(absolute_shift_mat, 5, "absolute_shift_mat")
    mas = _table(max_allow_shift, 6, "max_allow_shift")
    tidx = _table(tile_idx, 4, "tile_idx")
    ni, nj = _sorted_pairs(_overlap(overlap_matrix, n))
    filt = mx[:, 2] < xcorr_thresh
    mxf = mx[~filt]
    s, t = _tiles(mxf[:, 0]), _tiles(mxf[:, 1])
    sw = np.sqrt(_grid_weights(tidx, s, t, axis_weight, n <= SPARSE_TILE_NUM))
    C = _differences(s, t, n, sw)
    D = sw[:, None] * ash[~filt, 2:5]
    return _round_to_first(_bounded_system(n, ni, nj, C, D, mas))


def grid_assignment_test(
    n_tiles: int,
    max_xcorr_mat,
    absolute_shift_mat,
    overlap_matrix,
    max_allow_shift,
    xcorr_thresh: float,
    tile_idx,
) -> np.ndarray:
    """``stitch_global_grid_assignment_test`` (shiftMethod 'test'), as MATLAB has it: the bounds
    rows in ``find(overlap_matrix)``'s column-major order, neighbors weighted 1, 10 and 0.1 along
    tile index x, y and z (all four columns compared), the pairs' values taken from rows of
    ``absolute_shift_mat`` at their pair index (``pair_index(i, j, n)``, an error past its end)
    and from its first three columns (MATLAB's 1-based tile numbers, then x), solved by plain
    ``lsqlin`` and rounded without the offset."""
    from petakit5d.matlabcompat import lsq
    from petakit5d.stitch.geometry import pair_index

    n = int(n_tiles)
    mx = _table(max_xcorr_mat, 3, "max_xcorr_mat")
    ash = _table(absolute_shift_mat, 5, "absolute_shift_mat")
    mas = _table(max_allow_shift, 6, "max_allow_shift")
    tidx = _table(tile_idx, 4, "tile_idx")
    ni, nj = _find(_overlap(overlap_matrix, n))
    mxf = mx[mx[:, 2] >= xcorr_thresh]
    s, t = _tiles(mxf[:, 0]), _tiles(mxf[:, 1])
    tdiff = np.abs(tidx[s] - tidx[t])
    nb = ((tdiff[:, 0] + tdiff[:, 1]) + tdiff[:, 2]) + tdiff[:, 3] == 1
    w = np.zeros(len(s))
    if nb.any():
        nz = tdiff[nb] != 0
        if np.any(nz.sum(axis=1) != 1):
            raise AssignmentError(
                "SWITCH expression must be a scalar or a character vector.",
                "MATLAB:badSwitchExpression",
            )
        w[nb] = np.array([1.0, 10.0, 0.1, 0.0])[np.argmax(nz, axis=1)]
    sw = np.sqrt(w)
    C = _differences(s, t, n, sw)
    c = np.array([pair_index(int(i), int(j), n) for i, j in zip(s, t, strict=True)], np.int64)
    if c.size and c.max() >= ash.shape[0]:
        raise AssignmentError(
            f"Index in position 1 exceeds array bounds. Index must not exceed {ash.shape[0]}.",
            "MATLAB:badsubscript",
        )
    rows = ash[c]
    vals = np.column_stack([rows[:, 0] + 1, rows[:, 1] + 1, rows[:, 2]])
    D = sw[:, None] * vals
    A = _bounds_matrix(ni, nj, n)
    B = _bound_rows(mas)
    x = np.zeros((n, 3))
    for k in range(3):
        x[:, k] = lsq.lsqlin(C, D[:, k], A, B[:, k])
    return _round_to_first(x, offset=False)


def group_assignment(
    n_tiles: int,
    max_xcorr_mat,
    absolute_shift_mat,
    overlap_matrix,
    max_allow_shift,
    xcorr_thresh: float,
    tile_idx,
    axis_weight,
    grp_idx,
    cuboid_overlap_ij_mat=None,
    xcorr_dir: str = "",
    parse_cluster: bool = False,
    node_factor: float = 1,
    mcc_mode: bool = False,
    config_file: str = "",
) -> np.ndarray:
    """``stitch_global_group_assignment``: step 1, the grid method in each group of tiles
    (``grp_idx``: each tile's group; with ``parse_cluster``, more than 100 tiles and more than one
    group, as tasks of ``global_grid_assignment`` on .mat files in
    ``<xcorr_dir>/xcorr_group_optimization/``); step 2, for each pair of groups (first group lower)
    with pairs at or above the threshold from the first to the second, two nodes at the centroids
    of the pairs' tiles, linked by the correlation-weighted mean shift of these pairs (weight 0.01)
    within their widest bounds, the nodes of a group at their centroids' distances (equalities),
    solved as ``stitch_assignment_lsqlin``; step 3, each group moved by its first node, groups
    without nodes at their stage positions relative to the tiles of the groups they overlap
    (``lsq.free_group_offsets``). Rounded and relative to tile 0; with one group, the grid
    method's shifts. ``cuboid_overlap_ij_mat`` is not used."""
    from petakit5d.core.dtypes import mround
    from petakit5d.matlabcompat import lsq
    from petakit5d.matlabcompat import reduce as mr

    del cuboid_overlap_ij_mat
    n = int(n_tiles)
    mx = _table(max_xcorr_mat, 3, "max_xcorr_mat")
    ash = _table(absolute_shift_mat, 5, "absolute_shift_mat")
    mas = _table(max_allow_shift, 6, "max_allow_shift")
    tidx = _table(tile_idx, 4, "tile_idx")
    om = _overlap(overlap_matrix, n).copy()
    grp = np.asarray(grp_idx, dtype=np.float64).ravel()
    uniq = np.unique(grp)
    ng = uniq.size
    parse_cluster = bool(parse_cluster) and n > 100 and ng > 1

    # step 1: the grid method in each group
    d_pre = np.zeros((n, 3))
    all_gp = np.zeros(mx.shape[0], bool)
    gind_cell = []
    inputs: list = [None] * ng
    outputs: list = [None] * ng
    tasks: list = [None] * ng
    if parse_cluster:
        from petakit5d.core.ids import get_uuid

        uuid = get_uuid()
        opt_dir = os.path.join(xcorr_dir, "xcorr_group_optimization")
        if not os.path.isdir(opt_dir):
            os.mkdir(opt_dir)
    for g in range(ng):
        ing = grp == uniq[g]
        gind = np.nonzero(ing)[0]
        gind_cell.append(gind)
        if gind.size == 1:
            continue
        # indices within the group
        new = np.cumsum(ing) - 1
        gp = np.isin(_tiles(mx[:, 0]), gind) & np.isin(_tiles(mx[:, 1]), gind)
        all_gp |= gp
        mx_g = mx[gp].copy()
        mx_g[:, 0:2] = new[_tiles(mx_g[:, 0:2])]
        ash_g = ash[gp].copy()
        ash_g[:, 0:2] = mx_g[:, 0:2]
        om_g = om[np.ix_(gind, gind)]
        if parse_cluster:
            inputs[g], outputs[g], tasks[g] = _group_task(
                opt_dir, g, gind.size, mx_g, ash_g, om_g, mas[gp], tidx[gind], xcorr_thresh,
                axis_weight, uuid,
            )  # fmt: skip
        else:
            d_pre[gind] = grid_assignment(
                gind.size, mx_g, ash_g, om_g, mas[gp], xcorr_thresh, tidx[gind], axis_weight
            )
    if parse_cluster:
        d_pre = _run_group_tasks(
            d_pre, grp, uniq, inputs, outputs, tasks, node_factor, mcc_mode, config_file
        )

    # with one group, the grid method
    if ng == 1:
        d = mround(d_pre)
        return d - d[0]

    # step 2: the groups, by their pairs at or above the threshold
    gn_i, gn_j = _tiles(mx[all_gp, 0]), _tiles(mx[all_gp, 1])
    om[gn_i, gn_j] = False
    keep = ~all_gp
    mx, ash, mas = mx[keep], ash[keep], mas[keep]
    filt = mx[:, 2] < xcorr_thresh
    mxf, ashf, masf = mx[~filt], ash[~filt], mas[~filt]
    first, second = _tiles(mxf[:, 0]), _tiles(mxf[:, 1])
    pair_sizes = np.zeros((ng, ng))
    dist = np.zeros((ng, ng, 3))
    cent = np.zeros((ng, ng, 6))
    sel_of = {}
    for gi in range(ng - 1):
        in_i = np.isin(first, gind_cell[gi])
        for gj in range(gi + 1, ng):
            sel = in_i & np.isin(second, gind_cell[gj])
            pair_sizes[gi, gj] = np.count_nonzero(sel)
            if not sel.any():
                continue
            sel_of[gi, gj] = sel
            w = mxf[sel, 2]
            num = mr.sum(ashf[sel, 2:5] * w[:, None], axis=0, matlab_layout=True)
            dist[gi, gj] = num / mr.sum(w)
            ci = np.unique(first[sel])
            cj = np.unique(second[sel])
            cent[gi, gj, 0:3] = mr.mean(d_pre[ci], axis=0, matlab_layout=True)
            cent[gi, gj, 3:6] = mr.mean(d_pre[cj], axis=0, matlab_layout=True)
    # the pairs of groups in find's (column-major) order, two nodes each
    gpj, gpi = np.nonzero(pair_sizes.T)
    n_p = gpi.size
    n_node = 2 * n_p
    first_nodes = np.arange(n_p)
    second_nodes = n_p + np.arange(n_p)
    sw = math.sqrt(0.01)
    C = _differences(first_nodes, second_nodes, n_node, np.full(n_p, sw))
    D = sw * dist[gpi, gpj]
    gmas = np.zeros((n_p, 6))
    for k in range(n_p):
        sel = sel_of[int(gpi[k]), int(gpj[k])]
        gmas[k, 0:3] = mr.min(masf[sel, 0:3], axis=0, matlab_layout=True)[0]
        gmas[k, 3:6] = mr.max(masf[sel, 3:6], axis=0, matlab_layout=True)[0]
    # equalities: the nodes of a group at their centroids' distances
    eq_rows, eq_vals = [], []
    node_cent = np.zeros((n_node, 3))
    nodes_of = []
    for g in range(ng):
        nodes = np.unique(np.r_[first_nodes[gpi == g], second_nodes[gpj == g]])
        nodes_of.append(nodes)
        if nodes.size == 0:
            continue
        cg = np.array(
            [cent[gpi[k], gpj[k], 0:3] if k < n_p else cent[gpi[k - n_p], gpj[k - n_p], 3:6]
             for k in nodes]
        )  # fmt: skip
        node_cent[nodes] = cg
        # find(triu(ones(m), 1)): column-major
        for c in range(1, nodes.size):
            for r in range(c):
                eq_rows.append((nodes[r], nodes[c]))
                eq_vals.append(cg[c] - cg[r])
    if eq_rows:
        e = np.array(eq_rows, dtype=np.int64)
        Aeq = _differences(e[:, 0], e[:, 1], n_node)
        BEQ = np.array(eq_vals)
    else:
        Aeq, BEQ = None, None
    A = _bounds_matrix(first_nodes, second_nodes, n_node)
    d_opt = lsq.assignment(C, D, A, _bound_rows(gmas), Aeq, BEQ)

    # step 3: each group moved by its first node; groups without nodes at their stage positions
    # relative to the tiles of the groups they overlap
    d = np.zeros((n, 3))
    gid = np.zeros(n, np.int64)
    placed = np.zeros(ng, bool)
    for g in range(ng):
        gind = gind_cell[g]
        gid[gind] = g
        if nodes_of[g].size == 0:
            continue
        placed[g] = True
        k = nodes_of[g][0]
        d[gind] = (d_pre[gind] - node_cent[k]) + d_opt[k]
    if not placed.all():
        free = ~placed[gid]
        offset = lsq.free_group_offsets(placed, gid, d, d_pre, om)
        d[free] = d_pre[free] + offset[gid[free]]
    return _round_to_first(d)


def _group_task(opt_dir, g, n_g, mx_g, ash_g, om_g, mas_g, tidx_g, xcorr_thresh, axis_weight, uuid):
    """a group's optimization file (MATLAB's form) and its task"""
    from petakit5d.core.text import through_text
    from petakit5d.io.matfile import write_mat
    from petakit5d.jobs.tasks import Task

    path = os.path.join(opt_dir, f"xcorr_group_{g + 1}_opt_information.mat")
    tmp = os.path.join(opt_dir, f"xcorr_group_{g + 1}_opt_information_{uuid}.mat")
    one = np.column_stack([mx_g[:, 0:2] + 1, mx_g[:, 2]])
    write_mat(
        tmp,
        {
            "nF_g": float(n_g),
            "max_xcorr_mat_g": one,
            "absolute_shift_mat_g": np.column_stack([one[:, 0:2], ash_g[:, 2:5]]),
            "overlap_matrix_g": np.asarray(om_g, bool),
            "max_allow_shift_g": mas_g,
            "tileIdx_g": tidx_g,
        },
    )
    os.replace(tmp, path)
    out = os.path.join(opt_dir, f"xcorr_group_{g + 1}_opt_result.mat")
    # MATLAB passes the threshold as %d (%e for a non-integer: 7 digits) and the weights as %.20d
    task = Task.of(
        global_grid_assignment,
        path,
        out,
        float(through_text(float(xcorr_thresh), "%d")),
        [float(v) for v in np.asarray(axis_weight, dtype=np.float64).ravel()],
        uuid,
    )
    return path, out, task


def _run_group_tasks(d_pre, grp, uniq, inputs, outputs, tasks, node_factor, mcc_mode, config):
    """step 1's tasks of stitch_global_group_assignment on the cluster"""
    from petakit5d.core.dtypes import mround
    from petakit5d.io.matfile import read_mat
    from petakit5d.jobs import framework

    ng = uniq.size
    if any(t is None for t in tasks):
        # a group of one tile has no task: MATLAB's call gets empty paths and errors
        raise AssignmentError(
            "A group of one tile has no optimization task (MATLAB's cluster path errors for it)."
        )
    _, counts = np.unique(grp, return_counts=True)
    mem = max(20.0, float(counts.max()) ** 2 * 0.0001) * node_factor
    min_task_job_num = 1
    if ng >= 4:
        min_task_job_num = max(2, int(mround(ng / 8)))
    done: list[bool] = []
    for i in range(3):
        done = framework.run_tasks(
            inputs,
            outputs,
            tasks,
            max_trial_num=2,
            parse_cluster=True,
            mem_allocate=mem * 2**i,
            min_task_job_num=min_task_job_num,
            mcc_mode=mcc_mode,
            config_file=config,
        )
        if all(done):
            break
    if not all(done):
        raise AssignmentError("group registartion cannot be fully done!")
    for g in range(ng):
        a = read_mat(outputs[g], "d_shift_g")
        d_pre[grp == uniq[g]] = np.asarray(a["d_shift_g"], dtype=np.float64).reshape(-1, 3)
    return d_pre


@register(
    "stitch_global_grid_assignment_wrapper",
    Param("inputFullpath", "@(x) ischar(x)", "The group's optimization information (.mat)."),
    Param("outputFullpath", "@(x) ischar(x)", "The result file (.mat, with d_shift_g)."),
    Param("xcorr_thresh", "@(x) isscalar(x)", "Pairs below this correlation are not fitted."),
    Param(
        "axisWeight",
        "@(x) isvector(x)",
        "Weights of the neighbors along tile index z, y and x.",
        convert="size_yxz",
    ),
    Param("uuid", "@(x) ischar(x)", "Id of the temporary result file."),
    task="S05",
)
def global_grid_assignment(
    input_fullpath: str,
    output_fullpath: str,
    xcorr_thresh: float,
    axis_weight: Sequence[float],
    uuid: str,
) -> None:
    """The grid assignment of a group of tiles as a task: reads the group's problem from a .mat
    file (MATLAB's form: ``nF_g``, ``max_xcorr_mat_g``, ``absolute_shift_mat_g``,
    ``overlap_matrix_g``, ``max_allow_shift_g``, ``tileIdx_g``) and writes ``d_shift_g`` to
    ``output_fullpath`` (through a temporary file); an existing result is kept."""
    from petakit5d.io.matfile import read_mat, write_mat

    if os.path.exists(output_fullpath):
        print(f"output {output_fullpath} already exists, skip it!")
        return
    names = (
        "nF_g",
        "max_xcorr_mat_g",
        "absolute_shift_mat_g",
        "overlap_matrix_g",
        "max_allow_shift_g",
        "tileIdx_g",
    )
    a = read_mat(input_fullpath, list(names))
    n = int(np.asarray(a["nF_g"]).ravel()[0])
    mx = _table(a["max_xcorr_mat_g"], 3, "max_xcorr_mat_g")
    ash = _table(a["absolute_shift_mat_g"], 5, "absolute_shift_mat_g")
    mx[:, 0:2] -= 1
    ash[:, 0:2] -= 1
    om = np.asarray(a["overlap_matrix_g"], dtype=bool).reshape(n, n)
    d = grid_assignment(
        n,
        mx,
        ash,
        om,
        a["max_allow_shift_g"],
        float(np.asarray(xcorr_thresh).ravel()[0]),
        a["tileIdx_g"],
        axis_weight,
    )
    tmp = f"{output_fullpath[:-4]}_{uuid}.mat"
    write_mat(tmp, {"d_shift_g": d})
    os.replace(tmp, output_fullpath)


# ------------------------------------------------------------------------------------------------
# stitch_shift_assignment
# ------------------------------------------------------------------------------------------------


def xcorr_result_paths(xcorr_dir: str, n_tiles: int) -> list[str]:
    """The registration result of each leading tile: ``<xcorr_dir>/xcorr_tile_<k>.mat`` (k from 1),
    above 10,000 tiles in subfolders ``xcorr_tiles_<first>_<last>/`` of 10,000 (made here)."""
    n = int(n_tiles)
    if n <= XCORR_BATCH_SIZE:
        return [os.path.join(xcorr_dir, f"xcorr_tile_{k}.mat") for k in range(1, n + 1)]
    out = []
    for b in range(math.ceil(n / XCORR_BATCH_SIZE)):
        s = b * XCORR_BATCH_SIZE + 1
        t = min((b + 1) * XCORR_BATCH_SIZE, n)
        sub = os.path.join(xcorr_dir, f"xcorr_tiles_{s:06d}_{t:06d}")
        if not os.path.isdir(sub):
            os.mkdir(sub)
        out.extend(os.path.join(sub, f"xcorr_tile_{k}.mat") for k in range(s, t + 1))
    return out


def _matlab_size(shape) -> np.ndarray:
    """MATLAB's getImageSize of an image of ``(z, y, x)`` shape: ``[ny nx nz]``"""
    s = [int(v) for v in shape]
    if len(s) >= 3:
        return np.array([s[1], s[2], s[0]] + s[3:], dtype=np.float64)
    return np.array(s, dtype=np.float64)


def _tile_mip_slabs(tile_paths, pool_size, parse_cluster, mcc_mode, config_file) -> str:
    """``compute_tile_mip_slabs`` (S04's): the MIP slabs of large tiles for their registration;
    the name of their folder"""
    try:
        from petakit5d.stitch.registration import compute_tile_mip_slabs
    except ImportError:
        compute_tile_mip_slabs = _compute_tile_mip_slabs
    return compute_tile_mip_slabs(
        tile_paths, pool_size, parse_cluster=parse_cluster, mcc_mode=mcc_mode,
        config_file=config_file,
    )  # fmt: skip


def _compute_tile_mip_slabs(tile_paths, pool_size, *, parse_cluster, mcc_mode, config_file):
    """stand-in for S04's ``compute_tile_mip_slabs`` (MATLAB's, on T02's ``mip_zarr``): slabs of
    every tile in ``<tile folder>/MIP_Slabs_pooling_<pool sizes, y x z>/``, existing ones kept"""
    from petakit5d.core.axes import to_yxz
    from petakit5d.core.dtypes import mround
    from petakit5d.core.paths import fileparts
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.tools.mip import mip_zarr

    ps = [int(v) for v in np.asarray(pool_size).ravel()]
    yxz = [int(v) for g in range(0, len(ps), 3) for v in to_yxz(ps[g : g + 3])]
    mip_dir = "MIP_Slabs_pooling_" + "_".join(str(v) for v in yxz)
    parts = [fileparts(p) for p in tile_paths]
    for d in sorted({p[0] for p in parts}):
        path = f"{d}/{mip_dir}"
        if not os.path.isdir(path):
            os.mkdir(path)
    tasks = [
        Task.of(
            mip_zarr, p, result_dir_name=mip_dir, axis=(1, 1, 1), batch_size=(1024, 1024, 1024),
            pool_size=tuple(ps), mip_slab=True, parse_cluster=parse_cluster, parse_parfor=False,
            mcc_mode=mcc_mode, config_file=config_file,
        )
        for p in tile_paths
    ]  # fmt: skip
    outputs = [f"{d}/{mip_dir}/{name}_z.zarr" for d, name, _ in parts]
    mem = 1024**3 * 4 / 2**30 * 3
    for i in range(3):
        done = framework.run_tasks(
            list(tile_paths),
            outputs,
            tasks,
            cpus_per_task=1,
            max_trial_num=2,
            parse_cluster=parse_cluster,
            mem_allocate=mem * 2**i,
            min_task_job_num=int(mround(len(tile_paths) / 2)),
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        if all(done):
            break
    return mip_dir


def _read_xcorr(path: str):
    """a registration result: 0-based pairs (k, 2), relative shifts (k, 3) x, y, z, max xcorr"""
    from petakit5d.io.matfile import read_mat

    a = read_mat(path, ["relative_shift_mat", "max_xcorr_mat", "pair_indices"])
    pairs = np.asarray(a["pair_indices"], dtype=np.float64).reshape(-1, 2)
    rel = np.asarray(a["relative_shift_mat"], dtype=np.float64).reshape(-1, 3)
    mxv = np.asarray(a["max_xcorr_mat"], dtype=np.float64).reshape(-1)
    return _tiles(pairs - 1), rel, mxv


def _bounds(lo, hi, order) -> np.ndarray:
    """``[l .* o + u .* (o - 1), u .* o + l .* (o - 1)]``: a pair whose second tile is not ahead
    along an axis gets the bounds turned around"""
    o = order.astype(np.float64)
    return np.hstack([lo * o + hi * (o - 1), hi * o + lo * (o - 1)])


def shift_assignment(
    tile_paths: Sequence[str],
    xcorr_dir: str,
    im_sizes,
    xyz,
    xyz_voxel_sizes: Sequence[float],
    data_order: str,
    overlap_matrix,
    overlap_regions,
    max_offset: Sequence[float],
    xcorr_downsample: Sequence[float],
    xcorr_thresh: float,
    tile_idx,
    assign_method: str = "grid",
    stitch_2d: bool = False,
    axis_weight: Sequence[float] | None = AXIS_WEIGHT,
    group_file: str = "",
    large_file: bool = False,
    pool_size: Sequence[float] | None = None,
    parse_cluster: bool = False,
    node_factor: float = 2,
    mcc_mode: bool = False,
    config_file: str = "",
    local_worker_num: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """``[xyz_shift, dxyz_shift] = stitch_shift_assignment(...)``: registers the overlapping
    tiles and assigns their shifts by ``assign_method`` (``METHODS``).

    The pairs of ``overlap_matrix`` (pairs that touch only at a corner left out for 'grid' and
    'test', and within each group for 'group'; groups from ``group_file``'s x, y, z[, did], g
    columns, else the tiles' locations) are registered by one task per leading tile (per pair with
    ``large_file``, on the tiles' MIP slabs) that writes ``xcorr_tile_<k>.mat`` in ``xcorr_dir``
    (existing results are kept), tried three times with 2, 4 and 8 CPUs. The pairs' shifts get
    the sign of their tiles' order in ``xyz`` along each axis, their bounds the overlap size less
    1 (less 3 for 'grid'; z 0 in 2D for 'grid' and 'group') within ``max_offset``, turned around
    where the second tile is not ahead, twice ``max_offset`` between groups. Returns the shifted
    positions ``xyz + dxyz_shift`` and ``dxyz_shift = d_shift * xyz_voxel_sizes`` (x, y, z um)."""
    from petakit5d.core.paths import mkdir_recursive
    from petakit5d.core.text import sprintf, through_mat2str
    from petakit5d.io import get_image_size
    from petakit5d.jobs import framework
    from petakit5d.matlabcompat.stats import prctile
    from petakit5d.stitch.geometry import pair_index, parse_data_order

    print("Compute cross-correlation based registration between overlap tiles...")
    t0 = time.perf_counter()
    if axis_weight is None or np.size(axis_weight) == 0:
        axis_weight = AXIS_WEIGHT
    if not os.path.isdir(xcorr_dir):
        mkdir_recursive(xcorr_dir, True)
    vs = np.asarray(xyz_voxel_sizes, dtype=np.float64).reshape(3)
    mo_xyz = np.asarray(max_offset, dtype=np.float64).reshape(3)[::-1].copy()
    ds = np.asarray(xcorr_downsample, dtype=np.float64).ravel()
    n = len(tile_paths)
    xyz = _table(xyz, 3, "xyz")
    sizes = _table(im_sizes, 3, "im_sizes")
    tidx = _table(tile_idx, 4, "tile_idx")
    regions = _table(overlap_regions, 6, "overlap_regions")
    om = _overlap(overlap_matrix, n).copy()
    method = assign_method

    # the pairs of each method
    ti, tj = _find(om)
    grp_idx = None
    if method in ("grid", "test"):
        if not np.all(tidx[:, 3] == tidx[0, 3]):
            _warn(
                "The tiles are not from the same data folder, please make sure they are "
                "continuous grids when using grid method."
            )
        corner = _sum3(np.abs(tidx[ti, 0:3] - tidx[tj, 0:3])) != 1
        om[ti[corner], tj[corner]] = False
    elif method == "group":
        grp_idx = _groups(tidx, group_file)
        for g in np.unique(grp_idx):
            ing = grp_idx == g
            ij = ing[ti] & ing[tj]
            ti_g, tj_g = ti[ij], tj[ij]
            td = np.abs(tidx[ti_g] - tidx[tj_g])
            corner = ((td[:, 0] + td[:, 1]) + td[:, 2]) + td[:, 3] != 1
            om[ti_g[corner], tj_g[corner]] = False

    mip_dir = ""
    if large_file:
        mip_dir = _tile_mip_slabs(tile_paths, pool_size, parse_cluster, mcc_mode, config_file)

    # the registration tasks
    out_paths = xcorr_result_paths(xcorr_dir, n)
    ti, tj = _find(om)
    dop = list(parse_data_order(data_order))
    cuboids = np.hstack((xyz, xyz + (sizes[:, dop] - 1) * vs))
    pinds = np.array([pair_index(int(a), int(b), n) for a, b in zip(ti, tj, strict=True)], int)
    cov = regions[pinds] if pinds.size else np.zeros((0, 6))
    # the cuboids as MATLAB's mat2str gives them to the tasks (15 digits)
    cub15 = through_mat2str(cuboids).tolist()
    reg15 = dict(zip(pinds.tolist(), through_mat2str(cov).tolist(), strict=True))
    common = _registration_options(
        vs, data_order, stitch_2d, ds, max_offset, large_file, mip_dir, pool_size,
        parse_cluster, mcc_mode, config_file,
    )  # fmt: skip
    inputs, outputs, tasks = [], [], []
    if large_file:
        for f in range(ti.size):
            a, b = int(ti[f]), int(tj[f])
            inputs.append(tile_paths[a])
            outputs.append(os.path.join(xcorr_dir, f"xcorr_tile_{a + 1}_{b + 1}.mat"))
            tasks.append(_registration_task(tile_paths, a, [b], outputs[-1], cub15, reg15, n,
                                            common))  # fmt: skip
        min_task_job_num = ti.size
    else:
        for f in range(n):
            tj_f = tj[ti == f]
            if tj_f.size == 0:
                continue
            inputs.append(tile_paths[f])
            outputs.append(out_paths[f])
            tasks.append(_registration_task(tile_paths, f, list(tj_f), out_paths[f], cub15,
                                            reg15, n, common))  # fmt: skip
        min_task_job_num = 1

    print("Compute pairwise cross correlation between overlap tiles...")
    sz = _matlab_size(get_image_size(tile_paths[0]))
    ext = (cov[:, 3:6] - cov[:, 0:3]) / vs + mo_xyz
    if sz.size != 3:
        raise AssignmentError(
            "Arrays have incompatible sizes for this operation.", "MATLAB:sizeDimensionsMustMatch"
        )
    m = np.fmin(ext, sz)
    raw = ((m[:, 0] * m[:, 1]) * m[:, 2]) * 4 / 1024**3
    ds_prod = (ds[1] * ds[2]) * ds[0] if ds.size == 3 else float(np.prod(ds))
    mem = float(prctile(np.ceil(raw) * (2 + 30 / ds_prod), 99)[0]) * node_factor
    done: list[bool] = []
    for i in range(3):
        done = framework.run_tasks(
            inputs,
            outputs,
            tasks,
            cpus_per_task=2 * 2**i,
            max_trial_num=2,
            parse_cluster=parse_cluster,
            mem_allocate=mem * 2**i if tasks else None,
            min_task_job_num=min_task_job_num,
            mcc_mode=mcc_mode,
            config_file=config_file,
            local_worker_num=local_worker_num,
        )
        if all(done):
            break

    # the results, by first, then second tile
    if not all(done):
        raise AssignmentError("Some xcorr files are missing!")
    n_all = n * (n - 1) // 2
    ash = np.zeros((n_all, 5))
    mx = np.zeros((n_all, 3))
    filled = np.zeros(n_all, bool)
    for path in outputs:
        pr, rel, mxv = _read_xcorr(path)
        a, b = pr[:, 0], pr[:, 1]
        p = np.array([pair_index(int(u), int(v), n) for u, v in zip(a, b, strict=True)], int)
        order_flag = (0.5 - (xyz[a] > xyz[b])) * 2
        ash[p] = np.column_stack([a, b, rel * order_flag])
        mx[p] = np.column_stack([a, b, mxv])
        filled[p] = True
    ash, mx = ash[filled], mx[filled]

    print(f"Optimize registration for all tiles with {method} method ...")
    print(
        sprintf(
            "xcorr threshold: %f , %d / %d pairs are below the threashold.",
            xcorr_thresh,
            int(np.count_nonzero(mx[:, 2] < xcorr_thresh)),
            mx.shape[0],
        )
    )
    # the overlap sizes (voxels) and bounds of the pairs in the order of max_xcorr_mat
    i_p, j_p = _tiles(mx[:, 0]), _tiles(mx[:, 1])
    pinds_p = np.array([pair_index(int(u), int(v), n) for u, v in zip(i_p, j_p, strict=True)], int)
    ol = regions[pinds_p] if pinds_p.size else np.zeros((0, 6))
    ol_size = (ol[:, 3:6] - ol[:, 0:3]) / vs
    neq = mx.shape[0]
    lo = -np.ones((neq, 1)) * mo_xyz
    order = xyz[j_p] >= xyz[i_p]
    if method == "local":
        d_shift = local_assignment(n, mx, ash, om, xcorr_thresh)
    elif method == "global":
        hi = np.fmin(ol_size - 1, mo_xyz)
        d_shift = global_assignment(n, mx, ash, om, _bounds(lo, hi, order), xcorr_thresh)
    elif method == "grid":
        hi = np.fmin(ol_size - 3, mo_xyz)
        if stitch_2d:
            hi[:, 2] = 0
        mas = _bounds(lo, hi, order)
        d_shift = grid_assignment(n, mx, ash, om, mas, xcorr_thresh, tidx, axis_weight)
    elif method == "group":
        hi = np.fmin(ol_size - 1, mo_xyz)
        if stitch_2d:
            hi[:, 2] = 0
        mas = _bounds(lo, hi, order)
        inter = grp_idx[i_p] != grp_idx[j_p]
        k = int(np.count_nonzero(inter))
        mas[inter] = np.hstack([-np.ones((k, 1)) * mo_xyz, np.ones((k, 1)) * mo_xyz]) * 2
        d_shift = group_assignment(
            n, mx, ash, om, mas, xcorr_thresh, tidx, axis_weight, grp_idx, cov, xcorr_dir,
            parse_cluster, node_factor, mcc_mode, config_file,
        )  # fmt: skip
    elif method == "test":
        hi = np.fmin((cov[:, 3:6] - cov[:, 0:3]) / vs - 1, mo_xyz)
        if hi.shape[0] != neq:
            raise AssignmentError(
                "Dimensions of arrays being concatenated are not consistent.",
                "MATLAB:catenate:dimensionMismatch",
            )
        mas = np.hstack([lo, hi])
        d_shift = grid_assignment_test(n, mx, ash, om, mas, xcorr_thresh, tidx)
    else:
        raise AssignmentError(
            "Unrecognized function or variable 'd_shift'.", "MATLAB:UndefinedFunction"
        )

    dxyz_shift = d_shift * vs
    xyz_shift = xyz + dxyz_shift
    print("xcorr registration is done!")
    print(f"Elapsed time is {time.perf_counter() - t0:.6f} seconds.")
    print("")
    return xyz_shift, dxyz_shift


def _groups(tidx: np.ndarray, group_file: str) -> np.ndarray:
    """each tile's group: of the nearest row of ``group_file`` (x, y, z[, did], g; 4 columns: did
    0), else the tile's location"""
    from petakit5d.matlabcompat.neighbors import knnsearch

    if not group_file:
        if np.all(tidx[:, 3] == tidx[0, 3]):
            _warn(
                "The input group file is not provided, and the tiles are from the same data "
                "folder, treat all tiles as one group."
            )
        g_idx, gid = tidx, tidx[:, 3]
    else:
        from petakit5d.stitch.imagelist import read_table

        t = read_table(group_file)
        col = lambda k: np.asarray(t[k], dtype=np.float64)  # noqa: E731
        did = np.zeros(t.height) if len(t.names) == 4 else col("did")
        g_idx = np.column_stack([col("x"), col("y"), col("z"), did])
        gid = col("g")
    kidx, d = knnsearch(g_idx, tidx, 1)
    if np.any(d != 0):
        _warn("Some tiles are not included in the groups")
    return gid[kidx[:, 0]]


def _registration_options(vs, data_order, stitch_2d, ds, max_offset, large_file, mip_dir,
                          pool_size, parse_cluster, mcc_mode, config_file) -> dict:  # fmt: skip
    """the options of the registration tasks, with the rounding of MATLAB's call text
    (``mat2str``: 15 digits for the offsets and pool sizes; ``%.20d`` keeps the rest)"""
    from petakit5d.core.text import through_mat2str

    mo = np.asarray(max_offset, dtype=np.float64).ravel()
    ps = None if pool_size is None or np.size(pool_size) == 0 else np.asarray(pool_size, float)
    return {
        "xyz_voxel_sizes": [float(v) for v in vs],
        "data_order": data_order,
        "stitch_2d": bool(stitch_2d),
        "down_sample": [float(v) for v in ds],
        "max_offset": [float(v) for v in np.atleast_1d(through_mat2str(mo))],
        "large_file": bool(large_file),
        "mip_dir_str": mip_dir,
        "pool_size": None if ps is None else [float(v) for v in through_mat2str(ps.ravel())],
        "parse_cluster": bool(parse_cluster),
        "mcc_mode": bool(mcc_mode),
        "config_file": config_file,
    }


def _registration_task(tile_paths, f, js, out_path, cub15, reg15, n, common):
    """the registration of tile f with tiles js (``cross_correlation_registration_wrapper``):
    0-based pairs, the cuboids (``cub15``, per tile) and overlap cuboids (``reg15``, per pair
    index) as MATLAB's ``mat2str`` gives them"""
    from petakit5d.jobs.tasks import Task
    from petakit5d.stitch.geometry import pair_index

    js = [int(j) for j in js]
    kw = dict(common)
    vs, dorder = kw.pop("xyz_voxel_sizes"), kw.pop("data_order")
    return Task(
        REGISTRATION_TASK,
        [
            tile_paths[f],
            [tile_paths[j] for j in js],
            out_path,
            [[f, j] for j in js],
            cub15[f],
            [cub15[j] for j in js],
            [reg15[pair_index(f, j, n)] for j in js],
            vs,
            dorder,
        ],
        kw,
    )
