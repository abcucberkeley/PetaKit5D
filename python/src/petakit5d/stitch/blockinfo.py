"""Stitching's block information (``stitch_process_block_info``, ``stitch_organize_block_info``):
for each batch of the stitched volume, the tiles it takes, where each tile's part lies in the
batch, in the stitched volume and in the tile, and the parts of that tile which other tiles
fill (its half-overlap regions), for the blending tasks (``stitch.block``).

The files keep MATLAB's names and values (boxes ``[ymin xmin zmin ymax xmax zmax]``, 1-based,
inclusive; sizes ``[ny nx nz]``; tiles and batches numbered from 1), so MATLAB and the port
read each other's:

- ``<stitchInfo>/<name>_block_info_task_size_<k>.mat`` (v7.3): ``tileFns``, ``int_xyz_shift``
  (columns in the tiles' axis order), ``imSizes``, ``overlap_matrix``, ``half_ol_region_cell``,
  ``ol_region_cell``, ``overlap_map_mat``, ``bboxes`` (each tile's box in the stitched volume,
  zeros for tiles outside it), ``bboxStart_mat``, ``bboxEnd_mat`` (the same box in the tile),
  ``nvSize``, ``blockSize`` (the batch size), ``bSubSz``, ``BorderSize``;
- ``<stitchInfo>/<name>_task_size_<k>/stitch_block_info_blocks_<s>_<t>.json``: the batches
  ``s`` to ``t`` of one blending task, as MATLAB's ``jsonencode(..., 'PrettyPrint', true)``
  writes them, in subfolders ``task_size_<k>_tasks_<a>_<b>/`` when there are more than
  ``max_file_num_per_folder`` tasks;
- ``<...>_task_size_<k>//block_flags/stitch_block_info_task_size_<k>_blocks_<s>_<t>.mat``: the
  flag of each organizing task (``blockInds``, ``t1``).

The port's forms (``TileBlock``, ``RegionBlock``): tiles numbered from 0; boxes ``(z0, y0, x0,
z1, y1, x1)``, 0-based, end-exclusive; batches numbered from 0 in MATLAB's order (y fastest,
then x, then z); sizes and the border ``(z, y, x)``; a tile's half-overlap regions are boxes
over its voxels (``stitch.geometry.half_overlap_regions``). Tile offsets ``int_xyz_shift`` keep
their x, y, z columns (exception 2 of the conventions). Paths are built as MATLAB builds them
(``'%s/...'`` after folders that end with ``/``: the per-block names have ``//``), so that the
stitch information files hold the same strings.

Quirks kept: a tile takes a batch when its box and the batch's overlap by MATLAB's test (each
axis: the batch's start or end in the tile, or the tile's end in the batch; that is any
overlap); the organizing tasks' flags are made after their JSON files; with ``max_file_num_per
_folder`` exceeded the organizing tasks hold ``max_file_num_per_folder * task_size`` blocks or
a power-of-two fraction of it.
"""

from __future__ import annotations

import json
import math
import os
import re
import time
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "BlockInfoError",
    "RegionBlock",
    "TileBlock",
    "block_info_bytes",
    "block_info_json",
    "block_task_size",
    "organize_block_info",
    "organize_blocks",
    "per_block_info_paths",
    "process_block_info",
    "read_block_info",
    "read_per_block_info",
    "tile_blocks_from_matlab",
    "tile_blocks_to_matlab",
    "write_block_info",
]

Box = tuple[int, int, int, int, int, int]


class BlockInfoError(RuntimeError):
    """An error MATLAB raises; ``identifier`` is MATLAB's error id (``''`` for none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


@dataclass(frozen=True)
class RegionBlock:
    """One half-overlap region of a tile in a batch (MATLAB's ``mblockInfo(k)``): the region's
    part of the batch in the batch, in the stitched volume and in the tile (None where MATLAB's
    is empty)."""

    batch_box: Box | None
    canvas_box: Box | None
    tile_box: Box | None


@dataclass(frozen=True)
class TileBlock:
    """One tile's part of a batch (MATLAB's ``stitchBlockInfo{b}(i)``): its box in the batch
    (``bCoords``; with a border, in the batch with its border), in the stitched volume
    (``wCoords``) and in the tile (``bboxCoords``), the half-overlap regions of the tile in
    the batch (``mblockInfo``) and the border (``(z, y, x)``; None for MATLAB's ``[]``)."""

    tile: int
    batch_box: Box | None
    canvas_box: Box | None
    tile_box: Box | None
    regions: tuple[RegionBlock, ...]
    border_size: tuple[int, int, int] | None


# ---------------------------------------------------------------------------------------------
# conversions to and from MATLAB's values
# ---------------------------------------------------------------------------------------------


def _row_box(v: Any) -> Box | None:
    """a MATLAB box row (as jsondecode or load gives it) as a box; None for an empty one"""
    from petakit5d.core.bbox import box_from_matlab

    a = np.asarray(v if v is not None else [], dtype=np.float64).reshape(-1)
    if a.size == 0:
        return None
    if a.size != 6:
        raise BlockInfoError(
            f"Unable to perform assignment because the size of the left side is 1-by-6 and the "
            f"size of the right side is 1-by-{a.size}.",
            "MATLAB:subsassigndimmismatch",
        )
    return box_from_matlab(a)  # type: ignore[return-value]


def _box_row(box: Box | None) -> np.ndarray:
    """a box as MATLAB's 1 x 6 double row; 0 x 6 for None"""
    from petakit5d.core.bbox import box_to_matlab

    if box is None:
        return np.zeros((0, 6))
    return np.asarray(box_to_matlab(box), dtype=np.float64)


def _yxz_row(v: Sequence[int] | None) -> np.ndarray | None:
    from petakit5d.core.axes import to_yxz

    if v is None:
        return None
    return np.asarray(to_yxz(tuple(int(x) for x in v)), dtype=np.float64)


def _zyx(v: Any) -> tuple[int, int, int] | None:
    """MATLAB's [y x z] (a row or column of 3 integers) as (z, y, x); None for []"""
    from petakit5d.core.axes import from_yxz

    a = np.asarray(v if v is not None else [], dtype=np.float64).reshape(-1)
    if a.size == 0:
        return None
    return tuple(int(x) for x in from_yxz(tuple(a)))  # type: ignore[return-value]


def _elements(v: Any) -> list[dict]:
    """the structs of a decoded MATLAB value: [] (none), a struct (one), a struct array or a
    list of structs (jsondecode's struct array, a cell of structs)"""
    from petakit5d.io.matfile import StructArray

    if v is None:
        return []
    if isinstance(v, dict):
        return [v]
    if isinstance(v, StructArray):
        return v.ravel()
    if isinstance(v, np.ndarray):
        if v.dtype == object:
            return [e for x in v.ravel(order="F") for e in _elements(x)]
        if v.size == 0:
            return []
        raise BlockInfoError(f"expected structs in the block information, got {v.dtype}")
    if isinstance(v, list | tuple):
        out: list[dict] = []
        for x in v:
            out.extend(_elements(x))
        return out
    raise BlockInfoError(f"expected structs in the block information, got {type(v).__name__}")


def _tile_from_matlab(d: dict) -> TileBlock:
    regions = tuple(
        RegionBlock(
            _row_box(m.get("bCoords")), _row_box(m.get("wCoords")), _row_box(m.get("bboxCoords"))
        )
        for m in _elements(d.get("mblockInfo"))
    )
    return TileBlock(
        tile=int(float(np.asarray(d["tileInd"]).reshape(-1)[0])) - 1,
        batch_box=_row_box(d.get("bCoords")),
        canvas_box=_row_box(d.get("wCoords")),
        tile_box=_row_box(d.get("bboxCoords")),
        regions=regions,
        border_size=_zyx(d.get("borderSize")),
    )


def tile_blocks_from_matlab(value: Any) -> list[list[TileBlock]]:
    """MATLAB's per-block information (decoded JSON: per batch ``[]``, a struct or a list of
    structs; a struct array of batches x tiles, as jsondecode makes of equal tile counts, or a
    cell, as ``load`` gives ``stitchBlockInfo_t``) as a list of each batch's ``TileBlock``s."""
    from petakit5d.io.matfile import StructArray

    if isinstance(value, StructArray) and len(value.shape) == 2:
        # (MATLAB: arrayfun(@(x) stitchBlockInfo(x, :), ...): a batch per row)
        n, m = value.shape
        return [[_tile_from_matlab(value[r, c]) for c in range(m)] for r in range(n)]
    if isinstance(value, dict):
        return [[_tile_from_matlab(value)]]
    if isinstance(value, np.ndarray) and value.dtype == object:
        value = list(value.ravel(order="F"))
    if not isinstance(value, list | tuple):
        raise BlockInfoError(f"expected the batches' block information, got {type(value)}")
    if value and all(isinstance(x, dict) for x in value):
        # (jsondecode of one tile per batch: a struct array, one batch per element)
        return [[_tile_from_matlab(x)] for x in value]
    return [[_tile_from_matlab(d) for d in _elements(b)] for b in value]


def _tile_to_matlab(t: TileBlock) -> dict:
    regions = [
        {
            "bCoords": _box_row(r.batch_box),
            "wCoords": _box_row(r.canvas_box),
            "bboxCoords": _box_row(r.tile_box),
        }
        for r in t.regions
    ]
    return {
        "tileInd": float(t.tile + 1),
        "bCoords": _box_row(t.batch_box),
        "wCoords": _box_row(t.canvas_box),
        "bboxCoords": _box_row(t.tile_box),
        # struct([]), a struct or a 1 x K struct array
        "mblockInfo": None if not regions else (regions[0] if len(regions) == 1 else regions),
        "borderSize": _yxz_row(t.border_size),
    }


def tile_blocks_to_matlab(blocks: Sequence[Sequence[TileBlock]]) -> list:
    """the batches' block information as MATLAB's cell of struct arrays (``[]`` for a batch
    without tiles, a struct for one tile), for ``jsonencode``"""
    out: list = []
    for tiles in blocks:
        if not tiles:
            out.append(None)
        elif len(tiles) == 1:
            out.append(_tile_to_matlab(tiles[0]))
        else:
            out.append([_tile_to_matlab(t) for t in tiles])
    return out


def block_info_json(blocks: Sequence[Sequence[TileBlock]]) -> str:
    """the text of a per-block information file: ``jsonencode(stitchBlockInfo_t,
    'PrettyPrint', true)``"""
    from petakit5d.core.jsonenc import jsonencode

    return jsonencode(tile_blocks_to_matlab(blocks), pretty_print=True)


def read_per_block_info(path: str | os.PathLike) -> list[list[TileBlock]]:
    """the block information of a blending task's batches from its file: MATLAB's JSON
    (``jsondecode``), or a ``.mat`` file with ``stitchBlockInfo_t``"""
    path = os.fspath(path)
    if os.path.splitext(path)[1] == ".mat":
        from petakit5d.io import read_mat

        return tile_blocks_from_matlab(read_mat(path, ["stitchBlockInfo_t"])["stitchBlockInfo_t"])
    with open(path, "rb") as f:
        text = f.read().decode("utf-8")
    return tile_blocks_from_matlab(json.loads(text))


# ---------------------------------------------------------------------------------------------
# the block information file
# ---------------------------------------------------------------------------------------------


def _object_column(items: Sequence[Any]) -> np.ndarray:
    a = np.empty((len(items), 1), dtype=object)
    for k, v in enumerate(items):
        a[k, 0] = v
    return a


def _regions_cell(regions: Sequence[Sequence[Box]]) -> np.ndarray:
    """pairs of half-overlap boxes as MATLAB's K x 1 cell of 1 x 2 cells of [s t] regions"""
    from petakit5d.stitch.geometry import region_to_matlab

    return _object_column([[region_to_matlab(a), region_to_matlab(b)] for a, b in regions])


def _regions_from_cell(v: Any) -> list[tuple[Box, Box]]:
    from petakit5d.stitch.geometry import region_from_matlab

    if isinstance(v, np.ndarray) and v.dtype == object:
        items = list(v.ravel(order="F"))
    elif isinstance(v, list | tuple):
        items = list(v)
    else:
        items = []
    out = []
    for p in items:
        q = list(p.ravel(order="F")) if isinstance(p, np.ndarray) and p.dtype == object else p
        out.append((region_from_matlab(q[0]), region_from_matlab(q[1])))  # type: ignore[arg-type]
    return out


def _data_order_rev(data_order: str) -> list[int]:
    """for each MATLAB dimension (y, x, z of the standard order), the column of x, y, z it
    holds (MATLAB's data_order_reverse_mat, 0-based)"""
    from petakit5d.stitch.geometry import data_order_axes_to_matlab, parse_data_order

    dims = data_order_axes_to_matlab(parse_data_order(data_order))
    return [int(k) for k in np.argsort(np.asarray(dims), kind="stable")]


def write_block_info(
    path: str,
    *,
    tile_fns: Sequence[str],
    int_xyz_shift,
    im_sizes,
    overlap_matrix,
    half_regions: Sequence[Sequence[Box]],
    full_regions: Sequence[Sequence[Box]],
    pairs,
    canvas_boxes,
    tile_boxes,
    valid,
    nv_size: Sequence[int],
    block_size: Sequence[int],
    border_size: Sequence[int] | None,
    data_order: str = "z,y,x",
) -> None:
    """the block information file (``save('-v7.3', ...)``) with MATLAB's values; see the
    module (``int_xyz_shift`` in x, y, z columns, permuted to the tiles' axis order)."""
    from petakit5d.core.axes import to_yxz
    from petakit5d.io import write_mat
    from petakit5d.stitch.geometry import overlap_map_to_matlab, tile_canvas_boxes_to_matlab

    sh = np.asarray(int_xyz_shift, dtype=np.float64).reshape(-1, 3)
    n = sh.shape[0]
    sizes = np.asarray(im_sizes, dtype=np.float64).reshape(-1, 3)
    bboxes, bstart, bend = tile_canvas_boxes_to_matlab(canvas_boxes, tile_boxes, valid)
    nv = [int(v) for v in nv_size]
    bs = [int(v) for v in block_size]
    b_sub = [-(-a // b) for a, b in zip(nv, bs, strict=True)]
    data = {
        "tileFns": _object_column([os.fspath(f) for f in tile_fns]),
        "int_xyz_shift": sh[:, _data_order_rev(data_order)],
        "imSizes": sizes[:, [1, 2, 0]],
        "overlap_matrix": np.asarray(overlap_matrix, dtype=bool).reshape(n, n),
        "half_ol_region_cell": _regions_cell(half_regions),
        "ol_region_cell": _regions_cell(full_regions),
        "overlap_map_mat": overlap_map_to_matlab(pairs, n).reshape(-1, 2),
        "bboxes": bboxes,
        "bboxStart_mat": bstart,
        "bboxEnd_mat": bend,
        "nvSize": np.asarray(to_yxz(tuple(nv)), dtype=np.float64),
        "blockSize": np.asarray(to_yxz(tuple(bs)), dtype=np.float64),
        "bSubSz": np.asarray(to_yxz(tuple(b_sub)), dtype=np.float64),
        "BorderSize": np.zeros((0, 0)) if border_size is None else _yxz_row(border_size),
    }
    write_mat(path, data, version="7.3")


def read_block_info(path: str | os.PathLike) -> dict[str, Any]:
    """the block information file in the port's forms: ``tile_fns``, ``overlap_matrix``,
    ``half_regions`` and ``full_regions`` (pairs of boxes), ``pairs`` (``(K, 2)``, from
    ``overlap_map_mat``), ``canvas_boxes``, ``tile_boxes``, ``valid``, ``nv_size``,
    ``block_size``, ``b_sub_sz`` (``(z, y, x)``), ``border_size`` and ``im_sizes``."""
    from petakit5d.core.bbox import boxes_from_matlab
    from petakit5d.io import read_mat

    a = read_mat(os.fspath(path))
    bboxes = np.asarray(a["bboxes"], dtype=np.float64).reshape(-1, 6)
    bstart = np.asarray(a["bboxStart_mat"], dtype=np.float64).reshape(-1, 3)
    bend = np.asarray(a["bboxEnd_mat"], dtype=np.float64).reshape(-1, 3)
    n = bboxes.shape[0]
    om = np.asarray(a["overlap_matrix"], dtype=bool).reshape(n, n)
    omap = np.asarray(a["overlap_map_mat"], dtype=np.float64).reshape(-1, 2)
    order = np.argsort(omap[:, 0], kind="stable")
    ind = omap[order, 1].astype(np.int64) - 1
    pairs = np.column_stack((ind % n, ind // n)).astype(np.int64)
    fns = a.get("tileFns", [])
    fns = list(fns.ravel(order="F")) if isinstance(fns, np.ndarray) else list(fns)
    sizes = np.asarray(a.get("imSizes", np.zeros((0, 3))), dtype=np.float64).reshape(-1, 3)
    return {
        "tile_fns": [str(f) for f in fns],
        "overlap_matrix": om,
        "half_regions": _regions_from_cell(a["half_ol_region_cell"]),
        "full_regions": _regions_from_cell(a.get("ol_region_cell", [])),
        "pairs": pairs,
        "canvas_boxes": boxes_from_matlab(bboxes),
        "tile_boxes": boxes_from_matlab(np.hstack((bstart, bend))),
        "valid": ~np.all(bstart == 0, axis=1),
        "nv_size": _zyx(a["nvSize"]),
        "block_size": _zyx(a["blockSize"]),
        "b_sub_sz": _zyx(a["bSubSz"]),
        "border_size": _zyx(a.get("BorderSize")),
        "im_sizes": sizes[:, [2, 0, 1]].astype(np.int64),
    }


# ---------------------------------------------------------------------------------------------
# organizing: the tiles of each batch
# ---------------------------------------------------------------------------------------------


def _closed(boxes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """boxes as closed 1-based (start, end) in (z, y, x); zero rows stay zero"""
    b = np.asarray(boxes, dtype=np.int64).reshape(-1, 6)
    lo, hi = b[:, :3] + 1, b[:, 3:].copy()
    zero = ~np.any(b != 0, axis=1)
    lo[zero] = 0
    hi[zero] = 0
    return lo, hi


def _matlab_overlap(lo1, hi1, lo2, hi2) -> np.ndarray:
    """MATLAB's test of stitch_organize_block_info (closed boxes; per axis ol_1 | ol_2 | ol_3 |
    ol_4, where ol_4 repeats ol_3): rows of the first set against one box of the second"""
    ol_1 = (lo1 <= lo2) & (lo2 <= hi1)
    ol_2 = (lo1 <= hi2) & (hi2 <= hi1)
    ol_3 = (lo2 <= hi1) & (hi1 <= hi2)
    ol_4 = (lo2 <= hi1) & (hi1 <= hi2)
    return np.all(ol_1 | ol_2 | ol_3 | ol_4, axis=1)


def _batch_boxes(first: int, last: int, nv_size, block_size) -> np.ndarray:
    """the boxes of batches first to last (0-based, MATLAB's order) in the stitched volume"""
    nv = np.asarray(nv_size, dtype=np.int64)
    bs = np.asarray(block_size, dtype=np.int64)
    sub = -(-nv // bs)
    k = np.arange(first, last + 1, dtype=np.int64)
    # MATLAB's ind2sub(bSubSz, k): y fastest, then x, then z
    sy = k % sub[1]
    sx = (k // sub[1]) % sub[2]
    sz = k // (sub[1] * sub[2])
    s = np.column_stack((sz, sy, sx))
    start = s * bs
    stop = np.minimum(nv, (s + 1) * bs)
    return np.hstack((start, stop))


def organize_blocks(
    block_range: tuple[int, int],
    canvas_boxes,
    tile_boxes,
    valid,
    overlap_matrix,
    half_regions: Sequence[Sequence[Box]],
    pairs,
    nv_size: Sequence[int],
    block_size: Sequence[int],
    border_size: Sequence[int] | None = (0, 0, 0),
) -> list[list[TileBlock]]:
    """The tiles of batches ``block_range[0]`` to ``block_range[1]`` (0-based, inclusive, as
    MATLAB's ``blockInds``; batches of ``block_size`` over ``nv_size``, ``(z, y, x)``), as
    ``stitch_organize_block_info`` computes them: each tile whose box in the stitched volume
    (``canvas_boxes``; ``tile_boxes`` the same box in the tile, ``valid`` False for tiles
    outside) overlaps the batch, in tile order, with its blocks (``core.blocks.bbox_to_blocks``
    with ``border_size``) and the parts of its half-overlap regions in the batch: the region of
    each pair the tile is in (``half_regions[k]`` for ``pairs[k]``, the first box for the
    pair's first tile) cut to the tile's box."""
    from petakit5d.core.blocks import bbox_to_blocks

    cb = np.asarray(canvas_boxes, dtype=np.int64).reshape(-1, 6)
    tb = np.asarray(tile_boxes, dtype=np.int64).reshape(-1, 6)
    ok = np.asarray(valid, dtype=bool).reshape(-1)
    n = cb.shape[0]
    om = np.asarray(overlap_matrix, dtype=bool).reshape(n, n)
    pair_of = {(int(i), int(j)): k for k, (i, j) in enumerate(np.asarray(pairs).reshape(-1, 2))}
    border = None if border_size is None else tuple(int(v) for v in border_size)
    first, last = int(block_range[0]), int(block_range[1])
    blk = _batch_boxes(first, last, nv_size, block_size)
    blo, bhi = _closed(blk)

    # map each block with the tiles (MATLAB's bboxes are zeros for tiles outside)
    clo, chi = _closed(np.where(ok[:, None], cb, 0))
    block_tiles = [
        np.nonzero(_matlab_overlap(clo, chi, blo[b], bhi[b]))[0] for b in range(len(blk))
    ]
    used = np.zeros(n, dtype=bool)
    for t in block_tiles:
        used[t] = True

    tile_info: dict[int, tuple[Any, list[Any], np.ndarray]] = {}
    for i in range(n):
        if not used[i] or not ok[i]:
            continue
        if border is None:
            # (bboxToBlocks indexes the border, borderSize(i), of [])
            raise BlockInfoError(
                "Index exceeds the number of array elements. Index must not exceed 0.",
                "MATLAB:badsubscript",
            )
        start = tuple(int(v) for v in tb[i, :3])
        stop = tuple(int(v) for v in tb[i, 3:])
        info = bbox_to_blocks(tuple(cb[i]), start, block_size, nv_size, border)
        minfo: list[Any] = []
        mboxes: list[tuple[int, ...]] = []
        for j in range(n):
            if not om[i, j] and not om[j, i]:
                continue
            k = pair_of.get((i, j) if i < j else (j, i))
            if k is None:
                raise BlockInfoError(
                    "Index exceeds the number of array elements. Index must not exceed 0.",
                    "MATLAB:badsubscript",
                )
            region = half_regions[k][0 if i < j else 1]
            lo = [max(int(region[q]), start[q]) for q in range(3)]
            hi = [min(int(region[q + 3]), stop[q]) for q in range(3)]
            if any(hi[q] <= lo[q] for q in range(3)):
                continue
            # the region in the stitched volume
            off = [start[q] - int(cb[i, q]) for q in range(3)]
            mbox = tuple(lo[q] - off[q] for q in range(3)) + tuple(hi[q] - off[q] for q in range(3))
            minfo.append(bbox_to_blocks(mbox, tuple(lo), block_size, nv_size, border))
            mboxes.append(mbox)
        tile_info[i] = (info, minfo, np.asarray(mboxes, dtype=np.int64).reshape(-1, 6))

    def rows(info, b: int) -> tuple:
        sel = np.nonzero(info.block_ind == b)[0]
        if sel.size == 0:
            return None, None, None
        r = int(sel[0])
        return (
            tuple(int(v) for v in info.b_coords[r]),
            tuple(int(v) for v in info.w_coords[r]),
            tuple(int(v) for v in info.bbox_coords[r]),
        )

    out: list[list[TileBlock]] = []
    for b in range(len(blk)):
        bind = first + b
        tiles = []
        for ti in block_tiles[b]:
            info, minfo, mboxes = tile_info[int(ti)]
            bc, wc, bbc = rows(info, bind)
            mlo, mhi = _closed(mboxes)
            sel = np.nonzero(_matlab_overlap(mlo, mhi, blo[b], bhi[b]))[0] if len(minfo) else []
            regions = tuple(RegionBlock(*rows(minfo[int(k)], bind)) for k in sel)
            tiles.append(TileBlock(int(ti), bc, wc, bbc, regions, border))  # type: ignore[arg-type]
        out.append(tiles)
    return out


def _save_flag(flag_fullname: str, block_inds: Sequence[int], t1: float, uuid: str) -> None:
    """``save(flatTmpFn, 'blockInds', 't1'); movefile(flatTmpFn, flagFullname)``"""
    from petakit5d.io import write_mat

    tmp = f"{flag_fullname[:-4]}_{uuid}.mat"
    write_mat(
        tmp,
        {"blockInds": np.asarray(block_inds, dtype=np.float64) + 1, "t1": float(t1)},
        version="7",
    )
    os.replace(tmp, flag_fullname)


def _block_inds(v: Any) -> list[int]:
    a = np.asarray(v, dtype=np.float64).reshape(-1)
    if a.size == 0 or not np.all(np.isfinite(a)) or np.any(a != np.floor(a)):
        raise BlockInfoError(f"block indices must be integers, got {v!r}")
    return [int(x) for x in a]


@register(
    "stitch_organize_block_info",
    Param(
        "blockInds", "@isnumeric", "First and last batch (0-based) of the task.", convert="index1"
    ),
    Param(
        "taskSize", "@isnumeric", "Batches per blending task (one JSON file each).", convert="same"
    ),
    Param("BlockInfoFullname", "@(x) ischar(x)", "The block information file."),
    Param("PerBlockInfoPath", "@(x) ischar(x)", "Folder of the per-block information files."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when the task is done."),
    Param("Overwrite", "@islogical", "Redo a task whose flag exists."),
    Param("BorderSize", "@isnumeric", "Border of the batches (z, y, x).", convert="size_yxz"),
    Param("uuid", "@ischar", "Id of the temporary files (empty: a new one)."),
    task="S07",
)
def organize_block_info(
    block_inds: Sequence[int],
    task_size: int,
    block_info_fullname: str,
    per_block_info_path: str,
    flag_fullname: str,
    *,
    overwrite: bool = False,
    border_size: Sequence[int] | None = None,
    uuid: str = "",
) -> bool:
    """``stitch_organize_block_info``: the per-block information of batches ``block_inds[0]``
    to ``block_inds[-1]`` (0-based) from the block information file, written as one JSON file
    per ``task_size`` batches (``stitch_block_info_blocks_<s>_<t>.json``, 1-based numbers; a
    temporary name, then a rename), then the task's flag (``blockInds``, ``t1``). Nothing is
    done if the flag exists (``overwrite`` removes it). Returns whether the task ran."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.jobs.files import exist

    t0 = time.perf_counter()
    if not uuid:
        uuid = get_uuid()
    inds = _block_inds(block_inds)
    flag_fullname = os.fspath(flag_fullname)
    if exist(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            print(
                f"The block files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!",
                flush=True,
            )
            return False
    a = read_block_info(block_info_fullname)
    first, last = inds[0], inds[-1]
    blocks = organize_blocks(
        (first, last),
        a["canvas_boxes"],
        a["tile_boxes"],
        a["valid"],
        a["overlap_matrix"],
        a["half_regions"],
        a["pairs"],
        a["nv_size"],
        a["block_size"],
        None if border_size is None else tuple(int(v) for v in border_size),
    )
    ts = int(task_size)
    nb = last - first + 1
    for t in range(math.ceil(nb / ts)):
        st, tt = t * ts, min((t + 1) * ts, nb)
        bs, bt = st + first + 1, tt + first
        fn = f"{per_block_info_path}/stitch_block_info_blocks_{bs}_{bt}.json"
        tmp = f"{fn[:-5]}_{uuid}.json"
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            f.write(block_info_json(blocks[st:tt]))
        os.replace(tmp, fn)
    _save_flag(flag_fullname, [first, last], time.perf_counter() - t0, uuid)
    return True


# ---------------------------------------------------------------------------------------------
# the block information of a volume, with its organizing tasks
# ---------------------------------------------------------------------------------------------


def block_task_size(num_blocks: int, task_size: int, max_file_num_per_folder: int) -> int:
    """the batches of one organizing task (MATLAB's ``blk_taskSize``): 500 blending tasks,
    within 1e5 to 1e6 batches, rounded up to whole blending tasks; with more blending tasks than
    ``max_file_num_per_folder``, at most the batches of one folder of JSON files ``F``, or if
    fewer, ``F / 2^(i - 1)`` for the first ``i`` (from 1 to ``ceil(log2(F / size))``) with
    ``size > F / 2^i`` or ``F / 2^i`` not whole, and the last ``i`` without one (MATLAB's loop:
    ``F / 2`` for a size of exactly ``F / 4``)."""
    ts = int(task_size)
    n_file = int(max_file_num_per_folder)
    blk = max(1e5, min(1e6, ts * 500))
    blk = math.ceil(blk / ts) * ts
    if math.ceil(num_blocks / ts) > n_file:
        blk = min(blk, n_file * ts)
        if blk < n_file * ts:
            i = 0
            for i in range(1, math.ceil(math.log2(n_file * ts / blk)) + 1):
                cur = n_file * ts / 2**i
                if blk > cur or math.floor(cur) != cur:
                    break
            blk = n_file * ts / 2 ** (i - 1)
    return int(blk)


def per_block_info_paths(
    per_block_info_path: str, task_size: int, num_blocks: int, max_file_num_per_folder: int
) -> tuple[list[str], list[str]]:
    """the folders of the JSON files and the JSON file of each blending task, as MATLAB names
    them (more tasks than ``max_file_num_per_folder``: subfolders
    ``task_size_<k>_tasks_<a>_<b>/``, which are not made here)"""
    ts, n_file = int(task_size), int(max_file_num_per_folder)
    num_tasks = math.ceil(num_blocks / ts)
    folders = [per_block_info_path]
    if num_tasks > n_file:
        folders = [
            f"{per_block_info_path}/task_size_{ts}_tasks_{t * n_file + 1}_"
            f"{min((t + 1) * n_file, num_tasks)}/"
            for t in range(math.ceil(num_tasks / n_file))
        ]
    paths = []
    for t in range(num_tasks):
        bs, bt = t * ts + 1, min((t + 1) * ts, num_blocks)
        folder = folders[t // n_file] if num_tasks > n_file else per_block_info_path
        paths.append(f"{folder}/stitch_block_info_blocks_{bs}_{bt}.json")
    return folders, paths


_JSON_START = re.compile(r"stitch_block_info_blocks_(\d+)_\d+.json")


def block_info_bytes(folders: Sequence[str]) -> list[int]:
    """the sizes of the per-block information files of the folders, in the order of their
    first batch (MATLAB's ``get_block_info_file_bytes``)"""
    from petakit5d.core.paths import matlab_dir

    starts: list[float] = []
    sizes: list[int] = []
    for folder in folders:
        for e in matlab_dir(f"{folder}*.json"):
            m = _JSON_START.search(e.name)
            if m is None:
                raise BlockInfoError(
                    "Index exceeds the number of array elements. Index must not exceed 0.",
                    "MATLAB:badsubscript",
                )
            starts.append(float(m.group(1)))
            sizes.append(os.path.getsize(os.path.join(e.folder, e.name)))
    order = np.argsort(np.asarray(starts), kind="stable")
    return [sizes[int(k)] for k in order]


def _mkdir(path: str) -> None:
    # MATLAB's mkdir: parents too, an existing folder is fine (jobs share these folders)
    os.makedirs(path, exist_ok=True)


def process_block_info(
    int_xyz_shift,
    im_sizes,
    nv_size: Sequence[int],
    block_size: Sequence[int],
    data_order: str,
    overlap_matrix,
    full_regions: Sequence[Sequence[Box]],
    half_regions: Sequence[Sequence[Box]],
    pairs,
    border_size: Sequence[int] | None,
    tile_fns: Sequence[str],
    stitch_info_path: str,
    nv_fsname: str,
    is_primary_ch: bool,
    *,
    stitch_info_fullpath: str = "",
    stitch_2d: bool = False,
    task_size: int = 10,
    max_file_num_per_folder: int = 20000,
    uuid: str = "",
    parse_cluster: bool = True,
    mcc_mode: bool = False,
    config_file: str = "",
) -> tuple[str, list[str], list[int] | None]:
    """``stitch_process_block_info``: the block information of a stitched volume (batches of
    ``block_size`` over ``nv_size``, ``(z, y, x)``), from the tiles' voxel offsets
    ``int_xyz_shift`` (x, y, z columns), their sizes ``im_sizes`` (``(n, 3)`` in array order,
    axes named by ``data_order``), the rechecked ``overlap_matrix`` and the full and half
    overlap regions of ``pairs`` (``stitch.geometry.half_overlap_regions``): the block
    information file, and for a primary channel the per-block JSON files of every blending
    task of ``task_size`` batches, made by organizing tasks (``organize_block_info``) through
    the job framework (2 CPUs, the master computes, up to three rounds). A secondary channel
    uses the primary's JSON files (the folder named after ``stitch_info_fullpath``, which must
    exist).

    Returns the block information file, the JSON file of each blending task and, with
    ``parse_cluster``, the sizes of the JSON files (in task order; else None)."""
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.stitch.geometry import tile_canvas_boxes

    nv = [int(v) for v in nv_size]
    bsz = [int(v) for v in block_size]
    ts = int(task_size)
    n_file = int(max_file_num_per_folder)
    print("Process stitch block information...", flush=True)
    b_sub = [-(-a // b) for a, b in zip(nv, bsz, strict=True)]
    num_blocks = math.prod(b_sub)
    cb, tb, valid = tile_canvas_boxes(int_xyz_shift, im_sizes, nv, data_order, stitch_2d)

    tmp_fn = f"{stitch_info_path}/{nv_fsname}_block_info_task_size_{ts}_{uuid}.mat"
    block_info_fullname = f"{stitch_info_path}/{nv_fsname}_block_info_task_size_{ts}.mat"
    write_block_info(
        tmp_fn,
        tile_fns=tile_fns,
        int_xyz_shift=int_xyz_shift,
        im_sizes=im_sizes,
        overlap_matrix=overlap_matrix,
        half_regions=half_regions,
        full_regions=full_regions,
        pairs=pairs,
        canvas_boxes=cb,
        tile_boxes=tb,
        valid=valid,
        nv_size=nv,
        block_size=bsz,
        border_size=border_size,
        data_order=data_order,
    )
    os.replace(tmp_fn, block_info_fullname)

    if is_primary_ch:
        pbi_path = f"{stitch_info_path}/{nv_fsname}_task_size_{ts}/"
        _mkdir(pbi_path)
        flag_path = f"{pbi_path}/block_flags/"
        _mkdir(flag_path)
    else:
        from petakit5d.core.paths import fileparts

        pstr, fsn, _ = fileparts(stitch_info_fullpath)
        pbi_path = f"{pstr}/{fsn}_task_size_{ts}/"
        if not os.path.isdir(pbi_path):
            raise BlockInfoError(
                f"The Per Block Info path for the primary channel {pbi_path} does not exist!"
            )
    folders, pbi_fullpaths = per_block_info_paths(pbi_path, ts, num_blocks, n_file)
    if len(folders) > 1 or folders[0] != pbi_path:
        for f in folders:
            _mkdir(f)

    if not is_primary_ch:
        return (
            block_info_fullname,
            pbi_fullpaths,
            block_info_bytes(folders) if parse_cluster else None,
        )

    num_tasks = math.ceil(num_blocks / ts)
    blk = block_task_size(num_blocks, ts, n_file)
    blk_num_tasks = math.ceil(num_blocks / blk)
    flags, tasks = [], []
    border = None if border_size is None else [int(v) for v in border_size]
    for t in range(blk_num_tasks):
        bs, bt = t * blk + 1, min((t + 1) * blk, num_blocks)
        flag = f"{flag_path}/stitch_block_info_task_size_{ts}_blocks_{bs}_{bt}.mat"
        flags.append(flag)
        folder = folders[math.ceil(bs / (n_file * ts)) - 1] if num_tasks > n_file else pbi_path
        tasks.append(
            Task.of(
                organize_block_info,
                [bs - 1, bt - 1],
                ts,
                block_info_fullname,
                folder,
                flag,
                border_size=border,
            )
        )
    done: list[bool] = []
    for _ in range(3):
        done = framework.run_tasks(
            [block_info_fullname] * blk_num_tasks,
            flags,
            tasks,
            cpus_per_task=2,
            mem_allocate=blk * 0.0001,
            min_task_job_num=2 if blk_num_tasks >= 4 else 1,
            max_trial_num=2,
            master_compute=True,
            parse_cluster=bool(parse_cluster),
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        if all(done):
            break
    if not all(done):
        raise BlockInfoError(
            f"Some block ({sum(not d for d in done)} / {len(done)}) info file cannot finished!"
        )
    return block_info_fullname, pbi_fullpaths, block_info_bytes(folders) if parse_cluster else None
