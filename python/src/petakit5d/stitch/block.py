"""The blending tasks of stitching (``processStitchBlock``), the merge of their batches' z MIPs
(``stitch_merge_block_mips``), the frame's memory decisions (``stitch_mip_memory``,
``stitch_in_process``) and the separate-tile output (``processStitchSeparteTiles``).

``process_stitch_block`` blends its batches of the stitched volume from the tiles' zarr files
and writes each batch into the stitched zarr (made by the frame), then its flag: with
``save_mip`` the batches' z MIPs (one per batch position in y and x, MATLAB's ``.mat`` with
``batchInds``, ``mipBboxes``, ``mips``), else an empty file. As processStitchBlock.m at the
oracle commit:

- nothing is done if the flag's folder is missing (all batches were done and it was removed)
  or the flag exists (``overwrite`` removes it);
- a batch without tiles is zeros; a single tile that covers the batch is written as it is;
  tiles that are not of ``dtype`` are cast to it (MATLAB's ``cast``);
- 'none', 'max', 'mean' and 'median' blend with ``kernels.stitch.tile_blending_3d``, each tile
  counted as 0 in its half-overlap regions (its ``RegionBlock``s; a region outside the tile's
  part of the batch is an error);
- 'feather' blends the tiles with their distance maps (``imdist_fullpaths``, the map of each
  tile by ``imdist_file_idx``; maps of at most two batches in all are read whole once): in the
  region where the tiles overlap (the box of the pairwise overlaps of their boxes, where a
  one-voxel overlap counts as none, widened to tiles that miss it), or in the whole batch when
  that region is 80% of it or more ("full"); outside it the last tile read wins. With pooled
  maps (``pool_size``, large files) the maps are resized to each tile's part
  (``kernels.stitch.feather_distance_map_resize_3d``), and with more than two tiles the
  "minor" ones (whose largest weight is below the smallest weight of a tile that covers the
  region by more than 1000^(1/10)) are left out unless the major tiles leave voxels without
  data; one major tile that covers the others is written as it is;
- the result is cut to the batch (borders, the volume's end) and written in the zarr's class;
  the MIP of a blended batch is ``cast(max(nv_block, [], 3), dtype)``.

Batches are numbered from 0 (MATLAB's order: y fastest, then x, then z); boxes are ``(z0, y0,
x0, z1, y1, x1)``, 0-based, end-exclusive; sizes, the border and the pool size are ``(z, y,
x)``. The flags keep MATLAB's values (1-based batch numbers, ``[ymin xmin ymax xmax]`` MIP
boxes). Reads go through one buffer per kind (tiles, distance maps), reused across batches.

With pooled maps the voxels the major tiles leave without data are read from the minor tiles
(their box in the batch's coordinates), as are the minor tiles' parts the major tiles don't
cover (with the major tile's weight); 'feather' is matched in any case; with a border
(``border_size`` > 0) the batch is cut to its own size.
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "StitchBlockError",
    "in_process",
    "merge_block_mips",
    "mip_memory",
    "process_separate_tiles",
    "process_stitch_block",
    "save_block_flag",
]

Box = tuple[int, int, int, int, int, int]

# MATLAB's 1000 .^ (1/10) (the feather blending's dominance factor dmtFactor .^ (1/10); its
# bits from MATLAB R2024b, 0x3fffec982d5bb8b0)
_DMT = float.fromhex("0x1.fec982d5bb8b0p+0")
# MATLAB's whos: bytes of a cell element's header (R2024b, 64-bit)
_CELL_HEADER_BYTES = 112


class StitchBlockError(RuntimeError):
    """An error MATLAB raises; ``identifier`` is MATLAB's error id (``''`` for none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _ints(v: Any, what: str = "values") -> list[int]:
    a = np.asarray(v if v is not None else [], dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(a)) or np.any(a != np.floor(a)):
        raise StitchBlockError(f"{what} must be integers, got {v!r}")
    return [int(x) for x in a]


def _ndtype(dtype: str) -> np.dtype:
    from petakit5d.core.arrays import numpy_dtype

    return numpy_dtype(dtype)


def _as_class(a: np.ndarray, dtype: str) -> np.ndarray:
    """MATLAB's ``if ~isa(a, dtype), a = cast(a, dtype); end``"""
    if a.dtype == _ndtype(dtype):
        return a
    from petakit5d.core import dtypes

    return np.ascontiguousarray(dtypes.cast(a, dtype))


def _vol3(a: np.ndarray, box: Sequence[int]) -> np.ndarray:
    """a read of a box as (z, y, x) (an image of one slice comes back 2D)"""
    return a.reshape(tuple(int(box[q + 3]) - int(box[q]) for q in range(3)))


def _read(path: str, box: Sequence[int], buf: Any = None, nt: int | None = None) -> np.ndarray:
    from petakit5d.io import read_zarr

    return _vol3(read_zarr(path, input_bbox=list(box), out=buf, num_threads=nt), box)


def _proper(box: Sequence[int], shape: Sequence[int]) -> bool:
    return all(0 <= box[q] <= box[q + 3] <= shape[q] for q in range(len(shape)))


def _m4(box: Sequence[int], j: int) -> list[float]:
    """our 3D box and a tile number as MATLAB's 4D [ymin xmin zmin t ymax xmax zmax t]"""
    return [box[1] + 1, box[2] + 1, box[0] + 1, j + 1, box[4], box[5], box[3], j + 1]


def _put(
    dst: np.ndarray,
    j: int,
    box: Sequence[int],
    src: np.ndarray,
    src_box: Sequence[int] | None = None,
    nt: int | None = None,
) -> np.ndarray:
    """``dst = indexing4d(dst, src, [box(1:3), j, box(4:6), j][, rbbox])`` (``dst`` ``(t, z,
    y, x)``, in place): the copy kernel for boxes that fit, MATLAB's indexing4d (K01) for the
    empty and reversed boxes MATLAB's arithmetic can give"""
    from petakit5d.kernels.crop import indexing4d_matlab, indexing_4d

    ext = [box[q + 3] - box[q] for q in range(3)]
    sb = [0, 0, 0, *src.shape] if src_box is None else list(src_box)
    if (
        src.dtype == dst.dtype
        and 0 <= j < dst.shape[0]
        and _proper(box, dst.shape[1:])
        and _proper(sb, src.shape)
        and ext == [sb[q + 3] - sb[q] for q in range(3)]
    ):
        rb = None if src_box is None else (0, *sb[:3], 1, *sb[3:])
        indexing_4d(dst, (j, *box[:3], j + 1, *box[3:]), src[np.newaxis], rb, num_threads=nt)
        return dst
    rbm = None if src_box is None else _m4(sb, 0)
    return indexing4d_matlab(dst, src[np.newaxis], _m4(box, j), rbm, num_threads=nt)


def _crop3(a: np.ndarray, box: Sequence[int], nt: int | None = None) -> np.ndarray:
    """``crop3d(a, box)`` (MATLAB's errors for a box beyond the array)"""
    from petakit5d.kernels.crop import crop3d_matlab, crop_3d

    if _proper(box, a.shape):
        return crop_3d(a, list(box), num_threads=nt)
    m = [box[1] + 1, box[2] + 1, box[0] + 1, box[4], box[5], box[3]]
    return _vol3(crop3d_matlab(a, m, num_threads=nt), box)


def _batch_sub(bi: int, b_sub: Sequence[int]) -> tuple[int, int, int]:
    """MATLAB's ind2sub(bSubSz, bi) of a 0-based batch number: (z, y, x), 0-based"""
    nz, ny, nx = b_sub
    return bi // (ny * nx), bi % ny, (bi // ny) % nx


def _pool_psz(pool_size: Sequence[int]) -> tuple[int, int, int] | None:
    """the pooling of the distance maps in z, y, x: MATLAB's psz = [1, 1, poolSize(3)] (3
    values) or poolSize([4, 5, 3]) (6 or 9); None for other counts (undefined in MATLAB)"""
    p = list(pool_size)
    if len(p) == 3:
        return (p[0], 1, 1)
    if len(p) in (6, 9):
        # (z, y, x) groups: MATLAB's [py1 px1 pz1 py2 px2 pz2]([4, 5, 3]) = [py2 px2 pz1]
        return (p[0], p[4], p[5])
    return None


def _mround(x: float) -> float:
    from petakit5d.core.dtypes import mround

    return float(mround(x))


def _closed(box: Sequence[int]) -> list[int]:
    """our box as a closed 1-based box in (z, y, x) (for check_cuboids_overlaps)"""
    return [box[0] + 1, box[1] + 1, box[2] + 1, box[3], box[4], box[5]]


def _overlap_region(boxes: list[Box], bsz: Sequence[int]) -> tuple[list[int], bool, list[int]]:
    """the region where the tiles of a batch overlap (MATLAB's ``bCoords_c``), whether the
    batch is full (that region is 80% of it or more, or empty) and the region's size"""
    from petakit5d.stitch.geometry import check_cuboids_overlaps

    n = len(boxes)
    if n == 2:
        c = [max(boxes[0][q], boxes[1][q]) for q in range(3)] + [
            min(boxes[0][q + 3], boxes[1][q + 3]) for q in range(3)
        ]
    else:
        # the pairs of find(triu(ones(n), 1)): column by column
        s = [a for t in range(n) for a in range(t)]
        t_ = [t for t in range(n) for _ in range(t)]
        ov: np.ndarray | None
        _, ov = check_cuboids_overlaps(
            np.asarray([_closed(boxes[a]) for a in s], dtype=np.float64),
            np.asarray([_closed(boxes[b]) for b in t_], dtype=np.float64),
            False,
        )
        rows = (
            np.asarray([_closed(b) for b in boxes], dtype=np.float64)
            if ov is None
            else ov[~np.any(ov == 0, axis=1)]
        )
        if rows.shape[0] == 0:
            rows = np.asarray([_closed(b) for b in boxes], dtype=np.float64)
        c = [int(v) - 1 for v in rows[:, :3].min(axis=0)] + [
            int(v) for v in rows[:, 3:].max(axis=0)
        ]
        # tiles that miss the region widen it (rare)
        miss = [
            k
            for k, b in enumerate(boxes)
            if any(b[q] >= c[q + 3] or b[q + 3] <= c[q] for q in range(3))
        ]
        if miss:
            c = [min([boxes[k][q] for k in miss] + [c[q]]) for q in range(3)] + [
                max([boxes[k][q + 3] for k in miss] + [c[q + 3]]) for q in range(3)
            ]
    dbsz = [c[q + 3] - c[q] for q in range(3)]
    # prod(dbsz ./ bsz) in MATLAB's order (y, x, z)
    r = [dbsz[q] / bsz[q] for q in (1, 2, 0)]
    full = (r[0] * r[1]) * r[2] >= 0.8 or any(v <= 0 for v in dbsz)
    if full:
        c = [0, 0, 0, *bsz]
        dbsz = list(bsz)
    return c, full, dbsz


def _cut(b: Sequence[int], c: Sequence[int]) -> list[int]:
    """the part of box b in box c: [max(c(1:3), b(1:3)), min(c(4:6), b(4:6))]"""
    return [max(c[q], b[q]) for q in range(3)] + [min(c[q + 3], b[q + 3]) for q in range(3)]


def _rel(b: Sequence[int], c: Sequence[int]) -> list[int]:
    """box b relative to the start of box c"""
    return [b[q] - c[q] for q in range(3)] + [b[q + 3] - c[q] for q in range(3)]


def _range_max_z(a: np.ndarray, nt: int | None) -> np.ndarray:
    """MATLAB's max(a, [], 3) of a (z, y, x) volume as a (y, x) image"""
    from petakit5d.tools.mip import axis_max

    return np.ascontiguousarray(axis_max(a, -3, num_threads=nt).reshape(a.shape[1:]))


class _Maps:
    """the distance maps of a blending task: sizes, the whole maps read once (cacheDist) and
    the buffer of the regions read per tile and batch"""

    def __init__(self, paths: Sequence[str], file_idx: Sequence[int], cache_limit: float | None):
        from petakit5d.io import ReadBuffer, get_image_size

        self.paths = [os.fspath(p) for p in paths]
        self.idx = list(file_idx)
        self.sizes = []
        for p in self.paths:
            s = [int(v) for v in get_image_size(p)]
            # (z, y, x); a 2D map is one slice
            self.sizes.append(tuple(([1] * (3 - len(s)) + s)[-3:]))
        total = sum(math.prod(s) for s in self.sizes)
        self.cache = cache_limit is not None and total <= cache_limit
        self.whole: dict[int, np.ndarray] = {}
        self.buf = ReadBuffer()

    def index(self, tile: int) -> int:
        if len(self.paths) == 1:
            return 0
        if tile >= len(self.idx):
            raise StitchBlockError(
                f"Index exceeds the number of array elements. Index must not exceed "
                f"{len(self.idx)}.",
                "MATLAB:badsubscript",
            )
        return int(self.idx[tile])

    def size(self, tile: int) -> tuple[int, int, int]:
        return self.sizes[self.index(tile)]

    def read(self, tile: int, box: Sequence[int], nt: int | None) -> np.ndarray:
        """the region box of the tile's map (a view of the buffer, valid until the next read)"""
        return _read(self.paths[self.index(tile)], box, self.buf, nt)

    def get_whole(self, tile: int, nt: int | None) -> np.ndarray:
        from petakit5d.io import read_zarr

        d = self.index(tile)
        if d not in self.whole:
            a = read_zarr(self.paths[d], num_threads=nt)
            self.whole[d] = a.reshape(self.sizes[d])
        return self.whole[d]


def _pooled_box(
    dbb: Sequence[int],
    dbsz_j: Sequence[int],
    psz: Sequence[int],
    dsz: Sequence[int],
    one_slice: bool,
) -> list[int]:
    """the box of a pooled distance map for a tile region (MATLAB's p_bboxCoords: start
    max(1, round(start / psz)), end min(dsz, start + ceil(size / psz) - 1), one z slice if
    one_slice), as our box"""
    s = [max(1.0, _mround((dbb[q] + 1) / psz[q])) for q in range(3)]
    e = [min(float(dsz[q]), s[q] + math.ceil(dbsz_j[q] / psz[q]) - 1) for q in range(3)]
    if one_slice:
        e[0] = s[0]
    return [int(v) - 1 for v in s] + [int(v) for v in e]


def _resize(dmat: np.ndarray, size: Sequence[int], wd: Any, nt: int | None) -> np.ndarray:
    from petakit5d.kernels.stitch import feather_distance_map_resize_3d

    return feather_distance_map_resize_3d(
        np.ascontiguousarray(dmat), [0, 0, 0, *size], wd, num_threads=nt
    )


def _weight_power(d: float, wd: Any) -> float:
    """MATLAB's d ^ wd of two doubles"""
    from petakit5d.matlabcompat.elementary import power

    return float(power(np.float64(d), np.float64(float(np.asarray(wd).reshape(-1)[0]))))


# ---------------------------------------------------------------------------------------------
# one batch
# ---------------------------------------------------------------------------------------------


def _feather_weights(
    tiles: list,
    boxes: list[Box],
    c: list[int],
    dbsz: list[int],
    maps: _Maps,
    psz: tuple[int, int, int] | None,
    wd: Any,
    nt: int | None,
) -> tuple[np.ndarray | None, list, np.ndarray, np.ndarray, np.ndarray]:
    """the weights of the tiles in the region c (MATLAB's tim_d_block, ``(t, z, y, x)``),
    with pooled maps and more than two tiles also the maps' regions, their ranges and the minor
    tiles (MATLAB's tim_d_block_cell, d_ranges_mat, minor_inds)"""
    n = len(tiles)
    pooled_many = psz is not None and n > 2
    tim_d = None if pooled_many else np.zeros((n, *dbsz), dtype=np.float32)
    cells: list = [None] * n
    ranges = np.zeros((n, 2))
    dbsz_mat = np.zeros((n, 3), dtype=np.int64)
    for j, tb in enumerate(tiles):
        b = boxes[j]
        db = _cut(b, c)
        t = tb.tile_box
        dbb = [t[q] - b[q] + db[q] for q in range(3)] + [
            t[q + 3] - b[q + 3] + db[q + 3] for q in range(3)
        ]
        dbsz_j = [db[q + 3] - db[q] for q in range(3)]
        dbsz_mat[j] = dbsz_j
        rel = _rel(db, c)
        if psz is None:
            if maps.cache:
                _put(tim_d, j, rel, maps.get_whole(tb.tile, nt), dbb, nt)  # type: ignore[arg-type]
            else:
                _put(tim_d, j, rel, maps.read(tb.tile, dbb, nt), None, nt)  # type: ignore[arg-type]
            continue
        # pooled maps (for now only pooled in z): the region of the map, resized to the tile's
        # (MATLAB: bboxCoords(6) == dbboxCoords(3))
        one_slice = t[3] == dbb[0] + 1
        pb = _pooled_box(dbb, dbsz_j, psz, maps.size(tb.tile), one_slice)
        from petakit5d.io import read_zarr

        im_d = _vol3(read_zarr(maps.paths[maps.index(tb.tile)], input_bbox=pb, num_threads=nt), pb)
        if n == 2:
            _put(tim_d, j, rel, _resize(im_d, dbsz_j, wd, nt), None, nt)  # type: ignore[arg-type]
        else:
            if im_d.size == 0:
                raise StitchBlockError(
                    "Unable to perform assignment because the size of the left side is 1-by-2 "
                    "and the size of the right side is 0-by-0.",
                    "MATLAB:subsassigndimmismatch",
                )
            cells[j] = im_d
            # max and min of a single map (NaN left out), as doubles
            ranges[j] = (
                float(np.fmax.reduce(im_d, axis=None)),
                float(np.fmin.reduce(im_d, axis=None)),
            )
    minor = np.zeros(n, dtype=bool)
    if pooled_many:
        tim_d = np.zeros((n, *dbsz), dtype=np.float32)
        full_rows = np.all(dbsz_mat == np.asarray(dbsz), axis=1)
        if full_rows.any():
            top = float(np.fmax.reduce(ranges[full_rows, 1]))
            minor = top > _DMT * ranges[:, 0]
        # (with one major tile, no map is needed)
        if np.sum(~minor) > 1:
            for j in range(n):
                if minor[j]:
                    continue
                db = _cut(boxes[j], c)
                _put(tim_d, j, _rel(db, c), _resize(cells[j], dbsz_mat[j], wd, nt), None, nt)
    return tim_d, cells, ranges, dbsz_mat, minor


def _minor_tiles(
    st: dict,
    tiles: list,
    boxes: list[Box],
    tile_fns: Sequence[str],
    dtype: str,
    wd: Any,
    psz: tuple[int, int, int],
    nt: int | None,
) -> None:
    """MATLAB's large-file feather step for more than two tiles with minor ones (processStitch
    Block.m lines 342-444): read the minor tiles where the major tiles leave voxels empty or
    don't cover them, with weights; one major tile that covers the others becomes the batch
    (st['num_tiles'] = 1). Updates the batch state ``st`` in place."""
    from petakit5d.kernels.stitch import check_major_tile_valid
    from petakit5d.stitch.geometry import check_major_tile_cover

    minor = st["minor"]
    c, full = st["c"], st["full"]
    tim_f, tim_d, tim_block = st["tim_f"], st["tim_d"], st["tim_block"]
    cells, ranges = st["cells"], st["ranges"]
    n = len(tiles)
    valid, unvalid = check_major_tile_valid(tim_f, ~minor, num_threads=nt)
    # (in the batch's coordinates: tim_f is the overlap region)
    if not valid:
        unvalid = tuple(int(unvalid[q]) + c[q % 3] for q in range(6))  # type: ignore[index]
    cover, uncover = check_major_tile_cover(np.asarray(boxes, dtype=np.int64), ~minor)
    orig_major = ~minor
    if not valid or not bool(np.all(cover)):
        for j in range(n):
            if orig_major[j] and np.sum(orig_major) > 1 and cover[j]:
                continue
            u: list[int] | None
            if not cover[j]:
                uc = [int(v) for v in uncover[j]]
                u = (
                    uc
                    if valid
                    else [min(unvalid[q], uc[q]) for q in range(3)]  # type: ignore[index]
                    + [max(unvalid[q + 3], uc[q + 3]) for q in range(3)]  # type: ignore[index]
                )
            elif not valid:
                u = list(unvalid)  # type: ignore[arg-type]
            elif orig_major[j]:
                # (the major tile's weight, for the minor tiles it does not cover)
                u = list(boxes[j])
            else:
                u = None
            if u is None:
                continue
            b = boxes[j]
            t = tiles[j].tile_box
            if any(b[q] >= u[q + 3] or b[q + 3] <= u[q] for q in range(3)):
                continue
            u = _cut(u, b)
            ubb = [t[q] - b[q] + u[q] for q in range(3)] + [
                t[q] - b[q] + u[q + 3] for q in range(3)
            ]
            db = _cut(b, c)
            dbsz_j = [db[q + 3] - db[q] for q in range(3)]
            udb = _cut(u, c)
            udbsz = [udb[q + 3] - udb[q] for q in range(3)]
            if minor[j]:
                block_j = _as_class(_read(tile_fns[tiles[j].tile], ubb, None, nt), dtype)
                if not full:
                    from petakit5d.kernels.crop import indexing_3d

                    indexing_3d(tim_block, u, block_j, num_threads=nt)
                    tim_f = _put(tim_f, j, _rel(udb, c), block_j, _rel(udb, u), nt)
                else:
                    tim_f = _put(tim_f, j, u, block_j, None, nt)
                if np.sum(minor) == 1:
                    w = np.float32(_weight_power(ranges[j, 0], wd))
                    im_d = np.full([max(0, v) for v in udbsz], w, dtype=np.float32)
                    tim_d = _put(tim_d, j, _rel(udb, c), im_d, None, nt)
                    minor[j] = False
                    continue
            if orig_major[j]:
                w = np.float32(_weight_power(ranges[j, 0], wd))
                im_d = np.full([max(0, v) for v in dbsz_j], w, dtype=np.float32)
                tim_d = _put(tim_d, j, _rel(db, c), im_d, None, nt)
            else:
                im_d = cells[j]
                ps = [max(1.0, _mround((udb[q] - db[q] + 1) / psz[q])) for q in range(3)]
                pe = [
                    min(float(im_d.shape[q]), ps[q] + math.ceil(udbsz[q] / psz[q]) - 1)
                    for q in range(3)
                ]
                if udb[3] == udb[0] + 1:
                    pe[0] = ps[0]
                pbox = [int(v) - 1 for v in ps] + [int(v) for v in pe]
                im_d = _resize(_crop3(im_d, pbox, nt), udbsz, wd, nt)
                tim_d = _put(tim_d, j, _rel(udb, c), im_d, None, nt)
            minor[j] = False
    # one major tile that covers the others: the batch is that tile
    if np.sum(~minor) == 1 and bool(np.all(cover)):
        if not full:
            tim_f = tim_block[np.newaxis]
        else:
            k = int(np.nonzero(~minor)[0][0])
            tim_f = tim_f[k : k + 1]
        st["num_tiles"] = 1
    st.update(tim_f=tim_f, tim_d=tim_d, minor=minor)


def _blend_batch(
    tiles: list,
    bsz: list[int],
    blend_method: str,
    dtype: str,
    tile_fns: Sequence[str],
    maps: _Maps | None,
    psz: tuple[int, int, int] | None,
    pooled: bool,
    wd: Any,
    tile_buf: Any,
    nt: int | None,
) -> np.ndarray:
    """the blended batch (``(z, y, x)`` of bsz) of tiles (MATLAB's lines 152-482)"""
    from petakit5d.kernels import stitch as K

    n = len(tiles)
    ndt = _ndtype(dtype)
    boxes: list[Box] = []
    for tb in tiles:
        if tb.batch_box is None or tb.tile_box is None:
            raise StitchBlockError(
                "Unable to perform assignment because the size of the left side is 1-by-6 and "
                "the size of the right side is 1-by-0.",
                "MATLAB:subsassigndimmismatch",
            )
        boxes.append(tb.batch_box)
    feather = blend_method.lower() == "feather"
    full_tile = n == 1 and [boxes[0][q + 3] - boxes[0][q] for q in range(3)] == list(bsz)
    mboxes: list[list[int]] = []
    st: dict[str, Any] = {
        "num_tiles": n,
        "full": False,
        "c": [0, 0, 0, *bsz],
        "dbsz": list(bsz),
        "tim_d": None,
        "tim_block": None,
        "minor": np.zeros(n, dtype=bool),
        "cells": [None] * n,
        "ranges": np.zeros((n, 2)),
    }
    if feather:
        if n > 1:
            c, full, dbsz = _overlap_region(boxes, bsz)
            st.update(c=c, full=full, dbsz=dbsz)
            tim_d, cells, ranges, _, minor = _feather_weights(
                tiles,
                boxes,
                c,
                dbsz,
                maps,
                psz,
                wd,
                nt,  # type: ignore[arg-type]
            )
            st.update(tim_d=tim_d, cells=cells, ranges=ranges, minor=minor)
        st["tim_f"] = None if full_tile else np.zeros((n, *st["dbsz"]), dtype=ndt)
        if n > 1 and not st["full"]:
            st["tim_block"] = np.zeros(bsz, dtype=ndt)
    else:
        st["tim_f"] = None if full_tile else np.zeros((n, *bsz), dtype=ndt)

    c = st["c"]
    for j, tb in enumerate(tiles):
        if feather and n > 2 and pooled and st["minor"][j]:
            continue
        block_j = _as_class(_read(tile_fns[tb.tile], tb.tile_box, tile_buf, nt), dtype)
        b = boxes[j]
        if feather:
            if n > 1:
                db = _cut(b, c)
                if not st["full"]:
                    from petakit5d.kernels.crop import indexing_3d

                    indexing_3d(st["tim_block"], b, block_j, num_threads=nt)
                    st["tim_f"] = _put(st["tim_f"], j, _rel(db, c), block_j, _rel(db, b), nt)
                else:
                    st["tim_f"] = _put(st["tim_f"], j, b, block_j, None, nt)
                continue
        else:
            # the regions of this tile that other tiles fill, in the batch
            t = tb.tile_box
            for r in tb.regions:
                if r.tile_box is None:
                    raise StitchBlockError(
                        "Arrays have incompatible sizes for this operation.",
                        "MATLAB:sizeDimensionsMustMatch",
                    )
                m = _rel(r.tile_box, t)
                nonempty = [m[q + 3] > m[q] for q in range(3)]
                if any(nonempty[q] and (m[q] < 0 or m[q + 3] > block_j.shape[q]) for q in range(3)):
                    mm = [m[1] + 1, m[2] + 1, m[0] + 1, m[4], m[5], m[3]]
                    sz = [block_j.shape[1], block_j.shape[2], block_j.shape[0]]
                    raise StitchBlockError(
                        f"The overlap region [{' '.join(map(str, mm))}] is outside the region "
                        f"[{' '.join(map(str, sz))}] of tile {tb.tile + 1}."
                    )
                if all(nonempty):
                    mboxes.append([j] + [m[q] + b[q % 3] for q in range(6)])
        if full_tile:
            st["tim_f"] = block_j[np.newaxis]
        else:
            st["tim_f"] = _put(st["tim_f"], j, b, block_j, None, nt)

    if feather and n > 2 and pooled and st["minor"].any():
        _minor_tiles(st, tiles, boxes, tile_fns, dtype, wd, psz, nt)  # type: ignore[arg-type]

    if st["num_tiles"] == 1:
        return st["tim_f"][0]
    if blend_method.lower() == "feather":
        if st["full"]:
            nv, mex = K.feather_blending_3d(st["tim_f"], st["tim_d"], num_threads=nt)
        else:
            nv, mex = K.feather_blending_3d(
                st["tim_f"], st["tim_d"], st["tim_block"], c, num_threads=nt
            )
        if not mex and nv.dtype.kind == "f":
            nv = K.replace_nan_with_value(nv, 0, inplace=True, num_threads=nt)
        return nv
    rows = np.asarray(mboxes, dtype=np.float64).reshape(-1, 7)
    return K.tile_blending_3d(st["tim_f"], rows, blend_method, num_threads=nt)


# ---------------------------------------------------------------------------------------------
# the blending task
# ---------------------------------------------------------------------------------------------


def save_block_flag(
    flag_fullname: str, batch_inds: Sequence[int], mip_boxes: Sequence[Sequence[int]], mips: list
) -> None:
    """the flag of a blending task with the z MIPs of its batches, as processStitchBlock saves
    it: ``batchInds`` (1-based), ``mipBboxes`` (MATLAB's ``[ymin xmin ymax xmax]``, unique rows
    sorted) and ``mips`` (a column cell: per row the max of the batches' MIPs there, in batch
    order), uncompressed (MAT v6 below 2 GB, else v7.3), under a temporary name and moved into
    place. ``mip_boxes`` are our ``(y0, x0, y1, x1)``, ``mips`` ``(y, x)`` images."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import write_mat
    from petakit5d.matlabcompat.reduce import max as mmax

    rows = np.asarray(
        [[b[0] + 1, b[1] + 1, b[2], b[3]] for b in mip_boxes], dtype=np.float64
    ).reshape(-1, 4)
    uniq, ic = np.unique(rows, axis=0, return_inverse=True)
    ic = np.asarray(ic).reshape(-1)
    merged = []
    for x in range(uniq.shape[0]):
        group = [mips[k] for k in range(len(mips)) if ic[k] == x]
        if len(group) == 1:
            merged.append(np.asarray(group[0]))
        else:
            merged.append(np.asarray(mmax(np.stack(group), axis=0)[0]))
    cell = np.empty((len(merged), 1), dtype=object)
    for k, m in enumerate(merged):
        cell[k, 0] = m
    nbytes = sum(int(m.nbytes) + _CELL_HEADER_BYTES for m in merged)
    tmp = f"{flag_fullname[:-4]}_{get_uuid()}.mat"
    # (MATLAB's save writes the variables in the order of its workspace)
    data = {
        "batchInds": np.asarray(batch_inds, dtype=np.float64) + 1,
        "mips": cell,
        "mipBboxes": uniq,
    }
    write_mat(tmp, data, version="6" if nbytes < 2**31 else "7.3", compression=False)
    os.replace(tmp, flag_fullname)


def _tile_fns_from_info(block_info_fullname: str) -> list[str]:
    from petakit5d.io import read_mat

    fns = read_mat(block_info_fullname, ["tileFns"])["tileFns"]
    fns = list(fns.ravel(order="F")) if isinstance(fns, np.ndarray) else list(fns)
    return [str(f) for f in fns]


def _stitch_block_info_from_matlab(v: Any) -> Any:
    """MATLAB's struct (array) of the batches' tiles as {k: tiles of the task's batch k}"""
    from petakit5d.stitch.blockinfo import tile_blocks_from_matlab

    if v is None or (isinstance(v, np.ndarray) and v.size == 0):
        return None
    return dict(enumerate(tile_blocks_from_matlab(v)))


def _stitch_block_info_to_matlab(v: Any) -> Any:
    from petakit5d.stitch.blockinfo import tile_blocks_to_matlab

    if v is None:
        return None
    return tile_blocks_to_matlab([v[k] for k in sorted(v)])


@register(
    "processStitchBlock",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based) of the task.", convert="index1"),
    Param("BlockInfoFullname", "@(x) ischar(x)", "The block information file (tile paths)."),
    Param("PerBlockInfoFullname", "@(x) ischar(x)", "The task's per-block information file."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param("stitchFullname", "@(x) ischar(x)", "The stitched zarr (made before), written into."),
    Param(
        "stitchBlockInfo",
        "@(x) isempty(x) || isstruct(x)",
        "The batches' block information (without a per-block file).",
        convert=(_stitch_block_info_from_matlab, _stitch_block_info_to_matlab),
    ),
    Param("tileFns", "@(x) isempty(x) || iscell(x)", "Tile zarr files (without a block file)."),
    Param("Overwrite", "@islogical", "Redo a task whose flag exists."),
    Param("imSize", "@isnumeric", "Size (z, y, x) of the stitched volume.", convert="size_yxz"),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x).", convert="size_yxz"),
    Param("dtype", "@ischar", "Class of the stitched volume."),
    Param("BlendMethod", "@ischar", "'none', 'max', 'mean', 'median' or 'feather'."),
    Param("BorderSize", "@isnumeric", "Border of the batches (z, y, x).", convert="size_yxz"),
    Param("imdistFullpaths", "@iscell", "Distance maps (feather)."),
    Param("imdistFileIdx", "@isnumeric", "The distance map of each tile.", convert="index1"),
    Param("poolSize", "@isnumeric", "Pooling of the distance maps (z, y, x).", convert="size_yxz"),
    Param("weightDegree", "@isnumeric", "Power of the distance weights (pooled maps)."),
    Param("saveMIP", "@islogical", "Save the batches' z MIPs in the flag."),
    Param(
        "numThreads",
        None,
        "Threads of the kernels, reads and writes (None: the CPUs this process may use).",
        py="num_threads",
        python_only=True,
    ),
    task="S07",
)
def process_stitch_block(
    batch_inds: Sequence[int],
    block_info_fullname: str,
    per_block_info_fullname: str,
    flag_fullname: str,
    stitch_fullname: str,
    stitch_block_info: Mapping[int, Sequence[Any]] | None = None,
    tile_fns: Sequence[str] | None = None,
    *,
    overwrite: bool = False,
    im_size: Sequence[int] | None = None,
    batch_size: Sequence[int] | None = None,
    dtype: str = "uint16",
    blend_method: str = "feather",
    border_size: Sequence[int] | None = None,
    imdist_fullpaths: Sequence[str] = (),
    imdist_file_idx: Sequence[int] | None = None,
    pool_size: Sequence[int] | None = None,
    weight_degree: float = 10,
    save_mip: bool = False,
    num_threads: int | None = None,
) -> list[bool] | None:
    """``processStitchBlock``: blend the batches ``batch_inds`` of the stitched volume (size
    ``im_size``, batches ``batch_size``) into ``stitch_fullname`` and write the task's flag (see
    the module). The batches' tiles come from ``per_block_info_fullname`` (JSON, or ``.mat``),
    else ``stitch_block_info`` (``{k: [TileBlock, ...]}`` for the task's k-th batch, MATLAB's
    struct); the tiles' paths from
    ``block_info_fullname``, else ``tile_fns``. Returns the batches' done flags (None if
    nothing was done)."""
    from petakit5d.io import ReadBuffer, write_zarr
    from petakit5d.jobs.files import exist
    from petakit5d.stitch.blockinfo import read_per_block_info

    nt = num_threads
    inds = _ints(batch_inds, "batch_inds")
    flag_fullname = os.fspath(flag_fullname)
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        print(f"The block directory {flag_path} does not exist, skip the processing!", flush=True)
        return None
    if exist(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            print(
                f"The block files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!",
                flush=True,
            )
            return None
    info: Any = None
    if stitch_block_info is not None:
        info = [list(stitch_block_info.get(k, ())) for k in range(len(inds))]
    if per_block_info_fullname:
        ext = os.path.splitext(os.fspath(per_block_info_fullname))[1]
        if ext in (".mat", ".json"):
            info = read_per_block_info(per_block_info_fullname)
    if block_info_fullname:
        tile_fns = _tile_fns_from_info(os.fspath(block_info_fullname))
    fns = [os.fspath(f) for f in (tile_fns or [])]
    border = [0, 0, 0] if border_size is None or len(border_size) == 0 else _ints(border_size)
    ims = _ints(im_size, "im_size")
    bs = _ints(batch_size, "batch_size")
    method = str(blend_method)
    feather = method.lower() == "feather"
    pooled = pool_size is not None and len(pool_size) > 0
    psz = _pool_psz(_ints(pool_size)) if pooled else None
    maps = None
    if feather:
        maps = _Maps(
            imdist_fullpaths,
            [] if imdist_file_idx is None else _ints(imdist_file_idx),
            None if pooled else 2 * math.prod(bs),
        )
    b_sub = [-(-a // b) for a, b in zip(ims, bs, strict=True)]
    ndt = _ndtype(dtype)

    done = [False] * len(inds)
    mips: list[np.ndarray] = []
    mip_boxes: list[tuple[int, int, int, int]] = []
    tile_buf = ReadBuffer()
    for i, bi in enumerate(inds):
        print(f"Process batch {bi + 1}... ", end="", flush=True)
        t0 = time.perf_counter()
        sub = _batch_sub(bi, b_sub)
        ob0 = [sub[q] * bs[q] for q in range(3)]
        ob1 = [min(ob0[q] + bs[q], ims[q]) for q in range(3)]
        nbsz = [ob1[q] - ob0[q] for q in range(3)]
        obox = (*ob0, *ob1)
        if info is None or i >= len(info):
            raise StitchBlockError(
                f"Index exceeds the number of array elements. Index must not exceed "
                f"{0 if info is None else len(info)}.",
                "MATLAB:badsubscript",
            )
        tiles = list(info[i])
        if not tiles:
            nv = np.zeros(nbsz, dtype=ndt)
            write_zarr(nv, stitch_fullname, bbox=obox, create=False, num_threads=nt)
            if save_mip:
                mips.append(np.zeros(nbsz[1:], dtype=ndt))
                mip_boxes.append((ob0[1], ob0[2], ob1[1], ob1[2]))
            done[i] = True
            print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)
            continue
        if any(v > 0 for v in border):
            inner = [sub[q] != 0 and sub[q] != b_sub[q] - 1 for q in range(3)]
            bsz = [min(bs[q], nbsz[q]) + border[q] * (inner[q] + 1) for q in range(3)]
        else:
            bsz = [min(bs[q], nbsz[q]) for q in range(3)]
        # (one tile: of dtype; blended: of dtype, or single for 'mean', 'median' and MATLAB's
        # feather code)
        nv = _blend_batch(
            tiles, bsz, method, dtype, fns, maps, psz, pooled, weight_degree, tile_buf, nt
        )
        if any(v > 0 for v in border):
            s = [(sub[q] != 0) * border[q] for q in range(3)]
            nv = _crop3(nv, [s[q] for q in range(3)] + [s[q] + nbsz[q] for q in range(3)], nt)
        if list(nv.shape) != nbsz:
            nv = _crop3(nv, [0, 0, 0, *nbsz], nt)
        if save_mip:
            # cast(max(nv_block, [], 3), dtype): the cast keeps the order of values. (Before the
            # write, where MATLAB takes it after: the same values, and the port's kernels run
            # one after another, without the zarr library's threads between them.)
            mips.append(_as_class(_range_max_z(nv, nt), dtype))
            mip_boxes.append((ob0[1], ob0[2], ob1[1], ob1[2]))
        # (written in the zarr's class, as MATLAB's writezarr casts)
        write_zarr(
            np.ascontiguousarray(nv), stitch_fullname, bbox=obox, create=False, num_threads=nt
        )
        done[i] = True
        print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)
    tile_buf.release()

    if all(done):
        print(f"Processing of batches {inds[0] + 1} - {inds[-1] + 1} are done!", flush=True)
        if save_mip:
            save_block_flag(flag_fullname, inds, mip_boxes, mips)
        else:
            open(flag_fullname, "w").close()
    return done


# ---------------------------------------------------------------------------------------------
# the MIPs of the blending tasks
# ---------------------------------------------------------------------------------------------


def merge_block_mips(
    flag_fullpaths: Sequence[str], im_size: Sequence[int], dtype: str
) -> np.ndarray | None:
    """``stitch_merge_block_mips``: the z MIP (``(y, x)``) of a stitched volume of ``im_size``
    (``(z, y, x)``) merged from the MIPs its blending tasks saved in their flags (``max`` with
    the image so far, which starts at NaN for single and double and at the class's smallest
    value otherwise, flag by flag and box by box); None (the caller computes the MIP from the
    volume) if a flag can't be read or has no MIPs."""
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import from_matlab
    from petakit5d.io import read_mat
    from petakit5d.matlabcompat.reduce import maximum

    ims = _ints(im_size, "im_size")
    ndt = _ndtype(dtype)
    if dtype in ("single", "double"):
        mip = np.full(ims[1:3], dtypes.NAN, dtype=ndt)
    else:
        mip = np.full(ims[1:3], dtypes.intmin(dtype), dtype=ndt)
    for fn in flag_fullpaths:
        try:
            a = read_mat(os.fspath(fn))
        except Exception:
            a = {}
        if "mips" not in a or "mipBboxes" not in a:
            return None
        items = a["mips"]
        items = list(items.ravel(order="F")) if isinstance(items, np.ndarray) else list(items)
        boxes = np.asarray(a["mipBboxes"], dtype=np.float64).reshape(-1, 4)
        for b in range(len(items)):
            y0, x0, y1, x1 = (int(v) for v in boxes[b])
            m = from_matlab(np.asarray(items[b]))
            region = mip[y0 - 1 : y1, x0 - 1 : x1]
            region[...] = maximum(region, m.reshape(region.shape))
    return mip


def mip_memory(cluster_mode: bool, config_file: str = "") -> float:
    """``stitch_mip_memory``: the memory (GB) the frame measures the stitched volume against for
    its MIPs: this process's, and for a stitching job of the cluster (``cluster_mode`` in a
    Slurm job) at least the largest job of the config (``maxCPUNum`` x ``memPerCPU``, else 24 x
    20.9)."""
    from petakit5d.core.memory import get_system_memory
    from petakit5d.jobs.config import cluster_config_value

    mem = get_system_memory()
    if cluster_mode and os.environ.get("SLURM_JOB_ID"):
        job = float(cluster_config_value(config_file, "maxCPUNum", 24)) * float(
            cluster_config_value(config_file, "memPerCPU", 20.9)
        )
        mem = job if math.isnan(mem) else max(mem, job)
    return mem


def in_process(
    num_batches: float,
    batch_size: Sequence[int],
    dtype: str,
    blend_method: str,
    im_sizes,
    parse_cluster: bool,
    config_file: str = "",
) -> bool:
    """``stitch_in_process``: whether the frame stitches a volume of ``num_batches`` batches in
    its own process instead of with jobs: fewer than 30 batches whose batches (20 B per voxel of
    a batch of ``dtype``) and, with feather blending, distance maps (40 B per voxel of the
    largest tile, ``im_sizes`` ``(n, 3)``) fit this process (``master_task_fits``); without the
    cluster always (below 30 batches)."""
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.jobs.config import cluster_config_value
    from petakit5d.jobs.coordinator import master_task_fits

    if not num_batches < 30:
        return False
    mem = math.prod(float(v) for v in batch_size) * data_type_to_byte_number(dtype) / 2**30 * 20
    if blend_method.lower() == "feather":
        s = np.asarray(im_sizes, dtype=np.float64).reshape(-1, 3)
        mem = max(mem, float(np.max(np.prod(s, axis=1))) * 4 / 2**30 * 10)
    if not parse_cluster:
        return True
    return master_task_fits(
        mem,
        1,
        cluster_config_value(config_file, "memPerCPU", 20.9),
        cluster_config_value(config_file, "maxCPUNum", 24),
        1,
    )


# ---------------------------------------------------------------------------------------------
# separate tiles
# ---------------------------------------------------------------------------------------------


def process_separate_tiles(
    zarr_fullpaths: Sequence[str],
    stitch_path: str,
    st_indices,
    p_im_sz: Sequence[int],
    *,
    overwrite: bool = False,
    save_mip: bool = True,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 1,
    uuid: str = "",
    debug: bool = False,
    mcc_mode: bool = False,
    config_file: str = "",
) -> list[bool]:
    """``processStitchSeparteTiles``: each tile on its own canvas of the stitched size
    ``p_im_sz`` (``(z, y, x)``) at its place (``st_indices``: the tiles' 0-based voxel offsets in
    x, y, z), as ``<stitch_path>/<tile name>``: a crop of the tile with the canvas' box (from
    ``-offset``, padded; MATLAB's data order [2, 1, 3] is assumed) by ``tools.crop.crop_zarr``
    tasks through the job framework, with MIPs, retried twice with twice the CPUs and memory.
    Returns the tasks' done flags."""
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.io import get_image_dtype
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.tools.crop import crop_zarr

    del overwrite, debug
    if not uuid:
        uuid = get_uuid()
    fns = [os.fspath(f) for f in zarr_fullpaths]
    st = np.asarray(st_indices, dtype=np.float64).reshape(-1, 3)
    psz = _ints(p_im_sz, "p_im_sz")
    outs, tasks = [], []
    for f, fn in enumerate(fns):
        _, name, ext = fileparts(fn)
        out = f"{stitch_path}/{name}{ext}"
        outs.append(out)
        # (z, y, x) start: minus the offset in z, y, x
        start = [-int(st[f, 2]), -int(st[f, 1]), -int(st[f, 0])]
        box = start + [start[q] + psz[q] for q in range(3)]
        tasks.append(
            Task.of(
                crop_zarr,
                fn,
                out,
                box,
                pad=True,
                batch_size=[int(v) for v in batch_size],
                block_size=[int(v) for v in block_size],
                save_mip=bool(save_mip),
                parse_cluster=bool(parse_cluster),
                parse_parfor=False,
                master_compute=bool(master_compute),
                cpus_per_task=cpus_per_task,
                uuid=uuid,
                mcc_mode=bool(mcc_mode),
                config_file=config_file,
            )
        )
    dtype = get_image_dtype(fns[0])
    mem = math.prod(int(v) for v in batch_size) * data_type_to_byte_number(dtype) / 1024**3 * 10
    done = [False]
    for i in range(3):
        if not all(done):
            done = framework.run_tasks(
                fns,
                outs,
                tasks,
                cpus_per_task=cpus_per_task * 2**i,
                mem_allocate=mem * 2**i,
                max_trial_num=2,
                master_compute=bool(master_compute),
                parse_cluster=bool(parse_cluster),
                mcc_mode=mcc_mode,
                config_file=config_file,
            )
    return done
