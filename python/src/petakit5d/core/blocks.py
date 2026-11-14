"""Batch and block coordinates shared by the large-file paths.

The Python forms take and return ``(z, y, x)`` sizes and our boxes (``(z0, y0, x0, z1, y1, x1)``,
0-based, end-exclusive), in MATLAB's batch order and numbering:

- ``zarr_chunk_coordinates``: ``XR_zarrChunkCoordinatesExtraction`` (batches of whole blocks with
  borders; ``same_batch_size``; batches numbered with y fastest, then x, then z, as MATLAB)
- ``sub_volume_coordinates``: ``XR_subVolumeCoordinatesExtraction`` (overlapping sub-volumes, x
  fastest, then y, then z, as its ``meshgrid`` gives)
- ``bbox_to_blocks``: ``bboxToBlocks`` (the blocks of a stitched image a box touches)
- ``find_good_factor_number``: ``findGoodFactorNumber`` (sizes with prime factors 2, 3, 5, 7)

The ``*_matlab`` forms take and return MATLAB's values (``[y x z]`` sizes, 1-based inclusive
``[ymin xmin zmin ymax xmax zmax]`` rows) for files in MATLAB's formats; they are the reference
implementation (checked on 2,100 parameter sets against MATLAB), and the Python forms convert
around them. Inputs must be integral.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Sequence

import numpy as np

__all__ = [
    "BlockInfo",
    "MatlabBlockInfo",
    "bbox_to_blocks",
    "bbox_to_blocks_matlab",
    "find_good_factor_number",
    "prime_factors",
    "sub_volume_coordinates",
    "sub_volume_coordinates_matlab",
    "zarr_chunk_coordinates",
    "zarr_chunk_coordinates_matlab",
]


def _ints(v, n: int | None, what: str, scalar_ok: bool = True) -> np.ndarray:
    a = np.asarray(v, dtype=np.float64).ravel()
    if n is not None:
        if a.size == 1 and scalar_ok:
            a = np.repeat(a, n)
        if a.size != n:
            raise ValueError(
                f"{what} must have {n} values{' (or 1)' if scalar_ok else ''}, got {a.size}"
            )
    if not np.all(np.isfinite(a)) or np.any(a != np.floor(a)):
        raise ValueError(f"{what} must be integers, got {v!r}")
    return a.astype(np.int64)


def _ceil_div(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    # ceil(a ./ b) for positive integers, exactly
    return -(-a // b)


def zarr_chunk_coordinates_matlab(
    vol_size: Sequence[int],
    *,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] | int = (256, 256, 256),
    bbox: Sequence[int] | None = None,
    same_batch_size: bool = True,
    border_size: Sequence[int] | int = 50,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``XR_zarrChunkCoordinatesExtraction`` with MATLAB's values: ``(batch, region, local)``.

    Batches are whole blocks (``batch_size`` rounded up to a multiple of ``block_size``, both
    first limited to the region), numbered column-major (y fastest). ``batch_bboxes`` are the
    input regions with ``border_size`` (limited to the volume; with ``same_batch_size`` the
    last batch on an axis starts earlier so all have the same size), ``region_bboxes`` the
    output regions relative to ``bbox``, ``local_bboxes`` the output regions inside each batch.
    """
    vs = _ints(vol_size, 3, "volSize", scalar_ok=False)
    bs = _ints(batch_size, 3, "batchSize", scalar_ok=False)
    if np.any(bs < 1):
        raise ValueError("batchSize must be at least 1")
    blk = _ints(block_size, 3, "BlockSize")
    bd = _ints(border_size, 3, "BorderSize")
    if bbox is None or np.size(bbox) == 0:
        wd_start = np.ones(3, dtype=np.int64)
        wd_end = vs.copy()
    else:
        b = _ints(bbox, 6, "bbox", scalar_ok=False)
        wd_start, wd_end = b[:3], b[3:]
    wd_size = wd_end - wd_start + 1

    bs = np.minimum(wd_size, bs)
    blk = np.minimum(wd_size, blk)
    # whole blocks per batch
    bs = _ceil_div(bs, blk) * blk

    n_sub = _ceil_div(wd_size, bs)
    # ndgrid order: y fastest
    gy, gx, gz = np.meshgrid(
        np.arange(1, n_sub[0] + 1),
        np.arange(1, n_sub[1] + 1),
        np.arange(1, n_sub[2] + 1),
        indexing="ij",
    )
    subs = np.stack([gy.ravel(order="F"), gx.ravel(order="F"), gz.ravel(order="F")], axis=1)

    region = np.empty((subs.shape[0], 6), dtype=np.int64)
    region[:, :3] = (subs - 1) * bs + wd_start
    region[:, 3:] = np.minimum(region[:, :3] + bs - 1, wd_end)

    batch = region.copy()
    if same_batch_size:
        batch[:, :3] = batch[:, 3:] - bs + 1
    batch[:, :3] = np.maximum(1, batch[:, :3] - bd)
    batch[:, 3:] = np.minimum(vs, batch[:, 3:] + bd)

    local = np.empty_like(region)
    local[:, :3] = region[:, :3] - batch[:, :3] + 1
    local[:, 3:] = local[:, :3] + (region[:, 3:] - region[:, :3])

    region = region - np.concatenate([wd_start, wd_start]) + 1
    return batch, region, local


def prime_factors(n: int) -> list[int]:
    """MATLAB's ``factor`` for a non-negative integer (``factor(0)`` is 0, ``factor(1)`` is 1)."""
    n = int(n)
    if n < 0:
        raise ValueError("N must be a real non-negative integer")
    if n < 4:
        return [n]
    out = []
    while n % 2 == 0:
        out.append(2)
        n //= 2
    f = 3
    while f * f <= n:
        while n % f == 0:
            out.append(f)
            n //= f
        f += 2
    if n > 1:
        out.append(n)
    return out


def find_good_factor_number(
    given_num: int, direction_flag: int = 1, allow_odd_seven: bool = False
) -> int:
    """``findGoodFactorNumber``: the nearest number (stepping by ``direction_flag``) whose prime
    factors are at most 7; unless ``allow_odd_seven``, odd numbers with 7 as their largest
    factor are skipped (cuFFT failed for some)."""
    good = int(given_num)
    if good != given_num:
        raise ValueError("givenNum must be an integer")
    pf = prime_factors(good)
    while pf[-1] > 7 or (not allow_odd_seven and pf[0] > 2 and pf[-1] == 7):
        good += direction_flag
        pf = prime_factors(good)
    return good


def sub_volume_coordinates_matlab(
    vol_size: Sequence[int],
    *,
    chunk_size: Sequence[int] = (1024, 1024, 1024),
    overlap_size: Sequence[int] | int = 50,
    max_sub_volume: float = 1.2e9,
    size_range: Sequence[float] = (0.75, 1.5),
    use_good_num: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """``XR_subVolumeCoordinatesExtraction`` with MATLAB's values: ``(xmin, xmax, ymin, ymax,
    zmin, zmax, nn)``.

    Overlapping sub-volumes of about ``chunk_size`` (within ``size_range`` of it) whose voxel
    count is at most ``max_sub_volume``; sizes made good FFT sizes with ``use_good_num``. The
    order is MATLAB's ``meshgrid`` order: x fastest, then y, then z.
    """
    vs = np.asarray(vol_size, dtype=np.float64).ravel()
    if vs.size != 3:
        raise ValueError("volSize must have 3 values")
    if np.any(vs < 2):
        raise ValueError(
            f"The volume size [{' '.join(f'{v:g}' for v in vs)}] is not valid for a 3d volume."
        )
    cs = np.asarray(chunk_size, dtype=np.float64).ravel()
    ol = np.asarray(overlap_size, dtype=np.float64).ravel()
    if ol.size == 1:
        ol = np.repeat(ol, 3)
    elif ol.size == 2:
        ol = np.array([ol[0], ol[0], ol[1]])
    sr = [float(x) for x in np.asarray(size_range, dtype=np.float64).ravel()]

    ratio = vs / cs
    min_p = np.ceil(ratio / sr[1])
    max_p = np.ceil(ratio / sr[0])

    chunk_cell, ratio_cell, part_cell = [], [], []
    for i in range(3):
        pr = np.arange(min_p[i], max_p[i] + 1)  # colon of integer-valued doubles
        c = np.ceil(vs[i] / pr + ol[i])
        chunk_cell.append(c)
        ratio_cell.append(c / cs[i])
        part_cell.append(pr)

    n_parts = [c.size for c in chunk_cell]
    inds = [0, 0, 0]

    def chosen():
        return [chunk_cell[k][inds[k]] for k in range(3)]

    ch = chosen()
    while ch[0] * ch[1] * ch[2] > max_sub_volume:
        ratios = [ratio_cell[k][min(inds[k] + 1, n_parts[k] - 1)] for k in range(3)]
        # sort(..., 'descend') is stable
        order = sorted(range(3), key=lambda k: -ratios[k])
        ct = 0
        while inds[order[ct]] >= n_parts[order[ct]] - 1:
            ct += 1
            if ct == 3:
                raise IndexError(
                    "Index exceeds the number of array elements (no axis can be split further)."
                )
        inds[order[ct]] += 1
        ch = chosen()

    parts = np.array([part_cell[k][inds[k]] for k in range(3)])
    csz = np.ceil((vs + ol * (parts - 1)) / parts)
    if use_good_num:
        csz = np.array([float(find_good_factor_number(int(c), 1)) for c in csz])

    with np.errstate(divide="ignore", invalid="ignore"):
        aol = np.ceil((csz * parts - vs) / (parts - 1))
    aol[parts == 1] = 0
    res = vs - ((csz - aol) * parts + aol)
    res[parts == 1] = 0

    starts, ends = [], []
    for i in range(3):
        step = csz[i] - aol[i]
        p = int(parts[i])
        if step <= 0:
            raise IndexError("Index exceeds the number of array elements (empty start positions).")
        n_avail = int(math.floor((vs[i] - 1) / step)) + 1
        if p > n_avail:
            raise IndexError(
                "Index exceeds the number of array elements (too few start positions)."
            )
        lo = 1 + step * np.arange(p)
        r = int(res[i])
        if r > 0:
            if r > p:
                raise IndexError(
                    "Index in position 1 is invalid (shift beyond the first sub-volume)."
                )
            lo[p - r :] = lo[p - r :] + np.arange(1, r + 1)
        hi = np.minimum(vs[i], lo + csz[i] - 1)
        starts.append(lo)
        ends.append(hi)

    yl, xl, zl = starts
    yr, xr, zr = ends

    # [Y, X, Z] = meshgrid(y, x, z) has size (nx, ny, nz); (:) runs x fastest, then y, then z
    def grid(v, axis):
        shape = (xl.size, yl.size, zl.size)
        idx = np.indices(shape)[axis]
        return v[idx].ravel(order="F")

    ymin, ymax = grid(yl, 1), grid(yr, 1)
    xmin, xmax = grid(xl, 0), grid(xr, 0)
    zmin, zmax = grid(zl, 2), grid(zr, 2)
    as_int = lambda v: v.astype(np.int64)  # noqa: E731
    return (
        as_int(xmin),
        as_int(xmax),
        as_int(ymin),
        as_int(ymax),
        as_int(zmin),
        as_int(zmax),
        int(ymin.size),
    )


@dataclasses.dataclass(frozen=True)
class MatlabBlockInfo:
    """``bboxToBlocks``' struct (MATLAB field names in snake_case; all values 1-based).

    Rows of ``block_sub``, ``block_ind``, ``b_coords``, ``w_coords`` and ``bbox_coords`` are the
    touched blocks in column-major order (y fastest).
    """

    block_sub: np.ndarray  # (n, 3) block subscripts
    block_ind: np.ndarray  # (n,) linear block index (sub2ind of b_sub_size)
    b_coords: np.ndarray  # (n, 6) region inside each block (with borders)
    w_coords: np.ndarray  # (n, 6) the same region in the stitched image
    bbox_coords: np.ndarray  # (n, 6) the same region relative to bbox_start
    bbox: tuple[int, ...]
    bbox_start: tuple[int, ...]
    block_size: tuple[int, ...]
    image_size: tuple[int, ...]
    b_sub_size: tuple[int, ...]
    border_size: tuple[int, ...]


def _single_axis_split(
    s: int, t: int, bsz: int, nb: int, bdsz: int
) -> tuple[np.ndarray, np.ndarray]:
    lo = max(1, -(-(s - bdsz) // bsz))
    hi = min(nb, -(-(t + bdsz) // bsz))
    subs = np.arange(lo, hi + 1, dtype=np.int64)
    if subs.size == 0:
        raise IndexError("Index exceeds the number of array elements (the box touches no block).")
    n = subs.size
    if bdsz > 0:
        # the part of [s, t] in each block with its borders, which spans (b - 1) * bsz + 1 -
        # bdsz (none before the first block) to b * bsz + bdsz
        block_s = (subs - 1) * bsz + 1 - (subs != 1).astype(np.int64) * bdsz
        cs = np.maximum(s, block_s) - block_s + 1
        ct = np.minimum(t, subs * bsz + bdsz) - block_s + 1
    else:
        cs = np.concatenate([[max(1, s - (subs[0] - 1) * bsz)], np.ones(n - 1, dtype=np.int64)])
        ct = np.concatenate(
            [np.full(n - 1, bsz, dtype=np.int64), [min(bsz, t - (subs[-1] - 1) * bsz)]]
        )
    return subs, np.stack([cs.astype(np.int64), ct.astype(np.int64)])


def bbox_to_blocks_matlab(
    bbox: Sequence[int],
    bbox_start: Sequence[int],
    block_size: Sequence[int],
    image_size: Sequence[int],
    border_size: Sequence[int] = (0, 0, 0),
) -> MatlabBlockInfo:
    """``bboxToBlocks`` with MATLAB's values: the blocks (of ``block_size``, tiling
    ``image_size``) a box touches."""
    b = _ints(bbox, 6, "bbox", scalar_ok=False)
    st = _ints(bbox_start, 3, "bboxStart", scalar_ok=False)
    bs = _ints(block_size, 3, "blockSize", scalar_ok=False)
    im = _ints(image_size, 3, "imageSize", scalar_ok=False)
    bd = _ints(border_size, 3, "borderSize", scalar_ok=False)
    n_sub = _ceil_div(im, bs)

    subs, coords = [], []
    for i in range(3):
        s_i, c_i = _single_axis_split(
            int(b[i]), int(b[3 + i]), int(bs[i]), int(n_sub[i]), int(bd[i])
        )
        subs.append(s_i)
        coords.append(c_i)

    def nd(vs):
        g = np.meshgrid(*vs, indexing="ij")
        return np.stack([x.ravel(order="F") for x in g], axis=1)

    block_sub = nd(subs)
    b_coords = np.concatenate([nd([c[0] for c in coords]), nd([c[1] for c in coords])], axis=1)
    off = (block_sub - 1) * bs - (block_sub != 1) * bd
    w_coords = np.tile(off, (1, 2)) + b_coords
    bbox_coords = w_coords - np.tile(b[:3] - st, 2)
    block_ind = (
        (block_sub[:, 0] - 1)
        + n_sub[0] * ((block_sub[:, 1] - 1) + n_sub[1] * (block_sub[:, 2] - 1))
        + 1
    )
    tup = lambda v: tuple(int(x) for x in v)  # noqa: E731
    return MatlabBlockInfo(
        block_sub=block_sub,
        block_ind=block_ind.astype(np.int64),
        b_coords=b_coords,
        w_coords=w_coords,
        bbox_coords=bbox_coords,
        bbox=tup(b),
        bbox_start=tup(st),
        block_size=tup(bs),
        image_size=tup(im),
        b_sub_size=tup(n_sub),
        border_size=tup(bd),
    )


# -- the Python forms -------------------------------------------------------------------------


def _zyx3(v, what: str) -> list:
    """A (z, y, x) value (or a scalar for all axes) as MATLAB's [y x z]."""
    from petakit5d.core.axes import to_yxz

    a = np.asarray(v, dtype=np.float64).ravel()
    if a.size == 1:
        return [float(a[0])] * 3
    if a.size != 3:
        raise ValueError(f"{what} must have 3 values (z, y, x) or 1, got {a.size}")
    return list(to_yxz(a))


def zarr_chunk_coordinates(
    vol_size: Sequence[int],
    *,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] | int = (256, 256, 256),
    bbox=None,
    same_batch_size: bool = True,
    border_size: Sequence[int] | int = 50,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``XR_zarrChunkCoordinatesExtraction``: ``(batch_boxes, region_boxes, local_boxes)``.

    Sizes are ``(nz, ny, nx)``-ordered (``block_size`` and ``border_size`` may be scalars),
    ``bbox`` is our box (None: the whole volume). Each result is an ``(n, 6)`` int64 array of our
    boxes, one row per batch in MATLAB's order (y fastest, then x, then z): the input region with
    borders, the output region relative to ``bbox``'s start, and the output region inside the
    batch.
    """
    from petakit5d.core.bbox import box_to_matlab, boxes_from_matlab, check_box

    vs = _zyx3(vol_size, "vol_size")
    mb = None
    if bbox is not None and np.size(bbox) != 0:
        mb = list(
            box_to_matlab(
                check_box(bbox, tuple(np.asarray(vol_size, dtype=np.int64)), allow_outside=True)
            )
        )
        if any(isinstance(v, float) for v in mb):
            mb = [vs[i % 3] if v == math.inf else v for i, v in enumerate(mb)]
    b, r, loc = zarr_chunk_coordinates_matlab(
        vs,
        batch_size=_zyx3(batch_size, "batch_size"),
        block_size=_zyx3(block_size, "block_size"),
        bbox=mb,
        same_batch_size=same_batch_size,
        border_size=_zyx3(border_size, "border_size"),
    )
    return boxes_from_matlab(b), boxes_from_matlab(r), boxes_from_matlab(loc)


def sub_volume_coordinates(
    vol_size: Sequence[int],
    *,
    chunk_size: Sequence[int] = (1024, 1024, 1024),
    overlap_size: Sequence[int] | int = 50,
    max_sub_volume: float = 1.2e9,
    size_range: Sequence[float] = (0.75, 1.5),
    use_good_num: bool = True,
) -> np.ndarray:
    """``XR_subVolumeCoordinatesExtraction``: an ``(n, 6)`` int64 array of our boxes, in MATLAB's
    order (x fastest, then y, then z). Sizes and ``overlap_size`` are ``(z, y, x)`` (or scalars)."""
    from petakit5d.core.bbox import boxes_from_matlab

    xmin, xmax, ymin, ymax, zmin, zmax, _ = sub_volume_coordinates_matlab(
        _zyx3(vol_size, "vol_size"),
        chunk_size=_zyx3(chunk_size, "chunk_size"),
        overlap_size=_zyx3(overlap_size, "overlap_size"),
        max_sub_volume=max_sub_volume,
        size_range=size_range,
        use_good_num=use_good_num,
    )
    return boxes_from_matlab(np.stack([ymin, xmin, zmin, ymax, xmax, zmax], axis=1))


@dataclasses.dataclass(frozen=True)
class BlockInfo:
    """``bboxToBlocks``' struct in our forms: block subscripts and boxes ``(z, y, x)``, 0-based,
    end-exclusive; rows in MATLAB's order (block y fastest, then x, then z).
    """

    block_sub: np.ndarray  # (n, 3) block subscripts (z, y, x), 0-based
    block_ind: np.ndarray  # (n,) MATLAB's linear block number minus 1 (y fastest)
    b_coords: np.ndarray  # (n, 6) region inside each block (with borders)
    w_coords: np.ndarray  # (n, 6) the same region in the stitched image
    bbox_coords: np.ndarray  # (n, 6) the same region relative to bbox_start
    bbox: tuple[int, ...]
    bbox_start: tuple[int, ...]  # (z, y, x), 0-based
    block_size: tuple[int, ...]  # (z, y, x)
    image_size: tuple[int, ...]  # (nz, ny, nx)
    b_sub_size: tuple[int, ...]  # blocks per axis (z, y, x)
    border_size: tuple[int, ...]  # (z, y, x)


def bbox_to_blocks(
    bbox,
    bbox_start: Sequence[int],
    block_size: Sequence[int],
    image_size: Sequence[int],
    border_size: Sequence[int] = (0, 0, 0),
) -> BlockInfo:
    """``bboxToBlocks``: the blocks (of ``block_size``, tiling ``image_size``) that our box
    ``bbox`` touches; ``bbox_start`` is the box's start in the original tile (0-based, z, y, x)."""
    from petakit5d.core.axes import from_yxz
    from petakit5d.core.bbox import box_to_matlab, boxes_from_matlab, check_box

    b = check_box(bbox, allow_outside=True)
    if len(b) != 6 or None in b:
        raise ValueError("bbox_to_blocks takes a closed 3D box")
    start = [int(v) + 1 for v in _zyx3(bbox_start, "bbox_start")]
    m = bbox_to_blocks_matlab(
        box_to_matlab(b),
        start,
        _zyx3(block_size, "block_size"),
        _zyx3(image_size, "image_size"),
        _zyx3(border_size, "border_size"),
    )
    tup = lambda v: tuple(int(x) for x in from_yxz(np.asarray(v)))  # noqa: E731
    return BlockInfo(
        block_sub=m.block_sub[:, [2, 0, 1]] - 1,
        block_ind=m.block_ind - 1,
        b_coords=boxes_from_matlab(m.b_coords),
        w_coords=boxes_from_matlab(m.w_coords),
        bbox_coords=boxes_from_matlab(m.bbox_coords),
        bbox=tuple(b),
        bbox_start=tuple(int(v) - 1 for v in from_yxz(np.asarray(start))),
        block_size=tup(m.block_size),
        image_size=tup(m.image_size),
        b_sub_size=tup(m.b_sub_size),
        border_size=tup(m.border_size),
    )
