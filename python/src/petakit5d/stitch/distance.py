"""Stitching's distance maps: the feather weights of the tiles (``compute_tile_bwdist``,
``compute_tile_bwdist_mip_slabs``, ``compute_tile_distance_transform``,
``distance_weight_single_axis``). The kernels that resize and blend with the maps
(``feather_distance_map_resize_3d``, ``feather_blending_3d``) are K04's
(``petakit5d.kernels.stitch``), re-exported here.

A tile's map, as MATLAB computes it (``tile_distance_map``):

- the empty voxels: the tile's zeros, its first and last rows and columns, and, for a volume,
  its first and last slices;
- per z slice the 2D Euclidean distance to the nearest empty voxel (``bwdist``, single), as
  ``fastPower((d + 1) / 10, weightDegree)``; with ``single_dist_map`` (the default) the centre
  slice's map serves every slice whose empty voxels are the same, else (``overlaps``) the map is
  1 outside the tile's overlap regions and, in each region in turn, the distance within the
  region padded by one voxel (regions whose padded crop has no empty voxel stay 1);
- times a z window, ``fastPower`` of the middle ``nz`` values of a Tukey window of
  ``max(round(1.1 nz), nz + 2)`` points (``tukeywin(n, 0.5)``; a slice-less tile: 1), in double
  and rounded to single, times 0 on the tile's zeros;
- with ``dist_bbox``, times ``max(fastPower(wy .* wx .* wz, weightDegree), 1e-45)`` (single):
  1 inside the box, Hann tapers of 100 voxels beyond its ends and a decay of 0.99 per voxel
  further out (``distance_weight_single_axis``).

The large-file variant (``tile_distance_map_mip_slab``, ``compute_tile_bwdist_mip_slabs``)
computes the map of a tile's z MIP slab (``MIP_Slabs_pooling_<pool>/<name>_z.zarr``) without the
powers: ``(d + 1) / 10`` times the Tukey window, times ``max(wy .* wx .* wz, 1e-50^(1/wd))``
with the box, tapers and decay scaled by the pool sizes; the blending resizes it to each tile
(``feather_distance_map_resize_3d``) and raises it to the power there.

Every voxel gets MATLAB's operations in MATLAB's order and class: the distances and powers in
single, the window in double and the single-with-double product rounded to single; the slices
are independent, so they are computed in threads (results do not depend on the thread count).

Arrays are ``(z, y, x)`` C order, boxes ``(z0, y0, x0, z1, y1, x1)`` 0-based and end-exclusive,
sizes and block sizes ``(z, y, x)``. ``data_order_axes`` is the array axis of x, y and z of the
tile (``stitch.geometry.parse_data_order``; MATLAB's ``dataOrderMat``, ``[2 1 3]`` being
``(2, 1, 0)``). Overlap regions are the block information's (its ``ol_region_cell``, rows x, y
and z), over the tile's voxels. Files keep MATLAB's names (``imdist/<first
tile>/<tile>_wr_<wd>.zarr``, ``..._z_pooling_<pool>_wr_1_for_<wd>.zarr``, with MATLAB's ``//``
after the folder) and its ``parameters.mat``; the maps are native zarr arrays. On Windows the
folder is ``imdist/`` itself, as MATLAB's ``ispc`` branch has it (``imdist//<tile>_wr_<wd>.zarr``),
so the stitched volumes of one result folder share it there.

The overlap regions (rows x, y, z) are read in the order of the tile's axes (``x, y, z`` of the
array turned to put z last, else the tile's own order); the large-file map of a tile whose z is
not its last axis is turned along z's dimension and back, and takes its size after turning back;
without location ids the tiles are of one location, and one distance box is for every
per-tile map.

MATLAB quirks kept (each with a test): ``compute_tile_bwdist`` ignores its ``compressor`` and
``shard_size`` (createzarr's defaults); the large-file variant picks the chunk-file separator
from ``sz / blockSize``, a matrix division (one number: the least-squares quotient), not
``sz ./ blockSize``; an output that exists when a task starts is left as it is (``overwrite``
only replaces one that appears while it runs); ``distance_weight_single_axis`` returns double
ones when the box covers the axis and grows its result when the box starts beyond it; on Windows
a call whose overlap regions differ from the shared folder's ``parameters.mat`` removes the
whole ``imdist/``, other volumes' maps included.
"""

from __future__ import annotations

import math
import os
import shutil
import sys
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from petakit5d.kernels.stitch import feather_blending_3d, feather_distance_map_resize_3d
from petakit5d.params import Param, register
from petakit5d.stitch.geometry import data_order_axes_from_matlab, data_order_axes_to_matlab

__all__ = [
    "DistanceError",
    "compute_tile_bwdist",
    "compute_tile_bwdist_mip_slabs",
    "compute_tile_distance_transform",
    "distance_weight_single_axis",
    "feather_blending_3d",
    "feather_distance_map_resize_3d",
    "mip_slab_name",
    "tile_distance_map",
    "tile_distance_map_mip_slab",
    "z_window",
]

# compute_tile_bwdist's box weights: taper length, decay per voxel, window
DIST_BUFFER_SIZE = 100
DIST_DFACTOR = 0.99
DIST_WIN_TYPE = "hann"
# max(im_dist_wt, 1e-45): single(1e-45) is the smallest subnormal
_TINY = np.float32(1e-45)
# the GB a task asks for: this many times the first tile's voxels in single (MATLAB's * 10)
TASK_MEMORY_FACTOR = 10


class DistanceError(ValueError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------

# Python axis of MATLAB's dimensions 1 (y), 2 (x), 3 (z)
_AX = {1: 1, 2: 2, 3: 0}


def _threads(n: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(n)


def _parallel(nt: int, fn: Callable[[int], Any], items: Sequence[int]) -> None:
    """fn(i) for every item, on ``nt`` threads (each item's work releases the GIL)"""
    items = list(items)
    if nt <= 1 or len(items) <= 1:
        for i in items:
            fn(i)
        return
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=min(nt, len(items))) as ex:
        for f in [ex.submit(fn, i) for i in items]:
            f.result()


def _mround(x: float) -> float:
    from petakit5d.core.dtypes import mround

    return float(mround(float(x)))


def _int(v: Any, what: str) -> int:
    f = float(v)
    if not math.isfinite(f) or f != math.floor(f):
        raise DistanceError(f"{what} must be an integer, got {v!r}")
    return int(f)


def _mpermute_axes(p: Sequence[int]) -> tuple[int, int, int]:
    """the np.transpose axes of MATLAB's ``permute(a, p)`` of a 3D image array"""
    q = [0, 0, 0]
    for k in range(3):
        q[_AX[k + 1]] = _AX[int(p[k])]
    return tuple(q)  # type: ignore[return-value]


def _mpermute(vol: np.ndarray, p: Sequence[int]) -> np.ndarray:
    q = _mpermute_axes(p)
    if q == (0, 1, 2):
        return vol
    return np.ascontiguousarray(np.transpose(vol, q))


def _as_volume(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a)
    if a.ndim == 2:
        return a.reshape((1, *a.shape))
    if a.ndim != 3:
        raise ValueError(f"a tile is an image (y, x) or a volume (z, y, x), not {a.ndim}D")
    return a


def _dom(data_order_axes: Sequence[int]) -> list[int]:
    """MATLAB's dataOrderMat of data_order_axes"""
    dom = [int(v) for v in data_order_axes_to_matlab(data_order_axes)]
    if sorted(dom) != [1, 2, 3]:
        raise ValueError(
            f"data_order_axes must be a permutation of (0, 1, 2), got {data_order_axes}"
        )
    return dom


def _fast_power_single(b: np.ndarray, e: int) -> np.ndarray | None:
    """PetaKit5D's fastPower of a single array: ``1 ./ b`` first for a negative exponent, the
    base for exponent 1, ``result .* base`` for each set bit (the double scalar 1 first, so the
    first product is the base itself) and ``base .* base`` per bit; None for exponent 0 (the
    double scalar 1). Overflow gives Inf and underflow 0, without warnings, as in MATLAB."""
    with np.errstate(over="ignore", under="ignore", divide="ignore"):
        if e < 0:
            b = np.divide(np.float32(1), b)
            e = -e
        if e == 1:
            return b
        result = None
        while e > 0:
            if e & 1:
                result = b.copy() if result is None else np.multiply(result, b, out=result)
            e >>= 1
            if e:
                b = np.multiply(b, b)
    return result


def _weight_degree(wd: Any) -> int:
    e = float(wd)
    if e != math.floor(e) or not math.isfinite(e):
        # (MATLAB: bitget of a non-integer in fastPower)
        raise DistanceError(
            "Double inputs must have integer values in the range of ASSUMEDTYPE.",
            "MATLAB:bitget:invalidDoubleType",
        )
    return int(e)


def z_window(nz: int, weight_degree: float | None = None) -> np.ndarray:
    """The z window of a tile of ``nz`` slices: the middle values (from ``max(round(0.05 nz),
    2)``) of ``tukeywin(max(round(1.1 nz), nz + 2), 0.5)``, raised to ``weight_degree`` by
    fastPower in double when it is given; ones for one slice (MATLAB's ``win_z = 1``)."""
    from petakit5d.matlabcompat.elementary import fast_power, tukeywin

    nz = int(nz)
    if nz <= 1:
        w = np.ones(max(nz, 1))
    else:
        nw = int(max(_mround(nz * 1.1), nz + 2))
        ws = int(max(_mround(nz * 0.05), 2))
        w = tukeywin(nw, 0.5)[ws - 1 : ws - 1 + nz]
    if weight_degree is None:
        return np.ascontiguousarray(w, dtype=np.float64)
    r = fast_power(w, _weight_degree(weight_degree))
    return np.broadcast_to(np.asarray(r, dtype=np.float64), w.shape).copy()


def _empty(
    vol: np.ndarray,
    z: int,
    three_d: bool,
    keep: np.ndarray | None = None,
    ys: slice = slice(None),
    xs: slice = slice(None),
) -> np.ndarray:
    """the empty voxels of slice z, or of its crop [ys, xs]: the tile's zeros (``keep``, the
    slice's data, if given), its first and last rows and columns, and in a volume its whole
    first and last slices"""
    nz, ny, nx = vol.shape
    y0, y1, _ = ys.indices(ny)
    x0, x1, _ = xs.indices(nx)
    if three_d and (z == 0 or z == nz - 1):
        return np.ones((max(y1 - y0, 0), max(x1 - x0, 0)), dtype=bool)
    e = np.logical_not(keep) if keep is not None else np.equal(vol[z, ys, xs], 0)
    if e.size:
        if y0 == 0:
            e[0, :] = True
        if y1 == ny:
            e[-1, :] = True
        if x0 == 0:
            e[:, 0] = True
        if x1 == nx:
            e[:, -1] = True
    return e


def _distance(e: np.ndarray, wd: int | None) -> np.ndarray | None:
    """``fastPower((bwdist(e) + 1) / 10, wd)`` of one slice (single; ``wd`` None: no power;
    None for wd 0, the double scalar 1)"""
    from petakit5d.matlabcompat.distance import bwdist

    d = bwdist(np.ascontiguousarray(e), threads=1)
    # single + double 1 and single / double 10 round as single arithmetic does
    np.add(d, np.float32(1), out=d)
    np.divide(d, np.float32(10), out=d)
    if wd is None:
        return d
    return _fast_power_single(d, wd)


class _BoxWeights:
    """the box weights ``max(f(wy .* wx' .* wz), floor)`` of each slice: f is fastPower with the
    weight degree (None: none); slices with the same wz share them"""

    def __init__(self, wy, wx, wz, power: int | None, floor: np.float32):
        self.p = np.multiply.outer(np.asarray(wy, np.float32), np.asarray(wx, np.float32))
        self.wz = np.asarray(wz, np.float32)
        self.power = power
        self.floor = floor
        self._cache: dict[int, np.ndarray] = {}

    def of(self, z: int) -> np.ndarray:
        key = int(self.wz[z].view(np.uint32))
        w = self._cache.get(key)
        if w is None:
            w = np.multiply(self.p, self.wz[z])
            if self.power is not None:
                w = _fast_power_single(w, self.power)
                if w is None:  # (weight degree 0: the double scalar 1)
                    w = np.ones_like(self.p)
            w = np.fmax(w, self.floor)
            self._cache[key] = w
        return w


def _finish_slice(
    out_z: np.ndarray, src: np.ndarray | None, w: float, keep: np.ndarray, box: np.ndarray | None
) -> None:
    """``out_z = single(src * w) .* keep`` (src single, w double; src None: the double scalar
    1), then ``.* box``; Inf .* 0 is MATLAB's NaN on every processor (the kernel's rule)"""
    from petakit5d.kernels.ops import scale_offset

    if src is None:
        out_z[...] = keep
    else:
        scale_offset(src, scale=float(w), mask=keep, out=out_z, num_threads=1)
    if box is not None:
        np.multiply(out_z, box, out=out_z)


def _check_out(out: np.ndarray | None, shape: tuple[int, ...]) -> np.ndarray:
    if out is None:
        return np.empty(shape, dtype=np.float32)
    if out.dtype != np.float32 or out.shape != shape or not out.flags.c_contiguous:
        raise ValueError(f"out must be a C-contiguous single array of shape {shape}")
    if not out.flags.writeable:
        raise ValueError("out must be writeable")
    return out


# ---------------------------------------------------------------------------------------------
# distance_weight_single_axis
# ---------------------------------------------------------------------------------------------


def _grow(w: np.ndarray, n: int) -> np.ndarray:
    """MATLAB's growth of a vector by an assignment beyond its end (zeros between)"""
    if n <= w.size:
        return w
    return np.concatenate([w, np.zeros(n - w.size)])


def _dwsa(
    sz: int, e1: float, e2: float, buffer_size: float, dfactor: float, win_type: str
) -> np.ndarray:
    """distance_weight_single_axis with MATLAB's values ([s t], 1-based, inclusive)"""
    from petakit5d.matlabcompat import elementary as el

    w = np.ones(sz)
    s = max(1.0, float(e1))
    t = min(float(sz), float(e2))
    if s == 1 and t == sz:
        return w
    if win_type.lower() != "hann":
        raise DistanceError(f"Unsupported window type {win_type}")
    s = _int(s, "the box start")
    t = _int(t, "the box end")
    b = _int(buffer_size, "bufferSize")
    if t < 1:
        raise DistanceError(
            "Index in position 1 is invalid. Array indices must be positive integers or logical "
            "values.",
            "MATLAB:badsubscript",
        )
    w = _grow(w, s)
    w[:s] = 0
    w[t - 1 :] = 0
    win_dist = np.arange(b + 1, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        win = el.power(el.cos(((0.5 * math.pi) * (b - win_dist)) / b), 2)
    win = np.fmax(win, 1e-3)
    s0 = max(1, s - b)
    w[s0 - 1 : s] = win[b - (s - max(1, s - b)) :]
    pos = dfactor > 0
    lb, ub = 1e-4 * pos, 1e-3 * pos
    if s0 > 1:
        k = np.arange(s0 - 1, 0, -1, dtype=np.float64)
        w[: s0 - 1] = el.power(dfactor, k) * min(max(lb, w[s0 - 1]), ub)
    win_f = win[::-1]
    t1 = min(sz, t + b)
    w[t - 1 : t1] = win_f[: min(sz, t + b) - t + 1]
    if t1 < sz:
        k = np.arange(1, sz - t1 + 1, dtype=np.float64)
        w[t1:sz] = el.power(dfactor, k) * min(max(lb, w[t1 - 1]), ub)
    if dfactor > 0:
        w = np.fmax(w, 1e-30)
    return w.astype(np.float32)


def distance_weight_single_axis(
    n: int,
    end_points: Sequence[int | None],
    buffer_size: float = 10,
    dfactor: float = 0.99,
    win_type: str = "hann",
) -> np.ndarray:
    """``distance_weight_single_axis(sz, endPoints, bufferSize, dfactor, winType)``: weights of
    the ``n`` voxels of an axis for the range ``end_points`` (``(start, stop)``, 0-based,
    end-exclusive; ``None``: open): 1 on the range, Hann tapers ``cos(pi/2 (b - k)/b)^2`` (at
    least 1e-3) from its first and last voxel over ``buffer_size`` voxels outward, then a decay
    ``dfactor^k`` times the taper's last value (kept in [1e-4, 1e-3]); at least 1e-30.
    Single; double ones when the range covers the axis (as MATLAB returns them). A range that
    starts beyond the axis gives a longer vector, as MATLAB's assignment grows it."""
    a, b = (list(end_points) + [None, None])[:2]
    e1 = -math.inf if a is None else float(a) + 1
    e2 = math.inf if b is None else float(b)
    return _dwsa(int(n), e1, e2, buffer_size, float(dfactor), win_type)


def _box_axis_weights(
    sz_m: Sequence[int], box_m: Sequence[float], buffers: Sequence[float], dfactors: Sequence[float]
) -> list[np.ndarray]:
    """the weights along MATLAB's dimensions 1, 2, 3 of a box ``[ymin xmin zmin ymax xmax
    zmax]``"""
    return [
        _dwsa(int(sz_m[k]), box_m[k], box_m[k + 3], buffers[k], dfactors[k], DIST_WIN_TYPE)
        for k in range(3)
    ]


def _check_weights(ws: list[np.ndarray], shape: tuple[int, int, int]) -> None:
    """MATLAB's implicit expansion of im_dist (Python shape) with wy .* wx' .* wz"""
    want = (shape[1], shape[2], shape[0])
    for k in range(3):
        if ws[k].size not in (1, want[k]):
            raise DistanceError(
                "Arrays have incompatible sizes for this operation.",
                "MATLAB:sizeDimensionsMustMatch",
            )


def _expand(ws: list[np.ndarray], shape: tuple[int, int, int]) -> list[np.ndarray]:
    """wy, wx, wz of length 1 broadcast to the array's sizes"""
    want = (shape[1], shape[2], shape[0])
    return [np.broadcast_to(w, (want[k],)) if w.size == 1 else w for k, w in enumerate(ws)]


# ---------------------------------------------------------------------------------------------
# a tile's map, in memory
# ---------------------------------------------------------------------------------------------


def tile_distance_map(
    tile: np.ndarray,
    weight_degree: float = 10,
    *,
    single_dist_map: bool = True,
    overlaps: Sequence[Sequence[int]] = (),
    dist_bbox: Sequence[int | None] | None = None,
    data_order_axes: Sequence[int] = (2, 1, 0),
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """The feather weights of a tile (``(z, y, x)`` or ``(y, x)``, any numeric class; nonzero is
    data), single, as ``compute_tile_bwdist`` computes them (module docstring). ``overlaps``:
    the tile's overlap regions in order (boxes over its voxels in the turned array, see the
    module docstring), used without ``single_dist_map``. ``dist_bbox``: the box of full weight
    (``None`` ends are open). ``out``: a single array of the result's shape to write into
    (returned). An image's result is ``(1, ny, nx)``."""
    wd = _weight_degree(weight_degree)
    dom = _dom(data_order_axes)
    vol = _as_volume(tile)
    if dom[2] != 3:
        vol = _mpermute(vol, dom)
    nt = _threads(num_threads)
    shape = vol.shape
    rev = None
    if dom[2] != 3:
        rev = [int(i) + 1 for i in np.argsort(dom, kind="stable")]
        q = _mpermute_axes(rev)
        shape = tuple(shape[q[k]] for k in range(3))
    res = _check_out(out, shape)
    m = res if rev is None else np.empty(vol.shape, dtype=np.float32)
    ws = None
    if dist_bbox is not None and len(dist_bbox):
        from petakit5d.core.bbox import box_to_matlab

        box_m = [float(v) for v in box_to_matlab(dist_bbox)]
        sz_m = (shape[1], shape[2], shape[0])
        n = DIST_BUFFER_SIZE
        ws = _box_axis_weights(sz_m, box_m, (n, n, n), (DIST_DFACTOR,) * 3)
        _check_weights(ws, shape)
    if rev is None:
        boxes = _BoxWeights(*_expand(ws, shape), wd, _TINY) if ws is not None else None
        _map_slices(vol, m, wd, single_dist_map, overlaps, boxes, nt, power=True)
        return res
    _map_slices(vol, m, wd, single_dist_map, overlaps, None, nt, power=True)
    res[...] = _mpermute(m, rev)
    if ws is not None:
        boxes = _BoxWeights(*_expand(ws, shape), wd, _TINY)
        _parallel(nt, lambda z: np.multiply(res[z], boxes.of(z), out=res[z]), range(shape[0]))
    return res


def _map_slices(vol, out, wd, single, overlaps, boxes, nt, *, power: bool) -> None:
    """the map of a (turned) tile into out: distances, z window, zeros, box weights"""
    nz = vol.shape[0]
    three_d = nz > 1
    pw = wd if power else None
    win = z_window(nz, wd if power else None)
    if power and wd == 0:
        # fastPower(x, 0) is the double scalar 1: the map is 1 on the tile's data
        def zero_power(z: int) -> None:
            keep = np.not_equal(vol[z], 0)
            _finish_slice(out[z], None, 1.0, keep, boxes.of(z) if boxes is not None else None)

        _parallel(nt, zero_power, range(nz))
        return
    if single:
        c = int(_mround((nz + 1) / 2)) - 1
        e_c = _empty(vol, c, three_d)
        d_c = _distance(e_c, pw)

        def one(z: int) -> None:
            keep = np.not_equal(vol[z], 0)
            e = _empty(vol, z, three_d, keep)
            src = d_c if np.array_equal(e, e_c) else _distance(e, pw)
            _finish_slice(out[z], src, win[z], keep, boxes.of(z) if boxes is not None else None)

        _parallel(nt, one, range(nz))
        return
    out.fill(1)
    for box in overlaps:
        _overlap_distance(vol, out, box, pw, three_d, nt)

    def rest(z: int) -> None:
        keep = np.not_equal(vol[z], 0)
        _finish_slice(out[z], out[z], win[z], keep, boxes.of(z) if boxes is not None else None)

    _parallel(nt, rest, range(nz))


def _overlap_distance(vol, out, box, pw, three_d, nt) -> None:
    """the distance in one overlap region (padded by a voxel; MATLAB's crop3d_mex) into out"""
    b = [int(v) for v in box]
    shape = vol.shape
    if any(b[k + 3] <= b[k] for k in range(3)):
        return
    lo = [max(b[k] - 1, 0) for k in range(3)]
    hi = [min(b[k + 3] + 1, shape[k]) for k in range(3)]
    ys, xs = slice(lo[1], hi[1]), slice(lo[2], hi[2])

    def crop(z: int) -> np.ndarray:
        return _empty(vol, z, three_d, ys=ys, xs=xs)

    # skip it when the padded crop holds no empty voxel
    if not any(crop(z).any() for z in range(lo[0], hi[0])):
        return
    if any(b[k] < 0 or b[k + 3] > shape[k] for k in range(3)):
        # (MATLAB: the crop of im_dist_ij at [s, t] indexes beyond it)
        raise DistanceError("Index in position 1 exceeds array bounds.", "MATLAB:badsubscript")
    iy = slice(b[1] - lo[1], b[4] - lo[1])
    ix = slice(b[2] - lo[2], b[5] - lo[2])

    def one(z: int) -> None:
        d = _distance(crop(z), pw)
        out[z, b[1] : b[4], b[2] : b[5]] = 1.0 if d is None else d[iy, ix]

    _parallel(nt, one, range(b[0], b[3]))


def tile_distance_map_mip_slab(
    slab: np.ndarray,
    weight_degree: float = 10,
    *,
    pool_size: Sequence[int],
    dist_bbox: Sequence[int | None] | None = None,
    data_order_axes: Sequence[int] = (2, 1, 0),
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """The weights of a tile's z MIP slab (``(z, y, x)``) as ``compute_tile_bwdist_mip_slabs``
    computes them (module docstring): no powers, the box (in the tile's full-size voxels) scaled
    by ``pool_size`` (``(z, y, x)``, 3 values or 6: along z, then across). A tile whose z is
    not its last axis is turned to put z's dimension last and back."""
    wd = float(weight_degree)
    dom = _dom(data_order_axes)
    vol = _as_volume(slab)
    # (the dimension of z)
    zdo = dom[2]
    if dom[2] != 3:
        vol = _mpermute(vol, [d for d in (1, 2, 3) if d != zdo] + [zdo])
    nt = _threads(num_threads)
    shape = vol.shape
    back = None
    if dom[2] != 3:
        inds = [1, 2, 3]
        first = [i for i in inds if i < zdo]
        back = first + [3] + [i for i in inds if i not in first + [3]]
        q = _mpermute_axes(back)
        shape = tuple(vol.shape[q[k]] for k in range(3))
    # (MATLAB's sizes, of the map turned back)
    sz_m = (shape[1], shape[2], shape[0])
    res = _check_out(out, shape)
    m = res if back is None else np.empty(vol.shape, dtype=np.float32)
    _map_slices(vol, m, None, True, (), None, nt, power=False)
    if back is not None:
        res[...] = _mpermute(m, back)
    if dist_bbox is not None and len(dist_bbox):
        from petakit5d.core.bbox import box_to_matlab
        from petakit5d.matlabcompat import elementary as el

        p = _pool_matlab(pool_size)
        div = [p[3], p[4], p[2], p[3], p[4], p[2]]
        box_m = [max(1.0, _mround(float(v) / d)) for v, d in zip(box_to_matlab(dist_bbox), div,
                                                                 strict=True)]  # fmt: skip
        bufs = [max(1.0, _mround(DIST_BUFFER_SIZE / d)) for d in (p[3], p[4], p[2])]
        dfs = [float(el.power(DIST_DFACTOR, float(d))) for d in (p[3], p[4], p[2])]
        ws = _box_axis_weights(sz_m, box_m, bufs, dfs)
        _check_weights(ws, shape)
        with np.errstate(divide="ignore"):
            floor = float(el.power(1e-50, 1.0 / wd)) if wd != 0 else 0.0
        boxes = _BoxWeights(*_expand(ws, shape), None, np.float32(floor))
        _parallel(nt, lambda z: np.multiply(res[z], boxes.of(z), out=res[z]), range(shape[0]))
    return res


def _pool_matlab(pool_size: Sequence[int]) -> list[float]:
    """MATLAB's poolSize ([y x z], [y x z y x z]) of a (z, y, x) pool size, padded with ones to
    six values as compute_tile_bwdist_mip_slabs pads it"""
    from petakit5d.core.axes import to_yxz

    ps = [float(v) for v in np.asarray(pool_size, dtype=np.float64).ravel()]
    if len(ps) not in (3, 6, 9):
        raise ValueError(f"pool_size has 3, 6 or 9 values, got {len(ps)}")
    m = [float(v) for g in range(0, len(ps), 3) for v in to_yxz(ps[g : g + 3])]
    return m + [1.0, 1.0, 1.0] if len(m) == 3 else m


def mip_slab_name(pool_size: Sequence[int]) -> str:
    """``MIP_Slabs_pooling_<p1>_<p2>_<p3>[_<p4>_<p5>_<p6>]`` (MATLAB's ``[y x z]`` pool sizes,
    ``%d``) of a ``(z, y, x)`` pool size"""
    from petakit5d.core.text import sprintf

    p = _pool_matlab(pool_size)
    n = 3 if len(np.asarray(pool_size).ravel()) == 3 else 6
    return "MIP_Slabs_pooling_" + "_".join(sprintf("%d", v) for v in p[:n])


# ---------------------------------------------------------------------------------------------
# the tasks: compute_tile_bwdist, compute_tile_bwdist_mip_slabs
# ---------------------------------------------------------------------------------------------

# the block information last read (MATLAB's persistent variables): path, mtime, size, values
_BLOCK_INFO: list[Any] = [None, None, None, None]


def _cells(v: Any) -> list[Any]:
    """the elements of a MATLAB cell read by read_mat, in MATLAB's linear order"""
    if v is None:
        return []
    if isinstance(v, np.ndarray) and v.dtype == object:
        return list(v.ravel(order="F"))
    if isinstance(v, (list, tuple)):
        return list(v)
    return [v]


def _read_block_info(path: str) -> dict[str, Any]:
    """tileFns, overlap_matrix, ol_region_cell and overlap_map_mat of a block information file,
    read again only when its path or time changes"""
    from petakit5d.io.matfile import read_mat

    st = os.stat(path)
    key = (st.st_mtime_ns, st.st_size)
    cached = _BLOCK_INFO
    if cached[0] is not None and cached[0].lower() == path.lower() and cached[1] == key:
        return cached[3]
    a = read_mat(path, ["tileFns", "overlap_matrix", "ol_region_cell", "overlap_map_mat"])
    tiles = [str(t) for t in _cells(a["tileFns"])]
    om = np.asarray(a["overlap_matrix"])
    regions = [[np.asarray(r, dtype=np.float64) for r in _cells(c)]
               for c in _cells(a["ol_region_cell"])]  # fmt: skip
    omap = np.asarray(a["overlap_map_mat"], dtype=np.float64).reshape(-1, 2, order="F") \
        if np.size(a["overlap_map_mat"]) else np.zeros((0, 2))  # fmt: skip
    info = {"tileFns": tiles, "overlap_matrix": om, "ol_region_cell": regions,
            "overlap_map_mat": omap}  # fmt: skip
    _BLOCK_INFO[:] = [path, key, None, info]
    return info


def _region_box(m, ol_order: Sequence[int]) -> tuple[int, ...]:
    """MATLAB's region ``[s t]`` (rows x, y, z) as a box of the (turned) array whose MATLAB
    dimensions hold the rows ``ol_order`` (MATLAB's ``[midx(olOrder, 1)', midx(olOrder, 2)']``)"""
    r = np.asarray(m, dtype=np.float64).reshape(3, 2)
    o = [int(v) - 1 for v in ol_order]
    return (
        int(r[o[2], 0]) - 1,
        int(r[o[0], 0]) - 1,
        int(r[o[1], 0]) - 1,
        int(r[o[2], 1]),
        int(r[o[0], 1]),
        int(r[o[1], 1]),
    )


def _tile_overlaps(
    info: dict[str, Any], i: int, ol_order: Sequence[int] = (2, 1, 3)
) -> list[tuple[int, ...]]:
    """the overlap regions of tile i (0-based) in the order of the other tiles: each pair's
    first region for the lower tile, its second for the higher, as boxes of the (turned) array
    whose MATLAB dimensions hold the regions' rows ``ol_order`` (x, y, z: 1, 2, 3)"""
    n = len(info["tileFns"])
    om = np.asarray(info["overlap_matrix"]).reshape(n, n, order="F")
    omap = info["overlap_map_mat"]
    out = []
    for j in range(n):
        if not om[i, j] and not om[j, i]:
            continue
        ind = (j * n + i + 1) if i < j else (i * n + j + 1)
        rows = np.nonzero(omap[:, 1] == ind)[0]
        if rows.size != 1:
            # (MATLAB: ol_region_cell{...}{1} of no or several cells)
            raise DistanceError(
                f"Intermediate brace '{{}}' indexing produced a comma-separated list with "
                f"{rows.size} values, but it must produce a single value when followed by "
                "subsequent indexing operations.",
                "MATLAB:index:expected_one_output_from_intermediate_indexing",
            )
        pair = info["ol_region_cell"][int(rows[0])]
        out.append(_region_box(pair[0 if i < j else 1], ol_order))
    return out


def _dim_separator(n_chunks: float) -> str:
    return "/" if n_chunks > 10000 else "."


def _publish(tmp: str, final: str, overwrite: bool) -> None:
    """MATLAB's end of the tasks: an output that appeared meanwhile is kept (the new one removed)
    unless overwrite; else the temporary zarr is renamed"""
    if os.path.isdir(final):
        if not overwrite:
            print(f"Distance Zarr result {final} already exists, skip it!", flush=True)
            shutil.rmtree(tmp)
            return
        shutil.rmtree(final)
    os.rename(tmp, final)


def _write_map(
    m: np.ndarray,
    final: str,
    block: Sequence[int],
    sep: str,
    nt: int,
    *,
    shard_size=None,
    compressor: str = "zstd",
    overwrite: bool = False,
) -> None:
    """the map written as MATLAB writes it: createzarr of its size under a temporary name,
    writezarr, then the rename"""
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import create_zarr, write_zarr

    tmp = f"{final}_{get_uuid()}"
    create_zarr(tmp, data_size=tuple(m.shape), block_size=block, shard_size=shard_size,
                dtype="single", compressor=compressor, dim_separator=sep)  # fmt: skip
    write_zarr(m, tmp, num_threads=nt, bbox=(0, 0, 0, *m.shape))
    _publish(tmp, final, overwrite)


_BLOCK_INFO_PARAM = Param("blockInfoFullname", "@ischar", "Block information file (stitching).")
_TILE_PARAM = Param(
    "tileInd", "@isscalar", "The tile (0-based) in the block information.", convert="index1"
)
_DATA_ORDER_PARAM = Param(
    "dataOrderMat",
    "@(x) isvector(x)",
    "Array axes of x, y and z of the tile (stitch.geometry.parse_data_order).",
    py="data_order_axes",
    convert=(data_order_axes_from_matlab, data_order_axes_to_matlab),
)
_THREADS_PARAM = Param(
    "numThreads",
    None,
    "Threads of the map and the zarr reads and writes (None: the CPUs available).",
    py="num_threads",
    python_only=True,
)


@register(
    "compute_tile_bwdist",
    _BLOCK_INFO_PARAM,
    _TILE_PARAM,
    Param("bwdistFullpath", "@ischar", "Distance map to write (zarr)."),
    Param("weightDegree", "@isscalar", "Power of the weights (an integer)."),
    Param("singleDistMap", "@islogical", "The map from the tile alone, not its overlap regions."),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of the map.", convert="size_yxz"),
    Param(
        "shardSize",
        "@(x) isempty(x) || isvector(x)",
        "Not used (as in MATLAB).",
        convert="size_yxz",
    ),  # fmt: skip
    Param("compressor", "@ischar", "Not used (as in MATLAB: zstd)."),
    Param(
        "distBbox",
        "@(x) isempty(x) || isvector(x)",
        "Box (z0, y0, x0, z1, y1, x1) of full weight (None: none).",
        convert="bbox",
    ),
    _DATA_ORDER_PARAM,
    Param("Overwrite", "@islogical", "Replace a map that appears while this one is computed."),
    _THREADS_PARAM,
    task="S06",
)
def compute_tile_bwdist(
    block_info_fullname: str,
    tile_ind: int,
    bwdist_fullpath: str,
    weight_degree: float,
    single_dist_map: bool,
    block_size: Sequence[int],
    shard_size: Sequence[int] | None,
    compressor: str,
    dist_bbox: Sequence[int | None] | None,
    data_order_axes: Sequence[int],
    overwrite: bool = False,
    *,
    num_threads: int | None = None,
) -> None:
    """``compute_tile_bwdist``: the distance map of tile ``tile_ind`` of the block information
    (``tile_distance_map`` with its overlap regions), written as a native single zarr (chunks
    ``block_size``, zstd; the separator ``/`` above 10,000 chunks) under a temporary name and
    renamed to ``bwdist_fullpath``. Nothing is done if that exists."""
    del shard_size, compressor  # (MATLAB's createzarr call takes neither)
    from petakit5d.io import read_zarr

    final = os.fspath(bwdist_fullpath)
    if os.path.isdir(final):
        return
    info = _read_block_info(os.fspath(block_info_fullname))
    i = int(tile_ind)
    nt = _threads(num_threads)
    tile = read_zarr(info["tileFns"][i], num_threads=nt)
    # the rows (x, y, z) of the overlap regions in the order of the turned array's axes: x, y, z
    # once turned (to put z last), else the tile's own order
    dom = _dom(data_order_axes)
    ol_order = (1, 2, 3) if dom[2] != 3 else tuple(dom)
    overlaps = () if single_dist_map else _tile_overlaps(info, i, ol_order)
    m = tile_distance_map(
        tile, weight_degree, single_dist_map=bool(single_dist_map), overlaps=overlaps,
        dist_bbox=dist_bbox, data_order_axes=data_order_axes, num_threads=nt,
    )  # fmt: skip
    del tile
    block = [int(v) for v in block_size]
    n_chunks = float(np.prod([math.ceil(s / b) for s, b in zip(m.shape, block, strict=True)]))
    _write_map(m, final, block, _dim_separator(n_chunks), nt, overwrite=bool(overwrite))


def _mrdivide_rows(a: Sequence[float], b: Sequence[float]) -> float:
    """MATLAB's ``a / b`` of two rows: the least-squares x of ``x * b = a``"""
    a64 = np.asarray(a, dtype=np.float64).ravel()
    b64 = np.asarray(b, dtype=np.float64).ravel()
    return float(np.dot(b64, a64) / np.dot(b64, b64))


@register(
    "compute_tile_bwdist_mip_slabs",
    _BLOCK_INFO_PARAM,
    _TILE_PARAM,
    Param("bwdistFullpath", "@ischar", "Distance map to write (zarr)."),
    Param("weightDegree", "@isscalar", "Power of the weights, applied in the blending."),
    Param("singleDistMap", "@islogical", "Not used (always the tile's own map)."),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of the map.", convert="size_yxz"),
    Param(
        "shardSize",
        "@(x) isempty(x) || isvector(x)",
        "Shard size (z, y, x) (None: none).",
        convert="size_yxz",
    ),  # fmt: skip
    Param("compressor", "@ischar", "Compressor of the map."),
    Param("poolSize", "@isvector", "Pool sizes (z, y, x) of the MIP slabs.", convert="size_yxz"),
    Param(
        "distBbox",
        "@isnumeric",
        "Box (z0, y0, x0, z1, y1, x1) of full weight in the tile.",
        convert="bbox",
    ),  # fmt: skip
    _DATA_ORDER_PARAM,
    Param("Overwrite", "@islogical", "Replace a map that appears while this one is computed."),
    _THREADS_PARAM,
    task="S06",
)
def compute_tile_bwdist_mip_slabs(
    block_info_fullname: str,
    tile_ind: int,
    bwdist_fullpath: str,
    weight_degree: float,
    single_dist_map: bool,
    block_size: Sequence[int],
    shard_size: Sequence[int] | None,
    compressor: str,
    pool_size: Sequence[int],
    dist_bbox: Sequence[int | None] | None,
    data_order_axes: Sequence[int],
    overwrite: bool = False,
    *,
    num_threads: int | None = None,
) -> None:
    """``compute_tile_bwdist_mip_slabs``: the map of tile ``tile_ind``'s z MIP slab
    (``<tile folder>/MIP_Slabs_pooling_<pool>/<name>_z.zarr``; ``tile_distance_map_mip_slab``),
    written as a native single zarr with ``shard_size`` and ``compressor`` under a temporary
    name and renamed to ``bwdist_fullpath``. Nothing is done if that exists."""
    del single_dist_map
    from petakit5d.core.paths import fileparts
    from petakit5d.io import read_zarr

    final = os.fspath(bwdist_fullpath)
    if os.path.isdir(final):
        return
    info = _read_block_info(os.fspath(block_info_fullname))
    data_path, fsn, _ = fileparts(info["tileFns"][int(tile_ind)])
    mip = f"{data_path}/{mip_slab_name(pool_size)}/{fsn}_z.zarr"
    nt = _threads(num_threads)
    slab = read_zarr(mip, num_threads=nt)
    m = tile_distance_map_mip_slab(
        slab, weight_degree, pool_size=pool_size, dist_bbox=dist_bbox,
        data_order_axes=data_order_axes, num_threads=nt,
    )  # fmt: skip
    del slab
    # the separator from ceil(sz / blockSize), a matrix division
    block = [int(v) for v in block_size]
    sz_m = [m.shape[1], m.shape[2], m.shape[0]]
    from petakit5d.core.axes import to_yxz

    q = _mrdivide_rows(sz_m, to_yxz(block))
    sep = _dim_separator(float(math.ceil(q)))
    shard = None if shard_size is None or len(shard_size) == 0 else [int(v) for v in shard_size]
    _write_map(m, final, block, sep, nt, shard_size=shard, compressor=compressor,
               overwrite=bool(overwrite))  # fmt: skip


# ---------------------------------------------------------------------------------------------
# compute_tile_distance_transform: the wrapper
# ---------------------------------------------------------------------------------------------


def _is_pc() -> bool:
    """MATLAB's ``ispc`` (asked at each call, so that tests can take the Windows branch)"""
    return sys.platform.startswith("win")


def _unique_first(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """MATLAB's ``[c, ia] = unique(v)``: sorted values and the first index of each"""
    c, ia = np.unique(v, return_index=True)
    return c, ia


def _same_setting(path: str, om, regions_raw, single: bool) -> bool:
    from petakit5d.io.matfile import read_mat

    if single:
        return True
    b = read_mat(path, ["overlap_matrix", "ol_region_cell", "overlap_map_mat"])
    om_b = np.asarray(b["overlap_matrix"])
    om_a = np.asarray(om)
    if om_a.shape != om_b.shape or not np.all(om_a == om_b):
        return False
    ca, cb = _cells(regions_raw), _cells(b["ol_region_cell"])
    if len(ca) != len(cb):
        # (MATLAB: cellfun of cells of other sizes)
        raise DistanceError("All of the input arguments must be of the same size and shape.",
                            "MATLAB:cellfun:AllInputsSameSize")  # fmt: skip
    return all(_isequal(x, y) for x, y in zip(ca, cb, strict=True))


def _isequal(a: Any, b: Any) -> bool:
    """MATLAB's isequal of two values read by read_mat (cells as lists or object arrays)"""
    ca = isinstance(a, (list, tuple)) or (isinstance(a, np.ndarray) and a.dtype == object)
    cb = isinstance(b, (list, tuple)) or (isinstance(b, np.ndarray) and b.dtype == object)
    if ca or cb:
        if not (ca and cb):
            return False
        xa, xb = _cells(a), _cells(b)
        return len(xa) == len(xb) and all(_isequal(x, y) for x, y in zip(xa, xb, strict=True))
    x, y = np.asarray(a), np.asarray(b)
    return x.shape == y.shape and bool(np.all(x == y))


def compute_tile_distance_transform(
    block_info_fullname: str,
    stitch_path: str,
    tile_fullpaths: Sequence[str],
    *,
    overwrite: bool = False,
    data_order_axes: Sequence[int] = (2, 1, 0),
    blend_weight_degree: float = 10,
    single_dist_map: bool = True,
    loc_ids: Sequence[float] | None = None,
    dist_bboxes: Sequence[Sequence[int | None]] | None = None,
    block_size: Sequence[int] = (500, 500, 500),
    shard_size: Sequence[int] | None = None,
    compressor: str = "lz4",
    large_file: bool = False,
    pool_size: Sequence[int] | None = None,
    parse_cluster: bool = True,
    mcc_mode: bool = False,
    config_file: str = "",
    uuid: str = "",
) -> tuple[list[str], list[int]]:
    """``[outputFullpaths, imdistFileIdx] = compute_tile_distance_transform(blockInfoFullname,
    stitchPath, tileFullpaths, ...)``: the distance maps of the tiles, one ``compute_tile_bwdist``
    task each (``compute_tile_bwdist_mip_slabs`` on MIP slabs for ``large_file``, made first),
    through the job framework, in ``<stitch_path>/imdist/<first tile's name>/`` (``imdist/`` on
    Windows) with the call's ``parameters.mat``. With ``single_dist_map`` there is one map per
    location (``loc_ids``), and locations whose first tiles have the same size share one; else
    one per tile. Returns the maps (MATLAB's paths) and the map of each tile (0-based; one value
    when all tiles share a map). An existing folder made with other overlap regions (without
    ``single_dist_map``) is removed first."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive
    from petakit5d.core.text import sprintf, through_text
    from petakit5d.io import get_image_size, read_mat, write_mat
    from petakit5d.jobs.framework import run_tasks
    from petakit5d.jobs.tasks import Task

    if not uuid:
        uuid = get_uuid()
    paths = [os.fspath(p) for p in tile_fullpaths]
    n_f = len(paths)
    _, fsn, _ = fileparts(paths[0])
    if _is_pc():
        dist_path = sprintf("%s/imdist/", stitch_path)
    else:
        dist_path = sprintf("%s/imdist/%s/", stitch_path, fsn)
    if not os.path.isdir(dist_path):
        mkdir_recursive(dist_path)
    a = read_mat(block_info_fullname, ["overlap_matrix", "ol_region_cell", "overlap_map_mat"])
    om, regions, omap = a["overlap_matrix"], a["ol_region_cell"], a["overlap_map_mat"]
    param_fullname = dist_path + "parameters.mat"
    if os.path.isfile(param_fullname):
        if not _same_setting(param_fullname, om, regions, bool(single_dist_map)):
            print(
                "The existing distance transforms are not in the same setting as current "
                "running, delete them!",
                flush=True,
            )
            shutil.rmtree(dist_path)
    os.makedirs(dist_path, exist_ok=True)
    pr = _results(block_info_fullname, stitch_path, paths, overwrite, data_order_axes,
                  blend_weight_degree, single_dist_map, loc_ids, dist_bboxes, block_size,
                  shard_size, compressor, large_file, pool_size, parse_cluster, mcc_mode,
                  config_file, uuid)  # fmt: skip
    tmp = sprintf("%s/parameters_%s.mat", dist_path, uuid)
    # (MATLAB saves the inputParser itself, which the port can't write: its Results)
    write_mat(tmp, {"ip": {"Results": pr}, "overlap_matrix": om, "ol_region_cell": regions,
                    "overlap_map_mat": omap}, version="7.3")  # fmt: skip
    os.replace(tmp, param_fullname)

    file_idx = list(range(1, n_f + 1))
    dist_tile_idx = list(range(1, n_f + 1))
    boxes_m = None
    if dist_bboxes is not None and len(dist_bboxes):
        from petakit5d.core.bbox import box_to_matlab

        boxes_m = np.array([box_to_matlab(b) for b in dist_bboxes], dtype=np.float64)
    if single_dist_map:
        lid = np.asarray([] if loc_ids is None else loc_ids, dtype=np.float64).ravel()
        # (without location ids, the tiles are of one location)
        if lid.size == 0:
            lid = np.ones(n_f)
        uniq, first = _unique_first(lid)
        paths = [paths[k] for k in first]
        file_idx = [0] * n_f
        sizes = [get_image_size(paths[k]) for k in range(uniq.size)]
        for k in range(uniq.size):
            for t in np.nonzero(lid == uniq[k])[0]:
                file_idx[int(t)] = k + 1
        dist_tile_idx = [int(k) + 1 for k in first]
        if uniq.size > 1:
            mapped = list(range(1, uniq.size + 1))
            for k in range(1, uniq.size):
                same = [j for j in range(k) if tuple(sizes[j]) == tuple(sizes[k])]
                if same:
                    mapped[k] = same[0] + 1
            um = sorted(set(mapped))
            paths = [paths[k - 1] for k in um]
            for k in range(uniq.size):
                for t in np.nonzero(lid == uniq[k])[0]:
                    file_idx[int(t)] = um.index(mapped[k]) + 1
            dist_tile_idx = [dist_tile_idx[k - 1] for k in um]
        if boxes_m is not None:
            if boxes_m.shape[0] == 1:
                if len(dist_tile_idx) > 1:
                    boxes_m = np.repeat(boxes_m, len(dist_tile_idx), axis=0)
            elif boxes_m.shape[0] != len(dist_tile_idx):
                raise DistanceError(
                    "Number of distance bounding box does not match the number of unique tile "
                    "image sizes!"
                )
    elif boxes_m is not None and boxes_m.shape[0] == 1 and n_f > 1:
        # (one box for the maps of all tiles)
        boxes_m = np.repeat(boxes_m, n_f, axis=0)
    if file_idx and all(v == file_idx[0] for v in file_idx):
        file_idx = file_idx[:1]
    n_f = len(paths)
    if n_f == 0:
        # (MATLAB: tileFullpaths{1} of an empty cell)
        raise DistanceError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    fsnames = [fileparts(p)[1] for p in paths]
    box_rows: list[Any] = [None] * n_f
    if boxes_m is not None:
        if boxes_m.shape[0] < n_f:
            raise DistanceError(
                f"Index in position 1 exceeds array bounds. Index must not exceed "
                f"{boxes_m.shape[0]}.",
                "MATLAB:badsubscript",
            )
        from petakit5d.core.bbox import box_from_matlab
        from petakit5d.core.text import through_mat2str

        box_rows = [list(box_from_matlab(through_mat2str(boxes_m[k]))) for k in range(n_f)]
    wd = float(through_text(blend_weight_degree, "%d"))
    wd_out = int(wd) if wd == math.floor(wd) else wd
    block = [int(v) for v in block_size]
    shard = None if shard_size is None or len(shard_size) == 0 else [int(v) for v in shard_size]
    dao = [int(v) for v in data_order_axes]
    if large_file:
        from petakit5d.stitch.assignment import _tile_mip_slabs

        pool = [int(v) for v in np.asarray(pool_size).ravel()]
        _tile_mip_slabs(paths, pool, parse_cluster, mcc_mode, config_file)
        pool_str = mip_slab_name(pool).removeprefix("MIP_Slabs_")
        outputs = [sprintf("%s/%s_z_%s_wr_1_for_%d.zarr", dist_path, x, pool_str,
                           blend_weight_degree) for x in fsnames]  # fmt: skip
        tasks = [
            Task.of(
                compute_tile_bwdist_mip_slabs,
                block_info_fullname,
                dist_tile_idx[k] - 1,
                outputs[k],
                wd_out,
                bool(single_dist_map),
                block,
                shard,
                compressor,
                pool,
                box_rows[k],
                dao,
                bool(overwrite),
            )  # fmt: skip
            for k in range(n_f)
        ]
    else:
        outputs = [sprintf("%s/%s_wr_%d.zarr", dist_path, x, blend_weight_degree)
                   for x in fsnames]  # fmt: skip
        tasks = [
            Task.of(
                compute_tile_bwdist,
                block_info_fullname,
                dist_tile_idx[k] - 1,
                outputs[k],
                wd_out,
                bool(single_dist_map),
                block,
                shard,
                compressor,
                box_rows[k],
                dao,
                bool(overwrite),
            )  # fmt: skip
            for k in range(n_f)
        ]
    im_size = get_image_size(paths[0])
    mem = float(np.prod([float(v) for v in im_size])) * 4 / 1024**3 * TASK_MEMORY_FACTOR
    min_task_job_num = len(tasks) if mem > 10 else 1
    run_tasks(paths, outputs, tasks, mem_allocate=mem, min_task_job_num=min_task_job_num,
              parse_cluster=parse_cluster, mcc_mode=mcc_mode, config_file=config_file)  # fmt: skip
    return outputs, [v - 1 for v in file_idx]


def _results(
    block_info_fullname,
    stitch_path,
    paths,
    overwrite,
    data_order_axes,
    blend_weight_degree,
    single_dist_map,
    loc_ids,
    dist_bboxes,
    block_size,
    shard_size,
    compressor,
    large_file,
    pool_size,
    parse_cluster,
    mcc_mode,
    config_file,
    uuid,
) -> dict[str, Any]:
    """the call's ip.Results with MATLAB's values, fields sorted as inputParser sorts them"""
    from petakit5d.core.axes import to_yxz
    from petakit5d.core.bbox import box_to_matlab

    def num(v):
        return np.zeros((0, 0)) if v is None or len(v) == 0 else np.asarray(v, dtype=np.float64)

    pool = None
    if pool_size is not None and len(pool_size):
        ps = [float(v) for v in np.asarray(pool_size).ravel()]
        pool = [v for g in range(0, len(ps), 3) for v in to_yxz(ps[g : g + 3])]
    boxes = None
    if dist_bboxes is not None and len(dist_bboxes):
        boxes = np.array([box_to_matlab(b) for b in dist_bboxes], dtype=np.float64)
    pr = {
        "Overwrite": bool(overwrite),
        "blendWeightDegree": float(blend_weight_degree),
        "blockInfoFullname": os.fspath(block_info_fullname),
        "blockSize": num(to_yxz(block_size)),
        "compressor": compressor,
        "configFile": config_file,
        "dataOrderMat": num(data_order_axes_to_matlab(data_order_axes)),
        "distBboxes": num(boxes),
        "largeFile": bool(large_file),
        "locIds": num(None if loc_ids is None else list(loc_ids)),
        "mccMode": bool(mcc_mode),
        "parseCluster": bool(parse_cluster),
        "poolSize": num(pool),
        "shardSize": num(None if shard_size is None or not len(shard_size) else to_yxz(shard_size)),
        "singleDistMap": bool(single_dist_map),
        "stitchPath": os.fspath(stitch_path),
        "tileFullpaths": np.array(list(paths), dtype=object).reshape(-1, 1),
        "uuid": uuid,
    }
    return dict(sorted(pr.items(), key=lambda kv: kv[0].lower()))
