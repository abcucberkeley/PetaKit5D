"""The stitching kernels: feather and tile blending, distance map resize, integral image, the
major-tile check and NaN replacement (``petakit5d._ext.stitch``, moved from the mex files of
``microscopeDataProcessing/stitch/mex`` and ``pointDetection/mex_functions``).

Arrays are the port's (``docs/conventions.md``): volumes ``(z, y, x)``, images ``(y, x)``, a
stack of tiles ``(tile, z, y, x)``, all C-contiguous (never copied; another layout raises).
Boxes are ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive; per-axis values are ``(z, y,
x)``. Every function takes ``num_threads`` (None: the CPUs this process may use); results do not
depend on it.

Two levels:

- the kernels (``feather_blend``, ``feather_blend_into``, ``resize_distance_map``,
  ``tile_blend``, ``integral_image``, ``major_tiles_valid``, ``replace_nan_with_value``): the
  mex's checks with its error identifiers (``StitchKernelError.identifier``) and its results
  bit for bit: each voxel gets the mex's operations in the mex's order (tiles in order, the
  integral image's passes along y, then x, then z, the interpolation along y, then x, then z),
  in single precision and without fused multiply-adds;
- the MATLAB wrappers (``feather_blending_3d``, ``feather_distance_map_resize_3d``,
  ``tile_blending_3d``, ``integral_image_3d``, ``check_major_tile_valid``): the kernel for the
  classes the mex takes, MATLAB's own code for the others (what the wrapper runs when the mex
  rejects a class). Where the MATLAB wrapper falls back to its code after a mex error that is
  not about a class (a box outside the output, sizes that do not match), these raise the mex's
  error instead.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np

__all__ = [
    "BLEND_METHODS",
    "StitchKernelError",
    "StitchKernelTypeError",
    "check_major_tile_valid",
    "feather_blend",
    "feather_blend_into",
    "feather_blending_3d",
    "feather_distance_map_resize_3d",
    "integral_image",
    "integral_image_3d",
    "major_tiles_valid",
    "replace_nan_with_value",
    "resize_distance_map",
    "tile_blend",
    "tile_blending_3d",
    "tile_boxes_from_matlab",
    "tile_boxes_to_matlab",
]

BLEND_METHODS = ("none", "mean", "max", "median")


class StitchKernelError(ValueError):
    """An error the mex raises; ``identifier`` is its MATLAB error identifier."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


class StitchKernelTypeError(TypeError):
    """A class the mex rejects; ``identifier`` is its MATLAB error identifier."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _ext():
    from petakit5d._ext import load

    return load("stitch")


def _threads(n: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(n)


def _array(a, name: str) -> np.ndarray:
    if not isinstance(a, np.ndarray):
        raise TypeError(f"{name} must be a NumPy array, got {type(a).__name__}")
    return a


def _mndims(a: np.ndarray) -> int:
    """MATLAB's ndims of the array of the same data"""
    from petakit5d.core.arrays import size_to_matlab

    return len(size_to_matlab(a.shape))


def _pad(shape: Sequence[int], n: int) -> list[int]:
    """a shape with leading singletons to n dimensions"""
    s = [int(v) for v in shape]
    return [1] * (n - len(s)) + s


def _flat(a: np.ndarray, name: str) -> np.ndarray:
    """a 1-D view of a C-contiguous array (never a copy)"""
    if not a.flags.c_contiguous:
        hint = " (it is Fortran-ordered: MATLAB's layout? use core.arrays.from_matlab)"
        raise ValueError(f"{name} must be C-contiguous{hint if a.flags.f_contiguous else ''}")
    return a.reshape(-1)


def _floor(v: float) -> float:
    """C's floor (NaN and Inf pass through)"""
    return float(np.floor(v))


def _trunc_u64(v: float) -> int:
    """(uint64_t) of a double as GCC converts it on x86-64: truncation, negative values wrap;
    NaN, +Inf and values from 2^64 give 0, values to -2^63 give 2^63 (every caller checks the
    result, so these only decide which error comes)"""
    if v != v or v >= 2.0**64:
        return 0
    if v <= -(2.0**63):
        return 1 << 63
    return math.trunc(v) % (1 << 64)


def _numbers(v, name: str) -> np.ndarray:
    """real numbers of an argument as float64"""
    a = np.asarray(v)
    if a.dtype.kind not in "biuf":
        raise TypeError(f"{name} must be real numbers, got {a.dtype}")
    return a.astype(np.float64).reshape(-1)


def _ints(v, n: int, name: str) -> list[int]:
    """n integers of a box or a size (Python ints, NumPy integers, integral floats)"""
    a = _numbers(v, name)
    if a.size != n:
        raise ValueError(f"{name} must have {n} values, got {a.size}")
    if not np.all(np.isfinite(a)) or np.any(a != np.floor(a)):
        raise ValueError(f"{name} must be integers, got {list(a)}")
    return [int(x) for x in a]


# ---------------------------------------------------------------------------------------------
# feather blending
# ---------------------------------------------------------------------------------------------

_FB = "feather_blending:inputError"
_FB_TYPE = "feather_blending:dataTypeError"
_FB_DMAT_SIZE = "The size of distance map does not matach that of the data"


def _feather_checks(f: np.ndarray, d: np.ndarray) -> list[int]:
    """the mex's checks of the tiles and the weights: their (t, z, y, x) sizes"""
    if _mndims(f) > 4:
        raise StitchKernelError(_FB, "The data must have at most 4 dimensions")
    if _mndims(d) != _mndims(f):
        raise StitchKernelError(_FB, _FB_DMAT_SIZE)
    return _pad(f.shape, 4)


def feather_blend(
    tiles: np.ndarray, weights: np.ndarray, *, num_threads: int | None = None
) -> np.ndarray:
    """``feather_blending_3d_mex``: the feather blend of the tiles (uint16 or single,
    ``(t, z, y, x)``; a volume is one tile) with the single weights of the same size, per voxel
    ``sum(f * w) / sum(w)`` over the tiles with ``f != 0``, in tile order and in single
    precision; uint16 results are rounded with + 0.5 (in double), voxels without a weight are
    0. Returns ``(z, y, x)`` of the tiles' class."""
    f, d = _array(tiles, "tiles"), _array(weights, "weights")
    dims = _feather_checks(f, d)
    if d.dtype != np.float32:
        raise StitchKernelError(_FB, "The data type of the distance map must be single!")
    if _pad(d.shape, 4) != dims:
        raise StitchKernelError(_FB, _FB_DMAT_SIZE)
    if f.dtype not in (np.uint16, np.float32):
        raise StitchKernelTypeError(_FB_TYPE, "Data type not suppported")
    out = np.zeros(dims[1:], dtype=f.dtype)
    _ext().feather_blending_3d(
        _flat(f, "tiles"), _flat(d, "weights"), _flat(out, "out"), dims, dims[1:], [0, 0, 0],
        _threads(num_threads),
    )  # fmt: skip
    return out


def feather_blend_into(
    tiles: np.ndarray,
    weights: np.ndarray,
    out: np.ndarray,
    box,
    *,
    num_threads: int | None = None,
) -> None:
    """``feather_blending_3d_with_indexing_mex``: the feather blend written into ``out``
    (``(z, y, x)``, the tiles' class) in place, in ``box`` (the size of a tile); voxels without
    a weight keep their value."""
    f, d, o = _array(tiles, "tiles"), _array(weights, "weights"), _array(out, "out")
    dims = _feather_checks(f, d)
    if _mndims(o) > 3:
        raise StitchKernelError(_FB, "The output must have at most 3 dimensions")
    if _pad(d.shape, 4) != dims:
        raise StitchKernelError(_FB, _FB_DMAT_SIZE)
    odims = _pad(o.shape, 3)
    b = _numbers(box, "box")
    if b.size != 6:
        raise StitchKernelError(_FB, "The size of bounding box must be 1x6.")
    b = _ints(b, 6, "box")
    for i in range(3):
        if not (0 <= b[i] and b[i + 3] <= odims[i] and b[i + 3] - b[i] == dims[i + 1]):
            raise StitchKernelError(_FB, "The bounding box is out of bound.")
    if d.dtype != np.float32:
        raise StitchKernelError(_FB, "The data type of the distance map must be single!")
    if f.dtype != o.dtype:
        raise StitchKernelError(
            _FB, "The data type of the overlap region does not match that of the data!"
        )
    if f.dtype not in (np.uint16, np.float32):
        raise StitchKernelTypeError(_FB_TYPE, "Data type not suppported")
    if not o.flags.writeable:
        raise ValueError("out must be writeable")
    _ext().feather_blending_3d(
        _flat(f, "tiles"), _flat(d, "weights"), _flat(o, "out"), dims, odims, b[:3],
        _threads(num_threads),
    )  # fmt: skip


def _sum_tiles(x: np.ndarray) -> np.ndarray:
    """MATLAB's sum(x, 4) of a (y, x, z, t) array: the tiles added in order, in x's class"""
    if x.shape[0] == 0:
        return np.zeros(x.shape[1:], dtype=x.dtype)
    s = x[0].copy()
    for t in range(1, x.shape[0]):
        s = s + x[t]
    return s


def _feather_matlab(f: np.ndarray, d: np.ndarray) -> np.ndarray:
    """the MATLAB code of feather_blending_3d: sum(single(f) .* w, 4) ./ sum(w, 4) with
    w = d .* (f ~= 0) (single; NaN where no tile has a weight)"""
    from petakit5d.core import dtypes

    if d.dtype.kind in "iub":
        raise StitchKernelTypeError(
            "MATLAB:mixedClasses",
            "Integers can only be combined with integers of the same class, or scalar doubles.",
        )
    if f.dtype.kind not in "biuf" or d.dtype.kind != "f":
        raise TypeError(f"feather blending of {f.dtype} tiles with {d.dtype} weights")
    f4 = f.reshape(_pad(f.shape, 4))
    d4 = d.reshape(_pad(d.shape, 4))
    if f4.shape != d4.shape:
        raise StitchKernelError(_FB, _FB_DMAT_SIZE)
    w = d4 * (f4 != 0).astype(d4.dtype)
    fs = dtypes.single(f4)
    p = fs * w if w.dtype == np.float32 else dtypes.times(fs, w)
    num, den = _sum_tiles(p), _sum_tiles(w)
    with np.errstate(divide="ignore", invalid="ignore"):
        if den.dtype == np.float32:
            return np.ascontiguousarray(num / den)
        return np.ascontiguousarray(dtypes.rdivide(num, den))


def feather_blending_3d(
    tim_f_block: np.ndarray,
    tim_d_block: np.ndarray,
    tim_block: np.ndarray | None = None,
    box=None,
    *,
    num_threads: int | None = None,
) -> tuple[np.ndarray, bool]:
    """``[tim_block, mex_compute] = feather_blending_3d(tim_f_block, tim_d_block[, tim_block,
    bbox])``: the feather blend of the tiles (``(t, z, y, x)``), or with ``tim_block`` and
    ``box`` written into ``tim_block`` in place (returned).

    uint16 and single tiles with single weights (and a ``tim_block`` of the tiles' class) run the
    kernel (``mex_compute`` True). Other classes run MATLAB's code: a single result with NaN
    where no tile has a weight; with ``tim_block``, that result cast to the block's class
    (MATLAB's ``cast``: rounded half away from zero, saturated, NaN to 0) and written at ``box``.
    ``indexing4d_mex``'s classes (logical, uint8, uint16, single, double) are written into
    ``tim_block`` in place; for other classes MATLAB writes into a copy, which is returned. A box
    outside the block or of another size than the tiles raises (MATLAB would grow the block or
    fail in its assignment).
    """
    f, d = _array(tim_f_block, "tim_f_block"), _array(tim_d_block, "tim_d_block")
    indexed = tim_block is not None or box is not None
    if indexed and (tim_block is None or box is None):
        raise TypeError("give both tim_block and box, or neither")
    # the classes the mex takes; with others it raises and MATLAB's code runs
    mex_ok = f.dtype in (np.uint16, np.float32) and d.dtype == np.float32
    if indexed:
        mex_ok = mex_ok and _array(tim_block, "tim_block").dtype == f.dtype
    if mex_ok:
        if not indexed:
            return feather_blend(f, d, num_threads=num_threads), True
        feather_blend_into(f, d, tim_block, box, num_threads=num_threads)
        return tim_block, True
    crop = _feather_matlab(f, d)
    if not indexed:
        return crop, False
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import matlab_class

    b = _ints(box, 6, "box")
    o = _array(tim_block, "tim_block")
    odims = _pad(o.shape, 3)
    if _mndims(o) > 3 or not all(0 <= b[i] <= b[i + 3] <= odims[i] for i in range(3)):
        raise StitchKernelError(_FB, "The bounding box is out of bound.")
    if [b[i + 3] - b[i] for i in range(3)] != list(crop.shape):
        raise StitchKernelError(_FB, "The bounding box is out of bound.")
    # tim_block = indexing4d(tim_block, cast(crop, class(tim_block)), [bbox(1:3), 1, bbox(4:6), 1])
    region = dtypes.cast(crop, matlab_class(o))
    if o.dtype not in (np.bool_, np.uint8, np.uint16, np.float32, np.float64):
        o = o.copy()  # indexing4d's MATLAB code: an assignment to its own copy, returned
    o.reshape(odims)[b[0] : b[3], b[1] : b[4], b[2] : b[5]] = region
    return o, False


# ---------------------------------------------------------------------------------------------
# distance map resize
# ---------------------------------------------------------------------------------------------

_DR = "distmap_resize:inputError"
_DR_SMALL = "The resized dimensions must not be smaller than those of the input distance map!"


def resize_distance_map(
    dmat: np.ndarray, size, wd, *, num_threads: int | None = None
) -> np.ndarray:
    """``feather_distance_map_resize_3d_mex(dmat, sz, wd)``: the single distance map
    ``(z, y, x)`` resized to ``size`` ``(nz, ny, nx)`` (not smaller than the map; values are
    truncated) by trilinear interpolation (the corners map onto each other), raised to the
    power ``int(single(wd))``."""
    dm = _array(dmat, "dmat")
    if _mndims(dm) > 3:
        raise StitchKernelError(_FB, "The distance map must have at most 3 dimensions")
    dims = _pad(dm.shape, 3)
    s = _numbers(size, "size")
    if s.size != 3:
        raise StitchKernelError("crop:inputError", "Input range for bbox is not 3")
    for i in (1, 2, 0):  # the mex checks MATLAB's y, x, z
        if not (_floor(s[i]) >= dims[i] and s[i] < 9.0e15):
            raise StitchKernelError(_DR, _DR_SMALL)
    if 0 in dims and _floor(s[0]) * _floor(s[1]) * _floor(s[2]) > 0:
        raise StitchKernelError(_DR, "The distance map is empty")
    rdims = [_trunc_u64(v) for v in s]
    a = np.asarray(wd)
    d = np.float32(float(a.reshape(-1)[0]) if a.size else 0.0)  # (float) mxGetScalar(wd)
    if dm.dtype != np.float32:
        raise StitchKernelTypeError(_FB_TYPE, "Data type not suppported")
    out = np.zeros(rdims, dtype=np.float32)
    _ext().feather_distance_map_resize_3d(
        _flat(dm, "dmat"), _flat(out, "out"), dims, rdims, float(d), _threads(num_threads)
    )
    return out


def feather_distance_map_resize_3d(
    dmat: np.ndarray, box, wd, *, num_threads: int | None = None
) -> np.ndarray:
    """``feather_distance_map_resize_3d(dmat, bbox, wd)``: ``dmat`` resized to the size of
    ``box`` and raised to the power ``wd``.

    Single maps run the kernel. MATLAB's code for other classes (``imresize``/``imresize3`` and
    ``fastPower``) is not ported: they raise ``NotImplementedError``.
    """
    b = _ints(box, 6, "box")
    if _array(dmat, "dmat").dtype != np.float32:
        raise NotImplementedError(
            "feather_distance_map_resize_3d of a non-single distance map (MATLAB's imresize3 "
            "fall-back) is not ported"
        )
    return resize_distance_map(dmat, [b[3] - b[0], b[4] - b[1], b[5] - b[2]], wd,
                               num_threads=num_threads)  # fmt: skip


# ---------------------------------------------------------------------------------------------
# tile blending
# ---------------------------------------------------------------------------------------------

_TB = "tile_blending:inputError"


def tile_boxes_from_matlab(rows) -> np.ndarray:
    """MATLAB's regions ``[tile, ymin, xmin, zmin, ymax, xmax, zmax]`` (1-based, inclusive; a
    ``K x 7`` array or ``[]``) as rows ``(t, z0, y0, x0, z1, y1, x1)`` (0-based, end-exclusive),
    float64 (values that are not integers stay so, for ``tile_blend``'s check)."""
    m = np.asarray(rows if rows is not None else np.zeros((0, 7)), dtype=np.float64)
    if m.size == 0:
        return np.zeros((0, 7))
    if m.ndim == 1:
        m = m.reshape(1, -1)
    if m.ndim != 2 or m.shape[1] != 7:
        raise StitchKernelError(
            _TB, "The regions must be a K x 7 array (t, z0, y0, x0, z1, y1, x1)"
        )
    out = m[:, [0, 3, 1, 2, 6, 4, 5]].copy()
    out[:, :4] -= 1
    return out


def tile_boxes_to_matlab(rows) -> np.ndarray:
    """rows ``(t, z0, y0, x0, z1, y1, x1)`` as MATLAB's ``[tile, ymin, xmin, zmin, ymax, xmax,
    zmax]`` (float64)"""
    b = np.asarray(rows, dtype=np.float64).reshape(-1, 7)
    out = b[:, [0, 2, 3, 1, 5, 6, 4]].copy()
    out[:, :4] += 1
    return out


def _tile_boxes(boxes, dims: list[int]) -> list[int]:
    """the mex's checks of the regions, row by row: flat (t, z0, y0, x0, z1, y1, x1) of the
    regions that are not empty"""
    if boxes is None:
        return []
    b = np.asarray(boxes)
    if b.size == 0:
        return []
    if b.dtype.kind not in "biuf" or b.ndim != 2 or b.shape[1] != 7:
        raise StitchKernelError(
            _TB, "The regions must be a K x 7 array (t, z0, y0, x0, z1, y1, x1)"
        )
    v = b.astype(np.float64)
    out: list[int] = []
    for k, row in enumerate(v):
        for x in row:
            if not (x == _floor(x)) or abs(x) > 9.0e15:
                raise StitchKernelError(_TB, "The regions must have integer coordinates")
        if row[4] <= row[1] or row[5] <= row[2] or row[6] <= row[3]:
            continue
        if not (0 <= row[0] < dims[0] and row[1] >= 0 and row[2] >= 0 and row[3] >= 0):
            raise StitchKernelError(_TB, f"Region {k + 1} is outside the tiles")
        if row[4] > dims[1] or row[5] > dims[2] or row[6] > dims[3]:
            raise StitchKernelError(_TB, f"Region {k + 1} is outside the tiles")
        out.extend(int(x) for x in row)
    return out


def tile_blend(
    tiles: np.ndarray, boxes, blend_method: str, *, num_threads: int | None = None
) -> np.ndarray:
    """``tile_blending_3d_mex``: the tiles (uint8, uint16 or single, ``(t, z, y, x)``, 0 = no
    data) blended per voxel ('none': the first nonzero in tile order, 'max', 'mean', 'median'
    of the nonzero values) from the tiles outside their regions ``boxes`` (rows ``(t, z0, y0,
    x0, z1, y1, x1)``; empty ones are skipped), or where that is 0 from all tiles. 'none' and
    'max' give ``(z, y, x)`` of the tiles' class, 'mean' and 'median' single."""
    f = _array(tiles, "tiles")
    if f.dtype not in (np.uint8, np.uint16, np.float32):
        raise StitchKernelTypeError(
            "tile_blending:dataTypeError", "Data type not supported (uint8, uint16 or single)"
        )
    if _mndims(f) > 4:
        raise StitchKernelError(_TB, "The tiles must be a 4D array (tiles in the 4th dimension)")
    dims = _pad(f.shape, 4)
    if dims[0] == 0:
        raise StitchKernelError(_TB, "There are no tiles")
    if not isinstance(blend_method, str):
        raise StitchKernelError(_TB, "The blend method must be a char array")
    if len(blend_method.encode("utf-8")) >= 16:
        raise StitchKernelError(_TB, "Unsupported blending method")
    if blend_method not in ("none", "max", "mean", "median"):
        raise StitchKernelError(_TB, f"Unsupported blending method: {blend_method}")
    rows = _tile_boxes(boxes, dims)
    ocls = f.dtype if blend_method in ("none", "max") else np.dtype(np.float32)
    out = np.empty(dims[1:], dtype=ocls)
    getattr(_ext(), f"tile_blending_3d_{blend_method}")(
        _flat(f, "tiles"), _flat(out, "out"), dims, rows, _threads(num_threads)
    )
    return out


def tile_blending_3d(
    tim_f_block: np.ndarray, boxes, blend_method: str, *, num_threads: int | None = None
) -> np.ndarray:
    """``nv_block = tile_blending_3d(tim_f_block, mboxes, BlendMethod)``: the tiles of a
    stitching batch blended with 'none', 'mean', 'max' or 'median' (see ``tile_blend``; boxes in
    the port's form). uint8, uint16 and single tiles run the kernel; other classes run MATLAB's
    code, which blends in single and gives a single result: the kernel on ``single(tiles)``
    (the same operations)."""
    if not isinstance(blend_method, str) or blend_method not in BLEND_METHODS:
        raise ValueError(f"Unsupported blending method: {blend_method}.")
    f = _array(tim_f_block, "tim_f_block")
    if f.dtype in (np.uint8, np.uint16, np.float32):
        return tile_blend(f, boxes, blend_method, num_threads=num_threads)
    if f.dtype.kind not in "biuf":
        raise TypeError(f"tile blending of {f.dtype} tiles")
    from petakit5d.core import dtypes

    fs = np.ascontiguousarray(dtypes.single(f))
    return tile_blend(fs, boxes, blend_method, num_threads=num_threads)


# ---------------------------------------------------------------------------------------------
# integral image
# ---------------------------------------------------------------------------------------------

_II = "integral_image_3d:inputError"


def integral_image(a: np.ndarray, size_t, *, num_threads: int | None = None) -> np.ndarray:
    """``integral_image_3d_mex(A, szT)``: the sums of ``a`` (uint16 or single, ``(z, y, x)``)
    over windows of size ``size_t`` ``(tz, ty, tx)`` (1 <= size_t <= the image's size;
    truncated), single, of size ``a.shape + size_t - 1``: running sums along y, then
    differences along x and z, in single precision, as the mex computes them."""
    a = _array(a, "A")
    if _mndims(a) > 3:
        raise StitchKernelError(_II, "The input must have at most 3 dimensions")
    dims = _pad(a.shape, 3)
    t = _numbers(size_t, "size_t")
    if t.size != 3:
        raise StitchKernelError(_II, "The number of elements of size T must be 3")
    for i in (1, 2, 0):  # the mex checks MATLAB's y, x, z
        if not (_floor(t[i]) >= 1):
            raise StitchKernelError(_II, "The elements of size T must be greater than zeros")
        if not (_floor(t[i]) <= dims[i]):
            raise StitchKernelError(_II, "The size of the input must not be smaller than size T")
    dims_t = [_trunc_u64(x) for x in t]
    if a.dtype not in (np.uint16, np.float32):
        raise StitchKernelTypeError("integral_image_3d:dataTypeError", "Data type not suppported")
    out = np.zeros([dims[i] + dims_t[i] - 1 for i in range(3)], dtype=np.float32)
    _ext().integral_image_3d(_flat(a, "A"), _flat(out, "out"), dims, dims_t, _threads(num_threads))
    return out


def _cumsum(b: np.ndarray, axis: int) -> np.ndarray:
    """MATLAB's cumsum along an axis, in the class (logical: double; uint8 saturating)"""
    if b.dtype.kind == "b":
        return np.cumsum(b, axis=axis, dtype=np.float64)
    if b.dtype.kind == "f":
        return np.cumsum(b, axis=axis)
    from petakit5d.core import dtypes

    # nonnegative steps: the running sum saturates once, at the maximum
    c = np.cumsum(b, axis=axis, dtype=np.float64)
    return dtypes.cast(np.minimum(c, np.iinfo(b.dtype).max), "uint8")


def _integral_image_matlab(a: np.ndarray, t: Sequence[int]) -> np.ndarray:
    """integralImage of integral_image_3d.m (adapted from MATLAB's normxcorr2): zero padding,
    cumsum along y, then x, then z, and differences of the shifted sums, in the class of A"""
    from petakit5d.core import dtypes

    # indexing3d_mex and crop3d_mex take logical, uint8, uint16, single and double
    if a.dtype not in (np.bool_, np.uint8, np.float64):
        raise StitchKernelTypeError("tiff:dataTypeError", "Data type not suppported")
    s = _pad(a.shape, 3)
    b = np.zeros([s[i] + 2 * t[i] - 1 for i in range(3)], dtype=a.dtype)
    b[t[0] : t[0] + s[0], t[1] : t[1] + s[1], t[2] : t[2] + s[2]] = a.reshape(s)
    c = b
    for ax in (1, 2, 0):  # MATLAB's dimensions 1, 2, 3
        cs = _cumsum(c, ax)
        n = cs.shape[ax]
        hi = np.take(cs, range(t[ax], n), axis=ax)
        lo = np.take(cs, range(n - t[ax]), axis=ax)
        c = dtypes.minus(hi, lo) if cs.dtype.kind in "iu" else hi - lo
    return np.ascontiguousarray(c)


def integral_image_3d(a: np.ndarray, size_t, *, num_threads: int | None = None) -> np.ndarray:
    """``IntA = integral_image_3d(A, szT)``: the integral image of ``a`` over windows ``size_t``
    ``(tz, ty, tx)``. uint16 and single run the kernel (a single result); double, uint8 and
    logical run MATLAB's code in their class (double for logical)."""
    a = _array(a, "A")
    t = _numbers(size_t, "size_t")
    s = _pad(a.shape, 3)
    if t.size != 3:
        raise ValueError("size_t must have 3 values")
    if any(s[i] < t[i] for i in range(3)):
        raise ValueError("The size of A must be not smaller than size T!")
    if a.dtype in (np.uint16, np.float32):
        return integral_image(a, t, num_threads=num_threads)
    if _mndims(a) > 3 or any(x < 1 or x != int(x) for x in t):
        raise StitchKernelError(_II, "size_t must be 3 integers of at least 1 for a 3D image")
    return _integral_image_matlab(a, [int(x) for x in t])


# ---------------------------------------------------------------------------------------------
# the major-tile check
# ---------------------------------------------------------------------------------------------

_A4 = "any_4th_dim_mex:inputError"
_A4_TYPE = "any_4th_dim_mex:dataTypeError"


def major_tiles_valid(
    tiles: np.ndarray, major, *, num_threads: int | None = None
) -> tuple[bool, tuple | None]:
    """``[is_valid, unvalid_bbox] = any_4th_dim_mex(fmat, major_inds)``: whether every voxel of
    the tiles (uint16 or single, ``(t, z, y, x)``) has data (a nonzero value; NaN counts) in one
    of the major tiles (``major``, bool per tile, at least one true), else the box of the voxels
    without. Returns ``(True, None)`` or ``(False, box)``."""
    f = _array(tiles, "tiles")
    if _mndims(f) > 4:
        raise StitchKernelError(_A4, "The first input must have at most 4 dimensions")
    dims = _pad(f.shape, 4)
    m = np.asarray(major)
    if m.size != dims[0]:
        raise StitchKernelTypeError(
            _A4_TYPE,
            "The second input must be a logical array with same number element of the 4-th "
            "dimension of the first input",
        )
    if m.dtype != np.bool_:
        raise StitchKernelTypeError(_A4_TYPE, "The second input must be logical type")
    flags = np.ascontiguousarray(m.reshape(-1))
    if not flags.any():
        raise StitchKernelError(
            "any_4th_dim_mex:dataRangeError",
            "All elements in the second input are false, no indices to check!",
        )
    if f.dtype not in (np.uint16, np.float32):
        raise StitchKernelTypeError(_A4_TYPE, "Data type not suppported")
    valid, box = _ext().any_4th_dim(_flat(f, "tiles"), flags, dims, _threads(num_threads))
    return (True, None) if valid else (False, tuple(int(v) for v in box))


def check_major_tile_valid(
    fmat: np.ndarray, major_inds, *, num_threads: int | None = None
) -> tuple[bool, tuple | None]:
    """``[is_valid, unvalid_bbox] = check_major_tile_valid(fmat, major_inds)``: whether every
    voxel has data in one of the major tiles, else the box of the voxels without (``None``,
    MATLAB's ``[]``, when valid). uint16 and single tiles with a major tile run the kernel (NaN
    counts as data); other classes run MATLAB's code, where NaN is no data (``any`` ignores
    NaN)."""
    f = _array(fmat, "fmat")
    m = np.asarray(major_inds)
    if m.dtype != np.bool_:
        raise TypeError("major_inds must be bool (one flag per tile)")
    if f.dtype in (np.uint16, np.float32) and m.any():
        return major_tiles_valid(f, m, num_threads=num_threads)
    f4 = f.reshape(_pad(f.shape, 4))
    sel = f4[m.reshape(-1)]
    has = sel != 0
    if sel.dtype.kind in "fc":
        has &= ~np.isnan(sel)
    covered = has.any(axis=0)
    if covered.all():
        return True, None
    z, y, x = np.nonzero(~covered)
    return False, (int(z.min()), int(y.min()), int(x.min()),
                   int(z.max()) + 1, int(y.max()) + 1, int(x.max()) + 1)  # fmt: skip


# ---------------------------------------------------------------------------------------------
# NaN replacement
# ---------------------------------------------------------------------------------------------


def _value_scalar(value) -> float | None:
    """mxGetScalar of the value MATLAB receives: its first element as a double; None for values
    the mex rejects (not numeric or logical, complex, empty)"""
    if isinstance(value, np.ndarray | np.generic):
        a = np.asarray(value)
    elif value is None or isinstance(value, str):
        return None
    elif isinstance(value, bool):
        a = np.asarray(value, dtype=np.bool_)
    elif isinstance(value, int | float):
        a = np.asarray(value, dtype=np.float64)
    elif isinstance(value, list | tuple) and all(
        isinstance(x, int | float) and not isinstance(x, bool) for x in value
    ):
        a = np.asarray(value, dtype=np.float64)
    else:
        return None
    if a.dtype.kind not in "biuf" or a.size == 0:
        return None
    return float(a.reshape(-1)[0])  # the first element (MATLAB's order and NumPy's agree)


def replace_nan_with_value(
    x: np.ndarray, value: Any, *, inplace: bool = False, num_threads: int | None = None
) -> np.ndarray:
    """``replace_nan_with_value(x, value)`` (pointDetection's mex): ``x`` (single or double, any
    shape) with NaN replaced by ``value`` in ``x``'s class: its first element as a double
    (``mxGetScalar``), rounded to single for single ``x``. ``inplace`` writes into ``x``."""
    x = _array(x, "x")
    if x.dtype not in (np.float32, np.float64):
        raise StitchKernelTypeError(
            "data:TypeError", "Only single and double types are supported!\n"
        )
    dval = _value_scalar(value)
    if dval is None:
        raise StitchKernelError("data:inputError", "The value must be a full real numeric scalar.")
    val = x.dtype.type(dval)  # (float) dval for single input
    if not (x.flags.c_contiguous or x.flags.f_contiguous):
        raise ValueError("x must be a contiguous array")
    if inplace:
        if not x.flags.writeable:
            raise ValueError("x must be writeable")
        out = x
    else:
        out = np.empty_like(x)  # the layout of x (element-wise)
    # 1-D views in memory order (x and out have the same layout)
    xf, of = x.ravel(order="K"), out.ravel(order="K")
    if out.size and not np.shares_memory(of, out):
        raise ValueError("out could not be viewed without a copy")
    _ext().replace_nan_with_value(xf, of, float(val), _threads(num_threads))
    return out
