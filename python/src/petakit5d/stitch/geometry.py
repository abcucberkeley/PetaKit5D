"""Stitching geometry: axis and data orders, the voxel sizes of the stitching space, stage
coordinates in stitching space, tile cuboids and their overlaps, the half-overlap split, the
canvas (origin, padding, size, output box), tile placements, the batch grid, the size estimate
without registration and the major-tile cover check.

Ports ``stitch_process_coordinates``, ``check_cuboids_overlaps``,
``compute_half_of_overlap_region``, ``stitch_volume_size_estimate``, ``check_major_tile_cover``
and the geometry of ``XR_stitching_frame_zarr_dev_v1`` (axis and data orders, voxel sizes,
cuboids, overlaps, half overlap regions, canvas, batches) and of ``stitch_process_block_info``
(the tiles' boxes in the canvas). ``check_major_tile_valid`` is K04's
(``petakit5d.kernels.stitch``), re-exported here.

Conventions (``docs/conventions.md``):

- Stage and physical coordinates keep their names, x, y, z: ``coordinates`` and ``xyz`` are
  ``(n, 3)`` arrays with columns x, y, z (um), cuboids are ``(n, 6)`` rows
  ``[min_x min_y min_z max_x max_y max_z]`` (closed, physical), and the per-tile voxel offsets
  ``int_xyz_shift`` are columns x, y, z, as in MATLAB's stitch information files.
- ``axis_order`` strings (``'x,y,z'``, ``'-x,y,z'``) are unchanged.
- ``data_order`` names the axes of the tile arrays in their order: ``'z,y,x'`` for ``(z, y, x)``
  arrays (MATLAB's ``'y,x,z'``; ``core.axes.axis_string_from_matlab`` converts). Its parsed form,
  ``data_order_axes``, is the array axis of x, y and z (``(2, 1, 0)`` for ``'z,y,x'``; MATLAB's
  ``data_order_mat`` is the dimension of each, ``[2 1 3]``).
- Sizes (tile sizes, the stitched size) and per-axis vectors (batch, block and pad sizes, flags)
  are ``(z, y, x)`` in the order of the arrays' axes; ``resample_factor`` is ``(z, y, x)``
  (MATLAB's ``[x y z]``).
- Boxes are ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive. A half-overlap region is a box
  over a tile's voxels in ``(z, y, x)`` physical axes (MATLAB's ``3 x 2`` ``[s t]`` rows x, y,
  z, 1-based, inclusive), which the blending code reads with the standard data order.
- Tiles are numbered from 0; a pair is ``(i, j)``, ``i < j``.

Each value is computed with MATLAB's operations in MATLAB's order: ``round`` half away from zero,
``min``/``max`` along a dimension skip NaN, ``max`` and ``median`` of a single tile's ``1 x 3``
row reduce the row (a scalar, as in MATLAB), and ``sin``, ``cos``, ``tan`` and ``sind`` are
``matlabcompat.elementary``'s (fdlibm's: MATLAB's bits on every C library).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from petakit5d.matlabcompat import elementary

__all__ = [
    "GeometryError",
    "OVERLAP_TYPES",
    "StitchSpace",
    "StitchLayout",
    "StitchPlacement",
    "check_cuboids_overlaps",
    "check_major_tile_cover",
    "check_major_tile_valid",
    "compute_half_of_overlap_region",
    "batch_grid",
    "data_order_axes_from_matlab",
    "data_order_axes_to_matlab",
    "find_overlaps",
    "half_order_from_matlab",
    "half_order_to_matlab",
    "half_overlap_regions",
    "overlap_map_to_matlab",
    "pair_index",
    "parse_axis_order",
    "parse_data_order",
    "place_tiles",
    "prepare_layout",
    "process_coordinates",
    "region_from_matlab",
    "region_to_matlab",
    "resolve_overlap_type",
    "sort_tiles",
    "stage_to_xyz",
    "stitch_2d_from_input_sizes",
    "stitch_canvas",
    "stitch_pixel_sizes",
    "stitch_volume_size_estimate",
    "tile_canvas_boxes",
    "tile_canvas_boxes_to_matlab",
    "tile_counts",
    "tile_cuboids",
]

OVERLAP_TYPES = {"zero": 0.0, "none": 0.5, "half": 0.75, "full": 1.0}

_XYZ = "xyz"


class GeometryError(ValueError):
    """An error MATLAB raises; ``identifier`` is MATLAB's error id (``''`` for none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


def check_major_tile_valid(fmat, major_inds, *, num_threads=None):
    """K04's ``check_major_tile_valid`` (``petakit5d.kernels.stitch``)."""
    from petakit5d.kernels.stitch import check_major_tile_valid as f

    return f(fmat, major_inds, num_threads=num_threads)


# ---------------------------------------------------------------------------------------------
# MATLAB's elementary functions and reductions
# ---------------------------------------------------------------------------------------------


def _sin(x: float) -> float:
    return float(elementary.sin(x))


def _cos(x: float) -> float:
    return float(elementary.cos(x))


def _tan(x: float) -> float:
    return float(elementary.tan(x))


def _cot(x: float) -> float:
    """MATLAB's ``cot``: ``1 / tan(x)``."""
    return float(elementary.cot(x))


def _sind(x: float) -> float:
    """MATLAB's ``sind`` (``matlabcompat.elementary``: S03's measured rule)."""
    return float(elementary.sind(x))


def _mround(a):
    from petakit5d.core.dtypes import mround

    return mround(np.asarray(a, dtype=np.float64))


def _col_min(a: np.ndarray) -> np.ndarray:
    """``min(a, [], 1)``: per column, NaN skipped"""
    return np.fmin.reduce(a, axis=0)


def _max_default(a: np.ndarray) -> np.ndarray:
    """``max(a)``: per column of a matrix, the whole row of a single row (a scalar); NaN
    skipped"""
    if a.shape[0] == 1:
        return np.fmax.reduce(a[0])
    return np.fmax.reduce(a, axis=0)


def _median_default(a: np.ndarray) -> np.ndarray | np.float64:
    """``median(a)`` of sizes: per column of a matrix, of the row of a single row (a scalar);
    NaN gives NaN; an even count the mean of the two middle values"""

    def med(v: np.ndarray) -> np.float64:
        if np.isnan(v).any():
            return np.float64(np.nan)
        s = np.sort(v, kind="stable")
        k = s.size
        if k % 2:
            return s[k // 2]
        return (s[k // 2 - 1] + s[k // 2]) / 2

    if a.shape[0] == 1:
        return med(a[0])
    return np.array([med(a[:, c]) for c in range(a.shape[1])])


def _rows(a, ncol: int, name: str) -> np.ndarray:
    m = np.asarray(a, dtype=np.float64)
    if m.ndim == 1:
        m = m.reshape(1, -1)
    if m.ndim != 2 or m.shape[1] != ncol:
        raise GeometryError(f"{name} must be an (n, {ncol}) array, got shape {m.shape}")
    return m


# ---------------------------------------------------------------------------------------------
# axis and data orders
# ---------------------------------------------------------------------------------------------


def _letters(s: str) -> str:
    return s.replace(",", "").replace("-", "")


def parse_axis_order(axis_order: str) -> tuple[tuple[int, int, int], tuple[float, float, float]]:
    """MATLAB's ``order_sign_mat`` of an ``axisOrder`` string (``'x,y,z'``, ``'-x,y,z'``,
    ``'y,x,z'``): for x, y and z, the stage coordinate column (0-based; MATLAB's first row minus
    1) and its sign (``-1`` where the string has ``'-x'``, ``'-y'``, ``'-z'``; case-sensitive).

    Spaces are removed first. MATLAB's pattern check of the string never raises (``if
    ~regexpi(...)`` is false for a match and for no match), so neither does this: the columns
    are the sort order of its letters."""
    s = axis_order.replace(" ", "")
    pure = _letters(s)
    if len(pure) != 3:
        raise GeometryError(
            "Unable to perform assignment because the size of the left side is 1-by-3 and the "
            f"size of the right side is 1-by-{len(pure)}.",
            "MATLAB:subsassigndimmismatch",
        )
    perm = tuple(sorted(range(3), key=lambda k: pure[k]))
    signs = tuple(-1.0 if f"-{c}" in s else 1.0 for c in _XYZ)
    return perm, signs  # type: ignore[return-value]


def parse_data_order(data_order: str = "z,y,x") -> tuple[int, int, int]:
    """The array axis of x, y and z for tiles whose axes are named ``data_order`` in their
    order (``'z,y,x'``: ``(2, 1, 0)``; ``'z,x,y'``: ``(1, 2, 0)``), as MATLAB's ``[~,
    data_order_mat] = sort(dataOrder_pure)`` does (``,`` and ``-`` removed)."""
    pure = _letters(data_order)
    if len(pure) != 3:
        raise GeometryError(f"data_order names 3 axes, got {data_order!r}")
    return tuple(sorted(range(3), key=lambda k: pure[k]))  # type: ignore[return-value]


def data_order_axes_to_matlab(axes) -> tuple[int, int, int]:
    """``data_order_axes`` (the array axes of x, y, z) as MATLAB's ``data_order_mat`` (their
    dimensions: axis 1 is y, dimension 1; axis 2 x, 2; axis 0 z, 3)."""
    return tuple({1: 1, 2: 2, 0: 3}[int(a)] for a in axes)  # type: ignore[return-value]


def data_order_axes_from_matlab(dims) -> tuple[int, int, int]:
    """MATLAB's ``data_order_mat`` as ``data_order_axes``."""
    return tuple({1: 1, 2: 2, 3: 0}[int(d)] for d in np.asarray(dims).ravel())  # type: ignore


def _reverse(dop: Sequence[int]) -> tuple[int, int, int]:
    """for each array axis, the position (x 0, y 1, z 2) it holds (MATLAB's
    ``data_order_reverse_mat``)"""
    rev = [0, 0, 0]
    for k, a in enumerate(dop):
        rev[a] = k
    return tuple(rev)  # type: ignore[return-value]


def half_order_from_matlab(v) -> tuple[int, ...]:
    """MATLAB's ``halfOrder`` (positions in x, y, z: 1 x, 2 y, 3 z; ``[3, 2, 1]`` splits z
    first) as axes of ``(z, y, x)``: ``(0, 1, 2)``."""
    return tuple(3 - int(k) for k in np.asarray(v).ravel())


def half_order_to_matlab(v) -> tuple[int, ...]:
    """``half_order`` (axes of ``(z, y, x)``) as MATLAB's ``halfOrder``."""
    return tuple(3 - int(k) for k in v)


def stage_to_xyz(coordinates, axis_order: str = "x,y,z") -> np.ndarray:
    """``xyz = coordinates(:, order_sign_mat(1, :)) .* order_sign_mat(2, :); xyz = xyz -
    min(xyz, [], 1)``: stage coordinates (columns x, y, z, um) in the axis order, from 0."""
    c = _rows(coordinates, 3, "coordinates")
    perm, signs = parse_axis_order(axis_order)
    xyz = c[:, list(perm)] * np.asarray(signs)
    return xyz - _col_min(xyz) if xyz.shape[0] else xyz


def sort_tiles(tile_idx) -> tuple[np.ndarray, np.ndarray]:
    """``[tileIdx, sinds] = sortrows(tileIdx, [4, 3, 2, 1])``: tiles by location, then z, y and
    x index (stable, ascending, NaN last); ``sinds`` 0-based."""
    t = _rows(tile_idx, 4, "tile_idx")
    order = np.lexsort((t[:, 0], t[:, 1], t[:, 2], t[:, 3]))
    return t[order], order.astype(np.intp)


def tile_counts(tile_idx, axis_order: str = "x,y,z") -> np.ndarray:
    """``tileNum``: the number of distinct indices of columns x, y, z of ``tile_idx`` (each NaN
    counts), in the axis order's columns."""
    t = np.asarray(tile_idx, dtype=np.float64)
    t = t.reshape(-1, t.shape[-1]) if t.size else np.zeros((0, 4))
    if t.shape[1] < 3:
        raise GeometryError(f"tile_idx must have 3 or 4 columns, got shape {t.shape}")
    counts = []
    for c in range(3):
        v = t[:, c]
        nan = np.isnan(v)
        counts.append(np.unique(v[~nan]).size + int(nan.sum()))
    perm, _ = parse_axis_order(axis_order)
    return np.asarray(counts, dtype=np.float64)[list(perm)]


def stitch_2d_from_input_sizes(input_sizes, data_order: str = "z,y,x") -> bool:
    """The frame's first 2D check on the input tiles: 2D unless a tile has more than one voxel
    along ``imSize(data_order_mat == 3)``. That is the dimension at the position of z's
    dimension in ``data_order_mat``, which is z only for orders whose permutation is its own
    inverse (MATLAB's quirk: for MATLAB's ``'y,z,x'`` and ``'z,x,y'`` it is another axis);
    sizes ``(n, 3)`` in ``(z, y, x)`` array order."""
    s = _rows(input_sizes, 3, "input_sizes")
    dop = parse_data_order(data_order)
    dom = data_order_axes_to_matlab(dop)
    k = [i for i in range(3) if dom[i] == 3]  # MATLAB dimension k + 1 of the size
    if len(k) != 1:
        raise GeometryError("data_order must name z once")
    ax = {1: 1, 2: 2, 3: 0}[k[0] + 1]
    return not bool(np.any(s[:, ax] > 1))


# ---------------------------------------------------------------------------------------------
# voxel sizes of the stitching space
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class StitchSpace:
    """The stitching space: voxel sizes ``xp``, ``yp``, ``zp`` (um), the resample factors
    ``xf``, ``yf``, ``zf``, the skew angle ``theta`` (rad), ``z_aniso`` and ``dsr`` (DSR and
    not an IO scan)."""

    xp: float
    yp: float
    zp: float
    xf: float
    yf: float
    zf: float
    theta: float
    z_aniso: float
    dsr: bool

    @property
    def pixel_sizes(self) -> np.ndarray:
        """``[xp, yp, zp]``"""
        return np.array([self.xp, self.yp, self.zp])


def stitch_pixel_sizes(
    xy_pixel_size: float = 0.108,
    dz: float = 0.5,
    skew_angle: float = 32.45,
    *,
    resample_type: str = "xy_isotropic",
    resample_factor=None,
    objective_scan: bool = False,
    io_scan: bool = False,
    dsr: bool = False,
) -> StitchSpace:
    """The voxel sizes of the stitching space (``XR_stitching_frame_zarr_dev_v1``,
    ``stitch_volume_size_estimate``): ``resample_factor`` (``(z, y, x)``; given, it wins over
    ``resample_type``), ``'method_1'``, ``'isotropic'`` or ``'xy_isotropic'``; outside DSR
    space (also for IO scans) z is scaled by ``z_aniso``."""
    dsr = bool(dsr) and not bool(io_scan)
    if objective_scan or io_scan:
        z_aniso = dz / xy_pixel_size
    else:
        z_aniso = _sind(float(skew_angle)) * dz / xy_pixel_size
    theta = skew_angle * math.pi / 180
    rtype = resample_type
    if resample_factor is not None and np.size(resample_factor) > 0:
        rtype = "given"
    if rtype == "given":
        rf = np.asarray(resample_factor, dtype=np.float64).ravel()
        if rf.size < 3:
            raise GeometryError(
                f"Index exceeds the number of array elements. Index must not exceed {rf.size}.",
                "MATLAB:badsubscript",
            )
        # MATLAB's [x y z] is (z, y, x)
        xf, yf, zf = float(rf[rf.size - 1]), float(rf[rf.size - 2]), float(rf[rf.size - 3])
    elif rtype == "method_1":
        zf = _cot(abs(theta))
        yf = 1.0
        xf = 1 / _cos(abs(theta))
    elif rtype == "isotropic":
        zf = yf = xf = 1.0
    elif rtype == "xy_isotropic":
        s, c = _sin(theta), _cos(theta)
        zf = math.sqrt(
            (s * s + z_aniso * z_aniso * (c * c)) / (c * c + z_aniso * z_aniso * (s * s))
        )
        yf = xf = 1.0
    else:
        raise GeometryError("Unrecognized function or variable 'xf'.", "MATLAB:UndefinedFunction")
    if not dsr:
        zf = z_aniso * zf
    return StitchSpace(
        xy_pixel_size * xf,
        xy_pixel_size * yf,
        xy_pixel_size * zf,
        xf,
        yf,
        zf,
        theta,
        z_aniso,
        dsr,
    )


# ---------------------------------------------------------------------------------------------
# coordinates in stitching space
# ---------------------------------------------------------------------------------------------


def process_coordinates(
    xyz,
    im_sizes,
    flipped_tile,
    axis_order: str,
    data_order: str,
    space: StitchSpace,
    *,
    io_scan: bool = False,
    objective_scan: bool = False,
    ds: bool = False,
) -> np.ndarray:
    """``stitch_process_coordinates``: the x of flipped tiles moved by their extent, stage
    coordinates in skewed space (not DSR, DS, objective or IO scan) or DS/objective-scan space,
    then snapped to whole voxels from the minimum (``round``). ``im_sizes`` ``(n, 3)`` in array
    order; ``space.dsr`` is the DSR flag (False for IO scans)."""
    xyz = _rows(xyz, 3, "xyz").copy()
    s = _rows(im_sizes, 3, "im_sizes")
    dop = parse_data_order(data_order)
    _, signs = parse_axis_order(axis_order)
    xp, zp, theta = space.xp, space.zp, space.theta
    dsr = space.dsr
    if flipped_tile is not None and np.size(flipped_tile) > 0:
        flipped = (np.asarray(flipped_tile, dtype=np.float64).ravel() > 0).astype(np.float64)
        if dsr:
            offsize = (s[:, dop[0]] - 1) * xp - s[:, dop[2]] * _cot(theta) * zp
        else:
            offsize = (s[:, dop[2]] - 1) * zp / _sin(theta)
        xyz[:, 0] = xyz[:, 0] - offsize * flipped * signs[0]
    if not io_scan and not objective_scan and not ds and not dsr:
        st, ct = _sin(theta), _cos(theta)
        xyz = np.column_stack((xyz[:, 2] / st, xyz[:, 1], -xyz[:, 0] * st + xyz[:, 2] * ct))
    if objective_scan or ds:
        st, ct = _sin(theta), _cos(theta)
        xyz = np.column_stack(
            (xyz[:, 0] * ct + xyz[:, 2] * st, xyz[:, 1], -xyz[:, 0] * st + xyz[:, 2] * ct)
        )
    p = space.pixel_sizes
    if xyz.shape[0] == 0:
        return xyz
    return _mround((xyz - _col_min(xyz)) / p) * p


def tile_cuboids(xyz, im_sizes, data_order: str, space: StitchSpace) -> np.ndarray:
    """``[xyz, xyz + (imSizes(:, data_order_mat) - 1) .* [xp, yp, zp]]``: each tile's closed
    cuboid in stitching space, ``[min_x min_y min_z max_x max_y max_z]``."""
    xyz = _rows(xyz, 3, "xyz")
    s = _rows(im_sizes, 3, "im_sizes")
    dop = parse_data_order(data_order)
    return np.hstack((xyz, xyz + (s[:, list(dop)] - 1) * space.pixel_sizes))


# ---------------------------------------------------------------------------------------------
# overlaps
# ---------------------------------------------------------------------------------------------


def check_cuboids_overlaps(
    cuboid_1, cuboid_2, stitch_2d: bool = False, *, regions: bool = True
) -> tuple[np.ndarray, np.ndarray | None]:
    """``[is_overlap_mat, cuboid_overlap_mat] = check_cuboids_overlaps(cuboid_1, cuboid_2,
    stitch2D)``: whether closed cuboids ``[min1 min2 min3 max1 max2 max3]`` (rows; one row is
    paired with every row of the other set) overlap, and with ``regions`` (MATLAB's second
    output) their overlaps, ``None`` (MATLAB's ``[]``) when none does.

    With ``regions``, an overlap that is a single point along an axis (3D; along the first two
    axes in 2D) counts as none, with a zero row; without it (MATLAB's one output) it counts."""
    c1 = np.asarray(cuboid_1, dtype=np.float64)
    c2 = np.asarray(cuboid_2, dtype=np.float64)
    c1 = c1.reshape(1, -1) if c1.ndim == 1 else c1
    c2 = c2.reshape(1, -1) if c2.ndim == 1 else c2
    if c1.ndim != 2 or c2.ndim != 2 or c1.shape[1] != 6 or c2.shape[1] != 6:
        raise GeometryError("The size of the cuboids must be n X 6!")
    n1, n2 = c1.shape[0], c2.shape[0]
    msg = "The number of cuboids in the two sets must equal or one set need to have one cuboid!"
    if n1 < n2:
        if n1 != 1:
            raise GeometryError(msg)
        n = n2
        c1 = np.repeat(c1, n, axis=0)
    elif n1 > n2:
        if n2 != 1:
            raise GeometryError(msg)
        n = n1
        c2 = np.repeat(c2, n, axis=0)
    else:
        n = n1
    nol = (c2[:, 3:] < c1[:, :3]) | (c1[:, 3:] < c2[:, :3])
    is_overlap = ~nol.any(axis=1)
    if not regions or not is_overlap.any():
        return is_overlap, None
    out = np.zeros((n, 6))
    for p in np.nonzero(is_overlap)[0]:
        # per axis the middle two of the four bounds (MATLAB's sort: stable, NaN last)
        b = np.sort(
            np.column_stack((c1[p].reshape(2, 3).T, c2[p].reshape(2, 3).T)), axis=1, kind="stable"
        )[:, 1:3]
        out[p] = b.T.reshape(-1)
        k = 2 if stitch_2d else 3
        if np.any(b[:k, 0] == b[:k, 1]):
            out[p] = 0
            is_overlap[p] = False
    return is_overlap, out


def pair_index(i: int, j: int, n: int) -> int:
    """The row of pair ``(i, j)`` (0-based, ``i < j``) in ``overlap_regions``: MATLAB's
    ``(i - 1) * nF - i * (i + 1) / 2 + j`` minus 1."""
    return i * n - (i + 1) * (i + 2) // 2 + j


def find_overlaps(cuboids, stitch_2d: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """The frame's first overlap pass: ``overlap_matrix`` (``(n, n)`` bool, ``[i, j]`` for
    ``i < j``) and ``overlap_regions`` (``(n (n - 1) / 2, 6)``, the overlap cuboid of pair
    ``(i, j)`` in row ``pair_index(i, j, n)``, zeros where none)."""
    c = _rows(cuboids, 6, "cuboids")
    n = c.shape[0]
    om = np.zeros((n, n), dtype=bool)
    regions = np.zeros((n * (n - 1) // 2, 6))
    for i in range(n - 1):
        is_ov, ov = check_cuboids_overlaps(c[i], c[i + 1 :], stitch_2d)
        for k in np.nonzero(is_ov)[0]:
            j = i + 1 + int(k)
            om[i, j] = True
            regions[pair_index(i, j, n)] = ov[k]
    return om, regions


def resolve_overlap_type(overlap_type: str = "", blend_method: str = "none") -> str:
    """``overlapType`` ``''``: ``'none'`` for ``'none'`` blending, else ``'half'``."""
    if overlap_type == "":
        return "none" if blend_method == "none" else "half"
    return overlap_type


def compute_half_of_overlap_region(
    cuboid_1,
    cuboid_2,
    xyz_voxel_sizes,
    *,
    overlap_type: str = "half",
    half_order: Sequence[int] = (1, 0, 2),
    stitch_2d: bool = False,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """``[mregion_1, mregion_2] = compute_half_of_overlap_region(cuboid_1, cuboid_2,
    xyz_voxelsizes, ...)``: the overlap of two tile cuboids as a box over each tile's voxels
    (``(z0, y0, x0, z1, y1, x1)`` in physical axes), with ``round(len * f)`` voxels taken off
    along the split axis (``f``: zero 0, none 0.5, half 0.75, full 1), at the end of the first
    tile and the start of the second unless the first tile ends beyond the second. The split
    axis is the first of ``half_order`` (axes of ``(z, y, x)``; MATLAB's default ``[2, 3, 1]``
    is ``(1, 0, 2)``) with more than one voxel of overlap (the first if none has).

    Cuboids that do not overlap at all raise (MATLAB indexes an empty overlap); a single-point
    overlap gives regions from a zero overlap cuboid, as in MATLAB. An unknown ``overlap_type``
    raises after that (MATLAB's ``fc`` is undefined)."""
    vs = [float(v) for v in np.asarray(xyz_voxel_sizes, dtype=np.float64).ravel()]
    c1 = [float(v) for v in np.asarray(cuboid_1, dtype=np.float64).ravel()]
    c2 = [float(v) for v in np.asarray(cuboid_2, dtype=np.float64).ravel()]
    if len(c1) != 6 or len(c2) != 6:
        raise GeometryError("The size of the cuboids must be n X 6!")
    ov = _pair_overlap(c1, c2, stitch_2d)
    if ov is None:
        # (MATLAB: cubiod_overlap is [], and indexing it errors)
        raise GeometryError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    s1 = [_round((ov[k] - c1[k]) / vs[k]) + 1 for k in range(3)]
    s2 = [_round((ov[k] - c2[k]) / vs[k]) + 1 for k in range(3)]
    t1 = [_round((ov[3 + k] - c1[k]) / vs[k]) + 1 for k in range(3)]
    t2 = [_round((ov[3 + k] - c2[k]) / vs[k]) + 1 for k in range(3)]
    ol_len = [t1[k] - s1[k] + 1 for k in range(3)]
    # positions x 0, y 1, z 2 of the axes of half_order
    order = [2 - int(a) for a in half_order]
    split = order[0]
    if ol_len[split] <= 1:
        if ol_len[order[1]] > 1:
            split = order[1]
        elif ol_len[order[2]] > 1:
            split = order[2]
    if overlap_type not in OVERLAP_TYPES:
        raise GeometryError("Unrecognized function or variable 'fc'.", "MATLAB:UndefinedFunction")
    rm = _round(ol_len[split] * OVERLAP_TYPES[overlap_type])
    e1, e2 = list(t1), list(t2)
    if t1[split] > t2[split]:
        s1[split] += rm
        e2[split] -= rm
    else:
        e1[split] -= rm
        s2[split] += rm
    return _region(s1, e1), _region(s2, e2)


def _round(v: float) -> float:
    """MATLAB's round of a double (half away from zero; the sign of a zero result is not kept)"""
    if not math.isfinite(v):
        return v
    t = float(math.trunc(v))
    if abs(v - t) >= 0.5:
        t += 1.0 if v > 0 else -1.0
    return t


def _region(s: list, t: list) -> tuple[int, ...]:
    """MATLAB's [s t] rows x, y, z as a box (z0, y0, x0, z1, y1, x1)"""
    return (int(s[2]) - 1, int(s[1]) - 1, int(s[0]) - 1, int(t[2]), int(t[1]), int(t[0]))


def _pair_overlap(c1: list, c2: list, stitch_2d: bool) -> list | None:
    """check_cuboids_overlaps of one pair with its overlap: None without one, zeros for a
    single-point overlap"""
    if any(math.isnan(v) for v in c1) or any(math.isnan(v) for v in c2):
        _, ov = check_cuboids_overlaps(c1, c2, stitch_2d)
        return None if ov is None else [float(v) for v in ov[0]]
    for k in range(3):
        if c2[3 + k] < c1[k] or c1[3 + k] < c2[k]:
            return None
    lo, hi = [0.0] * 3, [0.0] * 3
    for k in range(3):
        b = sorted((c1[k], c1[3 + k], c2[k], c2[3 + k]))  # stable, as MATLAB's sort
        lo[k], hi[k] = b[1], b[2]
    m = 2 if stitch_2d else 3
    if any(lo[k] == hi[k] for k in range(m)):
        return [0.0] * 6
    return lo + hi


def region_from_matlab(m) -> tuple[int, ...]:
    """MATLAB's region ``[s t]`` (``3 x 2``, rows x, y, z, 1-based, inclusive) as a box
    ``(z0, y0, x0, z1, y1, x1)``."""
    r = np.asarray(m, dtype=np.float64).reshape(3, 2)
    return (
        int(r[2, 0]) - 1,
        int(r[1, 0]) - 1,
        int(r[0, 0]) - 1,
        int(r[2, 1]),
        int(r[1, 1]),
        int(r[0, 1]),
    )


def region_to_matlab(box) -> np.ndarray:
    """A half-overlap box as MATLAB's ``3 x 2`` region ``[s t]`` (rows x, y, z)."""
    z0, y0, x0, z1, y1, x1 = (int(v) for v in box)
    return np.array([[x0 + 1, x1], [y0 + 1, y1], [z0 + 1, z1]], dtype=np.float64)


def half_overlap_regions(
    cuboids,
    overlap_matrix,
    space: StitchSpace,
    *,
    overlap_type: str,
    half_order: Sequence[int] = (0, 1, 2),
    stitch_2d: bool = False,
) -> tuple[np.ndarray, list, list, np.ndarray]:
    """The frame's second overlap pass on the shifted cuboids: pairs of ``overlap_matrix`` that
    no longer overlap (MATLAB's one-output check: a single-point overlap still counts) are
    dropped; for the others the half-overlap regions (``overlap_type``) and the full ones
    (``'zero'``). Returns the new ``overlap_matrix``, the half and full regions (a pair of boxes
    per overlapping pair, in loop order) and the pairs ``(K, 2)`` ``(i, j)``."""
    c = _rows(cuboids, 6, "cuboids")
    om = np.array(overlap_matrix, dtype=bool, copy=True)
    n = c.shape[0]
    half, full, pairs = [], [], []
    vs = space.pixel_sizes
    for i in range(n - 1):
        js = np.nonzero(om[i, i + 1 :])[0] + i + 1
        if js.size == 0:
            continue
        is_ov, _ = check_cuboids_overlaps(c[i], c[js], stitch_2d, regions=False)
        for j1, j in enumerate(js):
            if not is_ov[j1]:
                om[i, j] = False
                continue
            half.append(
                compute_half_of_overlap_region(
                    c[i],
                    c[j],
                    vs,
                    overlap_type=overlap_type,
                    half_order=half_order,
                    stitch_2d=stitch_2d,
                )
            )
            full.append(
                compute_half_of_overlap_region(
                    c[i],
                    c[j],
                    vs,
                    overlap_type="zero",
                    half_order=half_order,
                    stitch_2d=stitch_2d,
                )
            )
            pairs.append((i, int(j)))
    return om, half, full, np.asarray(pairs, dtype=np.int64).reshape(-1, 2)


def overlap_map_to_matlab(pairs, n: int) -> np.ndarray:
    """``overlap_map_mat``: ``[counter, (j - 1) * nF + i]`` (1-based) of the pairs."""
    p = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    k = np.arange(1, p.shape[0] + 1, dtype=np.float64)
    return np.column_stack((k, (p[:, 1] * n + p[:, 0] + 1).astype(np.float64)))


# ---------------------------------------------------------------------------------------------
# canvas, tile boxes, batches
# ---------------------------------------------------------------------------------------------


def stitch_canvas(
    xyz_shift,
    im_sizes,
    data_order: str,
    space: StitchSpace,
    tile_num,
    *,
    xcorr_shift: bool = True,
    stitch_pad_size: Sequence[float] = (2, 2, 2),
    stitch_mip=None,
    out_bbox=None,
    p_im_sz=None,
) -> dict:
    """The stitched canvas: padding (``2`` voxels per tile along each axis by default, only with
    registration and more than one tile) on both sides, the origin, each tile's offset from it,
    the stitched size (``round(extent / voxel + 2)``, 1 along z for MIPs and 2D tiles; a
    secondary channel's ``p_im_sz``), an output box (``out_bbox``, open ends at the size) and the
    voxel offsets ``round(xyz_shift / voxel)``.

    Returns a dict: ``pad_size`` and ``origin`` (x, y, z, um), ``xyz_shift`` (relative,
    ``(n, 3)``), ``nv_size`` (``(z, y, x)`` array order), ``int_xyz_shift`` (``(n, 3)``, x, y,
    z voxels)."""
    xs = _rows(xyz_shift, 3, "xyz_shift")
    s = _rows(im_sizes, 3, "im_sizes")
    dop = parse_data_order(data_order)
    rev = _reverse(dop)
    p = space.pixel_sizes
    n = xs.shape[0]
    if n > 1 and xcorr_shift:
        pad = np.asarray(stitch_pad_size, dtype=np.float64).ravel()
        pad_size = p * pad[list(dop)] * np.asarray(tile_num, dtype=np.float64).ravel()
    else:
        pad_size = np.zeros(3)
    origin = _col_min(xs) - pad_size
    rel = xs - origin
    if p_im_sz is None:
        dxyz = _max_default(xs + s[:, list(dop)] * p) + pad_size - origin
        vxyz = _mround(dxyz / p + 2)
        nv = vxyz[list(rev)]
        if (stitch_mip is not None and np.any(stitch_mip)) or bool(np.any(s[:, dop[2]] == 1)):
            nv[dop[2]] = 1
    else:
        nv = np.asarray(p_im_sz, dtype=np.float64).ravel().copy()
    if out_bbox is not None and np.size(out_bbox) > 0:
        b = [None if v is None else float(v) for v in out_bbox]
        start = np.array(b[:3], dtype=np.float64)
        end = np.array([nv[a] if b[3 + a] is None else b[3 + a] for a in range(3)])
        pad_l = pad_size - start[list(dop)] * p
        origin = _col_min(xs) - pad_l
        rel = xs - origin
        nv = end - start
    return {
        "pad_size": pad_size,
        "origin": origin,
        "xyz_shift": rel,
        "nv_size": nv,
        "int_xyz_shift": _mround(rel / p),
    }


def tile_canvas_boxes(
    int_xyz_shift, im_sizes, nv_size, data_order: str = "z,y,x", stitch_2d: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Each tile's box in the canvas and in the tile (``stitch_process_block_info``): the part
    of the tile at its offset that lies in the canvas. Returns ``(canvas_boxes, tile_boxes,
    valid)``, ``(n, 6)`` int64 boxes and a flag per tile; rows of tiles outside the canvas are
    zeros. In 2D both boxes are ``(0, 1)`` along z."""
    sh = _rows(int_xyz_shift, 3, "int_xyz_shift")
    s = _rows(im_sizes, 3, "im_sizes")
    nv = np.asarray(nv_size, dtype=np.float64).ravel()
    dop = parse_data_order(data_order)
    rev = _reverse(dop)
    n = sh.shape[0]
    cb = np.zeros((n, 6), dtype=np.int64)
    tb = np.zeros((n, 6), dtype=np.int64)
    valid = np.ones(n, dtype=bool)
    for i in range(n):
        st = sh[i, list(rev)]
        for a in range(3):
            # MATLAB: ridx = max(1, 1 - st) : min(imSize, nvSize - st); idx = st + ridx
            r0 = max(1.0, 1.0 - st[a])
            r1 = min(s[i, a], nv[a] - st[a])
            if stitch_2d and a == dop[2]:
                r0 = r1 = 1.0
                st_a = 0.0
            else:
                st_a = st[a]
            nr = math.floor(r1 - r0) + 1 if r1 >= r0 else 0
            if nr <= 0:
                cb[i] = 0
                tb[i] = 0
                valid[i] = False
                break
            last = r0 + (nr - 1)
            cb[i, a], cb[i, 3 + a] = int(st_a + r0) - 1, int(st_a + last)
            tb[i, a], tb[i, 3 + a] = int(r0) - 1, int(last)
    return cb, tb, valid


def tile_canvas_boxes_to_matlab(canvas_boxes, tile_boxes, valid) -> tuple[np.ndarray, ...]:
    """MATLAB's ``bboxes`` (``[ymin xmin zmin ymax xmax zmax]`` in the canvas), ``bboxStart_mat``
    and ``bboxEnd_mat`` (first and last voxel in the tile), zeros for tiles outside."""
    from petakit5d.core.bbox import boxes_to_matlab

    v = np.asarray(valid, dtype=bool)
    cb = boxes_to_matlab(np.asarray(canvas_boxes)).astype(np.float64)
    tb = boxes_to_matlab(np.asarray(tile_boxes)).astype(np.float64)
    cb[~v] = 0
    tb[~v] = 0
    return cb, tb[:, :3].copy(), tb[:, 3:].copy()


def batch_grid(nv_size, im_sizes, batch_size=(500, 500, 500), block_size=(500, 500, 500)):
    """The frame's batches: ``batch_size`` at most the stitched size and twice the median tile
    size (of the row for a single tile: a scalar), ``block_size`` at most the batch size, the
    batch size rounded up to whole blocks, the batch grid. All ``(z, y, x)``. Returns
    ``(batch_size, block_size, b_sub_sz, num_batches)``; ``num_batches`` is MATLAB's double
    (NaN, as in MATLAB, for a size with a zero)."""
    nv = np.asarray(nv_size, dtype=np.float64).ravel()
    s = _rows(im_sizes, 3, "im_sizes")
    with np.errstate(invalid="ignore", divide="ignore"):
        bs = np.fmin(nv, np.asarray(batch_size, dtype=np.float64))
        bs = np.fmin(_median_default(s) * 2, bs)
        bl = np.fmin(bs, np.asarray(block_size, dtype=np.float64))
        bs = np.ceil(bs / bl) * bl
        sub = np.ceil(nv / bs)
    return bs, bl, sub, float(np.prod(sub))


def _batch_number(nv, im_size, batch_size, block_size) -> float:
    """``batch_number`` of stitch_volume_size_estimate (twice the tile size, per axis)"""
    with np.errstate(invalid="ignore", divide="ignore"):
        bs = np.fmin(nv, np.asarray(batch_size, dtype=np.float64))
        bs = np.fmin(np.asarray(im_size, dtype=np.float64) * 2, bs)
        bl = np.fmin(bs, np.asarray(block_size, dtype=np.float64))
        bs = np.ceil(bs / bl) * bl
        return float(np.prod(np.ceil(nv / bs)))


# ---------------------------------------------------------------------------------------------
# the frame's geometry in two steps (before and after the registration)
# ---------------------------------------------------------------------------------------------


@dataclass
class StitchLayout:
    """The geometry before registration: ``xyz_orig`` (stage coordinates from 0), the tile
    order ``sinds`` (0-based; the sorted ``tile_idx``), ``tile_num`` (x, y, z), the space, the
    tile sizes ``im_sizes`` in that order, the snapped coordinates ``xyz``, the ``cuboids``,
    ``overlap_matrix`` and ``overlap_regions``."""

    data_order: str
    axis_order: str
    xyz_orig: np.ndarray
    sinds: np.ndarray
    tile_idx: np.ndarray
    tile_num: np.ndarray
    space: StitchSpace
    im_sizes: np.ndarray
    xyz: np.ndarray
    cuboids: np.ndarray
    overlap_matrix: np.ndarray
    overlap_regions: np.ndarray
    stitch_2d: bool


@dataclass
class StitchPlacement:
    """The geometry after registration (``xyz_shift``): the rechecked ``overlap_matrix``, the
    half and full overlap regions per pair (``pairs``), the canvas (``pad_size``, ``origin``,
    relative ``xyz_shift``, ``nv_size``, ``int_xyz_shift``), the batch grid and the tiles'
    boxes in the canvas."""

    overlap_type: str
    overlap_matrix: np.ndarray
    half_regions: list
    full_regions: list
    pairs: np.ndarray
    pad_size: np.ndarray
    origin: np.ndarray
    xyz_shift: np.ndarray
    nv_size: np.ndarray
    int_xyz_shift: np.ndarray
    batch_size: np.ndarray
    block_size: np.ndarray
    b_sub_sz: np.ndarray
    num_batches: float
    canvas_boxes: np.ndarray = field(repr=False)
    tile_boxes: np.ndarray = field(repr=False)
    tile_valid: np.ndarray = field(repr=False)


def prepare_layout(
    coordinates,
    im_sizes,
    *,
    tile_idx=None,
    flipped_tile=None,
    axis_order: str = "x,y,z",
    data_order: str = "z,y,x",
    xy_pixel_size: float = 0.108,
    dz: float = 0.5,
    skew_angle: float = 32.45,
    objective_scan: bool = False,
    io_scan: bool = False,
    ds: bool = False,
    dsr: bool = False,
    resample_type: str = "xy_isotropic",
    resample_factor=None,
    stitch_2d: bool = False,
    coord_sizes=None,
    xyz_orig=None,
) -> StitchLayout:
    """The frame's geometry up to the registration, for tiles listed with stage
    ``coordinates`` and ``im_sizes`` (``(n, 3)`` array order, in the listed order; the
    processed tiles' sizes): the axis order, the tile order by ``tile_idx`` (``(n, 4)`` x, y, z,
    location), the voxel sizes, the coordinates in stitching space (with ``coord_sizes``, a
    secondary's primary sizes, if given; ``xyz_orig`` replaces the stage coordinates from 0, as
    ``usePrimaryCoords`` does), cuboids and overlaps. ``flipped_tile`` keeps its listed order, as
    in MATLAB. ``stitch_2d`` is the frame's first check (``stitch_2d_from_input_sizes`` of the
    input tiles; None: False); it is also set when every processed tile is one voxel along z."""
    xyz0 = stage_to_xyz(coordinates, axis_order)
    xyz = xyz0 if xyz_orig is None else _rows(xyz_orig, 3, "xyz_orig")
    n = xyz.shape[0]
    if tile_idx is not None and np.size(tile_idx) > 0:
        tidx, sinds = sort_tiles(tile_idx)
        xyz = xyz[sinds]
    else:
        tidx, sinds = np.zeros((n, 4)), np.arange(n, dtype=np.intp)
    tnum = tile_counts(tidx, axis_order)
    space = stitch_pixel_sizes(
        xy_pixel_size,
        dz,
        skew_angle,
        resample_type=resample_type,
        resample_factor=resample_factor,
        objective_scan=objective_scan,
        io_scan=io_scan,
        dsr=dsr,
    )
    s = _rows(im_sizes, 3, "im_sizes")[sinds]
    dop = parse_data_order(data_order)
    stitch_2d = bool(stitch_2d) or bool(np.all(s[:, dop[2]] == 1))
    cs = s if coord_sizes is None else _rows(coord_sizes, 3, "coord_sizes")
    pxyz = process_coordinates(
        xyz,
        cs,
        flipped_tile,
        axis_order,
        data_order,
        space,
        io_scan=io_scan,
        objective_scan=objective_scan,
        ds=ds,
    )
    cub = tile_cuboids(pxyz, s, data_order, space)
    om, regions = find_overlaps(cub, stitch_2d)
    return StitchLayout(
        data_order,
        axis_order,
        xyz0,
        sinds,
        tidx,
        tnum,
        space,
        s,
        pxyz,
        cub,
        om,
        regions,
        stitch_2d,
    )


def place_tiles(
    layout: StitchLayout,
    xyz_shift=None,
    *,
    blend_method: str = "none",
    overlap_type: str = "",
    half_order: Sequence[int] = (0, 1, 2),
    xcorr_shift: bool = True,
    stitch_pad_size: Sequence[float] = (2, 2, 2),
    stitch_mip=None,
    out_bbox=None,
    p_im_sz=None,
    batch_size: Sequence[float] = (500, 500, 500),
    block_size: Sequence[float] = (500, 500, 500),
) -> StitchPlacement:
    """The frame's geometry after the registration: ``xyz_shift`` (``(n, 3)``, the shifted
    coordinates; None: ``layout.xyz``, as without registration). Overlaps rechecked, half-overlap
    regions, canvas, batch grid and tile boxes (see the functions of each step)."""
    s = layout.im_sizes
    xs = layout.xyz if xyz_shift is None else _rows(xyz_shift, 3, "xyz_shift")
    otype = resolve_overlap_type(overlap_type, blend_method)
    cub = tile_cuboids(xs, s, layout.data_order, layout.space)
    om, half, full, pairs = half_overlap_regions(
        cub,
        layout.overlap_matrix,
        layout.space,
        overlap_type=otype,
        half_order=half_order,
        stitch_2d=layout.stitch_2d,
    )
    cv = stitch_canvas(
        xs,
        s,
        layout.data_order,
        layout.space,
        layout.tile_num,
        xcorr_shift=xcorr_shift,
        stitch_pad_size=stitch_pad_size,
        stitch_mip=stitch_mip,
        out_bbox=out_bbox,
        p_im_sz=p_im_sz,
    )
    bs, bl, sub, nb = batch_grid(cv["nv_size"], s, batch_size, block_size)
    cb, tb, valid = tile_canvas_boxes(
        cv["int_xyz_shift"], s, cv["nv_size"], layout.data_order, layout.stitch_2d
    )
    return StitchPlacement(
        otype,
        om,
        half,
        full,
        pairs,
        cv["pad_size"],
        cv["origin"],
        cv["xyz_shift"],
        cv["nv_size"],
        cv["int_xyz_shift"],
        bs,
        bl,
        sub,
        nb,
        cb,
        tb,
        valid,
    )


# ---------------------------------------------------------------------------------------------
# the size estimate
# ---------------------------------------------------------------------------------------------


def stitch_volume_size_estimate(
    coordinates,
    im_size,
    *,
    tile_idx=None,
    flipped_tile=None,
    axis_order: str = "x,y,z",
    data_order: str = "z,y,x",
    xy_pixel_size: float = 0.108,
    dz: float = 0.5,
    skew_angle: float = 32.45,
    objective_scan: bool = False,
    io_scan: bool = False,
    ds: bool = False,
    dsr: bool = False,
    resample_type: str = "xy_isotropic",
    resample_factor=None,
    input_bbox=None,
    tile_out_bbox=None,
    xcorr_shift: bool = True,
    xy_max_offset=300,
    z_max_offset=50,
    stitch_pad_size: Sequence[float] = (2, 2, 2),
    out_bbox=None,
    stitch_mip=None,
    batch_size: Sequence[float] = (500, 500, 500),
    block_size: Sequence[float] = (500, 500, 500),
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """``[nvSize, numBatches, nvSizeRange, numBatchesRange] = stitch_volume_size_estimate(
    coordinates, imSize, ...)``: the stitched size (``(z, y, x)``) of tiles of size ``im_size``
    (as read, all alike; ``(ny, nx)`` is one slice) at the stage ``coordinates`` (``(n, 3)``,
    um) without registration shifts, its number of batches, and with registration
    (``xcorr_shift`` and more than one tile, no ``out_bbox``) the smallest and largest sizes
    (rows of ``nv_size_range``, ``(2, 3)``) and batch counts, the tiles moving by up to the max
    offsets at both ends.

    Boxes are Python boxes, sizes and per-axis values ``(z, y, x)``, ``resample_factor``
    ``(z, y, x)`` (MATLAB's ``[x y z]``; for a skewed-space resample MATLAB divides the tile
    size ``[ny nx nz]`` by ``[xf yf zf]`` element by element, so ny by xf and nx by yf: kept).
    The batches use twice the tile size per axis (the frame uses twice the median size, of the
    row for a single tile)."""
    c = _rows(coordinates, 3, "coordinates")
    nf = c.shape[0]
    sz = [float(v) for v in np.asarray(im_size, dtype=np.float64).ravel()]
    if len(sz) < 3:
        sz = [1.0] * (3 - len(sz)) + sz  # (ny, nx) is (1, ny, nx); (n,) (1, 1, n)
    if len(sz) > 3:
        raise GeometryError("im_size has 2 or 3 values")
    im = np.array(sz)
    dop = parse_data_order(data_order)
    rev = _reverse(dop)
    xyz = stage_to_xyz(c, axis_order)
    tidx = np.zeros((nf, 4)) if tile_idx is None or np.size(tile_idx) == 0 else tile_idx
    tnum = tile_counts(tidx, axis_order)
    space = stitch_pixel_sizes(
        xy_pixel_size,
        dz,
        skew_angle,
        resample_type=resample_type,
        resample_factor=resample_factor,
        objective_scan=objective_scan,
        io_scan=io_scan,
        dsr=dsr,
    )
    rs = (
        np.ones(3)
        if resample_factor is None or np.size(resample_factor) == 0
        else np.asarray(resample_factor, dtype=np.float64).ravel()
    )
    if input_bbox is not None and np.size(input_bbox) == 6:
        b = [int(v) for v in input_bbox]
        im = np.array([b[3] - b[0], b[4] - b[1], b[5] - b[2]], dtype=np.float64)
    if tile_out_bbox is not None and np.size(tile_out_bbox) == 6:
        b = [int(v) for v in tile_out_bbox]
        im = np.array([b[3] - b[0], b[4] - b[1], b[5] - b[2]], dtype=np.float64)
    if not space.dsr and not np.all(rs == 1):
        # MATLAB's [r1 r2 r3] = [x y z] factors, padded from [xy z] or [r], divide [ny nx nz]
        m = rs[::-1]  # [x y z]
        m = np.concatenate((np.ones(4 - m.size) * m[0], m[1:]))
        im = _mround(im / m[[2, 0, 1]])
    sizes = np.repeat(im.reshape(1, 3), nf, axis=0)
    xyz = process_coordinates(
        xyz,
        sizes,
        flipped_tile,
        axis_order,
        data_order,
        space,
        io_scan=io_scan,
        objective_scan=objective_scan,
        ds=ds,
    )
    p = space.pixel_sizes
    if nf > 1 and xcorr_shift:
        pad = np.asarray(stitch_pad_size, dtype=np.float64).ravel()
        pad_size = p * pad[list(dop)] * tnum
    else:
        pad_size = np.zeros(3)
    origin = _col_min(xyz) - pad_size
    dxyz = np.fmax.reduce(xyz + sizes[:, list(dop)] * p, axis=0) + pad_size - origin
    vxyz = _mround(dxyz / p + 2)
    nv = vxyz[list(rev)]
    mip = stitch_mip is not None and bool(np.any(stitch_mip))
    if mip or im[dop[2]] == 1:
        nv[dop[2]] = 1
    has_out = out_bbox is not None and np.size(out_bbox) == 6
    if has_out:
        b = list(out_bbox)
        end = np.array([nv[a] if b[3 + a] is None else float(b[3 + a]) for a in range(3)])
        nv = end - np.array(b[:3], dtype=np.float64)
    rng = np.vstack((nv, nv))
    if nf > 1 and xcorr_shift and not has_out:
        xy = np.asarray(xy_max_offset, dtype=np.float64).ravel()
        zo = np.asarray(z_max_offset, dtype=np.float64).ravel()
        max_offset = np.array([xy[0], xy[-1], zo[0]])
        margin = 2 * max_offset[list(rev)]
        rng = np.vstack((np.fmax(nv - margin, 1), nv + margin))
        if mip or im[dop[2]] == 1:
            rng[:, dop[2]] = 1
    nb = _batch_number(nv, im, batch_size, block_size)
    nbr = np.array(
        [
            _batch_number(rng[0], im, batch_size, block_size),
            _batch_number(rng[1], im, batch_size, block_size),
        ],
        dtype=np.float64,
    )
    return nv, nb, rng, nbr


# ---------------------------------------------------------------------------------------------
# major tiles
# ---------------------------------------------------------------------------------------------


def _cuboid_diff(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``cuboid_diff``: A - B as at most 26 cuboids split on the bounds of A and of their
    overlap (inclusive MATLAB boxes, ``[ymin xmin zmin ymax xmax zmax]``)"""
    if np.all(a[:3] >= b[:3]) and np.all(a[3:] <= b[3:]):
        return np.zeros((0, 6))
    ov = np.sort(np.column_stack((a.reshape(2, 3).T, b.reshape(2, 3).T)), axis=1, kind="stable")
    ov = ov[:, 1:3].T.reshape(-1)
    d = np.vstack((ov.reshape(2, 3), a.reshape(2, 3)))
    cuts = []
    for col in range(3):
        u = np.unique(d[:, col])
        if u.size == 1:
            u = np.array([u[0], u[0]])
        cuts.append(u)
    ys, xs, zs = cuts
    rows = []
    # MATLAB's ndgrid(1 : ny-1, 1 : nx-1, 1 : nz-1)(:): y fastest, then x, then z
    for k in range(zs.size - 1):
        for j in range(xs.size - 1):
            for i in range(ys.size - 1):
                rows.append((ys[i], xs[j], zs[k], ys[i + 1], xs[j + 1], zs[k + 1]))
    cub = np.asarray(rows, dtype=np.float64).reshape(-1, 6)
    return cub[~np.all(cub == ov, axis=1)]


def check_major_tile_cover(coords, major_inds) -> tuple[np.ndarray, np.ndarray]:
    """``[major_cover_mat, uncover_bbox_mat] = check_major_tile_cover(coords_mat, major_inds)``:
    for each minor tile, whether the major tiles' boxes cover its box (``coords``, ``(n, 6)``
    Python boxes), and if not a box of the part left uncovered: MATLAB's ``[min of the starts,
    min of the ends]`` of the pieces, computed on MATLAB's inclusive boxes. Rows of covered and
    major tiles are zeros."""
    from petakit5d.core.bbox import boxes_from_matlab, boxes_to_matlab

    m = boxes_to_matlab(np.asarray(coords, dtype=np.int64).reshape(-1, 6)).astype(np.float64)
    major = np.asarray(major_inds, dtype=bool).ravel()
    n = major.size
    cover = np.ones(n, dtype=bool)
    unc = np.zeros((n, 6))
    for i in range(n):
        if major[i]:
            continue
        ci = m[i].reshape(1, 6)
        cov = np.zeros(1, dtype=bool)
        for j in range(n):
            if not major[j]:
                continue
            if np.all(cov):
                break
            cj = m[j]
            is_cover = np.all((cj[:3] <= ci[:, :3]) & (ci[:, 3:] <= cj[3:]), axis=1)
            cov = cov | is_cover
            nol = (ci[:, 3:] < cj[:3]) | (cj[3:] < ci[:, :3])
            is_ov = ~np.any(nol, axis=1)
            diffs = [
                _cuboid_diff(ci[k], cj) for k in range(is_ov.size) if not is_cover[k] and is_ov[k]
            ]
            dm = np.vstack(diffs) if diffs else np.zeros((0, 6))
            if dm.shape[0]:
                keep = ~cov & ~is_ov
                ci = np.vstack((ci[keep], dm))
                cov = np.concatenate((cov[keep], np.zeros(dm.shape[0], dtype=bool)))
        cover[i] = bool(np.all(cov))
        if not cover[i]:
            unc[i] = np.concatenate((_col_min(ci[:, :3]), _col_min(ci[:, 3:])))
    out = np.zeros((n, 6), dtype=np.int64)
    if (~cover).any():
        out[~cover] = boxes_from_matlab(unc[~cover])
    return cover, out
