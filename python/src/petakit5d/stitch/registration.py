"""The stitching registration of overlapping tile pairs by normalized cross-correlation.

=================================================  =============================================
Python                                             MATLAB
=================================================  =============================================
``cross_correlation_registration``                 cross_correlation_registration_wrapper: the
                                                   pairs of a leading tile (the task of
                                                   stitch_shift_assignment), their result file
``cross_correlation_registration_3d``              cross_correlation_registration_3d
``cross_correlation_registration_2d``              cross_correlation_registration_2d
``cross_correlation_registration_3d_mip_slabs``    cross_correlation_registration_3d_mip_slabs
``compute_tile_mip_slabs``                         compute_tile_mip_slabs
=================================================  =============================================

A pair: the overlap region of tile 2 is the template (cropped to its brightest part along axes
longer than ``dim_num_thrsh``, max-pooled by ``down_sample``, then cut to its rows, columns and
slices with data), the overlap of tile 1 widened by the max offsets the image; the template's
offset at the first maximum of their normalized cross-correlation in the window of the allowed
shifts (``stitch.xcorr``: the 3D NCC of ``kernels.xcorr``, ``normxcorr2`` in 2D) gives the shift of
tile 2 against its stage position. The large-file variant registers the tiles' three max-pooled
MIP slabs (each one as a 3D pair, in tasks) and averages the shifts along each axis from the two
slabs that see it, weighted by their squared peaks.

Values: positions and shifts are physical, x, y, z by name, as MATLAB's: cuboids
``[x0 y0 z0 x1 y1 z1]`` (um), ``xyz_voxelsizes`` (x, y, z), relative shifts (x, y, z, voxels;
negative toward, positive away). Per-axis vectors are ``(z, y, x)``: ``down_sample`` (MATLAB's
``downSample``, in the tiles' array order), ``max_offset`` (MATLAB's ``MaxOffset``, ``[x y z]``)
and ``pool_size`` (MATLAB's ``poolSize``, ``[y x z]`` in groups of three). The tiles' data order
is the array axes of x, y and z (``stitch.geometry.parse_data_order``, the Python form of MATLAB's
``data_order_mat``); the pair functions also take the data order's text (``'z,y,x'``). Pairs are
0-based.

Files keep MATLAB's form: the wrapper's ``relative_shift_mat`` (pairs x 3), ``max_xcorr_mat``
(pairs x 1), 1-based ``pair_indices`` (and with ``large_file`` the cells of each pair's slab
results); the pair functions' ``relative_shift`` and ``max_xcorr`` (the slabs' also
``relative_shift_mat`` and ``max_xcorr_mat``), each written under a temporary name and renamed.
Tiles are tiff files or zarr arrays (the port's native ones and MATLAB's).

The order of the axes follows MATLAB's code: the 3D pair takes x, y and z as
``data_order_mat(1..3)`` (6a8db79); the 2D pair and the slabs keep MATLAB's own axis choices
(the same for the data orders with z last).
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Sequence

import numpy as np

from petakit5d.params import Param, register
from petakit5d.stitch.geometry import data_order_axes_from_matlab, data_order_axes_to_matlab
from petakit5d.stitch.xcorr import RegistrationError, _bad_index

__all__ = [
    "RegistrationError",
    "compute_tile_mip_slabs",
    "cross_correlation_registration",
    "cross_correlation_registration_2d",
    "cross_correlation_registration_3d",
    "cross_correlation_registration_3d_mip_slabs",
]

# the default max offsets: MATLAB's [300, 300, 50] (x, y, z)
MAX_OFFSET = (50, 300, 300)
# compute_tile_mip_slabs: the MIP task's batches and the major downsampling of the slabs
_MIP_BATCH = (1024, 1024, 1024)
_MAJOR_DOWNSAMPLE = 8


# ---------------------------------------------------------------------------------------------
# helpers: axes, files
# ---------------------------------------------------------------------------------------------


def _axes(data_order) -> tuple[int, ...]:
    """the array axes of x, y and z of a data order string or of the axes themselves"""
    from petakit5d.stitch.geometry import parse_data_order

    if isinstance(data_order, str):
        return parse_data_order(data_order)
    return tuple(int(v) for v in np.asarray(data_order).ravel())


def _matlab_order(data_order) -> tuple[list[int], list[int]]:
    """MATLAB's ``data_order_mat`` (the dimension, 1 to 3, of x, y and z) and
    ``data_order_reverse_mat`` (for each dimension, which of x, y, z it holds)"""
    dom = [int(d) for d in data_order_axes_to_matlab(_axes(data_order))]
    rev = [0, 0, 0]
    for k, d in enumerate(dom):
        rev[d - 1] = k + 1
    return dom, rev


def _vec(v, n: int = 3) -> np.ndarray:
    a = np.asarray(v, dtype=np.float64).ravel()
    if a.size < n:
        raise _bad_index(a.size)
    return a


def _down_sample_m(down_sample) -> np.ndarray:
    """MATLAB's downSample ([d1 d2 d3], the data's order) of a (z, y, x) down_sample"""
    from petakit5d.core.axes import to_yxz

    return np.asarray(to_yxz(tuple(float(v) for v in _vec(down_sample))), dtype=np.float64)


def _max_offset_xyz(max_offset) -> np.ndarray:
    from petakit5d.core.axes import to_xyz

    return np.asarray(to_xyz(tuple(float(v) for v in _vec(max_offset))), dtype=np.float64)


def _ext(path: str) -> str:
    from petakit5d.core.paths import fileparts

    return fileparts(path)[2]


def _image_size(path: str, ext: str) -> list[int]:
    """blockedImage's Size of a tile: tiff [ny nx nz]; zarr the array's [ny nx nz] (2 values for
    a 2D array)"""
    from petakit5d.io import get_image_size

    if ext not in (".tif", ".zarr"):
        raise RegistrationError(
            "Unrecognized function or variable 'bim_1'.", "MATLAB:UndefinedFunction"
        )
    s = [int(v) for v in get_image_size(path)]
    if ext == ".tif":
        return [s[1], s[2], s[0]]
    return [s[1], s[2], s[0]] if len(s) == 3 else s


def _read_region(
    path: str,
    ext: str,
    start: Sequence[float],
    stop: Sequence[float],
    num_threads: int | None = None,
) -> np.ndarray:
    """getIORegion(start, stop) of the tile's adapter: MATLAB's 1-based, inclusive [d1 d2 d3]
    bounds, as a (z, y, x) volume"""
    from petakit5d.io import read_tiff, read_zarr

    s = [int(v) for v in start]
    e = [int(v) for v in stop]
    if ext == ".tif":
        im = read_tiff(path, range=(s[2] - 1, e[2]), num_threads=num_threads)
        if im.ndim == 2:
            im = im[np.newaxis]
        return np.ascontiguousarray(im[:, s[0] - 1 : e[0], s[1] - 1 : e[1]])
    a = read_zarr(
        path, input_bbox=(s[2] - 1, s[0] - 1, s[1] - 1, e[2], e[0], e[1]), num_threads=num_threads
    )
    return a[np.newaxis] if a.ndim == 2 else a


def _start_stop(cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes):
    """the overlap's first and last voxel (1-based, x, y, z) in tile 1 and tile 2"""
    from petakit5d.core.dtypes import mround

    vs = _vec(xyz_voxelsizes)[:3]
    c1 = _vec(cuboid_1, 6)
    c2 = _vec(cuboid_2, 6)
    co = _vec(cuboid_overlap_12, 6)
    s1 = mround((co[0:3] - c1[0:3]) / vs) + 1
    s2 = mround((co[0:3] - c2[0:3]) / vs) + 1
    t1 = mround((co[3:6] - c1[0:3]) / vs) + 1
    t2 = mround((co[3:6] - c2[0:3]) / vs) + 1
    return s1, s2, t1, t2


def _pool(region: np.ndarray, ds_m: np.ndarray, num_threads: int | None = None) -> np.ndarray:
    """max_pooling_3d(region, downSample) (MATLAB's [d1 d2 d3] pool sizes)"""
    from petakit5d.kernels.mip import max_pooling_3d_matlab

    out = max_pooling_3d_matlab(region, ds_m.reshape(1, 3), num_threads=num_threads)
    return out[np.newaxis] if out.ndim == 2 else out


def _save(path: str, variables: dict) -> None:
    """save('-v7.3', <path>_<uuid>.mat, ...), fileattrib +w g, movefile onto path"""
    from petakit5d.core.ids import get_uuid
    from petakit5d.io.matfile import write_mat
    from petakit5d.io.tiff import group_writable

    tmp = f"{path[:-4]}_{get_uuid()}.mat"
    write_mat(tmp, variables)
    group_writable(tmp)
    os.replace(tmp, path)


def _no_shift(xcorr_fullpath: str, done: bool) -> tuple[np.ndarray, float]:
    print(
        "The non-empty area of the region 2 is too small, skip the xcorr computing, set "
        "relative_shift as [0, 0, 0]\n. Save results ..."
    )
    rel = np.zeros(3)
    if xcorr_fullpath:
        _save(xcorr_fullpath, {"relative_shift": rel.reshape(1, 3), "max_xcorr": 1.0})
        if done:
            print("Done!")
    return rel, 1.0


def _bounds_zyx(lb_d: np.ndarray, ub_d: np.ndarray) -> np.ndarray:
    """MATLAB's maxShifts [lb; ub] in its dimensions [d1 d2 d3] as (2, 3) rows (z, y, x)"""
    return np.array([[lb_d[2], lb_d[0], lb_d[1]], [ub_d[2], ub_d[0], ub_d[1]]])


# ---------------------------------------------------------------------------------------------
# cross_correlation_registration_3d
# ---------------------------------------------------------------------------------------------


def cross_correlation_registration_3d(
    img_fullpath_1: str,
    img_fullpath_2: str,
    xcorr_fullpath: str,
    cuboid_1,
    cuboid_2,
    cuboid_overlap_12,
    xyz_voxelsizes,
    data_order,
    *,
    down_sample: Sequence[float] = (1, 1, 1),
    max_offset: Sequence[float] = MAX_OFFSET,
    dim_num_thrsh: float = 10000,
    blank_num_trsh: float | None = None,
    num_threads: int | None = None,
    backend: str | None = None,
) -> tuple[np.ndarray, np.floating | float] | None:
    """``[relative_shift, max_xcorr] = cross_correlation_registration_3d(...)``: the shift of
    tile 2 (``img_fullpath_2``) against its stage position, from the 3D NCC of the two tiles'
    overlap (``cuboid_overlap_12``; ``cuboid_1``, ``cuboid_2`` the tiles' cuboids, um).

    Returns ``relative_shift`` (x, y, z voxels, float64) and the NCC's maximum (single), or
    ``[0, 0, 0]`` and 1 when the template has too little data. With ``xcorr_fullpath`` they are
    also saved there (MATLAB's names); an existing file is kept and None returned.
    ``num_threads``: threads of the reads, the pooling and the NCC (None: the CPUs available);
    ``backend``: the CPU FFT backend of the NCC (None: ``petakit5d.fft``'s default).
    """
    from petakit5d.core.bbox import box_to_matlab
    from petakit5d.stitch.xcorr import (
        cap_region_inds,
        crop_subregion_by_intensity,
        normxcorr3_max_shift,
        template_indices,
    )

    ds_m = _down_sample_m(down_sample)
    mo = _max_offset_xyz(max_offset)
    dom, rev = _matlab_order(data_order)
    dom0 = [d - 1 for d in dom]
    rev0 = [k - 1 for k in rev]
    xdo, ydo, zdo = dom0
    if xcorr_fullpath and os.path.exists(xcorr_fullpath):
        print(f"The xcorr result {xcorr_fullpath} already exists!")
        return None
    print(f"Compute xcorr shifts for tiles:\n    {img_fullpath_1}\n    {img_fullpath_2}")
    ext = _ext(img_fullpath_1)
    sz_1 = _image_size(img_fullpath_1, ext)
    s1, s2, t1, t2 = _start_stop(cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes)

    region_2 = _read_region(img_fullpath_2, ext, s2[rev0], t2[rev0], num_threads)
    # crop region 2 to its part with the most signal along axes over dim_num_thrsh
    region_2, box = crop_subregion_by_intensity(
        region_2, dim_num_thrsh, blank_num_trsh, num_threads=num_threads
    )
    if region_2.shape[0] == 1:
        # MATLAB's size() of one slice has 2 values, crop_bbox 5
        raise _bad_index(5)
    crop = np.asarray(box_to_matlab(box), dtype=np.float64)
    s1 = s1 + (crop[dom0] - 1)
    t1 = s1 + (crop[[d + 3 for d in dom0]] - crop[dom0])
    s2 = s2 + (crop[dom0] - 1)
    pooled = bool(np.any(ds_m != 1))
    if pooled:
        region_2 = _pool(region_2, ds_m, num_threads)

    # the image: the overlap of tile 1 and the max offsets around it
    maxoff = np.array(
        [
            math.ceil(mo[0] / ds_m[xdo]) * ds_m[xdo],
            math.ceil(mo[1] / ds_m[ydo]) * ds_m[ydo],
            math.ceil(mo[2] / ds_m[zdo]) * ds_m[zdo],
        ]
    )
    sr1 = np.maximum(s1 - maxoff, 1)
    if len(sz_1) < 3:
        raise _bad_index(len(sz_1))
    tr1 = np.minimum(t1 + maxoff, np.asarray([sz_1[d] for d in dom0], dtype=np.float64))
    region_1 = _read_region(img_fullpath_1, ext, sr1[rev0], tr1[rev0], num_threads)
    if pooled:
        region_1 = _pool(region_1, ds_m, num_threads)

    # the template: the rows, columns and slices with data in a quarter (else an eighth) of
    # their area, less the margins cap_region_2_inds cuts
    sz_2 = [region_2.shape[1], region_2.shape[2], region_2.shape[0]]
    nonzero = region_2 != 0
    y, x, z = template_indices(nonzero, 4)
    if not (y.any() and x.any() and z.any()):
        y, x, z = template_indices(nonzero, 8)
    maxoff_d = maxoff[rev0] / ds_m
    y_inds = cap_region_inds(y, sz_2[0], maxoff_d[0])
    x_inds = cap_region_inds(x, sz_2[1], maxoff_d[1])
    z_inds = cap_region_inds(z, sz_2[2], maxoff_d[2])
    if not (z_inds.size and x_inds.size and y_inds.size):
        return _no_shift(xcorr_fullpath, True)
    if (
        y_inds[0] != 1
        or y_inds[-1] != sz_2[0]
        or x_inds[0] != 1
        or x_inds[-1] != sz_2[1]
        or z_inds[0] != 1
        or z_inds[-1] != sz_2[2]
    ):
        region_2 = region_2[
            z_inds[0] - 1 : z_inds[-1], y_inds[0] - 1 : y_inds[-1], x_inds[0] - 1 : x_inds[-1]
        ]
    src2_d = np.array([y_inds[0], x_inds[0], z_inds[0]], dtype=np.float64)
    src2 = src2_d[dom0]

    # the allowed shifts, the NCC and the shift
    ds_xyz = ds_m[dom0]
    sp2_down = (sr1 - s1) / ds_xyz - (src2 - 1)
    maxoff_ds = maxoff / ds_xyz
    lb = (-sp2_down - maxoff_ds)[rev0]
    ub = (-sp2_down + maxoff_ds)[rev0]
    off, max_xcorr, _ = normxcorr3_max_shift(
        np.ascontiguousarray(region_2, dtype=np.float32),
        np.ascontiguousarray(region_1, dtype=np.float32),
        _bounds_zyx(lb, ub),
        threads=num_threads,
        backend=backend,
    )
    off_d = np.array([off[1], off[2], off[0]], dtype=np.float64)
    sp2 = (sr1 - s1) - (src2 - 1) * ds_xyz
    rel = off_d[dom0] * ds_xyz + sp2
    rel = (rel * ((s1 >= s2) - 0.5)) * 2
    if xcorr_fullpath:
        print("Save results ...", end="")
        _save(xcorr_fullpath, {"relative_shift": rel.reshape(1, 3), "max_xcorr": max_xcorr})
        print("Done!")
    return rel, max_xcorr


# ---------------------------------------------------------------------------------------------
# cross_correlation_registration_2d
# ---------------------------------------------------------------------------------------------


def _matlab_form(region: np.ndarray, dom: Sequence[int] | None = None) -> np.ndarray:
    """a (z, y, x) region as MATLAB's array (permuted by dom), as the port's image for a 2D one"""
    from petakit5d.core.arrays import to_matlab

    m = to_matlab(region)
    if dom is not None:
        m = np.transpose(m, [d - 1 for d in dom])
    if m.ndim == 3 and m.shape[2] == 1:
        m = m[:, :, 0]
    return np.ascontiguousarray(m)


def cross_correlation_registration_2d(
    img_fullpath_1: str,
    img_fullpath_2: str,
    xcorr_fullpath: str,
    cuboid_1,
    cuboid_2,
    cuboid_overlap_12,
    xyz_voxelsizes,
    data_order,
    *,
    down_sample: Sequence[float] = (1, 1, 1),
    max_offset: Sequence[float] = MAX_OFFSET,
    dim_num_thrsh: float = 10000,
    num_threads: int | None = None,
) -> tuple[np.ndarray, float] | None:
    """``[relative_shift, max_xcorr] = cross_correlation_registration_2d(...)``: as
    ``cross_correlation_registration_3d`` for 2D tiles, with ``normxcorr2`` (``stitch.xcorr``).

    The template is not cropped by intensity (``dim_num_thrsh`` is not used, as in MATLAB). The
    max offsets are taken as MATLAB takes them: x's along the data's first dimension, y's along
    its second (the same for 'z,y,x' tiles); tiles whose z is not last are permuted to x, y, z.
    """
    del dim_num_thrsh
    from petakit5d.stitch.xcorr import normxcorr2_max_shift, template_indices

    ds_m = _down_sample_m(down_sample)
    mo = _max_offset_xyz(max_offset)
    dom, rev = _matlab_order(data_order)
    dom0 = [d - 1 for d in dom]
    rev0 = [k - 1 for k in rev]
    xdo = np.array([d == 1 for d in dom])
    ydo = np.array([d == 2 for d in dom])
    max_x = mo[xdo]
    max_y = mo[ydo]
    if xcorr_fullpath and os.path.exists(xcorr_fullpath):
        print(f"The xcorr result {xcorr_fullpath} already exists!")
        return None
    print(f"Compute xcorr shifts for tiles:\n    {img_fullpath_1}\n    {img_fullpath_2}")
    ext = _ext(img_fullpath_1)
    sz_1 = _image_size(img_fullpath_1, ext)
    s1, s2, t1, t2 = _start_stop(cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes)

    region_2 = _read_region(img_fullpath_2, ext, s2[rev0], t2[rev0], num_threads)
    region_2 = region_2.astype(np.float32)
    pooled = bool(np.any(ds_m != 1))
    if pooled:
        region_2 = _pool(region_2, ds_m, num_threads)
    maxoff_x = float((np.ceil(max_x / ds_m[xdo]) * ds_m[xdo])[0])
    maxoff_y = float((np.ceil(max_y / ds_m[ydo]) * ds_m[ydo])[0])
    border = np.array([maxoff_x, maxoff_y, 0.0])
    sr1 = np.maximum(s1 - border, 1)
    if len(sz_1) < 3:
        raise _bad_index(len(sz_1))
    tr1 = np.minimum(t1 + border, np.asarray([sz_1[d] for d in dom0], dtype=np.float64))
    region_1 = _read_region(img_fullpath_1, ext, sr1[rev0], tr1[rev0], num_threads)
    region_1 = region_1.astype(np.float32)
    if pooled:
        region_1 = _pool(region_1, ds_m, num_threads)

    # the template: the columns, then rows, of its first slice with data in a quarter (else an
    # eighth) of their length
    nonzero = region_2 != 0
    y, x, _ = template_indices(nonzero, 4, two_d=True)
    if not (x.any() and y.any()):
        y, x, _ = template_indices(nonzero, 8, two_d=True)
    if not (x.any() and y.any()):
        rel, mx = _no_shift(xcorr_fullpath, False)
        print("Done!")
        return rel, mx
    # MATLAB's ~all(y_inds([1, end]) & x_inds([1, end])): a column and a row, so a 2 x 2
    # matrix and an if on two values: the template is cut to its rows and columns with data
    # when its first or last row has none, or both its first and last columns
    if not (y[0] and y[-1]) or not (x[0] or x[-1]):
        region_2 = np.ascontiguousarray(region_2[:1][:, y][:, :, x])
    src2_d = np.array([np.flatnonzero(y)[0] + 1, np.flatnonzero(x)[0] + 1, 1], dtype=np.float64)
    src2 = src2_d[dom0]

    ds_xyz = ds_m[dom0]
    sp2_down = (sr1 - s1) / ds_xyz - (src2 - 1)
    maxoff_ds = np.array([maxoff_x, maxoff_y, 1.0]) / ds_xyz
    lb = (-sp2_down - maxoff_ds)[rev0]
    ub = (-sp2_down + maxoff_ds)[rev0]
    permute = dom[2] != 3
    if permute:
        lb, ub = lb[dom0], ub[dom0]
    T = _matlab_form(region_2, dom if permute else None)
    A = _matlab_form(region_1, dom if permute else None)
    off, max_xcorr, _ = normxcorr2_max_shift(T, A, _bounds_zyx(lb, ub), threads=num_threads)
    off_d = np.array([off[1], off[2], off[0]], dtype=np.float64)
    if permute:
        off_d = off_d[rev0]
    sp2 = (sr1 - s1) - (src2 - 1) * ds_xyz
    rel = off_d[dom0] * ds_xyz + sp2
    rel = (rel * ((s1 >= s2) - 0.5)) * 2
    max_xcorr = float(max_xcorr)
    if xcorr_fullpath:
        print("Save results ...")
        _save(xcorr_fullpath, {"relative_shift": rel.reshape(1, 3), "max_xcorr": max_xcorr})
    print("Done!")
    return rel, max_xcorr


# ---------------------------------------------------------------------------------------------
# cross_correlation_registration_3d_mip_slabs
# ---------------------------------------------------------------------------------------------


def _pool_matlab(pool_size) -> list[float]:
    """MATLAB's poolSize ([y x z] groups) of a (z, y, x) pool size"""
    from petakit5d.core.axes import to_yxz

    if pool_size is None:
        return []
    p = [float(v) for v in np.asarray(pool_size, dtype=np.float64).ravel()]
    return [float(v) for g in range(0, len(p) - len(p) % 3, 3) for v in to_yxz(p[g : g + 3])]


def _matlab_size(path: str) -> list[int]:
    """getImageSize: [ny nx nz]"""
    from petakit5d.io import get_image_size

    s = [int(v) for v in get_image_size(path)]
    return [s[1], s[2], s[0]] if len(s) >= 3 else s


def cross_correlation_registration_3d_mip_slabs(
    img_fullpath_1: str,
    img_fullpath_2: str,
    xcorr_fullpath: str,
    cuboid_1,
    cuboid_2,
    cuboid_overlap_12,
    xyz_voxelsizes,
    data_order,
    *,
    down_sample: Sequence[float] = (1, 1, 1),
    max_offset: Sequence[float] = MAX_OFFSET,
    mip_dir_str: str = "",
    pool_size: Sequence[float] | None = None,
    major_downsample: float = _MAJOR_DOWNSAMPLE,
    dim_num_thrsh: float = 10000,
    blank_num_trsh: float = 200,
    parse_cluster: bool = True,
    mcc_mode: bool = False,
    config_file: str = "",
    num_threads: int | None = None,
):
    """``[relative_shift, max_xcorr, relative_shift_mat, max_xcorr_mat] =
    cross_correlation_registration_3d_mip_slabs(...)``: the registration of two large tiles from
    their MIP slabs (``<tile folder>/<mip_dir_str>/<name>_<axis>.zarr``, pooled by ``pool_size``).

    Each slab pair (in MATLAB's axis order of the data order: x, y, z for 'z,y,x' tiles) is
    registered by a ``cross_correlation_registration`` task (``<xcorr>_<axis>.mat``; three tries
    with growing memory), its axis downsampled by ``major_downsample``, the others by
    ``down_sample / pool_size[3:6]`` (rounded). A slab whose peak is 1 at no shift while its peak is
    0.2 over the median of the three is given weight 0; each axis's shift is the average of its two
    slabs' shifts (in voxels of the tile) weighted by their squared peaks, rounded; the peak is
    ``sum(p^3) / sum(p^2)``. Saved to ``xcorr_fullpath`` (an existing file is read instead).
    ``dim_num_thrsh`` and ``blank_num_trsh`` are not used (as in MATLAB); ``num_threads`` is the
    slab tasks'.
    """
    del dim_num_thrsh, blank_num_trsh
    from petakit5d.core.dtypes import mround
    from petakit5d.core.paths import fileparts
    from petakit5d.core.text import through_mat2str, through_text
    from petakit5d.io.matfile import read_mat
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.matlabcompat import elementary, stats
    from petakit5d.stitch.geometry import check_cuboids_overlaps

    names = ["relative_shift", "max_xcorr", "relative_shift_mat", "max_xcorr_mat"]
    if xcorr_fullpath and os.path.exists(xcorr_fullpath):
        print(f"The xcorr result {xcorr_fullpath} already exists!")
        a = read_mat(xcorr_fullpath, names)
        return tuple(a.get(k) for k in names)
    ps = _pool_matlab(pool_size)
    pool_1 = np.ones(3)
    if len(ps) >= 6:
        pool_1 = np.asarray(ps[3:6])
        ps = ps[:3]
    if len(ps) < 3:
        raise _bad_index(len(ps))
    pool = np.asarray(ps, dtype=np.float64)
    ds_m = _down_sample_m(down_sample)
    mo = _max_offset_xyz(max_offset)
    vs = _vec(xyz_voxelsizes)[:3]
    c1 = _vec(cuboid_1, 6)[:6]
    c2 = _vec(cuboid_2, 6)[:6]
    co = _vec(cuboid_overlap_12, 6)[:6]
    dom, _ = _matlab_order(data_order)
    dom0 = [d - 1 for d in dom]
    print(
        "Compute xcorr shifts for tiles with MIP slab method:\n"
        f"    {img_fullpath_1}\n    {img_fullpath_2}"
    )
    axis_strs = [("x", "y", "z")[d] for d in dom0]
    inputs, outputs, tasks = [], [], []
    overlap_sizes = np.zeros((3, 3))
    for i in range(3):
        slabs = []
        for path in (img_fullpath_1, img_fullpath_2):
            data_path, fsn, _ = fileparts(path)
            slabs.append(f"{data_path}/{mip_dir_str}/{fsn}_{axis_strs[i]}.zarr")
        vs_i = vs * pool_1[dom0]
        vs_i[dom0[i]] = vs[dom0[i]] * pool[i]
        mo_i = np.ceil(mo / pool_1)
        mo_i[i] = math.ceil(mo[i] / pool[i])
        ds_i = mround(ds_m / pool_1)
        ds_i[i] = major_downsample
        sz_i1 = _matlab_size(slabs[0])
        sz_i2 = _matlab_size(slabs[1])
        c1i = c1.copy()
        c1i[3:6] = c1i[0:3] + (np.asarray([sz_i1[d] for d in dom0]) - 2) * vs_i
        c2i = c2.copy()
        c2i[3:6] = c2i[0:3] + (np.asarray([sz_i2[d] for d in dom0]) - 2) * vs_i
        is_ov, ov = check_cuboids_overlaps(c1i, c2i, False)
        if bool(np.asarray(is_ov).ravel()[0]):
            coi = np.asarray(ov, dtype=np.float64).reshape(-1)[:6]
        else:
            c1i, c2i, coi = c1, c2, co
        overlap_sizes[i] = np.ceil((coi[3:6] - coi[0:3]) / vs_i) + mo_i
        out_i = f"{xcorr_fullpath[:-4]}_{axis_strs[i]}.mat"
        inputs.append(slabs[0])
        outputs.append(out_i)
        # the values as MATLAB's function string gives them to the task (mat2str: 15 digits;
        # downSample with '%.20d')
        ds_t = through_text(ds_i, "%.20d")
        mo_t = through_mat2str(mo_i)
        tasks.append(
            Task.of(
                cross_correlation_registration,
                slabs[0],
                [slabs[1]],
                out_i,
                [[i, i + 1]],
                [float(v) for v in through_mat2str(c1i)],
                [[float(v) for v in through_mat2str(c2i)]],
                [[float(v) for v in through_mat2str(coi)]],
                [float(v) for v in through_mat2str(vs_i)],
                [int(v) for v in _axes(data_order)],
                stitch_2d=False,
                down_sample=[float(ds_t[2]), float(ds_t[0]), float(ds_t[1])],
                max_offset=[float(mo_t[2]), float(mo_t[1]), float(mo_t[0])],
                large_file=False,
                mip_dir_str="",
                pool_size=None,
                **({} if num_threads is None else {"num_threads": int(num_threads)}),
            )
        )
    mem = float(np.max((overlap_sizes[:, 0] * overlap_sizes[:, 1]) * overlap_sizes[:, 2])) * 4
    mem = mem / 2**30 * 50
    done: list[bool] = []
    for k in range(3):
        done = framework.run_tasks(
            inputs,
            outputs,
            tasks,
            cpus_per_task=2,
            max_trial_num=2,
            parse_cluster=parse_cluster,
            mem_allocate=mem * 2**k,
            min_task_job_num=3,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        if all(done):
            break
    if not all(done):
        raise RegistrationError("Some pairs of xcorr results are missing!")

    # the slabs' results; a failed one (peak 1 at no shift, 0.2 over the median) weighs 0
    rsm = np.zeros((3, 3))
    mxm = np.zeros((3, 1))
    for i in range(3):
        a = read_mat(outputs[i], ["relative_shift_mat", "max_xcorr_mat"])
        rsm[i] = np.asarray(a["relative_shift_mat"], dtype=np.float64).ravel()[:3]
        mxm[i, 0] = float(np.asarray(a["max_xcorr_mat"], dtype=np.float64).ravel()[0])
    failed = (mxm[:, 0] == 1) & np.all(rsm == 0, axis=1)
    if failed.any():
        if float(mxm.max()) - float(stats.median(mxm[:, 0])) > 0.2:
            mxm[failed, 0] = 0
    # each axis from the two slabs that see it, weighted by their squared peaks
    rel = np.zeros(3)
    for i in range(3):
        diff = [k for k in range(3) if k != dom0[i]]
        w = elementary.power(mxm[diff, 0], 2)
        rsm[diff, i] = rsm[diff, i] * pool_1[dom0[i]]
        prod = rsm[diff, i] * w
        rel[i] = (prod[0] + prod[1]) / (w[0] + w[1])
    rel = mround(rel)
    p3 = elementary.power(mxm[:, 0], 3)
    p2 = elementary.power(mxm[:, 0], 2)
    max_xcorr = float(((p3[0] + p3[1]) + p3[2]) / ((p2[0] + p2[1]) + p2[2]))
    if xcorr_fullpath:
        print("Save results ...")
        _save(
            xcorr_fullpath,
            {
                "relative_shift": rel.reshape(1, 3),
                "max_xcorr": max_xcorr,
                "relative_shift_mat": rsm,
                "max_xcorr_mat": mxm,
            },
        )
        print("Done!")
    return rel, max_xcorr, rsm, mxm


# ---------------------------------------------------------------------------------------------
# cross_correlation_registration_wrapper
# ---------------------------------------------------------------------------------------------


@register(
    "cross_correlation_registration_wrapper",
    Param("imgFullpath_1", "@ischar", "The leading tile (tiff or zarr)."),
    Param("imgFullpath_2", "@iscell", "The tiles paired with it."),
    Param("xcorrFullpath", "@ischar", "Result file (.mat)."),
    Param(
        "pair_indices", "@isnumeric", "The pairs (k x 2, 0-based tile numbers).", convert="index1"
    ),
    Param(
        "cuboid_1",
        "@isnumeric",
        "Cuboid [x0 y0 z0 x1 y1 z1] (um) of the leading tile.",
        convert="same",
    ),
    Param("cuboid_2", "@isnumeric", "Cuboids (k x 6) of the paired tiles.", convert="same"),
    Param(
        "cuboid_overlap_12", "@isnumeric", "Overlap cuboids (k x 6) of the pairs.", convert="same"
    ),
    Param("xyz_voxelsizes", "@isnumeric", "Voxel sizes (x, y, z, um).", convert="same"),
    Param(
        "data_order_mat",
        "@isnumeric",
        "Array axes of x, y and z of the tiles (stitch.geometry.parse_data_order).",
        py="data_order_axes",
        convert=(data_order_axes_from_matlab, data_order_axes_to_matlab),
    ),
    Param("Stitch2D", "@islogical", "2D tiles (normxcorr2)."),
    Param("downSample", "@isnumeric", "Max pooling (z, y, x) of the regions.", convert="size_yxz"),
    Param("MaxOffset", "@isnumeric", "Max offsets (z, y, x, voxels).", convert="per_axis_xyz"),
    Param("largeFile", "@islogical", "Register the tiles' MIP slabs."),
    Param("mipDirStr", "@ischar", "Folder of the MIP slabs next to the tiles."),
    Param("poolSize", "@isnumeric", "Pool sizes (z, y, x) of the MIP slabs.", convert="size_yxz"),
    Param("dimNumThrsh", "@isnumeric", "Longest template along an axis (cropped by intensity)."),
    Param("parseCluster", "@islogical", "Run the slab tasks on the cluster."),
    Param("mccMode", "@islogical", "Compiled tasks."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param(
        "numThreads",
        None,
        "Threads of each pair's reads and NCC (None: the CPUs available).",
        py="num_threads",
        python_only=True,
    ),
    task="S04",
)
def cross_correlation_registration(
    img_fullpath_1: str,
    img_fullpath_2: Sequence[str],
    xcorr_fullpath: str,
    pair_indices,
    cuboid_1,
    cuboid_2,
    cuboid_overlap_12,
    xyz_voxelsizes,
    data_order_axes: Sequence[int],
    *,
    stitch_2d: bool = False,
    down_sample: Sequence[float] = (1, 1, 1),
    max_offset: Sequence[float] = MAX_OFFSET,
    large_file: bool = False,
    mip_dir_str: str = "",
    pool_size: Sequence[float] | None = None,
    dim_num_thrsh: float = 10000,
    parse_cluster: bool = True,
    mcc_mode: bool = False,
    config_file: str = "",
    num_threads: int | None = None,
) -> None:
    """The registration of a leading tile with the tiles of its pairs, saved to
    ``xcorr_fullpath`` (``relative_shift_mat``, ``max_xcorr_mat``, 1-based ``pair_indices``; with
    ``large_file`` also each pair's slab results): 2D pairs with ``stitch_2d``, large tiles on
    their MIP slabs (``<xcorr>_mip_slabs.mat`` for one pair, ``<xcorr>_<second tile>_mip_slabs.mat``
    for more), else 3D pairs. An existing result is kept."""
    if os.path.exists(xcorr_fullpath):
        print(f"The xcorr result {xcorr_fullpath} already exists!")
        return
    pairs = np.asarray(pair_indices, dtype=np.int64).reshape(-1, 2)
    n_f = pairs.shape[0]
    c2 = np.asarray(cuboid_2, dtype=np.float64).reshape(n_f, -1)
    co = np.asarray(cuboid_overlap_12, dtype=np.float64).reshape(n_f, -1)
    rsm = np.zeros((n_f, 3))
    mxm = np.zeros((n_f, 1))
    rs_cells = np.empty((n_f, 1), dtype=object)
    mx_cells = np.empty((n_f, 1), dtype=object)
    for i in range(n_f):
        t0 = time.perf_counter()
        args = (img_fullpath_1, img_fullpath_2[i])
        kw = {"down_sample": down_sample, "max_offset": max_offset, "num_threads": num_threads}
        if stitch_2d:
            rel, mx = cross_correlation_registration_2d(
                *args,
                "",
                cuboid_1,
                c2[i],
                co[i],
                xyz_voxelsizes,
                data_order_axes,
                dim_num_thrsh=dim_num_thrsh,
                **kw,
            )
        elif large_file:
            if n_f == 1:
                path_i = f"{xcorr_fullpath[:-4]}_mip_slabs.mat"
            else:
                path_i = f"{xcorr_fullpath[:-4]}_{pairs[i, 1] + 1}_mip_slabs.mat"
            rel, mx, rs_cells[i, 0], mx_cells[i, 0] = cross_correlation_registration_3d_mip_slabs(
                *args,
                path_i,
                cuboid_1,
                c2[i],
                co[i],
                xyz_voxelsizes,
                data_order_axes,
                mip_dir_str=mip_dir_str,
                pool_size=pool_size,
                dim_num_thrsh=dim_num_thrsh,
                parse_cluster=parse_cluster,
                mcc_mode=mcc_mode,
                config_file=config_file,
                **kw,
            )
        else:
            rel, mx = cross_correlation_registration_3d(
                *args,
                "",
                cuboid_1,
                c2[i],
                co[i],
                xyz_voxelsizes,
                data_order_axes,
                dim_num_thrsh=dim_num_thrsh,
                **kw,
            )
        rsm[i] = np.asarray(rel, dtype=np.float64).ravel()
        mxm[i, 0] = float(mx)
        print(f"Elapsed time is {time.perf_counter() - t0:.6f} seconds.")
    print("Save results ... ", end="")
    variables = {
        "relative_shift_mat": rsm,
        "max_xcorr_mat": mxm,
        "pair_indices": (pairs + 1).astype(np.float64),
    }
    if large_file:
        variables["relative_shift_mat_cell"] = rs_cells
        variables["max_xcorr_mat_cell"] = mx_cells
    _save(xcorr_fullpath, variables)
    print("Done!")


# ---------------------------------------------------------------------------------------------
# compute_tile_mip_slabs
# ---------------------------------------------------------------------------------------------


def compute_tile_mip_slabs(
    zarr_fullpaths: Sequence[str],
    pool_size: Sequence[float],
    *,
    parse_cluster: bool = True,
    mcc_mode: bool = False,
    config_file: str = "",
) -> str:
    """``mipDirStr = compute_tile_mip_slabs(zarrFullpaths, poolSize, ...)``: the max-pooled MIP
    slabs of large tiles (``tools.mip.mip_zarr`` tasks with ``mip_slab``, batches of 1024) in
    ``<tile folder>/MIP_Slabs_pooling_<pool sizes, [y x z]>/``, existing slabs kept; three tries
    with growing memory. Returns the folder's name."""
    from petakit5d.core.dtypes import mround
    from petakit5d.core.paths import fileparts
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.stitch.distance import mip_slab_name
    from petakit5d.tools.mip import mip_zarr

    ps = [float(v) for v in np.asarray(pool_size, dtype=np.float64).ravel()]
    mip_dir = mip_slab_name(ps)
    parts = [fileparts(p) for p in zarr_fullpaths]
    for d in sorted({p[0] for p in parts}):
        path = f"{d}/{mip_dir}"
        if not os.path.isdir(path):
            os.mkdir(path)
    tasks = [
        Task.of(
            mip_zarr,
            p,
            result_dir_name=mip_dir,
            axis=(1, 1, 1),
            batch_size=_MIP_BATCH,
            pool_size=tuple(ps),
            mip_slab=True,
            parse_cluster=parse_cluster,
            parse_parfor=False,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        for p in zarr_fullpaths
    ]
    outputs = [f"{d}/{mip_dir}/{name}_z.zarr" for d, name, _ in parts]
    mem = 1024**3 * 4 / 2**30 * 3
    for i in range(3):
        done = framework.run_tasks(
            list(zarr_fullpaths),
            outputs,
            tasks,
            cpus_per_task=1,
            max_trial_num=2,
            parse_cluster=parse_cluster,
            mem_allocate=mem * 2**i,
            min_task_job_num=int(mround(len(zarr_fullpaths) / 2)),
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        if all(done):
            break
    return mip_dir
