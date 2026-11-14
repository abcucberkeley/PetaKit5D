"""Stitching of one volume end to end (``XR_stitching_frame_zarr_dev_v1``): ``stitch_volume``.

The steps, in MATLAB's order:

1. **setup**: the tiles (or a tile information file), the result ``<data>/matlab_stitch/<name>
   .zarr`` (an existing result is kept), the axis order check, the stitching space's voxel sizes
   and its marker file ``px<xp>_py<yp>_pz<zp>``, the edge-erosion function file
   ``<result folder>/tmp/processFunction_<uuid>.mat``;
2. **tile preparation** (``stitch_process_tiles``, S02's): tiff tiles to zarr through T01's
   ``tiff_to_zarr_dataset``, zarr tiles that need processing through ``process_zarr_tile``
   (a stand-in of S02's ``stitch_process_zarr_tile`` until it lands);
3. **geometry** (S03's ``geometry.prepare_layout``): coordinates in stitching space, overlaps;
4. **registration and assignment** of a primary channel (S05's ``assignment.shift_assignment``,
   whose registration tasks are S04's) into ``xcorr/<name>.mat``, or that file if it exists; a
   secondary channel takes the primary's shifts, overlaps and stitched size from its stitch
   information file, and makes its own block information when its tile sizes or its stitched
   size differ from the primary's (with the primary's tile positions);
5. **placement** (S03's ``geometry.place_tiles``): half-overlap regions, canvas, batches; the
   block information (S07's ``blockinfo.process_block_info``);
6. **the temporary stitched zarr** ``<name>.zarr_<uuid>``, resumed when one of the same size and
   chunks exists (a master whose job ran out of time leaves it for the next one);
7. **feather distance maps** (S06's ``distance.compute_tile_distance_transform``) and the stitch
   information file ``stitchInfo/<name>.mat`` of a primary channel;
8. **blending** (S07's ``block.process_stitch_block``) as tasks of the job framework, flags in
   ``zarr_flags/<name>_<uuid>/``, retried three times with twice the CPUs and time each;
   'separate' writes each tile on a canvas of its own instead (``block.process_separate_tiles``);
9. **the z MIP** (merged from the tasks' flags when the volume fits half the memory, else from
   the volume, or by the large-file MIPs), then the result is published by one rename.

Conventions (``docs/conventions.md``; S03's for the geometry): arrays and tile sizes ``(z, y,
x)`` with ``data_order`` naming the tiles' axes (``'z,y,x'``; MATLAB's ``'y,x,z'``); boxes ``(z0,
y0, x0, z1, y1, x1)``, 0-based, end-exclusive (``None`` for MATLAB's ``Inf`` end); sizes and
per-axis vectors ``(z, y, x)`` (``batch_size``, ``block_size``, ``shard_size``, ``pool_size``,
``border_size``, ``stitch_pad_size``, ``xcorr_downsample``, ``axis_weight``, the
``resample_factor`` in the tiles' axis order, the ``stitch_mip`` flags); ``half_order`` axes of
``(z, y, x)`` (MATLAB's default ``[3, 2, 1]`` is ``(0, 1, 2)``). Stage coordinates
``coordinates`` (columns x, y, z, um), ``tile_idx`` (columns x, y, z, location) and
``flipped_tile`` keep MATLAB's values (exception 2).

Files keep MATLAB's names and values, so that MATLAB and the port read each other's (a
secondary channel of one stitched with the primary's information of the other):

- ``xcorr/<name>.mat`` (v7.3): ``xyz_shift``, ``dxyz_shift`` (``(n, 3)``, x, y, z, um, tiles in
  the stitching order);
- ``stitchInfo/<name>.mat`` (v7.3): ``ip``, ``flippedTile``, ``overlap_regions``,
  ``overlap_matrix``, ``dxyz_shift``, ``shiftMethod``, ``pImSz`` (``[ny nx nz]`` in the tiles'
  dimensions), ``pTileSizes``, ``xp``, ``yp``, ``zp``, ``data_order_mat``,
  ``PerBlockInfoFullpaths``, ``imdistFullpaths``, ``imdistFileIdx`` (1-based), ``xyz_orig``
  (``write_stitch_info``, ``read_stitch_info``). MATLAB saves its inputParser object as ``ip``;
  the port saves a struct of its ``Results`` (MATLAB's values) instead.

Not ported: ``saveMultires`` (the multi-resolution result through Python's zarr in MATLAB;
PLAN 4.3). Waiting on other tasks: the registration (S04: without it a primary channel with
``xcorr_shift`` needs its ``xcorr/<name>.mat``), MATLAB-code tile functions such as the edge
erosion (S02's parser) and resampled tiles (T04, N01).
"""

from __future__ import annotations

import math
import os
import shutil
import time
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "StitchVolumeError",
    "StitchVolumeWarning",
    "blend_task_size",
    "pixel_info_name",
    "process_zarr_tile",
    "read_stitch_info",
    "read_xcorr_result",
    "stitch_volume",
    "write_stitch_info",
    "write_xcorr_result",
]

# the variables of the stitch information file, in MATLAB's save order
STITCH_INFO_VARS = (
    "ip",
    "flippedTile",
    "overlap_regions",
    "overlap_matrix",
    "dxyz_shift",
    "shiftMethod",
    "pImSz",
    "pTileSizes",
    "xp",
    "yp",
    "zp",
    "data_order_mat",
    "PerBlockInfoFullpaths",
    "imdistFullpaths",
    "imdistFileIdx",
    "xyz_orig",
)
# the registration's task module (S04's)
REGISTRATION_MODULE = "petakit5d.stitch.registration"
# a volume of more input than this (bytes, as single) is big data (fewer CPUs per task)
BIG_DATA_BYTES = 100 * 2**30


class StitchVolumeError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


class StitchVolumeWarning(UserWarning):
    """A warning MATLAB gives here; ``identifier`` is its warning id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _warn(message: str, identifier: str = "") -> None:
    warnings.warn(StitchVolumeWarning(message, identifier), stacklevel=3)


def _print(msg: str = "", end: str = "\n") -> None:
    print(msg, end=end, flush=True)


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    return np.size(np.asarray(v, dtype=object)) == 0


def _mround(x: float) -> float:
    """MATLAB's round of a double"""
    from petakit5d.core.dtypes import mround

    return float(mround(np.float64(x)))


def _ints(v: Any) -> list[int]:
    return [int(x) for x in np.asarray(v, dtype=np.float64).ravel()]


def _half(v: Any) -> list[int] | None:
    """``round(v / 2)`` of a size (None for MATLAB's ``[]``)"""
    if _empty(v):
        return None
    return [int(_mround(float(x) / 2)) for x in np.asarray(v, dtype=np.float64).ravel()]


def _isdir(path: str) -> bool:
    """MATLAB's ``exist(path, 'dir')``"""
    return os.path.isdir(path)


def _rmtree(path: str) -> None:
    """MATLAB's ``rmdir(path, 's')`` (a missing folder is not an error here)"""
    if os.path.isdir(path):
        shutil.rmtree(path)


def _mkdir_shared(path: str) -> None:
    """MATLAB's ``mkdir`` and ``fileattrib(path, '+w', 'g')`` of a folder other jobs share"""
    from petakit5d.io.tiff import group_writable

    os.makedirs(path, exist_ok=True)
    group_writable(path)


def _msize(sz: Sequence[float]) -> list[float]:
    """a size (z, y, x) in MATLAB's order [d1 d2 d3]"""
    from petakit5d.core.axes import to_yxz

    return [float(v) for v in to_yxz([float(x) for x in sz])]


def _msizes(sizes: np.ndarray) -> np.ndarray:
    """rows of sizes (z, y, x) as MATLAB's rows [d1 d2 d3]"""
    s = np.asarray(sizes, dtype=np.float64).reshape(-1, 3)
    return np.ascontiguousarray(s[:, [1, 2, 0]])


def _psizes(msizes: Any) -> np.ndarray:
    """MATLAB's rows of sizes [d1 d2 d3] as (z, y, x) rows"""
    s = np.asarray(msizes, dtype=np.float64).reshape(-1, 3)
    return np.ascontiguousarray(s[:, [2, 0, 1]])


def _size_text(sz: Sequence[float]) -> str:
    """MATLAB's ``strjoin(arrayfun(@num2str, sz, 'unif', 0), 'x')`` of a size (z, y, x)"""
    from petakit5d.core.text import num2str

    return "x".join(num2str(v) for v in _msize(sz))


def _cell_list(v: Any) -> list:
    """a MATLAB cell (or one char) read from a .mat file, as a list (column-major)"""
    if v is None:
        return []
    if isinstance(v, str):
        return [v] if v else []
    if isinstance(v, np.ndarray):
        return list(v.ravel(order="F"))
    return list(v)


def _column_cell(items: Sequence[Any]) -> Any:
    """a list as MATLAB's column cell (``cell(n, 1)``); ``{}`` (0 x 0) when empty"""
    if not items:
        return []
    a = np.empty((len(items), 1), dtype=object)
    for i, v in enumerate(items):
        a[i, 0] = v
    return a


def _mat_has(path: str, name: str) -> bool:
    """``any(strcmp(who('-file', path), name))``"""
    from petakit5d.io import read_mat

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return name in read_mat(path, [name])


# ---------------------------------------------------------------------------------------------
# the files of a stitched volume
# ---------------------------------------------------------------------------------------------


def write_xcorr_result(path: str, xyz_shift: Any, dxyz_shift: Any) -> None:
    """``save('-v7.3', path, 'xyz_shift', 'dxyz_shift')``: the shifted tile positions and their
    shifts (``(n, 3)``, x, y, z, um, tiles in the stitching order)."""
    from petakit5d.io import write_mat

    write_mat(
        path,
        {
            "xyz_shift": np.asarray(xyz_shift, dtype=np.float64).reshape(-1, 3),
            "dxyz_shift": np.asarray(dxyz_shift, dtype=np.float64).reshape(-1, 3),
        },
    )


def read_xcorr_result(path: str) -> tuple[np.ndarray, np.ndarray]:
    """``xyz_shift`` and ``dxyz_shift`` of a registration result file (MATLAB's or the port's)"""
    from petakit5d.io import read_mat

    a = read_mat(path, ["xyz_shift", "dxyz_shift"])
    return (
        np.asarray(a["xyz_shift"], dtype=np.float64).reshape(-1, 3),
        np.asarray(a["dxyz_shift"], dtype=np.float64).reshape(-1, 3),
    )


def _flipped_value(flipped_tile: Any) -> Any:
    """flippedTile as MATLAB holds it: [] when empty, logical flags or doubles"""
    if _empty(flipped_tile):
        return None
    a = np.asarray(flipped_tile)
    if a.dtype != np.bool_:
        a = a.astype(np.float64)
    return a


def write_stitch_info(
    path: str,
    *,
    params: dict,
    flipped_tile: Any,
    overlap_regions: Any,
    overlap_matrix: Any,
    dxyz_shift: Any,
    shift_method: str,
    p_im_sz: Sequence[float],
    p_tile_sizes: Any,
    xyz_voxel_sizes: Sequence[float],
    data_order: str,
    per_block_info_fullpaths: Sequence[str],
    imdist_fullpaths: Sequence[str],
    imdist_file_idx: Sequence[int] | None,
    xyz_orig: Any,
) -> None:
    """The stitch information file of a primary channel (v7.3, MATLAB's variables and values):
    ``params`` is the call's ``ip.Results`` (MATLAB's values) saved as the struct ``ip``;
    ``p_im_sz`` (``(z, y, x)``) and ``p_tile_sizes`` (``(n, 3)`` in array order) are saved in
    the tiles' dimensions, ``imdist_file_idx`` 1-based."""
    from petakit5d.io import write_mat
    from petakit5d.stitch.geometry import data_order_axes_to_matlab, parse_data_order

    idx = None
    if imdist_file_idx is not None and len(imdist_file_idx) > 0:
        v = np.asarray([float(i) + 1 for i in imdist_file_idx], dtype=np.float64)
        idx = v[0] if v.size == 1 else v.reshape(-1, 1)
    xp, yp, zp = (float(v) for v in xyz_voxel_sizes)
    n = np.asarray(dxyz_shift).reshape(-1, 3).shape[0]
    values = {
        "ip": {"Results": params},
        "flippedTile": _flipped_value(flipped_tile),
        "overlap_regions": np.asarray(overlap_regions, dtype=np.float64).reshape(-1, 6),
        "overlap_matrix": np.asarray(overlap_matrix, dtype=bool).reshape(n, n),
        "dxyz_shift": np.asarray(dxyz_shift, dtype=np.float64).reshape(-1, 3),
        "shiftMethod": str(shift_method),
        "pImSz": np.asarray(_msize(p_im_sz), dtype=np.float64),
        "pTileSizes": _msizes(p_tile_sizes),
        "xp": xp,
        "yp": yp,
        "zp": zp,
        "data_order_mat": np.asarray(
            data_order_axes_to_matlab(parse_data_order(data_order)), dtype=np.float64
        ),
        "PerBlockInfoFullpaths": _column_cell([os.fspath(p) for p in per_block_info_fullpaths]),
        "imdistFullpaths": _column_cell([os.fspath(p) for p in imdist_fullpaths]),
        "imdistFileIdx": idx,
        "xyz_orig": np.asarray(xyz_orig, dtype=np.float64).reshape(-1, 3),
    }
    write_mat(path, values)


def read_stitch_info(path: str, names: Sequence[str] | None = None) -> dict[str, Any]:
    """The variables ``names`` (MATLAB's; all by default) of a stitch information file written
    by MATLAB or the port, in the port's forms, by MATLAB name: ``pImSz`` a size ``(z, y, x)``,
    ``pTileSizes`` ``(n, 3)`` in array order, ``imdistFileIdx`` 0-based ints, ``data_order_mat``
    the array axes of x, y and z (``data_order_axes``), the cells lists of paths,
    ``overlap_matrix`` ``(n, n)`` bool, ``flippedTile`` as MATLAB holds it (None for ``[]``);
    ``ip`` as read (MATLAB's object is not decoded). A variable the file lacks is left out with
    a warning, as MATLAB's ``load`` does."""
    from petakit5d.io import read_mat
    from petakit5d.stitch.geometry import data_order_axes_from_matlab

    a = read_mat(path, None if names is None else list(names))
    out: dict[str, Any] = {}
    for k, v in a.items():
        if k == "pImSz":
            m = np.asarray(v, dtype=np.float64).ravel()
            out[k] = tuple(int(x) for x in (m[2], m[0], m[1])) if m.size == 3 else tuple(m)
        elif k == "pTileSizes":
            out[k] = _psizes(v)
        elif k == "imdistFileIdx":
            out[k] = [int(x) - 1 for x in np.asarray(v, dtype=np.float64).ravel(order="F")]
        elif k == "data_order_mat":
            out[k] = data_order_axes_from_matlab(v)
        elif k in ("PerBlockInfoFullpaths", "imdistFullpaths"):
            out[k] = [str(p) for p in _cell_list(v)]
        elif k == "overlap_matrix":
            b = np.asarray(v, dtype=bool)
            n = math.isqrt(b.size)
            out[k] = b.reshape(n, n, order="F") if b.size else np.zeros((0, 0), bool)
        elif k in ("dxyz_shift", "xyz_orig", "overlap_regions"):
            # (rows; one row is a 1 x c array, a scalar read for none)
            c = 6 if k == "overlap_regions" else 3
            a2 = np.asarray(v, dtype=np.float64)
            out[k] = a2.reshape(-1, c) if a2.size % c == 0 else a2
        elif k == "flippedTile":
            out[k] = None if _empty(v) else np.asarray(v)
        elif k in ("xp", "yp", "zp"):
            out[k] = float(np.asarray(v).ravel()[0])
        elif k == "shiftMethod":
            out[k] = str(v)
        else:
            out[k] = v
    return out


def pixel_info_name(xp: float, yp: float, zp: float) -> str:
    """``sprintf('px%0.5g_py%0.5g_pz%0.5g', xp, yp, zp)``: the name of the marker file of the
    stitching space's voxel sizes"""
    from petakit5d.core.text import sprintf

    return sprintf("px%0.5g_py%0.5g_pz%0.5g", xp, yp, zp)


def _write_pixel_info(result_root: str, name: str) -> None:
    """the marker file of the voxel sizes in ``result_root``; other markers there are deleted
    (MATLAB's ``fclose(fopen(..., 'w'))`` fails where the folder is missing)"""
    from petakit5d.core.paths import matlab_dir
    from petakit5d.io.tiff import group_writable

    for e in matlab_dir(f"{result_root}/px*_py*_pz*"):
        if e.name != name:
            os.remove(f"{e.folder}/{e.name}")
    path = f"{result_root}/{name}"
    if not os.path.exists(path):
        if not os.path.isdir(result_root):
            raise StitchVolumeError(
                "Invalid file identifier. Use fopen to generate a valid file identifier.",
                "MATLAB:FileIO:InvalidFid",
            )
        open(path, "w").close()
        group_writable(path)


def blend_task_size(
    num_batches: float,
    batch_size: Sequence[float],
    parse_cluster: bool,
    local_worker_num: float = 1,
) -> int:
    """The number of batches per blending task: 20 of 512^3 per cluster job (inversely
    proportional to the batch volume, at least 1, and enough for at most 5,000 jobs, within
    ``round(50 * ratio)``); locally all batches, at most 500 of 512^3, and about 4 tasks per
    local worker."""
    ratio = max(0.1, 512.0**3 / math.prod(float(v) for v in batch_size))
    if parse_cluster:
        ts = max(1.0, _mround(20 * ratio))
        ts = max(ts, min(_mround(50 * ratio), math.ceil(num_batches / 5000)))
    else:
        ts = min(num_batches, _mround(500 * ratio))
        if local_worker_num > 1:
            ts = min(ts, math.ceil(num_batches / (4 * local_worker_num)))
    return int(ts)


def _task_order(block_info_bytes: Sequence[float]) -> list[int]:
    """the cluster's task order: the largest block information first (a stable sort), and above
    10,000 tasks each run of four from the four quarters of that order"""
    n = len(block_info_bytes)
    sinds = sorted(range(n), key=lambda i: -float(block_info_bytes[i]))
    if n > 10000:
        m = n // 4 * 4
        q = m // 4
        # reshape(sinds(1 : m), [], 4)' then (:): column k // 4, row k % 4
        head = [sinds[(k % 4) * q + k // 4] for k in range(m)]
        sinds[:m] = head
    return sinds


# ---------------------------------------------------------------------------------------------
# tile preparation (stitch_process_tiles, S02's; a stand-in until it lands)
# ---------------------------------------------------------------------------------------------


def _matlab_code_check(process_fun_path: Any) -> None:
    """MATLAB code in a process function file needs S02's tile-function parser: fail before
    any tile is converted"""
    from petakit5d.io.convert import _fun_rows, _read_function_file, user_function

    for row in _fun_rows(process_fun_path):
        for fn in row:
            if fn:
                fun, _ = _read_function_file(fn)
                if fun:
                    user_function(fun)


def _process_tiles(
    input_fullpaths: Sequence[str],
    *,
    result_dir_name: str,
    zarr_file: bool,
    loc_ids: Sequence[float],
    block_size: Sequence[int],
    shard_size: Sequence[int] | None,
    flipped_tile: Sequence[bool],
    resample: Sequence[float],
    partial_file: bool,
    input_bbox: Any,
    tile_out_bbox: Any,
    process_fun_path: Any,
    parse_cluster: bool,
    master_compute: bool,
    big_data: bool,
    mcc_mode: bool,
    config_file: str,
) -> None:
    """``stitch_process_tiles``: tiff tiles to zarr (``XR_tiffToZarr_wrapper`` with the channel
    pattern '.'), zarr tiles through ``process_zarr_tile`` tasks when they need processing
    (``result_dir_name`` given)."""
    try:  # S02's, once it lands
        from petakit5d.stitch.tiles import process_tiles  # type: ignore[import-not-found]
    except ImportError:
        process_tiles = None
    kw = {
        "result_dir_name": result_dir_name,
        "zarr_file": zarr_file,
        "loc_ids": loc_ids,
        "block_size": block_size,
        "shard_size": shard_size,
        "flipped_tile": flipped_tile,
        "resample": resample,
        "partial_file": partial_file,
        "input_bbox": input_bbox,
        "tile_out_bbox": tile_out_bbox,
        "process_fun_path": process_fun_path,
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "big_data": big_data,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    if process_tiles is not None:
        process_tiles(list(input_fullpaths), **kw)
        return
    _matlab_code_check(process_fun_path)
    if not zarr_file:
        from petakit5d.io.convert import tiff_to_zarr_dataset

        _print("\nConvert tiles from Tiff to Zarr...")
        tiff_to_zarr_dataset(
            "",
            tiff_fullpaths=list(input_fullpaths),
            result_dir_name=result_dir_name,
            loc_ids=list(loc_ids),
            block_size=tuple(block_size),
            shard_size=None if shard_size is None else tuple(shard_size),
            flipped_tile=list(flipped_tile),
            resample_factor=tuple(resample),
            partial_file=partial_file,
            channel_patterns=(".",),
            input_bbox=input_bbox,
            tile_out_bbox=tile_out_bbox,
            process_fun_path=process_fun_path,
            parse_cluster=parse_cluster,
            master_compute=master_compute,
            big_data=big_data,
            cpus_per_task=1,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        return
    if not result_dir_name:
        return
    _process_zarr_tiles(list(input_fullpaths), **kw)


def _process_zarr_tiles(input_fullpaths: list[str], **kw: Any) -> None:
    """the zarr branch of ``stitch_process_tiles``: one ``process_zarr_tile`` task per tile, the
    function of its channel (the pattern '.') by keywords, boxes per location, two rounds"""
    from petakit5d.core.paths import fileparts
    from petakit5d.io import get_image_size
    from petakit5d.io.convert import _box_rows, _fun_rows, _keywords, _read_function_file
    from petakit5d.jobs import framework
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task

    _print("\nPreprocess Zarr tiles...")
    parse_cluster = bool(kw["parse_cluster"])
    data_path = fileparts(input_fullpaths[0])[0]
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(data_path, "../job_logs").parse_cluster)
    rows = _fun_rows(kw["process_fun_path"])
    n_cols = max(1, max((len(r) for r in rows), default=0))
    usr = [[""] * n_cols for _ in range(max(1, len(rows)))]
    kws = [[[""] for _ in range(n_cols)] for _ in range(max(1, len(rows)))]
    _print("Process function paths:")
    for r in rows:
        _print("    " + "    ".join(r))
    for i, r in enumerate(rows):
        for j, fn in enumerate(r):
            if not fn:
                continue
            fun, k = _read_function_file(fn)
            if fun:
                usr[i][j] = fun
            if k is not None:
                kws[i][j] = _keywords(k)
    in_rows, out_rows = _box_rows(kw["input_bbox"]), _box_rows(kw["tile_out_bbox"])
    loc = [float(v) for v in kw["loc_ids"]]
    uniq = sorted(set(loc))
    loc_specific = len(uniq) > 1 and (len(in_rows) == len(uniq) or len(out_rows) == len(uniq))
    if loc_specific:
        for what, r in (("InputBbox", in_rows), ("tileOutBbox", out_rows)):
            if len(r) > 1 and len(r) != len(uniq):
                raise StitchVolumeError(
                    f"The number of {what} ({len(r)}) does not match the number of tile "
                    f"locations ({len(uniq)})!"
                )
        in_rows = in_rows if len(in_rows) > 1 else [in_rows[0] if in_rows else None] * len(uniq)
        out_rows = (
            out_rows if len(out_rows) > 1 else [out_rows[0] if out_rows else None] * len(uniq)
        )
    flips = list(kw["flipped_tile"])
    resample = [float(v) for v in kw["resample"]]
    outs, tasks = [], []
    for i, p in enumerate(input_fullpaths):
        dpath, fsn, _ = fileparts(p)
        zpath = f"{dpath}/{kw['result_dir_name']}"
        if not os.path.isdir(zpath):
            _mkdir_shared(zpath)
        outs.append(f"{dpath}/{kw['result_dir_name']}/{fsn}.zarr/")
        if loc_specific:
            k = uniq.index(loc[i])
            in_i, out_i = in_rows[k], out_rows[k]
        else:
            in_i = in_rows[0] if in_rows else None
            out_i = out_rows[0] if out_rows else None
        # the channel's functions ('.' matches every tile), column-major, picked by keywords
        strs = [usr[0][j] for j in range(n_cols)]
        kw_j = [kws[0][j] for j in range(n_cols)]
        usr_i = strs[0]
        for j in range(len(strs)):
            if any(k in p for k in kw_j[j]):
                usr_i = strs[j]
                break
        tasks.append(
            Task.of(
                process_zarr_tile,
                [p],
                outs[-1],
                None,
                block_size=list(kw["block_size"]),
                shard_size=None if kw["shard_size"] is None else list(kw["shard_size"]),
                flip_zstack=bool(flips[i]) if flips else False,
                resample=resample,
                input_bbox=None if in_i is None else list(in_i),
                tile_out_bbox=None if out_i is None else list(out_i),
                compressor="zstd",
                usr_fcn=usr_i,
            )
        )
    sz = get_image_size(input_fullpaths[-1])
    mem = float(sz[1]) * float(sz[2]) * float(sz[0]) * 4 / 1024**3 * 2.5
    if not kw["big_data"]:
        mem *= 2
    common = {
        "parse_cluster": parse_cluster,
        "master_compute": bool(kw["master_compute"]),
        "max_trial_num": 2,
        "job_log_dir": "../job_logs",
        "mcc_mode": kw["mcc_mode"],
        "config_file": kw["config_file"],
    }
    done = framework.run_tasks(
        input_fullpaths, outs, tasks, cpus_per_task=1, mem_allocate=mem, **common
    )
    if not all(done):
        framework.run_tasks(
            input_fullpaths, outs, tasks, cpus_per_task=2, mem_allocate=mem * 2, **common
        )


def process_zarr_tile(
    input_filename: str | Sequence[str],
    zarr_filename: str,
    frame: np.ndarray | None = None,
    *,
    overwrite: bool = False,
    block_size: Sequence[int] = (250, 500, 500),
    shard_size: Sequence[int] | None = None,
    expand_2d_dim: bool = True,
    flip_zstack: bool = False,
    resample: Sequence[float] | None = None,
    input_bbox: Sequence[int] | None = None,
    tile_out_bbox: Sequence[int] | None = None,
    compressor: str = "lz4",
    usr_fcn: Any = "",
    uuid: str = "",
) -> None:
    """A stand-in of S02's ``stitch_process_zarr_tile``: a zarr tile (its box ``input_bbox``)
    flipped in z, cropped to ``tile_out_bbox``, passed through the user function and cast back
    to its class, as a native zarr (an existing result is kept). Resampling (MATLAB's
    ``imresize3``) is not ported yet."""
    from petakit5d.core.arrays import matlab_class
    from petakit5d.core.dtypes import cast
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import create_zarr, read_zarr, write_zarr
    from petakit5d.io.convert import user_function
    from petakit5d.jobs.publish import publish_result

    zarr_filename = os.fspath(zarr_filename).rstrip("/")
    if os.path.isdir(zarr_filename):
        if not overwrite:
            _print(f"Zarr result {zarr_filename} already exists, skip it!")
            return
        shutil.rmtree(zarr_filename)
    if not uuid:
        uuid = get_uuid()
    t0 = time.perf_counter()
    first = input_filename if isinstance(input_filename, str) else list(input_filename)[0]
    _print(f"Process zarr tile for {first if first else zarr_filename} ...")
    block = [int(v) for v in block_size]
    if frame is None or np.size(frame) == 0:
        frame = read_zarr(first, input_bbox=None if _empty(input_bbox) else list(input_bbox))
    dtype = matlab_class(frame)
    if frame.ndim == 2:
        frame = frame.reshape(1, *frame.shape)
    if any(s == 0 for s in frame.shape):
        _print("The input image is empty, skip it!")
        return
    if flip_zstack:
        frame = np.ascontiguousarray(frame[::-1])
    if not _empty(tile_out_bbox):
        from petakit5d.core.bbox import box_slices

        frame = np.ascontiguousarray(frame[box_slices(list(tile_out_bbox))])
    rs = [] if _empty(resample) else [float(v) for v in np.asarray(resample).ravel()]
    if rs and not all(v == 1 for v in rs):
        err = NotImplementedError(
            "resampled zarr tiles need imresize3 (task N01) and S02's stitch_process_zarr_tile"
        )
        err.identifier = "petakit5d:stitch:resampleNotPorted"  # type: ignore[attr-defined]
        raise err
    fn = user_function(usr_fcn)
    if fn is not None:
        _print(usr_fcn if isinstance(usr_fcn, str) else repr(usr_fcn))
        frame = fn(frame)
    sz = list(frame.shape)
    block = [min(a, b) for a, b in zip(sz, block[-3:], strict=True)]
    tmp = f"{zarr_filename}_{uuid}"
    if not os.path.isdir(tmp):
        sep = (
            "/"
            if math.prod(math.ceil(a / b) for a, b in zip(sz, block, strict=True)) > 10000
            else "."
        )
        create_zarr(
            tmp,
            data_size=sz,
            block_size=block,
            shard_size=shard_size,
            dtype=dtype,
            expand_2d_dim=expand_2d_dim,
            compressor=compressor,
            dim_separator=sep,
        )
    write_zarr(np.ascontiguousarray(cast(frame, dtype)), tmp, create=False)
    publish_result(tmp, zarr_filename)
    _print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.")
    _print("Done!")


# ---------------------------------------------------------------------------------------------
# the frame
# ---------------------------------------------------------------------------------------------


def _bbox_rows_in(v: Any) -> Any:
    """one box, or one box per tile location (rows)"""
    from petakit5d.io.convert import _bbox_rows_in as f

    return f(v)


def _bbox_rows_out(v: Any) -> Any:
    from petakit5d.io.convert import _bbox_rows_out as f

    return f(v)


def _out_bbox_in(v: Any) -> Any:
    """outBbox: a box (or MATLAB's 3 x 2 [starts ends]) as a Python box"""
    from petakit5d.core.bbox import box_from_matlab

    if _empty(v):
        return None
    return box_from_matlab(np.asarray(v, dtype=np.float64).ravel(order="F"))


def _out_bbox_out(v: Any) -> Any:
    from petakit5d.core.bbox import box_to_matlab

    return None if _empty(v) else box_to_matlab(v)


def _half_order_in(v: Any) -> Any:
    from petakit5d.stitch.geometry import half_order_from_matlab

    return None if _empty(v) else half_order_from_matlab(v)


def _half_order_out(v: Any) -> Any:
    from petakit5d.stitch.geometry import half_order_to_matlab

    return None if _empty(v) else half_order_to_matlab(v)


def _stitch_mip_flags(stitch_mip: Any) -> tuple[bool, bool, bool] | None:
    """stitchMIP expanded to three flags (z, y, x): one value is z; of several the first in
    MATLAB's order (y, x, z) is kept"""
    if _empty(stitch_mip):
        return None
    v = [bool(x) for x in np.asarray(stitch_mip).ravel()]
    if not any(v):
        return None
    if len(v) == 1:
        return (v[0], False, False)
    f = list(v)
    if sum(f) > 1:
        first = next(a for a in (1, 2, 0) if f[a])
        f = [a == first for a in range(3)]
    return (f[0], f[1], f[2])


def _edge_function(dsr: bool, stitch_mip, tile_offset: float, edge_artifacts: float) -> str:
    """the edge-erosion function of the tiles (MATLAB code, as MATLAB writes it)"""
    from petakit5d.core.text import sprintf

    if dsr:
        if stitch_mip is None:
            return sprintf("@(x)erodeVolumeBy2DProjection(x,%d)", edge_artifacts)
        return ""
    # (skewed space: TileOffset counts are added, so that no data voxel is 0)
    if tile_offset != 0:
        return sprintf("@(x)erodeVolumeBy2DProjection(x+%d,%d)", tile_offset, edge_artifacts)
    return sprintf("@(x)erodeVolumeBy2DProjection(x,%d)", edge_artifacts)


def _process_fun_rows(process_fun_path: Any) -> list[list[str]] | None:
    """processFunPath as MATLAB's cell: rows of function files (None: MATLAB's ``{''}``); a flat
    sequence is one row (the wrapper passes ``1 x k``)"""
    if _empty(process_fun_path):
        return None
    if isinstance(process_fun_path, str):
        # MATLAB: cellfun(@isempty, processFunPath) of a char fails
        raise StitchVolumeError(
            "Input #2 expected to be a cell array, was char instead.", "MATLAB:cellfun:NotACell"
        )
    items = list(process_fun_path)
    if items and all(isinstance(x, str) for x in items):
        rows = [[str(x) for x in items]]
    else:
        rows = [[str(x) for x in r] for r in items]
    if all(not f for r in rows for f in r):
        return None
    return rows


def _check_registration() -> None:
    import importlib.util

    if importlib.util.find_spec(REGISTRATION_MODULE) is None:
        err = NotImplementedError(
            "the registration of the tiles (cross_correlation_registration_wrapper, task S04) "
            "is not ported yet: give the volume's registration result xcorr/<name>.mat, or "
            "set xcorr_shift=False"
        )
        err.identifier = "petakit5d:stitch:registrationNotPorted"  # type: ignore[attr-defined]
        raise err


@register(
    "XR_stitching_frame_zarr_dev_v1",
    Param(
        "tileFullpaths",
        "@iscell",
        "The tiles (tiff or zarr files; raw names with processed tiles).",
    ),
    Param(
        "coordinates",
        "@isnumeric",
        "Stage coordinates (x, y, z, um), one row per tile.",
        convert="same",
    ),
    Param("ResultPath", "@ischar", "The stitched zarr ('': <data>/matlab_stitch/<name>.zarr)."),
    Param(
        "tileInfoFullpath",
        "@ischar",
        "A .mat file with the tiles, coordinates, tile indices and flips.",
    ),
    Param("stitchInfoDir", "@ischar", "Folder of the stitch information in the result folder."),
    Param(
        "stitchInfoFullpath", "@ischar", "The primary channel's stitch information (secondaries)."
    ),
    Param(
        "ProcessedDirstr", "@ischar", "Folder of processed tiles in the data folder (e.g. 'DSR')."
    ),
    Param("Overwrite", "@islogical", "Not used."),
    Param("SkewAngle", "@isscalar", "Skew angle (degrees)."),
    Param(
        "axisOrder",
        "@ischar",
        "Axis order of the stage coordinates ('x,y,z', '-x,y,z', ...).",
        convert="same",
    ),
    Param("dataOrder", "@ischar", "Axes of the tile arrays ('z,y,x').", convert="data_order"),
    Param(
        "flippedTile",
        "@(x) isempty(x) || all(islogical(x) | isnumeric(x))",
        "Per tile: flipped (bidirectional scan).",
        convert="same",
    ),
    Param("dz", "@isscalar", "Scan step (um)."),
    Param("xyPixelSize", "@isscalar", "Pixel size (um)."),
    Param("objectiveScan", "@islogical", "Objective scan."),
    Param("IOScan", "@islogical", "Scan of the sample stage along the objective (no deskew)."),
    Param("sCMOSCameraFlip", "@islogical", "Not used."),
    Param("Reverse", "@islogical", "Not used."),
    Param("Crop", "@islogical", "Not used."),
    Param(
        "InputBbox",
        "@isnumeric",
        "Box of each tile read (z0, y0, x0, z1, y1, x1) (or one per tile location).",
        convert=(_bbox_rows_in, _bbox_rows_out),
    ),
    Param(
        "tileOutBbox",
        "@isnumeric",
        "Box of each tile kept after the flip (or one per tile location).",
        convert=(_bbox_rows_in, _bbox_rows_out),
    ),
    Param(
        "TileOffset",
        "@isnumeric",
        "Counts added to the tiles (skewed space, edge erosion).",
        convert="same",
    ),
    Param("df", "@isnumeric", "Not used."),
    Param("save16bit", "@islogical", "Save the result as uint16.", convert="same"),
    Param("EdgeArtifacts", "@isnumeric", "Edge erosion of the tiles (voxels; 0: none)."),
    Param("Decon", "@islogical", "Not used."),
    Param("DS", "@islogical", "Tiles deskewed."),
    Param("DSR", "@islogical", "Tiles deskewed and rotated."),
    Param(
        "resampleType",
        "@ischar",
        "'xy_isotropic', 'isotropic', 'method_1' or 'given'.",
        convert="same",
    ),
    Param(
        "resampleFactor",
        "@isnumeric",
        "Resampling of the tiles (z, y, x), in the tiles' axis order (None: none).",
        convert="resample_yxz",
    ),
    Param("deconRotate", "@islogical", "Not used."),
    Param("BlendMethod", "@ischar", "'none', 'feather', 'max', 'mean', 'median' or 'separate'."),
    Param("blendWeightDegree", "@isnumeric", "Power of the feather weights."),
    Param(
        "halfOrder",
        "@isnumeric",
        "Order of the axes (z, y, x) in which overlaps are split.",
        convert=(_half_order_in, _half_order_out),
    ),
    Param("overlapType", "@ischar", "'' (by the blend method), 'none', 'half', 'full'."),
    Param("xcorrShift", "@islogical", "Register the tiles."),
    Param(
        "isPrimaryCh", "@islogical", "Primary channel (registered); else the primary's information."
    ),
    Param(
        "usePrimaryCoords",
        "@islogical",
        "A secondary takes the primary's stage coordinates.",
        convert="same",
    ),
    Param(
        "stitchPadSize",
        "@(x) isnumeric(x) && numel(x) == 3",
        "Padding per tile along each axis (z, y, x) with registration.",
        convert="size_yxz",
    ),
    Param(
        "outBbox",
        "@(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6)",
        "Box of the stitched volume kept (None ends: to the end).",
        convert=(_out_bbox_in, _out_bbox_out),
    ),
    Param("zNormalize", "@islogical", "Not used."),
    Param(
        "xcorrDownsample",
        "@isnumeric",
        "Downsampling of the registration (z, y, x).",
        convert="size_yxz",
    ),
    Param("xcorrThresh", "@isnumeric", "Pairs below this correlation are not used."),
    Param("xyMaxOffset", "@isnumeric", "Largest shift in x and y (voxels)."),
    Param("zMaxOffset", "@isnumeric", "Largest shift in z (voxels)."),
    Param("shiftMethod", "@ischar", "'grid', 'local', 'global', 'group' or 'test'."),
    Param(
        "axisWeight",
        "@isnumeric",
        "Weights of the tile index axes (z, y, x) in the assignment.",
        convert="size_yxz",
    ),
    Param("groupFile", "@ischar", "Groups of tiles (csv) for 'group'."),
    Param("singleDistMap", "@islogical", "One distance map for the tiles of a location."),
    Param("distBboxes", "@isnumeric", "Boxes of the distance maps.", convert="bboxes"),
    Param("zarrFile", "@islogical", "Tiles are zarr."),
    Param("largeFile", "@islogical", "Large tiles (registration and maps on MIP slabs)."),
    Param(
        "poolSize",
        "@isnumeric",
        "Pooling of the MIP slabs (z, y, x) for large files.",
        convert="size_yxz",
    ),
    Param("blockSize", "@isnumeric", "Chunk size of the result (z, y, x).", convert="size_yxz"),
    Param("batchSize", "@isnumeric", "Batch size of the blending (z, y, x).", convert="size_yxz"),
    Param("shardSize", "@isnumeric", "Shard size (z, y, x) (Zarr v3).", convert="size_yxz"),
    Param("saveMultires", "@islogical", "Not ported (a multi-resolution result)."),
    Param("resLevel", "@isnumeric", "Not ported (levels of the multi-resolution result)."),
    Param(
        "BorderSize", "@isnumeric", "Border of the blending batches (z, y, x).", convert="size_yxz"
    ),
    Param("saveMIP", "@islogical", "Save the z MIP of the result."),
    Param(
        "tileIdx",
        "@isnumeric",
        "Tile indices (x, y, z, location), one row per tile.",
        convert="same",
    ),
    Param(
        "processFunPath",
        "@(x) isempty(x) || iscell(x) || ischar(x)",
        "Function files (.mat or .txt) applied to the tiles.",
    ),
    Param(
        "stitchMIP",
        "@(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3))",
        "Stitch the tiles' MIPs (flags z, y, x; one value: z).",
        convert="flags_yxz",
    ),
    Param("stitch2D", "@(x)islogical(x)", "Not used (2D is decided from the tiles)."),
    Param("bigStitchData", "@(x)islogical(x)", "Big data: one CPU per blending task."),
    Param(
        "maxFileNumPerFolder", "@(x)isscalar(x)", "Files per folder of block information and flags."
    ),
    Param("parseCluster", "@islogical", "Run the tasks as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "Not used."),
    Param("uuid", "@ischar", "Id of the temporary results ('': a new one)."),
    Param("maxTrialNum", "@isnumeric", "Not used (two tries per round)."),
    Param("unitWaitTime", "@isnumeric", "Not used."),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only)."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param("localWorkerNum", "@isnumeric", "Local worker processes without the cluster."),
    Param("debug", "@islogical", "Keep the temporary results and flags."),
    task="S08",
)
def stitch_volume(
    tile_fullpaths: Sequence[str],
    coordinates: Any,
    *,
    result_path: str = "",
    tile_info_fullpath: str = "",
    stitch_info_dir: str = "stitchInfo",
    stitch_info_fullpath: str = "",
    processed_dirstr: str = "",
    overwrite: bool = False,
    skew_angle: float = 32.45,
    axis_order: str = "x,y,z",
    data_order: str = "z,y,x",
    flipped_tile: Any = None,
    dz: float = 0.5,
    xy_pixel_size: float = 0.108,
    objective_scan: bool = False,
    io_scan: bool = False,
    s_cmos_camera_flip: bool = False,
    reverse: bool = False,
    crop: bool = False,
    input_bbox: Any = None,
    tile_out_bbox: Any = None,
    tile_offset: float = 0,
    df: Any = None,
    save_16bit: bool = True,
    edge_artifacts: float = 0,
    decon: bool = False,
    ds: bool = False,
    dsr: bool = False,
    resample_type: str = "xy_isotropic",
    resample_factor: Sequence[float] | None = None,
    decon_rotate: bool = False,
    blend_method: str = "none",
    blend_weight_degree: float = 10,
    half_order: Sequence[int] = (0, 1, 2),
    overlap_type: str = "",
    xcorr_shift: bool = True,
    is_primary_ch: bool = True,
    use_primary_coords: bool = False,
    stitch_pad_size: Sequence[float] = (2, 2, 2),
    out_bbox: Sequence[int | None] | None = None,
    z_normalize: bool = False,
    xcorr_downsample: Sequence[float] = (1, 2, 2),
    xcorr_thresh: float = 0.25,
    xy_max_offset: float = 300,
    z_max_offset: float = 50,
    shift_method: str = "grid",
    axis_weight: Sequence[float] = (10, 1, 0.1),
    group_file: str = "",
    single_dist_map: bool = True,
    dist_bboxes: Any = None,
    zarr_file: bool = False,
    large_file: bool = False,
    pool_size: Sequence[float] | None = None,
    block_size: Sequence[int] = (500, 500, 500),
    batch_size: Sequence[int] = (500, 500, 500),
    shard_size: Sequence[int] | None = None,
    save_multires: bool = False,
    res_level: float = 4,
    border_size: Sequence[int] = (0, 0, 0),
    save_mip: bool = True,
    tile_idx: Any = None,
    process_fun_path: Any = "",
    stitch_mip: Any = None,
    stitch_2d: bool = False,
    big_stitch_data: bool = False,
    max_file_num_per_folder: float = 20000,
    parse_cluster: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 2,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 30,
    mcc_mode: bool = False,
    config_file: str = "",
    local_worker_num: float = 1,
    debug: bool = False,
) -> None:
    """Stitch one volume from its tiles (``XR_stitching_frame_zarr_dev_v1``; the steps and files
    are in the module's docstring): the tiles (``tile_fullpaths``, at stage ``coordinates``)
    are prepared, registered (a primary channel with ``xcorr_shift``) or placed with a primary
    channel's shifts (``is_primary_ch`` False, ``stitch_info_fullpath``), blended by
    ``blend_method`` in batches and published as ``result_path`` (by default
    ``<data>/matlab_stitch/<name>.zarr``, the first tile's name without its last 21
    characters), with its z MIP in ``MIPs/``. An existing result is kept."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    if save_multires:
        raise NotImplementedError(
            "saveMultires (a multi-resolution result) is not ported (PLAN.md 4.3)"
        )
    # (arrays as nested lists: the registry's matlab_results compares values with ())
    values = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in locals().items()}
    params = {k: to_plain(v) for k, v in matlab_results(stitch_volume, values).items()}

    tile_fullpaths = [os.fspath(p) for p in tile_fullpaths]
    # tile paths and coordinates from a file
    if tile_info_fullpath:
        from petakit5d.io import read_mat

        a = read_mat(tile_info_fullpath)
        tile_fullpaths = [str(p) for p in _cell_list(a["tile_fullpaths"])]
        coordinates = np.asarray(a["xyz"], dtype=np.float64)
        tile_idx = a["tileIdx"]
        flipped_tile = a["flippedTile"]
    data_path, fsname_first, _ = fileparts(tile_fullpaths[0])
    if not result_path:
        result_path = f"{data_path}/matlab_stitch/{fsname_first[:-21]}.zarr"
    n_f = len(tile_fullpaths)

    mip_flags = _stitch_mip_flags(stitch_mip)
    # a 2D (MIP) stitch of fewer than 50 tiles runs in this process
    if mip_flags is not None and n_f < 50:
        parse_cluster = False
    if parse_cluster:
        from petakit5d.jobs.slurm import check_slurm_cluster

        parse_cluster = bool(check_slurm_cluster(data_path).parse_cluster)
    cluster_mode = parse_cluster
    if not uuid:
        uuid = get_uuid()

    nv_fullname = result_path
    pstr, nv_fsname, _ = fileparts(result_path)
    result_dir = fileparts(pstr)[1]
    result_root = f"{data_path}/{result_dir}"
    if _isdir(nv_fullname):
        _print(f"Stitched result {nv_fullname} is already exist, skip it!")
        return
    _print("Stitch for tiles:\n" + "\n".join(tile_fullpaths))
    if io_scan:
        dsr = False

    from petakit5d.stitch import geometry as G

    G.check_axis_order(axis_order)
    dop = G.parse_data_order(data_order)
    xyz_override = None

    # a secondary channel or time point takes the primary's flipped tiles (and coordinates)
    if not is_primary_ch:
        if not os.path.exists(stitch_info_fullpath):
            raise StitchVolumeError(
                f"The stitch information filename {stitch_info_fullpath} does not exist!"
            )
        if _empty(flipped_tile):
            flipped_tile = read_stitch_info(stitch_info_fullpath, ["flippedTile"])["flippedTile"]
        if use_primary_coords:
            xyz_override = read_stitch_info(stitch_info_fullpath, ["xyz_orig"])["xyz_orig"]

    tidx_in = None if _empty(tile_idx) else np.asarray(tile_idx, dtype=np.float64).reshape(n_f, -1)
    if tidx_in is not None:
        tidx_sorted, sinds = G.sort_tiles(tidx_in)
    else:
        tidx_sorted, sinds = np.zeros((n_f, 4)), np.arange(n_f, dtype=np.intp)
    tile_fullpaths = [tile_fullpaths[k] for k in sinds]
    flipped_sorted = None if _empty(flipped_tile) else np.asarray(flipped_tile).ravel()
    if flipped_sorted is not None and flipped_sorted.size == n_f:
        flipped_sorted = flipped_sorted[sinds]

    space = G.stitch_pixel_sizes(
        xy_pixel_size,
        dz,
        skew_angle,
        resample_type=resample_type,
        resample_factor=resample_factor,
        objective_scan=objective_scan,
        io_scan=io_scan,
        dsr=dsr,
        data_order=data_order,
    )
    xp, yp, zp = space.xp, space.yp, space.zp

    # the marker file of the voxel sizes
    _write_pixel_info(result_root, pixel_info_name(xp, yp, zp))

    if not ds and not dsr and mip_flags is None:
        single_dist_map = True

    fun_rows = _process_fun_rows(process_fun_path)
    if fun_rows is None:
        fun_rows = [[""]]
        if edge_artifacts > 0:
            from petakit5d.io import write_mat

            usr_fun = _edge_function(dsr, mip_flags, float(tile_offset), float(edge_artifacts))
            # a name of its own: the jobs of other volumes write to the same folder at the same
            # time (the run's uuid is shared by them)
            tmp_dir = f"{pstr}/tmp/"
            fn = f"{tmp_dir}/processFunction_{get_uuid()}.mat"
            os.makedirs(tmp_dir, exist_ok=True)
            write_mat(fn, {"usrFun": usr_fun})
            if usr_fun:
                fun_rows = [[fn]]

    # DSR tiles are stitched as they are (resampling happens later in the pipeline)
    rf = (1.0, 1.0, 1.0) if _empty(resample_factor) else tuple(float(v) for v in resample_factor)
    zarr_flipped: list[bool] = [False] * (0 if flipped_sorted is None else flipped_sorted.size)
    stitch_resample = (1.0, 1.0, 1.0)
    if not dsr:
        stitch_resample = rf
        if not processed_dirstr and flipped_sorted is not None:
            zarr_flipped = [bool(v > 0) for v in flipped_sorted.astype(np.float64)]
    # partial files of tiles in skewed space
    partial_file = not ds and not dsr
    has_fun = any(f for r in fun_rows for f in r)
    process_tiles = not zarr_file or (
        not _empty(input_bbox)
        or not _empty(tile_out_bbox)
        or has_fun
        or any(zarr_flipped)
        or any(v != 1 for v in stitch_resample)
    )
    from petakit5d.core.axes import to_yxz
    from petakit5d.stitch.imagelist import process_filenames

    input_fullpaths, zarr_fullpaths, _, zarr_pathstr = process_filenames(
        tile_fullpaths,
        processed_dirstr,
        mip_flags is not None,
        [float(v) for v in to_yxz(list(stitch_resample))],
        zarr_file,
        process_tiles,
    )

    # 2D tiles (the first 3D tile ends the check; its size decides big data)
    from petakit5d.io import get_image_dtype, get_image_size, get_image_size_batch

    im_size: tuple[int, ...] = ()
    for f in range(n_f):
        im_size = tuple(get_image_size(input_fullpaths[f]))
        s3 = im_size if len(im_size) >= 3 else (1, *im_size)
        if s3[dop[2]] > 1:
            stitch_2d = False
            break
        stitch_2d = True
    if stitch_2d and n_f < 50:
        parse_cluster = False

    if large_file:
        if _empty(pool_size):
            raise StitchVolumeError("poolSize must not be empty when largeFile is true!")
    else:
        pool_size = None

    node_factor = 2
    compressor = "zstd"
    if (
        not stitch_2d
        and not big_stitch_data
        and n_f > 4
        and math.prod(float(v) for v in im_size) * n_f * 4 > BIG_DATA_BYTES
    ) or large_file:
        big_stitch_data = True
    if big_stitch_data:
        node_factor = 1

    loc_ids = [float(v) for v in tidx_sorted[:, 3]]
    _process_tiles(
        input_fullpaths,
        result_dir_name=zarr_pathstr,
        zarr_file=zarr_file,
        loc_ids=loc_ids,
        block_size=_half(block_size) or [],
        shard_size=_half(shard_size),
        flipped_tile=zarr_flipped,
        resample=stitch_resample,
        partial_file=partial_file,
        input_bbox=input_bbox,
        tile_out_bbox=tile_out_bbox,
        process_fun_path=fun_rows,
        parse_cluster=parse_cluster,
        master_compute=master_compute,
        big_data=big_stitch_data,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )

    # the processed tiles' sizes (z, y, x)
    im_sizes = np.asarray(get_image_size_batch(zarr_fullpaths), dtype=np.float64)
    if mip_flags is not None and im_sizes.shape[1] == 2:
        im_sizes = np.column_stack([np.ones(n_f), im_sizes])
    dtype = get_image_dtype(zarr_fullpaths[0])
    if save_16bit:
        dtype = "uint16"

    # a secondary channel or time point stitches with the primary's per-block information,
    # which is made for the primary's tile sizes and stitched size. If its tiles (e.g., a tile
    # whose acquisition stopped early) or its stitched size differ, it makes its own, with the
    # primary's tile positions (its coordinates, here, and its shifts) and task size
    own_block_info = False
    coord_sizes = None
    if not is_primary_ch:
        p_info = read_stitch_info(stitch_info_fullpath, ["pImSz", "pTileSizes"])
        p_im_sz = tuple(float(v) for v in p_info["pImSz"])
        p_tile_sizes = p_info["pTileSizes"]
        sv = list(p_im_sz)
        if out_bbox is not None and not _empty(out_bbox):
            b = list(out_bbox)
            sv = [(sv[a] if b[3 + a] is None else float(b[3 + a])) - float(b[a]) for a in range(3)]
        diff_tiles: list[int] = []
        if im_sizes.shape == p_tile_sizes.shape:
            diff_tiles = [int(k) for k in np.nonzero(np.any(im_sizes != p_tile_sizes, axis=1))[0]]
        if diff_tiles:
            texts = [
                f"Tile {k + 1} ({fileparts(tile_fullpaths[k])[1]}) is {_size_text(im_sizes[k])}, "
                f"the primary's {_size_text(p_tile_sizes[k])}"
                for k in diff_tiles
            ]
            _warn(
                "; ".join(texts) + ": so this volume gets its own block information, with the "
                "primary's tile positions.",
                "XR_stitching_frame_zarr_dev_v1:tileSizes",
            )
        if list(sv) != list(p_im_sz):
            _warn(
                f"The stitched size {_size_text(sv)} differs from the primary's "
                f"{_size_text(p_im_sz)}, so this volume gets its own block information, with "
                "the primary's tile positions.",
                "XR_stitching_frame_zarr_dev_v1:stitchedSize",
            )
        own_block_info = bool(diff_tiles) or list(sv) != list(p_im_sz)
        if own_block_info and _mat_has(stitch_info_fullpath, "xyz_orig"):
            xyz_override = read_stitch_info(stitch_info_fullpath, ["xyz_orig"])["xyz_orig"]
            coord_sizes = p_tile_sizes

    # the geometry: coordinates in stitching space, cuboids, overlaps (im_sizes back in the
    # listed order, as prepare_layout sorts them with the tiles)
    listed = np.empty_like(im_sizes)
    listed[sinds] = im_sizes
    lay = G.prepare_layout(
        coordinates,
        listed,
        tile_idx=tidx_in,
        flipped_tile=None if _empty(flipped_tile) else flipped_tile,
        axis_order=axis_order,
        data_order=data_order,
        xy_pixel_size=xy_pixel_size,
        dz=dz,
        skew_angle=skew_angle,
        objective_scan=objective_scan,
        io_scan=io_scan,
        ds=ds,
        dsr=dsr,
        resample_type=resample_type,
        resample_factor=resample_factor,
        stitch_2d=stitch_2d,
        coord_sizes=coord_sizes,
        xyz_orig=xyz_override,
    )
    stitch_2d = lay.stitch_2d
    vs = (xp, yp, zp)

    # the shifts: registration of a primary channel, or the primary's
    dxyz_shift = np.zeros((n_f, 3))
    xyz_shift = lay.xyz
    p_im_sz_pl = None
    imdist_fullpaths: list[str] = []
    imdist_file_idx: list[int] = []
    p_tile_sizes_si = None
    if xcorr_shift and is_primary_ch:
        xcorr_dir = f"{result_root}/xcorr/{nv_fsname}/"
        xcorr_final = f"{result_root}/xcorr/{nv_fsname}.mat"
        if not os.path.exists(xcorr_final):
            from petakit5d.stitch.assignment import shift_assignment

            if np.any(lay.overlap_matrix):
                _check_registration()
            xcorr_tmp = f"{result_root}/xcorr/{nv_fsname}_{uuid}.mat"
            xyz_shift, dxyz_shift = shift_assignment(
                zarr_fullpaths,
                xcorr_dir,
                lay.im_sizes,
                lay.xyz,
                vs,
                data_order,
                lay.overlap_matrix,
                lay.overlap_regions,
                (z_max_offset, xy_max_offset, xy_max_offset),
                xcorr_downsample,
                xcorr_thresh,
                lay.tile_idx,
                shift_method,
                stitch_2d,
                axis_weight,
                group_file,
                large_file,
                pool_size,
                parse_cluster,
                node_factor,
                mcc_mode,
                config_file,
                int(local_worker_num),
            )
            write_xcorr_result(xcorr_tmp, xyz_shift, dxyz_shift)
            os.replace(xcorr_tmp, xcorr_final)
        else:
            xyz_shift, dxyz_shift = read_xcorr_result(xcorr_final)
    elif not is_primary_ch:
        si = read_stitch_info(
            stitch_info_fullpath,
            ["overlap_matrix", "dxyz_shift", "pImSz", "pTileSizes", "imdistFullpaths",
             "imdistFileIdx"],
        )  # fmt: skip
        om_p = si["overlap_matrix"]
        om = lay.overlap_matrix
        if om_p.shape != om.shape:
            raise StitchVolumeError(
                "Arrays have incompatible sizes for this operation.",
                "MATLAB:sizeDimensionsMustMatch",
            )
        if not np.all(om_p == om):
            _warn(
                "The overlap matrix of current overlap matrix is different from the one in the "
                "stitch info, use the common ones!"
            )
            lay.overlap_matrix = om & om_p
        dxyz_shift = np.asarray(si["dxyz_shift"], dtype=np.float64).reshape(-1, 3)
        xyz_shift = lay.xyz + dxyz_shift
        p_im_sz_pl = si["pImSz"]
        p_tile_sizes_si = si["pTileSizes"]
        imdist_fullpaths = list(si.get("imdistFullpaths", []))
        imdist_file_idx = list(si.get("imdistFileIdx", []))

    # placement: the overlaps after the shifts, the half-overlap regions, the canvas, batches
    pl = G.place_tiles(
        lay,
        xyz_shift,
        blend_method=blend_method,
        overlap_type=overlap_type,
        half_order=tuple(half_order),
        xcorr_shift=xcorr_shift,
        stitch_pad_size=tuple(stitch_pad_size),
        stitch_mip=mip_flags,
        out_bbox=None if _empty(out_bbox) else list(out_bbox),
        p_im_sz=p_im_sz_pl,
        batch_size=tuple(batch_size),
        block_size=tuple(block_size),
    )
    nv_size = [int(v) for v in pl.nv_size]
    bsz = [int(v) for v in pl.batch_size]
    blk = [int(v) for v in pl.block_size]
    num_batches = int(pl.num_batches)

    stitch_info_path = f"{data_path}/{result_dir}/{stitch_info_dir}"
    if not os.path.isdir(stitch_info_path):
        _mkdir_shared(stitch_info_path)

    from petakit5d.stitch import block as SB
    from petakit5d.stitch import blockinfo as SBI

    # a small volume is stitched in this process, if its blocks and distance maps fit it
    if SB.in_process(num_batches, bsz, dtype, blend_method, im_sizes, parse_cluster, config_file):
        parse_cluster = False
    task_size = blend_task_size(num_batches, bsz, parse_cluster, local_worker_num)
    # secondary channels use the per-task block info of the primary channel, so they take its
    # task size, which depends on the mode and the local workers
    if not is_primary_ch and _mat_has(stitch_info_fullpath, "PerBlockInfoFullpaths"):
        import re

        pbi = read_stitch_info(stitch_info_fullpath, ["PerBlockInfoFullpaths"])
        tok = re.search(r"_task_size_(\d+)", pbi["PerBlockInfoFullpaths"][0])
        if tok:
            task_size = int(tok.group(1))

    block_info_fullname, per_block_info_fullpaths, block_info_bytes = SBI.process_block_info(
        pl.int_xyz_shift,
        im_sizes,
        nv_size,
        bsz,
        data_order,
        pl.overlap_matrix,
        pl.full_regions,
        pl.half_regions,
        pl.pairs,
        tuple(border_size),
        zarr_fullpaths,
        stitch_info_path,
        nv_fsname,
        is_primary_ch or own_block_info,
        stitch_info_fullpath=stitch_info_fullpath,
        stitch_2d=stitch_2d,
        task_size=task_size,
        max_file_num_per_folder=int(max_file_num_per_folder),
        uuid=uuid,
        parse_cluster=parse_cluster,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )

    # the temporary stitched zarr next to the final one, named after it: it is published by a
    # rename, and the large-file MIPs are named after the final one
    if os.name == "nt" and len(uuid) > 4:
        nv_tmp_fullname = f"{nv_fullname}_{uuid[:4]}"
    else:
        nv_tmp_fullname = f"{nv_fullname}_{uuid}"
    nv_tmp_raw_fullname = nv_tmp_fullname
    fresh_stitch = True
    if _isdir(nv_tmp_raw_fullname):
        from petakit5d.io import get_zarr_info

        try:
            info = get_zarr_info(nv_tmp_raw_fullname)
            same = (
                len(info.size) == len(nv_size)
                and list(info.size) == nv_size
                and list(info.block_size) == blk
                # (the port writes only into its own layout)
                and info.axes == "zyx"
            )
        except (OSError, ValueError, KeyError):
            same = False
        if same:
            fresh_stitch = False
        else:
            shutil.rmtree(nv_tmp_raw_fullname)
    # an existing temporary result is kept: a master that ran out of time leaves it for the next
    if fresh_stitch:
        from petakit5d.io import create_zarr

        sep = (
            "/"
            if math.prod(math.ceil(a / b) for a, b in zip(nv_size, blk, strict=True)) > 10000
            else "."
        )
        create_zarr(
            nv_tmp_raw_fullname,
            data_size=nv_size,
            block_size=blk,
            shard_size=None if _empty(shard_size) else _ints(shard_size),
            dtype=dtype,
            compressor=compressor,
            dim_separator=sep,
        )

    # feather blending: the tiles' distance maps
    stitch_path = f"{data_path}/{result_dir}/"
    if blend_method.lower() == "feather":
        # one map per location only if the tiles of each location have one size
        uniq: dict[float, int] = {}
        for k, v in enumerate(loc_ids):
            uniq.setdefault(v, k)
        for v in sorted(uniq):
            same = np.all(im_sizes[np.asarray(loc_ids) == v] == im_sizes[uniq[v]])
            single_dist_map = bool(single_dist_map and same)
            if not single_dist_map:
                break
        dist_kw = {
            "data_order_axes": dop,
            "blend_weight_degree": blend_weight_degree,
            "single_dist_map": single_dist_map,
            "loc_ids": loc_ids,
            "dist_bboxes": None if _empty(dist_bboxes) else dist_bboxes,
            "block_size": _half(blk),
            "shard_size": _half(shard_size),
            "compressor": compressor,
            "large_file": large_file,
            "pool_size": pool_size,
            "parse_cluster": parse_cluster,
            "mcc_mode": mcc_mode,
            "config_file": config_file,
        }
        from petakit5d.stitch.distance import compute_tile_distance_transform

        use_primary = False
        if not is_primary_ch:
            try:
                np.broadcast_shapes(im_sizes.shape, p_tile_sizes_si.shape)
            except ValueError:
                raise StitchVolumeError(
                    "Arrays have incompatible sizes for this operation.",
                    "MATLAB:sizeDimensionsMustMatch",
                ) from None
            use_primary = bool(np.all(im_sizes == p_tile_sizes_si))
        if is_primary_ch or not use_primary:
            _mkdir_shared(f"{data_path}/{result_dir}/imdist/")
            imdist_fullpaths, imdist_file_idx = compute_tile_distance_transform(
                block_info_fullname, stitch_path, zarr_fullpaths, **dist_kw
            )
    else:
        imdist_fullpaths, imdist_file_idx = [], []

    # the stitch information of a primary channel
    if is_primary_ch:
        stitch_info_tmp = f"{result_root}/{stitch_info_dir}/{nv_fsname}_{uuid}.mat"
        stitch_info_fullpath = f"{result_root}/{stitch_info_dir}/{nv_fsname}.mat"
        p_im_sz_si = nv_size
        write_stitch_info(
            stitch_info_tmp,
            params=params,
            flipped_tile=flipped_tile,
            overlap_regions=lay.overlap_regions,
            overlap_matrix=pl.overlap_matrix,
            dxyz_shift=dxyz_shift,
            shift_method=shift_method,
            p_im_sz=p_im_sz_si,
            p_tile_sizes=im_sizes,
            xyz_voxel_sizes=vs,
            data_order=data_order,
            per_block_info_fullpaths=per_block_info_fullpaths,
            imdist_fullpaths=imdist_fullpaths,
            imdist_file_idx=imdist_file_idx,
            xyz_orig=lay.xyz_orig,
        )
        os.replace(stitch_info_tmp, stitch_info_fullpath)
    else:
        p_im_sz_si = [int(v) for v in p_im_sz_pl]

    # 'separate': each tile on a canvas of its own at its place in the stitched volume
    if blend_method.lower() == "separate":
        from petakit5d.jobs.publish import publish_result

        SB.process_separate_tiles(
            zarr_fullpaths,
            stitch_path,
            pl.int_xyz_shift,
            p_im_sz_si,
            batch_size=[max(2048, b) for b in blk],
            block_size=blk,
            parse_cluster=parse_cluster,
            master_compute=master_compute,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        nv_fullname = f"{data_path}/{result_dir}/{nv_fsname}.zarr"
        publish_result(nv_tmp_fullname, nv_fullname)
        return

    _blend(
        nv_fullname=nv_fullname,
        nv_tmp_fullname=nv_tmp_fullname,
        result_root=result_root,
        nv_fsname=nv_fsname,
        uuid=uuid,
        fresh_stitch=fresh_stitch,
        num_batches=num_batches,
        task_size=task_size,
        max_file_num_per_folder=int(max_file_num_per_folder),
        block_info_fullname=block_info_fullname,
        per_block_info_fullpaths=per_block_info_fullpaths,
        block_info_bytes=block_info_bytes,
        nv_size=nv_size,
        batch_size=bsz,
        dtype=dtype,
        blend_method=blend_method,
        border_size=[int(v) for v in border_size],
        imdist_fullpaths=imdist_fullpaths,
        imdist_file_idx=imdist_file_idx,
        pool_size=pool_size,
        blend_weight_degree=blend_weight_degree,
        save_mip=save_mip,
        cluster_mode=cluster_mode,
        parse_cluster=parse_cluster,
        master_compute=master_compute,
        node_factor=node_factor,
        job_log_dir=job_log_dir,
        mcc_mode=mcc_mode,
        config_file=config_file,
        local_worker_num=local_worker_num,
        mip_flags=mip_flags,
        debug=debug,
    )


def _blend(**k: Any) -> None:
    """the blending tasks, the MIP and the publication of the stitched volume"""
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.text import num2str, through_text
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.files import batch_file_exist
    from petakit5d.jobs.publish import publish_result
    from petakit5d.jobs.tasks import Task
    from petakit5d.stitch import block as SB

    _print("Process blocks for stitched result...")
    nv_fullname, nv_tmp = k["nv_fullname"], k["nv_tmp_fullname"]
    num_batches, task_size = k["num_batches"], k["task_size"]
    n_file = k["max_file_num_per_folder"]
    num_tasks = math.ceil(num_batches / task_size)
    debug = bool(k["debug"])

    # the flag folder of this run; flags of an earlier run are removed unless it is resumed
    flag_root = f"{k['result_root']}/zarr_flags/"
    flag_path = f"{k['result_root']}/zarr_flags/{k['nv_fsname']}_{k['uuid']}/"
    if _isdir(flag_path) and k["fresh_stitch"]:
        shutil.rmtree(flag_path)
    os.makedirs(flag_path, exist_ok=True)
    group_writable(flag_root)
    group_writable(flag_path)
    flag_paths = [flag_path]
    if num_tasks > n_file:
        flag_paths = []
        for f in range(math.ceil(num_tasks / n_file)):
            sub = f"{flag_path}/tasks_{f * n_file + 1}_{(f + 1) * n_file}/"
            _mkdir_shared(sub)
            flag_paths.append(sub)

    dtype = k["dtype"]
    save_mip = bool(k["save_mip"])
    nv_size = k["nv_size"]
    mip_mem = math.nan
    fuse_mip = False
    # data up to half of the memory: the tasks save the z MIPs of their batches, merged into the
    # result's MIP, instead of reading the result again (one slice: the MIP is the result)
    if save_mip:
        mip_mem = SB.mip_memory(k["cluster_mode"], k["config_file"])
        fuse_mip = (
            len(nv_size) == 3
            and nv_size[0] > 1
            and math.prod(nv_size) * data_type_to_byte_number(dtype) / 2**30 < mip_mem / 2
        )

    imdist = list(k["imdist_fullpaths"])
    idx = list(k["imdist_file_idx"])
    pool = None if _empty(k["pool_size"]) else _ints(k["pool_size"])
    # the weight degree as MATLAB's task text has it (%d)
    wd = float(through_text(float(k["blend_weight_degree"]), "%d"))
    flags, tasks = [], []
    for t in range(num_tasks):
        b0, b1 = t * task_size, min((t + 1) * task_size, num_batches)
        folder = flag_paths[math.ceil((t + 1) / n_file) - 1] if num_tasks > n_file else flag_path
        flags.append(f"{folder}/blocks_{b0 + 1}_{b1}.mat")
        tasks.append(
            Task.of(
                SB.process_stitch_block,
                list(range(b0, b1)),
                k["block_info_fullname"],
                k["per_block_info_fullpaths"][t],
                flags[-1],
                nv_tmp,
                None,
                None,
                im_size=list(nv_size),
                batch_size=list(k["batch_size"]),
                dtype=dtype,
                blend_method=k["blend_method"],
                border_size=list(k["border_size"]),
                imdist_fullpaths=imdist,
                imdist_file_idx=idx,
                pool_size=pool,
                weight_degree=wd,
                save_mip=fuse_mip,
            )
        )
    inputs = list(k["per_block_info_fullpaths"])
    outputs = list(flags)
    if k["parse_cluster"]:
        # the tasks with the most block information first
        order = _task_order(k["block_info_bytes"])
        inputs = [inputs[i] for i in order]
        outputs = [outputs[i] for i in order]
        tasks = [tasks[i] for i in order]

    byte_num = data_type_to_byte_number(dtype)
    cpus = 1 * k["node_factor"]
    mem = math.prod(k["batch_size"]) * byte_num / 1024**3 * 20
    jtl = task_size * (3 / 60)
    common = {
        "final_out_fullpath": nv_fullname,
        "max_trial_num": 2,
        "master_compute": k["master_compute"],
        "parse_cluster": k["parse_cluster"],
        "job_log_dir": k["job_log_dir"],
        "mcc_mode": k["mcc_mode"],
        "config_file": k["config_file"],
    }
    done = None
    if not _isdir(nv_fullname):
        done = framework.run_tasks(
            inputs, outputs, tasks, cpus_per_task=cpus, mem_allocate=mem, job_time_limit=jtl,
            local_worker_num=k["local_worker_num"], **common,
        )  # fmt: skip
    # retries with more CPUs (memory) and time
    for i in range(1, 4):
        if not _isdir(nv_fullname) and not all(done or [False]):
            done = framework.run_tasks(
                inputs, outputs, tasks, cpus_per_task=cpus * 2**i, job_time_limit=jtl * 2**i,
                **common,
            )  # fmt: skip
    if not _isdir(nv_fullname) and done is not None and not all(done):
        if not debug:
            _rmtree(nv_tmp)
            _rmtree(flag_path)
        if _isdir(nv_fullname):
            _print(f"Stitched result {nv_fullname} is finished by another job, skip it!")
            return
        failed = [i + 1 for i, d in enumerate(done) if not d]
        raise StitchVolumeError(
            f"The block processing for blocks [{num2str(failed)}] cannot be finished!"
        )

    # another run may have published the result while this one ran (the job framework then
    # reports all tasks as done): that result and its MIP are kept
    if _isdir(nv_fullname) and (
        k["mip_flags"] is None or not all(batch_file_exist(flags, None, True))
    ):
        if not debug:
            _rmtree(nv_tmp)
            _rmtree(flag_path)
        _print(f"Stitched result {nv_fullname} is finished by another job, skip it!")
        return

    # the z MIP from the blending tasks' flags (None if a flag has none)
    fused = None
    if fuse_mip and _isdir(nv_tmp):
        fused = SB.merge_block_mips(flags, nv_size, dtype)
    if not debug and _isdir(flag_path):
        shutil.rmtree(flag_path)

    # the MIP from the temporary result: the caller counts the stitching as done once the
    # result exists
    if save_mip and _isdir(nv_tmp):
        t0 = time.perf_counter()
        mip_path = f"{k['result_root']}/MIPs/"
        if not os.path.isdir(mip_path):
            _mkdir_shared(mip_path)
        mip_name = f"{mip_path}{k['nv_fsname']}_MIP_z.tif"
        if fused is not None:
            from petakit5d.io import write_tiff

            tmp = f"{mip_name[:-4]}_{get_uuid()}.tif"
            write_tiff(fused, tmp)
            os.replace(tmp, mip_name)
        elif math.prod(nv_size) * byte_num / 2**30 < mip_mem / 2:
            from petakit5d.tools.mip import save_mip_zarr

            save_mip_zarr(nv_tmp, mip_name, dtype=dtype, axis=(1, 0, 0))
        else:
            # data greater than half of the memory: the large-file MIPs (cluster)
            from petakit5d.tools.mip import mip_zarr

            mip_zarr(
                nv_tmp,
                axis=(1, 1, 1),
                parse_cluster=k["parse_cluster"],
                master_compute=k["master_compute"],
                job_log_dir=k["job_log_dir"],
                mcc_mode=k["mcc_mode"],
                config_file=k["config_file"],
            )
        _print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.")

    # publish (a result another writer published meanwhile is kept, also for MIP stitching,
    # where MATLAB replaced it)
    if _isdir(nv_tmp):
        publish_result(nv_tmp, nv_fullname)
    if _isdir(nv_tmp) and not debug:
        shutil.rmtree(nv_tmp)
