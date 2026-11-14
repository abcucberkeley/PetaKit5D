"""Stitching of datasets (``XR_matlab_stitching_wrapper``): every volume of an image list,
stitched by ``stitch.volume.stitch_volume`` (XR_stitching_frame_zarr_dev_v1) here or in Slurm
jobs.

``stitch_dataset`` is XR_matlab_stitching_wrapper.m at the oracle commit:

- **Folders and parameters**: ``<data path>/<result_dir_name>/`` (of the first data path with
  ``multi_loc``) with ``stitchInfo/`` and ``tmp/``, made group-writable, and parameters.mat
  (v7.3) and parameters.json of the call (MATLAB's ``ip.Results``: the values as given).
- **Volumes**: the image list (S01's ``imagelist.parse_image_list``, or
  ``parse_multi_loc_image_list`` for one list per folder) gives the tiles; a volume is a time
  point (fullIter), camera, stack, channel and, for online stitching while streaming, z layer.
  Its result is ``<result folder>/<prefix>Scan_Iter_<iter>[_Cam<c>]_ch<n>_CAM1_stack<s>_<laser>
  nm_<msec>msec_<msecAbs>msecAbs[_<z>z][_<n>ntile].zarr``, named after its first tile; a
  result that exists is kept. Tiles may be processed ones (``processed_dir_str``, e.g. DSR
  tiles) of the raw names.
- **Primaries and secondaries** (``xcorr_mode``): 'primaryFirst' registers the primary
  channel (``primary_ch``, else the first found) at the first time point, and every other volume
  waits for its stitch information file and is placed with its shifts, so a time series is
  registered once and every time point has the same canvas; volumes with another number of
  tiles than the primary are skipped. 'primary' registers the primary channel at each time
  point, the other channels use it. 'all' registers every volume. A given stitch information
  file (``stitch_info_fullpath``) pins every volume to it ('stitchInfo').
- **Streaming** (``streaming``): tiles that don't exist yet are waited for (5 s rounds), up to
  100 rounds per tile of the list without a new tile (5 with ``online_stitch``); online
  stitching stitches the tiles there are, per z layer, names the result by its tile count and
  removes results of fewer tiles. Missing processed tiles count as tries, except in cluster
  mode while the pipeline's marker ``deskew_in_progress.tmp`` in their folder is less than 10
  minutes old (``processed_tiles_coming``).
- **Values in transit**: MATLAB hands each volume to the frame as a function string that it
  evaluates, here and in its jobs, so every number goes through ``sprintf`` (``%0.10f`` for the
  pixel size, z step, angle, offsets, threshold and edge erosion, ``%d`` for the tile offset,
  ``mat2str`` with 15 digits for the coordinates, sizes, boxes and weights, ``%.10d`` for the
  resample factor). ``frame_task_string`` writes that string and ``frame_task`` reads it back
  as MATLAB's parser does (``params.literal``) and binds it to ``stitch_volume``'s parameters,
  so the frame gets MATLAB's values in both modes. Without the cluster and with
  ``local_worker_num`` > 1 the string ends in ``'localWorkerNum',<n>``; in cluster mode without
  ``master_compute`` in ``'masterCompute',false``.
- **Local mode** (no Slurm, or ``parse_cluster`` False): each volume in turn, here, with
  MATLAB's ``.tmp`` markers (a marker younger than ``unit_wait_time`` minutes per tile is
  another process's volume), up to ``max_trial_num`` tries. ``local_worker_num`` None: 1, or
  for datasets with 32 or more tiles per volume to register or 32 G voxels in all, up to 8 local
  workers (one per 3 CPUs of this process) for the frame's tasks.
- **Cluster mode**: one single job per volume (``jobs.single``: never waiting, task runners
  from the second volume on, the call's pools ended with it) with ``coordinator_nice``, from
  ``config_file``. A volume of 30 batches or more (unless 2D or MIP with fewer than 50 tiles)
  runs its tile, map and block tasks as jobs it waits for, so its job is sized as a
  coordinator (``coordinator_resources``: one CPU and its MIP's memory, or with the config's
  ``masterCompute`` a block worker's 2 CPUs, 1 for big data); a smaller one is stitched in its
  job (``cpus_per_task``, 8 bytes per voxel of its tiles). With ``master_compute`` this process
  stitches the last unfinished volume itself if no job of it runs and it fits this process's
  job (``master_task_fits``), and the registered primary of 'primaryFirst' first. More than 10
  tiles go to the frame in ``stitchInfo/<name>_tile_info.mat``. A volume is done when its
  result exists (``batch_file_exist``); a volume's last try is waited for, and between rounds
  the wait returns within a second of the results, tiles or stitch information still to come.

The MATLAB function's quirks are kept: the axis-order check never raises; the job sizes use
the first tile size seen for every volume; with ``multi_loc`` a volume's number in the job's
array index, the last-volume rule and the job lookup is its tile count (the tile loops reuse
MATLAB's variable ``f``); 'primary' secondaries without camera names take the first stack's
stitch information; a missing primary volume in 'primary' mode skips the first time point's
volumes of every camera and channel (and its secondary then fails, as MATLAB's index error);
tiles of one volume whose msec fields differ give a garbled name (MATLAB's sprintf of the
unique values), whose ``.tmp`` marker can't be made; ``process_fun_path`` entries are chosen by
the first tile's path (``contains``) of each volume. Where MATLAB's function string would be
garbled (a vector where it formats a scalar, a 3 x 2 ``out_bbox``), the port raises
``StitchDatasetError``.

Waiting on others: the registration (S04, through the frame: without it a registered primary
with overlapping tiles needs its ``xcorr/<name>.mat``), MATLAB-code tile functions (S02).
"""

from __future__ import annotations

import math
import os
import re
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "FRAME_FORMAT",
    "StitchDatasetError",
    "frame_task",
    "frame_task_string",
    "processed_tiles_coming",
    "split_call",
    "stitch_dataset",
]

# the jobs' function (XR_stitching_frame_zarr_dev_v1)
FRAME_FUNC = "petakit5d.stitch.volume:stitch_volume"
FRAME_NAME = "XR_stitching_frame_zarr_dev_v1"
ARRAY_ID_MODULUS = 5000
# tiles above which cluster jobs get them in a tile information file
TILE_INFO_MIN = 10
# the pipeline's marker of deskew jobs still making tiles, and its age limit (minutes)
DESKEW_MARKER = "deskew_in_progress.tmp"
DESKEW_MARKER_MINUTES = 10
# a round's wait (s) for the results, tiles and stitch information still to come
ROUND_WAIT = 5
# a stitching job is a coordinator from this many batches on
COORDINATOR_BATCHES = 30
# 2D and MIP volumes of fewer tiles are stitched in their job
SMALL_2D_TILES = 50
BIG_DATA_BYTES = 100 * 2**30
# local workers without localWorkerNum: tiles per volume to register, voxels in all
LOCAL_WORKER_TILES = 32
LOCAL_WORKER_VOXELS = 32e9

# XR_matlab_stitching_wrapper's function string of a volume (lines 669-693), MATLAB's quotes
# undone
FRAME_FORMAT = (
    "XR_stitching_frame_zarr_dev_v1(%s,%s,'axisOrder','%s','dataOrder','%s',"
    "'xyPixelSize',%0.10f,'dz',%0.10f,'skewAngle',%0.10f,'reverse',%s,'objectiveScan',%s,"
    "'IOScan',%s,'resultPath','%s','tileInfoFullpath','%s','stitchInfoDir','%s',"
    "'stitchInfoFullpath','%s','processedDirStr','%s','DS',%s,'DSR',%s,'zarrFile',%s,"
    "'largeFile',%s,'poolSize',%s,'batchSize',%s,'blockSize',%s,'shardSize',%s,"
    "'resampleType','%s',"
    "'resampleFactor',[%s],'inputBbox',%s,'tileOutBbox',%s,'tileOffset',%d,'blendMethod','%s',"
    "'overlapType','%s','xcorrShift',%s,'xyMaxOffset',%0.10f,'zMaxOffset',%0.10f,"
    "'xcorrDownsample',%s,"
    "'xcorrThresh',%0.10f,'shiftMethod','%s','axisWeight',[%s],'groupFile','%s',"
    "'isPrimaryCh',%s,"
    "'usePrimaryCoords',%s,'outBbox',[%s],'save16bit',%s,"
    "'tileIdx',%s,'flippedTile',[%s],'processFunPath',%s,'stitchMIP',%s,'bigStitchData',%s,"
    "'edgeArtifacts',%0.10f,'distBboxes',%s,'saveMIP',%s,'parseCluster',%s,'uuid','%s',"
    "'mccMode',%s,'configFile','%s')"
)
# the values of the conversions, in order (s: text, l: string() of a logical, n: a number)
FRAME_FIELDS = (
    ("tile_fullpaths", "s"),
    ("xyz", "s"),
    ("axis_order", "s"),
    ("data_order", "s"),
    ("xy_pixel_size", "n"),
    ("dz", "n"),
    ("skew_angle", "n"),
    ("reverse", "l"),
    ("objective_scan", "l"),
    ("io_scan", "l"),
    ("result_path", "s"),
    ("tile_info_fullpath", "s"),
    ("stitch_info_dir", "s"),
    ("stitch_info_fullpath", "s"),
    ("processed_dir_str", "s"),
    ("ds", "l"),
    ("dsr", "l"),
    ("zarr_file", "l"),
    ("large_file", "l"),
    ("pool_size", "s"),
    ("batch_size", "s"),
    ("block_size", "s"),
    ("shard_size", "s"),
    ("resample_type", "s"),
    ("resample_factor", "s"),
    ("input_bbox", "s"),
    ("tile_out_bbox", "s"),
    ("tile_offset", "n"),
    ("blend_method", "s"),
    ("overlap_type", "s"),
    ("xcorr_shift", "l"),
    ("xy_max_offset", "n"),
    ("z_max_offset", "n"),
    ("xcorr_downsample", "s"),
    ("xcorr_thresh", "n"),
    ("shift_method", "s"),
    ("axis_weight", "s"),
    ("group_file", "s"),
    ("is_primary_ch", "l"),
    ("use_primary_coords", "l"),
    ("out_bbox", "s"),
    ("save_16bit", "l"),
    ("tile_idx", "s"),
    ("flipped_tile", "s"),
    ("process_fun_path", "s"),
    ("stitch_mip", "s"),
    ("big_stitch_data", "l"),
    ("edge_artifacts", "n"),
    ("dist_bboxes", "s"),
    ("save_mip", "l"),
    ("parse_cluster", "l"),
    ("uuid", "s"),
    ("mcc_mode", "l"),
    ("config_file", "s"),
)
# the volume names (the result, its .tmp marker and the stitch information files)
_NAME_CAM = "%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s"
_NAME_NOCAM = "%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s"
_CONVERSION = re.compile(r"%(?:%|[-+ 0#]*\d*(?:\.\d+)?[dfsgeiuxXoc])")


class StitchDatasetError(ValueError):
    """An error where MATLAB raises one, with its identifier ('' where MATLAB's has none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


def _print(msg: str = "") -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------------------------------------
# the function string of a volume and its values
# ---------------------------------------------------------------------------------------------


def _pieces(fmt: str) -> list[str]:
    """the format cut after each conversion's literal text: one conversion per piece"""
    out, start = [], 0
    convs = [m for m in _CONVERSION.finditer(fmt) if m.group(0) != "%%"]
    for k, m in enumerate(convs):
        end = convs[k + 1].start() if k + 1 < len(convs) else len(fmt)
        out.append(fmt[start:end])
        start = end
        del m
    return out


_FRAME_PIECES = _pieces(FRAME_FORMAT)


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    return np.size(np.asarray(v, dtype=object)) == 0


def _matlab_string(v: Any) -> str:
    """``string(x)`` of a logical (or of [], which prints nothing)"""
    if _empty(v):
        return ""
    a = np.asarray(v)
    if a.size != 1:
        raise StitchDatasetError(
            f"A logical value {v!r} with several elements garbles MATLAB's function string"
        )
    return "true" if bool(a.reshape(-1)[0]) else "false"


def _scalar_number(name: str, v: Any) -> float | int | bool:
    a = np.asarray(v)
    if a.size != 1 or a.dtype.kind not in "biuf":
        raise StitchDatasetError(
            f"{name}: {v!r} is not one number; in MATLAB's function string it would shift the "
            "values after it (MATLAB's volumes would get other parameters)"
        )
    x = a.reshape(-1)[0]
    return bool(x) if a.dtype.kind == "b" else x.item()


def frame_task_string(values: dict[str, Any]) -> str:
    """XR_matlab_stitching_wrapper's function string of a volume, ``sprintf`` of
    ``FRAME_FORMAT``: ``values`` has the keys of ``FRAME_FIELDS`` (texts already in MATLAB's
    form: the tiles' cell text, ``mat2str`` of the matrices and sizes, ``num2str`` of the resample
    factor, the out box and the flips). Each conversion is formatted with its own value, as
    MATLAB does for scalars (an empty value prints nothing and keeps its place)."""
    from petakit5d.core.text import sprintf

    parts = []
    for piece, (key, kind) in zip(_FRAME_PIECES, FRAME_FIELDS, strict=True):
        v = values[key]
        if kind == "l":
            arg: Any = _matlab_string(v)
        elif kind == "s":
            arg = os.fspath(v) if isinstance(v, os.PathLike) else str(v)
        else:
            arg = _scalar_number(key, v)
        parts.append(sprintf(piece, arg))
    return "".join(parts)


def split_call(func_str: str) -> tuple[str, list[Any]]:
    """The function name and the arguments (MATLAB model values) of a function string, read as
    MATLAB's parser reads it (``eval``)."""
    from petakit5d.params import literal
    from petakit5d.params.values import MArray

    m = re.match(r"\s*([A-Za-z]\w*)\((.*)\)\s*$", func_str, re.S)
    if m is None:
        raise StitchDatasetError(f"Not a function call: {func_str}")
    args = literal.parse("{" + m.group(2) + "}")
    # (char rows as text: parameter names, paths)
    out = [
        a.text()
        if isinstance(a, MArray) and a.mclass == "char" and len(a.shape) == 2 and a.shape[0] <= 1
        else a
        for a in args
    ]
    return m.group(1), out


def frame_task(func_str: str) -> Any:
    """The job task of a function string of XR_stitching_frame_zarr_dev_v1: ``stitch_volume``
    with the values MATLAB's evaluation of the string gives (names, case and conversions
    through the registry)."""
    from petakit5d.jobs.tasks import Task
    from petakit5d.params.registry import REGISTRY

    name, args = split_call(func_str)
    spec = REGISTRY.find(name)
    kwargs = spec.bind_matlab(args)
    hook = spec.matlab_in
    if hook is not None:
        kwargs = hook(kwargs)
    return Task(FRAME_FUNC, kwargs=kwargs)


def _mat2str_commas(v: Any) -> str:
    """``strrep(mat2str(v), ' ', ',')`` of a MATLAB value (None is [])"""
    from petakit5d.core.text import mat2str

    if v is None:
        return "[]"
    a = np.asarray(v)
    if a.dtype.kind not in "biuf":
        raise StitchDatasetError(f"{v!r} is not numeric or logical")
    if a.dtype.kind == "f" and a.dtype != np.float64:
        a = a.astype(np.float64)
    return mat2str(a).replace(" ", ",")


def _num2str_commas(v: Any, fmt: str) -> str:
    """``strrep(num2str(v, fmt), ' ', '')`` of a MATLAB row (None is [], which prints '')"""
    from petakit5d.core.text import num2str

    if v is None:
        return ""
    a = np.asarray(v)
    if a.size == 0:
        return ""
    if a.ndim > 1 and a.shape[0] > 1:
        raise StitchDatasetError(
            f"{v!r}: num2str of a matrix gives several rows, which garble MATLAB's function string"
        )
    return num2str(a.reshape(1, -1), fmt).replace(" ", "")


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _as_list(x: Any) -> list[Any]:
    if x is None:
        return []
    if isinstance(x, (str, os.PathLike)):
        return [os.fspath(x)]
    return [os.fspath(e) if isinstance(e, os.PathLike) else e for e in x]


def _mtrue(v: Any) -> bool:
    """MATLAB's ``if v``: not empty and all nonzero"""
    a = np.asarray(v)
    return bool(a.size > 0 and np.all(a))


def _unique(a: np.ndarray) -> np.ndarray:
    """MATLAB's unique of a double column (sorted; NaNs kept, each one)"""
    a = np.asarray(a, dtype=np.float64).ravel()
    fin = np.unique(a[~np.isnan(a)])
    return np.concatenate([fin, a[np.isnan(a)]])


def _mkdir_shared(path: str) -> None:
    from petakit5d.core.paths import mkdir_recursive
    from petakit5d.io.tiff import group_writable

    if not os.path.isdir(path):
        mkdir_recursive(path)
        group_writable(path)


def _save_parameters(root: str, pr: dict[str, Any]) -> None:
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.io.tiff import group_writable

    write_mat(root + "/parameters.mat", {"pr": pr})
    write_json_file(pr, root + "/parameters.json")
    group_writable(root + "/parameters.mat")
    group_writable(root + "/parameters.json")


def _touch(path: str) -> None:
    """``fclose(fopen(path, 'w'))``"""
    try:
        with open(path, "w"):
            pass
    except (OSError, ValueError):
        # (fopen gives -1, which fclose refuses)
        raise StitchDatasetError(
            "Invalid file identifier. Use fopen to generate a valid file identifier.",
            "MATLAB:FileIO:InvalidFid",
        ) from None


def _delete(path: str) -> None:
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def _exist_file(path: str) -> bool:
    """``exist(path, 'file')`` (a file or folder)"""
    from petakit5d.jobs.files import exist

    return exist(path)


def _batch_exist(paths: Sequence[str]) -> list[bool]:
    from petakit5d.jobs import files

    return files.batch_file_exist(list(paths))


def _wait(groups: list[list[str]], wait_time: float) -> None:
    from petakit5d.jobs import files

    files.wait_for_files(groups, wait_time)


def processed_tiles_coming(tile_fullpaths: Sequence[str]) -> bool:
    """Whether jobs still make these processed tiles: the pipeline (XR_microscopeAutomatic
    Processing) keeps the marker ``deskew_in_progress.tmp`` in their folder, written again each
    minute, while deskew jobs of the dataset are pending or running. A marker 10 minutes old is
    left from a pipeline that ended."""
    from petakit5d.core.paths import fileparts

    folders = sorted({fileparts(p)[0] for p in tile_fullpaths})
    for folder in folders:
        marker = f"{folder}/{DESKEW_MARKER}"
        # (batch_file_exist sees a marker another node just made)
        if _batch_exist([marker])[0]:
            try:
                age = (time.time() - os.path.getmtime(marker)) / 60
            except OSError:
                continue
            if age < DESKEW_MARKER_MINUTES:
                return True
    return False


def _get_image_size(path: str) -> tuple[int, ...]:
    from petakit5d.io import get_image_size

    return tuple(int(v) for v in get_image_size(path))


def _sprintf(fmt: str, *args: Any) -> str:
    from petakit5d.core.text import sprintf

    return sprintf(fmt, *args)


# ---------------------------------------------------------------------------------------------
# the wrapper
# ---------------------------------------------------------------------------------------------


@register(
    "XR_matlab_stitching_wrapper",
    Param(
        "dataPaths",
        "@(x) ischar(x) || iscell(x)",
        "Folder(s) of the tiles (one per group with multiLoc).",
    ),
    Param(
        "imageListFullpaths",
        "@(x) ischar(x) || iscell(x)",
        "Image list(s), in the order of the folders.",
    ),
    Param("resultDirName", "@ischar", "Result folder under the (first) data folder."),
    Param("streaming", "@islogical", "Wait for tiles that don't exist yet (real-time processing)."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the tiles to stitch."),
    Param("multiLoc", "@islogical", "Tiles in several folders, a group per folder."),
    Param(
        "processedDirStr",
        "@ischar",
        "Folder of processed tiles under the data folder (e.g. 'DSR').",
    ),
    Param(
        "stitchInfoFullpath", "@ischar", "A stitch information file to stitch every volume with."
    ),
    Param("DS", "@islogical", "Tiles deskewed."),
    Param("DSR", "@islogical", "Tiles deskewed and rotated (stage coordinates)."),
    Param("xyPixelSize", "@isnumeric", "Pixel size (um)."),
    Param("dz", "@isnumeric", "Scan step (um)."),
    Param("skewAngle", "@isnumeric", "Skew angle (degrees)."),
    Param("reverse", "@islogical", "Reverse the z axis."),
    Param(
        "parseSettingFile",
        "@islogical",
        "Flip tiles whose setting files have a negative stage interval.",
    ),
    Param(
        "axisOrder",
        "@ischar",
        "Axis order of the stage coordinates ('x,y,z', '-x,y,z', ...).",
        convert="same",
    ),
    Param("dataOrder", "@ischar", "Axes of the tile arrays ('z,y,x').", convert="data_order"),
    Param("objectiveScan", "@islogical", "Objective scan (DS space)."),
    Param("IOScan", "@islogical", "Inverted objective scan (stage coordinates, DSR space)."),
    Param("zarrFile", "@islogical", "Tiles are zarr."),
    Param("largeFile", "@islogical", "Large tiles (registration and maps on MIP slabs)."),
    Param(
        "poolSize",
        "@isnumeric",
        "Pooling of the MIP slabs (z, y, x) for large files.",
        convert="size_yxz",
    ),
    Param(
        "batchSize", "@isnumeric", "Batch size (z, y, x) of the blending tasks.", convert="size_yxz"
    ),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the results.", convert="size_yxz"),
    Param(
        "shardSize",
        "@isnumeric",
        "Shard size (z, y, x) of the results (Zarr v3).",
        convert="size_yxz",
    ),
    Param("resampleType", "@ischar", "'given', 'isotropic' or 'xy_isotropic'.", convert="same"),
    Param(
        "resampleFactor",
        "@isnumeric",
        "Resampling (z, y, x) in the tiles' axis order (None: none).",
        convert="resample_yxz",
    ),
    Param(
        "inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) of each tile read.", convert="bbox"
    ),
    Param(
        "tileOutBbox", "@isnumeric", "Box of each tile kept after its processing.", convert="bbox"
    ),
    Param(
        "tileOffset", "@isnumeric", "Counts added to the tiles in their processing.", convert="same"
    ),
    Param("blendMethod", "@ischar", "'none', 'feather', 'mean', 'median', 'max' or 'separate'."),
    Param("overlapType", "@ischar", "'' (by the blend method), 'none', 'half' or 'full'."),
    Param("xcorrShift", "@islogical", "Register the tiles."),
    Param(
        "xyMaxOffset",
        "@(x) isvector(x) && numel(x) <= 2",
        "Largest shift in x and y (voxels).",
        convert="same",
    ),
    Param("zMaxOffset", "@isnumeric", "Largest shift in z (voxels)."),
    Param(
        "xcorrDownsample",
        "@isnumeric",
        "Downsampling (z, y, x) of the registration.",
        convert="size_yxz",
    ),
    Param("xcorrThresh", "@isnumeric", "Pairs below this correlation are not used."),
    Param(
        "outBbox",
        "@(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6)",
        "Box of the stitched volume kept (None ends: to the end).",
        convert="bbox",
    ),
    Param(
        "xcorrMode",
        "@(x) ismember(x, {'primary', 'primaryFirst', 'all'})",
        "'primaryFirst', 'primary' or 'all' (see the module).",
    ),
    Param("shiftMethod", "@ischar", "'grid', 'local', 'global' or 'group'."),
    Param(
        "axisWeight",
        "@isnumeric",
        "Weights of the tile index axes (z, y, x) in the assignment.",
        convert="size_yxz",
    ),
    Param("groupFile", "@ischar", "Groups of tiles (csv) for 'group'."),
    Param(
        "primaryCh",
        "@(x) isempty(x) || ischar(x)",
        "The primary channel ('Cam<c>_ch<n>' or 'ch<n>'; '': the first found).",
    ),
    Param(
        "usePrimaryCoords",
        "@islogical",
        "Secondaries take the primary's stage coordinates.",
        convert="same",
    ),
    Param("save16bit", "@islogical", "Save the results as uint16.", convert="same"),
    Param("edgeArtifacts", "@isnumeric", "Edge erosion of the tiles (voxels; 0: none)."),
    Param("distBboxes", "@isnumeric", "Boxes of the feather distance maps.", convert="bboxes"),
    Param("saveMIP", "@islogical", "Save the z MIPs of the results."),
    Param(
        "stitchMIP",
        "@(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3))",
        "Stitch the tiles' MIPs (flags z, y, x; one value: z).",
        convert="flags_yxz",
    ),
    Param(
        "onlineStitch",
        "@(x) islogical(x)",
        "Stitch the tiles there are (with streaming, per z layer).",
    ),
    Param(
        "processFunPath",
        "@(x) isempty(x) || ischar(x) || iscell(x)",
        "Function files (.mat or .txt) for the tiles, one per channel pattern.",
    ),
    Param("parseCluster", "@islogical", "Stitch the volumes in Slurm jobs."),
    Param("masterCompute", "@islogical", "This process stitches volumes too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a volume's job (a volume stitched in its job)."),
    Param("uuid", "@ischar", "Id of the temporary names (empty: a new one)."),
    Param("maxTrialNum", "@isnumeric", "Tries per volume."),
    Param(
        "unitWaitTime", "@isnumeric", "Minutes per tile a local .tmp marker counts as another's."
    ),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; passed on, ignored)."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param(
        "localWorkerNum",
        "@isnumeric",
        "Local workers of the frame's tasks without the cluster (None: by the dataset's size).",
    ),
    task="S09",
)
def stitch_dataset(
    data_paths,
    image_list_fullpaths,
    *,
    result_dir_name="matlab_stitch",
    streaming=False,
    channel_patterns=("CamA_ch0", "CamA_ch1", "CamB_ch0"),
    multi_loc=False,
    processed_dir_str="",
    stitch_info_fullpath="",
    ds=False,
    dsr=False,
    xy_pixel_size=0.108,
    dz=0.5,
    skew_angle=32.45,
    reverse=False,
    parse_setting_file=False,
    axis_order="x,y,z",
    data_order="z,y,x",
    objective_scan=False,
    io_scan=False,
    zarr_file=False,
    large_file=False,
    pool_size=None,
    batch_size=(512, 512, 512),
    block_size=(256, 256, 256),
    shard_size=None,
    resample_type="xy_isotropic",
    resample_factor=None,
    input_bbox=None,
    tile_out_bbox=None,
    tile_offset=0,
    blend_method="feather",
    overlap_type="",
    xcorr_shift=True,
    xy_max_offset=300,
    z_max_offset=50,
    xcorr_downsample=(1, 2, 2),
    xcorr_thresh=0.25,
    out_bbox=None,
    xcorr_mode="primaryFirst",
    shift_method="grid",
    axis_weight=(10, 1, 0.1),
    group_file="",
    primary_ch="",
    use_primary_coords=False,
    save_16bit=True,
    edge_artifacts=0,
    dist_bboxes=None,
    save_mip=True,
    stitch_mip=None,
    online_stitch=False,
    process_fun_path="",
    parse_cluster=True,
    master_compute=True,
    job_log_dir="../job_logs",
    cpus_per_task=8,
    uuid="",
    max_trial_num=3,
    unit_wait_time=0.1,
    mcc_mode=False,
    config_file="",
    local_worker_num=None,
) -> None:
    """Stitch every volume of the image list(s) (XR_matlab_stitching_wrapper): see the module.
    Sizes and per-axis values in (z, y, x), boxes 0-based and end-exclusive."""
    values = dict(locals())
    from petakit5d.core.ids import get_uuid
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(stitch_dataset, values).items()}

    # the root folder
    if multi_loc:
        paths = [str(p) for p in _as_list(data_paths)]
        lists = [str(p) for p in _as_list(image_list_fullpaths)]
        stitching_rt = f"{paths[0]}/{result_dir_name}"
        data_path_1: Any = paths
    else:
        paths = _as_list(data_paths)
        lists = _as_list(image_list_fullpaths)
        data_path_1 = str(paths[0]) if paths else ""
        lists = [str(lists[0])] if lists else [""]
        stitching_rt = f"{data_path_1}/{result_dir_name}"
    _mkdir_shared(stitching_rt)
    _save_parameters(stitching_rt, pr)

    # the axis order (MATLAB's check never raises: `~regexpi(...)` is false for a match and
    # for no match)
    axis_order = str(axis_order).replace(" ", "")

    stitch_info_dir = "stitchInfo"
    stitch_info_path = f"{stitching_rt}/{stitch_info_dir}"
    _mkdir_shared(stitch_info_path)
    _mkdir_shared(f"{stitching_rt}/tmp")

    job_log_fname = job_error_fname = ""
    if parse_cluster:
        from petakit5d.jobs.slurm import check_slurm_cluster

        check = check_slurm_cluster(paths[0] if multi_loc else data_path_1, job_log_dir)
        parse_cluster = bool(check.parse_cluster)
        job_log_fname, job_error_fname = check.job_log_fname, check.job_log_error_fname

    if not uuid:
        uuid = get_uuid()
    mip_flags = _mtrue_any(stitch_mip)
    if mip_flags:
        cpus_per_task = 1

    use_processed = bool(processed_dir_str)
    from petakit5d.stitch import imagelist as IL

    kw = {
        "streaming": bool(streaming),
        "online_stitch": bool(online_stitch),
        "stitch_info_fullpath": str(stitch_info_fullpath),
        "stitch_info_path": stitch_info_path,
        "zarr_file": bool(zarr_file),
        "channel_patterns": [str(c) for c in _as_list(channel_patterns)],
        "use_processed_data": use_processed,
        "processed_dir_str": str(processed_dir_str),
        "xcorr_mode": str(xcorr_mode),
        "primary_ch": "" if primary_ch is None else str(primary_ch),
    }
    if multi_loc:
        info = IL.parse_multi_loc_image_list(paths, lists, **kw)
    else:
        info = IL.parse_image_list(data_path_1, lists[0], **kw)

    channels = [str(c) for c in _as_list(channel_patterns)]
    if isinstance(process_fun_path, str) or process_fun_path is None:
        fun_rows = [process_fun_path or ""] * len(channels)
    else:
        fun_rows = ["" if _empty(p) else str(p) for p in process_fun_path]

    s = {
        "multi_loc": bool(multi_loc),
        "paths": paths if multi_loc else [data_path_1],
        "stitching_rt": stitching_rt,
        "stitch_info_path": stitch_info_path,
        "stitch_info_dir": stitch_info_dir,
        "info": info,
        "use_processed": use_processed,
        "processed_dir_str": str(processed_dir_str),
        "streaming": bool(streaming),
        "online_stitch": bool(online_stitch),
        "parse_setting_file": bool(parse_setting_file),
        "channels": channels,
        "fun_rows": fun_rows,
        "parse_cluster": parse_cluster,
        "master_compute": bool(master_compute),
        "job_log_dir": job_log_dir,
        "job_log_fname": job_log_fname,
        "job_error_fname": job_error_fname,
        "cpus_per_task": cpus_per_task,
        "config_file": config_file,
        "mcc_mode": mcc_mode,
        "uuid": uuid,
        "unit_wait_time": float(unit_wait_time),
        "save_16bit": bool(save_16bit),
        "large_file": bool(large_file),
        "xcorr_shift": bool(xcorr_shift),
        "stitch_mip": mip_flags,
        "pr": pr,
        "axis_order": axis_order,
        # (the size estimate's Python values)
        "estimate": {
            "axis_order": axis_order,
            "data_order": data_order,
            "xy_pixel_size": xy_pixel_size,
            "dz": dz,
            "skew_angle": skew_angle,
            "objective_scan": objective_scan,
            "io_scan": io_scan,
            "ds": ds,
            "dsr": dsr,
            "resample_type": resample_type,
            "resample_factor": resample_factor,
            "input_bbox": input_bbox,
            "tile_out_bbox": tile_out_bbox,
            "xcorr_shift": xcorr_shift,
            "xy_max_offset": xy_max_offset,
            "z_max_offset": z_max_offset,
            "out_bbox": out_bbox,
            "stitch_mip": stitch_mip,
            "batch_size": batch_size,
            "block_size": block_size,
        },
    }
    s["local_worker_num"] = _local_worker_num(local_worker_num, s)
    loop = _Loop(s, max_trial_num=float(max_trial_num))
    loop.run()


def _mtrue_any(v: Any) -> bool:
    """``any(stitchMIP)``"""
    if v is None:
        return False
    return bool(np.any(np.asarray(v, dtype=bool)))


def _local_worker_num(local_worker_num: Any, s: dict[str, Any]) -> float:
    """the frame's local workers: the given number, or 1, or for large datasets without the
    cluster up to 8, one per 3 CPUs of this process (MATLAB's: of maxNumCompThreads)"""
    if not _empty(local_worker_num):
        v = float(np.asarray(local_worker_num, dtype=np.float64).ravel()[0])
        return int(v) if v.is_integer() else v
    info = s["info"]
    tab = info.tab
    n = 1
    if s["parse_cluster"] or tab.height == 0:
        return n
    did = int(tab["did"][0]) if s["multi_loc"] else 1
    folder = s["paths"][did - 1]
    tile = f"{folder}/{s['processed_dir_str']}/{tab['Filename'][0]}"
    mask = (
        (tab["ch"] == info.ch[0])
        & (tab["camera"] == info.cam[0])
        & (tab["fullIter"] == info.full_iter[0])
        & (tab["stack"] == info.stackn[0])
    )
    if _mtrue(info.zlayer_stitch):
        mask &= tab["z"] == info.nz[0]
    tile_num = int(np.sum(mask))
    if _exist_file(tile):
        voxels = math.prod(float(v) for v in _get_image_size(tile))
        if (s["xcorr_shift"] and tile_num >= LOCAL_WORKER_TILES) or (
            tab.height * voxels >= LOCAL_WORKER_VOXELS
        ):
            from petakit5d.core.threads import available_cpus

            n = max(1, min(8, available_cpus() // 3))
    return n


# ---------------------------------------------------------------------------------------------
# the rounds over the volumes
# ---------------------------------------------------------------------------------------------


class _Loop:
    """XR_matlab_stitching_wrapper's rounds over the volumes (time points x cameras x stacks x
    channels x z layers, in MATLAB's linear order for its numbers)"""

    def __init__(self, s: dict[str, Any], *, max_trial_num: float) -> None:
        info = s["info"]
        self.s = s
        self.info = info
        self.dims = (
            len(info.full_iter),
            len(info.cam),
            len(info.stackn),
            len(info.ch),
            len(info.nz),
        )
        shape = self.dims
        self.row_exist = np.ones(shape, dtype=bool)
        self.done = np.zeros(shape, dtype=bool)
        self.trials = np.zeros(shape, dtype=np.float64)
        self.max_trial = float(max_trial_num)
        self.active = np.zeros(shape, dtype=bool)
        self.job_ids = -np.ones(shape, dtype=np.int64)
        self.cluster = bool(s["parse_cluster"])
        self.im_size: tuple[int, ...] | None = None
        self.stream_counter = 0
        tab = info.tab
        if s["online_stitch"]:
            self.stream_max = 5
        else:
            self.stream_max = 100 * tab.height
        self.big_stitch_data = self._big_stitch_data()
        self.nice = 0
        self.job_mem_per_cpu = 20.9
        self.job_max_cpu_num = 24.0
        self.master_mem: float | None = None
        self._known: dict[tuple[int, ...], tuple[int, int]] = {}
        self._known_time = -math.inf
        self._task_ids: dict[tuple[int, ...], int] = {}
        self._tasks: dict[str, Any] = {}

    def _big_stitch_data(self) -> bool:
        s, info = self.s, self.info
        tab = info.tab
        if s["online_stitch"] or s["stitch_mip"] or s["streaming"] or tab.height <= 20:
            return False
        fn = tab["Filename"][0]
        did = int(tab["did"][0]) if s["multi_loc"] else 1
        folder = s["paths"][did - 1]
        if s["use_processed"]:
            tile = f"{folder}/{s['processed_dir_str']}/{fn}"
        else:
            tile = f"{folder}/{fn}"
        sz = _get_image_size(tile)
        return math.prod(float(v) for v in sz) * 4 * tab.height > 100 * 1024**3

    def _f(self, sub: tuple[int, ...]) -> int:
        """MATLAB's sub2ind (1-based linear index, the first subscript fastest)"""
        return int(np.ravel_multi_index(sub, self.dims, order="F")) + 1

    def _settled(self) -> bool:
        return bool(np.all(self.done | ((self.trials >= self.max_trial) & ~self.active)))

    def _last_f(self) -> int:
        cand = (~self.done & (self.trials < self.max_trial)).ravel(order="F")
        idx = np.flatnonzero(cand)
        return int(idx[-1]) + 1 if idx.size else 0

    def run(self) -> None:
        s = self.s
        mark = None
        if self.cluster:
            from petakit5d.jobs import single
            from petakit5d.jobs.config import cluster_config_value

            cfg = s["config_file"]
            self.job_mem_per_cpu = float(cluster_config_value(cfg, "memPerCPU", 20.9))
            self.job_max_cpu_num = float(cluster_config_value(cfg, "maxCPUNum", 24))
            if s["master_compute"] and os.environ.get("SLURM_JOB_ID"):
                from petakit5d.core import memory

                self.master_mem = float(memory.get_system_memory())
            from petakit5d.jobs.coordinator import coordinator_nice

            # the stitching jobs wait for jobs they submit: these go first
            self.nice = coordinator_nice()
            # the pools of task runners that the call's single jobs make end with the call
            mark = single.mark()
        try:
            self._rounds()
        finally:
            if mark is not None:
                from petakit5d.jobs import single

                single.end(mark)

    # -- one round ---------------------------------------------------------------------------
    def _rounds(self) -> None:
        s = self.s
        while not self._settled():
            # exit if no new images are transferred
            if s["streaming"] and self.stream_counter > self.stream_max:
                break
            last_f = self._last_f()
            # the files volumes wait for (cells of paths), and the results still to come
            awaited: list[list[str]] = []
            results: list[str] = []
            self._known = {}
            self._known_time = -math.inf
            nn, ncam, ns, nc, nzz = self.dims
            for n in range(nn):
                for cam in range(ncam):
                    for st in range(ns):
                        for c in range(nc):
                            for z in range(nzz):
                                self._volume((n, cam, st, c, z), last_f, awaited, results)
            if s["streaming"]:
                self.stream_counter += 1
            # wait if some volumes are still to come
            if not bool(np.all(self.done | (self.trials >= self.max_trial))) or (
                self.cluster and bool(np.any(self.active & ~self.done))
            ):
                if not s["streaming"]:
                    _wait([results, *awaited], ROUND_WAIT)
                else:
                    _wait([], ROUND_WAIT)

    def _volume(self, sub, last_f: int, awaited: list, results: list) -> None:  # noqa: C901
        s, info = self.s, self.info
        n, ncam, st, c, z = sub
        if self.done[sub] or (self.trials[sub] >= self.max_trial and not self.active[sub]):
            return
        if not self.row_exist[sub]:
            self.done[sub] = True
            return
        f = self._f(sub)
        tab = info.tab
        cam = info.cam[ncam]
        mask = (
            (tab["ch"] == info.ch[c])
            & (tab["camera"] == cam)
            & (tab["fullIter"] == info.full_iter[n])
            & (tab["stack"] == info.stackn[st])
        )
        zlayer = _mtrue(info.zlayer_stitch)
        if zlayer:
            mask &= tab["z"] == info.nz[z]
        cur = tab[np.asarray(mask, dtype=bool)]
        xcorr_mode = info.xcorr_mode
        if cur.height == 0:
            self.row_exist[sub] = False
            self.done[sub] = True
            # if the primary channel is missing, the other channels are skipped in 'primary'
            # mode (MATLAB marks the first time point, stack and z layer)
            if ncam == 0 and c == 0 and xcorr_mode == "primary":
                self.row_exist[0, :, 0, :, 0] = False
                self.done[0, :, 0, :, 0] = True
            return

        # 'primaryFirst': the number of tiles must be the primary channel's
        if (
            xcorr_mode.lower() == "primaryfirst"
            and info.primary_tab is not None
            and cur.height != info.primary_tab.height
            and not s["online_stitch"]
        ):
            _print("The number of tiles is different from the number of tiles in the primary "
                   "channel, skip it!")  # fmt: skip
            for fn in cur["Filename"]:
                _print(f"    {{'{fn}'}}")
            self.done[sub] = True
            return

        laser = _unique(cur["laser"])
        abstime = _unique(cur["abstime"])
        fpgatime = float(cur["fpgatime"][0])
        xyz = np.column_stack(
            [
                np.asarray(cur[k], dtype=np.float64)
                for k in ("StageX_um_", "StageY_um_", "StageZ_um_")
            ]
        )
        # the data path's number is part of the tile indices
        tile_idx = np.column_stack(
            [np.asarray(cur[k], dtype=np.float64) for k in ("x", "y", "z", "did")]
        )
        if laser.size > 1:
            laser = laser[:1]

        fnames = [str(x) for x in cur["Filename"]]
        nf = len(fnames)
        if s["multi_loc"]:
            did = [int(v) for v in cur["did"]]
            tile_fullpaths = [f"{s['paths'][d - 1]}/{x}" for d, x in zip(did, fnames, strict=True)]
            # (MATLAB's tile loop reuses the volume's number f: it is nF from here on)
            f = nf
        else:
            tile_fullpaths = [f"{s['paths'][0]}/{x}" for x in fnames]

        flipped: np.ndarray | None = None
        if s["parse_setting_file"]:
            from petakit5d.pipeline.settings import parse_image_setting_files

            setting = parse_image_setting_files(tile_fullpaths)
            flipped = np.array([float(x["StageInterval"]) < 0 for x in setting], dtype=bool)

        # streaming: the tiles must be there; processed tiles must be there
        if s["streaming"]:
            # (batch_file_exist sees tiles that arrived from another machine)
            if not all(_batch_exist(tile_fullpaths)):
                self.stream_counter += 1
                return
            self.stream_counter = 0
        processed_dir_str = s["processed_dir_str"]
        tile_processed = []
        if s["use_processed"]:
            if s["multi_loc"]:
                tile_processed = [
                    f"{s['paths'][d - 1]}/{processed_dir_str}/{x}"
                    for d, x in zip(did, fnames, strict=True)
                ]
                f = nf
            else:
                tile_processed = [f"{s['paths'][0]}/{processed_dir_str}/{x}" for x in fnames]
            exist = _batch_exist(tile_processed)
            if s["streaming"]:
                if not all(exist):
                    self.stream_counter += 1
                    return
                self.stream_counter = 0
            elif not all(exist):
                # (cluster) not a try while jobs of the pipeline still make the tiles
                missing = [p for p, e in zip(tile_processed, exist, strict=True) if not e]
                if self.cluster and processed_tiles_coming(missing):
                    awaited.append(missing)
                else:
                    self.trials[sub] += 1
                return
        else:
            processed_dir_str = ""

        z_str = _sprintf("_%03dz", info.nz[z]) if zlayer else ""
        if s["online_stitch"]:
            z_str = _sprintf("%s_%03dntile", z_str, nf)

        specify_cam = _mtrue(info.specify_cam)
        prefix = info.prefix
        fi = info.full_iter[n]
        sip = s["stitch_info_path"]
        stitch_info_fullpath = info.stitch_info_fullpath

        # a secondary channel waits for its primary's stitch information
        is_primary = True
        if xcorr_mode == "primaryFirst":
            if sub != (0, 0, 0, 0, 0):
                is_primary = False
        elif xcorr_mode == "stitchInfo":
            is_primary = False
        elif xcorr_mode in ("primary", "all"):
            stitch_info_fullpath = self._name(
                sip, fi, cam, info.ch[c], info.stackn[st], laser, abstime, fpgatime, z_str, ".mat"
            )
        if xcorr_mode.lower() == "primary" and not (ncam == 0 and c == 0):
            is_primary = False
            pmask = (
                (tab["ch"] == info.ch[0])
                & (tab["camera"] == info.cam[0])
                & (tab["fullIter"] == fi)
                & (tab["stack"] == info.stackn[st])
            )
            ptab = tab[np.asarray(pmask, dtype=bool)]
            if ptab.height == 0:
                raise StitchDatasetError(
                    "Index exceeds the number of array elements. Index must not exceed 0.",
                    "MATLAB:badsubscript",
                )
            p_laser = _unique(ptab["laser"])
            p_abstime = _unique(ptab["abstime"])
            p_fpgatime = float(ptab["fpgatime"][0])
            p_z = _unique(ptab["z"])[0]
            if p_laser.size > 1:
                p_laser = p_laser[:1]
            pz_str = _sprintf("_%03dz", p_z) if zlayer else ""
            if specify_cam:
                stitch_info_fullpath = _sprintf(
                    _NAME_CAM + ".mat", sip, prefix, fi, info.cam[0], info.ch[0],
                    info.stackn[st], p_laser, p_abstime, p_fpgatime, pz_str,
                )  # fmt: skip
            else:
                # (MATLAB takes the first stack here)
                stitch_info_fullpath = _sprintf(
                    _NAME_NOCAM + ".mat", sip, prefix, fi, info.ch[0], info.stackn[0],
                    p_laser, p_abstime, p_fpgatime, pz_str,
                )  # fmt: skip
        # (batch_file_exist sees a file another node just made)
        if not is_primary and not _batch_exist([stitch_info_fullpath])[0]:
            awaited.append([stitch_info_fullpath])
            return

        rt = s["stitching_rt"]
        save_fsname = self._name(
            rt, fi, cam, info.ch[c], info.stackn[st], laser, abstime, fpgatime, z_str, ""
        )
        tmp_fname = save_fsname + ".tmp"
        save_fname = save_fsname + ".zarr"
        from petakit5d.core.paths import fileparts

        fsname = fileparts(save_fname)[1]

        # online stitching: results of fewer tiles are removed
        if s["online_stitch"]:
            self._remove_fewer_tiles(save_fsname, z_str, nf, sub)

        # first check if the volume is done
        if _batch_exist([save_fname])[0]:
            self.done[sub] = True
            if os.path.isfile(tmp_fname):
                _delete(tmp_fname)
            return

        func_str = self._func_str(
            sub, tile_fullpaths, xyz, tile_idx, flipped, cur, save_fname, fsname,
            stitch_info_fullpath, processed_dir_str, is_primary,
        )  # fmt: skip

        job_running = False
        master_compute_f = s["master_compute"]
        stitch_cpus = s["cpus_per_task"]
        stitch_mem = 0.0
        if self.cluster:
            master_compute_f, stitch_cpus, stitch_mem = self._job_size(
                tile_fullpaths, tile_processed, xyz, tile_idx, flipped
            )
        primary_first_master = (
            master_compute_f and s["xcorr_shift"] and xcorr_mode.lower() == "primaryfirst"
            and is_primary
        )  # fmt: skip
        if os.path.isfile(tmp_fname) or (self.cluster and not primary_first_master):
            if self.cluster:
                if not (master_compute_f and f == last_f and self._job_id_at(f) == 0):
                    last_file = master_compute_f and f == last_f
                    job_running = self._submit(
                        sub, func_str, f, last_file, master_compute_f, stitch_cpus, stitch_mem,
                        save_fname,
                    )  # fmt: skip
            else:
                age = (time.time() - os.path.getmtime(tmp_fname)) / 60
                if age < nf * s["unit_wait_time"]:
                    # (another process stitches it)
                    results.append(save_fname)
                    return
                _touch(tmp_fname)
        else:
            _touch(tmp_fname)

        # the master takes over the last unfinished volume, unless its job runs (it would
        # write the same files at the same time; a pending job of it was cancelled)
        if not self.cluster or (
            master_compute_f
            and (
                (f == last_f and not job_running)
                or (primary_first_master and self.job_ids[sub] == -1)
            )
        ):
            self.trials[sub] += 1
            from petakit5d.jobs.tasks import run_record, task_record

            t0 = time.perf_counter()
            run_record(task_record(self._task(func_str), "python"))
            _print(f"\nStitching result {fsname} is generated. Elapsed time is "
                   f"{time.perf_counter() - t0:f} seconds.\n")  # fmt: skip

        if _batch_exist([save_fname])[0]:
            self.done[sub] = True
            if os.path.isfile(tmp_fname):
                _delete(tmp_fname)
        elif self.cluster and self.active[sub]:
            # (its job still makes it)
            results.append(save_fname)

    # -- names -------------------------------------------------------------------------------
    def _name(self, folder, fi, cam, ch, stack, laser, abstime, fpgatime, z_str, ext) -> str:
        prefix = self.info.prefix
        if _mtrue(self.info.specify_cam):
            return _sprintf(
                _NAME_CAM + ext, folder, prefix, fi, cam, ch, stack, laser, abstime, fpgatime,
                z_str,
            )  # fmt: skip
        return _sprintf(
            _NAME_NOCAM + ext, folder, prefix, fi, ch, stack, laser, abstime, fpgatime, z_str
        )

    def _remove_fewer_tiles(self, save_fsname: str, z_str: str, nf: int, sub) -> None:
        from petakit5d.core.paths import fileparts, matlab_dir
        from petakit5d.jobs.publish import remove_path

        rt = self.s["stitching_rt"]
        # (for bidirectional scans the fpga time may change)
        pattern = f"{save_fsname[:-32]}_*msecAbs_{z_str[1:5]}*ntile.zarr"
        names = [e.name for e in matlab_dir(pattern)]
        for name in names:
            m = re.search(r"_(\d+)ntile", name)
            num = float(m.group(1)) if m else math.nan
            full = f"{rt}/{name}"
            mip = f"{rt}/MIPs/{fileparts(full)[1]}_MIP_z.tif"
            if num < nf and _exist_file(full):
                if os.path.isdir(full):
                    remove_path(full)
                else:
                    _delete(full)
                _delete(mip)
                self.done[sub] = False

    # -- the function string -----------------------------------------------------------------
    def _func_str(
        self, sub, tile_fullpaths, xyz, tile_idx, flipped, cur, save_fname, fsname,
        stitch_info_fullpath, processed_dir_str, is_primary,
    ) -> str:  # fmt: skip
        s = self.s
        pr = s["pr"]
        tiles_str = "{'" + "','".join(tile_fullpaths) + "'}"
        xyz_str = _mat2str_commas(xyz)
        tile_idx_str = _mat2str_commas(tile_idx)
        tile_info_fullpath = ""
        flipped_str = "" if flipped is None else _num2str_commas(flipped, "%d,")

        # the process function of the first tile's channel pattern
        hits = [p in tile_fullpaths[0] for p in s["channels"]]
        if not any(hits) and "mappedFilename" in cur:
            mapped = str(cur["mappedFilename"][0])
            hits = [p in mapped for p in s["channels"]]
        rows = s["fun_rows"]
        for i, h in enumerate(hits):
            if h and i >= len(rows):
                raise StitchDatasetError(
                    f"Index in position 1 exceeds array bounds. Index must not exceed {len(rows)}.",
                    "MATLAB:badsubscript",
                )
        fun_str = "{'" + "','".join(r for r, h in zip(rows, hits, strict=False) if h) + "'}"

        # more than 10 tiles: the tiles' paths and coordinates in a file for the job
        if len(tile_fullpaths) > TILE_INFO_MIN and self.cluster:
            tile_info_fullpath = f"{s['stitching_rt']}/stitchInfo/{fsname}_tile_info.mat"
            if not os.path.exists(tile_info_fullpath):
                from petakit5d.io import write_mat

                _print("Save tile paths and coordinates to disk...")
                tmp = f"{s['stitching_rt']}/stitchInfo/{fsname}_tile_info_{s['uuid']}.mat"
                col = np.empty((len(tile_fullpaths), 1), dtype=object)
                for i, p in enumerate(tile_fullpaths):
                    col[i, 0] = p
                write_mat(
                    tmp,
                    {
                        "tile_fullpaths": col,
                        "xyz": xyz,
                        "tileIdx": tile_idx,
                        "flippedTile": None if flipped is None else flipped.reshape(1, -1),
                    },
                )
                os.replace(tmp, tile_info_fullpath)
            tiles_str, xyz_str, tile_idx_str, flipped_str = "{}", "[]", "[]", ""

        out_bbox = pr["outBbox"]
        if out_bbox is not None and np.ndim(out_bbox) == 2 and np.shape(out_bbox)[0] > 1:
            raise StitchDatasetError(
                "outBbox: a 3 x 2 box makes num2str give several rows, which garble MATLAB's "
                "function string; give it as [ymin xmin zmin ymax xmax zmax]"
            )
        v = {
            "tile_fullpaths": tiles_str,
            "xyz": xyz_str,
            "axis_order": s["axis_order"],
            "data_order": pr["dataOrder"],
            "xy_pixel_size": pr["xyPixelSize"],
            "dz": pr["dz"],
            "skew_angle": pr["skewAngle"],
            "reverse": pr["reverse"],
            "objective_scan": pr["objectiveScan"],
            "io_scan": pr["IOScan"],
            "result_path": save_fname,
            "tile_info_fullpath": tile_info_fullpath,
            "stitch_info_dir": s["stitch_info_dir"],
            "stitch_info_fullpath": stitch_info_fullpath,
            "processed_dir_str": processed_dir_str,
            "ds": pr["DS"],
            "dsr": pr["DSR"],
            "zarr_file": pr["zarrFile"],
            "large_file": pr["largeFile"],
            "pool_size": _mat2str_commas(pr["poolSize"]),
            "batch_size": _mat2str_commas(pr["batchSize"]),
            "block_size": _mat2str_commas(pr["blockSize"]),
            "shard_size": _mat2str_commas(pr["shardSize"]),
            "resample_type": pr["resampleType"],
            "resample_factor": _num2str_commas(pr["resampleFactor"], "%.10d,"),
            "input_bbox": _mat2str_commas(pr["inputBbox"]),
            "tile_out_bbox": _mat2str_commas(pr["tileOutBbox"]),
            "tile_offset": pr["tileOffset"],
            "blend_method": pr["blendMethod"],
            "overlap_type": pr["overlapType"],
            "xcorr_shift": pr["xcorrShift"],
            "xy_max_offset": pr["xyMaxOffset"],
            "z_max_offset": pr["zMaxOffset"],
            "xcorr_downsample": _mat2str_commas(pr["xcorrDownsample"]),
            "xcorr_thresh": pr["xcorrThresh"],
            "shift_method": pr["shiftMethod"],
            "axis_weight": _mat2str_commas(pr["axisWeight"]),
            "group_file": pr["groupFile"],
            "is_primary_ch": is_primary,
            "use_primary_coords": pr["usePrimaryCoords"],
            "out_bbox": _num2str_commas(out_bbox, "%d,"),
            "save_16bit": pr["save16bit"],
            "tile_idx": tile_idx_str,
            "flipped_tile": flipped_str,
            "process_fun_path": fun_str,
            "stitch_mip": _mat2str_commas(pr["stitchMIP"]),
            "big_stitch_data": self.big_stitch_data,
            "edge_artifacts": pr["edgeArtifacts"],
            "dist_bboxes": _mat2str_commas(pr["distBboxes"]),
            "save_mip": pr["saveMIP"],
            "parse_cluster": self.cluster,
            "uuid": s["uuid"],
            "mcc_mode": pr["mccMode"],
            "config_file": pr["configFile"],
        }
        func_str = frame_task_string(v)
        # (only then, so that the jobs of before the option still run)
        if not self.cluster and s["local_worker_num"] > 1:
            func_str = _sprintf("%s,'localWorkerNum',%d)", func_str[:-1], s["local_worker_num"])
        # (and masterCompute only when off, as the frame's default is true)
        if self.cluster and not s["master_compute"]:
            func_str = func_str[:-1] + ",'masterCompute',false)"
        return func_str

    # -- cluster jobs ------------------------------------------------------------------------
    def _job_size(self, tile_fullpaths, tile_processed, xyz, tile_idx, flipped):
        """(master_compute_f, cpus, memory) of the volume's job (XR_matlab_stitching_wrapper
        lines 713-776)"""
        s = self.s
        from petakit5d.jobs.config import cluster_config_value
        from petakit5d.jobs.coordinator import coordinator_resources, master_task_fits

        if self.im_size is None:
            # (MATLAB measures the first volume's first tile once, for every volume)
            src = tile_processed[0] if s["use_processed"] else tile_fullpaths[0]
            self.im_size = _get_image_size(src)
        datasize = math.prod(float(v) for v in self.im_size)
        total_size = datasize * len(tile_fullpaths)
        # assume for double
        total_dsize = total_size * 4 / 1024**3
        mem_factor = 2
        master_compute_f = bool(s["master_compute"]) and master_task_fits(
            total_dsize * mem_factor, s["cpus_per_task"], self.job_mem_per_cpu,
            self.job_max_cpu_num, None, self.master_mem,
        )  # fmt: skip

        cpus: float = s["cpus_per_task"]
        mem = total_dsize * mem_factor
        from petakit5d.stitch.geometry import stitch_volume_size_estimate

        est = dict(s["estimate"])
        _, _, nv_range, nb_range = stitch_volume_size_estimate(
            xyz, self.im_size, tile_idx=tile_idx, flipped_tile=flipped, **est
        )
        nv = [float(v) for v in nv_range[1]]
        small_2d = (s["stitch_mip"] or any(v == 1 for v in nv)) and len(
            tile_fullpaths
        ) < SMALL_2D_TILES
        if nb_range[0] >= COORDINATOR_BATCHES and not small_2d:
            dtype = "uint16"
            if not s["save_16bit"]:
                from petakit5d.io import get_image_dtype

                dtype = get_image_dtype(tile_fullpaths[0])
            from petakit5d.core.arrays import data_type_to_byte_number

            nbytes = data_type_to_byte_number(dtype)
            nz_, ny_, nx_ = nv
            mip_mem = ny_ * nx_ * nbytes / 2**30
            if nz_ * ny_ * nx_ * nbytes / 2**30 >= self.job_max_cpu_num * self.job_mem_per_cpu / 2:
                mip_mem = (ny_ * nx_ + ny_ * nz_ + nx_ * nz_) * nbytes / 2**30
            # (a block's job: 2 CPUs, 1 for big data)
            big = (
                self.big_stitch_data
                or s["large_file"]
                or (len(tile_fullpaths) > 4 and datasize * len(tile_fullpaths) * 4 > BIG_DATA_BYTES)
            )
            bsz = [float(v) for v in np.asarray(s["estimate"]["batch_size"]).ravel()]
            block_mem = math.prod(min(b, v) for b, v in zip(bsz, nv, strict=True)) * nbytes
            block_mem = block_mem / 2**30 * 20
            computes = cluster_config_value(s["config_file"], "masterCompute", s["master_compute"])
            cpus, mem = coordinator_resources(mip_mem, bool(computes), 2 - int(big), block_mem)
        return master_compute_f, cpus, mem

    def _job_id_at(self, f: int) -> int:
        """``job_ids(f)``: the job of linear index f (0 when out of range, never compared true
        then, as MATLAB's index would fail only where the rule can't reach it)"""
        flat = self.job_ids.ravel(order="F")
        return int(flat[f - 1]) if 1 <= f <= flat.size else -1

    def _known_status(self, sub, last_file: bool) -> int | None:
        """the job's state from one squeue per round for the volumes' jobs (made again after
        10 s); a pending job of the master's last volume is looked at by itself (it is
        cancelled on it)"""
        if last_file:
            return None
        from petakit5d.jobs import single

        now = time.monotonic()
        if now - self._known_time > 10:
            subs = [tuple(int(x) for x in t) for t in np.argwhere((self.job_ids > 0) & ~self.done)]
            self._known = {}
            if subs:
                ids = [int(self.job_ids[t]) for t in subs]
                arr = [self._task_ids.get(t, 0) % ARRAY_ID_MODULUS for t in subs]
                for t, j, st in zip(subs, ids, single.job_states(ids, arr), strict=True):
                    self._known[t] = (j, st)
            self._known_time = time.monotonic()
        hit = self._known.pop(tuple(sub), None)
        if hit is None or hit[0] != int(self.job_ids[sub]):
            return None
        return hit[1]

    def _task(self, func_str: str) -> Any:
        """frame_task of a function string (read once per string)"""
        task = self._tasks.get(func_str)
        if task is None:
            task = self._tasks[func_str] = frame_task(func_str)
        return task

    def _submit(self, sub, func_str, f, last_file, master_compute_f, cpus, mem, save_fname) -> bool:
        from petakit5d.jobs import single

        s = self.s
        job_id = int(self.job_ids[sub])
        task_id = f % ARRAY_ID_MODULUS
        known = self._known_status(sub, last_file) if job_id > 0 else None
        if self.trials[sub] < self.max_trial:
            job_id, job_status, submit_status = single.submit_single_job(
                self._task(func_str),
                job_id,
                task_id,
                job_log_dir=s["job_log_dir"],
                job_log_fname=s["job_log_fname"],
                job_error_fname=s["job_error_fname"],
                master_compute=master_compute_f,
                cpus_per_task=cpus,
                last_file=last_file,
                mem_allocate=mem,
                mcc_mode=s["mcc_mode"],
                config_file=s["config_file"],
                nice=self.nice,
                known_status=known,
                output=save_fname,
            )
        else:
            # the last try: only whether its job (or a task runner of another job that took
            # its task) is still there
            job_status = single.check_single_job_status(job_id, task_id, known_status=known)
            submit_status = 0
        self.job_ids[sub] = job_id
        self._task_ids[tuple(sub)] = task_id
        self.trials[sub] += submit_status
        self.active[sub] = submit_status > 0 or job_status >= 0
        return job_status > 0
