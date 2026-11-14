"""The microscope pipeline (``XR_microscopeAutomaticProcessing``): every file of the datasets
deskewed and rotated (with the flat-field correction) by ``deskew_rotate.file.deskew_rotate_file``
(XR_deskewRotateFrame), and each dataset stitched by ``stitch.dataset.stitch_dataset``
(XR_matlab_stitching_wrapper), here or in Slurm jobs, on existing data or while it is acquired.

``automatic_processing`` is XR_microscopeAutomaticProcessing.m at the oracle commit (decon is
no longer a step of it: ``decon.dataset.decon_dataset`` is its own wrapper):

- **Folders**: ``<data path>/DS/`` (without the combined deskew and rotation), ``<data
  path>/DSR/`` and, for stitching, ``<data path>//<stitch_result_dir_name>/`` (by default
  ``matlab_stitch`` or ``matlab_stitch_<blend method>``), each removed first with its
  ``overwrite`` flag (DS, DSR, stitch) and made group-writable. Stitching turns deskew and
  rotation on; without rotation the combined path is off, and the combined path turns both on.
  While streaming, data folders that don't exist yet are made.
- **Files**: XR_parseImageFilenames (``pipeline.filenames.parse_dataset_filenames``): the files of
  the channel patterns with their partial files, sizes and z flips (setting files with
  ``parse_setting_file``); in streaming mode files changed in the last ``min_modify_time``
  minutes are left for a later round.
- **Values in transit**: MATLAB hands each file to XR_deskewRotateFrame and each dataset to
  XR_matlab_stitching_wrapper as a function string that it evaluates, here and in its jobs,
  so every value goes through ``sprintf`` (``%.20d`` for the pixel size, z step, angle and
  lower limit, ``num2str`` with ``%0.10f`` for the offset, ``mat2str`` for the resample factor
  and the box, ``%.10f`` for the offsets of the registration, ``%d`` for the edge erosion and
  the tries, ``num2str`` with ``%.10d,`` for the stitching's resample factor). ``task_string``
  writes these strings and ``task_of`` reads them back as MATLAB's parser does
  (``params.literal``) and binds them to the port's functions through the registry, so the
  steps get MATLAB's values in both modes.
- **Two flags per file**: deskew (and rotation) and stitching; a file is finished when both
  are done (or the stitching is off). A file deleted during the run is finished. A deskew
  result is seen by ``batch_file_exist`` (results another node just made included); a
  stitched result by the pattern ``<stitch folder>/<name less 39 characters>*Abs.zarr``.
- **Local mode** (no Slurm, or ``parse_cluster`` False): each file in turn, here, with MATLAB's
  ``.tmp`` markers (a marker younger than ``unit_wait_time`` minutes is another process's file);
  after each deskew the dataset's stitching is called (it stitches the volumes whose tiles are
  there and leaves the others), and the round ends with a look for the stitched results.
- **Cluster mode**: one single job per file (``jobs.single``: never waiting, task runners from the
  second task of a kind, the call's pools ended with it) with MATLAB's memory requests (the
  combined path: ``size / 2^30 * 2 * (7 + 3 / prod(resample) + 8 * tries)`` GB; the separate
  path: XR_estimateComputingMemory's deskew estimate), and one stitching job per dataset,
  submitted after its deskew jobs (``size / 2^30 * 2 * 20`` GB, 2 instead of 20 for MIP
  stitching), whose tries and state are those of all the dataset's files. While a dataset's
  deskew jobs are pending or running, ``<DSR folder>/deskew_in_progress.tmp`` tells its
  stitching job that tiles are coming (written at the start of each round and again once it is
  a minute old, removed when they are done and when the call ends). Each job's last try is
  waited for. The wait between rounds is 0.3 s per unfinished file (3 to 30 s; 30 s locally),
  returning within a second of the results still to come.
- **Streaming** (``streaming``): new files are added each round; the run ends when every file
  is done, the newest file of each folder is ``max_modify_time`` minutes old and
  ``max_wait_loop_num`` idle rounds (30 s each) have passed. With ``generate_image_list``
  (``'from_encoder'``, ``'from_sqlite'``) a dataset's image list is made again before its
  stitching.

The MATLAB function's quirks are kept (each has a test): the data paths are sorted (unique)
but the image lists keep their order; the axis-order check never raises; the stitching
folders are made (and deskew and rotation forced on) even when a missing image list turns the
stitching off; local stitching sets ``max_trial_num`` to 1 for the rest of the call (the deskew
tries too); a file's z flip becomes the default flip of files found later while streaming; a new
file whose name is part of a file name already listed (in any folder) is not added; the
stitching's resample factor is the deskew's ``[x y z]`` factor, which the stitching reads in the
tiles' axis order (``[y x z]``); a single ``save_16bit`` value fails when the stitching's string
is made (MATLAB's ``save16bit(2)``); the job log folder is checked with the last data path; the
stitched result's ``.tmp`` marker is ``<name>..tmp``. In local mode a dataset whose stitching
never makes its result is called again in every round (MATLAB counts no stitching tries
there).
"""

from __future__ import annotations

import math
import os
import time
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "DESKEW_FORMAT",
    "STITCH_FORMAT",
    "AutomaticProcessingError",
    "automatic_processing",
    "check_resample_setting",
    "task_of",
    "task_string",
]

ARRAY_ID_MODULUS = 5000
DESKEW_MARKER = "deskew_in_progress.tmp"
# the job functions of the strings' MATLAB functions
TASK_FUNCS = {
    "XR_deskewRotateFrame": "petakit5d.deskew_rotate.file:deskew_rotate_file",
    "XR_matlab_stitching_wrapper": "petakit5d.stitch.dataset:stitch_dataset",
}

# XR_microscopeAutomaticProcessing's function strings (lines 527-532 and 649-655), MATLAB's
# quotes undone
DESKEW_FORMAT = (
    "XR_deskewRotateFrame(%s,%.20d,%.20d,'skewAngle',%.20d,"
    "'objectiveScan',%s,'zStageScan',%s,'reverse',%s,'inputAxisOrder','%s',"
    "'outputAxisOrder','%s','FFCorrection',%s,'BKRemoval',%s,'lowerLimit',%.20d,"
    "'constOffset',[%s],'FFImage','%s','BackgroundImage','%s','rotate',%s,"
    "'resampleFactor',[%s],'DSRCombined',%s,'inputBbox',%s,'flipZstack',%s,"
    "'save16bit',%s,'save3DStack',%s,'saveZarr',%s,'uuid','%s')"
)
STITCH_FORMAT = (
    "XR_matlab_stitching_wrapper('%s','%s','ResultDir','%s',"
    "'streaming',%s,'DS',%s,'DSR',%s,'channelPatterns',%s,'ProcessedDirStr','%s',"
    "'axisOrder','%s','resampleType','%s','resampleFactor',[%s],'reverse',%s,"
    "'parseSettingFile',%s,'xcorrShift',%s,'xcorrMode','%s','xyMaxOffset',%.10f,"
    "'zMaxOffset',%.10f,'edgeArtifacts',%d,'blendMethod','%s',"
    "'outBbox',[%s],'save16bit',%s,'primaryCh','%s','stitchMIP',%s,'onlineStitch',%s,"
    "'zarrFile',%s,'maxTrialNum',%d,'parseCluster',%s,'mccMode',%s,'configFile','%s','uuid','%s')"
)


class AutomaticProcessingError(ValueError):
    """An error where MATLAB raises one, with its identifier ('' where MATLAB's has none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


def _warn(message: str, identifier: str = "") -> None:
    w = UserWarning(message)
    w.identifier = identifier
    warnings.warn(w, stacklevel=3)


def _ispc() -> bool:
    # (asked at each call, so that tests can take MATLAB's Windows branch on any system)
    return os.name == "nt"


# ---------------------------------------------------------------------------------------------
# values in MATLAB's strings
# ---------------------------------------------------------------------------------------------


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    return np.size(np.asarray(v, dtype=object)) == 0


def _string(v: Any) -> str:
    """``string(x)`` of a logical (or of [], which prints nothing)"""
    if _empty(v):
        return ""
    a = np.asarray(v)
    if a.size != 1:
        raise AutomaticProcessingError(
            f"A logical value {v!r} with several elements garbles MATLAB's function string"
        )
    return "true" if bool(a.reshape(-1)[0]) else "false"


def _number(v: Any) -> Any:
    """a numeric MATLAB value for sprintf (a scalar as a Python number, else its array)"""
    if v is None:
        return np.zeros((0, 0))
    a = np.asarray(v)
    if a.dtype.kind == "b":
        a = a.astype(np.float64)
    return a.reshape(-1)[0].item() if a.size == 1 else a


def _mat2str_comma(v: Any) -> str:
    """``mat2str_comma(v)``: ``strrep(mat2str(v), ' ', ',')`` ([] for None)"""
    from petakit5d.core.text import mat2str

    if v is None:
        return "[]"
    a = np.asarray(v)
    if a.dtype.kind == "f" and a.dtype != np.float64:
        a = a.astype(np.float64)
    return mat2str(a).replace(" ", ",")


def _num2str(v: Any, fmt: str) -> str:
    """``num2str(v, fmt)`` of a row ('' for [])"""
    from petakit5d.core.text import num2str

    if v is None or np.size(v) == 0:
        return ""
    a = np.asarray(v, dtype=np.float64)
    if a.ndim > 1 and a.shape[0] > 1:
        raise AutomaticProcessingError(
            f"{v!r}: num2str of a matrix gives several rows, which garble MATLAB's function string"
        )
    return num2str(a.reshape(1, -1), fmt)


def _cell_text(items: Sequence[str]) -> str:
    """``sprintf('{''%s''}', strjoin(items, ''','''))``"""
    return "{'" + "','".join(items) + "'}"


def task_string(fmt: str, *args: Any) -> str:
    """``sprintf(fmt, args...)`` as MATLAB formats the pipeline's strings: one call, so a value
    of several elements takes the conversions after its own, as in MATLAB."""
    from petakit5d.core.text import sprintf

    return sprintf(fmt, *args)


def _split_call(func_str: str) -> tuple[str, list[Any]]:
    from petakit5d.decon.dataset import split_call

    return split_call(func_str)


def task_of(func_str: str) -> Any:
    """The job task of one of the pipeline's function strings: the port's function of the
    string's MATLAB function, with the values MATLAB's evaluation of the string gives (names,
    case, prefixes and conversions through the registry)."""
    from petakit5d.jobs.tasks import Task
    from petakit5d.params.registry import REGISTRY

    name, args = _split_call(func_str)
    spec = REGISTRY.find(name)
    kwargs = spec.bind_matlab(args)
    hook = spec.matlab_in
    if hook is not None:
        kwargs = hook(kwargs)
    return Task(TASK_FUNCS[spec.matlab_name], kwargs=kwargs)


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def check_resample_setting(
    resample_type: str,
    resample: Any,
    objective_scan: bool,
    skew_angle: float,
    xy_pixel_size: float,
    dz: float,
) -> tuple[Any, float]:
    """``XR_checkResampleSetting``: the resample factor (MATLAB's ``[x y z]`` of the deskew,
    as a list; None for a type it doesn't know and no factor) and ``zAniso``. 'given' completes
    one value to three and two to ``[r1 r1 r2]``; 'isotropic' is ``[1 1 1]``; 'xy_isotropic'
    ``[1 1 zf]`` with ``zf`` the ratio that makes the voxels of the rotation isotropic."""
    from petakit5d.matlabcompat import elementary as E

    if objective_scan:
        z_aniso = dz / xy_pixel_size
    else:
        z_aniso = float(E.sind(skew_angle)) * dz / xy_pixel_size
    r = None if resample is None else [float(x) for x in np.asarray(resample).reshape(-1)]
    if resample_type == "given":
        if not r:
            raise AutomaticProcessingError(
                'For resample Type "given", the parameter resample must not be empty!'
            )
        if len(r) == 1:
            r = [r[0]] * 3
        elif len(r) == 2:
            r = [r[0], r[0], r[1]]
    elif resample_type == "isotropic":
        r = [1.0, 1.0, 1.0]
    elif resample_type == "xy_isotropic":
        theta = skew_angle * math.pi / 180
        st = float(E.sin(theta))
        ct = float(E.cos(theta))
        zf = float(E.sqrt((st**2 + z_aniso**2 * ct**2) / (ct**2 + z_aniso**2 * st**2)))
        r = [1.0, 1.0, zf]
    return r, z_aniso


def _batch_exist(path: str) -> bool:
    from petakit5d.jobs.files import batch_file_exist

    return batch_file_exist([path])[0]


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _touch(path: str) -> None:
    """``fclose(fopen(path, 'w'))``"""
    try:
        with open(path, "w"):
            pass
    except OSError:
        raise AutomaticProcessingError(
            "Invalid file identifier. Use fopen to generate a valid file identifier.",
            "MATLAB:FileIO:InvalidFid",
        ) from None


def _delete(path: str) -> None:
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


def _make_folder(path: str) -> None:
    from petakit5d.core.paths import mkdir_recursive
    from petakit5d.io.tiff import group_writable

    if not os.path.isdir(path):
        mkdir_recursive(path)
        group_writable(path)


def _remove_folder(path: str) -> None:
    from petakit5d.jobs.publish import remove_path

    if os.path.isdir(path):
        remove_path(path)


def _dir_names(pattern: str) -> list[str]:
    from petakit5d.core.paths import matlab_dir

    return [e.name for e in matlab_dir(pattern)]


def _stitch_prefix(fsname: str) -> str:
    """``fsname(1 : end - 39)``: the name without its msecAbs value and tile fields"""
    return fsname[: max(0, len(fsname) - 39)]


def _pause(seconds: float) -> None:
    from petakit5d.jobs import files

    files._pause(seconds)


def _wait_for_files(groups: list[list[str]], wait_time: float) -> None:
    from petakit5d.jobs import files

    files.wait_for_files(groups, wait_time)


def update_deskew_markers(marker_fullpaths: Sequence[str], fd_inds, deskewing) -> None:
    """The marker of each dataset with a file whose deskewing is pending (written again once
    it is a minute old, so that it stays recent while the jobs wait in the queue); the
    others' removed."""
    ds = [False] * len(marker_fullpaths)
    for d, on in zip(fd_inds, deskewing, strict=True):
        if on:
            ds[int(d)] = True
    for d, fn in enumerate(marker_fullpaths):
        try:
            mtime = os.path.getmtime(fn)
        except OSError:
            mtime = None
        if ds[d]:
            if mtime is None or (time.time() - mtime) / 60 >= 1:
                try:
                    with open(fn, "w"):
                        pass
                except OSError:
                    pass
        elif mtime is not None:
            _delete(fn)


# ---------------------------------------------------------------------------------------------
# the wrapper
# ---------------------------------------------------------------------------------------------


@register(
    "XR_microscopeAutomaticProcessing",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folder(s) of the same settings."),
    Param(
        "overwrite",
        "@(x) (numel(x) == 1 || numel(x) == 3) && islogical(x)",
        "Remove the results first (one value, or DS, DSR and stitch).",
        convert="same",
    ),
    Param("streaming", "@islogical", "Real-time processing: wait for files still to come."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files."),
    Param("skewAngle", "@isscalar", "Skew angle (degrees)."),
    Param("dz", "@isscalar", "Scan step (um)."),
    Param("xyPixelSize", "@isscalar", "Pixel size (um)."),
    Param("reverse", "@islogical", "Reverse the z axis."),
    Param("objectiveScan", "@islogical", "Objective scan."),
    Param("zStageScan", "@islogical", "Z stage scan."),
    Param(
        "save16bit",
        "@(x) (numel(x) == 1 || numel(x) == 2) && islogical(x)",
        "uint16 results of deskew/rotation and of stitching.",
        convert="same",
    ),
    Param("dzFromEncoder", "@islogical", "Estimate dz from the encoder positions."),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("saveZarr", "@islogical", "Save the results as zarr."),
    Param("save3DStack", "@islogical", "Save the DS and DSR volumes."),
    Param("deskew", "@islogical", "Deskew the files."),
    Param("rotate", "@islogical", "Rotate the deskewed files."),
    Param("stitch", "@islogical", "Stitch the deskewed and rotated files."),
    Param("parseSettingFile", "@islogical", "Flip z stacks after their setting files."),
    Param("flipZstack", "@islogical", "Flip the z stacks."),
    Param("inputAxisOrder", "@ischar", "Axis order of the input.", convert="axis_string"),
    Param("outputAxisOrder", "@ischar", "Axis order of the output.", convert="axis_string"),
    Param("DSRCombined", "@islogical", "Combined deskew and rotation."),
    Param("FFCorrection", "@islogical", "Flat-field correction."),
    Param("BKRemoval", "@islogical", "Background removal in the flat-field correction."),
    Param("lowerLimit", "@isnumeric", "Lowest value of the normalized flat field."),
    Param(
        "constOffset",
        "@(x) isnumeric(x)",
        "Constant offset after the flat field (None: the background).",
        convert="same",
    ),
    Param("FFImagePaths", "@iscell", "Flat-field images, in the order of channel_patterns."),
    Param("backgroundPaths", "@iscell", "Background images, in the order of channel_patterns."),
    Param("resampleType", "@ischar", "'given', 'isotropic' or 'xy_isotropic'.", convert="same"),
    Param(
        "resampleFactor",
        "@isnumeric",
        "Resampling (z, y, x) of the deskew (None: by resample_type).",
        convert="resample_xyz",
    ),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) of the input.", convert="bbox"),
    Param("stitchResultDirName", "@ischar", "Stitching result folder (empty: by blend method)."),
    Param("imageListFullpaths", "@(x) ischar(x) || iscell(x)", "Image list(s), one per dataset."),
    Param(
        "axisOrder",
        "@(x) ischar(x)",
        "Axis order of the stage coordinates ('x,y,z', '-x,y,z', ...).",
        convert="same",
    ),
    Param("blendMethod", "@ischar", "'none', 'feather', 'mean', 'median' or 'max'."),
    Param("xcorrShift", "@islogical", "Register the tiles."),
    Param(
        "xcorrMode",
        "@(x) ismember(lower(x), {'primary', 'primaryfirst', 'all'})",
        "'primaryFirst', 'primary' or 'all'.",
    ),
    Param("xyMaxOffset", "@isnumeric", "Largest shift in x and y (voxels).", convert="same"),
    Param("zMaxOffset", "@isnumeric", "Largest shift in z (voxels)."),
    Param("edgeArtifacts", "@isnumeric", "Edge erosion of the tiles (voxels)."),
    Param("primaryCh", "@ischar", "The primary channel of the registration."),
    Param(
        "stitchMIP",
        "@(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3))",
        "Stitch the tiles' MIPs (flags z, y, x; one value: z).",
        convert="flags_yxz",
    ),
    Param("onlineStitch", "@(x) islogical(x)", "Stitch the tiles there are."),
    Param("generateImageList", "@(x) ischar(x)", "'', 'from_encoder' or 'from_sqlite'."),
    Param("parseCluster", "@islogical", "Process in Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a job."),
    Param("uuid", "@ischar", "Id of the temporary names."),
    Param("maxTrialNum", "@isnumeric", "Tries per file."),
    Param("unitWaitTime", "@isnumeric", "Minutes a local .tmp marker counts as another's."),
    Param("minModifyTime", "@isnumeric", "Minutes a file must be unchanged (streaming)."),
    Param("maxModifyTime", "@isnumeric", "Minutes to wait for new files (streaming)."),
    Param("maxWaitLoopNum", "@isnumeric", "Idle rounds before the end (streaming)."),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; passed on, ignored)."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="P02",
)
def automatic_processing(
    data_paths,
    *,
    overwrite=False,
    streaming=True,
    channel_patterns=("CamA_ch0", "CamA_ch1", "CamB_ch0"),
    skew_angle=32.45,
    dz=0.5,
    xy_pixel_size=0.108,
    reverse=True,
    objective_scan=False,
    z_stage_scan=False,
    save_16bit=(True, True),
    dz_from_encoder=False,
    zarr_file=False,
    save_zarr=False,
    save_3d_stack=True,
    deskew=True,
    rotate=True,
    stitch=False,
    parse_setting_file=False,
    flip_zstack=False,
    input_axis_order="zyx",
    output_axis_order="zyx",
    dsr_combined=True,
    ff_correction=False,
    bk_removal=False,
    lower_limit=0.4,
    const_offset=None,
    ff_image_paths=("", "", ""),
    background_paths=("", "", ""),
    resample_type="isotropic",
    resample_factor=None,
    input_bbox=None,
    stitch_result_dir_name="",
    image_list_fullpaths="",
    axis_order="xyz",
    blend_method="none",
    xcorr_shift=False,
    xcorr_mode="primaryFirst",
    xy_max_offset=300,
    z_max_offset=50,
    edge_artifacts=2,
    primary_ch="",
    stitch_mip=None,
    online_stitch=False,
    generate_image_list="",
    parse_cluster=True,
    master_compute=True,
    job_log_dir="../job_logs",
    cpus_per_task=1,
    uuid="",
    max_trial_num=3,
    unit_wait_time=1,
    min_modify_time=1,
    max_modify_time=10,
    max_wait_loop_num=10,
    mcc_mode=False,
    config_file="",
) -> None:
    """Deskew/rotate and stitch the datasets' files (XR_microscopeAutomaticProcessing): see the
    module. Python values: boxes (z0, y0, x0, z1, y1, x1), resample factors (z, y, x), axis
    orders 'zyx', MIP flags (z, y, x)."""
    from petakit5d.core.paths import matlab_sort_key
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    m = {k: to_plain(v) for k, v in matlab_results(automatic_processing, dict(locals())).items()}
    paths = m["dataPaths"]
    paths = [paths] if isinstance(paths, str) else [str(p) for p in paths]
    # unique: sorted, each once
    paths = sorted(set(paths), key=matlab_sort_key)
    lists = m["imageListFullpaths"]
    lists = [lists] if isinstance(lists, str) else [str(p) for p in lists]

    s = _State(m, paths, lists)
    s.setup()
    s.run()


class _State:
    """XR_microscopeAutomaticProcessing's variables and rounds (MATLAB's names in snake_case;
    values in MATLAB's form, file indices 0-based)"""

    def __init__(self, m: dict[str, Any], paths: list[str], lists: list[str]) -> None:
        self.m = m
        self.paths = paths
        self.lists = lists
        self.nd = len(paths)
        ow = [bool(x) for x in np.atleast_1d(np.asarray(m["overwrite"])).reshape(-1)]
        self.overwrite = ow * 3 if len(ow) == 1 else ow
        self.streaming = bool(m["streaming"])
        self.channels = [str(c) for c in m["channelPatterns"]]
        self.save16 = [bool(x) for x in np.atleast_1d(np.asarray(m["save16bit"])).reshape(-1)]
        self.deskew = bool(m["deskew"])
        self.rotate = bool(m["rotate"])
        self.stitch = bool(m["stitch"])
        self.dsr_combined = bool(m["DSRCombined"])
        self.flip_zstack = bool(m["flipZstack"])
        self.parse_cluster = bool(m["parseCluster"])
        self.max_trial = float(_number(m["maxTrialNum"]))
        self.mcc_mode = m["mccMode"]
        self.stitch_dir = m["stitchResultDirName"]
        self.axis_order = m["axisOrder"]
        self.resample_type = m["resampleType"]
        self.ext = ".zarr" if m["saveZarr"] else ".tif"
        self.ds_paths: list[str] = []
        self.dsr_paths: list[str] = []
        self.stch_paths: list[str] = []
        self.job_log_dir = m["jobLogDir"]
        self.job_log_fname = ""
        self.job_error_fname = ""

    # -- folders and settings (lines 204-391) ------------------------------------------------
    def setup(self) -> None:
        m = self.m
        nd = self.nd
        last = ""
        for d in range(nd):
            last = self.paths[d]
            if not self.paths[d].endswith("/"):
                self.paths[d] = self.paths[d] + "/"
        if self.streaming:
            for d in range(nd):
                last = self.paths[d]
                _make_folder(last)
        if _ispc():
            self.paths = [p.replace("\\", "/") for p in self.paths]
            self.lists = [p.replace("\\", "/") for p in self.lists]

        if self.parse_cluster:
            from petakit5d.jobs.slurm import check_slurm_cluster

            chk = check_slurm_cluster(last, self.job_log_dir)
            self.parse_cluster = bool(chk.parse_cluster)
            self.job_log_fname = chk.job_log_fname
            self.job_error_fname = chk.job_log_error_fname

        if self.stitch:
            if not self.streaming:
                if len(self.lists) != nd:
                    _warn(
                        "The number of image list files does not match that of the data paths, "
                        "please make sure image list files are provided for each dataset!"
                    )
                    self.stitch = False
                else:
                    for d in range(nd):
                        if not _exist(self.lists[d]):
                            _warn(
                                f"Image list filename {self.lists[d]} does not exist, set stitch "
                                "option as False"
                            )
                            self.stitch = False
                            break
            self.deskew = True
            self.rotate = True
            for d in range(nd):
                if not self.stitch_dir:
                    if m["blendMethod"] == "none":
                        self.stitch_dir = "matlab_stitch"
                    else:
                        self.stitch_dir = f"matlab_stitch_{m['blendMethod']}"
                stch = f"{self.paths[d]}/{self.stitch_dir}/"
                if self.overwrite[2]:
                    _remove_folder(stch)
                _make_folder(stch)
                self.stch_paths.append(stch)
            # (MATLAB's check: ~regexpi(...) is false for a match and for none, so it never
            # raises)
            self.axis_order = self.axis_order.replace(" ", "")

        if not self.rotate:
            self.dsr_combined = False
        if self.dsr_combined:
            self.deskew = True
            self.rotate = True

        if self.deskew:
            if not self.dsr_combined:
                for d in range(nd):
                    ds = f"{self.paths[d]}/DS/"
                    if self.overwrite[0]:
                        _remove_folder(ds)
                    _make_folder(ds)
                    self.ds_paths.append(ds)
            if m["FFCorrection"]:
                ffs = list(m["FFImagePaths"])
                bgs = list(m["backgroundPaths"])
                if len(self.channels) != len(ffs) or len(self.channels) != len(bgs):
                    raise AutomaticProcessingError(
                        "The number of channels in channelPatterns does not match the number of "
                        "LS files or Background Files!"
                    )
                for c in range(len(ffs)):
                    if not _exist(ffs[c]):
                        raise AutomaticProcessingError(f"LS Image file {ffs[c]} does not exist!")
                    if not _exist(bgs[c]):
                        # (MATLAB's message names the background an LS image too)
                        raise AutomaticProcessingError(f"LS Image file {bgs[c]} does not exist!")
        else:
            self.rotate = False

        if self.rotate:
            for d in range(nd):
                dsr = f"{self.paths[d]}/DSR/"
                if self.overwrite[1]:
                    _remove_folder(dsr)
                _make_folder(dsr)
                self.dsr_paths.append(dsr)

        rtype = self.resample_type
        if not _empty(m["resampleFactor"]):
            rtype = "given"
        self.resample, _ = check_resample_setting(
            rtype,
            m["resampleFactor"],
            bool(m["objectiveScan"]),
            float(_number(m["skewAngle"])),
            float(_number(m["xyPixelSize"])),
            float(_number(m["dz"])),
        )

        dz = _number(m["dz"])
        if m["dzFromEncoder"] and not self.streaming:
            from petakit5d.stitch.imagelist_generate import estimate_actual_step_size_from_encoder

            self.dz_all = [
                estimate_actual_step_size_from_encoder(self.paths[d], dz=dz) for d in range(nd)
            ]
        else:
            self.dz_all = [dz] * nd

    # -- files ---------------------------------------------------------------------------
    def _parse(self):
        from petakit5d.pipeline.filenames import parse_dataset_filenames

        return parse_dataset_filenames(
            self.paths,
            self.channels,
            bool(self.m["parseSettingFile"]),
            self.flip_zstack,
            False,
            [],
            self.streaming,
            float(_number(self.m["minModifyTime"])),
            bool(self.m["zarrFile"]),
        )

    # -- the rounds (lines 393-795) --------------------------------------------------------
    def run(self) -> None:
        parsed = self._parse()
        self.fnames = list(parsed.fnames)
        self.fd_inds = [int(v) for v in parsed.fd_inds]
        self.gfnames = [list(g) for g in parsed.gfnames]
        self.partialvols = [bool(v) for v in parsed.partialvols]
        self.data_sizes = [float(v) for v in parsed.data_sizes]
        self.flips = [bool(v) for v in parsed.flip_zstack]
        self.latest = [float(v) for v in parsed.latest_modify_times]
        n = len(self.fnames)
        self.done = [[False, False] for _ in range(n)]
        self.trials = [[0, 0] for _ in range(n)]
        self.active = [[False, False] for _ in range(n)]
        self.job_ids = [[-1, -1] for _ in range(n)]
        self.wait_loop = 0
        if not self.streaming and not self.stitch:
            for f in range(n):
                self.done[f][1] = True

        self.markers: list[str] = []
        mark = None
        if self.parse_cluster:
            from petakit5d.jobs import single

            mark = single.mark()
        if self.parse_cluster and self.stitch:
            self.markers = [p + DESKEW_MARKER for p in self.dsr_paths]
        try:
            self._rounds()
        finally:
            if self.markers:
                update_deskew_markers(self.markers, [], [])
            if mark is not None:
                from petakit5d.jobs import single

                single.end(mark)

    def settled(self) -> bool:
        return all(
            self.done[f][k] or (self.trials[f][k] >= self.max_trial and not self.active[f][k])
            for f in range(len(self.fnames))
            for k in (0, 1)
        )

    def _deskewing(self) -> list[bool]:
        return [
            not self.done[f][0] and (self.trials[f][0] < self.max_trial or self.active[f][0])
            for f in range(len(self.fnames))
        ]

    def _rounds(self) -> None:
        m = self.m
        ts = time.perf_counter()
        marker_time = 0.0
        max_modify = float(_number(m["maxModifyTime"]))
        max_wait_loops = float(_number(m["maxWaitLoopNum"]))
        while not self.settled() or (
            self.streaming
            and (
                len(self.fnames) == 0
                or any(t < max_modify for t in self.latest)
                or self.wait_loop < max_wait_loops
            )
        ):
            if self.markers:
                update_deskew_markers(self.markers, self.fd_inds, self._deskewing())
                marker_time = time.perf_counter() - ts
            for f in range(len(self.fnames)):
                if self.markers and time.perf_counter() - ts - marker_time >= 60:
                    update_deskew_markers(self.markers, self.fd_inds, self._deskewing())
                    marker_time = time.perf_counter() - ts
                if all(self.done[f]):
                    continue
                self._file(f)

            if not self.parse_cluster and self.stitch:
                # a later file's stitching can finish the dataset of files checked before it
                for f in range(len(self.fnames)):
                    if self.done[f][1]:
                        continue
                    pat = self._stitch_pattern(f)
                    if _dir_names(pat):
                        self.done[f][1] = True

            if not self.settled():
                self.wait_loop = 0
                wait_time = 30.0
                if self.parse_cluster:
                    left = sum(not all(x) for x in self.done)
                    wait_time = min(30.0, max(3.0, 0.3 * left))
                _wait_for_files([self._awaited()], wait_time)
                if not self.streaming:
                    continue
            elif self.streaming:
                if self.wait_loop < max_wait_loops:
                    self.wait_loop += 1
                    _pause(30)
            else:
                break
            self._new_files()

    def _stitch_pattern(self, f: int) -> str:
        from petakit5d.core.paths import fileparts

        fsname = fileparts(self.fnames[f])[1]
        return f"{self.stch_paths[self.fd_inds[f]]}/{_stitch_prefix(fsname)}*Abs.zarr"

    def _awaited(self) -> list[str]:
        """the results of the steps still to come, as a round looks for them"""
        from petakit5d.core.paths import fileparts

        out = []
        for f in range(len(self.fnames)):
            coming = [
                not self.done[f][k]
                and not (self.trials[f][k] >= self.max_trial and not self.active[f][k])
                for k in (0, 1)
            ]
            if not any(coming):
                continue
            fsname = fileparts(self.fnames[f])[1]
            fd = self.fd_inds[f]
            if coming[0] and not self.dsr_combined and self.ds_paths:
                out.append(f"{self.ds_paths[fd]}{fsname}{self.ext}")
            if coming[0] and self.rotate and self.dsr_paths:
                out.append(f"{self.dsr_paths[fd]}{fsname}{self.ext}")
            if coming[1] and self.stch_paths:
                out.append(self._stitch_pattern(f))
        return out

    def _new_files(self) -> None:
        """the files that came since the last look (streaming)"""
        parsed = self._parse()
        new_names = list(parsed.fnames)
        if not new_names:
            return
        # (MATLAB's contains(fnames, x): a name that is part of a listed name is not new)
        new = [not any(x in g for g in self.fnames) for x in new_names]
        if not any(new):
            return
        for i, k in enumerate(new):
            if not k:
                continue
            self.fnames.append(new_names[i])
            self.fd_inds.append(int(parsed.fd_inds[i]))
            self.partialvols.append(bool(parsed.partialvols[i]))
            self.gfnames.append(list(parsed.gfnames[i]))
            self.flips.append(bool(parsed.flip_zstack[i]))
            self.data_sizes.append(float(parsed.data_sizes[i]))
            self.done.append([False, False])
            self.trials.append([0, 0])
            self.active.append([False, False])
            self.job_ids.append([-1, -1])
        self.latest = [
            max(a, float(b)) for a, b in zip(self.latest, parsed.latest_modify_times, strict=True)
        ]

    # -- one file --------------------------------------------------------------------------
    def _file(self, f: int) -> None:
        from petakit5d.core.paths import fileparts

        m = self.m
        fname = self.fnames[f]
        fsname = fileparts(fname)[1]
        fd = self.fd_inds[f]
        data_path = self.paths[fd]
        frame = data_path + fname
        if not _exist(frame):
            self.done[f] = [True, True]
            return
        task_id = (f + 1) % ARRAY_ID_MODULUS
        ds_full = dsr_full = tmp = ""

        if self.deskew:
            if not self.dsr_combined:
                ds_full = f"{self.ds_paths[fd]}{fsname}{self.ext}"
                tmp = f"{self.ds_paths[fd]}/{fsname}.tmp"
            if self.rotate:
                dsr_full = f"{self.dsr_paths[fd]}{fsname}{self.ext}"
                tmp = f"{self.dsr_paths[fd]}/{fsname}.tmp"
            if self._deskewed(ds_full, dsr_full):
                self.done[f][0] = True
                if os.path.isfile(tmp):
                    _delete(tmp)
        else:
            self.done[f][0] = True

        if not self.done[f][0]:
            if not self._deskew(f, frame, data_path, task_id, ds_full, dsr_full, tmp):
                # (a young .tmp marker: another process's file)
                return

        # then stitching
        stch_full = stch_tmp = ""
        found = False
        if self.stitch:
            if self.done[f][1]:
                return
            pattern = self._stitch_pattern(f)
            names = _dir_names(pattern)
            if names:
                found = True
                stch_full = self.stch_paths[fd] + names[0]
                stch_tmp = stch_full[:-4] + ".tmp"
                if _exist(stch_full):
                    self.done[f][1] = True
                    if os.path.isfile(stch_tmp):
                        _delete(stch_tmp)
            gen = m["generateImageList"]
            if self.streaming and gen:
                firsts = [
                    i
                    for i in range(len(self.fnames))
                    if not self.done[i][1] and self.fd_inds[i] == fd
                ]
                if firsts and f == firsts[0]:
                    from petakit5d.stitch import imagelist_generate as G

                    if gen == "from_encoder":
                        self._set_list(
                            fd,
                            G.generate_image_list_from_encoder(
                                data_path, self.dz_all[fd], self.channels
                            ),
                        )
                    elif gen == "from_sqlite":
                        self._set_list(fd, G.generate_image_list_from_sqlite(data_path))
        else:
            self.done[f][1] = True

        if self.done[f][1]:
            return
        if not self.parse_cluster and (not _exist(frame) or not _exist(dsr_full)):
            return
        if not self.parse_cluster:
            self.mcc_mode = False
            # (MATLAB's maxTrialNum = 1: for the rest of the call, the deskew tries too)
            self.max_trial = 1.0
        func_str = self._stitch_string(fd)
        if self.parse_cluster:
            ds_inds = [i for i in range(len(self.fnames)) if self.fd_inds[i] == fd]
            if any(not self.done[i][0] and self.trials[i][0] == 0 for i in ds_inds):
                return
            first = ds_inds[0]
            mem_factor = 2 if _any_true(m["stitchMIP"]) else 20
            est = self.data_sizes[f] / 2**30 * 2 * mem_factor
            job_id, job_status, submit_status = self._submit(
                func_str, self.job_ids[first][1], (first + 1) % ARRAY_ID_MODULUS, est,
                self.trials[first][1], None,
            )  # fmt: skip
            # one stitching job per dataset: its tries and state are those of all its files
            trials = self.trials[first][1] + submit_status
            for i in ds_inds:
                self.job_ids[i][1] = job_id
                self.trials[i][1] = trials
                self.active[i][1] = submit_status > 0 or job_status >= 0
        else:
            self._run_here(func_str)
        if found and _exist(stch_full):
            self.done[f][1] = True
            if os.path.isfile(stch_tmp):
                _delete(stch_tmp)

    def _set_list(self, fd: int, fn: str) -> None:
        while len(self.lists) <= fd:
            self.lists.append("")
        self.lists[fd] = fn

    def _deskewed(self, ds_full: str, dsr_full: str) -> bool:
        # batch_file_exist sees a result another node just made, which exist() can miss
        return (self.dsr_combined or _batch_exist(ds_full)) and (
            not self.rotate or _batch_exist(dsr_full)
        )

    def _deskew(self, f, frame, data_path, task_id, ds_full, dsr_full, tmp) -> bool:
        """the file's deskew (and rotation) here or in its job; False if another process has
        the file (its young .tmp marker)"""
        m = self.m
        ff_image = bg_image = ""
        if m["FFCorrection"]:
            hits = [c in frame for c in self.channels]
            sel = [i for i, h in enumerate(hits) if h]
            if not sel:
                raise AutomaticProcessingError(
                    "Insufficient number of outputs from right hand side of equal sign to "
                    "satisfy assignment.",
                    "MATLAB:needMoreRhsOutputs",
                )
            ff_image = m["FFImagePaths"][sel[0]]
            bg_image = m["backgroundPaths"][sel[0]]
        # (MATLAB's flipZstack: also the default flip of the files found later)
        self.flip_zstack = self.flips[f]
        inputs = [frame]
        if self.partialvols[f]:
            inputs = [data_path + x for x in self.gfnames[f]]
        func_str = task_string(
            DESKEW_FORMAT,
            _cell_text(inputs),
            _number(m["xyPixelSize"]),
            _number(self.dz_all[self.fd_inds[f]]),
            _number(m["skewAngle"]),
            _string(m["objectiveScan"]),
            _string(m["zStageScan"]),
            _string(m["reverse"]),
            m["inputAxisOrder"],
            m["outputAxisOrder"],
            _string(m["FFCorrection"]),
            _string(m["BKRemoval"]),
            _number(m["lowerLimit"]),
            _num2str(m["constOffset"], "%0.10f"),
            ff_image,
            bg_image,
            _string(self.rotate),
            _mat2str_comma(self.resample),
            _string(self.dsr_combined),
            _mat2str_comma(m["inputBbox"]),
            _string(self.flip_zstack),
            _string(self.save16[0]),
            _string(m["save3DStack"]),
            _string(m["saveZarr"]),
            m["uuid"],
        )
        if os.path.exists(tmp) or self.parse_cluster:
            if self.parse_cluster:
                if not self.dsr_combined:
                    from petakit5d.tools.memory_estimate import estimate_computing_memory

                    est = float(np.asarray(estimate_computing_memory(frame, ("deskew",))[0])[0])
                else:
                    r = 1.0 if self.resample is None else math.prod(self.resample)
                    est = self.data_sizes[f] / 2**30 * 2 * (7 + 3 / r + self.trials[f][0] * 8)
                job_id, job_status, submit_status = self._submit(
                    func_str, self.job_ids[f][0], task_id, est, self.trials[f][0],
                    dsr_full or ds_full,
                )  # fmt: skip
                self.job_ids[f][0] = job_id
                self.trials[f][0] += submit_status
                self.active[f][0] = submit_status > 0 or job_status >= 0
            else:
                try:
                    age = (time.time() - os.path.getmtime(tmp)) / 60
                except OSError:
                    age = math.inf
                if age < float(_number(m["unitWaitTime"])):
                    return False
                _touch(tmp)
        else:
            _touch(tmp)
        if not self.parse_cluster:
            self._run_here(func_str)
            self.trials[f][0] += 1
        if self._deskewed(ds_full, dsr_full):
            self.done[f][0] = True
            if os.path.isfile(tmp):
                _delete(tmp)
        return True

    def _stitch_string(self, fd: int) -> str:
        m = self.m
        if len(self.save16) < 2:
            raise AutomaticProcessingError(
                "Index exceeds the number of array elements. Index must not exceed "
                f"{len(self.save16)}.",
                "MATLAB:badsubscript",
            )
        if fd >= len(self.lists):
            raise AutomaticProcessingError(
                "Index exceeds the number of array elements. Index must not exceed "
                f"{len(self.lists)}.",
                "MATLAB:badsubscript",
            )
        mip = m["stitchMIP"]
        return task_string(
            STITCH_FORMAT,
            self.paths[fd],
            self.lists[fd],
            self.stitch_dir,
            _string(self.parse_cluster and self.streaming),
            _string(False),
            _string(True),
            _cell_text(self.channels),
            "DSR",
            self.axis_order,
            "isotropic",
            _num2str(self.resample, "%.10d,").replace(" ", ""),
            _string(m["reverse"]),
            _string(m["parseSettingFile"]),
            _string(m["xcorrShift"]),
            m["xcorrMode"],
            _number(m["xyMaxOffset"]),
            _number(m["zMaxOffset"]),
            _number(m["edgeArtifacts"]),
            m["blendMethod"],
            "",
            _string(self.save16[1]),
            m["primaryCh"],
            _mat2str_comma(mip),
            _string(m["onlineStitch"]),
            _string(m["saveZarr"]),
            _number(self.max_trial),
            _string(self.parse_cluster),
            _string(self.mcc_mode),
            m["configFile"],
            m["uuid"],
        )

    def _run_here(self, func_str: str) -> None:
        """MATLAB's ``feval(str2func(['@()', func_str]))``"""
        from petakit5d.jobs.tasks import run_record, task_record

        t0 = time.perf_counter()
        run_record(task_record(task_of(func_str), "python"))
        print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.\n", flush=True)

    def _submit(self, func_str, job_id, task_id, mem, trials, output):
        """generic_single_job_submit_wrapper while tries are left, else the last try's
        check_single_job_status; returns (job_id, job_status, submit_status)"""
        from petakit5d.jobs import single

        if trials < self.max_trial:
            return single.submit_single_job(
                task_of(func_str),
                job_id,
                task_id,
                job_log_dir=self.job_log_dir,
                job_log_fname=self.job_log_fname,
                job_error_fname=self.job_error_fname,
                master_compute=bool(self.m["masterCompute"]),
                cpus_per_task=_number(self.m["cpusPerTask"]),
                mem_allocate=mem,
                mcc_mode=self.mcc_mode,
                config_file=self.m["configFile"],
                output=output,
            )
        return job_id, single.check_single_job_status(job_id, task_id), 0


def _any_true(v: Any) -> bool:
    if v is None:
        return False
    return bool(np.any(np.asarray(v)))
