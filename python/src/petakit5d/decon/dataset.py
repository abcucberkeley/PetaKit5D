"""Deconvolution of datasets (``XR_decon_data_wrapper``): every file of the datasets that
matches a channel pattern, deconvolved by ``decon.frame.decon_file`` (XR_RLdeconFrame3D) here or
in Slurm jobs.

``decon_dataset`` is XR_decon_data_wrapper.m at the oracle commit:

- **Folders and parameters**: the result folder ``<data path>/<result_dir_name>/`` of each
  dataset (removed first with ``overwrite``), made group-writable, with parameters.mat (v7.3) and
  parameters.json of the call (MATLAB's ``ip.Results``: the values as given).
- **PSFs**: with ``psf_gen``, each PSF is cleaned and resampled to ``dz`` (psf_gen_new 'masked',
  1.5) and, for methods other than OMW, cropped to its support (``crop_psf_support``), then
  written to ``<result folder>/psfgen/<psf>_<RLMethod>.tif`` of each dataset (a temporary name,
  then a rename). A PSF file that doesn't exist is an error.
- **Files**: XR_parseImageFilenames (``pipeline.filenames.parse_dataset_filenames``); a file's PSF,
  background, scale factor and (OMW) Wiener alpha and OTF threshold are those of the first channel
  pattern it matches (``regexpi``); iterations may be given per dataset. Results are
  ``<result folder>/<name>.tif`` (``.zarr`` with ``save_zarr``, forced by the in-place large-file
  method); an existing result is kept.
- **Values in transit**: MATLAB hands each file to XR_RLdeconFrame3D as a function string that it
  evaluates, here and in its jobs, so every number goes through ``sprintf`` (``%.10f`` for the
  pixel size and the z steps, ``%d`` for the background, angle, erosion, iterations, damp and
  scale factors and offset, which keeps 7 significant digits of a non-integer, ``%.20f`` for the
  OMW alpha and threshold, ``mat2str`` for the bounds and sizes). ``frame_task_string`` writes
  that string conversion by conversion and ``frame_task`` reads it back as MATLAB's parser does
  (``params.literal``) and binds it to ``decon_file``'s parameters, so the files get MATLAB's
  values in both modes.
- **First time point** (``edge_erosion`` > 0 with ``erode_by_ftp``, not for large files): the
  first file of each dataset saves its eroded mask (``Masks/<name>_eroded.zarr``; a finished
  file without its mask is deconvolved again), and the others wait for it and use it (their own
  erosion if its size differs, with MATLAB's warning).
- **Local mode** (no Slurm, or ``parse_cluster`` False): each file in turn, with MATLAB's
  ``.tmp`` markers (a marker younger than ``unit_wait_time`` minutes is another process's file),
  up to ``max_trial_num`` tries per file.
- **Cluster mode**: one single job per file (``jobs.single``; never waiting, task runners from
  the second task of a kind, the call's pools ended with it), 10 x 4 bytes per voxel of memory,
  ``config_file``; GPU jobs use ``gpu_config_file`` with at most 4 CPUs (the config's
  ``minCPUNum`` can raise it). Files are done when the result folder's listing has them; a
  file's job is looked at once per round (one squeue for all), its last try is waited for, and
  the wait between rounds is 0.3 s per unfinished file (3 to 30 s), returning within a second of
  the results still to come or of a mask files wait for.
- **Files per job**: with ``task_batch_num`` > 1, or for GPU jobs without it (up to 8 files, as
  long as there are 16 jobs: ``max(1, min(8, files // 16))``), files of the same settings (the
  function string after the file) and size share a job; files submitted again get a job each.
- **OMW on the cluster** (not large files): the back projector of each PSF is made once here,
  before the jobs, from the first dataset with a file of its channel (no figure), and copied to
  the other datasets.
- **In-place large files on the cluster**: a job per volume runs its block jobs (R07's
  ``decon_large_in_place``), with ``in_place_controller_resources`` for its CPUs and memory,
  ``coordinator_nice`` for its priority, the caller's ``master_compute`` (when False) and
  ``job_log_dir``; a volume another call's controller still makes is waited for, without a try.

The MATLAB function's quirks are kept (each has a test): the background defaults to 100;
``skew_angle`` reaches the files through ``%d``; ``mask_fullpaths`` empty reaches them as
``{''}``; the job log names come from the last dataset's folder; ``.zarr`` results have the marker
``<name>..tmp``; a finished zarr result without its first-time-point mask is not deleted (MATLAB's
``delete`` of a folder only warns). Where MATLAB's function string would be garbled (a vector
where it formats a scalar, a file of two channels), the port raises ``DeconDatasetError``.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import time
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "FRAME_FORMAT",
    "DeconDatasetError",
    "decon_dataset",
    "frame_task",
    "frame_task_string",
    "split_call",
]

# the job's function (XR_RLdeconFrame3D)
FRAME_FUNC = "petakit5d.decon.frame:decon_file"
FRAME_NAME = "XR_RLdeconFrame3D"
ARRAY_ID_MODULUS = 5000
# files per GPU job without task_batch_num: at most 8, while there are 16 jobs
AUTO_BATCH_MAX = 8
AUTO_BATCH_JOBS = 16
# CPUs of a GPU job (with a GPU config file)
GPU_JOB_CPUS = 4

# XR_decon_data_wrapper's function string of a file (lines 547-555), MATLAB's quotes undone
FRAME_FORMAT = (
    "XR_RLdeconFrame3D('%s',%.10f,%.10f,'%s','PSFfile','%s',"
    "'dzPSF',%.10f,'background',[%d],'skewAngle',%d,'flipZstack',%s,"
    "'edgeErosion',%d,'ErodeMaskfile','%s','SaveMaskfile',%s,'Rotate',%s,"
    "'deconIter',%d,'RLMethod','%s','wienerAlpha',%.20f,'OTFCumThresh',%.20f,"
    "'hannWinBounds',%s,'skewed',[%s],'debug',%s,'saveStep',%d,'psfGen',%s,'saveZarr',%s,"
    "'parseCluster',%s,'parseParfor',%s,'GPUJob',%s,'save16bit',%s,"
    "'largeFile',%s,'largeMethod','%s','batchSize',%s,'blockSize',%s,"
    "'dampFactor',%d,'scaleFactor',[%d],'deconOffset',%d,'maskFullpaths',%s,"
    "'uuid','%s','cpusPerTask',%d,'mccMode',%s,'configFile','%s','GPUConfigFile','%s')"
)
# the values of the conversions, in order (s: text, l: string() of a logical or [])
FRAME_FIELDS = (
    ("frame", "s"),
    ("xy_pixel_size", "n"),
    ("dz", "n"),
    ("decon_path", "s"),
    ("psf", "s"),
    ("dz_psf", "n"),
    ("background", "n"),
    ("skew_angle", "n"),
    ("flip_zstack", "l"),
    ("edge_erosion", "n"),
    ("erode_mask", "s"),
    ("save_maskfile", "l"),
    ("rotate", "l"),
    ("decon_iter", "n"),
    ("rl_method", "s"),
    ("wiener_alpha", "n"),
    ("otf_cum_thresh", "n"),
    ("hann_win_bounds", "s"),
    ("skewed", "l"),
    ("debug", "l"),
    ("save_step", "n"),
    ("psf_gen", "l"),
    ("save_zarr", "l"),
    ("parse_cluster", "l"),
    ("parse_parfor", "l"),
    ("gpu_job", "l"),
    ("save_16bit", "l"),
    ("large_file", "l"),
    ("large_method", "s"),
    ("batch_size", "s"),
    ("block_size", "s"),
    ("damp_factor", "n"),
    ("scale_factor", "n"),
    ("decon_offset", "n"),
    ("mask_fullpaths", "s"),
    ("uuid", "s"),
    ("cpus_per_task", "n"),
    ("mcc_mode", "l"),
    ("config_file", "s"),
    ("gpu_config_file", "s"),
)
_CONVERSION = re.compile(r"%[-+ 0#]*\d*(?:\.\d+)?[a-zA-Z]")


class DeconDatasetError(ValueError):
    """An error where MATLAB raises one, with its identifier ('' where MATLAB's has none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


def _warn(message: str, identifier: str = "") -> None:
    w = UserWarning(message)
    w.identifier = identifier
    warnings.warn(w, stacklevel=3)


# ---------------------------------------------------------------------------------------------
# the function string of a file and its values
# ---------------------------------------------------------------------------------------------


def _pieces(fmt: str) -> list[str]:
    """the format cut after each conversion's literal text: one conversion per piece"""
    out, start = [], 0
    for m in _CONVERSION.finditer(fmt):
        nxt = _CONVERSION.search(fmt, m.end())
        end = nxt.start() if nxt else len(fmt)
        out.append(fmt[start:end])
        start = end
    return out


_FRAME_PIECES = _pieces(FRAME_FORMAT)


def _matlab_string(v: Any) -> str:
    """``string(x)`` of a logical (or of [], which prints nothing)"""
    if v is None or (isinstance(v, (list, tuple, np.ndarray)) and np.size(v) == 0):
        return ""
    a = np.asarray(v)
    if a.size != 1:
        raise DeconDatasetError(
            f"A logical value {v!r} with several elements garbles MATLAB's function string"
        )
    return "true" if bool(a.reshape(-1)[0]) else "false"


def _scalar_number(name: str, v: Any) -> float | int | bool:
    a = np.asarray(v)
    if a.size != 1 or a.dtype.kind not in "biuf":
        raise DeconDatasetError(
            f"{name}: {v!r} is not one number; in MATLAB's function string it would shift the "
            "values after it (MATLAB's files would get other parameters)"
        )
    x = a.reshape(-1)[0]
    return bool(x) if a.dtype.kind == "b" else x.item()


def frame_task_string(values: dict[str, Any]) -> str:
    """XR_decon_data_wrapper's function string of a file, ``sprintf`` of ``FRAME_FORMAT``:
    ``values`` has the keys of ``FRAME_FIELDS`` (texts already in MATLAB's form: mat2str of the
    bounds and sizes, the masks' cell text). Each conversion is formatted with its own value, as
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
        raise DeconDatasetError(f"Not a function call: {func_str}")
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
    """The job task of a function string of XR_RLdeconFrame3D: ``decon_file`` with the values
    MATLAB's evaluation of the string gives (names, case and conversions through the registry)."""
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
    """``strrep(mat2str(v), ' ', ',')``"""
    from petakit5d.core.text import mat2str

    return mat2str(np.asarray(v, dtype=np.float64)).replace(" ", ",")


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _as_list(x: Any) -> list[Any]:
    if x is None:
        return []
    if isinstance(x, (str, os.PathLike)):
        return [os.fspath(x)]
    return [os.fspath(e) if isinstance(e, os.PathLike) else e for e in x]


def _vector(x: Any) -> np.ndarray:
    return np.atleast_1d(np.asarray(x, dtype=np.float64)).reshape(-1)


def _per_channel(x: Any, n_channels: int) -> np.ndarray:
    """MATLAB's ``x * ones(1, numel(channelPatterns))`` for a scalar with several channels"""
    v = _vector(x)
    if n_channels > 1 and v.size == 1:
        v = np.full(n_channels, v[0])
    return v


def _pick(v: np.ndarray, mapping: list[bool], what: str):
    """``v(psfMapping)``: the value of the file's channel (one, or MATLAB's string garbles)"""
    sel = [v[i] for i, m in enumerate(mapping) if m]
    if len(mapping) > v.size and any(mapping[v.size :]):
        raise DeconDatasetError(
            f"Index exceeds the number of array elements. Index must not exceed {v.size}.",
            "MATLAB:badsubscript",
        )
    if len(sel) != 1:
        raise DeconDatasetError(
            f"The file matches {len(sel)} channel patterns: its {what} would be several values "
            "in MATLAB's function string"
        )
    x = float(sel[0])
    return int(x) if x.is_integer() and abs(x) < 2**53 else x


def _listing(folder: str) -> set[str]:
    """the names in a folder (MATLAB's dir: results are moved into place under their final
    names, so a listed name is a finished result)"""
    try:
        return set(os.listdir(folder))
    except OSError:
        return set()


def _in_listing(path: str) -> bool:
    from petakit5d.core.paths import fileparts

    folder, name, ext = fileparts(path)
    return (name + ext) in _listing(folder or ".")


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _touch(path: str) -> None:
    with open(path, "w"):
        pass


def _delete(path: str) -> None:
    """MATLAB's ``delete``: a file is removed; a folder stays, with MATLAB's warning"""
    if os.path.isdir(path):
        _warn(
            f"'{path}' is a directory.  Use rmdir to delete directories.",
            "MATLAB:DELETE:DirectoryDeletion",
        )
    elif os.path.lexists(path):
        os.remove(path)


def _msize(shape: Sequence[int]) -> list[int]:
    """getImageSize's [y x z] of a (z, y, x) size"""
    s = [int(v) for v in shape]
    if len(s) >= 3:
        return [s[1], s[2], s[0], *s[3:]]
    return s


def _channel_hits(path: str, patterns: Sequence[str]) -> list[bool]:
    """``~cellfun(@isempty, regexpi(path, patterns))``"""
    from petakit5d.pipeline.filenames import matlab_regex

    return [matlab_regex(p, True).search(path) is not None for p in patterns]


def _save_parameters(decon_path: str, pr: dict[str, Any]) -> None:
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat

    write_mat(decon_path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, decon_path + "/parameters.json")


def _write_psfgen(psf: np.ndarray, final: str, uuid_name: str) -> None:
    """writetiff to the temporary name, then movefile onto the final name"""
    from petakit5d.io import write_tiff

    write_tiff(np.ascontiguousarray(psf), uuid_name)
    os.replace(uuid_name, final)


def _copy_into_place(src: str, dst: str, uuid: str) -> None:
    root, ext = os.path.splitext(dst)
    tmp = f"{root}_{uuid}{ext}"
    shutil.copyfile(src, tmp)
    os.replace(tmp, dst)


def _read_single(fn: str) -> np.ndarray:
    from petakit5d.core import dtypes
    from petakit5d.io import read_tiff

    return np.ascontiguousarray(dtypes.single(read_tiff(fn)))


def _in_place_resources(size, kw: dict[str, Any]) -> tuple[float, float]:
    from petakit5d import decon

    try:
        fn = decon.in_place_controller_resources
    except AttributeError:
        raise NotImplementedError(
            "the in-place large-file deconvolution (petakit5d.decon.decon_large_in_place, task "
            "R07) is not ported yet"
        ) from None
    return fn(size, kw["batch_size"], kw["block_size"], **kw["opts"])


# ---------------------------------------------------------------------------------------------
# the wrapper
# ---------------------------------------------------------------------------------------------


@register(
    "XR_decon_data_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folder(s)."),
    Param("resultDirName", "@(x) ischar(x)", "Result folder under each dataset."),
    Param("overwrite", "@(x) islogical(x)", "Remove the result folders first (overwrite(1))."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to deconvolve."),
    Param("skewAngle", "@isscalar", "Skew angle (degrees)."),
    Param("dz", "@isscalar", "Scan step (um)."),
    Param("xyPixelSize", "@(x) isvector(x) && numel(x) <= 2", "Pixel size (um).", convert="same"),
    Param("save16bit", "@(x) numel(x) == 1 && islogical(x)", "uint16 results.", convert="same"),
    Param("parseSettingFile", "@islogical", "Flip z stacks after their setting files."),
    Param("flipZstack", "@islogical", "Flip the z stacks."),
    Param("background", "@isnumeric", "Background, or one per channel ([]: 100).", convert="same"),
    Param("dzPSF", "@isnumeric", "z step of the PSFs (um)."),
    Param("edgeErosion", "@isnumeric", "Erosion of the edges after the decon (voxels)."),
    Param("erodeByFTP", "@islogical", "Erode with the mask of the first time point."),
    Param("psfFullpaths", "@iscell", "PSF files, in the order of channel_patterns."),
    Param("deconIter", "@isnumeric", "Iterations, or one per dataset."),
    Param("RLMethod", "@ischar", "'simplified', 'omw' or 'original'."),
    Param("wienerAlpha", "@isnumeric", "Wiener alpha (OMW), or one per channel."),
    Param("OTFCumThresh", "@isnumeric", "OTF cumulative sum threshold (OMW), or one per channel."),
    Param("hannWinBounds", "@isnumeric", "Apodization range of the distance matrix (OMW)."),
    Param("skewed", "@(x) isempty(x) || islogical(x)", "PSF in skewed space (None: automatic)."),
    Param("debug", "@islogical", "Save the iterations' estimates."),
    Param("saveStep", "@isnumeric", "Save an estimate every save_step iterations (debug)."),
    Param("psfGen", "@islogical", "Generate clean PSFs (psfgen/)."),
    Param("GPUJob", "@islogical", "Deconvolve on GPUs."),
    Param("deconRotate", "@islogical", "Deskew and rotate the results (DSR/)."),
    Param("batchSize", "@isvector", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of zarr results.", convert="size_yxz"),
    Param("largeFile", "@islogical", "Deconvolve the files as large files."),
    Param("largeMethod", "@ischar", "Large-file method: 'inmemory' or 'inplace'."),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("saveZarr", "@islogical", "Save the results as zarr."),
    Param("dampFactor", "@isnumeric", "Damp factor of the results."),
    Param("scaleFactor", "@isnumeric", "Scale factor, or one per channel ([]: 1)."),
    Param("deconOffset", "@isnumeric", "Offset of the results.", convert="same"),
    Param("maskFullpaths", "@iscell", "2D masks (xy, xz, yz) of the regions to deconvolve."),
    Param("parseCluster", "@islogical", "Deconvolve the files in Slurm jobs."),
    Param("parseParfor", "@islogical", "Use a local pool (large files)."),
    Param("masterCompute", "@islogical", "The in-place controllers deconvolve blocks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a file's job."),
    Param("taskBatchNum", "@isnumeric", "Files per job (None: automatic for GPU jobs)."),
    Param("uuid", "@ischar", "Id of the temporary names (empty: a new one)."),
    Param("unitWaitTime", "@isnumeric", "Minutes a local .tmp marker counts as another's."),
    Param("maxTrialNum", "@isnumeric", "Tries per file."),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; passed on, ignored)."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param("GPUConfigFile", "@ischar", "Cluster configuration file of the GPU jobs."),
    keep_unmatched=True,
    task="R06",
)
def decon_dataset(
    data_paths,
    *,
    result_dir_name="matlab_decon",
    overwrite=False,
    channel_patterns=("CamA_ch0", "CamA_ch1", "CamB_ch0"),
    skew_angle=32.45,
    dz=0.5,
    xy_pixel_size=0.108,
    save_16bit=True,
    parse_setting_file=False,
    flip_zstack=False,
    background=None,
    dz_psf=0.1,
    edge_erosion=0,
    erode_by_ftp=True,
    psf_fullpaths=("", "", ""),
    decon_iter=15,
    rl_method="simplified",
    wiener_alpha=0.005,
    otf_cum_thresh=0.9,
    hann_win_bounds=(0.8, 1.0),
    skewed=None,
    debug=False,
    save_step=5,
    psf_gen=True,
    gpu_job=False,
    decon_rotate=False,
    batch_size=(1024, 1024, 1024),
    block_size=(256, 256, 256),
    large_file=False,
    large_method="inmemory",
    zarr_file=False,
    save_zarr=False,
    damp_factor=1,
    scale_factor=None,
    decon_offset=0,
    mask_fullpaths=(),
    parse_cluster=True,
    parse_parfor=False,
    master_compute=True,
    job_log_dir="../job_logs",
    cpus_per_task=2,
    task_batch_num=None,
    uuid="",
    unit_wait_time=1,
    max_trial_num=3,
    mcc_mode=False,
    config_file="",
    gpu_config_file="",
) -> None:
    """Deconvolve the datasets' files (XR_decon_data_wrapper): see the module. Sizes in
    (z, y, x)."""
    values = dict(locals())
    from petakit5d.core.ids import get_uuid
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs.publish import remove_path
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain
    from petakit5d.pipeline.filenames import parse_dataset_filenames

    pr = {k: to_plain(v) for k, v in matlab_results(decon_dataset, values).items()}
    if not uuid:
        uuid = get_uuid()
    paths = _as_list(data_paths)
    channels = [str(c) for c in _as_list(channel_patterns)]
    psfs = [str(p) for p in _as_list(psf_fullpaths)]
    masks = [str(m) for m in _as_list(mask_fullpaths)]
    nd = len(paths)
    nch = len(channels)
    sep = os.sep
    last_path = paths[-1] if paths else ""
    paths = [p if p.endswith(sep) else p + sep for p in paths]

    iters = _vector(decon_iter)
    iters = iters if iters.size == nd else np.full(nd, iters[0] if iters.size else np.nan)
    ow = [bool(v) for v in np.atleast_1d(np.asarray(overwrite)).reshape(-1)]
    ow = ow * 2 if len(ow) == 1 else ow
    alphas = _per_channel(wiener_alpha, nch)
    threshs = _per_channel(otf_cum_thresh, nch)
    bgs = _per_channel(100 if background is None or np.size(background) == 0 else background, nch)
    scales = _per_channel(
        1 if scale_factor is None or np.size(scale_factor) == 0 else scale_factor, nch
    )
    method = str(rl_method)

    job_log_fname = job_error_fname = ""
    if parse_cluster:
        from petakit5d.jobs.slurm import check_slurm_cluster

        check = check_slurm_cluster(last_path, job_log_dir)
        parse_cluster = bool(check.parse_cluster)
        job_log_fname, job_error_fname = check.job_log_fname, check.job_log_error_fname

    decon_name = result_dir_name or "matlab_decon"
    decon_paths = []
    for d in range(nd):
        decon_path = f"{paths[d]}{decon_name}{sep}"
        if ow[0] and os.path.isdir(decon_path):
            remove_path(decon_path)
        if not os.path.isdir(decon_path):
            os.makedirs(decon_path, exist_ok=True)
            group_writable(decon_path)
        decon_paths.append(decon_path)
        _save_parameters(decon_path, pr)

    _generate_psfs(psfs, paths, decon_name, method, dz_psf, dz, psf_gen, uuid)

    if large_file and large_method == "inplace":
        save_zarr = True

    parsed = parse_dataset_filenames(
        paths, channels, parse_setting_file, flip_zstack, True, decon_paths, False, 1, zarr_file
    )
    fnames = list(parsed.fnames)
    fd_inds = [int(v) for v in parsed.fd_inds]
    flips = [bool(v) for v in parsed.flip_zstack]
    ftp_inds = [int(v) for v in parsed.ftp_inds]
    mask_paths = list(parsed.mask_fullpaths)
    n_files = len(fnames)

    if parse_cluster and method == "omw" and not large_file:
        _premake_backprojectors(
            psfs, channels, paths, fnames, fd_inds, decon_name, method, alphas, threshs,
            skewed, hann_win_bounds, psf_gen, uuid,
        )  # fmt: skip

    loop = _Loop(
        n_files=n_files,
        max_trial_num=max_trial_num,
        parse_cluster=parse_cluster,
    )
    static = {
        "xy_pixel_size": xy_pixel_size,
        "dz": dz,
        "dz_psf": dz_psf,
        "skew_angle": skew_angle,
        "edge_erosion": edge_erosion,
        "rotate": decon_rotate,
        "rl_method": method,
        "hann_win_bounds": _hann_text(hann_win_bounds),
        "skewed": skewed,
        "debug": debug,
        "save_step": save_step,
        "psf_gen": psf_gen,
        "save_zarr": save_zarr,
        "parse_cluster": parse_cluster,
        "parse_parfor": parse_parfor,
        "gpu_job": gpu_job,
        "save_16bit": save_16bit,
        "large_file": large_file,
        "large_method": large_method,
        "batch_size": _mat2str_commas(_size_matlab(batch_size)),
        "block_size": _mat2str_commas(_size_matlab(block_size)),
        "damp_factor": damp_factor,
        "decon_offset": decon_offset,
        "mask_fullpaths": "{'" + "','".join(masks) + "'}",
        "uuid": uuid,
        "cpus_per_task": cpus_per_task,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
        "gpu_config_file": gpu_config_file,
    }
    settings = {
        "paths": paths,
        "decon_paths": decon_paths,
        "fnames": fnames,
        "fd_inds": fd_inds,
        "flips": flips,
        "ftp_inds": ftp_inds,
        "mask_paths": mask_paths,
        "channels": channels,
        "psfs": psfs,
        "iters": iters,
        "alphas": alphas,
        "threshs": threshs,
        "bgs": bgs,
        "scales": scales,
        "method": method,
        "static": static,
        "edge_erosion": float(_vector(edge_erosion)[0]),
        "erode_by_ftp": bool(erode_by_ftp),
        "large_file": bool(large_file),
        "in_place": bool(large_file) and large_method == "inplace",
        "save_zarr": bool(save_zarr),
        "zarr_file": bool(zarr_file),
        "gpu_job": bool(gpu_job),
        "master_compute": bool(master_compute),
        "job_log_dir": job_log_dir,
        "job_log_fname": job_log_fname,
        "job_error_fname": job_error_fname,
        "cpus_per_task": cpus_per_task,
        "config_file": config_file,
        "gpu_config_file": gpu_config_file,
        "task_batch_num": task_batch_num,
        "unit_wait_time": float(unit_wait_time),
        "uuid": uuid,
        "mcc_mode": mcc_mode,
        "save_16bit": bool(save_16bit),
        "batch_size": batch_size,
        "block_size": block_size,
    }
    loop.run(settings)


def _size_matlab(v: Any) -> np.ndarray:
    """a (z, y, x) size (or a scalar) as MATLAB's [y x z]"""
    a = _vector(v)
    if a.size == 3:
        return a[[1, 2, 0]]
    return a


def _hann_text(v: Any) -> str:
    from petakit5d.core.text import mat2str_comma

    return mat2str_comma(np.asarray(v, dtype=np.float64))


def _generate_psfs(psfs, paths, decon_name, method, dz_psf, dz, psf_gen, uuid) -> None:
    """the generated PSFs of each dataset (psfgen/<psf>_<RLMethod>.tif)"""
    from petakit5d.core.paths import fileparts
    from petakit5d.decon.psf import crop_psf_support
    from petakit5d.decon.psf import psf_gen as gen

    for fn in psfs:
        if not _exist(fn):
            raise DeconDatasetError(f"PSF file {fn} does not exist!")
        if not psf_gen:
            continue
        print(f"PSF generation for {fn} ...", flush=True)
        fsn = fileparts(fn)[1]
        psf = gen(_read_single(fn), dz_psf, dz, 1.5, "masked")
        if method != "omw":
            psf = crop_psf_support(psf)
        for data_path in paths:
            final = f"{data_path}/{decon_name}/psfgen/{fsn}_{method}.tif"
            _write_psfgen(psf, final, f"{data_path}/{decon_name}/psfgen/{fsn}_{uuid}.tif")
        print("Done!\n", flush=True)


def _premake_backprojectors(
    psfs, channels, paths, fnames, fd_inds, decon_name, method, alphas, threshs, skewed,
    hann_win_bounds, psf_gen, uuid,
) -> None:  # fmt: skip
    """the OMW back projector of each PSF, made once before the jobs (as the jobs get the
    parameters from the function string), copied to the other datasets"""
    from petakit5d.core.paths import fileparts
    from petakit5d.core.text import mat2str_comma, sprintf, through_text
    from petakit5d.decon.backprojector import omw_backprojector_files
    from petakit5d.io import get_image_size
    from petakit5d.params.literal import str2num
    from petakit5d.params.values import to_python

    hann = to_python(str2num(mat2str_comma(np.asarray(hann_win_bounds, dtype=np.float64)))[0])
    hann = list(hann) if isinstance(hann, tuple) else [hann]
    for p in range(min(len(psfs), len(channels))):
        psf_fn = psfs[p]
        fsn = fileparts(psf_fn)[1]
        alpha = float(through_text(float(alphas[p]), "%.20f"))
        thresh = float(through_text(float(threshs[p]), "%.20f"))
        names = [
            sprintf("%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f.tif", fsn, alpha, thresh),
            f"{fsn}_OTF.tif",
            sprintf("%s_OTF_mask_otf_cum_thresh_%0.6f.tif", fsn, thresh),
        ]
        src = ""
        for d, data_path in enumerate(paths):
            full = [data_path + fnames[i] for i in range(len(fnames)) if fd_inds[i] == d]
            hits = [_channel_hits(f, [channels[p]])[0] for f in full]
            if not any(hits):
                continue
            first = full[hits.index(True)]
            psfgen_path = f"{data_path}/{decon_name}/psfgen/"
            psf = _read_single(f"{psfgen_path}/{fsn}_{method}.tif" if psf_gen else psf_fn)
            # RLdecon crops a PSF larger than the data first
            if any(
                a > b for a, b in zip(_psf_size3(psf), _size3(get_image_size(first)), strict=False)
            ):
                continue
            os.makedirs(psfgen_path, exist_ok=True)
            if not src:
                # without the figure, as in the cluster jobs
                omw_backprojector_files(
                    psf, psfgen_path, fsn, alpha, thresh, skewed, hann, psf_fn, False
                )
                src = psfgen_path
            elif not os.path.exists(psfgen_path + names[0]):
                # the same PSF and parameters as in the first dataset
                for name in names:
                    if os.path.exists(src + name):
                        _copy_into_place(src + name, psfgen_path + name, uuid)


def _psf_size3(psf: np.ndarray) -> tuple[int, int, int]:
    """size(psf, 1 : 3) as (z, y, x)"""
    s = psf.shape
    return (1, s[0], s[1]) if psf.ndim == 2 else (s[0], s[1], s[2])


def _size3(size: Sequence[int]) -> tuple[int, ...]:
    s = tuple(int(v) for v in size)
    return s[:3]


class _Loop:
    """XR_decon_data_wrapper's rounds over the files"""

    def __init__(self, *, n_files: int, max_trial_num: float, parse_cluster: bool) -> None:
        self.n = n_files
        self.max_trial = float(max_trial_num)
        self.cluster = parse_cluster
        self.done = [False] * n_files
        self.trials = [0] * n_files
        self.job_ids = [-1] * n_files
        self.job_task_ids = [0] * n_files
        self.active = [False] * n_files
        self.other_writer = [False] * n_files

    def settled(self) -> bool:
        return all(
            self.done[i] or (self.trials[i] >= self.max_trial and not self.active[i])
            for i in range(self.n)
        )

    def run(self, s: dict[str, Any]) -> None:
        from petakit5d.jobs import single

        mark = single.mark() if self.cluster else None
        try:
            self._rounds(s)
        finally:
            if mark is not None:
                single.end(mark)

    # -- one round ---------------------------------------------------------------------------
    def _rounds(self, s: dict[str, Any]) -> None:
        from petakit5d.core.paths import fileparts
        from petakit5d.jobs import files as F

        n_done = 0
        t0 = time.perf_counter()
        in_place = s["in_place"]
        nice = 0
        controller = None
        if self.cluster and in_place:
            from petakit5d.jobs.coordinator import coordinator_nice

            # (the controllers wait for their block jobs: these go first)
            nice = coordinator_nice()
            controller = {
                "batch_size": s["batch_size"],
                "block_size": s["block_size"],
                "opts": {
                    "save_16bit": s["save_16bit"],
                    "gpu_job": s["gpu_job"],
                    "master_compute": s["master_compute"],
                    "config_file": s["config_file"],
                    "cpus_per_task": s["cpus_per_task"],
                },
            }
        auto_batch = s["task_batch_num"] is None or np.size(s["task_batch_num"]) == 0
        batch_num = 1
        if not auto_batch:
            batch_num = max(1, int(_mround(float(_vector(s["task_batch_num"])[0]))))
        batch_files = (
            self.cluster
            and not s["large_file"]
            and ((auto_batch and s["gpu_job"]) or (not auto_batch and batch_num > 1))
        )
        out_ext = ".zarr" if s["save_zarr"] else ".tif"

        while not self.settled():
            cands: list[dict[str, Any]] = []
            listings = [_listing(p) for p in s["decon_paths"]] if self.cluster else []
            awaited: list[str] = []
            for f in range(self.n):
                if self.done[f]:
                    continue
                fname = s["fnames"][f]
                fsname = fileparts(fname)[1]
                fd = s["fd_inds"][f]
                frame = s["paths"][fd] + fname
                if not _exist(frame):
                    self.done[f] = True
                    continue
                task_id = (f + 1) % ARRAY_ID_MODULUS
                decon_path = s["decon_paths"][fd]
                result = f"{decon_path}/{fsname}{out_ext}"
                tmp_marker = result[:-4] + ".tmp"
                rname = fsname + out_ext
                if _exist(result) or (self.cluster and rname in listings[fd]):
                    self.done[f] = True
                    if os.path.isfile(tmp_marker):
                        os.remove(tmp_marker)

                erode_mask = ""
                save_mask = False
                if s["edge_erosion"] > 0 and s["erode_by_ftp"] and not s["large_file"]:
                    mask_fn = s["mask_paths"][fd] or ""
                    if f == s["ftp_inds"][fd]:
                        save_mask = True
                        # a result without the mask is computed again, to save the mask
                        if (
                            self.done[f]
                            and not os.path.isdir(mask_fn)
                            and not (self.cluster and _in_listing(mask_fn))
                        ):
                            self.done[f] = False
                            _delete(result)
                    elif not self.done[f]:
                        erode_mask = mask_fn
                        if not F.batch_file_exist([erode_mask])[0]:
                            awaited.append(erode_mask)
                            continue
                        from petakit5d.io import get_image_size

                        mask_sz = _msize(get_image_size(erode_mask))
                        img_sz = _msize(get_image_size(frame))
                        if mask_sz != img_sz:
                            from petakit5d.core.text import num2str

                            _warn(
                                f"The image size [{num2str(img_sz, '%d ')}] does not match the "
                                f"defined mask size [{num2str(mask_sz, '%d ')}], use its own mask "
                                "for edge erosion..."
                            )
                            erode_mask = ""
                if self.done[f]:
                    continue

                func_str = self._func_str(s, f, frame, decon_path, erode_mask, save_mask)
                if os.path.exists(tmp_marker) or self.cluster:
                    if self.cluster:
                        if in_place:
                            from petakit5d.jobs import writers

                            # the volume's controller of another call that may still make it:
                            # wait for it, without a try
                            live = writers.live(result, s["uuid"])
                            if live:
                                if not self.other_writer[f]:
                                    print(
                                        f"{result} is made by a job of another call ({live[0]}), "
                                        "wait for it...",
                                        flush=True,
                                    )
                                    self.other_writer[f] = True
                                self.active[f] = True
                                continue
                        cands.append(
                            self._candidate(s, f, frame, func_str, task_id, rname, controller)
                        )
                        continue
                    age = (time.time() - os.path.getmtime(tmp_marker)) / 60
                    if age < s["unit_wait_time"]:
                        continue
                    _touch(tmp_marker)
                else:
                    _touch(tmp_marker)

                # local mode: the file here, from its function string
                from petakit5d.jobs.tasks import run_record, task_record

                t1 = time.perf_counter()
                run_record(task_record(frame_task(func_str), "python"))
                print(f"Elapsed time is {time.perf_counter() - t1:f} seconds.\n", flush=True)
                self.trials[f] += 1
                if _exist(result):
                    self.done[f] = True
                    if os.path.isfile(tmp_marker):
                        os.remove(tmp_marker)

            if self.cluster and cands:
                self._submit(s, cands, batch_files, auto_batch, batch_num, nice)

            if not self.settled():
                left = sum(not x for x in self.done)
                wait = 30.0
                if self.cluster:
                    wait = min(30.0, max(3.0, 0.3 * left))
                coming = [
                    i
                    for i in range(self.n)
                    if not self.done[i]
                    and not (self.trials[i] >= self.max_trial and not self.active[i])
                ]
                results = [
                    f"{s['decon_paths'][s['fd_inds'][i]]}/{fileparts(s['fnames'][i])[1]}{out_ext}"
                    for i in coming
                ]
                groups = [results] + [[m] for m in sorted(set(awaited))]
                F.wait_for_files(groups, wait)
            if n_done < sum(self.done):
                n_done = sum(self.done)
                print(
                    f"Time {time.perf_counter() - t0:0.2f} s: {n_done} / {self.n} "
                    f"({n_done / self.n * 100:0.2f}%) are finished!",
                    flush=True,
                )
        if all(self.done):
            print(f"All output files ({self.n} / {self.n}) are finished!\n", flush=True)

    def _func_str(self, s, f, frame, decon_path, erode_mask, save_mask) -> str:
        fd = s["fd_inds"][f]
        hits = _channel_hits(frame, s["channels"])
        psf_idx = [i for i, h in enumerate(hits) if h]
        if not psf_idx:
            raise DeconDatasetError(
                "Insufficient number of outputs from right hand side of equal sign to satisfy "
                "assignment.",
                "MATLAB:needMoreRhsOutputs",
            )
        if psf_idx[0] >= len(s["psfs"]):
            raise DeconDatasetError(
                f"Index exceeds the number of array elements. Index must not exceed "
                f"{len(s['psfs'])}.",
                "MATLAB:badsubscript",
            )
        v = dict(s["static"])
        v["frame"] = frame
        v["decon_path"] = decon_path
        v["psf"] = s["psfs"][psf_idx[0]]
        v["background"] = _pick(s["bgs"], hits, "background")
        it = float(s["iters"][fd])
        v["decon_iter"] = int(it) if it.is_integer() else it
        if s["method"] == "omw":
            v["wiener_alpha"] = _pick(s["alphas"], hits, "Wiener alpha")
            v["otf_cum_thresh"] = _pick(s["threshs"], hits, "OTF threshold")
        else:
            v["wiener_alpha"] = float(s["alphas"][0])
            v["otf_cum_thresh"] = float(s["threshs"][0])
        v["scale_factor"] = _pick(s["scales"], hits, "scale factor")
        v["flip_zstack"] = s["flips"][f]
        v["erode_mask"] = erode_mask
        v["save_maskfile"] = save_mask
        func_str = frame_task_string(v)
        if s["in_place"] and not s["master_compute"]:
            # (only then, as the controller's default is true)
            func_str = func_str[:-1] + ",'masterCompute',false)"
        if s["in_place"]:
            # (and the job log folder, so that the controller's block jobs log where the
            # caller's jobs do)
            func_str = func_str[:-1] + f",'jobLogDir','{s['job_log_dir']}')"
        return func_str

    def _candidate(self, s, f, frame, func_str, task_id, rname, controller) -> dict[str, Any]:
        from petakit5d.io import get_image_size

        size = tuple(int(v) for v in get_image_size(frame))
        est = float(np.prod(np.asarray(size, dtype=np.float64))) * 4 / 2**30
        mem = est * 10
        cfg = s["config_file"]
        cpus = s["cpus_per_task"]
        if s["gpu_job"] and not s["in_place"]:
            cfg = s["gpu_config_file"]
            # the GPU deconvolves and the CPUs only read and write the file: at most 4 CPUs,
            # so that every GPU of a node can work
            if not s["large_file"] and s["gpu_config_file"]:
                cpus = min(cpus, GPU_JOB_CPUS)
        if controller is not None:
            cpus, mem = _in_place_resources(size, controller)
        fd = s["fd_inds"][f]
        return {
            "f": f,
            "frame": frame,
            "func_str": func_str,
            "task_id": task_id,
            "size": size,
            "out_name": rname,
            "decon_path": s["decon_paths"][fd],
            "fd": fd,
            "mem": mem,
            "cpus": cpus,
            "config_file": cfg,
        }

    # -- cluster jobs ------------------------------------------------------------------------
    def _job_params(self, s: dict[str, Any], c: dict[str, Any]) -> dict[str, Any]:
        return {
            "job_log_dir": s["job_log_dir"],
            "job_log_fname": s["job_log_fname"],
            "job_error_fname": s["job_error_fname"],
            "last_file": False,
            "cpus_per_task": c["cpus"],
            "mem_allocate": c["mem"],
            "mcc_mode": s["mcc_mode"],
            "config_file": c["config_file"],
        }

    def _submit(self, s, cands, batch_files, auto_batch, batch_num, nice) -> None:
        from petakit5d.jobs import single

        # one squeue for the jobs of the round's files
        keyed = sorted(
            {
                (self.job_ids[c["f"]], self.job_task_ids[c["f"]])
                for c in cands
                if self.job_ids[c["f"]] > 0
            }
        )
        states = (
            dict(
                zip(
                    keyed,
                    single.job_states(
                        [k[0] for k in keyed], [k[1] % ARRAY_ID_MODULUS for k in keyed]
                    ),
                    strict=True,
                )
            )
            if keyed
            else {}
        )

        if not batch_files:
            for c in cands:
                f = c["f"]
                known = states.get((self.job_ids[f], self.job_task_ids[f]))
                submit_status = 0
                if self.trials[f] < self.max_trial:
                    job_id, job_status, submit_status = single.submit_single_job(
                        frame_task(c["func_str"]),
                        self.job_ids[f],
                        c["task_id"],
                        nice=nice,
                        known_status=known,
                        output=f"{c['decon_path']}/{c['out_name']}",
                        **self._job_params(s, c),
                    )
                else:
                    # the last try: only whether its job (or a task runner of another job that
                    # took its task) is still there
                    job_id = self.job_ids[f]
                    job_status = single.check_single_job_status(
                        job_id, c["task_id"], known_status=known
                    )
                self.job_ids[f] = job_id
                self.job_task_ids[f] = c["task_id"]
                self.trials[f] += submit_status
                self.active[f] = submit_status > 0 or job_status >= 0
            return

        # files of several per job: those whose job is pending or running wait for it
        jobs = sorted({(self.job_ids[c["f"]], self.job_task_ids[c["f"]]) for c in cands})
        live = {}
        for j in jobs:
            live[j] = (
                j[0] > 0
                and single.check_single_job_status(j[0], j[1], known_status=states.get(j)) >= 0
            )
        subs = []
        for c in cands:
            f = c["f"]
            lv = live[(self.job_ids[f], self.job_task_ids[f])]
            self.active[f] = lv
            if not lv:
                subs.append(c)
        # a job may have finished since the listing at the start of the round
        if any(self.trials[c["f"]] > 0 for c in subs):
            listings = [_listing(p) for p in s["decon_paths"]]
            keep = []
            for c in subs:
                if self.trials[c["f"]] > 0 and c["out_name"] in listings[c["fd"]]:
                    self.done[c["f"]] = True
                else:
                    keep.append(c)
            subs = keep
        # files whose jobs ended max_trial_num times are not submitted again
        subs = [c for c in subs if self.trials[c["f"]] < self.max_trial]
        # files submitted before (their job failed) get a job each
        retry = [self.trials[c["f"]] > 0 for c in subs]
        if auto_batch:
            batch_num = max(1, min(AUTO_BATCH_MAX, sum(not r for r in retry) // AUTO_BATCH_JOBS))
        # a job takes files of the same settings (the function string after the file) and size
        prefix = [f"{FRAME_NAME}('{c['frame']}'" for c in subs]
        param_strs = [c["func_str"][len(p) :] for c, p in zip(subs, prefix, strict=True)]
        keys = [f"{ps}|{c['size']}" for ps, c in zip(param_strs, subs, strict=True)]
        groups: dict[str, int] = {}
        param_inds = []
        for k in keys:
            groups.setdefault(k, len(groups))
            param_inds.append(groups[k])
        n_groups = len(groups)
        for i, r in enumerate(retry):
            if r:
                param_inds[i] = n_groups
                n_groups += 1
        for p in range(n_groups):
            pos = [i for i, q in enumerate(param_inds) if q == p]
            for b in range(math.ceil(len(pos) / batch_num)):
                bpos = pos[b * batch_num : (b + 1) * batch_num]
                bfs = [subs[i] for i in bpos]
                func_str = bfs[0]["func_str"]
                if len(bfs) > 1:
                    frames = ",".join(f"'{c['frame']}'" for c in bfs)
                    func_str = f"{FRAME_NAME}({{{frames}}}{param_strs[bpos[0]]}"
                first = bfs[0]
                # (a file submitted again runs as before, without task runners)
                job_id, _, submit_status = single.submit_single_job(
                    frame_task(func_str),
                    -1,
                    first["task_id"],
                    task_runner=not retry[bpos[0]],
                    **{**self._job_params(s, first), "mem_allocate": max(c["mem"] for c in bfs)},
                )
                for c in bfs:
                    self.job_ids[c["f"]] = job_id
                    self.job_task_ids[c["f"]] = first["task_id"]
                    self.trials[c["f"]] += submit_status
                    self.active[c["f"]] = submit_status > 0


def _mround(x: float) -> float:
    from petakit5d.core.dtypes import mround

    return float(mround(x))
