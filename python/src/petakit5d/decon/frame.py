"""Deconvolution of the files of a job (``XR_RLdeconFrame3D``): each file in one piece
(``decon.file.rl_decon_file``) or, when it is large, with the large-file functions (task R07).

``decon_file`` is XR_RLdeconFrame3D.m at the oracle commit:

- **Files and results**: each file of ``frame_fullpaths`` (missing ones are skipped with a
  warning) goes to ``<decon_path>/<name>.tif`` (``.zarr`` with ``save_zarr``), ``decon_path``
  defaulting to ``<file folder>/matlab_decon``; DSR results to ``<decon_path>/DSR/``. An existing
  result is kept (the file is skipped), or with ``overwrite`` removed with its DSR result first
  (the decon skips existing results).
- **Large files**: a file needs ``10 x 4`` bytes per voxel; above ``cpu_max_mem`` GB it, and every
  later file of the call, goes to the large-file functions (``large_file``), as does a file whose
  one-piece decon wrote neither result. ``large_method`` 'inplace' or 'inmemory' then returns
  from the call, as MATLAB's; another method does nothing. These functions are task R07's
  (``petakit5d.decon.decon_large_in_place`` and ``decon_large_in_memory``): until it lands they
  raise ``NotImplementedError``.
- **First-time-point masks**: with ``edge_erosion`` > 0 and ``save_maskfile``, the mask of the
  file's nonzero voxels (flipped with ``flip_zstack``) eroded by ``edge_erosion`` is saved as
  ``<decon_path>/Masks/<name>_eroded.zarr`` (uint8, replacing an older one) and used as the eroded
  mask file of the decon (and of the files after it).
- **GPU jobs** (``gpu_job``): the decon runs with ``use_gpu`` True; a file that doesn't fit on the
  GPU (``GPUOutOfMemoryError`` or another error with one of MATLAB's GPU capacity identifiers) is
  deconvolved on the CPU instead, other errors are raised. With several files and room on the
  GPU for the decon twice over (MATLAB's 42 bytes per voxel, 38 for OMW: R04's peak for uint16
  data), the next file is read while the GPU deconvolves (``overlap_fcn``), if it will be
  deconvolved in one piece, and handed to its own decon.

The background defaults to 100, the method name is lowercased, and the decon gets
``DSRCombined`` true and the defaults of RLdecon for what this function doesn't pass on.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
import warnings
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = ["decon_file"]

# MATLAB's GPU capacity errors (memory of an array, of an operation, of an FFT plan; the gpuArray
# size limit), which XR_RLdeconFrame3D catches to deconvolve the file on the CPU instead
GPU_CAPACITY_IDS = (
    "parallel:gpu:array:OOM",
    "parallel:gpu:array:OOMForOperation",
    "parallel:gpu:array:FFTPlanInternalError",
    "parallel:gpu:array:pmaxsize",
)


class DeconFileError(ValueError):
    """An error MATLAB raises, with its identifier ('' where MATLAB's error has none)."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _warn(message: str, identifier: str = "") -> None:
    w = UserWarning(message)
    w.identifier = identifier
    warnings.warn(w, stacklevel=3)


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _slurm_available() -> bool:
    """MATLAB's ``system('sinfo')`` status is 0"""
    exe = shutil.which("sinfo")
    if exe is None:
        return False
    try:
        r = subprocess.run([exe], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return False
    return r.returncode == 0


def _mem_required(size) -> float:
    """``prod(imSize) * 4 / 2^30 * 10`` (GB)"""
    return float(np.prod(np.asarray(size, dtype=np.float64))) * 4 / 2**30 * 10


def _result_path(decon_path: str, fsname: str, save_zarr: bool) -> str:
    return f"{decon_path}/{fsname}{'.zarr' if save_zarr else '.tif'}"


def is_gpu_capacity_error(err: BaseException) -> bool:
    """the GPU ran out of memory (an array, an operation or an FFT plan) or the data is over the
    gpuArray size limit: the errors a GPU job deconvolves on the CPU instead"""
    try:
        from petakit5d.kernels.decon_gpu import GPUOutOfMemoryError
    except ImportError:  # pragma: no cover (the module is pure Python)
        GPUOutOfMemoryError = ()  # noqa: N806
    if GPUOutOfMemoryError and isinstance(err, GPUOutOfMemoryError):
        return True
    return getattr(err, "identifier", None) in GPU_CAPACITY_IDS


def gpu_free_bytes() -> int | None:
    """free memory of the GPU the job would use (MATLAB's ``gpuDevice().AvailableMemory``), None
    without a usable GPU"""
    try:
        from petakit5d.gpu.device import select_device
        from petakit5d.gpu.memory import available_bytes

        select_device()
        return int(available_bytes())
    except Exception:  # noqa: BLE001 (MATLAB's try/catch around gpuDevice)
        return None


def gpu_has_room(im_size, method: str) -> bool:
    """whether the GPU has the memory of a decon of im_size free twice over (MATLAB's loops on the
    GPU take about 38 (OMW) or 42 bytes per voxel; R04's peaks for uint16 data are the same)"""
    per_voxel = 38 if method == "omw" else 42
    free = gpu_free_bytes()
    if free is None:
        return False
    return free >= 2 * per_voxel * float(np.prod(np.asarray(im_size, dtype=np.float64)))


def next_frame_needs_decon(
    fn: str, decon_path: str, save_zarr: bool, overwrite: bool, cpu_max_mem: float
) -> bool:
    """whether the next file of the job will be deconvolved in one piece, as the loop decides
    (it exists, has no result or overwrite, and isn't large)"""
    from petakit5d.core.paths import fileparts
    from petakit5d.io import get_image_size

    if not _exist(fn):
        return False
    pathstr, fsname, _ = fileparts(fn)
    if not decon_path:
        decon_path = f"{pathstr}/matlab_decon"
    if not overwrite and _exist(_result_path(decon_path, fsname, save_zarr)):
        return False
    return _mem_required(get_image_size(fn)) <= cpu_max_mem


def read_frame(fn: str) -> np.ndarray | None:
    """the raw data of a file as RLdecon reads it (None if it can't be read: RLdecon then reads it
    itself)"""
    from petakit5d import io

    try:
        ext = os.path.splitext(fn.rstrip("/\\"))[1]
        if ext in (".tif", ".tiff"):
            return io.read_tiff(fn)
        if ext == ".zarr":
            return io.read_zarr(fn)
    except Exception as e:  # noqa: BLE001 (MATLAB's catch: disp(ME))
        print(e, flush=True)
    return None


def _remove(path: str) -> None:
    from petakit5d.jobs.publish import remove_path

    if os.path.lexists(path.rstrip("/\\") or path):
        remove_path(path)


def _save_mask(fn: str, ext: str, decon_path: str, fsname: str, kw: dict) -> str:
    """the eroded mask of the file's nonzero voxels as <decon_path>/Masks/<name>_eroded.zarr"""
    from petakit5d import io
    from petakit5d.decon.erosion import decon_mask_edge_erosion
    from petakit5d.decon.file import _flip_z, _read
    from petakit5d.jobs.publish import remove_path

    print("Create eroded masks using raw data...", flush=True)
    im = _read(fn, None) if ext in (".tif", ".tiff", ".zarr") else None
    if im is None:
        raise DeconFileError(_UNDEFINED, "Unrecognized function or variable 'im_raw'.")
    if kw["flip_zstack"]:
        im = _flip_z(im, True)
    bw = decon_mask_edge_erosion(im > 0, kw["edge_erosion"])
    del im
    mask_path = f"{decon_path}/Masks"
    os.makedirs(mask_path, exist_ok=True)
    final = f"{mask_path}/{fsname}_eroded.zarr"
    tmp = f"{mask_path}/{fsname}_eroded_{kw['uuid']}.zarr"
    io.write_zarr(
        np.asarray(bw).view(np.uint8) if bw.dtype == np.bool_ else bw.astype(np.uint8),
        tmp,
        block_size=kw["block_size"],
    )
    if os.path.isdir(final):
        remove_path(final)
    os.replace(tmp, final)
    return final


_UNDEFINED = "MATLAB:UndefinedFunction"


def _large(name: str):
    """R07's large-file function (through the decon package)"""
    from petakit5d import decon

    try:
        return getattr(decon, name)
    except AttributeError:
        raise NotImplementedError(
            f"the large-file deconvolution (petakit5d.decon.{name}, task R07) is not ported yet"
        ) from None


@register(
    "XR_RLdeconFrame3D",
    Param("frameFullpaths", "@(x) ischar(x) || iscell(x)", "File(s) to deconvolve."),
    Param("xyPixelSize", "@isnumeric", "Pixel size (um)."),
    Param("dz", "@isnumeric", "Scan step (um)."),
    Param(
        "deconPath",
        "@(x) ischar(x) || isempty(x)",
        "Result folder ('': <file folder>/matlab_decon).",
    ),
    Param("PSFfile", "@ischar", "PSF file.", py="psf_file"),
    Param("Overwrite", "@islogical", "Recompute existing results."),
    Param("save16bit", "@islogical", "Results in uint16 (else single).", convert="same"),
    Param("Rotate", "@islogical", "Rotate the result (with the deskew for sample scans)."),
    Param("Deskew", "@islogical", "Deskew the result."),
    Param("SkewAngle", "@isnumeric", "Skew angle (degrees)."),
    Param("flipZstack", "@islogical", "Flip the data in z first."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("Background", "@isnumeric", "Background ([]: 100).", convert="same"),
    Param("EdgeSoften", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("zEdgeSoften", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("dzPSF", "@isnumeric", "z step of the PSF (um)."),
    Param("DeconIter", "@isnumeric", "Number of iterations."),
    Param("EdgeErosion", "@isnumeric", "Erosion of the edges of the mask of nonzero data."),
    Param("ErodeMaskfile", "@ischar", "Eroded mask zarr used instead of the erosion."),
    Param("SaveMaskfile", "@islogical", "Save the eroded mask (for the other time points)."),
    Param("RLMethod", "@ischar", "'simplified', 'omw', 'original' or 'cudagen' (any case)."),
    Param("wienerAlpha", "@isnumeric", "Wiener alpha of the OMW back projector."),
    Param("OTFCumThresh", "@isnumeric", "OTF cumulative sum threshold (OMW)."),
    Param("hannWinBounds", "@isnumeric", "Apodization range of the distance matrix (OMW)."),
    Param("skewed", "@(x) isempty(x) || islogical(x)", "Decon in skewed space (None: automatic)."),
    Param("fixIter", "@islogical", "Fixed number of iterations (saved in the info file)."),
    Param("errThresh", "@isnumeric", "Error threshold (saved in the info file)."),
    Param("CPUMaxMem", "@isnumeric", "Memory (GB) above which a file is decon'd as a large file."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("blockSize", "@isnumeric", "Block size (z, y, x) of zarr results.", convert="size_yxz"),
    Param("largeFile", "@islogical", "Deconvolve the files as large files."),
    Param("largeMethod", "@ischar", "Large-file method: 'inmemory' or 'inplace'."),
    Param("saveZarr", "@islogical", "Save results as zarr."),
    Param(
        "save3Dstack",
        "@islogical",
        "Save the (decon, DS, DSR) stacks.",
        convert="same",
        py="save_3d_stack",
    ),
    Param("mipAxis", "@isnumeric", "MIP axes (z, y, x).", convert="flags_yxz"),
    Param("dampFactor", "@isnumeric", "Damp factor of the result."),
    Param("scaleFactor", "@isnumeric", "Scale factor of the result ([]: 1)."),
    Param("deconOffset", "@isnumeric", "Offset of the result.", convert="same"),
    Param(
        "maskFullpaths", "@iscell", "2D masks (xy, xz, yz) of the regions to decon (large files)."
    ),
    Param("parseCluster", "@islogical", "Use the Slurm cluster (large files)."),
    Param("parseParfor", "@islogical", "Use parfor (large files)."),
    Param("masterCompute", "@islogical", "The master runs tasks too (large files)."),
    Param("masterCPU", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("GPUJob", "@islogical", "Deconvolve on a GPU (the CPU for files that don't fit)."),
    Param("jobLogDir", "@ischar", "Job log folder (large files)."),
    Param("cpusPerTask", "@isnumeric", "CPUs per task (large files)."),
    Param("uuid", "@ischar", "uuid of the temporary names."),
    Param("maxTrialNum", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("unitWaitTime", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("debug", "@islogical", "Save the iterations' estimates."),
    Param("saveStep", "@isnumeric", "Save an estimate every saveStep iterations (debug)."),
    Param("psfGen", "@islogical", "Generate the PSF (psf_gen) or use it as it is."),
    Param("mccMode", "@islogical", "Compiled mode (large files)."),
    Param("configFile", "@ischar", "Cluster config file (large files)."),
    Param("GPUConfigFile", "@ischar", "GPU cluster config file (large files)."),
    task="R05",
)
def decon_file(
    frame_fullpaths,
    xy_pixel_size,
    dz,
    decon_path="",
    *,
    psf_file="",
    overwrite=False,
    save_16bit=True,
    rotate=False,
    deskew=False,
    skew_angle=-32.45,
    flip_zstack=False,
    reverse=True,
    background=None,
    edge_soften=5,
    z_edge_soften=2,
    dz_psf=0.1,
    decon_iter=15,
    edge_erosion=8,
    erode_maskfile="",
    save_maskfile=False,
    rl_method="simplified",
    wiener_alpha=0.005,
    otf_cum_thresh=0.9,
    hann_win_bounds=(0.8, 1),
    skewed=None,
    fix_iter=False,
    err_thresh=None,
    cpu_max_mem=500,
    batch_size=(1024, 1024, 1024),
    block_size=(256, 256, 256),
    large_file=False,
    large_method="inmemory",
    save_zarr=False,
    save_3d_stack=(True, True, True),
    mip_axis=(1, 0, 0),
    damp_factor=1,
    scale_factor=None,
    decon_offset=0,
    mask_fullpaths=(),
    parse_cluster=True,
    parse_parfor=False,
    master_compute=True,
    master_cpu=False,
    gpu_job=False,
    job_log_dir="../job_logs",
    cpus_per_task=4,
    uuid="",
    max_trial_num=3,
    unit_wait_time=2,
    debug=False,
    save_step=5,
    psf_gen=True,
    mcc_mode=False,
    config_file="",
    gpu_config_file="",
) -> None:
    """Deconvolve each file of a job (XR_RLdeconFrame3D): see the module. Sizes and MIP axes in
    (z, y, x)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.decon import file as F
    from petakit5d.io import get_image_size

    _ = edge_soften, z_edge_soften, master_cpu, max_trial_num, unit_wait_time
    files = (
        [frame_fullpaths] if isinstance(frame_fullpaths, (str, os.PathLike)) else frame_fullpaths
    )
    files = [os.fspath(f) for f in files]
    if not psf_file:
        raise DeconFileError("", "You should provide a PSF file for the frame...\\n")
    method = str(rl_method).lower()
    bg = 100 if background is None or np.size(background) == 0 else background
    if not uuid:
        uuid = get_uuid()
    if parse_cluster and not _slurm_available():
        _warn("A slurm-based computing cluster is not exist. Set parseCluster as false.")
        parse_cluster = False
    save3 = [bool(v) for v in np.asarray(save_3d_stack).reshape(-1)]
    if len(save3) == 1:
        save3 = save3 * 3
    save3[1] = save3[1] and bool(deskew)
    save3[2] = save3[2] and bool(rotate)
    mips = [float(v) for v in np.asarray(mip_axis).reshape(-1)]
    if not large_file and not any(save3) and not any(v > 0 for v in mips):
        raise DeconFileError(
            "", "Either the 3D stack or the MIP for at least one axis needs to be saved!"
        )
    large = bool(large_file)
    mask_file = erode_maskfile
    n = len(files)
    prefetched: tuple[str, Any] | None = None
    erosion = float(np.asarray(edge_erosion, dtype=np.float64).reshape(-1)[0])
    for f, fn in enumerate(files):
        if not _exist(fn):
            _warn(f"{fn} does not exist! Skip it!")
            continue
        pathstr, fsname, ext = fileparts(fn)
        dpath = decon_path
        if not dpath:
            dpath = f"{pathstr}/matlab_decon"
            os.makedirs(dpath, exist_ok=True)
        dsr_folder = f"{dpath}/DSR/"
        if save3[2]:
            os.makedirs(dsr_folder, exist_ok=True)
        im_size = tuple(int(v) for v in get_image_size(fn))
        if _mem_required(im_size) > float(cpu_max_mem):
            large = True
        result = _result_path(dpath, fsname, save_zarr)
        dsr_result = f"{dsr_folder}/{fsname}{'.zarr' if save_zarr else '.tif'}"
        if _exist(result):
            if not overwrite:
                print("Deconvolution results already exist, skip it!", flush=True)
                continue
            # the decon functions skip existing results, so remove them
            _remove(result)
            if save3[2]:
                _remove(dsr_result)
        if not large:
            if erosion > 0 and save_maskfile:
                mask_file = _save_mask(
                    fn,
                    ext,
                    dpath,
                    fsname,
                    {
                        "flip_zstack": flip_zstack,
                        "edge_erosion": edge_erosion,
                        "uuid": uuid,
                        "block_size": block_size,
                    },
                )
            t0 = time.perf_counter()
            rawdata = None
            if prefetched is not None and prefetched[0] == fn:
                rawdata = prefetched[1]
            prefetched = None
            # with several files on the GPU, the next file is read during the decon, if the GPU
            # has room for the decon twice over
            overlap_fcn = None
            next_fn = None
            if gpu_job and f < n - 1 and gpu_has_room(im_size, method):
                next_fn = files[f + 1]
                if next_frame_needs_decon(next_fn, decon_path, save_zarr, overwrite, cpu_max_mem):
                    overlap_fcn = lambda p=next_fn: read_frame(p)  # noqa: E731
            args = (fn, result, psf_file, xy_pixel_size, dz, dz_psf)
            kw = {
                "rawdata": rawdata,
                "save_16bit": save_16bit,
                "skew_angle": skew_angle,
                "deskew": deskew,
                "rotate": rotate,
                "dsr_combined": True,
                "reverse": reverse,
                "background": bg,
                "decon_iter": decon_iter,
                "rl_method": method,
                "skewed": skewed,
                "wiener_alpha": wiener_alpha,
                "otf_cum_thresh": otf_cum_thresh,
                "hann_win_bounds": hann_win_bounds,
                "save_zarr": save_zarr,
                "block_size": block_size,
                "fix_iter": fix_iter,
                "damp_factor": damp_factor,
                "scale_factor": scale_factor,
                "decon_offset": decon_offset,
                "edge_erosion": edge_erosion,
                "erode_maskfile": mask_file,
                "err_thresh": err_thresh,
                "save_step": save_step,
                "psf_gen": psf_gen,
                "debug": debug,
                "save_3d_stack": tuple(save3),
                "mip_axis": mip_axis,
                "flip_zstack": flip_zstack,
            }
            del rawdata
            next_data = None
            try:
                if overlap_fcn is None:
                    F.rl_decon_file(*args, **kw, use_gpu=bool(gpu_job))
                else:
                    r = F.rl_decon_file(*args, **kw, use_gpu=bool(gpu_job), overlap_fcn=overlap_fcn)
                    next_data = r.overlap_out if r is not None else None
            except Exception as e:
                # a file that doesn't fit on the GPU (memory or the gpuArray size limit) is
                # deconvolved on the CPU instead of failing
                if not (gpu_job and is_gpu_capacity_error(e)):
                    raise
                _warn(f"{fn} does not fit on the GPU ({e}), deconvolve it on the CPU instead...")
                F.rl_decon_file(*args, **kw, use_gpu=False)
            kw = None
            if next_data is not None:
                prefetched = (next_fn, next_data)
            next_data = None
            if (save3[0] and _exist(result)) or (save3[2] and _exist(dsr_result)):
                print(f"Done! Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)
                continue

        # large-file deconvolution: the file is too large, or the decon in one piece wrote no
        # result
        common = {
            "save_16bit": save_16bit,
            "deskew": deskew,
            "skew_angle": skew_angle,
            "flip_zstack": flip_zstack,
            "background": bg,
            "dz_psf": dz_psf,
            "decon_iter": decon_iter,
            "rl_method": method,
            "skewed": skewed,
            "wiener_alpha": wiener_alpha,
            "otf_cum_thresh": otf_cum_thresh,
            "hann_win_bounds": hann_win_bounds,
            "fix_iter": fix_iter,
            "batch_size": batch_size,
            "block_size": block_size,
            "damp_factor": damp_factor,
            "scale_factor": scale_factor,
            "decon_offset": decon_offset,
            "edge_erosion": edge_erosion,
            "mask_fullpaths": mask_fullpaths,
            "uuid": uuid,
            "debug": debug,
            "psf_gen": psf_gen,
        }
        if str(large_method).lower() == "inplace":
            _large("decon_large_in_place")(
                frame_fullpaths=fn,
                pixel_size=xy_pixel_size,
                dz=dz,
                decon_path=dpath,
                psf=psf_file,
                **common,
                parse_cluster=parse_cluster,
                parse_parfor=parse_parfor,
                master_compute=master_compute,
                job_log_dir=job_log_dir,
                cpus_per_task=cpus_per_task,
                gpu_job=gpu_job,
                mcc_mode=mcc_mode,
                config_file=config_file,
                gpu_config_file=gpu_config_file,
            )
            return
        if str(large_method).lower() == "inmemory":
            _large("decon_large_in_memory")(
                frame_fullpaths=fn,
                psf_fullpath=psf_file,
                decon_fullpath=result,
                xy_pixel_size=xy_pixel_size,
                dz=dz,
                **common,
                save_zarr=save_zarr,
                use_gpu=gpu_job,
            )
            return
