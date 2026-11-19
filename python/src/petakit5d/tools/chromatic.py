"""Chromatic shift correction (MATLAB's ``microscopeDataProcessing/chromatic_shift_correction``):
the shifts between channels measured from PSFs (``XR_chromatic_shift_estimation``) and the
datasets' channels aligned by integer voxel offsets
(``XR_chromatic_shift_correction_data_wrapper``).

The estimation: each PSF after the first (cropped by ``crop_length`` on each side) is placed in
the first by the normalized cross-correlation of ``kernels.xcorr.normxcorr3_max_shift`` (single
data) within ``max_offset`` of the crop; its offset less the crop is the channel's shift, and the
shifts, the bounds and the peaks are saved in MATLAB's ``.mat`` (``-v7.3``, values in
``[y x z]``) when a result file is named.

The correction: channel ``c``'s voxel ``i`` shows the scene at ``i + offset[c]``; each file is
cropped (and padded with ``pad_value``) by ``tools.crop.crop_file`` to one region of the scene,
one task per file through the job framework, into ``<data path>/<result_dir_name>/`` with the
input's name (``.zarr`` with ``save_zarr``). The region by ``mode``:

- ``'valid'``: the overlap of the channels (no padding);
- ``'same'``: the input's size, at ``new_origin`` (0-based, in the scene), or without it centered
  on the overlap (``overlap start + floor((overlap size - image size) / 2)``; the origin of this
  mode is the one the MATLAB tests define, tests/README.md);
- ``'full'``: the union of the channels.

Modes are taken in any case; others are refused (MATLAB's parser:
``MATLAB:InputParser:ArgumentFailedValidation``). Large zarr files (``zarr_file`` and
``large_file``) are cropped in batches (``crop_file`` with ``large_file``, MATLAB's
``XR_crop_zarr``) into zarr results, their jobs sized as ``XR_crop_dataset`` sizes them; their
padding is 0 (another ``pad_value`` is refused where they are padded). Offsets must be a row per
channel pattern (``...:noOffsets`` without offsets and PSFs, ``...:offsetSize`` otherwise).

Arrays are ``(z, y, x)``; offsets are rows ``(z, y, x)`` per channel (in ``channel_patterns``'
order), shifts, bounds and crops ``(z, y, x)``; the MATLAB-name layer takes MATLAB's ``[y x z]``
rows and a 1-based ``newOrigin``. The NCC runs on the process's CPU FFT backend
(``petakit5d.fft``): MATLAB's bits with ``exact``, within decision record 0006 with
``permissive``.

The estimation takes one path as one PSF (too few). MATLAB quirks kept: ``padValue`` goes to the
tasks through ``%d`` (a non-integer value as MATLAB's ``%e`` text).

Differences: offsets and origins must be integers (MATLAB passes others to XR_crop_frame, whose
mex truncate them).
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "ChromaticError",
    "chromatic_shift_correction_dataset",
    "chromatic_shift_estimation",
    "correction_boxes",
]

_DIMS = "MATLAB:sizeDimensionsMustMatch"
_DIMS_MSG = "Arrays have incompatible sizes for this operation."


class ChromaticError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _floats(v: Any) -> np.ndarray:
    return np.asarray(v, dtype=np.float64).ravel()


def _per_axis(v: Any) -> np.ndarray:
    """a (z, y, x) vector, one number for every axis"""
    a = _floats(v)
    return np.broadcast_to(a, (3,)).copy() if a.size == 1 else a


def _yxz_value(v: Any) -> np.ndarray:
    """a per-axis value as MATLAB holds it (a double row; one number stays one)"""
    a = _floats(v)
    return a.reshape(1, -1) if a.size != 3 else a[[1, 2, 0]].reshape(1, 3)


def _empty(v: Any) -> bool:
    return v is None or np.asarray(v).size == 0


def _offsets_in(v: Any) -> Any:
    """MATLAB's rows of [y x z] offsets as rows (z, y, x)"""
    if _empty(v):
        return None
    a = np.asarray(v, dtype=np.float64)
    a = a.reshape(1, -1) if a.ndim < 2 else a
    if a.shape[1] != 3:
        return tuple(tuple(r) for r in a.tolist())
    return tuple(tuple(r) for r in a[:, [2, 0, 1]].tolist())


def _offsets_out(v: Any) -> Any:
    """rows (z, y, x) as MATLAB's rows of [y x z]"""
    if _empty(v):
        return None
    a = np.asarray(v, dtype=np.float64)
    a = a.reshape(1, -1) if a.ndim < 2 else a
    if a.shape[1] != 3:
        return a
    return a[:, [1, 2, 0]]


def _integers(a: np.ndarray, what: str) -> np.ndarray:
    if not np.all(np.isfinite(a)) or np.any(a != np.trunc(a)):
        raise ValueError(f"{what} must be integers (voxels)")
    return a.astype(np.int64)


def _read(path: str, zarr_file: bool) -> np.ndarray:
    from petakit5d.io import read_tiff, read_zarr

    im = read_zarr(path) if zarr_file else read_tiff(path)
    return im if im.ndim == 3 else im.reshape(1, *im.shape)


# ---------------------------------------------------------------------------------------------
# XR_chromatic_shift_estimation
# ---------------------------------------------------------------------------------------------


def _two_channels(
    im_s: np.ndarray, im_t: np.ndarray, max_offset: np.ndarray, crop: np.ndarray, nt: int | None
) -> tuple[np.ndarray, Any]:
    """chromatic_shift_estimation_for_two_channels: the shift (z, y, x) of ``im_t`` (cropped)
    in ``im_s`` and the peak of the NCC"""
    from petakit5d.kernels.xcorr import normxcorr3_max_shift

    c = _integers(crop, "crop_length")
    n = im_t.shape
    t1 = im_t[c[0] : n[0] - c[0], c[1] : n[1] - c[1], c[2] : n[2] - c[2]]
    lower, upper = crop - max_offset, crop + max_offset
    bounds = np.stack([lower, upper])
    off, peak, _ = normxcorr3_max_shift(
        np.ascontiguousarray(t1, dtype=np.float32),
        np.ascontiguousarray(im_s, dtype=np.float32),
        bounds,
        threads=nt,
    )
    return np.asarray(off, dtype=np.float64) - crop, peak


def _matlab_out(shifts: Any) -> Any:
    """the shifts as MATLAB's rows of [y x z]"""
    return np.asarray(shifts, dtype=np.float64)[:, [1, 2, 0]]


@register(
    "XR_chromatic_shift_estimation",
    Param("psfFullpaths", "@(x) ischar(x) || iscell(x)", "PSF files, one per channel."),
    Param(
        "resultFullpath",
        "@(x) ischar(x) || iscell(x)",
        "Result .mat file of the measurement (empty: not saved).",
    ),
    Param(
        "maxOffset",
        "@isnumeric",
        "Largest shift (z, y, x) searched, around the crop.",
        convert="size_yxz",
    ),
    Param(
        "cropLength",
        "@isnumeric",
        "Voxels (z, y, x) cut from each side of the PSFs after the first.",
        convert="size_yxz",
    ),
    Param("zarrFile", "@islogical", "The PSFs are zarr."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Debug mode (not used)."),
    Param("num_threads", None, "Threads of the NCC (None: all available).", python_only=True),
    task="B03",
    matlab_out=_matlab_out,
)
def chromatic_shift_estimation(
    psf_fullpaths: str | Sequence[str],
    *,
    result_fullpath: str = "",
    max_offset: Sequence[float] | float = (30, 30, 30),
    crop_length: Sequence[float] | float = (0, 0, 0),
    zarr_file: bool = False,
    uuid: str = "",
    debug: bool = False,
    num_threads: int | None = None,
) -> np.ndarray:
    """The shifts (rows ``(z, y, x)``, the first channel's 0) of the channels' PSFs against the
    first (``XR_chromatic_shift_estimation``), each the offset of the cropped PSF in the first at
    the peak of their normalized cross-correlation within ``max_offset`` of the crop, less the
    crop. With ``result_fullpath`` the paths, bounds, shifts and peaks are saved there (MATLAB's
    ``.mat``, ``[y x z]`` values; published by one rename)."""
    del debug
    if isinstance(psf_fullpaths, (str, os.PathLike)):
        # (one path is one PSF)
        psf_fullpaths = [os.fspath(psf_fullpaths)]
    paths = [os.fspath(p) for p in psf_fullpaths]
    nf = len(paths)
    if nf < 2:
        raise ChromaticError("At least two PSF files are needed!")
    max_off = _per_axis(max_offset)
    crop = _per_axis(crop_length)

    im_s = _read(paths[0], zarr_file)
    shifts = np.zeros((nf, 3))
    peaks = np.zeros((nf, 1))
    for i in range(1, nf):
        im_t = _read(paths[i], zarr_file)
        shifts[i], peak = _two_channels(im_s, im_t, max_off, crop, num_threads)
        peaks[i, 0] = float(peak)

    if result_fullpath:
        from petakit5d.core.ids import get_uuid
        from petakit5d.core.paths import fileparts
        from petakit5d.io import write_mat

        if not uuid:
            uuid = get_uuid()
        result_fullpath = os.fspath(result_fullpath)
        path, fsn, ext = fileparts(result_fullpath)
        tmp = f"{path}/{fsn}_{uuid}{ext}"
        write_mat(
            tmp,
            {
                "psfFullpaths": list(psf_fullpaths),
                "maxOffset": _yxz_value(max_offset),
                "cropLength": _yxz_value(crop_length),
                "relative_shift_yxz_mat": _matlab_out(shifts),
                "max_xcorr_mat": peaks,
            },
        )
        os.replace(tmp, result_fullpath)
    return shifts


# ---------------------------------------------------------------------------------------------
# the regions
# ---------------------------------------------------------------------------------------------


def correction_boxes(
    im_size: Sequence[int],
    offsets: Any,
    mode: str = "valid",
    new_origin: Sequence[int] | None = None,
) -> tuple[np.ndarray, bool]:
    """The box ``(z0, y0, x0, z1, y1, x1)`` of each channel (rows) that the wrapper crops from
    its files for ``mode``, in the channel's coordinates, and whether the crops are padded.
    Channel ``c`` covers the scene box ``offsets[c] : offsets[c] + im_size``; ``new_origin``
    (0-based, in the scene) places the ``'same'`` region."""
    n = np.asarray(im_size, dtype=np.int64)[:3]
    if _empty(offsets):
        # (MATLAB's [1 1 1 imSize] + [])
        raise ChromaticError(_DIMS_MSG, _DIMS)
    o = _integers(np.atleast_2d(np.asarray(offsets, dtype=np.float64)), "offsets")
    if o.shape[1] != 3:
        raise ChromaticError(_DIMS_MSG, _DIMS)
    start, stop = o, o + n
    ov_start, ov_stop = start.max(axis=0), stop.min(axis=0)
    if mode == "valid":
        return np.concatenate([ov_start - start, ov_stop - start], axis=1), False
    if mode == "same":
        if _empty(new_origin):
            origin = ov_start + (ov_stop - ov_start - n) // 2
        else:
            origin = _integers(_floats(new_origin), "new_origin")
        return np.concatenate([origin - start, origin - start + n], axis=1), True
    if mode == "full":
        origin = start.min(axis=0)
        size = stop.max(axis=0) - origin
        return np.concatenate([origin - start, origin - start + size], axis=1), True
    return np.concatenate([start, stop], axis=1), True


# ---------------------------------------------------------------------------------------------
# XR_chromatic_shift_correction_data_wrapper
# ---------------------------------------------------------------------------------------------


@register(
    "XR_chromatic_shift_correction_data_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param(
        "chromaticOffset",
        "@(x) isnumeric(x)",
        "Offsets in voxels, a row (z, y, x) per channel (None: measured from the PSFs).",
        convert=(_offsets_in, _offsets_out),
    ),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param(
        "mode",
        "@(x) ischar(x) && any(strcmpi(x, {'valid', 'same', 'full'}))",
        "Output region, in any case: 'valid' (overlap), 'same' (input size) or 'full'.",
    ),
    Param("padValue", "@isnumeric", "Value of the voxels without data.", convert="same"),
    Param(
        "newOrigin",
        "@isnumeric",
        "Origin (z, y, x) of the 'same' region in the scene (None: centered on the overlap).",
        convert="point1_yxz",
    ),
    Param("channelPatterns", "@iscell", "Channel identifiers, one per channel."),
    Param("psfFullpaths", "@iscell", "PSF files, one per channel (for measured offsets)."),
    Param(
        "maxOffset",
        "@isnumeric",
        "Largest offset (z, y, x) of the measurement.",
        convert="size_yxz",
    ),
    Param(
        "cropLength",
        "@isnumeric",
        "Voxels (z, y, x) cut from each side of the PSFs for the measurement.",
        convert="size_yxz",
    ),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param(
        "largeFile",
        "@islogical",
        "Large files: zarr files cropped in batches into zarr, tiff files in memory.",
    ),
    Param("saveZarr", "@islogical", "Save the results as zarr."),
    Param("batchSize", "@isvector", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of zarr results.", convert="size_yxz"),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Debug mode (not used)."),
    task="B03",
)
def chromatic_shift_correction_dataset(
    data_paths: str | Sequence[str],
    *,
    chromatic_offset: Any = None,
    result_dir_name: str = "Chromatic_Shift_Corrected",
    mode: str = "valid",
    pad_value: float = 0,
    new_origin: Sequence[int] | None = None,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    psf_fullpaths: Sequence[str] = (),
    max_offset: Sequence[float] | float = (20, 20, 20),
    crop_length: Sequence[float] | float = (0, 0, 0),
    zarr_file: bool = False,
    large_file: bool = False,
    save_zarr: bool = False,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    parse_cluster: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 3,
    config_file: str = "",
    mcc_mode: bool = False,
    uuid: str = "",
    debug: bool = False,
) -> None:
    """Align the channels of the datasets by integer voxel offsets
    (``XR_chromatic_shift_correction_data_wrapper``): with ``chromatic_offset`` (rows
    ``(z, y, x)``, one per channel pattern) or, without, the offsets measured from
    ``psf_fullpaths`` (``chromatic_shift_estimation``, saved in
    ``<first result folder>/measurement/``). Each file is cropped to the ``mode``'s region by a
    ``crop_file`` task through the job framework (large zarr files in batches), its result in
    ``<data path>/<result_dir_name>/``."""
    import warnings

    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import simplify_path
    from petakit5d.core.text import through_mat2str, through_text
    from petakit5d.io import get_image_size
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_image_filenames
    from petakit5d.tools.crop import crop_file

    del debug
    if not isinstance(mode, str) or mode.lower() not in ("valid", "same", "full"):
        raise ChromaticError(
            "The value of 'mode' is invalid. It must satisfy the function: "
            "@(x)ischar(x)&&any(strcmpi(x,{'valid','same','full'})).",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )
    mode = mode.lower()
    paths = [data_paths] if isinstance(data_paths, str | os.PathLike) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    # (large zarr files are cropped in batches, into zarr results)
    large_zarr = bool(zarr_file) and bool(large_file)
    if large_zarr:
        save_zarr = True
    nd = len(paths)
    result_paths = []
    for d in range(nd):
        p = paths[d]
        if not p.endswith(("/", os.sep)):
            p = p + os.sep
        paths[d] = simplify_path(p) + "/"
        result_path = f"{paths[d]}/{result_dir_name}/"
        result_paths.append(result_path)
        os.makedirs(result_path, exist_ok=True)
        group_writable(result_path)
    if not uuid:
        uuid = get_uuid()

    parsed = parse_image_filenames(paths, zarr_file, channel_patterns)
    fnames, fsns = parsed.fnames, parsed.fsns
    fd_inds = [int(v) for v in parsed.fd_inds]
    ch_inds = [int(v) for v in parsed.ch_inds]
    nf = len(fnames)
    nc = len(channel_patterns)
    if len(set(ch_inds)) != nc:
        raise ChromaticError(
            "Number of channels of the images does not match the number of channel patterns!"
        )

    offsets = chromatic_offset
    if _empty(offsets) and len(psf_fullpaths) > 0:
        if len(psf_fullpaths) != nc:
            raise ChromaticError(
                "The number of PSF files does not match the number of the channels!"
            )
        print("Measure chromatic offsets across channels with PSF files...", flush=True)
        print("\n".join(f"    {{'{os.fspath(p)}'}}" for p in psf_fullpaths), flush=True)
        measure_path = f"{result_paths[0]}measurement/"
        os.makedirs(measure_path, exist_ok=True)
        # (the PSFs are tiff files; the data keeps zarr_file)
        shifts = chromatic_shift_estimation(
            psf_fullpaths,
            result_fullpath=f"{measure_path}/chromatic_shift_measurement_result.mat",
            max_offset=max_offset,
            crop_length=crop_length,
            zarr_file=False,
            uuid=uuid,
        )
        offsets = shifts - shifts.min(axis=0)
    if _empty(offsets):
        raise ChromaticError(
            "Give chromaticOffset (a [y, x, z] row per channel pattern) or psfFullpaths to "
            "measure it.",
            "XR_chromatic_shift_correction_data_wrapper:noOffsets",
        )
    shape = np.shape(np.atleast_2d(np.asarray(offsets, dtype=np.float64)))
    if shape != (nc, 3):
        raise ChromaticError(
            f"chromaticOffset must have a [y, x, z] row per channel pattern ({nc} x 3), not "
            f"[{shape[0]} {shape[1]}].",
            "XR_chromatic_shift_correction_data_wrapper:offsetSize",
        )
    print("Chromatic offset across channels:", flush=True)
    print("" if _empty(offsets) else np.asarray(offsets), flush=True)

    im_sizes = np.zeros((nd, 3), dtype=np.int64)
    boxes: list[np.ndarray] = []
    pad = True
    for d in range(nd):
        f_d = fd_inds.index(d)
        im_size = get_image_size(f"{paths[d]}{fnames[f_d]}")[:3]
        im_sizes[d] = im_size
        b, pad = correction_boxes(im_size, offsets, mode, new_origin)
        boxes.append(b)

    if large_zarr and pad and float(pad_value) != 0:
        from petakit5d.core.text import sprintf

        raise ChromaticError(
            sprintf(
                "Large Zarr files are padded with 0 (XR_crop_zarr): padValue %g is not supported.",
                float(pad_value),
            ),
            "XR_chromatic_shift_correction_data_wrapper:padValue",
        )
    ext = ".zarr" if save_zarr else ".tif"
    # (the values as the task strings give them: padValue %d, the sizes through mat2str)
    pv = float(through_text(float(pad_value), "%d"))
    batch = tuple(int(v) for v in np.asarray(through_mat2str(batch_size)).ravel())
    block = tuple(int(v) for v in np.asarray(through_mat2str(block_size)).ravel())
    cpus = float(through_text(float(cpus_per_task), "%d"))
    frame_fullpaths, result_fullpaths, tasks = [], [], []
    for f in range(nf):
        frame = f"{paths[fd_inds[f]]}/{fnames[f]}"
        result = f"{result_paths[fd_inds[f]]}/{fsns[f]}{ext}"
        b = boxes[fd_inds[f]]
        frame_fullpaths.append(frame)
        result_fullpaths.append(result)
        tasks.append(
            Task.of(
                crop_file,
                frame,
                result,
                tuple(int(v) for v in b[ch_inds[f]]),
                pad=pad,
                pad_value=pv,
                zarr_file=bool(zarr_file),
                save_zarr=bool(save_zarr),
                batch_size=batch,
                block_size=block,
                uuid=uuid,
                parse_cluster=bool(parse_cluster),
                master_compute=bool(master_compute),
                cpus_per_task=cpus,
                mcc_mode=bool(mcc_mode),
                config_file=config_file,
                # (only then, as MATLAB passes it)
                **({"large_file": True} if large_zarr else {}),
            )
        )

    job_cpus = cpus_per_task
    coord_nice = 0
    if large_zarr:
        # each file's job (crop_zarr) runs its batches as jobs: it holds the process and the
        # merge of the MIP slabs, and if it crops batches itself cpus_per_task and a batch's
        # memory, as crop_dataset sizes them
        from petakit5d.io import get_image_dtype
        from petakit5d.jobs.config import cluster_config_value
        from petakit5d.jobs.coordinator import coordinator_nice, coordinator_resources
        from petakit5d.tools.mip import merge_memory

        allb = np.concatenate(boxes, axis=0)
        out_size = [int(v) for v in (allb[:, 3:6] - allb[:, 0:3]).max(axis=0)]
        batch_b = [min(o, int(b)) for o, b in zip(out_size, batch, strict=True)]
        block_b = [min(o, int(b)) for o, b in zip(out_size, block, strict=True)]
        merge_mem = merge_memory(out_size, batch_b, block_b, get_image_dtype(frame_fullpaths[0]))
        job_cpus, mem_allocate = coordinator_resources(
            merge_mem,
            cluster_config_value(config_file, "masterCompute", master_compute),
            cpus_per_task,
            float(np.prod(np.asarray(batch_b, dtype=np.float64))) * 4 * 2 / 1024**3,
        )
        # (and they wait for their batch jobs: these go first)
        if parse_cluster:
            coord_nice = coordinator_nice()
    else:
        # (a large tiff file is cropped in memory)
        mem_allocate = float(np.prod(im_sizes.max(axis=0).astype(np.float64))) / 1024**3 * 4 * 2.5
    common = {
        "parse_cluster": parse_cluster,
        "job_log_dir": job_log_dir,
        "master_compute": master_compute,
        "cpus_per_task": job_cpus,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
        "nice": coord_nice,
    }
    done = framework.run_tasks(
        frame_fullpaths, result_fullpaths, tasks, mem_allocate=mem_allocate, **common
    )
    if not all(done):
        done = framework.run_tasks(
            frame_fullpaths, result_fullpaths, tasks, mem_allocate=mem_allocate * 2, **common
        )
    if not all(done):
        warnings.warn("Some files are not processed!", stacklevel=2)
