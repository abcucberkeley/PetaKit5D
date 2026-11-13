"""Cropping (MATLAB's ``microscopeDataProcessing/crop``): files (``XR_crop_frame``), datasets
(``XR_crop_dataset``), large zarr files in blocks with MIPs (``XR_crop_zarr``,
``XR_crop_block``) and the border trim of volumes (``trimBorder``).

Arrays are ``(z, y, x)`` C order; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive,
and may reach outside the image (with ``pad``, the part outside is filled); sizes, batch and
block sizes and borders are ``(z, y, x)``. The MATLAB-name layer converts MATLAB's values
(``[ymin xmin zmin ymax xmax zmax]``, 1-based, inclusive; ``[y x z]`` sizes; 1-based batch
numbers). Boxes are integral: MATLAB's non-integer boxes (truncated by the mex, rounded by its
indexing) are not taken.

Files keep MATLAB's names and contents: ``<data>/<result_dir_name>/<name>.tif`` or ``.zarr``
(the format follows ``save_zarr``), ``parameters.mat`` of the call, the MIPs of large crops in
``<result>/MIPs/`` (``mip_zarr``, merged from the MIPs the block tasks save in their flags), and
the flags ``zarr_flag/<name>_<uuid>/blocks_<i>_<j>.mat`` of the block tasks (MATLAB's format).
Zarr outputs are native arrays, as every zarr the port writes.

Results are published by one rename of a temporary result (``<name>_<uuid>.tif``,
``<name>_<uuid>.zarr``, ``<name>.zarr_<uuid>`` for large crops), after everything else the task
writes (a large crop's MIPs): the job framework counts a task as done once its result exists.
A large crop records its writer (``jobs.writers``), waits for a writer of another call of the
same result instead of making it a second time, keeps a result another writer published, and
its blocks skip their work once its writer has ended.

MATLAB quirks kept (each has a test):

- ``crop_zarr`` refuses a box past the image's end without ``pad``, but pads a box that starts
  before the image (its blocks always pad);
- ``crop_dataset`` passes neither ``cpus_per_task`` nor ``job_log_dir`` to the file tasks (the
  blocks of large crops ask for 1 CPU and log to ``../job_logs``);
- a moving crop of several datasets takes the time points of the last dataset's files,
  numbered over all files (MATLAB's ``ch_tps(f)``): an error when the last has fewer files;
- ``trim_border`` takes a scalar border only for ``'post'``, as MATLAB indexes ``s(2)`` of a
  scalar otherwise.

Differences: ``overwrite`` of ``crop_file`` removes an existing zarr result before the new one is
moved into place (MATLAB's movefile moved the new folder into the old one).
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "BLOCK_MEMORY_FACTOR",
    "CropError",
    "block_memory",
    "crop_block",
    "crop_dataset",
    "crop_file",
    "crop_zarr",
    "trim_border",
]

# GB a block task of a large crop asks for: this many times its batch (in the data's class).
# MATLAB asks for 8 bytes per voxel (prod(batchSize) * 4 * 2), 4 times a uint16 batch. The
# port's block task holds one batch and, for a batch that leaves the image, the part of it that
# is read before it is placed, then the writer's chunk buffers. Peak RSS of one task with a
# 1024^3 uint16 batch (2 GiB, random data, chunks 256^3, with the blocks' MIPs; c24 node,
# 2026-11-12): inside the image 1.27 / 1.78 times the batch at 8 / 24 CPUs, cut by one plane
# (the worst case: the batch and its part) 2.02 at 24 CPUs, half outside 1.28, wholly outside
# 0.40. 2.2 covers the worst case.
BLOCK_MEMORY_FACTOR = 2.2

_BADSUBSCRIPT = "MATLAB:badsubscript"


class CropError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _ints(v: Any) -> list[int]:
    out = []
    for x in np.asarray(v).ravel():
        f = float(x)
        if not math.isfinite(f) or f != math.floor(f):
            raise ValueError(f"expected integers, got {v!r}")
        out.append(int(f))
    return out


def _box(box: Any) -> list[int]:
    b = _ints(box)
    if len(b) != 6:
        raise ValueError(f"a box has 6 values (z0, y0, x0, z1, y1, x1), got {len(b)}")
    return b


def _size3(v: Any, name: str) -> list[int]:
    s = _ints(v)
    if len(s) != 3:
        raise ValueError(f"{name} has 3 values (z, y, x), got {len(s)}")
    return s


def _vol3(a: np.ndarray) -> np.ndarray:
    """an image (y, x) as the volume (1, y, x)"""
    return a[np.newaxis] if a.ndim == 2 else a


def _image_size(path: str) -> list[int]:
    from petakit5d.io import get_image_size

    return [int(v) for v in get_image_size(path)[:3]]


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _numpy_dtype(cls: str) -> np.dtype:
    from petakit5d.core.arrays import numpy_dtype

    return numpy_dtype(cls)


def _fill_value(value: float, dtype: np.dtype) -> Any:
    """padarray's constant: ``cast(value, class(im))`` (round half away from zero, saturate)"""
    from petakit5d.core.arrays import matlab_class
    from petakit5d.core.dtypes import cast

    return np.asarray(cast(np.asarray([float(value)]), matlab_class(dtype))).ravel()[0]


def _toc(t0: float) -> None:
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)


def block_memory(batch_size: Sequence[int], dtype: str) -> float:
    """GB a block task of ``crop_zarr`` asks for: ``BLOCK_MEMORY_FACTOR`` times a batch of
    ``batch_size`` (``(z, y, x)``) in the class ``dtype``."""
    from petakit5d.core.arrays import data_type_to_byte_number

    voxels = float(np.prod([float(v) for v in batch_size], dtype=np.float64))
    return voxels * data_type_to_byte_number(dtype) / 2**30 * BLOCK_MEMORY_FACTOR


# ---------------------------------------------------------------------------------------------
# trimBorder
# ---------------------------------------------------------------------------------------------


def trim_border(vol: np.ndarray, border_size: Any, method: str = "both") -> np.ndarray:
    """``trimBorder``: ``vol`` (``(z, y, x)``) without ``border_size`` (``(z, y, x)``) voxels at
    the start (``'pre'``), the end (``'post'``) or both ends (``'both'``) of each axis, as a new
    C-contiguous array; ``vol`` itself if the border is all zero. A border larger than the
    volume leaves an empty axis. MATLAB's errors: a scalar border for ``'pre'`` or ``'both'``
    (``s(2)`` of a scalar), a negative border that reaches outside, an unknown method."""
    from petakit5d.kernels.crop import crop_3d

    b = np.asarray(border_size, dtype=np.float64).ravel()
    if b.size and np.all(b == 0):
        return vol
    a = np.asarray(vol)
    if method not in ("pre", "post", "both"):
        raise CropError("Unrecognized function or variable 's'.", "MATLAB:UndefinedFunction")
    sz = list(a.shape)
    if a.ndim != 3:
        raise CropError(
            "Arrays have incompatible sizes for this operation.", "MATLAB:sizeDimensionsMustMatch"
        )
    if b.size == 1 and method != "post":
        raise CropError(
            "Index exceeds the number of array elements. Index must not exceed 1.", _BADSUBSCRIPT
        )
    if b.size not in (1, 3):
        raise CropError(
            "Arrays have incompatible sizes for this operation.", "MATLAB:sizeDimensionsMustMatch"
        )
    bb = _ints(np.broadcast_to(b, (3,)))
    lo = bb if method in ("pre", "both") else [0, 0, 0]
    hi = [n - (k if method in ("post", "both") else 0) for n, k in zip(sz, bb, strict=True)]
    # MATLAB indexes y, then x, then z (positions 1, 2, 3); an empty range is never checked
    box = [0] * 6
    for pos, p in enumerate((1, 2, 0), 1):
        s, t = lo[p], hi[p]
        if t > s:
            if s < 0:
                raise CropError(
                    f"Index in position {pos} is invalid. Array indices must be positive "
                    "integers or logical values.",
                    _BADSUBSCRIPT,
                )
            if t > sz[p]:
                raise CropError(
                    f"Index in position {pos} exceeds array bounds. Index must not exceed {sz[p]}.",
                    _BADSUBSCRIPT,
                )
        else:
            s, t = 0, 0
        box[p], box[p + 3] = s, t
    return crop_3d(a, box)


# ---------------------------------------------------------------------------------------------
# files
# ---------------------------------------------------------------------------------------------


def _read_box(path: str, box: list[int], zarr_file: bool) -> np.ndarray:
    """the box of a tiff or zarr file (MATLAB's readzarr with inputBbox; readtiff of the box's
    slices and crop3d of its rows and columns, with crop3d's errors)"""
    from petakit5d.core.bbox import box_to_matlab
    from petakit5d.io import read_tiff, read_zarr
    from petakit5d.kernels.crop import crop3d_matlab

    if zarr_file:
        return _vol3(read_zarr(path, input_bbox=box))
    im = _vol3(read_tiff(path, range=(box[0], box[3])))
    m = list(box_to_matlab(box))
    m[2], m[5] = 1, im.shape[0]
    return _vol3(np.ascontiguousarray(crop3d_matlab(im, m)))


def _padded(im: np.ndarray, box: list[int], im_size: list[int], pad_value: float) -> np.ndarray:
    """padarray of the read part to the box (``pad_value`` before and after)"""
    from petakit5d.kernels.crop import indexing_3d

    pre = [max(0, -box[q]) for q in range(3)]
    post = [max(0, box[q + 3] - im_size[q]) for q in range(3)]
    if not any(pre) and not any(post):
        return im
    shape = tuple(im.shape[q] + pre[q] + post[q] for q in range(3))
    fill = _fill_value(pad_value, im.dtype)
    if fill == 0 and not np.signbit(fill):
        out = np.zeros(shape, dtype=im.dtype)
    else:
        out = np.full(shape, fill, dtype=im.dtype)
    place = [*pre, *(pre[q] + im.shape[q] for q in range(3))]
    return indexing_3d(out, place, im)


@register(
    "XR_crop_frame",
    Param("dataFullpath", "@(x) ischar(x) || iscell(x)", "Input tiff or zarr file."),
    Param("saveFullpath", "@(x) ischar(x) || iscell(x)", "Result file (.tif or .zarr)."),
    Param(
        "bbox",
        "@isnumeric",
        "Box (z0, y0, x0, z1, y1, x1); its end is cut at the image's end without pad.",
        convert="bbox",
    ),
    Param("overwrite", "@islogical", "Replace an existing result."),
    Param("pad", "@islogical", "Pad the part of the box outside the image."),
    Param("padValue", "@isnumeric", "Value of the padding.", convert="same"),
    Param("zarrFile", "@islogical", "The input is zarr."),
    Param("largeFile", "@islogical", "Large zarr input: crop in blocks (crop_zarr)."),
    Param("saveZarr", "@islogical", "Save the result as zarr."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of zarr results.", convert="size_yxz"),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("parseCluster", "@islogical", "Run the blocks of large files as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes blocks too."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a block's job."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T03",
)
def crop_file(
    data_fullpath: str,
    save_fullpath: str,
    bbox: Sequence[int],
    *,
    overwrite: bool = False,
    pad: bool = False,
    pad_value: float = 0,
    zarr_file: bool = False,
    large_file: bool = False,
    save_zarr: bool = False,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    uuid: str = "",
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 1,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Crop a tiff or zarr file with a box (``XR_crop_frame``) into ``save_fullpath`` (tiff, or
    zarr with ``save_zarr``, chunks ``block_size``). Without ``pad`` the box's end is cut at
    the image's end; with it, the part outside the image is ``pad_value`` (cast to the image's
    class). Large zarr files are cropped in blocks (``crop_zarr``). Nothing is done if the
    input is missing (a warning) or the result exists (unless ``overwrite``)."""
    import warnings

    from petakit5d.core.ids import get_uuid
    from petakit5d.io import write_tiff, write_zarr
    from petakit5d.jobs.publish import publish_result, remove_path

    if not uuid:
        uuid = get_uuid()
    data_fullpath, save_fullpath = os.fspath(data_fullpath), os.fspath(save_fullpath)
    if not _exist(data_fullpath):
        warnings.warn(f"The file {data_fullpath} does not exist!", stacklevel=2)
        return
    if _exist(save_fullpath) and not overwrite:
        print(f"The cropped file {save_fullpath} is already exist, skip it!", flush=True)
        return
    box = _box(bbox)
    print(f"Crop frame {data_fullpath} with bounding box {tuple(box)}...", flush=True)

    im_size = _image_size(data_fullpath)
    outside = any(box[q] < 0 or box[q + 3] > im_size[q] for q in range(3))
    if pad and outside:
        box_1 = [max(0, v) for v in box[:3]] + [
            min(n, v) for n, v in zip(im_size, box[3:], strict=True)
        ]
    else:
        box_1 = box[:3] + [min(n, v) for n, v in zip(im_size, box[3:], strict=True)]

    if zarr_file and large_file:
        crop_zarr(
            data_fullpath,
            save_fullpath,
            box,
            pad=pad,
            batch_size=batch_size,
            block_size=block_size,
            uuid=uuid,
            parse_cluster=parse_cluster,
            parse_parfor=False,
            master_compute=master_compute,
            cpus_per_task=cpus_per_task,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        return
    im = _read_box(data_fullpath, box_1, zarr_file)
    if pad:
        im = _padded(im, box, im_size, pad_value)

    if save_zarr:
        tmp = f"{save_fullpath[:-5]}_{uuid}.zarr"
        write_zarr(im, tmp, block_size=_size3(block_size, "block_size"))
        if overwrite and os.path.isdir(save_fullpath):
            remove_path(save_fullpath)
        publish_result(tmp, save_fullpath)
    else:
        tmp = f"{save_fullpath[:-4]}_{uuid}.tif"
        write_tiff(im, tmp)
        os.replace(tmp, save_fullpath)


# ---------------------------------------------------------------------------------------------
# large zarr files
# ---------------------------------------------------------------------------------------------


@register(
    "XR_crop_block",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based).", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("cropFullpath", "@(x) ischar(x)", "Result zarr (made before), written into."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch in the input.",
        convert="bboxes",
        py="batch_bboxes",
    ),
    Param(
        "RegionBBoxes",
        "@isnumeric",
        "Box of each batch in the result.",
        convert="bboxes",
        py="region_bboxes",
    ),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param(
        "mipBatchSize",
        "@isnumeric",
        "Save the MIPs of the blocks in the flag, for MIP slabs of this batch size (z, y, x).",
        convert="size_yxz",
    ),
    Param("uuid", "@ischar", "Id of the call."),
    Param("debug", "@islogical", "Not used."),
    task="T03",
)
def crop_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    crop_fullpath: str,
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    region_bboxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    mip_batch_size: Sequence[int] | None = None,
    uuid: str = "",
    debug: bool = False,
) -> list[bool] | None:
    """``XR_crop_block``: each batch (``batch_bboxes``, boxes of the input zarr that may reach
    outside it: zeros there) written into its region (``region_bboxes``) of the result zarr
    ``crop_fullpath``; then the flag: empty, or with the MIPs of the blocks for MIP slabs of
    ``mip_batch_size`` (``save_block_mips``). Nothing is done if the flag exists (unless
    ``overwrite``) or if the writer of the crop has ended (``writers.orphaned``)."""
    from petakit5d.io import ReadBuffer, get_zarr_info, read_zarr, write_zarr
    from petakit5d.jobs import writers
    from petakit5d.kernels.crop import indexing_3d
    from petakit5d.tools.mip import block_mips, save_block_mips

    del uuid, debug
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise CropError(f"The block directory {flag_path} does not exist, skip the processing!")
    inds = _ints(batch_inds)
    if _exist(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            print(
                f"The batch files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!",
                flush=True,
            )
            return None
    # a block of a task whose job ended (e.g., was cancelled): nothing would use it
    if writers.orphaned(flag_fullname):
        return None
    if not os.path.isdir(zarr_fullpath):
        raise CropError(f"The input zarr file {zarr_fullpath} doesnot exist!")

    im_size = _image_size(zarr_fullpath)
    dtype = get_zarr_info(crop_fullpath).dtype
    ndtype = _numpy_dtype(dtype)
    in_boxes = np.asarray(batch_bboxes, dtype=np.int64).reshape(-1, 6)
    out_boxes = np.asarray(region_bboxes, dtype=np.int64).reshape(-1, 6)
    mips: list[list[np.ndarray]] = []
    mbs = None if mip_batch_size is None or len(mip_batch_size) == 0 else list(mip_batch_size)

    done = [False] * len(inds)
    # (one buffer for the blocks' reads; each block is written before the next is read)
    buf = ReadBuffer()
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        ib = [int(v) for v in in_boxes[i]]
        ob = [int(v) for v in out_boxes[i]]
        shape = tuple(ob[q + 3] - ob[q] for q in range(3))
        overlap = not (
            any(ib[q + 3] <= 0 for q in range(3)) or any(ib[q] >= im_size[q] for q in range(3))
        )
        if not overlap:
            out = np.zeros(shape, dtype=ndtype)
        else:
            cut = [max(0, v) for v in ib[:3]] + [
                min(n, v) for n, v in zip(im_size, ib[3:], strict=True)
            ]
            if cut != ib:
                part = _vol3(read_zarr(zarr_fullpath, input_bbox=cut, out=buf))
                out = np.zeros(shape, dtype=ndtype)
                place = [cut[q] - ib[q] for q in range(3)] + [cut[q + 3] - ib[q] for q in range(3)]
                indexing_3d(out, place, part)
                del part
            else:
                out = _vol3(read_zarr(zarr_fullpath, input_bbox=cut, out=buf))
        write_zarr(out, crop_fullpath, bbox=ob, create=False)
        if mbs is not None:
            mips.append(block_mips(out, ob, mbs, dtype))
        # (one batch in memory)
        del out
        done[i] = True
        _toc(t0)

    if all(done):
        if mbs is None:
            open(flag_fullname, "w").close()
        else:
            save_block_mips(flag_fullname, [list(map(int, r)) for r in out_boxes], mips, mbs)
    return done


@register(
    "XR_crop_zarr",
    Param("zarrFullpath", "@(x) ischar(x) || iscell(x)", "Input zarr."),
    Param("cropFullpath", "@(x) ischar(x) || iscell(x)", "Result zarr."),
    Param("bbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) to crop.", convert="bbox"),
    Param("pad", "@islogical", "Pad the part of the box past the image's end."),
    Param("batchSize", "@isvector", "Batch size (z, y, x) of a block.", convert="size_yxz"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the result.", convert="size_yxz"),
    Param("saveMIP", "@islogical", "Save the MIPs of the result (MIPs/)."),
    Param("parseCluster", "@islogical", "Run the blocks as Slurm jobs."),
    Param("parseParfor", "@islogical", "Run the blocks on a local process pool."),
    Param("masterCompute", "@islogical", "This process computes blocks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a block's job."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T03",
)
def crop_zarr(
    zarr_fullpath: str,
    crop_fullpath: str,
    bbox: Sequence[int],
    *,
    pad: bool = False,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    save_mip: bool = True,
    parse_cluster: bool = True,
    parse_parfor: bool = False,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 1,
    uuid: str = "",
    debug: bool = False,
    mcc_mode: bool = False,
    config_file: str = "",
) -> list[bool] | None:
    """Crop a large zarr file in blocks (``XR_crop_zarr``): the batches of the box
    (``batch_size``, whole chunks of ``block_size``) as ``crop_block`` tasks through the job
    framework into ``<crop>_<uuid>``, its MIPs (``mip_zarr``, merged from the MIPs of the
    blocks), then the result moved into place. A box past the image's end needs ``pad``.
    Nothing is done if the result exists; a writer of another call of the same result is
    waited for, and a result another writer published is kept. Returns the tasks' done
    flags (None if nothing was done)."""
    from petakit5d.core.blocks import zarr_chunk_coordinates
    from petakit5d.core.dtypes import mround
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive, strip_trailing_sep
    from petakit5d.io import create_zarr, get_zarr_info
    from petakit5d.jobs import framework, writers
    from petakit5d.jobs.publish import keep_other_result, publish_result
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.tools.mip import mip_zarr, mip_zarr_batch_size

    if not uuid:
        uuid = get_uuid()
    zarr_fullpath = os.fspath(zarr_fullpath)
    crop_fullpath = strip_trailing_sep(os.fspath(crop_fullpath))
    crop_path, fsname, _ext = fileparts(crop_fullpath)
    if os.path.isdir(crop_fullpath):
        print("The output result exists, skip it!", flush=True)
        return None
    print(f"Start Large-file Cropping for {fsname}...", flush=True)
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(crop_path, job_log_dir).parse_cluster)

    flag_path = f"{crop_path}/zarr_flag/{fsname}_{uuid}/"
    if not os.path.isdir(flag_path):
        mkdir_recursive(flag_path)
    tmp_path = f"{crop_fullpath}_{uuid}"
    # a writer of another call that may still make the crop: wait for it; the record is
    # marked ended when this call ends but by publishing (its blocks are skipped then)
    with writers.recorded_writer(crop_fullpath, uuid, flag_path):
        if writers.wait(crop_fullpath, uuid):
            keep_other_result(crop_fullpath, tmp_path, flag_path)
            return None

        info = get_zarr_info(zarr_fullpath)
        im_size = [int(v) for v in info.size[:3]]
        dtype = info.dtype
        box = _box(bbox)
        if not pad and any(box[q + 3] > im_size[q] for q in range(3)):
            raise CropError(
                "The bounding box is out of bound! Please check the bbox setting or turn on "
                "pad option."
            )
        out_size = [box[q + 3] - box[q] for q in range(3)]
        batch = [min(o, b) for o, b in zip(out_size, _size3(batch_size, "batch_size"), strict=True)]
        block = [min(o, b) for o, b in zip(out_size, _size3(block_size, "block_size"), strict=True)]
        batch_boxes, region_boxes, _ = zarr_chunk_coordinates(
            out_size,
            batch_size=batch,
            block_size=block,
            same_batch_size=False,
            border_size=0,
        )
        # the batches in the input
        batch_boxes = batch_boxes + np.asarray(box[:3] * 2, dtype=np.int64)

        if not os.path.isdir(tmp_path):
            n_chunks = 1
            for o, b in zip(out_size, block, strict=True):
                n_chunks *= math.ceil(o / b)
            create_zarr(
                tmp_path,
                data_size=out_size,
                block_size=block,
                dtype=dtype,
                dim_separator="/" if n_chunks > 10000 else ".",
            )

        # the tasks save the MIPs of their blocks for the MIP slabs (mip_zarr merges them
        # instead of reading the result again); the option only then, as MATLAB passes it
        mip_kw: dict[str, Any] = {}
        if save_mip:
            zi = get_zarr_info(tmp_path)
            mbs = mip_zarr_batch_size(zi.size, zi.block_size, zi.dtype)
            mip_kw["mip_batch_size"] = [int(v) for v in mbs]

        num_batch = batch_boxes.shape[0]
        task_size = max(5, int(mround(num_batch / 5000)))
        num_tasks = math.ceil(num_batch / task_size)
        tasks, outputs = [], []
        for t in range(num_tasks):
            inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
            flag = f"{flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
            outputs.append(flag)
            tasks.append(
                Task.of(
                    crop_block,
                    inds,
                    zarr_fullpath,
                    tmp_path,
                    flag,
                    batch_boxes[inds].tolist(),
                    region_boxes[inds].tolist(),
                    uuid=uuid,
                    debug=bool(debug),
                    **mip_kw,
                )
            )
        inputs = [zarr_fullpath] * num_tasks
        if parse_cluster or not parse_parfor:
            # (the blocks stop once the result exists: another writer made it)
            mem = block_memory(batch, dtype)
            done: list[bool] = []
            for k in range(2):
                done = framework.run_tasks(
                    inputs,
                    outputs,
                    tasks,
                    cpus_per_task=cpus_per_task,
                    max_job_num=math.inf,
                    mem_allocate=mem * 2**k,
                    task_batch_num=1,
                    master_compute=master_compute,
                    parse_cluster=parse_cluster,
                    job_log_dir=job_log_dir,
                    mcc_mode=mcc_mode,
                    config_file=config_file,
                    final_out_fullpath=crop_fullpath,
                )
                if all(done):
                    break
        else:
            done = framework.run_tasks(
                inputs,
                outputs,
                tasks,
                cluster_type="parfor",
                task_batch_num=1,
                gpu_job=False,
                local_worker_num=12,
                uuid=uuid,
            )
        if not all(done):
            raise CropError("Some blocks are not finished!")

        # another writer may have published the result: keep it and its MIPs (the same
        # result), and remove this one's temporary result, which can be incomplete
        if keep_other_result(crop_fullpath, tmp_path, flag_path):
            return done

        # the MIPs before the result is moved into place (the framework counts the task as
        # done once it exists), merged from the MIPs of the blocks in the flags
        if save_mip and os.path.isdir(tmp_path):
            mip_zarr(
                tmp_path,
                axis=(1, 1, 1),
                parse_cluster=parse_cluster,
                master_compute=master_compute,
                job_log_dir=job_log_dir,
                mcc_mode=mcc_mode,
                config_file=config_file,
                block_flag_fullpaths=outputs,
            )
        # (one may also have published it during the MIPs)
        publish_result(tmp_path, crop_fullpath, flag_path)
    return done


# ---------------------------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------------------------


def _save_parameters(result_path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) of the call (MATLAB's ``pr``)"""
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(crop_dataset, values).items()}
    write_mat(result_path + "/parameters.mat", {"pr": pr})


def _list_files(data_path: str, zarr_file: bool) -> list[str]:
    from petakit5d.core.paths import matlab_dir

    return [e.name for e in matlab_dir(f"{data_path}/*.{'zarr' if zarr_file else 'tif'}")]


@register(
    "XR_crop_dataset",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) to crop.", convert="bbox"),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param(
        "cropType",
        "@ischar",
        "'fixed' (the box), 'moving' (moved linearly to last_start_coords over the time "
        "points) or 'center' (a box of the box's size at the image's center).",
    ),
    Param("pad", "@islogical", "Pad the part of the box outside the image."),
    Param(
        "lastStartCoords",
        "@isnumeric",
        "Start (z, y, x) of the box at the last time point ('moving').",
        convert="point1_yxz",
    ),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to crop."),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("largeFile", "@islogical", "Large zarr files: crop in blocks (crop_zarr)."),
    Param("saveZarr", "@islogical", "Save the results as zarr."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of zarr results.", convert="size_yxz"),
    Param("save16bit", "@islogical", "Not used.", convert="same"),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T03",
)
def crop_dataset(
    data_paths: str | Sequence[str],
    input_bbox: Sequence[int],
    *,
    result_dir_name: str = "Cropped",
    crop_type: str = "fixed",
    pad: bool = False,
    last_start_coords: Sequence[int] | None = None,
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0"),
    zarr_file: bool = False,
    large_file: bool = False,
    save_zarr: bool = False,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    save_16bit: bool = True,
    parse_cluster: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 2,
    uuid: str = "",
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Crop every file of the datasets that matches a channel pattern (``XR_crop_dataset``)
    into ``<data path>/<result_dir_name>/<name>.tif`` (or ``.zarr`` with ``save_zarr``; large
    zarr files always give zarr), one ``crop_file`` task per file through the job framework.
    ``crop_type``: ``'fixed'`` crops ``input_bbox``; ``'moving'`` moves it linearly from its
    start at the first time point to ``last_start_coords`` at the last (per channel, rounded
    as MATLAB rounds); ``'center'`` crops a box of its size at each image's center.
    parameters.mat of the call goes into each result folder."""
    values = {k: v for k, v in locals().items() if k in _DATASET_ARGS}
    from petakit5d.core.dtypes import mround
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, simplify_path
    from petakit5d.io import get_image_dtype
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.config import cluster_config_value
    from petakit5d.jobs.coordinator import coordinator_nice, coordinator_resources
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import matches_channel_pattern
    from petakit5d.tools.mip import merge_memory

    if not uuid:
        uuid = get_uuid()
    paths = [data_paths] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    if large_file:
        save_zarr = True
    box = _box(input_bbox)

    result_paths = []
    for d, data_path in enumerate(paths):
        result_path = f"{data_path}/{result_dir_name}/"
        os.makedirs(result_path, exist_ok=True)
        paths[d] = simplify_path(data_path)
        result_paths.append(simplify_path(result_path))
        group_writable(result_path)
        _save_parameters(result_path, values)

    moving = crop_type == "moving"
    fnames: list[str] = []
    fd_inds: list[int] = []
    max_tps: list[int] = []
    ch_tps: list[int] = []
    sz: list[int] = []
    for d, data_path in enumerate(paths):
        names = _list_files(data_path, zarr_file)
        if not names:
            raise CropError(
                "Index exceeds the number of array elements. Index must not exceed 0.",
                _BADSUBSCRIPT,
            )
        sz = _image_size(f"{data_path}/{names[0]}")
        names = sorted(names)
        ch = [
            matches_channel_pattern(names, p, True) for p in channel_patterns
        ]  # (files x channels)
        keep = [any(c[i] for c in ch) for i in range(len(names))]
        fnames += [n for n, k in zip(names, keep, strict=True) if k]
        fd_inds += [d] * sum(keep)
        if moving:
            # the time point of each file in its channel (cumsum over the files kept)
            rows = [[c[i] for c in ch] for i in range(len(names)) if keep[i]]
            count = [0] * len(channel_patterns)
            ch_tps = []
            for r in rows:
                count = [n + int(v) for n, v in zip(count, r, strict=True)]
                ch_tps.append(sum(n * int(v) for n, v in zip(count, r, strict=True)))
            max_tps.append(max(ch_tps) if ch_tps else 0)
            last = _ints(last_start_coords if last_start_coords is not None else [])
            ext = [box[q + 3] - box[q] for q in range(3)]
            if len(last) != 3:
                raise CropError(
                    "Arrays have incompatible sizes for this operation.",
                    "MATLAB:sizeDimensionsMustMatch",
                )
            if any(ext[q] + last[q] > sz[q] for q in range(3)):
                raise CropError(
                    "The moving bounding box moves out of the image for the last time point!"
                )
    n_files = len(fnames)

    frame_fullpaths, crop_fullpaths, tasks = [], [], []
    for f in range(n_files):
        d = fd_inds[f]
        frame = f"{paths[d]}/{fnames[f]}"
        frame_fullpaths.append(frame)
        # the output format follows save_zarr, whatever the input format
        fsn = fileparts(fnames[f])[1]
        crop = f"{result_paths[d]}/{fsn}.{'zarr' if save_zarr else 'tif'}"
        crop_fullpaths.append(crop)
        ext = [box[q + 3] - box[q] for q in range(3)]
        if moving:
            if f >= len(ch_tps):
                raise CropError(
                    "Index exceeds the number of array elements. Index must not exceed "
                    f"{len(ch_tps)}.",
                    _BADSUBSCRIPT,
                )
            last = _ints(last_start_coords)
            start = [
                box[q] + int(mround((last[q] - box[q]) / (max_tps[d] - 1) * (ch_tps[f] - 1)))
                for q in range(3)
            ]
            box_f = start + [s + e for s, e in zip(start, ext, strict=True)]
        elif crop_type == "center":
            im_size = _image_size(frame)
            # center = round((imSize + 1) / 2) (1-based), start = center - floor(size / 2)
            start = [
                int(mround((n + 1) / 2)) - e // 2 - 1 for n, e in zip(im_size, ext, strict=True)
            ]
            box_f = start + [s + e for s, e in zip(start, ext, strict=True)]
        else:
            box_f = list(box)
        kw: dict[str, Any] = {
            "pad": bool(pad),
            "zarr_file": bool(zarr_file),
            "large_file": bool(large_file),
            "save_zarr": bool(save_zarr),
            "batch_size": _size3(batch_size, "batch_size"),
            "block_size": _size3(block_size, "block_size"),
            "uuid": uuid,
            "parse_cluster": bool(parse_cluster),
            "mcc_mode": bool(mcc_mode),
            "config_file": config_file,
        }
        # (only then, as the default of crop_file is true)
        if zarr_file and large_file and not master_compute:
            kw["master_compute"] = False
        tasks.append(Task.of(crop_file, frame, crop, box_f, **kw))

    coord_nice = 0
    if zarr_file and large_file:
        # each file's job (crop_zarr) runs its blocks as jobs: it holds the process and the
        # merge of the MIP slabs, and if it crops blocks itself (master_compute, the config's
        # if it sets it), cpus_per_task and a block's estimate
        out_size = [box[q + 3] - box[q] for q in range(3)]
        batch_b = [
            min(o, b) for o, b in zip(out_size, _size3(batch_size, "batch_size"), strict=True)
        ]
        block_b = [
            min(o, b) for o, b in zip(out_size, _size3(block_size, "block_size"), strict=True)
        ]
        dtype = get_image_dtype(frame_fullpaths[0])
        merge_mem = merge_memory(out_size, batch_b, block_b, dtype)
        computes = cluster_config_value(config_file, "masterCompute", master_compute)
        cpus_per_task, mem_allocate = coordinator_resources(
            merge_mem, computes, cpus_per_task, block_memory(batch_b, dtype)
        )
        # (and they wait for their block jobs: these go first)
        if parse_cluster:
            coord_nice = coordinator_nice()
    else:
        # (MATLAB's 32 bytes per voxel of the last dataset's first file)
        mem_allocate = float(np.prod(sz, dtype=np.float64)) * 4 / 1024**3 * 8
    framework.run_tasks(
        frame_fullpaths,
        crop_fullpaths,
        tasks,
        parse_cluster=parse_cluster,
        job_log_dir=job_log_dir,
        master_compute=master_compute,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        mcc_mode=mcc_mode,
        config_file=config_file,
        nice=coord_nice,
    )


_DATASET_ARGS = (
    "data_paths",
    "input_bbox",
    "result_dir_name",
    "crop_type",
    "pad",
    "last_start_coords",
    "channel_patterns",
    "zarr_file",
    "large_file",
    "save_zarr",
    "batch_size",
    "block_size",
    "save_16bit",
    "parse_cluster",
    "master_compute",
    "job_log_dir",
    "cpus_per_task",
    "uuid",
    "mcc_mode",
    "config_file",
)
