"""The block task of the in-place large-file decon (``RLdecon_for_zarr_block``).

``decon_zarr_block`` deconvolves its batches of a zarr file and writes each batch's region into
the result zarr that ``decon.large.decon_large_in_place`` made (``<result>_<uuid>``), then its
flag ``blocks_<first>_<last>.mat``: with ``mip_batch_size`` the MIPs of its blocks for the MIP
slabs (``tools.mip.save_block_mips``, MATLAB's format), else the time it took (``t``).

As RLdecon_for_zarr_block.m at the oracle commit:

- nothing is done if the flag exists (unless ``overwrite``, which removes it), or if the writer
  of the result has ended (``jobs.writers.orphaned``: its record says so, or its Slurm job is
  gone); the flag's folder, the input and the result must exist;
- a batch whose region is empty in one of the 2D masks (``mask_fullpaths``: xy, xz, yz) is
  written as zeros without a decon; the others are read with their border, deconvolved as
  RLdecon does (``decon.file.rl_decon_file`` with the batch as ``rawdata``, its region as
  ``decon_bbox``, no deskew or rotation, iterations fixed, ``use_gpu`` as given: a GPU job's
  block fails without a usable GPU, so that the framework retries it) and written into the
  result in its class;
- every fifth batch of a CPU task looks for its flag, and stops if another job made it.

Batch numbers are 0-based (the flag names keep MATLAB's 1-based numbers), boxes are ours
(``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive). ``flip_zstack`` and ``fix_iter`` are
taken and not used, as in MATLAB.
"""

from __future__ import annotations

import os
import time
from collections.abc import Sequence

import numpy as np

from petakit5d.params import Param, register

__all__ = ["decon_zarr_block"]


def _ints(v) -> list[int]:
    a = np.asarray(v, dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(a)) or np.any(a != np.floor(a)):
        raise ValueError(f"expected integers, got {v!r}")
    return [int(x) for x in a]


def _boxes(v) -> np.ndarray:
    a = np.asarray(v, dtype=np.int64)
    return a.reshape(-1, 6)


def _save_time_flag(flag: str, t: float) -> None:
    """MATLAB's ``save(flagFullname, 't')`` (a temporary name, then moved into place: the flag
    marks the task as done)"""
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import write_mat

    tmp = f"{flag[:-4]}_{get_uuid()}.mat"
    write_mat(tmp, {"t": float(t)}, version="7")
    os.replace(tmp, flag)


@register(
    "RLdecon_for_zarr_block",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based).", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("psfFullpath", "@(x) ischar(x)", "PSF file."),
    Param("deconFullpath", "@(x) ischar(x)", "Result zarr (made before), written into."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch (with its border) in the input.",
        convert="bboxes",
        py="batch_bboxes",
    ),
    Param(
        "RegionBBoxes",
        "@isnumeric",
        "Box of each batch's region in the result.",
        convert="bboxes",
        py="region_bboxes",
    ),
    Param("xyPixelSize", "@isnumeric", "Pixel size (um)."),
    Param("dz", "@isnumeric", "Scan step (um)."),
    Param("save16bit", "@islogical", "Result in uint16 (else single).", convert="same"),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("SkewAngle", "@isnumeric", "Skew angle (degrees)."),
    Param("flipZstack", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("Background", "@isnumeric", "Background.", convert="same"),
    Param("dzPSF", "@isnumeric", "z step of the PSF (um)."),
    Param("DeconIter", "@isnumeric", "Number of iterations."),
    Param("dampFactor", "@isnumeric", "Damp factor of the result."),
    Param("scaleFactor", "@isnumeric", "Scale factor of the result ([]: 1)."),
    Param("deconOffset", "@isnumeric", "Offset of the result.", convert="same"),
    Param("EdgeErosion", "@isnumeric", "Erosion of the edges of each batch's mask."),
    Param("maskFullpaths", "@iscell", "2D masks (xy, xz, yz) of the regions to decon."),
    Param("RLMethod", "@ischar", "'simplified', 'omw', 'original' or 'cudagen'."),
    Param("wienerAlpha", "@isnumeric", "Wiener alpha of the OMW back projector."),
    Param("OTFCumThresh", "@isnumeric", "OTF cumulative sum threshold (OMW)."),
    Param("hannWinBounds", "@isnumeric", "Apodization range of the distance matrix (OMW)."),
    Param("skewed", "@(x) isempty(x) || islogical(x)", "Decon in skewed space (None: automatic)."),
    Param("fixIter", "@islogical", "Accepted; the blocks fix the iterations (as in MATLAB)."),
    Param("useGPU", "@islogical", "Deconvolve on a GPU (an error if none is usable)."),
    Param(
        "mipBatchSize",
        "@isnumeric",
        "Save the MIPs of the blocks in the flag, for MIP slabs of this batch size (z, y, x).",
        convert="size_yxz",
    ),
    Param("uuid", "@ischar", "Id of the call."),
    Param("debug", "@islogical", "Not used."),
    Param("psfGen", "@islogical", "Generate the PSF (psf_gen) or use it as it is."),
    task="R07",
)
def decon_zarr_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    psf_fullpath: str,
    decon_fullpath: str,
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    region_bboxes: Sequence[Sequence[int]],
    xy_pixel_size,
    dz,
    *,
    save_16bit: bool = True,
    overwrite: bool = False,
    skew_angle=-32.45,
    flip_zstack: bool = False,
    background=None,
    dz_psf=0.1,
    decon_iter=15,
    damp_factor=1,
    scale_factor=1.0,
    decon_offset=0,
    edge_erosion=0,
    mask_fullpaths: Sequence[str] = (),
    rl_method: str = "simplified",
    wiener_alpha=0.005,
    otf_cum_thresh=0.9,
    hann_win_bounds=(0.8, 1),
    skewed=None,
    fix_iter: bool = False,
    use_gpu: bool = False,
    mip_batch_size: Sequence[int] | None = None,
    uuid: str = "",
    debug: bool = False,
    psf_gen: bool = True,
) -> list[bool] | None:
    """``RLdecon_for_zarr_block``: deconvolve the batches of a large zarr file into the result
    zarr and write the task's flag (see the module). Returns the batches' done flags (None if
    nothing was done)."""
    from petakit5d.core.arrays import numpy_dtype
    from petakit5d.decon.large import (
        DeconLargeError,
        _local_box,
        _rl_batch,
        _use_masks,
        masks_empty,
    )
    from petakit5d.io import ReadBuffer, get_image_dtype, read_zarr, write_zarr
    from petakit5d.jobs import writers
    from petakit5d.jobs.files import exist
    from petakit5d.tools.mip import block_mips, save_block_mips

    del flip_zstack, fix_iter, uuid, debug
    t0 = time.perf_counter()
    zarr_fullpath = os.fspath(zarr_fullpath)
    decon_fullpath = os.fspath(decon_fullpath)
    flag_fullname = os.fspath(flag_fullname)
    # the flag's folder must exist (else all blocks were done and it was removed)
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise DeconLargeError(
            f"The block directory {flag_path} does not exist, skip the processing!"
        )
    inds = _ints(batch_inds)
    if exist(flag_fullname):
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
        raise DeconLargeError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    if not os.path.isdir(decon_fullpath):
        raise DeconLargeError(f"The output zarr file {decon_fullpath} doesnot exist!")

    dtype = get_image_dtype(decon_fullpath, True)
    ndtype = numpy_dtype(dtype)
    in_boxes = _boxes(batch_bboxes)
    out_boxes = _boxes(region_bboxes)
    mbs = None if mip_batch_size is None or len(mip_batch_size) == 0 else _ints(mip_batch_size)
    masks = _use_masks(mask_fullpaths)
    kw = {
        "save_16bit": bool(save_16bit),
        "skew_angle": skew_angle,
        "deskew": False,
        "rotate": False,
        "background": background,
        "decon_iter": decon_iter,
        "rl_method": str(rl_method),
        "skewed": skewed,
        "wiener_alpha": wiener_alpha,
        "otf_cum_thresh": otf_cum_thresh,
        "hann_win_bounds": hann_win_bounds,
        "fix_iter": True,
        "damp_factor": damp_factor,
        "scale_factor": scale_factor,
        "decon_offset": decon_offset,
        "edge_erosion": edge_erosion,
        "use_gpu": bool(use_gpu),
        "psf_gen": bool(psf_gen),
    }

    mips: list[list[np.ndarray]] = []
    done = [False] * len(inds)
    # (one buffer for the batches' reads; a batch is deconvolved before the next is read)
    buf = ReadBuffer()
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1}... ", end="", flush=True)
        tb = time.perf_counter()
        ib = [int(v) for v in in_boxes[i]]
        ob = [int(v) for v in out_boxes[i]]
        if masks and masks_empty(masks, ob):
            out = np.zeros(tuple(ob[q + 3] - ob[q] for q in range(3)), dtype=ndtype)
            write_zarr(out, decon_fullpath, bbox=ob, create=False)
            if mbs is not None:
                mips.append(block_mips(out, ob, mbs, dtype))
            del out
            done[i] = True
            print(f"Elapsed time is {time.perf_counter() - tb:f} seconds.", flush=True)
            continue
        in_batch = read_zarr(zarr_fullpath, input_bbox=ib, out=buf)
        out = _rl_batch(
            in_batch,
            decon_fullpath,
            psf_fullpath,
            xy_pixel_size,
            dz,
            dz_psf,
            _local_box(ob, ib),
            kw,
        )
        del in_batch
        out = np.asarray(out)
        write_zarr(out, decon_fullpath, bbox=ob, create=False)
        if mbs is not None:
            mips.append(block_mips(out, ob, mbs, dtype))
        del out
        done[i] = True
        print(f"Elapsed time is {time.perf_counter() - tb:f} seconds.", flush=True)
        # the task may have been finished by another job: look for the flag every five
        # batches (not on GPUs)
        if not use_gpu and (i + 1) % 5 == 0 and exist(flag_fullname):
            break
    buf.release()

    if all(done):
        t = time.perf_counter() - t0
        if mbs is None:
            _save_time_flag(flag_fullname, t)
        else:
            save_block_mips(flag_fullname, [list(map(int, r)) for r in out_boxes], mips, mbs, t)
    return done
