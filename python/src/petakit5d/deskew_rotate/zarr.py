"""Deskew and rotation of large zarr files in blocks (``XR_deskewRotateZarr``,
``XR_deskewRotateBlock``) and the boxes of MIP masks (``XR_getDeskeRotateBoxesFromMasks``).

``deskew_rotate_zarr`` splits the combined deskew/rotation of one zarr file (sample scan only)
into batches of output rows along y, each across the whole x and z, with one more input row on
each side (MATLAB's borders). ``deskew_rotate_block`` tasks run through the job framework: each
reads its batches' input rows, deskews and rotates them with ``deskew_rotate_volume`` on the
output box of the batch, writes them into the temporary result ``<DSR>/<name>.zarr_<uuid>``
and saves its flag ``<DSR>/zarr_flag/<name>_<uuid>/blocks_<i>_<j>.mat`` (with the MIPs of its
blocks, for ``tools.mip.mip_zarr`` to merge into the MIP slabs instead of reading the result
again). Then the MIPs ``<DSR>/MIPs/<name>_MIP_{y,x,z}.tif``, and last the one rename that
publishes the result: the job framework counts a file as done once its result exists.

A call records its writer (``jobs.writers``): it waits for a live writer of another call of
the same result instead of making it a second time, keeps a result another writer published
(before or during its MIPs) and removes its own temporary result; its blocks skip their work
once its writer has ended (a cancelled job). Without all blocks done nothing is published.

With MIP masks (``mask_fullpaths``: the xy, xz and yz MIPs of a mask), the input box and the
output region come from the masks (``deskew_rotate_boxes_from_masks``): the mask's box with a
border, and the deskewed and rotated y MIP's box with a border.

Arrays are ``(z, y, x)``, boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive; sizes
``(z, y, x)``; ``resample_factor`` ``(z, y, x)`` (``deskew_rotate_volume``'s). The batch and
block geometry is computed in MATLAB's numbers (``[y x z]``, 1-based) as MATLAB computes it
and converted at the end, so that the batches, their borders and the tasks are MATLAB's. Each
block runs the deskew kernels on MATLAB's layout (``deskew_rotate_volume_matlab_layout``) and
is written from it: the output voxels are MATLAB's.

MATLAB quirks kept (each has a test):

- the batches read the resample factor as ``[y x z]`` (``rs(1)`` for the rows of a batch, the
  result's size divided by ``rs``), where ``deskewRotateFrame3D`` reads it as ``[x y z]``: only
  factors with equal x and y give the result the blocks write; a factor of one or two values
  raises MATLAB's index error;
- ``flip_zstack`` is accepted and not used (the blocks never flip);
- ``parse_parfor`` without a cluster raises MATLAB's error for its undefined ``GPUJob``.
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
    "DeskewRotateZarrError",
    "deskew_rotate_block",
    "deskew_rotate_boxes_from_masks",
    "deskew_rotate_zarr",
    "large_file_geometry",
]

_INTERP = "@(x) any(strcmpi(x, {'cubic', 'linear'}))"
_OBJECTIVE = "Currently objective scan is not supported for large scale deskew/rotation!\\n"


class DeskewRotateZarrError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _logical(x: Any) -> bool:
    return bool(np.asarray(x).reshape(-1)[0]) if np.size(x) else False


def _toc(t0: float) -> None:
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)


def _mround(x: float) -> float:
    """MATLAB's round of a double (half away from zero)"""
    from petakit5d.deskew_rotate.geometry import _round

    return _round(float(x))


def _resample_matlab(resample_factor: Any) -> list[float] | None:
    """the MATLAB vector of a resample factor (``[x y z]``, ``[xy z]`` or one value), None for
    none"""
    if resample_factor is None or np.size(resample_factor) == 0:
        return None
    from petakit5d.params.kinds import get_kind

    v = np.asarray(resample_factor, dtype=np.float64).reshape(-1)
    if v.size == 1:
        return [float(v[0])]
    return [float(x) for x in get_kind("resample_xyz").to_matlab(tuple(v.tolist()))]


def _index_error(n: int) -> DeskewRotateZarrError:
    return DeskewRotateZarrError(
        f"Index exceeds the number of array elements. Index must not exceed {n}.",
        "MATLAB:badsubscript",
    )


def _rs3(rs: list[float] | None) -> list[float]:
    """MATLAB's ``rs([1, 2, 3])`` of ``rs = resampleFactor`` (``[1, 1, 1]`` if empty)"""
    if rs is None:
        return [1.0, 1.0, 1.0]
    if len(rs) < 3:
        raise _index_error(len(rs))
    return [float(v) for v in rs[:3]]


# ---------------------------------------------------------------------------------------------
# XR_getDeskeRotateBoxesFromMasks
# ---------------------------------------------------------------------------------------------


def _read_mask(path: str) -> np.ndarray:
    from petakit5d.io import read_tiff, read_zarr
    from petakit5d.jobs.files import fileparts

    ext = fileparts(path)[2]
    if ext in (".tif", ".tiff"):
        im = read_tiff(path)
    elif ext == ".zarr":
        im = read_zarr(path)
    else:
        raise DeskewRotateZarrError(
            "Unrecognized function or variable 'im'.", "MATLAB:UndefinedFunction"
        )
    im = np.asarray(im)
    # (a MIP is an image; MATLAB's 2D indexing of a volume would take its first slice's view)
    return im if im.ndim == 2 else im.reshape(im.shape[-2:]) if im.shape[0] == 1 else im[0]


def _box2_matlab(im: np.ndarray) -> list[int]:
    """``getImageBoundingBox`` of an image in MATLAB's values ``[r0 c0 r1 c1]``"""
    from petakit5d.io import get_image_bounding_box

    b = get_image_bounding_box(im)
    if not b:
        raise DeskewRotateZarrError(
            "Unable to perform assignment because the size of the left side is 1-by-4 and the "
            "size of the right side is 0-by-0.",
            "MATLAB:subsassigndimmismatch",
        )
    return [b[0] + 1, b[1] + 1, b[2], b[3]]


def deskew_rotate_boxes_from_masks(
    mask_fullpaths: Sequence[str],
    border_size: Sequence[int],
    xy_pixel_size: float,
    dz: float,
    skew_angle: float,
    reverse: bool,
    resample: Sequence[float] | None = None,
    input_bbox: Sequence[int] | None = None,
) -> tuple[tuple[int, ...], tuple[int, ...], np.ndarray]:
    """``XR_getDeskeRotateBoxesFromMasks``: the input box and the output region of a large-file
    deskew/rotation from the MIPs of a mask (``mask_fullpaths``: the z, y and x MIP files, in
    MATLAB's orientation: ``ny x nx``, ``nx x nz`` and ``ny x nz``).

    The input box is the mask's box (from the three MIPs) widened by ``border_size`` (``(z, y,
    x)``) within ``input_bbox`` (None: the whole image). The output region: the box of the y
    MIP's crop deskewed and rotated (as a volume of one row), widened by the border, and all
    rows of the input box (divided by the resample factor's first value, as MATLAB does).
    ``resample``: MATLAB's ``[y x z]`` vector of ``deskew_rotate_zarr`` (None: none).

    Returns the input box, the output region (Python boxes of the image and of the full
    result) and the deskewed y MIP (MATLAB's ``nx x nz`` image).
    """
    from petakit5d.core import dtypes
    from petakit5d.deskew_rotate.volume import deskew_rotate_volume

    bd = [int(v) for v in np.asarray(border_size).reshape(-1)]
    bd_yxz = [bd[1], bd[2], bd[0]]
    boxes = []
    sz = [0, 0, 0]
    im_y = None
    for i in range(3):
        im = _read_mask(os.fspath(mask_fullpaths[i]))
        boxes.append(_box2_matlab(im))
        if i == 0:
            sz[0], sz[1] = int(im.shape[0]), int(im.shape[1])
        elif i == 1:
            sz[2] = int(im.shape[1])
            im_y = im
    b1, b2, b3 = boxes
    bbox = [
        min(b1[0], b3[0]),
        min(b1[1], b2[0]),
        min(b2[1], b3[1]),
        max(b1[2], b3[2]),
        max(b1[3], b2[2]),
        max(b2[3], b3[3]),
    ]
    if input_bbox is None or np.size(input_bbox) == 0:
        in_m = [1, 1, 1, *sz]
    else:
        from petakit5d.core.bbox import box_to_matlab

        in_m = [int(v) for v in box_to_matlab(tuple(int(v) for v in input_bbox))]
    in_box = [max(in_m[k], bbox[k] - bd_yxz[k]) for k in range(3)] + [
        min(in_m[k + 3], bbox[k + 3] + bd_yxz[k]) for k in range(3)
    ]
    crop = np.asarray(im_y)[in_box[1] - 1 : in_box[4], in_box[2] - 1 : in_box[5]]
    # permute(im_y, [3, 1, 2]): a volume of one row, x along the crop's rows, z its columns
    frame = np.ascontiguousarray(dtypes.to_class(crop, "uint16").T[:, np.newaxis, :])
    rs = [1.0, 1.0, 1.0]
    if resample is not None and len(resample):
        if len(resample) < 3:
            raise _index_error(len(resample))
        rs[0], rs[2] = float(resample[1]), float(resample[2])
    dsr = deskew_rotate_volume(
        frame,
        skew_angle,
        dz,
        xy_pixel_size,
        reverse=reverse,
        bbox=None,
        objective_scan=False,
        resample_factor=(rs[2], rs[1], rs[0]),
        interp_method="linear",
        x_step_thresh=2,
        save_16bit=all(r == 1 for r in rs),
    )
    # squeeze: MATLAB's nx x nz image
    dsr2 = np.ascontiguousarray(dsr[:, 0, :].T)
    ob = _box2_matlab(dsr2)
    r1 = 1.0 if resample is None or not len(resample) else float(resample[0])
    out_m = [
        1,
        max(1, ob[0] - bd_yxz[1]),
        max(1, ob[1] - bd_yxz[2]),
        int(_mround((in_box[3] - in_box[0] + 1) / r1)),
        min(dsr2.shape[0], ob[2] + bd_yxz[1]),
        min(dsr2.shape[1], ob[3] + bd_yxz[2]),
    ]
    from petakit5d.core.bbox import box_from_matlab

    return tuple(box_from_matlab(in_box)), tuple(box_from_matlab(out_m)), dsr2


# ---------------------------------------------------------------------------------------------
# XR_deskewRotateBlock
# ---------------------------------------------------------------------------------------------


def _turn_pages(vol: np.ndarray, n: int) -> np.ndarray:
    from petakit5d.deskew_rotate.file import _turn_pages as turn

    return turn(vol, n)


def _block_mips_zxy(out_zxy: np.ndarray, box: Sequence[int], mbs: Sequence[int], dtype: str):
    """``MIP_zarr_block_mips`` of a block in MATLAB's layout (``(z, x, y)`` C): the MIPs of the
    ``(z, y, x)`` block, from the layout itself (a maximum does not depend on the order of its
    values), the MIP planes turned back"""
    from petakit5d.tools.mip import block_mips

    b = [int(v) for v in box]
    m = [int(v) for v in mbs]
    zb, xb, yb = block_mips(
        out_zxy, (b[0], b[2], b[1], b[3], b[5], b[4]), (m[0], m[2], m[1]), dtype
    )
    # (z, x, y) parts -> (z, y, x): the z MIP (p, x, y), the y MIP (z, x, p), the x MIP (z, p, y)
    return [
        np.ascontiguousarray(np.swapaxes(zb, 1, 2)),
        np.ascontiguousarray(np.swapaxes(yb, 1, 2)),
        np.ascontiguousarray(np.swapaxes(xb, 1, 2)),
    ]


@register(
    "XR_deskewRotateBlock",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based).", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("dsrFullpath", "@(x) ischar(x)", "Result zarr (made before), written into."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "inBatchBBoxes",
        "@isnumeric",
        "Box of each batch in the input (with its border rows).",
        convert="bboxes",
        py="in_batch_bboxes",
    ),
    Param(
        "outBatchBBoxes",
        "@isnumeric",
        "Box of each batch in the full result (with its border rows).",
        convert="bboxes",
        py="out_batch_bboxes",
    ),
    Param(
        "outRegionBBoxes",
        "@isnumeric",
        "Box of each batch in the result zarr.",
        convert="bboxes",
        py="out_region_bboxes",
    ),
    Param(
        "outLocalBboxes",
        "@isnumeric",
        "Box of each batch's region inside the batch.",
        convert="bboxes",
    ),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Scan interval in um."),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("SkewAngle", "@isnumeric", "Skew angle in degrees."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("flipZstack", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("save16bit", "@islogical", "uint16 result (else single).", convert="same"),
    Param("interpMethod", _INTERP, "Interpolation: 'linear' or 'cubic'."),
    Param(
        "resampleFactor",
        "@(x) isempty(x) || isnumeric(x)",
        "Resampling after rotation (z, y, x).",
        convert="resample_xyz",
    ),
    Param(
        "mipBatchSize",
        "@isnumeric",
        "Save the MIPs of the blocks in the flag, for MIP slabs of this batch size (z, y, x).",
        convert="size_yxz",
    ),
    Param("uuid", "@ischar", "Id of the call."),
    Param("debug", "@islogical", "Not used."),
    Param(
        "numThreads",
        None,
        "Threads of the reads, kernels and writes (None: the CPUs available).",
        py="num_threads",
        python_only=True,
    ),
    task="D04",
)
def deskew_rotate_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    dsr_fullpath: str,
    flag_fullname: str,
    in_batch_bboxes: Sequence[Sequence[int]],
    out_batch_bboxes: Sequence[Sequence[int]],
    out_region_bboxes: Sequence[Sequence[int]],
    out_local_bboxes: Sequence[Sequence[int]],
    xy_pixel_size: float,
    dz: float,
    *,
    overwrite: bool = False,
    skew_angle: float = 32.45,
    reverse: bool = True,
    flip_zstack: bool = False,
    save_16bit: bool = True,
    interp_method: str = "linear",
    resample_factor: Sequence[float] | None = None,
    mip_batch_size: Sequence[int] | None = None,
    uuid: str = "",
    debug: bool = False,
    num_threads: int | None = None,
) -> list[bool] | None:
    """``XR_deskewRotateBlock``: each batch (rows of the input zarr, ``in_batch_bboxes``)
    deskewed and rotated on the output box its region has in the full result
    (``out_local_bboxes`` moved by ``out_batch_bboxes``' x and z starts), written into its
    region (``out_region_bboxes``) of the result zarr ``dsr_fullpath``; then the flag: ``t``,
    or with the MIPs of the blocks for MIP slabs of ``mip_batch_size``. Nothing is done if the
    flag exists (unless ``overwrite``) or if the writer of the result has ended
    (``writers.orphaned``). Returns the batches' done flags (None if nothing was done)."""
    from petakit5d.core import dtypes
    from petakit5d.core.threads import resolve_threads
    from petakit5d.deskew_rotate.volume import deskew_rotate_volume_matlab_layout
    from petakit5d.io import ReadBuffer, get_image_dtype, read_zarr, write_mat, write_zarr
    from petakit5d.jobs import writers
    from petakit5d.tools.mip import save_block_mips

    t_start = time.perf_counter()
    del flip_zstack, uuid, debug
    n = resolve_threads(num_threads)
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise DeskewRotateZarrError(
            f"The block directory {flag_path} does not exist, skip the processing!"
        )
    inds = [int(v) for v in np.asarray(batch_inds).reshape(-1)]
    if os.path.isfile(flag_fullname):
        if _logical(overwrite):
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
        raise DeskewRotateZarrError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    if not os.path.isdir(dsr_fullpath):
        raise DeskewRotateZarrError(f"The output zarr file {dsr_fullpath} doesnot exist!")

    save_16bit = _logical(save_16bit)
    reverse = _logical(reverse)
    dtype = get_image_dtype(zarr_fullpath)
    mbs = None if mip_batch_size is None or np.size(mip_batch_size) == 0 else list(mip_batch_size)
    out_dtype = get_image_dtype(dsr_fullpath, True) if mbs is not None else ""
    in_boxes = np.asarray(in_batch_bboxes, dtype=np.int64).reshape(-1, 6)
    ob_boxes = np.asarray(out_batch_bboxes, dtype=np.int64).reshape(-1, 6)
    or_boxes = np.asarray(out_region_bboxes, dtype=np.int64).reshape(-1, 6)
    ol_boxes = np.asarray(out_local_bboxes, dtype=np.int64).reshape(-1, 6)
    mips: list[list[np.ndarray]] = []
    done = [False] * len(inds)
    # (one buffer for the batches' reads; each is written before the next is read)
    buf = ReadBuffer()
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        ib = [int(v) for v in in_boxes[i]]
        ob = [int(v) for v in or_boxes[i]]
        in_batch = read_zarr(zarr_fullpath, input_bbox=ib, num_threads=n, out=buf)
        if in_batch.ndim == 2:
            in_batch = in_batch[np.newaxis]
        if not (save_16bit and dtype == "uint16"):
            in_batch = np.ascontiguousarray(dtypes.to_class(in_batch, "single"))
        # the box of the batch's region in the full result: the local box moved by the
        # batch's x and z starts (its y is the batch's own rows)
        bs = ob_boxes[i]
        dsr_box = [int(v) for v in ol_boxes[i] + np.array([bs[0], 0, bs[2], bs[0], 0, bs[2]])]
        # the kernels on MATLAB's layout: the input's pages turned in place, the result kept
        out_zxy = deskew_rotate_volume_matlab_layout(
            _turn_pages(in_batch, n),
            skew_angle,
            dz,
            xy_pixel_size,
            reverse=reverse,
            bbox=dsr_box,
            objective_scan=False,
            resample_factor=resample_factor,
            save_16bit=save_16bit,
            interp_method=interp_method,
            num_threads=n,
        )
        del in_batch
        if save_16bit and dtype != "uint16" and out_zxy.dtype != np.uint16:
            # uint16 rounds to the nearest integer
            out_zxy = np.ascontiguousarray(dtypes.to_class(out_zxy, "uint16"))
        # written from MATLAB's layout through a (z, y, x) view
        write_zarr(_zyx_view(out_zxy), dsr_fullpath, bbox=ob, create=False, num_threads=n)
        if mbs is not None:
            mips.append(_block_mips_zxy(out_zxy, ob, mbs, out_dtype))
        del out_zxy
        done[i] = True
        _toc(t0)

    if all(done):
        t = time.perf_counter() - t_start
        if mbs is None:
            from petakit5d.core.ids import get_uuid

            tmp = f"{flag_fullname[:-4]}_{get_uuid()}.mat"
            write_mat(tmp, {"t": float(t)}, version="7")
            os.replace(tmp, flag_fullname)
        else:
            save_block_mips(flag_fullname, [list(map(int, r)) for r in or_boxes], mips, mbs, t)
    return done


# ---------------------------------------------------------------------------------------------
# XR_deskewRotateZarr
# ---------------------------------------------------------------------------------------------


class LargeFileGeometry:
    """The batches of a large-file deskew/rotation in MATLAB's values (``[y x z]`` sizes,
    rows of 1-based inclusive boxes), as ``XR_deskewRotateZarr`` computes them."""

    def __init__(self, **kw: Any):
        self.__dict__.update(kw)


def large_file_geometry(
    im_size_m: Sequence[int],
    xy_pixel_size: float,
    dz: float,
    skew_angle: float,
    batch_size_m: Sequence[int],
    block_size_m: Sequence[int],
    resample_m: list[float] | None,
    input_bbox_m: Sequence[int] | None,
    output_bbox_m: Sequence[int] | None,
) -> LargeFileGeometry:
    """``XR_deskewRotateZarr``'s geometry in MATLAB's numbers: ``im_size_m`` the input's
    ``[ny nx nz]``, sizes ``[y x z]``, boxes MATLAB's, ``resample_m`` MATLAB's vector."""
    from petakit5d.core.blocks import zarr_chunk_coordinates_matlab
    from petakit5d.deskew_rotate.geometry import _cos, _sin

    rs = _rs3(resample_m)
    bim = [int(v) for v in im_size_m]
    if input_bbox_m is not None and len(input_bbox_m):
        ib = [int(v) for v in input_bbox_m]
        wd_start = ib[:3]
        im_size = [ib[k + 3] - ib[k] + 1 for k in range(3)]
    else:
        wd_start = [1, 1, 1]
        im_size = list(bim)
    ny, nx, nz = im_size
    theta = float(skew_angle) * math.pi / 180
    st = _sin(abs(theta))
    z_aniso = st * float(dz) / float(xy_pixel_size)
    dsr_out = [
        _mround(ny),
        _mround((nx - 1) * _cos(theta) + (nz - 1) * z_aniso / st),
        _mround((nx - 1) * st - 4),
    ]
    dsr_out = [int(_mround(dsr_out[k] / rs[k])) for k in range(3)]
    out_size = list(dsr_out)
    if output_bbox_m is None or not len(output_bbox_m):
        out_bbox = [1, 1, 1, *dsr_out]
    else:
        out_bbox = [int(v) for v in output_bbox_m]
    out_size[1] = out_bbox[4] - out_bbox[1] + 1
    out_size[2] = out_bbox[5] - out_bbox[2] + 1
    bsz = [int(v) for v in batch_size_m]
    blk = [int(v) for v in block_size_m]
    out_batch = [min(o, b) for o, b in zip(out_size, bsz, strict=True)]
    out_batch[1:] = out_size[1:]
    out_block = [min(a, b) for a, b in zip(im_size, blk, strict=True)]
    out_block = [min(a, b) for a, b in zip(out_batch, out_block, strict=True)]
    batch_b, region_b, local_b = zarr_chunk_coordinates_matlab(
        dsr_out,
        batch_size=out_batch,
        block_size=out_block,
        bbox=out_bbox,
        same_batch_size=False,
        border_size=[1, 0, 0],
    )
    in_b = batch_b.astype(np.int64).copy()
    first = np.array([_mround((v - 1) * rs[0]) for v in batch_b[:, 0]], dtype=np.int64)
    in_b[:, 0] = first + wd_start[0]
    in_b[:, 1] = wd_start[1]
    in_b[:, 2] = wd_start[2]
    rows = np.array(
        [
            _mround((e - s + 1) * rs[0] - 1)
            for s, e in zip(batch_b[:, 0], batch_b[:, 3], strict=True)
        ],
        dtype=np.int64,
    )
    in_b[:, 3] = np.minimum(in_b[:, 0] + rows, bim[0])
    in_b[:, 4] = in_b[:, 1] + im_size[1] - 1
    in_b[:, 5] = in_b[:, 2] + im_size[2] - 1
    return LargeFileGeometry(
        im_size=im_size,
        dsr_out_size=dsr_out,
        out_size=out_size,
        out_bbox=out_bbox,
        out_batch_size=out_batch,
        out_block_size=out_block,
        in_batch_bboxes=in_b,
        out_batch_bboxes=batch_b.astype(np.int64),
        out_region_bboxes=region_b.astype(np.int64),
        out_local_bboxes=local_b.astype(np.int64),
    )


@register(
    "XR_deskewRotateZarr",
    Param("frameFullpath", "@(x) ischar(x) || iscell(x)", "Input zarr (or a cell of one)."),
    Param("xyPixelSize", "@isscalar", "Pixel size in um."),
    Param("dz", "@isscalar", "Scan interval in um."),
    Param("resultDirStr", "@ischar", "Result folder under the input's folder."),
    Param("objectiveScan", "@islogical", "Objective scan (not supported: MATLAB's error)."),
    Param("overwrite", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("crop", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("skewAngle", "@isscalar", "Skew angle in degrees."),
    Param("reverse", "@islogical", "Reverse scan direction."),
    Param("DSRCombined", "@(x) islogical(x)", "Accepted; not used (always combined)."),
    Param("flipZstack", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("save16bit", "@islogical", "uint16 result (else single).", convert="same"),
    Param("saveMIP", "@islogical", "Save the MIPs of the result (MIPs/)."),
    Param("saveZarr", "@islogical", "Accepted; not used (the result is zarr)."),
    Param(
        "batchSize",
        "@isvector",
        "Batch size (z, y, x); a batch is all of x and z.",
        convert="size_yxz",
    ),  # fmt: skip
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of the result.", convert="size_yxz"),
    Param(
        "inputBbox",
        "@(x) isempty(x) || isvector(x)",
        "Box (z0, y0, x0, z1, y1, x1) of the input to deskew and rotate.",
        convert="bbox",
    ),
    Param("taskSize", "@isnumeric", "Batches per task (empty: MATLAB's rule).", convert="same"),
    Param(
        "resampleFactor",
        "@(x) isempty(x) || isnumeric(x)",
        "Resampling after rotation (z, y, x).",
        convert="resample_xyz",
    ),
    Param("interpMethod", _INTERP, "Interpolation: 'linear' or 'cubic'."),
    Param("maskFullpaths", "@iscell", "MIP masks (xy, xz, yz) that set the boxes."),
    Param("parseCluster", "@islogical", "Run the blocks as Slurm jobs."),
    Param("parseParfor", "@islogical", "Local process pool (MATLAB's error without a cluster)."),
    Param("masterCompute", "@islogical", "This process computes blocks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a block's job."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="D04",
)
def deskew_rotate_zarr(
    frame_fullpath: str | Sequence[str],
    xy_pixel_size: float,
    dz: float,
    *,
    result_dir_str: str = "DSR/",
    objective_scan: bool = False,
    overwrite: bool = False,
    crop: bool = False,
    skew_angle: float = 32.45,
    reverse: bool = False,
    dsr_combined: bool = True,
    flip_zstack: bool = False,
    save_16bit: bool = True,
    save_mip: bool = True,
    save_zarr: bool = False,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    input_bbox: Sequence[int] | None = None,
    task_size: Any = None,
    resample_factor: Sequence[float] | None = None,
    interp_method: str = "linear",
    mask_fullpaths: Sequence[str] = (),
    parse_cluster: bool = True,
    parse_parfor: bool = False,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 8,
    uuid: str = "",
    debug: bool = False,
    mcc_mode: bool = False,
    config_file: str = "",
) -> list[bool] | None:
    """Deskew and rotate a large zarr file in blocks (``XR_deskewRotateZarr``): the batches of
    output rows (``batch_size``'s y; all of x and z), with a border row on each side, as
    ``deskew_rotate_block`` tasks through the job framework into
    ``<dir>/<result_dir_str>/<name>.zarr_<uuid>``, then its MIPs (``mip_zarr``, merged from the
    MIPs of the blocks) and the result moved into place. Sample scan only. Nothing is done if
    the result exists; a writer of another call of the same result is waited for, and a
    result another writer published is kept. Returns the tasks' done flags (None if nothing
    was done)."""
    from petakit5d.core.bbox import box_to_matlab, boxes_from_matlab
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.deskew_rotate.memory import estimate_deskew_rotate_memory
    from petakit5d.io import create_zarr, get_image_dtype, get_image_size, get_zarr_info
    from petakit5d.jobs import framework, writers
    from petakit5d.jobs.publish import keep_other_result, publish_result, remove_path
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.tools.mip import mip_zarr, mip_zarr_batch_size

    del overwrite, crop, dsr_combined, save_zarr, flip_zstack
    if not uuid:
        uuid = get_uuid()
    # the block processing only supports sample scan
    if _logical(objective_scan):
        raise DeskewRotateZarrError(_OBJECTIVE)
    if not isinstance(frame_fullpath, str | os.PathLike):
        frame_fullpath = list(frame_fullpath)[0]
    frame_fullpath = os.fspath(frame_fullpath)
    if not os.path.isdir(frame_fullpath):
        raise DeskewRotateZarrError(f"The input zarr file {frame_fullpath} does not exist!")
    save_16bit = _logical(save_16bit)
    save_mip = _logical(save_mip)

    data_path, fsname, _ = fileparts(frame_fullpath)
    dsr_path = f"{data_path}/{result_dir_str}"
    if not os.path.isdir(dsr_path):
        os.makedirs(dsr_path, exist_ok=True)
    dsr_fullpath = f"{dsr_path}/{fsname}.zarr"
    if os.path.isdir(dsr_fullpath):
        print("The output result exists, skip it!", flush=True)
        return None
    print(f"Start Large-file Deskew and Rotate for {frame_fullpath}...", flush=True)
    dsr_tmp = f"{dsr_fullpath}_{uuid}"
    dtype = "uint16" if save_16bit else "single"
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(data_path, job_log_dir).parse_cluster)

    t_start = time.perf_counter()
    flag_path = f"{dsr_path}/zarr_flag/{fsname}_{uuid}/"
    # (the jobs of other files make zarr_flag/ at the same time: exist_ok, not mkdir_recursive)
    os.makedirs(flag_path, exist_ok=True)
    # a writer of another call that may still make the result: wait for it; the record is
    # marked ended when this call ends but by publishing (its blocks are skipped then)
    with writers.recorded_writer(dsr_fullpath, uuid, flag_path):
        if writers.wait(dsr_fullpath, uuid):
            keep_other_result(dsr_fullpath, dsr_tmp, flag_path)
            return None

        rs_m = _resample_matlab(resample_factor)
        rs = _rs3(rs_m)
        nz_i, ny_i, nx_i = (int(v) for v in get_image_size(frame_fullpath))
        bim = [ny_i, nx_i, nz_i]

        box_m = None
        if input_bbox is not None and np.size(input_bbox):
            box_m = [int(v) for v in box_to_matlab(tuple(int(v) for v in input_bbox))]
        out_box_m = None
        masks = [os.fspath(p) for p in (mask_fullpaths or ())]
        if masks and masks[0]:
            print("Compute input and out bounding boxes with MIP masks...", flush=True)
            for p in masks:
                print(f"    {{'{p}'}}", flush=True)  # (disp of the cell column)
            border = [max(10, int(_mround(100 / r))) for r in rs]  # [y x z]
            in_py, out_py, _ = deskew_rotate_boxes_from_masks(
                masks,
                (border[2], border[0], border[1]),
                xy_pixel_size,
                dz,
                skew_angle,
                reverse,
                rs_m if rs_m is not None else [1.0, 1.0, 1.0],
                None if box_m is None else input_bbox,
            )
            box_m = [int(v) for v in box_to_matlab(in_py)]
            out_box_m = [int(v) for v in box_to_matlab(out_py)]
            print(f"Input bounding box: [{' '.join(map(str, box_m))}]", flush=True)
            print(f"Output bounding box: [{' '.join(map(str, out_box_m))}]", flush=True)

        g = large_file_geometry(
            bim,
            xy_pixel_size,
            dz,
            skew_angle,
            _yxz(batch_size),
            _yxz(block_size),
            rs_m,
            box_m,
            out_box_m,
        )
        out_size = _zyx(g.out_size)
        out_block = _zyx(g.out_block_size)

        # a temporary result of another size or chunks (an earlier call with other settings)
        if os.path.isdir(dsr_tmp):
            info = get_zarr_info(dsr_tmp)
            if tuple(info.block_size) != out_block or tuple(info.size) != out_size:
                remove_path(dsr_tmp)
                remove_path(flag_path, folder_only=True)
                os.makedirs(flag_path, exist_ok=True)
        if not os.path.isdir(dsr_tmp):
            n_chunks = math.prod(math.ceil(o / b) for o, b in zip(out_size, out_block, strict=True))
            create_zarr(
                dsr_tmp,
                data_size=out_size,
                block_size=out_block,
                dtype=dtype,
                dim_separator="/" if n_chunks > 10000 else ".",
            )

        # the tasks save the MIPs of their blocks for the MIP slabs (mip_zarr merges them
        # instead of reading the result again); the option only then, as MATLAB passes it
        mip_kw: dict[str, Any] = {}
        if save_mip:
            zi = get_zarr_info(dsr_tmp)
            mip_kw["mip_batch_size"] = [
                int(v) for v in mip_zarr_batch_size(zi.size, zi.block_size, zi.dtype)
            ]

        in_b = boxes_from_matlab(g.in_batch_bboxes)
        ob_b = boxes_from_matlab(g.out_batch_bboxes)
        or_b = boxes_from_matlab(g.out_region_bboxes)
        ol_b = boxes_from_matlab(g.out_local_bboxes)
        num_batch = len(in_b)
        if task_size is None or np.size(task_size) == 0:
            # the number of batches a job processes
            tsz = int(max(2, min(10, _mround(num_batch / 5000))))
        else:
            tsz = int(np.asarray(task_size).reshape(-1)[0])
        num_tasks = math.ceil(num_batch / tsz)
        rs_task = None if rs_m is None else resample_factor
        tasks, outputs = [], []
        for t in range(num_tasks):
            inds = list(range(t * tsz, min((t + 1) * tsz, num_batch)))
            flag = f"{flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
            outputs.append(flag)
            tasks.append(
                Task.of(
                    deskew_rotate_block,
                    inds,
                    frame_fullpath,
                    dsr_tmp,
                    flag,
                    [list(map(int, in_b[k])) for k in inds],
                    [list(map(int, ob_b[k])) for k in inds],
                    [list(map(int, or_b[k])) for k in inds],
                    [list(map(int, ol_b[k])) for k in inds],
                    _through(xy_pixel_size, "%0.20d"),
                    _through(dz, "%0.20d"),
                    overwrite=False,
                    skew_angle=_through(skew_angle, "%0.20d"),
                    reverse=_logical(reverse),
                    flip_zstack=False,
                    save_16bit=save_16bit,
                    interp_method=interp_method,
                    resample_factor=None if rs_task is None else _through_mat2str(rs_task),
                    uuid=uuid,
                    debug=_logical(debug),
                    **mip_kw,
                )
            )
        inputs = [frame_fullpath] * num_tasks

        if parse_cluster or not parse_parfor:
            # memory of the largest batch: its input with border rows, the deskew and its
            # output region, plus the process
            in_sizes = g.in_batch_bboxes[:, 3:] - g.in_batch_bboxes[:, :3] + 1
            out_sizes = g.out_region_bboxes[:, 3:] - g.out_region_bboxes[:, :3] + 1
            ind = int(np.argmax(np.prod(in_sizes, axis=1) + np.prod(out_sizes, axis=1)))
            mem = estimate_deskew_rotate_memory(
                _zyx(in_sizes[ind]),
                xy_pixel_size,
                dz,
                skew_angle=skew_angle,
                save_16bit=save_16bit,
                interp_method=interp_method,
                resample_factor=resample_factor,
                out_size=_zyx(out_sizes[ind]),
                large_file=True,
                block_size=out_block,
                data_type=get_image_dtype(frame_fullpath),
            )[0]
            done: list[bool] = []
            for _ in range(2):
                # (the blocks stop once the result exists: another writer made it)
                done = framework.run_tasks(
                    inputs,
                    outputs,
                    tasks,
                    cpus_per_task=cpus_per_task,
                    max_job_num=math.inf,
                    task_batch_num=1,
                    master_compute=master_compute,
                    parse_cluster=parse_cluster,
                    mem_allocate=mem,
                    job_log_dir=job_log_dir,
                    mcc_mode=mcc_mode,
                    config_file=config_file,
                    final_out_fullpath=dsr_fullpath,
                )
                if all(done):
                    break
        else:
            # MATLAB passes a GPUJob it never defines
            raise DeskewRotateZarrError(
                "Unrecognized function or variable 'GPUJob'.", "MATLAB:UndefinedFunction"
            )
        if not all(done):
            raise DeskewRotateZarrError("Some blocks are not finished!")

        # another writer may have published the result: keep it and its MIPs (the same
        # result), and remove this one's temporary result, which can be incomplete
        if keep_other_result(dsr_fullpath, dsr_tmp, flag_path):
            return done

        # the MIPs before the result is moved into place (the framework counts the task as
        # done once it exists), named after the final result, merged from the MIPs of the
        # blocks in the flags
        if save_mip and os.path.isdir(dsr_tmp):
            mip_dir = f"{dsr_path}/MIPs/"
            if not os.path.isdir(mip_dir):
                os.makedirs(mip_dir, exist_ok=True)
                if os.name != "nt":
                    from petakit5d.io.tiff import group_writable

                    group_writable(mip_dir)
            mip_zarr(
                dsr_tmp,
                parse_cluster=parse_cluster,
                master_compute=master_compute,
                job_log_dir=job_log_dir,
                mcc_mode=mcc_mode,
                config_file=config_file,
                block_flag_fullpaths=outputs,
            )
        # (one may also have published it during the MIPs)
        publish_result(dsr_tmp, dsr_fullpath, flag_path)
        print(f"Elapsed time is {time.perf_counter() - t_start:f} seconds.", flush=True)
    return done


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _zyx_view(out_zxy: np.ndarray) -> np.ndarray:
    """the (z, y, x) view of a C array (z, x, y); copied where the view would read as a
    Fortran array (one slice or one row), which the zarr writer refuses"""
    v = out_zxy.transpose(0, 2, 1)
    if v.flags.f_contiguous and not v.flags.c_contiguous:
        return np.ascontiguousarray(v)
    return v


def _yxz(v: Any) -> list[int]:
    """a (z, y, x) size as MATLAB's [y x z]"""
    a = [int(x) for x in np.asarray(v).reshape(-1)]
    if len(a) == 1:
        return a * 3
    return [a[1], a[2], a[0]]


def _zyx(v: Sequence[Any]) -> tuple[int, int, int]:
    """MATLAB's [y x z] as (z, y, x)"""
    a = [int(x) for x in v]
    return (a[2], a[0], a[1])


def _through(value: float, fmt: str) -> float:
    """a number after MATLAB's task string (``sprintf(fmt)``, read back)"""
    from petakit5d.core.text import through_text

    return float(through_text(float(value), fmt))


def _through_mat2str(value: Any) -> Any:
    """a vector after ``mat2str`` in a task string, in its own form (tuple)"""
    from petakit5d.core.text import through_mat2str

    v = np.asarray(value, dtype=np.float64).reshape(-1)
    return tuple(float(x) for x in np.asarray(through_mat2str(v)).reshape(-1))
