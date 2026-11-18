"""Resampling (MATLAB's ``microscopeDataProcessing/tools/resample``): files
(``XR_resampleFrame``), datasets (``XR_resample_dataset``), large zarr files in batches
(``XR_resampleSingleZarr``, ``resampleZarrBlock``) and the resample setting of the pipeline
(``XR_checkResampleSetting``).

The resizing is ``matlabcompat.resize`` (MATLAB's ``imresize3``/``imresize`` bit for bit, the
lab's ``imresize3_fast`` and ``imresize3_average``, re-exported here); max pooling is
``kernels.mip``.

Arrays are ``(z, y, x)`` C order, boxes ``(z0, y0, x0, z1, y1, x1)`` (0-based,
end-exclusive), sizes and batch, block and border sizes ``(z, y, x)``. A resample factor is a
number, ``(z, xy)`` or ``(z, y, x)``; a factor above 1 shrinks. The MATLAB-name layer converts
MATLAB's values (``[y x z]`` or ``[xy z]`` factors, ``[ymin xmin zmin ymax xmax zmax]`` boxes,
1-based batch numbers).

Sizes and arithmetic are MATLAB's: an output size is ``round(size ./ factor)`` (half away from
zero), the batches of the zarr path are placed in double arithmetic as ``resampleZarrBlock``
places them, and the factors reach the tasks through MATLAB's text: ``%0.10f`` from the
dataset wrapper to its file tasks, ``%d`` (``%e`` for a non-integer: 7 digits) from
``XR_resampleSingleZarr`` to its block tasks.

Files keep MATLAB's names and contents: ``<data>/<result_dir_name>/<name>.tif`` or ``.zarr``,
``parameters.mat`` and ``parameters.json`` of the dataset call, the flags
``zarr_flag/<name>_<uuid>/blocks_<i>_<j>.mat`` of the block tasks. Zarr outputs are native
arrays. Results are published by one rename of a temporary result (``<name>_<uuid>.tif``,
``<name>_<uuid>.zarr``, ``<name>.zarr_<uuid>`` for large files); a result another writer
published is kept.

MATLAB quirks kept (each has a test or an oracle case):

- ``resample_file`` reads the image as single and resamples a volume with ``imresize3_fast``,
  an image (one plane) with ``imresize`` (``'linear'`` as ``'bilinear'``); its method switch is
  case-sensitive and an unknown method (``'box'``, ``'Linear'``) saves the image unresampled;
- an image takes a factor's MATLAB values ``[y x]``: of a ``(z, y, x)`` factor its y and x
  values, of a ``(z, xy)`` factor (MATLAB's ``[xy z]``, read there as ``[y x]``) the xy value
  for y and the z value for x;
- ``'max'`` and ``'mean'`` give ``ceil(size ./ factor)`` (edge blocks of the voxels inside);
  ``'max'`` truncates a non-integer factor (the mex), ``'mean'`` errors (``imresize3_average``
  needs integers), and both error for an image;
- the zarr path resizes in the data's class (rounded after every pass), the file path in
  single, so the two differ by rounding for integer data; its ``'mean'`` is rounded to the
  data's class when written;
- the zarr path's result has ``round(size ./ factor)`` voxels for every method, so ``'max'``
  and ``'mean'`` drop the partial edge blocks where that is less than ``ceil``;
- a zarr batch is resized with the scale ``outSize / size(batch)``: its voxels equal the whole
  image resized at once only where that is the image's scale and the border (rounded up to a
  multiple of the factor) holds the kernel's reach (shrinking: the factor for linear, twice
  it for cubic). Integer factors that divide the image do (borders 4 and 5 with factor 2);
  factors that don't divide it, non-integer ones and too small borders don't;
- ``resample_zarr_block``'s method switch is case-sensitive: ``'Linear'``, which
  ``resample_zarr``'s check takes, makes its tasks fail (MATLAB's undefined ``out_batch``)
  unless every factor is 1;
- with a box, a zarr batch's end is cut at the image's end before the box's start is added, so
  the last batches may read and resample voxels past the box's end (and keep them where they
  fall into the output);
- the dataset wrapper doesn't pass ``block_size`` to its file tasks (zarr results have chunks
  of at most 256) and gives large files ``.tif`` names unless ``save_zarr``.
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.matlabcompat.resize import imresize3_average, imresize3_fast
from petakit5d.params import Param, register

__all__ = [
    "ResampleError",
    "check_resample_setting",
    "imresize3_average",
    "imresize3_fast",
    "resample_dataset",
    "resample_file",
    "resample_zarr",
    "resample_zarr_block",
]

_DIMS = "MATLAB:sizeDimensionsMustMatch"
_DIMS_MSG = "Arrays have incompatible sizes for this operation."


class ResampleError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _floats(v: Any) -> list[float]:
    return [float(x) for x in np.asarray(v, dtype=np.float64).ravel()]


def _ints(v: Any) -> list[int]:
    out = []
    for x in _floats(v):
        if not math.isfinite(x) or x != math.floor(x):
            raise ValueError(f"expected integers, got {v!r}")
        out.append(int(x))
    return out


def _empty(v: Any) -> bool:
    return v is None or np.asarray(v).size == 0


def _matlab_factor(rf: Any) -> list[float]:
    """a factor's MATLAB values: ``r``, ``(z, xy)`` -> ``[xy z]``, ``(z, y, x)`` -> ``[y x z]``"""
    r = _floats(rf)
    if len(r) == 2:
        return [r[1], r[0]]
    if len(r) == 3:
        return [r[1], r[2], r[0]]
    return r


def _zyx_factor(rf: Any) -> np.ndarray:
    """a factor as ``(z, y, x)`` (a number for every axis, ``(z, xy)`` -> ``(z, xy, xy)``)"""
    r = _floats(rf)
    if len(r) == 1:
        r = r * 3
    elif len(r) == 2:
        r = [r[0], r[1], r[1]]
    return np.asarray(r, dtype=np.float64)


def _mround(v: Any) -> np.ndarray:
    from petakit5d.core.dtypes import mround

    return np.asarray(mround(np.asarray(v, dtype=np.float64)), dtype=np.float64)


def _toc(t0: float) -> None:
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _crop3d(a: np.ndarray, mbox: Sequence[float]) -> np.ndarray:
    """``crop3d_mex(a, mbox)`` with the indexing fall-back of its callers for the classes the
    mex doesn't take (``crop:dataTypeError``); other mex errors are raised"""
    from petakit5d.kernels.crop import crop3d_matlab

    try:
        return np.asarray(crop3d_matlab(a, list(mbox), mex_only=True))
    except Exception as err:  # the mex's errors carry its identifier
        if getattr(err, "identifier", "") != "crop:dataTypeError":
            raise
    return np.asarray(crop3d_matlab(a, list(mbox)))


def _max_pooling(a: np.ndarray, rs_matlab: Sequence[float]) -> np.ndarray:
    """MATLAB's ``max_pooling_3d(a, rs)`` with MATLAB's ``rs`` values"""
    from petakit5d.kernels.mip import max_pooling_3d_matlab

    return np.asarray(max_pooling_3d_matlab(a, np.asarray(rs_matlab, dtype=np.float64)))


def _average(a: np.ndarray, rs_matlab: Sequence[float]) -> np.ndarray:
    """MATLAB's ``imresize3_average(a, rs)`` with MATLAB's ``rs`` values, with its errors: an
    image's two values (``rem(sz, rs)``), a padding that is not an integer (padarray), one
    value or a factor that is not an integer (reshape)"""
    rs = list(rs_matlab)
    if len(rs) not in (1, 3):
        raise ResampleError(_DIMS_MSG, _DIMS)
    v = a.reshape((1,) * (3 - a.ndim) + a.shape) if a.ndim < 3 else a
    sz = [v.shape[1], v.shape[2], v.shape[0]]
    rr = rs * 3 if len(rs) == 1 else rs
    rem = [math.fmod(s, r) if r != 0 else math.nan for s, r in zip(sz, rr, strict=True)]
    if any(x != 0 for x in rem):
        pad = [
            math.ceil(s / r) * r - s if math.isfinite(s / r) else math.nan
            for s, r in zip(sz, rr, strict=True)
        ]
        if any(not math.isfinite(p) or p != math.floor(p) for p in pad):
            raise ResampleError(
                "Expected input number 2, PADSIZE, to be integer-valued.",
                "MATLAB:padarray:expectedInteger",
            )
    if len(rs) == 1:
        # (reshape indexes rs(2))
        raise ResampleError(
            "Index exceeds the number of array elements. Index must not exceed 1.",
            "MATLAB:badsubscript",
        )
    if any(not math.isfinite(r) or r < 1 or r != math.floor(r) for r in rs):
        raise ResampleError(
            "Size arguments must be real integers.", "MATLAB:getReshapeDims:notRealInt"
        )
    return np.asarray(imresize3_average(a, (int(rs[2]), int(rs[0]), int(rs[1]))))


# ---------------------------------------------------------------------------------------------
# XR_checkResampleSetting
# ---------------------------------------------------------------------------------------------


def check_resample_setting(
    resample_type: str,
    resample: Any,
    objective_scan: bool,
    skew_angle: float,
    xy_pixel_size: float,
    dz: float,
) -> tuple[Any, float]:
    """``XR_checkResampleSetting``: the resample factor of the deskew and rotation (MATLAB's
    ``[x y z]``) and ``zAniso`` (``dz / xy_pixel_size``, times ``sind(skew_angle)`` unless
    ``objective_scan``). ``'given'`` completes one value to three and two to ``[r1 r1 r2]``
    (an error without a factor); ``'isotropic'`` is ``[1 1 1]``; ``'xy_isotropic'`` is
    ``[1 1 zf]``, ``zf`` the ratio that makes the voxels of the rotation isotropic; another
    type returns the factor unchanged. (The pipeline's implementation,
    ``pipeline.automatic.check_resample_setting``.)"""
    from petakit5d.pipeline.automatic import check_resample_setting as check

    return check(resample_type, resample, objective_scan, skew_angle, xy_pixel_size, dz)


# ---------------------------------------------------------------------------------------------
# XR_resampleFrame
# ---------------------------------------------------------------------------------------------


def _read_frame(fn: str, box: list[int] | None, zarr_file: bool) -> np.ndarray:
    """readzarr with the box, or readtiff of the box's slices and crop3d of its rows and
    columns (MATLAB's crop3d_mex, with the indexing fall-back for other classes)"""
    from petakit5d.core.bbox import box_to_matlab
    from petakit5d.io import read_tiff, read_zarr

    if zarr_file:
        return read_zarr(fn, input_bbox=box)
    if box is None:
        return read_tiff(fn)
    im = read_tiff(fn, range=(box[0], box[3]))
    m = [float(v) for v in box_to_matlab(box)]
    return _crop3d(im, [m[0], m[1], 1, m[3], m[4], m[5] - m[2] + 1])


def _resample_frame(im: np.ndarray, resample_factor: Any, interp_method: str) -> np.ndarray:
    """XR_resampleFrame's resampling of a single array (``(z, y, x)`` or an image ``(y, x)``)"""
    from petakit5d.matlabcompat.resize import imresize

    rs = _matlab_factor(resample_factor)
    if not rs:
        raise ResampleError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    if im.ndim == 2:
        # a [y x z] factor uses its y and x factors
        if len(rs) == 3:
            rs = rs[:2]
        rs = [rs[0]] * max(0, 3 - len(rs)) + rs[1:]
        msize = [im.shape[0], im.shape[1]]
    else:
        # completed to 3D
        rs = [rs[0]] * max(0, 4 - len(rs)) + rs[1:]
        msize = [im.shape[1], im.shape[2], im.shape[0]]
    if len(rs) != len(msize):
        raise ResampleError(_DIMS_MSG, _DIMS)
    out_size = [float(v) for v in _mround(np.asarray(msize, np.float64) / np.asarray(rs))]

    # (MATLAB's switch: case-sensitive; other methods leave the image)
    if interp_method in ("nearest", "linear", "cubic"):
        if im.ndim == 2:
            method = "bilinear" if interp_method == "linear" else interp_method
            return imresize(im, out_size, method)
        return imresize3_fast(im, (out_size[2], out_size[0], out_size[1]), interp_method)
    if interp_method == "max":
        return _max_pooling(im, rs)
    if interp_method == "mean":
        return _average(im, rs)
    return im


@register(
    "XR_resampleFrame",
    Param("fn", "@ischar", "Input tiff or zarr file."),
    Param("fnout", "@ischar", "Result file (.tif, or .zarr with save_zarr)."),
    Param(
        "resampleFactor",
        "@isnumeric",
        "Factor (z, y, x), (z, xy) or one number; above 1 shrinks.",
        convert="resample_yxz",
    ),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) to resample.", convert="bbox"),
    Param(
        "interpMethod",
        "@ischar",
        "'linear', 'cubic', 'nearest', 'max' or 'mean' (another: not resampled).",
    ),
    Param("save16bit", "@islogical", "Save as uint16 (else single).", convert="same"),
    Param("zarrFile", "@islogical", "The input is zarr."),
    Param("saveZarr", "@islogical", "Save the result as zarr."),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of a zarr result.", convert="size_yxz"),
    Param("uuid", "@ischar", "Id of the temporary result (empty: a new one)."),
    task="T04",
)
def resample_file(
    fn: str,
    fnout: str,
    resample_factor: Any,
    *,
    input_bbox: Sequence[int] | None = None,
    interp_method: str = "linear",
    save_16bit: bool = True,
    zarr_file: bool = False,
    save_zarr: bool = False,
    block_size: Sequence[int] = (256, 256, 256),
    uuid: str = "",
) -> None:
    """Resample one tiff or zarr file (``XR_resampleFrame``) into ``fnout``: the image (or its
    box) in single, resized to ``round(size ./ factor)`` (``imresize3_fast``, or ``imresize``
    for one plane) or pooled by the factor (``'max'``, ``'mean'``: ``ceil(size ./ factor)``),
    saved as uint16 (``save_16bit``) or single, as tiff or zarr (chunks ``block_size``).
    Nothing is done if the result exists."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.io import write_tiff, write_zarr
    from petakit5d.jobs.publish import publish_result
    from petakit5d.kernels.ops import convert

    if not uuid:
        uuid = get_uuid()
    fn, fnout = os.fspath(fn), os.fspath(fnout)
    if _exist(fnout):
        print(f"The result file {fnout} already exists, skip it!", end="", flush=True)
        return
    from petakit5d.core.bbox import refuse_nan_end

    refuse_nan_end(input_bbox)  # (MATLAB indexes with it)
    box = None if _empty(input_bbox) else _ints(input_bbox)
    # (MATLAB's single() and uint16(), threaded)
    im = convert(_read_frame(fn, box, zarr_file), "single")
    im = _resample_frame(im, resample_factor, interp_method)
    if save_16bit:
        im = convert(im, "uint16")

    pstr, fsname, ext = fileparts(fnout)
    tmp = f"{pstr}/{fsname}_{uuid}{ext}"
    if save_zarr:
        write_zarr(im, tmp, block_size=_ints(block_size))
    else:
        write_tiff(im, tmp)
    publish_result(tmp, fnout)


# ---------------------------------------------------------------------------------------------
# resampleZarrBlock
# ---------------------------------------------------------------------------------------------


def batch_place(
    sub: Sequence[int],
    batch_size: Sequence[float],
    resample_factor: Sequence[float],
    border_size: Sequence[float],
    im_size: Sequence[int],
    out_size: Sequence[int],
    box_start: Sequence[int] | None = None,
) -> dict[str, list[int]]:
    """Where ``resampleZarrBlock`` places the batch ``sub`` (``(z, y, x)``, 0-based): MATLAB's
    1-based values per axis ``(z, y, x)``, in its double arithmetic. ``border_size`` is the
    border already rounded to the factor; ``box_start`` the box's 0-based start (None: none).

    Returns ``in_start``/``in_end`` (the input read, with the border), ``in_start_orig`` (the
    batch's start without it), ``out_start``/``out_end`` (the region of the result) and
    ``ba_start``/``ba_end`` (that region in the batch's resized output)."""
    s = np.asarray(sub, dtype=np.float64) + 1
    bs = np.asarray(batch_size, dtype=np.float64)
    rf = np.asarray(resample_factor, dtype=np.float64)
    bd = np.asarray(border_size, dtype=np.float64)
    sz = np.asarray(im_size, dtype=np.float64)
    osz = np.asarray(out_size, dtype=np.float64)
    ib_start_orig = _mround((s - 1) * bs * rf) + 1
    ib_end_orig = np.minimum(_mround(s * bs * rf), sz)
    ob_start = _mround((ib_start_orig - 1) / rf) + 1
    ob_start = _mround((ob_start - 1) / bs) * bs + 1
    if box_start is not None:
        off = np.asarray(box_start, dtype=np.float64)
        ib_start_orig = ib_start_orig + off
        ib_end_orig = ib_end_orig + off
    ib_start = np.maximum(ib_start_orig - bd, 1)
    ib_end = np.minimum(ib_end_orig + bd, sz)
    ob_end = np.minimum(ob_start + bs - 1, osz)
    ba_start = _mround((ib_start_orig - ib_start) / rf) + 1
    ba_end = ba_start + (ob_end - ob_start + 1) - 1

    def ints(v: np.ndarray) -> list[int]:
        return [int(x) for x in v]

    return {
        "in_start": ints(ib_start),
        "in_end": ints(ib_end),
        "in_start_orig": ints(ib_start_orig),
        "out_start": ints(ob_start),
        "out_end": ints(ob_end),
        "ba_start": ints(ba_start),
        "ba_end": ints(ba_end),
    }


def _batch_sub(bi: int, n_sub: Sequence[int]) -> tuple[int, int, int]:
    """MATLAB's ``ind2sub([ny nx nz], bi + 1)`` (y fastest), 0-based, as ``(z, y, x)``"""
    nz, ny, nx = n_sub
    del nz
    return bi // (ny * nx), bi % ny, (bi // ny) % nx


@register(
    "resampleZarrBlock",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based).", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("dsFullpath", "@(x) ischar(x)", "Result zarr (made before), written into."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "resampleFactor",
        "@isnumeric",
        "Factor (z, y, x) or one number; above 1 shrinks.",
        convert="resample_yxz",
    ),
    Param(
        "inputBbox",
        "@isnumeric",
        "Box (z0, y0, x0, z1, y1, x1) of the input that was resampled.",
        convert="bbox",
    ),
    Param(
        "batchSize",
        "@isnumeric",
        "Batch size (z, y, x) in the result (None: twice block_size).",
        convert="size_yxz",
    ),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the result.", convert="size_yxz"),
    Param(
        "borderSize",
        "@isnumeric",
        "Border (z, y, x) read around a batch, rounded up to a multiple of the factor.",
        convert="size_yxz",
    ),
    Param("overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("interpMethod", "@ischar", "'linear', 'cubic', 'nearest', 'max' or 'mean'."),
    task="T04",
)
def resample_zarr_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    ds_fullpath: str,
    flag_fullname: str,
    resample_factor: Any,
    *,
    input_bbox: Sequence[int] | None = None,
    batch_size: Sequence[int] | None = None,
    block_size: Sequence[int] | None = None,
    border_size: Sequence[float] | None = None,
    overwrite: bool = False,
    interp_method: str = "linear",
) -> list[bool] | None:
    """``resampleZarrBlock``: each batch of the result zarr ``ds_fullpath`` (``batch_size``,
    numbered with y fastest, then x, then z), resampled from the input region it maps to plus
    the border, in the data's class, then the flag. Nothing is done if the flag exists (unless
    ``overwrite``)."""
    from petakit5d.io import ReadBuffer, get_image_size, read_zarr, write_zarr

    # (z, y, x) for the arithmetic (a number for every axis); MATLAB's values for the pooling
    rs_m = _matlab_factor(resample_factor)
    if len(rs_m) not in (1, 3):
        raise ResampleError(_DIMS_MSG, _DIMS)
    rf = _zyx_factor(resample_factor)
    border = np.zeros(3) if _empty(border_size) else np.asarray(_floats(border_size))
    # a border that is a multiple of the resample factor keeps each batch on the grid of the
    # whole output
    border = _mround(np.ceil(border / rf) * rf)

    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise ResampleError(f"The block directory {flag_path} does not exist, skip the processing!")
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
    if not os.path.isdir(zarr_fullpath):
        raise ResampleError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    if not os.path.isdir(ds_fullpath):
        raise ResampleError(f"The output zarr file {ds_fullpath} doesnot exist!")

    sz = [int(v) for v in get_image_size(zarr_fullpath)[:3]]
    o_sz = [int(v) for v in get_image_size(ds_fullpath)[:3]]
    batch = (
        np.asarray(_floats(block_size)) * 2
        if _empty(batch_size)
        else np.asarray(_floats(batch_size))
    )
    n_sub = [int(v) for v in np.ceil(np.asarray(o_sz) / batch)]
    box_start = None if _empty(input_bbox) else _ints(input_bbox)[:3]
    every_one = bool(np.all(rf == 1))

    done = [False] * len(inds)
    buf = ReadBuffer()
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1}... ", end="", flush=True)
        t0 = time.perf_counter()
        p = batch_place(_batch_sub(bi, n_sub), batch, rf, border, sz, o_sz, box_start)
        in_box = [v - 1 for v in p["in_start"]] + p["in_end"]
        in_batch = read_zarr(zarr_fullpath, input_bbox=in_box, out=buf)
        if in_batch.ndim != 3:
            # (MATLAB's size of one plane has two values)
            raise ResampleError(_DIMS_MSG, _DIMS)
        # MATLAB's region of the resized batch, [ymin xmin zmin ymax xmax zmax]
        ba_s, ba_e = p["ba_start"], p["ba_end"]
        region = [ba_s[1], ba_s[2], ba_s[0], ba_e[1], ba_e[2], ba_e[0]]
        if every_one:
            out = _crop3d(in_batch, region)
        elif interp_method in ("nearest", "linear", "cubic"):
            # (the whole batch, then its region: computing the region alone would copy the
            # input region it needs, which costs more than the border's voxels)
            fit = np.maximum(_mround(np.asarray(in_batch.shape) / rf), np.asarray(ba_e))
            out = imresize3_fast(in_batch, [float(v) for v in fit], interp_method)
            out = _crop3d(out, region)
        elif interp_method == "max":
            out = _crop3d(_max_pooling(in_batch, rs_m), region)
        elif interp_method == "mean":
            out = _crop3d(_average(in_batch, rs_m), region)
        else:
            raise ResampleError(
                "Unrecognized function or variable 'out_batch'.", "MATLAB:UndefinedFunction"
            )
        del in_batch
        out_box = [v - 1 for v in p["out_start"]] + p["out_end"]
        write_zarr(out, ds_fullpath, bbox=out_box, create=False)
        del out
        done[i] = True
        _toc(t0)

    if all(done):
        open(flag_fullname, "w").close()
    return done


# ---------------------------------------------------------------------------------------------
# XR_resampleSingleZarr
# ---------------------------------------------------------------------------------------------


@register(
    "XR_resampleSingleZarr",
    Param("zarrFullpath", "@ischar", "Input zarr."),
    Param("rsFullpath", "@ischar", "Result zarr."),
    Param(
        "resampleFactor",
        "@isnumeric",
        "Factor (z, y, x), (z, xy) or one number; above 1 shrinks.",
        convert="resample_yxz",
    ),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) to resample.", convert="bbox"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the result.", convert="size_yxz"),
    Param(
        "batchSize",
        "@isnumeric",
        "Batch size (z, y, x) of a block task, in the result.",
        convert="size_yxz",
    ),
    Param("borderSize", "@isnumeric", "Border (z, y, x) read around a batch.", convert="size_yxz"),
    Param(
        "interpMethod",
        "@(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest', 'max', 'mean'}))",
        "'linear', 'cubic', 'nearest', 'max' or 'mean'.",
    ),
    Param("parseCluster", "@islogical", "Run the blocks as Slurm jobs."),
    Param("cpusPerTask", "@isscalar", "CPUs of a block's job."),
    Param("uuid", "@ischar", "Id of the temporary result (empty: a new one)."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T04",
)
def resample_zarr(
    zarr_fullpath: str,
    rs_fullpath: str,
    resample_factor: Any,
    *,
    input_bbox: Sequence[int] | None = None,
    block_size: Sequence[int] = (256, 256, 256),
    batch_size: Sequence[int] = (512, 512, 512),
    border_size: Sequence[float] = (5, 5, 5),
    interp_method: str = "linear",
    parse_cluster: bool = True,
    cpus_per_task: float = 1,
    uuid: str = "",
    mcc_mode: bool = False,
    config_file: str = "",
) -> list[bool] | None:
    """Resample a large zarr file in batches (``XR_resampleSingleZarr``): the result
    (``round(size ./ factor)``, chunks ``block_size``) made in ``<result>_<uuid>`` by
    ``resample_zarr_block`` tasks of whole chunks through the job framework, then moved into
    place. Nothing is done if the result exists; a result another writer published is kept.
    Returns the tasks' done flags (None if nothing was done)."""
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive
    from petakit5d.core.text import through_text
    from petakit5d.io import create_zarr, get_zarr_info
    from petakit5d.jobs import framework
    from petakit5d.jobs.publish import publish_result, remove_path
    from petakit5d.jobs.tasks import Task

    if not uuid:
        uuid = get_uuid()
    rf = _zyx_factor(resample_factor)
    zarr_fullpath, rs_fullpath = os.fspath(zarr_fullpath), os.fspath(rs_fullpath)
    if os.path.isdir(rs_fullpath):
        print("The output result exists, skip it!", flush=True)
        return None
    rs_path, rsfsname, _ = fileparts(rs_fullpath)
    tmp_path = f"{rs_fullpath[:-1] if rs_fullpath.endswith('/') else rs_fullpath}_{uuid}"
    if not os.path.isdir(zarr_fullpath):
        raise ResampleError(f"The input file {zarr_fullpath} does not exist, pleast double check!")
    info = get_zarr_info(zarr_fullpath)
    dtype = info.dtype
    in_size = np.asarray(info.size[:3], dtype=np.float64)
    from petakit5d.core.bbox import refuse_nan_end

    refuse_nan_end(input_bbox, size=True)  # (MATLAB's batches of the size NaN)
    box = None if _empty(input_bbox) else _ints(input_bbox)
    if box is not None:
        in_size = np.asarray([box[q + 3] - box[q] for q in range(3)], dtype=np.float64)
    rs_size = _mround(in_size / rf)

    block = np.minimum(rs_size, np.asarray(_floats(block_size)))
    batch = np.minimum(rs_size, np.asarray(_floats(batch_size)))
    batch = np.ceil(batch / block) * block
    num_batch = int(np.prod(np.ceil(rs_size / batch)))

    print("Process blocks for resampled data...", flush=True)
    flag_path = f"{rs_path}/zarr_flag/{rsfsname}_{uuid}/"
    if not os.path.isdir(flag_path):
        mkdir_recursive(flag_path)
    out_size = [int(v) for v in rs_size]
    chunks = [int(v) for v in block]
    if os.path.isdir(tmp_path):
        ti = get_zarr_info(tmp_path)
        if list(ti.block_size[:3]) != chunks or list(ti.size[:3]) != out_size:
            remove_path(tmp_path)
            remove_path(flag_path)
            os.makedirs(flag_path, exist_ok=True)
    if not os.path.isdir(tmp_path):
        n_chunks = math.prod(math.ceil(o / b) for o, b in zip(out_size, chunks, strict=True))
        create_zarr(
            tmp_path,
            data_size=out_size,
            block_size=chunks,
            dtype=dtype,
            dim_separator="/" if n_chunks > 10000 else ".",
        )

    task_size = max(5, math.ceil(num_batch / 5000))
    num_tasks = math.ceil(num_batch / task_size)
    # the values as MATLAB's function strings give them to the tasks (%d: %e for fractions)
    rf_t = [float(v) for v in through_text(rf, "%d")]
    border_t = None if _empty(border_size) else [float(v) for v in through_text(border_size, "%d")]
    batch_t = [float(v) for v in through_text(batch, "%d")]
    kw: dict[str, Any] = {
        "interp_method": interp_method,
        "input_bbox": box,
        "border_size": border_t,
        "batch_size": batch_t,
        "block_size": chunks,
    }
    tasks, outputs = [], []
    for t in range(num_tasks):
        inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
        flag = f"{flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
        outputs.append(flag)
        tasks.append(Task.of(resample_zarr_block, inds, zarr_fullpath, tmp_path, flag, rf_t, **kw))
    inputs = [zarr_fullpath] * num_tasks

    byte = data_type_to_byte_number(dtype)
    mem = float(np.prod(batch)) * byte / 2**30 * (float(np.prod(rf)) + 1) * 3
    common: dict[str, Any] = {
        "final_out_fullpath": rs_fullpath,
        "cpus_per_task": cpus_per_task,
        "mem_allocate": mem,
        "parse_cluster": parse_cluster,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = framework.run_tasks(inputs, outputs, tasks, **common)
    if not all(done):
        done = framework.run_tasks(inputs, outputs, tasks, **common)
    if not all(done):
        raise ResampleError(f"{sum(not d for d in done)} / {len(done)} tasks failed!")

    # another copy of this task may have published the result: keep it (the job framework then
    # reports all blocks as done, so this copy's temporary result can be incomplete)
    publish_result(tmp_path, rs_fullpath, flag_path)
    return done


# ---------------------------------------------------------------------------------------------
# XR_resample_dataset
# ---------------------------------------------------------------------------------------------


def _save_parameters(result_path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) and parameters.json of the call (MATLAB's ``pr``)"""
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(resample_dataset, values).items()}
    write_mat(result_path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, result_path + "/parameters.json")


@register(
    "XR_resample_dataset",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param(
        "resampleFactor",
        "@isnumeric",
        "Factor (z, y, x), (z, xy) or one number; above 1 shrinks.",
        convert="resample_yxz",
    ),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to resample."),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) to resample.", convert="bbox"),
    Param(
        "interpMethod",
        "@(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest', 'max', 'mean'}))",
        "'linear', 'cubic', 'nearest', 'max' or 'mean' (integer factors only).",
    ),
    Param("save16bit", "@islogical", "Save uint16 results (else single).", convert="same"),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("largeFile", "@islogical", "Large zarr files: resample in batches (resample_zarr)."),
    Param("saveZarr", "@islogical", "Save the results as zarr."),
    Param(
        "blockSize",
        "@isnumeric",
        "Chunk size (z, y, x) of the results of large files.",
        convert="size_yxz",
    ),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param(
        "borderSize",
        "@isnumeric",
        "Border (z, y, x) of a batch of large files.",
        convert="size_yxz",
    ),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("cpusPerTask", "@isscalar", "CPUs of a task's job."),
    Param("uuid", "@ischar", "Id of the temporary results (empty: a new one)."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T04",
)
def resample_dataset(
    data_paths: str | Sequence[str],
    resample_factor: Any,
    *,
    result_dir_name: str = "resampled",
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamA_ch1", "CamB_ch0", "CamB_ch1"),
    input_bbox: Sequence[int] | None = None,
    interp_method: str = "linear",
    save_16bit: bool = True,
    zarr_file: bool = False,
    large_file: bool = False,
    save_zarr: bool = False,
    block_size: Sequence[int] = (256, 256, 256),
    batch_size: Sequence[int] = (512, 512, 512),
    border_size: Sequence[float] = (5, 5, 5),
    parse_cluster: bool = True,
    job_log_dir: str = "../job_logs",
    master_compute: bool = True,
    cpus_per_task: float = 2,
    uuid: str = "",
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Resample every file of the datasets that matches a channel pattern
    (``XR_resample_dataset``) into ``<data path>/<result_dir_name>/<name>.tif`` (``.zarr`` with
    ``save_zarr``), one ``resample_file`` task per file (``resample_zarr`` for large zarr
    files) through the job framework. parameters.mat and .json of the call go into each result
    folder."""
    values = {k: v for k, v in locals().items() if k in _DATASET_ARGS}
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import simplify_path
    from petakit5d.core.text import through_text
    from petakit5d.io import get_image_dtype, get_image_size
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_image_filenames

    if not uuid:
        uuid = get_uuid()
    paths = [data_paths] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    result_paths = []
    for d, p in enumerate(paths):
        paths[d] = simplify_path(p)
        result_path = f"{paths[d]}/{result_dir_name}/"
        result_paths.append(result_path)
        os.makedirs(result_path, exist_ok=True)
        group_writable(result_path)
        _save_parameters(result_path, values)

    _, fsns, fd_inds, filepaths, _ = parse_image_filenames(paths, zarr_file, channel_patterns)
    n_f = len(fsns)
    if n_f == 0:
        raise ResampleError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    sz = [float(v) for v in get_image_size(filepaths[0])[:3]]
    byte = data_type_to_byte_number(get_image_dtype(filepaths[0]))

    ext = ".zarr" if save_zarr else ".tif"
    results = [f"{result_paths[int(fd_inds[i])]}/{fsns[i]}{ext}" for i in range(n_f)]
    # the factors as the function strings give them (mat2str_comma(resampleFactor, 10))
    rf_t = [float(v) for v in through_text(_floats(resample_factor), "%0.10f")]
    from petakit5d.core.bbox import refuse_nan_end

    # a NaN end: MATLAB's file tasks fail where they index with it (large files: where the
    # size of the box makes their batches)
    refuse_nan_end(input_bbox, size=large_file)
    box = None if _empty(input_bbox) else _ints(input_bbox)
    rf = _floats(resample_factor)
    if large_file:
        cpus = float(through_text(float(cpus_per_task), "%d"))
        kw: dict[str, Any] = {
            "input_bbox": box,
            "block_size": _ints(block_size),
            "batch_size": _ints(batch_size),
            "border_size": None if _empty(border_size) else _floats(border_size),
            "interp_method": interp_method,
            "parse_cluster": bool(parse_cluster),
            "cpus_per_task": int(cpus) if cpus == int(cpus) else cpus,
            "uuid": uuid,
            "mcc_mode": bool(mcc_mode),
            "config_file": str(config_file),
        }
        tasks = [Task.of(resample_zarr, filepaths[i], results[i], rf_t, **kw) for i in range(n_f)]
        voxels = float(np.prod(np.asarray(_floats(batch_size))))
        mem_allocate = voxels * byte / 2**30 * (float(np.prod(rf)) + 1) * 3
    else:
        kw = {
            "input_bbox": box,
            "zarr_file": bool(zarr_file),
            "save_zarr": bool(save_zarr),
            "interp_method": interp_method,
            "save_16bit": bool(save_16bit),
            "uuid": uuid,
        }
        tasks = [Task.of(resample_file, filepaths[i], results[i], rf_t, **kw) for i in range(n_f)]
        mem_allocate = float(np.prod(sz)) * byte / 1024**3 * (2 + 2 / float(np.prod(rf)))

    framework.run_tasks(
        list(filepaths),
        results,
        tasks,
        parse_cluster=parse_cluster,
        job_log_dir=job_log_dir,
        master_compute=master_compute,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )


_DATASET_ARGS = (
    "data_paths",
    "resample_factor",
    "result_dir_name",
    "channel_patterns",
    "input_bbox",
    "interp_method",
    "save_16bit",
    "zarr_file",
    "large_file",
    "save_zarr",
    "block_size",
    "batch_size",
    "border_size",
    "parse_cluster",
    "job_log_dir",
    "master_compute",
    "cpus_per_task",
    "uuid",
    "mcc_mode",
    "config_file",
)
