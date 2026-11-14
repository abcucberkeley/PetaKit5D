"""Channel unmixing (MATLAB's ``microscopeDataProcessing/unmix``): the weighted channel sum
(``unmix_channels_sum`` and its mex), frames (``XR_unmix_channels_frame``), large zarr files in
blocks (``XR_unmix_channels_zarr``, ``unmix_channels_block``, ``unmix_channels_gaussian_block``)
and datasets (``XR_unmix_channels_data_wrapper``).

The unmixed image of channels ``x_c`` with factors ``f_c``::

    out = cast(max(0, sum_c (single(x_c) - background_c) * f_c) + offset, class)

in ``unmix_channels_sum_mex``'s arithmetic (``petakit5d._ext.unmix``): each single-double
operation in double and rounded to single once, the sum in single in channel order,
``max(0, NaN)`` and ``max(0, -0)`` are +0, the offset added only when it is not 0, and MATLAB's
``cast`` (round half away from zero, saturate, NaN to 0) for uint8 and uint16. In 'gaussian'
mode a channel with a nonzero sigma is filtered first (``imgaussfilt3`` in single, after its
background); with flat-field correction each channel is corrected first
(``XR_LSFlatFieldCorrection`` with offset 0). The sum is per element, so it runs on the port's
``(z, y, x)`` arrays as they are.

As MATLAB, the channels are summed in one pass where the channels as read and the result take
at most 8 bytes per voxel (for frames and gaussian blocks; linear blocks: at most 4 bytes more
than the largest), otherwise one at a time into a single partial sum (here in place, so no more
memory than the partial sum and one channel). Both give the same values.

Arrays are ``(z, y, x)`` C order; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive;
batch, block and border sizes ``(z, y, x)``; ``channel_ind`` and batch numbers 0-based; the
sigmas one entry per channel, a scalar or ``(z, y, x)``. The MATLAB-name layer converts
MATLAB's values (``[y x z]`` sizes and sigma rows, 1-based numbers and boxes).

Files keep MATLAB's names and contents: ``<data>/<result_dir_name>/<name>.tif`` or ``.zarr``
(the channel ``channel_ind``'s name), the flags ``zarr_flag/<name>_<uuid>/blocks_<i>_<j>.mat``
of the block tasks. Zarr results are native arrays. Results are published by one rename of a
temporary result (``<name>_<uuid>.tif`` / ``.zarr``); a result another writer published is
kept (MATLAB's movefile replaced a tiff).

Values that MATLAB passes to its tasks as text keep MATLAB's rounding: the dataset's factors
with ``%0.10f``, its sigmas with 10 significant digits (``mat2str``), its backgrounds and offset
with 15; the zarr function's factors and sigmas with 15.

MATLAB quirks kept (each has a test):

- ``unmix_channels_zarr`` writes the class of the channel it unmixes; ``save_16bit`` is not
  used there;
- ``unmix_channels_dataset`` with ``large_file`` passes neither the backgrounds, the offset nor
  the flat-field correction to the zarr tasks, and names their results ``.tif`` unless
  ``save_zarr`` is set, so that without it the job framework waits for results that never come
  (it runs each task three times, minutes in local mode);
- the linear blocks include channels with factor 0 (read and multiplied); the frames skip them.

Differences: integer-class factors, backgrounds or offsets are taken as doubles (MATLAB's code
would compute in that integer class).
"""

from __future__ import annotations

import functools
import math
import os
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "UnmixError",
    "extension_available",
    "unmix_channels_block",
    "unmix_channels_dataset",
    "unmix_channels_file",
    "unmix_channels_gaussian_block",
    "unmix_channels_sum",
    "unmix_channels_sum_ref",
    "unmix_channels_zarr",
]

# the classes of the kernel (unmix_channels_sum_mex's)
_KERNEL = ("uint8", "uint16", "single", "double")
_F32 = np.dtype(np.float32)


class UnmixError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


@functools.cache
def _ext():
    """the extension module, or None (not built: the NumPy reference runs)"""
    from petakit5d import _ext as ext

    try:
        return ext.load("unmix")
    except ImportError:
        return None


def extension_available() -> bool:
    """Whether the C++ kernel is built (otherwise the NumPy reference runs)."""
    return _ext() is not None


def _threads(num_threads: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(num_threads)


def _class(a: np.ndarray) -> str:
    from petakit5d.core.arrays import matlab_class

    return matlab_class(a)


def _dtype(cls: str) -> np.dtype:
    from petakit5d.core.arrays import numpy_dtype

    return numpy_dtype(cls)


def _bytes(cls: str) -> int:
    from petakit5d.core.arrays import data_type_to_byte_number

    return data_type_to_byte_number(cls)


def _numbers(v: Any) -> np.ndarray:
    """a MATLAB numeric value as a flat float64 array in MATLAB's order (None and [] empty)"""
    if v is None:
        return np.zeros(0)
    a = np.asarray(v)
    if a.dtype.kind == "b":
        a = a.astype(np.float64)
    return np.asarray(a, dtype=np.float64).reshape(-1, order="F")


def _empty(v: Any) -> bool:
    """MATLAB's isempty of a value (None, '', an empty array or sequence)"""
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    if isinstance(v, np.ndarray):
        return v.size == 0
    if isinstance(v, list | tuple):
        return len(v) == 0
    return False


def _flat(a: np.ndarray) -> np.ndarray:
    """a C-contiguous 1-D view (a copy only for other layouts)"""
    return np.ascontiguousarray(a).reshape(-1)


def _toc(t0: float) -> None:
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)


# ---------------------------------------------------------------------------------------------
# the sum (unmix_channels_sum, unmix_channels_sum_mex)
# ---------------------------------------------------------------------------------------------


def _cast_ref(v: np.ndarray, cls: str) -> np.ndarray:
    """the mex's cast of single values"""
    if cls == "single":
        return v.copy()
    if cls == "double":
        return v.astype(np.float64)
    hi = 65535.0 if cls == "uint16" else 255.0
    with np.errstate(invalid="ignore"):
        d = v.astype(np.float64) + 0.5
        r = np.where(v >= np.float32(0.5), np.where(v >= np.float32(hi), hi, np.floor(d)), 0.0)
    return r.astype(_dtype(cls))


def unmix_channels_sum_ref(
    ims: Sequence[np.ndarray],
    factors: Sequence[float],
    backgrounds: Sequence[float],
    offset: float,
    dtype: str,
    im: np.ndarray | None = None,
) -> np.ndarray:
    """The kernel's arithmetic in NumPy (the reference of the tests, and the path without the
    extension): ``ims`` of the kernel's classes and one shape, one double factor and
    background per channel, ``dtype`` ``''`` for the partial sum."""
    shape = ims[0].shape if len(ims) else im.shape
    with np.errstate(all="ignore"):
        acc = np.zeros(shape, _F32) if im is None else np.array(im, dtype=_F32)
        for x, f, b in zip(ims, factors, backgrounds, strict=True):
            v = np.asarray(x).astype(_F32)
            if b != 0:
                v = (v.astype(np.float64) - float(b)).astype(_F32)
            acc = acc + (v.astype(np.float64) * float(f)).astype(_F32)
        if not dtype:
            return acc
        v = np.where(acc > 0, acc, np.float32(0))
        if offset != 0:
            v = (v.astype(np.float64) + float(offset)).astype(_F32)
        return _cast_ref(v, dtype)


def _sum(
    ims: Sequence[np.ndarray],
    factors: Sequence[float],
    backgrounds: Sequence[float],
    offset: float,
    dtype: str,
    im: np.ndarray | None,
    out: np.ndarray | None,
    num_threads: int | None,
) -> np.ndarray:
    """the kernel on arrays of its classes (one shape); ``out`` may be ``im``"""
    shape = ims[0].shape if len(ims) else im.shape
    odt = _F32 if not dtype else _dtype(dtype)
    if out is None:
        out = np.empty(shape, dtype=odt)
    elif out.dtype != odt or out.shape != tuple(shape) or not out.flags.c_contiguous:
        raise ValueError(f"out must be a C-contiguous {odt} array of shape {tuple(shape)}")
    ext = _ext()
    if ext is None:
        out[...] = unmix_channels_sum_ref(ims, factors, backgrounds, offset, dtype, im)
        return out
    ext.unmix_sum(
        [_flat(x) for x in ims],
        np.ascontiguousarray(factors, dtype=np.float64).reshape(-1),
        np.ascontiguousarray(backgrounds, dtype=np.float64).reshape(-1),
        float(offset),
        bool(dtype),
        None if im is None else _flat(im),
        out.reshape(-1),
        _threads(num_threads),
    )
    return out


def _subtract(x: np.ndarray, b: float, num_threads: int | None) -> np.ndarray:
    """``x - b`` of a single array and a double, in place (MATLAB's single minus double)"""
    ext = _ext()
    if ext is None or not x.flags.c_contiguous:
        with np.errstate(all="ignore"):
            x[...] = (x.astype(np.float64) - float(b)).astype(_F32)
        return x
    f = x.reshape(-1)
    ext.subtract(f, float(b), f, _threads(num_threads))
    return x


def _single(x: np.ndarray, num_threads: int | None = None) -> np.ndarray:
    """``single(x)``: x itself for single data"""
    if x.dtype == _F32:
        return x
    from petakit5d.kernels.ops import convert

    return convert(x, "single", num_threads=num_threads)


def unmix_channels_sum(
    ims: Sequence[np.ndarray],
    unmix_factors: Any,
    const_background: Any = None,
    const_offset: Any = None,
    dtype: str = "",
    im: np.ndarray | None = None,
    *,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """``unmix_channels_sum``: the unmixed image of the channel images ``ims`` (one shape)::

        cast(max(0, im + sum_c (single(ims[c]) - const_background[c]) * unmix_factors[c])
             + const_offset, dtype)

    ``const_background`` None or one value per channel; ``const_offset`` None or a scalar, added
    when it is not 0; ``dtype`` a MATLAB class, or ``''`` for the partial sum (single, not
    clipped and without the offset); ``im`` None (0) or a single partial sum to add the channels
    to. uint8, uint16, single and double channels and results in the mex's arithmetic; other
    channel classes are made single first and other results cast from the single result, as
    MATLAB's code computes them (the same values). ``out``: an array of the result's class and
    shape to write into (it may be ``im``: the channels are added in place).
    """
    from petakit5d.core import dtypes

    arrs = [np.asarray(x) for x in ims]
    nc = len(arrs)
    if nc == 0 and im is None:
        raise UnmixError(
            "The channel images must be a non-empty cell array (or empty with a partial sum)",
            "unmix_channels:inputError",
        )
    shape = arrs[0].shape if nc else np.shape(im)
    for c, x in enumerate(arrs):
        if x.shape != shape:
            raise UnmixError(
                f"Channel {c + 1} has a different size than channel 1", "unmix_channels:inputError"
            )
        if np.iscomplexobj(x):
            raise UnmixError(
                f"Channel {c + 1}: data type not supported (uint8, uint16, single or double)",
                "unmix_channels:dataTypeError",
            )
    if im is not None:
        im = np.asarray(im)
        if im.dtype != _F32 or im.shape != shape:
            raise UnmixError(
                "The partial sum must be a real single array of the channels' size",
                "unmix_channels:dataTypeError",
            )
    f = _numbers(unmix_factors)
    if f.size < nc:
        raise UnmixError(
            f"Index exceeds the number of array elements. Index must not exceed {f.size}.",
            "MATLAB:badsubscript",
        )
    bg = np.zeros(nc)
    if not _empty(const_background):
        b = _numbers(const_background)
        if b.size < nc:
            raise UnmixError(
                f"Index exceeds the number of array elements. Index must not exceed {b.size}.",
                "MATLAB:badsubscript",
            )
        bg = b[:nc]
    offset = 0.0
    if not _empty(const_offset):
        o = _numbers(const_offset)
        if o.size != 1:
            raise UnmixError("The offset must be a scalar", "unmix_channels:inputError")
        offset = float(o[0])
    dtype = "" if dtype is None else str(dtype)
    # other channel classes: single() first, as MATLAB's code does
    arrs = [x if _class(x) in _KERNEL else np.asarray(dtypes.to_class(x, "single")) for x in arrs]
    if not dtype or dtype in _KERNEL:
        return _sum(arrs, f[:nc], bg, offset, dtype, im, out, num_threads)
    r = np.asarray(
        dtypes.cast(_sum(arrs, f[:nc], bg, offset, "single", im, None, num_threads), dtype)
    )
    if out is None:
        return r
    out[...] = r
    return out


# ---------------------------------------------------------------------------------------------
# sigmas, factors and texts
# ---------------------------------------------------------------------------------------------


def _sigma_rows(unmix_sigmas: Any) -> list[np.ndarray]:
    """one entry per channel: a scalar (isotropic) or (z, y, x)"""
    if _empty(unmix_sigmas):
        return []
    a = np.asarray(unmix_sigmas, dtype=object)
    if a.ndim == 0:
        return [np.asarray([float(a)])]
    rows = []
    for r in unmix_sigmas:
        v = np.atleast_1d(np.asarray(r, dtype=np.float64)).reshape(-1)
        rows.append(v)
    return rows


def _sigma(rows: list[np.ndarray], c: int) -> np.ndarray:
    if c >= len(rows):
        raise UnmixError(
            f"Index in position 1 exceeds array bounds. Index must not exceed {len(rows)}.",
            "MATLAB:badsubscript",
        )
    return rows[c]


def _filtered(rows: list[np.ndarray], c: int) -> bool:
    """``~all(unmixSigmas(c, :) == 0)``"""
    return not bool(np.all(_sigma(rows, c) == 0))


def _sigmas_in(v: Any) -> Any:
    """MATLAB's sigmas (one row per channel, [sy sx sz] or one value) as the Python form; a
    vector is a column (one value per channel), a NumPy matrix keeps its rows"""
    if _empty(v):
        return None
    a = np.asarray(v, dtype=np.float64)
    if a.ndim == 0:
        return (float(a),)
    if a.ndim == 1 or (a.ndim == 2 and a.shape[1] == 1 and type(v).__module__ != "numpy"):
        return tuple(float(x) for x in a.reshape(-1))
    if a.ndim == 2 and a.shape[1] == 1:
        return tuple(float(x) for x in a[:, 0])
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(f"unmixSigmas: one value or [y x z] per channel, got shape {a.shape}")
    return tuple((float(r[2]), float(r[0]), float(r[1])) for r in a)


def _sigmas_out(v: Any) -> Any:
    """the Python sigmas as MATLAB's matrix (a column, or rows [sy sx sz])"""
    if _empty(v):
        return None
    rows = _sigma_rows(v)
    if all(r.size == 1 for r in rows):
        return tuple(float(r[0]) for r in rows)
    out = []
    for r in rows:
        z, y, x = (float(r[0]),) * 3 if r.size == 1 else (float(r[0]), float(r[1]), float(r[2]))
        out.append((y, x, z))
    return tuple(out)


def _gauss_sigma(row: np.ndarray) -> Any:
    """the sigma of imgaussfilt3 (a scalar or (z, y, x))"""
    return float(row[0]) if row.size == 1 else tuple(float(v) for v in row)


def _through_text(v: Any, fmt: str) -> Any:
    """values after ``sprintf(fmt, ...)`` and MATLAB's parser (NaN is MATLAB's NaN)"""
    from petakit5d.core.text import through_text

    if _empty(v):
        return v
    a = np.asarray(v, dtype=np.float64)
    r = np.asarray(through_text(a.reshape(-1), fmt), dtype=np.float64).reshape(a.shape)
    return _like(v, r)


def _through_mat2str(v: Any, precision: int) -> Any:
    """values after ``mat2str(v, precision)`` and MATLAB's parser (NaN is MATLAB's NaN)"""
    from petakit5d.core.text import through_mat2str

    if _empty(v):
        return v
    a = np.asarray(v, dtype=np.float64)
    r = np.asarray(through_mat2str(a.reshape(-1), precision), dtype=np.float64)
    return _like(v, r.reshape(a.shape))


def _like(v: Any, r: np.ndarray) -> Any:
    """``r`` in the form of ``v`` (a number, nested tuples or an array)"""
    if type(v).__module__ == "numpy":
        return r
    if np.ndim(v) == 0:
        return float(r)
    return tuple(_like(e, re) for e, re in zip(v, r, strict=True))


def _sigmas_text(v: Any, precision: int) -> Any:
    """sigmas after ``mat2str(unmixSigmas, precision)``"""
    if _empty(v):
        return v
    return tuple(_through_mat2str(r, precision) for r in v)


# ---------------------------------------------------------------------------------------------
# frames
# ---------------------------------------------------------------------------------------------


def _read(path: str, zarr_file: bool) -> np.ndarray:
    from petakit5d.io import read_tiff, read_zarr

    return read_zarr(path) if zarr_file else read_tiff(path)


@register(
    "XR_unmix_channels_frame",
    Param("frameFullpaths", "@(x) ischar(x) || iscell(x)", "Files of the channels (tiff or zarr)."),
    Param("unmixFactors", "@(x) isnumeric(x)", "Unmixing factor of each channel."),
    Param("mode", "@ischar", "'linear', or 'gaussian' (each channel filtered with its sigma)."),
    Param(
        "unmixSigmas",
        "@isnumeric",
        "Sigmas of the gaussian mode, one per channel: a scalar or (z, y, x).",
        convert=(_sigmas_in, _sigmas_out),
    ),
    Param("resultDirName", "@ischar", "Result folder under the data path."),
    Param("channelInd", "@isnumeric", "The channel to unmix (0-based).", convert="index1"),
    Param("FFCorrection", "@islogical", "Flat-field correct the channels first."),
    Param("lowerLimit", "@isnumeric", "Lowest value of the normalized flat fields."),
    Param("FFImagePaths", "@iscell", "Flat-field image of each channel."),
    Param("backgroundPaths", "@iscell", "Background image of each channel."),
    Param("constBackground", "@isnumeric", "Background of each channel, subtracted first."),
    Param("constOffset", "@(x) isnumeric(x)", "If set, added after the unmixing.", convert="same"),
    Param("save16bit", "@islogical", "Save uint16; otherwise single.", convert="same"),
    Param("zarrFile", "@islogical", "The channels are zarr."),
    Param("saveZarr", "@islogical", "Save the result as zarr."),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of a zarr result.", convert="size_yxz"),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    task="B02",
)
def unmix_channels_file(
    frame_fullpaths: str | Sequence[str],
    unmix_factors: Any,
    *,
    mode: str = "linear",
    unmix_sigmas: Any = None,
    result_dir_name: str = "Unmixed",
    channel_ind: int = 0,
    ff_correction: bool = False,
    lower_limit: Any = 0.4,
    ff_image_paths: Sequence[str] = ("", "", ""),
    background_paths: Sequence[str] = ("", "", ""),
    const_background: Any = None,
    const_offset: Any = None,
    save_16bit: bool = True,
    zarr_file: bool = False,
    save_zarr: bool = False,
    block_size: Sequence[int] = (256, 256, 256),
    uuid: str = "",
    debug: bool = False,
) -> None:
    """Unmix one frame (``XR_unmix_channels_frame``): the channels' files (one per channel,
    tiff or zarr with ``zarr_file``) into ``<data path of channel channel_ind>/
    <result_dir_name>/<its name>.tif`` (or ``.zarr``, chunks ``block_size``), uint16 with
    ``save_16bit`` else single. Channels with factor 0 are not read. Nothing is done if the
    result exists."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, strip_trailing_sep
    from petakit5d.io import get_image_dtype, get_image_size, read_tiff, write_tiff, write_zarr
    from petakit5d.jobs.files import exist
    from petakit5d.jobs.publish import publish_result

    del debug
    if not uuid:
        uuid = get_uuid()
    paths = [frame_fullpaths] if isinstance(frame_fullpaths, str | os.PathLike) else frame_fullpaths
    paths = [strip_trailing_sep(os.fspath(p)) for p in paths]
    ci = int(channel_ind)
    data_path, fsname, _ext_name = fileparts(paths[ci])
    unmix_path = f"{data_path}/{result_dir_name}/"
    os.makedirs(unmix_path, exist_ok=True)
    if save_zarr:
        full = f"{unmix_path}/{fsname}.zarr"
        tmp = f"{unmix_path}/{fsname}_{uuid}.zarr"
    else:
        full = f"{unmix_path}/{fsname}.tif"
        tmp = f"{unmix_path}/{fsname}_{uuid}.tif"
    if exist(full):
        print("The output results exist, skip it!", flush=True)
        return
    print(f"Start channel unmixing for {fsname}...", flush=True)

    im_size = get_image_size(paths[ci])
    nc = len(paths)
    factors = _numbers(unmix_factors)
    if factors.size < nc:
        raise UnmixError(
            f"Index exceeds the number of array elements. Index must not exceed {factors.size}.",
            "MATLAB:badsubscript",
        )
    rows = _sigma_rows(unmix_sigmas)
    gaussian = mode == "gaussian"
    used = [bool(factors[c] != 0) for c in range(nc)]
    filtered = [used[c] and gaussian and _filtered(rows, c) for c in range(nc)]
    processed = [used[c] and (bool(ff_correction) or filtered[c]) for c in range(nc)]
    dtype = "uint16" if save_16bit else "single"
    bgs_given = not _empty(const_background)
    bgv = _numbers(const_background)

    def background(c: int) -> float:
        if c >= bgv.size:
            raise UnmixError(
                f"Index exceeds the number of array elements. Index must not exceed {bgv.size}.",
                "MATLAB:badsubscript",
            )
        return float(bgv[c])

    # the channels as read and the result in one pass when they take at most 8 bytes per voxel
    # (MATLAB's rule); otherwise one at a time into the partial sum, in place
    at_once = any(used) and not any(processed) and sum(used) + _bytes(dtype) <= 8
    if at_once:
        dts = [get_image_dtype(paths[c], zarr_file) for c in range(nc) if used[c]]
        at_once = (
            all(d in _KERNEL for d in dts) and sum(_bytes(d) for d in dts) + _bytes(dtype) <= 8
        )

    ims: list[np.ndarray] = []
    fs: list[float] = []
    bgs: list[float] = []
    im: np.ndarray | None = None
    for c in range(nc):
        if not used[c]:
            continue
        im_c = _read(paths[c], zarr_file)
        bg_c = 0.0
        if processed[c] or im_c.dtype == _F32:
            im_c = _single(im_c)
            if ff_correction:
                from petakit5d.flatfield import flatfield_correct_volume

                ff_image, bk_image = ff_image_paths[c], background_paths[c]
                print(
                    f"Flat-field correction for frame {paths[c]}...\n"
                    f"  Flat-field image: {ff_image}\n  Background image: {bk_image}",
                    flush=True,
                )
                im_c = flatfield_correct_volume(
                    im_c,
                    read_tiff(ff_image),
                    read_tiff(bk_image),
                    lower_limit=lower_limit,
                    const_offset=0,
                )
            if bgs_given:
                im_c = _subtract(np.ascontiguousarray(im_c), background(c), None)
            if filtered[c]:
                from petakit5d.matlabcompat.filters import imgaussfilt3

                im_c = imgaussfilt3(im_c, _gauss_sigma(_sigma(rows, c)))
        elif bgs_given:
            bg_c = background(c)
        if at_once:
            ims.append(im_c)
            fs.append(float(factors[c]))
            bgs.append(bg_c)
        else:
            im = _sum([im_c], [factors[c]], [bg_c], 0.0, "", im, im, None)
        del im_c

    if at_once:
        im = unmix_channels_sum(ims, fs, bgs, const_offset, dtype)
        del ims
    else:
        if im is None:
            # no channel contributes
            im = np.zeros(tuple(int(v) for v in im_size[:3]), dtype=_F32)
        im = unmix_channels_sum([], [], None, const_offset, dtype, im)

    if save_zarr:
        write_zarr(im, tmp, block_size=[int(v) for v in block_size], create=True)
    else:
        write_tiff(im, tmp)
    if exist(tmp):
        publish_result(tmp, full)


# ---------------------------------------------------------------------------------------------
# zarr blocks
# ---------------------------------------------------------------------------------------------


def _block_start(
    batch_inds: Sequence[int], flag_fullname: str, overwrite: bool, zarr_fullpaths: Sequence[str]
) -> list[int] | None:
    """the checks of the block functions: the batch numbers, or None if the flag exists"""
    from petakit5d.jobs.files import exist

    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise UnmixError(f"The block directory {flag_path} does not exist, skip the processing!")
    inds = [int(v) for v in np.asarray(batch_inds).reshape(-1)]
    if exist(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            print(
                f"The batch files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!",
                flush=True,
            )
            return None
    if not os.path.isdir(zarr_fullpaths[0]):
        raise UnmixError(f"The input zarr file {zarr_fullpaths[0]} doesnot exist!")
    return inds


def _check_outputs(factors: np.ndarray, nc: int, unmix_fullpath: str) -> str:
    from petakit5d.io import get_image_dtype

    if factors.size != nc:
        raise UnmixError("The number of unmix factors must equal to the number of channesl!")
    if not os.path.isdir(unmix_fullpath):
        raise UnmixError(f"The output zarr file {unmix_fullpath} does not exist!")
    return get_image_dtype(unmix_fullpath, True)


class _Buffers:
    """read buffers, one per channel, and one result buffer, reused over the batches"""

    def __init__(self, n: int):
        from petakit5d.io import ReadBuffer

        self.reads = [ReadBuffer() for _ in range(n)]
        self.result = ReadBuffer()

    def read(self, c: int, path: str, box: Sequence[int]) -> np.ndarray:
        from petakit5d.io import read_zarr

        a = read_zarr(path, input_bbox=list(box), out=self.reads[c])
        return a[np.newaxis] if a.ndim == 2 else a

    def out(self, shape: Sequence[int], dtype: str) -> np.ndarray:
        return self.result.view(shape, _dtype(dtype))


def _box_shape(box: Sequence[int]) -> tuple[int, int, int]:
    return (int(box[3] - box[0]), int(box[4] - box[1]), int(box[5] - box[2]))


@register(
    "unmix_channels_block",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based).", convert="index1"),
    Param("zarrFullpaths", "@(x) iscell(x)", "Zarr files of the channels."),
    Param("unmixFullpath", "@(x) ischar(x)", "Result zarr (made before), written into."),
    Param("unmixFactors", "@(x) isnumeric(x)", "Unmixing factor of each channel."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch.",
        convert="bboxes",
        py="batch_bboxes",
    ),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("uuid", "@ischar", "Id of the call."),
    Param("debug", "@islogical", "Not used."),
    task="B02",
)
def unmix_channels_block(
    batch_inds: Sequence[int],
    zarr_fullpaths: Sequence[str],
    unmix_fullpath: str,
    unmix_factors: Any,
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    uuid: str = "",
    debug: bool = False,
) -> list[bool] | None:
    """``unmix_channels_block``: the linear unmixing of each batch (``batch_bboxes``) of the
    channels' zarr files, written into the same box of the result zarr ``unmix_fullpath`` in
    its class; then the flag. Nothing is done if the flag exists (unless ``overwrite``)."""
    from petakit5d.io import get_image_dtype, write_zarr

    del uuid, debug
    paths = [os.fspath(p) for p in zarr_fullpaths]
    inds = _block_start(batch_inds, flag_fullname, overwrite, paths)
    if inds is None:
        return None
    nc = len(paths)
    factors = _numbers(unmix_factors)
    dtype = _check_outputs(factors, nc, unmix_fullpath)
    # the channels in one pass when they and the result take no more than the single sum and
    # the largest of them (MATLAB's rule); otherwise one at a time into the sum
    dts = [get_image_dtype(p, True) for p in paths]
    at_once = all(d in _KERNEL for d in [*dts, dtype])
    if at_once:
        b = [_bytes(d) for d in [*dts, dtype]]
        at_once = sum(b) <= 4 + max(b)

    boxes = np.asarray(batch_bboxes, dtype=np.int64).reshape(-1, 6)
    buf = _Buffers(nc if at_once else 1)
    zero_bgs = [0.0] * nc
    done = [False] * len(inds)
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        box = [int(v) for v in boxes[i]]
        shape = _box_shape(box)
        out = buf.out(shape, dtype)
        if at_once:
            ims = [buf.read(c, paths[c], box) for c in range(nc)]
            _sum(ims, factors, zero_bgs, 0.0, dtype, None, out, None)
            del ims
        else:
            acc = np.zeros(shape, dtype=_F32)
            for c in range(nc):
                im_c = buf.read(0, paths[c], box)
                _sum([im_c], [factors[c]], [0.0], 0.0, "", acc, acc, None)
            _sum([], [], [], 0.0, dtype, acc, out, None)
            del acc
        write_zarr(out, unmix_fullpath, bbox=box, create=False)
        done[i] = True
        _toc(t0)
    if all(done):
        open(flag_fullname, "w").close()
    return done


@register(
    "unmix_channels_gaussian_block",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based).", convert="index1"),
    Param("zarrFullpaths", "@(x) iscell(x)", "Zarr files of the channels."),
    Param("unmixFullpath", "@(x) ischar(x)", "Result zarr (made before), written into."),
    Param("unmixFactors", "@(x) isnumeric(x)", "Unmixing factor of each channel."),
    Param(
        "unmixSigmas",
        "@(x) isnumeric(x)",
        "Sigma of each channel: a scalar or (z, y, x); 0 for no filter.",
        convert=(_sigmas_in, _sigmas_out),
    ),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch with its border.",
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
    Param(
        "localBBoxes",
        "@isnumeric",
        "Box of each region inside its batch.",
        convert="bboxes",
        py="local_bboxes",
    ),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("uuid", "@ischar", "Id of the call."),
    Param("debug", "@islogical", "Not used."),
    task="B02",
)
def unmix_channels_gaussian_block(
    batch_inds: Sequence[int],
    zarr_fullpaths: Sequence[str],
    unmix_fullpath: str,
    unmix_factors: Any,
    unmix_sigmas: Any,
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    region_bboxes: Sequence[Sequence[int]],
    local_bboxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    uuid: str = "",
    debug: bool = False,
) -> list[bool] | None:
    """``unmix_channels_gaussian_block``: each batch with its border (``batch_bboxes``), each
    channel with a nonzero sigma filtered (``imgaussfilt3`` in single), unmixed, its region
    (``local_bboxes``) written into ``region_bboxes`` of the result zarr in its class; then the
    flag. Nothing is done if the flag exists (unless ``overwrite``)."""
    from petakit5d.io import get_image_dtype, write_zarr
    from petakit5d.kernels.crop import crop_3d
    from petakit5d.matlabcompat.filters import imgaussfilt3

    del uuid, debug
    paths = [os.fspath(p) for p in zarr_fullpaths]
    inds = _block_start(batch_inds, flag_fullname, overwrite, paths)
    if inds is None:
        return None
    nc = len(paths)
    factors = _numbers(unmix_factors)
    dtype = _check_outputs(factors, nc, unmix_fullpath)
    rows = _sigma_rows(unmix_sigmas)
    filtered = [_filtered(rows, c) for c in range(nc)]
    dts = [get_image_dtype(p, True) for p in paths]
    at_once = (
        not any(filtered)
        and all(d in _KERNEL for d in [*dts, dtype])
        and sum(_bytes(d) for d in dts) + _bytes(dtype) <= 8
    )

    boxes = np.asarray(batch_bboxes, dtype=np.int64).reshape(-1, 6)
    regions = np.asarray(region_bboxes, dtype=np.int64).reshape(-1, 6)
    locals_ = np.asarray(local_bboxes, dtype=np.int64).reshape(-1, 6)
    buf = _Buffers(nc if at_once else 1)
    zero_bgs = [0.0] * nc
    done = [False] * len(inds)
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        box = [int(v) for v in boxes[i]]
        shape = _box_shape(box)
        batch = buf.out(shape, dtype)
        if at_once:
            ims = [buf.read(c, paths[c], box) for c in range(nc)]
            _sum(ims, factors, zero_bgs, 0.0, dtype, None, batch, None)
            del ims
        else:
            acc: np.ndarray | None = None
            for c in range(nc):
                im_c = buf.read(0, paths[c], box)
                if filtered[c]:
                    im_c = imgaussfilt3(_single(im_c), _gauss_sigma(rows[c]))
                acc = _sum([im_c], [factors[c]], [0.0], 0.0, "", acc, acc, None)
                del im_c
            _sum([], [], [], 0.0, dtype, acc, batch, None)
            del acc
        out = crop_3d(batch, [int(v) for v in locals_[i]])
        write_zarr(out, unmix_fullpath, bbox=[int(v) for v in regions[i]])
        del out
        done[i] = True
        _toc(t0)
    if all(done):
        open(flag_fullname, "w").close()
    return done


# ---------------------------------------------------------------------------------------------
# large zarr files
# ---------------------------------------------------------------------------------------------


@register(
    "XR_unmix_channels_zarr",
    Param("zarrFullpaths", "@(x) ischar(x) || iscell(x)", "Zarr files of the channels."),
    Param("unmixFactors", "@(x) isnumeric(x)", "Unmixing factor of each channel."),
    Param("mode", "@ischar", "'linear', or 'gaussian' (each channel filtered with its sigma)."),
    Param(
        "unmixSigmas",
        "@isnumeric",
        "Sigmas of the gaussian mode, one per channel: a scalar or (z, y, x).",
        convert=(_sigmas_in, _sigmas_out),
    ),
    Param("resultDirName", "@ischar", "Result folder under the data path."),
    Param("channelInd", "@isnumeric", "The channel to unmix (0-based).", convert="index1"),
    Param(
        "save16bit", "@islogical", "Not used (the result has the channel's class).", convert="same"
    ),
    Param("batchSize", "@isvector", "Batch size (z, y, x) of a block.", convert="size_yxz"),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of the result.", convert="size_yxz"),
    Param("borderSize", "@isvector", "Border (z, y, x) of each batch.", convert="size_yxz"),
    Param("parseCluster", "@islogical", "Run the blocks as Slurm jobs."),
    Param("parseParfor", "@islogical", "Run the blocks on a local process pool."),
    Param("masterCompute", "@islogical", "This process computes blocks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a block's job."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    task="B02",
)
def unmix_channels_zarr(
    zarr_fullpaths: str | Sequence[str],
    unmix_factors: Any,
    *,
    mode: str = "linear",
    unmix_sigmas: Any = None,
    result_dir_name: str = "Unmixed",
    channel_ind: int = 0,
    save_16bit: bool = True,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    border_size: Sequence[int] = (0, 0, 0),
    parse_cluster: bool = True,
    parse_parfor: bool = False,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 3,
    config_file: str = "",
    mcc_mode: bool = False,
    uuid: str = "",
    debug: bool = False,
) -> list[bool] | None:
    """Unmix large zarr files in blocks (``XR_unmix_channels_zarr``): the batches of the image
    (``batch_size``, whole chunks of ``block_size``, with ``border_size`` in 'gaussian' mode) as
    ``unmix_channels_block`` or ``unmix_channels_gaussian_block`` tasks through the job
    framework into ``<data path>/<result_dir_name>/<name>_<uuid>.zarr``, then moved into place
    as ``<name>.zarr``, in the class of the channel ``channel_ind``. Nothing is done if the
    result exists; a result another writer published is kept. Returns the tasks' done flags
    (None if nothing was done)."""
    from petakit5d.core.blocks import zarr_chunk_coordinates
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive, strip_trailing_sep
    from petakit5d.io import create_zarr, get_image_dtype, get_image_size
    from petakit5d.jobs import framework
    from petakit5d.jobs.files import exist
    from petakit5d.jobs.publish import keep_other_result, publish_result
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task

    del save_16bit
    if not uuid:
        uuid = get_uuid()
    paths = [zarr_fullpaths] if isinstance(zarr_fullpaths, str | os.PathLike) else zarr_fullpaths
    paths = [strip_trailing_sep(os.fspath(p)) for p in paths]
    ci = int(channel_ind)
    data_path, fsname, _ext_name = fileparts(paths[ci])
    unmix_path = f"{data_path}/{result_dir_name}/"
    os.makedirs(unmix_path, exist_ok=True)
    full = f"{unmix_path}/{fsname}.zarr"
    tmp = f"{unmix_path}/{fsname}_{uuid}.zarr"
    if exist(full):
        print("The output results exist, skip it!", flush=True)
        return None
    print(f"Start Large-file channel unmixing for {fsname}...", flush=True)
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(data_path, job_log_dir).parse_cluster)

    flag_path = f"{unmix_path}/zarr_flag/{fsname}_{uuid}/"
    if not os.path.isdir(flag_path):
        mkdir_recursive(flag_path)
    im_size = [int(v) for v in get_image_size(paths[ci])[:3]]
    dtype = get_image_dtype(paths[ci], True)
    batch = [min(n, int(b)) for n, b in zip(im_size, batch_size, strict=True)]
    block = [min(b, int(k)) for b, k in zip(batch, block_size, strict=True)]
    batch_boxes, region_boxes, local_boxes = zarr_chunk_coordinates(
        im_size,
        batch_size=batch,
        block_size=block,
        same_batch_size=False,
        border_size=[int(v) for v in border_size],
    )
    if not os.path.isdir(tmp):
        n_chunks = 1
        for n, b in zip(im_size, block, strict=True):
            n_chunks *= math.ceil(n / b)
        create_zarr(
            tmp,
            data_size=im_size,
            block_size=block,
            dtype=dtype,
            dim_separator="/" if n_chunks > 10000 else ".",
        )

    num_batch = batch_boxes.shape[0]
    task_size = max(1, min(math.ceil(num_batch / 10000), 20))
    num_tasks = math.ceil(num_batch / task_size)
    # (MATLAB's task strings: mat2str of the factors and sigmas)
    factors = _through_mat2str(tuple(_numbers(unmix_factors).tolist()), 15)
    sigmas = _sigmas_text(unmix_sigmas, 15)
    tasks, outputs = [], []
    for t in range(num_tasks):
        inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
        flag = f"{flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
        outputs.append(flag)
        if mode == "linear":
            tasks.append(
                Task.of(
                    unmix_channels_block,
                    inds,
                    paths,
                    tmp,
                    list(factors),
                    flag,
                    batch_boxes[inds].tolist(),
                    uuid=uuid,
                    debug=bool(debug),
                )
            )
            mem_factor = len(paths) + 1
        elif mode == "gaussian":
            tasks.append(
                Task.of(
                    unmix_channels_gaussian_block,
                    inds,
                    paths,
                    tmp,
                    list(factors),
                    None if _empty(sigmas) else [list(np.atleast_1d(r)) for r in sigmas],
                    flag,
                    batch_boxes[inds].tolist(),
                    region_boxes[inds].tolist(),
                    local_boxes[inds].tolist(),
                    uuid=uuid,
                    debug=bool(debug),
                )
            )
            mem_factor = (len(paths) + 1) * 1.5
        else:
            raise UnmixError(
                "Unrecognized function or variable 'memFactor'.", "MATLAB:UndefinedFunction"
            )

    mem_allocate = float(np.prod(batch, dtype=np.float64)) * 4 / 1024**3 * mem_factor
    inputs = [paths[ci]] * num_tasks
    done: list[bool] = []
    if parse_cluster or not parse_parfor:
        for _ in range(2):
            done = framework.run_tasks(
                inputs,
                outputs,
                tasks,
                final_out_fullpath=full,
                parse_cluster=parse_cluster,
                master_compute=master_compute,
                task_batch_num=1,
                cpus_per_task=cpus_per_task,
                mem_allocate=mem_allocate,
                job_log_dir=job_log_dir,
                mcc_mode=mcc_mode,
                config_file=config_file,
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
        raise UnmixError("Some blocks are not finished!")

    # another copy of this task may have published the result: keep it (this copy's
    # temporary result can be incomplete then)
    if not keep_other_result(full, tmp, flag_path):
        publish_result(tmp, full, flag_path)
    return done


# ---------------------------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------------------------


@register(
    "XR_unmix_channels_data_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param("unmixFactors", "@(x) isnumeric(x)", "Unmixing factor of each channel."),
    Param("mode", "@ischar", "'linear', or 'gaussian' (each channel filtered with its sigma)."),
    Param(
        "unmixSigmas",
        "@isnumeric",
        "Sigmas of the gaussian mode, one per channel: a scalar or (z, y, x).",
        convert=(_sigmas_in, _sigmas_out),
    ),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files."),
    Param("channelInd", "@isnumeric", "The channel to unmix (0-based).", convert="index1"),
    Param("FFCorrection", "@islogical", "Flat-field correct the channels first (not large files)."),
    Param("lowerLimit", "@isnumeric", "Lowest value of the normalized flat fields."),
    Param("FFImagePaths", "@iscell", "Flat-field image of each channel."),
    Param("backgroundPaths", "@iscell", "Background image of each channel."),
    Param(
        "constBackground", "@isnumeric", "Background (scalar, or one per channel) subtracted first."
    ),
    Param("constOffset", "@(x) isnumeric(x)", "If set, added after the unmixing.", convert="same"),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("largeFile", "@islogical", "Large zarr files: unmix in blocks (unmix_channels_zarr)."),
    Param("saveZarr", "@islogical", "Save the results as zarr."),
    Param("save16bit", "@islogical", "Save uint16; otherwise single.", convert="same"),
    Param("batchSize", "@isvector", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of zarr results.", convert="size_yxz"),
    Param(
        "borderSize",
        "@isvector",
        "Border (z, y, x) of the batches of large files.",
        convert="size_yxz",
    ),
    Param("parseCluster", "@islogical", "Run the frames as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    task="B02",
)
def unmix_channels_dataset(
    data_paths: str | Sequence[str],
    *,
    unmix_factors: Any = None,
    mode: str = "linear",
    unmix_sigmas: Any = None,
    result_dir_name: str = "Unmixed",
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    channel_ind: int = 0,
    ff_correction: bool = False,
    lower_limit: Any = 0.4,
    ff_image_paths: Sequence[str] = ("", ""),
    background_paths: Sequence[str] = ("", ""),
    const_background: Any = None,
    const_offset: Any = None,
    zarr_file: bool = False,
    large_file: bool = False,
    save_zarr: bool = False,
    save_16bit: bool = True,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    border_size: Sequence[int] = (0, 0, 0),
    parse_cluster: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 3,
    config_file: str = "",
    mcc_mode: bool = False,
    uuid: str = "",
    debug: bool = False,
) -> None:
    """Unmix a channel of every time point of the datasets (``XR_unmix_channels_data_wrapper``):
    the files are grouped by channel (``channel_patterns``, in each folder, or over the folders
    when each holds one channel), groups without every channel are left out (a warning), and each
    group is one ``unmix_channels_file`` task (``unmix_channels_zarr`` with ``large_file``)
    through the job framework, its result in ``<data path>/<result_dir_name>/``."""
    import warnings

    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import simplify_path
    from petakit5d.io import get_image_size
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import (
        group_image_filenames_by_channels,
        parse_image_filenames,
    )

    del debug
    paths = [data_paths] if isinstance(data_paths, str | os.PathLike) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
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
    nc = len(channel_patterns)
    factors = _numbers(unmix_factors)
    if len(set(ch_inds)) != nc:
        raise UnmixError(
            "Number of channels of the images does not match the number of channel patterns!"
        )
    if nc != factors.size:
        raise UnmixError("The number of images does not match that of the unmixing factors!")
    if ff_correction and (nc != len(ff_image_paths) or nc != len(background_paths)):
        raise UnmixError(
            "The number of flat field or background images does not match number of channels!"
        )
    bg = const_background
    if not _empty(bg):
        b = _numbers(bg)
        if b.size == 1:
            bg = tuple([float(b[0])] * nc)
        elif b.size != nc:
            raise UnmixError(
                "The number of constBackground values does not match number of channels!"
            )

    # group the files by channel: over the folders if each holds one channel, else per folder
    single_per_folder = (
        nd > 1
        and nc == nd
        and all(
            len({ch_inds[i] for i in range(len(fnames)) if fd_inds[i] == d}) == 1 for d in range(nd)
        )
    )
    if single_per_folder:
        groups = group_image_filenames_by_channels(fnames, channel_patterns, ch_inds)
    else:
        parts = []
        for d in range(nd):
            idx = [i for i in range(len(fnames)) if fd_inds[i] == d]
            g = group_image_filenames_by_channels(
                [fnames[i] for i in idx], channel_patterns, [ch_inds[i] for i in idx]
            )
            if g.size:
                # (-1: no file of the channel in the group)
                parts.append(np.where(g >= 0, np.asarray(idx, np.int64)[np.maximum(g, 0)], -1))
        groups = np.concatenate(parts, axis=0) if parts else np.zeros((0, nc), np.int64)
    incomplete = np.any(groups < 0, axis=1)
    if incomplete.any():
        warnings.warn(
            "Some images do not have all channels and are excluded for unmixing!", stacklevel=2
        )
        groups = groups[~incomplete]

    ci = int(channel_ind)
    ext = ".zarr" if save_zarr else ".tif"
    im_size = get_image_size(f"{paths[fd_inds[0]]}/{fnames[0]}")

    # (MATLAB's task strings: the factors with %0.10f, the sigmas with mat2str(..., 10), the
    # backgrounds and the offset with mat2str, lowerLimit with %.20d)
    factors_t = _through_text(tuple(factors.tolist()), "%0.10f")
    sigmas_t = _sigmas_text(unmix_sigmas, 10)
    bg_t = _through_mat2str(bg, 15)
    off_t = _through_mat2str(const_offset, 15)
    lower_t = _through_text(lower_limit, "%.20d")
    frame_fullpaths, result_fullpaths, tasks = [], [], []
    for row in groups:
        files = [f"{paths[fd_inds[x]]}/{fnames[x]}" for x in row]
        frame_fullpaths.append(files[ci])
        p = int(row[ci])
        result_fullpaths.append(f"{result_paths[fd_inds[p]]}/{fsns[p]}{ext}")
        sig = None if _empty(sigmas_t) else [list(np.atleast_1d(r)) for r in sigmas_t]
        if large_file:
            tasks.append(
                Task.of(
                    unmix_channels_zarr,
                    files,
                    list(factors_t),
                    mode=mode,
                    unmix_sigmas=sig,
                    result_dir_name=result_dir_name,
                    channel_ind=ci,
                    save_16bit=bool(save_16bit),
                    batch_size=[int(v) for v in batch_size],
                    block_size=[int(v) for v in block_size],
                    border_size=[int(v) for v in border_size],
                    parse_cluster=bool(parse_cluster),
                    master_compute=bool(master_compute),
                    job_log_dir=job_log_dir,
                    cpus_per_task=int(cpus_per_task),
                    config_file=config_file,
                    mcc_mode=bool(mcc_mode),
                    uuid=uuid,
                )
            )
        else:
            tasks.append(
                Task.of(
                    unmix_channels_file,
                    files,
                    list(factors_t),
                    mode=mode,
                    unmix_sigmas=sig,
                    result_dir_name=result_dir_name,
                    channel_ind=ci,
                    ff_correction=bool(ff_correction),
                    lower_limit=lower_t,
                    ff_image_paths=list(ff_image_paths),
                    background_paths=list(background_paths),
                    const_background=None if _empty(bg_t) else list(np.atleast_1d(bg_t)),
                    const_offset=None if _empty(off_t) else float(np.asarray(off_t).reshape(-1)[0]),
                    save_16bit=bool(save_16bit),
                    zarr_file=bool(zarr_file),
                    save_zarr=bool(save_zarr),
                    block_size=[int(v) for v in block_size],
                    uuid=uuid,
                )
            )

    if large_file:
        mem_allocate = float(np.prod([float(v) for v in batch_size])) / 1024**3 * 4 * 2.5
    else:
        mem_allocate = float(np.prod([float(v) for v in im_size])) / 1024**3 * 4 * 2.5
    framework.run_tasks(
        frame_fullpaths,
        result_fullpaths,
        tasks,
        parse_cluster=parse_cluster,
        job_log_dir=job_log_dir,
        master_compute=master_compute,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )
