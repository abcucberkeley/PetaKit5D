"""Flat-field correction of light-sheet frames: ``XR_LSFlatFieldCorrection`` and
``processFFCorrectionFrame``.

- ``prepare_flat_field(ls_image, background, plane_shape, ...)``: the flat field and the
  background as ``XR_LSFlatFieldCorrection`` makes them before the correction: a stack of
  flat-field images averaged over z (MATLAB's ``mean(LSImage, 3)``), both images center-cropped
  to the data's plane, the background subtracted, the flat field divided by its maximum and
  raised to ``lower_limit``;
- ``flatfield_correct_volume(rawdata, ls_image, background, ...)``: the whole function, the
  correction by ``petakit5d.kernels.flatfield`` (the mex for uint16 and single data, MATLAB's
  code for the other classes);
- ``flatfield_correct_file(frame, ls_image, background_image, const_offset, ...)``:
  ``processFFCorrectionFrame``, the images read from their tiff files.

Each step keeps MATLAB's classes and arithmetic: the z mean as MATLAB adds (in double for
integer stacks, in single for single ones), ``single(LS) - single(BG)``,
``single(LS) / single(max(LS(:)))`` (NaN ignored), ``max(LS, LowerLimit)`` in MATLAB's class
of both (NaN gives ``LowerLimit``), then ``single`` of both images. The crop takes rows and
columns ``D + 1 : D + n`` (1-based) with ``D = ceil((m - n) / 2)``, as MATLAB indexes; images
smaller than the data raise MATLAB's index errors, images with a third dimension MATLAB's size
error, with MATLAB's identifiers. uint16 results above 65535 (or below 0, with a negative
``const_offset``) wrap as the mex converts them; MATLAB's code for the other classes saturates
(known issue, kept).

Arrays are the port's: data ``(z, y, x)`` (or an image ``(y, x)``), images ``(y, x)``, stacks
of flat-field images ``(z, y, x)``. ``out``: a C-contiguous array of the result's class and the
data's shape to write into; frame after frame, a reused buffer saves the first touch of new
pages (a third to a half of the kernel's time). ``num_threads`` (None: the CPUs available)
never changes a result.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "flatfield_correct_file",
    "flatfield_correct_volume",
    "prepare_flat_field",
]

_F32 = np.dtype(np.float32)


class MatlabError(ValueError):
    """An error MATLAB raises, with its identifier."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


SIZE_ERROR = "MATLAB:sizeDimensionsMustMatch"
INDEX_ERROR = "MATLAB:badsubscript"


def _if(x: Any) -> bool:
    """MATLAB's ``if x``: true when x is not empty and all its elements are nonzero"""
    a = np.asarray(x)
    return bool(a.size) and bool(np.all(a != 0))


def _first(x: Any) -> bool:
    """the mex's ``mxGetScalar(x) != 0``: the first element (MATLAB's order) of a flag; an empty
    one is false"""
    a = np.asarray(x)
    if a.size == 0:
        return False
    return bool(a.reshape(-1, order="F")[0] != 0)  # NaN is true


def _double(x: Any) -> Any:
    """a Python number as MATLAB's double (the MATLAB-name layer gives integer-valued doubles as
    int; ``np.asarray`` would make them int64)"""
    if isinstance(x, int | float) and not isinstance(x, bool):
        return np.float64(x)
    return x


def _given(const_offset: Any) -> bool:
    """``~isempty(constOffset)``"""
    return const_offset is not None and np.asarray(const_offset).size > 0


def _plane_size(raw: np.ndarray) -> tuple[int, int]:
    """``size(squeeze(Rawdata(:, :, 1)))``: MATLAB's rows and columns of the data's planes"""
    from petakit5d.core.arrays import size_to_matlab

    ms = size_to_matlab(np.shape(raw))
    return int(ms[0]), int(ms[1])


def _crop(img: np.ndarray, im_size: tuple[int, int]) -> np.ndarray:
    """``img(D(1)+1:D(1)+ImSize(1), D(2)+1:D(2)+ImSize(2))`` with ``D = ceil((size(img) -
    ImSize) / 2)`` (a view), with MATLAB's errors"""
    from petakit5d.core.arrays import size_to_matlab

    ms = size_to_matlab(img.shape)
    if len(ms) != 2:
        raise MatlabError(SIZE_ERROR, "Arrays have incompatible sizes for this operation.")
    img2 = img.reshape(ms)
    sl = []
    for pos in range(2):
        m, n = int(ms[pos]), int(im_size[pos])
        d = (m - n + 1) // 2  # ceil((m - n) / 2) of integers
        if n > 0 and d + 1 < 1:
            raise MatlabError(
                INDEX_ERROR,
                f"Index in position {pos + 1} is invalid. Array indices must be positive "
                "integers or logical values.",
            )
        if n > 0 and d + n > m:
            raise MatlabError(
                INDEX_ERROR,
                f"Index in position {pos + 1} exceeds array bounds. Index must not exceed {m}.",
            )
        sl.append(slice(d, d + n) if n > 0 else slice(0, 0))
    return img2[sl[0], sl[1]]


def prepare_flat_field(
    ls_image: np.ndarray,
    background: np.ndarray,
    plane_shape: tuple[int, int],
    *,
    lower_limit: Any = 0.4,
    remove_ff_im_background: Any = True,
    ls_rescale: Any = True,
) -> tuple[np.ndarray, np.ndarray]:
    """The flat field and the background of ``XR_LSFlatFieldCorrection`` for data whose planes
    are ``plane_shape`` ``(ny, nx)``, as MATLAB has them before they reach the correction (each
    in its MATLAB class; the correction makes both single)::

        if ndims(LSImage) == 3, LSImage = squeeze(mean(LSImage, 3)); end
        LSImage = <center crop to the plane>; background = <center crop to the plane>;
        if removeFFImBackground, LSImage = single(LSImage) - single(background); end
        if LSRescale, LSImage = single(LSImage) / single(max(LSImage(:))); end
        LSImage = max(LSImage, LowerLimit);

    ``ls_image``: a ``(y, x)`` image or a ``(z, y, x)`` stack; ``background``: a ``(y, x)``
    image. Raises ``MatlabError`` with MATLAB's identifiers where MATLAB raises.
    """
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import size_to_matlab
    from petakit5d.matlabcompat import reduce

    ls = np.asarray(ls_image)
    bg = np.asarray(background)
    if len(size_to_matlab(ls.shape)) == 3:
        ls = np.asarray(reduce.mean(ls, axis=-3))
    im_size = (int(plane_shape[0]), int(plane_shape[1]))
    ls = _crop(ls, im_size)
    bg = _crop(bg, im_size)
    if _if(remove_ff_im_background):
        ls = _single_op(np.subtract, dtypes.minus, _single(ls), _single(bg))
    if _if(ls_rescale):
        m, _ = reduce.max(ls)
        ls = _single_op(np.divide, dtypes.rdivide, _single(ls), _single(np.asarray(m)))
    return _lower_limit(ls, lower_limit), bg


def _single_op(op, matlab_op, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``matlab_op(a, b)`` (``core.dtypes``) of single arrays: NumPy's IEEE operation when no
    result is NaN (the same bits), ``matlab_op`` for MATLAB's NaN payloads otherwise"""
    with np.errstate(all="ignore"):
        r = op(a, b)
    if np.isnan(r).any():
        return np.asarray(matlab_op(a, b))
    return r


def _lower_limit(ls: np.ndarray, lower_limit: Any) -> np.ndarray:
    """``max(LSImage, LowerLimit)`` (``reduce.maximum``); for single data and a nonzero real
    scalar limit, NumPy's fmax with the limit rounded to single (the same: a NaN gives the limit,
    and only a zero limit could tie with a value of the other sign)"""
    from petakit5d.matlabcompat import reduce

    lower_limit = _double(lower_limit)
    lim = np.asarray(lower_limit)
    if ls.dtype == _F32 and lim.ndim == 0 and lim.dtype in (np.float64, _F32):
        lo = np.float32(lim)
        if lo != 0 and not np.isnan(lo):
            return np.fmax(ls, lo)
    return np.asarray(reduce.maximum(ls, lower_limit))


def _single(x: Any) -> np.ndarray:
    from petakit5d.core.dtypes import to_class

    a = np.asarray(x)
    return a if a.dtype == _F32 else np.asarray(to_class(a, "single"))


def _correct(
    raw: np.ndarray,
    flat: np.ndarray,
    bg: np.ndarray,
    const_offset: Any,
    cast_data_type: Any,
    out: np.ndarray | None,
    num_threads: int | None,
) -> np.ndarray:
    """XR_LSFlatFieldCorrection's ``try`` mex ``catch`` MATLAB code, into ``out`` if given"""
    from petakit5d.kernels import flatfield as ff

    f, b = _single(flat), _single(bg)
    off = _double(const_offset) if _given(const_offset) else None
    try:
        return ff.flat_field_correction_3d(
            raw,
            f,
            b,
            cast_data_type=_first(cast_data_type),
            const_offset=off,
            out=out,
            num_threads=num_threads,
        )
    except ff.MatlabError:
        pass
    # the classes and sizes the mex refuses: MATLAB's code (its `if castDataType`)
    r = ff.apply_flat_field(raw, f, b, const_offset=off, cast_data_type=_if(cast_data_type))
    if out is None:
        return r
    if out.dtype != r.dtype or out.shape != r.shape:
        raise ValueError(
            f"out must be a {r.dtype} array of shape {r.shape}, got {out.dtype} {out.shape}"
        )
    out[...] = r
    return out


def _to_matlab_result(result: Any) -> Any:
    from petakit5d.core.arrays import to_matlab

    return to_matlab(result)


@register(
    "XR_LSFlatFieldCorrection",
    Param("Rawdata", None, "Data, (z, y, x) or (y, x).", convert="volume"),
    Param(
        "LSImage", None, "Flat-field image (y, x), or a stack of them (z, y, x).", convert="volume"
    ),
    Param("background", None, "Camera background image (y, x).", convert="volume"),
    Param("LowerLimit", "@isnumeric", "Lowest value of the normalized flat field."),
    Param(
        "constOffset",
        "@(x) isempty(x) || isnumeric(x)",
        "If set, added after the correction instead of the background.",
        convert="same",
    ),
    Param(
        "removeFFImBackground",
        "@(x) islogical(x)",
        "Subtract the background from the flat-field image.",
    ),
    Param("LSRescale", "@(x) islogical(x)", "Divide the flat field by its maximum."),
    Param(
        "castDataType",
        "@(x) islogical(x)",
        "Cast the result to the data's class; otherwise single.",
    ),
    Param("out", None, "Array to write the result into (reused over frames).", python_only=True),
    Param(
        "numThreads",
        None,
        "Threads of the correction (None: the CPUs available).",
        py="num_threads",
        python_only=True,
    ),
    task="D01",
    matlab_out=_to_matlab_result,
)
def flatfield_correct_volume(
    rawdata: np.ndarray,
    ls_image: np.ndarray,
    background: np.ndarray,
    *,
    lower_limit: Any = 0.4,
    const_offset: Any = None,
    remove_ff_im_background: Any = True,
    ls_rescale: Any = True,
    cast_data_type: Any = True,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """Flat-field correction of light-sheet data (``XR_LSFlatFieldCorrection``).

    The flat field ``ls_image`` (an image, or a stack averaged over z) and the camera
    ``background`` are prepared as ``prepare_flat_field`` describes, then every plane of
    ``rawdata`` is corrected::

        out = max(single(raw) - background, 0) ./ flat + (const_offset or background)

    uint16 and single data in the mex's arithmetic (uint16 results wrap above 65535), other
    classes in MATLAB's code (saturated). The result has the data's class with
    ``cast_data_type``, else single. ``out``: an array of the result's class and the data's
    shape, written and returned (reuse it frame after frame).
    """
    raw = np.asarray(rawdata)
    flat, bg = prepare_flat_field(
        ls_image,
        background,
        _plane_size(raw),
        lower_limit=lower_limit,
        remove_ff_im_background=remove_ff_im_background,
        ls_rescale=ls_rescale,
    )
    return _correct(raw, flat, bg, const_offset, cast_data_type, out, num_threads)


@register(
    "processFFCorrectionFrame",
    Param("frame", None, "Data, (z, y, x) or (y, x).", convert="volume"),
    Param("LSImage", None, "Tiff file of the flat-field image (or stack).", convert="same"),
    Param("BackgroundImage", None, "Tiff file of the camera background image."),
    Param(
        "constOffset",
        None,
        "If not empty, added after the correction instead of the background.",
        convert="same",
    ),
    Param("LowerLimit", None, "Lowest value of the normalized flat field."),
    Param("removeFFImBackground", None, "Subtract the background from the flat-field image."),
    Param("LSRescale", None, "Divide the flat field by its maximum."),
    Param("castDataType", None, "Cast the result to the data's class; otherwise single."),
    Param("out", None, "Array to write the result into (reused over frames).", python_only=True),
    Param(
        "numThreads",
        None,
        "Threads of the correction (None: the CPUs available).",
        py="num_threads",
        python_only=True,
    ),
    task="D01",
    matlab_out=_to_matlab_result,
)
def flatfield_correct_file(
    frame: np.ndarray,
    ls_image: str,
    background_image: str,
    const_offset: Any,
    lower_limit: Any = 0.4,
    remove_ff_im_background: Any = True,
    ls_rescale: Any = True,
    cast_data_type: Any = True,
    *,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """Flat-field correction of a frame in memory with the flat field and the background read
    from tiff files (``processFFCorrectionFrame``); the rest as ``flatfield_correct_volume``.
    The files are read at every call, as MATLAB reads them."""
    from petakit5d.io import read_tiff

    ls = read_tiff(ls_image)
    bk = read_tiff(background_image)
    return flatfield_correct_volume(
        frame,
        ls,
        bk,
        lower_limit=lower_limit,
        const_offset=const_offset,
        remove_ff_im_background=remove_ff_im_background,
        ls_rescale=ls_rescale,
        cast_data_type=cast_data_type,
        out=out,
        num_threads=num_threads,
    )
