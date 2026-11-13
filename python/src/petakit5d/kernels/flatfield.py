"""Flat-field correction of a ``(z, y, x)`` volume, the kernel of ``flat_field_correction_3d_mex``.

The Python side of ``petakit5d._ext.flatfield`` (``python/cpp/flatfield``):

- ``flat_field_correction_3d(vol, flat, background, cast_data_type=, const_offset=)``: per voxel
  ``max((vol - background) / flat, 0) + (const_offset or background)`` in single precision, with
  the ``(y, x)`` single images ``flat`` and ``background`` applied to every plane;
- ``flat_field_correction_3d_mex(im, flat, background, cast_data_type[, const_offset])``: the
  mex's call with MATLAB's argument values (the cast flag and the offset read as
  ``mxGetScalar`` reads them) and its checks, in its order, with its error identifiers;
- ``apply_flat_field(raw, flat, background, ...)``: the correction step of
  ``XR_LSFlatFieldCorrection`` once the flat field is prepared (``processFFCorrectionFrame``'s
  path): the images made single, then the mex, or where the mex refuses (other classes), the
  MATLAB code after its ``catch``;
- ``flat_field_correction_3d_ref``: a NumPy reference of the kernel (tests; used when the
  extension is not built).

The mex takes uint16 and single data. uint16 data gives uint16 (``cast_data_type`` true) or
single; single data gives single. Rounding of uint16 output as the mex's: ``(uint16)(v + 0.5)``
computed in double and converted through 32 bits, so values above 65535 wrap (65536 is 0) and
negative ones (a negative ``const_offset``) wrap from 65535 down, and NaN and values beyond the
int32 range give 0. MATLAB's code for the other classes casts with saturation instead
(``cast``). Both are kept as MATLAB has them (known issue). Arrays are C-contiguous ``(z, y, x)``
volumes or ``(y, x)`` images (other layouts are copied); results have the input's shape.
``num_threads`` (None or 0: the CPUs available) never changes a result.
"""

from __future__ import annotations

from types import ModuleType
from typing import Any

import numpy as np

__all__ = [
    "MatlabError",
    "apply_flat_field",
    "extension_available",
    "flat_field_correction_3d",
    "flat_field_correction_3d_mex",
    "flat_field_correction_3d_ref",
]

INPUT_ERROR = "flat_field_correction_3d_mex:inputError"
TYPE_ERROR = "flat_field_correction_3d_mex:dataTypeError"

_U16 = np.dtype(np.uint16)
_F32 = np.dtype(np.float32)
# below this many voxels one thread (a team costs more than the loop)
_SERIAL_VOXELS = 1 << 16

_ext_module: ModuleType | None = None
_ext_missing = False
_available: int | None = None


class MatlabError(ValueError):
    """An error MATLAB raises, with its identifier
    (``'flat_field_correction_3d_mex:inputError'``, ...)."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _ext() -> ModuleType | None:
    global _ext_module, _ext_missing
    if _ext_module is None and not _ext_missing:
        from petakit5d import _ext

        try:
            _ext_module = _ext.load("flatfield")
        except ImportError:
            _ext_missing = True
    return _ext_module


def extension_available() -> bool:
    """Whether the C++ kernel is built (otherwise the NumPy reference runs)."""
    return _ext() is not None


def _threads(num_threads: int | None) -> int:
    """``core.threads.resolve_threads`` with the CPUs counted once per process"""
    global _available
    if _available is None:
        from petakit5d.core.threads import available_cpus

        _available = available_cpus()
    if num_threads is None or num_threads == 0:
        return _available
    n = int(num_threads)
    if n < 0:
        raise ValueError(f"thread count must be positive, got {n}")
    return min(n, _available)


# ---------------------------------------------------------------------------------------------
# the mex's checks
# ---------------------------------------------------------------------------------------------


def _msize(a: np.ndarray) -> tuple[int, ...]:
    from petakit5d.core.arrays import size_to_matlab

    return size_to_matlab(a.shape)


def _check(a: np.ndarray, flat: np.ndarray, bg: np.ndarray) -> None:
    """the mex's checks of the arrays, in its order"""
    ms = _msize(a)
    if len(ms) > 3:
        raise MatlabError(INPUT_ERROR, "The data must have at most 3 dimensions")
    # MATLAB's arrays have at least 2 dimensions, so only the first two sizes are compared
    for img, what in ((flat, "flat field"), (bg, "background")):
        mi = _msize(img)
        if len(mi) > 2:
            raise MatlabError(INPUT_ERROR, f"The {what} image must be 1d or 2d")
        if mi[:2] != ms[:2]:
            raise MatlabError(
                INPUT_ERROR, f"The size of {what} image does not matach that of the data"
            )


def _check_classes(a: np.ndarray, flat: np.ndarray, bg: np.ndarray) -> None:
    if flat.dtype != _F32 or bg.dtype != _F32:
        raise MatlabError(
            INPUT_ERROR,
            "The data type of the flat field and background images must be single!",
        )
    if a.dtype not in (_U16, _F32):
        raise MatlabError(TYPE_ERROR, "Data type not suppported")


def _mx_scalar(x: Any) -> float:
    """``mxGetScalar``: the first element of a numeric, logical or char array as a double (the
    real part of a complex one; 0 for an empty array)"""
    a = np.asarray(x)
    if a.size == 0:
        return 0.0
    v = a.reshape(-1, order="F")[0]
    if a.dtype.kind == "U":
        return float(ord(str(v)[0])) if str(v) else 0.0
    if a.dtype.kind == "c":
        v = v.real
    return float(v)


# ---------------------------------------------------------------------------------------------
# the kernel
# ---------------------------------------------------------------------------------------------


def _vol3(a: np.ndarray) -> np.ndarray:
    """the (z, y, x) volume the kernel reads: an image is (1, ny, nx), a row (1, 1, n); leading
    singleton axes beyond 3 dimensions are dropped (checked: at most 3 MATLAB dimensions)"""
    if a.ndim < 3:
        a = a.reshape((1,) * (3 - a.ndim) + a.shape)
    elif a.ndim > 3:
        a = a.reshape(a.shape[-3:])
    return np.ascontiguousarray(a)


def _img2(img: np.ndarray, ny: int, nx: int) -> np.ndarray:
    return np.ascontiguousarray(img.reshape(ny, nx))


def _round_output_ref(w: np.ndarray) -> np.ndarray:
    """the mex's ``(uint16)(w + 0.5)`` of float32 values ``w``: truncation through int32,
    low 16 bits; NaN and values outside int32's range give 0"""
    d = w.astype(np.float64) + 0.5
    ok = (d > -2147483649.0) & (d < 2147483648.0)
    t = np.where(ok, d, -2147483648.0).astype(np.int64)
    return (t & 0xFFFF).astype(np.uint16)


def _ref(v: np.ndarray, f: np.ndarray, b: np.ndarray, out_dtype, offset) -> np.ndarray:
    # each step a float32 operation; a NaN v is replaced by 0 before it reaches a result
    from petakit5d.core import dtypes

    with np.errstate(all="ignore"):
        w = (v.astype(np.float32) - b) / f
        w = np.where(w >= 0, w, np.float32(0))
        add = b if offset is None else offset
        w = dtypes.x86_nan(w + add, w, add)  # Inf - Inf: MATLAB's NaN on every processor
    if out_dtype == _U16:
        return _round_output_ref(w)
    return w.astype(np.float32, copy=False)


def _offset(const_offset: Any) -> np.float32 | None:
    if const_offset is None:
        return None
    return np.float32(_mx_scalar(const_offset))


def _run(
    a: np.ndarray,
    flat: np.ndarray,
    bg: np.ndarray,
    cast: bool,
    offset: np.float32 | None,
    num_threads: int | None,
    reference: bool = False,
    out: np.ndarray | None = None,
) -> np.ndarray:
    """the kernel on checked arrays: the result in a's shape (written into out if given)"""
    v = _vol3(a)
    nz, ny, nx = v.shape
    f = _img2(flat, ny, nx)
    b = _img2(bg, ny, nx)
    out_dtype = _U16 if (a.dtype == _U16 and cast) else _F32
    if out is not None:
        if out.dtype != out_dtype or out.shape != a.shape or not out.flags.c_contiguous:
            raise ValueError(
                f"out must be a C-contiguous {out_dtype} array of shape {a.shape}, got "
                f"{out.dtype} {out.shape}"
            )
        if not out.flags.writeable:
            raise ValueError("out must be writeable")
        if any(np.may_share_memory(out, x) for x in (a, flat, bg)):
            raise ValueError("out must not overlap the inputs")
    ext = None if reference else _ext()
    if ext is None:
        r = _ref(v, f, b, out_dtype, offset).reshape(a.shape)
        if out is None:
            return r
        out[...] = r
        return out
    res = np.empty(v.shape, dtype=out_dtype) if out is None else out.reshape(v.shape)
    if res.size:
        nt = 1 if v.size < _SERIAL_VOXELS else _threads(num_threads)
        has = offset is not None
        ext.flat_field_correct(v, f, b, res, has, float(offset) if has else 0.0, nt)
    return res.reshape(a.shape) if out is None else out


def flat_field_correction_3d(
    vol: np.ndarray,
    flat: np.ndarray,
    background: np.ndarray,
    *,
    cast_data_type: bool = True,
    const_offset: float | None = None,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """Flat-field correction of each plane of ``vol`` (``(z, y, x)`` or ``(y, x)``, uint16 or
    single) by the single ``(y, x)`` images ``flat`` and ``background``, as
    ``flat_field_correction_3d_mex``:

        out = max((single(vol) - background) / flat, 0) + (const_offset or background)

    in single precision (NaN gives 0 before the offset is added). uint16 ``vol`` gives uint16
    with the mex's rounding (``cast_data_type``; values beyond 65535 and negative ones wrap) or
    single; single ``vol`` gives single. ``const_offset`` is rounded to single. Raises
    ``MatlabError`` with the mex's identifiers for bad sizes and classes. ``out``: a
    C-contiguous array of the result's class and ``vol``'s shape to write into (a buffer reused
    over frames saves the first touch of new pages, a third of the time at one thread).
    """
    a, f, b = np.asarray(vol), np.asarray(flat), np.asarray(background)
    _check(a, f, b)
    _check_classes(a, f, b)
    return _run(a, f, b, bool(cast_data_type), _offset(const_offset), num_threads, out=out)


def flat_field_correction_3d_ref(
    vol: np.ndarray,
    flat: np.ndarray,
    background: np.ndarray,
    *,
    cast_data_type: bool = True,
    const_offset: float | None = None,
) -> np.ndarray:
    """``flat_field_correction_3d`` in NumPy (the reference of the kernel's tests)."""
    a, f, b = np.asarray(vol), np.asarray(flat), np.asarray(background)
    _check(a, f, b)
    _check_classes(a, f, b)
    return _run(a, f, b, bool(cast_data_type), _offset(const_offset), None, reference=True)


def flat_field_correction_3d_mex(*args: Any, num_threads: int | None = None) -> np.ndarray:
    """``flat_field_correction_3d_mex(im, flat, background, castDataType[, constOffset])`` with
    ``(z, y, x)`` arrays: the mex's checks in its order (number of arguments, at most 3
    dimensions, 2-D images of the data's plane size, then the classes: single images, uint16 or
    single data), ``castDataType`` and ``constOffset`` read as ``mxGetScalar`` reads them (the
    first element as a double; a nonzero or NaN flag is true; the offset rounded to single).
    """
    if len(args) not in (4, 5):
        raise MatlabError(INPUT_ERROR, "Number of input arguments must be 4 or 5")
    a, f, b = (np.asarray(x) for x in args[:3])
    _check(a, f, b)
    cast = _mx_scalar(args[3]) != 0  # (bool) of a double: NaN is true
    offset = np.float32(_mx_scalar(args[4])) if len(args) == 5 else None
    _check_classes(a, f, b)
    return _run(a, f, b, cast, offset, num_threads)


# ---------------------------------------------------------------------------------------------
# XR_LSFlatFieldCorrection's correction step
# ---------------------------------------------------------------------------------------------


def _single(x: Any) -> np.ndarray:
    from petakit5d.core.dtypes import to_class

    a = np.asarray(x)
    return a if a.dtype == _F32 else np.asarray(to_class(a, "single"))


def _matlab_code(
    raw: np.ndarray, flat: np.ndarray, bg: np.ndarray, const_offset: Any, cast_data_type: bool
) -> np.ndarray:
    """XR_LSFlatFieldCorrection's MATLAB code after the mex's ``catch``:

    Temp = single(Rawdata) - single(background);
    Temp = max(Temp, 0);
    Temp = Temp ./ single(LSImage) + constOffset;   % or + single(background)
    if castDataType, Temp = cast(Temp, class(Rawdata)); end
    """
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import matlab_class
    from petakit5d.kernels.ops import clamp

    cls = matlab_class(raw)
    try:
        t = dtypes.minus(_single(raw), bg)
        t = clamp(np.asarray(t, dtype=np.float32), lo=0.0)
        t = dtypes.rdivide(t, flat)
        t = dtypes.plus(t, bg if const_offset is None else np.asarray(const_offset))
    except TypeError as err:
        # an integer offset with single data
        if "MATLAB:mixedClasses" not in str(err):
            raise
        raise MatlabError(
            "MATLAB:mixedClasses",
            "Integers can only be combined with integers of the same class, or scalar doubles.",
        ) from None
    except ValueError as err:
        if "broadcast" not in str(err):
            raise
        raise MatlabError(
            "MATLAB:sizeDimensionsMustMatch", "Arrays have incompatible sizes for this operation."
        ) from None
    t = np.asarray(t)
    if cast_data_type:
        if cls == "logical" and np.isnan(t).any():
            raise MatlabError("MATLAB:nologicalnan", "NaN's cannot be converted to logicals.")
        t = np.asarray(dtypes.cast(t, cls))
    return t


def apply_flat_field(
    raw: np.ndarray,
    flat: np.ndarray,
    background: np.ndarray,
    *,
    const_offset: Any = None,
    cast_data_type: bool = True,
    mex_only: bool = False,
    num_threads: int | None = None,
) -> np.ndarray:
    """The correction of ``XR_LSFlatFieldCorrection`` once ``flat`` (``LSImage``) is prepared,
    with a ``(z, y, x)`` ``raw`` and ``(y, x)`` images:

        LSImage = single(LSImage); background = single(background);
        try
            Temp = flat_field_correction_3d_mex(Rawdata, LSImage, background, ...
                                                castDataType[, constOffset]);
        catch
            <MATLAB's code>

    The mex for uint16 and single data; for every other class (and wherever the mex raises),
    MATLAB's code: ``single(raw) - background``, ``max(., 0)``, ``./ flat``, ``+ const_offset``
    (any numeric, as MATLAB adds it) or ``+ background``, then ``cast`` to the class of ``raw``
    (rounding half away from zero with saturation) when ``cast_data_type``. MATLAB's errors carry
    its identifiers. ``mex_only``: the mex alone (its ``MatlabError``).
    """
    f = _single(flat)
    b = _single(background)
    a = np.asarray(raw)
    args: list[Any] = [a, f, b, bool(cast_data_type)]
    if const_offset is not None and np.asarray(const_offset).size:
        args.append(const_offset)
    try:
        return flat_field_correction_3d_mex(*args, num_threads=num_threads)
    except MatlabError:
        if mex_only:
            raise
    off = const_offset if (const_offset is not None and np.asarray(const_offset).size) else None
    return _matlab_code(a, f, b, off, bool(cast_data_type))
