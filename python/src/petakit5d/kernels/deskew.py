"""The deskew/rotate kernels: the three mex files of ``deskew_rotate/mex`` (K03).

- ``volume_deskew_rotate_warp_mex(vol, tmat, bbox, save16bit)``:
  ``volume_deskew_rotate_warp(vol, tmat, box, save_16bit)``
- ``skewed_space_interp_volume_deskew_rotate_warp_mex(vol, xstep, stepsize, Reverse, tmat,
  bbox, save16bit)``: ``skewed_space_interp_volume_deskew_rotate_warp(vol, xstep, stepsize,
  reverse, tmat, box, save_16bit)``
- ``skewed_space_interp_defined_stepsize_mex(vol, xstep, stepsize, Reverse, save16bit)``:
  ``skewed_space_interp_defined_stepsize(vol, xstep, stepsize, reverse, save_16bit)``

Volumes are ``(z, y, x)`` C-ordered single or uint16 arrays (an image ``(y, x)`` is one slice);
results are ``(z, y, x)``, single, or uint16 for uint16 input with ``save_16bit``. ``box`` is the
output region ``(z0, y0, x0, z1, y1, x1)``, 0-based and end-exclusive, of integers (MATLAB's
``bbox`` ``[ymin, xmin, zmin, ymax, xmax, zmax]`` is ``(zmin - 1, ymin - 1, xmin - 1, zmax, ymax,
xmax)``); rows y past the input stay 0, as in the mex. ``tmat`` keeps MATLAB's convention: the
4 x 4 backward transform the MATLAB code builds (``M(i, j)`` is ``tmat[i - 1, j - 1]``, axes
x, y, z of MATLAB's array as there). ``xstep`` and ``stepsize`` are the mex's (raw slices
``xstep`` pixels apart in x; the grid step in raw slices).

Every voxel has the mex's arithmetic: the same operands, weights and order of operations,
including the uint16 rules (pairwise sums in the first ``floor(n / 16) * 16`` values of a y
line, ``(uint16)(v + 0.5f)`` output, not MATLAB's ``uint16``). As in the mex, only a logical true
(``True``, ``np.True_``, a 1-element bool array) counts as true for ``reverse`` and
``save_16bit``. Argument errors are the mex's, with MATLAB's identifier on the exception
(``identifier``): ``volume_warp:inputError`` (``ValueError``) and
``feather_blending:dataTypeError`` (``TypeError``). ``num_threads`` (None: all CPUs the process
may use) splits the work; results do not depend on it.

Two implementations give the same bits, and the data's layout picks one: the functions above
run kernels made for ``(z, y, x)`` volumes; the ``*_matlab_layout`` functions, for code that
holds MATLAB's layout already (a tiff read y-fastest, a write from that layout), take and return
C-ordered ``(z, x, y)`` arrays (MATLAB's ``(y, x, z)`` in memory) and run the mex's loops, at the
mex's speed. Neither turns a volume around.

The compiled module is ``petakit5d._ext.deskew``: AVX2 on x86-64 (as the mex), elsewhere the same
operations in the same order (``simd()`` says which), so every platform gives the x86 mex's bits.
"""

from __future__ import annotations

import math
import numbers
from typing import Any

import numpy as np

from petakit5d.core.threads import resolve_threads

__all__ = [
    "interp_grid_size",
    "skewed_space_interp_defined_stepsize",
    "skewed_space_interp_defined_stepsize_matlab_layout",
    "skewed_space_interp_volume_deskew_rotate_warp",
    "skewed_space_interp_volume_deskew_rotate_warp_matlab_layout",
    "volume_deskew_rotate_warp",
    "volume_deskew_rotate_warp_matlab_layout",
]

INPUT_ERROR = "volume_warp:inputError"
TYPE_ERROR = "feather_blending:dataTypeError"
_LIMIT = 9.0e15  # the mex's bound for finite values


def _ext(module: str = "deskew"):
    from petakit5d._ext import load

    return load(module)


def _error(message: str, identifier: str = INPUT_ERROR, cls: type = ValueError) -> Exception:
    err = cls(message)
    err.identifier = identifier
    return err


def _logical_true(x: Any) -> bool:
    """mxIsLogicalScalarTrue: a logical scalar that is true (a 1 is not)."""
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, np.ndarray) and x.dtype == np.bool_ and x.size == 1:
        return bool(x.reshape(-1)[0])
    return False


def _double_scalar(x: Any) -> float | None:
    """x if it is a double scalar as the mex requires (Python int or float, float64), else None"""
    if isinstance(x, (bool, np.bool_)):
        return None
    if isinstance(x, (int, float)) and not isinstance(x, np.generic):
        return float(x)
    if isinstance(x, (np.ndarray, np.generic)):
        a = np.asarray(x)
        if a.dtype == np.float64 and a.size == 1:
            return float(a.reshape(-1)[0])
    return None


def _volume(vol: Any) -> np.ndarray:
    """the input as a 3D (z, y, x) C-ordered view; errors as the mex"""
    if not isinstance(vol, np.ndarray):
        raise TypeError(f"vol must be a NumPy array, got {type(vol).__name__}")
    if vol.ndim > 3:
        raise _error("The input must have at most 3 dimensions")
    if not vol.flags.c_contiguous:
        raise ValueError("vol must be C-ordered (z, y, x); use core.arrays.as_volume")
    if vol.ndim < 3:
        vol = vol.reshape((1,) * (3 - vol.ndim) + vol.shape)
    return vol


def _check_dtype(vol: np.ndarray) -> None:
    if vol.dtype not in (np.float32, np.uint16):
        raise _error("Data type not suppported", TYPE_ERROR, TypeError)


def _scalars(xstep: Any, stepsize: Any) -> tuple[float, float]:
    vals = [_double_scalar(v) for v in (xstep, stepsize)]
    if vals[0] is None or vals[1] is None:
        raise _error("xstep and the step size must be double scalars")
    return vals[0], vals[1]


def _check_steps(xstep: float, stepsize: float) -> None:
    # the column ranges of the interpolated slices need xstep > 0
    if not (0 < xstep < _LIMIT) or not (0 < stepsize < _LIMIT):
        raise _error("xstep and the step size must be positive and finite")


def _tmat(tmat: Any) -> np.ndarray:
    a = tmat if isinstance(tmat, np.ndarray) else np.asarray(tmat, dtype=np.float64)
    if a.ndim != 2 or a.shape != (4, 4):
        raise _error("Input transformation matrix must be 4x4!")
    if a.dtype != np.float64:
        raise _error("Input transformation matrix must be a full double array")
    # NaN coordinates would pass the range checks of the kernels (finite and bounded, so that
    # the coordinates can't overflow to Inf - Inf)
    if not np.all(np.abs(a) < _LIMIT):
        raise _error("Input transformation matrix must be finite")
    return np.ascontiguousarray(a)


def _box(box: Any, ny: int) -> tuple[tuple[float, float, float], tuple[int, int, int], int]:
    """the mex's single precision geometry of a box: the starts (y, x, z: the mex's startX,
    startY, startZ), the output shape (z, y, x) and the rows y computed per line"""
    v = list(box) if not isinstance(box, np.ndarray) else box.reshape(-1).tolist()
    if len(v) != 6:
        raise _error("Input range for bbox is not 6")
    if not all(isinstance(b, numbers.Real) and not isinstance(b, bool) for b in v):
        raise _error("Input bbox must be finite")
    if not all(abs(b) < _LIMIT for b in v):
        raise _error("Input bbox must be finite")
    if any(b != math.floor(b) for b in v):
        raise _error("The box values must be integers")
    z0, y0, x0, z1, y1, x1 = (int(b) for b in v)
    # as the mex: MATLAB's 1-based start in single precision, minus 1; the end in single
    one = np.float32(1)
    start = [np.float32(b + 1) - one for b in (y0, x0, z0)]
    end = [np.float32(b) for b in (y1, x1, z1)]
    if any(s < 0 for s in start) or any(s > e for s, e in zip(start, end, strict=True)):
        raise _error(
            "Input bbox is invalid, the numbers must be nonnegative and the end ranges must be "
            "nonless than the start ranges."
        )
    oy, ox, oz = (int(e - s) for s, e in zip(start, end, strict=True))
    # output row i reads input row y0 + i: the rows past the input are outside and stay 0
    ys = math.floor(float(start[0]))
    line_len = min(oy, ny - ys) if ys < ny else 0
    return (float(start[0]), float(start[1]), float(start[2])), (oz, oy, ox), line_len


def _out_dtype(vol: np.ndarray, save_16bit: bool) -> type:
    return np.uint16 if vol.dtype == np.uint16 and save_16bit else np.float32


def interp_grid_size(nz: int, stepsize: float) -> int:
    """Slices of the skewed space interpolation grid of ``nz`` raw slices at ``stepsize``:
    ``floor(round((nz - 1) / stepsize * 100000) / 100000) + 1`` (0 for no slices)."""
    return int(_ext().interp_grid_size(int(nz), float(stepsize)))


# -- running a kernel ---------------------------------------------------------------------------
#
# Two implementations with the same bits; the data's layout picks one, so no volume is turned
# around:
# - (z, y, x) volumes: the C-layout kernels (lanes along y in 8 x 8-transposed scratch).
#   Measured faster on (z, y, x) data than the MATLAB-layout kernels with a threaded transpose
#   before and after, at 1, 4 and 24 threads (2.1 s against 2.7 s and 0.17 s against 0.51 s
#   for the fused uint16 kernel on a demo tile), since the transposes add two passes and two
#   fresh volumes.
# - MATLAB-layout data ((z, x, y) C arrays, y fastest: a tiff read y-fastest, a write from
#   that layout): the mex's loops (ml/ of the extension), at the mex's speed.
# ``_kernel_layout`` forces one ("c" or "matlab", with transposes as needed: for the tests).


def _swap(a: np.ndarray, num_threads: int) -> np.ndarray:
    """a with its last two axes swapped, C-contiguous (in threads)"""
    out = np.empty((a.shape[0], a.shape[2], a.shape[1]), dtype=a.dtype)
    _ext().swap_axes(a, out, num_threads)
    return out


def _run(name, vol, matlab_in, pre, post, out_zyx, dtype, num_threads, matlab_out, layout):
    """kernel `name` on vol ((z, y, x), or (z, x, y) with matlab_in); the output (z, y, x), or
    (z, x, y) with matlab_out"""
    ext = _ext()
    nt = resolve_threads(num_threads)
    oz, oy, ox = out_zyx
    if layout == "c":  # (z, y, x) data only: the *_matlab_layout functions never force it
        out = np.zeros(out_zyx, dtype=dtype)
        getattr(ext, name)(vol, *pre, out, *post, nt)
        return out
    v = vol if matlab_in else _swap(vol, nt)
    out = np.zeros((oz, ox, oy), dtype=dtype)
    getattr(ext, name + "_ml")(v, *pre, out, *post, nt)
    return out if matlab_out else _swap(out, nt)


def _layout(kernel_layout: str | None, matlab_in: bool) -> str:
    """the kernels to run: those of the data's layout unless forced"""
    if kernel_layout not in (None, "matlab", "c"):
        raise ValueError(f"kernel layout 'matlab' or 'c', not {kernel_layout!r}")
    return kernel_layout or ("matlab" if matlab_in else "c")


def _matlab_volume(vol: Any) -> tuple[np.ndarray, tuple[int, int, int]]:
    """MATLAB-layout data, the C array (z, x, y), and its (z, y, x) shape"""
    if not isinstance(vol, np.ndarray):
        raise TypeError(f"vol must be a NumPy array, got {type(vol).__name__}")
    if vol.ndim != 3:
        raise _error("The input must have at most 3 dimensions")
    if not vol.flags.c_contiguous:
        raise ValueError("vol must be C-contiguous (z, x, y): MATLAB's layout")
    return vol, (vol.shape[0], vol.shape[2], vol.shape[1])


def _warp(vol, shape, matlab_in, matlab_out, tmat, box, save_16bit, num_threads, layout):
    t = _tmat(tmat)
    start, out_shape, line_len = _box(box, shape[1])
    s16 = _logical_true(save_16bit)
    _check_dtype(vol)
    return _run(
        "warp",
        vol,
        matlab_in,
        (t,),
        (*start, line_len),
        out_shape,
        _out_dtype(vol, s16),
        num_threads,
        matlab_out,
        _layout(layout, matlab_in),
    )


def _fused(
    vol,
    shape,
    matlab_in,
    matlab_out,
    xstep,
    stepsize,
    reverse,
    tmat,
    box,
    save_16bit,
    num_threads,
    layout,
):
    xs, step = _scalars(xstep, stepsize)
    rev = _logical_true(reverse)
    s16 = _logical_true(save_16bit)
    _check_steps(xs, step)
    # the columns of an interpolated slice come from columns up to ceil(xstep) + 1 apart
    if shape[2] < math.ceil(xs) + 2:
        raise _error("The second dimension of the input must be at least ceil(xstep) + 2")
    if not _ext().interp_grid_size(shape[0], step) < _LIMIT:
        raise _error("The step size is too small")
    t = _tmat(tmat)
    start, out_shape, line_len = _box(box, shape[1])
    _check_dtype(vol)
    return _run(
        "skewed_warp",
        vol,
        matlab_in,
        (xs, step, rev, t),
        (*start, line_len),
        out_shape,
        _out_dtype(vol, s16),
        num_threads,
        matlab_out,
        _layout(layout, matlab_in),
    )


def _interp(
    vol, shape, matlab_in, matlab_out, xstep, stepsize, reverse, save_16bit, num_threads, layout
):
    xs, step = _scalars(xstep, stepsize)
    xs32 = float(np.float32(xs))
    rev = _logical_true(reverse)
    s16 = _logical_true(save_16bit)
    _check_steps(xs32, step)
    if shape[2] < math.ceil(xs32) + 2:
        raise _error("The second dimension of the input must be at least ceil(xstep) + 2")
    nint = _ext().interp_grid_size(shape[0], step)
    if not nint < _LIMIT:
        raise _error("The step size is too small")
    _check_dtype(vol)
    out_shape = (int(nint), shape[1], shape[2])
    return _run(
        "interp_defined_stepsize",
        vol,
        matlab_in,
        (xs32, step, rev),
        (),
        out_shape,
        _out_dtype(vol, s16),
        num_threads,
        matlab_out,
        _layout(layout, matlab_in),
    )


# -- the kernels on (z, y, x) volumes -------------------------------------------------------------


def volume_deskew_rotate_warp(
    vol: np.ndarray,
    tmat: Any,
    box: Any,
    save_16bit: Any = False,
    *,
    num_threads: int | None = None,
    _kernel_layout: str | None = None,
) -> np.ndarray:
    """``volume_deskew_rotate_warp_mex``: deskew/rotate warp of ``vol`` (no interpolation in
    skewed space) into the output region ``box``, with the backward transform ``tmat``."""
    vol3 = _volume(vol)
    return _warp(vol3, vol3.shape, False, False, tmat, box, save_16bit, num_threads, _kernel_layout)


def skewed_space_interp_volume_deskew_rotate_warp(
    vol: np.ndarray,
    xstep: Any,
    stepsize: Any,
    reverse: Any,
    tmat: Any,
    box: Any,
    save_16bit: Any = False,
    *,
    num_threads: int | None = None,
    _kernel_layout: str | None = None,
) -> np.ndarray:
    """``skewed_space_interp_volume_deskew_rotate_warp_mex``: skewed space interpolation at
    ``stepsize`` and the deskew/rotate warp in one pass, without an interpolated volume (the
    fused kernel)."""
    vol3 = _volume(vol)
    return _fused(
        vol3,
        vol3.shape,
        False,
        False,
        xstep,
        stepsize,
        reverse,
        tmat,
        box,
        save_16bit,
        num_threads,
        _kernel_layout,
    )


def skewed_space_interp_defined_stepsize(
    vol: np.ndarray,
    xstep: Any,
    stepsize: Any,
    reverse: Any,
    save_16bit: Any = False,
    *,
    num_threads: int | None = None,
    _kernel_layout: str | None = None,
) -> np.ndarray:
    """``skewed_space_interp_defined_stepsize_mex``: skewed space interpolation of ``vol`` to
    ``interp_grid_size(nz, stepsize)`` slices (xstep in single precision, as in the mex)."""
    vol3 = _volume(vol)
    return _interp(
        vol3,
        vol3.shape,
        False,
        False,
        xstep,
        stepsize,
        reverse,
        save_16bit,
        num_threads,
        _kernel_layout,
    )


# -- the same on MATLAB's layout (file-level code that reads and writes that layout) -------------


def volume_deskew_rotate_warp_matlab_layout(
    vol_zxy: np.ndarray,
    tmat: Any,
    box: Any,
    save_16bit: Any = False,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """``volume_deskew_rotate_warp`` on MATLAB's layout: ``vol_zxy`` and the result are C arrays
    ``(z, x, y)`` (y fastest, MATLAB's ``(y, x, z)`` in memory); ``box`` is ours, ``(z0, y0, x0,
    z1, y1, x1)``."""
    v, shape = _matlab_volume(vol_zxy)
    return _warp(v, shape, True, True, tmat, box, save_16bit, num_threads, None)


def skewed_space_interp_volume_deskew_rotate_warp_matlab_layout(
    vol_zxy: np.ndarray,
    xstep: Any,
    stepsize: Any,
    reverse: Any,
    tmat: Any,
    box: Any,
    save_16bit: Any = False,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """``skewed_space_interp_volume_deskew_rotate_warp`` on MATLAB's layout (``(z, x, y)`` C
    arrays in and out)."""
    v, shape = _matlab_volume(vol_zxy)
    return _fused(
        v, shape, True, True, xstep, stepsize, reverse, tmat, box, save_16bit, num_threads, None
    )


def skewed_space_interp_defined_stepsize_matlab_layout(
    vol_zxy: np.ndarray,
    xstep: Any,
    stepsize: Any,
    reverse: Any,
    save_16bit: Any = False,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """``skewed_space_interp_defined_stepsize`` on MATLAB's layout (``(z, x, y)`` C arrays in
    and out)."""
    v, shape = _matlab_volume(vol_zxy)
    return _interp(v, shape, True, True, xstep, stepsize, reverse, save_16bit, num_threads, None)
