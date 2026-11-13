"""Richardson-Lucy deconvolution on the CPU: ``decon_lucy`` (decon_lucy_function, the
'simplified' RL with Biggs-Andrews acceleration) and ``decon_lucy_omw`` (decon_lucy_omw_function,
RL with an OTF-masked Wiener back projector).

Arrays are the port's: volumes ``(z, y, x)``, images ``(y, x)``, C-contiguous; a volume with one z
slice is MATLAB's 2D image. Results equal MATLAB's (decon_lucy_function.m and
decon_lucy_omw_function.m at the head of ``test``, on the CPU):

- **The iterations** run in the half-spectrum kernels (``petakit5d.kernels.decon.LucyRFFT``,
  MATLAB's decon_lucy_rfft_mex, its default path on Linux) for data MATLAB's mex takes: real, not
  double, 2D or 3D with at least 2 elements along y and x, a scalar background. On the exact FFT
  backend they are MATLAB's bits (the kernels keep RL's arrays in MATLAB's layout); on the
  permissive one (the default) within decision record 0006's tolerance. Other data (double, a
  single row or column, a background per voxel) or ``use_mex=False`` run MATLAB's loop on full
  spectra (``matlabcompat.fft``), with MATLAB's mixed-class arithmetic.
- **RL's acceleration**: lambda from sums in double in a fixed order (blocks of 4,096, four lanes:
  ``kernels.decon.accel_lambda``) of the kernels' internal arrays: MATLAB's linear order in MATLAB's
  layout (the exact backend: MATLAB's bits), the arrays' own order on the permissive backend
  (faster; last bits as the backend's), clamped to [0, 1]; the prediction ``max(J_2 + lambda *
  (J_2 - J_3), 0)``; their arrays are reused across iterations. MATLAB's loop (data the kernels
  don't take) sums in MATLAB's order. The PSF is normalized by its fixed-order
  double sum (decon_sum_double). Nothing depends on the thread count.
- **The OTF cache**: the last OTF(s) of each method are kept for the same PSF bits (and back
  projector bits), data size, path and FFT backend and layout (``clear_otf_cache`` empties it).
- **Post-processing** as MATLAB: ``dampFactor`` (> 1), the mask of nonzero input voxels eroded by
  ``EdgeErosion`` (``decon.erosion``), ``scaleFactor`` and ``deconOffset`` (the offset and the
  erosion multiply by the mask), uint16 output (``save16bit``), the box ``bbox``.
- **Debug output** (``debug``): the error metric per iteration (``err_mat``: iteration, mean
  squared difference of the estimate and the data in MATLAB's order at one thread, its ratio to
  the first one, MATLAB's stop criterion) and the estimate as ``<debug_folder>/Iter_%04d.tif``
  every ``save_step`` iterations, with MATLAB's bookkeeping and its errors (a step past the
  iterations or one that isn't an integer can leave ``istp`` undefined).

``use_gpu``: False runs on the CPU; None (MATLAB's ``[]``) on a GPU if one is usable, else on the
CPU with a ``GPUFallbackWarning`` (once per new OTF, as MATLAB warns); True requires a usable GPU
(``NoUsableGPUError``, identifier ``decon:noUsableGPU``). The GPU path itself is
``petakit5d.decon.rl_gpu`` (task R04); where it isn't installed, None runs on the CPU with a
warning and True raises.
"""

from __future__ import annotations

import math
import os
import threading
import warnings
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any, NamedTuple

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "DeconError",
    "OMWResult",
    "RLResult",
    "clear_otf_cache",
    "decon_lucy",
    "decon_lucy_omw",
    "lucy",
    "lucy_omw",
    "use_rfft_kernels",
]

_DBL_EPS = 2.220446049250313e-16
_UNDEFINED_ISTP = "MATLAB:UndefinedFunction"
_NONINTEGER_SIZE = "MATLAB:NonIntegerInput"


class DeconError(ValueError):
    """An error MATLAB raises, with its identifier ('' where MATLAB's error has none)."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


class RLResult(NamedTuple):
    """decon_lucy_function's outputs: the result, err_mat (NUMIT x 4; with debug output the rows
    MATLAB keeps), the number of iterations run, and the output of ``overlap_fcn`` (or None)."""

    decon: np.ndarray
    err_mat: np.ndarray
    iterations: int
    overlap_out: Any = None


class OMWResult(NamedTuple):
    """decon_lucy_omw_function's outputs: the result, err_mat, the output of ``overlap_fcn``."""

    decon: np.ndarray
    err_mat: np.ndarray
    overlap_out: Any = None


# ---------------------------------------------------------------------------------------------
# small MATLAB rules
# ---------------------------------------------------------------------------------------------


def _threads(threads: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(threads)


def _matlab_2d(a: np.ndarray) -> np.ndarray:
    """a (1, ny, nx) volume as MATLAB's 2D array (y, x): MATLAB's size has no third dimension"""
    return a.reshape(a.shape[1:]) if a.ndim == 3 and a.shape[0] == 1 else a


def _data(i: Any) -> np.ndarray:
    a = np.asarray(i)
    if a.dtype.kind not in "iuf":
        raise TypeError(f"the data must be a real numeric array, not {a.dtype}")
    if a.ndim not in (2, 3):
        raise ValueError("the data must be a (z, y, x) volume or a (y, x) image")
    return _matlab_2d(a)


def _single(x: Any) -> np.ndarray:
    """MATLAB's single(x) of a small array (PSFs)"""
    from petakit5d.core import dtypes

    a = np.asarray(x)
    if a.dtype == np.float32:
        return np.ascontiguousarray(a)
    return np.ascontiguousarray(dtypes.single(a))


def _psf(p: Any) -> np.ndarray:
    a = _single(p)
    if a.ndim == 0:
        a = a.reshape(1, 1)
    if a.ndim == 1:
        a = a.reshape(1, -1)  # a row
    if a.ndim > 3:
        raise ValueError("the PSF must be a (z, y, x) volume or a (y, x) image")
    return np.ascontiguousarray(_matlab_2d(a))


def _iterations(numit: Any) -> int:
    """the iterations of ``for k = 1 : NUMIT`` with MATLAB's error of ``zeros(NUMIT, 4)``"""
    a = np.asarray(numit)
    if a.size != 1 or a.dtype.kind not in "biuf":
        raise DeconError("MATLAB:NonScalarInput", "Size inputs must be scalar.")
    x = float(a.reshape(-1)[0])
    if not math.isfinite(x) or x != math.floor(x):
        raise DeconError(_NONINTEGER_SIZE, "Size inputs must be integers.")
    return max(0, int(x))


def _background(background: Any) -> tuple[bool, Any]:
    """(scalar, value): MATLAB's [] is 0; a real scalar is a float (a MATLAB double) or the NumPy
    scalar of its class; anything else the array (MATLAB's form of a 2D one)"""
    if background is None:
        return True, 0.0
    if isinstance(background, (int, float)) and not isinstance(background, bool):
        return True, float(background)
    b = np.asarray(background)
    if b.size == 0:
        return True, 0.0
    if b.size == 1 and not np.iscomplexobj(b):
        v = b.reshape(-1)[0]
        return True, float(v) if b.dtype == np.float64 else v
    return False, _matlab_2d(b)


def use_rfft_kernels(
    use_mex: bool | None, method: str, data: np.ndarray, debug: bool, background: Any
) -> tuple[bool, bool]:
    """(use, can): whether the iterations run in the half-spectrum kernels (decon_use_rfft_mex on
    the CPU): ``can`` for real, non-empty, non-double 2D or 3D data with at least 2 elements
    along y and x and a real scalar background; None uses them where they can, True requires
    them (else an error), False runs MATLAB's loop."""
    a = _matlab_2d(np.asarray(data))
    scalar_bg = background is None or (
        np.size(background) == 1 and not np.iscomplexobj(np.asarray(background))
    )
    can = (
        a.dtype != np.float64
        and a.ndim in (2, 3)
        and a.shape[-2] > 1
        and a.shape[-1] > 1
        and not np.iscomplexobj(a)
        and a.size > 0
        and scalar_bg
    )
    if method not in ("simplified", "omw"):
        raise ValueError(f"unknown method {method!r}")
    if use_mex is None:
        return can, can
    if use_mex and not can:
        from petakit5d.core.arrays import matlab_class, size_to_matlab

        msize = "[" + " ".join(str(v) for v in size_to_matlab(a.shape)) + "]"
        raise DeconError(
            "",
            f"decon_lucy_rfft_mex does not support the data (class {matlab_class(a)}, size "
            f"{msize}) or the options (debug on the GPU or a non-scalar background).",
        )
    return bool(use_mex), can


# ---------------------------------------------------------------------------------------------
# threaded element-wise steps
# ---------------------------------------------------------------------------------------------


def _chunks(n: int, threads: int) -> list[tuple[int, int]]:
    """pieces of 0 .. n for the threads (at least 1 Mi elements each)"""
    k = max(1, min(4 * threads, n // (1 << 20)))
    step = -(-n // k)
    return [(s, min(n, s + step)) for s in range(0, n, step)]


def _parallel(fn: Callable[[int, int], None], n: int, threads: int) -> None:
    """fn(start, stop) over pieces of 0 .. n in threads (NumPy releases the interpreter lock)"""
    parts = _chunks(n, threads)
    if len(parts) <= 1 or threads <= 1:
        for s, e in parts:
            fn(s, e)
        return
    with ThreadPoolExecutor(min(threads, len(parts))) as ex:
        list(ex.map(lambda p: fn(*p), parts))


def _damp(j: np.ndarray, i: np.ndarray, damp: float, threads: int) -> np.ndarray:
    """``J .* (J <= d * I) + min(J, d * I .* (J >= d * I))``: d * I in double rounded to the class
    (single data), MATLAB's min (NaN gives the other value, a tie J)"""
    from petakit5d.kernels import ops

    if j.dtype == np.float32 and i.dtype == np.float32:
        di = ops.scale_offset(i, float(damp), num_threads=threads)
    else:
        from petakit5d.core import dtypes

        di = np.asarray(dtypes.times(float(damp), i))
    jf, df = j.reshape(-1), np.ascontiguousarray(di).reshape(-1)
    out = np.empty_like(jf)

    def piece(s: int, e: int) -> None:
        a, b = jf[s:e], df[s:e]
        with np.errstate(invalid="ignore"):
            lo = a * (a <= b).astype(a.dtype)
            hi = b * (a >= b).astype(a.dtype)
            pick = (hi < a) | (np.isnan(a) & ~np.isnan(hi))
            out[s:e] = lo + np.where(pick, hi, a)

    _parallel(piece, jf.size, threads)
    return out.reshape(j.shape)


def _subtract_background(
    raw: np.ndarray, scalar: bool, bg: Any, threads: int, inplace: bool = False
) -> np.ndarray:
    """MATLAB's ``I = max(I - Background, 0)`` of the data (single for every class but double):
    a new array, or ``raw`` itself with ``inplace`` (a single array of the caller's own)"""
    from petakit5d.kernels import ops

    if raw.dtype == np.float64 or raw.dtype == np.float32:
        x = raw
    else:
        x = ops.convert(raw, "single", num_threads=threads)
        inplace = True
    if scalar and isinstance(bg, float):
        return ops.scale_offset(
            x, None, -bg, lo=0.0, out=x if inplace else None, num_threads=threads
        )
    from petakit5d.core import dtypes

    return ops.clamp(np.asarray(dtypes.minus(x, bg)), lo=0.0, num_threads=threads)


def _single_input(raw: np.ndarray, threads: int) -> np.ndarray:
    """MATLAB's single(I) for data that isn't double"""
    from petakit5d.kernels import ops

    if raw.dtype in (np.float32, np.float64):
        return raw
    return ops.convert(raw, "single", num_threads=threads)


# ---------------------------------------------------------------------------------------------
# MATLAB's loop on full spectra (data the kernels don't take, or use_mex=False)
# ---------------------------------------------------------------------------------------------


def _cmul(h: np.ndarray, s: np.ndarray, conj: bool = False) -> np.ndarray:
    """MATLAB's ``h .* s`` (``conj(h) .* s``) of spectra: each product and sum one IEEE operation
    (no fused multiply-add), in single when both are single, else in double and rounded to single
    if either is single; a real operand multiplies both parts of the other"""
    singles = (np.float32, np.complex64)
    both = h.dtype in singles and s.dtype in singles
    ft = np.float32 if both else np.float64
    ct = np.complex64 if h.dtype in singles or s.dtype in singles else np.complex128
    out = np.empty(np.broadcast_shapes(h.shape, s.shape), np.complex64 if both else np.complex128)
    hc, sc = np.iscomplexobj(h), np.iscomplexobj(s)
    hr = h.real.astype(ft)
    sr = s.real.astype(ft)
    if hc and sc:
        hi, si = h.imag.astype(ft), s.imag.astype(ft)
        if conj:
            out.real = hr * sr + hi * si
            out.imag = hr * si - hi * sr
        else:
            out.real = hr * sr - hi * si
            out.imag = hr * si + hi * sr
    elif hc:
        hi = h.imag.astype(ft)
        out.real = hr * sr
        out.imag = -(hi * sr) if conj else hi * sr
    elif sc:
        out.real = hr * s.real.astype(ft)
        out.imag = hr * s.imag.astype(ft)
    else:
        return (hr * sr).astype(ct)
    return out.astype(ct) if out.dtype != ct else out


def _psf2otf(psf: np.ndarray, shape: Sequence[int], backend: str | None, threads: int):
    """decon_psf2otf(psf, outSize) (``decon.otf.psf2otf``)"""
    from petakit5d.decon.otf import psf2otf

    return psf2otf(psf, shape, backend=backend, threads=threads)


def _conv(x: np.ndarray, h: np.ndarray, backend, threads, conj: bool = False) -> np.ndarray:
    """``real(ifftn(h .* fftn(x)))`` (``conj(h)``) as MATLAB computes it"""
    from petakit5d.matlabcompat import fft as F

    s = F.fftn(x, backend=backend, threads=threads)
    r = F.ifftn(_cmul(h, s, conj), backend=backend, threads=threads)
    return np.ascontiguousarray(r.real) if np.iscomplexobj(r) else r


def _like(x: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """MATLAB's cast(x, 'like', ref)"""
    if x.dtype == ref.dtype:
        return x
    from petakit5d.core import dtypes

    return np.asarray(dtypes.cast(x, "single" if ref.dtype == np.float32 else "double"))


def _clamp0(x: np.ndarray, threads: int, lo: float = 0.0) -> np.ndarray:
    from petakit5d.kernels import ops

    return ops.clamp(x, lo=lo, num_threads=threads)


def _lambda_any(j2, y, j4, threads):
    """decon_lucy_accel('lambda'): the kernels for arrays of one class, else MATLAB's code"""
    from petakit5d.kernels import decon as K

    if j2.dtype == y.dtype == j4.dtype and j2.dtype in (np.float32, np.float64):
        r = K.accel_lambda(j2, y, j4, threads=threads)
        return r.lam, r.j4
    from petakit5d.core import dtypes

    d = np.asarray(dtypes.minus(j2, y))
    s0 = K.sum_double(np.asarray(dtypes.times(d, j4)), threads=threads)
    s1 = K.sum_double(np.asarray(dtypes.times(j4, j4)), threads=threads)
    with np.errstate(all="ignore"):
        lam = float(np.float64(s0) / (np.float64(s1) + _DBL_EPS))
    return K._clamp(lam), d


def _predict_any(j2, j3, lam, threads, out=None):
    """decon_lucy_accel('predict'): the kernel for arrays of one class, else MATLAB's code"""
    from petakit5d.kernels import decon as K

    if j2.dtype == j3.dtype and j2.dtype in (np.float32, np.float64):
        return K.accel_predict(j2, j3, lam, threads=threads, out=out)
    from petakit5d.core import dtypes

    t = dtypes.plus(j2, dtypes.times(lam, dtypes.minus(j2, j3)))
    return _clamp0(np.asarray(t), threads)


# ---------------------------------------------------------------------------------------------
# debug output
# ---------------------------------------------------------------------------------------------


def _rem(k: int, step: float) -> float:
    """MATLAB's rem(k, step)"""
    if step == 0 or math.isnan(step):
        return math.nan
    return math.fmod(k, step)


def _mmin(a: float, b: float) -> float:
    if math.isnan(a):
        return b
    if math.isnan(b):
        return a
    return a if a <= b else b


def _ratio(a: float, b: float) -> float:
    with np.errstate(all="ignore"):
        return float(np.float64(a) / np.float64(b))


def _mse(j: np.ndarray, i: np.ndarray, matlab_layout: bool, threads: int) -> float:
    """fixed_order_stat('mean', (double(J) - double(I)) .^ 2): the sum in MATLAB's linear order at
    one thread (decon_sum_double) divided by the count"""
    from petakit5d.kernels import decon as K

    d = j.astype(np.float64)
    np.subtract(d, i, out=d)
    np.multiply(d, d, out=d)
    return _ratio(K.sum_double(d, matlab_layout=matlab_layout, threads=threads), d.size)


class _Debug:
    """err_mat and the iteration files of the MATLAB loops' debug output"""

    def __init__(self, n_iter: int, step: Any, folder: str, threads: int):
        self.err = np.zeros((n_iter, 4))
        self.step = float(np.asarray(step).reshape(-1)[0])
        self.folder = folder
        self.istp: float | None = None
        self.threads = threads

    def _undefined(self) -> DeconError:
        return DeconError(_UNDEFINED_ISTP, "Unrecognized function or variable 'istp'.")

    def _row(self, k: float) -> int:
        """a 1-based row index of err_mat as MATLAB checks it"""
        if k != math.floor(k) or k < 1:
            raise DeconError(
                "MATLAB:badsubscript",
                "Index in position 1 is invalid. Array indices must be positive integers or "
                "logical values.",
            )
        if k > self.err.shape[0]:
            raise DeconError(
                "MATLAB:badsubscript",
                "Index in position 1 exceeds array bounds. Index must not exceed "
                f"{self.err.shape[0]}.",
            )
        return int(k) - 1

    def bookkeeping(self, k: int, err: float | None) -> None:
        e, s = self.err, self.step
        if err is not None:
            e[k - 1, 0:2] = (k, err)
        else:
            e[k - 1, 0] = k
        e[k - 1, 2] = _ratio(e[k - 1, 1], e[0, 1])
        if _rem(k, s) == 0 and k > s * 2:
            self.istp = k / s
            i = self._row(self.istp)
            i1, i2 = self._row(self.istp - 1), self._row(self.istp - 2)
            e[k - 1, 3] = _mmin(
                abs(e[i, 2] - e[i1, 2]) / 10, abs(e[i, 2] + e[i2, 2] - 2 * e[i1, 2]) / 2
            )
        elif k < 3 * s:
            e[k - 1, 3] = 1
        else:
            if self.istp is None:
                raise self._undefined()
            e[k - 1, 3] = e[self._row(self.istp * s), 3]

    def files(self, k: int, estimate: Callable[[], np.ndarray]) -> None:
        if k == 1:
            os.makedirs(self.folder, exist_ok=True)
        if _rem(k, self.step) == 0:
            from petakit5d.io import write_tiff
            from petakit5d.kernels import ops

            u = ops.convert(estimate(), "uint16", num_threads=self.threads)
            write_tiff(u, f"{self.folder}/Iter_{k:04d}.tif")

    def rows_to_istp(self) -> np.ndarray:
        """decon_lucy_function's ``err_mat(1 : istp, :)``"""
        if self.istp is None:
            raise self._undefined()
        self._row(self.istp)
        return self.err[: int(self.istp)].copy()

    def rows_kept(self, n_iter: int) -> np.ndarray:
        """decon_lucy_omw_function's ``err_mat(1 : floor(NUMIT / estep), :)``"""
        with np.errstate(all="ignore"):
            m = np.floor(np.float64(n_iter) / np.float64(self.step))
        if not m >= 1:
            return np.zeros((0, 4))
        self._row(m)
        return self.err[: int(m)].copy()


# ---------------------------------------------------------------------------------------------
# the OTF cache
# ---------------------------------------------------------------------------------------------


class _Entry(NamedTuple):
    key: tuple  # (data shape, kernels or loop, backend, layout, whether a GPU was checked)
    psfs: tuple[np.ndarray, ...]
    otfs: tuple[np.ndarray, ...]


_CACHE: dict[str, _Entry] = {}
_CACHE_LOCK = threading.Lock()


def clear_otf_cache() -> None:
    """Forget the cached OTFs (MATLAB's ``clear decon_lucy_function decon_lucy_omw_function``)."""
    with _CACHE_LOCK:
        _CACHE.clear()


def _same_psfs(a: tuple[np.ndarray, ...], b: tuple[np.ndarray, ...]) -> bool:
    """the same sizes and bits (-0 differs from 0, a NaN equals itself)"""
    return len(a) == len(b) and all(
        x.shape == y.shape and x.tobytes() == y.tobytes() for x, y in zip(a, b, strict=True)
    )


def _gpu_device(use_gpu: bool | None):
    """the GPU to run on, or None (the CPU): gpu_for's warning or NoUsableGPUError (with MATLAB's
    identifier)"""
    from petakit5d.gpu.device import NoUsableGPUError, gpu_for

    try:
        return gpu_for(use_gpu, purpose="the deconvolution")
    except NoUsableGPUError as e:
        e.identifier = "decon:noUsableGPU"
        raise


def _normalized(psf: np.ndarray, threads: int) -> np.ndarray:
    """``PSF ./ decon_sum_double(PSF)``: in double, rounded to single"""
    from petakit5d.core import dtypes
    from petakit5d.kernels import decon as K

    s = K.sum_double(psf, threads=threads)
    return np.ascontiguousarray(dtypes.rdivide(psf, s), dtype=np.float32)


def _otfs(
    method: str,
    psfs: tuple[np.ndarray, ...],
    shape: tuple[int, ...],
    kernels: bool,
    use_gpu: bool | None,
    backend: str,
    layout: str | None,
    threads: int,
):
    """(OTFs, GPU device): the cached OTFs of the same PSFs, data size, path, backend and layout,
    else new ones; the GPU check of a new OTF as MATLAB makes it (None: the CPU)"""
    want_gpu = use_gpu is not False
    key = (shape, kernels, backend, layout, want_gpu)
    with _CACHE_LOCK:
        e = _CACHE.get(method)
        hit = e is not None and e.key == key and _same_psfs(e.psfs, psfs)
        if hit and use_gpu is not True:
            return e.otfs, None
        _CACHE.pop(method, None)
    dev = _gpu_device(use_gpu) if want_gpu else None
    if dev is not None:
        return None, dev
    norm = [_normalized(p, threads) for p in psfs]
    if kernels:
        from petakit5d.kernels.decon import LucyRFFT

        k = LucyRFFT(shape, backend=backend, threads=threads, layout=layout)
        otfs = tuple(k.otf(p) for p in norm)
    else:
        otfs = tuple(_psf2otf(p, shape, backend, threads) for p in norm)
    with _CACHE_LOCK:
        _CACHE[method] = _Entry(key, tuple(p.copy() for p in psfs), otfs)
    return otfs, None


def _gpu_arguments(values: dict[str, Any]) -> dict[str, Any]:
    """the arguments of a call for the GPU path (the registered function's)"""
    return {k: v for k, v in values.items() if k not in ("threads", "backend", "layout")}


def _gpu_path(method: str, use_gpu: bool | None, kwargs: dict) -> Any:
    """the GPU path (task R04), or None for the CPU where it isn't installed"""
    try:
        from petakit5d.decon import rl_gpu  # type: ignore[attr-defined]
    except ImportError:
        msg = "The GPU path of the deconvolution (petakit5d.decon.rl_gpu) is not installed"
        if use_gpu:
            from petakit5d.gpu.device import NoUsableGPUError

            err = NoUsableGPUError(f"{msg}, but the deconvolution is to run on a GPU.")
            err.identifier = "decon:noUsableGPU"
            raise err from None
        from petakit5d.gpu.device import GPUFallbackWarning

        warnings.warn(f"{msg}, so the deconvolution runs on the CPU.", GPUFallbackWarning, 3)
        return None
    return getattr(rl_gpu, method)(**kwargs)


# ---------------------------------------------------------------------------------------------
# post-processing
# ---------------------------------------------------------------------------------------------


def _crop(j: np.ndarray, bbox: Sequence[Any]) -> np.ndarray:
    """crop3d: the box of the result (a 2D result takes a 6-value box whose z range is 0 .. 1)"""
    from petakit5d.kernels.crop import CROP_INPUT, MatlabError, crop_3d

    box = list(bbox)
    if j.ndim == 2 and len(box) == 6:
        if (box[0] or 0) != 0 or box[3] not in (None, 1):
            raise MatlabError(CROP_INPUT, f"bbox {tuple(box)} is outside the image")
        box = [box[1], box[2], box[4], box[5]]
    return crop_3d(j, box)


def _post(
    j: np.ndarray,
    damp_input: Callable[[], np.ndarray] | None,
    mask: np.ndarray | None,
    *,
    damp_factor: Any,
    scale_factor: Any,
    decon_offset: Any,
    edge_erosion: Any,
    save_16bit: bool,
    bbox: Any,
    shape_in: tuple[int, ...],
    threads: int,
) -> np.ndarray:
    """dampFactor, EdgeErosion, scaleFactor, deconOffset, save16bit and bbox as MATLAB"""
    from petakit5d.kernels import ops

    damp, scale, offset, erosion = (
        float(np.asarray(v).reshape(-1)[0])
        for v in (damp_factor, scale_factor, decon_offset, edge_erosion)
    )
    if damp > 1:
        j = _damp(j, _like(damp_input(), j), damp, threads)
    if erosion > 0:
        from petakit5d.decon.erosion import decon_mask_edge_erosion

        mask = decon_mask_edge_erosion(mask, edge_erosion, threads=threads)
    if scale != 1 or offset != 0 or erosion > 0:
        if offset != 0 or erosion > 0:
            j = ops.scale_offset(j, scale, offset, mask=mask, num_threads=threads)
        else:
            j = ops.scale_offset(j, scale, num_threads=threads)
    if save_16bit:
        j = ops.convert(j, "uint16", num_threads=threads)
    if len(shape_in) == 3 and j.ndim == 2:
        j = j.reshape((1, *j.shape))
    if bbox is not None and len(bbox) > 0:
        j = _crop(j, bbox)
    return j


def _mask_of(i: np.ndarray, decon_offset: Any, edge_erosion: Any, threads: int):
    """MATLAB's ``I_mask = I ~= 0`` where the offset or the erosion needs it"""
    from petakit5d.kernels import ops

    off = float(np.asarray(decon_offset).reshape(-1)[0])
    ero = float(np.asarray(edge_erosion).reshape(-1)[0])
    if off != 0 or ero > 0:
        return ops.compare(i, "ne", 0.0, num_threads=threads)
    return None


# ---------------------------------------------------------------------------------------------
# decon_lucy_function
# ---------------------------------------------------------------------------------------------


def _internal_zeros(k, n: int) -> np.ndarray:
    """zeros in a kernel object's internal layout"""
    return np.zeros(n, np.float32).reshape(k.internal_shape, order="F" if k.matlab_layout else "C")


def _lambda_step(j2, y, j4, matlab_layout: bool, threads: int) -> float:
    """decon_lucy_accel('lambda') on the kernels' internal arrays, the new J_4 written into J_4: the
    sums in the arrays' memory order, a fixed order at any thread count. In MATLAB's layout (the
    exact backend) that is MATLAB's linear order (MATLAB's bits); on the permissive backend's
    (z, y, x) arrays it is their own order, whose sums take a third of those in MATLAB's order
    (17-21% of an RL iteration against 6-13%, R03's measurement), and whose last bits are of the
    size of the backend's FFT differences (decision record 0006)."""
    from petakit5d.kernels import decon as K

    if matlab_layout:
        return K.accel_lambda(j2, y, j4, matlab_layout=True, threads=threads, out=j4).lam
    f4 = j4.reshape(-1)
    return K.accel_lambda(j2.reshape(-1), y.reshape(-1), f4, threads=threads, out=f4).lam


def lucy(
    i: np.ndarray,
    psf: np.ndarray,
    numit: Any,
    *,
    background: Any = None,
    use_gpu: bool | None = False,
    save_16bit: bool = True,
    damp_factor: Any = 1,
    scale_factor: Any = 1,
    decon_offset: Any = 0,
    edge_erosion: Any = 0,
    bbox: Sequence[Any] | None = None,
    debug: bool = False,
    debug_folder: str = "./debug/",
    save_step: Any = 1,
    use_mex: bool | None = None,
    overlap_fcn: Callable[[], Any] | None = None,
    threads: int | None = None,
    backend: str | None = None,
    layout: str | None = None,
) -> RLResult:
    """``decon_lucy`` with the thread count (None: the CPUs this process may use), the CPU FFT
    backend (None: ``petakit5d.fft.backend()``) and the kernels' layout (None: theirs) as
    arguments; ``use_gpu`` False by default (the CPU). Results don't depend on ``threads``."""
    call = _gpu_arguments(locals())
    from petakit5d import fft
    from petakit5d.kernels import decon as K

    t = _threads(threads)
    be = fft.backend() if backend is None else backend
    shape_in = np.shape(i)
    data = _data(i)
    scalar_bg, bg = _background(background)
    p = _psf(psf)
    use, _ = use_rfft_kernels(use_mex, "simplified", data, debug, None if scalar_bg else bg)
    otfs, dev = _otfs("simplified", (p,), data.shape, use, use_gpu, be, layout, t)
    if dev is not None:
        r = _gpu_path("decon_lucy", use_gpu, call)
        if r is not None:
            return r
        otfs, _ = _otfs("simplified", (p,), data.shape, use, False, be, layout, t)
    h = otfs[0]

    x = _single_input(data, t)
    mask = _mask_of(x, decon_offset, edge_erosion, t)
    im = _subtract_background(x, scalar_bg, bg, t, inplace=x is not data)
    del x
    n_iter = _iterations(numit)
    dbg = _Debug(n_iter, save_step, debug_folder, t) if debug else None

    if use:
        k = K.LucyRFFT(data.shape, backend=be, threads=t, layout=layout)
        ml = k.matlab_layout
        i_int = k.to_internal(im)
        if ml:
            del im
        j2 = i_int
        j3 = _internal_zeros(k, data.size)
        j4 = _internal_zeros(k, data.size)
        y = None
        lam = 0.0
        for it in range(1, n_iter + 1):
            if it > 2:
                lam = _lambda_step(j2, y, j4, ml, t)
            elif it == 2:
                # J_4 = J_2(:) - Y(:): the lambda step writes it (its factor is not used here)
                _lambda_step(j2, y, j4, ml, t)
            y = K.accel_predict(j2, j3, lam, threads=t, out=y)
            j3 = j2
            j2 = k.rlstep(y, i_int, h, internal=True)
            if dbg is not None:
                dbg.bookkeeping(it, _mse(j2, i_int, ml, t))
                dbg.files(it, lambda j=j2: k.from_internal(j))
        del j3, j4, y
        j = k.from_internal(j2)
        del j2

        def damp_input():
            return k.from_internal(i_int)

    else:
        j2 = im
        j3 = np.zeros(data.shape, dtype=im.dtype)
        j4 = np.zeros(data.shape, dtype=im.dtype)
        y = None
        lam = 0.0
        for it in range(1, n_iter + 1):
            if it > 2:
                lam, j4 = _lambda_any(j2, y, j4, t)
            elif it == 2:
                from petakit5d.core import dtypes

                j4 = np.asarray(dtypes.minus(j2, y))
            y = _predict_any(j2, j3, lam, t)
            r = _clamp0(_conv(y, h, be, t), t, _DBL_EPS)
            r = _divide(_like(im, j2), r)
            j3 = j2
            r = _conv(r, h, be, t, conj=True)
            j2 = _clamp0(_times(y, r), t)
            if dbg is not None:
                dbg.bookkeeping(it, _mse(j2, im, False, t))
                dbg.files(it, lambda j=j2: j)
        del j3, j4, y
        j = j2

        def damp_input():
            return im

    err_mat = dbg.rows_to_istp() if dbg is not None else np.zeros((n_iter, 4))
    overlap_out = overlap_fcn() if overlap_fcn is not None else None
    out = _post(
        j,
        damp_input,
        mask,
        damp_factor=damp_factor,
        scale_factor=scale_factor,
        decon_offset=decon_offset,
        edge_erosion=edge_erosion,
        save_16bit=save_16bit,
        bbox=bbox,
        shape_in=shape_in,
        threads=t,
    )
    return RLResult(out, err_mat, n_iter, overlap_out)


def _divide(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """MATLAB's ``a ./ b`` (single division for single arrays)"""
    if a.dtype == b.dtype:
        with np.errstate(all="ignore"):
            return a / b
    from petakit5d.core import dtypes

    return np.asarray(dtypes.rdivide(a, b))


def _times(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """MATLAB's ``a .* b``"""
    if a.dtype == b.dtype:
        with np.errstate(all="ignore"):
            return a * b
    from petakit5d.core import dtypes

    return np.asarray(dtypes.times(a, b))


# ---------------------------------------------------------------------------------------------
# decon_lucy_omw_function
# ---------------------------------------------------------------------------------------------


def lucy_omw(
    i: np.ndarray,
    psf_f: np.ndarray,
    psf_b: np.ndarray,
    numit: Any,
    *,
    background: Any = None,
    use_gpu: bool | None = False,
    save_16bit: bool = True,
    damp_factor: Any = 1,
    scale_factor: Any = 1,
    decon_offset: Any = 0,
    edge_erosion: Any = 0,
    bbox: Sequence[Any] | None = None,
    debug: bool = False,
    debug_folder: str = "./debug/",
    save_step: Any = 1,
    use_mex: bool | None = None,
    overlap_fcn: Callable[[], Any] | None = None,
    threads: int | None = None,
    backend: str | None = None,
    layout: str | None = None,
) -> OMWResult:
    """``decon_lucy_omw`` with the thread count, the CPU FFT backend and the kernels' layout as
    arguments (as ``lucy``); ``use_gpu`` False by default."""
    call = _gpu_arguments(locals())
    from petakit5d import fft
    from petakit5d.kernels import decon as K

    t = _threads(threads)
    be = fft.backend() if backend is None else backend
    shape_in = np.shape(i)
    data = _data(i)
    scalar_bg, bg = _background(background)
    pf, pb = _psf(psf_f), _psf(psf_b)
    use, _ = use_rfft_kernels(use_mex, "omw", data, debug, None if scalar_bg else bg)
    otfs, dev = _otfs("omw", (pf, pb), data.shape, use, use_gpu, be, layout, t)
    if dev is not None:
        r = _gpu_path("decon_lucy_omw", use_gpu, call)
        if r is not None:
            return r
        otfs, _ = _otfs("omw", (pf, pb), data.shape, use, False, be, layout, t)
    hf, hb = otfs

    if use:
        k = K.LucyRFFT(data.shape, backend=be, threads=t, layout=layout)
        raw = data if data.dtype in (np.uint16, np.float32) else _single_input(data, t)
        mask = _mask_of(raw, decon_offset, edge_erosion, t)
        bgv = float(bg)
        if debug:
            j, err_mat = _omw_kernel_debug(k, raw, hf, hb, numit, bg, save_step, debug_folder, t)
        else:
            j = k.omw(raw, hf, hb, numit, bgv)
            err_mat = np.zeros((_iterations(numit), 4))

        def damp_input():
            return k.omw(raw, hf, hb, 0, bgv)

    else:
        x = _single_input(data, t)
        mask = _mask_of(x, decon_offset, edge_erosion, t)
        im = _subtract_background(x, scalar_bg, bg, t, inplace=x is not data)
        del x
        n_iter = _iterations(numit)
        dbg = _Debug(n_iter, save_step, debug_folder, t) if debug else None
        j = im
        for it in range(1, n_iter + 1):
            cx = _clamp0(_conv(j, hf, be, t), t, _DBL_EPS)
            cx = _conv(_divide(_like(im, j), cx), hb, be, t)
            j = _clamp0(_times(cx, j), t)
            if dbg is not None:
                dbg.bookkeeping(it, _mse(j, im, False, t))
                dbg.files(it, lambda v=j: v)
        err_mat = dbg.rows_kept(n_iter) if dbg is not None else np.zeros((n_iter, 4))

        def damp_input():
            return im

    overlap_out = overlap_fcn() if overlap_fcn is not None else None
    out = _post(
        j,
        damp_input,
        mask,
        damp_factor=damp_factor,
        scale_factor=scale_factor,
        decon_offset=decon_offset,
        edge_erosion=edge_erosion,
        save_16bit=save_16bit,
        bbox=bbox,
        shape_in=shape_in,
        threads=t,
    )
    return OMWResult(out, err_mat, overlap_out)


def _omw_kernel_debug(k, raw, hf, hb, numit, bg, step, folder, threads):
    """omw_mex_debug: the iterations in the kernels in pieces from the last result to the
    iterations the debug outputs need (the error metric of the rows kept, the saved estimates),
    the rest of err_mat as MATLAB's loop computes it"""
    n_iter = _iterations(numit)
    dbg = _Debug(n_iter, step, folder, threads)
    s = _subtract_background(raw, True, bg, threads)
    bgv = float(bg)
    with np.errstate(all="ignore"):
        n_err = np.floor(np.float64(n_iter) / np.float64(dbg.step))
    j, kd = None, 0
    for it in range(1, n_iter + 1):
        if it <= n_err or _rem(it, dbg.step) == 0 or it == n_iter:
            j = k.omw(raw, hf, hb, it, bgv) if j is None else k.omw(raw, hf, hb, it - kd, bgv, j)
            kd = it
        dbg.bookkeeping(it, _mse(j, s, False, threads) if it <= n_err else None)
        dbg.files(it, lambda v=j: v)
    if j is None:
        j = k.omw(raw, hf, hb, 0, bgv)
    return j, dbg.rows_kept(n_iter)


# ---------------------------------------------------------------------------------------------
# the registered functions
# ---------------------------------------------------------------------------------------------


def _matlab_out_rl(r: RLResult):
    from petakit5d.core.arrays import to_matlab

    k = float(r.iterations) if r.iterations > 0 else np.zeros((0, 0))
    return (to_matlab(r.decon), r.err_mat, k, r.overlap_out)


def _matlab_out_omw(r: OMWResult):
    from petakit5d.core.arrays import to_matlab

    return (to_matlab(r.decon), r.err_mat, r.overlap_out)


def _background_in(v: Any) -> Any:
    """MATLAB's Background: a scalar as it is (a double stays a float), an image array in
    Python's form"""
    if v is None or isinstance(v, (int, float)) or np.ndim(v) < 2:
        return v
    from petakit5d.core.arrays import from_matlab

    return from_matlab(v)


def _background_out(v: Any) -> Any:
    if v is None or isinstance(v, (int, float)) or np.ndim(v) < 2:
        return v
    from petakit5d.core.arrays import to_matlab

    return to_matlab(v)


@register(
    "decon_lucy_function",
    Param("I", "@isnumeric", "Data, a (z, y, x) volume or a (y, x) image.", convert="volume"),
    Param("PSF", "@isnumeric", "PSF (z, y, x), at most the data's size.", convert="volume"),
    Param("NUMIT", "@isnumeric", "Number of iterations."),
    Param(
        "Background",
        "@isnumeric",
        "Background subtracted from the data: a scalar ([]: 0), or one per voxel.",
        convert=(_background_in, _background_out),
    ),
    Param(
        "useGPU",
        "@(x) isempty(x) || islogical(x)",
        "True: on a GPU, an error (decon:noUsableGPU) if none is usable; False: on the CPU; "
        "None: a GPU if one is usable, else the CPU with a warning.",
    ),
    Param("save16bit", "@islogical", "Return uint16 (else single).", convert="same"),
    Param("dampFactor", "@isnumeric", "Damp factor (> 1: results at most dampFactor x data)."),
    Param("scaleFactor", "@isnumeric", "Scale factor of the result."),
    Param("deconOffset", "@isnumeric", "Offset added to the result.", convert="same"),
    Param("EdgeErosion", "@isnumeric", "Erosion of the edges of the mask of nonzero data."),
    Param(
        "bbox",
        "@isnumeric",
        "Box (z0, y0, x0, z1, y1, x1) of the result to return (None: all).",
        convert="bbox",
    ),
    Param("debug", "@islogical", "Error metric and estimates of the iterations."),
    Param("debug_folder", "@ischar", "Folder of the debug estimates."),
    Param("saveStep", "@isnumeric", "Save an estimate every saveStep iterations (debug)."),
    Param(
        "useMex",
        "@(x) isempty(x) || islogical(x)",
        "Iterations in the half-spectrum kernels: None where they can, True required, False "
        "MATLAB's loop on full spectra.",
    ),
    Param(
        "overlapFcn",
        "@(x) isempty(x) || isa(x, 'function_handle')",
        "Function run after the iterations (on the GPU while they run); its output is returned.",
    ),
    task="R03",
    matlab_out=_matlab_out_rl,
)
def decon_lucy(
    i,
    psf,
    numit,
    *,
    background=None,
    use_gpu=None,
    save_16bit=True,
    damp_factor=1,
    scale_factor=1,
    decon_offset=0,
    edge_erosion=0,
    bbox=None,
    debug=False,
    debug_folder="./debug/",
    save_step=1,
    use_mex=None,
    overlap_fcn=None,
) -> RLResult:
    """Richardson-Lucy deconvolution with acceleration (decon_lucy_function): the result, err_mat,
    the number of iterations run and the output of ``overlap_fcn`` (see the module)."""
    return lucy(
        i,
        psf,
        numit,
        background=background,
        use_gpu=use_gpu,
        save_16bit=save_16bit,
        damp_factor=damp_factor,
        scale_factor=scale_factor,
        decon_offset=decon_offset,
        edge_erosion=edge_erosion,
        bbox=bbox,
        debug=debug,
        debug_folder=debug_folder,
        save_step=save_step,
        use_mex=use_mex,
        overlap_fcn=overlap_fcn,
    )


@register(
    "decon_lucy_omw_function",
    Param("I", "@isnumeric", "Data, a (z, y, x) volume or a (y, x) image.", convert="volume"),
    Param("PSF_f", "@isnumeric", "Forward PSF (z, y, x).", convert="volume"),
    Param("PSF_b", "@isnumeric", "Back projector (z, y, x).", convert="volume"),
    Param("NUMIT", "@isnumeric", "Number of iterations."),
    Param(
        "Background",
        "@isnumeric",
        "Background subtracted from the data: a scalar ([]: 0), or one per voxel.",
        convert=(_background_in, _background_out),
    ),
    Param(
        "useGPU",
        "@(x) isempty(x) || islogical(x)",
        "True: on a GPU, an error (decon:noUsableGPU) if none is usable; False: on the CPU; "
        "None: a GPU if one is usable, else the CPU with a warning.",
    ),
    Param("save16bit", "@islogical", "Return uint16 (else single).", convert="same"),
    Param("dampFactor", "@isnumeric", "Damp factor (> 1: results at most dampFactor x data)."),
    Param("scaleFactor", "@isnumeric", "Scale factor of the result."),
    Param("deconOffset", "@isnumeric", "Offset added to the result.", convert="same"),
    Param("EdgeErosion", "@isnumeric", "Erosion of the edges of the mask of nonzero data."),
    Param(
        "bbox",
        "@isnumeric",
        "Box (z0, y0, x0, z1, y1, x1) of the result to return (None: all).",
        convert="bbox",
    ),
    Param("debug", "@islogical", "Error metric and estimates of the iterations."),
    Param("debug_folder", "@ischar", "Folder of the debug estimates."),
    Param("saveStep", "@isnumeric", "Save an estimate every saveStep iterations (debug)."),
    Param(
        "useMex",
        "@(x) isempty(x) || islogical(x)",
        "Iterations in the half-spectrum kernels: None where they can, True required, False "
        "MATLAB's loop on full spectra.",
    ),
    Param(
        "overlapFcn",
        "@(x) isempty(x) || isa(x, 'function_handle')",
        "Function run after the iterations (on the GPU while they run); its output is returned.",
    ),
    task="R03",
    matlab_out=_matlab_out_omw,
)
def decon_lucy_omw(
    i,
    psf_f,
    psf_b,
    numit,
    *,
    background=None,
    use_gpu=None,
    save_16bit=True,
    damp_factor=1,
    scale_factor=1,
    decon_offset=0,
    edge_erosion=0,
    bbox=None,
    debug=False,
    debug_folder="./debug/",
    save_step=1,
    use_mex=None,
    overlap_fcn=None,
) -> OMWResult:
    """Richardson-Lucy deconvolution with the OMW back projector (decon_lucy_omw_function): the
    result, err_mat and the output of ``overlap_fcn`` (see the module)."""
    return lucy_omw(
        i,
        psf_f,
        psf_b,
        numit,
        background=background,
        use_gpu=use_gpu,
        save_16bit=save_16bit,
        damp_factor=damp_factor,
        scale_factor=scale_factor,
        decon_offset=decon_offset,
        edge_erosion=edge_erosion,
        bbox=bbox,
        debug=debug,
        debug_folder=debug_folder,
        save_step=save_step,
        use_mex=use_mex,
        overlap_fcn=overlap_fcn,
    )
