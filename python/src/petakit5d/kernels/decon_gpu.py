"""The deconvolution kernels on the GPU, moved from the CUDA mex files of
``microscopeDataProcessing/decon/mex``: the half-spectrum Richardson-Lucy steps of
decon_lucy_rfft_gpu_mex (``LucyRFFTGPU``: 'otf', 'omw', 'rlstep') and the mask erosion of
decon_mask_edge_erosion_gpu_mex (``mask_edge_erosion``, and ``decon_mask_edge_erosion_gpu`` for
decon_mask_edge_erosion.m's gpuArray masks).

The CUDA sources (``cpp/cuda/decon/rl_rfft.cu``, ``edge_erosion.cu``) are compiled at run time
(``petakit5d.gpu.kernels``); the FFTs are single precision real-to-complex cuFFT plans with one
work area (``petakit5d.gpu.fft.RealFFTPlan``). CuPy is imported only when a GPU is used; the shapes,
plan sizes and argument checks need neither CuPy nor a GPU.

**Layouts** (``layout=``):

- ``"matlab"`` (the default): the mex's plans on MATLAB's layout, with the mex's bits on the same
  GPU model with cuFFT 11.0.8 (the version MATLAB R2024b bundles). MATLAB's ``(y, x, z)`` Fortran
  volume is the C array ``(z, x, y)``: the plans transform that, halved along y (cuFFT's
  row-major ``(nz, nx, ny)``, size-1 axes other than the last left out, as the mex's
  ``fft_dims``). Half OTFs are MATLAB's ``H``: Fortran arrays ``(ny // 2 + 1, nx, nz)`` (2D
  ``(ny // 2 + 1, nx)``). Volumes in and out are ``(z, y, x)``; kernels that read or write them
  next to MATLAB-layout arrays go through 32 x 32 tiles in shared memory, so no volume is turned
  around: the input stays as it is (uint16 as uint16), and only RL's ``Y`` is transposed into the
  scratch the mex has anyway. ``internal=True`` takes and returns MATLAB-layout arrays instead
  (``internal_shape``: Fortran ``(ny, nx, nz)``; ``to_internal``/``from_internal``), as RL's loop
  keeps them across iterations.
- ``"native"``: plans of the ``(z, y, x)`` arrays themselves, halved along x (half OTFs
  ``(nz, ny, nx // 2 + 1)``, C order), faster; the element-wise steps are the same, the FFTs round
  differently (measured against the mex in the K09 report, per output type).

**Arithmetic** (the mex's): the background subtracted as MATLAB's GPU code does (uint16: the
saturating ``I - uint16(background)``, with the mex's rounding ``floor(b + 0.5)``; single:
``max(single(double(I) - background), 0)``), MATLAB's complex product without fused multiply-adds,
the ratio ``I ./ max(r / N, eps)`` with double ``eps`` as single and ``1 / N`` as
``single(1 / double(N))``, the update ``max(r / N .* J, 0)`` (-0 to +0). Every kernel is
element-wise (no atomics, no reductions) and cuFFT's plans are deterministic, so results don't
change from run to run or with the grid.

**Memory** as the mex's: an 'omw' call holds the input (uint16 stays uint16), the two half OTFs,
J, one real scratch, one half spectrum and the plans' shared work area (about another half
spectrum); 'rlstep' Y, I, H, its result, a half spectrum and the work area. Where a result goes
to a ``(z, y, x)`` volume, it is allocated after the half spectrum is freed. With
``keep_plans=True`` (the default) the plans come from ``petakit5d.gpu.fft``'s plan cache, whose
work area stays allocated between calls (``petakit5d.gpu.memory.release_cached(plans=True)``
frees it); ``keep_plans=False`` makes and destroys them in each call, as the mex does. Out of GPU
memory raises ``GPUOutOfMemoryError`` (a ``MemoryError``) with MATLAB's identifiers, as the mex:
``parallel:gpu:array:OOM`` (arrays and the work area), ``parallel:gpu:array:OOMForOperation`` and
``parallel:gpu:array:FFTPlanInternalError`` (cuFFT planning); a failed call keeps no arrays.

Arrays may be NumPy arrays (copied to the GPU) or CuPy arrays on the current GPU; results are CuPy
arrays on the current stream. Errors the mex raises carry its identifier
(``DeconKernelError.identifier``).
"""

from __future__ import annotations

import math
import os
import threading
import traceback
from collections.abc import Callable, Sequence
from typing import Any, NamedTuple

import numpy as np

from petakit5d.kernels.decon import DeconKernelError

__all__ = [
    "LAYOUTS",
    "GPUOutOfMemoryError",
    "Geometry",
    "LucyRFFTGPU",
    "background_uint16",
    "decon_mask_edge_erosion_gpu",
    "geometry",
    "mask_edge_erosion",
    "nhood_offsets",
    "omw",
    "otf",
    "rlstep",
]

LAYOUTS = ("matlab", "native")

_RFFT = "decon_lucy_rfft_gpu_mex:inputError"
_EROSION = "decon_mask_edge_erosion_gpu_mex:inputError"
_OOM = "parallel:gpu:array:OOM"
_GPU_BAD_INDEX = "parallel:array:SubsasgnInvalidSubscript"  # gpuArray assignments' index error
_OOM_PLAN = "parallel:gpu:array:OOMForOperation"
_PLAN_FAILED = "parallel:gpu:array:FFTPlanInternalError"
_MAX_DIM = 2147483647
_MAX_OFFSETS = 4096
_THREADS = 256
_TILE, _TILE_ROWS = 32, 8
_MAX_GRID_YZ = 65535

_U16, _F32 = "unsigned short", "float"
_RL_NAMES = tuple(
    f"{k}<{t}>"
    for k in ("init_kernel", "ratio_kernel", "zyx_to_int", "int_to_zyx", "init_t", "ratio_t")
    for t in (_U16, _F32)
)


class GPUOutOfMemoryError(MemoryError):
    """Out of GPU memory, with the identifier MATLAB's GPU code (and the mex) gives."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _err(message: str, identifier: str = _RFFT) -> DeconKernelError:
    return DeconKernelError(identifier, message)


# ---------------------------------------------------------------------------------------------
# geometry: shapes, plans and layouts (no GPU needed)
# ---------------------------------------------------------------------------------------------


class Geometry(NamedTuple):
    """The arrays of the Richardson-Lucy steps for data of one shape in one layout.

    shape: the data, (nz, ny, nx) or (ny, nx). fft_shape: the real C arrays the plans transform
    (MATLAB's layout: (nz, nx, ny) or (nx, ny); native: shape). half_shape: their half spectra,
    C order. plan_dims: cuFFT's dimensions (the mex's fft_dims in MATLAB's layout).
    internal_shape, spectrum_shape: an internal array and a half OTF as callers get them
    (MATLAB's layout: Fortran arrays of MATLAB's sizes; native: C arrays). zyx: (nz, ny, nx).
    """

    shape: tuple[int, ...]
    layout: str
    fft_shape: tuple[int, ...]
    half_shape: tuple[int, ...]
    plan_dims: tuple[int, ...]
    internal_shape: tuple[int, ...]
    spectrum_shape: tuple[int, ...]
    zyx: tuple[int, int, int]

    @property
    def n(self) -> int:
        """voxels"""
        return math.prod(self.zyx)

    @property
    def nh(self) -> int:
        """elements of a half spectrum"""
        return math.prod(self.half_shape)

    @property
    def matlab_layout(self) -> bool:
        return self.layout == "matlab"

    def psf_strides(self) -> tuple[int, int, int]:
        """the strides (z, y, x) of a real array of fft_shape, in elements"""
        nz, ny, nx = self.zyx
        if self.matlab_layout:
            return ny * nx, 1, ny
        return ny * nx, nx, 1


def geometry(shape: Sequence[int], layout: str = "matlab") -> Geometry:
    """The Geometry of data of ``shape`` ((nz, ny, nx) or (ny, nx)) in ``layout``; the mex's
    error (identifier decon_lucy_rfft_gpu_mex:inputError) for other sizes."""
    from petakit5d.gpu.fft import matlab_shape, plan_dims

    if layout not in LAYOUTS:
        raise ValueError(f"layout must be one of {LAYOUTS}, not {layout!r}")
    s = tuple(int(v) for v in shape)
    if len(s) not in (2, 3) or any(v < 1 or v > _MAX_DIM for v in s):
        raise _err("Only 2D and 3D data with dimensions below 2^31 are supported.")
    nz, ny, nx = (1, *s) if len(s) == 2 else s
    if -(-ny // _TILE) > _MAX_GRID_YZ:
        raise ValueError(f"at most {_TILE * _MAX_GRID_YZ} voxels along y, not {ny}")
    if layout == "matlab":
        fs = matlab_shape(s)
        internal = (ny, nx, nz) if len(s) == 3 else (ny, nx)
        spectrum = (ny // 2 + 1, nx, nz) if len(s) == 3 else (ny // 2 + 1, nx)
    else:
        fs = s
        internal = s
        spectrum = (*s[:-1], nx // 2 + 1)
    half = (*fs[:-1], fs[-1] // 2 + 1)
    return Geometry(s, layout, fs, half, plan_dims(fs), internal, spectrum, (nz, ny, nx))


def _grid(n: int) -> int:
    """the mex's grid of 256-thread blocks for n elements (grid-stride loops)"""
    b = -(-n // _THREADS)
    return max(1, min(b, 65535 * 32))


def _tile_grid(g: Geometry) -> tuple[int, int, int]:
    nz, ny, nx = g.zyx
    return (-(-nx // _TILE), -(-ny // _TILE), max(1, min(nz, _MAX_GRID_YZ)))


# ---------------------------------------------------------------------------------------------
# arguments (the mex's checks, in its order; no GPU needed)
# ---------------------------------------------------------------------------------------------


def _is_cupy(a: Any) -> bool:
    t = type(a)
    return t.__name__ == "ndarray" and t.__module__.split(".")[0] == "cupy"


def _is_array(a: Any) -> bool:
    return isinstance(a, np.ndarray) or _is_cupy(a)


def _count(v: Any, name: str) -> int:
    """a count as a MATLAB loop 1 : v runs it (0 for v < 1 or NaN), the mex's scalar_int"""
    if _is_cupy(v):
        raise _err(f"{name} must be a real numeric scalar.")
    a = np.asarray(v)
    if a.size != 1 or a.dtype.kind not in "biuf":
        raise _err(f"{name} must be a real numeric scalar.")
    x = float(a.reshape(-1)[0])
    if not x >= 1:
        return 0
    if x > 2147483647.0:
        raise _err(f"{name} is too large.")
    return int(x)


def _background(v: Any) -> float:
    if _is_cupy(v):
        raise _err("background must be a real numeric scalar.")
    a = np.asarray(v)
    if a.size != 1 or a.dtype.kind not in "biuf":
        raise _err("background must be a real numeric scalar.")
    return float(a.reshape(-1)[0])


def background_uint16(bg: float) -> int:
    """the mex's uint16(background): rounded half away from zero as floor(bg + 0.5) in double
    (so 0.49999999999999994 gives 1, unlike MATLAB's uint16), saturated, NaN 0"""
    if bg != bg or bg <= 0:
        return 0
    if bg >= 65535:
        return 65535
    return int(math.floor(bg + 0.5))


def _real_class(a: Any) -> str | None:
    """'uint16' or 'single' for a real non-empty NumPy or CuPy array of that class, else None"""
    if not _is_array(a) or a.size == 0:
        return None
    return {"uint16": "uint16", "float32": "single"}.get(a.dtype.name)


# ---------------------------------------------------------------------------------------------
# GPU plumbing
# ---------------------------------------------------------------------------------------------

_modules: dict[tuple[str, int, str], Any] = {}
_modules_lock = threading.Lock()


def _cp() -> Any:
    from petakit5d.gpu._cupy import import_cupy

    return import_cupy()


def _module(filename: str) -> Any:
    """the compiled kernels of cpp/cuda/decon/<filename> on the current GPU"""
    cp = _cp()
    key = (filename, cp.cuda.Device().id, os.environ.get("PETAKIT5D_GPU_TARGET", ""))
    with _modules_lock:
        m = _modules.get(key)
        if m is None:
            from petakit5d.gpu import kernels

            exprs = _RL_NAMES if filename == "rl_rfft.cu" else ()
            m = kernels.load_module(kernels.cuda_source("decon", filename), name_expressions=exprs)
            _modules[key] = m
        return m


def _run(fn: Callable[..., Any], *args: Any) -> Any:
    """fn(*args), with CuPy's and cuFFT's out-of-memory errors as MATLAB's (the mex's mapping);
    the failed call's arrays are released before the error reaches the caller"""
    cp = _cp()
    from petakit5d.gpu.fft import CuFFTMemoryError

    err: GPUOutOfMemoryError | None = None
    try:
        return fn(*args)
    except cp.cuda.memory.OutOfMemoryError as e:
        err = GPUOutOfMemoryError(_OOM, f"Out of memory on device. {e}")
        traceback.clear_frames(e.__traceback__)
    except CuFFTMemoryError as e:
        if e.code == 5:
            err = GPUOutOfMemoryError(_PLAN_FAILED, f"cuFFT failed to make a plan ({e}).")
        else:
            err = GPUOutOfMemoryError(
                _OOM_PLAN, f"Out of memory on device for the cuFFT plans ({e})."
            )
        traceback.clear_frames(e.__traceback__)
    raise err


def _upload(a: Any, cp: Any) -> Any:
    """a C-contiguous CuPy array of a (a copy of a NumPy array or of a strided CuPy array)"""
    if _is_cupy(a):
        if a.device.id != cp.cuda.Device().id:
            raise ValueError(f"the array is on GPU {a.device.id}, not the current GPU")
        return a if a.flags.c_contiguous else cp.ascontiguousarray(a)
    return cp.asarray(np.ascontiguousarray(a))


def _same_device(cp: Any, *arrays: Any) -> None:
    dev = cp.cuda.Device().id
    for a in arrays:
        if a.device.id != dev:
            raise ValueError(f"an OTF is on GPU {a.device.id}, not the current GPU {dev}")


class _Timer:
    """the mex's t: seconds in [planning, r2c FFTs, c2r FFTs, element-wise kernels], each step
    synchronized (only when timing)"""

    def __init__(self, on: bool, cp: Any):
        self.on = on
        self.t = [0.0, 0.0, 0.0, 0.0]
        if on:
            self._cp = cp
            self._e0, self._e1 = cp.cuda.Event(), cp.cuda.Event()

    def start(self) -> None:
        if self.on:
            self._e0.record()

    def stop(self, k: int) -> None:
        if self.on:
            self._e1.record()
            self._e1.synchronize()
            self.t[k] += self._cp.cuda.get_elapsed_time(self._e0, self._e1) / 1000.0


# ---------------------------------------------------------------------------------------------
# Richardson-Lucy steps (decon_lucy_rfft_gpu_mex)
# ---------------------------------------------------------------------------------------------


class LucyRFFTGPU:
    """The half-spectrum Richardson-Lucy steps of decon_lucy_rfft_gpu_mex for data of one shape.

    shape: (nz, ny, nx) or (ny, nx). layout: 'matlab' (the mex's plans and bits) or 'native' (see
    the module). keep_plans: take the plans from the plan cache (True) or make them in each call
    (False, the mex's memory between calls). timing: fill ``last_times`` (seconds in planning,
    r2c FFTs, c2r FFTs and element-wise kernels of the last call; every step is synchronized).
    The methods take the mex's arguments in the port's form and run on the current GPU.
    """

    def __init__(
        self,
        shape: Sequence[int],
        *,
        layout: str = "matlab",
        keep_plans: bool = True,
        timing: bool = False,
    ):
        self.geometry = geometry(shape, layout)
        g = self.geometry
        self.shape = g.shape
        self.layout = g.layout
        self.internal_shape = g.internal_shape
        self.spectrum_shape = g.spectrum_shape
        self.keep_plans = bool(keep_plans)
        self.timing = bool(timing)
        self.last_times = (0.0, 0.0, 0.0, 0.0)
        self._inv_n = np.float32(1.0 / float(g.n))

    @property
    def matlab_layout(self) -> bool:
        """whether the internal arrays and half OTFs are in MATLAB's layout"""
        return self.geometry.matlab_layout

    # -- arrays ------------------------------------------------------------------------------------

    def _expose(self, c: Any) -> Any:
        """a C array of fft_shape (or half_shape) as callers get it: MATLAB's Fortran array"""
        return c.T if self.matlab_layout else c

    def _check(self, a: Any, name: str, internal: bool) -> None:
        """a (z, y, x) volume of the data's shape, or an internal array (MATLAB's layout: of its
        Fortran shape and order)"""
        g = self.geometry
        if not internal:
            if tuple(a.shape) != self.shape:
                raise _err(f"{name} must have the data's shape {self.shape}, not {tuple(a.shape)}.")
            return
        if tuple(a.shape) != g.internal_shape:
            raise _err(f"{name} must be an internal array of shape {g.internal_shape}.")
        if self.matlab_layout and not a.flags.f_contiguous:
            raise ValueError(f"{name}: an internal array in MATLAB's layout is Fortran-contiguous")

    def _c(self, a: Any, internal: bool, cp: Any) -> Any:
        """the C array on the GPU of a checked volume, or of fft_shape of an internal array"""
        return _upload(a.T if internal and self.matlab_layout else a, cp)

    def _spectrum(self, h: Any) -> Any:
        """the C array of half_shape of a half OTF (the mex's check_otf)"""
        g = self.geometry
        if not _is_cupy(h) or h.dtype != np.complex64:
            raise _err("The OTFs must be complex single half OTF gpuArrays.")
        if tuple(h.shape) != g.spectrum_shape:
            raise _err(
                "The OTF size must be [floor(ny / 2) + 1, nx, nz] of the data."
                if self.matlab_layout
                else f"The OTF size must be the native half spectrum {g.spectrum_shape}."
            )
        c = h.T if self.matlab_layout else h
        if not c.flags.c_contiguous:
            raise ValueError(
                "a half OTF in MATLAB's layout is Fortran-contiguous"
                if self.matlab_layout
                else "a native half OTF is C-contiguous"
            )
        return c

    def to_internal(self, vol: Any) -> Any:
        """A (z, y, x) volume (float32 or uint16, NumPy or CuPy) as this object's internal array on
        the GPU: MATLAB's Fortran array (by a tiled transpose), or the volume itself (native)."""
        self._check(vol, "The volume", False)
        if vol.dtype.name not in ("float32", "uint16"):
            raise ValueError("to_internal: a float32 or uint16 volume")
        a = self._c(vol, False, _cp())
        if not self.matlab_layout:
            return a
        return _run(self._transpose, a, "zyx_to_int").T

    def from_internal(self, a: Any) -> Any:
        """The (z, y, x) volume of an internal array, on the GPU."""
        self._check(a, "The internal array", True)
        if a.dtype.name not in ("float32", "uint16"):
            raise ValueError("from_internal: a float32 or uint16 array")
        c = self._c(a, True, _cp())
        if not self.matlab_layout:
            return c
        return _run(self._transpose, c, "int_to_zyx")

    def _transpose(self, src: Any, kernel: str) -> Any:
        cp = _cp()
        g = self.geometry
        nz, ny, nx = g.zyx
        dst = cp.empty(g.fft_shape if kernel == "zyx_to_int" else g.shape, dtype=src.dtype)
        t = _U16 if src.dtype == np.uint16 else _F32
        _module("rl_rfft.cu")[f"{kernel}<{t}>"](
            _tile_grid(g), (_TILE, _TILE_ROWS), (dst, src, np.int32(ny), np.int32(nx), np.int32(nz))
        )
        return dst

    def _plan(self, *, inverse: bool = True) -> Any:
        """the R2C and C2R plans: cached (with the inverse, which the iterations share), or made
        for the call (the OTF's without the inverse, as the mex makes them)"""
        from petakit5d.gpu import fft

        if self.keep_plans:
            return fft.get_plan(self.geometry.fft_shape)
        return fft.RealFFTPlan(self.geometry.fft_shape, inverse=inverse)

    def _done(self, plan: Any, timer: _Timer) -> None:
        if not self.keep_plans:
            plan.close()
        self.last_times = tuple(timer.t)

    # -- 'otf' -------------------------------------------------------------------------------------

    def otf(self, psf: Any) -> Any:
        """``decon_lucy_rfft_gpu_mex('otf', psf, outSize)``: the half OTF of a normalized single
        PSF ((z, y, x) or (y, x), at most the data's size along each axis): decon_psf2otf's OTF
        (the PSF zero-padded and shifted so that its center floor(size / 2) is at the origin)
        halved, MATLAB's H (``spectrum_shape``; Fortran) or the native half spectrum."""
        if _real_class(psf) != "single":
            raise _err("The PSF must be a real, non-empty single array.")
        if psf.ndim > 3 or any(v > _MAX_DIM for v in psf.shape):
            raise _err("Only 2D and 3D data with dimensions below 2^31 are supported.")
        ps = (1,) * (3 - psf.ndim) + tuple(int(v) for v in psf.shape)
        if any(p > d for p, d in zip(ps, self.geometry.zyx, strict=True)):
            raise _err(
                "outSizeIsSmallerThanPsfSize", "decon_lucy_rfft_gpu_mex:outSizeIsSmallerThanPsfSize"
            )
        return _run(self._otf, psf, ps)

    def _otf(self, psf: Any, ps: tuple[int, int, int]) -> Any:
        cp = _cp()
        g = self.geometry
        p = _upload(psf, cp)
        timer = _Timer(self.timing, cp)
        r = cp.zeros(g.fft_shape, dtype=cp.float32)
        s = cp.empty(g.half_shape, dtype=cp.complex64)
        timer.start()
        plan = self._plan(inverse=False)
        timer.stop(0)
        try:
            nz, ny, nx = g.zyx
            sz, sy, sx = g.psf_strides()
            npsf = math.prod(ps)
            timer.start()
            _module("rl_rfft.cu")["scatter_psf_kernel"](
                _grid(npsf),
                _THREADS,
                (
                    r,
                    p,
                    *(np.int32(v) for v in ps),
                    *(np.int32(v) for v in (nz, ny, nx)),
                    *(np.int64(v) for v in (sz, sy, sx)),
                ),
            )
            timer.stop(3)
            timer.start()
            plan.forward(r, s)
            timer.stop(1)
        finally:
            self._done(plan, timer)
        return self._expose(s)

    # -- 'omw' -------------------------------------------------------------------------------------

    def omw(
        self, raw: Any, h_f: Any, h_b: Any, n_iter: Any, background: Any, *, internal: bool = False
    ) -> Any:
        """``decon_lucy_rfft_gpu_mex('omw', I, H_f, H_b, nIter, background)``: ``n_iter`` OMW
        iterations ``J <- max(real(ifftn(fftn(I ./ max(real(ifftn(fftn(J) .* H_f)), eps)) .* H_b))
        .* J, 0)`` from J = I, I the background-subtracted raw data (uint16 or single, (z, y, x),
        or an internal array with ``internal=True``); the half OTFs from ``otf`` on the GPU. The
        result: single, (z, y, x) or internal."""
        n = _count(n_iter, "nIter")
        bg = _background(background)
        cls = _real_class(raw)
        if cls is None or raw.dtype.kind == "c":
            raise _err("The data must be a real, non-empty uint16 or single array.")
        self._check(raw, "The data", internal)
        hf, hb = self._spectrum(h_f), self._spectrum(h_b)
        bgv = float(background_uint16(bg)) if cls == "uint16" else bg
        cp = _cp()
        a = self._c(raw, internal, cp)
        _same_device(cp, hf, hb)
        out = _run(self._omw, a, hf, hb, n, bgv, internal)
        return self._expose(out) if internal else out

    def _omw(self, raw: Any, hf: Any, hb: Any, n_iter: int, bg: float, internal: bool) -> Any:
        cp = _cp()
        g = self.geometry
        mod = _module("rl_rfft.cu")
        t = _U16 if raw.dtype == np.uint16 else _F32
        nz, ny, nx = g.zyx
        dims = (np.int32(ny), np.int32(nx), np.int32(nz))
        linear = internal or not self.matlab_layout
        n, nh, inv_n = np.int64(g.n), np.int64(g.nh), self._inv_n
        timer = _Timer(self.timing, cp)
        bgd = np.float64(bg)
        if n_iter == 0:
            # the background-subtracted input in the input's layout
            j = cp.empty(raw.shape, dtype=cp.float32)
            timer.start()
            mod[f"init_kernel<{t}>"](_grid(g.n), _THREADS, (j, raw, bgd, n))
            timer.stop(3)
            self.last_times = tuple(timer.t)
            return j
        j = cp.empty(g.fft_shape, dtype=cp.float32)
        timer.start()
        if linear:
            mod[f"init_kernel<{t}>"](_grid(g.n), _THREADS, (j, raw, bgd, n))
        else:
            mod[f"init_t<{t}>"](_tile_grid(g), (_TILE, _TILE_ROWS), (j, raw, bgd, *dims))
        timer.stop(3)
        r = cp.empty(g.fft_shape, dtype=cp.float32)
        s = cp.empty(g.half_shape, dtype=cp.complex64)
        timer.start()
        plan = self._plan()
        timer.stop(0)
        try:
            for k in range(1, n_iter + 1):
                for h, last in ((hf, False), (hb, True)):
                    timer.start()
                    plan.forward(j if not last else r, s)
                    timer.stop(1)
                    timer.start()
                    mod["multiply_kernel"](_grid(g.nh), _THREADS, (s, h, nh, np.int32(0)))
                    timer.stop(3)
                    timer.start()
                    plan.inverse(s, r)
                    timer.stop(2)
                    if last:
                        break
                    timer.start()
                    if linear:
                        mod[f"ratio_kernel<{t}>"](_grid(g.n), _THREADS, (r, raw, bgd, inv_n, n))
                    else:
                        mod[f"ratio_t<{t}>"](
                            _tile_grid(g), (_TILE, _TILE_ROWS), (r, raw, bgd, inv_n, *dims)
                        )
                    timer.stop(3)
                timer.start()
                if k == n_iter and not linear:
                    # the result as a (z, y, x) volume, in the half spectrum's memory
                    del s
                    out = cp.empty(g.shape, dtype=cp.float32)
                    mod["update_out_t"](
                        _tile_grid(g), (_TILE, _TILE_ROWS), (out, r, j, inv_n, *dims)
                    )
                    timer.stop(3)
                    return out
                mod["update_kernel"](_grid(g.n), _THREADS, (j, r, j, inv_n, n))
                timer.stop(3)
            return j
        finally:
            self._done(plan, timer)

    # -- 'rlstep' ----------------------------------------------------------------------------------

    def rlstep(self, y: Any, i: Any, h: Any, *, internal: bool = False) -> Any:
        """``decon_lucy_rfft_gpu_mex('rlstep', Y, I, H)``: one update of decon_lucy_function,
        ``max(Y .* real(ifftn(conj(H) .* fftn(I ./ max(real(ifftn(H .* fftn(Y))), eps)))), 0)``
        (the ratio's numerator ``max(I, 0)``, as in the mex). Y: single; I: background-subtracted
        data, uint16 or single; both (z, y, x), or internal arrays with ``internal=True``, which
        the result then is too."""
        if _real_class(y) != "single" or _real_class(i) is None:
            raise _err("Y must be a real single and I a real uint16 or single array.")
        if tuple(y.shape) != tuple(i.shape):
            raise _err("Y and I must have the same size.")
        self._check(y, "Y", internal)
        self._check(i, "I", internal)
        hc = self._spectrum(h)
        cp = _cp()
        yc, ic = self._c(y, internal, cp), self._c(i, internal, cp)
        _same_device(cp, hc)
        out = _run(self._rlstep, yc, ic, hc, internal)
        return self._expose(out) if internal else out

    def _rlstep(self, y: Any, raw: Any, h: Any, internal: bool) -> Any:
        cp = _cp()
        g = self.geometry
        mod = _module("rl_rfft.cu")
        t = _U16 if raw.dtype == np.uint16 else _F32
        nz, ny, nx = g.zyx
        dims = (np.int32(ny), np.int32(nx), np.int32(nz))
        linear = internal or not self.matlab_layout
        n, nh, inv_n = np.int64(g.n), np.int64(g.nh), self._inv_n
        zero = np.float64(0)
        timer = _Timer(self.timing, cp)
        r = cp.empty(g.fft_shape, dtype=cp.float32)
        s = cp.empty(g.half_shape, dtype=cp.complex64)
        timer.start()
        plan = self._plan()
        timer.stop(0)
        try:
            if linear:
                src = y  # an out-of-place R2C leaves its input as it is
            else:
                timer.start()
                mod["zyx_to_int<float>"](_tile_grid(g), (_TILE, _TILE_ROWS), (r, y, *dims))
                timer.stop(3)
                src = r
            for conj in (0, 1):
                timer.start()
                plan.forward(src, s)
                timer.stop(1)
                timer.start()
                mod["multiply_kernel"](_grid(g.nh), _THREADS, (s, h, nh, np.int32(conj)))
                timer.stop(3)
                timer.start()
                plan.inverse(s, r)
                timer.stop(2)
                if conj:
                    break
                timer.start()
                if linear:
                    mod[f"ratio_kernel<{t}>"](_grid(g.n), _THREADS, (r, raw, zero, inv_n, n))
                else:
                    mod[f"ratio_t<{t}>"](
                        _tile_grid(g), (_TILE, _TILE_ROWS), (r, raw, zero, inv_n, *dims)
                    )
                timer.stop(3)
                src = r
            timer.start()
            if linear:
                mod["update_kernel"](_grid(g.n), _THREADS, (r, r, y, inv_n, n))
                timer.stop(3)
                return r
            # the result as a (z, y, x) volume, in the half spectrum's memory
            del s
            out = cp.empty(g.shape, dtype=cp.float32)
            mod["update_y_t"](_tile_grid(g), (_TILE, _TILE_ROWS), (out, r, y, inv_n, *dims))
            timer.stop(3)
            return out
        finally:
            self._done(plan, timer)


def otf(psf: Any, shape: Sequence[int], **kw: Any) -> Any:
    """``decon_lucy_rfft_gpu_mex('otf', psf, outSize)`` for data of ``shape`` (``LucyRFFTGPU``
    keywords)."""
    return LucyRFFTGPU(shape, **kw).otf(psf)


def omw(raw: Any, h_f: Any, h_b: Any, n_iter: Any, background: Any, **kw: Any) -> Any:
    """``decon_lucy_rfft_gpu_mex('omw', I, H_f, H_b, nIter, background)`` (``LucyRFFTGPU.omw``)."""
    if not _is_array(raw):
        raise _err("The data must be a real, non-empty uint16 or single array.")
    return LucyRFFTGPU(raw.shape, **kw).omw(raw, h_f, h_b, n_iter, background)


def rlstep(y: Any, i: Any, h: Any, **kw: Any) -> Any:
    """``decon_lucy_rfft_gpu_mex('rlstep', Y, I, H)`` (``LucyRFFTGPU.rlstep``)."""
    if not _is_array(y):
        raise _err("Y must be a real single and I a real uint16 or single array.")
    return LucyRFFTGPU(y.shape, **kw).rlstep(y, i, h)


# ---------------------------------------------------------------------------------------------
# mask erosion (decon_mask_edge_erosion_gpu_mex)
# ---------------------------------------------------------------------------------------------


def _scalar_arg(v: Any, name: str) -> float:
    """the mex's scalar_arg: a real numeric or logical scalar on the host"""
    if _is_cupy(v):
        raise _err(f"{name} must be a real numeric scalar.", _EROSION)
    a = np.asarray(v)
    if a.size != 1 or a.dtype.kind not in "biuf":
        raise _err(f"{name} must be a real numeric scalar.", _EROSION)
    return float(a.reshape(-1)[0])


def _is_double(v: Any) -> bool:
    """MATLAB's class of v is double: a float64 array or scalar, a Python number or a list or
    tuple of Python numbers (the oracle's mapping)"""
    if isinstance(v, (list, tuple)):
        return all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in v)
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return True
    return isinstance(v, (np.ndarray, np.floating)) and np.asarray(v).dtype == np.float64


def _pad_arg(v: Any, name: str) -> tuple[int, int]:
    """the mex's pad_arg: a double vector of 2 nonnegative integers"""
    if not _is_double(v) or np.asarray(v).size != 2:
        raise _err(f"{name} must be a double vector of 2 elements.", _EROSION)
    out = []
    for x in np.asarray(v, dtype=np.float64).reshape(-1).tolist():
        if not (x >= 0) or x > 1e9 or x != math.floor(x):
            raise _err(f"{name} must contain nonnegative integers.", _EROSION)
        out.append(int(x))
    return out[0], out[1]


def nhood_offsets(nhoods: Sequence[Any]) -> tuple[list[int], list[int], list[int]]:
    """(start, rows, columns): the offsets (row = y, column = x) of the nonzero elements of each
    neighborhood about its center floor((size + 1) / 2) (1-based), in MATLAB's column-major order,
    and where each neighborhood's begin; numeric neighborhoods are converted as MATLAB's
    ``logical`` (NaN: MATLAB's error)"""
    start, dr, dc = [0], [], []
    for e in nhoods:
        if _is_cupy(e):
            raise _err("Each neighborhood must be a real 2D array.", _EROSION)
        a = np.asarray(e)
        if a.ndim > 2 or a.size == 0 or a.dtype.kind not in "biuf":
            raise _err("Each neighborhood must be a real 2D array.", _EROSION)
        a = a.reshape((1,) * (2 - a.ndim) + a.shape)
        if a.dtype.kind == "f" and np.isnan(a).any():
            raise _err("NaN values cannot be converted to logicals.", "MATLAB:nologicalnan")
        m, k = a.shape
        cy, cx = (m + 1) // 2 - 1, (k + 1) // 2 - 1
        cols, rows = np.nonzero((a != 0).T)  # column-major: by column, then row
        if rows.size == 0:
            raise _err("Each neighborhood must have a nonzero element.", _EROSION)
        dr.extend(int(r - cy) for r in rows.tolist())
        dc.extend(int(c - cx) for c in cols.tolist())
        start.append(len(dr))
    return start, dr, dc


def mask_edge_erosion(
    mask: Any,
    nhoods: Sequence[Any],
    pad_ul: Any,
    pad_lr: Any,
    rz: Any,
    clear_border: Any,
) -> Any:
    """``decon_mask_edge_erosion_gpu_mex(A, nhoods, padUL, padLR, rz, clearBorder)``: each z slice
    of the bool mask ((z, y, x) or (y, x), NumPy or CuPy) padded with true by ``pad_ul`` and
    ``pad_lr`` ((rows, columns)), eroded by the 2D neighborhoods in turn (rows y, columns x; logical
    or numeric, nonzero is true) as imerode applies a decomposed structuring element (pixels
    outside the padded slice count as true) and cropped back, then (3D) eroded along z by a line of
    ``2 * rz + 1`` voxels (voxels outside count as true); ``clear_border`` sets the first and last
    rows and columns (3D: and slices) false first. Returns a bool CuPy array of the mask's
    shape, identical to decon_mask_edge_erosion_mex's and imerode's result."""
    if not isinstance(nhoods, (list, tuple)):
        raise _err("nhoods must be a cell array.", _EROSION)
    pu, pl = _pad_arg(pad_ul, "padUL"), _pad_arg(pad_lr, "padLR")
    rzv = _scalar_arg(rz, "rz")
    if not (rzv >= 0) or rzv > 1e9 or rzv != math.floor(rzv):
        raise _err("rz must be a nonnegative integer.", _EROSION)
    clear = _scalar_arg(clear_border, "clearBorder") != 0
    start, dr, dc = nhood_offsets(nhoods)
    if len(dr) > _MAX_OFFSETS:
        raise _err(f"The neighborhoods have more than {_MAX_OFFSETS} elements.", _EROSION)
    if not _is_array(mask) or mask.dtype != np.bool_ or mask.ndim > 3:
        raise _err("The mask must be a full 2D or 3D logical array.", _EROSION)
    return _run(_erode, mask, start, dr, dc, pu, pl, int(rzv), clear)


def _erode(mask, start, dr, dc, pu, pl, rz: int, clear: bool) -> Any:
    cp = _cp()
    a = _upload(mask, cp)
    d, h, w = (1,) * (3 - a.ndim) + tuple(int(v) for v in a.shape)
    hp, wp = h + pu[0] + pl[0], w + pu[1] + pl[1]
    out = cp.empty(a.shape, dtype=cp.bool_)
    if d * h * w == 0:
        return out
    mod = _module("edge_erosion.cu")
    i64 = np.int64
    off = cp.asarray(
        np.array([dr, dc], dtype=np.int32).T.copy() if dr else np.zeros((1, 2), np.int32)
    )
    p = cp.empty((d, hp, wp), dtype=cp.uint8)
    q = cp.empty((d, hp, wp), dtype=cp.uint8)
    block = (32, 8)
    grid = (min(-(-wp // 32), 65535), min(-(-hp // 8), 65535), min(d, 65535))
    geo = (i64(h), i64(w), i64(d), i64(pu[0]), i64(pu[1]), i64(hp), i64(wp))
    mod["load_kernel"](grid, block, (p, a.view(cp.uint8), *geo, np.int32(clear)))
    for k in range(len(start) - 1):
        mod["erode_kernel"](
            grid,
            block,
            (
                q,
                p,
                off,
                np.int32(start[k]),
                np.int32(start[k + 1] - start[k]),
                i64(d),
                i64(hp),
                i64(wp),
            ),
        )
        p, q = q, p
    zgrid = (min(-(-w // 32), 65535), min(-(-h // 8), 65535), 1)
    mod["z_kernel"](zgrid, block, (out.view(cp.uint8), p, i64(rz if d > 1 else 0), *geo))
    return out


def decon_mask_edge_erosion_gpu(mask: Any, edge_erosion: Any) -> Any:
    """decon_mask_edge_erosion.m for gpuArray masks: the mask (bool or numeric, (z, y, x) or (y,
    x), NumPy or CuPy; numeric counts its positive values) with its edges eroded by
    ``edge_erosion`` (taken as a double), a bool CuPy array. A full mask gets its first and last
    ``edge_erosion`` rows, columns and (3D) slices set false directly (within the mask: a
    dimension shorter than the erosion is all edge); others go through ``mask_edge_erosion``
    with the decomposition of ``strel('disk', edge_erosion - 1)`` and the line along z, for
    erosions the mex takes (an integer of at least 1), else MATLAB's imerode code with its errors
    (computed on the CPU: MATLAB's GPU imerode gives the CPU's result). As
    ``petakit5d.decon.erosion.decon_mask_edge_erosion`` on the CPU, with the same results."""
    from petakit5d.decon import erosion as E

    cp = _cp()
    m = _upload(mask, cp)
    if m.dtype != cp.bool_:
        if m.dtype.kind == "c":
            m = m.real
        if m.dtype.kind not in "biuf":
            raise TypeError(f"the mask must be logical or numeric, not {m.dtype}")
        m = m > 0
    if m.ndim not in (2, 3):
        raise ValueError("decon_mask_edge_erosion_gpu: a (z, y, x) volume or a (y, x) image")
    e = E._scalar(edge_erosion)
    ndims = 3 if m.ndim == 3 and m.shape[0] > 1 else 2
    if bool(m.all()):
        # MATLAB's edge indices and their check (petakit5d.decon.erosion), with a gpuArray's error
        try:
            rows = E._edge_rows(m.shape[-2], e, 1)
            cols = E._edge_rows(m.shape[-1], e, 2)
            planes = E._edge_rows(m.shape[0], e, 3) if ndims == 3 else None
        except E.MatlabError as err:
            if err.identifier != E._BAD_INDEX:
                raise
            raise E.MatlabError(
                _GPU_BAD_INDEX,
                "Subscript indices must either be real positive integers or logicals.",
            ) from None
        out = m.copy()
        out[..., cp.asarray(rows), :] = False
        out[..., cp.asarray(cols)] = False
        if planes is not None:
            out[cp.asarray(planes)] = False
        return out
    if E._mex_takes(e):
        nhoods, pad_ul, pad_lr = E.disk_decomposition(int(e) - 1)
        rz = int(e) - 1 if ndims == 3 else 0
        return mask_edge_erosion(m, nhoods, pad_ul, pad_lr, rz, True)
    return cp.asarray(E._imerode_code(cp.asnumpy(m), e, ndims, None))
