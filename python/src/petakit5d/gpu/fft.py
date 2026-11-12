"""Single precision cuFFT plans with the work area under our control.

RealFFTPlan, as the GPU decon mex (decon_lucy_rfft_gpu_mex): a (y, x, z) Fortran-ordered volume
is cuFFT's row-major (z, x, y), so the half spectrum is along y, (ny // 2 + 1, nx, nz)
Fortran-ordered; size-1 dimensions other than y are left out of the plan; the forward (R2C) and
inverse (C2R) plans share one work area, allocated from CuPy's memory pool, the larger of the two.
ComplexFFTPlan: the same for complex-to-complex transforms, which is what MATLAB's gpuArray fftn
and ifftn use (fftn of real data is the forward transform of complex(x)).

cuFFT is called through ctypes on the library CuPy loads (cuda-pathfinder: the cuFFT wheel),
since CuPy's own plans each allocate a work area. PETAKIT5D_CUFFT_LIBRARY names another
libcufft to use instead (testing, e.g. the one MATLAB bundles).
"""

from __future__ import annotations

import contextlib
import os
import sys
import threading
import weakref
from collections import OrderedDict
from collections.abc import Sequence
from typing import Any

from petakit5d.gpu._cupy import import_cupy

__all__ = [
    "ComplexFFTPlan",
    "CuFFTError",
    "CuFFTMemoryError",
    "RealFFTPlan",
    "clear_plan_cache",
    "cufft_info",
    "get_complex_plan",
    "get_plan",
    "plan_dims",
    "set_plan_cache_size",
]

_R2C, _C2R, _C2C = 0x2A, 0x2C, 0x29
_FORWARD, _INVERSE = -1, 1
_ALLOC_FAILED, _INTERNAL_ERROR, _INVALID_SIZE = 2, 5, 8
_NAMES = {
    1: "CUFFT_INVALID_PLAN",
    2: "CUFFT_ALLOC_FAILED",
    3: "CUFFT_INVALID_TYPE",
    4: "CUFFT_INVALID_VALUE",
    5: "CUFFT_INTERNAL_ERROR",
    6: "CUFFT_EXEC_FAILED",
    7: "CUFFT_SETUP_FAILED",
    8: "CUFFT_INVALID_SIZE",
    9: "CUFFT_UNALIGNED_DATA",
    11: "CUFFT_INVALID_DEVICE",
    13: "CUFFT_NO_WORKSPACE",
    14: "CUFFT_NOT_IMPLEMENTED",
    16: "CUFFT_NOT_SUPPORTED",
}

_lib_lock = threading.Lock()
_lib: Any = None
_lib_path: str | None = None

_plans: OrderedDict[tuple[Any, ...], _FFTPlan] = OrderedDict()
_plans_lock = threading.Lock()
_plan_cache_size = 2


class CuFFTError(RuntimeError):
    """A cuFFT call failed; code is cuFFT's result."""

    def __init__(self, code: int, what: str) -> None:
        super().__init__(f"cuFFT error {code} ({_NAMES.get(code, 'unknown')}) in {what}.")
        self.code = code


class CuFFTMemoryError(CuFFTError, MemoryError):
    """Out of GPU memory for a plan (or a plan that fails to make, as the mex counts it)."""


def _check(code: int, what: str, planning: bool = False) -> None:
    if code == 0:
        return
    if code == _ALLOC_FAILED or (planning and code == _INTERNAL_ERROR):
        raise CuFFTMemoryError(code, what)
    raise CuFFTError(code, what)


def _library() -> Any:
    global _lib, _lib_path
    with _lib_lock:
        if _lib is not None:
            return _lib
        import ctypes
        from ctypes import POINTER, c_int, c_longlong, c_size_t, c_void_p

        path = os.environ.get("PETAKIT5D_CUFFT_LIBRARY") or None
        if path is None:
            import_cupy()
            try:
                from cuda import pathfinder

                path = pathfinder.load_nvidia_dynamic_lib("cufft").abs_path
            except Exception:  # (an older CuPy, a system toolkit)
                path = None
        if path is None:
            path = "cufft64_11.dll" if sys.platform == "win32" else "libcufft.so.11"
        lib = ctypes.CDLL(path)
        sigs = {
            "cufftCreate": [POINTER(c_int)],
            "cufftSetAutoAllocation": [c_int, c_int],
            "cufftMakePlanMany": [
                c_int,
                c_int,
                POINTER(c_int),
                POINTER(c_int),
                c_int,
                c_int,
                POINTER(c_int),
                c_int,
                c_int,
                c_int,
                c_int,
                POINTER(c_size_t),
            ],
            "cufftMakePlanMany64": [
                c_int,
                c_int,
                POINTER(c_longlong),
                POINTER(c_longlong),
                c_longlong,
                c_longlong,
                POINTER(c_longlong),
                c_longlong,
                c_longlong,
                c_int,
                c_longlong,
                POINTER(c_size_t),
            ],
            "cufftSetWorkArea": [c_int, c_void_p],
            "cufftSetStream": [c_int, c_void_p],
            "cufftExecR2C": [c_int, c_void_p, c_void_p],
            "cufftExecC2R": [c_int, c_void_p, c_void_p],
            "cufftExecC2C": [c_int, c_void_p, c_void_p, c_int],
            "cufftDestroy": [c_int],
            "cufftGetVersion": [POINTER(c_int)],
            "cufftGetProperty": [c_int, POINTER(c_int)],
        }
        for name, args in sigs.items():
            f = getattr(lib, name)
            f.argtypes = args
            f.restype = c_int
        _lib, _lib_path = lib, path
        return lib


def cufft_info() -> dict[str, Any]:
    """cuFFT's version ('11.4.1') and the library file in use."""
    import ctypes

    lib = _library()
    v = [ctypes.c_int() for _ in range(3)]
    for i in range(3):
        _check(lib.cufftGetProperty(i, ctypes.byref(v[i])), "cufftGetProperty")
    n = ctypes.c_int()
    _check(lib.cufftGetVersion(ctypes.byref(n)), "cufftGetVersion")
    return {
        "version": ".".join(str(x.value) for x in v),
        "version_code": n.value,
        "library": _lib_path,
    }


def _yxz(shape: Sequence[int]) -> tuple[int, int, int]:
    s = tuple(int(n) for n in shape)
    if not 1 <= len(s) <= 3 or any(n < 1 for n in s):
        raise ValueError(
            f"The FFT size must have 1 to 3 positive dimensions (y, x, z), not {tuple(shape)}"
        )
    return s + (1,) * (3 - len(s))


def plan_dims(shape: Sequence[int]) -> tuple[int, ...]:
    """cuFFT's row-major dimensions of a (y, x, z) Fortran-ordered size, as the mex makes them."""
    ny, nx, nz = _yxz(shape)
    return (*(n for n in (nz, nx) if n > 1), ny)


def _destroy(lib: Any, handles: list[int], work: list[Any]) -> None:
    for h in handles:
        with contextlib.suppress(Exception):  # (at exit CUDA may be gone)
            lib.cufftDestroy(h)
    handles.clear()
    work.clear()


class _FFTPlan:
    # cuFFT plans of the given types for one size on the current GPU, sharing one work area

    def __init__(self, shape: Sequence[int], types: Sequence[int]) -> None:
        import ctypes

        cp = import_cupy()
        self._cp = cp
        self.shape = tuple(int(n) for n in shape)
        self.dims = plan_dims(self.shape)
        self.device = cp.cuda.Device().id
        self._lock = threading.Lock()
        lib = _library()
        self._lib = lib
        handles: list[int] = []
        work: list[Any] = []
        self._handles, self._work = handles, work
        self._finalizer = weakref.finalize(self, _destroy, lib, handles, work)
        sizes = []
        for t in types:
            h = ctypes.c_int()
            _check(lib.cufftCreate(ctypes.byref(h)), "cufftCreate")
            handles.append(h.value)
            _check(lib.cufftSetAutoAllocation(h.value, 0), "cufftSetAutoAllocation")
            sizes.append(self._make(cp, lib, h.value, t))
        self.work_sizes = tuple(sizes)  # what each plan needs
        self.work_bytes = max(sizes)  # the shared area
        if self.work_bytes:
            mem = cp.cuda.alloc(self.work_bytes)
            work.append(mem)
            for h in handles:
                _check(lib.cufftSetWorkArea(h, ctypes.c_void_p(mem.ptr)), "cufftSetWorkArea")

    def _make(self, cp: Any, lib: Any, h: int, t: int) -> int:
        import ctypes

        rank = len(self.dims)
        ws = ctypes.c_size_t()
        for retry in (False, True):
            n = (ctypes.c_int * rank)(*self.dims)
            r = lib.cufftMakePlanMany(h, rank, n, None, 1, 0, None, 1, 0, t, 1, ctypes.byref(ws))
            if r == _INVALID_SIZE:  # more elements than the int interface takes
                n64 = (ctypes.c_longlong * rank)(*self.dims)
                r = lib.cufftMakePlanMany64(
                    h, rank, n64, None, 1, 0, None, 1, 0, t, 1, ctypes.byref(ws)
                )
            if r == _ALLOC_FAILED and not retry:
                cp.get_default_memory_pool().free_all_blocks()
                continue
            _check(r, "cufftMakePlanMany", planning=True)
            break
        return int(ws.value)

    @property
    def closed(self) -> bool:
        return not self._finalizer.alive

    def close(self) -> None:
        """Destroy the plans and return the work area to the pool (also when garbage collected)."""
        self._finalizer()

    def __enter__(self) -> _FFTPlan:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def _ready(self, a: Any, dtype: Any, shape: tuple[int, ...], what: str) -> None:
        if self.closed:
            raise ValueError("The FFT plan is closed.")
        if a.dtype != dtype or a.shape != shape or not a.flags.f_contiguous:
            import numpy as np

            raise ValueError(
                f"{what} must be a {np.dtype(dtype).name} array of shape {shape} in Fortran order, "
                f"not {a.dtype.name} {a.shape} ({'F' if a.flags.f_contiguous else 'not F'})"
            )
        if a.device.id != self.device:
            raise ValueError(f"{what} is on GPU {a.device.id}, the plan on GPU {self.device}.")

    def _exec(self, f: Any, h: int, a: Any, b: Any, what: str, *direction: int) -> None:
        import ctypes

        _check(
            self._lib.cufftSetStream(h, ctypes.c_void_p(self._cp.cuda.get_current_stream().ptr)),
            "cufftSetStream",
        )
        _check(f(h, ctypes.c_void_p(a.data.ptr), ctypes.c_void_p(b.data.ptr), *direction), what)


class RealFFTPlan(_FFTPlan):
    """R2C and C2R plans for one size on the current GPU, sharing one work area.

    shape: (ny, nx, nz), (ny, nx) or (ny,), the size of the real arrays (float32, Fortran
    order); their spectra are complex64 of half_shape, Fortran order. inverse: also make the
    C2R plan. A plan runs on the current CuPy stream and may be used by one thread at a time
    (calls are serialized). work_sizes: (R2C, C2R).
    """

    def __init__(self, shape: Sequence[int], *, inverse: bool = True) -> None:
        super().__init__(shape, (_R2C, _C2R)[: 2 if inverse else 1])
        self.half_shape = (_yxz(self.shape)[0] // 2 + 1, *self.shape[1:])
        self.inverse_plan = inverse

    def forward(self, x: Any, out: Any = None) -> Any:
        """The half spectrum of x (R2C, out of place; x is unchanged)."""
        cp = self._cp
        self._ready(x, cp.float32, self.shape, "x")
        if out is None:
            out = cp.empty(self.half_shape, dtype=cp.complex64, order="F")
        self._ready(out, cp.complex64, self.half_shape, "out")
        with self._lock:
            self._exec(self._lib.cufftExecR2C, self._handles[0], x, out, "cufftExecR2C")
        return out

    def inverse(self, s: Any, out: Any = None) -> Any:
        """The real array of half spectrum s, unscaled (C2R; cuFFT overwrites s)."""
        cp = self._cp
        if not self.inverse_plan:
            raise ValueError("The plan was made without the inverse (inverse=False).")
        self._ready(s, cp.complex64, self.half_shape, "s")
        if out is None:
            out = cp.empty(self.shape, dtype=cp.float32, order="F")
        self._ready(out, cp.float32, self.shape, "out")
        with self._lock:
            self._exec(self._lib.cufftExecC2R, self._handles[1], s, out, "cufftExecC2R")
        return out


class ComplexFFTPlan(_FFTPlan):
    """A C2C plan (forward and inverse) for one size on the current GPU, with its work area.

    Arrays are complex64 of shape (Fortran order), in place allowed (out=z). With the cuFFT
    version MATLAB bundles, forward(complex(x)) is MATLAB's gpuArray fftn(x) bit for bit (see
    docs/decisions/0003-gpu.md); inverse is unscaled.
    """

    def __init__(self, shape: Sequence[int]) -> None:
        super().__init__(shape, (_C2C,))

    def _run(self, z: Any, out: Any, direction: int) -> Any:
        cp = self._cp
        self._ready(z, cp.complex64, self.shape, "z")
        if out is None:
            out = cp.empty(self.shape, dtype=cp.complex64, order="F")
        self._ready(out, cp.complex64, self.shape, "out")
        with self._lock:
            self._exec(self._lib.cufftExecC2C, self._handles[0], z, out, "cufftExecC2C", direction)
        return out

    def forward(self, z: Any, out: Any = None) -> Any:
        """The spectrum of z."""
        return self._run(z, out, _FORWARD)

    def inverse(self, z: Any, out: Any = None) -> Any:
        """The inverse transform of z, unscaled."""
        return self._run(z, out, _INVERSE)


def _cached(key: tuple[Any, ...], make: Any) -> Any:
    with _plans_lock:
        p = _plans.get(key)
        if p is not None and not p.closed:
            _plans.move_to_end(key)
            return p
    p = make()
    with _plans_lock:
        _plans[key] = p
        _plans.move_to_end(key)
        while len(_plans) > _plan_cache_size:
            _plans.popitem(last=False)  # (destroyed once no caller holds it)
    return p


def get_plan(shape: Sequence[int], *, inverse: bool = True) -> RealFFTPlan:
    """A cached RealFFTPlan for the size on the current GPU (the last plan_cache_size are kept)."""
    cp = import_cupy()
    shape = tuple(int(n) for n in shape)
    return _cached(
        (cp.cuda.Device().id, "real", shape, inverse), lambda: RealFFTPlan(shape, inverse=inverse)
    )


def get_complex_plan(shape: Sequence[int]) -> ComplexFFTPlan:
    """A cached ComplexFFTPlan for the size on the current GPU (sharing the cache with get_plan)."""
    cp = import_cupy()
    shape = tuple(int(n) for n in shape)
    return _cached((cp.cuda.Device().id, "complex", shape), lambda: ComplexFFTPlan(shape))


def clear_plan_cache() -> None:
    """Drop the cached plans (each is destroyed once no caller holds it)."""
    with _plans_lock:
        _plans.clear()


def set_plan_cache_size(n: int) -> None:
    """Keep the last n plans (default 2: a decon job's sizes repeat from file to file)."""
    global _plan_cache_size
    if int(n) < 0:
        raise ValueError("The plan cache size must be 0 or more.")
    with _plans_lock:
        _plan_cache_size = int(n)
        while len(_plans) > _plan_cache_size:
            _plans.popitem(last=False)


def cached_plan_bytes() -> int:
    """The work areas the cached plans hold."""
    with _plans_lock:
        return sum(p.work_bytes for p in _plans.values() if not p.closed)
