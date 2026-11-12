"""Single precision cuFFT plans with the work area under our control.

Plans work on C-contiguous arrays: a volume is (z, y, x) as everywhere in the port.

- RealFFTPlan(shape): R2C and C2R of a C array, the half spectrum along the last axis (x of a
  volume), (..., nx // 2 + 1); the two plans share one work area, the larger of the two,
  allocated from CuPy's memory pool. Size-1 axes other than the last are left out of the plan.
- ComplexFFTPlan(shape): C2C of a C array, the same rules.
- MatlabFFTPlan(shape): C2C of (z, y, x) volumes with the bits of MATLAB's gpuArray fftn and
  ifftn. MATLAB's plan follows its (y, x, z) Fortran layout, cuFFT's row-major (z, x, y), which
  no cuFFT plan over (z, y, x) memory can express (and 1D transforms in MATLAB's order don't give
  its bits at every size), so the plan keeps MATLAB's layout inside: its spectra are C arrays
  (z, x, y), the layout of MATLAB's own spectra; the forward transform transposes each z plane
  of its (real or complex) input into it, and the inverse transposes back and divides by N, as
  ifftn. The transposes replace the casts a real input needs anyway (real to complex, real part
  and scaling after the inverse), so they add no pass and no buffer; the work area is that of
  the same plan.

A plan in MATLAB's layout: ComplexFFTPlan or RealFFTPlan of the shape (nz, nx, ny) applied to
the C array np.ascontiguousarray(core.arrays.matlab_order(vol)) is the plan of the GPU decon mex
and of MATLAB's GPU code (cuFFT's row-major (z, x, y)); with cuFFT 11.0.8 it gives MATLAB's bits.

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
    "MatlabFFTPlan",
    "RealFFTPlan",
    "clear_plan_cache",
    "cufft_info",
    "get_complex_plan",
    "get_matlab_plan",
    "get_plan",
    "matlab_shape",
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


def _checked(shape: Sequence[int]) -> tuple[int, ...]:
    s = tuple(int(n) for n in shape)
    if not 1 <= len(s) <= 3 or any(n < 1 for n in s):
        raise ValueError(f"The FFT size must have 1 to 3 positive dimensions, not {tuple(shape)}")
    return s


def plan_dims(shape: Sequence[int]) -> tuple[int, ...]:
    """cuFFT's row-major dimensions of a C array of shape: size-1 axes other than the last left
    out, as the mex makes them (its (y, x, z) Fortran volume is the C array (z, x, y))."""
    s = _checked(shape)
    return (*(n for n in s[:-1] if n > 1), s[-1])


def matlab_shape(shape: Sequence[int]) -> tuple[int, ...]:
    """The C shape of MATLAB's layout of a (z, y, x) volume or (y, x) image: (nz, nx, ny), (nx,
    ny)."""
    s = _checked(shape)
    if len(s) < 2:
        raise ValueError(f"MATLAB's layout is that of a (y, x) image or (z, y, x) volume, not {s}")
    return (*s[:-2], s[-1], s[-2])


def _check_array(a: Any, dtype: Any, shape: tuple[int, ...], what: str) -> None:
    if a.dtype != dtype or a.shape != shape or not a.flags.c_contiguous:
        import numpy as np

        raise ValueError(
            f"{what} must be a C-contiguous {np.dtype(dtype).name} array of shape {shape}, not "
            f"{a.dtype.name} {a.shape}{'' if a.flags.c_contiguous else ' (not C-contiguous)'}"
        )


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
        _check_array(a, dtype, shape, what)
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

    shape: the size of the real arrays (float32, C order; a volume (nz, ny, nx)); their spectra
    are complex64 of half_shape, C order, halved along the last axis ((nz, ny, nx // 2 + 1)).
    inverse: also make the C2R plan. A plan runs on the current CuPy stream and may be used by
    one thread at a time (calls are serialized). work_sizes: (R2C, C2R).
    """

    def __init__(self, shape: Sequence[int], *, inverse: bool = True) -> None:
        super().__init__(shape, (_R2C, _C2R)[: 2 if inverse else 1])
        self.half_shape = (*self.shape[:-1], self.shape[-1] // 2 + 1)
        self.inverse_plan = inverse

    def forward(self, x: Any, out: Any = None) -> Any:
        """The half spectrum of x (R2C, out of place; x is unchanged)."""
        cp = self._cp
        self._ready(x, cp.float32, self.shape, "x")
        if out is None:
            out = cp.empty(self.half_shape, dtype=cp.complex64)
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
            out = cp.empty(self.shape, dtype=cp.float32)
        self._ready(out, cp.float32, self.shape, "out")
        with self._lock:
            self._exec(self._lib.cufftExecC2R, self._handles[1], s, out, "cufftExecC2R")
        return out


class ComplexFFTPlan(_FFTPlan):
    """A C2C plan (forward and inverse) for one size on the current GPU, with its work area.

    Arrays are complex64 of shape, C order, in place allowed (out=z); inverse is unscaled. Of
    the shape matlab_shape(s) on MATLAB's layout of a volume of shape s it is MATLAB's plan
    (see the module and MatlabFFTPlan).
    """

    def __init__(self, shape: Sequence[int]) -> None:
        super().__init__(shape, (_C2C,))

    def _run(self, z: Any, out: Any, direction: int) -> Any:
        cp = self._cp
        self._ready(z, cp.complex64, self.shape, "z")
        if out is None:
            out = cp.empty(self.shape, dtype=cp.complex64)
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


# Transposes of each z plane between our (z, y, x) and MATLAB's (z, x, y) C layout, through
# 32 x 32 tiles in shared memory (reads and writes coalesced), fused with the casts around the FFT.
# div 1: each part divided by n in double and rounded to single, as MATLAB's ifftn of single data
# (the same bits as a division in single when n is a single; only this gives MATLAB's at N =
# 181 * 1031 * 101, which isn't); div 0: no division.
_LAYOUT_SOURCE = r"""
#define TILE 32
#define ROWS 8
typedef long long idx_t;

__device__ __forceinline__ float scaled(const float v, const double n, const int div) {
    return div ? (float)((double)v / n) : v;
}

// out (z, x, y) <- in (z, y, x); real input gets a zero imaginary part
__device__ __forceinline__ float2 cplx(const float v) { return make_float2(v, 0.0f); }
__device__ __forceinline__ float2 cplx(const float2 v) { return v; }

template <typename In>
__global__ void to_matlab(const In* __restrict__ in, float2* __restrict__ out,
                          const int ny, const int nx) {
    __shared__ In tile[TILE][TILE + 1];
    const idx_t plane = (idx_t)blockIdx.z * ny * nx;
    const int x = blockIdx.x * TILE + threadIdx.x, y = blockIdx.y * TILE + threadIdx.y;
    for (int j = 0; j < TILE; j += ROWS)
        if (x < nx && y + j < ny)
            tile[threadIdx.y + j][threadIdx.x] = in[plane + (idx_t)(y + j) * nx + x];
    __syncthreads();
    const int yo = blockIdx.y * TILE + threadIdx.x, xo = blockIdx.x * TILE + threadIdx.y;
    for (int j = 0; j < TILE; j += ROWS)
        if (yo < ny && xo + j < nx)
            out[plane + (idx_t)(xo + j) * ny + yo] = cplx(tile[threadIdx.x][threadIdx.y + j]);
}

// out (z, y, x) <- in (z, x, y), divided by n (div); a real output keeps the real part
__device__ __forceinline__ void put(float* o, const float2 v, const double n, const int div) {
    *o = scaled(v.x, n, div);
}
__device__ __forceinline__ void put(float2* o, const float2 v, const double n, const int div) {
    *o = make_float2(scaled(v.x, n, div), scaled(v.y, n, div));
}

template <typename Out>
__global__ void from_matlab(const float2* __restrict__ in, Out* __restrict__ out,
                            const int ny, const int nx, const double n, const int div) {
    __shared__ float2 tile[TILE][TILE + 1];
    const idx_t plane = (idx_t)blockIdx.z * ny * nx;
    const int y = blockIdx.y * TILE + threadIdx.x, x = blockIdx.x * TILE + threadIdx.y;
    for (int j = 0; j < TILE; j += ROWS)
        if (y < ny && x + j < nx)
            tile[threadIdx.y + j][threadIdx.x] = in[plane + (idx_t)(x + j) * ny + y];
    __syncthreads();
    const int xo = blockIdx.x * TILE + threadIdx.x, yo = blockIdx.y * TILE + threadIdx.y;
    for (int j = 0; j < TILE; j += ROWS)
        if (xo < nx && yo + j < ny)
            put(out + plane + (idx_t)(yo + j) * nx + xo, tile[threadIdx.x][threadIdx.y + j],
                n, div);
}
"""
_LAYOUT_NAMES = (
    "to_matlab<float>",
    "to_matlab<float2>",
    "from_matlab<float>",
    "from_matlab<float2>",
)


class MatlabFFTPlan:
    """C2C of (z, y, x) volumes or (y, x) images with the bits of MATLAB's gpuArray fftn, ifftn.

    The plan is MATLAB's (ComplexFFTPlan of matlab_shape(shape)); its spectra are C arrays of
    spectrum_shape = matlab_shape(shape), (kz, kx, ky), the layout of MATLAB's spectra:
    element-wise products of spectra (an OTF and a transform) need no transposes. forward takes
    a float32 or complex64 volume (MATLAB's fftn of real data is the transform of complex(x));
    inverse is ifftn (divided by N) and overwrites the spectrum. to_native and from_native
    transpose a spectrum between MATLAB's layout and (kz, ky, kx).
    """

    def __init__(self, shape: Sequence[int]) -> None:
        import numpy as np

        from petakit5d.gpu import kernels

        cp = import_cupy()
        self._cp = cp
        self.shape = _checked(shape)
        self.spectrum_shape = matlab_shape(self.shape)
        nz = 1 if len(self.shape) == 2 else self.shape[0]
        if nz > 65535:
            raise ValueError(f"At most 65535 z planes, not {nz}")
        self._plan = ComplexFFTPlan(self.spectrum_shape)
        self.dims, self.device = self._plan.dims, self._plan.device
        self.work_sizes, self.work_bytes = self._plan.work_sizes, self._plan.work_bytes
        self._k = kernels.load_source(
            _LAYOUT_SOURCE, name="fft_layout.cu", name_expressions=_LAYOUT_NAMES
        )
        ny, nx = self.shape[-2], self.shape[-1]
        self._grid = ((nx + 31) // 32, (ny + 31) // 32, nz)
        self._nyx = (np.int32(ny), np.int32(nx))
        self._n = np.float64(np.prod(self.shape, dtype=np.float64))
        self._one, self._div = np.float64(1), (np.int32(0), np.int32(1))

    @property
    def closed(self) -> bool:
        return self._plan.closed

    def close(self) -> None:
        """Destroy the plan and return its work area to the pool."""
        self._plan.close()

    def __enter__(self) -> MatlabFFTPlan:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def _check(self, a: Any, dtype: Any, shape: tuple[int, ...], what: str) -> Any:
        if self.closed:
            raise ValueError("The FFT plan is closed.")
        if a is None and what == "out":
            return self._cp.empty(shape, dtype=dtype)
        _check_array(a, dtype, shape, what)
        if a.device.id != self.device:
            raise ValueError(f"{what} is on GPU {a.device.id}, the plan on GPU {self.device}.")
        return a

    def from_native(self, s: Any, out: Any = None) -> Any:
        """A (kz, ky, kx) spectrum in MATLAB's layout (kz, kx, ky)."""
        cp = self._cp
        s = self._check(s, cp.complex64, self.shape, "s")
        out = self._check(out, cp.complex64, self.spectrum_shape, "out")
        self._k["to_matlab<float2>"](self._grid, (32, 8), (s, out, *self._nyx))
        return out

    def forward(self, x: Any, out: Any = None) -> Any:
        """MATLAB's fftn(x) of a float32 or complex64 volume, in MATLAB's layout."""
        cp = self._cp
        if x.dtype not in (cp.float32, cp.complex64):
            raise ValueError(f"x must be float32 or complex64, not {x.dtype.name}")
        x = self._check(x, x.dtype, self.shape, "x")
        out = self._check(out, cp.complex64, self.spectrum_shape, "out")
        name = "to_matlab<float>" if x.dtype == cp.float32 else "to_matlab<float2>"
        self._k[name](self._grid, (32, 8), (x, out, *self._nyx))
        return self._plan.forward(out, out)

    def inverse(self, s: Any, out: Any = None, *, real: bool = False, scale: bool = True) -> Any:
        """MATLAB's ifftn(s) as a (z, y, x) volume (real: its real part, float32), from a spectrum
        in MATLAB's layout, which is overwritten. scale False: unscaled (no division by N)."""
        cp = self._cp
        s = self._check(s, cp.complex64, self.spectrum_shape, "s")
        out = self._check(out, cp.float32 if real else cp.complex64, self.shape, "out")
        self._plan.inverse(s, s)
        name = "from_matlab<float>" if real else "from_matlab<float2>"
        div = self._div[1 if scale else 0]
        self._k[name](self._grid, (32, 8), (s, out, *self._nyx, self._n, div))
        return out

    def to_native(self, s: Any, out: Any = None) -> Any:
        """A spectrum in MATLAB's layout (kz, kx, ky) as (kz, ky, kx)."""
        cp = self._cp
        s = self._check(s, cp.complex64, self.spectrum_shape, "s")
        out = self._check(out, cp.complex64, self.shape, "out")
        self._k["from_matlab<float2>"](
            self._grid, (32, 8), (s, out, *self._nyx, self._one, self._div[0])
        )
        return out


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


def get_matlab_plan(shape: Sequence[int]) -> MatlabFFTPlan:
    """A cached MatlabFFTPlan for the (z, y, x) size on the current GPU (sharing the cache)."""
    cp = import_cupy()
    shape = tuple(int(n) for n in shape)
    return _cached((cp.cuda.Device().id, "matlab", shape), lambda: MatlabFFTPlan(shape))


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
