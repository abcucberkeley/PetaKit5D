"""MATLAB's ``fft``, ``ifft``, ``fftn`` and ``ifftn``, bit for bit.

The results equal MATLAB R2024b's on Linux x86-64 (``fftw-3.3.8-sse2-avx``). The transforms run on
our own FFTW build, which makes the same plans and gives the same bits as MATLAB's library, through
a driver that reproduces MATLAB's: 1D FFTW plans per dimension; lines in groups of four and the
lines MATLAB transforms alone; the real-to-complex and complex-to-real stages; the Hermitian fill of
``fftn`` (and of ``ifftn`` with ``'symmetric'``); MATLAB's exact test for conjugate-symmetric input
of ``ifft`` and ``ifftn``; the division by the length after every inverse stage; a complex result
whose imaginary parts are all zero returned as real. See ``python/cpp/fft/driver.h`` and
``python/docs/decisions/0002-fft.md``.

Arrays are taken and returned in Fortran order (MATLAB's memory layout); ``dim`` is 0-based. The
results never depend on the thread count. MATLAB itself has two modes, fixed at a session's first
FFT: with several computational threads (the default here, as production MATLAB runs) or with one
(``matlab_parallel=False``); they differ in the last bits for some lengths.

Test-only: ``PETAKIT5D_TEST_FFTW`` set to MATLAB's ``bin/glnxa64`` folder (or to its single and
double precision libraries, separated by ``os.pathsep``) runs the same driver on MATLAB's own FFTW.
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence

import numpy as np

__all__ = ["fft", "ifft", "fftn", "ifftn", "fftw_info", "use_fftw"]

# lines per piece of MATLAB's split of a stage among its threads
_GRAIN_PARALLEL = 32
_MATLAB_LIBS = ("libmwfftw3f.so.3", "libmwfftw3.so.3")
_state = {"ext": None}


def _ext():
    """The extension module, with the FFTW of PETAKIT5D_TEST_FFTW loaded once if set."""
    if _state["ext"] is None:
        from petakit5d._ext import load

        ext = load("fft")
        spec = os.environ.get("PETAKIT5D_TEST_FFTW", "")
        if spec:
            ext.set_backend(*_library_paths(spec))
        _state["ext"] = ext
    return _state["ext"]


def _library_paths(spec: str) -> tuple[str, str]:
    if os.path.isdir(spec):
        return tuple(os.path.join(spec, f) for f in _MATLAB_LIBS)
    parts = spec.split(os.pathsep)
    if len(parts) != 2:
        raise ValueError(
            "PETAKIT5D_TEST_FFTW must be a folder or two library paths separated by os.pathsep"
        )
    return parts[0], parts[1]


def use_fftw(libf: str = "", libd: str = "") -> None:
    """Runs the transforms on the FFTW libraries ``libf`` (single) and ``libd`` (double); without
    arguments on ours."""
    _ext().set_backend(libf, libd)


def fftw_info() -> dict[str, str]:
    """The FFTW in use: its name (``builtin`` or the loaded files) and version strings."""
    name, vf, vd = _ext().backend_info()
    return {"name": name, "single": vf, "double": vd}


def _threads() -> int:
    """The CPUs this process may use."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def _input(x) -> np.ndarray:
    """MATLAB's input classes: single and double, real or complex; integers and logicals become
    double."""
    a = np.asarray(x)
    if a.dtype in (np.float32, np.float64, np.complex64, np.complex128):
        pass
    elif a.dtype == np.bool_ or np.issubdtype(a.dtype, np.integer):
        if a.dtype.itemsize > 4:
            raise TypeError("Invalid data type: 64-bit integers are not supported (as in MATLAB)")
        a = a.astype(np.float64)
    else:
        raise TypeError(
            f"Invalid data type {a.dtype}: the input must be single, double, an integer or logical"
        )
    if a.ndim == 0:
        a = a.reshape(1)
    return np.asfortranarray(a)


def _first_nonsingleton(shape: tuple[int, ...]) -> int:
    """MATLAB's default dimension: the first whose length is not 1; a scalar counts as a row."""
    for k, s in enumerate(shape):
        if s != 1:
            return k
    return 1


def _count(v, what: str) -> int:
    """A non-negative integer value (MATLAB takes one of any numeric class: 8.0 is 8, 2.5 an
    error)."""
    f = float(v)
    if not (math.isfinite(f) and f >= 0 and f == math.floor(f)):
        raise ValueError(f"{what} must be a non-negative integer")
    return int(f)


def _options(threads: int | None, matlab_parallel: bool) -> tuple[int, int]:
    return (threads if threads else _threads()), (_GRAIN_PARALLEL if matlab_parallel else 0)


def _shape_back(y: np.ndarray, ndim: int, keep: int) -> np.ndarray:
    """Drops the trailing singleton dimensions the driver added beyond the input's (not below
    keep)."""
    shape = list(y.shape)
    while len(shape) > max(ndim, keep) and shape[-1] == 1:
        shape.pop()
    return y.reshape(shape, order="F") if len(shape) != y.ndim else y


def fft(
    x,
    n: int | None = None,
    dim: int | None = None,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
) -> np.ndarray:
    """``Y = fft(X, n, dim)``: the FFT along ``dim`` (0-based; default the first non-singleton
    dimension) of ``X`` zero-padded or truncated to ``n`` points (default: the size along
    ``dim``)."""
    a = _input(x)
    d = _first_nonsingleton(a.shape) if dim is None else _count(dim, "dim")
    t, g = _options(threads, matlab_parallel)
    y = _ext().fft(a, -1 if n is None else _count(n, "n"), d, t, g)
    return _shape_back(y, a.ndim, d + 1 if n is not None else 0)


def ifft(
    x,
    n: int | None = None,
    dim: int | None = None,
    symmetric: bool = False,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
) -> np.ndarray:
    """``X = ifft(Y, n, dim)``, or ``ifft(Y, n, dim, 'symmetric')`` with ``symmetric=True``: the
    inverse FFT along ``dim``; conjugate-symmetric input (MATLAB's exact test) gives a real
    result."""
    a = _input(x)
    d = _first_nonsingleton(a.shape) if dim is None else _count(dim, "dim")
    t, g = _options(threads, matlab_parallel)
    y = _ext().ifft(a, -1 if n is None else _count(n, "n"), d, bool(symmetric), t, g)
    return _shape_back(y, a.ndim, d + 1 if n is not None else 0)


def _siz(a: np.ndarray, siz: Sequence[int] | None) -> list[int]:
    """MATLAB's siz: empty for the input's size, else non-negative integers, at least ndims(X)
    of them (MATLAB's ndims: 2 at least, so a vector or a scalar needs two)."""
    if siz is None:
        return []
    s = [_count(v, "siz") for v in np.asarray(siz).ravel()]
    if not s:
        return []
    nd = a.ndim
    while nd > 2 and a.shape[nd - 1] == 1:
        nd -= 1
    if len(s) < max(nd, 2):
        raise ValueError("siz must have an element for every dimension of the input, 2 at least")
    return s


def fftn(
    x,
    siz: Sequence[int] | None = None,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
) -> np.ndarray:
    """``Y = fftn(X, siz)``: the N-D FFT of ``X`` zero-padded or truncated to ``siz`` (``None`` or
    empty: its size)."""
    a = _input(x)
    s = _siz(a, siz)
    t, g = _options(threads, matlab_parallel)
    y = _ext().fftn(a, s, t, g)
    return _shape_back(y, a.ndim, len(s))


def ifftn(
    x,
    siz: Sequence[int] | None = None,
    symmetric: bool = False,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
) -> np.ndarray:
    """``X = ifftn(Y, siz)``, or ``ifftn(Y, siz, 'symmetric')`` with ``symmetric=True``: the N-D
    inverse FFT; conjugate-symmetric input (MATLAB's exact test) gives a real result."""
    a = _input(x)
    s = _siz(a, siz)
    t, g = _options(threads, matlab_parallel)
    y = _ext().ifftn(a, s, bool(symmetric), t, g)
    return _shape_back(y, a.ndim, len(s))
