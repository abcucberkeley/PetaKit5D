"""MATLAB's ``fft``, ``ifft``, ``fftn`` and ``ifftn`` on the port's (..., z, y, x) C arrays.

A call takes its array in the Python form and returns MATLAB's result in the Python form: the
transform along MATLAB's dimension k runs along the mapped axis (dimension 1, y, is axis -2;
dimension 2, x, is axis -1; dimension k >= 3 is axis -k), ``shape`` is the Python shape of the
padded or truncated array (MATLAB's ``siz``, mapped), and the result is
``from_matlab(MATLAB's result of to_matlab(x))``. A 1-D array is a row.

The transforms run MATLAB's driver (``python/cpp/fft/driver.h``): 1D FFTW-style plans per
dimension; the lines MATLAB transforms in groups of four and alone, numbered in MATLAB's order and
visited in the C array's memory order; the real-to-complex and complex-to-real stages; the
Hermitian fill of ``fftn`` (and of ``ifftn`` with ``symmetric``); MATLAB's exact test for
conjugate-symmetric input of ``ifft`` and ``ifftn``; the division by the length after every inverse
stage; a complex result whose imaginary parts are all zero returned as real. On the exact backend
(``petakit5d.fft``, FFTW 3.3.8 with MATLAB's codelets, the optional GPL add-on) the results are
MATLAB's bit for bit; on the permissive one (pocketfft, the default) they differ in the last bits.

MATLAB itself has two FFT modes, fixed at a session's first FFT: with several computational threads
(the default here, as production MATLAB runs) or with one (``matlab_parallel=False``); they differ
in the last bits for some lengths. The results never depend on our thread count.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from petakit5d import fft as _fft

__all__ = ["fft", "ifft", "fftn", "ifftn"]

_GRAIN_PARALLEL = 32  # lines per piece of MATLAB's split of a stage among its threads


def _input(x) -> np.ndarray:
    """MATLAB's input classes: single and double, real or complex; integers and logicals become
    double. C order."""
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
    return np.ascontiguousarray(a)


def _count(v, what: str) -> int:
    """A non-negative integer value (MATLAB takes one of any numeric class: 8.0 is 8, 2.5 an
    error)."""
    f = float(v)
    if not (math.isfinite(f) and f >= 0 and f == math.floor(f)):
        raise ValueError(f"{what} must be a non-negative integer")
    return int(f)


def _integer(v, what: str) -> int:
    f = float(v)
    if not (math.isfinite(f) and f == math.floor(f)):
        raise ValueError(f"{what} must be an integer")
    return int(f)


def _matlab_ndim(shape: tuple[int, ...]) -> int:
    """MATLAB's ndims of an array of this Python shape: its leading axes of length 1 beyond the last
    two are MATLAB's trailing singletons; 2 at least."""
    nd = len(shape)
    k = 0
    while nd - k > 2 and shape[k] == 1:
        k += 1
    return max(nd - k, 2)


def _dim_of(axis: int, ndim: int) -> int:
    """MATLAB's dimension (1-based) of a NumPy axis; a negative axis beyond the array's is a
    dimension of length 1 (MATLAB's trailing singletons)."""
    if axis >= 0:
        if axis >= ndim:
            raise ValueError(f"axis {axis} is out of range for an array of {ndim} dimensions")
        axis -= ndim
    return {-2: 1, -1: 2}.get(axis, -axis)


def _first_dim(shape: tuple[int, ...]) -> int:
    """MATLAB's default dimension: the first whose length is not 1, in MATLAB's order (y, x, z,
    ...); 2 (x) when there is none (a scalar counts as a row)."""
    msize = (shape[-2], shape[-1], *shape[-3::-1]) if len(shape) >= 2 else (1, *shape)
    for k, s in enumerate(msize):
        if s != 1:
            return k + 1
    return 2


def _options(threads: int | None, matlab_parallel: bool, backend: str | None):
    return (
        _fft._threads(threads),
        _GRAIN_PARALLEL if matlab_parallel else 0,
        _fft._backend_id(backend),
    )


def _with_dims(a: np.ndarray, ndim: int) -> np.ndarray:
    """``a`` with leading axes of length 1 up to ``ndim`` dimensions (MATLAB's trailing
    singletons)."""
    return a.reshape((1,) * (ndim - a.ndim) + a.shape) if a.ndim < ndim else a


def _back(y: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    return y.reshape(shape) if y.shape != shape else y


def fft(
    x,
    n: int | None = None,
    axis: int | None = None,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
    backend: str | None = None,
) -> np.ndarray:
    """``Y = fft(X, n, dim)``: the FFT along ``axis`` (default: MATLAB's first non-singleton
    dimension) of ``X`` zero-padded or truncated to ``n`` points (default: the length along it)."""
    a = _input(x)
    dim = _first_dim(a.shape) if axis is None else _dim_of(_integer(axis, "axis"), a.ndim)
    nn = -1 if n is None else _count(n, "n")
    t, g, b = _options(threads, matlab_parallel, backend)
    out = _out_shape_along(a, dim, nn)
    y = _fft._ext().fft(a, nn, dim - 1, t, g, b)
    return _back(y, out)


def ifft(
    x,
    n: int | None = None,
    axis: int | None = None,
    symmetric: bool = False,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
    backend: str | None = None,
) -> np.ndarray:
    """``X = ifft(Y, n, dim)``, or ``ifft(Y, n, dim, 'symmetric')`` with ``symmetric=True``: the
    inverse FFT along ``axis``; conjugate-symmetric input (MATLAB's exact test) gives a real
    result."""
    a = _input(x)
    dim = _first_dim(a.shape) if axis is None else _dim_of(_integer(axis, "axis"), a.ndim)
    nn = -1 if n is None else _count(n, "n")
    t, g, b = _options(threads, matlab_parallel, backend)
    out = _out_shape_along(a, dim, nn)
    y = _fft._ext().ifft(a, nn, dim - 1, bool(symmetric), t, g, b)
    return _back(y, out)


def _out_shape_along(a: np.ndarray, dim: int, n: int) -> tuple[int, ...]:
    """The Python shape of the result of a transform along MATLAB's dimension ``dim`` with ``n``
    points: the input's; with ``n`` (n >= 0) along the dimension's axis, which may be a new leading
    one (MATLAB's trailing singleton)."""
    if n < 0:
        return a.shape
    ax = {1: -2, 2: -1}.get(dim, -dim)
    shape = list(_with_dims(a, max(a.ndim, -ax)).shape)
    shape[ax] = n
    return tuple(shape)


def _siz(a: np.ndarray, shape: Sequence[int] | None) -> tuple[list[int], tuple[int, ...]]:
    """MATLAB's siz (its order) and the result's Python shape for ``shape``: None or empty for the
    input's size; otherwise non-negative integers, as MATLAB's siz at least ndims(X) of them (a
    shape of one length is a row, MATLAB's [1 n]; a leading axis of length 1 beyond the last two is
    MATLAB's trailing singleton)."""
    if shape is None:
        return [], a.shape
    s = [_count(v, "shape") for v in np.asarray(shape).ravel()]
    if not s:
        return [], a.shape
    if max(len(s), 2) < _matlab_ndim(a.shape):
        raise ValueError("shape must have a length for every dimension of the input")
    if len(s) == 1:
        msiz = [1, s[0]]
    else:
        msiz = [s[-2], s[-1], *s[-3::-1]]
    return msiz, tuple(s)


def fftn(
    x,
    shape: Sequence[int] | None = None,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
    backend: str | None = None,
) -> np.ndarray:
    """``Y = fftn(X, siz)``: the N-D FFT of ``X`` zero-padded or truncated to ``shape`` (the Python
    shape; ``None`` or empty: its own)."""
    a = _input(x)
    msiz, out = _siz(a, shape)
    t, g, b = _options(threads, matlab_parallel, backend)
    y = _fft._ext().fftn(a, msiz, t, g, b)
    return _back(y, out)


def ifftn(
    x,
    shape: Sequence[int] | None = None,
    symmetric: bool = False,
    *,
    threads: int | None = None,
    matlab_parallel: bool = True,
    backend: str | None = None,
) -> np.ndarray:
    """``X = ifftn(Y, siz)``, or ``ifftn(Y, siz, 'symmetric')`` with ``symmetric=True``: the N-D
    inverse FFT; conjugate-symmetric input (MATLAB's exact test) gives a real result."""
    a = _input(x)
    msiz, out = _siz(a, shape)
    t, g, b = _options(threads, matlab_parallel, backend)
    y = _fft._ext().ifftn(a, msiz, bool(symmetric), t, g, b)
    return _back(y, out)
