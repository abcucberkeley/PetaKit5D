"""The registration's normalized cross-correlation in the window of the allowed shifts:
``normxcorr3_fast``, ``normxcorr3_max_shift`` and ``normxcorr2_max_shift`` of PetaKit5D (MATLAB
4a3bef3), with ``normxcorr3_window_mex`` behind them (``petakit5d._ext.xcorr``).

``normxcorr3_fast`` is the lab's adaptation of Daniel Eaton's normxcorr3 (MATLAB File Exchange
73946, BSD license: ``THIRD_PARTY_NOTICES.md``), which follows J. P. Lewis, "Fast Normalized
Cross-Correlation", Vision Interface 1995: the correlation by FFTs, the window sums of the image
from running sums.

The NCC is MATLAB's: the template and the image each less its mean (in their class), their
correlation by FFTs in the class of the data, the window sums (by cumulative sums along y, then x,
then z), the normalization and the NCC in double; 0 where a window's variance is at most
``eps(class)`` times the image's (``var(A(:), 1)``) or where the template is flat; the result in
the class of the image.

Arrays are the port's: volumes ``(z, y, x)``, images ``(y, x)`` (a volume of one slice),
C-contiguous float32 or float64 (another layout is copied, another class raises). Boxes are
``(z0, y0, x0, z1, y1, x1)`` in the full output of shape ``T.shape + A.shape - 1``, 0-based and
end-exclusive (a 2D box ``(y0, x0, y1, x1)`` for images); shifts and offsets are ``(z, y, x)``.
Inside, the kernel works in MATLAB's layout (y fastest), in which MATLAB's FFT stages, groups of
lines and cumulative sums are defined: the inputs are moved there at entry and the window back at
exit.

FFTs run on the CPU FFT backends of ``petakit5d.fft`` (``backend=``: ``"permissive"``, pocketfft,
the default; ``"exact"``, the FFTW add-on), in the process's one FFT state: with ``"exact"`` the
results are MATLAB's bit for bit (MATLAB's FFT mode with several threads;
``matlab_parallel=False`` for a MATLAB that started on one thread), with ``"permissive"`` they
differ from MATLAB's by the transforms' rounding (``docs/decisions/0006-fft-backends.md``). The
template's and the image's means, sum, standard deviation and variance are MATLAB's over ``x(:)``
in its order at one thread, as MATLAB's ``fixed_order_stat`` computes them. Results do not depend
on ``threads``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

__all__ = [
    "XcorrError",
    "max_shift_box",
    "normxcorr2_max_shift",
    "normxcorr3_fast",
    "normxcorr3_max_shift",
    "normxcorr3_window",
]

_MEX_ID = "normxcorr3_window_mex:inputError"
_SHAPES = ("full", "same", "valid")
_state: dict = {"ext": None}


class XcorrError(ValueError):
    """An error MATLAB raises; ``identifier`` is its MATLAB error identifier (``''`` for
    normxcorr3_fast's errors, which have none)."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _ext():
    if _state["ext"] is None:
        from petakit5d._ext import load

        _state["ext"] = load("xcorr")
    return _state["ext"]


def _backend(name: str | None) -> int:
    """The backend id for ``name`` (None: petakit5d.fft's default); the exact one is loaded on
    first use into the process's FFT state, which this module shares with petakit5d._ext.fft."""
    from petakit5d import fft

    name = fft.backend() if name is None else name
    if name not in fft.BACKENDS:
        raise ValueError(f"unknown FFT backend {name!r}: use one of {fft.BACKENDS}")
    return fft._backend_id(name)


def _threads(threads: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(threads)


def _grain(matlab_parallel: bool) -> int:
    return int(_ext().GRAIN_PARALLEL) if matlab_parallel else 0


def _array(x, name: str, identifier: str = "") -> np.ndarray:
    a = np.asarray(x)
    if a.ndim > 3:
        raise XcorrError(identifier, "A and T must be no more than 3 dimensional matrices")
    if a.dtype not in (np.float32, np.float64):
        raise TypeError(f"{name} must be float32 or float64 (MATLAB's single or double)")
    if a.ndim < 2:
        a = a.reshape(1, -1)  # a vector is a row
    return np.ascontiguousarray(a)


def _shape3(a: np.ndarray) -> tuple[int, int, int]:
    return (1, *a.shape) if a.ndim == 2 else a.shape


def _full_shape(t: np.ndarray, a: np.ndarray) -> tuple[int, int, int]:
    return tuple(p + q - 1 for p, q in zip(_shape3(t), _shape3(a), strict=True))


def _win(box, full: tuple[int, int, int], identifier: str, message: str) -> list[int]:
    """MATLAB's window [s1 s2 s3 t1 t2 t3] (y, x, z; 1-based, inclusive) of a box of the full
    output (``None``: all of it)."""
    if box is None:
        z0, y0, x0, (z1, y1, x1) = 0, 0, 0, full
    else:
        b = list(np.asarray(box, dtype=object).ravel())
        if len(b) == 4:
            b = [0, b[0], b[1], 1, b[2], b[3]]
        ok = len(b) == 6 and all(
            isinstance(v, int | np.integer) or (isinstance(v, float) and v.is_integer()) for v in b
        )
        if not ok:
            raise XcorrError(identifier, message)
        z0, y0, x0, z1, y1, x1 = (int(v) for v in b)
        lo, hi = (z0, y0, x0), (z1, y1, x1)
        if not all(0 <= s < e <= n for s, e, n in zip(lo, hi, full, strict=True)):
            raise XcorrError(identifier, message)
    return [y0 + 1, x0 + 1, z0 + 1, y1, x1, z1]


def _window_shape(win: Sequence[int], two_d: bool) -> tuple[int, ...]:
    wy, wx, wz = (win[3 + i] - win[i] + 1 for i in range(3))
    return (wy, wx) if two_d else (wz, wy, wx)


def normxcorr3_window(
    T0,
    A0,
    box=None,
    fast_fft: bool = False,
    *,
    threads: int | None = None,
    backend: str | None = None,
    matlab_parallel: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``[corr, sA, sA2, sT] = normxcorr3_window_mex('ncc', T0, A0, bbox, nThreads, grain, fast)``.

    In the box of the full output: ``corr``, the correlation of ``T0`` and ``A0`` (single), and in
    double the sums of ``A0`` and ``A0 .* A0`` over the image's elements that the template overlaps
    at each position, and of ``T0`` over the template's elements that overlap the image. ``T0``
    and ``A0`` are float32, ``T0`` not larger than ``A0``, the full output at least 2 long along y
    and x. ``fast_fft``: smaller FFTs (2, 3, 5, 7-smooth sizes without wrap-around in the box),
    which changes ``corr`` by rounding.
    """
    a0, t0 = np.asarray(A0), np.asarray(T0)
    if t0.dtype != np.float32 or a0.dtype != np.float32 or t0.size == 0 or a0.size == 0:
        raise XcorrError(_MEX_ID, "T and A must be real, non-empty single arrays.")
    t = _array(t0, "T0", _MEX_ID)
    a = _array(a0, "A0", _MEX_ID)
    if any(p > q for p, q in zip(_shape3(t), _shape3(a), strict=True)):
        raise XcorrError(_MEX_ID, "The template must not be larger than the image.")
    full = _full_shape(t, a)
    if full[1] < 2 or full[2] < 2:
        raise XcorrError(_MEX_ID, "The full output must have at least 2 elements along y and x.")
    win = _win(
        box,
        full,
        _MEX_ID,
        "bbox must hold integers with 1 <= s <= t <= size of the full output.",
    )
    two_d = t.ndim == 2 and a.ndim == 2
    out = _ext().window(
        t, a, win, _threads(threads), _grain(matlab_parallel), bool(fast_fft), _backend(backend)
    )
    shape = _window_shape(win, two_d)
    return tuple(v.reshape(shape) for v in out)


def _ncc(T, A, box, fast_fft, threads, backend, matlab_parallel):
    """normxcorr3_fast in the box (the full output for None): C (window, class of A) and the
    0-based index of its first maximum in MATLAB's order."""
    t = _array(T, "T")
    a = _array(A, "A")
    if any(p > q for p, q in zip(_shape3(t), _shape3(a), strict=True)):
        raise XcorrError("", "template must be smaller than image")
    if t.dtype != a.dtype:
        raise TypeError("T and A must have the same class")
    full = _full_shape(t, a)
    win = _win(
        box,
        full,
        "",
        "BBOX must hold integers [s1, s2, s3, t1, t2, t3] with 1 <= s <= t <= the size of the "
        "full output",
    )
    nt = _threads(threads)
    be = _backend(backend)
    grain = _grain(matlab_parallel)
    t3 = t.reshape(_shape3(t))
    a3 = a.reshape(_shape3(a))
    ext = _ext()
    if a.dtype == np.float32 and full[1] > 1 and full[2] > 1:
        # the mex's path
        C, imax, _, _ = ext.ncc(t3, a3, win, nt, grain, bool(fast_fft), be)
    else:
        # MATLAB's code: the correlation by fftn and ifftn of the whole output, in the class
        C, imax = _ncc_by_fftn(ext, t3, a3, full, win, nt, be, matlab_parallel, backend)
    two_d = t.ndim == 2 and a.ndim == 2
    return C.reshape(_window_shape(win, two_d)), int(imax)


def _times(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """MATLAB's a .* b of spectra, each complex or real (fftn returns an all-real result as
    real): the complex product without fused multiply-adds, a real factor times each part"""
    ac, bc = np.iscomplexobj(a), np.iscomplexobj(b)
    if not ac and not bc:
        return np.multiply(a, b)
    if ac and bc:
        out = np.empty(a.shape, dtype=a.dtype)
        out.real = np.subtract(np.multiply(a.real, b.real), np.multiply(a.imag, b.imag))
        out.imag = np.add(np.multiply(a.real, b.imag), np.multiply(a.imag, b.real))
        return out
    r, c = (b, a) if ac else (a, b)
    out = np.empty(c.shape, dtype=c.dtype)
    out.real = np.multiply(r, c.real)
    out.imag = np.multiply(r, c.imag)
    return out


def _ncc_by_fftn(ext, t3, a3, full, win, nt, be, matlab_parallel, backend):
    from petakit5d.matlabcompat.fft import fftn, ifftn

    T0, A0, stats = ext.prepare(t3, a3, nt)
    kw = {"threads": nt, "matlab_parallel": matlab_parallel, "backend": backend}
    fa = fftn(A0, full, **kw)
    ft = fftn(np.ascontiguousarray(T0[::-1, ::-1, ::-1]), full, **kw)
    corr = ifftn(_times(fa, ft), **kw)
    if np.iscomplexobj(corr):
        corr = corr.real
    corr = np.ascontiguousarray(
        corr[win[2] - 1 : win[5], win[0] - 1 : win[3], win[1] - 1 : win[4]], dtype=a3.dtype
    )
    return ext.finish(corr, T0, A0, win, stats, nt)


def normxcorr3_fast(
    T,
    A,
    shape: str = "full",
    box=None,
    fast_fft: bool = False,
    *,
    threads: int | None = None,
    backend: str | None = None,
    matlab_parallel: bool = True,
) -> np.ndarray:
    """``C = normxcorr3_fast(T, A, shape, bbox, fastFFT)``: the NCC of the template ``T`` with the
    image ``A``.

    ``shape`` ``'full'`` (size ``T + A - 1``), ``'same'`` (size of ``A``) or ``'valid'`` (size
    ``A - T + 1``); with ``box`` (``shape`` must be ``'full'``) only that part of the full output.
    Single data whose full output has at least 2 elements along y and x take the mex's path;
    double data and thin outputs MATLAB's code (``fftn`` and ``ifftn`` of the whole output), with
    the same results. ``fast_fft`` lets the mex's path use smaller FFTs.
    """
    t = _array(T, "T")
    a = _array(A, "A")
    t3, a3 = _shape3(t), _shape3(a)
    if any(p > q for p, q in zip(t3, a3, strict=True)):
        raise XcorrError("", "template must be smaller than image")
    if box is not None and str(shape).lower() != "full":
        raise XcorrError("", "BBOX requires SHAPE 'full'")
    C, _ = _ncc(t, a, box, fast_fft and box is not None, threads, backend, matlab_parallel)
    if box is not None:
        return C
    s = str(shape).lower()
    if s not in _SHAPES:
        raise XcorrError("", f"unknown SHAPE {shape}, assuming FULL by default")
    if s == "full":
        return C
    if s == "same":
        lo = [(n - 1) // 2 for n in t3]
        sl = tuple(slice(p, p + m) for p, m in zip(lo, a3, strict=True))
    else:
        # MATLAB's C(szT(1) : end - szT(1) + 1, ...) along each dimension
        sl = tuple(slice(n - 1, f - n + 1) for n, f in zip(t3, _full_shape_of(t3, a3), strict=True))
    C3 = C.reshape((1, *C.shape)) if C.ndim == 2 else C
    out = np.ascontiguousarray(C3[sl])
    return out.reshape(out.shape[1:]) if C.ndim == 2 else out


def _full_shape_of(t3, a3):
    return tuple(p + q - 1 for p, q in zip(t3, a3, strict=True))


def _per_axis_bounds(max_shifts, n: int) -> tuple[list[float], list[float]]:
    """(lower, upper) per axis of the last n axes from a vector of n or 3 values (``[-m, m]``) or
    a (2, n or 3) array of lower and upper rows, in (z, y, x) order."""
    m = np.asarray(max_shifts, dtype=np.float64)
    if m.ndim == 1:
        m = np.stack([-m, m])
    if m.ndim != 2 or m.shape[0] != 2 or m.shape[1] < n:
        raise ValueError(
            f"max_shifts must be {n} or 3 values or a (2, {n}) array of lower and upper bounds"
        )
    m = m[:, -n:]
    return [float(v) for v in m[0]], [float(v) for v in m[1]]


def max_shift_box(t_shape: Sequence[int], a_shape: Sequence[int], max_shifts) -> tuple[int, ...]:
    """The box of the full output of the shifts within ``max_shifts`` (``(z, y, x)``: one value per
    axis for ``[-m, m]``, or rows of lower and upper bounds), as normxcorr3_max_shift computes it:
    from ``max(1, ceil(lower + szT))`` to ``min(floor(upper + szT), szT + szA - 1)`` (MATLAB's
    1-based window)."""
    t3 = (1, *t_shape) if len(t_shape) == 2 else tuple(t_shape)
    a3 = (1, *a_shape) if len(a_shape) == 2 else tuple(a_shape)
    lower, upper = _per_axis_bounds(max_shifts, 3)
    full = _full_shape_of(t3, a3)
    start = [max(1, math.ceil(lo + n)) - 1 for lo, n in zip(lower, t3, strict=True)]
    stop = [min(math.floor(hi + n), f) for hi, n, f in zip(upper, t3, full, strict=True)]
    return (*start, *stop)


def normxcorr3_max_shift(
    T,
    A,
    max_shifts,
    fast_fft: bool = False,
    *,
    threads: int | None = None,
    backend: str | None = None,
    matlab_parallel: bool = True,
) -> tuple[tuple[int, int, int], np.floating, np.ndarray]:
    """``[max_off, max_corr, C] = normxcorr3_max_shift(T, A, maxShifts, fastFFT)``.

    ``max_shifts`` in ``(z, y, x)``: one value per axis (shifts from ``-m`` to ``m``) or a (2, 3)
    array of lower and upper bounds. Returns the offset ``(z, y, x)`` of the template in the image
    at the first maximum of the NCC in MATLAB's order (``max_off``), that maximum, and the NCC in
    the window of the allowed shifts (``C``, ``(z, y, x)``; ``(y, x)`` for images).
    """
    t = _array(T, "T")
    a = _array(A, "A")
    t3, a3 = _shape3(t), _shape3(a)
    box = max_shift_box(t3, a3, max_shifts)
    C, imax = _ncc(t, a, box, fast_fft, threads, backend, matlab_parallel)
    from petakit5d.core.arrays import matlab_index_to_zyx

    wshape = tuple(box[3 + i] - box[i] for i in range(3))
    idx = matlab_index_to_zyx(imax, wshape)
    off = tuple(int(i + s + 1 - n) for i, s, n in zip(idx, box[:3], t3, strict=True))
    return off, C.reshape(-1)[_c_index(idx, wshape)], C


def _c_index(idx, shape) -> int:
    return int(np.ravel_multi_index(idx, shape))


def normxcorr2_max_shift(
    T, A, max_shifts, *, normxcorr2
) -> tuple[tuple[int, int, int], float, np.ndarray]:
    """``[max_off, max_corr, C] = normxcorr2_max_shift(T, A, maxShifts)`` for images ``(y, x)``.

    ``normxcorr2`` is the function that gives MATLAB's ``normxcorr2(T, A)`` (the full output, of
    size ``T + A - 1``, as an image); ``max_shifts`` gives the bounds along y and x (``(y, x)`` or
    ``(z, y, x)`` values, or rows of lower and upper bounds). Returns the offset ``(0, y, x)`` at
    the first maximum in MATLAB's order (NaN ignored), the maximum and the window of ``C``.
    """
    t = np.asarray(T)
    C = np.asarray(normxcorr2(T, A))
    lower, upper = _per_axis_bounds(max_shifts, 2)
    sz_t = t.shape[-2:]
    s = [max(1, math.ceil(lo + n)) for lo, n in zip(lower, sz_t, strict=True)]
    e = [min(math.floor(hi + n), c) for hi, n, c in zip(upper, sz_t, C.shape, strict=True)]
    Cw = np.ascontiguousarray(C[s[0] - 1 : e[0], s[1] - 1 : e[1]])
    # MATLAB's order: y fastest
    v = Cw.T.reshape(-1)
    k = 0 if np.all(np.isnan(v)) else int(np.nanargmax(v))
    y, x = k % Cw.shape[0], k // Cw.shape[0]
    off = (0, y + s[0] - sz_t[0], x + s[1] - sz_t[1])
    return off, Cw[y, x], Cw
