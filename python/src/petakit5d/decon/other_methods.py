"""RLdecon's 'original' method: MATLAB's ``deconvlucy`` (the Image Processing Toolbox's accelerated,
damped Richardson-Lucy deconvolution), reproduced from its documentation, the published method
and black-box probes of MATLAB R2024b (task B11; no toolbox source was read).

The method (Richardson 1972, Lucy 1974; acceleration of Biggs and Andrews, Appl. Opt. 36, 1766,
1997; damping of White, "Image restoration using the damped Richardson-Lucy method", 1994): with
the OTF ``H`` of the PSF, the data ``I``, the weights ``W`` and the read-out noise ``R``, each
iteration predicts ``Y = max(J + lambda (J - J_prev), 0)``, re-blurs it, ``B = real(ifftn(H .*
fftn(Y))) + R``, takes the ratio ``W (I + R) / B`` (damped where the fit is within ``DAMPAR``),
correlates it with the PSF and normalizes by the correlation of the weights:
``J_new = max(Y .* real(ifftn(conj(H) .* fftn(ratio))) ./ (real(ifftn(conj(H) .* fftn(W))) +
sqrt(eps)), 0)``. ``lambda`` is the dot product of the last two changes ``J_new - Y`` over the
squared norm of the older one (``+ eps``), clamped to [0, 1]; it is 0 until two changes exist.

What the probes showed (``tests/decon/test_other_methods.py`` has a case for each):

- **Precision**: the iterations run in double for every class; the estimate is returned in the
  class of the data. Integer data is scaled to [0, 1] first (MATLAB's ``im2double``: uint8 / 255,
  uint16 / 65535, int16 (x + 32768) / 65535) and back with rounding and saturation; ``DAMPAR``
  and ``READOUT`` must have the data's class and are scaled the same way.
- **The OTF** is MATLAB's ``psf2otf`` of the PSF in double: the PSF padded at the end to the
  estimate's size, its center ``floor(size / 2)`` shifted to the origin, ``fftn``; real when its
  imaginary parts are round-off, ``max|imag| <= n log2(n) eps max|OTF|`` with ``n`` the number of
  PSF elements (bounded by probes; ``psf2otf_real``). The PSF is not normalized: the weights'
  correlation normalizes the update.
- **Zeros**: a re-blurred value of exactly 0 becomes ``eps`` and the ratio is ``W (I + R) / B +
  eps``; ``W (I + R)`` is clamped at 0 (negative data counts as 0, the first estimate is the data).
- **Damping** (``DAMPAR`` > 0, N = 10): ``g = min((W I log(r) + B - W I) / (DAMPAR^2 / 2), 1)`` and
  the ratio ``r`` becomes ``1 + g^9 (10 - 9 g) (r - 1)``. The log is NumPy's (the C library's last
  bits; MATLAB's ``log`` is not written out in ``matlabcompat``); the power is MATLAB's.
- **Subsampling** (``SUBSMPL`` s > 1, a PSF s times finer than the data): the estimate lives on the
  grid s times finer along every axis, starts as the data replicated into s^d blocks, is re-blurred
  there and averaged over the blocks; the ratio is replicated back. The result has the fine size.
- **The cell form** (``deconvlucy_state``): MATLAB's ``J = deconvlucy({I}, ...)`` returns ``{I, the
  last estimate, the one before (double), the last two changes as numel x 2 (double)}``, and a call
  with that cell continues the iterations exactly (the acceleration included).
- **Arguments**: ``NUMIT`` default 10 (also for []), a positive scalar, fractions run
  ``floor(NUMIT)`` iterations; ``DAMPAR``, ``WEIGHT`` and ``READOUT`` [] are their defaults (0,
  ones, 0); ``WEIGHT`` is double of the data's size; ``SUBSMPL`` a positive integer; the PSF not
  all zero and not larger than the (fine) estimate. Errors carry MATLAB's identifiers.

**Where MATLAB's bits are not reached**: MATLAB's dot products for ``lambda`` are BLAS's (their
order follows the CPU and the thread count); here they are fixed-order double sums in the arrays'
memory order (``kernels.decon.accel_lambda``), the same on every CPU and at every thread count. A
``lambda`` one unit in the last place away moves double results by a few units in the last place
(at most 26 ulp after 5 iterations in the probes); single and integer results were bit-identical to
MATLAB's in every probe on the exact FFT backend. The permissive backend is within decision record
0006's tolerance.

``rl_decon_original`` is RLdecon's branch for the method (``deconvlucy`` of the background-free
single data with the normalized PSF, then ``deconBbox`` and ``scaleFactor``). Arrays are the
port's: volumes ``(z, y, x)`` and images ``(y, x)``, C-contiguous; a volume with one z slice is
MATLAB's 2D image. The iterations run on half spectra (``petakit5d.fft.HalfSpectrum``; on the exact
backend MATLAB's ``real(ifftn(H .* fftn(Y)))`` bit for bit) and in threads; nothing depends on the
thread count. MATLAB runs this method on the CPU only (RLdecon passes no GPU array), and so does the
port.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any, NamedTuple

import numpy as np

__all__ = [
    "DeconvlucyError",
    "LucyState",
    "deconvlucy",
    "deconvlucy_state",
    "psf2otf_real",
    "rl_decon_original",
]

_EPS = 2.220446049250313e-16
_SQRT_EPS = math.sqrt(_EPS)  # 1.4901161193847656e-08, exact (eps is a power of two)
_DAMP_N = 10  # White's exponent
# MATLAB's classes of deconvlucy's data and their im2double ranges
_RANGES = {
    np.dtype(np.uint8): (0.0, 255.0),
    np.dtype(np.uint16): (0.0, 65535.0),
    np.dtype(np.int16): (-32768.0, 32767.0),
    np.dtype(np.float32): None,
    np.dtype(np.float64): None,
}
_CLASS = {
    np.dtype(np.uint8): "uint8",
    np.dtype(np.uint16): "uint16",
    np.dtype(np.int16): "int16",
    np.dtype(np.float32): "single",
    np.dtype(np.float64): "double",
}


class DeconvlucyError(ValueError):
    """An error deconvlucy raises, with MATLAB's identifier."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


class LucyState(NamedTuple):
    """MATLAB's cell form of deconvlucy: the data as given, the last estimate and the one before
    (double, the fine size with subsampling), and ``internal``, the last two changes of the
    estimate ``(2, *estimate.shape)`` (double; the latest first; MATLAB's numel x 2 array is
    ``internal`` reshaped in MATLAB's linear order)."""

    image: np.ndarray
    estimate: np.ndarray
    previous: np.ndarray
    internal: np.ndarray


# ---------------------------------------------------------------------------------------------
# arguments
# ---------------------------------------------------------------------------------------------


def _empty(v: Any) -> bool:
    return v is None or (isinstance(v, np.ndarray) and v.size == 0)


def _matlab_2d(a: np.ndarray) -> np.ndarray:
    return a.reshape(a.shape[1:]) if a.ndim == 3 and a.shape[0] == 1 else a


def _data(i: Any) -> np.ndarray:
    a = np.asarray(i)
    if a.dtype not in _RANGES:
        raise DeconvlucyError(
            "MATLAB:deconvlucy:invalidType",
            "Expected input number 1, I, to be one of these types: uint8, uint16, int16, "
            f"single, double; it is {a.dtype}.",
        )
    if a.ndim not in (2, 3):
        raise DeconvlucyError(
            "MATLAB:deconvlucy:invalidType", "I must be a (z, y, x) volume or a (y, x) image."
        )
    return np.ascontiguousarray(_matlab_2d(a))


def _to_double(x: np.ndarray, dt: np.dtype) -> np.ndarray:
    """MATLAB's im2double of an array of class ``dt`` (integers scaled to [0, 1]), a new array"""
    r = _RANGES[dt]
    d = np.array(x, dtype=np.float64)
    if r is None:
        return d
    lo, hi = r
    if lo != 0.0:
        d -= lo
    d /= hi - lo
    return d


def _from_double(j: np.ndarray, dt: np.dtype) -> np.ndarray:
    """the estimate in the data's class (MATLAB's im2uint8/im2uint16/im2int16, single, double)"""
    from petakit5d.core import dtypes

    r = _RANGES[dt]
    if r is None:
        return j if dt == np.float64 else np.asarray(dtypes.cast(j, "single"))
    lo, hi = r
    v = j * (hi - lo) if lo == 0.0 else j * (hi - lo) + lo
    return np.asarray(dtypes.cast(v, _CLASS[dt]))


def _iterations(numit: Any) -> int:
    if _empty(numit):
        return 10
    v = np.asarray(numit)
    if v.size != 1 or v.dtype.kind not in "iuf":
        raise DeconvlucyError("MATLAB:deconvlucy:expectedScalar", "NUMIT must be a scalar.")
    f = float(v.reshape(-1)[0])
    if not f > 0:
        raise DeconvlucyError(
            "MATLAB:deconvlucy:expectedPositive", "Expected input number 3, NUMIT, to be positive."
        )
    if not math.isfinite(f):
        raise DeconvlucyError(
            "MATLAB:deconvlucy:expectedFinite", "Expected input number 3, NUMIT, to be finite."
        )
    return int(math.floor(f))


def _same_class(v: Any, dt: np.dtype, name: str, identifier: str) -> float:
    """DAMPAR or READOUT: a scalar of the data's class (a Python number is MATLAB's double),
    scaled as the data"""
    if _empty(v):
        return 0.0
    a = np.asarray(v)
    cls = np.dtype(np.float64) if isinstance(v, (int, float)) else a.dtype
    if cls != dt:
        raise DeconvlucyError(identifier, f"{name} must be the same class as the input image.")
    if a.size != 1:
        raise DeconvlucyError("MATLAB:deconvlucy:expectedScalar", f"{name} must be a scalar.")
    return float(_to_double(a.reshape(1), dt)[0])


def _weight(w: Any, shape: tuple[int, ...]) -> np.ndarray | None:
    if _empty(w):
        return None
    a = np.asarray(w)
    if a.dtype != np.float64:
        raise DeconvlucyError(
            "MATLAB:deconvlucy:invalidType",
            "Expected input number 5, WEIGHT, to be one of these types: double.",
        )
    a = _matlab_2d(a)
    if a.shape != shape:
        raise DeconvlucyError(
            "images:deconvlucy:weightMustBeSameSizeAsImage",
            "WEIGHT must be the same size as the input image.",
        )
    return np.ascontiguousarray(a)


def _subsampling(s: Any) -> int:
    if _empty(s):
        return 1
    v = np.asarray(s)
    if v.size != 1 or v.dtype.kind not in "iuf":
        raise DeconvlucyError("MATLAB:deconvlucy:expectedScalar", "SUBSMPL must be a scalar.")
    f = float(v.reshape(-1)[0])
    if not f > 0:
        raise DeconvlucyError(
            "MATLAB:deconvlucy:expectedPositive",
            "Expected input number 7, SUBSMPL, to be positive.",
        )
    if f != math.floor(f) or not math.isfinite(f):
        raise DeconvlucyError(
            "MATLAB:deconvlucy:expectedInteger",
            "Expected input number 7, SUBSMPL, to be integer-valued.",
        )
    return int(f)


def _psf(psf: Any, fine: tuple[int, ...]) -> np.ndarray:
    p = np.asarray(psf)
    if p.dtype.kind not in "iuf" or p.size == 0:
        raise DeconvlucyError("MATLAB:deconvlucy:invalidType", "PSF must be a numeric array.")
    p = _matlab_2d(p.astype(np.float64))
    if not np.any(p != 0):
        raise DeconvlucyError(
            "images:deconvlucy:psfMustNotBeZeroEverywhere", "PSF must not be zero everywhere."
        )
    ps = (1,) * (len(fine) - p.ndim) + p.shape if p.ndim <= len(fine) else p.shape
    if len(ps) != len(fine) or any(a > b for a, b in zip(ps, fine, strict=True)):
        raise DeconvlucyError(
            "images:deconvlucy:psfMustBeSmallerThanImage",
            "PSF must be smaller than or equal to the image size.",
        )
    return np.ascontiguousarray(p.reshape(ps))


# ---------------------------------------------------------------------------------------------
# the OTF
# ---------------------------------------------------------------------------------------------


def _round_off(psf_size: int) -> float:
    """psf2otf's bound for round-off in the imaginary part, relative to max|OTF|"""
    return psf_size * math.log2(psf_size) * _EPS if psf_size > 1 else 0.0


def _is_round_off(otf: np.ndarray, psf_size: int) -> bool:
    if not np.iscomplexobj(otf):
        return True
    im = float(np.max(np.abs(otf.imag)))
    return im <= _round_off(psf_size) * float(np.max(np.abs(otf)))


def psf2otf_real(
    psf: Any,
    out_size: Sequence[int],
    *,
    backend: str | None = None,
    threads: int | None = None,
) -> np.ndarray:
    """MATLAB's ``psf2otf(psf, out_size)``: the OTF in double (``decon.otf.psf2otf`` of the PSF in
    double), real where the imaginary parts are round-off (see the module)."""
    from petakit5d.decon.otf import psf2otf

    p = np.asarray(psf, dtype=np.float64)
    o = psf2otf(p, out_size, backend=backend, threads=threads)
    if np.iscomplexobj(o) and _is_round_off(o, p.size):
        return np.ascontiguousarray(o.real)
    return o


def _half_otf(hs, psf: np.ndarray) -> np.ndarray:
    """the half spectrum of the OTF; its imaginary parts zeroed where they are round-off (MATLAB's
    real OTF: a real factor multiplies both parts of the spectrum alike)"""
    from petakit5d.decon.otf import pad_shift_psf

    h = hs.forward(pad_shift_psf(psf, hs.shape))
    if _is_round_off(h, psf.size):
        h.imag[...] = 0.0
    return h


# ---------------------------------------------------------------------------------------------
# element-wise steps in threads
# ---------------------------------------------------------------------------------------------


def _parallel(fn, n: int, threads: int) -> None:
    from petakit5d.decon.rl import _parallel as par

    par(fn, n, threads)


def _max0(x: np.ndarray) -> None:
    """MATLAB's ``max(x, 0)`` in place: NaN and negative values to +0, -0 kept"""
    np.copyto(x, 0.0, where=~(x >= 0))


def _ratio(rb, wi, readout, d22, threads) -> np.ndarray:
    """``W (I + R) / B + eps`` with ``B + R`` (0 as eps), damped, in place in ``rb``"""
    from petakit5d.matlabcompat.elementary import power

    fb, fw = rb.reshape(-1), wi.reshape(-1)

    def piece(s: int, e: int) -> None:
        b, w = fb[s:e], fw[s:e]
        if readout:
            b += readout
        np.copyto(b, _EPS, where=b == 0)
        with np.errstate(all="ignore"):
            r = w / b
            r += _EPS
            if d22:
                g = w * np.log(r)
                g += b
                g -= w
                g /= d22
                np.copyto(g, 1.0, where=~(g <= 1))  # MATLAB's min(g, 1): NaN gives 1
                f = power(g, _DAMP_N - 1) * (_DAMP_N - (_DAMP_N - 1) * g)
                r = 1 + f * (r - 1)
        b[...] = r

    _parallel(piece, fb.size, threads)
    return rb


def _update(y, cc, scale, threads) -> np.ndarray:
    """``max(Y .* CC ./ scale, 0)`` in place in ``cc``"""
    fy, fc, fs = y.reshape(-1), cc.reshape(-1), scale.reshape(-1)

    def piece(s: int, e: int) -> None:
        c = fc[s:e]
        with np.errstate(all="ignore"):
            np.multiply(fy[s:e], c, out=c)
            np.divide(c, fs[s:e], out=c)
        _max0(c)

    _parallel(piece, fc.size, threads)
    return cc


def _replicate(a: np.ndarray, s: int) -> np.ndarray:
    """each element as an s x ... x s block (MATLAB's estimate on the fine grid)"""
    out = a
    for ax in range(a.ndim):
        out = np.repeat(out, s, axis=ax)
    return np.ascontiguousarray(out)


def _block_mean(a: np.ndarray, s: int) -> np.ndarray:
    """the mean of each s x ... x s block (the fine grid back on the data's)"""
    if a.ndim == 2:
        ny, nx = a.shape
        v = a.reshape(ny // s, s, nx // s, s)
        return v.sum(axis=(1, 3)) / (s * s)
    nz, ny, nx = a.shape
    v = a.reshape(nz // s, s, ny // s, s, nx // s, s)
    return v.sum(axis=(1, 3, 5)) / (s * s * s)


# ---------------------------------------------------------------------------------------------
# deconvlucy
# ---------------------------------------------------------------------------------------------


def _lucy(state, psf, numit, dampar, weight, readout, subsmpl, backend, threads) -> LucyState:
    from petakit5d import fft
    from petakit5d.core.threads import resolve_threads
    from petakit5d.kernels import decon as K

    t = resolve_threads(threads)
    be = fft.backend() if backend is None else backend
    fresh = not isinstance(state, LucyState)
    if fresh:
        given = state
        image = _data(state)
        est = prev = g = None
    else:
        given = state.image
        image = _data(state.image)
        est = np.ascontiguousarray(_matlab_2d(np.asarray(state.estimate, dtype=np.float64)))
        prev = np.ascontiguousarray(_matlab_2d(np.asarray(state.previous, dtype=np.float64)))
        g = np.asarray(state.internal, dtype=np.float64).reshape((2, *est.shape))
    dt = image.dtype
    n = _iterations(numit)
    d = _same_class(dampar, dt, "DAMPAR", "images:deconvlucy:damparMustBeSameClassAsInputImage")
    r = _same_class(readout, dt, "READOUT", "images:deconvlucy:readoutMustBeSameClassAsInputImage")
    w = _weight(weight, image.shape)
    s = _subsampling(subsmpl)
    fine = tuple(v * s for v in image.shape)
    p = _psf(psf, fine)

    i_d = _to_double(image, dt)  # a new array
    if fresh:
        est = _replicate(i_d, s) if s > 1 else i_d
        prev = np.zeros(fine)
        g = np.zeros((2, *fine))
    elif est.shape != fine:
        raise DeconvlucyError(
            "images:deconvlucy:invalidCellArray", "The estimates of the cell don't fit the image."
        )
    # max(WEIGHT .* (READOUT + I), 0)
    wi = i_d + r if r else i_d.copy()
    if w is not None:
        np.multiply(w, wi, out=wi)
    _max0(wi)
    hs = fft.HalfSpectrum(fine, np.float64, backend=be, threads=t)
    h = _half_otf(hs, p)
    ones = np.ones(fine) if w is None else (_replicate(w, s) if s > 1 else w)
    scale = hs.convolve(ones, h, conjugate=True)
    scale += _SQRT_EPS
    d22 = d * d / 2

    j2, j3 = est, prev
    # buf: the change before the latest (MATLAB's J{4}(:, 2) after the call)
    buf = np.zeros(fine) if fresh else np.ascontiguousarray(g[1])
    y = None
    for k in range(n):
        if k == 0:
            # a fresh call has no changes yet (lambda 0); a resumed one has the last two
            lam = 0.0
            if not fresh:
                g0 = np.ascontiguousarray(g[0])
                lam = K.accel_lambda(g0.reshape(-1), np.zeros(g0.size), buf.reshape(-1)).lam
                buf = g0.copy()
        else:
            # lambda of the changes J - Y (written into buf) and the one before
            fb = buf.reshape(-1)
            lam = K.accel_lambda(j2.reshape(-1), y.reshape(-1), fb, threads=t, out=fb).lam
        y = K.accel_predict(j2, j3, lam, threads=t, out=y)
        rb = hs.convolve(y, h)
        if s > 1:
            rb = np.ascontiguousarray(_block_mean(rb, s))
        rb = _ratio(rb, wi, r, d22, t)
        if s > 1:
            rb = _replicate(rb, s)
        cc = hs.convolve(rb, h, conjugate=True)
        j3 = j2
        j2 = _update(y, cc, scale, t)
    internal = np.stack([j2 - y, buf]) if n else np.asarray(g)
    return LucyState(np.asarray(given), j2, j3, internal)


def deconvlucy_state(
    i: Any,
    psf: Any,
    numit: Any = 10,
    dampar: Any = None,
    weight: Any = None,
    readout: Any = None,
    subsmpl: Any = None,
    *,
    backend: str | None = None,
    threads: int | None = None,
) -> LucyState:
    """MATLAB's ``J = deconvlucy({I}, PSF, ...)`` (``i`` an array) or ``deconvlucy(J, PSF, ...)``
    (``i`` a ``LucyState``): the state after ``numit`` more iterations (see the module)."""
    return _lucy(i, psf, numit, dampar, weight, readout, subsmpl, backend, threads)


def deconvlucy(
    i: Any,
    psf: Any,
    numit: Any = 10,
    dampar: Any = None,
    weight: Any = None,
    readout: Any = None,
    subsmpl: Any = None,
    *,
    backend: str | None = None,
    threads: int | None = None,
) -> np.ndarray:
    """MATLAB's ``deconvlucy(I, PSF, NUMIT, DAMPAR, WEIGHT, READOUT, SUBSMPL)``: Richardson-Lucy
    deconvolution with acceleration and damping, the result in the class of ``i`` (see the
    module). None (or an empty array) is a parameter's default: NUMIT 10, DAMPAR 0, WEIGHT ones,
    READOUT 0, SUBSMPL 1. ``backend``: the CPU FFT backend (``petakit5d.fft``; None: the default);
    ``threads``: None for the CPUs this process may use (results don't depend on it)."""
    st = _lucy(i, psf, numit, dampar, weight, readout, subsmpl, backend, threads)
    return _from_double(st.estimate, np.asarray(st.image).dtype)


# ---------------------------------------------------------------------------------------------
# RLdecon's branch
# ---------------------------------------------------------------------------------------------


def rl_decon_original(
    rawdata: np.ndarray,
    psf: np.ndarray,
    n_iter: Any,
    *,
    decon_bbox: Sequence[Any] | None = None,
    scale_factor: Any = 1,
    backend: str | None = None,
    threads: int | None = None,
) -> np.ndarray:
    """RLdecon's ``'original'`` branch: ``deconvlucy(rawdata, psf, nIter)`` (the data already
    ``max(single(raw) - Background, 0)``, the PSF already normalized), cropped to ``decon_bbox``
    (a Python box) and times ``scale_factor`` where it isn't 1 (MATLAB's single times double)."""
    from petakit5d.core import dtypes
    from petakit5d.decon.rl import _crop

    out = deconvlucy(rawdata, psf, n_iter, backend=backend, threads=threads)
    if decon_bbox is not None and len(decon_bbox):
        out = _crop(out, decon_bbox)
    if scale_factor != 1:
        out = np.asarray(dtypes.times(out, scale_factor))
    return out
