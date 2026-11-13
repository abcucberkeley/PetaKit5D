"""The deconvolution kernels (``petakit5d._ext.decon``, moved from the mex files of
``microscopeDataProcessing/decon/mex``): the half-spectrum Richardson-Lucy steps of
decon_lucy_rfft_mex, the acceleration step and fixed-order sums of decon_lucy_accel_mex
(decon_lucy_accel, decon_sum_double), the mask erosion of decon_mask_edge_erosion_mex and the steps
of the OMW back projector (omw_backprojector_mex).

Arrays are the port's (``docs/conventions.md``): volumes ``(z, y, x)``, images ``(y, x)``,
C-contiguous. Every function takes ``threads`` (None: the CPUs this process may use); results do
not depend on it. Errors the mex raises carry its identifier (``DeconKernelError.identifier``).

**Richardson-Lucy steps** (``LucyRFFT``; ``otf``, ``omw``, ``rlstep``, ``fft1``), on the CPU FFT
backends of ``petakit5d.fft``, in one of two layouts:

- ``layout="matlab"`` (the default on the exact backend): the mex's loops in MATLAB's layout
  (MATLAB's ``(y, x, z)`` array, y fastest), with MATLAB's 1D transforms, groups of four lines,
  lines run alone, Hermitian fill and division after every inverse stage, on the backend's FFT
  table. On the exact backend the results are the mex's (MATLAB's ``fftn``/``ifftn`` code) bit
  for bit; ``matlab_parallel`` is MATLAB's FFT mode (its split of a stage's lines among threads,
  decon_lucy_rfft's ``grain``). OTFs are MATLAB's half OTFs: Fortran arrays
  ``(ny, nx, nz // 2 + 1)`` (2D ``(ny, nx // 2 + 1)``), exactly MATLAB's ``H``. Volumes in and
  out are ``(z, y, x)``: one transpose in and one out per call; ``to_internal``/``from_internal``
  and ``internal=True`` keep an iteration's arrays in MATLAB's layout across calls (RL's loop),
  and the acceleration functions take them with ``matlab_layout=True``.
- ``layout="native"`` (the default on the permissive backend): the native half spectra of
  ``petakit5d.fft.HalfSpectrum`` (pocketfft, halved along x: ``(..., nx // 2 + 1)``) with the
  mex's element-wise steps per voxel; within decision record 0006's tolerances of MATLAB. No
  transposes; the internal arrays are the ``(z, y, x)`` arrays themselves.

The element-wise arithmetic is the mex's: the background subtracted in double and rounded to
single, ``max(., 0)``, ``max(v, eps)`` with double ``eps``, MATLAB's complex product without fused
multiply-adds, the update ``max(v .* J, 0)`` (which makes -0 +0).

**Acceleration and sums**: ``accel_lambda`` and ``accel_predict`` (decon_lucy_accel's 'lambda' and
'predict'), ``sum_double`` (decon_sum_double): from 89,000 elements the mex's fixed order (blocks of
4,096 in MATLAB's linear order, four lanes each, the blocks in order), below MATLAB's own
(``petakit5d.kernels.ops.sum_double``). ``(z, y, x)`` arrays are read in MATLAB's linear order;
``matlab_layout=True`` takes MATLAB's arrays (Fortran-contiguous, or vectors).

**Mask erosion**: ``mask_edge_erosion`` (decon_mask_edge_erosion_mex): numeric neighborhoods are
converted as MATLAB's ``logical`` does (the mex called back into MATLAB for it).

**OMW back projector steps** (omw_backprojector_mex's modes): ``omw_convhull`` (with N03's
``bwconvhull_slices``, the mex's algorithm), ``omw_window`` (the window bounds), ``omw_components``
(``bwareaopen`` and the largest components), ``omw_sortdesc``, ``omw_medfilt3`` (N03's
``medfilt3``, the mex's selection network); the last two return None where the mex returns
empty (NaN or negative values; NaN or -0).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any, NamedTuple

import numpy as np

__all__ = [
    "SUM_BLOCKS_FROM",
    "AccelLambda",
    "DeconKernelError",
    "LucyRFFT",
    "accel_lambda",
    "accel_predict",
    "fft1",
    "from_matlab_layout",
    "mask_edge_erosion",
    "omw",
    "omw_components",
    "omw_convhull",
    "omw_medfilt3",
    "omw_sortdesc",
    "omw_window",
    "otf",
    "rlstep",
    "sum_double",
    "to_matlab_layout",
]

# from this many elements MATLAB's sums depend on its thread count; the mex sums them in blocks
SUM_BLOCKS_FROM = 89000
_DBL_EPS = 2.220446049250313e-16
_RFFT = "decon_lucy_rfft_mex"
_ACCEL = "decon_lucy_accel_mex:inputError"
_EROSION = "decon_mask_edge_erosion_mex:inputError"
_OMW = "omw_backprojector_mex:invalidInput"


class DeconKernelError(ValueError):
    """An error the mex raises; ``identifier`` is its MATLAB error identifier."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _ext():
    from petakit5d._ext import load

    return load("decon")


def _threads(n: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(n)


def _err(identifier: str, message: str) -> DeconKernelError:
    return DeconKernelError(identifier, message)


def _backend(name: str | None) -> tuple[str, int]:
    """the backend's name and its id; the exact one is loaded on first use into the process's FFT
    state, which this module shares with petakit5d._ext.fft"""
    from petakit5d import fft

    name = fft.backend() if name is None else name
    if name not in fft.BACKENDS:
        raise ValueError(f"unknown FFT backend {name!r}: use one of {fft.BACKENDS}")
    return name, fft._backend_id(name)


def _count(v: Any, name: str) -> int:
    """a count as a MATLAB loop 1 : v runs it (0 for v < 1 or NaN), the mex's scalar_int"""
    a = np.asarray(v)
    if a.size != 1 or a.dtype.kind not in "biuf":
        raise _err(f"{_RFFT}:inputError", f"{name} must be a real numeric scalar.")
    x = float(a.reshape(-1)[0])
    if not x >= 1:
        return 0
    if x > 2147483647.0:
        raise _err(f"{_RFFT}:inputError", f"{name} is too large.")
    return int(x)


# ---------------------------------------------------------------------------------------------
# layouts
# ---------------------------------------------------------------------------------------------


def to_matlab_layout(vol: np.ndarray, *, threads: int | None = None) -> np.ndarray:
    """MATLAB's array of a (z, y, x) volume or (y, x) image: a Fortran array (y, x, z) (or (y, x)),
    ``core.arrays.to_matlab`` by a threaded blocked transpose."""
    a = np.ascontiguousarray(vol)
    if a.ndim not in (2, 3):
        raise ValueError("to_matlab_layout: a (z, y, x) volume or a (y, x) image")
    nz, ny, nx = (1, *a.shape) if a.ndim == 2 else a.shape
    out = np.empty((nz, nx, ny), dtype=a.dtype)
    _ext().transpose_planes(a.reshape(-1), out.reshape(-1), ny, nx, _threads(threads))
    return out[0].T if a.ndim == 2 else np.transpose(out, (2, 1, 0))


def from_matlab_layout(a: np.ndarray, *, threads: int | None = None) -> np.ndarray:
    """The (z, y, x) volume (or (y, x) image) of MATLAB's array (y, x, z) (or (y, x)), Fortran-
    contiguous: ``core.arrays.from_matlab`` by a threaded blocked transpose."""
    a = np.asarray(a)
    if a.ndim not in (2, 3) or not a.flags.f_contiguous:
        raise ValueError("from_matlab_layout: a Fortran-contiguous (y, x, z) or (y, x) array")
    ny, nx = a.shape[0], a.shape[1]
    nz = a.shape[2] if a.ndim == 3 else 1
    out = np.empty((nz, ny, nx), dtype=a.dtype)
    _ext().transpose_planes(a.reshape(-1, order="F"), out.reshape(-1), nx, ny, _threads(threads))
    return out[0] if a.ndim == 2 else out


# ---------------------------------------------------------------------------------------------
# Richardson-Lucy steps (decon_lucy_rfft_mex)
# ---------------------------------------------------------------------------------------------


def _volume(a, name: str, dtypes: tuple, shape: tuple[int, ...] | None = None) -> np.ndarray:
    if not isinstance(a, np.ndarray) or np.iscomplexobj(a) or a.dtype not in dtypes:
        raise _err(f"{_RFFT}:inputError", f"{name} must be a real array of class {dtypes}.")
    if shape is not None and a.shape != shape:
        raise _err(
            f"{_RFFT}:inputError", f"{name} must have the data's shape {shape}, not {a.shape}."
        )
    return np.ascontiguousarray(a)


class LucyRFFT:
    """The half-spectrum Richardson-Lucy steps of decon_lucy_rfft_mex for data of one shape.

    shape: (nz, ny, nx) or (ny, nx), ny and nx at least 2 (the mex's rule). backend: the CPU FFT
    backend (``petakit5d.fft``; None: its default). layout: 'matlab' or 'native' (None: 'matlab'
    on the exact backend, 'native' on the permissive one; see the module docstring).
    matlab_parallel: MATLAB's FFT mode for the 'matlab' layout (the parallel split of its
    multithreaded sessions, or a session started on one thread). ``last_times``: the seconds of
    the last call in planning, forward transforms, inverse transforms and other work.
    """

    def __init__(
        self,
        shape: Sequence[int],
        *,
        backend: str | None = None,
        threads: int | None = None,
        matlab_parallel: bool = True,
        layout: str | None = None,
    ):
        self.shape = tuple(int(v) for v in shape)
        if len(self.shape) not in (2, 3) or any(v < 1 for v in self.shape):
            raise _err(f"{_RFFT}:inputError", "Only 2D and 3D data are supported.")
        if self.shape[-2] < 2 or self.shape[-1] < 2:
            raise _err(
                f"{_RFFT}:inputError", "The data must have at least 2 elements along y and x."
            )
        if any(v > 2147483647 for v in self.shape):
            raise _err(f"{_RFFT}:inputError", "Dimensions must be smaller than 2^31.")
        self.backend, self._be = _backend(backend)
        self.layout = layout or ("matlab" if self.backend == "exact" else "native")
        if self.layout not in ("matlab", "native"):
            raise ValueError("layout must be 'matlab' or 'native'")
        if self.layout == "native" and self.backend != "permissive":
            raise ValueError("the native layout runs on the permissive backend's half spectra")
        self.threads = threads
        self.matlab_parallel = bool(matlab_parallel)
        self._grain = int(_ext().GRAIN_PARALLEL) if matlab_parallel else 0
        nz = self.shape[0] if len(self.shape) == 3 else 1
        ny, nx = self.shape[-2], self.shape[-1]
        self._d = [ny, nx, nz]
        if self.layout == "matlab":
            self.spectrum_shape = (
                (ny, nx, nz // 2 + 1) if len(self.shape) == 3 else (ny, nx // 2 + 1)
            )
            self.internal_shape = (ny, nx, nz) if len(self.shape) == 3 else (ny, nx)
        else:
            self.spectrum_shape = (*self.shape[:-1], nx // 2 + 1)
            self.internal_shape = self.shape
        self.last_times = (0.0, 0.0, 0.0, 0.0)

    # -- layouts ---------------------------------------------------------------------------------

    @property
    def matlab_layout(self) -> bool:
        """whether the internal arrays are in MATLAB's layout (for the acceleration functions)"""
        return self.layout == "matlab"

    def to_internal(self, vol) -> np.ndarray:
        """a (z, y, x) array as this object's internal array (MATLAB's Fortran array, or itself)"""
        a = np.asarray(vol)
        if a.shape != self.shape:
            raise ValueError(f"expected an array of shape {self.shape}, got {a.shape}")
        if self.layout == "native":
            return np.ascontiguousarray(a)
        return to_matlab_layout(a, threads=self.threads)

    def from_internal(self, a) -> np.ndarray:
        """the (z, y, x) array of an internal array"""
        a = np.asarray(a)
        if a.shape != self.internal_shape:
            raise ValueError(f"expected an internal array of shape {self.internal_shape}")
        if self.layout == "native":
            return np.ascontiguousarray(a)
        if not a.flags.f_contiguous:
            from petakit5d.core.arrays import from_matlab

            return from_matlab(a)
        return from_matlab_layout(a, threads=self.threads)

    def _flat_internal(self, a, name: str) -> np.ndarray:
        if not isinstance(a, np.ndarray) or a.dtype != np.float32 or a.shape != self.internal_shape:
            raise _err(
                f"{_RFFT}:inputError",
                f"{name} must be a single internal array of shape {self.internal_shape}.",
            )
        if self.layout == "native":
            return np.ascontiguousarray(a)
        if not a.flags.f_contiguous:
            raise ValueError(f"{name}: an internal array in MATLAB's layout is Fortran-contiguous")
        return a.reshape(-1, order="F")

    def _spectrum(self, h, name: str = "H") -> np.ndarray:
        if not isinstance(h, np.ndarray) or h.dtype != np.complex64:
            raise _err(f"{_RFFT}:inputError", "The OTFs must be complex single half OTFs.")
        if h.shape != self.spectrum_shape:
            raise _err(
                f"{_RFFT}:inputError", "The OTF size must be that of the half OTF of the data size."
            )
        if self.layout == "native":
            if not h.flags.c_contiguous:
                raise ValueError(f"{name}: a native half spectrum is C-contiguous")
            return h
        if not h.flags.f_contiguous:
            raise ValueError(f"{name}: MATLAB's half OTF is a Fortran array")
        return h.reshape(-1, order="F")

    def _out(self, flat: np.ndarray, internal: bool) -> np.ndarray:
        if self.layout == "native" or not internal:
            return flat.reshape(self.shape)
        return flat.reshape(self.internal_shape, order="F")

    # -- steps -----------------------------------------------------------------------------------

    def otf(self, psf) -> np.ndarray:
        """The half OTF of a (normalized) PSF: decon_psf2otf's OTF (the PSF zero-padded to the
        data's shape, shifted so that its center floor(size / 2) is at the origin, transformed),
        its half in this object's layout. The PSF has the data's number of axes and fits in it."""
        p = np.asarray(psf)
        if p.dtype != np.float32 or np.iscomplexobj(p) or p.size == 0:
            raise _err(f"{_RFFT}:inputError", "The PSF must be a real, non-empty single array.")
        if p.ndim > len(self.shape):
            raise _err(f"{_RFFT}:inputError", "Only 2D and 3D data are supported.")
        p = np.ascontiguousarray(p.reshape((1,) * (len(self.shape) - p.ndim) + p.shape))
        if any(a > b for a, b in zip(p.shape, self.shape, strict=True)):
            raise _err(f"{_RFFT}:outSizeIsSmallerThanPsfSize", "outSizeIsSmallerThanPsfSize")
        m, t = _ext(), _threads(self.threads)
        if self.layout == "native":
            h, tm = m.native_otf(p, list(self.shape), t)
            self.last_times = tuple(tm)
            return h
        p3 = p if p.ndim == 3 else p.reshape((1, *p.shape))
        h, tm = m.matlab_otf(p3, self._d, t, self._grain, self._be)
        self.last_times = tuple(tm)
        return h.reshape(self.spectrum_shape, order="F")

    def omw(self, raw, h_f, h_b, n_iter, background, j0=None) -> np.ndarray:
        """``n_iter`` OMW iterations from ``j0`` (default ``max(raw - background, 0)``):
        ``J <- max(real(ifftn(fftn(I ./ max(real(ifftn(fftn(J) .* H_f)), eps)) .* H_b)) .* J, 0)``
        with ``I = max(raw - background, 0)``. raw: uint16 or single (z, y, x); j0: single; the
        result: single (z, y, x)."""
        raw = _volume(raw, "The data", (np.dtype(np.uint16), np.dtype(np.float32)), self.shape)
        hf, hb = self._spectrum(h_f, "H_f"), self._spectrum(h_b, "H_b")
        n = _count(n_iter, "nIter")
        bg = np.asarray(background)
        if bg.size != 1 or np.iscomplexobj(bg) or bg.dtype.kind not in "biuf":
            raise _err(f"{_RFFT}:inputError", "background must be a real numeric scalar.")
        bgv = float(bg.reshape(-1)[0])
        if j0 is not None:
            if not isinstance(j0, np.ndarray) or j0.dtype != np.float32 or j0.shape != self.shape:
                raise _err(
                    f"{_RFFT}:inputError", "J0 must be a real single array of the data size."
                )
            j0 = np.ascontiguousarray(j0)
        m, t = _ext(), _threads(self.threads)
        kind = "uint16" if raw.dtype == np.uint16 else "single"
        if self.layout == "native":
            r, tm = getattr(m, f"native_omw_{kind}")(raw, bgv, hf, hb, n, j0, t)
        else:
            r, tm = getattr(m, f"matlab_omw_{kind}")(
                raw.reshape(-1),
                bgv,
                hf,
                hb,
                n,
                None if j0 is None else j0.reshape(-1),
                self._d,
                t,
                self._grain,
                self._be,
                True,
            )
        self.last_times = tuple(tm)
        return r.reshape(self.shape)

    def rlstep(self, y, i, h, *, internal: bool = False) -> np.ndarray:
        """One update of decon_lucy_function:
        ``max(Y .* real(ifftn(conj(H) .* fftn(I ./ max(real(ifftn(H .* fftn(Y))), eps)))), 0)``
        (the ratio's numerator is ``max(I, 0)``, as in the mex). Y, I: single; (z, y, x), or with
        ``internal=True`` this object's internal arrays, which the result then is too."""
        m, t = _ext(), _threads(self.threads)
        hs = self._spectrum(h)
        if internal:
            yv, iv = self._flat_internal(y, "Y"), self._flat_internal(i, "I")
        else:
            yv = _volume(y, "Y", (np.dtype(np.float32),), self.shape)
            iv = _volume(i, "I", (np.dtype(np.float32),), self.shape)
        if self.layout == "native":
            r, tm = m.native_rlstep(yv, iv, hs, t)
        else:
            if not internal:
                yv, iv = yv.reshape(-1), iv.reshape(-1)
            r, tm = m.matlab_rlstep(yv, iv, hs, self._d, t, self._grain, self._be, not internal)
        self.last_times = tuple(tm)
        return self._out(r.reshape(-1), internal)


def otf(psf, shape: Sequence[int], **kw) -> np.ndarray:
    """``decon_lucy_rfft('otf', psf, outSize)``: the half OTF of ``psf`` for data of ``shape``
    (``LucyRFFT`` keywords)."""
    return LucyRFFT(shape, **kw).otf(psf)


def omw(raw, h_f, h_b, n_iter, background, j0=None, **kw) -> np.ndarray:
    """``decon_lucy_rfft('omw', I, H_f, H_b, nIter, background[, J0])`` (``LucyRFFT.omw``)."""
    return LucyRFFT(np.shape(raw), **kw).omw(raw, h_f, h_b, n_iter, background, j0)


def rlstep(y, i, h, **kw) -> np.ndarray:
    """``decon_lucy_rfft('rlstep', Y, I, H)`` (``LucyRFFT.rlstep``)."""
    return LucyRFFT(np.shape(y), **kw).rlstep(y, i, h)


def fft1(
    x,
    *,
    backend: str | None = None,
    threads: int | None = None,
    matlab_parallel: bool = True,
) -> np.ndarray:
    """``decon_lucy_rfft_mex('fft1', X, ...)``: rows 0 .. n // 2 of MATLAB's ``fft(X, [], 1)`` of a
    single (n, m) matrix (n, m >= 2), with MATLAB's groups of lines in the given FFT mode (the
    lines are X's columns, numbered in order)."""
    a = np.asarray(x)
    if a.dtype != np.float32 or a.ndim != 2 or a.shape[0] < 2 or a.shape[1] < 2:
        raise _err(
            f"{_RFFT}:inputError",
            "X must be a real single matrix with at least 2 rows and columns.",
        )
    _, be = _backend(backend)
    n, mcols = a.shape
    grain = int(_ext().GRAIN_PARALLEL) if matlab_parallel else 0
    cols = np.ascontiguousarray(a.T)  # the columns contiguous: MATLAB's layout
    y = _ext().matlab_fft1(cols.reshape(-1), n, mcols, _threads(threads), grain, be)
    return np.ascontiguousarray(y.reshape((mcols, n // 2 + 1)).T)


# ---------------------------------------------------------------------------------------------
# acceleration and sums (decon_lucy_accel_mex, decon_lucy_accel, decon_sum_double)
# ---------------------------------------------------------------------------------------------


class AccelLambda(NamedTuple):
    """decon_lucy_accel('lambda'): the clamped factor, the new J_4 (J_2 - Y, in J_2's shape and
    layout) and the two sums."""

    lam: float
    j4: np.ndarray
    sums: tuple[float, float]


def _accel_arrays(arrays: Sequence[Any], names: Sequence[str], matlab_layout: bool):
    """the arrays as flat views in MATLAB's linear order of one class; (ny, nx) of (z, y, x) arrays
    (None: flat already)"""
    a0 = arrays[0]
    if (
        not isinstance(a0, np.ndarray)
        or a0.dtype not in (np.float32, np.float64)
        or np.iscomplexobj(a0)
    ):
        raise _err(_ACCEL, f"{names[0]} must be a real single or double array.")
    for a, nm in zip(arrays[1:], names[1:], strict=True):
        if not isinstance(a, np.ndarray) or a.dtype not in (np.float32, np.float64):
            raise _err(_ACCEL, f"{nm} must be a real single or double array.")
        if a.dtype != a0.dtype or a.size != a0.size:
            raise _err(_ACCEL, f"{nm} must have the class and number of elements of J_2.")
    flats = []
    yx = None
    if matlab_layout or a0.ndim < 2:
        for a, nm in zip(arrays, names, strict=True):
            if a.ndim >= 2 and not a.flags.f_contiguous:
                raise ValueError(f"{nm}: an array in MATLAB's layout is Fortran-contiguous")
            flats.append(a.reshape(-1, order="F"))
    else:
        for a, nm in zip(arrays, names, strict=True):
            if a.shape != a0.shape:
                raise ValueError(f"{nm}: the (z, y, x) arrays have J_2's shape")
            flats.append(np.ascontiguousarray(a).reshape(-1))
        yx = (a0.shape[-2], a0.shape[-1])
    return flats, yx


def _clamp(v: float) -> float:
    """max(min(v, 1), 0) as MATLAB: min(NaN, 1) is 1, -0 is kept"""
    m = 1.0 if math.isnan(v) else (v if v < 1.0 else 1.0)
    return m if m >= 0.0 else 0.0


def accel_lambda(
    j2, y, j4, *, matlab_layout: bool = False, threads: int | None = None, out=None
) -> AccelLambda:
    """decon_lucy_accel('lambda', J_2, Y, J_4): ``lambda = sum((J_2 - Y) .* J_4) / (sum(J_4 .* J_4)
    + eps)`` with the products in the data's class and the sums in double (``sum_double``'s
    order), clamped to [0, 1] as MATLAB's max(min(lambda, 1), 0), and ``J_4 = J_2 - Y``. The
    arrays: real single or double, of one class and number of elements; (z, y, x) arrays of one
    shape read in MATLAB's linear order, or MATLAB's arrays with ``matlab_layout=True``. ``out``:
    an array of J_2's class, shape and layout for the new J_4 (it may be J_4 itself)."""
    arrays = [j2, y, j4]
    (f2, fy, f4), yx = _accel_arrays(arrays, ["J_2", "Y", "J_4"], matlab_layout)
    n = f2.size
    kind = "single" if f2.dtype == np.float32 else "double"

    def shaped(flat):
        if yx is None:
            return flat.reshape(j2.shape, order="F") if j2.ndim >= 2 else flat
        return flat.reshape(j2.shape)

    if n >= SUM_BLOCKS_FROM:
        ny, nx = yx if yx is not None else (1, 1)
        fo = None if out is None else _out_flat(out, j2, yx is None)
        lam, j4n, s0, s1 = getattr(_ext(), f"accel_lambda_{kind}")(
            f2, fy, f4, yx is not None, ny, nx, _threads(threads), fo
        )
        return AccelLambda(
            float(lam), out if out is not None else shaped(j4n), (float(s0), float(s1))
        )
    # MATLAB's code: the products in the class, MATLAB's sums
    with np.errstate(all="ignore"):
        d = f2 - fy
        p1, p2 = d * f4, f4 * f4
    s0 = sum_double(_matlab_vector(p1, yx, j2.shape), matlab_layout=True, threads=threads)
    s1 = sum_double(_matlab_vector(p2, yx, j2.shape), matlab_layout=True, threads=threads)
    if out is not None:
        _out_flat(out, j2, yx is None)[...] = d
        return AccelLambda(_clamp(s0 / (s1 + _DBL_EPS)), out, (s0, s1))
    return AccelLambda(_clamp(s0 / (s1 + _DBL_EPS)), shaped(d), (s0, s1))


def _out_flat(out, like: np.ndarray, matlab_layout: bool) -> np.ndarray:
    """a flat view of an output array of like's class, shape and layout"""
    if not isinstance(out, np.ndarray) or out.dtype != like.dtype or out.shape != like.shape:
        raise ValueError("out must have the class and shape of J_2")
    if matlab_layout and like.ndim >= 2:
        if not out.flags.f_contiguous:
            raise ValueError("out: an array in MATLAB's layout is Fortran-contiguous")
        return out.reshape(-1, order="F")
    if not out.flags.c_contiguous:
        raise ValueError("out must be C-contiguous")
    return out.reshape(-1)


def _matlab_vector(flat: np.ndarray, yx, shape) -> np.ndarray:
    """values of a flat array in MATLAB's linear order (flat already, or of a (z, y, x) array)"""
    if yx is None:
        return flat
    from petakit5d.core.arrays import matlab_order

    return np.ascontiguousarray(matlab_order(flat.reshape(shape))).reshape(-1)


def accel_predict(j2, j3, lam: float, *, threads: int | None = None, out=None) -> np.ndarray:
    """decon_lucy_accel('predict', J_2, J_3, lambda): ``max(J_2 + lambda * (J_2 - J_3), 0)``, the
    product in double rounded to the class, NaN to 0, -0 kept; in J_2's shape and layout (the
    arithmetic is per element, so any layout). ``out``: an array of J_2's class, shape and layout
    for the result (it may be J_2 or J_3)."""
    for a, nm in ((j2, "J_2"), (j3, "J_3")):
        if not isinstance(a, np.ndarray) or a.dtype not in (np.float32, np.float64):
            raise _err(_ACCEL, f"{nm} must be a real single or double array.")
    if j3.dtype != j2.dtype or j3.size != j2.size:
        raise _err(_ACCEL, "J_3 must have the class and number of elements of J_2.")
    lv = np.asarray(lam)
    if lv.size != 1 or lv.dtype != np.float64 and not isinstance(lam, int | float):
        raise _err(_ACCEL, "lambda must be a real double scalar.")
    # MATLAB's layout (Fortran arrays) or the port's: the result in J_2's
    matlab_layout = j2.ndim >= 2 and j2.flags.f_contiguous and not j2.flags.c_contiguous
    if matlab_layout and not j3.flags.f_contiguous:
        raise ValueError("J_3: J_2 is in MATLAB's layout, so is J_3")
    if not matlab_layout and j3.shape != j2.shape and j3.ndim >= 2:
        raise ValueError("J_3: the (z, y, x) arrays have J_2's shape")
    if matlab_layout:
        f2, f3 = j2.reshape(-1, order="F"), j3.reshape(-1, order="F")
    else:
        f2, f3 = np.ascontiguousarray(j2).reshape(-1), np.ascontiguousarray(j3).reshape(-1)
    kind = "single" if j2.dtype == np.float32 else "double"
    fo = None if out is None else _out_flat(out, j2, matlab_layout)
    y = getattr(_ext(), f"accel_predict_{kind}")(
        f2, f3, float(lv.reshape(-1)[0]), _threads(threads), fo
    )
    if out is not None:
        return out
    return y.reshape(j2.shape, order="F") if matlab_layout else y.reshape(j2.shape)


def sum_double(x, *, matlab_layout: bool = False, threads: int | None = None) -> float:
    """decon_sum_double(x): ``sum(x(:), 'double')`` in an order that doesn't depend on the number
    of threads. Below 89,000 elements (and for classes other than single and double) MATLAB's own
    sum (``petakit5d.kernels.ops.sum_double``); from there decon_lucy_accel_mex's: blocks of 4,096
    elements in MATLAB's linear order, each in four lanes, the block sums in order (a NaN sum with
    the mex's payload). (z, y, x) arrays are read in MATLAB's linear order; MATLAB's arrays with
    ``matlab_layout=True``."""
    from petakit5d.kernels import ops

    a = np.asarray(x)
    if a.size < SUM_BLOCKS_FROM or a.dtype not in (np.float32, np.float64):
        return ops.sum_double(a, matlab_layout=matlab_layout, num_threads=threads)
    kind = "single" if a.dtype == np.float32 else "double"
    if matlab_layout or a.ndim < 2:
        # MATLAB's linear order of MATLAB's array (a view when it is Fortran-contiguous)
        v = np.ascontiguousarray(a.reshape(-1, order="F"))
        return float(getattr(_ext(), f"accel_sum_{kind}")(v, False, 1, 1, _threads(threads)))
    a = np.ascontiguousarray(a)
    return float(
        getattr(_ext(), f"accel_sum_{kind}")(
            a.reshape(-1), True, a.shape[-2], a.shape[-1], _threads(threads)
        )
    )


# ---------------------------------------------------------------------------------------------
# mask erosion (decon_mask_edge_erosion_mex)
# ---------------------------------------------------------------------------------------------


def _nonneg_int(v, name: str, kind: str = "scalar") -> int:
    x = float(v)
    if not (x >= 0) or x > 1e9 or x != math.floor(x):
        if kind == "pad":
            raise _err(_EROSION, f"{name} must contain nonnegative integers.")
        raise _err(_EROSION, f"{name} must be a nonnegative integer.")
    return int(x)


def _pad(v, name: str) -> list[int]:
    a = np.asarray(v)
    if a.dtype.kind not in "iuf" or a.size != 2 or np.iscomplexobj(a):
        raise _err(_EROSION, f"{name} must be a double vector of 2 elements.")
    return [_nonneg_int(t, name, "pad") for t in a.reshape(-1).tolist()]


def _scalar(v, name: str) -> float:
    a = np.asarray(v)
    if a.size != 1 or np.iscomplexobj(a) or a.dtype.kind not in "biuf":
        raise _err(_EROSION, f"{name} must be a real numeric scalar.")
    return float(a.reshape(-1)[0])


def _nhood_offsets(nhoods) -> tuple[list[int], list[int], list[int]]:
    """the offsets (row = y, column = x) of the nonzero elements of each neighborhood about its
    center floor((size + 1) / 2) (1-based), in MATLAB's column-major order"""
    if isinstance(nhoods, np.ndarray) or not isinstance(nhoods, list | tuple):
        raise _err(_EROSION, "nhoods must be a cell array.")
    start, dr, dc = [0], [], []
    for e in nhoods:
        a = np.asarray(e)
        if a.ndim > 2 or a.size == 0 or np.iscomplexobj(a) or a.dtype.kind not in "biuf":
            raise _err(_EROSION, "Each neighborhood must be a real 2D array.")
        a = a.reshape((1,) * (2 - a.ndim) + a.shape)
        if a.dtype.kind == "f" and np.isnan(a).any():
            # the mex's call of MATLAB's logical() raises MATLAB's own error
            raise _err("MATLAB:nologicalnan", "NaN values cannot be converted to logicals.")
        m, n = a.shape
        cy, cx = (m + 1) // 2 - 1, (n + 1) // 2 - 1
        cols, rows = np.nonzero((a != 0).T)  # column-major: by column, then row
        if rows.size == 0:
            raise _err(_EROSION, "Each neighborhood must have a nonzero element.")
        dr.extend(int(r - cy) for r in rows.tolist())
        dc.extend(int(c - cx) for c in cols.tolist())
        start.append(len(dr))
    return start, dr, dc


def mask_edge_erosion(
    mask,
    nhoods: Sequence[Any],
    pad_ul: Sequence[int],
    pad_lr: Sequence[int],
    rz: int,
    clear_border: bool,
    *,
    threads: int | None = None,
) -> np.ndarray:
    """decon_mask_edge_erosion_mex(A, nhoods, padUL, padLR, rz, clearBorder): each z slice of the
    bool mask ((z, y, x) or (y, x)) eroded by the neighborhoods in turn as imerode applies a
    decomposed structuring element (2D arrays (rows y, columns x), logical or numeric; numeric
    nonzero is true), on the slice padded with true by pad_ul and pad_lr ((rows, columns)), then
    along z by 2 * rz + 1 voxels (3D); clear_border sets the first and last rows and columns (3D:
    and slices) false first. Pixels outside count as true."""
    a = np.asarray(mask)
    if a.dtype != np.bool_ or a.ndim > 3:
        raise _err(_EROSION, "The mask must be a full 2D or 3D logical array.")
    start, dr, dc = _nhood_offsets(nhoods)
    pu, pl = _pad(pad_ul, "padUL"), _pad(pad_lr, "padLR")
    rzv = _scalar(rz, "rz")
    if not (rzv >= 0) or rzv > 1e9 or rzv != math.floor(rzv):
        raise _err(_EROSION, "rz must be a nonnegative integer.")
    clear = _scalar(clear_border, "clearBorder") != 0
    v = np.ascontiguousarray(a.reshape((1,) * (3 - a.ndim) + a.shape))
    out = _ext().mask_edge_erosion(
        v.view(np.uint8), start, dr, dc, pu, pl, int(rzv), clear, _threads(threads)
    )
    return out.reshape(a.shape)


# ---------------------------------------------------------------------------------------------
# OMW back projector steps (omw_backprojector_mex)
# ---------------------------------------------------------------------------------------------


def _slice_indices(v, n: int, name: str) -> np.ndarray:
    a = np.asarray(v if v is not None else [], dtype=np.float64).reshape(-1)
    if a.size and not np.all((a >= 0) & (a <= n - 1) & (a == np.floor(a))):
        raise _err(_OMW, f"{name} must be indices from 0 to {n - 1}.")
    idx = a.astype(np.int64)
    if np.unique(idx).size != idx.size:
        raise _err(_OMW, f"{name} must not repeat indices.")
    return idx


def omw_convhull(mask, zinds, yinds, xinds, *, threads: int | None = None) -> np.ndarray:
    """omw_backprojector_mex('convhull'): the (z, y, x) bool mask (every size above 1) with its
    z slices ``vol[z]`` at zinds, then the y slices (MATLAB's ``squeeze(A(y, :, :))``, rows x,
    columns z) at yinds, then the x slices (``squeeze(A(:, x, :))``, rows y, columns z) at xinds
    (0-based, each slice once) replaced by their ``bwconvhull(slice, 'union')``."""
    from petakit5d.matlabcompat.morphology import bwconvhull_slices

    a = np.asarray(mask)
    if a.dtype != np.bool_ or a.ndim != 3 or min(a.shape) < 2:
        raise _err(_OMW, "A must be a 3D logical array with all dimensions > 1.")
    zi = _slice_indices(zinds, a.shape[0], "zinds")
    yi = _slice_indices(yinds, a.shape[1], "yinds")
    xi = _slice_indices(xinds, a.shape[2], "xinds")
    out = np.ascontiguousarray(a)
    for axis, inds in ((0, zi), (1, yi), (2, xi)):
        if inds.size:
            out = bwconvhull_slices(out, axis, inds, threads=threads)
    return out.copy() if out is a else out


def omw_window(shape: Sequence[int], e, r, bounds, *, threads: int | None = None):
    """omw_backprojector_mex('window'): where the apodization window of a grid of ``shape``
    ((nz, ny, nx), center (shape + 1) / 2 in MATLAB's 1-based coordinates) is 1 or 0 whatever the
    nearest of the edge points are (unit vectors ``e``, (n, 3) with columns (z, y, x), n >= 2,
    radii ``r``), for ``bounds = (l, u)``, 0 < l < u < inf. Returns the (z, y, x) float64 mask (1
    where the window is 1) and the (z, y, x) indices of the voxels whose window has to be computed,
    in MATLAB's linear order."""
    from petakit5d.core.arrays import matlab_index_to_zyx

    sz = [float(v) for v in np.asarray(shape, dtype=np.float64).reshape(-1)]
    if len(sz) != 3 or not all(v >= 1 and v == math.floor(v) for v in sz):
        raise _err(_OMW, "psz must be positive integers.")
    ev = np.asarray(e)
    rv = np.asarray(r)
    if (
        ev.dtype != np.float64
        or ev.ndim != 2
        or ev.shape[1] != 3
        or ev.shape[0] < 2
        or rv.dtype != np.float64
        or rv.size != ev.shape[0]
    ):
        raise _err(_OMW, "E must be an n x 3 double array (n >= 2) and R n doubles.")
    if not np.all(np.abs(ev) <= 1.5):
        raise _err(_OMW, "E must be unit vectors.")
    b = np.asarray(bounds, dtype=np.float64).reshape(-1)
    if b.size != 2:
        raise _err(_OMW, "bounds must be [l, u].")
    lo, hi = float(b[0]), float(b[1])
    if not (lo > 0 and lo < hi and hi < math.inf):
        raise _err(_OMW, "bounds must be [l, u] with 0 < l < u < Inf.")
    nz, ny, nx = (int(v) for v in sz)
    # the mex's columns (y, x, z), column-major
    cols = np.ascontiguousarray(ev[:, [1, 2, 0]].T).reshape(-1)
    mask, inds = _ext().omw_window(
        [ny, nx, nz], cols, np.ascontiguousarray(rv.reshape(-1)), lo, hi, _threads(threads)
    )
    return mask, matlab_index_to_zyx(inds, (nz, ny, nx))


def omw_components(mask, mode: str, n, *, threads: int | None = None) -> np.ndarray:
    """omw_backprojector_mex('components'): the 26-connected components of the (z, y, x) bool mask
    that ``bwareaopen(A, n)`` keeps (mode 'minsize': at least n voxels), or the n largest ones
    ('largest'; ties go to the component whose first voxel comes first in MATLAB's order)."""
    a = np.asarray(mask)
    if a.dtype != np.bool_ or a.ndim != 3:
        raise _err(_OMW, "A must be a 3D logical array.")
    if a.size >= 2147483647:
        raise _err(_OMW, "A must have fewer than 2^31 - 1 elements.")
    if mode not in ("minsize", "largest"):
        raise _err(_OMW, "The mode must be 'minsize' or 'largest'.")
    nv = np.asarray(n)
    if nv.size != 1 or nv.dtype.kind not in "iuf":
        raise _err(_OMW, "n must be a nonnegative integer.")
    x = float(nv.reshape(-1)[0])
    if not (x >= 0) or x != math.floor(x) or x > 9e15:
        raise _err(_OMW, "n must be a nonnegative integer.")
    v = np.ascontiguousarray(a)
    return _ext().omw_components(v.view(np.uint8), mode == "largest", int(x), _threads(threads))


def omw_sortdesc(a, *, threads: int | None = None) -> np.ndarray | None:
    """omw_backprojector_mex('sortdesc'): ``sort(double(a(:)), 'descend')`` of single values, or
    None when a has NaN or negative values (-0 included), which the mex leaves to MATLAB's sort.
    The sorted values don't depend on a's order or layout."""
    v = np.asarray(a)
    if v.dtype != np.float32:
        raise _err(
            _OMW, "Usage: v = omw_backprojector_mex('sortdesc', a, nThreads) with real single a"
        )
    return _ext().omw_sortdesc(np.ascontiguousarray(v).reshape(-1), _threads(threads))


def omw_medfilt3(a, *, threads: int | None = None) -> np.ndarray | None:
    """omw_backprojector_mex('medfilt3'): ``medfilt3(A, [3, 3, 3])`` with symmetric padding of a
    single (z, y, x) volume (or (y, x) image), or None when A has NaN or -0 values, whose order
    the mex's network doesn't keep as medfilt3 does (the mex returns empty)."""
    from petakit5d.matlabcompat.morphology import medfilt3

    v = np.asarray(a)
    if v.dtype != np.float32 or v.ndim > 3:
        raise _err(
            _OMW, "Usage: B = omw_backprojector_mex('medfilt3', A, nThreads) with real single 3D A"
        )
    if v.size == 0:
        return np.zeros((0, 0))
    if np.isnan(v).any() or np.any((v == 0) & np.signbit(v)):
        return None
    v = v.reshape((1,) * (2 - v.ndim) + v.shape) if v.ndim < 2 else v
    return medfilt3(np.ascontiguousarray(v), (3, 3, 3), "symmetric", threads=threads)
