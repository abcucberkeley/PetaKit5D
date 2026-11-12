"""The CPU FFT backends and the half-spectrum transforms FFT-based kernels convolve with.

Two backends run MATLAB's fft driver (``petakit5d.matlabcompat.fft``) and the half spectra below:

- ``"permissive"`` (the default): pocketfft (BSD-3), built in. Results differ from MATLAB's in the
  last bits (decision record 0006: tolerances and measurements); they do not depend on the thread
  count, the SIMD width or the platform's vector unit.
- ``"exact"``: FFTW 3.3.8 with MATLAB R2024b's codelet sets (GPL-2.0-or-later), an optional add-on
  of two shared libraries (the ``petakit5d-fftw`` package, ``python/cpp/fft/exact``). With it the
  results are MATLAB's bit for bit on Linux x86-64. petakit5d contains no GPL code; the add-on is
  found through ``PETAKIT5D_FFTW`` (its folder, or the single and double precision libraries
  separated by ``os.pathsep``), the installed ``petakit5d_fftw`` package, or ``PETAKIT5D_TEST_FFTW``
  (MATLAB's ``bin/glnxa64`` folder, whose libmwfftw3f / libmwfftw3 are the same library).

``backend()`` gives the default backend of the calls that name none, ``set_backend(name)`` sets it;
``PETAKIT5D_FFT`` (``permissive``, ``exact`` or ``auto``: exact when the add-on is there) sets it at
import. Both backends can be used in one process: every call may name its own.
"""

from __future__ import annotations

import os
import sys
from collections.abc import Iterator, Sequence
from contextlib import contextmanager

import numpy as np

__all__ = [
    "BACKENDS",
    "HalfSpectrum",
    "available",
    "backend",
    "exact_libraries",
    "info",
    "set_backend",
    "using",
]

BACKENDS = ("permissive", "exact")
_MATLAB_LIBS = ("libmwfftw3f.so.3", "libmwfftw3.so.3")
_state: dict = {"ext": None, "default": None, "exact": None}


def _ext():
    if _state["ext"] is None:
        from petakit5d._ext import load

        _state["ext"] = load("fft")
    return _state["ext"]


def _pair(spec: str, names: tuple[str, str]) -> tuple[str, str] | None:
    if os.path.isdir(spec):
        paths = tuple(os.path.join(spec, n) for n in names)
        return paths if all(os.path.exists(p) for p in paths) else None
    parts = spec.split(os.pathsep)
    if len(parts) != 2:
        raise ValueError(
            f"an FFTW location must be a folder or two library paths separated by {os.pathsep!r}"
        )
    return parts[0], parts[1]


def _addon_names() -> tuple[str, str]:
    if sys.platform == "win32":
        return ("pk5d_fftw3f.dll", "pk5d_fftw3.dll")
    if sys.platform == "darwin":
        return ("libpk5d_fftw3f.dylib", "libpk5d_fftw3.dylib")
    return ("libpk5d_fftw3f.so", "libpk5d_fftw3.so")


def exact_libraries() -> tuple[str, str] | None:
    """The exact backend's single and double precision libraries, or None: ``PETAKIT5D_FFTW``,
    the ``petakit5d_fftw`` package, ``PETAKIT5D_TEST_FFTW`` (MATLAB's own libraries)."""
    spec = os.environ.get("PETAKIT5D_FFTW", "")
    if spec:
        found = _pair(spec, _addon_names())
        if found is None:
            raise FileNotFoundError(f"PETAKIT5D_FFTW={spec}: no FFTW libraries there")
        return found
    try:
        import petakit5d_fftw

        paths = petakit5d_fftw.library_paths()
        if all(os.path.exists(p) for p in paths):
            return paths
    except ImportError:
        pass
    spec = os.environ.get("PETAKIT5D_TEST_FFTW", "")
    if spec:
        found = _pair(spec, _MATLAB_LIBS)
        if found is None:
            raise FileNotFoundError(f"PETAKIT5D_TEST_FFTW={spec}: no FFTW libraries there")
        return found
    return None


def _load_exact() -> bool:
    """Loads the exact backend once; whether it is there."""
    if _state["exact"] is None:
        paths = exact_libraries()
        if paths is None:
            _state["exact"] = False
        else:
            _ext().load_exact(*paths)
            _state["exact"] = True
    return bool(_state["exact"])


def available() -> tuple[str, ...]:
    """The backends this installation has: ``("permissive",)``, or both with the add-on."""
    return BACKENDS if _load_exact() else BACKENDS[:1]


def _initial() -> str:
    v = os.environ.get("PETAKIT5D_FFT", "").strip().lower() or "permissive"
    if v == "auto":
        return "exact" if _load_exact() else "permissive"
    if v not in BACKENDS:
        raise ValueError(f"PETAKIT5D_FFT={v!r}: use 'permissive', 'exact' or 'auto'")
    return v


def backend() -> str:
    """The backend of the calls that name none."""
    if _state["default"] is None:
        _state["default"] = _initial()
    return _state["default"]


def set_backend(name: str) -> None:
    """Makes ``name`` (``"permissive"`` or ``"exact"``) the default backend."""
    _backend_id(name)
    _state["default"] = name


@contextmanager
def using(name: str) -> Iterator[None]:
    """The default backend set to ``name`` inside the block."""
    old = backend()
    set_backend(name)
    try:
        yield
    finally:
        _state["default"] = old


def _backend_id(name: str | None) -> int:
    name = backend() if name is None else name
    if name == "permissive":
        return int(_ext().PERMISSIVE)
    if name == "exact":
        if not _load_exact():
            raise RuntimeError(
                "the exact FFT backend is not installed: install the petakit5d-fftw add-on "
                "(python/cpp/fft/exact, GPL) or set PETAKIT5D_FFTW to its libraries"
            )
        return int(_ext().EXACT)
    raise ValueError(f"unknown FFT backend {name!r}: use one of {BACKENDS}")


def info(name: str | None = None) -> dict[str, str]:
    """A backend's name, where it comes from and its version strings."""
    n, where, vf, vd = _ext().backend_info(_backend_id(name))
    return {"name": n, "where": where, "single": vf, "double": vd}


def _threads(threads: int | None) -> int:
    if threads:
        return int(threads)
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(None)


class HalfSpectrum:
    """Half-spectrum transforms of real (..., z, y, x) arrays of one shape, for FFT convolutions.

    ``forward(x)`` is the part of ``fftn(x)`` that the inverse of a conjugate-symmetric spectrum
    needs, ``multiply(s, h)`` the product with another spectrum of the same form (an OTF from
    ``otf(psf)``), ``inverse(s)`` the real result, ``real(ifftn(...))`` of the full spectrum.

    - exact backend: MATLAB's spectra, ``fftn(x)`` with index <= n / 2 along MATLAB's last
      transformed dimension (z of a volume, x of an image), and MATLAB's ifftn of the product: the
      results equal MATLAB's ``real(ifftn(H .* fftn(Y)))`` bit for bit, its stages, groups of lines
      and fill included;
    - permissive backend: plain real-to-complex transforms on the native layout, halved along x.

    The spectra's layout is the backend's (``spectrum_shape``): keep OTFs and spectra of one object
    together. ``matlab_parallel`` selects MATLAB's FFT mode (the exact backend only).
    """

    def __init__(
        self,
        shape: Sequence[int],
        dtype=np.float32,
        *,
        backend: str | None = None,
        threads: int | None = None,
        matlab_parallel: bool = True,
    ):
        self.shape = tuple(int(v) for v in shape)
        if not self.shape or any(v < 1 for v in self.shape):
            raise ValueError(f"a half spectrum needs lengths >= 1, got {self.shape}")
        self.dtype = np.dtype(dtype)
        if self.dtype not in (np.float32, np.float64):
            raise TypeError("the half spectrum takes float32 or float64 data")
        self.backend = globals()["backend"]() if backend is None else backend
        self._id = _backend_id(self.backend)
        self.threads = threads
        self._grain = int(_ext().GRAIN_PARALLEL) if matlab_parallel else 0
        if self.backend == "exact" and len(self.shape) > 1:
            self.spectrum_shape = tuple(_ext().half_shape(list(self.shape)))
        else:
            # permissive, and a 1-D array (a row: MATLAB's last dimension is its own)
            self.spectrum_shape = (*self.shape[:-1], self.shape[-1] // 2 + 1)

    def _real(self, x) -> np.ndarray:
        x = np.ascontiguousarray(x, dtype=self.dtype)
        if x.shape != self.shape:
            raise ValueError(f"expected an array of shape {self.shape}, got {x.shape}")
        return x

    def _spectrum(self, s) -> np.ndarray:
        cdt = np.complex64 if self.dtype == np.float32 else np.complex128
        if not (isinstance(s, np.ndarray) and s.dtype == cdt and s.flags.c_contiguous):
            raise TypeError(f"a spectrum of this object: C-contiguous {np.dtype(cdt)}")
        if s.shape != self.spectrum_shape:
            raise ValueError(f"expected a spectrum of shape {self.spectrum_shape}, got {s.shape}")
        return s

    def forward(self, x) -> np.ndarray:
        """The half spectrum of real data of this object's shape."""
        x = self._real(x)
        t = _threads(self.threads)
        if self.backend == "exact":
            return _ext().half_forward(x, t, self._grain, self._id).reshape(self.spectrum_shape)
        return _ext().native_forward(x, t)

    def inverse(self, s) -> np.ndarray:
        """The real result of the spectrum whose half is ``s`` (``s`` is overwritten)."""
        s = self._spectrum(s)
        t = _threads(self.threads)
        if self.backend == "exact":
            return (
                _ext()
                .half_inverse(s, list(self.shape), t, self._grain, self._id)
                .reshape(self.shape)
            )
        return _ext().native_inverse(s, list(self.shape), t)

    def multiply(self, s, h, conjugate: bool = False) -> np.ndarray:
        """``s .* h`` (``conjugate``: ``conj(h) .* s``) in place in ``s``: MATLAB's complex
        product."""
        s = self._spectrum(s)
        h = self._spectrum(h)
        _ext().half_multiply(s, h, bool(conjugate), _threads(self.threads))
        return s

    def otf(self, psf) -> np.ndarray:
        """The OTF of a PSF for this shape, as decon_psf2otf makes it: the PSF padded at the end to
        the shape, shifted so that its center (floor(size / 2) per axis) is at the origin, and
        transformed."""
        psf = np.asarray(psf, dtype=self.dtype)
        if psf.ndim != len(self.shape) or any(
            a > b for a, b in zip(psf.shape, self.shape, strict=True)
        ):
            raise ValueError(f"the PSF {psf.shape} must fit in {self.shape}")
        p = np.zeros(self.shape, dtype=self.dtype)
        p[tuple(slice(0, n) for n in psf.shape)] = psf
        p = np.roll(p, tuple(-(n // 2) for n in psf.shape), axis=tuple(range(p.ndim)))
        return self.forward(p)

    def convolve(self, x, h, conjugate: bool = False) -> np.ndarray:
        """``real(ifftn(h .* fftn(x)))`` (``conjugate``: ``conj(h)``) for a spectrum ``h`` of
        this object."""
        s = self.forward(x)
        self.multiply(s, h, conjugate)
        return self.inverse(s)
