"""Tests of petakit5d.fft: the two CPU FFT backends, the license boundary of the default build, and
the half-spectrum transforms FFT-based kernels convolve with."""

import os
import re
from pathlib import Path

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.matlabcompat import fft as F

BACKENDS = P.available()
EXACT = "exact" in BACKENDS
needs_exact = pytest.mark.skipif(
    not EXACT, reason="the exact FFT backend (the FFTW add-on) is not installed"
)


def module_file():
    return P._ext().__file__


def test_the_default_build_contains_no_fftw():
    # FFTW (GPL) is reached only through the backend table: its code, codelets and version strings
    # are not in the module, and the module links no FFTW library
    data = open(module_file(), "rb").read()
    for marker in (b"fftw-3.3", b"n1fv_", b"t3bv_", b"r2cf_64", b"Matteo Frigo"):
        assert marker not in data, marker
    assert b"pocketfft-cpp-" in data
    if os.name == "posix" and os.path.exists("/proc/self/maps"):
        # Linux: FFTW libraries are mapped only when the exact backend was loaded
        libs = {ln.split()[-1] for ln in open("/proc/self/maps") if ln.rstrip().endswith(".so")}
        fftw = [p for p in libs if re.search(r"fftw", os.path.basename(p))]
        if not EXACT:
            assert not fftw, fftw


def test_only_the_exact_component_includes_fftw3_h():
    # the exact backend is reached only through the backend table (cpp/fft/backend.h)
    cpp = Path(__file__).resolve().parents[2] / "cpp"
    if not cpp.is_dir():
        pytest.skip("the C++ sources are not here")
    include = re.compile(r"^\s*#\s*include\s*[<\"][^>\"]*fftw3[^>\"]*\.h[>\"]")
    found = []
    for p in cpp.rglob("*"):
        if p.suffix not in (".c", ".cc", ".cpp", ".h", ".hpp", ".cu", ".cuh", ".inl"):
            continue
        if p.relative_to(cpp).parts[:2] == ("fft", "exact"):
            continue
        for k, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if include.search(line):
                found.append(f"{p.relative_to(cpp)}:{k}")
    assert not found, found


def test_the_default_backend_is_permissive():
    env = os.environ.get("PETAKIT5D_FFT", "").strip().lower()
    if env and env != "permissive":
        pytest.skip(f"PETAKIT5D_FFT={env} sets another default")
    assert P.backend() == "permissive"
    info = P.info()
    assert info["name"] == "permissive" and info["single"].startswith("pocketfft-cpp-")


def test_unknown_or_missing_backends_raise():
    with pytest.raises(ValueError):
        P.set_backend("fftw")
    if not EXACT:
        with pytest.raises(RuntimeError):
            P.set_backend("exact")
        assert P.available() == ("permissive",)


def test_using_restores_the_default():
    before = P.backend()
    with P.using("permissive"):
        assert P.backend() == "permissive"
    assert P.backend() == before


@needs_exact
def test_both_backends_run_in_one_process_with_their_own_results():
    rng = np.random.default_rng(3)
    x = rng.random((40, 64, 37)).astype(np.float32)
    a = F.fftn(x, backend="exact")
    b = F.fftn(x, backend="permissive")
    assert a.tobytes() != b.tobytes()
    assert np.max(np.abs(a - b)) <= 1e-5 * np.max(np.abs(a))
    # each again after the other: the same bits
    assert F.fftn(x, backend="exact").tobytes() == a.tobytes()
    assert F.fftn(x, backend="permissive").tobytes() == b.tobytes()
    with P.using("exact"):
        assert F.fftn(x).tobytes() == a.tobytes()
        assert P.info()["single"].startswith("fftw-3.3")


# --- half spectra ---------------------------------------------------------------------------------


def matlab_product(s, h, conjugate=False):
    """MATLAB's complex product (no fused multiply-add), h .* s or conj(h) .* s"""
    sr, si, hr, hi = s.real, s.imag, h.real, h.imag
    out = np.empty(s.shape, dtype=s.dtype)
    if conjugate:
        out.real = sr * hr + si * hi
        out.imag = si * hr - sr * hi
    else:
        out.real = sr * hr - si * hi
        out.imag = sr * hi + si * hr
    return out


@pytest.mark.parametrize("be", BACKENDS)
@pytest.mark.parametrize("shape", [(6, 10, 12), (5, 9, 7), (12, 10), (1, 16, 20), (7, 1, 9), (64,)])
def test_half_spectrum_shapes_and_round_trip(be, shape):
    hs = P.HalfSpectrum(shape, np.float32, backend=be)
    rng = np.random.default_rng(len(shape))
    x = rng.random(shape).astype(np.float32)
    s = hs.forward(x)
    assert s.shape == hs.spectrum_shape and s.dtype == np.complex64
    y = hs.inverse(s.copy())
    assert y.shape == shape and y.dtype == np.float32
    assert np.max(np.abs(y - x)) <= 1e-5 * np.max(np.abs(x))
    if be == "permissive":
        assert hs.spectrum_shape == (*shape[:-1], shape[-1] // 2 + 1)
        ref = np.fft.rfftn(x.astype(np.float64))
        assert np.max(np.abs(s - ref)) <= 1e-5 * np.max(np.abs(ref))


@needs_exact
@pytest.mark.parametrize(
    "shape", [(6, 10, 12), (5, 9, 7), (12, 10), (1, 16, 20), (64, 50, 40), (2, 64, 33)]
)
def test_exact_half_spectrum_is_matlabs_fftn_and_ifftn(shape):
    # forward: fftn's elements with index <= n / 2 along MATLAB's last transformed dimension;
    # inverse: MATLAB's ifftn of the conjugate-symmetric product, bit for bit
    rng = np.random.default_rng(7)
    hs = P.HalfSpectrum(shape, np.float32, backend="exact")
    x = rng.random(shape).astype(np.float32)
    psf = rng.random(tuple(max(1, n // 3) for n in shape)).astype(np.float32)
    full = F.fftn(x, backend="exact")
    s = hs.forward(x)
    last = _matlab_last_axis(shape)
    half = [slice(None)] * len(shape)
    half[last] = slice(0, shape[last] // 2 + 1)
    assert s.tobytes() == np.ascontiguousarray(full[tuple(half)]).tobytes()
    # the OTF as decon_psf2otf makes it, and real(ifftn(H .* fftn(x)))
    p = np.zeros(shape, np.float32)
    p[tuple(slice(0, n) for n in psf.shape)] = psf
    p = np.roll(p, tuple(-(n // 2) for n in psf.shape), axis=tuple(range(len(shape))))
    H = F.fftn(p, backend="exact")
    assert hs.otf(psf).tobytes() == np.ascontiguousarray(H[tuple(half)]).tobytes()
    ref = F.ifftn(matlab_product(full, H), backend="exact")
    assert ref.dtype == np.float32
    assert hs.convolve(x, hs.otf(psf)).tobytes() == ref.tobytes()
    ref = F.ifftn(matlab_product(full, H, True), backend="exact")
    assert hs.convolve(x, hs.otf(psf), conjugate=True).tobytes() == ref.tobytes()


def _matlab_last_axis(shape):
    """The NumPy axis of MATLAB's last transformed dimension of a (z, y, x) shape (x at least)."""
    msize = (shape[-2], shape[-1], *shape[-3::-1]) if len(shape) >= 2 else (1, shape[0])
    last = 2
    for k, n in enumerate(msize):
        if n > 1:
            last = max(last, k + 1)
    if len(shape) == 1:
        return 0
    return {1: -2, 2: -1}.get(last, -last) % len(shape)


@pytest.mark.parametrize("be", BACKENDS)
def test_half_spectrum_convolution_is_close_to_double_precision(be):
    rng = np.random.default_rng(11)
    shape = (16, 40, 33)
    x = rng.random(shape).astype(np.float32) * 1000
    psf = rng.random((7, 9, 9)).astype(np.float32)
    psf /= psf.sum()
    hs = P.HalfSpectrum(shape, np.float32, backend=be)
    y = hs.convolve(x, hs.otf(psf))
    p = np.zeros(shape)
    p[:7, :9, :9] = psf
    p = np.roll(p, (-3, -4, -4), axis=(0, 1, 2))
    ref = np.fft.irfftn(
        np.fft.rfftn(x.astype(np.float64)) * np.fft.rfftn(p), s=shape, axes=(0, 1, 2)
    )
    d = np.abs(y - ref)
    assert d.max() <= 1e-5 * np.abs(ref).max()
    assert np.linalg.norm(d) <= 1e-6 * np.linalg.norm(ref)


@pytest.mark.parametrize("be", BACKENDS)
def test_half_spectra_do_not_depend_on_the_thread_count(be):
    rng = np.random.default_rng(5)
    x = rng.random((20, 64, 37)).astype(np.float32)
    ref = P.HalfSpectrum(x.shape, backend=be, threads=1)
    s1 = ref.forward(x)
    y1 = ref.inverse(s1.copy())
    for t in (2, 3, 8):
        hs = P.HalfSpectrum(x.shape, backend=be, threads=t)
        assert hs.forward(x).tobytes() == s1.tobytes()
        assert hs.inverse(s1.copy()).tobytes() == y1.tobytes()


def test_half_spectrum_checks_its_arrays():
    hs = P.HalfSpectrum((8, 10), np.float32)
    with pytest.raises(ValueError):
        hs.forward(np.zeros((8, 11), np.float32))
    with pytest.raises(TypeError):
        hs.inverse(np.zeros(hs.spectrum_shape, np.complex128))
    with pytest.raises(ValueError):
        hs.otf(np.ones((9, 3), np.float32))
    with pytest.raises(TypeError):
        P.HalfSpectrum((8, 10), np.int16)
