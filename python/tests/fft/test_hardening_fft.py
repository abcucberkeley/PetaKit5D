"""Tests the mutation campaign asked for (E06, docs/hardening.md): the FFT interface."""

import numpy as np
import pytest

from petakit5d import fft


@pytest.mark.parametrize("shape", [(), (0, 4), (3, 0, 2)])
def test_half_spectrum_needs_lengths_of_at_least_one(shape):
    with pytest.raises(ValueError, match="lengths >= 1"):
        fft.HalfSpectrum(shape, np.float32)


@pytest.mark.parametrize("exact", [False, True])
def test_auto_takes_the_exact_backend_when_it_is_there(monkeypatch, exact):
    monkeypatch.setenv("PETAKIT5D_FFT", "auto")
    monkeypatch.setattr(fft, "_load_exact", lambda: exact)
    assert fft._initial() == ("exact" if exact else "permissive")
    monkeypatch.setenv("PETAKIT5D_FFT", "fftw")
    with pytest.raises(ValueError, match="PETAKIT5D_FFT"):
        fft._initial()
