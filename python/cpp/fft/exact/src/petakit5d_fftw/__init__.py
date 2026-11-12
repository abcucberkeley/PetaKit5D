"""FFTW 3.3.8 with MATLAB R2024b's codelet sets (GPL-2.0-or-later): the exact CPU FFT backend of
petakit5d. Installing this package lets ``petakit5d.fft.set_backend("exact")`` give MATLAB's bits;
petakit5d itself contains no GPL code and runs without it (pocketfft, the permissive backend).
"""

from __future__ import annotations

import os
import sys

__all__ = ["library_paths", "VERSION"]

VERSION = "3.3.8"


def library_paths() -> tuple[str, str]:
    """The single and double precision libraries of this package."""
    here = os.path.dirname(os.path.abspath(__file__))
    if sys.platform == "win32":
        names = ("pk5d_fftw3f.dll", "pk5d_fftw3.dll")
    elif sys.platform == "darwin":
        names = ("libpk5d_fftw3f.dylib", "libpk5d_fftw3.dylib")
    else:
        names = ("libpk5d_fftw3f.so", "libpk5d_fftw3.so")
    return tuple(os.path.join(here, n) for n in names)
