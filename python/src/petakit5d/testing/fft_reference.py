"""Whether results recorded from MATLAB's FFTs can be reproduced bit for bit here.

The exact FFT backend (``petakit5d.fft``: FFTW 3.3.8 with MATLAB R2024b's codelets) gives MATLAB's
bits where FFTW makes MATLAB's plans and MATLAB's twiddle factors:

- the plans: FFTW picks its codelets from the CPU, so they are MATLAB's on x86-64 with AVX;
- the twiddle factors come from the C library (FFTW's kernel/trig.c, in double for both
  precisions), so their last bits, and with them the FFTs', follow its version and the code path
  it takes on the CPU. GCC merges real_cexp's cos and sin calls into one call of ``sincos``, in
  MATLAB's ``libmwfftw3f`` / ``libmwfftw3`` as in our build, and glibc's ``sincos`` is not always
  its ``sin`` and ``cos``: on glibc 2.28 and 2.35 they differ in the last bit at 944 of FFTW's
  arguments for lengths up to 2048, while glibc 2.39's ``sincos`` gives ``sin``'s and ``cos``'s
  values. MATLAB's references (digests of its outputs, oracle cases with FFTs) were recorded on
  glibc 2.28 (Intel x86-64); where the C library gives FFTW other values, MATLAB's own FFTW
  would not reproduce them either.

``digests_hold_here()`` gives ``""`` where MATLAB's recorded FFT bits apply, else the reason (for a
test's skip): no exact backend, other plans, or other twiddle factors (``twiddle_fingerprint()``
against ``REFERENCE_TWIDDLE_FINGERPRINT``).
"""

from __future__ import annotations

import ctypes
import ctypes.util
import functools
import hashlib
import math
import platform
import sys

import numpy as np

__all__ = [
    "K2PI",
    "REFERENCE_TWIDDLE_FINGERPRINT",
    "c_sincos",
    "digests_hold_here",
    "exact_plan",
    "fftw_sincos",
    "imports_sincos",
    "matlab_platform",
    "separate_sincos",
    "twiddle_arguments",
    "twiddle_fingerprint",
]

# twiddle_fingerprint() where MATLAB's references were recorded (glibc 2.28's sincos, Intel x86-64)
REFERENCE_TWIDDLE_FINGERPRINT = "9c27d5b80b58341c8c9f1ee72499b969dc7673edb6a30793fd23fd87761327a6"
K2PI = 6.2831853071795864769252867665590057683943388  # FFTW's K2PI, rounded to double as C does


def twiddle_arguments(n: int) -> np.ndarray:
    """The arguments FFTW's real_cexp takes the cosine and sine of for 2 pi m / n, m = 0 .. n - 1:
    4m and 4n reduced to the first octant, then (K2PI * m) / n in double."""
    m = 4 * np.arange(n, dtype=np.int64)
    n4 = 4 * n
    m = np.where(m > n4 - m, n4 - m, m)
    m = np.where(m - n > 0, m - n, m)
    m = np.where(m > n - m, n - m, m)
    return (K2PI * np.unique(m).astype(np.float64)) / float(n4)


def separate_sincos(v: float) -> tuple[float, float]:
    """(sin(v), cos(v)) from the C library's sin and cos."""
    return math.sin(v), math.cos(v)


def c_sincos():
    """The C library's sincos(x, &s, &c) as a function of x giving (s, c), or None."""
    if sys.platform == "win32":
        return None
    try:
        f = ctypes.CDLL(ctypes.util.find_library("m") or "libm.so.6").sincos
    except (OSError, AttributeError):
        return None
    f.argtypes = [ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
    f.restype = None
    s, c = ctypes.c_double(), ctypes.c_double()

    def sincos(v):
        f(v, ctypes.byref(s), ctypes.byref(c))
        return s.value, c.value

    return sincos


def imports_sincos(path: str) -> bool:
    """Whether a library names sincos among its symbols (its dynamic string table)."""
    with open(path, "rb") as fh:
        return b"\0sincos\0" in fh.read()


def _exact() -> bool:
    from petakit5d import fft

    return "exact" in fft.available()


def fftw_sincos():
    """The function FFTW's twiddle factors come from: the C library's sincos where the exact
    backend's libraries import it (built with GCC), else its sin and cos."""
    from petakit5d import fft

    paths = fft.exact_libraries() if _exact() else None
    f = c_sincos()
    if paths and f is not None and all(imports_sincos(q) for q in paths):
        return f
    return separate_sincos


def twiddle_fingerprint(nmax: int = 2048, sincos=None) -> str:
    """sha256 of the cosines and sines of every argument of twiddle_arguments(n), n = 1 .. nmax
    (the lengths of the oracle packs and of FFTW's sub-plans for them), from sincos (default: the
    function FFTW calls here)."""
    sincos = fftw_sincos() if sincos is None else sincos
    h = hashlib.sha256()
    for n in range(1, nmax + 1):
        sc = [sincos(v) for v in twiddle_arguments(n).tolist()]
        h.update(np.array([c for _, c in sc] + [s for s, _ in sc], dtype=np.float64).tobytes())
    return h.hexdigest()


def exact_plan(kind: str, n: int, howmany: int, double: bool = False) -> str:
    """The exact backend's plan for kind r2c, c2c-1, c2c+1 or c2r, whitespace collapsed."""
    from petakit5d import fft

    ext = fft._ext()
    return " ".join(ext.plan_string(fft._backend_id("exact"), double, kind, n, howmany).split())


def matlab_platform() -> bool:
    """Whether the exact backend's plans can be MATLAB's here: x86-64 with AVX (FFTW picks its
    codelets from the CPU)."""
    if not _exact() or platform.machine().lower() not in ("x86_64", "amd64"):
        return False
    return exact_plan(
        "r2c", 64, 4
    ) == '(rdft2-r2hc-direct-64-x4 "r2cf_64")' and "avx" in exact_plan("c2c-1", 64, 4)


@functools.cache
def digests_hold_here() -> str:
    """``""`` where MATLAB's recorded FFT bits apply, else why not (no exact backend, or FFTW's
    plans or twiddle factors differ from MATLAB's)."""
    if not _exact():
        return "the exact FFT backend (the FFTW add-on) is not installed"
    if not matlab_platform():
        return "FFTW's plans differ from MATLAB's here (not x86-64 with AVX)"
    fp = twiddle_fingerprint()
    if fp != REFERENCE_TWIDDLE_FINGERPRINT:
        return (
            "the C library gives FFTW other twiddle factors than where MATLAB's references were "
            f"recorded (fingerprint {fp[:16]}, not {REFERENCE_TWIDDLE_FINGERPRINT[:16]}): MATLAB's "
            "own FFTW would differ from them here too"
        )
    return ""
