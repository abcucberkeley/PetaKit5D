"""The optical transfer function of a PSF (``decon_psf2otf``, MATLAB's ``psf2otf``).

``psf2otf(psf, out_size)`` pads the PSF with zeros at the end to ``out_size``, shifts it
circularly so that its center ``floor(size / 2)`` (per axis, 0-based) is at the origin, and
transforms it with MATLAB's ``fftn`` (``matlabcompat.fft``): the exact FFT backend gives MATLAB's
bits, the permissive one (the default) the last bits of its library (decision record 0006). The
padding and the shift are one copy into the zero array; the transform is what MATLAB computes.

Arrays are the port's: ``(z, y, x)`` volumes and ``(y, x)`` images, C-contiguous; sizes are
Python shapes. As in MATLAB, sizes count three dimensions (a (y, x) image is one z plane), the
output size must not be smaller than the PSF along any of them, an all-zero PSF gives real
zeros of its class, and the result is complex (real where every imaginary part is zero, as
MATLAB returns it).
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import product

import numpy as np

__all__ = ["pad_shift_psf", "psf2otf"]


def _size3(shape: Sequence[int]) -> tuple[int, int, int]:
    """a Python shape as three (z, y, x) lengths (MATLAB's size padded with ones)"""
    s = tuple(int(v) for v in shape)
    if len(s) > 3:
        raise ValueError("psf2otf: sizes of at most three dimensions (z, y, x)")
    return (1,) * (3 - len(s)) + s


def _error(message: str):
    from petakit5d.decon.rl import DeconError

    return DeconError("", message)


def pad_shift_psf(psf: np.ndarray, out_shape: Sequence[int]) -> np.ndarray:
    """``circshift(padarray(psf, out - size, 'post'), -floor(size / 2))``: the PSF in a zero
    array of ``out_shape`` with element ``i`` at ``(i - size // 2) mod out`` per axis (MATLAB's
    padding with +0 and shift), in one copy. ``psf`` and ``out_shape`` as (z, y, x) or (y, x);
    the PSF must fit."""
    a = np.asarray(psf)
    ps, os_ = _size3(a.shape), _size3(out_shape)
    src = a.reshape(ps)
    out = np.zeros(os_, dtype=a.dtype)
    # per axis the two runs of the shift: source [h, p) to [0, p - h), [0, h) to [n - h, n)
    runs = []
    for p, n in zip(ps, os_, strict=True):
        h = p // 2
        r = [(slice(h, p), slice(0, p - h))]
        if h:
            r.append((slice(0, h), slice(n - h, n)))
        runs.append(r)
    for combo in product(*runs):
        out[tuple(d for _, d in combo)] = src[tuple(s for s, _ in combo)]
    return out.reshape(tuple(int(v) for v in out_shape))


def psf2otf(
    psf,
    out_size: Sequence[int] | None = None,
    *,
    backend: str | None = None,
    threads: int | None = None,
) -> np.ndarray:
    """``decon_psf2otf(psf, outSize)``: the OTF of ``psf`` for an array of ``out_size`` (a
    Python shape, (z, y, x) or (y, x); default the PSF's): the PSF padded at the end with
    zeros, shifted so that its center ``floor(size / 2)`` is at the origin, and transformed by
    MATLAB's ``fftn`` on ``backend`` (``petakit5d.fft``; None: the default backend). An all-zero
    (or empty) PSF gives zeros of its class; a size smaller than the PSF's along an axis raises
    MATLAB's error ``outSizeIsSmallerThanPsfSize`` (identifier '')."""
    from petakit5d.matlabcompat import fft as F

    a = np.asarray(psf)
    shape = a.shape if out_size is None else tuple(int(v) for v in np.asarray(out_size).ravel())
    if a.size and any(o < p for o, p in zip(_size3(shape), _size3(a.shape), strict=True)):
        raise _error("outSizeIsSmallerThanPsfSize")
    if not np.any(a != 0):
        return np.zeros(shape, dtype=a.dtype)
    return F.fftn(pad_shift_psf(a, shape), backend=backend, threads=threads)
