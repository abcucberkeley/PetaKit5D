"""The translation of a volume by a fractional shift (``imtranslate_function`` in
``third_parties/sim3D``, a modified copy of MathWorks' ``imtranslate``), as ``simRecon``'s Wiener
filter of ``wVer`` 0 uses it: linear interpolation, the output view of the input, 0 outside.

Written from the documented behaviour of ``imtranslate`` and measured MATLAB outputs (the
function's code is MathWorks' and not read): the output voxel ``(z, y, x)`` is the input at
``(z - tz, y - ty, x - tx)``, interpolated linearly along each axis, and 0 where that position
is outside the input's voxel grid along any axis (a sample between the last voxel and the
outside is 0, not blended with the fill value: measured). A translation is the same at every
voxel, so its weights are per axis.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

__all__ = ["imtranslate"]


def _shift_axis(a, axis: int, t: float, xp):
    """``a`` moved by ``t`` along ``axis``: ``out[i] = (1 - w) a[i - k] + w a[i - k - 1]``
    (``t = k + w``, ``0 <= w < 1``) where both samples are inside, else 0"""
    n = a.shape[axis]
    k = math.floor(t)
    w = t - k
    # output range whose samples are all inside
    lo = max(0, k + (1 if w > 0 else 0))
    hi = min(n, n + k)
    out = xp.zeros_like(a)
    if hi <= lo:
        return out

    def view(arr, start, stop):
        sl = [slice(None)] * arr.ndim
        sl[axis] = slice(start, stop)
        return arr[tuple(sl)]

    dst = view(out, lo, hi)
    dst += view(a, lo - k, hi - k) * (1.0 - w)
    if w > 0:
        dst += view(a, lo - k - 1, hi - k - 1) * w
    return out


def imtranslate(vol, translation: Sequence[float], *, _xp=np):
    """``imtranslate_function(vol, [tx ty tz])``: the ``(z, y, x)`` volume moved by ``tx``,
    ``ty``, ``tz`` voxels (MATLAB's [x y z] order of the translation), double."""
    xp = _xp
    tx, ty, tz = (float(v) for v in translation)
    a = xp.asarray(vol).astype(np.float64)
    a = _shift_axis(a, 2, tx, xp)
    a = _shift_axis(a, 1, ty, xp)
    return _shift_axis(a, 0, tz, xp)
