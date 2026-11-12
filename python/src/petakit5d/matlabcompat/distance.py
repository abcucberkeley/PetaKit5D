"""MATLAB's ``bwdist`` (Euclidean), exactly: the distance of every pixel to the nearest nonzero
pixel, as single. MATLAB computes the squared distance (exact integers) and returns
``sqrt(D)`` of it in single; the squared distance is unique, so any exact algorithm gives the
same bits, whatever the order of the axes. Images without nonzero pixels give Inf everywhere.

Images are NumPy's ``(y, x)`` or ``(z, y, x)``; the result has the image's shape."""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = ["bwdist"]


def bwdist(bw: Any, method: str = "euclidean", *, threads: int | None = None) -> np.ndarray:
    """MATLAB's ``D = bwdist(BW)`` of a 2D or 3D image (any class: nonzero is a feature):
    single, the Euclidean distance to the nearest nonzero pixel. Only the Euclidean method
    (and no second output) is supported."""
    from petakit5d import _ext
    from petakit5d.core.threads import resolve_threads

    m = method.lower()
    if not m or not "euclidean".startswith(m):
        raise NotImplementedError(f"bwdist: method {method!r} is not supported (euclidean only)")
    a = np.asarray(bw)
    if a.dtype.kind == "c":
        raise ValueError("bwdist: BW must be real")
    a = a if a.dtype == np.bool_ else a != 0
    if a.ndim > 3:
        raise NotImplementedError("bwdist: images of more than 3 dimensions are not supported")
    if a.ndim < 2:
        a = a.reshape((1, -1) if a.ndim == 1 else (1, 1))
    v = np.ascontiguousarray(a.reshape((1,) * (3 - a.ndim) + a.shape))
    out = np.empty(v.shape, dtype=np.float32)
    _ext.load("morph").bwdist_edt(v, out, resolve_threads(threads))
    return out.reshape(a.shape)
