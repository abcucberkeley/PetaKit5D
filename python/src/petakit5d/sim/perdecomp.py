"""The periodic-plus-smooth decomposition of a volume (``perdecomp_3D`` in ``third_parties/sim3D``).

A new implementation from the published method (L. Moisan, "Periodic plus smooth image
decomposition", J. Math. Imaging Vis. 39(2), 2011): ``u = p + s`` with ``s`` smooth and ``p``
periodic. The boundary image ``v`` holds, at each face of the volume, the jump across the
periodic boundary (``v[0] = u[n-1] - u[0]``, ``v[n-1] = u[0] - u[n-1]`` along each axis, summed
over the axes); ``s`` solves the periodic discrete Poisson problem ``laplacian(s) = v``:

    s^(q, r, t) = v^(q, r, t) / (2 cos(2 pi q/M) + 2 cos(2 pi r/N) + 2 cos(2 pi t/P) - 6),

with ``s^(0) = 0``, and ``p = u - s``. ``perdecomp_3d`` returns ``p``, which has no jumps at
the borders, so its spectrum has no cross of edge artifacts.

The MATLAB function's code (another author's) is not read; its outputs were measured and agree
with this within the tolerance of ``docs/parity/sim.md``.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["perdecomp_3d"]


def _boundary(u, xp):
    """the boundary image: the jumps across the periodic boundary at each face, summed"""
    v = xp.zeros(u.shape, dtype=u.dtype)
    for axis in range(u.ndim):
        n = u.shape[axis]
        if n < 2:
            continue
        first = [slice(None)] * u.ndim
        last = [slice(None)] * u.ndim
        first[axis] = slice(0, 1)
        last[axis] = slice(n - 1, n)
        jump = u[tuple(last)] - u[tuple(first)]
        v[tuple(first)] += jump
        v[tuple(last)] -= jump
    return v


def _denominator(shape, dtype, xp):
    """``2 cos(2 pi q/M) + 2 cos(2 pi r/N) + 2 cos(2 pi t/P) - 2 ndim``, 1 at the origin"""
    from petakit5d.matlabcompat import elementary as E

    total = None
    for axis, n in enumerate(shape):
        c = 2.0 * np.asarray(E.cos(2.0 * math.pi * np.arange(n, dtype=np.float64) / n))
        sh = [1] * len(shape)
        sh[axis] = n
        term = (c - 2.0).reshape(sh)
        total = term if total is None else total + term
    total = np.broadcast_to(total, shape).copy()
    total.reshape(-1)[0] = 1.0
    return xp.asarray(total.astype(dtype))


def perdecomp_3d(u, *, use_gpu: bool = False, gpu_precision: str = "double", _ops=None):
    """The periodic component ``p`` of ``u`` (``perdecomp_3D(u, useGPU, gpuPrecision)``): an
    array of ``u``'s shape (``(z, y, x)`` or an image) of ``gpu_precision``, on the CPU as on
    the GPU."""
    from petakit5d.sim.recon import _div_parts, _Ops, _precision

    ops = _ops or _Ops(bool(use_gpu), gpu_precision, 1, None)
    xp = ops.xp
    # (computed in gpuPrecision on the CPU too: measured)
    dt = _precision(gpu_precision)
    a = xp.asarray(u).astype(dt)
    v = _boundary(a, xp)
    s_hat = _div_parts(ops.fftn(v), _denominator(a.shape, dt, xp))
    s_hat.reshape(-1)[0] = 0
    s = ops.ifftn(s_hat)
    s = s.real if s.dtype.kind == "c" else s
    return (a - s.astype(dt)).astype(dt)
