"""Synthetic lattice SIM data with a known pattern, for tests, examples and benchmarks.

Beads under a stationary lattice pattern along y, imaged with a 3D Gaussian PSF: the phase
image ``n`` of the bead at ``r_b`` is ``A_b g(r - r_b) I_n(y_b, z - z_b)``, with the
illumination ``I_n(y, dz) = 1 + sum_m 2 c_m(dz) cos(2 pi m (y / L - n s / L))`` of the lattice
period ``L`` and phase step ``s`` (um) and axial envelopes ``c_m`` (Gaussians, one per order m).
A bead image (``phase_psf``) under the same model gives the OTF that reconstructs it. Pages are
``(z planes * nphases, y, x)``, the phases of each plane in a row, uint16 after a background and
noise. Bit for bit the same on every platform (``testing.synth``'s exp, sine and cosine and
counter-based noise; sums in a fixed order).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = ["phase_psf", "phase_stack"]


def _pattern(y_um: float, dz_um: np.ndarray, n: int, *, lattice_period, phase_step, amps, sigmas):
    """``I_n(y, dz)`` for one bead's y and the axial offsets ``dz_um``"""
    from petakit5d.testing.synth import stable_exp, stable_sincos_2pi

    total = np.ones(dz_um.shape)
    for m, (a, s) in enumerate(zip(amps, sigmas, strict=True), start=1):
        env = a * stable_exp(-(dz_um * dz_um) / (2 * s * s))
        u = m * (y_um / lattice_period - n * phase_step / lattice_period)
        _, c = stable_sincos_2pi(np.array([u - np.floor(u)]))
        total = total + 2 * env * float(c[0])
    return total


def phase_stack(
    shape: Sequence[int],
    beads: Sequence[Sequence[float]],
    *,
    nphases: int = 5,
    pixel: Sequence[float] = (0.16074803849369898, 0.11, 0.11),
    sigma: Sequence[float] = (2.0, 1.3, 1.3),
    lattice_period: float = 1.2021,
    phase_step: float = 0.232,
    order_amps: Sequence[float] = (0.45, 0.2),
    order_sigmas: Sequence[float] = (0.9, 0.6),
    background: float = 100.0,
    noise: float = 3.0,
    seed: int = 0,
    dtype=np.uint16,
) -> np.ndarray:
    """Raw phases of beads ``(z, y, x, amplitude)`` (0-based voxel positions) in a volume of
    ``shape`` ``(nz, ny, nx)``: ``(nz * nphases, ny, nx)``; ``pixel`` and ``sigma`` are
    ``(z, y, x)`` (um and voxels)."""
    from petakit5d.core.dtypes import mround
    from petakit5d.testing.synth import CounterRNG, stable_exp

    nz, ny, nx = (int(v) for v in shape)
    sz, sy, sx = (float(v) for v in sigma)
    pz, py, _ = (float(v) for v in pixel)
    zz = np.arange(nz, dtype=np.float64)
    yy = np.arange(ny, dtype=np.float64)
    xx = np.arange(nx, dtype=np.float64)
    out = np.zeros((nz, nphases, ny, nx))
    for b in beads:
        bz, by, bx, amp = (float(v) for v in b)
        gz = stable_exp(-((zz - bz) ** 2) / (2 * sz * sz))
        gy = stable_exp(-((yy - by) ** 2) / (2 * sy * sy))
        gx = stable_exp(-((xx - bx) ** 2) / (2 * sx * sx))
        g2 = gy[:, None] * gx[None, :]
        for n in range(nphases):
            ill = _pattern(
                by * py, (zz - bz) * pz, n, lattice_period=lattice_period,
                phase_step=phase_step, amps=order_amps, sigmas=order_sigmas,
            )  # fmt: skip
            for z in range(nz):
                out[z, n] = out[z, n] + (amp * gz[z] * ill[z]) * g2
    rng = CounterRNG(seed, "sim/phase_stack")
    vals = out.reshape(nz * nphases, ny, nx) + background
    if noise:
        vals = vals + noise * rng.normal(vals.shape)
    vals = np.maximum(vals, 0)
    if np.dtype(dtype) == np.float32:
        return vals.astype(np.float32)
    return np.asarray(mround(vals)).astype(dtype)


def phase_psf(
    shape: Sequence[int] = (11, 32, 30),
    *,
    nphases: int = 5,
    pixel: Sequence[float] = (0.10716535899579933, 0.11, 0.11),
    offset: Sequence[float] = (0.3, -0.4, 0.25),
    amplitude: float = 20000.0,
    **kw,
) -> np.ndarray:
    """The phase-stepped image of one bead near the center (``offset`` voxels from
    ``(nz // 2, ny // 2, nx // 2)``), as a PSF for ``sim_psf_to_otf``."""
    nz, ny, nx = (int(v) for v in shape)
    bead = (nz // 2 + offset[0], ny // 2 + offset[1], nx // 2 + offset[2], amplitude)
    return phase_stack(shape, [bead], nphases=nphases, pixel=pixel, **kw)
