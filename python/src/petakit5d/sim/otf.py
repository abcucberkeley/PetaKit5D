"""The SIM OTF from a phase-stepped bead image (``sim_PSFtoOTF_gen`` of ``third_parties/sim3D``)
and the bead's 3D Gaussian fit (``Lsq_GaussianFit_3D``).

``sim_psf_to_otf``: the background subtracted, the widefield PSF (the sum of the phase images
of each plane) fitted with a 3D Gaussian to find the bead's center, every phase image moved so
that the bead sits at the center that ``fftshift`` takes as the origin (no phase ramp in the
OTF), transformed, the orders separated with the pseudo-inverse of the separation matrix and
each orientation's zero order normalized to unit energy. Arrays are ``(z, y, x)``; the bead
image is ``(nz * nphases * norientations, y, x)`` with the phases interleaved per plane.

The fit is ``fminsearch`` (``sim.neldermead``) on the sum of squared residuals. Its objective
sums in MATLAB's order at one thread (MATLAB sums this many values on all its threads) with the
C library's exp, so the fit, and the OTF, are MATLAB's within the tolerances of
``docs/parity/sim.md``.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register
from petakit5d.sim.recon import ARRAY_OR_PATH

__all__ = ["gaussian_fit_3d", "sim_psf_to_otf"]


def gaussian_fit_3d(peak_img, display_fit: bool = False, *, num_threads: int | None = None):
    """``Lsq_GaussianFit_3D``: ``[amplitude, cy, cx, cz, sigma_y, sigma_x, sigma_z]`` of the 3D
    Gaussian that best fits ``peak_img`` ``(nz, ny, nx)`` in least squares (``fminsearch`` from
    the maximum voxel and sigmas of 1.5, ``TolX`` and ``TolFun`` 1e-6, at most 10,000 iterations
    and evaluations). The center is 1-based, in MATLAB's subscripts (y, x, z), as MATLAB's
    function returns it; ``display_fit`` (MATLAB's figures) is not supported."""
    from petakit5d.core.threads import resolve_threads
    from petakit5d.kernels.ops import sum_double
    from petakit5d.matlabcompat import reduce
    from petakit5d.sim.neldermead import fminsearch

    if display_fit:
        raise NotImplementedError("display_fit (MATLAB's figures of the fit) is not ported")
    img = np.ascontiguousarray(np.asarray(peak_img, dtype=np.float64))
    if img.ndim == 2:
        img = img[None]
    nz, ny, nx = img.shape
    nt = resolve_threads(num_threads)
    # [maxval, ind] = max(PeakImg(:)): the first maximum in MATLAB's order
    m, idx = reduce.max(img)
    iz, iy, ix = (int(v) for v in idx)
    guess = [float(m), iy + 1.0, ix + 1.0, iz + 1.0, 1.5, 1.5, 1.5]
    from petakit5d.sim.recon import _par

    yy = np.arange(1, ny + 1, dtype=np.float64)
    xx = np.arange(1, nx + 1, dtype=np.float64)
    zz = np.arange(1, nz + 1, dtype=np.float64)
    r2 = np.empty(img.shape)

    def objective(p: np.ndarray) -> float:
        # (data - p1 exp(-(dy^2/(2 sy^2) + dx^2/(2 sx^2) + dz^2/(2 sz^2)))).^2 per voxel, the
        # planes on the threads, then the sum in MATLAB's order
        ty = (yy - p[1]) ** 2 / (2 * p[4] ** 2)
        tx = (xx - p[2]) ** 2 / (2 * p[5] ** 2)
        tz = (zz - p[3]) ** 2 / (2 * p[6] ** 2)
        tyx = ty[:, None] + tx[None, :]

        def run(lo, hi):
            with np.errstate(over="ignore", invalid="ignore"):
                r = img[lo:hi] - p[0] * np.exp(-(tyx[None] + tz[lo:hi, None, None]))
            r2[lo:hi] = r * r

        _par(nz, run, nt)
        return sum_double(r2, num_threads=nt)

    res = fminsearch(
        objective, guess, tol_x=1e-6, tol_fun=1e-6, max_iter=10000, max_fun_evals=10000
    )
    return res.x


def _read_psf(psf: Any) -> np.ndarray:
    if isinstance(psf, (str, os.PathLike)):
        from petakit5d.io import read_tiff

        return np.asarray(read_tiff(os.fspath(psf)))
    a = np.asarray(psf)
    return a if a.ndim == 3 else a[None]


def _otf_to_matlab(o):
    from petakit5d.core.arrays import to_matlab

    return to_matlab(np.asarray(o))


@register(
    "sim_PSFtoOTF_gen",
    Param(
        "PSF",
        None,
        "Phase-stepped bead image (z * phases * orientations, y, x), or its tiff file.",
        convert=ARRAY_OR_PATH,
    ),
    Param("nphases", "@isnumeric", "Phases."),
    Param("norders", "@isnumeric", "Information orders."),
    Param("norientations", "@isnumeric", "Pattern orientations."),
    Param("lattice_period", "@isnumeric", "Lattice period in um (the coarsest)."),
    Param("phase_step", "@isnumeric", "Phase step in um."),
    Param("pxl_dim_PSF", "@isnumeric", "Voxel size (z, y, x) in um.", convert="size_yxz"),
    Param("Background", "@isnumeric", "Background subtracted.", convert="same"),
    Param("displayFit", "@islogical", "Show the Gaussian fit (MATLAB's figures; not ported)."),
    Param("Save16bit", "@islogical", "Not used."),
    Param("gpuPrecision", "@ischar", "Precision of the phase ramps and GPU arrays."),
    Param("useGPU", "@islogical", "Run on a GPU if one is usable."),
    Param("saveOTF", "@islogical", "Save the OTF (fails, as MATLAB's: see the docstring)."),
    Param("Overwrite", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    Param("fft_backend", None, "CPU FFT backend (None: the default).", python_only=True),
    task="C06",
    matlab_out=_otf_to_matlab,
)
def sim_psf_to_otf(
    psf,
    *,
    nphases: float = 5,
    norders: float = 5,
    norientations: float = 1,
    lattice_period: float = 1.2021,
    phase_step: float = 0.232,
    pxl_dim_psf: Sequence[float] = (0.10716535899579933, 0.11, 0.11),
    background: float = 105,
    display_fit: bool = False,
    save_16bit: bool = False,
    gpu_precision: str = "double",
    use_gpu: bool = True,
    save_otf: bool = False,
    overwrite: bool = False,
    num_threads: int | None = None,
    fft_backend: str | None = None,
):
    """The OTF orders ``(norders, z, y, x)`` (``(norientations, norders, z, y, x)`` for several
    orientations), complex double, of a phase-stepped bead image or tiff file
    (``sim_PSFtoOTF_gen``). ``save_otf`` fails as MATLAB's does (its file name is taken from the
    image, not the file name)."""
    from petakit5d.core.threads import resolve_threads
    from petakit5d.sim.recon import (
        SimError,
        _div_parts,
        _Ops,
        _subtract_background,
        _sum_vector,
        _yxz_of,
        fourier_shift,
        inverse_separation_matrix,
        separation_matrix,
    )

    del save_16bit, overwrite
    gpu = False
    if use_gpu:
        from petakit5d.gpu.device import usable_gpu_count

        gpu = usable_gpu_count()[0] > 0
    nt = resolve_threads(num_threads)
    ops = _Ops(gpu, gpu_precision, nt, fft_backend)
    xp, p = ops.xp, ops.p
    a = _subtract_background(_read_psf(psf), float(background))
    nph, nor, nori = int(nphases), int(norders), int(norientations)
    nimgs, ny, nx = a.shape
    nz_f = nimgs / (float(nphases) * float(norientations))
    pxp = _yxz_of(pxl_dim_psf)
    dk_psf = 1.0 / (np.array([ny, nx, nz_f], dtype=np.float64) * pxp)
    if gpu:
        dk_psf = dk_psf.astype(p)
    if nz_f != int(nz_f):
        raise SimError("Size inputs must be integers.", "MATLAB:NonIntegerInput")
    nz = int(nz_f)
    step = nph * nori

    # the widefield PSF: the sum of each plane's phase images (in index order, double)
    wide = np.float64 if not gpu or p == np.float64 else np.float32
    wf = np.empty((nz, ny, nx), dtype=wide)
    for ii in range(nz):
        s = a[ii * step].astype(wide)
        for q in range(1, step):
            s = s + a[ii * step + q].astype(wide)
        wf[ii] = s
    fit = gaussian_fit_3d(wf.astype(np.float64), bool(display_fit), num_threads=nt)

    inv_sep = inverse_separation_matrix(
        separation_matrix(nph, nor, float(lattice_period), float(phase_step)), gpu_precision
    )
    store = ops.pc if gpu else np.complex128
    o = xp.zeros((nori, nor, nz, ny, nx), dtype=store)
    shift = [
        np.ceil((ny + 1) / 2) - fit[1],
        np.ceil((nx + 1) / 2) - fit[2],
        np.ceil((nz + 1) / 2) - fit[3],
    ]
    m = xp.asarray(inv_sep if gpu else inv_sep.astype(np.complex128))
    pdk = np.prod(dk_psf)
    for jj in range(nori):
        dks = []
        for ii in range(nph):
            dr = xp.asarray(a[ii + jj * nph :: step].astype(p if gpu else np.float64))
            sh = fourier_shift(dr, shift, gpu_precision, _ops=ops)
            dk = ops.fftshift(ops.ifftn(ops.ifftshift(sh.astype(store))))
            dks.append(_div_parts(dk.astype(store), pdk))
            del dr, sh, dk
        for j in range(nor):
            acc = dks[0] * m[j, 0]
            for kk in range(1, nph):
                acc = acc + dks[kk] * m[j, kk]
            o[jj, j] = acc if gpu or p == np.float64 else acc.astype(np.complex64)
            del acc
        del dks
    zero = int(np.ceil(nor / 2)) - 1
    for jj in range(nori):
        o0 = o[jj, zero]
        energy = _sum_vector(((o0.real * o0.real + o0.imag * o0.imag) * pdk).ravel(), gpu)
        o[jj] = _div_parts(o[jj], float(np.sqrt(energy)))
    o = ops.host(o)
    if save_otf:
        # MATLAB's fileparts(PSF) of the image array
        raise SimError(
            "Input must be a row vector of characters, or a string scalar.",
            "MATLAB:fileparts:MustBeChar",
        )
    return o[0] if nori == 1 else o
