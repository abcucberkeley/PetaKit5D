"""The PSF sigma estimate of 3-D point detection: MATLAB's ``GU_estimateSigma3D(file, [])``.

Three rounds of the 'speedup' detector on the bead image (``detect.point_source_detection_speedup``
in mode ``'xyzAsrc'``, the sigmas free, with the fitGaussian3D mex's fits), starting from
``[1.2 2.5]`` and then from the previous round's result; each round takes the brightest
detection's sigmas and position. The demos then divide ``sigma_z`` by the condition's
``zAniso`` (``z_aniso_sigma``).

As in MATLAB, the detector's t threshold is computed once per session and alpha and shared with
every later 'speedup' detection (``detect.clear_speedup_cache`` starts a new session).

``estimate_sigma3d`` is registered with the lab's interface (``rt`` and ``FileName`` by position,
no validators); the MATLAB-name layer gives MATLAB's outputs (1-based positions, the mask in
MATLAB's layout).
"""

from __future__ import annotations

import os

import numpy as np

from petakit5d.params import Param, register

__all__ = ["estimate_sigma3d", "psf_volume", "z_aniso_sigma"]


def psf_volume(path_or_volume) -> np.ndarray:
    """the bead image as the estimate reads it: double, zeros as NaN (MATLAB's), (z, y, x)"""
    from petakit5d.core.dtypes import NAN

    if isinstance(path_or_volume, (str, os.PathLike)):
        from petakit5d.io import read_tiff

        v = read_tiff(str(path_or_volume)).astype(np.float64)
    else:
        v = np.array(path_or_volume, dtype=np.float64)
    v[v == 0] = NAN
    return v


def _sigma_to_matlab(result: tuple) -> tuple:
    """MATLAB's outputs: XYZ and the detections with 1-based positions, the mask (y, x, z)"""
    from petakit5d.core.arrays import to_matlab
    from petakit5d.detection.run import to_matlab_fields

    sxy, sz, xyz, ps, mask = result
    return sxy, sz, {k: v + 1 for k, v in xyz.items()}, to_matlab_fields(ps), to_matlab(mask)


@register(
    "GU_estimateSigma3D",
    Param("rt", None, "The bead image's file (or, in the port, its volume).", py="path_or_volume"),
    Param("FileName", None, "A file name appended to rt ([]: rt is the file)."),
    Param("rounds", None, "Detection rounds (3, as MATLAB).", python_only=True),
    Param(
        "start", None, "The first round's sigmas [xy z] ([1.2 2.5], as MATLAB).", python_only=True
    ),
    Param("numThreads", None, "Threads (None: the CPUs available).", python_only=True),
    task="C02",
    matlab_out=_sigma_to_matlab,
)
def estimate_sigma3d(
    path_or_volume,
    file_name=None,
    *,
    rounds: int = 3,
    start=(1.2, 2.5),
    num_threads: int | None = None,
) -> tuple[float, float, dict, dict, np.ndarray]:
    """``[sigmaXY, sigmaZ, XYZ, pstruct, mask] = GU_estimateSigma3D(rt, FileName)``.

    Returns the brightest detection's ``sigma_xy`` and ``sigma_z`` (voxels), its position
    ``{'x', 'y', 'z'}`` (0-based, x along the last axis), and the last round's detections and
    mask. A tie of brightest detections, which MATLAB would carry on as vectors, raises. A
    ``file_name`` is appended to the path as MATLAB's ``[rt FileName]``."""
    from petakit5d.detection import detect

    if file_name is not None and np.size(file_name) != 0:
        path_or_volume = str(path_or_volume) + str(file_name)
    vol = psf_volume(path_or_volume)
    sigma = (float(start[0]), float(start[1]))
    ps: dict = {}
    mask = np.zeros(vol.shape, dtype=bool)
    xyz: dict = {}
    for _ in range(int(rounds)):
        ps, mask = detect.point_source_detection_speedup(
            vol,
            sigma,
            mode="xyzAsrc",
            alpha=0.05,
            remove_redundant_fits=True,
            window=None,
            fit_method="original",
            num_threads=num_threads,
        )
        A = np.asarray(ps["A"], dtype=np.float64)
        if A.size == 0:
            raise ValueError("GU_estimateSigma3D: no detection in the PSF image")
        idx = np.nonzero(A == np.max(A))[0]
        if idx.size != 1:
            raise ValueError("GU_estimateSigma3D: several detections share the largest amplitude")
        i = int(idx[0])
        s = np.asarray(ps["s"], dtype=np.float64)
        sigma = (float(s[0, i]), float(s[1, i]))
        xyz = {k: float(np.asarray(ps[k]).ravel()[i]) for k in ("x", "y", "z")}
    return sigma[0], sigma[1], xyz, ps, mask


def z_aniso_sigma(sigma_xy: float, sigma_z: float, z_aniso: float) -> tuple[float, float]:
    """the demos' sigma for the deskewed frames: ``[sigmaXY, sigmaZ ./ zAniso]`` (the bead image
    is about isotropic, the frames are not)"""
    return float(sigma_xy), float(sigma_z) / float(z_aniso)
