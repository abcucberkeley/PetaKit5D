"""Deskew and rotation of a volume in memory (``deskewRotateFrame3D``, ``deskewFrame3D``,
``rotateFrame3D``) and of tiff files (``XR_RotateFrame3D``).

``deskew_rotate_volume`` is the lab's ``deskewRotateFrame3D``: the geometry of
``deskew_rotate.geometry``, then

- the deskew kernels (``kernels.deskew``, the mex functions) for a sample scan with linear
  interpolation and no resampling on the CPU: the warp alone, or for x steps above
  ``xStepThresh`` the skewed space interpolation and warp in one pass (the fused kernel);
- ``imwarp`` (``matlabcompat.warp``, MATLAB's own) otherwise (objective scan, resampling, cubic
  interpolation, ``gpuProcess``), after the skewed space interpolation (its kernel) where the x
  step is above the threshold;
- MATLAB's code paths where a mex rejects its input: classes other than uint16 and single, and
  inputs narrower than the x step, go to the interpolation in MATLAB code
  (``skewed_space_interp_matlab``, the lab's ``skewed_space_interp_defined_stepsize.m``) and
  ``imwarp`` as in MATLAB (without MATLAB's ``PetaKit5D:mexFallback`` warning: the port has no
  mex to miss).

Results are bit for bit MATLAB's on every path. ``gpu_process`` runs MATLAB's GPU code path
(``imwarp`` instead of the kernels) on the CPU, bit for bit MATLAB's ``imwarp`` on the CPU;
MATLAB's own GPU result differs from that by 1 count in 0.008 to 0.016 % of uint16 voxels
(single: relative L2 error 2e-6 to 3e-6), measured on an A100.

``deskew_volume`` and ``rotate_volume`` replace ``deskewFrame3D`` and ``rotateFrame3D``, another
author's GPL functions, which the port does not read: the shear of x by z and the rotation about
y with the z anisotropy are affine transforms through ``imwarp``, with the output sizes, options
and classes found from MATLAB's behaviour (oracle cases on ramps and random volumes).

Volumes are ``(z, y, x)`` C arrays; boxes ``(z0, y0, x0, z1, y1, x1)``; resample factors ``(z, y,
x)``. ``deskew_rotate_volume_matlab_layout`` takes and returns MATLAB's layout (C arrays
``(z, x, y)``, y fastest: a tiff read y-fastest) for file-level code, which then runs the deskew
kernels on that layout at the mex's speed without turning the volume around (the ``imwarp``
paths turn the input around once, and the result back).
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "deskew_rotate_volume",
    "deskew_rotate_volume_matlab_layout",
    "deskew_volume",
    "rotate_file",
    "rotate_volume",
    "skewed_space_interp_matlab",
]

_KERNEL_IDS = ("volume_warp:inputError", "feather_blending:dataTypeError")
_INTERP = "@(x) any(strcmpi(x, {'cubic', 'linear'}))"


def _volume_to_matlab(v):
    """a returned (z, y, x) volume in MATLAB's form (the MATLAB-name layer)"""
    from petakit5d.core.arrays import to_matlab

    return to_matlab(v)


def _error(cls: type, message: str, identifier: str) -> Exception:
    err = cls(message)
    err.identifier = identifier
    return err


def _check_interp_method(method: Any) -> str:
    if not isinstance(method, str) or method.lower() not in ("cubic", "linear"):
        raise _error(
            ValueError,
            "The value of 'interpMethod' is invalid. It must satisfy the function: "
            "@(x)any(strcmpi(x,{'cubic','linear'})).",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )
    return method


def _kernel_error(err: Exception) -> bool:
    """a rejection by a deskew kernel (the mex's errors), after which MATLAB's code runs"""
    return getattr(err, "identifier", None) in _KERNEL_IDS


# -- layouts -----------------------------------------------------------------------------------


def _to_zyx(vol: np.ndarray, matlab_layout: bool, num_threads: int | None = None) -> np.ndarray:
    """the (z, y, x) C array of the data (MATLAB-layout data turned around in threads)"""
    if not matlab_layout:
        return vol
    from petakit5d.core.arrays import from_matlab

    return from_matlab(vol.T, num_threads=num_threads)


def _from_zyx(out: np.ndarray, matlab_layout: bool, num_threads: int | None = None) -> np.ndarray:
    if not matlab_layout:
        return out
    from petakit5d.core.arrays import to_matlab

    return to_matlab(out, num_threads=num_threads).T


# -- MATLAB's skewed space interpolation (skewed_space_interp_defined_stepsize.m) ---------------


def _colon(a: float, d: float, b: float) -> np.ndarray:
    """MATLAB's ``a : d : b`` for doubles (d > 0): the count with MATLAB's tolerance, values
    from both ends, the middle one the mean of the ends"""
    if d <= 0 or b < a or not (math.isfinite(a) and math.isfinite(d) and math.isfinite(b)):
        return np.zeros(0)
    from petakit5d.core.dtypes import mround

    tol = 2.0 * np.finfo(np.float64).eps * max(abs(a), abs(b))
    if a == math.floor(a) and d == math.floor(d):
        n = int(math.floor((b - a) / d))
    else:
        n = int(mround(np.float64((b - a) / d)))
        if a + n * d - b > tol:
            n -= 1
    last = a + n * d
    if last - b >= -tol:
        last = b
    out = np.empty(n + 1)
    for k in range(n // 2 + 1):
        out[k] = a + k * d
        out[n - k] = last - k * d
    if n % 2 == 0:
        out[n // 2] = (a + last) / 2
    return out


def _times(a: np.ndarray, w: float) -> np.ndarray:
    """single .* double: computed in double, rounded to single"""
    return (a.astype(np.float64) * w).astype(np.float32)


def skewed_space_interp_matlab(
    vol: np.ndarray, xstep: float, stepsize: float, *, reverse: bool = True
) -> np.ndarray:
    """``skewed_space_interp_defined_stepsize(im, xstep, stepsize, 'Reverse', reverse)``, the
    MATLAB function (the code path where the mex rejects its input), on a ``(z, y, x)`` volume.

    The slices of the grid ``1 : stepsize : nz`` (MATLAB's colon), each from the two raw slices
    around it, shifted by their x step: linear along x in each, then between the two; single
    arithmetic as MATLAB's (a single value with a double weight computed in double and
    rounded). A grid slice on a raw slice copies it, and the last slice is the last raw slice.
    ``stepsize`` 1 returns ``vol`` itself; otherwise the result is single.
    """
    v = vol if vol.ndim == 3 else vol.reshape((1,) * (3 - vol.ndim) + vol.shape)
    nz, ny, nx = v.shape
    if stepsize == 1:
        return vol
    xa = math.ceil(xstep)
    zout = _colon(1.0, float(stepsize), float(nz))
    nint = zout.size
    zs = np.floor(zout)
    out = np.zeros((nint, ny, nx), dtype=np.float32)
    if reverse:
        szm = xa - (zout - zs) * xstep + 1
        tzm = szm - (xa - xstep)
    else:
        szm = (zout - zs) * xstep + 1
        tzm = szm + (xa - xstep)
    s_mat = np.floor(szm)
    t_mat = np.floor(tzm)
    sw_mat = szm - s_mat
    tw_mat = tzm - t_mat
    zw_mat = zout - zs
    for z in range(1, nz + 1):
        if z == nz:
            out[-1] = v[z - 1]
            break
        im_s = np.zeros((ny, nx + xa), dtype=np.float32)
        im_t = np.zeros((ny, nx + xa), dtype=np.float32)
        if reverse:
            im_s[:, xa : xa + nx] = v[z - 1]
            im_t[:, :nx] = v[z]
        else:
            im_s[:, :nx] = v[z - 1]
            im_t[:, xa : xa + nx] = v[z]
        for zint in np.flatnonzero(zs == z):
            s, t = int(s_mat[zint]), int(t_mat[zint])
            sw, tw, zw = float(sw_mat[zint]), float(tw_mat[zint]), float(zw_mat[zint])
            if (reverse and sw == 0) or (not reverse and tw == 0):
                out[zint] = v[z - 1]
                continue
            im_sq = _times(im_s[:, s - 1 : s - 1 + nx], 1 - sw) + _times(im_s[:, s : s + nx], sw)
            im_tq = _times(im_t[:, t - 1 : t - 1 + nx], 1 - tw) + _times(im_t[:, t : t + nx], tw)
            out[zint] = _times(im_sq, 1 - zw) + _times(im_tq, zw)
    return out


# -- deskewRotateFrame3D --------------------------------------------------------------------------


def _imwarp(vol, geo, interp_method, box, matlab_layout, num_threads):
    from petakit5d.matlabcompat.warp import imwarp

    v = _to_zyx(vol, matlab_layout, num_threads)
    if box is None:
        out = imwarp(
            v, geo.affine, interp_method, output_size=geo.out_size, num_threads=num_threads
        )
    else:
        out = imwarp(v, geo.affine, interp_method, output_box=box, num_threads=num_threads)
    return _from_zyx(out, matlab_layout, num_threads)


def _interp(vol, x_step, step, reverse, save_16bit, matlab_layout, num_threads):
    """skewed_space_interp_defined_stepsize_mex, or the MATLAB function where it rejects vol"""
    from petakit5d.kernels import deskew

    try:
        if matlab_layout:
            return deskew.skewed_space_interp_defined_stepsize_matlab_layout(
                vol, x_step, step, reverse, save_16bit, num_threads=num_threads
            )
        return deskew.skewed_space_interp_defined_stepsize(
            vol, x_step, step, reverse, save_16bit, num_threads=num_threads
        )
    except (ValueError, TypeError) as err:
        if not _kernel_error(err):
            raise
    v = _to_zyx(vol, matlab_layout, num_threads)
    out = skewed_space_interp_matlab(v, x_step, step, reverse=reverse)
    return _from_zyx(out, matlab_layout, num_threads)


def _warp(vol, geo, box, save_16bit, interp_method, matlab_layout, num_threads):
    """volume_deskew_rotate_warp_mex, or imwarp where it rejects vol"""
    from petakit5d.kernels import deskew

    kbox = box if box is not None else (0, 0, 0, *geo.out_size)
    try:
        if matlab_layout:
            return deskew.volume_deskew_rotate_warp_matlab_layout(
                vol, geo.tmat, kbox, save_16bit, num_threads=num_threads
            )
        return deskew.volume_deskew_rotate_warp(
            vol, geo.tmat, kbox, save_16bit, num_threads=num_threads
        )
    except (ValueError, TypeError) as err:
        if not _kernel_error(err):
            raise
    return _imwarp(vol, geo, interp_method, box, matlab_layout, num_threads)


def _fused(vol, geo, reverse, box, save_16bit, interp_method, matlab_layout, num_threads):
    """skewed_space_interp_volume_deskew_rotate_warp_mex, or its two steps where it rejects vol"""
    from petakit5d.kernels import deskew

    kbox = box if box is not None else (0, 0, 0, *geo.out_size)
    try:
        f = (
            deskew.skewed_space_interp_volume_deskew_rotate_warp_matlab_layout
            if matlab_layout
            else deskew.skewed_space_interp_volume_deskew_rotate_warp
        )
        return f(
            vol,
            geo.x_step,
            geo.int_stepsize,
            reverse,
            geo.tmat,
            kbox,
            save_16bit,
            num_threads=num_threads,
        )
    except (ValueError, TypeError) as err:
        if not _kernel_error(err):
            raise
    vol_1 = _interp(
        vol, geo.x_step, geo.int_stepsize, reverse, save_16bit, matlab_layout, num_threads
    )
    return _warp(vol_1, geo, box, save_16bit, interp_method, matlab_layout, num_threads)


def _shape_zyx(vol: np.ndarray, matlab_layout: bool) -> tuple[int, int, int]:
    if matlab_layout:
        return vol.shape[0], vol.shape[2], vol.shape[1]
    return vol.shape


def _check(vol: Any, matlab_layout: bool) -> np.ndarray:
    if not isinstance(vol, np.ndarray):
        raise TypeError(f"vol must be a NumPy array, got {type(vol).__name__}")
    if matlab_layout:
        if vol.ndim != 3 or not vol.flags.c_contiguous:
            raise ValueError("vol must be a C-contiguous (z, x, y) array: MATLAB's layout")
        return vol
    from petakit5d.core.arrays import check_volume

    v = check_volume(vol, name="vol")
    return v if v.ndim == 3 else v[np.newaxis]


def _logical(x: Any) -> bool:
    return bool(np.asarray(x).reshape(-1)[0]) if np.size(x) else False


def _deskew_rotate(
    vol,
    angle,
    dz,
    xy_pixel_size,
    *,
    reverse,
    bbox,
    objective_scan,
    x_step_thresh,
    resample_factor,
    gpu_process,
    save_16bit,
    interp_method,
    matlab_layout,
    num_threads,
):
    from petakit5d.deskew_rotate.geometry import _check_size, deskew_rotate_geometry

    method = _check_interp_method(interp_method)
    v = _check(vol, matlab_layout)
    reverse = _logical(reverse)
    save_16bit = _logical(save_16bit)
    box = None if bbox is None or np.size(bbox) == 0 else tuple(bbox)
    kw = {
        "reverse": reverse,
        "objective_scan": _logical(objective_scan),
        "x_step_thresh": float(x_step_thresh),
        "resample_factor": resample_factor,
        "gpu_process": _logical(gpu_process),
        "interp_method": method,
    }
    shape = _shape_zyx(v, matlab_layout)
    geo = deskew_rotate_geometry(shape, float(angle), float(dz), float(xy_pixel_size), **kw)
    to_matlab_layout = False
    if matlab_layout and not geo.use_fast_method:
        # imwarp works on (z, y, x): turn the input around (smaller than an interpolated volume)
        # and the result back, rather than the interpolated volume
        v = _to_zyx(v, True, num_threads)
        matlab_layout, to_matlab_layout = False, True
    vol_1 = v
    if geo.do_interp:
        print(
            "The step size is greater than the threshold, use skewed space interpolation "
            "for combined deskew rotate...",
            flush=True,
        )
        print(f"Input dz: {float(dz):f} , interpolated dz: {geo.dz_out:f}", flush=True)
        if not geo.use_fast_method:
            import time

            t0 = time.perf_counter()
            vol_1 = _interp(
                v, geo.x_step, geo.int_stepsize, reverse, save_16bit, matlab_layout, num_threads
            )
            print(f"Skewed space interpolation time: {time.perf_counter() - t0:f} s", flush=True)
            n1 = _shape_zyx(vol_1, matlab_layout)[0]
            if n1 != geo.nz:
                geo = deskew_rotate_geometry(
                    shape, float(angle), float(dz), float(xy_pixel_size), interp_nz=n1, **kw
                )
    # MATLAB makes the output views here (imref3d of the size, and of the box): their errors
    _check_size(geo.out_size)
    if box is not None:
        _check_size([box[i + 3] - box[i] for i in range(3)])
    if not geo.use_fast_method:
        out = _imwarp(vol_1, geo, method, box, False, num_threads)
        return _from_zyx(out, to_matlab_layout, num_threads)
    if geo.do_interp:
        return _fused(v, geo, reverse, box, save_16bit, method, matlab_layout, num_threads)
    return _warp(v, geo, box, save_16bit, method, matlab_layout, num_threads)


@register(
    "deskewRotateFrame3D",
    Param("vol", None, "The raw volume (z, y, x).", convert="volume"),
    Param("angle", None, "Skew angle in degrees (typical value: 32.8)."),
    Param("dz", None, "Scan step (typical value: 0.2-0.5)."),
    Param("xyPixelSize", None, "Pixel size (typical value: 0.1)."),
    Param("reverse", "@islogical", "Reverse scan direction."),
    Param("Crop", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("bbox", "@isnumeric", "Output box (z0, y0, x0, z1, y1, x1).", convert="bbox"),
    Param("objectiveScan", "@islogical", "Objective scan (a rotation, no deskew)."),
    Param("xStepThresh", "@isnumeric", "x step above which to interpolate in skewed space."),
    Param("resampleFactor", "@isnumeric", "Resample factor (z, y, x).", convert="resample_xyz"),
    Param("gpuProcess", "@islogical", "MATLAB's GPU code path (imwarp), run on the CPU."),
    Param(
        "save16bit",
        "@islogical",
        "uint16 output for uint16 input on the kernel paths.",
        convert="same",
    ),
    Param("interpMethod", _INTERP, "Interpolation: 'linear' or 'cubic'."),
    task="D02",
    matlab_out=_volume_to_matlab,
)
def deskew_rotate_volume(
    vol: np.ndarray,
    angle: float,
    dz: float,
    xy_pixel_size: float,
    *,
    reverse: bool = False,
    crop: bool = True,
    bbox: Sequence[int] | None = None,
    objective_scan: bool = False,
    x_step_thresh: float = 2.0,
    resample_factor: Sequence[float] | None = None,
    gpu_process: bool = False,
    save_16bit: bool = False,
    interp_method: str = "linear",
) -> np.ndarray:
    """Deskew, rotate and resample a raw light-sheet volume (``deskewRotateFrame3D``).

    ``vol``: ``(z, y, x)`` (uint16 and single run the deskew kernels; other classes MATLAB's code
    paths). ``angle`` in degrees, ``dz`` the scan step and ``xy_pixel_size`` the pixel size (the
    output voxels are isotropic at that size). ``bbox``: the output region ``(z0, y0, x0, z1,
    y1, x1)`` of the full output (None: all of it). ``resample_factor``: ``(z, y, x)``, applied
    after the rotation. ``save_16bit``: uint16 output for uint16 input where the kernels run
    (single otherwise; ``imwarp`` keeps the class of what it warps). ``crop`` is accepted and not
    used, as in MATLAB. The kernels use all CPUs this process may use.

    Returns the ``(nz, ny, nx)`` output (the box's shape with ``bbox``).
    """
    _ = crop
    return _deskew_rotate(
        vol,
        angle,
        dz,
        xy_pixel_size,
        reverse=reverse,
        bbox=bbox,
        objective_scan=objective_scan,
        x_step_thresh=x_step_thresh,
        resample_factor=resample_factor,
        gpu_process=gpu_process,
        save_16bit=save_16bit,
        interp_method=interp_method,
        matlab_layout=False,
        num_threads=None,
    )


def deskew_rotate_volume_matlab_layout(
    vol_zxy: np.ndarray,
    angle: float,
    dz: float,
    xy_pixel_size: float,
    *,
    reverse: bool = False,
    bbox: Sequence[int] | None = None,
    objective_scan: bool = False,
    x_step_thresh: float = 2.0,
    resample_factor: Sequence[float] | None = None,
    gpu_process: bool = False,
    save_16bit: bool = False,
    interp_method: str = "linear",
    num_threads: int | None = None,
) -> np.ndarray:
    """``deskew_rotate_volume`` on MATLAB's layout: ``vol_zxy`` and the result are C arrays
    ``(z, x, y)`` (y fastest; MATLAB's ``(y, x, z)`` in memory, as a tiff read y-fastest gives
    it). The kernel paths run the mex's loops on that layout; the others turn the volume around
    for ``imwarp`` and back. Boxes and factors are in the Python form."""
    return _deskew_rotate(
        vol_zxy,
        angle,
        dz,
        xy_pixel_size,
        reverse=reverse,
        bbox=bbox,
        objective_scan=objective_scan,
        x_step_thresh=x_step_thresh,
        resample_factor=resample_factor,
        gpu_process=gpu_process,
        save_16bit=save_16bit,
        interp_method=interp_method,
        matlab_layout=True,
        num_threads=num_threads,
    )


# -- deskewFrame3D and rotateFrame3D: new implementations from MATLAB's behaviour ---------------


def _imwarp_view(vol, t, method, size, box):
    from petakit5d.matlabcompat.warp import imwarp

    if box is None:
        return imwarp(vol, t, method, output_size=size)
    return imwarp(vol, t, method, output_box=box)


def _as_input(vol: Any) -> np.ndarray:
    from petakit5d.core.arrays import check_volume

    return check_volume(vol, name="vol")


@register(
    "deskewFrame3D",
    Param("vol", None, "The raw volume (z, y, x).", convert="volume"),
    Param("angle", None, "Skew angle in degrees."),
    Param("dz", None, "Scan step."),
    Param("xyPixelSize", None, "Pixel size."),
    Param("reverse", "@islogical", "Reverse scan direction."),
    Param("Crop", "@islogical", "Crop half of the empty triangles in x."),
    Param("Interp", _INTERP, "Interpolation: 'linear' or 'cubic'."),
    task="D02",
    matlab_out=_volume_to_matlab,
)
def deskew_volume(
    vol: np.ndarray,
    angle: float,
    dz: float,
    xy_pixel_size: float,
    reverse: bool = False,
    *,
    crop: bool = False,
    interp: str = "linear",
) -> np.ndarray:
    """Deskew a raw light-sheet volume: the shear of x by z (``deskewFrame3D``).

    Slice ``k`` moves by ``k dx`` in x, ``dx = cos(angle) dz / xy_pixel_size`` (``-dx`` for
    reverse scans); the output is ``(nz, ny, ceil((nz - 1) dx) + nx)``, of the input's class,
    through ``imwarp`` (``interp`` 'linear' or 'cubic', 0 outside). ``crop`` takes
    ``floor(min(nx, ceil((nz - 1) dx)) / 2)`` columns off each side, as the output view of the
    warp. A new implementation from the math and MATLAB's behaviour: bit for bit MATLAB's
    ``deskewFrame3D`` in the oracle cases.
    """
    from petakit5d.deskew_rotate.geometry import deskew_geometry

    method = _check_interp_method(interp)
    v = _as_input(vol)
    shape = v.shape if v.ndim == 3 else (1, *v.shape)
    t, size, box = deskew_geometry(
        shape,
        float(angle),
        float(dz),
        float(xy_pixel_size),
        reverse=_logical(reverse),
        crop=_logical(crop),
    )
    return _imwarp_view(v, t, method, size, box)


@register(
    "rotateFrame3D",
    Param("vol", None, "The deskewed volume (z, y, x).", convert="volume"),
    Param("angle", None, "Skew angle in degrees."),
    Param("zxRatio", None, "Slice spacing in pixels (z anisotropy).", convert="same"),
    Param("reverse", "@islogical", "Reverse scan direction."),
    Param("Crop", "@islogical", "Accepted; does not change the result (as in MATLAB)."),
    Param("resample", "@isnumeric", "Resample factor (z, y, x).", convert="resample_xyz"),
    Param("ObjectiveScan", "@islogical", "Accepted; does not change the result (as in MATLAB)."),
    Param(
        "outSize", "@isnumeric", "Output size (nz, ny, nx) before resampling.", convert="size_yxz"
    ),
    Param("Interp", _INTERP, "Interpolation: 'linear' or 'cubic'."),
    task="D02",
    matlab_out=_volume_to_matlab,
)
def rotate_volume(
    vol: np.ndarray,
    angle: float,
    zx_ratio: float,
    reverse: bool = False,
    *,
    crop: bool = True,
    resample: Sequence[float] | None = None,
    objective_scan: bool = False,
    out_size: Sequence[int] | None = None,
    interp: str = "linear",
) -> np.ndarray:
    """Rotate a deskewed volume about y into the coverslip frame (``rotateFrame3D``).

    z is scaled by ``zx_ratio`` (the slice spacing in pixels), the volume rotated by ``angle``
    (``-angle`` for reverse scans) about its center and centered in the output, then resampled
    by ``resample`` (``(z, y, x)``) about the center. The output size is ``out_size`` or
    ``round(nx cos(angle) + nz zx / sin|angle|)`` in x and ``round(nz zx cos(angle) + nx
    sin|angle|)`` in z, each divided by its factor and rounded with ``resample``; the class is
    the input's (``imwarp``, ``interp`` 'linear' or 'cubic', 0 outside). ``crop`` and
    ``objective_scan`` are accepted and change nothing, as in MATLAB. A new implementation from
    the math and MATLAB's behaviour: bit for bit MATLAB's ``rotateFrame3D`` in the oracle cases.
    """
    from petakit5d.deskew_rotate.geometry import rotate_geometry

    _ = crop, objective_scan
    method = _check_interp_method(interp)
    v = _as_input(vol)
    shape = v.shape if v.ndim == 3 else (1, *v.shape)
    rs = None if resample is None or np.size(resample) == 0 else resample
    osz = None if out_size is None or np.size(out_size) == 0 else out_size
    t, size = rotate_geometry(
        shape, float(angle), float(zx_ratio), reverse=_logical(reverse), resample=rs, out_size=osz
    )
    return _imwarp_view(v, t, method, size, None)


# -- XR_RotateFrame3D ------------------------------------------------------------------------------


def _crop_box(a: np.ndarray, box) -> np.ndarray:
    """``a(b(1) : b(4), b(2) : b(5), b(3) : b(6))`` for the box ``b`` in MATLAB's values: MATLAB's
    error where the indexing leaves the array (an empty range is empty)"""
    z0, y0, x0, z1, y1, x1 = (int(v) for v in box)
    for pos, (lo, hi, n) in enumerate(
        ((y0, y1, a.shape[1]), (x0, x1, a.shape[2]), (z0, z1, a.shape[0])), start=1
    ):
        if hi > lo and (lo < 0 or hi > n):
            raise _error(
                IndexError,
                f"Index in position {pos} exceeds array bounds. Index must not exceed {n}.",
                "MATLAB:badsubscript",
            )
    sl = tuple(slice(lo, max(hi, lo)) for lo, hi in ((z0, z1), (y0, y1), (x0, x1)))
    return np.ascontiguousarray(a[sl])


@register(
    "XR_RotateFrame3D",
    Param("framePaths", "@(x) ischar(x) || iscell(x)", "Tiff file or files."),
    Param("xyPixelSize", "@isscalar", "Pixel size."),
    Param("dz", "@isscalar", "Scan step."),
    Param("objectiveScan", "@islogical", "Objective scan (zAniso dz / px)."),
    Param("Overwrite", "@islogical", "Rotate again where the result exists."),
    Param("Crop", "@islogical", "Passed to the rotation."),
    Param(
        "bbox", "@(x) isnumeric(x)", "Box (z0, y0, x0, z1, y1, x1) of the result.", convert="bbox"
    ),
    Param(
        "resample",
        "@(x) isnumeric(x)",
        "Resampling after rotation (z, y, x).",
        convert="resample_xyz",
    ),
    Param("SkewAngle", "@isscalar", "Skew angle in degrees."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("sCMOSCameraFlip", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("save16bit", "@islogical", "Save uint16 (else single).", convert="same"),
    Param("uuid", "@ischar", "Id in the temporary name (empty: a new one)."),
    task="D02",
)
def rotate_file(
    frame_paths,
    xy_pixel_size: float,
    dz: float,
    *,
    objective_scan: bool = False,
    overwrite: bool = False,
    crop: bool = True,
    bbox: Sequence[int] | None = None,
    resample: Sequence[float] | None = None,
    skew_angle: float = 31.5,
    reverse: bool = False,
    s_cmos_camera_flip: bool = False,
    save_16bit: bool = True,
    uuid: str = "",
) -> None:
    """Rotate tiff files of deskewed data (``XR_RotateFrame3D``).

    Each file is read as double, rotated (``rotate_volume`` with ``zAniso = sin(angle) dz / px``,
    ``dz / px`` for objective scans), cut to ``bbox`` and written as uint16 (or single) to
    ``<folder>/Rotated/<name>.tif``, first as ``<name>_<uuid>.tif``, then renamed. A missing
    file is skipped with a warning, an existing result unless ``overwrite``.
    """
    import warnings

    from petakit5d.core import dtypes
    from petakit5d.core.ids import get_uuid
    from petakit5d.deskew_rotate.geometry import _sin
    from petakit5d.io import read_tiff, write_tiff

    _ = s_cmos_camera_flip
    paths = [frame_paths] if isinstance(frame_paths, str) else list(frame_paths)
    uid = uuid if uuid else get_uuid()
    if _logical(objective_scan):
        z_aniso = dz / xy_pixel_size
    else:
        theta = skew_angle * math.pi / 180
        dz0 = _sin(theta) * dz
        z_aniso = dz0 / xy_pixel_size
    box = None if bbox is None or np.size(bbox) == 0 else tuple(bbox)
    for frame_path in paths:
        if not os.path.exists(frame_path):
            warnings.warn(f"{frame_path} does not exist!", stacklevel=2)
            continue
        # MATLAB's fileparts: the folder without its last separator, the name without extension
        rt, base = os.path.split(frame_path)
        fsname = os.path.splitext(base)[0]
        rt_path = rt + "/Rotated/"
        os.makedirs(rt_path, exist_ok=True)
        rt_fullname = rt_path + fsname + ".tif"
        if os.path.exists(rt_fullname) and not _logical(overwrite):
            continue
        print(f"Rotate frame {frame_path}...", flush=True)
        im = dtypes.double(read_tiff(frame_path))
        im_rt = rotate_volume(
            im,
            skew_angle,
            z_aniso,
            _logical(reverse),
            crop=crop,
            resample=resample,
            objective_scan=objective_scan,
        )
        if box is not None:
            im_rt = _crop_box(im_rt, box)
        tmp = f"{rt_path}{fsname}_{uid}.tif"
        out = dtypes.uint16(im_rt) if _logical(save_16bit) else dtypes.single(im_rt)
        write_tiff(np.ascontiguousarray(out), tmp)
        os.replace(tmp, rt_fullname)
