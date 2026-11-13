"""The geometry of the deskew and rotation of light-sheet volumes, value for value.

For the combined deskew and rotation (``deskewRotateFrame3D.m``, the lab's): output sizes, the
skewed space interpolation step, and the transforms, the 4 x 4 matrix that ``imwarp`` gets
(``affine3d``'s ``T``) and the backward matrix of the deskew kernels (``tmat``). For the separate
deskew and rotation (``deskew_geometry``, ``rotate_geometry``): the same transforms on their own,
with the output sizes and crops that MATLAB's ``deskewFrame3D`` and ``rotateFrame3D`` give
(another author's GPL code, not read: found from their behaviour on 2,500 geometries each and
bit for bit on 240 random volumes each).
The matrices keep MATLAB's convention (the conventions' exception 4): a point is a row
``[x y z 1]`` of MATLAB's array axes (x the column, y the row, z the plane, in pixels) and ``p M``
maps it; entries are compared with MATLAB's.

The transforms are the textbook ones, composed as the lab's code composes them:

1. deskew, a shear of x by z: slice ``k`` moves by ``k * xstep`` in x (forward scans: ``xstep =
   dx`` and a shift of ``-dx``; reverse scans: ``xstep = -dx`` and a shift of
   ``dx + ceil((nz - 1) dx)``, the width the shear adds), ``dx = cos(theta) dz / px`` (an
   objective scan has no shear);
2. translation of the deskewed volume's center to the origin;
3. scaling of z by the slice spacing in pixels (``zAniso``);
4. rotation about y by theta in the x-z plane (``-theta`` for reverse scans);
5. translation to the center of the output, then (with resampling) to the center, the scale
   ``1 / rs`` per axis, and to the center of the resampled output.

MATLAB's arithmetic, made explicit: ``cos`` and ``sin`` are fdlibm's (``matlabcompat.elementary``,
the same bits on every platform), ``round`` is half away from zero, every expression keeps
MATLAB's order of operations, a 4 x 4 product sums each entry from +0 over k = 1..4 without
fused multiply-adds (MATLAB's, measured: 3,000 random and sparse products and 3,000 three-term
chains identical on two CPU types), and the right division ``eye(4) / M'`` is MATLAB's
(``matlabcompat.warp``: an LU of ``M'`` with partial pivoting and the transposed solve).

Sizes and boxes are in the Python form here: ``(nz, ny, nx)`` and ``(z0, y0, x0, z1, y1, x1)``;
resample factors ``(z, y, x)`` (MATLAB's ``[x y z]``).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

__all__ = [
    "DeskewRotateGeometry",
    "deskew_geometry",
    "deskew_rotate_geometry",
    "interp_step",
    "matlab_matmul",
    "output_size",
    "rotate_geometry",
]

Matrix = list[list[float]]


def _cos(x: float) -> float:
    from petakit5d.matlabcompat import elementary

    return float(np.asarray(elementary.cos(np.array([x], dtype=np.float64))).reshape(-1)[0])


def _sin(x: float) -> float:
    from petakit5d.matlabcompat import elementary

    return float(np.asarray(elementary.sin(np.array([x], dtype=np.float64))).reshape(-1)[0])


def _round(x: float) -> float:
    """MATLAB's ``round`` of a double (half away from zero)"""
    from petakit5d.core.dtypes import mround

    return float(mround(np.float64(x)))


def _ceil(x: float) -> float:
    return float(math.ceil(x)) if math.isfinite(x) else x


def _floor(x: float) -> float:
    return float(math.floor(x)) if math.isfinite(x) else x


def _eye() -> Matrix:
    return [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]


def _translation(row: Sequence[float]) -> Matrix:
    """the identity with ``[tx ty tz 1]`` as its last row"""
    m = _eye()
    m[3][:3] = [float(v) for v in row]
    return m


def matlab_matmul(a: Matrix, b: Matrix) -> Matrix:
    """MATLAB's ``a * b`` of two 4 x 4 doubles: each entry summed from +0 over k in order,
    one rounding per product and per addition (no fused multiply-add)."""
    out = []
    for i in range(4):
        row = []
        for j in range(4):
            s = 0.0
            for k in range(4):
                s = s + a[i][k] * b[k][j]
            row.append(s)
        out.append(row)
    return out


def _chain(*ms: Matrix) -> Matrix:
    """``m1 * m2 * ...``, left to right as MATLAB evaluates it"""
    out = ms[0]
    for m in ms[1:]:
        out = matlab_matmul(out, m)
    return out


def _add(a: Matrix, b: Matrix, sign: float = 1.0) -> Matrix:
    """``a + b`` (``sign`` 1) or ``a - b`` (-1), elementwise"""
    if sign > 0:
        return [[a[i][j] + b[i][j] for j in range(4)] for i in range(4)]
    return [[a[i][j] - b[i][j] for j in range(4)] for i in range(4)]


def _transpose(a: Matrix) -> Matrix:
    return [[a[j][i] for j in range(4)] for i in range(4)]


def _mrdivide(a: Matrix, t: Matrix) -> Matrix:
    from petakit5d.matlabcompat.warp import mrdivide

    return mrdivide(a, t)


def _resample_xyz(resample_factor) -> list[float] | None:
    """MATLAB's ``rs`` (``[x y z]``) from the Python (z, y, x) factor; None for ``[]``"""
    if resample_factor is None:
        return None
    v = [float(x) for x in np.asarray(resample_factor, dtype=np.float64).reshape(-1)]
    if not v:
        return None
    if len(v) != 3:
        # MATLAB indexes rs(2) and rs(3) (RS and rs([2, 1, 3])): any other length errors there
        err = IndexError(
            f"Index exceeds the number of array elements. Index must not exceed {len(v)}."
        )
        err.identifier = "MATLAB:badsubscript"
        raise err
    return v[::-1]


def output_size(
    raw_size: Sequence[int],
    angle: float,
    dz: float,
    xy_pixel_size: float,
    *,
    objective_scan: bool = False,
) -> tuple[int, int, int]:
    """The output size ``(nz, ny, nx)`` of the deskew and rotation of a raw ``(nz, ny, nx)``
    volume before resampling (``outSize``): sample scans span ``(nx - 1) cos(theta) + (nz - 1)
    dz / px`` in x and ``(nx - 1) sin|theta| - 4`` in z (the first and last two planes have
    interpolation artifacts); objective scans the rotated box exactly."""
    o = _out_size_matlab(raw_size, angle, dz, xy_pixel_size, objective_scan)
    return int(o[2]), int(o[0]), int(o[1])


def _out_size_matlab(raw_size, angle, dz, px, objective_scan) -> list[float]:
    nz, ny, nx = (float(v) for v in raw_size)
    theta = angle * math.pi / 180
    st = _sin(abs(theta))
    ct = _cos(theta)
    z_aniso = dz / px if objective_scan else st * dz / px
    if not objective_scan:
        return [
            _round(ny),
            _round((nx - 1) * ct + (nz - 1) * z_aniso / st),
            _round((nx - 1) * st - 4),
        ]
    return [_round(ny), _round(nx * ct + nz * z_aniso * st), _round(nz * z_aniso * ct + nx * st)]


def interp_step(dz: float, dx: float, theta: float, xy_pixel_size: float, x_step_thresh: float):
    """``(dzout, int_stepsize)`` of the skewed space interpolation for an x step ``dx`` above the
    threshold: below 1.5 thresholds the largest ``dz / ceil(dz / dzout_thresh)``, refined by
    fractions ``1 / (counter + 1)`` while it is below 0.4 of the threshold's dz; from 1.5
    thresholds ``dz / ceil(|dx| / x_step_thresh)``."""
    if abs(dx) / x_step_thresh < 1.5:
        dzout_thresh = xy_pixel_size * x_step_thresh / _cos(theta)
        dzout = dz / (_ceil((dz / dzout_thresh) / (1 / 1)) * (1 / 1))
        counter = 1
        while dzout / dzout_thresh < 0.4:
            dzout = dz / (_ceil((dz / dzout_thresh) / (1 / (counter + 1))) * (1 / (counter + 1)))
            counter = counter + 1
    else:
        ndiv = _ceil(abs(dx) / x_step_thresh)
        dzout = dz / ndiv
    return dzout, dzout / dz


@dataclass(frozen=True)
class DeskewRotateGeometry:
    """What ``deskewRotateFrame3D`` computes before it warps (Python sizes and boxes; matrices
    in MATLAB's convention, as ``np.ndarray`` 4 x 4 float64).

    - ``out_size``: ``(nz, ny, nx)`` of the full output (after resampling); ``out_size_matlab``
      MATLAB's ``outSize`` as it computes it (``[ny nx nz]`` doubles: -0 where a size rounds
      to zero from below);
    - ``do_interp``: skewed space interpolation first (sample scan, ``|dx| > xStepThresh``);
      ``dz_out``, ``int_stepsize`` its dz and its step in raw slices (dz and 1 without);
      ``x_step``: ``|dx|`` of the raw slices (the kernels' ``xstep``);
    - ``use_fast_method``: the deskew kernels (sample scan, linear, no resampling, CPU), else
      ``imwarp``;
    - ``nz``: the slices the transforms are built for (interpolated slices with ``do_interp``);
    - ``affine``: ``ds_S * (T1 * S * R * T2) * (RT1 * RS * RT2)``, ``affine3d``'s ``T`` for
      ``imwarp``;
    - ``tmat``: the kernels' backward matrix (0-based subscripts: the shear shifted by one
      step, the centers by one; ``eye(4) / M'`` with the first two rows and columns swapped).
    """

    raw_size: tuple[int, int, int]
    out_size: tuple[int, int, int]
    out_size_matlab: tuple[float, float, float]
    theta: float
    dx: float
    z_aniso: float
    x_step: float
    do_interp: bool
    use_fast_method: bool
    dz_out: float
    int_stepsize: float
    nz: int
    nx_ds: float
    xshift: float
    xstep: float
    affine: np.ndarray
    tmat: np.ndarray


def deskew_rotate_geometry(
    raw_size: Sequence[int],
    angle: float,
    dz: float,
    xy_pixel_size: float,
    *,
    reverse: bool = False,
    objective_scan: bool = False,
    x_step_thresh: float = 2.0,
    resample_factor=None,
    gpu_process: bool = False,
    interp_method: str = "linear",
    interp_nz: int | None = None,
) -> DeskewRotateGeometry:
    """The geometry of ``deskewRotateFrame3D(vol, angle, dz, xyPixelSize, ...)`` for a raw
    volume of size ``raw_size`` ``(nz, ny, nx)``.

    ``resample_factor``: ``(z, y, x)`` (MATLAB's ``[x y z]``), None for none. ``interp_nz``: the
    slices of the interpolated volume where the code interpolates first and then warps with
    ``imwarp`` (MATLAB takes them from the volume: the mex gives the kernels' grid size, the
    MATLAB function the count of ``1 : stepsize : nz``); None: the kernels' grid size,
    ``floor(round((nz - 1) / stepsize * 100000) / 100000) + 1``.
    """
    nz, ny, nx = (int(v) for v in raw_size)
    rs = _resample_xyz(resample_factor)
    theta = angle * math.pi / 180
    dx = _cos(theta) * dz / xy_pixel_size
    if objective_scan:
        z_aniso = dz / xy_pixel_size
    else:
        z_aniso = _sin(abs(theta)) * dz / xy_pixel_size
    out_size = _out_size_matlab((nz, ny, nx), angle, dz, xy_pixel_size, objective_scan)

    do_interp = (not objective_scan) and abs(dx) > x_step_thresh
    use_fast = not (
        objective_scan
        or gpu_process
        or (rs is not None and any(r != 1 for r in rs))
        or interp_method.lower() != "linear"
    )
    x_step = abs(dx)
    dz_out, step = dz, 1.0
    nzf = float(nz)
    if do_interp:
        dz_out, step = interp_step(dz, dx, theta, xy_pixel_size, x_step_thresh)
        if use_fast or interp_nz is None:
            nzf = _floor(_round((nzf - 1) / step * 100000) / 100000) + 1
        else:
            nzf = float(interp_nz)
        dz = dz_out
        dx = _cos(theta) * dz / xy_pixel_size
        z_aniso = _sin(abs(theta)) * dz / xy_pixel_size

    # 1. deskew
    if not reverse:
        xshift = -dx
        xstep = dx
    else:
        xshift = dx + _ceil((nzf - 1) * dx)
        xstep = -dx
    nx_ds = _ceil((nzf - 1) * dx) + nx
    if objective_scan:
        nx_ds = float(nx)
        ds_s = _eye()
    else:
        ds_s = _eye()
        ds_s[2][0] = xstep
        ds_s[3][0] = xshift

    # 2-5. rotate
    th = -theta if reverse else theta
    center = [(ny + 1) / 2, (nx_ds + 1) / 2, (nzf + 1) / 2]
    t1 = _translation([-center[1], -center[0], -center[2]])
    s = _eye()
    s[2][2] = z_aniso
    c, sn = _cos(th), _sin(th)
    r = [[c, 0.0, -sn, 0.0], [0.0, 1.0, 0.0, 0.0], [sn, 0.0, c, 0.0], [0.0, 0.0, 0.0, 1.0]]
    t2 = _translation([(out_size[1] + 1) / 2, (out_size[0] + 1) / 2, (out_size[2] + 1) / 2])

    # resampling after deskew and rotate
    if rs is not None:
        rt1 = _translation([-(out_size[1] + 1) / 2, -(out_size[0] + 1) / 2, -(out_size[2] + 1) / 2])
        rsc = _eye()
        for i in range(3):
            rsc[i][i] = 1 / rs[i]
        ryxz = (rs[1], rs[0], rs[2])
        out_size = [_round(out_size[i] / ryxz[i]) for i in range(3)]
        rt2 = _translation([(out_size[1] + 1) / 2, (out_size[0] + 1) / 2, (out_size[2] + 1) / 2])
    else:
        rt1 = rsc = rt2 = _eye()

    resample = _chain(rt1, rsc, rt2)
    affine = _chain(ds_s, _chain(t1, s, r, t2), resample)

    # the kernels' form: backward, for 0-based subscripts
    offset = [[0.0] * 4 for _ in range(4)]
    offset[3][:3] = [1.0, 1.0, 1.0]
    ds_s_c = [list(row) for row in ds_s]
    ds_s_c[3][0] = ds_s_c[3][0] + xstep
    mc = _chain(ds_s_c, _chain(_add(t1, offset), s, r, _add(t2, offset, -1.0)), resample)
    tm = _mrdivide(_eye(), _transpose(mc))
    order = (1, 0, 2, 3)
    tmat = [[tm[i][j] for j in order] for i in order]

    return DeskewRotateGeometry(
        raw_size=(nz, ny, nx),
        out_size=(int(out_size[2]), int(out_size[0]), int(out_size[1])),
        out_size_matlab=tuple(float(v) for v in out_size),
        theta=theta,
        dx=dx,
        z_aniso=z_aniso,
        x_step=x_step,
        do_interp=do_interp,
        use_fast_method=use_fast,
        dz_out=dz_out,
        int_stepsize=step,
        nz=int(nzf),
        nx_ds=nx_ds,
        xshift=xshift,
        xstep=xstep,
        affine=np.array(affine, dtype=np.float64),
        tmat=np.array(tmat, dtype=np.float64),
    )


# -- the separate deskew and rotation (deskewFrame3D and rotateFrame3D, from their behaviour) ------


def _check_size(size: Sequence[float]) -> None:
    if any(not (v > 0) for v in size):
        err = ValueError("Expected ImageSize to be positive.")
        err.identifier = "MATLAB:imref3d:expectedPositive"
        raise err


def deskew_geometry(
    raw_size: Sequence[int],
    angle: float,
    dz: float,
    xy_pixel_size: float,
    *,
    reverse: bool = False,
    crop: bool = False,
) -> tuple[np.ndarray, tuple[int, int, int], tuple[int, int, int, int, int, int] | None]:
    """The deskew of a raw ``(nz, ny, nx)`` volume as ``deskewFrame3D`` does it: ``(T, size,
    box)``, ``T`` the shear for ``imwarp`` (the deskew step above), ``size`` ``(nz, ny, nxDs)``
    with ``nxDs = ceil((nz - 1) dx) + nx``, and ``box`` the output view with ``crop`` (None
    without): ``floor(min(nx, ceil((nz - 1) dx)) / 2)`` columns off each side in x, half of the
    empty triangles' width."""
    nz, ny, nx = (int(v) for v in raw_size)
    theta = angle * math.pi / 180
    dx = _cos(theta) * dz / xy_pixel_size
    if not reverse:
        xshift, xstep = -dx, dx
    else:
        xshift, xstep = dx + _ceil((nz - 1) * dx), -dx
    w = _ceil((nz - 1) * dx)
    nx_ds = int(w + nx)
    t = np.eye(4)
    t[2, 0] = xstep
    t[3, 0] = xshift
    _check_size((ny, nx_ds, nz))
    if not crop:
        return t, (nz, ny, nx_ds), None
    m = int(_floor(min(nx, w) / 2))
    return t, (nz, ny, nx_ds), (0, 0, m, nz, ny, nx_ds - m)


def rotate_geometry(
    size: Sequence[int],
    angle: float,
    zx_ratio: float,
    *,
    reverse: bool = False,
    resample=None,
    out_size: Sequence[int] | None = None,
) -> tuple[np.ndarray, tuple[int, int, int]]:
    """The rotation of a deskewed ``(nz, ny, nx)`` volume as ``rotateFrame3D`` does it: ``(T,
    size)``, ``T = (T1 * S * R * T2) * (RT1 * RS * RT2)`` (steps 2-5 above, z scaled by
    ``zx_ratio``, centered on the input's center) and the ``(nz, ny, nx)`` output size:
    ``out_size`` or ``round([ny, nx cos(theta) + nz zx / sin|theta|, nz zx cos(theta) + nx
    sin|theta|])`` in MATLAB's ``[y x z]``, divided by the factors and rounded with
    ``resample`` (``(z, y, x)``; MATLAB's ``[x y z]``, which must have 3 values)."""
    nz, ny, nx = (int(v) for v in size)
    rs = _resample_xyz(resample)
    theta = angle * math.pi / 180
    if out_size is None:
        c, sa = _cos(theta), _sin(abs(theta))
        o = [float(ny), _round(nx * c + nz * zx_ratio / sa), _round(nz * zx_ratio * c + nx * sa)]
    else:
        oz, oy, ox = (float(v) for v in out_size)
        o = [oy, ox, oz]
    th = -theta if reverse else theta
    center = [(ny + 1) / 2, (nx + 1) / 2, (nz + 1) / 2]
    t1 = _translation([-center[1], -center[0], -center[2]])
    s = _eye()
    s[2][2] = float(zx_ratio)
    c, sn = _cos(th), _sin(th)
    r = [[c, 0.0, -sn, 0.0], [0.0, 1.0, 0.0, 0.0], [sn, 0.0, c, 0.0], [0.0, 0.0, 0.0, 1.0]]
    t2 = _translation([(o[1] + 1) / 2, (o[0] + 1) / 2, (o[2] + 1) / 2])
    if rs is not None:
        rt1 = _translation([-(o[1] + 1) / 2, -(o[0] + 1) / 2, -(o[2] + 1) / 2])
        rsc = _eye()
        for i in range(3):
            rsc[i][i] = 1 / rs[i]
        ryxz = (rs[1], rs[0], rs[2])
        o = [_round(o[i] / ryxz[i]) for i in range(3)]
        rt2 = _translation([(o[1] + 1) / 2, (o[0] + 1) / 2, (o[2] + 1) / 2])
    else:
        rt1 = rsc = rt2 = _eye()
    t = matlab_matmul(_chain(t1, s, r, t2), _chain(rt1, rsc, rt2))
    _check_size(o)
    return np.array(t, dtype=np.float64), (int(o[2]), int(o[0]), int(o[1]))
