"""The deskew/rotate kernels (petakit5d.kernels.deskew) against the three mex files.

Oracle cases (MATLAB runs the mex, Python the port, bit for bit): the warp alone
(``volume_deskew_rotate_warp_mex``), the fused interpolation and warp
(``skewed_space_interp_volume_deskew_rotate_warp_mex``) and the interpolation alone
(``skewed_space_interp_defined_stepsize_mex``), for uint16 with and without save16bit and
single input, forward and reverse scans, step sizes around the interpolation threshold, odd
sizes and output boxes (inside, past the input's rows, entirely past them); the argument
checks of the mex (MexInputCheckTest); and the interpolation against the MATLAB function
``skewed_space_interp_defined_stepsize.m`` (DeskewRotateKernelExactTest's
interpolationMatchesTheMatlabFunction, with its tolerance).

Tests without MATLAB: a NumPy model of the mex arithmetic in its order (the port of
DeskewRotateKernelExactTest's referenceKernel), the scalar build against the AVX2 build, and
results independent of the thread count. The transforms are the sample-scan deskew and
rotation, derived below with IEEE operations only: both sides receive the same matrix.

Volumes and boxes of the cases are MATLAB's ((y, x, z) Fortran arrays, 1-based inclusive
bboxes); the harness's roles (``Vol``, ``Box``, ``out_roles=('vol',)``) give the Python side the
port's (z, y, x) C arrays and 0-based boxes and map its results back. Each case's Python side
runs both kernels (C layout on (z, y, x), the mex's loops on MATLAB's layout) and requires them
to be equal.
"""

from __future__ import annotations

import dataclasses
import math
from functools import partial

import numpy as np
import pytest

from petakit5d._ext import available, load
from petakit5d.core.arrays import from_matlab, size_from_matlab, to_matlab
from petakit5d.core.bbox import box_from_matlab
from petakit5d.core.dtypes import mround
from petakit5d.kernels import deskew
from petakit5d.testing import oracle
from petakit5d.testing.cases import Box, Case, In, Value, Vol
from petakit5d.testing.compare import Custom, MArray
from petakit5d.testing.synth import CounterRNG, stable_sincos_2pi

# ------------------------------------------------------------------------------------------
# the cases' MATLAB forms and the port's (by hand until the oracle's roles land, F14)
# ------------------------------------------------------------------------------------------


def _py_box(b):
    """a case's MATLAB bbox in the Python form; one that is not 6 values as given"""
    if isinstance(b, np.ndarray):
        v = [float(x) for x in b.reshape(-1)]
    else:
        v = list(b)
    if len(v) != 6 or not all(math.isfinite(x) and x == math.floor(x) for x in v):
        # what the mex rejects: given as it is (shifted to 0-based starts, MATLAB's order kept
        # apart from the axes), so that the port rejects it too
        return v if len(v) != 6 else (v[2] - 1, v[0] - 1, v[1] - 1, v[5], v[3], v[4])
    return box_from_matlab(v)


def _zxy(vol):
    """a (z, y, x) volume (an image: one slice) as MATLAB's layout, the (z, x, y) C array"""
    v = vol if vol.ndim == 3 else vol.reshape((1,) * (3 - vol.ndim) + vol.shape)
    return np.ascontiguousarray(v.transpose(0, 2, 1))


def _both(c_call, ml_call):
    """the port on (z, y, x) data (the C-layout kernels) and on MATLAB's layout (the mex's
    loops): both must give the same bits, which are then compared with MATLAB's"""
    out = c_call()
    if isinstance(out, np.ndarray) and out.ndim == 3:
        ml = ml_call().transpose(0, 2, 1)
        if ml.dtype != out.dtype or not np.array_equal(ml, out):
            raise AssertionError("the MATLAB-layout kernel differs from the C-layout kernel")
    return out


# the cases' Python sides: the port's (z, y, x) form in and out (the roles convert)


def p_warp(vol, tmat, box, save16bit):
    """volume_deskew_rotate_warp_mex's call through the port (both kernels)"""
    return _both(
        lambda: deskew.volume_deskew_rotate_warp(vol, tmat, box, save16bit),
        lambda: deskew.volume_deskew_rotate_warp_matlab_layout(_zxy(vol), tmat, box, save16bit),
    )


def p_fused(vol, xstep, stepsize, reverse, tmat, box, save16bit):
    """skewed_space_interp_volume_deskew_rotate_warp_mex's call through the port (both
    kernels)"""
    return _both(
        lambda: deskew.skewed_space_interp_volume_deskew_rotate_warp(
            vol, xstep, stepsize, reverse, tmat, box, save16bit
        ),
        lambda: deskew.skewed_space_interp_volume_deskew_rotate_warp_matlab_layout(
            _zxy(vol), xstep, stepsize, reverse, tmat, box, save16bit
        ),
    )


def p_interp(vol, xstep, stepsize, reverse, save16bit):
    """skewed_space_interp_defined_stepsize_mex's call through the port (both kernels)"""
    return _both(
        lambda: deskew.skewed_space_interp_defined_stepsize(
            vol, xstep, stepsize, reverse, save16bit
        ),
        lambda: deskew.skewed_space_interp_defined_stepsize_matlab_layout(
            _zxy(vol), xstep, stepsize, reverse, save16bit
        ),
    )


# the same on MATLAB's form, for the tests without MATLAB


def m_warp(vol, tmat, bbox, save16bit):
    return to_matlab(p_warp(from_matlab(vol), tmat, _py_box(bbox), save16bit))


def m_fused(vol, xstep, stepsize, reverse, tmat, bbox, save16bit):
    out = p_fused(from_matlab(vol), xstep, stepsize, reverse, tmat, _py_box(bbox), save16bit)
    return to_matlab(out)


def m_interp(vol, xstep, stepsize, reverse, save16bit):
    return to_matlab(p_interp(from_matlab(vol), xstep, stepsize, reverse, save16bit))


SKEW_ANGLE = 32.45
XY_PIXEL_SIZE = 0.108
X_STEP_THRESH = 2.0

pytestmark = pytest.mark.skipif(not available("deskew"), reason="petakit5d._ext.deskew not built")

# ------------------------------------------------------------------------------------------
# geometry of a sample-scan light-sheet stack, derived here, with IEEE operations only (the same
# bits on every platform: the matrix is an input of the cases, so both sides use the same one)
#
# Points are rows p = [x, y, z, 1] in the axes of MATLAB's array (x: dimension 2, y: dimension
# 1, z: dimension 3, in pixels), and a 4 x 4 matrix M maps p to p M (the convention of the
# kernels' tmat and of MATLAB's affine3d). Raw slice k of a sample scan is the plane z = k,
# imaged after the stage moved by dz along the scan, which is tilted by the skew angle theta
# against the sheet: in raw pixels, slice k is displaced by k * dx in x, dx = cos(theta) * dz /
# px, and the slices are sin(theta) * dz / px pixels apart along the normal of the sheet.
# ------------------------------------------------------------------------------------------


def _sincos_deg(deg: float) -> tuple[float, float]:
    s, c = stable_sincos_2pi(np.array([abs(deg) / 360.0]))
    return math.copysign(float(s[0]), deg), float(c[0])


def sum_fixed(values) -> float:
    s = 0.0
    for v in values:
        s = s + v
    return s


def _matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    """a b, each entry summed over k = 0..3 in order"""
    return [[sum_fixed(a[i][k] * b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def _inverse(m: list[list[float]]) -> list[list[float]]:
    """Gauss-Jordan elimination with partial pivoting, in a fixed order"""
    a = [list(r) + [1.0 if i == j else 0.0 for j in range(4)] for i, r in enumerate(m)]
    for c in range(4):
        p = max(range(c, 4), key=lambda r: abs(a[r][c]))
        a[c], a[p] = a[p], a[c]
        piv = a[c][c]
        a[c] = [v / piv for v in a[c]]
        for r in range(4):
            if r != c and a[r][c] != 0.0:
                f = a[r][c]
                a[r] = [a[r][j] - f * a[c][j] for j in range(8)]
    return [r[4:] for r in a]


def _identity() -> list[list[float]]:
    return [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]


def _translation(tx: float, ty: float, tz: float) -> list[list[float]]:
    """p -> p + (tx, ty, tz): the offset in the last row"""
    m = _identity()
    m[3][:3] = [tx, ty, tz]
    return m


def _interp_step(dz: float, cos_t: float) -> float:
    """the grid step (in raw slices) for an x step above the threshold: the largest step
    1 / n (n slices per raw step) with an x step at most X_STEP_THRESH, made finer while the
    interpolated dz is below 0.4 of the threshold's dz (then 1 / n is searched over fractions
    of the next denominator), or, for x steps from 1.5 thresholds on, 1 / ceil(dx / threshold)"""
    dx = cos_t * dz / XY_PIXEL_SIZE
    if abs(dx) / X_STEP_THRESH < 1.5:
        dz_thresh = XY_PIXEL_SIZE * X_STEP_THRESH / cos_t  # the dz whose dx is the threshold
        dz_out = dz / math.ceil(dz / dz_thresh)
        k = 1
        while dz_out / dz_thresh < 0.4:
            dz_out = dz / (math.ceil((dz / dz_thresh) / (1 / (k + 1))) * (1 / (k + 1)))
            k += 1
    else:
        dz_out = dz / math.ceil(abs(dx) / X_STEP_THRESH)
    return dz_out


def geometry(raw_size, dz: float, reverse: bool) -> dict:
    """the backward transform of the sample-scan deskew and rotation (no resampling) for the
    kernels, the output size, and the interpolation step"""
    ny, nx, nz = raw_size
    sin_t, cos_t = _sincos_deg(SKEW_ANGLE)
    dx = cos_t * dz / XY_PIXEL_SIZE  # x displacement per raw slice
    z_aniso = abs(sin_t) * dz / XY_PIXEL_SIZE  # slice spacing normal to the sheet, in pixels
    # output extents: x spans the deskewed slices along the coverslip, (nx - 1) cos + (nz - 1)
    # dz / px (z_aniso / sin is dz / px); z spans the slice height (nx - 1) sin, less the first
    # and last 2 planes (interpolation artifacts)
    out_size = [
        ny,
        int(mround((nx - 1) * cos_t + (nz - 1) * z_aniso / abs(sin_t))),
        int(mround((nx - 1) * abs(sin_t) - 4)),
    ]

    # x steps above the threshold interpolate in skewed space onto a finer grid of slices
    stepsize = 1.0
    xstep_raw = abs(dx)
    interp = abs(dx) > X_STEP_THRESH
    if interp:
        dz_out = _interp_step(dz, cos_t)
        stepsize = dz_out / dz
        nz = deskew.interp_grid_size(nz, stepsize)
        dz = dz_out
        dx = cos_t * dz / XY_PIXEL_SIZE
        z_aniso = abs(sin_t) * dz / XY_PIXEL_SIZE

    # 1. deskew, a shear of x by z: x' = x + shear * z + shift. Forward scans move slices to
    # +x (shear dx, the first slice at -dx so that the 1-based first slice lands at 0);
    # reverse scans to -x, offset by the width they add, ceil((nz - 1) dx). In 0-based
    # coordinates (z one less) the offset grows by one shear step.
    if not reverse:
        shift, shear = -dx, dx
    else:
        shift, shear = dx + math.ceil((nz - 1) * dx), -dx
    nx_ds = math.ceil((nz - 1) * dx) + nx  # width of the deskewed volume
    deskew_m = _identity()
    deskew_m[2][0] = shear
    deskew_m[3][0] = shift + shear
    # 2. to the center of the deskewed volume (1-based center, 0-based coordinates: + 1)
    center = [(ny + 1) / 2, (nx_ds + 1) / 2, (nz + 1) / 2]
    to_center = _translation(-center[1] + 1, -center[0] + 1, -center[2] + 1)
    # 3. slices to isotropic pixels: z scaled by the slice spacing
    scale = _identity()
    scale[2][2] = z_aniso
    # 4. rotation about y by the skew angle (its negative for reverse scans) in the x-z plane
    sin_r = -sin_t if reverse else sin_t
    rotate = [
        [cos_t, 0.0, -sin_r, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [sin_r, 0.0, cos_t, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    # 5. to the center of the output (0-based)
    from_center = _translation(
        (out_size[1] + 1) / 2 - 1, (out_size[0] + 1) / 2 - 1, (out_size[2] + 1) / 2 - 1
    )
    # the forward map raw -> output, then the backward map output -> raw. The kernels read
    # the matrix by columns (their tmat[i][j] is the backward map's (j, i) entry, the
    # transpose of the inverse) with the first two axes in array order (y, x, z): rows and
    # columns 1 and 2 swapped.
    forward = _matmul(deskew_m, _matmul(_matmul(_matmul(to_center, scale), rotate), from_center))
    backward_t = [list(col) for col in zip(*_inverse(forward), strict=True)]
    order = [1, 0, 2, 3]
    tmat = np.array([[backward_t[i][j] for j in order] for i in order], dtype=np.float64)
    return {
        "out_size": out_size,
        "tmat": tmat,
        "xstep": xstep_raw,
        "stepsize": stepsize,
        "interp": interp,
    }


# ------------------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------------------

KINDS = ("uint16", "uint16ToSingle", "single")


def make_volume(size, kind: str, seed: int = 0) -> np.ndarray:
    """the test's volumes in MATLAB's form, of MATLAB size ``size``: uint16 over the full
    range, single rand * 1000 (drawn in the port's form, which fills MATLAB's linear order, so
    the arrays and the cases' input hashes are those made before the convention change)"""
    rng = CounterRNG(seed, f"deskew/{kind}/{'x'.join(map(str, size))}")
    shape = size_from_matlab(size)
    if kind.startswith("uint16"):
        return to_matlab(rng.integers(0, 65535, shape, dtype=np.uint16))
    return to_matlab((rng.uniform(shape) * 1000.0).astype(np.float32))


def _save16(kind: str) -> bool:
    return kind != "uint16ToSingle"


def _vol(size, kind, seed=0):
    return Value(partial(make_volume, tuple(size), kind, seed))


def _mat(m):
    return Value(lambda m=m: np.array(m, dtype=np.float64))


# ------------------------------------------------------------------------------------------
# oracle cases
# ------------------------------------------------------------------------------------------

WARP = p_warp
FUSED = p_fused
INTERP = p_interp

# dz: 0.2 and 0.255 below the x step threshold (dx 1.56 and 1.99); 0.257 just above it (2.008,
# step 1/2), 0.3 (2.34, one interpolated slice per step, the default path), 0.38 (ratio 1.48,
# the search for a step), 0.39 (ratio 1.52: ndiv), 0.8 (three per step, neighbors both
# interpolated), 1.026 (a step just below 0.2)
DZ_DIRECT = (0.2, 0.255)
DZ_INTERP = (0.257, 0.3, 0.38, 0.39, 0.8, 1.026)
RAW_SIZES = {"y61": (61, 40, 30), "y48": (48, 37, 26)}


def _dz_id(dz: float) -> str:
    return f"dz{dz:g}".replace(".", "")


def _dir_id(reverse: bool) -> str:
    return "reverse" if reverse else "forward"


def _bboxes(raw, g) -> dict:
    o = g["out_size"]
    inside = [5, 7, 3, 57, o[1] - 9, o[2] - 2]
    past = [5, 7, 3, raw[0] + 9, o[1] - 9, o[2] - 2]
    beyond = [raw[0] + 2, 7, 3, raw[0] + 9, o[1] - 9, o[2] - 2]
    return {"bbox": inside, "bboxpast": past, "bboxbeyond": beyond}


def _kernel_cases() -> list[Case]:
    cases = []
    for dz in DZ_DIRECT + DZ_INTERP:
        for reverse in (True, False):
            for kind in KINDS:
                for sname, raw in RAW_SIZES.items():
                    g = geometry(raw, dz, reverse)
                    boxes = {"whole": [1, 1, 1, *g["out_size"]]}
                    if sname == "y61":
                        boxes.update(_bboxes(raw, g))
                    for bname, bbox in boxes.items():
                        kernel = "fused" if g["interp"] else "warp"
                        cid = f"{kernel}.{_dz_id(dz)}.{_dir_id(reverse)}.{kind}.{sname}.{bname}"
                        inputs = {"vol": _vol(raw, kind), "tmat": _mat(g["tmat"])}
                        if g["interp"]:
                            cases.append(
                                Case(
                                    cid,
                                    matlab="skewed_space_interp_volume_deskew_rotate_warp_mex",
                                    args=[
                                        In("vol"),
                                        g["xstep"],
                                        g["stepsize"],
                                        reverse,
                                        In("tmat"),
                                        bbox,
                                        _save16(kind),
                                    ],
                                    python=FUSED,
                                    inputs=inputs,
                                )
                            )
                        else:
                            cases.append(
                                Case(
                                    cid,
                                    matlab="volume_deskew_rotate_warp_mex",
                                    args=[In("vol"), In("tmat"), bbox, _save16(kind)],
                                    python=WARP,
                                    inputs=inputs,
                                )
                            )
    # the fused kernel's interpolated grid fed to the warp alone (deskewRotateFrame3D's two steps)
    for dz in (0.3, 0.8):
        for kind in KINDS:
            raw = RAW_SIZES["y61"]
            g = geometry(raw, dz, True)
            cases.append(
                Case(
                    f"warp.grid.{_dz_id(dz)}.{kind}",
                    code=(
                        "v = skewed_space_interp_defined_stepsize_mex(in{1}, in{2}, in{3}, true, "
                        "in{6}); out = {volume_deskew_rotate_warp_mex(v, in{4}, in{5}, in{6})};"
                    ),
                    args=[
                        In("vol"),
                        g["xstep"],
                        g["stepsize"],
                        In("tmat"),
                        [1, 1, 1, *g["out_size"]],
                        _save16(kind),
                    ],
                    python=_two_steps,
                    inputs={"vol": _vol(raw, kind), "tmat": _mat(g["tmat"])},
                )
            )
    # odd, larger sizes (lines with remainders of 8 and 16; many threads' worth of slices)
    for kind in KINDS:
        for dz, reverse in ((0.2, True), (0.3, False), (0.8, True)):
            raw = (101, 67, 53)
            g = geometry(raw, dz, reverse)
            bbox = [1, 1, 1, *g["out_size"]]
            inputs = {"vol": _vol(raw, kind, 1), "tmat": _mat(g["tmat"])}
            cid = f"odd.{_dz_id(dz)}.{_dir_id(reverse)}.{kind}"
            if g["interp"]:
                cases.append(
                    Case(
                        cid,
                        matlab="skewed_space_interp_volume_deskew_rotate_warp_mex",
                        args=[
                            In("vol"),
                            g["xstep"],
                            g["stepsize"],
                            reverse,
                            In("tmat"),
                            bbox,
                            _save16(kind),
                        ],
                        python=FUSED,
                        inputs=inputs,
                    )
                )
            else:
                cases.append(
                    Case(
                        cid,
                        matlab="volume_deskew_rotate_warp_mex",
                        args=[In("vol"), In("tmat"), bbox, _save16(kind)],
                        python=WARP,
                        inputs=inputs,
                    )
                )
    return cases


def _two_steps(vol, xstep, stepsize, tmat, box, save16bit):
    v = deskew.skewed_space_interp_defined_stepsize(vol, xstep, stepsize, True, save16bit)
    return deskew.volume_deskew_rotate_warp(v, tmat, box, save16bit)


def _interp_cases() -> list[Case]:
    cases = []
    sizes = {"y61": (61, 40, 30), "y8z100": (8, 48, 100), "y37": (37, 23, 17)}
    for dz in DZ_INTERP + (0.6,):
        for reverse in (True, False):
            for kind in KINDS:
                for sname, raw in sizes.items():
                    g = geometry(raw, dz, reverse)
                    cases.append(
                        Case(
                            f"interp.{_dz_id(dz)}.{_dir_id(reverse)}.{kind}.{sname}",
                            matlab="skewed_space_interp_defined_stepsize_mex",
                            args=[In("vol"), g["xstep"], g["stepsize"], reverse, _save16(kind)],
                            python=INTERP,
                            inputs={"vol": _vol(raw, kind)},
                        )
                    )
    # steps 1 (all raw slices), 1/2 and 1/3 exactly; xstep just below and above an integer
    for stepsize, xstep in ((1.0, 2.5), (0.5, 3.0), (1 / 3, 2.9999999), (1 / 3, 3.0000001)):
        for kind in KINDS:
            sid = f"s{stepsize:.4g}_x{xstep:.8g}".replace(".", "p")
            cases.append(
                Case(
                    f"interp.{sid}.{kind}",
                    matlab="skewed_space_interp_defined_stepsize_mex",
                    args=[In("vol"), xstep, stepsize, True, _save16(kind)],
                    python=INTERP,
                    inputs={"vol": _vol((19, 9, 7), kind)},
                )
            )
            cases.append(
                Case(
                    f"fused.{sid}.{kind}",
                    matlab="skewed_space_interp_volume_deskew_rotate_warp_mex",
                    args=[
                        In("vol"),
                        xstep,
                        stepsize,
                        False,
                        In("tmat"),
                        [1, 1, 1, 19, 13, 9],
                        _save16(kind),
                    ],
                    python=FUSED,
                    inputs={"vol": _vol((19, 9, 7), kind), "tmat": _mat(_GRID_TMAT)},
                )
            )
    return cases


# output (y, z) reads grid slice y and row z: the interpolated grid read out (the MATLAB test's
# transform for the combined mex), with a half-pixel shift so that the weights are not 0 and 1
_GRID_TMAT = [[1.0, 0, 0, 0], [0, 0, 0.5, 0.25], [0, 0.75, 0, 0.5], [0, 0, 0, 1]]


def _edge_cases() -> list[Case]:
    """2D input, a single row or slice, empty inputs, a bbox of 0 sizes, the last-row and
    last-slice nudge (ty + 1 == ny exactly)"""
    cases = []
    eye = [[1.0, 0, 0, 0], [0, 1.0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1]]
    for kind in KINDS:
        cases += [
            Case(
                f"edge.2d.{kind}",
                matlab="volume_deskew_rotate_warp_mex",
                args=[In("vol"), In("tmat"), [1, 1, 1, 17, 9, 1], _save16(kind)],
                python=WARP,
                inputs={"vol": _vol((17, 9), kind), "tmat": _mat(eye)},
            ),
            Case(
                f"edge.identity.{kind}",
                matlab="volume_deskew_rotate_warp_mex",
                args=[In("vol"), In("tmat"), [1, 1, 1, 17, 9, 5], _save16(kind)],
                python=WARP,
                inputs={"vol": _vol((17, 9, 5), kind), "tmat": _mat(eye)},
            ),
            Case(
                f"edge.lastrow.{kind}",
                matlab="volume_deskew_rotate_warp_mex",
                args=[In("vol"), In("tmat"), [3, 1, 1, 17, 4, 4], _save16(kind)],
                python=WARP,
                inputs={
                    "vol": _vol((17, 9, 5), kind),
                    # ty = 2 z + 2 hits ny - 1 = 8 at z = 3; tz = (4/3) y hits nz - 1 = 4 at y = 3
                    "tmat": _mat(
                        [[1.0, 0, 0, 0], [0, 0, 2.0, 2.0], [0, 4 / 3, 0, 0], [0, 0, 0, 1]]
                    ),
                },
            ),
            Case(
                f"edge.onerow.{kind}",
                matlab="volume_deskew_rotate_warp_mex",
                args=[In("vol"), In("tmat"), [1, 1, 1, 17, 1, 3], _save16(kind)],
                python=WARP,
                inputs={"vol": _vol((17, 1, 3), kind), "tmat": _mat(eye)},
            ),
            Case(
                f"edge.emptybox.{kind}",
                matlab="volume_deskew_rotate_warp_mex",
                args=[In("vol"), In("tmat"), [3, 2, 2, 2, 5, 4], _save16(kind)],
                python=WARP,
                inputs={"vol": _vol((17, 9, 5), kind), "tmat": _mat(eye)},
            ),
            Case(
                f"edge.emptyinput.{kind}",
                matlab="skewed_space_interp_volume_deskew_rotate_warp_mex",
                args=[In("vol"), 2.5, 0.5, True, In("tmat"), [1, 1, 1, 16, 12, 10], _save16(kind)],
                python=FUSED,
                inputs={"vol": _vol((16, 12, 0), kind), "tmat": _mat(eye)},
            ),
            Case(
                f"edge.emptyinterp.{kind}",
                matlab="skewed_space_interp_defined_stepsize_mex",
                args=[In("vol"), 2.5, 0.5, True, _save16(kind)],
                python=INTERP,
                inputs={"vol": _vol((16, 12, 0), kind)},
            ),
        ]
    return cases


# ------------------------------------------------------------------------------------------
# argument checks (MexInputCheckTest)
# ------------------------------------------------------------------------------------------

ID = "volume_warp:inputError"


def _err(cid, matlab, args, python, inputs, identifier=ID, py_error=ValueError):
    return Case(
        cid,
        matlab=matlab,
        args=args,
        python=python,
        inputs=inputs,
        expect_error=identifier,
        py_error=py_error,
        error_message=True,
    )


def _error_cases() -> list[Case]:
    eye = np.eye(4)
    bbox = [1, 1, 1, 16, 12, 10]
    vin = {"vol": _vol((16, 12, 10), "single")}
    w, f, i = (
        "volume_deskew_rotate_warp_mex",
        "skewed_space_interp_volume_deskew_rotate_warp_mex",
        "skewed_space_interp_defined_stepsize_mex",
    )
    nan41 = eye.copy()
    nan41[3, 0] = np.nan
    inf34 = eye.copy()
    inf34[2, 3] = np.inf
    v4d = {"vol": Value(lambda: np.ones((16, 12, 10, 2), dtype=np.float32, order="F"))}
    v = In("vol")
    cases = [
        # deskewRotateWarpChecksItsInputs
        _err("err.warp.4d", w, [v, eye, bbox, False], WARP, v4d),
        _err("err.warp.single_tmat", w, [v, eye.astype(np.float32), bbox, False], WARP, vin),
        _err("err.warp.int64_tmat", w, [v, eye.astype(np.int64), bbox, False], WARP, vin),
        _err("err.warp.int64_bbox", w, [v, eye, np.array(bbox, dtype=np.int64), False], None, vin),
        _err("err.warp.nan_tmat", w, [v, nan41, bbox, False], WARP, vin),
        _err("err.warp.inf_tmat", w, [v, inf34, bbox, False], WARP, vin),
        _err("err.warp.3x4_tmat", w, [v, eye[:3], bbox, False], WARP, vin),
        _err(
            "err.warp.single_bbox", w, [v, eye, np.array(bbox, dtype=np.float32), False], None, vin
        ),
        _err("err.warp.empty_bbox", w, [v, eye, np.zeros((0, 6)), False], None, vin),
        _err(
            "err.warp.column_bbox",
            w,
            [v, eye, np.array(bbox, dtype=np.float64).reshape(6, 1), False],
            None,
            vin,
        ),
        _err("err.warp.bbox5", w, [v, eye, bbox[:5], False], WARP, vin),
        _err("err.warp.nan_bbox", w, [v, eye, [np.nan, 1, 1, 16, 12, 10], False], WARP, vin),
        _err("err.warp.bbox0", w, [v, eye, [0, 1, 1, 16, 12, 10], False], WARP, vin),
        _err("err.warp.bbox_reversed", w, [v, eye, [1, 5, 1, 16, 3, 10], False], WARP, vin),
        _err(
            "err.warp.double",
            w,
            [v, eye, bbox, False],
            WARP,
            {"vol": Value(lambda: np.ones((16, 12, 10), order="F"))},
            identifier="feather_blending:dataTypeError",
            py_error=TypeError,
        ),
        # skewedSpaceInterpWarpChecksItsInputs
        _err("err.fused.4d", f, [v, 2.5, 0.5, True, eye, bbox, False], FUSED, v4d),
        _err("err.fused.xstep0", f, [v, 0, 0.5, True, eye, bbox, False], FUSED, vin),
        _err("err.fused.xstep_vector", f, [v, [2.5, 1], 0.5, True, eye, bbox, False], FUSED, vin),
        _err(
            "err.fused.single_xstep",
            f,
            [v, np.float32(2.5), 0.5, True, eye, bbox, False],
            FUSED,
            vin,
        ),
        _err("err.fused.int64_xstep", f, [v, np.int64(2), 0.5, True, eye, bbox, False], FUSED, vin),
        _err("err.fused.int64_step", f, [v, 2.5, np.int64(1), True, eye, bbox, False], FUSED, vin),
        _err(
            "err.fused.int64_tmat",
            f,
            [v, 2.5, 0.5, True, eye.astype(np.int64), bbox, False],
            FUSED,
            vin,
        ),
        _err(
            "err.fused.int64_bbox",
            f,
            [v, 2.5, 0.5, True, eye, np.array(bbox, dtype=np.int64), False],
            None,
            vin,
        ),
        _err("err.fused.step0", f, [v, 2.5, 0, True, eye, bbox, False], FUSED, vin),
        _err("err.fused.step_negative", f, [v, 2.5, -0.5, True, eye, bbox, False], FUSED, vin),
        _err("err.fused.step_nan", f, [v, 2.5, np.nan, True, eye, bbox, False], FUSED, vin),
        _err("err.fused.step_tiny", f, [v, 2.5, 1e-20, True, eye, bbox, False], FUSED, vin),
        _err("err.fused.nan_tmat", f, [v, 2.5, 0.5, True, nan41, bbox, False], FUSED, vin),
        _err(
            "err.fused.nan_bbox",
            f,
            [v, 2.5, 0.5, True, eye, [1, np.nan, 1, 16, 12, 10], False],
            FUSED,
            vin,
        ),
        _err(
            "err.fused.int32_bbox",
            f,
            [v, 2.5, 0.5, True, eye, np.array(bbox, dtype=np.int32), False],
            None,
            vin,
        ),
        _err("err.fused.narrow", f, [v, 10.5, 0.5, True, eye, bbox, False], FUSED, vin),
        _err(
            "err.fused.double",
            f,
            [v, 2.5, 0.5, True, eye, bbox, False],
            FUSED,
            {"vol": Value(lambda: np.ones((16, 12, 10), order="F"))},
            identifier="feather_blending:dataTypeError",
            py_error=TypeError,
        ),
        # skewedSpaceInterpChecksItsInputs
        _err("err.interp.4d", i, [v, 2.5, 0.5, True, False], INTERP, v4d),
        _err("err.interp.xstep0", i, [v, 0, 0.5, True, False], INTERP, vin),
        _err("err.interp.xstep_negative", i, [v, -0.5, 0.5, True, False], INTERP, vin),
        _err("err.interp.xstep_vector", i, [v, [2.5, 1], 0.5, True, False], INTERP, vin),
        _err("err.interp.xstep_underflow", i, [v, 1e-50, 0.5, True, False], INTERP, vin),
        _err("err.interp.int64_xstep", i, [v, np.int64(2), 0.5, True, False], INTERP, vin),
        _err("err.interp.int64_step", i, [v, 2.5, np.int64(1), True, False], INTERP, vin),
        _err("err.interp.step0", i, [v, 2.5, 0, True, False], INTERP, vin),
        _err("err.interp.step_inf", i, [v, 2.5, np.inf, True, False], INTERP, vin),
        _err("err.interp.step_tiny", i, [v, 2.5, 1e-20, True, False], INTERP, vin),
        _err("err.interp.narrow", i, [v, 10.5, 0.5, True, False], INTERP, vin),
        _err(
            "err.interp.double",
            i,
            [v, 2.5, 0.5, True, False],
            INTERP,
            {"vol": Value(lambda: np.ones((16, 12, 10), order="F"))},
            identifier="feather_blending:dataTypeError",
            py_error=TypeError,
        ),
    ]
    return cases


# ------------------------------------------------------------------------------------------
# the interpolation against skewed_space_interp_defined_stepsize.m
# (interpolationMatchesTheMatlabFunction)
# ------------------------------------------------------------------------------------------


def _matches_matlab_function(xstep: float, stepsize: float, is_uint16_out: bool):
    """the test's comparison: interior columns within 1e-5 of the maximum (uint16: 0.51), and
    slices at raw slice positions equal to those raw slices (checked in test_reference_model)"""

    def fn(ref: MArray, act: MArray) -> list[str]:
        if ref.data is None:
            pytest.skip("a tolerance comparison: needs the golden store (the manifest has a hash)")
        r = np.asarray(ref.data, dtype=np.float64)
        a = np.asarray(act.data, dtype=np.float64)
        if r.shape != a.shape:
            return [f"size {a.shape}, MATLAB {r.shape}"]
        tol = 0.51 if is_uint16_out else 1e-5 * float(r.max())
        cols = slice(math.ceil(xstep) + 2, r.shape[1] - math.ceil(xstep) - 2)
        d = np.abs(a[:, cols, :] - r[:, cols, :])
        m = float(d.max()) if d.size else 0.0
        return [] if m <= tol else [f"largest difference {m:.4g} > {tol:.4g}"]

    return Custom(fn, f"interior columns within the test's tolerance (step {stepsize:.4g})")


def _function_cases() -> list[Case]:
    cases = []
    raw = (8, 48, 500)
    for dz in (0.6, 0.8, 1.026):
        for reverse in (True, False):
            for kind in KINDS:
                g = geometry(raw, dz, reverse)
                nz = deskew.interp_grid_size(raw[2], g["stepsize"])
                cmp = _matches_matlab_function(g["xstep"], g["stepsize"], kind == "uint16")
                common = {
                    "matlab": "skewed_space_interp_defined_stepsize",
                    "args": [In("vol"), g["xstep"], g["stepsize"], "Reverse", reverse],
                    "inputs": {"vol": _vol(raw, kind)},
                    "compare": cmp,
                }
                cases.append(
                    Case(
                        f"mfunction.separate.{_dz_id(dz)}.{_dir_id(reverse)}.{kind}",
                        python=partial(_separate, save16bit=_save16(kind)),
                        **common,
                    )
                )
                cases.append(
                    Case(
                        f"mfunction.combined.{_dz_id(dz)}.{_dir_id(reverse)}.{kind}",
                        python=partial(_combined, save16bit=_save16(kind), nz=nz),
                        **common,
                    )
                )
    return cases


def _separate(vol, xstep, stepsize, _name, reverse, *, save16bit):
    return p_interp(vol, xstep, stepsize, reverse, save16bit)


def _combined(vol, xstep, stepsize, _name, reverse, *, save16bit, nz):
    """the fused kernel with a transform that reads grid slice y, row z: its (x, y, grid)
    output turned to the interpolated grid's (grid, y, x)"""
    tmat = np.array([[1.0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    box = (0, 0, 0, vol.shape[2], vol.shape[1], nz)  # MATLAB's [1, 1, 1, ny, nz, nx]
    out = p_fused(vol, xstep, stepsize, reverse, tmat, box, save16bit)
    return np.ascontiguousarray(out.transpose(2, 1, 0))


def _box_arg(b):
    """a case's bbox with its role; one the mex rejects for its values (not 6 of them, not
    finite, not integers) is given to Python as it is, so that the port rejects it too"""
    if isinstance(b, np.ndarray) or len(b) != 6:
        return b
    if not all(isinstance(v, (int, float)) and math.isfinite(v) and v == math.floor(v) for v in b):
        return b
    return Box(b)


# the positions of the bbox, by Python side (the volume is the first argument of all)
_ROLES = {p_warp: (2,), p_fused: (5,), p_interp: (), _two_steps: (4,)}


def _with_roles(case: Case) -> Case:
    """the case's MATLAB arguments with their roles (the key does not change)"""
    fn = case.python
    if fn is None:
        return case
    boxes = () if isinstance(fn, partial) else _ROLES[fn]
    args = list(case.args)
    args[0] = Vol(args[0])
    for i in boxes:
        args[i] = _box_arg(args[i])
    return dataclasses.replace(case, args=args, out_roles=("vol",))


CASES = [
    _with_roles(c)
    for c in _kernel_cases() + _interp_cases() + _edge_cases() + _error_cases() + _function_cases()
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# the mex arithmetic in its order, in NumPy (DeskewRotateKernelExactTest.referenceKernel)
# ------------------------------------------------------------------------------------------

F32 = np.float32


def _sum2(a, b, wa, wb, kind):
    r = a * wa + b * wb
    if kind == "uint16":
        return np.floor(r + F32(0.5)).astype(np.uint16)
    return r


def _sum4(a, b, c, d, wa, wb, wc, wd, kind):
    """single lines left to right; uint16 lines (and the uint16/single mixes) as the SIMD code
    for the first floor(n / 16) * 16 values"""
    n = a.shape[0]
    nb = n - n % 16 if kind != "single" else 0
    r = np.empty(n, dtype=F32)
    i = slice(0, nb)
    r[i] = (a[i] * wa + b[i] * wb) + (c[i] * wc + d[i] * wd)
    i = slice(nb, n)
    r[i] = ((wa * a[i] + wb * b[i]) + wc * c[i]) + wd * d[i]
    if kind == "uint16":
        return np.floor(r + F32(0.5)).astype(np.uint16)
    return r


def _grid(vol, xstep, stepsize, reverse, kind):
    """the interpolation grid: is_original, index, the interpolated slices (all rows)"""
    ny, nx, nz = vol.shape
    ntotal = deskew.interp_grid_size(nz, stepsize)
    is_orig = np.zeros(ntotal, dtype=bool)
    idx = np.zeros(ntotal, dtype=np.int64)
    nint = 0
    for i in range(ntotal):
        cs = i * stepsize
        if abs(cs - float(mround(cs))) < 1e-3:
            is_orig[i] = True
            idx[i] = int(mround(cs))
        else:
            idx[i] = nint
            nint += 1
    im_int = np.zeros((ny, nx, nint), dtype=np.uint16 if kind == "uint16" else F32, order="F")
    xa = math.ceil(xstep)

    def raw(z, r):
        return vol[:, r, z].astype(F32)

    for i in range(ntotal):
        if is_orig[i]:
            continue
        zout = stepsize * i
        zw = zout - math.floor(zout)
        if reverse:
            s, t = math.ceil(xstep) - xstep * zw, xstep * (1 - zw)
        else:
            s, t = xstep * zw, math.ceil(xstep) - xstep * (1 - zw)
        sw, s = s - math.floor(s), math.floor(s)
        tw, t = t - math.floor(t), math.floor(t)
        swf, twf, zwf = (float(F32(v)) for v in (sw, tw, zw))
        c00, c01 = F32((1 - swf) * (1 - zwf)), F32(swf * (1 - zwf))
        c10, c11 = F32((1 - twf) * zwf), F32(twf * zwf)
        zorg = math.floor(zout)
        lim = nx - t - 1
        for y in range(nx):
            if reverse:
                if y < xa - s:
                    v = _sum2(raw(zorg + 1, y + t), raw(zorg + 1, y + t + 1), c10, c11, kind)
                elif y < lim:
                    v = _sum4(
                        raw(zorg, y + s - xa),
                        raw(zorg, y + s - xa + 1),
                        raw(zorg + 1, y + t),
                        raw(zorg + 1, y + t + 1),
                        c00,
                        c01,
                        c10,
                        c11,
                        kind,
                    )
                else:
                    v = _sum2(raw(zorg, y + s - xa), raw(zorg, y + s - xa + 1), c00, c01, kind)
            else:
                if y < xa - s:
                    v = _sum2(raw(zorg, y + s), raw(zorg, y + s + 1), c00, c01, kind)
                elif y < lim:
                    v = _sum4(
                        raw(zorg, y + s),
                        raw(zorg, y + s + 1),
                        raw(zorg + 1, y + t - xa),
                        raw(zorg + 1, y + t - xa + 1),
                        c00,
                        c01,
                        c10,
                        c11,
                        kind,
                    )
                else:
                    v = _sum2(
                        raw(zorg + 1, y + t - xa), raw(zorg + 1, y + t - xa + 1), c10, c11, kind
                    )
            im_int[:, y, idx[i]] = v
    return is_orig, idx, im_int


def _kind(vol, save16bit):
    if vol.dtype != np.uint16:
        return "single"
    return "uint16" if save16bit else "mixed"


def reference_kernel(vol, tmat, bbox, save16bit, xstep=None, stepsize=1.0, reverse=False):
    """the mex's output, computed with the mex's single precision arithmetic in its order"""
    ny, nx, nz = vol.shape
    kind = _kind(vol, save16bit)
    if xstep is None:
        ntotal = nz
        is_orig, idx = np.ones(nz, dtype=bool), np.arange(nz)
        im_int = None
    else:
        is_orig, idx, im_int = _grid(vol, xstep, stepsize, reverse, kind)
        ntotal = len(is_orig)
    start = [F32(b) - F32(1) for b in bbox[:3]]
    out_shape = [int(F32(e) - s) for e, s in zip(bbox[3:], start, strict=True)]
    x0 = math.floor(float(start[0]))
    line_len = min(out_shape[0], ny - x0) if x0 < ny else 0
    xr = slice(x0, x0 + line_len)
    out = np.zeros(out_shape, dtype=np.uint16 if kind == "uint16" else F32, order="F")

    def line(zq, r):
        if is_orig[zq]:
            return vol[xr, r, idx[zq]].astype(F32)
        return im_int[xr, r, idx[zq]].astype(F32)

    for z in range(out_shape[2]):
        zf = F32(z) + start[2]
        ty = tmat[1, 2] * float(zf) + tmat[1, 3]
        if ty < 0 or ty + 1 >= nx:
            if ty + 1 == nx and nx > 1:
                counter = 1
                while math.floor(ty) + 1 == nx:
                    ty -= 1e-9 * counter
                    counter *= 2
            else:
                continue
        y0 = math.floor(ty)
        wy = F32(ty - y0)
        for y in range(out_shape[1]):
            tz = (tmat[2, 1] * float(F32(y) + start[1]) + tmat[2, 2] * float(zf)) + tmat[2, 3]
            if tz < 0 or tz + 1 >= ntotal:
                if tz + 1 == ntotal and ntotal > 1:
                    counter = 1
                    while math.floor(tz) + 1 == ntotal:
                        tz -= 1e-9 * counter
                        counter *= 2
                else:
                    continue
            z0 = math.floor(tz)
            wz = F32(tz - z0)
            one = F32(1)
            w00, w01 = (one - wy) * (one - wz), (one - wy) * wz
            w10, w11 = wy * (one - wz), wy * wz
            if line_len == 0:
                continue
            l00, l01, l10, l11 = (
                line(z0, y0),
                line(z0 + 1, y0),
                line(z0, y0 + 1),
                line(z0 + 1, y0 + 1),
            )
            if kind == "mixed":
                if is_orig[z0] and not is_orig[z0 + 1]:
                    v = _sum4(l00, l01, l10, l11, w00, w01, w10, w11, "mixed")
                elif not is_orig[z0] and is_orig[z0 + 1]:
                    v = _sum4(l01, l00, l11, l10, w01, w00, w11, w10, "mixed")
                elif is_orig[z0]:
                    v = _sum4(l00, l01, l10, l11, w00, w01, w10, w11, "mixed")
                else:
                    v = _sum4(l00, l01, l10, l11, w00, w01, w10, w11, "single")
            else:
                v = _sum4(l00, l01, l10, l11, w00, w01, w10, w11, kind)
            out[:line_len, y, z] = v
    return out


MODEL_CASES = [
    (raw, dz, reverse, kind, box)
    for raw in ((61, 40, 30), (48, 37, 26), (17, 23, 11))
    for dz in (0.2, 0.3, 0.8, 1.026)
    for reverse in (True, False)
    for kind in KINDS
    for box in ("whole", "bbox", "bboxpast")
    if box == "whole" or raw[0] == 61
]


@pytest.mark.parametrize(
    "raw, dz, reverse, kind, box",
    MODEL_CASES,
    ids=[f"{r[0]}x{r[1]}x{r[2]}-{_dz_id(d)}-{_dir_id(v)}-{k}-{b}" for r, d, v, k, b in MODEL_CASES],
)
def test_reference_model(raw, dz, reverse, kind, box):
    """the kernels equal the mex arithmetic written out in NumPy (and the grid's raw slices are
    those raw slices)"""
    vol = make_volume(raw, kind, 2)
    g = geometry(raw, dz, reverse)
    bbox = [1, 1, 1, *g["out_size"]] if box == "whole" else _bboxes(raw, g)[box]
    s16 = _save16(kind)
    if g["interp"]:
        got = m_fused(vol, g["xstep"], g["stepsize"], reverse, g["tmat"], bbox, s16)
        exp = reference_kernel(vol, g["tmat"], bbox, s16, g["xstep"], g["stepsize"], reverse)
        np.testing.assert_array_equal(got, exp)
        assert got.dtype == exp.dtype
        # the separate interpolation: xstep in single precision
        sep = m_interp(vol, g["xstep"], g["stepsize"], reverse, s16)
        xs32 = float(F32(g["xstep"]))
        is_orig, idx, im_int = _grid(vol, xs32, g["stepsize"], reverse, _kind(vol, s16))
        assert sep.shape == (raw[0], raw[1], len(is_orig))
        for zq in range(len(is_orig)):
            src = vol[:, :, idx[zq]] if is_orig[zq] else im_int[:, :, idx[zq]]
            np.testing.assert_array_equal(
                sep[:, :, zq], src.astype(sep.dtype), err_msg=f"slice {zq}"
            )
    else:
        got = m_warp(vol, g["tmat"], bbox, s16)
        exp = reference_kernel(vol, g["tmat"], bbox, s16)
        np.testing.assert_array_equal(got, exp)
        assert got.dtype == exp.dtype


# ------------------------------------------------------------------------------------------
# platforms and threads
# ------------------------------------------------------------------------------------------


def _run_all(module, vol_m, g, bbox, s16, reverse, threads):
    """the three kernels through one compiled module (the Python wrapper's checks skipped), on
    the (z, y, x) form of MATLAB's volume vol_m"""
    ext = load(module)
    vol = from_matlab(vol_m)
    start = [float(F32(b) - F32(1)) for b in bbox[:3]]  # y, x, z
    oy, ox, oz = (int(F32(e) - F32(s)) for e, s in zip(bbox[3:], start, strict=True))
    line_len = min(oy, vol.shape[1] - math.floor(start[0]))
    dt = np.uint16 if vol.dtype == np.uint16 and s16 else F32
    outs = []
    a = np.zeros((oz, oy, ox), dtype=dt)
    ext.warp(vol, g["tmat"], a, *start, line_len, threads)
    outs.append(a)
    if g["interp"]:
        b = np.zeros((oz, oy, ox), dtype=dt)
        ext.skewed_warp(
            vol, g["xstep"], g["stepsize"], reverse, g["tmat"], b, *start, line_len, threads
        )
        outs.append(b)
        n = ext.interp_grid_size(vol.shape[0], g["stepsize"])
        c = np.zeros((int(n), vol.shape[1], vol.shape[2]), dtype=dt)
        xs32 = float(F32(g["xstep"]))
        ext.interp_defined_stepsize(vol, xs32, g["stepsize"], reverse, c, threads)
        outs.append(c)
    # the kernels on MATLAB's layout: (z, x, y) C arrays, MATLAB's memory
    vm = np.ascontiguousarray(vol_m.transpose(2, 1, 0))
    a = np.zeros((oz, ox, oy), dtype=dt)
    ext.warp_ml(vm, g["tmat"], a, *start, line_len, threads)
    outs.append(a)
    if g["interp"]:
        b = np.zeros((oz, ox, oy), dtype=dt)
        ext.skewed_warp_ml(
            vm, g["xstep"], g["stepsize"], reverse, g["tmat"], b, *start, line_len, threads
        )
        outs.append(b)
        c = np.zeros((int(n), vol.shape[2], vol.shape[1]), dtype=dt)
        ext.interp_defined_stepsize_ml(vm, xs32, g["stepsize"], reverse, c, threads)
        outs.append(c)
    return outs


PLATFORM_CASES = [
    (dz, rev, kind) for dz in (0.2, 0.3, 0.8) for rev in (True, False) for kind in KINDS
]


@pytest.mark.skipif(not available("deskew_scalar"), reason="no scalar build (x86-64 only)")
@pytest.mark.parametrize("dz, reverse, kind", PLATFORM_CASES)
def test_scalar_build_equals_avx2(dz, reverse, kind):
    """the code other processors run (no AVX2) gives the bits of the AVX2 code"""
    assert load("deskew").simd() == "avx2"
    assert load("deskew_scalar").simd() == "scalar"
    # lines of 101 and 61 rows: remainders 5 and 13 of 16, 5 of 8
    for raw in ((101, 67, 53), (61, 40, 30)):
        vol = make_volume(raw, kind, 3)
        g = geometry(raw, dz, reverse)
        oy, ox, oz = g["out_size"]
        for bbox in ([1, 1, 1, oy, ox, oz], [4, 3, 2, oy, ox - 1, oz - 1]):
            a = _run_all("deskew", vol, g, bbox, _save16(kind), reverse, 4)
            b = _run_all("deskew_scalar", vol, g, bbox, _save16(kind), reverse, 4)
            for x, y in zip(a, b, strict=True):
                np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("dz, reverse, kind", PLATFORM_CASES)
def test_threads_do_not_change_results(dz, reverse, kind):
    raw = (61, 40, 30)
    vol = make_volume(raw, kind, 4)
    g = geometry(raw, dz, reverse)
    bbox = [1, 1, 1, *g["out_size"]]
    ref = _run_all("deskew", vol, g, bbox, _save16(kind), reverse, 1)
    for n in (2, 3, 8):
        for x, y in zip(
            ref, _run_all("deskew", vol, g, bbox, _save16(kind), reverse, n), strict=True
        ):
            np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("dz, reverse, kind", PLATFORM_CASES)
def test_kernel_layouts_agree(dz, reverse, kind):
    """the kernels on MATLAB's layout (the default) and on ours give the same bits, and the
    *_matlab_layout entry points equal the (z, y, x) ones"""
    for raw in ((61, 40, 30), (37, 23, 17)):
        vol = from_matlab(make_volume(raw, kind, 7))
        vm = np.ascontiguousarray(vol.transpose(0, 2, 1))
        g = geometry(raw, dz, reverse)
        oy, ox, oz = g["out_size"]
        s16 = _save16(kind)
        for box in ((0, 0, 0, oz, oy, ox), (2, 3, 1, oz - 1, oy - 2, ox - 3)):
            if g["interp"]:
                f = partial(
                    deskew.skewed_space_interp_volume_deskew_rotate_warp,
                    vol,
                    g["xstep"],
                    g["stepsize"],
                    reverse,
                    g["tmat"],
                    box,
                    s16,
                )
                fm = partial(
                    deskew.skewed_space_interp_volume_deskew_rotate_warp_matlab_layout,
                    vm,
                    g["xstep"],
                    g["stepsize"],
                    reverse,
                    g["tmat"],
                    box,
                    s16,
                )
            else:
                f = partial(deskew.volume_deskew_rotate_warp, vol, g["tmat"], box, s16)
                fm = partial(
                    deskew.volume_deskew_rotate_warp_matlab_layout, vm, g["tmat"], box, s16
                )
            a = f(_kernel_layout="matlab")
            np.testing.assert_array_equal(a, f(_kernel_layout="c"))
            np.testing.assert_array_equal(a, fm().transpose(0, 2, 1))
        if g["interp"]:
            a = deskew.skewed_space_interp_defined_stepsize(
                vol, g["xstep"], g["stepsize"], reverse, s16
            )
            c = deskew.skewed_space_interp_defined_stepsize(
                vol, g["xstep"], g["stepsize"], reverse, s16, _kernel_layout="c"
            )
            m = deskew.skewed_space_interp_defined_stepsize_matlab_layout(
                vm, g["xstep"], g["stepsize"], reverse, s16
            )
            np.testing.assert_array_equal(a, c)
            np.testing.assert_array_equal(a, m.transpose(0, 2, 1))


@pytest.mark.parametrize("dtype", (np.uint16, np.float32))
def test_swap_axes(dtype):
    """the threaded transpose of the last two axes, for sizes with and without remainders of 8"""
    ext = load("deskew")
    for shape in (
        (3, 17, 9),
        (1, 64, 128),
        (2, 1800 // 50, 512 // 16),
        (4, 8, 8),
        (0, 5, 3),
        (2, 1, 7),
    ):
        a = (CounterRNG(9).uniform(shape) * 60000).astype(dtype)
        a = np.ascontiguousarray(a)
        out = np.empty((shape[0], shape[2], shape[1]), dtype=dtype)
        ext.swap_axes(a, out, 3)
        np.testing.assert_array_equal(out, a.transpose(0, 2, 1))


# ------------------------------------------------------------------------------------------
# the Python side
# ------------------------------------------------------------------------------------------


# input row (x) depends on output z, input slice on output y: this one swaps them
_SWAP = np.array([[1.0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])


@pytest.mark.parametrize("kind", KINDS)
def test_swap_transform_copies_values(kind):
    """weights 1 and 0 reproduce the input: output (z, y, x) is input (x, y, z) away from the
    last column and slice (which are interpolated from just below)"""
    vol = from_matlab(make_volume((17, 9, 5), kind))  # (5, 17, 9)
    out = deskew.volume_deskew_rotate_warp(vol, _SWAP, (0, 0, 0, 8, 17, 4), _save16(kind))
    np.testing.assert_array_equal(out, vol[:4, :, :8].transpose(2, 1, 0).astype(out.dtype))


def test_logical_flags_as_the_mex():
    """only a logical true is true (mxIsLogicalScalarTrue): a 1 is not"""
    vol = from_matlab(make_volume((17, 9, 5), "uint16"))
    run = partial(deskew.volume_deskew_rotate_warp, vol, _SWAP, (0, 0, 0, 8, 17, 4))
    assert run(True).dtype == np.uint16
    assert run(np.True_).dtype == np.uint16
    assert run(np.array([True])).dtype == np.uint16
    assert run(1).dtype == np.float32
    assert run().dtype == np.float32
    a = deskew.skewed_space_interp_defined_stepsize(vol, 2.5, 0.5, True, True)
    b = deskew.skewed_space_interp_defined_stepsize(vol, 2.5, 0.5, 1, True)
    c = deskew.skewed_space_interp_defined_stepsize(vol, 2.5, 0.5, False, True)
    np.testing.assert_array_equal(b, c)
    assert not np.array_equal(a, c)


def test_layout_and_identifiers():
    vol = from_matlab(make_volume((17, 9, 5), "single"))
    with pytest.raises(ValueError, match="C-ordered"):
        deskew.volume_deskew_rotate_warp(np.asfortranarray(vol), np.eye(4), (0, 0, 0, 2, 2, 2))
    with pytest.raises(ValueError) as e:
        deskew.volume_deskew_rotate_warp(vol, np.eye(4), (0, 0, 0, 2, 2))
    assert e.value.identifier == "volume_warp:inputError"
    with pytest.raises(ValueError) as e:
        deskew.volume_deskew_rotate_warp(vol, np.eye(4), (0, 0, 0.5, 2, 2, 2))
    assert e.value.identifier == "volume_warp:inputError"
    with pytest.raises(TypeError) as e:
        deskew.skewed_space_interp_defined_stepsize(vol.astype(np.float64), 2.5, 0.5, True)
    assert e.value.identifier == "feather_blending:dataTypeError"
    # the output is always 3D, also for a box of one slice
    out = deskew.volume_deskew_rotate_warp(vol, np.eye(4), (0, 0, 0, 1, 17, 9))
    assert out.shape == (1, 17, 9)


def test_two_dimensional_and_row_inputs():
    """an image (y, x) is a volume of one slice; a 1-D array a row (1 x n). Interpolation needs
    two slices, so the output is all outside (0)."""
    img = from_matlab(make_volume((17, 9), "single"))
    out = deskew.volume_deskew_rotate_warp(img, _SWAP, (0, 0, 0, 8, 17, 1))
    assert out.shape == (8, 17, 1) and not out.any()
    row = np.arange(5, dtype=np.float32)
    assert deskew.volume_deskew_rotate_warp(row, _SWAP, (0, 0, 0, 1, 1, 5)).shape == (1, 1, 5)


BOX_STARTS = ((0, 0, 0), (3, 5, 2), (17, 1, 4), (8, 9, 1), (21, 0, 0))


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("dz", (0.2, 0.3, 0.8))
def test_box_offsets(kind, dz):
    """boxes that start anywhere (lines that start off a multiple of 8 or 16 keep the arithmetic
    of their positions in the line) equal the model"""
    raw = (61, 40, 30)
    vol_m = make_volume(raw, kind, 5)
    g = geometry(raw, dz, True)
    oy, ox, oz = g["out_size"]
    extra = (g["xstep"], g["stepsize"], True) if g["interp"] else ()
    for y0, x0, z0 in BOX_STARTS:
        bbox = [y0 + 1, x0 + 1, z0 + 1, oy - 1, ox - 2, oz - 1]
        if g["interp"]:
            got = m_fused(vol_m, g["xstep"], g["stepsize"], True, g["tmat"], bbox, _save16(kind))
        else:
            got = m_warp(vol_m, g["tmat"], bbox, _save16(kind))
        ref = reference_kernel(vol_m, g["tmat"], bbox, _save16(kind), *extra)
        np.testing.assert_array_equal(got, ref, err_msg=f"start {(y0, x0, z0)}")


def test_interp_grid_size():
    assert deskew.interp_grid_size(0, 0.5) == 0
    assert deskew.interp_grid_size(30, 1 / 3) == 88
    assert deskew.interp_grid_size(500, 0.2) == 2496
    # the rounding to 5 decimals absorbs the representation error of the step: (0.101 / 3) /
    # 0.101 is 0.33333333333333337, and 29 / it is just below 87
    assert deskew.interp_grid_size(11, 0.1) == 101
    assert deskew.interp_grid_size(30, (0.101 / 3) / 0.101) == 88


# steps whose rules the cases above do not reach: a grid size that needs the rounding to 5
# decimals; a slice 5e-4 from a raw slice (a copy: within 1e-3); a step that is not a short
# binary fraction
GRID_STEPS = ((0.101 / 3) / 0.101, 0.50025, 0.2873)
# xsteps whose single precision value changes the separate interpolation's weights
XSTEPS = (2.3456789123, 2.9876543211, 4.1234567891)


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("stepsize", GRID_STEPS)
def test_grid_rules_model(kind, stepsize):
    """both interpolations at these steps equal the model (the separate one with xstep in
    single precision, which matters for these xsteps)"""
    raw = (37, 23, 30)
    vol = make_volume(raw, kind, 6)
    s16 = _save16(kind)
    sensitive = False
    for xstep in XSTEPS:
        for reverse in (True, False):
            sep = m_interp(vol, xstep, stepsize, reverse, s16)
            xs32 = float(F32(xstep))
            is_orig, idx, im_int = _grid(vol, xs32, stepsize, reverse, _kind(vol, s16))
            assert sep.shape == (raw[0], raw[1], len(is_orig))
            exp = np.stack(
                [
                    vol[:, :, idx[q]] if is_orig[q] else im_int[:, :, idx[q]]
                    for q in range(len(is_orig))
                ],
                axis=2,
            ).astype(sep.dtype)
            np.testing.assert_array_equal(sep, exp)
            _, _, im64 = _grid(vol, xstep, stepsize, reverse, _kind(vol, s16))
            sensitive = sensitive or not np.array_equal(im64, im_int)
            n = len(is_orig)
            tmat = np.array([[1.0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
            bbox = [1, 1, 1, raw[0], n, raw[1]]
            got = m_fused(vol, xstep, stepsize, reverse, tmat, bbox, s16)
            ref = reference_kernel(vol, tmat, bbox, s16, xstep, stepsize, reverse)
            np.testing.assert_array_equal(got, ref)
    assert sensitive, "no xstep here where single precision changes the interpolation"
