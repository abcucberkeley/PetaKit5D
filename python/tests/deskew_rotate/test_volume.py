"""Deskew and rotation of volumes (petakit5d.deskew_rotate.volume, .geometry) against MATLAB.

Oracle cases (MATLAB runs the oracle tree's functions, Python the port; bit for bit):

- ``deskewRotateFrame3D`` (``deskew_rotate_volume``): sample and objective scans, forward and
  reverse, dz on both sides of the x step threshold and of its 1.5x rule, other thresholds,
  boxes (inside, past the input's rows), resample factors, linear and cubic, uint16 with and
  without save16bit, single, and double and uint8 (MATLAB's code paths after the mex rejects
  them), inputs narrower than the x step (the MATLAB interpolation), argument errors; through
  the snake_case API (on both layouts: ``(z, y, x)`` and MATLAB's) and the MATLAB-name layer;
- its geometry for 1,000 random geometries (MATLAB's own lines of ``deskewRotateFrame3D`` in
  the case: the output size, the slices after interpolation and the step, ``imwarp``'s matrix
  and the kernels' ``tmat``), and the output sizes of ``deskewRotateFrame3D`` itself for 1,000
  more;
- ``deskewFrame3D`` and ``rotateFrame3D`` (``deskew_volume``, ``rotate_volume``: another
  author's GPL functions, used as black boxes only): output sizes for 500 geometries each, and
  voxels of random and bead volumes over their options and classes;
- ``XR_RotateFrame3D`` (``rotate_file``): the files in ``Rotated/``;
- ``skewed_space_interp_defined_stepsize.m`` (``skewed_space_interp_matlab``) and MATLAB's colon
  that it uses;
- GoldenOutputTest's deskew/rotate scenarios: the port on each scenario's inputs against the
  oracle tree's result and the stored reference, both bit for bit (golden store only);
- ``gpu_process``: MATLAB's gpuProcess code path (the interpolation, then imwarp) computed on
  the CPU, which the port runs (MATLAB's own GPU result differs from it in the last bits).

Tests without MATLAB: DeskewRotateFrame3DTest in intent (bead grids: fitted raw-to-output maps,
canvas origins, intensities, classes, boxes, the fall-back paths), the layouts and thread
counts, MATLAB's 4 x 4 product and the MATLAB-name layer.

Volumes and boxes of the cases are MATLAB's ((y, x, z) arrays, 1-based bboxes, ``[x y z]``
resample factors); the roles give the Python side its forms.
"""

from __future__ import annotations

import math
from functools import partial
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, size_from_matlab, to_matlab
from petakit5d.deskew_rotate import geometry as geo
from petakit5d.deskew_rotate import volume as dv
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import (
    Box,
    Case,
    In,
    InPath,
    PerAxis,
    Tree,
    Value,
    Vol,
    Work,
    mcell,
)
from petakit5d.testing.synth import CounterRNG

SKEW = 32.45
PX = 0.108

# ------------------------------------------------------------------------------------------
# inputs (MATLAB's form)
# ------------------------------------------------------------------------------------------


def make_volume(size, kind: str, seed: int = 0) -> np.ndarray:
    """a test volume of MATLAB size ``size`` ([ny nx nz]) in MATLAB's form: uint16 over the
    full range, uint8, single or double rand * 1000, or a grid of single beads on a
    background"""
    rng = CounterRNG(seed, f"dsr/{kind}/{'x'.join(map(str, size))}")
    shape = size_from_matlab(size)
    if kind == "uint16":
        return to_matlab(rng.integers(0, 65535, shape, dtype=np.uint16))
    if kind == "uint8":
        return to_matlab(rng.integers(0, 255, shape, dtype=np.uint8))
    if kind == "beads":
        nz, ny, nx = shape
        pts = [
            (z, 2 + (k * 5) % max(ny - 4, 1), x)
            for k, (z, x) in enumerate(
                (z, x) for z in range(3, nz - 2, 6) for x in range(4, nx - 3, 9)
            )
        ]
        return to_matlab(synth.beads(shape, pts, (2.0, 1.5), 1000.0, 100.0))
    u = rng.uniform(shape) * 1000.0
    return to_matlab(u.astype(np.float32 if kind == "single" else np.float64))


def _vol(size, kind, seed=0):
    return Value(partial(make_volume, tuple(size), kind, seed))


# ------------------------------------------------------------------------------------------
# the Python sides
# ------------------------------------------------------------------------------------------


def p_dsr(vol, angle, dz, px, **kw):
    """deskew_rotate_volume, required to equal the same on MATLAB's layout"""
    out = dv.deskew_rotate_volume(vol, angle, dz, px, **kw)
    if vol.ndim == 3:
        kw.pop("crop", None)
        ml = dv.deskew_rotate_volume_matlab_layout(
            np.ascontiguousarray(vol.transpose(0, 2, 1)), angle, dz, px, **kw
        )
        assert np.array_equal(ml.transpose(0, 2, 1), out), "the layouts differ"
    return out


def p_matlab(name, *args):
    """the MATLAB-name layer (MATLAB's arguments and values, the result in MATLAB's form)"""
    from petakit5d import matlab

    return getattr(matlab, name)(*args)


# ------------------------------------------------------------------------------------------
# deskewRotateFrame3D: cases
# ------------------------------------------------------------------------------------------

# dz: 0.2 and 0.255 below the threshold 2 (dx 1.56 and 1.99, the warp alone); 0.3 (2.34: one
# interpolated slice per step), 0.38 (ratio 1.48 to the threshold: the search for a step), 0.39
# (1.52: ndiv), 0.8 (6.25: 4 per step)
DZ = (0.2, 0.255, 0.3, 0.38, 0.39, 0.8)
SIZE = (21, 47, 29)  # MATLAB [ny nx nz]
KINDS = {
    "uint16": ("uint16", True),
    "uint16ToSingle": ("uint16", False),
    "single": ("single", False),
}


def _id(v) -> str:
    return f"{v:g}".replace(".", "").replace("-", "m")


def dsr_case(cid, size, kind, angle, dz, px, *, seed=0, cls=None, **opts):
    """deskewRotateFrame3D(vol, angle, dz, px, name-value options) and the port"""
    names = {
        "reverse": "reverse",
        "bbox": "bbox",
        "objective_scan": "objectiveScan",
        "x_step_thresh": "xStepThresh",
        "resample_factor": "resampleFactor",
        "save_16bit": "save16bit",
        "interp_method": "interpMethod",
        "crop": "Crop",
    }
    args = [In("vol"), angle, dz, px]
    kw = {}
    for k, v in opts.items():
        args += [names[k], v]
        if k == "bbox":
            kw[k] = Box(v)
        elif k == "resample_factor":
            kw[k] = PerAxis(v, "xyz")
        else:
            kw[k] = v
    return Case(
        cid,
        matlab="deskewRotateFrame3D",
        args=args,
        python=p_dsr,
        py_args=[Vol(In("vol")), angle, dz, px],
        py_kwargs=kw,
        out_roles=("vol",),
        inputs={"vol": _vol(size, cls or kind, seed)},
    )


def _dsr_cases() -> list[Case]:
    cases = []
    # the kernel paths: dz x direction x class
    for dz in DZ:
        for reverse in (True, False):
            for kname, (cls, s16) in KINDS.items():
                cases.append(
                    dsr_case(
                        f"dsr.{_id(dz)}.{'rev' if reverse else 'fwd'}.{kname}",
                        SIZE,
                        cls,
                        SKEW,
                        dz,
                        PX,
                        reverse=reverse,
                        save_16bit=s16,
                    )
                )
    # boxes: inside, past the input's rows (y), each with the warp alone and the fused kernel
    for dz in (0.2, 0.3):
        for kname, (cls, s16) in KINDS.items():
            for bname, bbox in (("inside", [3, 9, 2, 17, 71, 13]), ("past", [5, 9, 2, 33, 71, 13])):
                cases.append(
                    dsr_case(
                        f"dsr.box.{bname}.{_id(dz)}.{kname}",
                        SIZE,
                        cls,
                        SKEW,
                        dz,
                        PX,
                        reverse=True,
                        bbox=bbox,
                        save_16bit=s16,
                    )
                )
    # other thresholds: 2.5 (dz 0.3 below it), 1.5 (ratio 1.56: ndiv), 1.0 (dz 0.2: ratio 1.56)
    for thr, dz in ((2.5, 0.3), (1.5, 0.3), (1.0, 0.2), (1.2, 0.2)):
        cases.append(
            dsr_case(
                f"dsr.thresh{_id(thr)}.{_id(dz)}",
                SIZE,
                "single",
                SKEW,
                dz,
                PX,
                reverse=True,
                x_step_thresh=thr,
            )
        )
    # other angles and pixel sizes (forward and reverse)
    for angle, px, dz in ((31.5, 0.104, 0.25), (30.1, 0.147, 0.5), (-32.45, 0.108, 0.3)):
        for reverse in (True, False):
            cases.append(
                dsr_case(
                    f"dsr.angle{_id(angle)}.px{_id(px)}.{'rev' if reverse else 'fwd'}",
                    SIZE,
                    "uint16",
                    angle,
                    dz,
                    px,
                    reverse=reverse,
                    save_16bit=True,
                )
            )
    # resampling (imwarp), with the interpolation first at dz 0.3; [1 1 1] stays on the kernels
    for rs in ([2, 2, 2], [1.5, 1, 2], [1, 1, 1], [0.8, 1.3, 1]):
        for dz in (0.2, 0.3):
            for kname, (cls, s16) in (("uint16", ("uint16", True)), ("single", ("single", False))):
                cases.append(
                    dsr_case(
                        f"dsr.rs{''.join(_id(r) for r in rs)}.{_id(dz)}.{kname}",
                        SIZE,
                        cls,
                        SKEW,
                        dz,
                        PX,
                        reverse=True,
                        resample_factor=rs,
                        save_16bit=s16,
                    )
                )
    for dz in (0.2, 0.3):
        cases.append(
            dsr_case(
                f"dsr.rs.box.{_id(dz)}",
                SIZE,
                "single",
                SKEW,
                dz,
                PX,
                reverse=True,
                resample_factor=[2, 2, 2],
                bbox=[2, 5, 2, 9, 30, 8],
            )
        )
    # objective scan (imwarp: a rotation), linear and cubic
    for dz in (0.2, 0.5):
        for reverse in (True, False):
            for method in ("linear", "cubic"):
                for kname, (cls, s16) in (
                    ("uint16", ("uint16", True)),
                    ("single", ("single", False)),
                ):
                    cases.append(
                        dsr_case(
                            f"dsr.os.{_id(dz)}.{'rev' if reverse else 'fwd'}.{method}.{kname}",
                            SIZE,
                            cls,
                            SKEW,
                            dz,
                            PX,
                            reverse=reverse,
                            objective_scan=True,
                            interp_method=method,
                            save_16bit=s16,
                        )
                    )
    # cubic sample scans (imwarp after the interpolation at dz 0.3)
    for dz in (0.2, 0.3):
        for reverse in (True, False):
            for kname, (cls, s16) in (("uint16", ("uint16", True)), ("single", ("single", False))):
                cases.append(
                    dsr_case(
                        f"dsr.cubic.{_id(dz)}.{'rev' if reverse else 'fwd'}.{kname}",
                        SIZE,
                        cls,
                        SKEW,
                        dz,
                        PX,
                        reverse=reverse,
                        interp_method="cubic",
                        save_16bit=s16,
                    )
                )
    # classes the kernels reject: double (imwarp, or the MATLAB interpolation and the warp
    # kernel), uint8 (imwarp in uint8, or the MATLAB interpolation in single)
    for cls in ("double", "uint8"):
        for dz in (0.2, 0.3):
            for extra in ({}, {"resample_factor": [2, 2, 2]}):
                tag = ".rs" if extra else ""
                cases.append(
                    dsr_case(
                        f"dsr.class.{cls}.{_id(dz)}{tag}",
                        SIZE,
                        cls,
                        SKEW,
                        dz,
                        PX,
                        reverse=True,
                        **extra,
                    )
                )
    # inputs narrower than the x step (dz 1.5: x step 11.7, 13 columns): the kernels reject them
    for kname, (cls, s16) in KINDS.items():
        for extra in ({}, {"resample_factor": [1, 1, 2]}):
            tag = ".rs" if extra else ""
            cases.append(
                dsr_case(
                    f"dsr.narrow.{kname}{tag}",
                    (8, 13, 12),
                    cls,
                    SKEW,
                    1.5,
                    PX,
                    reverse=True,
                    save_16bit=s16,
                    **extra,
                )
            )
    # Crop is accepted and not used; one slice; odd larger sizes
    cases.append(dsr_case("dsr.cropfalse", SIZE, "single", SKEW, 0.3, PX, reverse=True, crop=False))
    for kname, (cls, s16) in KINDS.items():
        cases.append(
            dsr_case(f"dsr.oneslice.{kname}", (9, 30, 1), cls, SKEW, 0.2, PX, save_16bit=s16)
        )
        for dz, reverse in ((0.2, True), (0.3, False), (0.8, True)):
            cases.append(
                dsr_case(
                    f"dsr.odd.{_id(dz)}.{'rev' if reverse else 'fwd'}.{kname}",
                    (67, 101, 53),
                    cls,
                    SKEW,
                    dz,
                    PX,
                    seed=1,
                    reverse=reverse,
                    save_16bit=s16,
                )
            )
    # bead volumes
    for dz in (0.2, 0.3, 0.6):
        for reverse in (True, False):
            cases.append(
                dsr_case(
                    f"dsr.beads.{_id(dz)}.{'rev' if reverse else 'fwd'}",
                    (24, 64, 40),
                    "beads",
                    SKEW,
                    dz,
                    PX,
                    reverse=reverse,
                )
            )
    # the MATLAB-name layer (MATLAB's values in and out)
    for dz, extra in (
        (0.3, ["reverse", True, "bbox", [3, 9, 2, 17, 71, 13]]),
        (0.2, ["resampleFactor", [1.5, 1, 2], "save16bit", True]),
        (0.3, ["objectiveScan", True, "interpMethod", "cubic"]),
    ):
        args = [In("vol"), SKEW, dz, PX, *extra]
        cases.append(
            Case(
                f"dsr.matlabname.{_id(dz)}.{extra[0]}",
                matlab="deskewRotateFrame3D",
                args=args,
                python=partial(p_matlab, "deskewRotateFrame3D"),
                inputs={"vol": _vol(SIZE, "uint16")},
            )
        )
    return cases


def _dsr_error_cases() -> list[Case]:
    v = {"vol": _vol(SIZE, "single")}
    small = {"vol": _vol((6, 5, 9), "single")}
    return [
        Case(
            "dsr.err.interp_method",
            matlab="deskewRotateFrame3D",
            args=[In("vol"), SKEW, 0.3, PX, "interpMethod", "nearest"],
            python=p_dsr,
            py_args=[Vol(In("vol")), SKEW, 0.3, PX],
            py_kwargs={"interp_method": "nearest"},
            inputs=v,
            expect_error="MATLAB:InputParser:ArgumentFailedValidation",
            py_error=ValueError,
        ),
        Case(
            "dsr.err.resample_two",
            matlab="deskewRotateFrame3D",
            args=[In("vol"), SKEW, 0.3, PX, "resampleFactor", [2, 2]],
            python=p_dsr,
            py_args=[Vol(In("vol")), SKEW, 0.3, PX],
            py_kwargs={"resample_factor": (2, 2)},
            inputs=v,
            expect_error="MATLAB:badsubscript",
            py_error=IndexError,
        ),
        # (nx - 1) sin(theta) - 4 < 1: no output
        Case(
            "dsr.err.no_output",
            matlab="deskewRotateFrame3D",
            args=[In("vol"), SKEW, 0.2, PX],
            python=p_dsr,
            py_args=[Vol(In("vol")), SKEW, 0.2, PX],
            inputs=small,
            expect_error="MATLAB:imref3d:expectedPositive",
            py_error=ValueError,
        ),
        # (nx - 1) sin(theta) - 4 rounds to -0
        Case(
            "dsr.err.zero_output",
            matlab="deskewRotateFrame3D",
            args=[In("vol"), SKEW, 0.2, PX],
            python=p_dsr,
            py_args=[Vol(In("vol")), SKEW, 0.2, PX],
            inputs={"vol": _vol((6, 8, 9), "single")},
            expect_error="MATLAB:imref3d:expectedPositive",
            py_error=ValueError,
        ),
    ]


# ------------------------------------------------------------------------------------------
# the geometry: 1,000 random geometries, MATLAB's lines and MATLAB's function
# ------------------------------------------------------------------------------------------

ANGLES = (31.5, 32.45, 30.0, 32.8, 30.1, 45.0, -32.45, 10.0)
PXS = (0.108, 0.1, 0.104, 0.147, 0.0973)

# deskewRotateFrame3D.m, lines 50-222 of the oracle tree, for a table of geometries: rows
# [ny nx nz angle dz xyPixelSize reverse objectiveScan xStepThresh rs(1:3) hasRs]
GEOMETRY_CODE = """
G = in{1}; n = size(G, 1);
S = zeros(n, 3); N = zeros(n, 4); MM = zeros(4, 4, n); TT = zeros(4, 4, n);
for i = 1 : n
    ny = G(i, 1); nx = G(i, 2); nz = G(i, 3); angle = G(i, 4); dz = G(i, 5);
    xyPixelSize = G(i, 6); Reverse = G(i, 7) ~= 0; objectiveScan = G(i, 8) ~= 0;
    xStepThresh = G(i, 9); rs = [];
    if G(i, 13), rs = G(i, 10 : 12); end
    theta = angle * pi/180;
    dx = cos(theta)*dz/xyPixelSize;
    if objectiveScan
        zAniso = dz / xyPixelSize;
    else
        zAniso = sin(abs(theta)) * dz / xyPixelSize;
    end
    if ~objectiveScan
        outSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), ...
            (nx-1)*sin(abs(theta))-4]);
    else
        outSize = round([ny, nx*cos(theta)+nz*zAniso*sin(abs(theta)), ...
            nz*zAniso*cos(theta)+nx*sin(abs(theta))]);
    end
    do_interp = ~objectiveScan && abs(dx) > xStepThresh;
    use_fast_method = ~(objectiveScan || (~isempty(rs) && any(rs ~= 1)));
    int_stepsize = 1;
    if do_interp
        if abs(dx) / xStepThresh < 1.5
            dzout_thresh = xyPixelSize * xStepThresh / cos(theta);
            dzout = dz / (ceil((dz / dzout_thresh) / (1 / 1)) *(1 / 1));
            counter = 1;
            while dzout / dzout_thresh < 0.4
                dzout = dz / (ceil((dz / dzout_thresh) / (1 / (counter + 1))) *(1 / (counter + 1)));
                counter = counter + 1;
            end
        else
            ndiv = ceil(abs(dx) / xStepThresh);
            dzout = dz / ndiv;
        end
        int_stepsize = dzout / dz;
        nz = floor(round((nz - 1) / int_stepsize * 100000) / 100000) + 1;
        dz = dzout;
        dx = cos(theta)*dz/xyPixelSize;
        zAniso = sin(abs(theta)) * dz / xyPixelSize;
    end
    if ~Reverse
        xshift = -dx;
        xstep = dx;
    else
        xshift = dx + ceil((nz-1)*dx);
        xstep = -dx;
    end
    nxDs = ceil((nz-1)*dx) + nx;
    if objectiveScan
        nxDs = nx;
        ds_S = eye(4);
    else
        ds_S = [1 0 0 0; 0 1 0 0; xstep 0 1 0; xshift 0 0 1];
    end
    if Reverse
        theta = -theta;
    end
    center = ([ny nxDs nz]+1)/2;
    T1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; -center([2 1 3]) 1];
    S1 = [1 0 0 0; 0 1 0 0; 0 0 zAniso 0; 0 0 0 1];
    R = [cos(theta) 0 -sin(theta) 0; 0 1 0 0; sin(theta) 0 cos(theta) 0; 0 0 0 1];
    T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; (outSize([2 1 3])+1)/2 1];
    if ~isempty(rs)
        RT1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; -(outSize([2,1,3])+1)/2 1];
        RS = [1/rs(1) 0 0 0; 0 1/rs(2) 0 0; 0 0 1/rs(3) 0; 0 0 0 1];
        outSize = round(outSize ./ rs([2,1,3]));
        RT2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; (outSize([2,1,3])+1)/2 1];
    else
        RT1 = eye(4); RS = eye(4); RT2 = eye(4);
    end
    M = ds_S*(T1*S1*R*T2)*(RT1*RS*RT2);
    offset = zeros(4, 4);
    offset(4, 1 : 3) = 1;
    ds_S_c = ds_S;
    ds_S_c(4, 1) = ds_S_c(4, 1) + xstep;
    tmat = eye(4) / (ds_S_c*((T1+offset)*S1*R*(T2-offset))*(RT1*RS*RT2))';
    tmat = tmat([2, 1, 3, 4], [2, 1, 3, 4]);
    S(i, :) = outSize; N(i, :) = [nz, int_stepsize, do_interp, use_fast_method];
    MM(:, :, i) = M; TT(:, :, i) = tmat;
end
out = {S, N, MM, TT};
"""

SIZES_CODE = """
G = in{1}; n = size(G, 1); S = zeros(n, 3);
for i = 1 : n
    rs = [];
    if G(i, 13), rs = G(i, 10 : 12); end
    v = zeros(G(i, 1), G(i, 2), G(i, 3), 'single');
    o = deskewRotateFrame3D(v, G(i, 4), G(i, 5), G(i, 6), 'reverse', G(i, 7) ~= 0, ...
        'objectiveScan', G(i, 8) ~= 0, 'xStepThresh', G(i, 9), 'resampleFactor', rs);
    S(i, :) = [size(o, 1), size(o, 2), size(o, 3)];
end
out = {S};
"""


def geometry_table(n: int, seed: int, small: bool) -> np.ndarray:
    """n random geometries (rows as GEOMETRY_CODE reads them), from exact operations: dz,
    pixel sizes, thresholds and factors are integers divided by powers of ten"""
    rng = CounterRNG(seed, f"dsr/geometry/{n}/{small}")
    k = rng.integers(0, 1 << 30, (n, 16), dtype=np.int64)
    rows = []
    for r in k.tolist():
        if small:
            ny, nx, nz = 4 + r[0] % 3, 60 + r[1] % 101, 2 + r[2] % 59
            angle = (31.5, 32.45, 30.0, 45.0, -32.45, 20.0, 60.0)[r[3] % 7]
            dz = (50 + r[4] % 1151) / 1000
            px = (0.108, 0.1, 0.147)[r[5] % 3] if r[6] % 2 else (800 + r[5] % 1200) / 10000
        else:
            ny, nx, nz = 1 + r[0] % 2000, 2 + r[1] % 599, 1 + r[2] % 1200
            angle = ANGLES[r[3] % len(ANGLES)] if r[6] % 2 else (r[3] % 12001 - 6000) / 100
            if abs(angle) < 1:
                angle = 31.5
            dz = (50 + r[4] % 1551) / 1000
            px = PXS[r[5] % len(PXS)] if r[7] % 10 < 7 else (500 + r[5] % 2001) / 10000
        rev = r[8] % 2
        osc = 1 if r[9] % 5 == 0 else 0
        thr = 2.0 if r[10] % 10 < 7 else (50 + r[10] % 301) / 100
        u = r[11] % 10
        if u < 5:
            rs, has = (0.0, 0.0, 0.0), 0
        elif u < 6:
            rs, has = (1.0, 1.0, 1.0), 1
        elif u < 8:
            rs, has = tuple(float(1 + v % 3) for v in r[12:15]), 1
        else:
            rs, has = tuple((50 + v % 251) / 100 for v in r[12:15]), 1
        rows.append([ny, nx, nz, angle, dz, px, rev, osc, thr, *rs, has])
    return np.array(rows, dtype=np.float64)


def _geometry_of_row(row) -> geo.DeskewRotateGeometry:
    ny, nx, nz, angle, dz, px, rev, osc, thr, r1, r2, r3, has = (float(v) for v in row)
    return geo.deskew_rotate_geometry(
        (int(nz), int(ny), int(nx)),
        angle,
        dz,
        px,
        reverse=bool(rev),
        objective_scan=bool(osc),
        x_step_thresh=thr,
        resample_factor=(r3, r2, r1) if has else None,
    )


def p_geometry(g: np.ndarray) -> list:
    n = g.shape[0]
    s = np.zeros((n, 3))
    nn = np.zeros((n, 4))
    mm = np.zeros((4, 4, n))
    tt = np.zeros((4, 4, n))
    for i in range(n):
        x = _geometry_of_row(g[i])
        s[i] = x.out_size_matlab
        nn[i] = (x.nz, x.int_stepsize, x.do_interp, x.use_fast_method)
        mm[:, :, i] = x.affine
        tt[:, :, i] = x.tmat
    return [s, nn, mm, tt]


def p_sizes(g: np.ndarray) -> np.ndarray:
    out = np.zeros((g.shape[0], 3))
    for i in range(g.shape[0]):
        x = _geometry_of_row(g[i])
        out[i] = (x.out_size[1], x.out_size[2], x.out_size[0])
    return out


def _geometry_cases() -> list[Case]:
    return [
        Case(
            "geometry.matrices",
            code=GEOMETRY_CODE,
            args=[In("g")],
            nargout=4,
            outputs=["outSize", "nzStep", "affine", "tmat"],
            python=p_geometry,
            inputs={"g": Value(partial(geometry_table, 1000, 0, False))},
        ),
        Case(
            "geometry.sizes",
            code=SIZES_CODE,
            args=[In("g")],
            python=p_sizes,
            inputs={"g": Value(partial(geometry_table, 1000, 1, True))},
        ),
    ]


# ------------------------------------------------------------------------------------------
# deskewFrame3D and rotateFrame3D (black boxes)
# ------------------------------------------------------------------------------------------

DESKEW_SIZES_CODE = """
G = in{1}; n = size(G, 1); S = zeros(n, 3);
for i = 1 : n
    v = zeros(G(i, 1), G(i, 2), G(i, 3), 'single');
    o = deskewFrame3D(v, G(i, 4), G(i, 5), G(i, 6), G(i, 7) ~= 0, 'Crop', G(i, 8) ~= 0);
    S(i, :) = [size(o, 1), size(o, 2), size(o, 3)];
end
out = {S};
"""

ROTATE_SIZES_CODE = """
G = in{1}; n = size(G, 1); S = zeros(n, 3);
for i = 1 : n
    rs = []; os = [];
    if G(i, 12), rs = G(i, 9 : 11); end
    if G(i, 16), os = G(i, 13 : 15); end
    v = zeros(G(i, 1), G(i, 2), G(i, 3), 'single');
    o = rotateFrame3D(v, G(i, 4), G(i, 5), G(i, 6) ~= 0, 'Crop', G(i, 7) ~= 0, ...
        'ObjectiveScan', G(i, 8) ~= 0, 'resample', rs, 'outSize', os);
    S(i, :) = [size(o, 1), size(o, 2), size(o, 3)];
end
out = {S};
"""


def deskew_size_table(n: int) -> np.ndarray:
    """rows [ny nx nz angle dz px reverse crop]"""
    rng = CounterRNG(2, f"dsr/deskew_sizes/{n}")
    rows = []
    for r in rng.integers(0, 1 << 30, (n, 9), dtype=np.int64).tolist():
        angle = (
            (31.5, 32.45, -32.45, 45.0, 10.0, 60.0)[r[3] % 6]
            if r[8] % 2
            else ((r[3] % 14001 - 7000) / 100)
        )
        rows.append(
            [
                1 + r[0] % 4,
                2 + r[1] % 199,
                2 + r[2] % 149,
                angle,
                (50 + r[4] % 1201) / 1000,
                (500 + r[5] % 2001) / 10000,
                r[6] % 2,
                r[7] % 2,
            ]
        )
    return np.array(rows, dtype=np.float64)


def rotate_size_table(n: int) -> np.ndarray:
    """rows [ny nx nz angle zx reverse crop os rs(1:3) hasRs outSize(1:3) hasOut]"""
    rng = CounterRNG(3, f"dsr/rotate_sizes/{n}")
    rows = []
    for r in rng.integers(0, 1 << 30, (n, 16), dtype=np.int64).tolist():
        ny = 3 + r[0] % 4
        angle = (
            (31.5, 32.45, -32.45, 45.0, 20.0, 60.0)[r[3] % 6]
            if r[15] % 2
            else ((15 + r[3] % 46) * (1 if r[3] % 2 else -1))
        )
        u = r[8] % 10
        if u < 5:
            rs, has = (0.0, 0.0, 0.0), 0
        elif u < 7:
            rs, has = (1.5, 1.5, 1.5), 1
        else:
            rs, has = tuple((60 + v % 191) / 100 for v in r[9:12]), 1
        if r[12] % 7 == 0:
            osz, has_o = (ny, 10 + r[13] % 300, 5 + r[14] % 200), 1
        else:
            osz, has_o = (0, 0, 0), 0
        rows.append(
            [
                ny,
                2 + r[1] % 119,
                2 + r[2] % 79,
                angle,
                (20 + r[4] % 281) / 100,
                r[5] % 2,
                r[6] % 2,
                1 if r[7] % 3 == 0 else 0,
                *rs,
                has,
                *osz,
                has_o,
            ]
        )
    return np.array(rows, dtype=np.float64)


def p_deskew_sizes(g: np.ndarray) -> np.ndarray:
    out = np.zeros((g.shape[0], 3))
    for i, (ny, nx, nz, angle, dz, px, rev, crop) in enumerate(g.tolist()):
        _, size, box = geo.deskew_geometry(
            (int(nz), int(ny), int(nx)), angle, dz, px, reverse=bool(rev), crop=bool(crop)
        )
        nx_out = size[2] if box is None else box[5] - box[2]
        out[i] = (size[1], nx_out, size[0])
    return out


def p_rotate_sizes(g: np.ndarray) -> np.ndarray:
    out = np.zeros((g.shape[0], 3))
    for i, r in enumerate(g.tolist()):
        ny, nx, nz, angle, zx = r[:5]
        rs = (r[10], r[9], r[8]) if r[11] else None
        osz = (int(r[14]), int(r[12]), int(r[13])) if r[15] else None
        _, size = geo.rotate_geometry(
            (int(nz), int(ny), int(nx)), angle, zx, reverse=bool(r[5]), resample=rs, out_size=osz
        )
        out[i] = (size[1], size[2], size[0])
    return out


def _deskew_rotate_frame_cases() -> list[Case]:
    cases = [
        Case(
            "deskew.sizes",
            code=DESKEW_SIZES_CODE,
            args=[In("g")],
            python=p_deskew_sizes,
            inputs={"g": Value(partial(deskew_size_table, 500))},
        ),
        Case(
            "rotate.sizes",
            code=ROTATE_SIZES_CODE,
            args=[In("g")],
            python=p_rotate_sizes,
            inputs={"g": Value(partial(rotate_size_table, 500))},
        ),
    ]
    size = (17, 43, 23)
    # deskewFrame3D(vol, angle, dz, px, reverse, 'Crop', crop, 'Interp', method)
    for cls in ("single", "uint16", "double", "beads"):
        for angle, dz, px in ((SKEW, 0.3, PX), (-31.5, 0.5, 0.104)):
            for reverse in (True, False):
                for crop in (False, True):
                    for method in ("linear", "cubic"):
                        if cls == "beads" and method == "cubic":
                            continue
                        cid = (
                            f"deskew.{cls}.{_id(angle)}.{'rev' if reverse else 'fwd'}."
                            f"{'crop' if crop else 'full'}.{method}"
                        )
                        cases.append(
                            Case(
                                cid,
                                matlab="deskewFrame3D",
                                args=[
                                    Vol(In("v")),
                                    angle,
                                    dz,
                                    px,
                                    reverse,
                                    "Crop",
                                    crop,
                                    "Interp",
                                    method,
                                ],
                                python=dv.deskew_volume,
                                py_args=[Vol(In("v")), angle, dz, px, reverse],
                                py_kwargs={"crop": crop, "interp": method},
                                out_roles=("vol",),
                                inputs={"v": _vol(size, cls)},
                            )
                        )
    # rotateFrame3D(vol, angle, zxRatio, reverse, 'Crop', ..., 'resample', ..., 'outSize', ...)
    opts = (
        ("plain", {}),
        ("nocrop", {"Crop": False}),
        ("os", {"ObjectiveScan": True}),
        ("rs", {"resample": [1.5, 1.5, 1.5]}),
        ("rsxyz", {"resample": [1.2, 1.7, 0.8]}),
        ("out", {"outSize": [17, 60, 30]}),
        ("outrs", {"outSize": [17, 60, 30], "resample": [2, 1, 1.5]}),
    )
    for cls in ("single", "uint16", "double", "beads"):
        for angle, zx in ((SKEW, 1.4993), (-31.5, 0.8)):
            for reverse in (True, False):
                for oname, o in opts:
                    for method in ("linear", "cubic"):
                        if method == "cubic" and (cls == "beads" or oname not in ("plain", "rs")):
                            continue
                        cid = (
                            f"rotate.{cls}.{_id(angle)}.{'rev' if reverse else 'fwd'}."
                            f"{oname}.{method}"
                        )
                        margs, kw = [], {"interp": method}
                        for k, val in o.items():
                            margs += [k, val]
                            if k == "resample":
                                kw["resample"] = PerAxis(val, "xyz")
                            elif k == "outSize":
                                kw["out_size"] = PerAxis(val, "yxz")
                            elif k == "Crop":
                                kw["crop"] = val
                            else:
                                kw["objective_scan"] = val
                        cases.append(
                            Case(
                                cid,
                                matlab="rotateFrame3D",
                                args=[Vol(In("v")), angle, zx, reverse, *margs, "Interp", method],
                                python=dv.rotate_volume,
                                py_args=[Vol(In("v")), angle, zx, reverse],
                                py_kwargs=kw,
                                out_roles=("vol",),
                                inputs={"v": _vol(size, cls)},
                            )
                        )
    cases.append(
        Case(
            "rotate.err.scalar_resample",
            matlab="rotateFrame3D",
            args=[Vol(In("v")), SKEW, 1.5, False, "resample", 2],
            python=dv.rotate_volume,
            py_args=[Vol(In("v")), SKEW, 1.5, False],
            py_kwargs={"resample": (2,)},
            inputs={"v": _vol(size, "single")},
            expect_error="MATLAB:badsubscript",
            py_error=IndexError,
        )
    )
    return cases


# ------------------------------------------------------------------------------------------
# XR_RotateFrame3D, skewed_space_interp_defined_stepsize.m, MATLAB's colon
# ------------------------------------------------------------------------------------------


def frame_tree(cls: str, names=("scan_ch0",), rotated=()):
    def make(folder):
        for k, n in enumerate(names):
            a = from_matlab(make_volume((19, 37, 21), cls, k))
            synth.write_tiff(Path(folder) / f"{n}.tif", a)
        if rotated:
            (Path(folder) / "Rotated").mkdir()
            for n in rotated:
                a = np.full((3, 4, 5), 7, dtype=np.uint16)
                synth.write_tiff(Path(folder) / "Rotated" / f"{n}.tif", a)

    return make


def rotate_file_case(cid, tree, frames, *args):
    """XR_RotateFrame3D(<d>/<frames>, ...) on a copy of the folder d"""
    code = (
        "f = cellfun(@(x) [in{1}, x], strsplit(in{2}, '|'), 'unif', 0); "
        "if numel(f) == 1, f = f{1}; end; XR_RotateFrame3D(f, in{3 : end}); out = {};"
    )
    return Case(
        cid,
        code=code,
        nargout=0,
        args=[InPath("d"), "|".join(frames), *args],
        python=p_rotate_file,
        inputs={"d": Tree(tree, "d")},
    )


def p_rotate_file(d, frames, *args):
    from petakit5d import matlab

    f = [d + x for x in frames.split("|")]
    matlab.XR_RotateFrame3D(f[0] if len(f) == 1 else f, *args)
    return []


def _rotate_file_cases() -> list[Case]:
    base = ["SkewAngle", SKEW, "uuid", "u0"]
    return [
        rotate_file_case(
            "rotate_file.default", frame_tree("uint16"), ["scan_ch0.tif"], PX, 0.3, *base
        ),
        rotate_file_case(
            "rotate_file.single",
            frame_tree("uint16"),
            ["scan_ch0.tif"],
            PX,
            0.3,
            *base,
            "save16bit",
            False,
        ),
        rotate_file_case(
            "rotate_file.reverse_os",
            frame_tree("single"),
            ["scan_ch0.tif"],
            PX,
            0.2,
            *base,
            "Reverse",
            True,
            "objectiveScan",
            True,
        ),
        rotate_file_case(
            "rotate_file.bbox_rs",
            frame_tree("uint16"),
            ["scan_ch0.tif"],
            PX,
            0.3,
            *base,
            "bbox",
            [2, 5, 3, 15, 40, 20],
            "resample",
            [1.5, 1, 2],
        ),
        rotate_file_case(
            "rotate_file.two",
            frame_tree("uint16", ("a_ch0", "b_ch0")),
            ["a_ch0.tif", "b_ch0.tif"],
            PX,
            0.3,
            *base,
        ),
        rotate_file_case(
            "rotate_file.exists",
            frame_tree("uint16", ("a_ch0", "b_ch0"), ("a_ch0",)),
            ["a_ch0.tif", "b_ch0.tif"],
            PX,
            0.3,
            *base,
        ),
        rotate_file_case(
            "rotate_file.overwrite",
            frame_tree("uint16", ("a_ch0",), ("a_ch0",)),
            ["a_ch0.tif"],
            PX,
            0.3,
            *base,
            "Overwrite",
            True,
        ),
        rotate_file_case(
            "rotate_file.missing",
            frame_tree("uint16"),
            ["none.tif", "scan_ch0.tif"],
            PX,
            0.3,
            *base,
        ),
    ]


def p_interp_matlab(vol, xstep, stepsize, _name, reverse):
    return dv.skewed_space_interp_matlab(vol, xstep, stepsize, reverse=reverse)


def p_colon(a, d, b):
    return dv._colon(a, d, b).reshape(1, -1)


def _interp_function_cases() -> list[Case]:
    cases = []
    # x steps of dz 0.3, 0.6 (step 1/3 in the deskew), 0.8 and 1.5 at 32.45 degrees
    for dz, step in ((0.3, 0.5), (0.6, 1 / 3), (0.8, 0.25), (1.5, 1 / 6), (0.39, 0.5)):
        xstep = math.cos(SKEW * math.pi / 180) * dz / PX
        for reverse in (True, False):
            for cls in ("single", "uint16", "double"):
                cases.append(
                    Case(
                        f"mfunction.{_id(dz)}.{'rev' if reverse else 'fwd'}.{cls}",
                        matlab="skewed_space_interp_defined_stepsize",
                        args=[Vol(In("v")), xstep, step, "Reverse", reverse],
                        python=p_interp_matlab,
                        out_roles=("vol",),
                        inputs={"v": _vol((9, 23, 13), cls)},
                    )
                )
    cases.append(
        Case(
            "mfunction.step1",
            matlab="skewed_space_interp_defined_stepsize",
            args=[Vol(In("v")), 2.5, 1, "Reverse", True],
            python=p_interp_matlab,
            out_roles=("vol",),
            inputs={"v": _vol((9, 23, 13), "uint16")},
        )
    )
    for a, d, b in (
        (1, 1 / 3, 30),
        (1, 0.5, 13),
        (1, 0.1, 7),
        (1, 1 / 7, 100),
        (1, 0.3, 1),
        (1, 1 / 6, 12),
        # an even count whose middle value (the mean of the ends) is not a + k d
        (1, 1 / 49, 7),
        (1, 1 / 49, 13),
    ):
        cases.append(
            Case(
                f"colon.{_id(d)}.{b}",
                code="out = {in{1} : in{2} : in{3}};",
                args=[float(a), d, float(b)],
                python=p_colon,
            )
        )
    return cases


# ------------------------------------------------------------------------------------------
# the golden scenarios of GoldenOutputTest (tests/matlab_map/test_golden.py's cases)
# ------------------------------------------------------------------------------------------

GOLDEN = {
    "deskew_rotate_sample_scan": (0.2, None, True),
    "deskew_rotate_skewed_interp": (0.3, None, True),
    "deskew_rotate_resampled": (0.2, (2, 2, 2), False),
}
GOLDEN_HELPERS = ["tests/matlab_map/matlab"]


# the inputs cases of test_golden.py (the same keys: their stored results)
GOLDEN_INPUTS = {
    s: Case(
        f"golden.{s}.inputs",
        matlab="pkpy_golden",
        args=["inputs", s, mcell("vol"), Work()],
        outputs=["vol"],
        helpers=GOLDEN_HELPERS,
        files=False,
    )
    for s in GOLDEN
}


def p_golden(scenario: str):
    """GoldenOutputTest.computeScenario of a deskew/rotate scenario with the port"""
    ref = oracle.reference(GOLDEN_INPUTS[scenario])
    if any(getattr(v, "data", None) is None for v in ref.outputs):
        pytest.skip("the scenario's inputs are in the golden store (the manifest has hashes)")
    vol = from_matlab(ref.outputs[0].data)
    dz, rs, rounded = GOLDEN[scenario]
    out = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=True, resample_factor=rs)
    if rounded:
        from petakit5d.core import dtypes

        out = dtypes.uint16(dtypes.mround(out))
    return [out]


def _golden_cases() -> list[Case]:
    cases = list(GOLDEN_INPUTS.values())
    for s in GOLDEN:
        for what in ("compute", "stored"):
            # the stored references in the MATLAB repo too: bit for bit (GoldenOutputTest allows a
            # relative L2 error of 1e-4 between platforms)
            cases.append(
                Case(
                    f"golden.{s}.{what}",
                    matlab="pkpy_golden",
                    args=[what, s, mcell("dsr"), Work()],
                    outputs=["dsr"],
                    python=p_golden,
                    py_args=[s],
                    py_out=list,
                    out_roles=("vol",),
                    helpers=GOLDEN_HELPERS,
                    files=False,
                )
            )
    return cases


# ------------------------------------------------------------------------------------------
# gpu_process: MATLAB's gpuProcess code path computed on the CPU, bit for bit
# ------------------------------------------------------------------------------------------

# deskewRotateFrame3D.m with gpuProcess (lines 50-213; sample scan, no resampling, linear): the
# interpolation mex where the x step is above the threshold, then imwarp with the combined
# matrix. MATLAB runs that imwarp on a gpuArray; here it runs on the CPU.
GPU_PATH_CODE = """
vol = in{1}; angle = in{2}; dz = in{3}; xyPixelSize = in{4}; Reverse = in{5}; save16bit = in{6};
xStepThresh = 2.0;
[ny, nx, nz] = size(vol);
theta = angle * pi/180;
dx = cos(theta)*dz/xyPixelSize;
zAniso = sin(abs(theta)) * dz / xyPixelSize;
outSize = round([ny, (nx-1)*cos(theta)+(nz-1)*zAniso/sin(abs(theta)), ...
    (nx-1)*sin(abs(theta))-4]);
vol_1 = vol;
if abs(dx) > xStepThresh
    if abs(dx) / xStepThresh < 1.5
        dzout_thresh = xyPixelSize * xStepThresh / cos(theta);
        dzout = dz / (ceil((dz / dzout_thresh) / (1 / 1)) *(1 / 1));
        counter = 1;
        while dzout / dzout_thresh < 0.4
            dzout = dz / (ceil((dz / dzout_thresh) / (1 / (counter + 1))) *(1 / (counter + 1)));
            counter = counter + 1;
        end
    else
        dzout = dz / ceil(abs(dx) / xStepThresh);
    end
    int_stepsize = dzout / dz;
    vol_1 = skewed_space_interp_defined_stepsize_mex(vol, abs(dx), int_stepsize, ...
        Reverse, save16bit);
    [ny, nx, nz] = size(vol_1);
    dz = dzout;
    dx = cos(theta)*dz/xyPixelSize;
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end
if ~Reverse
    xshift = -dx; xstep = dx;
else
    xshift = dx + ceil((nz-1)*dx); xstep = -dx;
end
nxDs = ceil((nz-1)*dx) + nx;
ds_S = [1 0 0 0; 0 1 0 0; xstep 0 1 0; xshift 0 0 1];
if Reverse
    theta = -theta;
end
center = ([ny nxDs nz]+1)/2;
T1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; -center([2 1 3]) 1];
S1 = [1 0 0 0; 0 1 0 0; 0 0 zAniso 0; 0 0 0 1];
R = [cos(theta) 0 -sin(theta) 0; 0 1 0 0; sin(theta) 0 cos(theta) 0; 0 0 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; (outSize([2 1 3])+1)/2 1];
M = ds_S*(T1*S1*R*T2)*(eye(4)*eye(4)*eye(4));
RA = imref3d(outSize, 1, 1, 1);
out = {imwarp(vol_1, affine3d(M), 'linear', 'FillValues', 0, 'OutputView', RA)};
"""


def _gpu_cases() -> list[Case]:
    cases = []
    for sname, size, seed in (("small", SIZE, 0), ("odd", (67, 101, 53), 1), ("beads", None, 0)):
        for kname, (cls, s16) in KINDS.items():
            if sname == "beads" and cls == "uint16":
                continue
            for dz in (0.2, 0.3):
                for reverse in (True, False):
                    if sname == "beads":
                        vin = _vol((24, 64, 40), "beads")
                    else:
                        vin = _vol(size, cls, seed)
                    cases.append(
                        Case(
                            f"gpu.{sname}.{kname}.{_id(dz)}.{'rev' if reverse else 'fwd'}",
                            code=GPU_PATH_CODE,
                            args=[In("vol"), SKEW, dz, PX, reverse, s16],
                            python=p_dsr,
                            py_args=[Vol(In("vol")), SKEW, dz, PX],
                            py_kwargs={"reverse": reverse, "save_16bit": s16, "gpu_process": True},
                            out_roles=("vol",),
                            inputs={"vol": vin},
                        )
                    )
    return cases


CASES = (
    _dsr_cases()
    + _dsr_error_cases()
    + _geometry_cases()
    + _deskew_rotate_frame_cases()
    + _rotate_file_cases()
    + _interp_function_cases()
    + _golden_cases()
    + _gpu_cases()
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_case_ids_and_counts():
    ids = [c.id for c in CASES]
    assert len(ids) == len(set(ids))
    assert geometry_table(1000, 0, False).shape == (1000, 13)
    assert deskew_size_table(500).shape == (500, 8)
    assert rotate_size_table(500).shape == (500, 16)


# ------------------------------------------------------------------------------------------
# DeskewRotateFrame3DTest in intent: a grid of beads with known raw positions
# ------------------------------------------------------------------------------------------

RAW_SIZE = (64, 96, 80)  # MATLAB [ny nx nz]


def _sind(x: float) -> float:
    from petakit5d.matlabcompat import elementary

    return float(np.asarray(elementary.sind(np.array([x]))).reshape(-1)[0])


def _cosd(x: float) -> float:
    from petakit5d.matlabcompat import elementary

    return float(np.asarray(elementary.cosd(np.array([x]))).reshape(-1)[0])


def bead_grid():
    """15 beads on a raw x-z grid, each at its own y (raw: MATLAB's 1-based [y x z]); the
    single volume (z, y, x)"""
    xs, zs = (30, 48, 66), (12, 26, 40, 54, 68)
    raw = np.array(
        [[8 + 3.4 * ((k * 7) % 15), xs[k % 3], zs[k // 3]] for k in range(15)], dtype=np.float64
    )
    pos = np.stack([raw[:, 2] - 1, raw[:, 0] - 1, raw[:, 1] - 1], axis=1)
    return raw, synth.beads(size_from_matlab(RAW_SIZE), pos, 2.0, 1000.0)


def sample_scan_out_size(dz: float) -> tuple[int, int, int]:
    """MATLAB's [ny nx nz]"""
    sz = RAW_SIZE
    theta = SKEW * math.pi / 180
    v = (
        sz[0],
        (sz[1] - 1) * math.cos(theta) + (sz[2] - 1) * dz / PX,
        (sz[1] - 1) * math.sin(theta) - 4,
    )
    return tuple(int(math.floor(x + 0.5)) for x in v)


def expected_sample_scan_map(dz: float, reverse: bool, raw: np.ndarray):
    """[y; x; z]_out = map [y; x; z; 1]_raw of the centered deskew/rotation convention"""
    sz = RAW_SIZE
    theta = SKEW * math.pi / 180
    q = dz / PX
    dx = math.cos(theta) * q
    f = math.ceil((sz[2] - 1) * dx) - (sz[2] - 1) * dx
    c = (np.array(sample_scan_out_size(dz), dtype=np.float64) + 1) / 2
    s = 1 - 2 * reverse
    m = np.array([[1, 0, 0, 0], [0, math.cos(theta), s * q, 0], [0, -s * math.sin(theta), 0, 0]])
    m[:, 3] = np.array([0, c[1] - s * math.cos(theta) * f / 2, c[2] + math.sin(theta) * f / 2]) - m[
        :, :3
    ] @ np.array([0, (sz[1] + 1) / 2, (sz[2] + 1) / 2])
    return m, np.c_[raw, np.ones(len(raw))] @ m.T


def fit_map(out: np.ndarray, raw: np.ndarray, approx: np.ndarray):
    """least-squares affine map from bead centroids near their approximate output positions
    (both in MATLAB's 1-based [y x z])"""
    c = np.array(
        [synth.centroid(out, (p[2] - 1, p[0] - 1, p[1] - 1), 6)[[1, 2, 0]] + 1 for p in approx]
    )
    a = np.c_[raw, np.ones(len(raw))]
    sol, *_ = np.linalg.lstsq(a, c, rcond=None)
    return sol.T, float(np.abs(a @ sol - c).max())


def total(x) -> float:
    return math.fsum(np.asarray(x, dtype=np.float64).ravel().tolist())


def rel_l2(a, e) -> float:
    a, e = np.asarray(a, dtype=np.float64), np.asarray(e, dtype=np.float64)
    return float(np.linalg.norm((a - e).ravel()) / max(np.linalg.norm(e.ravel()), 2.0**-52))


def _reject(*_a, **_k):
    err = ValueError("stubbed kernel")
    err.identifier = "volume_warp:inputError"
    raise err


DZ_REV = [(dz, r) for dz in (0.2, 0.3) for r in (True, False)]


@pytest.mark.parametrize("dz,reverse", DZ_REV)
def test_output_size_matches_geometry(dz, reverse):
    vol = np.zeros(size_from_matlab(RAW_SIZE), dtype=np.float32)
    out = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=reverse)
    assert to_matlab(out).shape == sample_scan_out_size(dz)


@pytest.mark.parametrize("dz,reverse", DZ_REV)
def test_geometry_matches_acquisition(dz, reverse):
    raw, vol = bead_grid()
    out = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=reverse)
    expected, pos = expected_sample_scan_map(dz, reverse, raw)
    actual, residual = fit_map(out, raw, pos)
    np.testing.assert_allclose(actual[:, :3], expected[:, :3], atol=0.005 * dz / PX, rtol=0)
    assert residual < 0.6


@pytest.mark.parametrize("dz,reverse", DZ_REV)
def test_canvas_origin_matches_convention(dz, reverse):
    raw, vol = bead_grid()
    out = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=reverse)
    _, pos = expected_sample_scan_map(dz, reverse, raw)
    actual, _ = fit_map(out, raw, pos)
    np.testing.assert_allclose(np.c_[raw, np.ones(len(raw))] @ actual.T, pos, atol=0.25, rtol=0)


@pytest.mark.parametrize("dz,reverse", DZ_REV)
def test_total_intensity_scales_with_voxel_volume(dz, reverse):
    _, vol = bead_grid()
    out = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=reverse)
    z_aniso = _sind(SKEW) * dz / PX
    assert total(out) / total(vol) == pytest.approx(z_aniso, rel=0.02)


def test_uint16_input_with_save16bit_gives_uint16():
    from petakit5d.core import dtypes

    _, vol = bead_grid()
    vol = dtypes.uint16(dtypes.plus(vol, 100))
    out = dv.deskew_rotate_volume(vol, SKEW, 0.2, PX, reverse=True, save_16bit=True)
    assert out.dtype == np.uint16
    out_single = dv.deskew_rotate_volume(vol.astype(np.float32), SKEW, 0.2, PX, reverse=True)
    assert out_single.dtype == np.float32
    d = np.abs(out.astype(np.float64) - dtypes.mround(out_single).astype(np.float64))
    assert d.max() <= 1


def test_resample_factor_scales_output_size():
    vol = np.zeros(size_from_matlab(RAW_SIZE), dtype=np.float32)
    out = dv.deskew_rotate_volume(vol, SKEW, 0.2, PX, reverse=True, resample_factor=(2, 2, 2))
    expected = tuple(int(math.floor(v / 2 + 0.5)) for v in sample_scan_out_size(0.2))
    assert to_matlab(out).shape == expected


def test_bbox_crops_the_full_output():
    from petakit5d.core.bbox import box_from_matlab, box_slices

    _, vol = bead_grid()
    full = dv.deskew_rotate_volume(vol, SKEW, 0.2, PX, reverse=True)
    # a y extent that isn't a multiple of 8 runs the scalar tail of the SIMD loops too
    box = box_from_matlab([9, 21, 5, 53, 150, 40])
    out = dv.deskew_rotate_volume(vol, SKEW, 0.2, PX, reverse=True, bbox=box)
    assert out.shape == full[box_slices(box)].shape
    assert rel_l2(out, full[box_slices(box)]) <= 1e-6


def test_separate_kernel_steps_replace_the_combined_kernel(monkeypatch):
    """where the fused kernel rejects its input, its two steps run (dz 0.3)"""
    from petakit5d.core.bbox import box_from_matlab
    from petakit5d.kernels import deskew

    _, vol = bead_grid()
    box = box_from_matlab([9, 21, 5, 53, 150, 40])
    combined = dv.deskew_rotate_volume(vol, SKEW, 0.3, PX, reverse=True, bbox=box)
    monkeypatch.setattr(deskew, "skewed_space_interp_volume_deskew_rotate_warp", _reject)
    separate = dv.deskew_rotate_volume(vol, SKEW, 0.3, PX, reverse=True, bbox=box)
    assert rel_l2(separate, combined) <= 1e-6


@pytest.mark.parametrize("dz", [0.2, 0.3])
def test_bbox_past_the_input_rows_matches_imwarp(dz, monkeypatch):
    from petakit5d.core.bbox import box_from_matlab
    from petakit5d.kernels import deskew

    _, vol = bead_grid()
    box = box_from_matlab([9, 21, 5, RAW_SIZE[0] + 12, 150, 40])
    out = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=True, bbox=box)
    monkeypatch.setattr(deskew, "volume_deskew_rotate_warp", _reject)
    monkeypatch.setattr(deskew, "skewed_space_interp_volume_deskew_rotate_warp", _reject)
    ref = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=True, bbox=box)
    n_in = RAW_SIZE[0] - 9 + 1
    assert out.shape == tuple(box[i + 3] - box[i] for i in range(3))
    assert not np.any(ref[:, n_in:, :])
    assert not np.any(out[:, n_in:, :]), "values past the input's rows"
    assert rel_l2(out, ref) <= 1e-6


def test_forward_interpolation_ends_at_the_last_raw_slice():
    """dz 0.6 (step 1/3): the last grid slice passes the last raw slice by a rounding error"""
    from petakit5d.kernels import deskew

    xstep = _cosd(SKEW) * 0.6 / PX
    v = (CounterRNG(4, "dsr/fwdend").uniform((30, 16, 40)) * 1000).astype(np.float32)
    for save_16bit in (True, False):
        for inp in (v, v.astype(np.uint16)):
            out = deskew.skewed_space_interp_defined_stepsize(inp, xstep, 1 / 3, False, save_16bit)
            assert out.shape[0] == 88
            assert np.array_equal(out[-1], inp[-1].astype(out.dtype))


def test_inputs_narrower_than_the_shift_run_the_matlab_code(monkeypatch):
    from petakit5d.kernels import deskew

    v = (CounterRNG(5, "dsr/narrow").uniform((12, 8, 6))).astype(np.float32)
    for call in (
        lambda: deskew.skewed_space_interp_volume_deskew_rotate_warp(
            v, 4.69, 1 / 3, True, np.eye(4), (0, 0, 0, 5, 8, 5), False
        ),
        lambda: deskew.skewed_space_interp_defined_stepsize(v, 4.69, 1 / 3, True, False),
    ):
        with pytest.raises(ValueError) as err:
            call()
        assert err.value.identifier == "volume_warp:inputError"
    # dz 1.5 (x step 11.7), 13 columns
    v = (CounterRNG(6, "dsr/narrow13").uniform((12, 8, 13)) * 1000).astype(np.float32)
    out = dv.deskew_rotate_volume(v, SKEW, 1.5, PX, reverse=True)
    monkeypatch.setattr(deskew, "skewed_space_interp_volume_deskew_rotate_warp", _reject)
    monkeypatch.setattr(deskew, "skewed_space_interp_defined_stepsize", _reject)
    ref = dv.deskew_rotate_volume(v, SKEW, 1.5, PX, reverse=True)
    assert np.count_nonzero(ref) > 0
    assert np.array_equal(out, ref)


def test_single_slice_or_column_reads_nothing_outside():
    from petakit5d.kernels import deskew

    v = (CounterRNG(7, "dsr/slice").uniform((6, 8, 40)) * 1000 + 1).astype(np.float32)
    box = (0, 0, 0, 5, 8, 5)
    eye = np.eye(4)
    zeros = np.zeros((5, 8, 5), dtype=np.float32)
    assert np.array_equal(deskew.volume_deskew_rotate_warp(v[:1], eye, box, False), zeros)
    col = np.ascontiguousarray(v[:, :, :1])
    assert np.array_equal(deskew.volume_deskew_rotate_warp(col, eye, box, False), zeros)
    u = v[:1].astype(np.uint16)
    assert np.array_equal(
        deskew.volume_deskew_rotate_warp(u, eye, box, True), zeros.astype(np.uint16)
    )
    fused = deskew.skewed_space_interp_volume_deskew_rotate_warp(
        v[:1], 4.69, 1 / 3, True, eye, box, False
    )
    assert np.array_equal(fused, zeros)
    # with two slices the points on them have values
    assert np.count_nonzero(deskew.volume_deskew_rotate_warp(v[:2], eye, box, False)) > 0


def test_combined_matches_separate_deskew_then_rotate():
    """the combined transform equals the deskew followed by the rotation (up to the extra
    interpolation of the separate path)"""
    dz = 0.2
    raw, vol = bead_grid()
    combined = dv.deskew_rotate_volume(vol, SKEW, dz, PX, reverse=True)
    theta = SKEW * math.pi / 180
    ds = dv.deskew_volume(vol, SKEW, dz, PX, True, crop=False)
    separate = dv.rotate_volume(ds, SKEW, geo._sin(theta) * dz / PX, True, out_size=combined.shape)
    _, pos = expected_sample_scan_map(dz, True, raw)
    np.testing.assert_allclose(
        fit_map(separate, raw, pos)[0], fit_map(combined, raw, pos)[0], atol=0.1, rtol=0
    )
    # the separate path interpolates twice, which blurs the beads
    assert rel_l2(separate, combined) <= 0.1


@pytest.mark.parametrize("reverse", [True, False])
def test_objective_scan_geometry_with_cubic_interpolation(reverse):
    """objective scan is a pure rotation (no deskew)"""
    dz = 0.2
    q = dz / PX
    raw, vol = bead_grid()
    out = dv.deskew_rotate_volume(
        vol, SKEW, dz, PX, objective_scan=True, reverse=reverse, interp_method="cubic"
    )
    sz = RAW_SIZE
    theta = SKEW * math.pi / 180
    o = tuple(
        int(math.floor(v + 0.5))
        for v in (
            sz[0],
            sz[1] * math.cos(theta) + sz[2] * q * math.sin(theta),
            sz[2] * q * math.cos(theta) + sz[1] * math.sin(theta),
        )
    )
    assert to_matlab(out).shape == o
    th = -theta if reverse else theta
    expected = np.array(
        [
            [1, 0, 0, 0],
            [0, math.cos(th), math.sin(th) * q, 0],
            [0, -math.sin(th), math.cos(th) * q, 0],
        ]
    )
    expected[:, 3] = np.array([0, (o[1] + 1) / 2, (o[2] + 1) / 2]) - expected[:, :3] @ np.array(
        [0, (sz[1] + 1) / 2, (sz[2] + 1) / 2]
    )
    pos = np.c_[raw, np.ones(len(raw))] @ expected.T
    actual, _ = fit_map(out, raw, pos)
    np.testing.assert_allclose(actual[:, :3], expected[:, :3], atol=0.005 * q, rtol=0)
    np.testing.assert_allclose(np.c_[raw, np.ones(len(raw))] @ actual.T, pos, atol=0.25, rtol=0)
    assert total(out) / total(vol) == pytest.approx(q, rel=0.02)


def test_objective_scan_preserves_content():
    dz = 0.2
    _, vol = bead_grid()
    out = dv.deskew_rotate_volume(vol, SKEW, dz, PX, objective_scan=True)
    assert total(out) / total(vol) == pytest.approx(dz / PX, rel=0.02)


def test_save16bit_rounds_to_nearest():
    from petakit5d.core import dtypes

    _, vol = bead_grid()
    vol = dtypes.uint16(dtypes.plus(vol, 100))
    out = dv.deskew_rotate_volume(vol, SKEW, 0.2, PX, reverse=True, save_16bit=True)
    ref = dv.deskew_rotate_volume(vol.astype(np.float32), SKEW, 0.2, PX, reverse=True)
    # values within 1e-3 of a tie can round either way
    r = ref.astype(np.float64)
    not_tie = np.abs(r - np.floor(r) - 0.5) > 1e-3
    mismatched = np.count_nonzero(out[not_tie] != dtypes.uint16(ref[not_tie])) / not_tie.sum()
    assert mismatched < 1e-3


# ------------------------------------------------------------------------------------------
# the port's own: layouts, threads, the MATLAB-name layer, MATLAB's product
# ------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "dz,opts",
    [
        (0.2, {}),
        (0.3, {"reverse": True}),
        (0.3, {"resample_factor": (2, 1, 1.5)}),
        (0.2, {"objective_scan": True}),
        (1.5, {"reverse": True}),
    ],
)
def test_layouts_and_thread_counts_give_the_same_result(dz, opts):
    v = from_matlab(make_volume((19, 41, 23), "uint16"))
    ref = dv.deskew_rotate_volume(v, SKEW, dz, PX, save_16bit=True, **opts)
    vm = np.ascontiguousarray(v.transpose(0, 2, 1))
    for n in (1, 2, 5):
        out = dv.deskew_rotate_volume_matlab_layout(
            vm, SKEW, dz, PX, save_16bit=True, num_threads=n, **opts
        )
        assert np.array_equal(out.transpose(0, 2, 1), ref)


def test_gpu_process_runs_matlabs_gpu_code_path_on_the_cpu():
    """gpu_process runs imwarp after the interpolation, as MATLAB's gpuProcess, with its output
    class (uint16 input without save16bit and interpolation stays uint16), on both layouts"""
    from petakit5d.matlabcompat.warp import imwarp

    v = from_matlab(make_volume((19, 41, 23), "uint16"))
    out = dv.deskew_rotate_volume(v, SKEW, 0.2, PX, reverse=True, gpu_process=True)
    g = geo.deskew_rotate_geometry(v.shape, SKEW, 0.2, PX, reverse=True, gpu_process=True)
    assert not g.use_fast_method and out.dtype == np.uint16
    assert np.array_equal(out, imwarp(v, g.affine, "linear", output_size=g.out_size))
    for dz in (0.2, 0.3):
        ref = dv.deskew_rotate_volume(v, SKEW, dz, PX, gpu_process=True)
        vm = np.ascontiguousarray(v.transpose(0, 2, 1))
        ml = dv.deskew_rotate_volume_matlab_layout(vm, SKEW, dz, PX, gpu_process=True)
        assert np.array_equal(ml.transpose(0, 2, 1), ref)


def test_the_matlab_name_layer_converts_boxes_factors_and_volumes():
    from petakit5d import matlab
    from petakit5d.core.bbox import box_from_matlab

    m = make_volume((19, 41, 23), "single")
    mbox = [3, 9, 2, 15, 50, 9]
    out = matlab.deskewRotateFrame3D(m, SKEW, 0.3, PX, "Reverse", True, "bbox", mbox)
    ref = dv.deskew_rotate_volume(
        from_matlab(m), SKEW, 0.3, PX, reverse=True, bbox=box_from_matlab(mbox)
    )
    assert np.array_equal(out, to_matlab(ref))
    out = matlab.rotateFrame3D(m, SKEW, 1.5, True, "resample", [1.2, 1.7, 0.8])
    ref = dv.rotate_volume(from_matlab(m), SKEW, 1.5, True, resample=(0.8, 1.7, 1.2))
    assert np.array_equal(out, to_matlab(ref))


def test_interp_method_is_not_a_parameter_of_the_separate_functions():
    """MATLAB's deskewFrame3D and rotateFrame3D reject 'interpMethod' (their parameter is
    'Interp'; measured: MATLAB:InputParser:UnmatchedParameter), as RLdecon.m passes it"""
    from petakit5d import matlab

    m = make_volume((4, 20, 10), "single")
    for call in (
        lambda: matlab.deskewFrame3D(m, SKEW, 0.3, PX, False, "Crop", False, "interpMethod", "x"),
        lambda: matlab.rotateFrame3D(m, SKEW, 1.5, False, "interpMethod", "linear"),
    ):
        with pytest.raises(Exception) as err:
            call()
        assert err.value.identifier == "MATLAB:InputParser:UnmatchedParameter"


def test_the_matlab_interpolation_sets_the_slices_of_the_transform(monkeypatch):
    """after MATLAB's interpolation (where the kernel rejects the input) the transform is made
    for the slices it returned (numel(1 : step : nz), which can differ from the kernels' grid)"""
    from petakit5d.matlabcompat.warp import imwarp

    v = from_matlab(make_volume((9, 30, 14), "double"))
    real = dv.skewed_space_interp_matlab

    def fewer(*a, **k):
        return real(*a, **k)[:-1]

    monkeypatch.setattr(dv, "skewed_space_interp_matlab", fewer)
    out = dv.deskew_rotate_volume(v, SKEW, 0.3, PX, reverse=True, resample_factor=(2, 2, 2))
    vol_1 = fewer(v, abs(geo._cos(SKEW * math.pi / 180) * 0.3 / PX), 0.5, reverse=True)
    g = geo.deskew_rotate_geometry(
        v.shape, SKEW, 0.3, PX, reverse=True, resample_factor=(2, 2, 2), interp_nz=vol_1.shape[0]
    )
    assert g.nz == vol_1.shape[0]
    assert np.array_equal(out, imwarp(vol_1, g.affine, "linear", output_size=g.out_size))


def test_registered_names():
    from petakit5d.params import find

    for mname, pname in (
        ("deskewRotateFrame3D", "deskew_rotate_volume"),
        ("deskewFrame3D", "deskew_volume"),
        ("rotateFrame3D", "rotate_volume"),
        ("XR_RotateFrame3D", "rotate_file"),
    ):
        name = find(mname).python_name
        assert name.startswith("petakit5d.deskew_rotate.") and name.endswith(f".{pname}")
        assert getattr(__import__("petakit5d.deskew_rotate").deskew_rotate, pname).__name__ == pname


def test_matlab_matmul_sums_from_positive_zero():
    """MATLAB's product of 4 x 4 doubles starts each sum at +0: products that are all -0 give
    +0 (measured)"""
    a = [[0.0, 0.0, 0.0, 0.0]] * 4
    b = [[-1.0, -2.0, -3.0, -4.0]] * 4
    out = geo.matlab_matmul(a, b)
    assert all(math.copysign(1.0, v) == 1.0 for row in out for v in row)


def test_interp_method_and_resample_errors():
    v = np.zeros((5, 6, 20), dtype=np.float32)
    with pytest.raises(ValueError) as err:
        dv.deskew_rotate_volume(v, SKEW, 0.3, PX, interp_method="nearest")
    assert err.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
    with pytest.raises(IndexError) as err:
        dv.rotate_volume(v, SKEW, 1.5, resample=(2, 2))
    assert err.value.identifier == "MATLAB:badsubscript"


def test_colon_counts_and_ends():
    assert dv._colon(1.0, 1 / 3, 30.0).size == 88
    z = dv._colon(1.0, 0.1, 7.0)
    assert z.size == 61 and z[0] == 1.0 and z[-1] == 7.0
    assert dv._colon(1.0, 0.3, 1.0).tolist() == [1.0]
    assert dv._colon(2.0, 0.5, 1.0).size == 0
