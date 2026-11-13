"""Deskew/rotate of files (petakit5d.deskew_rotate.file, .memory) against MATLAB.

Oracle cases (MATLAB runs the oracle tree's functions on a copy of an input folder, Python the
port through the MATLAB-name layer; the folders they leave are compared, files bit for bit):

- ``XR_deskewRotateFrame`` (``deskew_rotate_file``): the option matrix of
  XRDeskewRotateFrameTest and more: the combined path (deskew kernels at dz 0.2, 0.3 and 0.5,
  forward and reverse; ``imwarp`` for objective scans, resampling and cubic interpolation),
  the separate path (deskew only, with rotation, objective and z stage scans, above 1000
  slices in blocks, a DS made earlier and rotated in a second call, overwrite), tiff and zarr
  input and output, uint16, single and uint8 data, ``save16bit`` on and off, flat field (with
  and without a constant offset), input boxes, flips, input and output axis orders, partial
  files, an existing result, ``save3DStack`` and ``saveMIP`` off, other folder names, the
  rescale (and its MATLAB error for 16-bit results), the background removal's known issue and
  a missing input;
- ``scaleContrast`` (another author's GPL function, used as a black box; ``scale_contrast``):
  classes, ranges, NaN and Inf, degenerate ranges, its errors and defaults;
- ``XR_estimateDeskewRotateMemory``: the output size of every geometry (the port's memory is
  its own, measured).

Tests without MATLAB: XRDeskewRotateFrameTest and DeskewRotateMemoryTest in intent (files equal
the in-memory functions' results, the publish order, private zarr writes, the known issues),
and the port's own rules (the flat-field cache, partial files read into one array, the read
buffer, thread counts, the memory model).

Inputs are made in the port's form and handed to MATLAB in its form (tiff files are the same
bytes; zarr inputs are written in MATLAB's layout).
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.deskew_rotate import file as df
from petakit5d.deskew_rotate import memory as dm
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, Vol
from petakit5d.testing.compare import ignore, tree
from petakit5d.testing.synth import CounterRNG

SKEW = 32.45
PX = 0.108
FSN = "Scan_Iter_0000_0000_CamA_ch0_stack0000_488nm_0000000msec_0106060251msecAbs"
# the raw volume of the cases (z, y, x): not cubic, x wide enough for the kernels
RAW = (29, 37, 61)

# ------------------------------------------------------------------------------------------
# inputs (the port's form)
# ------------------------------------------------------------------------------------------


def make_frame(shape=RAW, kind: str = "uint16", seed: int = 0) -> np.ndarray:
    """beads on a background of 100 with noise, ``(z, y, x)``: uint16, single or uint8"""
    nz, ny, nx = shape
    rng = CounterRNG(seed, f"dsrfile/{kind}/{'x'.join(map(str, shape))}")
    pos = []
    for k in range(9):
        pos.append(
            (2 + (k * 7) % max(1, nz - 4), 3 + (k * 5) % max(1, ny - 6), 4 + (k * 11) % (nx - 8))
        )
    vol = synth.beads(shape, pos, 1.5, 1000.0, 100.0)
    noise = rng.integers(0, 50, shape, dtype=np.uint16)
    if kind == "single":
        return np.ascontiguousarray(vol + noise.astype(np.float32))
    v = np.floor(vol).astype(np.uint16) + noise
    if kind == "uint8":
        return np.ascontiguousarray((v // 8).astype(np.uint8))
    return np.ascontiguousarray(v)


def frame_tree(
    kind: str = "uint16",
    shape=RAW,
    *,
    parts: int = 1,
    zarr: bool = False,
    ff: bool = False,
    existing: tuple[str, ...] = (),
    sub: str = "",
):
    """a folder with the frame ``<FSN>.tif`` (``parts`` partial files ``<FSN>_part0001.tif``,
    ...; or ``<FSN>.zarr`` in MATLAB's layout), flat-field and background images ``ff.tif``,
    ``bg.tif`` and existing results (folder-relative file names, uint16 7s)"""

    def make(folder):
        d = Path(folder) / sub
        d.mkdir(parents=True, exist_ok=True)
        vol = make_frame(shape, kind)
        if zarr:
            synth.write_zarr(d / f"{FSN}.zarr", vol, (16, 32, 32), layout="matlab")
        else:
            edges = np.linspace(0, shape[0], parts + 1).astype(int)
            for k in range(parts):
                name = FSN if k == 0 else f"{FSN}_part{k:04d}"
                synth.write_tiff(d / f"{name}.tif", vol[edges[k] : edges[k + 1]])
        if ff:
            ny, nx = shape[1], shape[2]
            bg = np.full((ny + 3, nx + 2), 100, dtype=np.float32)
            ramp = np.linspace(500, 1000, nx + 2, dtype=np.float64).astype(np.float32)
            synth.write_tiff(d / "ff.tif", np.ascontiguousarray(bg + ramp[np.newaxis]))
            synth.write_tiff(d / "bg.tif", bg)
        for e in existing:
            p = Path(folder) / e
            p.parent.mkdir(parents=True, exist_ok=True)
            synth.write_tiff(p, np.full((3, 4, 5), 7, dtype=np.uint16))

    return make


# ------------------------------------------------------------------------------------------
# XR_deskewRotateFrame cases
# ------------------------------------------------------------------------------------------

# runs{r} is one call's arguments after the file; '@d/x' is the input folder's file x
_CODE = (
    "f = cellfun(@(x) [in{1}, x], strsplit(in{2}, '|'), 'unif', 0); "
    "if numel(f) == 1, f = f{1}; end; runs = in{3}; "
    "for r = 1 : numel(runs), a = runs{r}; "
    "for k = 1 : numel(a), if ischar(a{k}) && startsWith(a{k}, '@d/'), "
    "a{k} = [in{1}, a{k}(4 : end)]; end, end; "
    "XR_deskewRotateFrame(f, a{:}); end; out = {};"
)


def p_frame(d, frames, runs):
    from petakit5d import matlab

    f = [d + x for x in frames.split("|")]
    for a in runs:
        a = [d + v[3:] if isinstance(v, str) and v.startswith("@d/") else v for v in a]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            matlab.XR_deskewRotateFrame(f[0] if len(f) == 1 else f, *a)
    return []


def frame_case(cid, make_tree, runs, frames=f"{FSN}.tif", **kw):
    return Case(
        f"frame.{cid}",
        code=_CODE,
        nargout=0,
        args=[InPath("d"), frames, [list(r) for r in runs]],
        python=p_frame,
        inputs={"d": Tree(make_tree, "d")},
        **kw,
    )


def _r(*extra, dz=0.3, reverse=True):
    """one call: pixel size, dz, skew angle, reverse and the extra name-value pairs"""
    return [PX, dz, "skewAngle", SKEW, "reverse", reverse, *extra]


def _frame_cases() -> list[Case]:
    t = frame_tree()
    sep = ("DSRCombined", False, "rotate", True)
    ds_only = ("DSRCombined", False, "rotate", False)
    ffa = ("FFCorrection", True, "FFImage", "@d/ff.tif", "backgroundImage", "@d/bg.tif")
    cases = [
        # combined
        frame_case("combined", t, [_r()]),
        frame_case("combined.forward", t, [_r(reverse=False)]),
        frame_case("combined.dz02", t, [_r(dz=0.2)]),
        frame_case("combined.dz05", t, [_r(dz=0.5)]),
        frame_case("combined.zarr", t, [_r("saveZarr", True)]),
        frame_case("combined.zarr_block", t, [_r("saveZarr", True, "blockSize", [16, 32, 8])]),
        frame_case("combined.single_out", t, [_r("save16bit", False)]),
        frame_case("combined.single_in", frame_tree("single"), [_r()]),
        frame_case("combined.single_in_out", frame_tree("single"), [_r("save16bit", False)]),
        frame_case("combined.uint8_in", frame_tree("uint8"), [_r()]),
        frame_case("combined.objective", t, [_r("objectiveScan", True, dz=0.2)]),
        frame_case("combined.zstage", t, [_r("zStageScan", True, reverse=False)]),
        frame_case("combined.resample", t, [_r("resampleFactor", [2, 2, 2])]),
        frame_case("combined.resample_xyz", t, [_r("resampleFactor", [1, 1, 1.5])]),
        frame_case("combined.resample_one", t, [_r("resampleFactor", 1.5)]),
        frame_case("combined.cubic", t, [_r("interpMethod", "cubic")]),
        frame_case("combined.thresh", t, [_r("xStepThresh", 2.5)]),
        frame_case("combined.flatfield", frame_tree(ff=True), [_r(*ffa)]),
        frame_case(
            "combined.flatfield_offset",
            frame_tree(ff=True),
            [_r(*ffa, "constOffset", 10, "lowerLimit", 0.6)],
        ),
        frame_case(
            "combined.bbox_flip",
            t,
            [_r("inputBbox", [5, 3, 2, 33, 58, 27], "flipZstack", True)],
        ),
        frame_case("combined.in_xyz", t, [_r("inputAxisOrder", "xyz")]),
        frame_case("combined.in_zyx", t, [_r("inputAxisOrder", "zyx")]),
        frame_case("combined.out_zyx", t, [_r("outputAxisOrder", "zyx")]),
        frame_case("combined.out_xzy", t, [_r("outputAxisOrder", "xzy", "saveZarr", True)]),
        frame_case(
            "combined.parts",
            frame_tree(parts=3),
            [_r()],
            frames=f"{FSN}.tif|{FSN}_part0001.tif|{FSN}_part0002.tif",
        ),
        frame_case("combined.zarr_in", frame_tree(zarr=True), [_r()], frames=f"{FSN}.zarr"),
        frame_case(
            "combined.exists",
            frame_tree(existing=(f"DSR/{FSN}.tif",)),
            [_r()],
        ),
        # without save3DStack the results are empty files (compared as such)
        frame_case("combined.no3d", t, [_r("save3DStack", False)]),
        frame_case("combined.no3d_zarr", t, [_r("save3DStack", False, "saveZarr", True)]),
        frame_case("combined.nomip", t, [_r("saveMIP", False)]),
        frame_case(
            "combined.names",
            t,
            [_r("DSRDirName", "rotated", "surffix", "_v2", "uuid", "0123abcd")],
        ),
        frame_case("combined.rescale_single", t, [_r("rescaleRotate", True, "save16bit", False)]),
        frame_case("combined.rescale_single_in", frame_tree("single"), [_r("rescaleRotate", True)]),
        frame_case(
            "combined.rescale_u16",
            t,
            [_r("rescaleRotate", True)],
            expect_error="MATLAB:mixedClasses",
            py_error=df.MatlabError,
        ),
        frame_case(
            "combined.bkremoval",
            frame_tree(ff=True),
            [_r("BKRemoval", True, "backgroundImage", "@d/bg.tif")],
            expect_error="MATLAB:UndefinedFunction",
            py_error=df.MatlabError,
        ),
        frame_case("combined.missing", t, [_r()], frames=f"none.tif|{FSN}.tif"),
        # separate
        frame_case("separate", t, [_r(*sep)]),
        frame_case("separate.zarr", t, [_r(*sep, "saveZarr", True)]),
        frame_case("separate.single_out", t, [_r(*sep, "save16bit", False)]),
        frame_case("separate.crop", t, [_r(*sep, "crop", True)]),
        frame_case("separate.deskew_only", t, [_r(*ds_only)]),
        frame_case("separate.objective", t, [_r(*sep, "objectiveScan", True, dz=0.2)]),
        frame_case("separate.zstage", t, [_r(*sep, "zStageScan", True, reverse=False)]),
        frame_case("separate.resample", t, [_r(*sep, "resampleFactor", [2, 2, 2])]),
        frame_case("separate.rescale", t, [_r(*sep, "rescaleRotate", True)]),
        frame_case("separate.out_zyx", t, [_r(*sep, "outputAxisOrder", "zyx")]),
        frame_case(
            "separate.two_runs",
            t,
            [
                _r(*ds_only, "inputAxisOrder", "xyz", "inputBbox", [5, 3, 2, 33, 58, 27]),
                _r(*sep, "inputAxisOrder", "xyz", "inputBbox", [5, 3, 2, 33, 58, 27]),
            ],
        ),
        frame_case(
            "separate.two_runs_single",
            t,
            [_r(*ds_only, "save16bit", False), _r(*sep, "save16bit", False)],
        ),
        frame_case(
            "separate.ds_exists",
            frame_tree(existing=(f"DS/{FSN}.tif",)),
            [_r(*ds_only)],
        ),
        frame_case(
            "separate.overwrite",
            frame_tree(existing=(f"DS/{FSN}.tif",)),
            [_r(*ds_only, "overwrite", True)],
        ),
        frame_case("separate.no3d", t, [_r(*ds_only, "save3DStack", False)]),
    ]
    # above 1000 slices: blocks of rows (XRDeskewRotateFrameTest/largeZSeparatePathWritesDSMIP);
    # MATLAB writes that DS file through blockedImage's MPageTiffAdapter, another tiff layout
    # than writetiff's (the port writes every tiff with cpp-tiff): voxels compared
    tz = frame_tree(shape=(1001, 16, 16))
    for cid, extra in (
        ("large_z", ()),
        ("large_z.single_out", ("save16bit", False)),
        ("large_z.rotate", ("rotate", True)),
    ):
        args = [PX, 0.01, "skewAngle", SKEW, "DSRCombined", False, "rotate", False, *extra]
        if "rotate" in extra:
            args = [PX, 0.01, "skewAngle", SKEW, "DSRCombined", False, *extra]
        cases.append(
            frame_case(f"separate.{cid}", tz, [args], marks=("slow",), tree=tree(tiff_bytes=False))
        )
    return cases


# ------------------------------------------------------------------------------------------
# scaleContrast cases
# ------------------------------------------------------------------------------------------


def _values(kind: str, n: int = 4000, seed: int = 0):
    def make():
        rng = CounterRNG(seed, f"scale/{kind}/{n}")
        if kind == "u16":
            return rng.integers(0, 3000, (n,), dtype=np.uint16)
        if kind == "s":
            return (rng.uniform((n,), np.float32) * np.float32(3000)).astype(np.float32)
        if kind == "s3":
            return to_matlab((rng.uniform((5, 7, 9), np.float32) - np.float32(0.5)) * 40)
        if kind == "snan":
            a = (rng.uniform((n,), np.float32) * np.float32(50)).astype(np.float32)
            a[3] = np.float32(np.inf)
            a[7] = -np.float32(np.inf)
            from petakit5d.core.dtypes import NAN

            a[11] = np.float32(NAN)
            return a
        return rng.uniform((n,), np.float64) * 2000 - 1000

    return Value(make)


def sc_case(cid, kind, *args, vol=False, **kw):
    a = Vol(In("a")) if vol else In("a")
    return Case(
        f"scale.{cid}",
        matlab="scaleContrast",
        args=[In("a"), *args],
        python=df.scale_contrast,
        py_args=[a, *args],
        out_roles=("vol",) if vol else (None,),
        inputs={"a": _values(kind)},
        **kw,
    )


def _f32(*v):
    return np.array(v, dtype=np.float32)


def _scale_cases() -> list[Case]:
    err = {"expect_error": "MATLAB:mixedClasses", "py_error": df.MatlabError}
    return [
        sc_case("u16_double", "u16", [100.0, 1300.0], [0, 65535]),
        sc_case("u16_double_odd", "u16", [37.5, 2999.25], [3, 70001]),
        sc_case("u16_single_range", "u16", _f32(100, 1300), [0, 65535]),
        sc_case("u16_u16_range", "u16", np.array([1, 3000], dtype=np.uint16), [0, 65535], **err),
        sc_case("u16_u16_const", "u16", np.array([3, 3], dtype=np.uint16), [0, 65535]),
        sc_case("s_single", "s", _f32(0.25, 2999.5), [0, 65535]),
        sc_case("s_inside", "s", _f32(100, 2000), [0, 65535]),
        sc_case("s_double_range", "s", [0.25, 2999.5], [0, 65535]),
        sc_case("s_out_1020", "s", _f32(0.25, 2999.5), [10, 20]),
        sc_case("s_out_rev", "s", _f32(0.25, 2999.5), [1, 0]),
        sc_case("s_out_single", "s", _f32(1, 40), _f32(0, 65535)),
        sc_case("s_reversed", "s", _f32(40, 10), [0, 65535]),
        sc_case("s_const", "s", _f32(3.5, 3.5), [10, 20]),
        sc_case("s_signed0", "s", np.array([-0.0, 0.0], dtype=np.float32), [0, 65535]),
        sc_case("s_inf_end", "s", _f32(0, np.inf), [0, 65535]),
        sc_case("s_inf_both", "s", _f32(np.inf, np.inf), [0, 65535]),
        sc_case("s_nan_range", "s", _f32(np.nan, 3), [0, 65535]),
        sc_case("s_naninf", "snan", _f32(0.5, 49), [0, 65535]),
        sc_case("s_3d", "s3", _f32(-19.5, 19.75), [0, 65535], vol=True),
        sc_case("d_double", "d", [-1000.0, 1000.0], [0, 65535]),
        sc_case("d_inside", "d", [-100.0, 100.0], [0, 65535]),
        sc_case("d_col_range", "d", np.array([[0.2], [0.7]]), np.array([[0.0], [65535.0]])),
        sc_case("d_const", "d", [5.0, 5.0], [10, 20]),
        # r2 - r1 saturates to 0 in uint16: the degenerate range, double zeros
        sc_case("u16_reversed_u16", "u16", np.array([5, 3], dtype=np.uint16), [0, 65535]),
        sc_case("s_default_out", "s", _f32(10, 50)),
        Case(
            "scale.s_defaults",
            matlab="scaleContrast",
            args=[In("a")],
            python=df.scale_contrast,
            inputs={"a": _values("s")},
        ),
    ]


# ------------------------------------------------------------------------------------------
# XR_estimateDeskewRotateMemory: the output size
# ------------------------------------------------------------------------------------------


def p_memory_out_size(im_size, px, dz, *nv):
    kw = {}
    names = {
        "skewAngle": "skew_angle",
        "objectiveScan": "objective_scan",
        "zStageScan": "z_stage_scan",
        "resampleFactor": "resample_factor",
        "interpMethod": "interp_method",
        "outSize": "out_size",
        "largeFile": "large_file",
        "DSRCombined": "dsr_combined",
        "rotate": "rotate",
    }
    for k in range(0, len(nv), 2):
        name, v = nv[k], nv[k + 1]
        if name in ("outSize",):
            v = (v[2], v[0], v[1])
        if name == "resampleFactor":
            v = tuple(v) if isinstance(v, list | tuple) else (v,)
            v = v[::-1] if len(v) == 3 else (v if len(v) == 1 else (v[1], v[0]))
        kw[names[name]] = v
    sz = (int(im_size[2]), int(im_size[0]), int(im_size[1]))
    mem, arr, osz = dm.estimate_deskew_rotate_memory(sz, px, dz, **kw)
    return [mem, arr, tuple(float(v) for v in osz)]


def mem_case(cid, im_size, dz, *nv):
    return Case(
        f"memory.{cid}",
        matlab="XR_estimateDeskewRotateMemory",
        args=[list(im_size), PX, dz, "skewAngle", SKEW, *nv],
        nargout=3,
        python=p_memory_out_size,
        out_roles=(None, None, "yxz"),
        compare={
            "out1": ignore("the port's memory is its own (measured), not MATLAB's"),
            "out2": ignore("the port's arrays are its own (measured), not MATLAB's"),
        },
    )


def _memory_cases() -> list[Case]:
    s = [64, 96, 80]
    return [
        mem_case("dz03", s, 0.3),
        mem_case("dz02", s, 0.2),
        mem_case("dz08", s, 0.8),
        mem_case("objective", s, 0.3, "objectiveScan", True),
        mem_case("resample222", s, 0.3, "resampleFactor", [2, 2, 2]),
        mem_case("resample1115", s, 0.3, "resampleFactor", [1, 1, 1.5]),
        mem_case("resample_scalar", s, 0.3, "resampleFactor", 2),
        mem_case("resample_two", s, 0.3, "resampleFactor", [1.5, 2]),
        mem_case("resample_xyz", s, 0.3, "resampleFactor", [1.2, 1.5, 2]),
        mem_case("cubic", s, 0.3, "interpMethod", "cubic"),
        mem_case("zstage", s, 0.3, "zStageScan", True),
        mem_case("demo", [1800, 512, 500], 0.3),
        mem_case("preresample", [16, 24, 1001], 0.3, "resampleFactor", [2, 2, 2]),
        mem_case("preresample_big", [1800, 512, 2000], 0.3, "resampleFactor", [2, 2, 2]),
        mem_case("outsize", [18, 96, 80], 0.3, "outSize", [16, 120, 47], "largeFile", True),
        mem_case("separate", [1800, 512, 500], 0.3, "DSRCombined", False),
    ]


# ------------------------------------------------------------------------------------------
# rotateFrame3D's Crop (another author's GPL function, a black box): the output size from the
# data's zero regions (deskew_rotate.volume._crop_out_size)
# ------------------------------------------------------------------------------------------


def make_band(shape, step, width, start, kind="single", reverse=False, seed=0, fill=None):
    """MATLAB's form of a (z, y, x) volume whose slice k holds data in x from
    start + floor(k * step) (from the other end with ``reverse``) over ``width`` columns, 0
    elsewhere: a deskewed volume's zero triangles"""

    def make():
        nz, ny, nx = shape
        rng = CounterRNG(seed, f"band/{kind}/{shape}/{step}/{width}/{start}/{reverse}")
        vals = rng.uniform(shape, np.float32) * np.float32(1000) + np.float32(1)
        v = np.zeros(shape, dtype=np.float32)
        for k in range(nz):
            x0 = start + int(np.floor(k * step))
            if reverse:
                x0 = nx - x0 - width
            v[k, :, x0 : x0 + width] = vals[k, :, x0 : x0 + width] if fill is None else fill
        if kind == "uint16":
            v = np.floor(v).astype(np.uint16)
        elif kind == "double":
            v = v.astype(np.float64)
        return to_matlab(v)

    return Value(make)


def make_fringe(mode, shape=(24, 19, 70), step=1.5627, width=30, start=2):
    """a band of positive values with a column of negative values on both ends of each slice,
    as a cubic deskew's ringing: in every row (``'neg'``: the columns' maximum along y is
    negative) or in every other row with zeros between (``'mixed'``: the maximum is 0)"""

    def make():
        nz, ny, nx = shape
        rng = CounterRNG(7, f"fringe/{mode}/{shape}")
        vals = rng.uniform(shape, np.float32) * np.float32(1000) + np.float32(1)
        neg = -(rng.uniform(shape, np.float32) * np.float32(0.5) + np.float32(0.01))
        v = np.zeros(shape, dtype=np.float32)
        rows = slice(None) if mode == "neg" else slice(0, None, 2)
        for k in range(nz):
            x0 = start + int(np.floor(k * step))
            v[k, :, x0 : x0 + width] = vals[k, :, x0 : x0 + width]
            for x in (x0 - 1, x0 + width):
                v[k, rows, x] = neg[k, rows, x]
        return to_matlab(v)

    return Value(make)


def make_box(shape, box, kind="single", seed=0):
    def make():
        rng = CounterRNG(seed, f"rbox/{shape}/{box}")
        v = np.zeros(shape, dtype=np.float32)
        z0, y0, x0, z1, y1, x1 = box
        v[z0:z1, y0:y1, x0:x1] = rng.uniform(shape, np.float32)[z0:z1, y0:y1, x0:x1] + 1
        return to_matlab(v.astype(np.uint16) if kind == "uint16" else v)

    return Value(make)


def rot_case(cid, make, angle, zx, reverse, *nv, py=None, **kw):
    return Case(
        f"rotcrop.{cid}",
        matlab="rotateFrame3D",
        args=[Vol(In("v")), angle, zx, reverse, *nv],
        python="petakit5d.deskew_rotate.volume:rotate_volume",
        py_args=[Vol(In("v")), angle, zx, reverse],
        py_kwargs=py or {},
        out_roles=("vol",),
        inputs={"v": make},
        **kw,
    )


def _rotcrop_cases() -> list[Case]:
    z2 = 0.9937  # sind(32.45) 0.2 / 0.108, rounded
    s2 = 1.5627  # the deskew's x step for dz 0.2
    sh = (24, 19, 70)
    cases = []
    for cls in ("single", "uint16", "double"):
        cases.append(
            rot_case(
                f"band.fwd.{cls}",
                make_band(sh, s2, 30, 0, cls),
                SKEW,
                z2,
                False,
                "Crop",
                True,
                py={"crop": True},
            )
        )
    for cid, make, ang, zx, rev in (
        ("band.rev", make_band(sh, s2, 30, 0, reverse=True), SKEW, z2, True),
        ("band.rev_fwd_data", make_band(sh, s2, 30, 0, reverse=True), SKEW, z2, False),
        ("band.interior", make_band(sh, s2, 30, 3), SKEW, z2, False),
        ("band.step23", make_band((20, 11, 80), 2.344, 26, 1), SKEW, 1.4906, False),
        ("band.angle", make_band(sh, s2, 30, 0), -31.5, z2, False),
        ("band.zx2", make_band(sh, s2, 30, 0), SKEW, 2.0, False),
        ("band.nan", make_band(sh, s2, 30, 2, fill=np.float32(np.nan)), SKEW, z2, False),
        ("box.inside", make_box((21, 13, 61), (3, 2, 9, 18, 13, 50)), SKEW, z2, False),
        ("box.edge", make_box((21, 13, 61), (0, 0, 0, 21, 13, 40)), SKEW, z2, False),
        ("box.off", make_box((21, 13, 61), (2, 4, 20, 15, 9, 61), "uint16"), SKEW, z2, False),
        ("full", make_box((21, 13, 61), (0, 0, 0, 21, 13, 61)), SKEW, z2, False),
        # a column is data where its maximum along y is not 0
        ("fringe.mixed", make_fringe("mixed"), SKEW, z2, False),
        ("fringe.neg", make_fringe("neg"), SKEW, z2, False),
    ):
        cases.append(rot_case(cid, make, ang, zx, rev, "Crop", True, py={"crop": True}))
    band = make_band(sh, s2, 30, 0)
    cases += [
        rot_case("band.nocrop", band, SKEW, z2, False, "Crop", False, py={"crop": False}),
        rot_case("band.default", band, SKEW, z2, False),
        rot_case(
            "band.resample",
            band,
            SKEW,
            z2,
            False,
            "Crop",
            True,
            "resample",
            [1.5, 1, 2],
            py={"crop": True, "resample": (2, 1, 1.5)},
        ),
        rot_case(
            "band.outsize",
            band,
            SKEW,
            z2,
            False,
            "Crop",
            True,
            "outSize",
            [19, 60, 20],
            py={"crop": True, "out_size": (20, 19, 60)},
        ),
        rot_case(
            "band.cubic",
            band,
            SKEW,
            z2,
            False,
            "Crop",
            True,
            "Interp",
            "cubic",
            py={"crop": True, "interp": "cubic"},
        ),
        # an objective scan is never cropped (RLdecon's objective-scan rotation)
        rot_case(
            "band.objective",
            band,
            SKEW,
            z2,
            False,
            "Crop",
            True,
            "objectiveScan",
            True,
            py={"crop": True, "objective_scan": True},
        ),
        rot_case(
            "band.objective_cubic",
            band,
            SKEW,
            z2,
            True,
            "Crop",
            True,
            "objectiveScan",
            True,
            "Interp",
            "cubic",
            py={"crop": True, "objective_scan": True, "interp": "cubic"},
        ),
        rot_case(
            "band.objective_false",
            band,
            SKEW,
            z2,
            False,
            "Crop",
            True,
            "objectiveScan",
            False,
            py={"crop": True, "objective_scan": False},
        ),
        rot_case(
            "err.voxel",
            make_box((21, 13, 61), (10, 6, 30, 11, 7, 31)),
            SKEW,
            z2,
            False,
            "Crop",
            True,
            py={"crop": True},
            expect_error="MATLAB:imref3d:expectedPositive",
            py_error=ValueError,
        ),
    ]

    def band_tree(folder):
        v = np.zeros((24, 19, 70), dtype=np.uint16)
        rng = CounterRNG(4, "rotfile")
        vals = rng.integers(1, 3000, v.shape, dtype=np.uint16)
        for k in range(24):
            x0 = int(np.floor(k * s2))
            v[k, :, x0 : x0 + 30] = vals[k, :, x0 : x0 + 30]
        synth.write_tiff(Path(folder) / "scan_ch0.tif", v)

    cases.append(
        Case(
            "rotcrop.rotate_file",
            code="XR_RotateFrame3D([in{1}, 'scan_ch0.tif'], in{2 : end}); out = {};",
            nargout=0,
            args=[InPath("d"), PX, 0.2, "SkewAngle", SKEW, "uuid", "u0"],
            python=p_rotate_file,
            inputs={"d": Tree(band_tree, "d")},
        )
    )
    return cases


def p_rotate_file(d, *args):
    from petakit5d import matlab

    matlab.XR_RotateFrame3D(d + "scan_ch0.tif", *args)
    return []


CASES = _frame_cases() + _scale_cases() + _memory_cases() + _rotcrop_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# XRDeskewRotateFrameTest in intent (no MATLAB)
# ------------------------------------------------------------------------------------------


def _write_frame(tmp_path, kind="uint16", shape=RAW, name=FSN) -> tuple[str, np.ndarray]:
    frame = make_frame(shape, kind)
    fn = tmp_path / f"{name}.tif"
    synth.write_tiff(fn, frame)
    return str(fn), frame


def _read(path):
    from petakit5d.io import read_tiff, read_zarr

    a = read_zarr(path) if str(path).endswith(".zarr") else read_tiff(path)
    return a


def _expected_dsr(frame, dz=0.3, **kw):
    from petakit5d.deskew_rotate import deskew_rotate_volume

    kw.setdefault("save_16bit", True)
    return deskew_rotate_volume(frame, SKEW, dz, PX, reverse=True, **kw)


def _mip(vol):
    return np.max(vol, axis=0)


def test_combined_writes_dsr_and_mip(tmp_path):
    fn, frame = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, skew_angle=SKEW, reverse=True)
    expected = _expected_dsr(frame)
    dsr = _read(tmp_path / "DSR" / f"{FSN}.tif")
    np.testing.assert_array_equal(dsr, expected)
    assert dsr.dtype == np.uint16
    np.testing.assert_array_equal(
        _read(tmp_path / "DSR" / "MIPs" / f"{FSN}_MIP_z.tif"), _mip(expected)
    )
    assert not (tmp_path / "DS").exists()


def test_save_zarr_matches_tiff(tmp_path):
    fn, _ = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True)
    os.rename(tmp_path / "DSR", tmp_path / "tif")
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, save_zarr=True)
    np.testing.assert_array_equal(
        _read(str(tmp_path / "DSR" / f"{FSN}.zarr")), _read(tmp_path / "tif" / f"{FSN}.tif")
    )
    np.testing.assert_array_equal(
        _read(str(tmp_path / "DSR" / "MIPs" / f"{FSN}_MIP_z.zarr")),
        _read(tmp_path / "tif" / "MIPs" / f"{FSN}_MIP_z.tif"),
    )


def test_save16bit_false_writes_single(tmp_path):
    fn, frame = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, save_16bit=False)
    dsr = _read(tmp_path / "DSR" / f"{FSN}.tif")
    assert dsr.dtype == np.float32
    np.testing.assert_array_equal(dsr, _expected_dsr(frame, save_16bit=False))


def test_separate_path_writes_deskewed_and_rotated(tmp_path):
    from petakit5d.core import dtypes
    from petakit5d.deskew_rotate import deskew_volume

    fn, frame = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, dsr_combined=False, rotate=True)
    ds = _read(tmp_path / "DS" / f"{FSN}.tif")
    expected = dtypes.uint16(deskew_volume(frame.astype(np.float32), SKEW, 0.3, PX, True))
    np.testing.assert_array_equal(ds, expected)
    assert (tmp_path / "DS" / "MIPs" / f"{FSN}_MIP_z.tif").is_file()
    # the same geometry as the combined transform (beads within 0.15 voxel)
    dsr = _read(tmp_path / "DSR" / f"{FSN}.tif").astype(np.float32)
    combined = _expected_dsr(frame).astype(np.float32)
    assert dsr.shape == combined.shape
    peaks = _peaks(combined)
    got = synth.centroids(dsr - 100, peaks, 3)
    ref = synth.centroids(combined - 100, peaks, 3)
    assert np.abs(got - ref).max() < 0.15


def _peaks(vol):
    """the brightest local maxima (z, y, x) of a bead volume"""
    from scipy import ndimage

    mx = ndimage.maximum_filter(vol, size=5)
    mask = (vol == mx) & (vol > 0.5 * vol.max())
    pos = np.argwhere(mask)
    inner = np.all((pos >= 4) & (pos < np.array(vol.shape) - 4), axis=1)
    return pos[inner]


def test_deskew_only_writes_only_ds(tmp_path):
    fn, _ = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, dsr_combined=False, rotate=False)
    assert (tmp_path / "DS" / f"{FSN}.tif").is_file()
    assert not (tmp_path / "DSR").exists()


def test_objective_scan_separate_path_runs(tmp_path):
    from petakit5d.io import get_image_size

    fn, _ = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.2, objective_scan=True, dsr_combined=False, rotate=True)
    theta = SKEW * np.pi / 180
    nz, ny, nx = RAW
    q = 0.2 / PX
    from petakit5d.core.dtypes import mround

    out = (
        int(mround(nz * q * np.cos(theta) + nx * np.sin(theta))),
        ny,
        int(mround(nx * np.cos(theta) + nz * q * np.sin(theta))),
    )
    assert tuple(get_image_size(tmp_path / "DSR" / f"{FSN}.tif")) == out


@pytest.mark.slow
def test_large_z_separate_path_writes_ds_mip(tmp_path):
    from petakit5d.core import dtypes
    from petakit5d.deskew_rotate import deskew_volume

    rng = CounterRNG(3, "largez")
    frame = (100 + rng.integers(0, 1000, (1001, 16, 16), dtype=np.uint16)).astype(np.uint16)
    fn = tmp_path / f"{FSN}.tif"
    synth.write_tiff(fn, frame)
    df.deskew_rotate_file(str(fn), PX, 0.01, dsr_combined=False, rotate=False)
    expected = dtypes.uint16(deskew_volume(frame.astype(np.float32), SKEW, 0.01, PX, False))
    np.testing.assert_array_equal(_read(tmp_path / "DS" / f"{FSN}.tif"), expected)
    np.testing.assert_array_equal(
        _read(tmp_path / "DS" / "MIPs" / f"{FSN}_MIP_z.tif"), _mip(expected)
    )


def test_flat_field_correction_is_applied_first(tmp_path):
    from petakit5d.flatfield import flatfield_correct_volume

    fn, frame = _write_frame(tmp_path)
    bg = np.zeros(RAW[1:], dtype=np.float32) + 100
    ff = bg + np.linspace(500, 1000, RAW[2]).astype(np.float32)
    synth.write_tiff(tmp_path / "ff.tif", ff)
    synth.write_tiff(tmp_path / "bg.tif", bg)
    df.deskew_rotate_file(
        fn,
        PX,
        0.3,
        reverse=True,
        ff_correction=True,
        ff_image=str(tmp_path / "ff.tif"),
        background_image=str(tmp_path / "bg.tif"),
        lower_limit=0.4,
    )
    corrected = flatfield_correct_volume(frame, ff, bg, lower_limit=0.4)
    np.testing.assert_array_equal(_read(tmp_path / "DSR" / f"{FSN}.tif"), _expected_dsr(corrected))


def test_input_bbox_and_flip_are_applied_first(tmp_path):
    fn, frame = _write_frame(tmp_path)
    box = (2, 8, 4, 27, 33, 58)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, input_bbox=box, flip_zstack=True)
    cropped = np.ascontiguousarray(frame[2:27, 8:33, 4:58][::-1])
    np.testing.assert_array_equal(_read(tmp_path / "DSR" / f"{FSN}.tif"), _expected_dsr(cropped))


@pytest.mark.parametrize("order", ["zxy", "xzy"])
def test_input_axis_order_is_permuted_to_zyx(tmp_path, order):
    from petakit5d.core.axes import axis_permutation

    frame = make_frame()
    (tmp_path / order).mkdir()
    fn = tmp_path / order / f"{FSN}.tif"
    synth.write_tiff(fn, np.ascontiguousarray(np.transpose(frame, axis_permutation("zyx", order))))
    df.deskew_rotate_file(str(fn), PX, 0.3, reverse=True, input_axis_order=order)
    np.testing.assert_array_equal(
        _read(tmp_path / order / "DSR" / f"{FSN}.tif"), _expected_dsr(frame)
    )


def test_output_axis_order_permutes_the_result(tmp_path):
    # a 3-cycle ('xzy' in our names is MATLAB's 'zyx'), which its inverse would not give
    fn, frame = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, output_axis_order="xzy")
    expected = np.transpose(_expected_dsr(frame), (2, 0, 1))
    np.testing.assert_array_equal(_read(tmp_path / "DSR" / f"{FSN}.tif"), expected)


def test_existing_output_is_not_recomputed(tmp_path):
    fn, frame = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True)
    first = _read(tmp_path / "DSR" / f"{FSN}.tif")
    synth.write_tiff(fn, frame * 2)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True)
    np.testing.assert_array_equal(_read(tmp_path / "DSR" / f"{FSN}.tif"), first)


def test_z_stage_scan_separate_matches_combined(tmp_path):
    from petakit5d.io import get_image_size

    fn, _ = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, z_stage_scan=True)
    combined = get_image_size(tmp_path / "DSR" / f"{FSN}.tif")
    df.deskew_rotate_file(
        fn, PX, 0.3, z_stage_scan=True, dsr_combined=False, rotate=True, dsr_dir_name="DSR2/"
    )
    assert get_image_size(tmp_path / "DSR2" / f"{FSN}.tif") == combined


def test_rotation_of_saved_deskewed_volume_keeps_crop_and_axis_order(tmp_path):
    frame = make_frame()
    (tmp_path / "xzy").mkdir()
    fn = str(tmp_path / "xzy" / f"{FSN}.tif")
    synth.write_tiff(fn, np.ascontiguousarray(np.transpose(frame, (2, 0, 1))))
    kw = {
        "reverse": True,
        "input_axis_order": "xzy",
        "input_bbox": (4, 2, 8, 58, 27, 33),
        "dsr_combined": False,
        "save_16bit": False,
    }
    df.deskew_rotate_file(fn, PX, 0.3, rotate=False, **kw)
    df.deskew_rotate_file(fn, PX, 0.3, rotate=True, **kw)
    two_runs = _read(tmp_path / "xzy" / "DSR" / f"{FSN}.tif")
    import shutil

    shutil.rmtree(tmp_path / "xzy" / "DS")
    shutil.rmtree(tmp_path / "xzy" / "DSR")
    df.deskew_rotate_file(fn, PX, 0.3, rotate=True, **kw)
    np.testing.assert_array_equal(two_runs, _read(tmp_path / "xzy" / "DSR" / f"{FSN}.tif"))


@pytest.mark.parametrize(
    ("extra", "in_place"),
    [
        ({}, True),
        ({"dsr_combined": False, "rotate": False}, True),
        ({"dsr_combined": False, "rotate": False, "save_16bit": False}, True),
        ({"uuid": "01234567-89ab-cdef-0123-456789abcdef"}, False),
    ],
)
def test_private_zarr_results_are_written_in_place(tmp_path, monkeypatch, extra, in_place):
    import petakit5d.io as pio
    from petakit5d.core import ids

    calls = []
    real = pio.write_zarr

    def spy(data, path, **kw):
        calls.append((str(path), kw.get("uuid", True)))
        return real(data, path, **kw)

    monkeypatch.setattr(pio, "write_zarr", spy)
    fn, _ = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, save_zarr=True, **extra)
    results = [c for c in calls if "MIPs" not in c[0]]
    assert len(results) == 1
    expect_private = in_place and len(ids.get_uuid()) >= 32
    assert results[0][1] is (not expect_private)


def test_background_removal_raises_matlabs_error(tmp_path):
    # KNOWN ISSUE (kept): MATLAB calls XR_CameraBackgroundRemoval, which does not exist
    fn, _ = _write_frame(tmp_path)
    synth.write_tiff(tmp_path / "bg.tif", np.zeros(RAW[1:], dtype=np.float32) + 100)
    with pytest.raises(df.MatlabError) as err:
        df.deskew_rotate_file(
            fn, PX, 0.3, reverse=True, bk_removal=True, background_image=str(tmp_path / "bg.tif")
        )
    assert err.value.identifier == "MATLAB:UndefinedFunction"
    assert not (tmp_path / "DSR" / f"{FSN}.tif").exists()


def test_objective_scan_separate_path_is_deskewed(tmp_path):
    # KNOWN ISSUE (kept as MATLAB behaves): the separate path deskews objective scans before
    # rotating them, so beads are sheared against the pure rotation of the combined cubic path
    from petakit5d.deskew_rotate import deskew_rotate_volume

    fn, frame = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.2, objective_scan=True, dsr_combined=False, rotate=True)
    dsr = _read(tmp_path / "DSR" / f"{FSN}.tif").astype(np.float32)
    expected = deskew_rotate_volume(
        frame.astype(np.float32), SKEW, 0.2, PX, objective_scan=True, interp_method="cubic"
    )
    peaks = _peaks(expected)
    got = synth.centroids(dsr - 100, peaks, 3)
    ref = synth.centroids(expected - 100, peaks, 3)
    assert np.abs(got - ref).max() > 0.5


# ------------------------------------------------------------------------------------------
# the port's own rules
# ------------------------------------------------------------------------------------------


def test_publish_order_and_temporary_names(tmp_path):
    from petakit5d.testing.publish_spy import PublishSpy, normalize

    fn, _ = _write_frame(tmp_path)
    with PublishSpy() as spy:
        df.deskew_rotate_file(fn, PX, 0.3, reverse=True, dsr_combined=False, rotate=True, uuid="u7")
    ds, dsr = tmp_path / "DS" / f"{FSN}.tif", tmp_path / "DSR" / f"{FSN}.tif"
    assert spy.was_moved(ds) and spy.was_moved(dsr)
    assert f"MIPs/{FSN}_MIP_z.tif" in spy.files_when_moved(ds)
    assert f"MIPs/{FSN}_MIP_z.tif" in spy.files_when_moved(dsr)
    srcs = [normalize(c.src) for c in spy.calls]
    assert normalize(tmp_path / "DS" / f"{FSN}_u7.tif") in srcs
    assert normalize(tmp_path / "DSR" / f"{FSN}_u7.tif") in srcs
    assert normalize(tmp_path / "DSR" / "MIPs" / f"{FSN}_MIP_z.tiff_u7") in srcs
    # the DS MIP is written under its final name (as MATLAB)
    assert not any("DS/MIPs" in s and "DSR" not in s for s in srcs)


def test_another_writers_result_is_kept(tmp_path, monkeypatch):
    from petakit5d.jobs import publish

    fn, _ = _write_frame(tmp_path)
    dsr = tmp_path / "DSR" / f"{FSN}.tif"
    real = publish.rename_no_replace

    def other_writer_first(src, dst):
        if os.fspath(dst).endswith(f"{FSN}.tif"):
            synth.write_tiff(dst, np.full((2, 2, 2), 5, dtype=np.uint16))
        return real(src, dst)

    monkeypatch.setattr(publish, "rename_no_replace", other_writer_first)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True)
    assert np.all(_read(dsr) == 5)
    assert sorted(p.name for p in (tmp_path / "DSR").iterdir()) == ["MIPs", f"{FSN}.tif"]


def test_save_3d_stack_false_writes_empty_results(tmp_path):
    fn, _ = _write_frame(tmp_path)
    df.deskew_rotate_file(fn, PX, 0.3, reverse=True, save_3d_stack=False)
    assert (tmp_path / "DSR" / f"{FSN}.tif").stat().st_size == 0
    assert (tmp_path / "DSR" / "MIPs" / f"{FSN}_MIP_z.tif").stat().st_size > 0
    # zarr: an empty file with the zarr's name (as MATLAB); the DS of the separate path too
    df.deskew_rotate_file(fn, PX, 0.3, save_3d_stack=False, save_zarr=True, dsr_dir_name="Z/")
    assert (tmp_path / "Z" / f"{FSN}.zarr").is_file()
    assert (tmp_path / "Z" / f"{FSN}.zarr").stat().st_size == 0
    df.deskew_rotate_file(fn, PX, 0.3, save_3d_stack=False, dsr_combined=False)
    assert (tmp_path / "DS" / f"{FSN}.tif").stat().st_size == 0


def test_missing_input_warns_and_writes_nothing(tmp_path):
    fn, _ = _write_frame(tmp_path)
    with pytest.warns(UserWarning, match="does not exist"):
        df.deskew_rotate_file([str(tmp_path / "none.tif"), fn], PX, 0.3)
    assert sorted(p.name for p in tmp_path.iterdir()) == [f"{FSN}.tif"]


def test_partial_files_are_read_into_one_array(tmp_path):
    frame = make_frame()
    names = []
    for k, (a, b) in enumerate(((0, 10), (10, 11), (11, 29))):
        name = FSN if k == 0 else f"{FSN}_part{k:04d}"
        synth.write_tiff(tmp_path / f"{name}.tif", frame[a:b])
        names.append(str(tmp_path / f"{name}.tif"))
    np.testing.assert_array_equal(df._read_frame(names, 2, None), frame)
    from petakit5d.io import ReadBuffer

    buf = ReadBuffer()
    np.testing.assert_array_equal(df._read_frame(names, 2, buf), frame)
    df.deskew_rotate_file(names, PX, 0.3, reverse=True, read_buffer=buf)
    np.testing.assert_array_equal(_read(tmp_path / "DSR" / f"{FSN}.tif"), _expected_dsr(frame))


def test_read_buffer_is_reused_over_files(tmp_path):
    from petakit5d.io import ReadBuffer

    buf = ReadBuffer()
    for sub in ("a", "b"):
        (tmp_path / sub).mkdir()
        _write_frame(tmp_path / sub)
        df.deskew_rotate_file(str(tmp_path / sub / f"{FSN}.tif"), PX, 0.3, read_buffer=buf)
        first = buf._buf
    assert buf._buf is first
    np.testing.assert_array_equal(
        _read(tmp_path / "a" / "DSR" / f"{FSN}.tif"), _read(tmp_path / "b" / "DSR" / f"{FSN}.tif")
    )


def test_flat_field_cache_follows_the_files(tmp_path):
    fn, frame = _write_frame(tmp_path)
    bg = np.zeros(RAW[1:], dtype=np.float32) + 100
    ff = bg + np.linspace(500, 1000, RAW[2]).astype(np.float32)
    synth.write_tiff(tmp_path / "ff.tif", ff)
    synth.write_tiff(tmp_path / "bg.tif", bg)
    a = df._flat_field(str(tmp_path / "ff.tif"), str(tmp_path / "bg.tif"), RAW[1:], 0.4)
    b = df._flat_field(str(tmp_path / "ff.tif"), str(tmp_path / "bg.tif"), RAW[1:], 0.4)
    assert a is b
    synth.write_tiff(tmp_path / "ff.tif", ff[:, ::-1].copy())
    os.utime(tmp_path / "ff.tif", ns=(1, 1))
    c = df._flat_field(str(tmp_path / "ff.tif"), str(tmp_path / "bg.tif"), RAW[1:], 0.4)
    assert c is not a
    np.testing.assert_array_equal(c[0], a[0][:, ::-1])
    d = df._flat_field(str(tmp_path / "ff.tif"), str(tmp_path / "bg.tif"), RAW[1:], 0.5)
    assert d is not c


def test_thread_counts_give_the_same_files(tmp_path):
    for n in (1, 3):
        (tmp_path / str(n)).mkdir()
        fn, _ = _write_frame(tmp_path / str(n))
        df.deskew_rotate_file(fn, PX, 0.3, reverse=True, num_threads=n, save_zarr=True)
    a = _read(str(tmp_path / "1" / "DSR" / f"{FSN}.zarr"))
    b = _read(str(tmp_path / "3" / "DSR" / f"{FSN}.zarr"))
    np.testing.assert_array_equal(a, b)


def test_matlab_layer_converts_boxes_orders_and_factors(tmp_path):
    from petakit5d import matlab

    fn, frame = _write_frame(tmp_path)
    matlab.XR_deskewRotateFrame(
        fn,
        PX,
        0.3,
        "reverse",
        True,
        "inputBbox",
        [9, 5, 3, 33, 58, 27],
        "BlockSize",
        [16, 32, 8],
        "saveZarr",
        True,
        "resampleFactor",
        [1, 1, 1],
    )
    np.testing.assert_array_equal(
        _read(str(tmp_path / "DSR" / f"{FSN}.zarr")),
        _expected_dsr(np.ascontiguousarray(frame[2:27, 8:33, 4:58]), resample_factor=(1, 1, 1)),
    )
    from petakit5d.io import get_zarr_info

    assert tuple(get_zarr_info(str(tmp_path / "DSR" / f"{FSN}.zarr")).block_size) == (8, 16, 32)


def test_resample_factors_are_completed_as_matlab():
    # MATLAB's [ones(1, 4 - numel(rs)) * rs(1), rs(2:end)] of [x y z], in (z, y, x)
    assert df._resample_completed(None) is None
    assert df._resample_completed((2,)) == (2.0, 2.0, 2.0)
    assert df._resample_completed((3, 1.5)) == (3.0, 1.5, 1.5)  # [xy z] = [1.5 3]
    assert df._resample_completed((3, 2, 1)) == (3.0, 2.0, 1.0)


def test_an_empty_zarr_result_is_made_again(tmp_path, monkeypatch):
    # MATLAB checks a zarr result with exist(..., 'dir'): the empty file that save3DStack false
    # leaves is not one, so the next call computes it again (and keeps the file)
    fn, _ = _write_frame(tmp_path)
    calls = []
    real = df._combined
    monkeypatch.setattr(df, "_combined", lambda *a, **k: calls.append(1) or real(*a, **k))
    for _ in range(2):
        df.deskew_rotate_file(fn, PX, 0.3, save_3d_stack=False, save_zarr=True)
    assert len(calls) == 2
    df.deskew_rotate_file(fn, PX, 0.3)
    df.deskew_rotate_file(fn, PX, 0.3)
    assert len(calls) == 3  # a tiff result exists: not computed again


def test_pre_resampling_is_not_ported_yet(tmp_path):
    frame = np.zeros((1001, 16, 24), dtype=np.uint16) + 100
    fn = tmp_path / f"{FSN}.tif"
    synth.write_tiff(fn, frame)
    with pytest.raises(NotImplementedError, match="imresize"):
        df.deskew_rotate_file(str(fn), PX, 0.3, reverse=True, resample_factor=(2, 2, 2))


# -- rotateFrame3D's Crop ----------------------------------------------------------------------


def test_rotate_crop_is_the_transform_with_a_smaller_centred_size():
    from petakit5d.deskew_rotate import volume as dv

    v = np.zeros((24, 9, 70), dtype=np.float32)
    for k in range(24):
        v[k, :, int(k * 1.5627) : int(k * 1.5627) + 30] = 1 + k
    size = dv._crop_out_size(v, SKEW, 0.9937, False)
    full = dv.rotate_volume(v, SKEW, 0.9937, crop=False)
    cropped = dv.rotate_volume(v, SKEW, 0.9937, crop=True)
    assert cropped.shape == size and size[0] < full.shape[0] and size[2] < full.shape[2]
    np.testing.assert_array_equal(
        cropped, dv.rotate_volume(v, SKEW, 0.9937, crop=False, out_size=size)
    )
    # out_size wins over crop
    assert dv.rotate_volume(v, SKEW, 0.9937, crop=True, out_size=(5, 9, 6)).shape == (5, 9, 6)
    # data without zeros: no crop; a single voxel: MATLAB's imref3d error
    assert dv._crop_out_size(np.ones((4, 3, 8), np.float32), SKEW, 1.0, False) is None
    one = np.zeros((9, 3, 9), np.float32)
    one[4, 1, 4] = 1
    with pytest.raises(ValueError) as err:
        dv.rotate_volume(one, SKEW, 1.0, crop=True)
    assert err.value.identifier == "MATLAB:imref3d:expectedPositive"
    # data reaching the x edge in every slice has no left edge: the width collapses
    edge = np.zeros((21, 3, 61), np.float32)
    edge[:, :, :40] = 1
    inner = np.zeros((21, 3, 61), np.float32)
    inner[:, :, 1:40] = 1
    assert (
        dv._crop_out_size(edge, SKEW, 1.0, False)[2] < dv._crop_out_size(inner, SKEW, 1.0, False)[2]
    )


# -- scale_contrast -----------------------------------------------------------------------------


def test_scale_contrast_rules():
    a = np.array([10, 20, 30], dtype=np.float32)
    out = df.scale_contrast(a, np.array([10, 30], dtype=np.float32), (0, 65535))
    assert out.dtype == np.float32
    np.testing.assert_array_equal(out, np.array([0, 32767.5, 65535], dtype=np.float32))
    assert df.scale_contrast(a, (10.0, 30.0), (0, 1)).dtype == np.float64
    z = df.scale_contrast(a, np.array([5, 5], dtype=np.float32), (10, 20))
    assert z.dtype == np.float64 and not z.any()
    with pytest.raises(df.MatlabError) as err:
        df.scale_contrast(np.array([1, 2], dtype=np.uint16), np.array([1, 2], dtype=np.uint16))
    assert err.value.identifier == "MATLAB:mixedClasses"
    # not clipped
    np.testing.assert_array_equal(df.scale_contrast(a, (20.0, 30.0), (0, 10)), [-10.0, 0.0, 10.0])


@pytest.mark.parametrize("n", [1, 3])
def test_chunked_conversions_equal_matlabs_cast(n):
    from petakit5d.core import dtypes

    # one element past whole chunks: the last chunk must not take the scalar rule
    rng = CounterRNG(9, "chunks")
    size = 2 * df._CHUNK + 1
    a = (rng.uniform((size,), np.float32) * np.float32(70000) - np.float32(1000)).astype(np.float32)
    a[:4] = np.float32([0.49999997, 2.5, -0.5, 65535.5])
    a[-1] = np.float32(0.49999997)
    vol = a.reshape(1, 1, size)
    got = df._to_class(vol, "uint16", n)
    np.testing.assert_array_equal(got, dtypes.to_class(vol, "uint16"))
    assert got[0, 0, -1] == 1  # floor(|x| + 0.5) in single, MATLAB's array rule
    np.testing.assert_array_equal(df._to_class(vol, "double", n), vol.astype(np.float64))


def test_scale_contrast_slabs_equal_one_pass():
    rng = CounterRNG(5, "slabs")
    a = (rng.uniform((40, 300, 400), np.float32) * 1000).astype(np.float32)
    r = (np.float32(a.min()), np.float32(a.max()))
    whole = df.scale_contrast(a, r, (0, 65535))
    small = np.concatenate([df.scale_contrast(a[k : k + 1], r, (0, 65535)) for k in range(40)])
    np.testing.assert_array_equal(whole, small)


# ------------------------------------------------------------------------------------------
# DeskewRotateMemoryTest in intent, and the port's memory model
# ------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("dz", "kw"),
    [
        (0.3, {}),
        (0.2, {}),
        (0.8, {}),
        (0.3, {"objective_scan": True}),
        (0.3, {"resample_factor": (2, 2, 2)}),
        (0.3, {"resample_factor": (1.5, 1, 1)}),
        (0.3, {"interp_method": "cubic"}),
    ],
)
def test_output_size_matches_deskew_rotate_volume(dz, kw):
    from petakit5d.deskew_rotate import deskew_rotate_volume

    vol = np.zeros((80, 64, 96), dtype=np.uint16)
    dsr = deskew_rotate_volume(vol, SKEW, dz, PX, reverse=True, save_16bit=True, **kw)
    _, _, out = dm.estimate_deskew_rotate_memory(vol.shape, PX, dz, skew_angle=SKEW, **kw)
    assert out == dsr.shape


def test_z_stage_scan_uses_the_complementary_angle():
    from petakit5d.deskew_rotate import deskew_rotate_volume

    vol = np.zeros((80, 64, 96), dtype=np.uint16)
    dsr = deskew_rotate_volume(vol, 90 - SKEW, 0.3, PX, reverse=True, save_16bit=True)
    _, _, out = dm.estimate_deskew_rotate_memory(vol.shape, PX, 0.3, z_stage_scan=True)
    assert out == dsr.shape


def test_kernel_path_holds_input_and_output():
    # the port's fused kernel interpolates in line caches: no interpolated volume
    sz = (500, 1800, 512)
    for dz in (0.2, 0.3, 0.8):
        _, arrays, out = dm.estimate_deskew_rotate_memory(sz, PX, dz)
        n_in, n_out = int(np.prod(sz)), int(np.prod(out))
        turn = max(2**26, max(2 * n_in, 2 * n_out) / 16)
        assert arrays == pytest.approx((2 * n_in + 2 * n_out + turn) / 2**30, rel=1e-12)


def test_estimate_grows_with_the_options():
    sz = (500, 1800, 512)
    base, _, _ = dm.estimate_deskew_rotate_memory(sz, PX, 0.3)
    for kw in (
        {"save_16bit": False},
        {"save_zarr": True},
        {"dsr_combined": False},
        {"interp_method": "cubic"},
        {"objective_scan": True},
        {"output_axis_order": "xzy"},
        {"rescale_rotate": True},
    ):
        assert dm.estimate_deskew_rotate_memory(sz, PX, 0.3, **kw)[0] > base, kw
    assert dm.estimate_deskew_rotate_memory(sz, PX, 0.3, read_size=(500, 7200, 512))[0] > base


# peak resident memory (GB) of deskew_rotate_file processes of the port (the demo tile cropped or
# repeated along z; job ids in the task report): size (z, y, x), dz, options, the largest peak
# over 1 and 24 threads (c24 nodes, 2026-11-13)
_Z = {k: (k, 1800, 512) for k in (100, 500, 1000, 2000, 5000)}
MEASURED = {
    "tiny": ((80, 64, 96), 0.3, {}, 0.05),
    "default100": (_Z[100], 0.3, {}, 0.91),
    "default500": (_Z[500], 0.3, {}, 2.56),
    "default2000": (_Z[2000], 0.3, {}, 8.92),
    "default5000": (_Z[5000], 0.3, {}, 21.66),
    "wide250": ((250, 1800, 1024), 0.3, {}, 3.87),
    "dz02": (_Z[2000], 0.2, {}, 7.24),
    "dz05": (_Z[500], 0.5, {}, 3.39),
    "dz08": (_Z[500], 0.8, {}, 4.66),
    "single": (_Z[2000], 0.3, {"save_16bit": False}, 14.36),
    "zarrOutput100": (_Z[100], 0.3, {"save_zarr": True}, 1.64),
    "zarrOutput2000": (_Z[2000], 0.3, {"save_zarr": True}, 9.17),
    "zarrOutputSingle": (_Z[500], 0.3, {"save_zarr": True, "save_16bit": False}, 7.91),
    "zarrInput500": (_Z[500], 0.3, {}, 2.57),
    "partialFiles500": (_Z[500], 0.3, {"partial_file": True}, 2.56),
    "partialFiles2000": (_Z[2000], 0.3, {"partial_file": True}, 8.92),
    "inputBbox": ((500, 900, 512), 0.3, {"read_size": _Z[500]}, 1.33),
    "flatField": (_Z[500], 0.3, {"preprocess": True}, 2.58),
    "flip": (_Z[500], 0.3, {"preprocess": True}, 2.55),
    "outputAxisOrder": (_Z[500], 0.3, {"output_axis_order": "xzy"}, 3.35),
    "resample222": (_Z[500], 0.3, {"resample_factor": (2, 2, 2)}, 2.82),
    "resample1115": (_Z[500], 0.3, {"resample_factor": (1.5, 1, 1)}, 3.71),
    "objectiveScan": (_Z[500], 0.3, {"objective_scan": True}, 6.62),
    "cubic100": (_Z[100], 0.3, {"interp_method": "cubic"}, 1.20),
    "cubic500": (_Z[500], 0.3, {"interp_method": "cubic"}, 4.27),
    "separate500": (_Z[500], 0.3, {"dsr_combined": False}, 9.80),
    "separate1000": (_Z[1000], 0.3, {"dsr_combined": False}, 29.57),
    "deskewOnly500": (_Z[500], 0.3, {"dsr_combined": False, "rotate": False}, 9.33),
    "deskewOnly2000": (_Z[2000], 0.3, {"dsr_combined": False, "rotate": False}, 107.27),
    "rescaleSingle": (_Z[500], 0.3, {"rescale_rotate": True, "save_16bit": False}, 7.60),
}


@pytest.mark.parametrize("name", sorted(MEASURED))
def test_estimate_covers_the_measured_peak(name):
    sz, dz, kw, peak = MEASURED[name]
    mem, _, _ = dm.estimate_deskew_rotate_memory(sz, PX, dz, **kw)
    assert mem >= peak
    assert mem <= peak * 1.3 + 2.5  # MATLAB's test: at most 1.3 x the peak + 2.5 GB


def test_partial_files_are_not_doubled():
    sz = (500, 1800, 512)
    assert dm.estimate_deskew_rotate_memory(sz, PX, 0.3, partial_file=True) == (
        dm.estimate_deskew_rotate_memory(sz, PX, 0.3)
    )
