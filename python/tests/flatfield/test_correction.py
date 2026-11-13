"""Flat-field correction of light-sheet frames (petakit5d.flatfield).

Oracle cases (MATLAB: XR_LSFlatFieldCorrection and processFFCorrectionFrame): the inputs of
tests/flatfield/FlatFieldCorrectionTest.m, then the options (LowerLimit, constOffset,
removeFFImBackground, LSRescale, castDataType) over flat fields and backgrounds of several
classes, 2-D images and z stacks (MATLAB's mean over z), center crops of larger images (even
and odd margins), MATLAB's errors for smaller images and extra dimensions, NaN and Inf in the
flat field, a zero or NaN maximum, every data class (the mex for uint16 and single, MATLAB's
code for the others), 2-D, row, 4-D and larger data, the MATLAB-name layer, and the tiff files
of processFFCorrectionFrame with 4 to 8 arguments. Inputs are MATLAB's arrays; the roles give
the Python calls (z, y, x) volumes. Then FlatFieldCorrectionTest's methods with their exact
expected values, and the reuse of output buffers (``out=``).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.core.arrays import to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.flatfield import correction as fc
from petakit5d.kernels import flatfield as kff
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, File, In, InPath, Value, Vol
from petakit5d.testing.synth import CounterRNG

DTYPES = {
    "logical": np.bool_,
    "uint8": np.uint8,
    "int8": np.int8,
    "uint16": np.uint16,
    "int16": np.int16,
    "uint32": np.uint32,
    "int32": np.int32,
    "uint64": np.uint64,
    "int64": np.int64,
    "single": np.float32,
    "double": np.float64,
}

# ---------------------------------------------------------------------------------------------
# input generators in Python's form ((z, y, x) data, (y, x) images); cases hand MATLAB to_matlab
# ---------------------------------------------------------------------------------------------


def _rng(name, shape):
    return CounterRNG(11, f"flatfield_d01/{name}/{'x'.join(map(str, shape))}")


def raw(shape, cls, name="raw", hi=4000):
    """data (z, y, x) or (y, x): integers 0..hi (some at 0 and at hi), for single and double the
    same plus a fraction of 1/8"""
    rng = _rng(f"{name}/{cls}", shape)
    v = rng.integers(0, hi, shape, dtype=np.int64)
    if v.size > 4:
        f = v.reshape(-1)
        f[::17] = 0
        f[5::23] = hi
    dt = np.dtype(DTYPES[cls])
    if cls == "logical":
        return v % 2 == 1
    if dt.kind == "f":
        frac = rng.spawn("frac").integers(0, 7, shape, dtype=np.int64)
        return (v.astype(dt) + frac.astype(dt) / dt.type(8)).astype(dt)
    return v.astype(dt)


def ls_image(shape, cls="single", name="ls"):
    """a flat-field image (or stack) above a background of about 100: 150..3000, single and
    double with a fraction of 1/4"""
    rng = _rng(f"{name}/{cls}", shape)
    v = rng.integers(150, 3000, shape, dtype=np.int64)
    dt = np.dtype(DTYPES[cls])
    if dt.kind == "f":
        q = rng.spawn("q").integers(0, 3, shape, dtype=np.int64)
        return (v.astype(dt) + q.astype(dt) / dt.type(4)).astype(dt)
    if cls == "uint8":
        return (v % 250 + 5).astype(dt)
    return v.astype(dt)


def bg_image(shape, cls="single", name="bg"):
    """a background image: 95 + 0..15 (single and double in quarters)"""
    rng = _rng(f"{name}/{cls}", shape)
    dt = np.dtype(DTYPES[cls])
    if dt.kind == "f":
        q = rng.integers(0, 60, shape, dtype=np.int64)
        return (dt.type(95) + q.astype(dt) / dt.type(4)).astype(dt)
    return (95 + rng.integers(0, 15, shape, dtype=np.int64)).astype(dt)


def border(img, pads, value_fn):
    """img with a border of pads ((top, bottom), (left, right)) filled by value_fn(shape)"""
    out = value_fn((img.shape[0] + sum(pads[0]), img.shape[1] + sum(pads[1]))).astype(img.dtype)
    out[pads[0][0] : pads[0][0] + img.shape[0], pads[1][0] : pads[1][0] + img.shape[1]] = img
    return out


def border_stack(st, pads, value):
    nz, ny, nx = st.shape
    out = np.full((nz, ny + sum(pads[0]), nx + sum(pads[1])), value, st.dtype)
    out[:, pads[0][0] : pads[0][0] + ny, pads[1][0] : pads[1][0] + nx] = st
    return out


def exact_inputs():
    """FlatFieldCorrectionTest's exactInputs in Python's form"""
    bg = np.full((16, 24), 100, np.float32)
    ff = (bg + np.repeat(np.array([1000, 500, 200], np.float32), 8)[None, :]).astype(np.float32)
    r = np.full((4, 16, 24), 300, np.uint16)
    r[1, 2, :] = 50  # MATLAB raw(3, :, 2)
    return r, ff, bg


def expected_output(r, offset):
    """FlatFieldCorrectionTest's expectedOutput"""
    ff_norm = np.repeat(np.array([1.0, 0.5, 0.4]), 8)[None, :]
    return np.asarray(
        dtypes.cast(np.maximum(r.astype(np.float64) - 100, 0) / ff_norm + offset, "uint16")
    )


# ---------------------------------------------------------------------------------------------
# the Python sides
# ---------------------------------------------------------------------------------------------

_KW = {
    "LowerLimit": "lower_limit",
    "constOffset": "const_offset",
    "removeFFImBackground": "remove_ff_im_background",
    "LSRescale": "ls_rescale",
    "castDataType": "cast_data_type",
}


def xr(rawdata, ls, bg, *nv):
    """XR_LSFlatFieldCorrection(raw, LS, bg, name, value, ...) through the native function"""
    kw = {_KW[nv[i]]: nv[i + 1] for i in range(0, len(nv), 2)}
    return fc.flatfield_correct_volume(rawdata, ls, bg, **kw)


def xr_layer(*args):
    """the MATLAB-name layer (MATLAB's arrays in and out)"""
    from petakit5d import matlab

    return matlab.XR_LSFlatFieldCorrection(*args)


def ff_layer(*args):
    from petakit5d import matlab

    return matlab.processFFCorrectionFrame(*args)


def m(make):
    """a Value input in MATLAB's form of a Python-form generator"""
    return Value(lambda: to_matlab(make()))


def xr_case(cid, r, ls, bg, opts=(), **kw):
    return Case(
        cid,
        matlab="XR_LSFlatFieldCorrection",
        args=[Vol(In("raw")), Vol(In("ls")), Vol(In("bg")), *opts],
        python=xr,
        out_roles=("vol",),
        inputs={"raw": m(r), "ls": m(ls), "bg": m(bg)},
        **kw,
    )


def err_case(cid, r, ls, bg, ident, opts=()):
    return xr_case(
        cid, r, ls, bg, opts, expect_error=ident, py_error=fc.MatlabError, error_message=True
    )


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

LL = ["LowerLimit", 0.4]


def _test_cases():
    """FlatFieldCorrectionTest's calls"""

    def r():
        return exact_inputs()[0]

    def ff():
        return exact_inputs()[1]

    def bg():
        return exact_inputs()[2]

    cases = [
        xr_case("test.formula", r, ff, bg, LL),
        xr_case("test.const_offset", r, ff, bg, [*LL, "constOffset", 5.0]),
        xr_case("test.single", lambda: r().astype(np.float32), ff, bg, LL),
        xr_case("test.nocast", r, ff, bg, [*LL, "castDataType", False]),
        xr_case("test.double_fallback", lambda: r().astype(np.float64), ff, bg, LL),
        xr_case(
            "test.stack",
            r,
            lambda: np.stack([ff() - np.float32(50), ff() + np.float32(50)]).astype(np.float32),
            bg,
            LL,
        ),
        xr_case(
            "test.center_crop",
            r,
            lambda: np.pad(ff(), 2, constant_values=5000),
            lambda: np.pad(bg(), 2, constant_values=0),
            LL,
        ),
        xr_case(
            "test.uint16_images",
            r,
            lambda: np.where(np.arange(24) < 12, 1, 2).astype(np.uint16)[None, :].repeat(16, 0),
            lambda: bg().astype(np.uint16),
            [*LL, "removeFFImBackground", False, "LSRescale", False],
        ),
    ]

    def write_ff(path):
        synth.write_tiff(path, ff())

    def write_bg(path):
        synth.write_tiff(path, bg())

    cases.append(
        Case(
            "test.process_frame",
            matlab="processFFCorrectionFrame",
            args=[Vol(In("raw")), InPath("ff"), InPath("bg"), None, 0.4],
            python=fc.flatfield_correct_file,
            out_roles=("vol",),
            inputs={"raw": m(r), "ff": File(write_ff, "ff.tif"), "bg": File(write_bg, "bg.tif")},
        )
    )
    return cases


SHAPE = (5, 37, 29)  # (z, y, x): MATLAB 37 x 29 x 5
PLANE = SHAPE[1:]


def _option_cases():
    """uint16 data, single images: every combination of the options"""
    cases = []
    for rm in (True, False):
        for rs in (True, False):
            for lname, lval in (
                ("ll04", 0.4),
                ("ll005", 0.05),
                ("llinf", -math.inf),
                ("ll25", 2.5),
            ):
                for cast in (True, False):
                    for oname, off in (("nooff", []), ("off", ["constOffset", 7.5])):
                        cid = (
                            f"opt.{'rm' if rm else 'norm'}.{'rs' if rs else 'nors'}.{lname}."
                            f"{'cast' if cast else 'nocast'}.{oname}"
                        )
                        cases.append(
                            xr_case(
                                cid,
                                lambda: raw(SHAPE, "uint16"),
                                lambda: ls_image(PLANE),
                                lambda: bg_image(PLANE),
                                [
                                    "LowerLimit",
                                    lval,
                                    "removeFFImBackground",
                                    rm,
                                    "LSRescale",
                                    rs,
                                    "castDataType",
                                    cast,
                                    *off,
                                ],
                            )
                        )
    cases.append(
        xr_case(
            "opt.empty_offset",
            lambda: raw(SHAPE, "uint16"),
            lambda: ls_image(PLANE),
            lambda: bg_image(PLANE),
            ["constOffset", None],
        )
    )
    cases.append(
        xr_case(
            "opt.defaults",
            lambda: raw(SHAPE, "uint16"),
            lambda: ls_image(PLANE),
            lambda: bg_image(PLANE),
        )
    )
    cases.append(
        xr_case(
            "opt.single_lower_limit",
            lambda: raw(SHAPE, "uint16"),
            lambda: ls_image(PLANE),
            lambda: bg_image(PLANE),
            ["LowerLimit", np.float32(0.4)],
        )
    )
    cases.append(
        xr_case(
            "opt.negative_offset",
            lambda: raw(SHAPE, "uint16"),
            lambda: ls_image(PLANE),
            lambda: bg_image(PLANE),
            ["constOffset", -40.0],
        )  # uint16 wraps below 0
    )
    return cases


FLAG_SETS = [
    ("rm_rs", True, True),
    ("norm_nors", False, False),
    ("rm_nors", True, False),
    ("norm_rs", False, True),
]


def _class_cases():
    """flat-field images and stacks, and backgrounds, of other classes"""
    cases = []
    for kind, shape in (("image", PLANE), ("stack", (3, *PLANE))):
        for cls in ("uint16", "single", "double", "uint8"):
            for fname, rm, rs in FLAG_SETS:
                cases.append(
                    xr_case(
                        f"cls.ls_{kind}_{cls}.{fname}",
                        lambda: raw(SHAPE, "uint16"),
                        lambda s=shape, c=cls: ls_image(s, c),
                        lambda: bg_image(PLANE),
                        ["LowerLimit", 0.4, "removeFFImBackground", rm, "LSRescale", rs],
                    )
                )
    for cls in ("uint16", "double"):
        for fname, rm, rs in FLAG_SETS[:2]:
            cases.append(
                xr_case(
                    f"cls.bg_{cls}.{fname}",
                    lambda: raw(SHAPE, "uint16"),
                    lambda: ls_image(PLANE),
                    lambda c=cls: bg_image(PLANE, c),
                    ["LowerLimit", 0.4, "removeFFImBackground", rm, "LSRescale", rs],
                )
            )
    # a stack of 37 planes and a long one: MATLAB's sums along z
    cases.append(
        xr_case(
            "cls.ls_stack37_single",
            lambda: raw(SHAPE, "uint16"),
            lambda: ls_image((37, *PLANE), "single", "deep"),
            lambda: bg_image(PLANE),
        )
    )
    cases.append(
        xr_case(
            "cls.ls_stack1500_single",
            lambda: raw((2, 9, 7), "uint16"),
            lambda: ls_image((1500, 9, 7), "single", "long"),
            lambda: bg_image((9, 7)),
        )
    )
    # the mean's last bits reach a single result: stacks of 24-bit fractions, the flat field as
    # it is (sums rounded in single, then single(double(sum) / n))
    for nz, ny, nx in ((37, *PLANE), (1500, 9, 7)):
        cases.append(
            xr_case(
                f"cls.ls_stack{nz}_fine_nocast",
                lambda s=(nz, ny, nx): raw((2, *s[1:]), "single", "fine"),
                lambda s=(nz, ny, nx): fine_stack(s),
                lambda s=(ny, nx): bg_image(s),
                ["removeFFImBackground", False, "LSRescale", False, "castDataType", False],
            )
        )
    return cases


def fine_stack(shape):
    """150 + 2850 * uniform [0, 1) in single (24-bit fractions)"""
    u = _rng("fine", shape).uniform(shape, np.float32)
    return (np.float32(150) + u * np.float32(2850)).astype(np.float32)


def _raw_cases():
    """every class of data (the mex for uint16 and single, MATLAB's code for the others)"""
    cases = []
    for cls in DTYPES:
        hi = {"logical": 1, "uint8": 255, "int8": 127}.get(cls, 4000)
        for cast in (True, False):
            cases.append(
                xr_case(
                    f"raw.{cls}.{'cast' if cast else 'nocast'}",
                    lambda c=cls, h=hi: raw(SHAPE, c, "rawcls", hi=h),
                    lambda: ls_image(PLANE),
                    lambda: bg_image(PLANE),
                    ["castDataType", cast],
                )
            )
    cases.append(
        xr_case(
            "raw.int16.negative",
            lambda: (raw(SHAPE, "int16", "neg") - np.int16(3000)).astype(np.int16),
            lambda: ls_image(PLANE),
            lambda: bg_image(PLANE),
        )
    )
    return cases


def _crop_cases():
    """larger images: center crops with even and odd margins"""
    cases = []

    def big(v):
        return lambda s: np.full(s, v)

    for name, pads in (
        ("y1", ((0, 1), (0, 0))),  # D = ceil(1/2) = 1: the extra row is the first
        ("x1", ((0, 0), (1, 0))),
        ("y2x3", ((1, 1), (2, 1))),
        ("y5x4", ((3, 2), (2, 2))),
        ("y3x3", ((2, 1), (2, 1))),
    ):
        cases.append(
            xr_case(
                f"crop.ls_{name}",
                lambda: raw(SHAPE, "uint16"),
                lambda p=pads: border(ls_image(PLANE), _ceil_pads(p), big(9000)),
                lambda: bg_image(PLANE),
            )
        )
    cases.append(
        xr_case(
            "crop.bg_y1x2",
            lambda: raw(SHAPE, "uint16"),
            lambda: ls_image(PLANE),
            lambda: border(bg_image(PLANE), _ceil_pads(((0, 1), (1, 1))), big(0)),
        )
    )
    cases.append(
        xr_case(
            "crop.both_y4x1",
            lambda: raw(SHAPE, "uint16"),
            lambda: border(ls_image(PLANE), _ceil_pads(((2, 2), (0, 1))), big(9000)),
            lambda: border(bg_image(PLANE), _ceil_pads(((1, 1), (1, 0))), big(0)),
        )
    )
    cases.append(
        xr_case(
            "crop.stack_y3x1",
            lambda: raw(SHAPE, "uint16"),
            lambda: border_stack(
                ls_image((3, *PLANE), "single", "cst"),
                _ceil_pads(((2, 1), (0, 1))),
                np.float32(9000),
            ),
            lambda: bg_image(PLANE),
        )
    )
    return cases


def _ceil_pads(p):
    """pads with MATLAB's center: the crop starts at ceil(margin / 2) (the extra row of an odd
    margin is on the top / left)"""
    (t, b), (lft, r) = p
    my, mx = t + b, lft + r
    return ((-(-my // 2), my - (-(-my // 2))), (-(-mx // 2), mx - (-(-mx // 2))))


def _error_cases():
    """smaller images and extra dimensions: MATLAB's errors (identifiers from MATLAB)"""

    def r():
        return raw(SHAPE, "uint16")

    ny, nx = PLANE
    pos = bounds = fc.INDEX_ERROR  # MATLAB has one identifier for both (the messages differ)
    return [
        err_case(
            "err.ls_y_minus1", r, lambda: ls_image((ny - 1, nx)), lambda: bg_image(PLANE), bounds
        ),
        err_case(
            "err.ls_y_minus2", r, lambda: ls_image((ny - 2, nx)), lambda: bg_image(PLANE), pos
        ),
        err_case(
            "err.ls_x_minus1", r, lambda: ls_image((ny, nx - 1)), lambda: bg_image(PLANE), bounds
        ),
        err_case(
            "err.ls_x_minus3", r, lambda: ls_image((ny, nx - 3)), lambda: bg_image(PLANE), pos
        ),
        err_case(
            "err.ls_y_plus2_x_minus1",
            r,
            lambda: ls_image((ny + 2, nx - 1)),
            lambda: bg_image(PLANE),
            bounds,
        ),
        err_case(
            "err.ls_y_minus1_x_minus2",
            r,
            lambda: ls_image((ny - 1, nx - 2)),
            lambda: bg_image(PLANE),
            bounds,
        ),
        err_case(
            "err.bg_y_minus2", r, lambda: ls_image(PLANE), lambda: bg_image((ny - 2, nx)), pos
        ),
        err_case(
            "err.bg_3d", r, lambda: ls_image(PLANE), lambda: bg_image((2, *PLANE)), fc.SIZE_ERROR
        ),
        err_case(
            "err.ls_4d", r, lambda: ls_image((2, 3, *PLANE)), lambda: bg_image(PLANE), fc.SIZE_ERROR
        ),
    ]


def _special_cases():
    """NaN and Inf in the flat field, a maximum of zero or NaN"""

    def r():
        return raw(SHAPE, "uint16")

    def nan_ls():
        a = ls_image(PLANE, "single", "nan")
        a[3, :] = NAN
        a[5, 7] = np.inf
        a[6, 2] = -np.inf
        return a

    def nan_only_ls():
        a = ls_image(PLANE, "single", "nan")
        a[3, :] = NAN
        a[0, 0] = 1e30
        return a

    def bg_flat():
        return np.full(PLANE, 100, np.float32)

    def zero_max_ls():  # LS - BG <= 0, its maximum +0
        a = np.full(PLANE, 100, np.float32)
        a[::3, :] = 60
        return a

    def neg_max_ls():  # LS - BG < 0 everywhere
        return (
            np.full(PLANE, 40, np.float32) + ls_image(PLANE, "single", "neg") / np.float32(1000)
        ).astype(np.float32)

    def all_nan_ls():
        return np.full(PLANE, NAN, np.float32)

    return [
        xr_case("special.nan_inf", r, nan_ls, lambda: bg_image(PLANE)),
        xr_case(
            "special.nan_inf_llinf", r, nan_ls, lambda: bg_image(PLANE), ["LowerLimit", -math.inf]
        ),
        xr_case("special.nan_big", r, nan_only_ls, lambda: bg_image(PLANE)),
        xr_case("special.zero_max", r, zero_max_ls, bg_flat),
        xr_case("special.negative_max", r, neg_max_ls, bg_flat),
        xr_case("special.negative_max_llinf", r, neg_max_ls, bg_flat, ["LowerLimit", -math.inf]),
        xr_case("special.all_nan", r, all_nan_ls, bg_flat),
        xr_case(
            "special.nan_bg",
            r,
            lambda: ls_image(PLANE),
            lambda: np.where(np.arange(PLANE[1]) == 4, NAN, bg_image(PLANE)).astype(np.float32),
        ),
    ]


def _shape_cases():
    return [
        xr_case(
            "shape.image",
            lambda: raw(PLANE, "uint16", "img"),
            lambda: ls_image(PLANE),
            lambda: bg_image(PLANE),
        ),
        xr_case(
            "shape.row",
            lambda: raw((1, 40), "uint16", "row"),
            lambda: ls_image((3, 44)),
            lambda: bg_image((1, 41)),
        ),
        xr_case(
            "shape.column",
            lambda: raw((40, 1), "single", "col"),
            lambda: ls_image((41, 1)),
            lambda: bg_image((40, 3)),
        ),
        xr_case(
            "shape.4d",
            lambda: raw((2, 3, 14, 11), "uint16", "4d"),
            lambda: ls_image((15, 13)),
            lambda: bg_image((14, 11)),
        ),
        xr_case(
            "shape.big",
            lambda: raw((7, 301, 257), "uint16", "big"),
            lambda: ls_image((3, 311, 263), "single", "big"),
            lambda: bg_image((303, 258)),
            ["constOffset", 10.0],
        ),
        xr_case(
            "shape.big_single",
            lambda: raw((7, 301, 257), "single", "big"),
            lambda: ls_image((311, 263), "uint16", "big"),
            lambda: bg_image((303, 258)),
            ["castDataType", False],
        ),
    ]


def _layer_cases():
    """the MATLAB-name layer: MATLAB's arrays and name-value pairs in, MATLAB's array out"""
    r, ls, bg = (
        (lambda: raw(SHAPE, "uint16")),
        (lambda: ls_image((3, 39, 30))),
        (lambda: bg_image(PLANE)),
    )
    ins = {"raw": m(r), "ls": m(ls), "bg": m(bg)}
    return [
        Case(
            "layer.xr_defaults",
            matlab="XR_LSFlatFieldCorrection",
            args=[In("raw"), In("ls"), In("bg")],
            python=xr_layer,
            inputs=ins,
        ),
        Case(
            "layer.xr_options",
            matlab="XR_LSFlatFieldCorrection",
            args=[
                In("raw"),
                In("ls"),
                In("bg"),
                "lowerlimit",
                0.2,
                "CONSTOFFSET",
                3.0,
                "LSRescale",
                False,
                "castDataType",
                False,
            ],
            python=xr_layer,
            inputs=ins,
        ),
        # integer-valued doubles (int on the Python side of the layer), MATLAB's code for uint8
        Case(
            "layer.xr_int_values_uint8",
            matlab="XR_LSFlatFieldCorrection",
            args=[In("raw8"), In("ls"), In("bg"), "LowerLimit", 1.0, "constOffset", 5.0],
            python=xr_layer,
            inputs={**ins, "raw8": m(lambda: raw(SHAPE, "uint8", "l8", hi=255))},
        ),
    ]


def _file_cases():
    """processFFCorrectionFrame with tiff files and 4 to 8 arguments"""

    def r():
        return raw(SHAPE, "uint16")

    def w(make):
        return lambda path: synth.write_tiff(path, make())

    files = {
        "raw": m(r),
        "ls": File(w(lambda: ls_image((39, 30))), "ls.tif"),
        "bg": File(w(lambda: bg_image((38, 29))), "bg.tif"),
        "lsst": File(w(lambda: ls_image((4, 37, 31), "uint16", "st")), "ls_stack.tif"),
        "bgu": File(w(lambda: bg_image(PLANE, "uint16")), "bg_u16.tif"),
    }
    out = []
    for cid, ls_key, bg_key, rest in (
        ("file.args4", "ls", "bg", [None]),
        ("file.args4_offset", "ls", "bg", [5.0]),
        ("file.args5", "ls", "bg", [None, 0.1]),
        ("file.args6", "lsst", "bgu", [None, 0.4, False]),
        ("file.args7", "lsst", "bgu", [2.0, 0.4, True, False]),
        ("file.args8", "lsst", "bgu", [None, 0.3, True, True, False]),
    ):
        ins = {"raw": files["raw"], ls_key: files[ls_key], bg_key: files[bg_key]}
        out.append(
            Case(
                cid,
                matlab="processFFCorrectionFrame",
                args=[Vol(In("raw")), InPath(ls_key), InPath(bg_key), *rest],
                python=fc.flatfield_correct_file,
                out_roles=("vol",),
                inputs=ins,
            )
        )
    ins = {"raw": files["raw"], "ls": files["ls"], "bg": files["bg"]}
    out.append(
        Case(
            "file.layer_args8",
            matlab="processFFCorrectionFrame",
            args=[In("raw"), InPath("ls"), InPath("bg"), None, 0.4, True, True, True],
            python=ff_layer,
            inputs=ins,
        )
    )
    return out


CASES = (
    _test_cases()
    + _option_cases()
    + _class_cases()
    + _raw_cases()
    + _crop_cases()
    + _error_cases()
    + _special_cases()
    + _shape_cases()
    + _layer_cases()
    + _file_cases()
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# FlatFieldCorrectionTest
# ---------------------------------------------------------------------------------------------


def test_corrects_with_expected_formula():
    r, ff, bg = exact_inputs()
    out = fc.flatfield_correct_volume(r, ff, bg, lower_limit=0.4)
    assert out.dtype == np.uint16
    np.testing.assert_array_equal(out, expected_output(r, 100))


def test_const_offset_replaces_background():
    r, ff, bg = exact_inputs()
    out = fc.flatfield_correct_volume(r, ff, bg, lower_limit=0.4, const_offset=5)
    np.testing.assert_array_equal(out, expected_output(r, 5))


def test_single_input_gives_single():
    r, ff, bg = exact_inputs()
    exp = expected_output(r, 100).astype(np.float32)
    out = fc.flatfield_correct_volume(r.astype(np.float32), ff, bg, lower_limit=0.4)
    assert out.dtype == np.float32
    np.testing.assert_allclose(out, exp, atol=1e-3, rtol=0)
    out = fc.flatfield_correct_volume(r, ff, bg, lower_limit=0.4, cast_data_type=False)
    assert out.dtype == np.float32
    np.testing.assert_allclose(out, exp, atol=1e-3, rtol=0)


def test_matlab_fallback_matches_mex():
    r, ff, bg = exact_inputs()
    out = fc.flatfield_correct_volume(r.astype(np.float64), ff, bg, lower_limit=0.4)
    assert out.dtype == np.float64
    np.testing.assert_allclose(out, expected_output(r, 100).astype(np.float64), atol=1e-3)


def test_flat_field_stack_is_averaged_over_z():
    r, ff, bg = exact_inputs()
    ff3 = np.stack([ff - np.float32(50), ff + np.float32(50)])
    out = fc.flatfield_correct_volume(r, ff3, bg, lower_limit=0.4)
    np.testing.assert_array_equal(out, expected_output(r, 100))


def test_larger_images_are_center_cropped():
    r, ff, bg = exact_inputs()
    out = fc.flatfield_correct_volume(
        r, np.pad(ff, 2, constant_values=5000), np.pad(bg, 2, constant_values=0), lower_limit=0.4
    )
    np.testing.assert_array_equal(out, expected_output(r, 100))


def test_non_single_flat_field_uses_fast_path(monkeypatch):
    if not kff.extension_available():
        pytest.skip("the flat-field extension is not built")
    r, _, bg = exact_inputs()
    ff = np.ones(bg.shape, np.uint16)
    ff[:, 12:] = 2

    def no_fallback(*a, **k):
        raise AssertionError("fell back to MATLAB's code")

    monkeypatch.setattr(kff, "apply_flat_field", no_fallback)
    out = fc.flatfield_correct_volume(
        r,
        ff,
        bg.astype(np.uint16),
        lower_limit=0.4,
        remove_ff_im_background=False,
        ls_rescale=False,
    )
    exp = dtypes.cast(np.maximum(r.astype(np.float64) - 100, 0) / ff + 100, "uint16")
    np.testing.assert_array_equal(out, exp)


def test_process_frame_reads_images_from_disk(tmp_path):
    r, ff, bg = exact_inputs()
    synth.write_tiff(tmp_path / "ff.tif", ff)
    synth.write_tiff(tmp_path / "bg.tif", bg)
    out = fc.flatfield_correct_file(
        r, str(tmp_path / "ff.tif"), str(tmp_path / "bg.tif"), None, 0.4
    )
    np.testing.assert_array_equal(out, expected_output(r, 100))


# ---------------------------------------------------------------------------------------------
# the preparation, the crop, out= and threads
# ---------------------------------------------------------------------------------------------


def test_crop_starts_at_ceil_of_half_the_margin():
    img = np.arange(7 * 9, dtype=np.float32).reshape(7, 9)
    c = fc._crop(img, (4, 6))  # margins 3 and 3: D = 2, rows/columns 3.. (1-based)
    np.testing.assert_array_equal(c, img[2:6, 2:8])
    c = fc._crop(img, (6, 9))  # margin 1: D = 1
    np.testing.assert_array_equal(c, img[1:7, :])
    with pytest.raises(fc.MatlabError) as e:
        fc._crop(img, (8, 9))  # D = 0, 8 > 7
    assert e.value.identifier == fc.INDEX_ERROR
    with pytest.raises(fc.MatlabError) as e:
        fc._crop(img, (7, 11))  # D = -1
    assert e.value.identifier == fc.INDEX_ERROR
    with pytest.raises(fc.MatlabError) as e:
        fc._crop(np.zeros((2, 7, 9), np.float32), (7, 9))
    assert e.value.identifier == fc.SIZE_ERROR


def test_prepare_keeps_matlab_classes():
    r, ff, bg = exact_inputs()
    f, b = fc.prepare_flat_field(ff, bg, r.shape[1:])
    assert f.dtype == np.float32 and b.dtype == np.float32
    np.testing.assert_array_equal(f, np.repeat(np.float32([1, 0.5, 0.4]), 8)[None, :].repeat(16, 0))
    u = ff.astype(np.uint16)
    f, _ = fc.prepare_flat_field(
        u, bg, r.shape[1:], remove_ff_im_background=False, ls_rescale=False
    )
    assert f.dtype == np.uint16  # max(uint16, 0.4) stays uint16
    st = np.stack([u, u + np.uint16(1)])
    f, _ = fc.prepare_flat_field(
        st, bg, r.shape[1:], remove_ff_im_background=False, ls_rescale=False
    )
    assert f.dtype == np.float64  # mean of uint16 is double
    np.testing.assert_array_equal(f, u + 0.5)


def test_nan_in_flat_field_becomes_the_lower_limit():
    r, ff, bg = exact_inputs()
    ff = ff.copy()
    ff[0, 0] = NAN
    f, _ = fc.prepare_flat_field(ff, bg, r.shape[1:], lower_limit=0.25)
    assert f[0, 0] == np.float32(0.25)
    assert f[0, 1] == np.float32(1)  # the maximum ignores NaN


@pytest.mark.parametrize("cls,cast", [("uint16", True), ("uint16", False), ("single", True)])
def test_out_is_reused_frame_after_frame(cls, cast):
    shape = (6, 33, 41)
    ls, bg = ls_image(shape[1:]), bg_image(shape[1:])
    first = raw(shape, cls, "f1")
    ref = fc.flatfield_correct_volume(first, ls, bg, cast_data_type=cast)
    buf = np.empty_like(ref)
    for k in range(3):
        fr = raw(shape, cls, f"f{k + 1}")
        exp = fc.flatfield_correct_volume(fr, ls, bg, cast_data_type=cast)
        got = fc.flatfield_correct_volume(fr, ls, bg, cast_data_type=cast, out=buf)
        assert got is buf
        assert got.tobytes() == exp.tobytes()


def test_out_on_the_fallback_path():
    shape = (3, 12, 10)
    fr = raw(shape, "double", "fb")
    ls, bg = ls_image(shape[1:]), bg_image(shape[1:])
    exp = fc.flatfield_correct_volume(fr, ls, bg)
    buf = np.empty_like(exp)
    assert fc.flatfield_correct_volume(fr, ls, bg, out=buf) is buf
    assert buf.tobytes() == exp.tobytes()
    with pytest.raises(ValueError):
        fc.flatfield_correct_volume(fr, ls, bg, out=np.empty(shape, np.float32))


def test_out_must_match_the_result():
    shape = (3, 12, 10)
    fr = raw(shape, "uint16", "om")
    ls, bg = ls_image(shape[1:]), bg_image(shape[1:])
    with pytest.raises(ValueError):
        fc.flatfield_correct_volume(fr, ls, bg, out=np.empty(shape, np.float32))
    with pytest.raises(ValueError):
        fc.flatfield_correct_volume(fr, ls, bg, out=np.empty((3, 12, 11), np.uint16))


def test_file_out_and_threads(tmp_path):
    shape = (9, 70, 101)
    fr = raw(shape, "uint16", "thr")
    synth.write_tiff(tmp_path / "ls.tif", ls_image((2, 72, 101)))
    synth.write_tiff(tmp_path / "bg.tif", bg_image(shape[1:]))
    a = [str(tmp_path / "ls.tif"), str(tmp_path / "bg.tif"), 3.0]
    ref = fc.flatfield_correct_file(fr, *a, num_threads=1)
    buf = np.empty_like(ref)
    for t in (2, 3, 8):
        assert fc.flatfield_correct_file(fr, *a, out=buf, num_threads=t) is buf
        assert buf.tobytes() == ref.tobytes()


def test_matlab_layer_equals_native():
    from petakit5d import matlab

    shape = (4, 21, 17)
    fr, ls, bg = raw(shape, "uint16", "lay"), ls_image((2, 23, 17)), bg_image(shape[1:])
    nat = fc.flatfield_correct_volume(fr, ls, bg, lower_limit=0.3, ls_rescale=False)
    lay = matlab.XR_LSFlatFieldCorrection(
        to_matlab(fr), to_matlab(ls), bg, "LowerLimit", 0.3, "LSRescale", False
    )
    assert lay.tobytes(order="F") == to_matlab(nat).tobytes(order="F")


def _special_singles(shape, name):
    a = ls_image(shape, "single", name) - np.float32(1500)
    f = a.reshape(-1)
    f[::7] = NAN
    f[1::11] = np.inf
    f[2::13] = -np.inf
    f[3::17] = -0.0
    f[4::19] = 0.0
    return a


def test_fast_paths_give_the_bits_of_matlab_helpers():
    from petakit5d.matlabcompat import reduce

    a, b = _special_singles((13, 17), "fa"), _special_singles((13, 17), "fb")
    for op, ref, x, y in (
        (np.subtract, dtypes.minus, a, b),
        (np.subtract, dtypes.minus, np.abs(a), np.abs(b)),  # no NaN: NumPy's operation
        (np.divide, dtypes.rdivide, a, np.float32(3.7)),
        (np.divide, dtypes.rdivide, a, np.float32(0)),
    ):
        got = fc._single_op(op, ref, x, y)
        exp = np.asarray(ref(x, y))
        assert got.dtype == exp.dtype == np.float32
        assert got.tobytes() == exp.tobytes()
    for lim in (0.4, 0.0, -0.0, -math.inf, math.inf, NAN, np.float32(0.3), 2, np.float64(-3)):
        got = fc._lower_limit(a, lim)
        exp = np.asarray(reduce.maximum(a, fc._double(lim)))  # a Python int is a double
        assert got.dtype == exp.dtype and got.tobytes() == exp.tobytes(), lim


def test_python_ints_are_doubles():
    # the MATLAB-name layer gives integer-valued doubles as int: LowerLimit 1, constOffset 5
    shape = (3, 12, 10)
    ls, bg = ls_image(shape[1:]), bg_image(shape[1:])
    for cls in ("uint16", "double", "uint8"):
        fr = raw(shape, cls, "ints", hi=200)
        a = fc.flatfield_correct_volume(fr, ls, bg, lower_limit=1, const_offset=5)
        b = fc.flatfield_correct_volume(fr, ls, bg, lower_limit=1.0, const_offset=5.0)
        assert a.dtype == b.dtype and a.tobytes() == b.tobytes()


def test_flags_read_as_matlab_reads_them():
    assert fc._if(True) and not fc._if(False)
    assert not fc._if(np.array([True, False])) and not fc._if(np.zeros((0,), bool))
    assert fc._first(np.array([[True, False]])) and not fc._first(np.zeros((0,), bool))
    assert fc._given(5) and not fc._given(None) and not fc._given(np.zeros((0, 0)))
