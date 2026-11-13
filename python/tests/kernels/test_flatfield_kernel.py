"""Flat-field correction kernel (petakit5d.kernels.flatfield).

Oracle cases (MATLAB: flat_field_correction_3d_mex, and XR_LSFlatFieldCorrection with the flat
field used as it is, which reaches the mex or, for the classes the mex refuses, its MATLAB code)
over uint16 and single data, the cast flag, constant offsets (none, positive, negative,
fractional, NaN, Inf, other classes), flat fields with lower limits, zeros, negative values,
NaN and Inf, uint16 results beyond 65535 and below 0 (wrapped by the mex), -0, 2-D, row, empty
and odd-sized volumes, and the input errors. Their inputs are MATLAB's (y, x, z) arrays; the
roles give the Python calls (z, y, x) volumes. Then tests/platform/MexInputCheckTest.m's
flatFieldCorrectionChecksItsInputs, and the kernel against its NumPy reference at several
thread counts, instruction sets and layouts.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.kernels import flatfield as ff
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Vol
from petakit5d.testing.synth import CounterRNG

MOD = "petakit5d.kernels.flatfield"
INPUT_ERROR = "flat_field_correction_3d_mex:inputError"
TYPE_ERROR = "flat_field_correction_3d_mex:dataTypeError"
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
# input generators: MATLAB's (y, x, z) arrays (the case keys are made from these bytes)
# ---------------------------------------------------------------------------------------------


def _rng(name, shape):
    return CounterRNG(7, f"flatfield/{name}/{'x'.join(map(str, shape))}")


def raw(shape, cls, name="raw", hi=4000):
    """data of MATLAB size shape: integers 0..hi (a few at 0 and at hi), or for single and
    double the same plus a fraction of 1/8"""
    pshape = tuple(shape[2:][::-1]) + tuple(shape[:2]) if len(shape) > 2 else tuple(shape)
    rng = _rng(f"{name}/{cls}", shape)
    v = rng.integers(0, hi, pshape, dtype=np.int64)
    if v.size > 4:
        flat = v.reshape(-1)
        flat[::17] = 0
        flat[5::23] = hi
    dt = np.dtype(DTYPES[cls])
    if cls == "logical":
        return to_matlab(v % 2 == 1)
    if dt.kind == "f":
        frac = rng.spawn("frac").integers(0, 7, pshape, dtype=np.int64)
        return to_matlab((v.astype(dt) + frac.astype(dt) / dt.type(8)).astype(dt))
    return to_matlab(v.astype(dt))


def flat(shape, name="flat", lo=0.3, hi=1.3):
    """a single flat field (y, x): uniform in [lo, hi), a column at the lower limit 0.4"""
    u = _rng(name, shape).uniform(shape, np.float32)
    f = (np.float32(lo) + u * np.float32(hi - lo)).astype(np.float32)
    if f.shape[1] > 1:
        f[:, 1] = np.float32(0.4)
    return f


def background(shape, name="bg", base=100.0):
    """a single background (y, x): base plus 0 to 15 in quarters"""
    q = _rng(name, shape).integers(0, 60, shape, dtype=np.int64)
    return (np.float32(base) + q.astype(np.float32) / np.float32(4)).astype(np.float32)


def special_flat(shape):
    """a flat field with zeros of both signs, negative values, NaN, Inf and tiny values"""
    f = flat(shape, "special")
    f[0, 0] = 0.0
    f[1, 0] = -0.0
    f[2, :3] = -0.5
    f[3, 0] = NAN
    f[3, 2] = np.inf
    f[4, 0] = -np.inf
    f[5, :] = np.float32(1e-6)  # results beyond 2^31
    f[6, :] = np.float32(0.01)  # results beyond 65535
    return f


def special_bg(shape):
    b = background(shape, "special")
    b[7, 0] = NAN
    b[7, 1] = np.inf
    b[7, 2] = -np.inf
    b[8, :] = 0.0
    b[9, :] = -50.0
    return b


def special_raw(shape, cls):
    """uint16 or single data with 0 and 65535; single also NaN, Inf, -Inf and -0"""
    a = raw(shape, cls, "special", hi=65535)
    if cls == "single":
        a[0, 1, :] = NAN
        a[1, 1, :] = np.inf
        a[2, 1, :] = -np.inf
        a[8, :, 0] = -0.0
        a[9, 2, :] = 1e30
    return a


# ---------------------------------------------------------------------------------------------
# oracle cases: flat_field_correction_3d_mex
# ---------------------------------------------------------------------------------------------

SHAPES = {
    "odd": (37, 23, 5),
    "image": (19, 31),
    "row": (1, 40),
    "column": (40, 1),
    "plane1": (13, 11, 1),
    "big": (301, 257, 17),  # several tasks per thread
    "longrow": (3, 9001, 2),  # rows longer than a block
}

# MATLAB's castDataType values: (name, value, read as true)
CASTS = [("true", True), ("false", False), ("one", 1.0), ("zero", 0.0)]

OFFSETS = [
    ("off5", 5.0),
    ("off0", 0.0),
    ("offneg", -3.25),  # uint16 results below 0 wrap from 65535 down
    ("offfrac", 0.1),  # rounded to single
    ("offbig", 1e10),  # beyond int32: 0 in uint16
]


def _mex_case(cid, im, fl, bg, extra, **kw):
    inputs = {"im": Value(im), "ff": Value(fl), "bg": Value(bg)}
    return Case(
        cid,
        matlab="flat_field_correction_3d_mex",
        args=[Vol(In("im")), Vol(In("ff")), Vol(In("bg")), *extra],
        python=f"{MOD}:flat_field_correction_3d_mex",
        out_roles=("vol",),
        inputs=inputs,
        **kw,
    )


def _plane(shape):
    return tuple(shape[:2])


def _mex_cases():
    cases = []
    for cls in ("uint16", "single"):
        for sname, shape in SHAPES.items():
            p = _plane(shape)

            def im(s=shape, c=cls):
                return raw(s, c)

            def fl(p=p):
                return flat(p)

            def bg(p=p):
                return background(p)

            for cname, cast in CASTS:
                if sname not in ("odd", "image") and cname not in ("true", "false"):
                    continue
                cases.append(_mex_case(f"mex.{cls}.{sname}.{cname}", im, fl, bg, [cast]))
            for oname, off in OFFSETS:
                if sname in ("big", "longrow", "column") and oname not in ("off5", "offneg"):
                    continue
                for cname, cast in CASTS[:2]:
                    cases.append(
                        _mex_case(f"mex.{cls}.{sname}.{cname}.{oname}", im, fl, bg, [cast, off])
                    )
        # flat fields and backgrounds with zeros, negatives, NaN and Inf; data with the
        # extremes (and NaN, Inf, -0 in single)
        shape = (12, 9, 4)
        p = _plane(shape)

        def im(s=shape, c=cls):
            return special_raw(s, c)

        for cname, cast in CASTS[:2]:
            for extra, oname in (([], "nooff"), ([7.0], "off7"), ([-1000.0], "offneg")):
                cases.append(
                    _mex_case(
                        f"mex.{cls}.specials.{cname}.{oname}",
                        im,
                        lambda p=p: special_flat(p),
                        lambda p=p: special_bg(p),
                        [cast, *extra],
                    )
                )
    # uint16 wrapping and the int32 range, value by value: f = 1, b = 0, offsets that put the
    # results at the edges
    edge = np.array([[0, 1, 2, 65534, 65535], [10, 20, 30, 40, 50]], dtype=np.uint16).reshape(
        2, 5, 1
    )
    ones = np.ones((2, 5), np.float32)
    zeros = np.zeros((2, 5), np.float32)
    for oname, off in [
        ("half", 0.5),
        ("below_half", 0.49999997),
        ("wrap", 65535.5),
        ("m1", -1.0),
        ("m05", -0.5),
        ("m051", -0.51),
        ("m70000", -70000.0),
        ("int32max", 2147483520.0),
        ("int32over", 2147483648.0),
        ("int32min", -2147483648.0),
        ("int32under", -2147483904.0),
        ("nan", math.nan),
        ("inf", math.inf),
        ("minf", -math.inf),
    ]:
        cases.append(
            _mex_case(
                f"mex.uint16.edge.{oname}",
                lambda: edge,
                lambda: ones,
                lambda: zeros,
                [True, off],
            )
        )
    # -0: (-0 - 0) / f is -0 and stays (-0 >= 0); -0 + -0 is -0
    negz = np.full((3, 4, 2), -0.0, np.float32)
    negz[0, 0, 0] = 2.0
    for oname, extra in (("negzero", [False, -0.0]), ("zero", [False, 0.0]), ("bg", [False])):
        cases.append(
            _mex_case(
                f"mex.single.negzero.{oname}",
                lambda: negz,
                lambda: np.ones((3, 4), np.float32),
                lambda: np.zeros((3, 4), np.float32),
                extra,
            )
        )
    # mxGetScalar's reading of the flag and the offset
    shape = (9, 7, 3)

    def im():
        return raw(shape, "uint16")

    def fl():
        return flat(shape[:2])

    def bg():
        return background(shape[:2])

    for cname, cast in [
        ("flag_2p5", 2.5),
        ("flag_nan", math.nan),
        ("flag_neg", -1.0),
        ("flag_int8_0", np.int8(0)),
        ("flag_uint16_3", np.uint16(3)),
        ("flag_row_0_1", np.array([[0.0, 1.0]])),
        ("flag_row_1_0", np.array([[1.0, 0.0]])),
        ("flag_single_tiny", np.float32(1e-30)),
    ]:
        cases.append(_mex_case(f"mex.uint16.{cname}", im, fl, bg, [cast]))
    for oname, off in [
        ("off_single", np.float32(2.75)),
        ("off_int16", np.int16(-7)),
        ("off_uint8", np.uint8(200)),
        ("off_logical", True),
        ("off_row", np.array([[3.0, 100.0]])),
        ("off_double_rounds", 1.0000000596046448),  # half way between two singles
    ]:
        cases.append(_mex_case(f"mex.uint16.{oname}", im, fl, bg, [True, off]))
        cases.append(
            _mex_case(f"mex.single.{oname}", lambda: raw(shape, "single"), fl, bg, [True, off])
        )
    # empty data
    for cls in ("uint16", "single"):
        cases.append(
            _mex_case(
                f"mex.{cls}.empty",
                lambda c=cls: np.zeros((0, 4, 3), DTYPES[c], order="F"),
                lambda: np.zeros((0, 4), np.float32),
                lambda: np.zeros((0, 4), np.float32),
                [True],
            )
        )
    cases += _error_cases()
    return cases


def _error_cases():
    """the mex's checks (MexInputCheckTest.flatFieldCorrectionChecksItsInputs and the rest)"""
    im = np.ones((4, 3, 2), np.uint16)
    fl = np.ones((4, 3), np.float32)
    bg = np.zeros((4, 3), np.float32)
    cases = []

    def err(cid, i, f, b, extra, ident):
        cases.append(
            _mex_case(
                f"mex.error.{cid}",
                lambda: i,
                lambda: f,
                lambda: b,
                extra,
                expect_error=ident,
                py_error=ff.MatlabError,
            )
        )

    err("4d_data", np.ones((4, 3, 2, 2), np.uint16), fl, bg, [True], INPUT_ERROR)
    err("3d_background", im, fl, np.zeros((4, 3, 2), np.float32), [True], INPUT_ERROR)
    err("3d_flat", im, np.ones((4, 3, 2), np.float32), bg, [True], INPUT_ERROR)
    err("flat_size", im, np.ones((4, 4), np.float32), bg, [True], INPUT_ERROR)
    err("flat_transposed", im, np.ones((3, 4), np.float32), bg, [True], INPUT_ERROR)
    err("background_size", im, fl, np.zeros((5, 3), np.float32), [True], INPUT_ERROR)
    err("flat_double", im, np.ones((4, 3)), bg, [True], INPUT_ERROR)
    err("background_uint16", im, fl, np.zeros((4, 3), np.uint16), [True], INPUT_ERROR)
    err("three_args", im, fl, bg, [], INPUT_ERROR)
    err("six_args", im, fl, bg, [True, 1.0, 2.0], INPUT_ERROR)
    # the sizes are checked before the classes, the image classes before the data's
    err("size_before_class", im.astype(np.float64), np.ones((4, 4)), bg, [True], INPUT_ERROR)
    err("flat_class_before_data", im.astype(np.float64), np.ones((4, 3)), bg, [True], INPUT_ERROR)
    for cls in ("double", "uint8", "int16", "uint32", "logical"):
        err(f"data_{cls}", im.astype(DTYPES[cls]), fl, bg, [True], TYPE_ERROR)
    # accepted: the flat field and background of a 2-D image, and the check test's identity
    cases.append(
        _mex_case(
            "mex.check.identity",
            lambda: im,
            lambda: fl,
            lambda: bg,
            [True],
        )
    )
    return cases


# ---------------------------------------------------------------------------------------------
# oracle cases: XR_LSFlatFieldCorrection with the flat field used as it is
# ---------------------------------------------------------------------------------------------


def xr_correct(rawdata, ls, bg, *opts):
    """XR_LSFlatFieldCorrection(raw, LS, bg, 'LowerLimit', -Inf, 'removeFFImBackground',
    false, 'LSRescale', false, ...) through apply_flat_field: with these options the flat field
    and the background reach the correction as they are (same sizes, no crop)"""
    o = {opts[i]: opts[i + 1] for i in range(0, len(opts), 2)}
    assert o.pop("LowerLimit") == -math.inf
    assert o.pop("removeFFImBackground") is False and o.pop("LSRescale") is False
    if ls.dtype.kind == "f":
        ls = np.where(np.isnan(ls), ls.dtype.type(-np.inf), ls)  # max(LS, -Inf): NaN is -Inf
    return ff.apply_flat_field(
        rawdata,
        ls,
        bg,
        const_offset=o.pop("constOffset", None),
        cast_data_type=o.pop("castDataType", True),
    )


XR_OPTS = ["LowerLimit", -math.inf, "removeFFImBackground", False, "LSRescale", False]


def _xr_case(cid, im, fl, bg, extra, **kw):
    return Case(
        cid,
        matlab="XR_LSFlatFieldCorrection",
        args=[Vol(In("im")), Vol(In("ff")), Vol(In("bg")), *XR_OPTS, *extra],
        python=xr_correct,
        out_roles=("vol",),
        inputs={"im": Value(im), "ff": Value(fl), "bg": Value(bg)},
        **kw,
    )


def _xr_cases():
    cases = []
    shape = (14, 11, 3)
    p = shape[:2]

    def fl():
        return flat(p)

    def bg():
        return background(p)

    for cls in DTYPES:
        hi = {"logical": 1, "uint8": 255, "int8": 127}.get(cls, 4000)

        def im(c=cls, h=hi):
            return raw(shape, c, "xr", hi=h)

        for cast in (True, False):
            for oname, extra in (
                ("nooff", []),
                ("off5", ["constOffset", 5.0]),
                ("offneg", ["constOffset", -2.5]),
            ):
                cases.append(
                    _xr_case(
                        f"xr.{cls}.{'cast' if cast else 'nocast'}.{oname}",
                        im,
                        fl,
                        bg,
                        ["castDataType", cast, *extra],
                    )
                )

    # signed classes: negative data, saturation below; int16 offset (integer result class)
    def neg():
        return (raw(shape, "int16", "xrneg") - np.int16(3000)).astype(np.int16)

    cases.append(_xr_case("xr.int16.negative", neg, fl, bg, ["castDataType", True]))
    cases.append(
        _xr_case(
            "xr.int16.offset_int16",
            neg,
            fl,
            bg,
            ["castDataType", True, "constOffset", np.int16(-30000)],
            # MATLAB's code adds the offset to single data: an error for an integer class
            expect_error="MATLAB:mixedClasses",
            py_error=ff.MatlabError,
        )
    )
    cases.append(
        _xr_case(
            "xr.double.specials",
            lambda: special_raw((12, 9, 4), "single").astype(np.float64),
            lambda: special_flat((12, 9)),
            lambda: special_bg((12, 9)),
            ["castDataType", True],
        )
    )

    # MATLAB's code saturates where the mex wraps
    def big():
        return raw(shape, "uint8", "xrbig", hi=255)

    cases.append(
        _xr_case(
            "xr.uint8.saturates",
            big,
            lambda: (flat(p) * np.float32(0.25)).astype(np.float32),
            bg,
            ["castDataType", True, "constOffset", 100.0],
        )
    )
    # 4-D data: the mex refuses, MATLAB's code broadcasts over z and t (with a fall-back warning)
    cases.append(
        _xr_case(
            "xr.uint16.4d",
            lambda: raw((14, 11, 3, 2), "uint16", "xr4d"),
            fl,
            bg,
            ["castDataType", True],
        )
    )
    # uint16 flat field and background: made single before the mex (the fast path)
    cases.append(
        _xr_case(
            "xr.uint16.uint16_images",
            lambda: raw(shape, "uint16", "xru16"),
            lambda: np.where(np.arange(p[1]) < 5, 1, 2).astype(np.uint16)[None, :].repeat(p[0], 0),
            lambda: np.full(p, 100, np.uint16),
            ["castDataType", True],
        )
    )
    return cases


CASES = _mex_cases() + _xr_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# helpers on (z, y, x) volumes
# ---------------------------------------------------------------------------------------------


def native(shape, cls, name="native", hi=4000):
    """a (z, y, x) volume of a class"""
    m = (shape[1], shape[2], shape[0]) if len(shape) == 3 else tuple(shape)
    return from_matlab(raw(m, cls, name, hi=hi))


def same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, a.shape, b.dtype, b.shape)
    return np.array_equal(
        np.ascontiguousarray(a).view(np.uint8), np.ascontiguousarray(b).view(np.uint8)
    )


def by_hand(v, f, b, cast, offset):
    """the mex's formula voxel by voxel in Python floats rounded to single at each step"""
    out_u16 = v.dtype == np.uint16 and cast
    res = np.empty(v.shape, np.uint16 if out_u16 else np.float32)
    f32 = np.float32
    with np.errstate(all="ignore"):
        for idx in np.ndindex(v.shape):
            yx = idx[-2:]
            w = (f32(v[idx]) - b[yx]) / f[yx]
            w = w if w >= 0 else f32(0)
            w = w + (b[yx] if offset is None else f32(offset))
            if out_u16:
                d = float(w) + 0.5
                t = math.trunc(d) if -2147483649.0 < d < 2147483648.0 else -(2**31)
                res[idx] = t & 0xFFFF
            else:
                res[idx] = w
    return res


# ---------------------------------------------------------------------------------------------
# ported MATLAB test: MexInputCheckTest.flatFieldCorrectionChecksItsInputs
# ---------------------------------------------------------------------------------------------


def test_flat_field_correction_checks_its_inputs():
    im = np.ones((2, 4, 3), np.uint16)
    f = np.ones((4, 3), np.float32)
    b = np.zeros((4, 3), np.float32)
    for args, what in [
        ((np.ones((2, 2, 4, 3), np.uint16), f, b, True), "4-D data"),
        ((im, f, np.zeros((2, 4, 3), np.float32), True), "3-D background"),
        ((im, np.ones((2, 4, 3), np.float32), b, True), "3-D flat field"),
        ((im, f.T.copy(), b, True), "flat field of another size"),
        ((im, f, b), "three arguments"),
        ((im, f.astype(np.float64), b, True), "double flat field"),
    ]:
        with pytest.raises(ff.MatlabError) as e:
            ff.flat_field_correction_3d_mex(*args)
        assert e.value.identifier == INPUT_ERROR, what
    with pytest.raises(ff.MatlabError) as e:
        ff.flat_field_correction_3d_mex(im.astype(np.float64), f, b, True)
    assert e.value.identifier == TYPE_ERROR
    assert same_bits(ff.flat_field_correction_3d_mex(im, f, b, True), im)
    with pytest.raises(ff.MatlabError) as e:
        ff.flat_field_correction_3d(im, f, b[:3])
    assert e.value.identifier == INPUT_ERROR


# ---------------------------------------------------------------------------------------------
# the kernel against its references
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("cls", ["uint16", "single"])
@pytest.mark.parametrize("cast", [True, False])
@pytest.mark.parametrize("offset", [None, 5.0, -3.25, 0.1])
def test_kernel_matches_formula(cls, cast, offset):
    v = native((3, 7, 5), cls, "hand", hi=65535)
    f = flat((7, 5), "hand", lo=0.01, hi=1.3)
    f[0, 0] = 0.0
    f[1, 1] = NAN
    b = background((7, 5), "hand")
    b[2, 2] = -40.0
    want = by_hand(v, f, b, cast, offset)
    got = ff.flat_field_correction_3d(v, f, b, cast_data_type=cast, const_offset=offset)
    assert same_bits(got, want)
    ref = ff.flat_field_correction_3d_ref(v, f, b, cast_data_type=cast, const_offset=offset)
    assert same_bits(ref, want)


def test_uint16_results_wrap_as_the_mex_converts():
    """(uint16)(w + 0.5) through 32 bits (GCC's cvttsd2si into a 32-bit register): the low 16
    bits of the truncation, 0 for NaN and outside int32 (known issue: the mex wraps; MATLAB's
    cast would saturate). With data 0, flat 1 and background 0 the result is the offset's."""
    cases = [
        (-3.2, 65534),
        (70000, 4464),
        (3e9, 0),
        (-3e9, 0),
        (1e20, 0),
        (65535.6, 0),
        (math.nan, 0),
        (-0.7, 0),
        (2147483520.0, 65408),
        (2147483648.0, 0),
        (-2147483648.0, 1),
        (-2147483904.0, 0),
        (4294967296.0, 0),
        (math.inf, 0),
        (-math.inf, 0),
        (65535.4, 65535),
        (0.5, 1),
        (-0.5, 0),
        (-1.5, 65535),
        (1.4999999, 1),
        (0.49999997, 0),
    ]
    v = np.zeros((1, 1, 1), np.uint16)
    f = np.ones((1, 1), np.float32)
    b = np.zeros((1, 1), np.float32)
    for w, want in cases:
        for fn in (ff.flat_field_correction_3d, ff.flat_field_correction_3d_ref):
            got = fn(v, f, b, const_offset=w)
            assert int(got[0, 0, 0]) == want, (w, fn)
    # through the background too (no offset): max((0 - b) / 1, 0) + b is b
    w = np.array([[c[0] for c in cases if c[0] >= 0 or math.isnan(c[0])]], np.float32)
    got = ff.flat_field_correction_3d(np.zeros((2, *w.shape), np.uint16), np.ones_like(w), w)
    want = [c[1] for c in cases if c[0] >= 0 or math.isnan(c[0])]
    assert got[1, 0].tolist() == want


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_results_do_not_depend_on_threads(cls, monkeypatch):
    monkeypatch.setattr(ff, "_threads", lambda n: 1 if n is None else n)
    monkeypatch.setattr(ff, "_SERIAL_VOXELS", 0)
    for shape in [(23, 97, 61), (1, 300, 7), (40, 2, 3), (5, 3, 9001)]:
        v = native(shape, cls, "thr", hi=65535)
        f = flat(shape[1:], "thr", lo=0.01)
        b = background(shape[1:], "thr")
        for cast, off in ((True, None), (False, None), (True, -3.0)):
            ref = ff.flat_field_correction_3d_ref(v, f, b, cast_data_type=cast, const_offset=off)
            for n in (1, 2, 3, 4, 7, 16, 64):
                got = ff.flat_field_correction_3d(
                    v, f, b, cast_data_type=cast, const_offset=off, num_threads=n
                )
                assert same_bits(got, ref), (shape, cast, off, n)


@pytest.mark.parametrize(
    "pair", [("uint16", np.uint16), ("uint16", np.float32), ("single", np.float32)]
)
def test_every_instruction_set_gives_the_baseline_bits(pair):
    cls, out_dtype = pair
    ext = ff._ext()
    isas = list(range(ext.best_isa() + 1))
    v = native((6, 13, 77), cls, "isa", hi=65535)
    if cls == "single":
        v[0, 0, :9] = [NAN, np.inf, -np.inf, -0.0, 1e30, -1e30, 0.0, 3e9, -3e9]
    f = special_flat((13, 77))
    b = special_bg((13, 77))
    for has, off in ((False, 0.0), (True, -3.25), (True, 1e10), (True, math.nan)):
        runs = []
        for isa in isas:
            out = np.empty(v.shape, out_dtype)
            ext.flat_field_correct(v, f, b, out, has, off, 3, isa)
            runs.append(out)
        want = ff._ref(v, f, b, np.dtype(out_dtype), np.float32(off) if has else None)
        for isa, out in zip(isas, runs, strict=True):
            assert same_bits(out, want), (isa, has, off)
    with pytest.raises(ValueError):
        ext.flat_field_correct(
            v, f, b, np.empty(v.shape, out_dtype), False, 0.0, 1, ext.best_isa() + 1
        )


def test_binding_refuses_other_layouts_and_classes():
    ext = ff._ext()
    v = np.zeros((2, 4, 3), np.uint16)
    f = np.ones((4, 3), np.float32)
    b = np.zeros((4, 3), np.float32)
    out = np.empty(v.shape, np.uint16)
    with pytest.raises(TypeError):
        ext.flat_field_correct(np.asfortranarray(v), f, b, out, False, 0.0, 1)
    with pytest.raises(TypeError):
        ext.flat_field_correct(v, f.astype(np.float64), b, out, False, 0.0, 1)
    with pytest.raises(TypeError):
        ext.flat_field_correct(v.astype(np.float32), f, b, out, False, 0.0, 1)
    with pytest.raises(ValueError):
        ext.flat_field_correct(v, f[:3], b, out, False, 0.0, 1)
    with pytest.raises(ValueError):
        ext.flat_field_correct(v, f, b, out[:1], False, 0.0, 1)
    with pytest.raises(ValueError):
        ext.flat_field_correct(v, f, b, out, False, 0.0, 0)


def test_layout_and_input_are_respected():
    v = native((9, 15, 14), "uint16", "layout")
    f = flat((15, 14), "layout")
    b = background((15, 14), "layout")
    keep = (v.copy(), f.copy(), b.copy())
    want = ff.flat_field_correction_3d(v, f, b)
    assert want.flags.c_contiguous and want.dtype == np.uint16 and want.shape == v.shape
    assert same_bits(
        ff.flat_field_correction_3d(np.asfortranarray(v), np.asfortranarray(f), b), want
    )
    sliced = native((9, 30, 14), "uint16", "layout2")[:, ::2]
    assert same_bits(
        ff.flat_field_correction_3d(sliced, f, b), ff.flat_field_correction_3d(sliced.copy(), f, b)
    )
    for a, k in zip((v, f, b), keep, strict=True):
        assert same_bits(a, k)
    # an image stays an image; a 4-D array with a leading 1 is its volume
    assert same_bits(ff.flat_field_correction_3d(v[0], f, b), want[0])
    assert same_bits(ff.flat_field_correction_3d(v[None], f, b), want[None])
    # a (1, ny, nx) flat field is a 2-D image to MATLAB
    assert same_bits(ff.flat_field_correction_3d(v, f[None], b[None]), want)
    # single output for uint16 data without the cast, single for single data
    assert ff.flat_field_correction_3d(v, f, b, cast_data_type=False).dtype == np.float32
    assert ff.flat_field_correction_3d(v.astype(np.float32), f, b).dtype == np.float32
    assert ff.flat_field_correction_3d(np.zeros((0, 15, 14), np.uint16), f, b).shape == (0, 15, 14)


def test_reference_keeps_negative_zero():
    """-0 passes max(., 0) (-0 >= 0) and -0 + -0 is -0, in the kernel and the reference"""
    v = np.full((2, 3, 4), -0.0, np.float32)
    f = np.ones((3, 4), np.float32)
    b = np.zeros((3, 4), np.float32)
    for off in (-0.0, 0.0, None):
        want = by_hand(v, f, b, True, off)
        assert same_bits(ff.flat_field_correction_3d(v, f, b, const_offset=off), want)
        assert same_bits(ff.flat_field_correction_3d_ref(v, f, b, const_offset=off), want)
    assert np.signbit(ff.flat_field_correction_3d_ref(v, f, b, const_offset=-0.0)).all()


def test_results_go_into_out():
    v = native((4, 9, 8), "uint16", "out")
    f = flat((9, 8), "out")
    b = background((9, 8), "out")
    want = ff.flat_field_correction_3d(v, f, b, const_offset=3.0)
    out = np.full(v.shape, 7, np.uint16)
    assert ff.flat_field_correction_3d(v, f, b, const_offset=3.0, out=out) is out
    assert same_bits(out, want)
    img = np.empty((9, 8), np.float32)
    ff.flat_field_correction_3d(v[1], f, b, cast_data_type=False, out=img)
    assert same_bits(img, ff.flat_field_correction_3d(v[1], f, b, cast_data_type=False))
    for bad in (
        np.empty(v.shape, np.float32),
        np.empty((4, 9, 7), np.uint16),
        np.empty((4, 8, 9), np.uint16).transpose(0, 2, 1),
    ):
        with pytest.raises(ValueError):
            ff.flat_field_correction_3d(v, f, b, out=bad)
    x = v.astype(np.float32)
    with pytest.raises(ValueError):
        ff.flat_field_correction_3d(x, f, b, out=x)


def test_reference_runs_without_the_extension(monkeypatch):
    v = native((4, 9, 8), "uint16", "noext")
    f = flat((9, 8), "noext")
    b = background((9, 8), "noext")
    want = ff.flat_field_correction_3d(v, f, b, const_offset=2.0)
    monkeypatch.setattr(ff, "_ext", lambda: None)
    assert not ff.extension_available()
    assert same_bits(ff.flat_field_correction_3d(v, f, b, const_offset=2.0), want)


def test_mx_scalar_reads_the_first_element():
    assert ff._mx_scalar(np.array([[0.0, 1.0]])) == 0.0
    assert ff._mx_scalar(np.array([[2.0], [3.0]])) == 2.0
    assert ff._mx_scalar(np.asfortranarray(np.array([[4.0, 5.0], [6.0, 7.0]]))) == 4.0
    assert ff._mx_scalar(True) == 1.0
    assert ff._mx_scalar(np.int16(-7)) == -7.0
    assert ff._mx_scalar(np.array([], np.float64)) == 0.0
    assert ff._mx_scalar(2 + 3j) == 2.0
    assert ff._mx_scalar("a") == 97.0
    assert math.isnan(ff._mx_scalar(math.nan))


# ---------------------------------------------------------------------------------------------
# XR_LSFlatFieldCorrection's step: the mex or MATLAB's code
# ---------------------------------------------------------------------------------------------


def test_apply_uses_the_mex_for_uint16_and_single(monkeypatch):
    v = native((3, 8, 6), "uint16", "apply")
    f = flat((8, 6), "apply")
    b = background((8, 6), "apply")
    want = ff.flat_field_correction_3d(v, f, b)
    # uint16 images are made single first (nonSingleFlatFieldUsesFastPath's regression)
    calls = []
    real = ff._matlab_code
    monkeypatch.setattr(ff, "_matlab_code", lambda *a: calls.append(1) or real(*a))
    got = ff.apply_flat_field(v, (np.floor(f * 10) + 1).astype(np.uint16), b.astype(np.uint16))
    assert got.dtype == np.uint16 and not calls
    assert same_bits(
        got,
        ff.flat_field_correction_3d(
            v, (np.floor(f * 10) + 1).astype(np.float32), b.astype(np.uint16).astype(np.float32)
        ),
    )
    assert same_bits(ff.apply_flat_field(v, f, b), want)
    assert same_bits(ff.apply_flat_field(v, f, b, const_offset=[]), want)  # [] is no offset
    ff.apply_flat_field(v.astype(np.float64), f, b)
    assert calls


def test_matlab_code_casts_with_saturation():
    """the MATLAB code (other classes) rounds half away from zero and saturates where the mex
    wraps: (250 - 100) / 0.25 + 100 is 700, uint8 255; uint16 through the mex wraps"""
    f = np.array([[0.5, 0.25, 1.0]], np.float32)
    b = np.full((1, 3), 100.0, np.float32)
    got = ff.apply_flat_field(np.full((1, 1, 3), 250, np.uint8), f, b)
    assert got.dtype == np.uint8 and got.tolist() == [[[255, 255, 250]]]
    v = np.full((1, 1, 3), 300, np.int16)
    got = ff.apply_flat_field(v, f, b, const_offset=-1000.25)
    assert got.dtype == np.int16 and got.tolist() == [[[-600, -200, -800]]]
    got = ff.apply_flat_field(v, f, b, const_offset=-1000.5)
    assert got.tolist() == [[[-601, -201, -801]]]  # half away from zero
    got = ff.apply_flat_field(v.astype(np.float64), f, b, cast_data_type=False)
    assert got.dtype == np.float32 and got.tolist() == [[[500.0, 900.0, 300.0]]]
    got = ff.apply_flat_field(v.astype(np.float64), f, b)
    assert got.dtype == np.float64 and got.tolist() == [[[500.0, 900.0, 300.0]]]
    u = np.full((1, 1, 3), 65535, np.uint16)
    tiny = np.full((1, 3), 0.25, np.float32)
    zero = np.zeros((1, 3), np.float32)
    assert ff.apply_flat_field(u, tiny, zero).tolist() == [[[65532, 65532, 65532]]]  # wraps
    with pytest.raises(ff.MatlabError) as e:
        ff.apply_flat_field(v.astype(np.float64), f, b, mex_only=True)
    assert e.value.identifier == TYPE_ERROR
    with pytest.raises(ff.MatlabError) as e:
        ff.apply_flat_field(v.astype(np.float64), f[:, :2], b[:, :2])
    assert e.value.identifier == "MATLAB:sizeDimensionsMustMatch"


def test_matlab_code_broadcasts_four_dimensions():
    v = native((3, 5, 4), "uint16", "four")
    f = flat((5, 4), "four")
    b = background((5, 4), "four")
    v4 = np.stack([v, v[::-1]])
    got = ff.apply_flat_field(v4, f, b)
    assert got.shape == v4.shape and got.dtype == np.uint16
    assert same_bits(got[0], ff.flat_field_correction_3d(v, f, b))
