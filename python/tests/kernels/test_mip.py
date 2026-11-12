"""Max pooling, range maximum and box minimum (petakit5d.kernels.mip).

Oracle cases (MATLAB: max_pooling_3d_mex, max_pooling_3d, range_max_3d_mex, min_bbox_3d_mex,
min_bbox_3d) over every class, sizes not divisible by the pool size, 2-D, 4-D, empty and
all-zero input, NaN, -0 and Inf, and the input errors. Their inputs are MATLAB's (y, x, z)
arrays and values; the roles give the Python calls (z, y, x) arrays, (pz, py, px) pools,
axes, 0-based starts and Python boxes, and the MATLAB-value forms (``*_matlab``) take MATLAB's
pool sizes and boxes where the mex and MATLAB's code read them as numbers. Then the MATLAB
tests' checks (tests/tools/MIPTest.m, BlockMIPTest.m, tests/platform/MexInputCheckTest.m) on
(z, y, x) data against NumPy references, and results at several thread counts.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, matlab_order, size_from_matlab, to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.kernels import mip
from petakit5d.testing import oracle
from petakit5d.testing.cases import Box, Case, Dim, In, Index1, PerAxis, Value, Vol, mcol
from petakit5d.testing.synth import CounterRNG

MOD = "petakit5d.kernels.mip"
CLASSES = ("logical", "uint8", "uint16", "single", "double")
FALLBACK_CLASSES = ("int8", "int16", "int32", "int64", "uint32", "uint64")
DTYPES = {
    "logical": np.bool_,
    "uint8": np.uint8,
    "uint16": np.uint16,
    "single": np.float32,
    "double": np.float64,
    "int8": np.int8,
    "int16": np.int16,
    "int32": np.int32,
    "int64": np.int64,
    "uint32": np.uint32,
    "uint64": np.uint64,
}
# value ranges of the random integer volumes (signed classes have negative values)
INT_RANGES = {
    "uint8": (0, 255),
    "uint16": (0, 65535),
    "int8": (-128, 127),
    "int16": (-1000, 1000),
    "int32": (-1000000, 1000000),
    "int64": (-(10**9), 10**9),
    "uint32": (0, 4 * 10**9),
    "uint64": (0, 4 * 10**9),
}

# ---------------------------------------------------------------------------------------------
# input generators: MATLAB's (y, x, z) arrays (the case keys are made from these bytes)
# ---------------------------------------------------------------------------------------------


def vol(shape, cls, name, positive=False):
    """a deterministic random array of a MATLAB class (Fortran order); floats are normal * 1000
    (half negative) or, with positive, uniform * 1000"""
    # CounterRNG fills (z, y, x) arrays in MATLAB's linear order: the MATLAB form of its array
    # of our shape is the array it gave for MATLAB's shape before the (z, y, x) change
    rng = CounterRNG(2, f"mip/{name}/{cls}/{'x'.join(map(str, shape))}")
    ps = size_from_matlab(shape)
    if cls == "logical":
        return to_matlab(rng.uniform(ps) > 0.5)
    dt = np.dtype(DTYPES[cls])
    if dt.kind == "f":
        v = rng.uniform(ps) * 1000.0 if positive else rng.normal(ps) * 1000.0
        return to_matlab(v.astype(dt))
    lo, hi = INT_RANGES[cls]
    return to_matlab(rng.integers(lo, hi, ps, dtype=dt))


def specials(cls, shape=(9, 8, 7)):
    """single or double with NaN of several payloads, -0, +0, Inf, -Inf, a block of negative
    values only and lines of -0 and 0"""
    dt = np.dtype(DTYPES[cls])
    a = vol(shape, cls, "specials")
    bits = np.uint32 if dt == np.float32 else np.uint64
    if dt == np.float32:
        nans = np.array([0xFFC00000, 0x7FC00000, 0xFFC00123], dtype=bits).view(dt)
    else:
        nans = np.array(
            [0xFFF8000000000000, 0x7FF8000000000000, 0xFFF8000000000123], dtype=bits
        ).view(dt)
    flat = a.reshape(-1, order="F")
    rng = CounterRNG(3, f"specials/{cls}")
    pick = rng.integers(0, 9, flat.size, dtype=np.int64).reshape(-1)
    flat[pick == 0] = nans[0]
    flat[pick == 1] = nans[1]
    flat[pick == 2] = nans[2]
    a[2, :, :] = nans[0]  # a line of NaN
    a[:, 3, :] = -0.0
    a[:, 3, ::3] = 0.0
    a[5, 5, :] = np.inf
    a[6, 6, :] = -np.inf
    a[:2, :2, :2] = -np.abs(a[:2, :2, :2]) - 1  # negative only
    a[4:6, 6:8, 2:4] = -0.0
    a[4, 6, 3] = 0.0
    return a


def _value(make):
    return Value(make)


def _mex(fn, **kw):
    """a case of a mex function, run by the MATLAB-value form alone (mex_only)"""
    return dict(python=f"{MOD}:{fn}_matlab", py_kwargs={"mex_only": True}, **kw)


# ---------------------------------------------------------------------------------------------
# oracle cases: max pooling
# ---------------------------------------------------------------------------------------------

POOL_SHAPES = [
    ((23, 19, 11), (2, 2, 2)),
    ((23, 19, 11), (3, 2, 1)),
    ((23, 19, 11), (4, 5, 5)),
    ((40, 36, 12), (2, 2, 2)),
    ((7, 5, 3), (30, 30, 30)),
    ((9, 6, 5), (1, 1, 1)),
    ((17, 13), (3, 4, 1)),
    ((17, 13), (2, 2, 3)),
    ((1, 13, 4), (2, 3, 2)),
    ((23, 19, 11), (2, 3, 2.5)),
    ((6, 5, 4, 2), (2, 2, 2)),
    ((0, 5, 3), (2, 2, 2)),
    ((4, 0, 0), (2, 2, 2)),
    ((5, 4, 0), (1, 1, 1)),
]


def _shape_id(shape, pool):
    s = "x".join(map(str, shape))
    p = "x".join(str(v).replace(".", "p") for v in pool)
    return f"{s}_p{p}"


def _pool_native(case_id, pool, make):
    """a mex case the native function runs: (z, y, x) volume, (pz, py, px) pool"""
    return Case(
        case_id,
        matlab="max_pooling_3d_mex",
        args=[Vol(In("im")), PerAxis(list(pool), "yxz")],
        python=f"{MOD}:max_pooling_3d",
        inputs={"im": _value(make)},
        out_roles=("vol",),
    )


def _pool_wrapper(case_id, pool, make, **kw):
    """a case of max_pooling_3d.m: MATLAB's pool values"""
    return Case(
        case_id,
        matlab="max_pooling_3d",
        args=[Vol(In("im")), pool],
        python=f"{MOD}:max_pooling_3d_matlab",
        inputs={"im": _value(make)},
        out_roles=("vol",),
        **kw,
    )


def _pool_cases():
    cases = []
    for cls in CLASSES:
        for shape, pool in POOL_SHAPES:
            sid = _shape_id(shape, pool)
            make = lambda s=shape, c=cls: vol(s, c, "pool")  # noqa: E731
            cases.append(_pool_native(f"max_pooling_3d_mex.{cls}.{sid}", pool, make))
            cases.append(_pool_wrapper(f"max_pooling_3d.{cls}.{sid}", list(pool), make))
        cases.append(
            _pool_native(
                f"max_pooling_3d_mex.{cls}.zeros",
                (3, 2, 2),
                lambda c=cls: np.zeros((10, 7, 5), DTYPES[c], order="F"),
            )
        )
        cases.append(
            _pool_native(
                f"max_pooling_3d_mex.{cls}.threads",
                (3, 4, 2),
                lambda c=cls: vol((301, 257, 33), c, "big", positive=True),
            )
        )
    for cls in ("single", "double"):
        for pool in ((2, 2, 2), (3, 1, 2), (1, 1, 1), (9, 8, 7)):
            cases.append(
                _pool_native(
                    f"max_pooling_3d_mex.{cls}.specials_p{'x'.join(map(str, pool))}",
                    pool,
                    lambda c=cls: specials(c),
                )
            )
    # the known issue: blocks of negative values give 0 (MIPTest.maxPoolingOfNegativeValues)
    cases.append(
        _pool_native(
            "max_pooling_3d_mex.single.negative_block",
            (2, 2, 2),
            lambda: np.asfortranarray(
                -np.arange(1, 9, dtype=np.float32).reshape((2, 2, 2), order="F")
            ),
        )
    )
    # MexInputCheckTest.mipMexCheckTheirInputs: a fractional pool size is truncated
    cases.append(
        _pool_native(
            "max_pooling_3d_mex.single.input_check_truncates",
            (2, 3, 2.5),
            lambda: np.asfortranarray(
                np.arange(1, 25, dtype=np.float32).reshape((4, 3, 2), order="F")
            ),
        )
    )
    # other integer classes: the MATLAB code of max_pooling_3d.m
    for cls in FALLBACK_CLASSES:
        for shape, pool in (
            ((23, 19, 11), (4, 5, 5)),
            ((23, 19, 11), (2, 2, 2)),
            ((8, 6, 4), (2, 3, 2)),
            ((5, 19, 11), (8, 5, 5)),
            ((7, 5, 3), (30, 30, 30)),
            ((0, 5, 3), (2, 2, 2)),
        ):
            cases.append(
                _pool_wrapper(
                    f"max_pooling_3d.{cls}.{_shape_id(shape, pool)}",
                    list(pool),
                    lambda s=shape, c=cls: vol(s, c, "pool"),
                )
            )
        cases.append(
            Case(
                f"max_pooling_3d_mex.{cls}.class_error",
                matlab="max_pooling_3d_mex",
                args=[Vol(In("im")), [2, 2, 2]],
                inputs={"im": _value(lambda c=cls: vol((6, 5, 4), c, "pool"))},
                expect_error="mat:dataTypeError",
                **_mex("max_pooling_3d"),
            )
        )
    # errors
    bad_pools = {
        "two": [2, 2],
        "four": [2, 2, 2, 2],
        "zero": [0, 1, 1],
        "half": [0.5, 1, 1],
        "nan": [2, NAN, 1],
        "negative": [-2, 2, 2],
        "huge": [9.5e15, 1, 1],
    }
    for name, pool in bad_pools.items():
        cases.append(
            Case(
                f"max_pooling_3d_mex.single.pool_{name}",
                matlab="max_pooling_3d_mex",
                args=[Vol(In("im")), pool],
                inputs={"im": _value(lambda: vol((6, 5, 4), "single", "pool"))},
                expect_error="pooling:inputError",
                **_mex("max_pooling_3d"),
            )
        )
    cases.append(
        Case(
            "max_pooling_3d_mex.single.pool_column",
            matlab="max_pooling_3d_mex",
            args=[Vol(In("im")), mcol([2, 2, 2])],
            py_args=[Vol(In("im")), np.array([[2.0], [2.0], [2.0]])],
            inputs={"im": _value(lambda: vol((6, 5, 4), "single", "pool"))},
            expect_error="pooling:inputError",
            **_mex("max_pooling_3d"),
        )
    )
    for name, (shape, pool) in {
        "int16_2d": ((6, 5), (2, 2, 1)),
        "int16_4d": ((6, 5, 4, 2), (2, 2, 2)),
        "int16_fraction_nopad": ((5, 4, 4), (2.5, 2, 2)),
        "int16_fraction_pad": ((6, 4, 4), (2.5, 2, 2)),
        "int16_pool_zero": ((6, 4, 4), (0, 2, 2)),
    }.items():
        cls = name.split("_")[0]
        cases.append(
            _pool_wrapper(
                f"max_pooling_3d.error_{name}",
                list(pool),
                lambda s=shape, c=cls: vol(s, c, "pool"),
                expect_error=True,
            )
        )
    for name, pool in {
        "two": [2, 2],
        "zero": [0, 1, 1],
        "half": [0.5, 1, 1],
        "negative": [-2, 2, 2],
        "negative_even": [-3, 2, 2],
        "inf": [float("inf"), 1, 1],
        "scalar": 2,
    }.items():
        cases.append(
            _pool_wrapper(
                f"max_pooling_3d.error_uint16_pool_{name}",
                pool,
                lambda: vol((6, 5, 4), "uint16", "pool"),
                expect_error=True,
            )
        )
    return cases


# ---------------------------------------------------------------------------------------------
# oracle cases: range maximum
# ---------------------------------------------------------------------------------------------


def _starts_sets(n):
    sets = {"one": [1], "quarters": sorted({1, min(2, n), -(-n // 2), n})}
    if n > 2:
        sets["odd"] = list(range(1, n + 1, 2))
        sets["each"] = list(range(1, n + 1))
    return sets


def _range_case(case_id, make, dim, starts, **kw):
    return Case(
        case_id,
        matlab="range_max_3d_mex",
        args=[Vol(In("im")), Dim(dim), Index1(starts)],
        python=f"{MOD}:range_max_3d",
        inputs={"im": _value(make)},
        out_roles=("vol",),
        **kw,
    )


def _range_cases():
    cases = []
    for cls in ("uint8", "uint16", "single", "double"):
        for shape in ((13, 11, 9), (17, 9), (1, 12, 7), (4100, 5, 3), (150, 131, 5)):
            sid = "x".join(map(str, shape))
            for dim in (1, 2, 3):
                n = shape[dim - 1] if dim <= len(shape) else 1
                for sname, starts in _starts_sets(n).items():
                    if shape[0] > 1000 and sname == "each":
                        continue
                    cases.append(
                        _range_case(
                            f"range_max_3d.{cls}.{sid}.d{dim}.{sname}",
                            lambda s=shape, c=cls: vol(s, c, "range"),
                            dim,
                            starts,
                        )
                    )
    for cls in ("single", "double"):
        for dim in (1, 2, 3):
            for sname, starts in {"one": [1], "two": [1, 3], "each": list(range(1, 8))}.items():
                cases.append(
                    _range_case(
                        f"range_max_3d.{cls}.specials.d{dim}.{sname}",
                        lambda c=cls: specials(c),
                        dim,
                        starts,
                    )
                )
    for cls in ("uint8", "uint16", "single", "double"):
        cases.append(
            _range_case(
                f"range_max_3d.{cls}.zeros",
                lambda c=cls: np.zeros((6, 5, 4), DTYPES[c], order="F"),
                3,
                [1, 3],
            )
        )
    magic4 = np.array(
        [[16, 2, 3, 13], [5, 11, 10, 8], [9, 7, 6, 12], [4, 14, 15, 1]], dtype=np.uint16, order="F"
    )
    # MATLAB's arguments, and the Python axis and starts they stand for
    errors = {
        "not_first": (magic4, 1, [2, 3], -2, [1, 2]),
        "not_increasing": (magic4, 1, [1, 3, 2], -2, [0, 2, 1]),
        "past_end": (magic4, 1, [1, 5], -2, [0, 4]),
        "fraction": (magic4, 1, [1, 2.5], -2, [0, 1.5]),
        "dim4": (magic4, 4, [1], -4, [0]),
        "dim_fraction": (magic4, 1.5, [1], 1.5, [0]),
        "dim_two": (magic4, [1, 2], [1], [-2, -1], [0]),
        "starts_empty": (magic4, 1, np.zeros((1, 0)), -2, []),
        "starts_nan": (magic4, 1, [1, NAN], -2, [0, NAN]),
        "int16": (magic4.astype(np.int16), 1, [1], -2, [0]),
        "logical": (magic4 > 5, 1, [1], -2, [0]),
        "empty": (np.zeros((0, 3), np.uint16), 2, [1], -1, [0]),
        "four_dims": (np.zeros((2, 2, 2, 2), np.uint16, order="F"), 1, [1], -2, [0]),
    }
    for name, (im, dim, starts, axis, starts0) in errors.items():
        cases.append(
            Case(
                f"range_max_3d.error_{name}",
                matlab="range_max_3d_mex",
                args=[Vol(In("im")), dim, starts],
                py_args=[Vol(In("im")), axis, starts0],
                python=f"{MOD}:range_max_3d",
                inputs={"im": _value(lambda a=im: np.asfortranarray(a))},
                expect_error="range_max_3d:inputError",
            )
        )
    return cases


# ---------------------------------------------------------------------------------------------
# oracle cases: minimum in a box
# ---------------------------------------------------------------------------------------------

BOXES = {
    "region": [3, 4, 2, 15, 25, 9],
    "empty": [5, 1, 1, 4, 30, 10],
    "full": [1, 1, 1, 20, 30, 10],
    "voxel": [7, 8, 3, 7, 8, 3],
    "fraction": [2.7, 3.2, 1.9, 10.5, 20.9, 8.1],
}
BAD_BOXES = {
    "start_zero": [0, 1, 1, 5, 5, 5],
    "past_end": [1, 1, 1, 21, 5, 5],
    "reversed": [5, 1, 1, 3, 30, 10],
    "nan_start": [NAN, 1, 1, 5, 5, 5],
    "inf_end": [1, 1, 1, float("inf"), 5, 5],
    "negative": [-1, 1, 1, 5, 5, 5],
    "fraction_start": [0.5, 1, 1, 5, 5, 5],
    "five": [1, 1, 1, 5, 5],
    "seven": [3, 4, 2, 15, 25, 9, 1],
}


def _min_input(cls, nan=False):
    def make():
        if cls == "logical":
            return np.asfortranarray(vol((20, 30, 10), "double", "min", positive=True) > 50.0)
        a = vol((20, 30, 10), "double", "min", positive=True) / 5.0 + 1.0
        if nan:
            a[1:9, :, :] = NAN
            a[11, 6, 2] = 0.5
            a[12, 9, 3] = -0.0
            a[13, 9, 3] = 0.0
        if np.dtype(DTYPES[cls]).kind == "f":
            return np.asfortranarray(a.astype(DTYPES[cls]))
        # integer classes: NaN as 0, the rest floored (exact integers below 2^16)
        a = np.floor(np.where(np.isnan(a), 0.0, a))
        return np.asfortranarray(a.astype(DTYPES[cls]))

    return make


COLON_ENDS = [
    (2.7, 10.5),
    (2.5, 10.5),
    (2.4, 9.6),
    (3.5, 9.5),
    (1.5, 20.4),
    (0.6, 5),
    (1.2, 1.3),
    (2.6, 2.5),
    (19.7, 20.6),
    (0.4, 3),
    (4, 7.5),
]
COLON_ERRORS = {(0.4, 3)}


def _num(v):
    return str(v).replace(".", "p")


def _ramp(values):
    return np.asfortranarray(np.broadcast_to(values.astype(np.int16)[:, None, None], (20, 3, 2)))


def _ramp_up():
    return _ramp(np.arange(1, 21))


def _ramp_down():
    return _ramp(100 - np.arange(1, 21))


def _integral(box):
    return all(np.isfinite(float(v)) and float(v) == int(v) for v in box)


def _min_case(fn, case_id, box, make, *, py_box=None, **kw):
    """min_bbox_3d_mex with a valid integral box: the native function with a Python box; else
    the MATLAB-value form (the mex alone for min_bbox_3d_mex)"""
    if fn == "min_bbox_3d_mex" and not kw.get("expect_error") and _integral(box):
        return Case(
            case_id,
            matlab=fn,
            args=[Vol(In("im")), Box(box)],
            python=f"{MOD}:min_bbox_3d",
            inputs={"im": _value(make)},
            **kw,
        )
    extra = (
        _mex("min_bbox_3d") if fn == "min_bbox_3d_mex" else {"python": f"{MOD}:min_bbox_3d_matlab"}
    )
    if py_box is not None:
        extra["py_args"] = [Vol(In("im")), py_box]
    return Case(
        case_id,
        matlab=fn,
        args=[Vol(In("im")), box],
        inputs={"im": _value(make)},
        **extra,
        **kw,
    )


def _min_cases():
    cases = []
    for cls in (*CLASSES, *FALLBACK_CLASSES):
        for bname, box in BOXES.items():
            for fn in ("min_bbox_3d_mex", "min_bbox_3d"):
                if fn == "min_bbox_3d_mex" and cls in FALLBACK_CLASSES:
                    continue
                if fn == "min_bbox_3d" and cls in FALLBACK_CLASSES and bname == "fraction":
                    continue  # MATLAB's indexing error, a case below
                cases.append(_min_case(fn, f"{fn}.{cls}.{bname}", box, _min_input(cls)))
    for cls in ("single", "double"):
        for bname, box in {
            "all": [1, 1, 1, 20, 30, 10],
            "nan_only": [2, 1, 1, 9, 30, 10],
            "zeros": [12, 9, 3, 14, 10, 4],
            "neg_zero_first": [13, 10, 4, 14, 10, 4],
        }.items():
            for fn in ("min_bbox_3d_mex", "min_bbox_3d"):
                cases.append(
                    _min_case(fn, f"{fn}.{cls}.nan_{bname}", box, _min_input(cls, nan=True))
                )
    for cls in ("uint16", "single"):
        cases.append(
            _min_case(
                "min_bbox_3d_mex",
                f"min_bbox_3d_mex.{cls}.image",
                [2, 3, 1, 9, 7, 1],
                lambda c=cls: vol((11, 9), c, "min2d", positive=True),
            )
        )
        cases.append(
            _min_case(
                "min_bbox_3d_mex",
                f"min_bbox_3d_mex.{cls}.four_dims",
                [1, 1, 1, 6, 5, 4],
                lambda c=cls: vol((6, 5, 4, 2), c, "min4d", positive=True),
            )
        )
        cases.append(
            _min_case(
                "min_bbox_3d_mex",
                f"min_bbox_3d_mex.{cls}.empty_array",
                [1, 1, 1, 0, 5, 3],
                lambda c=cls: np.zeros((0, 5, 3), DTYPES[c], order="F"),
            )
        )
        column = np.array([[3.0], [4], [2], [15], [25], [9]])
        cases.append(
            _min_case(
                "min_bbox_3d",
                f"min_bbox_3d.{cls}.column_box",
                mcol([3, 4, 2, 15, 25, 9]),
                _min_input(cls, nan=True),
                py_box=column,
            )
        )
        cases.append(
            _min_case(
                "min_bbox_3d_mex",
                f"min_bbox_3d_mex.{cls}.column_box",
                mcol([3, 4, 2, 15, 25, 9]),
                _min_input(cls, nan=True),
                py_box=column,
                expect_error="minbbox:inputError",
            )
        )
    # boxes the mex rejects that the MATLAB code (crop3d's indexing) takes
    returns = {"reversed", "fraction_start", "seven"}
    for cls in ("uint16", "single", "int16"):
        for bname, box in BAD_BOXES.items():
            # (uint64_t)Inf is 0 in the mex: an empty box; MATLAB's colon refuses Inf
            ok = bname in returns or (bname == "inf_end" and cls != "int16")
            cases.append(
                _min_case(
                    "min_bbox_3d",
                    f"min_bbox_3d.{cls}.bad_{bname}",
                    box,
                    _min_input(cls),
                    expect_error=not ok,
                )
            )
            if cls != "int16":
                cases.append(
                    _min_case(
                        "min_bbox_3d_mex",
                        f"min_bbox_3d_mex.{cls}.bad_{bname}",
                        box,
                        _min_input(cls),
                        expect_error=bname != "inf_end",
                    )
                )
    cases.append(
        _min_case(
            "min_bbox_3d", "min_bbox_3d.int16.fraction", BOXES["fraction"], _min_input("int16")
        )
    )
    # MATLAB's colon in an index with non-integer ends (a warning, then rounded): the first
    # index of a range from values increasing along y, the last from values decreasing
    for lo, hi in COLON_ENDS:
        for name, make in (("first", _ramp_up), ("last", _ramp_down)):
            cases.append(
                _min_case(
                    "min_bbox_3d",
                    f"min_bbox_3d.int16.colon_{name}_{_num(lo)}_{_num(hi)}",
                    [lo, 1, 1, hi, 3, 2],
                    make,
                    expect_error=(lo, hi) in COLON_ERRORS,
                )
            )
    for cls in ("uint16", "single"):
        cases.append(
            _min_case(
                "min_bbox_3d",
                f"min_bbox_3d.{cls}.column_box_neg_zero_first",
                mcol([13, 10, 4, 14, 10, 4]),
                _min_input(cls, nan=True),
                py_box=np.array([[13.0], [10], [4], [14], [10], [4]]),
            )
        )
    cases.append(
        _min_case(
            "min_bbox_3d_mex",
            "min_bbox_3d_mex.int16.class_error",
            BOXES["region"],
            _min_input("int16"),
            expect_error="mat:dataTypeError",
        )
    )
    cases.append(
        _min_case(
            "min_bbox_3d_mex",
            "min_bbox_3d_mex.single.threads",
            [3, 2, 5, 299, 250, 31],
            lambda: vol((301, 257, 33), "single", "big"),
        )
    )
    return cases


CASES = _pool_cases() + _range_cases() + _min_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# NumPy references on (z, y, x) volumes
# ---------------------------------------------------------------------------------------------


def native(shape, cls, name, positive=False):
    """a deterministic random (z, y, x) volume (C order)"""
    return from_matlab(vol(tuple(shape[1:]) + (shape[0],), cls, name, positive))


def native_specials(cls, shape):
    """specials() as a (z, y, x) volume"""
    return from_matlab(specials(cls, tuple(shape[1:]) + (shape[0],)))


def block_max(v, pool):
    """MIPTest's blockMax with the mex's start at 0: max(0, max of each block); NaN skipped"""
    p = [int(x) for x in pool]
    m = [-(-s // q) for s, q in zip(v.shape, p, strict=True)]
    out = np.zeros(m, dtype=v.dtype)
    for k in range(m[0]):
        for i in range(m[1]):
            for j in range(m[2]):
                blk = v[
                    k * p[0] : (k + 1) * p[0], i * p[1] : (i + 1) * p[1], j * p[2] : (j + 1) * p[2]
                ]
                best = out[k, i, j]
                for x in blk.ravel():
                    if x > best:
                        best = x
                out[k, i, j] = best
    return out


def range_max(v, axis, starts):
    """BlockMIPTest's rangeMax, with MATLAB max's rules folded in order along the axis: NaN
    ignored (the first NaN of a range of NaN), the first of equal values"""
    n = v.shape[axis]
    bounds = [*starts, n]
    parts = []
    for r in range(len(starts)):
        m = np.take(v, [bounds[r]], axis=axis).copy()
        for i in range(bounds[r] + 1, bounds[r + 1]):
            x = np.take(v, [i], axis=axis)
            take = x > m
            if v.dtype.kind == "f":
                take |= np.isnan(m) & ~np.isnan(x)
            m = np.where(take, x, m)
        parts.append(m)
    return np.concatenate(parts, axis=axis)


def serial_min(v, box):
    """min_bbox_3d_mex's loop on a (z, y, x) box: from the largest value, take a value only if
    smaller, in MATLAB's order (y fastest, then x, then z)"""
    z0, y0, x0, z1, y1, x1 = box
    dt = v.dtype
    m = np.finfo(dt).max if dt.kind == "f" else (True if dt == np.bool_ else np.iinfo(dt).max)
    for x in matlab_order(v[z0:z1, y0:y1, x0:x1]).ravel():
        if x < m:
            m = x
    return dt.type(m)


def same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, a.shape, b.dtype, b.shape)
    return np.array_equal(
        np.ascontiguousarray(a).view(np.uint8), np.ascontiguousarray(b).view(np.uint8)
    )


# ---------------------------------------------------------------------------------------------
# ported MATLAB tests
# ---------------------------------------------------------------------------------------------

# MIPTest's pool sizes [py px pz] as (pz, py, px)
POOL_SIZES = {"cube": (2, 2, 2), "anisotropic": (1, 3, 2), "nonDivisible": (5, 4, 5)}


@pytest.mark.parametrize("pool_type", ["uint16", "single"])
@pytest.mark.parametrize("pool_size", list(POOL_SIZES), ids=list(POOL_SIZES))
def test_max_pooling_matches_blockwise_max(pool_size, pool_type):
    """MIPTest.maxPoolingMatchesBlockwiseMax (rand * 1000: non-negative)"""
    v = native((11, 23, 19), pool_type, "blockwise", positive=True)
    pool = POOL_SIZES[pool_size]
    expected = block_max(v, pool)
    assert same_bits(mip.max_pooling_3d(v, pool), expected)
    py, px, pz = pool[1], pool[2], pool[0]
    assert same_bits(mip.max_pooling_3d_matlab(v, [py, px, pz]), expected)


def test_max_pooling_of_negative_values_gives_zero_known_issue():
    """MIPTest.maxPoolingOfNegativeValues (KNOWN ISSUE, kept): a block of negative values gives
    0, not its maximum (-1), as max_pooling_3d_mex"""
    v = -np.arange(1, 9, dtype=np.float32).reshape((2, 2, 2))
    out = mip.max_pooling_3d(v, (2, 2, 2))
    assert same_bits(out, np.zeros((1, 1, 1), np.float32))
    assert same_bits(mip.max_pooling_3d_matlab(v, [2, 2, 2]), out)
    # the MATLAB code of other classes: whole blocks give their maximum, edge blocks are padded
    # with 0. MATLAB's 3 x 2 x 2 of -1 .. -12 pools to 2 x 1 x 1, which MATLAB's code returns
    # as the row [-1 0] (max gives 1 x 2 x 1 x 1 x 1 x 1)
    m = -np.arange(1, 13, dtype=np.int16).reshape((3, 2, 2), order="F")
    v16 = from_matlab(m)
    assert np.array_equal(mip.max_pooling_3d(v16, (2, 2, 2)), np.array([[[-1], [0]]], np.int16))
    assert np.array_equal(mip.max_pooling_3d_matlab(v16, [2, 2, 2]), np.array([[-1, 0]], np.int16))


def test_mip_mex_check_their_inputs():
    """MexInputCheckTest.mipMexCheckTheirInputs (Python has no sparse arrays or classes of
    pool-size lists: those checks are not ported)"""
    # MATLAB's reshape(single(1 : 24), 4, 3, 2)
    v = from_matlab(np.arange(1, 25, dtype=np.float32).reshape((4, 3, 2), order="F"))
    for pool in ([1, 2, NAN], [1, 0.5, 1], [2, 2], [1, 2, 2, 1]):
        with pytest.raises(mip.MatlabError) as e:
            mip.max_pooling_3d(v, pool)
        assert e.value.identifier == "pooling:inputError"
    out = mip.max_pooling_3d_matlab(v, [2, 3, 2.5], mex_only=True)
    assert same_bits(out, np.array([22, 24], np.float32).reshape((1, 2, 1)))
    for args in ((v[0], -2, [0, 4]), (np.zeros((0, 3), np.float32), -2, [0])):
        with pytest.raises(mip.MatlabError) as e:
            mip.range_max_3d(*args)
        assert e.value.identifier == "range_max_3d:inputError"


@pytest.mark.parametrize("dtype", ["uint16", "single"])
def test_range_max_matches_max_of_each_range(dtype):
    """BlockMIPTest.rangeMaxMatchesMaxOfEachRange (sizes 13 x 11 x 9, 17 x 9, 1 x 12 x 7), with
    NaN, -0 and 0 for single"""
    for shape in ((9, 13, 11), (17, 9), (7, 1, 12)):
        v = (
            native(shape, dtype, "rangemax")
            if len(shape) == 3
            else from_matlab(vol(shape, dtype, "rangemax"))
        )
        if dtype == "single" and len(shape) == 3 and min(shape) > 7:
            v = native_specials("single", shape)
        v3 = v.reshape((1,) * (3 - v.ndim) + v.shape)
        for axis in (-2, -1, -3):
            n = v3.shape[axis]
            for starts in ([0], sorted({0, min(1, n - 1), -(-n // 2) - 1, n - 1})):
                out = mip.range_max_3d(v, axis, starts)
                want = range_max(v3, axis + 3, starts)
                assert same_bits(out, want if v.ndim == 3 else want[0]), (shape, axis, starts)


def test_range_max_rejects_invalid_input():
    """BlockMIPTest.rangeMaxRejectsInvalidInput (MATLAB's uint16(magic(4)), dim 1 is axis -2)"""
    v = np.array([[16, 2, 3, 13], [5, 11, 10, 8], [9, 7, 6, 12], [4, 14, 15, 1]], np.uint16)
    for args in (
        (v, -2, [1, 2]),
        (v, -2, [0, 2, 1]),
        (v, -2, [0, 4]),
        (v, -2, [0, 1.5]),
        (v, -4, [0]),
        (v.astype(np.int16), -2, [0]),
    ):
        with pytest.raises(mip.MatlabError) as e:
            mip.range_max_3d(*args)
        assert e.value.identifier == "range_max_3d:inputError"


@pytest.mark.parametrize("min_type", ["logical", "uint8", "uint16", "single", "double", "int16"])
def test_minimum_of_a_region(min_type):
    """MIPTest.minimumOfARegion: MATLAB's box [3 4 2 15 25 9] is (1, 2, 3, 9, 15, 25); the
    MATLAB code without the mex is the int16 path"""
    v = from_matlab(_min_input(min_type)())
    expected = v[1:9, 2:15, 3:25].min()
    out = mip.min_bbox_3d(v, (1, 2, 3, 9, 15, 25))
    assert out.dtype == v.dtype and out == expected
    assert mip.min_bbox_3d_matlab(v, [3, 4, 2, 15, 25, 9]) == expected
    largest = mip.min_bbox_3d(v, (0, 4, 0, 10, 4, 30))
    dt = v.dtype
    want = True if dt == np.bool_ else (np.finfo(dt).max if dt.kind == "f" else np.iinfo(dt).max)
    assert largest.dtype == dt and largest == want


def test_minimum_of_a_region_skips_nan():
    """MIPTest.minimumOfARegionSkipsNaN"""
    v = (native((10, 20, 30), "single", "skipnan", positive=True) / np.float32(1000)) + np.float32(
        1
    )
    v[:, 1:9, :] = NAN
    v[2, 11, 6] = 0.5
    big = np.finfo(np.float32).max
    for box, want in (
        ((0, 0, 0, 10, 20, 30), np.float32(0.5)),
        ((0, 1, 0, 10, 9, 30), big),
        ((0, 4, 0, 10, 4, 30), big),
    ):
        assert same_bits(mip.min_bbox_3d(v, box), np.float32(want))


# ---------------------------------------------------------------------------------------------
# threads, layouts, inputs
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("cls", ["logical", "uint16", "single", "double"])
def test_results_do_not_depend_on_threads(cls, monkeypatch):
    # thread counts beyond the job's CPUs too
    monkeypatch.setattr(mip, "_threads", lambda n: 1 if n is None else n)
    shape = (23, 97, 61)
    v = native_specials(cls, shape) if cls in ("single", "double") else native(shape, cls, "thr")
    box = (3, 1, 2, 20, 90, 60)
    ref_pool = mip.max_pooling_3d(v, (2, 3, 4), num_threads=1)
    assert same_bits(ref_pool, block_max(v, (2, 3, 4)))
    ref_min = mip.min_bbox_3d(v, box, num_threads=1)
    assert same_bits(ref_min, serial_min(v, box))
    for n in (2, 3, 4, 7, 16):
        assert same_bits(mip.max_pooling_3d(v, (2, 3, 4), num_threads=n), ref_pool)
        assert same_bits(mip.min_bbox_3d(v, box, num_threads=n), ref_min)
    if cls == "logical":
        return
    for axis in (0, 1, 2):
        starts = list(range(0, v.shape[axis], 5))
        ref = mip.range_max_3d(v, axis, starts, num_threads=1)
        assert same_bits(ref, range_max(v, axis, starts))
        for n in (2, 3, 4, 7, 16):
            assert same_bits(mip.range_max_3d(v, axis, starts, num_threads=n), ref)


@pytest.mark.parametrize("shape", [(5, 9, 700), (4, 40, 1500), (4, 8, 5000), (4, 300, 9)])
def test_range_max_parts_do_not_change_results(shape, monkeypatch):
    """rows along x are cut into parts by the thread count for the y maximum"""
    monkeypatch.setattr(mip, "_threads", lambda n: 1 if n is None else n)
    v = native_specials("single", shape)
    for axis in (0, 1):
        starts = [0] if shape[axis] < 3 else [0, 1]
        want = range_max(v, axis, starts)
        for n in (1, 2, 3, 7, 24, 64):
            assert same_bits(mip.range_max_3d(v, axis, starts, num_threads=n), want), (axis, n)


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.uint16])
def test_range_max_along_x_keeps_matlabs_rules(dtype):
    """long runs along x (the vectorized fold): rows of NaN only, of -Inf, of negatives with
    zeros of both signs, a NaN first, and plain rows, in runs of several lengths"""
    rng = CounterRNG(9, f"xrows/{np.dtype(dtype).name}")
    v = (rng.uniform((3, 8, 301)) * 100).astype(dtype)
    if np.dtype(dtype).kind == "f":
        v -= 50
        v[0, 0] = NAN
        v[0, 1] = -np.inf
        v[0, 2] = -np.abs(v[0, 2]) - 1
        v[0, 2, 41] = -0.0  # the first zero; a later one in a lower lane
        v[0, 2, 72] = 0.0
        v[0, 3] = v[0, 2]
        v[0, 3, 41] = 0.0
        v[0, 3, 72] = -0.0
        v[0, 4, 0] = NAN
        v[0, 5] = NAN
        v[0, 5, 77] = -np.inf
        v[1, 6, ::2] = -0.0
        v[1, 6, 1::2] = NAN
    for starts in ([0], [0, 33], [0, 40, 41, 200], [0, 1, 2, 3, 100, 300]):
        assert same_bits(mip.range_max_3d(v, -1, starts), range_max(v, 2, starts)), starts


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_min_bbox_keeps_the_first_zero_in_matlab_order(dtype, monkeypatch):
    """a zero minimum is the box's first zero (-0 or 0) in MATLAB's order (y fastest), at every
    thread count and in every lane of the vectorized loop"""
    monkeypatch.setattr(mip, "_threads", lambda n: 1 if n is None else n)
    rng = CounterRNG(7, "zeros")
    for t in range(16):
        v = (rng.uniform((5, 9, 37)) + 0.25).astype(dtype)
        flat = v.reshape(-1)
        pos = rng.integers(0, flat.size - 1, 6, dtype=np.int64).ravel()
        flat[pos[:3]] = -0.0
        flat[pos[3:]] = 0.0
        flat[pos[0] if t % 2 else pos[3]] = NAN
        if t >= 12:
            # a zero of each sign in one plane: MATLAB's order (y fastest) reaches x = 3 before
            # x = 10, C order reaches y = 4 before y = 6
            v[:] = 1
            v[2, 6, 3] = -0.0 if t % 2 else 0.0
            v[2, 4, 10] = 0.0 if t % 2 else -0.0
        box = (0, t % 2, t % 3, 5 - t % 2, 9, 37 - t)
        want = serial_min(v, box)
        for n in (1, 2, 3, 5, 16):
            assert same_bits(mip.min_bbox_3d(v, box, num_threads=n), want), (t, n)


def test_max_pooling_matches_reference_on_odd_shapes():
    rng = CounterRNG(5, "odd")
    for t in range(30):
        shape = tuple(int(x) for x in rng.integers(1, 13, 3).ravel())
        pool = tuple(int(x) for x in rng.integers(1, 6, 3).ravel())
        v = (
            native_specials("single", shape)
            if min(shape) > 8
            else native(shape, "single", f"odd{t}")
        )
        assert same_bits(mip.max_pooling_3d(v, pool), block_max(v, pool)), (shape, pool)


def test_other_integer_classes_pool_as_matlabs_code():
    v = native((7, 11, 9), "int16", "fallback")
    out = mip.max_pooling_3d(v, (2, 3, 4))
    m = [-(-s // p) for s, p in zip(v.shape, (2, 3, 4), strict=True)]
    want = np.zeros(m, np.int16)
    for k in range(m[0]):
        for i in range(m[1]):
            for j in range(m[2]):
                blk = v[2 * k : 2 * k + 2, 3 * i : 3 * i + 3, 4 * j : 4 * j + 4]
                full = blk.shape == (2, 3, 4)
                want[k, i, j] = blk.max() if full else max(blk.max(), 0)
    assert np.array_equal(out, want)


def test_layout_and_input_are_respected():
    v = native((9, 15, 14), "uint16", "layout")
    keep = v.copy()
    f = np.asfortranarray(v)
    assert same_bits(mip.max_pooling_3d(f, (2, 2, 3)), mip.max_pooling_3d(v, (2, 2, 3)))
    assert same_bits(mip.range_max_3d(f, 1, [0, 3]), mip.range_max_3d(v, 1, [0, 3]))
    box = (1, 1, 1, 9, 9, 9)
    assert mip.min_bbox_3d(f, box) == mip.min_bbox_3d(v, box)
    sliced = native((9, 30, 14), "uint16", "layout")[:, ::2]
    assert same_bits(mip.max_pooling_3d(sliced, (2, 2, 2)), block_max(sliced.copy(), (2, 2, 2)))
    assert same_bits(v, keep)
    out = mip.max_pooling_3d(v, (2, 2, 2))
    assert out.flags.c_contiguous
    # an image stays an image; a 2D box works on it
    img = v[0]
    assert mip.max_pooling_3d(img, (1, 2, 2)).shape == (8, 7)
    assert mip.range_max_3d(img, -2, [0]).shape == (1, 14)
    assert mip.min_bbox_3d(img, (1, 2, 5, 9)) == img[1:5, 2:9].min()
    with pytest.raises(mip.MatlabError):
        mip.min_bbox_3d(v, (0, 0, 0, 10, 15, 14))
    assert mip.min_bbox_3d(v, (0, 0, 0, None, None, None)) == v.min()


def test_binding_refuses_other_layouts_and_classes():
    ext = mip._ext()
    vf = np.asfortranarray(np.zeros((4, 4, 4), np.uint16))
    out = np.zeros((2, 2, 2), np.uint16)
    with pytest.raises(TypeError):
        ext.max_pooling_3d(vf, out, 2, 2, 2, 1)
    with pytest.raises(TypeError):
        ext.max_pooling_3d(np.zeros((4, 4, 4), np.uint16), out.astype(np.float32), 2, 2, 2, 1)
    with pytest.raises(ValueError):
        ext.max_pooling_3d(np.zeros((4, 4, 4), np.uint16), out, 2, 2, 2, 0)


@pytest.mark.parametrize("dtype", [np.bool_, np.uint8, np.uint16, np.float32, np.float64])
def test_every_instruction_set_gives_the_baseline_bits(dtype):
    """the kernels run the best instruction set of the CPU (SSE4.1, AVX2 on x86-64); each one
    the CPU has gives the baseline's results, bit for bit"""
    ext = mip._ext()
    isas = list(range(ext.best_isa() + 1))
    shape = (13, 41, 133)
    if np.dtype(dtype).kind == "f":
        v = native_specials(np.dtype(dtype) == np.float32 and "single" or "double", shape)
        v[3, 7] = -0.0
        v[3, 7, 50] = 0.0
        v[4] = np.where(np.isnan(v[4]), v[4], -np.abs(v[4]))
        v[4, :, 9] = 0.0
        v[4, :, 3] = -0.0
    elif dtype == np.bool_:
        v = native(shape, "logical", "isa")
    else:
        v = native(shape, np.dtype(dtype).name, "isa")
    runs = []
    for isa in isas:
        r = {}
        for pool in ((1, 1, 1), (2, 3, 4), (1, 1, 7), (13, 41, 133)):
            m = [-(-s // p) for s, p in zip(shape, pool, strict=True)]
            out = np.empty(m, dtype=v.dtype)
            ext.max_pooling_3d(v, out, *pool, 3, isa)
            r["pool", pool] = out
        for box in ((0, 0, 0, 13, 41, 133), (2, 3, 4, 11, 40, 100), (4, 0, 0, 5, 41, 133)):
            r["min", box] = ext.min_box_3d(v, *box, 3, isa)
        if dtype != np.bool_:
            for axis in (0, 1, 2):
                for starts in ([0], list(range(0, shape[axis], 6))):
                    osz = list(shape)
                    osz[axis] = len(starts)
                    out = np.empty(osz, dtype=v.dtype)
                    ext.range_max_3d(v, out, axis, starts, 3, isa)
                    r["range", axis, len(starts)] = out
        runs.append(r)
    for isa, r in zip(isas[1:], runs[1:], strict=True):
        for k, val in r.items():
            assert same_bits(np.asarray(val), np.asarray(runs[0][k])), (isa, k)
    with pytest.raises(ValueError):
        ext.min_box_3d(v, 0, 0, 0, 1, 1, 1, 1, ext.best_isa() + 1)


def test_mex_uint64_conversion():
    u = mip._mex_u64
    assert u(2.9) == 2 and u(-0.5) == 0 and u(-1.0) == 2**64 - 1
    assert u(float("nan")) == 2**63 and u(float("inf")) == 0 and u(-float("inf")) == 2**63
    assert u(2.0**63) == 2**63 and u(2.0**64) == 0
