"""Crop and indexing kernels (petakit5d.kernels.crop).

Oracle cases (MATLAB: crop3d_mex, crop4d_mex, indexing3d_mex, indexing4d_mex,
indexing4d_crop_mex and the wrappers crop3d, crop4d, indexing4d) over the mex's classes
(logical, uint8, uint16, single, double), 1-D to 4-D arrays, boxes at every border, whole rows,
planes and volumes (the kernel's merged runs), empty ranges, non-integer bounds, arrays large
enough for several threads, complex data, the classes and boxes the mex rejects (the wrappers'
MATLAB code: indexing, assignment with conversion, growth) and the input errors. Inputs are
MATLAB's (y, x, z, t) arrays; the roles give the Python calls (t, z, y, x) arrays, and the
boxes stay MATLAB's values (the *_matlab functions take them). The indexing functions write
into their first argument in place: their cases return it.

Then tests/crop/CropMexTest.m and the crop and indexing methods of
tests/platform/MexInputCheckTest.m, and the kernels against NumPy at several thread counts,
layouts and classes.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.kernels import crop as c
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Vol
from petakit5d.testing.synth import CounterRNG, random_image

MOD = "petakit5d.kernels.crop"
CROP_INPUT = "crop:inputError"
CROP_TYPE = "crop:dataTypeError"
INDEXING_INPUT = "indexing:inputError"
TIFF_TYPE = "tiff:dataTypeError"
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
MEX_CLASSES = ("logical", "uint8", "uint16", "single", "double")
OTHER_CLASSES = ("int8", "int16", "uint32", "int32", "uint64", "int64")
IM_SIZE = (23, 19, 11)  # CropMexTest's ImSize (MATLAB's y, x, z)
IM_SIZE_4D = (13, 11, 7, 3)

# ---------------------------------------------------------------------------------------------
# input generators: MATLAB's arrays (the case keys are made from these bytes)
# ---------------------------------------------------------------------------------------------


def _python_size(msize):
    """our shape of MATLAB size msize"""
    m = list(msize)
    return tuple(m[2:][::-1]) + tuple(m[:2])


def image(msize, cls, name="im"):
    """CropMexTest's randomImage of MATLAB size msize (single and double: integral values
    times 1000 / 8, exact everywhere)"""
    rng = CounterRNG(11, f"crop/{name}/{cls}/{'x'.join(map(str, msize))}")
    shape = _python_size(msize)
    dt = np.dtype(DTYPES[cls])
    if dt.kind == "f":
        v = rng.integers(-8000, 8000, shape).astype(dt) / dt.type(8)
        return to_matlab(v.astype(dt))
    if dt.itemsize == 8:
        # 64-bit integers from two 32-bit words, beyond 2^53
        hi = rng.uint32(shape).astype(np.uint64) >> np.uint64(1 if dt.kind == "i" else 0)
        lo = rng.spawn("lo").uint32(shape).astype(np.uint64)
        return to_matlab(((hi << np.uint64(32)) | lo).astype(dt))
    return to_matlab(random_image(shape, cls, rng))


def zeros(msize, cls):
    return to_matlab(np.zeros(_python_size(msize), DTYPES[cls]))


def complex_image(msize, cls, name="cim"):
    re = image(msize, cls, name + "re")
    im = image(msize, cls, name + "im")
    return np.asfortranarray(re + 1j * im).astype(np.complex64 if cls == "single" else complex)


def _b(v):
    """a box for a case id"""
    return "_".join(str(x).replace(".", "p").replace("-", "m") for x in v)


# ---------------------------------------------------------------------------------------------
# oracle cases: crop3d_mex and crop3d
# ---------------------------------------------------------------------------------------------

BOXES_3D = {
    # CropMexTest's bbox3d
    "interior": [3, 5, 2, 20, 17, 9],
    "full": [1, 1, 1, 23, 19, 11],
    "singleVoxel": [7, 8, 9, 7, 8, 9],
    "singleSlice": [1, 1, 11, 23, 19, 11],
    "lastRow": [23, 1, 1, 23, 19, 11],
    # borders, whole Python rows (x) and planes (y and x), corners
    "firstColumn": [1, 1, 1, 23, 1, 11],
    "lastColumn": [1, 19, 1, 23, 19, 11],
    "firstPlane": [1, 1, 1, 23, 19, 1],
    "corner": [23, 19, 11, 23, 19, 11],
    "wholeRows": [5, 1, 2, 9, 19, 10],
    "wholePlanes": [1, 1, 3, 23, 19, 8],
    "lowCorner": [1, 1, 1, 4, 3, 2],
    # empty ranges (end = start - 1) and non-integer bounds (truncated)
    "emptyY": [4, 1, 1, 3, 19, 11],
    "emptyX": [1, 20, 1, 23, 19, 11],
    "emptyZ": [1, 1, 6, 23, 19, 5],
    "fraction": [1.5, 1, 2, 6.9, 5, 4],
    "fraction2": [2.999, 3.5, 1.01, 2.5, 19.99, 11.9],
}
SHAPES_3D = {
    "image": ((23, 19), {"in": [3, 5, 1, 20, 17, 1], "full": [1, 1, 1, 23, 19, 1]}),
    "row": ((1, 40), {"in": [1, 3, 1, 1, 38, 1], "full": [1, 1, 1, 1, 40, 1]}),
    "column": ((40, 1), {"in": [3, 1, 1, 38, 1, 1], "full": [1, 1, 1, 40, 1, 1]}),
    # crop3d_mex reads the first volume of a 4-D array
    "first_volume": (IM_SIZE_4D, {"in": [2, 3, 2, 12, 9, 6], "full": [1, 1, 1, 13, 11, 7]}),
    "plane1": ((13, 11, 1, 2), {"in": [2, 3, 1, 12, 9, 1]}),
    # several threads (over 1 MB)
    "big": ((301, 257, 17), {"in": [2, 3, 1, 300, 250, 17], "full": [1, 1, 1, 301, 257, 17]}),
    "longrow": ((3, 90001, 2), {"in": [1, 2, 1, 3, 90000, 2], "rows": [2, 1, 1, 3, 90001, 2]}),
}


def _crop_case(cid, fn, data, bbox, **kw):
    py = {"crop3d_mex": "crop3d_matlab", "crop3d": "crop3d_matlab"}.get(fn)
    py = py or ("crop4d_matlab" if fn.startswith("crop4d") else None)
    kwargs = {"mex_only": True} if fn.endswith("_mex") else {}
    return Case(
        cid,
        matlab=fn,
        args=[Vol(In("im")), bbox],
        python=f"{MOD}:{py}",
        py_kwargs=kwargs,
        out_roles=("vol",),
        inputs={"im": Value(data)},
        **kw,
    )


def _err(cid, fn, data, bbox, ident, **kw):
    return _crop_case(
        cid, fn, data, bbox, expect_error=ident, py_error=c.MatlabError, error_message=True, **kw
    )


def _crop3d_cases():
    cases = []
    for cls in MEX_CLASSES:
        for bname, b in BOXES_3D.items():
            cases.append(
                _crop_case(
                    f"crop3d_mex.{cls}.{bname}", "crop3d_mex", lambda c_=cls: image(IM_SIZE, c_), b
                )
            )
        for bname in ("interior", "full", "emptyY", "fraction"):
            cases.append(
                _crop_case(
                    f"crop3d.{cls}.{bname}",
                    "crop3d",
                    lambda c_=cls: image(IM_SIZE, c_),
                    BOXES_3D[bname],
                )  # fmt: skip
            )
        for sname, (msize, boxes) in SHAPES_3D.items():
            if sname in ("big", "longrow") and cls not in ("uint16", "logical", "double"):
                continue
            for bname, b in boxes.items():
                cases.append(
                    _crop_case(
                        f"crop3d_mex.{cls}.{sname}.{bname}",
                        "crop3d_mex",
                        lambda m=msize, c_=cls: image(m, c_),
                        b,
                    )
                )
    # the classes the mex rejects: crop3d's MATLAB code
    for cls in OTHER_CLASSES:
        data = lambda c_=cls: image(IM_SIZE, c_)  # noqa: E731
        cases.append(_crop_case(f"crop3d.{cls}.interior", "crop3d", data, BOXES_3D["interior"]))
        cases.append(
            _err(f"crop3d_mex.{cls}.interior", "crop3d_mex", data, BOXES_3D["interior"], CROP_TYPE)
        )
    # complex data: the mex's API reads the real part
    for cls in ("single", "double"):
        cases.append(
            _crop_case(
                f"crop3d_mex.complex_{cls}.interior",
                "crop3d_mex",
                lambda c_=cls: complex_image(IM_SIZE, c_),
                BOXES_3D["interior"],
            )
        )
    # boxes the mex rejects that crop3d's MATLAB code takes or refuses as MATLAB does
    fallback = {
        "zero_start": ([0, 1, 1, 5, 5, 5], True),
        "half_start": ([0.5, 1, 1, 5, 5, 5], False),
        "reversed": ([4, 1, 1, 2, 2, 2], False),
        "beyond": ([1, 1, 1, 24, 19, 11], True),
        "five": ([1, 1, 1, 5, 5], True),
        "seven": ([1, 1, 1, 5, 5, 5, 9], False),
        "column": (np.array([[2.0], [3], [1], [9], [10], [4]]), False),
        "single": (np.array([2, 3, 1, 9, 10, 4], np.float32), False),
        "int32": (np.array([2, 3, 1, 9, 10, 4], np.int32), False),
        "nan_start": ([math.nan, 1, 1, 2, 2, 2], True),
        "inf_end": ([1, 1, 1, math.inf, 2, 2], True),
        "two_rows": (np.array([[2.0, 1, 9, 4, 10, 5], [3, 1, 10, 5, 11, 6]]), False),
    }
    for cls in ("uint16", "int16"):
        for bname, (b, fails) in fallback.items():
            cases.append(
                _crop_case(
                    f"crop3d.{cls}.bad_{bname}",
                    "crop3d",
                    lambda c_=cls: image(IM_SIZE, c_),
                    b,
                    expect_error=fails,
                    py_error=c.MatlabError,
                )
            )
    # the mex's checks (MexInputCheckTest.crop3dChecksItsInputs, CropMexTest.crop3dRejects...)
    im = lambda: image((6, 5, 4), "uint16", "checks")  # noqa: E731
    big = lambda: image(IM_SIZE, "uint16")  # noqa: E731
    errors = {
        "five_values": (big, [1, 1, 1, 5, 5]),
        "zero_start": (big, [0, 1, 1, 5, 5, 5]),
        "beyond_z": (big, [1, 1, 1, 23, 19, 12]),
        "single_bbox": (im, np.array([1, 1, 1, 2, 2, 2], np.float32)),
        "int32_bbox": (im, np.array([1, 1, 1, 2, 2, 2], np.int32)),
        "int64_bbox": (im, np.array([1, 1, 1, 2, 2, 2], np.float64).view(np.int64)),
        "empty_bbox": (im, np.zeros((0, 6))),
        "start_m1": (im, [-1, 1, 1, 2, 2, 2]),
        "nan_start": (im, [math.nan, 1, 1, 2, 2, 2]),
        "inf_end": (im, [1, 1, 1, math.inf, 2, 2]),
        "end_below": (im, [4, 1, 1, 2, 2, 2]),
        "nan_end": (im, [1, 1, 1, 2, math.nan, 2]),
        "seven": (im, [1, 1, 1, 2, 2, 2, 2]),
        "column": (im, np.array([[1.0], [1], [1], [2], [2], [2]])),
        "none": (im, None),
    }
    for ename, (data, b) in errors.items():
        cases.append(_err(f"crop3d_mex.error.{ename}", "crop3d_mex", data, b, CROP_INPUT))
    cases.append(
        Case(
            "crop3d_mex.error.sparse",
            code="out = {crop3d_mex(sparse(6, 5), [1, 1, 1, 2, 2, 1])};",
            python=lambda: c.crop3d_matlab(_sparse(6, 5), [1, 1, 1, 2, 2, 1], mex_only=True),
            expect_error=CROP_INPUT,
            py_error=c.MatlabError,
            error_message=True,
        )
    )
    return cases


def _sparse(m, n):
    from scipy import sparse

    return sparse.csc_matrix((m, n))


# ---------------------------------------------------------------------------------------------
# oracle cases: crop4d_mex and crop4d
# ---------------------------------------------------------------------------------------------

BOXES_4D = {
    # CropMexTest's boxes
    "interior": [2, 3, 2, 1, 12, 9, 6, 3],
    "full": [1, 1, 1, 1, 13, 11, 7, 3],
    "voxel": [5, 5, 5, 2, 5, 5, 5, 2],
    # whole Python rows, planes and volumes, the last time point, borders
    "wholeRows": [2, 1, 2, 1, 12, 11, 6, 3],
    "wholePlanes": [1, 1, 2, 1, 13, 11, 6, 3],
    "wholeVolumes": [1, 1, 1, 2, 13, 11, 7, 3],
    "lastTime": [1, 1, 1, 3, 13, 11, 7, 3],
    "corner": [13, 11, 7, 3, 13, 11, 7, 3],
    "lastColumn": [1, 11, 1, 1, 13, 11, 7, 2],
    # non-integer bounds (MexInputCheckTest)
    "fraction": [2, 1, 1, 2, 3.5, 3, 2, 2],
    "fraction2": [1.9, 2.2, 1, 1, 4.99, 3.01, 7.5, 3.9],
}


def _crop4d_cases():
    cases = []
    for cls in MEX_CLASSES:
        for bname, b in BOXES_4D.items():
            cases.append(
                _crop_case(
                    f"crop4d_mex.{cls}.{bname}",
                    "crop4d_mex",
                    lambda c_=cls: image(IM_SIZE_4D, c_),
                    b,
                )  # fmt: skip
            )
        for bname in ("interior", "full", "fraction"):
            cases.append(
                _crop_case(
                    f"crop4d.{cls}.{bname}",
                    "crop4d",
                    lambda c_=cls: image(IM_SIZE_4D, c_),
                    BOXES_4D[bname],
                )  # fmt: skip
            )
        # 3-D and 2-D data, and many threads
        cases.append(
            _crop_case(
                f"crop4d_mex.{cls}.volume",
                "crop4d_mex",
                lambda c_=cls: image(IM_SIZE, c_),
                [3, 5, 2, 1, 20, 17, 9, 1],
            )  # fmt: skip
        )
        cases.append(
            _crop_case(
                f"crop4d_mex.{cls}.image",
                "crop4d_mex",
                lambda c_=cls: image((23, 19), c_),
                [3, 5, 1, 1, 20, 17, 1, 1],
            )  # fmt: skip
        )
        if cls in ("uint16", "double"):
            for bname, b in {
                "in": [2, 3, 2, 1, 100, 95, 33, 3],
                "volumes": [1, 1, 1, 2, 101, 97, 33, 3],
                "planes": [1, 1, 5, 1, 101, 97, 30, 3],
            }.items():
                cases.append(
                    _crop_case(
                        f"crop4d_mex.{cls}.big.{bname}",
                        "crop4d_mex",
                        lambda c_=cls: image((101, 97, 33, 3), c_),
                        b,
                    )  # fmt: skip
                )
    for cls in OTHER_CLASSES:
        data = lambda c_=cls: image(IM_SIZE_4D, c_)  # noqa: E731
        cases.append(_crop_case(f"crop4d.{cls}.interior", "crop4d", data, BOXES_4D["interior"]))
        cases.append(
            _err(f"crop4d_mex.{cls}.interior", "crop4d_mex", data, BOXES_4D["interior"], TIFF_TYPE)
        )
    for cls in ("single", "double"):
        cases.append(
            _crop_case(
                f"crop4d_mex.complex_{cls}.interior",
                "crop4d_mex",
                lambda c_=cls: complex_image(IM_SIZE_4D, c_),
                BOXES_4D["interior"],
            )
        )
    fallback = {
        "reversed": ([3, 1, 1, 1, 2, 2, 2, 1], False),
        "zero_start": ([0, 1, 1, 1, 2, 2, 2, 1], True),
        "beyond": ([1, 1, 1, 1, 14, 11, 7, 3], True),
        "seven": ([1, 1, 1, 1, 2, 2, 2], True),
        "nine": ([1, 1, 1, 1, 2, 2, 2, 2, 5], False),
        "single": (np.array([2, 3, 2, 1, 12, 9, 6, 3], np.float32), False),
        "half_start": ([0.5, 1, 1, 1, 2, 2, 2, 2], False),
    }
    for cls in ("uint16", "int16"):
        for bname, (b, fails) in fallback.items():
            cases.append(
                _crop_case(
                    f"crop4d.{cls}.bad_{bname}",
                    "crop4d",
                    lambda c_=cls: image(IM_SIZE_4D, c_),
                    b,
                    expect_error=fails,
                    py_error=c.MatlabError,
                )
            )
    # the mex's checks (MexInputCheckTest.crop4dChecksItsInputs, CropMexTest.crop4dRejects...)
    im = lambda: image((4, 3, 2, 2), "uint16", "checks4")  # noqa: E731
    big = lambda: image(IM_SIZE_4D, "uint16")  # noqa: E731
    errors = {
        "seven_values": (big, [1, 1, 1, 1, 2, 2, 2]),
        "beyond_y": (big, [1, 1, 1, 1, 14, 11, 7, 3]),
        "reversed": (big, [3, 1, 1, 1, 2, 2, 2, 1]),
        "five_d": (lambda: zeros((2, 2, 2, 2, 2), "uint16"), [1, 1, 1, 1, 2, 2, 2, 2]),
        "single_bbox": (im, np.array([1, 1, 1, 1, 2, 2, 2, 2], np.float32)),
        "int64_bbox": (im, np.array([1, 1, 1, 1, 2, 2, 2, 2], np.float64).view(np.int64)),
        "empty_bbox": (im, np.zeros((0, 8))),
        "nan_start": (im, [math.nan, 1, 1, 1, 2, 2, 2, 2]),
        "nan_end": (im, [1, 1, 1, 1, 2, math.nan, 2, 2]),
        "inf_end": (im, [1, 1, 1, 1, 2, 2, math.inf, 2]),
        "zero_start": (im, [1, 1, 0, 1, 2, 2, 2, 2]),
        "empty_range": (im, [2, 1, 1, 1, 1, 3, 2, 2]),
        "fraction_empty": (im, [3.5, 1, 1, 1, 2.9, 3, 2, 2]),
        "huge_start": (im, [1e300, 1, 1, 1, 2, 2, 2, 2]),
    }
    for ename, (data, b) in errors.items():
        cases.append(_err(f"crop4d_mex.error.{ename}", "crop4d_mex", data, b, CROP_INPUT))
    cases.append(
        Case(
            "crop4d_mex.error.sparse",
            code="out = {crop4d_mex(sparse(4, 3), [1, 1, 1, 1, 2, 2, 1, 1])};",
            python=lambda: c.crop4d_matlab(_sparse(4, 3), [1, 1, 1, 1, 2, 2, 1, 1], mex_only=True),
            expect_error=CROP_INPUT,
            py_error=c.MatlabError,
            error_message=True,
        )
    )
    return cases


# ---------------------------------------------------------------------------------------------
# oracle cases: indexing3d_mex, indexing4d_mex, indexing4d_crop_mex, indexing4d
# ---------------------------------------------------------------------------------------------

# the mex writes into A in place; the case returns it
_IDX3 = "A = in{1}; indexing3d_mex(A, in{2}, in{3}); out = {A};"
_IDX4 = "A = in{1}; indexing4d_mex(A, in{2}, in{3}); out = {A};"
_IDX4C3 = "A = in{1}; indexing4d_crop_mex(A, in{2}, in{3}); out = {A};"
_IDX4C4 = "A = in{1}; indexing4d_crop_mex(A, in{2}, in{3}, in{4}); out = {A};"


def _indexing4d_mex(out, bbox, region):
    """indexing4d_mex(A, bbox, R): the mex alone, in its argument order"""
    return c.indexing4d_matlab(out, region, bbox, mex_only=True)


def _idx_case(cid, code, py, a, b, r, rb=None, **kw):
    args = [Vol(In("A")), b, Vol(In("R"))] + ([] if rb is None else [rb])
    return Case(
        cid,
        code=code,
        args=args,
        python=_indexing4d_mex if py == "indexing4d_mex" else f"{MOD}:{py}",
        out_roles=("vol",),
        inputs={"A": Value(a), "R": Value(r)},
        **kw,
    )


def _idx_err(cid, code, py, a, b, r, ident, rb=None):
    return _idx_case(
        cid,
        code,
        py,
        a,
        b,
        r,
        rb,
        expect_error=ident,
        py_error=c.MatlabError,
        error_message=True,
        nargout=1,
    )


def _wrapper_case(cid, a, r, b, rb=None, **kw):
    """A = indexing4d(A, R, bbox[, rbbox])"""
    args = [Vol(In("A")), Vol(In("R")), b] + ([] if rb is None else [rb])
    return Case(
        cid,
        matlab="indexing4d",
        args=args,
        python=f"{MOD}:indexing4d_matlab",
        out_roles=("vol",),
        inputs={"A": Value(a), "R": Value(r)},
        **kw,
    )


def _size(b, n):
    """MATLAB's size of the region of box b (truncated as the mex truncates)"""
    return [int(b[i + n]) - (int(b[i]) - 1) for i in range(n)]


def _indexing3d_cases():
    cases = []
    boxes = {
        # CropMexTest.indexing3dAssignsRegionInPlace
        "interior": [3, 5, 2, 20, 17, 9],
        "full": [1, 1, 1, 23, 19, 11],
        "wholeRows": [5, 1, 2, 9, 19, 10],
        "wholePlanes": [1, 1, 3, 23, 19, 8],
        "corner": [23, 19, 11, 23, 19, 11],
        "lastRow": [23, 1, 1, 23, 19, 11],
        "emptyY": [4, 1, 1, 3, 19, 11],
        "fraction": [2.5, 1, 1, 3.5, 2, 2],
    }
    for cls in MEX_CLASSES:
        for bname, b in boxes.items():
            rsize = _size(b, 3)
            cases.append(
                _idx_case(
                    f"indexing3d_mex.{cls}.{bname}",
                    _IDX3,
                    "indexing3d_matlab",
                    lambda c_=cls: image(IM_SIZE, c_, "A3"),
                    b,
                    lambda c_=cls, s=rsize: image(s, c_, "R3"),
                )
            )
        cases.append(
            _idx_case(
                f"indexing3d_mex.{cls}.image",
                _IDX3,
                "indexing3d_matlab",
                lambda c_=cls: zeros((23, 19), c_),
                [3, 5, 1, 20, 17, 1],
                lambda c_=cls: image((18, 13), c_, "R2"),
            )
        )
        if cls in ("uint16", "single"):
            cases.append(
                _idx_case(
                    f"indexing3d_mex.{cls}.big",
                    _IDX3,
                    "indexing3d_matlab",
                    lambda c_=cls: zeros((301, 257, 17), c_),
                    [2, 3, 1, 300, 250, 17],
                    lambda c_=cls: image((299, 248, 17), c_, "Rbig"),
                )
            )
    for cls in ("single", "double"):
        cases.append(
            _idx_case(
                f"indexing3d_mex.complex_{cls}.interior",
                _IDX3,
                "indexing3d_matlab",
                lambda c_=cls: complex_image(IM_SIZE, c_, "cA"),
                [3, 5, 2, 20, 17, 9],
                lambda c_=cls: complex_image((18, 13, 8), c_, "cR"),
            )
        )
    # the mex's checks (MexInputCheckTest.indexing3dChecksItsInputs, CropMexTest)
    a = lambda: zeros((6, 5, 4), "uint16")  # noqa: E731
    r = lambda: to_matlab(np.ones((2, 2, 2), np.uint16))  # noqa: E731
    b = [1, 1, 1, 2, 2, 2]
    bi = [3, 5, 2, 20, 17, 9]
    big = lambda: zeros(IM_SIZE, "uint16")  # noqa: E731
    errors = {
        "four_d_data": (lambda: zeros((6, 5, 4, 2), "uint16"), b, r),
        "four_d_region": (a, b, lambda: zeros((2, 2, 2, 2), "uint16")),
        "single_bbox": (a, np.array(b, np.float32), r),
        "int64_bbox": (a, np.array(b, np.float64).view(np.int64), r),
        "empty_bbox": (a, np.zeros((0, 6)), r),
        "start_m1": (a, [-1, 1, 1, 0, 2, 2], r),
        "nan_start": (a, [math.nan, 1, 1, 2, 2, 2], r),
        "beyond": (a, [1, 1, 1, 2, 2, 5], r),
        "end_below": (a, [3, 1, 1, 1, 2, 2], r),
        "five_values": (a, [1, 1, 1, 2, 2], r),
        "mismatched_region": (big, bi, lambda: zeros((3, 3, 3), "uint16")),
        "mismatched_class": (big, bi, lambda: zeros(_size(bi, 3), "single")),
        "int16": (lambda: zeros((6, 5, 4), "int16"), b, lambda: zeros((2, 2, 2), "int16")),
    }
    for ename, (data, bb, reg) in errors.items():
        ident = TIFF_TYPE if ename == "int16" else INDEXING_INPUT
        cases.append(
            _idx_err(
                f"indexing3d_mex.error.{ename}", _IDX3, "indexing3d_matlab", data, bb, reg, ident
            )
        )
    for ename, code, py in (
        ("sparse_data", "indexing3d_mex(sparse(6, 5), [1, 1, 1, 2, 2, 1], sparse(2, 2));",
         lambda: c.indexing3d_matlab(_sparse(6, 5), [1, 1, 1, 2, 2, 1], _sparse(2, 2))),
        ("sparse_region", "indexing3d_mex(zeros(6, 5), [1, 1, 1, 2, 2, 1], sparse(2, 2));",
         lambda: c.indexing3d_matlab(np.zeros((6, 5)), [1, 1, 1, 2, 2, 1], _sparse(2, 2))),
    ):  # fmt: skip
        cases.append(
            Case(
                f"indexing3d_mex.error.{ename}",
                code=code + " out = {};",
                nargout=0,
                python=py,
                expect_error=INDEXING_INPUT,
                py_error=c.MatlabError,
                error_message=True,
            )
        )
    return cases


B4 = [2, 3, 2, 1, 12, 9, 6, 3]  # CropMexTest's indexing4dInputs


def _indexing4d_cases():
    cases = []
    boxes = {
        "interior": B4,
        "full": [1, 1, 1, 1, 13, 11, 7, 3],
        "wholeRows": [2, 1, 2, 1, 12, 11, 6, 3],
        "wholePlanes": [1, 1, 2, 1, 13, 11, 6, 3],
        "wholeVolumes": [1, 1, 1, 2, 13, 11, 7, 3],
        "corner": [13, 11, 7, 3, 13, 11, 7, 3],
        "fraction": [2.5, 1, 1, 1.2, 3.5, 2, 2, 2.9],
    }
    for cls in MEX_CLASSES:
        for bname, b in boxes.items():
            rsize = _size(b, 4)
            r = lambda c_=cls, s=rsize: image(s, c_, "R4")  # noqa: E731
            a = lambda c_=cls: image(IM_SIZE_4D, c_, "A4")  # noqa: E731
            cases.append(
                _idx_case(f"indexing4d_mex.{cls}.{bname}", _IDX4, "indexing4d_mex", a, b, r)
            )
            cases.append(_idx_case(f"indexing4d_crop_mex.{cls}.{bname}", _IDX4C3,
                                   "indexing4d_crop_matlab", a, b, r))  # fmt: skip
            if bname in ("interior", "full", "corner"):
                cases.append(_wrapper_case(f"indexing4d.{cls}.{bname}", a, r, b))
        # indexing4d_mex's empty ranges (indexing4d_crop_mex refuses them)
        cases.append(
            _idx_case(
                f"indexing4d_mex.{cls}.emptyZ",
                _IDX4,
                "indexing4d_mex",
                lambda c_=cls: image(IM_SIZE_4D, c_, "A4"),
                [1, 1, 4, 1, 13, 11, 3, 3],
                lambda c_=cls: to_matlab(np.zeros((3, 0, 13, 11), DTYPES[c_])),
            )
        )
        # CropMexTest.indexing4dCropAssignsSubregion and other region boxes
        rbig = lambda c_=cls: image((15, 15, 9, 3), c_, "Rc")  # noqa: E731
        for rname, rb in {
            "crop": [2, 3, 1, 1, 12, 9, 5, 3],
            "corner": [5, 7, 5, 1, 15, 13, 9, 3],
            "fraction": [2.9, 3.2, 1.5, 1, 12.7, 9.9, 5.01, 3.5],
        }.items():
            cases.append(
                _idx_case(
                    f"indexing4d_crop_mex.{cls}.rbox_{rname}",
                    _IDX4C4,
                    "indexing4d_crop_matlab",
                    lambda c_=cls: zeros(IM_SIZE_4D, c_),
                    B4,
                    rbig,
                    rb,
                )
            )
            cases.append(
                _wrapper_case(
                    f"indexing4d.{cls}.rbox_{rname}",
                    lambda c_=cls: zeros(IM_SIZE_4D, c_),
                    rbig,
                    B4,
                    rb,
                )
            )
        # whole volumes from a region box: the merged copies
        cases.append(
            _idx_case(
                f"indexing4d_crop_mex.{cls}.rbox_volumes",
                _IDX4C4,
                "indexing4d_crop_matlab",
                lambda c_=cls: zeros(IM_SIZE_4D, c_),
                [1, 1, 1, 2, 13, 11, 7, 3],
                lambda c_=cls: image((13, 11, 7, 5), c_, "Rv"),
                [1, 1, 1, 3, 13, 11, 7, 4],
            )
        )
    # many threads
    for cls in ("uint16", "double"):
        cases.append(
            _idx_case(
                f"indexing4d_crop_mex.{cls}.big",
                _IDX4C4,
                "indexing4d_crop_matlab",
                lambda c_=cls: zeros((101, 97, 33, 3), c_),
                [2, 1, 2, 1, 100, 97, 33, 3],
                lambda c_=cls: image((103, 99, 35, 3), c_, "Rb4"),
                [3, 2, 1, 1, 101, 98, 32, 3],
            )
        )
        cases.append(
            _idx_case(
                f"indexing4d_mex.{cls}.big",
                _IDX4,
                "indexing4d_mex",
                lambda c_=cls: zeros((101, 97, 33, 3), c_),
                [1, 1, 1, 1, 101, 97, 33, 3],
                lambda c_=cls: image((101, 97, 33, 3), c_, "Rf4"),
            )
        )
    cases += _indexing4d_fallback_cases()
    cases += _indexing4d_error_cases()
    return cases


def _indexing4d_fallback_cases():
    """indexing4d's MATLAB code after the mex's catch"""
    cases = []
    a16 = lambda: image(IM_SIZE_4D, "uint16", "A4")  # noqa: E731
    sz = _size(B4, 4)
    for cls in OTHER_CLASSES:
        cases.append(
            _wrapper_case(
                f"indexing4d.{cls}.interior",
                lambda c_=cls: image(IM_SIZE_4D, c_, "A4"),
                lambda c_=cls: image(sz, c_, "R4"),
                B4,
            )
        )
        cases.append(
            _wrapper_case(
                f"indexing4d.{cls}.rbox_crop",
                lambda c_=cls: zeros(IM_SIZE_4D, c_),
                lambda c_=cls: image((15, 15, 9, 3), c_, "Rc"),
                B4,
                [2, 3, 1, 1, 12, 9, 5, 3],
            )
        )

    def fractions(cls):
        """values with halves, negatives, NaN, Inf and beyond the integer ranges"""
        v = np.arange(math.prod(sz), dtype=np.float64).reshape(_python_size(sz)) * 37.25 - 900
        f = v.reshape(-1)
        f[:8] = [0.5, 1.5, -0.5, -1.5, 2.5, 65535.5, 65536, -70000]
        f[8:12] = [math.nan, math.inf, -math.inf, 1e20]
        return to_matlab(v.astype(DTYPES[cls]))

    # classes that differ: R converted to A's class
    for acls, rcls in (
        ("uint16", "double"),
        ("uint16", "single"),
        ("uint8", "double"),
        ("single", "double"),
        ("double", "single"),
        ("double", "uint16"),
        ("uint16", "int32"),
        ("uint16", "logical"),
        ("logical", "uint16"),
        ("int16", "double"),
    ):
        r = (lambda c_=rcls: fractions(c_)) if DTYPES[rcls] in (np.float32, np.float64) else (
            lambda c_=rcls: image(sz, c_, "Rx")
        )  # fmt: skip
        cases.append(
            _wrapper_case(
                f"indexing4d.{acls}_from_{rcls}",
                lambda c_=acls: image(IM_SIZE_4D, c_, "A4"),
                r,
                B4,
            )
        )
    # NaN into logical
    cases.append(
        _wrapper_case(
            "indexing4d.logical_from_nan",
            lambda: zeros(IM_SIZE_4D, "logical"),
            lambda: fractions("double"),
            B4,
            expect_error=True,
            py_error=c.MatlabError,
        )
    )
    # a scalar fills the box; sizes that match without singletons; sizes that do not
    cases.append(
        _wrapper_case("indexing4d.uint16.scalar", a16, lambda: np.array([[7]], np.uint16), B4)
    )
    cases.append(
        _wrapper_case(
            "indexing4d.uint16.squeezed",
            a16,
            lambda: image([11, 1, 5], "uint16", "Rsq"),
            [2, 3, 2, 1, 12, 3, 6, 1],
        )
    )
    cases.append(
        _wrapper_case(
            "indexing4d.uint16.permuted",
            a16,
            lambda: image([1, 11, 5], "uint16", "Rsq"),
            [2, 3, 2, 1, 12, 3, 6, 1],
        )
    )
    cases.append(
        _wrapper_case(
            "indexing4d.uint16.transposed",
            a16,
            lambda: image([5, 11], "uint16", "Rsq"),
            [2, 3, 2, 1, 12, 3, 6, 1],
            expect_error=True,
            py_error=c.MatlabError,
        )
    )
    cases.append(
        _wrapper_case(
            "indexing4d.uint16.mismatch",
            a16,
            lambda: image([3, 3, 3], "uint16", "Rm"),
            B4,
            expect_error=True,
            py_error=c.MatlabError,
        )
    )
    # boxes beyond A: MATLAB grows it
    for gname, b in (
        ("grow_y", [12, 3, 2, 1, 15, 9, 6, 3]),
        ("grow_t", [2, 3, 2, 3, 12, 9, 6, 4]),
        ("grow_xz", [2, 10, 6, 1, 12, 13, 9, 1]),
        ("grow_from_3d", [2, 3, 2, 1, 12, 9, 6, 2]),
    ):
        rsize = _size(b, 4)
        a = a16 if gname != "grow_from_3d" else (lambda: image(IM_SIZE, "uint16", "A3"))
        cases.append(
            _wrapper_case(
                f"indexing4d.uint16.{gname}",
                a,
                lambda s=rsize: image(s, "uint16", "Rg"),
                b,
            )
        )
    # other boxes the mex rejects
    for bname, b, fails in (
        ("zero_start", [0, 3, 2, 1, 10, 9, 6, 3], True),
        ("single", np.array(B4, np.float32), False),
        ("seven", B4[:7], True),
        ("half_start", [1.5, 3, 2, 1, 11.5, 9, 6, 3], False),
    ):
        cases.append(
            _wrapper_case(
                f"indexing4d.uint16.bad_{bname}",
                a16,
                lambda: image(sz, "uint16", "R4"),
                b,
                expect_error=fails,
                py_error=c.MatlabError,
            )
        )
    # a region box beyond R
    cases.append(
        _wrapper_case(
            "indexing4d.uint16.rbox_beyond",
            lambda: zeros(IM_SIZE_4D, "uint16"),
            lambda: image((15, 15, 9, 3), "uint16", "Rc"),
            B4,
            [6, 7, 5, 1, 16, 13, 9, 3],
            expect_error=True,
            py_error=c.MatlabError,
        )
    )
    return cases


def _indexing4d_error_cases():
    """MexInputCheckTest.indexing4dChecksItsInputs and the other checks of both mex functions"""
    cases = []
    a = lambda: zeros((4, 3, 2, 2), "uint16")  # noqa: E731
    r = lambda: to_matlab(np.ones((1, 1, 2, 2), np.uint16))  # noqa: E731
    b = [1, 1, 1, 1, 2, 2, 1, 1]
    common = {
        "five_d_data": (lambda: zeros((4, 3, 2, 2, 2), "uint16"), b, r),
        "five_d_region": (a, b, lambda: zeros((2, 2, 1, 1, 2), "uint16")),
        "single_bbox": (a, np.array(b, np.float32), r),
        "int64_bbox": (a, np.array(b, np.float64).view(np.int64), r),
        "empty_bbox": (a, np.zeros((0, 8)), r),
        "nan_start": (a, [math.nan, 1, 1, 1, 2, 2, 1, 1], r),
        "seven": (a, b[:7], r),
        "beyond": (a, [1, 1, 1, 1, 2, 2, 1, 3], r),
        "class": (a, b, lambda: zeros((2, 2), "double")),
        "int16": (lambda: zeros((4, 3, 2, 2), "int16"), b, lambda: zeros((2, 2), "int16")),
        "size": (a, b, lambda: zeros((2, 3), "uint16")),
    }
    for fn, code, py in (
        ("indexing4d_mex", _IDX4, "indexing4d_mex"),
        ("indexing4d_crop_mex", _IDX4C3, "indexing4d_crop_matlab"),
    ):
        for ename, (data, bb, reg) in common.items():
            ident = TIFF_TYPE if ename == "int16" else INDEXING_INPUT
            case = _idx_err(f"{fn}.error.{ename}", code, py, data, bb, reg, ident)
            cases.append(case)
    cases.append(
        _idx_err(
            "indexing4d_mex.error.start_m1",
            _IDX4,
            "indexing4d_mex",
            a,
            [-1, 1, 1, 1, 0, 2, 1, 1],
            r,
            INDEXING_INPUT,
        )
    )
    cases.append(
        _idx_err(
            "indexing4d_mex.error.end_below",
            _IDX4,
            "indexing4d_mex",
            a,
            [3, 1, 1, 1, 1, 2, 1, 1],
            r,
            INDEXING_INPUT,
        )
    )
    # indexing4d_crop_mex: empty ranges, and its region box
    rr = lambda: image((15, 15, 9, 3), "uint16", "Rc")  # noqa: E731
    big = lambda: zeros(IM_SIZE_4D, "uint16")  # noqa: E731
    cases.append(
        _idx_err(
            "indexing4d_crop_mex.error.empty_range",
            _IDX4C3,
            "indexing4d_crop_matlab",
            big,
            [1, 1, 4, 1, 13, 11, 3, 3],
            lambda: to_matlab(np.zeros((3, 0, 13, 11), np.uint16)),
            INDEXING_INPUT,
        )
    )
    rboxes = {
        "single": np.array(b, np.float32),
        "int64": np.array(b, np.float64).view(np.int64),
        "nan": [math.nan, 1, 1, 1, 2, 2, 1, 1],
        "seven": b[:7],
        "zero_start": [0, 1, 1, 1, 1, 2, 1, 1],
        "beyond": [1, 1, 1, 1, 16, 9, 5, 3],
        "reversed": [3, 1, 1, 1, 2, 9, 5, 3],
        "other_range": [1, 1, 1, 1, 11, 8, 5, 3],
        "nan_end": [1, 1, 1, 1, 11, math.nan, 5, 3],
    }
    for rname, rb in rboxes.items():
        cases.append(
            _idx_err(
                f"indexing4d_crop_mex.error.rbox_{rname}",
                _IDX4C4,
                "indexing4d_crop_matlab",
                big,
                B4,
                rr,
                INDEXING_INPUT,
                rb,
            )
        )
    for ename, code, py in (
        ("sparse_data",
         "indexing4d_mex(sparse(4, 3), [1, 1, 1, 1, 2, 2, 1, 1], sparse(2, 2));",
         lambda: c.indexing4d_matlab(_sparse(4, 3), _sparse(2, 2), [1, 1, 1, 1, 2, 2, 1, 1],
                                     mex_only=True)),
        ("crop_sparse_data",
         "indexing4d_crop_mex(sparse(4, 3), [1, 1, 1, 1, 2, 2, 1, 1], sparse(2, 2));",
         lambda: c.indexing4d_crop_matlab(_sparse(4, 3), [1, 1, 1, 1, 2, 2, 1, 1],
                                          _sparse(2, 2))),
    ):  # fmt: skip
        cases.append(
            Case(
                f"indexing4d.error.{ename}",
                code=code + " out = {};",
                nargout=0,
                python=py,
                expect_error=INDEXING_INPUT,
                py_error=c.MatlabError,
                error_message=True,
            )
        )
    return cases


CASES = _crop3d_cases() + _crop4d_cases() + _indexing3d_cases() + _indexing4d_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# tests/crop/CropMexTest.m, in intent
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("cls", MEX_CLASSES)
@pytest.mark.parametrize("bname", ["interior", "full", "singleVoxel", "singleSlice", "lastRow"])
def test_crop3d_matches_indexing(cls, bname):
    im = image(IM_SIZE, cls)
    b = BOXES_3D[bname]
    expected = im[b[0] - 1 : b[3], b[1] - 1 : b[4], b[2] - 1 : b[5]]
    v = from_matlab(im)
    for got in (c.crop3d_matlab(v, b, mex_only=True), c.crop3d_matlab(v, b)):
        assert got.dtype == expected.dtype
        np.testing.assert_array_equal(to_matlab(got).reshape(expected.shape), expected)
    pb = (b[2] - 1, b[0] - 1, b[1] - 1, b[5], b[3], b[4])
    np.testing.assert_array_equal(c.crop_3d(v, pb), v[pb[0] : pb[3], pb[1] : pb[4], pb[2] : pb[5]])


def test_crop3d_rejects_invalid_bbox():
    v = from_matlab(image(IM_SIZE, "uint16"))
    for b in ([1, 1, 1, 5, 5], [0, 1, 1, 5, 5, 5], [1, 1, 1, 23, 19, 12]):
        with pytest.raises(c.MatlabError) as e:
            c.crop3d_matlab(v, b, mex_only=True)
        assert e.value.identifier == CROP_INPUT
    with pytest.raises(c.MatlabError) as e:
        c.crop3d_matlab(v.astype(np.int16), [1, 1, 1, 5, 5, 5], mex_only=True)
    assert e.value.identifier == CROP_TYPE
    # the native boxes
    for b in (
        (0, 0, 0, 12, 23, 19),
        (0, 0, -1, 11, 23, 19),
        (0, 0, 3, 11, 23, 2),
        (0, 0, 0.5, 1, 1, 1),
    ):
        with pytest.raises(c.MatlabError) as e:
            c.crop_3d(v, b)
        assert e.value.identifier == CROP_INPUT


@pytest.mark.parametrize("cls", MEX_CLASSES)
def test_crop4d_matches_indexing(cls):
    im = image(IM_SIZE_4D, cls)
    v = from_matlab(im)
    for b in ([2, 3, 2, 1, 12, 9, 6, 3], [1, 1, 1, 1, *IM_SIZE_4D], [5, 5, 5, 2, 5, 5, 5, 2]):
        expected = im[b[0] - 1 : b[4], b[1] - 1 : b[5], b[2] - 1 : b[6], b[3] - 1 : b[7]]
        for got in (c.crop4d_matlab(v, b, mex_only=True), c.crop4d_matlab(v, b)):
            np.testing.assert_array_equal(to_matlab(got).reshape(expected.shape), expected)
        pb = [b[3] - 1, b[2] - 1, b[0] - 1, b[1] - 1, b[7], b[6], b[4], b[5]]
        np.testing.assert_array_equal(
            c.crop_4d(v, pb), v[pb[0] : pb[4], pb[1] : pb[5], pb[2] : pb[6], pb[3] : pb[7]]
        )


def test_crop4d_rejects_invalid_bbox():
    v = from_matlab(image(IM_SIZE_4D, "uint16"))
    for b in ([1, 1, 1, 1, 2, 2, 2], [1, 1, 1, 1, 14, 11, 7, 3], [3, 1, 1, 1, 2, 2, 2, 1]):
        with pytest.raises(c.MatlabError) as e:
            c.crop4d_matlab(v, b, mex_only=True)
        assert e.value.identifier == CROP_INPUT


@pytest.mark.parametrize("cls", MEX_CLASSES)
def test_indexing3d_assigns_region_in_place(cls):
    a = from_matlab(zeros(IM_SIZE, cls))
    b = [3, 5, 2, 20, 17, 9]
    r = image(_size(b, 3), cls)
    expected = zeros(IM_SIZE, cls)
    expected[b[0] - 1 : b[3], b[1] - 1 : b[4], b[2] - 1 : b[5]] = r
    assert c.indexing3d_matlab(a, b, from_matlab(r)) is a
    np.testing.assert_array_equal(to_matlab(a), expected)
    a2 = from_matlab(zeros(IM_SIZE, cls))
    assert c.indexing_3d(a2, (1, 2, 4, 9, 20, 17), from_matlab(r)) is a2
    np.testing.assert_array_equal(a2, a)


def test_indexing3d_rejects_mismatched_region():
    a = from_matlab(zeros(IM_SIZE, "uint16"))
    b = [3, 5, 2, 20, 17, 9]
    for r in (np.zeros((3, 3, 3), np.uint16), from_matlab(zeros(_size(b, 3), "single"))):
        with pytest.raises(c.MatlabError) as e:
            c.indexing3d_matlab(a, b, r)
        assert e.value.identifier == INDEXING_INPUT
    assert not a.any()


@pytest.mark.parametrize("cls", MEX_CLASSES)
def test_indexing4d_wrapper_matches_indexing(cls):
    a = from_matlab(zeros(IM_SIZE_4D, cls))
    r = image(_size(B4, 4), cls)
    expected = zeros(IM_SIZE_4D, cls)
    expected[1:12, 2:9, 1:6, 0:3] = r
    out = c.indexing4d_matlab(a, from_matlab(r), B4)
    assert out is a
    np.testing.assert_array_equal(to_matlab(out), expected)


def test_indexing4d_crop_assigns_subregion():
    a = from_matlab(zeros(IM_SIZE_4D, "uint16"))
    r = image((15, 15, 9, 3), "uint16")
    rb = [2, 3, 1, 1, *[B4[4 + i] - B4[i] + [2, 3, 1, 1][i] for i in range(4)]]
    c.indexing4d_crop_matlab(a, B4, from_matlab(r), rb)
    expected = zeros(IM_SIZE_4D, "uint16")
    expected[1:12, 2:9, 1:6, 0:3] = r[
        rb[0] - 1 : rb[4], rb[1] - 1 : rb[5], rb[2] - 1 : rb[6], rb[3] - 1 : rb[7]
    ]
    np.testing.assert_array_equal(to_matlab(a), expected)
    a2 = from_matlab(zeros(IM_SIZE_4D, "uint16"))
    c.indexing_4d(a2, (0, 1, 1, 2, 3, 6, 12, 9), from_matlab(r), (0, 0, 1, 2, 3, 5, 12, 9))
    np.testing.assert_array_equal(a2, a)


def test_indexing4d_mex_supports_logical():
    a = from_matlab(zeros(IM_SIZE_4D, "logical"))
    r = image(_size(B4, 4), "logical")
    c.indexing4d_matlab(a, from_matlab(r), B4, mex_only=True)
    expected = zeros(IM_SIZE_4D, "logical")
    expected[1:12, 2:9, 1:6, 0:3] = r
    np.testing.assert_array_equal(to_matlab(a), expected)


# ---------------------------------------------------------------------------------------------
# tests/platform/MexInputCheckTest.m: the crop and indexing methods, in intent
# ---------------------------------------------------------------------------------------------


def _as_int64(b):
    """an int64 array with the bits of double(b) (MexInputCheckTest.asInt64)"""
    return np.asarray(b, np.float64).view(np.int64)


def _raises(fn, ident):
    with pytest.raises(c.MatlabError) as e:
        fn()
    assert e.value.identifier == ident, (e.value.identifier, str(e.value))


def test_crop3d_checks_its_inputs():
    im = from_matlab(np.arange(1, 6 * 5 * 4 + 1, dtype=np.uint16).reshape(6, 5, 4, order="F"))
    for b in (
        np.array([1, 1, 1, 2, 2, 2], np.float32),
        np.array([1, 1, 1, 2, 2, 2], np.int32),
        _as_int64([1, 1, 1, 2, 2, 2]),
        np.zeros((0, 6)),
        [-1, 1, 1, 2, 2, 2],
        [math.nan, 1, 1, 2, 2, 2],
        [1, 1, 1, math.inf, 2, 2],
        [4, 1, 1, 2, 2, 2],
    ):
        _raises(lambda b=b: c.crop3d_matlab(im, b, mex_only=True), CROP_INPUT)
    _raises(lambda: c.crop3d_matlab(_sparse(6, 5), [1, 1, 1, 2, 2, 1], mex_only=True), CROP_INPUT)
    m = to_matlab(im)
    np.testing.assert_array_equal(c.crop3d_matlab(im, [4, 1, 1, 3, 5, 4]), from_matlab(m[3:3]))
    np.testing.assert_array_equal(
        to_matlab(c.crop3d_matlab(im, [1.5, 1, 2, 6.9, 5, 4])), m[0:6, :, 1:4]
    )


def test_crop4d_checks_its_inputs():
    im = from_matlab(
        np.arange(1, 4 * 3 * 2 * 2 + 1, dtype=np.uint16).reshape(4, 3, 2, 2, order="F")
    )
    _raises(
        lambda: c.crop4d_matlab(np.zeros((2,) * 5, np.uint16), [1, 1, 1, 1, 2, 2, 2, 2],
                                mex_only=True),
        CROP_INPUT,
    )  # fmt: skip
    _raises(
        lambda: c.crop4d_matlab(_sparse(4, 3), [1, 1, 1, 1, 2, 2, 1, 1], mex_only=True), CROP_INPUT
    )
    for b in (
        np.array([1, 1, 1, 1, 2, 2, 2, 2], np.float32),
        _as_int64([1, 1, 1, 1, 2, 2, 2, 2]),
        np.zeros((0, 8)),
        [math.nan, 1, 1, 1, 2, 2, 2, 2],
        [1, 1, 1, 1, 2, math.nan, 2, 2],
    ):
        _raises(lambda b=b: c.crop4d_matlab(im, b, mex_only=True), CROP_INPUT)
    got = c.crop4d_matlab(im, [2, 1, 1, 2, 3.5, 3, 2, 2], mex_only=True)
    np.testing.assert_array_equal(to_matlab(got), to_matlab(im)[1:3, :, :, 1])


def test_indexing3d_checks_its_inputs():
    a = np.zeros((4, 6, 5), np.uint16)
    r = np.ones((2, 2, 2), np.uint16)
    b = [1, 1, 1, 2, 2, 2]
    _raises(lambda: c.indexing3d_matlab(np.zeros((2, 4, 6, 5), np.uint16), b, r), INDEXING_INPUT)
    _raises(lambda: c.indexing3d_matlab(a, b, np.ones((2, 2, 2, 2), np.uint16)), INDEXING_INPUT)
    _raises(
        lambda: c.indexing3d_matlab(_sparse(6, 5), [1, 1, 1, 2, 2, 1], _sparse(2, 2)),
        INDEXING_INPUT,
    )
    _raises(
        lambda: c.indexing3d_matlab(np.zeros((6, 5)), [1, 1, 1, 2, 2, 1], _sparse(2, 2)),
        INDEXING_INPUT,
    )
    for bb in (np.array(b, np.float32), _as_int64(b), np.zeros((0, 6)), [-1, 1, 1, 0, 2, 2],
               [math.nan, 1, 1, 2, 2, 2]):  # fmt: skip
        _raises(lambda bb=bb: c.indexing3d_matlab(a, bb, r), INDEXING_INPUT)
    assert not a.any()
    c.indexing3d_matlab(a, [2.5, 1, 1, 3.5, 2, 2], r)
    expected = np.zeros((6, 5, 4), np.uint16, order="F")
    expected[1:3, 0:2, 0:2] = 1
    np.testing.assert_array_equal(to_matlab(a), expected)


@pytest.mark.parametrize("fn", ["indexing4d_mex", "indexing4d_crop_mex"])
def test_indexing4d_checks_its_inputs(fn):
    def call(a, b, r, rb=None):
        if fn == "indexing4d_mex":
            return c.indexing4d_matlab(a, r, b, mex_only=True)
        return c.indexing4d_crop_matlab(a, b, r, rb)

    a = np.zeros((2, 2, 4, 3), np.uint16)
    r = np.ones((1, 1, 2, 2), np.uint16)
    b = [1, 1, 1, 1, 2, 2, 1, 1]
    _raises(lambda: call(np.zeros((2, 2, 2, 4, 3), np.uint16), b, r), INDEXING_INPUT)
    _raises(lambda: call(a, b, np.ones((2, 1, 1, 2, 2), np.uint16)), INDEXING_INPUT)
    _raises(lambda: call(_sparse(4, 3), b, _sparse(2, 2)), INDEXING_INPUT)
    for bb in (
        np.array(b, np.float32),
        _as_int64(b),
        np.zeros((0, 8)),
        [math.nan, 1, 1, 1, 2, 2, 1, 1],
    ):
        _raises(lambda bb=bb: call(a, bb, r), INDEXING_INPUT)
    if fn == "indexing4d_mex":
        _raises(lambda: call(a, [-1, 1, 1, 1, 0, 2, 1, 1], r), INDEXING_INPUT)
    else:
        for rb in (np.array(b, np.float32), _as_int64(b), [math.nan, 1, 1, 1, 2, 2, 1, 1]):
            _raises(lambda rb=rb: call(a, b, r, rb), INDEXING_INPUT)
    assert not a.any()


# ---------------------------------------------------------------------------------------------
# the kernels against NumPy
# ---------------------------------------------------------------------------------------------

ALL_DTYPES = [
    np.bool_, np.uint8, np.int8, np.uint16, np.int16, np.uint32, np.int32, np.uint64, np.int64,
    np.float16, np.float32, np.float64, np.complex64, np.complex128, "M8[s]", "V3",
]  # fmt: skip


def _filled(shape, dtype, name="k"):
    dt = np.dtype(dtype)
    n = math.prod(shape)
    raw = CounterRNG(5, f"crop/kern/{name}/{shape}").integers(0, 255, n * dt.itemsize)
    return raw.astype(np.uint8).view(dt).reshape(shape)


def _same_bytes(a, b):
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def _boxes(shape):
    """boxes that cover the kernel's run merges: interior, whole x, whole (y, x), whole (z, y,
    x), everything, a single voxel, borders, empty"""
    n = len(shape)
    full = list(shape)
    out = [
        [0] * n + full,
        [1] * n + [s - 1 for s in full],
        [0] * n + [1] * n,
        [s - 1 for s in full] + full,
    ]
    for k in range(n):
        b = [1] * n + [s - 1 for s in full]
        for j in range(k, n):
            b[j], b[j + n] = 0, full[j]
        out.append(b)
        e = [0] * n + full
        e[k + n] = e[k]
        out.append(e)
    return [b for b in out if all(b[i] <= b[i + n] for i in range(n))]


def _np_crop(v, b):
    n = v.ndim
    return v[tuple(slice(b[i], b[i + n]) for i in range(n))]


@pytest.mark.parametrize("dtype", ALL_DTYPES, ids=str)
def test_kernel_copies_every_class(dtype):
    v = _filled((5, 7, 9), dtype)
    for b in _boxes(v.shape):
        got = c.crop_3d(v, b)
        assert got.flags.c_contiguous
        assert _same_bytes(got, np.ascontiguousarray(_np_crop(v, b))), b
    w = _filled((3, 4, 5, 6), dtype, "4")
    for b in _boxes(w.shape):
        assert _same_bytes(c.crop_4d(w, b), np.ascontiguousarray(_np_crop(w, b))), b
        out = np.zeros_like(w)
        reg = np.ascontiguousarray(_np_crop(w, b))
        c.indexing_4d(out, b, reg)
        want = np.zeros_like(w)
        want[tuple(slice(b[i], b[i + 4]) for i in range(4))] = reg
        assert _same_bytes(out, want), b


@pytest.mark.parametrize("threads", [1, 2, 3, 4, 7, 16, 64])
def test_results_do_not_depend_on_threads(threads, monkeypatch):
    monkeypatch.setattr(c, "_threads", lambda n: 1 if n is None else n)
    monkeypatch.setattr(c, "_SERIAL_BYTES", 0)
    monkeypatch.setattr(c, "_BYTES_PER_THREAD", 1)
    monkeypatch.setattr(c, "_PIECE", 7 if threads % 2 else 1 << 18)
    for shape, dtype in (((23, 97, 61), np.uint16), ((1, 300, 7), np.float64),
                         ((40, 2, 3), np.bool_), ((3, 5, 9001), np.float32)):  # fmt: skip
        v = _filled(shape, dtype, "thr")
        for b in _boxes(shape):
            want = np.ascontiguousarray(_np_crop(v, b))
            assert _same_bytes(c.crop_3d(v, b, num_threads=threads), want), (shape, b)
            out = np.zeros_like(v)
            c.indexing_3d(out, b, want, num_threads=threads)
            ref = np.zeros_like(v)
            ref[tuple(slice(b[i], b[i + 3]) for i in range(3))] = want
            assert _same_bytes(out, ref), (shape, b)


def test_big_copy_splits_runs_between_threads():
    """a box of one run (a whole volume) and one of few long runs, at many threads"""
    v = _filled((4, 512, 640), np.uint16, "big")
    for b in ((0, 0, 0, 4, 512, 640), (1, 0, 0, 3, 512, 640), (0, 5, 0, 4, 500, 640),
              (0, 5, 3, 4, 500, 637)):  # fmt: skip
        for t in (1, 2, 5, 24):
            assert _same_bytes(c.crop_3d(v, b, num_threads=t), np.ascontiguousarray(_np_crop(v, b)))


def test_layouts_and_open_ends():
    v = _filled((6, 9, 11), np.int16, "lay")
    b = (1, 2, 3, 5, 8, 10)
    want = np.ascontiguousarray(_np_crop(v, b))
    for src in (v, np.asfortranarray(v), v[:, :, :], np.flip(np.flip(v, 0), 0)):
        assert _same_bytes(c.crop_3d(src, b), want)
    s = _filled((6, 18, 11), np.int16, "lay")[:, ::2, :]
    assert _same_bytes(c.crop_3d(s, b), np.ascontiguousarray(_np_crop(s, b)))
    assert _same_bytes(
        c.crop_3d(v, (1, None, 3, None, 8, None)), np.ascontiguousarray(v[1:, :8, 3:])
    )
    img = v[0]
    assert _same_bytes(c.crop_3d(img, (2, 3, 8, 10)), np.ascontiguousarray(img[2:8, 3:10]))
    assert _same_bytes(c.crop_3d(v, np.array(b, np.int64)), want)
    assert _same_bytes(c.crop_3d(v, [float(x) for x in b]), want)


def test_results_go_into_out():
    v = _filled((6, 9, 11), np.float32, "out")
    b = (1, 2, 3, 5, 8, 10)
    out = np.full((4, 6, 7), -1, np.float32)
    assert c.crop_3d(v, b, out=out) is out
    assert _same_bytes(out, np.ascontiguousarray(_np_crop(v, b)))
    for bad in (np.empty((4, 6, 8), np.float32), np.empty((4, 6, 7), np.float64),
                np.empty((7, 6, 4), np.float32).T):  # fmt: skip
        with pytest.raises(ValueError):
            c.crop_3d(v, b, out=bad)
    with pytest.raises(ValueError):
        c.crop_3d(v, (0, 0, 0, 1, 1, 1), out=v[:1, :1, :1])


def test_indexing_writes_only_out():
    out = np.zeros((6, 9, 11), np.uint8)
    region = _filled((4, 6, 7), np.uint8, "reg")
    keep = region.copy()
    got = c.indexing_3d(out, (1, 2, 3, 5, 8, 10), np.asfortranarray(region))
    assert got is out
    assert _same_bytes(region, keep)
    assert _same_bytes(out[1:5, 2:8, 3:10], region)
    assert out.sum() == region.astype(np.int64).sum()
    # out is never copied: other layouts and read-only arrays are refused
    for bad in (
        np.asfortranarray(np.zeros((6, 9, 11), np.uint8)),
        np.zeros((6, 9, 22), np.uint8)[:, :, ::2],
    ):
        with pytest.raises(ValueError):
            c.indexing_3d(bad, (1, 2, 3, 5, 8, 10), region)
    ro = np.zeros((6, 9, 11), np.uint8)
    ro.flags.writeable = False
    with pytest.raises(ValueError):
        c.indexing_3d(ro, (1, 2, 3, 5, 8, 10), region)
    with pytest.raises(ValueError):
        c.indexing3d_matlab(np.asfortranarray(np.zeros((6, 9, 11), np.uint8)), [1, 1, 1, 2, 2, 2],
                            np.zeros((2, 2, 2), np.uint8))  # fmt: skip


def test_indexing_with_a_region_inside_out(monkeypatch):
    v = _filled((4, 10, 12), np.uint16, "alias")
    want = v.copy()
    want[0:3, 2:10, 1:12] = v[1:4, 0:8, 0:11].copy()
    c.indexing_3d(v, (0, 2, 1, 3, 10, 12), v[1:4, 0:8, 0:11])
    assert _same_bytes(v, want)
    # the box after the region in memory: a copy in place would read what it wrote
    v = _filled((4, 10, 12), np.uint16, "alias")
    want = v.copy()
    want[1:4, 0:8, 0:11] = v[0:3, 2:10, 1:12].copy()
    c.indexing_3d(v, (1, 0, 0, 4, 8, 11), v[0:3, 2:10, 1:12])
    assert _same_bytes(v, want)
    # a C-contiguous region (whole planes) that the kernel itself would copy over itself, in
    # small pieces (glibc's memcpy of one call is a memmove, so one piece would hide it)
    monkeypatch.setattr(c, "_PIECE", 2)
    v = _filled((5, 10, 12), np.uint16, "alias")
    want = v.copy()
    want[1:5] = v[0:4].copy()
    c.indexing_3d(v, (1, 0, 0, 5, 10, 12), v[0:4])
    assert _same_bytes(v, want)
    w = _filled((2, 4, 10, 12), np.uint16, "alias4")
    want = w.copy()
    want[:, 0:3, 2:10, 1:12] = w[:, 1:4, 0:8, 0:11].copy()
    c.indexing_4d(w, (0, 0, 2, 1, 2, 3, 10, 12), w, (0, 1, 0, 0, 2, 4, 8, 11))
    assert _same_bytes(w, want)


def test_native_checks():
    v = np.zeros((3, 4, 5), np.uint16)
    _raises(lambda: c.crop_3d(v, (0, 0, 0, 3, 4)), CROP_INPUT)
    _raises(lambda: c.crop_3d(v, (0, 0, 0, 3, 4, 6)), CROP_INPUT)
    _raises(lambda: c.crop_3d(v, (2, 0, 0, 1, 4, 5)), CROP_INPUT)
    _raises(lambda: c.crop_3d(v, (True, 0, 0, 1, 4, 5)), CROP_INPUT)
    _raises(lambda: c.crop_3d(v[None], (0,) * 4 + (1, 3, 4, 5)), CROP_INPUT)
    _raises(lambda: c.crop_4d(v, (0,) * 6), CROP_INPUT)
    _raises(
        lambda: c.indexing_3d(v, (0, 0, 0, 1, 1, 1), np.zeros((1, 1, 1), np.int16)), INDEXING_INPUT
    )
    _raises(
        lambda: c.indexing_3d(v, (0, 0, 0, 1, 1, 1), np.zeros((1, 1, 2), np.uint16)), INDEXING_INPUT
    )
    _raises(
        lambda: c.indexing_3d(v, (0, 0, 0, 1, 1, 1), np.zeros((1, 1), np.uint16)), INDEXING_INPUT
    )
    _raises(
        lambda: c.indexing_4d(v[None], (0,) * 4 + (1,) * 4, np.zeros((1, 2, 2, 2), np.uint16),
                              (0,) * 4 + (1, 3, 1, 1)),
        INDEXING_INPUT,
    )  # fmt: skip
    with pytest.raises(TypeError):
        c.crop_3d(np.array([[[None]]], dtype=object), (0, 0, 0, 1, 1, 1))
    with pytest.raises(TypeError):
        c.crop_3d([[[1]]], (0, 0, 0, 1, 1, 1))
    # empty boxes give empty arrays and write nothing
    assert c.crop_3d(v, (1, 1, 1, 1, 4, 5)).shape == (0, 3, 4)
    assert c.indexing_3d(v, (1, 1, 1, 1, 4, 5), np.zeros((0, 3, 4), np.uint16)) is v


def test_binding_refuses_bad_sizes():
    ext = c._ext()
    if ext is None:
        pytest.skip("the extension is not built")
    src = np.zeros(4 * 5 * 6, np.uint8)
    dst = np.zeros(2 * 2 * 2, np.uint8)
    ok = {
        "src_shape": [1, 4, 5, 6],
        "src_origin": [0, 0, 0, 0],
        "dst_shape": [1, 2, 2, 2],
        "dst_origin": [0, 0, 0, 0],
        "extent": [1, 2, 2, 2],
        "itemsize": 1,
        "num_threads": 1,
    }
    ext.copy_box(src, dst=dst, **ok)
    bad = [
        ("extent", [1, 2, 2, 3]),
        ("src_origin", [0, 3, 0, 0]),
        ("src_shape", [1, 4, 5, 7]),
        ("dst_shape", [1, 2, 2, 3]),
        ("num_threads", 0),
        ("itemsize", 0),
        ("extent", [1, 2, 2]),
    ]
    for key, val in bad:
        with pytest.raises(ValueError):
            ext.copy_box(src, dst=dst, **{**ok, key: val})
    with pytest.raises(TypeError):
        ext.copy_box(src.astype(np.uint16), dst=dst, **ok)
    ro = dst.copy()
    ro.flags.writeable = False
    with pytest.raises(TypeError):
        ext.copy_box(src, dst=ro, **ok)


def test_numpy_copy_without_the_extension(monkeypatch):
    v = _filled((5, 7, 9), np.float64, "noext")
    b = (1, 2, 3, 4, 6, 8)
    want = c.crop_3d(v, b)
    monkeypatch.setattr(c, "_ext", lambda: None)
    assert _same_bytes(c.crop_3d(v, b), want)
    out = np.zeros_like(v)
    c.indexing_3d(out, b, want)
    assert _same_bytes(out[1:4, 2:6, 3:8], want)


def test_matlab_layer_keeps_the_mex_quirks():
    v = from_matlab(image(IM_SIZE_4D, "uint16"))
    # crop3d_mex crops the first volume of a 4-D array; crop3d's MATLAB code would read z and t
    got = c.crop3d_matlab(v, [1, 1, 1, 13, 11, 7], mex_only=True)
    np.testing.assert_array_equal(got, v[0])
    # a box of two rows: the first 6 values in MATLAB's (column) order
    b2 = np.array([[2.0, 1, 9, 4, 10, 5], [3, 1, 10, 5, 11, 6]])
    np.testing.assert_array_equal(
        c.crop3d_matlab(from_matlab(image(IM_SIZE, "uint16")), b2),
        c.crop3d_matlab(from_matlab(image(IM_SIZE, "uint16")), [2, 3, 1, 1, 9, 10]),
    )
    # complex data: the real part
    z = from_matlab(complex_image(IM_SIZE, "double"))
    got = c.crop3d_matlab(z, [3, 5, 2, 20, 17, 9], mex_only=True)
    assert got.dtype == np.float64
    np.testing.assert_array_equal(got, z.real[1:9, 2:20, 4:17])


def test_wrapper_fallback_converts_and_grows():
    a = from_matlab(zeros((4, 3, 2, 2), "uint16"))
    r = np.array([[0.5, 1.5], [-3.0, 70000.0]])  # (y, x): MATLAB 2 x 2
    out = c.indexing4d_matlab(a, r, [1, 1, 1, 1, 2, 2, 1, 1])
    assert out is a
    assert to_matlab(a)[:2, :2, 0, 0].tolist() == [[1, 2], [0, 65535]]
    g = c.indexing4d_matlab(a, np.ones((1, 1), np.uint16), [5, 1, 1, 3, 5, 1, 1, 3])
    assert g is not a and g.shape == (3, 2, 5, 3) and to_matlab(g)[4, 0, 0, 2] == 1
    assert to_matlab(g)[:4, :3, :2, :2].tolist() == to_matlab(a).tolist()
    with pytest.raises(c.MatlabError) as e:
        c.indexing4d_matlab(from_matlab(zeros((4, 3, 2, 2), "logical")), np.array([[math.nan]]),
                            [1, 1, 1, 1, 1, 1, 1, 1])  # fmt: skip
    assert e.value.identifier == "MATLAB:nologicalnan"


def test_import_is_light():
    import subprocess
    import sys

    code = (
        "import sys; import petakit5d.kernels.crop; "
        "print(sorted(m for m in ('scipy', 'skimage', 'cupy') if m in sys.modules))"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "[]"
