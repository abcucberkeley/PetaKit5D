"""MAT-file reading and writing (petakit5d.io.matfile).

Every case is a MATLAB expression, the value read_mat returns for it, and the other Python values
write_mat maps to the same MATLAB value. Against MATLAB (the oracle, cluster only): MATLAB saves
each expression with -v7.3, -v7 and -v6 and read_mat must return the value; write_mat writes the
value in each version and MATLAB's load must give the expression's value with the same class,
size, complexity, sparsity and bits. Without MATLAB the cases run Python round trips and check the
HDF5 layout MATLAB uses.
"""

from __future__ import annotations

import dataclasses
import math
import os
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

VERSIONS = ("7.3", "7", "6")
# MATLAB's NaN on x86-64 has the sign bit set (0xFFF8...), unlike Python's float("nan")
NAN = np.array([0xFFF8000000000000], dtype=np.uint64).view(np.float64)[0]


@dataclasses.dataclass(frozen=True)
class MatCase:
    """name: variable name in the file; matlab: expression of the value; python: what read_mat
    returns; inputs: more values write_mat must store as the same MATLAB value; write: False for
    values only MATLAB can write (objects, strings, function handles); versions: the versions
    write_mat writes it in; ml: the versions MATLAB saves it in without loss"""

    name: str
    matlab: str
    python: Callable[[], Any]
    inputs: tuple[Callable[[], Any], ...] = ()
    write: bool = True
    versions: tuple[str, ...] = VERSIONS
    ml: tuple[str, ...] = VERSIONS


# MATLAB's -v6 writes char through the session's default encoding: in a C locale every
# character outside ASCII becomes 26 (SUB), and a surrogate pair is an error
_NO_V6 = ("7.3", "7")


def _f(a, dtype=np.float64, shape=None):
    """an F-ordered array, filled column-major like MATLAB's reshape"""
    a = np.asarray(a, dtype=dtype)
    return np.asfortranarray(a.reshape(shape, order="F") if shape is not None else a)


def _row(*v, dtype=np.float64):
    return _f(v, dtype, (1, len(v)))


def _col(*v, dtype=np.float64):
    return _f(v, dtype, (len(v), 1))


def _empty(*shape, dtype=np.float64):
    return np.zeros(shape, dtype=dtype, order="F")


def _chars(s, shape):
    """a char array that read_mat returns as a numpy array (not a row): column-major fill"""
    return _f(list(s) if s else [], "<U1", shape)


def _cell(items, shape):
    a = np.empty(shape, dtype=object, order="F")
    flat = a.reshape(-1, order="F")
    for i, v in enumerate(items):
        flat[i] = v
    return a


def _structs(shape, fields, rows):
    from petakit5d.io.matfile import StructArray

    return StructArray.from_list(rows, shape=shape, fields=fields)


def _ramp(n, dtype, shape):
    """values 0..n-1 times a prime, modulo the class range: covers every byte"""
    v = (np.arange(n, dtype=np.uint64) * np.uint64(40503)) % np.uint64(65521)
    return _f(v.astype(dtype), dtype, shape)


_INT_LIMITS = [
    (t, np.iinfo(t))
    for t in (np.int8, np.uint8, np.int16, np.uint16, np.int32, np.uint32, np.int64, np.uint64)
]


def _int_cases():
    out = []
    for t, info in _INT_LIMITS:
        c = np.dtype(t).name
        out.append(
            MatCase(
                f"{c}_limits",
                f"[intmin('{c}') 0 intmax('{c}')]",
                lambda t=t, info=info: _f([info.min, 0, info.max], t, (1, 3)),
            )
        )
        out.append(MatCase(f"{c}_scalar", f"{c}(7)", lambda t=t: t(7)))
        out.append(MatCase(f"{c}_empty", f"zeros(0, 0, '{c}')", lambda t=t: _empty(0, 0, dtype=t)))
    return out


CASES: list[MatCase] = [
    # double
    MatCase("d_scalar", "1.5", lambda: np.float64(1.5), (lambda: 1.5, lambda: np.array(1.5))),
    MatCase("d_int", "3", lambda: np.float64(3), (lambda: 3,)),
    MatCase("d_zero", "0", lambda: np.float64(0), (lambda: 0,)),
    MatCase("d_negzero", "-0", lambda: np.float64(-0.0), (lambda: -0.0,)),
    MatCase("d_row", "[1 2 3]", lambda: _row(1, 2, 3), (lambda: np.array([1.0, 2, 3]),)),
    MatCase("d_col", "[1; 2; 3]", lambda: _col(1, 2, 3)),
    MatCase(
        "d_mat",
        "reshape(1 : 6, 2, 3)",
        lambda: _f(range(1, 7), shape=(2, 3)),
        (lambda: np.ascontiguousarray(_f(range(1, 7), shape=(2, 3))),),
    ),
    MatCase("d_3d", "reshape(1 : 24, 2, 3, 4)", lambda: _f(range(1, 25), shape=(2, 3, 4))),
    MatCase("d_4d", "reshape(1 : 48, 2, 3, 2, 4)", lambda: _f(range(1, 49), shape=(2, 3, 2, 4))),
    MatCase(
        "d_trailing_one", "zeros(2, 3, 1, 1)", lambda: _empty(2, 3), (lambda: _empty(2, 3, 1, 1),)
    ),
    MatCase("d_middle_one", "zeros(2, 1, 3)", lambda: _empty(2, 1, 3)),
    MatCase(
        "d_special",
        "[NaN Inf -Inf -0 realmin realmax eps pi 2^-1074]",
        lambda: _row(
            NAN,
            math.inf,
            -math.inf,
            -0.0,
            2.2250738585072014e-308,
            1.7976931348623157e308,
            2.220446049250313e-16,
            math.pi,
            5e-324,
        ),
    ),
    MatCase(
        "d_nan_payload",
        "typecast(uint8([1 0 0 0 0 0 248 127]), 'double')",
        lambda: np.array([1, 0, 0, 0, 0, 0, 248, 127], dtype=np.uint8).view(np.float64)[0],
    ),
    MatCase("d_empty", "[]", lambda: _empty(0, 0), (lambda: None, lambda: np.zeros((0, 0)))),
    MatCase("d_empty_0x3", "zeros(0, 3)", lambda: _empty(0, 3)),
    MatCase("d_empty_3x0", "zeros(3, 0)", lambda: _empty(3, 0)),
    MatCase("d_empty_1x0", "zeros(1, 0)", lambda: _empty(1, 0), (lambda: np.zeros(0),)),
    MatCase("d_empty_2x0x3", "zeros(2, 0, 3)", lambda: _empty(2, 0, 3)),
    MatCase(
        "d_large",
        "reshape(mod((0 : 299999) * 40503, 65521), 300, 1000)",
        lambda: _ramp(300000, np.float64, (300, 1000)),
    ),
    # single
    MatCase("s_scalar", "single(1.25)", lambda: np.float32(1.25)),
    MatCase(
        "s_mat",
        "single(reshape([1 2.5 -3 NaN Inf 1e-40], 2, 3))",
        lambda: _f([1, 2.5, -3, NAN, math.inf, 1e-40], np.float32, (2, 3)),
    ),
    MatCase("s_empty", "single([])", lambda: _empty(0, 0, dtype=np.float32)),
    MatCase(
        "s_3d", "single(reshape(1 : 60, 3, 4, 5))", lambda: _f(range(1, 61), np.float32, (3, 4, 5))
    ),
    # integers
    *_int_cases(),
    MatCase(
        "uint16_volume",
        "reshape(uint16(mod((0 : 59999) * 40503, 65521)), 30, 40, 50)",
        lambda: _ramp(60000, np.uint16, (30, 40, 50)),
    ),
    # complex
    MatCase("z_scalar", "complex(1, 2)", lambda: np.complex128(1 + 2j), (lambda: 1 + 2j,)),
    MatCase(
        "z_zero_imag", "complex(1, 0)", lambda: np.complex128(1 + 0j), (lambda: complex(1, 0),)
    ),
    MatCase(
        "z_row", "[1+2i, 3-4i, -0.5i]", lambda: _f([1 + 2j, 3 - 4j, -0.5j], np.complex128, (1, 3))
    ),
    MatCase("z_single", "single([1+2i; 3])", lambda: _f([1 + 2j, 3], np.complex64, (2, 1))),
    # logical
    MatCase("l_true", "true", lambda: np.bool_(True), (lambda: True,)),
    MatCase("l_false", "false", lambda: np.bool_(False), (lambda: False,)),
    MatCase(
        "l_row",
        "[true false true]",
        lambda: _row(1, 0, 1, dtype=np.bool_),
        (lambda: np.array([True, False, True]),),
    ),
    MatCase("l_eye", "logical(eye(3))", lambda: _f(np.eye(3), np.bool_)),
    MatCase("l_empty", "logical([])", lambda: _empty(0, 0, dtype=np.bool_)),
    MatCase("l_empty_0x3", "false(0, 3)", lambda: _empty(0, 3, dtype=np.bool_)),
    # char
    MatCase("c_row", "'abc'", lambda: "abc", (lambda: np.str_("abc"),)),
    MatCase("c_one", "'x'", lambda: "x"),
    MatCase("c_empty", "''", lambda: ""),
    MatCase("c_empty_1x0", "char(zeros(1, 0))", lambda: _chars("", (1, 0))),
    MatCase("c_col", "['a'; 'b']", lambda: _chars("ab", (2, 1))),
    MatCase(
        "c_matrix",
        "['abc'; 'def']",
        lambda: _chars("adbecf", (2, 3)),
        (lambda: np.array([list("abc"), list("def")]),),
    ),
    MatCase("c_unicode", "char([233 8364 937 32 955])", lambda: "é€Ω λ", ml=_NO_V6),
    MatCase("c_surrogates", "char([55357 56832 97])", lambda: "\U0001f600a", ml=_NO_V6),
    MatCase("c_lone_surrogate", "char([97 55357])", lambda: "a\ud83d", ml=_NO_V6),
    MatCase("c_null", "char([97 0 98])", lambda: "a\x00b"),
    MatCase(
        "c_path",
        "'/mnt/projects/embryo_2024/Scan_Iter_0000_CamA_ch0_000x_000y_000z_0000t.tif'",
        lambda: "/mnt/projects/embryo_2024/Scan_Iter_0000_CamA_ch0_000x_000y_000z_0000t.tif",
    ),
    # cell
    MatCase("cell_empty", "{}", list, (lambda: (),)),
    MatCase(
        "cell_row",
        "{1, 'a', [1 2]}",
        lambda: [np.float64(1), "a", _row(1, 2)],
        (lambda: (1, "a", np.array([1.0, 2])),),
    ),
    MatCase("cell_col", "{1; 'b'}", lambda: _cell([np.float64(1), "b"], (2, 1))),
    MatCase("cell_cellstr", "{'CamA', 'CamB'}", lambda: ["CamA", "CamB"]),
    MatCase("cell_empties", "cell(2, 3)", lambda: _cell([_empty(0, 0)] * 6, (2, 3))),
    MatCase("cell_1x0", "cell(1, 0)", lambda: _cell([], (1, 0))),
    MatCase("cell_0x3", "cell(0, 3)", lambda: _cell([], (0, 3))),
    MatCase("cell_nested", "{{}, {1, {'x'}}, {''}}", lambda: [[], [np.float64(1), ["x"]], [""]]),
    MatCase(
        "cell_3d",
        "reshape(num2cell(1 : 8), 2, 2, 2)",
        lambda: _cell([np.float64(i) for i in range(1, 9)], (2, 2, 2)),
    ),
    MatCase(
        "cell_mixed",
        "{uint8(5), single([1 2]), true, struct('a', 1), {}}",
        lambda: [
            np.uint8(5),
            _row(1, 2, dtype=np.float32),
            np.bool_(True),
            {"a": np.float64(1)},
            [],
        ],
    ),
    MatCase(
        "cell_mips",
        "{reshape(uint16(1 : 12), 3, 4); reshape(uint16(13 : 18), 2, 3)}",
        lambda: _cell(
            [_f(range(1, 13), np.uint16, (3, 4)), _f(range(13, 19), np.uint16, (2, 3))], (2, 1)
        ),
    ),
    # struct
    MatCase(
        "st_scalar",
        "struct('a', 1, 'b', 'x')",
        lambda: {"a": np.float64(1), "b": "x"},
        (lambda: {"a": 1, "b": "x"},),
    ),
    MatCase(
        "st_order",
        "struct('zeta', 1, 'alpha', 2, 'Mid', 3)",
        lambda: {"zeta": np.float64(1), "alpha": np.float64(2), "Mid": np.float64(3)},
    ),
    MatCase("st_nofields", "struct()", dict),
    MatCase(
        "st_nested",
        "struct('p', struct('q', struct('r', [1 2])), 's', {{}})",
        lambda: {"p": {"q": {"r": _row(1, 2)}}, "s": []},
    ),
    MatCase(
        "st_cell_field",
        "struct('c', {{'CamA', 'CamB'}}, 'e', {{}})",
        lambda: {"c": ["CamA", "CamB"], "e": []},
    ),
    MatCase(
        "st_empty_fields",
        "struct('a', [], 'b', '', 'c', zeros(0, 3), 'd', {{}}, 'e', false(0, 0))",
        lambda: {
            "a": _empty(0, 0),
            "b": "",
            "c": _empty(0, 3),
            "d": [],
            "e": _empty(0, 0, dtype=np.bool_),
        },
        (
            lambda: {
                "a": None,
                "b": "",
                "c": np.zeros((0, 3)),
                "d": [],
                "e": np.zeros((0, 0), bool),
            },
        ),
    ),
    MatCase(
        "st_array_row",
        "struct('a', {1, 'two', [3 3]}, 'b', {[], true, {}})",
        lambda: _structs(
            (1, 3),
            ("a", "b"),
            [
                {"a": np.float64(1), "b": _empty(0, 0)},
                {"a": "two", "b": np.bool_(True)},
                {"a": _row(3, 3), "b": []},
            ],
        ),
    ),
    MatCase(
        "st_array_col",
        "struct('a', {1; 2})",
        lambda: _structs((2, 1), ("a",), [{"a": np.float64(1)}, {"a": np.float64(2)}]),
    ),
    MatCase(
        "st_array_2x2",
        "reshape(struct('v', {1, 2, 3, 4}), 2, 2)",
        lambda: _structs((2, 2), ("v",), [{"v": np.float64(i)} for i in range(1, 5)]),
    ),
    MatCase(
        "st_empty_with_fields", "struct('a', {}, 'b', {})", lambda: _structs((0, 0), ("a", "b"), [])
    ),
    MatCase("st_empty_1x0", "repmat(struct('a', 1), 1, 0)", lambda: _structs((1, 0), ("a",), [])),
    MatCase("st_empty_nofields", "struct([])", lambda: _structs((0, 0), (), [])),
    MatCase("st_array_nofields", "repmat(struct(), 1, 3)", lambda: _structs((1, 3), (), [{}] * 3)),
    MatCase(
        "st_long_field",
        "struct('a23456789012345678901234567890123456789012345678901234567890123', 5)",
        lambda: {"a23456789012345678901234567890123456789012345678901234567890123": np.float64(5)},
    ),
    MatCase(
        "st_in_cell",
        "{struct('a', {1, 2})}",
        lambda: [_structs((1, 2), ("a",), [{"a": np.float64(1)}, {"a": np.float64(2)}])],
    ),
    # sparse
    MatCase(
        "sp_double", "sparse([1 0 0; 0 0 2.5])", lambda: _sparse([[1, 0, 0], [0, 0, 2.5]], float)
    ),
    MatCase("sp_logical", "sparse(logical([0 1; 1 0]))", lambda: _sparse([[0, 1], [1, 0]], bool)),
    MatCase("sp_empty", "sparse(3, 4)", lambda: _sparse(np.zeros((3, 4)), float)),
    MatCase("sp_complex", "sparse([1i 0; 0 2])", lambda: _sparse([[1j, 0], [0, 2]], complex)),
    # shaped like PetaKit5D's files: a parameters.mat 'pr' struct, the MIPs of a MIP flag file
    MatCase(
        "prod_parameters",
        "struct('dataPaths', {{'/mnt/projects/embryo_2024/'}}, 'outDirStr', '', "
        "'channelPatterns', {{'CamA', 'CamB'}}, 'dz', 0.5, 'skewAngle', 32.45, "
        "'flipZstack', false, 'inputBbox', [], 'blockSize', [256 256 256], "
        "'FFImagePaths', {{''}}, 'parseCluster', true, 'cpusPerTask', 1, 'uuid', 'a1b2', "
        "'resampleFactor', [], 'maskFullpaths', {{}})",
        lambda: {
            "dataPaths": ["/mnt/projects/embryo_2024/"],
            "outDirStr": "",
            "channelPatterns": ["CamA", "CamB"],
            "dz": np.float64(0.5),
            "skewAngle": np.float64(32.45),
            "flipZstack": np.bool_(False),
            "inputBbox": _empty(0, 0),
            "blockSize": _row(256, 256, 256),
            "FFImagePaths": [""],
            "parseCluster": np.bool_(True),
            "cpusPerTask": np.float64(1),
            "uuid": "a1b2",
            "resampleFactor": _empty(0, 0),
            "maskFullpaths": [],
        },
        (
            lambda: {
                "dataPaths": ["/mnt/projects/embryo_2024/"],
                "outDirStr": "",
                "channelPatterns": ("CamA", "CamB"),
                "dz": 0.5,
                "skewAngle": 32.45,
                "flipZstack": False,
                "inputBbox": None,
                "blockSize": np.array([256.0, 256, 256]),
                "FFImagePaths": [""],
                "parseCluster": True,
                "cpusPerTask": 1,
                "uuid": "a1b2",
                "resampleFactor": None,
                "maskFullpaths": [],
            },
        ),
    ),
    MatCase(
        "prod_mip_flag",
        "{reshape(uint16(1 : 2048), 1, 64, 32); reshape(uint16(1 : 1024), 1, 32, 32)}",
        lambda: _cell(
            [
                _f(range(1, 2049), np.uint16, (1, 64, 32)),
                _f(range(1, 1025), np.uint16, (1, 32, 32)),
            ],
            (2, 1),
        ),
    ),
    # MATLAB-only values: read, never written
    MatCase("obj_string", '"abc"', lambda: "abc", write=False),
    MatCase(
        "obj_string_row",
        '["CamA", "", "Camé", string(missing)]',
        lambda: _f(["CamA", "", "Camé", None], object, (1, 4)),
        write=False,
    ),
    MatCase("obj_string_col", '["a"; "bc"]', lambda: _f(["a", "bc"], object, (2, 1)), write=False),
    MatCase(
        "obj_string_in_struct",
        "struct('f', \"x\", 'g', {{\"y\", 1}})",
        lambda: {"f": "x", "g": ["y", np.float64(1)]},
        write=False,
    ),
    MatCase("obj_inputparser", "inputParser", lambda: _object("inputParser"), write=False),
    MatCase("obj_function_handle", "@(x) x + 1", lambda: _object("function_handle"), write=False),
]


def _sparse(a, dtype):
    import scipy.sparse

    return scipy.sparse.csc_matrix(np.asarray(a, dtype=dtype))


def _object(cls):
    from petakit5d.io.matfile import MatlabObject

    return MatlabObject(cls)


# ------------------------------------------------------------------------------------------------
# comparison of Python values: types, classes, shapes and bits
# ------------------------------------------------------------------------------------------------


def assert_same(got: Any, want: Any, where: str = "") -> None:
    """exact equality: same Python type, dtype, shape, memory order irrelevant, values bit for bit
    (NaN payloads, -0), dict key order"""
    import scipy.sparse

    from petakit5d.io.matfile import MatlabObject, StructArray

    w = where or "value"
    assert type(got) is type(want), f"{w}: type {type(got).__name__} != {type(want).__name__}"
    if isinstance(want, dict):
        assert list(got) == list(want), f"{w}: fields {list(got)} != {list(want)}"
        for k in want:
            assert_same(got[k], want[k], f"{w}.{k}")
    elif isinstance(want, list):
        assert len(got) == len(want), f"{w}: length {len(got)} != {len(want)}"
        for i, (g, v) in enumerate(zip(got, want, strict=True)):
            assert_same(g, v, f"{w}[{i}]")
    elif isinstance(want, StructArray):
        assert got.shape == want.shape, f"{w}: shape {got.shape} != {want.shape}"
        assert got.fields == want.fields, f"{w}: fields {got.fields} != {want.fields}"
        for i, (g, v) in enumerate(zip(got.ravel(), want.ravel(), strict=True)):
            assert_same(g, v, f"{w}({i + 1})")
    elif isinstance(want, MatlabObject):
        assert got.class_name == want.class_name, (
            f"{w}: class {got.class_name} != {want.class_name}"
        )
    elif scipy.sparse.issparse(want):
        assert got.format == "csc" and got.dtype == want.dtype and got.shape == want.shape, (
            f"{w}: {got.format} {got.dtype} {got.shape} != csc {want.dtype} {want.shape}"
        )
        assert (got != want).nnz == 0, f"{w}: sparse values differ"
    elif isinstance(want, np.ndarray):
        assert got.dtype == want.dtype, f"{w}: dtype {got.dtype} != {want.dtype}"
        assert got.shape == want.shape, f"{w}: shape {got.shape} != {want.shape}"
        if want.dtype == object:
            for i, (g, v) in enumerate(
                zip(got.ravel(order="F"), want.ravel(order="F"), strict=True)
            ):
                assert_same(g, v, f"{w}{{{i + 1}}}")
        elif want.dtype.kind in "fc":
            gb = np.ascontiguousarray(got).view(np.uint8)
            wb = np.ascontiguousarray(want).view(np.uint8)
            assert np.array_equal(gb, wb), f"{w}: values differ: {got!r} != {want!r}"
        else:
            assert np.array_equal(got, want), f"{w}: values differ: {got!r} != {want!r}"
    elif isinstance(want, np.generic):
        assert got.dtype == want.dtype, f"{w}: dtype {got.dtype} != {want.dtype}"
        assert np.asarray(got).tobytes() == np.asarray(want).tobytes(), f"{w}: {got!r} != {want!r}"
    else:
        assert got == want, f"{w}: {got!r} != {want!r}"


# ------------------------------------------------------------------------------------------------
# Python round trips (every system, no MATLAB)
# ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("case", [c for c in CASES if c.write], ids=lambda c: c.name)
def test_round_trip(case: MatCase, version: str, tmp_path: Path) -> None:
    from petakit5d.io.matfile import read_mat, write_mat

    want = case.python()
    for i, make in enumerate((case.python, *case.inputs)):
        fn = tmp_path / f"{case.name}_{i}.mat"
        write_mat(fn, {case.name: make()}, version=version)
        got = read_mat(fn)
        assert list(got) == [case.name]
        assert_same(got[case.name], want, f"{case.name} (input {i}, -v{version})")


def test_many_variables_and_selection(tmp_path: Path) -> None:
    from petakit5d.io.matfile import read_mat, write_mat

    values = {c.name: c.python() for c in CASES if c.write}
    for version in VERSIONS:
        fn = tmp_path / f"all_{version}.mat"
        write_mat(fn, values, version=version)
        got = read_mat(fn)
        assert sorted(got) == sorted(values)
        for k, v in values.items():
            assert_same(got[k], v, f"{k} (-v{version})")
        part = read_mat(fn, ["st_scalar", "d_row"])
        assert list(part) == ["st_scalar", "d_row"]
        with pytest.warns(UserWarning, match="not found"):
            assert read_mat(fn, ["no_such_variable"]) == {}


@pytest.mark.parametrize("version", VERSIONS)
def test_version_detection(version: str, tmp_path: Path) -> None:
    from petakit5d.io.matfile import mat_version, write_mat

    fn = tmp_path / "v.mat"
    write_mat(fn, {"a": 1}, version=version)
    assert mat_version(fn) == ("7.3" if version == "7.3" else "5")
    head = fn.read_bytes()[:128]
    if version == "7.3":
        assert head.startswith(b"MATLAB 7.3 MAT-file, Platform: ")
        assert b"HDF5 schema 1.00 ." in head[:116]
        assert head[124:128] == b"\x00\x02IM"
    else:
        assert head.startswith(b"MATLAB 5.0 MAT-file")
        assert head[126:128] == b"IM"


def test_v73_layout_as_matlab_writes_it(tmp_path: Path) -> None:
    """the attributes, types and shapes MATLAB R2024b writes (seen in its files)"""
    import h5py

    from petakit5d.io.matfile import write_mat

    fn = tmp_path / "layout.mat"
    write_mat(
        fn,
        {
            "d": np.arange(6.0).reshape(2, 3),
            "c": "ab",
            "e": "",
            "l": True,
            "z": None,
            "cl": [1, "x"],
            "s": {"b": 1, "a": 2},
            "sa": _structs((1, 2), ("f",), [{"f": 1.0}, {"f": 2.0}]),
            "u": np.uint16(3),
        },
    )
    with h5py.File(fn, "r") as f:
        assert f.userblock_size == 512

        def cls(o):
            a = h5py.h5a.open(o.id, b"MATLAB_class")
            t = a.get_type()
            assert t.get_strpad() == h5py.h5t.STR_NULLTERM and t.get_cset() == h5py.h5t.CSET_ASCII
            assert a.get_space().get_simple_extent_type() == h5py.h5s.SCALAR
            return o.attrs["MATLAB_class"].decode()

        assert cls(f["d"]) == "double" and f["d"].shape == (3, 2)
        assert np.array_equal(f["d"][()], np.arange(6.0).reshape(2, 3).T)
        assert cls(f["c"]) == "char" and f["c"].dtype == np.uint16 and f["c"].shape == (2, 1)
        assert (
            f["c"].attrs["MATLAB_int_decode"] == 2
            and f["c"].attrs["MATLAB_int_decode"].dtype == np.int32
        )
        for k, c in (("e", "char"), ("z", "double")):
            assert cls(f[k]) == c and f[k].dtype == np.uint64
            assert list(f[k][()]) == [0, 0] and f[k].attrs["MATLAB_empty"] == 1
            assert f[k].attrs["MATLAB_empty"].dtype == np.uint8
        assert cls(f["l"]) == "logical" and f["l"].dtype == np.uint8
        assert f["l"].attrs["MATLAB_int_decode"] == 1
        assert (
            cls(f["cl"]) == "cell" and f["cl"].dtype == h5py.ref_dtype and f["cl"].shape == (2, 1)
        )
        assert f[f["cl"][0, 0]].name.startswith("/#refs#/")
        assert cls(f["#refs#/a"]) == "canonical empty"
        assert isinstance(f["s"], h5py.Group) and cls(f["s"]) == "struct"
        assert [x.tobytes().decode() for x in f["s"].attrs["MATLAB_fields"]] == ["b", "a"]
        assert cls(f["sa"]) == "struct" and f["sa/f"].dtype == h5py.ref_dtype
        assert f["sa/f"].shape == (2, 1)
        assert cls(f["u"]) == "uint16" and f["u"].shape == (1, 1)


def _v5_elements(buf: bytes, pos: int = 0, end: int | None = None) -> list[tuple[int, bytes, bool]]:
    """(type, data, small format) of the data elements in buf (MAT 5, little-endian)"""
    import struct
    import zlib

    out = []
    end = len(buf) if end is None else end
    while pos < end:
        t, n = struct.unpack_from("<II", buf, pos)
        if t >> 16:
            out.append((t & 0xFFFF, buf[pos + 4 : pos + 4 + (t >> 16)], True))
            pos += 8
        elif t == 15:
            out.extend(_v5_elements(zlib.decompress(buf[pos + 8 : pos + 8 + n])))
            pos += 8 + n
        else:
            out.append((t, buf[pos + 8 : pos + 8 + n], False))
            pos += 8 + n + (-n % 8)
    return out


@pytest.mark.parametrize("version", ("7", "6"))
def test_v5_layout_as_matlab_writes_it(version: str, tmp_path: Path) -> None:
    """the element types and flags MATLAB R2024b writes (seen in its -v6 and -v7 files)"""
    import struct

    from petakit5d.io.matfile import write_mat

    def matrix(value, name="abcd"):
        fn = tmp_path / "x.mat"
        write_mat(fn, {name: value}, version=version)
        buf = fn.read_bytes()
        assert buf[116:128] == bytes(8) + b"\x00\x01IM"
        top = _v5_elements(buf, 128)
        assert top[0][0] == 14
        return _v5_elements(top[0][1])

    def flags(els):
        w0, nzmax = struct.unpack("<II", els[0][1])
        return w0 & 0xFF, (w0 >> 8) & 0xFF, nzmax

    els = matrix(3.0)
    assert flags(els) == (6, 0, 0)
    assert els[2] == (1, b"abcd", True)  # names of up to 4 bytes in the small format
    assert flags(matrix(True)) == (9, 0x02, 0)
    assert flags(matrix(np.float32(1) + 1j))[:2] == (7, 0x08)
    import scipy.sparse

    assert flags(matrix(scipy.sparse.csc_matrix((3, 4)))) == (5, 0x10, 1)
    assert flags(matrix(scipy.sparse.csc_matrix(np.eye(2, dtype=bool)))) == (5, 0x12, 2)
    # char: UTF-8 for ASCII in v7, uint16 in v6, UTF-16 for anything else
    assert matrix("abcdef")[3][0] == (16 if version == "7" else 4)
    assert matrix("abé")[3][0] == 17
    assert matrix("a\ud83d")[3][0] == 17
    # struct field names: one length for all, the longest name plus its terminating zero
    els = matrix({"ab": 1, "c": 2})
    assert els[3] == (5, struct.pack("<i", 3), True)
    assert els[4][1] == b"ab\x00c\x00\x00"
    assert matrix({})[3][1] == struct.pack("<i", 1)


@pytest.mark.parametrize("version", VERSIONS)
def test_write_errors(version: str, tmp_path: Path) -> None:
    from petakit5d.io.matfile import MatlabObject, write_mat

    fn = tmp_path / "e.mat"
    bad = [
        ({"1abc": 1}, ValueError),
        ({"a b": 1}, ValueError),
        ({"a" * 64: 1}, ValueError),
        ({"s": {"bad key": 1}}, ValueError),
        ({"s": {"_x": 1}}, ValueError),
        ({"x": 2**53 + 1}, ValueError),
        ({"x": {1, 2}}, TypeError),
        ({"x": object()}, TypeError),
        ({"x": np.array(["ab", "c"])}, ValueError),
        ({"x": MatlabObject("inputParser")}, TypeError),
        ({"x": np.datetime64("2026-01-01")}, TypeError),
    ]
    for variables, err in bad:
        with pytest.raises(err):
            write_mat(fn, variables, version=version)
    with pytest.raises(ValueError):
        write_mat(fn, {"a": 1}, version="5")


def test_large_ints_and_bytes(tmp_path: Path) -> None:
    from petakit5d.io.matfile import read_mat, write_mat

    for version in VERSIONS:
        fn = tmp_path / f"i{version}.mat"
        write_mat(fn, {"a": 2**53, "b": -(2**53), "c": b"\x00\xff"}, version=version)
        got = read_mat(fn)
        assert_same(got["a"], np.float64(2**53))
        assert_same(got["b"], np.float64(-(2**53)))
        assert_same(got["c"], _row(0, 255, dtype=np.uint8))


def test_compression_options(tmp_path: Path) -> None:
    import h5py

    from petakit5d.io.matfile import read_mat, write_mat

    big = _ramp(300000, np.float64, (300, 1000))
    write_mat(tmp_path / "c.mat", {"x": big, "s": np.arange(1000.0)})
    write_mat(tmp_path / "n.mat", {"x": big}, compression=False)
    with h5py.File(tmp_path / "c.mat", "r") as f:
        assert f["x"].compression == "gzip" and f["x"].compression_opts == 3
        assert f["x"].chunks == (27, 300)
        assert f["s"].id.get_create_plist().get_layout() == h5py.h5d.COMPACT
    with h5py.File(tmp_path / "n.mat", "r") as f:  # -nocompression: chunked, deflate level 0
        assert f["x"].compression == "gzip" and f["x"].compression_opts == 0
    for n in ("c", "n"):
        assert_same(read_mat(tmp_path / f"{n}.mat")["x"], big)
    write_mat(tmp_path / "v6.mat", {"x": big}, version="6")
    write_mat(tmp_path / "v7.mat", {"x": big}, version="7")
    assert (tmp_path / "v7.mat").stat().st_size < (tmp_path / "v6.mat").stat().st_size


# chunk shapes MATLAB R2024b chose (HDF5 order: MATLAB's size reversed), from its own files
MATLAB_CHUNKS = [
    ((9000, 1), 8, (4500, 1)),
    ((100000, 1), 8, (7692, 1)),
    ((1, 100000), 8, (1, 7692)),
    ((10000, 10), 8, (769, 10)),
    ((70000, 3), 8, (2692, 3)),
    ((3, 70000), 1, (3, 17500)),
    ((3, 300, 300), 8, (3, 300, 9)),
    ((3, 300, 300), 4, (3, 300, 17)),
    ((3, 300, 300), 1, (3, 300, 60)),
    ((50, 40, 30, 20), 8, (1, 13, 30, 20)),
    ((50, 40, 30, 20), 4, (1, 20, 30, 20)),
    ((50, 40, 30, 20), 1, (2, 40, 30, 20)),
    ((1000, 300), 8, (27, 300)),
    ((50, 40, 30), 2, (25, 40, 30)),
    ((100, 100), 16, (100, 33)),
    ((300, 300), 1, (300, 150)),
    ((8192, 1), 8, (8192, 1)),
    ((5000, 1), 8, (5000, 1)),
    ((40000, 1), 2, (20000, 1)),
    ((100000, 1), 4, (14285, 1)),
]


@pytest.mark.parametrize(("shape", "itemsize", "chunks"), MATLAB_CHUNKS)
def test_chunks_as_matlab(shape, itemsize, chunks) -> None:
    from petakit5d.io.matfile import _chunks

    assert _chunks(shape, itemsize) == chunks


def test_overwrite_is_whole_file(tmp_path: Path) -> None:
    from petakit5d.io.matfile import read_mat, write_mat

    fn = tmp_path / "o.mat"
    for version in VERSIONS:
        write_mat(fn, {"a": 1, "b": 2}, version=version)
        write_mat(fn, {"c": 3}, version="7.3" if version != "7.3" else "6")
        assert list(read_mat(fn)) == ["c"]


def test_structarray_api() -> None:
    from petakit5d.io.matfile import StructArray

    s = StructArray.from_list([{"a": 1, "b": 2}, {"a": 3, "b": 4}], shape=(1, 2))
    assert s.shape == (1, 2) and s.fields == ("a", "b") and len(s) == 2 and s.size == 2
    assert [e["a"] for e in s] == [1, 3]
    assert s[1]["b"] == 4 and s[0, 1]["b"] == 4
    assert s.field("a") == [1, 3]
    with pytest.raises(ValueError):
        StructArray.from_list([{"a": 1}, {"b": 2}])
    e = StructArray.empty(("x", "y"))
    assert e.shape == (0, 0) and e.fields == ("x", "y") and len(e) == 0
    assert repr(s).startswith("StructArray(1x2")


def test_read_not_a_mat_file(tmp_path: Path) -> None:
    from petakit5d.io.matfile import read_mat

    fn = tmp_path / "x.mat"
    fn.write_bytes(b"hello")
    with pytest.raises(ValueError, match="not a MAT-file"):
        read_mat(fn)
    with pytest.raises(FileNotFoundError):
        read_mat(tmp_path / "missing.mat")


# ------------------------------------------------------------------------------------------------
# save -ascii
# ------------------------------------------------------------------------------------------------

ASCII_CASES = [
    # (name, MATLAB expression, Python value, options); MATLAB: save(fn, 'x', '-ascii', opts{:})
    ("a_int", "12", 12, ()),
    ("a_frac", "0.105", 0.105, ()),
    ("a_neg_small", "-3.25e-10", -3.25e-10, ()),
    ("a_big", "1e100", 1e100, ()),
    ("a_zero", "0", 0, ()),
    ("a_negzero", "-0", -0.0, ()),
    ("a_special", "[NaN Inf -Inf]", np.array([math.nan, math.inf, -math.inf]), ()),
    ("a_matrix", "[1 2; 3 4.5]", np.array([[1, 2], [3, 4.5]]), ()),
    ("a_col", "[1; -2]", np.array([[1.0], [-2.0]]), ()),
    ("a_round", "[0.123456789 9.99999999]", np.array([0.123456789, 9.99999999]), ()),
    ("a_pi_double", "pi", math.pi, ("-double",)),
    ("a_tabs", "[1 2 3]", np.array([1.0, 2, 3]), ("-tabs",)),
    ("a_double_tabs", "[pi; -1e-300]", np.array([[math.pi], [-1e-300]]), ("-double", "-tabs")),
    ("a_single", "single(pi)", np.float32(math.pi), ()),
    ("a_uint16", "uint16([5 65535])", np.array([5, 65535], dtype=np.uint16), ()),
    ("a_logical", "true", True, ()),
    ("a_char", "'ab'", "ab", ()),
    ("a_empty", "[]", np.zeros((0, 0)), ()),
    ("a_3d", "reshape(1 : 8, 2, 2, 2)", _f(range(1, 9), shape=(2, 2, 2)), ()),
    ("a_complex", "[1+2i 3]", np.array([1 + 2j, 3]), ()),
    ("a_double_special", "[NaN -Inf 0.1]", np.array([math.nan, -math.inf, 0.1]), ("-double",)),
    ("a_1x0", "zeros(1, 0)", np.zeros((1, 0)), ()),
]


@pytest.mark.parametrize(
    ("name", "expr", "value", "opts"), ASCII_CASES, ids=[a[0] for a in ASCII_CASES]
)
def test_save_ascii_format(name, expr, value, opts, tmp_path: Path) -> None:
    """format rules of MATLAB's save -ascii (checked against MATLAB by the oracle cases)"""
    from petakit5d.io.matfile import save_ascii

    fn = tmp_path / f"{name}.txt"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        save_ascii(fn, value, double="-double" in opts, tabs="-tabs" in opts)
    assert bool(w) == (name in _ASCII_UNSUPPORTED)
    text = fn.read_bytes().decode("ascii")
    want = _ASCII_EXPECTED.get(name)
    if want is not None:
        assert text == want


# MATLAB warns for these: it writes nothing for the first four and the real part of a complex
_ASCII_UNSUPPORTED = {"a_single", "a_uint16", "a_logical", "a_3d", "a_complex"}

_ASCII_EXPECTED = {
    "a_int": "   1.2000000e+01\n",
    "a_matrix": "   1.0000000e+00   2.0000000e+00\n   3.0000000e+00   4.5000000e+00\n",
    "a_empty": "",
    "a_1x0": "\n",
    "a_complex": "   1.0000000e+00   3.0000000e+00\n",
    "a_tabs": "   1.0000000e+00\t   2.0000000e+00\t   3.0000000e+00\t\n",
    "a_special": "             NaN             Inf            -Inf\n",
}


def test_save_ascii_several_variables(tmp_path: Path) -> None:
    from petakit5d.io.matfile import save_ascii

    fn = tmp_path / "two.txt"
    save_ascii(fn, 12, np.array([1.0, 2.0]))
    assert fn.read_text() == "   1.2000000e+01\n   1.0000000e+00   2.0000000e+00\n"
    for bad in ({"a": 1}, [1, 2], [1, "a"], (1,)):
        with pytest.raises(TypeError):
            save_ascii(fn, bad)


# ------------------------------------------------------------------------------------------------
# against MATLAB (the oracle; references made on the cluster)
# ------------------------------------------------------------------------------------------------


def oracle_spec() -> dict:
    """the cases as data for the MATLAB side: expressions, versions, ASCII cases"""
    return {
        "cases": [
            {
                "name": c.name,
                "matlab": c.matlab,
                "write": c.write,
                "versions": list(c.versions) if c.write else list(VERSIONS),
            }
            for c in CASES
        ],
        "ascii": [{"name": n, "matlab": e, "opts": list(o)} for n, e, _, o in ASCII_CASES],
    }


def write_python_files(out: Path) -> None:
    """the files write_mat and save_ascii make for MATLAB to load: <out>/py<version>/<case>_<i>.mat
    (i: 0 for the value read_mat returns, then each other input) and <out>/ascii/<case>.txt"""
    from petakit5d.io.matfile import save_ascii, write_mat

    for v in VERSIONS:
        (out / f"py{v}").mkdir(parents=True, exist_ok=True)
    for c in CASES:
        if not c.write:
            continue
        for v in c.versions:
            for i, make in enumerate((c.python, *c.inputs)):
                write_mat(out / f"py{v}" / f"{c.name}_{i}.mat", {c.name: make()}, version=v)
    (out / "ascii").mkdir(parents=True, exist_ok=True)
    for name, _, value, opts in ASCII_CASES:
        save_ascii(
            out / "ascii" / f"{name}.txt", value, double="-double" in opts, tabs="-tabs" in opts
        )


def _matlab_ref_dir() -> Path | None:
    """MATLAB's results for these cases (stand-in for the F04 store until it lands): the folder
    PETAKIT5D_MATFILE_REF with ml<version>/<case>.mat, ascii/<case>.txt and load_results.txt"""
    ref = os.environ.get("PETAKIT5D_MATFILE_REF")
    return Path(ref) if ref and Path(ref).is_dir() else None


@pytest.mark.oracle
@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_matlab_reads(case: MatCase, version: str) -> None:
    """read_mat reads what MATLAB R2024b saves"""
    from petakit5d.io.matfile import read_mat

    ref = _matlab_ref_dir()
    if ref is None:
        pytest.skip("no MATLAB references (cluster only)")
    if version not in case.ml:
        pytest.skip(f"MATLAB's -v{version} does not keep this value")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got = read_mat(ref / f"ml{version}" / f"{case.name}.mat")
    assert list(got) == [case.name]
    assert_same(got[case.name], case.python(), f"{case.name} (MATLAB -v{version})")


@pytest.mark.oracle
@pytest.mark.parametrize("case", [c for c in CASES if c.write], ids=lambda c: c.name)
def test_matlab_loads(case: MatCase) -> None:
    """MATLAB loads what write_mat writes: the MATLAB side's verdicts (load_results.txt)"""
    ref = _matlab_ref_dir()
    if ref is None:
        pytest.skip("no MATLAB references (cluster only)")
    verdicts = {}
    for line in (ref / "load_results.txt").read_text().splitlines():
        key, _, verdict = line.partition(" ")
        verdicts[key] = verdict
    n = 1 + len(case.inputs)
    for v in case.versions:
        for i in range(n):
            key = f"py{v}/{case.name}_{i}"
            assert verdicts.get(key) == "ok", f"{key}: {verdicts.get(key, 'not loaded')}"


@pytest.mark.oracle
@pytest.mark.parametrize(
    ("name", "expr", "value", "opts"), ASCII_CASES, ids=[a[0] for a in ASCII_CASES]
)
def test_matlab_ascii(name, expr, value, opts, tmp_path: Path) -> None:
    """save_ascii writes the bytes MATLAB's save -ascii writes"""
    from petakit5d.io.matfile import save_ascii

    ref = _matlab_ref_dir()
    if ref is None:
        pytest.skip("no MATLAB references (cluster only)")
    fn = tmp_path / f"{name}.txt"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        save_ascii(fn, value, double="-double" in opts, tabs="-tabs" in opts)
    assert fn.read_bytes() == (ref / "ascii" / f"{name}.txt").read_bytes()
