"""The io wrappers' arguments-block declarations against MATLAB: conversions copied and refused.

``readtiff``, ``writetiff``, ``readzarr``, ``writezarr`` and ``createzarr`` declare their 36
inputs in arguments blocks (11 declarations: ``DECLS``). MATLAB converts a value to the declared
class and size before the function runs, then validates it. The MATLAB-name layer does the same
(``params.declarations``), and refuses the conversions that only hide mistakes, with the port's
error (``petakit5d:params:conversionRefused``; decided by the user on 2026-11-22).

MATLAB R2024b decides: each oracle case calls a fake function with one declaration
(``matlab_decl/pkpy_decl_<key>.m``: the lab's text and default, the argument positional or
name-value as the lab's are) on a battery of values of every class and many sizes, and records
for each value 'accepted <class> <size> <value>' (the value as MATLAB converted it) or MATLAB's
error identifier. The Python side is the registry's declaration of the lab's parameters.

- ``decl.<key>``: the values whose conversion the port copies: the same outcome for each, in
  every form the registry takes (a NumPy array, an ``MArray``, a tuple, a Python scalar);
- ``decl.<key>.refused``: the values whose conversion the port refuses (numbers, a cell of char
  arrays and a function handle for char, characters for logical): MATLAB's outcome is recorded
  (each converted: none is ``MATLAB:validation:UnableToConvert``), and the test checks that the
  port refuses each, in every form.

The battery holds no string array of other than one element: the registry's values have no form
for them.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from petakit5d.params import (
    ConversionRefusedError,
    MatlabArgumentError,
    MatlabValidationError,
    Param,
    Registry,
)
from petakit5d.params import values as mv
from petakit5d.params.declarations import REFUSED_ID, Declaration, read, value
from petakit5d.params.registry import REGISTRY
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, mcell, mcol, mempty, mstring, mtyped
from petakit5d.testing.compare import MArray, MCell, MFunc, MString, MStruct, MValue, to_mvalue

HELPERS = Path(__file__).resolve().parent / "matlab_decl"

COMPRESSION = "(1, :) char {mustBeMember(options.Compression, {'none', 'lzw', 'zstd'})}"
MODE = "(1, :) char {mustBeMember(options.Mode, {'parallel', 'libtiff', 'imwrite'})}"
# key: (declaration, how the fake takes the argument (p1, p2: its first or second positional
# argument; n: a name-value argument), the argument's name, the lab's parameters)
DECLS = {
    "char_p1": (
        "char",
        "p1",
        "filepath",
        [("readtiff", "filepath"), ("readzarr", "filepath"), ("createzarr", "filepath")],
    ),
    "char_p2": ("char", "p2", "filepath", [("writetiff", "filepath"), ("writezarr", "filepath")]),
    "char": (
        "char",
        "n",
        "dtype",
        [
            ("createzarr", "dtype"),
            ("createzarr", "order"),
            ("createzarr", "compressor"),
            ("createzarr", "dimSeparator"),
            ("writezarr", "dimSeparator"),
        ],
    ),
    "flag": (
        "(1, 1) logical",
        "n",
        "groupWrite",
        [
            ("writetiff", "groupWrite"),
            ("writezarr", "overwrite"),
            ("writezarr", "expand2dDim"),
            ("writezarr", "groupWrite"),
            ("writezarr", "uuid"),
            ("createzarr", "overwrite"),
            ("createzarr", "expand2dDim"),
            ("createzarr", "groupWrite"),
        ],
    ),
    "row": (
        "(1, :) {mustBeNumeric}",
        "n",
        "dataSize",
        [
            ("createzarr", "dataSize"),
            ("createzarr", "blockSize"),
            ("createzarr", "shardSize"),
            ("createzarr", "zarrSubSize"),
            ("writezarr", "blockSize"),
            ("writezarr", "shardSize"),
            ("writezarr", "zarrSubSize"),
            ("writezarr", "bbox"),
        ],
    ),
    "scalar": ("(1, 1) {mustBeNumeric}", "n", "clevel", [("createzarr", "clevel")]),
    "finite_row": (
        "(1, :) {mustBeNumeric, mustBeFinite}",
        "n",
        "range",
        [("readtiff", "range"), ("readzarr", "inputBbox")],
    ),
    "numeric_or_logical_row": (
        "(1, :) {mustBeNumericOrLogical}",
        "n",
        "sparseData",
        [("readzarr", "sparseData")],
    ),
    "numeric": ("{mustBeNumeric}", "p1", "img", [("writetiff", "img"), ("writezarr", "data")]),
    "compression": (COMPRESSION, "n", "Compression", [("writetiff", "Compression")]),
    "mode": (MODE, "n", "Mode", [("writetiff", "Mode")]),
    "logical_row": ("(1, :) logical", "n", "sparseData", [("writezarr", "sparseData")]),
    "logical": ("logical", "n", "create", [("writezarr", "create")]),
}


def _chars(*rows: str) -> np.ndarray:
    """a char matrix (rows of one length)"""
    return np.array([list(r) for r in rows], dtype="<U1")


def _nd(shape: tuple[int, ...]) -> np.ndarray:
    return np.arange(1, math.prod(shape) + 1, dtype=np.float64).reshape(shape, order="F")


def _cells(size: tuple[int, ...], items: list) -> MCell:
    return MCell(size, [x if isinstance(x, MValue) else to_mvalue(x) for x in items])


def _structs(size: tuple[int, ...], rows: list[dict]) -> MStruct:
    fields = list(rows[0]) if rows else []
    return MStruct(size, fields, [{k: to_mvalue(v) for k, v in r.items()} for r in rows])


# (name, group, value in MATLAB's form): the groups say which conversion a value meets
BATTERY = [
    ("true", "logical", True),
    ("false", "logical", False),
    ("logical_row2", "logical", np.array([True, False])),
    ("logical_col2", "logical", mcol([True, False], "logical")),
    ("logical_2x2", "logical", np.ones((2, 2), dtype=bool)),
    ("logical_1x1x2", "logical", np.ones((1, 1, 2), dtype=bool)),
    ("logical_empty", "logical", mempty("logical")),
    ("logical_1x0", "logical", mempty("logical", (1, 0))),
    ("zero", "number", 0.0),
    ("minus_zero", "number", -0.0),
    ("one", "number", 1.0),
    ("two", "number", 2.0),
    ("half", "number", 0.5),
    ("two_half", "number", 2.5),
    ("negative", "number", -3.0),
    ("nan", "number", math.nan),
    ("inf", "number", math.inf),
    ("minus_inf", "number", -math.inf),
    ("big", "number", 1e10),
    ("c65", "number", 65.0),
    ("c65_6", "number", 65.6),
    ("c_big", "number", 70000.0),
    ("tiny", "number", 1e-300),
    ("empty", "number", None),
    ("e1x0", "number", mempty("double", (1, 0))),
    ("e0x1", "number", mempty("double", (0, 1))),
    ("e0x3", "number", mempty("double", (0, 3))),
    ("e3x0", "number", mempty("double", (3, 0))),
    ("e1x1x0", "number", mempty("double", (1, 1, 0))),
    ("e0x0x3", "number", mempty("double", (0, 0, 3))),
    ("e1x0x2", "number", mempty("double", (1, 0, 2))),
    ("row3", "number", [1.0, 2.0, 3.0]),
    ("col3", "number", mcol([1.0, 2.0, 3.0])),
    ("nd_1x1x3", "number", _nd((1, 1, 3))),
    ("nd_1x1x1x3", "number", _nd((1, 1, 1, 3))),
    ("nd_1x2x3", "number", _nd((1, 2, 3))),
    ("m2x2", "number", _nd((2, 2))),
    ("m3x2", "number", _nd((3, 2))),
    ("nd_2x2x2", "number", _nd((2, 2, 2))),
    ("row_nan", "number", [1.0, math.nan]),
    ("col_nan", "number", mcol([math.nan, 2.0])),
    ("row_inf", "number", [1.0, 2.0, math.inf]),
    ("row_minus_inf", "number", [-math.inf, 3.0]),
    ("col_zero", "number", mcol([0.0, 1.0, 2.0])),
    ("uint8", "number", mtyped(3, "uint8")),
    ("uint8_zero", "number", mtyped(0, "uint8")),
    ("int8_col", "number", mcol([1, -2], "int8")),
    ("uint16_row3", "number", mtyped([1, 2, 3], "uint16")),
    ("int64_big", "number", mtyped(np.int64(9007199254740993), "int64")),
    ("single", "number", mtyped(1.5, "single")),
    ("single_nan", "number", mtyped(math.nan, "single")),
    ("single_inf", "number", mtyped(math.inf, "single")),
    ("single_col", "number", mcol([1.0, 2.0], "single")),
    ("codes_lzw", "number", [108.0, 122.0, 119.0]),
    ("codes_lzw_col", "number", mcol([108.0, 122.0, 119.0])),
    ("codes_lzw_u8", "number", mtyped([108, 122, 119], "uint8")),
    ("complex", "complex", 1 + 2j),
    ("complex_row", "complex", np.array([1 + 2j, 3 + 0j])),
    ("complex_inf", "complex", complex(math.inf, 1.0)),
    ("complex_zero", "complex", 0j),
    ("char", "char", "abc"),
    ("char_empty", "char", ""),
    ("char_a", "char", "a"),
    ("char_1x0", "char", mempty("char", (1, 0))),
    ("char_0x3", "char", mempty("char", (0, 3))),
    ("char_col", "char", _chars("a", "b")),
    ("char_col_lzw", "char", _chars("l", "z", "w")),
    ("char_matrix", "char", _chars("ab", "cd")),
    ("char_nd", "char", np.full((1, 1, 2), "a", dtype="<U1")),
    ("char_zero", "char", "0"),
    ("lzw", "char", "lzw"),
    ("upper_lzw", "char", "LZW"),
    ("zstd", "char", "zstd"),
    ("none", "char", "none"),
    ("lzw_blank", "char", "lzw "),
    ("blank_lzw", "char", " lzw"),
    ("lz", "char", "lz"),
    ("parallel", "char", "parallel"),
    ("upper_parallel", "char", "Parallel"),
    ("libtiff", "char", "libtiff"),
    ("imwrite", "char", "imwrite"),
    ("string", "string", mstring("abc")),
    ("string_empty", "string", mstring("")),
    ("string_lzw", "string", mstring("lzw")),
    ("string_upper_lzw", "string", mstring("LZW")),
    ("string_parallel", "string", mstring("parallel")),
    ("cell_empty", "cell_char", MCell((0, 0), [])),
    ("cell_1x0", "cell_char", mcell()),
    ("cell_a", "cell_char", mcell("a")),
    ("cell_lzw", "cell_char", mcell("lzw")),
    ("cell_upper_lzw", "cell_char", mcell("LZW")),
    ("cell_ab", "cell_char", mcell("a", "b")),
    ("cell_col", "cell_char", _cells((2, 1), ["lzw", "zstd"])),
    ("cell_empty_char", "cell_char", mcell("")),
    ("cell_char_col", "cell_char", mcell(_chars("a", "b"))),
    ("cell_num", "cell", mcell(1.0)),
    ("cell_nums", "cell", mcell(1.0, 2.0)),
    ("cell_true", "cell", mcell(True)),
    ("cell_string", "cell", mcell(mstring("lzw"))),
    ("cell_mixed", "cell", mcell("a", 1.0)),
    ("cell_cell", "cell", mcell(mcell("a"))),
    ("struct", "struct", {"a": 1.0}),
    ("struct_row", "struct", _structs((1, 2), [{"a": 1.0}, {"a": 2.0}])),
    ("struct_empty", "struct", MStruct((0, 0), [], [])),
    ("handle", "handle", MFunc("@sin")),
    ("handle_anon", "handle", MFunc("@(t)t+1")),
]
GROUP = {n: g for n, g, _ in BATTERY}
# the groups whose conversion the port refuses, by the declared class
REFUSED_GROUPS = {"char": ("number", "cell_char", "handle"), "logical": ("char",)}
# of the numbers, those the refused cases record (the refusal is the same for every number)
REFUSED_NUMBERS = (
    "zero one two_half negative nan inf minus_inf c65 c65_6 c_big empty e1x0 e0x1 e0x3 "
    "e1x1x0 row3 col3 nd_1x1x3 nd_1x2x3 m2x2 row_nan uint8 int64_big single single_col "
    "codes_lzw codes_lzw_col codes_lzw_u8"
).split()


def _anon(t):
    return t + 1


def _py(v: MValue) -> object:
    """A MATLAB value of the battery in the registry's form: char rows as str ('' 0x0), other
    char arrays and numeric and logical arrays as the registry's MArray (MATLAB's class and
    size), complex numbers as Python complex or NumPy arrays, cells as Cell, string scalars as
    MString, a 1x1 struct as a dict, other structs as StructArray, a handle as a callable"""
    from petakit5d.io.matfile import StructArray

    if isinstance(v, MArray):
        size = tuple(v.size)
        flat = v.data.reshape(-1, order="F")
        if v.mclass == "char":
            text = "".join(chr(int(c)) for c in flat)
            if size == (0, 0) or (len(size) == 2 and size[0] == 1 and size[1] > 0):
                return text
            return mv.MArray("char", size, list(text))
        if v.complex:
            return complex(flat[0]) if size == (1, 1) else v.data
        items = [bool(x) for x in flat] if v.mclass == "logical" else [x.item() for x in flat]
        return mv.MArray(v.mclass, size, items)
    if isinstance(v, MCell):
        return mv.Cell([_py(x) for x in v.items], tuple(v.size))
    if isinstance(v, MString):
        (text,) = v.items
        return mv.MString(text)
    if isinstance(v, MStruct):
        els = [{k: _py(x) for k, x in e.items()} for e in v.items]
        if tuple(v.size) == (1, 1):
            return els[0]
        if not els:
            return StructArray.empty(tuple(v.fields), tuple(v.size))
        return StructArray(tuple(v.size), tuple(v.fields), els)
    if isinstance(v, MFunc):
        return math.sin if v.text == "@sin" else _anon
    raise TypeError(f"no registry form for {v!r}")


_NP = {"double": np.float64, "single": np.float32, "logical": np.bool_, "char": "<U1"}


def _forms(v: object) -> list:
    """the battery value in each form the registry takes: its registry form, and the NumPy
    array, Python scalar, tuple (a column as rows of one), None or tuple of texts it can be"""
    port = _py(to_mvalue(v))
    out = [port]
    if isinstance(port, mv.MArray):
        dt = _NP.get(port.mclass) or np.dtype(port.mclass)
        out.append(np.array(port.data, dtype=dt).reshape(port.shape, order="F"))
        if port.mclass in ("double", "logical") and len(port.shape) == 2:
            vals = [bool(e) if port.mclass == "logical" else float(e) for e in port.data]
            if port.shape == (1, 1):
                out.append(vals[0])
            elif port.shape[0] == 1 and port.shape[1]:
                out.append(tuple(vals))
            elif port.shape[1] == 1 and port.shape[0] > 1:
                out.append(tuple((e,) for e in vals))
            elif port.shape == (0, 0) and port.mclass == "double":
                out.append(None)
    elif type(port) is str and port:  # (an MString is a str too: a string scalar)
        out.append(np.array(list(port), dtype="<U1"))
    elif (
        isinstance(port, mv.Cell)
        and port.shape[0] == 1
        and len(port)
        and all(type(e) is str for e in port)
    ):
        out.append(tuple(port))  # (an empty tuple is MATLAB's [], not a cell)
    return out


FORMS = {n: _forms(v) for n, _, v in BATTERY}


def _mat2str(a: np.ndarray) -> str:
    from petakit5d.core.text import mat2str

    return mat2str(a)


def describe(v) -> str:
    """an accepted value as pkpy_decl_describe.m records it"""
    s = f"accepted {v.mclass} {_mat2str(np.array([v.size], dtype=np.float64))}"
    if v.kind != "array" or len(v.size) != 2:
        return s
    a = v.array
    if a.size <= 12 and not v.complex and (v.mclass in mv.NUMERIC_CLASSES or v.mclass == "logical"):
        return f"{s} {_mat2str(a)}"
    if v.mclass == "char" and a.size <= 24:
        codes = np.vectorize(ord, otypes=[np.float64])(a) if a.size else np.zeros(a.shape)
        return f"{s} chars {_mat2str(codes)}"
    return s


def outcome(decl: Declaration, x: object, position: int | None) -> str:
    """the port's outcome of a value: as MATLAB's side records it; 'refused' for the port's
    refusal of a conversion"""
    try:
        v = decl.convert(x, "x", position)
    except ConversionRefusedError as err:
        assert err.identifier == REFUSED_ID
        return "refused"
    except MatlabValidationError as err:
        return err.identifier
    return describe(v)


def _position(key: str) -> int | None:
    form = DECLS[key][1]
    return {"p1": 1, "p2": 2}.get(form)


def _declaration(key: str) -> Declaration:
    text, _, _, params = DECLS[key]
    function, param = params[0]
    d = REGISTRY.find(function).param(param).declaration
    assert d is not None and d.text == text
    return d


def outcomes(key: str, names: list[str]) -> dict[str, str]:
    """the port's outcome of each value, the same in every form of it"""
    decl, position = _declaration(key), _position(key)
    out = {}
    for n in names:
        got = [outcome(decl, f, position) for f in FORMS[n]]
        if len(set(got)) != 1:
            raise AssertionError(f"{key} {n}: the forms of the value give {got}")
        out[n] = got[0]
    return out


def _refused(key: str, name: str) -> bool:
    cls = read(DECLS[key][0]).cls
    return GROUP[name] in REFUSED_GROUPS.get(cls, ())


def copied_names(key: str) -> list[str]:
    return [n for n, _, _ in BATTERY if not _refused(key, n)]


def refused_names(key: str) -> list[str]:
    return [
        n for n, g, _ in BATTERY if _refused(key, n) and (g != "number" or n in REFUSED_NUMBERS)
    ]


_CODE = (
    "v = in{1}; n = in{2}; f = in{3}; k = in{4}; a = in{5}; r = struct(); "
    "for i = 1 : numel(v); try; "
    "if k == 1, o = feval(f, v{i}); elseif k == 2, o = feval(f, 0, v{i}); "
    "else, o = feval(f, a, v{i}); end; "
    "r.(n{i}) = pkpy_decl_describe(o); "
    "catch ME; r.(n{i}) = ME.identifier; end; end; out = {r};"
)
_VALUES = {n: v for n, _, v in BATTERY}


def _case(key: str, names: list[str], refused: bool) -> Case:
    _, form, arg, _ = DECLS[key]
    return Case(
        f"decl.{key}" + (".refused" if refused else ""),
        code=_CODE,
        args=[
            mcell(*(_VALUES[n] for n in names)),
            mcell(*names),
            f"pkpy_decl_{key}",
            {"p1": 1.0, "p2": 2.0, "n": 3.0}[form],
            arg,
        ],
        python=None if refused else (lambda *_, key=key, names=names: outcomes(key, names)),
        helpers=("tests/params/matlab_decl",),
    )


CASES = [_case(k, copied_names(k), False) for k in DECLS]
REFUSED_CASES = [_case(k, refused_names(k), True) for k in DECLS if refused_names(k)]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def _texts(ref_value: MValue) -> dict[str, str] | None:
    """MATLAB's outcome per value, from a reference's struct (None where a manifest keeps its
    hash only)"""
    if not isinstance(ref_value, MStruct) or ref_value.items is None:
        return None
    (row,) = ref_value.items
    out = {}
    for k, x in row.items():
        if not isinstance(x, MArray) or x.data is None:
            return None
        out[k] = x.text()
    return out


@oracle.parametrize(REFUSED_CASES)
def test_matlab_refused(case, tmp_path):
    """MATLAB converts each of these values (its outcome recorded: accepted, or refused by a
    later check, never UnableToConvert); the port refuses each, in every form"""
    oracle.check(case, tmp_path)
    key = case.id.split(".")[1]
    names = refused_names(key)
    ref = oracle.reference(case)
    texts = _texts(ref.outputs[0])
    if texts is None:
        pytest.skip(f"case {case.id}: MATLAB's outcomes are not kept in full by the manifest")
    assert list(texts) == names
    unable = [n for n, t in texts.items() if t == "MATLAB:validation:UnableToConvert"]
    assert not unable, f"MATLAB does not convert {unable}: the port should copy its error"
    got = outcomes(key, names)
    assert all(o == "refused" for o in got.values()), got


# --- without MATLAB -------------------------------------------------------------------------------


def test_decls_lists_every_declaration_of_the_registry():
    """every registered arguments-block declaration is read whole and listed with its parameters
    (a new one needs a row, its fake function and its MATLAB references); every other registered
    text that is no inputParser validator has none"""
    found, unread = {}, []
    for s in REGISTRY.specs():
        for p in s.params:
            if p.declaration is not None:
                found[(s.matlab_name, p.matlab)] = p.declaration.text
            elif p.validator and not p.validator.lstrip().startswith("@"):
                unread.append((s.matlab_name, p.matlab, p.validator))
    assert not unread, unread
    want = {fp: text for text, _, _, params in DECLS.values() for fp in params}
    assert found == want
    assert len(found) == 36 and len({t for t, *_ in DECLS.values()}) == 11


@pytest.mark.parametrize("key", list(DECLS))
def test_the_fake_functions_hold_the_labs_texts(key):
    text, form, arg, _ = DECLS[key]
    src = (HELPERS / f"pkpy_decl_{key}.m").read_text(encoding="utf-8")
    line = f"options.{arg} {text} = " if form == "n" else f"    {arg} {text}\n"
    assert line in src
    head = {"p1": f"({arg})", "p2": f"(first, {arg})", "n": "(options)"}[form]
    assert src.startswith(f"function out = pkpy_decl_{key}{head}\n")


def test_the_declarations_are_read_once():
    for text, *_ in DECLS.values():
        assert read(text) is read(text) and read(text).text == text


@pytest.mark.parametrize(
    ("text", "size", "cls", "validators"),
    [
        ("char", None, "char", ()),
        ("(1, 1) logical", (1, 1), "logical", ()),
        ("logical", None, "logical", ()),
        ("(1, :) logical", (1, None), "logical", ()),
        ("{mustBeNumeric}", None, None, (("mustBeNumeric", ()),)),
        (
            "(1, :) {mustBeNumeric, mustBeFinite}",
            (1, None),
            None,
            (("mustBeNumeric", ()), ("mustBeFinite", ())),
        ),
        ("(1,:){mustBeNumericOrLogical}", (1, None), None, (("mustBeNumericOrLogical", ()),)),
        (COMPRESSION, (1, None), "char", (("mustBeMember", ("none", "lzw", "zstd")),)),
        (
            "char {mustBeMember(x, {'it''s', 'b'})}",
            None,
            "char",
            (("mustBeMember", ("it's", "b")),),
        ),
    ],
)
def test_reading(text, size, cls, validators):
    d = read(text)
    assert (d.size, d.cls, d.validators) == (size, cls, validators)


@pytest.mark.parametrize(
    "text",
    [
        None,
        "",
        "@ischar",
        "@(x) numel(x) == 2",
        "(1, 3) double",  # sizes other than (1, 1) and (1, :)
        "(:, 1) {mustBeNumeric}",
        "(1, 1, :) logical",
        "string",  # classes other than char and logical
        "(1, :) double",
        "{mustBePositive}",  # other validation functions
        "(1, :) {mustBeNumeric, mustBeInteger}",
        "char {mustBeMember(x, {1, 2})}",
        "char {mustBeMember(x, {})}",
        "(1, :",
        "(1, 1) logical {mustBeNumeric} extra",
    ],
)
def test_texts_not_read(text):
    """sizes, classes and validation functions the reading does not know give no declaration
    (the class alternatives check such a text, as before)"""
    assert read(text) is None


@pytest.mark.parametrize(("name", "group", "v"), BATTERY)
def test_the_groups_are_the_ports_classes(name, group, v):
    """each battery value's group is what the port sees in every form of it"""
    want = {
        "logical": ("array", "logical"),
        "complex": ("array", None),
        "char": ("array", "char"),
        "string": ("string", "string"),
        "cell_char": ("cell", "cell"),
        "cell": ("cell", "cell"),
        "struct": ("struct", "struct"),
        "handle": ("handle", "function_handle"),
    }
    for f in FORMS[name]:
        got = value(f)
        if group == "number":
            assert got.kind == "array" and got.mclass in mv.NUMERIC_CLASSES and not got.complex
        elif group == "complex":
            assert got.kind == "array" and got.complex
        else:
            assert (got.kind, got.mclass) == want[group]
        if group in ("cell_char", "cell"):
            chars = [value(e).kind == "array" and value(e).mclass == "char" for e in got.items]
            assert all(chars) == (group == "cell_char")


# each rule on both sides: (key, value taken (with its outcome), value refused, identifier)
SIDES = [
    ("row", "col3", "accepted double [1 3] [1 2 3]", "m3x2", "MATLAB:validation:IncompatibleSize"),
    (
        "row",
        "nd_1x1x3",
        "accepted double [1 3] [1 2 3]",
        "nd_1x2x3",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "row",
        "empty",
        "accepted double [1 0] zeros(1,0)",
        "e0x3",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "row",
        "e0x1",
        "accepted double [1 0] zeros(1,0)",
        "e3x0",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "row",
        "e1x1x0",
        "accepted double [1 0] zeros(1,0)",
        "e1x0x2",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "row",
        "int8_col",
        "accepted int8 [1 2] [1 -2]",
        "logical_row2",
        "MATLAB:validators:mustBeNumeric",
    ),
    ("row", "complex", "accepted double [1 1]", "char", "MATLAB:validators:mustBeNumeric"),
    (
        "row",
        "nan",
        "accepted double [1 1] NaN",
        "logical_2x2",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "finite_row",
        "row3",
        "accepted double [1 3] [1 2 3]",
        "row_inf",
        "MATLAB:validators:mustBeFinite",
    ),
    (
        "finite_row",
        "big",
        "accepted double [1 1] 10000000000",
        "col_nan",
        "MATLAB:validators:mustBeFinite",
    ),
    (
        "finite_row",
        "complex",
        "accepted double [1 1]",
        "complex_inf",
        "MATLAB:validators:mustBeFinite",
    ),
    (
        "finite_row",
        "single",
        "accepted single [1 1] 1.5",
        "single_nan",
        "MATLAB:validators:mustBeFinite",
    ),
    (
        "finite_row",
        "col_zero",
        "accepted double [1 3] [0 1 2]",
        "char_col",
        "MATLAB:validators:mustBeNumeric",
    ),
    ("scalar", "uint8", "accepted uint8 [1 1] 3", "row3", "MATLAB:validation:IncompatibleSize"),
    ("scalar", "nan", "accepted double [1 1] NaN", "empty", "MATLAB:validation:IncompatibleSize"),
    ("scalar", "one", "accepted double [1 1] 1", "true", "MATLAB:validators:mustBeNumeric"),
    (
        "numeric_or_logical_row",
        "logical_col2",
        "accepted logical [1 2] [true false]",
        "char",
        "MATLAB:validators:mustBeNumericOrLogical",
    ),
    (
        "numeric",
        "nd_2x2x2",
        "accepted double [2 2 2]",
        "logical_2x2",
        "MATLAB:validators:mustBeNumeric",
    ),
    ("flag", "two_half", "accepted logical [1 1] true", "nan", "MATLAB:validation:UnableToConvert"),
    (
        "flag",
        "minus_zero",
        "accepted logical [1 1] false",
        "logical_row2",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "flag",
        "int64_big",
        "accepted logical [1 1] true",
        "row_nan",
        "MATLAB:validation:UnableToConvert",
    ),
    (
        "flag",
        "uint8_zero",
        "accepted logical [1 1] false",
        "complex_zero",
        "MATLAB:validation:UnableToConvert",
    ),
    ("flag", "true", "accepted logical [1 1] true", "string", "MATLAB:validation:UnableToConvert"),
    (
        "flag",
        "false",
        "accepted logical [1 1] false",
        "cell_true",
        "MATLAB:validation:UnableToConvert",
    ),
    (
        "logical_row",
        "col_zero",
        "accepted logical [1 3] [false true true]",
        "m2x2",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "logical_row",
        "empty",
        "accepted logical [1 0] zeros(1,0)",
        "single_nan",
        "MATLAB:validation:UnableToConvert",
    ),
    (
        "logical",
        "m3x2",
        "accepted logical [3 2] [true true;true true;true true]",
        "handle",
        "MATLAB:validation:UnableToConvert",
    ),
    (
        "compression",
        "char_col_lzw",
        "accepted char [1 3] chars [108 122 119]",
        "upper_lzw",
        "MATLAB:validators:mustBeMember",
    ),
    (
        "compression",
        "string_lzw",
        "accepted char [1 3] chars [108 122 119]",
        "lzw_blank",
        "MATLAB:validators:mustBeMember",
    ),
    (
        "compression",
        "zstd",
        "accepted char [1 4] chars [122 115 116 100]",
        "char_matrix",
        "MATLAB:validation:IncompatibleSize",
    ),
    (
        "compression",
        "none",
        "accepted char [1 4] chars [110 111 110 101]",
        "char_empty",
        "MATLAB:validators:mustBeMember",
    ),
    (
        "mode",
        "parallel",
        "accepted char [1 8] chars [112 97 114 97 108 108 101 108]",
        "upper_parallel",
        "MATLAB:validators:mustBeMember",
    ),
    (
        "mode",
        "string_parallel",
        "accepted char [1 8] chars [112 97 114 97 108 108 101 108]",
        "char_nd",
        "MATLAB:validators:mustBeMember",
    ),
    (
        "char",
        "string_empty",
        "accepted char [0 0] chars []",
        "true",
        "MATLAB:validation:UnableToConvert",
    ),
    (
        "char",
        "char_matrix",
        "accepted char [2 2] chars [97 98;99 100]",
        "cell_mixed",
        "MATLAB:validation:UnableToConvert",
    ),
    (
        "char",
        "char_0x3",
        "accepted char [0 3] chars zeros(0,3)",
        "struct",
        "MATLAB:validation:UnableToConvert",
    ),
    (
        "char_p1",
        "char_1x0",
        "accepted char [1 0] chars zeros(1,0)",
        "complex",
        "MATLAB:validation:UnableToConvert",
    ),
]


@pytest.mark.parametrize(("key", "took", "want", "refused", "ident"), SIDES)
def test_each_rule_on_both_sides(key, took, want, refused, ident):
    got = outcomes(key, [took, refused])
    assert got == {took: want, refused: ident}


@pytest.mark.parametrize("key", list(DECLS))
def test_the_ports_refusals(key):
    """the conversions the port refuses: numbers, cells of char arrays and function handles for
    char; characters for logical; none for the declarations without a class"""
    cls = read(DECLS[key][0]).cls
    got = outcomes(key, [n for n, _, _ in BATTERY])
    refused = {n for n, o in got.items() if o == "refused"}
    assert refused == {n for n, g, _ in BATTERY if g in REFUSED_GROUPS.get(cls, ())}


def test_the_refusals_messages():
    """the port's own words, naming the argument, with its identifier"""
    p = REGISTRY.find("writetiff").param("filepath")
    for x, what in [
        (65.0, "not numbers"),
        (("a.tif",), "not a cell"),
        (math.sin, "not a function"),
    ]:
        with pytest.raises(ConversionRefusedError) as err:
            p.take(x, "writetiff", 2)
        assert err.value.identifier == "petakit5d:params:conversionRefused"
        assert str(err.value).startswith("The argument at position 2 ('filepath') must be text")
        assert what in str(err.value) and "docs/parity/io.md" in str(err.value)
    p = REGISTRY.find("writetiff").param("groupWrite")
    with pytest.raises(ConversionRefusedError) as err:
        p.take("a", "writetiff")
    assert str(err.value).startswith("The value of 'groupWrite' must be logical")
    assert isinstance(err.value, MatlabArgumentError) and isinstance(err.value, TypeError)


@pytest.mark.parametrize(
    ("function", "param", "x", "position", "message"),
    [
        (
            "readtiff",
            "range",
            (1, math.inf),
            None,
            "Invalid value for 'range' argument. Value must be finite.",
        ),
        (
            "readzarr",
            "inputBbox",
            np.ones((2, 3)),
            None,
            "Invalid value for 'inputBbox' argument. Value must be a vector.",
        ),
        (
            "createzarr",
            "clevel",
            (1, 2),
            None,
            "Invalid value for 'clevel' argument. Value must be a scalar.",
        ),
        (
            "createzarr",
            "clevel",
            "a",
            None,
            "Invalid value for 'clevel' argument. Value must be numeric.",
        ),
        (
            "readzarr",
            "sparseData",
            "a",
            None,
            "Invalid value for 'sparseData' argument. Value must be numeric or logical.",
        ),
        (
            "writetiff",
            "img",
            np.ones((2, 2), bool),
            1,
            "Invalid argument at position 1. Value must be numeric.",
        ),
        (
            "writetiff",
            "filepath",
            True,
            2,
            "Invalid argument at position 2. Value must be of type char or be convertible to char.",
        ),
        (
            "writetiff",
            "groupWrite",
            math.nan,
            None,
            "Invalid value for 'groupWrite' argument. Value must be of type logical or be "
            "convertible to logical.",
        ),
        (
            "writetiff",
            "Compression",
            "LZW",
            None,
            "Invalid value for 'Compression' argument. Value must be a member of this set:\n"
            "    'none'\n    'lzw'\n    'zstd'\n",
        ),
    ],
)
def test_matlabs_messages(function, param, x, position, message):
    """MATLAB's messages (measured with the lab's wrappers): 'Invalid value for '<name>'
    argument.' or 'Invalid argument at position <n>.', then the rule's sentence"""
    with pytest.raises(MatlabValidationError) as err:
        REGISTRY.find(function).param(param).take(x, function, position)
    assert str(err.value) == message
    assert isinstance(err.value, ValueError) and isinstance(err.value, MatlabArgumentError)


# normal calls: the value itself, unchanged
NORMAL = [
    ("readtiff", "filepath", "/mnt/projects/embryo_2024/a.tif"),
    ("readtiff", "range", [2, 5]),
    ("readtiff", "range", (2.0, 5.0)),
    ("readtiff", "range", np.array([2, 5])),
    ("readtiff", "range", None),
    ("writetiff", "Compression", "zstd"),
    ("writetiff", "Mode", "parallel"),
    ("writetiff", "groupWrite", False),
    ("writetiff", "groupWrite", np.True_),
    ("writetiff", "img", np.zeros((4, 5, 3), np.uint16)),
    ("createzarr", "dataSize", (500, 500, 500)),
    ("createzarr", "shardSize", ()),
    ("createzarr", "clevel", 1),
    ("createzarr", "dtype", "uint16"),
    ("writezarr", "create", True),
    ("writezarr", "sparseData", False),
    ("writezarr", "bbox", [1, 1, 1, 20, 11, 6]),
    ("readzarr", "sparseData", 1),
    ("readzarr", "inputBbox", mv.MArray("double", (1, 6), [1.0, 1.0, 1.0, 5.0, 5.0, 5.0])),
]


@pytest.mark.parametrize(("function", "param", "x"), NORMAL)
def test_normal_calls_are_unchanged(function, param, x):
    assert REGISTRY.find(function).param(param).take(x, function) is x


# what a converted value becomes for the function: (function, param, value, what it gets)
PASSED = [
    ("readtiff", "range", np.array([[2.0], [5.0]]), np.array([2.0, 5.0])),
    ("readtiff", "range", np.ones((1, 1, 2), np.int16), np.ones(2, np.int16)),
    (
        "readtiff",
        "range",
        mv.MArray("double", (2, 1), [2.0, 5.0]),
        mv.MArray("double", (1, 2), [2.0, 5.0]),
    ),
    ("readtiff", "range", ((2,), (5,)), (2, 5)),
    ("writezarr", "bbox", mv.MArray("double", (0, 0), []), mv.MArray("double", (1, 0), [])),
    ("writetiff", "Compression", np.array([["z"], ["s"], ["t"], ["d"]]), "zstd"),
    ("writetiff", "Compression", mv.MArray("char", (3, 1), list("lzw")), "lzw"),
    ("writetiff", "Compression", mv.MString("none"), "none"),
    ("createzarr", "dtype", mv.MString("uint8"), "uint8"),
    ("writetiff", "groupWrite", 1, True),
    ("writetiff", "groupWrite", 0.0, False),
    ("writetiff", "groupWrite", np.float32(-2.5), True),
    ("writetiff", "groupWrite", mv.MArray("uint8", (1, 1), [0]), False),
    ("writezarr", "sparseData", (1, 0, 2), (True, False, True)),
    ("writezarr", "sparseData", None, None),
    ("writezarr", "sparseData", np.array([[0.0], [3.0]]), np.array([False, True])),
    ("writezarr", "create", ((1, 0), (2, 3)), ((True, False), (True, True))),
    ("writezarr", "create", np.zeros((1, 1, 2)), np.zeros((1, 1, 2), bool)),
    (
        "writezarr",
        "create",
        mv.MArray("double", (2, 1), [1.0, 0.0]),
        mv.MArray("logical", (2, 1), [True, False]),
    ),
]


@pytest.mark.parametrize(("function", "param", "x", "want"), PASSED)
def test_the_converted_value_is_passed_on(function, param, x, want):
    got = REGISTRY.find(function).param(param).take(x, function)
    assert type(got) is type(want)
    if isinstance(want, np.ndarray):
        assert got.dtype == want.dtype and got.shape == want.shape and np.array_equal(got, want)
    else:
        assert got == want


def test_the_matlab_name_layer_passes_the_converted_values_on(tmp_path):
    """bind_matlab gives the functions MATLAB's converted values (in the Python forms of their
    kinds), and refuses what MATLAB refuses, with its errors"""
    spec = REGISTRY.find("createzarr")
    got = spec.bind_matlab(
        ["/mnt/projects/a.zarr", "dataSize", np.array([[20.0], [11.0], [6.0]]), "blockSize",
         np.full((1, 1, 3), 8.0), "groupWrite", 0, "expand2", 1.0, "shardSize", None]
    )  # fmt: skip
    assert got["data_size"].tolist() == [6.0, 20.0, 11.0]  # [y x z] -> (z, y, x)
    assert got["block_size"].tolist() == [8.0, 8.0, 8.0]
    assert got["group_write"] is False and got["expand_2d_dim"] is True
    assert got["shard_size"] is None
    got = REGISTRY.find("readtiff").bind_matlab(["a.tif", "range", ((3,), (5,))])
    assert got["range"] == (2, 5)
    got = REGISTRY.find("readzarr").bind_matlab(
        ["a.zarr"], {"inputBbox": mv.MArray("double", (6, 1), [2.0, 3.0, 4.0, 10.0, 11.0, 6.0])}
    )
    assert got["input_bbox"] == (3, 1, 2, 6, 10, 11)
    with pytest.raises(MatlabValidationError) as err:
        REGISTRY.find("readzarr").bind_matlab(["a.zarr", "inputBbox", [1, 1, 1, math.inf, 5, 5]])
    assert err.value.identifier == "MATLAB:validators:mustBeFinite"
    with pytest.raises(ConversionRefusedError):
        REGISTRY.find("writezarr").bind_matlab([np.zeros((2, 2)), "a.zarr", "dimSeparator", 47])
    with pytest.raises(MatlabValidationError) as err:
        REGISTRY.find("writetiff").bind_matlab([np.zeros((2, 2)), "a.tif", "comp", "LZW"])
    assert err.value.identifier == "MATLAB:validators:mustBeMember"


def test_the_positions_in_messages():
    """a positional argument's place among the function's positional ones, as MATLAB's
    messages name it; name-value arguments by their declared names"""
    with pytest.raises(MatlabValidationError, match="^Invalid argument at position 1. "):
        REGISTRY.find("writetiff").bind_matlab([np.zeros((2, 2), bool), "a.tif"])
    with pytest.raises(MatlabValidationError, match="^Invalid argument at position 2. "):
        REGISTRY.find("writezarr").bind_matlab([np.zeros((2, 2)), True])
    with pytest.raises(MatlabValidationError, match="^Invalid value for 'groupWrite' argument"):
        REGISTRY.find("writetiff").bind_matlab([np.zeros((2, 2)), "a.tif", "GROUP", (1, 1)])


def test_the_command_line_applies_the_declarations():
    """bind_strings: the declarations' rules on the converted strings (these functions have no
    mcc parser file: their strings stay text, which a numeric declaration refuses, and which the
    port refuses for a logical one)"""
    spec = REGISTRY.find("createzarr")
    assert spec.bind_strings(["a.zarr", "dtype", "uint8"]) == {
        "filepath": "a.zarr",
        "dtype": "uint8",
    }
    with pytest.raises(MatlabValidationError) as err:
        spec.bind_strings(["a.zarr", "dataSize", "[1 2 3]"])
    assert err.value.identifier == "MATLAB:validators:mustBeNumeric"
    with pytest.raises(ConversionRefusedError):
        spec.bind_strings(["a.zarr", "groupWrite", "true"])
    with pytest.raises(MatlabValidationError, match="^Invalid argument at position 1. "):
        REGISTRY.find("writetiff").bind_strings(["img", "a.tif"])


def test_the_snake_case_api_is_unchanged(tmp_path):
    """the snake_case API checks a declaration's classes as before: its values are the port's
    forms (an open range end None is MATLAB's Inf, which mustBeFinite would refuse)"""
    spec = REGISTRY.find("readtiff")
    assert spec.bind_python(("a.tif",), {"range": (2, None)}) == {
        "filepath": "a.tif",
        "range": (2, None),
    }
    spec = REGISTRY.find("writetiff")
    got = spec.bind_python(
        (np.zeros((2, 2), np.uint16), tmp_path / "a.tif"), {"compression": "LZW"}
    )
    assert got["compression"] == "LZW"  # (write_tiff refuses it itself)
    with pytest.raises(MatlabArgumentError) as err:
        spec.bind_python((np.zeros((2, 2), np.uint16), 65), {})
    assert err.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"


def _fake(text: str):
    reg = Registry()

    @reg.register("XR_decl_fake", Param("x", None), Param("y", text))
    def fake(x, y=None):
        return y

    return reg.find("XR_decl_fake")


def test_an_optional_argument_with_a_declaration():
    """an optional positional argument with a declaration is converted too (no registered one
    has it)"""
    spec = _fake("(1, :) {mustBeNumeric}")
    assert spec.bind_matlab([1, ((1,), (2,))]) == {"x": 1, "y": (1, 2)}
    with pytest.raises(MatlabValidationError, match="^Invalid argument at position 2. "):
        spec.bind_matlab([1, ((1, 2), (3, 4))])


def test_a_struct_array_and_other_values():
    d = read("(1, 1) {mustBeNumeric}")
    from petakit5d.io.matfile import StructArray

    rows = StructArray.from_list([{"a": 1}, {"a": 2}])
    with pytest.raises(MatlabValidationError, match="scalar"):
        d.convert(rows, "x")
    with pytest.raises(MatlabValidationError, match="numeric"):
        d.convert(object(), "x")
    with pytest.raises(MatlabValidationError, match="numeric"):
        d.convert(np.uint16, "x")  # a type is no MATLAB value
    assert read("char").convert(np.array(["abc"]), "x").size == (1, 3)  # NumPy text: its text
    assert value(10**400).array[0, 0] == math.inf
    assert value(np.float16(1)).kind == "other"
