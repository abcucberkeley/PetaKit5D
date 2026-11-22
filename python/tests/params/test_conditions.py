"""The lab's validator texts with conditions beyond class predicates, against MATLAB.

Until 2026-11-22 the registry checked a validator's class predicates only, so a text such as
``@(x) isvector(x) && numel(x) <= 2`` took a 1x3 vector that MATLAB refuses. Now a text with other
terms is read whole and evaluated as MATLAB evaluates it (``params.conditions``). ``LAB`` lists the
38 lab parameters whose validators have such terms (26 texts). MATLAB R2024b decides: each oracle
case runs inputParser with one parameter's own validator text (the lab's, as the registry and the
stored extract hold it) on a battery of values of every class and of many sizes, on both sides of
each condition, and records for each value 'accepted', inputParser's identifier, or
'raised(<identifier>)' for an error that the validator raises itself (inputParser passes it on
with that identifier). The Python side is the registry's check of that parameter on the same
values in the registry's forms. The checks follow inputParser only: a function may still fail
further on with a value they take, in MATLAB as in the port.

The battery holds no complex number and no string array of other than one element: the
registry's values have no form for them.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.params import MatlabArgumentError, MatlabCheckError, Param, Registry
from petakit5d.params import spec as S
from petakit5d.params import values as mv
from petakit5d.params.conditions import read
from petakit5d.params.registry import REGISTRY
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, mcell, mcol, mempty, mstring, mtyped
from petakit5d.testing.compare import MArray, MCell, MFunc, MString, MStruct, MValue, to_mvalue

# the lab's validator texts with other terms than class predicates, and their parameters
T_SAVE16 = "@(x) (numel(x) == 1 || numel(x) == 2) && islogical(x)"
T_OVERWRITE3 = "@(x) (numel(x) == 1 || numel(x) == 3) && islogical(x)"
T_OVERWRITE1 = "@(x) (numel(x) == 1) && islogical(x)"
T_RLMETHOD = "@(x) any(strcmp(x, {'original', 'simplified', 'omw', 'cudagen'}))"
T_INTERP = "@(x) any(strcmpi(x, {'cubic', 'linear'}))"
T_RESAMPLE = "@(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest', 'max', 'mean'}))"
T_MODE = "@(x) ischar(x) && any(strcmpi(x, {'valid', 'same', 'full'}))"
T_OVERLAP = "@(x) ischar(x) && ismember(lower(x), {'pixel', 'um'})"
T_STITCH_MIP = "@(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3))"
T_INTERVAL = "@(x) isempty(x) || (isnumeric(x) && size(x, 1) == 1 && size(x, 2) == 3)"
T_INDICES = "@(x) isempty(x) || (isnumeric(x) && size(x, 2) == 5)"
T_FLIPPED = "@(x) isempty(x) || all(islogical(x) | isnumeric(x))"
T_SAVE3D = "@(x) islogical(x) && (numel(x) == 1 || numel(x) == 3)"
T_XCORR_LOWER = "@(x) ismember(lower(x), {'primary', 'primaryfirst', 'all'})"
T_XCORR = "@(x) ismember(x, {'primary', 'primaryFirst', 'all'})"
T_OUT_BBOX = "@(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6)"
T_DATA_SIZE = "@(x) isnumeric(x) && (isempty(x) || numel(x) == 1)"
T_IM_SIZE = "@(x) isnumeric(x) && (isempty(x) || numel(x) == 3)"
T_TWO = "@(x) isnumeric(x) && numel(x) == 2"
T_THREE = "@(x) isnumeric(x) && numel(x) == 3"
T_PIXEL = "@(x) isvector(x) && numel(x) <= 2"
T_AXIS = "@(x) isvector(x) || numel(x) == 3"
T_LOGICAL1 = "@(x) numel(x) == 1 && islogical(x)"
T_NUMEL1 = "@(x) numel(x)==1"
T_MIXTURES = "@(x) numel(x)==1 && x>0 && round(x)==x"
T_NUMEL2 = "@(x) numel(x)==2"
LAB = [
    ("XR_microscopeAutomaticProcessing", "save16bit", "parameter", T_SAVE16),
    ("XR_microscopeAutomaticProcessing", "overwrite", "parameter", T_OVERWRITE3),
    ("XR_fftSpectrumComputingWrapper", "overwrite", "parameter", T_OVERWRITE1),
    ("RLdecon", "RLMethod", "parameter", T_RLMETHOD),
    ("deskewRotateFrame3D", "interpMethod", "parameter", T_INTERP),
    ("RLdecon", "interpMethod", "parameter", T_INTERP),
    ("XR_deskew_rotate_data_wrapper", "interpMethod", "parameter", T_INTERP),
    ("XR_deskewRotateBlock", "interpMethod", "parameter", T_INTERP),
    ("XR_deskewRotateFrame", "interpMethod", "parameter", T_INTERP),
    ("XR_deskewRotateZarr", "interpMethod", "parameter", T_INTERP),
    ("XR_resample_dataset", "interpMethod", "parameter", T_RESAMPLE),
    ("XR_resampleSingleZarr", "interpMethod", "parameter", T_RESAMPLE),
    ("XR_chromatic_shift_correction_data_wrapper", "mode", "parameter", T_MODE),
    ("XR_generate_image_list_wrapper", "overlapSizeType", "parameter", T_OVERLAP),
    ("XR_matlab_stitching_wrapper", "stitchMIP", "parameter", T_STITCH_MIP),
    ("XR_microscopeAutomaticProcessing", "stitchMIP", "parameter", T_STITCH_MIP),
    ("XR_stitching_frame_zarr_dev_v1", "stitchMIP", "parameter", T_STITCH_MIP),
    ("XR_generate_image_list_wrapper", "tileInterval", "parameter", T_INTERVAL),
    ("XR_generate_image_list_wrapper", "tileIndices", "parameter", T_INDICES),
    ("XR_stitching_frame_zarr_dev_v1", "flippedTile", "parameter", T_FLIPPED),
    ("RLdecon", "save3Dstack", "parameter", T_SAVE3D),
    ("XR_microscopeAutomaticProcessing", "xcorrMode", "parameter", T_XCORR_LOWER),
    ("XR_matlab_stitching_wrapper", "xcorrMode", "parameter", T_XCORR),
    ("XR_matlab_stitching_wrapper", "outBbox", "parameter", T_OUT_BBOX),
    ("XR_stitching_frame_zarr_dev_v1", "outBbox", "parameter", T_OUT_BBOX),
    ("XR_estimateComputingMemory", "dataSize", "parameter", T_DATA_SIZE),
    ("XR_estimateComputingMemory", "imSize", "parameter", T_IM_SIZE),
    ("omw_backprojector_generation", "hannWinBounds", "parameter", T_TWO),
    ("RLdecon", "mipAxis", "parameter", T_THREE),
    ("XR_stitching_frame_zarr_dev_v1", "stitchPadSize", "parameter", T_THREE),
    ("XR_decon_data_wrapper", "xyPixelSize", "parameter", T_PIXEL),
    ("XR_deskew_rotate_data_wrapper", "xyPixelSize", "parameter", T_PIXEL),
    ("XR_matlab_stitching_wrapper", "xyMaxOffset", "parameter", T_PIXEL),
    ("saveMIP_zarr", "axis", "parameter", T_AXIS),
    ("XR_decon_data_wrapper", "save16bit", "parameter", T_LOGICAL1),
    ("XR_runDetTrack3d", "TrackingGapLength", "parameter", T_NUMEL1),
    ("XR_runDetTrack3d", "MaxMixtures", "parameter", T_MIXTURES),
    ("XR_runDetTrack3d", "TrackingRadius", "parameter", T_NUMEL2),
]


def _chars(*rows: str) -> np.ndarray:
    """a char matrix (rows of one length)"""
    return np.array([list(r) for r in rows], dtype="<U1")


def _structs(size: tuple[int, ...], rows: list[dict]) -> MStruct:
    fields = list(rows[0]) if rows else []
    return MStruct(size, fields, [{k: to_mvalue(v) for k, v in r.items()} for r in rows])


def _cells(size: tuple[int, ...], items: list) -> MCell:
    return MCell(size, [x if isinstance(x, MValue) else to_mvalue(x) for x in items])


def _nd(shape: tuple[int, ...], mclass: str = "double") -> np.ndarray:
    a = np.arange(1, math.prod(shape) + 1, dtype=np.float64).reshape(shape, order="F")
    return a.astype(bool) if mclass == "logical" else a


# the values: classes and sizes on both sides of each condition (MATLAB's form)
BATTERY: dict[str, object] = {
    # logical
    "true": True,
    "false": False,
    "logical_row2": np.array([True, False]),
    "logical_row3": np.array([True, False, True]),
    "logical_col3": mcol([True, False, True], "logical"),
    "logical_2x2": np.ones((2, 2), dtype=bool),
    "logical_1x1x3": np.ones((1, 1, 3), dtype=bool),
    "logical_1x1x2": np.ones((1, 1, 2), dtype=bool),
    "logical_empty": mempty("logical"),
    "logical_1x0": mempty("logical", (1, 0)),
    "logical_0x1": mempty("logical", (0, 1)),
    # double scalars
    "zero": 0.0,
    "minus_zero": -0.0,
    "one": 1.0,
    "two": 2.0,
    "half": 2.5,
    "point5": 0.5,
    "minus_half": -2.5,
    "negative": -3.0,
    "tiny": 1e-300,
    "nan": math.nan,
    "inf": math.inf,
    "minus_inf": -math.inf,
    "big_odd": 2.0**52 + 1,  # MATLAB's round keeps it
    "big_half": 2.0**51 + 0.5,
    # empty arrays
    "empty": None,
    "empty_1x0": mempty("double", (1, 0)),
    "empty_0x1": mempty("double", (0, 1)),
    "empty_0x3": mempty("double", (0, 3)),
    "empty_1x1x0": mempty("double", (1, 1, 0)),
    # rows, columns, matrices, n-D arrays
    "row2": [1.0, 2.0],
    "row3": [1.0, 2.0, 3.0],
    "row4": [1.0, 2.0, 3.0, 4.0],
    "row5": [1.0, 2.0, 3.0, 4.0, 5.0],
    "row6": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    "row2_nan": [1.0, math.nan],
    "nan_row2": [math.nan, math.nan],
    "row3_inf": [1.0, 2.0, math.inf],
    "row2_half": [0.5, 1.0],
    "col2": mcol([1.0, 2.0]),
    "col3": mcol([1.0, 2.0, 3.0]),
    "col6": mcol([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
    "m2x2": _nd((2, 2)),
    "m3x2": _nd((3, 2)),
    "m2x3": _nd((2, 3)),
    "m2x5": _nd((2, 5)),
    "m3x5": _nd((3, 5)),
    "nd_1x1x2": _nd((1, 1, 2)),
    "nd_1x1x3": _nd((1, 1, 3)),
    "nd_1x2x3": _nd((1, 2, 3)),
    "nd_1x3x2": _nd((1, 3, 2)),
    "nd_3x2x2": _nd((3, 2, 2)),
    "nd_2x5x2": _nd((2, 5, 2)),
    "nd_3x2x1x2": _nd((3, 2, 1, 2)),
    # other numeric classes
    "uint8": mtyped(3, "uint8"),
    "uint8_zero": mtyped(0, "uint8"),
    "uint8_row": mtyped([1, 2], "uint8"),
    "int8_five": mtyped(5, "int8"),
    "int8_negative": mtyped(-2, "int8"),
    "int16_row": mtyped([1, 2], "int16"),
    "int16_row3": mtyped([1, 2, 3], "int16"),
    "uint16_row3": mtyped([1, 2, 3], "uint16"),
    "single": mtyped(1.5, "single"),
    "single_whole": mtyped(2.0, "single"),
    "single_row2": mtyped([1.0, 2.0], "single"),
    "single_nan": mtyped(math.nan, "single"),
    "single_inf": mtyped(math.inf, "single"),
    # char
    "char": "abc",
    "char_empty": "",
    "char_one": "a",
    "char_1x0": mempty("char", (1, 0)),
    "char_0x3": mempty("char", (0, 3)),
    "char_col": _chars("a", "b"),
    "char_col6": _chars(*"linear"),
    "char_matrix": _chars("ab", "cd"),
    "char_3x2": _chars("ab", "cd", "ef"),
    "char_3rows": _chars("a", "b", "c"),
    "char_pad": _chars("cubic ", "xxxxxx"),  # rows end in blanks
    "char_rows_match": _chars("cubic ", "linear"),
    "char_rows_full": _chars("valid", "same ", "full "),
    "char_rows_um": _chars("um", "ab"),
    "char_nd": np.full((1, 1, 2), "a", dtype="<U1"),
    # the words of the texts, their case, prefixes and blanks
    **{
        w: w
        for w in (
            "linear cubic nearest max mean valid same full pixel um primary primaryFirst "
            "primaryfirst all original simplified omw cudagen none lzw"
        ).split()
    },
    "upper_linear": "LINEAR",
    "mixed_cubic": "Cubic",
    "upper_mean": "MEAN",
    "prefix_lin": "lin",
    "linear_blank": "linear ",
    "upper_full": "FULL",
    "upper_um": "UM",
    "mixed_pixel": "Pixel",
    "upper_primary": "PRIMARY",
    "upper_all": "ALL",
    "mixed_original": "Original",
    # string scalars
    "string": mstring("abc"),
    "string_empty": mstring(""),
    "string_x": mstring("x"),
    "string_linear": mstring("linear"),
    "string_um": mstring("um"),
    "string_upper_um": mstring("UM"),
    "string_primary": mstring("primary"),
    "string_original": mstring("original"),
    # cells
    "cell_empty": MCell((0, 0), []),
    "cell_1x0": mcell(),
    "cell_char": mcell("a"),
    "cell_empty_char": mcell(""),
    "cell_linear": mcell("linear"),
    "cell_upper_um": mcell("UM"),
    "cell_um": mcell("um"),
    "cell_primary": mcell("primary"),
    "cell_upper_primary": mcell("PRIMARY"),
    "cell_original": mcell("original"),
    "cell_linear_cubic": mcell("linear", "cubic"),
    "cell_x_linear": mcell("x", "linear"),
    "cell_upper_x_linear": mcell("X", "LINEAR"),
    "cell_abc_def": mcell("abc", "def"),
    "cell_pixel_um": mcell("pixel", "um"),
    "cell_primary_all": mcell("primary", "all"),
    "cell_primary_x": mcell("primary", "x"),
    "cell_original_omw": mcell("original", "omw"),
    "cell_three": mcell("linear", "cubic", "abc"),
    "cell_x_same_y": mcell("x", "same", "y"),
    "cell_four": mcell("a", "b", "c", "d"),
    "cell_four_match": mcell("x", "simplified", "y", "z"),
    "cell_five_match": mcell("x", "y", "nearest", "z", "w"),
    "cell_col2": _cells((2, 1), ["linear", "cubic"]),
    "cell_2x2": _cells((2, 2), ["a", "c", "b", "d"]),
    "cell_1x1x2": _cells((1, 1, 2), ["linear", "cubic"]),
    "cell_num": mcell(1.0),
    "cell_nums": mcell(1.0, 2.0),
    "cell_true": mcell(True),
    "cell_string": mcell(mstring("linear")),
    "cell_strings": mcell(mstring("cubic"), mstring("linear")),
    "cell_num_linear": mcell(1.0, "linear"),
    "cell_linear_num": mcell("linear", 1.0),
    "cell_char_matrix": mcell(_chars("ab", "cd")),
    # structs and a function handle
    "struct": {"a": 1.0},
    "struct_row": _structs((1, 2), [{"a": 1.0}, {"a": 2.0}]),
    "struct_3": _structs((1, 3), [{"a": 1.0}, {"a": 2.0}, {"a": 3.0}]),
    "struct_empty": MStruct((0, 0), [], []),
    "struct_empty_field": MStruct((0, 0), ["a"], []),
    "handle": MFunc("@sin"),
}


def _py(v: MValue) -> object:
    """A MATLAB value of the battery in the registry's form: char rows as str ('' 0x0), other
    char arrays and numeric and logical arrays as the registry's MArray (MATLAB's class and
    size), cells as Cell, string scalars as MString, a 1x1 struct as a dict, other structs as
    StructArray, a handle as a callable"""
    from petakit5d.io.matfile import StructArray

    if isinstance(v, MArray):
        size = tuple(v.size)
        flat = v.data.reshape(-1, order="F")
        if v.mclass == "char":
            text = "".join(chr(int(c)) for c in flat)
            if size == (0, 0) or (len(size) == 2 and size[0] == 1 and size[1] > 0):
                return text
            return mv.MArray("char", size, list(text))
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
        return math.sin
    raise TypeError(f"no registry form for {v!r}")


BATTERY_PORT = {k: _py(to_mvalue(v)) for k, v in BATTERY.items()}

# inputParser with the parameter's validator text: 'accepted', inputParser's own identifier, or
# raised(<identifier>) for an error of the validator's own
_ADD = {"parameter": "ip.addParameter('x', [], %s);"}


def _parse_code(kind: str, validator: str) -> str:
    return (
        "ip = inputParser; " + _ADD[kind] % validator + " v = in{1}; n = in{2}; r = struct(); "
        "for i = 1 : numel(v); try; ip.parse('x', v{i}); r.(n{i}) = 'accepted'; "
        "catch ME; if startsWith(ME.identifier, 'MATLAB:InputParser:'); "
        "r.(n{i}) = ME.identifier; else; r.(n{i}) = ['raised(' ME.identifier ')']; end; end; end; "
        "out = {r};"
    )


def outcomes(function: str, param: str, values: dict) -> dict[str, str]:
    """the registry's check of one parameter on each value, as the MATLAB side records it"""
    p = REGISTRY.find(function).param(param)
    out = {}
    for name, v in values.items():
        try:
            p.check(v, function)
        except MatlabCheckError as err:
            out[name] = f"raised({err.identifier})"
        except MatlabArgumentError as err:
            out[name] = err.identifier
        else:
            out[name] = "accepted"
    return out


def _lab_case(function: str, param: str, kind: str, validator: str) -> Case:
    return Case(
        f"lab.{function}.{param}",
        code=_parse_code(kind, validator),
        args=[mcell(*BATTERY.values()), mcell(*BATTERY)],
        python=lambda *_: outcomes(function, param, BATTERY_PORT),
    )


CASES = [_lab_case(*row) for row in LAB]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- without MATLAB -------------------------------------------------------------------------------


def _class_terms_only(validator: str) -> bool:
    """the reading before 2026-11-22: whether every term of an anonymous text is a class
    predicate (its alternatives covered it whole)"""
    import re

    m = re.match(r"^@\s*\(\s*x\s*\)\s*(.+)$", validator.strip(), re.S)
    if not m:
        return True
    for alt in S._split_top(m.group(1), ("||", "|")):
        for term in S._split_top(S._strip_parens(alt), ("&&", "&")):
            term = S._strip_parens(term)
            pm = S._PRED_RE.match(term)
            am = S._ISA_RE.match(term)
            if not (pm and pm.group(1) in S._KNOWN_PREDS) and not (
                am and am.group(1) == "function_handle"
            ):
                return False
    return True


def test_lab_lists_every_parameter_whose_validator_has_a_condition():
    """the table is the registry's: each lab parameter whose validator has terms beyond class
    predicates, with its kind and text, is read whole (a new one needs a row and its MATLAB
    reference), and every other text is not"""
    found, unread = {}, []
    for s in REGISTRY.specs():
        for p in s.params:
            if p.condition is not None:
                found[(s.matlab_name, p.matlab)] = (p.kind, p.validator)
            elif p.validator and not _class_terms_only(p.validator):
                unread.append((s.matlab_name, p.matlab, p.validator))
    assert not unread, unread
    assert found == {(f, p): (k, v) for f, p, k, v in LAB}


@pytest.mark.parametrize("text", sorted({v for *_, v in LAB}))
def test_each_text_is_read(text):
    cond = read(text)
    assert cond is not None and cond.text == text
    assert read(text) is cond  # read once


@pytest.mark.parametrize(
    "text",
    [
        None,
        "",
        "@ischar",
        "@(x) ischar(x) || iscell(x)",
        "@(x) isempty(x) || isa(x,'function_handle') || ischar(x)",
        "@(x)islogical(x)",
        "(1, :) {mustBeNumeric}",
        "@(x) myCheck(x) && numel(x) == 2",  # a function this reading does not know
        "@(x) numel(x) + 1 == 2",  # arithmetic
        "@(x) strcmp(x, 'a')",  # strcmp of one text
        "@(x) ismember(x, [1 2])",
        "@(x) size(x, 0) == 1",
        "@(x) numel(x)' == 2",  # a transpose
        "@(x) numel(x) == 2 && x",  # && of the value itself
        "@(x) lower(x) == 2",  # lower() only inside ismember, as the lab's texts have it
        "@(x) numel(lower(x)) == 2",
        "@(x, y) numel(x) == 2",
        "@(x) numel(y) == 2",
        "@(x) numel(x) == 2)",
    ],
)
def test_texts_not_read(text):
    """class predicates only (their alternatives check them) and constructs this reading does
    not know give no condition"""
    assert read(text) is None


def test_matlabs_precedence():
    """MATLAB's operator precedence (its documentation): || below &&, below |, below &, below the
    comparisons; the lab's texts set parentheses wherever they mix them"""
    cond = read("@(x) isempty(x) || ischar(x) && numel(x) == 2")  # isempty || (ischar && ...)
    assert [cond.outcome(v)[0] for v in ("", "ab", "abc", 5.0)] == [True, True, False, False]
    cond = read("@(x) numel(x) == 1 | numel(x) == 2 & islogical(x)")  # == 1 | (== 2 & ...)
    assert [cond.outcome(v)[0] for v in (5.0, (True, False), (1.0, 2.0))] == [True, True, False]
    cond = read("@(x) numel(x) == 2 == 1")  # left to right: (numel == 2) == 1
    assert [cond.outcome(v)[0] for v in ("ab", "a")] == [True, False]


def test_a_quote_after_a_name_is_no_text():
    """MATLAB's transpose: a quote right after a name, a number or a closing bracket"""
    assert read("@(x) numel(x') == numel('ab')") is None
    assert read("@(x) numel(x) == numel(x) && any(strcmp(x, {'a'}))") is not None


def test_a_negation_is_read():
    """the class alternatives have no negation: ~ is a condition"""
    cond = read("@(x) ~isempty(x) && ischar(x)")
    assert cond is not None
    assert cond.outcome("a") == (True, None)
    assert cond.outcome("") == (False, None)


def test_the_class_alternatives_stay_for_hints_and_the_command_line():
    """the reading whole replaces the check only: a text's class alternatives, its command-line
    conversion and cell hint are as before"""
    p = REGISTRY.find("XR_decon_data_wrapper").param("xyPixelSize")
    assert p.alternatives == [frozenset({"isvector"})]
    assert p.condition is not None and p.mcc == "str2num"
    p = REGISTRY.find("XR_runDetTrack3d").param("TrackingRadius")
    assert p.alternatives is None and p.condition is not None


@pytest.mark.parametrize(("function", "param", "kind", "validator"), LAB)
def test_each_parameter_takes_both_sides_of_its_condition(function, param, kind, validator):
    """one value on each side of each parameter's condition (the cases compare all of them)"""
    took = outcomes(function, param, BATTERY_PORT)
    side = {
        T_SAVE16: ("logical_row2", "logical_row3"),
        T_OVERWRITE3: ("logical_col3", "logical_row2"),
        T_OVERWRITE1: ("true", "logical_row2"),
        T_RLMETHOD: ("cell_four_match", "cell_x_linear"),
        T_INTERP: ("cell_x_linear", "cell_linear_cubic"),
        T_RESAMPLE: ("upper_mean", "char_rows_match"),
        T_MODE: ("char_rows_full", "char_pad"),
        T_OVERLAP: ("upper_um", "char_rows_um"),
        T_STITCH_MIP: ("logical_1x1x3", "logical_2x2"),
        T_INTERVAL: ("nd_1x3x2", "col3"),
        T_INDICES: ("nd_2x5x2", "row4"),
        T_FLIPPED: ("nd_3x2x2", "char"),
        T_SAVE3D: ("logical_row3", "logical_row2"),
        T_XCORR_LOWER: ("cell_upper_primary", "cell_primary_all"),
        T_XCORR: ("primaryFirst", "primaryfirst"),
        T_OUT_BBOX: ("m3x2", "nd_1x2x3"),
        T_DATA_SIZE: ("empty_1x1x0", "row2"),
        T_IM_SIZE: ("nd_1x1x3", "row2"),
        T_TWO: ("nd_1x1x2", "row3"),
        T_THREE: ("col3", "row4"),
        T_PIXEL: ("col2", "row3"),
        T_AXIS: ("nd_1x1x3", "m2x2"),
        T_LOGICAL1: ("false", "logical_1x0"),
        T_NUMEL1: ("handle", "struct_row"),
        T_MIXTURES: ("big_odd", "big_half"),
        T_NUMEL2: ("char_nd", "char"),
    }[validator]
    assert took[side[0]] == "accepted"
    assert took[side[1]] != "accepted"


# MATLAB's errors (measured; the cases compare every value)
RAISED = [
    (T_OUT_BBOX, "nd_1x2x3", "MATLAB:sizeDimensionsMustMatch"),
    (T_OUT_BBOX, "nd_3x2x1x2", "MATLAB:sizeDimensionsMustMatch"),
    (T_INTERP, "cell_empty", "MATLAB:strcmp:InputsSizeMismatch"),
    (T_INTERP, "cell_col2", "MATLAB:strcmp:InputsSizeMismatch"),
    (T_INTERP, "char_3rows", "MATLAB:strcmp:RowsElementsMismatch"),
    (T_RLMETHOD, "char_col", "MATLAB:strcmp:RowsElementsMismatch"),
    (T_OVERLAP, "char_matrix", "MATLAB:nonLogicalConditional"),
    (T_OVERLAP, "char_nd", "MATLAB:cellstr:InputShape"),
    (T_XCORR, "one", "MATLAB:ISMEMBER:InputClass"),
    (T_XCORR, "cell_nums", "MATLAB:ISMEMBER:InputClass"),
    (T_XCORR, "handle", "MATLAB:ISMEMBER:InvalidInputsDataType"),
    (T_XCORR_LOWER, "cell_string", "MATLAB:lower:CellsMustContainChars"),
    (T_XCORR_LOWER, "struct", "MATLAB:ISMEMBER:InvalidInputsDataType"),
    (T_MIXTURES, "string", "MATLAB:string:ComparisonNotDefined"),
    (T_MIXTURES, "cell_num", "MATLAB:UndefinedFunction"),
    (T_MIXTURES, "struct", "MATLAB:UndefinedFunction"),
]


@pytest.mark.parametrize(("text", "value", "ident"), RAISED)
def test_errors_of_the_validators_own(text, value, ident):
    assert read(text).outcome(BATTERY_PORT[value]) == (False, ident)


def test_a_logical_array_is_a_refusal():
    """inputParser takes one true element only: ismember of a cell, even all true, refuses"""
    cond = read(T_XCORR)
    assert cond.outcome(mv.Cell(["primary", "all"])) == (False, None)
    assert cond.outcome(mv.Cell(["primary"])) == (True, None)
    assert cond.outcome(mv.Cell([], (0, 0))) == (False, None)


# the port's own forms (snake_case values, NumPy arrays): (text, value, outcome)
FORMS = [
    (T_PIXEL, (0.1, 0.2), (True, None)),
    (T_PIXEL, (0.1, 0.2, 0.3), (False, None)),
    (T_PIXEL, np.array([0.1, 0.2, 0.3]), (False, None)),
    (T_PIXEL, np.zeros((2, 1)), (True, None)),
    (T_PIXEL, np.zeros((2, 1, 1)), (True, None)),  # MATLAB's size drops trailing singletons
    (T_PIXEL, 0.1, (True, None)),
    (T_PIXEL, None, (False, None)),
    (T_THREE, (1, 0, 0), (True, None)),
    (T_THREE, np.array([1, 0, 0], dtype=np.uint8), (True, None)),
    (T_THREE, (True, False, False), (False, None)),
    (T_DATA_SIZE, None, (True, None)),
    (T_DATA_SIZE, np.float32(5.0), (True, None)),
    (T_OUT_BBOX, np.zeros((3, 2)), (True, None)),
    (T_OUT_BBOX, np.zeros((2, 3)), (True, None)),
    (T_OUT_BBOX, np.zeros((1, 2, 3)), (False, "MATLAB:sizeDimensionsMustMatch")),
    (T_OUT_BBOX, tuple(range(6)), (True, None)),
    (T_INTERVAL, (1.0, 2.0, 3.0), (True, None)),
    (T_INTERVAL, ((1.0, 2.0, 3.0),), (True, None)),
    (T_INDICES, np.zeros((4, 5)), (True, None)),
    (T_INDICES, np.zeros(5), (True, None)),
    (T_INTERP, "Linear", (True, None)),
    (T_INTERP, ("x", "linear"), (True, None)),  # element by element: 'linear' is the second
    (T_INTERP, ("x", "cubic"), (False, None)),
    (T_INTERP, ("cubic",), (True, None)),
    (T_INTERP, ("a", "b", "c"), (False, "MATLAB:strcmp:InputsSizeMismatch")),
    (T_INTERP, np.array(["l", "i", "n", "e", "a", "r"]), (True, None)),  # a char row
    (T_INTERP, 5, (False, None)),
    (T_RESAMPLE, mv.MString("linear"), (False, None)),  # ischar first
    (T_OVERLAP, "Pixel", (True, None)),
    (T_OVERLAP, ("um",), (False, None)),
    (T_XCORR, ("primary",), (True, None)),
    (T_XCORR, ("primary", "all"), (False, None)),
    (T_XCORR, 1, (False, "MATLAB:ISMEMBER:InputClass")),
    (T_XCORR, None, (False, "MATLAB:ISMEMBER:InputClass")),
    (T_XCORR, {"a": 1}, (False, "MATLAB:ISMEMBER:InvalidInputsDataType")),
    (T_XCORR, len, (False, "MATLAB:ISMEMBER:InvalidInputsDataType")),
    (T_XCORR_LOWER, ("Primary",), (True, None)),
    (T_XCORR_LOWER, (1,), (False, "MATLAB:ISMEMBER:InputClass")),  # (1,) is a numeric row
    (T_XCORR_LOWER, ("a", 1), (False, "MATLAB:lower:CellsMustContainChars")),
    (T_MIXTURES, 2, (True, None)),
    (T_MIXTURES, True, (True, None)),
    (T_MIXTURES, np.int16(3), (True, None)),
    (T_MIXTURES, 2.5, (False, None)),
    (T_MIXTURES, "a", (True, None)),  # 'a' is 97
    (T_MIXTURES, (2,), (True, None)),
    (T_MIXTURES, ("a",), (False, "MATLAB:UndefinedFunction")),
    (T_MIXTURES, mv.MString("a"), (False, "MATLAB:string:ComparisonNotDefined")),
    (T_MIXTURES, len, (False, "MATLAB:UndefinedFunction")),
    (T_NUMEL2, (3, 6), (True, None)),
    (T_NUMEL2, np.array([[3, 6]]), (True, None)),
    (T_NUMEL2, "ab", (True, None)),
    (T_NUMEL1, {"a": 1, "b": 2}, (True, None)),
    (T_FLIPPED, np.zeros((2, 2), dtype=bool), (True, None)),
    (T_FLIPPED, (), (True, None)),
    (T_STITCH_MIP, (True, False, True), (True, None)),
    (T_STITCH_MIP, np.array([True, False]), (False, None)),
    (T_SAVE16, np.array([True]), (True, None)),
]


@pytest.mark.parametrize(("text", "value", "want"), FORMS)
def test_the_ports_forms(text, value, want):
    assert read(text).outcome(value) == want


def test_check_raises_matlabs_identifiers():
    p = REGISTRY.find("XR_matlab_stitching_wrapper").param("xcorrMode")
    p.check("primary", "XR_matlab_stitching_wrapper")
    with pytest.raises(MatlabCheckError) as err:
        p.check(1.0, "XR_matlab_stitching_wrapper")
    assert err.value.identifier == "MATLAB:ISMEMBER:InputClass"
    assert "The function's validator (@(x)ismember(x,{'primary','primaryFirst','all'}))" in str(
        err.value
    )
    with pytest.raises(MatlabArgumentError) as err:
        p.check("Primary", "XR_matlab_stitching_wrapper")
    assert not isinstance(err.value, MatlabCheckError)
    assert err.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"


def test_the_snake_case_api_checks_matlabs_form():
    """bind_python checks a value after its kind's conversion, as MATLAB's validator sees it"""
    spec = REGISTRY.find("RLdecon")
    p = spec.param("mipAxis")
    assert p.conv is not None and p.conv.name == "flags_yxz"
    p.check(S._as_matlab_sees(p, (1, 0, 0)), "RLdecon")
    with pytest.raises(MatlabArgumentError):
        p.check(S._as_matlab_sees(p, (1, 0, 0, 0)), "RLdecon")
    p = REGISTRY.find("XR_decon_data_wrapper").param("xyPixelSize")
    with pytest.raises(MatlabArgumentError):
        p.check(S._as_matlab_sees(p, (0.108, 0.108, 0.5)), "XR_decon_data_wrapper")


def _spec(validator: str, *, optional: bool = False):
    reg = Registry()
    params = (Param("x", validator), Param("Flag", "@islogical"))
    if optional:

        @reg.register("XR_conditions_fake", *params)
        def fake_opt(x=None, *, flag=False):
            return x

    else:

        @reg.register("XR_conditions_fake", *params)
        def fake(x, *, flag=False):
            return x

    return reg.find("XR_conditions_fake")


def test_an_optional_argument_with_a_condition():
    """inputParser at an optional argument: a value its validator takes is the argument, one it
    refuses fails (a name followed by a value starts the name-value part), one that makes it
    raise fails with that error"""
    spec = _spec(T_XCORR, optional=True)
    assert spec.bind_matlab(["all"]) == {"x": "all"}
    assert spec.bind_matlab(["Flag", True]) == {"flag": True}
    with pytest.raises(MatlabArgumentError) as err:
        spec.bind_matlab(["xyz"])
    assert err.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
    with pytest.raises(MatlabCheckError) as err:
        spec.bind_matlab([1.0])
    assert err.value.identifier == "MATLAB:ISMEMBER:InputClass"
    spec = _spec(T_PIXEL, optional=True)
    assert spec.bind_matlab([(0.1, 0.2), "Flag", True]) == {"x": (0.1, 0.2), "flag": True}
    with pytest.raises(MatlabArgumentError):
        spec.bind_matlab([(0.1, 0.2, 0.3)])


def test_the_three_ways_in_check_the_condition():
    spec = _spec(T_PIXEL)
    with pytest.raises(MatlabArgumentError):
        spec.bind_python(((0.1, 0.2, 0.3),), {})
    with pytest.raises(MatlabArgumentError):
        spec.bind_matlab([[0.1, 0.2, 0.3]])
    with pytest.raises(MatlabArgumentError):
        spec.bind_strings(["[0.1 0.2 0.3]"])
    assert spec.bind_strings(["[0.1 0.2]"]) == {"x": (0.1, 0.2)}
