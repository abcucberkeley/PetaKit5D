"""MATLAB's isstruct in the parameter registry, in the port's forms of MATLAB values.

A dict (a 1x1 struct), an ``io.matfile.StructArray`` of any size (empty ones, with fields or
without, included) and the classes that hold MATLAB structs (bwconncomp's ``ConnComp``, a
condition's movies) are structs; objects (``io.matfile.MatlabObject``), cells (lists and
tuples, of structs too) and all other values are not. MATLAB R2024b decides the same values
(oracle cases): ``isstruct`` itself, and inputParser's ``@isstruct`` validator, which the
registry's two ways in (the MATLAB-name layer and the snake_case API) must follow, so that
registered functions take struct arrays (a condition of several movies).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.io.matfile import MatlabObject, StructArray
from petakit5d.params import MatlabArgumentError, MString, Param, Registry
from petakit5d.params.spec import _pred
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, mstring, mtyped
from petakit5d.testing.compare import MCell, MStruct, to_mvalue


def isstruct(x) -> bool:
    """the registry's isstruct (the Python side of the cases)"""
    return _pred("isstruct", x)


def mstruct_array(size, fields, rows) -> MStruct:
    """MATLAB's struct array of ``rows`` (Python values, column-major) as the harness passes it"""
    return MStruct(size, list(fields), [{k: to_mvalue(v) for k, v in r.items()} for r in rows])


ROWS3 = [{"a": 1.0}, {"a": 2.0}, {"a": 3.0}]
# struct arrays: name -> (size, fields, rows), MATLAB's struct array and the port's StructArray
ARRAYS = {
    "empty_0x0": ((0, 0), (), []),
    "empty_0x0_fields": ((0, 0), ("b", "a"), []),
    "empty_0x3_fields": ((0, 3), ("a",), []),
    "empty_3x0_fields": ((3, 0), ("a",), []),
    "empty_0x2x3_fields": ((0, 2, 3), ("a",), []),
    "empty_0x3": ((0, 3), (), []),
    "no_fields_1x3": ((1, 3), (), [{}, {}, {}]),
    "row_1x3": ((1, 3), ("a",), ROWS3),
    "column_3x1": ((3, 1), ("a",), ROWS3),
    "matrix_3x2": ((3, 2), ("b", "a"), [{"b": float(i), "a": "x"} for i in range(6)]),
    "nd_2x2x2": ((2, 2, 2), ("a",), [{"a": float(i)} for i in range(8)]),
}
# other values: name -> (MATLAB's argument, the port's value)
VALUES = {
    "scalar_no_fields": ({}, {}),
    "scalar": ({"a": 1.0}, {"a": 1.0}),
    "scalar_with_empty_field": (
        MStruct((1, 1), ["t"], [{"t": MStruct((0, 0), ["a"], [])}]),
        {"t": StructArray.empty(("a",))},
    ),
    "cell_of_struct": ([{"a": 1.0}], [{"a": 1.0}]),
    "cell_of_empty_struct": (
        MCell((1, 1), [MStruct((0, 0), ["a"], [])]),
        [StructArray.empty(("a",))],
    ),
    "empty_cell": (MCell((0, 0), []), []),
    "empty_double": (None, None),
    "char": ("abc", "abc"),
    "string": (mstring("abc"), MString("abc")),
    "double": (1.0, 1.0),
    "logical": (True, True),
    "int8": (mtyped(1, "int8"), np.int8(1)),
}


def _conncomp():
    from petakit5d.matlabcompat.label import bwconncomp

    return bwconncomp(np.ones((3, 4), dtype=bool))


# MATLAB objects and handles: name -> (MATLAB expression, the port's value of it)
OBJECTS = {
    "function_handle": ("@sin", lambda: math.sin),
    "inputParser": ("inputParser", lambda: MatlabObject("inputParser")),
    "inputParser_empty": ("inputParser.empty", lambda: MatlabObject("inputParser")),
    "containers_Map": ("containers.Map({'a'}, {1})", lambda: MatlabObject("containers.Map")),
    "datetime": ("datetime(2026, 1, 1)", lambda: MatlabObject("datetime")),
    "table": ("table([1; 2])", lambda: MatlabObject("table")),
    "dictionary": ("dictionary", lambda: MatlabObject("dictionary")),
    "struct_of_object": ("struct('p', inputParser)", lambda: {"p": MatlabObject("inputParser")}),
    "bwconncomp": ("bwconncomp(true(3, 4))", _conncomp),
}

ISSTRUCT = (
    [
        Case(
            f"isstruct.{name}",
            matlab="isstruct",
            args=[mstruct_array(size, fields, rows)],
            python=isstruct,
            py_args=[StructArray(size, fields, rows)],
        )
        for name, (size, fields, rows) in ARRAYS.items()
    ]
    + [
        Case(f"isstruct.{name}", matlab="isstruct", args=[m], python=isstruct, py_args=[p])
        for name, (m, p) in VALUES.items()
    ]
    + [
        Case(
            f"isstruct.{name}",
            code=f"out = {{isstruct({expr})}};",
            python=lambda make=make: isstruct(make()),
        )
        for name, (expr, make) in OBJECTS.items()
    ]
)

# inputParser's @isstruct validator (addRequired), and a registered function with the same one
PARSE = "ip = inputParser; ip.addRequired('data', @isstruct); ip.parse({}); out = {{}};"
REFUSED = "MATLAB:InputParser:ArgumentFailedValidation"
REG = Registry()


@REG.register("pkpy_isstruct_parse", Param("data", "@isstruct", "A struct."))
def takes_a_struct(data):
    return None


def parse(x) -> None:
    """``x`` through the registry's two ways in, which must agree; their error if they refuse"""
    spec = REG.find("pkpy_isstruct_parse")
    errors = []
    for call in (spec.call_matlab, spec.call_python):
        try:
            call(x)
            errors.append(None)
        except MatlabArgumentError as err:
            errors.append(err)
    if (errors[0] is None) != (errors[1] is None):
        raise AssertionError(f"the MATLAB-name layer and the snake_case API disagree: {errors}")
    if errors[0] is not None:
        raise errors[0]


PARSED = {  # name -> (MATLAB's argument, the port's value, refused)
    **{
        n: (mstruct_array(*ARRAYS[n]), StructArray(*ARRAYS[n]), False)
        for n in ("row_1x3", "empty_0x0_fields", "empty_0x0", "no_fields_1x3")
    },
    "scalar": (*VALUES["scalar"], False),
    "cell_of_struct": (*VALUES["cell_of_struct"], True),
    "empty_double": (*VALUES["empty_double"], True),
}
PARSE_CASES = [
    Case(
        f"parse.{name}",
        code=PARSE.format("in{1}"),
        args=[m],
        python=parse,
        py_args=[p],
        nargout=0,
        expect_error=REFUSED if refused else False,
        py_error=MatlabArgumentError,
        error_message=True,
    )
    for name, (m, p, refused) in PARSED.items()
] + [
    Case(
        "parse.inputParser",
        code=PARSE.format("inputParser"),
        python=lambda: parse(MatlabObject("inputParser")),
        nargout=0,
        expect_error=REFUSED,
        py_error=MatlabArgumentError,
        error_message=True,
    )
]

CASES = ISSTRUCT + PARSE_CASES


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- the port's forms, without MATLAB -------------------------------------------------------------


def movie2d(source: str):
    from petakit5d.detection.condition import Movie2D

    return Movie2D(source, ["ch1/"], None, 1.0, (3, 4, 5), 2, ["gfp"], [["f.tif"]], ["m.tif"],
                   1.4, 1.0, 0.104)  # fmt: skip


def movie(source: str):
    from petakit5d.detection.condition import Movie

    return Movie(source, ["ch1/"], None, 1.0, (3, 4, 5), 2, ["gfp"], [["f.tif"]], ["m.tif"])


@pytest.mark.parametrize("name", list(ARRAYS))
def test_struct_arrays_of_any_size_are_structs(name):
    assert isstruct(StructArray(*ARRAYS[name]))


def test_struct_arrays_of_one_element_and_dicts_are_structs():
    assert isstruct(StructArray((1, 1), ("a",), [{"a": 1.0}]))
    assert isstruct({}) and isstruct({"a": 1.0})
    import types

    assert isstruct(types.MappingProxyType({"a": 1.0}))  # another mapping


def test_the_classes_that_hold_matlab_structs_are_structs():
    from petakit5d.detection.condition import Condition

    assert isstruct(_conncomp())  # bwconncomp's struct
    assert isstruct(movie2d("m/")) and isstruct(movie("m/"))  # Movie is a Movie2D
    assert isstruct(Condition([movie("a/"), movie("b/")]))


def test_the_classes_are_known_by_module_and_name():
    # the rule imports nothing: a class counts by its module and name (and its subclasses)
    fake = type("Condition", (list,), {"__module__": "petakit5d.detection.condition"})
    assert isstruct(fake()) and isstruct(type("Sub", (fake,), {})())
    other = type("Condition", (list,), {"__module__": "elsewhere"})
    assert not isstruct(other())


NOT_STRUCTS = {
    "object": MatlabObject("inputParser"),
    "object_named_struct": MatlabObject("struct"),  # an object, whatever its class name says
    "cell_of_struct": [{"a": 1.0}],
    "cell_of_struct_tuple": ({"a": 1.0},),
    "cell_of_empty_struct": [StructArray.empty(("a",))],
    "empty_cell": [],
    "empty_tuple": (),
    "empty_double": None,
    "char": "abc",
    "string": MString("abc"),
    "double": 1.0,
    "logical": True,
    "int8": np.int8(1),
    "empty_array": np.zeros((0, 0)),
    "numpy_structured": np.zeros(1, dtype=[("a", "f8")]),  # not the port's form of a struct
    "function_handle": math.sin,
    "the_class_itself": StructArray,
}


@pytest.mark.parametrize("name", list(NOT_STRUCTS))
def test_objects_cells_and_other_values_are_not_structs(name):
    assert not isstruct(NOT_STRUCTS[name])


def test_validators_with_isstruct_take_struct_arrays():
    from petakit5d.params.spec import validator_alternatives

    empty_or_struct = validator_alternatives("@(x) isempty(x) || isstruct(x)")
    sa = StructArray(*ARRAYS["row_1x3"])
    assert any(all(_pred(q, sa) for q in alt) for alt in empty_or_struct)
    assert not _pred("iscell", sa) and not _pred("isempty", sa)
    assert _pred("isempty", StructArray.empty(("a",))) and isstruct(StructArray.empty(("a",)))


def test_registered_functions_take_a_condition_of_several_movies():
    # the gap: XR_runTracking3D (and every @isstruct parameter) refused a 1 x n struct array
    from petakit5d.detection.condition import struct_array
    from petakit5d.params import find

    data = struct_array([movie("a/"), movie("b/")])
    assert isinstance(data, StructArray) and data.shape == (1, 2)
    for name in ("XR_runTracking3D", "XR_runDetTrack3d"):
        spec = find(name)
        assert spec.bind_python((data,), {})["data"] is data
    assert find("XR_runTracking3D").bind_matlab([data])["data"] is data
    got = find("XR_runDetTrack3d").bind_matlab([data])["data"]
    assert [m.source for m in got] == ["a/", "b/"]
    with pytest.raises(MatlabArgumentError) as err:
        find("XR_runTracking3D").bind_python(([m.to_matlab() for m in got],), {})
    assert err.value.identifier == REFUSED  # a cell of structs, as in MATLAB
