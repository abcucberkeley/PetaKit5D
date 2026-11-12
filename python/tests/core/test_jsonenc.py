import math

import numpy as np
import pytest

from petakit5d.core.jsonenc import format_json_number, jsonencode


@pytest.mark.parametrize(
    "x,expected",
    [
        (0.1 + 0.2, "0.30000000000000004"),
        (1 / 3, "0.33333333333333331"),
        (math.pi, "3.1415926535897931"),
        (1e20, "1.0E+20"),
        (1e15, "1.0E+15"),
        (123456789012345678, "1.2345678901234568E+17"),
        (2.0**53, "9.007199254740992E+15"),
        (1e-7, "1E-7"),
        (1e-5, "1E-5"),
        (1.5e-12, "1.5E-12"),
        (1e-4, "0.0001"),
        (999999.0, "999999"),
        (1e6, "1.0E+6"),
        (-0.0, "-0"),
        (0.5, "0.5"),
        (5e-324, "4.94065645841247E-324"),
        (100.0, "100"),
    ],
)
def test_double_format(x, expected):
    assert format_json_number(x) == expected


@pytest.mark.parametrize(
    "x,expected",
    [
        (0.1, "0.1"),
        (1 / 3, "0.333333343"),
        (1e-12, "1E-12"),
        (1234567.88, "1.23456788E+6"),
        (-0.999999881, "-0.999999881"),
    ],
)
def test_single_format(x, expected):
    assert format_json_number(np.float32(x), single=True) == expected


def test_values():
    assert jsonencode(np.array([1.0, 2.0, 3.0])) == "[1,2,3]"
    assert jsonencode(np.array([[1.0], [2.0], [3.0]])) == "[1,2,3]"  # vectors are flat
    assert jsonencode(np.array([[1, 2], [3, 4]], dtype=float)) == "[[1,2],[3,4]]"
    assert (
        jsonencode(np.arange(1.0, 9.0).reshape((2, 2, 2), order="F"))
        == "[[[1,5],[3,7]],[[2,6],[4,8]]]"
    )
    assert jsonencode(np.array([5.0])) == "5" and jsonencode([5]) == "[5]"
    assert (
        jsonencode(None) == "[]"
        and jsonencode([]) == "[]"
        and jsonencode("") == '""'
        and jsonencode({}) == "{}"
    )
    assert jsonencode([math.nan, math.inf]) == "[null,null]"
    assert jsonencode([math.nan, -math.inf], convert_inf_and_nan=False) == "[NaN,-Infinity]"
    assert jsonencode(np.uint64(2**64 - 1)) == "18446744073709551615"
    assert jsonencode(np.array([True, False])) == "[true,false]"
    assert jsonencode('q"\\/\té\x01') == '"q\\"\\\\/\\té\\u0001"'


def test_pretty_print():
    s = jsonencode(
        {"a": 1, "c": np.array([1.0, 2.0]), "e": None, "f": {"g": True}, "h": []}, pretty_print=True
    )
    assert (
        s == '{\n  "a": 1,\n  "c": [\n    1,\n    2\n  ],\n  "e": [],\n'
        '  "f": {\n    "g": true\n  },\n  "h": []\n}'
    )
    assert (
        jsonencode([{"a": 1}, {"a": 2}], pretty_print=True)
        == '[\n  {\n    "a": 1\n  },\n  {\n    "a": 2\n  }\n]'
    )


def test_rejects_objects():
    with pytest.raises(TypeError):
        jsonencode(object())
