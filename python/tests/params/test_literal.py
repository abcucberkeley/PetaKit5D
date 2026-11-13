"""The MATLAB literal parser and formatter (unit tests; MATLAB's answers are in
test_literal_oracle.py)."""

from __future__ import annotations

import math

import pytest

from petakit5d.params import Cell, LiteralError, MArray, MString, format_literal, parse, str2num
from petakit5d.params.values import describe, from_python, same, to_python


@pytest.mark.parametrize(
    ("text", "py"),
    [
        ("1", 1),
        ("-32.45", -32.45),
        ("1e-3", 0.001),
        (".5", 0.5),
        ("true", True),
        ("false", False),
        ("[]", None),
        ("{}", ()),
        ("''", ""),
        ("'it''s'", "it's"),
        ('"a""b"', 'a"b'),
        ("[1,2,3]", (1, 2, 3)),
        ("[1 2 3]", (1, 2, 3)),
        ("[1 -2]", (1, -2)),
        ("[1 - 2]", -1),
        ("[1 , -2]", (1, -2)),
        ("[1;2]", (1, 2)),
        ("[1,2;3,4]", ((1, 2), (3, 4))),
        ("[1 2\n3 4]", ((1, 2), (3, 4))),
        ("[1, 2, ...\n 3]", (1, 2, 3)),
        ("{'a','b'}", ("a", "b")),
        ("{'CamA_ch0', 'CamB_ch0'}", ("CamA_ch0", "CamB_ch0")),
        ("{{'a'},{'b','c'}}", (("a",), ("b", "c"))),
        ("{'', '', ''}", ("", "", "")),
        ("{1, 'a'; [1 2], {}}", ((1, "a"), ((1, 2), ()))),
        ("[true, false]", (True, False)),
        ("[true, 1]", (1, 1)),
        ("0:3", (0, 1, 2, 3)),
        ("1:0", None),
        ("[0:2:6]'", (0, 2, 4, 6)),
        ("pi/2", math.pi / 2),
        ("2^-1", 0.5),
        ("-2^2", -4),
        ("[1,2]*2", (2, 4)),
        ("zeros(1,0)", None),
        ("ones(1,2)", (1, 1)),
        ("cell(1,2)", (None, None)),
        ("['ab', 'cd']", "abcd"),
        ("['ab'; 'cd']", ("ab", "cd")),
        ("['' 65]", "A"),
        ("[{1}, 2]", (1, 2)),
        ("uint16([1.5 2.5 -1])", (2, 3, 0)),
        ("int8(-2.5)", -3),
    ],
)
def test_parse(text, py):
    got = to_python(parse(text))
    assert got == py and type(got) is type(py)


@pytest.mark.parametrize(
    "text",
    ["[1,2;3]", "[1 2", "{1 2", "'abc", "1i", "x", "a.b", "foo(1)", "[1 2]*[3 4]", "", "1 +"],
)
def test_parse_rejects(text):
    with pytest.raises(LiteralError):
        parse(text)


def test_classes_and_sizes():
    assert parse("[1;2;3]").shape == (3, 1)
    assert parse("[]").shape == (0, 0)
    assert parse("zeros(0,3)").shape == (0, 3)
    assert parse("{'a';'b'}").shape == (2, 1)
    assert parse("single(1)").mclass == "single"
    assert parse("[single(1), 2]").mclass == "single"
    assert parse("[uint8(1), 2.7]").data == [1, 3]
    assert parse("logical([1 0 2])").data == [True, False, True]
    assert isinstance(parse('"s"'), MString)
    assert isinstance(parse("{1}"), Cell)
    assert parse("''").shape == (0, 0)


# MATLAB R2024b's answers (eval; probes of task X01, 2026-11-13). The literal corpus holds more,
# up to 2^53 (test_literal_oracle.py); these include the 64-bit ends.
@pytest.mark.parametrize(
    ("text", "mclass", "data"),
    [
        ("0x100", "uint16", [256]),
        ("0x0001", "uint16", [1]),  # the digits' bits decide the class, not the value
        ("0x000000FF", "uint32", [255]),
        ("0xFFFFFFFFFFFFFFFF", "uint64", [2**64 - 1]),
        ("0x8000000000000000", "uint64", [2**63]),
        ("0b" + "1" * 64, "uint64", [2**64 - 1]),
        ("0b100000000", "uint16", [256]),
        ("0X1F", "uint8", [31]),
        ("0x1e3", "uint16", [483]),
        ("0xFFFFFFFFFFFFFFFFu64", "uint64", [2**64 - 1]),
        ("0x8000000000000000s64", "int64", [-(2**63)]),
        ("0x7FFFFFFFFFFFFFFFs64", "int64", [2**63 - 1]),
        ("0xFFFFFFFFFFFFFFFFs64", "int64", [-1]),
        ("0xFFs16", "int16", [255]),
        ("0xFFS8", "int8", [-1]),
        ("0b1s64", "int64", [1]),
        ("-0xFFFFFFFFFFFFFFFF", "uint64", [0]),
        ("--0x10", "uint8", [0]),
        ("-0x80s8", "int8", [127]),
        ("int8(-0x1s16)", "int8", [-1]),
        ("[uint64(5), 0xFFFFFFFFFFFFFFFF]", "uint64", [5, 2**64 - 1]),
        ("[0xFFFFFFFFFFFFFFFF, int8(-1)]", "uint64", [2**64 - 1, 0]),
        ("[int64(-1), uint64(18446744073709551615)]", "int64", [-1, 2**63 - 1]),
        ("uint64(0xFFFFFFFFFFFFFFFFs64)", "uint64", [0]),
        ("int64(0xFFFFFFFFFFFFFFFF)", "int64", [2**63 - 1]),
        ("-int64(-9223372036854775808)", "int64", [2**63 - 1]),
        ("[0x41 'a']", "char", ["A", "a"]),
        ("~0x10", "logical", [False]),
    ],
)
def test_hex_and_binary_literals(text, mclass, data):
    v = parse(text)
    assert (v.mclass, v.data) == (mclass, data)
    assert same(parse(format_literal(v)), v)


@pytest.mark.parametrize(
    "text",
    [
        "0x10000000000000000",  # 17 digits: more bits than uint64 has
        "0x00000000000000000",
        "0b1" + "0" * 64,
        "0x00FFu8",  # the digits' bits, not the value, must fit the suffix
        "0x0080s8",
        "0b000000001u8",
        "0x1FFFFFFFFFFFFFFFFu64",
        "0xFFu", "0xFFs", "0x1u128", "0xFFi8", "0x1us8", "0x", "0b", "0xG", "0b2", "0b102",
        "0x1.5", "0x1p3", "0x1F.", "0x10i", "0xFF_FF", "00x10", "0b1e3", "1x10",
    ],
)  # fmt: skip
def test_hex_and_binary_literals_matlab_rejects(text):
    with pytest.raises(LiteralError):
        parse(text)


@pytest.mark.parametrize("text", ["0x10 + 1", "2 * 0x10", "0x10:0x12", "abs(0x10)"])
def test_integer_arithmetic_is_not_supported(text):
    """valid MATLAB (uint8 17, 32, [16 17 18], 16) that the parser refuses"""
    with pytest.raises(LiteralError):
        parse(text)


def test_str2num_reads_hex_literals():
    x, ok, _ = str2num("512, 0x100")
    assert ok and (x.mclass, x.data) == ("uint16", [512, 256])


def test_nan_has_matlabs_bits():
    assert describe(parse("NaN"))["d"] == ["fff8000000000000"]
    assert describe(parse("-NaN"))["d"] == ["7ff8000000000000"]
    assert describe(parse("0/0"))["d"] == ["fff8000000000000"]


def test_colon_is_symmetric_like_matlab():
    v = parse("0:0.1:1").data
    assert len(v) == 11 and v[-1] == 1.0 and v[3] == 0.30000000000000004
    assert parse("1:0").shape == (1, 0)
    assert parse("10:-3:1").data == [10.0, 7.0, 4.0, 1.0]


def test_str2num():
    assert str2num("1 2 3")[0].data == [1.0, 2.0, 3.0]
    assert str2num("true")[0].mclass == "logical"
    for text in ("yes", "'abc'", "{1}", ""):
        x, ok, _ = str2num(text)
        assert x.shape == (0, 0) and not ok


@pytest.mark.parametrize(
    "value",
    [
        1,
        -0.0,
        0.1,
        1e-310,
        1.7976931348623157e308,
        math.inf,
        -math.inf,
        True,
        "",
        "it's",
        "a\nb",
        (1, 2, 3),
        ((1, 2), (3, 4)),
        (True, False),
        ("a", "b"),
        ("a", 1, (1, 2), ()),
        None,
        MString('q"'),
        Cell(),
        MArray("single", (1, 2), [0.5, 0.10000000149011612]),
        MArray("uint16", (1, 1), [7]),
        MArray("double", (0, 3), []),
        MArray("logical", (0, 0), []),
        MArray("char", (2, 2), list("acbd")),
        Cell([MArray.scalar("double", 1.0), MArray.char("x")], (2, 1)),
    ],
)
def test_format_round_trip(value):
    hint = "cell" if isinstance(value, tuple) and any(isinstance(e, str) for e in value) else None
    m = from_python(value, hint)
    assert same(parse(format_literal(m)), m), format_literal(m)


def test_format_texts():
    assert format_literal(("CamA_ch0", "CamB_ch0")) == "{'CamA_ch0','CamB_ch0'}"
    assert format_literal((1024, 1024, 1024)) == "[1024,1024,1024]"
    assert format_literal(None) == "[]"
    assert format_literal((), "cell") == "{}"
    assert format_literal(0.108) == "0.108"
    assert format_literal(False) == "false"
    assert format_literal(math.nan) == "-NaN"  # Python's NaN has the sign bit clear
    assert format_literal(-math.nan) == "NaN"


def test_python_matrices_are_rows():
    assert same(from_python(((1, 2), (3, 4))), parse("[1,2;3,4]"))
    assert same(from_python(((True, False),)), parse("[true,false]"))
    assert format_literal(((1, 2, 3), (4, 5, 6))) == "[1,2,3;4,5,6]"
    assert same(from_python(("a", ("b",)), "cell"), parse("{'a',{'b'}}"))


def test_parse_takes_only_text():
    with pytest.raises(TypeError):
        parse(1)  # type: ignore[arg-type]
