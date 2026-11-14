import math

import numpy as np
import pytest

from petakit5d.core import dtypes, text


@pytest.mark.parametrize(
    "fmt,args,expected",
    [
        ("%d", (1.5,), "1.500000e+00"),
        ("%.20d", (0.108,), "1.07999999999999998779e-01"),
        ("%.20d", (32,), "00000000000000000032"),
        ("%.20d", (-3,), "-00000000000000000003"),
        ("%d", (2.0**63,), "9.223372e+18"),
        ("%d", (-(2.0**63),), "-9223372036854775808"),
        ("%d", (np.uint64(2**64 - 1),), "1.844674e+19"),
        ("%u", (np.uint64(2**64 - 1),), "18446744073709551615"),
        ("%x", (-3,), "-3.000000e+00"),
        ("%x", (255,), "ff"),
        ("%#o", (8,), "010"),
        ("%d", (-0.0,), "0"),
        ("%g", (-0.0,), "-0"),
        ("%05d", (math.nan,), "  NaN"),
        ("%+d", (math.inf,), "+Inf"),
        ("%+d", (math.nan,), "NaN"),
        ("%-8d|", (32,), "32      |"),
        ("%s", (65,), "A"),
        ("%s", (70000,), "7.000000e+04"),
        ("%s", (1.5,), "1.500000e+00"),
        ("%c%c", (72, 105.5), "H1.055000e+02"),
        ("%d,", ([1, 2, 3],), "1,2,3,"),
        ("%d %s;", (1, "ab", 2), "1 ab;2 "),
        ("%s,", ("ab", "c"), "ab,c,"),
        ("%d,", ("ab",), "97,98,"),
        ("%d %d", (np.array([[1, 2], [3, 4]]),), "1 32 4"),
        ("abc%ddef", (), "abc"),
        ("a%db", ([],), "ab"),
        ("%d,", ([],), ","),
        ("%s\\n", ("",), "\n"),
        ("a\\nb\\\\c\\tq", (), "a\nb\\c\tq"),
        ('a\\"b', (), 'a"b'),
        ("a\\qb", (), "a"),
        ("a%qb", (1,), "a"),
        ("\\x41\\101", (), "AA"),
        ("100%%", (), "100%"),
        ("%*d|", (5, 3), "    3|"),
        ("%.f", (2.5,), "2"),
        ("%.0f", (1.5,), "2"),
        ("%.20f", (0.108,), "0.10799999999999999878"),
        ("%.10f", (-0.0,), "-0.0000000000"),
        ("%f", (np.int64(9007199254740993),), "9007199254740992.000000"),
        ("%d", (np.int64(9007199254740993),), "9007199254740993"),
        ("%.20g", (np.float32(0.1),), "0.10000000149011611938"),
        ("%d", (np.float32(3e9),), "3000000000"),
    ],
)
def test_sprintf(fmt, args, expected):
    assert text.sprintf(fmt, *args) == expected


# MATLAB R2024b's answers (probes of task X01, 2026-11-13): each empty argument holds one
# conversion's place and prints the field of empty text; text pads with zeros under the 0 flag
E = np.zeros
EMPTY_ARGUMENT_CASES = [
    ("a[%d]b%d", ([], 5), "a[]b5"),
    ("a[%d]b%d|", ([], 5, 6), "a[]b5|a[6]b"),
    ("%d,", (1, [], 2), "1,,2,"),
    ("%d,", (1, []), "1,,"),
    ("%d,", ([], [], 3), ",,3,"),
    ("%d,", ([], 3, []), ",3,,"),
    ("%d,", ([], 3, [], []), ",3,,,"),
    ("[%5d]", ([], 1), "[     ][    1]"),
    ("[%-5d]", ([], 1), "[     ][1    ]"),
    ("[%05.2f]", ([], 1), "[00000][01.00]"),
    ("[%+d]", ([], 1), "[][+1]"),
    ("[%s][%c]", ([], "x"), "[][x]"),
    ("[%5s]", ([], "x"), "[     ][    x]"),
    ("[%d][%s]", ("", 5), "[][\x05]"),
    ("[%c]%d", ("", 5), "[]5"),
    ("[%5c]%d", ("", 5), "[     ]5"),
    ("[%5s]", ("", "x"), "[     ][    x]"),
    ("[%.2s]", ("", "xyz"), "[][xy]"),
    ("a%db%dc", ([],), "ab"),
    ("a%db%dc", ([], []), "abc"),
    ("a%db%dc", ([], [], []), "abcab"),
    ("a%db%dc%dd", ([], 5), "ab5c"),
    ("a%db%dc%dd", (5, []), "a5bc"),
    ("%d %d %d", ([],), " "),
    ("%d %s|", ([], "x", 3, ""), " x|3 |"),
    ("%*d|", ([], 3, 4), ""),
    ("%*d|", (3, [], 4), "   |"),
    ("%.*f|", ([], 2.5), ""),
    ("%d|%d", (E(0, bool), True), "|1"),
    ("%d|", (E((0, 3)), 7), "|7|"),
    ("%s", ([], "ab"), "ab"),
    ("%d", ([], "ab"), "9798"),
    ("x%sy", ("", ""), "xyxy"),
    ("x%sy", ("", "a"), "xyxay"),
    ("x%sy", ("a", ""), "xayxy"),
    ("x%dy", ([], []), "xyxy"),
    ("x%dy", (1, []), "x1yxy"),
    ("x%dy", ([], 1), "xyx1y"),
    ("%s\\n", ("", ""), "\n\n"),
    ("<%d>", ([], [], 1), "<><><1>"),
    ("a%%%d", ([], 2), "a%a%2"),
    ("%d,%s;", (1, ""), "1,;"),
    ("%d,%s;", (1, "", 2, "x"), "1,;2,x;"),
    ("%s=%d;", ("a", [], "b", 2), "a=;b=2;"),
    ("%c%c|", ("a", "", "b"), "a|b"),
    ("%5.1f|%s", ([], "x"), "     |x"),
    ("%g %g", ([], [1, 2]), " 12 "),
    ("%d%d%d", ([], 1, [], 2), "12"),
    ("%d\\n", ([], 1), "\n1\n"),
    ("abc", ([],), "abc"),
    ("abc", ([], 1), "abc"),
    ("abc", ("",), "abc"),
    ("%d", (E(0, np.float32),), ""),
    ("%d", (E(0, np.int8),), ""),
    ("%x|", ([], 255), "|ff|"),
    ("%e|", ([], 1), "|1.000000e+00|"),
    ("%d %d", ([1, 2], [], 3), "1 2 3"),
    ("%d %d\\n", ([1, 2], [], [3, 4]), "1 2\n 3\n4 "),
    ("%d", (E((1, 0)), E((0, 1))), ""),
    ("%s-%s", ("", E((1, 0))), "-"),
    ("a%sb", ([],), "ab"),
    ("a%sb", ("",), "ab"),
    ("a\\qb%d", ([], 1), "a"),
    ("%d %y %d", ([], 1), " "),
    ("[%05s]", ("",), "[00000]"),
    ("[%05s]", ([], "ab"), "[00000][000ab]"),
    ("[%05c]", ("",), "[00000]"),
    ("[%05c]", ("", "a"), "[00000][0000a]"),
    ("[%05d]", ([],), "[00000]"),
    ("[%05d]", ([], 1), "[00000][00001]"),
    ("[%-05d]", ([], 1), "[00000][1    ]"),
    ("[% 5d]", ([], 1), "[     ][    1]"),
    ("[% d]", ([], 1), "[][ 1]"),
    ("[%#5x]", ([], 1), "[     ][  0x1]"),
    ("[%#o]", ([], 8), "[][010]"),
    ("[%5.3d]", ([], 1), "[     ][  001]"),
    ("[%5%]", ([], 1), "["),
    ("ab%*d|", ([], 3), "ab"),
    ("%d|%*d|", (1, [], 3), "1|"),
    ("%d|%.*f|", (1, [], 2.5), "1|"),
    ("%d|%*.*f|", (1, 4, [], 2.5), "1|"),
    ("%d|%*d|", (1, 3, []), "1|   |"),
    ("[%*s]", (4, ""), "[    ]"),
    ("[%-*d]", (4, []), "[    ]"),
    ("[%0*d]", (-5, []), "[00000]"),
    ("[%s]", ([],), "[]"),
    ("[%d]", ("",), "[]"),
    ("%s|%s", ("ab", [], "c"), "ab|c|"),
    ("%d", ([],), ""),
    ("", ([],), ""),
    ("", ([], 1), ""),
]
TEXT_PADDING_CASES = [
    ("[%05s]", ("ab",), "[000ab]"),
    ("[%-05s]", ("ab",), "[ab000]"),
    ("[%0-5s]", ("ab",), "[ab000]"),
    ("[%05s]", ("-ab",), "[00-ab]"),
    ("[%05s]", (" ab",), "[00 ab]"),
    ("[%05s]", (65,), "[0000A]"),
    ("[%05c]", (66,), "[0000B]"),
    ("[%-05c]", ("a",), "[a0000]"),
    ("[%05s]", ("abcdef",), "[abcdef]"),
    ("[%05.1s]", ("abc",), "[0000a]"),
    ("[%015s]", (3.5,), "[0003.500000e+00]"),
    ("[%+05s]", ("ab",), "[000ab]"),
    ("[% 05s]", ("ab",), "[000ab]"),
    ("[%#05s]", ("ab",), "[000ab]"),
    ("[%*s]", (-5, "ab"), "[ab   ]"),
    ("[%0*s]", (-5, "ab"), "[ab000]"),
    ("[%05s|%-5s]", ("a", "b"), "[0000a|b    ]"),
]


@pytest.mark.parametrize("fmt,args,expected", EMPTY_ARGUMENT_CASES + TEXT_PADDING_CASES)
def test_sprintf_empty_arguments_and_text_padding(fmt, args, expected):
    assert text.sprintf(fmt, *args) == expected


# MATLAB R2024b's answers (a probe of 78 calls): %s of a numeric argument takes the run of
# character codes (integer values in [0, 65535], logical too) that starts at the current
# element of that argument as one text field; an element that is no code goes alone, as %e
NUMERIC_TEXT_CASES = [
    ("[%s]", ([65, 66, 67],), "[ABC]"),
    ("[%s]", (np.array([[65], [66], [67]]),), "[ABC]"),
    ("[%s]", (np.array([[65, 67], [66, 68]]),), "[ABCD]"),
    ("[%s]|", ([65, 66], [67, 68]), "[AB]|[CD]|"),
    ("[%d][%s]", ([1, 65, 66],), "[1][AB]"),
    ("[%s][%d]", ([65, 66], 5), "[AB][5]"),
    ("[%s][%d]", ([65, 66, 5],), "[AB\x05]["),
    ("[%s][%d]", ([65, 66, 5], 7), "[AB\x05][7]"),
    ("[%5s]", ([65, 66],), "[   AB]"),
    ("[%-5s]", ([65, 66],), "[AB   ]"),
    ("[%05s]", ([65, 66],), "[000AB]"),
    ("[%-05s]", ([65, 66],), "[AB000]"),
    ("[%.1s]", ([65, 66, 67],), "[A]"),
    ("[%5.2s]", ([65, 66, 67],), "[   AB]"),
    ("[%.0s]", ([65, 66, 67],), "[]"),
    ("[%s]", (np.array([65, 66], np.uint8),), "[AB]"),
    ("[%s]", (np.array([65, 66], np.int16),), "[AB]"),
    ("[%s]", (np.array([65.0, 66.0], np.float32),), "[AB]"),
    ("[%s]", (np.array([True, False, True]),), "[\x01\x00\x01]"),
    ("[%s]", (np.array([True, True]),), "[\x01\x01]"),
    ("[%s]", (np.array([-1, 65], np.int8),), "[-1.000000e+00][A]"),
    ("[%s]", (np.array([65, 65535], np.uint16),), "[A\uffff]"),
    ("[%s]", ([0, 65],), "[\x00A]"),
    ("[%s]", ([65, 66], "cd"), "[AB][cd]"),
    ("[%s]", ("ab", [67, 68]), "[ab][CD]"),
    ("%s=%d;", ([65, 66], 3), "AB=3;"),
    ("[%c]", ([65, 66],), "[A][B]"),
    ("[%s|%c]", ([65, 66, 67],), "[ABC|"),
    ("[%c|%s]", ([65, 66, 67],), "[A|BC]"),
    ("%s,", ([65, 66], 67), "AB,C,"),
    ("%s", ("a", 66), "aB"),
    ("[%*s]", ([3, 65, 66],), "[ AB]"),
    ("[%*s]", (3, [65, 66]), "[ AB]"),
    ("[%-*s]", (4, [65, 66]), "[AB  ]"),
    ("[%s]", (65,), "[A]"),
    ("[%s]", (True,), "[\x01]"),
    ("[%s]", ([65, 66.5],), "[A][6.650000e+01]"),
    ("[%s]", ([66.5, 65],), "[6.650000e+01][A]"),
    ("[%s]", ([65, 66.5, 67],), "[A][6.650000e+01][C]"),
    ("[%s]", ([1.5, 2.5],), "[1.500000e+00][2.500000e+00]"),
    ("[%s]", ([65, -1],), "[A][-1.000000e+00]"),
    ("[%s]", ([-1, 65],), "[-1.000000e+00][A]"),
    ("[%s]", ([65, 70000],), "[A][7.000000e+04]"),
    ("[%s]", ([65, 65535],), "[A\uffff]"),
    ("[%s]", ([65, 65536],), "[A][6.553600e+04]"),
    ("[%s]", ([65, math.nan],), "[A][NaN]"),
    ("[%s]", ([math.nan, 65],), "[NaN][A]"),
    ("[%s]", ([65, math.inf],), "[A][Inf]"),
    ("[%s]", ([-math.inf, 65],), "[-Inf][A]"),
    ("[%s]", ([math.nan, math.nan],), "[NaN][NaN]"),
    ("[%10s]", ([65, 66.5],), "[         A][6.650000e+01]"),
    ("[%-14s|]", ([1.5, 2],), "[1.500000e+00  |][\x02             |]"),
    ("[%.2s]", ([65, 66.5],), "[A][6.65e+01]"),
    ("[%05s]", ([65, 66.5],), "[0000A][6.650000e+01]"),
    ("[%s|%s]", ([65, 66.5], 67), "[A|6.650000e+01][C|"),
    ("[%d|%s]", ([1, 65, 66.5],), "[1|A][6.650000e+01|"),
    ("[%s]", (np.array([65, 1099511627776], np.int64),), "[A][1.099512e+12]"),
    ("[%s]", (np.uint64(18446744073709551615),), "[1.844674e+19]"),
    ("[%s]", (-0.0,), "[\x00]"),
    ("[%s]", ([-0.0, 65],), "[\x00A]"),
    ("[%s]", (np.array([65.0, 66.5], np.float32),), "[A][6.650000e+01]"),
    ("[%s]", (np.array([[65.5, 66], [67, 68]]),), "[6.550000e+01][CBD]"),
    ("[%s]", ([65, 66.25, 10000000000],), "[A][6.625000e+01][1.000000e+10]"),
    ("[%+s]", ([1.5, 65],), "[+1.500000e+00][A]"),
    ("[%s]", ([1e-05, 65],), "[1.000000e-05][A]"),
    ("[%.3s]", ([65, 66, 67, 68, 69],), "[ABC]"),
    ("[%.3s]", ([65, 66.5, 67],), "[A][6.650e+01][C]"),
    ("[%s][%s]", ([65, 66], []), "[AB][]"),
    ("[%s][%s]", ([], [65, 66]), "[][AB]"),
    ("[%s]\\n", (np.array([[65, 66], [67, 68]]),), "[ACBD]\n"),
    ("%s %s", ([72, 105], [33, 33]), "Hi !!"),
    ("%s", ([104, 105], " ", 33), "hi !"),
    ("[%s]", ([65, 66.0000000001],), "[A][6.600000e+01]"),
    ("[%s]", (np.float32(16777216.0),), "[1.677722e+07]"),
    ("[%s]", ([65, 66],), "[AB]"),
    ("[%d,%s]", ([1, 65.5, 2, 66],), "[1,6.550000e+01][2,B]"),
    ("[%s,%d]", ([65.5, 2, 66],), "[6.550000e+01,2][B,"),
    ("[%s,%d]", ([65, 66, 2.5],), "[AB,2.500000e+00]"),
]


@pytest.mark.parametrize("fmt,args,expected", NUMERIC_TEXT_CASES)
def test_sprintf_text_of_numeric_arguments(fmt, args, expected):
    assert text.sprintf(fmt, *args) == expected


def test_sprintf_text_runs_longer_than_the_look_ahead():
    """the run of character codes is looked for in growing chunks; runs end anywhere in them
    (the first two are oracle cases too)"""
    assert (
        text.sprintf("[%s]", np.array([65.0] * 100 + [66.5])) == "[" + "A" * 100 + "][6.650000e+01]"
    )
    w = [65] * 17 + [1.5] + [66] * 40
    assert text.sprintf("%s|", w) == "A" * 17 + "|1.500000e+00|" + "B" * 40 + "|"
    assert text.sprintf("[%s]", np.full(1000, 66, np.uint16)) == "[" + "B" * 1000 + "]"
    assert text.sprintf("%s,", np.arange(40) + 0.5) == "".join(f"{v + 0.5:e}," for v in range(40))


def test_sprintf_rejects_cells():
    with pytest.raises(TypeError, match="invalidInputType"):
        text.sprintf("%s", ["a", "b"])


@pytest.mark.parametrize(
    "x,expected",
    [
        (0.108, "0.108"),
        (math.pi, "3.1416"),
        ([1.5, -2.25, 3], "1.5       -2.25           3"),
        (123456.789, "123456.789"),
        (1e-7, "1e-07"),
        (32, "32"),
        ([1, 2, 10], "1   2  10"),
        (1e20, "1e+20"),
        (-0.0, "0"),
        (math.nan, "NaN"),
        (True, "1"),
        (np.int16(-5), "-5"),
        ([], ""),
    ],
)
def test_num2str(x, expected):
    assert text.num2str(x) == expected


def test_num2str_with_format_and_precision():
    assert text.num2str([1.5, 22.25], "%.10f,") == "1.5000000000,22.2500000000,"
    assert text.num2str(np.pi, 8) == "3.1415927"
    assert text.num2str("text") == "text"
    with pytest.raises(ValueError):
        text.num2str(1, "no percent")


def test_int2str():
    assert text.int2str(2.5) == "3" and text.int2str(-2.5) == "-3" and text.int2str(-0.4) == "0"
    assert text.int2str([1, 22, -3]) == "1  22  -3"
    with pytest.raises(ValueError, match="subsassigndimmismatch"):
        text.int2str([0.001, 0.002])  # MATLAB's own int2str errors here


def test_mat2str_and_comma():
    assert text.mat2str(0.1) == "0.1"
    assert text.mat2str(1 / 3) == "0.333333333333333"
    assert text.mat2str(np.float32(np.pi)) == "3.14159274101257"
    assert text.mat2str([1, math.nan, -math.inf]) == "[1 NaN -Inf]"
    assert text.mat2str(np.array([[1, 2], [3, 4]])) == "[1 2;3 4]"
    assert text.mat2str(np.zeros((1, 0))) == "zeros(1,0)" and text.mat2str(np.zeros((0, 0))) == "[]"
    assert text.mat2str(np.array([True, False])) == "[true false]"
    assert text.mat2str("it's") == "'it''s'"
    assert text.mat2str(np.int64(2**62)) == "4611686018427387904"
    assert text.mat2str(np.int64(2**62), 6) == "4.61169e+18"
    assert text.mat2str(np.array([1.0, 2.0], np.float32), cls=True) == "single([1 2])"
    assert text.mat2str_comma([1, 2, 3]) == "[1,2,3]"
    assert text.mat2str_comma(5) == "5"
    assert text.mat2str_comma([0.5, 1.25], 3) == "[0.500,1.250]"


def test_str2double():
    assert text.str2double("1,000") == 1000 and text.str2double(" 2.5 ") == 2.5
    assert text.str2double("0x10") == 16 and text.str2double("-Inf") == -math.inf
    assert text.str2double("1e400") == math.inf
    for s in ("abc", "", "Infinity", "1 2", "+-1", "NaN"):
        v = text.str2double(s)
        assert math.isnan(v)
        # MATLAB's NaN, with the sign bit set
        assert (
            np.float64(v).view(np.uint64)
            == np.float64(dtypes.NAN).view(np.uint64)
            == 0xFFF8000000000000
        )
    assert math.isnan(text.str2double(5))


def test_through_text():
    # a pixel size through '%.10f' and '%.20d' (MATLAB's wrappers)
    assert text.through_text(0.1234567890123, "%.10f") == 0.123456789
    assert text.through_text(0.108, "%.20d") == 0.108
    np.testing.assert_array_equal(
        text.through_text(np.array([1 / 3, 2.0]), "%.10f"), [0.3333333333, 2.0]
    )
    assert text.through_mat2str(1 / 3) == 0.333333333333333
    assert math.isinf(text.through_mat2str(math.inf))
