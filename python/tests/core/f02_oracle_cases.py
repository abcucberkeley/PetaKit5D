# ruff: noqa: E501, E731
"""Oracle cases of petakit5d.core: the MATLAB snippet, the Python call and their inputs.

Stand-in for F04's case objects (petakit5d.testing.cases is not in this tree yet): each case
names a MATLAB snippet that sets ``out`` from its inputs (MATLAB helpers in ``matlab/``), the
Python callable with the same inputs, and how the output is compared (``array``: class, size and
bytes; ``list``: a list of strings or numbers; ``text``: one string). References are the
SHA-256 values in ``oracle_f02.json`` (from MATLAB R2024b at the oracle commit). The pack is
run on the cluster by the task's runner; see docs/conventions.md, "How the rules were measured".
"""

from __future__ import annotations

import dataclasses
import math
import os
from collections.abc import Callable
from fractions import Fraction
from typing import Any

import _f02_gen as _gen  # this folder is on sys.path (the test module and the pack runner add it)
import numpy as np

from petakit5d.core import arrays, axes, blocks, dtypes, jsonenc, paths, text, textfiles
from petakit5d.core import bbox as bbox_

SEED = 20261111


@dataclasses.dataclass
class OCase:
    id: str
    matlab: str
    python: Callable[..., Any]
    inputs: dict[str, Callable[[], Any]] = dataclasses.field(default_factory=dict)
    kind: str = "array"  # array | list | text
    setup: str = ""  # 'tree': inputs include a file tree; the callable gets its root
    # roles (petakit5d.testing.cases): the inputs are MATLAB's; roles[name] ('vol', 'box', ...)
    # gives the Python callable an input in Python's form, out_role maps its result back
    roles: dict[str, str] = dataclasses.field(default_factory=dict)
    out_role: Any = None

    def run_python(self, inputs: dict[str, Any]) -> Any:
        """the Python side: inputs through their roles, the result through out_role"""
        from petakit5d.testing.cases import Role, out_to_matlab, role_to_python

        args = {
            k: role_to_python(Role(self.roles[k], v), v) if k in self.roles else v
            for k, v in inputs.items()
        }
        return out_to_matlab(self.out_role, self.python(**args))


CASES: list[OCase] = []


def case(id, matlab, python, kind="array", setup="", roles=None, out_role=None, **inputs):
    CASES.append(OCase(id, matlab, python, inputs, kind, setup, dict(roles or {}), out_role))


# -- conversions -------------------------------------------------------------------------------

FLOAT_SRC = {"double": np.float64, "single": np.float32}
TARGETS = (*_gen.INT_CLASSES, "single", "double")

for _src, _dt in FLOAT_SRC.items():
    _vals = lambda dt=_dt: _gen.float_values(SEED, 200000, dt)
    for _cls in TARGETS:
        case(
            f"conv.{_src}.{_cls}.direct",
            f"out = pkpy_f02_conv(a, '{_cls}');",
            lambda a, c=_cls: dtypes.to_class(a, c),
            a=_vals,
        )
        case(
            f"conv.{_src}.{_cls}.cast",
            f"out = cast(a, '{_cls}');",
            lambda a, c=_cls: dtypes.cast(a, c),
            a=_vals,
        )
    case(
        f"conv.{_src}.logical",
        "out = pkpy_f02_conv(a(~isnan(a)), 'logical');",
        lambda a: dtypes.to_class(a[~np.isnan(a)], "logical"),
        a=_vals,
    )
    case(f"round.{_src}", "out = round(a);", lambda a: dtypes.mround(a), a=_vals)
    for _cls in ("int8", "uint8", "int16", "uint16", "int32"):
        # a scalar as the argument of a function file takes the array rounding ...
        case(
            f"conv.{_src}.{_cls}.scalars",
            f"out = arrayfun(@(v) pkpy_f02_conv(v, '{_cls}'), a);",
            lambda a, c=_cls: np.array(
                [dtypes.to_class(v, c, scalar_fast=True) for v in a], dtype=c
            ),
            a=(lambda dt=_dt: _gen.float_values(SEED + 7, 4000, dt)),
        )
        # ... and so does one in an anonymous function made in eval'd code (a local scalar in
        # a function or script rounds exactly; see docs/conventions.md, 3.1)
        case(
            f"conv.{_src}.{_cls}.scalars_eval_anon",
            f"out = arrayfun(@(v) {_cls}(v), a);",
            lambda a, c=_cls: np.array(
                [dtypes.to_class(v, c, scalar_fast=True) for v in a], dtype=c
            ),
            a=(lambda dt=_dt: _gen.float_values(SEED + 7, 4000, dt)),
        )

for _src in _gen.INT_CLASSES:
    _vals = lambda s=_src: _gen.int_values(SEED + 1, 20000, s)
    for _cls in (*TARGETS, "logical"):
        case(
            f"conv.{_src}.{_cls}.direct",
            f"out = pkpy_f02_conv(a, '{_cls}');",
            lambda a, c=_cls: dtypes.to_class(a, c),
            a=_vals,
        )

# -- integer arithmetic ------------------------------------------------------------------------

OPS = ("plus", "minus", "times", "rdivide")
_PYOPS = {
    "plus": dtypes.plus,
    "minus": dtypes.minus,
    "times": dtypes.times,
    "rdivide": dtypes.rdivide,
    "ldivide": dtypes.ldivide,
}


def _table(a, op):
    A, B = np.meshgrid(a, a, indexing="ij")
    return _PYOPS[op](A, B)


def _scalar_table(a, d, op, side):
    out = np.zeros((a.size, d.size), dtype=a.dtype)
    for k in range(d.size):
        dk = float(d[k])
        out[:, k] = _PYOPS[op](a, dk) if side == "r" else _PYOPS[op](dk, a)
    return out


for _cls in _gen.INT_CLASSES:
    _edges = lambda c=_cls: _gen.int_edges(c)
    for _op in OPS:
        case(
            f"arith.{_cls}.{_op}.table",
            f"out = pkpy_f02_table(a, '{_op}');",
            lambda a, op=_op: _table(a, op),
            a=_edges,
        )
        for _side in "rl":
            case(
                f"arith.{_cls}.{_op}.double_{_side}",
                f"out = pkpy_f02_scalar_table(a, d, '{_op}', '{_side}');",
                lambda a, d, op=_op, s=_side: _scalar_table(a, d, op, s),
                a=_edges,
                d=lambda: _gen.DOUBLE_SCALARS,
            )
    for _opt in ("fix", "round", "floor", "ceil"):
        case(
            f"idivide.{_cls}.{_opt}",
            f"out = pkpy_f02_idivide_table(a, '{_opt}');",
            lambda a, o=_opt: dtypes.idivide(*np.meshgrid(a, a, indexing="ij"), o),
            a=_edges,
        )
    case(f"uminus.{_cls}", "out = -a;", lambda a: dtypes.uminus(a), a=_edges)


NEAR_HALVES = np.array(
    sorted(
        {
            float(x)
            for k in range(-3, 4)
            for x in (k + 0.5, np.nextafter(k + 0.5, -np.inf), np.nextafter(k + 0.5, np.inf))
        }
        | {
            1 / 3,
            -1 / 3,
            0.1,
            -0.1,
            2.0**-60,
            -(2.0**-60),
            1.0 - 2.0**-53,
            2.0**52 + 0.5,
            2.0**53,
            -(2.0**53),
        }
    )
)
for _cls in ("int64", "uint64", "int32", "uint16", "int8"):
    _small = lambda c=_cls: np.unique(
        np.concatenate(
            [_gen.int_edges(c), np.array([0, 1, 2, 3, 4, 5, 6, 7, 9, 10], dtype=c)]
            + ([np.array([-1, -2, -3, -4, -9, -10], dtype=c)] if c[0] == "i" else [])
        )
    )
    for _op in OPS:
        for _side in "rl":
            case(
                f"arith.{_cls}.{_op}.halves_{_side}",
                f"out = pkpy_f02_scalar_table(a, d, '{_op}', '{_side}');",
                lambda a, d, op=_op, s=_side: _scalar_table(a, d, op, s),
                a=_small,
                d=lambda: NEAR_HALVES,
            )


def _near_half_pairs(kind: str):
    """uint16 n and double d where the double result of n*d (or n/d) is k + 1/2 and the exact
    result is not: these tell whether MATLAB rounds the double result or the exact one."""
    u = _gen.uniform(SEED + 3, 40000)
    ns, ds = [], []
    for i in range(0, u.size, 2):
        n = 1 + int(u[i] * 60000)
        k = int(u[i + 1] * n * 3)
        d = (k + 0.5) / n if kind == "times" else n / (k + 0.5)
        if d == 0:
            continue
        exact = Fraction(n) * Fraction(d) if kind == "times" else Fraction(n) / Fraction(d)
        fl = n * d if kind == "times" else n / d
        if fl == k + 0.5 and exact != Fraction(2 * k + 1, 2):
            ns.append(n)
            ds.append(d)
        if len(ns) >= 300:
            break
    return np.array(ns, dtype=np.uint16), np.array(ds, dtype=np.float64)


for _op in ("times", "rdivide"):
    case(
        f"arith.uint16.{_op}.near_half",
        f"out = zeros(size(a), 'uint16'); for i = 1:numel(a); out(i) = pkpy_f02_arith(a(i), d(i), '{_op}'); end",
        lambda a, d, op=_op: np.array(
            [_PYOPS[op](a[i], float(d[i])) for i in range(a.size)], dtype=np.uint16
        ),
        a=lambda op=_op: _near_half_pairs(op)[0],
        d=lambda op=_op: _near_half_pairs(op)[1],
    )

# -- single with double ------------------------------------------------------------------------

_sv = lambda: _gen.float_values(SEED + 4, 20000, np.float32)
_dv = lambda: (
    _gen.uniform(SEED + 5, _gen.float_values(SEED + 4, 20000, np.float32).size) * 3.0 + 0.01
)
for _op in ("plus", "minus", "times", "rdivide", "ldivide"):
    for _side in "rl":
        case(
            f"mixed.single.{_op}.double_scalar_{_side}",
            f"out = pkpy_f02_scalar_table(a, d, '{_op}', '{_side}');",
            lambda a, d, op=_op, s=_side: _scalar_table(a, d, op, s),
            a=_sv,
            d=lambda: _gen.DOUBLE_SCALARS,
        )
    case(
        f"mixed.single.{_op}.double_array",
        f"out = pkpy_f02_arith(a, b, '{_op}');",
        lambda a, b, op=_op: _PYOPS[op](a, b),
        a=_sv,
        b=_dv,
    )
    case(
        f"mixed.double_array.{_op}.single",
        f"out = pkpy_f02_arith(b, a, '{_op}');",
        lambda a, b, op=_op: _PYOPS[op](b, a),
        a=_sv,
        b=_dv,
    )
    case(
        f"mixed.single.{_op}.single",
        f"out = pkpy_f02_arith(a, single(b), '{_op}');",
        lambda a, b, op=_op: _PYOPS[op](a, b.astype(np.float32)),
        a=_sv,
        b=_dv,
    )
case(
    "mixed.single.plus.logical",
    "out = pkpy_f02_arith(a, true, 'plus');",
    lambda a: dtypes.plus(a, True),
    a=_sv,
)
case(
    "mixed.logical.plus.logical",
    "out = pkpy_f02_arith(true(3, 1), false(3, 1), 'plus');",
    lambda: dtypes.plus(np.ones(3, bool), np.zeros(3, bool)),
)
case(
    "mixed.uint16.plus.char",
    "out = pkpy_f02_arith(a, 'a', 'plus');",
    lambda a: dtypes.plus(a, "a"),
    a=lambda: _gen.int_edges("uint16"),
)

# -- formatting --------------------------------------------------------------------------------


def fmt_values() -> np.ndarray:
    v = [
        0.0,
        -0.0,
        1.0,
        -1.0,
        0.5,
        1.5,
        2.5,
        -2.5,
        0.1,
        0.108,
        32.45,
        1 / 3,
        math.pi,
        1e-5,
        1e-20,
        1e15,
        1e16,
        1e20,
        2.0**53,
        2.0**53 + 2,
        2.0**62,
        2.0**63,
        -(2.0**63),
        2.0**63 - 1024,
        -(2.0**63) - 2048,
        2.0**64,
        2.0**64 - 2048,
        1e300,
        5e-324,
        2.2250738585072014e-308,
        123456789.123,
        0.005,
        0.9,
        65535.0,
        65535.5,
        70000.0,
        math.nan,
        math.inf,
        -math.inf,
        0.1 + 0.2,
        0.1 + 0.7,
        999999.5,
        1e-4,
        0.0001234,
        100.0,
        99999.0,
        123456.0,
        1234567.0,
        0.3,
        0.45,
        1e5,
        1e6,
        1e-6,
        0.00012,
        5e-5,
        1.2e-5,
        31.999999999999996,
        2.675,
        1.005,
        0.125,
        0.0625,
        7.0,
        10.0,
        1e10,
        12345678901.5,
    ]
    for e in range(-12, 26):
        p = _gen.pow10(e)
        v += [p, 1.5 * p, -2.25 * p, p / 3, p - 1, 123456789 * _gen.pow10(e - 8)]
    u = _gen.uniform(SEED + 6, 900)
    v += list((u[:300] * 2 - 1) * _gen.pow10_array(np.floor(u[300:600] * 40 - 15)))
    v += list(np.floor(u[600:750] * 2e6) - 1e6)
    v += list(np.floor(u[750:900] * 2000) / 8 - 125)
    return np.array(v, dtype=np.float64)


SPRINTF_FORMATS = (
    "%d",
    "%.20d",
    "%5d",
    "%-8d|",
    "%+d",
    "%05d",
    "%x",
    "%i",
    "%u",
    "%.20f",
    "%.10f",
    "%g",
    "%0.5g",
    "%.f",
    "%f",
    "%0.6f",
    "%e",
    "%s",
    "%04d",
    "%010d",
    "%5.1f",
    "%.3d",
    "%.15g",
    "%G",
    "%E",
    "%c",
    "%o",
)


def _matlab_str(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def _to_single(v):
    with np.errstate(over="ignore"):
        return v.astype(np.float32)


_SRC_VALUES = {
    "double": fmt_values,
    "single": lambda: _to_single(fmt_values()),
    "int16": lambda: _gen.int_edges("int16"),
    "uint8": lambda: _gen.int_edges("uint8"),
    "int64": lambda: _gen.int_edges("int64"),
    "uint64": lambda: _gen.int_edges("uint64"),
    "logical": lambda: np.array([True, False]),
}
for _src, _mk in _SRC_VALUES.items():
    for _i, _f in enumerate(SPRINTF_FORMATS):
        case(
            f"sprintf.{_src}.{_i:02d}",
            f"out = arrayfun(@(v) sprintf({_matlab_str(_f)}, v), a, 'UniformOutput', false);",
            lambda a, f=_f: [text.sprintf(f, v) for v in a],
            kind="list",
            a=_mk,
        )
    case(
        f"num2str.{_src}",
        "out = arrayfun(@(v) num2str(v), a, 'UniformOutput', false);",
        lambda a: [text.num2str(v) for v in a],
        kind="list",
        a=_mk,
    )
    case(
        f"num2str.{_src}.fmt10f",
        "out = arrayfun(@(v) num2str(v, '%.10f,'), a, 'UniformOutput', false);",
        lambda a: [text.num2str(v, "%.10f,") for v in a],
        kind="list",
        a=_mk,
    )
    case(
        f"num2str.{_src}.prec8",
        "out = arrayfun(@(v) num2str(v, 8), a, 'UniformOutput', false);",
        lambda a: [text.num2str(v, 8) for v in a],
        kind="list",
        a=_mk,
    )
    case(
        f"int2str.{_src}",
        "out = arrayfun(@(v) int2str(v), a, 'UniformOutput', false);",
        lambda a: [text.int2str(v) for v in a],
        kind="list",
        a=_mk,
    )
    case(
        f"mat2str.{_src}",
        "out = arrayfun(@(v) mat2str(v), a, 'UniformOutput', false);",
        lambda a: [text.mat2str(v) for v in a],
        kind="list",
        a=_mk,
    )
    case(
        f"mat2str.{_src}.prec6",
        "out = arrayfun(@(v) mat2str(v, 6), a, 'UniformOutput', false);",
        lambda a: [text.mat2str(v, 6) for v in a],
        kind="list",
        a=_mk,
    )
    case(
        f"jsonencode.{_src}",
        "out = arrayfun(@(v) jsonencode(v), a, 'UniformOutput', false);",
        lambda a: [jsonenc.jsonencode(v) for v in a],
        kind="list",
        a=_mk,
    )
    case(
        f"jsonencode.{_src}.vector",
        "out = jsonencode(a);",
        lambda a: jsonenc.jsonencode(a),
        kind="text",
        a=_mk,
    )
    case(
        f"mat2str.{_src}.vector",
        "out = mat2str(a');",
        lambda a: text.mat2str(a.reshape(1, -1)),
        kind="text",
        a=_mk,
    )


def _rows(seed, n, ints: bool):
    u = _gen.uniform(seed, n * 5)
    m = (u[: n * 4].reshape(n, 4) * 2 - 1) * _gen.pow10_array(
        np.floor(u[n * 4 :].reshape(n, 1) * 12 - 4)
    )
    if ints:
        m = np.floor(m)
    m[::7, 1] = np.nan
    m[::11, 2] = np.inf
    m[::13, 0] = -0.0
    return m


for _kind in ("ints", "floats"):
    _mk = lambda k=_kind: _rows(SEED + 8, 300, k == "ints")
    case(
        f"num2str.rows.{_kind}",
        "out = arrayfun(@(i) num2str(a(i, :)), (1:size(a, 1))', 'UniformOutput', false);",
        lambda a: [text.num2str(a[i, :]) for i in range(a.shape[0])],
        kind="list",
        a=_mk,
    )
    case(
        f"num2str.rows.{_kind}.fmt10f",
        "out = arrayfun(@(i) num2str(a(i, :), '%.10f,'), (1:size(a, 1))', 'UniformOutput', false);",
        lambda a: [text.num2str(a[i, :], "%.10f,") for i in range(a.shape[0])],
        kind="list",
        a=_mk,
    )
    case(
        f"num2str.rows.{_kind}.fmt20d",
        "out = arrayfun(@(i) strrep(num2str(a(i, :), '%.20d,'), ' ', ''), (1:size(a, 1))', 'UniformOutput', false);",
        lambda a: [text.num2str(a[i, :], "%.20d,").replace(" ", "") for i in range(a.shape[0])],
        kind="list",
        a=_mk,
    )
    case(
        f"int2str.rows.{_kind}",
        "out = arrayfun(@(i) int2str(a(i, :)), (1:size(a, 1))', 'UniformOutput', false);",
        lambda a: [text.int2str(a[i, :]) for i in range(a.shape[0])],
        kind="list",
        a=_mk,
    )
    case(
        f"mat2str.rows.{_kind}",
        "out = arrayfun(@(i) mat2str(a(i, :)), (1:size(a, 1))', 'UniformOutput', false);",
        lambda a: [text.mat2str(a[i, :]) for i in range(a.shape[0])],
        kind="list",
        a=_mk,
    )
    case(
        f"mat2str_comma.rows.{_kind}",
        "out = arrayfun(@(i) mat2str_comma(a(i, :)), (1:size(a, 1))', 'UniformOutput', false);",
        lambda a: [text.mat2str_comma(a[i, :]) for i in range(a.shape[0])],
        kind="list",
        a=_mk,
    )
    case(
        f"mat2str_comma.rows.{_kind}.sn3",
        "out = arrayfun(@(i) mat2str_comma(a(i, :), 3), (1:size(a, 1))', 'UniformOutput', false);",
        lambda a: [text.mat2str_comma(a[i, :], 3) for i in range(a.shape[0])],
        kind="list",
        a=_mk,
    )
    case(
        f"mat2str.matrix.{_kind}",
        "out = mat2str(a(1:5, :));",
        lambda a: text.mat2str(a[:5, :]),
        kind="text",
        a=_mk,
    )
    case(
        f"jsonencode.matrix.{_kind}",
        "out = jsonencode(a(1:5, :));",
        lambda a: jsonenc.jsonencode(a[:5, :]),
        kind="text",
        a=_mk,
    )

# sprintf with several arguments, escapes and empty data: (MATLAB call, Python args)
_PI = math.pi
SPRINTF_CALLS = [
    ("sprintf('%d,', [1 2 3])", ("%d,", [1, 2, 3])),
    ("sprintf('%d %s;', 1, 'ab', 2)", ("%d %s;", 1, "ab", 2)),
    ("sprintf('%s', [72 105])", ("%s", [72, 105])),
    ("sprintf('%d', [1 2.5 3])", ("%d", [1, 2.5, 3])),
    ("sprintf('abc%ddef')", ("abc%ddef",)),
    ("sprintf('abc%sdef%d')", ("abc%sdef%d",)),
    ("sprintf('a%db', [])", ("a%db", [])),
    ("sprintf('%d,', [])", ("%d,", [])),
    ("sprintf('%s\\n', '')", ("%s\\n", "")),
    ("sprintf('%s,', 'ab', 'c')", ("%s,", "ab", "c")),
    ("sprintf('%d,', 'ab')", ("%d,", "ab")),
    ("sprintf('%d %d', [1 2; 3 4])", ("%d %d", np.array([[1.0, 2.0], [3.0, 4.0]]))),
    ("sprintf('%*d|', 5, 3)", ("%*d|", 5, 3)),
    ("sprintf('%-*.*f|', 8, 2, pi)", ("%-*.*f|", 8, 2, _PI)),
    ("sprintf('a\\\"b')", ('a\\"b',)),
    ("sprintf('a\\nb\\tc\\\\d')", ("a\\nb\\tc\\\\d",)),
    ("sprintf('%5%|')", ("%5%|",)),
    ("sprintf('a%')", ("a%",)),
    ("sprintf('\\x41\\101')", ("\\x41\\101",)),
    ("sprintf('%c%c', 72, 105.5)", ("%c%c", 72, 105.5)),
    ("sprintf('%.20d,%.20d,%.20d', 0.108, 0.108, 0.3)", ("%.20d,%.20d,%.20d", 0.108, 0.108, 0.3)),
    ("sprintf('%s=%d', 'abc', 7)", ("%s=%d", "abc", 7)),
    ("sprintf('%s|', 'abc', 66, 'd')", ("%s|", "abc", 66, "d")),
    ("sprintf('%d|%s|', 'ab', 'cd')", ("%d|%s|", "ab", "cd")),
    ("sprintf('%c', 'hello')", ("%c", "hello")),
    ("sprintf('%5.2s|', 'hello')", ("%5.2s|", "hello")),
    ("sprintf('%d%%|', [10 20])", ("%d%%|", [10, 20])),
    ("sprintf('[%s]', '')", ("[%s]", "")),
    ("sprintf('%d', int8([-1 2]))", ("%d", np.array([-1, 2], dtype=np.int8))),
    ("sprintf('%f|', single([0.1 1e20]))", ("%f|", np.array([0.1, 1e20], dtype=np.float32))),
    ("sprintf('%s', 70000)", ("%s", 70000)),
    ("sprintf('%s', 65535)", ("%s", 65535)),
    ("sprintf('%+5s|', 'ab')", ("%+5s|", "ab")),
    ("sprintf('%#x %#o %#g', 255, 8, 1)", ("%#x %#o %#g", 255, 8, 1)),
    ("sprintf('%d', true)", ("%d", True)),
    ("sprintf('%.2f %.2f\\n', [1 2; 3 4])", ("%.2f %.2f\\n", np.array([[1.0, 2.0], [3.0, 4.0]]))),
    ("sprintf('x%sy', 'a', 'b')", ("x%sy", "a", "b")),
    ("sprintf('%s %s', 'one')", ("%s %s", "one")),
    ("sprintf('%d %d %d', 1, 2)", ("%d %d %d", 1, 2)),
    ("sprintf('no conversions')", ("no conversions",)),
    ("sprintf('no conversions', 5)", ("no conversions", 5)),
    (
        "sprintf('%.0f %.0f %.0f %.0f', 0.5, 1.5, 2.5, -0.5)",
        ("%.0f %.0f %.0f %.0f", 0.5, 1.5, 2.5, -0.5),
    ),
    ("sprintf('%.3e', 12345.6789)", ("%.3e", 12345.6789)),
    ("sprintf('% d|% d', 5, -5)", ("% d|% d", 5, -5)),
    ("sprintf('%-6.2f|%6.2f|', -1.005, 1.005)", ("%-6.2f|%6.2f|", -1.005, 1.005)),
]
for _i, (_m, _p) in enumerate(SPRINTF_CALLS):
    case(f"sprintf.call.{_i:02d}", f"out = {_m};", lambda p=_p: text.sprintf(*p), kind="text")

# str2double
STR2DOUBLE_INPUTS = [
    "1",
    " 2.5 ",
    "1,000",
    "1e3",
    "-Inf",
    "inf",
    "NaN",
    "nan",
    "abc",
    "",
    "1.2.3",
    "0x10",
    ".5",
    "5.",
    "+3",
    "1e",
    "1 2",
    "32768\n",
    "Infinity",
    "1,2,3",
    "  -0  ",
    "1E+20",
    "-.5e-3",
    "\t7",
    "8\t",
    "1,,0",
    "3.0000000000000000001",
    "123456789012345678901234567890",
    "1e400",
    "-1e-400",
    "+-1",
    "--1",
    "1+",
]
case(
    "str2double.list",
    "out = cellfun(@str2double, a);",
    lambda a: np.array([text.str2double(s) for s in a]),
    a=lambda: list(STR2DOUBLE_INPUTS),
)

# -- jsonencode of structures, writeJsonFile, text files ---------------------------------------

_K_STR = 'q"\\/\t\u00e9%'
JSON_STRUCTS = [
    (
        "struct('a', 1, 'b', 'x', 'c', [1 2], 'd', {{'p', 'q'}}, 'e', [], 'f', struct('g', true), 'h', {{}}, 'k', ['q\"\\/', char(9), char(233), '%'])",
        lambda: {
            "a": 1,
            "b": "x",
            "c": np.array([1.0, 2.0]),
            "d": ["p", "q"],
            "e": None,
            "f": {"g": True},
            "h": [],
            "k": _K_STR,
        },
    ),
    ("struct('a', {1, 2})", lambda: [{"a": 1}, {"a": 2}]),
    (
        "struct('m', [1 2; 3 4], 'z', zeros(1, 0), 'c1', {{5}}, 'emptyS', struct([]), 'col', [1; 2; 3], 'n', NaN, 'nv', [NaN 1 Inf])",
        lambda: {
            "m": np.array([[1.0, 2.0], [3.0, 4.0]]),
            "z": np.zeros((1, 0)),
            "c1": [5],
            "emptyS": [],
            "col": np.array([[1.0], [2.0], [3.0]]),
            "n": math.nan,
            "nv": np.array([math.nan, 1.0, math.inf]),
        },
    ),
    (
        "struct('dataPaths', {{'/mnt/data/a/'}}, 'xyPixelSize', 0.108, 'dz', 0.3, 'skewAngle', 32.45, 'reverse', true, "
        "'ChannelPatterns', {{'CamA_ch0', 'CamB_ch0'}}, 'bbox', [], 'save16bit', true, 'LowerLimit', 0.4, 'FFImagePaths', {{''}}, "
        "'BlockSize', [256, 256, 256], 'jobLogDir', '../job_logs', 'cpusPerTask', 1, 'configFile', '', 'resampleFactor', [], "
        "'xStepThresh', 2.35, 'tileOutBbox', [1 1 1 100 200 300], 'single', single(0.1), 'u16', uint16(7), 'ints', int32([1 -2]), "
        "'mixed', {{1, 'a', [1 2]}}, 'emptyStr', \"\", 'nested', struct('a', struct('b', {{}})))",
        lambda: {
            "dataPaths": ["/mnt/data/a/"],
            "xyPixelSize": 0.108,
            "dz": 0.3,
            "skewAngle": 32.45,
            "reverse": True,
            "ChannelPatterns": ["CamA_ch0", "CamB_ch0"],
            "bbox": None,
            "save16bit": True,
            "LowerLimit": 0.4,
            "FFImagePaths": [""],
            "BlockSize": np.array([256.0, 256.0, 256.0]),
            "jobLogDir": "../job_logs",
            "cpusPerTask": 1,
            "configFile": "",
            "resampleFactor": None,
            "xStepThresh": 2.35,
            "tileOutBbox": np.array([1.0, 1, 1, 100, 200, 300]),
            "single": np.float32(0.1),
            "u16": np.uint16(7),
            "ints": np.array([1, -2], dtype=np.int32),
            "mixed": [1, "a", np.array([1.0, 2.0])],
            "emptyStr": "",
            "nested": {"a": {"b": []}},
        },
    ),
    (
        "struct('esc', [char(0:31), char(127), char(128), '\\', '/', '\"', ''''])",
        lambda: {"esc": "".join(chr(c) for c in range(32)) + "\x7f\x80\\/\"'"},
    ),
    ("struct()", dict),
    (
        "struct('x', reshape(1:8, 2, 2, 2))",
        lambda: {"x": np.arange(1.0, 9.0).reshape((2, 2, 2), order="F")},
    ),
    ("struct('t', true(2, 2), 'c', {{{'a'}}})", lambda: {"t": np.ones((2, 2), bool), "c": [["a"]]}),
]
for _i, (_m, _p) in enumerate(JSON_STRUCTS):
    case(
        f"jsonencode.struct.{_i}",
        f"out = jsonencode({_m});",
        lambda p=_p: jsonenc.jsonencode(p()),
        kind="text",
    )
    case(
        f"jsonencode.struct.{_i}.pretty",
        f"out = jsonencode({_m}, PrettyPrint=true);",
        lambda p=_p: jsonenc.jsonencode(p(), pretty_print=True),
        kind="text",
    )
case(
    "jsonencode.nan.noconvert",
    "out = jsonencode({NaN, Inf, -Inf, [1 NaN]}, 'ConvertInfAndNaN', false);",
    lambda: jsonenc.jsonencode(
        [math.nan, math.inf, -math.inf, np.array([1.0, math.nan])], convert_inf_and_nan=False
    ),
    kind="text",
)

_READ = "f = fopen(fn, 'r'); out = fread(f, inf, '*uint8'); fclose(f); delete(fn);"


def _bytes_of(write: Callable[[str], None]) -> np.ndarray:
    import tempfile

    fd, fn = tempfile.mkstemp()
    os.close(fd)
    try:
        write(fn)
        with open(fn, "rb") as f:
            return np.frombuffer(f.read(), dtype=np.uint8).copy()
    finally:
        os.remove(fn)


WRITE_JSON = [
    (
        "struct('a', 'C:\\data\\new', 'b', 'x%03d.tif', 'c', ['tab', char(9), 'here'], 'd', ['quote\"', char(233)], 'e', 1)",
        lambda: {
            "a": "C:\\data\\new",
            "b": "x%03d.tif",
            "c": "tab\there",
            "d": 'quote"\u00e9',
            "e": 1,
        },
    ),
    (
        "struct('a', 'plain', 'b', [1 2 3], 'c', 0.108)",
        lambda: {"a": "plain", "b": np.array([1.0, 2, 3]), "c": 0.108},
    ),
    ("struct('a', 'per%cent', 'b', 2)", lambda: {"a": "per%cent", "b": 2}),
    ("struct('u', ['caf', char(233), ' ', char(8364)])", lambda: {"u": "caf\u00e9 \u20ac"}),
    ("struct('p', '100%% done', 'q', 'a%%b')", lambda: {"p": "100%% done", "q": "a%%b"}),
    ("struct('ctl', [char(1), 'x'])", lambda: {"ctl": "\x01x"}),
    ("struct('nl', ['a', newline, 'b'], 'pct', '50 %')", lambda: {"nl": "a\nb", "pct": "50 %"}),
]
for _i, (_m, _p) in enumerate(WRITE_JSON):
    case(
        f"writeJsonFile.{_i}",
        f"fn = [tempname, '.json']; writeJsonFile({_m}, fn); {_READ}",
        lambda p=_p: _bytes_of(lambda fn: textfiles.write_json_file(p(), fn)),
    )

WRITE_TEXT = [
    (
        "{'line1', ['l2 ', char(233)], '', 'x\\ny %d'}",
        lambda: ["line1", "l2 \u00e9", "", "x\\ny %d"],
        "",
    ),
    ("'raw text no newline'", lambda: "raw text no newline", ""),
    ("{}", list, ""),
    (
        "arrayfun(@(i) sprintf('n%d', i), 1:25, 'UniformOutput', false)",
        lambda: [f"n{i}" for i in range(1, 26)],
        ", 10",
    ),
    ("{['eur ', char(8364)]}", lambda: ["eur \u20ac"], ""),
    ("['two', newline, 'lines']", lambda: "two\nlines", ""),
]
for _i, (_m, _p, _extra) in enumerate(WRITE_TEXT):
    case(
        f"writeTextFile.{_i}",
        f"fn = [tempname, '.txt']; writeTextFile({_m}, fn{_extra}); {_READ}",
        lambda p=_p: _bytes_of(lambda fn: textfiles.write_text_file(p(), fn)),
    )

READ_TEXT = [
    b"a\r\nb\n\nc",
    b"a\n",
    b"",
    bytes([99, 97, 102, 195, 169, 10, 226, 130, 172]),
    bytes([233, 10, 120]),
    b"no newline",
    b"\n",
    b"\n\n",
    b"x\ry\n",
    b"tab\there\n  spaced  \n",
]
for _i, _b in enumerate(READ_TEXT):
    case(
        f"readTextFile.{_i}",
        "fn = [tempname, '.txt']; f = fopen(fn, 'w'); fwrite(f, a); fclose(f); out = readTextFile(fn); delete(fn);",
        lambda a: _read_back(a),
        kind="list",
        a=lambda b=_b: np.frombuffer(b, dtype=np.uint8).copy(),
    )


def _read_back(a: np.ndarray):
    import tempfile

    fd, fn = tempfile.mkstemp()
    with os.fdopen(fd, "wb") as f:
        f.write(a.tobytes())
    try:
        return textfiles.read_text_file(fn)
    finally:
        os.remove(fn)


# -- paths -------------------------------------------------------------------------------------

TREE = [
    "f|B",
    "f|a",
    "f|A",
    "f|b",
    "f|_x",
    "f|1",
    "f|10",
    "f|2",
    "f|.hidden",
    "f|a.b",
    "f|a_b",
    "f|a-b",
    "f|Z",
    "f|z",
    "f|~t",
    "f| sp",
    "f|a b",
    "f|#h",
    "f|a.B",
    "f|ab",
    "f|aB",
    "f|[x]",
    "f|{y}",
    "f|a,b",
    "f|x.tar.gz",
    "d|sub",
    "d|Sub2",
    "d|emptydir",
    "d|sub/inner",
    "f|sub/inner/deep.txt",
    "f|sub/s1",
    "f|sub/S0",
    "d|.hdir",
    "f|.hdir/h1",
    "d|sub/inner/deeper",
    "f|sub/inner/deeper/x.y",
    "l|link_to_sub|sub",
    "l|link_to_file|a.b",
    "l|broken|nonexistent_target",
]
TREE_UNICODE = ["f|\u00e9", "f|e", "f|\u00c9", "f|\u4e2d", "f|z", "f|\u00df"]


def _make_tree(root: str, spec: list[str]) -> None:
    for item in spec:
        t = item.split("|")
        p = os.path.join(root, t[1])
        if t[0] == "d":
            os.mkdir(p)
        elif t[0] == "f":
            with open(p, "wb") as f:
                f.write(t[1].encode("utf-8", "replace"))
        else:
            os.symlink(t[2], p)


def _with_tree(fn, spec_key: str = "tree"):
    def run(tree, **kw):
        import tempfile

        root = os.path.realpath(tempfile.mkdtemp())
        try:
            _make_tree(root, tree)
            return fn(root, **kw)
        finally:
            import shutil

            shutil.rmtree(root, ignore_errors=True)

    return run


def _rel(s: str, root: str) -> str:
    return s.replace(root, "<root>")


def _dir_names(root, sub=""):
    return [
        f"{e.name}|{int(e.isdir)}"
        for e in paths.matlab_dir(os.path.join(root, sub) if sub else root)
    ]


_TREE_M = "root = pkpy_f02_tree(tree); "
_TREE_END = " rmdir(root, 's');"
case(
    "dir.order",
    _TREE_M
    + "d = dir(root); out = arrayfun(@(e) sprintf('%s|%d', e.name, e.isdir), d, 'UniformOutput', false);"
    + _TREE_END,
    _with_tree(lambda root: _dir_names(root)),
    kind="list",
    tree=lambda: list(TREE),
)
case(
    "dir.order.sub",
    _TREE_M
    + "d = dir([root, '/sub']); out = arrayfun(@(e) sprintf('%s|%d', e.name, e.isdir), d, 'UniformOutput', false);"
    + _TREE_END,
    _with_tree(lambda root: _dir_names(root, "sub")),
    kind="list",
    tree=lambda: list(TREE),
)
case(
    "dir.order.unicode",
    _TREE_M
    + "d = dir(root); out = arrayfun(@(e) sprintf('%s|%d', e.name, e.isdir), d, 'UniformOutput', false);"
    + _TREE_END,
    _with_tree(lambda root: _dir_names(root)),
    kind="list",
    tree=lambda: list(TREE_UNICODE),
)
case(
    "dir.pattern",
    _TREE_M
    + "d = dir([root, '/a*']); out = arrayfun(@(e) sprintf('%s|%d', e.name, e.isdir), d, 'UniformOutput', false);"
    + _TREE_END,
    _with_tree(
        lambda root: [
            f"{e.name}|{int(e.isdir)}" for e in paths.matlab_dir(os.path.join(root, "a*"))
        ]
    ),
    kind="list",
    tree=lambda: list(TREE),
)
for _mode in ("file", "folder", "other"):
    case(
        f"dir_recursive.{_mode}",
        _TREE_M + f"out = strrep(dir_recursive(root, '{_mode}'), root, '<root>');" + _TREE_END,
        _with_tree(lambda root, m=_mode: [_rel(s, root) for s in paths.dir_recursive(root, m)]),
        kind="list",
        tree=lambda: list(TREE),
    )
case(
    "dir_recursive.trailing",
    _TREE_M + "out = strrep(dir_recursive([root, '/']), root, '<root>');" + _TREE_END,
    _with_tree(lambda root: [_rel(s, root) for s in paths.dir_recursive(root + "/")]),
    kind="list",
    tree=lambda: list(TREE),
)
case(
    "dir_recursive.missing",
    "out = dir_recursive([tempname, '/nothing']);",
    lambda: paths.dir_recursive("/nonexistent_pk5d_f02/nothing"),
    kind="list",
)

SIMPLIFY = [
    "sub",
    "sub/",
    "./sub",
    "sub/../Sub2",
    "a.b",
    "sub/inner/deep.txt",
    "link_to_sub",
    "link_to_file",
    "broken",
    "nonexist",
    "*.b",
    "a*",
    "sub//inner",
    ".",
    "emptydir",
    "sub/inner/",
    "x.tar.gz",
    "link_to_sub/s1",
    "./",
    "z*",
]
case(
    "simplifyPath",
    _TREE_M
    + "old = cd(root); out = cellfun(@(p) strrep(simplifyPath(p), root, '<root>'), plist, 'UniformOutput', false); cd(old);"
    + _TREE_END,
    _with_tree(lambda root, plist: _simplify_all(root, plist)),
    kind="list",
    tree=lambda: list(TREE),
    plist=lambda: list(SIMPLIFY),
)


def _simplify_all(root, ps):
    old = os.getcwd()
    os.chdir(root)
    try:
        return [_rel(paths.simplify_path(p), root) for p in ps]
    finally:
        os.chdir(old)


FILEPARTS = [
    "",
    "a",
    "a.b",
    "a.b.c",
    ".bashrc",
    "a.",
    "/",
    "/a",
    "/a/",
    "/a/b.c",
    "a/b",
    "a/b/",
    "a//b",
    "./a.txt",
    "../x",
    "/a/b/c.tar.gz",
    "dir.d/file",
    "dir.d/",
    "C:\\x\\y.z",
    "a\\b.c",
    "/.x",
    "~/f.m",
    "x/.",
    "x/..",
    "..",
    ".",
    "/a/b/.hidden.txt",
    "no_ext",
    "trailing.dot.",
    " sp ace.t xt",
    "a/b.c/d",
]
case(
    "fileparts.list",
    "out = cellfun(@(p) pkpy_f02_fileparts(p), plist, 'UniformOutput', false);",
    lambda plist: ["|".join(paths.fileparts(p)) for p in plist],
    kind="list",
    plist=lambda: list(FILEPARTS),
)

# -- geometry ----------------------------------------------------------------------------------


# The geometry cases run the Python forms (our boxes and (z, y, x) sizes) and convert MATLAB's
# parameters in and the results back with the core helpers, so they check the conversions and
# MATLAB's batch order against the stored MATLAB references (the MATLAB calls are unchanged).


def _zarr_chunks_py(P: np.ndarray) -> np.ndarray:
    out = []
    for i, p in enumerate(P, start=1):
        bbox = bbox_.box_from_matlab(p[14:20]) if p[13] else None
        try:
            b, r, loc = blocks.zarr_chunk_coordinates(
                axes.from_yxz(p[0:3]),
                batch_size=axes.from_yxz(p[3:6]),
                block_size=axes.from_yxz(p[6:9]),
                border_size=axes.from_yxz(p[9:12]),
                same_batch_size=bool(p[12]),
                bbox=bbox,
            )
            mb, mr, ml = (bbox_.boxes_to_matlab(x) for x in (b, r, loc))
            out.append(
                np.concatenate(
                    [
                        [i, mb.shape[0]],
                        mb.ravel(order="F"),
                        mr.ravel(order="F"),
                        ml.ravel(order="F"),
                    ]
                )
            )
        except (ValueError, IndexError):
            out.append(np.array([i, -1]))
    return np.concatenate(out).astype(np.float64)


def _subvol_py(P: np.ndarray) -> np.ndarray:
    out = []
    for i, p in enumerate(P, start=1):
        try:
            boxes = blocks.sub_volume_coordinates(
                axes.from_yxz(p[0:3]),
                chunk_size=axes.from_yxz(p[3:6]),
                overlap_size=p[6],
                max_sub_volume=p[7],
                use_good_num=bool(p[8]),
            )
            m = bbox_.boxes_to_matlab(boxes)  # [ymin xmin zmin ymax xmax zmax]
            out.append(
                np.concatenate(
                    [[i, m.shape[0]], m[:, 1], m[:, 4], m[:, 0], m[:, 3], m[:, 2], m[:, 5]]
                )
            )
        except (ValueError, IndexError):
            out.append(np.array([i, -1]))
    return np.concatenate(out).astype(np.float64)


def _bbox_blocks_py(P: np.ndarray) -> np.ndarray:
    out = []
    for i, p in enumerate(P, start=1):
        try:
            s = blocks.bbox_to_blocks(
                bbox_.box_from_matlab(p[0:6]),
                np.asarray(axes.from_yxz(p[6:9])) - 1,
                axes.from_yxz(p[9:12]),
                axes.from_yxz(p[12:15]),
                axes.from_yxz(p[15:18]),
            )
            out.append(
                np.concatenate(
                    [
                        [i, s.block_sub.shape[0]],
                        (s.block_sub[:, [1, 2, 0]] + 1).ravel(order="F"),
                        s.block_ind + 1,
                        bbox_.boxes_to_matlab(s.b_coords).ravel(order="F"),
                        bbox_.boxes_to_matlab(s.w_coords).ravel(order="F"),
                        bbox_.boxes_to_matlab(s.bbox_coords).ravel(order="F"),
                        axes.to_yxz(np.asarray(s.b_sub_size)),
                    ]
                )
            )
        except (ValueError, IndexError):
            out.append(np.array([i, -1]))
    return np.concatenate(out).astype(np.float64)


case(
    "zarrChunkCoordinates",
    "out = pkpy_f02_zarr_chunks(P);",
    _zarr_chunks_py,
    P=lambda: _gen.bbox_params(SEED + 10, 800),
)
case(
    "subVolumeCoordinates",
    "out = pkpy_f02_subvol(P);",
    _subvol_py,
    P=lambda: _gen.subvol_params(SEED + 11, 500),
)
case(
    "bboxToBlocks",
    "out = pkpy_f02_bbox_blocks(P);",
    _bbox_blocks_py,
    P=lambda: _gen.blocks_params(SEED + 12, 800),
)
case(
    "findGoodFactorNumber.up",
    "out = arrayfun(@(n) findGoodFactorNumber(n, 1), a);",
    lambda a: np.array([blocks.find_good_factor_number(int(n), 1) for n in a], dtype=np.float64),
    a=lambda: np.arange(1.0, 4001.0),
)
case(
    "findGoodFactorNumber.down",
    "out = arrayfun(@(n) findGoodFactorNumber(n, -1), a);",
    lambda a: np.array([blocks.find_good_factor_number(int(n), -1) for n in a], dtype=np.float64),
    a=lambda: np.arange(1.0, 4001.0),
)
case(
    "findGoodFactorNumber.odd7",
    "out = arrayfun(@(n) findGoodFactorNumber(n, 1, true), a);",
    lambda a: np.array(
        [blocks.find_good_factor_number(int(n), 1, True) for n in a], dtype=np.float64
    ),
    a=lambda: np.arange(1.0, 4001.0),
)

AXIS_ORDERS = ["xyz", "yxz", "zyx", "zxy", "yzx", "xzy", "XYZ", "xy", "xxy", "abc", "yz"]


def _axis_py(o):
    out = []
    for i in o:
        for j in o:
            try:
                out.append(
                    f"{i}|{j}|"
                    + text.mat2str(np.array(arrays.axis_order_mapping(i, j), dtype=np.float64))
                )
            except ValueError:
                out.append(f"{i}|{j}|ERR")
    return out


case(
    "axis_order_mapping",
    "out = {}; for i = 1:numel(o); for j = 1:numel(o); try; r = mat2str(axis_order_mapping(o{i}, o{j})); catch; r = 'ERR'; end; out{end + 1} = sprintf('%s|%s|%s', o{i}, o{j}, r); end; end",
    _axis_py,
    kind="list",
    o=lambda: list(AXIS_ORDERS),
)

SHAPES = [
    [3, 4, 1],
    [3, 1, 1, 1],
    [1, 1, 5],
    [1, 5],
    [0, 3],
    [3, 0, 1],
    [1, 1, 1],
    [2, 3, 4, 1, 1],
    [1, 3, 1, 4],
    [5, 1],
    [1, 1, 1, 1, 7],
    [2, 1, 3],
]


def _size_py(shapes):
    out = []
    for s in shapes:
        z = np.zeros(s)
        out.append(
            f"{text.mat2str(np.array(arrays.msize(z), dtype=float))}|{arrays.mndims(z)}|{text.mat2str(np.array(arrays.msize(arrays.msqueeze(z)), dtype=float))}|{arrays.msize(z, 3)}"
        )
    return out


case(
    "size_squeeze",
    "out = cellfun(@(s) sprintf('%s|%d|%s|%d', mat2str(size(zeros(s(:)'))), ndims(zeros(s(:)')), mat2str(size(squeeze(zeros(s(:)')))), size(zeros(s(:)'), 3)), shapes, 'UniformOutput', false);",
    _size_py,
    kind="list",
    shapes=lambda: [list(s) for s in SHAPES],
)

# the field order of inputParser's Results (parameters.json follows it)
case(
    "inputParser.order",
    "ip = inputParser; ip.CaseSensitive = false; for k = 1:numel(names); ip.addParameter(names{k}, k); end; ip.parse(); out = fieldnames(ip.Results);",
    lambda names: sorted(names, key=lambda s: (s.lower(), s)),
    kind="list",
    names=lambda: [
        "zeta",
        "Alpha",
        "beta",
        "BKRemoval",
        "bbox",
        "ChannelPatterns",
        "dz",
        "DSDirName",
        "xyPixelSize",
        "Z",
        "a1",
        "A2",
    ],
)

# -- F12: the mapping between MATLAB's arrays and ours ----------------------------------------
# Inputs are MATLAB-form arrays (NumPy shape = MATLAB size); their role 'vol' gives the Python
# side our (z, y, x) arrays, and out_role 'vol' maps a result that is one of our arrays back for
# the comparison (the MATLAB snippets and references are unchanged).

F12_CLASSES = (
    "double",
    "single",
    "int8",
    "uint8",
    "int16",
    "uint16",
    "int32",
    "uint32",
    "int64",
    "uint64",
    "logical",
)
F12_SHAPES = {
    "2d": (5, 7),
    "3d": (4, 5, 3),
    "4d": (3, 4, 2, 3),
    "5d": (2, 3, 4, 2, 2),
    "3d_y1": (1, 5, 3),
    "3d_x1": (4, 1, 3),
}


def f12_array(cls: str, shape, seed: int) -> np.ndarray:
    """A MATLAB-form array of a class from the bit-stable stream (distinct values)."""
    n = int(np.prod(shape))
    r = _gen.raw(seed, n)
    if cls == "logical":
        a = (r & np.uint64(1)).astype(bool)
    elif cls in ("double", "single"):
        a = ((r >> np.uint64(11)).astype(np.float64) * 2.0**-53 * 200 - 100).astype(cls)
    else:
        # the stream's bytes as values of the class (little-endian), the first n of them
        a = np.frombuffer(r.astype("<u8").tobytes(), dtype=np.dtype(cls).newbyteorder("<"))[
            :n
        ].astype(cls)
    return np.asarray(a).reshape(shape, order="F")


for _cls in F12_CLASSES:
    for _nd, _shape in F12_SHAPES.items():
        if _nd.startswith("3d_") and _cls != "double":
            continue
        _mk = lambda c=_cls, sh=_shape: f12_array(c, sh, SEED + 20)
        case(
            f"f12.permute.{_cls}.{_nd}",
            "out = permute(a, [ndims(a):-1:3, 1, 2]);",
            lambda a: a,  # MATLAB's permute is our array itself
            roles={"a": "vol"},
            a=_mk,
        )
    case(
        f"f12.colon.{_cls}",
        "out = a(:);",
        lambda a: arrays.matlab_order(a).reshape(-1),
        roles={"a": "vol"},
        a=(lambda c=_cls: f12_array(c, F12_SHAPES["4d"], SEED + 21)),
    )

F12_PAIRS = [
    ((4, 5, 3), (4, 5)),
    ((4, 5, 3), (4, 1)),
    ((4, 5, 3), (1, 5)),
    ((4, 5, 3), (1, 1, 3)),
    ((4, 5, 3), (1, 5, 3)),
    ((4, 5, 3), (4, 1, 3)),
    ((4, 1, 3), (1, 5)),
    ((4, 5, 3, 2), (4, 5)),
    ((4, 5, 3, 2), (1, 1, 3, 2)),
    ((4, 5, 3, 2), (1, 1, 1, 2)),
    ((4, 1), (1, 5)),
    ((1, 1, 3), (4, 5)),
]
for _k, (_sa, _sb) in enumerate(F12_PAIRS):
    _ma = lambda sh=_sa: f12_array("double", sh, SEED + 22)
    _mb = lambda sh=_sb: f12_array("double", sh, SEED + 23)
    case(
        f"f12.broadcast.{_k:02d}.plus",
        "out = pkpy_f02_arith(a, b, 'plus');",
        lambda a, b: dtypes.plus(a, b),
        roles={"a": "vol", "b": "vol"},
        out_role="vol",
        a=_ma,
        b=_mb,
    )
    case(
        f"f12.broadcast.{_k:02d}.numpy_times",
        "out = a .* b;",
        lambda a, b: a * b,
        roles={"a": "vol", "b": "vol"},
        out_role="vol",
        a=_ma,
        b=_mb,
    )
for _k in (0, 3, 7):
    _sa, _sb = F12_PAIRS[_k]
    case(
        f"f12.broadcast.{_k:02d}.single_double",
        "out = pkpy_f02_arith(a, b, 'times');",
        lambda a, b: dtypes.times(a, b),
        roles={"a": "vol", "b": "vol"},
        out_role="vol",
        a=(lambda sh=_sa: f12_array("single", sh, SEED + 24)),
        b=(lambda sh=_sb: f12_array("double", sh, SEED + 25)),
    )

CASE_IDS = [c.id for c in CASES]
assert len(set(CASE_IDS)) == len(CASE_IDS), "duplicate case ids"
