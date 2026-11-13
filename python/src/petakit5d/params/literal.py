"""MATLAB literals: a safe parser (no eval) and a formatter.

``parse`` evaluates the subset of MATLAB expressions that parameter values use: numbers
(``1e-3``, ``.5``, ``Inf``, ``NaN``, ``pi``, ``eps``), ``true``/``false``, char vectors
(``'it''s'``), string scalars (``"a"``), matrices (``[1,2;3,4]``, ``[1 -2]``, newlines as row
breaks, ``...`` continuations), cells (``{'a',{'b'}}``), colon ranges, unary and element-wise
arithmetic, transposes, the constructors ``zeros ones true false Inf NaN cell char double
single logical int8 ... uint64`` and ``sind cosd sqrt abs floor ceil round fix``.
Concatenation follows MATLAB's class rules (char wins, then integer classes, single, double;
logical only if all are logical; ``[]`` is dropped). Anything else (variables, other
functions, struct fields, complex numbers) raises ``LiteralError``. MATLAB's answers for the
corpus of literals are in tests/params/data/literal_oracle.json.

``str2num`` and ``matlab_eval`` reproduce what the mcc parser functions do with string
arguments (``str2num(s)`` evaluates ``['[' s ']']`` and gives ``[]`` on any error).
``format_literal`` writes text that MATLAB evaluates back to the same class, size and bits.
"""

from __future__ import annotations

import math
import struct
from typing import Any

from petakit5d.params.values import (
    FLOAT_CLASSES,
    INT_RANGES,
    Cell,
    MArray,
    MString,
    from_python,
    round_half_away,
    to_int_class,
    to_single,
)

__all__ = ["LiteralError", "parse", "str2num", "matlab_eval", "format_literal"]

EPS = 2.0**-52
# MATLAB's NaN (and the x86-64 default NaN of invalid operations): sign bit set
NAN = struct.unpack(">d", bytes.fromhex("fff8000000000000"))[0]


class LiteralError(ValueError):
    """Text that is not a supported MATLAB literal (or that MATLAB rejects too)."""


_CONSTANTS = {
    "true": lambda: MArray.scalar("logical", True),
    "false": lambda: MArray.scalar("logical", False),
    "Inf": lambda: MArray.scalar("double", math.inf),
    "inf": lambda: MArray.scalar("double", math.inf),
    "NaN": lambda: MArray.scalar("double", NAN),
    "nan": lambda: MArray.scalar("double", NAN),
    "pi": lambda: MArray.scalar("double", math.pi),
    "eps": lambda: MArray.scalar("double", EPS),
}

_FILLS = {
    "zeros": ("double", 0.0),
    "ones": ("double", 1.0),
    "true": ("logical", True),
    "false": ("logical", False),
    "Inf": ("double", math.inf),
    "inf": ("double", math.inf),
    "NaN": ("double", NAN),
    "nan": ("double", NAN),
}


def _is_ident_start(c: str) -> bool:
    return c.isascii() and c.isalpha()


def _is_ident_char(c: str) -> bool:
    return c.isascii() and (c.isalnum() or c == "_")


class _Parser:
    def __init__(self, text: str):
        self.s = text
        self.n = len(text)
        self.i = 0
        # True while directly inside [] or {} (whitespace separates elements there)
        self.mstack: list[bool] = [False]

    # -- low level ------------------------------------------------------------------------
    def err(self, msg: str) -> LiteralError:
        return LiteralError(f"{msg} at position {self.i + 1} of {self.s!r}")

    @property
    def in_matrix(self) -> bool:
        return self.mstack[-1]

    def peek(self, k: int = 0) -> str:
        j = self.i + k
        return self.s[j] if j < self.n else ""

    def skip_ws(self, newlines: bool | None = None) -> bool:
        """Skip blanks, ``...`` continuations and comments; newlines too outside matrices."""
        if newlines is None:
            newlines = not self.in_matrix
        start = self.i
        while self.i < self.n:
            c = self.s[self.i]
            if c in " \t\r":
                self.i += 1
            elif c == "\n" and newlines:
                self.i += 1
            elif self.s.startswith("...", self.i):
                j = self.s.find("\n", self.i)
                self.i = self.n if j < 0 else j + 1
            elif c == "%":
                j = self.s.find("\n", self.i)
                self.i = self.n if j < 0 else j
            else:
                break
        return self.i > start

    # -- grammar --------------------------------------------------------------------------
    def parse_all(self) -> Any:
        self.skip_ws(True)
        if self.i >= self.n:
            raise self.err("empty expression")
        v = self.expr()
        self.skip_ws(True)
        while self.peek() in (";", ","):
            self.i += 1
            self.skip_ws(True)
        if self.i < self.n:
            raise self.err(f"unexpected {self.peek()!r}")
        return v

    def expr(self) -> Any:
        return self.colon()

    def colon(self) -> Any:
        a = self.additive()
        save = self.i
        self.skip_ws(False)
        if self.peek() != ":":
            self.i = save
            return a
        self.i += 1
        self.skip_ws(False)
        b = self.additive()
        save = self.i
        self.skip_ws(False)
        if self.peek() == ":":
            self.i += 1
            self.skip_ws(False)
            c = self.additive()
            return _colon(a, b, c)
        self.i = save
        return _colon(a, MArray.scalar("double", 1.0), b)

    def additive(self) -> Any:
        v = self.multiplicative()
        while True:
            save = self.i
            spaced = self.skip_ws(False)
            c = self.peek()
            if c not in ("+", "-"):
                self.i = save
                return v
            if self.in_matrix and spaced and self.peek(1) not in (" ", "\t"):
                # "[1 -2]": a new element
                self.i = save
                return v
            self.i += 1
            self.skip_ws(False)
            rhs = self.multiplicative()
            v = _arith(c, v, rhs)

    def multiplicative(self) -> Any:
        v = self.unary()
        while True:
            save = self.i
            self.skip_ws(False)
            if self.s.startswith((".*", "./"), self.i):
                op = self.s[self.i : self.i + 2]
                self.i += 2
            elif self.peek() in ("*", "/"):
                op = self.peek()
                self.i += 1
            elif self.peek() == "\\" or self.s.startswith(".\\", self.i):
                raise self.err("left division is not supported")
            else:
                self.i = save
                return v
            self.skip_ws(False)
            rhs = self.unary()
            v = _arith(op, v, rhs)

    def unary(self) -> Any:
        c = self.peek()
        if c in ("-", "+", "~") and not self.s.startswith("~=", self.i):
            self.i += 1
            self.skip_ws(False)
            v = self.unary()
            if c == "-":
                return _negate(v)
            if c == "+":
                return _uplus(v)
            return _not(v)
        return self.power()

    def power(self) -> Any:
        v = self.postfix()
        while True:
            save = self.i
            self.skip_ws(False)
            if self.s.startswith(".^", self.i):
                self.i += 2
            elif self.peek() == "^":
                self.i += 1
            else:
                self.i = save
                return v
            self.skip_ws(False)
            rhs = self.power_operand()
            v = _arith("^", v, rhs)

    def power_operand(self) -> Any:
        c = self.peek()
        if c in ("-", "+", "~"):
            self.i += 1
            self.skip_ws(False)
            v = self.power_operand()
            return _negate(v) if c == "-" else _uplus(v) if c == "+" else _not(v)
        return self.postfix()

    def postfix(self) -> Any:
        v = self.primary()
        while True:
            if self.s.startswith(".'", self.i):
                self.i += 2
                v = _transpose(v)
            elif self.peek() == "'":
                self.i += 1
                v = _transpose(v)
            else:
                return v

    def primary(self) -> Any:
        c = self.peek()
        if not c:
            raise self.err("unexpected end")
        if c.isdigit() or (c == "." and self.peek(1).isdigit()):
            return self.number()
        if c == "'":
            return self.char_literal()
        if c == '"':
            return self.string_literal()
        if c == "[":
            return self.matrix("[", "]")
        if c == "{":
            return self.matrix("{", "}")
        if c == "(":
            self.i += 1
            self.mstack.append(False)
            self.skip_ws(True)
            v = self.expr()
            self.skip_ws(True)
            self.mstack.pop()
            if self.peek() != ")":
                raise self.err("missing )")
            self.i += 1
            return v
        if _is_ident_start(c):
            return self.identifier()
        raise self.err(f"unexpected {c!r}")

    def number(self) -> MArray:
        j = self.i
        s = self.s
        while j < self.n and s[j].isdigit():
            j += 1
        if j < self.n and s[j] == "." and not s.startswith((".*", "./", ".^", ".'"), j):
            j += 1
            while j < self.n and s[j].isdigit():
                j += 1
        if j < self.n and s[j] in "eE":
            k = j + 1
            if k < self.n and s[k] in "+-":
                k += 1
            if k < self.n and s[k].isdigit():
                while k < self.n and s[k].isdigit():
                    k += 1
                j = k
        text = s[self.i : j]
        if j < self.n and _is_ident_char(s[j]):
            if s[j] in "ij" and not (j + 1 < self.n and _is_ident_char(s[j + 1])):
                raise self.err("complex numbers are not supported")
            raise self.err(f"bad number {s[self.i : j + 1]!r}")
        self.i = j
        return MArray.scalar("double", float(text))

    def char_literal(self) -> MArray:
        j = self.i + 1
        out = []
        while True:
            if j >= self.n or self.s[j] == "\n":
                raise self.err("unterminated char vector")
            if self.s[j] == "'":
                if j + 1 < self.n and self.s[j + 1] == "'":
                    out.append("'")
                    j += 2
                    continue
                break
            out.append(self.s[j])
            j += 1
        self.i = j + 1
        return MArray.char("".join(out))

    def string_literal(self) -> MString:
        j = self.i + 1
        out = []
        while True:
            if j >= self.n or self.s[j] == "\n":
                raise self.err("unterminated string")
            if self.s[j] == '"':
                if j + 1 < self.n and self.s[j + 1] == '"':
                    out.append('"')
                    j += 2
                    continue
                break
            out.append(self.s[j])
            j += 1
        self.i = j + 1
        return MString("".join(out))

    def identifier(self) -> Any:
        j = self.i
        while j < self.n and _is_ident_char(self.s[j]):
            j += 1
        name = self.s[self.i : j]
        self.i = j
        if self.peek() == "(":
            args = self.call_args()
            return _call(name, args, self)
        if name in _CONSTANTS:
            return _CONSTANTS[name]()
        if name in _CALLABLE:
            return _call(name, [], self)
        if self.peek() == "." and _is_ident_start(self.peek(1)):
            raise self.err(f"field access {name}. is not supported")
        raise self.err(f"unknown name {name!r} (variables and functions are not literals)")

    def call_args(self) -> list[Any]:
        self.i += 1  # (
        self.mstack.append(False)
        args: list[Any] = []
        self.skip_ws(True)
        if self.peek() == ")":
            self.i += 1
            self.mstack.pop()
            return args
        while True:
            self.skip_ws(True)
            args.append(self.expr())
            self.skip_ws(True)
            c = self.peek()
            if c == ",":
                self.i += 1
                continue
            if c == ")":
                self.i += 1
                break
            raise self.err("expected , or )")
        self.mstack.pop()
        return args

    def matrix(self, open_: str, close: str) -> Any:
        self.i += 1
        self.mstack.append(True)
        rows: list[list[Any]] = []
        row: list[Any] = []
        expect_value = True  # nothing since the last separator
        while True:
            spaced = self.skip_ws(False)
            c = self.peek()
            if not c:
                raise self.err(f"missing {close}")
            if c == close:
                self.i += 1
                break
            if c in (";", "\n"):
                self.i += 1
                if row:
                    rows.append(row)
                row = []
                expect_value = True
                continue
            if c == ",":
                if expect_value and not row:
                    raise self.err("unexpected ,")
                self.i += 1
                expect_value = True
                continue
            if not expect_value and not spaced and c not in ("'", '"'):
                raise self.err(f"unexpected {c!r}")
            if c in (")", "]", "}"):
                raise self.err(f"unexpected {c!r}")
            row.append(self.expr())
            expect_value = False
        if row:
            rows.append(row)
        self.mstack.pop()
        if open_ == "{":
            rows = [[Cell([e], (1, 1)) for e in r] for r in rows]
            if not rows:
                return Cell()
        return _vcat([_hcat(r, self) for r in rows], self)


# ----------------------------------------------------------------------------------------------
# operations on model values


def _as_array(v: Any, what: str) -> MArray:
    if isinstance(v, MArray):
        return v
    kind = "cell" if isinstance(v, Cell) else "string" if isinstance(v, MString) else "struct"
    raise LiteralError(f"{what} of a {kind} is not supported")


def _numeric_data(a: MArray) -> list[float]:
    if a.mclass == "char":
        return [float(ord(ch)) for ch in a.data]
    if a.mclass == "logical":
        return [1.0 if x else 0.0 for x in a.data]
    if a.mclass in INT_RANGES:
        raise LiteralError("arithmetic on integer classes is not supported")
    return [float(x) for x in a.data]


def _result_class(*arrays: MArray) -> str:
    return "single" if any(a.mclass == "single" for a in arrays) else "double"


def _finish(mclass: str, shape: tuple[int, ...], data: list[float]) -> MArray:
    if mclass == "single":
        data = [to_single(x) for x in data]
    return MArray(mclass, shape, data)


def _negate(v: Any) -> MArray:
    a = _as_array(v, "unary minus")
    return _finish(_result_class(a), a.shape, [-x for x in _numeric_data(a)])


def _uplus(v: Any) -> MArray:
    a = _as_array(v, "unary plus")
    if a.mclass in FLOAT_CLASSES:
        return a
    return _finish("double", a.shape, _numeric_data(a))


def _not(v: Any) -> MArray:
    a = _as_array(v, "not")
    data = _numeric_data(a)
    if any(math.isnan(x) for x in data):
        raise LiteralError("NaN cannot be converted to logical")
    return MArray("logical", a.shape, [x == 0 for x in data])


def _pow(x: float, y: float) -> float:
    """MATLAB's ``x ^ y`` of doubles (``matlabcompat.elementary.power``: fdlibm's pow)"""
    if x < 0 and math.isfinite(y) and y != math.floor(y):  # complex in MATLAB
        raise LiteralError("complex results are not supported")
    from petakit5d.matlabcompat import elementary

    return float(elementary.power(x, y))


def _div(x: float, y: float) -> float:
    if y == 0:
        if x == 0 or math.isnan(x):
            return NAN
        return math.copysign(math.inf, x) * math.copysign(1.0, y)
    return x / y


_BIN = {
    "+": lambda x, y: x + y,
    "-": lambda x, y: x - y,
    "*": lambda x, y: x * y,
    ".*": lambda x, y: x * y,
    "/": _div,
    "./": _div,
    "^": _pow,
}


def _arith(op: str, u: Any, w: Any) -> MArray:
    a = _as_array(u, f"operator {op}")
    b = _as_array(w, f"operator {op}")
    if op in ("*", "^") and not (a.is_scalar() or b.is_scalar()):
        raise LiteralError(f"matrix operator {op} is not supported (use .{op})")
    if op == "/" and not b.is_scalar():
        raise LiteralError("matrix division is not supported")
    x, y = _numeric_data(a), _numeric_data(b)
    if a.is_scalar() and not b.is_scalar():
        shape, x = b.shape, x * b.numel
    elif b.is_scalar() and not a.is_scalar():
        shape, y = a.shape, y * a.numel
    elif a.shape == b.shape:
        shape = a.shape
    else:
        raise LiteralError(f"sizes {a.shape} and {b.shape} do not match")
    f = _BIN[op]
    mclass = _result_class(a, b)
    if mclass == "single":
        x = [to_single(t) for t in x]
        y = [to_single(t) for t in y]
    return _finish(mclass, shape, [_x86_nan(f(p, q), p, q) for p, q in zip(x, y, strict=True)])


def _x86_nan(r: float, p: float, q: float) -> float:
    """r, or MATLAB's NaN for an invalid operation (Inf - Inf, Inf * 0) of numbers, as x86 gives
    it also where the CPU's default NaN is another (ARM's has the sign bit clear)"""
    return NAN if r != r and p == p and q == q else r


def _transpose(v: Any) -> Any:
    if isinstance(v, MString | dict):
        return v
    if len(v.shape) > 2:
        raise LiteralError("transpose of an N-D array")
    r, c = v.shape
    data = [
        v[i + j * r] if isinstance(v, Cell) else v.data[i + j * r]
        for i in range(r)
        for j in range(c)
    ]
    if isinstance(v, Cell):
        return Cell(data, (c, r))
    return MArray(v.mclass, (c, r), data)


def _scalar_value(v: Any, what: str) -> float:
    a = _as_array(v, what)
    if not a.is_scalar():
        raise LiteralError(f"{what} needs scalars")
    return _numeric_data(a)[0]


def _colon(a: Any, d: Any, b: Any) -> MArray:
    """MATLAB's colon operator for doubles (values symmetric about the midpoint)."""
    for v in (a, d, b):
        if isinstance(v, MArray) and v.is_empty():
            return MArray.empty("double", (1, 0))
    base = _scalar_value(a, "colon")
    inc = _scalar_value(d, "colon")
    limit = _scalar_value(b, "colon")
    if any(isinstance(v, MArray) and v.mclass == "single" for v in (a, d, b)):
        raise LiteralError("colon of single values is not supported")
    if not (math.isfinite(base) and math.isfinite(inc) and math.isfinite(limit)):
        if math.isnan(base) or math.isnan(inc) or math.isnan(limit):
            return MArray.scalar("double", NAN)
        raise LiteralError("colon with Inf is not supported")
    if inc == 0 or (base < limit and inc < 0) or (limit < base and inc > 0):
        return MArray.empty("double", (1, 0))
    tol = 2.0 * EPS * max(abs(base), abs(limit))
    sig = 1.0 if inc > 0 else -1.0
    if base == math.floor(base) and inc == 1:
        n = int(math.floor(limit) - base)
    elif base == math.floor(base) and inc == math.floor(inc):
        n = int(math.floor((limit - base) / inc))
    else:
        n = int(round_half_away((limit - base) / inc))
        if sig * (base + n * inc - limit) > tol:
            n -= 1
    last = base + n * inc
    if sig * (last - limit) >= -tol:
        last = limit
    out = [0.0] * (n + 1)
    for k in range(n // 2 + 1):
        out[k] = base + k * inc
        out[n - k] = last - k * inc
    if n % 2 == 0:
        out[n // 2] = (base + last) / 2
    return MArray("double", (1, n + 1), out)


def _class_of_concat(arrays: list[MArray]) -> str:
    classes = [a.mclass for a in arrays]
    if "char" in classes:
        return "char"
    for c in classes:
        if c in INT_RANGES:
            return c
    if "single" in classes:
        return "single"
    if classes and all(c == "logical" for c in classes):
        return "logical"
    return "double"


def _convert(a: MArray, mclass: str) -> list[Any]:
    if a.mclass == mclass:
        return a.data
    if mclass == "char":
        out = []
        for x in _numeric_data(a) if a.mclass != "char" else a.data:
            if isinstance(x, str):
                out.append(x)
                continue
            if x != math.floor(x) or x < 0 or x > 0x10FFFF:
                raise LiteralError("only whole character codes are supported")
            out.append(chr(int(x)))
        return out
    if mclass in INT_RANGES:
        return [to_int_class(x, mclass) for x in _numeric_data(a)]
    data = _numeric_data(a)
    return [to_single(x) for x in data] if mclass == "single" else data


def _concat(items: list[Any], axis: int, p: _Parser | None) -> Any:
    if not items:
        return MArray.empty()
    if any(isinstance(v, MString) for v in items):
        if len(items) == 1:
            return items[0]
        raise LiteralError("string arrays are not supported")
    if any(isinstance(v, dict) for v in items):
        raise LiteralError("struct arrays are not supported")
    if any(isinstance(v, Cell) for v in items):
        cells = []
        for v in items:
            if isinstance(v, Cell):
                cells.append(v)
            elif v.shape == (0, 0):
                continue
            else:
                cells.append(Cell([v], (1, 1)))
        parts = [c for c in cells if c.shape != (0, 0)]
        if not parts:
            return Cell()
        return _join(parts, axis, None, lambda c: list(c), lambda d, s: Cell(d, s))
    arrays = [_as_array(v, "concatenation") for v in items]
    mclass = _class_of_concat(arrays)
    parts = [a for a in arrays if a.shape != (0, 0)]
    if not parts:
        return MArray.empty(mclass)
    conv = [MArray(mclass, a.shape, _convert(a, mclass)) for a in parts]
    return _join(conv, axis, mclass, lambda a: a.data, lambda d, s: MArray(mclass, s, d))


def _join(parts, axis, mclass, data_of, make):
    for a in parts:
        if len(a.shape) > 2:
            raise LiteralError("N-D concatenation is not supported")
    if len(parts) == 1:
        return parts[0]
    other = 1 - axis
    if len({a.shape[other] for a in parts}) != 1:
        raise LiteralError("Dimensions of arrays being concatenated are not consistent.")
    if axis == 1:  # horizontal: column-major data just append
        data: list[Any] = []
        for a in parts:
            data.extend(data_of(a))
        return make(data, (parts[0].shape[0], sum(a.shape[1] for a in parts)))
    rows_total = sum(a.shape[0] for a in parts)
    ncol = parts[0].shape[1]
    data = []
    for j in range(ncol):
        for a in parts:
            r = a.shape[0]
            d = data_of(a)
            data.extend(d[j * r : (j + 1) * r])
    return make(data, (rows_total, ncol))


def _hcat(items: list[Any], p: _Parser | None = None) -> Any:
    return _concat(items, 1, p)


def _vcat(items: list[Any], p: _Parser | None = None) -> Any:
    if len(items) == 1:
        return items[0]
    return _concat(items, 0, p)


def _dims(args: list[Any], name: str) -> tuple[int, ...]:
    if not args:
        return (1, 1)
    if len(args) == 1:
        a = _as_array(args[0], name)
        vals = _numeric_data(a)
        if a.is_scalar():
            dims = (vals[0], vals[0])
        elif a.is_vector():
            dims = tuple(vals)
        else:
            raise LiteralError(f"bad size for {name}")
    else:
        dims = tuple(_scalar_value(a, name) for a in args)
    out = []
    for x in dims:
        if math.isnan(x) or x != math.floor(x):
            raise LiteralError(f"size arguments of {name} must be integers")
        out.append(max(int(x), 0))
    if len(out) < 2:
        out.append(out[0])
    return tuple(out)


def _cast(name: str, args: list[Any]) -> Any:
    if len(args) != 1:
        raise LiteralError(f"{name} takes one argument here")
    v = args[0]
    if name == "cell":
        raise LiteralError("cell needs sizes")
    if isinstance(v, Cell | MString | dict):
        raise LiteralError(f"{name} of a {type(v).__name__} is not supported")
    if name == "logical":
        data = _numeric_data(v)
        if any(math.isnan(x) for x in data):
            raise LiteralError("NaN cannot be converted to logical")
        return MArray("logical", v.shape, [x != 0 for x in data])
    if name == "char":
        return MArray("char", v.shape, _convert(v, "char"))
    if name in INT_RANGES:
        return MArray(name, v.shape, _convert(v, name))
    if name == "single":
        if v.mclass in INT_RANGES:
            return MArray("single", v.shape, [to_single(float(x)) for x in v.data])
        return MArray("single", v.shape, [to_single(x) for x in _numeric_data(v)])
    if v.mclass in INT_RANGES:
        return MArray("double", v.shape, [float(x) for x in v.data])
    return MArray("double", v.shape, _numeric_data(v))


def _call(name: str, args: list[Any], p: _Parser) -> Any:
    if name in _FILLS:
        mclass, fill = _FILLS[name]
        shape = _dims(args, name)
        return MArray(mclass, shape, [fill] * math.prod(shape))
    if name == "cell":
        if not args:
            return Cell()
        shape = _dims(args, name)
        return Cell([MArray.empty() for _ in range(math.prod(shape))], shape)
    if name in ("double", "single", "logical", "char", *INT_RANGES):
        if not args and name in INT_RANGES:
            return MArray.scalar(name, 0)
        return _cast(name, args)
    if name in _MATH:
        if len(args) != 1:
            raise p.err(f"{name} takes one argument here")
        a = _as_array(args[0], name)
        if a.mclass not in ("double", "logical", "char"):
            raise LiteralError(f"{name} of {a.mclass} is not supported")
        return MArray("double", a.shape, [_MATH[name](x) for x in _numeric_data(a)])
    raise p.err(f"function {name!r} is not supported")


def _sincosd(x: float, cos: bool) -> float:
    """MATLAB's sind and cosd (``matlabcompat.elementary``, MATLAB's bits on every platform)"""
    from petakit5d.matlabcompat import elementary

    return float(elementary.cosd(x) if cos else elementary.sind(x))


def _sqrt(x: float) -> float:
    if x < 0:
        raise LiteralError("complex results are not supported")
    return math.sqrt(x)


def _rounder(f):
    # keeps the sign of a zero result (ceil(-0.5) is -0 in IEEE arithmetic)
    def g(x: float) -> float:
        if not math.isfinite(x):
            return x
        return math.copysign(float(f(x)), x)

    return g


_MATH = {
    "sind": lambda x: _sincosd(x, False),
    "cosd": lambda x: _sincosd(x, True),
    "sqrt": _sqrt,
    "abs": abs,
    "floor": _rounder(math.floor),
    "ceil": _rounder(math.ceil),
    "round": round_half_away,
    "fix": _rounder(math.trunc),
}

_CALLABLE = {"zeros", "ones", "cell"}


# ----------------------------------------------------------------------------------------------
# public functions


def parse(text: str) -> Any:
    """Evaluate a MATLAB literal expression without eval; returns a model value."""
    if not isinstance(text, str):
        raise TypeError("parse takes a str")
    return _Parser(text).parse_all()


def matlab_eval(text: str) -> Any:
    """``eval(text)`` of a literal expression, as the mcc parsers run it on cells."""
    return parse(text)


def str2num(text: str) -> tuple[Any, bool, str]:
    """MATLAB's ``[x, ok, msg] = str2num(text)``: ``['[' text ']']`` evaluated; [] on error."""
    if not text:
        return MArray.empty(), False, ""
    text = text.replace("\0", " ")
    try:
        x = parse("[" + text + "]")
    except LiteralError as err:
        return MArray.empty(), False, str(err)
    if isinstance(x, Cell | MString) or (isinstance(x, MArray) and x.mclass == "char"):
        # str2num gives [] for char, cell and string results
        return MArray.empty(), False, "not a number"
    return x, True, ""


# ----------------------------------------------------------------------------------------------
# formatting


def _fmt_double(x: float) -> str:
    if math.isnan(x):
        # MATLAB's NaN has the sign bit set; the other one is -NaN
        return "NaN" if math.copysign(1.0, x) < 0 else "-NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "-0" if math.copysign(1.0, x) < 0 else "0"
    if x == math.floor(x) and abs(x) < 1e15:
        return str(int(x))
    return repr(x)


def _fmt_char_row(text: str) -> str:
    if all(c.isprintable() for c in text):
        return "'" + text.replace("'", "''") + "'"
    parts: list[str] = []
    run = ""
    for c in text:
        if c.isprintable():
            run += c
        else:
            if run:
                parts.append("'" + run.replace("'", "''") + "'")
                run = ""
            parts.append(f"char({ord(c)})")
    if run:
        parts.append("'" + run.replace("'", "''") + "'")
    return "[" + ",".join(parts) + "]"


def _fmt_elem(x: Any, mclass: str) -> str:
    if mclass == "logical":
        return "true" if x else "false"
    if mclass in INT_RANGES:
        return str(int(x))
    if mclass == "single":
        return _fmt_double(x)
    return _fmt_double(x)


def _fmt_array(a: MArray) -> str:
    if len(a.shape) > 2:
        raise ValueError("N-D arrays have no literal form")
    r, c = a.shape
    if a.mclass == "char":
        if r == 0 and c == 0:
            return "''"
        if r == 1:
            return _fmt_char_row(a.text())
        if a.is_empty():
            return f"char(zeros({r},{c}))"
        return "[" + ";".join(_fmt_char_row("".join(row)) for row in a.rows()) + "]"
    if a.is_empty():
        if a.mclass == "logical":
            return f"false({r},{c})"
        body = "[]" if (r, c) == (0, 0) else f"zeros({r},{c})"
        return body if a.mclass == "double" else f"{a.mclass}({body})"
    if a.is_scalar():
        body = _fmt_elem(a.data[0], a.mclass)
    else:
        body = (
            "[" + ";".join(",".join(_fmt_elem(x, a.mclass) for x in row) for row in a.rows()) + "]"
        )
    if a.mclass in ("double", "logical"):
        return body
    return f"{a.mclass}({body})"


def format_literal(value: Any, hint: str | None = None) -> str:
    """MATLAB text that evaluates (``eval``) to the same class, size and bits.

    Plain Python values are converted with ``values.from_python`` first (``hint`` as there).
    """
    v = from_python(value, hint)
    if isinstance(v, MString):
        return '"' + str(v).replace('"', '""') + '"'
    if isinstance(v, dict):
        raise ValueError("structs have no literal form here")
    if isinstance(v, Cell):
        if len(v.shape) > 2:
            raise ValueError("N-D cells have no literal form")
        r, c = v.shape
        if not len(v):
            return "{}" if (r, c) == (0, 0) else f"cell({r},{c})"
        return "{" + ";".join(",".join(format_literal(e) for e in row) for row in v.rows()) + "}"
    return _fmt_array(v)
