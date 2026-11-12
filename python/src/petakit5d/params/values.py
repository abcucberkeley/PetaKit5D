"""MATLAB values as the literal parser, the formatter and the registry see them.

The model is exact (class, size, bits) and needs no NumPy:

- ``MArray(mclass, shape, data)``: a numeric, logical or char array; ``data`` in column-major
  order (floats for double and single, ints for the integer classes, bools, one-character
  strings for char);
- ``Cell``: a cell array (a list in column-major order with a ``shape``);
- ``MString``: a string scalar (``"abc"``); a Python ``str`` is a char row vector;
- ``dict``: a scalar struct.

``to_python`` turns a model value into the plain Python form the snake_case API takes (str,
bool, int or float, tuples, None for an empty numeric array), ``from_python`` goes back, and
``describe`` gives the canonical JSON form that the MATLAB side of the oracle prints too.
"""

from __future__ import annotations

import math
import struct
from collections.abc import Mapping, Sequence
from typing import Any

INT_RANGES = {
    "int8": (-(2**7), 2**7 - 1),
    "uint8": (0, 2**8 - 1),
    "int16": (-(2**15), 2**15 - 1),
    "uint16": (0, 2**16 - 1),
    "int32": (-(2**31), 2**31 - 1),
    "uint32": (0, 2**32 - 1),
    "int64": (-(2**63), 2**63 - 1),
    "uint64": (0, 2**64 - 1),
}
FLOAT_CLASSES = ("double", "single")
NUMERIC_CLASSES = (*FLOAT_CLASSES, *INT_RANGES)
ARRAY_CLASSES = (*NUMERIC_CLASSES, "logical", "char")

_MAX_EXACT_INT = 2**53


def to_single(x: float) -> float:
    """Round a double to the nearest single (as MATLAB's single() does), kept as a float."""
    if math.isnan(x) or math.isinf(x):
        return x
    try:
        return struct.unpack("<f", struct.pack("<f", x))[0]
    except OverflowError:
        return math.copysign(math.inf, x)


def round_half_away(x: float) -> float:
    """MATLAB's round: halves away from zero."""
    if math.isnan(x) or math.isinf(x):
        return x
    r = math.floor(abs(x) + 0.5)
    # floor(|x| + 0.5) is off by one for the largest double below 0.5
    if r - abs(x) > 0.5:
        r -= 1
    return math.copysign(r, x)


def to_int_class(x: float, mclass: str) -> int:
    """MATLAB's conversion of a double to an integer class: round half away, saturate, NaN 0."""
    lo, hi = INT_RANGES[mclass]
    if math.isnan(x):
        return 0
    if math.isinf(x):
        return hi if x > 0 else lo
    r = int(round_half_away(x))
    return min(max(r, lo), hi)


class MArray:
    """A MATLAB numeric, logical or char array (column-major data)."""

    __slots__ = ("mclass", "shape", "data")

    def __init__(self, mclass: str, shape: tuple[int, ...], data: Sequence[Any]):
        if mclass not in ARRAY_CLASSES:
            raise ValueError(f"unknown MATLAB class {mclass!r}")
        shape = tuple(int(n) for n in shape)
        if len(shape) < 2:
            raise ValueError("MATLAB arrays have at least two dimensions")
        while len(shape) > 2 and shape[-1] == 1:
            shape = shape[:-1]
        if math.prod(shape) != len(data):
            raise ValueError(f"{len(data)} values for size {shape}")
        self.mclass = mclass
        self.shape = shape
        self.data = list(data)

    @classmethod
    def scalar(cls, mclass: str, value: Any) -> MArray:
        return cls(mclass, (1, 1), [value])

    @classmethod
    def empty(cls, mclass: str = "double", shape: tuple[int, ...] = (0, 0)) -> MArray:
        return cls(mclass, shape, [])

    @classmethod
    def char(cls, text: str) -> MArray:
        """A char row vector; '' is 0x0 as in MATLAB."""
        return cls("char", (1, len(text)) if text else (0, 0), list(text))

    @property
    def numel(self) -> int:
        return len(self.data)

    def is_empty(self) -> bool:
        return not self.data

    def is_scalar(self) -> bool:
        return len(self.data) == 1

    def is_vector(self) -> bool:
        return len(self.shape) == 2 and (self.shape[0] == 1 or self.shape[1] == 1)

    def rows(self) -> list[list[Any]]:
        """The rows of a 2-D array."""
        r, c = self.shape[0], math.prod(self.shape[1:])
        return [[self.data[i + j * r] for j in range(c)] for i in range(r)]

    def text(self) -> str:
        """The characters of a char array in column-major order."""
        return "".join(self.data)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MArray):
            return NotImplemented
        return (
            self.mclass == other.mclass
            and self.shape == other.shape
            and all(_same_element(a, b) for a, b in zip(self.data, other.data, strict=True))
        )

    def __hash__(self) -> int:
        return hash((self.mclass, self.shape, len(self.data)))

    def __repr__(self) -> str:
        if self.mclass == "char" and len(self.shape) == 2 and self.shape[0] <= 1:
            return f"MArray.char({self.text()!r})"
        return f"MArray({self.mclass!r}, {self.shape!r}, {self.data!r})"


class Cell(list):
    """A MATLAB cell array: a list in column-major order with a shape."""

    __slots__ = ("shape",)

    def __init__(self, items: Sequence[Any] = (), shape: tuple[int, ...] | None = None):
        super().__init__(items)
        if shape is None:
            shape = (1, len(self)) if len(self) else (0, 0)
        shape = tuple(int(n) for n in shape)
        while len(shape) > 2 and shape[-1] == 1:
            shape = shape[:-1]
        if math.prod(shape) != len(self):
            raise ValueError(f"{len(self)} elements for cell size {shape}")
        self.shape = shape

    def rows(self) -> list[list[Any]]:
        r, c = self.shape[0], math.prod(self.shape[1:])
        return [[self[i + j * r] for j in range(c)] for i in range(r)]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Cell):
            return NotImplemented
        return self.shape == other.shape and all(
            same(a, b) for a, b in zip(self, other, strict=True)
        )

    def __ne__(self, other: object) -> bool:
        eq = self.__eq__(other)
        return eq if eq is NotImplemented else not eq

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f"Cell({list(self)!r}, shape={self.shape!r})"


class MString(str):
    """A MATLAB string scalar ("abc"); a plain str is a char row vector."""

    __slots__ = ()

    def __repr__(self) -> str:
        return f"MString({str(self)!r})"


def _same_element(a: Any, b: Any) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        return a == b and math.copysign(1.0, a) == math.copysign(1.0, b)
    return type(a) is type(b) and a == b


def same(a: Any, b: Any) -> bool:
    """Exact equality of two model values (class, size and bits; NaN equals NaN)."""
    if isinstance(a, MString) or isinstance(b, MString):
        return isinstance(a, MString) and isinstance(b, MString) and str(a) == str(b)
    if isinstance(a, dict) and isinstance(b, dict):
        return list(a) == list(b) and all(same(a[k], b[k]) for k in a)
    if type(a) is not type(b):
        return False
    return a == b


# ----------------------------------------------------------------------------------------------
# describe: the canonical JSON form shared with the MATLAB side of the oracle


def float_hex(x: float, mclass: str = "double") -> str:
    """The bits of a double (16 hex digits) or a single (8), as MATLAB's num2hex prints them."""
    if mclass == "single":
        return struct.pack(">f", x).hex()
    return struct.pack(">d", x).hex()


def describe(v: Any) -> dict[str, Any]:
    """Canonical description: {"c": class, "s": size, "d": data} (column-major data)."""
    v = from_python(v)
    if isinstance(v, MString):
        return {"c": "string", "s": [1, 1], "d": [str(v)]}
    if isinstance(v, Cell):
        return {"c": "cell", "s": list(v.shape), "d": [describe(e) for e in v]}
    if isinstance(v, dict):
        return {"c": "struct", "s": [1, 1], "f": {k: describe(e) for k, e in v.items()}}
    if v.mclass in FLOAT_CLASSES:
        data: list[Any] = [float_hex(x, v.mclass) for x in v.data]
    elif v.mclass == "logical":
        data = [int(x) for x in v.data]
    elif v.mclass == "char":
        data = [v.text()]
    else:
        data = [int(x) for x in v.data]
    return {"c": v.mclass, "s": list(v.shape), "d": data}


# ----------------------------------------------------------------------------------------------
# conversion to and from plain Python values


def _py_number(x: Any, mclass: str) -> Any:
    if mclass == "logical":
        return bool(x)
    if mclass in INT_RANGES:
        return int(x)
    if mclass == "char":
        return x
    if math.isfinite(x) and x == math.floor(x) and abs(x) < _MAX_EXACT_INT:
        # an integer-valued double is an int in Python (counts, sizes, indices); -0 stays 0
        return int(x)
    return float(x)


def to_python(v: Any) -> Any:
    """A model value as the snake_case API takes it.

    char -> str (a char matrix -> tuple of row strings); string -> str; logical, numeric
    scalars -> bool, int (integer-valued) or float; empty numeric or logical arrays -> None;
    vectors -> tuple; matrices -> tuple of row tuples; cells -> tuple (rows likewise);
    structs -> dict. Plain Python values pass through (lists become tuples).
    """
    if isinstance(v, MString):
        return str(v)
    if isinstance(v, Cell):
        if len(v.shape) == 2 and (v.shape[0] <= 1 or v.shape[1] == 1):
            return tuple(to_python(e) for e in v)
        return tuple(tuple(to_python(e) for e in row) for row in v.rows())
    if isinstance(v, MArray):
        if v.mclass == "char":
            if len(v.shape) == 2 and v.shape[0] <= 1:
                return v.text()
            return tuple("".join(row) for row in v.rows())
        if v.is_empty():
            return None
        vals = [_py_number(x, v.mclass) for x in v.data]
        if v.is_scalar():
            return vals[0]
        if v.is_vector():
            return tuple(vals)
        if len(v.shape) > 2:
            raise ValueError(f"no plain Python form for a {v.shape} array")
        r = v.shape[0]
        return tuple(tuple(vals[i + j * r] for j in range(v.shape[1])) for i in range(r))
    if isinstance(v, dict):
        return {k: to_python(e) for k, e in v.items()}
    if isinstance(v, list | tuple):
        return tuple(to_python(e) for e in v)
    return v


def _is_number(x: Any) -> bool:
    return isinstance(x, int | float) and not isinstance(x, bool)


def _numpy_to_model(x: Any) -> MArray | None:
    import numpy as np

    a = np.asarray(x)
    kind = a.dtype.kind
    if kind == "b":
        mclass = "logical"
    elif kind == "f":
        if a.dtype == np.float32:
            mclass = "single"
        elif a.dtype == np.float64:
            mclass = "double"
        else:
            raise ValueError(f"no MATLAB class for {a.dtype}")
    elif kind in "iu":
        mclass = f"{'u' if kind == 'u' else ''}int{a.dtype.itemsize * 8}"
    elif kind == "U":
        if a.ndim <= 1:
            return MArray.char("".join(a.reshape(-1).tolist()))
        raise ValueError("no MATLAB class for a string array")
    else:
        raise ValueError(f"no MATLAB class for {a.dtype}")
    shape = a.shape if a.ndim >= 2 else (1, a.size) if a.ndim == 1 else (1, 1)
    data = a.reshape(-1, order="F").tolist()
    if mclass in FLOAT_CLASSES:
        data = [float(d) for d in data]
    return MArray(mclass, shape, data)


def from_python(x: Any, hint: str | None = None) -> Any:
    """A Python value as a model value (MString, Cell and MArray pass through).

    str -> char; bool -> logical; int, float -> double; None -> [] (0x0 double); a sequence
    of bools -> logical row; of numbers (bools allowed among them) -> double row; of equal
    length number sequences -> double matrix; else -> cell row. ``hint="cell"`` makes every
    sequence (and an empty one) a cell; ``hint="char"`` makes a sequence of str a char matrix.
    NumPy arrays and scalars keep their class; dicts are structs.
    """
    if isinstance(x, MArray | Cell | MString):
        return x
    if isinstance(x, str):
        return MArray.char(x)
    if isinstance(x, bool):
        return MArray.scalar("logical", x)
    if _is_number(x):
        return MArray.scalar("double", float(x))
    if x is None:
        return MArray.empty()
    if isinstance(x, Mapping):
        return {str(k): from_python(e) for k, e in x.items()}
    if type(x).__module__ == "numpy":
        return _numpy_to_model(x)
    if isinstance(x, list | tuple):
        if hint == "cell":
            return Cell([from_python(e) for e in x])
        if not x:
            return MArray.empty()
        if all(isinstance(e, bool) for e in x):
            return MArray("logical", (1, len(x)), list(x))
        if all(isinstance(e, bool) or _is_number(e) for e in x):
            return MArray("double", (1, len(x)), [float(e) for e in x])
        if (
            all(isinstance(e, list | tuple) and e for e in x)
            and len({len(e) for e in x}) == 1
            and all(isinstance(f, bool) or _is_number(f) for e in x for f in e)
        ):
            r, c = len(x), len(x[0])
            allbool = all(isinstance(f, bool) for e in x for f in e)
            data = [x[i][j] for j in range(c) for i in range(r)]
            if allbool:
                return MArray("logical", (r, c), data)
            return MArray("double", (r, c), [float(d) for d in data])
        return Cell([from_python(e) for e in x])
    raise TypeError(f"no MATLAB value for a {type(x).__name__}")


def matlab_class(v: Any) -> str:
    """MATLAB's class() of a value (model or plain Python)."""
    if callable(v) and not isinstance(v, type):
        return "function_handle"
    m = from_python(v)
    if isinstance(m, MString):
        return "string"
    if isinstance(m, Cell):
        return "cell"
    if isinstance(m, dict):
        return "struct"
    return m.mclass


def matlab_size(v: Any) -> tuple[int, ...]:
    """MATLAB's size() of a value (model or plain Python)."""
    m = from_python(v)
    if isinstance(m, MString | dict):
        return (1, 1)
    return m.shape
