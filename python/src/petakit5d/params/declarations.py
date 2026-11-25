"""The lab's arguments-block declarations, applied as MATLAB's arguments block applies them.

The io wrappers ``readtiff``, ``writetiff``, ``readzarr``, ``writezarr`` and ``createzarr``
declare their inputs in ``arguments`` blocks (``range (1, :) {mustBeNumeric, mustBeFinite}``,
``groupWrite (1, 1) logical``). MATLAB R2024b converts a value to the declared class and size
before the function runs, then validates it, in this order (measured with fake functions of each
of the lab's declarations, positional and name-value, on values of every class and many sizes;
``tests/params/test_declarations.py``):

1. **the class** (``char``, ``logical``): char takes char arrays of any size, a string scalar as
   its text; ``logical`` takes logical arrays, and numbers of any size (0 false, any other
   number true); a NaN anywhere, a complex number, a string, a cell, a struct or a function
   handle raises ``MATLAB:validation:UnableToConvert``, as logical values and complex numbers
   do for char;
2. **the size**: ``(1, 1)`` takes one element only; ``(1, :)`` takes a vector of any
   orientation (one dimension other than 1: a column, 1x1x3), made a row, and ``[]`` (0x0),
   made 1x0; anything else raises ``MATLAB:validation:IncompatibleSize`` (a matrix, 0x3, 1x2x3);
3. **the validation functions**, in their order: ``mustBeNumeric`` (complex numbers too; not
   logical values or characters), ``mustBeNumericOrLogical``, ``mustBeFinite`` (NaN and Inf
   refused, of either part of a complex number) and ``mustBeMember``, whose list is
   case-sensitive (``'LZW'`` is refused where ``'lzw'`` is a member); each raises
   ``MATLAB:validators:<name>``.

The MATLAB-name layer and the command line give the function the converted value (a column made
a row, numbers made logical values), and raise MATLAB's errors with its identifiers and messages
(``spec.MatlabValidationError``).

Some of MATLAB's class conversions only hide mistakes. The port refuses them instead, with an
error of its own (``spec.ConversionRefusedError``, identifier
``petakit5d:params:conversionRefused``; decided by the user on 2026-11-22, docs/parity/io.md):

- for ``char``: numbers (MATLAB takes the characters of their codes: ``writetiff(img, 65)``
  writes a file named ``A``, ``[]`` is ``''``), a cell of char arrays (MATLAB makes its texts a
  char matrix, ``{'a.tif'}`` is ``'a.tif'``, ``{}`` is ``''``) and a function handle (its name);
- for ``logical``: characters (each is true but ``char(0)``: ``'a'`` is true, so is ``'0'``).

The snake_case API checks a declaration's classes as before (``spec.validator_alternatives``):
its values are the port's own forms (``range=(2, None)``: to the end). A declaration with a size,
class or validation function this reading does not know gives no ``Declaration`` (``read``
returns None) and is checked by its classes as before; a test fails for a registered one.
"""

from __future__ import annotations

import functools
import math
import re
from collections.abc import Mapping
from typing import Any

import numpy as np

from petakit5d.params import spec as _s
from petakit5d.params import values as _v

__all__ = ["REFUSED_ID", "Declaration", "Value", "read", "value"]

REFUSED_ID = "petakit5d:params:conversionRefused"
_SIZE_ID = "MATLAB:validation:IncompatibleSize"
_CONVERT_ID = "MATLAB:validation:UnableToConvert"
_VALIDATORS = ("mustBeNumeric", "mustBeNumericOrLogical", "mustBeFinite", "mustBeMember")
_CLASSES = ("char", "logical")
_NUMERIC = frozenset(_v.NUMERIC_CLASSES)
_DTYPES = {
    "double": np.float64,
    "single": np.float32,
    "logical": np.bool_,
    "char": "<U1",
    **{c: np.dtype(c) for c in _v.INT_RANGES},
}
_DECL_RE = re.compile(
    r"^\s*(?:\((?P<size>[^()]*)\))?\s*(?P<cls>[A-Za-z]\w*)?\s*(?:\{(?P<val>.*)\})?\s*$", re.S
)
_MEMBER_RE = re.compile(r"^mustBeMember\s*\(\s*[\w.]+\s*,\s*\{(?P<words>[^{}]*)\}\s*\)$", re.S)
_WORD_RE = re.compile(r"\s*'((?:[^']|'')*)'\s*(?:,|$)")


# ----------------------------------------------------------------------------------------------
# reading


def _split_top(text: str) -> list[str]:
    """the comma-separated items of a validator list (commas in brackets and quotes kept)"""
    out, depth, start, quote = [], 0, 0, False
    for i, c in enumerate(text):
        if c == "'":
            quote = not quote
        elif quote:
            continue
        elif c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == "," and depth == 0:
            out.append(text[start:i])
            start = i + 1
    out.append(text[start:])
    return [s.strip() for s in out if s.strip()]


def _validator(item: str) -> tuple[str, tuple[str, ...]] | None:
    name = re.sub(r"\s*\(\s*\)$", "", item)
    if name in _VALIDATORS[:3]:
        return name, ()
    m = _MEMBER_RE.match(item)
    if m is None:
        return None
    words, rest = [], m.group("words")
    while rest.strip():
        w = _WORD_RE.match(rest)
        if w is None:
            return None
        words.append(w.group(1).replace("''", "'"))
        rest = rest[w.end() :]
    return ("mustBeMember", tuple(words)) if words else None


@functools.cache
def read(text: str | None) -> Declaration | None:
    """The declaration of an arguments-block text (``"(1, :) {mustBeNumeric}"``,
    ``"(1, 1) logical"``, ``"char"``); None for an inputParser validator (``@...``), no text,
    or a size (other than ``(1, 1)`` and ``(1, :)``), class or validation function this reading
    does not know."""
    if not text or text.lstrip().startswith("@"):
        return None
    m = _DECL_RE.match(text)
    if m is None:
        return None
    size = None
    if m.group("size") is not None:
        dims = tuple(d.strip() for d in m.group("size").split(","))
        if dims not in (("1", "1"), ("1", ":")):
            return None
        size = (1, 1) if dims[1] == "1" else (1, None)
    cls = m.group("cls")
    if cls is not None and cls not in _CLASSES:
        return None
    validators = []
    if m.group("val") is not None:
        for item in _split_top(m.group("val")):
            v = _validator(item)
            if v is None:
                return None
            validators.append(v)
    if size is None and cls is None and not validators:
        return None
    return Declaration(text, size, cls, tuple(validators))


# ----------------------------------------------------------------------------------------------
# values as the arguments block sees them


def _msize(shape: tuple[int, ...]) -> tuple[int, ...]:
    """MATLAB's size of an array of this NumPy shape (0-D 1x1, 1-D a row, no trailing ones)"""
    s = tuple(int(n) for n in shape)
    if len(s) == 0:
        return (1, 1)
    if len(s) == 1:
        return (1, s[0])
    while len(s) > 2 and s[-1] == 1:
        s = s[:-1]
    return s


class Value:
    """A value as MATLAB's arguments block sees it: ``kind`` ('array', 'cell', 'struct',
    'string', 'handle' or 'other'), MATLAB's ``mclass`` and ``size``; for arrays ``array``, a
    NumPy array of that size (char as '<U1') and ``complex``; for cells their ``items``, for a
    string scalar its text (``items[0]``)."""

    __slots__ = ("kind", "mclass", "size", "array", "complex", "items")

    def __init__(
        self,
        kind: str,
        mclass: str,
        size: tuple[int, ...],
        array: np.ndarray | None = None,
        complex_: bool = False,
        items: list | None = None,
    ):
        self.kind = kind
        self.mclass = mclass
        self.size = size
        self.array = array
        self.complex = complex_
        self.items = items

    def __repr__(self) -> str:
        return f"Value({self.kind}, {self.mclass}, {self.size})"


def _array(a: np.ndarray, mclass: str, complex_: bool = False) -> Value:
    size = _msize(a.shape)
    return Value("array", mclass, size, a.reshape(size), complex_)


def _numpy_value(x: Any) -> Value:
    a = np.asarray(x)
    k = a.dtype.kind
    if k == "b":
        return _array(a, "logical")
    if k in "iu":
        return _array(a, a.dtype.name)
    if k == "f" and a.dtype in (np.float32, np.float64):
        return _array(a, "single" if a.dtype == np.float32 else "double")
    if k == "c" and a.dtype in (np.complex64, np.complex128):
        return _array(a, "single" if a.dtype == np.complex64 else "double", True)
    if k == "U":
        if a.dtype.itemsize <= 4:  # one character per element: a char array of that size
            return _array(a.astype("<U1"), "char")
        if a.ndim <= 1:  # as the registry reads it (spec._text_of): the text of its elements
            return value("".join(a.reshape(-1).tolist()))
    return Value("other", "", (1, 1))


def value(x: Any) -> Value:
    """``x`` in any of the registry's forms (Python values, NumPy arrays, ``MArray``, ``Cell``,
    ``MString``, dicts and struct arrays, callables) as MATLAB sees it."""
    if isinstance(x, _v.MString):
        return Value("string", "string", (1, 1), items=[str(x)])
    if isinstance(x, str):
        size = (1, len(x)) if x else (0, 0)
        return Value("array", "char", size, np.array(list(x), dtype="<U1").reshape(size))
    if isinstance(x, bool):
        return Value("array", "logical", (1, 1), np.full((1, 1), x))
    if isinstance(x, int | float):
        try:
            f = float(x)
        except OverflowError:  # (an int beyond the doubles: MATLAB's literal is Inf)
            f = math.inf if x > 0 else -math.inf
        return Value("array", "double", (1, 1), np.full((1, 1), f))
    if isinstance(x, complex):
        return Value("array", "double", (1, 1), np.full((1, 1), x), True)
    if x is None:
        return Value("array", "double", (0, 0), np.zeros((0, 0)))
    if isinstance(x, _v.MArray):
        a = np.array(x.data, dtype=_DTYPES[x.mclass]).reshape(x.shape, order="F")
        return Value("array", x.mclass, tuple(x.shape), a)
    if isinstance(x, _v.Cell):
        return Value("cell", "cell", tuple(x.shape), items=list(x))
    if type(x).__module__ == "numpy" and hasattr(x, "dtype"):
        return _numpy_value(x)
    if isinstance(x, list | tuple):
        try:
            m = _v.from_python(x)
        except (TypeError, ValueError):
            return Value("other", "", (1, 1))
        return value(m)
    if _s._is_struct(x):
        return Value("struct", "struct", (1, 1) if isinstance(x, Mapping) else _s._size(x))
    if callable(x) and not isinstance(x, type):
        return Value("handle", "function_handle", (1, 1))
    return Value("other", "", (1, 1))


# ----------------------------------------------------------------------------------------------
# the declaration


class Declaration:
    """An arguments-block declaration read whole: ``size`` ((1, 1), (1, None) for (1, :), or
    None), ``cls`` ('char', 'logical' or None) and ``validators`` ((name, words) in their
    order). ``convert`` is MATLAB's conversion and validation of a value, ``apply`` the value
    the function gets."""

    __slots__ = ("text", "size", "cls", "validators")

    def __init__(
        self,
        text: str,
        size: tuple[int, int | None] | None,
        cls: str | None,
        validators: tuple[tuple[str, tuple[str, ...]], ...],
    ):
        self.text = text
        self.size = size
        self.cls = cls
        self.validators = validators

    def __repr__(self) -> str:
        return f"Declaration({self.text!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Declaration) and other.text == self.text

    def __hash__(self) -> int:
        return hash(self.text)

    def convert(self, x: Any, name: str, position: int | None = None) -> Value:
        """``x`` converted to the declared class and size and validated, as MATLAB does it:
        ``MatlabValidationError`` with MATLAB's identifier and message (``position``: the
        argument's place among the positional ones, else it is a name-value argument), or
        ``ConversionRefusedError`` for the conversions the port refuses."""
        return self._convert(value(x), _Where(name, position))

    def _convert(self, v: Value, where: _Where) -> Value:
        # (each step gives back the same Value where it changes nothing)
        if self.cls == "char":
            v = _to_char(v, where)
        elif self.cls == "logical":
            v = _to_logical(v, where)
        if self.size is not None:
            v = _fit(v, self.size, where)
        for fn, words in self.validators:
            _validate(v, fn, words, where)
        return v

    def apply(self, x: Any, name: str, position: int | None = None) -> Any:
        """The value the function gets: ``x`` itself where MATLAB keeps it as it is (a normal
        call is unchanged), else the converted value in the registry's plain forms, as ``x``
        came in: a row of numbers as a tuple, a NumPy vector or an ``MArray`` row; text as
        str; logical values as bool, a tuple, a NumPy array or an ``MArray``."""
        before = value(x)
        v = self._convert(before, _Where(name, position))
        if v is before:
            return x
        if v.kind != "array":
            return _v.Cell(v.items, v.size) if v.kind == "cell" else x
        if v.mclass == "char":
            if len(v.size) == 2 and v.size[0] <= 1:
                return "".join(v.array.reshape(-1).tolist())
            return x
        if v.mclass != before.mclass:  # numbers made logical values
            return _as_logical(x, v)
        return _as_row(x, v)


class _Where:
    """the argument, as MATLAB's messages name it"""

    __slots__ = ("name", "position")

    def __init__(self, name: str, position: int | None):
        self.name = name
        self.position = position

    def matlab(self) -> str:
        if self.position is not None:
            return f"Invalid argument at position {self.position}."
        return f"Invalid value for '{self.name}' argument."

    def port(self) -> str:
        if self.position is not None:
            return f"The argument at position {self.position} ('{self.name}')"
        return f"The value of '{self.name}'"


def _matlab_error(where: _Where, text: str, identifier: str) -> _s.MatlabValidationError:
    return _s.MatlabValidationError(f"{where.matlab()} {text}", identifier)


def _refused(where: _Where, what: str, matlab: str) -> _s.ConversionRefusedError:
    return _s.ConversionRefusedError(
        f"{where.port()} {what}: MATLAB's arguments block would {matlab}, a conversion the "
        "port refuses because it hides mistakes (docs/parity/io.md).",
        REFUSED_ID,
    )


def _is_char_array(e: Any) -> bool:
    v = value(e)
    return v.kind == "array" and v.mclass == "char"


def _to_char(v: Value, where: _Where) -> Value:
    if v.kind == "array":
        if v.mclass == "char":
            return v
        if v.mclass != "logical" and not v.complex:
            raise _refused(
                where,
                "must be text (char), not numbers",
                "take the characters of their codes (65 is 'A', [] is '')",
            )
    elif v.kind == "string":
        return value(v.items[0])  # its text ("" is '', 0x0)
    elif v.kind == "cell" and all(_is_char_array(e) for e in v.items):
        raise _refused(
            where,
            "must be text (char), not a cell",
            "make the cell's texts a char matrix ({'a.tif'} is 'a.tif', {} is '')",
        )
    elif v.kind == "handle":
        raise _refused(
            where,
            "must be text (char), not a function handle",
            "take its name ('sin' for @sin)",
        )
    raise _matlab_error(where, "Value must be of type char or be convertible to char.", _CONVERT_ID)


def _to_logical(v: Value, where: _Where) -> Value:
    if v.kind == "array":
        if v.mclass == "logical":
            return v
        if v.mclass == "char":
            raise _refused(
                where,
                "must be logical (or numbers), not text",
                "take each character as true ('a' is true, and so is '0')",
            )
        a = v.array
        if not v.complex and not (a.dtype.kind == "f" and bool(np.isnan(a).any())):
            return Value("array", "logical", v.size, a != 0)
    raise _matlab_error(
        where, "Value must be of type logical or be convertible to logical.", _CONVERT_ID
    )


def _fit(v: Value, size: tuple[int, int | None], where: _Where) -> Value:
    if size == (1, 1):
        if v.size != (1, 1):
            raise _matlab_error(where, "Value must be a scalar.", _SIZE_ID)
        return v
    s = v.size
    if s == (0, 0):
        row = (1, 0)
    elif sum(d != 1 for d in s) <= 1:
        row = (1, math.prod(s))
    else:
        raise _matlab_error(where, "Value must be a vector.", _SIZE_ID)
    if row == s:
        return v
    if v.kind == "array":
        return Value("array", v.mclass, row, v.array.reshape(row), v.complex)
    return Value(v.kind, v.mclass, row, items=v.items)


def _validate(v: Value, fn: str, words: tuple[str, ...], where: _Where) -> None:
    array = v.kind == "array"
    if fn == "mustBeNumeric":
        if array and v.mclass in _NUMERIC:
            return
        raise _matlab_error(where, "Value must be numeric.", "MATLAB:validators:mustBeNumeric")
    if fn == "mustBeNumericOrLogical":
        if array and (v.mclass in _NUMERIC or v.mclass == "logical"):
            return
        raise _matlab_error(
            where, "Value must be numeric or logical.", "MATLAB:validators:mustBeNumericOrLogical"
        )
    if fn == "mustBeFinite":
        # (the lab's declarations reach it after mustBeNumeric: numbers only)
        if array and v.mclass != "char" and bool(np.isfinite(v.array).all()):
            return
        raise _matlab_error(where, "Value must be finite.", "MATLAB:validators:mustBeFinite")
    # mustBeMember: one of the words, case-sensitive (the lab's take a char row only)
    if array and v.mclass == "char" and "".join(v.array.reshape(-1).tolist()) in words:
        return
    listed = "".join(f"    '{w}'\n" for w in words)
    raise _matlab_error(
        where, f"Value must be a member of this set:\n{listed}", "MATLAB:validators:mustBeMember"
    )


def _flatten(x: Any) -> list:
    out = []
    for e in x:
        out.extend(_flatten(e) if isinstance(e, list | tuple) else [e])
    return out


def _as_row(x: Any, v: Value) -> Any:
    """a vector or [] of numbers or logical values made MATLAB's row, in the form ``x`` came in"""
    if x is None:
        return None
    if isinstance(x, _v.MArray):
        return _v.MArray(x.mclass, v.size, list(x.data))
    if type(x).__module__ == "numpy":
        return np.asarray(x).reshape(-1)
    if isinstance(x, list | tuple):
        return tuple(_flatten(x))
    return x


def _as_logical(x: Any, v: Value) -> Any:
    """numbers made MATLAB's logical values (of MATLAB's size ``v.size``), in the form ``x`` came
    in: one value as bool; a NumPy array or an ``MArray`` as one of class logical (a row 1-D);
    a sequence as a tuple of bools (a matrix as a tuple of rows); an empty one as None"""
    a = v.array
    if a.size == 1:
        return bool(a.reshape(-1)[0])
    if isinstance(x, _v.MArray):
        return _v.MArray("logical", v.size, [bool(e) for e in a.reshape(-1, order="F")])
    row = len(v.size) == 2 and v.size[0] == 1
    if type(x).__module__ == "numpy":
        return a.reshape(-1) if row else a
    if a.size == 0:  # (None, ())
        return None
    if isinstance(x, list | tuple) and len(v.size) == 2:
        rows = tuple(tuple(bool(e) for e in r) for r in a)
        return rows[0] if row else rows
    return a.reshape(-1) if row else a
