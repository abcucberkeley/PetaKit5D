"""The lab's validator texts read whole: the conditions beyond their class predicates.

``spec.validator_alternatives`` reads a validator text as alternatives of class predicates; the
command-line conversion and the cell hint come from them, and a text of class predicates only is
checked by them. A text with other terms (``numel(x) <= 2``, ``size(x, 2) == 5``,
``all(size(x) == [3, 2])``, ``any(strcmpi(x, {'cubic', 'linear'}))``,
``ismember(lower(x), {...})``, ``x > 0 && round(x) == x``) is read here as an expression, and its
check is the expression as MATLAB evaluates it (``read``). What MATLAB R2024b does, measured
with inputParser and each of the lab's texts on values of every class and of many sizes
(``tests/params/test_conditions.py``):

- inputParser takes a value only when the validator returns one true element: a logical array
  (``ismember`` of a cell, even all true) is a refusal;
- ``&&`` and ``||`` evaluate their right side only when needed, and an operand of other than one
  element raises ``MATLAB:nonLogicalConditional``;
- the comparisons, ``&`` and ``|`` work element by element with implicit expansion; sizes that
  do not fit raise ``MATLAB:sizeDimensionsMustMatch`` (``size(x) == [3, 2]`` of an n-D array);
- ``x`` compared with a number: numbers, logicals and characters by their values; a string
  raises ``MATLAB:string:ComparisonNotDefined``, a cell, a struct or a function handle
  ``MATLAB:UndefinedFunction``; ``round`` is MATLAB's (half away from zero, exact);
- ``strcmp`` and ``strcmpi`` of ``x`` and a row of N words: a char row (``''`` too) is one
  text; a char matrix of R > 1 rows compares its rows, trailing blanks removed, with the words
  in turn when R is N, and raises ``MATLAB:strcmp:RowsElementsMismatch`` otherwise; a string
  is one text; a cell of one element or of the words' size compares element by element (an
  element that is no char row is unequal), another cell raises
  ``MATLAB:strcmp:InputsSizeMismatch``; anything else (a char array of more dimensions too) is
  unequal;
- ``ismember`` of ``x`` (or ``lower(x)``) in a row of words: a char row (``''`` too) or a string
  is one element; a char matrix gives one element per row and a cell of char arrays one per
  element (so only a cell of one element can pass); a char array of more dimensions raises
  ``MATLAB:cellstr:InputShape``, a cell with an element that is no char array
  ``MATLAB:ISMEMBER:InputClass``, as a number or a logical does, and a struct or a function
  handle ``MATLAB:ISMEMBER:InvalidInputsDataType``;
- ``lower`` changes text (char, string, a cell of char arrays) and passes other values on; a cell
  with an element that is no char array raises ``MATLAB:lower:CellsMustContainChars``.

An error raised by a term reaches the caller as ``spec.MatlabCheckError`` with MATLAB's
identifier, as inputParser passes it on. A text with a construct this reading does not know
gives no condition (``read`` returns None): its class predicates are checked and the rest is not,
as before; the tests fail for a registered text of that kind.
"""

from __future__ import annotations

import functools
import re
from typing import Any

import numpy as np

from petakit5d.params import spec as _s
from petakit5d.params import values as _v

__all__ = ["Condition", "read"]

# the functions this reading knows: the class predicates (spec._pred), and these
_CLASS = frozenset(_s._KNOWN_PREDS - {"function_handle"})
_FUNCTIONS = frozenset(
    {"isa", "numel", "size", "all", "any", "strcmp", "strcmpi", "ismember", "lower", "round"}
)
_COMPARE = ("==", "~=", "<=", ">=", "<", ">")
_TOKEN = re.compile(
    r"\s*(?:(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)|(?P<name>[A-Za-z]\w*)"
    r"|(?P<str>'(?:[^']|'')*')|(?P<op>\|\||&&|==|~=|<=|>=|[<>|&~()\[\]{},]))"
)
# a text of class predicates joined by || && | & only: the class alternatives check it
_ONLY_CLASS = re.compile(
    r"^@\s*\(\s*(\w+)\s*\)"
    r"(?:[\s()|&]|is\w+\s*\(\s*\1\s*\)|isa\s*\(\s*\1\s*,\s*'function_handle'\s*\))*$"
)


class _Unread(Exception):
    """a construct this reading does not know"""


class _Raised(Exception):
    """MATLAB raises an error with this identifier"""

    def __init__(self, identifier: str):
        super().__init__(identifier)
        self.identifier = identifier


# ----------------------------------------------------------------------------------------------
# reading


def _tokens(body: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    i = 0
    while body[i:].strip():
        m = _TOKEN.match(body, i)
        if m is None:
            raise _Unread(body[i:])
        kind = m.lastgroup or ""
        # a quote right after a name, a number or a closing bracket is MATLAB's transpose
        if (
            kind == "str"
            and out
            and m.start(kind) == i
            and (out[-1][0] in ("name", "num") or out[-1][1] in (")", "]", "}"))
        ):
            raise _Unread("transpose")
        out.append((kind, m.group(kind)))
        i = m.end()
    return out


class _Parser:
    """MATLAB's precedence: || below &&, below |, below &, below the comparisons, below ~"""

    def __init__(self, tokens: list[tuple[str, str]], arg: str):
        self.t = tokens
        self.i = 0
        self.arg = arg

    def peek(self) -> str | None:
        return self.t[self.i][1] if self.i < len(self.t) else None

    def take(self, want: str | None = None) -> tuple[str, str]:
        if self.i >= len(self.t) or (want is not None and self.t[self.i][1] != want):
            raise _Unread(f"expected {want}")
        self.i += 1
        return self.t[self.i - 1]

    def expr(self) -> tuple:
        return self._chain("||", "or", self.andand)

    def andand(self) -> tuple:
        return self._chain("&&", "and", self.elor)

    def elor(self) -> tuple:
        return self._chain("|", "|", self.eland)

    def eland(self) -> tuple:
        return self._chain("&", "&", self.compare)

    def _chain(self, op: str, node: str, sub: Any) -> tuple:
        a = sub()
        while self.peek() == op:
            self.take()
            a = (node, a, sub())
        return a

    def compare(self) -> tuple:
        a = self.unary()
        while self.peek() in _COMPARE:
            op = self.take()[1]
            a = ("cmp", op, a, self.unary())
        return a

    def unary(self) -> tuple:
        if self.peek() == "~":
            self.take()
            return ("not", self.unary())
        return self.primary()

    def primary(self) -> tuple:
        kind, text = self.take()
        if kind == "num":
            return ("num", float(text))
        if kind == "str":
            return ("text", text[1:-1].replace("''", "'"))
        if text == "(":
            e = self.expr()
            self.take(")")
            return e
        if text == "[":
            nums = []
            while self.peek() != "]":
                k, n = self.take()
                if k != "num":
                    raise _Unread("[...] of numbers only")
                nums.append(float(n))
                if self.peek() == ",":
                    self.take()
            self.take("]")
            return ("row", tuple(nums))
        if text == "{":
            words = []
            while self.peek() != "}":
                k, w = self.take()
                if k != "str":
                    raise _Unread("{...} of char rows only")
                words.append(w[1:-1].replace("''", "'"))
                if self.peek() == ",":
                    self.take()
            self.take("}")
            return ("words", tuple(words))
        if kind != "name":
            raise _Unread(text)
        if text == self.arg:
            return ("arg",)
        if self.peek() != "(":
            raise _Unread(text)
        self.take("(")
        args = []
        while self.peek() != ")":
            args.append(self.expr())
            if self.peek() == ",":
                self.take()
        self.take(")")
        return _call(text, args)


_BOOLEAN = ("or", "and", "|", "&", "cmp", "not")


def _boolean(n: tuple) -> bool:
    """whether a node gives a logical value (the operands of && and || and the text's result)"""
    return n[0] in _BOOLEAN or (
        n[0] == "call" and n[1] in (*_CLASS, "isa", "all", "any", "strcmp", "strcmpi", "ismember")
    )


def _numeric(n: tuple) -> bool:
    """whether a node gives numbers (the operands of the comparisons)"""
    return n[0] in ("arg", "num", "row") or (n[0] == "call" and n[1] in ("numel", "size", "round"))


def _call(name: str, args: list[tuple]) -> tuple:
    """a call node, its arguments checked against what this reading knows"""
    if name not in _CLASS and name not in _FUNCTIONS:
        raise _Unread(name)
    arg = [a[0] == "arg" for a in args]
    if name in _CLASS or name in ("numel", "lower"):
        ok = arg == [True]
    elif name == "isa":
        ok = arg == [True, False] and args[1] == ("text", "function_handle")
    elif name == "size":
        ok = arg == [True] or (
            arg == [True, False]
            and args[1][0] == "num"
            and args[1][1] >= 1
            and args[1][1] == int(args[1][1])
        )
    elif name in ("strcmp", "strcmpi"):
        ok = arg == [True, False] and args[1][0] == "words"
    elif name == "ismember":
        ok = (
            len(args) == 2
            and (args[0] == ("arg",) or args[0] == ("call", "lower", (("arg",),)))
            and args[1][0] == "words"
        )
    elif name == "round":
        ok = len(args) == 1 and _numeric(args[0])
    else:  # all, any
        ok = len(args) == 1 and _boolean(args[0])
    if not ok:
        raise _Unread(f"{name}(...)")
    return ("call", name, tuple(args))


def _conditions(n: tuple) -> bool:
    """whether a node has a term other than a class predicate (a negation counts: the class
    alternatives have none)"""
    if n[0] in ("or", "and", "|", "&"):
        return _conditions(n[1]) or _conditions(n[2])
    return not (n[0] == "call" and (n[1] in _CLASS or n[1] == "isa"))


class Condition:
    """A validator text read whole; ``outcome(x)`` is its check of a value."""

    __slots__ = ("text", "_node")

    def __init__(self, text: str, node: tuple):
        self.text = text
        self._node = node

    def outcome(self, x: Any) -> tuple[bool, str | None]:
        """(taken, identifier): taken when the validator gives one true element; refused with
        MATLAB's identifier where a term raises an error (None for a plain refusal)"""
        try:
            r = _eval(self._node, x)
        except _Raised as err:
            return False, err.identifier
        return bool(r.size == 1 and r.reshape(-1)[0]), None

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Condition) and other.text == self.text

    def __hash__(self) -> int:
        return hash(self.text)

    def __repr__(self) -> str:
        return f"Condition({self.text!r})"


@functools.cache
def read(text: str | None) -> Condition | None:
    """The condition of an anonymous validator text with terms beyond class predicates; None for
    a text of class predicates only (the class alternatives check it), for one that is no
    anonymous function of one argument, and for one with a construct this reading does not
    know."""
    if not text:
        return None
    v = text.strip()
    if _ONLY_CLASS.match(v):
        return None
    m = re.match(r"^@\s*\(\s*(\w+)\s*\)(.+)$", v, re.S)
    if not m:
        return None
    try:
        p = _Parser(_tokens(m.group(2)), m.group(1))
        node = p.expr()
        if p.i != len(p.t) or not _boolean(node) or not _conditions(node):
            return None
        _check_operands(node)
    except _Unread:
        return None
    return Condition(text, node)


def _check_operands(n: tuple) -> None:
    """&&, ||, &, | and ~ take logical values, the comparisons numbers or logical values
    (anything else is not read)"""
    if n[0] in ("or", "and", "|", "&", "not"):
        subs = n[1:]
        if not all(_boolean(s) for s in subs):
            raise _Unread(f"{n[0]} of logical values only")
    elif n[0] == "cmp":
        subs = n[2:]
        if not all(_numeric(s) or _boolean(s) for s in subs):
            raise _Unread("comparisons of numbers and logical values only")
    elif n[0] == "call":
        subs = tuple(a for a in n[2] if a[0] not in ("text", "words"))
    else:
        subs = ()
    for s in subs:
        _check_operands(s)


# ----------------------------------------------------------------------------------------------
# evaluating


class _Arg:
    """the value checked (or lower() of it)"""

    __slots__ = ("value",)

    def __init__(self, value: Any):
        self.value = value


def _eval(n: tuple, x: Any) -> Any:
    kind = n[0]
    if kind == "arg":
        return _Arg(x)
    if kind == "num":
        return np.array([[n[1]]])
    if kind == "row":
        return np.array([list(n[1])], dtype=np.float64).reshape(1, len(n[1]))
    if kind == "or":
        return np.array([[_scalar(_eval(n[1], x)) or _scalar(_eval(n[2], x))]])
    if kind == "and":
        return np.array([[_scalar(_eval(n[1], x)) and _scalar(_eval(n[2], x))]])
    if kind in ("|", "&"):
        a, b = _expand(_logical(_eval(n[1], x)), _logical(_eval(n[2], x)))
        return (a | b) if kind == "|" else (a & b)
    if kind == "not":
        return ~_logical(_eval(n[1], x))
    if kind == "cmp":
        a, b = _expand(_numbers(_eval(n[2], x)), _numbers(_eval(n[3], x)))
        with np.errstate(invalid="ignore"):
            return {
                "==": np.equal,
                "~=": np.not_equal,
                "<": np.less,
                "<=": np.less_equal,
                ">": np.greater,
                ">=": np.greater_equal,
            }[n[1]](a, b)
    name, args = n[1], n[2]
    if name in _CLASS:
        return np.array([[_s._pred(name, x)]])
    if name == "isa":
        return np.array([[_s._pred("function_handle", x)]])
    if name == "numel":
        return np.array([[float(_s._numel(x))]])
    if name == "size":
        sz = _s._size(x)
        if len(args) == 1:
            return np.array([[float(d) for d in sz]])
        d = int(args[1][1])
        return np.array([[float(sz[d - 1] if d <= len(sz) else 1)]])
    if name in ("all", "any"):
        return _reduce(name, _logical(_eval(args[0], x)))
    if name in ("strcmp", "strcmpi"):
        return _strcmp(x, args[1][1], name == "strcmpi")
    if name == "ismember":
        a = _eval(args[0], x)
        return _ismember(a.value, args[1][1])
    if name == "lower":
        return _Arg(_lower(x))
    if name == "round":
        from petakit5d.core.dtypes import mround

        return mround(_numbers(_eval(args[0], x)))
    raise _s.RegistryError(f"no reading of {name}")  # (read refuses it)


def _scalar(v: Any) -> bool:
    """an operand of && or ||: one element"""
    a = np.asarray(v)
    if a.size != 1:
        raise _Raised("MATLAB:nonLogicalConditional")
    return bool(a.reshape(-1)[0])


def _logical(v: Any) -> np.ndarray:
    if isinstance(v, _Arg):
        return _numbers(v) != 0
    return np.asarray(v, dtype=bool)


def _expand(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """MATLAB's implicit expansion: dimensions equal or one of them 1"""
    n = max(a.ndim, b.ndim)
    sa = a.shape + (1,) * (n - a.ndim)
    sb = b.shape + (1,) * (n - b.ndim)
    if any(p != q and p != 1 and q != 1 for p, q in zip(sa, sb, strict=True)):
        raise _Raised("MATLAB:sizeDimensionsMustMatch")
    return a.reshape(sa), b.reshape(sb)


def _reduce(name: str, a: np.ndarray) -> np.ndarray:
    """MATLAB's all and any: along the first dimension that is not 1 (0x0: one element)"""
    if a.shape == (0, 0):
        return np.array([[name == "all"]])
    dim = next((i for i, d in enumerate(a.shape) if d != 1), 0)
    f = np.all if name == "all" else np.any
    return f(a, axis=dim, keepdims=True)


def _model(x: Any) -> Any:
    """the value as the registry's model: MArray (numbers, logicals, char), Cell, MString,
    'struct', 'handle' or 'object'"""
    if isinstance(x, _v.MArray | _v.Cell | _v.MString):
        return x
    if _s._is_struct(x):
        return "struct"
    if callable(x) and not isinstance(x, type):
        return "handle"
    if type(x).__module__ == "numpy" and getattr(x, "dtype", None) is not None:
        if x.dtype.kind == "U" and x.ndim >= 2:  # a char matrix
            return _v.MArray("char", x.shape, np.asarray(x).reshape(-1, order="F").tolist())
    try:
        return _v.from_python(x)
    except (TypeError, ValueError):
        return "object"


def _numbers(v: Any) -> np.ndarray:
    """numbers of a computed value, or of the value checked: MATLAB's values of numbers,
    logicals and char, of MATLAB's size; MATLAB's errors for the other classes"""
    if not isinstance(v, _Arg):
        return np.asarray(v, dtype=np.float64)
    m = _model(v.value)
    if isinstance(m, _v.MString):
        raise _Raised("MATLAB:string:ComparisonNotDefined")
    if not isinstance(m, _v.MArray):
        raise _Raised("MATLAB:UndefinedFunction")
    vals = [float(ord(c)) for c in m.data] if m.mclass == "char" else [float(e) for e in m.data]
    return np.array(vals, dtype=np.float64).reshape(m.shape, order="F")


def _char_rows(m: _v.MArray) -> list[str] | None:
    """the texts of a 2-D char array's rows (None for more dimensions)"""
    if len(m.shape) != 2:
        return None
    return ["".join(r) for r in m.rows()]


def _same(a: str, b: str, fold: bool) -> bool:
    return a.lower() == b.lower() if fold else a == b


def _element_text(e: Any) -> str | None:
    """a cell element's text when it is a char row (a string or a char matrix is none)"""
    if isinstance(e, _v.MString):
        return None
    m = _model(e)
    if isinstance(m, _v.MArray) and m.mclass == "char" and len(m.shape) == 2 and m.shape[0] <= 1:
        return m.text()
    return None


def _strcmp(x: Any, words: tuple[str, ...], fold: bool) -> np.ndarray:
    n = len(words)
    m = _model(x)
    if isinstance(m, _v.MString):
        return np.array([[_same(str(m), w, fold) for w in words]])
    if isinstance(m, _v.MArray) and m.mclass == "char":
        rows = _char_rows(m)
        if rows is None:
            return np.zeros((1, n), dtype=bool)
        if len(rows) <= 1:
            t = rows[0] if rows else ""
            return np.array([[_same(t, w, fold) for w in words]])
        if len(rows) != n:
            raise _Raised("MATLAB:strcmp:RowsElementsMismatch")
        rows = [r.rstrip(" \t\n\v\f\r") for r in rows]
        return np.array([[_same(r, w, fold) for r, w in zip(rows, words, strict=True)]])
    if isinstance(m, _v.Cell):
        if m.shape == (1, 1):
            t = _element_text(m[0])
            return np.array([[t is not None and _same(t, w, fold) for w in words]])
        if m.shape != (1, n):
            raise _Raised("MATLAB:strcmp:InputsSizeMismatch")
        texts = [_element_text(e) for e in m]
        return np.array(
            [[t is not None and _same(t, w, fold) for t, w in zip(texts, words, strict=True)]]
        )
    return np.zeros((1, n), dtype=bool)


def _is_char(e: Any) -> bool:
    m = _model(e)
    return isinstance(m, _v.MArray) and m.mclass == "char"


def _ismember(x: Any, words: tuple[str, ...]) -> np.ndarray:
    m = _model(x)
    if isinstance(m, _v.MString):
        return np.array([[str(m) in words]])
    if isinstance(m, _v.MArray) and m.mclass == "char":
        rows = _char_rows(m)
        if rows is None:
            raise _Raised("MATLAB:cellstr:InputShape")
        if len(rows) <= 1:
            return np.array([[(rows[0] if rows else "") in words]])
        return np.array([[r.rstrip(" \t\n\v\f\r") in words] for r in rows])
    if isinstance(m, _v.Cell):
        if not all(_is_char(e) for e in m):
            raise _Raised("MATLAB:ISMEMBER:InputClass")
        texts = [_element_text(e) for e in m]
        flags = [t is not None and t in words for t in texts]
        return np.array(flags, dtype=bool).reshape(m.shape, order="F")
    if isinstance(m, _v.MArray):
        raise _Raised("MATLAB:ISMEMBER:InputClass")
    raise _Raised("MATLAB:ISMEMBER:InvalidInputsDataType")


def _lower(x: Any) -> Any:
    m = _model(x)
    if isinstance(m, _v.MString):
        return _v.MString(str(m).lower())
    if isinstance(m, _v.MArray) and m.mclass == "char":
        return _v.MArray("char", m.shape, [_lower_char(c) for c in m.data])
    if isinstance(m, _v.Cell):
        if not all(_is_char(e) for e in m):
            raise _Raised("MATLAB:lower:CellsMustContainChars")
        return _v.Cell([_lower(e) for e in m], m.shape)
    return x


def _lower_char(c: str) -> str:
    low = c.lower()
    return low if len(low) == 1 else c
