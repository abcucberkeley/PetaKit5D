"""Numbers as MATLAB formats and parses them: ``sprintf``, ``num2str``, ``int2str``, ``mat2str``.

The MATLAB wrappers hand parameters to their tasks as text (``%.10f``, ``%.20d``, ``%d`` of a
non-integer, ``mat2str`` with 15 digits, ``num2str(x, '%.10f,')``), and folder names hold
``num2str`` and ``%g`` forms. ``through_text`` applies the same rounding at the same boundary.

MATLAB's rules that differ from C and Python (measured on R2024b, see docs/conventions.md):

- ``%d``, ``%i`` print integer values in ``[-2^63, 2^63)`` as integers; any other value (a
  fraction, a larger value, a uint64 above ``intmax('int64')``) is printed with ``%e`` and the
  same flags, width and precision (``sprintf('%.20d', 0.108)`` is ``1.07999999999999998779e-01``).
  ``%u``, ``%x``, ``%o`` take ``[0, 2^64)``; negative values switch to ``%e``.
- ``%s`` and ``%c`` of a number print the character for integers in ``[0, 65535]``, else ``%e``.
- ``Inf``, ``-Inf``, ``NaN`` print as such for every numeric conversion (width and ``-`` kept,
  ``0`` ignored, ``+`` kept for ``Inf`` only).
- The format is applied again while data remain; it stops at the first conversion without
  data. Without any data argument it stops at the first conversion; with only empty data it
  prints the format once with empty conversions.
- Escapes ``\\n \\t \\\\ \\" \\a \\b \\f \\r \\v \\xHH \\NNN`` are processed; an unknown
  escape or an invalid conversion ends the output there (MATLAB warns).
"""

from __future__ import annotations

import math
import re

import numpy as np

__all__ = [
    "int2str",
    "mat2str",
    "mat2str_comma",
    "matlab_strtrim",
    "num2str",
    "sprintf",
    "str2double",
    "through_mat2str",
    "through_text",
]

_I64_MIN = -(2**63)
_I64_LIM = 2**63
_U64_LIM = 2**64

# -- data items --------------------------------------------------------------------------------


class _Char:
    """A char array argument: consumed whole by %s, element-wise by the other conversions."""

    __slots__ = ("text",)

    def __init__(self, text: str):
        self.text = text


def _items(args) -> list:
    items = []
    for a in args:
        if isinstance(a, str):
            items.append(_Char(a))
        elif isinstance(a, dict | set) or (
            isinstance(a, list | tuple)
            and any(not isinstance(v, int | float | bool | np.number | np.bool_) for v in a)
        ):
            # a cell array in MATLAB (a list of numbers is a numeric vector here)
            raise TypeError(
                "Only numeric, logical and text arguments (MATLAB:printf:invalidInputType)"
            )
        else:
            arr = np.asarray(a)
            if arr.dtype.kind == "U":
                items.append(_Char("".join(arr.ravel(order="F").tolist())))
                continue
            if arr.dtype.kind not in "biuf":
                raise TypeError(
                    "Only numeric, logical and text arguments (MATLAB:printf:invalidInputType)"
                )
            if arr.size:
                items.append(arr.ravel(order="F"))
    return items


def _num(v):
    """A data element as int (exact integer classes, logical) or float."""
    if isinstance(v, np.bool_ | bool):
        return int(v)
    if isinstance(v, np.integer):
        return int(v)
    return float(v)


# -- one conversion ----------------------------------------------------------------------------

_SPEC = re.compile(r"%([-+ 0#]*)(\*|\d+)?(?:\.(\*|\d*))?([diouxXfFeEgGcs%]?)")


def _pad(s: str, flags: str, width: int | None, zero_ok: bool = False) -> str:
    if width is None or len(s) >= width:
        return s
    if "-" in flags:
        return s + " " * (width - len(s))
    if zero_ok and "0" in flags:
        sign = s[0] if s[:1] in "+- " else ""
        return sign + "0" * (width - len(s)) + s[len(sign) :]
    return " " * (width - len(s)) + s


def _nonfinite(v: float, flags: str, width: int | None) -> str:
    if math.isnan(v):
        s = "NaN"
    else:
        s = "Inf" if v > 0 else "-Inf"
        if v > 0 and "+" in flags:
            s = "+Inf"
        elif v > 0 and " " in flags:
            s = " Inf"
    return _pad(s, flags.replace("0", ""), width)


def _pyfmt(flags: str, width: int | None, prec: int | None, conv: str, v) -> str:
    spec = (
        "%"
        + flags
        + (str(width) if width is not None else "")
        + (f".{prec}" if prec is not None else "")
        + conv
    )
    return spec % v


def _format_float(flags, width, prec, conv, v: float) -> str:
    if not math.isfinite(v):
        return _nonfinite(v, flags, width)
    return _pyfmt(flags, width, prec, conv, v)


def _format_int(flags, width, prec, conv, iv: int) -> str:
    if conv == "o" and "#" in flags:
        # C prints a leading 0, Python '0o'
        s = _pyfmt(flags.replace("#", ""), None, prec, "o", iv)
        if not s.startswith("0"):
            s = "0" + s
        return _pad(s, flags, width, zero_ok=prec is None)
    return _pyfmt(flags, width, prec, conv, iv)


def _convert(flags: str, width: int | None, prec: int | None, conv: str, v) -> str:
    """One numeric value through a conversion, with MATLAB's switch to %e."""
    if conv in "cs":
        if isinstance(v, int) or (math.isfinite(v) and v == math.floor(v)):
            iv = int(v)
            if 0 <= iv <= 65535:
                s = chr(iv)
                if prec is not None and conv == "s":
                    s = s[:prec]
                return _pad(s, flags, width)
        if isinstance(v, float) and not math.isfinite(v):
            return _nonfinite(v, flags, width)
        return _format_float(flags, width, prec, "e", float(v))
    if conv in "di":
        if isinstance(v, float):
            if not math.isfinite(v):
                return _nonfinite(v, flags, width)
            if v != math.floor(v) or not (_I64_MIN <= v < _I64_LIM):
                return _format_float(flags, width, prec, "e", v)
            v = int(v)
        if not (_I64_MIN <= v < _I64_LIM):
            return _format_float(flags, width, prec, "e", float(v))
        return _format_int(flags, width, prec, "d", v)
    if conv in "uoxX":
        if isinstance(v, float):
            if not math.isfinite(v):
                return _nonfinite(v, flags, width)
            if v != math.floor(v) or not (0 <= v < _U64_LIM):
                return _format_float(flags, width, prec, "e", v)
            v = int(v)
        if not (0 <= v < _U64_LIM):
            return _format_float(flags, width, prec, "e", float(v))
        return _format_int(flags, width, prec, "d" if conv == "u" else conv, v)
    return _format_float(flags, width, prec, conv, float(v))


# -- escapes -----------------------------------------------------------------------------------

_SIMPLE_ESC = {
    "n": "\n",
    "t": "\t",
    "\\": "\\",
    '"': '"',
    "a": "\a",
    "b": "\b",
    "f": "\f",
    "r": "\r",
    "v": "\v",
}


def _unescape(fmt: str) -> tuple[list, bool]:
    """Split a format into literal text and conversion specs, processing escapes.

    Returns ``(parts, complete)``: parts are ``str`` or match objects of ``_SPEC``; ``complete``
    is False when an unknown escape or an invalid conversion ended the format early.
    """
    parts: list = []
    buf: list[str] = []
    i, n = 0, len(fmt)
    while i < n:
        c = fmt[i]
        if c == "\\":
            if i + 1 >= n:
                buf.append("\\")
                i += 1
                continue
            e = fmt[i + 1]
            if e in _SIMPLE_ESC:
                buf.append(_SIMPLE_ESC[e])
                i += 2
            elif e == "x":
                m = re.match(r"[0-9a-fA-F]+", fmt[i + 2 :])
                if not m:
                    return parts + (["".join(buf)] if buf else []), False
                buf.append(chr(int(m.group(0), 16)))
                i += 2 + len(m.group(0))
            elif e in "01234567":
                m = re.match(r"[0-7]{1,3}", fmt[i + 1 :])
                buf.append(chr(int(m.group(0), 8)))
                i += 1 + len(m.group(0))
            else:
                return parts + (["".join(buf)] if buf else []), False
        elif c == "%":
            m = _SPEC.match(fmt, i)
            if fmt.startswith("%%", i):
                buf.append("%")
                i += 2
                continue
            if m is None or m.group(4) in ("", "%"):
                return parts + (["".join(buf)] if buf else []), False
            if buf:
                parts.append("".join(buf))
                buf = []
            parts.append(m)
            i = m.end()
        else:
            buf.append(c)
            i += 1
    if buf:
        parts.append("".join(buf))
    return parts, True


# -- sprintf -----------------------------------------------------------------------------------


class _Data:
    def __init__(self, items: list):
        self.items = items
        self.i = 0  # item
        self.j = 0  # element inside the item

    def empty(self) -> bool:
        return self.i >= len(self.items)

    def _advance(self, n: int, size: int):
        self.j += n
        if self.j >= size:
            self.i += 1
            self.j = 0

    def take_number(self):
        it = self.items[self.i]
        if isinstance(it, _Char):
            v = ord(it.text[self.j])
            self._advance(1, len(it.text))
            return v
        v = _num(it[self.j])
        self._advance(1, it.size)
        return v

    def take_text(self):
        """For %s: a whole char item (from the current position), or one number."""
        it = self.items[self.i]
        if isinstance(it, _Char):
            s = it.text[self.j :]
            self.i += 1
            self.j = 0
            return s
        return self.take_number()

    def take_char(self):
        it = self.items[self.i]
        if isinstance(it, _Char):
            s = it.text[self.j]
            self._advance(1, len(it.text))
            return s
        return self.take_number()


def sprintf(fmt: str, *args) -> str:
    """MATLAB's ``sprintf`` for numeric, logical and char arguments (arrays column-major)."""
    parts, complete = _unescape(fmt)
    has_args = len(args) > 0
    data = _Data(_items(args))
    convs = [p for p in parts if not isinstance(p, str)]
    out: list[str] = []

    def spec_values(m, need_data: bool):
        flags, width, prec = m.group(1), m.group(2), m.group(3)
        w = p = None
        if width == "*":
            w = int(data.take_number()) if need_data else None
        elif width:
            w = int(width)
        if prec == "*":
            p = int(data.take_number()) if need_data else None
        elif prec is not None:
            p = int(prec) if prec else 0
        if w is not None and w < 0:
            flags, w = flags + "-", -w
        return flags, w, p

    if not has_args or (data.empty() and not convs):
        # no data: up to the first conversion (the whole format if it has none)
        for p in parts:
            if not isinstance(p, str):
                break
            out.append(p)
        return "".join(out)
    if data.empty():
        # only empty data: the format once, conversions print nothing
        out.extend(p for p in parts if isinstance(p, str))
        return "".join(out)

    while True:
        for p in parts:
            if isinstance(p, str):
                out.append(p)
                continue
            if data.empty():
                return "".join(out)
            flags, w, pr = spec_values(p, True)
            if data.empty():
                return "".join(out)
            conv = p.group(4)
            if conv == "s":
                v = data.take_text()
                if isinstance(v, str):
                    out.append(_pad(v if pr is None else v[:pr], flags, w))
                else:
                    out.append(_convert(flags, w, pr, "s", v))
            elif conv == "c":
                v = data.take_char()
                out.append(
                    _pad(v, flags, w) if isinstance(v, str) else _convert(flags, w, pr, "c", v)
                )
            else:
                out.append(_convert(flags, w, pr, conv, data.take_number()))
        if not complete or data.empty() or not convs:
            return "".join(out)


# -- trimming ----------------------------------------------------------------------------------

_WS = " \t\n\v\f\r\0"


def matlab_strtrim(s: str) -> str:
    """``strtrim``: leading and trailing whitespace and null characters removed."""
    return s.strip(_WS)


# -- int2str and num2str -----------------------------------------------------------------------


def _as_matrix(x) -> np.ndarray:
    a = np.asarray(x)
    if a.ndim == 0:
        a = a.reshape(1, 1)
    elif a.ndim == 1:
        a = a.reshape(1, -1)
    elif a.ndim > 2:
        a = a.reshape(a.shape[0], -1, order="F")
    return a


def _round_half_away(a: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        t = np.trunc(a)
        return np.where(np.abs(a - t) >= 0.5, t + np.sign(a), t)


def _width_max(a: np.ndarray) -> float:
    w = a.astype(np.float64) if a.dtype.kind == "f" else a
    if a.dtype.kind == "f":
        w = np.where(np.isfinite(w), w, 314.0)
    return float(np.max(np.abs(w.astype(np.float64)))) if w.size else 0.0


def _rows_text(rows: list[str]) -> str:
    if len(rows) != 1:
        raise NotImplementedError(
            "a matrix gives a char matrix in MATLAB; only scalars and row vectors here"
        )
    return rows[0]


def int2str(x) -> str:
    """``int2str``: rounded (half away from zero), printed with ``%.0f`` (``%d`` for integers)."""
    a = _as_matrix(x)
    if a.size == 0:
        return ""
    if a.dtype.kind == "f":
        a = _round_half_away(a.astype(np.float64)) + 0.0
        conv = ".0f"
    elif a.dtype.kind == "b":
        a = a.astype(np.uint8)
        conv = "d"
    else:
        conv = "u" if a.dtype == np.uint64 else "d"
    if a.size == 1:
        return sprintf("%" + conv, a)
    wmax = _width_max(np.asarray(x) if np.asarray(x).dtype.kind == "f" else a)
    width = 3 if wmax == 0 else int(math.floor(math.log10(wmax))) + 3
    f = f"%{width}{conv}"
    rows = [sprintf(f, a[r, :]) for r in range(a.shape[0])]
    if any(len(r) != width * a.shape[1] for r in rows):
        # MATLAB assigns each row into a char matrix of width*cols columns and errors
        # (for example when all values are below 0.01 and the width is negative)
        raise ValueError(
            "Unable to perform assignment because the size of the left side differs"
            " (MATLAB:subsassigndimmismatch)"
        )
    return matlab_strtrim(_rows_text(rows))


def _recycled(a: np.ndarray, d: int) -> str:
    f = f"%{d + 7}.{d}g"
    rows = []
    for r in range(a.shape[0]):
        s = sprintf(f, a[r, :])
        if a.shape[1] > 1 and np.min(a[r, :]) < 0:
            # a space before each '-' that follows a non-space other than e/E
            s = re.sub(r"([^\sEe])-", r"\1 -", s)
        rows.append(s)
    return matlab_strtrim(_rows_text(rows))


def num2str(x, f=None) -> str:
    """``num2str`` with no format, a precision (number) or a format (text).

    Scalars and row vectors (a matrix makes a char matrix in MATLAB and raises here).
    """
    if isinstance(x, str):
        return x
    a = _as_matrix(x)
    if a.size == 0:
        return ""
    if a.dtype.kind not in "biuf":
        raise TypeError("num2str takes numeric or logical input")
    if a.dtype.kind == "f":
        a = a + 0.0  # no negative zero
    is_int = a.dtype.kind in "iub"
    if f is None or (is_int and not isinstance(f, str)):
        xmax = _width_max(a)
        if is_int or (
            bool(np.all((a == np.trunc(a)) | np.isnan(a))) and np.spacing(np.float64(xmax)) <= 1
        ):
            return int2str(a)
        d = 1 if xmax == 0 else min(12, max(1, int(math.floor(math.log10(xmax))) + 1)) + 4
        return _recycled(a, d)
    if not isinstance(f, str):
        d = int(_round_half_away(np.asarray(float(f))))
        return _recycled(a, d)
    if "%" not in f:
        raise ValueError(f"Invalid format: {f}")
    # each element formatted alone; columns right-justified to the widest
    cols = [sprintf(f, a[0, j]) for j in range(a.shape[1])] if a.shape[0] == 1 else None
    if cols is None:
        raise NotImplementedError("num2str of a matrix with a format is not supported")
    left = "-" in f.split("%", 1)[1][:4]
    width = max(len(c) for c in cols)
    cols = [(c.ljust(width) if left else c.rjust(width)) if len(c) < width else c for c in cols]
    return matlab_strtrim("".join(cols))


# -- mat2str -----------------------------------------------------------------------------------


def _fmt_g(v, precision: int) -> str:
    if isinstance(v, int):
        return str(v)
    if not math.isfinite(v):
        return "NaN" if math.isnan(v) else ("Inf" if v > 0 else "-Inf")
    return "%.*g" % (precision, v)  # noqa: UP031 (C formatting, as MATLAB's)


def mat2str(a, precision: int | None = None, *, cls: bool = False) -> str:
    """``mat2str`` (15 significant digits by default) of a scalar, vector or matrix."""
    prec_given = precision is not None
    precision = 15 if precision is None else int(precision)
    if isinstance(a, str):
        if a == "":
            return "''"
        dangerous = re.compile(r"[\0\n-\r]")
        s = a.replace("'", "''")
        if dangerous.search(s):
            s = dangerous.sub(lambda m: f"' char({ord(m.group(0))}) '", s)
            return "['" + s + "']"
        return "'" + s + "'"
    arr = np.asarray(a)
    if arr.ndim > 2:
        raise ValueError("Input matrix must be 2-D.")
    m = _as_matrix(arr)
    rows, cols = m.shape
    cname = None
    if cls:
        from petakit5d.core.arrays import matlab_class

        cname = matlab_class(arr if arr.ndim else a)
    if m.size == 0:
        s = "[]" if rows == 0 and cols == 0 else f"zeros({rows},{cols})"
        return f"{cname}({s})" if cname else s
    if m.dtype.kind == "b":
        vals = [["true" if v else "false" for v in row] for row in m]
    else:
        if m.dtype.kind == "f":
            m = m + 0.0  # no negative zero
        exact_int = False
        if m.dtype in (np.dtype(np.int64), np.dtype(np.uint64)):
            # %d / %u unless a given precision has fewer digits than the largest value
            # (MATLAB counts the digits as the first k with |x| <= 10^k)
            big = int(np.max(np.abs(m.astype(object))))
            ndig = next(k for k in range(1, 21) if big <= 10**k or k == 20)
            exact_int = (not prec_given) or precision >= ndig
        if exact_int:
            vals = [[str(int(v)) for v in row] for row in m]
        else:
            # %.<precision>g of the value as a double (integers too)
            vals = [[_fmt_g(float(v), precision) for v in row] for row in m]
    s = ";".join(" ".join(r) for r in vals)
    if m.size != 1:
        s = "[" + s + "]"
    if cname and m.dtype.kind != "b":
        s = f"{cname}({s})"
    return s


def mat2str_comma(a, sn: int | None = None) -> str:
    """``mat2str_comma``: ``mat2str`` with commas, or ``sn`` decimals with ``%0.<sn>f,``."""
    if sn is None:
        s = mat2str(a).replace(" ", ",")
    else:
        s = "[" + num2str(a, f"%0.{int(sn)}f,").strip(",") + "]"
    return s.replace(" ", "")


# -- parsing -----------------------------------------------------------------------------------

_FLOAT = re.compile(
    r"\s*([+-]?(?:0[xX](?:[0-9a-fA-F]+\.?[0-9a-fA-F]*|\.[0-9a-fA-F]+)(?:[pP][+-]?\d+)?"
    r"|(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|[iI][nN][fF]|[nN][aA][nN]))"
)


def _nan() -> float:
    from petakit5d.core.dtypes import NAN

    return NAN


def str2double(s) -> float:
    """``str2double`` of a real number; MATLAB's NaN (``dtypes.NAN``) for anything else.

    As MATLAB: leading and trailing whitespace and commas are dropped (``'1,000'`` is 1000),
    hexadecimal reads (``'0x10'`` is 16), ``'Inf'`` and ``'NaN'`` in any case, but not
    ``'Infinity'``; complex text gives NaN here (MATLAB returns a complex number).
    """
    if not isinstance(s, str):
        return _nan()
    for text in (s, s.rstrip(_WS).replace(",", "")):
        m = _FLOAT.match(text)
        if m and m.end() == len(text):
            t = m.group(1)
            low = t.lower().lstrip("+-")
            if low == "inf":
                return -math.inf if t.startswith("-") else math.inf
            if low == "nan":
                return _nan()
            if low.startswith("0x"):
                v = float.fromhex(low if "p" in low or "." in low else low)
                return -v if t.startswith("-") else v
            return float(t)
    return _nan()


def through_text(value, fmt: str):
    """A value after MATLAB's trip through text: ``str2double(sprintf(fmt, v))`` per element.

    Use where a wrapper passes a number to its tasks in a formatted string (``%.10f``,
    ``%.20d``, ``%d``): ``through_text(dz, '%.10f')``. Arrays keep their shape (float64).
    """
    a = np.asarray(value)
    if a.ndim == 0:
        return str2double(sprintf(fmt, a))
    flat = [str2double(sprintf(fmt, v)) for v in a.ravel(order="F")]
    return np.asarray(flat, dtype=np.float64).reshape(a.shape, order="F")


def through_mat2str(value, precision: int = 15):
    """A value after ``mat2str`` and MATLAB's parser: 15 significant digits per element."""
    a = np.asarray(value)
    conv = lambda v: float(_fmt_g(float(v), precision).replace("Inf", "inf"))  # noqa: E731
    if a.ndim == 0:
        return conv(a)
    return np.asarray([conv(v) for v in a.ravel(order="F")], dtype=np.float64).reshape(
        a.shape, order="F"
    )
