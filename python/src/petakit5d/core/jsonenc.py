"""``jsonencode`` as MATLAB R2024b writes it (compact and ``PrettyPrint``).

Python values stand for MATLAB values as follows:

| Python | MATLAB | JSON |
|---|---|---|
| ``dict`` | 1x1 struct (field order kept) | object (``{}`` without fields) |
| ``list``, ``tuple`` | 1xN cell | array, each element encoded alone |
| ``str`` | char row vector | string |
| ``bool``, ``np.bool_`` | logical scalar | ``true``/``false`` |
| ``int``, ``float`` | double scalar | number |
| NumPy scalar | scalar of its class | number |
| ``np.ndarray`` | array of its class (1-D: a row) | see below |
| ``None`` | ``[]`` | ``[]`` |

An array is ``[]`` when empty, its value with 1 element, a flat array for a vector (either
orientation), else nested rows (``a.tolist()`` order of the MATLAB-sized array).

Numbers: double values take ``%.15g`` if it reads back as the same double, else ``%.17g``
(single: ``%.6g``, else ``%.9g``, read back in single). Decimal exponents -4 to 5 print plainly,
others as ``dE-n``, ``d.0E+n`` or ``d.ddE±n``. Integers print in full, NaN and Inf as ``null``
(``NaN``, ``Infinity``, ``-Infinity`` with ``convert_inf_and_nan=False``). Strings escape ``"``,
``\\`` and control characters (``\\b \\t \\n \\f \\r``, others ``\\u00xx``); ``/`` and non-ASCII
characters are written as they are.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = ["format_json_number", "jsonencode"]

_SHORT = {
    0x08: "\\b",
    0x09: "\\t",
    0x0A: "\\n",
    0x0C: "\\f",
    0x0D: "\\r",
    0x22: '\\"',
    0x5C: "\\\\",
}


def _string(s: str) -> str:
    out = ['"']
    for ch in s:
        c = ord(ch)
        if c in _SHORT:
            out.append(_SHORT[c])
        elif c < 0x20:
            out.append("\\u%04x" % c)  # noqa: UP031 (C formatting, as MATLAB's)
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def format_json_number(x, single: bool = False) -> str:
    """One finite floating-point number as ``jsonencode`` prints it."""
    x = float(x)
    p_short, p_long = (6, 9) if single else (15, 17)
    s = "%.*g" % (p_short, x)  # noqa: UP031 (C formatting, as MATLAB's)
    back = float(np.float32(float(s))) if single else float(s)
    p = p_short if back == x else p_long
    mant, ex = ("%.*e" % (p - 1, x)).split("e")  # noqa: UP031 (C formatting, as MATLAB's)
    ex = int(ex)
    if -5 < ex < 6:
        return "%.*g" % (p, x)  # noqa: UP031 (C formatting, as MATLAB's)
    sign = "-" if mant.startswith("-") else ""
    digits = mant.lstrip("-").replace(".", "").rstrip("0") or "0"
    if len(digits) > 1:
        m = digits[0] + "." + digits[1:]
    else:
        m = digits + (".0" if ex >= 0 else "")
    return f"{sign}{m}E{'+' if ex >= 0 else '-'}{abs(ex)}"


def _scalar(v, kind: str, single: bool, convert: bool) -> str:
    if kind == "b":
        return "true" if v else "false"
    if kind in "iu":
        return str(int(v))
    f = float(v)
    if not math.isfinite(f):
        if convert:
            return "null"
        return "NaN" if math.isnan(f) else ("Infinity" if f > 0 else "-Infinity")
    return format_json_number(f, single)


class _Writer:
    def __init__(self, pretty: bool, convert: bool):
        self.pretty = pretty
        self.convert = convert

    def seq(self, items: list[str], depth: int, open_: str, close: str) -> str:
        if not items:
            return open_ + close
        if not self.pretty:
            return open_ + ",".join(items) + close
        ind = "  " * (depth + 1)
        return open_ + "\n" + ",\n".join(ind + it for it in items) + "\n" + "  " * depth + close

    def value(self, v, depth: int) -> str:
        if v is None:
            return "[]"
        if isinstance(v, dict):
            items = []
            for k, x in v.items():
                sep = ": " if self.pretty else ":"
                items.append(_string(str(k)) + sep + self.value(x, depth + 1))
            return self.seq(items, depth, "{", "}")
        if isinstance(v, list | tuple):
            return self.seq([self.value(x, depth + 1) for x in v], depth, "[", "]")
        if isinstance(v, str):
            return _string(v)
        if isinstance(v, bool | np.bool_):
            return "true" if v else "false"
        if isinstance(v, int) and not isinstance(v, np.generic):
            return _scalar(float(v), "f", False, self.convert)
        if isinstance(v, float) and not isinstance(v, np.generic):
            return _scalar(v, "f", False, self.convert)
        a = np.asarray(v)
        if a.dtype.kind == "U":
            if a.ndim == 0:
                return _string(str(a))
            return self.seq(
                [self.value(str(x), depth + 1) for x in a.ravel(order="F")], depth, "[", "]"
            )
        if a.dtype.kind not in "biuf":
            raise TypeError(f"jsonencode cannot encode {type(v).__name__} of dtype {a.dtype}")
        single = a.dtype == np.float32
        kind = a.dtype.kind
        if a.size == 0:
            return "[]"
        if a.size == 1:
            return _scalar(a.reshape(-1)[0], kind, single, self.convert)
        from petakit5d.core.arrays import msize

        shape = msize(a)
        m = a.reshape(shape, order="F") if a.ndim != len(shape) or a.shape != shape else a
        if len(shape) == 2 and (shape[0] == 1 or shape[1] == 1):
            flat = m.ravel(order="F")
            return self.seq([_scalar(x, kind, single, self.convert) for x in flat], depth, "[", "]")
        return self._nested(m, depth, kind, single)

    def _nested(self, m: np.ndarray, depth: int, kind: str, single: bool) -> str:
        if m.ndim == 1:
            return self.seq([_scalar(x, kind, single, self.convert) for x in m], depth, "[", "]")
        return self.seq(
            [self._nested(m[i], depth + 1, kind, single) for i in range(m.shape[0])],
            depth,
            "[",
            "]",
        )


def jsonencode(value, *, pretty_print: bool = False, convert_inf_and_nan: bool = True) -> str:
    """MATLAB's ``jsonencode(value, 'PrettyPrint', pretty_print, 'ConvertInfAndNaN', ...)``."""
    return _Writer(pretty_print, convert_inf_and_nan).value(value, 0)
