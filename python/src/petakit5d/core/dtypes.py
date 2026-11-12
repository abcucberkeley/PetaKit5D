"""MATLAB's class rules, exactly: conversions, rounding, integer-class and single/double arithmetic.

NumPy does none of these by default: ``np.round`` rounds half to even, ``astype`` truncates and
wraps, a float64 scalar or array promotes a float32 array, and a Python float is converted to
float32 *before* the operation. Port MATLAB code through these functions (see
docs/conventions.md for the measurements behind each rule):

- ``cast(x, cls)``: MATLAB's ``cast(x, cls)``, ``feval(cls, x)``, assignment into an array of
  class ``cls``, and the class function on a scalar. Rounds half away from zero (exactly),
  saturates, NaN gives 0.
- ``to_class(x, cls)`` (and ``uint16(x)``, ``uint8(x)``, ... below): the class function called
  directly, as ``uint16(x)`` in MATLAB code. Equal to ``cast`` except one case: for single
  arrays of 2 or more elements to int8, uint8, int16 or uint16, MATLAB rounds as
  ``floor(|x| + 0.5)`` in single precision, so ``0.49999997`` (the single below 0.5) gives 1
  (``-1`` for the signed classes) where ``cast`` gives 0. A scalar rounds exactly, unless it is
  the argument of a function file whose size MATLAB doesn't know there (``scalar_fast=True``).
- ``plus``, ``minus``, ``times``, ``rdivide``, ``ldivide``, ``uminus``: MATLAB's ``+ - .* ./ .\\``.
  Single with double (a scalar or an array) is computed in double and rounded to single;
  integers combine only with the same class or a double scalar (else ``TypeError``, MATLAB's
  ``MATLAB:mixedClasses``, also for logical and single operands); 8 to 32-bit results are
  computed in double, then rounded half away from zero and saturated (so
  ``uint16(1) + 0.49999999999999994`` is 2, as in MATLAB); 64-bit results are exact (rational)
  and then rounded and saturated, as MATLAB does.
- NaN: ``NAN`` is MATLAB's NaN (sign bit set); a NaN operand keeps its payload (the first NaN
  operand, quieted), as MATLAB on x86.
- ``mround``: MATLAB's ``round`` (half away from zero), class kept.
- ``EPS`` and ``eps(x)``: MATLAB's ``eps`` is the double ``2^-52``, also with single data.
"""

from __future__ import annotations

import operator
from fractions import Fraction

import numpy as np

from petakit5d.core.arrays import CLASS_DTYPES, matlab_class, numpy_dtype

__all__ = [
    "EPS",
    "NAN",
    "INT_CLASSES",
    "cast",
    "double",
    "eps",
    "idivide",
    "int8",
    "int16",
    "int32",
    "int64",
    "intmax",
    "intmin",
    "ldivide",
    "logical",
    "minus",
    "mround",
    "plus",
    "rdivide",
    "result_class",
    "single",
    "times",
    "to_class",
    "uint8",
    "uint16",
    "uint32",
    "uint64",
    "uminus",
]

# eps, eps('double'): MATLAB code's bare eps is this, also next to single data
EPS = 2.220446049250313e-16

# MATLAB's NaN is the x86 default NaN, with the sign bit set (0xFFF8000000000000); Python's
# float('nan') and np.nan are 0x7FF8000000000000. Use NAN where MATLAB code writes NaN.
NAN: float = float(np.array([0xFFF8000000000000], dtype=np.uint64).view(np.float64)[0])

INT_CLASSES = ("int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64")
_SMALL_INTS = ("int8", "uint8", "int16", "uint16", "int32", "uint32")
_FAST_PATH = ("int8", "uint8", "int16", "uint16")
_PRED_HALF_SINGLE = np.nextafter(np.float32(0.5), np.float32(0))


def intmax(cls: str) -> int:
    return int(np.iinfo(numpy_dtype(cls)).max)


def intmin(cls: str) -> int:
    return int(np.iinfo(numpy_dtype(cls)).min)


def _cls(x) -> str:
    if isinstance(x, str):
        return "char"
    return matlab_class(x)


def _arr(x) -> np.ndarray:
    if isinstance(x, str):
        return (
            np.array([ord(c) for c in x], dtype=np.float64).reshape(1, -1)
            if len(x) != 1
            else np.float64(ord(x))
        )
    if isinstance(x, bool):
        return np.asarray(x)
    if isinstance(x, int | float) and not isinstance(x, np.generic):
        return np.asarray(float(x))
    return np.asarray(x)


# -- rounding ----------------------------------------------------------------------------------


def mround(x):
    """MATLAB's ``round``: half away from zero, exactly, in the class of ``x``."""
    a = _arr(x)
    if a.dtype.kind != "f":
        return a.copy()
    with np.errstate(invalid="ignore"):
        t = np.trunc(a)
        frac = a - t  # exact (NaN for Inf, which stays Inf below)
        up = np.abs(frac) >= a.dtype.type(0.5)
        # (t + 0 would turn round(-0.3), which is -0, into +0)
        return np.where(up, t + np.sign(a), t).astype(a.dtype)


def eps(x=None):
    """MATLAB's ``eps``: ``eps()`` and ``eps('double')`` are 2^-52, ``eps('single')`` 2^-23,
    ``eps(x)`` the distance from ``|x|`` to the next larger number of x's class (NaN for
    Inf and NaN)."""
    if x is None or (isinstance(x, str) and x == "double"):
        return EPS
    if isinstance(x, str) and x == "single":
        return np.float32(2.0**-23)
    a = _arr(x)
    if a.dtype.kind != "f":
        raise TypeError("eps takes 'double', 'single' or a floating-point value")
    mant_bits, min_exp = (24, -149) if a.dtype == np.float32 else (53, -1074)
    m, e = np.frexp(np.abs(a))
    with np.errstate(invalid="ignore"):
        r = np.ldexp(np.ones_like(a), np.maximum(e - mant_bits, min_exp).astype(np.int32))
    r = np.where(np.isfinite(a), r, np.nan).astype(a.dtype)
    r = np.where(a == 0, np.ldexp(a.dtype.type(1), min_exp), r).astype(a.dtype)
    return r[()] if r.ndim == 0 else r


# -- conversions -------------------------------------------------------------------------------


def _float_to_int(r: np.ndarray, cls: str) -> np.ndarray:
    """Rounded float values (integral, or NaN/Inf) to an integer class with saturation."""
    dt = CLASS_DTYPES[cls]
    info = np.iinfo(dt)
    r = r.astype(np.float64)
    r = np.where(np.isnan(r), 0.0, r)
    if cls in _SMALL_INTS:
        return np.clip(r, info.min, info.max).astype(dt)
    out = np.empty(r.shape, dtype=dt)
    hi = r >= float(2**63 if cls == "int64" else 2**64)
    lo = r <= (float(-(2**63)) if cls == "int64" else 0.0)
    mid = ~(hi | lo)
    out[hi] = info.max
    out[lo] = info.min
    out[mid] = r[mid].astype(dt)
    return out


def _int_to_int(a: np.ndarray, cls: str) -> np.ndarray:
    dt = CLASS_DTYPES[cls]
    info = np.iinfo(dt)
    if a.dtype == np.uint64:
        if dt == np.uint64:
            return a.copy()
        return np.minimum(a, np.uint64(info.max)).astype(dt)
    w = a.astype(np.int64)
    if dt == np.uint64:
        return np.maximum(w, 0).astype(dt)
    return np.clip(w, info.min, info.max).astype(dt)


def _convert(x, cls: str, fast_path: bool, scalar_fast: bool = False):
    a = _arr(x)
    cls = matlab_class(cls) if not isinstance(cls, str) else cls
    if cls not in CLASS_DTYPES:
        raise ValueError(f"unknown class {cls!r}")
    dt = CLASS_DTYPES[cls]
    kind = a.dtype.kind
    if cls == "logical":
        if kind == "f" and np.any(np.isnan(a)):
            raise ValueError("NaN's cannot be converted to logicals. (MATLAB:nologicalnan)")
        return (a != 0) if kind != "b" else a.copy()
    if cls in ("double", "single"):
        with np.errstate(over="ignore", invalid="ignore"):
            return a.astype(dt)
    # integer target
    if kind == "b":
        return a.astype(dt)
    if kind in "iu":
        return _int_to_int(a, cls)
    if kind != "f":
        raise TypeError(f"cannot convert {a.dtype} to {cls}")
    if fast_path and a.dtype == np.float32 and (a.size >= 2 or scalar_fast) and cls in _FAST_PATH:
        # MATLAB's direct call on a single array: floor(|x| + 0.5) in single precision
        with np.errstate(invalid="ignore"):
            r = np.floor(np.abs(a) + np.float32(0.5)) * np.sign(a)
    else:
        r = mround(a)
    return _float_to_int(np.asarray(r), cls)


def cast(x, cls):
    """MATLAB's ``cast(x, cls)`` (also ``feval(cls, x)``, assignment into a ``cls`` array and the
    class function on a scalar): round half away from zero, saturate, NaN to 0."""
    return _convert(x, cls, fast_path=False)


def to_class(x, cls, *, scalar_fast: bool = False):
    """MATLAB's class function called directly, ``cls(x)`` (``uint16(x)`` in MATLAB code).

    Equals ``cast`` except for single arrays of 2 or more elements to int8, uint8, int16 and
    uint16, where MATLAB rounds as ``floor(|x| + 0.5)`` in single precision (``0.49999997``
    gives 1). Scalars round exactly; ``scalar_fast=True`` gives the array rounding to a scalar
    too, as MATLAB does for a scalar that is the argument of a function file.
    """
    return _convert(x, cls, fast_path=True, scalar_fast=scalar_fast)


def int8(x):
    return to_class(x, "int8")


def uint8(x):
    return to_class(x, "uint8")


def int16(x):
    return to_class(x, "int16")


def uint16(x):
    return to_class(x, "uint16")


def int32(x):
    return to_class(x, "int32")


def uint32(x):
    return to_class(x, "uint32")


def int64(x):
    return to_class(x, "int64")


def uint64(x):
    return to_class(x, "uint64")


def single(x):
    """``single(x)``: IEEE round to nearest (Inf beyond the single range)."""
    return to_class(x, "single")


def double(x):
    return to_class(x, "double")


def logical(x):
    """``logical(x)``: nonzero is true; NaN raises, as in MATLAB."""
    return to_class(x, "logical")


# -- arithmetic --------------------------------------------------------------------------------


_MIXED = (
    "Integers can only be combined with integers of the same class, or scalar doubles"
    " ({}, {}) (MATLAB:mixedClasses)"
)


def result_class(a, b) -> str:
    """The class of ``a op b`` for ``+ - .* ./ .\\`` (raises ``TypeError`` where MATLAB errors)."""
    ca, cb = _cls(a), _cls(b)
    ia, ib = ca in INT_CLASSES, cb in INT_CLASSES
    if ia or ib:
        if ia and ib:
            if ca != cb:
                raise TypeError(_MIXED.format(ca, cb))
            return ca
        other, oc = (b, cb) if ia else (a, ca)
        if oc in ("double", "char") and np.size(_arr(other)) == 1:
            return ca if ia else cb
        raise TypeError(_MIXED.format(ca, cb))
    if "single" in (ca, cb):
        return "single"
    return "double"


def _broadcast_pair(a: np.ndarray, b: np.ndarray):
    """MATLAB's implicit expansion: missing dimensions are trailing (a 1-D array is a row)."""
    if a.ndim == 0 or b.ndim == 0 or a.ndim == b.ndim:
        return a, b
    nd = max(a.ndim, b.ndim)

    def m(x):
        if x.ndim == 1:
            x = x.reshape(1, -1)
        return x.reshape(x.shape + (1,) * (nd - x.ndim))

    return m(a), m(b)


def _exact_int_op(op, a: np.ndarray, b: np.ndarray, cls: str) -> np.ndarray:
    """64-bit integer arithmetic, exactly (rational), rounded half away, saturated."""
    info = np.iinfo(CLASS_DTYPES[cls])
    a, b = np.broadcast_arrays(a, b)
    out = np.empty(a.shape, dtype=CLASS_DTYPES[cls])
    fa, fb = a.ravel(), b.ravel()
    res = out.reshape(-1)
    for k in range(fa.size):
        x, y = fa[k], fb[k]
        if (isinstance(x, np.floating) and not np.isfinite(x)) or (
            isinstance(y, np.floating) and not np.isfinite(y)
        ):
            with np.errstate(all="ignore"):
                r = op(np.float64(x), np.float64(y))
            res[k] = _float_to_int(np.asarray([r]), cls)[0]
            continue
        fx = Fraction(int(x)) if isinstance(x, np.integer) else Fraction(float(x))
        fy = Fraction(int(y)) if isinstance(y, np.integer) else Fraction(float(y))
        if op is operator.truediv and fy == 0:
            r = 0 if fx == 0 else (info.max if fx > 0 else info.min)
        else:
            q = op(fx, fy)
            n = abs(q.numerator) * 2 + q.denominator  # round half away: floor(|q| + 1/2)
            r = (n // (2 * q.denominator)) * (1 if q >= 0 else -1)
            if cls == "int64" and _int64_quirk(op, x, y):
                r += 1
        res[k] = min(max(r, info.min), info.max)
    return out


_PRED_HALF_D = float(np.nextafter(0.5, 0.0))


def _int64_quirk(op, x, y) -> bool:
    """MATLAB R2024b's int64 bug: a positive integer plus -0.49999999999999994 (as ``a + d``,
    ``d + a``, ``a - d`` or ``d - a``) gives one more than the exact rounding; nothing else
    differs in the measured grids (uint64, times and rdivide are exact)."""
    if op is operator.add:
        i, d = (x, y) if isinstance(x, np.integer) else (y, x)
        return isinstance(i, np.integer) and int(i) > 0 and float(d) == -_PRED_HALF_D
    if op is operator.sub:
        if isinstance(x, np.integer) and not isinstance(y, np.integer):
            return int(x) > 0 and float(y) == _PRED_HALF_D
        if isinstance(y, np.integer) and not isinstance(x, np.integer):
            return int(y) < 0 and float(x) == -_PRED_HALF_D
    return False


def _binop(op, a, b):
    cls = result_class(a, b)
    A, B = _broadcast_pair(_arr(a), _arr(b))
    if A.dtype.kind == "b":
        A = A.astype(np.float64)
    if B.dtype.kind == "b":
        B = B.astype(np.float64)
    if cls in ("int64", "uint64"):
        r = _exact_int_op(op, A, B, cls)
        return r[()] if r.ndim == 0 else r
    with np.errstate(invalid="ignore"):
        A64, B64 = A.astype(np.float64), B.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        r = op(A64, B64)
    if cls in ("double", "single"):
        r = _nan_operands(r, A64, B64)
    if cls == "double":
        return r
    if cls == "single":
        with np.errstate(over="ignore", invalid="ignore"):
            return r.astype(np.float32)
    return _float_to_int(np.asarray(mround(r)), cls)


_QUIET = np.uint64(0x0008000000000000)


def _nan_operands(r, a: np.ndarray, b: np.ndarray):
    """x86's rule, which MATLAB shows: a NaN result of a NaN operand is the first NaN operand,
    quieted (NumPy may swap the operands of + and *); other NaNs are the default NaN."""
    na, nb = np.isnan(a), np.isnan(b)
    if not (na.any() or nb.any()):
        return r
    r = np.array(r, dtype=np.float64, copy=True)
    a, b, r = np.broadcast_arrays(a, b, r)
    r = r.copy()
    na, nb = (
        np.broadcast_to(na, r.shape),
        np.broadcast_to(nb, r.shape) & ~np.broadcast_to(na, r.shape),
    )
    r[na] = (a[na].view(np.uint64) | _QUIET).view(np.float64)
    r[nb] = (b[nb].view(np.uint64) | _QUIET).view(np.float64)
    return r[()] if r.ndim == 0 else r


def plus(a, b):
    """MATLAB's ``a + b``."""
    return _binop(operator.add, a, b)


def minus(a, b):
    """MATLAB's ``a - b``."""
    return _binop(operator.sub, a, b)


def times(a, b):
    """MATLAB's ``a .* b`` (and ``a * b`` with a scalar)."""
    return _binop(operator.mul, a, b)


def rdivide(a, b):
    """MATLAB's ``a ./ b`` (and ``a / b`` with a scalar ``b``)."""
    return _binop(operator.truediv, a, b)


def ldivide(a, b):
    """MATLAB's ``a .\\ b`` (``b ./ a``)."""
    return _binop(operator.truediv, b, a)


def uminus(a):
    """MATLAB's ``-a`` (saturating for integers: ``-int8(-128)`` is 127)."""
    x = _arr(a)
    cls = _cls(a)
    if cls in INT_CLASSES:
        if cls.startswith("u"):
            return np.zeros_like(x)
        info = np.iinfo(x.dtype)
        return np.where(x == info.min, info.max, -np.where(x == info.min, 0, x)).astype(x.dtype)
    if x.dtype.kind == "b":
        return -x.astype(np.float64)
    return -x


def idivide(a, b, opt: str = "fix"):
    """MATLAB's ``idivide(a, b, opt)``: integer division rounded by 'fix' (default), 'round',
    'floor' or 'ceil'; at least one argument is an integer, both of the same class or one a
    double scalar."""
    cls = result_class(a, b)
    if cls not in INT_CLASSES:
        raise TypeError(
            "At least one argument must belong to an integer class. (MATLAB:idivide:mustBeInteger)"
        )
    if opt == "round":
        return rdivide(a, b)
    A, B = _broadcast_pair(_arr(a), _arr(b))
    if cls in ("int64", "uint64"):
        info = np.iinfo(CLASS_DTYPES[cls])
        A, B = np.broadcast_arrays(A, B)
        out = np.empty(A.shape, dtype=CLASS_DTYPES[cls])
        for idx in np.ndindex(A.shape):
            x, y = (
                Fraction(int(A[idx])),
                Fraction(int(mround(B[idx])) if B.dtype.kind == "f" else int(B[idx])),
            )
            if y == 0:
                r = 0 if x == 0 else (info.max if x > 0 else info.min)
            else:
                q = x / y
                r = {"fix": int(q), "floor": q.__floor__(), "ceil": q.__ceil__()}[opt]
            out[idx] = min(max(r, info.min), info.max)
        return out[()] if out.ndim == 0 else out
    # 8 to 32 bits: the double quotient is exact enough that floor/ceil/fix are exact
    if B.dtype.kind == "f":
        B = mround(B)
    if A.dtype.kind == "f":
        A = mround(A)
    with np.errstate(divide="ignore", invalid="ignore"):
        q = A.astype(np.float64) / B.astype(np.float64)
    f = {"fix": np.trunc, "floor": np.floor, "ceil": np.ceil}.get(opt)
    if f is None:
        raise ValueError(f"unknown rounding option {opt!r}")
    return _float_to_int(f(q), cls)
