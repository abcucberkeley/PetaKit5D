"""Threaded element-wise primitives and the fixed-order double sum (petakit5d._ext.core).

MATLAB runs element-wise expressions on whole volumes on all its threads; NumPy runs them on
one. These primitives give MATLAB's bits at MATLAB's speed:

- ``convert(x, cls)``: a class conversion. ``mode='class'`` is the class function called
  directly (``uint16(x)``), ``mode='cast'`` is ``cast(x, 'uint16')`` (see
  docs/conventions.md 3.1: they differ for single arrays to 8 and 16-bit integers).
- ``scale_offset(x, scale, offset, hi=, lo=, mask=, out_class=)``: MATLAB's
  ``x * scale + offset``, then ``min(., hi)``, ``max(., lo)``, ``.* mask`` and a conversion,
  each step only when given, each rounded to the class of ``x`` as MATLAB rounds it.
- ``clamp(x, lo, hi)``: ``max(min(x, hi), lo)``.
- ``compare(x, op, value)``: ``x == value`` (``'eq'``, ``'ne'``, ``'lt'``, ``'le'``, ``'gt'``,
  ``'ge'``) as a logical array.
- ``sum_double(x)``: ``sum(x(:), 'double')`` in MATLAB's order at one thread, at any thread
  count (decon_sum_double).

Arrays are the port's ``(z, y, x)`` C arrays of any shape (docs/conventions.md 1). The
element-wise primitives run in memory order, so a MATLAB-layout Fortran array that a kernel
keeps inside works too; results have the layout of ``x`` (a non-contiguous ``x`` is copied to C
order first). ``sum_double`` adds in MATLAB's linear order (y fastest), read in place through
a transposing scratch; ``sum_double(v, matlab_layout=True)`` sums MATLAB's own array
contiguously, for kernels that keep MATLAB's layout inside (faster: no transpose).
Each primitive has a NumPy reference built on ``petakit5d.core.dtypes`` (``convert_ref``,
``scale_offset_ref``, ``clamp_ref``, ``compare_ref``, ``sum_double_ref``), used when the
extension is not built, and in the tests. ``num_threads`` (None or 0: all CPUs available)
never changes a result.
"""

from __future__ import annotations

import functools
import math

import numpy as np

from petakit5d.core import dtypes
from petakit5d.core.arrays import CLASS_DTYPES, matlab_class, matlab_order
from petakit5d.core.threads import resolve_threads

__all__ = [
    "SUM_BLOCK",
    "SUM_BLOCK_FROM",
    "SUM_SEQUENTIAL_UPTO",
    "clamp",
    "clamp_ref",
    "compare",
    "compare_ref",
    "convert",
    "convert_ref",
    "extension_available",
    "scale_offset",
    "scale_offset_ref",
    "sum_double",
    "sum_double_ref",
]

# MATLAB's sum(x, 'double') at one thread: index order up to 43 elements, one block of four
# lanes up to 88,999, from 89,000 on blocks of 4,096 (each in four lanes) added in index order
SUM_SEQUENTIAL_UPTO = 43
SUM_BLOCK_FROM = 89000
SUM_BLOCK = 4096

_NUMERIC = ("int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64")
_ARITH = ("int8", "uint8", "int16", "uint16", "int32", "uint32", "single", "double")
_MODES = ("class", "cast")
_OPS = ("eq", "ne", "lt", "le", "gt", "ge")


@functools.cache
def _core():
    """The extension module with the ops, or None (not built, or built without them)."""
    from petakit5d import _ext

    if not _ext.available("core"):
        return None
    m = _ext.load("core")
    return m if hasattr(m, "ops_convert") else None


def extension_available() -> bool:
    """Whether the compiled primitives are there (else the NumPy references run)."""
    return _core() is not None


# -- layout ------------------------------------------------------------------------------------


def _contig(x: np.ndarray) -> np.ndarray:
    """x itself when contiguous (C order, or a MATLAB-layout Fortran array inside a kernel),
    else a C-ordered copy"""
    if x.flags.c_contiguous or x.flags.f_contiguous:
        return x
    return np.ascontiguousarray(x)


def _flat(x: np.ndarray) -> np.ndarray:
    """the 1-D view of a contiguous array in memory order"""
    return x.ravel(order="K")


def _empty_like(x: np.ndarray, dtype) -> np.ndarray:
    """an array of the shape and memory order of x"""
    return np.empty_like(x, dtype=dtype)


def _like_layout(m: np.ndarray, x: np.ndarray) -> np.ndarray:
    """m in the memory order of x (a copy only when the layouts differ)"""
    if (x.flags.c_contiguous and m.flags.c_contiguous) or (
        x.flags.f_contiguous and m.flags.f_contiguous
    ):
        return m
    r = np.empty_like(x, dtype=m.dtype)
    r[...] = m
    return r


def _out_array(x: np.ndarray, dtype, out: np.ndarray | None) -> np.ndarray:
    if out is None:
        return _empty_like(x, dtype)
    if not isinstance(out, np.ndarray) or out.dtype != np.dtype(dtype) or out.shape != x.shape:
        raise ValueError(f"out must be a {np.dtype(dtype)} array of shape {x.shape}")
    if not out.flags.writeable:
        raise ValueError("out is read-only")
    fx, cx = x.flags.f_contiguous, x.flags.c_contiguous
    if not ((fx and out.flags.f_contiguous) or (cx and out.flags.c_contiguous)):
        raise ValueError("out must be contiguous in the memory order of x")
    if np.may_share_memory(x, out) and not (
        out.dtype == x.dtype
        and out.__array_interface__["data"][0] == x.__array_interface__["data"][0]
    ):
        raise ValueError("out overlaps x (only out=x, of the same class, may share memory)")
    return out


def _class_name(cls) -> str:
    try:
        name = matlab_class(cls)
    except TypeError:
        name = None
    if name not in CLASS_DTYPES or name == "logical":
        raise ValueError(f"unsupported target class {cls!r}")
    return name


def _array(x) -> np.ndarray:
    a = np.asarray(x)
    if a.dtype.kind not in "biuf":
        raise TypeError(f"real numeric or logical array required, got {a.dtype}")
    return a


def _mode(mode: str) -> str:
    if mode not in _MODES:
        raise ValueError(f"mode must be 'class' or 'cast', got {mode!r}")
    return mode


# -- convert -----------------------------------------------------------------------------------


def _exact_rounding(a: np.ndarray, mode: str, scalar_fast: bool) -> bool:
    """whether a conversion rounds exactly (else floor(|x| + 0.5) in single, for single data to
    8 and 16-bit classes, and through double for 64-bit integers to single): to_class rounds a
    scalar exactly unless scalar_fast"""
    return mode == "cast" or (a.size < 2 and not scalar_fast)


def convert_ref(x, cls, *, mode: str = "class", scalar_fast: bool = False) -> np.ndarray:
    """NumPy reference of ``convert``: ``dtypes.to_class`` (mode 'class') or ``dtypes.cast``."""
    a = _array(x)
    cls = _class_name(cls)
    if _mode(mode) == "cast":
        r = dtypes.cast(a, cls)
    elif cls == "single" and a.dtype.kind in "iu" and a.dtype.itemsize == 8 and a.size >= 2:
        # MATLAB's single() of a 64-bit integer array rounds through double (cast doesn't)
        r = dtypes.to_class(dtypes.to_class(a, "double"), "single")
    else:
        r = dtypes.to_class(a, cls, scalar_fast=scalar_fast)
    return np.asarray(r)


def convert(
    x,
    cls,
    *,
    mode: str = "class",
    scalar_fast: bool = False,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """MATLAB's ``cls(x)`` (mode 'class', ``uint16(x)`` in MATLAB code) or ``cast(x, cls)``
    (mode 'cast'): round half away from zero (single data to int8, uint8, int16 and uint16 in
    mode 'class': ``floor(|x| + 0.5)`` in single), saturate, NaN to 0. In mode 'class' a 64-bit
    integer array becomes single through double (rounded twice, as MATLAB's ``single(x)``; its
    ``cast`` rounds once). ``cls`` is a MATLAB class name or a NumPy dtype; logical and real
    numeric inputs."""
    a = _contig(_array(x))
    cls = _class_name(cls)
    exact = _exact_rounding(a, _mode(mode), scalar_fast)
    m = _core()
    if m is None:
        r = convert_ref(a, cls, mode="cast" if exact else "class", scalar_fast=True)
        if out is None:
            return r
        out = _out_array(a, CLASS_DTYPES[cls], out)
        out[...] = r
        return out
    o = _out_array(a, CLASS_DTYPES[cls], out)
    m.ops_convert(_flat(a), _flat(o), exact, resolve_threads(num_threads))
    return o


# -- scale, offset, clamp, mask ----------------------------------------------------------------


def _bound(v, cls: str):
    """a min/max bound in the class of the data (None: no step; NaN: no step, as MATLAB's min
    and max ignore NaN)"""
    if v is None:
        return None
    v = float(v)
    if math.isnan(v):
        return None
    if cls in ("single", "double"):
        return v if cls == "double" else float(np.float32(v))
    return float(dtypes.cast(v, cls))


def _mmin(y: np.ndarray, b: float) -> np.ndarray:
    """MATLAB's min(y, b) for a bound b of y's class: y on a tie, b for NaN y"""
    bb = y.dtype.type(b)
    return np.where(y <= bb, y, bb).astype(y.dtype)


def _mmax(y: np.ndarray, b: float) -> np.ndarray:
    bb = y.dtype.type(b)
    return np.where(y >= bb, y, bb).astype(y.dtype)


def _check_arith(a: np.ndarray, mask) -> str:
    cls = matlab_class(a)
    if cls not in _ARITH:
        raise TypeError(
            f"scale_offset: {cls} data is not supported (single, double, 8 to 32-bit integers)"
        )
    if mask is not None:
        if cls not in ("single", "double"):
            raise TypeError(
                "Integers can only be combined with integers of the same class, or scalar "
                f"doubles ({cls}, logical) (MATLAB:mixedClasses)"
            )
        mk = np.asarray(mask)
        if mk.dtype != np.bool_ or mk.shape != a.shape:
            raise ValueError(f"mask must be a logical array of shape {a.shape}")
    return cls


def scale_offset_ref(
    x,
    scale=None,
    offset=None,
    *,
    hi=None,
    lo=None,
    mask=None,
    out_class=None,
    mode: str = "class",
    scalar_fast: bool = False,
) -> np.ndarray:
    """NumPy reference of ``scale_offset`` (``dtypes.times``, ``plus``, MATLAB's min and max)."""
    a = _array(x)
    cls = _check_arith(a, mask)
    y = a
    if scale is not None:
        y = np.asarray(dtypes.times(y, float(scale)))
    if offset is not None:
        y = np.asarray(dtypes.plus(y, float(offset)))
    b = _bound(hi, cls)
    if b is not None:
        y = _mmin(y, b)
    b = _bound(lo, cls)
    if b is not None:
        y = _mmax(y, b)
    if mask is not None:
        y = np.asarray(dtypes.times(y, np.asarray(mask)))
    if y is a:
        y = a.copy()
    if out_class is not None:
        y = convert_ref(y, out_class, mode=mode, scalar_fast=scalar_fast)
    return y


def scale_offset(
    x,
    scale: float | None = None,
    offset: float | None = None,
    *,
    hi: float | None = None,
    lo: float | None = None,
    mask=None,
    out_class=None,
    mode: str = "class",
    scalar_fast: bool = False,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """MATLAB's ``y = x * scale + offset; y = min(y, hi); y = max(y, lo); y = y .* mask;
    y = out_class(y)`` in one pass, each step only when its argument is given.

    ``scale``, ``offset``, ``hi`` and ``lo`` are double scalars. Single data with a double
    scalar is computed in double and rounded to single at every step; 8 to 32-bit integer data
    is computed in double, rounded half away from zero and saturated at every step (MATLAB's
    rules). A NaN bound is no step (MATLAB's min and max ignore NaN); on a tie the data is kept,
    and a NaN value gives the bound. ``mask`` (logical, the shape of ``x``) needs single or
    double data, as in MATLAB. ``out_class`` converts at the end as ``convert`` with ``mode``.
    """
    a = _contig(_array(x))
    cls = _check_arith(a, mask)
    if out_class is not None:
        ocls = _class_name(out_class)
    else:
        ocls = cls
    exact = _exact_rounding(a, _mode(mode), scalar_fast)
    m = _core()
    if m is None:
        r = scale_offset_ref(
            a,
            scale,
            offset,
            hi=hi,
            lo=lo,
            mask=mask,
            out_class=ocls if out_class is not None else None,
            mode="cast" if exact else "class",
            scalar_fast=True,
        )
        if out is None:
            return r
        out = _out_array(a, CLASS_DTYPES[ocls], out)
        out[...] = r
        return out
    mk = None
    if mask is not None:
        mk = _flat(_like_layout(np.asarray(mask), a))
    bh, bl = _bound(hi, cls), _bound(lo, cls)
    o = _out_array(a, CLASS_DTYPES[ocls], out)
    m.ops_affine(
        _flat(a),
        _flat(o),
        scale is not None,
        float(scale) if scale is not None else 1.0,
        offset is not None,
        float(offset) if offset is not None else 0.0,
        bh is not None,
        bh if bh is not None else 0.0,
        bl is not None,
        bl if bl is not None else 0.0,
        mk,
        exact,
        resolve_threads(num_threads),
    )
    return o


def clamp_ref(x, lo=None, hi=None) -> np.ndarray:
    """NumPy reference of ``clamp``."""
    return scale_offset_ref(x, hi=hi, lo=lo)


def clamp(
    x,
    lo: float | None = None,
    hi: float | None = None,
    *,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """MATLAB's ``max(min(x, hi), lo)`` (each only when given) with double scalar bounds, in the
    class of ``x``: NaN gives ``hi`` (else ``lo``), ties keep ``x`` (``max(-0, 0)`` is -0), a NaN
    bound changes nothing."""
    return scale_offset(x, hi=hi, lo=lo, out=out, num_threads=num_threads)


# -- compare -----------------------------------------------------------------------------------


def _op(op: str) -> int:
    if op not in _OPS:
        raise ValueError(f"op must be one of {_OPS}, got {op!r}")
    return _OPS.index(op)


def _int64_bounds(value: float, signed: bool) -> tuple[int, int, bool]:
    """floor and ceil of a finite double, clipped to just outside the class's range, and whether
    it is integral"""
    lo, hi = (-(2**63), 2**63 - 1) if signed else (0, 2**64 - 1)
    f, c = math.floor(value), math.ceil(value)
    return max(lo - 1, min(hi + 1, f)), max(lo - 1, min(hi + 1, c)), f == c


def compare_ref(x, op: str, value: float) -> np.ndarray:
    """NumPy reference of ``compare``: the exact comparison of each value with ``value``."""
    a = _array(x)
    k = _op(op)
    v = float(value)
    if math.isnan(v):
        return np.full(a.shape, op == "ne", dtype=bool)
    if a.dtype.kind in "iu" and a.dtype.itemsize == 8:
        if math.isinf(v):
            r = {
                "eq": False,
                "ne": True,
                "lt": v > 0,
                "le": v > 0,
                "gt": v < 0,
                "ge": v < 0,
            }[op]
            return np.full(a.shape, r, dtype=bool)
        f, c, integral = _int64_bounds(v, a.dtype.kind == "i")
        w = a.astype(object)
        if k < 2:
            eq = np.asarray(w == f, dtype=bool) & integral
            return eq if k == 0 else ~eq
        r = [None, None, w < c, w <= f, w > f, w >= c][k]
        return np.asarray(r, dtype=bool).reshape(a.shape)
    d = a.astype(np.float64)
    return [d == v, d != v, d < v, d <= v, d > v, d >= v][k]


def compare(x, op: str, value: float, *, num_threads: int | None = None) -> np.ndarray:
    """MATLAB's ``x == value``, ``~=``, ``<``, ``<=``, ``>``, ``>=`` (op 'eq', 'ne', 'lt', 'le',
    'gt', 'ge') with a double scalar: each value of ``x`` compared exactly with ``value``
    (NaN: only 'ne' is true). A logical array of the shape and layout of ``x``."""
    a = _contig(_array(x))
    k = _op(op)
    m = _core()
    if m is None:
        return compare_ref(a, op, value)
    o = _empty_like(a, np.bool_)
    m.ops_compare(_flat(a), _flat(o), k, float(value), resolve_threads(num_threads))
    return o


# -- sum ---------------------------------------------------------------------------------------


def _seq_sum(v: np.ndarray) -> float:
    """((v1 + v2) + v3) + ... in double (cumsum adds in order; -0 for an empty lane)"""
    if not v.size:
        return -0.0
    with np.errstate(invalid="ignore", over="ignore"):
        return float(np.cumsum(v)[-1])


def _add_first(a: float, b: float) -> float:
    """a + b with the NaN of a, else of b, quieted (x86's addition of NaNs)"""
    if a != a:
        return a + a
    if b != b:
        return b + b
    return a + b


def _lanes_exact(v: list[float], tail_first: bool) -> float:
    acc = [-0.0, -0.0, -0.0, -0.0]
    main = len(v) - len(v) % 4
    for t, x in enumerate(v):
        k = t & 3
        acc[k] = _add_first(x, acc[k]) if t >= main and tail_first else _add_first(acc[k], x)
    return _add_first(_add_first(_add_first(acc[0], acc[1]), acc[2]), acc[3])


def _nan_exact_sum(v: np.ndarray, single: bool) -> float:
    """the sum with MATLAB's operand order, which decides the payload of a NaN result: the
    running sums are the first operand, except in the lanes' tails (the last mod(m, 4) elements
    of the four-lane sum, and of the last block of single data) and for the block sums of
    double data (measured in R2024b)"""
    n = v.size
    w = v.tolist()
    if n <= SUM_SEQUENTIAL_UPTO:
        s = -0.0
        for x in w:
            s = _add_first(s, x)
        return s
    if n < SUM_BLOCK_FROM:
        return _lanes_exact(w, True)
    s = -0.0
    for i0 in range(0, n, SUM_BLOCK):
        b = _lanes_exact(w[i0 : i0 + SUM_BLOCK], single)
        s = _add_first(s, b) if single else _add_first(b, s)
    return s


def _matlab_values(a: np.ndarray, matlab_layout: bool) -> np.ndarray:
    """the values of a in MATLAB's linear order (``vol(:)``)"""
    if matlab_layout:
        return a.reshape(-1, order="F")
    return matlab_order(a).reshape(-1)


def _check_sum_class(a: np.ndarray) -> None:
    if a.dtype.kind in "iu" and a.dtype.itemsize == 8:
        raise TypeError("sum_double: 64-bit integers are not supported")


def sum_double_ref(x, *, matlab_layout: bool = False) -> float:
    """NumPy reference of ``sum_double``."""
    a = _array(x)
    _check_sum_class(a)
    with np.errstate(invalid="ignore"):
        v = _matlab_values(a, matlab_layout).astype(np.float64)
    s = _fast_sum(v)
    if s != s:
        return _nan_exact_sum(v, a.dtype == np.float32)
    return s


def _fast_sum(v: np.ndarray) -> float:
    n = v.size
    if n == 0:
        return 0.0
    if n <= SUM_SEQUENTIAL_UPTO:
        return _seq_sum(v)
    if n < SUM_BLOCK_FROM:
        lanes = [_seq_sum(v[k::4]) for k in range(4)]
        with np.errstate(invalid="ignore", over="ignore"):
            return float(((np.float64(lanes[0]) + lanes[1]) + lanes[2]) + lanes[3])
    nb = -(-n // SUM_BLOCK)
    p = np.full(nb * SUM_BLOCK, -0.0)
    p[:n] = v
    with np.errstate(invalid="ignore", over="ignore"):
        lanes = np.cumsum(p.reshape(4, SUM_BLOCK // 4, nb, order="F"), axis=1)[:, -1, :]
        bs = ((lanes[0] + lanes[1]) + lanes[2]) + lanes[3]
    return _seq_sum(bs)


def sum_double(x, *, matlab_layout: bool = False, num_threads: int | None = None) -> float:
    """``sum(v(:), 'double')`` in MATLAB's order at one thread, whatever the thread count
    (decon_sum_double), where ``v`` is MATLAB's array of ``x`` (``to_matlab(x)``).

    The elements are added in MATLAB's linear order ``t = y + ny*(x + nx*z)`` (y fastest; the C
    order of ``matlab_order(x)``): in index order up to 43 elements; in four lanes (element t in
    lane mod(t, 4), each lane in increasing t, added as ((l1 + l2) + l3) + l4) up to 88,999;
    from 89,000 on in blocks of 4,096 elements (block t // 4096), each in four lanes, the block
    sums added in index order. A sum of -0s is -0, an empty sum 0; a NaN sum has MATLAB's
    payload (that of the first NaN in MATLAB's operand order, ``_nan_exact_sum``). Logical,
    single, double and 8 to 32-bit integer data.

    A C-contiguous ``x`` is read in place (strided along y); ``matlab_layout=True`` takes
    MATLAB's own array ``(y, x, z, ...)`` instead (as a kernel that works in MATLAB's layout
    holds it; a Fortran-contiguous one is read in memory order).
    """
    a = _array(x)
    _check_sum_class(a)
    m = _core()
    if m is None:
        return sum_double_ref(a, matlab_layout=matlab_layout)
    nt = resolve_threads(num_threads)
    if a.size == 0:
        return 0.0
    if matlab_layout:
        v = a.reshape(-1, order="F")  # a view when a is Fortran-contiguous
        return float(m.ops_sum_double(np.ascontiguousarray(v), nt))
    if a.ndim < 2:
        return float(m.ops_sum_double(np.ascontiguousarray(a.reshape(-1)), nt))
    a = np.ascontiguousarray(a)
    return float(m.ops_sum_double_yx(a.reshape(-1, a.shape[-2], a.shape[-1]), nt))
