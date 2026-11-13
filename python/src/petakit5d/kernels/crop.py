"""Crop a box out of an array and write a region into the box of an array: the kernels of
``crop3d_mex``, ``crop4d_mex``, ``indexing3d_mex``, ``indexing4d_mex`` and
``indexing4d_crop_mex`` (``microscopeDataProcessing/crop/mex``), and the MATLAB wrappers
``crop3d``, ``crop4d`` and ``indexing4d``.

The Python side of ``petakit5d._ext.crop`` (``python/cpp/crop``). Two levels:

- the kernels, with the port's arrays and boxes (``docs/conventions.md``): volumes ``(z, y,
  x)``, images ``(y, x)``, 4-D arrays ``(t, z, y, x)``; boxes ``(z0, y0, x0, z1, y1, x1)`` (2-D
  ``(y0, x0, y1, x1)``, 4-D ``(t0, z0, y0, x0, t1, z1, y1, x1)``), 0-based and end-exclusive,
  ``None`` for an open end:

  - ``crop_3d(vol, box)``, ``crop_4d(vol, box)``: a new C-contiguous array with the box of
    ``vol`` (``vol`` may have any layout; ``out=`` writes into a given array instead, to reuse
    a buffer over frames);
  - ``indexing_3d(out, box, region)``, ``indexing_4d(out, box, region, region_box=None)``:
    write ``region`` (or its box ``region_box``) into the box of ``out``, **in place**, and
    return ``out``. ``out`` must be a writeable C-contiguous array: it is never copied, and no
    other array is written. A ``region`` that shares memory with ``out`` is copied first.

- the MATLAB functions with MATLAB's box values (``[ymin xmin zmin ymax xmax zmax]``, 4-D
  ``[ymin xmin zmin tmin ymax xmax zmax tmax]``, 1-based, inclusive, any numbers) and
  ``(z, y, x)`` arrays, for the MATLAB-name layer and the oracle: ``crop3d_matlab``,
  ``crop4d_matlab`` and ``indexing4d_matlab`` are the ``.m`` wrappers (the mex's checks, in
  its order, with its error identifiers and messages, its truncation of non-integer bounds and
  its results; where the mex raises, the wrapper's MATLAB indexing code with MATLAB's results
  and errors; ``mex_only=True``: the mex alone), ``indexing3d_matlab`` and
  ``indexing4d_crop_matlab`` the mex functions that have no wrapper. Like the mex functions,
  the ``indexing*`` ones write into their first argument in place (MATLAB's mex writes into
  the input's memory, which ``feather_blending_3d``'s fall-back relies on); where the
  wrapper's MATLAB code grows the array (a box beyond it), the grown array is a new one.

Every NumPy class is copied by the kernels (the bytes of each element, as memcpy copies them
in the mex); the MATLAB functions take the mex's classes (logical, uint8, uint16, single,
double) and, where the wrapper falls back, every class MATLAB indexing takes. A complex array
reaches the mex as its real part (the mex's separate-complex API reads only the real part): a
crop of complex data is real, and an indexing mex leaves complex data unchanged (MATLAB passes
the mex a copy; known issue, kept). ``num_threads``
(None or 0: the CPUs this process may use) never changes a result. Errors raise
``MatlabError`` with MATLAB's identifier (``identifier``).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from types import ModuleType
from typing import Any

import numpy as np

__all__ = [
    "MatlabError",
    "crop3d_matlab",
    "crop4d_matlab",
    "crop_3d",
    "crop_4d",
    "extension_available",
    "indexing3d_matlab",
    "indexing4d_crop_matlab",
    "indexing4d_matlab",
    "indexing_3d",
    "indexing_4d",
]

CROP_INPUT = "crop:inputError"
CROP_TYPE = "crop:dataTypeError"
INDEXING_INPUT = "indexing:inputError"
# crop4d_mex and the indexing mex functions name their class error so
TIFF_TYPE = "tiff:dataTypeError"

# the classes of the mex functions
_MEX_DTYPES = frozenset(
    np.dtype(t) for t in (np.bool_, np.uint8, np.uint16, np.float32, np.float64)
)
_COMPLEX = frozenset(np.dtype(t) for t in (np.complex64, np.complex128))
_U64 = (1 << 64) - 1
# below this many bytes one thread copies (a team costs more than the copy); above, at least
# this many bytes per thread
_SERIAL_BYTES = 1 << 20
_BYTES_PER_THREAD = 1 << 18
# bytes of one memcpy at most: one large call uses non-temporal stores, which into a fresh
# array (its pages faulted in by the copy) measured 1.3 times slower at one thread than calls
# of 16 KB to 1 MB (K01's micro benchmark)
_PIECE = 1 << 18
# per-region thread binding: off. Bound teams measured the same at 1 and 24 threads and 5%
# slower at 4 threads for a crop (micro jobs 4572759, 4572760); the mex's team is not bound
_BIND = False

_ext_module: ModuleType | None = None
_ext_missing = False
_available: int | None = None


class MatlabError(ValueError):
    """An error MATLAB raises, with its identifier (``'crop:inputError'``, ...)."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _ext() -> ModuleType | None:
    global _ext_module, _ext_missing
    if _ext_module is None and not _ext_missing:
        from petakit5d import _ext as ext

        try:
            _ext_module = ext.load("crop")
        except ImportError:
            _ext_missing = True
    return _ext_module


def extension_available() -> bool:
    """Whether the C++ kernel is built (otherwise NumPy copies)."""
    return _ext() is not None


def _threads(num_threads: int | None) -> int:
    """``core.threads.resolve_threads`` with the CPUs counted once per process"""
    global _available
    if _available is None:
        from petakit5d.core.threads import available_cpus

        _available = available_cpus()
    if num_threads is None or num_threads == 0:
        return _available
    n = int(num_threads)
    if n < 0:
        raise ValueError(f"thread count must be positive, got {n}")
    return min(n, _available)


# ---------------------------------------------------------------------------------------------
# the copy
# ---------------------------------------------------------------------------------------------


def _pad4(shape: Sequence[int]) -> list[int]:
    s = [int(v) for v in shape]
    return [1] * (4 - len(s)) + s


def _copy(
    src: np.ndarray,
    s0: Sequence[int],
    dst: np.ndarray,
    d0: Sequence[int],
    ext: Sequence[int],
    num_threads: int | None,
    reference: bool = False,
) -> None:
    """dst[d0 : d0 + ext] = src[s0 : s0 + ext] for arrays of up to 4 dimensions (boxes checked;
    dst C-contiguous, src any layout, the same dtype; they do not overlap)"""
    n = len(ext)
    nbytes = math.prod(int(e) for e in ext) * src.dtype.itemsize
    if nbytes == 0:
        return
    k = None if reference or not src.flags.c_contiguous else _ext()
    if k is None or src.dtype.itemsize == 0:
        ss = tuple(slice(int(a), int(a) + int(e)) for a, e in zip(s0, ext, strict=True))
        ds = tuple(slice(int(a), int(a) + int(e)) for a, e in zip(d0, ext, strict=True))
        np.copyto(dst[ds], src[ss], casting="no")
        return
    nt = 1
    if nbytes >= _SERIAL_BYTES:
        nt = max(1, min(_threads(num_threads), nbytes // _BYTES_PER_THREAD))
    pad = [0] * (4 - n)
    k.copy_box(
        src.reshape(-1).view(np.uint8),
        _pad4(src.shape),
        pad + [int(v) for v in s0],
        dst.reshape(-1).view(np.uint8),
        _pad4(dst.shape),
        pad + [int(v) for v in d0],
        [1] * (4 - n) + [int(v) for v in ext],
        src.dtype.itemsize,
        nt,
        _BIND,
        _PIECE,
    )


def _array(x: Any, name: str) -> np.ndarray:
    if not isinstance(x, np.ndarray):
        raise TypeError(f"{name} must be a NumPy array, got {type(x).__name__}")
    if x.dtype.hasobject:
        raise TypeError(f"{name} must not hold Python objects ({x.dtype})")
    return x


def _writeable(out: Any, name: str = "out") -> np.ndarray:
    o = _array(out, name)
    if not o.flags.c_contiguous:
        hint = " (it is Fortran-ordered: MATLAB's layout? use core.arrays.from_matlab)"
        raise ValueError(f"{name} must be C-contiguous{hint if o.flags.f_contiguous else ''}")
    if not o.flags.writeable:
        raise ValueError(f"{name} must be writeable")
    return o


def _int_value(v: Any, what: str, ident: str) -> int:
    if isinstance(v, bool | np.bool_):
        raise MatlabError(ident, f"{what} must be integers, got {v!r}")
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise MatlabError(ident, f"{what} must be integers, got {v!r}") from None
    if not math.isfinite(f) or f != math.floor(f):
        raise MatlabError(ident, f"{what} must be integers, got {v!r}")
    return int(v) if isinstance(v, int | np.integer) else int(f)


def _box(box: Any, shape: Sequence[int], ident: str, what: str = "box") -> list[int]:
    """a Python box of an array of shape: starts then stops, integers (None: open), within"""
    n = len(shape)
    vals = list(box) if not isinstance(box, np.ndarray) else list(box.reshape(-1))
    if len(vals) != 2 * n:
        raise MatlabError(ident, f"{what} must have {2 * n} values for {n} dimensions")
    b = []
    for i, v in enumerate(vals):
        if v is None:
            v = 0 if i < n else shape[i - n]
        b.append(_int_value(v, f"the values of {what}", ident))
    for i in range(n):
        if not 0 <= b[i] <= b[i + n] <= shape[i]:
            raise MatlabError(
                ident, f"{what} {tuple(b)} is outside the array of shape {tuple(shape)}"
            )
    return b


def _result(out: np.ndarray | None, shape: tuple[int, ...], dtype: np.dtype) -> np.ndarray:
    if out is None:
        return np.empty(shape, dtype=dtype)
    o = _writeable(out)
    if o.dtype != dtype or o.shape != shape:
        raise ValueError(f"out must be a {dtype} array of shape {shape}, got {o.dtype} {o.shape}")
    return o


# ---------------------------------------------------------------------------------------------
# the kernels
# ---------------------------------------------------------------------------------------------


def _crop(
    vol: Any,
    box: Any,
    ndims: tuple[int, ...],
    out: np.ndarray | None,
    num_threads: int | None,
) -> np.ndarray:
    a = _array(vol, "vol")
    if a.ndim not in ndims:
        raise MatlabError(CROP_INPUT, f"vol must have {' or '.join(map(str, ndims))} dimensions")
    b = _box(box, a.shape, CROP_INPUT)
    n = a.ndim
    ext = tuple(b[i + n] - b[i] for i in range(n))
    res = _result(out, ext, a.dtype)
    if out is not None and np.may_share_memory(res, a):
        raise ValueError("out must not overlap vol")
    _copy(a, b[:n], res, (0,) * n, ext, num_threads)
    return res


def crop_3d(
    vol: np.ndarray,
    box: Sequence[int | None],
    *,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """``vol[z0:z1, y0:y1, x0:x1]`` of a ``(z, y, x)`` volume (an image ``(y, x)`` takes
    ``(y0, x0, y1, x1)``) as a new C-contiguous array, the copy of ``crop3d_mex``. Empty ranges
    give empty arrays. ``out``: a C-contiguous array of the box's shape and ``vol``'s class to
    write into (returned). ``MatlabError`` ``crop:inputError`` for a box that is not integral,
    reversed or outside ``vol``."""
    return _crop(vol, box, (2, 3), out, num_threads)


def crop_4d(
    vol: np.ndarray,
    box: Sequence[int | None],
    *,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """``vol[t0:t1, z0:z1, y0:y1, x0:x1]`` of a ``(t, z, y, x)`` array, box ``(t0, z0, y0,
    x0, t1, z1, y1, x1)``, as a new C-contiguous array, the copy of ``crop4d_mex`` (empty
    ranges allowed here; ``crop4d_matlab`` keeps the mex's refusal). ``out`` as in
    ``crop_3d``."""
    return _crop(vol, box, (4,), out, num_threads)


def _indexing(
    out: Any,
    box: Any,
    region: Any,
    region_box: Any,
    ndims: tuple[int, ...],
    num_threads: int | None,
) -> np.ndarray:
    o = _writeable(out)
    if o.ndim not in ndims:
        raise MatlabError(
            INDEXING_INPUT, f"out must have {' or '.join(map(str, ndims))} dimensions"
        )
    r = _array(region, "region")
    if r.dtype != o.dtype:
        raise MatlabError(
            INDEXING_INPUT, "The data type of the region does not match that of the data!"
        )
    n = o.ndim
    b = _box(box, o.shape, INDEXING_INPUT)
    ext = tuple(b[i + n] - b[i] for i in range(n))
    if r.ndim != n:
        raise MatlabError(INDEXING_INPUT, f"region must have {n} dimensions")
    if region_box is None:
        rb = [0] * n + list(r.shape)
    else:
        rb = _box(region_box, r.shape, INDEXING_INPUT, "region_box")
    if tuple(rb[i + n] - rb[i] for i in range(n)) != ext:
        raise MatlabError(INDEXING_INPUT, "Subregion size does not match the bounding box size")
    if np.may_share_memory(o, r):
        r = r.copy()
    _copy(r, rb[:n], o, b[:n], ext, num_threads)
    return o


def indexing_3d(
    out: np.ndarray,
    box: Sequence[int | None],
    region: np.ndarray,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """``out[z0:z1, y0:y1, x0:x1] = region``, in place (the copy of ``indexing3d_mex``):
    ``out`` is a writeable C-contiguous ``(z, y, x)`` volume or ``(y, x)`` image (never
    copied), ``region`` an array of its class and of the box's shape (any layout). Returns
    ``out``. ``MatlabError`` ``indexing:inputError`` for a bad box, size or class."""
    return _indexing(out, box, region, None, (2, 3), num_threads)


def indexing_4d(
    out: np.ndarray,
    box: Sequence[int | None],
    region: np.ndarray,
    region_box: Sequence[int | None] | None = None,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """``out[t0:t1, z0:z1, y0:y1, x0:x1] = region[region_box]``, in place (the copy of
    ``indexing4d_mex`` and ``indexing4d_crop_mex``): ``out`` a writeable C-contiguous ``(t, z,
    y, x)`` array, ``region`` a 4-D array of its class, ``region_box`` (default: all of
    ``region``) a box of ``region`` of the same extent as ``box``. Returns ``out``."""
    return _indexing(out, box, region, region_box, (4,), num_threads)


# ---------------------------------------------------------------------------------------------
# MATLAB's values: the mex's arguments
# ---------------------------------------------------------------------------------------------


def _is_sparse(x: Any) -> bool:
    return type(x).__module__.startswith("scipy.sparse")


def _mx(x: Any) -> tuple[str, tuple[int, ...], np.ndarray, bool]:
    """MATLAB's view of a parameter value: class, size, the values in memory order (the real
    part, as double) and whether it is sparse. Lists of numbers are double rows (of bools,
    logical rows); NumPy arrays keep their class, a 1-D array is a row, matrices keep their
    shape (M(i, j) = P[i - 1, j - 1])."""
    from petakit5d.core.arrays import matlab_class

    if _is_sparse(x):
        d = np.asarray(x.toarray())
        return matlab_class(d), tuple(d.shape), d.real.astype(np.float64).ravel(order="F"), True
    if x is None:
        return "double", (0, 0), np.zeros(0), False
    if isinstance(x, list | tuple):
        if x and all(isinstance(v, bool | np.bool_) for v in x):
            return "logical", (1, len(x)), np.asarray(x, dtype=np.float64), False
        if all(isinstance(v, int | float | np.integer | np.floating) for v in x):
            return "double", (1, len(x)) if x else (0, 0), np.asarray(x, np.float64), False
        raise TypeError(f"not a MATLAB number array: {x!r}")
    if isinstance(x, bool | int | float) and not isinstance(x, np.generic):
        return ("logical" if isinstance(x, bool) else "double"), (1, 1), np.array([x], float), False
    a = np.asarray(x)
    cls = matlab_class(a.dtype)
    size = (1, 1) if a.ndim == 0 else (1, a.size) if a.ndim == 1 else tuple(a.shape)
    return cls, size, a.real.astype(np.float64).ravel(order="F"), False


def _floor(v: float) -> float:
    """C's floor (NaN and Inf pass through)"""
    return float(np.floor(v))


def _mex_u64(v: float) -> int:
    """``(uint64_t)v`` of a double as GCC's x86-64 code converts it (values beyond the range
    as cvttsd2si gives them)"""
    two63 = 9223372036854775808.0
    if v >= two63:
        w = v - two63
        return (int(w) if w < two63 else 1 << 63) ^ (1 << 63)
    if v != v or v < -two63:
        return 1 << 63
    return int(v) & _U64


def _mex_bbox(bbox: Any, n: int, ident: str, name: str = "bbox") -> np.ndarray:
    """the mex's checks of a box argument: mxGetN (columns) of n, a full double array of at
    least n elements; its first n values in memory order"""
    cls, size, vals, sparse = _mx(bbox)
    if math.prod(size[1:]) != n:
        raise MatlabError(ident, f"Input range for {name} is not {n}")
    if cls != "double" or sparse or vals.size < n:
        raise MatlabError(ident, f"The {name} must be a full double array")
    return vals[:n]


def _mex_data(a: Any, what: str, ident: str) -> np.ndarray:
    if _is_sparse(a):
        raise MatlabError(ident, what)
    return _array(a, "the data")


def _msize(a: np.ndarray) -> list[int]:
    from petakit5d.core.arrays import size_to_matlab

    return list(size_to_matlab(a.shape))


def _python_shape(msize: Sequence[int]) -> tuple[int, ...]:
    """our shape of a MATLAB array of size msize (trailing singletons beyond 2 dropped)"""
    from petakit5d.core.arrays import size_from_matlab

    m = [int(v) for v in msize]
    while len(m) > 2 and m[-1] == 1:
        m.pop()
    return size_from_matlab(m)


def _view(a: np.ndarray, n: int) -> np.ndarray:
    """the n-D (..., z, y, x) view of an array of at most n MATLAB dimensions"""
    if a.ndim <= n:
        return a.reshape((1,) * (n - a.ndim) + a.shape)
    return a.reshape(a.shape[a.ndim - n :])


def _mex_class(a: np.ndarray, ident: str) -> np.ndarray:
    """the array the mex copies: the mex's classes as they are, a complex single or double
    array as its real part (the mex reads mxGetPr); other classes raise"""
    v = _real(a)
    if v.dtype not in _MEX_DTYPES:
        raise MatlabError(ident, "Data type not suppported")
    return v


def _crop3d_mex(vol: Any, bbox: Any, num_threads: int | None) -> np.ndarray:
    a = _mex_data(vol, "The input must be a full array", CROP_INPUT)
    ms = _msize(a)
    b = _mex_bbox(bbox, 6, CROP_INPUT)
    dims = [ms[0], ms[1], ms[2] if len(ms) > 2 else 1]
    for i in range(3):
        if not b[i] >= 1:
            raise MatlabError(CROP_INPUT, "Lower bounds must be at least 1")
        if not _floor(b[i + 3]) <= dims[i]:
            raise MatlabError(CROP_INPUT, "Upper bound is invalid")
        if not _floor(b[i + 3]) >= _floor(b[i]) - 1:
            raise MatlabError(CROP_INPUT, "Upper bound is smaller than the lower bound minus 1")
    lo = [_mex_u64(float(v)) - 1 for v in b[:3]]  # y, x, z
    hi = [_mex_u64(float(v)) for v in b[3:]]
    src = _mex_class(a, CROP_TYPE)
    # the mex reads MATLAB's first three dimensions: the first (z, y, x) volume
    v = src.reshape((1,) * (3 - src.ndim) + src.shape) if src.ndim <= 3 else src
    while v.ndim > 3:
        v = v[0]
    ext = (hi[2] - lo[2], hi[0] - lo[0], hi[1] - lo[1])
    res = np.empty(ext, dtype=v.dtype)
    _copy(v, (lo[2], lo[0], lo[1]), res, (0, 0, 0), ext, num_threads)
    return res.reshape(_python_shape((ext[1], ext[2], ext[0])))


def _crop4d_mex(vol: Any, bbox: Any, num_threads: int | None) -> np.ndarray:
    a = _mex_data(vol, "The input must be a full array", CROP_INPUT)
    ms = _msize(a)
    if len(ms) > 4:
        raise MatlabError(CROP_INPUT, "The input must have at most 4 dimensions")
    dims = ms + [1] * (4 - len(ms))
    b = _mex_bbox(bbox, 8, CROP_INPUT)
    for i in range(4):
        if not b[i] >= 1:
            raise MatlabError(CROP_INPUT, "Lower bounds must be at least 1")
        if not _floor(b[i + 4]) <= dims[i]:
            raise MatlabError(CROP_INPUT, "Upper bound for bbox is invalid")
    lo = [(_mex_u64(float(v)) - 1) & _U64 for v in b[:4]]  # y, x, z, t
    hi = [_mex_u64(float(v)) for v in b[4:]]
    for i in range(4):
        if hi[i] <= lo[i]:
            raise MatlabError(CROP_INPUT, "Upper bound for bbox is smaller than lower bound")
        if hi[i] > dims[i]:
            raise MatlabError(CROP_INPUT, "Upper bound for bbox is invalid")
    src = _mex_class(a, TIFF_TYPE)
    v = _view(src, 4)
    ext = (hi[3] - lo[3], hi[2] - lo[2], hi[0] - lo[0], hi[1] - lo[1])
    res = np.empty(ext, dtype=v.dtype)
    _copy(v, (lo[3], lo[2], lo[0], lo[1]), res, (0, 0, 0, 0), ext, num_threads)
    return res.reshape(_python_shape((ext[2], ext[3], ext[1], ext[0])))


def _same_class(o: np.ndarray, r: np.ndarray) -> None:
    """mxGetClassID of both arrays (a complex array has the class of its real part)"""
    if _real(o).dtype != _real(r).dtype:
        raise MatlabError(
            INDEXING_INPUT, "The data type of the region does not match that of the data!"
        )


def _real(a: np.ndarray) -> np.ndarray:
    """what the mex reads of an array: its real part (mxGetPr)"""
    return a.real if a.dtype in _COMPLEX else a


def _indexing_copy(
    o: np.ndarray,
    lo: Sequence[int],
    r: np.ndarray,
    rlo: Sequence[int],
    ext: Sequence[int],
    n: int,
    num_threads: int | None,
) -> None:
    """the box rlo + ext of r into the box lo + ext of o, both seen as (..., z, y, x) of n
    dimensions; the mex's classes, else its class error. A complex o is left as it is: MATLAB
    hands the mex's separate-complex API a copy of a complex array, so the mex's write is
    lost (oracle cases indexing*_mex.complex_*); a complex r gives its real part."""
    ov, rv = _real(o), _real(r)
    if ov.dtype not in _MEX_DTYPES:
        raise MatlabError(TIFF_TYPE, "Data type not suppported")
    if o.dtype in _COMPLEX:
        return
    o4, r4 = _view(ov, n), _view(rv, n)
    if np.may_share_memory(o4, r4):
        r4 = r4.copy()
    _copy(r4, rlo, o4, lo, ext, num_threads)


def _mex_indexing(
    out: Any,
    bbox: Any,
    region: Any,
    n: int,
    num_threads: int | None,
) -> np.ndarray:
    """indexing3d_mex (n = 3) and indexing4d_mex (n = 4)"""
    if _is_sparse(out) or _is_sparse(region):
        raise MatlabError(INDEXING_INPUT, "The data and the region must be full arrays")
    o = _array(out, "the data")
    r = _array(region, "the region")
    ms, mr = _msize(o), _msize(r)
    if len(ms) > n or len(mr) > n:
        raise MatlabError(
            INDEXING_INPUT, f"The data and the region must have at most {n} dimensions"
        )
    dims = ms + [1] * (n - len(ms))
    b = _mex_bbox(bbox, 2 * n, INDEXING_INPUT)
    for i in range(n):
        if not b[i] >= 1:
            raise MatlabError(INDEXING_INPUT, "Lower bounds must be at least 1")
        if not _floor(b[i + n]) <= dims[i]:
            raise MatlabError(INDEXING_INPUT, "Upper bound is invalid")
        if not _floor(b[i + n]) >= _floor(b[i]) - 1:
            raise MatlabError(INDEXING_INPUT, "Upper bound is smaller than the lower bound minus 1")
    lo = [_mex_u64(float(v)) - 1 for v in b[:n]]  # y, x, z(, t)
    hi = [_mex_u64(float(v)) for v in b[n:]]
    rdims = mr + [1] * (n - len(mr))
    if rdims != [hi[i] - lo[i] for i in range(n)]:
        raise MatlabError(INDEXING_INPUT, "Subregion size does not match the bounding box size")
    _same_class(o, r)
    _writeable(o, "the data")
    order = (3, 2, 0, 1) if n == 4 else (2, 0, 1)
    plo = [lo[i] for i in order]
    ext = [hi[i] - lo[i] for i in order]
    _indexing_copy(o, plo, r, [0] * n, ext, n, num_threads)
    return o


def _mex_indexing_crop(
    out: Any, bbox: Any, region: Any, rbbox: Any, num_threads: int | None
) -> np.ndarray:
    """indexing4d_crop_mex(A, bbox, R[, rbbox])"""
    if _is_sparse(out) or _is_sparse(region):
        raise MatlabError(INDEXING_INPUT, "The data and the region must be full arrays")
    o = _array(out, "the data")
    r = _array(region, "the region")
    ms, mr = _msize(o), _msize(r)
    if len(ms) > 4 or len(mr) > 4:
        raise MatlabError(INDEXING_INPUT, "The data and the region must have at most 4 dimensions")
    dims = ms + [1] * (4 - len(ms))
    b = _mex_bbox(bbox, 8, INDEXING_INPUT)
    for i in range(4):
        if not b[i] >= 1:
            raise MatlabError(INDEXING_INPUT, "Lower bounds must be at least 1")
        if not _floor(b[i + 4]) <= dims[i]:
            raise MatlabError(INDEXING_INPUT, "Upper bound for bbox is invalid")
    lo = [(_mex_u64(float(v)) - 1) & _U64 for v in b[:4]]
    hi = [_mex_u64(float(v)) for v in b[4:]]
    for i in range(4):
        if hi[i] <= lo[i]:
            raise MatlabError(INDEXING_INPUT, "Upper bound for bbox is smaller than lower bound")
        if hi[i] > dims[i]:
            raise MatlabError(INDEXING_INPUT, "Upper bound for bbox is invalid")
    rdims = mr + [1] * (4 - len(mr))
    rlo, rhi = [0, 0, 0, 0], list(rdims)
    if rbbox is not None:
        rb = _mex_bbox(rbbox, 8, INDEXING_INPUT, "region bbox")
        for i in range(4):
            if not rb[i] >= 1:
                raise MatlabError(
                    INDEXING_INPUT, "Lower bounds of the region bbox must be at least 1"
                )
            if not _floor(rb[i + 4]) <= rdims[i]:
                raise MatlabError(INDEXING_INPUT, "Upper bound for region bbox is invalid")
        rlo = [(_mex_u64(float(v)) - 1) & _U64 for v in rb[:4]]
        rhi = [_mex_u64(float(v)) for v in rb[4:]]
    for i in range(4):
        if rhi[i] <= rlo[i]:
            raise MatlabError(
                INDEXING_INPUT, "Upper bound for region bbox is smaller than lower bound"
            )
        if rhi[i] > rdims[i]:
            raise MatlabError(INDEXING_INPUT, "Upper bound for region bbox is invalid")
        if hi[i] - lo[i] != rhi[i] - rlo[i]:
            raise MatlabError(
                INDEXING_INPUT, "The range of region bbox is different from that of bbox"
            )
    _same_class(o, r)
    _writeable(o, "the data")
    order = (3, 2, 0, 1)
    _indexing_copy(
        o,
        [lo[i] for i in order],
        r,
        [rlo[i] for i in order],
        [hi[i] - lo[i] for i in order],
        4,
        num_threads,
    )
    return o


# ---------------------------------------------------------------------------------------------
# MATLAB's code of the wrappers (after the mex's catch)
# ---------------------------------------------------------------------------------------------


def _colon(lo: float, hi: float) -> np.ndarray:
    """``lo : hi`` as an index: MATLAB refuses non-finite ends and rounds the values of a range
    with non-integer ends (after its warning MATLAB:colon:nonIntegerIndex, not given here)"""
    from petakit5d.core.dtypes import mround

    if not (math.isfinite(lo) and math.isfinite(hi)):
        raise MatlabError(
            "MATLAB:colon:nonFiniteEndpoint",
            "Nonfinite endpoints or increment for colon operator in index.",
        )
    if hi < lo:
        return np.zeros(0)
    tol = 2.0 * 2.220446049250313e-16 * max(abs(lo), abs(hi))
    n = math.floor((hi - lo) + tol)
    return mround(lo + np.arange(n + 1, dtype=np.float64))


def _subscripts(bbox: Any, n: int) -> list[np.ndarray]:
    """``bbox(1):bbox(n + 1), ..., bbox(n):bbox(2n)`` (linear indexing into bbox)"""
    try:
        vals = _mx(bbox)[2]
    except TypeError:
        raise MatlabError("MATLAB:badsubscript", "bbox is not a numeric array") from None
    if vals.size < 2 * n:
        raise MatlabError(
            "MATLAB:badsubscript",
            f"Index exceeds the number of array elements. Index must not exceed {vals.size}.",
        )
    return [_colon(float(vals[i]), float(vals[i + n])) for i in range(n)]


def _check_subscript(v: np.ndarray, n: int | None, position: int) -> np.ndarray:
    if v.size and v.min() < 1:
        raise MatlabError(
            "MATLAB:badsubscript",
            f"Index in position {position} is invalid. Array indices must be positive integers "
            "or logical values.",
        )
    if n is not None and v.size and v.max() > n:
        raise MatlabError(
            "MATLAB:badsubscript",
            f"Index in position {position} exceeds array bounds. Index must not exceed {n}.",
        )
    return v.astype(np.intp) - 1


def _msubs_view(a: np.ndarray, n: int) -> tuple[np.ndarray, list[int]]:
    """the array as MATLAB indexes it with n subscripts: (k, ..., y, x) with the last subscript
    over MATLAB's dimension n and the ones after it (k), and MATLAB's sizes per subscript"""
    ms = _msize(a)
    ms = ms + [1] * (n - len(ms))
    lead = math.prod(ms[n - 1 :])
    msz = ms[: n - 1] + [lead]
    # (y, x[, z], k) in MATLAB: ours is (k[, z], y, x)
    pshape = tuple(reversed(msz[2:])) + (msz[0], msz[1])
    return a.reshape(pshape), msz


def _to_python_index(subs: Sequence[np.ndarray]) -> tuple[np.ndarray, ...]:
    """MATLAB's subscripts (y, x, z, t, ...) as an np.ix_ index of the (..., z, y, x) view"""
    order = list(reversed(range(2, len(subs)))) + [0, 1]
    return np.ix_(*[subs[i] for i in order])


def _crop_code(vol: Any, bbox: Any, n: int) -> np.ndarray:
    """``A(bbox(1) : bbox(n + 1), ...)`` (crop3d.m and crop4d.m after the mex's catch)"""
    a = _array(vol, "the data")
    subs = _subscripts(bbox, n)
    v, msz = _msubs_view(a, n)
    idx = [_check_subscript(subs[i], msz[i], i + 1) for i in range(n)]
    res = np.ascontiguousarray(v[_to_python_index(idx)])
    return res.reshape(_python_shape([len(s) for s in idx]))


def _assigned(o: np.ndarray, r: np.ndarray) -> np.ndarray:
    """the values of r converted to o's class, as MATLAB's indexed assignment converts them"""
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import matlab_class

    if r.dtype == o.dtype:
        return r
    if r.dtype.kind == "c" and o.dtype.kind != "c":
        raise MatlabError(
            "MATLAB:assignment:complexToReal", "Complex values cannot be assigned to a real array."
        )
    try:
        return np.asarray(dtypes.cast(r, matlab_class(o.dtype)))
    except ValueError as err:
        if "nologicalnan" not in str(err):
            raise
        raise MatlabError("MATLAB:nologicalnan", "NaN's cannot be converted to logicals.") from None


def _size_text(s: Sequence[int]) -> str:
    return "-by-".join(str(int(v)) for v in s)


def _assign_code(out: Any, region: Any, bbox: Any, rbbox: Any) -> np.ndarray:
    """indexing4d.m after the mex's catch: ``A(bbox(1) : bbox(5), ...) = R`` (with rbbox:
    ``= R(rbbox(1) : rbbox(5), ...)``). The array is written in place unless the subscripts
    reach beyond it: then a grown copy (MATLAB grows the array, the new elements 0)."""
    o = _array(out, "the data")
    r = _array(region, "the region")
    subs = _subscripts(bbox, 4)
    if rbbox is not None:
        r = _crop_code(r, rbbox, 4)
    idx = [_check_subscript(subs[i], None, i + 1) for i in range(4)]
    lens = [len(s) for s in idx]
    # the right side: a scalar fills the box; otherwise its size without singletons must equal
    # the box's without singletons
    if r.size != 1:
        mr = [v for v in _msize(r) if v != 1]
        ml = [v for v in lens if v != 1]
        if mr != ml:
            raise MatlabError(
                "MATLAB:subsassigndimmismatch",
                "Unable to perform assignment because the size of the left side is "
                f"{_size_text(_trim(lens))} and the size of the right side is "
                f"{_size_text(_msize(r))}.",
            )
    vals = _assigned(o, r)
    ms = _msize(o)
    ms = ms + [1] * (4 - len(ms))
    need = [max(ms[i], int(idx[i].max()) + 1 if idx[i].size else 0) for i in range(4)]
    if all(lens) and need != ms:
        if len(_msize(o)) > 4:
            raise MatlabError(
                "MATLAB:indexed:ambiguousGrowth",
                "Attempt to grow array along ambiguous dimension.",
            )
        grown = np.zeros(_python_shape(need), dtype=o.dtype)
        _view(grown, 4)[: ms[3], : ms[2], : ms[0], : ms[1]] = _view(o, 4)
        o = grown
    if not all(lens):
        return o
    v, _ = _msubs_view(o, 4)
    block = vals.reshape(()) if vals.size == 1 else _msubs_ordered(vals, lens)
    v[_to_python_index(idx)] = block
    return o


def _trim(s: Sequence[int]) -> list[int]:
    m = [int(v) for v in s]
    while len(m) > 2 and m[-1] == 1:
        m.pop()
    return m


def _msubs_ordered(vals: np.ndarray, lens: Sequence[int]) -> np.ndarray:
    """the right side's values in MATLAB's order laid out as the (t, z, y, x) block of the
    subscripts' lengths (MATLAB assigns element by element in linear order)"""
    from petakit5d.core.arrays import matlab_order

    flat = np.ascontiguousarray(matlab_order(vals)).reshape(-1)
    # MATLAB's linear order of the block (y fastest, then x, z, t) is the C order of (t, z, x, y)
    blk = flat.reshape(lens[3], lens[2], lens[1], lens[0])
    return np.swapaxes(blk, -1, -2)


# ---------------------------------------------------------------------------------------------
# the MATLAB functions
# ---------------------------------------------------------------------------------------------


def crop3d_matlab(
    vol: np.ndarray, bbox: Any, *, mex_only: bool = False, num_threads: int | None = None
) -> np.ndarray:
    """``crop3d(A, bbox)`` (``crop3d_mex`` with ``mex_only``) with MATLAB's ``bbox`` (``[ymin
    xmin zmin ymax xmax zmax]``, 1-based, inclusive) and a ``(z, y, x)`` array.

    The mex: a full double array with 6 columns, bounds truncated (``(uint64_t)``) after the
    checks ``min >= 1``, ``floor(max) <= size`` and ``floor(max) >= floor(min) - 1`` (empty
    ranges allowed), then the class (logical, uint8, uint16, single, double; complex: the real
    part). Of an array of more than 3 MATLAB dimensions it crops the first volume. The result
    has MATLAB's size (an image for one plane). Where the mex raises, crop3d.m's MATLAB code:
    ``A(bbox(1):bbox(4), bbox(2):bbox(5), bbox(3):bbox(6))`` (the third subscript over z and the
    dimensions after it; non-integer range ends rounded; MATLAB's index errors).
    """
    try:
        return _crop3d_mex(vol, bbox, num_threads)
    except MatlabError:
        if mex_only:
            raise
    return _crop_code(vol, bbox, 3)


def crop4d_matlab(
    vol: np.ndarray, bbox: Any, *, mex_only: bool = False, num_threads: int | None = None
) -> np.ndarray:
    """``crop4d(A, bbox)`` (``crop4d_mex`` with ``mex_only``) with MATLAB's ``bbox`` (``[ymin
    xmin zmin tmin ymax xmax zmax tmax]``) and a ``(t, z, y, x)`` array (of at most 4 MATLAB
    dimensions for the mex).

    The mex: a full double array with 8 columns; ``min >= 1`` and ``floor(max) <= size``, then
    the truncated bounds must give a range of at least one element; class errors are
    ``tiff:dataTypeError``. Where the mex raises, crop4d.m's ``A(bbox(1):bbox(5), ...,
    bbox(4):bbox(8))``.
    """
    try:
        return _crop4d_mex(vol, bbox, num_threads)
    except MatlabError:
        if mex_only:
            raise
    return _crop_code(vol, bbox, 4)


def indexing3d_matlab(
    out: np.ndarray, bbox: Any, region: np.ndarray, *, num_threads: int | None = None
) -> np.ndarray:
    """``indexing3d_mex(A, bbox, R)``: ``R`` written into the box of ``A`` (``out``, a
    ``(z, y, x)`` array of at most 3 MATLAB dimensions) **in place**; returns ``out``.

    The mex's checks in its order (full arrays of at most 3 dimensions, a full double box with
    6 columns, its bounds as ``crop3d_mex``'s, ``R`` of the box's size, the same class, then
    the class). ``out`` must be a writeable C-contiguous array.
    """
    return _mex_indexing(out, bbox, region, 3, num_threads)


def indexing4d_crop_matlab(
    out: np.ndarray,
    bbox: Any,
    region: np.ndarray,
    rbbox: Any = None,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """``indexing4d_crop_mex(A, bbox, R[, rbbox])``: the box ``rbbox`` of ``R`` (default: all
    of it) written into the box ``bbox`` of ``A`` (``out``, ``(t, z, y, x)``) **in place**;
    returns ``out``. Both boxes must hold at least one element per dimension and have the same
    extent."""
    return _mex_indexing_crop(out, bbox, region, rbbox, num_threads)


def indexing4d_matlab(
    out: np.ndarray,
    region: np.ndarray,
    bbox: Any,
    rbbox: Any = None,
    *,
    mex_only: bool = False,
    num_threads: int | None = None,
) -> np.ndarray:
    """``A = indexing4d(A, R, bbox[, rbbox])`` with ``A`` = ``out`` (``(t, z, y, x)``): the
    mex (``indexing4d_mex(A, bbox, R)``, or ``indexing4d_crop_mex(A, bbox, R, rbbox)``) writes
    into ``out`` **in place**, as in MATLAB, and ``out`` is returned.

    Where the mex raises (another class, classes that differ, a box beyond the array, sizes
    that do not match), indexing4d.m's MATLAB code: ``A(bbox(1):bbox(5), ...) = R`` (or
    ``R(rbbox(1):rbbox(5), ...)``), with MATLAB's conversion of ``R`` to ``A``'s class (round
    half away from zero, saturation, NaN to 0; NaN refused for logical), a scalar ``R`` filling
    the box, and MATLAB's errors. That code also writes into ``out`` in place, unless the box
    reaches beyond ``out``: MATLAB grows the array then, and the grown array (a new one, the
    new elements 0) is returned. ``mex_only``: the mex alone.
    """
    try:
        if rbbox is None:
            return _mex_indexing(out, bbox, region, 4, num_threads)
        return _mex_indexing_crop(out, bbox, region, rbbox, num_threads)
    except MatlabError:
        if mex_only:
            raise
    _writeable(out)
    return _assign_code(out, region, bbox, rbbox)
