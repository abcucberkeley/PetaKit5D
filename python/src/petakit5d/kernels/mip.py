"""Max pooling, range maximum (MIPs) and the minimum in a box, on ``(z, y, x)`` volumes.

The Python side of ``petakit5d._ext.mip`` (``python/cpp/mip``), the kernels of MATLAB's
``max_pooling_3d_mex``, ``range_max_3d_mex`` and ``min_bbox_3d_mex``:

- ``max_pooling_3d(vol, pool_size)``: pools ``(pz, py, px)``;
- ``range_max_3d(vol, axis, starts)``: a NumPy axis, 0-based range starts;
- ``min_bbox_3d(vol, box)``: a box ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive.

Arrays are C-contiguous ``(z, y, x)`` volumes or ``(y, x)`` images (other layouts are copied);
results have the input's layout and number of dimensions. Arrays of more than 3 dimensions
are read as the mex reads MATLAB's: their first ``(z, y, x)`` volume. Errors carry MATLAB's
identifier (``MatlabError.identifier``).

``max_pooling_3d_matlab`` and ``min_bbox_3d_matlab`` take the pool size and the box in
MATLAB's values (``[py px pz]``, ``[ymin xmin zmin ymax xmax zmax]``, any numbers) with
``(z, y, x)`` arrays, for the MATLAB-name layer: the mex's checks and truncation, and where
the mex rejects the class or the value, the MATLAB code of ``max_pooling_3d.m`` and
``min_bbox_3d.m`` with its results and errors (``mex_only=True``: the mex alone).

Known issue kept from MATLAB: max pooling starts each block's maximum at 0, so a block of
negative values gives 0 (the mex; MATLAB's code for the other integer classes pads edge blocks
with 0, so its partial blocks do the same).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from types import ModuleType
from typing import Any

import numpy as np

__all__ = [
    "MatlabError",
    "max_pooling_3d",
    "max_pooling_3d_matlab",
    "min_bbox_3d",
    "min_bbox_3d_matlab",
    "range_max_3d",
]

# classes the kernels (and the mex functions) take
_POOL_DTYPES = frozenset(
    np.dtype(t) for t in (np.bool_, np.uint8, np.uint16, np.float32, np.float64)
)
_RANGE_DTYPES = frozenset(np.dtype(t) for t in (np.uint8, np.uint16, np.float32, np.float64))
_INT_DTYPES = frozenset(
    np.dtype(t)
    for t in (np.int8, np.uint8, np.int16, np.uint16, np.int32, np.uint32, np.int64, np.uint64)
)

_U64 = 0xFFFFFFFFFFFFFFFF
_ext_module: ModuleType | None = None


class MatlabError(ValueError):
    """An error MATLAB raises, with its identifier (``'pooling:inputError'``, ...)."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _ext() -> ModuleType:
    global _ext_module
    if _ext_module is None:
        from petakit5d import _ext

        _ext_module = _ext.load("mip")
    return _ext_module


_available: int | None = None


def _threads(num_threads: int | None) -> int:
    """``core.threads.resolve_threads`` with the CPUs counted once per process (counting reads
    the cgroup files, about 0.1 ms: 3% of a 200-million-voxel MIP at 24 threads)"""
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


def _array(x: Any) -> np.ndarray:
    a = np.asarray(x)
    if a.dtype.kind not in "biuf":
        raise MatlabError("mat:dataTypeError", "Data type not suppported")
    return a


def _vol3(a: np.ndarray) -> np.ndarray:
    """the (z, y, x) volume a kernel reads, C-contiguous: an image is (1, ny, nx), a row
    (1, 1, n); of more dimensions the first volume (the mex reads MATLAB's first 3)"""
    if a.ndim < 3:
        a = a.reshape((1,) * (3 - a.ndim) + a.shape)
    elif a.ndim > 3:
        a = a.reshape((-1, *a.shape[-3:]))[0]
    return np.ascontiguousarray(a)


def _shaped(out: np.ndarray, a: np.ndarray) -> np.ndarray:
    """a kernel's (z, y, x) result with the input's number of dimensions (2 or 3)"""
    return out[0] if a.ndim < 3 else out


def _largest(dt: np.dtype) -> np.generic:
    if dt == np.bool_:
        return np.bool_(True)
    if dt.kind == "f":
        return np.finfo(dt).max
    return dt.type(np.iinfo(dt).max)


def _mex_u64(v: float) -> int:
    """``(uint64_t)v`` of a double as x86-64 code compiled by GCC computes it (values beyond
    the range wrap as cvttsd2si gives them)"""
    two63 = 9223372036854775808.0
    if v >= two63:
        # subtract 2^63, convert, flip the top bit (Inf and 2^64 or more give 0)
        w = v - two63
        return (int(w) if w < two63 else 1 << 63) ^ (1 << 63)
    # cvttsd2si: NaN and values below -2^63 give 0x8000000000000000
    if v != v or v < -two63:
        return 1 << 63
    return int(v) & _U64


# ---------------------------------------------------------------------------------------------
# max pooling
# ---------------------------------------------------------------------------------------------


def _check_pool(vals: Sequence[float]) -> tuple[int, int, int]:
    """the mex's check of each pool size (at least 1 after truncation, NaN refused)"""
    for v in vals:
        if not (np.floor(v) >= 1 and v < 9.0e15):
            raise MatlabError("pooling:inputError", "Pooling size must be at least 1")
    return int(vals[0]), int(vals[1]), int(vals[2])


def _pool_blocks(v: np.ndarray, pool: tuple[int, int, int]) -> np.ndarray:
    """MATLAB's code for the integer classes the mex does not take (``padarray`` with 0 to
    whole blocks, then the maximum of each block): whole blocks give their true maximum"""
    m = [-(-n // p) for n, p in zip(v.shape, pool, strict=True)]
    padded = np.zeros([mi * p for mi, p in zip(m, pool, strict=True)], dtype=v.dtype)
    padded[: v.shape[0], : v.shape[1], : v.shape[2]] = v
    blocks = padded.reshape((m[0], pool[0], m[1], pool[1], m[2], pool[2]))
    return np.ascontiguousarray(blocks.max(axis=(1, 3, 5)))


def max_pooling_3d(
    vol: np.ndarray, pool_size: Sequence[float], *, num_threads: int | None = None
) -> np.ndarray:
    """The maximum of each block of ``pool_size`` (``(pz, py, px)``, truncated) of a volume;
    edge blocks are cut at the volume's end, the result is ``ceil(shape / pool_size)``.

    Logical, uint8, uint16, single and double as ``max_pooling_3d_mex``: each maximum starts at
    0 (a block of negative values gives 0; NaN are skipped). Other integer classes as
    ``max_pooling_3d.m``'s MATLAB code: whole blocks give their maximum, partial edge blocks are
    padded with 0. Raises ``MatlabError`` ``pooling:inputError`` for a bad pool size.
    """
    a = _array(vol)
    p = np.asarray(pool_size, dtype=np.float64).ravel()
    if p.size != 3:
        raise MatlabError("pooling:inputError", "Input range for pooling size is not 3")
    pool = _check_pool([float(x) for x in p])
    if a.dtype not in _POOL_DTYPES and a.dtype not in _INT_DTYPES:
        raise MatlabError("mat:dataTypeError", "Data type not suppported")
    v = _vol3(a)
    if a.dtype not in _POOL_DTYPES:
        return _shaped(_pool_blocks(v, pool), a)
    out = np.empty([-(-n // q) for n, q in zip(v.shape, pool, strict=True)], dtype=a.dtype)
    if out.size:
        _ext().max_pooling_3d(v, out, *pool, _threads(num_threads))
    return _shaped(out, a)


def _max_pooling_matlab_code(a: np.ndarray, pool_size: Any) -> np.ndarray:
    """max_pooling_3d.m's MATLAB code on MATLAB's pool values, with the errors MATLAB raises:

    sz = size(orig);
    orig = padarray(orig, ceil(sz ./ poolsize) .* poolsize - sz, 0, 'post');
    sz = size(orig);
    orig = reshape(orig, poolsize(1), sz(1) ./ poolsize(1), poolsize(2), ...);
    out = squeeze(max(orig, [], [1, 3, 5]));
    """
    from petakit5d.core.arrays import from_matlab, size_to_matlab, to_matlab

    p = np.asarray(pool_size, dtype=np.float64)
    if p.ndim == 2 and p.shape[0] != 1 and p.size != 1:
        # a column expands against size(orig) into a matrix: not reproduced
        raise MatlabError("pooling:inputError", "Input range for pooling size is not 3")
    pv = [float(x) for x in p.ravel()]
    sz = size_to_matlab(a.shape)
    if len(pv) not in (1, len(sz)):
        raise MatlabError(
            "MATLAB:sizeDimensionsMustMatch", "Arrays have incompatible sizes for this operation."
        )
    full = pv * len(sz) if len(pv) == 1 else pv
    with np.errstate(all="ignore"):
        q = np.asarray(sz, dtype=np.float64) / np.asarray(full)
        pad = np.ceil(q) * np.asarray(full) - np.asarray(sz, dtype=np.float64)
    if np.isnan(pad).any():
        raise MatlabError(
            "MATLAB:padarray:expectedNonNaN", "Expected input number 2, PADSIZE, to be non-NaN."
        )
    if (pad < 0).any():
        raise MatlabError(
            "MATLAB:padarray:expectedNonnegative",
            "Expected input number 2, PADSIZE, to be nonnegative.",
        )
    if (pad != np.floor(pad)).any() or not np.isfinite(pad).all():
        raise MatlabError(
            "MATLAB:padarray:expectedInteger",
            "Expected input number 2, PADSIZE, to be integer-valued.",
        )
    if len(pv) < 3:
        raise MatlabError(
            "MATLAB:badsubscript",
            f"Index exceeds the number of array elements. Index must not exceed {len(pv)}.",
        )
    psz = [int(s + d) for s, d in zip(sz, pad, strict=True)]
    dims = []
    for i in range(3):
        dims += [pv[i], psz[i] / pv[i]]
    if any(d != np.floor(d) or not np.isfinite(d) for d in dims):
        raise MatlabError(
            "MATLAB:getReshapeDims:notRealInt", "Size arguments must be real integers."
        )
    if any(d < 0 for d in dims):
        raise MatlabError(
            "MATLAB:checkDimCommon:nonnegativeSize", "Size vector elements should be nonnegative."
        )
    if len(sz) != 3:
        raise MatlabError(
            "MATLAB:getReshapeDims:notSameNumel",
            "Number of elements must not change. Use [] as one of the size inputs to "
            "automatically calculate the appropriate size for that dimension.",
        )
    out = _pool_blocks(_vol3(a), (int(pv[2]), int(pv[0]), int(pv[1])))
    # MATLAB's max gives size [1 my 1 mx 1 mz] (trailing ones dropped), then squeeze: a row
    # 1 x my when mx = mz = 1
    mz, my, mx = out.shape
    s = [1, my, 1, mx, 1, mz]
    while len(s) > 2 and s[-1] == 1:
        s.pop()
    if len(s) > 2:
        s = [d for d in s if d != 1]
        s += [1] * (2 - len(s))
    return from_matlab(to_matlab(out).reshape(s, order="F"))


def max_pooling_3d_matlab(
    vol: np.ndarray,
    pool_size: Any,
    *,
    mex_only: bool = False,
    num_threads: int | None = None,
) -> np.ndarray:
    """``max_pooling_3d(orig, poolSize)`` with MATLAB's ``poolSize`` (``[py px pz]``) and a
    ``(z, y, x)`` volume.

    The mex's checks (a 1 x 3 double row, each at least 1 after truncation) and its result for
    logical, uint8, uint16, single and double. Where the mex rejects the pool size or the class,
    ``max_pooling_3d.m`` runs MATLAB's code (MATLAB warns ``PetaKit5D:mexFallback`` first for a
    rejected pool size; not here): its result (squeezed as MATLAB squeezes) or its error.
    ``mex_only``: the mex alone (``pooling:inputError``, ``mat:dataTypeError``).
    """
    a = _array(vol)
    p = np.asarray(pool_size, dtype=np.float64)
    try:
        if p.size != 3 or p.ndim > 2 or (p.ndim == 2 and p.shape[0] != 1):
            raise MatlabError("pooling:inputError", "Input range for pooling size is not 3")
        py, px, pz = _check_pool([float(x) for x in p.ravel()])
        if a.dtype in _POOL_DTYPES:
            return max_pooling_3d(a, (pz, py, px), num_threads=num_threads)
        raise MatlabError("mat:dataTypeError", "Data type not suppported")
    except MatlabError as err:
        if mex_only:
            raise
        if err.identifier == "mat:dataTypeError" and a.dtype not in _INT_DTYPES:
            raise
    return _max_pooling_matlab_code(a, pool_size)


# ---------------------------------------------------------------------------------------------
# range maximum
# ---------------------------------------------------------------------------------------------


def range_max_3d(
    vol: np.ndarray,
    axis: int,
    starts: Sequence[int] | np.ndarray,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """The maximum of ``vol`` along ``axis`` over consecutive ranges that begin at ``starts``
    (0-based, increasing, the first 0; a range ends where the next begins, the last at the end).

    The result has ``vol``'s shape with ``len(starts)`` along ``axis``; range r is
    ``vol.take(range(starts[r], starts[r + 1]), axis).max(axis)`` with MATLAB's rules: NaN are
    ignored (a range of NaN gives its first NaN), and of equal values (0 and -0) the first is
    kept. A MIP is ``range_max_3d(vol, axis, [0])``. uint8, uint16, single and double, at most 3
    dimensions, not empty (``range_max_3d_mex``); an image takes axis -3 as a z of size 1.
    ``MatlabError`` ``range_max_3d:inputError`` otherwise.
    """
    ident = "range_max_3d:inputError"
    a = np.asarray(vol)
    if a.dtype not in _RANGE_DTYPES:
        raise MatlabError(ident, "im must be a real uint8, uint16, single or double array.")
    if a.ndim > 3:
        raise MatlabError(ident, "im must have at most 3 dimensions.")
    v = _vol3(a)
    if v.size == 0:
        raise MatlabError(ident, "im must not be empty.")
    if isinstance(axis, (bool, np.bool_)) or not isinstance(axis, (int, np.integer)):
        raise MatlabError(ident, "axis must be an axis of im.")
    neg = int(axis) - a.ndim if axis >= 0 else int(axis)
    if not -3 <= neg <= -1 or (axis >= 0 and axis >= a.ndim):
        raise MatlabError(ident, "axis must be an axis of im.")
    kax = 3 + neg  # 0: z, 1: y, 2: x
    st = np.asarray(starts)
    if st.size == 0 or st.dtype.kind not in "iuf":
        raise MatlabError(ident, "starts must be a nonempty vector of numbers.")
    n = v.shape[kax]
    vals = [float(x) for x in st.astype(np.float64).ravel()]
    prev = -1.0
    for r, x in enumerate(vals):
        if x <= prev or x >= float(n) or x != math.floor(x) or (r == 0 and x != 0):
            raise MatlabError(
                ident, "starts must be increasing integers from 0 to below the size along axis."
            )
        prev = x
    begin = [int(x) for x in vals]
    shape = list(v.shape)
    shape[kax] = len(begin)
    out = np.empty(shape, dtype=a.dtype)
    _ext().range_max_3d(v, out, kax, begin, _threads(num_threads))
    return _shaped(out, a)


# ---------------------------------------------------------------------------------------------
# minimum in a box
# ---------------------------------------------------------------------------------------------


def _min_of(v: np.ndarray, box: Sequence[int], num_threads: int | None) -> np.generic:
    """the minimum of v[z0:z1, y0:y1, x0:x1] (a valid box) from the largest value of the class"""
    if any(box[i + 3] == box[i] for i in range(3)):
        return _largest(v.dtype)
    if v.dtype in _POOL_DTYPES:
        return v.dtype.type(_ext().min_box_3d(v, *box, _threads(num_threads)))
    region = v[box[0] : box[3], box[1] : box[4], box[2] : box[5]]
    return v.dtype.type(min(region.min(), _largest(v.dtype)))


def min_bbox_3d(
    vol: np.ndarray, box: Sequence[int | None], *, num_threads: int | None = None
) -> np.generic:
    """The minimum of ``vol`` in ``box`` (``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive;
    an image also takes ``(y0, x0, y1, x1)``; a ``None`` end is the volume's end), starting from
    the largest value of the class as ``min_bbox_3d_mex``: NaN are skipped, and an empty box or
    a box of NaN gives that value. Of zeros (0 and -0) it is the first in MATLAB's order.

    Logical and every integer and float class. ``MatlabError`` ``minbbox:inputError`` for a box
    that is not integral, reversed or outside the volume.
    """
    a = _array(vol)
    v = _vol3(a)
    vals = list(box)
    if len(vals) == 4 and a.ndim == 2:
        vals = [0, vals[0], vals[1], 1, vals[2], vals[3]]
    if len(vals) != 6:
        raise MatlabError("minbbox:inputError", "The length of bbox is not 6")
    b = []
    for i, x in enumerate(vals):
        if x is None and i >= 3:
            x = v.shape[i - 3]
        if x is None or isinstance(x, (bool, np.bool_)) or float(x) != math.floor(float(x)):
            raise MatlabError("minbbox:inputError", "Bbox is not valid")
        b.append(int(x))
    if any(b[i] < 0 or b[i + 3] < b[i] or b[i + 3] > v.shape[i] for i in range(3)):
        raise MatlabError("minbbox:inputError", "Bbox is not valid")
    return _min_of(v, b, num_threads)


def _colon_values(lo: float, hi: float) -> np.ndarray:
    """``lo : hi`` in an index: MATLAB refuses non-finite ends, and rounds the values of a
    range with non-integer ends (warning MATLAB:colon:nonIntegerIndex; not given here)"""
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


def _check_index(v: np.ndarray, n: int, position: int) -> np.ndarray:
    if v.size and v.min() < 1:
        raise MatlabError(
            "MATLAB:badsubscript",
            f"Index in position {position} is invalid. Array indices must be positive integers "
            "or logical values.",
        )
    if v.size and v.max() > n:
        raise MatlabError(
            "MATLAB:badsubscript",
            f"Index in position {position} exceeds array bounds. Index must not exceed {n}.",
        )
    return v.astype(np.intp) - 1


def _min_bbox_matlab_code(a: np.ndarray, bbox: Any, num_threads: int | None) -> np.generic:
    """min_bbox_3d.m's MATLAB code: ``crop3d`` (``A(b(1):b(4), b(2):b(5), b(3):b(6))``, the
    third subscript over z and the dimensions after it) and the region's minimum, starting
    from the largest value of the class"""
    b = np.asarray(bbox, dtype=np.float64).ravel(order="F")
    if b.size < 6:
        raise MatlabError(
            "MATLAB:badsubscript",
            f"Index exceeds the number of array elements. Index must not exceed {b.size}.",
        )
    if a.ndim < 3:
        a = a.reshape((1,) * (3 - a.ndim) + a.shape)
    v = a.reshape((-1, *a.shape[-2:]))  # (k, y, x): MATLAB's third subscript is k
    sizes = (v.shape[1], v.shape[2], v.shape[0])  # MATLAB's subscripts y, x, k
    subs = [_colon_values(float(b[i]), float(b[i + 3])) for i in range(3)]
    iy, ix, ik = (_check_index(subs[i], sizes[i], i + 1) for i in range(3))
    region = np.ascontiguousarray(v[np.ix_(ik, iy, ix)])
    if region.size == 0:
        return _largest(a.dtype)
    return _min_of(region, [0, 0, 0, *region.shape], num_threads)


def min_bbox_3d_matlab(
    vol: np.ndarray,
    bbox: Any,
    *,
    mex_only: bool = False,
    num_threads: int | None = None,
) -> np.generic:
    """``min_bbox_3d(orig, bbox)`` with MATLAB's ``bbox`` (``[ymin xmin zmin ymax xmax zmax]``,
    1-based, inclusive) and a ``(z, y, x)`` volume.

    The mex's checks (a 1 x 6 row; values truncated as ``(uint64_t)``; ``max >= min - 1``, within
    the volume) and its result for logical, uint8, uint16, single and double. Where the mex
    rejects the class or the box, ``min_bbox_3d.m`` runs MATLAB's code (MATLAB warns
    ``PetaKit5D:mexFallback`` first for a rejected box; not here): the box as MATLAB indexing
    reads ``b(1):b(4)`` and so on. A reversed range is empty, the values of a range with
    non-integer ends are rounded (``0.6:5`` is 1 to 5), any 6 or more elements in any shape are
    taken, and ranges outside the array or with non-finite ends raise MATLAB's errors.
    ``mex_only``: the mex alone (``minbbox:inputError``, ``mat:dataTypeError``).
    """
    a = _array(vol)
    b = np.asarray(bbox, dtype=np.float64)
    try:
        if b.size != 6 or b.ndim > 2 or (b.ndim == 2 and b.shape[0] != 1):
            raise MatlabError("minbbox:inputError", "The length of bbox is not 6")
        v = _vol3(a)
        u = [_mex_u64(float(x)) for x in b.ravel()]
        lo = [(x - 1) & _U64 for x in u[:3]]  # y, x, z
        hi = u[3:]
        msz = (v.shape[1], v.shape[2], v.shape[0])
        if any(h < s for h, s in zip(hi, lo, strict=True)) or any(
            h > n for h, n in zip(hi, msz, strict=True)
        ):
            raise MatlabError("minbbox:inputError", "Bbox is not valid")
        if a.dtype not in _POOL_DTYPES:
            raise MatlabError("mat:dataTypeError", "Data type not suppported")
        return _min_of(v, [lo[2], lo[0], lo[1], hi[2], hi[0], hi[1]], num_threads)
    except MatlabError:
        if mex_only:
            raise
    return _min_bbox_matlab_code(a, bbox, num_threads)
