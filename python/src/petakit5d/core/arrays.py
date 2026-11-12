"""Array conventions: (z, y, x) volumes in C order, and the mapping to and from MATLAB's arrays.

A volume is a NumPy array ``(z, y, x)``, C-contiguous (x fastest, as in a tiff page); an image is
``(y, x)``; extra dimensions lead (``(tile, z, y, x)``). MATLAB's array of the same data is
``(y, x, z, d4, ...)`` in Fortran order, which is NOT the same bytes: every z plane is
transposed. ``from_matlab`` and ``to_matlab`` map image arrays of any number of dimensions;
matrices and vectors of parameters keep MATLAB's indices (``M(i, j)`` is ``P[i - 1, j - 1]``,
a 1-D array is a row) and are not mapped.

MATLAB's linear order (y fastest, then x, then z) is the C order of ``matlab_order(vol)``, the
view with axes ``(..., z, x, y)``: the voxel ``(z, y, x)`` has MATLAB's 0-based linear index
``y + ny*(x + nx*z)``.

The functions ``msize``, ``mndims``, ``mnumel``, ``msqueeze``, ``sub2ind``, ``ind2sub`` and
``mpermute`` are MATLAB's rules on arrays in MATLAB's form (parameter values, ``.mat``
contents, the MATLAB-name layer).
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

__all__ = [
    "CLASS_DTYPES",
    "as_volume",
    "axis_order_mapping",
    "check_volume",
    "data_type_to_byte_number",
    "from_matlab",
    "ind2sub",
    "matlab_class",
    "matlab_index_to_zyx",
    "matlab_order",
    "mndims",
    "mnumel",
    "mpermute",
    "msize",
    "msqueeze",
    "numpy_dtype",
    "size_from_matlab",
    "size_to_matlab",
    "sub2ind",
    "to_matlab",
    "zyx_to_matlab_index",
]

# MATLAB class name -> NumPy dtype
CLASS_DTYPES: dict[str, np.dtype] = {
    "double": np.dtype(np.float64),
    "single": np.dtype(np.float32),
    "int8": np.dtype(np.int8),
    "uint8": np.dtype(np.uint8),
    "int16": np.dtype(np.int16),
    "uint16": np.dtype(np.uint16),
    "int32": np.dtype(np.int32),
    "uint32": np.dtype(np.uint32),
    "int64": np.dtype(np.int64),
    "uint64": np.dtype(np.uint64),
    "logical": np.dtype(np.bool_),
}
_DTYPE_CLASSES = {v: k for k, v in CLASS_DTYPES.items()}

# above this size from_matlab and to_matlab copy in several threads (NumPy's copy releases the GIL)
_THREADED_BYTES = 64 << 20


def numpy_dtype(cls: str | np.dtype | type) -> np.dtype:
    """The NumPy dtype of a MATLAB class name ('uint16', 'single', ...); dtypes pass through."""
    if isinstance(cls, str) and cls in CLASS_DTYPES:
        return CLASS_DTYPES[cls]
    dt = np.dtype(cls)
    if dt not in _DTYPE_CLASSES:
        raise TypeError(f"no MATLAB class for dtype {dt}")
    return dt


def matlab_class(x) -> str:
    """MATLAB's ``class`` of an array, a dtype or a Python scalar (float: 'double')."""
    if isinstance(x, bool | np.bool_):
        return "logical"
    if isinstance(x, int | float) and not isinstance(x, np.generic):
        return "double"
    if isinstance(x, str) and x in CLASS_DTYPES:
        return x
    dt = x.dtype if hasattr(x, "dtype") else np.dtype(x)
    try:
        return _DTYPE_CLASSES[np.dtype(dt)]
    except KeyError:
        raise TypeError(f"no MATLAB class for dtype {dt}") from None


def data_type_to_byte_number(dtype: str) -> int:
    """``dataTypeToByteNumber``: bytes per element of 'uint8', 'uint16', 'single' or 'double'."""
    table = {"uint8": 1, "uint16": 2, "single": 4, "double": 8}
    if not isinstance(dtype, str) or dtype not in table:
        raise ValueError(f"Unsupported data type {dtype}")
    return table[dtype]


# -- the mapping between MATLAB's arrays and ours ----------------------------------------------


def _from_matlab_axes(n: int) -> tuple[int, ...]:
    # MATLAB (y, x, z, d4, ..., dn) -> (dn, ..., d4, z, y, x)
    return (*range(n - 1, 1, -1), 0, 1) if n >= 2 else tuple(range(n))


def _to_matlab_axes(n: int) -> tuple[int, ...]:
    # (dn, ..., d4, z, y, x) -> MATLAB (y, x, z, d4, ..., dn)
    return (n - 2, n - 1, *range(n - 3, -1, -1)) if n >= 2 else tuple(range(n))


def _copy_threaded(dst: np.ndarray, src: np.ndarray, num_threads: int | None) -> None:
    """``dst[...] = src`` split over the first axis of ``dst`` in threads (bit-exact copy)."""
    from petakit5d.core.threads import resolve_threads

    n = resolve_threads(num_threads)
    if dst.nbytes < _THREADED_BYTES or n == 1 or dst.ndim == 0 or dst.shape[0] < 2:
        dst[...] = src
        return
    from concurrent.futures import ThreadPoolExecutor

    edges = np.linspace(0, dst.shape[0], min(n, dst.shape[0]) * 4 + 1).astype(int)
    parts = [(a, b) for a, b in zip(edges[:-1], edges[1:], strict=True) if b > a]

    def copy(ab):
        dst[ab[0] : ab[1]] = src[ab[0] : ab[1]]

    with ThreadPoolExecutor(max_workers=min(n, len(parts))) as ex:
        list(ex.map(copy, parts))


def from_matlab(a, *, num_threads: int | None = None) -> np.ndarray:
    """MATLAB's image array ``(y, x, z, d4, ...)`` as ours: ``(..., d4, z, y, x)``, C-contiguous.

    ``np.ascontiguousarray(np.transpose(a, (n-1, ..., 2, 0, 1)))``: a 2D image keeps its shape
    (``(y, x)``), a 1-D array is returned as it is (a row). Large arrays are copied in threads.
    """
    a = np.asarray(a)
    if a.ndim <= 1:
        return np.ascontiguousarray(a)
    t = np.transpose(a, _from_matlab_axes(a.ndim))
    if t.flags.c_contiguous:
        return t
    out = np.empty(t.shape, dtype=t.dtype)
    _copy_threaded(out, t, num_threads)
    return out


def to_matlab(a, *, num_threads: int | None = None) -> np.ndarray:
    """Our image array ``(..., d4, z, y, x)`` as MATLAB's ``(y, x, z, d4, ...)`` (Fortran order).

    The inverse of ``from_matlab`` (``np.transpose(a, (n-2, n-1, n-3, ..., 0))``).
    """
    a = np.asarray(a)
    if a.ndim <= 1:
        return np.asfortranarray(a)
    t = np.transpose(a, _to_matlab_axes(a.ndim))
    if t.flags.f_contiguous:
        return t
    # fill a Fortran array through its C-ordered transpose (the axes of `a` reversed)
    out = np.empty(t.shape[::-1], dtype=t.dtype)
    _copy_threaded(out, t.T, num_threads)
    return out.T


def matlab_order(a: np.ndarray) -> np.ndarray:
    """The view ``(..., z, x, y)`` (last two axes swapped) whose C order is MATLAB's linear order.

    ``matlab_order(vol).ravel()`` is MATLAB's ``vol(:)``; "first in MATLAB's order" means first
    in this view.
    """
    a = np.asarray(a)
    return np.swapaxes(a, -1, -2) if a.ndim >= 2 else a


def size_to_matlab(shape: Sequence[int]) -> tuple[int, ...]:
    """MATLAB's size of the array of ``shape`` (``(nz, ny, nx)`` -> ``(ny, nx, nz)``).

    MATLAB's size rules apply: at least 2 values, no trailing singletons beyond the second (an
    image or a volume with nz = 1 is ``(ny, nx)``); a 1-D shape is a row ``(1, n)``.
    """
    s = [int(v) for v in shape]
    if len(s) == 0:
        return (1, 1)
    if len(s) == 1:
        return (1, s[0])
    m = [s[i] for i in _to_matlab_axes(len(s))]
    while len(m) > 2 and m[-1] == 1:
        m.pop()
    return tuple(m)


def size_from_matlab(msize: Sequence[int], ndim: int | None = None) -> tuple[int, ...]:
    """Our shape of a MATLAB array of size ``msize`` (``(ny, nx, nz)`` -> ``(nz, ny, nx)``).

    ``ndim`` asks for that many dimensions (leading singletons added: ``size_from_matlab((ny,
    nx), 3)`` is ``(1, ny, nx)``); by default the shape has as many as ``msize``.
    """
    m = [int(v) for v in msize]
    if len(m) == 1:
        m.append(1)
    if ndim is not None:
        if ndim < 2:
            raise ValueError("ndim must be at least 2")
        trailing = m[ndim:]
        if any(v != 1 for v in trailing):
            raise ValueError(f"MATLAB size {list(msize)} has more than {ndim} dimensions")
        m = m[:ndim] + [1] * (ndim - len(m))
    return tuple(m[i] for i in _from_matlab_axes(len(m)))


def zyx_to_matlab_index(idx, shape: Sequence[int]):
    """MATLAB's 0-based linear index of the voxel ``idx`` (a tuple of index arrays or ints, in the
    axes of ``shape``): ``y + ny*(x + nx*z)`` for ``(z, y, x)`` (MATLAB's 1-based index is this
    plus 1)."""
    shape = tuple(int(v) for v in shape)
    idx = tuple(np.asarray(i) for i in idx)
    if len(idx) != len(shape):
        raise ValueError(f"{len(idx)} indices for a shape of {len(shape)} dimensions")
    for i, n in zip(idx, shape, strict=True):
        if np.any((i < 0) | (i >= n)):
            raise IndexError("index out of range")
    m = np.ravel_multi_index(tuple(idx[k] for k in _swap_last(len(shape))), _swapped(shape))
    return int(m) if np.ndim(m) == 0 else m


def matlab_index_to_zyx(k, shape: Sequence[int]) -> tuple:
    """The voxel (indices in the axes of ``shape``, e.g. ``(z, y, x)``) of MATLAB's 0-based linear
    index ``k``; the inverse of ``zyx_to_matlab_index``."""
    shape = tuple(int(v) for v in shape)
    sw = np.unravel_index(np.asarray(k), _swapped(shape))
    out = [sw[j] for j in _swap_last(len(shape))]
    return tuple(int(v) if np.ndim(v) == 0 else v for v in out)


def _swap_last(n: int) -> list[int]:
    p = list(range(n))
    if n >= 2:
        p[-1], p[-2] = p[-2], p[-1]
    return p


def _swapped(shape: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(shape[i] for i in _swap_last(len(shape)))


# -- the volume convention at the public boundary ----------------------------------------------


def as_volume(a, dtype=None, *, copy: bool = False) -> np.ndarray:
    """``a`` as a C-contiguous ``(z, y, x)`` array (an image ``(y, x)`` gets ``nz = 1``).

    A view when ``a`` already has that layout (and the dtype, if given), else a C-ordered copy. A
    value conversion, if any, is NumPy's ``astype``: pass integer dtypes only for data that are
    already integral (``dtypes.cast`` is MATLAB's conversion).
    """
    a = np.asarray(a)
    if a.ndim == 2:
        a = a[np.newaxis]
    elif a.ndim != 3:
        raise ValueError(f"a volume has 2 or 3 dimensions (z, y, x), got shape {a.shape}")
    dt = a.dtype if dtype is None else numpy_dtype(dtype)
    if copy or dt != a.dtype or not a.flags.c_contiguous:
        return np.array(a, dtype=dt, order="C", copy=True)
    return a


def check_volume(
    a,
    *,
    name: str = "volume",
    ndim: Iterable[int] = (2, 3),
    dtypes: Iterable | None = None,
    contiguous: bool = True,
) -> np.ndarray:
    """Check a volume at a public boundary: dimensions, class and C order (never converts).

    Raises ``TypeError`` or ``ValueError``; callers decide whether a copy is acceptable
    (``as_volume``). A Fortran-ordered array is refused with a hint: it is likely MATLAB's layout.
    """
    if not isinstance(a, np.ndarray):
        raise TypeError(f"{name} must be a NumPy array, got {type(a).__name__}")
    ndim = tuple(ndim)
    if a.ndim not in ndim:
        raise ValueError(f"{name} must have {' or '.join(map(str, ndim))} dimensions, got {a.ndim}")
    if dtypes is not None:
        allowed = tuple(numpy_dtype(d) for d in dtypes)
        if a.dtype not in allowed:
            names = ", ".join(matlab_class(d) for d in allowed)
            raise TypeError(f"{name} must be of class {names}, got {matlab_class(a)}")
    if contiguous and not a.flags.c_contiguous:
        hint = (
            " (it is Fortran-ordered: MATLAB's layout? use core.arrays.from_matlab)"
            if a.flags.f_contiguous
            else ""
        )
        raise ValueError(f"{name} must be C-contiguous (z, y, x){hint}")
    return a


# -- MATLAB's size rules (on arrays in MATLAB's form) ------------------------------------------


def _mshape(a) -> tuple[int, ...]:
    shape = np.shape(a)
    if len(shape) == 0:
        return (1, 1)
    if len(shape) == 1:
        return (1, shape[0])
    n = len(shape)
    while n > 2 and shape[n - 1] == 1:
        n -= 1
    return tuple(int(s) for s in shape[:n])


def msize(a, dim: int | None = None):
    """MATLAB's ``size`` of an array in MATLAB's form: at least 2 values, no trailing singletons;
    ``dim`` is 1-based. A NumPy scalar is 1 x 1, a 1-D array 1 x n. (Our shape's MATLAB size is
    ``size_to_matlab(a.shape)``.)"""
    s = _mshape(a)
    if dim is None:
        return s
    if dim < 1:
        raise ValueError("Dimension argument must be a positive integer scalar")
    return s[dim - 1] if dim <= len(s) else 1


def mndims(a) -> int:
    """MATLAB's ``ndims`` (at least 2) of an array in MATLAB's form."""
    return len(_mshape(a))


def mnumel(a) -> int:
    return int(np.prod(_mshape(a)))


def msqueeze(a) -> np.ndarray:
    """MATLAB's ``squeeze`` of an array in MATLAB's form: 2D unchanged, else singletons removed,
    at least 2 dimensions kept (1 x 1 x n gives n x 1); element order is MATLAB's linear order."""
    a = np.asarray(a)
    s = _mshape(a)
    if len(s) <= 2:
        return a.reshape(s, order="F") if a.ndim != len(s) else a
    kept = [n for n in s if n != 1]
    while len(kept) < 2:
        kept.append(1)
    return a.reshape(kept, order="F")


def sub2ind(size: Sequence[int], *subs) -> np.ndarray | int:
    """MATLAB's ``sub2ind`` (1-based subscripts in, 1-based linear index out, column-major), for
    MATLAB-valued sizes and subscripts. Out-of-range subscripts raise ``IndexError``."""
    size = [int(s) for s in size]
    if len(subs) > len(size):
        size = size + [1] * (len(subs) - len(size))
    subs_a = [np.asarray(s) for s in subs]
    dims = size[: len(subs)]
    dims[-1] = int(np.prod(size[len(subs) - 1 :]))
    ind = np.zeros(np.broadcast(*subs_a).shape, dtype=np.int64)
    stride = 1
    for s, n in zip(subs_a, dims, strict=True):
        if np.any((s < 1) | (s > n)) or np.any(s != np.floor(s)):
            raise IndexError("Out of range subscript.")
        ind = ind + (s.astype(np.int64) - 1) * stride
        stride *= n
    ind = ind + 1
    return int(ind) if ind.ndim == 0 else ind


def ind2sub(size: Sequence[int], ind, nout: int | None = None) -> tuple:
    """MATLAB's ``ind2sub``: 1-based subscripts of 1-based linear indices (column-major)."""
    size = [int(s) for s in size]
    nout = len(size) if nout is None else nout
    dims = size[:nout]
    if nout < len(size):
        dims[-1] = int(np.prod(size[nout - 1 :]))
    elif nout > len(size):
        dims = dims + [1] * (nout - len(size))
    ind = np.asarray(ind)
    if np.any((ind < 1) | (ind > int(np.prod(dims)))):
        raise IndexError("Index exceeds number of array elements.")
    rest = ind.astype(np.int64) - 1
    out = []
    for n in dims[:-1]:
        out.append(rest % n + 1)
        rest = rest // n
    out.append(rest + 1)
    return tuple(int(o) if o.ndim == 0 else o for o in out)


def mpermute(a, order: Sequence[int]) -> np.ndarray:
    """MATLAB's ``permute`` with a 1-based ``order`` on an array in MATLAB's form (a
    Fortran-ordered copy, MATLAB's layout of the result)."""
    a = np.asarray(a)
    order = [int(o) for o in order]
    if sorted(order) != list(range(1, len(order) + 1)) or len(order) < a.ndim:
        raise ValueError("ORDER must be a permutation of 1:N with N at least ndims(A)")
    if a.ndim < len(order):
        a = a.reshape(a.shape + (1,) * (len(order) - a.ndim), order="A")
    return np.asfortranarray(a.transpose([o - 1 for o in order]))


def axis_order_mapping(input_axis_order: str, output_axis_order: str = "yxz") -> tuple[int, ...]:
    """``axis_order_mapping`` on MATLAB's axis strings (the MATLAB-name layer): the 1-based
    positions of ``input_axis_order``'s axes in the output. ``mpermute(vol, axis_order_mapping(
    order, 'yxz'))`` stores MATLAB's yxz data in ``order``. For our arrays and strings
    (``'zyx'``) use ``core.axes.axis_permutation``."""
    orders = ("xyz", "yxz", "zyx", "zxy", "yzx", "xzy")
    inp = input_axis_order.lower()
    if len(inp) != 3 or not any(inp in o for o in orders):
        raise ValueError("Input axis order must contains x, y, and z.")
    out = output_axis_order
    if len(out) != 3 or not any(out in o for o in orders):
        raise ValueError("Output axis order must contains x, y, and z.")
    return tuple(out.index(c) + 1 for c in inp)
