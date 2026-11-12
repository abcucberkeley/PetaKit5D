"""Array conventions: (y, x, z) volumes in Fortran order, MATLAB's size rules and class names.

A volume is a NumPy array of shape ``(ny, nx, nz)``, Fortran-ordered (y fastest in memory), as
in MATLAB and in the zarr files PetaKit5D writes. A 2D image is ``(ny, nx)``, and counts as a
volume with ``nz == 1``. Linear indices are MATLAB's minus one, in column-major order
(``a.ravel(order="F")``). A 1-D NumPy array stands for a MATLAB row vector (``1 x n``).
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
    "fempty",
    "ffull",
    "fones",
    "fzeros",
    "ind2sub",
    "matlab_class",
    "mndims",
    "mnumel",
    "mpermute",
    "msize",
    "msqueeze",
    "numpy_dtype",
    "sub2ind",
    "yxz_to_zyx",
    "zyx_to_yxz",
]

# MATLAB class name -> NumPy dtype (char is UTF-16 code units in MATLAB; uint16 here)
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


# -- allocation in Fortran order ---------------------------------------------------------------


def fempty(shape: Sequence[int] | int, dtype="single") -> np.ndarray:
    """An uninitialized Fortran-ordered array; dtype may be a MATLAB class name."""
    return np.empty(shape, dtype=numpy_dtype(dtype), order="F")


def fzeros(shape: Sequence[int] | int, dtype="double") -> np.ndarray:
    """MATLAB's ``zeros(shape, class)``, Fortran-ordered."""
    return np.zeros(shape, dtype=numpy_dtype(dtype), order="F")


def fones(shape: Sequence[int] | int, dtype="double") -> np.ndarray:
    return np.ones(shape, dtype=numpy_dtype(dtype), order="F")


def ffull(shape: Sequence[int] | int, value, dtype="double") -> np.ndarray:
    return np.full(shape, value, dtype=numpy_dtype(dtype), order="F")


# -- the volume convention at the public boundary ----------------------------------------------


def as_volume(a, dtype=None, *, copy: bool = False) -> np.ndarray:
    """``a`` as a 3D (y, x, z) Fortran-ordered array (2D gets nz = 1).

    Returns a view when ``a`` already has that layout (and the dtype, if given), else a copy
    in Fortran order. A value conversion, if any, is NumPy's ``astype``: pass integer dtypes
    only for data that are already integral (use ``dtypes.cast`` for MATLAB's conversion).
    """
    a = np.asarray(a)
    if a.ndim == 2:
        a = a[:, :, np.newaxis]
    elif a.ndim != 3:
        raise ValueError(f"a volume has 2 or 3 dimensions (y, x, z), got shape {a.shape}")
    dt = a.dtype if dtype is None else numpy_dtype(dtype)
    if copy or dt != a.dtype or not a.flags.f_contiguous:
        return np.array(a, dtype=dt, order="F", copy=True)
    return a


def check_volume(
    a,
    *,
    name: str = "volume",
    ndim: Iterable[int] = (2, 3),
    dtypes: Iterable | None = None,
    fortran: bool = True,
) -> np.ndarray:
    """Check a volume at a public boundary: number of dimensions, class, Fortran order.

    Raises ``TypeError`` or ``ValueError``; never converts (callers decide whether a copy is
    acceptable, with ``as_volume``).
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
    if fortran and a.ndim > 1 and not a.flags.f_contiguous:
        raise ValueError(f"{name} must be Fortran-ordered (y fastest); use core.arrays.as_volume")
    return a


def zyx_to_yxz(a: np.ndarray) -> np.ndarray:
    """A C-ordered (z, y, x) stack (tiff pages first) as a (y, x, z) Fortran-ordered copy."""
    a = np.asarray(a)
    if a.ndim == 2:
        return np.asfortranarray(a)
    if a.ndim != 3:
        raise ValueError(f"expected (z, y, x), got shape {a.shape}")
    return np.asfortranarray(a.transpose(1, 2, 0))


def yxz_to_zyx(a: np.ndarray) -> np.ndarray:
    """A (y, x, z) volume as a C-ordered (z, y, x) copy (the layout of tiff pages)."""
    a = np.asarray(a)
    if a.ndim == 2:
        return np.ascontiguousarray(a)
    if a.ndim != 3:
        raise ValueError(f"expected (y, x, z), got shape {a.shape}")
    return np.ascontiguousarray(a.transpose(2, 0, 1))


# -- MATLAB's size rules -----------------------------------------------------------------------


def _mshape(a) -> tuple[int, ...]:
    shape = np.shape(a)
    if len(shape) == 0:
        return (1, 1)
    if len(shape) == 1:
        return (1, shape[0])
    # trailing singleton dimensions beyond the second do not exist in MATLAB
    n = len(shape)
    while n > 2 and shape[n - 1] == 1:
        n -= 1
    return tuple(int(s) for s in shape[:n])


def msize(a, dim: int | None = None):
    """MATLAB's ``size``: a tuple of at least 2, without trailing singletons; ``dim`` is 1-based.

    A NumPy scalar is 1 x 1, a 1-D array 1 x n; ``msize(a, 3)`` of a 2D image is 1.
    """
    s = _mshape(a)
    if dim is None:
        return s
    if dim < 1:
        raise ValueError("Dimension argument must be a positive integer scalar")
    return s[dim - 1] if dim <= len(s) else 1


def mndims(a) -> int:
    """MATLAB's ``ndims`` (at least 2)."""
    return len(_mshape(a))


def mnumel(a) -> int:
    return int(np.prod(_mshape(a)))


def msqueeze(a) -> np.ndarray:
    """MATLAB's ``squeeze``: 2D arrays unchanged (1 x n stays a row), else singletons removed.

    The result keeps at least 2 dimensions as MATLAB does (1 x 1 x n gives n x 1), in Fortran
    order (a view when possible).
    """
    a = np.asarray(a)
    s = _mshape(a)
    if len(s) <= 2:
        return a.reshape(s, order="F") if a.ndim != len(s) else a
    kept = [n for n in s if n != 1]
    while len(kept) < 2:
        kept.append(1)
    return a.reshape(kept, order="F")


def sub2ind(size: Sequence[int], *subs) -> np.ndarray | int:
    """MATLAB's ``sub2ind`` (1-based subscripts in, 1-based linear index out, column-major).

    Out-of-range subscripts raise ``IndexError``, as MATLAB errors.
    """
    size = [int(s) for s in size]
    if len(subs) > len(size):
        size = size + [1] * (len(subs) - len(size))
    subs_a = [np.asarray(s) for s in subs]
    # the last given dimension absorbs the rest, as in MATLAB
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
    """MATLAB's ``permute`` with a 1-based ``order``; returns a Fortran-ordered copy."""
    a = np.asarray(a)
    order = [int(o) for o in order]
    if sorted(order) != list(range(1, len(order) + 1)) or len(order) < a.ndim:
        raise ValueError("ORDER must be a permutation of 1:N with N at least ndims(A)")
    if a.ndim < len(order):
        a = a.reshape(a.shape + (1,) * (len(order) - a.ndim), order="A")
    return np.asfortranarray(a.transpose([o - 1 for o in order]))


def axis_order_mapping(input_axis_order: str, output_axis_order: str = "yxz") -> tuple[int, ...]:
    """``axis_order_mapping``: the 1-based positions of ``input_axis_order``'s axes in the output.

    ``mpermute(vol, axis_order_mapping(order, 'yxz'))`` stores yxz data in ``order`` (as
    XR_deskewRotateFrame applies outputAxisOrder). MATLAB's checks are kept: the orders must
    be 3 letters and a substring of one of the six orders (so 'xx' style repeats fail, but the
    check uses ``contains``, as MATLAB does).
    """
    orders = ("xyz", "yxz", "zyx", "zxy", "yzx", "xzy")
    inp = input_axis_order.lower()
    if len(inp) != 3 or not any(inp in o for o in orders):
        raise ValueError("Input axis order must contains x, y, and z.")
    out = output_axis_order
    if len(out) != 3 or not any(out in o for o in orders):
        raise ValueError("Output axis order must contains x, y, and z.")
    # strfind(outputAxisOrder, c): the first position (1-based); MATLAB errors if absent
    return tuple(out.index(c) + 1 for c in inp)
