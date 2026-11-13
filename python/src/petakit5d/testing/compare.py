"""Comparison library of the oracle harness: MATLAB values in Python, value files, diffs.

MATLAB values are held in small classes that keep MATLAB's class and size:

- ``MArray``: numeric, logical and char arrays. ``data`` has MATLAB's size (two or more
  dimensions, no trailing singleton dimensions after the second); char data are UTF-16 code
  units (uint16); complex double and single arrays have a complex dtype.
- ``MCell``, ``MStruct`` (field order kept, compared without regard to it), ``MString``
  (string arrays, ``None`` for missing), ``MSparse``, ``MFunc`` (function handles, by text),
  ``MMap`` (containers.Map), ``MObject`` (other objects: class name, the value MATLAB gave for
  it, display text).

``to_mvalue`` converts Python values (guided by the MATLAB reference when one is given),
``read_pkv``/``write_pkv`` read and write the value files exchanged with MATLAB, ``vhash`` and
``describe`` give the SHA-256 and the JSON descriptor kept in the repo manifests, ``compare``
and ``compare_trees`` compare and explain. Formats and rules: ``python/docs/testing.md``.
"""

from __future__ import annotations

import base64
import contextlib
import contextvars
import dataclasses
import fnmatch
import hashlib
import itertools
import json
import math
import os
import re
import struct
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

HASH_VERSION = b"pk1"
INLINE_BYTES = 2048  # arrays up to this size are kept inline in descriptors
INLINE_TEXT = 4096  # text and json files up to this size are kept inline
DESCRIPTOR_BYTES = 8192  # a container whose descriptor would be larger keeps only its hash
MAX_PROBLEMS = 40  # problems listed in a report (all are counted)
HINT_NUMEL = 20_000_000  # arrays up to this size get layout hints

_DTYPE = {
    "double": np.dtype("<f8"),
    "single": np.dtype("<f4"),
    "int8": np.dtype("i1"),
    "uint8": np.dtype("u1"),
    "int16": np.dtype("<i2"),
    "uint16": np.dtype("<u2"),
    "int32": np.dtype("<i4"),
    "uint32": np.dtype("<u4"),
    "int64": np.dtype("<i8"),
    "uint64": np.dtype("<u8"),
    "logical": np.dtype("?"),
    "char": np.dtype("<u2"),
}
_CDTYPE = {"double": np.dtype("<c16"), "single": np.dtype("<c8")}
NUMERIC_CLASSES = (
    "double",
    "single",
    "int8",
    "uint8",
    "int16",
    "uint16",
    "int32",
    "uint32",
    "int64",
    "uint64",
)
INT_CLASSES = NUMERIC_CLASSES[2:]
_CLASS_OF_DTYPE = {np.dtype(v).newbyteorder("="): k for k, v in _DTYPE.items() if k != "char"}
_CLASS_OF_DTYPE.update({np.dtype("c16"): "double", np.dtype("c8"): "single"})
_UINT = {1: np.uint8, 2: np.uint16, 4: np.uint32, 8: np.uint64}
AXES = ("y", "x", "z", "t")

MAGIC = b"PKV1\x00\x00\x00\x00"
END_MAGIC = b"PKV1END\x00"


def msize(shape: Iterable[int]) -> tuple[int, ...]:
    """MATLAB size: at least two dimensions, no trailing singletons after the second"""
    s = [int(n) for n in shape]
    if len(s) < 2:
        raise ValueError(f"a MATLAB size has two or more dimensions, got {tuple(s)}")
    while len(s) > 2 and s[-1] == 1:
        s.pop()
    return tuple(s)


def size_text(size: Sequence[int]) -> str:
    return "[" + " ".join(str(n) for n in size) + "]"


def subscripts(index: int, size: Sequence[int]) -> tuple[int, ...]:
    """MATLAB subscripts (1-based) of a 0-based column-major index"""
    return tuple(int(i) + 1 for i in np.unravel_index(int(index), tuple(size), order="F"))


# what the value being compared is: "vol" for an image array that Python holds as (z, y, x)
# (an output with out role 'vol', a zarr in Python's layout); set with ``layout(kind)``
_LAYOUT: contextvars.ContextVar[str | None] = contextvars.ContextVar("pk5d_layout", default=None)


@contextlib.contextmanager
def layout(kind: str | None) -> Iterator[None]:
    """compare values of this kind ("vol": diff reports give (z, y, x) indices too)"""
    tok = _LAYOUT.set(kind)
    try:
        yield
    finally:
        _LAYOUT.reset(tok)


def _perm_from_matlab(n: int) -> tuple[int, ...]:
    # MATLAB (y, x, z, d4, ..., dn) -> ours (dn, ..., d4, z, y, x), as core.arrays.from_matlab
    return (*range(n - 1, 1, -1), 0, 1) if n >= 2 else tuple(range(n))


def _perm_to_matlab(n: int) -> tuple[int, ...]:
    # ours (dn, ..., d4, z, y, x) -> MATLAB (y, x, z, d4, ..., dn), as core.arrays.to_matlab
    return (n - 2, n - 1, *range(n - 3, -1, -1)) if n >= 2 else tuple(range(n))


def zyx_text(sub: Sequence[int]) -> str:
    """the 0-based index in Python's form of MATLAB's 1-based subscripts ``sub`` (an array of
    3 or more dimensions, or any image while comparing a volume), else ''"""
    n = len(sub)
    if n < 2 or (n == 2 and _LAYOUT.get() != "vol"):
        return ""
    p = [int(v) - 1 for v in sub]
    if n == 2:
        return f"(y, x) ({p[0]}, {p[1]})"
    names = [f"d{k + 1}" for k in range(n - 1, 2, -1)] + ["z", "y", "x"]
    vals = p[2:][::-1] + [p[0], p[1]]
    return f"({', '.join(names)}) ({', '.join(map(str, vals))})"


def _at(i: int, size: Sequence[int]) -> str:
    """MATLAB's subscripts of a 0-based column-major index, with Python's index where it differs"""
    sub = subscripts(i, size)
    z = zyx_text(sub)
    return f"({','.join(map(str, sub))})" + (f" = {z}" if z else "")


def _sha(*parts: bytes | str) -> hashlib._Hash:
    h = hashlib.sha256(HASH_VERSION)
    for p in parts:
        h.update(b":")
        h.update(p.encode("utf-8") if isinstance(p, str) else p)
    return h


def _size_bytes(size: Sequence[int]) -> bytes:
    return ",".join(str(n) for n in size).encode()


# ---------------------------------------------------------------------------------------------
# value classes
# ---------------------------------------------------------------------------------------------


class MValue:
    """base of the MATLAB value classes"""

    kind = ""
    __slots__ = ()

    def hash(self) -> str:
        return vhash(self)


class MArray(MValue):
    """numeric, logical or char array with MATLAB's class and size

    ``data`` is None for an array known only by its descriptor (``size``, ``complex`` and
    ``hash`` are then set).
    """

    kind = "array"
    __slots__ = ("mclass", "size", "complex", "data", "_hash")

    def __init__(
        self,
        mclass: str,
        data: np.ndarray | None = None,
        *,
        size: Sequence[int] | None = None,
        complex: bool = False,
        hash: str | None = None,
    ):
        if mclass not in _DTYPE:
            raise ValueError(f"unknown MATLAB array class {mclass!r}")
        if data is not None:
            data = np.asarray(data)
            if data.ndim < 2:
                raise ValueError("MArray data must have MATLAB size (2 or more dimensions)")
            if data.dtype.kind == "c":
                if mclass not in _CDTYPE or data.dtype.newbyteorder("=") != _CDTYPE[
                    mclass
                ].newbyteorder("="):
                    raise ValueError(
                        f"complex {mclass} needs dtype {_CDTYPE.get(mclass)}, got {data.dtype}"
                    )
                complex = True
            else:
                if data.dtype.newbyteorder("=") != _DTYPE[mclass].newbyteorder("="):
                    raise ValueError(f"{mclass} needs dtype {_DTYPE[mclass]}, got {data.dtype}")
                complex = False
            s = msize(data.shape)
            if s != data.shape:
                data = data.reshape(s, order="F")
            size = s
        elif size is None:
            raise ValueError("MArray needs data or a size")
        self.mclass = mclass
        self.size = msize(size)
        self.complex = bool(complex)
        self.data = data
        self._hash = hash

    @property
    def numel(self) -> int:
        return math.prod(self.size)

    @property
    def nbytes(self) -> int:
        return self.numel * _DTYPE[self.mclass].itemsize * (2 if self.complex else 1)

    def text(self) -> str:
        """a char row as str"""
        if self.mclass != "char" or self.data is None:
            raise ValueError("not a char array with data")
        return _units_to_str(self.data.reshape(-1, order="F"))

    def is_char_row(self) -> bool:
        return (
            self.mclass == "char" and len(self.size) == 2 and (self.size[0] == 1 or self.numel == 0)
        )

    def __repr__(self) -> str:
        c = " complex" if self.complex else ""
        if self.is_char_row() and self.data is not None:
            return f"MArray(char {size_text(self.size)} {self.text()!r})"
        return f"MArray({self.mclass}{c} {size_text(self.size)})"


class MCell(MValue):
    kind = "cell"
    __slots__ = ("size", "items", "_hash")

    def __init__(self, size: Sequence[int], items: list | None, hash: str | None = None):
        self.size = msize(size)
        if items is not None and len(items) != math.prod(self.size):
            raise ValueError(
                f"cell of size {size_text(self.size)} needs {math.prod(self.size)} items, got "
                f"{len(items)}"
            )
        self.items = items
        self._hash = hash

    def __repr__(self) -> str:
        return f"MCell({size_text(self.size)})"


class MStruct(MValue):
    """struct array: ``items`` holds one dict per element, in column-major order"""

    kind = "struct"
    __slots__ = ("size", "fields", "items", "_hash")

    def __init__(
        self,
        size: Sequence[int],
        fields: Sequence[str],
        items: list[dict] | None,
        hash: str | None = None,
    ):
        self.size = msize(size)
        self.fields = list(fields)
        if items is not None and len(items) != math.prod(self.size):
            raise ValueError("struct items do not match its size")
        self.items = items
        self._hash = hash

    def __repr__(self) -> str:
        return f"MStruct({size_text(self.size)} fields {self.fields})"


class MString(MValue):
    kind = "string"
    __slots__ = ("size", "items", "_hash")

    def __init__(self, size: Sequence[int], items: list[str | None], hash: str | None = None):
        self.size = msize(size)
        self.items = list(items) if items is not None else None
        self._hash = hash


class MSparse(MValue):
    """sparse matrix; rows and cols are 0-based, in MATLAB's (column-major) order"""

    kind = "sparse"
    __slots__ = ("mclass", "size", "rows", "cols", "values", "_hash")

    def __init__(
        self, mclass: str, size: Sequence[int], rows, cols, values, hash: str | None = None
    ):
        self.mclass = mclass
        self.size = msize(size)
        self.rows = None if rows is None else np.asarray(rows, dtype=np.int64)
        self.cols = None if cols is None else np.asarray(cols, dtype=np.int64)
        self.values = None if values is None else np.asarray(values)
        self._hash = hash


class MFunc(MValue):
    kind = "func"
    __slots__ = ("text", "_hash")

    def __init__(self, text: str, hash: str | None = None):
        self.text = text
        self._hash = hash


class MMap(MValue):
    kind = "map"
    __slots__ = ("key_type", "value_type", "keys", "values", "_hash")

    def __init__(
        self, key_type: str, value_type: str, keys: MCell, values: MCell, hash: str | None = None
    ):
        self.key_type, self.value_type, self.keys, self.values = key_type, value_type, keys, values
        self._hash = hash


class MObject(MValue):
    """an object; ``value`` is what MATLAB gave for it (an inputParser's Results), or None"""

    kind = "object"
    __slots__ = ("mclass", "value", "text", "_hash")

    def __init__(self, mclass: str, value: MValue | None, text: str = "", hash: str | None = None):
        self.mclass, self.value, self.text = mclass, value, text
        self._hash = hash


class MHashed(MValue):
    """a value known only by its kind, size and hash (a large cell, struct, string array or map
    in a manifest); ``info`` holds what else is known: struct fields, map key and value types,
    and ``sig``, the reading signature (see ``reading_signature``)"""

    kind = "hashed"
    __slots__ = ("of", "size", "info", "_hash")

    def __init__(self, of: str, size: Sequence[int] | None, hash: str, info: dict | None = None):
        self.of, self.size, self._hash = of, tuple(size) if size else None, hash
        self.info = info or {}

    def skeleton(self) -> MValue | None:
        """a value without data, of this kind, size and element kinds and classes, to read
        Python results as (as the full value would read them, when the signature is known)"""
        if self.info.get("sig") is not None:
            return _like_of_signature(self.info["sig"])
        if self.of == "cell" and self.size:
            return MCell(self.size, None)
        if self.of == "struct" and self.size:
            return MStruct(self.size, self.info.get("fields", []), None)
        if self.of == "string" and self.size:
            return MString(self.size, None)
        if self.of == "map":
            empty = MCell((1, 0), [])
            return MMap(self.info.get("kt", "char"), self.info.get("vt", "any"), empty, empty)
        return None

    @property
    def decides(self) -> bool:
        """whether a Python value is read as the full value would read it (so that a hash
        mismatch is a difference): cells, structs and maps need their reading signature"""
        return (
            self.of not in ("cell", "struct", "map", "object") or self.info.get("sig") is not None
        )


class MPlaceholder(MValue):
    """an argument filled in when a case runs: an input value (``input``), the path of an input
    file or folder in the work folder (``path``) or the work folder itself (``work``)"""

    kind = "placeholder"
    __slots__ = ("what", "name", "sub", "_hash")

    def __init__(self, what: str, name: str = "", sub: str = ""):
        if what not in ("input", "path", "work"):
            raise ValueError(what)
        self.what, self.name, self.sub, self._hash = what, name, sub, None

    def __repr__(self) -> str:
        name = {"input": "In", "path": "InPath", "work": "Work"}[self.what]
        return f"{name}({self.name or self.sub!r})"


class MRef(MValue):
    """an argument read from a value file (an input in the golden store)"""

    kind = "ref"
    __slots__ = ("path", "_hash")

    def __init__(self, path: str | os.PathLike, hash: str):
        self.path, self._hash = str(path), hash


# ---------------------------------------------------------------------------------------------
# char helpers
# ---------------------------------------------------------------------------------------------


def _str_to_units(s: str) -> np.ndarray:
    return np.frombuffer(s.encode("utf-16-le", "surrogatepass"), dtype="<u2").copy()


def _units_to_str(u: np.ndarray) -> str:
    return np.ascontiguousarray(u, dtype="<u2").tobytes().decode("utf-16-le", "surrogatepass")


def mchar(s: str) -> MArray:
    """MATLAB char row of a str ('' is 0x0, as MATLAB's '')"""
    u = _str_to_units(s)
    return MArray("char", u.reshape(1, -1) if u.size else np.zeros((0, 0), "<u2"))


# ---------------------------------------------------------------------------------------------
# raw data, hashing
# ---------------------------------------------------------------------------------------------


def _flat(a: np.ndarray) -> np.ndarray:
    f = a.reshape(-1, order="F")
    if f.dtype.byteorder == ">":
        f = f.astype(f.dtype.newbyteorder("<"))
    return f


def _real_parts(a: np.ndarray) -> list[np.ndarray]:
    f = _flat(a)
    return [f.real, f.imag] if f.dtype.kind == "c" else [f]


def _iter_raw(a: np.ndarray, chunk: int = 1 << 23) -> Iterator[memoryview]:
    """the raw little-endian bytes in column-major order (real block, then imaginary block)"""
    for p in _real_parts(a):
        for i in range(0, p.size, chunk):
            c = np.ascontiguousarray(p[i : i + chunk])
            if c.dtype == np.bool_:
                c = c.view(np.uint8)
            yield memoryview(c).cast("B")


def _raw_bytes(a: np.ndarray) -> bytes:
    return b"".join(bytes(m) for m in _iter_raw(a))


def _array_hash(a: MArray) -> str:
    h = _sha(b"a", a.mclass, _size_bytes(a.size), b"c" if a.complex else b"r", b"")
    for m in _iter_raw(a.data):
        h.update(m)
    return h.hexdigest()


def vhash(v: MValue) -> str:
    """SHA-256 of a value's class, size and contents (struct field order is not part of it)"""
    if getattr(v, "_hash", None):
        return v._hash
    if isinstance(v, MArray):
        if v.data is None:
            raise ValueError("array without data or hash")
        h = _array_hash(v)
    elif isinstance(v, MCell):
        h = _sha(b"cell", _size_bytes(v.size), ";".join(vhash(x) for x in v.items)).hexdigest()
    elif isinstance(v, MStruct):
        fs = sorted(v.fields)
        parts = ";".join(",".join(vhash(el[f]) for f in fs) for el in v.items)
        h = _sha(b"struct", _size_bytes(v.size), ",".join(fs), parts).hexdigest()
    elif isinstance(v, MString):
        h = _sha(b"string", _size_bytes(v.size), json.dumps(v.items)).hexdigest()
    elif isinstance(v, MSparse):
        vals = MArray(v.mclass, v.values.reshape(-1, 1))
        h = _sha(
            b"sparse",
            v.mclass,
            _size_bytes(v.size),
            v.rows.astype("<i8").tobytes(),
            v.cols.astype("<i8").tobytes(),
            vhash(vals),
        ).hexdigest()
    elif isinstance(v, MFunc):
        h = _sha(b"func", v.text).hexdigest()
    elif isinstance(v, MMap):
        h = _sha(b"map", v.key_type, v.value_type, vhash(v.keys), vhash(v.values)).hexdigest()
    elif isinstance(v, MObject):
        h = _sha(b"object", v.mclass, vhash(v.value) if v.value is not None else "-").hexdigest()
    elif isinstance(v, MPlaceholder):
        h = _sha(b"ph", v.what, v.name, v.sub).hexdigest()
    else:
        raise TypeError(f"cannot hash {type(v).__name__}")
    try:
        v._hash = h
    except AttributeError:
        pass
    return h


def sha256_file(path: str | os.PathLike, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------
# Python values to MATLAB values
# ---------------------------------------------------------------------------------------------


def _class_of(dtype: np.dtype) -> str:
    try:
        return _CLASS_OF_DTYPE[np.dtype(dtype).newbyteorder("=")]
    except KeyError:
        raise TypeError(f"no MATLAB class for dtype {dtype}") from None


def _array_from_numpy(x: np.ndarray, like: MValue | None = None) -> MValue:
    x = np.asarray(x)
    if x.dtype.kind in "US":
        if x.dtype.kind == "U" and x.dtype.itemsize == 4 and x.ndim >= 2:
            u = x.view(np.uint32)
            if np.all(u < 0x10000):
                return MArray("char", u.astype("<u2"))
        items = [str(s) for s in x.reshape(-1, order="F")]
        cells = [mchar(s) for s in items]
        return MCell(_size_like(x.shape, like), cells)
    if x.dtype.kind == "O":
        items = list(x.reshape(-1, order="F"))
        litems = (
            like.items
            if isinstance(like, MCell) and like.items and len(like.items) == len(items)
            else None
        )
        return MCell(
            _size_like(x.shape, like),
            [to_mvalue(it, litems[i] if litems else None) for i, it in enumerate(items)],
        )
    if x.dtype.kind == "f" and x.dtype.itemsize == 2:
        raise TypeError("float16 has no MATLAB class")
    cls = _class_of(x.dtype)
    if x.dtype.byteorder == ">":
        x = x.astype(x.dtype.newbyteorder("<"))
    return MArray(cls, x.reshape(_size_like(x.shape, like), order="F"))


def _size_like(shape: tuple, like: MValue | None) -> tuple[int, ...]:
    """MATLAB size of a NumPy shape; a 1-D shape becomes a row, or a column or empty size if the
    reference has one"""
    if len(shape) >= 2:
        return msize(shape)
    if len(shape) == 0:
        return (1, 1)
    n = int(shape[0])
    ls = getattr(like, "size", None)
    if ls and len(ls) == 2:
        if ls == (n, 1):
            return (n, 1)
        if n == 0 and math.prod(ls) == 0:
            return tuple(ls)
    return (1, n)


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float, np.integer, np.floating)) and not isinstance(
        x, (bool, np.bool_)
    )


def _int_fits(v: int, cls: str) -> bool:
    if cls == "double":
        return abs(v) <= 2**53
    if cls == "single":
        return abs(v) <= 2**24
    info = np.iinfo(_DTYPE[cls])
    return info.min <= v <= info.max


def _scalar(cls: str, v: Any) -> MArray:
    return MArray(
        cls, np.array([[v]], dtype=_CDTYPE[cls] if isinstance(v, complex) else _DTYPE[cls])
    )


def to_mvalue(x: Any, like: MValue | None = None) -> MValue:
    """the MATLAB value of a Python value

    Without ``like`` (arguments): None is [] (0x0 double); bool is logical; int and float are
    double; complex is complex double; str is a char row; numpy arrays and scalars keep their
    class (1-D arrays are rows); lists and tuples of numbers or of bools are double or logical
    rows, other lists are cell rows; dicts, dataclasses and namespaces are 1x1 structs.

    With ``like`` (a MATLAB reference), Python results are read as that reference where Python
    has no counterpart: a 1-D array or list matches a row or a column, a Python int takes an
    integer class it fits, '' and an empty 1-D array match any empty size, a list of dicts is a
    struct array, a dict is a containers.Map, and a value compared with a MATLAB object is
    compared with the object's value. Classes are never changed otherwise.
    """
    if isinstance(x, MValue):
        return x
    if isinstance(like, MObject) and like.value is not None:
        # compared with the object's value: read as that object, so that hashes agree too
        return MObject(like.mclass, to_mvalue(x, like.value), READ_AS_OBJECT)
    if isinstance(like, MHashed):
        like = like.skeleton()
    if x is None:
        if isinstance(like, MStruct) and math.prod(like.size) == 0:
            return MStruct(like.size, like.fields, [])
        if isinstance(like, MCell) and math.prod(like.size) == 0:
            return MCell(like.size, [])
        return MArray("double", np.zeros((0, 0)))
    if isinstance(x, (bool, np.bool_)):
        return _scalar("logical", bool(x))
    if isinstance(x, np.generic):
        return _array_from_numpy(np.asarray(x), like)
    if isinstance(x, int):
        cls = (
            like.mclass if isinstance(like, MArray) and like.mclass in NUMERIC_CLASSES else "double"
        )
        if not _int_fits(x, cls):
            if cls != "double" or not _int_fits(x, "double"):
                raise ValueError(f"Python int {x} does not fit MATLAB {cls}; pass a NumPy scalar")
            cls = "double"
        return _scalar(cls, x)
    if isinstance(x, float):
        return _scalar("double", x)
    if isinstance(x, complex):
        return _scalar("double", x)
    if isinstance(x, (str, os.PathLike)):
        s = os.fspath(x)
        if isinstance(like, MString):
            return MString((1, 1), [s])
        m = mchar(s)
        if not s and isinstance(like, MArray) and like.mclass == "char" and like.numel == 0:
            return MArray("char", np.zeros(like.size, "<u2"))
        return m
    if isinstance(x, bytes):
        return MArray("uint8", np.frombuffer(x, np.uint8).reshape(1, -1).copy())
    if isinstance(x, np.ndarray):
        return _array_from_numpy(x, like)
    if _is_sparse(x):
        return _sparse_from_scipy(x)
    if isinstance(x, (list, tuple)):
        return _list_to_mvalue(list(x), like)
    if isinstance(x, Mapping):
        if isinstance(like, MMap):
            return _dict_to_map(x, like)
        return _dict_to_struct(x, like)
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return _dict_to_struct({f.name: getattr(x, f.name) for f in dataclasses.fields(x)}, like)
    if hasattr(x, "__dict__") and type(x).__name__ == "SimpleNamespace":
        return _dict_to_struct(vars(x), like)
    raise TypeError(f"no MATLAB value for a Python {type(x).__name__}")


def _list_to_mvalue(xs: list, like: MValue | None) -> MValue:
    n = len(xs)
    if isinstance(like, MStruct) and all(isinstance(x, Mapping) for x in xs):
        size = _size_like((n,), like)
        litems = like.items if like.items is not None and len(like.items) == n else None
        els = [
            _dict_to_struct(x, MStruct((1, 1), like.fields, [litems[i]]) if litems else None)
            for i, x in enumerate(xs)
        ]
        fields = els[0].fields if els else like.fields
        return MStruct(size, fields, [e.items[0] for e in els])
    if isinstance(like, MString) and all(isinstance(x, str) or x is None for x in xs):
        return MString(_size_like((n,), like), xs)
    if isinstance(like, MCell):
        litems = like.items if like.items is not None and len(like.items) == n else None
        return MCell(
            _size_like((n,), like),
            [to_mvalue(x, litems[i] if litems else None) for i, x in enumerate(xs)],
        )
    if n == 0:
        if isinstance(like, MArray) and like.numel == 0:
            return MArray(like.mclass, np.zeros(like.size, _DTYPE[like.mclass]))
        return MArray("double", np.zeros((0, 0)))
    if all(isinstance(x, (bool, np.bool_)) for x in xs):
        return MArray(
            "logical", np.array(xs, dtype=bool).reshape(_size_like((n,), like), order="F")
        )
    if all(_is_number(x) for x in xs):
        if (
            all(isinstance(x, int) for x in xs)
            and isinstance(like, MArray)
            and like.mclass in NUMERIC_CLASSES
            and all(_int_fits(x, like.mclass) for x in xs)
        ):
            cls = like.mclass
        elif all(isinstance(x, (int, float)) for x in xs):
            cls = "double"
            for x in xs:
                if isinstance(x, int) and not _int_fits(x, "double"):
                    raise ValueError(f"Python int {x} does not fit MATLAB double")
        else:
            return MCell(_size_like((n,), like), [to_mvalue(x) for x in xs])
        return MArray(
            cls, np.array(xs, dtype=_DTYPE[cls]).reshape(_size_like((n,), like), order="F")
        )
    if all(isinstance(x, complex) or _is_number(x) for x in xs) and any(
        isinstance(x, complex) for x in xs
    ):
        return MArray(
            "double", np.array(xs, dtype=np.complex128).reshape(_size_like((n,), like), order="F")
        )
    return MCell(_size_like((n,), like), [to_mvalue(x) for x in xs])


def _dict_to_struct(d: Mapping, like: MValue | None) -> MStruct:
    lel = like.items[0] if isinstance(like, MStruct) and like.items and len(like.items) == 1 else {}
    items = {}
    for k, v in d.items():
        if not isinstance(k, str):
            raise TypeError(f"struct field names must be str, got {k!r}")
        items[k] = to_mvalue(v, lel.get(k))
    return MStruct((1, 1), list(d.keys()), [items])


def _dict_to_map(d: Mapping, like: MMap) -> MMap:
    keys = sorted(d.keys())
    kl = (
        like.keys.items
        if like.keys.items and len(like.keys.items) == len(keys)
        else [None] * len(keys)
    )
    vl = (
        like.values.items
        if like.values.items and len(like.values.items) == len(keys)
        else [None] * len(keys)
    )
    return MMap(
        like.key_type,
        like.value_type,
        MCell((1, len(keys)), [to_mvalue(k, kl[i]) for i, k in enumerate(keys)]),
        MCell((1, len(keys)), [to_mvalue(d[k], vl[i]) for i, k in enumerate(keys)]),
    )


def _is_sparse(x: Any) -> bool:
    return type(x).__module__.startswith("scipy.sparse")


def _sparse_from_scipy(x: Any) -> MSparse:
    c = x.tocsc()
    c.sort_indices()
    rows = c.indices.astype(np.int64)
    cols = np.repeat(np.arange(c.shape[1], dtype=np.int64), np.diff(c.indptr))
    vals = c.data
    cls = "logical" if vals.dtype == np.bool_ else _class_of(vals.dtype)
    return MSparse(cls, c.shape, rows, cols, vals)


# ---------------------------------------------------------------------------------------------
# value files (.pkv), exchanged with MATLAB (+pkpy_oracle/read_value.m, write_value.m)
# ---------------------------------------------------------------------------------------------


def write_pkv(path: str | os.PathLike, value: MValue | Any) -> None:
    """write a value file (written to a temporary name, then renamed)"""
    path = Path(path)
    tmp = path.with_name(f"{path.name}.tmp{os.getpid()}_{id(value)}")
    v = to_mvalue(value)
    try:
        with open(tmp, "wb") as f:
            f.write(MAGIC)
            root = _enc(v, f)
            js = json.dumps(
                {"pkv": 1, "root": root}, separators=(",", ":"), ensure_ascii=False
            ).encode("utf-8")
            f.write(js)
            f.write(struct.pack("<Q", len(js)))
            f.write(END_MAGIC)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def _align(f) -> int:
    o = f.tell()
    pad = (-o) % 8
    if pad:
        f.write(b"\x00" * pad)
    return o + pad


def _enc(v: MValue, f) -> dict:
    if isinstance(v, MArray):
        o = _align(f)
        n = 0
        for m in _iter_raw(v.data):
            f.write(m)
            n += m.nbytes
        return {"k": "a", "c": v.mclass, "s": list(v.size), "o": o, "n": n, "x": int(v.complex)}
    if isinstance(v, MCell):
        return {"k": "cell", "s": list(v.size), "v": [_enc(x, f) for x in v.items]}
    if isinstance(v, MStruct):
        return {
            "k": "struct",
            "s": list(v.size),
            "f": v.fields,
            "v": [_enc(el[fn], f) for el in v.items for fn in v.fields],
        }
    if isinstance(v, MString):
        return {"k": "string", "s": list(v.size), "v": v.items}
    if isinstance(v, MSparse):
        return {
            "k": "sparse",
            "c": v.mclass,
            "s": list(v.size),
            "i": _enc(MArray("double", (v.rows + 1).astype(np.float64).reshape(-1, 1)), f),
            "j": _enc(MArray("double", (v.cols + 1).astype(np.float64).reshape(-1, 1)), f),
            "v": _enc(MArray(v.mclass, v.values.reshape(-1, 1)), f),
        }
    if isinstance(v, MFunc):
        return {"k": "func", "t": v.text}
    if isinstance(v, MMap):
        return {
            "k": "map",
            "kt": v.key_type,
            "vt": v.value_type,
            "keys": _enc(v.keys, f),
            "values": _enc(v.values, f),
        }
    if isinstance(v, MObject):
        return {
            "k": "object",
            "c": v.mclass,
            "v": None if v.value is None else _enc(v.value, f),
            "t": v.text,
        }
    if isinstance(v, MPlaceholder):
        return {"k": "ph", "t": v.what, "n": v.name, "sub": v.sub}
    if isinstance(v, MRef):
        return {"k": "ref", "p": v.path, "h": v._hash}
    raise TypeError(f"cannot write {type(v).__name__}")


def read_pkv(path: str | os.PathLike) -> MValue:
    """read a value file written by MATLAB or by write_pkv"""
    with open(path, "rb") as f:
        if f.read(8) != MAGIC:
            raise ValueError(f"{path}: not a value file")
        f.seek(-16, os.SEEK_END)
        n = struct.unpack("<Q", f.read(8))[0]
        if f.read(8) != END_MAGIC:
            raise ValueError(f"{path}: value file is incomplete")
        f.seek(-16 - n, os.SEEK_END)
        hdr = json.loads(f.read(n).decode("utf-8"))
        return _dec(hdr["root"], f)


def _seq(v: Any) -> list:
    # MATLAB's jsonencode writes a one-element cell as [x]; a missing list is null
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


def _dec(node: dict, f) -> MValue:
    k = node["k"]
    if k == "a":
        cls, size = node["c"], msize(_seq(node["s"]))
        cx = bool(node.get("x"))
        n = int(node["n"])
        f.seek(int(node["o"]))
        dt = _DTYPE[cls]
        cnt = math.prod(size)
        raw = np.fromfile(f, dtype=dt if dt != np.bool_ else np.uint8, count=cnt * (2 if cx else 1))
        if raw.nbytes != n:
            raise ValueError(f"value file: array of {n} bytes has {raw.nbytes}")
        if dt == np.bool_:
            raw = raw.astype(bool)
        if cx:
            data = np.empty(cnt, dtype=_CDTYPE[cls])
            data.real, data.imag = raw[:cnt], raw[cnt:]
        else:
            data = raw
        return MArray(cls, data.reshape(size, order="F"))
    if k == "cell":
        size = msize(_seq(node["s"]))
        return MCell(size, [_dec(x, f) for x in _seq(node.get("v"))])
    if k == "struct":
        size = msize(_seq(node["s"]))
        fields = [str(x) for x in _seq(node.get("f"))]
        vals = [_dec(x, f) for x in _seq(node.get("v"))]
        nf = len(fields)
        items = (
            [
                dict(zip(fields, vals[i * nf : (i + 1) * nf], strict=True))
                for i in range(math.prod(size))
            ]
            if nf
            else [{} for _ in range(math.prod(size))]
        )
        return MStruct(size, fields, items)
    if k == "string":
        size = msize(_seq(node["s"]))
        items = node.get("v")
        items = [] if items is None else (items if isinstance(items, list) else [items])
        items = [None if (x is None or x == []) else str(x) for x in items]
        return MString(size, items)
    if k == "sparse":
        i = _dec(node["i"], f).data.reshape(-1).astype(np.int64) - 1
        j = _dec(node["j"], f).data.reshape(-1).astype(np.int64) - 1
        vals = _dec(node["v"], f).data.reshape(-1)
        return MSparse(node["c"], _seq(node["s"]), i, j, vals)
    if k == "func":
        return MFunc(node["t"])
    if k == "map":
        return MMap(node["kt"], node["vt"], _dec(node["keys"], f), _dec(node["values"], f))
    if k == "object":
        v = node.get("v")
        return MObject(node["c"], _dec(v, f) if isinstance(v, dict) else None, node.get("t") or "")
    if k == "ph":
        return MPlaceholder(node["t"], node.get("n") or "", node.get("sub") or "")
    if k == "ref":
        return read_pkv(node["p"])
    raise ValueError(f"value file: unknown node kind {k!r}")


# ---------------------------------------------------------------------------------------------
# descriptors (repo manifests)
# ---------------------------------------------------------------------------------------------

HASH_ONLY = (
    f"the manifest keeps the data of arrays up to {INLINE_BYTES} bytes and the elements of "
    f"containers up to {DESCRIPTOR_BYTES // 1024} KB, else hashes: the readable element-by-element "
    "diff is in the golden store (run this test on the cluster)"
)


_MIN_NODE = 40  # bytes of the smallest element descriptor
READ_AS_OBJECT = "(a Python value, read as this object's value)"


def _runs(sigs: list) -> list:
    out: list = []
    for s in sigs:
        if out and out[-1][1] == s:
            out[-1][0] += 1
        else:
            out.append([1, s])
    return out


def _unruns(runs: list) -> list:
    return [s for n, s in runs for _ in range(n)]


def reading_signature(v: MValue) -> Any:
    """what reading a Python value as ``v`` uses (``to_mvalue``'s rules), for hash-only
    descriptors: kinds, array classes, the sizes of empty arrays and of columns, struct fields,
    map types; element lists run-length encoded"""
    if isinstance(v, MArray):
        s = v.size
        if v.numel == 0 or (len(s) == 2 and s[1] == 1 and s[0] > 1):
            return f"{v.mclass}:{'x'.join(str(n) for n in s)}"
        return v.mclass
    if isinstance(v, MCell):
        return {"k": "cell", "s": list(v.size), "v": _runs([reading_signature(x) for x in v.items])}
    if isinstance(v, MStruct):
        cols = {f: _runs([reading_signature(el[f]) for el in v.items]) for f in v.fields}
        return {"k": "struct", "s": list(v.size), "f": v.fields, "v": cols}
    if isinstance(v, MString):
        return {"k": "string", "s": list(v.size)}
    if isinstance(v, MMap):
        return {
            "k": "map",
            "kt": v.key_type,
            "vt": v.value_type,
            "keys": reading_signature(v.keys),
            "values": reading_signature(v.values),
        }
    if isinstance(v, MObject):
        return {
            "k": "object",
            "c": v.mclass,
            "v": None if v.value is None else reading_signature(v.value),
        }
    if isinstance(v, MHashed) and v.info.get("sig") is not None:
        return v.info["sig"]
    return {"k": v.kind}


def _like_of_signature(sig: Any) -> MValue | None:
    """a value without data that reads Python values as the value of this signature does"""
    if isinstance(sig, str):
        cls, _, size = sig.partition(":")
        return MArray(cls, None, size=tuple(int(n) for n in size.split("x")) if size else (1, 1))
    k = sig.get("k")
    if k == "cell":
        return MCell(sig["s"], [_like_of_signature(x) for x in _unruns(sig["v"])])
    if k == "struct":
        cols = {f: _unruns(sig["v"][f]) for f in sig["f"]}
        n = math.prod(sig["s"])
        return MStruct(
            sig["s"],
            sig["f"],
            [{f: _like_of_signature(cols[f][i]) for f in sig["f"]} for i in range(n)],
        )
    if k == "string":
        return MString(sig["s"], None)
    if k == "map":
        keys, values = _like_of_signature(sig["keys"]), _like_of_signature(sig["values"])
        empty = MCell((1, 0), [])
        return MMap(
            sig["kt"],
            sig["vt"],
            keys if isinstance(keys, MCell) else empty,
            values if isinstance(values, MCell) else empty,
        )
    if k == "object":
        return MObject(sig["c"], None if sig.get("v") is None else _like_of_signature(sig["v"]))
    return None


def _with_signature(d: dict, v: MValue, budget: int) -> dict:
    """a hash-only container descriptor with its reading signature, when that fits"""
    if "v" not in d and "keys" not in d:
        sig = reading_signature(v)
        if _fits(sig, budget):
            d["r"] = sig
    return d


def _fits(obj: Any, budget: int) -> bool:
    return len(json.dumps(obj, separators=(",", ":"))) <= budget


def _derived(d: dict) -> dict:
    """an element's descriptor without the hash that its data gives again"""
    k = d["k"]
    if (
        (k == "a" and ("b64" in d or "t" in d))
        or (k in ("cell", "struct", "string") and "v" in d)
        or (k == "map" and "keys" in d)
    ):
        d.pop("h", None)
    return d


def describe(v: MValue, inline: int = INLINE_BYTES, budget: int = DESCRIPTOR_BYTES) -> dict:
    """JSON descriptor: kind, class, size, hash, and the data of small arrays (up to ``inline``
    bytes; char rows as text); a cell, struct, string array or map whose descriptor would take
    more than ``budget`` bytes keeps only its kind, size (struct fields, map types) and hash.
    Elements whose data is kept carry no hash of their own (it is computed again)."""
    h = vhash(v)
    if isinstance(v, MArray):
        d = {"k": "a", "c": v.mclass, "s": list(v.size), "h": h}
        if v.complex:
            d["x"] = 1
        if v.data is not None and v.nbytes <= inline:
            if v.mclass == "char" and v.is_char_row():
                d["t"] = v.text()
            else:
                d["b64"] = base64.b64encode(_raw_bytes(v.data)).decode("ascii")
        return d
    if isinstance(v, MCell):
        d = {"k": "cell", "s": list(v.size), "h": h}
        if v.items is not None and len(v.items) * _MIN_NODE <= budget:
            kids = [_derived(describe(x, inline, budget)) for x in v.items]
            if _fits(kids, budget):
                d["v"] = kids
        return _with_signature(d, v, budget)
    if isinstance(v, MStruct):
        d = {"k": "struct", "s": list(v.size), "f": v.fields, "h": h}
        if v.items is not None and len(v.items) * max(len(v.fields), 1) * _MIN_NODE <= budget:
            kids = [
                {fn: _derived(describe(el[fn], inline, budget)) for fn in v.fields}
                for el in v.items
            ]
            if _fits(kids, budget):
                d["v"] = kids
        return _with_signature(d, v, budget)
    if isinstance(v, MString):
        d = {"k": "string", "s": list(v.size), "h": h}
        if v.items is not None and _fits(v.items, budget):
            d["v"] = v.items
        return d
    if isinstance(v, MObject):
        return {
            "k": "object",
            "c": v.mclass,
            "h": h,
            "v": None if v.value is None else describe(v.value, inline, budget),
        }
    if isinstance(v, MFunc):
        return {"k": "func", "t": v.text, "h": h}
    if isinstance(v, MSparse):
        return {"k": "sparse", "c": v.mclass, "s": list(v.size), "h": h, "nnz": int(v.rows.size)}
    if isinstance(v, MMap):
        d = {"k": "map", "h": h, "kt": v.key_type, "vt": v.value_type}
        kd, vd = describe(v.keys, inline, budget), describe(v.values, inline, budget)
        if "v" in kd and "v" in vd and _fits([kd, vd], budget):
            d["keys"], d["values"] = _derived(kd), _derived(vd)
        return _with_signature(d, v, budget)
    if isinstance(v, MHashed):
        d = {"k": v.of, "s": list(v.size) if v.size else None, "h": h}
        if "fields" in v.info:
            d["f"] = v.info["fields"]
        if "kt" in v.info:
            d["kt"], d["vt"] = v.info["kt"], v.info.get("vt")
        if v.info.get("sig") is not None:
            d["r"] = v.info["sig"]
        return d
    raise TypeError(f"cannot describe {type(v).__name__}")


def from_descriptor(d: dict) -> MValue:
    """the value a descriptor stands for (arrays without inline data keep only their hash);
    reads descriptors of every manifest layout (elements with or without their own hash)"""
    k = d["k"]
    if k == "a":
        cls, size, cx = d["c"], msize(d["s"]), bool(d.get("x"))
        if "t" in d:
            return MArray("char", _str_to_units(d["t"]).reshape(size, order="F"), hash=d.get("h"))
        if "b64" in d:
            raw = base64.b64decode(d["b64"])
            dt = _DTYPE[cls]
            cnt = math.prod(size)
            a = np.frombuffer(raw, dtype=dt if dt != np.bool_ else np.uint8).copy()
            if dt == np.bool_:
                a = a.astype(bool)
            if cx:
                c = np.empty(cnt, dtype=_CDTYPE[cls])
                c.real, c.imag = a[:cnt], a[cnt:]
                a = c
            return MArray(cls, a.reshape(size, order="F"), hash=d.get("h"))
        return MArray(cls, None, size=size, complex=cx, hash=d["h"])
    if k == "cell":
        if "v" not in d:
            return MHashed("cell", d["s"], d["h"], {"sig": d.get("r")})
        return MCell(d["s"], [from_descriptor(x) for x in d["v"]], hash=d.get("h"))
    if k == "struct":
        if "v" not in d:
            return MHashed("struct", d["s"], d["h"], {"fields": list(d["f"]), "sig": d.get("r")})
        return MStruct(
            d["s"],
            d["f"],
            [{fn: from_descriptor(el[fn]) for fn in d["f"]} for el in d["v"]],
            hash=d.get("h"),
        )
    if k == "string":
        if "v" not in d:
            return MHashed("string", d["s"], d["h"])
        return MString(d["s"], d["v"], hash=d.get("h"))
    if k == "object":
        return MObject(d["c"], None if d.get("v") is None else from_descriptor(d["v"]), hash=d["h"])
    if k == "func":
        return MFunc(d["t"], hash=d["h"])
    if k == "map" and "keys" in d:
        keys, values = from_descriptor(d["keys"]), from_descriptor(d["values"])
        if isinstance(keys, MCell) and isinstance(values, MCell):
            return MMap(d["kt"], d["vt"], keys, values, hash=d.get("h"))
    if k == "map":
        info = {"kt": d.get("kt", "char"), "vt": d.get("vt", "any"), "sig": d.get("r")}
        return MHashed("map", None, d["h"], info)
    return MHashed(k, d.get("s"), d["h"])


# ---------------------------------------------------------------------------------------------
# masks (volatile text)
# ---------------------------------------------------------------------------------------------

UUID_RE = r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
STAMP_RE = r"(?<!\d)(?:19|20)\d{12}(?!\d)"
TIME_RE = (
    r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:?\d{2})?"
    r"|\d{2}-[A-Z][a-z]{2}-\d{4} \d{2}:\d{2}:\d{2}"
)


@dataclasses.dataclass(frozen=True)
class Masks:
    """text substitutions applied to both sides before comparing: (regex, replacement) pairs"""

    subs: tuple[tuple[str, str], ...] = (
        (UUID_RE, "<uuid>"),
        (STAMP_RE, "<timestamp>"),
        (TIME_RE, "<time>"),
    )

    def text(self, s: str) -> str:
        for pat, rep in self.subs:
            s = re.sub(pat, rep, s)
        return s

    def plus(self, *subs: tuple[str, str]) -> Masks:
        return Masks(tuple(subs) + self.subs)

    def with_paths(
        self, work: str | None = None, host: str | None = None, **roots: str | None
    ) -> Masks:
        """mask a side's own work folder (``<WORK>``), host name and other roots"""
        subs = []
        for tag, p in [("WORK", work)] + [(k.upper(), v) for k, v in roots.items()]:
            if not p:
                continue
            p = str(p).rstrip("/\\")
            for form in {p, p.replace("\\", "/"), p.replace("/", "\\")}:
                subs.append((re.escape(form), f"<{tag}>"))
        if host:
            subs.append((r"(?<![\w.-])" + re.escape(host) + r"(?![\w-])", "<host>"))
        return Masks(tuple(subs) + self.subs)

    def key(self) -> str:
        return hashlib.sha256(json.dumps(self.subs).encode()).hexdigest()[:16]


NO_MASKS = Masks(())


def mask_value(v: MValue, masks: Masks) -> MValue:
    """the value with masks applied to its char rows and strings"""
    if not masks.subs:
        return v
    if isinstance(v, MArray):
        if v.mclass == "char" and v.data is not None and v.is_char_row():
            t = v.text()
            m = masks.text(t)
            return v if m == t else mchar(m)
        return v
    if isinstance(v, MCell) and v.items is not None:
        items = [mask_value(x, masks) for x in v.items]
        return (
            v if all(a is b for a, b in zip(items, v.items, strict=True)) else MCell(v.size, items)
        )
    if isinstance(v, MStruct) and v.items is not None:
        items = [{k: mask_value(x, masks) for k, x in el.items()} for el in v.items]
        same = all(items[i][k] is v.items[i][k] for i in range(len(items)) for k in v.fields)
        return v if same else MStruct(v.size, v.fields, items)
    if isinstance(v, MString) and v.items is not None:
        items = [None if s is None else masks.text(s) for s in v.items]
        return v if items == v.items else MString(v.size, items)
    if isinstance(v, MObject) and v.value is not None:
        nv = mask_value(v.value, masks)
        return v if nv is v.value else MObject(v.mclass, nv, v.text)
    if isinstance(v, MMap):
        return MMap(
            v.key_type, v.value_type, mask_value(v.keys, masks), mask_value(v.values, masks)
        )
    return v


def drop_fields(v: MValue, names: Iterable[str]) -> MValue:
    """the value without struct fields of these names, at any depth"""
    names = set(names)
    if not names:
        return v
    if isinstance(v, MStruct) and v.items is not None:
        keep = [f for f in v.fields if f not in names]
        return MStruct(
            v.size, keep, [{f: drop_fields(el[f], names) for f in keep} for el in v.items]
        )
    if isinstance(v, MCell) and v.items is not None:
        return MCell(v.size, [drop_fields(x, names) for x in v.items])
    if isinstance(v, MObject) and v.value is not None:
        return MObject(v.mclass, drop_fields(v.value, names), v.text)
    return v


# ---------------------------------------------------------------------------------------------
# comparison specifications
# ---------------------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Exact:
    """bit for bit: class, size and every bit of every value (NaN payloads and -0 included)"""

    def describe(self) -> str:
        return "exact"


@dataclasses.dataclass(frozen=True)
class Tolerance:
    """a declared tolerance; ``reason`` names the parity exception that allows it

    A differing value is within the tolerance if its absolute difference is at most the
    largest of ``abs``, ``rel`` times the reference value and ``rel_max`` times the largest
    reference magnitude, or if it is at most ``ulp`` units in the last place away. The values
    that differ at all may be at most ``max_count`` and at most ``max_frac`` of all values, and
    the relative L2 error at most ``rel_l2``. Class and size must still be equal; NaN must be
    NaN where MATLAB has NaN.
    """

    reason: str
    abs: float | None = None
    rel: float | None = None
    rel_max: float | None = None
    ulp: int | None = None
    rel_l2: float | None = None
    max_frac: float | None = None
    max_count: int | None = None

    def __post_init__(self):
        if not self.reason or not self.reason.strip():
            raise ValueError("a tolerance needs a reason: the parity exception that allows it")

    def describe(self) -> str:
        parts = [
            f"{f.name}={getattr(self, f.name)}"
            for f in dataclasses.fields(self)
            if f.name != "reason" and getattr(self, f.name) is not None
        ]
        return f"tolerance ({', '.join(parts)}; {self.reason})"

    @property
    def elementwise(self) -> bool:
        return any(x is not None for x in (self.abs, self.rel, self.rel_max, self.ulp))


@dataclasses.dataclass(frozen=True)
class Ignore:
    reason: str

    def describe(self) -> str:
        return f"not compared ({self.reason})"


@dataclasses.dataclass(frozen=True)
class Custom:
    """``fn(ref, act)`` returns a list of problems (strings); empty when equal"""

    fn: Callable[[MValue, MValue], list[str]]
    name: str = "custom"

    def describe(self) -> str:
        return self.name


Spec = Exact | Tolerance | Ignore | Custom
EXACT = Exact()


def exact() -> Exact:
    return EXACT


def tol(
    *,
    reason: str,
    abs: float | None = None,
    rel: float | None = None,
    rel_max: float | None = None,
    ulp: int | None = None,
    rel_l2: float | None = None,
    max_frac: float | None = None,
    max_count: int | None = None,
) -> Tolerance:
    return Tolerance(reason, abs, rel, rel_max, ulp, rel_l2, max_frac, max_count)


def ignore(reason: str) -> Ignore:
    return Ignore(reason)


# ---------------------------------------------------------------------------------------------
# reports
# ---------------------------------------------------------------------------------------------


@dataclasses.dataclass
class Problem:
    path: str
    what: str  # class, size, kind, values, fields, error, missing, extra, text, layout, meta, ...
    lines: list[str]


@dataclasses.dataclass
class Report:
    problems: list[Problem] = dataclasses.field(default_factory=list)
    notes: list[str] = dataclasses.field(default_factory=list)
    undecided: list[str] = dataclasses.field(default_factory=list)  # needs the golden store

    @property
    def ok(self) -> bool:
        return not self.problems

    def add(self, path: str, what: str, *lines: str) -> None:
        self.problems.append(Problem(path, what, [ln for ln in lines if ln]))

    def note(self, path: str, text: str) -> None:
        self.notes.append(f"{path}: {text}" if path else text)

    def undecide(self, path: str, why: str) -> None:
        self.undecided.append(f"{path}: {why}")

    def extend(self, other: Report) -> None:
        self.problems.extend(other.problems)
        self.notes.extend(other.notes)
        self.undecided.extend(other.undecided)

    def text(self, title: str = "") -> str:
        out = []
        if title:
            out.append(f"{title}: " + ("equal" if self.ok else f"{len(self.problems)} problem(s)"))
        for p in self.problems[:MAX_PROBLEMS]:
            head, *rest = p.lines or [""]
            out.append(f"  {p.path}: {p.what}: {head}")
            out.extend(f"      {ln}" for ln in rest)
        if len(self.problems) > MAX_PROBLEMS:
            out.append(f"  ... and {len(self.problems) - MAX_PROBLEMS} more problems")
        if self.undecided:
            out.append("  not decided without the golden store (run on the cluster):")
            out.extend(f"    {u}" for u in self.undecided[:MAX_PROBLEMS])
        if self.notes:
            out.append("  notes:")
            out.extend(f"    {n}" for n in self.notes[:MAX_PROBLEMS])
            if len(self.notes) > MAX_PROBLEMS:
                out.append(f"    ... and {len(self.notes) - MAX_PROBLEMS} more notes")
        return "\n".join(out)

    def to_json(self) -> dict:
        return {
            "ok": self.ok,
            "problems": [dataclasses.asdict(p) for p in self.problems],
            "notes": self.notes,
            "undecided": self.undecided,
        }


# ---------------------------------------------------------------------------------------------
# comparison of values
# ---------------------------------------------------------------------------------------------


def _kind_name(v: MValue) -> str:
    if isinstance(v, MArray):
        return f"{v.mclass}{' complex' if v.complex else ''} array {size_text(v.size)}"
    size = getattr(v, "size", None)
    kind = v.of if isinstance(v, MHashed) else v.kind
    name = {
        "cell": "cell array",
        "struct": "struct",
        "string": "string array",
        "sparse": "sparse matrix",
        "func": "function handle",
        "map": "containers.Map",
        "object": f"{getattr(v, 'mclass', '')} object",
    }
    return f"{name.get(kind, kind)}{' ' + size_text(size) if size else ''}"


def compare(
    ref: MValue,
    act: Any,
    spec: Spec = EXACT,
    path: str = "value",
    *,
    masks: Masks | None = None,
    act_masks: Masks | None = None,
) -> Report:
    """compare a Python value (converted with ``to_mvalue(act, like=ref)``) with a MATLAB value"""
    rep = Report()
    a = to_mvalue(act, like=ref)
    if masks is not None:
        ref = mask_value(ref, masks)
    if act_masks is not None or masks is not None:
        a = mask_value(a, act_masks if act_masks is not None else masks)
    _cmp(ref, a, spec, path, rep)
    return rep


def _cmp(r: MValue, a: MValue, spec: Spec, path: str, rep: Report) -> None:
    if isinstance(spec, Ignore):
        return
    if isinstance(spec, Custom):
        if isinstance(r, MHashed) or (isinstance(r, MArray) and r.data is None):
            if vhash(a) != r._hash:
                rep.undecide(
                    path, f"{spec.name} needs MATLAB's data; the manifest has its hash only"
                )
            return
        for p in spec.fn(r, a):
            rep.add(path, spec.name, p)
        return
    if isinstance(r, MObject) and not isinstance(a, MObject):
        if r.value is None:
            rep.add(
                path,
                "kind",
                f"MATLAB returned a {r.mclass} object that the harness cannot read, Python a "
                f"{_kind_name(a)}",
            )
            return
        rep.note(path, f"MATLAB returned a {r.mclass} object; compared its value")
        r = r.value
    if isinstance(r, MHashed):
        asize = getattr(a, "size", None)
        if getattr(a, "kind", None) != r.of:
            rep.add(path, "kind", f"MATLAB: {_kind_name(r)}, Python: {_kind_name(a)}")
        elif r.size is not None and asize is not None and tuple(asize) != tuple(r.size):
            rep.add(
                path,
                "size",
                f"{r.of} size: MATLAB {size_text(r.size)}, Python {size_text(asize)}",
                _size_hint(r.size, asize),
            )
        elif r.of == "struct" and set(r.info.get("fields", [])) != set(a.fields):
            fr, fa = set(r.info["fields"]), set(a.fields)
            rep.add(
                path,
                "fields",
                f"struct fields differ: only in MATLAB {sorted(fr - fa)}, only in Python "
                f"{sorted(fa - fr)}",
            )
        elif vhash(a) != r._hash:
            if not isinstance(spec, Exact):
                rep.undecide(
                    path, f"{spec.describe()} needs MATLAB's data; the manifest has its hash only"
                )
            elif not r.decides:
                rep.undecide(path, "the manifest lacks the reading signature of this large value")
            else:
                rep.add(
                    path,
                    "values",
                    f"the {_kind_name(r)} differs from MATLAB's (SHA-256 of all its elements)",
                    HASH_ONLY,
                )
        return
    if r.kind != a.kind:
        hint = ""
        if isinstance(r, MArray) and r.mclass == "char" and isinstance(a, MString):
            hint = "MATLAB has a char array, Python a string: return a str"
        elif isinstance(r, MString) and isinstance(a, MArray) and a.mclass == "char":
            hint = "MATLAB has a string array, Python a char array"
        elif isinstance(r, MSparse):
            hint = "MATLAB has a sparse matrix: return a scipy.sparse matrix"
        rep.add(path, "kind", f"MATLAB: {_kind_name(r)}, Python: {_kind_name(a)}", hint)
        return
    if r._hash is not None and isinstance(spec, Exact) and not isinstance(r, MArray):
        if vhash(a) == r._hash:
            return
    if isinstance(r, MArray):
        _cmp_array(r, a, spec, path, rep)
    elif isinstance(r, MCell):
        if r.size != a.size:
            rep.add(
                path,
                "size",
                f"cell size: MATLAB {size_text(r.size)}, Python {size_text(a.size)}",
                _size_hint(r.size, a.size),
            )
            return
        for i, (x, y) in enumerate(zip(r.items, a.items, strict=True)):
            _cmp(x, y, spec, f"{path}{{{_sub_text(i, r.size)}}}", rep)
    elif isinstance(r, MStruct):
        fr, fa = set(r.fields), set(a.fields)
        if fr != fa:
            rep.add(
                path,
                "fields",
                f"struct fields differ: only in MATLAB {sorted(fr - fa)}, only in Python "
                f"{sorted(fa - fr)}",
            )
        elif r.fields != a.fields and len(r.size) == 2 and r.size == a.size:
            rep.note(
                path, f"field order differs (MATLAB {r.fields}, Python {a.fields}); not compared"
            )
        if r.size != a.size:
            rep.add(
                path,
                "size",
                f"struct array size: MATLAB {size_text(r.size)}, Python {size_text(a.size)}",
            )
            return
        common = [f for f in r.fields if f in fa]
        n = len(r.items)
        for i in range(n):
            el = "" if n == 1 else f"({_sub_text(i, r.size)})"
            for f in common:
                _cmp(r.items[i][f], a.items[i][f], spec, f"{path}{el}.{f}", rep)
    elif isinstance(r, MString):
        if r.size != a.size:
            rep.add(
                path,
                "size",
                f"string array size: MATLAB {size_text(r.size)}, Python {size_text(a.size)}",
            )
            return
        bad = [i for i, (x, y) in enumerate(zip(r.items, a.items, strict=True)) if x != y]
        if bad:
            i = bad[0]
            rep.add(
                path,
                "values",
                f"{len(bad)} of {len(r.items)} strings differ; first at ({_sub_text(i, r.size)}): "
                f"MATLAB {r.items[i]!r}, Python {a.items[i]!r}",
            )
    elif isinstance(r, MSparse):
        if (r.mclass, r.size) != (a.mclass, a.size):
            rep.add(
                path,
                "class",
                f"sparse: MATLAB {r.mclass} {size_text(r.size)}, Python {a.mclass} "
                f"{size_text(a.size)}",
            )
        elif not (np.array_equal(r.rows, a.rows) and np.array_equal(r.cols, a.cols)):
            rep.add(
                path,
                "values",
                f"sparse pattern differs: MATLAB {r.rows.size} nonzeros, Python {a.rows.size}",
            )
        else:
            _cmp_array(
                MArray(r.mclass, r.values.reshape(-1, 1)),
                MArray(a.mclass, a.values.reshape(-1, 1)),
                spec,
                f"{path} (nonzeros)",
                rep,
            )
    elif isinstance(r, MFunc):
        if r.text != a.text:
            rep.add(path, "values", f"function handle: MATLAB {r.text!r}, Python {a.text!r}")
    elif isinstance(r, MMap):
        _cmp(r.keys, a.keys, EXACT, f"{path}.keys", rep)
        _cmp(r.values, a.values, spec, f"{path}.values", rep)
    elif isinstance(r, MObject):
        if a.text == READ_AS_OBJECT and r.mclass == a.mclass:
            rep.note(path, f"MATLAB returned a {r.mclass} object; compared its value")
        if r.mclass != a.mclass:
            rep.add(path, "class", f"object class: MATLAB {r.mclass}, Python {a.mclass}")
        elif r.value is not None and a.value is not None:
            _cmp(r.value, a.value, spec, path, rep)


def _sub_text(i: int, size: Sequence[int]) -> str:
    if len(size) == 2 and (size[0] == 1 or size[1] == 1):
        return str(i + 1)
    return ",".join(str(s) for s in subscripts(i, size))


MISSING_MAP = (
    "Python's result is MATLAB's in (z, y, x) form: a volume not mapped back to MATLAB's "
    "(out role 'vol' missing, or to_matlab left out of py_out)?"
)
DOUBLE_MAP = (
    "Python's result is MATLAB's mapped once more (to_matlab of MATLAB's (y, x, z) array): a "
    "volume mapped twice (out role 'vol' and to_matlab in py_out), or a (y, x, z) argument "
    "without its Vol role?"
)
YX_SWAP = (
    "Python's result is MATLAB's with y and x swapped (dimensions 1 and 2): x and y exchanged "
    "in an index, a per-axis vector or a transpose?"
)


def _size_hint(rs: Sequence[int], as_: Sequence[int]) -> str:
    rs, as_ = tuple(rs), tuple(as_)
    if sorted(rs) == sorted(as_) and len(rs) == len(as_) and rs != as_:
        n = len(rs)
        named = ""
        if n >= 3 and as_ == tuple(rs[i] for i in _perm_from_matlab(n)):
            named = MISSING_MAP
        elif n >= 3 and as_ == tuple(rs[i] for i in _perm_to_matlab(n)):
            named = DOUBLE_MAP
        elif n >= 2 and min(rs[:2]) > 1 and as_ == (rs[1], rs[0], *rs[2:]):
            named = YX_SWAP
        return (
            "the sizes are a permutation of each other: axes in another order ((y, x, z) is "
            "MATLAB's)?" + (f" {named}" if named else "")
        )
    if len(rs) == len(as_):
        d = [b - a for a, b in zip(rs, as_, strict=True)]
        nz = [i for i, x in enumerate(d) if x]
        if len(nz) == 1:
            return (
                f"dimension {nz[0] + 1} ({AXES[nz[0]] if nz[0] < 4 else nz[0] + 1}) differs by "
                f"{d[nz[0]]:+d}: an off-by-one in a size or range?"
            )
    if len(rs) == 2 and len(as_) == 2 and rs == tuple(reversed(as_)) and 1 in rs:
        return "row against column vector"
    return ""


def _fmt(x: Any, cls: str) -> str:
    if cls == "char":
        try:
            return repr(chr(int(x)))
        except (ValueError, OverflowError):
            return str(int(x))
    if isinstance(x, (np.floating, float)):
        if np.isnan(x):
            return f"NaN(0x{_bits_hex(x)})"
        if cls == "single":
            return (
                repr(float(np.float32(x)))
                if np.float32(x) == float(np.float32(x))
                else str(np.float32(x))
            )
        return repr(float(x))
    if isinstance(x, (np.complexfloating, complex)):
        return (
            f"{_fmt(x.real, cls)}{'+' if not np.signbit(x.imag) else '-'}{_fmt(abs(x.imag), cls)}i"
        )
    if isinstance(x, (np.bool_, bool)):
        return "true" if x else "false"
    return str(x)


def _bits_hex(x: Any) -> str:
    a = np.asarray(x)
    return f"{int(a.view(_UINT[a.dtype.itemsize])):0{a.dtype.itemsize * 2}X}"


def _neq_mask(rf: np.ndarray, af: np.ndarray) -> np.ndarray:
    """differing bits, per element (inputs: flat arrays of one dtype)"""
    m = np.zeros(rf.shape, dtype=bool)
    for pr, pa in zip(_real_parts(rf), _real_parts(af), strict=True):
        pr = np.ascontiguousarray(pr)
        pa = np.ascontiguousarray(pa)
        if pr.dtype == np.bool_:
            m |= pr != pa
        else:
            u = _UINT[pr.dtype.itemsize]
            m |= pr.view(u) != pa.view(u)
    return m


def _bits_equal(x: np.ndarray, y: np.ndarray) -> bool:
    if x.shape != y.shape or x.dtype != y.dtype:
        return False
    return not _neq_mask(_flat(x), _flat(y)).any()


def _ordered(x: np.ndarray) -> np.ndarray:
    """int64 keys in the order of the float values (for ULP distances)"""
    if x.dtype.itemsize == 4:
        i = x.view(np.int32).astype(np.int64)
        return np.where(i < 0, np.int64(-(2**31)) - i, i)
    i = x.view(np.int64)
    return np.where(i < 0, np.int64(-(2**63)) - i, i)


def _ulp_dist(r: np.ndarray, a: np.ndarray) -> np.ndarray:
    kr, ka = _ordered(np.ascontiguousarray(r)), _ordered(np.ascontiguousarray(a))
    same_sign = (kr >= 0) == (ka >= 0)
    d = np.where(
        same_sign,
        np.abs(ka - np.where(same_sign, kr, 0)).astype(np.float64),
        np.abs(ka.astype(np.float64) - kr.astype(np.float64)),
    )
    return d


def _histogram(absd: np.ndarray, ulp: np.ndarray | None, is_int: bool) -> str:
    if is_int:
        edges = [(1, 1), (2, 2), (3, 4), (5, 16), (17, 256), (257, math.inf)]
        lab = ["1", "2", "3-4", "5-16", "17-256", ">256"]
        vals = absd
        name = "|difference|"
    else:
        edges = [(0, 1), (2, 4), (5, 16), (17, 256), (257, 65536), (65537, math.inf)]
        lab = ["<=1 ulp", "2-4 ulp", "5-16 ulp", "17-256 ulp", "257-65536 ulp", ">65536 ulp"]
        vals = ulp
        name = "distance"
    parts = []
    for (lo, hi), lb in zip(edges, lab, strict=True):
        c = int(np.count_nonzero((vals >= lo) & (vals <= hi)))
        if c:
            parts.append(f"{lb}: {c}")
    return f"{name} histogram: " + ", ".join(parts) if parts else ""


def _cmp_array(r: MArray, a: MArray, spec: Spec, path: str, rep: Report) -> None:
    if r.mclass != a.mclass or r.complex != a.complex:
        line = (
            f"MATLAB {r.mclass}{' complex' if r.complex else ''}, Python "
            f"{a.mclass}{' complex' if a.complex else ''}"
        )
        extra = ""
        if r.data is not None and r.size == a.size and r.numel:
            rd, ad = (
                r.data.astype(np.complex128 if (r.complex or a.complex) else np.float64),
                a.data.astype(np.complex128 if (r.complex or a.complex) else np.float64),
            )
            ndiff = int(np.count_nonzero(~((rd == ad) | (np.isnan(rd) & np.isnan(ad)))))
            extra = (
                "values are equal after conversion to double"
                if ndiff == 0
                else f"{ndiff} of {r.numel} values differ even after conversion to double"
            )
        rep.add(
            path,
            "class",
            line,
            f"size {size_text(r.size)}"
            + ("" if r.size == a.size else f" (Python {size_text(a.size)})"),
            extra,
        )
        if r.size != a.size:
            rep.add(
                path,
                "size",
                f"MATLAB {size_text(r.size)}, Python {size_text(a.size)}",
                _size_hint(r.size, a.size),
            )
        return
    if r.mclass == "char" and r.data is not None and r.is_char_row() and a.is_char_row():
        rt, at = r.text(), a.text()
        if rt != at:
            i = next(
                (k for k, (x, y) in enumerate(zip(rt, at, strict=False)) if x != y),
                min(len(rt), len(at)),
            )
            rep.add(
                path,
                "text",
                f"MATLAB {_clip(rt)!r}",
                f"Python {_clip(at)!r}",
                f"first difference at character {i + 1}",
            )
            return
    if r.size != a.size:
        lines = [
            f"MATLAB {size_text(r.size)}, Python {size_text(a.size)} (class {r.mclass})",
            _size_hint(r.size, a.size),
        ]
        if r.data is not None and a.data is not None:
            lines.extend(_size_data_hints(r.data, a.data))
        rep.add(path, "size", *lines)
        return
    if r.data is None:
        if vhash(a) != r._hash and not isinstance(spec, Exact):
            rep.undecide(
                path, f"{spec.describe()} needs MATLAB's data; the manifest has its hash only"
            )
        elif vhash(a) != r._hash:
            rep.add(
                path,
                "values",
                f"values differ from MATLAB's ({r.mclass} {size_text(r.size)}; SHA-256 mismatch)",
                HASH_ONLY,
            )
        return
    rf, af = _flat(r.data), _flat(a.data)
    neq = _neq_mask(rf, af)
    nd = int(np.count_nonzero(neq))
    if nd == 0:
        return
    lines, within = _value_details(r, a, rf, af, neq, nd, spec)
    if within:
        rep.note(path, f"{nd} of {r.numel} values differ, within {spec.describe()}")
        return
    rep.add(path, "values", *lines)


def _clip(s: str, n: int = 200) -> str:
    return s if len(s) <= n else s[:n] + "..."


def _value_details(r: MArray, a: MArray, rf, af, neq, nd, spec) -> tuple[list[str], bool]:
    n = r.numel
    idx = np.flatnonzero(neq)
    rv, av = rf[idx], af[idx]
    cls = r.mclass
    is_float = rf.dtype.kind in "fc"
    lines = [
        f"{nd} of {n} values differ ({100.0 * nd / n:.3g}%); class "
        f"{cls}{' complex' if r.complex else ''}, size {size_text(r.size)}"
    ]
    i0 = int(idx[0])
    lines.append(
        f"first difference at index {i0} (0-based, column-major) = MATLAB {_at(i0, r.size)}: "
        f"MATLAB {_fmt(rv[0], cls)}, Python {_fmt(av[0], cls)}"
    )
    if nd > 1:
        more = ", ".join(
            f"{_at(int(j), r.size)} {_fmt(rf[j], cls)} vs {_fmt(af[j], cls)}" for j in idx[1:5]
        )
        lines.append(f"next: {more}" + (" ..." if nd > 5 else ""))
    within = False
    if is_float:
        rr = np.concatenate([np.real(rv), np.imag(rv)]) if r.complex else rv
        aa = np.concatenate([np.real(av), np.imag(av)]) if r.complex else av
        rn, an = np.isnan(rr), np.isnan(aa)
        both_nan = rn & an
        one_nan = rn ^ an
        zero_sign = (rr == 0) & (aa == 0)
        if both_nan.any():
            lines.append(
                f"{int(both_nan.sum())} differ only in NaN bits (MATLAB NaN is 0xFFF8000000000000 "
                "for double; np.nan is 0x7FF8000000000000)"
            )
        if one_nan.any():
            lines.append(f"{int(one_nan.sum())} are NaN on one side only")
        if zero_sign.any():
            lines.append(f"{int(zero_sign.sum())} differ only in the sign of zero (-0 vs +0)")
        fin = ~(rn | an) & ~zero_sign
        absd = np.abs(aa.astype(np.float64) - rr.astype(np.float64))
        if fin.any():
            ulp = _ulp_dist(rr[fin], aa[fin])
            ad = absd[fin]
            k = int(np.argmax(ad))
            lines.append(f"largest |difference| {ad[k]:.6g}; largest distance {int(ulp.max())} ulp")
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = ad / np.abs(rr[fin].astype(np.float64))
            rel = rel[np.isfinite(rel)]
            if rel.size:
                lines.append(f"largest relative difference {rel.max():.6g}")
            h = _histogram(ad, ulp, False)
            if h:
                lines.append(h)
    else:
        absd = np.abs(av.astype(np.float64) - rv.astype(np.float64))
        k = int(np.argmax(absd))
        j = int(idx[k])
        lines.append(
            f"largest |difference| {absd[k]:.6g} at {_at(j, r.size)}: "
            f"MATLAB {_fmt(rf[j], cls)}, Python {_fmt(af[j], cls)}"
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = absd / np.abs(rv.astype(np.float64))
        rel = rel[np.isfinite(rel)]
        if rel.size:
            lines.append(f"largest relative difference {rel.max():.6g}")
        if cls != "logical":
            h = _histogram(absd, None, True)
            if h:
                lines.append(h)
    if isinstance(spec, Tolerance):
        ok, why = _within(spec, r, rf, af, idx, rv, av, nd, n)
        if ok:
            return lines, True
        lines.append(f"outside {spec.describe()}: {why}")
    else:
        lines.append("comparison: exact (bit for bit)")
    if n <= HINT_NUMEL:
        lines.extend(f"hint: {h}" for h in _layout_hints(r.data, a.data))
    return lines, within


def _within(t: Tolerance, r, rf, af, idx, rv, av, nd, n) -> tuple[bool, str]:
    if t.max_count is not None and nd > t.max_count:
        return False, f"{nd} values differ, more than {t.max_count}"
    if t.max_frac is not None and nd > t.max_frac * n:
        return False, f"{nd / n:.3g} of values differ, more than {t.max_frac}"
    if rf.dtype.kind in "fc":
        rn, an = np.isnan(rv), np.isnan(av)
        if (rn ^ an).any():
            return False, "NaN on one side only"
        sel = ~(rn & an)
    else:
        sel = np.ones(rv.shape, bool)
    rvd, avd = (
        rv[sel].astype(np.complex128 if r.complex else np.float64),
        av[sel].astype(np.complex128 if r.complex else np.float64),
    )
    d = np.abs(avd - rvd)
    if t.elementwise:
        allow = np.zeros(d.shape)
        if t.abs is not None:
            allow = np.maximum(allow, t.abs)
        if t.rel is not None:
            allow = np.maximum(allow, t.rel * np.abs(rvd))
        if t.rel_max is not None:
            mx = _max_abs(rf)
            allow = np.maximum(allow, t.rel_max * mx)
        ok = d <= allow
        if t.ulp is not None and rf.dtype.kind == "f":
            ok |= _ulp_dist(rv[sel], av[sel]) <= t.ulp
        if not ok.all():
            k = int(np.flatnonzero(~ok)[0])
            return False, f"{int((~ok).sum())} values beyond it; first: |difference| {d[k]:.6g}"
    if t.rel_l2 is not None:
        num = math.sqrt(float(np.sum(d.astype(np.float64) ** 2)))
        den = math.sqrt(_sum_sq(rf))
        e = num / den if den else (0.0 if num == 0 else math.inf)
        if e > t.rel_l2:
            return False, f"relative L2 error {e:.3g} > {t.rel_l2}"
    return True, ""


def _max_abs(f: np.ndarray, chunk: int = 1 << 23) -> float:
    m = 0.0
    for i in range(0, f.size, chunk):
        c = np.abs(f[i : i + chunk].astype(np.complex128 if f.dtype.kind == "c" else np.float64))
        c = c[np.isfinite(c)]
        if c.size:
            m = max(m, float(c.max()))
    return m


def _sum_sq(f: np.ndarray, chunk: int = 1 << 23) -> float:
    s = 0.0
    for i in range(0, f.size, chunk):
        c = np.abs(f[i : i + chunk].astype(np.complex128 if f.dtype.kind == "c" else np.float64))
        s += float(np.sum(c[np.isfinite(c)] ** 2))
    return s


def _mapping_hints(r: np.ndarray, a: np.ndarray) -> list[str]:
    """the mapping's mistakes that explain Python's array from MATLAB's: not mapped back,
    mapped twice (or an argument without its role), y and x swapped"""
    out = []
    n = r.ndim
    if r.size < 2:
        return out
    tests = []
    if n >= 3:
        tests += [(_perm_from_matlab(n), MISSING_MAP), (_perm_to_matlab(n), DOUBLE_MAP)]
    if n >= 2 and r.shape[0] > 1 and r.shape[1] > 1:
        tests.append(((1, 0, *range(2, n)), YX_SWAP))
    for perm, text in tests:
        t = np.transpose(r, perm)
        if t.shape != a.shape or (t.shape == r.shape and _bits_equal(t, r)):
            continue
        if _bits_equal(t, a):
            out.append(text)
    return out


def _layout_hints(r: np.ndarray, a: np.ndarray) -> list[str]:
    """explanations for values in the wrong place: the mapping, order, permutation, shift,
    offset"""
    nd = r.ndim
    if r.size < 2:
        return []
    out: list[str] = _mapping_hints(r, a)
    rc = r.reshape(-1, order="F").reshape(r.shape, order="C")
    if not _bits_equal(rc, r) and _bits_equal(a, rc):
        out.append(
            "Python has MATLAB's values in C (row-major) order: a reshape or ravel without "
            "order='F'?"
        )
    if 2 <= nd <= 4:
        for p in itertools.permutations(range(nd)):
            if p == tuple(range(nd)) or tuple(r.shape[i] for i in p) != a.shape:
                continue
            t = np.transpose(r, p)
            if not _bits_equal(t, r) and _bits_equal(a, t):
                out.append(
                    f"Python equals permute(MATLAB, [{' '.join(str(i + 1) for i in p)}]): axes in "
                    "another order"
                )
                break
    for ax in range(nd):
        for k in (1, -1, 2, -2):
            if r.shape[ax] <= abs(k) + 1:
                continue
            sa, sr = [slice(None)] * nd, [slice(None)] * nd
            sa[ax], sr[ax] = (
                (slice(k, None), slice(None, -k)) if k > 0 else (slice(None, k), slice(-k, None))
            )
            ra, rr_ = r[tuple(sa)], r[tuple(sr)]
            if not _bits_equal(ra, rr_) and _bits_equal(a[tuple(sa)], rr_):
                name = AXES[ax] if ax < 4 else str(ax + 1)
                out.append(
                    f"Python equals MATLAB shifted by {k:+d} along dimension {ax + 1} ({name}), "
                    "away from the border: an off-by-one in an index or box?"
                )
                break
    if r.dtype.kind in "iuf" and r.dtype.kind == a.dtype.kind:
        d = a.astype(np.float64) - r.astype(np.float64)
        if np.all(np.isfinite(d)) and d.size and np.all(d == d.flat[0]) and d.flat[0] != 0:
            out.append(f"Python = MATLAB {d.flat[0]:+g} everywhere")
    if r.dtype.kind in "iub" or (r.dtype.kind == "f" and not np.isnan(r).any()):
        if np.array_equal(np.sort(r, axis=None), np.sort(a, axis=None)):
            out.append("Python has the same values as MATLAB, in other positions")
    return out[:4]


def _size_data_hints(r: np.ndarray, a: np.ndarray) -> list[str]:
    out = []
    if r.ndim == a.ndim and sorted(r.shape) == sorted(a.shape):
        out += [f"hint: {h} (the values agree)" for h in _mapping_hints(r, a)]
        for p in itertools.permutations(range(r.ndim)):
            if tuple(r.shape[i] for i in p) == a.shape and _bits_equal(np.transpose(r, p), a):
                out.append(
                    f"hint: Python equals permute(MATLAB, [{' '.join(str(i + 1) for i in p)}])"
                )
                break
    if (
        r.size == a.size
        and r.dtype == a.dtype
        and _bits_equal(r.reshape(-1, order="F"), a.reshape(-1, order="F"))
    ):
        out.append("hint: same values in column-major order: a reshape with another size")
    if (
        r.size == a.size
        and r.dtype == a.dtype
        and _bits_equal(r.reshape(-1, order="F"), a.reshape(-1, order="C"))
    ):
        out.append("hint: MATLAB's column-major values are Python's in C order: order='F' missing")
    if r.ndim == a.ndim:
        d = [y - x for x, y in zip(r.shape, a.shape, strict=True)]
        nz = [i for i, x in enumerate(d) if x]
        if len(nz) == 1 and r.dtype == a.dtype:
            ax = nz[0]
            k = abs(d[ax])
            big, small = (r, a) if d[ax] < 0 else (a, r)
            for side, sl in (("last", slice(None, big.shape[ax] - k)), ("first", slice(k, None))):
                s = [slice(None)] * r.ndim
                s[ax] = sl
                if _bits_equal(big[tuple(s)], small):
                    who = "Python" if d[ax] < 0 else "MATLAB"
                    out.append(
                        f"hint: {who} lacks the {side} {k} element(s) along dimension {ax + 1}; "
                        "otherwise equal"
                    )
                    break
    return out


# ---------------------------------------------------------------------------------------------
# MAT-files (the Python side; MATLAB decodes its own files in the runner)
# ---------------------------------------------------------------------------------------------


def read_mat(path: str | os.PathLike) -> MStruct:
    """the variables of a MAT-file as a 1x1 struct (v7.3 through h5py, older through SciPy)"""
    with open(path, "rb") as f:
        head = f.read(128)
    if b"MATLAB 7.3 MAT-file" in head or head[:8] == b"\x89HDF\r\n\x1a\n":
        return _read_mat73(path)
    from scipy.io import loadmat

    d = loadmat(
        path, squeeze_me=False, struct_as_record=True, chars_as_strings=False, mat_dtype=True
    )
    names = [k for k in d if not k.startswith("__")]
    return MStruct((1, 1), names, [{k: _from_scipy(d[k]) for k in names}])


def _from_scipy(x: Any) -> MValue:
    if _is_sparse(x):
        return _sparse_from_scipy(x)
    if isinstance(x, np.ndarray):
        if x.dtype.names:
            fields = list(x.dtype.names)
            els = x.reshape(-1, order="F")
            return MStruct(
                msize(x.shape) if x.ndim >= 2 else (1, x.size),
                fields,
                [{f: _from_scipy(el[f]) for f in fields} for el in els],
            )
        if x.dtype.kind == "O":
            return MCell(
                msize(x.shape) if x.ndim >= 2 else (1, x.size),
                [_from_scipy(e) for e in x.reshape(-1, order="F")],
            )
        if x.dtype.kind == "U":
            if x.ndim < 2:
                x = x.reshape(1, -1) if x.size else x.reshape(0, 0)
            if x.dtype.itemsize == 4:
                cp = x.view(np.uint32) if x.size else np.zeros(x.shape, np.uint32)
                if np.all(cp < 0x10000):
                    return MArray("char", cp.astype("<u2"))
            return mchar("".join(x.reshape(-1, order="F")))
        if x.ndim < 2:
            x = x.reshape(1, -1)
        return _array_from_numpy(x)
    name = type(x).__name__
    if name == "MatlabFunction":
        return MFunc("<function handle>")
    if name in ("MatlabOpaque", "MatlabObject"):
        cls = getattr(x, "classname", None) or (x.dtype.names and "object") or "object"
        return MObject(str(cls), None, "(MATLAB object, not decoded)")
    return to_mvalue(x)


def _read_mat73(path) -> MStruct:
    import h5py

    with h5py.File(path, "r") as f:
        names = [k for k in f.keys() if not k.startswith("#")]
        return MStruct((1, 1), names, [{k: _h5_value(f[k], f) for k in names}])


def _attr(obj, name: str) -> str | None:
    v = obj.attrs.get(name)
    if v is None:
        return None
    return v.decode() if isinstance(v, bytes) else (str(v) if not isinstance(v, np.ndarray) else v)


def _h5_value(obj, f) -> MValue:
    import h5py

    cls = _attr(obj, "MATLAB_class")
    if isinstance(obj, h5py.Group):
        if "MATLAB_sparse" in obj.attrs:
            m = int(obj.attrs["MATLAB_sparse"])
            jc = obj["jc"][()].astype(np.int64)
            n = jc.size - 1
            rows = obj["ir"][()].astype(np.int64) if "ir" in obj else np.zeros(0, np.int64)
            data = obj["data"][()] if "data" in obj else np.zeros(0)
            if data.dtype.names:
                data = data["real"] + 1j * data["imag"]
            cols = np.repeat(np.arange(n, dtype=np.int64), np.diff(jc))
            if cls == "logical":
                data = data.astype(bool)
            return MSparse(cls or "double", (m, n), rows, cols, data)
        if cls == "function_handle":
            return MFunc("<function handle>")
        fields = obj.attrs.get("MATLAB_fields")
        names = (
            [b"".join(x).decode() if isinstance(x, np.ndarray) else str(x) for x in fields]
            if fields is not None
            else list(obj.keys())
        )
        if not names:
            return MStruct((1, 1), [], [{}])
        first = obj[names[0]]
        if (
            isinstance(first, h5py.Dataset)
            and first.dtype == h5py.ref_dtype
            and _attr(first, "MATLAB_class") is None
        ):
            shape = tuple(reversed(first.shape))
            refs = {k: obj[k][()].T.reshape(-1, order="F") for k in names}
            n = math.prod(shape)
            items = [{k: _h5_value(f[refs[k][i]], f) for k in names} for i in range(n)]
            return MStruct(shape, names, items)
        return MStruct((1, 1), names, [{k: _h5_value(obj[k], f) for k in names}])
    if _attr(obj, "MATLAB_object_decode") is not None or cls in (
        "string",
        "inputParser",
        "datetime",
        "table",
    ):
        return MObject(cls or "object", None, "(MATLAB object, not decoded)")
    data = obj[()]
    if obj.attrs.get("MATLAB_empty"):
        size = tuple(int(n) for n in np.asarray(data).reshape(-1))
        if cls == "cell":
            return MCell(size, [])
        if cls == "struct":
            return MStruct(size, [], [])
        return MArray(cls, np.zeros(size, _DTYPE[cls]))
    if cls == "cell":
        refs = data.T.reshape(-1, order="F") if data.ndim else data.reshape(-1)
        return MCell(tuple(reversed(data.shape)), [_h5_value(f[r], f) for r in refs])
    if isinstance(data, np.ndarray) and data.dtype.names and "real" in data.dtype.names:
        data = data["real"] + 1j * data["imag"]
        data = data.astype(_CDTYPE[cls])
    a = np.asarray(data).T
    if cls == "logical":
        a = a.astype(bool)
    elif cls == "char":
        a = a.astype("<u2")
    elif cls in _DTYPE:
        a = a.astype(_CDTYPE[cls] if a.dtype.kind == "c" else _DTYPE[cls], copy=False)
    else:
        return MObject(cls or "unknown", None, "(not decoded)")
    if a.ndim < 2:
        a = a.reshape(1, -1)
    return MArray(cls, a)


# ---------------------------------------------------------------------------------------------
# file trees
# ---------------------------------------------------------------------------------------------

TEXT_EXT = {
    ".txt",
    ".log",
    ".csv",
    ".m",
    ".sh",
    ".xml",
    ".out",
    ".err",
    ".md",
    ".py",
    ".tsv",
    ".ini",
    ".cfg",
}


def file_kind(rel: str, is_dir: bool = False) -> str:
    low = rel.lower()
    if is_dir:
        return "zarr" if low.endswith(".zarr") else "dir"
    ext = os.path.splitext(low)[1]
    if ext in (".tif", ".tiff"):
        return "tiff"
    if ext == ".mat":
        return "mat"
    if ext == ".json":
        return "json"
    if ext == ".png":
        return "png"
    if ext == ".pdf":
        return "pdf"
    if ext in TEXT_EXT:
        return "text"
    return "bytes"


def normalize_name(rel: str, masks: Masks | None = None) -> str:
    """the name a file is matched by: uuids and time stamps masked, forward slashes"""
    rel = rel.replace("\\", "/")
    return (masks or Masks()).text(rel)


@dataclasses.dataclass
class FileItem:
    rel: str
    kind: str
    path: Path | None = None  # the file itself (None: known by descriptor only)
    desc: dict | None = None  # descriptor (manifests)
    decoded: Path | None = None  # MATLAB's own decoding of a .mat file (value file)


@dataclasses.dataclass
class FileTree:
    items: dict[str, list[FileItem]]
    deleted: list[str] = dataclasses.field(default_factory=list)

    def names(self) -> list[str]:
        return sorted(self.items)


def scan_tree(
    root: str | os.PathLike,
    rels: Iterable[str] | None = None,
    decoded_root: str | os.PathLike | None = None,
    deleted: Iterable[str] = (),
    masks: Masks | None = None,
) -> FileTree:
    """the files under root, or only ``rels`` (files, empty folders, zarr folders)"""
    root = Path(root)
    if rels is None:
        rels = list_tree(root)
    items: dict[str, list[FileItem]] = {}
    for rel in sorted(rels):
        p = root / rel
        k = file_kind(rel, p.is_dir())
        dec = None
        if (
            k == "mat"
            and decoded_root is not None
            and (Path(decoded_root) / (rel + ".pkv")).is_file()
        ):
            dec = Path(decoded_root) / (rel + ".pkv")
        items.setdefault(normalize_name(rel, masks), []).append(FileItem(rel, k, p, None, dec))
    return FileTree(items, list(deleted))


def list_tree(root: str | os.PathLike) -> list[str]:
    """relative paths of the files, the empty folders and the zarr folders under root"""
    root = Path(root)
    out = []
    for dp, dns, fns in os.walk(root):
        d = Path(dp)
        rel = d.relative_to(root).as_posix()
        if rel != "." and d.name.lower().endswith(".zarr"):
            out.append(rel)
            dns[:] = []
            continue
        if rel != "." and not dns and not fns:
            out.append(rel)
        out.extend((d / f).relative_to(root).as_posix() for f in fns)
    return sorted(out)


def tree_from_descriptor(d: dict) -> FileTree:
    items: dict[str, list[FileItem]] = {}
    for name, entries in d.get("files", {}).items():
        items[name] = [FileItem(e.get("rel", name), e["k"], None, e) for e in entries]
    return FileTree(items, list(d.get("deleted", [])))


@dataclasses.dataclass(frozen=True)
class TreeSpec:
    """how output files are compared

    values: spec for decoded voxels and .mat contents; per_file: {glob: spec} for matching
    names; ignore: regexes of names not compared; ignore_vars: .mat variables not compared;
    ignore_fields: struct fields and JSON keys not compared, at any depth; tiff_bytes: tiff
    files must be equal byte for byte (else voxels and tags only); zarr_chunks: 'note' reports
    differing chunk-file sets as a note (zarr is compared by decoded data), 'strict' as a
    problem; text_format: JSON and text must be equal as text (else by content);
    masks: extra (regex, replacement) pairs for text, names and values.
    """

    values: Spec = EXACT
    per_file: tuple[tuple[str, Spec], ...] = ()
    ignore: tuple[str, ...] = ()
    ignore_vars: tuple[str, ...] = ("t", "t1")
    ignore_fields: tuple[str, ...] = ()
    tiff_bytes: bool = True
    zarr_chunks: str = "note"
    text_format: bool = True
    masks: tuple[tuple[str, str], ...] = ()

    def spec_for(self, name: str) -> Spec:
        for g, s in self.per_file:
            if fnmatch.fnmatch(name, g):
                return s
        return self.values

    def key(self) -> str:
        return hashlib.sha256(
            repr((self.ignore, self.ignore_vars, self.ignore_fields, self.masks)).encode()
        ).hexdigest()[:16]


def tree(
    *,
    values: Spec = EXACT,
    per_file: Mapping[str, Spec] | None = None,
    ignore: Sequence[str] = (),
    ignore_vars: Sequence[str] = ("t", "t1"),
    ignore_fields: Sequence[str] = (),
    tiff_bytes: bool = True,
    zarr_chunks: str = "note",
    text_format: bool = True,
    masks: Sequence[tuple[str, str]] = (),
) -> TreeSpec:
    return TreeSpec(
        values,
        tuple((per_file or {}).items()),
        tuple(ignore),
        tuple(ignore_vars),
        tuple(ignore_fields),
        tiff_bytes,
        zarr_chunks,
        text_format,
        tuple(masks),
    )


def _read_tiff(path: Path) -> MArray:
    try:
        import cpptiff

        a = cpptiff.read_tiff(str(path))
    except ImportError:
        import tifffile

        a = tifffile.imread(str(path))
    a = np.asarray(a)
    if a.ndim == 3:
        a = np.transpose(a, (1, 2, 0))
    return _array_from_numpy(a if a.ndim >= 2 else a.reshape(1, -1))


def _read_zarr(path: Path) -> MArray:
    import cppzarr

    a = np.asarray(cppzarr.read_zarr(str(path)))
    return _array_from_numpy(a if a.ndim >= 2 else a.reshape(1, -1))


ZARR_META = (".zarray", ".zattrs", ".zgroup", "zarr.json")


def _zarr_parts(path: Path) -> tuple[dict[str, bytes], list[str]]:
    meta, chunks = {}, []
    for dp, _, fns in os.walk(path):
        for fn in fns:
            p = Path(dp) / fn
            rel = p.relative_to(path).as_posix()
            if fn in ZARR_META:
                meta[rel] = p.read_bytes()
            else:
                chunks.append(rel)
    return meta, sorted(chunks)


def _png_payload(data: bytes) -> bytes:
    i, out = 8, [data[:8]]
    while i + 8 <= len(data):
        n = struct.unpack(">I", data[i : i + 4])[0]
        t = data[i + 4 : i + 8]
        body = data[i + 8 : i + 8 + n]
        if not (t == b"tIME" or (t == b"tEXt" and body.startswith(b"Creation Time"))):
            out.append(t + body)
        i += 12 + n
    return b"".join(out)


class _Side:
    """lazy access to one file of one side, with that side's masks"""

    def __init__(self, item: FileItem, masks: Masks, spec: TreeSpec):
        self.it, self.masks, self.spec = item, masks, spec
        self._cache: dict[str, Any] = {}

    def _get(self, key: str, fn):
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    @property
    def has_file(self) -> bool:
        return self.it.path is not None

    def d(self, key: str, default=None):
        return (self.it.desc or {}).get(key, default)

    def sha(self) -> str:
        if not self.has_file:
            return self.d("sha")
        return self._get("sha", lambda: sha256_file(self.it.path))

    def text(self) -> str | None:
        if not self.has_file:
            return self.d("text")
        return self._get(
            "text", lambda: self.masks.text(self.it.path.read_bytes().decode("utf-8", "replace"))
        )

    def text_sha(self) -> str:
        if not self.has_file:
            return self.d("sha")
        return self._get("tsha", lambda: hashlib.sha256(self.text().encode("utf-8")).hexdigest())

    def png_sha(self) -> str:
        if not self.has_file:
            return self.d("sha")
        return self._get(
            "psha", lambda: hashlib.sha256(_png_payload(self.it.path.read_bytes())).hexdigest()
        )

    def pdf_sha(self) -> str:
        """PDF bytes without the creation and modification dates (as isolation_compare)"""
        if not self.has_file:
            return self.d("sha")
        pat = rb"/(CreationDate|ModDate) \([^)]*\)"
        return self._get(
            "dsha", lambda: hashlib.sha256(re.sub(pat, b"", self.it.path.read_bytes())).hexdigest()
        )

    def voxels(self) -> MValue:
        if not self.has_file:
            return from_descriptor(self.d("data"))
        fn = _read_tiff if self.it.kind == "tiff" else _read_zarr
        return self._get("vox", lambda: fn(self.it.path))

    def voxels_to_matlab(self) -> MValue:
        """the voxels of a zarr in Python's layout ((..., z, y, x)) as MATLAB's (y, x, z, ...)"""
        if not self.has_file:
            return self.voxels()

        def load():
            import cppzarr

            x = np.asarray(cppzarr.read_zarr(str(self.it.path)))
            if x.ndim >= 2:
                x = np.transpose(x, _perm_to_matlab(x.ndim))
            return _array_from_numpy(x if x.ndim >= 2 else x.reshape(1, -1))

        return self._get("voxm", load)

    def mat(self) -> MValue:
        if not self.has_file:
            return from_descriptor(self.d("vars"))

        def load():
            v = read_pkv(self.it.decoded) if self.it.decoded is not None else read_mat(self.it.path)
            if isinstance(v, MStruct) and self.spec.ignore_vars:
                keep = [f for f in v.fields if f not in self.spec.ignore_vars]
                v = MStruct(v.size, keep, [{f: el[f] for f in keep} for el in v.items])
            return mask_value(drop_fields(v, self.spec.ignore_fields), self.masks)

        return self._get("mat", load)

    def mat_error(self) -> str | None:
        """why the file is not a readable MAT-file (an empty flag file named .mat), or None"""
        if not self.has_file:
            return self.d("unreadable")

        def check():
            try:
                self.mat()
            except Exception as e:  # any reader's failure: the file is compared as bytes
                return f"{type(e).__name__}: {str(e)[:200]}"
            return None

        return self._get("materr", check)

    def zarr_meta(self) -> dict[str, str]:
        if not self.has_file:
            return self.d("meta", {})
        return self._get(
            "zmeta",
            lambda: {
                k: self.masks.text(v.decode("utf-8", "replace"))
                for k, v in _zarr_parts(self.it.path)[0].items()
            },
        )

    def zarr_chunks(self) -> tuple[int, str]:
        if not self.has_file:
            return self.d("nchunks", 0), self.d("chunks_sha", "")
        ch = self._get("zch", lambda: _zarr_parts(self.it.path)[1])
        return len(ch), hashlib.sha256("\n".join(ch).encode()).hexdigest()

    def chunk_names(self) -> list[str] | None:
        return self._get("zch", lambda: _zarr_parts(self.it.path)[1]) if self.has_file else None

    def describe(self) -> dict:
        k = self.it.kind
        d: dict[str, Any] = {"k": k, "rel": normalize_name(self.it.rel, self.masks)}
        if k == "tiff":
            d["sha"] = self.sha()
            d["data"] = describe(self.voxels())
        elif k == "zarr":
            d["meta"] = self.zarr_meta()
            d["nchunks"], d["chunks_sha"] = self.zarr_chunks()
            d["data"] = describe(self.voxels())
        elif k == "mat":
            err = self.mat_error()
            if err is None:
                d["vars"] = describe(self.mat())
            else:
                d["sha"], d["unreadable"] = self.sha(), err
        elif k in ("json", "text"):
            d["sha"] = self.text_sha()
            t = self.text()
            if t is not None and len(t.encode("utf-8")) <= INLINE_TEXT:
                d["text"] = t
        elif k == "png":
            d["sha"] = self.png_sha()
        elif k == "pdf":
            d["sha"] = self.pdf_sha()
        elif k == "bytes":
            d["sha"] = self.sha()
        return d


def describe_tree(t: FileTree, spec: TreeSpec, masks: Masks) -> dict:
    files = {}
    for name in t.names():
        if any(re.search(p, name) for p in spec.ignore):
            continue
        files[name] = [_Side(it, masks, spec).describe() for it in t.items[name]]
    return {"files": files, "deleted": sorted(normalize_name(x, masks) for x in t.deleted)}


def compare_trees(
    ref: FileTree,
    act: FileTree,
    spec: TreeSpec | None = None,
    *,
    ref_masks: Masks | None = None,
    act_masks: Masks | None = None,
    path: str = "files",
) -> Report:
    """compare the files Python wrote with those MATLAB wrote"""
    spec = spec or TreeSpec()
    rep = Report()
    rm = (ref_masks or Masks()).plus(*spec.masks)
    am = (act_masks or Masks()).plus(*spec.masks)
    rnames = {n for n in ref.items if not any(re.search(p, n) for p in spec.ignore)}
    anames = {n for n in act.items if not any(re.search(p, n) for p in spec.ignore)}
    for n in sorted(rnames - anames):
        kinds = ", ".join(sorted({it.kind for it in ref.items[n]}))
        rep.add(f"{path}/{n}", "missing", f"MATLAB wrote this {kinds}, Python did not")
    for n in sorted(anames - rnames):
        rep.add(f"{path}/{n}", "extra", "Python wrote this, MATLAB did not")
    rdel = sorted(normalize_name(x, rm) for x in ref.deleted)
    adel = sorted(normalize_name(x, am) for x in act.deleted)
    if rdel != adel:
        rep.add(
            path, "deleted", f"files removed from the work folder: MATLAB {rdel}, Python {adel}"
        )
    for n in sorted(rnames & anames):
        ri, ai = ref.items[n], act.items[n]
        if len(ri) != len(ai):
            rep.add(
                f"{path}/{n}",
                "count",
                f"{len(ri)} entries match this name in MATLAB, {len(ai)} in Python",
            )
            continue
        for r_it, a_it in zip(ri, ai, strict=True):
            _cmp_file(
                _Side(r_it, rm, spec),
                _Side(a_it, am, spec),
                spec.spec_for(n),
                spec,
                f"{path}/{n}",
                rep,
            )
    return rep


def _cmp_file(r: _Side, a: _Side, vspec: Spec, spec: TreeSpec, path: str, rep: Report) -> None:
    if isinstance(vspec, Ignore):
        return
    rk, ak = r.it.kind, a.it.kind
    if rk != ak:
        rep.add(path, "kind", f"MATLAB wrote a {rk}, Python a {ak}")
        return
    if rk == "dir":
        return
    if rk == "tiff":
        if r.sha() == a.sha():
            return
        sub = Report()
        _cmp(r.voxels(), a.voxels(), vspec, path + " (voxels)", sub)
        rep.extend(sub)
        if sub.ok and not sub.undecided:
            if spec.tiff_bytes:
                rep.add(
                    path,
                    "layout",
                    "voxels are equal but the tiff bytes differ (layout, tags or compression)",
                    f"MATLAB {r.sha()[:16]}, Python {a.sha()[:16]}",
                )
            else:
                rep.note(path, "tiff bytes differ, voxels equal")
        return
    if rk == "zarr":
        rmeta, ameta = r.zarr_meta(), a.zarr_meta()
        rz, az = _zarray_of(rmeta), _zarray_of(ameta)
        if rz is not None and az is not None and rz.get("order") == "F" and az.get("order") == "C":
            _cmp_zarr_layouts(r, a, rz, az, vspec, spec, path, rep)
            return
        for k in sorted(set(rmeta) | set(ameta)):
            if rmeta.get(k) != ameta.get(k):
                rep.add(f"{path}/{k}", "meta", _text_diff(rmeta.get(k), ameta.get(k)))
        (rn, rh), (an, ah) = r.zarr_chunks(), a.zarr_chunks()
        if rh != ah:
            line = f"chunk files differ: MATLAB {rn}, Python {an}"
            rc, ac = r.chunk_names(), a.chunk_names()
            if rc is not None and ac is not None:
                only_r, only_a = sorted(set(rc) - set(ac)), sorted(set(ac) - set(rc))
                line += (
                    f" (only MATLAB: {only_r[:5]}{'...' if len(only_r) > 5 else ''}; only Python: "
                    f"{only_a[:5]}{'...' if len(only_a) > 5 else ''})"
                )
            if spec.zarr_chunks == "strict":
                rep.add(path, "layout", line)
            else:
                rep.note(path, line + "; zarr is compared by decoded data")
        _cmp(r.voxels(), a.voxels(), vspec, path + " (voxels)", rep)
        return
    if rk == "mat":
        rerr, aerr = r.mat_error(), a.mat_error()
        if rerr is None and aerr is None:
            _cmp(r.mat(), a.mat(), vspec, path, rep)
        elif rerr is not None and aerr is not None:
            # neither side is a MAT-file (an empty flag file named .mat): byte for byte
            if r.sha() != a.sha():
                rep.add(
                    path,
                    "bytes",
                    f"not a readable MAT-file on either side ({aerr}); the bytes differ: MATLAB "
                    f"{r.sha()[:16]}, Python {a.sha()[:16]}",
                )
            else:
                rep.note(path, f"not a readable MAT-file ({aerr}); compared byte for byte")
        else:
            who, err = ("MATLAB's", rerr) if rerr is not None else ("Python's", aerr)
            rep.add(path, "mat", f"{who} file is not a readable MAT-file ({err}); the other's is")
        return
    if rk in ("json", "text"):
        if r.text_sha() == a.text_sha():
            return
        rt, at = r.text(), a.text()
        if rk == "json" and rt is not None and at is not None:
            try:
                rj, aj = json.loads(rt), json.loads(at)
            except ValueError:
                rj = aj = None
            if rj is not None or aj is not None:
                for f in spec.ignore_fields:
                    rj, aj = _drop_json(rj, f), _drop_json(aj, f)
                if rj == aj:
                    if spec.text_format:
                        rep.add(
                            path,
                            "format",
                            "JSON content is equal but the text differs",
                            _text_diff(rt, at),
                        )
                    return
                rep.add(path, "values", _json_diff(rj, aj, ""), _text_diff(rt, at))
                return
        if rt is None or at is None:
            rep.add(
                path,
                "text",
                f"text differs from MATLAB's (SHA-256 mismatch; the manifest keeps texts up to "
                f"{INLINE_TEXT} bytes, else hashes: the readable diff is in the golden store, run "
                "this test on the cluster)",
            )
        else:
            rep.add(path, "text", _text_diff(rt, at))
        return
    if rk == "png":
        if r.png_sha() != a.png_sha():
            rep.add(path, "bytes", "png content differs (time stamps not compared)")
        return
    if rk == "pdf":
        if r.pdf_sha() != a.pdf_sha():
            rep.add(path, "bytes", "pdf differs (creation and modification dates not compared)")
        return
    if r.sha() != a.sha():
        rep.add(path, "bytes", f"bytes differ: MATLAB {r.sha()[:16]}, Python {a.sha()[:16]}")


def _zarray_of(meta: Mapping[str, str]) -> dict | None:
    """the array's .zarray (zarr v2) as a dict, or None"""
    t = meta.get(".zarray")
    if t is None:
        return None
    try:
        d = json.loads(t)
    except ValueError:
        return None
    return d if isinstance(d, dict) and isinstance(d.get("shape"), list) else None


def zarray_to_matlab(z: Mapping[str, Any]) -> dict:
    """a .zarray of Python's layout (order C, shape and chunks (..., z, y, x)) as MATLAB's
    (order F, (y, x, z, ...)): what MATLAB writes for the same array"""
    out = dict(z)
    n = len(z["shape"])
    perm = _perm_to_matlab(n)
    out["shape"] = [z["shape"][i] for i in perm]
    if isinstance(z.get("chunks"), list) and len(z["chunks"]) == n:
        out["chunks"] = [z["chunks"][i] for i in perm]
    out["order"] = "F"
    return out


def _chunk_names_to_matlab(names: list[str], n: int, sep: str) -> list[str]:
    perm = _perm_to_matlab(n)
    out = []
    for nm in names:
        parts = nm.split("/" if sep == "/" else ".")
        if len(parts) == n and all(x.isdigit() for x in parts):
            nm = (sep or ".").join(parts[i] for i in perm)
        out.append(nm)
    return sorted(out)


def _cmp_zarr_layouts(
    r: _Side, a: _Side, rz: dict, az: dict, vspec: Spec, spec: TreeSpec, path: str, rep: Report
) -> None:
    """MATLAB's zarr (order F, (y, x, z)) against Python's native one (order C, (z, y, x)):
    .zarray modulo the permutation of shape and chunks and the order flag, the voxels after
    the mapping (core.arrays.to_matlab)"""
    rep.note(
        path,
        "Python's zarr is in Python's layout (order C, (z, y, x)), MATLAB's in its own (order F, "
        "(y, x, z)): .zarray and voxels compared after the mapping",
    )
    mz = zarray_to_matlab(az)
    for k in sorted(set(rz) | set(mz)):
        if rz.get(k) == mz.get(k):
            continue
        lines = [f"{k}: MATLAB {json.dumps(rz.get(k))}, Python {json.dumps(az.get(k))}"]
        if k in ("shape", "chunks"):
            lines[0] += f" (in MATLAB's axes {json.dumps(mz.get(k))})"
            if az.get(k) == rz.get(k):
                lines.append(
                    f"hint: Python wrote MATLAB's (y, x, z) {k} with order C: a volume written "
                    "without its mapping (from_matlab)?"
                )
        rep.add(f"{path}/.zarray", "meta", *lines)
    rmeta, ameta = r.zarr_meta(), a.zarr_meta()
    for k in sorted((set(rmeta) | set(ameta)) - {".zarray"}):
        if rmeta.get(k) != ameta.get(k):
            rep.add(f"{path}/{k}", "meta", _text_diff(rmeta.get(k), ameta.get(k)))
    n = len(az["shape"])
    rc, ac = r.chunk_names(), a.chunk_names()
    if rc is not None and ac is not None:
        sep = az.get("dimension_separator", ".") or "."
        mapped = _chunk_names_to_matlab(ac, n, sep)
        rsep = rz.get("dimension_separator", ".") or "."
        rmapped = sorted(nm if sep == rsep else nm.replace(rsep, sep) for nm in rc)
        if mapped != rmapped:
            line = (
                f"chunk files differ after the mapping: MATLAB {len(rc)}, Python {len(ac)}"
                f" (only MATLAB: {sorted(set(rmapped) - set(mapped))[:5]}; only Python: "
                f"{sorted(set(mapped) - set(rmapped))[:5]})"
            )
            if spec.zarr_chunks == "strict":
                rep.add(path, "layout", line)
            else:
                rep.note(path, line + "; zarr is compared by decoded data")
    av = a.voxels_to_matlab()
    with layout("vol"):
        _cmp(r.voxels(), av, vspec, path + " (voxels)", rep)


def _drop_json(x: Any, name: str) -> Any:
    if isinstance(x, dict):
        return {k: _drop_json(v, name) for k, v in x.items() if k != name}
    if isinstance(x, list):
        return [_drop_json(v, name) for v in x]
    return x


def _json_diff(r: Any, a: Any, path: str) -> str:
    if type(r) is not type(a):
        return f"JSON {path or 'root'}: MATLAB {json.dumps(r)[:120]}, Python {json.dumps(a)[:120]}"
    if isinstance(r, dict):
        if set(r) != set(a):
            return (
                f"JSON {path or 'root'}: keys only in MATLAB {sorted(set(r) - set(a))}, only in "
                f"Python {sorted(set(a) - set(r))}"
            )
        for k in r:
            if r[k] != a[k]:
                return _json_diff(r[k], a[k], f"{path}.{k}")
    if isinstance(r, list):
        if len(r) != len(a):
            return f"JSON {path or 'root'}: MATLAB has {len(r)} items, Python {len(a)}"
        for i, (x, y) in enumerate(zip(r, a, strict=True)):
            if x != y:
                return _json_diff(x, y, f"{path}[{i}]")
    return f"JSON {path or 'root'}: MATLAB {json.dumps(r)[:120]}, Python {json.dumps(a)[:120]}"


def _text_diff(rt: str | None, at: str | None) -> str:
    if rt is None or at is None:
        return f"present only in {'Python' if rt is None else 'MATLAB'}"
    rl, al = rt.split("\n"), at.split("\n")
    for i in range(max(len(rl), len(al))):
        x = rl[i] if i < len(rl) else None
        y = al[i] if i < len(al) else None
        if x != y:
            return (
                f"first difference on line {i + 1}: MATLAB {_clip(x or '', 160)!r}, Python "
                f"{_clip(y or '', 160)!r}"
                + (f" ({len(rl)} vs {len(al)} lines)" if len(rl) != len(al) else "")
            )
    return "texts differ"
