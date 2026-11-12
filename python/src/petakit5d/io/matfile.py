"""MAT-files as MATLAB R2024b reads and writes them.

``read_mat`` reads v7.3 files (HDF5, through h5py) and v6/v7 files (MAT 5 format, our own
reader); ``write_mat`` writes v7.3, v7 (compressed MAT 5) and v6 (uncompressed MAT 5) files
that MATLAB's ``load`` returns with the classes, sizes and values given; ``save_ascii`` writes
what ``save(fn, 'x', '-ascii')`` writes.

Python to MATLAB (``write_mat``):

=========================================  ==================================================
Python                                     MATLAB
=========================================  ==================================================
``bool``, ``numpy.bool_``                  logical 1x1
``int``                                    double 1x1 (an int a double can't hold exactly is
                                           an error)
``float``                                  double 1x1
``complex``                                complex double 1x1
``str``, ``numpy.str_``                    char 1xN, UTF-16 code units (``''`` is 0x0)
``bytes``                                  uint8 1xN (``b''`` is 0x0)
``None``                                   ``[]`` (double 0x0)
``list``, ``tuple``                        cell 1xN (empty: ``{}``, 0x0)
``dict``                                   struct 1x1 (keys are the fields, in order; ``{}``
                                           is ``struct()``)
``StructArray``                            struct array of its shape (empty ones too)
numpy scalar, 0-d array                    1x1 of its class
``numpy.ndarray``                          its class (below) and shape; 1-D arrays are rows
                                           (1xN); trailing singleton dimensions after the
                                           second are dropped, as in MATLAB
  float64, float32                         double, single
  int8 ... uint64                          int8 ... uint64
  bool                                     logical
  complex128, complex64                    complex double, complex single
  ``<U1`` (one character each)             char array of that shape
  object                                   cell array of that shape
``scipy.sparse`` matrix or array           sparse double, logical or complex double
=========================================  ==================================================

Anything else (sets, float16, datetimes, wider string arrays, ``MatlabObject``) is a TypeError
or ValueError. Variable and field names must be MATLAB names (a letter, then letters, digits
and underscores, 63 at most).

MATLAB to Python (``read_mat``), chosen so that ``write_mat`` of what ``read_mat`` returns gives
the same MATLAB values back:

- numeric and logical: 1x1 is a numpy scalar of the class's dtype (``numpy.float64``,
  ``numpy.uint16``, ``numpy.bool_``, ...); any other size an F-ordered ``numpy.ndarray`` of
  MATLAB's size (at least 2-D, ``[]`` is shape ``(0, 0)``).
- char: 1xN (N >= 1) and 0x0 are ``str`` (UTF-16 surrogate pairs joined); other sizes a
  ``<U1`` array of the size.
- cell: 1xN (N >= 1) and 0x0 are a ``list``; other sizes an object array of the size.
- struct: 1x1 is a ``dict`` (fields in MATLAB's order); other sizes a ``StructArray``.
- sparse: ``scipy.sparse.csc_matrix`` (float64, bool or complex128).
- string: a scalar is a ``str``, other sizes an object array of ``str`` (``None`` for
  ``<missing>``). Written back it becomes char (or a cell), not string.
- objects (``inputParser``, function handles, ``datetime``, ...): a ``MatlabObject`` that names
  the class; their contents are not read, and ``write_mat`` refuses them.
- complex integer arrays (which numpy lacks) are read as complex128.
"""

from __future__ import annotations

import math
import os
import re
import struct
import sys
import time
import warnings
import zlib
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = ["MatlabObject", "StructArray", "mat_version", "read_mat", "save_ascii", "write_mat"]

PathLike = str | os.PathLike

_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,62}")
_NUMERIC = {
    "double": np.float64,
    "single": np.float32,
    "int8": np.int8,
    "uint8": np.uint8,
    "int16": np.int16,
    "uint16": np.uint16,
    "int32": np.int32,
    "uint32": np.uint32,
    "int64": np.int64,
    "uint64": np.uint64,
}
_DTYPE_CLASS = {np.dtype(t): c for c, t in _NUMERIC.items()}
_DTYPE_CLASS[np.dtype(np.bool_)] = "logical"
_DTYPE_CLASS[np.dtype(np.complex128)] = "double"
_DTYPE_CLASS[np.dtype(np.complex64)] = "single"


class MatlabObject:
    """A MATLAB object read_mat does not decode (an inputParser, a function handle, ...)."""

    __slots__ = ("class_name", "note")

    def __init__(self, class_name: str, note: str = "contents not read") -> None:
        self.class_name = class_name
        self.note = note

    def __repr__(self) -> str:
        return f"MatlabObject({self.class_name!r}: {self.note})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, MatlabObject) and other.class_name == self.class_name

    def __hash__(self) -> int:
        return hash(self.class_name)


class StructArray:
    """A MATLAB struct array of any size but 1x1 (a 1x1 struct is a dict).

    Elements are dicts with the same fields, kept in column-major order. ``s[i]`` is the
    element at linear index i (0-based, column-major), ``s[i, j]`` the one at those 0-based
    subscripts; ``len(s)`` is the number of elements and iteration goes column-major.
    """

    __slots__ = ("_elements", "_shape", "fields")

    def __init__(
        self, shape: Sequence[int], fields: Sequence[str], elements: Iterable[Mapping] = ()
    ):
        self._shape = tuple(int(n) for n in shape)
        self.fields = tuple(fields)
        els = [dict(e) for e in elements]
        if len(els) != math.prod(self._shape):
            raise ValueError(f"StructArray: {len(els)} elements for size {self._shape}")
        for e in els:
            if set(e) != set(self.fields) or len(e) != len(self.fields):
                raise ValueError(f"StructArray: element fields {list(e)} != {list(self.fields)}")
        self._elements = [{f: e[f] for f in self.fields} for e in els]

    @classmethod
    def from_list(
        cls,
        rows: Sequence[Mapping],
        shape: Sequence[int] | None = None,
        fields: Sequence[str] | None = None,
    ) -> StructArray:
        """a struct array from dicts in column-major order (default size 1xN)"""
        if fields is None:
            fields = tuple(rows[0]) if rows else ()
        return cls((1, len(rows)) if shape is None else shape, fields, rows)

    @classmethod
    def empty(cls, fields: Sequence[str] = (), shape: Sequence[int] = (0, 0)) -> StructArray:
        return cls(shape, fields, [])

    @property
    def shape(self) -> tuple[int, ...]:
        return self._shape

    @property
    def size(self) -> int:
        return len(self._elements)

    @property
    def ndim(self) -> int:
        return len(self._shape)

    def __len__(self) -> int:
        return len(self._elements)

    def __iter__(self) -> Iterator[dict]:
        return iter(self._elements)

    def ravel(self) -> list[dict]:
        """the elements, column-major"""
        return list(self._elements)

    def __getitem__(self, index: int | tuple[int, ...]) -> dict:
        if isinstance(index, tuple):
            index = int(np.ravel_multi_index(index, self._shape, order="F"))
        return self._elements[index]

    def field(self, name: str) -> list:
        """one field of every element, column-major (MATLAB's {s.name})"""
        if name not in self.fields:
            raise KeyError(name)
        return [e[name] for e in self._elements]

    def __repr__(self) -> str:
        size = "x".join(str(n) for n in self._shape)
        return f"StructArray({size}, fields={list(self.fields)})"


# ------------------------------------------------------------------------------------------------
# Python values to MATLAB values
# ------------------------------------------------------------------------------------------------


@dataclass
class _Node:
    """a MATLAB value: cls the class; shape MATLAB's size; data an F-ordered array of that size
    (char: uint16 code units; logical: uint8), a column-major list of _Node (cell), a
    column-major list of dicts of _Node (struct) or a CSC matrix (sparse)"""

    cls: str
    shape: tuple[int, ...]
    data: Any = None
    fields: tuple[str, ...] = ()
    sparse: bool = False

    @property
    def numel(self) -> int:
        return math.prod(self.shape)


def _mshape(shape: Sequence[int]) -> tuple[int, ...]:
    """MATLAB's size: at least two dimensions, trailing singletons after the second dropped"""
    s = [int(n) for n in shape]
    if len(s) == 0:
        s = [1, 1]
    elif len(s) == 1:
        s = [1, s[0]]
    while len(s) > 2 and s[-1] == 1:
        s.pop()
    return tuple(s)


def _check_name(name: Any, what: str) -> str:
    if not isinstance(name, str) or not _NAME.fullmatch(name):
        raise ValueError(
            f"{what} {name!r} is not a MATLAB name (a letter, then letters, digits "
            "or underscores, at most 63)"
        )
    return name


def _char_node(s: str) -> _Node:
    units = np.frombuffer(s.encode("utf-16-le", "surrogatepass"), dtype="<u2")
    if units.size == 0:
        return _Node("char", (0, 0), np.zeros((0, 0), np.uint16))
    return _Node("char", (1, units.size), units.reshape(1, -1))


def _array_node(a: np.ndarray, where: str) -> _Node:
    shape = _mshape(a.shape)
    kind = a.dtype.kind
    if kind == "U":
        if a.ndim == 0:
            return _char_node(str(a[()]))
        if a.dtype.itemsize != 4:
            raise ValueError(
                f"{where}: string arrays with more than one character per element "
                "become char only as <U1 arrays; use a list for a cell of char"
            )
        flat = a.reshape(-1, order="F")
        units = np.empty(flat.shape, np.uint16)
        for i, ch in enumerate(flat):
            if len(ch) != 1 or ord(ch) > 0xFFFF:
                raise ValueError(
                    f"{where}: a char array element must be one UTF-16 code unit, not {ch!r}"
                )
            units[i] = ord(ch)
        return _Node("char", shape, units.reshape(shape, order="F"))
    if kind == "O":
        flat = a.reshape(-1, order="F")
        return _Node(
            "cell", shape, [_to_node(v, f"{where}{{{i + 1}}}") for i, v in enumerate(flat)]
        )
    cls = _DTYPE_CLASS.get(a.dtype.newbyteorder("=")) if kind in "biufc" else None
    if cls is None:
        raise TypeError(f"{where}: numpy dtype {a.dtype} has no MATLAB class")
    a = np.asarray(a, dtype=a.dtype.newbyteorder("="))
    if kind == "b":
        a = a.view(np.uint8)
    return _Node(cls, shape, np.reshape(a, shape, order="F"))


def _sparse_node(m: Any, where: str) -> _Node:
    import scipy.sparse

    m = scipy.sparse.csc_matrix(m)
    if m.dtype == np.bool_:
        cls = "logical"
    elif m.dtype in (np.float64, np.complex128):
        cls = "double"
    else:
        raise TypeError(f"{where}: MATLAB sparse matrices are double or logical, not {m.dtype}")
    m.sum_duplicates()
    m.sort_indices()
    return _Node(cls, (int(m.shape[0]), int(m.shape[1])), m, sparse=True)


def _to_node(v: Any, where: str) -> _Node:
    if isinstance(v, _Node):
        return v
    if isinstance(v, MatlabObject):
        raise TypeError(f"{where}: a MATLAB {v.class_name} object can't be written")
    if isinstance(v, (bool, np.bool_)):
        return _Node("logical", (1, 1), np.full((1, 1), bool(v), np.uint8))
    if isinstance(v, int) and not isinstance(v, np.integer):
        try:
            f = float(v)
        except OverflowError:
            f = math.inf
        if math.isinf(f) or int(f) != v:
            raise ValueError(
                f"{where}: the int {v} is not exactly a double; pass a numpy "
                "integer for an integer class"
            )
        return _Node("double", (1, 1), np.full((1, 1), f))
    if isinstance(v, float) and not isinstance(v, np.floating):
        return _Node("double", (1, 1), np.full((1, 1), v))
    if isinstance(v, complex) and not isinstance(v, np.complexfloating):
        return _Node("double", (1, 1), np.full((1, 1), v, np.complex128))
    if isinstance(v, str):
        return _char_node(str(v))
    if isinstance(v, (bytes, bytearray)):
        b = np.frombuffer(bytes(v), np.uint8)
        return _Node(
            "uint8", (1, b.size) if b.size else (0, 0), b.reshape((1, -1) if b.size else (0, 0))
        )
    if v is None:
        return _Node("double", (0, 0), np.zeros((0, 0)))
    if isinstance(v, (list, tuple)):
        if not v:
            return _Node("cell", (0, 0), [])
        return _Node(
            "cell", (1, len(v)), [_to_node(x, f"{where}{{{i + 1}}}") for i, x in enumerate(v)]
        )
    if isinstance(v, Mapping):
        fields = tuple(_check_name(k, f"{where}: field") for k in v)
        return _Node(
            "struct", (1, 1), [{k: _to_node(v[k], f"{where}.{k}") for k in fields}], fields
        )
    if isinstance(v, StructArray):
        fields = tuple(_check_name(k, f"{where}: field") for k in v.fields)
        els = [
            {k: _to_node(e[k], f"{where}({i + 1}).{k}") for k in fields} for i, e in enumerate(v)
        ]
        return _Node("struct", _mshape(v.shape), els, fields)
    if isinstance(v, (np.ndarray, np.generic)):
        return _array_node(np.asarray(v), where)
    if type(v).__module__.startswith("scipy.sparse") and hasattr(v, "tocsc"):
        return _sparse_node(v, where)
    raise TypeError(f"{where}: {type(v).__name__} has no MATLAB equivalent")


# ------------------------------------------------------------------------------------------------
# MATLAB values to Python values (shared by both readers)
# ------------------------------------------------------------------------------------------------


def _fortran(flat: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    return np.reshape(flat, shape, order="F")


def _py_numeric(a: np.ndarray) -> Any:
    if a.shape == (1, 1):
        return a.reshape(-1)[0]
    return a


def _py_char(units: np.ndarray, shape: tuple[int, ...]) -> Any:
    """units: uint16 code units, column-major"""
    units = np.asarray(units, dtype="<u2").reshape(-1, order="F")
    if shape == (0, 0):
        return ""
    if len(shape) == 2 and shape[0] == 1 and shape[1] >= 1:
        return units.tobytes().decode("utf-16-le", "surrogatepass")
    chars = np.array([chr(u) for u in units.tolist()], dtype="<U1")
    return _fortran(chars, shape)


def _py_cell(items: list, shape: tuple[int, ...]) -> Any:
    if shape == (0, 0):
        return []
    if len(shape) == 2 and shape[0] == 1 and shape[1] >= 1:
        return items
    a = np.empty(math.prod(shape), dtype=object)
    for i, x in enumerate(items):
        a[i] = x
    return _fortran(a, shape)


def _py_struct(elements: list[dict], shape: tuple[int, ...], fields: tuple[str, ...]) -> Any:
    if shape == (1, 1):
        return elements[0]
    return StructArray(shape, fields, elements)


def _py_strings(data: np.ndarray) -> Any:
    """MATLAB string array data (the 'any' property of a string object): uint64 words
    [1, ndims, dims..., one character count per element, UTF-16 units packed 4 per word]"""
    w = np.asarray(data, dtype="<u8").reshape(-1, order="F")
    nd = int(w[1])
    shape = tuple(int(x) for x in w[2 : 2 + nd])
    n = math.prod(shape)
    counts = w[2 + nd : 2 + nd + n]
    units = w[2 + nd + n :].view("<u2")
    out, pos = [], 0
    for c in counts.tolist():
        if c == 0xFFFFFFFFFFFFFFFF:
            out.append(None)
            continue
        out.append(units[pos : pos + c].tobytes().decode("utf-16-le", "surrogatepass"))
        pos += c
    if shape == (1, 1) and out[0] is not None:
        return out[0]
    a = np.empty(n, dtype=object)
    a[:] = out
    return _fortran(a, shape)


class _Mcos:
    """MATLAB's object store (the FileWrapper__ cell of the subsystem): enough of it to read
    string objects. cell(i) returns the raw value of cell element i."""

    def __init__(self, blob: bytes, cell: Any) -> None:
        self.cell = cell
        u32 = np.frombuffer(blob, "<u4")
        nnames = int(u32[1])
        off = [int(x) for x in u32[2:10]]
        names = blob[40 : off[0]].split(b"\x00")[:nnames]
        self.names = [n.decode("ascii", "replace") for n in names]
        # one row per object, the object id being the row number: class index, two unknown
        # words, property list index in segment 2, in segment 4, and a dependency id
        objs = np.frombuffer(blob[off[2] : off[3]], "<u4").reshape(-1, 6)
        self.objects = {i: (int(r[0]), int(r[3]), int(r[4])) for i, r in enumerate(objs) if i}
        self.props = [self._props(blob[off[1] : off[2]]), self._props(blob[off[3] : off[4]])]

    def _name(self, i: int) -> str:
        return self.names[i - 1] if 0 < i <= len(self.names) else ""

    @staticmethod
    def _props(seg: bytes) -> list[list[tuple[int, int, int]]]:
        """per-object property lists of one segment: (name index, kind, value)"""
        u = np.frombuffer(seg, "<u4")
        out, pos = [], 0
        while pos < u.size:
            n = int(u[pos])
            out.append(
                [tuple(int(x) for x in u[pos + 1 + 3 * k : pos + 4 + 3 * k]) for k in range(n)]
            )
            pos += 1 + 3 * n
            pos += pos % 2  # entries are padded to 8 bytes
        return out

    def value(self, meta: np.ndarray, class_name: str) -> Any:
        """the value of an object variable: meta is its uint32 data
        [0xDD000000, ndims, dims..., object ids..., class id]"""
        m = np.asarray(meta, dtype=np.uint32).reshape(-1, order="F")
        if class_name != "string" or m.size < 5 or m[0] != 0xDD000000:
            return MatlabObject(class_name)
        nd = int(m[1])
        ids = m[2 + nd : 2 + nd + math.prod(int(x) for x in m[2 : 2 + nd])]
        if ids.size != 1 or int(ids[0]) not in self.objects:
            return MatlabObject(class_name, "string object not found")
        _, p2, p4 = self.objects[int(ids[0])]
        props = self.props[0][p2] if p2 else self.props[1][p4] if p4 < len(self.props[1]) else []
        for name, kind, val in props:
            if self._name(name) == "any" and kind == 1:
                return _py_strings(self.cell(val + 2))
        return MatlabObject(class_name, "string data not found")


# ------------------------------------------------------------------------------------------------
# v7.3: HDF5
# ------------------------------------------------------------------------------------------------

_V73_TEXT = "MATLAB 7.3 MAT-file, Platform: {}, Created on: {} HDF5 schema 1.00 ."
# MATLAB stores small variables in the object header (compact), larger ones chunked and
# deflated at level 3 (level 0 with -nocompression), in chunks of at most 64 KiB
_COMPACT_MAX = 9000
_CHUNK_BYTES = 65536
_DEFLATE = 3


def _platform() -> str:
    if sys.platform.startswith("win"):
        return "PCWIN64"
    if sys.platform == "darwin":
        import platform

        return "MACA64" if platform.machine() == "arm64" else "MACI64"
    return "GLNXA64"


def _ctime() -> str:
    """time.ctime-like text in English whatever the locale ('Tue Nov  3 19:24:23 2026')"""
    t = time.localtime()
    day = "Mon Tue Wed Thu Fri Sat Sun".split()[t.tm_wday]
    mon = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()[t.tm_mon - 1]
    return f"{day} {mon} {t.tm_mday:2d} {t.tm_hour:02d}:{t.tm_min:02d}:{t.tm_sec:02d} {t.tm_year}"


def _chunks(shape: tuple[int, ...], itemsize: int) -> tuple[int, ...]:
    """MATLAB's chunk shape for an HDF5 shape: the number of 64 KiB pieces the data needs,
    taken from the largest dimension (the last of equal ones); when that dimension is too short
    it becomes 1 and the rest is taken from the next"""
    ch = list(shape)
    n = math.ceil(math.prod(ch) * itemsize / _CHUNK_BYTES)
    while n > 1:
        d = max(range(len(ch)), key=lambda i: (ch[i], i))
        if ch[d] <= 1:
            break
        if ch[d] >= n:
            ch[d] //= n
            break
        n = math.ceil(n / ch[d])
        ch[d] = 1
    return tuple(ch)


class _H5Writer:
    def __init__(self, f: Any, compression: bool) -> None:
        import h5py

        self.h5py = h5py
        self.f = f
        self.compression = compression
        self.refs = None
        self.count = 0

    # attributes as MATLAB writes them: fixed-length NULLTERM ASCII class names, int32 decode
    # hints, uint8 empty flags, a variable-length list of one-character strings for fields
    def _attr_str(self, obj: Any, name: str, s: str) -> None:
        h5t = self.h5py.h5t
        tid = h5t.C_S1.copy()
        tid.set_size(len(s))
        tid.set_strpad(h5t.STR_NULLTERM)
        aid = self.h5py.h5a.create(
            obj.id, name.encode(), tid, self.h5py.h5s.create(self.h5py.h5s.SCALAR)
        )
        aid.write(np.array(s.encode("ascii"), dtype=f"S{len(s)}"), mtype=tid)

    def _attrs(
        self,
        obj: Any,
        cls: str,
        *,
        decode: int = 0,
        empty: bool = False,
        fields: tuple[str, ...] | None = None,
        sparse: int | None = None,
    ) -> None:
        self._attr_str(obj, "MATLAB_class", cls)
        if decode:
            obj.attrs.create("MATLAB_int_decode", np.int32(decode))
        if empty:
            obj.attrs.create("MATLAB_empty", np.uint8(1))
        if sparse is not None:
            obj.attrs.create("MATLAB_sparse", np.uint64(sparse))
        if fields:
            # MATLAB's own files have NULLTERM one-character strings here; HDF5 would truncate
            # those to empty when converting, and MATLAB reads the NULLPAD ones h5py writes
            vt = self.h5py.vlen_dtype(np.dtype("S1"))
            data = np.empty(len(fields), dtype=vt)
            for i, k in enumerate(fields):
                data[i] = np.frombuffer(k.encode("ascii"), dtype="S1")
            obj.attrs.create("MATLAB_fields", data, dtype=vt)

    def _dataset(self, parent: Any, name: str, data: np.ndarray) -> Any:
        """data in HDF5 order (MATLAB's dimensions reversed), C-contiguous"""
        h5p, h5d = self.h5py.h5p, self.h5py.h5d
        kw: dict[str, Any] = {"track_times": False}
        if data.nbytes <= _COMPACT_MAX:
            dcpl = h5p.create(h5p.DATASET_CREATE)
            dcpl.set_layout(h5d.COMPACT)
            kw["dcpl"] = dcpl
        else:
            kw.update(
                chunks=_chunks(data.shape, data.dtype.itemsize),
                compression="gzip",
                compression_opts=_DEFLATE if self.compression else 0,
            )
        return parent.create_dataset(name, data=data, **kw)

    def _refs_group(self) -> Any:
        if self.refs is None:
            self.refs = self.f.create_group("#refs#")
            ds = self._dataset(self.refs, "a", np.zeros(2, np.uint64))
            self._attrs(ds, "canonical empty", empty=True)
        return self.refs

    def _ref(self, node: _Node) -> Any:
        """write a value into #refs# and return a reference to it"""
        g = self._refs_group()
        self.count += 1
        i, name = self.count, ""
        while True:  # b, c, ..., z, ba, bb, ...: 'a' is the canonical empty
            name = "abcdefghijklmnopqrstuvwxyz"[i % 26] + name
            i //= 26
            if i == 0:
                break
        return self.write(g, name, node).ref

    def write(self, parent: Any, name: str, n: _Node) -> Any:
        h5py = self.h5py
        hshape = tuple(reversed(n.shape))
        if n.sparse:
            m = n.data
            g = parent.create_group(name, track_order=False)
            self._attrs(g, n.cls, decode=1 if n.cls == "logical" else 0, sparse=n.shape[0])
            if m.nnz:
                d = m.data.astype(np.uint8) if n.cls == "logical" else m.data
                self._dataset(g, "data", self._complex(d) if np.iscomplexobj(d) else d)
                self._dataset(g, "ir", m.indices.astype(np.uint64))
            self._dataset(g, "jc", m.indptr.astype(np.uint64))
            return g
        if n.cls == "struct" and n.fields and n.numel:
            g = parent.create_group(name)
            self._attrs(g, "struct", fields=n.fields)
            if n.shape == (1, 1):
                for k in n.fields:
                    self.write(g, k, n.data[0][k])
            else:
                for k in n.fields:
                    refs = np.array([self._ref(e[k]) for e in n.data], dtype=h5py.ref_dtype)
                    self._dataset(g, k, refs.reshape(hshape))
            return g
        if n.numel == 0 or n.cls == "struct":
            ds = self._dataset(parent, name, np.array(n.shape, dtype=np.uint64))
            self._attrs(ds, n.cls, empty=True, fields=n.fields if n.cls == "struct" else None)
            return ds
        if n.cls == "cell":
            refs = np.array([self._ref(x) for x in n.data], dtype=h5py.ref_dtype)
            ds = self._dataset(parent, name, refs.reshape(hshape))
            self._attrs(ds, "cell")
            return ds
        a = n.data.T  # C-ordered, HDF5's shape
        if np.iscomplexobj(a):
            a = self._complex(a)
        ds = self._dataset(parent, name, np.ascontiguousarray(a))
        self._attrs(ds, n.cls, decode={"logical": 1, "char": 2}.get(n.cls, 0))
        return ds

    @staticmethod
    def _complex(a: np.ndarray) -> np.ndarray:
        t = a.real.dtype
        c = np.empty(a.shape, dtype=[("real", t), ("imag", t)])
        c["real"] = a.real
        c["imag"] = a.imag
        return c


def _write_v73(path: PathLike, nodes: dict[str, _Node], compression: bool) -> None:
    import h5py

    # no HDF5 file locks (they fail on some network file systems); results are published by
    # rename, so nobody reads a file while it is written
    with h5py.File(path, "w", userblock_size=512, libver="earliest", locking=False) as f:
        w = _H5Writer(f, compression)
        for name, n in nodes.items():
            w.write(f, name, n)
    text = _V73_TEXT.format(_platform(), _ctime()).encode("ascii").ljust(116)
    with open(path, "r+b") as fh:
        fh.write(text + bytes(8) + b"\x00\x02IM")


class _H5Reader:
    def __init__(self, f: Any) -> None:
        import h5py

        self.h5py = h5py
        self.f = f
        self._mcos: _Mcos | None = None

    def mcos(self) -> _Mcos | None:
        if self._mcos is None and "#subsystem#/MCOS" in self.f:
            refs = self.f["#subsystem#/MCOS"][()].reshape(-1)
            blob = self.f[refs[0]][()].tobytes()
            self._mcos = _Mcos(blob, lambda i: self.f[refs[i]][()].T)
        return self._mcos

    def read(self, obj: Any) -> Any:
        h5py = self.h5py
        a = obj.attrs
        cls = a.get("MATLAB_class", b"")
        cls = cls.decode("ascii") if isinstance(cls, bytes) else str(cls)
        if "MATLAB_object_decode" in a:
            if cls == "string" and isinstance(obj, h5py.Dataset):
                m = self.mcos()
                return m.value(obj[()].T, cls) if m else MatlabObject(cls, "no object store")
            return MatlabObject(cls)
        if isinstance(obj, h5py.Group):
            if "MATLAB_sparse" in a:
                return self._sparse(obj, cls, int(a["MATLAB_sparse"]))
            if cls != "struct":
                return MatlabObject(cls or "group")
            fields = self._fields(obj)
            first = obj[fields[0]] if fields else None
            if (
                isinstance(first, h5py.Dataset)
                and first.dtype == h5py.ref_dtype
                and "MATLAB_class" not in first.attrs
            ):
                shape = tuple(reversed(first.shape))
                cols = {k: [self.read(self.f[r]) for r in obj[k][()].reshape(-1)] for k in fields}
                els = [{k: cols[k][i] for k in fields} for i in range(math.prod(shape))]
                return _py_struct(els, _mshape(shape), fields)
            return {k: self.read(obj[k]) for k in fields}
        if a.get("MATLAB_empty", 0):
            shape = _mshape(obj[()].reshape(-1).tolist()) if obj.size else (0, 0)
            if cls == "struct":
                fields = self._fields(obj) if "MATLAB_fields" in a else ()
                els = [dict.fromkeys(fields) for _ in range(math.prod(shape))]
                return _py_struct(els, shape, fields)
            return self._empty(cls, shape)
        data = obj[()]
        if cls == "cell":
            shape = _mshape(tuple(reversed(data.shape)))
            return _py_cell([self.read(self.f[r]) for r in data.reshape(-1)], shape)
        if cls not in _NUMERIC and cls not in ("char", "logical"):
            # a dataset without a MATLAB class (an HDF5 file not written by MATLAB) is read as is
            return _py_numeric(np.asarray(data).T) if cls == "" else MatlabObject(cls)
        a2 = np.asarray(data).T
        shape = _mshape(a2.shape)
        if a2.dtype.names:
            a2 = a2["real"] + 1j * a2["imag"].astype(a2["real"].dtype)
            if cls not in ("double", "single"):
                a2 = a2.astype(np.complex128)
        if cls == "char":
            return _py_char(a2, shape)
        if cls == "logical":
            return _py_numeric(np.reshape(a2.astype(np.bool_), shape, order="F"))
        if not a2.dtype.names and a2.dtype.kind != "c":
            a2 = a2.astype(_NUMERIC[cls], copy=False)
        return _py_numeric(np.reshape(a2, shape, order="F"))

    @staticmethod
    def _fields(obj: Any) -> tuple[str, ...]:
        if "MATLAB_fields" in obj.attrs:
            return tuple(
                b"".join(np.asarray(x).tolist()).decode("ascii") for x in obj.attrs["MATLAB_fields"]
            )
        return tuple(obj.keys())

    @staticmethod
    def _empty(cls: str, shape: tuple[int, ...]) -> Any:
        if cls == "char":
            return _py_char(np.zeros(0, np.uint16), shape)
        if cls == "cell":
            return _py_cell([], shape)
        if cls == "logical":
            return np.zeros(shape, np.bool_, order="F")
        if cls == "canonical empty":
            return np.zeros((0, 0))
        if cls in _NUMERIC:
            return np.zeros(shape, _NUMERIC[cls], order="F")
        return MatlabObject(cls)

    def _sparse(self, g: Any, cls: str, nrows: int) -> Any:
        import scipy.sparse

        jc = g["jc"][()].astype(np.int64)
        ir = g["ir"][()].astype(np.int64) if "ir" in g else np.zeros(0, np.int64)
        if "data" in g:
            d = g["data"][()]
            if d.dtype.names:
                d = d["real"] + 1j * d["imag"]
        else:
            d = np.zeros(0)
        d = (
            d.astype(np.bool_)
            if cls == "logical"
            else d.astype(np.complex128 if np.iscomplexobj(d) else np.float64)
        )
        return scipy.sparse.csc_matrix((d, ir, jc), shape=(nrows, jc.size - 1))


def _read_v73(path: PathLike, names: Sequence[str] | None) -> dict[str, Any]:
    import h5py

    out = {}
    with h5py.File(path, "r", locking=False) as f:
        r = _H5Reader(f)
        keys = [k for k in f if not k.startswith("#")]
        for k in keys if names is None else names:
            if k in keys:
                out[k] = r.read(f[k])
            else:
                warnings.warn(f"Variable '{k}' not found in {os.fspath(path)}", stacklevel=3)
    return out


# ------------------------------------------------------------------------------------------------
# v6 and v7: MAT 5
# ------------------------------------------------------------------------------------------------

_MI = {
    1: "i1",
    2: "u1",
    3: "i2",
    4: "u2",
    5: "i4",
    6: "u4",
    7: "f4",
    9: "f8",
    12: "i8",
    13: "u8",
    16: "utf8",
    17: "utf16",
    18: "utf32",
}
_MI_OF = {
    "double": 9,
    "single": 7,
    "int8": 1,
    "uint8": 2,
    "int16": 3,
    "uint16": 4,
    "int32": 5,
    "uint32": 6,
    "int64": 12,
    "uint64": 13,
    "logical": 2,
    "char": 4,
}
_MX_OF = {
    "cell": 1,
    "struct": 2,
    "char": 4,
    "double": 6,
    "single": 7,
    "int8": 8,
    "uint8": 9,
    "int16": 10,
    "uint16": 11,
    "int32": 12,
    "uint32": 13,
    "int64": 14,
    "uint64": 15,
    "logical": 9,
}
_MX_CLASS = {v: k for k, v in _MX_OF.items() if k != "logical"}
_V5_MAX = 2**31 - 1


class _V5Writer:
    """MAT 5 data elements, as a list of byte strings and array buffers (no copies of data)"""

    def __init__(self, utf8: bool) -> None:
        self.parts: list[Any] = []
        self.utf8 = utf8

    @staticmethod
    def _tag(mi: int, n: int) -> bytes:
        return struct.pack("<II", mi, n)

    def _elem(self, mi: int, data: Any, n: int) -> int:
        """append one element; data is bytes or a contiguous array of n bytes"""
        if 0 < n <= 4:
            raw = bytes(data) if not isinstance(data, np.ndarray) else data.tobytes()
            self.parts.append(struct.pack("<HH", mi, n) + raw.ljust(4, b"\x00"))
            return 8
        self.parts.append(self._tag(mi, n))
        if n:
            self.parts.append(data)
        pad = -n % 8
        if pad:
            self.parts.append(bytes(pad))
        return 8 + n + pad

    def matrix(self, n: _Node, name: str) -> int:
        """append a miMATRIX element; returns its size"""
        self.parts.append(None)  # the tag, filled in below
        at = len(self.parts) - 1
        size = 0
        flags = 0
        nzmax = 0
        if n.sparse:
            cls_code, flags, nzmax = 5, 0x10, max(n.data.nnz, 1)
            if n.cls == "logical":
                flags |= 0x02
            if np.iscomplexobj(n.data.data):
                flags |= 0x08
        else:
            cls_code = _MX_OF[n.cls]
            if n.cls == "logical":
                flags = 0x02
            elif n.cls not in ("cell", "struct", "char") and np.iscomplexobj(n.data):
                flags = 0x08
        size += self._elem(6, struct.pack("<II", cls_code | flags << 8, nzmax), 8)
        size += self._elem(5, struct.pack(f"<{len(n.shape)}i", *n.shape), 4 * len(n.shape))
        nb = name.encode("ascii")
        size += self._elem(1, nb, len(nb))
        if n.sparse:
            m = n.data
            ir = m.indices.astype("<i4")
            jc = m.indptr.astype("<i4")
            size += self._elem(5, ir, ir.nbytes) + self._elem(5, jc, jc.nbytes)
            if n.cls == "logical":
                d = m.data.astype(np.uint8)
                size += self._elem(2, d, d.nbytes)
            else:
                re_ = np.ascontiguousarray(m.data.real, "<f8")
                size += self._elem(9, re_, re_.nbytes)
                if flags & 0x08:
                    im = np.ascontiguousarray(m.data.imag, "<f8")
                    size += self._elem(9, im, im.nbytes)
        elif n.cls == "cell":
            for x in n.data:
                size += self.matrix(x, "")
        elif n.cls == "struct":
            ln = max([len(k) for k in n.fields] + [0]) + 1
            size += self._elem(5, struct.pack("<i", ln), 4)
            names = b"".join(k.encode("ascii").ljust(ln, b"\x00") for k in n.fields)
            size += self._elem(1, names, len(names))
            for e in n.data:
                for k in n.fields:
                    size += self.matrix(e[k], "")
        else:
            a = n.data
            mi = _MI_OF[n.cls]
            if n.cls == "char" and n.numel:
                # as MATLAB: UTF-8 bytes for ASCII text (v7), UTF-16 otherwise; MATLAB reads
                # uint16 char data through the session's default encoding, losing non-ASCII
                ascii = bool((a < 128).all())
                if ascii and self.utf8:
                    b = np.ascontiguousarray(a.reshape(-1, order="F").astype(np.uint8))
                    size += self._elem(16, b, b.nbytes)
                    self.parts[at] = self._tag(14, size)
                    return size + 8
                if not ascii:
                    mi = 17
            parts = (a.real, a.imag) if np.iscomplexobj(a) else (a,)
            for p in parts:
                p = np.asarray(p, dtype=p.dtype.newbyteorder("<")).reshape(-1, order="F")
                p = np.ascontiguousarray(p)
                size += self._elem(mi, p, p.nbytes)
        self.parts[at] = self._tag(14, size)
        return size + 8


def _write_v5(path: PathLike, nodes: dict[str, _Node], compress: bool) -> None:
    text = f"MATLAB 5.0 MAT-file, Platform: {_platform()}, Created on: {_ctime()}"
    with open(path, "wb") as fh:
        fh.write(text.encode("ascii").ljust(116) + bytes(8) + b"\x00\x01IM")
        for name, n in nodes.items():
            w = _V5Writer(utf8=compress)
            size = w.matrix(n, name)
            if size - 8 > _V5_MAX:
                raise ValueError(
                    f"variable {name!r} needs {size} bytes; -v6 and -v7 hold at "
                    "most 2 GB per variable (use version '7.3')"
                )
            if compress:
                z = zlib.compressobj(_DEFLATE)
                chunks = [
                    z.compress(memoryview(p).cast("B") if isinstance(p, np.ndarray) else p)
                    for p in w.parts
                ]
                chunks.append(z.flush())
                body = b"".join(chunks)
                fh.write(struct.pack("<II", 15, len(body)))
                fh.write(body)
            else:
                for p in w.parts:
                    fh.write(memoryview(p).cast("B") if isinstance(p, np.ndarray) else p)


class _V5Reader:
    def __init__(self, buf: bytes | memoryview, order: str, mcos: _Mcos | None = None) -> None:
        self.buf = buf
        self.o = order
        self.mcos = mcos

    def elem(self, pos: int) -> tuple[int, memoryview, int]:
        """(type, data, position of the next element)"""
        t, n = struct.unpack_from(self.o + "II", self.buf, pos)
        if t >> 16:  # small data element: the byte count in the upper half of the type word
            n, t = t >> 16, t & 0xFFFF
            return t, memoryview(self.buf)[pos + 4 : pos + 4 + n], pos + 8
        end = pos + 8 + n
        return t, memoryview(self.buf)[pos + 8 : end], end + (-n % 8 if t != 15 else 0)

    def numbers(self, t: int, data: memoryview) -> np.ndarray:
        code = _MI.get(t)
        if code is None or code.startswith("utf"):
            raise ValueError(f"unexpected MAT-file data type {t}")
        return np.frombuffer(data, dtype=self.o + code)

    def text(self, data: memoryview) -> str:
        return bytes(data).split(b"\x00", 1)[0].decode("ascii", "replace")

    def matrix(self, data: memoryview, want: Iterable[str] | None = None) -> tuple[str, Any]:
        """the name and value of the payload of a miMATRIX element"""
        sub = _V5Reader(data, self.o, self.mcos)
        if len(data) == 0:  # an empty element: []
            return "", np.zeros((0, 0))
        _, fl, pos = sub.elem(0)
        w0, nzmax = struct.unpack(self.o + "II", bytes(fl[:8]))
        mx, flags = w0 & 0xFF, (w0 >> 8) & 0xFF
        if mx == 17:  # opaque: name, type, class, data
            _, nm, pos = sub.elem(pos)
            _, _, pos = sub.elem(pos)
            _, cn, pos = sub.elem(pos)
            name, cname = sub.text(nm), sub.text(cn)
            if want is not None and name not in want:
                return name, None
            _, d, pos = sub.elem(pos)
            _, meta = sub.matrix(d)
            if cname == "FileWrapper__":
                return name, meta
            if self.mcos is None:
                return name, MatlabObject(cname, "no object store")
            return name, self.mcos.value(meta, cname)
        _, dm, pos = sub.elem(pos)
        shape = _mshape(sub.numbers(5, dm).tolist())
        _, nm, pos = sub.elem(pos)
        name = sub.text(nm)
        if want is not None and name not in want:
            return name, None
        numel = math.prod(shape)
        if mx == 1:
            items = []
            for _ in range(numel):
                _, d, pos = sub.elem(pos)
                items.append(sub.matrix(d)[1])
            return name, _py_cell(items, shape)
        if mx in (2, 3):
            cname = ""
            if mx == 3:
                _, cn, pos = sub.elem(pos)
                cname = sub.text(cn)
            _, ln, pos = sub.elem(pos)
            ln = int(sub.numbers(5, ln)[0])
            _, fn, pos = sub.elem(pos)
            raw = bytes(fn)
            fields = (
                tuple(
                    raw[i : i + ln].split(b"\x00", 1)[0].decode("ascii")
                    for i in range(0, len(raw), ln)
                )
                if ln
                else ()
            )
            els = []
            for _ in range(numel):
                e = {}
                for k in fields:
                    _, d, pos = sub.elem(pos)
                    e[k] = sub.matrix(d)[1]
                els.append(e)
            if mx == 3:
                return name, MatlabObject(cname)
            return name, _py_struct(els, shape, fields)
        if mx == 16:
            return name, MatlabObject("function_handle")
        if mx == 5:
            import scipy.sparse

            t, d, pos = sub.elem(pos)
            ir = sub.numbers(t, d).astype(np.int64)
            t, d, pos = sub.elem(pos)
            jc = sub.numbers(t, d).astype(np.int64)
            nnz = int(jc[-1]) if jc.size else 0
            t, d, pos = sub.elem(pos)
            pr = sub.numbers(t, d)[:nnz]
            if flags & 0x02:
                vals = pr.astype(np.bool_)
            else:
                vals = pr.astype(np.float64)
                if flags & 0x08:
                    t, d, pos = sub.elem(pos)
                    vals = vals + 1j * sub.numbers(t, d)[:nnz].astype(np.float64)
            return name, scipy.sparse.csc_matrix((vals, ir[:nnz], jc), shape=shape)
        cls = "char" if mx == 4 else _MX_CLASS.get(mx)
        if cls is None:
            return name, MatlabObject(f"class {mx}")
        t, d, pos = sub.elem(pos)
        if cls == "char":
            if t == 16:
                units = np.frombuffer(
                    bytes(d).decode("utf-8", "surrogatepass").encode("utf-16-le", "surrogatepass"),
                    "<u2",
                )
            elif t == 17:
                units = np.frombuffer(d, self.o + "u2").astype(np.uint16)
            elif t == 18:
                units = np.frombuffer(
                    bytes(d)
                    .decode("utf-32-be" if self.o == ">" else "utf-32-le")
                    .encode("utf-16-le", "surrogatepass"),
                    "<u2",
                )
            else:
                units = sub.numbers(t, d).astype(np.uint16)
            return name, _py_char(units, shape)
        a = sub.numbers(t, d) if numel else np.zeros(0)
        if flags & 0x02:
            return name, _py_numeric(_fortran(a.astype(np.bool_), shape))
        dt = np.dtype(_NUMERIC[cls])
        if flags & 0x08:
            t, d, pos = sub.elem(pos)
            im = sub.numbers(t, d) if numel else np.zeros(0)
            ct = np.complex64 if cls == "single" else np.complex128
            v = np.empty(numel, ct)
            v.real = a
            v.imag = im
            return name, _py_numeric(_fortran(v, shape))
        return name, _py_numeric(_fortran(a.astype(dt), shape))

    def top(self, pos: int) -> tuple[int, memoryview | None, int]:
        """the next top-level variable: (type, its miMATRIX payload, next position)"""
        t, d, nxt = self.elem(pos)
        if t == 15:
            raw = zlib.decompressobj().decompress(d)
            t2, d2, _ = _V5Reader(raw, self.o).elem(0)
            return t2, d2, nxt
        return t, d, nxt


def _read_v5(path: PathLike, names: Sequence[str] | None) -> dict[str, Any]:
    with open(path, "rb") as fh:
        buf = fh.read()
    order = ">" if buf[126:128] == b"MI" else "<"
    sub = struct.unpack_from(order + "Q", buf, 116)[0]
    if buf[116:124] in (b" " * 8, b"\x00" * 8):
        sub = 0
    mcos = None
    if 128 <= sub < len(buf):
        mcos = _v5_mcos(buf, sub, order)
    r = _V5Reader(buf, order, mcos)
    out: dict[str, Any] = {}
    want = None if names is None else set(names)
    pos = 128
    while pos + 8 <= len(buf):
        if pos == sub:
            break
        t, d, nxt = r.top(pos)
        if t == 14 and d is not None:
            name, value = r.matrix(d, want)
            if want is None or name in want:
                out[name] = value
        pos = nxt
    if names is not None:
        for k in names:
            if k not in out:
                warnings.warn(f"Variable '{k}' not found in {os.fspath(path)}", stacklevel=3)
        out = {k: out[k] for k in names if k in out}
    return out


def _v5_mcos(buf: bytes, pos: int, order: str) -> _Mcos | None:
    """the object store in the subsystem data of a MAT 5 file: a uint8 array holding a small
    MAT stream whose first variable is a struct with an MCOS field (a FileWrapper__ cell)"""
    try:
        t, d, _ = _V5Reader(buf, order).top(pos)
        _, raw = _V5Reader(buf, order).matrix(d)
        raw = np.asarray(raw).tobytes()
        o = ">" if raw[2:4] == b"MI" else "<"
        r = _V5Reader(raw, o)
        t, d, _ = r.top(8)
        _, s = r.matrix(d)
        cell = s["MCOS"] if isinstance(s, dict) else None
        if not isinstance(cell, (list, np.ndarray)):
            return None
        items = (
            list(np.asarray(cell, dtype=object).reshape(-1, order="F"))
            if isinstance(cell, np.ndarray)
            else cell
        )
        blob = np.asarray(items[0]).tobytes()
        return _Mcos(blob, lambda i: items[i])
    except (ValueError, KeyError, IndexError, struct.error):
        return None


# ------------------------------------------------------------------------------------------------
# public functions
# ------------------------------------------------------------------------------------------------

_HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"


def mat_version(path: PathLike) -> str:
    """'7.3' for an HDF5 MAT-file, '5' for a v6 or v7 one (the two share a format); ValueError
    for anything else (v4 files are not supported)"""
    with open(path, "rb") as fh:
        head = fh.read(520)
    if head[512:520] == _HDF5_SIGNATURE or head[:8] == _HDF5_SIGNATURE:
        return "7.3"
    if (
        len(head) >= 128
        and head[126:128] in (b"IM", b"MI")
        and head[124:126] in (b"\x00\x01", b"\x01\x00")
    ):
        return "5"
    raise ValueError(f"{os.fspath(path)} is not a MAT-file (v6, v7 or v7.3)")


def read_mat(path: PathLike, variables: str | Sequence[str] | None = None) -> dict[str, Any]:
    """the variables of a MAT-file (all, or the ones named), as MATLAB's load returns them,
    with the types of the module's table. A named variable the file lacks is left out with a
    warning, as MATLAB's load does."""
    if isinstance(variables, str):
        variables = [variables]
    if not os.path.exists(path):
        raise FileNotFoundError(f"no such file: {os.fspath(path)}")
    v = mat_version(path)
    return _read_v73(path, variables) if v == "7.3" else _read_v5(path, variables)


def write_mat(
    path: PathLike, variables: Mapping[str, Any], *, version: str = "7.3", compression: bool = True
) -> None:
    """write variables to a MAT-file (replacing the file), as MATLAB's save does.

    version: '7.3' (HDF5; PetaKit5D's default), '7' (MAT 5, compressed) or '6' (MAT 5,
    uncompressed). compression=False is save's -nocompression (v7.3 and v7)."""
    if version not in ("7.3", "7", "6"):
        raise ValueError(f"version must be '7.3', '7' or '6', not {version!r}")
    nodes = {_check_name(k, "variable name"): _to_node(v, k) for k, v in variables.items()}
    if version == "7.3":
        _write_v73(path, nodes, compression)
    else:
        _write_v5(path, nodes, compression and version == "7")


def _ascii_number(x: float, width: int, prec: int) -> str:
    if math.isnan(x):
        return "NaN".rjust(width)
    if math.isinf(x):
        return ("Inf" if x > 0 else "-Inf").rjust(width)
    return f"{x:{width}.{prec}e}"


def save_ascii(path: PathLike, *values: Any, double: bool = False, tabs: bool = False) -> None:
    """MATLAB's save(path, names..., '-ascii'[, '-double'][, '-tabs']): each value's rows, values
    printed with 8 (or 16) significant digits. As in MATLAB, only real 2-D double and char
    values are written (char as its codes); others are skipped with a warning."""
    width, prec = (25, 16) if double else (16, 7)
    lines = []
    for i, v in enumerate(values):
        if v is not None and not isinstance(
            v, (str, bytes, bool, int, float, complex, np.ndarray, np.generic)
        ):
            raise TypeError(
                f"save_ascii: value {i + 1} is a {type(v).__name__}; give numbers, "
                "numpy arrays or str"
            )
        n = _to_node(v, f"value {i + 1}")
        if n.cls not in ("double", "char") or n.sparse or len(n.shape) > 2:
            warnings.warn(
                "Attempt to write an unsupported data type to an ASCII file. "
                f"Value {i + 1} not written to file.",
                stacklevel=2,
            )
            continue
        a = n.data
        if np.iscomplexobj(a):
            warnings.warn(
                f"Imaginary part of complex value {i + 1} not saved to ASCII file.", stacklevel=2
            )
            a = a.real
        a = a.astype(np.float64)
        for r in range(n.shape[0]):  # a 1x0 value is one empty line
            row = [_ascii_number(float(x), width, prec) + ("\t" if tabs else "") for x in a[r]]
            lines.append("".join(row) + "\n")
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("".join(lines))
