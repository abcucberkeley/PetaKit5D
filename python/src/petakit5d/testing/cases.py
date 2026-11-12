"""Oracle cases: the MATLAB call, the Python call, their inputs and how results are compared.

A case names a MATLAB function (or a MATLAB snippet) with its arguments, the Python callable
with its arguments, the inputs both sides read, and the comparison of each output and of the
files written. Inputs are generated once by deterministic Python functions; both sides read the
same bytes (MATLAB from the golden store, Python from the same generator, checked by hash).

    Case('crop.uint16_odd', matlab='crop3d_mex', args=[In('im'), [2, 3, 1, 9, 8, 5]],
         python='petakit5d.tools.crop:crop_3d', inputs={'im': Value(lambda: synth.texture(...))})

Placeholders in arguments: ``In(name)`` (the value of a ``Value`` input), ``InPath(name)`` (the
path of a ``File``, ``Tree`` or ``Data`` input in the case's work folder) and ``Work(sub)``
(the work folder, with a trailing slash). Both sides run in a fresh work folder; the files a
call writes there (new or changed) are its file output.

Roles: arguments are written as MATLAB takes them; a role says what one means, and the Python
call gets it in Python's form (``Vol``: a ``(y, x, z)`` array as ``(z, y, x)``; ``Box``: a
1-based inclusive box as ``(z0, y0, x0, z1, y1, x1)``; ``PerAxis``, ``Index1``, ``Range1``,
``Dim``, ``LinIndex``). ``out_roles`` map Python's results back to MATLAB's form before they
are compared and hashed. Roles change neither the case key nor MATLAB's arguments.
"""

from __future__ import annotations

import atexit
import copy
import dataclasses
import hashlib
import importlib
import json
import os
import shutil
import tempfile
import time
import traceback
import uuid
import warnings
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from petakit5d.testing.compare import (
    EXACT,
    FileTree,
    MArray,
    MCell,
    MPlaceholder,
    MRef,
    MString,
    MStruct,
    MValue,
    Spec,
    TreeSpec,
    mchar,
    msize,
    scan_tree,
    to_mvalue,
    vhash,
)

FORMAT = 2  # case key format; bump when the meaning of a key or of a stored result changes


# ---------------------------------------------------------------------------------------------
# placeholders and MATLAB literals
# ---------------------------------------------------------------------------------------------


def In(name: str) -> MPlaceholder:
    """the value of the ``Value`` input ``name``"""
    return MPlaceholder("input", name)


def InPath(name: str) -> MPlaceholder:
    """the path of the ``File``, ``Tree`` or ``Data`` input ``name`` (folders end with '/')"""
    return MPlaceholder("path", name)


def Work(sub: str = "") -> MPlaceholder:
    """the case's work folder with a trailing '/', followed by ``sub``"""
    return MPlaceholder("work", "", sub)


# ---------------------------------------------------------------------------------------------
# roles: what a MATLAB argument means, so that Python gets it in Python's form
# ---------------------------------------------------------------------------------------------

ROLE_KINDS = ("vol", "box", "peraxis", "index1", "range1", "dim", "linindex")
OUT_KINDS = ("vol", "box", "boxes", "yxz", "xyz", "index1", "range1", "dim", "linindex", "each")


class Role:
    """a MATLAB argument (as MATLAB takes it) with its meaning; MATLAB and the case key see the
    value, the Python call the value in Python's form"""

    __slots__ = ("kind", "value", "param")

    def __init__(self, kind: str, value: Any, param: Any = None):
        if kind not in ROLE_KINDS:
            raise ValueError(f"unknown role {kind!r}")
        self.kind, self.value, self.param = kind, value, param

    def __repr__(self) -> str:
        name = {
            "vol": "Vol",
            "box": "Box",
            "peraxis": "PerAxis",
            "index1": "Index1",
            "range1": "Range1",
            "dim": "Dim",
            "linindex": "LinIndex",
        }[self.kind]
        extra = "" if self.param is None else f", {self.param!r}"
        return f"{name}({self.value!r}{extra})"


def Vol(x: Any) -> Role:
    """an image array, MATLAB's ``(y, x, z, ...)``; Python gets ``(..., z, y, x)``
    (``arrays.from_matlab``)"""
    return Role("vol", x)


def Box(x: Any) -> Role:
    """MATLAB's ``[ymin xmin zmin ymax xmax zmax]`` (1-based, inclusive; 2D ``[ymin xmin ymax
    xmax]``); Python gets ``(z0, y0, x0, z1, y1, x1)`` (``bbox.box_from_matlab``)"""
    return Role("box", x)


def PerAxis(x: Any, order: str = "yxz") -> Role:
    """a per-axis vector in MATLAB's ``[y x z]`` (``'yxz'``) or ``[x y z]`` (``'xyz'``) order;
    Python gets ``(z, y, x)`` (``axes.from_yxz`` / ``from_xyz``)"""
    if order not in ("yxz", "xyz"):
        raise ValueError(f"PerAxis order is 'yxz' or 'xyz', got {order!r}")
    return Role("peraxis", x, order)


def Index1(x: Any) -> Role:
    """a 1-based index (or indices); Python gets them 0-based"""
    return Role("index1", x)


def Range1(x: Any) -> Role:
    """MATLAB's ``[first last]`` (1-based, inclusive); Python gets ``(start, stop)``"""
    return Role("range1", x)


def Dim(x: Any, ndim: int | None = None) -> Role:
    """MATLAB's dimension (1: y, 2: x, 3: z); Python gets the NumPy axis (``axes.dim_to_axis``;
    negative unless ``ndim`` is given)"""
    return Role("dim", x, ndim)


def LinIndex(x: Any, size: Any) -> Role:
    """MATLAB's 1-based linear index (or indices) into an array of MATLAB size ``size`` (a size,
    or ``In(name)`` / the name of a ``Value`` input); Python gets the ``(z, y, x)`` index
    tuple (``arrays.matlab_index_to_zyx``)"""
    return Role("linindex", x, size)


def _strip_roles(x: Any) -> Any:
    """the arguments as MATLAB takes them (roles removed; a container without roles is
    returned as it is, so that its reading as a MATLAB value cannot change)"""
    if isinstance(x, Role):
        return _strip_roles(x.value)
    if type(x) in (list, tuple):
        items = [_strip_roles(y) for y in x]
        if all(a is b for a, b in zip(items, x, strict=True)):
            return x
        return items if type(x) is list else tuple(items)
    if type(x) is dict:
        items = {k: _strip_roles(y) for k, y in x.items()}
        return x if all(items[k] is y for k, y in x.items()) else items
    return x


def _plain(v: Any) -> Any:
    """a value for the conversions: a MATLAB array as its data"""
    if isinstance(v, MArray):
        return v.data
    return v


def _as_int(v: Any) -> Any:
    a = np.asarray(_plain(v))
    if a.dtype.kind == "f":
        if not np.all(a == np.trunc(a)):
            raise ValueError(f"an index must be integral, got {v!r}")
        a = a.astype(np.int64)
    return int(a) if a.ndim == 0 else a


def _shift(v: Any, d: int) -> Any:
    """MATLAB's indices as integers moved by d, keeping the container (list, tuple or array)"""
    if isinstance(v, (list, tuple)):
        return type(v)(_shift(x, d) for x in v)
    return _as_int(v) + d


def _shift_keep(v: Any, d: int) -> Any:
    """Python's indices moved by d, keeping the container and an array's class (a float array
    of places stays float, as MATLAB's double indices)"""
    if isinstance(v, (list, tuple)):
        return type(v)(_shift_keep(x, d) for x in v)
    if isinstance(v, (np.ndarray, np.generic)):
        return v + np.asarray(d, dtype=v.dtype)
    return v + d


class _Sizes:
    """the MATLAB size a LinIndex role or output refers to"""

    def __init__(self, case: Case):
        self.case = case

    def __call__(self, param: Any) -> tuple[int, ...]:
        if isinstance(param, MPlaceholder):
            param = param.name
        if isinstance(param, str):
            x = self.case.inputs.get(param)
            if not isinstance(x, Value):
                raise ValueError(
                    f"case {self.case.id}: LinIndex size {param!r} names no Value input"
                )
            return tuple(x.mvalue().size)
        return tuple(int(n) for n in np.asarray(_plain(param)).ravel())


def role_to_python(role: Role, v: Any, sizes: Callable[[Any], tuple] | None = None) -> Any:
    """a MATLAB argument value in Python's form, as its role says"""
    from petakit5d.core import arrays, axes, bbox

    k = role.kind
    if k == "vol":
        return arrays.from_matlab(np.asarray(_plain(v)))
    if k == "box":
        return bbox.box_from_matlab(_plain(v))
    if k == "peraxis":
        f = axes.from_yxz if role.param == "yxz" else axes.from_xyz
        return f(_plain(v))
    if k == "index1":
        return _shift(v, -1)
    if k == "range1":
        return bbox.range_from_matlab(_plain(v))
    if k == "dim":
        return axes.dim_to_axis(int(_as_int(v)), role.param)
    if k == "linindex":
        size = sizes(role.param) if sizes is not None else tuple(role.param)
        return arrays.matlab_index_to_zyx(np.asarray(_as_int(v)) - 1, arrays.size_from_matlab(size))
    raise ValueError(k)


def out_to_matlab(role: Any, v: Any, sizes: Callable[[Any], tuple] | None = None) -> Any:
    """a Python result in MATLAB's form, as its out role says

    None: unchanged. 'vol': ``to_matlab``; 'box': ``box_to_matlab``; 'boxes': rows of boxes;
    'yxz' / 'xyz': a ``(z, y, x)`` vector as MATLAB's ``[y x z]`` / ``[x y z]``; 'index1': 0-based
    indices plus 1; 'range1': ``(start, stop)`` as ``[first last]``; ('dim', ndim): a NumPy
    axis as MATLAB's dimension; ('linindex', size): a ``(z, y, x)`` index tuple as MATLAB's
    1-based linear index into an array of MATLAB size ``size`` (or of the ``Value`` input of
    that name); ('each', role): the role on every item of a list (a cell); a list of roles:
    applied in order; a dict: by field; a callable: applied.
    """
    from petakit5d.core import arrays, axes, bbox

    if role is None:
        return v
    if callable(role) and not isinstance(role, str):
        return role(v)
    if isinstance(role, list):
        for r in role:
            v = out_to_matlab(r, v, sizes)
        return v
    if isinstance(role, dict):
        if not isinstance(v, Mapping):
            raise TypeError(f"out role {role!r} needs a dict result, got {type(v).__name__}")
        return {f: out_to_matlab(role.get(f), x, sizes) for f, x in v.items()}
    kind, param = role if isinstance(role, tuple) else (role, None)
    if kind not in OUT_KINDS:
        raise ValueError(f"unknown out role {kind!r}")
    if kind == "each":
        if not isinstance(v, (list, tuple)):
            raise TypeError(f"out role ('each', ...) needs a list result, got {type(v).__name__}")
        return [out_to_matlab(param, x, sizes) for x in v]
    if kind == "vol":
        return arrays.to_matlab(np.asarray(v))
    if kind == "box":
        return bbox.box_to_matlab(v)
    if kind == "boxes":
        return bbox.boxes_to_matlab(v)
    if kind == "yxz":
        return axes.to_yxz(v)
    if kind == "xyz":
        return axes.to_xyz(v)
    if kind == "index1":
        return _shift_keep(v, 1)
    if kind == "range1":
        return bbox.range_to_matlab(v)
    if kind == "dim":
        if param is None:
            raise ValueError("out role 'dim' needs ('dim', ndim)")
        return axes.axis_to_dim(int(v), param)
    if kind == "linindex":
        if param is None:
            raise ValueError("out role 'linindex' needs ('linindex', MATLAB size or input name)")
        size = sizes(param) if sizes is not None else tuple(param)
        k = arrays.zyx_to_matlab_index(tuple(v), arrays.size_from_matlab(size))
        return k + 1
    raise ValueError(kind)


def role_kind(role: Any) -> str | None:
    """the kind of an out role ('vol', 'box', ...; ('each', r) and lists by their last role)"""
    if role is None or callable(role) and not isinstance(role, str):
        return None
    if isinstance(role, list):
        return role_kind(role[-1]) if role else None
    if isinstance(role, tuple):
        return role_kind(role[1]) if role and role[0] == "each" else role[0]
    return role if isinstance(role, str) else None


def _numbers(v: MValue) -> np.ndarray | None:
    if isinstance(v, MArray) and v.data is not None and v.mclass not in ("char",):
        return np.asarray(v.data, dtype=np.float64).reshape(-1, order="F")
    return None


def _same(x: Any, y: np.ndarray) -> bool:
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    return x.shape == y.shape and bool(np.all((x == y) | (np.isnan(x) & np.isnan(y))))


def role_hints(case: Case, role: Any, ref: MValue, act: MValue) -> list[str]:
    """what a role's mistakes explain of an output that differs from MATLAB's (the array
    mapping's are in compare's hints): a box converted twice or not converted, a per-axis
    vector in another order, a first-index rule taken in C order"""
    from petakit5d.core import axes, bbox

    kind = role_kind(role)
    r, a = _numbers(ref), _numbers(act)
    if r is None or a is None:
        return []
    out = []
    if kind in ("box", None) and r.size in (4, 6) and a.size == r.size:
        try:
            twice = bbox.box_to_matlab(tuple(int(x) for x in r))
            once = bbox.box_from_matlab(r)
        except (ValueError, TypeError, OverflowError):
            twice = once = None
        if twice is not None and _same(twice, a):
            out.append(
                "Python's box is MATLAB's converted once more (box_to_matlab of a MATLAB box): "
                "a box converted twice (out role 'box' and box_to_matlab in py_out, or a Box "
                "role on a box already in Python's form)?"
            )
        elif once is not None and _same([np.inf if x is None else x for x in once], a):
            out.append(
                "Python's box is in Python's form ((z0, y0, x0, z1, y1, x1), 0-based, end "
                "exclusive): out role 'box' missing?"
            )
    if kind in ("yxz", "xyz", None) and r.size == 3 and a.size == 3 and not _same(r, a):
        if _same(r[[1, 0, 2]], a):
            out.append(
                "y and x swapped: Python's per-axis vector is MATLAB's in the other of [y x z] "
                "and [x y z] (PerAxis(..., 'xyz') for 'yxz', or out role 'xyz' for 'yxz')?"
            )
        else:
            frm = axes.from_yxz if kind != "xyz" else axes.from_xyz
            to = axes.to_yxz if kind != "xyz" else axes.to_xyz
            if _same(frm(r), a):
                out.append(
                    "Python's per-axis vector is in (z, y, x) order: not mapped back (out role "
                    f"'{kind or 'yxz'}' missing)?"
                )
            elif _same(to(r), a):
                out.append(
                    "Python's per-axis vector is MATLAB's mapped once more: converted twice?"
                )
    if kind == "linindex" and isinstance(role, tuple) and len(role) > 1:
        out += _first_index_hints(case, role[1], r, a)
    return out


def _first_index_hints(case: Case, param: Any, r: np.ndarray, a: np.ndarray) -> list[str]:
    """a linear index result that names another voxel of the same value: a first-index rule
    (first maximum, first nonzero) taken in C order"""
    from petakit5d.core import arrays

    name = param.name if isinstance(param, MPlaceholder) else param
    x = case.inputs.get(name) if isinstance(name, str) else None
    if not isinstance(x, Value) or r.shape != a.shape:
        return []
    vol = np.asarray(_plain(x.python_value()))
    flat = vol.reshape(-1, order="F")
    bad = np.flatnonzero(r != a)
    if not bad.size:
        return []
    kr, ka = int(r[bad[0]]) - 1, int(a[bad[0]]) - 1
    if not (0 <= kr < flat.size and 0 <= ka < flat.size) or flat[ka] != flat[kr]:
        return []
    zyx = arrays.from_matlab(vol)
    hits = np.flatnonzero(zyx.reshape(-1) == flat[kr])
    first_c = (
        arrays.zyx_to_matlab_index(np.unravel_index(int(hits[0]), zyx.shape), zyx.shape)
        if hits.size
        else None
    )
    if first_c == ka:
        return [
            f"Python's index {ka + 1} is the first voxel of value {flat[kr]!r} in C order of "
            "(z, y, x) (x fastest); MATLAB's first is in its linear order (y fastest, then x, "
            "then z): a first-index rule (first maximum, first nonzero) taken in C order? "
            "Search arrays.matlab_order(vol) instead"
        ]
    return [
        f"Python's index {ka + 1} names another voxel of the same value {flat[kr]!r} as "
        f"MATLAB's {kr + 1} (a tie): MATLAB takes the first in its linear order (y fastest, then "
        "x, then z)"
    ]


def mcell(*items: Any) -> MCell:
    """a 1xN MATLAB cell array (a Python list of numbers would be a numeric row)"""
    return MCell(
        (1, len(items)),
        [x if isinstance(x, (MValue, MPlaceholder)) else to_mvalue(x) for x in items],
    )


def mcol(values: Sequence[Any] | np.ndarray, mclass: str = "double") -> MArray:
    """a MATLAB column vector"""
    a = np.asarray(values, dtype=np.dtype("?") if mclass == "logical" else None)
    return mtyped(a.reshape(-1, 1), mclass)


def mtyped(value: Any, mclass: str) -> MArray:
    """a numeric or logical array of a given MATLAB class (MATLAB's uint16(x) for integral x)"""
    from petakit5d.testing.compare import _DTYPE

    a = np.asarray(value)
    if mclass not in _DTYPE or mclass == "char":
        raise ValueError(mclass)
    if (
        a.dtype.kind == "f"
        and mclass in ("int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64")
        and not np.all(a == np.trunc(a))
    ):
        raise ValueError("mtyped converts integral values only; round in the generator")
    a = a.astype(_DTYPE[mclass])
    if a.ndim == 0:
        a = a.reshape(1, 1)
    elif a.ndim == 1:
        a = a.reshape(1, -1)
    return MArray(mclass, a)


def mstring(*items: str | None) -> MString:
    """a MATLAB string array (1xN; a single item gives a string scalar)"""
    return MString((1, len(items)), list(items))


def mstruct(**fields: Any) -> MStruct:
    return MStruct(
        (1, 1),
        list(fields),
        [
            {
                k: x if isinstance(x, (MValue, MPlaceholder)) else to_mvalue(x)
                for k, x in fields.items()
            }
        ],
    )


def mempty(mclass: str = "double", size: Sequence[int] = (0, 0)) -> MArray:
    from petakit5d.testing.compare import _DTYPE

    return MArray(mclass, np.zeros(msize(size), _DTYPE[mclass]))


# ---------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------

_CACHE: Path | None = None


def _cache_dir() -> Path:
    global _CACHE
    if _CACHE is None:
        _CACHE = Path(tempfile.mkdtemp(prefix="pkpy_inputs_"))
        atexit.register(shutil.rmtree, _CACHE, True)
    return _CACHE


def tree_hash(root: str | os.PathLike, content: bool = True) -> str:
    """SHA-256 of a folder: relative paths with the bytes (or only the sizes) of every file"""
    root = Path(root)
    h = hashlib.sha256(b"pk1:tree:" if content else b"pk1:stat:")
    for dp, dns, fns in sorted(os.walk(root)):
        dns.sort()
        d = Path(dp)
        if not dns and not fns and d != root:
            h.update(d.relative_to(root).as_posix().encode() + b"/\n")
        for fn in sorted(fns):
            p = d / fn
            rel = p.relative_to(root).as_posix()
            if content:
                from petakit5d.testing.compare import sha256_file

                h.update(f"{rel}\0{sha256_file(p)}\n".encode())
            else:
                h.update(f"{rel}\0{p.stat().st_size}\n".encode())
    return h.hexdigest()


class Input:
    """base of the input kinds"""

    kind = ""

    def hash(self) -> str:
        raise NotImplementedError


class Value(Input):
    """an in-memory input: ``make()`` returns an array or value (deterministically)

    MATLAB receives ``to_mvalue(make())``; Python receives a fresh copy of what make returned.
    """

    kind = "value"

    def __init__(self, make: Callable[[], Any]):
        self.make = make
        self._v: Any = None
        self._mv: MValue | None = None
        self._done = False

    def _value(self) -> Any:
        if not self._done:
            self._v = self.make()
            self._done = True
        return self._v

    def python_value(self) -> Any:
        v = self._value()
        return v.copy(order="K") if isinstance(v, np.ndarray) else copy.deepcopy(v)

    def mvalue(self) -> MValue:
        if self._mv is None:
            self._mv = to_mvalue(self._value())
        return self._mv

    def hash(self) -> str:
        return vhash(self.mvalue())


class File(Input):
    """a generated file: ``make(path)`` writes it; it is copied to ``<work>/<name>``"""

    kind = "file"

    def __init__(self, make: Callable[[Path], Any], name: str):
        self.make, self.name = make, name
        self._local: Path | None = None
        self._h: str | None = None

    def local(self) -> Path:
        if self._local is None:
            d = _cache_dir() / uuid.uuid4().hex
            d.mkdir()
            p = d / Path(self.name).name
            self.make(p)
            if not p.is_file():
                raise RuntimeError(f"input generator for {self.name!r} did not write {p}")
            self._local = p
        return self._local

    def hash(self) -> str:
        if self._h is None:
            from petakit5d.testing.compare import sha256_file

            self._h = hashlib.sha256(b"pk1:file:" + sha256_file(self.local()).encode()).hexdigest()
        return self._h


class Tree(Input):
    """a generated folder: ``make(folder)`` fills it; it is copied to ``<work>/<name>/``"""

    kind = "tree"

    def __init__(self, make: Callable[[Path], Any], name: str):
        self.make, self.name = make, name
        self._local: Path | None = None
        self._h: str | None = None

    def local(self) -> Path:
        if self._local is None:
            d = _cache_dir() / uuid.uuid4().hex / self.name
            d.mkdir(parents=True)
            self.make(d)
            self._local = d
        return self._local

    def hash(self) -> str:
        if self._h is None:
            self._h = tree_hash(self.local())
        return self._h


class Data(Input):
    """existing data (a demo dataset or a real-data copy), named by a path that may use
    environment variables (``$PETAKIT5D_DEMO_DATA/...``)

    copy: copied into the work folder as ``name`` (functions may write next to their input);
    else the path itself is passed. content: hash file contents (else names and sizes).
    """

    kind = "data"

    def __init__(
        self, path: str, name: str | None = None, *, copy: bool = True, content: bool = False
    ):
        self.path, self.copy, self.content = path, copy, content
        self.name = name or Path(path.rstrip("/\\")).name
        self._h: str | None = None

    def local(self) -> Path:
        p = Path(os.path.expandvars(self.path))
        if "$" in str(p) or not p.exists():
            raise FileNotFoundError(
                f"data input {self.path!r} not found (is its environment variable set?)"
            )
        return p

    def hash(self) -> str:
        if self._h is None:
            p = self.local()
            if p.is_dir():
                self._h = tree_hash(p, self.content)
            else:
                from petakit5d.testing.compare import sha256_file

                self._h = hashlib.sha256(
                    b"pk1:file:"
                    + (sha256_file(p) if self.content else str(p.stat().st_size)).encode()
                ).hexdigest()
        return self._h


# ---------------------------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------------------------


@dataclasses.dataclass
class PyResult:
    outputs: list | None
    error: BaseException | None
    error_text: str
    warnings: list[dict]
    files: FileTree | None
    work: Path
    elapsed: float
    raw_outputs: list | None = None  # before out_roles (Python's form)


@dataclasses.dataclass(eq=False)
class Case:
    """one oracle case (see the module docstring and python/docs/testing.md)

    id: unique in its test module. matlab: function name (``feval``), or None with ``code``: a
    MATLAB snippet that reads ``in`` (cell of the arguments) and ``work`` and sets ``out`` (cell
    of the outputs). args: MATLAB arguments (Python values, see ``to_mvalue``, and
    placeholders). nargout: outputs requested. python: 'package.module:function' or a
    callable; py_args/py_kwargs: its arguments (default: args); py_out: maps the Python
    result to the list of outputs. inputs: name -> Value, File, Tree or Data. compare: a spec
    for every output, or {output name or index: spec}. files: compare the files written
    (default: when the case uses a work folder or file inputs); tree: how. expect_error:
    MATLAB must raise (an identifier string: that one). py_error: exception classes Python
    may raise when MATLAB raises. error_message: compare error messages too. helpers:
    folders of MATLAB helper files (relative to python/) put on MATLAB's path. device: 'cpu',
    'gpu' or 'gpu:<model>'. marks: extra pytest marks. outputs: output names. out_roles: one
    out role per output (see ``out_to_matlab``): Python's results in MATLAB's form.
    warnings: 'record' or 'ids' (Python must give warnings with the same identifiers).
    Arguments in args, py_args and py_kwargs may carry roles (Vol, Box, PerAxis, Index1,
    Range1, Dim, LinIndex): MATLAB and the key see the value, Python its Python form.
    """

    id: str
    matlab: str | None = None
    args: Sequence[Any] = ()
    nargout: int = 1
    python: str | Callable | None = None
    py_args: Sequence[Any] | None = None
    py_kwargs: Mapping[str, Any] | None = None
    py_out: Callable[[Any], Sequence[Any]] | None = None
    inputs: Mapping[str, Input] = dataclasses.field(default_factory=dict)
    compare: Spec | Mapping[int | str, Spec] = EXACT
    files: bool | None = None
    tree: TreeSpec = TreeSpec()
    expect_error: bool | str = False
    py_error: type | tuple = Exception
    error_message: bool = False
    code: str | None = None
    helpers: Sequence[str] = ()
    device: str = "cpu"
    marks: Sequence[str] = ()
    outputs: Sequence[str] | None = None
    out_roles: Sequence[Any] | None = None
    warnings: str = "record"
    note: str = ""

    def __post_init__(self):
        if not self.id or any(c in self.id for c in " \t\n/\\"):
            raise ValueError(f"case id {self.id!r}: no spaces or slashes")
        if (self.matlab is None) == (self.code is None):
            raise ValueError(f"case {self.id}: give either matlab (a function name) or code")
        if self.outputs is not None and len(self.outputs) != self.nargout:
            raise ValueError(
                f"case {self.id}: {len(self.outputs)} output names for nargout={self.nargout}"
            )
        if self.out_roles is not None and len(self.out_roles) != self.nargout:
            raise ValueError(
                f"case {self.id}: {len(self.out_roles)} out roles for nargout={self.nargout}"
            )
        for ph in _placeholders(self.args):
            if ph.what in ("input", "path") and ph.name not in self.inputs:
                raise ValueError(f"case {self.id}: argument uses unknown input {ph.name!r}")
            if ph.what == "input" and not isinstance(self.inputs[ph.name], Value):
                raise ValueError(
                    f"case {self.id}: In({ph.name!r}) needs a Value input; use InPath for files"
                )
            if ph.what == "path" and isinstance(self.inputs[ph.name], Value):
                raise ValueError(
                    f"case {self.id}: InPath({ph.name!r}) needs a File, Tree or Data input"
                )
        self._key: str | None = None
        self._module_globals: dict | None = None

    # names, specs, marks

    def output_names(self) -> list[str]:
        return (
            list(self.outputs)
            if self.outputs is not None
            else [f"out{i + 1}" for i in range(self.nargout)]
        )

    def spec_for(self, i: int) -> Spec:
        if isinstance(self.compare, Mapping):
            name = self.output_names()[i]
            return self.compare.get(name, self.compare.get(i, EXACT))
        return self.compare

    @property
    def captures_files(self) -> bool:
        if self.files is not None:
            return self.files
        return any(p.what in ("work", "path") for p in _placeholders(self.args)) or any(
            not isinstance(x, Value) for x in self.inputs.values()
        )

    def all_marks(self) -> list[str]:
        m = list(self.marks)
        if self.device.startswith("gpu") and "gpu" not in m:
            m.append("gpu")
        if any(isinstance(x, Data) for x in self.inputs.values()) and "demodata" not in m:
            m.append("demodata")
        return m

    @property
    def module_file(self) -> str | None:
        return (self._module_globals or {}).get("__file__")

    # the key: what determines MATLAB's result

    def key(self, harness: Path | None = None) -> str:
        if self._key is None:
            spec = {
                "format": FORMAT,
                "matlab": self.matlab or "",
                "code": self.code or "",
                "nargout": self.nargout,
                "args": vhash(self._args_mvalue(lambda ph: MRef("", self.inputs[ph.name].hash()))),
                "inputs": {n: f"{x.kind}:{x.hash()}" for n, x in sorted(self.inputs.items())},
                "names": {n: getattr(x, "name", "") for n, x in sorted(self.inputs.items())},
                "files": self.captures_files,
                "device": self.device,
                "helpers": _helpers_hash(self.helpers),
                "harness": _harness_file_hash(self.matlab, harness),
            }
            # InPath(key) of an input whose name differs from its key: MATLAB resolved the key
            # until the harness resolved the name (results made before are not reused)
            renamed = {
                p.name: self.inputs[p.name].name
                for p in _placeholders(self.args)
                if p.what == "path" and getattr(self.inputs[p.name], "name", p.name) != p.name
            }
            if renamed:
                spec["paths"] = renamed
            self._key = hashlib.sha256(
                json.dumps(spec, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
        return self._key

    def _args_mvalue(self, input_ref: Callable[[MPlaceholder], MValue]) -> MCell:
        args = _strip_roles(list(self.args))  # roles change neither MATLAB's call nor the key
        mv = to_mvalue(args) if args else MCell((1, 0), [])
        if not isinstance(mv, MCell):
            mv = MCell((1, len(args)), [to_mvalue(a) for a in args])
        if len(mv.items) != len(args):
            mv = MCell((1, len(args)), [to_mvalue(a) for a in args])
        return _subst(mv, input_ref)

    def matlab_args(self, input_path: Callable[[str, Input], str]) -> MCell:
        """the arguments for MATLAB: inputs as references to value files, data paths filled in"""

        def ref(ph: MPlaceholder) -> MValue:
            return MRef(input_path(ph.name, self.inputs[ph.name]), self.inputs[ph.name].hash())

        mv = self._args_mvalue(ref)

        def data_paths(v: MValue) -> MValue:
            # InPath(key): MATLAB gets the path of the input's copy, which has the input's name
            if isinstance(v, MPlaceholder) and v.what == "path":
                x = self.inputs[v.name]
                if isinstance(x, Data) and not x.copy:
                    p = str(x.local())
                    return mchar(p + ("/" if x.local().is_dir() else ""))
                return MPlaceholder("path", x.name)
            return v

        return _walk(mv, data_paths)

    def check_inputs(self) -> list[str]:
        """problems of the generated inputs: each is made a second time and must give the
        same bytes (otherwise a reference cannot be matched when they are made again, as CI does)"""
        out = []
        for name, x in sorted(self.inputs.items()):
            if isinstance(x, Value):
                h = x.hash()
                again = vhash(to_mvalue(x.make()))
            elif isinstance(x, (File, Tree)):
                h = x.hash()
                fresh = type(x)(x.make, x.name)
                again = fresh.hash()
            else:
                continue
            if again != h:
                out.append(
                    f"case {self.id}: input {name!r} is not deterministic (made twice, the bytes "
                    "differ; "
                    "files must not hold times, uuids or random data)"
                )
        return out

    # the Python side

    def python_callable(self) -> Callable:
        if callable(self.python):
            return self.python
        if not isinstance(self.python, str) or ":" not in self.python:
            raise ValueError(
                f'case {self.id}: python must be "package.module:function" or a callable'
            )
        mod, _, name = self.python.partition(":")
        obj = importlib.import_module(mod)
        for part in name.split("."):
            obj = getattr(obj, part)
        return obj

    def materialize(self, work: Path) -> None:
        """copy the file, folder and data inputs into the work folder"""
        for x in self.inputs.values():
            if isinstance(x, File):
                shutil.copyfile(x.local(), work / x.name)
            elif isinstance(x, Tree):
                shutil.copytree(x.local(), work / x.name)
            elif isinstance(x, Data) and x.copy:
                src = x.local()
                if src.is_dir():
                    shutil.copytree(src, work / x.name)
                else:
                    shutil.copyfile(src, work / x.name)

    def python_args(self, work: Path) -> tuple[list, dict]:
        sizes = _Sizes(self)

        def fill(v: Any) -> Any:
            if isinstance(v, Role):
                return role_to_python(v, fill(v.value), sizes)
            if isinstance(v, MPlaceholder):
                if v.what == "input":
                    return self.inputs[v.name].python_value()
                if v.what == "work":
                    return work.as_posix() + "/" + v.sub
                x = self.inputs[v.name]
                if isinstance(x, Data) and not x.copy:
                    p = x.local()
                else:
                    p = work / x.name
                return p.as_posix() + ("/" if p.is_dir() else "")
            if isinstance(v, list):
                return [fill(y) for y in v]
            if isinstance(v, tuple):
                return tuple(fill(y) for y in v)
            if isinstance(v, dict):
                return {k: fill(y) for k, y in v.items()}
            return v

        args = self.args if self.py_args is None else self.py_args
        return [fill(a) for a in args], {k: fill(v) for k, v in (self.py_kwargs or {}).items()}

    def run_python(self, tmp_dir: str | os.PathLike) -> PyResult:
        """run the Python side in a fresh work folder under tmp_dir"""
        work = Path(tmp_dir) / "work"
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        self.materialize(work)
        fn = self.python_callable()
        args, kwargs = self.python_args(work)
        before = snapshot(work)
        cwd = os.getcwd()
        err, etext, res = None, "", None
        t0 = time.perf_counter()
        os.chdir(work)
        try:
            with warnings.catch_warnings(record=True) as wl:
                warnings.simplefilter("always")
                try:
                    res = fn(*args, **kwargs)
                except Exception as e:  # the case's own exception is a result
                    err, etext = e, traceback.format_exc()
        finally:
            os.chdir(cwd)
        elapsed = time.perf_counter() - t0
        wlist = [
            {
                "identifier": getattr(w.message, "identifier", "") or "",
                "message": str(w.message),
                "category": w.category.__name__,
            }
            for w in wl
        ]
        outs = raw = None
        if err is None:
            raw = self._outputs(res)
            outs = self.map_outputs(raw)
        files = None
        if self.captures_files:
            rels, deleted = changed_files(before, snapshot(work))
            files = scan_tree(work, rels, deleted=deleted)
        return PyResult(outs, err, etext, wlist, files, work, elapsed, raw)

    def map_outputs(self, outs: list) -> list:
        """Python's results in MATLAB's form, as out_roles say"""
        if self.out_roles is None:
            return outs
        sizes = _Sizes(self)
        return [out_to_matlab(r, v, sizes) for r, v in zip(self.out_roles, outs, strict=True)]

    def _outputs(self, res: Any) -> list:
        if self.py_out is not None:
            outs = list(self.py_out(res))
        elif self.nargout == 0:
            outs = []
        elif self.nargout == 1:
            outs = [res]
        else:
            if not isinstance(res, (tuple, list)):
                raise TypeError(
                    f"case {self.id}: nargout={self.nargout} needs a tuple from Python, got "
                    f"{type(res).__name__}"
                )
            outs = list(res)
        if len(outs) != self.nargout:
            raise TypeError(
                f"case {self.id}: Python gave {len(outs)} outputs, nargout is {self.nargout}"
            )
        return outs


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _placeholders(x: Any) -> list[MPlaceholder]:
    out = []
    if isinstance(x, Role):
        out.extend(_placeholders(x.value))
    elif isinstance(x, MPlaceholder):
        out.append(x)
    elif isinstance(x, (list, tuple)):
        for y in x:
            out.extend(_placeholders(y))
    elif isinstance(x, dict):
        for y in x.values():
            out.extend(_placeholders(y))
    elif isinstance(x, MCell) and x.items:
        for y in x.items:
            out.extend(_placeholders(y))
    elif isinstance(x, MStruct) and x.items:
        for el in x.items:
            for y in el.values():
                out.extend(_placeholders(y))
    return out


def _walk(v: MValue, fn: Callable[[MValue], MValue]) -> MValue:
    if isinstance(v, MCell) and v.items is not None:
        return MCell(v.size, [_walk(x, fn) for x in v.items])
    if isinstance(v, MStruct) and v.items is not None:
        return MStruct(
            v.size, v.fields, [{k: _walk(x, fn) for k, x in el.items()} for el in v.items]
        )
    return fn(v)


def _subst(v: MValue, input_ref: Callable[[MPlaceholder], MValue]) -> MValue:
    return _walk(
        v, lambda x: input_ref(x) if isinstance(x, MPlaceholder) and x.what == "input" else x
    )


def _python_root() -> Path:
    """the python/ folder of the repo (oracle.python_root: also when petakit5d is installed)"""
    from petakit5d.testing.oracle import python_root

    return python_root()


def _helpers_hash(helpers: Sequence[str]) -> dict[str, str]:
    """{folder: tree_hash} of a case's helper folders (relative to python/); a folder that is
    not there is an error, not the hash of an empty folder"""
    out = {}
    for h in sorted(helpers):
        p = _python_root() / h
        if not p.is_dir():
            raise FileNotFoundError(f"helper folder {h} not found at {p}")
        out[h] = tree_hash(p)
    return out


def _harness_file_hash(func: str | None, harness: Path | None) -> str:
    """hash of a harness MATLAB function a case calls (pkpy_oracle.selftest.*), else ''"""
    if not func or not func.startswith("pkpy_oracle."):
        return ""
    from petakit5d.testing.oracle import harness_dir

    root = harness or harness_dir()
    parts = func.split(".")
    p = root.joinpath(*[f"+{x}" for x in parts[:-1]], parts[-1] + ".m")
    if not p.is_file():
        raise FileNotFoundError(f"harness function {func} not found at {p}")
    return hashlib.sha256(p.read_bytes()).hexdigest()


def snapshot(root: Path) -> dict[str, tuple[bool, int, int]]:
    """relative path -> (is folder, size, mtime in ns) of everything under root"""
    out = {}
    for dp, dns, fns in os.walk(root):
        d = Path(dp)
        for n in dns:
            p = d / n
            out[p.relative_to(root).as_posix()] = (True, 0, 0)
        for n in fns:
            p = d / n
            st = p.stat()
            out[p.relative_to(root).as_posix()] = (False, st.st_size, st.st_mtime_ns)
    return out


def _zarr_root(rel: str) -> str | None:
    parts = rel.split("/")
    for i, p in enumerate(parts):
        if p.lower().endswith(".zarr"):
            return "/".join(parts[: i + 1])
    return None


def changed_files(before: dict, after: dict) -> tuple[list[str], list[str]]:
    """what a call wrote: new or changed files (a zarr folder as a whole), new empty folders;
    and the files it removed"""
    out, deleted = set(), set()
    parents = {rel[: rel.rfind("/")] for rel in after if "/" in rel}
    for rel, st in after.items():
        z = _zarr_root(rel)
        if st[0]:
            if z == rel:
                if rel not in before:
                    out.add(rel)
            elif z is None and rel not in before and rel not in parents:
                out.add(rel)
            continue
        if before.get(rel) != st:
            out.add(z or rel)
    for rel, st in before.items():
        if rel in after:
            continue
        z = _zarr_root(rel)
        if z is not None:
            if z in after:
                out.add(z)
            elif z == rel:
                deleted.add(rel)
        elif not st[0]:
            deleted.add(rel)
    return sorted(out), sorted(deleted)
