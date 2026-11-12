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
    'gpu' or 'gpu:<model>'. marks: extra pytest marks. outputs: output names.
    warnings: 'record' or 'ids' (Python must give warnings with the same identifiers).
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
                "helpers": {h: tree_hash(_python_root() / h) for h in sorted(self.helpers)},
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
        mv = to_mvalue(list(self.args)) if self.args else MCell((1, 0), [])
        if not isinstance(mv, MCell):
            mv = MCell((1, len(self.args)), [to_mvalue(a) for a in self.args])
        if len(mv.items) != len(self.args):
            mv = MCell((1, len(self.args)), [to_mvalue(a) for a in self.args])
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
        def fill(v: Any) -> Any:
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
        outs = None
        if err is None:
            outs = self._outputs(res)
        files = None
        if self.captures_files:
            rels, deleted = changed_files(before, snapshot(work))
            files = scan_tree(work, rels, deleted=deleted)
        return PyResult(outs, err, etext, wlist, files, work, elapsed)

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
    if isinstance(x, MPlaceholder):
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
    """the python/ folder of the repo"""
    return Path(__file__).resolve().parents[3]


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
