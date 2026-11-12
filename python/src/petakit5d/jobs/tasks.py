"""Tasks of the job framework as records, and the call folder that holds them.

A task is a record, never code text:

* a Python task names a callable by its import path and gives JSON arguments,
  ``{"func": "package.module:function", "args": [...], "kwargs": {...}}``, with an optional
  ``"path"``: folders to add to ``sys.path`` if the module can't be imported without them
  (a user module outside an installed package);
* a bash task is a command, ``{"bash": "command"}`` (``language='bash'``).

Arguments are made canonical when a task is made (tuples become lists, paths strings, NumPy
values Python numbers through ``tolist``) and pass through a JSON round trip, so that a task
gets the same values in the caller's process (local mode) as in a cluster job, which reads
them from the call folder. Floats keep every bit (shortest repr), and inf and NaN are kept
(Python's JSON extension).

The call folder (written by the cluster backend; local runs need none)::

    <folder>/call.json     metadata, settings and the byte offset of each task line
    <folder>/tasks.ndjson  one line per task: {"index", "input", "output", "task"}

``call.json`` is written last, by rename, so a job that sees it sees every task. A job runs
``python -m petakit5d.jobs.exec <folder> <batch>``; batch ``b`` (1-based) holds tasks
``(b - 1) * task_batch_num + 1`` to ``min(b * task_batch_num, n)``.

This module imports only the standard library: the job entry point reads it at start-up.
"""

from __future__ import annotations

import importlib
import json
import math
import os
import sys
import time
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any

CALL_FORMAT = "petakit5d.jobs.call"
CALL_VERSION = 1
CALL_FILE = "call.json"
TASKS_FILE = "tasks.ndjson"

LANGUAGES = ("python", "bash")


class Task:
    """A Python task: ``func`` is ``'package.module:qualified.name'``; ``args`` and ``kwargs``
    are JSON-serializable (made canonical here); ``path`` lists import roots to add if the
    module isn't importable without them."""

    __slots__ = ("func", "args", "kwargs", "path")

    def __init__(
        self,
        func: str,
        args: Sequence[Any] = (),
        kwargs: Mapping[str, Any] | None = None,
        path: Sequence[str] = (),
    ) -> None:
        if not isinstance(func, str) or func.count(":") != 1:
            raise ValueError(f"func must be 'package.module:function', not {func!r}")
        module, qualname = func.split(":")
        if not module or not qualname or any(not p.isidentifier() for p in module.split(".")):
            raise ValueError(f"func must be 'package.module:function', not {func!r}")
        if any(not p.isidentifier() for p in qualname.split(".")):
            raise ValueError(f"{func!r}: {qualname!r} is not importable by name")
        if isinstance(args, (str, bytes, Mapping)):
            raise TypeError("args must be a sequence of arguments")
        self.func = func
        self.args = canonical(list(args), "args")
        self.kwargs = canonical(dict(kwargs or {}), "kwargs")
        if isinstance(path, str):
            path = [path]
        self.path = [os.fspath(p) for p in path]

    @classmethod
    def of(cls, func: Callable[..., Any], /, *args: Any, **kwargs: Any) -> Task:
        """The task that calls ``func(*args, **kwargs)``; ``func`` must be importable by its
        module and qualified name (not a lambda or a nested function). A function of the
        script run as ``__main__`` is imported by the script's name, which runs the script's
        top-level code in the job (guard it with ``if __name__ == "__main__":``, as for
        multiprocessing)."""
        import functools

        if isinstance(func, functools.partial):
            return cls.of(func.func, *func.args, *args, **{**func.keywords, **kwargs})
        module = getattr(func, "__module__", None)
        qualname = getattr(func, "__qualname__", None)
        if not module or not qualname or "<" in qualname:
            raise ValueError(f"{func!r} can't be imported by name; define it at module level")
        path: list[str] = []
        if module == "__main__":
            module, root = _main_module_name()
            path.append(root)
        elif module.split(".")[0] != "petakit5d":
            root = _module_root(module)
            if root:
                path.append(root)
        return cls(f"{module}:{qualname}", args, kwargs, path)

    def record(self) -> dict[str, Any]:
        rec: dict[str, Any] = {"func": self.func, "args": self.args, "kwargs": self.kwargs}
        if self.path:
            rec["path"] = list(self.path)
        return rec

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Task) and self.record() == other.record()

    def __repr__(self) -> str:
        return f"Task({self.func!r}, {self.args!r}, {self.kwargs!r})"


def canonical(value: Any, where: str = "value") -> Any:
    """``value`` as JSON data: None, bool, int, float, str, lists and dicts with str keys.
    Tuples become lists, ``os.PathLike`` strings, objects with ``tolist`` (NumPy arrays and
    scalars) their ``tolist()``; anything else is a TypeError naming where it was."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return float(value)
    if isinstance(value, (list, tuple)):
        return [canonical(v, f"{where}[{i}]") for i, v in enumerate(value)]
    if isinstance(value, Mapping):
        out = {}
        for k, v in value.items():
            if not isinstance(k, str):
                raise TypeError(f"{where}: dict keys must be str, not {type(k).__name__}")
            out[k] = canonical(v, f"{where}[{k!r}]")
        return out
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        return canonical(tolist(), where)
    raise TypeError(f"{where}: {type(value).__name__} is not JSON data (task arguments must be)")


def task_record(task: Any, language: str) -> dict[str, Any]:
    """The record of one task given to run_tasks: a Task, a record dict, a callable (no
    arguments) or a functools.partial for Python; a command string for bash. The record has
    been through a JSON round trip."""
    if language == "bash":
        if isinstance(task, Mapping) and set(task) == {"bash"}:
            task = task["bash"]
        if not isinstance(task, str):
            raise TypeError(f"a bash task is a command string, not {type(task).__name__}")
        rec: dict[str, Any] = {"bash": task}
    elif isinstance(task, Task):
        rec = task.record()
    elif isinstance(task, Mapping):
        if "bash" in task:
            raise TypeError("a bash task needs language='bash'")
        unknown = set(task) - {"func", "args", "kwargs", "path"}
        if unknown:
            raise ValueError(f"unknown task record keys {sorted(unknown)}")
        rec = Task(task["func"], task.get("args", ()), task.get("kwargs"), task.get("path", ()))
        rec = rec.record()
    elif isinstance(task, str):
        raise TypeError(f"a str task is a bash command and needs language='bash': {task!r}")
    elif callable(task):
        rec = Task.of(task).record()
    else:
        raise TypeError(f"a task is a Task, a record or a callable, not {type(task).__name__}")
    return json.loads(json.dumps(rec, allow_nan=True))


def describe(rec: Mapping[str, Any]) -> str:
    """A one-line text of a task record, as the framework prints its function strings."""
    if "bash" in rec:
        return rec["bash"]
    parts = [json.dumps(a) for a in rec.get("args", ())]
    parts += [f"{k}={json.dumps(v)}" for k, v in rec.get("kwargs", {}).items()]
    return f"{rec['func']}({', '.join(parts)})"


def resolve(rec: Mapping[str, Any]) -> Callable[..., Any]:
    """The callable of a Python task record. The record's import roots are appended to
    sys.path only if the module can't be imported without them."""
    module, qualname = rec["func"].split(":")
    try:
        obj: Any = importlib.import_module(module)
    except ModuleNotFoundError as err:
        roots = [p for p in rec.get("path", ()) if p not in sys.path]
        if not roots or not (err.name and module.startswith(err.name)):
            raise
        sys.path.extend(roots)
        obj = importlib.import_module(module)
    for name in qualname.split("."):
        obj = getattr(obj, name)
    return obj


def run_record(rec: Mapping[str, Any]) -> Any:
    """Calls a Python task record in this process and returns what the callable returns."""
    return resolve(rec)(*rec.get("args", ()), **rec.get("kwargs", {}))


def batch_indices(batch: int, task_batch_num: int, n_tasks: int) -> range:
    """0-based indices of the tasks of 1-based ``batch``."""
    return range((batch - 1) * task_batch_num, min(batch * task_batch_num, n_tasks))


class Call:
    """One call of run_tasks: the tasks with their inputs and outputs, and the effective
    settings (python parameter names; config file applied). This is what the backends
    (local loop, process pool, cluster loop) receive.

    ``inputs[i]`` is a path or a list of paths, ``outputs[i]`` a path (trailing separators
    removed, as MATLAB does), ``records[i]`` a task record.
    """

    def __init__(
        self,
        inputs: Sequence[str | Sequence[str]],
        outputs: Sequence[str],
        records: Sequence[Mapping[str, Any]],
        settings: Mapping[str, Any],
    ) -> None:
        if not len(inputs) == len(outputs) == len(records):
            raise ValueError(
                f"{len(inputs)} inputs, {len(outputs)} outputs and {len(records)} tasks: "
                "each task needs one input entry and one output path"
            )
        self.inputs = [p if isinstance(p, str) else [os.fspath(q) for q in p] for p in inputs]
        self.outputs = list(outputs)
        self.records = list(records)
        self.settings = dict(settings)
        self.folder: str | None = None

    @property
    def n(self) -> int:
        return len(self.records)

    @property
    def language(self) -> str:
        return self.settings.get("language", "python")

    @property
    def task_batch_num(self) -> int:
        return int(self.settings.get("task_batch_num", 1))

    @property
    def batch_num(self) -> int:
        return math.ceil(self.n / self.task_batch_num) if self.n else 0

    def batch(self, batch: int) -> range:
        return batch_indices(batch, self.task_batch_num, self.n)


def new_call_folder(parent: str, uuid: str) -> str:
    """Makes and returns ``<parent>/pk5d_call_<yyyyMMdd_HHmmss_SSS>_<uuid[:6]>`` (absolute)."""
    now = time.time()
    stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime(now)) + f"_{int(now * 1000) % 1000:03d}"
    folder = os.path.join(os.path.abspath(parent), f"pk5d_call_{stamp}_{uuid[:6]}")
    os.makedirs(folder)
    return folder


def write_call(call: Call, folder: str, extra: Mapping[str, Any] | None = None) -> str:
    """Writes the call's tasks and metadata into ``folder`` (made if needed) and returns its
    absolute path. ``extra`` adds keys to call.json (e.g. ``"runner"``: the import path of a
    batch runner that the entry point calls instead of running the batch itself)."""
    folder = os.path.abspath(folder)
    os.makedirs(folder, exist_ok=True)
    offsets = []
    pos = 0
    tmp = os.path.join(folder, f".{TASKS_FILE}.{os.getpid()}.tmp")
    with open(tmp, "wb") as fh:
        for i, rec in enumerate(call.records):
            line = {"index": i + 1, "input": call.inputs[i], "output": call.outputs[i], "task": rec}
            data = (json.dumps(line, allow_nan=True) + "\n").encode()
            offsets.append(pos)
            pos += len(data)
            fh.write(data)
    os.replace(tmp, os.path.join(folder, TASKS_FILE))
    meta = {
        "format": CALL_FORMAT,
        "version": CALL_VERSION,
        "n_tasks": call.n,
        "language": call.language,
        "task_batch_num": call.task_batch_num,
        "final_out_fullpath": call.settings.get("final_out_fullpath", ""),
        "bash_launch_str": call.settings.get("bash_launch_str", ""),
        "uuid": call.settings.get("uuid", ""),
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "python": sys.executable,
        "settings": canonical_settings(call.settings),
        "offsets": offsets,
    }
    if extra:
        meta.update(canonical(dict(extra), "extra"))
    tmp = os.path.join(folder, f".{CALL_FILE}.{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, allow_nan=True)
    os.replace(tmp, os.path.join(folder, CALL_FILE))
    call.folder = folder
    return folder


def canonical_settings(settings: Mapping[str, Any]) -> dict[str, Any]:
    """The settings that are JSON data (others, such as objects, are left out)."""
    out = {}
    for k, v in settings.items():
        try:
            out[k] = canonical(v, k)
        except TypeError:
            continue
    return out


class CallFolder:
    """Reads a call folder: ``meta`` (call.json) and task lines by index."""

    def __init__(self, folder: str) -> None:
        self.folder = folder
        with open(os.path.join(folder, CALL_FILE), encoding="utf-8") as fh:
            self.meta: dict[str, Any] = json.load(fh)
        if self.meta.get("format") != CALL_FORMAT:
            raise ValueError(f"{folder} is not a call folder of the job framework")
        if self.meta.get("version", 0) > CALL_VERSION:
            raise ValueError(
                f"{folder}: call folder version {self.meta['version']} is newer "
                f"than this petakit5d reads ({CALL_VERSION})"
            )

    @property
    def n(self) -> int:
        return int(self.meta["n_tasks"])

    @property
    def language(self) -> str:
        return self.meta["language"]

    def batch(self, batch: int) -> range:
        return batch_indices(batch, int(self.meta["task_batch_num"]), self.n)

    def lines(self, indices: Iterable[int]) -> list[dict[str, Any]]:
        """The task lines (``index``, ``input``, ``output``, ``task``) at 0-based indices."""
        offsets = self.meta["offsets"]
        out = []
        with open(os.path.join(self.folder, TASKS_FILE), "rb") as fh:
            for i in indices:
                fh.seek(offsets[i])
                out.append(json.loads(fh.readline()))
        return out


def read_call(folder: str) -> Call:
    """The Call written to ``folder`` (all its tasks)."""
    cf = CallFolder(folder)
    lines = cf.lines(range(cf.n))
    call = Call(
        [x["input"] for x in lines],
        [x["output"] for x in lines],
        [x["task"] for x in lines],
        cf.meta.get("settings", {}),
    )
    call.folder = folder
    return call


def _main_module_name() -> tuple[str, str]:
    # the import name and root of the script run as __main__
    main = sys.modules.get("__main__")
    spec = getattr(main, "__spec__", None)
    if spec is not None and spec.name and spec.name != "__main__":
        return spec.name, _module_root(spec.name, getattr(main, "__file__", None)) or ""
    file = getattr(main, "__file__", None)
    if not file:
        raise ValueError(
            "a function defined interactively can't be imported by a job; put it in a module"
        )
    name = os.path.splitext(os.path.basename(file))[0]
    if not name.isidentifier():
        raise ValueError(f"the script {file} has no importable module name")
    return name, os.path.dirname(os.path.abspath(file))


def _module_root(module: str, file: str | None = None) -> str | None:
    # the sys.path folder a module was imported from (None for namespace packages)
    if file is None:
        file = getattr(sys.modules.get(module), "__file__", None)
    if not file:
        return None
    root = os.path.dirname(os.path.abspath(file))
    depth = module.count(".") + (os.path.basename(file).startswith("__init__.") and 1)
    for _ in range(depth):
        root = os.path.dirname(root)
    return root
