"""Tasks that share a process: what a task can leave behind for the next one, given back after
each task (the task runners), and a harness that shows that tasks run one after another in one
process give what fresh processes give, bit for bit (MATLAB's tests/isolation and the reset of
task_runner.m).

``ProcessState.capture()`` takes the state of this process before its first task;
``state.restore()`` gives it back after a task and returns the names of what the task had
changed. It covers, MATLAB's reset item by item where Python has one:

  folder          the current folder
  environment     variables added, changed or removed (also for C code: os.environ calls putenv)
  path            sys.path, sys.meta_path, sys.path_hooks (a task record's import roots)
  warnings        the warning filters (and the "once" registries, so that warnings show again
                  as in a new process), showwarning and formatwarning
  threads         the CPUs the process may use (affinity); the thread variables of the numeric
                  libraries are in the environment
  fp              the floating-point environment of this thread (rounding mode, flush to zero,
                  exception masks; C's fegetenv), NumPy's error state and its print options
                  (with NumPy 2 these and decimal's context are context variables: a task runs
                  in a copy of the process's context, which goes with it)
  random          Python's and NumPy's global random numbers, seeded anew as in a new process
                  (not given back: two tasks would then draw the same numbers)
  umask, locale, signal handlers, a pending alarm, stdin/stdout/stderr, the recursion limit,
  the switch interval, gc settings, exception hooks, tempfile's folder, decimal's context,
  the root logger (handlers, level, logging.disable), and CuPy's memory pools (freed)

What no reset undoes is left to the runner's guards, which end the process: memory a task
leaves (resident memory above the warm process), threads and open files it leaves, a slow
build-up (a process runs tasks for 10 minutes), and a task that failed. Modules stay
imported (the point of running tasks in one process); their caches must be exact (keyed by
contents), which the harness checks for the task functions it is given.

The harness (``run_harness``, or ``python -m petakit5d.jobs.isolation run <spec> <work>
<results>``): a spec names task records that write into ``<work>/out`` and a task that
prepares the inputs. Each run starts from a new work folder (always the same absolute path,
since task records hold absolute paths): first every task in a process of its own ("fresh"),
then all tasks in one process in random orders with the runner's restore between them; each
run's ``out`` folder is kept in ``<results>/<run>`` and compared with the fresh one (files,
tiff and zarr data and metadata through ``petakit5d.testing.compare``).
"""

from __future__ import annotations

import contextvars
import gc
import json
import os
import signal
import sys
import threading
import time
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    "HarnessResult",
    "ProcessState",
    "TaskResult",
    "main",
    "run_harness",
    "run_in_process",
]

_LINUX = sys.platform.startswith("linux")
_SIGNALS = tuple(
    getattr(signal, n)
    for n in ("SIGTERM", "SIGINT", "SIGHUP", "SIGUSR1", "SIGUSR2", "SIGALRM")
    if hasattr(signal, n)
)
# NumPy's error state in a new process
_NUMPY_ERR = {"divide": "warn", "over": "warn", "under": "ignore", "invalid": "warn"}


# --- the floating-point environment of this thread (C's fenv) ---------------------------------


class _FPEnv:
    """fegetenv/fesetenv through ctypes; ``control(buf)`` is the part that changes results
    (rounding mode, flush to zero, exception masks), without the sticky exception flags"""

    _lib: Any = None
    _ok: bool | None = None

    @classmethod
    def _load(cls) -> bool:
        if cls._ok is None:
            cls._ok = False
            try:
                import ctypes

                name = {"linux": "libm.so.6", "darwin": "libSystem.B.dylib"}.get(sys.platform)
                if name:
                    cls._lib = ctypes.CDLL(name)
                    cls._lib.fegetenv.argtypes = [ctypes.c_void_p]
                    cls._lib.fesetenv.argtypes = [ctypes.c_void_p]
                    cls._ok = True
            except (OSError, AttributeError):
                cls._ok = False
        return bool(cls._ok)

    @classmethod
    def get(cls) -> bytes | None:
        if not cls._load():
            return None
        import ctypes

        buf = ctypes.create_string_buffer(64)
        if cls._lib.fegetenv(buf) != 0:
            return None
        return buf.raw

    @classmethod
    def set(cls, raw: bytes) -> None:
        if not cls._load():
            return
        import ctypes

        buf = ctypes.create_string_buffer(raw, 64)
        cls._lib.fesetenv(buf)

    @staticmethod
    def control(raw: bytes) -> bytes:
        import platform

        m = platform.machine().lower()
        if m in ("x86_64", "amd64"):
            if sys.platform == "darwin":
                # {u16 control, u16 status, u32 mxcsr, ...}
                mx = int.from_bytes(raw[4:8], "little") & ~0x3F
                return raw[0:2] + mx.to_bytes(4, "little")
            # glibc: x87 control word first, mxcsr at byte 28
            mx = int.from_bytes(raw[28:32], "little") & ~0x3F
            return raw[0:2] + mx.to_bytes(4, "little")
        if m in ("aarch64", "arm64"):
            return raw[8:16] if sys.platform == "darwin" else raw[0:4]
        return raw[:32]


# --- the state --------------------------------------------------------------------------------


@dataclass
class ProcessState:
    """What a task can change in this process (see the module docstring)."""

    cwd: str | None
    environ: dict[str, str]
    path: list[str]
    meta_path: list[Any]
    path_hooks: list[Any]
    filters: list[Any]
    showwarning: Any
    formatwarning: Any
    affinity: set[int] | None
    umask: int | None
    locale: str | None
    handlers: dict[int, Any]
    stdio: tuple[Any, Any, Any]
    recursion: int
    switch: float
    gc_state: tuple[bool, tuple[int, ...], int]
    hooks: tuple[Any, Any, Any]
    tempdir: Any
    fpenv: bytes | None
    logging_state: tuple[Any, ...] | None = None
    decimal_context: Any = None
    numpy_state: tuple[Any, ...] | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def capture(cls) -> ProcessState:
        main = threading.current_thread() is threading.main_thread()
        try:
            cwd = os.getcwd()
        except OSError:
            cwd = None
        affinity = None
        if hasattr(os, "sched_getaffinity"):
            try:
                affinity = set(os.sched_getaffinity(0))
            except OSError:
                pass
        umask = None
        if hasattr(os, "umask"):
            umask = os.umask(0o022)
            os.umask(umask)
        handlers = {}
        if main:
            for s in _SIGNALS:
                try:
                    handlers[s] = signal.getsignal(s)
                except (ValueError, OSError):
                    pass
        tempfile = sys.modules.get("tempfile")
        st = cls(
            cwd=cwd,
            environ=dict(os.environ),
            path=list(sys.path),
            meta_path=list(sys.meta_path),
            path_hooks=list(sys.path_hooks),
            filters=list(warnings.filters),
            showwarning=warnings.showwarning,
            formatwarning=warnings.formatwarning,
            affinity=affinity,
            umask=umask,
            locale=_get_locale(),
            handlers=handlers,
            stdio=(sys.stdin, sys.stdout, sys.stderr),
            recursion=sys.getrecursionlimit(),
            switch=sys.getswitchinterval(),
            gc_state=(gc.isenabled(), gc.get_threshold(), gc.get_debug()),
            hooks=(sys.excepthook, threading.excepthook, sys.displayhook),
            tempdir=getattr(tempfile, "tempdir", None),
            fpenv=_FPEnv.get(),
        )
        st.logging_state = _logging_state()
        st.decimal_context = _decimal_context()
        st.numpy_state = _numpy_state()
        # (a timer already pending is not a task's: the process's own)
        st.extra["timer"] = (
            hasattr(signal, "getitimer") and signal.getitimer(signal.ITIMER_REAL)[0] > 0
        )
        return st

    def restore(self) -> list[str]:
        """gives the state back; returns what had changed (in the order of the docstring)"""
        changed: list[str] = []
        main = threading.current_thread() is threading.main_thread()
        if main:
            # (a pending alarm or timer of the task would go off in the next one)
            if hasattr(signal, "setitimer") and not self.extra.get("timer"):
                try:
                    if signal.setitimer(signal.ITIMER_REAL, 0)[0] > 0:
                        changed.append("alarm")
                except (OSError, ValueError):
                    pass
            for s, h in self.handlers.items():
                try:
                    if signal.getsignal(s) is not h and h is not None:
                        signal.signal(s, h)
                        changed.append(f"signal {s}")
                except (ValueError, OSError, TypeError):
                    pass
        if self.cwd is not None:
            try:
                now = os.getcwd()
            except OSError:
                now = None
            if now != self.cwd:
                try:
                    os.chdir(self.cwd)
                except OSError:
                    pass
                changed.append("folder")
        if _restore_environ(self.environ):
            changed.append("environment")
        if sys.path != self.path:
            sys.path[:] = self.path
            changed.append("path")
        if sys.meta_path != self.meta_path or sys.path_hooks != self.path_hooks:
            sys.meta_path[:] = self.meta_path
            sys.path_hooks[:] = self.path_hooks
            sys.path_importer_cache.clear()
            changed.append("import system")
        if (
            warnings.filters != self.filters
            or warnings.showwarning is not self.showwarning
            or warnings.formatwarning is not self.formatwarning
        ):
            warnings.filters[:] = self.filters
            warnings.showwarning = self.showwarning
            warnings.formatwarning = self.formatwarning
            changed.append("warnings")
        # (the registries of warnings shown "once" or by "default" start empty, as in a new
        # process: a change of the filters' version empties them)
        mutated = getattr(warnings, "_filters_mutated", None) or getattr(
            warnings, "_filters_mutated_lock_held", None
        )
        if mutated is not None:
            try:
                mutated()
            except Exception:
                pass
        if self.affinity is not None:
            try:
                if set(os.sched_getaffinity(0)) != self.affinity:
                    os.sched_setaffinity(0, self.affinity)
                    changed.append("threads")
            except OSError:
                changed.append("threads (not restored)")
        if self.fpenv is not None:
            now_env = _FPEnv.get()
            if now_env is not None and _FPEnv.control(now_env) != _FPEnv.control(self.fpenv):
                _FPEnv.set(self.fpenv)
                changed.append("fp")
        if self.umask is not None:
            u = os.umask(self.umask)
            if u != self.umask:
                changed.append("umask")
        loc = _get_locale()
        if loc != self.locale and self.locale is not None:
            try:
                import locale

                locale.setlocale(locale.LC_ALL, self.locale)
            except Exception:
                pass
            changed.append("locale")
        if (sys.stdin, sys.stdout, sys.stderr) != self.stdio:
            sys.stdin, sys.stdout, sys.stderr = self.stdio
            changed.append("stdio")
        if sys.getrecursionlimit() != self.recursion:
            sys.setrecursionlimit(self.recursion)
            changed.append("recursion limit")
        if sys.getswitchinterval() != self.switch:
            sys.setswitchinterval(self.switch)
            changed.append("switch interval")
        if (gc.isenabled(), gc.get_threshold(), gc.get_debug()) != self.gc_state:
            (gc.enable if self.gc_state[0] else gc.disable)()
            gc.set_threshold(*self.gc_state[1])
            gc.set_debug(self.gc_state[2])
            changed.append("gc")
        if (sys.excepthook, threading.excepthook, sys.displayhook) != self.hooks:
            sys.excepthook, threading.excepthook, sys.displayhook = self.hooks
            changed.append("hooks")
        tempfile = sys.modules.get("tempfile")
        if tempfile is not None and tempfile.tempdir != self.tempdir:
            tempfile.tempdir = self.tempdir
            changed.append("tempfile")
        if _restore_logging(self.logging_state):
            changed.append("logging")
        if _restore_decimal(self.decimal_context):
            changed.append("decimal")
        changed += _restore_numpy(self.numpy_state)
        _reseed_random()
        _free_gpu_pools()
        return changed


def _get_locale() -> str | None:
    import locale

    try:
        return locale.setlocale(locale.LC_ALL)
    except Exception:
        return None


def _restore_environ(saved: Mapping[str, str]) -> bool:
    now = os.environ
    if now == saved:
        return False
    for k in [k for k in now if k not in saved]:
        del now[k]
    for k, v in saved.items():
        if now.get(k) != v:
            now[k] = v
    return True


def _logging_state() -> tuple[Any, ...] | None:
    logging = sys.modules.get("logging")
    if logging is None:
        return None
    root = logging.root
    return (list(root.handlers), root.level, root.manager.disable)


def _restore_logging(saved: tuple[Any, ...] | None) -> bool:
    logging = sys.modules.get("logging")
    if logging is None:
        return False
    root = logging.root
    if saved is None:
        # (imported by a task: as logging starts)
        saved = ([], logging.WARNING, 0)
    if (root.handlers, root.level, root.manager.disable) == (saved[0], saved[1], saved[2]):
        return False
    root.handlers[:] = saved[0]
    root.setLevel(saved[1])
    logging.disable(saved[2])
    return True


def _decimal_context() -> Any:
    decimal = sys.modules.get("decimal")
    return None if decimal is None else decimal.getcontext().copy()


def _restore_decimal(saved: Any) -> bool:
    decimal = sys.modules.get("decimal")
    if decimal is None:
        return False
    if saved is None:
        saved = decimal.Context()
    now = decimal.getcontext()

    def settings(c: Any) -> tuple[Any, ...]:
        # (not the flags: they are signals of the last operations, cleared anyway)
        return (c.prec, c.rounding, c.Emin, c.Emax, c.capitals, c.clamp, dict(c.traps))

    changed = settings(now) != settings(saved)
    decimal.setcontext(saved.copy())
    return changed


def _numpy_state() -> tuple[Any, ...] | None:
    np = sys.modules.get("numpy")
    if np is None or not hasattr(np, "geterr"):
        return None
    return (np.geterr(), np.geterrcall(), np.get_printoptions())


def _numpy_default_printoptions(np: Any) -> dict[str, Any] | None:
    # NumPy's print options in a new process (as np.get_printoptions gives them)
    mod = sys.modules.get("numpy._core.printoptions")
    d = getattr(mod, "default_format_options_dict", None)
    if not isinstance(d, dict):
        return None
    out = dict(d)
    if out.get("legacy") == sys.maxsize:
        out["legacy"] = False
    return out


def _restore_numpy(saved: tuple[Any, ...] | None) -> list[str]:
    np = sys.modules.get("numpy")
    if np is None or not hasattr(np, "geterr"):
        return []
    if saved is None:
        # (imported by a task: as NumPy starts)
        saved = (dict(_NUMPY_ERR), None, _numpy_default_printoptions(np))
    out = []
    if np.geterr() != saved[0] or np.geterrcall() is not saved[1]:
        np.seterr(**saved[0])
        np.seterrcall(saved[1])
        out.append("numpy errors")
    if saved[2] is not None:
        now = np.get_printoptions()
        if now != saved[2]:
            keys = [k for k in saved[2] if k in now]
            np.set_printoptions(**{k: saved[2][k] for k in keys})
            out.append("numpy print options")
    return out


def _reseed_random() -> None:
    # as in a new process: seeded from the operating system (if a task imported them)
    random = sys.modules.get("random")
    if random is not None:
        random.seed()
    npr = sys.modules.get("numpy.random")
    if npr is not None and hasattr(npr, "seed"):
        try:
            npr.seed()
        except Exception:
            pass


def _free_gpu_pools() -> None:
    cp = sys.modules.get("cupy")
    if cp is None:
        return
    try:
        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass


# --- one task in this process ----------------------------------------------------------------


@dataclass
class TaskResult:
    """a task run in this process: its error (None if it returned; a SystemExit with code 0 or
    None counts as returned, as for a process of its own), its seconds and what it changed"""

    error: BaseException | None
    seconds: float
    changed: list[str]
    text: str = ""


def run_in_process(
    record: Mapping[str, Any],
    baseline: ProcessState | None,
    before: Callable[[], Any] | None = None,
    after: Callable[[Any, BaseException | None], None] | None = None,
) -> TaskResult:
    """Runs a Python task record in this process and gives the process state back to
    ``baseline`` after it (no restore without one; the task then also runs in the process's
    own context). ``before()`` runs just before the task and ``after(token, error)`` just
    after it, before the restore (the task metrics). A task's exception is printed with its
    traceback (stderr) and returned; BaseExceptions other than SystemExit (SIGTERM's
    Terminated, KeyboardInterrupt) propagate after ``after``."""
    import traceback

    from petakit5d.jobs.tasks import run_record

    token = before() if before else None
    t0 = time.perf_counter()
    err: BaseException | None = None
    text = ""
    try:
        if baseline is None:
            run_record(record)
        else:
            # (in a context of its own: what the task sets in context variables, such as
            # NumPy's print options and error state and decimal's context, stays there)
            contextvars.copy_context().run(run_record, record)
    except SystemExit as e:
        if e.code not in (None, 0):
            err = e
            text = traceback.format_exc()
    except Exception as e:
        err = e
        text = traceback.format_exc()
    except BaseException as e:
        if after:
            after(token, e)
        raise
    seconds = time.perf_counter() - t0
    if after:
        after(token, err)
    if err is not None:
        _flush()
        print(text, end="", file=sys.stderr, flush=True)
    changed = baseline.restore() if baseline is not None else []
    return TaskResult(err, seconds, changed, text)


def _flush() -> None:
    for s in (sys.stdout, sys.stderr):
        try:
            s.flush()
        except (AttributeError, OSError, ValueError):
            pass


# --- the harness ----------------------------------------------------------------------------


@dataclass
class HarnessResult:
    """per run (order): the order, whether its outputs equal the fresh ones, the comparison
    text, what each task changed in the process; the folder of the kept outputs"""

    results: str
    runs: list[dict[str, Any]]

    @property
    def ok(self) -> bool:
        return all(r["equal"] for r in self.runs)

    def text(self) -> str:
        lines = []
        for r in self.runs:
            state = "equal" if r["equal"] else "DIFFERENT"
            lines.append(f"{r['name']}: {state} to fresh processes ({r['files']} files)")
            if not r["equal"]:
                lines.append(r["report"])
        return "\n".join(lines)


def _spec(spec: str | Mapping[str, Any]) -> dict[str, Any]:
    if isinstance(spec, (str, os.PathLike)):
        with open(spec, encoding="utf-8") as fh:
            return json.load(fh)
    return dict(spec)


def _python_env() -> dict[str, str]:
    env = dict(os.environ)
    import petakit5d

    root = os.path.dirname(os.path.dirname(os.path.abspath(petakit5d.__file__)))
    old = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join([root, old]) if old else root
    return env


def _child(args: Sequence[str], cwd: str, python: str | None = None) -> Any:
    # a process of the harness, started in the work folder
    import subprocess

    return subprocess.run(
        [python or sys.executable, "-u", "-m", "petakit5d.jobs.isolation", *args],
        env=_python_env(),
        capture_output=True,
        text=True,
        cwd=cwd,
    )


def run_harness(
    spec: str | Mapping[str, Any],
    work: str,
    results: str,
    orders: int | Sequence[Sequence[int]] = 5,
    seed: int = 0,
    jobs: int = 1,
    python: str | None = None,
    reset: bool = True,
) -> HarnessResult:
    """Runs the tasks of ``spec`` in fresh processes and in one process in ``orders`` random
    orders (seeded; or the given orders of 0-based task indices), keeps each run's
    ``<work>/out`` in ``<results>/<run>`` and compares each order with the fresh run. spec:
    ``{"tasks": [record, ...], "prepare": record or None}`` (a JSON file or a dict); the
    processes start in the work folder. ``jobs``: fresh processes at once; ``reset`` False:
    the one-process runs without the restore between tasks (to show what tasks leave)."""
    import random
    import shutil

    s = _spec(spec)
    n = len(s["tasks"])
    work = os.path.abspath(work)
    results = os.path.abspath(results)
    os.makedirs(results, exist_ok=True)
    spec_file = os.path.join(results, "spec.json")
    with open(spec_file, "w", encoding="utf-8") as fh:
        json.dump(s, fh, allow_nan=True)
    rng = random.Random(seed)
    plans: list[tuple[str, list[int] | None]] = [("fresh", None)]
    if isinstance(orders, int):
        for k in range(1, orders + 1):
            order = list(range(n))
            rng.shuffle(order)
            plans.append((f"order_{k}", order))
    else:
        plans += [(f"order_{k}", list(o)) for k, o in enumerate(orders, 1)]
    runs: list[dict[str, Any]] = []
    for name, order in plans:
        _new_work(work, s.get("prepare"), spec_file, python)
        t0 = time.perf_counter()
        if order is None:
            _fresh(spec_file, n, jobs, work, python)
            records: list[dict[str, Any]] = []
        else:
            rec_file = os.path.join(results, f"{name}.ndjson")
            args = ["sequence", spec_file, ",".join(map(str, order)), rec_file]
            r = _child(args + ([] if reset else ["noreset"]), work, python)
            if r.returncode != 0:
                raise RuntimeError(f"isolation harness: {name} failed:\n{r.stdout}\n{r.stderr}")
            with open(rec_file, encoding="utf-8") as fh:
                records = [json.loads(x) for x in fh if x.strip()]
        seconds = time.perf_counter() - t0
        dest = os.path.join(results, name)
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.move(os.path.join(work, "out"), dest)
        runs.append({"name": name, "order": order, "seconds": seconds, "tasks": records})
    shutil.rmtree(work, ignore_errors=True)
    out = HarnessResult(results, [])
    from petakit5d.testing.compare import compare_trees, scan_tree

    ref = scan_tree(os.path.join(results, "fresh"))
    for r in runs[1:]:
        act = scan_tree(os.path.join(results, r["name"]))
        rep = compare_trees(ref, act)
        r["equal"] = rep.ok
        r["files"] = len(ref.items)
        r["report"] = "" if rep.ok else rep.text(r["name"])
        out.runs.append(r)
    with open(os.path.join(results, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"fresh_seconds": runs[0]["seconds"], "runs": out.runs}, fh, indent=1)
    return out


def _new_work(work: str, prepare: Any, spec_file: str, python: str | None) -> None:
    import shutil

    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(os.path.join(work, "out"))
    if prepare:
        r = _child(["prepare", spec_file], work, python)
        if r.returncode != 0:
            raise RuntimeError(f"isolation harness: prepare failed:\n{r.stdout}\n{r.stderr}")


def _fresh(spec_file: str, n: int, jobs: int, work: str, python: str | None) -> None:
    from concurrent.futures import ThreadPoolExecutor

    def one(i: int) -> None:
        r = _child(["task", spec_file, str(i)], work, python)
        if r.returncode != 0:
            raise RuntimeError(f"isolation harness: task {i + 1} failed:\n{r.stdout}\n{r.stderr}")

    with ThreadPoolExecutor(max(1, int(jobs))) as ex:
        list(ex.map(one, range(n)))


def _job_like_start() -> None:
    # as the job entry point starts: the numeric libraries' thread pools sized to our CPUs
    from petakit5d.jobs.exec import limit_library_threads

    limit_library_threads()


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    usage = (
        "usage: python -m petakit5d.jobs.isolation run <spec> <work> <results> [orders] [seed]"
        " [jobs]\n       ... task <spec> <i> | sequence <spec> <i,j,...> <records> [noreset]"
        " | prepare <spec>"
    )
    if not args:
        print(usage, file=sys.stderr)
        return 2
    cmd = args[0]
    if cmd == "run" and len(args) >= 4:
        nums = [int(x) for x in args[4:7]]
        res = run_harness(args[1], args[2], args[3], *nums)
        print(res.text())
        return 0 if res.ok else 1
    if cmd in ("task", "sequence", "prepare") and len(args) >= 2:
        _job_like_start()
        s = _spec(args[1])
        if cmd == "prepare":
            r = run_in_process(s["prepare"], None)
            return 0 if r.error is None else 1
        if cmd == "task" and len(args) == 3:
            r = run_in_process(s["tasks"][int(args[2])], None)
            return 0 if r.error is None else 1
        if cmd == "sequence" and len(args) in (4, 5):
            order = [int(x) for x in args[2].split(",") if x]
            baseline = None if args[4:] == ["noreset"] else ProcessState.capture()
            status = 0
            with open(args[3], "w", encoding="utf-8") as fh:
                for i in order:
                    r = run_in_process(s["tasks"][i], baseline)
                    status = status or (r.error is not None)
                    line = {"task": i, "seconds": r.seconds, "changed": r.changed}
                    line["error"] = None if r.error is None else repr(r.error)
                    fh.write(json.dumps(line) + "\n")
                    fh.flush()
            return int(status)
    print(usage, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
