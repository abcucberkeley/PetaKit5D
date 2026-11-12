"""The job framework's local process pool: MATLAB's matlab_pool_generic_computing_wrapper
(``localWorkerNum > 1`` without cluster) and matlab_parfor_generic_computing_wrapper
(``clusterType 'parfor'``).

Tasks run in worker processes, never in threads of this interpreter. A worker is a fresh
interpreter (``sys.executable``) started for the call, with this process's ``sys.path``,
current folder and environment, and runs task records as a cluster job does (by import
path, with their JSON arguments), one at a time. Each worker gets a share of the CPUs, as
MATLAB gives each pool worker ``floor(maxNumCompThreads / workers)`` threads: the thread
variables of OpenMP and the BLAS libraries are set to it and, on Linux, the worker is
bound to its own block of this process's CPUs, so ``core.threads.available_cpus()`` in a
task gives that share. A task's output (Python's and the C libraries') is printed with
the task when it ends, as MATLAB prints a future's diary. A worker restores its current
folder, environment and warning filters after each task.

The pool ends with a normal shutdown: each worker's task stream is closed, the worker
finishes its task and exits, and the pool waits for it. Nothing is killed: when another
run makes the final output, tasks that are running finish (MATLAB cancels them).

Seams of run_tasks (framework.py): ``run_pool(call, worker_num=...)`` runs first with
``local_worker_num > 1`` and leaves what it didn't finish to the sequential loop;
``run_parfor(call)`` is ``cluster_type='parfor'``.
"""

from __future__ import annotations

import collections
import json
import math
import os
import queue
import subprocess
import sys
import threading
import time
import traceback
from collections.abc import Callable, Sequence
from typing import IO, Any, NamedTuple

from petakit5d.jobs.files import batch_file_exist, exist
from petakit5d.jobs.tasks import Call, task_record

# MATLAB: a worker takes about 2 GB besides its task (the pool may use half the memory)
WORKER_MEMORY_GB = 2

# environment of a worker: its 1-based id, its threads, its CPUs (Linux), its GPU
ENV_WORKER = "PETAKIT5D_POOL_WORKER"
ENV_THREADS = "PETAKIT5D_POOL_THREADS"
ENV_CPUS = "PETAKIT5D_POOL_CPUS"
ENV_GPU = "PETAKIT5D_POOL_GPU"
# thread counts of the libraries a task may load (read when they load)
THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
# how often the pool's wait wakes up (Ctrl-C is seen then)
_POLL = 0.5

__all__ = [
    "PoolWorker",
    "pool_worker_num",
    "run_on_pool",
    "run_parfor",
    "run_pool",
    "worker_info",
]


class PoolWorker(NamedTuple):
    """The pool worker this process is: 1-based id, threads, GPU (CUDA_VISIBLE_DEVICES)."""

    id: int
    threads: int
    gpu: str | None


def worker_info() -> PoolWorker | None:
    """The pool worker this process (or its parent task process) is, None outside a pool
    (MATLAB's getCurrentTask)."""
    wid = os.environ.get(ENV_WORKER, "")
    if not wid.isdigit():
        return None
    threads = os.environ.get(ENV_THREADS, "")
    return PoolWorker(int(wid), int(threads) if threads.isdigit() else 1, os.environ.get(ENV_GPU))


def pool_worker_num(
    worker_num: float, mem_allocate: float | None, memory: float | None = None
) -> float:
    """The workers of the local pool, as MATLAB's framework sizes it: ``worker_num``, at most
    ``floor(memory / 2 / (mem_allocate + 2))`` when the task memory (GB) is known (half of
    the memory, a worker taking 2 GB besides its task). ``memory`` is getSystemMemory's
    (GB; NaN when unknown, which sets no limit, as MATLAB's min ignores NaN)."""
    if mem_allocate is None:
        return worker_num
    if memory is None:
        memory = _system_memory()
    denom = mem_allocate + WORKER_MEMORY_GB
    if math.isnan(memory) or math.isnan(denom):
        return worker_num
    cap = memory / 2 / denom if denom else math.copysign(math.inf, memory)
    return min(worker_num, math.floor(cap) if math.isfinite(cap) else cap)


def run_pool(call: Call, worker_num: float | None = None) -> list[bool] | None:
    """Runs the call's tasks on a local pool of ``worker_num`` workers (the call's
    ``local_worker_num`` by default), limited by memory (``pool_worker_num``), before the
    sequential loop. Returns the pool's done flags, or None if no pool ran (one worker or
    fewer fit). Failed tasks are left unfinished, for the sequential loop to retry."""
    s = call.settings
    if worker_num is None:
        worker_num = s.get("local_worker_num", 1)
    workers = pool_worker_num(worker_num, s.get("mem_allocate"))
    if not workers > 1:
        return None
    return run_on_pool(
        call.inputs,
        call.outputs,
        call.records,
        final_out_fullpath=s.get("final_out_fullpath", ""),
        worker_num=workers,
        language=call.language,
    )


def run_on_pool(
    input_paths: Sequence[str | Sequence[str]],
    output_paths: Sequence[str],
    tasks: Sequence[Any],
    *,
    final_out_fullpath: str = "",
    worker_num: float = 2,
    language: str = "python",
) -> list[bool]:
    """MATLAB's matlab_pool_generic_computing_wrapper: runs the tasks at most
    ``worker_num`` at a time on a pool of worker processes, in task order, and returns
    whether each task's output exists. Tasks whose output exists or whose input is missing
    are skipped (not reported done), and none starts once the final output exists. A task
    that raises, or whose worker ends, is reported and left unfinished. If no worker can
    be started, no task runs."""
    records = [task_record(r, language) for r in tasks]
    outputs = list(output_paths)
    n = len(outputs)
    if final_out_fullpath and batch_file_exist([final_out_fullpath])[0]:
        return [True] * n
    done = list(batch_file_exist(outputs, None, True))
    if all(done):
        return done
    todo = [i for i, ok in enumerate(_input_exist(input_paths)) if ok and not done[i]]
    if not todo:
        return done

    workers = _worker_count(worker_num, len(todo))
    pool = _Pool(_slots(workers), language)
    try:
        if not pool.start():
            return done
        print(f"Run {len(todo)} tasks on the local pool ({workers} at a time)...")
        progress = _Progress(outputs, done, final_out_fullpath)
        pool.run(todo, records, progress.task_done)
        if progress.final_exists:
            done[:] = [True] * n
    finally:
        pool.close()
    return done


def run_parfor(call: Call) -> list[bool]:
    """``cluster_type='parfor'``: the tasks run on a local pool only (never in this
    process), in up to ``max_trial_num`` rounds (MATLAB's 3 loops), and the done flags are
    returned as run_tasks gives them (a task with a missing input counts as done).

    Workers: for a GPU job one per usable GPU, worker k on GPU ``(k - 1) mod n`` of the
    GPUs this process sees (its CUDA_VISIBLE_DEVICES set to that one); without a usable
    GPU, or for CPU jobs, ``local_worker_num`` if above 1, else the CPUs this process may
    use, at most ``max_job_num``, within the memory limit of the pool (at least 1).
    (MATLAB's framework passes this wrapper parameters it doesn't take, so the 'parfor'
    path fails there; this ports what it was for.)"""
    s = call.settings
    n = call.n
    final_out = s.get("final_out_fullpath", "")
    if final_out and batch_file_exist([final_out])[0]:
        return [True] * n
    gpus: list[str] = []
    if s.get("gpu_job"):
        gpus = _parfor_gpus(_usable_gpu_num())
        if not gpus:
            print(
                "Warning: no usable GPU for the GPU job; the pool workers get no GPU of their own.",
                file=sys.stderr,
            )
    if gpus:
        workers = len(gpus)
    else:
        workers = s.get("local_worker_num", 1)
        if not workers > 1:
            workers = min(s.get("max_job_num", math.inf), _available_cpus())
        workers = pool_worker_num(workers, s.get("mem_allocate"))

    in_exist = _input_exist(call.inputs)
    out_exist = list(batch_file_exist(call.outputs, None, True))
    todo = [i for i in range(n) if in_exist[i] and not out_exist[i]]
    if todo:
        pool = _Pool(_slots(_worker_count(workers, len(todo)), gpus), call.language)
        try:
            if pool.start():
                for g, slot in enumerate(pool.slots, 1):
                    if slot.gpu is not None:
                        print(f"Pool worker {g} uses GPU {slot.gpu}")
                progress = _Progress(call.outputs, out_exist, final_out)
                for _ in range(int(s.get("max_trial_num", 3))):
                    print(
                        f"Run {len(todo)} tasks on the local pool ({len(pool.slots)} at a time)..."
                    )
                    pool.run(todo, call.records, progress.task_done)
                    if progress.final_exists:
                        return [True] * n
                    out_exist = list(batch_file_exist(call.outputs, None, True))
                    todo = [i for i in todo if not out_exist[i]]
                    if not todo or not pool.alive():
                        break
        finally:
            pool.close()
    return [out_exist[i] or not in_exist[i] for i in range(n)]


# --- the pool -------------------------------------------------------------------------------


class _Slot(NamedTuple):
    id: int  # 1-based
    threads: int
    cpus: tuple[int, ...] | None  # CPUs to bind to (Linux), None: no binding
    gpu: str | None  # CUDA_VISIBLE_DEVICES of the worker


def _slots(workers: int, gpus: Sequence[str] = ()) -> list[_Slot]:
    # floor(CPUs / workers) threads each (at least 1); on Linux disjoint blocks of the CPUs
    avail = _available_cpus()
    threads = max(1, avail // workers)
    cpus = _affinity()
    bind = cpus is not None and len(cpus) >= threads * workers
    out = []
    for k in range(workers):
        block = tuple(cpus[k * threads : (k + 1) * threads]) if bind and cpus else None
        gpu = gpus[k % len(gpus)] if gpus else None
        out.append(_Slot(k + 1, threads, block, gpu))
    return out


class _Worker:
    """A worker process and the thread that reads its messages."""

    def __init__(self, slot: _Slot, proc: subprocess.Popen, events: queue.Queue) -> None:
        self.slot = slot
        self.proc = proc
        self.ready = False
        self.alive = True
        self.task: int | None = None
        self.reader = threading.Thread(target=_read_messages, args=(self, events), daemon=True)
        self.reader.start()

    def send(self, msg: dict[str, Any]) -> bool:
        try:
            assert self.proc.stdin is not None
            self.proc.stdin.write((json.dumps(msg, allow_nan=True) + "\n").encode())
            self.proc.stdin.flush()
            return True
        except (OSError, ValueError):
            return False

    def end_stream(self) -> None:
        try:
            if self.proc.stdin is not None:
                self.proc.stdin.close()
        except OSError:
            pass


def _read_messages(w: _Worker, events: queue.Queue) -> None:
    assert w.proc.stdout is not None
    for line in w.proc.stdout:
        try:
            msg = json.loads(line)
        except ValueError:
            continue
        events.put((w, msg))
    events.put((w, None))


class _Pool:
    """Worker processes of one call: start, run tasks (several rounds), close."""

    def __init__(self, slots: list[_Slot], language: str) -> None:
        self.slots = slots
        self.language = language
        self.events: queue.Queue = queue.Queue()
        self.workers: list[_Worker] = []
        self.restarts = 0

    def start(self) -> bool:
        for slot in self.slots:
            w = self._spawn(slot)
            if w is None:
                break
            self.workers.append(w)
        if not self.workers:
            return False
        # wait until a worker can take tasks, or all ended
        while not any(w.ready for w in self.workers):
            if not self.alive():
                print("The local pool could not start its workers, no task runs on it.")
                return False
            self._handle(self._next_event(), None, None)
        return True

    def _spawn(self, slot: _Slot) -> _Worker | None:
        try:
            proc = subprocess.Popen(
                [_python_executable(), "-c", _WORKER_CODE],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                env=_worker_env(slot),
            )
        except OSError as err:
            print(f"The local pool can't start a worker: {err}")
            return None
        w = _Worker(slot, proc, self.events)
        w.send({"path": list(sys.path)})
        return w

    def alive(self) -> bool:
        return any(w.alive for w in self.workers)

    def _next_event(self) -> tuple[_Worker, dict[str, Any] | None]:
        while True:
            try:
                return self.events.get(timeout=_POLL)
            except queue.Empty:
                continue

    def run(
        self,
        todo: Sequence[int],
        records: Sequence[dict[str, Any]],
        task_done: Callable[[int, float, str | None, str], bool],
    ) -> None:
        """Runs the tasks ``todo`` (0-based, in order) until all ended or ``task_done``
        returns True (stop: no more tasks start); tasks no worker can take stay undone."""
        pending = collections.deque(todo)
        state = {"running": 0, "stop": False}

        def dispatch(w: _Worker) -> None:
            while pending and not state["stop"] and w.alive and w.task is None:
                i = pending.popleft()
                if w.send({"i": i, "language": self.language, "task": records[i]}):
                    w.task = i
                    state["running"] += 1
                else:
                    pending.appendleft(i)
                    return

        for w in self.workers:
            if w.ready:
                dispatch(w)
        # (after a stop, tasks that are running finish while the pool closes)
        while not state["stop"] and (pending or state["running"]) and self.alive():
            w, msg = self._next_event()
            if self._handle((w, msg), state, task_done) and pending:
                self._restart(w)
            for other in self.workers:
                if other.ready:
                    dispatch(other)
        if pending and not state["stop"]:
            print(f"The local pool has no workers left: {len(pending)} tasks did not run on it.")

    def _handle(self, event, state, task_done) -> bool:
        # one message of a worker; True if the worker ended
        w, msg = event
        if msg is not None and "ready" in msg:
            w.ready = True
            return False
        if msg is None:
            w.alive = False
            status = _exit_status(w.proc)
            i = w.task
            w.task = None
            if i is not None and state is not None:
                state["running"] -= 1
                text = f"the pool worker {w.slot.id} ended (exit status {status})"
                if task_done(i, math.nan, text, ""):
                    state["stop"] = True
            return True
        i = w.task
        w.task = None
        if i is None or state is None:
            return False
        state["running"] -= 1
        if task_done(i, msg.get("t", math.nan), msg.get("error"), msg.get("output", "")):
            state["stop"] = True
        return False

    def _restart(self, w: _Worker) -> None:
        # a worker that ended while tasks wait is replaced, at most once per slot on average
        if self.restarts >= len(self.slots):
            return
        self.restarts += 1
        new = self._spawn(w.slot)
        if new is not None:
            print(f"Pool worker {w.slot.id} ended; started a new one.")
            self.workers.append(new)

    def close(self) -> None:
        """Normal shutdown: ends each worker's task stream and waits for it to exit."""
        for w in self.workers:
            w.end_stream()
        for w in self.workers:
            w.proc.wait()
            w.reader.join()
            if w.proc.stdout is not None:
                w.proc.stdout.close()
            w.alive = False


def _exit_status(proc: subprocess.Popen) -> str:
    try:
        return str(proc.wait(timeout=10))
    except subprocess.TimeoutExpired:
        return "unknown"


class _Progress:
    """What the MATLAB pool prints after each task, and its done flags."""

    def __init__(self, outputs: Sequence[str], done: list[bool], final_out: str) -> None:
        self.outputs = outputs
        self.done = done
        self.final_out = final_out
        self.n_done = sum(done)
        self.final_exists = False
        self.ts = time.monotonic()

    def task_done(self, i: int, t: float, error: str | None, output: str) -> bool:
        if output:
            # (a text stream on Windows wrote \r\n; print writes \n as the platform does)
            output = output.replace("\r\n", "\n")
            print(output, end="" if output.endswith("\n") else "\n")
        if error:
            print(f"Task {i + 1} failed: {error}")
        self.done[i] = batch_file_exist([self.outputs[i]], None, True)[0]
        print(f"Done! Elapsed time is {_fmt_f(t)} seconds.")
        n = len(self.done)
        if self.n_done < sum(self.done):
            self.n_done = sum(self.done)
            print(
                f"Time {time.monotonic() - self.ts:0.2f} s: {self.n_done} / {n} "
                f"({self.n_done / n * 100:0.2f}%) are finished!"
            )
        # another run made the final output: stop, as the framework does
        if self.final_out and batch_file_exist([self.final_out])[0]:
            self.final_exists = True
            print(
                f"Time {time.monotonic() - self.ts:0.2f} s: The final output file "
                f"{self.final_out} already exists!"
            )
            return True
        return False


# --- the worker process -------------------------------------------------------------------

_WORKER_CODE = "from petakit5d.jobs.local import _worker_main; _worker_main()"


def _worker_env(slot: _Slot) -> dict[str, str]:
    env = dict(os.environ)
    root = _package_root()
    if root:
        old = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = os.pathsep.join([root, old]) if old else root
    for name in THREAD_VARIABLES:
        env[name] = str(slot.threads)
    env[ENV_WORKER] = str(slot.id)
    env[ENV_THREADS] = str(slot.threads)
    env.pop(ENV_CPUS, None)
    env.pop(ENV_GPU, None)
    if slot.cpus:
        env[ENV_CPUS] = ",".join(str(c) for c in slot.cpus)
    if slot.gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = slot.gpu
        env[ENV_GPU] = slot.gpu
    return env


def _worker_main() -> None:
    """A pool worker: reads tasks (JSON lines) from stdin, runs each, writes its result
    (JSON line) to the original stdout; exits when stdin ends."""
    cpus = os.environ.get(ENV_CPUS, "")
    if cpus and hasattr(os, "sched_setaffinity"):
        try:
            os.sched_setaffinity(0, [int(c) for c in cpus.split(",")])
        except (OSError, ValueError):
            pass
    # the protocol keeps its own copies of stdin and stdout; tasks read nothing, and what
    # is printed between tasks goes to stderr
    proto_in = os.fdopen(os.dup(0), "rb")
    proto_out = os.fdopen(os.dup(1), "wb")
    null = os.open(os.devnull, os.O_RDONLY)
    os.dup2(null, 0)
    os.close(null)
    os.dup2(2, 1)

    def send(msg: dict[str, Any]) -> None:
        proto_out.write((json.dumps(msg, allow_nan=True) + "\n").encode())
        proto_out.flush()

    hello = proto_in.readline()
    if not hello:
        return
    sys.path[:] = json.loads(hello)["path"]
    send({"ready": os.getpid()})
    import tempfile

    with tempfile.TemporaryFile() as cap:
        for line in proto_in:
            msg = json.loads(line)
            send(_run_in_worker(msg, cap))


def _run_in_worker(msg: dict[str, Any], cap: IO[bytes]) -> dict[str, Any]:
    import warnings

    from petakit5d.jobs.tasks import run_record

    cwd = os.getcwd()
    environ = dict(os.environ)
    filters = list(warnings.filters)
    cap.seek(0)
    cap.truncate()
    _flush_all()
    saved = (os.dup(1), os.dup(2))
    os.dup2(cap.fileno(), 1)
    os.dup2(cap.fileno(), 2)
    err = None
    t0 = time.perf_counter()
    try:
        if msg["language"] == "bash":
            from petakit5d.jobs.framework import run_shell

            run_shell(msg["task"]["bash"])  # its exit status is ignored, as in MATLAB
        else:
            run_record(msg["task"])
    except (Exception, SystemExit):
        err = traceback.format_exc()
    finally:
        t = time.perf_counter() - t0
        _flush_all()
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(saved[0])
        os.close(saved[1])
    # the next task starts from the same state
    if os.getcwd() != cwd:
        os.chdir(cwd)
    if dict(os.environ) != environ:
        os.environ.clear()
        os.environ.update(environ)
    warnings.filters[:] = filters
    cap.seek(0)
    output = cap.read().decode("utf-8", "replace")
    return {"i": msg["i"], "t": t, "error": err, "output": output}


def _flush_all() -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.flush()
        except (AttributeError, OSError, ValueError):
            pass
    if os.name == "posix":
        # C stdio buffers of extension modules
        try:
            import ctypes

            ctypes.CDLL(None).fflush(None)
        except (OSError, AttributeError):
            pass


# --- helpers --------------------------------------------------------------------------------


def _input_exist(inputs: Sequence[str | Sequence[str]]) -> list[bool]:
    # MATLAB: batch_file_exist for paths; every path of each entry otherwise
    if all(isinstance(p, str) for p in inputs):
        return list(batch_file_exist(list(inputs), None, True))  # type: ignore[arg-type]
    return [
        exist(p) if isinstance(p, str) else all(exist(q) for q in p)  # type: ignore[union-attr]
        for p in inputs
    ]


def _parfor_gpus(n_gpu: int) -> list[str]:
    # the GPUs this process sees, as CUDA_VISIBLE_DEVICES values, the first n_gpu of them
    if n_gpu <= 0:
        return []
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None and visible.strip():
        ids = [x.strip() for x in visible.split(",") if x.strip()]
    else:
        ids = [str(k) for k in range(n_gpu)]
    return ids[:n_gpu]


def _usable_gpu_num() -> int:
    """MATLAB's usable_gpu_num (gpu.device)."""
    from petakit5d.gpu.device import usable_gpu_count

    return usable_gpu_count()[0]


def _system_memory() -> float:
    from petakit5d.core import memory

    return memory.get_system_memory()


def _available_cpus() -> int:
    from petakit5d.core import threads

    return threads.available_cpus()


def _affinity() -> list[int] | None:
    if not hasattr(os, "sched_getaffinity"):
        return None
    try:
        return sorted(os.sched_getaffinity(0))
    except OSError:
        return None


def _python_executable() -> str:
    return sys.executable


def _package_root() -> str:
    import petakit5d

    return os.path.dirname(os.path.dirname(os.path.abspath(petakit5d.__file__)))


def _worker_count(worker_num: float, n_tasks: int) -> int:
    # MATLAB: min(max(1, round(workerNum)), nT); round(NaN) is NaN, and max(1, NaN) is 1
    if math.isnan(worker_num):
        return 1
    if math.isinf(worker_num):
        return n_tasks if worker_num > 0 else 1
    t = math.trunc(worker_num)
    r = int(t + (math.copysign(1, worker_num) if abs(worker_num - t) >= 0.5 else 0))
    return min(max(1, r), n_tasks)


def _fmt_f(t: float) -> str:
    # MATLAB's %f (NaN, Inf)
    if math.isnan(t):
        return "NaN"
    if math.isinf(t):
        return "Inf" if t > 0 else "-Inf"
    return f"{t:f}"
