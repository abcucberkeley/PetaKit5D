"""Task runners: a job's process runs the tasks of its batch one after another, then unstarted
tasks of the call's other batches (MATLAB's task_runner, task_runner_run and task_runner.sh, for
Python's jobs). On by default in cluster mode (run_tasks' ``task_runner``; Linux only, as in
MATLAB); bash tasks and jobs with in-job slots have none.

Job side. call.json names ``petakit5d.jobs.runner:run`` as its runner, so the job entry point
(``python -m petakit5d.jobs.exec <call folder> <batch>``) hands the batch to ``run``. The process
claims a task (``jobs.claims``), runs it in the process, gives the process state back after it
(``jobs.isolation``) and claims the next one: first the tasks of its own batch in order, then
(stealing) unclaimed tasks of the call's other batches, from the last task backwards, of
batches with the key of its own (the CPUs and GPU of their jobs), never of the batch the
master runs itself (key 'master'), of a batch the master submitted again, or a task whose
output existed or whose input was missing when the call started (it stays in its batch).
Tasks of other batches are taken only while the job has the time for them: twice its longest
task and a minute before its end (SLURM_JOB_END_TIME; before its first task, while less than
half of its time is used). A claimed task whose output (or the call's final output) exists is
skipped. A batch the master submitted again runs as without runners (its tasks in order; the
first error ends it), as do calls without runner settings.

The job goes on in a new process (``os.execv`` of the same command, with the runner's state
in PETAKIT5D_RUNNER_STATE, which tasks never see) before its next task when: a task failed (it
is marked failed, and the job then takes no tasks of other batches); the process holds more
memory than a warm process (after its first task up to twice PETAKIT5D_RUNNER_MEM_GB, 0.5 GB,
more than at its start, later up to the margin more than after its first task); it has 64
more open files or 256 more threads than after its first task; or it has run tasks for
PETAKIT5D_RUNNER_MINUTES (10) minutes.

Master side. ``cluster_runners(call, shape)`` is the object the cluster loop (``jobs.cluster``)
asks: the runner settings for call.json (``call_extra``); before it submits a batch again,
whether to leave it (``keep``: while live runners of other jobs hold its tasks, or, for a
first submission, all of them; claims of owners that ended are removed, and a batch submitted
again is marked so); which tasks of its own batch the master runs (``master_tasks``: those it
can claim). While jobs' runners take tasks of other batches, it also keeps the call's jobs in
the queue to maxJobNum (a job that runs other batches' tasks counts, though its own batch is
done) and cancels the pending jobs of batches whose tasks are all done or taken by live
runners (``scancel --state=PENDING``), with one squeue per round of the loop.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from collections import deque
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from typing import Any

from petakit5d.jobs import claims as _claims

__all__ = [
    "MAX_MINUTES",
    "MEMORY_MARGIN_GB",
    "MORE_FILES",
    "MORE_THREADS",
    "PROTOCOL_VERSION",
    "RUNNER",
    "STATE_VARIABLE",
    "ClusterRunners",
    "NoRunners",
    "batch_keys",
    "cluster_runners",
    "run",
    "run_plain",
    "time_left_ok",
]

RUNNER = "petakit5d.jobs.runner:run"
PROTOCOL_VERSION = 1

# the runner's state for the next process of the job (removed before any task sees it)
STATE_VARIABLE = "PETAKIT5D_RUNNER_STATE"
# memory a process may hold above a warm one (GB), and how long it runs tasks (minutes)
MEMORY_VARIABLE = "PETAKIT5D_RUNNER_MEM_GB"
MINUTES_VARIABLE = "PETAKIT5D_RUNNER_MINUTES"
MEMORY_MARGIN_GB = 0.5
MAX_MINUTES = 10.0
# open files and threads a process may have above those after its first task
MORE_FILES = 64
MORE_THREADS = 256
# tasks of other batches only while the job has 2 x its longest task + 60 s left
TIME_LEFT_FACTOR = 2.0
TIME_LEFT_MARGIN = 60.0
# failed claims of other batches' tasks in a row before the claims are listed again
RELIST_MISSES = 32

_GB_KB = 1048576.0


def batch_keys(n_batches: int, cpus: int, gpu: bool, master_batch: int = 0) -> list[str]:
    """per batch (1-based ``master_batch``) the key of the processes that run its tasks: their
    CPUs and GPU (``'4'``, ``'8gpu'``), ``'master'`` for the batch the master runs itself"""
    key = f"{int(cpus)}{'gpu' if gpu else ''}"
    keys = [key] * int(n_batches)
    if 1 <= master_batch <= n_batches:
        keys[master_batch - 1] = "master"
    return keys


# --- the job side ------------------------------------------------------------------------------


@dataclass
class _JobState:
    """what a job's next process needs of its runner's state"""

    t0: float  # time.time() when the job's first process started
    est: float = 0.0  # longest task so far (s)
    steal: bool = True  # (False after a failed task)
    status: int = 0  # 1 once a task failed
    ran: int = 0  # tasks run
    processes: int = 1
    stealer: bool = False  # the job's stealer mark is made
    next: int = -1  # a task claimed for the next process (0-based), -1 none
    next_kind: str = ""

    @classmethod
    def take(cls) -> _JobState:
        """the state the previous process left (removed from the environment), else a new one"""
        text = os.environ.pop(STATE_VARIABLE, "")
        if text:
            try:
                d = json.loads(text)
                st = cls(**{k: d[k] for k in d if k in cls.__dataclass_fields__})
                st.processes += 1
                return st
            except (ValueError, TypeError):
                pass
        return cls(time.time())


def run(call_folder: str, batch: int) -> int:
    """The job entry point's runner (call.json ``"runner"``): runs batch ``batch`` (1-based)
    and other batches' tasks; returns the exit status (1 if a task failed, 2 for a batch that
    isn't there)."""
    from petakit5d.jobs.tasks import CallFolder

    state = _JobState.take()
    cf = CallFolder(call_folder)
    own = list(cf.batch(batch))
    if batch < 1 or not own:
        nb = -(-cf.n // int(cf.meta["task_batch_num"]))
        print(f"petakit5d.jobs.runner: batch {batch} is not in 1..{nb}", file=sys.stderr)
        return 2
    cfg = cf.meta.get("task_runner") or {}
    store = _claims.Claims(call_folder)
    if cfg.get("version") != PROTOCOL_VERSION or cf.language != "python":
        return run_plain(cf, own)
    if not store.usable():
        print(
            f"Task runner {_name()}: the claims folder {store.folder} takes no symbolic links; "
            "the batch runs without task runners."
        )
        return run_plain(cf, own)
    if store.resubmitted(batch):
        if state.processes > 1:
            # (the master submitted the batch again while this job ran: its job runs it)
            print(f"Task runner {_name()}: batch {batch} was submitted again; this job stops.")
            _elapsed(state.t0)
            return state.status
        return run_plain(cf, own)
    return _Runner(cf, batch, own, cfg, store, state).run()


def run_plain(cf: Any, indices: Sequence[int]) -> int:
    """the tasks at ``indices`` in this process as the job entry point runs a batch without
    runners: in order, the first exception ends the batch (exit status 1); a task metrics
    record each (mode slurm_python)"""
    import traceback

    from petakit5d.jobs.tasks import run_record

    metrics = _metrics(cf)
    t0 = time.perf_counter()
    status = 0
    for i, x in zip(indices, cf.lines(indices), strict=True):
        tm = metrics.task_start(cf, [i], "slurm_python") if metrics else None
        try:
            run_record(x["task"])
        except Exception as err:
            if metrics:
                metrics.task_stop(tm, err)
            sys.stdout.flush()
            print(f"Error in task {x['index']} ({x['task']['func']}):", file=sys.stderr)
            traceback.print_exc()
            status = 1
            break
        except BaseException as err:
            if metrics:
                metrics.task_stop(tm, err)
            raise
        if metrics:
            metrics.task_stop(tm, None)
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)
    return status


def time_left_ok(est: float, now: float | None = None) -> bool:
    """whether the job has the time for a task of about ``est`` s: 2 x est + 60 s before its end
    (est 0: while less than half of its time is used); True outside Slurm or without an end"""
    end = _env_float("SLURM_JOB_END_TIME")
    if end is None:
        return True
    now = time.time() if now is None else now
    if est > 0:
        need = TIME_LEFT_FACTOR * est + TIME_LEFT_MARGIN
    else:
        start = _env_float("SLURM_JOB_START_TIME")
        need = (end - (now if start is None else start)) / 2
    return end - now >= need


class _Picker:
    """the next task of a job's runner: its own batch's in order, then other batches' (last
    task first) that the job's key allows; each claimed before it is returned"""

    def __init__(
        self,
        n: int,
        bs: int,
        batch: int,
        keys: Sequence[str],
        skip: set[int],
        store: _claims.Claims,
        me: str,
        start: Sequence[int],
    ) -> None:
        self.bs = bs
        self.store = store
        self.me = me
        self.own: deque[int] = deque(start)
        mykey = keys[batch - 1] if 0 < batch <= len(keys) else ""
        self.others = [
            i
            for i in range(n - 1, -1, -1)
            if i // bs + 1 != batch
            and i not in skip
            and i // bs < len(keys)
            and keys[i // bs] == mykey
        ]
        self.queue: deque[int] = deque()
        self.known: set[int] = set()  # tasks seen taken (claimed, or of resubmitted batches)

    def next_own(self) -> int | None:
        while self.own:
            i = self.own.popleft()
            if self.store.claim(i, self.me):
                return i
        return None

    def next_other(self) -> int | None:
        """the next unclaimed task of another batch, claimed; None if there is none. Each task
        is tried once: a claim that fails (another process has it) is not tried again"""
        misses = 0
        listed = False
        while True:
            if not self.queue:
                if listed:
                    return None
                names = self.store.names()
                self.known |= _claims.claimed_set(names)
                again = _claims.resubmitted_set(names)
                self.queue = deque(
                    i for i in self.others if i not in self.known and i // self.bs + 1 not in again
                )
                listed = True
                misses = 0
                continue
            i = self.queue.popleft()
            self.known.add(i)
            if self.store.claim(i, self.me):
                return i
            misses += 1
            if misses >= RELIST_MISSES:
                # (others took these since the listing: list again rather than try each)
                self.queue.clear()
                listed = False


class _Runner:
    """one process of a job's runner (see the module docstring)"""

    def __init__(
        self,
        cf: Any,
        batch: int,
        own: list[int],
        cfg: dict[str, Any],
        store: _claims.Claims,
        state: _JobState,
    ) -> None:
        self.cf = cf
        self.batch = batch
        self.own = own
        self.store = store
        self.state = state
        self.bs = int(cf.meta["task_batch_num"])
        self.keys = list(cfg.get("keys") or [])
        self.me = _claims.owner(1)
        self.name = _name()
        key = self.keys[batch - 1] if batch <= len(self.keys) else ""
        self.steal = state.steal and key not in ("", "master") and _cpus_match(key)
        self.picker = _Picker(
            cf.n, self.bs, batch, self.keys, set(cfg.get("skip") or ()), store, self.me, own
        )
        self.final_out = cf.meta.get("final_out_fullpath") or ""
        self.margin = _env_float(MEMORY_VARIABLE)
        self.margin = MEMORY_MARGIN_GB if self.margin is None else self.margin
        minutes = _env_float(MINUTES_VARIABLE)
        self.minutes = MAX_MINUTES if minutes is None else minutes

    def _next(self) -> tuple[int, str] | None:
        """the next task, claimed: ('o' own, 'a' another batch's)"""
        st = self.state
        if st.next >= 0:
            i, kind = st.next, st.next_kind or "o"
            st.next, st.next_kind = -1, ""
            return i, kind
        i = self.picker.next_own()
        if i is not None:
            return i, "o"
        if not self.steal or not time_left_ok(st.est):
            return None
        while True:
            i = self.picker.next_other()
            if i is None:
                return None
            if self.store.resubmitted(i // self.bs + 1):
                self.store.give_back(i, self.me)
                continue
            return i, "a"

    def run(self) -> int:
        from petakit5d.jobs.files import batch_file_exist
        from petakit5d.jobs.isolation import ProcessState, run_in_process

        st = self.state
        metrics = _metrics(self.cf)
        baseline = ProcessState.capture()
        start_rss = _process_numbers()[0]
        t_proc = time.monotonic()
        warm: tuple[float, int, int] | None = None
        ran_here = 0
        while True:
            nxt = self._next()
            if nxt is None:
                break
            i, kind = nxt
            line = self.cf.lines([i])[0]
            out = line["output"]
            paths = [out, self.final_out] if self.final_out else [out]
            seen = batch_file_exist(paths)
            if self.final_out and seen[-1]:
                print(f"Task runner {self.name}: the final output {self.final_out} exists.")
                break
            if seen[0]:
                print(f'Task runner {self.name}: output "{out}" already exists')
                continue
            if kind == "a" and not st.stealer:
                key = _claims.job_key_of(self.me)
                if key:
                    self.store.mark_stealer(key, self.me)
                st.stealer = True
            where = "" if kind == "o" else f" of batch {i // self.bs + 1}"
            print(f"Task runner {self.name}: task {i + 1}{where}", flush=True)

            def before(i: int = i) -> Any:
                return metrics.task_start(self.cf, [i], "slurm_python") if metrics else None

            def after(tm: Any, err: BaseException | None) -> None:
                if metrics:
                    metrics.task_stop(tm, err)

            res = run_in_process(line["task"], baseline, before, after)
            ran_here += 1
            st.ran += 1
            st.est = max(st.est, res.seconds)
            rss, fds, thr = _process_numbers()
            if res.changed:
                print(f"Task runner {self.name}: task {i + 1} changed {', '.join(res.changed)}")
            if res.error is not None:
                print(
                    f"Task runner {self.name}: task {i + 1} failed after {res.seconds:.3f} s; "
                    "it is marked failed, and this job takes no tasks of other batches."
                )
                self.store.mark_failed(i, self.me)
                st.status, st.steal = 1, False
                self.steal = False
                self._go_on("a task failed")
                continue
            print(
                f"Task runner {self.name}: task {i + 1} done in {res.seconds:.3f} s, memory "
                f"{rss / _GB_KB:.2f} GB",
                flush=True,
            )
            if warm is None:
                warm = (rss, fds, thr)
                limit = start_rss + 2 * self.margin * _GB_KB
            else:
                limit = warm[0] + self.margin * _GB_KB
            why = ""
            if rss > limit:
                rss = _process_numbers(trim=True)[0]
            if rss > limit:
                why = f"memory {rss / _GB_KB:.2f} GB ({start_rss / _GB_KB:.2f} GB at its start)"
            elif ran_here > 1 and fds > warm[1] + MORE_FILES:
                why = f"{fds} open files, {warm[1]} after its first task"
            elif ran_here > 1 and thr > warm[2] + MORE_THREADS:
                why = f"{thr} threads, {warm[2]} after its first task"
            elif time.monotonic() - t_proc >= self.minutes * 60:
                why = f"tasks for {(time.monotonic() - t_proc) / 60:.1f} min"
            if why:
                self._go_on(why)
        if st.ran == 0:
            print(f"Task runner {self.name}: no task left to run for batch {self.batch}")
        _elapsed(st.t0)
        return st.status

    def _go_on(self, why: str) -> None:
        """the job's next task in a new process (exec), if there is one; returns only if the
        new process can't be started (then this one goes on)"""
        nxt = self._next()
        if nxt is None:
            return
        st = self.state
        st.next, st.next_kind = nxt
        print(f"Task runner {self.name}: {why}, a new process runs task {nxt[0] + 1}", flush=True)
        try:
            _reexec(self.cf.folder, self.batch, st)
        except OSError as err:
            print(f"Task runner {self.name}: no new process ({err}); this one goes on")
            st.next, st.next_kind = nxt


def _reexec(folder: str, batch: int, st: _JobState) -> None:
    """replaces this process by a new one of the same job entry point (os.execv; the pid, and
    so the claims' owner, stays)"""
    import atexit

    argv = list(getattr(sys, "orig_argv", None) or [])
    k = argv.index("-m") if "-m" in argv else -1
    if not (0 <= k < len(argv) - 1 and argv[k + 1] == "petakit5d.jobs.exec"):
        argv = [sys.executable, "-u", "-m", "petakit5d.jobs.exec", folder, str(batch)]
    if not sys.executable or not hasattr(os, "execv") or not os.access(sys.executable, os.X_OK):
        raise OSError("this Python can't start itself again")
    os.environ[STATE_VARIABLE] = json.dumps(asdict(st))
    try:
        atexit._run_exitfuncs()
    except Exception:
        pass
    for s in (sys.stdout, sys.stderr):
        try:
            s.flush()
        except (AttributeError, OSError, ValueError):
            pass
    try:
        os.execv(sys.executable, argv)
    finally:
        os.environ.pop(STATE_VARIABLE, None)


def _process_numbers(trim: bool = False) -> tuple[float, int, int]:
    """(resident memory in kB, open files, threads) of this process (Linux; zeros elsewhere);
    with ``trim``, after the C heap gave back the free memory it can (malloc_trim)"""
    if not sys.platform.startswith("linux"):
        return 0.0, 0, 0
    if trim:
        _malloc_trim()
    rss, thr = 0.0, 0
    try:
        with open("/proc/self/status", encoding="ascii", errors="replace") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    rss = float(line.split()[1])
                elif line.startswith("Threads:"):
                    thr = int(line.split()[1])
    except OSError:
        pass
    try:
        fds = len(os.listdir("/proc/self/fd"))
    except OSError:
        fds = 0
    return rss, fds, thr


_libc: list[Any] = []


def _malloc_trim() -> None:
    if not _libc:
        try:
            import ctypes

            _libc.append(ctypes.CDLL("libc.so.6").malloc_trim)
        except (OSError, AttributeError):
            _libc.append(None)
    if _libc[0] is not None:
        try:
            _libc[0](0)
        except Exception:
            pass


def _cpus_match(key: str) -> bool:
    # a job with other CPUs than its batch's key (another place) takes no other batches' tasks
    digits = ""
    for c in key:
        if not c.isdigit():
            break
        digits += c
    have = os.environ.get("SLURM_CPUS_PER_TASK", "")
    return not (digits and have.isdigit() and int(have) != int(digits))


def _metrics(cf: Any) -> Any:
    folder = cf.meta.get("settings", {}).get("task_metrics")
    if not isinstance(folder, str) or not folder:
        return None
    try:
        from petakit5d.jobs import metrics
    except ImportError:
        print("Warning: task metrics are not available (petakit5d.jobs.metrics)", file=sys.stderr)
        return None
    return metrics


def _name() -> str:
    from petakit5d.jobs.writers import own_job

    job, arr = own_job()
    return f"{job}_{arr}" if arr else str(job)


def _elapsed(t0: float) -> None:
    print(f"Elapsed time is {time.time() - t0:f} seconds.", flush=True)


def _env_float(name: str) -> float | None:
    try:
        x = float(os.environ.get(name, ""))
    except ValueError:
        return None
    return x if math.isfinite(x) else None


# --- the master side ---------------------------------------------------------------------------


class NoRunners:
    """the cluster loop's runner object without runners (bash tasks, in-job slots)"""

    def call_extra(self) -> dict[str, Any]:
        return {}

    def waiting(self, flags: Sequence[bool], timestamp: float) -> None:
        pass

    def keep(self, indices: Sequence[int], retry: bool) -> bool:
        return False

    def master_tasks(self, indices: Sequence[int]) -> list[int]:
        return list(indices)

    def close(self) -> None:
        pass


def cluster_runners(call: Any, shape: Any) -> Any:
    """the master side of the task runners for the cluster loop (jobs.cluster): a
    ``ClusterRunners``, or ``NoRunners`` for bash tasks and jobs with in-job slots"""
    if call.language != "python" or getattr(shape, "slot_mode", False):
        return NoRunners()
    return ClusterRunners(call, shape)


class ClusterRunners:
    """The master side (see the module docstring). Claims are read once per round of the loop
    (``waiting`` starts a round). Slurm is asked only when a claim's owner is another job (one
    ``squeue -j`` of the owners' jobs not asked yet in the round), and while jobs' runners take
    tasks of other batches (one ``squeue --me`` per round: the call's jobs, by their
    comments)."""

    def __init__(self, call: Any, shape: Any, slurm: Any = None) -> None:
        from petakit5d.jobs import framework
        from petakit5d.jobs.files import batch_file_exist

        self.call = call
        s = call.settings
        n = call.n
        self.bs = int(shape.task_batch_num)
        nb = -(-n // self.bs)
        out = batch_file_exist(call.outputs, None, True)
        inp = framework.input_exist(call, out)
        self.skip = [i for i in range(n) if out[i] or not inp[i]]
        master = 0
        if s.get("master_compute", True):
            undone = [i for i in range(n) if not out[i]]
            if undone:
                master = undone[-1] // self.bs + 1
        self.keys = batch_keys(nb, shape.cpus_per_task, bool(s.get("gpu_job", False)), master)
        self.max_job = s.get("max_job_num", math.inf)
        self.me = _claims.owner(0)
        self._slurm = slurm
        self._store: _claims.Claims | None = None
        self._off = False
        self._names: set[str] | None = None
        self._view: dict[str, tuple[str, str]] | None | bool = False  # False: not asked
        self._live: dict[str, bool] = {}  # (this round) job key: running or pending
        self._approved = 0
        self.cancelled: list[str] = []

    # -- the protocol of the cluster loop ------------------------------------------------------

    def call_extra(self) -> dict[str, Any]:
        return {
            "runner": RUNNER,
            "task_runner": {"version": PROTOCOL_VERSION, "keys": self.keys, "skip": self.skip},
        }

    def waiting(self, flags: Sequence[bool], timestamp: float) -> None:
        """a new round of the loop (after its query): claims and Slurm are read again; pending
        jobs of batches that runners did are cancelled"""
        self._names = None
        self._view = False
        self._live = {}
        self._approved = 0
        if self.store is None:
            return
        if any(n.startswith("s") for n in self.names()):
            self._cancel_covered()

    def keep(self, indices: Sequence[int], retry: bool) -> bool:
        """before batch of ``indices`` (its tasks not done) is submitted again (``retry``) or
        for the first time in a later round: True to leave it this round"""
        if not indices or self.store is None:
            return False
        names = self.names()
        claimed, live = self.store.states(list(indices), names, self._listed)
        if any(c and lv is None for c, lv in zip(claimed, live, strict=True)):
            return True  # (Slurm didn't answer: nothing is released or submitted on it)
        alive = [bool(x) for x in live]
        if any(alive) and (retry or all(alive)):
            return True
        if self._at_job_limit():
            return True
        if retry:
            # (before the release: a runner that sees the claims gone sees the mark too)
            self.store.mark_resubmitted(indices[0] // self.bs + 1)
        dead = [i for i, c, a in zip(indices, claimed, alive, strict=True) if c and not a]
        if dead:
            self.store.release(dead)
            self._names = None
        self._approved += 1
        return False

    def master_tasks(self, indices: Sequence[int]) -> list[int]:
        """the tasks of the master's batch it runs itself: those without output that it can
        claim (claims of owners that ended are removed first)"""
        from petakit5d.jobs.files import batch_file_exist

        if self.store is None:
            return list(indices)
        found = batch_file_exist([self.call.outputs[i] for i in indices], None, True)
        todo = [i for i, f in zip(indices, found, strict=True) if not f]
        if not todo:
            return []
        self._names = None  # (the claims as they are now)
        names = self.names()
        claimed, live = self.store.states(todo, names, self._listed)
        dead = [i for i, c, lv in zip(todo, claimed, live, strict=True) if c and lv is False]
        if dead:
            self.store.release(dead)
            self._names = None
        out = []
        for i, c, lv in zip(todo, claimed, live, strict=True):
            if c and lv is not False:
                continue
            if self.store.claim(i, self.me):
                out.append(i)
        return out

    def close(self) -> None:
        pass

    # -- helpers ---------------------------------------------------------------------------------

    @property
    def store(self) -> _claims.Claims | None:
        """the call's claims (None before the call folder exists, or if it takes no symbolic
        links: then the master does as without runners, and so do the jobs)"""
        if self._store is None and not self._off and getattr(self.call, "folder", None):
            store = _claims.Claims(self.call.folder)
            if store.usable():
                self._store = store
            else:
                self._off = True
                print(f"Task runner: {store.folder} takes no symbolic links; no task runners.")
        return self._store

    @property
    def slurm(self) -> Any:
        if self._slurm is None:
            from petakit5d.jobs.slurm import Slurm

            self._slurm = Slurm()
        return self._slurm

    def names(self) -> set[str]:
        if self._names is None:
            self._names = self.store.names() if self.store is not None else set()
        return self._names

    def view(self) -> dict[str, tuple[str, str]] | None:
        """{job key: (state, comment)} of this user's jobs (squeue --me), once per round; None
        if Slurm didn't answer"""
        if self._view is False:
            r = self.slurm.run(["squeue", "--me", "-h", "-r", "-o", "%i|%t|%k"])
            if not r.ok:
                self._view = None
            else:
                v: dict[str, tuple[str, str]] = {}
                for line in r.stdout.splitlines():
                    parts = line.strip().split("|", 2)
                    if len(parts) == 3 and parts[0]:
                        v[parts[0].strip()] = (parts[1].strip(), parts[2].strip())
                self._view = v
        return self._view

    def _listed(self, keys: set[str]) -> set[str] | None:
        """those of ``keys`` (owners' jobs) that Slurm has running or pending (a completing
        job has ended): one squeue -j of those not asked this round; None if it didn't answer"""
        new = sorted(k for k in keys if k not in self._live)
        if new:
            found = self._states_of(new)
            if found is None:
                return None
            for k in new:
                self._live[k] = found.get(k, "CG") != "CG"
        return {k for k in keys if self._live.get(k)}

    def _states_of(self, keys: Sequence[str]) -> dict[str, str] | None:
        # squeue -h -r -o '%i %t' -j <keys> (1000 per call): {key: state}, None unknown
        out: dict[str, str] = {}
        for b in range(0, len(keys), 1000):
            part = list(keys[b : b + 1000])
            r = self.slurm.run(["squeue", "-h", "-r", "-o", "%i %t", "-j", ",".join(part)])
            if "Invalid job id" in r.output:
                continue
            if not r.ok:
                return None
            for line in r.stdout.splitlines():
                f = line.split()
                if len(f) >= 2:
                    out[f[0]] = f[1]
        return out

    def _view_listed(self, keys: set[str]) -> set[str] | None:
        # as _listed, from this round's view of all the user's jobs
        v = self.view()
        if v is None:
            return None
        return {k for k in keys if v.get(k, ("CG", ""))[0] != "CG"}

    def _call_jobs(self) -> dict[str, tuple[str, int]] | None:
        """{job key: (state, batch)} of this call's jobs that Slurm has (by their comment
        ``<call folder name>:<batch>``)"""
        v = self.view()
        if v is None:
            return None
        base = os.path.basename(os.path.normpath(self.call.folder)) + ":"
        out = {}
        for k, (state, comment) in v.items():
            if comment.startswith(base) and comment[len(base) :].isdigit() and state != "CG":
                out[k] = (state, int(comment[len(base) :]))
        return out

    def _at_job_limit(self) -> bool:
        """while jobs' runners take other batches' tasks, the loop doesn't count those jobs
        (their batches are done): the call's jobs in the queue (and those submitted this round)
        must stay below maxJobNum"""
        if not math.isfinite(self.max_job) or not any(n.startswith("s") for n in self.names()):
            return False
        jobs = self._call_jobs()
        if jobs is None:
            return True
        return len(jobs) + self._approved >= self.max_job

    def _cancel_covered(self) -> None:
        """cancels the pending jobs of batches whose tasks are all done or held by live
        runners of other jobs"""
        from petakit5d.jobs.files import batch_file_exist

        jobs = self._call_jobs()
        if not jobs:
            return
        pending = {k: b for k, (state, b) in jobs.items() if state == "PD"}
        if not pending:
            return
        n = self.call.n
        tasks = {
            k: list(range((b - 1) * self.bs, min(b * self.bs, n)))
            for k, b in pending.items()
            if 1 <= b
        }
        flat = sorted({i for v in tasks.values() for i in v})
        seen = batch_file_exist([self.call.outputs[i] for i in flat], None, True)
        found = dict(zip(flat, seen, strict=True))
        open_tasks = [i for i in flat if not found[i]]
        names = self.names()
        claimed, live = self.store.states(open_tasks, names, self._view_listed)
        held = {i for i, c, lv in zip(open_tasks, claimed, live, strict=True) if c and lv}
        for k, idx in tasks.items():
            if all(found[i] or i in held for i in idx):
                r = self.slurm.run(["scancel", "--state=PENDING", k])
                if r.ok:
                    self.cancelled.append(k)
                    print(
                        f"Task runner: the tasks of batch {pending[k]} are done or taken by "
                        f"other jobs; its pending job {k} is cancelled."
                    )
