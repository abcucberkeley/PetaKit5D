"""Pools of single jobs: the single jobs this process submits with the same function,
PythonLaunchStr, CPUs, SlurmParam, SlurmConstraint, time limit and job log folder run each
other's tasks (MATLAB's task_runner_run 'single_*' actions and task_runner.sh's steal mode 3).

A pool is a folder of single jobs (``jobs.single.TaskFolder``),
``<job log folder>/task_runner_single_<session>_<n>/``, made with the second task of a kind
(the first runs as before). Each job is a batch of one task, its own; tasks are added as their
jobs are submitted. Besides the task lines the pool has:

    runnable.txt  the tasks' numbers, each written after its line (a job reads only those)
    jobs.txt      '<task> <job>_<array task>' per submitted task (written after sbatch)
    claims/       the claims of jobs.claims: t<task> (the claimer), f<task> (it failed there),
                  r<task> (taken back by this process: it runs or submits the task itself)

Job side (``run_job``, the job entry point's runner of a pool): the job claims its own task
(one symbolic link, no listing; if another process has it, the job ends at once) and runs it in
its process, then the pool's unstarted tasks in submission order, each only while Slurm keeps
its job pending with a reason (resources, priority, a limit, a hold; not reason None, which
may start at once, nor a task whose job is not recorded yet), and cancels that job while it is
pending (``scancel --state=PENDING``). A task runs once in a runner: a failed one is marked
failed and the job takes no more. The process state is given back after each task, and the
runner's guards (memory, open files, threads, 10 minutes) go on in a new process, as
``jobs.runner`` does; tasks of other jobs only while the job has the time for them. The job
leaves the pool when it ends.

Master side (``jobs.single`` calls these): ``task``, ``add``, ``job``, ``live``, ``take``,
``status``, ``end``. A task is taken back (``take``) before it is submitted again or run by the
caller: its claim becomes this process's (a claim of an owner that ended is replaced in one
rename), so no runner takes it. A claim's owner may still run its task while its Slurm job is
in the queue (or, without a job, its process runs on this node); a failed claim is not live; a
query without an answer counts as live. Only this process removes the claims of owners that
ended, so no task runs in two processes at once.
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from petakit5d.jobs import claims as _claims
from petakit5d.jobs import single
from petakit5d.jobs.slurm import GONE, PENDING, RUNNING, job_key
from petakit5d.jobs.writers import node_name, own_job

__all__ = [
    "CLAIM_CACHE_SECONDS",
    "JOBS_FILE",
    "POOL_RUNNER",
    "RUNNABLE_FILE",
    "SLURM_TIMEOUT",
    "Pool",
    "PoolTask",
    "add",
    "end",
    "job",
    "live",
    "pending",
    "reset",
    "run_job",
    "status",
    "take",
    "task",
]

POOL_RUNNER = "petakit5d.jobs.pools:run_job"
RUNNABLE_FILE = "runnable.txt"
JOBS_FILE = "jobs.txt"
# the claims of a pool and their owners' jobs, read at most this often per pool while its
# tasks' jobs are pending (s)
CLAIM_CACHE_SECONDS = 10.0
# the job's squeue and scancel calls (MATLAB's timeout 10)
SLURM_TIMEOUT = 10.0

_LINUX = sys.platform.startswith("linux")


# --- the master side -------------------------------------------------------------------------


@dataclass
class PoolTask:
    """a task of a pool: its kind, number, state ('new', 'submitted' or 'taken': the master
    took it back) and the key of its job"""

    kind: str
    i: int
    state: str = "new"
    job_key: str = ""
    folder: str = ""

    def command(self, launch: str) -> str:
        return single.job_command(self.folder, self.i, launch)


@dataclass
class Pool:
    """a pool of this process: its folder, whether it takes tasks (False after it couldn't
    be made or written), and the claims read while its tasks' jobs are pending"""

    folder: single.TaskFolder | None
    on: bool = True
    claim_cache: dict[str, Any] | None = None
    claim_time: float = 0.0

    @property
    def claims(self) -> _claims.Claims:
        return _claims.Claims(self.folder.path)


@dataclass
class _Session:
    pools: dict[str, Pool] = field(default_factory=dict)  # kind -> pool
    tasks: dict[str, PoolTask] = field(default_factory=dict)  # task key -> task
    jobs: dict[str, str] = field(default_factory=dict)  # job key -> task key


_SESSION = _Session()


def reset() -> None:
    """forget the pools of this process (tests; jobs.single.reset calls it)"""
    global _SESSION
    _SESSION = _Session()


def _slurm(timeout: float | None = None) -> Any:
    return single._slurm(timeout)


def task(key: str) -> PoolTask | None:
    """the pool task of a task (its key in jobs.single), None if it has none or this process
    took it back"""
    t = _SESSION.tasks.get(key)
    return None if t is None or t.state == "taken" else t


def _make_pool(parent: str, cpus: int) -> single.TaskFolder:
    """a new pool folder: the claims folder (a probe that it takes symbolic links), the task
    files, the mark of this process"""
    f = single.TaskFolder.make(
        parent, single.POOL_PREFIX, POOL_RUNNER, {"kind": "pool", "cpus": int(cpus)}
    )
    store = _claims.Claims(f.path)
    try:
        if not store.usable():
            raise OSError(f"the claims folder {store.folder} takes no symbolic links")
        for name in (RUNNABLE_FILE, JOBS_FILE):
            with open(os.path.join(f.path, name), "ab"):
                pass
    except OSError:
        # (this process leaves it, so that it goes)
        f.leave_master(set())
        raise
    return f


def add(key: str, kind: str, parent: str, line: Mapping[str, Any], cpus: int) -> PoolTask | None:
    """adds the task (its line) to the pool of its kind, which its first task makes in parent;
    None if it can't be added (it then runs as before, with a warning once per pool)"""
    if key in _SESSION.tasks:
        return task(key)
    p = _SESSION.pools.get(kind)
    if p is not None and not p.on:
        return None
    try:
        if p is None:
            p = Pool(_make_pool(parent, cpus))
            _SESSION.pools[kind] = p
        i = p.folder.add(line)
        # (its number last: a job takes only tasks in runnable.txt, whose lines are there)
        with open(os.path.join(p.folder.path, RUNNABLE_FILE), "a", encoding="ascii") as fh:
            fh.write(f"{i}\n")
    except OSError as err:
        import warnings

        warnings.warn(
            single.SingleJobWarning(
                f"The single job runs without task runners: {err}", "task_runner_run:single"
            ),
            stacklevel=3,
        )
        # (the pool's files may be incomplete: its later tasks run as before)
        if p is None:
            p = Pool(None, on=False)
        p.on = False
        _SESSION.pools[kind] = p
        return None
    t = PoolTask(kind, i, folder=p.folder.path)
    _SESSION.tasks[key] = t
    return t


def job(key: str, job_id: int | None, array_id: int | None) -> None:
    """records the job of a pool task (jobs.txt: a runner that takes the task cancels it while
    it is pending); without one (its submission failed) the task is taken back"""
    if job_id is None:
        take(key)
        return
    t = _SESSION.tasks[key]
    t.state = "submitted"
    t.job_key = job_key(job_id, array_id)
    _SESSION.jobs[t.job_key] = key
    p = _SESSION.pools[t.kind]
    p.folder.jobs[t.i] = t.job_key
    try:
        with open(os.path.join(p.folder.path, JOBS_FILE), "a", encoding="ascii") as fh:
            fh.write(f"{t.i} {t.job_key}\n")
    except OSError as err:
        # (its job then starts and ends at once if a runner took its task)
        import warnings

        warnings.warn(
            single.SingleJobWarning(
                f"The job of a single job's task is not recorded: {err}", "task_runner_run:single"
            ),
            stacklevel=3,
        )


def _read_claim(store: _claims.Claims, i: int) -> tuple[str, bool]:
    """the owner of task i's claim ('' if none) and whether it is marked failed; opening the
    claims folder first revalidates it on NFS, so that claims other nodes just made are seen"""
    try:
        with os.scandir(store.folder):
            pass
    except OSError:
        pass
    owner = store.owner_of(i - 1) or ""
    failed = os.path.islink(os.path.join(store.folder, f"f{i}"))
    return owner, failed


def _my_job() -> tuple[int, int]:
    job_id, arr = own_job()
    return job_id, int(arr) if arr.isdigit() else 0


def _owner_live(owner: str, failed: bool, own_job_key: str) -> bool:
    """whether the owner of a claim may still run its task: its Slurm job is in the queue
    (squeue gave no answer: it may be) or, without a job or in this process's job, its process
    runs on this node; not a claim marked failed, of the task's own job (the caller has its
    state) or of this process"""
    if not owner or failed:
        return False
    p = _claims.parse_owner(owner)
    if p is None:
        return True  # (an owner the runners didn't write: taken as live)
    job_id, arr, node, pid, _ = p
    if f"{job_id}_{arr}" == own_job_key:
        return False
    if job_id > 0 and (job_id, arr) != _my_job():
        st = _slurm().batch_status([job_id], [arr])[0]
        return st is None or st >= PENDING
    return node == node_name() and pid != os.getpid() and _LINUX and os.path.isdir(f"/proc/{pid}")


def live(key: str) -> tuple[bool, bool]:
    """whether a process other than the task's own job has its claim and may still run it, and
    whether the task is claimed"""
    t = _SESSION.tasks[key]
    p = _SESSION.pools[t.kind]
    owner, failed = _read_claim(p.claims, t.i)
    return _owner_live(owner, failed, t.job_key), bool(owner)


def take(key: str) -> bool:
    """this process takes the task back: it replaces a claim whose owner ended (and its failed
    mark) by its own and marks the task taken (r<task>), so that no runner takes it; False if
    another process has it"""
    t = _SESSION.tasks[key]
    if t.state == "taken":
        return True
    p = _SESSION.pools[t.kind]
    store = p.claims
    owner, failed = _read_claim(store, t.i)
    if _owner_live(owner, failed, t.job_key):
        return False
    me = _claims.owner(0)
    c = store.folder
    ok = True
    try:
        try:
            os.remove(os.path.join(c, f"f{t.i}"))
        except FileNotFoundError:
            pass
        if not owner:
            # (fails if a runner just claimed it)
            os.symlink(me, os.path.join(c, f"t{t.i}"))
        else:
            # a claim whose owner ended is replaced at once (only its owner, which ended, and
            # this process remove claims), so that no runner claims the task in between
            tmp = os.path.join(c, f".take_t{t.i}_{os.getpid()}")
            try:
                os.remove(tmp)
            except FileNotFoundError:
                pass
            os.symlink(me, tmp)
            os.replace(tmp, os.path.join(c, f"t{t.i}"))
        with open(os.path.join(c, f"r{t.i}"), "w"):
            pass
    except OSError:
        ok = False
    if not ok:
        # another process just claimed it, or (if no live process has it now) the claims
        # folder failed us, and then no runner can claim it either
        owner, failed = _read_claim(store, t.i)
        ok = not _owner_live(owner, failed, t.job_key)
    if ok:
        t.state = "taken"
    return ok


def status(job_id: int, array_id: int, st: int) -> int:
    """the state of a pool task's own job whose Slurm state is st (check_single_job_status):
    after the job ended, 1 while another process has the task, else -1 (the task is taken
    back); while the job is pending, -1 if a runner of another job failed the task (it is taken
    back and the job cancelled), and the claim of a runner that ended without finishing it is
    removed, so that the job runs the task when it starts"""
    if st == RUNNING:
        return st
    key = job_key(job_id, array_id)
    tk = _SESSION.jobs.get(key)
    if tk is None:
        return st
    t = _SESSION.tasks[tk]
    if t.state == "taken" or t.job_key != key:
        return st
    if st == PENDING:
        if pending(tk) == "failed" and take(tk):
            _slurm().cancel([key])
            return GONE
        return st
    if not take(tk):
        return RUNNING
    return st


def pending(key: str) -> str:
    """while the task's own job is pending: 'live' if a runner of another job has the task,
    'failed' if one failed it, 'released' if one had it and ended without finishing it (its
    claim is removed, so that the job runs the task when it starts, or another runner takes
    it), else 'none'. The pool's claims and their owners' jobs are read at most every 10 s."""
    t = _SESSION.tasks[key]
    p = _SESSION.pools[t.kind]
    now = time.monotonic()
    if p.claim_cache is None or now - p.claim_time > CLAIM_CACHE_SECONDS:
        p.claim_cache = {"names": p.claims.names(), "jobs": {}}
        p.claim_time = now
    cache = p.claim_cache
    names = cache["names"]
    if f"t{t.i}" not in names:
        return "none"
    if f"f{t.i}" in names:
        return "failed"
    owner = p.claims.owner_of(t.i - 1) or ""
    q = _claims.parse_owner(owner)
    if q is None:
        alive = True
    elif q[0] > 0 and (q[0], q[1]) != _my_job():
        k = f"{q[0]}_{q[1]}"
        if k not in cache["jobs"]:
            st = _slurm().batch_status([q[0]], [q[1]])[0]
            cache["jobs"][k] = RUNNING if st is None else st
        alive = cache["jobs"][k] >= PENDING
    else:
        alive = q[2] == node_name() and q[3] != os.getpid() and os.path.isdir(f"/proc/{q[3]}")
    if alive:
        return "live"
    try:
        os.remove(os.path.join(p.claims.folder, f"t{t.i}"))
    except OSError:
        pass
    names.discard(f"t{t.i}")
    return "released"


def end(mark: int) -> list[single.TaskFolder]:
    """the end of a call: the pools made after mark take no more tasks (a later task of their
    kind makes a new pool) and are forgotten; returns their folders, which jobs.single.end
    leaves (a pool that couldn't be made stays off for the session, as in MATLAB)"""
    kinds = [
        kind for kind, p in _SESSION.pools.items() if p.folder is not None and p.folder.index > mark
    ]
    folders = [_SESSION.pools.pop(kind).folder for kind in kinds]
    gone = set(kinds)
    for tk in [k for k, t in _SESSION.tasks.items() if t.kind in gone]:
        t = _SESSION.tasks.pop(tk)
        if t.job_key:
            _SESSION.jobs.pop(t.job_key, None)
    return folders


# --- the job side ----------------------------------------------------------------------------


@dataclass
class _JobState:
    """what a pool job's next process needs (PETAKIT5D_RUNNER_STATE, as jobs.runner's)"""

    t0: float
    est: float = 0.0  # longest task so far (s)
    steal: bool = True
    status: int = 0
    ran: int = 0
    processes: int = 1
    own_done: bool = False
    next: int = -1  # a task claimed for the next process (1-based), -1 none

    @classmethod
    def take(cls) -> _JobState:
        from petakit5d.jobs.runner import STATE_VARIABLE

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


def _read_numbers(path: str) -> list[int]:
    try:
        with open(path, encoding="ascii", errors="replace") as fh:
            return [int(x) for x in fh.read().split() if x.isdigit()]
    except OSError:
        return []


def _read_jobs(folder: str) -> dict[int, str]:
    """task -> job key of jobs.txt (the last line of a task counts)"""
    out: dict[int, str] = {}
    try:
        with open(os.path.join(folder, JOBS_FILE), encoding="ascii", errors="replace") as fh:
            for row in fh:
                f = row.split()
                if len(f) == 2 and f[0].isdigit():
                    out[int(f[0])] = f[1]
    except OSError:
        pass
    return out


class _Runner:
    """one process of a pool job (see the module docstring)"""

    def __init__(self, folder: str, own: int, state: _JobState) -> None:
        from petakit5d.jobs.tasks import CallFolder

        self.folder = folder
        self.own = own
        self.state = state
        self.store = _claims.Claims(folder)
        self.me = _claims.owner(1)
        job_id, arr = own_job()
        self.name = f"{job_id}_{arr}" if arr else str(job_id)
        meta = CallFolder(folder).meta.get("single", {})
        cpus = os.environ.get("SLURM_CPUS_PER_TASK", "")
        # (a job with other CPUs than its pool's, at another place, takes only its own task)
        self.steal = state.steal and not (
            cpus.isdigit() and int(cpus) != int(meta.get("cpus", cpus))
        )
        from petakit5d.jobs import runner

        self.margin = runner._env_float(runner.MEMORY_VARIABLE)
        self.margin = runner.MEMORY_MARGIN_GB if self.margin is None else self.margin
        minutes = runner._env_float(runner.MINUTES_VARIABLE)
        self.minutes = runner.MAX_MINUTES if minutes is None else minutes

    def _waiting(self) -> set[int]:
        """the tasks whose jobs Slurm keeps pending with a reason (not None, which may start at
        once), of the jobs recorded in jobs.txt"""
        jobs = _read_jobs(self.folder)
        if not jobs:
            return set()
        reasons = _slurm(SLURM_TIMEOUT).pending_reasons()
        if not reasons:
            return set()
        keep = {k for k, r in reasons.items() if r != "None"}
        return {i for i, k in jobs.items() if k in keep}

    def _pick(self) -> int | None:
        """the next unstarted task of the pool's other jobs, in submission order, claimed (its
        pending job is cancelled); None if there is none or the job has no time for it"""
        if self.state.next > 0:
            i, self.state.next = self.state.next, -1
            return i
        from petakit5d.jobs.runner import time_left_ok

        if not self.steal:
            return None
        names = self.store.names()
        todo = [
            i
            for i in _read_numbers(os.path.join(self.folder, RUNNABLE_FILE))
            if i != self.own and f"t{i}" not in names and f"r{i}" not in names
        ]
        if not todo or not time_left_ok(self.state.est):
            return None
        waiting: set[int] | None = None
        for i in todo:
            if waiting is None:
                waiting = self._waiting()
            if i not in waiting:
                continue
            if self.store.claim(i - 1, self.me):
                self._cancel_job_of(i)
                return i
        return None

    def _cancel_job_of(self, i: int) -> None:
        # the job of a task this job took: cancelled while it is pending (it would only start
        # to end at once, at the front of the queue); if it just started, its own claim failed
        # and it ends anyway
        k = _read_jobs(self.folder).get(i)
        if k:
            _slurm(SLURM_TIMEOUT).run(["scancel", "--state=PENDING", k], SLURM_TIMEOUT)

    def run(self) -> int:
        from petakit5d.jobs import runner
        from petakit5d.jobs.isolation import ProcessState, run_in_process

        st = self.state
        baseline = ProcessState.capture()
        start_rss = runner._process_numbers()[0]
        t_proc = time.monotonic()
        warm: tuple[float, int, int] | None = None
        ran_here = 0
        if not st.own_done:
            # its own task first, claimed without listing the pool; a job whose task another
            # process has ends at once
            st.own_done = True
            if not self.store.claim(self.own - 1, self.me):
                print(f"Task runner: no task left to run for {self.own}", flush=True)
                return st.status
            nxt: int | None = self.own
        else:
            nxt = self._pick()
        while nxt is not None:
            i = nxt
            line = single.read_task_line(self.folder, i)
            where = "" if i == self.own else " of another job"
            print(f"Task runner {self.name}: task {i}{where}: {_describe(line)}", flush=True)
            res = self._run_task(line, baseline, run_in_process)
            ran_here += 1
            st.ran += 1
            st.est = max(st.est, res.seconds)
            rss, fds, thr = runner._process_numbers()
            if res.changed:
                print(f"Task runner {self.name}: task {i} changed {', '.join(res.changed)}")
            if res.error is not None:
                print(
                    f"Task runner {self.name}: task {i} failed after {res.seconds:.3f} s; it is "
                    "marked failed, and this job takes no more tasks."
                )
                self.store.mark_failed(i - 1, self.me)
                st.status = 1
                break
            print(
                f"Task runner {self.name}: task {i} done in {res.seconds:.3f} s, memory "
                f"{rss / runner._GB_KB:.2f} GB",
                flush=True,
            )
            if warm is None:
                warm = (rss, fds, thr)
                limit = start_rss + 2 * self.margin * runner._GB_KB
            else:
                limit = warm[0] + self.margin * runner._GB_KB
            why = ""
            if rss > limit:
                rss = runner._process_numbers(trim=True)[0]
            if rss > limit:
                gb, gb0 = rss / runner._GB_KB, start_rss / runner._GB_KB
                why = f"memory {gb:.2f} GB ({gb0:.2f} GB at its start)"
            elif ran_here > 1 and fds > warm[1] + runner.MORE_FILES:
                why = f"{fds} open files, {warm[1]} after its first task"
            elif ran_here > 1 and thr > warm[2] + runner.MORE_THREADS:
                why = f"{thr} threads, {warm[2]} after its first task"
            elif time.monotonic() - t_proc >= self.minutes * 60:
                why = f"tasks for {(time.monotonic() - t_proc) / 60:.1f} min"
            nxt = self._pick()
            if why and nxt is not None:
                self._go_on(why, nxt)
        if st.ran == 0:
            print(f"Task runner {self.name}: no task left to run for {self.own}")
        return st.status

    def _run_task(self, line: Mapping[str, Any], baseline: Any, run_in_process: Any) -> Any:
        metrics = None
        if line.get("metrics"):
            try:
                from petakit5d.jobs import metrics
            except ImportError:
                metrics = None
        call = single._OneTask(line)

        def before() -> Any:
            return metrics.task_start(call, [call.index], "slurm_python") if metrics else None

        def after(tm: Any, err: BaseException | None) -> None:
            if metrics:
                metrics.task_stop(tm, err)

        return run_in_process(line["task"], baseline, before, after)

    def _go_on(self, why: str, nxt: int) -> None:
        """the job's next task (claimed) in a new process (exec); returns only if the new
        process can't be started (then this one goes on with it)"""
        from petakit5d.jobs import runner

        st = self.state
        st.next = nxt
        print(f"Task runner {self.name}: {why}, a new process runs task {nxt}", flush=True)
        try:
            runner._reexec(self.folder, self.own, st)
        except OSError as err:
            print(f"Task runner {self.name}: no new process ({err}); this one goes on")
        st.next = -1


def _describe(line: Mapping[str, Any]) -> str:
    from petakit5d.jobs.tasks import describe

    return f"(task {line.get('task_id')}) {describe(line['task'])}"


def run_job(folder: str, i: int) -> int:
    """The job entry point's runner of a pool (call.json "runner"): the job of pool task i
    (see the module docstring); returns the exit status (1 if a task failed)."""
    state = _JobState.take()
    t0 = state.t0
    try:
        status = _Runner(folder, int(i), state).run()
        print(f"Elapsed time is {time.time() - t0:f} seconds.", flush=True)
        return status
    finally:
        single.leave(folder, i)
