"""Single jobs: one task per Slurm job, submitted and looked at by the caller's own loop
(MATLAB's generic_single_job_submit_wrapper, check_single_job_status, single_job_pools and the
loops of the wrappers that call them).

``submit_single_job(task, job_id, task_id, **params) -> (job_id, job_status, submit_status)``
never waits: a task whose job is pending or running returns at once (job_status 0 or 1,
submit_status 0); a task without a job in the queue is submitted (submit_status 1 with the new
job id), up to max_trial_num sbatch tries in the call. With master_compute and last_file the
caller runs the task itself: a pending job of it is cancelled, and nothing is submitted.
``check_single_job_status(job_id, task_id)`` is the callers' look at a job of a last try.
States: 1 running, 0 pending, -1 gone; a Slurm query without an answer counts as running
(MATLAB's rule), so nothing is submitted again on it.

The task never enters the job's command line. Each task is a line of a folder of single jobs
in the job log folder (tmp_dir if given), written by this process, and the job runs
``<PythonLaunchStr> -u -m petakit5d.jobs.exec <folder> <line>``:

- the first task of a kind, tasks submitted again, bash tasks, and all tasks with task
  runners off or off Linux go to ``pk5d_single_<session>_<n>/`` and run as before: the job
  runs its task and ends;
- from the second task of a kind on (the same function, PythonLaunchStr, CPUs, SlurmParam,
  SlurmConstraint, time limit and job log folder), first submissions join a pool,
  ``task_runner_single_<session>_<n>/`` (``jobs.pools``): a job runs its own task, then
  unstarted tasks of the pool's other jobs while Slurm keeps those pending, and cancels them.

A folder's ``live/`` holds a mark of this process (``m``) and of each task (``t<line>``). A
job removes its mark when it ends; ``end(mark)`` (or ``SingleJobs`` at the end of a call)
leaves the folders made since ``mark()``: this process's mark and those of tasks whose jobs are
no longer in the queue go, and the last one to leave removes the folder. A process that exits
leaves its folders the same way.

``SingleJobs`` is the callers' side (MATLAB's wrappers' loops): per item a job id, a trial
counter and whether its job is active; one squeue per round for all jobs (``query``) instead of
one per job; the last try's status check; the wait for results between rounds; the pools of
the call ended with it.

Session state (MATLAB's persistent variables: the Slurm check of the job log folder, the
kinds submitted, the folders) is this process's; ``reset()`` forgets it (tests).
"""

from __future__ import annotations

import json
import math
import numbers
import os
import shlex
import shutil
import sys
import time
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from petakit5d.jobs.config import default_python_launch_str
from petakit5d.jobs.shapes import ARRAY_ID_MODULUS, TIME_FORMAT_SINGLE, time_str
from petakit5d.jobs.slurm import (
    GONE,
    PENDING,
    RUNNING,
    Slurm,
    check_slurm_cluster,
    job_key,
    split_options,
)
from petakit5d.jobs.tasks import CALL_FILE, CALL_FORMAT, CALL_VERSION, describe, task_record

__all__ = [
    "JOB_RUNNER",
    "LIVE_FOLDER",
    "POOL_PREFIX",
    "SINGLE_JOB_PARAMS",
    "TASK_FOLDER_PREFIX",
    "SingleJobWarning",
    "SingleJobs",
    "TaskFolder",
    "check_single_job_status",
    "end",
    "job_command",
    "job_states",
    "leave",
    "mark",
    "pools_used",
    "read_task_line",
    "reset",
    "run_job",
    "submit_single_job",
]

# the job entry point's runner for tasks that run as before (call.json "runner")
JOB_RUNNER = "petakit5d.jobs.single:run_job"
TASK_FOLDER_PREFIX = "pk5d_single_"
POOL_PREFIX = "task_runner_single_"
FOLDER_VERSION = 1
LIVE_FOLDER = "live"
TASKS_FILE = "tasks.ndjson"
INDEX_FILE = "tasks.idx"
INDEX_WIDTH = 12  # digits of a byte offset in tasks.idx (then a line break)

DEFAULT_JOB_LOG_DIR = "../job_logs"
DEFAULT_JOB_LOG_FNAME = "../job_logs/job_%A_%a.out"
DEFAULT_JOB_ERROR_FNAME = "../job_logs/job_%A_%a.err"
# the callers' wait between rounds: 0.3 s per item left, between 3 and 30 s (MATLAB's wrappers)
WAIT_PER_ITEM = 0.3
WAIT_MIN = 3.0
WAIT_MAX = 30.0

_LINUX = sys.platform.startswith("linux")


class SingleJobWarning(UserWarning):
    """A single job runs without what it would have had (its pool, its task file)."""

    def __init__(self, message: str, identifier: str) -> None:
        super().__init__(message)
        self.identifier = identifier


@dataclass(frozen=True)
class Param:
    """One parameter: MATLAB name, Python name, default, kind (bool, num, num_or_none, str,
    any), meaning. STAND-IN for the parameter registry of F03 (as jobs.framework's)."""

    matlab: str
    python: str
    default: Any
    kind: str
    doc: str


# generic_single_job_submit_wrapper's parameters after funcStr, job_id and task_id, in MATLAB's
# order with its names and defaults (site defaults are empty: they come from the config file)
SINGLE_JOB_PARAMS: tuple[Param, ...] = (
    Param("clusterType", "cluster_type", "slurm", "str", "'slurm' (the only one)"),
    Param("parseCluster", "parse_cluster", True, "bool", "accepted (the caller decides)"),
    Param(
        "masterCompute",
        "master_compute",
        False,
        "bool",
        "the caller runs tasks too (with last_file: this one)",
    ),
    Param(
        "lastFile",
        "last_file",
        False,
        "bool",
        "with master_compute: the caller runs this task itself (a pending job of it is "
        "cancelled; nothing is submitted)",
    ),
    Param("jobLogDir", "job_log_dir", DEFAULT_JOB_LOG_DIR, "str", "folder of the job logs"),
    Param("jobLogFname", "job_log_fname", DEFAULT_JOB_LOG_FNAME, "str", "sbatch -o"),
    Param("jobErrorFname", "job_error_fname", DEFAULT_JOB_ERROR_FNAME, "str", "sbatch -e"),
    Param(
        "tmpDir",
        "tmp_dir",
        "",
        "str",
        "folder of the task files (empty: the folder of the job logs)",
    ),
    Param("maxCPUNum", "max_cpu_num", 24, "num", "most CPUs of the job"),
    Param("minCPUNum", "min_cpu_num", 1, "num", "fewest CPUs of the job"),
    Param("cpusPerTask", "cpus_per_task", 1, "num", "CPUs of the job"),
    Param("memPerCPU", "mem_per_cpu", 20.9, "num", "GB of memory per CPU of the job"),
    Param(
        "memAllocate",
        "mem_allocate",
        None,
        "num_or_none",
        "GB of memory the task needs (None: unknown); sets the CPUs of its job",
    ),
    Param("uuid", "uuid", "", "str", "accepted (unused, as in MATLAB)"),
    Param("maxTrialNum", "max_trial_num", 3, "num", "sbatch tries in this call"),
    Param("unitWaitTime", "unit_wait_time", 30, "num", "accepted (unused, as in MATLAB)"),
    Param("maxJobNum", "max_job_num", math.inf, "num", "accepted (unused, as in MATLAB)"),
    Param("taskBatchNum", "task_batch_num", 1, "num", "accepted (unused, as in MATLAB)"),
    Param(
        "PythonLaunchStr",
        "python_launch_str",
        "",
        "str",
        "how the job starts Python (empty: the interpreter of this process)",
    ),
    Param("BashLaunchStr", "bash_launch_str", "", "str", "shell commands before a bash task"),
    Param(
        "SlurmParam",
        "slurm_param",
        "",
        "str",
        "sbatch options of the job (partition, QOS, memory per CPU); from the config",
    ),
    Param("SlurmConstraint", "slurm_constraint", "", "str", "sbatch constraint of the job"),
    Param("jobTimeLimit", "job_time_limit", 24, "num_or_none", "hours of the job (None: none)"),
    Param(
        "language",
        "language",
        "python",
        "str",
        "'python' (a task record; 'matlab' is taken for it) or 'bash' (a command)",
    ),
    Param("GNUparallel", "gnu_parallel", False, "bool", "accepted (unused, as in MATLAB)"),
    Param("paraJobNum", "para_job_num", 1, "num", "accepted (unused, as in MATLAB)"),
    Param(
        "ConfigFile",
        "config_file",
        "",
        "str",
        "cluster configuration file (JSON or .mat) whose values override these",
    ),
    Param(
        "taskMetrics",
        "task_metrics",
        False,
        "bool",
        "record the task's resources in <job log folder>/task_metrics/",
    ),
    Param(
        "taskRunner",
        "task_runner",
        True,
        "bool",
        "from the second task of a kind on, jobs share a pool and run each other's tasks "
        "(off if the caller or the config says so)",
    ),
    Param(
        "SlurmAlternatives",
        "slurm_alternatives",
        None,
        "any",
        "other places for the job (see jobs.alternatives)",
    ),
    Param(
        "gpuHeadroom",
        "gpu_headroom",
        None,
        "num_or_none",
        "with SlurmAlternatives: CPUs and GB per idle GPU that CPU jobs leave free",
    ),
    Param("nice", "nice", 0, "num", "lower the job's priority by this (sbatch --nice)"),
)

# MATLAB-only parameters and config keys: accepted and ignored
MATLAB_ONLY_PARAMS: tuple[Param, ...] = (
    Param("mccMode", "mcc_mode", False, "any", "MATLAB only (compiled jobs); ignored"),
    Param("MatlabLaunchStr", "matlab_launch_str", "", "any", "MATLAB only; ignored"),
    Param("MCRCacheRoot", "mcr_cache_root", "", "any", "MATLAB only; ignored"),
    Param("MCRParam", "mcr_param", "", "any", "MATLAB only; ignored"),
    Param("MCCMasterStr", "mcc_master_str", "", "any", "MATLAB only; ignored"),
)

_BY_PYTHON = {p.python: p for p in SINGLE_JOB_PARAMS + MATLAB_ONLY_PARAMS}


# --- the session ---------------------------------------------------------------------------


class _Session:
    """this process's state of single jobs (MATLAB's persistent variables)"""

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        # the Slurm check and the job log names of the last job log folder: (key, -o, -e)
        self.log_check: tuple[str, str, str] | None = None
        # the kinds of single jobs this process submitted a task of
        self.kinds_seen: set[str] = set()
        # whether a task joined a pool (jobs.pools is imported only then)
        self.pools_used = False
        self.session_id = ""
        self.n_folders = 0  # folders made (plain and pools; mark/end count them)
        self.folders: dict[str, TaskFolder] = {}  # plain folders by parent folder
        self.pid = os.getpid()

    def new_folder_index(self) -> tuple[str, int]:
        """the session id and the index of a new folder"""
        if not self.session_id:
            # the time (also its microseconds) and this process's id, as MATLAB names pools
            now = time.time()
            stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime(now))
            self.session_id = f"{stamp}_{os.getpid()}_{int(now * 1e6) % 1000000:06d}"
            _register_exit()
        self.n_folders += 1
        return self.session_id, self.n_folders


_SESSION = _Session()
_EXIT_REGISTERED = [False]


def reset() -> None:
    """forget this process's single jobs: the Slurm check, the kinds, the folders and pools
    (MATLAB: clear generic_single_job_submit_wrapper task_runner_run single_job_pools)"""
    _SESSION.reset()
    pools = sys.modules.get("petakit5d.jobs.pools")
    if pools is not None:
        pools.reset()


def _register_exit() -> None:
    if not _EXIT_REGISTERED[0]:
        import atexit

        atexit.register(_leave_at_exit)
        _EXIT_REGISTERED[0] = True


def _leave_at_exit() -> None:
    # (a forked child that exits must not leave its parent's folders)
    if os.getpid() == _SESSION.pid:
        try:
            end(0)
        except Exception:
            pass


def _slurm(timeout: float | None = None) -> Slurm:
    """the Slurm commands (tests replace this)"""
    return Slurm() if timeout is None else Slurm(query_timeout=timeout)


def pools_used() -> bool:
    """whether a task of this process joined a pool (single_job_pools())"""
    return _SESSION.pools_used


def mark() -> int:
    """the folders this process made so far (single_job_pools('mark')), for ``end``"""
    return _SESSION.n_folders


def end(mark: int) -> None:
    """the end of a call whose tasks are done or given up (single_job_pools('end', mark)): the
    folders made after mark take no more tasks (a later task makes a new one), and this process
    leaves them (see the module docstring); a folder a job still in the queue has stays until
    that job ends"""
    ended = [f for f in _SESSION.folders.values() if f.index > mark]
    for f in ended:
        _SESSION.folders.pop(f.parent, None)
    if _SESSION.pools_used:
        from petakit5d.jobs import pools

        ended += pools.end(mark)
    if not ended:
        return
    # (one squeue for the jobs of all of them; MATLAB's unique() sorts the keys)
    keys = sorted({k for f in ended for k in f.jobs.values() if k})
    left: set[str] | None = set()
    if keys:
        left = _slurm().live_jobs(keys)
    for f in ended:
        f.leave_master(left)


# --- the folders of single jobs ------------------------------------------------------------


def _append(path: str, data: bytes) -> None:
    with open(path, "ab") as fh:
        fh.write(data)


def _touch(path: str) -> None:
    with open(path, "ab"):
        pass


class TaskFolder:
    """A folder of single jobs' tasks, written by this process and read by the jobs:

    call.json     the job entry point's call file (its "runner" runs a task by its line)
    tasks.ndjson  one JSON line per task: index, task (record), task_id, language,
                  bash_launch_str, metrics (the task metrics folder, or '')
    tasks.idx     the byte offset of each line (12 digits and a line break), so that a job
                  reads its line without the others; written after the line
    live/         m (this process) and t<line> per task, until its job ends
    """

    def __init__(self, path: str, parent: str, index: int) -> None:
        self.path = path
        self.parent = parent
        self.index = index  # this session's folder number (mark/end)
        self.n = 0
        self.bytes = 0
        self.jobs: dict[int, str] = {}  # line -> job key of its job

    @classmethod
    def make(
        cls,
        parent: str,
        prefix: str = TASK_FOLDER_PREFIX,
        runner: str = JOB_RUNNER,
        extra: Mapping[str, Any] | None = None,
    ) -> TaskFolder:
        """makes a new folder in parent (raises OSError if it can't)"""
        session, index = _SESSION.new_folder_index()
        path = os.path.join(parent, f"{prefix}{session}_{index}")
        os.makedirs(os.path.join(path, LIVE_FOLDER))
        for name in (TASKS_FILE, INDEX_FILE):
            _touch(os.path.join(path, name))
        _touch(os.path.join(path, LIVE_FOLDER, "m"))
        meta: dict[str, Any] = {
            "format": CALL_FORMAT,
            "version": CALL_VERSION,
            "n_tasks": 0,
            "language": "python",
            "task_batch_num": 1,
            "runner": runner,
            "offsets": [],
            "settings": {},
            "single": {"version": FOLDER_VERSION},
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),
            "python": sys.executable,
        }
        if extra:
            meta["single"].update(extra)
        tmp = os.path.join(path, f".{CALL_FILE}.{os.getpid()}.tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(meta, fh)
        os.replace(tmp, os.path.join(path, CALL_FILE))
        return cls(path, parent, index)

    def add(self, line: Mapping[str, Any]) -> int:
        """appends a task line; returns its 1-based number (raises OSError)"""
        i = self.n + 1
        data = (json.dumps({"index": i, **line}, allow_nan=True) + "\n").encode()
        _append(os.path.join(self.path, TASKS_FILE), data)
        _append(os.path.join(self.path, INDEX_FILE), f"{self.bytes:0{INDEX_WIDTH}d}\n".encode())
        _touch(os.path.join(self.path, LIVE_FOLDER, f"t{i}"))
        self.n = i
        self.bytes += len(data)
        return i

    def command(self, i: int, launch: str) -> str:
        """the job's command for task line i"""
        return job_command(self.path, i, launch)

    def leave_master(self, left: set[str] | None) -> None:
        """this process leaves the folder: its mark and those of tasks whose jobs Slurm no
        longer has (left: the job keys it has; None: no answer, then only its own) go, and if
        live/ is then empty the folder goes"""
        live = os.path.join(self.path, LIVE_FOLDER)
        marks = ["m"]
        if left is not None:
            marks += [f"t{i}" for i in range(1, self.n + 1) if self.jobs.get(i, "") not in left]
        for m in marks:
            try:
                os.remove(os.path.join(live, m))
            except OSError:
                pass
        _remove_if_left(self.path)


def job_command(folder: str, i: int, launch: str) -> str:
    """the command of the job of task line i of a folder of single jobs"""
    return f"{launch} -u -m petakit5d.jobs.exec {shlex.quote(folder)} {int(i)}"


def _remove_if_left(path: str) -> None:
    """removes a folder of single jobs whose live/ is empty (rmdir succeeds only then, which
    the file server decides, so one of the last to leave removes it)"""
    name = os.path.basename(os.path.normpath(path))
    if not name.startswith((TASK_FOLDER_PREFIX, POOL_PREFIX)):
        return
    try:
        os.rmdir(os.path.join(path, LIVE_FOLDER))
    except OSError:
        return
    shutil.rmtree(path, ignore_errors=True)


def leave(folder: str, i: int) -> None:
    """the job of task line i leaves its folder: its mark goes, and the last one to leave
    removes the folder"""
    try:
        os.remove(os.path.join(folder, LIVE_FOLDER, f"t{int(i)}"))
    except OSError:
        pass
    _remove_if_left(folder)


def read_task_line(folder: str, i: int) -> dict[str, Any]:
    """the line of task i (1-based) of a folder of single jobs"""
    width = INDEX_WIDTH + 1
    with open(os.path.join(folder, INDEX_FILE), "rb") as fh:
        fh.seek((int(i) - 1) * width)
        raw = fh.read(width)
    if len(raw) != width or int(i) < 1:
        raise ValueError(f"{folder} has no task {i}")
    with open(os.path.join(folder, TASKS_FILE), "rb") as fh:
        fh.seek(int(raw[:INDEX_WIDTH]))
        return json.loads(fh.readline())


def _plain_folder(parent: str) -> TaskFolder:
    """this session's folder of tasks that run as before in parent (made on first use)"""
    f = _SESSION.folders.get(parent)
    if f is None or not os.path.isdir(os.path.join(f.path, LIVE_FOLDER)):
        os.makedirs(parent, exist_ok=True)
        f = TaskFolder.make(parent)
        _SESSION.folders[parent] = f
    return f


# --- the job side ----------------------------------------------------------------------------


class _OneTask:
    """what metrics.task_start reads of a call, for one single job's task: its record under
    its task id (MATLAB records a single job's task as task task_id)"""

    def __init__(self, line: Mapping[str, Any]) -> None:
        self.settings = {"task_metrics": line.get("metrics") or ""}
        self.records = {int(line.get("task_id", 1)) - 1: line["task"]}
        self.index = int(line.get("task_id", 1)) - 1


def run_line(line: Mapping[str, Any], mode: str = "slurm_python") -> tuple[int, float]:
    """runs a task line in this process as a job runs its task; returns (exit status, seconds):
    a Python task's exception is printed with its traceback (status 1); a bash task runs in
    bash after the line's BashLaunchStr (its exit status)"""
    import traceback

    rec = line["task"]
    metrics = None
    if line.get("metrics"):
        try:
            from petakit5d.jobs import metrics
        except ImportError:
            metrics = None
    call = _OneTask(line)
    tm = metrics.task_start(call, [call.index], mode) if metrics else None
    t0 = time.perf_counter()
    if "bash" in rec:
        from petakit5d.jobs.framework import run_shell

        cmd = rec["bash"]
        launch = line.get("bash_launch_str") or ""
        status = run_shell(f"{launch}; {cmd}" if launch else cmd)
        if metrics:
            metrics.task_stop(tm, status)
        return status, time.perf_counter() - t0
    from petakit5d.jobs.tasks import run_record

    status = 0
    try:
        run_record(rec)
    except Exception as err:
        if metrics:
            metrics.task_stop(tm, err)
        sys.stdout.flush()
        print(f"Error in task {line.get('task_id')} ({rec.get('func')}):", file=sys.stderr)
        traceback.print_exc()
        status = 1
    except BaseException as err:
        if metrics:
            metrics.task_stop(tm, err)
        raise
    else:
        if metrics:
            metrics.task_stop(tm, None)
    return status, time.perf_counter() - t0


def run_job(folder: str, i: int) -> int:
    """The job entry point's runner for a task that runs as before (call.json "runner"): runs
    task line i in this process, then leaves the folder; returns the exit status."""
    t0 = time.perf_counter()
    try:
        line = read_task_line(folder, i)
        print(f"Task {line.get('task_id')}: {describe(line['task'])}", flush=True)
        status, _ = run_line(line)
        print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)
        return status
    finally:
        leave(folder, i)


# --- submission ------------------------------------------------------------------------------


def _check_types(settings: Mapping[str, Any]) -> None:
    # MATLAB's inputParser validators (islogical, isnumeric, ischar)
    for name, value in settings.items():
        kind = _BY_PYTHON[name].kind
        ok = True
        if kind == "bool":
            ok = isinstance(value, bool) or type(value).__name__ == "bool_"
        elif kind in ("num", "num_or_none"):
            ok = (value is None and kind == "num_or_none") or (
                isinstance(value, numbers.Real) and not isinstance(value, bool)
            )
        elif kind == "str":
            ok = isinstance(value, str)
        if not ok:
            raise TypeError(
                f"The value of '{_BY_PYTHON[name].matlab}' ({name}) is invalid: expected "
                f"{kind.replace('_or_none', ' or None')}, got {value!r}"
            )


def _apply_config(settings: dict[str, Any], config_file: str) -> dict[str, Any]:
    """the settings with the config file's values (keys matched without regard to case;
    MATLAB-only and unknown keys have no effect); the file is read again only when it
    changed, and then printed"""
    from petakit5d.jobs import config

    cfg, read_now = config._CACHE.load(config_file)
    if read_now:
        print(f"Set up cluster configuration according to {config_file}...")
        print("Cluster configuration:")
        for key, value in cfg.entries:
            print(f"    {key}: {value}")
        print()
    by_matlab = {p.matlab: p for p in SINGLE_JOB_PARAMS}
    out = dict(settings)
    for name, value in cfg.overrides(by_matlab).items():
        param = by_matlab[name]
        # (MATLAB's jsonencode writes an empty value as [], which jsondecode reads back)
        if isinstance(value, list) and not value and param.kind in ("num_or_none", "any"):
            value = None
        out[param.python] = value
    return out


def _round_half_away(x: float) -> int:
    return int(math.copysign(math.floor(abs(x) + 0.5), x))


def _task_key(rec: Mapping[str, Any]) -> str:
    """a task's identity in this session (MATLAB's function string)"""
    return json.dumps(rec, sort_keys=True, separators=(",", ":"), allow_nan=True)


def _alternatives_move(job_id: int, array_id: int, pending: bool) -> None:
    # SlurmAlternatives: a pending job may move, a gone one is forgotten (nothing to do unless
    # the module recorded jobs of this process, which needs it imported)
    alt = sys.modules.get("petakit5d.jobs.alternatives")
    if alt is None:
        return
    if pending:
        alt.move([job_id], [array_id])
    else:
        alt.forget([job_id], [array_id])


def _status_of(job_id: int, array_id: int, known: int | None) -> int:
    """the job's state: known (from a query of the caller's round), else one squeue (no
    answer: running, as MATLAB's check_slurm_job_status)"""
    if known is not None:
        return int(known)
    st = _slurm().job_status(job_id, array_id)
    return RUNNING if st is None else st


def _pools() -> Any:
    from petakit5d.jobs import pools

    return pools


def submit_single_job(
    task: Any,
    job_id: int,
    task_id: int,
    *,
    cluster_type: str = "slurm",
    parse_cluster: bool = True,
    master_compute: bool = False,
    last_file: bool = False,
    job_log_dir: str = DEFAULT_JOB_LOG_DIR,
    job_log_fname: str = DEFAULT_JOB_LOG_FNAME,
    job_error_fname: str = DEFAULT_JOB_ERROR_FNAME,
    tmp_dir: str = "",
    max_cpu_num: float = 24,
    min_cpu_num: float = 1,
    cpus_per_task: float = 1,
    mem_per_cpu: float = 20.9,
    mem_allocate: float | None = None,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 30,
    max_job_num: float = math.inf,
    task_batch_num: float = 1,
    python_launch_str: str = "",
    bash_launch_str: str = "",
    slurm_param: str = "",
    slurm_constraint: str = "",
    job_time_limit: float | None = 24,
    language: str = "python",
    gnu_parallel: bool = False,
    para_job_num: float = 1,
    config_file: str = "",
    task_metrics: bool = False,
    task_runner: bool = True,
    slurm_alternatives: Any = None,
    gpu_headroom: float | None = None,
    nice: float = 0,
    mcc_mode: Any = False,
    matlab_launch_str: Any = "",
    mcr_cache_root: Any = "",
    mcr_param: Any = "",
    mcc_master_str: Any = "",
    known_status: int | None = None,
    output: str | Sequence[str] | None = None,
) -> tuple[int, int, int]:
    """Submits a task as a Slurm job of its own unless its job is pending or running
    (generic_single_job_submit_wrapper); returns (job_id, job_status, submit_status).

    task        : a ``Task`` (or its record, a callable, a functools.partial); with
                  ``language='bash'`` a command string
    job_id      : the task's job of before (<= 0: none)
    task_id     : the caller's number of the task; the job's array index is task_id mod 5000

    job_status is 1 (running), 0 (pending) or -1 (no job in the queue); submit_status 1 if
    a job was submitted in this call (job_id is then the new job's), else 0. It never waits.
    With master_compute and last_file the caller runs the task itself: a pending job of it is
    cancelled, and nothing is submitted (job_status -1 or 0). A task that a task runner of
    another job has counts as running.

    The other parameters are those of MATLAB's wrapper in snake_case with its defaults
    (SINGLE_JOB_PARAMS lists them); a config file's values override them. Not in MATLAB:

    known_status : the job's state from a query of the caller's round (``SingleJobs.query``),
                   used instead of a squeue of the job
    output       : the task's result (a path, or paths that must all exist): once Slurm says
                   the task's job, and the job of a runner that had it, are gone, it is looked
                   at, and if it exists the task is done: nothing is submitted (job_status -1,
                   submit_status 0). A job writes its result before it leaves the queue, so
                   this closes the window in which a caller that looked for the result before
                   the job ended submits a finished task again (MATLAB's callers have it).
    """
    del mcc_mode, matlab_launch_str, mcr_cache_root, mcr_param, mcc_master_str
    given = dict(locals())
    s = {p.python: given[p.python] for p in SINGLE_JOB_PARAMS}
    _check_types(s)
    if s["config_file"]:
        s = _apply_config(s, s["config_file"])
        _check_types(s)
    if s["cluster_type"] != "slurm":
        raise ValueError(f"cluster_type must be 'slurm', not {s['cluster_type']!r}")
    lang = s["language"].lower()
    lang = "python" if lang == "matlab" else lang
    if lang not in ("python", "bash"):
        raise ValueError(f"language must be 'python' or 'bash', not {s['language']!r}")
    rec = task_record(task, lang)
    key = _task_key(rec)

    mc_last = bool(s["master_compute"]) and bool(s["last_file"])
    log_dir = s["job_log_dir"] or DEFAULT_JOB_LOG_DIR
    log_out = s["job_log_fname"] or DEFAULT_JOB_LOG_FNAME
    log_err = s["job_error_fname"] or DEFAULT_JOB_ERROR_FNAME
    main_param = s["slurm_param"]
    param = main_param
    if s["nice"] > 0:
        param = f"{param} --nice={_round_half_away(s['nice'])}"
    constraint = s["slurm_constraint"]
    # (the caller's False wins over the config's True)
    runners = bool(task_runner) and bool(s["task_runner"])
    launch = s["python_launch_str"] or default_python_launch_str()

    # the Slurm check and the job log names, kept for the job log folder while Slurm is there
    if os.path.isdir(log_dir):
        check_key = f"{os.getcwd()}|{log_dir}"
        if _SESSION.log_check is not None and _SESSION.log_check[0] == check_key:
            log_out, log_err = _SESSION.log_check[1], _SESSION.log_check[2]
        else:
            import tempfile

            check = check_slurm_cluster(tempfile.gettempdir(), log_dir, _slurm())
            # (without Slurm the job keeps the log names it has)
            if check.parse_cluster:
                log_out, log_err = check.job_log_fname, check.job_log_error_fname
                _SESSION.log_check = (check_key, log_out, log_err)

    cpus = s["cpus_per_task"]
    max_trial = s["max_trial_num"]
    places: Any = None
    places_made = False
    k = 0
    plain: tuple[TaskFolder, int] | None = None
    pools: Any = None
    pool_task: Any = None
    job_id = int(job_id)
    job_status = GONE
    submit_status = 0
    trial = 0
    while submit_status == 0 and trial < max_trial:
        array_id = int(task_id) % ARRAY_ID_MODULUS
        # (a task without a job yet: squeue could only say so)
        job_status = GONE
        if job_id > 0:
            job_status = _status_of(job_id, array_id, known_status if trial == 0 else None)
            if job_status == PENDING and not mc_last:
                _alternatives_move(job_id, array_id, True)
            elif job_status < 0:
                _alternatives_move(job_id, array_id, False)

        # task runners: a task of a pool whose job ended may still run in a runner of another
        # job (then it counts as running); before it is submitted again or run by the caller,
        # it is taken back from the pool, so that no runner takes it. While its job is
        # pending, a task a runner failed is taken back (the job cancelled) and submitted
        # again, and one whose runner ended without finishing it is left to the job
        pool_task = None
        if job_status < RUNNING and _SESSION.pools_used:
            pools = _pools()
            pool_task = pools.task(key)
        if pool_task is not None and job_status == PENDING and not mc_last:
            if (
                pool_task.job_key == job_key(job_id, array_id)
                and pools.status(job_id, array_id, PENDING) == GONE
            ):
                job_status = GONE
            pool_task = None
        if pool_task is not None:
            live, claimed = pools.live(key)
            if live:
                return job_id, RUNNING, submit_status
            # (a task added in this call is submitted again after a failed submission)
            if pool_task.state != "new" or claimed or mc_last:
                if not pools.take(key):
                    return job_id, RUNNING, submit_status
                pool_task = None

        # a task whose job (and runner) ended may have made its result just now
        if job_status == GONE and job_id > 0 and output and _outputs_exist(output):
            return job_id, GONE, submit_status
        # the caller runs the last task itself: its pending job is cancelled
        if job_status == PENDING and mc_last:
            r = _slurm().cancel([job_key(job_id, array_id)])
            if r.output:
                print(r.output, end="" if r.output.endswith("\n") else "\n")
        # a pending or running job is left to the caller's next round
        if job_status >= PENDING or mc_last:
            return job_id, job_status, submit_status

        # no job: submit one
        mem = s["mem_allocate"]
        if mem is not None and cpus * s["mem_per_cpu"] < mem:
            cpus = math.ceil(mem / s["mem_per_cpu"])
        cpus = min(s["max_cpu_num"], max(s["min_cpu_num"], cpus))
        time_option = time_str(s["job_time_limit"], param, TIME_FORMAT_SINGLE)
        if not places_made and s["slurm_alternatives"] not in (None, [], {}):
            places = _single_job_places(s, main_param, cpus, time_option)
            places_made = True
        if places is not None:
            k, cpus, param, constraint, time_option = places.choose()

        # the metrics folder of this process's single jobs ('' without metrics)
        run_dir = ""
        if s["task_metrics"]:
            from petakit5d.jobs import metrics

            run_dir = metrics.session_dir(os.path.dirname(log_out))
        line = {
            "task": rec,
            "task_id": int(task_id),
            "language": lang,
            "bash_launch_str": s["bash_launch_str"],
            "metrics": run_dir,
        }

        # task runners (Linux, Python tasks): a first submission of a task joins the pool of
        # this process's single jobs of its kind, made with the second task of the kind
        pool_dir = os.path.dirname(log_out)
        if (
            pool_task is None
            and runners
            and job_id <= 0
            and _LINUX
            and lang == "python"
            and os.path.isabs(pool_dir)
            and os.path.isdir(pool_dir)
        ):
            kind = "|".join(
                [rec["func"], launch, str(int(cpus)), param, constraint, time_option, pool_dir]
            )
            if kind not in _SESSION.kinds_seen:
                _SESSION.kinds_seen.add(kind)
            else:
                _SESSION.pools_used = True
                pools = _pools()
                pool_task = pools.add(key, kind, pool_dir, line, int(cpus))

        if pool_task is not None:
            cmd = pool_task.command(launch)
        else:
            if plain is None:
                plain = _add_plain(s["tmp_dir"] or pool_dir, line)
            if plain is None:
                print("Job submission failed, try again!")
                trial += 1
                continue
            cmd = plain[0].command(plain[1], launch)
        wrap = f"date; echo Python command: {shlex.quote(cmd)}; {cmd}"
        argv = [
            "sbatch",
            f"--array={array_id}",
            "-o",
            log_out,
            "-e",
            log_err,
            f"--cpus-per-task={int(cpus)}",
            *split_options(param),
            *split_options(constraint),
            *split_options(time_option),
            f"--wrap={wrap}",
        ]
        print(shlex.join(argv))
        r = _slurm().sbatch(argv)
        if r.result.output:
            out = r.result.output
            print(out, end="" if out.endswith("\n") else "\n")
        if r.submitted:
            job_id = int(r.job_id)
            submit_status = 1
            if places is not None:
                places.submitted(job_id, array_id, k)
            if pool_task is not None:
                pools.job(key, job_id, array_id)
            elif plain is not None:
                plain[0].jobs[plain[1]] = job_key(job_id, array_id)
            if run_dir:
                from petakit5d.jobs import metrics

                metrics.record(
                    run_dir,
                    "submit",
                    {
                        "task": int(task_id),
                        "jobId": job_id,
                        "functionStr": describe(rec),
                        "cpusPerTask": cpus,
                        "memAllocate": s["mem_allocate"],
                        "memPerCPU": s["mem_per_cpu"],
                        "mccMode": False,
                        "SlurmParam": param,
                    },
                )
        else:
            print("Job submission failed, try again!")
        trial += 1

    # a pool task whose job couldn't be submitted: taken back (the caller submits it again)
    if submit_status == 0 and pool_task is not None:
        pools.job(key, None, None)
    return job_id, job_status, submit_status


def _outputs_exist(output: str | Sequence[str]) -> bool:
    from petakit5d.jobs.files import batch_file_exist

    paths = [output] if isinstance(output, (str, os.PathLike)) else list(output)
    return all(batch_file_exist([os.fspath(x)])[0] for x in paths)


def _add_plain(parent: str, line: Mapping[str, Any]) -> tuple[TaskFolder, int] | None:
    """the task as a line of this session's folder in parent (made if needed), or None with a
    warning if it can't be written"""
    try:
        folder = _plain_folder(os.path.abspath(parent))
        return folder, folder.add(line)
    except OSError as err:
        warnings.warn(
            SingleJobWarning(
                f"The task of a single job can't be written to {parent}: {err}",
                "generic_single_job_submit_wrapper:taskFile",
            ),
            stacklevel=3,
        )
        return None


def _single_job_places(s: Mapping[str, Any], main_param: str, cpus: float, time_option: str):
    # SlurmAlternatives: the places of the job (None for one place), made at its first
    # submission in the call
    from petakit5d.jobs import alternatives

    return alternatives.single_job_places(
        alternatives=s["slurm_alternatives"],
        slurm_param=main_param,
        slurm_constraint=s["slurm_constraint"],
        mem_per_cpu=s["mem_per_cpu"],
        max_cpu_num=s["max_cpu_num"],
        min_cpu_num=s["min_cpu_num"],
        gpu_headroom=s["gpu_headroom"],
        requested_cpus=s["cpus_per_task"],
        mem_allocate=s["mem_allocate"],
        main_cpus=int(cpus),
        job_time_limit=s["job_time_limit"],
        main_time_option=time_option,
        nice=s["nice"],
    )


def check_single_job_status(job_id: int, task_id: int, *, known_status: int | None = None) -> int:
    """The state of a single job of task task_id (check_single_job_status): 1 running, 0
    pending, -1 gone. A task that a task runner of another job has counts as running; one that
    a runner failed while its job was pending counts as gone (the job is cancelled). With -1 no
    job or runner has the task, and none will take it: the caller may run or submit it again.
    A pending job submitted to one of several places (SlurmAlternatives) may move to another
    one where it can start now. known_status: as for submit_single_job."""
    job_id = int(job_id)
    array_id = int(task_id) % ARRAY_ID_MODULUS
    if job_id <= 0 and known_status is None:
        status = GONE
    else:
        status = _status_of(job_id, array_id, known_status)
    if status == PENDING:
        _alternatives_move(job_id, array_id, True)
    elif status < 0:
        _alternatives_move(job_id, array_id, False)
    # (a session without pools doesn't import them)
    if status < RUNNING and _SESSION.pools_used:
        status = _pools().status(job_id, array_id, status)
    return status


# --- the callers' side ---------------------------------------------------------------------


def job_states(
    job_ids: Sequence[int], array_ids: Sequence[int], slurm: Slurm | None = None
) -> list[int]:
    """The states of single jobs as check_slurm_job_status gives them (1 running, 0 pending or
    another state Slurm lists, -1 not listed or no job), from one ``squeue -j <keys>`` per 5000
    jobs; jobs of a query that got no answer count as running."""
    slurm = _slurm() if slurm is None else slurm
    keys = [job_key(j, a) if int(j) > 0 else None for j, a in zip(job_ids, array_ids, strict=True)]
    order = sorted({k for k in keys if k is not None}, key=lambda x: tuple(map(int, x.split("_"))))
    found: dict[str, int] = {}
    for b in range(0, len(order), 5000):
        part = order[b : b + 5000]
        r = slurm.run(["squeue", "-j", ",".join(part)])
        if "Invalid job id" in r.output:
            states: dict[str, int] | None = {}
        elif not r.ok:
            states = None
        else:
            states = {}
            wanted = set(part)
            # (the header first; a job is running if a line of it says R, as MATLAB's check
            # looks for ' R ' in the output of its one job)
            for row in r.stdout.splitlines()[1:]:
                f = row.split()
                if f and f[0] in wanted:
                    if f" {row} ".find(" R ") >= 0:
                        states[f[0]] = RUNNING
                    else:
                        states.setdefault(f[0], PENDING)
        for x in part:
            found[x] = RUNNING if states is None else states.get(x, GONE)
    return [GONE if k is None else found[k] for k in keys]


@dataclass
class SingleJobs:
    """The callers' side of single jobs (the loops of MATLAB's wrappers), for n items:

    ``query(done)``                one squeue for the jobs of the items not done (a round's
                                   start); their states are used by the next submit/check
    ``submit(i, task, task_id, **params)``
                                   submit_single_job while item i has tries left, else the last
                                   try's check_single_job_status; updates job_ids, trials and
                                   active and returns the job's state
    ``ran_here(i)``                the caller ran item i itself (a try)
    ``settled(done)``              every item done, or given up with no active job
    ``coming(done)``               the items not done and not given up
    ``wait(file_groups, n_left)``  wait_for_files between rounds: min(30, max(3, 0.3 n_left)) s
    ``with SingleJobs(n) as jobs`` the pools the call makes end with it (mark/end)
    """

    n: int
    max_trial_num: float = 3
    job_ids: list[int] = field(default_factory=list)
    task_ids: list[int] = field(default_factory=list)
    trials: list[int] = field(default_factory=list)
    active: list[bool] = field(default_factory=list)
    _known: dict[int, tuple[int, int]] = field(default_factory=dict, repr=False)
    _mark: int | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        n = int(self.n)
        self.job_ids = self.job_ids or [-1] * n
        self.task_ids = self.task_ids or list(range(1, n + 1))
        self.trials = self.trials or [0] * n
        self.active = self.active or [False] * n

    def __enter__(self) -> SingleJobs:
        self._mark = mark()
        return self

    def __exit__(self, *exc: Any) -> None:
        if self._mark is not None:
            end(self._mark)
            self._mark = None

    def query(self, done: Sequence[bool] | None = None) -> None:
        """one squeue for the jobs of the items not done (their states hold until used)"""
        items = [
            i for i in range(self.n) if self.job_ids[i] > 0 and not (done is not None and done[i])
        ]
        self._known = {}
        if not items:
            return
        states = job_states(
            [self.job_ids[i] for i in items],
            [self.task_ids[i] % ARRAY_ID_MODULUS for i in items],
        )
        for i, st in zip(items, states, strict=True):
            self._known[i] = (self.job_ids[i], st)

    def _take_known(self, i: int) -> int | None:
        hit = self._known.pop(i, None)
        if hit is None or hit[0] != self.job_ids[i]:
            return None
        return hit[1]

    def submit(self, i: int, task: Any, task_id: int | None = None, **params: Any) -> int:
        """submit_single_job for item i while it has tries left, else the last try's
        check_single_job_status; returns the job's state"""
        if task_id is not None:
            self.task_ids[i] = int(task_id)
        tid = self.task_ids[i]
        known = self._take_known(i)
        submit_status = 0
        if self.trials[i] < self.max_trial_num:
            job_id, job_status, submit_status = submit_single_job(
                task, self.job_ids[i], tid, known_status=known, **params
            )
            self.job_ids[i] = job_id
        else:
            job_status = check_single_job_status(self.job_ids[i], tid, known_status=known)
        self.trials[i] += submit_status
        self.active[i] = submit_status > 0 or job_status >= PENDING
        return job_status

    def check(self, i: int) -> int:
        """check_single_job_status for item i (also updates active)"""
        job_status = check_single_job_status(
            self.job_ids[i], self.task_ids[i], known_status=self._take_known(i)
        )
        self.active[i] = job_status >= PENDING
        return job_status

    def ran_here(self, i: int) -> None:
        self.trials[i] += 1

    def settled(self, done: Sequence[bool]) -> bool:
        return all(
            done[i] or (self.trials[i] >= self.max_trial_num and not self.active[i])
            for i in range(self.n)
        )

    def coming(self, done: Sequence[bool]) -> list[int]:
        return [
            i
            for i in range(self.n)
            if not done[i] and not (self.trials[i] >= self.max_trial_num and not self.active[i])
        ]

    def wait(
        self,
        file_groups: Sequence[Sequence[str]] | None,
        n_left: int,
        wait_time: float | None = None,
    ) -> bool:
        """wait_for_files between rounds (True if a group's files all exist): wait_time, or
        0.3 s per item left between 3 and 30 s"""
        from petakit5d.jobs.files import wait_for_files

        if wait_time is None:
            wait_time = min(WAIT_MAX, max(WAIT_MIN, WAIT_PER_ITEM * n_left))
        return wait_for_files(file_groups, wait_time)
