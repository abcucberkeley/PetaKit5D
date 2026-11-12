"""The cluster loop of the job framework: tasks as Slurm jobs.

``run_cluster(call)`` is run_tasks' cluster mode, the parseCluster branch of MATLAB's
slurm_cluster_generic_computing_wrapper, with the batch shaping of
generic_computing_frameworks_wrapper and, for jobs of whole nodes or with GNU parallel, the
variants of mcc_slurm_cluster_generic_computing_wrapper (in-job slots of
``petakit5d.jobs.slots`` replace GNU parallel). The decisions and timings are MATLAB's:

- one ``sbatch --array=<batch mod 5000> -o/-e <jobLogDir>/job_%A_%a.out/.err
  --cpus-per-task=N <SlurmParam> <SlurmConstraint> -t H:MM:00 --wrap=...`` per batch; the wrap
  names only the call folder and the batch (``python -m petakit5d.jobs.exec <folder> <b>``);
- a round of the loop per queryInterval: one batched ``squeue -j`` of the jobs, the outputs
  looked for by folder listing, then the batches in order;
- a batch whose job left the queue without its outputs is submitted again after 45 s (30 s
  with master compute), up to maxTrialNum tries, CPU jobs with CPUs x (try + 1);
- at most maxJobNum batches with a job at once (a batch counts until all its tasks are done);
- master compute: the master runs the last batch that waits (cancelling its pending job);
- a pending job whose outputs all exist is cancelled; the final output ends the run; the
  outputs are looked at once more after the loop; the call's jobs are cancelled at the end.

What differs from MATLAB, on purpose:

- A Slurm query that fails or times out is "unknown": the jobs keep their last known states
  (MATLAB read it as "all jobs gone" and submitted them again after the grace period).
- A submission without an answer (a timeout, or a lost connection to the controller, after
  which Slurm may have made the job) is "unknown": the batch is left alone until a query finds
  its job (by its ``--comment``, which every job of a call gets) or finds that there is none.
- A refused submission is tried again in the next round, not at once with the tasks saved to a
  file (the wrap no longer holds the tasks, which was why MATLAB retried that way).
- The call's jobs are cancelled also when the loop ends with an exception, and on SIGTERM
  (the master then exits with status 143); MATLAB left them in the queue.
- In-job slots give GPU slot k the k-th device of the job's CUDA_VISIBLE_DEVICES when it lists
  enough (MATLAB: device k - 1, another job's GPU where Slurm doesn't renumber them).
- Timing only: a pause ends once the outputs of all tasks not done exist (MATLAB saw them at
  its next query, up to two query intervals later); and while only jobs of last tries are left,
  the loop pauses until its next query (MATLAB looped without a pause).
- The last waiting batch is left out of the submissions only when the master takes it (MATLAB
  left out batch 1 with runExtraTasks too, which then never ran again).
"""

from __future__ import annotations

import importlib
import math
import os
import shlex
import shutil
import signal
import sys
import threading
import time
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from petakit5d.jobs import framework
from petakit5d.jobs.config import default_python_launch_str
from petakit5d.jobs.files import batch_file_exist, exist, fileparts
from petakit5d.jobs.shapes import (
    ARRAY_ID_MODULUS,
    TIME_FORMAT_BATCH,
    TIME_FORMAT_MCC,
    job_cpus,
    retry_cpus,
    time_str,
)
from petakit5d.jobs.slurm import GONE, PENDING, RUNNING, Slurm, check_slurm_cluster, split_options
from petakit5d.jobs.tasks import Call, describe, new_call_folder, write_call

__all__ = [
    "GRACE_AFTER_JOB",
    "GRACE_MASTER_COMPUTE",
    "JobShape",
    "job_shape",
    "run_cluster",
]

# --- timing constants (MATLAB's) ---------------------------------------------------------------

# s a batch whose job left the queue without its outputs waits for them (NFS shows a file made
# on another node late) before it is submitted again; with master compute 15 s less
GRACE_AFTER_JOB = 45.0
GRACE_MASTER_COMPUTE = 15.0
# the master takes a batch whose job left the queue (while no job is pending) only after
# max(MASTER_WAIT_MIN, MASTER_WAIT_FACTOR x the time of its last batch) s, at most MASTER_WAIT_MAX
MASTER_WAIT_MIN = 30.0
MASTER_WAIT_MAX = 180.0
MASTER_WAIT_FACTOR = 0.5
# the shortest pause of a round (s), and MATLAB's pause margin
PAUSE_MIN = 0.001
PAUSE_MARGIN = 0.001
# the outputs and jobs are looked at again every this many batches within a round (n_status_check)
STATUS_CHECK_BATCHES = 10000
# a pause ends once the outputs of all tasks not done exist: they are looked at every
# OUTPUT_LOOK_INTERVAL s, or OUTPUT_LOOK_COST_FACTOR x the time of a look if that is longer
# (MATLAB saw them only at its next query, up to two query intervals later)
OUTPUT_LOOK_INTERVAL = 0.5
OUTPUT_LOOK_COST_FACTOR = 10.0
# a submission prints its batch's task with the function text every this many batches
PRINT_FUNCTION_EVERY = 50

# --- Slurm answers -----------------------------------------------------------------------------

# errors of a failed sbatch after which Slurm may have made the job (the reply was lost)
SUBMIT_UNKNOWN_ERRORS = (
    "Socket timed out",
    "Unable to contact slurm controller",
    "Zero Bytes were transmitted or received",
    "Connection timed out",
    "Connection refused",
    "Transport endpoint is not connected",
    "Communication connection failure",
)
# the job script's clock (MATLAB's date_cmd_str): "Start time: " first and "End time: " last
DATE_CMD = 'date "+%Y-%m-%d %H:%M:%S.%6N"'
# exit status of a master stopped by a signal (after it cancelled its jobs)
SIGNAL_EXIT_BASE = 128

# the loop's clock and pause (tests replace them)
_clock = time.monotonic
_sleep = time.sleep


# --- the shape of the jobs ---------------------------------------------------------------------


@dataclass(frozen=True)
class JobShape:
    """CPUs and batches of a call's jobs (generic_computing_frameworks_wrapper)

    cpus_per_task  : CPUs of a batch's first try
    task_batch_num : tasks per batch (one job per batch)
    para_job_num   : tasks at once in a job (its in-job slots)
    slot_mode      : jobs of whole nodes or with GNU parallel (MATLAB's compiled mode): the
                     jobs run their tasks with petakit5d.jobs.slots, one process each
    parallel       : GNUparallel and para_job_num > 1 (several tasks at once)
    """

    cpus_per_task: int
    task_batch_num: int
    para_job_num: int
    slot_mode: bool
    parallel: bool


def _round(x: float) -> int:
    # MATLAB's round of a scalar (half away from zero)
    return int(math.copysign(math.floor(abs(x) + 0.5), x))


def _rem(x: float, y: float) -> float:
    # MATLAB's rem (rem(x, Inf) is x)
    if math.isinf(y):
        return x
    return math.fmod(x, y) if y != 0 else math.nan


def job_shape(settings: dict[str, Any], n_tasks: int) -> JobShape:
    """The CPUs, batch size and slots of the jobs, as generic_computing_frameworks_wrapper
    shapes them for cluster jobs: job_cpus (rounded to a factor of maxCPUNum unless jobs take
    whole nodes); whole-node jobs run floor(maxCPUNum / cpusPerTask) tasks at once (at most
    GPUNum GPU tasks) in maxCPUNum CPUs; batches are a multiple of the tasks at once, at least
    minBatchNum and paraJobNum x paraBatchNum tasks, split into minTaskJobNum jobs or more
    (whole-node jobs: at most twice the tasks at once), and with runExtraTasks a multiple of
    maxJobNum (+ 1 with master compute) jobs."""
    s = settings
    whole = bool(s.get("whole_node_job", False))
    gpu = bool(s.get("gpu_job", False))
    max_cpu = s.get("max_cpu_num", 24)
    cpus = job_cpus(
        s.get("cpus_per_task", 1),
        s.get("mem_allocate"),
        s.get("mem_per_cpu", 20.9),
        s.get("min_cpu_num", 1),
        max_cpu,
        not whole,
    )
    para = s.get("para_job_num", 1)
    if whole:
        para = math.floor(max_cpu / cpus)
        cpus = max_cpu
        if gpu:
            para = min(para, s.get("gpu_num", 0))
    if para < 1:
        raise ValueError(
            f"jobs would run {para} tasks at once (whole-node GPU jobs need gpu_num of 1 or more)"
        )
    n = n_tasks
    master = 1 if s.get("master_compute", True) else 0
    max_job = s.get("max_job_num", math.inf)
    min_task_job = s.get("min_task_job_num", 1)
    tb = max(s.get("task_batch_num", 1), s.get("min_batch_num", 1))
    tb = math.ceil(tb / para) * para
    tb = max(tb, para * s.get("para_batch_num", 1))
    if n / tb < min_task_job / para:
        tb = max(math.ceil(n / min_task_job), 1)
    if whole:
        jobs = max_job + master
        if s.get("run_extra_tasks", False) and _rem(math.ceil(n / tb), jobs) != 0:
            rounds = max(1, math.floor(n / (jobs * tb)))
            tb1 = math.ceil((n + tb * 0.1 * master) / (rounds * jobs)) if math.isfinite(jobs) else 0
            tb = max(min(tb, n), tb1)
        if gpu:
            tb = math.ceil(tb / para) * para
        else:
            tb = min(2 * tb, math.ceil(tb / para) * para)
            if math.ceil(n / tb) < min_task_job:
                tb = math.ceil(n / min_task_job)
    parallel = bool(s.get("gnu_parallel", False)) and para > 1
    return JobShape(int(cpus), max(int(tb), 1), int(para), whole or parallel, parallel)


# --- the seams of other tasks ------------------------------------------------------------------


def _plugin(module: str) -> Any | None:
    # an optional module of the framework (task runners, SlurmAlternatives), None if absent
    name = f"petakit5d.jobs.{module}"
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as err:
        if err.name == name:
            return None
        raise


class _NoRunners:
    """The master side of the task runners (J06) when there are none. A runner object has:

    call_extra() -> dict              keys for call.json (e.g. "runner")
    waiting(flags, timestamp)         after each query: per task, whether it waits for a job
    keep(indices, retry) -> bool      before a batch is submitted again: True while runners
                                      of other jobs hold its tasks
    master_tasks(indices) -> list     the tasks of its batch the master runs itself
    close()                           at the end
    """

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


class _OnePlace:
    """The places of the jobs (SlurmAlternatives, J09) when there is one. A places object has:

    choose(cpus, trial, batch) -> (cpus, slurm_param, slurm_constraint, time_str)
    move(job_keys)                    pending jobs that can start now elsewhere move there
    submitted(job_id, array_id, trial, batch)
    forget(job_keys)                  at the end
    """

    def __init__(self, slurm_param: str, slurm_constraint: str, time_option: str) -> None:
        self.place = (slurm_param, slurm_constraint, time_option)

    def choose(self, cpus: int, trial: int, batch: int = 0) -> tuple[int, str, str, str]:
        return (cpus, *self.place)

    def move(self, job_keys: Sequence[str]) -> None:
        pass

    def submitted(self, job_id: int, array_id: int, trial: int, batch: int = 0) -> None:
        pass

    def forget(self, job_keys: Sequence[str]) -> None:
        pass


# --- signals -----------------------------------------------------------------------------------


class MasterTerminated(SystemExit):
    """The master got SIGTERM; its jobs were cancelled (exit status 128 + the signal)."""


class _TermGuard:
    """SIGTERM raises MasterTerminated in the master (main thread only), so that the loop's
    finally cancels the call's jobs; further signals are ignored until it is done."""

    def __init__(self) -> None:
        self.previous: Any = None
        self.active = False

    def __enter__(self) -> _TermGuard:
        if threading.current_thread() is threading.main_thread() and hasattr(signal, "SIGTERM"):
            try:
                self.previous = signal.signal(signal.SIGTERM, self._handler)
                self.active = True
            except (ValueError, OSError):
                self.active = False
        return self

    def _handler(self, signum: int, frame: Any) -> None:
        signal.signal(signum, signal.SIG_IGN)
        raise MasterTerminated(SIGNAL_EXIT_BASE + signum)

    def __exit__(self, *exc: Any) -> None:
        if self.active:
            signal.signal(signal.SIGTERM, self.previous)
            self.active = False


# --- the loop ----------------------------------------------------------------------------------


def run_cluster(call: Call, slurm: Slurm | None = None) -> list[bool]:
    """Runs a call's tasks as Slurm jobs until their outputs exist; returns the done flags.

    Without Slurm (sinfo fails) the tasks run here (framework.run_local), with MATLAB's
    warning. The jobs of the call are cancelled at the end, also after an exception and on
    SIGTERM (then MasterTerminated, exit status 143, after the jobs are cancelled).
    """
    loop = _ClusterLoop(call, Slurm() if slurm is None else slurm)
    return loop.run()


class _ClusterLoop:
    """The state machine of slurm_cluster_generic_computing_wrapper (parseCluster true).

    Per task (0-based): done, out_exist, in_exist, trials, job id, status (now and in the
    previous query: RUNNING 1, PENDING 0, GONE -1) and the time of the last query that saw its
    job. Batches are 1-based (b), as in MATLAB; a batch's array id is b mod 5000.
    """

    def __init__(self, call: Call, slurm: Slurm) -> None:
        self.call = call
        self.slurm = slurm
        self.s = call.settings
        self.n = call.n
        self.job_ids: list[int] = []
        self.unknown: dict[int, str] = {}  # batch -> marker of a submission without answer
        self.cancelled = False

    # -- set-up ----------------------------------------------------------------------------

    def run(self) -> list[bool]:
        call, s, n = self.call, self.s, self.n
        final_out = s.get("final_out_fullpath", "")
        if final_out and batch_file_exist([final_out])[0]:
            print(f"The final output file {final_out} already exists!\n")
            return [True] * n

        first_in = call.inputs[0] if isinstance(call.inputs[0], str) else call.inputs[0][0]
        data_path = fileparts(first_in)[0]
        check = check_slurm_cluster(data_path, s.get("job_log_dir") or "../job_logs", self.slurm)
        if not check.parse_cluster:
            s["parse_cluster"] = False
            return framework.run_local(call)

        self.shape = job_shape(s, n)
        # (the places shape their jobs from the CPUs asked for, as MATLAB's do)
        requested_cpus = s.get("cpus_per_task", 1)
        s["cpus_per_task"] = self.shape.cpus_per_task
        s["task_batch_num"] = self.shape.task_batch_num
        s["para_job_num"] = self.shape.para_job_num
        self.bs = self.shape.task_batch_num
        self.nb = math.ceil(n / self.bs)
        fmt = TIME_FORMAT_MCC if self.shape.slot_mode else TIME_FORMAT_BATCH
        self.time_option = time_str(s.get("job_time_limit"), s.get("slurm_param", ""), fmt)
        self.log_out, self.log_err = check.job_log_fname, check.job_log_error_fname
        self.log_dir = check.job_log_dir

        outputs = call.outputs
        self.out_exist = batch_file_exist(outputs, None, True)
        if all(self.out_exist):
            print(f"All output files ({n} / {n}) already exist(s)!\n")
            return [True] * n
        self.done = list(self.out_exist)
        self.in_exist = framework.input_exist(call, self.out_exist)

        self.master_compute = bool(s.get("master_compute", True))
        self.gpu = bool(s.get("gpu_job", False))
        self.max_trial = s.get("max_trial_num", 3)
        self.max_job = s.get("max_job_num", math.inf)
        self.query_interval = s.get("query_interval", 3)
        self.extra_tasks = self.shape.slot_mode and bool(s.get("run_extra_tasks", False))
        self.launch = s.get("python_launch_str") or default_python_launch_str()

        self.runners = self._runners()
        self.places = self._places(requested_cpus)
        parent = s.get("tmp_dir") or self.log_dir
        self.folder = new_call_folder(parent, s.get("uuid") or "call")
        self.marker_base = os.path.basename(self.folder)
        extra = dict(self.runners.call_extra())
        if self.shape.slot_mode:
            extra["slots"] = {
                "extra_tasks": self.extra_tasks,
                "outputs_existed": any(self.out_exist),
            }
        write_call(call, self.folder, extra)

        self.trials = [0] * n
        self.job_ids = [-1] * n
        self.cur = [GONE] * n
        self.prev = [GONE] * n
        self.job_ts = [0.0] * n
        self.task_ids = [((i // self.bs) + 1) % ARRAY_ID_MODULUS for i in range(n)]
        if self.shape.slot_mode:
            print(
                f"Task number : {n}, task batch size : {self.bs}, task parallel number : "
                f"{self.shape.para_job_num}, job number : {self.nb} "
            )
        else:
            print(f"Task number : {n}, task batch size : {self.bs}, job number : {self.nb}")

        with _TermGuard():
            try:
                self._loop()
                self._finish()
            finally:
                if not self.cancelled:
                    self._cancel_all()
                self.runners.close()
        return self.done

    def _runners(self) -> Any:
        s = self.s
        if not (s.get("task_runner", True) and sys.platform.startswith("linux")):
            return _NoRunners()
        mod = _plugin("runners")
        if mod is None:
            return _NoRunners()
        return mod.cluster_runners(self.call, self.shape)

    def _places(self, requested_cpus: float = 1) -> Any:
        s = self.s
        one = _OnePlace(s.get("slurm_param", ""), s.get("slurm_constraint", ""), self.time_option)
        if s.get("slurm_alternatives") in (None, [], {}) or s.get("whole_node_job"):
            return one
        mod = _plugin("alternatives")
        if mod is None:
            warnings.warn(
                "SlurmAlternatives needs petakit5d.jobs.alternatives; the jobs go to the main "
                "place",
                stacklevel=3,
            )
            return one
        places = mod.cluster_places(self.call, self.shape, self.time_option, requested_cpus)
        return places or one

    # -- the rounds ------------------------------------------------------------------------

    def _loop(self) -> None:
        n, s = self.n, self.s
        final_out = s.get("final_out_fullpath", "")
        self.t_start = _clock()
        ts = self.t_start
        self.qts = ts
        self.timestamp = 0.0
        self.pending_flag = False
        self.t1 = 0.0
        self.look_every = OUTPUT_LOOK_INTERVAL
        self.loop_counter = 0
        n_done = 0
        while self._more():
            query = self.loop_counter == 0 or _clock() - self.qts > self.query_interval
            if query:
                self._query()
            if self.loop_counter > 0 and query:
                self._look_at_outputs([i for i in range(n) if not self.done[i]])
            self._round()

            # (MATLAB pauses only while a task has tries left, and otherwise loops without a
            # pause until the jobs of the last tries leave the queue)
            if self._more():
                wt = max(PAUSE_MIN, self.query_interval - (_clock() - self.qts) - PAUSE_MARGIN)
                # all jobs pending: wait another interval
                if not self.master_compute and all(
                    self.cur[i] == PENDING for i in range(n) if not self.done[i]
                ):
                    wt += self.query_interval
                self._pause(wt)
            if n_done < sum(self.done):
                n_done = sum(self.done)
                if n_done == n:
                    print()
                print(
                    f"Time {_clock() - ts:0.2f} s: {n_done} / {n} "
                    f"({n_done / n * 100:0.2f}%) are finished!"
                )
            # batch_file_exist sees a final output another node just made
            if final_out and batch_file_exist([final_out])[0]:
                self.done = [True] * n
                print(
                    f"Time {_clock() - ts:0.2f} s: The final output file {final_out} "
                    "already exists!"
                )
                break
            self.loop_counter += 1

    def _pause(self, seconds: float) -> None:
        """MATLAB's pause of a round; it ends early once the outputs of all tasks not done
        (or the final output) exist, which then count as done"""
        end = _clock() + seconds
        while True:
            left = end - _clock()
            if left <= 0:
                return
            _sleep(min(left, self.look_every))
            if _clock() >= end:
                return
            t0 = time.perf_counter()
            if self._all_outputs_now():
                return
            cost = time.perf_counter() - t0
            self.look_every = max(OUTPUT_LOOK_INTERVAL, OUTPUT_LOOK_COST_FACTOR * cost)

    def _all_outputs_now(self) -> bool:
        final_out = self.s.get("final_out_fullpath", "")
        if final_out and batch_file_exist([final_out])[0]:
            return True
        undone = [i for i in range(self.n) if not self.done[i]]
        if not undone:
            return True
        found = batch_file_exist([self.call.outputs[i] for i in undone], None, True)
        if not all(found):
            return False
        for i in undone:
            self.out_exist[i] = True
            self.done[i] = True
        return True

    def _more(self) -> bool:
        # a task is not done and has tries left, or its job of the last try is in the queue
        return not all(
            self.done[i] or (self.trials[i] >= self.max_trial and self.cur[i] < 0)
            for i in range(self.n)
        )

    def _query(self, indices: Sequence[int] | None = None) -> None:
        """squeue -j of the jobs of the tasks (default: those not done); an unknown answer
        keeps their last states"""
        n = self.n
        undone = [i for i in range(n) if not self.done[i]] if indices is None else list(indices)
        if indices is None:
            for i in undone:
                self.prev[i] = self.cur[i]
        states = self.slurm.batch_status(
            [self.job_ids[i] for i in undone], [self.task_ids[i] for i in undone]
        )
        lost = False
        for i, st in zip(undone, states, strict=True):
            if st is None:
                lost = True
            else:
                self.cur[i] = st
        if lost:
            print(
                "Warning: squeue did not answer; the jobs keep their last known states.",
                file=sys.stderr,
            )
        if self.unknown:
            self._find_unknown()
        self.timestamp = _clock() - self.t_start
        if indices is None:
            for i in range(n):
                if self.cur[i] >= 0:
                    self.job_ts[i] = self.timestamp
            undone_now = [i for i in range(n) if not self.done[i]]
        else:
            for i in undone:
                self.job_ts[i] = self.timestamp
            undone_now = [i for i in range(n) if not self.done[i]]
        self.pending_flag = any(self.cur[i] == PENDING for i in undone_now)
        self.qts = _clock()
        if indices is None:
            if self.pending_flag:
                keys, seen = [], set()
                for i in undone_now:
                    if self.cur[i] == PENDING and self.job_ids[i] > 0:
                        key = f"{self.job_ids[i]}_{self.task_ids[i]}"
                        if key not in seen:
                            seen.add(key)
                            keys.append(key)
                if keys:
                    self.places.move(keys)
            self.runners.waiting(
                [
                    not self.done[i] and (self.cur[i] == PENDING or self.job_ids[i] < 0)
                    for i in range(n)
                ],
                self.timestamp,
            )

    def _look_at_outputs(self, indices: Sequence[int]) -> None:
        if not indices:
            return
        found = batch_file_exist([self.call.outputs[i] for i in indices], None, True)
        for i, ok in zip(indices, found, strict=True):
            self.out_exist[i] = ok
            self.done[i] = ok

    def _batch(self, b: int) -> range:
        return range((b - 1) * self.bs, min(b * self.bs, self.n))

    def _batch_done(self, b: int) -> bool:
        return all(self.done[f] for f in self._batch(b))

    def _last_waiting_batch(self) -> int:
        # the batch of the last task not done, with tries left, whose job doesn't run
        for i in range(self.n - 1, -1, -1):
            if (
                not self.done[i]
                and self.trials[i] < self.max_trial
                and self.cur[i] < RUNNING
                and (i // self.bs) + 1 not in self.unknown
            ):
                return (i // self.bs) + 1
        return -1

    def _round(self) -> None:
        n, call = self.n, self.call
        lastp = self._last_waiting_batch()
        fsnames = [""] * n
        active = None
        if self.nb > self.max_job:
            active = self._active_counter()
        for b in range(1, self.nb + 1):
            fs = self._batch(b)
            task_id = b % ARRAY_ID_MODULUS
            undone = [f for f in fs if not self.done[f]]
            if undone:
                mc = min(self.cur[f] for f in undone)
                mp = min(self.prev[f] for f in undone)
                for f in fs:
                    self.cur[f] = mc
                    self.prev[f] = mp
                if active is not None:
                    active.refresh(b)

            if b % STATUS_CHECK_BATCHES == 0 and b + STATUS_CHECK_BATCHES * 0.5 < self.nb:
                lastp = self._status_check(b, lastp)
                if active is not None:
                    active = self._active_counter()

            if all(self.done[f] or self.trials[f] >= self.max_trial for f in fs):
                continue

            for f in fs:
                inp = call.inputs[f]
                if isinstance(inp, list):
                    if not all(exist(x) for x in inp):
                        self.done[f] = True
                elif not self.in_exist[f]:
                    self.done[f] = True
                _, fsname, ext = fileparts(call.outputs[f])
                fsnames[f] = fsname if fsname or not ext else ext
                if self.out_exist[f]:
                    self.done[f] = True
                    # a pending job whose outputs all exist is cancelled
                    if (
                        (len(fs) == 1 or all(self.out_exist[x] for x in fs))
                        and self.job_ids[f] > 0
                        and self.cur[f] == PENDING
                    ):
                        self._scancel([f"{self.job_ids[f]}_{task_id}"])
            if active is not None:
                active.refresh(b)

            if all(self.done[f] or self.trials[f] >= self.max_trial for f in fs):
                continue
            f = fs[-1]
            names = [fsnames[x] for x in fs]
            master_batch = (
                self.master_compute
                and b == lastp
                and not (self.extra_tasks and self.nb > 1 and b == 1)
            )

            skip_submission = active is not None and active.count >= self.max_job
            if not skip_submission:
                # the master takes the last waiting batch: its pending job is cancelled
                if self.cur[f] == PENDING and master_batch and self.job_ids[f] > 0:
                    self._scancel([f"{self.job_ids[f]}_{task_id}"])
                    self.pending_flag = True
                    for x in fs:
                        self.trials[x] -= 1
                        self.cur[x] = GONE
                        self.prev[x] = GONE
                    if active is not None:
                        active.refresh(b)
                if self.cur[f] == RUNNING:
                    continue
                # a job that just left the queue: its outputs may show late
                grace = GRACE_AFTER_JOB - self.master_compute * GRACE_MASTER_COMPUTE
                if (
                    self.loop_counter > 0
                    and self.trials[f] > 0
                    and self.cur[f] == GONE
                    and self.timestamp - self.job_ts[f] < grace
                ):
                    continue
                # (MATLAB leaves out the last waiting batch also when the master doesn't take
                # it, batch 1 with runExtraTasks, which then never ran again)
                if (
                    self.cur[f] == GONE
                    and self.prev[f] == GONE
                    and not master_batch
                    and b not in self.unknown
                ):
                    if self.loop_counter > 0 and self.runners.keep(
                        [x for x in fs if not self.done[x]], self.trials[f] > 0
                    ):
                        continue
                    self._submit(b, fs, task_id, names)
                    if active is not None:
                        active.refresh(b)

            if master_batch:
                if self.loop_counter > 0 and not self.pending_flag:
                    wait = min(max(MASTER_WAIT_MIN, self.t1 * MASTER_WAIT_FACTOR), MASTER_WAIT_MAX)
                    if self.timestamp - self.job_ts[f] < wait:
                        continue
                    if all(self.done[x] for x in fs):
                        continue
                self._master_run(b, fs, names)

    def _active_counter(self) -> _ActiveJobs:
        return _ActiveJobs(self)

    def _status_check(self, b: int, lastp: int) -> int:
        """the look at outputs and jobs within a round (every STATUS_CHECK_BATCHES batches),
        as MATLAB has it (its batch numbers index the tasks there)"""
        if not _clock() - self.qts > self.query_interval:
            return lastp
        n = self.n
        self.prev = list(self.cur)
        job_inds = [i for i in range(n) if not self.done[i] and self.cur[i] > GONE]
        self._look_at_outputs(job_inds)
        self._query(job_inds)
        if self.master_compute:
            if b + STATUS_CHECK_BATCHES * 1.5 > self.nb:
                check = range(b, self.nb + 1)
            else:
                check = range(b, min(self.nb, b + STATUS_CHECK_BATCHES - 1) + 1)
            for k in reversed(check):
                i = k - 1
                if (
                    i < n
                    and not self.done[i]
                    and self.trials[i] < self.max_trial
                    and self.cur[i] < PENDING
                ):
                    return k
        return lastp

    # -- submission --------------------------------------------------------------------------

    def _job_cpus(self, trial: int) -> int:
        cpus, max_cpu = self.shape.cpus_per_task, self.s.get("max_cpu_num", 24)
        if self.shape.slot_mode:
            # compiled mode: CPU jobs keep their CPUs and run fewer tasks at once, GPU jobs
            # get more CPUs
            return int(min(max_cpu, cpus * (trial + 1) if self.gpu else cpus))
        return int(retry_cpus(cpus, max_cpu, trial, self.gpu))

    def _job_slots(self, trial: int) -> tuple[int, int, bool]:
        """(tasks at once, tries per task, GPU per slot) of a slot job's try"""
        para = self.shape.para_job_num
        if not self.shape.parallel or trial > 1:
            return 1, 1, False
        if self.gpu:
            return para, 2, True
        return max(1, _round(para / (trial + 1))), 2, False

    def _marker(self, b: int) -> str:
        return f"{self.marker_base}:{b}"

    def _wrap(self, b: int, trial: int) -> str:
        folder = shlex.quote(self.folder)
        if self.shape.slot_mode:
            slots, tries, gpu = self._job_slots(trial)
            cmd = (
                f"{self.launch} -u -m petakit5d.jobs.slots {folder} {b} --slots {slots} "
                f"--tries {tries}" + (" --gpu" if gpu else "")
            )
        else:
            cmd = f"{self.launch} -u -m petakit5d.jobs.exec {folder} {b}"
        where = "echo $PWD; " if self.shape.slot_mode or self.call.language == "bash" else ""
        return (
            f'printf "Start time: "; {DATE_CMD}; {where}echo Python command: '
            f'{shlex.quote(cmd)}; {cmd}; echo; printf "End time: "; {DATE_CMD}'
        )

    def _submit(self, b: int, fs: range, task_id: int, names: list[str]) -> None:
        f = fs[-1]
        trial = self.trials[f]
        text = ";".join(describe(self.call.records[x]) for x in fs)
        if b % PRINT_FUNCTION_EVERY == 0 or b == 1:
            print(f"Task {b: 4d}:    Process {', '.join(names)} with function {text}... ")
        else:
            print(f"Task {b: 4d}:    Process {', '.join(names)} ... ")
        place = self.places.choose(self._job_cpus(trial), trial, b)
        cpus, slurm_param, constraint, time_option = place
        argv = [
            "sbatch",
            f"--array={task_id}",
            "-o",
            self.log_out,
            "-e",
            self.log_err,
            f"--cpus-per-task={int(cpus)}",
        ]
        if self.shape.slot_mode and self.shape.parallel and trial <= 1:
            argv.append("--ntasks=1")
        argv += split_options(slurm_param) + split_options(constraint) + split_options(time_option)
        argv += [f"--comment={self._marker(b)}", f"--wrap={self._wrap(b, trial)}"]
        r = self.slurm.sbatch(argv)
        out = r.result.output
        if out:
            print(out, end="" if out.endswith("\n") else "\n")
        state = r.state
        if state == "failed" and any(e in out for e in SUBMIT_UNKNOWN_ERRORS):
            state = "unknown"
        if state == "submitted":
            self._submitted(b, fs, int(r.job_id), PENDING, trial)
        elif state == "unknown":
            print(
                f"Warning: sbatch gave no answer for batch {b}; its job is looked for before "
                "it is submitted again.",
                file=sys.stderr,
            )
            self.unknown[b] = self._marker(b)
        else:
            print("Unable to submit the job; it is tried again in the next round.")

    def _submitted(self, b: int, fs: range, job_id: int, state: int, trial: int) -> None:
        for x in fs:
            self.job_ids[x] = job_id
            self.cur[x] = state
            self.trials[x] += 1
        self.places.submitted(job_id, b % ARRAY_ID_MODULUS, trial, b)

    def _find_unknown(self) -> None:
        """looks for the jobs of submissions without answer by their comments"""
        found = self._lookup()
        if found is None:
            return
        for b in sorted(self.unknown):
            fs = self._batch(b)
            hit = found.get(self.unknown[b])
            del self.unknown[b]
            if hit is not None:
                print(f"Batch {b}: its job {hit[0]} was submitted.")
                self._submitted(b, fs, hit[0], hit[1], self.trials[fs[-1]])

    def _lookup(self) -> dict[str, tuple[int, int]] | None:
        # squeue --me -h -r -o '%i|%t|%k': {comment: (job id, state)} of this user's jobs
        r = self.slurm.run(["squeue", "--me", "-h", "-r", "-o", "%i|%t|%k"])
        if not r.ok:
            return None
        found: dict[str, tuple[int, int]] = {}
        for line in r.stdout.splitlines():
            parts = line.strip().split("|", 2)
            if len(parts) != 3 or not parts[0].split("_")[0].isdigit():
                continue
            st = RUNNING if parts[1].strip() == "R" else PENDING
            found[parts[2].strip()] = (int(parts[0].split("_")[0]), st)
        return found

    def _scancel(self, keys: Sequence[str]) -> None:
        r = self.slurm.cancel(keys)
        if r.output:
            print(r.output, end="" if r.output.endswith("\n") else "\n")

    # -- the master's batch ------------------------------------------------------------------

    def _master_run(self, b: int, fs: range, names: list[str]) -> None:
        # (with task runners, the tasks of the batch that no job's runner has)
        indices = self.runners.master_tasks(list(fs))
        if not indices:
            return
        if self.shape.slot_mode:
            from petakit5d.jobs import slots

            text = ";".join(describe(self.call.records[x]) for x in fs)
            print(f"\nProcess task batch {b} ... ")
            print(f"\nProcess {', '.join(names)} with function {text}... ")
            para, factor = self.shape.para_job_num, self.s.get("master_para_factor", 1)
            if not self.shape.parallel:
                n_slots, tries = 1, 1
            elif self.gpu:
                n_slots, tries = max(1, min(para, _round(para * factor))), 2
            else:
                trial = self.trials[fs[-1]]
                n_slots = max(1, min(para - 1, _round(para * factor / (trial + 1))))
                tries = 2
            t0 = time.perf_counter()
            slots.run_batch(self.folder, b, n_slots, tries=tries, gpu=False)
            self.t1 = time.perf_counter() - t0
        else:
            self.t1 = framework.run_batch_in_process(
                self.call, indices, [names[x - fs[0]] for x in indices], mode="master"
            )
        for x in fs:
            self.trials[x] += 1
        print(f"Done! Elapsed time is {self.t1:f} seconds.")

    # -- the end -----------------------------------------------------------------------------

    def _finish(self) -> None:
        n = self.n
        # a batch the master ran in its last try ran after the round looked at its outputs
        if not all(self.done):
            self._look_at_outputs([i for i in range(n) if not self.done[i]])
        self._cancel_all()
        pairs = zip(self.job_ids, self.task_ids, strict=True)
        keys = sorted({f"{j}_{a}" for j, a in pairs if j > 0})
        self.places.forget(keys)
        if all(self.done):
            print(
                f"Time {_clock() - self.t_start:0.2f} s: All output files ({n} / {n}) are "
                "finished!\n"
            )
            try:
                shutil.rmtree(self.folder)
            except OSError as err:
                warnings.warn(
                    f"The task files in {self.folder} could not be removed: {err}", stacklevel=3
                )

    def _cancel_all(self) -> None:
        """scancel of every job of the call (MATLAB: at a normal end only)"""
        self.cancelled = True
        if self.unknown:
            found = self._lookup()
            for b, marker in sorted(self.unknown.items()):
                hit = None if found is None else found.get(marker)
                if hit is not None:
                    for x in self._batch(b):
                        self.job_ids[x] = hit[0]
                elif found is None:
                    print(
                        f"Warning: a job of batch {b} may be left in the queue (comment {marker}).",
                        file=sys.stderr,
                    )
            self.unknown.clear()
        ids = sorted({j for j in self.job_ids if j > 0})
        if ids:
            self._scancel([str(j) for j in ids])


class _ActiveJobs:
    """the count of batches with a job and tasks to do (MATLAB counts it for every batch over
    all batches; here each batch's part is updated when the round changes it)"""

    def __init__(self, loop: _ClusterLoop) -> None:
        self.loop = loop
        self.part = [False] * (loop.nb + 1)
        for b in range(1, loop.nb + 1):
            self.part[b] = self._active(b)
        self.count = sum(self.part)

    def _active(self, b: int) -> bool:
        lp = self.loop
        first = (b - 1) * lp.bs
        return (lp.cur[first] >= 0 and not lp._batch_done(b)) or b in lp.unknown

    def refresh(self, b: int) -> None:
        now = self._active(b)
        if now != self.part[b]:
            self.count += 1 if now else -1
            self.part[b] = now
