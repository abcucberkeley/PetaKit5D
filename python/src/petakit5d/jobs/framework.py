"""run_tasks: the job framework's entry point (MATLAB's generic_computing_frameworks_wrapper),
and its sequential local mode (slurm_cluster_generic_computing_wrapper with parseCluster
false).

Area wrappers build lists of inputs, outputs and tasks (``tasks.Task`` records, or bash
commands with ``language='bash'``) and call run_tasks; they never call Slurm themselves.
A task is done when its output exists; run_tasks returns one done flag per task.

Backends: the sequential loop here; the local process pool (``jobs.local.run_pool``, task
J04), which runs first when ``local_worker_num > 1`` and leaves what it didn't finish to
the loop; the cluster loop (``jobs.cluster.run_cluster``, task J03). Each gets a
``tasks.Call``.
"""

from __future__ import annotations

import math
import numbers
import os
import subprocess
import sys
import time
import traceback
from collections.abc import Mapping, Sequence
from typing import Any, NamedTuple

from petakit5d.jobs.files import batch_file_exist, exist, fileparts
from petakit5d.jobs.tasks import Call, Task, describe, run_record, task_record

# the local loop sleeps this long after a round that ran nothing (it waits for another
# caller's .tmp marker); MATLAB loops without a pause
LOCAL_IDLE_WAIT = 0.2
# the local loop also looks at the outputs every this many batches (MATLAB's n_status_check)
STATUS_CHECK_BATCHES = 10000
# the default job log folder (relative to the current folder; see J01 for the fallback)
DEFAULT_JOB_LOG_DIR = "../job_logs"


class Param(NamedTuple):
    """One parameter: MATLAB name, Python name, default, kind (bool, num, num_or_none,
    int, str, any), meaning. STAND-IN for the parameter registry of F03."""

    matlab: str
    python: str
    default: Any
    kind: str
    doc: str


# the parameters of generic_computing_frameworks_wrapper that keep a meaning, MATLAB's
# order, names and defaults (site defaults are empty: they come from the config file)
RUN_TASKS_PARAMS: tuple[Param, ...] = (
    Param(
        "finalOutFullpath",
        "final_out_fullpath",
        "",
        "str",
        "stop (all tasks done) once this path exists",
    ),
    Param(
        "clusterType",
        "cluster_type",
        "slurm",
        "str",
        "'slurm' (cluster jobs or local) or 'parfor' (local process pool)",
    ),
    Param(
        "parseCluster",
        "parse_cluster",
        True,
        "bool",
        "run the tasks as Slurm jobs (if there is a cluster); else in this process",
    ),
    Param(
        "masterCompute",
        "master_compute",
        True,
        "bool",
        "this process runs tasks too (in cluster mode: the last waiting batch)",
    ),
    Param(
        "jobLogDir",
        "job_log_dir",
        DEFAULT_JOB_LOG_DIR,
        "str",
        "folder of the job logs (empty: the default)",
    ),
    Param(
        "tmpDir",
        "tmp_dir",
        "",
        "str",
        "folder of the .tmp markers of local runs and the call folders (empty: the "
        "output folder / the job log folder)",
    ),
    Param("maxCPUNum", "max_cpu_num", 24, "num", "most CPUs of a job"),
    Param("minCPUNum", "min_cpu_num", 1, "num", "fewest CPUs of a job"),
    Param("cpusPerTask", "cpus_per_task", 1, "num", "CPUs of a task's job"),
    Param("memPerCPU", "mem_per_cpu", 20.9, "num", "GB of memory per CPU of the jobs"),
    Param(
        "memAllocate",
        "mem_allocate",
        None,
        "num_or_none",
        "GB of memory a task needs (None: unknown); sets the CPUs of its job",
    ),
    Param(
        "wholeNodeJob",
        "whole_node_job",
        False,
        "bool",
        "jobs take whole nodes and run several tasks at once (in-job slots)",
    ),
    Param("uuid", "uuid", "", "str", "id of this call (empty: a new one)"),
    Param("GPUJob", "gpu_job", False, "bool", "the tasks use a GPU"),
    Param("GPUNum", "gpu_num", 0, "num", "GPUs of a whole-node job"),
    Param("maxTrialNum", "max_trial_num", 3, "num", "tries of a task"),
    Param(
        "unitWaitTime",
        "unit_wait_time",
        3,
        "num",
        "minutes a local run leaves a task to another caller whose .tmp marker is newer",
    ),
    Param("maxJobNum", "max_job_num", math.inf, "num", "most jobs pending or running"),
    Param(
        "minTaskJobNum",
        "min_task_job_num",
        1,
        "num",
        "split the tasks into at least this many jobs (whole-node jobs)",
    ),
    Param("taskBatchNum", "task_batch_num", 1, "num", "tasks per batch (one job per batch)"),
    Param("minBatchNum", "min_batch_num", 1, "num", "fewest tasks per batch (whole-node jobs)"),
    Param("paraBatchNum", "para_batch_num", 1, "num", "batches per slot of a whole-node job"),
    Param(
        "runExtraTasks",
        "run_extra_tasks",
        False,
        "bool",
        "whole-node jobs: shape batches to a multiple of maxJobNum",
    ),
    Param(
        "PythonLaunchStr",
        "python_launch_str",
        "",
        "str",
        "how a job starts Python: module loads, environment activation, exports, then "
        "the interpreter (empty: the interpreter of this process)",
    ),
    Param(
        "BashLaunchStr",
        "bash_launch_str",
        "",
        "str",
        "shell commands a job runs before its bash tasks",
    ),
    Param(
        "SlurmParam",
        "slurm_param",
        "",
        "str",
        "sbatch options of the jobs (partition, QOS, memory per CPU); from the config",
    ),
    Param("SlurmConstraint", "slurm_constraint", "", "str", "sbatch constraint of the jobs"),
    Param("jobTimeLimit", "job_time_limit", 24, "num_or_none", "hours of a job (None: no limit)"),
    Param("queryInterval", "query_interval", 3, "num", "seconds between looks at the jobs"),
    Param(
        "language",
        "language",
        "python",
        "str",
        "'python' (task records; 'matlab' is taken for it) or 'bash' (commands)",
    ),
    Param(
        "GNUparallel",
        "gnu_parallel",
        False,
        "bool",
        "run several tasks at once in a job (in-job slots; GNU parallel in MATLAB)",
    ),
    Param("paraJobNum", "para_job_num", 1, "num", "tasks at once in a job"),
    Param(
        "masterParaFactor",
        "master_para_factor",
        1,
        "num",
        "fraction of the in-job slots on the master",
    ),
    Param(
        "ConfigFile",
        "config_file",
        "",
        "str",
        "cluster configuration file (JSON or .mat) whose values override these",
    ),
    Param(
        "localWorkerNum",
        "local_worker_num",
        1,
        "num",
        "without cluster, run up to this many tasks at once in a local process pool",
    ),
    Param(
        "taskMetrics",
        "task_metrics",
        False,
        "bool",
        "record each task's resources in <jobLogDir>/task_metrics/",
    ),
    Param(
        "taskRunner",
        "task_runner",
        True,
        "bool",
        "cluster jobs run several tasks per process and take waiting tasks of others",
    ),
    Param(
        "SlurmAlternatives",
        "slurm_alternatives",
        None,
        "any",
        "other places for the jobs (see jobs.alternatives)",
    ),
    Param(
        "gpuHeadroom",
        "gpu_headroom",
        None,
        "num_or_none",
        "with SlurmAlternatives: CPUs and GB per idle GPU that CPU jobs leave free",
    ),
    Param("nice", "nice", 0, "num", "lower the jobs' priority by this (sbatch --nice)"),
)

# MATLAB-only parameters and config keys: accepted and ignored
MATLAB_ONLY_PARAMS: tuple[Param, ...] = (
    Param("mccMode", "mcc_mode", False, "any", "MATLAB only (compiled jobs); ignored"),
    Param("MatlabLaunchStr", "matlab_launch_str", "", "any", "MATLAB only; ignored"),
    Param("MCRCacheRoot", "mcr_cache_root", "", "any", "MATLAB only; ignored"),
    Param("MCRParam", "mcr_param", "", "any", "MATLAB only; ignored"),
    Param("MCCMasterStr", "mcc_master_str", "", "any", "MATLAB only; ignored"),
)

_BY_PYTHON = {p.python: p for p in RUN_TASKS_PARAMS + MATLAB_ONLY_PARAMS}
_BY_MATLAB_LOWER = {p.matlab.lower(): p for p in RUN_TASKS_PARAMS + MATLAB_ONLY_PARAMS}


def run_tasks(
    input_paths: Sequence[str | Sequence[str]] | str,
    output_paths: Sequence[str] | str,
    tasks: Sequence[Any] | Any,
    *,
    final_out_fullpath: str = "",
    cluster_type: str = "slurm",
    parse_cluster: bool = True,
    master_compute: bool = True,
    job_log_dir: str = DEFAULT_JOB_LOG_DIR,
    tmp_dir: str = "",
    max_cpu_num: float = 24,
    min_cpu_num: float = 1,
    cpus_per_task: float = 1,
    mem_per_cpu: float = 20.9,
    mem_allocate: float | None = None,
    whole_node_job: bool = False,
    uuid: str = "",
    gpu_job: bool = False,
    gpu_num: float = 0,
    max_trial_num: float = 3,
    unit_wait_time: float = 3,
    max_job_num: float = math.inf,
    min_task_job_num: float = 1,
    task_batch_num: float = 1,
    min_batch_num: float = 1,
    para_batch_num: float = 1,
    run_extra_tasks: bool = False,
    python_launch_str: str = "",
    bash_launch_str: str = "",
    slurm_param: str = "",
    slurm_constraint: str = "",
    job_time_limit: float | None = 24,
    query_interval: float = 3,
    language: str = "python",
    gnu_parallel: bool = False,
    para_job_num: float = 1,
    master_para_factor: float = 1,
    config_file: str = "",
    local_worker_num: float = 1,
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
) -> list[bool]:
    """Runs tasks until their outputs exist; returns one done flag per task.

    input_paths  : per task a path, or a list of paths; a task whose input is missing is
                   not run and counts as done (MATLAB's behavior).
    output_paths : per task the path whose existence means the task is done.
    tasks        : per task a ``Task`` (or its record dict, a callable without arguments,
                   a functools.partial), or with ``language='bash'`` a command string.

    The other parameters are those of MATLAB's generic_computing_frameworks_wrapper in
    snake_case, with its defaults (RUN_TASKS_PARAMS lists them with their meaning).
    A config file's values override them, its keys matched without regard to case. The
    MATLAB-only parameters (mcc_mode, matlab_launch_str, mcr_cache_root, mcr_param,
    mcc_master_str) are accepted and ignored.

    Local mode (no cluster, or one task with master_compute) runs the tasks in order in
    this process, as MATLAB does: a ``<output name>.tmp`` marker in the output folder (or
    tmp_dir) leaves a task to another caller for unit_wait_time minutes; a task that raises
    is retried once at once, and a second exception propagates; a task that ends without
    its output is tried again in the next round, up to max_trial_num times; outputs are
    looked at once more at the end.
    """
    del mcc_mode, matlab_launch_str, mcr_cache_root, mcr_param, mcc_master_str
    given = dict(locals())
    settings = {p.python: given[p.python] for p in RUN_TASKS_PARAMS}
    inputs, outputs, task_list = _as_lists(input_paths, output_paths, tasks)
    _check_types(settings)
    if settings["config_file"]:
        settings = apply_config(settings, settings["config_file"])
        _check_types(settings)
    settings = _derived_settings(settings)
    language = settings["language"]

    records = [task_record(t, language) for t in task_list]
    outputs = [_strip_output(p) for p in outputs]
    call = Call(inputs, outputs, records, settings)
    if call.n == 0:
        return []

    # in cluster mode, a master in a Slurm job runs tasks itself only if they fit its memory
    if settings["parse_cluster"] and settings["master_compute"]:
        fits, task_mem, master_mem = _master_task_fits(settings)
        if not fits:
            print(
                f"The tasks ({task_mem:0.1f} GB) do not fit this job ({master_mem:0.1f} GB), "
                "so the master computes none of them."
            )
            settings["master_compute"] = False
    # one task with master compute runs here
    if call.n == 1 and settings["master_compute"]:
        settings["parse_cluster"] = False
    call.settings = settings

    cluster_type = settings["cluster_type"]
    if cluster_type == "parfor":
        from petakit5d.jobs import local  # J04

        return local.run_parfor(call)
    if cluster_type != "slurm":
        raise ValueError(f"cluster_type must be 'slurm' or 'parfor', not {cluster_type!r}")

    if settings["parse_cluster"]:
        if _slurm_available():
            from petakit5d.jobs import cluster  # J03

            return cluster.run_cluster(call)
        print(
            "Warning: A slurm-based computing cluster is not exist. Set parseCluster as false.",
            file=sys.stderr,
        )
        settings["parse_cluster"] = False

    # without cluster, a local process pool first (it leaves what it didn't finish here)
    if settings["local_worker_num"] > 1 and call.n > 1 and language == "python":
        try:
            from petakit5d.jobs import local  # J04
        except ImportError:
            local = None
        if local is not None:
            local.run_pool(call, worker_num=settings["local_worker_num"])
    return run_local(call)


def apply_config(settings: dict[str, Any], config_file: str) -> dict[str, Any]:
    """The settings with the values of a cluster config file: its keys match the MATLAB
    parameter names regardless of case (MemPerCPU sets mem_per_cpu); MATLAB-only and
    unknown keys are ignored."""
    conf = _read_config(config_file)
    out = dict(settings)
    for key, value in conf.items():
        param = _BY_MATLAB_LOWER.get(key.lower()) or _BY_PYTHON.get(key)
        if param is None or param in MATLAB_ONLY_PARAMS:
            continue
        # MATLAB's jsonencode writes an empty value as [] (jsondecode reads it back as [])
        if isinstance(value, list) and not value and param.kind in ("num_or_none", "any"):
            value = None
        out[param.python] = value
    return out


def run_local(call: Call) -> list[bool]:
    """The sequential local mode of MATLAB's framework, for a Call (also the cluster
    backend's fallback when there is no cluster). Returns the done flags."""
    s = call.settings
    n = call.n
    outputs = call.outputs
    final_out = s.get("final_out_fullpath", "")
    max_trial = s.get("max_trial_num", 3)
    unit_wait = s.get("unit_wait_time", 3)
    tmp_dir = s.get("tmp_dir", "")
    bs = call.task_batch_num

    if final_out and batch_file_exist([final_out])[0]:
        print(f"The final output file {final_out} already exists!\n")
        return [True] * n
    out_exist = batch_file_exist(outputs, None, True)
    if all(out_exist):
        print(f"All output files ({n} / {n}) already exist(s)!\n")
        return [True] * n
    done = list(out_exist)
    in_exist = input_exist(call, out_exist)

    print(f"Task number : {n}")
    trials = [0] * n
    loop = 0
    n_done = 0
    ts = time.monotonic()
    fsnames = [""] * n

    def recheck() -> None:
        undone = [i for i in range(n) if not done[i]]
        for i, ok in zip(
            undone, batch_file_exist([outputs[i] for i in undone], None, True), strict=True
        ):
            out_exist[i] = ok
            done[i] = ok

    while not all(done[i] or trials[i] >= max_trial for i in range(n)):
        if loop > 0:
            recheck()
        ran = False
        # (MATLAB loops over nF batches in local mode; those past the last are empty)
        for b in range(1, n + 1):
            fs = range((b - 1) * bs, min(b * bs, n))
            if b % STATUS_CHECK_BATCHES == 0 and b + STATUS_CHECK_BATCHES * 0.5 < n:
                recheck()
            if all(done[f] or trials[f] >= max_trial for f in fs):
                continue
            for f in fs:
                inp = call.inputs[f]
                if isinstance(inp, list):
                    if not all(exist(x) for x in inp):
                        done[f] = True
                elif not in_exist[f]:
                    done[f] = True
                out_dir, fsname, ext = fileparts(outputs[f])
                if not fsname and ext:
                    fsname = ext
                fsnames[f] = fsname
                cur_tmp_dir = tmp_dir or out_dir
                tmp = f"{cur_tmp_dir}/{fsname}.tmp"
                if out_exist[f]:
                    done[f] = True
                    _remove(tmp)
            if all(done[f] or trials[f] >= max_trial for f in fs):
                continue

            # the marker of the batch is that of its last task
            if os.path.exists(tmp):
                age = _marker_age(tmp)
                if age is not None and age / 60 < unit_wait:
                    continue
                _touch(tmp)
            elif os.path.isdir(cur_tmp_dir):
                _touch(tmp)

            ran = True
            t1 = run_batch_in_process(call, fs, [fsnames[f] for f in fs])
            for f in fs:
                trials[f] += 1
            if os.path.exists(outputs[fs[-1]]) and os.path.exists(tmp):
                _remove(tmp)
            print(f"Done! Elapsed time is {t1:f} seconds.")

        if n_done < sum(done):
            n_done = sum(done)
            if n_done == n:
                print()
            print(
                f"Time {time.monotonic() - ts:0.2f} s: {n_done} / {n} "
                f"({n_done / n * 100:0.2f}%) are finished!"
            )

        # batch_file_exist sees a final output another node just made
        if final_out and batch_file_exist([final_out])[0]:
            done = [True] * n
            print(
                f"Time {time.monotonic() - ts:0.2f} s: The final output file {final_out} "
                "already exists!"
            )
            break
        loop += 1
        if not ran and not all(done[i] or trials[i] >= max_trial for i in range(n)):
            time.sleep(LOCAL_IDLE_WAIT)

    # a batch run in its last try ran after this round looked at its outputs: look once more
    if not all(done):
        recheck()
    if all(done):
        print(f"Time {time.monotonic() - ts:0.2f} s: All output files ({n} / {n}) are finished!\n")
    return done


def run_batch_in_process(
    call: Call, indices: Sequence[int], names: Sequence[str] = (), mode: str = "local"
) -> float:
    """Runs the tasks at ``indices`` (0-based) in this process as MATLAB's master does and
    returns the seconds they took. Python tasks run in order; if one raises, the batch runs
    once more from its first task, and a second exception propagates. Bash tasks run as one
    shell command (joined by ';'), their exit status ignored. ``mode`` ('local' or 'master')
    goes to the task metrics."""
    recs = [call.records[i] for i in indices]
    text = ";".join(describe(r) for r in recs)
    print(f"\nProcess {', '.join(names)} with function {text}... ")
    tm = _metrics_start(call, indices, mode)
    t0 = time.perf_counter()
    if call.language == "bash":
        status = run_shell(text)
        t1 = time.perf_counter() - t0
        _metrics_stop(tm, status)
        return t1
    try:
        for rec in recs:
            run_record(rec)
    except Exception:
        traceback.print_exc(file=sys.stdout)
        t0 = time.perf_counter()
        try:
            for rec in recs:
                run_record(rec)
        except Exception as err:
            _metrics_stop(tm, err)
            raise
    t1 = time.perf_counter() - t0
    _metrics_stop(tm, None)
    return t1


_metrics_warned = False


def _metrics_start(call: Call, indices: Sequence[int], mode: str) -> Any:
    # task metrics (J10): metrics.task_start(call, indices, mode) -> record, and
    # metrics.task_stop(record, outcome) with None, the exception or the shell's status
    global _metrics_warned
    if not call.settings.get("task_metrics"):
        return None
    try:
        from petakit5d.jobs import metrics
    except ImportError:
        if not _metrics_warned:
            print(
                "Warning: task metrics are not available (petakit5d.jobs.metrics)", file=sys.stderr
            )
            _metrics_warned = True
        return None
    return (metrics, metrics.task_start(call, list(indices), mode))


def _metrics_stop(tm: Any, outcome: Any) -> None:
    if tm is not None:
        tm[0].task_stop(tm[1], outcome)


def run_shell(command: str) -> int:
    """Runs a command in bash (cmd.exe on Windows) with this process's output; returns its
    exit status."""
    sys.stdout.flush()
    sys.stderr.flush()
    if os.name == "nt":
        return subprocess.run(command, shell=True).returncode
    import shutil

    bash = shutil.which("bash")
    return subprocess.run(command, shell=True, executable=bash).returncode


def input_exist(call: Call, out_exist: Sequence[bool]) -> list[bool]:
    """Whether each task's input exists, looked at for the tasks without output (the others
    count as having it), as MATLAB does; for a list of inputs its first path."""
    flags = list(out_exist)
    todo = [i for i in range(call.n) if not out_exist[i]]
    if not todo:
        return flags
    paths = [call.inputs[i] if isinstance(call.inputs[i], str) else call.inputs[i][0] for i in todo]
    for i, ok in zip(todo, batch_file_exist(paths, None, True), strict=True):
        flags[i] = ok
    return flags


def _as_lists(input_paths, output_paths, tasks):
    def paths(x, what):
        if isinstance(x, (str, os.PathLike)):
            return [os.fspath(x)]
        out = []
        for p in x:
            if isinstance(p, (str, os.PathLike)):
                out.append(os.fspath(p))
            elif what == "input" and isinstance(p, (list, tuple)):
                out.append([os.fspath(q) for q in p])
            else:
                raise TypeError(
                    f"an {what} path must be a str or os.PathLike, not {type(p).__name__}"
                )
        return out

    if isinstance(tasks, (str, Mapping, Task)) or callable(tasks):
        tasks = [tasks]
    return paths(input_paths, "input"), paths(output_paths, "output"), list(tasks)


def _strip_output(path: str) -> str:
    # MATLAB: strip(..., 'right', filesep), then strip(..., 'right', '/')
    return path.rstrip(os.sep).rstrip("/")


def _check_types(settings: dict[str, Any]) -> None:
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
            matlab = _BY_PYTHON[name].matlab
            raise TypeError(
                f"The value of '{matlab}' ({name}) is invalid: expected "
                f"{kind.replace('_or_none', ' or None')}, got {value!r}"
            )


def _derived_settings(settings: dict[str, Any]) -> dict[str, Any]:
    s = dict(settings)
    s["parse_cluster"] = bool(s["parse_cluster"])
    s["master_compute"] = bool(s["master_compute"])
    # an empty job log folder (as "jobLogDir": "" in a config file) is the default one
    if not s["job_log_dir"]:
        s["job_log_dir"] = DEFAULT_JOB_LOG_DIR
    lang = s["language"].lower()
    s["language"] = "python" if lang == "matlab" else lang
    if s["language"] not in ("python", "bash"):
        raise ValueError(f"language must be 'python' or 'bash', not {settings['language']!r}")
    if s["nice"] > 0:
        s["slurm_param"] = f"{s['slurm_param']} --nice={_round_half_away(s['nice'])}"
    s["task_batch_num"] = int(s["task_batch_num"])
    if s["task_batch_num"] < 1:
        raise ValueError("task_batch_num must be at least 1")
    if not s["uuid"]:
        s["uuid"] = _get_uuid()
    return s


def _remove(path: str) -> None:
    if os.path.isfile(path) or os.path.islink(path):
        try:
            os.remove(path)
        except OSError:
            pass


def _touch(path: str) -> None:
    # MATLAB's fclose(fopen(path, 'w')): makes the file empty and new
    try:
        with open(path, "w"):
            pass
    except OSError:
        pass


def _marker_age(path: str) -> float | None:
    # seconds since the marker's modification time, whole seconds as MATLAB's dir().date;
    # None if it went away
    try:
        mtime = os.stat(path).st_mtime
    except OSError:
        return None
    return time.time() - math.floor(mtime)


# --- STAND-INS (replaced when their owners' modules land) --------------------------------


def _read_config(config_file: str) -> dict[str, Any]:
    """STAND-IN for J01's cluster config reader (jobs.config): JSON, or .mat through F07."""
    if not os.path.isfile(config_file):
        raise FileNotFoundError(
            f"The cluster configuration file {config_file} does not exist, please check and "
            "provide the right one!"
        )
    print(f"Set up cluster configuration according to {config_file}...")
    if config_file.lower().endswith(".mat"):
        from petakit5d.io import matfile  # F07

        conf = dict(matfile.read_mat(config_file))
    else:
        import json

        with open(config_file, encoding="utf-8") as fh:
            conf = json.load(fh)
    if not isinstance(conf, dict):
        raise ValueError(f"The cluster configuration file {config_file} is not a JSON object")
    return conf


def _master_task_fits(settings: dict[str, Any]) -> tuple[bool, float, float]:
    """STAND-IN for J08's master_task_fits: outside a Slurm job always; in one, the task's
    memory against 1.05 x the job's memory from Slurm's variables."""
    mem_allocate = settings["mem_allocate"]
    max_cpu, min_cpu = settings["max_cpu_num"], settings["min_cpu_num"]
    mem_per_cpu = settings["mem_per_cpu"]
    if mem_allocate is None:
        task_mem = min(max_cpu, max(min_cpu, settings["cpus_per_task"])) * mem_per_cpu
    else:
        task_mem = min(mem_allocate, max_cpu * mem_per_cpu)
    if not os.environ.get("SLURM_JOB_ID"):
        return True, task_mem, math.inf
    env = os.environ
    if env.get("SLURM_MEM_PER_NODE"):
        master_mem = float(env["SLURM_MEM_PER_NODE"]) / 1024
    elif env.get("SLURM_MEM_PER_CPU") and env.get("SLURM_CPUS_ON_NODE"):
        master_mem = float(env["SLURM_MEM_PER_CPU"]) * float(env["SLURM_CPUS_ON_NODE"]) / 1024
    else:
        return False, task_mem, math.nan
    return task_mem <= master_mem * 1.05, task_mem, master_mem


def _slurm_available() -> bool:
    """STAND-IN for J01's cluster check (checkSlurmCluster): sinfo runs and exits with 0."""
    import shutil

    sinfo = shutil.which("sinfo")
    if not sinfo:
        return False
    try:
        return subprocess.run([sinfo], capture_output=True, timeout=60).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def _round_half_away(x: float) -> int:
    """STAND-IN for F02's MATLAB round."""
    return int(math.copysign(math.floor(abs(x) + 0.5), x))


def _get_uuid() -> str:
    """STAND-IN for F02's get_uuid."""
    import uuid

    return str(uuid.uuid4())
