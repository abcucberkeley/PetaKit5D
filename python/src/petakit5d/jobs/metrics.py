"""Task metrics of the job framework (MATLAB's task_metrics, task_metrics_run and
read_task_metrics): with ``task_metrics=True`` every task run leaves a record of what it used in
``<jobLogDir>/task_metrics/<run>/``, and read_task_metrics collects the records into tables.

Off by default, and then nothing here is imported. Errors never stop or fail a task: a metrics
error is one warning per process (``task_metrics:failed``), and the task then has no record.

Writer (the framework's hooks):
  ``tm = task_start(call, indices, mode)``  the state before a batch of tasks (0-based indices
                                            of ``call``; mode 'local', 'master', 'pool' or
                                            'slurm_python'); None without metrics
  ``task_stop(tm, outcome)``                appends the record; outcome None (finished), the
                                            exception, or a shell's exit status
  ``run_start(call, backend)``              makes the run folder and its run_start record
  ``run_end(call, done)``                   the run_end record (elapsed time, tasks done)
  ``session_dir(job_log_dir)``, ``record(run_dir, 'submit', info)``  single jobs

Files, as MATLAB writes them (all single-writer, so safe on NFS):
  ``task_<node>_<pid>_<process start ticks>.ndjson``  one JSON line per task, one file per process
  ``run.ndjson``                                      the master's run_start, run_end, submit lines

Where: ``<jobLogDir>/task_metrics/<yyyyMMdd_HHmmss>_<function>_<id8>/`` when jobLogDir is an
existing absolute folder, else ``<first input folder>/job_logs/task_metrics/...``; single jobs
one folder per process and job log folder (``..._single_<id8>``).

Record format (version 1, MATLAB's; MATLAB's read_task_metrics reads our folders and we read
MATLAB's). A task line is one JSON object, keys in this order:

  kind "task"; runId (the run folder's name); func (the first task's function,
  "package.module:name", or the bash command, up to its first "("); mode ("local", "master",
  "pool", "slurm_python"; MATLAB also writes "slurm_matlab", "slurm_mcc", "mcc"); status ("ok",
  "error": an exception, "failed": a nonzero exit status); error ("<identifier or exception
  class>: <message>", at most 1000 characters); node; jobId, arrayJobId, arrayTaskId (Slurm's
  variables); gpu (name); cudaDevices; parent ("<runId>/<task>" of the task this framework call
  ran in, or ""); version 1;
  task (1-based index of the first task of the batch); batchSize; exitCode; pid; slot, seq
  (PARALLEL_JOBSLOT, PARALLEL_SEQ); start, end (seconds since 1970, UTC); wallSec; cpuSec (this
  process, all threads); userSec, sysSec; childCpuSec (waited-for children); rssStartGB,
  rssEndGB, rssPeakGB (the task's peak: VmHWM is reset at its start), rssPeakBeforeGB;
  readBytes, writeBytes (storage); rchar, wchar (all read/write calls); minorFaults,
  majorFaults; threads (the thread count our kernels use by default; MATLAB's
  maxNumCompThreads); cpus (allowed CPUs); slurmCpus (SLURM_CPUS_PER_TASK); osThreads;
  gpuMemGB, gpuFreeGB; processStart (seconds since 1970); memTotalGB (node); metricsSec (the
  metrics' own time). Numbers as %.16g, null where unknown (memory, I/O and faults are from
  /proc and null off Linux; user, system and child CPU come from os.times() there).

run.ndjson lines (MATLAB's jsonencode layout): kind ("run_start", "run_end", "submit"),
version, runId, time, node, pid, jobId, parent, then run_start: the parameters of run_tasks by
MATLAB name, inputFullpaths, outputFullpaths, functionStrs, backend ("python"; MATLAB
"matlab", "mcc"), parseCluster, masterCompute, cpusPerTask (effective), mccMode (false),
taskBatchNum, paraJobNum, poolWorkerNum, elapsedSec; run_end: nTasks, nDone, elapsedSec;
submit: task, jobId, functionStr, cpusPerTask, memAllocate, memPerCPU, mccMode, SlurmParam.
"""

from __future__ import annotations

import json
import math
import os
import re
import struct
import sys
import time
import warnings
from collections.abc import Mapping, Sequence
from typing import Any

__all__ = [
    "TaskMetricsWarning",
    "process_info",
    "read_task_metrics",
    "record",
    "run_dir",
    "run_end",
    "run_start",
    "session_dir",
    "task_start",
    "task_stop",
    "worker_task_start",
]

# MATLAB's NaN (sign bit set), for missing values of the tables
NAN: float = struct.unpack("<d", struct.pack("<Q", 0xFFF8000000000000))[0]
# the modes of tasks that ran in a process of their own (start-up costs)
OWN_PROCESS_MODES = ("slurm_matlab", "slurm_mcc", "mcc", "slurm_python")
MATLAB_OWN_PROCESS_MODES = ("slurm_matlab", "slurm_mcc", "mcc")
# the numbers of a task record, in MATLAB's order
TASK_NUM_FIELDS = (
    "task", "batchSize", "exitCode", "pid", "slot", "seq", "start", "end", "wallSec", "cpuSec",
    "userSec", "sysSec", "childCpuSec", "rssStartGB", "rssEndGB", "rssPeakGB", "rssPeakBeforeGB",
    "readBytes", "writeBytes", "rchar", "wchar", "minorFaults", "majorFaults", "threads", "cpus",
    "slurmCpus", "osThreads", "gpuMemGB", "gpuFreeGB", "processStart", "memTotalGB", "metricsSec",
)  # fmt: skip
TASK_STR_FIELDS = (
    "runId", "func", "mode", "status", "error", "node", "jobId", "arrayJobId", "arrayTaskId",
    "gpu", "cudaDevices", "parent",
)  # fmt: skip
TASK_COLUMNS = (
    "runId", "task", "func", "mode", "status", "exitCode", "error", "node", "pid", "jobId",
    "arrayJobId", "arrayTaskId", "slot", "seq", "start", "end", "wallSec", "cpuSec", "userSec",
    "sysSec", "childCpuSec", "rssStartGB", "rssEndGB", "rssPeakGB", "rssPeakBeforeGB",
    "readBytes", "writeBytes", "rchar", "wchar", "minorFaults", "majorFaults", "threads", "cpus",
    "slurmCpus", "osThreads", "gpu", "gpuMemGB", "gpuFreeGB", "cudaDevices", "processStart",
    "memTotalGB", "batchSize", "parent", "metricsSec",
)  # fmt: skip
RUN_COLUMNS = (
    "runId", "kind", "func", "backend", "node", "start", "end", "elapsedSec", "nTasks", "nDone",
    "nRecords", "parseCluster", "mccMode", "cpusPerTask", "memAllocate", "taskBatchNum",
    "paraJobNum", "GPUJob", "parent", "folder",
)  # fmt: skip
SACCT_COLUMNS = (
    "jobSubmit", "jobStart", "jobEnd", "jobState", "queueSec", "jobHeadSec", "jobTailSec",
)  # fmt: skip


class TaskMetricsWarning(UserWarning):
    """A metrics problem (the task itself is not affected)."""

    def __init__(self, message: str, identifier: str = "task_metrics:failed") -> None:
        super().__init__(message)
        self.identifier = identifier


# --- the process ---------------------------------------------------------------------------

_LINUX = sys.platform.startswith("linux")
_GB_KB = 1048576.0  # kB per GB
_GB_B = 1073741824.0  # bytes per GB
_process: dict[str, Any] | None = None
_current = ""  # the task running in this process ("<runId>/<task>"), the parent of nested calls
_warned = False
_sessions: dict[str, str] = {}  # single jobs: job log folder -> metrics folder
_run_t0: dict[str, float] = {}  # run folder -> perf_counter at run_start


def _env(name: str) -> str:
    return os.environ.get(name, "")


def _read_text(path: str, n: int = -1) -> str:
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read(n)
    except OSError:
        return ""


def _field(s: str, name: str) -> float:
    # a number after name in /proc/self/status, NaN if missing
    k = s.find(name)
    if k < 0:
        return NAN
    m = re.match(r"\s*([-+]?[0-9.]+(?:[eE][-+]?[0-9]+)?)", s[k + len(name) :])
    return float(m.group(1)) if m else NAN


def _snapshot() -> list[float]:
    """From /proc/self: resident memory and its peak (kB), threads, allowed CPUs, minor and
    major page faults, user, system and waited-for children's CPU (clock ticks), the process
    start (ticks after boot), rchar, wchar, read_bytes, write_bytes; NaN elsewhere."""
    v = [NAN] * 15
    if not _LINUX:
        return v
    s = _read_text("/proc/self/status")
    v[0] = _field(s, "VmRSS:")
    v[1] = _field(s, "VmHWM:")
    v[2] = _field(s, "Threads:")
    k = s.find("Cpus_allowed:")
    if k >= 0:
        # the CPUs set in the hexadecimal mask
        words = s[k + 13 :].split(None, 1)
        try:
            v[3] = float(bin(int(words[0].replace(",", ""), 16)).count("1"))
        except (IndexError, ValueError):
            pass
    s = _read_text("/proc/self/stat")
    k = s.rfind(")")
    if k >= 0:
        x = s[k + 2 :].split()  # from the state (field 3) on
        if len(x) >= 20:
            try:
                v[4:11] = [float(x[i]) for i in (7, 9, 11, 12, 13, 14, 19)]
            except ValueError:
                pass
    io = dict(
        line.split(":", 1) for line in _read_text("/proc/self/io").splitlines() if ":" in line
    )
    try:
        v[11:15] = [float(io[n]) for n in ("rchar", "wchar", "read_bytes", "write_bytes")]
    except (KeyError, ValueError):
        pass
    return v


def _cpu_times() -> tuple[float, ...]:
    return tuple(os.times())


def _hz() -> float:
    try:
        return float(os.sysconf("SC_CLK_TCK"))
    except (AttributeError, ValueError, OSError):
        return 100.0


def process_info() -> dict[str, Any]:
    """This process's node, id, start (seconds since 1970), node memory, record file name and
    the task running in it (``current``); found once per process."""
    global _process
    p = _process
    if p is None or p["pid"] != os.getpid():
        p = _new_process_info()
        _process = p
    return dict(p, current=_current)


def _new_process_info() -> dict[str, Any]:
    p: dict[str, Any] = {"linux": _LINUX, "pid": os.getpid(), "hz": _hz()}
    host = _env("SLURMD_NODENAME")
    if not host and _LINUX:
        host = _read_text("/proc/sys/kernel/hostname").strip()
    host = host or _env("COMPUTERNAME") or _env("HOSTNAME")
    if not host:
        import socket

        host = socket.gethostname()
    p["host"] = host
    now = time.time()
    p["processStart"] = NAN
    p["memTotalGB"] = NAN
    token = f"{round(now * 1000):d}"
    if _LINUX:
        ticks = _snapshot()[10]
        up = _read_text("/proc/uptime").split()
        if not math.isnan(ticks) and up:
            p["processStart"] = now - (float(up[0]) - ticks / p["hz"])
            token = f"{int(ticks):d}"
        p["memTotalGB"] = _field(_read_text("/proc/meminfo", 256), "MemTotal:") / _GB_KB
    try:
        from petakit5d.core.threads import available_cpus

        p["threads"] = float(available_cpus())
    except Exception:
        p["threads"] = float(os.cpu_count() or 1)
    p["file"] = f"task_{_safe_name(host)}_{p['pid']:d}_{token}.ndjson"
    return p


def _safe_name(s: str) -> str:
    # for file names: letters, digits, '_', '-' and '.'
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", s)


def _warn_once(err: BaseException) -> None:
    global _warned
    if not _warned:
        _warned = True
        warnings.warn(TaskMetricsWarning(f"Task metrics are incomplete: {err}"), stacklevel=3)


# --- task records ----------------------------------------------------------------------------


def task_start(call: Any, indices: Sequence[int], mode: str) -> dict[str, Any] | None:
    """The state of this process before the tasks at 0-based ``indices`` of ``call`` (a
    ``tasks.Call`` or ``tasks.CallFolder``) run as one batch; None if the call has no metrics
    folder. ``mode``: 'local', 'master', 'pool' or 'slurm_python'. A call with
    ``task_metrics=True`` but no folder yet gets one here in the local modes."""
    try:
        t0 = time.perf_counter()
        settings = _settings(call)
        folder = settings.get("task_metrics")
        if folder is True and mode in ("local", "master"):
            folder = run_start(call)
        if not isinstance(folder, str) or not folder:
            return None
        indices = list(indices)
        text = _first_text(call, indices)
        return _begin(folder, indices, text, mode, bool(settings.get("gpu_job")), t0)
    except Exception as err:
        _warn_once(err)
        return None


def worker_task_start(
    record: Mapping[str, Any],
    index: int,
    metrics_dir: str,
    worker: Any = None,
    gpu_job: bool = False,
) -> dict[str, Any] | None:
    """task_start for a process that has a task record but not its Call (a local pool
    worker): ``record`` the task record, ``index`` its 0-based index in the call,
    ``metrics_dir`` the run's metrics folder (``call.settings['task_metrics']``) or the call
    folder (its call.json names the metrics folder), ``worker`` the pool worker
    (``local.worker_info()``: its id is the record's slot, its threads the record's threads,
    a GPU of its own makes it a GPU task). Mode 'pool'. None without a metrics folder. End
    it with ``task_stop(tm, outcome)``."""
    try:
        t0 = time.perf_counter()
        folder = _metrics_folder(metrics_dir)
        if not folder:
            return None
        text = record.get("bash", "") if "bash" in record else record.get("func", "")
        gpu = bool(gpu_job) or bool(getattr(worker, "gpu", None))
        tm = _begin(folder, [int(index)], text, "pool", gpu, t0)
        if worker is not None:
            tm["slot"] = float(worker.id)
            tm["threads"] = float(worker.threads)
        return tm
    except Exception as err:
        _warn_once(err)
        return None


_call_folders: dict[str, str] = {}  # call folder -> its metrics folder


def _metrics_folder(path: str) -> str:
    # a metrics folder as it is, or the one a call folder's call.json names
    if not path:
        return ""
    if path in _call_folders:
        return _call_folders[path]
    fn = os.path.join(path, "call.json")
    if not os.path.isfile(fn):
        return path
    with open(fn, encoding="utf-8") as fh:
        folder = json.load(fh).get("settings", {}).get("task_metrics")
    folder = folder if isinstance(folder, str) else ""
    _call_folders[path] = folder
    return folder


def _begin(
    folder: str, indices: list[int], text: str, mode: str, gpu: bool, t0: float
) -> dict[str, Any]:
    # the state before a task, and this process's current task
    global _current
    v0 = _snapshot()
    p = process_info()
    d = folder.replace("\\", "/").rstrip("/")
    tm: dict[str, Any] = {
        "runId": d[d.rfind("/") + 1 :],
        "runDir": folder,
        "task": [i + 1 for i in indices],
        "func": _func_name(text),
        "mode": mode,
        "gpu": gpu,
        "threads": p["threads"],
        "slot": None,
        "parent": _current,
        "proc0": v0,
    }
    if p["linux"]:
        # the peak resident memory from here on is the task's
        try:
            with open("/proc/self/clear_refs", "w") as fh:
                fh.write("5")
        except OSError:
            pass
    _current = f"{tm['runId']}/{tm['task'][0] if tm['task'] else 'NaN'}"
    tm["metricsSec"] = time.perf_counter() - t0
    tm["cpu0"] = time.process_time()
    tm["times0"] = _cpu_times()
    tm["start"] = time.time()
    tm["wall0"] = time.perf_counter()
    return tm


def task_stop(tm: Mapping[str, Any] | None, outcome: Any = None) -> None:
    """Appends the record of the task(s) ``task_start`` began to this process's file in the
    run folder. ``outcome``: None or 0 (finished), an exit status, or the exception."""
    global _current
    if tm is None:
        return
    try:
        wall = time.perf_counter() - tm["wall0"]
        cpu = time.process_time() - tm["cpu0"]
        t0 = time.perf_counter()
        _current = tm["parent"]
        times1 = _cpu_times()
        p = process_info()
        _append_line(
            os.path.join(tm["runDir"], p["file"]), _task_line(tm, outcome, wall, cpu, p, times1, t0)
        )
    except Exception as err:
        _warn_once(err)


def _task_line(tm, outcome, wall, cpu, p, times1, t0) -> str:
    v0 = tm["proc0"]
    v1 = _snapshot()
    hz = p["hz"]
    status, exit_code, err_str = "ok", 0.0, ""
    if isinstance(outcome, BaseException):
        status, exit_code = "error", 1.0
        ident = getattr(outcome, "identifier", "") or type(outcome).__qualname__
        err_str = f"{ident}: {outcome}"[:1000]
    elif outcome is not None and outcome != 0:
        status, exit_code = "failed", float(outcome)
    if p["linux"]:
        user = (v1[6] - v0[6]) / hz
        sys_ = (v1[7] - v0[7]) / hz
        child = (v1[8] + v1[9] - v0[8] - v0[9]) / hz
    else:
        t0s = tm["times0"]
        user, sys_ = times1[0] - t0s[0], times1[1] - t0s[1]
        child = times1[2] + times1[3] - t0s[2] - t0s[3]
    gpu_name, gpu_mem = "", (NAN, NAN)
    cuda = _env("CUDA_VISIBLE_DEVICES")
    if tm["gpu"] or cuda:
        gpu_name, gpu_mem = _gpu_info()
    task = tm["task"]
    start = tm["start"]
    slot = _env_num("PARALLEL_JOBSLOT") if tm.get("slot") is None else tm["slot"]
    nums = (
        task[0] if task else NAN, len(task), exit_code, p["pid"], slot,
        _env_num("PARALLEL_SEQ"), start, start + wall, wall, cpu, user, sys_, child,
        v0[0] / _GB_KB, v1[0] / _GB_KB, v1[1] / _GB_KB, v0[1] / _GB_KB, v1[13] - v0[13],
        v1[14] - v0[14], v1[11] - v0[11], v1[12] - v0[12], v1[4] - v0[4], v1[5] - v0[5],
        tm["threads"], v0[3], _env_num("SLURM_CPUS_PER_TASK"), v1[2], gpu_mem[0], gpu_mem[1],
        p["processStart"], p["memTotalGB"], tm["metricsSec"] + time.perf_counter() - t0,
    )  # fmt: skip
    head = (
        f'{{"kind":"task","runId":"{_json_str(tm["runId"])}","func":"{_json_str(tm["func"])}",'
        f'"mode":"{_json_str(tm["mode"])}","status":"{status}","error":"{_json_str(err_str)}",'
        f'"node":"{_json_str(p["host"])}","jobId":"{_json_str(_env("SLURM_JOB_ID"))}",'
        f'"arrayJobId":"{_json_str(_env("SLURM_ARRAY_JOB_ID"))}",'
        f'"arrayTaskId":"{_json_str(_env("SLURM_ARRAY_TASK_ID"))}",'
        f'"gpu":"{_json_str(gpu_name)}","cudaDevices":"{_json_str(cuda)}",'
        f'"parent":"{_json_str(tm["parent"])}","version":1'
    )
    body = ",".join(f'"{k}":{_num(x)}' for k, x in zip(TASK_NUM_FIELDS, nums, strict=True))
    return f"{head},{body}}}"


def _num(x: float) -> str:
    # MATLAB's %.16g, null for NaN and Inf
    x = float(x)
    return f"{x:.16g}" if math.isfinite(x) else "null"


def _env_num(name: str) -> float:
    # MATLAB's sscanf(getenv(name), '%f', 1), NaN if none
    m = re.match(r"\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)", _env(name))
    return float(m.group(1)) if m else NAN


def _json_str(s: str) -> str:
    # a string for JSON, escaped as MATLAB's task_metrics escapes it
    if any(c in s for c in '\\"') or any(ord(c) < 32 for c in s):
        s = s.replace("\\", "\\\\").replace('"', '\\"')
        s = s.replace("\n", "\\n").replace("\r", "\\r").replace("\t", "\\t")
        s = "".join(" " if ord(c) < 32 else c for c in s)
    return s


def _gpu_info() -> tuple[str, tuple[float, float]]:
    # the name, memory and free memory (GB) of the GPU this process uses, without making a
    # CUDA context (only if CuPy is imported and has one)
    cp = sys.modules.get("cupy")
    if cp is None:
        return "", (NAN, NAN)
    try:
        if not cp.cuda.driver.ctxGetCurrent():
            return "", (NAN, NAN)
        dev = cp.cuda.runtime.getDevice()
        name = cp.cuda.runtime.getDeviceProperties(dev)["name"]
        free, total = cp.cuda.runtime.memGetInfo()
        name = name.decode() if isinstance(name, bytes) else str(name)
        return name, (total / _GB_B, free / _GB_B)
    except Exception:
        return "", (NAN, NAN)


def _append_line(path: str, line: str) -> None:
    # appends a line to a file only this process writes
    try:
        fh = open(path, "a", encoding="utf-8", newline="\n")
    except FileNotFoundError:
        # the folder is gone (or not yet seen on this node)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fh = open(path, "a", encoding="utf-8", newline="\n")
    with fh:
        fh.write(line + "\n")


def _settings(call: Any) -> Mapping[str, Any]:
    s = getattr(call, "settings", None)
    if s is None:
        s = getattr(call, "meta", {}).get("settings", {})
    return s


def _first_text(call: Any, indices: list[int]) -> str:
    # the first task's text: a Python task's function, a bash command
    if not indices:
        return ""
    recs = getattr(call, "records", None)
    rec = recs[indices[0]] if recs is not None else call.lines(indices[:1])[0]["task"]
    return rec.get("bash", "") if "bash" in rec else rec.get("func", "")


def _func_name(s: str) -> str:
    # the function a task string calls: the text before its first '(', trimmed
    k = s.find("(")
    return (s if k < 0 else s[:k]).strip()


# --- runs ------------------------------------------------------------------------------------


def run_dir(job_log_dir: str, input_paths: Any, func_strs: Any) -> str:
    """A new metrics folder for a framework call,
    ``<base>/task_metrics/<date>_<time>_<function>_<id>/``: base is ``job_log_dir`` if it is
    an existing absolute folder, else ``<first input folder>/job_logs``. '' if it can't be made
    (with a warning)."""
    try:
        first = input_paths
        while isinstance(first, (list, tuple)) and first:
            first = first[0]
        data_path = os.path.dirname(os.fspath(first) if first else "").replace("\\", "/")
        if not _is_absolute(data_path):
            data_path = os.getcwd().replace("\\", "/") + "/" + data_path
        if (
            isinstance(job_log_dir, str)
            and _is_absolute(job_log_dir)
            and os.path.isdir(job_log_dir)
        ):
            base = job_log_dir.replace("\\", "/")
        else:
            base = data_path + "/job_logs"
        func = func_strs
        while isinstance(func, (list, tuple)) and func:
            func = func[0]
        name = _safe_name(_func_name(func if isinstance(func, str) else ""))
        return _new_dir(base, name[:40])
    except Exception as err:
        _warn_run(err, "run_dir")
        return ""


def _token() -> str:
    import uuid

    return uuid.uuid4().hex[:8]


def _new_dir(base: str, name: str) -> str:
    run_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{name}_{_token()}"
    folder = f"{base.rstrip('/')}/task_metrics/{run_id}/"
    os.makedirs(folder)
    return folder


def _is_absolute(p: str) -> bool:
    return bool(p) and (p[0] in "/\\" or (len(p) > 1 and p[1] == ":"))


def _warn_run(err: BaseException, action: str) -> None:
    warnings.warn(
        TaskMetricsWarning(
            f"Task metrics ({action}) are incomplete: {err}", "task_metrics_run:failed"
        ),
        stacklevel=3,
    )


def run_start(call: Any, backend: str = "python", **effective: Any) -> str:
    """Makes the metrics folder of a framework call and its run_start record (the settings by
    MATLAB name, the tasks, and ``effective`` values such as cpusPerTask, taskBatchNum,
    paraJobNum, poolWorkerNum), and sets ``call.settings['task_metrics']`` to the folder, so
    that the backends and jobs record into it. Returns the folder ('' if it can't be made)."""
    try:
        from petakit5d.jobs.tasks import describe

        settings = call.settings
        folder = settings.get("task_metrics")
        if isinstance(folder, str) and folder:
            return folder
        t0 = time.perf_counter()
        func_strs = [describe(r) for r in call.records]
        folder = run_dir(settings.get("job_log_dir", ""), call.inputs, func_strs)
        if not folder:
            settings["task_metrics"] = False
            return ""
        _run_t0[folder] = t0
        settings["task_metrics"] = folder
        info = _matlab_settings(settings)
        info["taskMetrics"] = True
        info.update(
            inputFullpaths=list(call.inputs),
            outputFullpaths=list(call.outputs),
            functionStrs=func_strs,
            backend=backend,
            parseCluster=bool(settings.get("parse_cluster")),
            masterCompute=bool(settings.get("master_compute")),
            cpusPerTask=_effective_cpus(settings),
            mccMode=False,
            taskBatchNum=settings.get("task_batch_num", 1),
            paraJobNum=settings.get("para_job_num", 1),
            poolWorkerNum=0,
        )
        info.update(effective)
        info["elapsedSec"] = time.perf_counter() - t0
        record(folder, "run_start", info)
        return folder
    except Exception as err:
        _warn_run(err, "run_start")
        return ""


def run_end(call: Any, done: Sequence[bool]) -> None:
    """The run_end record of a call that has a metrics folder: tasks, tasks done, elapsed
    seconds since run_start."""
    folder = _settings(call).get("task_metrics")
    if not isinstance(folder, str) or not folder:
        return
    t0 = _run_t0.pop(folder, None)
    elapsed = NAN if t0 is None else time.perf_counter() - t0
    record(
        folder,
        "run_end",
        {"nTasks": len(done), "nDone": sum(bool(x) for x in done), "elapsedSec": elapsed},
    )


def _matlab_settings(settings: Mapping[str, Any]) -> dict[str, Any]:
    # the settings by MATLAB name, as MATLAB's run_start has ip.Results
    try:
        from petakit5d.jobs.framework import RUN_TASKS_PARAMS
    except ImportError:
        return {}
    return {p.matlab: settings.get(p.python, p.default) for p in RUN_TASKS_PARAMS}


def _effective_cpus(settings: Mapping[str, Any]) -> Any:
    try:
        from petakit5d.jobs.shapes import job_cpus

        return job_cpus(
            settings["cpus_per_task"],
            settings.get("mem_allocate"),
            settings["mem_per_cpu"],
            settings["min_cpu_num"],
            settings["max_cpu_num"],
            bool(settings.get("parse_cluster")) and not settings.get("whole_node_job"),
        )
    except Exception:
        return settings.get("cpus_per_task")


def session_dir(job_log_dir: str, task_metrics: bool = True) -> str:
    """The metrics folder of the single jobs this process submits to ``job_log_dir`` (one per
    process and folder, ``..._single_<id>``); '' without metrics or if the folder doesn't
    exist."""
    if not task_metrics:
        return ""
    try:
        folder = _sessions.get(job_log_dir)
        if folder and os.path.isdir(folder):
            return folder
        base = job_log_dir.replace("\\", "/")
        if not _is_absolute(base):
            base = os.getcwd().replace("\\", "/") + "/" + base
        if not os.path.isdir(base):
            return ""
        folder = _new_dir(base, "single")
        _sessions[job_log_dir] = folder
        return folder
    except Exception as err:
        _warn_run(err, "session_dir")
        return ""


def record(folder: str, kind: str, info: Mapping[str, Any]) -> None:
    """Appends a record of the run (master only) of the given kind ('run_start', 'run_end',
    'submit') to ``<folder>/run.ndjson``, with the master's node, pid, Slurm job and the task
    it runs in (``parent``)."""
    if not folder:
        return
    try:
        p = process_info()
        d = folder.replace("\\", "/").rstrip("/")
        r: dict[str, Any] = {
            "kind": kind,
            "version": 1,
            "runId": d[d.rfind("/") + 1 :],
            "time": time.time(),
            "node": p["host"],
            "pid": p["pid"],
            "jobId": _env("SLURM_JOB_ID"),
            "parent": p["current"],
        }
        r.update(info)
        text = json.dumps(_jsonable(r), ensure_ascii=False, separators=(",", ":"), allow_nan=False)
        with open(os.path.join(d, "run.ndjson"), "a", encoding="utf-8", newline="\n") as fh:
            fh.write(text + "\n")
    except Exception as err:
        _warn_run(err, "record")


def _jsonable(x: Any) -> Any:
    # JSON data as MATLAB's jsonencode writes it: NaN and Inf as null
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, Mapping):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if x is None or isinstance(x, (bool, int, str)):
        return x
    if isinstance(x, os.PathLike):
        return os.fspath(x)
    try:
        return _jsonable(x.tolist())  # NumPy values
    except AttributeError:
        return str(x)


# --- the reader ------------------------------------------------------------------------------


def read_task_metrics(
    metrics_path: str, *, slurm_accounting: bool = False
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """The task metrics of a run folder, a task_metrics folder or a job log folder (all runs in
    it) as two tables, lists of rows (dicts): a row per task run (a task that ran twice has
    two) and a row per run (framework call, or a process's single jobs). Reads MATLAB's
    folders and ours alike; MATLAB's read_task_metrics with lists for tables.

    Besides the recorded fields (times in seconds since 1970, memory in GB; compiled MATLAB
    tasks the master ran itself have mode 'master', or 'local' without a cluster), task rows
    have startupSec (process start to its first task, for tasks in a process of their own),
    lineSec and exitSec (MATLAB's GNU parallel joblogs; task lines without a record get a row
    with status 'skipped', 'failed' or 'killed'), input, output, cpusPerTask, memAllocate (from
    the run), and with ``slurm_accounting`` jobSubmit, jobStart, jobEnd, jobState, queueSec,
    jobHeadSec and jobTailSec from sacct (whole seconds). Missing numbers are NaN, missing
    texts ''."""
    if not isinstance(metrics_path, str):
        raise TypeError("metrics_path must be a str")
    if not isinstance(slurm_accounting, bool):
        raise TypeError("slurm_accounting must be a bool")
    return _read(metrics_path, slurm_accounting, OWN_PROCESS_MODES)


def _read(path: str, slurm_accounting: bool, own_modes: Sequence[str]):
    tasks: list[Any] = []
    runs: list[dict[str, Any]] = []
    lines: list[dict[str, Any]] = []
    for d in _find_run_dirs(path):
        run_id = d[:-1]
        run_id = run_id[run_id.rfind("/") + 1 :]
        run_recs = _read_records(d + "run.ndjson")
        task_recs: list[Any] = []
        names = _listdir(d)
        for name in names:
            if _match(name, "task_", ".ndjson"):
                task_recs += _read_records(d + name)
        for name in names:
            if _match(name, "parallel_", ".joblog"):
                lines += _read_joblog(d + name, run_id)
        tasks += task_recs
        runs.append({"runId": run_id, "folder": d, "records": run_recs, "nRecords": len(task_recs)})
    t = _task_table(tasks)
    _add_startup(t, own_modes)
    t = _add_joblog(t, lines)
    run_table = _add_run_settings(t, runs)
    if slurm_accounting:
        _add_slurm_accounting(t)
    return t, run_table


def _listdir(d: str) -> list[str]:
    try:
        return sorted(os.listdir(d))
    except OSError:
        return []


def _match(name: str, prefix: str, suffix: str) -> bool:
    return (
        name.startswith(prefix) and name.endswith(suffix) and len(name) >= len(prefix) + len(suffix)
    )


def _find_run_dirs(path: str) -> list[str]:
    p = re.sub("//", "/", path.replace("\\", "/") + "/")
    names = _listdir(p)
    if "run.ndjson" in names or any(_match(n, "task_", ".ndjson") for n in names):
        return [p]
    if os.path.isdir(p + "task_metrics"):
        p += "task_metrics/"
    return [p + n + "/" for n in _listdir(p) if not n.startswith(".") and os.path.isdir(p + n)]


def _no_constant(name: str) -> Any:
    raise ValueError(f"not JSON: {name}")


def _read_records(fn: str) -> list[Any]:
    # the JSON records of a file, one per line; lines that don't decode are left out
    try:
        with open(fn, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return []
    recs = []
    for line in re.split(r"\r\n|\n|\r", text):
        if not line.strip():
            continue
        try:
            recs.append(json.loads(line, parse_constant=_no_constant))
        except ValueError:
            pass
    return recs


def _str2double(s: str) -> float:
    try:
        x = float(s.strip())
    except ValueError:
        return NAN
    return NAN if math.isnan(x) else x


def _read_joblog(fn: str, run_id: str) -> list[dict[str, Any]]:
    # the lines of a GNU parallel joblog with the task index from their metrics environment
    rows = []
    name = os.path.splitext(os.path.basename(fn))[0]
    m = re.match(r"parallel_(\d+)_", name)
    job_id = m.group(1) if m else ""
    try:
        with open(fn, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return rows
    for line in re.split(r"\r\n|\n|\r", text):
        c = line.split("\t")
        if len(c) < 9 or c[0] == "Seq":
            continue
        k = re.search(r"PETAKIT5D_TASK=(\d+)", c[8])
        if not k:
            continue
        rows.append(
            {
                "runId": run_id,
                "task": _str2double(k.group(1)),
                "jobId": job_id,
                "lineStart": _str2double(c[2]),
                "lineSec": _str2double(c[3]),
                "exitCode": _str2double(c[6]),
                "signal": _str2double(c[7]),
            }
        )
    return rows


def _is_num(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _first_num(x: Any) -> float | None:
    # the first value of a JSON number or numeric array, as MATLAB's jsondecode gives it
    # (null in an array is NaN); None if it is not numeric or empty
    if _is_num(x):
        return float(x)
    if isinstance(x, list) and x:
        if all(_is_num(v) or v is None for v in x):
            return NAN if x[0] is None else float(x[0])
        if all(isinstance(v, list) for v in x):
            return _first_num(x[0])
    return None


def _task_table(recs: list[Any]) -> list[dict[str, Any]]:
    rows = []
    for r in recs:
        if not (isinstance(r, dict) and r.get("kind") == "task"):
            continue
        row: dict[str, Any] = {}
        for name in TASK_COLUMNS:
            if name in TASK_NUM_FIELDS:
                x = _first_num(r.get(name))
                row[name] = NAN if x is None else x
            else:
                v = r.get(name)
                row[name] = v if isinstance(v, str) else ""
        rows.append(row)
    return rows


def _empty_task_row() -> dict[str, Any]:
    return {n: (NAN if n in TASK_NUM_FIELDS else "") for n in TASK_COLUMNS}


def _round(x: float) -> float:
    # MATLAB's round: half away from zero
    if not math.isfinite(x):
        return x
    return math.copysign(math.floor(abs(x) + 0.5), x)


def _add_startup(t: list[dict[str, Any]], own_modes: Sequence[str]) -> None:
    # process start to its first task, for tasks in a process of their own
    groups: dict[tuple, list[int]] = {}
    for i, r in enumerate(t):
        r["startupSec"] = NAN
        if r["mode"] in own_modes and not math.isnan(r["processStart"]):
            key = (r["node"], _key_num(r["pid"]), _key_num(_round(r["processStart"])))
            groups.setdefault(key, []).append(i)
    for inds in groups.values():
        k = _nanargmin([t[i]["start"] for i in inds])
        r = t[inds[k]]
        r["startupSec"] = r["start"] - r["processStart"]


def _key_num(x: float) -> Any:
    return "NaN" if math.isnan(x) else x


def _nanargmin(xs: list[float]) -> int:
    # MATLAB's [~, k] = min(x): the first smallest, NaN ignored (all NaN: the first)
    k = 0
    best = math.inf
    found = False
    for i, x in enumerate(xs):
        if not math.isnan(x) and (not found or x < best):
            best, k, found = x, i, True
    return k


def _add_joblog(t: list[dict[str, Any]], lines: list[dict[str, Any]]) -> list[dict[str, Any]]:
    # GNU parallel's line run time and exit for compiled tasks; lines without a record get a row
    for r in t:
        r["lineSec"] = NAN
        r["exitSec"] = NAN
    extra = []
    for ln in lines:
        line_end = ln["lineStart"] + ln["lineSec"]
        for r in t:
            if (
                r["runId"] == ln["runId"]
                and r["task"] == ln["task"]
                and (r["jobId"] == ln["jobId"] or ln["jobId"] == "")
                and r["start"] >= ln["lineStart"] - 1
                and r["end"] <= line_end + 1
                and math.isnan(r["lineSec"])
            ):
                r["lineSec"] = ln["lineSec"]
                r["exitSec"] = line_end - r["end"]
                break
        else:
            row = _empty_task_row()
            row.update(startupSec=NAN, lineSec=ln["lineSec"], exitSec=NAN)
            row.update(task=ln["task"], runId=ln["runId"], jobId=ln["jobId"])
            row.update(start=ln["lineStart"], end=line_end, exitCode=ln["exitCode"])
            if ln["signal"] > 0:
                row["status"] = "killed"
            elif ln["exitCode"] != 0:
                row["status"] = "failed"
            else:
                row["status"] = "skipped"
            extra.append(row)
    return t + extra


def _get_str(r: dict, name: str) -> str:
    v = r.get(name)
    return v if isinstance(v, str) else ""


def _get_num(r: dict, name: str) -> float:
    v = r.get(name)
    if isinstance(v, bool):
        return 1.0 if v else 0.0
    if isinstance(v, list) and v and all(isinstance(x, bool) for x in v):
        return 1.0 if v[0] else 0.0
    x = _first_num(v)
    return NAN if x is None else x


def _get_cell(r: dict, name: str) -> list:
    # a list of a run record as MATLAB's get_cell has it: jsondecode gives a cell for a list
    # of texts or lists, but an array (one element of the cell) for a list of numbers, of
    # bools or of objects with the same fields, and for an empty list
    if name not in r:
        return []
    c = r[name]
    if isinstance(c, list) and c and not _decodes_to_array(c):
        return c
    return [c]


def _decodes_to_array(c: list) -> bool:
    if all(_is_num(x) or x is None for x in c) or all(isinstance(x, bool) for x in c):
        return True
    return all(isinstance(x, dict) for x in c) and len({tuple(x) for x in c}) == 1


def _first_str(x: Any) -> str:
    while isinstance(x, list) and x:
        x = x[0]
    return x if isinstance(x, str) else ""


def _strtok(s: str, delim: str = "(") -> str:
    s = s.lstrip(delim)
    k = s.find(delim)
    return s if k < 0 else s[:k]


def _sprintf_d(x: float) -> str:
    # MATLAB's sprintf('%d', x) of a double
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == int(x):
        return str(int(x))
    return f"{x:e}"


def _nanmin(a: float, b: float) -> float:
    if math.isnan(a):
        return b
    if math.isnan(b):
        return a
    return min(a, b)


def _add_run_settings(t: list[dict[str, Any]], runs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    for r in t:
        r.update(input="", output="", cpusPerTask=NAN, memAllocate=NAN)
    table = []
    for run in runs:
        rt: dict[str, Any] = dict.fromkeys(RUN_COLUMNS, NAN)
        rt.update(runId=run["runId"], kind="framework", func="", backend="", node="", parent="")
        rt.update(folder=run["folder"], nRecords=float(run["nRecords"]))
        in_run = [r for r in t if r["runId"] == run["runId"]]
        for rec in run["records"]:
            if not isinstance(rec, dict) or "kind" not in rec:
                continue
            kind = rec["kind"]
            if kind == "run_start":
                rt["node"] = _get_str(rec, "node")
                el = _get_num(rec, "elapsedSec")
                rt["start"] = _get_num(rec, "time") - (0.0 if math.isnan(el) else max(0.0, el))
                rt["backend"] = _get_str(rec, "backend")
                rt["parent"] = _get_str(rec, "parent")
                for name in ("parseCluster", "mccMode", "cpusPerTask", "memAllocate",
                             "taskBatchNum", "paraJobNum", "GPUJob"):  # fmt: skip
                    rt[name] = _get_num(rec, name)
                inputs = _get_cell(rec, "inputFullpaths")
                outputs = _get_cell(rec, "outputFullpaths")
                funcs = _get_cell(rec, "functionStrs")
                rt["nTasks"] = float(len(outputs))
                if funcs:
                    rt["func"] = _strtok(_first_str(funcs[0]))
                for r in in_run:
                    f = r["task"]
                    if f >= 1 and f <= len(outputs):
                        r["output"] = _first_str(outputs[int(f) - 1])
                    if f >= 1 and f <= len(inputs):
                        r["input"] = _first_str(inputs[int(f) - 1])
                    r["cpusPerTask"] = rt["cpusPerTask"]
                    r["memAllocate"] = rt["memAllocate"]
                # compiled tasks in the master's Slurm job (or, like the master, in none) are
                # the ones the master ran itself
                job = _get_str(rec, "jobId")
                for r in in_run:
                    if r["mode"] in ("slurm_mcc", "mcc") and r["jobId"] == job:
                        r["mode"] = "local" if rt["parseCluster"] == 0 else "master"
            elif kind == "run_end":
                rt["end"] = _get_num(rec, "time")
                rt["elapsedSec"] = _get_num(rec, "elapsedSec")
                rt["nTasks"] = _get_num(rec, "nTasks")
                rt["nDone"] = _get_num(rec, "nDone")
            elif kind == "submit":
                # a single job: the tasks that ran in it, by its array job id
                rt["kind"] = "single"
                rt["node"] = _get_str(rec, "node")
                rt["start"] = _nanmin(rt["start"], _get_num(rec, "time"))
                rt["nTasks"] = 1.0 if math.isnan(rt["nTasks"]) else rt["nTasks"] + 1
                job = _sprintf_d(_get_num(rec, "jobId"))
                rows = [r for r in in_run if r["arrayJobId"] == job]
                for r in rows:
                    r["cpusPerTask"] = _get_num(rec, "cpusPerTask")
                    r["memAllocate"] = _get_num(rec, "memAllocate")
                if not rows:
                    continue
                func = _get_str(rec, "functionStr")
                if func:
                    rt["func"] = _strtok(func)
        table.append(rt)
    return table


def _sacct_time(s: str) -> float:
    try:
        return time.mktime(time.strptime(s.strip(), "%Y-%m-%dT%H:%M:%S"))
    except (ValueError, OverflowError):
        return NAN


def _add_slurm_accounting(t: list[dict[str, Any]]) -> None:
    # each task's Slurm job from sacct: submit, start, end, state, queue wait, and the job's
    # time before its first and after its last task
    import subprocess

    for r in t:
        r.update(dict.fromkeys(SACCT_COLUMNS, NAN))
        r["jobState"] = ""
    ids = sorted({r["jobId"] for r in t if r["jobId"]})
    by_job: dict[str, list[dict[str, Any]]] = {}
    for r in t:
        by_job.setdefault(r["jobId"], []).append(r)
    for b in range(0, len(ids), 500):
        bids = ids[b : b + 500]
        cmd = [
            "sacct",
            "-n",
            "-P",
            "-X",
            "--format=JobIDRaw,Submit,Start,End,State",
            "-j",
            ",".join(bids),
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            status, out = res.returncode, res.stdout + res.stderr
        except (OSError, subprocess.TimeoutExpired) as err:
            status, out = 1, str(err)
        if status != 0:
            warnings.warn(
                TaskMetricsWarning(f"sacct failed: {out.strip()}", "read_task_metrics:sacct"),
                stacklevel=3,
            )
            break
        for line in res.stdout.strip().splitlines():
            c = line.split("|")
            if len(c) < 5:
                continue
            times = [_sacct_time(x) for x in c[1:4]]
            state = c[4].split()[0] if c[4].split() else ""
            for r in by_job.get(c[0], []):
                r["jobSubmit"], r["jobStart"], r["jobEnd"] = times
                r["jobState"] = state
                r["queueSec"] = times[1] - times[0]
    for job in ids:
        rows = by_job[job]
        starts = [r["start"] for r in rows if not math.isnan(r["start"])]
        ends = [r["end"] for r in rows if not math.isnan(r["end"])]
        head = (min(starts) if starts else NAN) - rows[0]["jobStart"]
        tail = rows[0]["jobEnd"] - (max(ends) if ends else NAN)
        for r in rows:
            r["jobHeadSec"], r["jobTailSec"] = head, tail
