"""``getSystemMemory``: the memory (GB, ``2^30`` bytes) the job may use; NaN when unknown.

The same decisions as MATLAB's ``getSystemMemory`` (so both choose the same code paths):

- macOS: ``sysctl -n hw.memsize``.
- Linux in a Slurm job (``SLURM_JOB_ID`` set): the smallest memory limit on this process's
  cgroup path (v1 ``memory.limit_in_bytes`` under ``/sys/fs/cgroup/memory``, v2 ``memory.max``),
  if that path is the job's own (``/job_<id>``); else ``SLURM_MEM_PER_NODE`` (MB), unless
  ``SLURM_MEM_PER_CPU`` is also set (then one of them is another job's); else the ``mem=``
  value of ``scontrol show job <id>[_<array id>]``.
- Linux without Slurm: ``MemTotal`` of ``/proc/meminfo``.
- Windows: MATLAB's ``memory().MaxPossibleArrayBytes``, taken as the smaller of the available
  commit (physical plus page file) and the physical memory (MATLAB's default array size
  limit). Not compared with MATLAB yet.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
import sys
from collections.abc import Callable, Mapping

from petakit5d.core.text import str2double

__all__ = ["get_system_memory", "job_cgroup_memory", "linux_memory"]


def _read_file(path: str) -> str:
    with open(path) as f:
        return f.read()


def _run(cmd: list[str]) -> str:
    """Command output (stdout and stderr, as MATLAB's ``system``); '' if it can't start."""
    try:
        r = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return r.stdout or ""


def _parent(p: str) -> str:
    # MATLAB's fileparts(p) on a cgroup path
    i = p.rfind("/")
    if i < 0:
        return ""
    return p[:i] or "/"


def job_cgroup_memory(job_id: str, read_file: Callable[[str], str] = _read_file) -> float:
    """The smallest memory limit (GB) on this process's cgroup path if it is job ``job_id``'s."""
    try:
        txt = read_file("/proc/self/cgroup")
    except OSError:
        return math.nan
    m = re.search(r"^\d+:[^:\n]*\bmemory\b[^:\n]*:(\S+)$", txt, re.M)
    if m:
        root, limit_name = "/sys/fs/cgroup/memory", "memory.limit_in_bytes"
    else:
        m = re.search(r"^0::(\S+)$", txt, re.M)
        root, limit_name = "/sys/fs/cgroup", "memory.max"
    if m is None:
        return math.nan
    try:
        if not re.search("/job_" + job_id + "(/|$)", m.group(1)):
            return math.nan
    except re.error:
        return math.nan
    cg = m.group(1)
    limit = math.inf
    while cg and cg != "/":
        try:
            v = str2double(read_file(root + cg + "/" + limit_name).strip())
            if v > 0:
                limit = min(limit, v)
        except OSError:
            pass
        cg = _parent(cg)
    # (cgroup v1 shows no limit as a huge number)
    return limit / 2**30 if limit < 2**60 else math.nan


def linux_memory(
    env: Mapping[str, str] | None = None,
    read_file: Callable[[str], str] = _read_file,
    run: Callable[[list[str]], str] = _run,
) -> float:
    """``getSystemMemory``'s Linux branch, with its inputs injectable for tests."""
    env = os.environ if env is None else env
    job_id = env.get("SLURM_JOB_ID", "")
    if job_id:
        total = job_cgroup_memory(job_id, read_file)
        if math.isnan(total) and not env.get("SLURM_MEM_PER_CPU", ""):
            total = str2double(env.get("SLURM_MEM_PER_NODE", "")) / 1024
        if math.isnan(total):
            array_id = env.get("SLURM_ARRAY_TASK_ID", "")
            job_array_id = job_id if not array_id else f"{job_id}_{array_id}"
            out = run(["scontrol", "show", "job", job_array_id])
            m = re.search(r"mem=([\d.]+[TGMK])", out)
            if m:
                v = m.group(1)
                denom = {"T": 1 / 1024, "G": 1, "M": 1024, "K": 1024**2}[v[-1]]
                total = str2double(v[:-1]) / denom
            else:
                print("No match found")
        return total
    try:
        txt = read_file("/proc/meminfo")
    except OSError:
        txt = ""
    m = re.search(r"MemTotal:\s+(\d+)\s+kB", txt)
    if m:
        return str2double(m.group(1)) / 1024**2
    print("No match found")
    return math.nan


def _windows_memory() -> float:
    import ctypes

    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    st = MEMORYSTATUSEX()
    st.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):  # type: ignore[attr-defined]
        return math.nan
    return min(st.ullAvailPageFile, st.ullTotalPhys, st.ullAvailVirtual) / 1024**3


def get_system_memory() -> float:
    """``getSystemMemory()``: memory in GB (2^30 bytes) for choosing code paths; NaN if unknown."""
    if sys.platform == "darwin":
        out = _run(["sysctl", "-n", "hw.memsize"])
        return str2double(out.strip()) / 1024**3
    if sys.platform.startswith("win"):
        return _windows_memory()
    return linux_memory()
