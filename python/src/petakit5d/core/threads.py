"""Thread counts from the machine: the CPUs this process may use, and per-call counts.

``available_cpus()`` is the smaller of the CPU affinity (Slurm's cpuset, ``taskset``) and the
cgroup CPU quota (``cpu.max`` in v2, ``cpu.cfs_quota_us / cpu.cfs_period_us`` in v1, rounded
up), never an environment variable. Kernels take their thread count as an argument:
``resolve_threads(n)`` turns a user's ``n`` (None or 0 for all) into that count. Nothing here
sets global state (no ``OMP_NUM_THREADS``).

The affinity is read once, when this module is first imported (importing any part of
``petakit5d.core`` imports it), before any kernel has run: with ``OMP_PROC_BIND`` set, libgomp
binds the main thread to one CPU at its first parallel region, and the process's affinity read
after that is that one CPU.

``available_cpus()`` is computed once per process (the cgroup files cost about 0.1 ms to read,
which showed in the per-call overhead of small kernels) and then returned from a cache.
``refresh_cpus()`` reads the affinity and the cgroup quota again and replaces the cached value:
for tests and for a process whose affinity or cgroup changed. Call it before any kernel has
run, or with ``OMP_PROC_BIND`` unset: after a bound OpenMP region it would read one CPU.
"""

from __future__ import annotations

import math
import os
import re
import sys

__all__ = ["available_cpus", "cgroup_cpu_limit", "refresh_cpus", "resolve_threads"]


def _affinity_now() -> int:
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    if hasattr(os, "process_cpu_count"):  # Python 3.13: affinity on Windows too
        n = os.process_cpu_count()
        if n:
            return n
    return os.cpu_count() or 1


# taken before any OpenMP region can bind the main thread (see the module's docstring)
_AFFINITY_AT_IMPORT = _affinity_now()


def _affinity_count() -> int:
    """CPUs in this process's affinity when petakit5d.core was first imported."""
    return _AFFINITY_AT_IMPORT


def _read(path: str) -> str | None:
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return None


def cgroup_cpu_limit(proc_cgroup: str | None = None, root: str = "/sys/fs/cgroup") -> float | None:
    """The CPU quota of this process's cgroups in CPUs (smallest along the path), or None.

    ``proc_cgroup`` is the text of ``/proc/self/cgroup`` (read when None); ``root`` the cgroup
    mount (tests point both at a fake tree).
    """
    if not sys.platform.startswith("linux") and proc_cgroup is None:
        return None
    txt = _read("/proc/self/cgroup") if proc_cgroup is None else proc_cgroup
    if not txt:
        return None
    limits = []
    m1 = re.search(r"^\d+:([^:\n]*\bcpu\b[^:\n]*):(\S+)$", txt, re.M)
    if m1:
        controllers, path = m1.group(1), m1.group(2)
        bases = [
            os.path.join(root, controllers),
            os.path.join(root, "cpu"),
            os.path.join(root, "cpu,cpuacct"),
        ]
        base = next((b for b in bases if os.path.isdir(b)), None)
        if base is None:
            return None
        p = path
        while True:
            q = _read(os.path.join(base, p.lstrip("/"), "cpu.cfs_quota_us"))
            t = _read(os.path.join(base, p.lstrip("/"), "cpu.cfs_period_us"))
            if q is not None and t is not None:
                try:
                    quota, period = int(q), int(t)
                    if quota > 0 and period > 0:
                        limits.append(quota / period)
                except ValueError:
                    pass
            if p in ("", "/"):
                break
            p = os.path.dirname(p)
    else:
        m2 = re.search(r"^0::(\S+)$", txt, re.M)
        if not m2:
            return None
        p = m2.group(1)
        while True:
            v = _read(os.path.join(root, p.lstrip("/"), "cpu.max"))
            if v is not None:
                fields = v.split()
                if len(fields) >= 1 and fields[0] != "max":
                    try:
                        quota = int(fields[0])
                        period = int(fields[1]) if len(fields) > 1 else 100000
                        if quota > 0 and period > 0:
                            limits.append(quota / period)
                    except ValueError:
                        pass
            if p in ("", "/"):
                break
            p = os.path.dirname(p)
    return min(limits) if limits else None


def _available_uncached() -> int:
    n = _affinity_count()
    q = cgroup_cpu_limit()
    if q is not None:
        n = min(n, max(1, math.ceil(q - 1e-9)))
    return max(1, n)


_AVAILABLE: int | None = None  # the cached value of available_cpus()


def available_cpus() -> int:
    """CPUs this process may use: affinity, limited by the cgroup CPU quota (at least 1).

    Computed on the first call and cached for the process (``refresh_cpus`` recomputes it).
    """
    global _AVAILABLE
    n = _AVAILABLE
    if n is None:
        n = _AVAILABLE = _available_uncached()
    return n


def refresh_cpus() -> int:
    """Read the affinity and the cgroup CPU quota again, cache and return ``available_cpus()``.

    For tests and for a process whose affinity or cgroup changed; call it before any OpenMP
    region has bound the main thread (see the module's docstring).
    """
    global _AFFINITY_AT_IMPORT, _AVAILABLE
    _AFFINITY_AT_IMPORT = _affinity_now()
    _AVAILABLE = None
    return available_cpus()


def resolve_threads(n: int | None = None) -> int:
    """The thread count for one call: all available CPUs for None or 0, else ``n`` limited to
    them (a negative ``n`` raises)."""
    avail = available_cpus()
    if n is None or n == 0:
        return avail
    n = int(n)
    if n < 0:
        raise ValueError(f"thread count must be positive, got {n}")
    return min(n, avail)
