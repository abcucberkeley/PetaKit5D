"""Jobs that coordinate jobs: what they ask for, what their master computes, and their nice.

A coordinating job is a job whose task runs tasks of its own as jobs (a large-file task and
its block jobs, a stitching job and its blending jobs) and waits for them. Three rules, as in
MATLAB's job framework (coordinator_resources.m, master_task_fits.m, coordinator_nice.m):

- ``coordinator_resources``: it asks for what it holds itself (the process and its merges,
  with one CPU), or, if it computes some of its tasks itself, for one of these tasks too.
- ``master_task_fits``: a master in a Slurm job runs a task itself only if the task fits the
  job's memory (up to 5% over); outside Slurm always.
- ``coordinator_nice``: it is submitted with ``--nice`` the part of a job's priority that age
  and size can add (Slurm's multifactor priority), so that the jobs it waits for go first.

One value differs from MATLAB: the memory of the process itself (``PROCESS_MEMORY``), 1 GB for
Python where MATLAB counts 2 GB for a MATLAB process (measured; see the constant).
"""

from __future__ import annotations

import math
import os
import re
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, NamedTuple

if TYPE_CHECKING:
    from petakit5d.jobs.slurm import Slurm

__all__ = [
    "FIT_SLACK",
    "MERGE_MEMORY_FACTOR",
    "PROCESS_MEMORY",
    "MasterFit",
    "coordinator_nice",
    "coordinator_resources",
    "master_fit",
    "master_task_fits",
    "settings_master_fit",
]

# GB a coordinating Python process holds besides its merges (MATLAB: 2 GB for MATLAB, whose
# coordinators peaked at 1.5-1.7 GB). Measured on a Linux compute node (2026-11): the
# interpreter with the job framework, NumPy, SciPy, scikit-image, cpp-zarr, cpp-tiff, h5py and
# every petakit5d module peaks at 0.12 GB; with a call of 10,000 task records (30 keyword
# arguments each) written and its outputs looked at, 0.18 GB; of 100,000, 0.74 GB (0.85 GB in
# the job's cgroup). 1 GB covers that.
PROCESS_MEMORY = 1.0
# a merge (XR_MIP_zarr's slab merge) takes at most half of the job's memory
MERGE_MEMORY_FACTOR = 2
# memPerCPU is often a rounded value of the jobs' memory per CPU: up to 5% over fits
FIT_SLACK = 1.05

# PriorityWeightAge and PriorityWeightJobSize lines of scontrol show config (MATLAB's regexp
# with 'lineanchors'; \s also matches line ends there)
_WEIGHTS = re.compile(r"^\s*PriorityWeight(Age|JobSize)\s*=\s*(\d+)", re.M)


def _mmax(a: float, b: float) -> float:
    """MATLAB's max of two doubles (a NaN is ignored, the first of equals is kept)"""
    if math.isnan(a):
        return b
    if math.isnan(b):
        return a
    return b if b > a else a


def _mmin(a: float, b: float) -> float:
    """MATLAB's min of two doubles (a NaN is ignored)"""
    if math.isnan(a):
        return b
    if math.isnan(b):
        return a
    return b if b < a else a


def _empty(x: Any) -> bool:
    # MATLAB's isempty for the optional numbers ([] is None here; an empty list is accepted)
    return x is None or (hasattr(x, "__len__") and len(x) == 0)


def coordinator_resources(
    merge_mem: float, computes: bool, task_cpus: float, task_mem: float | None
) -> tuple[float, float]:
    """coordinator_resources.m: (cpus_per_task, mem_allocate) of a coordinating job.

    merge_mem : GB of its merges (XR_MIP_zarr merges its MIP slabs in up to half of the job).
    computes  : whether it computes some of its tasks itself (master_compute, for tasks that
                run on the same kind of node).
    task_cpus : CPUs the caller asks for its jobs (e.g. cpus_per_task, the blocks' CPUs).
    task_mem  : GB one of its tasks gets as a job (None: the memory of its CPUs).

    The job holds the process (PROCESS_MEMORY GB) and its merges, with one CPU:
    max(PROCESS_MEMORY + merge_mem, 2 x merge_mem) GB. If it computes, it gets task_cpus and
    at least task_mem, so that it is one of the workers instead of an idle job next to them.
    Its master runs a task only if the task fits the job (master_task_fits).
    """
    merge_mem = float(merge_mem)
    mem_allocate = _mmax(PROCESS_MEMORY + merge_mem, MERGE_MEMORY_FACTOR * merge_mem)
    cpus_per_task: float = 1
    if computes:
        cpus_per_task = task_cpus
        if not _empty(task_mem):
            mem_allocate = _mmax(mem_allocate, float(task_mem))  # type: ignore[arg-type]
    return cpus_per_task, mem_allocate


class MasterFit(NamedTuple):
    """master_fit's answer: whether the master may run the task, the task's memory and the
    master's (GB; inf outside a Slurm job, NaN if unknown)"""

    fits: bool
    task_mem: float
    master_mem: float


def master_fit(
    mem_allocate: float | None,
    cpus_per_task: float,
    mem_per_cpu: float,
    max_cpu_num: float | None = None,
    min_cpu_num: float | None = None,
    master_mem: float | None = None,
) -> MasterFit:
    """master_task_fits.m with all three outputs (see master_task_fits).

    The task's memory is mem_allocate (at most max_cpu_num x mem_per_cpu, what its job can
    get), or without it the memory of its job's CPUs (cpus_per_task within min_cpu_num and
    max_cpu_num, times mem_per_cpu). Outside a Slurm job (no SLURM_JOB_ID) every task fits
    (the user's process has all of the machine). In a job, the task fits if it needs at most
    FIT_SLACK x the job's memory: master_mem if given (callers that check many tasks read it
    once), else get_system_memory(); an unknown (NaN) memory doesn't fit.
    """
    max_cpu = math.inf if _empty(max_cpu_num) else float(max_cpu_num)  # type: ignore[arg-type]
    min_cpu = 1.0 if _empty(min_cpu_num) else float(min_cpu_num)  # type: ignore[arg-type]
    if _empty(mem_allocate):
        task_mem = _mmin(max_cpu, _mmax(min_cpu, float(cpus_per_task))) * float(mem_per_cpu)
    else:
        task_mem = _mmin(float(mem_allocate), max_cpu * float(mem_per_cpu))  # type: ignore[arg-type]

    if not os.environ.get("SLURM_JOB_ID"):
        return MasterFit(True, task_mem, math.inf)
    if _empty(master_mem):
        from petakit5d.core import memory

        master_mem = memory.get_system_memory()
    master_mem = float(master_mem)  # type: ignore[arg-type]
    return MasterFit(task_mem <= master_mem * FIT_SLACK, task_mem, master_mem)


def master_task_fits(
    mem_allocate: float | None,
    cpus_per_task: float,
    mem_per_cpu: float,
    max_cpu_num: float | None = None,
    min_cpu_num: float | None = None,
    master_mem: float | None = None,
) -> bool:
    """master_task_fits.m: whether the master of a cluster-mode call may run a task itself.

    Outside a Slurm job always; in a Slurm job (a job that coordinates jobs, or a pipeline's
    job) only if the task's memory fits the job's, so that the job isn't killed for memory with
    a task its own job would have run. Memories in GB; max_cpu_num defaults to inf and
    min_cpu_num to 1. ``master_fit`` gives the task's and the master's memory too.
    """
    return master_fit(
        mem_allocate, cpus_per_task, mem_per_cpu, max_cpu_num, min_cpu_num, master_mem
    ).fits


def settings_master_fit(settings: Mapping[str, Any]) -> MasterFit:
    """The fit run_tasks checks in cluster mode with master_compute (MATLAB's
    generic_computing_frameworks_wrapper): the call's mem_allocate, cpus_per_task, mem_per_cpu,
    max_cpu_num and min_cpu_num (config file applied)."""
    return master_fit(
        settings.get("mem_allocate"),
        settings["cpus_per_task"],
        settings["mem_per_cpu"],
        settings["max_cpu_num"],
        settings["min_cpu_num"],
    )


def coordinator_nice(slurm: Slurm | None = None) -> int:
    """coordinator_nice.m: the nice (sbatch --nice) of a job that submits jobs and waits for
    them, so that the jobs it waits for go first.

    With Slurm's multifactor priority a job's priority grows with its age (up to
    PriorityWeightAge) and its size (up to PriorityWeightJobSize), while fair share, partition
    and QOS are the same for a user's jobs of one partition and QOS. Lowered by the sum of the
    two weights (from ``scontrol show config``), a coordinating job stays below the jobs it
    submits, whatever their ages. 0 without Slurm (scontrol fails or doesn't answer in time),
    or with another priority type.
    """
    if slurm is None:
        from petakit5d.jobs.slurm import Slurm

        slurm = Slurm()
    r = slurm.run(["scontrol", "show", "config"])
    out = r.output  # (stdout and stderr, as MATLAB's system gives them)
    if not r.ok or "priority/multifactor" not in out:
        return 0
    return sum(int(m.group(2)) for m in _WEIGHTS.finditer(out))
