"""Jobs that coordinate jobs: master_task_fits, coordinator_resources, coordinator_nice.

Ported from tests/framework/MasterComputeFitTest.m and CoordinatorNiceTest.m; the test names
keep the MATLAB names in snake_case. The framework tests run run_tasks in cluster mode against
the fake Slurm in two ways: with a recording stand-in for the cluster loop (what run_tasks hands
it: master_compute, the number of tasks, SlurmParam), and through the cluster loop itself
(petakit5d.jobs.cluster, task J03; skipped until it exists), as the MATLAB tests do. The oracle
cases compare master_task_fits and coordinator_resources with MATLAB on grids of their inputs.
The cluster test submits a coordinating job and its block jobs to a real Slurm.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
import types
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.dtypes import NAN
from petakit5d.jobs import coordinator, framework
from petakit5d.jobs.coordinator import (
    PROCESS_MEMORY,
    coordinator_nice,
    coordinator_resources,
    master_fit,
    master_task_fits,
)
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.tasks import Task
from petakit5d.testing import oracle, slurm_fake
from petakit5d.testing.cases import Case, In, Value

fake_slurm = slurm_fake.fake_slurm  # (the fixture)

posix_only = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")

# MATLAB's allowance for a MATLAB process (coordinator_resources.m); Python's is PROCESS_MEMORY
MATLAB_PROCESS_MEMORY = 2


# --- helpers --------------------------------------------------------------------------------


def framework_task(input_fn: str, output_fn: str, log_fn: str) -> None:
    """a task (pktest.frameworkTask 'copy'): logs its output path, copies its input there"""
    fd = os.open(log_fn, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        os.write(fd, (output_fn + "\n").encode())
    finally:
        os.close(fd)
    Path(output_fn).write_bytes(Path(input_fn).read_bytes())


def read_log(log_fn: Path) -> list[str]:
    """the output paths of all task calls this process made, in order"""
    return log_fn.read_text().split() if log_fn.exists() else []


@pytest.fixture
def slurm_job(monkeypatch):
    """set_job(job_id, memory): the tests' Slurm job ('' for none) and its memory in GB
    (get_system_memory stubbed, as the MATLAB tests stub getSystemMemory)"""

    def set_job(job_id: str, memory: float = math.nan) -> None:
        for name in slurm_fake.SLURM_JOB_VARIABLES:
            monkeypatch.delenv(name, raising=False)
        if job_id:
            monkeypatch.setenv("SLURM_JOB_ID", job_id)
        from petakit5d.core import memory as core_memory

        monkeypatch.setattr(core_memory, "get_system_memory", lambda: memory)

    return set_job


@pytest.fixture
def fits_in_framework(monkeypatch):
    # run_tasks decides with coordinator.settings_master_fit; until J02's stand-in
    # (framework._master_task_fits) is replaced by it, the tests put it in its place
    current = getattr(framework, "_master_task_fits", None)
    if current is not None and getattr(current, "__module__", "") == framework.__name__:
        monkeypatch.setattr(framework, "_master_task_fits", coordinator.settings_master_fit)


class RecordingCluster:
    """stand-in for petakit5d.jobs.cluster: records the calls run_tasks hands to the cluster
    loop; its jobs end without results (as the fake Slurm's jobs)"""

    def __init__(self) -> None:
        self.calls: list[tuple[int, dict]] = []

    def run_cluster(self, call, slurm=None):
        self.calls.append((call.n, dict(call.settings)))
        return [False] * call.n


@pytest.fixture
def recording_cluster(monkeypatch, fake_slurm, fits_in_framework):
    import petakit5d.jobs as jobs_pkg

    rec = RecordingCluster()
    module = types.ModuleType("petakit5d.jobs.cluster")
    module.run_cluster = rec.run_cluster
    monkeypatch.setitem(sys.modules, "petakit5d.jobs.cluster", module)
    monkeypatch.setattr(jobs_pkg, "cluster", module, raising=False)
    return rec


class VirtualClock:
    """time of the cluster loop that passes only in its pauses (as J03's tests use it)"""

    def __init__(self) -> None:
        self.now = 1000.0

    def clock(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.now += seconds + 0.01


@pytest.fixture
def cluster_loop(monkeypatch, fake_slurm, fits_in_framework):
    """the real cluster loop (J03) with the fake Slurm, its pauses virtual when it has them"""
    cluster = pytest.importorskip("petakit5d.jobs.cluster", reason="needs the cluster loop (J03)")
    if hasattr(cluster, "_clock") and hasattr(cluster, "_sleep"):
        vc = VirtualClock()
        monkeypatch.setattr(cluster, "_clock", vc.clock)
        monkeypatch.setattr(cluster, "_sleep", vc.sleep)
    return fake_slurm


def run_framework(folder: Path, n: int, **params):
    """run_tasks in cluster mode with master_compute on n copy tasks, one try (MATLAB's
    runFramework); returns the log file and the outputs"""
    (folder / "out").mkdir(parents=True)
    (folder / "logs").mkdir()
    ins = [str(folder / f"in_{k}.txt") for k in range(1, n + 1)]
    outs = [str(folder / "out" / f"out_{k}.txt") for k in range(1, n + 1)]
    for k, fn in enumerate(ins, 1):
        Path(fn).write_text(f"task {k}")
    log_fn = folder / "calls.log"
    tasks = [Task.of(framework_task, i, o, str(log_fn)) for i, o in zip(ins, outs, strict=True)]
    args = {
        "parse_cluster": True,
        "master_compute": True,
        "max_trial_num": 1,
        "query_interval": 0,
        "job_log_dir": str(folder / "logs"),
        "task_runner": False,
    }
    args.update(params)
    run_tasks(ins, outs, tasks, **args)
    return log_fn, outs


# --- MasterComputeFitTest: master_task_fits and coordinator_resources -----------------------


def test_fit_is_the_tasks_memory_against_the_jobs(slurm_job):
    slurm_job("4242", 16)
    # memAllocate, up to 5% over (memPerCPU is often rounded)
    assert master_task_fits(16.5, 1, 20.9, 24, 1) is True
    assert master_task_fits(17, 1, 20.9, 24, 1) is False
    # at most what a job gets (maxCPUNum x memPerCPU)
    assert master_task_fits(100, 1, 10, 1, 1) is True
    # without memAllocate, the memory of the job's CPUs (at least minCPUNum, at most maxCPUNum)
    assert master_task_fits(None, 2, 10, 24, 1) is False
    assert master_task_fits(None, 2, 10, 1, 1) is True
    assert master_task_fits(None, 1, 10, 24, 2) is False
    fit = master_fit(None, 3, 5, 24, 1)
    assert (fit.task_mem, fit.master_mem) == (15, 16)
    # (the defaults: maxCPUNum inf, minCPUNum 1; [] is None)
    assert master_fit(None, 3, 5).task_mem == 15
    assert master_fit([], 0, 5, [], []).task_mem == 5
    assert master_fit(None, 50, 1).task_mem == 50


def test_unknown_memory_does_not_fit(slurm_job):
    # (get_system_memory gives NaN when it can't find the memory)
    slurm_job("4242", math.nan)
    assert master_task_fits(1, 1, 20.9, 24, 1) is False


def test_every_task_fits_outside_a_job(slurm_job):
    slurm_job("", 1)
    fit = master_fit(1000, 24, 20.9, 24, 1)
    assert fit.fits is True
    assert fit.master_mem == math.inf


def test_the_jobs_memory_is_read_once_given(slurm_job, monkeypatch):
    # callers that check many tasks pass the job's memory (MATLAB's sixth argument)
    slurm_job("4242", 1)
    assert master_fit(30, 1, 20.9, 24, 1, 32).fits is True
    assert master_fit(30, 1, 20.9, 24, 1, 32).master_mem == 32
    # outside a job the given memory doesn't matter
    slurm_job("", 1)
    assert master_fit(30, 1, 20.9, 24, 1, 1) == (True, 30, math.inf)


def test_master_memory_comes_from_get_system_memory(monkeypatch):
    # in a job, without the job's memory: the memory Slurm gives the job (getSystemMemory)
    from petakit5d.core import memory as core_memory

    for name in slurm_fake.SLURM_JOB_VARIABLES:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "4242")
    asked = []
    monkeypatch.setattr(core_memory, "get_system_memory", lambda: asked.append(1) or 8.0)
    assert master_fit(8.4, 1, 20.9).fits is True
    assert master_fit(8.41, 1, 20.9).fits is False
    assert len(asked) == 2


def test_coordinator_asks_for_what_it_holds():
    # the process (PROCESS_MEMORY GB) and the merges with one CPU; a merge needs up to half of
    # the job's memory (XR_MIP_zarr)
    assert coordinator_resources(0.5, False, 24, 40) == (1, PROCESS_MEMORY + 0.5)
    assert coordinator_resources(3, False, 24, 40) == (1, 6)
    # one of its tasks too, if it computes them
    assert coordinator_resources(0.5, True, 24, 40) == (24, 40)
    assert coordinator_resources(0.5, True, 24, None) == (24, PROCESS_MEMORY + 0.5)
    assert coordinator_resources(30, True, 2, 40) == (2, 60)


def test_coordinator_asks_for_what_it_holds_with_matlabs_allowance(monkeypatch):
    # MATLAB's values (2 GB for a MATLAB process): only the allowance differs
    monkeypatch.setattr(coordinator, "PROCESS_MEMORY", MATLAB_PROCESS_MEMORY)
    assert coordinator_resources(0.5, False, 24, 40) == (1, 2.5)
    assert coordinator_resources(3, False, 24, 40) == (1, 6)
    assert coordinator_resources(0.5, True, 24, 40) == (24, 40)
    assert coordinator_resources(0.5, True, 24, None) == (24, 2.5)
    assert coordinator_resources(30, True, 2, 40) == (2, 60)


def test_process_allowance_is_the_measured_one():
    # 1 GB for a Python coordinator (MATLAB: 2 GB); changed only with a new measurement
    # (0.12 GB with every import, 0.74 GB with a call of 100,000 tasks; see the constant)
    assert PROCESS_MEMORY == 1.0


# --- MasterComputeFitTest: the framework (the master's decision) ----------------------------


@posix_only
@pytest.mark.parametrize("n", [1, 3])
def test_master_in_a_job_runs_no_task_that_does_not_fit_it(
    n, tmp_path, slurm_job, recording_cluster, capsys
):
    # a master in a Slurm job ran the last task (or a single task) itself whatever its
    # estimate: a user's 4-CPU, 64 GB pipeline job ran an Imaris conversion and was killed for
    # memory. Now tasks that don't fit its memory all run as jobs
    slurm_job("4242", 16)
    log_fn, _ = run_framework(tmp_path, n, mem_allocate=40)
    assert read_log(log_fn) == [], f"{n} task(s): the master ran a task."
    assert len(recording_cluster.calls) == 1, f"{n} task(s) were not all handed to the jobs."
    n_tasks, settings = recording_cluster.calls[0]
    assert n_tasks == n
    assert settings["master_compute"] is False
    assert "do not fit this job" in capsys.readouterr().out


@posix_only
def test_master_in_a_job_runs_tasks_that_fit_it(tmp_path, slurm_job, recording_cluster):
    # with the memory of a task, the master runs the last one as before (master compute stays
    # on for the cluster loop)
    slurm_job("4242", 64)
    run_framework(tmp_path, 3, mem_allocate=40)
    assert [(n, s["master_compute"]) for n, s in recording_cluster.calls] == [(3, True)]


@posix_only
def test_master_in_a_job_runs_a_single_task_that_fits_it(tmp_path, slurm_job, recording_cluster):
    # one task that fits runs in the master, without a job
    slurm_job("4242", 64)
    log_fn, outs = run_framework(tmp_path, 1, mem_allocate=40)
    assert read_log(log_fn) == outs
    assert recording_cluster.calls == []


@posix_only
@pytest.mark.parametrize(("cpus", "computes"), [(4, False), (1, True)])
def test_tasks_without_an_estimate_need_the_memory_of_their_jobs(
    cpus, computes, tmp_path, slurm_job, recording_cluster
):
    # without memAllocate, a task's job has the memory of its CPUs (4 x 10 GB here), which is
    # what the master needs
    slurm_job("4242", 16)
    run_framework(tmp_path, 2, cpus_per_task=cpus, mem_per_cpu=10)
    assert [(n, s["master_compute"]) for n, s in recording_cluster.calls] == [(2, computes)]


@posix_only
def test_master_outside_a_job_runs_tasks_as_before(tmp_path, slurm_job, recording_cluster):
    # the user's process (no Slurm job) has all of the machine: it computes whatever the
    # estimate, as before
    slurm_job("", 1)
    run_framework(tmp_path, 3, mem_allocate=400)
    assert [(n, s["master_compute"]) for n, s in recording_cluster.calls] == [(3, True)]


def test_local_runs_run_all_tasks_as_before(tmp_path, slurm_job, fits_in_framework):
    # without a cluster the master runs every task, also in a job
    slurm_job("4242", 1)
    log_fn, outs = run_framework(tmp_path, 3, mem_allocate=400, parse_cluster=False)
    assert sorted(read_log(log_fn)) == sorted(outs)


@posix_only
def test_config_values_decide_the_fit(tmp_path, slurm_job, recording_cluster):
    # the fit uses the call's values after the config file (memPerCPU, maxCPUNum)
    slurm_job("4242", 16)
    conf = tmp_path / "config.json"
    conf.write_text(json.dumps({"MemPerCPU": 4, "maxCPUNum": 4}))
    (tmp_path / "a").mkdir()
    run_framework(tmp_path / "a", 2, cpus_per_task=8, mem_per_cpu=20.9, config_file=str(conf))
    assert [s["master_compute"] for _, s in recording_cluster.calls] == [True]  # 4 x 4 GB


# the same through the cluster loop and the fake Slurm (as the MATLAB tests run)


@posix_only
@pytest.mark.parametrize("n", [1, 3])
def test_master_in_a_job_runs_no_task_that_does_not_fit_it_cluster_loop(
    n, tmp_path, slurm_job, cluster_loop, capsys
):
    slurm_job("4242", 16)
    log_fn, _ = run_framework(tmp_path, n, mem_allocate=40)
    assert read_log(log_fn) == [], f"{n} task(s): the master ran a task."
    assert cluster_loop.call_count("sbatch") == n, f"{n} task(s) were not all submitted."
    assert "do not fit this job" in capsys.readouterr().out


@posix_only
def test_master_in_a_job_runs_tasks_that_fit_it_cluster_loop(tmp_path, slurm_job, cluster_loop):
    slurm_job("4242", 64)
    log_fn, outs = run_framework(tmp_path, 3, mem_allocate=40)
    assert read_log(log_fn) == [outs[2]]
    assert cluster_loop.call_count("sbatch") == 2


@posix_only
@pytest.mark.parametrize(("cpus", "master_runs", "jobs"), [(4, 0, 2), (1, 1, 1)])
def test_tasks_without_an_estimate_need_the_memory_of_their_jobs_cluster_loop(
    cpus, master_runs, jobs, tmp_path, slurm_job, cluster_loop
):
    slurm_job("4242", 16)
    log_fn, _ = run_framework(tmp_path, 2, cpus_per_task=cpus, mem_per_cpu=10)
    assert len(read_log(log_fn)) == master_runs, f"{cpus} CPUs per task."
    assert cluster_loop.call_count("sbatch") == jobs, f"{cpus} CPUs per task."


@posix_only
def test_master_outside_a_job_runs_tasks_as_before_cluster_loop(tmp_path, slurm_job, cluster_loop):
    slurm_job("", 1)
    log_fn, outs = run_framework(tmp_path, 3, mem_allocate=400)
    assert read_log(log_fn) == [outs[2]]
    assert cluster_loop.call_count("sbatch") == 2


# --- MasterComputeFitTest: oracle cases -----------------------------------------------------


def _with_job(job_id: str):
    # SLURM_JOB_ID set ('' unset) for the duration of a call; restored afterwards
    class _Job:
        def __enter__(self):
            self.orig = os.environ.get("SLURM_JOB_ID")
            if job_id:
                os.environ["SLURM_JOB_ID"] = job_id
            else:
                os.environ.pop("SLURM_JOB_ID", None)

        def __exit__(self, *exc):
            if self.orig is None:
                os.environ.pop("SLURM_JOB_ID", None)
            else:
                os.environ["SLURM_JOB_ID"] = self.orig

    return _Job()


def py_master_task_fits_grid(grid, job_id):
    """master_fit on each row (mem_allocate, cpus, mem_per_cpu, max, min, master_mem)"""
    tf, task_mem, master_mem = [], [], []
    with _with_job(job_id):
        for row in grid:
            fit = master_fit(*row)
            tf.append(fit.fits)
            task_mem.append(fit.task_mem)
            master_mem.append(fit.master_mem)
    return np.array(tf), np.array(task_mem), np.array(master_mem)


def py_coordinator_resources_grid(grid):
    """coordinator_resources on each row, with MATLAB's process allowance (2 GB)"""
    orig = coordinator.PROCESS_MEMORY
    coordinator.PROCESS_MEMORY = MATLAB_PROCESS_MEMORY
    try:
        out = [coordinator_resources(*row) for row in grid]
    finally:
        coordinator.PROCESS_MEMORY = orig
    return (np.array([float(c) for c, _ in out]), np.array([float(m) for _, m in out]))


# MATLAB: each row of the cell in{1} through the function (rows of numbers come as double rows)
_GRID_CODE = (
    "o = getenv('SLURM_JOB_ID'); if isempty(in{2}), unsetenv('SLURM_JOB_ID'); "
    "else, setenv('SLURM_JOB_ID', in{2}); end; g = in{1}; n = numel(g); "
    "tf = false(1, n); tm = zeros(1, n); mm = zeros(1, n); "
    "for k = 1 : n, c = g{k}; if ~iscell(c), c = num2cell(c); end; "
    "[tf(k), tm(k), mm(k)] = master_task_fits(c{:}); end; "
    "if isempty(o), unsetenv('SLURM_JOB_ID'); else, setenv('SLURM_JOB_ID', o); end; "
    "out = {tf, tm, mm};"
)
_RESOURCES_CODE = (
    "g = in{1}; n = numel(g); cpus = zeros(1, n); mem = zeros(1, n); "
    "for k = 1 : n, c = g{k}; if ~iscell(c), c = num2cell(c); end; "
    "[cpus(k), mem(k)] = coordinator_resources(c{:}); end; out = {cpus, mem};"
)


def fit_grid() -> list[list]:
    rows = []
    for mem in (None, 0.5, 15, 16.5, 16.8, 17, 100, NAN):
        for cpus in (1, 3, 30):
            for mem_per_cpu in (5, 20.9):
                for max_cpu in (None, 1, 24):
                    for min_cpu in (None, 2):
                        rows.append([mem, cpus, mem_per_cpu, max_cpu, min_cpu, 16])
    for master_mem in (NAN, 0, 15.9, 16, 1e6):
        for mem in (None, 16.8, 16.81, NAN):
            rows.append([mem, 3, 5.6, 24, 1, master_mem])
    return rows


def outside_grid() -> list[list]:
    return [[m, c, 20.9, x, 1, 1] for m in (None, 1000, NAN) for c in (1, 24) for x in (None, 24)]


def resources_grid() -> list[list]:
    return [
        [merge, computes, cpus, task_mem]
        for merge in (0, 0.25, 0.5, 2, 3, 30, NAN)
        for computes in (False, True)
        for cpus in (2, 24)
        for task_mem in (None, 0.5, 2.5, 40, NAN)
    ]


ORACLE_CASES = [
    Case(
        "master_task_fits_in_a_job",
        code=_GRID_CODE,
        nargout=3,
        args=[In("grid"), "4242"],
        python=py_master_task_fits_grid,
        inputs={"grid": Value(fit_grid)},
    ),
    Case(
        "master_task_fits_outside_a_job",
        code=_GRID_CODE,
        nargout=3,
        args=[In("grid"), ""],
        python=py_master_task_fits_grid,
        inputs={"grid": Value(outside_grid)},
    ),
    Case(
        "coordinator_resources_grid",
        code=_RESOURCES_CODE,
        nargout=2,
        args=[In("grid")],
        python=py_coordinator_resources_grid,
        inputs={"grid": Value(resources_grid)},
    ),
]


@oracle.parametrize(ORACLE_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- CoordinatorNiceTest --------------------------------------------------------------------


def stub_scontrol(monkeypatch, folder: Path, lines: list[str], status: int, delay: float = 0):
    """an scontrol first on the PATH that prints lines and exits with status"""
    folder.mkdir(parents=True, exist_ok=True)
    fn = folder / "scontrol"
    body = "".join(f'echo "{line}"\n' for line in lines)
    sleep = f"exec sleep {delay}\n" if delay else ""
    fn.write_text(f"#!/bin/bash\n{sleep}{body}exit {status}\n")
    fn.chmod(0o755)
    monkeypatch.setenv("PATH", f"{folder}{os.pathsep}{os.environ.get('PATH', '')}")


@posix_only
@pytest.mark.parametrize(
    ("ptype", "weights", "nice"),
    [
        ("priority/multifactor", ["PriorityWeightAge       = 1000", "PriorityWeightJobSize   = 0"],
         1000),
        ("priority/multifactor", ["PriorityWeightAge = 1000", "PriorityWeightJobSize = 500"], 1500),
        ("priority/basic", ["PriorityWeightAge = 1000"], 0),
        ("priority/multifactor", ["PriorityWeightAge = 1000", "PriorityWeightQOS = 7"], 1000),
        ("priority/multifactor", ["PriorityWeightFairShare = 1000"], 0),
    ],
)  # fmt: skip
def test_nice_is_the_weights_of_age_and_size(ptype, weights, nice, tmp_path, monkeypatch):
    # priority type and weights in `scontrol show config`, nice
    lines = ["PriorityDecayHalfLife   = 1-00:00:00", f"PriorityType            = {ptype}"]
    stub_scontrol(monkeypatch, tmp_path, [*lines, *weights, "PriorityWeightFairShare = 100000"], 0)
    assert coordinator_nice() == nice


@posix_only
def test_nice_is_0_without_slurm(tmp_path, monkeypatch):
    # (scontrol fails: no Slurm)
    stub_scontrol(monkeypatch, tmp_path, [], 1)
    assert coordinator_nice() == 0
    # (scontrol fails after printing the weights)
    stub_scontrol(
        monkeypatch,
        tmp_path / "b",
        ["PriorityType = priority/multifactor", "PriorityWeightAge = 9"],
        1,
    )
    assert coordinator_nice() == 0


@posix_only
def test_nice_is_0_when_scontrol_does_not_answer(tmp_path, monkeypatch):
    # (MATLAB would wait; a timed-out query is unknown: no nice, and no wait for the stub's
    # 600 s: the bound keeps a busy machine's slowness apart from that)
    from petakit5d.jobs.slurm import Slurm

    stub_scontrol(
        monkeypatch,
        tmp_path,
        ["PriorityType = priority/multifactor", "PriorityWeightAge = 9"],
        0,
        600,
    )
    t0 = time.monotonic()
    assert coordinator_nice(Slurm(query_timeout=0.5)) == 0
    assert time.monotonic() - t0 < 60


@posix_only
def test_framework_jobs_get_their_nice(tmp_path, slurm_job, recording_cluster):
    # the framework's nice option lowers the priority of its jobs by it (sbatch --nice)
    slurm_job("", 1)
    (tmp_path / "a").mkdir()
    run_framework(tmp_path / "a", 2, master_compute=False, nice=1000)
    (tmp_path / "b").mkdir()
    run_framework(tmp_path / "b", 2, master_compute=False)
    with_nice, without = (s["slurm_param"] for _, s in recording_cluster.calls)
    assert with_nice.endswith(" --nice=1000")
    assert "--nice" not in without, "A job without nice got one."


@posix_only
def test_framework_jobs_get_their_nice_cluster_loop(tmp_path, slurm_job, cluster_loop):
    # (the jobs don't run: the framework gives up after its try)
    slurm_job("", 1)
    run_framework(tmp_path, 2, master_compute=False, nice=1000, unit_wait_time=0)
    calls = cluster_loop.calls("sbatch")
    assert calls
    assert all(" --nice=1000" in c for c in calls), "A job of the call went without its nice."


# --- the real cluster: a coordinating job and its block jobs ---------------------------------

BLOCK_JOBS = 5
CHILD_DELAY = 600  # s; block jobs wait this long (--begin), so their priorities can be read


def coordinator_job_task(scratch: str, child_params: str, child_time: str) -> None:
    """the task of the coordinating job: submits BLOCK_JOBS one-CPU block jobs as the job
    framework submits them (one array job per batch), writes their ids and its own view of its
    job (memory, CPUs) to <scratch>/children.json, waits for <scratch>/release (at most 10
    minutes) and writes its peak memory to <scratch>/coordinator_done.json"""
    import resource

    from petakit5d.core import memory
    from petakit5d.jobs.slurm import Slurm

    slurm = Slurm()
    ids = []
    for b in range(1, BLOCK_JOBS + 1):
        r = slurm.submit(
            array_id=b,
            out=os.path.join(scratch, "block_%A_%a.out"),
            err=os.path.join(scratch, "block_%A_%a.err"),
            cpus_per_task=1,
            wrap="sleep 1",
            slurm_param=f"{child_params} --begin=now+{CHILD_DELAY}",
            time_str=child_time,
        )
        if r.job_id is not None:
            ids.append(f"{r.job_id}_{b}")
        if not r.submitted:
            break
    info = {
        "children": ids,
        "memory": memory.get_system_memory(),
        "cpus": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "job": os.environ.get("SLURM_JOB_ID"),
    }
    tmp = os.path.join(scratch, f"children.json.{os.getpid()}")
    with open(tmp, "w") as fh:
        json.dump(info, fh)
    os.replace(tmp, os.path.join(scratch, "children.json"))
    t0 = time.monotonic()
    while time.monotonic() - t0 < 600 and not os.path.exists(os.path.join(scratch, "release")):
        os.listdir(scratch)
        time.sleep(1)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
    with open(os.path.join(scratch, "coordinator_done.json"), "w") as fh:
        json.dump({"peak_rss_gb": peak}, fh)


def _fields(text: str | None) -> dict[str, str]:
    out = {}
    for word in (text or "").split():
        if "=" in word:
            k, v = word.split("=", 1)
            out.setdefault(k, v)
    return out


def _wait_for(path: Path, seconds: float) -> bool:
    t0 = time.monotonic()
    while time.monotonic() - t0 < seconds:
        os.listdir(path.parent)
        if path.exists():
            return True
        time.sleep(2)
    return False


def _gb(text: str) -> float:
    units = {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1, "T": 1024}
    return float(text[:-1]) * units[text[-1]] if text and text[-1] in units else float(text)


@pytest.mark.cluster
@pytest.mark.slow
def test_coordinator_on_the_cluster(tmp_path):
    """A coordinating job sized by coordinator_resources (merges of one CPU's memory, without
    computing: 2 CPUs by memory) and submitted with coordinator_nice's nice gets that request,
    and the block jobs it submits outrank it.

    The cluster comes from the config file PETAKIT5D_TEST_SLURM_CONFIG (SlurmParam with the
    partition, account, QOS and --mem-per-cpu; memPerCPU, maxCPUNum, minCPUNum, jobTimeLimit,
    PythonLaunchStr); the shared scratch folder from PETAKIT5D_TEST_SLURM_SCRATCH (else
    pytest's temp folder, which the jobs must see). One coordinating job and BLOCK_JOBS block
    jobs at a time, all cancelled at the end by their ids."""
    from petakit5d.jobs import shapes
    from petakit5d.jobs.config import ClusterConfig, default_python_launch_str
    from petakit5d.jobs.slurm import Slurm
    from petakit5d.jobs.tasks import Call, task_record, write_call

    conf = ClusterConfig.load(os.environ["PETAKIT5D_TEST_SLURM_CONFIG"])
    slurm_param = conf.get("SlurmParam", "")
    mem_per_cpu = float(conf.get("memPerCPU", 20.9))
    max_cpu = float(conf.get("maxCPUNum", 24))
    min_cpu = float(conf.get("minCPUNum", 1))
    time_limit = conf.get("jobTimeLimit", 0.25)
    python = conf.get("PythonLaunchStr", "") or default_python_launch_str()
    base = os.environ.get("PETAKIT5D_TEST_SLURM_SCRATCH") or str(tmp_path)
    scratch = Path(base) / f"coordinator_{os.getpid()}_{int(time.time())}"
    scratch.mkdir(parents=True)

    merge_mem = mem_per_cpu  # with the process: more than one CPU's memory
    cpus_per_task, mem_allocate = coordinator_resources(merge_mem, False, 4, None)
    cpus = int(shapes.job_cpus(cpus_per_task, mem_allocate, mem_per_cpu, min_cpu, max_cpu, True))
    assert cpus == 2 or max_cpu < 2, (cpus_per_task, mem_allocate, cpus)
    nice = coordinator_nice()
    t_str = shapes.time_str(time_limit, slurm_param)

    task = Task.of(coordinator_job_task, str(scratch), slurm_param, t_str)
    call = Call(
        [str(scratch)], [str(scratch / "coordinator_done.json")], [task_record(task, "python")],
        {"task_batch_num": 1, "language": "python", "uuid": "coordinator-test"},
    )  # fmt: skip
    folder = write_call(call, str(scratch / "call"))
    slurm = Slurm()
    coord_param = f"{slurm_param} --nice={nice}" if nice > 0 else slurm_param
    r = slurm.submit(
        array_id=1,
        out=str(scratch / "coordinator_%A_%a.out"),
        err=str(scratch / "coordinator_%A_%a.err"),
        cpus_per_task=cpus,
        wrap=f"{python} -m petakit5d.jobs.exec {folder} 1",
        slurm_param=coord_param,
        time_str=t_str,
    )
    assert r.submitted, r.result.output
    coord = f"{r.job_id}_1"
    children: list[str] = []
    report = {"coordinator": coord, "cpus": cpus, "mem_allocate": mem_allocate, "nice": nice}
    try:
        report["coordinator_sprio"] = slurm.run(["sprio", "-l", "-j", str(r.job_id)]).output
        job = _fields(slurm.show_job(coord))
        report["coordinator_show_job"] = {
            k: job.get(k) for k in ("JobState", "NumCPUs", "MinMemoryCPU", "Nice", "Priority")
        }
        assert job.get("NumCPUs") == str(cpus)
        assert job.get("Nice") == str(nice)
        assert _gb(job.get("MinMemoryCPU", "0")) * cpus >= mem_allocate

        assert _wait_for(scratch / "children.json", 900), "The coordinating job didn't start."
        info = json.loads((scratch / "children.json").read_text())
        children = info["children"]
        report["coordinator_view"] = info
        assert len(children) == BLOCK_JOBS
        assert info["memory"] >= mem_allocate  # (its cgroup: what master_task_fits sees)

        job = _fields(slurm.show_job(coord))
        coord_priority = int(job["Priority"])
        report["coordinator_priority"] = coord_priority
        prios = {}
        for c in children:
            cj = _fields(slurm.show_job(c))
            assert cj.get("Nice") == "0" and cj.get("JobState") == "PENDING", cj
            prios[c] = int(cj["Priority"])
        report["block_priorities"] = prios
        sprio = slurm.run(["sprio", "-l", "-j", ",".join(c.split("_")[0] for c in children)])
        report["sprio"] = sprio.output.splitlines()
        squeue = slurm.run(
            ["squeue", "-h", "-j", ",".join([coord, *children]), "-o", "%i %Q %y %T"]
        )
        report["squeue"] = squeue.output.splitlines()
        assert all(p > coord_priority for p in prios.values()), report
    finally:
        (scratch / "release").touch()
        if children:
            slurm.cancel(children)
        _wait_for(scratch / "coordinator_done.json", 120)
        if slurm.job_state(coord):  # (still listed: cancel it)
            slurm.cancel([coord])
        report_fn = scratch / "report.json"
        if (scratch / "coordinator_done.json").exists():
            report["coordinator_done"] = json.loads((scratch / "coordinator_done.json").read_text())
        report_fn.write_text(json.dumps(report, indent=1))
        print(json.dumps(report, indent=1))
    for _ in range(30):  # (cancelled jobs can be listed for a few seconds)
        left = slurm.live_jobs([coord, *children])
        if left is not None and not left:
            break
        time.sleep(2)
    assert left == set(), f"jobs of the test still in the queue: {left}"
