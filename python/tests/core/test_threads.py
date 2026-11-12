import os

import pytest

from petakit5d.core import threads


def _tree(tmp_path, files: dict):
    for rel, content in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content)
    return str(tmp_path)


def test_cgroup_v2_quota(tmp_path):
    root = _tree(
        tmp_path,
        {
            "a/cpu.max": "max 100000\n",
            "a/job_1/cpu.max": "200000 100000\n",
            "a/job_1/step/cpu.max": "350000 100000\n",
        },
    )
    assert threads.cgroup_cpu_limit("0::/a/job_1/step\n", root) == 2.0
    assert threads.cgroup_cpu_limit("0::/a\n", root) is None
    assert threads.cgroup_cpu_limit("garbage", root) is None


def test_cgroup_v1_quota(tmp_path):
    root = _tree(
        tmp_path,
        {
            "cpu,cpuacct/slurm/job_1/cpu.cfs_quota_us": "150000",
            "cpu,cpuacct/slurm/job_1/cpu.cfs_period_us": "100000",
            "cpu,cpuacct/slurm/cpu.cfs_quota_us": "-1",
            "cpu,cpuacct/slurm/cpu.cfs_period_us": "100000",
        },
    )
    assert threads.cgroup_cpu_limit("4:cpu,cpuacct:/slurm/job_1\n", root) == 1.5


def test_available_and_resolved_threads(monkeypatch):
    n = threads.available_cpus()
    assert 1 <= n <= (os.cpu_count() or 1)
    assert threads.resolve_threads() == n == threads.resolve_threads(0)
    assert threads.resolve_threads(1) == 1
    assert threads.resolve_threads(10**6) == n
    with pytest.raises(ValueError):
        threads.resolve_threads(-2)
    monkeypatch.setattr(threads, "_affinity_count", lambda: 8)
    monkeypatch.setattr(threads, "cgroup_cpu_limit", lambda: 1.5)
    assert threads.available_cpus() == 2
    monkeypatch.setattr(threads, "cgroup_cpu_limit", lambda: None)
    assert threads.available_cpus() == 8


def test_agrees_with_the_compiled_count():
    from petakit5d import _ext

    if not _ext.available("core"):
        pytest.skip("core extension not built")
    core = _ext.load("core")
    assert threads._affinity_count() == core.affinity_cpu_count()


def test_no_global_thread_state(monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    # the environment is not a source of thread counts
    assert threads.available_cpus() == threads.resolve_threads(None)
    assert os.environ["OMP_NUM_THREADS"] == "1"


_BOUND_CHILD = """
import os, petakit5d.core.threads as t
from petakit5d import _ext
pre = t.available_cpus()
mask = len(os.sched_getaffinity(0))
_ext.load("core").parallel_threads(max(2, mask))  # libgomp binds the main thread here
print(pre, mask, len(os.sched_getaffinity(0)), t.available_cpus(), t.resolve_threads())
"""


def test_omp_proc_bind_does_not_shrink_the_cpus():
    # K03: with OMP_PROC_BIND set, libgomp binds the main thread to one CPU at its first parallel
    # region; the affinity read after that is 1 CPU, and every kernel ran single-threaded
    import subprocess
    import sys

    from petakit5d import _ext

    if not hasattr(os, "sched_getaffinity") or not _ext.available("core"):
        pytest.skip("needs Linux affinity masks and the core extension")
    if len(os.sched_getaffinity(0)) < 2:
        pytest.skip("needs 2 or more CPUs")
    env = dict(os.environ, OMP_PROC_BIND="true")
    env.pop("OMP_NUM_THREADS", None)
    out = subprocess.run(
        [sys.executable, "-c", _BOUND_CHILD], env=env, capture_output=True, text=True, check=True
    ).stdout
    pre, mask, bound, after, resolved = (int(v) for v in out.split())
    if bound == mask:
        pytest.skip("this OpenMP runtime did not bind the main thread")
    assert bound == 1  # what the old code took as the process's CPUs
    assert after == pre == resolved and pre > 1
