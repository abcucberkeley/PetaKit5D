"""Tests the mutation campaign asked for (E06, docs/hardening.md): the job framework."""

import math
import sys
import types
from types import SimpleNamespace

import pytest

from petakit5d.jobs import cluster, config, metrics, runner, slurm, tasks


def test_a_config_counts_its_keys():
    c = config.ClusterConfig.from_mapping({"maxCPUNum": 4, "MemPerCPU": 1, "SlurmParam": ""})
    assert len(c) == 3 and bool(c)
    assert len(config.ClusterConfig.empty()) == 0


def test_a_missing_status_field_is_nan():
    text = "Name:\tpython\nVmRSS:\t  2048 kB\nThreads:\t5\n"
    assert metrics._field(text, "VmRSS:") == 2048
    assert metrics._field(text, "Threads:") == 5
    assert math.isnan(metrics._field(text, "VmHWM:"))
    assert math.isnan(metrics._field("VmHWM:\tnone", "VmHWM:"))
    # (a missing name is no position: not the number at its length from the start)
    assert math.isnan(metrics._field("Name:  42\n", "Threads:"))


class _Slurm:
    def __init__(self, stdout=""):
        self.argv = []
        self.stdout = stdout

    def run(self, argv, *a, **k):
        self.argv.append(list(argv))
        return SimpleNamespace(ok=True, stdout=self.stdout, output=self.stdout)


def test_job_states_are_asked_1000_jobs_at_a_time():
    r = runner.ClusterRunners.__new__(runner.ClusterRunners)
    r._slurm = _Slurm()
    keys = [f"{100 + k}_1" for k in range(2500)]
    assert r._states_of(keys) == {}
    sizes = [len(a[a.index("-j") + 1].split(",")) for a in r._slurm.argv]
    assert sizes == [1000, 1000, 500]


def test_pending_reasons_are_of_the_login_user_without_user(monkeypatch):
    s = slurm.Slurm()
    fake = _Slurm("12_3 Priority\n")
    monkeypatch.setattr(s, "run", fake.run)
    monkeypatch.setattr(slurm, "_login_name", lambda: "login_user")
    monkeypatch.delenv("USER", raising=False)
    assert s.pending_reasons() == {"12_3": "Priority"}
    monkeypatch.setenv("USER", "env_user")
    s.pending_reasons()
    users = [a[a.index("-u") + 1] for a in fake.argv]
    assert users == ["login_user", "env_user"]


def test_cancelling_finds_the_jobs_of_lost_submissions_by_their_comment():
    loop = cluster._ClusterLoop.__new__(cluster._ClusterLoop)
    loop.unknown = {2: "mk_b2", 3: "mk_b3"}
    loop.job_ids = [11, 0, 0, 0]
    loop.cancelled = False
    loop._batch = lambda b: range(b - 1, b)
    loop._lookup = lambda: {"mk_b2": (42, 1)}
    cancelled = []
    loop._scancel = lambda keys: cancelled.extend(keys)
    loop._cancel_all()
    assert loop.cancelled and loop.unknown == {}
    assert loop.job_ids == [11, 42, 0, 0]
    assert cancelled == ["11", "42"]


def main_task():
    """(a module-level function that the next test presents as one of the script run as
    __main__)"""


def _as_main(monkeypatch, file=None, spec_name=None):
    main = types.ModuleType("__main__")
    if file is not None:
        main.__file__ = str(file)
    main.__spec__ = SimpleNamespace(name=spec_name) if spec_name else None
    monkeypatch.setitem(sys.modules, "__main__", main)
    monkeypatch.setattr(main_task, "__module__", "__main__")


def test_a_function_of_the_main_script_is_imported_by_the_scripts_name(monkeypatch, tmp_path):
    _as_main(monkeypatch, tmp_path / "my_script.py")
    t = tasks.Task.of(main_task)
    assert t.func == "my_script:main_task" and t.path == [str(tmp_path)]
    # (python -m pkg.mod: the module's name, and the folder its package is in)
    _as_main(monkeypatch, tmp_path / "pkg" / "mod.py", "pkg.mod")
    t = tasks.Task.of(main_task)
    assert t.func == "pkg.mod:main_task" and t.path == [str(tmp_path)]
    _as_main(monkeypatch, tmp_path / "my-script.py")
    with pytest.raises(ValueError, match="has no importable module name"):
        tasks.Task.of(main_task)
    _as_main(monkeypatch)
    with pytest.raises(ValueError, match="defined interactively"):
        tasks.Task.of(main_task)


def test_a_failing_task_of_a_job_is_reported_by_its_index_and_function(tmp_path):
    # (a job without task metrics: the error is the task's, named, and the job's status is 1)
    import os
    import subprocess
    from pathlib import Path

    import petakit5d

    pkg = tmp_path / "lib" / "hardtasks"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "jobs.py").write_text("def fail():\n    raise RuntimeError('on purpose')\n")
    task = tasks.Task("hardtasks.jobs:fail", [], path=[str(tmp_path / "lib")])
    call = tasks.Call(["a"], ["b"], [tasks.task_record(task, "python")], {"task_batch_num": 1})
    folder = tasks.write_call(call, tmp_path / "c")
    env = dict(os.environ, PYTHONPATH=str(Path(petakit5d.__file__).resolve().parents[1]))
    cmd = [sys.executable, "-m", "petakit5d.jobs.exec", str(folder), "1"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
    assert r.returncode == 1
    assert "RuntimeError: on purpose" in r.stderr
    assert "(hardtasks.jobs:fail):" in r.stderr and "Error in task" in r.stderr
    assert "AttributeError" not in r.stderr
