"""Task metrics of the job framework (petakit5d.jobs.metrics). Ported from
tests/framework/TaskMetricsTest.m; the test names keep the MATLAB names in snake_case. The
compiled-environment tests (MATLAB Runtime task lines, GNU parallel joblogs written by jobs,
mccMaster) have no counterpart: Python jobs have no compiled mode, no GNU parallel and no task
line environment. Two oracle cases check that MATLAB's read_task_metrics reads run folders
that Python wrote: one written by run_tasks itself (clocks and /proc made deterministic), one
by hand with every kind of value the reader meets.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import textwrap
import time
import types
from pathlib import Path

import pytest

import petakit5d
from petakit5d.jobs import framework, metrics
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.metrics import TaskMetricsWarning, read_task_metrics
from petakit5d.jobs.tasks import Call, CallFolder, Task, run_record, write_call
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, InPath, Tree

# importing this module at collection must not hide a stand-in metrics module that other
# tests put in sys.modules ("from petakit5d.jobs import metrics" takes the package attribute
# first); the module stays in sys.modules
del sys.modules["petakit5d.jobs"].metrics
SRC_ROOT = str(Path(petakit5d.__file__).resolve().parents[1])
LOCAL = {"parse_cluster": False, "master_compute": True}
LINUX = sys.platform.startswith("linux")


class TaskFailed(Exception):
    """An error with MATLAB's identifier."""

    identifier = "pktest:frameworkTask:failed"


def framework_task(input_fn: str, output_fn: str, log_fn: str, mode: str = "copy") -> None:
    """Appends output_fn to log_fn, then 'copy' the input to the output or 'error' raise."""
    with open(log_fn, "a") as fh:
        fh.write(output_fn + "\n")
    if mode == "copy":
        Path(output_fn).write_bytes(Path(input_fn).read_bytes())
    elif mode == "error":
        raise TaskFailed(f"Task failed on purpose for {input_fn}")
    else:
        raise ValueError(mode)


def alloc_task(output_fn: str, n: int) -> None:
    """Fills n bytes (every page touched), then writes the output."""
    b = b"\x01" * n
    Path(output_fn).write_text(str(len(b)), newline="\n")


def compute_task(output_fn: str) -> None:
    """About half a second of arithmetic on small values (user time)."""
    t0 = time.process_time()
    s = 0.0
    while time.process_time() - t0 < 0.5:
        s += sum(math.sqrt(k) for k in range(10000))
    Path(output_fn).write_text(repr(s), newline="\n")


def nested_task(input_fn: str, output_fn: str, log_fn: str, metrics_on: bool, n: int = 2) -> None:
    """A task that runs framework calls of its own, one after another (as large-file tasks
    do): the last one writes output_fn. (The inner calls skip the wait for a fresh .tmp
    marker: the last one has the outer task's output, and so its marker.)"""
    for k in range(1, n + 1):
        out = output_fn if k == n else f"{output_fn}.{k}"
        run_tasks(
            [input_fn],
            [out],
            [Task(FUNC, [input_fn, out, log_fn, "copy"])],
            task_metrics=metrics_on,
            unit_wait_time=0,
            **LOCAL,
        )


FUNC = f"{__name__}:framework_task"


def make_tasks(tmp_path: Path, n: int):
    (tmp_path / "out").mkdir()
    in_fns = [str(tmp_path / f"in_{i}.txt") for i in range(1, n + 1)]
    out_fns = [str(tmp_path / "out" / f"out_{i}.txt") for i in range(1, n + 1)]
    for i, fn in enumerate(in_fns, 1):
        Path(fn).write_text(f"task {i}", newline="\n")
    return in_fns, out_fns, str(tmp_path / "calls.log")


def task_list(in_fns, out_fns, log_fn, mode="copy"):
    return [Task(FUNC, [i, o, log_fn, mode]) for i, o in zip(in_fns, out_fns, strict=True)]


def read_log(log_fn: str) -> list[str]:
    if not os.path.exists(log_fn):
        return []
    return [x for x in Path(log_fn).read_text().splitlines() if x]


def by_task(rows):
    return sorted(rows, key=lambda r: r["task"])


def col(rows, name):
    return [r[name] for r in rows]


def write_lines(fn: Path, lines) -> None:
    fn.write_text("".join(x + "\n" for x in lines), encoding="utf-8", newline="\n")


def the_run_folder(job_log_dir: Path) -> Path:
    runs = [p for p in (job_log_dir / "task_metrics").iterdir() if p.is_dir()]
    assert len(runs) == 1, "There should be one metrics folder."
    return runs[0]


def job_call(folder: str, n: int = 1, settings=None) -> Call:
    # a call whose tasks record into folder, as a job's call folder has it
    s = {"task_metrics": folder, "task_batch_num": 1, "language": "python"}
    s.update(settings or {})
    return Call(
        [f"/in_{i}" for i in range(n)],
        [f"/out_{i}" for i in range(n)],
        [{"func": "shutil:copyfile", "args": []}] * n,
        s,
    )


@pytest.fixture(autouse=True)
def fresh_metrics(monkeypatch):
    """task_metrics keeps this process's file name, the current task, the single-job folders
    and its one warning (MATLAB: clear task_metrics task_metrics_run); the metrics
    environment of the MATLAB suite is cleared."""
    monkeypatch.setattr(metrics, "_process", None)
    monkeypatch.setattr(metrics, "_current", "")
    monkeypatch.setattr(metrics, "_warned", False)
    monkeypatch.setattr(metrics, "_sessions", {})
    monkeypatch.setattr(metrics, "_run_t0", {})
    for name in ("SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "PARALLEL_JOBSLOT", "PARALLEL_SEQ"):
        monkeypatch.delenv(name, raising=False)


# --- TaskMetricsTest ---------------------------------------------------------------------------


def test_local_tasks_are_recorded(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    t0 = time.time()
    flags = run_tasks(
        in_fns, out_fns, task_list(in_fns, out_fns, log_fn), task_metrics=True, **LOCAL
    )
    t1 = time.time()
    assert flags == [True] * 3
    assert [Path(o).read_text() for o in out_fns] == [Path(i).read_text() for i in in_fns]
    # without a job log folder the metrics go to <data>/job_logs
    t, runs = read_task_metrics(str(tmp_path / "job_logs"))
    assert len(t) == 3
    t = by_task(t)
    assert col(t, "task") == [1, 2, 3]
    assert col(t, "func") == [FUNC] * 3
    assert col(t, "mode") == ["local"] * 3
    assert col(t, "status") == ["ok"] * 3
    assert col(t, "exitCode") == [0] * 3
    assert col(t, "output") == out_fns
    assert col(t, "input") == in_fns
    assert col(t, "pid") == [os.getpid()] * 3
    assert all(col(t, "node")), "No node name."
    for r in t:
        assert t0 - 1 <= r["start"] <= r["end"] <= t1 + 1, "Task times are not seconds since 1970."
        assert abs(r["end"] - r["start"] - r["wallSec"]) < 1e-3
        assert r["cpuSec"] >= 0
        assert r["threads"] >= 1
        assert r["metricsSec"] < 0.5, "The metrics took longer than expected."
    assert len(runs) == 1
    run = runs[0]
    assert (run["nTasks"], run["nRecords"]) == (3, 3)
    assert run["backend"] == "python"
    assert run["func"] == FUNC
    assert run["parseCluster"] == 0
    # the master's Slurm job (if any) tells the tasks it runs itself
    folder = Path(run["folder"])
    rec = json.loads((folder / "run.ndjson").read_text().splitlines()[0])
    assert rec["jobId"] == os.environ.get("SLURM_JOB_ID", "")
    assert rec["kind"] == "run_start"
    assert rec["taskMetrics"] is True
    # strict JSON for other readers: null, not NaN, for missing values
    text = next(folder.glob("task_*.ndjson")).read_text()
    assert "NaN" not in text and "Infinity" not in text
    for line in text.splitlines():
        json.loads(line, parse_constant=metrics._no_constant)


def test_local_run_records_its_end(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn), task_metrics=True, **LOCAL)
    _, runs = read_task_metrics(str(tmp_path / "job_logs"))
    assert (runs[0]["nTasks"], runs[0]["nDone"]) == (3, 3)


def test_metrics_are_off_by_default(tmp_path, recwarn):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn), **LOCAL)
    assert flags == [True] * 2
    assert not [w for w in recwarn if issubclass(w.category, TaskMetricsWarning)]
    assert not (tmp_path / "job_logs").exists(), "A metrics folder was made."


def test_metrics_off_imports_nothing(tmp_path):
    # without metrics the metrics module isn't even imported (no cost), in a fresh process
    in_fns, out_fns, _ = make_tasks(tmp_path, 2)
    code = textwrap.dedent(f"""
        import sys
        from petakit5d.jobs.framework import run_tasks
        from petakit5d.jobs.tasks import Task
        ins, outs = {in_fns!r}, {out_fns!r}
        run_tasks(ins, outs, [Task("shutil:copyfile", [i, o]) for i, o in zip(ins, outs)],
                  parse_cluster=False)
        print("petakit5d.jobs.metrics" in sys.modules)
    """)
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "False"


def test_config_file_turns_metrics_on(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"taskMetrics": True}), newline="\n")
    run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn), config_file=str(config), **LOCAL)
    t, _ = read_task_metrics(str(tmp_path / "job_logs"))
    assert sorted(col(t, "task")) == [1, 2]


def test_failed_local_task_is_recorded_and_still_fails(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    with pytest.raises(TaskFailed):
        run_tasks(
            in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "error"), task_metrics=True, **LOCAL
        )
    t, _ = read_task_metrics(str(tmp_path / "job_logs"))
    assert len(t) == 1
    assert t[0]["status"] == "error"
    assert t[0]["exitCode"] == 1
    assert "pktest:frameworkTask:failed" in t[0]["error"]
    # the task ran in try and again in except, in one record
    assert len(read_log(log_fn)) == 2


def test_failed_task_without_identifier_names_its_class(tmp_path):
    in_fns, out_fns, _ = make_tasks(tmp_path, 1)
    missing = str(tmp_path / "missing.txt")
    with pytest.raises(FileNotFoundError):
        run_tasks(in_fns, out_fns, [Task("os:remove", [missing])], task_metrics=True, **LOCAL)
    t, _ = read_task_metrics(str(tmp_path / "job_logs"))
    assert t[0]["error"].startswith("FileNotFoundError: ")
    assert repr(missing)[1:-1] in t[0]["error"]  # the OSError's text quotes it (\ as \\)


def test_bash_task_exit_status_is_recorded(tmp_path):
    in_fns, out_fns, _ = make_tasks(tmp_path, 1)
    run_tasks(
        in_fns,
        out_fns,
        ["exit 3"],
        language="bash",
        max_trial_num=1,
        unit_wait_time=0,
        task_metrics=True,
        **LOCAL,
    )
    t, _ = read_task_metrics(str(tmp_path / "job_logs"))
    assert len(t) == 1
    assert (t[0]["status"], t[0]["func"]) == ("failed", "exit 3")
    assert t[0]["exitCode"] == 3


@pytest.mark.skipif(not LINUX, reason="Resource use beyond time is recorded on Linux only.")
def test_memory_and_cpu_on_linux(tmp_path):
    # the peak resident memory is the task's own: it is reset at the start of each task (a
    # task that fills 512 MB, then one that doesn't); user and system CPU add up to the CPU
    # time; a task that computes on small values is mostly user time
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    tasks = [
        Task(f"{__name__}:alloc_task", [out_fns[0], 512 * 2**20]),
        Task(FUNC, [in_fns[1], out_fns[1], log_fn, "copy"]),
        Task(f"{__name__}:compute_task", [out_fns[2]]),
    ]
    run_tasks(in_fns, out_fns, tasks, max_trial_num=1, unit_wait_time=0, task_metrics=True, **LOCAL)
    t = by_task(read_task_metrics(str(tmp_path / "job_logs"))[0])
    assert len(t) == 3
    assert t[0]["rssPeakGB"] - t[0]["rssStartGB"] > 0.45, "The 512 MB are not in the peak."
    assert t[1]["rssPeakGB"] - t[1]["rssStartGB"] < 0.2, (
        "The peak was not reset for the second task."
    )
    for r in t:
        assert abs(r["userSec"] + r["sysSec"] - r["cpuSec"]) < 0.03 + 0.05 * r["cpuSec"]
    assert t[2]["userSec"] > 0.5 * t[2]["cpuSec"], "The computing task should be mostly user time."
    # (transparent huge pages fault 2 MB at a time)
    assert t[0]["minorFaults"] > 200, "The page faults of filling 512 MB are missing."
    assert all(r["cpus"] >= 1 and r["osThreads"] >= 1 and r["memTotalGB"] > 0 for r in t)
    assert all(r["processStart"] < r["start"] for r in t)
    assert all(r["rchar"] >= 0 and r["wchar"] > 0 for r in t)


def test_unwritable_metrics_do_not_fail_tasks(tmp_path):
    # a file where the metrics folder would be: no metrics, but the tasks run as usual
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    (tmp_path / "job_logs").write_text("", newline="\n")
    with pytest.warns(TaskMetricsWarning) as rec:
        flags = run_tasks(
            in_fns, out_fns, task_list(in_fns, out_fns, log_fn), task_metrics=True, **LOCAL
        )
    assert flags == [True] * 2
    assert sorted(read_log(log_fn)) == sorted(out_fns)
    assert rec[0].message.identifier == "task_metrics_run:failed"


def test_failed_record_write_only_warns(tmp_path):
    (tmp_path / "notAFolder").write_text("", newline="\n")
    tm = metrics.task_start(job_call(str(tmp_path / "notAFolder" / "run") + "/"), [0], "local")
    assert tm is not None
    with pytest.warns(TaskMetricsWarning) as rec:
        metrics.task_stop(tm)
    assert rec[0].message.identifier == "task_metrics:failed"
    # one warning per process
    tm = metrics.task_start(job_call(str(tmp_path / "notAFolder" / "run") + "/"), [0], "local")
    with warnings_as_list() as ws:
        metrics.task_stop(tm)
    assert not ws


def test_task_start_never_raises():
    with pytest.warns(TaskMetricsWarning):
        assert metrics.task_start(job_call("/x/run/"), [5], "local") is None  # no task 6
    metrics.task_stop(None, RuntimeError("x"))


def test_job_tasks_are_recorded(tmp_path):
    # (MATLAB: matlabJobsRunEachTaskWithMetrics) the master makes the run folder and writes it
    # into the call; a job reads the call folder and records each task it runs
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    call = Call(
        in_fns,
        out_fns,
        [t.record() for t in task_list(in_fns, out_fns, log_fn)],
        {
            "task_metrics": True,
            "job_log_dir": str(tmp_path),
            "task_batch_num": 2,
            "parse_cluster": True,
            "language": "python",
        },
    )
    folder = metrics.run_start(call)
    assert folder and call.settings["task_metrics"] == folder
    cf = CallFolder(write_call(call, str(tmp_path / "call")))
    assert cf.meta["settings"]["task_metrics"] == folder
    for i in cf.batch(1):
        tm = metrics.task_start(cf, [i], "slurm_python")
        run_record(cf.lines([i])[0]["task"])
        metrics.task_stop(tm, None)
    assert [Path(o).read_text() for o in out_fns] == [Path(i).read_text() for i in in_fns]
    t, runs = read_task_metrics(str(tmp_path / "task_metrics"))
    assert sorted(col(t, "task")) == [1, 2]
    assert col(t, "mode") == ["slurm_python"] * 2
    assert col(t, "status") == ["ok"] * 2
    assert col(t, "batchSize") == [1, 1]
    assert sorted(col(t, "output")) == sorted(out_fns)
    assert runs[0]["backend"] == "python"
    assert runs[0]["parseCluster"] == 1
    # the first task of a job process is its start-up
    assert t[0]["startupSec"] > 0 or math.isnan(t[0]["processStart"])


def test_job_entry_point_records_tasks(tmp_path):
    in_fns, out_fns, _ = make_tasks(tmp_path, 1)
    call = Call(
        in_fns,
        out_fns,
        [Task("shutil:copyfile", [in_fns[0], out_fns[0]]).record()],
        {"task_metrics": True, "job_log_dir": str(tmp_path), "language": "python"},
    )
    metrics.run_start(call)
    folder = write_call(call, str(tmp_path / "call"))
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    subprocess.run([sys.executable, "-m", "petakit5d.jobs.exec", folder, "1"], env=env, check=False)
    t, _ = read_task_metrics(str(tmp_path / "task_metrics"))
    assert col(t, "mode") == ["slurm_python"]


def test_master_computed_tasks_are_recorded_as_master(tmp_path):
    # in cluster mode with master compute the master runs the last batch itself (here task 2)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    call = Call(
        in_fns,
        out_fns,
        [t.record() for t in task_list(in_fns, out_fns, log_fn)],
        {
            "task_metrics": True,
            "job_log_dir": str(tmp_path),
            "language": "python",
            "parse_cluster": True,
        },
    )
    # (no run_start: the master's first task makes the run folder)
    framework.run_batch_in_process(call, [1], ["2"], mode="master")
    t, _ = read_task_metrics(str(tmp_path / "task_metrics"))
    assert len(t) == 1
    assert (t[0]["task"], t[0]["batchSize"], t[0]["mode"]) == (2, 1, "master")


def test_tasks_run_without_the_metrics_module(tmp_path, monkeypatch, capsys):
    # (MATLAB: matlabJobStatementRunsWithoutTaskMetrics) a process that can't import the
    # metrics still runs its tasks
    monkeypatch.delattr(sys.modules["petakit5d.jobs"], "metrics", raising=False)
    monkeypatch.setitem(sys.modules, "petakit5d.jobs.metrics", None)
    monkeypatch.setattr(framework, "_metrics_warned", False)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    flags = run_tasks(
        in_fns, out_fns, task_list(in_fns, out_fns, log_fn), task_metrics=True, **LOCAL
    )
    assert flags == [True]
    assert Path(out_fns[0]).read_text() == Path(in_fns[0]).read_text()
    assert "task metrics are not available" in capsys.readouterr().err


def test_job_task_errors_are_recorded(tmp_path):
    # (MATLAB: matlabJobStatementRethrowsTaskErrors) a failed task fails its job as before,
    # after its record
    folder = tmp_path / "metrics"
    folder.mkdir()
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 4)
    call = Call(
        in_fns,
        out_fns,
        [t.record() for t in task_list(in_fns, out_fns, log_fn, "error")],
        {"task_metrics": str(folder) + "/"},
    )
    tm = metrics.task_start(call, [3], "slurm_python")
    with pytest.raises(TaskFailed):
        try:
            run_record(call.records[3])
        except Exception as err:
            metrics.task_stop(tm, err)
            raise
    t, _ = read_task_metrics(str(folder))
    assert [(r["status"], r["mode"], r["task"]) for r in t] == [("error", "slurm_python", 4)]


def test_call_without_metrics_writes_no_metrics(tmp_path):
    # (MATLAB: matlabJobCommandIsUnchangedWithoutMetrics) the call folder of a call without
    # metrics says so, and no metrics folder is made
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    call = Call(
        in_fns,
        out_fns,
        [t.record() for t in task_list(in_fns, out_fns, log_fn)],
        {"task_metrics": False, "job_log_dir": str(tmp_path)},
    )
    cf = CallFolder(write_call(call, str(tmp_path / "call")))
    assert cf.meta["settings"]["task_metrics"] is False
    assert metrics.task_start(cf, [0], "slurm_python") is None
    assert metrics.task_start(call, [0], "local") is None
    assert not (tmp_path / "task_metrics").exists()


def test_nested_call_without_metrics_records_nothing(tmp_path):
    # (MATLAB: callWithoutMetricsClearsTaskEnvironment, singleJobWithoutMetricsClearsTask-
    # Environment) a task with metrics runs a framework call without metrics: that call's tasks
    # don't record into the task's run
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    task = Task(f"{__name__}:nested_task", [in_fns[0], out_fns[0], log_fn, False, 1])
    run_tasks(in_fns, out_fns, [task], task_metrics=True, **LOCAL)
    t, runs = read_task_metrics(str(tmp_path / "job_logs"))
    assert len(runs) == 1
    assert col(t, "func") == [f"{__name__}:nested_task"]
    assert metrics._current == ""


def test_reader_labels_tasks_the_master_ran(tmp_path):
    # compiled MATLAB tasks in the master's own Slurm job ran in the master
    run = tmp_path / "run_2"
    run.mkdir()
    write_lines(
        run / "run.ndjson",
        [
            json.dumps(
                {
                    "kind": "run_start",
                    "runId": "run_2",
                    "time": 100,
                    "jobId": "77",
                    "parseCluster": True,
                    "outputFullpaths": ["/d/a", "/d/b"],
                }
            )
        ],
    )
    write_lines(
        run / "task_n1_1_1.ndjson",
        [
            json.dumps(
                {"kind": "task", "runId": "run_2", "task": 1, "mode": "slurm_mcc", "jobId": "77"}
            ),
            json.dumps(
                {"kind": "task", "runId": "run_2", "task": 2, "mode": "slurm_mcc", "jobId": "78"}
            ),
        ],
    )
    t = by_task(read_task_metrics(str(run))[0])
    assert col(t, "mode") == ["master", "slurm_mcc"]


def test_single_jobs_record_submission_and_task(tmp_path, monkeypatch):
    # a single job's submission (in the submitting process's session folder) and its task
    # (in the job, with Slurm's array job and task ids) are joined by the array job id
    folder = metrics.session_dir(str(tmp_path))
    assert folder and metrics.session_dir(str(tmp_path)) == folder, "One folder per process."
    func_str = 'shutil:copyfile("/a", "/b")'
    metrics.record(
        folder,
        "submit",
        {
            "task": 7,
            "jobId": 1000,
            "functionStr": func_str,
            "cpusPerTask": 2,
            "memAllocate": 30,
            "memPerCPU": 20.9,
            "mccMode": False,
            "SlurmParam": "",
        },
    )
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "1000")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "7")
    call = job_call(folder, 7)
    metrics.task_stop(metrics.task_start(call, [6], "slurm_python"))
    t, runs = read_task_metrics(str(tmp_path / "task_metrics"))
    assert len(t) == 1
    assert (t[0]["mode"], t[0]["status"], t[0]["arrayJobId"]) == ("slurm_python", "ok", "1000")
    assert t[0]["task"] == 7
    assert (t[0]["cpusPerTask"], t[0]["memAllocate"]) == (2, 30), "The submission was not joined."
    assert len(runs) == 1
    assert (runs[0]["kind"], runs[0]["func"]) == ("single", "shutil:copyfile")
    assert runs[0]["nTasks"] == 1


def test_session_dir_without_metrics_or_folder():
    assert metrics.session_dir("/no/such/folder") == ""
    assert metrics.session_dir("/tmp", task_metrics=False) == ""


def test_nested_run_records_its_parent(tmp_path):
    # a task that runs framework calls of its own: the inner calls record the outer task as
    # their parent
    in_fns, _, log_fn = make_tasks(tmp_path, 1)
    inner_out = str(tmp_path / "out" / "inner.txt")
    task = Task(f"{__name__}:nested_task", [in_fns[0], inner_out, log_fn, True])
    run_tasks(in_fns, [inner_out], [task], task_metrics=True, **LOCAL)
    t, runs = read_task_metrics(str(tmp_path / "job_logs"))
    assert len(runs) == 3
    outer = [r for r in runs if r["func"] == f"{__name__}:nested_task"]
    inner = [r for r in runs if r["func"] == FUNC]
    assert (len(outer), len(inner)) == (1, 2)
    assert outer[0]["parent"] == ""
    assert col(inner, "parent") == [outer[0]["runId"] + "/1"] * 2
    assert sorted(col(t, "func")) == sorted([FUNC, FUNC, f"{__name__}:nested_task"])
    assert [r["parent"] for r in t if r["func"] == FUNC] == [outer[0]["runId"] + "/1"] * 2
    assert [r["parent"] for r in t if r["func"] != FUNC] == [""]


def make_joblog_run(run: Path) -> None:
    """MATLAB's readerJoinsJoblogAndRunSettings folder: a compiled run with GNU parallel's
    joblog, a killed line and a line whose output existed."""
    run.mkdir(parents=True)
    write_lines(
        run / "run.ndjson",
        [
            json.dumps(
                {
                    "kind": "run_start",
                    "runId": "run_1",
                    "time": 100,
                    "elapsedSec": 1,
                    "backend": "mcc",
                    "cpusPerTask": 2,
                    "memAllocate": 5,
                    "mccMode": True,
                    "inputFullpaths": ["/d/in_1.tif", "/d/in_2.tif", "/d/in_3.tif"],
                    "outputFullpaths": ["/d/out_1.tif", "/d/out_2.tif", "/d/out_3.tif"],
                    "functionStrs": ["f(1)", "f(2)", "f(3)"],
                }
            ),
            json.dumps({"kind": "run_end", "time": 140, "elapsedSec": 41, "nTasks": 3, "nDone": 2}),
        ],
    )
    write_lines(
        run / "task_n1_11_5.ndjson",
        [
            json.dumps(
                {
                    "kind": "task",
                    "runId": "run_1",
                    "task": 1,
                    "func": "f",
                    "mode": "slurm_mcc",
                    "status": "ok",
                    "node": "n1",
                    "pid": 11,
                    "jobId": "55",
                    "start": 110,
                    "end": 112,
                    "wallSec": 2,
                    "processStart": 106.5,
                }
            )
        ],
    )
    write_lines(
        run / "parallel_55_1.joblog",
        [
            "Seq\tHost\tStarttime\tJobRuntime\tSend\tReceive\tExitval\tSignal\tCommand",
            '1\t:\t106.400\t6.100\t0\t0\t0\t0\tOUT="/d/out_1.tif"; PETAKIT5D_TASK=1 run',
            '2\t:\t106.600\t20.000\t0\t0\t137\t9\tOUT="/d/out_2.tif"; PETAKIT5D_TASK=2 run',
            '3\t:\t106.800\t0.010\t0\t0\t0\t0\tOUT="/d/out_3.tif"; PETAKIT5D_TASK=3 run',
        ],
    )


def test_reader_joins_joblog_and_run_settings(tmp_path):
    # a compiled MATLAB run: GNU parallel's joblog gives the time from the end of a task to the
    # end of its process, and lines without a record (the process was killed, or the output
    # existed)
    run = tmp_path / "run_1"
    make_joblog_run(run)
    t, runs = read_task_metrics(str(run))
    t = by_task(t)
    assert len(t) == 3
    assert col(t, "status") == ["ok", "killed", "skipped"]
    assert t[0]["startupSec"] == pytest.approx(3.5, abs=1e-9)
    assert t[0]["lineSec"] == pytest.approx(6.1, abs=1e-9)
    assert t[0]["exitSec"] == pytest.approx(0.5, abs=1e-9)
    assert t[1]["exitCode"] == 137
    assert col(t, "output") == ["/d/out_1.tif", "/d/out_2.tif", "/d/out_3.tif"]
    assert t[0]["input"] == "/d/in_1.tif"
    assert (t[0]["cpusPerTask"], t[0]["memAllocate"]) == (2, 5)
    r = runs[0]
    assert [r["start"], r["end"], r["elapsedSec"], r["nDone"]] == [99, 140, 41, 2]


def test_reader_joins_a_line_up_to_a_second_around_its_record(tmp_path):
    # a record joins a joblog line if it starts at most 1 s before the line and ends at most
    # 1 s after it (clocks of the record and of GNU parallel differ a little)
    run = tmp_path / "run_1"
    run.mkdir()
    rec = {"kind": "task", "runId": "run_1", "task": 1, "jobId": "55"}
    write_lines(
        run / "task_n1_1_1.ndjson",
        [json.dumps({**rec, "start": 99.5, "end": 110.5}), json.dumps({**rec, "start": 98})],
    )
    write_lines(
        run / "parallel_55_1.joblog",
        [
            "1\t:\t100\t10\t0\t0\t0\t0\tPETAKIT5D_TASK=1 a",
            "2\t:\t200\t1\t0\t0\t0\t0\tPETAKIT5D_TASK=1 b",
        ],
    )
    t, _ = read_task_metrics(str(run))
    assert [r["lineSec"] for r in t[:2]] == [10, pytest.approx(math.nan, nan_ok=True)]
    assert len(t) == 3 and t[2]["status"] == "skipped" and t[2]["start"] == 200


def test_reader_reads_lists_as_jsondecode_does(tmp_path):
    # a list of numbers is one numeric array (one element), a list of texts a cell; a task
    # string's function is the text up to its first '(' after leading ones
    run = tmp_path / "run_1"
    run.mkdir()
    write_lines(
        run / "run.ndjson",
        [
            json.dumps(
                {"kind": "run_start", "outputFullpaths": [1, 2, 3], "functionStrs": ["((f(1)"]}
            )
        ],
    )
    _, runs = read_task_metrics(str(run))
    assert (runs[0]["nTasks"], runs[0]["func"]) == (1, "f")


def test_pool_tasks_are_recorded_in_their_processes(tmp_path):
    # (MATLAB: poolTasksAreRecordedOnTheWorkers) tasks other processes run record into their
    # own files of the run, which the reader joins
    folder = tmp_path / "task_metrics" / "run_p"
    folder.mkdir(parents=True)
    code = textwrap.dedent(f"""
        from petakit5d.jobs import metrics
        from petakit5d.jobs.tasks import Call
        call = Call(["/a"] * 4, ["/b"] * 4, [{{"func": "m:f"}}] * 4,
                    {{"task_metrics": {str(folder) + "/"!r}}})
        for i in {{i}}:
            metrics.task_stop(metrics.task_start(call, [i], "pool"))
    """)
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    for part in ("(0, 2)", "(1, 3)"):
        subprocess.run([sys.executable, "-c", code.replace("{i}", part)], env=env, check=True)
    t, _ = read_task_metrics(str(tmp_path))
    assert sorted(col(t, "task")) == [1, 2, 3, 4]
    assert col(t, "mode") == ["pool"] * 4
    assert os.getpid() not in col(t, "pid"), "The tasks did not run in other processes."
    assert len(set(col(t, "pid"))) == 2
    assert len(list(folder.glob("task_*.ndjson"))) == 2


def test_worker_hook_records_a_task_from_its_record(tmp_path):
    # a pool worker has the task record, its index, its worker info and the call folder
    in_fns, out_fns, _ = make_tasks(tmp_path, 3)
    recs = [Task("shutil:copyfile", [i, o]).record() for i, o in zip(in_fns, out_fns, strict=True)]
    call = Call(in_fns, out_fns, recs, {"task_metrics": True, "job_log_dir": str(tmp_path)})
    folder = metrics.run_start(call)
    call_folder = write_call(call, str(tmp_path / "call"))
    worker = types.SimpleNamespace(id=2, threads=3, gpu=None)
    for i, where in ((1, call_folder), (2, folder)):
        tm = metrics.worker_task_start(recs[i], i, where, worker)
        run_record(recs[i])
        metrics.task_stop(tm, None)
    assert metrics.worker_task_start(recs[0], 0, "", worker) is None
    t, _ = read_task_metrics(folder)
    assert [(r["task"], r["mode"], r["slot"], r["threads"]) for r in by_task(t)] == [
        (2, "pool", 2, 3),
        (3, "pool", 2, 3),
    ]
    assert col(t, "func") == ["shutil:copyfile"] * 2
    # without a metrics folder in the call folder: nothing
    off = Call(in_fns, out_fns, recs, {"task_metrics": False})
    assert metrics.worker_task_start(recs[0], 0, write_call(off, str(tmp_path / "c2"))) is None


def test_pool_tasks_are_recorded_on_the_workers(tmp_path):
    pytest.importorskip("petakit5d.jobs.local", reason="the local process pool (J04)")
    in_fns, out_fns, _ = make_tasks(tmp_path, 4)
    tasks = [Task("shutil:copyfile", [i, o]) for i, o in zip(in_fns, out_fns, strict=True)]
    flags = run_tasks(
        in_fns, out_fns, tasks, parse_cluster=False, local_worker_num=2, task_metrics=True
    )
    assert flags == [True] * 4
    t, runs = read_task_metrics(str(tmp_path / "job_logs"))
    assert sorted(col(t, "task")) == [1, 2, 3, 4]
    assert col(t, "mode") == ["pool"] * 4
    assert os.getpid() not in col(t, "pid"), "The tasks did not run on the workers."
    assert len(runs) == 1


# --- more of the writer and the reader ---------------------------------------------------------


def test_record_has_matlabs_fields_in_order(tmp_path):
    folder = str(tmp_path) + "/"
    metrics.task_stop(metrics.task_start(job_call(folder), [0], "local"), 2)
    (fn,) = tmp_path.glob("task_*.ndjson")
    p = metrics.process_info()
    assert fn.name == p["file"]
    assert fn.name.startswith(f"task_{metrics._safe_name(p['host'])}_{os.getpid()}_")
    rec = json.loads(fn.read_text())
    assert list(rec) == [
        "kind",
        "runId",
        "func",
        "mode",
        "status",
        "error",
        *metrics.TASK_STR_FIELDS[5:],
        "version",
        *metrics.TASK_NUM_FIELDS,
    ]
    assert (rec["status"], rec["exitCode"], rec["func"]) == ("failed", 2, "shutil:copyfile")
    assert rec["runId"] == tmp_path.name
    assert rec["version"] == 1


def test_error_texts_are_escaped_and_cut(tmp_path):
    folder = str(tmp_path) + "/"
    err = ValueError('a "quoted"\npath C:\\x\ttab\x01 é' + "x" * 2000)
    metrics.task_stop(metrics.task_start(job_call(folder), [0], "local"), err)
    t, _ = read_task_metrics(folder)
    e = t[0]["error"]
    assert e.startswith('ValueError: a "quoted"\npath C:\\x\ttab  é')
    assert len(e) == 1000


def test_slot_and_slurm_variables_are_recorded(tmp_path, monkeypatch):
    monkeypatch.setenv("PARALLEL_JOBSLOT", "3")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "6")
    monkeypatch.setenv("SLURM_JOB_ID", "4242")
    folder = str(tmp_path) + "/"
    metrics.task_stop(metrics.task_start(job_call(folder), [0], "slurm_python"))
    t, _ = read_task_metrics(folder)
    assert (t[0]["slot"], t[0]["slurmCpus"], t[0]["jobId"]) == (3, 6, "4242")
    assert math.isnan(t[0]["seq"])


def test_run_folder_location(tmp_path, monkeypatch):
    # in the job log folder if it is an existing absolute folder, else <first input>/job_logs
    (tmp_path / "logs").mkdir()
    d = metrics.run_dir(str(tmp_path / "logs"), [[str(tmp_path / "x" / "a.tif")]], ["pkg.m:f(1)"])
    assert d.startswith(str(tmp_path / "logs").replace("\\", "/") + "/task_metrics/")
    assert d.endswith("/") and "_pkg.m_f_" in d
    d = metrics.run_dir("../job_logs", [str(tmp_path / "x" / "a.tif")], "f")
    assert d.startswith((tmp_path / "x").as_posix() + "/job_logs/task_metrics/")
    monkeypatch.chdir(tmp_path)
    d = metrics.run_dir("", ["rel/a.tif"], ["f"])
    assert d.startswith(os.getcwd().replace("\\", "/") + "/rel/job_logs/task_metrics/")


def test_reader_leaves_out_broken_lines_and_other_kinds(tmp_path):
    run = tmp_path / "r"
    run.mkdir()
    good = json.dumps(
        {
            "kind": "task",
            "runId": "r",
            "task": 1,
            "pid": True,
            "start": [5, 6],
            "wallSec": None,
            "node": 7,
            "func": "f",
        }
    )
    (run / "task_a.ndjson").write_text(
        good + '\n{"kind":"task","runId":"r","ta\n{"kind":"x"}\n{"kind":"task","task":NaN}\n\n',
        newline="\n",
    )
    t, runs = read_task_metrics(str(run))
    assert len(t) == 1
    assert math.isnan(t[0]["pid"]) and t[0]["start"] == 5 and math.isnan(t[0]["wallSec"])
    assert t[0]["node"] == ""
    assert runs[0]["nRecords"] == 2
    # MATLAB's NaN (sign bit set) for missing values
    assert math.copysign(1, t[0]["pid"]) == -1


def test_reader_finds_runs_in_a_job_log_folder(tmp_path):
    make_joblog_run(tmp_path / "task_metrics" / "run_1")
    (tmp_path / "task_metrics" / ".hidden").mkdir()
    _, runs = read_task_metrics(str(tmp_path))
    assert col(runs, "runId") == ["run_1"]
    assert runs[0]["folder"] == (tmp_path / "task_metrics" / "run_1").as_posix() + "/"
    assert read_task_metrics(str(tmp_path / "nothing")) == ([], [])
    with pytest.raises(TypeError):
        read_task_metrics(str(tmp_path), slurm_accounting=1)


@pytest.mark.skipif(os.name == "nt", reason="the fake sacct is a shell script")
def test_reader_adds_slurm_accounting(tmp_path, monkeypatch):
    run = tmp_path / "r"
    run.mkdir()
    t0 = time.mktime(time.strptime("2026-11-11T10:00:00", "%Y-%m-%dT%H:%M:%S"))
    write_lines(
        run / "task_x.ndjson",
        [
            json.dumps(
                {
                    "kind": "task",
                    "runId": "r",
                    "task": 1,
                    "jobId": "55",
                    "start": t0 + 70,
                    "end": t0 + 80,
                }
            ),
            json.dumps(
                {
                    "kind": "task",
                    "runId": "r",
                    "task": 2,
                    "jobId": "55",
                    "start": t0 + 81,
                    "end": t0 + 90,
                }
            ),
            json.dumps(
                {"kind": "task", "runId": "r", "task": 3, "jobId": "", "start": t0, "end": t0}
            ),
        ],
    )
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "sacct").write_text(
        '#!/bin/sh\necho "$*" >> "$(dirname "$0")/sacct.log"\n'
        "echo '55|2026-11-11T09:59:00|2026-11-11T10:01:00|2026-11-11T10:01:35|CANCELLED by 1'\n",
        newline="\n",
    )
    (bin_dir / "sacct").chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    t, _ = read_task_metrics(str(run), slurm_accounting=True)
    assert "-j 55" in (bin_dir / "sacct.log").read_text()
    a, b, c = by_task(t)
    assert (a["queueSec"], a["jobState"], a["jobHeadSec"], a["jobTailSec"]) == (
        120,
        "CANCELLED",
        10,
        5,
    )
    assert b["jobHeadSec"] == 10 and math.isnan(c["queueSec"]) and c["jobState"] == ""
    (bin_dir / "sacct").write_text("#!/bin/sh\necho no >&2\nexit 1\n", newline="\n")
    with pytest.warns(TaskMetricsWarning) as rec:
        read_task_metrics(str(run), slurm_accounting=True)
    assert rec[0].message.identifier == "read_task_metrics:sacct"


def test_env_numbers_as_matlabs_sscanf(monkeypatch):
    for text, want in (("3", 3.0), (" 2.5x", 2.5), ("1e3", 1000.0), ("x", None), ("", None)):
        monkeypatch.setenv("PK5D_TEST_NUM", text)
        x = metrics._env_num("PK5D_TEST_NUM")
        assert math.isnan(x) if want is None else x == want, text


class warnings_as_list:
    def __enter__(self):
        import warnings

        self._cm = warnings.catch_warnings(record=True)
        ws = self._cm.__enter__()
        warnings.simplefilter("always")
        return ws

    def __exit__(self, *exc):
        return self._cm.__exit__(*exc)


# --- oracle: MATLAB's read_task_metrics reads folders Python wrote ---------------------------


class _Clocks:
    """Deterministic clocks, /proc and environment for the metrics module."""

    def __init__(self, env):
        self.t, self.p, self.c, self.k, self.n = 1_800_000_000.0, 100.0, 5.0, 0, 0
        self.env = env

    def time(self):
        self.t += 0.25
        return self.t

    def perf_counter(self):
        self.p += 0.125
        return self.p

    def process_time(self):
        self.c += 0.0625
        return self.c

    def strftime(self, fmt, *args):
        return "20261111_120000"

    def snapshot(self):
        self.k += 1
        k = float(self.k)
        return [
            262144 + 1024 * k,
            300000 + 2048 * k,
            4 + k % 3,
            8,
            1000 * k,
            k % 2,
            10 * k,
            2 * k,
            k,
            0,
            12345,
            4096 * k,
            8192 * k,
            0,
            65536 * k,
        ]

    def times(self):
        return (0.5, 0.25, 0.0, 0.0, 1.0)

    def token(self):
        self.n += 1
        return f"{self.n:08x}"

    def process(self):
        return {
            "linux": True,
            "pid": 4242,
            "hz": 100.0,
            "host": "node12",
            "processStart": 1_799_999_990.5,
            "memTotalGB": 512.0,
            "threads": 4.0,
            "file": "task_node12_4242_12345.ndjson",
        }


def _patched(clocks):
    fake_time = types.SimpleNamespace(
        time=clocks.time,
        perf_counter=clocks.perf_counter,
        process_time=clocks.process_time,
        strftime=clocks.strftime,
    )
    return {
        "time": fake_time,
        "_snapshot": clocks.snapshot,
        "_cpu_times": clocks.times,
        "_token": clocks.token,
        "_new_process_info": clocks.process,
        "_env": lambda name: clocks.env.get(name, ""),
        "_process": None,
        "_current": "",
        "_warned": False,
        "_sessions": {},
        "_run_t0": {},
    }


def make_python_runs(root: Path) -> None:
    """Run folders written by run_tasks with metrics on, with clocks, /proc, host and Slurm's
    variables made deterministic: three copy tasks, a failing task, a failing bash task, a
    task that runs two framework calls of its own, and a single job's submission with its
    task in a job."""
    env: dict[str, str] = {}
    clocks = _Clocks(env)
    saved = {k: getattr(metrics, k) for k in _patched(clocks)}
    cwd = os.getcwd()
    try:
        for k, v in _patched(clocks).items():
            setattr(metrics, k, v)
        os.chdir(root)
        (root / "data" / "out").mkdir(parents=True)
        ins = [f"data/in_{i}.txt" for i in range(1, 5)]
        outs = [f"data/out/out_{i}.txt" for i in range(1, 5)]
        for i, fn in enumerate(ins, 1):
            Path(fn).write_bytes(f"task {i}\n".encode())
        # (unit_wait_time 0: the failed task leaves its .tmp marker, which the later runs on
        # the same output would wait for)
        opts = {"task_metrics": True, "uuid": "0123abcd", "unit_wait_time": 0, **LOCAL}
        run_tasks(
            ins[:3],
            outs[:3],
            [Task("shutil:copyfile", [i, o]) for i, o in zip(ins, outs[:3], strict=False)],
            **opts,
        )
        try:
            run_tasks(ins[3:], outs[3:], [Task("os:remove", ["data/missing.txt"])], **opts)
        except FileNotFoundError:
            pass
        run_tasks(
            ins[3:],
            outs[3:],
            ["exit 3"],
            language="bash",
            max_trial_num=1,
            **opts,
        )
        inner = {"func": "shutil:copyfile", "args": [ins[3], outs[3]]}
        nested = Task(
            "petakit5d.jobs.framework:run_tasks",
            [[ins[3]], [outs[3]], [inner]],
            {"task_metrics": True, "uuid": "fedcba98", "unit_wait_time": 0, **LOCAL},
        )
        run_tasks(ins[3:], outs[3:], [nested], **opts)
        # a single job: the submitting process's record, then the job's task
        folder = metrics.session_dir(str(root / "data" / "job_logs"))
        metrics.record(
            folder,
            "submit",
            {
                "task": 7,
                "jobId": 1000,
                "functionStr": "shutil:copyfile(1)",
                "cpusPerTask": 2,
                "memAllocate": 30,
                "memPerCPU": 20.9,
                "mccMode": False,
                "SlurmParam": "-p compute",
            },
        )
        env.update(
            SLURM_JOB_ID="1001",
            SLURM_ARRAY_JOB_ID="1000",
            SLURM_ARRAY_TASK_ID="7",
            PARALLEL_JOBSLOT="2",
            SLURM_CPUS_PER_TASK="2",
        )
        call = Call(
            ["/a"] * 7, ["/b"] * 7, [{"func": "shutil:copyfile"}] * 7, {"task_metrics": folder}
        )
        metrics.task_stop(metrics.task_start(call, [6], "slurm_python"))
    finally:
        os.chdir(cwd)
        for k, v in saved.items():
            setattr(metrics, k, v)


def make_reader_runs(root: Path) -> None:
    """Run folders by hand, with what the reader meets: MATLAB's joblog run, a run with
    numbers as texts and bools, arrays and nulls, broken lines, texts with escapes, a single
    job session, a nested parent, records of other kinds."""
    tm = root / "task_metrics"
    tm.mkdir()
    make_joblog_run(tm / "run_1")
    run = tm / "run_2"
    run.mkdir()
    write_lines(
        run / "run.ndjson",
        [
            json.dumps(
                {
                    "kind": "run_start",
                    "runId": "run_2",
                    "time": 200.5,
                    "jobId": "",
                    "parseCluster": False,
                    "backend": "python",
                    "cpusPerTask": [3, 4],
                    "memAllocate": None,
                    "GPUJob": True,
                    "taskBatchNum": 2,
                    "inputFullpaths": [["/d/a_1", "/d/b_1"], ["/d/a_2"]],
                    "outputFullpaths": "/d/only",
                    "functionStrs": ["(g(1)", "g(2)"],
                    "parent": "run_0/3",
                }
            ),
            "not json",
            json.dumps(
                {"kind": "submit", "time": 150, "jobId": 9, "cpusPerTask": 1, "memAllocate": 2}
            ),
            json.dumps({"kind": "submit", "time": 210, "jobId": 2.5, "functionStr": "h(1)"}),
        ],
    )
    write_lines(
        run / "task_b_2_3.ndjson",
        [
            json.dumps(
                {
                    "kind": "task",
                    "runId": "run_2",
                    "task": 1,
                    "mode": "mcc",
                    "jobId": "",
                    "func": "g",
                    "error": 'x "y"\n\tz é',
                    "start": 201,
                    "end": 202,
                    "pid": 3,
                    "node": "b",
                    "processStart": 199.25,
                    "arrayJobId": "9",
                }
            ),
            json.dumps(
                {
                    "kind": "task",
                    "runId": "run_2",
                    "task": 2,
                    "mode": "slurm_python",
                    "pid": 3,
                    "node": "b",
                    "start": 203,
                    "processStart": 199.25,
                    "gpu": ["x"],
                    "exitCode": [[4, 5]],
                }
            ),
            json.dumps(
                {
                    "kind": "task",
                    "runId": "run_2",
                    "task": 2,
                    "mode": "slurm_python",
                    "pid": 3,
                    "node": "b",
                    "start": 202.5,
                    "processStart": 199.25,
                }
            ),
            json.dumps(
                {
                    "kind": "task",
                    "runId": "run_9",
                    "task": 1,
                    "mode": "slurm_matlab",
                    "pid": 5,
                    "node": "c",
                    "start": 300,
                    "end": 301,
                    "processStart": 290.4,
                }
            ),
            '{"kind":"task","runId":"run_2","task":',
            json.dumps({"kind": "other", "task": 3}),
        ],
    )


def read_as_matlab(path: str, sub: str = ""):
    """The tables as MATLAB's table2struct(t, 'ToScalar', true) has them (columns), with
    MATLAB's start-up modes (MATLAB's reader doesn't know slurm_python, which only we write)."""
    t, r = metrics._read(path + sub, False, metrics.MATLAB_OWN_PROCESS_MODES)
    return _columns(t), _columns(r)


def _columns(rows):
    # (the table column 'end' is the struct field 'xEnd' in MATLAB: 'end' is a keyword)
    name = {"end": "xEnd"}
    return {name.get(k, k): [row[k] for row in rows] for k in rows[0]} if rows else {}


MATLAB_READ = (
    "[t, r] = read_task_metrics([in{1}, in{2}]); "
    "out = {table2struct(t, 'ToScalar', true), table2struct(r, 'ToScalar', true)};"
)

CASES = [
    Case(
        "read_python_run_folders",
        code=MATLAB_READ,
        args=[InPath("pyrun"), "data/job_logs"],
        nargout=2,
        python=read_as_matlab,
        inputs={"pyrun": Tree(make_python_runs, "pyrun")},
        outputs=["tasks", "runs"],
    ),
    Case(
        "read_hand_written_run_folders",
        code=MATLAB_READ,
        args=[InPath("runs"), ""],
        nargout=2,
        python=read_as_matlab,
        inputs={"runs": Tree(make_reader_runs, "runs")},
        outputs=["tasks", "runs"],
    ),
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_python_run_folders_are_deterministic(tmp_path):
    # the oracle input is made twice the same; and what it holds
    for name in ("a", "b"):
        (tmp_path / name).mkdir()
        make_python_runs(tmp_path / name)
    files = sorted(
        p.relative_to(tmp_path / "a").as_posix() for p in (tmp_path / "a").rglob("*") if p.is_file()
    )
    for f in files:
        assert (tmp_path / "a" / f).read_bytes() == (tmp_path / "b" / f).read_bytes(), f
    t, runs = read_task_metrics(str(tmp_path / "a" / "data" / "job_logs"))
    assert sorted(col(t, "status")) == ["error", "failed"] + ["ok"] * 6
    assert len(runs) == 6
    assert sorted(col(runs, "kind")) == ["framework"] * 5 + ["single"]
    assert [r["mode"] for r in t if r["arrayJobId"] == "1000"] == ["slurm_python"]
