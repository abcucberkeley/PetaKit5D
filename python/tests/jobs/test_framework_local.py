"""run_tasks in local mode (parse_cluster=False), task records, call folders and the job entry
point. Ported from tests/framework/GenericComputingFrameworkTest.m (local-mode tests) and
GenericComputingFrameworkDemoTest.m; the test names keep the MATLAB names in snake_case.

Tasks call framework_task, which logs every call (so runs can be counted), as
pktest.frameworkTask does.
"""

from __future__ import annotations

import functools
import json
import math
import os
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

import petakit5d
from petakit5d.jobs import framework, tasks
from petakit5d.jobs.framework import RUN_TASKS_PARAMS, run_tasks
from petakit5d.jobs.tasks import Call, CallFolder, Task, read_call, task_record, write_call

SRC_ROOT = str(Path(petakit5d.__file__).resolve().parents[1])


class TaskFailed(Exception):
    """pktest:frameworkTask:failed"""


def framework_task(input_fn: str, output_fn: str, log_fn: str, mode: str = "copy") -> None:
    """Appends output_fn to log_fn, then by mode: 'copy' the input to the output, 'error'
    raise, 'noop' nothing, 'secondRun' as noop in the first run and copy from the second,
    'sleep:<s>' copy after s seconds (by a temporary name and a rename)."""
    _append_line(log_fn, output_fn)
    if mode == "copy":
        Path(output_fn).write_bytes(Path(input_fn).read_bytes())
    elif mode == "error":
        raise TaskFailed(f"Task failed on purpose for {input_fn}")
    elif mode == "noop":
        pass
    elif mode == "secondRun":
        if Path(log_fn).read_text().splitlines().count(output_fn) >= 2:
            Path(output_fn).write_bytes(Path(input_fn).read_bytes())
    elif mode.startswith("sleep:"):
        time.sleep(float(mode[6:]))
        tmp = f"{output_fn}_{os.getpid()}"
        Path(tmp).write_bytes(Path(input_fn).read_bytes())
        os.replace(tmp, output_fn)
    else:
        raise ValueError(mode)


def _append_line(log_fn: str, line: str) -> None:
    # one write with O_APPEND: lines of tasks in other processes don't mix
    fd = os.open(log_fn, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        os.write(fd, (line + "\n").encode())
    finally:
        os.close(fd)


FUNC = f"{__name__}:framework_task"


def make_tasks(tmp_path: Path, n: int):
    """n input text files, their output paths in out/, the call log."""
    (tmp_path / "out").mkdir()
    in_fns = [str(tmp_path / f"in_{i}.txt") for i in range(1, n + 1)]
    out_fns = [str(tmp_path / "out" / f"out_{i}.txt") for i in range(1, n + 1)]
    for i, fn in enumerate(in_fns, 1):
        Path(fn).write_text(f"task {i}")
    return in_fns, out_fns, str(tmp_path / "calls.log")


def task_list(in_fns, out_fns, log_fn, mode):
    return [Task(FUNC, [i, o, log_fn, mode]) for i, o in zip(in_fns, out_fns, strict=True)]


def read_log(log_fn: str) -> list[str]:
    if not os.path.exists(log_fn):
        return []
    return [x for x in Path(log_fn).read_text().splitlines() if x]


def tmp_files(folder: Path) -> list[str]:
    return sorted(p.name for p in folder.glob("*.tmp"))


LOCAL = {"parse_cluster": False, "master_compute": True}


# --- GenericComputingFrameworkTest (local mode) ---------------------------------------------


def test_runs_every_task_once_and_writes_outputs(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert flags == [True] * 3
    for i, o in zip(in_fns, out_fns, strict=True):
        assert Path(o).read_text() == Path(i).read_text()
    assert sorted(read_log(log_fn)) == sorted(out_fns), "Each task should run exactly once."
    assert tmp_files(tmp_path / "out") == [], "Task flag files should be removed after success."


def test_skips_tasks_whose_outputs_exist(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    for i, o in zip(in_fns, out_fns, strict=True):
        Path(o).write_bytes(Path(i).read_bytes())
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "error"), **LOCAL)
    assert flags == [True] * 3
    assert read_log(log_fn) == [], "No task should run when all outputs exist."


def test_runs_only_tasks_with_missing_outputs(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    for k in range(2):
        Path(out_fns[k]).write_bytes(Path(in_fns[k]).read_bytes())
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert flags == [True] * 3
    assert read_log(log_fn) == out_fns[2:]


def test_failing_task_is_retried_once_then_rethrown(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    with pytest.raises(TaskFailed):
        run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "error"), **LOCAL)
    # one run in try, one more in except, then the exception propagates
    assert len(read_log(log_fn)) == 2


def test_task_without_output_is_retried_up_to_max_trial_num(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    # unit_wait_time = 0, otherwise each retry waits for minutes for the task's own marker
    flags = run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "noop"),
        **LOCAL,
        max_trial_num=2,
        unit_wait_time=0,
    )
    assert flags == [False, False]
    assert sorted(read_log(log_fn)) == sorted(out_fns * 2), (
        "Each task should run maxTrialNum times."
    )
    # the flag files of failed tasks stay behind
    assert tmp_files(tmp_path / "out") == ["out_1.tmp", "out_2.tmp"]


def test_task_done_in_its_last_try_counts_as_done(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    flags = run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "secondRun"),
        **LOCAL,
        max_trial_num=2,
        unit_wait_time=0,
    )
    assert flags == [True, True], "A task done in its last try counts as unfinished."
    assert all(os.path.isfile(o) for o in out_fns)
    assert sorted(read_log(log_fn)) == sorted(out_fns * 2)


def test_stale_task_flag_from_crashed_run_is_ignored(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    # flag file left by a worker that died an hour ago
    tmp_fn = tmp_path / "out" / "out_2.tmp"
    tmp_fn.touch()
    old = time.time() - 3600
    os.utime(tmp_fn, (old, old))
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert flags == [True, True]
    assert os.path.isfile(out_fns[1])
    assert not tmp_fn.exists(), "The stale flag file should be removed."


def test_missing_input_is_skipped_and_reported_done(tmp_path):
    # characterizes MATLAB: tasks whose input is missing are not run, and their done flag
    # is true even though no output exists
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    os.remove(in_fns[1])
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert flags == [True, True]
    assert read_log(log_fn) == out_fns[:1]
    assert not os.path.exists(out_fns[1])


def test_existing_final_output_skips_all_tasks(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    final_fn = tmp_path / "final.txt"
    final_fn.touch()
    flags = run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "error"),
        **LOCAL,
        final_out_fullpath=str(final_fn),
    )
    assert flags == [True, True]
    assert read_log(log_fn) == []


def test_final_output_made_by_a_task_ends_the_run(tmp_path):
    # the first task writes the final output; the second never writes its own and would be
    # retried, but the run ends after the first pass once the final output exists
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    out_fns[0] = str(tmp_path / "final.txt")
    task = task_list(in_fns[:1], out_fns[:1], log_fn, "copy")
    task += task_list(in_fns[1:], out_fns[1:], log_fn, "noop")
    flags = run_tasks(
        in_fns, out_fns, task, **LOCAL, unit_wait_time=0, final_out_fullpath=out_fns[0]
    )
    assert flags == [True, True]
    assert read_log(log_fn) == out_fns, "Each task should run once."


def test_task_batches_run_each_task_once(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 4)
    flags = run_tasks(
        in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL, task_batch_num=2
    )
    assert flags == [True] * 4
    assert sorted(read_log(log_fn)) == sorted(out_fns)


def test_config_file_overrides_parameters(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    config_fn = tmp_path / "config.json"
    config_fn.write_text(json.dumps({"maxTrialNum": 1}))
    # max_trial_num = 3 here, the config file sets it to 1
    run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "noop"),
        **LOCAL,
        max_trial_num=3,
        unit_wait_time=0,
        config_file=str(config_fn),
    )
    assert sorted(read_log(log_fn)) == sorted(out_fns), (
        "With maxTrialNum = 1 from the config file, each task should run once."
    )


def test_config_file_keys_ignore_case(tmp_path):
    # (local part of configFileKeysIgnoreCase: MemPerCPU sets memPerCPU; the CPUs of the
    # jobs are J03's)
    config_fn = tmp_path / "config.json"
    config_fn.write_text(
        json.dumps(
            {
                "MemPerCPU": 5,
                "MAXTRIALNUM": 1,
                "MCRParam": "/opt/mcr",
                "mccMode": True,
                "someUnknownKey": 1,
            }
        )
    )
    settings = {p.python: p.default for p in RUN_TASKS_PARAMS}
    out = framework.apply_config(settings, str(config_fn))
    assert out["mem_per_cpu"] == 5
    assert out["max_trial_num"] == 1
    assert set(out) == set(settings), "MATLAB-only and unknown keys are ignored."
    # MATLAB's jsonencode writes empty values as []
    config_fn.write_text(json.dumps({"memAllocate": [], "jobTimeLimit": [], "jobLogDir": ""}))
    out = framework.apply_config(settings, str(config_fn))
    assert out["mem_allocate"] is None and out["job_time_limit"] is None
    assert out["job_log_dir"] == ""


def test_missing_config_file_errors(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    with pytest.raises(FileNotFoundError, match="does not exist"):
        run_tasks(
            in_fns,
            out_fns,
            task_list(in_fns, out_fns, log_fn, "copy"),
            parse_cluster=False,
            config_file=str(tmp_path / "missing.json"),
        )
    assert read_log(log_fn) == []


def test_master_para_factor_accepts_fractions(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    flags = run_tasks(
        in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL, master_para_factor=0.5
    )
    assert flags == [True, True]


# --- more local-mode behavior of MATLAB's loop ----------------------------------------------


def test_failing_task_in_a_batch_reruns_the_whole_batch(tmp_path):
    # MATLAB evals the batch's joined function string again from its first task
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    task = task_list(in_fns[:1], out_fns[:1], log_fn, "copy")
    task += task_list(in_fns[1:], out_fns[1:], log_fn, "error")
    with pytest.raises(TaskFailed):
        run_tasks(in_fns, out_fns, task, **LOCAL, task_batch_num=2)
    assert read_log(log_fn) == [out_fns[0], out_fns[1], out_fns[0], out_fns[1]]
    # the batch's marker (that of its last task) stays
    assert tmp_files(tmp_path / "out") == ["out_2.tmp"]


def test_batch_runs_all_its_tasks_even_if_some_are_done(tmp_path):
    # MATLAB runs the batch's whole function string; the task functions skip done work
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    Path(out_fns[0]).write_text("done")
    flags = run_tasks(
        in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL, task_batch_num=2
    )
    assert flags == [True, True]
    assert read_log(log_fn) == out_fns


def test_fresh_marker_of_another_caller_leaves_the_task_to_it(tmp_path, monkeypatch):
    # a marker newer than unitWaitTime minutes: the task is left; here the "other caller"
    # writes the output during the wait
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    (tmp_path / "out" / "out_1.tmp").touch()
    sleeps = []

    def fake_sleep(t):
        sleeps.append(t)
        Path(out_fns[0]).write_text("by the other caller")
        os.remove(tmp_path / "out" / "out_1.tmp")

    monkeypatch.setattr(framework.time, "sleep", fake_sleep)
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert flags == [True, True]
    assert read_log(log_fn) == out_fns[1:]
    assert sleeps == [framework.LOCAL_IDLE_WAIT]
    assert tmp_files(tmp_path / "out") == []


def test_a_round_that_ran_nothing_waits(tmp_path, monkeypatch):
    # one task left to another caller: the loop sleeps between its looks (MATLAB spins)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    (tmp_path / "out" / "out_1.tmp").touch()
    sleeps, looks = [], []
    real_bfe = framework.batch_file_exist

    def counting_bfe(*args):
        looks.append(len(sleeps))
        assert len(looks) < 50, "the loop looks without waiting"
        return real_bfe(*args)

    def fake_sleep(t):
        sleeps.append(t)
        if len(sleeps) == 3:
            Path(out_fns[0]).write_text("by the other caller")

    monkeypatch.setattr(framework, "batch_file_exist", counting_bfe)
    monkeypatch.setattr(framework.time, "sleep", fake_sleep)
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert flags == [True] and read_log(log_fn) == []
    assert sleeps == [framework.LOCAL_IDLE_WAIT] * 3


def test_markers_of_folder_outputs_with_trailing_separators(tmp_path):
    # the output's trailing separators are removed first: the marker is <folder>/x.tmp
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    out = [str(tmp_path / "out" / "x.zarr") + "/"]
    flags = run_tasks(in_fns, out, task_list(in_fns, out, log_fn, "noop"), **LOCAL, max_trial_num=1)
    assert flags == [False]
    assert tmp_files(tmp_path / "out") == ["x.tmp"]


def test_markers_go_to_tmp_dir_and_take_the_name_without_extension(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    out_fns = [str(tmp_path / "out" / "a.b.tif"), str(tmp_path / "out" / ".hidden")]
    marks = tmp_path / "marks"
    marks.mkdir()
    flags = run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "noop"),
        **LOCAL,
        max_trial_num=1,
        tmp_dir=str(marks),
    )
    assert flags == [False, False]
    assert tmp_files(marks) == [".hidden.tmp", "a.b.tmp"]
    assert tmp_files(tmp_path / "out") == []


def test_no_marker_without_the_output_folder(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    out_fns = [str(tmp_path / "later" / "out_1.txt")]
    flags = run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "noop"),
        **LOCAL,
        max_trial_num=2,
        unit_wait_time=3,
    )
    # no marker, so the second try runs at once
    assert flags == [False]
    assert read_log(log_fn) == out_fns * 2


def test_list_inputs_are_all_checked(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    inputs = [[in_fns[0], in_fns[1]], [in_fns[2], str(tmp_path / "missing.txt")], in_fns[1]]
    flags = run_tasks(inputs, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert flags == [True, True, True]
    assert read_log(log_fn) == [out_fns[0], out_fns[2]]


def test_single_task_with_a_list_input(tmp_path):
    # (MATLAB fails here: batch_file_exist gets a cell in a cell)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    flags = run_tasks(
        [in_fns], out_fns[:1], task_list(in_fns[:1], out_fns[:1], log_fn, "copy"), **LOCAL
    )
    assert flags == [True]


def test_folder_outputs_with_trailing_separators(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    out = tmp_path / "out" / "result.zarr"
    out.mkdir()
    flags = run_tasks(
        in_fns, [str(out) + "//"], task_list(in_fns, [str(out)], log_fn, "error"), **LOCAL
    )
    assert flags == [True]


def test_bash_tasks(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    if os.name == "nt":
        pytest.skip("POSIX shell commands")
    cmds = [f"cp '{i}' '{o}'; echo {o} >> '{log_fn}'" for i, o in zip(in_fns, out_fns, strict=True)]
    cmds[1] = f"echo {out_fns[1]} >> '{log_fn}'; exit 3"
    flags = run_tasks(
        in_fns,
        out_fns,
        cmds,
        **LOCAL,
        language="bash",
        max_trial_num=2,
        unit_wait_time=0,
        task_batch_num=1,
    )
    # a failing command is not retried at once (its status is ignored), only in later rounds
    assert flags == [True, False, True]
    assert read_log(log_fn) == [out_fns[0], out_fns[1], out_fns[2], out_fns[1]]


def test_bash_batches_run_as_one_shell(tmp_path):
    if os.name == "nt":
        pytest.skip("POSIX shell commands")
    in_fns, out_fns, _ = make_tasks(tmp_path, 2)
    cmds = [f"X=hello; touch '{out_fns[0]}'", f"echo $X > '{out_fns[1]}'"]
    flags = run_tasks(in_fns, out_fns, cmds, **LOCAL, language="bash", task_batch_num=2)
    assert flags == [True, True]
    assert Path(out_fns[1]).read_text() == "hello\n"


def test_language_matlab_means_python_tasks(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    flags = run_tasks(
        in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL, language="Matlab"
    )
    assert flags == [True]


def test_single_task_with_master_compute_runs_here(tmp_path, monkeypatch):
    # one task with masterCompute never goes to the cluster (no Slurm check)
    def no_slurm_check():
        raise AssertionError("Slurm was checked")

    monkeypatch.setattr(framework, "_slurm_available", no_slurm_check)
    # (outside a Slurm job: in one, the task must also fit the job's memory)
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    flags = run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"))
    assert flags == [True]


def test_without_a_cluster_tasks_run_here(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(framework, "_slurm_available", lambda: False)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    flags = run_tasks(
        in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), master_compute=False
    )
    assert flags == [True, True]
    assert "Set parseCluster as false" in capsys.readouterr().err


def test_matlab_only_parameters_are_accepted(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    flags = run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "copy"),
        **LOCAL,
        mcc_mode=True,
        matlab_launch_str="matlab -batch",
        mcr_param="/opt/mcr",
        mcr_cache_root="/tmp",
        mcc_master_str="/opt/run_mccMaster.sh",
    )
    assert flags == [True]


def test_parameter_types_are_checked(tmp_path):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 1)
    t = task_list(in_fns, out_fns, log_fn, "copy")
    with pytest.raises(TypeError, match="parseCluster"):
        run_tasks(in_fns, out_fns, t, parse_cluster=1)
    with pytest.raises(TypeError, match="maxTrialNum"):
        run_tasks(in_fns, out_fns, t, **LOCAL, max_trial_num="3")
    with pytest.raises(ValueError, match="language"):
        run_tasks(in_fns, out_fns, t, **LOCAL, language="perl")
    with pytest.raises(ValueError, match="each task"):
        run_tasks(in_fns, out_fns * 2, t, **LOCAL)
    assert read_log(log_fn) == []


def test_no_tasks():
    assert run_tasks([], [], [], **LOCAL) == []


def test_parameters_match_matlab():
    # names and defaults of generic_computing_frameworks_wrapper (site defaults empty)
    p = {x.matlab: x.default for x in RUN_TASKS_PARAMS}
    assert p["maxCPUNum"] == 24 and p["memPerCPU"] == 20.9 and p["maxTrialNum"] == 3
    assert p["unitWaitTime"] == 3 and p["maxJobNum"] == math.inf and p["jobTimeLimit"] == 24
    assert p["jobLogDir"] == "../job_logs" and p["queryInterval"] == 3
    assert p["taskRunner"] is True and p["masterCompute"] is True and p["parseCluster"] is True
    assert p["SlurmParam"] == "" and p["memAllocate"] is None
    import inspect

    sig = inspect.signature(run_tasks).parameters
    for x in RUN_TASKS_PARAMS:
        assert sig[x.python].default == x.default, x.python


# --- task records ---------------------------------------------------------------------------


def module_level_function(a, b=0):
    return a + b


def test_task_records_are_canonical_json():
    rec = task_record(
        Task(FUNC, (1, 2.5, (3, 4), Path("/a/b")), {"x": None, "y": [True]}), "python"
    )
    assert rec == {
        "func": FUNC,
        "args": [1, 2.5, [3, 4], "/a/b"],
        "kwargs": {"x": None, "y": [True]},
    }
    values = [0.1, 1 / 3, -0.0, 5e-324, 1.7976931348623157e308, math.inf, -math.inf]
    back = task_record(Task(FUNC, values), "python")["args"]
    assert [v.hex() for v in back] == [v.hex() for v in values]
    assert math.isnan(task_record(Task(FUNC, [math.nan]), "python")["args"][0])
    with pytest.raises(TypeError, match=r"args\[0\]"):
        Task(FUNC, [object()])
    with pytest.raises(TypeError, match="keys"):
        Task(FUNC, [], {1: 2})


def test_task_records_take_numpy_values():
    np = pytest.importorskip("numpy")
    t = Task(FUNC, [np.array([1, 2], dtype=np.uint16), np.float32(0.1), np.bool_(True)])
    assert t.args == [[1, 2], float(np.float32(0.1)), True]
    assert type(t.args[2]) is bool


def test_task_of_a_function_and_a_partial():
    t = Task.of(module_level_function, 1, b=2)
    assert t.func == f"{__name__}:module_level_function"
    assert tasks.run_record(t.record()) == 3
    p = Task.of(functools.partial(module_level_function, 5), b=1)
    assert p.args == [5] and p.kwargs == {"b": 1}
    with pytest.raises(ValueError, match="module level"):
        Task.of(lambda: None)


def test_task_names_are_checked():
    for bad in ("module", "a.b:", ":f", "a b:f", "a.b:f-g", "a:b:c"):
        with pytest.raises(ValueError):
            Task(bad)
    with pytest.raises(TypeError, match="language='bash'"):
        task_record("echo hi", "python")
    with pytest.raises(TypeError, match="language='bash'"):
        task_record({"bash": "echo hi"}, "python")
    assert task_record("echo hi", "bash") == {"bash": "echo hi"}
    assert task_record({"func": FUNC, "args": [1]}, "python")["kwargs"] == {}


def test_call_folder_round_trip(tmp_path):
    recs = [task_record(Task(FUNC, [i, [i] * i, {"k": i / 7}]), "python") for i in range(7)]
    inputs = [f"/data/in_{i}.tif" for i in range(6)] + [["/data/a.tif", "/data/b.tif"]]
    call = Call(
        inputs,
        [f"/data/out/out_{i}.tif" for i in range(7)],
        recs,
        {"task_batch_num": 3, "language": "python", "uuid": "abc", "max_job_num": math.inf},
    )
    folder = write_call(call, tmp_path / "call")
    cf = CallFolder(folder)
    assert cf.n == 7 and list(cf.batch(3)) == [6] and list(cf.batch(2)) == [3, 4, 5]
    assert [x["task"] for x in cf.lines([5, 0])] == [recs[5], recs[0]]
    back = read_call(folder)
    assert back.records == recs and back.inputs == inputs and back.outputs == call.outputs
    assert back.settings["max_job_num"] == math.inf and back.batch_num == 3
    assert sorted(os.listdir(folder)) == ["call.json", "tasks.ndjson"]


# --- the job entry point and user modules ---------------------------------------------------

USER_MODULE = """
import os, sys

def make(inp, out, factor, log):
    with open(log, "a") as fh:
        fh.write(f"{out} {os.getpid()}\\n")
    data = open(inp).read()
    tmp = out + ".part"
    with open(tmp, "w") as fh:
        fh.write(data * factor)
    os.replace(tmp, out)

def fail(log):
    with open(log, "a") as fh:
        fh.write("fail\\n")
    raise RuntimeError("on purpose")

def runner(folder, batch):
    print(f"runner {batch}")
    return 7
"""


@pytest.fixture
def user_module(tmp_path, monkeypatch):
    """A user package in a folder that is not on sys.path."""
    root = tmp_path / "userlib"
    (root / "usertasks").mkdir(parents=True)
    (root / "usertasks" / "__init__.py").write_text("")
    (root / "usertasks" / "jobs.py").write_text(USER_MODULE)
    monkeypatch.setattr(sys, "path", list(sys.path))
    yield root
    for name in [m for m in sys.modules if m.startswith("usertasks")]:
        del sys.modules[name]


def run_exec(folder, batch):
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    return subprocess.run(
        [sys.executable, "-m", "petakit5d.jobs.exec", str(folder), str(batch)],
        capture_output=True,
        text=True,
        env=env,
        timeout=120,
    )


def test_user_function_in_a_user_module_runs_locally(tmp_path, user_module):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    t = [
        Task("usertasks.jobs:make", [i, o, 2, log_fn], path=[str(user_module)])
        for i, o in zip(in_fns, out_fns, strict=True)
    ]
    flags = run_tasks(in_fns, out_fns, t, **LOCAL)
    assert flags == [True] * 3
    assert Path(out_fns[2]).read_text() == "task 3task 3"


def test_task_of_a_user_function_records_its_import_root(tmp_path, user_module):
    sys.path.insert(0, str(user_module))
    from usertasks import jobs as user_jobs

    t = Task.of(user_jobs.make, "a", "b", 1, "c")
    assert t.func == "usertasks.jobs:make" and t.path == [str(user_module)]


def test_user_function_runs_through_the_job_entry_point(tmp_path, user_module):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 5)
    recs = [
        task_record(
            Task("usertasks.jobs:make", [i, o, 3, log_fn], path=[str(user_module)]), "python"
        )
        for i, o in zip(in_fns, out_fns, strict=True)
    ]
    call = Call(in_fns, out_fns, recs, {"task_batch_num": 2, "language": "python"})
    folder = write_call(call, tmp_path / "call")
    for b in (1, 2, 3):
        r = run_exec(folder, b)
        assert r.returncode == 0, r.stderr
        assert r.stdout.strip().splitlines()[-1].startswith("Elapsed time is ")
    assert [Path(o).read_text() for o in out_fns] == [f"task {i}" * 3 for i in range(1, 6)]
    # each batch in one process
    pids = [line.split()[1] for line in read_log(log_fn)]
    assert pids[0] == pids[1] and pids[2] == pids[3] and len(set(pids)) == 3


def test_job_entry_point_stops_at_the_first_exception(tmp_path, user_module):
    log_fn = str(tmp_path / "log.txt")
    recs = [
        task_record(Task("usertasks.jobs:fail", [log_fn], path=[str(user_module)]), "python"),
        task_record(Task("usertasks.jobs:fail", [log_fn], path=[str(user_module)]), "python"),
    ]
    folder = write_call(Call(["a", "b"], ["c", "d"], recs, {"task_batch_num": 2}), tmp_path / "c")
    r = run_exec(folder, 1)
    assert r.returncode == 1
    assert "RuntimeError: on purpose" in r.stderr
    assert read_log(log_fn) == ["fail"]


def test_job_entry_point_bash_batches(tmp_path):
    if os.name == "nt":
        pytest.skip("POSIX shell commands")
    out = tmp_path / "out.txt"
    recs = [{"bash": f"echo $GREETING one >> '{out}'"}, {"bash": f"echo two >> '{out}'"}]
    call = Call(
        ["a", "b"],
        ["c", "d"],
        recs,
        {"task_batch_num": 2, "language": "bash", "bash_launch_str": "export GREETING=hi"},
    )
    r = run_exec(write_call(call, tmp_path / "c"), 1)
    assert r.returncode == 0, r.stderr
    assert out.read_text() == "hi one\ntwo\n"


def test_job_entry_point_errors(tmp_path, user_module):
    folder = write_call(
        Call(["a"], ["b"], [{"bash": "true"}], {"language": "bash"}), tmp_path / "c"
    )
    assert run_exec(folder, 2).returncode == 2
    assert run_exec(folder, 0).returncode == 2
    assert run_exec(tmp_path / "nothing", 1).returncode == 2
    r = subprocess.run(
        [sys.executable, "-m", "petakit5d.jobs.exec", str(folder)],
        capture_output=True,
        text=True,
        env=dict(os.environ, PYTHONPATH=SRC_ROOT),
    )
    assert r.returncode == 2 and "usage" in r.stderr
    # a runner named in call.json takes the batch
    runner = f"{user_module}"
    folder = write_call(
        Call(["a"], ["b"], [{"bash": "true"}], {"language": "bash"}),
        tmp_path / "r",
        extra={"runner": "usertasks.jobs:runner"},
    )
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([SRC_ROOT, runner]))
    r = subprocess.run(
        [sys.executable, "-m", "petakit5d.jobs.exec", folder, "1"],
        capture_output=True,
        text=True,
        env=env,
    )
    assert r.returncode == 7 and r.stdout.strip() == "runner 1"


def test_job_entry_point_imports_little(tmp_path):
    # start-up: neither NumPy nor the framework's local loop for a Python batch
    folder = write_call(
        Call(["a"], ["b"], [task_record(Task("os.path:join", ["x"]), "python")], {}), tmp_path / "c"
    )
    code = (
        "import sys, runpy; sys.argv = ['exec', sys.argv[1], '1']\n"
        "try:\n    runpy.run_module('petakit5d.jobs.exec', run_name='__main__')\n"
        "except SystemExit:\n    pass\n"
        "print(sorted(m for m in sys.modules if m.split('.')[0] in "
        "('numpy', 'scipy', 'petakit5d')))"
    )
    r = subprocess.run(
        [sys.executable, "-c", code, folder],
        capture_output=True,
        text=True,
        env=dict(os.environ, PYTHONPATH=SRC_ROOT),
    )
    assert r.returncode == 0, r.stderr
    loaded = eval(r.stdout.strip().splitlines()[-1])
    assert not [m for m in loaded if m.startswith(("numpy", "scipy"))]
    assert "petakit5d.jobs.framework" not in loaded


# --- two callers on the same data -----------------------------------------------------------

CALLER = """
import sys
sys.path.insert(0, sys.argv[1])
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.tasks import Task
ins, outs, log = sys.argv[2].split(","), sys.argv[3].split(","), sys.argv[4]
t = [Task("usertasks.jobs:make", [i, o, 1, log]) for i, o in zip(ins, outs)]
flags = run_tasks(ins, outs, t, parse_cluster=False)
print("FLAGS", all(flags), len(flags))
"""


def test_two_callers_cooperate_through_tmp_markers(tmp_path, user_module):
    # tasks of about 1 s; the second caller starts 0.5 s after the first: each task is left
    # to the caller whose marker is fresh, so every task runs once, and both callers wait
    # for the other's tasks and return all done
    n = 6
    in_fns, out_fns, log_fn = make_tasks(tmp_path, n)
    slow = USER_MODULE + textwrap.dedent("""
        import time
        _make = make

        def make(inp, out, factor, log):
            time.sleep(1.0)
            _make(inp, out, factor, log)
        """)
    (user_module / "usertasks" / "jobs.py").write_text(slow)
    script = tmp_path / "caller.py"
    script.write_text(CALLER)
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    args = [
        sys.executable,
        str(script),
        str(user_module),
        ",".join(in_fns),
        ",".join(out_fns),
        log_fn,
    ]
    p1 = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    time.sleep(0.5)
    p2 = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    outs = [p.communicate(timeout=240) for p in (p1, p2)]
    for p, (out, err) in zip((p1, p2), outs, strict=True):
        assert p.returncode == 0, err
        assert f"FLAGS True {n}" in out, out
    lines = read_log(log_fn)
    runs = [x.split()[0] for x in lines]
    assert sorted(runs) == sorted(out_fns), f"each task should run once:\n{lines}"
    assert len({x.split()[1] for x in lines}) == 2, "both callers should run tasks"
    assert tmp_files(tmp_path / "out") == []


# --- GenericComputingFrameworkDemoTest ------------------------------------------------------

DEMO_MODULE = """
import os
import numpy as np

def crop_flip_demo_function(input_fullpath, output_fullpath, scale, in_bbox):
    # the user function of the demo: skip a done output, compute, write by a temporary
    # name and a rename
    if os.path.exists(output_fullpath):
        return
    im = np.load(input_fullpath)
    y0, x0, z0, y1, x1, z1 = in_bbox
    im = im[y0 - 1:y1, x0 - 1:x1, z0 - 1:z1][::-1] * scale
    tmp = output_fullpath[:-4] + "_" + str(os.getpid()) + ".npy"
    np.save(tmp, im)
    os.replace(tmp, output_fullpath)
"""


def test_user_function_runs_through_framework(tmp_path, user_module):
    # the demo's workflow (demos/demo_generic_computing_framework.m) with a user function in
    # a user module; the demo's crop + deskew/rotate becomes crop + flip until the deskew
    # port exists
    np = pytest.importorskip("numpy")
    (user_module / "usertasks" / "demo.py").write_text(DEMO_MODULE)
    out_path = tmp_path / "Cropped_DSR"
    out_path.mkdir()
    fsns = [
        f"Scan_Iter_0000_0000_Cam{c}_ch0_CAM1_stack0000_488nm_0000000msec_000x.npy" for c in "AB"
    ]
    rng = np.random.default_rng(3)
    ins = [str(tmp_path / f) for f in fsns]
    outs = [str(out_path / f) for f in fsns]
    for fn in ins:
        np.save(fn, rng.integers(0, 1000, (64, 48, 32)).astype(np.uint16))
    bbox = [9, 1, 1, 56, 48, 30]
    t = [
        Task("usertasks.demo:crop_flip_demo_function", [i, o, 2, bbox], path=[str(user_module)])
        for i, o in zip(ins, outs, strict=True)
    ]
    run_tasks(
        ins,
        outs,
        t,
        parse_cluster=False,
        master_compute=True,
        mem_allocate=20,
        config_file="",
        mcc_mode=False,
    )
    for i, o in zip(ins, outs, strict=True):
        expected = np.load(i)[8:56, 0:48, 0:30][::-1] * 2
        np.testing.assert_array_equal(np.load(o), expected)
    assert tmp_files(out_path) == []


def test_demo_config_file_round_trips_and_runs_locally(tmp_path):
    a = {
        "MCRParam": "/opt/MATLAB_Runtime/R2024b",
        "MCCMasterStr": "/opt/PetaKit5D/mcc/linux/run_mccMaster.sh",
        "BashLaunchStr": "module load parallel/20190822",
        "SlurmParam": "-p compute --qos normal -n1 --mem-per-cpu=21418M",
        "jobLogDir": "/mnt/projects/job_logs",
        "maxCPUNum": 24,
        "minCPUNum": 1,
        "wholeNodeJob": True,
        "memPerCPU": 20.9,
        "GNUparallel": True,
        "jobTimeLimit": 12,
        "maxJobNum": 5000,
    }
    config_file = tmp_path / "test_slurm_cpu_config_file.json"
    config_file.write_text(json.dumps(a, indent=2))
    assert json.loads(config_file.read_text()) == a
    # with the config file and no cluster, tasks still run locally
    in_fn = tmp_path / "in.txt"
    in_fn.touch()
    log_fn = str(tmp_path / "calls.log")
    outs = [str(tmp_path / "out_1.txt"), str(tmp_path / "out_2.txt")]
    t = [Task(FUNC, [str(in_fn), o, log_fn, "copy"]) for o in outs]
    flags = run_tasks(
        [str(in_fn)] * 2,
        outs,
        t,
        parse_cluster=False,
        master_compute=True,
        config_file=str(config_file),
        mcc_mode=False,
    )
    assert flags == [True, True]


def test_run_tasks_prints_matlab_progress(tmp_path, capsys):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    out = capsys.readouterr().out
    assert "Task number : 2" in out
    assert f"Process out_1 with function {FUNC}(" in out
    assert "Done! Elapsed time is " in out
    assert "2 / 2 (100.00%) are finished!" in out
    assert "All output files (2 / 2) are finished!" in out
    run_tasks(in_fns, out_fns, task_list(in_fns, out_fns, log_fn, "copy"), **LOCAL)
    assert "All output files (2 / 2) already exist(s)!" in capsys.readouterr().out


def test_task_metrics_hooks(tmp_path, monkeypatch):
    # the seam of the task metrics (J10): task_start(call, indices, mode), task_stop(rec, outcome)
    import types

    events = []
    fake = types.ModuleType("petakit5d.jobs.metrics")
    fake.task_start = lambda call, indices, mode: events.append(("start", indices, mode)) or 7
    fake.task_stop = lambda rec, outcome: events.append(("stop", rec, type(outcome).__name__))
    monkeypatch.setitem(sys.modules, "petakit5d.jobs.metrics", fake)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = task_list(in_fns[:1], out_fns[:1], log_fn, "copy")
    t += task_list(in_fns[1:], out_fns[1:], log_fn, "error")
    with pytest.raises(TaskFailed):
        run_tasks(in_fns, out_fns, t, **LOCAL, task_metrics=True)
    assert events == [
        ("start", [0], "local"),
        ("stop", 7, "NoneType"),
        ("start", [1], "local"),
        ("stop", 7, "TaskFailed"),
    ]
    events.clear()
    x = [str(tmp_path / "x.txt")]
    assert run_tasks(in_fns[:1], x, task_list(in_fns[:1], x, log_fn, "copy"), **LOCAL) == [True]
    assert events == [], "no metrics unless asked for"


def test_call_settings_reach_the_backend(tmp_path, monkeypatch):
    # the seam of the cluster loop (J03): run_cluster gets a Call with the effective settings
    import types

    seen = {}
    fake = types.ModuleType("petakit5d.jobs.cluster")

    def run_cluster(call):
        seen["call"] = call
        return [False] * call.n

    fake.run_cluster = run_cluster
    monkeypatch.setitem(sys.modules, "petakit5d.jobs.cluster", fake)
    monkeypatch.setattr(framework, "_slurm_available", lambda: True)
    config_fn = tmp_path / "c.json"
    config_fn.write_text(json.dumps({"SlurmParam": "-p compute", "jobLogDir": ""}))
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    flags = run_tasks(
        in_fns,
        out_fns,
        task_list(in_fns, out_fns, log_fn, "copy"),
        nice=2.5,
        config_file=str(config_fn),
        task_batch_num=2,
    )
    call = seen["call"]
    assert flags == [False] * 3 and read_log(log_fn) == []
    assert call.settings["slurm_param"] == "-p compute --nice=3"
    assert call.settings["job_log_dir"] == "../job_logs" and call.settings["uuid"]
    assert call.batch_num == 2 and call.records[0]["func"] == FUNC
    # one task without master compute goes to the cluster too; no nice, no --nice
    seen.clear()
    run_tasks(
        in_fns[:1],
        out_fns[:1],
        task_list(in_fns, out_fns, log_fn, "copy")[:1],
        master_compute=False,
        slurm_param="-p compute",
    )
    assert seen["call"].n == 1 and seen["call"].settings["slurm_param"] == "-p compute"
