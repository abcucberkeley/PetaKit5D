"""In-job slots (petakit5d.jobs.slots), which replace GNU parallel in jobs of whole nodes:
one process per task, up to P at once, starts 0.2 s apart, two runs of a failing task, a GPU
per slot, tasks whose outputs exist skipped, and the extra tasks of runExtraTasks in MATLAB's
order."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import pytest

import petakit5d
from petakit5d.jobs import slots
from petakit5d.jobs.slots import batch_tasks, extra_batches, gpu_for_slot, run_batch, run_slots
from petakit5d.jobs.tasks import Call, Task, task_record, write_call

SRC_ROOT = str(Path(petakit5d.__file__).resolve().parents[1])

TASKS = """
import json
import os
import time


def record(path, seconds=0.0, fail_first=False):
    tries = path + ".tries"
    n = int(open(tries).read()) + 1 if os.path.exists(tries) else 1
    with open(tries, "w") as f:
        f.write(str(n))
    start = time.time()
    time.sleep(seconds)
    if fail_first and n == 1:
        raise RuntimeError("first try fails on purpose")
    info = {
        "start": start,
        "end": time.time(),
        "slot": os.environ.get("PARALLEL_JOBSLOT"),
        "seq": os.environ.get("PARALLEL_SEQ"),
        "gpu": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "pid": os.getpid(),
    }
    tmp = f"{path}_{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(info, f)
    os.replace(tmp, path)


def fail(path):
    tries = path + ".tries"
    n = int(open(tries).read()) + 1 if os.path.exists(tries) else 1
    with open(tries, "w") as f:
        f.write(str(n))
    raise RuntimeError("always fails")
"""


@pytest.fixture
def task_module(tmp_path, monkeypatch):
    root = tmp_path / "userlib"
    root.mkdir()
    (root / "j03_slot_tasks.py").write_text(TASKS)
    monkeypatch.setenv("PYTHONPATH", SRC_ROOT)
    return root


def make_call(tmp_path, root, specs, task_batch_num=None, language="python", **settings):
    """a call folder of tasks: specs are (function, kwargs) per task"""
    (tmp_path / "out").mkdir(exist_ok=True)
    outs = [str(tmp_path / "out" / f"out_{i}.json") for i in range(1, len(specs) + 1)]
    if language == "bash":
        recs = [
            task_record(cmd.format(out=o), "bash") for (cmd, _), o in zip(specs, outs, strict=True)
        ]
    else:
        recs = [
            task_record(Task(f"j03_slot_tasks:{fn}", [o], kw, path=[str(root)]), "python")
            for (fn, kw), o in zip(specs, outs, strict=True)
        ]
    extra = settings.pop("extra", None)
    s = {"task_batch_num": task_batch_num or len(specs), "language": language}
    s.update(settings)
    call = Call(["in"] * len(specs), outs, recs, s)
    folder = write_call(call, tmp_path / "call", extra)
    return folder, outs


def info(path):
    return json.loads(Path(path).read_text())


def tries(path):
    p = Path(path + ".tries")
    return int(p.read_text()) if p.exists() else 0


def test_tasks_run_in_slots_at_most_p_at_once(tmp_path, task_module):
    folder, outs = make_call(tmp_path, task_module, [("record", {"seconds": 1.0})] * 6)
    assert run_slots(folder, range(6), 3) == 0
    runs = [info(o) for o in outs]
    # never more than 3 at once, each task in its own process, slots 1..3
    for r in runs:
        at_once = sum(x["start"] <= r["start"] < x["end"] for x in runs)
        assert at_once <= 3
    assert len({r["pid"] for r in runs}) == 6
    assert {r["slot"] for r in runs} == {"1", "2", "3"}
    assert [r["seq"] for r in runs] == [str(i) for i in range(1, 7)]


def test_starts_are_slot_delay_apart(tmp_path, task_module, monkeypatch):
    # with free slots, each start after the first waits until 0.2 s after the one before
    # (parallel --delay 0.2)
    folder, outs = make_call(tmp_path, task_module, [("record", {})] * 4)
    waits = []

    def sleep(seconds):
        waits.append(seconds)
        time.sleep(seconds)

    monkeypatch.setattr(slots, "_sleep", sleep)
    t0 = datetime.now().timestamp()
    run_slots(folder, range(4), 4)
    assert datetime.now().timestamp() - t0 >= 3 * 0.2 * 0.95
    assert len(waits) == 3 and all(0.1 < w <= 0.2 for w in waits)


def test_a_failing_task_runs_twice(tmp_path, task_module):
    folder, outs = make_call(
        tmp_path, task_module, [("record", {"fail_first": True}), ("fail", {}), ("record", {})]
    )
    assert run_slots(folder, range(3), 2) == 1
    assert Path(outs[0]).exists() and tries(outs[0]) == 2
    assert not Path(outs[1]).exists() and tries(outs[1]) == 2
    assert tries(outs[2]) == 1
    # one run each without retries
    for o in outs:
        Path(o + ".tries").unlink(missing_ok=True)
        Path(o).unlink(missing_ok=True)
    assert run_slots(folder, range(3), 2, tries=1) == 2
    assert [tries(o) for o in outs] == [1, 1, 1]


def test_tasks_whose_outputs_exist_are_skipped(tmp_path, task_module, capfd):
    folder, outs = make_call(tmp_path, task_module, [("record", {})] * 3)
    Path(outs[1]).write_text("{}")
    assert run_slots(folder, range(3), 2) == 0
    assert tries(outs[1]) == 0 and tries(outs[0]) == tries(outs[2]) == 1
    assert f'output "{outs[1]}" already exists' in capfd.readouterr().out


def test_nothing_runs_once_the_final_output_exists(tmp_path, task_module):
    final = tmp_path / "final.txt"
    folder, outs = make_call(
        tmp_path, task_module, [("record", {})] * 2, final_out_fullpath=str(final)
    )
    final.write_text("done")
    assert run_slots(folder, range(2), 2) == 0
    assert [tries(o) for o in outs] == [0, 0]


def test_gpu_slots_get_one_gpu_each(tmp_path, task_module, monkeypatch, capfd):
    folder, outs = make_call(tmp_path, task_module, [("record", {"seconds": 0.5})] * 4)

    def seen():
        got = {(info(o)["slot"], info(o)["gpu"]) for o in outs}
        for o in outs:
            Path(o).unlink()
        return got

    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "4,6")
    run_slots(folder, range(4), 2, gpu=True)
    assert seen() == {("1", "4"), ("2", "6")}
    # fewer devices listed than slots: fewer slots, one of the job's devices each, never
    # another job's GPU (MATLAB since its fix)
    capfd.readouterr()
    run_slots(folder, range(4), 3, gpu=True)
    assert seen() == {("1", "4"), ("2", "6")}
    assert "GPUs 4,6: 2 slots" in capfd.readouterr().out
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "5")
    run_slots(folder, range(4), 3, gpu=True)
    assert seen() == {("1", "5")}
    assert "GPUs 5: 1 slots" in capfd.readouterr().out
    # without a list: devices 0 to P - 1
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES")
    run_slots(folder, range(4), 3, gpu=True)
    assert seen() == {("1", "0"), ("2", "1"), ("3", "2")}


def test_gpu_for_slot():
    assert gpu_for_slot(1, 2, "4,6") == "4" and gpu_for_slot(2, 2, "4, 6") == "6"
    assert gpu_for_slot(1, 1, None) == "0" and gpu_for_slot(2, 2, "") == "1"
    assert gpu_for_slot(2, 2, "GPU-a,GPU-b") == "GPU-b"
    with pytest.raises(ValueError):
        gpu_for_slot(2, 2, "3")
    assert [slots.gpu_slots(3, v) for v in ("4,6", "5", None, "", "1,2,3,4")] == [2, 1, 3, 3, 3]


def test_bash_tasks_run_in_slots_with_bash_launch_str(tmp_path, task_module):
    specs = [("echo $GREETING $PARALLEL_JOBSLOT > {out}", None)] * 2
    folder, outs = make_call(
        tmp_path, task_module, specs, language="bash", bash_launch_str="export GREETING=hi"
    )
    assert run_slots(folder, range(2), 1) == 0
    assert [Path(o).read_text().split() for o in outs] == [["hi", "1"], ["hi", "1"]]


def test_extra_batches_in_matlabs_order():
    # mcc_slurm_cluster_generic_computing_wrapper: the later batches from the last down,
    # reordered b by b, the first max(2, ceil(500 / batch size)) of them for b > 1
    assert extra_batches(1, 5, False, 4) == [5, 4, 3, 2]
    assert extra_batches(2, 6, False, 4) == [5, 3, 4]
    assert extra_batches(2, 6, True, 4) == [5, 3, 4, 6]
    assert extra_batches(3, 12, False, 300) == [10, 7]
    assert extra_batches(3, 12, False, 10) == [10, 7, 4, 9, 6, 8, 5]
    assert extra_batches(5, 5, False, 4) == []
    assert extra_batches(2, 3, True, 4) == [3]
    # batches over 500 tasks: still 2 extra batches (MATLAB's task files, extra_tasks_big_batches)
    assert extra_batches(2, 6, False, 600) == [5, 3]
    assert extra_batches(1, 6, False, 600) == [6, 5, 4, 3, 2]


def test_batch_tasks_with_extra_tasks():
    meta = {"n_tasks": 10, "task_batch_num": 4}
    assert batch_tasks(meta, 3) == [8, 9]
    meta["slots"] = {"extra_tasks": True, "outputs_existed": False}
    # batch 1: its tasks, then batches 3 and 2, each last task first
    assert batch_tasks(meta, 1) == [0, 1, 2, 3, 9, 8, 7, 6, 5, 4]
    assert batch_tasks(meta, 2) == [4, 5, 6, 7]


def test_run_batch_runs_the_batch_and_its_extra_tasks(tmp_path, task_module):
    folder, outs = make_call(
        tmp_path,
        task_module,
        [("record", {})] * 5,
        task_batch_num=2,
        extra={"slots": {"extra_tasks": True, "outputs_existed": False}},
    )
    assert run_batch(folder, 1, 2) == 0
    assert all(Path(o).exists() for o in outs)
    # (PARALLEL_SEQ: the place in the batch's list, own tasks first, batch 3 then batch 2)
    assert [info(o)["seq"] for o in outs] == ["1", "2", "5", "4", "3"]


def run_cli(*args, env=None):
    return subprocess.run(
        [sys.executable, "-m", "petakit5d.jobs.slots", *map(str, args)],
        capture_output=True,
        text=True,
        env=env,
        timeout=120,
    )


def test_job_entry_point(tmp_path, task_module):
    folder, outs = make_call(
        tmp_path, task_module, [("record", {}), ("fail", {}), ("record", {})], task_batch_num=3
    )
    r = run_cli(folder, 1, "--slots", 2, "--tries", 1)
    assert r.returncode == 1, r.stderr
    assert "always fails" in r.stderr and Path(outs[0]).exists() and Path(outs[2]).exists()
    # one task in this process
    Path(outs[0]).unlink()
    assert run_cli("--task", folder, 0).returncode == 0 and Path(outs[0]).exists()
    assert run_cli("--task", folder, 1).returncode == 1
    # wrong arguments
    assert run_cli().returncode == 2
    assert run_cli(folder).returncode == 2
    assert run_cli(folder, 1, "--slots").returncode == 2
    assert run_cli(folder, "x").returncode == 2
    assert run_cli(tmp_path / "nothing", 1).returncode == 2


def test_the_job_entry_point_imports_little(tmp_path, task_module):
    code = (
        "import sys\n"
        "import petakit5d.jobs.slots\n"
        "import petakit5d.jobs.tasks\n"
        "assert 'numpy' not in sys.modules, sorted(sys.modules)\n"
        "print('ok')"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    assert r.stdout.strip() == "ok", r.stderr


def test_task_processes_start_fast(tmp_path, task_module):
    # job start to task start: the slot's process start-up (the gate: job overhead under 1 s)
    folder, outs = make_call(tmp_path, task_module, [("record", {})])
    t0 = datetime.now().timestamp()
    run_slots(folder, [0], 1)
    assert info(outs[0])["start"] - t0 < 1.0
