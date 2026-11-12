"""The local process pool of the job framework (jobs.local). Ported from
tests/framework/LocalPoolFrameworkTest.m (8 tests) and the two parfor tests of
GenericComputingFrameworkTest.m; the test names keep the MATLAB names in snake_case
("in this MATLAB" is "in this process").

Tasks call pooltasks.framework_task, a user module written for each test (pool workers
import it by its record, as a cluster job does). It logs every call, as
pktest.frameworkTask does; 'where' writes 'worker' or 'client'.
"""

from __future__ import annotations

import json
import math
import os
import sys
import textwrap
import time
from pathlib import Path

import pytest

from petakit5d.jobs import local
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.tasks import Call, Task, task_record

POOL_TASKS = '''
import json
import os
import sys
import time

from petakit5d.jobs import local


class TaskFailed(Exception):
    """pktest:frameworkTask:failed"""


def framework_task(input_fn, output_fn, log_fn, mode="copy"):
    append_line(log_fn, output_fn)
    on_worker = local.worker_info() is not None
    if mode == "copy":
        write(output_fn, open(input_fn, "rb").read())
    elif mode == "error":
        raise TaskFailed(f"Task failed on purpose for {input_fn}")
    elif mode == "noop":
        pass
    elif mode == "secondRun":
        if open(log_fn).read().splitlines().count(output_fn) >= 2:
            write(output_fn, open(input_fn, "rb").read())
    elif mode in ("where", "clientonly"):
        if on_worker and mode == "clientonly":
            raise TaskFailed(f"Task failed on purpose on a worker for {input_fn}")
        write(output_fn, b"worker" if on_worker else b"client")
    elif mode.startswith("sleep:"):
        time.sleep(float(mode[6:]))
        write(output_fn, open(input_fn, "rb").read())
    else:
        raise ValueError(mode)


def write(fn, data):
    # by a temporary name and one rename, as tasks publish
    tmp = f"{fn}_{os.getpid()}"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, fn)


def append_line(log_fn, line):
    # under a lock folder: on Windows two appends at once can overwrite each other
    lock = log_fn + ".lock"
    t0 = time.monotonic()
    locked = False
    while time.monotonic() - t0 < 30:
        try:
            os.mkdir(lock)
            locked = True
            break
        except OSError:
            time.sleep(0.002)
    for _ in range(1000):
        try:
            with open(log_fn, "a") as fh:
                fh.write(line + "\\n")
            break
        except PermissionError:
            time.sleep(0.002)
    if locked:
        for _ in range(1000):
            try:
                os.rmdir(lock)
                break
            except OSError:
                time.sleep(0.002)


def info(output_fn, mode=""):
    """Writes what the task sees of its process as JSON."""
    w = local.worker_info()
    rec = {
        "pid": os.getpid(),
        "worker": None if w is None else w.id,
        "threads": None if w is None else w.threads,
        "gpu": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "omp": os.environ.get("OMP_NUM_THREADS"),
        "cwd": os.getcwd(),
        "marker": os.environ.get("POOL_TEST_MARKER"),
        "affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
    }
    from petakit5d.core.threads import available_cpus

    rec["available"] = available_cpus()
    if mode == "dirty":
        os.environ["POOL_TEST_MARKER"] = "changed"
        os.chdir(os.path.dirname(output_fn))
    if mode == "sleep":
        time.sleep(0.5)
    write(output_fn, json.dumps(rec).encode())


def compute(input_fn, output_fn, n):
    """Numbers from the input, written as bytes (elementwise math, an FFT, sums)."""
    import numpy as np

    seed = sum(open(input_fn, "rb").read())
    x = np.random.default_rng(seed).standard_normal(n).astype(np.float32)
    y = np.sin(x) * np.float32(1.5) + np.exp(-x * x)
    f = np.fft.rfft(y.astype(np.float64))
    y64 = y.astype(np.float64)
    out = np.concatenate([y64, f.real, f.imag, [np.add.reduce(y64)]])
    write(output_fn, out.tobytes())


def say(output_fn, text):
    print(text)
    sys.stdout.flush()
    os.write(2, (text + " on fd 2\\n").encode())
    write(output_fn, b"x")


def die(output_fn, code):
    os._exit(code)


def leave(output_fn):
    sys.exit(3)


def make_final(output_fn, final_fn):
    write(output_fn, b"x")
    write(final_fn, b"final")
'''


@pytest.fixture
def pooltasks(tmp_path, monkeypatch):
    """The task module, in a folder that is not on sys.path."""
    root = tmp_path / "userlib"
    root.mkdir()
    (root / "pooltasks.py").write_text(POOL_TASKS)
    monkeypatch.setattr(sys, "path", list(sys.path))
    yield str(root)
    sys.modules.pop("pooltasks", None)


def make_tasks(tmp_path: Path, n: int):
    """n input text files, their output paths in out/, the call log."""
    (tmp_path / "out").mkdir()
    in_fns = [str(tmp_path / f"in_{i}.txt") for i in range(1, n + 1)]
    out_fns = [str(tmp_path / "out" / f"out_{i}.txt") for i in range(1, n + 1)]
    for i, fn in enumerate(in_fns, 1):
        Path(fn).write_text(f"task {i}")
    return in_fns, out_fns, str(tmp_path / "calls.log")


def task_list(root, in_fns, out_fns, log_fn, mode):
    return [
        Task("pooltasks:framework_task", [i, o, log_fn, mode], path=[root])
        for i, o in zip(in_fns, out_fns, strict=True)
    ]


def task(root, name, *args):
    return Task(f"pooltasks:{name}", list(args), path=[root])


def read_log(log_fn: str) -> list[str]:
    if not os.path.exists(log_fn):
        return []
    return [x for x in Path(log_fn).read_text().splitlines() if x]


def read_texts(fns):
    return [Path(f).read_text() for f in fns]


LOCAL = {"parse_cluster": False, "master_compute": True}


# --- LocalPoolFrameworkTest -----------------------------------------------------------------


def test_tasks_run_on_the_local_pool(tmp_path, pooltasks):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 4)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "where")
    flags = run_tasks(in_fns, out_fns, t, local_worker_num=2, **LOCAL)
    assert flags == [True] * 4
    assert read_texts(out_fns) == ["worker"] * 4, "The tasks should run on the pool workers."
    assert sorted(read_log(log_fn)) == sorted(out_fns), "Each task should run once."


def test_call_log_keeps_tasks_that_run_at_once(tmp_path, pooltasks):
    # the call log the framework tests count runs with: tasks on two workers append to it at
    # the same time (a plain append lost lines on Windows in MATLAB's CI)
    n = 200
    in_fns, out_fns, log_fn = make_tasks(tmp_path, n)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "noop")
    flags = local.run_on_pool(in_fns, out_fns, t, worker_num=2)
    assert flags == [False] * n
    assert sorted(read_log(log_fn)) == sorted(out_fns), "The log lost calls."
    assert not os.path.isdir(log_fn + ".lock"), "The lock folder was left behind."


def test_failed_tasks_are_retried_in_this_process(tmp_path, pooltasks):
    # the tasks fail on the workers; the framework runs them again
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "clientonly")
    flags = run_tasks(in_fns, out_fns, t, local_worker_num=2, **LOCAL)
    assert flags == [True] * 3
    assert read_texts(out_fns) == ["client"] * 3
    assert sorted(read_log(log_fn)) == sorted(out_fns * 2), (
        "Each task should run once on a worker and once in this process."
    )


def test_failed_tasks_are_reported_unfinished(tmp_path, pooltasks, capfd):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "clientonly")
    flags = local.run_on_pool(in_fns, out_fns, t, worker_num=2)
    assert flags == [False, False]
    assert sorted(read_log(log_fn)) == sorted(out_fns)
    out = capfd.readouterr().out
    assert "Task 1 failed:" in out and "Task 2 failed:" in out
    assert "TaskFailed: Task failed on purpose on a worker" in out


def test_tasks_that_do_not_fit_in_memory_run_in_this_process(tmp_path, pooltasks, monkeypatch):
    # 8 GB of memory: half of it holds no worker with a 10 GB task
    from petakit5d.core import memory

    monkeypatch.setattr(memory, "get_system_memory", lambda: 8.0)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "where")
    flags = run_tasks(in_fns, out_fns, t, local_worker_num=2, mem_allocate=10, **LOCAL)
    assert flags == [True] * 3
    assert read_texts(out_fns) == ["client"] * 3


def test_existing_final_output_skips_all_tasks(tmp_path, pooltasks):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    final_fn = tmp_path / "final.txt"
    final_fn.write_text("")
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "where")
    flags = local.run_on_pool(in_fns, out_fns, t, worker_num=2, final_out_fullpath=str(final_fn))
    assert flags == [True, True]
    assert read_log(log_fn) == []


def test_tasks_with_missing_inputs_are_skipped(tmp_path, pooltasks):
    # as in the framework, which reports them as done
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    os.remove(in_fns[1])
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "where")
    flags = run_tasks(in_fns, out_fns, t, local_worker_num=2, **LOCAL)
    assert flags == [True] * 3
    assert sorted(read_log(log_fn)) == sorted([out_fns[0], out_fns[2]])
    assert not os.path.exists(out_fns[1])


def test_thread_pool_is_not_used(tmp_path, pooltasks, monkeypatch, capfd):
    # MATLAB: thread-based workers can't run mex functions, so no task goes to them. Here
    # tasks never run in threads of this process: when no worker process can start, no task
    # runs on the pool, and the framework runs the tasks in this process
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "where")
    monkeypatch.setattr(local, "_python_executable", lambda: str(tmp_path / "no_python"))
    flags = local.run_on_pool(in_fns, out_fns, t, worker_num=2)
    assert flags == [False, False]
    assert "tasks on the local pool" not in capfd.readouterr().out
    assert read_log(log_fn) == []

    flags = run_tasks(in_fns, out_fns, t, local_worker_num=2, **LOCAL)
    assert flags == [True, True]
    assert read_texts(out_fns) == ["client"] * 2


# --- GenericComputingFrameworkTest, parfor --------------------------------------------------


def test_parfor_wrapper_runs_every_task(tmp_path, pooltasks):
    # through run_tasks: MATLAB's framework passes the parfor wrapper parameters it doesn't
    # take (maxJobNum, taskBatchNum), so this path fails there
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    flags = run_tasks(in_fns, out_fns, t, cluster_type="parfor", local_worker_num=2, gpu_job=False)
    assert flags == [True, True]
    assert read_texts(out_fns) == read_texts(in_fns)


@pytest.mark.parametrize("visible", [None, "4,6"])
def test_parfor_gpu_job_has_a_worker_per_usable_gpu(
    tmp_path, pooltasks, monkeypatch, visible, capfd
):
    # usable_gpu_num is a stub (2 GPUs); each worker gets one of the GPUs this process sees
    monkeypatch.setattr(local, "_usable_gpu_num", lambda: 2)
    if visible is None:
        monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    else:
        monkeypatch.setenv("CUDA_VISIBLE_DEVICES", visible)
    in_fns, out_fns, _ = make_tasks(tmp_path, 6)
    t = [task(pooltasks, "info", o, "sleep") for o in out_fns]
    flags = run_tasks(in_fns, out_fns, t, cluster_type="parfor", gpu_job=True)
    assert flags == [True] * 6
    out = capfd.readouterr().out
    assert out.count(" uses GPU ") == 2, out
    recs = [json.loads(Path(o).read_text()) for o in out_fns]
    gpus = ["0", "1"] if visible is None else ["4", "6"]
    assert {r["worker"] for r in recs} == {1, 2}
    for r in recs:
        assert r["gpu"] == gpus[r["worker"] - 1]


# --- new tests ------------------------------------------------------------------------------


@pytest.mark.parametrize("workers", [2, 3])
def test_outputs_equal_the_sequential_run_bit_for_bit(tmp_path, pooltasks, workers):
    pytest.importorskip("numpy")
    n = 7
    seq_dir = tmp_path / "seq"
    pool_dir = tmp_path / "pool"
    seq_dir.mkdir()
    pool_dir.mkdir()
    in_fns, _, _ = make_tasks(tmp_path, n)
    seq = [str(seq_dir / f"out_{i}.bin") for i in range(n)]
    par = [str(pool_dir / f"out_{i}.bin") for i in range(n)]
    sizes = [1, 2, 17, 1000, 4097, 65536, 3]
    run_tasks(
        in_fns,
        seq,
        [task(pooltasks, "compute", i, o, k) for i, o, k in zip(in_fns, seq, sizes, strict=True)],
        **LOCAL,
    )
    calls = [
        task(pooltasks, "compute", i, o, k) for i, o, k in zip(in_fns, par, sizes, strict=True)
    ]
    flags = run_tasks(in_fns, par, calls, local_worker_num=workers, **LOCAL)
    assert flags == [True] * n
    for a, b in zip(seq, par, strict=True):
        assert Path(a).read_bytes() == Path(b).read_bytes()


def test_workers_share_the_cpus(tmp_path, pooltasks):
    # floor(CPUs / workers) threads each; on Linux each worker has its own block of CPUs
    from petakit5d.core.threads import available_cpus

    in_fns, out_fns, _ = make_tasks(tmp_path, 6)
    t = [task(pooltasks, "info", o, "sleep") for o in out_fns]
    workers = 2
    flags = local.run_on_pool(in_fns, out_fns, t, worker_num=workers)
    assert flags == [True] * 6
    recs = [json.loads(Path(o).read_text()) for o in out_fns]
    threads = max(1, available_cpus() // workers)
    assert {r["threads"] for r in recs} == {threads}
    assert {r["omp"] for r in recs} == {str(threads)}
    assert all(r["pid"] != os.getpid() for r in recs), "a task ran in this process"
    assert len({r["pid"] for r in recs}) == workers
    if hasattr(os, "sched_getaffinity") and len(os.sched_getaffinity(0)) >= threads * workers:
        blocks = {r["worker"]: tuple(r["affinity"]) for r in recs}
        assert all(len(b) == threads for b in blocks.values())
        assert len(set(blocks[1]) & set(blocks[2])) == 0
        assert {r["available"] for r in recs} == {threads}


def test_workers_import_what_this_process_imports(tmp_path, pooltasks):
    # a task named only by its import path (no import roots in the record) runs on the pool
    # if this process can import it: the workers get this process's sys.path
    sys.path.insert(0, pooltasks)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = [
        Task("pooltasks:framework_task", [i, o, log_fn, "where"])
        for i, o in zip(in_fns, out_fns, strict=True)
    ]
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=2) == [True, True]
    assert read_texts(out_fns) == ["worker"] * 2


def test_worker_count_is_at_most_the_task_count(tmp_path, pooltasks, capfd):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=5) == [True, True]
    assert "Run 2 tasks on the local pool (2 at a time)..." in capfd.readouterr().out


def test_tasks_start_in_task_order(tmp_path, pooltasks):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 5)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=1) == [True] * 5
    assert read_log(log_fn) == out_fns


def test_done_tasks_are_not_run(tmp_path, pooltasks):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    Path(out_fns[0]).write_text("before")
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=2) == [True] * 3
    assert sorted(read_log(log_fn)) == sorted(out_fns[1:])
    assert Path(out_fns[0]).read_text() == "before"


def test_missing_inputs_are_not_reported_done_by_the_pool(tmp_path, pooltasks):
    # MATLAB's pool reports only outputs; the framework counts missing inputs as done
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    os.remove(in_fns[0])
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=2) == [False, True, True]
    # every path of a list is looked at (the missing one is not the first)
    lists = [[in_fns[1], in_fns[0]], [in_fns[1]], [in_fns[2], in_fns[1]]]
    assert local.run_on_pool(lists, out_fns, t, worker_num=2) == [False, True, True]
    assert sorted(read_log(log_fn)) == sorted(out_fns[1:]), "a task with a missing input ran"


def test_final_output_made_during_the_run_stops_the_pool(tmp_path, pooltasks, capfd):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 6)
    final_fn = str(tmp_path / "final.txt")
    t = [task(pooltasks, "make_final", out_fns[0], final_fn)]
    t += task_list(pooltasks, in_fns[1:], out_fns[1:], log_fn, "copy")
    flags = local.run_on_pool(in_fns, out_fns, t, worker_num=1, final_out_fullpath=final_fn)
    assert flags == [True] * 6
    assert read_log(log_fn) == [], "no task should start once the final output exists"
    assert f"The final output file {final_fn} already exists!" in capfd.readouterr().out


def test_a_worker_that_ends_fails_its_task_and_is_replaced(tmp_path, pooltasks, capfd):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 4)
    t = [task(pooltasks, "die", out_fns[0], 3)]
    t += task_list(pooltasks, in_fns[1:], out_fns[1:], log_fn, "copy")
    flags = local.run_on_pool(in_fns, out_fns, t, worker_num=1)
    assert flags == [False, True, True, True]
    out = capfd.readouterr().out
    assert "Task 1 failed: the pool worker 1 ended (exit status 3)" in out
    assert "Done! Elapsed time is NaN seconds." in out


def test_sys_exit_in_a_task_is_a_failure_of_the_task(tmp_path, pooltasks, capfd):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    t = [task(pooltasks, "leave", out_fns[0])]
    t += task_list(pooltasks, in_fns[1:], out_fns[1:], log_fn, "copy")
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=1) == [False, True, True]
    out = capfd.readouterr().out
    assert "SystemExit: 3" in out and "started a new one" not in out


def test_each_task_starts_from_the_same_process_state(tmp_path, pooltasks, monkeypatch):
    monkeypatch.setenv("POOL_TEST_MARKER", "original")
    in_fns, out_fns, _ = make_tasks(tmp_path, 3)
    t = [task(pooltasks, "info", o, "dirty") for o in out_fns]
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=1) == [True] * 3
    recs = [json.loads(Path(o).read_text()) for o in out_fns]
    assert len({r["pid"] for r in recs}) == 1
    assert {r["marker"] for r in recs} == {"original"}
    assert {r["cwd"] for r in recs} == {os.getcwd()}


def test_task_output_is_printed_with_its_task(tmp_path, pooltasks, capfd):
    in_fns, out_fns, _ = make_tasks(tmp_path, 2)
    t = [task(pooltasks, "say", o, f"hello {k}") for k, o in enumerate(out_fns, 1)]
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=2) == [True, True]
    out = capfd.readouterr().out
    for k in (1, 2):
        assert f"hello {k}\n" in out and f"hello {k} on fd 2\n" in out
    blocks = out.split("Done! Elapsed time is")
    assert len(blocks) == 3 and all("hello" in b for b in blocks[:2])


def test_the_pool_ends_its_workers(tmp_path, pooltasks):
    in_fns, out_fns, _ = make_tasks(tmp_path, 4)
    t = [task(pooltasks, "info", o) for o in out_fns]
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=2) == [True] * 4
    pids = {json.loads(Path(o).read_text())["pid"] for o in out_fns}
    if os.name == "posix":
        for pid in pids:
            with pytest.raises(ProcessLookupError):
                os.kill(pid, 0)


def test_bash_tasks_run_on_the_pool(tmp_path):
    in_fns, out_fns, _ = make_tasks(tmp_path, 3)
    cmds = [f'echo done > "{o}"' for o in out_fns]
    flags = run_tasks(
        in_fns, out_fns, cmds, cluster_type="parfor", local_worker_num=2, language="bash"
    )
    assert flags == [True] * 3


def test_parfor_retries_unfinished_tasks_in_rounds_on_the_pool(tmp_path, pooltasks):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "secondRun")
    flags = run_tasks(in_fns, out_fns, t, cluster_type="parfor", local_worker_num=2)
    assert flags == [True, True]
    assert sorted(read_log(log_fn)) == sorted(out_fns * 2)
    assert read_texts(out_fns) == read_texts(in_fns)


def test_parfor_never_runs_tasks_in_this_process(tmp_path, pooltasks):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 2)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "clientonly")
    flags = run_tasks(
        in_fns, out_fns, t, cluster_type="parfor", local_worker_num=2, max_trial_num=2
    )
    assert flags == [False, False]
    assert sorted(read_log(log_fn)) == sorted(out_fns * 2)


def test_parfor_reports_missing_inputs_as_done(tmp_path, pooltasks):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    os.remove(in_fns[2])
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    flags = run_tasks(in_fns, out_fns, t, cluster_type="parfor", local_worker_num=2)
    assert flags == [True] * 3
    assert sorted(read_log(log_fn)) == sorted(out_fns[:2])


def test_parfor_gpu_job_without_a_usable_gpu_uses_cpu_workers(
    tmp_path, pooltasks, monkeypatch, capfd
):
    monkeypatch.setattr(local, "_usable_gpu_num", lambda: 0)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "7")
    in_fns, out_fns, _ = make_tasks(tmp_path, 2)
    t = [task(pooltasks, "info", o) for o in out_fns]
    flags = run_tasks(in_fns, out_fns, t, cluster_type="parfor", gpu_job=True, local_worker_num=2)
    assert flags == [True, True]
    assert "no usable GPU" in capfd.readouterr().err
    assert {json.loads(Path(o).read_text())["gpu"] for o in out_fns} == {"7"}


def test_parfor_workers_by_default_are_the_cpus_at_most_max_job_num(tmp_path, pooltasks, capfd):
    from petakit5d.core.threads import available_cpus

    in_fns, out_fns, log_fn = make_tasks(tmp_path, 40)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    assert run_tasks(in_fns, out_fns, t, cluster_type="parfor") == [True] * 40
    n = min(available_cpus(), 40)
    assert f"({n} at a time)" in capfd.readouterr().out
    for o in out_fns:
        os.remove(o)
    assert run_tasks(in_fns, out_fns, t, cluster_type="parfor", max_job_num=1) == [True] * 40
    assert "(1 at a time)" in capfd.readouterr().out


def test_run_pool_follows_the_call_settings(tmp_path, pooltasks, monkeypatch):
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 3)
    recs = [
        task_record(x, "python") for x in task_list(pooltasks, in_fns, out_fns, log_fn, "where")
    ]
    call = Call(in_fns, out_fns, recs, {"local_worker_num": 2, "mem_allocate": None})
    assert local.run_pool(call) == [True] * 3
    assert read_texts(out_fns) == ["worker"] * 3
    from petakit5d.core import memory

    monkeypatch.setattr(memory, "get_system_memory", lambda: 47.9)
    call.settings["mem_allocate"] = 10  # floor(47.9 / 2 / 12) = 1 worker: no pool
    assert local.run_pool(call, worker_num=4) is None


@pytest.mark.parametrize(
    "worker_num, mem_allocate, memory, expected",
    [
        (2, None, 8.0, 2),
        (4, 10, 48.0, 2),  # floor(48 / 2 / 12)
        (4, 10, 47.9, 1),
        (2, 10, 8.0, 0),
        (4, 10, math.nan, 4),  # memory unknown: no limit (MATLAB's min ignores NaN)
        (3, 0, 100.0, 3),
        (30, 0, 100.0, 25),
        (2.5, 1, 1000.0, 2.5),
    ],
)
def test_pool_worker_num_follows_matlab(worker_num, mem_allocate, memory, expected):
    assert local.pool_worker_num(worker_num, mem_allocate, memory) == expected


@pytest.mark.parametrize(
    "worker_num, n_tasks, expected",
    [
        (2, 5, 2),
        (2.5, 5, 3),
        (2.4999, 5, 2),
        (0, 5, 1),
        (math.nan, 5, 1),
        (math.inf, 5, 5),
        (9, 3, 3),
    ],
)
def test_worker_count_rounds_as_matlab(worker_num, n_tasks, expected):
    assert local._worker_count(worker_num, n_tasks) == expected


def test_worker_info_is_none_outside_a_pool(monkeypatch):
    monkeypatch.delenv(local.ENV_WORKER, raising=False)
    assert local.worker_info() is None


def test_a_pool_of_quick_tasks_starts_and_ends_quickly(tmp_path, pooltasks):
    # start-up and shutdown of a pool of 4 workers (spawn on every platform)
    in_fns, out_fns, log_fn = make_tasks(tmp_path, 8)
    t = task_list(pooltasks, in_fns, out_fns, log_fn, "copy")
    t0 = time.perf_counter()
    assert local.run_on_pool(in_fns, out_fns, t, worker_num=4) == [True] * 8
    assert time.perf_counter() - t0 < 60


def test_worker_imports_only_what_tasks_need(tmp_path):
    # a worker starts with the standard library and petakit5d.jobs (no NumPy)
    code = textwrap.dedent("""
        import sys
        from petakit5d.jobs import local
        print("numpy" in sys.modules)
    """)
    import subprocess

    from petakit5d.jobs.local import _package_root

    env = dict(os.environ, PYTHONPATH=_package_root())
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.stdout.strip() == "False", r.stderr
