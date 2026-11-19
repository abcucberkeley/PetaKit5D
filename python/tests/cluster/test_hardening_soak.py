"""Soak of the job framework on a real Slurm under crashes, memory kills, time limits and a
second caller of the same tasks (E06; J11's test_soak covers raised errors, exits, scancel and
outside cancels).

PETAKIT5D_TEST_SLURM_HARDENING=<tasks>x<jobs>[x<max jobs>] (e.g. 400x40x20) runs that many
Python tasks in batches of tasks / jobs, with the task runners on and off (two runs). The
cluster config must give jobs ``--mem-per-cpu=1G`` (``memPerCPU`` 1) and allow 2 CPUs. On the
first run of a task, in a job of the test, failures are injected by batch:

- a crash: the task's process dies of SIGSEGV (batches 1, 9, 17, ...) or SIGABRT (3, 11, ...);
- memory: the task touches 1.3 GiB on every run (5, 13, ...): the 1-CPU job (1 GiB) is killed
  by Slurm's memory limit, the framework's second try has 2 CPUs (2 GiB);
- a time limit: every run that starts in the first 2 minutes sleeps past its job's limit of 3
  minutes (7, 15, ...), so that both callers' jobs of those batches end at their limit.

A second process calls run_tasks with the same tasks, inputs and outputs 20 s after the first
(its own call folder and jobs). At the end both calls are done, every output is complete and was
published once, every planned failure happened, the memory tasks were published by a job of 2
CPUs or more after a 1-CPU run, the time-limit tasks after a run that hung, no call folder is
left and no job of the test is in the queue.
"""

from __future__ import annotations

import collections
import json
import os
import subprocess
import sys
import time

import pytest

pytestmark = [pytest.mark.cluster, pytest.mark.slow]

VARIABLE = "PETAKIT5D_TEST_SLURM_HARDENING"
JOB_TIME_LIMIT = 0.05  # h: 3 minutes
FAILURES = {1: "segv", 3: "abort", 5: "memory", 7: "timeout"}  # batch mod 8 -> failure
SIZE = 32768
SECOND_CALLER_DELAY = 20.0
HANG_WINDOW = 120.0  # s from the start in which time-limit tasks hang

TASKS = '''"""the hardening soak's task"""

import hashlib
import os
import socket
import time
import uuid

from petakit5d.jobs.publish import publish_result, temp_path

SIZE = {size}


def payload(index):
    seed = hashlib.sha256(f"hardening soak {{index}}".encode()).digest()
    return (seed * (SIZE // len(seed) + 1))[:SIZE]


def _note(path, text=""):
    with open(path, "w") as f:
        f.write(text)


def _first_run(folder, index):
    path = os.path.join(folder, "inject", str(index))
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    os.close(fd)
    return True


def work(folder, output, index, failure, job_name, seconds):
    env = os.environ
    job = f"{{env.get('SLURM_ARRAY_JOB_ID', '')}}_{{env.get('SLURM_ARRAY_TASK_ID', '')}}"
    cpus = env.get("SLURM_CPUS_PER_TASK", "")
    run = f"{{index}}_{{uuid.uuid4().hex}}"
    _note(os.path.join(folder, "runs", run), f"{{job}} {{cpus}} {{socket.gethostname()}}\\n")
    in_test_job = env.get("SLURM_JOB_NAME") == job_name and job != "_"
    if failure == "memory" and in_test_job:
        # (every run: only a job of 2 CPUs or more holds it)
        import numpy as np

        if _first_run(folder, index):
            _note(os.path.join(folder, "injected", str(index)), f"memory {{job}} {{cpus}}\\n")
        block = np.ones(int(1.3 * 2**30), dtype=np.uint8)
        time.sleep(2)
        del block
    elif failure == "timeout" and in_test_job:
        with open(os.path.join(folder, "hang_until")) as f:
            hang = time.time() < float(f.read())
        if hang:
            if _first_run(folder, index):
                _note(os.path.join(folder, "injected", str(index)), f"timeout {{job}} {{cpus}}\\n")
            _note(os.path.join(folder, "hung", run), job)
            time.sleep(3600)  # (until the job's limit)
    elif failure and in_test_job and _first_run(folder, index):
        _note(os.path.join(folder, "injected", str(index)), f"{{failure}} {{job}} {{cpus}}\\n")
        if failure == "segv":
            import ctypes

            ctypes.string_at(0)
        if failure == "abort":
            os.abort()
    time.sleep(seconds)
    tmp = temp_path(output)
    with open(tmp, "wb") as f:
        f.write(payload(index))
    if publish_result(tmp, output):
        _note(os.path.join(folder, "published", run), f"{{job}} {{cpus}}")
'''

SECOND_CALLER = """
import json, sys
from petakit5d.jobs.framework import run_tasks
c = json.load(open(sys.argv[1]))
done = run_tasks(c["inputs"], c["outputs"], c["tasks"], **c["settings"])
json.dump(done, open(sys.argv[2], "w"))
"""


def soak_size() -> tuple[int, int, int]:
    spec = os.environ.get(VARIABLE, "")
    if not spec:
        pytest.skip(f"{VARIABLE} is not set (<tasks>x<jobs>[x<max jobs>], e.g. 400x40x20)")
    parts = [int(x) for x in spec.lower().split("x")]
    n, jobs = parts[0], parts[1]
    max_jobs = parts[2] if len(parts) > 2 else 20
    assert n > 0 and jobs > 0 and n % jobs == 0, f"{VARIABLE}={spec}: tasks / jobs"
    return n, jobs, max_jobs


@pytest.mark.parametrize("task_runner", [True, False], ids=["runners", "no_runners"])
def test_soak_with_crashes_memory_kills_time_limits_and_a_second_caller(cluster, task_runner):
    from petakit5d.jobs.tasks import Task, task_record

    n, n_jobs, max_jobs = soak_size()
    size = n // n_jobs
    folder = cluster.shared_folder()
    for sub in ("out", "runs", "inject", "injected", "published", "hung"):
        (folder / sub).mkdir()
    module = f"e06_soak_{folder.name[-8:]}"
    (folder / f"{module}.py").write_text(TASKS.format(size=SIZE))
    sys.path.insert(0, str(folder))
    try:
        payload = __import__(module).payload
    finally:
        sys.path.remove(str(folder))
        sys.modules.pop(module, None)

    planned: dict[int, str] = {}
    outputs, tasks = [], []
    for i in range(1, n + 1):
        b = (i - 1) // size + 1
        failure = FAILURES.get(b % 8) if (i - 1) % size == size // 2 else None
        if failure:
            planned[i] = failure
        out = str(folder / "out" / f"out_{i:05d}.bin")
        args = [str(folder), out, i, failure, cluster.job_name, 0.2 + 0.1 * (i % 4)]
        outputs.append(out)
        tasks.append(Task(f"{module}:work", args, path=[str(folder)]))
    inputs = [str(folder / f"{module}.py")] * n
    params = {
        "task_batch_num": size, "max_trial_num": 3, "max_job_num": max_jobs,
        "job_time_limit": JOB_TIME_LIMIT, "query_interval": 3, "task_runner": task_runner,
        "cpus_per_task": 1,
    }  # fmt: skip

    # the second caller: the same call in another process, 20 s later
    settings = {
        "parse_cluster": True, "master_compute": False, "language": "python",
        "job_log_dir": str(folder / "job_logs"), "config_file": cluster.config_file,
        "slurm_constraint": f"-J {cluster.job_name}", **params,
    }  # fmt: skip
    spec = folder / "second_call.json"
    spec.write_text(json.dumps({"inputs": inputs, "outputs": outputs, "settings": settings,
                                "tasks": [task_record(t, "python") for t in tasks]}))  # fmt: skip
    (folder / "second_caller.py").write_text(SECOND_CALLER)
    second_done = folder / "second_done.json"
    log = open(folder / "second_caller.log", "w")
    env = dict(os.environ)
    (folder / "hang_until").write_text(repr(time.time() + HANG_WINDOW))
    t0 = time.monotonic()
    second = subprocess.Popen(
        ["bash", "-c", f"sleep {SECOND_CALLER_DELAY:g}; exec {sys.executable} "
         f"{folder / 'second_caller.py'} {spec} {second_done}"],
        stdout=log, stderr=subprocess.STDOUT, env=env,
    )  # fmt: skip
    try:
        done, fired = cluster.run_framework(
            folder, inputs, outputs, tasks, language="python", watchdog_minutes=60, **params
        )
        rc = second.wait(timeout=3600)
    finally:
        if second.poll() is None:
            second.terminate()
            second.wait(timeout=120)
        log.close()
    elapsed = time.monotonic() - t0
    left = cluster.wait_until_gone()

    def notes(sub):
        return {p.name: p.read_text().split() for p in (folder / sub).iterdir()}

    runs = notes("runs")
    published = notes("published")
    injected = {int(k): v[0] for k, v in notes("injected").items()}
    runs_of = collections.defaultdict(list)
    for name, v in runs.items():
        runs_of[int(name.split("_")[0])].append(v)
    pub_of = collections.Counter(int(name.split("_")[0]) for name in published)
    complete = [os.path.isfile(o) and open(o, "rb").read() == payload(i)
                for i, o in enumerate(outputs, start=1)]  # fmt: skip
    memory = [i for i, f in planned.items() if f == "memory"]
    hung = collections.Counter(int(p.name.split("_")[0]) for p in (folder / "hung").iterdir())
    timeouts = [i for i, f in planned.items() if f == "timeout"]
    memory_cpus = {i: sorted(int(v[1] or 0) for v in runs_of[i]) for i in memory}
    memory_pub = {i: [int(v[1] or 0) for k, v in published.items() if int(k.split("_")[0]) == i]
                  for i in memory}  # fmt: skip
    second_out = json.loads(second_done.read_text()) if second_done.is_file() else None
    summary = {
        "tasks": n, "jobs": n_jobs, "max_jobs": max_jobs, "runners": task_runner,
        "seconds": round(elapsed), "done": sum(done), "second_caller_exit": rc,
        "second_done": None if second_out is None else sum(second_out),
        "complete": sum(complete), "runs": len(runs),
        "published_twice": sorted(i for i, k in pub_of.items() if k > 1),
        "never_published": n - len(pub_of),
        "planned": collections.Counter(planned.values()),
        "injected": collections.Counter(injected.values()),
        "memory_run_cpus": memory_cpus, "memory_published_cpus": memory_pub,
        "hung_runs": {i: hung[i] for i in timeouts},
        "watchdog": fired, "jobs_left": left,
    }  # fmt: skip
    print(f"hardening soak summary: {json.dumps(summary, default=dict)}")

    assert not fired, "Jobs did not finish within the watchdog's limit and were cancelled."
    assert done == [True] * n, f"{n - sum(done)} tasks not done"
    assert rc == 0 and second_out == [True] * n, "The second caller did not finish its tasks."
    assert all(complete), f"{n - sum(complete)} outputs missing or incomplete"
    assert summary["published_twice"] == [], "Results were published more than once."
    assert summary["never_published"] == 0
    assert injected == planned, "Not every planned failure happened (or more did)."
    for i in memory:
        assert memory_cpus[i][0] == 1, f"task {i} first ran with {memory_cpus[i][0]} CPUs"
        assert memory_pub[i] and min(memory_pub[i]) >= 2, f"task {i}: {memory_pub[i]}"
    assert all(hung[i] >= 1 for i in timeouts), "A time-limit task never hung."
    assert all(len(runs_of[i]) >= 2 for i in planned), "A failed task was not run again."
    assert not list((folder / "job_logs").glob("pk5d_call_*")), "A call folder is left."
    assert left == [], f"Jobs of the test are left in the queue: {left}"
