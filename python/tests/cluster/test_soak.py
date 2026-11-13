"""Soak test of the job framework on a real Slurm, with injected failures (task J11).

PETAKIT5D_TEST_SLURM_SOAK=<tasks>x<jobs>[x<max jobs>] (e.g. 2000x200) runs that many Python
tasks in batches of tasks / jobs, one job per batch, at most <max jobs> (default 40) in the
queue at once, with the task runners (the default) and without master compute. Each task sleeps
a little and publishes a 64 KB result by a temporary name and one rename
(``jobs.publish.publish_result``). Failures are injected on the first run of a task, and only
in a job of the test:

- a non-zero exit: the task raises (batches 1, 11, 21, ...), or its process ends with status 7
  at once (``os._exit``; batches 3, 13, ...);
- a time limit: the task sleeps past its job's limit of 3 minutes (batches 5, 15, ...);
- scancel: the task cancels its own job (batches 7, 17, ...);
- and from outside: the test cancels running jobs of other batches, one every 15 s, at most
  jobs / 20 of them (at least 1), each batch once.

At the end every output is complete (the expected bytes), none was published twice (each
publish_result that published leaves a record), every planned failure happened, the call
folder (claims) is gone and no job of the test is left in the queue.
"""

from __future__ import annotations

import collections
import json
import os
import random
import sys
import threading
import time

import pytest

pytestmark = [pytest.mark.cluster, pytest.mark.slow]

SOAK_VARIABLE = "PETAKIT5D_TEST_SLURM_SOAK"
DEFAULT_MAX_JOBS = 40
JOB_TIME_LIMIT = 0.05  # h: 3 minutes
FAILURES = {1: "raise", 3: "exit", 5: "timeout", 7: "scancel"}  # batch mod 10 -> failure
CANCEL_PERIOD = 15.0  # s between the test's own cancels
SIZE = 65536

SOAK_TASKS = '''"""the soak test's task"""

import hashlib
import os
import socket
import subprocess
import time
import uuid

from petakit5d.jobs.publish import publish_result, temp_path

SIZE = {size}


def payload(index):
    seed = hashlib.sha256(f"soak {{index}}".encode()).digest()
    return (seed * (SIZE // len(seed) + 1))[:SIZE]


def _note(path, text=""):
    with open(path, "w") as f:
        f.write(text)


def _first_run(folder, index):
    """True for the first run of the task (an exclusive create)"""
    try:
        path = os.path.join(folder, "inject", str(index))
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    os.close(fd)
    return True


def work(folder, output, index, failure, job_name, seconds):
    env = os.environ
    job = f"{{env.get('SLURM_ARRAY_JOB_ID', '')}}_{{env.get('SLURM_ARRAY_TASK_ID', '')}}"
    run = f"{{index}}_{{uuid.uuid4().hex}}"
    # one file per run (appends from several nodes are not atomic on network file systems)
    # (and whether the output existed then: a batch submitted again runs all its tasks)
    existed = int(os.path.exists(output))
    note = f"{{existed}} {{job}} {{socket.gethostname()}} {{os.getpid()}}\\n"
    _note(os.path.join(folder, "runs", run), note)
    in_test_job = os.environ.get("SLURM_JOB_NAME") == job_name and job != "_"
    if failure and in_test_job and _first_run(folder, index):
        _note(os.path.join(folder, "injected", str(index)), f"{{failure}} {{job}}\\n")
        if failure == "raise":
            raise RuntimeError(f"task {{index}} fails on purpose")
        if failure == "exit":
            os._exit(7)
        if failure == "scancel":
            subprocess.run(["scancel", job], timeout=60)
        time.sleep(3600)  # (timeout: until the time limit; scancel: until the job ends)
    time.sleep(seconds)
    tmp = temp_path(output)
    data = payload(index)
    with open(tmp, "wb") as f:
        f.write(data[: SIZE // 2])
        f.flush()
        f.write(data[SIZE // 2 :])
    if publish_result(tmp, output):
        _note(os.path.join(folder, "published", run), job)
'''


def soak_size() -> tuple[int, int, int]:
    spec = os.environ.get(SOAK_VARIABLE, "")
    if not spec:
        pytest.skip(f"{SOAK_VARIABLE} is not set (<tasks>x<jobs>, e.g. 2000x200)")
    parts = [int(x) for x in spec.lower().split("x")]
    n, jobs = parts[0], parts[1]
    max_jobs = parts[2] if len(parts) > 2 else DEFAULT_MAX_JOBS
    assert n > 0 and jobs > 0 and n % jobs == 0, f"{SOAK_VARIABLE}={spec}: tasks / jobs"
    return n, jobs, max_jobs


class Canceller:
    """cancels a running job of the test now and then, of a batch without a planned failure,
    each batch at most once"""

    def __init__(self, cluster, batches_without_failure: set[int], limit: int) -> None:
        self.cluster = cluster
        self.allowed = set(batches_without_failure)
        self.limit = limit
        self.cancelled: list[str] = []
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _running(self) -> list[tuple[str, int]]:
        c = self.cluster
        argv = ["squeue", "-h", "-r", "-u", c.user, f"--name={c.job_name}", "-t", "R",
                "-o", "%i|%k"]  # fmt: skip
        r = c.slurm.run(argv)
        if not r.ok:
            return []
        jobs = []
        for line in r.stdout.splitlines():
            key, _, comment = line.strip().partition("|")
            batch = comment.rpartition(":")[2]
            if key and batch.isdigit():
                jobs.append((key, int(batch)))
        return jobs

    def _run(self) -> None:
        rng = random.Random(11)
        while not self._stop.wait(CANCEL_PERIOD) and len(self.cancelled) < self.limit:
            choices = [(k, b) for k, b in self._running() if b in self.allowed]
            if not choices:
                continue
            key, batch = rng.choice(sorted(choices))
            self.allowed.discard(batch)
            if self.cluster.slurm.cancel([key]).ok:
                self.cancelled.append(key)
                print(f"(soak: cancelled job {key} of batch {batch})", flush=True)

    def stop(self) -> None:
        self._stop.set()
        self._thread.join(timeout=CANCEL_PERIOD + 70)


def test_soak_with_injected_failures(cluster):
    from petakit5d.jobs.tasks import Task

    n, n_jobs, max_jobs = soak_size()
    size = n // n_jobs
    folder = cluster.shared_folder()
    for sub in ("out", "runs", "inject", "injected", "published"):
        (folder / sub).mkdir()
    module = f"j11_soak_{folder.name[-8:]}"
    (folder / f"{module}.py").write_text(SOAK_TASKS.format(size=SIZE))
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
        failure = FAILURES.get(b % 10) if (i - 1) % size == size // 2 else None
        if failure:
            planned[i] = failure
        out = str(folder / "out" / f"out_{i:05d}.bin")
        seconds = 0.2 + 0.2 * (i % 5)
        args = [str(folder), out, i, failure, cluster.job_name, seconds]
        outputs.append(out)
        tasks.append(Task(f"{module}:work", args, path=[str(folder)]))
    inputs = [str(folder / f"{module}.py")] * n
    planned_batches = {(i - 1) // size + 1 for i in planned}
    canceller = Canceller(
        cluster, set(range(1, n_jobs + 1)) - planned_batches, max(1, n_jobs // 20)
    )
    tee = cluster.Tee(sys.stdout)
    sys.stdout = tee
    t0 = time.monotonic()
    try:
        done, fired = cluster.run_framework(
            folder, inputs, outputs, tasks, language="python", task_batch_num=size,
            max_trial_num=3, max_job_num=max_jobs, job_time_limit=JOB_TIME_LIMIT,
            query_interval=3, watchdog_minutes=90,
        )  # fmt: skip
    finally:
        sys.stdout = tee.stream
        canceller.stop()
    elapsed = time.monotonic() - t0
    left = cluster.wait_until_gone()

    submitted = cluster.job_ids_of(tee.text())
    run_files = list((folder / "runs").iterdir())
    runs = collections.Counter(int(p.name.split("_")[0]) for p in run_files)
    after_output = sum(p.read_text().startswith("1 ") for p in run_files)
    published = collections.Counter(
        int(p.name.split("_")[0]) for p in (folder / "published").iterdir()
    )
    injected = {int(p.name): p.read_text().split()[0] for p in (folder / "injected").iterdir()}
    complete = [os.path.isfile(o) and open(o, "rb").read() == payload(i)
                for i, o in enumerate(outputs, start=1)]  # fmt: skip
    summary = {
        "tasks": n, "jobs": n_jobs, "max_jobs": max_jobs, "seconds": round(elapsed),
        "submissions": len(submitted), "done": sum(done), "complete": sum(complete),
        "published_twice": sorted(i for i, k in published.items() if k > 1),
        "never_published": n - len(published), "runs": sum(runs.values()),
        "tasks_run_more_than_once": sum(k > 1 for k in runs.values()),
        "runs_after_the_output_existed": after_output,
        "planned_failures": collections.Counter(planned.values()),
        "injected_failures": collections.Counter(injected.values()),
        "cancelled_by_the_test": canceller.cancelled, "watchdog": fired, "jobs_left": left,
    }  # fmt: skip
    print(f"soak summary: {json.dumps(summary, default=dict)}")

    assert not fired, "Jobs did not finish within the watchdog's limit and were cancelled."
    assert done == [True] * n, f"{n - sum(done)} tasks not done"
    assert all(complete), f"{n - sum(complete)} outputs missing or incomplete"
    assert summary["published_twice"] == [], "Results were published more than once."
    assert summary["never_published"] == 0
    assert injected == planned, "Not every planned failure happened (or more did)."
    assert all(runs[i] >= 2 for i in planned), "A task that failed was not run again."
    assert not list((folder / "job_logs").glob("pk5d_call_*")), "The call folder is left."
    assert left == [], f"Jobs of the test are left in the queue: {left}"
