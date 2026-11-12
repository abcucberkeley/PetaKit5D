"""Task runners (jobs.runner): a job's process runs its batch, then unstarted tasks of other
batches, with the process state given back after each task and guards that start a new
process; the master side that the cluster loop asks (MATLAB's TaskRunnerTest).

The job side runs through the real job entry point (``python -m petakit5d.jobs.exec``) on call
folders written as the cluster loop writes them; the master side runs against stubs of Slurm
and, through the cluster loop, the fake Slurm of the MATLAB suite."""

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

from petakit5d.jobs import claims, runner
from petakit5d.jobs.claims import Claims
from petakit5d.jobs.runner import RUNNER, ClusterRunners, batch_keys
from petakit5d.jobs.slurm import CommandResult
from petakit5d.jobs.tasks import Call, Task, task_record, write_call
from petakit5d.testing.slurm_fake import SLURM_JOB_VARIABLES, FakeSlurm

pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="task runners run on Linux only"
)

SRC_ROOT = str(Path(runner.__file__).resolve().parents[2])

TASKS = '''
import json
import os
import threading
import time
import uuid
import warnings

TOKEN = uuid.uuid4().hex  # one per interpreter: a new process has another
KEPT = []


def work(out, log, mode="ok", seconds=0.0):
    """appends what it sees to log, then writes out (by a temporary name and a rename)"""
    rec = {"out": out, "pid": os.getpid(), "token": TOKEN, "mode": mode,
           "job": os.environ.get("SLURM_ARRAY_JOB_ID", "")}
    with open(log, "a") as fh:
        fh.write(json.dumps(rec) + "\\n")
    if mode == "error":
        raise RuntimeError("on purpose")
    if mode == "hold":
        KEPT.append(bytearray(b"x") * int(0.4 * 2**30))
    if mode == "files":
        KEPT.extend(open(os.devnull) for _ in range(100))
    if mode == "threads":
        stop = threading.Event()
        KEPT.append(stop)
        for _ in range(300):
            threading.Thread(target=stop.wait, daemon=True).start()
    if mode == "exit0":
        raise SystemExit(0)
    if seconds:
        time.sleep(seconds)
    tmp = f"{out}_{uuid.uuid4().hex}"
    with open(tmp, "w") as fh:
        fh.write(mode)
    os.replace(tmp, out)


def leak(out, log):
    """changes what a task can change in its process"""
    import numpy as np
    import signal
    import sys

    os.chdir("/")
    os.environ["PK_LEAKED"] = "1"
    os.environ.pop("PK_KEPT", None)
    warnings.simplefilter("error")
    np.set_printoptions(precision=2)
    np.seterr(all="raise")
    sys.path.insert(0, "/nowhere")
    os.umask(0o077)
    signal.signal(signal.SIGUSR1, lambda *a: None)
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    with open(out, "w") as fh:
        fh.write("leaked")


def probe(out, log):
    """writes what it sees of the process state"""
    import numpy as np
    import signal
    import sys

    from petakit5d.core.threads import available_cpus

    seen = {
        "cwd": os.getcwd(),
        "leaked": os.environ.get("PK_LEAKED"),
        "kept": os.environ.get("PK_KEPT"),
        "state": os.environ.get("PETAKIT5D_RUNNER_STATE"),
        "error filter": any(f[0] == "error" for f in warnings.filters),
        "print": np.get_printoptions()["precision"],
        "err": np.geterr()["divide"],
        "path0": sys.path[0],
        "umask": os.umask(0o022),
        "usr1": str(signal.getsignal(signal.SIGUSR1)),
        "cpus": available_cpus(),
    }
    os.umask(seen["umask"])
    with open(out, "w") as fh:
        json.dump(seen, fh)
    with open(log, "a") as fh:
        fh.write(json.dumps({"out": out, "pid": os.getpid(), "token": TOKEN}) + "\\n")
'''


@pytest.fixture
def no_job(monkeypatch):
    for name in SLURM_JOB_VARIABLES + ("SLURM_CPUS_PER_TASK", "PETAKIT5D_RUNNER_STATE"):
        monkeypatch.delenv(name, raising=False)
    for name in (runner.MEMORY_VARIABLE, runner.MINUTES_VARIABLE):
        monkeypatch.delenv(name, raising=False)


class Calls:
    """a call folder as the cluster loop writes it, with fake tasks (module rtasks)"""

    def __init__(
        self, tmp_path, n, bs, modes=None, keys=None, skip=None, settings=None, seconds=0.0
    ):
        self.root = tmp_path
        tmp_path.mkdir(parents=True, exist_ok=True)
        (tmp_path / "rtasks.py").write_text(TASKS)
        (tmp_path / "out").mkdir(exist_ok=True)
        self.log = str(tmp_path / "runs.log")
        self.outs = [str(tmp_path / "out" / f"o_{k}.txt") for k in range(1, n + 1)]
        ins = [str(tmp_path / f"in_{k}.txt") for k in range(1, n + 1)]
        for p in ins:
            Path(p).write_text("x")
        modes = modes or {}
        recs = []
        for k, o in enumerate(self.outs, 1):
            m = modes.get(k, "ok")
            func, args = ("rtasks:work", [o, self.log, m, seconds])
            if m in ("leak", "probe"):
                func, args = (f"rtasks:{m}", [o, self.log])
            recs.append(task_record(Task(func, args, path=[str(tmp_path)]), "python"))
        s = {"task_batch_num": bs, "language": "python"}
        s.update(settings or {})
        call = Call(ins, self.outs, recs, s)
        nb = -(-n // bs)
        cfg = {"version": 1, "keys": keys or batch_keys(nb, 1, False), "skip": skip or []}
        self.folder = write_call(call, tmp_path / "call", {"runner": RUNNER, "task_runner": cfg})
        self.claims = Claims(self.folder)

    def job(self, batch, env=None, timeout=180):
        e = dict(os.environ, PYTHONPATH=SRC_ROOT)
        e.update(env or {})
        return subprocess.run(
            [sys.executable, "-u", "-m", "petakit5d.jobs.exec", self.folder, str(batch)],
            capture_output=True,
            text=True,
            env=e,
            timeout=timeout,
        )

    def runs(self):
        if not os.path.exists(self.log):
            return []
        return [json.loads(x) for x in Path(self.log).read_text().splitlines() if x.strip()]

    def tasks_run(self):
        return [self.outs.index(r["out"]) + 1 for r in self.runs()]

    def names(self):
        return sorted(self.claims.names())


def tnames(ks):
    return sorted(f"t{k}" for k in ks)


# --- the job side ------------------------------------------------------------------------------


def test_batch_keys():
    assert batch_keys(3, 4, False, 3) == ["4", "4", "master"]
    assert batch_keys(2, 8, True, 0) == ["8gpu", "8gpu"]


def test_job_runs_each_task_once_in_one_process(tmp_path, no_job):
    c = Calls(tmp_path, 8, 8)
    r = c.job(1)
    assert r.returncode == 0, r.stdout + r.stderr
    assert sorted(c.tasks_run()) == list(range(1, 9))
    assert all(os.path.isfile(o) for o in c.outs)
    assert len({x["token"] for x in c.runs()}) == 1, "each task ran in a process of its own"
    assert c.names() == tnames(range(1, 9))
    assert r.stdout.strip().splitlines()[-1].startswith("Elapsed time is ")


def test_job_takes_unstarted_tasks_of_other_batches_last_first(tmp_path, no_job):
    c = Calls(tmp_path, 12, 4)
    r = c.job(1)
    assert r.returncode == 0, r.stderr
    # its batch in order, then the others from the last task backwards
    assert c.tasks_run() == [1, 2, 3, 4, 12, 11, 10, 9, 8, 7, 6, 5]
    assert "task 12 of batch 3" in r.stdout
    # (the job's stealer mark is only made in a Slurm job)
    assert not any(n.startswith("s") for n in c.names())


def test_job_takes_only_tasks_of_batches_with_its_key(tmp_path, no_job):
    c = Calls(tmp_path, 16, 4, keys=["4", "4", "master", "8"])
    assert c.job(1).returncode == 0
    assert sorted(c.tasks_run()) == list(range(1, 9))
    # the master's batch: its own tasks only
    c = Calls(tmp_path / "m", 16, 4, keys=["4", "4", "master", "8"])
    assert c.job(3).returncode == 0
    assert sorted(c.tasks_run()) == [9, 10, 11, 12]


def test_job_with_other_cpus_than_its_key_takes_only_its_batch(tmp_path, no_job):
    c = Calls(tmp_path, 8, 4, keys=["4", "4"])
    assert c.job(1, {"SLURM_CPUS_PER_TASK": "8"}).returncode == 0
    assert sorted(c.tasks_run()) == [1, 2, 3, 4]


def test_job_leaves_claimed_tasks_resubmitted_batches_and_outputs(tmp_path, no_job):
    c = Calls(tmp_path, 12, 4)
    c.claims.claim(2, "99_1:othernode:1:1")
    Path(c.outs[3]).write_text("x\n")
    c.claims.mark_resubmitted(2)
    r = c.job(1)
    assert r.returncode == 0, r.stderr
    assert sorted(c.tasks_run()) == [1, 2, 9, 10, 11, 12]
    assert f'output "{c.outs[3]}" already exists' in r.stdout
    assert Path(c.outs[3]).read_text() == "x\n", "the existing output was rewritten"


def test_tasks_without_inputs_stay_in_their_batch(tmp_path, no_job):
    # a task whose input was missing when the call started runs only in its own batch's job
    c = Calls(tmp_path, 4, 2, skip=[3])
    assert c.job(1).returncode == 0
    assert sorted(c.tasks_run()) == [1, 2, 3]


def test_failed_task_is_marked_and_the_job_takes_no_more_of_other_batches(tmp_path, no_job):
    c = Calls(tmp_path, 12, 8, modes={2: "error"})
    r = c.job(1)
    assert r.returncode == 1
    assert "RuntimeError: on purpose" in r.stderr
    runs = c.runs()
    assert sorted(c.tasks_run()) == list(range(1, 9)), "own tasks run, no other batch's"
    assert "f2" in c.names() and "f1" not in c.names()
    # the tasks after the failed one ran in a new process (the same pid: exec)
    assert runs[0]["token"] == runs[1]["token"] != runs[2]["token"]
    assert len({x["pid"] for x in runs}) == 1
    assert "a new process runs task 3" in r.stdout


def test_job_with_no_task_left_runs_none(tmp_path, no_job):
    c = Calls(tmp_path, 4, 4)
    for i in range(4):
        c.claims.claim(i, "99_1:othernode:1:1")
    r = c.job(1)
    assert r.returncode == 0
    assert c.runs() == []
    assert "no task left" in r.stdout


def test_stealing_stops_near_the_jobs_end(tmp_path, no_job):
    now = int(time.time())
    env = {"SLURM_JOB_START_TIME": str(now - 3600), "SLURM_JOB_END_TIME": str(now + 30)}
    c = Calls(tmp_path, 8, 4)
    assert c.job(1, env).returncode == 0
    assert sorted(c.tasks_run()) == [1, 2, 3, 4]
    env["SLURM_JOB_END_TIME"] = str(now + 3600)
    c = Calls(tmp_path / "later", 8, 4)
    assert c.job(1, env).returncode == 0
    assert sorted(c.tasks_run()) == list(range(1, 9))


def test_time_left():
    os.environ.pop("SLURM_JOB_END_TIME", None)
    assert runner.time_left_ok(1e6)
    env = {"SLURM_JOB_START_TIME": "1000", "SLURM_JOB_END_TIME": "2000"}
    with pytest.MonkeyPatch.context() as mp:
        for k, v in env.items():
            mp.setenv(k, v)
        assert runner.time_left_ok(10, now=1879) and not runner.time_left_ok(10, now=1921)
        # (before any task: less than half of the job's time used)
        assert runner.time_left_ok(0, now=1499) and not runner.time_left_ok(0, now=1501)


def test_runner_resets_the_process_between_tasks(tmp_path, no_job):
    modes = {1: "probe", 2: "leak", 3: "probe", 4: "leak", 5: "probe"}
    c = Calls(tmp_path, 5, 5, modes=modes)
    r = c.job(1, {"PK_KEPT": "yes"})
    assert r.returncode == 0, r.stdout + r.stderr
    seen = [json.loads(Path(c.outs[k]).read_text()) for k in (0, 2, 4)]
    assert seen[0] == seen[1] == seen[2], seen
    first = seen[0]
    assert first["leaked"] is None and first["kept"] == "yes" and first["state"] is None
    assert first["print"] == 8 and first["err"] == "warn" and first["path0"] != "/nowhere"
    assert len({x["token"] for x in c.runs()}) == 1, "the tasks did not share a process"
    # (NumPy's print options and error state stay in the task's own context)
    line = next(x for x in r.stdout.splitlines() if "task 2 changed" in x)
    for item in ("folder", "environment", "path", "warnings", "umask", "signal"):
        assert item in line, line


def test_memory_guard_after_the_first_task(tmp_path, no_job):
    c = Calls(tmp_path, 3, 3, modes={1: "hold"})
    r = c.job(1, {runner.MEMORY_VARIABLE: "0.1"})
    assert r.returncode == 0, r.stderr
    t = [x["token"] for x in c.runs()]
    assert t[0] != t[1] == t[2], r.stdout
    why = r"memory [\d.]+ GB \([\d.]+ GB at its start\), a new process runs task 2"
    assert re.search(why, r.stdout), r.stdout


def test_memory_guard_after_a_later_task(tmp_path, no_job):
    c = Calls(tmp_path, 3, 3, modes={2: "hold"})
    r = c.job(1, {runner.MEMORY_VARIABLE: "0.1"})
    assert r.returncode == 0, r.stderr
    t = [x["token"] for x in c.runs()]
    assert t[0] == t[1] != t[2], r.stdout


def test_runner_ends_after_its_time(tmp_path, no_job):
    c = Calls(tmp_path, 3, 3)
    r = c.job(1, {runner.MINUTES_VARIABLE: "0"})
    assert r.returncode == 0, r.stderr
    assert len({x["token"] for x in c.runs()}) == 3
    assert sorted(c.tasks_run()) == [1, 2, 3]


def test_tasks_never_see_the_runner_state(tmp_path, no_job):
    # (the state a process leaves for the next one of the job is not the tasks')
    c = Calls(tmp_path, 3, 3, modes={1: "probe", 2: "probe", 3: "probe"})
    r = c.job(1, {runner.MINUTES_VARIABLE: "0"})
    assert r.returncode == 0, r.stderr
    assert len({x["token"] for x in c.runs()}) == 3
    assert [json.loads(Path(o).read_text())["state"] for o in c.outs] == [None] * 3


def test_process_with_many_more_open_files_is_renewed(tmp_path, no_job):
    c = Calls(tmp_path, 3, 3, modes={2: "files"})
    r = c.job(1)
    t = [x["token"] for x in c.runs()]
    assert r.returncode == 0 and t[0] == t[1] != t[2], r.stdout
    assert "open files" in r.stdout


def test_process_with_many_more_threads_is_renewed(tmp_path, no_job):
    c = Calls(tmp_path, 3, 3, modes={2: "threads"})
    r = c.job(1)
    t = [x["token"] for x in c.runs()]
    assert r.returncode == 0 and t[0] == t[1] != t[2], r.stdout
    assert "threads" in r.stdout


def test_a_task_that_exits_with_status_0_is_done(tmp_path, no_job):
    c = Calls(tmp_path, 2, 2, modes={1: "exit0"})
    r = c.job(1)
    assert r.returncode == 0, r.stderr
    assert c.tasks_run() == [1, 2]


def test_task_metrics_of_runner_tasks(tmp_path, no_job):
    from petakit5d.jobs.metrics import read_task_metrics

    run_dir = tmp_path / "metrics" / "run_1"
    run_dir.mkdir(parents=True)
    c = Calls(tmp_path, 6, 3, settings={"task_metrics": str(run_dir)})
    assert c.job(1).returncode == 0
    tasks, _ = read_task_metrics(str(run_dir))
    assert sorted(t["task"] for t in tasks) == [1, 2, 3, 4, 5, 6]
    assert {t["mode"] for t in tasks} == {"slurm_python"}
    assert {t["status"] for t in tasks} == {"ok"}
    assert len({t["pid"] for t in tasks}) == 1


def test_resubmitted_batch_runs_as_before(tmp_path, no_job):
    # a batch the master submitted again: its tasks in order without claims, the first error
    # ends the batch (as the job entry point runs a batch without runners)
    c = Calls(tmp_path, 6, 3, modes={2: "error"})
    c.claims.mark_resubmitted(1)
    r = c.job(1)
    assert r.returncode == 1
    assert c.tasks_run() == [1, 2]
    assert c.names() == ["r1"]
    assert r.stdout.strip().splitlines()[-1].startswith("Elapsed time is ")


@pytest.mark.parametrize("settings", [{"version": 2}, None])
def test_runner_settings_of_another_version_run_the_batch_as_before(tmp_path, no_job, settings):
    # (another version of the runner settings, or none: the batch as without runners)
    c = Calls(tmp_path, 4, 2)
    meta = json.loads(Path(c.folder, "call.json").read_text())
    if settings is None:
        del meta["task_runner"]
    else:
        meta["task_runner"].update(settings)
    Path(c.folder, "call.json").write_text(json.dumps(meta))
    assert c.job(1).returncode == 0
    assert c.tasks_run() == [1, 2] and c.names() == []


def test_outputs_given_with_a_tilde_that_exist_are_skipped(tmp_path, no_job):
    home = tmp_path / "home"
    home.mkdir()
    c = Calls(tmp_path, 2, 2)
    meta = json.loads(Path(c.folder, "call.json").read_text())
    lines = Path(c.folder, "tasks.ndjson").read_text().splitlines()
    for k, line in enumerate(lines):
        x = json.loads(line)
        x["output"] = f"~/o_{k + 1}.txt"
        (home / f"o_{k + 1}.txt").write_text("x")
        lines[k] = json.dumps(x)
    data = "".join(x + "\n" for x in lines).encode()
    meta["offsets"] = [0]
    for line in lines[:-1]:
        meta["offsets"].append(meta["offsets"][-1] + len(line.encode()) + 1)
    Path(c.folder, "tasks.ndjson").write_bytes(data)
    Path(c.folder, "call.json").write_text(json.dumps(meta))
    r = c.job(1, {"HOME": str(home)})
    assert r.returncode == 0
    assert c.runs() == []
    assert 'output "~/o_1.txt" already exists' in r.stdout


def test_claims_folder_without_symbolic_links_runs_the_batch_as_before(
    tmp_path, no_job, monkeypatch
):
    c = Calls(tmp_path, 4, 2)
    monkeypatch.setattr(Claims, "usable", lambda self: False)
    monkeypatch.syspath_prepend(str(tmp_path))
    assert runner.run(c.folder, 1) == 0
    assert c.tasks_run() == [1, 2] and c.names() == []


def test_two_jobs_share_the_tasks(tmp_path, no_job):
    # two jobs on 2 batches of 20 tasks, the second starting late (its job waited): the first
    # takes tasks of the second batch from its end, the second its first ones; each task runs
    # once, and the job that took tasks of another batch leaves its mark for the master
    c = Calls(tmp_path, 40, 20, seconds=0.1)
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    ps = []
    for b in (1, 2):
        if b == 2:
            time.sleep(2.5)
        e = dict(env, SLURM_ARRAY_JOB_ID=str(1000 + b), SLURM_ARRAY_TASK_ID=str(b))
        ps.append(
            subprocess.Popen(
                [sys.executable, "-u", "-m", "petakit5d.jobs.exec", c.folder, str(b)],
                env=e,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
        )
    outs = [p.communicate(timeout=180)[0] for p in ps]
    assert [p.returncode for p in ps] == [0, 0], outs
    done = c.tasks_run()
    assert sorted(done) == list(range(1, 41)), "a task ran twice or not at all"
    by_job = {
        j: [k for x, k in zip(c.runs(), done, strict=True) if x["job"] == j]
        for j in ("1001", "1002")
    }
    assert by_job["1001"][:20] == list(range(1, 21)) and 40 in by_job["1001"], outs
    assert by_job["1002"] and by_job["1002"][0] == 21, outs
    assert {n for n in c.names() if n.startswith("s")} == {"s1001_1"}


# --- the master side ---------------------------------------------------------------------------


class StubSlurm:
    """squeue answers from a table {key: (state, comment)}: -j <keys> with -o '%i %t', --me
    with -o '%i|%t|%k'; scancel recorded"""

    def __init__(self, jobs=None, fail=False):
        self.jobs = dict(jobs or {})
        self.fail = fail
        self.calls = []

    def run(self, argv, timeout=None):
        argv = tuple(argv)
        self.calls.append(argv)
        if argv[0] == "squeue":
            if self.fail:
                return CommandResult(argv, 1, "", "slurm_load_jobs error")
            if "-j" in argv:
                keys = argv[argv.index("-j") + 1].split(",")
                text = "".join(f"{k} {self.jobs[k][0]}\n" for k in keys if k in self.jobs)
            else:
                text = "".join(f"{k}|{s}|{c}\n" for k, (s, c) in self.jobs.items())
            return CommandResult(argv, 0, text)
        return CommandResult(argv, 0)


class Shape:
    def __init__(self, bs=1, cpus=1, slot_mode=False):
        self.task_batch_num = bs
        self.cpus_per_task = cpus
        self.slot_mode = slot_mode


def master_call(tmp_path, n, bs=1, settings=None, exist=(), missing=()):
    (tmp_path / "out").mkdir(exist_ok=True)
    ins = [str(tmp_path / f"in_{k}.txt") for k in range(n)]
    outs = [str(tmp_path / "out" / f"o_{k}.txt") for k in range(n)]
    for k, p in enumerate(ins):
        if k not in missing:
            Path(p).write_text("x")
    for k in exist:
        Path(outs[k]).write_text("x")
    recs = [task_record(Task("os:getcwd"), "python") for _ in range(n)]
    s = {"task_batch_num": bs, "master_compute": True, "max_job_num": float("inf")}
    s.update(settings or {})
    return Call(ins, outs, recs, s)


def started(call, tmp_path, slurm):
    rc = ClusterRunners(call, Shape(call.task_batch_num), slurm)
    folder = tmp_path / "pk5d_call_x"
    write_call(call, folder, rc.call_extra())
    rc.waiting([], 0.0)
    return rc


def test_runner_settings_for_the_jobs(tmp_path, no_job):
    # tasks done or without input at the start stay in their batches; the master's batch (the
    # last with a missing output) has its own key
    call = master_call(tmp_path, 6, bs=2, exist=[5], missing=[1])
    rc = ClusterRunners(call, Shape(2, cpus=4))
    extra = rc.call_extra()
    assert extra["runner"] == RUNNER
    assert extra["task_runner"] == {"version": 1, "keys": ["4", "4", "master"], "skip": [1, 5]}
    call.settings["master_compute"] = False
    call.settings["gpu_job"] = True
    rc = ClusterRunners(call, Shape(2, cpus=4))
    assert rc.call_extra()["task_runner"]["keys"] == ["4gpu"] * 3
    # bash tasks and in-job slots: none
    assert isinstance(runner.cluster_runners(call, Shape(2, slot_mode=True)), runner.NoRunners)
    call.settings["language"] = "bash"
    assert runner.cluster_runners(call, Shape(2)).call_extra() == {}


def test_master_leaves_batches_that_live_runners_hold(tmp_path, no_job):
    slurm = StubSlurm({"77_1": ("R", ""), "79_1": ("PD", ""), "80_1": ("CG", "")})
    call = master_call(tmp_path, 6, settings={"master_compute": False})
    rc = started(call, tmp_path, slurm)
    c = Claims(call.folder)
    c.claim(0, "77_1:n1:5:1")
    c.claim(1, "78_1:n1:5:1")  # its job is gone
    c.claim(3, "80_1:n1:5:1")  # completing: gone
    c.claim(4, "79_1:n2:6:1")
    rc.waiting([], 1.0)
    # first submission: left while all its tasks are held, a retry while any is
    assert rc.keep([0], retry=False)
    assert not rc.keep([0, 2], retry=False)
    assert rc.keep([0, 2], retry=True)
    assert rc.keep([4], retry=True)
    # a claim of an owner that ended is removed, and a retry marks its batch
    assert not rc.keep([1], retry=True)
    assert not rc.keep([3], retry=False)
    names = c.names()
    assert "t2" not in names and "r2" in names and "t4" not in names and "r4" not in names
    # (a job is asked about once per round)
    asked = [x[-1] for x in slurm.calls if x[0] == "squeue"]
    assert asked == ["77_1", "79_1", "78_1", "80_1"]


def test_master_leaves_batches_while_slurm_does_not_answer(tmp_path, no_job):
    slurm = StubSlurm(fail=True)
    call = master_call(tmp_path, 2, settings={"master_compute": False})
    rc = started(call, tmp_path, slurm)
    Claims(call.folder).claim(0, "77_1:n1:5:1")
    rc.waiting([], 1.0)
    assert rc.keep([0], retry=True) and rc.keep([0], retry=False)
    assert "t1" in Claims(call.folder).names()
    assert not rc.keep([1], retry=False)  # (no claims: nothing to ask)


def test_master_runs_the_tasks_it_can_claim(tmp_path, no_job):
    slurm = StubSlurm({"77_1": ("R", "")})
    call = master_call(tmp_path, 4, bs=4, exist=[2])
    rc = started(call, tmp_path, slurm)
    c = Claims(call.folder)
    c.claim(0, "77_1:n1:5:1")
    c.claim(1, "78_1:n1:5:1")
    assert rc.master_tasks([0, 1, 2, 3]) == [1, 3]
    assert claims.parse_owner(c.owner_of(1))[3] == os.getpid()


def test_job_limit_counts_jobs_that_run_other_batches_tasks(tmp_path, no_job):
    call = master_call(tmp_path, 8, settings={"master_compute": False, "max_job_num": 4})
    rc = ClusterRunners(call, Shape(1), StubSlurm())
    folder = tmp_path / "pk5d_call_x"
    write_call(call, folder, rc.call_extra())
    base = folder.name
    rc._slurm = StubSlurm({f"100{k}_{k}": ("R", f"{base}:{k}") for k in (1, 2, 3)})
    rc.waiting([], 0.0)
    # without a job that takes other batches' tasks the loop's count stands
    assert not rc.keep([4], retry=False) and not rc.keep([5], retry=False)
    Claims(folder).mark_stealer("1001_1", "1001_1:n1:5:1")
    rc.waiting([], 1.0)
    assert not rc.keep([6], retry=False)
    assert rc.keep([7], retry=False), "a fifth job of the call in the queue"
    rc._slurm.fail = True
    rc.waiting([], 2.0)
    assert rc.keep([7], retry=False)


def test_pending_jobs_of_batches_that_runners_did_are_cancelled(tmp_path, no_job):
    call = master_call(tmp_path, 6, bs=2, settings={"master_compute": False}, exist=[2])
    rc = ClusterRunners(call, Shape(2), StubSlurm())
    folder = tmp_path / "pk5d_call_x"
    write_call(call, folder, rc.call_extra())
    base = folder.name
    jobs = {
        "1000_1": ("R", f"{base}:1"),
        "1001_2": ("PD", f"{base}:2"),
        "1002_3": ("PD", f"{base}:3"),
        "999_2": ("PD", "pk5d_call_other:2"),
    }
    rc._slurm = StubSlurm(jobs)
    c = Claims(folder)
    c.mark_stealer("1000_1", "1000_1:n1:5:1")
    c.claim(3, "1000_1:n1:5:1")  # batch 2: task 3 done, task 4 taken by a live job
    c.claim(4, "1000_1:n1:5:1")  # batch 3: task 6 not taken
    rc.waiting([], 1.0)
    cancels = [x for x in rc._slurm.calls if x[0] == "scancel"]
    assert cancels == [("scancel", "--state=PENDING", "1001_2")]
    assert rc.cancelled == ["1001_2"]


# --- with the cluster loop (fake Slurm) ----------------------------------------------------------


class VirtualClock:
    def __init__(self):
        self.now = 1000.0
        self.hooks = []

    def clock(self):
        return self.now

    def sleep(self, s):
        self.now += s + 0.001
        while self.hooks and self.hooks[0][0] <= self.now:
            self.hooks.pop(0)[1]()


@pytest.fixture
def loop_with_runners(monkeypatch, no_job):
    """the cluster loop on a virtual clock, with this module's runners as its plugin (the loop
    looks for them as petakit5d.jobs.runners until its name is changed)"""
    import importlib

    # (the module run_tasks will use: another test module's fixture may have dropped the one
    # imported when this module was)
    loop = importlib.import_module("petakit5d.jobs.cluster")
    vc = VirtualClock()
    monkeypatch.setattr(loop, "_clock", vc.clock)
    monkeypatch.setattr(loop, "_sleep", vc.sleep)
    monkeypatch.setitem(sys.modules, "petakit5d.jobs.runners", runner)
    return vc


def fake_tasks(tmp_path, n):
    (tmp_path / "out").mkdir()
    (tmp_path / "logs").mkdir()
    ins = [str(tmp_path / f"in_{k}.txt") for k in range(1, n + 1)]
    outs = [str(tmp_path / "out" / f"out_{k}.txt") for k in range(1, n + 1)]
    for p in ins:
        Path(p).write_text("x")
    return ins, outs, [Task("os:getcwd") for _ in outs]


def run_cluster_call(fake, ins, outs, tasks, tmp_path, **params):
    from petakit5d.jobs.framework import run_tasks

    old = os.environ["PATH"]
    os.environ["PATH"] = fake.path_with(old)
    try:
        return run_tasks(
            ins, outs, tasks, job_log_dir=str(tmp_path / "logs"), slurm_param="-p queue", **params
        )
    finally:
        os.environ["PATH"] = old


def call_json_of(tmp_path):
    (p,) = (tmp_path / "logs").glob("pk5d_call_*/call.json")
    return json.loads(p.read_text())


def test_cluster_jobs_get_the_runner(tmp_path, loop_with_runners):
    # (the fake jobs never run: the call ends with its tasks not done, and its folder stays)
    import shutil

    fake = FakeSlurm(tmp_path / "fake").install()
    ins, outs, tasks = fake_tasks(tmp_path, 2)
    run_cluster_call(fake, ins, outs, tasks, tmp_path, master_compute=False, max_trial_num=1)
    meta = call_json_of(tmp_path)
    assert meta["runner"] == RUNNER
    assert meta["task_runner"]["keys"] == ["1", "1"]
    wraps = [x.split("--wrap=", 1)[1] for x in fake.calls("sbatch")]
    assert len(wraps) == 2 and all(" -u -m petakit5d.jobs.exec " in w for w in wraps)
    # taskRunner false: no runner
    for p in (tmp_path / "logs").glob("pk5d_call_*"):
        shutil.rmtree(p)
    fake2 = FakeSlurm(tmp_path / "fake2").install()
    run_cluster_call(
        fake2, ins, outs, tasks, tmp_path, master_compute=False, max_trial_num=1, task_runner=False
    )
    assert "runner" not in call_json_of(tmp_path)


def test_master_leaves_a_task_to_a_running_runner_until_its_job_is_gone(
    tmp_path, loop_with_runners
):
    # MATLAB's masterLeavesTasksToRunningRunnersAndResubmitsStaleOnes: job 1000's runner takes
    # task 3 before batch 3 is submitted (maxJobNum 2); batch 3 is submitted once, after 1000
    # left the queue without task 3's output (its claim is then removed)
    fake = FakeSlurm(tmp_path / "fake").install()
    ins, outs, tasks = fake_tasks(tmp_path, 3)
    fake.set_job_states("1000_1", ["R"] * 6 + ["done"], [outs[0]])
    fake.set_job_states("1001_2", ["R", "done"], [outs[1]])
    fake.set_job_states("1002_3", ["R"], [outs[2]])

    def claim_task_3():
        for d in (tmp_path / "logs").glob("pk5d_call_*"):
            Claims(d).claim(2, "1000_1:n1:5:1")

    loop_with_runners.hooks.append((1000.5, claim_task_3))
    flags = run_cluster_call(
        fake,
        ins,
        outs,
        tasks,
        tmp_path,
        master_compute=False,
        max_trial_num=2,
        max_job_num=2,
        query_interval=1,
    )
    assert flags == [True] * 3
    hist = fake.history()
    submit3 = [k for k, x in enumerate(hist) if x.startswith("sbatch") and "--array=3 " in x]
    assert len(submit3) == 1, "batch 3 should be submitted once"
    running = [k for k, x in enumerate(hist) if x.startswith("squeue") and "1000_1" in x]
    assert len(running) >= 6
    assert submit3[0] > running[5], "batch 3 was submitted while the job that took it ran"
