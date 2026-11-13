"""Pools of single jobs (jobs.pools): the single jobs a process submits with the same function
and resources share a pool, so a job, after its own task, runs the unstarted tasks of the
pool's other jobs that Slurm keeps pending, and cancels those jobs; the callers' status checks
count such a task as running, and the master takes a task back before it submits it again or
runs it (MATLAB's SingleJobRunnerTest).

A kind starts with a first task (start_kind), so that tasks 1, 2, ... are the pool's tasks 1, 2,
... with jobs 1000, 1001, ... The jobs' commands, from the fake Slurm of the MATLAB suite, run
in a shell as Slurm would run them, with Slurm's job variables. MATLAB's Slurm command histories
of the same tests (oracle tree PetaKit5D_7ed28c0) are compared with ours: the commands are the
same, apart from the job's command line (the wrap) and the differences each test names."""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

import petakit5d
from petakit5d.jobs import claims, config, pools, single
from petakit5d.jobs.single import SingleJobWarning, check_single_job_status, submit_single_job
from petakit5d.jobs.tasks import Task, task_record
from petakit5d.testing import slurm_fake
from petakit5d.testing.slurm_fake import SLURM_JOB_VARIABLES, FakeSlurm, hide_slurm

fake_slurm = slurm_fake.fake_slurm  # (the fixture)

pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="Task runners run on Linux only."
)

SRC_ROOT = str(Path(petakit5d.__file__).resolve().parents[1])
ABC = "-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M"
WAIT = "squeue -r -h -u <USER> -t PENDING -o %i %r"
OTHER = "99_1:othernode:1:1"
TASKS = '''
import json
import os
import shutil
import time
import uuid

TOKEN = uuid.uuid4().hex  # one per interpreter: a new process has another


def work(out, log, mode="ok", seconds=0.0):
    """appends what it sees to log, then writes out (by a temporary name and a rename)"""
    rec = {"out": out, "pid": os.getpid(), "token": TOKEN, "mode": mode,
           "job": os.environ.get("SLURM_ARRAY_JOB_ID", "")}
    with open(log, "a") as fh:
        fh.write(json.dumps(rec) + "\\n")
    if mode == "error":
        raise RuntimeError("on purpose")
    if seconds:
        time.sleep(seconds)
    tmp = f"{out}_{uuid.uuid4().hex}"
    with open(tmp, "w") as fh:
        fh.write(mode)
    os.replace(tmp, out)


def other(out, log, mode="ok", seconds=0.0):
    work(out, log, mode, seconds)


def copy(src, out, log):
    work(out, log)
    tmp = f"{out}_{uuid.uuid4().hex}"
    shutil.copyfile(src, tmp)
    os.replace(tmp, out)
'''

# --- MATLAB's histories (oracle tree PetaKit5D_7ed28c0, the recorded fake's history.log) -------
#
# Compact: "sub <array> [<cpus>]" an sbatch of these tests (<T1>/logs, {ABC}, -t 24:0:00, the
# wrap set aside), "wait" a job's squeue of the user's pending jobs (task_runner.sh's
# waiting_tasks), "take <key>" its scancel of a pending job whose task it took. A test with
# several fake Slurms has a history per fake, in the order of their teardown (the last first).

MATLAB_HISTORIES: dict[str, list[list[str]]] = {
    "compiledJobsRunTheUnstartedTasksOfTheirPool": [
        ["sub 1", "sub 2", "sub 3", "sub 4", "wait", "take 1001_2", "wait", "take 1002_3"]
        + ["wait", "take 1003_4"],
    ],
    "firstTaskOfAKindRunsAsBefore": [["sinfo ", "sub 1", "squeue -j 1000_1", "sub 2"]],
    "runnersTakeOnlyTasksWhoseJobsWaitInTheQueue": [
        ["sub 1", "sub 2", "sub 3", "sub 4", "sub 5", "sub 6", "wait", "take 1002_3"]
        + ["wait", "take 1005_6", "wait"],
    ],
    "jobWhoseTaskARunnerHasEndsAtOnce": [["sub 1", "sub 2", "sub 3"]],
    "slotStartsWithItsOwnTask": [["sub 1", "sub 2", "sub 3"]],
    "matlabJobGoesOnAfterAProcessEnded": [
        ["sub 1", "sub 2", "sub 3", "wait", "take 1001_2", "wait", "take 1002_3"],
    ],
    "jobWithNothingToTakeRunsItsTaskAsBefore": [["sub 1"]],
    "tasksSubmittedWhileAJobRunsAreTaken": [
        ["sub 1", "sub 2", "sub 3", "wait", "take 1001_2", "wait", "take 1002_3"],
    ],
    "poolsKeepKindsApart": [
        ["sinfo ", "sub 1", "sub 2 2", "sub 3", "sub 4", "sub 5 2", "sub 6", "sub 7 2", "sub 8"],
    ],
    "offSwitchOldBuildsAndRelativeLogFoldersSubmitAsBefore": [
        ["sinfo ", "sub 1", "squeue -j 555_2", "sub 2"],
        ["sinfo ", "sub 1", "sub 2"],
        [
            "sbatch --array=1 -o job_%A_%a.out -e <T1>/logs/job_%A_%a.err --cpus-per-task=1 "
            "{ABC} -t 24:0:00 --wrap=<W>",
            "sbatch --array=2 -o job_%A_%a.out -e <T1>/logs/job_%A_%a.err --cpus-per-task=1 "
            "{ABC} -t 24:0:00 --wrap=<W>",
        ],
        ["sinfo ", "sub 1", "sub 2"],
        ["sinfo ", "sub 1", "sub 2"],
        ["sinfo ", "sub 1", "sub 2"],
    ],
    "taskTakenBetweenSubmissionTriesIsSubmittedAsBefore": [
        ["squeue -j 999_9", "squeue -j 999_9", "sub 1"],
    ],
    "poolThatCannotBeMadeLeavesJobsAsBefore": [["sub 1", "sub 2"]],
    "claimsFolderThatFailsLeavesTasksToTheMaster": [["sub 1", "squeue -j 1000_1", "sub 1"]],
    "matlabRunnerTakesTaskStringsFromThePool": [
        ["sub 1", "sub 2", "sub 3", "wait", "take 1001_2", "wait", "take 1002_3"],
    ],
    "matlabPoolTasksRecordTheirMetrics": [["sub 7 2", "sub 8 2", "wait", "take 1001_8"]],
    "liveRunnerKeepsTheTaskOfAnEndedJob": [
        ["sub 1", "sub 2", "squeue -j 1001_2", "squeue -j 1000_1", "squeue -j 1001_2"]
        + ["squeue -j 1000_1", "squeue -j 1001_2", "squeue -j 1000_1", "squeue -j 1000_1"]
        + ["sub 2"],
    ],
    "runnerKeepsItsTaskWhileSqueueFails": [
        ["sub 1", "sub 2", "squeue -j 1000_1", "squeue -j 1000_1"],
    ],
    "pendingJobGetsBackTheTaskOfARunnerThatEnded": [
        ["sub 1", "sub 2", "squeue -j 1001_2", "squeue -j 1000_1"],
    ],
    "pendingJobOfATaskARunnerFailedIsReplaced": [
        ["sub 1", "sub 2", "sub 3", "squeue -j 1001_2", "scancel 1001_2", "sub 2"]
        + ["squeue -j 1002_3", "scancel 1002_3"],
    ],
    "failedTaskIsSubmittedAgainAsBefore": [["sub 1", "squeue -j 1000_1", "sub 1"]],
    "singleJobStatusCountsTaskRunners": [
        ["sub 1", "sub 2", "sub 3", "squeue -j 1001_2", "squeue -j 1000_1", "squeue -j 1001_2"]
        + ["squeue -j 1000_1", "squeue -j 1002_3", "squeue -j 1000_1", "squeue -j 77_4"]
        + ["squeue -j 77_4"],
    ],
    "masterTakesTheLastFileFromThePool": [
        ["sub 1", "sub 2", "squeue -j 1001_2", "squeue -j 1000_1", "squeue -j 1001_2"]
        + ["scancel 1001_2"],
    ],
    "failedSubmissionLeavesTheTaskToTheMaster": [["sub 1", "sub 1", "sub 1"]],
    "endedPoolGoesWithItsLastJob": [
        ["sub 1", "sub 2", "wait", "squeue -h -r -o %i -j 1000_1,1001_2", "sub 3"],
    ],
    "endedPoolWithoutJobsGoesAtOnce": [
        ["sub 1", "sub 2", "sub 3", "wait", "take 1001_2", "wait", "take 1002_3"]
        + ["squeue -h -r -o %i -j 1000_1,1001_2,1002_3"],
    ],
    "pendingJobKeepsAnEndedPool": [
        ["sub 1", "sub 2", "wait", "take 1001_2", "squeue -h -r -o %i -j 1000_1,1001_2"],
    ],
    "poolOfAnEnclosingCallStays": [["sub 1", "sub 2"]],
}


def expand(lines: list[str]) -> list[str]:
    """the commands of a compact history"""
    out = []
    for x in lines:
        f = x.split()
        if f[0] == "sub":
            cpus = f[2] if len(f) > 2 else "1"
            out.append(
                f"sbatch --array={f[1]} -o <T1>/logs/job_%A_%a.out -e <T1>/logs/job_%A_%a.err "
                f"--cpus-per-task={cpus} {ABC} -t 24:0:00 --wrap=<W>"
            )
        elif x == "wait":
            out.append(WAIT)
        elif f[0] == "take":
            out.append(f"scancel --state=PENDING {f[1]}")
        else:
            out.append(x.replace("{ABC}", ABC))
    return out


def matlab_history(name: str, k: int = 0) -> list[str]:
    return expand(MATLAB_HISTORIES[name][k])


@pytest.fixture(autouse=True)
def fresh_session(tmp_path_factory, monkeypatch):
    """a new session (MATLAB: clear task_runner_run generic_single_job_submit_wrapper
    single_job_pools), Slurm hidden unless a test asks for the fake, no job variables and none
    of the runner's settings"""
    hide_slurm(monkeypatch, tmp_path_factory.mktemp("no_slurm"))
    for name in (*SLURM_JOB_VARIABLES, "SLURM_CPUS_PER_TASK", "PETAKIT5D_RUNNER_STATE"):
        monkeypatch.delenv(name, raising=False)
    for name in ("PETAKIT5D_RUNNER_MEM_GB", "PETAKIT5D_RUNNER_MINUTES"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("USER", os.environ.get("USER") or "tester")
    single.reset()
    config.clear_cache()
    yield
    single.reset()
    config.clear_cache()


def quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


class Setup:
    """a job log folder, task outputs o_1 ... o_10 and the tasks module (MATLAB's
    compiledSetup and its helpers)"""

    def __init__(self, root: Path, fake: FakeSlurm) -> None:
        self.root = root
        self.fake = fake
        self.logs = root / "logs"
        self.logs.mkdir(parents=True, exist_ok=True)
        (root / "out").mkdir(exist_ok=True)
        (root / "ptasks.py").write_text(TASKS)
        self.log = root / "runs.log"
        self.outs = [""] + [str(root / "out" / f"o_{k}.txt") for k in range(1, 11)]

    def task(self, k: int, mode: str = "ok", func: str = "ptasks:work", seconds: float = 0):
        return Task(func, [self.outs[k], str(self.log), mode, seconds], path=[str(self.root)])

    def key(self, k: int, mode: str = "ok") -> str:
        return single._task_key(task_record(self.task(k, mode), "python"))

    def submit(
        self,
        k: int,
        mode: str = "ok",
        job_id: int = -1,
        func: str = "ptasks:work",
        seconds: float = 0,
        task: Task | None = None,
        **kw,
    ):
        """the wrapper for task k (array task k) and its job job_id, one sbatch try"""
        kw.setdefault("job_log_dir", str(self.logs))
        kw.setdefault("job_log_fname", f"{self.logs}/job_%A_%a.out")
        kw.setdefault("job_error_fname", f"{self.logs}/job_%A_%a.err")
        kw.setdefault("max_trial_num", 1)
        kw.setdefault("slurm_param", ABC)
        task = self.task(k, mode, func, seconds) if task is None else task
        return quiet(submit_single_job, task, job_id, k, **kw)

    def start_kind(self, **kw) -> None:
        """the kind's first task (task 10): it runs as before, without a pool; then the fake
        forgets it, so that tasks 1, 2, ... are the pool's, with jobs 1000, 1001, ..."""
        assert self.submit(10, **kw)[2] == 1
        assert not self.pools(), "The first task of a kind made a pool."
        self.fake.clear_calls()

    def pools(self) -> list[Path]:
        return sorted(self.logs.glob(f"{single.POOL_PREFIX}*"))

    def pool(self) -> Path:
        ps = self.pools()
        assert len(ps) == 1, ps
        return ps[0]

    def job_env(self, job_id: int, array_id: int, env=None) -> dict:
        e = dict(os.environ, PYTHONPATH=SRC_ROOT, SLURM_JOB_ID=str(job_id))
        e.update(SLURM_ARRAY_JOB_ID=str(job_id), SLURM_ARRAY_TASK_ID=str(array_id))
        e.update(SLURMD_NODENAME="testnode")
        e.pop("SLURM_JOB_END_TIME", None)
        e.update(env or {})
        return e

    def run_job(self, k: int, env=None, may_fail: bool = False, background: bool = False):
        """the k-th submitted job's command (array job 999 + k) in a shell, as Slurm runs it"""
        cmd = self.fake.calls("sbatch")[k - 1]
        array_id = int(re.search(r"--array=(\d+)", cmd).group(1))
        argv = ["sh", "-c", cmd.split("--wrap=", 1)[1]]
        e = self.job_env(999 + k, array_id, env)
        if background:
            return subprocess.Popen(
                argv, env=e, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
            )
        r = subprocess.run(argv, capture_output=True, text=True, env=e, timeout=180)
        if not may_fail:
            assert r.returncode == 0, r.stdout + r.stderr
        return r

    def runs(self) -> list[dict]:
        if not self.log.exists():
            return []
        return [json.loads(x) for x in self.log.read_text().splitlines() if x.strip()]

    def tasks_run(self) -> list[int]:
        return [self.outs.index(r["out"]) for r in self.runs()]

    def claim_of(self, i: int, pool: Path | None = None) -> tuple[str, bool, bool]:
        """the owner of task i's claim, and its taken-back (r) and failed marks"""
        c = (self.pool() if pool is None else pool) / "claims"
        try:
            owner = os.readlink(c / f"t{i}")
        except OSError:
            owner = ""
        return owner, (c / f"r{i}").exists(), os.path.islink(c / f"f{i}")

    def claim(self, i: int, owner: str, failed: bool = False) -> None:
        os.symlink(owner, self.pool() / "claims" / f"t{i}")
        if failed:
            os.symlink(owner, self.pool() / "claims" / f"f{i}")

    def history(self) -> list[str]:
        """the fake's commands with the test's folder as <T1>, the user as <USER> and the wrap
        of submissions as <W>"""
        out = []
        for line in self.fake.history():
            line = line.replace(str(self.root), "<T1>")
            line = re.sub(r"(?<= -u )\S+", "<USER>", line)
            if line.startswith("sbatch ") and "--test-only" not in line:
                line = re.sub(r" --wrap=.*$", " --wrap=<W>", line)
            out.append(line)
        return out


@pytest.fixture
def c(tmp_path, fake_slurm) -> Setup:
    return Setup(tmp_path, fake_slurm)


def wrap(cmd: str) -> str:
    return cmd.split("--wrap=", 1)[1]


def in_pool(cmd: str) -> bool:
    return f"/{single.POOL_PREFIX}" in wrap(cmd)


def as_before(cmd: str) -> bool:
    return f"/{single.TASK_FOLDER_PREFIX}" in wrap(cmd) and not in_pool(cmd)


# --- SingleJobRunnerTest ----------------------------------------------------------------------


def test_jobs_run_the_unstarted_tasks_of_their_pool(c):
    # (MATLAB: compiledJobsRunTheUnstartedTasksOfTheirPool) four single jobs of one kind: the
    # first one to start runs its task, then the others' in their order in its process, and
    # cancels their jobs while they are pending; if they start, they end without a task. The
    # jobs have the resources of before
    c.start_kind()
    for k in range(1, 5):
        assert c.submit(k)[2] == 1
    for k in range(2, 5):
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    cmds = c.fake.calls("sbatch")
    assert len(cmds) == 4
    for k in range(1, 5):
        assert re.match(
            rf"^--array={k} -o \S+ -e \S+ --cpus-per-task=1 -p cpu_queue .*-t 24:0:00 "
            r"--wrap=date; echo Python command: .* -u -m petakit5d\.jobs\.exec "
            rf"\S+/{single.POOL_PREFIX}\S+ {k}$",
            cmds[k - 1],
        ), cmds[k - 1]
    assert len(c.pools()) == 1, "The jobs are not in one pool."
    r = c.run_job(1)
    runs = c.runs()
    assert c.tasks_run() == [1, 2, 3, 4], r.stdout
    assert len({x["token"] for x in runs}) == 1, "The tasks did not run in the job's process."
    assert {x["job"] for x in runs} == {"1000"}
    assert c.fake.calls("scancel") == [
        "--state=PENDING 1001_2",
        "--state=PENDING 1002_3",
        "--state=PENDING 1003_4",
    ]
    for k in range(2, 5):
        assert "no task left" in c.run_job(k).stdout
    assert len(c.runs()) == 4, "A task ran more than once."
    assert all(os.path.isfile(c.outs[k]) for k in range(1, 5))
    assert c.history() == matlab_history("compiledJobsRunTheUnstartedTasksOfTheirPool")


FIRST_KIND = """
import contextlib, io, json, os, sys
sys.path.insert(0, sys.argv[1])
from petakit5d.jobs import single
from petakit5d.jobs.tasks import Task

root = sys.argv[2]
logs = root + "/logs"
task = lambda k: Task("ptasks:work", [f"{root}/out/o_{k}.txt", f"{root}/runs.log", "ok", 0],
                      path=[root])
kw = dict(job_log_dir=logs, job_log_fname=logs + "/job_%A_%a.out",
          job_error_fname=logs + "/job_%A_%a.err", max_trial_num=1, slurm_param=sys.argv[3])
with contextlib.redirect_stdout(io.StringIO()):
    single.submit_single_job(task(1), -1, 1, **kw)
    pools_after_first = sorted(x for x in os.listdir(logs) if x.startswith("task_runner_single_"))
    status = single.check_single_job_status(1000, 1)
    loaded = "petakit5d.jobs.pools" in sys.modules
    single.submit_single_job(task(2), -1, 2, **kw)
pool = [x for x in os.listdir(logs) if x.startswith("task_runner_single_")]
files = {n: open(f"{logs}/{pool[0]}/{n}").read() for n in ("runnable.txt", "jobs.txt")}
print(json.dumps({"pools": pools_after_first, "status": status, "loaded": loaded, **files}))
"""


def test_first_task_of_a_kind_runs_as_before(c):
    # a session's first task of a kind runs as before: a session with one task pays nothing
    # (no pool; neither its submission nor its status checks import the pools); the second
    # task makes the pool, with itself as the pool's task 1
    script = c.root / "first.py"
    script.write_text(FIRST_KIND)
    r = subprocess.run(
        [sys.executable, str(script), SRC_ROOT, str(c.root), ABC],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, r.stderr
    out = json.loads(r.stdout.strip().splitlines()[-1])
    assert out == {
        "pools": [],
        "status": -1,
        "loaded": False,
        "runnable.txt": "1\n",
        "jobs.txt": "1 1001_2\n",
    }
    cmds = c.fake.calls("sbatch")
    assert len(cmds) == 2
    assert as_before(cmds[0]), cmds[0]
    assert re.search(rf"/{single.POOL_PREFIX}\S+ 1$", cmds[1]), cmds[1]
    # (the process left its folders when it ended: Slurm has none of their jobs, so they went)
    assert not c.pools()
    exit_query = "squeue -h -r -o %i -j 1000_1,1001_2"
    assert c.history() == [*matlab_history("firstTaskOfAKindRunsAsBefore"), exit_query]


def test_runners_take_only_tasks_whose_jobs_wait_in_the_queue(c):
    # a job's runner takes another job's task only while Slurm keeps that job in the queue
    # (pending with a reason: resources, priority, a limit), not a job its scheduler hasn't
    # looked at yet (reason None: it may start at once), a job that runs, or a task whose job
    # is not recorded yet
    c.start_kind()
    for k in range(1, 7):
        c.submit(k)
    # 2: just submitted; 3: pending (resources); 4: its job runs; 5: pending, its job not
    # recorded; 6: pending (priority)
    states = ["PD:None", "PD", "R", "PD:Priority", "PD:Priority"]
    for k, st in zip(range(2, 7), states, strict=True):
        c.fake.set_job_states(f"{999 + k}_{k}", [st])
    jobs = (c.pool() / "jobs.txt").read_text().splitlines()
    (c.pool() / "jobs.txt").write_text("".join(x + "\n" for x in jobs if not x.startswith("5 ")))
    r = c.run_job(1)
    assert c.tasks_run() == [1, 3, 6], r.stdout
    assert c.fake.calls("scancel") == ["--state=PENDING 1002_3", "--state=PENDING 1005_6"]
    assert c.history() == matlab_history("runnersTakeOnlyTasksWhoseJobsWaitInTheQueue")


def test_job_whose_task_a_runner_has_ends_at_once(c):
    # a job whose own task another job's runner took ends without running a task, and takes
    # no other task (those have jobs of their own)
    c.start_kind()
    for k in range(1, 4):
        c.submit(k)
    c.claim(2, OTHER)
    assert "no task left" in c.run_job(2).stdout
    assert c.runs() == []
    assert c.history() == matlab_history("jobWhoseTaskARunnerHasEndsAtOnce")


def test_job_starts_with_its_own_task(c):
    # (MATLAB: slotStartsWithItsOwnTask) the first task of a single job is its own: if another
    # process took it, the job takes no other task, also while others wait in the queue
    c.start_kind()
    for k in range(1, 4):
        c.submit(k)
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    c.claim(2, OTHER)
    assert "no task left" in c.run_job(2).stdout
    assert c.runs() == [], "The job took a task of another job first."
    assert c.fake.call_count("squeue") == 0
    assert c.history() == matlab_history("slotStartsWithItsOwnTask")


def test_job_goes_on_after_a_process_ended(c):
    # (MATLAB: matlabJobGoesOnAfterAProcessEnded) a job whose process ends early (its memory,
    # time or open-file guard) goes on in a new process with the pool's next tasks
    c.start_kind()
    for k in range(1, 4):
        c.submit(k)
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    r = c.run_job(1, {"PETAKIT5D_RUNNER_MINUTES": "0"})
    assert c.tasks_run() == [1, 2, 3], r.stdout
    runs = c.runs()
    assert len({x["token"] for x in runs}) == 3, "The tasks did not run in new processes."
    assert len({x["pid"] for x in runs}) == 1, "The job did not go on in its own process."
    assert r.stdout.count("a new process runs task") == 2
    assert c.history() == matlab_history("matlabJobGoesOnAfterAProcessEnded")


def test_job_with_nothing_to_take_runs_its_task_as_before(c):
    # with no unstarted task in its pool, a job runs its task and asks Slurm nothing
    c.start_kind()
    c.submit(1)
    c.run_job(1)
    assert c.tasks_run() == [1]
    assert c.fake.call_count("squeue") == 0
    assert c.history() == matlab_history("jobWithNothingToTakeRunsItsTaskAsBefore")


def test_tasks_submitted_while_a_job_runs_are_taken(c):
    # (streaming) a task added while a job runs is taken by it, after the tasks submitted
    # before it
    c.start_kind()
    c.submit(1, seconds=3)
    c.submit(2)
    c.fake.set_job_states("1001_2", ["PD"])
    p = c.run_job(1, background=True)
    t0 = time.monotonic()
    while not c.runs() and time.monotonic() - t0 < 60:
        time.sleep(0.05)
    c.submit(3)
    c.fake.set_job_states("1002_3", ["PD"])
    out, _ = p.communicate(timeout=120)
    assert p.returncode == 0, out
    assert c.tasks_run() == [1, 2, 3], out
    assert len({x["token"] for x in c.runs()}) == 1
    assert c.history() == matlab_history("tasksSubmittedWhileAJobRunsAreTaken")


def test_pools_keep_kinds_apart(c):
    # a pool has the tasks of one function and its resources; another CPU count or function
    # makes another (each with its second task)
    kinds = [{}, {"cpus_per_task": 2}, {"func": "ptasks:other"}, {}, {"cpus_per_task": 2}, {}]
    kinds += [{"cpus_per_task": 2}, {"func": "ptasks:other"}]
    for k, kw in enumerate(kinds, 1):
        c.submit(k, **kw)
    pools_ = c.pools()
    assert len(pools_) == 3
    n = sorted(len((p / "runnable.txt").read_text().split()) for p in pools_)
    assert n == [1, 2, 2]
    assert c.history() == matlab_history("poolsKeepKindsApart")


def test_pools_keep_job_log_folders_apart(c):
    other = c.root / "logs2"
    other.mkdir()
    kw = {
        "job_log_dir": str(other),
        "job_log_fname": f"{other}/job_%A_%a.out",
        "job_error_fname": f"{other}/job_%A_%a.err",
    }
    for k in (1, 2):
        c.submit(k)
        c.submit(k + 4, **kw)
    assert len(c.pools()) == 1
    assert len(list(other.glob(f"{single.POOL_PREFIX}*"))) == 1


def test_off_switch_and_relative_log_folders_submit_as_before(c, monkeypatch):
    # (MATLAB: offSwitchOldBuildsAndRelativeLogFoldersSubmitAsBefore) task_runner False (as an
    # option or a config key; the caller's False wins over the config's True), a job log folder
    # that is not an absolute path, or a job of before: the second task of a kind runs as
    # before too, and there is no pool
    cases = [
        ({}, {"task_runner": False}, 5),
        ({"taskRunner": False}, {}, 3),
        ({}, {"job_log_fname": "job_%A_%a.out", "job_log_dir": "no_job_log_folder"}, 2),
        ({"taskRunner": True}, {"task_runner": False}, 1),
        ({}, {"job_id": 555}, 0),
    ]
    for n, (cfg, args, h) in enumerate(cases):
        root = c.root / f"case{n}"
        s = Setup(root, c.fake)
        s.fake.clear_calls()
        monkeypatch.chdir(root)
        args = dict(args)
        job_id = args.pop("job_id", -1)
        if cfg:
            (root / "config.json").write_text(json.dumps(cfg))
            args["config_file"] = str(root / "config.json")
        s.submit(1, **args)
        s.submit(2, job_id=job_id, **args)
        cmds = s.fake.calls("sbatch")
        assert len(cmds) == 2, n
        assert all(as_before(x) for x in cmds), (n, cmds)
        assert not s.pools(), f"Case {n} made a pool."
        assert s.history() == matlab_history(
            "offSwitchOldBuildsAndRelativeLogFoldersSubmitAsBefore", h
        ), n
        single.reset()
        config.clear_cache()


def test_task_taken_between_submission_tries_is_submitted_as_before(c, monkeypatch):
    # sbatch fails once, and meanwhile a runner (whose job ended since) claimed the new task:
    # the next try submits it as before (its own job could not run it), and the master has it
    c.start_kind()
    stubs = c.root / "stub"
    stubs.mkdir()
    (stubs / "sbatch").write_text(
        "#!/bin/sh\n"
        f'if [ ! -f "{stubs}/failed" ]; then\n'
        f'    touch "{stubs}/failed"\n'
        f'    for d in "{c.logs}"/{single.POOL_PREFIX}*/claims; do '
        'ln -s "999_9:othernode:5:1" "$d/t1"; done\n'
        "    exit 1\n"
        "fi\n"
        f'exec "{c.fake.folder}/sbatch" "$@"\n'
    )
    (stubs / "sbatch").chmod(0o755)
    monkeypatch.setenv("PATH", f"{stubs}{os.pathsep}{os.environ['PATH']}")
    assert c.submit(1, max_trial_num=2)[2] == 1
    cmds = c.fake.calls("sbatch")
    assert len(cmds) == 1
    assert as_before(cmds[0]), cmds[0]
    assert c.claim_of(1)[0].endswith(":0"), "The master did not take the task."
    assert c.history() == matlab_history("taskTakenBetweenSubmissionTriesIsSubmittedAsBefore")


def test_pool_that_cannot_be_made_leaves_jobs_as_before(c, monkeypatch):
    # a job log folder where no pool can be made (its claims take no symbolic links): a
    # warning, and the jobs are submitted as before
    c.start_kind()
    monkeypatch.setattr(claims.Claims, "usable", lambda self: False)
    with pytest.warns(SingleJobWarning) as rec:
        c.submit(1)
    assert [w.message.identifier for w in rec] == ["task_runner_run:single"]
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        c.submit(2)
    cmds = c.fake.calls("sbatch")
    assert len(cmds) == 2
    assert all(as_before(x) for x in cmds), "A job was submitted with task runners."
    assert not c.pools(), "The pool that could not be made is left."
    assert c.history() == matlab_history("poolThatCannotBeMadeLeavesJobsAsBefore")


def test_read_only_job_log_folder_leaves_jobs_as_before_with_their_tasks_elsewhere(c):
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        pytest.skip("root writes into read-only folders")
    c.start_kind(tmp_dir=str(c.root / "tasks"))
    c.logs.chmod(0o555)
    try:
        with pytest.warns(SingleJobWarning, match="without task runners"):
            assert c.submit(1, tmp_dir=str(c.root / "tasks"))[2] == 1
        cmd = c.fake.calls("sbatch")[0]
        assert as_before(cmd) and str(c.root / "tasks") in wrap(cmd)
    finally:
        c.logs.chmod(0o755)


def test_claims_folder_that_fails_leaves_tasks_to_the_master(c):
    # a claims folder that takes no new claims (here read-only): the master can't claim the
    # task of an ended job either, but as no process has it, it is submitted again as before
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        pytest.skip("root writes into read-only folders")
    c.start_kind()
    c.submit(1)
    claims_dir = c.pool() / "claims"
    claims_dir.chmod(0o555)
    try:
        job_id, _, submitted = c.submit(1, job_id=1000)
    finally:
        claims_dir.chmod(0o755)
    assert (job_id, submitted) == (1001, 1)
    assert as_before(c.fake.calls("sbatch")[-1])
    assert c.history() == matlab_history("claimsFolderThatFailsLeavesTasksToTheMaster")


def test_task_with_line_breaks_runs_from_its_pool(c):
    # (MATLAB: taskWithALineBreakIsSubmittedAsBefore) a pool keeps a task's line in a line of
    # its task file: a Python task record is one JSON line whatever its arguments hold, so a
    # task with line breaks joins the pool and runs from it
    c.start_kind()
    c.submit(1, mode="o\nk")
    cmd = c.fake.calls("sbatch")[0]
    assert in_pool(cmd) and "\n" not in cmd
    c.run_job(1)
    assert Path(c.outs[1]).read_text() == "o\nk"


def test_runner_takes_task_records_from_the_pool(c):
    # (MATLAB: matlabRunnerTakesTaskStringsFromThePool) a job's runner runs its task, then the
    # pool's other tasks in their order, each from the record its job got
    c.start_kind(task=Task("ptasks:copy", ["/nowhere", c.outs[10], str(c.log)], path=[str(c.root)]))
    ins = []
    for k in range(1, 4):
        ins.append(c.root / f"in_{k}.txt")
        ins[-1].write_text(f"task {k}\n")
        task = Task("ptasks:copy", [str(ins[-1]), c.outs[k], str(c.log)], path=[str(c.root)])
        c.submit(k, task=task)
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    c.run_job(1)
    assert c.tasks_run() == [1, 2, 3]
    assert [Path(c.outs[k]).read_text() for k in (1, 2, 3)] == [x.read_text() for x in ins]
    assert c.history() == matlab_history("matlabRunnerTakesTaskStringsFromThePool")


def test_pool_tasks_record_their_metrics(c):
    # (MATLAB: matlabPoolTasksRecordTheirMetrics) with task metrics, a task a runner runs
    # records itself as the task of its single job, joined to its submission
    from petakit5d.jobs.metrics import read_task_metrics

    kw = {"task_metrics": True, "cpus_per_task": 2, "mem_allocate": 30}
    c.start_kind(**kw)
    c.submit(7, **kw)
    c.submit(8, **kw)
    assert len(c.fake.calls("sbatch")) == 2
    c.fake.set_job_states("1001_8", ["PD"])
    c.run_job(1)
    t, _ = read_task_metrics(str(c.logs / "task_metrics"))
    t = sorted((x for x in t if x["task"] in (7, 8)), key=lambda x: x["task"])
    assert [x["task"] for x in t] == [7, 8]
    assert [(x["mode"], x["status"]) for x in t] == [("slurm_python", "ok")] * 2
    assert [(x["cpusPerTask"], x["memAllocate"]) for x in t] == [(2, 30)] * 2, (
        "The submissions were not joined."
    )
    assert c.history() == matlab_history("matlabPoolTasksRecordTheirMetrics")


def test_live_runner_keeps_the_task_of_an_ended_job(c):
    # after the own job of a task ended, a runner of another job that took the task may still
    # run it: the task counts as running (no new job) until that job ends; then it is
    # submitted again as before, and the master has it in the pool (no runner takes it)
    c.start_kind()
    c.submit(1)
    c.submit(2)
    c.claim(2, "1000_1:othernode:5:1")
    c.fake.set_job_states("1000_1", ["R", "R", "done"])
    for i in range(2):
        assert c.submit(2, job_id=1001) == (1001, 1, 0), f"Check {i + 1}."
    job_id, _, submitted = c.submit(2, job_id=1001)
    assert (job_id, submitted) == (1002, 1)
    cmds = c.fake.calls("sbatch")
    assert len(cmds) == 3
    assert as_before(cmds[2]), "The task was submitted again with task runners."
    owner, marked, _ = c.claim_of(2)
    assert owner.endswith(":0"), "The master did not take the task."
    assert marked, "The task is not marked taken."
    assert c.history() == matlab_history("liveRunnerKeepsTheTaskOfAnEndedJob")


def test_task_a_runner_finished_is_not_submitted_again(c):
    # (not in MATLAB) the runner that had the task finished it and its job ended: with the
    # task's output, the master sees the result after Slurm says both jobs are gone
    c.start_kind()
    c.submit(1)
    c.submit(2)
    c.claim(2, "1000_1:othernode:5:1")
    Path(c.outs[2]).write_text("ok")
    assert c.submit(2, job_id=1001, output=c.outs[2]) == (1001, -1, 0)
    assert c.fake.call_count("sbatch") == 2
    # (the result is looked at after Slurm's answers about both jobs)
    assert c.fake.calls("squeue") == ["-j 1001_2", "-j 1000_1", "-j 1000_1"]


def test_claim_of_an_owner_the_runners_did_not_write_is_live(c):
    c.start_kind()
    c.submit(1)
    c.submit(2)
    c.claim(2, "someone else")
    assert pools.live(c.key(2)) == (True, True)
    assert not pools.take(c.key(2))
    assert c.submit(2, job_id=1001) == (1001, 1, 0)
    assert c.fake.calls("squeue") == ["-j 1001_2"]


def test_task_whose_own_job_ended_unfinished_is_submitted_again(c):
    # the task's own job claimed it and ended without finishing it (killed, no failed mark):
    # the master takes it back and submits it as before, asking Slurm only about that job
    c.start_kind()
    c.submit(1)
    c.claim(1, "1000_1:testnode:4242:1")
    job_id, _, submitted = c.submit(1, job_id=1000)
    assert (job_id, submitted) == (1001, 1)
    assert as_before(c.fake.calls("sbatch")[-1])
    assert c.fake.calls("squeue") == ["-j 1000_1"]
    assert c.claim_of(1)[0].endswith(":0")


def test_runner_whose_job_is_pending_again_keeps_the_task(c):
    # (a runner's job that Slurm requeued) a claim of a job in the queue is live, pending too
    c.start_kind()
    c.submit(1)
    c.submit(2)
    c.claim(2, "1000_1:othernode:5:1")
    c.fake.set_job_states("1000_1", ["PD"])
    c.fake.set_job_states("1001_2", ["PD"])
    assert c.submit(2, job_id=1001) == (1001, 0, 0)
    assert c.claim_of(2)[0] == "1000_1:othernode:5:1", "The claim of a live runner was removed."


def test_a_new_process_runs_task_1_of_its_pool(c):
    # job 2's own task, then task 1 of the pool (its job pending) in a new process
    c.start_kind()
    for k in (1, 2):
        c.submit(k)
    c.fake.set_job_states("1000_1", ["PD"])
    r = c.run_job(2, {"PETAKIT5D_RUNNER_MINUTES": "0"})
    assert c.tasks_run() == [2, 1], r.stdout
    assert len({x["token"] for x in c.runs()}) == 2
    assert c.fake.calls("scancel") == ["--state=PENDING 1000_1"]


def test_runner_keeps_its_task_while_squeue_fails(c):
    # a failed squeue (a timeout of a busy controller) read as the job of a runner that took a
    # pool task gone, so the master took the task back; the runner keeps it now
    c.start_kind()
    c.submit(1)
    c.submit(2)
    c.claim(2, "1000_1:othernode:5:1")
    c.fake.set_job_states("1000_1", ["R"])
    c.fake.fail("squeue")
    assert pools.live(c.key(2)) == (True, True)
    assert not pools.take(c.key(2)), (
        "The master took the task of a runner Slurm didn't answer about."
    )
    assert c.claim_of(2) == ("1000_1:othernode:5:1", False, False)
    assert c.history() == matlab_history("runnerKeepsItsTaskWhileSqueueFails")


def test_pending_job_gets_back_the_task_of_a_runner_that_ended(c):
    # a runner of another job took the task of a pending job and ended without finishing it
    # (its job was killed): the master removes the claim, so that the pending job runs the
    # task when it starts (as without runners); nothing is submitted
    c.start_kind()
    c.submit(1)
    c.submit(2)
    c.claim(1, "77_1:othernode:5:1")  # (task 1 runs elsewhere)
    c.claim(2, "1000_1:othernode:5:1")
    c.fake.set_job_states("1001_2", ["PD"])
    assert c.submit(2, job_id=1001) == (1001, 0, 0)
    assert c.claim_of(2)[0] == "", "The claim of the runner that ended is left."
    c.run_job(2)
    assert c.tasks_run() == [2]
    assert c.fake.call_count("sbatch") == 2
    assert c.fake.call_count("scancel") == 0
    assert c.history() == matlab_history("pendingJobGetsBackTheTaskOfARunnerThatEnded")


def test_pending_job_of_a_task_a_runner_failed_is_replaced(c):
    # a runner of another job failed the task of a pending job: the master takes the task
    # back, cancels the pending job and submits the task again as before
    c.start_kind()
    for k in range(1, 4):
        c.submit(k)
    c.fake.set_job_states("1000_1", ["R"])  # (the runner's job still runs)
    for k in (2, 3):
        c.claim(k, "1000_1:othernode:5:1", failed=True)
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    job_id, _, submitted = c.submit(2, job_id=1001)
    assert (job_id, submitted) == (1003, 1)
    assert c.fake.calls("scancel") == ["1001_2"]
    assert as_before(c.fake.calls("sbatch")[-1])
    # the callers' status check of a last try: the task has ended
    assert quiet(check_single_job_status, 1002, 3) == -1
    assert c.fake.calls("scancel") == ["1001_2", "1002_3"]
    assert c.claim_of(3)[0].endswith(":0"), "The master did not take the task."
    assert c.history() == matlab_history("pendingJobOfATaskARunnerFailedIsReplaced")


def test_failed_task_is_submitted_again_as_before(c):
    c.start_kind()
    c.submit(1, mode="error")
    r = c.run_job(1, may_fail=True)
    assert r.returncode == 1
    assert c.pool().is_dir()
    owner, _, failed = c.claim_of(1)
    assert owner.startswith("1000_1:") and failed, "The failed task is not marked."
    job_id, _, submitted = c.submit(1, mode="error", job_id=1000)
    assert (job_id, submitted) == (1001, 1)
    cmds = c.fake.calls("sbatch")
    assert as_before(cmds[1]), cmds[1]
    # (the job of before runs it again; the task fails again)
    assert c.run_job(2, may_fail=True).returncode == 1
    assert c.tasks_run() == [1, 1]
    assert c.history() == matlab_history("failedTaskIsSubmittedAgainAsBefore")


def test_single_job_status_counts_task_runners(c):
    # check_single_job_status: a task a runner of another job has counts as running; once none
    # has, the master takes it from the pool (-1), so that no job's runner takes it after
    c.start_kind()
    for k in range(1, 4):
        c.submit(k)
    c.claim(2, "1000_1:othernode:5:1")
    c.fake.set_job_states("1000_1", ["R", "done"])
    c.fake.set_job_states("1002_3", ["PD"])
    assert quiet(check_single_job_status, 1001, 2) == 1
    assert quiet(check_single_job_status, 1001, 2) == -1
    assert quiet(check_single_job_status, 1002, 3) == 0
    # a job's task the master took: its job (or a runner) leaves it
    assert quiet(check_single_job_status, 1000, 1) == -1
    assert "no task left" in c.run_job(1).stdout
    assert c.runs() == []
    # a job of no pool: as check_slurm_job_status
    c.fake.set_job_states("77_4", ["R", "done"])
    assert [quiet(check_single_job_status, 77, 4) for _ in range(2)] == [1, -1]
    assert c.history() == matlab_history("singleJobStatusCountsTaskRunners")


def test_master_takes_the_last_file_from_the_pool(c):
    # master compute with the last file: the master cancels the file's pending job and runs it
    # itself, after it took the task from the pool; a task a runner of another job has stays
    # with it
    c.start_kind()
    c.submit(1)
    c.submit(2)
    c.fake.set_job_states("1001_2", ["PD"])
    c.fake.set_job_states("1000_1", ["R"])
    c.claim(2, "1000_1:othernode:5:1")
    last = {"master_compute": True, "last_file": True}
    assert c.submit(2, job_id=1001, **last)[1] == 1
    assert c.fake.call_count("scancel") == 0
    os.remove(c.pool() / "claims" / "t2")
    assert c.submit(2, job_id=1001, **last)[1] == 0
    assert c.fake.calls("scancel") == ["1001_2"]
    assert c.claim_of(2)[0].endswith(":0"), "The master did not take the task."
    assert c.fake.call_count("sbatch") == 2
    assert c.history() == matlab_history("masterTakesTheLastFileFromThePool")


def test_failed_submission_leaves_the_task_to_the_master(c):
    # a task whose job couldn't be submitted (sbatch failed each try) is taken back from the
    # pool; the caller's next try submits it as before
    c.start_kind()
    c.fake.reject_submissions(2)
    job_id, _, submitted = c.submit(1, max_trial_num=2)
    assert (job_id, submitted) == (-1, 0)
    assert c.claim_of(1)[0].endswith(":0"), "The master did not take the task."
    assert c.submit(1)[2] == 1
    assert as_before(c.fake.calls("sbatch")[-1])
    assert c.history() == matlab_history("failedSubmissionLeavesTheTaskToTheMaster")


# (the end's squeue also names the job of the kind's first task, whose line is in the folder of
# tasks that ran as before; MATLAB has no file for such a task)
FIRST_JOB = "1000_10"


def with_first_job(lines: list[str]) -> list[str]:
    out = []
    for x in lines:
        if x.startswith("squeue -h -r -o %i -j "):
            keys = sorted([*x.split()[-1].split(","), FIRST_JOB])
            x = "squeue -h -r -o %i -j " + ",".join(keys)
        out.append(x)
    return out


def test_ended_pool_goes_with_its_last_job(c):
    # at the end of a call (end(mark)) the master leaves the pools the call made, which take no
    # more tasks: a job of a pool still in the queue keeps it, and the last of its jobs to end
    # removes it; a later task of the kind makes a new pool
    mark = single.mark()
    c.start_kind()
    c.submit(1)
    c.submit(2)
    d = c.pool()
    c.fake.set_job_states("1001_2", ["R"])
    out = c.run_job(1).stdout
    quiet(single.end, mark)
    assert d.is_dir(), "The pool went while a job of it ran."
    out += c.run_job(2).stdout
    assert not d.exists(), out
    assert c.tasks_run() == [1, 2], out
    c.submit(3)
    assert c.pool() != d, "The ended pool took a task."
    assert c.history() == with_first_job(matlab_history("endedPoolGoesWithItsLastJob"))


def test_ended_pool_without_jobs_goes_at_once(c):
    # a pool whose jobs are no longer in the queue goes when the master leaves it, also with
    # the marks of jobs that never ran (a runner took their tasks and cancelled them)
    mark = single.mark()
    c.start_kind()
    for k in range(1, 4):
        c.submit(k)
    for k in (2, 3):
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    d = c.pool()
    c.run_job(1)
    assert c.tasks_run() == [1, 2, 3]
    for k in (2, 3):  # (cancelled)
        c.fake.set_job_states(f"{999 + k}_{k}", ["done"])
    assert d.is_dir(), "The pool went before the master left it."
    quiet(single.end, mark)
    assert not d.exists(), "The pool is still there."
    assert c.history() == with_first_job(matlab_history("endedPoolWithoutJobsGoesAtOnce"))


def test_pending_job_keeps_an_ended_pool(c):
    # a job whose task a runner took, but whose cancel was missed, is still pending when the
    # master leaves the pool: the pool stays until the job starts, ends at once (its task is
    # taken) and removes it
    mark = single.mark()
    c.start_kind()
    for k in (1, 2):
        c.submit(k)
    c.fake.set_job_states("1001_2", ["PD"])
    d = c.pool()
    c.run_job(1)
    assert c.tasks_run() == [1, 2]
    quiet(single.end, mark)
    assert d.is_dir(), "The pool went while a job of it was pending."
    out = c.run_job(2).stdout
    assert "no task left" in out
    assert not d.exists(), out
    assert len(c.runs()) == 2, "A task ran more than once."
    assert c.history() == with_first_job(matlab_history("pendingJobKeepsAnEndedPool"))


def test_pool_of_an_enclosing_call_stays(c):
    # a call that ends (an inner one, as a stitching wrapper the pipeline's master runs
    # itself) leaves only the pools it made
    c.start_kind()
    c.submit(1)
    d = c.pool()
    mark = single.mark()
    quiet(single.end, mark)
    assert (d / "live").is_dir()
    c.submit(2)
    assert (d / "runnable.txt").read_text() == "1\n2\n", "The pool took no more tasks."
    assert c.history() == matlab_history("poolOfAnEnclosingCallStays")


# --- more of the job side ---------------------------------------------------------------------


def test_failed_task_of_another_job_is_marked_and_the_job_takes_no_more(c):
    c.start_kind()
    c.submit(1)
    c.submit(2, mode="error")
    c.submit(3)
    for k in (2, 3):
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    r = c.run_job(1, may_fail=True)
    assert r.returncode == 1
    assert "RuntimeError: on purpose" in r.stdout + r.stderr
    assert c.tasks_run() == [1, 2]
    owner, _, failed = c.claim_of(2)
    assert owner.startswith("1000_1:") and failed
    assert c.claim_of(3)[0] == "", "The job took a task after a failed one."


def test_job_with_other_cpus_than_its_pool_takes_only_its_own_task(c):
    # (a job moved to another place with other CPUs: its results could differ from the
    # pool's other jobs')
    c.start_kind()
    for k in (1, 2):
        c.submit(k)
    c.fake.set_job_states("1001_2", ["PD"])
    c.run_job(1, {"SLURM_CPUS_PER_TASK": "8"})
    assert c.tasks_run() == [1]


def test_job_takes_tasks_of_others_only_while_it_has_the_time(c):
    c.start_kind()
    for k in (1, 2):
        c.submit(k)
    c.fake.set_job_states("1001_2", ["PD"])
    now = int(time.time())
    c.run_job(1, {"SLURM_JOB_START_TIME": str(now - 3600), "SLURM_JOB_END_TIME": str(now + 60)})
    assert c.tasks_run() == [1], "The job took a task without the time for it."
    assert c.fake.call_count("squeue") == 0


def test_job_leaves_tasks_the_master_took_back(c):
    c.start_kind()
    for k in (1, 2):
        c.submit(k)
    c.fake.set_job_states("1001_2", ["PD"])
    (c.pool() / "claims" / "r2").write_text("")
    c.run_job(1)
    assert c.tasks_run() == [1]


def test_two_jobs_share_the_tasks_of_a_pool(c):
    # tasks claimed by one job's runner are not run by another's: each task runs once
    c.start_kind()
    n = 8
    for k in range(1, n + 1):
        c.submit(k, seconds=0.2)
    for k in range(3, n + 1):
        c.fake.set_job_states(f"{999 + k}_{k}", ["PD"])
    procs = [c.run_job(k, background=True) for k in (1, 2)]
    for p in procs:
        out, _ = p.communicate(timeout=180)
        assert p.returncode == 0, out
    assert sorted(c.tasks_run()) == list(range(1, n + 1))
    assert {x["job"] for x in c.runs()} == {"1000", "1001"}


def test_job_state_survives_a_new_process(c):
    # the guards' new process keeps the job's tries: a failed task in a later process still
    # marks the job's exit status
    c.start_kind()
    c.submit(1)
    c.submit(2, mode="error")
    c.fake.set_job_states("1001_2", ["PD"])
    r = c.run_job(1, {"PETAKIT5D_RUNNER_MINUTES": "0"}, may_fail=True)
    assert r.returncode == 1
    assert c.tasks_run() == [1, 2]
    assert len({x["token"] for x in c.runs()}) == 2


# --- a real Slurm (marker cluster): PETAKIT5D_TEST_SLURM_CONFIG names a cluster config with
# the site's SlurmParam (partition, account, a small --mem-per-cpu, a job name), jobTimeLimit and
# PythonLaunchStr; PETAKIT5D_TEST_SLURM_SCRATCH a folder all nodes see -----------------------

real_slurm = slurm_fake.real_slurm  # (the fixture)


@pytest.fixture
def real_folder(real_slurm):
    import shutil
    import uuid

    scratch = os.environ.get("PETAKIT5D_TEST_SLURM_SCRATCH")
    if not scratch or not os.path.isdir(scratch):
        pytest.skip("PETAKIT5D_TEST_SLURM_SCRATCH names no folder")
    folder = Path(scratch) / f"pk5d_single_test_{uuid.uuid4().hex[:12]}"
    (folder / "out").mkdir(parents=True)
    (folder / "logs").mkdir()
    (folder / "ptasks.py").write_text(TASKS)
    yield folder
    shutil.rmtree(folder, ignore_errors=True)


@pytest.mark.cluster
@pytest.mark.slow
def test_real_cluster_single_jobs_with_half_held(real_folder):
    # 40 single jobs of one kind, submitted held; the first 20 released: they run their tasks
    # and take those of the held jobs (pending, reason JobHeldUser), which they cancel. Every
    # task runs once, no held job runs, none is left in the queue, and the pool goes with the
    # call
    from petakit5d.jobs.files import batch_file_exist
    from petakit5d.jobs.single import SingleJobs
    from petakit5d.jobs.slurm import Slurm

    conf = json.loads(Path(os.environ["PETAKIT5D_TEST_SLURM_CONFIG"]).read_text())
    conf["SlurmParam"] = f"{conf.get('SlurmParam', '')} --hold"
    cfg = real_folder / "held.json"
    cfg.write_text(json.dumps(conf))
    n, released = 40, 20
    logs, log = real_folder / "logs", real_folder / "runs.log"
    outs = [str(real_folder / "out" / f"o_{k}.txt") for k in range(1, n + 1)]
    tasks = [Task("ptasks:work", [o, str(log), "ok", 2.0], path=[str(real_folder)]) for o in outs]
    kw = {"job_log_dir": str(logs), "config_file": str(cfg), "max_trial_num": 1}
    slurm = Slurm()
    ids: set[str] = set()
    try:
        with SingleJobs(n) as jobs:
            for i in range(n):
                quiet(jobs.submit, i, tasks[i], i + 1, output=outs[i], **kw)
            assert all(j > 0 for j in jobs.job_ids), jobs.job_ids
            ids |= {f"{j}_{i + 1}" for i, j in enumerate(jobs.job_ids)}
            pool = list(logs.glob(f"{single.POOL_PREFIX}*"))
            assert len(pool) == 1
            first = {f"{jobs.job_ids[i]}_{i + 1}" for i in range(released)}
            r = slurm.run(["scontrol", "release", *sorted(first)])
            assert r.ok, r.output
            done = [False] * n
            t0 = time.monotonic()
            while not jobs.settled(done) and time.monotonic() - t0 < 900:
                found = batch_file_exist(outs, None, True)
                done = [d or f for d, f in zip(done, found, strict=True)]
                jobs.query(done)
                for i in range(n):
                    if not done[i]:
                        quiet(jobs.submit, i, tasks[i], output=outs[i], **kw)
                ids |= {f"{j}_{i + 1}" for i, j in enumerate(jobs.job_ids) if j > 0}
                coming = jobs.coming(done)
                if coming:
                    jobs.wait([[outs[i] for i in coming]], len(coming))
            assert all(done), f"{n - sum(done)} tasks not done"
            assert len(ids) == n, f"{len(ids) - n} tasks were submitted again"
        runs = [json.loads(x) for x in log.read_text().splitlines()]
        assert sorted(x["out"] for x in runs) == sorted(outs), "A task ran more than once."
        ran_in = {x["job"] for x in runs}
        assert ran_in <= {k.split("_")[0] for k in first}, "A held job ran."
        # (the held jobs were cancelled while pending, the released ones ended)
        deadline = time.monotonic() + 180
        while time.monotonic() < deadline:
            left = slurm.live_jobs(sorted(ids))
            if left is not None and not left and not pool[0].exists():
                break
            time.sleep(5)
        assert slurm.live_jobs(sorted(ids)) == set(), "Jobs of the test are left in the queue."
        assert not pool[0].exists(), "The pool is left."
        print(f"tasks ran in {len(ran_in)} of the {released} released jobs")
    finally:
        left = slurm.live_jobs(sorted(ids)) or set()
        if left:
            slurm.cancel(sorted(left))
