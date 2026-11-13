"""Single jobs (jobs.single): the submission of one task as a job of its own, the callers' status
checks and the callers' side of their loops (MATLAB's SingleJobSubmitWrapperTest, the single-job
test of CoordinatorNiceTest, and the single job of SlurmAlternativesTest through the wrapper).

The tests run against the fake Slurm of the MATLAB suite (petakit5d.testing.slurm_fake); jobs
whose commands run here run them as Slurm would, in a shell with Slurm's job variables."""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

import petakit5d
from petakit5d.jobs import alternatives as alt
from petakit5d.jobs import config, single
from petakit5d.jobs.single import (
    SINGLE_JOB_PARAMS,
    SingleJobs,
    SingleJobWarning,
    check_single_job_status,
    job_states,
    read_task_line,
    submit_single_job,
)
from petakit5d.jobs.tasks import Task
from petakit5d.testing import slurm_fake
from petakit5d.testing.slurm_fake import SLURM_JOB_VARIABLES, FakeSlurm, hide_slurm

fake_slurm = slurm_fake.fake_slurm  # (the fixture)

pytestmark = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")

SRC_ROOT = str(Path(petakit5d.__file__).resolve().parents[1])
# MATLAB's default SlurmParam (site values replaced, as in the fake's recordings)
ABC = "-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M"
GPU = "-p gpu8 --qos qos_normal -n1 --mem-per-cpu=7800M"
TASKS = '''
import os
import uuid


def work(out, log, mode="ok"):
    """appends its output's name to log, then writes out (by a temporary name and a rename)"""
    with open(log, "a") as fh:
        fh.write(out + "\\n")
    if mode == "error":
        raise RuntimeError("on purpose")
    tmp = f"{out}_{uuid.uuid4().hex}"
    with open(tmp, "w") as fh:
        fh.write(mode)
    os.replace(tmp, out)
'''


@pytest.fixture(autouse=True)
def fresh_session(tmp_path_factory, monkeypatch):
    """a new session (MATLAB: clear generic_single_job_submit_wrapper task_runner_run
    single_job_pools), Slurm hidden unless a test asks for the fake, no job variables"""
    hide_slurm(monkeypatch, tmp_path_factory.mktemp("no_slurm"))
    for name in (*SLURM_JOB_VARIABLES, "SLURM_CPUS_PER_TASK", "PETAKIT5D_RUNNER_STATE"):
        monkeypatch.delenv(name, raising=False)
    single.reset()
    alt.reset()
    config.clear_cache()
    alt.default_session().set_timing(view_age=0)
    yield
    single.reset()
    alt.reset()
    config.clear_cache()


def quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def make_task(folder: Path, k: int, mode: str = "ok") -> Task:
    (folder / "stasks.py").write_text(TASKS)
    (folder / "out").mkdir(exist_ok=True)
    return Task("stasks:work", [str(folder / "out" / f"o_{k}.txt"), str(folder / "runs.log"), mode],
                path=[str(folder)])  # fmt: skip


def log_folder(folder: Path) -> Path:
    """a new job log folder of the test's (MATLAB's withLogDir: createTempFolder per call)"""
    n = len(list(folder.glob("t[0-9]*"))) + 1
    d = folder / f"t{n}"
    d.mkdir()
    return d


def call_wrapper(folder: Path, job_id: int, task_id: int, **kwargs):
    """the wrapper for task 1 of folder (MATLAB's callWrapper: one sbatch try, MATLAB's default
    SlurmParam, a new job log folder unless one is given)"""
    if "job_log_dir" not in kwargs:
        kwargs["job_log_dir"] = str(log_folder(folder))
    kwargs.setdefault("max_trial_num", 1)
    kwargs.setdefault("slurm_param", ABC)
    return quiet(submit_single_job, make_task(folder, 1), job_id, task_id, **kwargs)


def submit_command(fake: FakeSlurm, folder: Path, **kwargs) -> str:
    """the sbatch command of a task without a job (MATLAB's submitCommand: the fake's log)"""
    call_wrapper(folder, -1, 7, **kwargs)
    cmds = fake.calls("sbatch")
    return cmds[-1] if cmds else ""


def wrap_of(cmd: str) -> str:
    return cmd.split("--wrap=", 1)[1]


def folder_of(cmd: str) -> tuple[str, int]:
    """the folder and task line a job's command names"""
    m = re.search(r"-m petakit5d\.jobs\.exec (\S+) (\d+)$", wrap_of(cmd))
    assert m, cmd
    return m.group(1).strip("'"), int(m.group(2))


def run_wrap(cmd: str, job_id: int, array_id: int, env=None) -> subprocess.CompletedProcess:
    """a job's --wrap command in a shell with Slurm's job variables"""
    e = dict(os.environ, PYTHONPATH=SRC_ROOT, SLURM_JOB_ID=str(job_id))
    e.update(SLURM_ARRAY_JOB_ID=str(job_id), SLURM_ARRAY_TASK_ID=str(array_id))
    e.update(SLURMD_NODENAME="testnode")
    e.update(env or {})
    return subprocess.run(
        ["sh", "-c", wrap_of(cmd)], capture_output=True, text=True, env=e, timeout=120
    )


def history(fake: FakeSlurm, folder: Path) -> list[str]:
    """the fake's commands with the test's folders (those in folder) as <T1>, <T2>, ... in order
    of first use, as in MATLAB's recordings, and the wrap as <W>"""
    tags: dict[str, str] = {}
    rx = re.compile(re.escape(str(folder)) + r"/[^/\s]+")

    def sub(m: re.Match) -> str:
        return tags.setdefault(m.group(0), f"<T{len(tags) + 1}>")

    out = []
    for line in fake.history():
        if line.startswith("sbatch ") and "--test-only" not in line:
            line = re.sub(r" --wrap=.*$", " --wrap=<W>", line)
        out.append(rx.sub(sub, line))
    return out


# --- MATLAB's histories (oracle tree PetaKit5D_7ed28c0, the recorded fake's history.log) -------
#
# The wrap of submissions is <W> (the job's command differs: Python's job entry point), the
# test's folders <T1>, <T2>, ...; MATLAB's default SlurmParam with the fake's made-up names
# ({ABC}). "sub <array> <folder>" is a submission of one CPU to the job log folder.

MATLAB_HISTORIES: dict[str, list[str]] = {
    "logPathsStayWhenSinfoFails": [
        "sbatch --array=7 -o <T1>/own_%A_%a.out -e <T1>/own_%A_%a.err --cpus-per-task=1 {ABC} "
        "-t 24:0:00 --wrap=<W>"
    ],
    "pendingJobReturnsWithoutWaiting": ["sinfo ", "squeue -j 123_7"],
    "runningJobReturnsWithoutWaiting": ["sinfo ", "squeue -j 123_7"],
    "jobNotInTheQueueIsSubmitted": ["sinfo ", "sub 7 <T1>"],
    "taskWithoutAJobIsSubmittedWithoutAQueueCheck": [
        "sinfo ",
        "sub 7 <T1>",
        "sinfo ",
        "squeue -j 1000_7",
    ],
    "slurmIsCheckedOncePerJobLogFolder": [
        "sinfo ",
        "sub 7 <T1>",
        "sub 8 <T1>",
        "sinfo ",
        "sub 9 <T2>",
    ],
    "sessionWithoutRunnersDoesNotLoadThem": ["sinfo ", "sub 7 <T1>", "sub 8 <T1>"],
    "failedQueryLeavesTheJob": ["sinfo ", "squeue -j 123_7", "sinfo ", "squeue -j 123_7"],
    "pendingLastJobGoesToTheMaster": ["sinfo ", "squeue -j 123_7", "scancel 123_7"],
    # CoordinatorNiceTest
    "singleJobsGetTheirNice": [
        "sinfo ",
        "sbatch --array=7 -o <T1>/job_%A_%a.out -e <T1>/job_%A_%a.err --cpus-per-task=1 {ABC} "
        "--nice=1000 -t 24:0:00 --wrap=<W>",
        "sub 8 <T1>",
    ],
}


def matlab_history(name: str) -> list[str]:
    out = []
    for x in MATLAB_HISTORIES[name]:
        f = x.split()
        if f[0] == "sub":
            x = (
                f"sbatch --array={f[1]} -o {f[2]}/job_%A_%a.out -e {f[2]}/job_%A_%a.err "
                "--cpus-per-task=1 {ABC} -t 24:0:00 --wrap=<W>"
            )
        out.append(x.replace("{ABC}", ABC))
    return out


def assert_matlab_history(fake: FakeSlurm, name: str, folder: Path) -> None:
    assert history(fake, folder) == matlab_history(name)


# --- SingleJobSubmitWrapperTest -----------------------------------------------------------------

CPU_CASES = {
    "cpusAboveMax": ({"cpus_per_task": 40}, 24),
    "minCPUNum": ({"min_cpu_num": 4}, 4),
    "memAllocate": ({"mem_allocate": 100}, 5),
    "memAllocateAboveMax": ({"mem_allocate": 1000}, 24),
}


@pytest.mark.parametrize("case", list(CPU_CASES))
def test_cpus_per_task_is_clamped_to_limits(case, tmp_path, fake_slurm):
    # fb6a98b regression: cpusPerTask/minCPUNum were only clamped to [minCPUNum, maxCPUNum]
    # when memAllocate was set
    kwargs, cpus = CPU_CASES[case]
    cmd = submit_command(fake_slurm, tmp_path, **kwargs)
    m = re.search(r"--cpus-per-task=(\d+)", cmd)
    assert m, "No sbatch command."
    assert int(m.group(1)) == cpus


def test_command_runs_the_task_from_its_folder(tmp_path, fake_slurm):
    # (MATLAB: commandWrapsTheTaskString) the job's command runs the job entry point on the
    # task's line in a folder of single jobs (the task never enters the command line), with
    # the launch string as given
    cmd = submit_command(
        fake_slurm, tmp_path, slurm_param="-p test_partition -n1", python_launch_str="pk_env/python"
    )
    assert cmd.startswith("--array=7 ")
    assert " -p test_partition -n1 " in cmd
    wrap = wrap_of(cmd)
    assert wrap.startswith("date; echo Python command: ")
    assert "stasks" not in wrap and "o_1.txt" not in wrap
    folder, i = folder_of(cmd)
    assert wrap.endswith(f"; pk_env/python -u -m petakit5d.jobs.exec {folder} {i}")
    assert os.path.basename(folder).startswith(single.TASK_FOLDER_PREFIX)
    line = read_task_line(folder, i)
    assert line["task"] == make_task(tmp_path, 1).record()
    assert (line["task_id"], line["language"]) == (7, "python")


def test_config_file_keys_ignore_case(tmp_path, fake_slurm):
    # the cluster configs in use spell MemPerCPU (it was matched case-sensitively)
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps({"MemPerCPU": 5, "PythonLaunchStr": "cfg_env/python"}))
    cmd = submit_command(fake_slurm, tmp_path, config_file=str(cfg), mem_allocate=50)
    assert "--cpus-per-task=10" in cmd
    assert " cfg_env/python -u -m petakit5d.jobs.exec " in cmd, "The config's launch was lost."


def test_empty_job_log_paths_in_a_config_file_are_the_defaults(tmp_path, fake_slurm, monkeypatch):
    # empty job log paths in a cluster config gave the job "-o  -e" (its log went to a file
    # named -e); they are the default ones (relative, as given: ../job_logs isn't there)
    work = tmp_path / "work" / "here"
    work.mkdir(parents=True)
    monkeypatch.chdir(work)
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps({"jobLogDir": "", "jobLogFname": "", "jobErrorFname": ""}))
    cmd = submit_command(
        fake_slurm, tmp_path, config_file=str(cfg), tmp_dir=str(tmp_path / "tasks")
    )
    m = re.search(r" -o (\S+) -e (\S+) ", cmd)
    assert m, cmd
    assert m.groups() == ("../job_logs/job_%A_%a.out", "../job_logs/job_%A_%a.err")
    assert folder_of(cmd)[0].startswith(str(tmp_path / "tasks"))


def test_log_paths_stay_when_sinfo_fails(tmp_path, fake_slurm, monkeypatch):
    # a job log folder that is there gives the job its log paths from the Slurm check, which
    # gives none if sinfo fails (as at times on a busy cluster): the job keeps its own
    stubs = tmp_path / "stub"
    stubs.mkdir()
    (stubs / "sinfo").write_text("#!/bin/sh\nexit 1\n")
    (stubs / "sinfo").chmod(0o755)
    monkeypatch.setenv("PATH", f"{stubs}{os.pathsep}{os.environ['PATH']}")
    logs = log_folder(tmp_path)
    with pytest.warns(UserWarning, match="slurm-based computing cluster"):
        call_wrapper(
            tmp_path,
            -1,
            7,
            job_log_dir=str(logs),
            job_log_fname=f"{logs}/own_%A_%a.out",
            job_error_fname=f"{logs}/own_%A_%a.err",
        )
    calls = fake_slurm.calls("sbatch")
    assert len(calls) == 1
    assert re.search(r" -o (\S+) -e (\S+) ", calls[0]).groups() == (
        f"{logs}/own_%A_%a.out",
        f"{logs}/own_%A_%a.err",
    )
    assert_matlab_history(fake_slurm, "logPathsStayWhenSinfoFails", tmp_path)


def test_pending_job_returns_without_waiting(tmp_path, fake_slurm):
    # the wrapper used to call squeue in a loop until a pending job started; callers check the
    # job again in their next sweep, as for running jobs
    fake_slurm.set_job_states("123_7", ["PD", "PD", "PD", "R"])
    assert call_wrapper(tmp_path, 123, 7) == (123, 0, 0)
    assert fake_slurm.call_count("squeue") == 1
    assert fake_slurm.call_count("sbatch") == 0
    assert_matlab_history(fake_slurm, "pendingJobReturnsWithoutWaiting", tmp_path)


def test_running_job_returns_without_waiting(tmp_path, fake_slurm):
    fake_slurm.set_job_states("123_7", ["R"])
    assert call_wrapper(tmp_path, 123, 7) == (123, 1, 0)
    assert fake_slurm.call_count("squeue") == 1
    assert_matlab_history(fake_slurm, "runningJobReturnsWithoutWaiting", tmp_path)


def test_job_not_in_the_queue_is_submitted(tmp_path, fake_slurm):
    job_id, _, submitted = call_wrapper(tmp_path, -1, 7)
    assert (job_id, submitted) == (1000, 1)
    assert fake_slurm.call_count("sbatch") == 1
    assert_matlab_history(fake_slurm, "jobNotInTheQueueIsSubmitted", tmp_path)


def test_task_without_a_job_is_submitted_without_a_queue_check(tmp_path, fake_slurm):
    # squeue ran also for a task without a job (job id -1), where it could only say so
    job_id, _, submitted = call_wrapper(tmp_path, -1, 7)
    assert (job_id, submitted) == (1000, 1)
    assert fake_slurm.call_count("squeue") == 0
    fake_slurm.set_job_states("1000_7", ["R"])
    assert call_wrapper(tmp_path, 1000, 7)[1] == 1
    assert fake_slurm.calls("squeue") == ["-j 1000_7"]
    assert_matlab_history(fake_slurm, "taskWithoutAJobIsSubmittedWithoutAQueueCheck", tmp_path)


def test_slurm_is_checked_once_per_job_log_folder(tmp_path, fake_slurm):
    # callers call the wrapper for each job in each round, and each call ran sinfo and listed
    # the job log folder: both are kept for the folder while Slurm is there
    log_dir, other = log_folder(tmp_path), log_folder(tmp_path)
    call_wrapper(tmp_path, -1, 7, job_log_dir=str(log_dir))
    call_wrapper(tmp_path, -1, 8, job_log_dir=str(log_dir))
    assert sum(x.startswith("sinfo") for x in fake_slurm.history()) == 1
    call_wrapper(tmp_path, -1, 9, job_log_dir=str(other))
    assert sum(x.startswith("sinfo") for x in fake_slurm.history()) == 2
    logs = [re.search(r"-o (\S+) -e (\S+)", x).groups() for x in fake_slurm.calls("sbatch")]
    assert len(logs) == 3
    assert logs[1] == logs[0] == (f"{log_dir}/job_%A_%a.out", f"{log_dir}/job_%A_%a.err")
    assert logs[2] == (f"{other}/job_%A_%a.out", f"{other}/job_%A_%a.err")
    assert_matlab_history(fake_slurm, "slurmIsCheckedOncePerJobLogFolder", tmp_path)


SESSION_SCRIPT = """
import contextlib, io, json, sys
sys.path.insert(0, sys.argv[1])
from petakit5d.jobs import single
from petakit5d.jobs.tasks import Task

folder, logs = sys.argv[2], sys.argv[3]
task = lambda k: Task("stasks:work", [f"{folder}/out/o_{k}.txt", f"{folder}/runs.log"],
                      path=[folder])
with contextlib.redirect_stdout(io.StringIO()):
    for k in (7, 8):
        single.submit_single_job(task(k), -1, k, job_log_dir=logs, max_trial_num=1,
                                 slurm_param=sys.argv[4], task_runner=False)
mods = sorted(m for m in sys.modules if m.startswith("petakit5d.jobs."))
print(json.dumps({"modules": mods}))
"""


def test_session_without_runners_does_not_load_them(tmp_path, fake_slurm):
    # with runners off, the submissions import neither the pools, the claims nor the runners
    make_task(tmp_path, 1)
    logs = log_folder(tmp_path)
    script = tmp_path / "session.py"
    script.write_text(SESSION_SCRIPT)
    r = subprocess.run(
        [sys.executable, str(script), SRC_ROOT, str(tmp_path), str(logs), ABC],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, r.stderr
    mods = set(json.loads(r.stdout.strip().splitlines()[-1])["modules"])
    assert fake_slurm.call_count("sbatch") == 2
    assert not {"petakit5d.jobs.pools", "petakit5d.jobs.claims", "petakit5d.jobs.runner"} & mods
    assert "petakit5d.jobs.alternatives" not in mods
    # (and the process left its folder when it ended: Slurm has none of its jobs)
    exit_query = "squeue -h -r -o %i -j 1000_7,1001_8"
    assert history(fake_slurm, tmp_path) == [
        *matlab_history("sessionWithoutRunnersDoesNotLoadThem"),
        exit_query,
    ]


def test_slurm_is_checked_in_each_call_without_slurm(tmp_path, monkeypatch):
    # without Slurm (sinfo fails) nothing is kept: each call checks again, as before
    stubs = tmp_path / "stub"
    stubs.mkdir()
    log = tmp_path / "sinfo.log"
    (stubs / "sinfo").write_text(f'#!/bin/sh\necho x >> "{log}"\nexit 1\n')
    (stubs / "sinfo").chmod(0o755)
    monkeypatch.setenv("PATH", f"{stubs}{os.pathsep}{os.environ['PATH']}")
    with pytest.warns(UserWarning):
        call_wrapper(tmp_path, -1, 7)
        call_wrapper(tmp_path, -1, 8)
    assert log.read_text().splitlines() == ["x", "x"]


def test_failed_query_leaves_the_job(tmp_path, fake_slurm):
    # a failed squeue (a timeout of a busy controller) read as the job gone: the task was
    # submitted again, or (master compute, last file) the caller ran it too; it counts as
    # running
    fake_slurm.set_job_states("123_7", ["R"])
    fake_slurm.fail("squeue", "slurm_load_jobs error: Socket timed out on send/recv operation")
    assert call_wrapper(tmp_path, 123, 7) == (123, 1, 0)
    assert call_wrapper(tmp_path, 123, 7, master_compute=True, last_file=True)[1] == 1
    assert fake_slurm.call_count("sbatch") == 0
    assert fake_slurm.call_count("scancel") == 0
    assert_matlab_history(fake_slurm, "failedQueryLeavesTheJob", tmp_path)


def test_pending_last_job_goes_to_the_master(tmp_path, fake_slurm):
    # with master compute the master cancels the last file's pending job and runs it itself
    fake_slurm.set_job_states("123_7", ["PD"])
    _, status, submitted = call_wrapper(tmp_path, 123, 7, master_compute=True, last_file=True)
    assert (status, submitted) == (0, 0)
    assert fake_slurm.calls("scancel") == ["123_7"]
    assert_matlab_history(fake_slurm, "pendingLastJobGoesToTheMaster", tmp_path)


def test_last_file_without_a_job_is_left_to_the_master(tmp_path, fake_slurm):
    assert call_wrapper(tmp_path, -1, 7, master_compute=True, last_file=True) == (-1, -1, 0)
    assert fake_slurm.call_count("sbatch") == 0


def test_array_index_is_the_task_id_mod_5000(tmp_path, fake_slurm):
    # (Slurm's MaxArraySize is 5001 by default)
    assert call_wrapper(tmp_path, -1, 5003)[2] == 1
    assert fake_slurm.calls("sbatch")[0].startswith("--array=3 ")
    fake_slurm.set_job_states("1000_3", ["R"])
    assert call_wrapper(tmp_path, 1000, 5003)[1] == 1
    assert quiet(check_single_job_status, 1000, 10003) == 1
    assert fake_slurm.calls("squeue") == ["-j 1000_3", "-j 1000_3"]


def test_finished_task_is_not_submitted_again(tmp_path, fake_slurm):
    # (not in MATLAB) a caller that looked for the result before the job ended would submit a
    # finished task again; with its output, the wrapper looks after Slurm says the job is gone
    out = tmp_path / "out" / "o_1.txt"
    assert call_wrapper(tmp_path, 1000, 7, output=str(out))[2] == 1, "No result: submitted."
    out.parent.mkdir(exist_ok=True)
    out.write_text("done")
    assert call_wrapper(tmp_path, 1000, 7, output=str(out)) == (1000, -1, 0)
    assert call_wrapper(tmp_path, 1000, 7, output=[str(out), str(out) + "x"])[2] == 1
    assert fake_slurm.call_count("sbatch") == 2
    fake_slurm.set_job_states("1002_7", ["PD"])
    assert call_wrapper(tmp_path, 1002, 7, output=str(out)) == (1002, 0, 0)


def test_failed_submissions_are_tried_max_trial_num_times(tmp_path, fake_slurm, capsys):
    fake_slurm.reject_submissions(5)
    out = submit_single_job(
        make_task(tmp_path, 1), -1, 7, job_log_dir=str(tmp_path), max_trial_num=3
    )
    assert out == (-1, -1, 0)
    assert fake_slurm.call_count("sbatch") == 3
    assert capsys.readouterr().out.count("Job submission failed, try again!") == 3


# --- CoordinatorNiceTest/singleJobsGetTheirNice ------------------------------------------------


def test_single_jobs_get_their_nice(tmp_path, fake_slurm):
    logs = str(log_folder(tmp_path))
    call_wrapper(tmp_path, -1, 7, job_log_dir=logs, nice=1000)
    call_wrapper(tmp_path, -1, 8, job_log_dir=logs)
    calls = fake_slurm.calls("sbatch")
    assert len(calls) == 2
    assert " --nice=1000" in calls[0]
    assert "--nice" not in calls[1], "A job without nice got one."
    assert_matlab_history(fake_slurm, "singleJobsGetTheirNice", tmp_path)


def test_nice_is_rounded(tmp_path, fake_slurm):
    call_wrapper(tmp_path, -1, 7, nice=2.5)
    assert " --nice=3 " in fake_slurm.calls("sbatch")[0]


# --- SlurmAlternativesTest/singleJobGoesToThePlaceWithRoomAndMoves, through the wrapper -------

ALT_HISTORY = [
    "sinfo ",
    "view",
    f"test --cpus-per-task=4 {ABC} -t 24:0:00",
    f"sbatch --array=3 -o <T1>/job_%A_%a.out -e <T1>/job_%A_%a.err --cpus-per-task=4 {ABC} "
    "-t 24:0:00 --wrap=<W>",
    "squeue -j 1000_3",
    "view",
    "test --cpus-per-task=4 -p gpu8 --qos qos_normal -n1 --mem-per-cpu=21418M -t 24:0:00",
    "update 1000_3 Partition=gpu8 QOS=qos_normal NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418",
    "squeue -j 1000_3",
    "view",
    f"test --cpus-per-task=4 {ABC} -t 24:0:00",
    "update 1000_3 Partition=cpu_queue QOS=qos_normal NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418",
]
VIEW = [
    "scontrol show node -o",
    "squeue --me -t PD -h -O JobArrayID:|,Partition:|,QOS:|,NumCPUs:|,tres-alloc:|,Reason:|",
]


def expand_alt(lines: list[str]) -> list[str]:
    out: list[str] = []
    for x in lines:
        if x == "view":
            out += VIEW
        elif x.startswith("test "):
            out.append(f"sbatch --test-only {x[5:]} --wrap=true")
        elif x.startswith("update "):
            out.append(f"scontrol update JobId={x[7:]}")
        else:
            out.append(x)
    return out


def node(name: str, partition: str, cpus: int, alloc: int, mb: int, mb_per: int, gpus=0):
    return FakeSlurm.node_line(name, partition, cpus, alloc, mb, alloc * mb_per, gpus, 0)


def test_single_job_goes_to_the_place_with_room_and_moves(tmp_path, fake_slurm):
    # the pipeline's and decon's single jobs (main place gpu8 without GPU headroom): submitted
    # to cpu_queue (gpu8 full) with a cpu_queue shape, then moved back to gpu8 once it has
    # room and cpu_queue not, by the wrapper's and check_single_job_status's looks at pending
    # jobs (the commands are MATLAB's)
    def nodes(cpu_alloc, gpu_alloc):
        fake_slurm.set_nodes(
            [
                node("n1", "cpu_queue", 24, cpu_alloc, 514048, 21418),
                FakeSlurm.node_line(
                    "g1", "gpu8", 128, gpu_alloc, 1031600, gpu_alloc * 8000, 8, 0, "MIXED", "gpu8"
                ),
            ]
        )

    nodes(0, 128)
    conf = {
        "SlurmParam": GPU,
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "gpuHeadroom": 0,
        "SlurmAlternatives": [{"SlurmParam": ABC, "memPerCPU": 20.9, "maxCPUNum": 24}],
    }
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(conf))
    args = {
        "job_log_dir": str(log_folder(tmp_path)),
        "mem_allocate": 70,
        "config_file": str(cfg),
        "task_runner": False,
    }
    task = make_task(tmp_path, 1)
    job_id = quiet(submit_single_job, task, -1, 3, **args)[0]
    jobs = [x for x in fake_slurm.calls("sbatch") if "--test-only" not in x]
    assert len(jobs) == 1
    assert f"--cpus-per-task=4 {ABC} " in jobs[0]
    fake_slurm.set_job_states(f"{job_id}_3", ["PD:Priority"])
    nodes(24, 0)
    alt.default_session().set_timing(min_wait=0, look_interval=0)
    assert quiet(submit_single_job, task, job_id, 3, **args)[2] == 0
    ups = [x for x in fake_slurm.calls("scontrol") if x.startswith("update")]
    assert ups == [
        f"update JobId={job_id}_3 Partition=gpu8 QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    ]
    nodes(0, 128)
    alt.default_session().set_timing(min_wait=0)
    quiet(check_single_job_status, job_id, 3)
    assert len([x for x in fake_slurm.calls("scontrol") if x.startswith("update")]) == 2
    assert history(fake_slurm, tmp_path) == expand_alt(ALT_HISTORY)


# --- the jobs of tasks that run as before -------------------------------------------------------


def test_job_runs_its_task_and_its_folder_goes_with_the_call(tmp_path, fake_slurm):
    # a job of a task that runs as before runs it in its process and leaves its folder; the
    # folder goes when the caller's call ended and its last job left
    mark = single.mark()
    call_wrapper(tmp_path, -1, 7)
    cmd = fake_slurm.calls("sbatch")[0]
    folder, i = folder_of(cmd)
    r = run_wrap(cmd, 1000, 7)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "out" / "o_1.txt").read_text() == "ok"
    assert "Task 7: stasks:work(" in r.stdout and "Elapsed time is " in r.stdout
    assert r.stdout.split("\n")[1].startswith("Python command: "), "The date is not first."
    assert os.path.isdir(folder), "The folder went before the caller left it."
    assert not os.path.exists(f"{folder}/live/t{i}")
    fake_slurm.set_job_states("1000_7", ["done"])
    single.end(mark)
    assert not os.path.exists(folder)


def test_folder_stays_while_a_job_of_it_is_in_the_queue(tmp_path, fake_slurm):
    # (MATLAB: matlabPoolGoesWithItsLastJob, for the folder of tasks that run as before) the
    # caller leaves the folder; its job, still in the queue, removes it when it ends
    mark = single.mark()
    call_wrapper(tmp_path, -1, 7)
    cmd = fake_slurm.calls("sbatch")[0]
    folder = folder_of(cmd)[0]
    fake_slurm.set_job_states("1000_7", ["R"])
    single.end(mark)
    assert os.path.isdir(folder)
    assert fake_slurm.calls("squeue")[-1] == "-h -r -o %i -j 1000_7"
    r = run_wrap(cmd, 1000, 7)
    assert r.returncode == 0, r.stderr
    assert not os.path.exists(folder), r.stdout
    # a later task makes a new folder
    call_wrapper(tmp_path, -1, 8, task_runner=False)
    assert folder_of(fake_slurm.calls("sbatch")[-1])[0] != folder


def test_folder_of_an_enclosing_call_stays(tmp_path, fake_slurm):
    # a call that ends (an inner one) leaves only the folders it made, not the one made just
    # before its mark
    logs = str(log_folder(tmp_path))
    call_wrapper(tmp_path, -1, 7, job_log_dir=logs, task_runner=False)
    folder = folder_of(fake_slurm.calls("sbatch")[0])[0]
    quiet(single.end, single.mark())
    assert os.path.exists(f"{folder}/live/m")
    call_wrapper(tmp_path, -1, 8, job_log_dir=logs, task_runner=False)
    assert folder_of(fake_slurm.calls("sbatch")[1]) == (folder, 2)


def test_failing_task_ends_its_job_with_status_1(tmp_path, fake_slurm):
    quiet(
        submit_single_job,
        make_task(tmp_path, 1, "error"),
        -1,
        7,
        job_log_dir=str(tmp_path),
        max_trial_num=1,
    )
    r = run_wrap(fake_slurm.calls("sbatch")[0], 1000, 7)
    assert r.returncode == 1
    assert "RuntimeError: on purpose" in r.stderr
    assert "Error in task 7 (stasks:work):" in r.stderr


def test_bash_task_runs_after_its_launch_string(tmp_path, fake_slurm):
    out = tmp_path / "bash.txt"
    quiet(
        submit_single_job,
        f'echo "$PK_X" > "{out}"',
        -1,
        2,
        language="bash",
        bash_launch_str="export PK_X=launched",
        job_log_dir=str(tmp_path),
    )
    cmd = fake_slurm.calls("sbatch")[0]
    assert "task_runner_single_" not in cmd
    r = run_wrap(cmd, 1000, 2)
    assert r.returncode == 0, r.stderr
    assert out.read_text() == "launched\n"


def test_bash_tasks_never_join_pools(tmp_path, fake_slurm):
    for k in (1, 2, 3):
        quiet(submit_single_job, f"true {k}", -1, k, language="bash", job_log_dir=str(tmp_path))
    assert not list(tmp_path.glob("task_runner_single_*"))
    assert len(fake_slurm.calls("sbatch")) == 3


def test_single_job_task_records_its_metrics(tmp_path, fake_slurm):
    # with task metrics, the submission is recorded in this process's folder and the task in
    # its job, and they are joined by the job
    from petakit5d.jobs.metrics import read_task_metrics

    logs = log_folder(tmp_path)
    call_wrapper(
        tmp_path, -1, 7, job_log_dir=str(logs), task_metrics=True, cpus_per_task=2, mem_allocate=30
    )
    r = run_wrap(fake_slurm.calls("sbatch")[0], 1000, 7)
    assert r.returncode == 0, r.stderr
    t, runs = read_task_metrics(str(logs / "task_metrics"))
    assert len(t) == 1
    assert (t[0]["task"], t[0]["mode"], t[0]["status"]) == (7, "slurm_python", "ok")
    assert (t[0]["cpusPerTask"], t[0]["memAllocate"]) == (2, 30)
    assert runs[0]["kind"] == "single"


def test_unwritable_task_folder_fails_the_submission(tmp_path, fake_slurm):
    # a task that can't be written anywhere is not submitted (the caller tries again); with
    # tmp_dir its folder goes there
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        pytest.skip("root writes into read-only folders")
    logs = tmp_path / "logs"
    logs.mkdir()
    logs.chmod(0o555)
    try:
        with pytest.warns(SingleJobWarning) as rec:
            out = call_wrapper(tmp_path, -1, 7, job_log_dir=str(logs), task_runner=False)
        assert out == (-1, -1, 0)
        assert rec[0].message.identifier == "generic_single_job_submit_wrapper:taskFile"
        assert fake_slurm.call_count("sbatch") == 0
        out = call_wrapper(tmp_path, -1, 7, job_log_dir=str(logs), tmp_dir=str(tmp_path / "tasks"))
        assert out[2] == 1
        assert folder_of(fake_slurm.calls("sbatch")[0])[0].startswith(str(tmp_path / "tasks"))
    finally:
        logs.chmod(0o755)


def test_wrong_parameters_raise(tmp_path):
    with pytest.raises(TypeError, match="'masterCompute'"):
        submit_single_job("x", -1, 1, master_compute=1)
    with pytest.raises(ValueError, match="cluster_type"):
        submit_single_job(make_task(tmp_path, 1), -1, 1, cluster_type="parfor")
    with pytest.raises(TypeError, match="needs language='bash'"):
        submit_single_job("echo x", -1, 1)


def test_parameters_are_matlabs_in_order():
    names = [p.matlab for p in SINGLE_JOB_PARAMS]
    assert names[:6] == [
        "clusterType",
        "parseCluster",
        "masterCompute",
        "lastFile",
        "jobLogDir",
        "jobLogFname",
    ]
    assert names[-4:] == ["taskRunner", "SlurmAlternatives", "gpuHeadroom", "nice"]
    defaults = {p.matlab: p.default for p in SINGLE_JOB_PARAMS}
    assert defaults["SlurmParam"] == "" and defaults["memPerCPU"] == 20.9
    assert defaults["masterCompute"] is False and defaults["jobTimeLimit"] == 24


# --- the callers' side (SingleJobs) -----------------------------------------------------------


def test_job_states_of_one_query(fake_slurm):
    fake_slurm.set_job_states("1000_1", ["R"])
    fake_slurm.set_job_states("1001_2", ["PD"])
    fake_slurm.set_job_states("1002_3", ["done"])
    out = job_states([1000, 1001, 1002, 1003, -1], [1, 2, 3, 4, 5])
    assert out == [1, 0, -1, -1, -1]
    assert fake_slurm.calls("squeue") == ["-j 1000_1,1001_2,1002_3,1003_4"]
    fake_slurm.fail("squeue")
    assert job_states([1000], [1]) == [1], "A query without an answer must count as running."
    assert job_states([], []) == []


def test_round_has_one_squeue_for_all_jobs(tmp_path, fake_slurm):
    # (the benchmark gate's rule) a round of a caller asks Slurm once for all its jobs, not
    # once per job, and the states it gets are those the wrapper would have got
    (tmp_path / "logs").mkdir()
    n = 6
    tasks = [make_task(tmp_path, k) for k in range(1, n + 1)]
    kw = {"job_log_dir": str(tmp_path / "logs"), "max_trial_num": 1, "task_runner": False}
    jobs = SingleJobs(n)
    done = [False] * n
    for i in range(n):
        assert quiet(jobs.submit, i, tasks[i], **kw) == -1
    assert jobs.job_ids == list(range(1000, 1000 + n))
    assert jobs.trials == [1] * n and all(jobs.active)
    assert fake_slurm.call_count("squeue") == 0
    for i in range(n):
        fake_slurm.set_job_states(f"{1000 + i}_{i + 1}", ["PD" if i % 2 else "R"])
    done[0] = True
    fake_slurm.clear_calls()
    jobs.query(done)
    states = [quiet(jobs.submit, i, tasks[i], **kw) for i in range(1, n)]
    assert states == [0, 1, 0, 1, 0]
    assert fake_slurm.calls("squeue") == ["-j 1001_2,1002_3,1003_4,1004_5,1005_6"]
    assert fake_slurm.call_count("sbatch") == 0
    # an item whose job changed since the query is looked at by itself
    fake_slurm.clear_calls()
    jobs.query(done)
    jobs.job_ids[1] = 1077
    fake_slurm.set_job_states("1077_2", ["R"])
    assert quiet(jobs.submit, 1, tasks[1], **kw) == 1
    assert fake_slurm.calls("squeue")[-1] == "-j 1077_2"


def test_last_try_only_checks_the_job(tmp_path, fake_slurm):
    jobs = SingleJobs(1, max_trial_num=1)
    task = make_task(tmp_path, 1)
    quiet(jobs.submit, 0, task, 7, job_log_dir=str(tmp_path), max_trial_num=1)
    assert jobs.trials == [1] and jobs.active == [True]
    fake_slurm.set_job_states("1000_7", ["PD", "done"])
    assert not jobs.settled([False])
    assert quiet(jobs.submit, 0, task, job_log_dir=str(tmp_path)) == 0
    assert jobs.coming([False]) == [0]
    assert quiet(jobs.submit, 0, task, job_log_dir=str(tmp_path)) == -1
    assert jobs.active == [False]
    assert jobs.settled([False]) and jobs.coming([False]) == []
    assert fake_slurm.call_count("sbatch") == 1
    jobs.ran_here(0)
    assert jobs.trials == [2]


def test_wait_between_rounds(monkeypatch):
    from petakit5d.jobs import files

    waits = []
    monkeypatch.setattr(files, "wait_for_files", lambda groups, t: waits.append(t) or False)
    jobs = SingleJobs(2)
    jobs.wait([["/nowhere"]], 1)
    jobs.wait([["/nowhere"]], 20)
    jobs.wait([["/nowhere"]], 1000)
    jobs.wait(None, 1, wait_time=0.5)
    assert waits == [3.0, 6.0, 30.0, 0.5]


def test_call_ends_its_pools_and_folders(tmp_path, fake_slurm):
    (tmp_path / "logs").mkdir()
    outer = single.mark()
    with SingleJobs(2) as jobs:
        for i in range(2):
            quiet(jobs.submit, i, make_task(tmp_path, i + 1), job_log_dir=str(tmp_path / "logs"))
        folders = sorted(p.name for p in (tmp_path / "logs").iterdir() if p.is_dir())
        assert len(folders) == 2
        assert folders[0].startswith(single.TASK_FOLDER_PREFIX)
        assert folders[1].startswith(single.POOL_PREFIX)
    # (the jobs never ran: Slurm no longer has them, so the folders went)
    assert not [p for p in (tmp_path / "logs").iterdir() if p.is_dir()]
    assert single.mark() > outer


def test_process_that_exits_leaves_its_folders(tmp_path, fake_slurm):
    # a process that ends without ending its calls leaves its folders as end() does
    (tmp_path / "logs").mkdir()
    make_task(tmp_path, 1)
    script = tmp_path / "exit.py"
    script.write_text(
        "import contextlib, io, sys\n"
        "sys.path.insert(0, sys.argv[1])\n"
        "from petakit5d.jobs import single\n"
        "from petakit5d.jobs.tasks import Task\n"
        "t = Task('stasks:work', [sys.argv[2] + '/out/o_1.txt', sys.argv[2] + '/runs.log'],\n"
        "         path=[sys.argv[2]])\n"
        "with contextlib.redirect_stdout(io.StringIO()):\n"
        "    single.submit_single_job(t, -1, 1, job_log_dir=sys.argv[2] + '/logs')\n"
    )
    r = subprocess.run(
        [sys.executable, str(script), SRC_ROOT, str(tmp_path)], capture_output=True, text=True
    )
    assert r.returncode == 0, r.stderr
    assert fake_slurm.call_count("sbatch") == 1
    assert not [p for p in (tmp_path / "logs").iterdir() if p.is_dir()]
