"""The cluster loop (petakit5d.jobs.cluster) against the fake Slurm, and on a real Slurm.

Ported from the MATLAB suite: the cluster-mode tests of GenericComputingFrameworkTest
(clusterModeSubmitsEachTaskOnceWithMaxJobNum, clusterModeCountsJobUntilAllItsTasksAreDone,
clusterModeFinishedShortLastBatchFreesItsJobSlot, configFileKeysIgnoreCase (its job part),
emptyJobLogDirInAConfigFileIsTheDefaultFolder) and the job-count and clean-up tests of
MccClusterWrapperTest (submitsEachTaskOnceWithMaxJobNum, countsJobUntilAllItsTasksAreDone,
finishedShortLastBatchFreesItsJobSlot, taskFilesThatCannotBeRemovedDoNotFailTheRun).

SCENARIOS are run by MATLAB's framework (oracle tree, the same fake Slurm, recorded by
PORT/agents/J03/scripts/j03_record.py) and by run_tasks here; the commands each gives the fake
(history.log, paths as placeholders, consecutive repeats collapsed, the --wrap text and the
Python job's --comment left out) must be the same (MATLAB_HISTORIES). Python runs them on a
virtual clock, so the 45 s grace periods take no time.

The real-cluster test (marker cluster) takes its partition, account and limits from the config
file PETAKIT5D_TEST_SLURM_CONFIG and a shared folder from PETAKIT5D_TEST_SLURM_SCRATCH.
"""

from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

import petakit5d
from petakit5d.jobs import cluster, framework
from petakit5d.jobs.cluster import job_shape
from petakit5d.jobs.framework import RUN_TASKS_PARAMS, run_tasks
from petakit5d.jobs.slurm import CommandResult, Slurm
from petakit5d.jobs.tasks import Task
from petakit5d.testing import slurm_fake
from petakit5d.testing.slurm_fake import SLURM_JOB_VARIABLES, FakeSlurm

fake_slurm = slurm_fake.fake_slurm  # (the fixture)

pytestmark = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")

SRC_ROOT = str(Path(petakit5d.__file__).resolve().parents[1])
MATLAB_TO_PYTHON = {p.matlab: p.python for p in RUN_TASKS_PARAMS}

# --- the scenarios (MATLAB names of the parameters; job states as FakeSlurm takes them, the
# files a job writes when it first runs as task numbers or "final") --------------------------

SLURM_PARAM = "-p cpu_queue -n1 --mem-per-cpu=1000M"
SCENARIOS: dict[str, dict] = {
    "first_submission_jobs_finish": {
        "n": 3,
        "params": {"masterCompute": False, "queryInterval": 0},
        "states": {
            "1000_1": (["PD", "R"], [1]),
            "1001_2": (["R"], [2]),
            "1002_3": (["PD", "PD", "R"], [3]),
        },
    },
    "pending_running_done_master": {
        "n": 4,
        "params": {"masterCompute": True, "queryInterval": 0},
        "states": {
            "1000_1": (["PD", "PD", "R"], [1]),
            "1001_2": (["PD", "R"], [2]),
            "1002_3": (["PD"], []),
        },
    },
    "retry_doubles_and_triples_cpus": {
        "n": 1,
        "params": {"masterCompute": False, "cpusPerTask": 2, "jobTimeLimit": 1.5},
        "states": {
            "1000_1": (["R", "done"], []),
            "1001_1": (["R", "done"], []),
            "1002_1": (["R"], [1]),
        },
    },
    "retry_gpu_jobs_keep_cpus": {
        "n": 2,
        "params": {"masterCompute": False, "cpusPerTask": 4, "GPUJob": True, "maxTrialNum": 2},
        "states": {
            "1000_1": (["R", "done"], [1]),
            "1001_2": (["R", "done"], []),
            "1002_2": (["PD", "R"], [2]),
        },
    },
    "max_job_num_holds_batches": {
        "n": 3,
        "params": {"masterCompute": False, "maxTrialNum": 1, "maxJobNum": 1, "queryInterval": 0},
        "states": {},
    },
    "job_counts_until_all_its_tasks_are_done": {
        "n": 4,
        "params": {
            "masterCompute": False,
            "maxTrialNum": 1,
            "maxJobNum": 1,
            "taskBatchNum": 2,
            "queryInterval": 0,
        },
        "states": {"1000_1": (["R", "R", "R", "done"], [1])},
    },
    "final_output_ends_the_run": {
        "n": 3,
        "final": True,
        "params": {"masterCompute": False, "queryInterval": 0},
        "states": {
            "1000_1": (["PD", "R"], ["final", 1]),
            "1001_2": (["PD"], []),
            "1002_3": (["PD"], []),
        },
    },
    "existing_outputs_and_missing_inputs": {
        "n": 4,
        "pre_outputs": [2],
        "missing_inputs": [3],
        "params": {"masterCompute": False, "queryInterval": 0},
        "states": {"1000_1": (["R"], [1]), "1001_4": (["R"], [4])},
    },
    "outputs_made_by_another_job": {
        "n": 2,
        "params": {"masterCompute": False, "queryInterval": 0},
        "states": {"1000_1": (["PD"], []), "1001_2": (["R"], [1, 2])},
    },
    "batches_of_tasks_with_master_compute": {
        "n": 5,
        "params": {"masterCompute": True, "taskBatchNum": 2, "queryInterval": 0},
        "states": {"1000_1": (["PD", "R"], [1, 2]), "1001_2": (["R"], [3, 4])},
    },
    "whole_node_jobs_retry_with_fewer_slots": {
        "n": 8,
        "mcc": True,
        "params": {
            "masterCompute": False,
            "wholeNodeJob": True,
            "GNUparallel": True,
            "cpusPerTask": 6,
            "jobTimeLimit": 1.5,
        },
        "states": {
            "1000_1": (["R", "done"], []),
            "1001_2": (["R", "done"], []),
            "1002_1": (["R", "done"], []),
            "1003_2": (["R", "done"], []),
            "1004_1": (["R"], [1, 2, 3, 4]),
            "1005_2": (["R"], [5, 6, 7, 8]),
        },
    },
    "whole_node_gpu_jobs": {
        "n": 4,
        "mcc": True,
        "params": {
            "masterCompute": False,
            "wholeNodeJob": True,
            "GNUparallel": True,
            "GPUJob": True,
            "GPUNum": 2,
            "cpusPerTask": 6,
            "SlurmParam": SLURM_PARAM + " --gres=gpu:2",
        },
        "states": {
            "1000_1": (["R", "done"], []),
            "1001_2": (["R", "done"], [3, 4]),
            "1002_1": (["R", "done"], []),
            "1003_1": (["R"], [1, 2]),
        },
    },
}

# The MATLAB framework's commands on the scenarios and its done flags (oracle tree
# PetaKit5D_7ed28c0, MATLAB R2024b, taskRunner false, recorded by j03_record.py of task J03):
# normalized, consecutive repeats collapsed, the job logs as <LOGS> and SLURM_PARAM as <P>
MATLAB_HISTORIES: dict[str, dict] = {
    "batches_of_tasks_with_master_compute": {
        "flags": [True, True, True, True, True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=2 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2",
            "scancel 1000_1",
            "squeue -j 1000_1",
            "scancel 1000 1001",
        ],
    },
    "existing_outputs_and_missing_inputs": {
        "flags": [True, True, True, True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=4 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_4",
            "scancel 1000 1001",
        ],
    },
    "final_output_ends_the_run": {
        "flags": [True, True, True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=2 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=3 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2,1002_3",
            "scancel 1000 1001 1002",
        ],
    },
    "first_submission_jobs_finish": {
        "flags": [True, True, True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=2 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=3 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2,1002_3",
            "squeue -j 1000_1,1002_3",
            "squeue -j 1002_3",
            "scancel 1000 1001 1002",
        ],
    },
    "job_counts_until_all_its_tasks_are_done": {
        "flags": [True, False, False, False],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1",
            "sbatch --array=2 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2",
            "scancel 1000 1001",
        ],
    },
    "max_job_num_holds_batches": {
        "flags": [False, False, False],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1",
            "sbatch --array=2 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2",
            "sbatch --array=3 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2,1002_3",
            "scancel 1000 1001 1002",
        ],
    },
    "outputs_made_by_another_job": {
        "flags": [True, True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=2 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2",
            "scancel 1000 1001",
        ],
    },
    "pending_running_done_master": {
        "flags": [True, True, True, True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=2 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=3 <LOGS> --cpus-per-task=1 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2,1002_3",
            "scancel 1002_3",
            "squeue -j 1000_1,1001_2,1002_3",
            "scancel 1000_1",
            "squeue -j 1000_1",
            "scancel 1000 1001 1002",
        ],
    },
    "retry_doubles_and_triples_cpus": {
        "flags": [True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=2 <P> -t 1:30:00 --wrap=<W>",
            "squeue -j 1000_1",
            "sbatch --array=1 <LOGS> --cpus-per-task=4 <P> -t 1:30:00 --wrap=<W>",
            "squeue -j 1001_1",
            "sbatch --array=1 <LOGS> --cpus-per-task=6 <P> -t 1:30:00 --wrap=<W>",
            "squeue -j 1002_1",
            "scancel 1002",
        ],
    },
    "retry_gpu_jobs_keep_cpus": {
        "flags": [True, True],
        "history": [
            "sinfo ",
            "sbatch --array=1 <LOGS> --cpus-per-task=4 <P> -t 24:00:00 --wrap=<W>",
            "sbatch --array=2 <LOGS> --cpus-per-task=4 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1000_1,1001_2",
            "squeue -j 1001_2",
            "sbatch --array=2 <LOGS> --cpus-per-task=4 <P> -t 24:00:00 --wrap=<W>",
            "squeue -j 1002_2",
            "scancel 1000 1002",
        ],
    },
    "whole_node_gpu_jobs": {
        "flags": [True, True, True, True],
        "history": [
            "sinfo ",
            (
                "sbatch --array=1 <LOGS> --cpus-per-task=24 --ntasks=1 <P> --gres=gpu:2 -t "
                "24:00:00 --wrap=<W slots=2 gpu>"
            ),
            (
                "sbatch --array=2 <LOGS> --cpus-per-task=24 --ntasks=1 <P> --gres=gpu:2 -t "
                "24:00:00 --wrap=<W slots=2 gpu>"
            ),
            "squeue -j 1000_1,1001_2",
            "squeue -j 1000_1",
            (
                "sbatch --array=1 <LOGS> --cpus-per-task=24 --ntasks=1 <P> --gres=gpu:2 -t "
                "24:00:00 --wrap=<W slots=2 gpu>"
            ),
            "squeue -j 1002_1",
            (
                "sbatch --array=1 <LOGS> --cpus-per-task=24 <P> --gres=gpu:2 -t 24:00:00 "
                "--wrap=<W slots=1>"
            ),
            "squeue -j 1003_1",
            "scancel 1001 1003",
        ],
    },
    "whole_node_jobs_retry_with_fewer_slots": {
        "flags": [True, True, True, True, True, True, True, True],
        "history": [
            "sinfo ",
            (
                "sbatch --array=1 <LOGS> --cpus-per-task=24 --ntasks=1 <P> -t 01:30:00 --wrap=<W "
                "slots=4>"
            ),
            (
                "sbatch --array=2 <LOGS> --cpus-per-task=24 --ntasks=1 <P> -t 01:30:00 --wrap=<W "
                "slots=4>"
            ),
            "squeue -j 1000_1,1001_2",
            (
                "sbatch --array=1 <LOGS> --cpus-per-task=24 --ntasks=1 <P> -t 01:30:00 --wrap=<W "
                "slots=2>"
            ),
            (
                "sbatch --array=2 <LOGS> --cpus-per-task=24 --ntasks=1 <P> -t 01:30:00 --wrap=<W "
                "slots=2>"
            ),
            "squeue -j 1002_1,1003_2",
            "sbatch --array=1 <LOGS> --cpus-per-task=24 <P> -t 01:30:00 --wrap=<W slots=1>",
            "sbatch --array=2 <LOGS> --cpus-per-task=24 <P> -t 01:30:00 --wrap=<W slots=1>",
            "squeue -j 1004_1,1005_2",
            "scancel 1004 1005",
        ],
    },
}


# --- running a scenario ------------------------------------------------------------------------


def touch(path: str) -> None:
    """a task: makes its output (by a temporary name and a rename)"""
    tmp = f"{path}_{os.getpid()}"
    Path(tmp).write_text("done")
    os.replace(tmp, path)


TOUCH = f"{__name__}:touch"


class VirtualClock:
    """the loop's clock and pause: time passes only in pauses, each a little longer than
    asked (as MATLAB's pause and the commands of a round take some time), so that a round
    after a pause of the query interval queries again, as MATLAB's rounds do"""

    OVERSHOOT = 0.01

    def __init__(self) -> None:
        self.now = 1000.0
        self.pauses: list[float] = []

    def clock(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.pauses.append(seconds)
        self.now += seconds + self.OVERSHOOT


@pytest.fixture(scope="module", autouse=True)
def forget_the_cluster_module():
    """after this module's tests, petakit5d.jobs as without the cluster module: a test that
    stands in a fake one through sys.modules (test_framework_local's seam test) gets it, as
    'from petakit5d.jobs import cluster' takes the package's attribute first"""
    yield
    import petakit5d.jobs as package

    if getattr(package, "cluster", None) is cluster:
        delattr(package, "cluster")
        if sys.modules.get("petakit5d.jobs.cluster") is cluster:
            del sys.modules["petakit5d.jobs.cluster"]


@pytest.fixture
def virtual_clock(monkeypatch):
    vc = VirtualClock()
    monkeypatch.setattr(cluster, "_clock", vc.clock)
    monkeypatch.setattr(cluster, "_sleep", vc.sleep)
    return vc


@pytest.fixture
def no_job_variables(monkeypatch):
    # the tests' own Slurm job must not count as the master's job (master_task_fits)
    for name in SLURM_JOB_VARIABLES:
        monkeypatch.delenv(name, raising=False)
    # (the stand-in Slurm check of run_tasks would call sinfo once more than MATLAB does)
    monkeypatch.setattr(framework, "_slurm_available", lambda: True, raising=False)


def scenario_files(spec: dict, folder: Path) -> tuple[list[str], list[str], str]:
    """the inputs and outputs of a scenario in folder, and its final output ('' without)"""
    n = spec["n"]
    (folder / "out").mkdir(parents=True, exist_ok=True)
    (folder / "logs").mkdir(exist_ok=True)
    ins = [str(folder / f"in_{k}.txt") for k in range(1, n + 1)]
    outs = [str(folder / "out" / f"out_{k}.txt") for k in range(1, n + 1)]
    for k, fn in enumerate(ins, 1):
        if k not in spec.get("missing_inputs", []):
            Path(fn).write_text(f"task {k}")
    for k in spec.get("pre_outputs", []):
        Path(outs[k - 1]).write_text("there")
    final = str(folder / "final.txt") if spec.get("final") else ""
    return ins, outs, final


def set_scenario_states(fake: FakeSlurm, spec: dict, outs: list[str], final: str) -> None:
    for key, (states, files) in spec["states"].items():
        paths = [final if f == "final" else outs[f - 1] for f in files]
        fake.set_job_states(key, states, paths or None)


def scenario_params(spec: dict, folder: Path) -> dict:
    """the scenario's parameters by Python name (MATLAB's names in spec)"""
    params = {"SlurmParam": SLURM_PARAM, "taskRunner": False, "jobLogDir": str(folder / "logs")}
    params.update(spec["params"])
    return {MATLAB_TO_PYTHON[k]: v for k, v in params.items()}


def collapse(lines: list[str]) -> tuple[list[str], list[int]]:
    """consecutive repeats as one line, and how many each was"""
    out: list[str] = []
    counts: list[int] = []
    for line in lines:
        if out and out[-1] == line:
            counts[-1] += 1
        else:
            out.append(line)
            counts.append(1)
    return out, counts


_WRAP = re.compile(r" --wrap=.*$", re.S)
_COMMENT = re.compile(r" --comment=\S+")


def normalize(line: str, folder: Path) -> str:
    """a history line with the folder as <D>, the wrap as <W> (slots and GPU of in-job
    slots kept: Python's --slots N --gpu, MATLAB's parallel --jobs N, or P=N; ... --jobs $P
    in GPU jobs, and its CUDA_VISIBLE_DEVICES, bash <task file> as one slot) and without the
    --comment"""
    line = line.replace(str(folder), "<D>")
    m = _WRAP.search(line)
    if m:
        wrap = m.group(0)
        tag = "<W>"
        slots = re.search(r"--slots (\d+)", wrap) or re.search(r"--jobs (\d+)", wrap)
        if slots is None and "--jobs $P " in wrap:
            slots = re.search(r"\bP=(\d+);", wrap)
        if slots:
            gpu = "--gpu" in wrap or "CUDA_VISIBLE_DEVICES" in wrap
            tag = f"<W slots={slots.group(1)}{' gpu' if gpu else ''}>"
        elif re.search(r"; bash \S*input_\d+\.txt", wrap):
            tag = "<W slots=1>"
        line = line[: m.start()] + f" --wrap={tag}"
    return _COMMENT.sub("", line)


def run_scenario(name: str, folder: Path) -> dict:
    """runs a scenario with run_tasks (fake Slurm first on PATH, a virtual clock); its flags
    and the collapsed, normalized history"""
    spec = SCENARIOS[name]
    folder.mkdir(parents=True, exist_ok=True)
    fake = FakeSlurm(folder / "fake").install()
    ins, outs, final = scenario_files(spec, folder)
    set_scenario_states(fake, spec, outs, final)
    params = scenario_params(spec, folder)
    if final:
        params["final_out_fullpath"] = final
    old_path = os.environ["PATH"]
    os.environ["PATH"] = fake.path_with(old_path)
    try:
        flags = run_tasks(ins, outs, [Task(TOUCH, [o]) for o in outs], **params)
    finally:
        os.environ["PATH"] = old_path
    lines, counts = collapse([normalize(x, folder) for x in fake.history()])
    return {"flags": flags, "history": lines, "counts": counts}


# --- the scenarios against MATLAB's recordings --------------------------------------------------


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_commands_equal_matlabs(name, tmp_path, virtual_clock, no_job_variables):
    if name not in MATLAB_HISTORIES:
        pytest.skip("no MATLAB recording of this scenario")
    got = run_scenario(name, tmp_path / name)
    want = MATLAB_HISTORIES[name]
    assert [short(x) for x in got["history"]] == want["history"]
    assert got["flags"] == want["flags"]


def short(line: str) -> str:
    """a history line with the job logs as <LOGS> and SLURM_PARAM as <P>"""
    logs = "-o <D>/logs/job_%A_%a.out -e <D>/logs/job_%A_%a.err"
    return line.replace(logs, "<LOGS>").replace(SLURM_PARAM, "<P>")


# --- the scenarios' own expectations (also without the recordings) ----------------------------


def sbatch_lines(history: list[str]) -> list[str]:
    return [x for x in history if x.startswith("sbatch ")]


def arrays(history: list[str]) -> list[int]:
    return [int(re.search(r"--array=(\d+)", x).group(1)) for x in sbatch_lines(history)]


def cpus(history: list[str]) -> list[int]:
    return [int(re.search(r"--cpus-per-task=(\d+)", x).group(1)) for x in sbatch_lines(history)]


def test_first_submission_one_job_per_batch(tmp_path, virtual_clock, no_job_variables):
    got = run_scenario("first_submission_jobs_finish", tmp_path)
    assert got["flags"] == [True] * 3
    assert got["history"][0] == "sinfo "
    first = sbatch_lines(got["history"])[0]
    assert first.startswith(
        "sbatch --array=1 -o <D>/logs/job_%A_%a.out -e <D>/logs/job_%A_%a.err "
        "--cpus-per-task=1 -p cpu_queue -n1 --mem-per-cpu=1000M -t 24:00:00 --wrap=<W>"
    )
    assert arrays(got["history"]) == [1, 2, 3]
    assert got["history"][-1] == "scancel 1000 1001 1002"
    # one batched squeue per round, of the jobs whose tasks are not done
    queries = [x for x in got["history"] if x.startswith("squeue")]
    assert queries == [
        "squeue -j 1000_1,1001_2,1002_3",
        "squeue -j 1000_1,1002_3",
        "squeue -j 1002_3",
    ]


def test_failed_job_is_retried_with_2x_and_3x_cpus(tmp_path, virtual_clock, no_job_variables):
    got = run_scenario("retry_doubles_and_triples_cpus", tmp_path)
    assert got["flags"] == [True]
    assert cpus(got["history"]) == [2, 4, 6]
    assert all(" -t 1:30:00 --wrap" in x for x in sbatch_lines(got["history"]))
    # the batch went again only 45 s after its job left the queue
    assert sum(virtual_clock.pauses) > 2 * cluster.GRACE_AFTER_JOB


def test_gpu_jobs_keep_their_cpus(tmp_path, virtual_clock, no_job_variables):
    got = run_scenario("retry_gpu_jobs_keep_cpus", tmp_path)
    assert got["flags"] == [True, True]
    assert cpus(got["history"]) == [4, 4, 4]
    assert arrays(got["history"]) == [1, 2, 2]


def test_master_takes_the_last_pending_batch(tmp_path, virtual_clock, no_job_variables):
    got = run_scenario("pending_running_done_master", tmp_path)
    assert got["flags"] == [True] * 4
    h = got["history"]
    # batch 4 was never submitted (the master ran it); the pending jobs of batches 3 and 1
    # were cancelled for the master to run them
    assert arrays(h) == [1, 2, 3]
    assert "scancel 1002_3" in h and "scancel 1000_1" in h
    assert h.index("scancel 1002_3") < h.index("scancel 1000_1")


def test_final_output_ends_the_run(tmp_path, virtual_clock, no_job_variables):
    got = run_scenario("final_output_ends_the_run", tmp_path)
    assert got["flags"] == [True] * 3
    assert got["history"][-1] == "scancel 1000 1001 1002"


def test_existing_outputs_and_missing_inputs_are_not_submitted(
    tmp_path, virtual_clock, no_job_variables
):
    got = run_scenario("existing_outputs_and_missing_inputs", tmp_path)
    assert got["flags"] == [True] * 4
    assert arrays(got["history"]) == [1, 4]


def test_outputs_made_by_another_job(tmp_path, virtual_clock, no_job_variables):
    # the job of batch 2 makes batch 1's output too: batch 1 is done, and its pending job is
    # cancelled at the end with the others (MATLAB's cancel of a pending job whose outputs
    # exist is never reached: such a batch counts as done before)
    got = run_scenario("outputs_made_by_another_job", tmp_path)
    assert got["flags"] == [True, True]
    assert got["history"][-1] == "scancel 1000 1001"
    assert "scancel 1000_1" not in got["history"]


def test_whole_node_jobs_retry_with_fewer_slots(tmp_path, virtual_clock, no_job_variables):
    got = run_scenario("whole_node_jobs_retry_with_fewer_slots", tmp_path)
    assert got["flags"] == [True] * 8
    s = sbatch_lines(got["history"])
    assert len(s) == 6 and cpus(got["history"]) == [24] * 6
    # 4 slots, then 2, then one task at a time (without --ntasks=1); -t as the compiled mode
    tags = [x.rsplit("--wrap=", 1)[1] for x in s]
    assert tags == ["<W slots=4>"] * 2 + ["<W slots=2>"] * 2 + ["<W slots=1>"] * 2
    assert all("--ntasks=1" in x for x in s[:4]) and not any("--ntasks" in x for x in s[4:])
    assert all(" -t 01:30:00 --wrap" in x for x in s)


def test_whole_node_gpu_jobs_have_a_slot_per_gpu(tmp_path, virtual_clock, no_job_variables):
    got = run_scenario("whole_node_gpu_jobs", tmp_path)
    assert got["flags"] == [True] * 4
    s = sbatch_lines(got["history"])
    assert [x.rsplit("--wrap=", 1)[1] for x in s] == ["<W slots=2 gpu>"] * 3 + ["<W slots=1>"]


# --- GenericComputingFrameworkTest (cluster mode) -------------------------------------------


def make_tasks(folder: Path, n: int) -> tuple[list[str], list[str]]:
    (folder / "out").mkdir(parents=True, exist_ok=True)
    ins = [str(folder / f"in_{i}.txt") for i in range(1, n + 1)]
    outs = [str(folder / "out" / f"out_{i}.txt") for i in range(1, n + 1)]
    for i, fn in enumerate(ins, 1):
        Path(fn).write_text(f"task {i}")
    return ins, outs


def noop_tasks(outs: list[str]) -> list[Task]:
    return [Task(f"{__name__}:noop", [o]) for o in outs]


def noop(path: str) -> None:
    """a task that makes nothing"""
    del path


def test_cluster_mode_submits_each_task_once_with_max_job_num(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    # with more jobs than maxJobNum the framework holds tasks back until earlier jobs end;
    # the fake jobs end at once without results
    ins, outs = make_tasks(tmp_path, 3)
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        max_job_num=1,
        job_log_dir=str(tmp_path),
        task_runner=False,
    )
    assert sorted(arrays(fake_slurm.history())) == [1, 2, 3]
    # the jobs start Python's job entry point, naming only the call folder and the batch
    wraps = [x.split("--wrap=", 1)[1] for x in fake_slurm.calls("sbatch")]
    assert all(re.search(r"-m petakit5d\.jobs\.exec \S+ [123]; echo; ", w) for w in wraps)


def test_cluster_mode_counts_job_until_all_its_tasks_are_done(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    ins, outs = make_tasks(tmp_path, 4)
    # the job of batch 1 (tasks 1 and 2) writes task 1's output at the first poll and is
    # gone at the fourth
    fake_slurm.set_job_states("1000_1", ["R", "R", "R", "done"], outs[:1])
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        max_job_num=1,
        task_batch_num=2,
        query_interval=0,
        job_log_dir=str(tmp_path),
    )
    cmds = fake_slurm.history()
    submit2 = [i for i, x in enumerate(cmds) if x.startswith("sbatch") and "--array=2 " in x]
    assert len(submit2) == 1, "Batch 2 should be submitted once."
    polls = [i for i, x in enumerate(cmds) if x.startswith("squeue") and "1000_1" in x]
    assert len(polls) >= 4, "The job of batch 1 was not polled until it ended."
    assert submit2[0] > polls[3], "Batch 2 was submitted while the job of batch 1 was running."


def finished_short_last_batch(fake_slurm, tmp_path, monkeypatch):
    # (Python: a refused submission is tried again in the next round, as the compiled mode
    # does; MATLAB jobs tried it once more at once, from a file)
    monkeypatch.setattr(cluster, "SUBMIT_UNKNOWN_ERRORS", ())
    ins, outs = make_tasks(tmp_path, 3)
    # Slurm refuses batch 1 (tasks 1 and 2), so batch 2 (task 3) is submitted first; its job
    # writes the output at the first poll and stays listed
    fake_slurm.reject_submissions(1)
    fake_slurm.set_job_states("1000_2", ["R"], outs[2:])
    final = tmp_path / "final.txt"
    # a framework that holds batch 1 back forever is stopped by a watchdog
    clock = {"rounds": 0}
    real_sleep = cluster._sleep

    def watchdog_sleep(seconds):
        clock["rounds"] += 1
        if clock["rounds"] > 200:
            final.write_text("watchdog")
        real_sleep(seconds)

    monkeypatch.setattr(cluster, "_sleep", watchdog_sleep)
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        max_job_num=1,
        task_batch_num=2,
        query_interval=0,
        job_log_dir=str(tmp_path),
        final_out_fullpath=str(final),
    )
    assert not final.exists(), "The finished batch 2 kept counting as a job."
    return arrays(fake_slurm.history())


def test_cluster_mode_finished_short_last_batch_frees_its_job_slot(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    assert finished_short_last_batch(fake_slurm, tmp_path, monkeypatch) == [1, 2, 1]


def test_existing_final_output_skips_all_tasks_in_cluster_mode(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    # as in local mode (existingFinalOutputSkipsAllTasks): nothing is submitted, not even
    # sinfo runs
    ins, outs = make_tasks(tmp_path, 3)
    final = tmp_path / "final.txt"
    final.write_text("done")
    flags = run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        job_log_dir=str(tmp_path),
        final_out_fullpath=str(final),
    )
    assert flags == [True] * 3 and fake_slurm.history() == []


def test_existing_outputs_skip_all_tasks_in_cluster_mode(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    ins, outs = make_tasks(tmp_path, 2)
    for o in outs:
        Path(o).write_text("there")
    flags = run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    assert flags == [True, True] and fake_slurm.history() == ["sinfo "]
    assert not list(tmp_path.glob("pk5d_call_*"))


def test_config_file_keys_ignore_case(fake_slurm, tmp_path, virtual_clock, no_job_variables):
    # 50 GB at 5 GB per CPU is 10 CPUs, rounded up to 12 to fit 24-CPU nodes
    ins, outs = make_tasks(tmp_path, 1)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"MemPerCPU": 5}))
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        mem_allocate=50,
        job_log_dir=str(tmp_path),
        config_file=str(config),
    )
    sbatch = fake_slurm.calls("sbatch")
    assert sbatch and "--cpus-per-task=12" in sbatch[0]


def test_empty_job_log_dir_in_a_config_file_is_the_default_folder(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # "jobLogDir": "" is the default folder: ../job_logs of the current folder if it is
    # there, else <data folder>/job_logs; the fake jobs make their logs as Slurm does
    fake_slurm.make_log_files()
    data = tmp_path / "data"
    data.mkdir()
    ins, outs = make_tasks(data, 2)
    config = data / "config.json"
    config.write_text(json.dumps({"jobLogDir": ""}))
    monkeypatch.chdir(data)
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        task_runner=False,
        config_file=str(config),
    )
    log_dir = str(data / "job_logs")
    sbatch = fake_slurm.calls("sbatch")
    assert len(sbatch) == 2
    for line in sbatch:
        m = re.search(r" -o (\S+) -e (\S+) ", line)
        assert m and m.group(1) == f"{log_dir}/job_%A_%a.out"
        assert m.group(2) == m.group(1).replace(".out", ".err")
    assert len(list((data / "job_logs").glob("job_*.out"))) == 2


# --- MccClusterWrapperTest (jobs of whole nodes) -------------------------------------------


SLOT_PARAMS = {"whole_node_job": True, "gnu_parallel": True, "cpus_per_task": 6}


def test_mcc_submits_each_task_once_with_max_job_num(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    ins, outs = make_tasks(tmp_path, 3)
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        max_job_num=1,
        job_log_dir=str(tmp_path),
        whole_node_job=True,
        cpus_per_task=24,
    )
    assert sorted(arrays(fake_slurm.history())) == [1, 2, 3]
    assert all("-m petakit5d.jobs.slots " in x for x in fake_slurm.calls("sbatch"))
    # one task at a time (24 CPUs per task): no --ntasks=1, no GNU parallel in MATLAB
    assert not any("--ntasks" in x for x in fake_slurm.calls("sbatch"))
    assert all("--slots 1 --tries 1" in x for x in fake_slurm.calls("sbatch"))


def test_mcc_counts_job_until_all_its_tasks_are_done(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    ins, outs = make_tasks(tmp_path, 4)
    fake_slurm.set_job_states("1000_1", ["R", "R", "R", "done"], outs[:1])
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        max_job_num=1,
        task_batch_num=2,
        query_interval=0,
        job_log_dir=str(tmp_path),
        whole_node_job=True,
        cpus_per_task=24,
    )
    cmds = fake_slurm.history()
    submit2 = [i for i, x in enumerate(cmds) if x.startswith("sbatch") and "--array=2 " in x]
    assert len(submit2) == 1
    polls = [i for i, x in enumerate(cmds) if x.startswith("squeue") and "1000_1" in x]
    assert len(polls) >= 4 and submit2[0] > polls[3]


def test_mcc_finished_short_last_batch_frees_its_job_slot(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    assert finished_short_last_batch(fake_slurm, tmp_path, monkeypatch) == [1, 2, 1]


def test_task_files_that_cannot_be_removed_do_not_fail_the_run(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # on NFS a job that still has a file of the call folder open keeps a .nfs file there,
    # so removing the folder fails; the run must not fail with it
    ins, outs = make_tasks(tmp_path, 2)
    fake_slurm.set_job_states("1000_1", ["R", "done"], outs[:1])
    fake_slurm.set_job_states("1001_2", ["R", "done"], outs[1:])
    removed = []
    real_rmtree = cluster.shutil.rmtree

    def refusing_rmtree(path, *args, **kwargs):
        if not Path(path).name.startswith("pk5d_call_"):
            return real_rmtree(path, *args, **kwargs)
        removed.append(path)
        raise OSError(39, "Directory not empty (a .nfs file)")

    monkeypatch.setattr(cluster.shutil, "rmtree", refusing_rmtree)
    with pytest.warns(UserWarning, match="could not be removed"):
        flags = run_tasks(
            ins,
            outs,
            noop_tasks(outs),
            master_compute=False,
            query_interval=0,
            job_log_dir=str(tmp_path),
        )
    assert flags == [True, True] and len(removed) == 1
    assert Path(removed[0]).name.startswith("pk5d_call_")


def test_the_call_folder_is_removed_when_all_tasks_are_done(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    ins, outs = make_tasks(tmp_path, 2)
    fake_slurm.set_job_states("1000_1", ["R"], outs[:1])
    fake_slurm.set_job_states("1001_2", ["R"], outs[1:])
    run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    assert not list(tmp_path.glob("pk5d_call_*"))


def submission_times(monkeypatch, vc) -> list[float]:
    """the virtual times of the sbatch calls from now on"""
    times: list[float] = []
    real = Slurm.sbatch

    def sbatch(self, argv):
        times.append(vc.now)
        return real(self, argv)

    monkeypatch.setattr(Slurm, "sbatch", sbatch)
    return times


def test_a_batch_goes_again_45_s_after_its_job_left(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # MATLAB's grace: 45 s after the last query that saw the job (the outputs of a job that
    # just ended may show late on NFS)
    ins, outs = make_tasks(tmp_path, 1)
    fake_slurm.set_job_states("1000_1", ["R", "done"], None)
    fake_slurm.set_job_states("1001_1", ["R"], outs)
    times = submission_times(monkeypatch, virtual_clock)
    run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    # submitted at 0, seen running at 6 (a doubled pause while all jobs pend), gone at 9
    assert len(times) == 2
    assert 45 + 6 <= times[1] - times[0] <= 45 + 6 + 3 + 0.5


def test_with_master_compute_the_grace_is_30_s_and_the_master_waits_30_s(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # the master ran batch 3 at once; the jobs of batches 1 and 2 end without outputs.
    # Batch 1 goes again 30 s after its job was last seen; the master takes batch 2, the
    # last waiting one, 30 s after its job left (max(30, half its last batch's time))
    ins, outs = make_tasks(tmp_path, 3)
    fake_slurm.set_job_states("1000_1", ["R", "done"], None)
    fake_slurm.set_job_states("1001_2", ["R", "done"], None)
    fake_slurm.set_job_states("1002_1", ["R"], outs[:1])
    times = submission_times(monkeypatch, virtual_clock)
    master_runs: list[tuple[float, list[int]]] = []
    real_run = framework.run_batch_in_process

    def run_batch(call, indices, names=(), mode="local"):
        master_runs.append((virtual_clock.now, list(indices)))
        return real_run(call, indices, names, mode)

    monkeypatch.setattr(framework, "run_batch_in_process", run_batch)
    flags = run_tasks(ins, outs, [Task(TOUCH, [o]) for o in outs], job_log_dir=str(tmp_path))
    assert flags == [True] * 3
    assert [x[1] for x in master_runs] == [[2], [1]]
    assert len(times) == 3
    # batch 1's job was last seen running at about 3 s
    assert 30 + 3 <= times[2] - times[0] <= 30 + 3 + 3 + 0.5
    assert 30 + 3 <= master_runs[1][0] - times[0] <= 30 + 3 + 3 + 0.5


def second_run(path: str, log_fn: str) -> None:
    """a task that makes its output only in its second run"""
    with open(log_fn, "a") as f:
        f.write("run\n")
    if Path(log_fn).read_text().count("run") >= 2:
        touch(path)


def test_the_masters_batch_done_in_its_last_try_counts_as_done(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    # (taskDoneInItsLastTryCountsAsDone in cluster mode) the master runs its batch again 30 s
    # later, in its last try; the loop ends right after it, and the last look at the outputs
    # counts it
    ins, outs = make_tasks(tmp_path, 2)
    fake_slurm.set_job_states("1000_1", ["R"], outs[:1])
    log_fn = str(tmp_path / "runs.log")
    t = [Task(f"{__name__}:noop", [outs[0]]), Task(f"{__name__}:second_run", [outs[1], log_fn])]
    flags = run_tasks(ins, outs, t, max_trial_num=2, job_log_dir=str(tmp_path))
    assert flags == [True, True]
    assert Path(log_fn).read_text().count("run") == 2
    assert fake_slurm.history()[-1] == "scancel 1000"


def test_the_pause_is_doubled_only_without_master_compute(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # all jobs pending: MATLAB waits another query interval only if the master doesn't compute
    ins, outs = make_tasks(tmp_path, 2)
    fake_slurm.set_job_states("1000_1", ["PD", "PD", "PD", "R"], outs[:1])
    queries: list[float] = []
    real = Slurm.batch_status

    def batch_status(self, job_ids, array_ids=None):
        queries.append(virtual_clock.now)
        return real(self, job_ids, array_ids)

    monkeypatch.setattr(Slurm, "batch_status", batch_status)
    flags = run_tasks(
        ins, outs, [Task(TOUCH, [o]) for o in outs], max_trial_num=1, job_log_dir=str(tmp_path)
    )
    assert flags == [True, True] and len(queries) == 5
    assert all(2.9 < b - a < 3.5 for a, b in zip(queries, queries[1:], strict=False))
    # without master compute: the pauses while all jobs pend are two intervals
    for o in outs:
        Path(o).unlink()
    fake_slurm.clear_calls()
    fake_slurm.set_job_states("1000_1", ["PD", "PD", "R"], outs[:1])
    fake_slurm.set_job_states("1001_2", ["PD", "PD", "R"], outs[1:])
    queries.clear()
    flags = run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        job_log_dir=str(tmp_path),
    )
    assert flags == [True, True]
    assert all(5.9 < b - a < 6.5 for a, b in zip(queries[:3], queries[1:3], strict=False))


def test_a_pause_ends_once_all_outputs_exist(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # the jobs make their outputs between two queries: the run ends without waiting for the
    # next query (MATLAB saw them only then)
    ins, outs = make_tasks(tmp_path, 2)
    for k in (1, 2):
        fake_slurm.set_job_states(f"{999 + k}_{k}", ["R"], None)
    marks = {}
    pause = virtual_clock.sleep

    def sleep(seconds):
        pause(seconds)
        if virtual_clock.now > 1010 and not marks:
            for o in outs:
                Path(o).write_text("done")
            marks["now"], marks["squeue"] = virtual_clock.now, fake_slurm.call_count("squeue")

    monkeypatch.setattr(cluster, "_sleep", sleep)
    flags = run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    assert flags == [True, True]
    assert fake_slurm.call_count("squeue") == marks["squeue"]
    assert virtual_clock.now - marks["now"] <= cluster.OUTPUT_LOOK_INTERVAL
    assert fake_slurm.history()[-1] == "scancel 1000 1001"


# --- unknown answers --------------------------------------------------------------------


def test_a_failed_query_keeps_the_job_states(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # squeue fails for a while: the running job is not taken as gone (MATLAB submitted it
    # again 45 s later)
    ins, outs = make_tasks(tmp_path, 1)
    fake_slurm.set_job_states("1000_1", ["R"], None)
    calls = {"n": 0}
    real = Slurm.batch_status

    def flaky(self, job_ids, array_ids=None):
        calls["n"] += 1
        if calls["n"] == 3:
            fake_slurm.fail("squeue")
        if calls["n"] == 40:
            fake_slurm.recover("squeue")
            Path(outs[0]).write_text("done")
        return real(self, job_ids, array_ids)

    monkeypatch.setattr(Slurm, "batch_status", flaky)
    flags = run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    assert flags == [True]
    assert fake_slurm.call_count("sbatch") == 1
    assert sum(virtual_clock.pauses) > 2 * cluster.GRACE_AFTER_JOB


def test_a_submission_without_answer_is_looked_for_before_it_is_submitted_again(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # sbatch times out, but Slurm made the job: a query finds it by its comment, and the
    # batch is not submitted again
    ins, outs = make_tasks(tmp_path, 1)
    made = {}
    real_run = Slurm.run

    def run(self, argv, timeout=None):
        argv = list(argv)
        if argv[0] == "sbatch" and not made:
            r = real_run(self, argv, timeout)
            made["id"] = int(r.stdout.split()[-1])
            made["comment"] = next(a for a in argv if a.startswith("--comment=")).split("=", 1)[1]
            return CommandResult(tuple(argv), None, "", "", timed_out=True)
        if argv[:2] == ["squeue", "--me"]:
            real_run(self, argv, timeout)  # (logged)
            line = f"{made['id']}_1|R|{made['comment']}\n"
            Path(outs[0]).write_text("done")
            return CommandResult(tuple(argv), 0, line, "")
        return real_run(self, argv, timeout)

    monkeypatch.setattr(Slurm, "run", run)
    flags = run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    assert flags == [True]
    assert fake_slurm.call_count("sbatch") == 1
    hist = fake_slurm.history()
    assert any(x.startswith("squeue --me -h -r -o %i|%t|%k") for x in hist)
    # its job is cancelled at the end with the others
    assert hist[-1] == f"scancel {made['id']}"


def test_nothing_is_submitted_while_a_lost_submission_cannot_be_looked_for(
    fake_slurm, tmp_path, virtual_clock, no_job_variables, monkeypatch
):
    # sbatch times out and then squeue fails for a while: the batch is not submitted again
    # until a query can say whether its job exists
    ins, outs = make_tasks(tmp_path, 1)
    made = {}
    lookups = {"n": 0}
    real_run = Slurm.run

    def run(self, argv, timeout=None):
        argv = list(argv)
        if argv[0] == "sbatch" and not made:
            r = real_run(self, argv, timeout)
            made["id"] = int(r.stdout.split()[-1])
            made["comment"] = next(a for a in argv if a.startswith("--comment=")).split("=", 1)[1]
            return CommandResult(tuple(argv), None, "", "", timed_out=True)
        if argv[:2] == ["squeue", "--me"]:
            lookups["n"] += 1
            if lookups["n"] <= 5:
                return CommandResult(tuple(argv), 1, "", "slurm_load_jobs error: timed out")
            Path(outs[0]).write_text("done")
            return CommandResult(tuple(argv), 0, f"{made['id']}_1|R|{made['comment']}\n", "")
        return real_run(self, argv, timeout)

    monkeypatch.setattr(Slurm, "run", run)
    flags = run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    assert flags == [True]
    assert fake_slurm.call_count("sbatch") == 1 and lookups["n"] == 6
    assert fake_slurm.history()[-1] == f"scancel {made['id']}"


def test_a_submission_without_answer_and_without_a_job_is_submitted_again(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    # the fake refuses with "Socket timed out" (after which Slurm may have made the job);
    # squeue doesn't list it, so the batch goes again
    ins, outs = make_tasks(tmp_path, 1)
    fake_slurm.reject_submissions(1)
    fake_slurm.set_job_states("1000_1", ["R"], outs)
    flags = run_tasks(ins, outs, noop_tasks(outs), master_compute=False, job_log_dir=str(tmp_path))
    assert flags == [True]
    hist = fake_slurm.history()
    lookups = [i for i, x in enumerate(hist) if x.startswith("squeue --me")]
    submits = [i for i, x in enumerate(hist) if x.startswith("sbatch")]
    assert len(submits) == 2 and lookups and submits[0] < lookups[0] < submits[1]


def test_every_job_gets_the_calls_comment(fake_slurm, tmp_path, virtual_clock, no_job_variables):
    ins, outs = make_tasks(tmp_path, 2)
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        job_log_dir=str(tmp_path),
    )
    comments = [re.search(r"--comment=(\S+)", x).group(1) for x in fake_slurm.calls("sbatch")]
    assert len(set(comments)) == 2
    assert all(re.fullmatch(r"pk5d_call_\d{8}_\d{6}_\d{3}_[0-9a-f-]{6}:[12]", c) for c in comments)


# --- cancelling the call's jobs --------------------------------------------------------------


def test_jobs_are_cancelled_when_the_master_raises(
    fake_slurm, tmp_path, virtual_clock, no_job_variables
):
    # the master's own batch raises twice: its jobs must not be left in the queue
    ins, outs = make_tasks(tmp_path, 3)
    fake_slurm.set_job_states("1000_1", ["PD"], None)
    fake_slurm.set_job_states("1001_2", ["PD"], None)
    t = noop_tasks(outs[:2]) + [Task(f"{__name__}:fail", [outs[2]])]
    with pytest.raises(RuntimeError, match="on purpose"):
        run_tasks(ins, outs, t, master_compute=True, job_log_dir=str(tmp_path))
    assert fake_slurm.history()[-1] == "scancel 1000 1001"


def fail(path: str) -> None:
    raise RuntimeError(f"on purpose: {path}")


MASTER_SCRIPT = """
import sys
from petakit5d.jobs import framework
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.tasks import Task

if hasattr(framework, "_slurm_available"):
    framework._slurm_available = lambda: True
folder = sys.argv[1]
ins = [f"{folder}/in_{i}.txt" for i in range(1, 4)]
outs = [f"{folder}/out/out_{i}.txt" for i in range(1, 4)]
print("master started", flush=True)
run_tasks(ins, outs, [Task("os:getcwd") for _ in outs], master_compute=False,
          job_log_dir=folder, query_interval=0.2)
"""


def test_jobs_are_cancelled_when_the_master_gets_sigterm(fake_slurm, tmp_path, no_job_variables):
    ins, outs = make_tasks(tmp_path, 3)
    for k in (1, 2, 3):
        fake_slurm.set_job_states(f"{999 + k}_{k}", ["PD"], None)
    script = tmp_path / "master.py"
    script.write_text(MASTER_SCRIPT)
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    for name in SLURM_JOB_VARIABLES:
        env.pop(name, None)
    proc = subprocess.Popen(
        [sys.executable, str(script), str(tmp_path)],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        deadline = time.monotonic() + 60
        while fake_slurm.call_count("squeue") < 2 and time.monotonic() < deadline:
            time.sleep(0.05)
        assert fake_slurm.call_count("sbatch") == 3
        proc.send_signal(signal.SIGTERM)
        out, err = proc.communicate(timeout=60)
    finally:
        if proc.poll() is None:
            proc.kill()
    assert proc.returncode == 128 + signal.SIGTERM, err
    assert fake_slurm.history()[-1] == "scancel 1000 1001 1002"


# --- job shapes, the wrap ----------------------------------------------------------------


def defaults(**kw) -> dict:
    s = {p.python: p.default for p in RUN_TASKS_PARAMS}
    s.update(kw)
    return s


# MATLAB's compiled mode on these (task number, its parameters): its task batch size, tasks at
# once and job number, and the --cpus-per-task of its jobs (masterCompute false; recorded by
# j03_shapes.py of task J03, oracle tree PetaKit5D_5df18b2)
W = {"wholeNodeJob": True, "GNUparallel": True}
MATLAB_SHAPES = {
    "whole_6": (100, dict(W, cpusPerTask=6), (4, 4, 25, 24)),
    "whole_6_min_jobs_50": (100, dict(W, cpusPerTask=6, minTaskJobNum=50), (2, 4, 50, 24)),
    "whole_5": (100, dict(W, cpusPerTask=5), (4, 4, 25, 24)),
    "whole_1_batch_10": (100, dict(W, cpusPerTask=1, taskBatchNum=10), (24, 24, 5, 24)),
    "whole_1_30_tasks": (30, dict(W, cpusPerTask=1), (24, 24, 2, 24)),
    "whole_6_mem_50": (40, dict(W, memAllocate=50, memPerCPU=10), (4, 4, 10, 24)),
    "whole_min_batch_7_para_batch_2": (
        50,
        dict(W, cpusPerTask=6, minBatchNum=7, paraBatchNum=2),
        (8, 4, 7, 24),
    ),
    "whole_max_cpu_16": (50, dict(W, cpusPerTask=3, maxCPUNum=16), (5, 5, 10, 16)),
    "whole_gpu_2": (10, dict(W, cpusPerTask=6, GPUJob=True, GPUNum=2), (2, 2, 5, 24)),
    "whole_gpu_4": (
        17,
        dict(W, cpusPerTask=4, GPUJob=True, GPUNum=4, taskBatchNum=5),
        (8, 4, 3, 24),
    ),
    "parallel_3_not_whole": (
        20,
        {"GNUparallel": True, "paraJobNum": 3, "cpusPerTask": 8, "taskBatchNum": 4},
        (6, 3, 4, 8),
    ),
    "extra_tasks_max_6": (
        100,
        dict(W, cpusPerTask=6, runExtraTasks=True, maxJobNum=6),
        (8, 4, 13, 24),
    ),
    "extra_tasks_order": (
        20,
        dict(W, cpusPerTask=12, runExtraTasks=True, maxJobNum=4),
        (4, 2, 5, 24),
    ),
    "extra_tasks_order_done": (
        24,
        dict(W, cpusPerTask=12, runExtraTasks=True, maxJobNum=3),
        (2, 2, 12, 24),
    ),
    "extra_tasks_small_batches": (
        60,
        dict(W, cpusPerTask=24, runExtraTasks=True, maxJobNum=50),
        (2, 1, 30, 24),
    ),
    "whole_1_10_tasks_extra_max_2": (
        10,
        dict(W, cpusPerTask=1, runExtraTasks=True, maxJobNum=2),
        (20, 24, 1, 24),
    ),
    "parallel_4_min_jobs_equal": (
        100,
        {"GNUparallel": True, "paraJobNum": 4, "cpusPerTask": 6, "minTaskJobNum": 100},
        (4, 4, 25, 6),
    ),
    "extra_tasks_big_batches": (
        3600,
        dict(W, cpusPerTask=24, taskBatchNum=600, runExtraTasks=True, maxJobNum=6),
        (600, 1, 6, 24),
    ),
}


@pytest.mark.parametrize("name", sorted(MATLAB_SHAPES))
def test_job_shape_equals_matlabs(name):
    n, params, (tb, para, jobs, cpus_) = MATLAB_SHAPES[name]
    s = defaults(**{MATLAB_TO_PYTHON[k]: v for k, v in params.items()})
    s["master_compute"] = False
    sh = job_shape(s, n)
    assert (sh.task_batch_num, sh.para_job_num, -(-n // sh.task_batch_num), sh.cpus_per_task) == (
        tb,
        para,
        jobs,
        cpus_,
    )


def test_job_shape_of_cpu_jobs():
    sh = job_shape(defaults(mem_allocate=50, mem_per_cpu=5), 100)
    assert (sh.cpus_per_task, sh.task_batch_num, sh.para_job_num, sh.slot_mode) == (12, 1, 1, False)
    # batches larger than the tasks shrink to them (minTaskJobNum 1)
    assert job_shape(defaults(task_batch_num=10), 4).task_batch_num == 4
    # GNU parallel with one task at a time is no slot mode
    assert not job_shape(defaults(gnu_parallel=True), 4).slot_mode


def test_job_shape_of_whole_node_jobs():
    s = defaults(whole_node_job=True, gnu_parallel=True, cpus_per_task=6)
    sh = job_shape(s, 100)
    assert (sh.cpus_per_task, sh.para_job_num, sh.task_batch_num, sh.slot_mode, sh.parallel) == (
        24,
        4,
        4,
        True,
        True,
    )
    # at least minTaskJobNum jobs
    s.update(min_task_job_num=50)
    assert job_shape(s, 100).task_batch_num == 2
    # GPU jobs: one task per GPU
    s.update(min_task_job_num=1, gpu_job=True, gpu_num=2)
    sh = job_shape(s, 100)
    assert (sh.para_job_num, sh.task_batch_num) == (2, 2)
    s.update(gpu_num=0)
    with pytest.raises(ValueError, match="gpu_num"):
        job_shape(s, 100)


def test_job_shape_with_extra_tasks():
    # runExtraTasks: as many jobs as maxJobNum + 1 (the master)
    s = defaults(whole_node_job=True, gnu_parallel=True, cpus_per_task=6, run_extra_tasks=True)
    sh = job_shape(dict(s, max_job_num=5), 100)
    assert sh.task_batch_num == 8 and -(-100 // sh.task_batch_num) == 13


def test_job_log_layout_of_the_wrap(fake_slurm, tmp_path, virtual_clock, no_job_variables):
    ins, outs = make_tasks(tmp_path, 1)
    run_tasks(
        ins,
        outs,
        noop_tasks(outs),
        master_compute=False,
        max_trial_num=1,
        job_log_dir=str(tmp_path),
        python_launch_str="export X=1; python3",
    )
    wrap = fake_slurm.calls("sbatch")[0].split("--wrap=", 1)[1]
    assert wrap.startswith('printf "Start time: "; date "+%Y-%m-%d %H:%M:%S.%6N"; ')
    assert wrap.endswith('; echo; printf "End time: "; date "+%Y-%m-%d %H:%M:%S.%6N"')
    assert "; export X=1; python3 -u -m petakit5d.jobs.exec " in wrap


# --- a real Slurm (marker cluster): PETAKIT5D_TEST_SLURM_CONFIG names a cluster config with
# the site's SlurmParam (partition, account, a small --mem-per-cpu, a job name), jobTimeLimit,
# maxJobNum and PythonLaunchStr; PETAKIT5D_TEST_SLURM_SCRATCH a folder all nodes see -----------

REAL_TASKS = '''
import os
import time
from datetime import datetime


def make(path, seconds=0.0):
    """writes its start and end time to path (by a temporary name and a rename)"""
    start = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")
    if seconds:
        time.sleep(seconds)
    tmp = f"{path}_{os.getpid()}"
    with open(tmp, "w") as f:
        f.write(start + "\\n" + datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f") + "\\n")
    os.replace(tmp, path)
'''

REAL_MASTER = """
import sys

sys.path.insert(0, sys.argv[1])
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.tasks import Task

folder, config, n, seconds = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
outs = [f"{folder}/out/out_{i}.txt" for i in range(1, n + 1)]
tasks = [Task("j03_real_tasks:make", [o, seconds], path=[folder]) for o in outs]
run_tasks([f"{folder}/j03_real_tasks.py"] * n, outs, tasks, task_batch_num=2,
          master_compute=False, job_log_dir=f"{folder}/logs", config_file=config)
"""


@pytest.fixture
def real_folder():
    import shutil
    import uuid

    scratch = os.environ.get("PETAKIT5D_TEST_SLURM_SCRATCH")
    if not scratch or not os.path.isdir(scratch):
        pytest.skip("PETAKIT5D_TEST_SLURM_SCRATCH names no folder")
    folder = Path(scratch) / f"pk5d_cluster_test_{uuid.uuid4().hex[:12]}"
    (folder / "out").mkdir(parents=True)
    (folder / "logs").mkdir()
    (folder / "j03_real_tasks.py").write_text(REAL_TASKS)
    yield folder
    shutil.rmtree(folder, ignore_errors=True)


def job_ids_of(text: str) -> list[int]:
    return [int(x) for x in re.findall(r"Submitted batch job (\d+)", text)]


def listed_jobs(job_ids: list[int]) -> set[int] | None:
    """the array jobs of job_ids that Slurm lists in any state (None: no answer)"""
    if not job_ids:
        return set()
    r = Slurm().run(["squeue", "-h", "-o", "%A", "-j", ",".join(str(j) for j in job_ids)])
    if not r.ok and "Invalid job id" not in r.output:
        return None
    return {int(x) for x in r.stdout.split() if x.isdigit()} & set(job_ids)


def wait_until_gone(job_ids: list[int], seconds: float = 180) -> list[int]:
    """the jobs Slurm still lists after waiting up to seconds for them to leave"""
    deadline = time.monotonic() + seconds
    while True:
        left = listed_jobs(job_ids)
        if left is not None and not left:
            return []
        if time.monotonic() > deadline:
            return sorted(left) if left is not None else [-1]
        time.sleep(5)


@pytest.mark.cluster
@pytest.mark.slow
def test_real_cluster_runs_tiny_tasks_in_batches(real_folder, capfd):
    # 200 tasks in batches of 10 (20 jobs; the config's maxJobNum caps the queue)
    config = os.environ["PETAKIT5D_TEST_SLURM_CONFIG"]
    n = 200
    outs = [str(real_folder / "out" / f"out_{i}.txt") for i in range(1, n + 1)]
    tasks = [Task("j03_real_tasks:make", [o], path=[str(real_folder)]) for o in outs]
    ins = [str(real_folder / "j03_real_tasks.py")] * n
    try:
        flags = run_tasks(
            ins,
            outs,
            tasks,
            task_batch_num=10,
            master_compute=False,
            job_log_dir=str(real_folder / "logs"),
            config_file=config,
        )
    finally:
        printed = capfd.readouterr().out
        ids = job_ids_of(printed)
    assert flags == [True] * n
    assert len(ids) == 20, printed[-3000:]
    # every job's log is in the job log folder, Start time first and End time last
    os.listdir(real_folder / "logs")
    for j in ids:
        logs = list((real_folder / "logs").glob(f"job_{j}_*.out"))
        assert len(logs) == 1, j
        lines = [x for x in logs[0].read_text().splitlines() if x.strip()]
        assert lines[0].startswith("Start time: ") and lines[-1].startswith("End time: "), lines
    # the call folder is gone, and so are the jobs
    assert not list((real_folder / "logs").glob("pk5d_call_*"))
    assert wait_until_gone(ids) == []


@pytest.mark.cluster
@pytest.mark.slow
def test_real_cluster_master_killed_with_sigterm_cancels_its_jobs(real_folder):
    config = os.environ["PETAKIT5D_TEST_SLURM_CONFIG"]
    script = real_folder / "master.py"
    script.write_text(REAL_MASTER)
    env = dict(os.environ, PYTHONPATH=SRC_ROOT)
    proc = subprocess.Popen(
        [sys.executable, "-u", str(script), str(real_folder), config, "40", "300"],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    printed: list[str] = []
    ids: list[int] = []
    try:
        deadline = time.monotonic() + 300
        while len(ids) < 20 and time.monotonic() < deadline and proc.poll() is None:
            line = proc.stdout.readline()
            printed.append(line)
            ids += job_ids_of(line)
        assert len(ids) == 20, "".join(printed[-50:])
        # let some of them start
        time.sleep(30)
        proc.send_signal(signal.SIGTERM)
        rest, _ = proc.communicate(timeout=180)
        printed.append(rest)
        assert proc.returncode == 128 + signal.SIGTERM, "".join(printed[-50:])
        assert wait_until_gone(ids) == []
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
        left = listed_jobs(ids)
        if left:
            # (only this test's own jobs, by their ids)
            Slurm().cancel([str(j) for j in sorted(left)])
