"""SlurmAlternatives (petakit5d.jobs.alternatives) against the fake Slurm, against MATLAB, and on a
real Slurm.

Ported from the MATLAB suite: SlurmAlternativesTest (all 35 methods, names in snake_case;
singleJobGoesToThePlaceWithRoomAndMoves at the level of the single jobs' places, which J07's
wrapper calls), SlurmAlternativesNiceTest (both methods) and
CoordinatorNiceTest/jobsAtAnotherPlaceKeepTheirNice. The tests on the fake Slurm also compare
its command history with MATLAB's for the same steps (MATLAB_HISTORIES: SlurmAlternativesTest
and SlurmAlternativesNiceTest of the oracle tree, recorded with an instrumented
FakeSlurmFixture (tasks J01, F04 and J09); the --wrap text of submissions as <W>, the test's
folder as <T1>; the --comment of submissions, MATLAB's and the cluster loop's call tags, is set
aside).

Oracle cases (ORACLE_CASES): the local functions of MATLAB's slurm_alternatives.m (extracted
from the oracle tree) on grids of inputs, and scenarios on a fake cluster (nodes, pending jobs
of the user, sbatch --test-only answers, choices, submissions, node changes and moves), run by
MATLAB's slurm_alternatives and by a Session against the same fake: the choices, the places,
the fake's command history and its jobs at the end must be equal.

Places of the tests: cpu_queue (24-CPU nodes, 20.9 GB per CPU) and gpu8 (128-CPU nodes with 8
GPUs, 7.6 GB per CPU; as an alternative it takes CPU-only jobs only with cpuJobs). A job of 70
GB is 4 CPUs at cpu_queue and 10 at gpu8 (16 in the batch framework, rounded to a factor of 32).
Place indices are 0-based (MATLAB's minus one).

The real-cluster test (marker cluster) takes its places from the config file
PETAKIT5D_TEST_SLURM_CONFIG (its SlurmParam and SlurmAlternatives) and submits only held jobs,
which it inspects with scontrol and cancels by their ids.
"""

from __future__ import annotations

import contextlib
import io
import json
import math
import os
import random
import re
import shlex
import subprocess
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from petakit5d.jobs import alternatives as alt
from petakit5d.jobs import cluster, config
from petakit5d.jobs.alternatives import AlternativesError, AlternativesWarning, Session
from petakit5d.jobs.cluster import JobShape
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.shapes import TIME_FORMAT_SINGLE, time_str
from petakit5d.jobs.slurm import Slurm, check_slurm_cluster
from petakit5d.jobs.tasks import Call, Task
from petakit5d.testing import oracle, slurm_fake
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.slurm_fake import SLURM_JOB_VARIABLES, FakeSlurm, hide_slurm

fake_slurm = slurm_fake.fake_slurm  # (the fixtures)
real_slurm = slurm_fake.real_slurm
POSIX = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")

ABC = "-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M"
GPU = "-p gpu8 --qos qos_normal -n1 --mem-per-cpu=7800M"
NOOP = f"{__name__}:noop"
TOUCH = f"{__name__}:touch"


def noop(path: str) -> None:
    """a task (the jobs of these tests never run)"""


def touch(path: str) -> None:
    """a task that makes its output (by a temporary name and a rename)"""
    tmp = f"{path}_{os.getpid()}"
    Path(tmp).write_text("done")
    os.replace(tmp, path)


# --- MATLAB's command histories (oracle tree PetaKit5D_7ed28c0) ---------------------------------
#
# A test's history is the commands of its last call: those that started after its last
# clearCalls. (The framework's normal end runs scancel in the background, which can write to the
# fake's log after the next clearCalls; recorded so, configWithoutAlternativesIsAsBefore started
# with its first call's scancel.) MATLAB's submissions carry its call tag
# (--comment=pk5d_<uuid12>_<batch>), the cluster loop's its own: the comparison sets both aside.
#
# Compact: "view" is scontrol show node -o and the squeue of the user's pending jobs, "test X"
# sbatch --test-only X --wrap=true, "update X" scontrol update JobId=X; {ABC}, {GPU}, {GPU_ABC}
# (gpu8 with cpu_queue's memory per CPU), {LOGS}, {M4}, and {CPU} and {GPU_CPU}
# (SlurmAlternativesNiceTest's places) as in ABBREV.

ABBREV = {
    "{ABC}": ABC,
    "{GPU_ABC}": "-p gpu8 --qos qos_normal -n1 --mem-per-cpu=21418M",
    "{GPU}": GPU,
    "{LOGS}": "-o <T1>/job_%A_%a.out -e <T1>/job_%A_%a.err",
    "{M4}": "NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418",
    "{CPU}": "-p cpu --qos normal -n1 --mem-per-cpu=21418M",
    "{GPU_CPU}": "-p gpu --qos normal -n1 --mem-per-cpu=21418M",
}
NODES_CMD = "scontrol show node -o"
PENDING_CMD = (
    "squeue --me -t PD -h -O JobArrayID:|,Partition:|,QOS:|,NumCPUs:|,tres-alloc:|,Reason:|"
)

MATLAB_HISTORIES: dict[str, list[str]] = {
    "memoryPerCPUOfAPlaceWithoutOneComesFromItsNodes": ["view"],
    "jobGoesToTheFirstPlaceWithRoom": [
        "view",
        "test --cpus-per-task=4 {ABC}",
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
    ],
    "jobKeepsItsShapeAtAPlaceWithLessMemoryPerCPU": ["view", "test --cpus-per-task=4 {GPU_ABC}"],
    "jobIsReshapedAtAPlaceWithMoreMemoryPerCPU": ["view", "test --cpus-per-task=4 {ABC}"],
    "cpuJobsDoNotGoToGPUNodesOfAnAlternativeWithoutCpuJobs": [
        "view",
        "view",
        "view",
        "test --cpus-per-task=4 {ABC}",
    ],
    "cpuJobsKeepOutOfGPUNodesWhereGPUJobsOfTheUserWait": [
        "view",
        "view",
        "update 1000_1 Partition=cpu_queue QOS=qos_normal {M4}",
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
    ],
    "mainPlaceOfGPUNodesTakesCPUJobs": ["view", "test --cpus-per-task=10 {GPU}"],
    "plannedNodesHaveNoRoom": ["view", "test --cpus-per-task=4 {GPU_ABC}"],
    "pendingJobsOfTheUserStartFirst": [
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "view",
        "test --cpus-per-task=4 {ABC}",
    ],
    "slurmsOwnCheckDecides": [
        "view",
        "test --cpus-per-task=4 {ABC}",
        "test --cpus-per-task=4 {GPU_ABC}",
    ],
    "placeSlurmRefusesIsNotUsed": ["view", "test --cpus-per-task=4 {GPU_ABC}", "view"],
    "mainPlaceIsNotLeftOutWhenSlurmRefusesItOnce": [
        "view",
        "test --cpus-per-task=4 {ABC}",
        "test --cpus-per-task=4 {GPU_ABC}",
        "view",
        "test --cpus-per-task=4 {ABC}",
    ],
    "burstOfJobsFillsTheRoomInOrder": [
        "view",
        "test --cpus-per-task=4 {ABC}",
        "test --cpus-per-task=4 {GPU_ABC}",
    ],
    "cpuJobsLeaveIdleGPUsTheirShare": [
        "view",
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
    ],
    "featuresAndNodeListsOfAPlaceCount": [
        "view",
        "view",
        "view",
        "view",
        "test --cpus-per-task=4 {GPU_ABC} -C gpu8 -x g1",
    ],
    "gpuJobsNeedAFreeGPU": [
        "view",
        "test --cpus-per-task=4 -p gpu9 --qos qos_normal -n1 --mem-per-cpu=8000M --gres=gpu:1",
    ],
    "limitThatPendingJobsWaitForClosesAPlace": ["view", "test --cpus-per-task=4 {GPU_ABC}"],
    "cpuJobsWaitOutsideGPUNodes": ["view"],
    "pendingJobMovesToAPlaceWithRoom": [
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "update 1000_1 Partition=gpu8 QOS=qos_normal {M4}",
    ],
    "onlyAsManyJobsMoveAsFit": [
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "update 1001_1 Partition=gpu8 QOS=qos_normal {M4}",
        "update 1002_2 Partition=gpu8 QOS=qos_normal {M4}",
    ],
    "jobsThatDoNotWaitForResourcesStay": ["view"],
    "jobWithRoomOnlyOutsideItsPlaceMoves": [
        "view",
        "test --cpus-per-task=4 {GPU_ABC} -x n1",
        "update 1000_1 Partition=gpu8 QOS=qos_normal {M4}",
    ],
    "pendingJobsAreFoundByTheirIds": [
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "update 1002_2 Partition=gpu8 QOS=qos_normal {M4}",
    ],
    "noSlurmCallsBeforeAJobHasWaited": [],
    "cpuJobLeavesTheQueueOfGPUNodes": [
        "view",
        "update 1000_1 Partition=cpu_queue QOS=qos_normal {M4}",
    ],
    "placesWithOtherOptionsDoNotTrade": ["view", "view"],
    "jobMovedBySomethingElseIsLeftToIt": ["view"],
    "failedMoveLeavesTheJob": [
        "view",
        "test --cpus-per-task=4 {GPU_ABC}",
        "update 1000_1 Partition=gpu8 QOS=qos_normal {M4}",
    ],
    "frameworkSubmitsToThePlaceWithRoom": [
        "sinfo ",
        "view",
        "test --cpus-per-task=4 {ABC} -t 2:00:00",
        "sbatch --array=1 --comment=pk5d_175c0ac1cd56_1 {LOGS} --cpus-per-task=4 {ABC} -t 2:00:00 "
        "--wrap=<W>",
        "view",
        "test --cpus-per-task=4 {ABC} -t 2:00:00",
        "sbatch --array=2 --comment=pk5d_175c0ac1cd56_2 {LOGS} --cpus-per-task=4 {ABC} -t 2:00:00 "
        "--wrap=<W>",
        "squeue -j 1000_1,1001_2",
        "scancel 1000 1001",
    ],
    "compiledFrameworkSubmitsToThePlaceWithRoom": [
        "sinfo ",
        "view",
        "test --cpus-per-task=4 {ABC} -t 24:00:00",
        "sbatch --array=1 --comment=pk5d_4fcba99a6897_1 {LOGS} --cpus-per-task=4 {ABC} "
        "-t 24:00:00 --wrap=<W>",
        "view",
        "test --cpus-per-task=4 {ABC} -t 24:00:00",
        "sbatch --array=2 --comment=pk5d_4fcba99a6897_2 {LOGS} --cpus-per-task=4 {ABC} "
        "-t 24:00:00 --wrap=<W>",
        "squeue -j 1000_1,1001_2",
        "scancel 1000 1001",
    ],
    "frameworkMovesPendingJobs": [
        "sinfo ",
        "view",
        "sbatch --array=1 --comment=pk5d_e782c4fd7c46_1 {LOGS} --cpus-per-task=4 {ABC} "
        "-t 24:00:00 --wrap=<W>",
        "squeue -j 1000_1",
        "view",
        "test --cpus-per-task=4 {GPU_ABC} -t 24:00:00",
        "update 1000_1 Partition=gpu8 QOS=qos_normal {M4}",
        *["squeue -j 1000_1", "view"] * 5,
        "squeue -j 1000_1",
        "scancel 1000",
    ],
    "singleJobGoesToThePlaceWithRoomAndMoves": [
        "sinfo ",
        "view",
        "test --cpus-per-task=4 {ABC} -t 24:0:00",
        "sbatch --array=3 {LOGS} --cpus-per-task=4 {ABC} -t 24:0:00 --wrap=<W>",
        "squeue -j 1000_3",
        "view",
        "test --cpus-per-task=4 {GPU_ABC} -t 24:0:00",
        "update 1000_3 Partition=gpu8 QOS=qos_normal {M4}",
        "squeue -j 1000_3",
        "view",
        "test --cpus-per-task=4 {ABC} -t 24:0:00",
        "update 1000_3 Partition=cpu_queue QOS=qos_normal {M4}",
    ],
    "pendingJobOfACallWithNiceMoves": [
        "sinfo ",
        "view",
        "sbatch --array=1 --comment=pk5d_6b16c0658c2f_1 {LOGS} --cpus-per-task=4 {CPU} --nice=1000 "
        "-t 24:00:00 --wrap=<W>",
        "squeue -j 1000_1",
        "view",
        "test --cpus-per-task=4 {GPU_CPU} --nice=1000 -t 24:00:00",
        "update 1000_1 Partition=gpu QOS=normal {M4}",
        *["squeue -j 1000_1", "view"] * 5,
        "squeue -j 1000_1",
        "scancel 1000",
    ],
    "jobsGetTheCallsNiceAtEveryPlace": [
        "sinfo ",
        "view",
        "test --cpus-per-task=4 {GPU_CPU} --nice=1000 -t 24:00:00",
        "sbatch --array=1 --comment=pk5d_3d7153916e5a_1 {LOGS} --cpus-per-task=4 {GPU_CPU} "
        "--nice=1000 -t 24:00:00 --wrap=<W>",
        "squeue -j 1000_1",
        "scancel 1000",
        "sinfo ",
        "view",
        "test --cpus-per-task=4 {GPU_CPU} --nice=5 -t 24:0:00",
        "sbatch --array=3 {LOGS} --cpus-per-task=4 {GPU_CPU} --nice=1000 -t 24:0:00 --wrap=<W>",
    ],
    "configWithoutAlternativesIsAsBefore": [
        "sinfo ",
        "sbatch --array=1 --comment=pk5d_u1_1 {LOGS} --cpus-per-task=4 {ABC} -t 24:00:00 "
        "--wrap=<W>",
        "sbatch --array=2 --comment=pk5d_u1_2 {LOGS} --cpus-per-task=4 {ABC} -t 24:00:00 "
        "--wrap=<W>",
        "squeue -j 1000_1,1001_2",
        "scancel 1000 1001",
    ],
}


def expand(lines: list[str]) -> list[str]:
    """the commands of a compact history"""
    out: list[str] = []
    for x in lines:
        if x == "view":
            out += [NODES_CMD, PENDING_CMD]
            continue
        for k, v in ABBREV.items():
            x = x.replace(k, v)
        if x.startswith("test "):
            x = f"sbatch --test-only {x[5:]} --wrap=true"
        elif x.startswith("update "):
            x = f"scontrol update JobId={x[7:]}"
        out.append(x)
    return out


_WRAP = re.compile(r" --wrap=.*$", re.S)
_COMMENT = re.compile(r" --comment=\S+")
_CALL_FOLDER = re.compile(r"pk5d_call_\d{8}_\d{6}_\d{3}_[0-9A-Za-z]+")


def history(fake: FakeSlurm, folder: Path | None = None) -> list[str]:
    """the fake's commands with the test's folder as <T1> and the wrap of submissions as <W>"""
    out = []
    for line in fake.history():
        if folder is not None:
            line = line.replace(str(folder), "<T1>")
        if line.startswith("sbatch ") and "--test-only" not in line:
            line = _WRAP.sub(" --wrap=<W>", line)
        out.append(line)
    return out


def comment_aside(lines: list[str]) -> list[str]:
    """commands without their --comment: a call's tag (MATLAB's pk5d_<uuid12>_<batch>, the
    cluster loop's pk5d_call_<...>:<batch>) differs in form and from call to call"""
    return [_COMMENT.sub("", x) for x in lines]


def assert_matlab_history(fake: FakeSlurm, name: str, folder: Path | None = None) -> None:
    assert comment_aside(history(fake, folder)) == comment_aside(expand(MATLAB_HISTORIES[name]))


# --- set-up and helpers (SlurmAlternativesTest's) ------------------------------------------------


@pytest.fixture(autouse=True)
def fresh_session(tmp_path_factory, monkeypatch):
    """a new session (MATLAB: clear slurm_alternatives and the wrappers), Slurm hidden unless a
    test asks for the fake, and each choice reading the cluster again (viewAge 0), as the MATLAB
    class sets up its tests"""
    hide_slurm(monkeypatch, tmp_path_factory.mktemp("no_slurm"))
    for name in SLURM_JOB_VARIABLES:
        monkeypatch.delenv(name, raising=False)
    alt.reset()
    config.clear_cache()
    alt.default_session().set_timing(view_age=0)
    yield
    alt.reset()
    config.clear_cache()


class VirtualClock:
    """the cluster loop's clock and pause: time passes only in pauses (a little longer than
    asked, as MATLAB's rounds take some time)"""

    def __init__(self, start: float = 1000.0) -> None:
        self.now = start

    def clock(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.now += seconds + 0.01


@pytest.fixture
def virtual_clock(monkeypatch):
    vc = VirtualClock()
    monkeypatch.setattr(cluster, "_clock", vc.clock)
    monkeypatch.setattr(cluster, "_sleep", vc.sleep)
    return vc


def S() -> Session:
    return alt.default_session()


def main_settings() -> dict:
    return {
        "SlurmParam": ABC,
        "SlurmConstraint": "",
        "memPerCPU": 20.9,
        "maxCPUNum": 24,
        "minCPUNum": 1,
        "gpuHeadroom": None,
    }


def abc_first(headroom, cpu_jobs: bool = True):
    """cpu_queue first, then gpu8 with the given GPU headroom, taking CPU jobs (unless
    cpu_jobs is False)"""
    entry = {
        "SlurmParam": GPU,
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "gpuHeadroom": headroom,
        "cpuJobs": cpu_jobs,
    }
    return S().places(main_settings(), entry)


def gpu_first():
    main = {
        "SlurmParam": GPU,
        "SlurmConstraint": "",
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "minCPUNum": 1,
        "gpuHeadroom": None,
    }
    return S().places(main, {"SlurmParam": ABC, "memPerCPU": 20.9, "maxCPUNum": 24})


def choose(places) -> int:
    """a job of 70 GB: 4 CPUs at cpu_queue, 10 at gpu8 (single-job shape)"""
    return S().choose(places, [4, 10], ["", ""])[0]


def submit_pending(fake: FakeSlurm, places, job: int, array: int, k: int) -> None:
    """a job submitted to place k (shaped for cpu_queue: 4 CPUs of 21418 MB) that waits for
    resources"""
    S().submitted(job, array, places, [4, 10], k, ["", ""])
    key = f"{job}_{array}"
    fake.set_job_states(key, ["PD:Priority"])
    fake.set_job_attributes(key, places[k].opts.partition, 4, 4 * 21418, 0, places[k].opts.qos)


def updates(fake: FakeSlurm) -> list[str]:
    return [x for x in fake.calls("scontrol") if x.startswith("update")]


def updated_jobs(fake: FakeSlurm) -> list[str]:
    return [re.sub(r"^update JobId=(\S+) .*$", r"\1", x) for x in updates(fake)]


def submissions(fake: FakeSlurm) -> list[str]:
    return [x for x in fake.calls("sbatch") if "--test-only" not in x]


def abc_node(name: str, cpu_alloc: int, state: str = "MIXED") -> str:
    return FakeSlurm.node_line(
        name, "cpu_queue", 24, cpu_alloc, 514048, cpu_alloc * 21418, 0, 0, state, "c24"
    )


def gpu_node(
    name: str, cpu_alloc: int, gpu_alloc: int, state: str = "MIXED", partition: str = "gpu8"
) -> str:
    return FakeSlurm.node_line(
        name, partition, 128, cpu_alloc, 1031600, cpu_alloc * 8000, 8, gpu_alloc, state, "gpu8"
    )


def abc_first_config(headroom) -> dict:
    return {
        "SlurmParam": ABC,
        "memPerCPU": 20.9,
        "maxCPUNum": 24,
        "SlurmAlternatives": [
            {
                "SlurmParam": GPU,
                "memPerCPU": 7.6,
                "maxCPUNum": 32,
                "gpuHeadroom": headroom,
                "cpuJobs": True,
            }
        ],
    }


def gpu_first_config() -> dict:
    """a config for CPU jobs on GPU nodes, with cpu_queue as an alternative"""
    return {
        "SlurmParam": GPU,
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "SlurmAlternatives": [{"SlurmParam": ABC, "memPerCPU": 20.9, "maxCPUNum": 24}],
    }


def write_config(folder: Path, conf: dict, name: str = "config.json") -> str:
    fn = folder / name
    fn.write_text(json.dumps(conf))
    return str(fn)


def make_tasks(folder: Path, n: int) -> tuple[list[str], list[str], list[Task]]:
    (folder / "out").mkdir(parents=True, exist_ok=True)
    ins = [str(folder / f"in_{k}.txt") for k in range(1, n + 1)]
    outs = [str(folder / "out" / f"out_{k}.txt") for k in range(1, n + 1)]
    for fn in ins:
        Path(fn).write_text("x")
    return ins, outs, [Task(NOOP, [o]) for o in outs]


def framework(ins, outs, tasks, folder: Path, config_file: str, **params) -> list[bool]:
    """run_tasks in cluster mode as SlurmAlternativesTest calls it"""
    settings = {
        "parse_cluster": True,
        "master_compute": False,
        "max_trial_num": 1,
        "mem_allocate": 70,
        "job_log_dir": str(folder),
        "task_runner": False,
        "config_file": config_file,
    }
    settings.update(params)
    with contextlib.redirect_stdout(io.StringIO()):
        return run_tasks(ins, outs, tasks, **settings)


def identifiers(record) -> list[str]:
    return [getattr(w.message, "identifier", "") for w in record]


# --- SlurmAlternativesTest -----------------------------------------------------------------------


def test_places_take_the_main_settings_for_keys_they_leave_out():
    main = main_settings()
    assert len(S().places(main, None)) == 1
    assert len(S().places(main, [])) == 1
    # (a JSON list of objects with other keys, as a config file gives it)
    alts = config.decode_json(
        f'{{"a": [{{"SlurmParam": "{GPU}", "memPerCPU": 7.6, "maxCPUNum": 32}}, '
        '{"slurmparam": "-p x --qos q --mem-per-cpu=32G"}]}'
    )["a"]
    places = S().places(main, alts)
    assert len(places) == 3
    assert [p.slurm_param for p in places] == [ABC, GPU, "-p x --qos q --mem-per-cpu=32G"]
    # (memPerCPU: an entry's, else its --mem-per-cpu)
    assert [p.mem_per_cpu for p in places] == [20.9, 7.6, 32]
    assert [p.max_cpu_num for p in places] == [24, 32, 24]
    assert [p.min_cpu_num for p in places] == [1, 1, 1]
    assert places[1].opts.partition == "gpu8"
    assert places[1].opts.qos == "qos_normal"
    assert places[1].opts.mem_per_cpu == 7800
    # (alternatives take CPU-only jobs at GPU nodes only if they say so)
    assert [p.cpu_jobs for p in places] == [True, False, False]
    # one entry (a JSON object) is the same; an unknown key is ignored with a warning
    with pytest.warns(AlternativesWarning) as record:
        places = S().places(main, {"SlurmParam": GPU, "maxCPU": 32})
    assert identifiers(record) == ["slurm_alternatives:key"]
    assert [p.max_cpu_num for p in places] == [24, 24]


@POSIX
def test_memory_per_cpu_of_a_place_without_one_comes_from_its_nodes(fake_slurm):
    # (the smallest RealMemory / CPUs of its nodes: 1031600 / 128 MB)
    fake_slurm.set_nodes(
        [
            abc_node("n1", 0),
            gpu_node("g1", 0, 0),
            FakeSlurm.node_line("g2", "gpu8", 64, 0, 1031600, 0, 4, 0),
        ]
    )
    places = S().places(main_settings(), {"SlurmParam": "-p gpu8 --qos qos_normal -n1 --mem=100G"})
    assert places[1].mem_per_cpu == 1031600 / 128 / 1024
    assert_matlab_history(fake_slurm, "memoryPerCPUOfAPlaceWithoutOneComesFromItsNodes")


def test_options_of_a_slurm_param_are_read():
    param = (
        "-A lab --partition=gpu4 -q high -n1 --mem=78G --gres=gpu:a100:2 -C c24 -x g3.cluster "
        "-J pk -t 1:00:00 --exclusive"
    )
    o = S().places(main_settings(), {"SlurmParam": param, "memPerCPU": 30})[1].opts
    assert (o.account, o.partition, o.qos, o.gres, o.constraint, o.exclude, o.job_name) == (
        "lab",
        "gpu4",
        "high",
        "gpu:a100:2",
        "c24",
        "g3.cluster",
        "pk",
    )
    assert (o.time, o.mem_token) == ("1:00:00", "--mem=78G")
    assert (o.mem, o.gpus, o.ntasks) == (78 * 1024, 2, 1)
    assert o.other == ("--exclusive",)


@POSIX
def test_job_goes_to_the_first_place_with_room(fake_slurm):
    places = abc_first(0)
    # cpu_queue has room: the main place
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    assert choose(places) == 0
    # cpu_queue is full: gpu8
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    assert choose(places) == 1
    # Slurm checked the place it chose
    assert any("--test-only" in x and "-p gpu8" in x for x in fake_slurm.calls("sbatch"))
    assert_matlab_history(fake_slurm, "jobGoesToTheFirstPlaceWithRoom")


@POSIX
def test_job_keeps_its_shape_at_a_place_with_less_memory_per_cpu(fake_slurm):
    # (shaped for cpu_queue: 4 CPUs of 21418 MB; at gpu8 it keeps them, not the 10 CPUs of
    # 7.8 GB a gpu8 config gives it)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    got = S().choose(abc_first(0), [4, 10], ["", ""])
    assert got == (1, 4, "-p gpu8 --qos qos_normal -n1 --mem-per-cpu=21418M")
    assert (
        "--cpus-per-task=4 -p gpu8 --qos qos_normal -n1 --mem-per-cpu=21418M"
        in fake_slurm.calls("sbatch")[-1]
    )
    assert_matlab_history(fake_slurm, "jobKeepsItsShapeAtAPlaceWithLessMemoryPerCPU")


@POSIX
def test_job_is_reshaped_at_a_place_with_more_memory_per_cpu(fake_slurm):
    # (shaped for gpu8: 10 CPUs of 7.8 GB; at cpu_queue 4 CPUs of 21418 MB)
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 128, 0)])
    assert S().choose(gpu_first(), [10, 4], ["", ""]) == (1, 4, ABC)
    assert_matlab_history(fake_slurm, "jobIsReshapedAtAPlaceWithMoreMemoryPerCPU")


@POSIX
def test_cpu_jobs_do_not_go_to_gpu_nodes_of_an_alternative_without_cpu_jobs(fake_slurm):
    # GPU nodes are for jobs that need a GPU: cpu_queue full, gpu8 idle (no headroom), but the
    # entry doesn't say cpuJobs: the job waits at cpu_queue, and a job waiting there doesn't
    # move to gpu8
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    places = abc_first(0, False)
    assert choose(places) == 0
    submit_pending(fake_slurm, places, 1000, 1, 0)
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1000], [1])
    assert updates(fake_slurm) == []
    # (an alternative of CPU and GPU nodes is not one of GPU nodes: with a GPU-node main place,
    # a CPU job goes there)
    fake_slurm.set_nodes(
        [
            FakeSlurm.node_line("n1", "cpu_queue", 24, 0, 514048, 0),
            FakeSlurm.node_line("t1", "cpu_queue", 20, 20, 515638, 0, 4, 4),
            gpu_node("g1", 128, 0),
        ]
    )
    assert S().choose(gpu_first(), [10, 4], ["", ""])[0] == 1
    assert_matlab_history(fake_slurm, "cpuJobsDoNotGoToGPUNodesOfAnAlternativeWithoutCpuJobs")


@POSIX
def test_cpu_jobs_keep_out_of_gpu_nodes_where_gpu_jobs_of_the_user_wait(fake_slurm):
    # (gpu8 takes CPU jobs and has room, but a GPU job of the user is pending there, held): new
    # CPU jobs wait at cpu_queue, and one waiting at gpu8 moves to cpu_queue, also where it
    # would fit
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    fake_slurm.set_job_states("900_1", ["PD:JobHeldUser"])
    fake_slurm.set_job_attributes("900_1", "gpu8", 4, 32000, 1, "qos_normal")
    places = abc_first(0)
    assert choose(places) == 0
    submit_pending(fake_slurm, places, 1000, 1, 1)
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1000], [1])
    assert updates(fake_slurm) == [
        "update JobId=1000_1 Partition=cpu_queue QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    ]
    # (without GPU jobs pending there it goes to gpu8, as before)
    fake_slurm.set_job_states("900_1", ["done"])
    assert choose(places) == 1
    assert_matlab_history(fake_slurm, "cpuJobsKeepOutOfGPUNodesWhereGPUJobsOfTheUserWait")


@POSIX
def test_main_place_of_gpu_nodes_takes_cpu_jobs(fake_slurm):
    # (the main place is the user's choice: used as before)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    main = {
        "SlurmParam": GPU,
        "SlurmConstraint": "",
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "minCPUNum": 1,
        "gpuHeadroom": 0,
    }
    places = S().places(main, {"SlurmParam": ABC, "memPerCPU": 20.9})
    assert S().choose(places, [10, 4], ["", ""])[0] == 0
    assert_matlab_history(fake_slurm, "mainPlaceOfGPUNodesTakesCPUJobs")


@POSIX
def test_planned_nodes_have_no_room(fake_slurm):
    # a node the backfill scheduler planned for a pending job keeps its free CPUs for it
    fake_slurm.set_nodes([abc_node("n1", 16, "MIXED+PLANNED"), gpu_node("g1", 0, 0)])
    assert choose(abc_first(0)) == 1
    assert_matlab_history(fake_slurm, "plannedNodesHaveNoRoom")


@POSIX
def test_pending_jobs_of_the_user_start_first(fake_slurm):
    # cpu_queue has 8 free CPUs, which a pending job of the user takes
    fake_slurm.set_nodes([abc_node("n1", 16), gpu_node("g1", 0, 0)])
    fake_slurm.set_job_states("900_1", ["PD:Priority"])
    fake_slurm.set_job_attributes("900_1", "cpu_queue", 8, 8 * 21418)
    assert choose(abc_first(0)) == 1
    # one that fits no node waits, and a new job would wait behind it (Slurm starts them by
    # priority): no room at cpu_queue either
    fake_slurm.set_job_attributes("900_1", "cpu_queue", 32, 1000)
    assert choose(abc_first(0)) == 1
    # (only waiting jobs: a held one doesn't count)
    fake_slurm.set_job_states("900_1", ["PD:JobHeldUser"])
    assert choose(abc_first(0)) == 0
    assert_matlab_history(fake_slurm, "pendingJobsOfTheUserStartFirst")


@POSIX
def test_slurms_own_check_decides(fake_slurm):
    # the nodes have room, Slurm says the job would start later
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    fake_slurm.set_test_only("cpu_queue", "later")
    assert choose(abc_first(0)) == 1
    assert_matlab_history(fake_slurm, "slurmsOwnCheckDecides")


@POSIX
def test_place_slurm_refuses_is_not_used(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    fake_slurm.set_test_only("gpu8", "error")
    places = abc_first(0)
    with pytest.warns(AlternativesWarning) as record:
        k = S().choose(places, [4, 10], ["", ""])[0]
    assert identifiers(record) == ["slurm_alternatives:testOnly"]
    assert k == 0
    n = fake_slurm.call_count("sbatch")
    assert choose(places) == 0
    assert fake_slurm.call_count("sbatch") == n, "Slurm was asked again about a place it refused."
    assert_matlab_history(fake_slurm, "placeSlurmRefusesIsNotUsed")


@POSIX
def test_main_place_is_not_left_out_when_slurm_refuses_it_once(fake_slurm):
    # (no room there then, but it is the user's place: asked again)
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    fake_slurm.set_test_only("cpu_queue", "error")
    places = abc_first(0)
    assert choose(places) == 1
    fake_slurm.set_test_only("cpu_queue", "now")
    assert choose(places) == 0
    assert_matlab_history(fake_slurm, "mainPlaceIsNotLeftOutWhenSlurmRefusesItOnce")


@POSIX
def test_burst_of_jobs_fills_the_room_in_order(fake_slurm):
    # (one view of the cluster for all jobs of a burst) cpu_queue has room for 2 jobs (4 CPUs
    # of 21418 MB), gpu8 for 3 of them (its free memory); the others wait at cpu_queue (not in
    # the GPU nodes' queue)
    fake_slurm.set_nodes([abc_node("n1", 16), gpu_node("g1", 96, 0)])
    places = abc_first(0)
    S().set_timing(view_age=60)
    assert [choose(places) for _ in range(7)] == [0, 0, 1, 1, 1, 0, 0]
    assert fake_slurm.call_count("scontrol") == 1, "The nodes were read more than once in a view."
    n_test = sum("--test-only" in x for x in fake_slurm.calls("sbatch"))
    assert n_test == 2, "Slurm was asked more than once per place and shape in a view."
    assert_matlab_history(fake_slurm, "burstOfJobsFillsTheRoomInOrder")


@POSIX
def test_cpu_jobs_leave_idle_gpus_their_share(fake_slurm):
    # an 8-GPU node with 6 idle GPUs keeps 16 CPUs and 1/8 of its memory per idle GPU (the
    # node's share), here all its free CPUs
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 32, 2)])
    assert choose(abc_first(None)) == 0
    # with the share of the run's GPU jobs (4 CPUs, 32 GB): gpu8
    assert choose(abc_first([4, 32])) == 1
    # 6 x (4 CPUs, 32 GB) kept: 578992 MB free for 6 jobs of 85672 MB
    S().set_timing(view_age=60)
    places = abc_first([4, 32])
    assert [choose(places) for _ in range(8)] == [1, 1, 1, 1, 1, 1, 0, 0]
    # (CPUs kept) 32 free CPUs and all memory, 4 idle GPUs: 16 CPUs for 4 jobs
    fake_slurm.set_nodes(
        [abc_node("n1", 24), FakeSlurm.node_line("g1", "gpu8", 128, 96, 1031600, 0, 8, 4)]
    )
    S().set_timing(view_age=60)
    assert [choose(places) for _ in range(5)] == [1, 1, 1, 1, 0]
    assert_matlab_history(fake_slurm, "cpuJobsLeaveIdleGPUsTheirShare")


@POSIX
def test_features_and_node_lists_of_a_place_count(fake_slurm):
    # a place's -C and -x: here no node of gpu8 is one of its nodes
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0), gpu_node("g2", 0, 0)])
    for param in (f"{GPU} -C gpu1", f"{GPU} -x g[1-2]", f"{GPU} -x g1,g2.cluster"):
        entry = {"SlurmParam": param, "memPerCPU": 7.6, "gpuHeadroom": 0, "cpuJobs": True}
        assert choose(S().places(main_settings(), entry)) == 0, param
    entry = {"SlurmParam": f"{GPU} -C gpu8 -x g1", "memPerCPU": 7.6, "gpuHeadroom": 0}
    entry["cpuJobs"] = True
    assert choose(S().places(main_settings(), entry)) == 1
    assert_matlab_history(fake_slurm, "featuresAndNodeListsOfAPlaceCount")


@POSIX
def test_gpu_jobs_need_a_free_gpu(fake_slurm):
    fake_slurm.set_nodes([gpu_node("g1", 8, 8), gpu_node("g2", 8, 7, "MIXED", "gpu9")])
    main = main_settings()
    main["SlurmParam"] = "-p gpu8 --qos qos_normal -n1 --mem-per-cpu=8000M --gres=gpu:1"
    main["memPerCPU"] = 7.8
    entry = {"SlurmParam": "-p gpu9 --qos qos_normal -n1 --mem-per-cpu=8000M --gres=gpu:1"}
    places = S().places(main, entry)
    assert S().choose(places, [4, 4], ["", ""])[0] == 1
    assert_matlab_history(fake_slurm, "gpuJobsNeedAFreeGPU")


@POSIX
def test_limit_that_pending_jobs_wait_for_closes_a_place(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    fake_slurm.set_job_states("900_1", ["PD:QOSMaxCpuPerUserLimit"])
    fake_slurm.set_job_attributes("900_1", "cpu_queue", 4, 4 * 21418, 0, "qos_normal")
    assert choose(abc_first(0)) == 1
    assert_matlab_history(fake_slurm, "limitThatPendingJobsWaitForClosesAPlace")


@POSIX
def test_cpu_jobs_wait_outside_gpu_nodes(fake_slurm):
    # main place gpu8 (GPU nodes), nothing has room: the job waits at cpu_queue (which has GPU
    # nodes too, not only)
    fake_slurm.set_nodes(
        [
            abc_node("n1", 24),
            gpu_node("g1", 128, 0),
            FakeSlurm.node_line("t1", "cpu_queue", 20, 20, 515638, 0, 4, 4),
        ]
    )
    assert S().choose(gpu_first(), [10, 4], ["", ""])[0] == 1
    assert_matlab_history(fake_slurm, "cpuJobsWaitOutsideGPUNodes")


@POSIX
def test_pending_job_moves_to_a_place_with_room(fake_slurm):
    # (at gpu8 it keeps its cpu_queue shape)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    places = abc_first(0)
    submit_pending(fake_slurm, places, 1000, 1, 0)
    S().set_timing(min_wait=0, look_interval=0)
    with contextlib.redirect_stdout(io.StringIO()) as out:
        S().move([1000], [1])
    assert updates(fake_slurm) == [
        "update JobId=1000_1 Partition=gpu8 QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    ]
    assert "Job 1000_1 moved from cpu_queue (4 CPUs) to gpu8 (4 CPUs)." in out.getvalue()
    assert "of 1 waiting jobs, 0 have room where they wait, 1 moved, 0 have no room" in (
        out.getvalue()
    )
    attrs = fake_slurm.job_attributes("1000_1")
    assert (attrs["partition"], attrs["cpus"], attrs["mem"]) == ("gpu8", 4, 4 * 21418)
    # not again within 2 minutes of the move
    S().set_timing(min_wait=120)
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 128, 0)])
    n = fake_slurm.call_count("scontrol")
    S().move([1000], [1])
    assert fake_slurm.call_count("scontrol") == n, "A job moved again within 2 minutes."
    assert_matlab_history(fake_slurm, "pendingJobMovesToAPlaceWithRoom")


@POSIX
def test_only_as_many_jobs_move_as_fit(fake_slurm):
    # gpu8 has room for 2 jobs of 4 CPUs of 21418 MB
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 106, 0)])
    places = abc_first(0)
    for a in (1, 2, 3):
        submit_pending(fake_slurm, places, 1000 + a, a, 0)
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1001, 1002, 1003], [1, 2, 3])
    assert updated_jobs(fake_slurm) == ["1001_1", "1002_2"]
    assert_matlab_history(fake_slurm, "onlyAsManyJobsMoveAsFit")


@POSIX
def test_jobs_that_do_not_wait_for_resources_stay(fake_slurm):
    # held, waiting for a start time or dependency, not seen by the scheduler yet, running, or
    # with room where they are (cpu_queue has room for 2: the one not seen yet and the last)
    fake_slurm.set_nodes([abc_node("n1", 16), gpu_node("g1", 0, 0)])
    places = abc_first(0)
    states = ["PD:JobHeldUser", "PD:BeginTime", "PD:Dependency", "PD:None", "R", "PD:Priority"]
    for a, state in enumerate(states, 1):
        submit_pending(fake_slurm, places, 1000 + a, a, 0)
        fake_slurm.set_job_states(f"{1000 + a}_{a}", [state])
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1000 + a for a in range(1, 7)], list(range(1, 7)))
    assert updates(fake_slurm) == []
    assert_matlab_history(fake_slurm, "jobsThatDoNotWaitForResourcesStay")


@POSIX
def test_job_with_room_only_outside_its_place_moves(fake_slurm):
    # (its place is cpu_queue without n1; the room at cpu_queue is on n1)
    fake_slurm.set_nodes([abc_node("n1", 0), abc_node("n2", 24), gpu_node("g1", 0, 0)])
    main = main_settings()
    main["SlurmParam"] = f"{ABC} -x n1"
    entry = {
        "SlurmParam": f"{GPU} -x n1",
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "gpuHeadroom": 0,
        "cpuJobs": True,
    }
    places = S().places(main, entry)
    submit_pending(fake_slurm, places, 1000, 1, 0)
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1000], [1])
    assert len(updates(fake_slurm)) == 1
    assert_matlab_history(fake_slurm, "jobWithRoomOnlyOutsideItsPlaceMoves")


@POSIX
def test_pending_jobs_are_found_by_their_ids(fake_slurm):
    # (the first pending job in squeue's list is held; the second one waits for resources and
    # moves)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    places = abc_first(0)
    submit_pending(fake_slurm, places, 1001, 1, 0)
    fake_slurm.set_job_states("1001_1", ["PD:JobHeldUser"])
    submit_pending(fake_slurm, places, 1002, 2, 0)
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1001, 1002], [1, 2])
    assert updated_jobs(fake_slurm) == ["1002_2"]
    assert_matlab_history(fake_slurm, "pendingJobsAreFoundByTheirIds")


@POSIX
def test_no_slurm_calls_before_a_job_has_waited(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    places = abc_first(0)
    submit_pending(fake_slurm, places, 1000, 1, 0)
    fake_slurm.clear_calls()
    S().move([1000], [1])
    assert fake_slurm.history() == [], "Slurm was asked about a job that just came."
    # jobs not recorded (of other calls or users) are never looked at
    S().set_timing(min_wait=0, look_interval=0)
    S().move([2000], [1])
    assert fake_slurm.history() == []
    assert_matlab_history(fake_slurm, "noSlurmCallsBeforeAJobHasWaited")


@POSIX
def test_cpu_job_leaves_the_queue_of_gpu_nodes(fake_slurm):
    # a CPU-only job waiting at gpu8 (main place) without room there or anywhere: it waits at
    # cpu_queue instead, shaped for cpu_queue
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 128, 0)])
    places = gpu_first()
    S().submitted(1000, 1, places, [10, 4], 0, ["", ""])
    fake_slurm.set_job_states("1000_1", ["PD:Resources"])
    fake_slurm.set_job_attributes("1000_1", "gpu8", 10, 78000, 0, "qos_normal")
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1000], [1])
    assert updates(fake_slurm) == [
        "update JobId=1000_1 Partition=cpu_queue QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    ]
    assert_matlab_history(fake_slurm, "cpuJobLeavesTheQueueOfGPUNodes")


@POSIX
def test_places_with_other_options_do_not_trade(fake_slurm):
    # options a pending job can't change (here --exclusive), and a place without a QOS for jobs
    # of another partition (Slurm keeps the job's QOS, also one the new partition doesn't allow)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    for param in (f"{GPU} --exclusive", "-p gpu8 -n1 --mem-per-cpu=7800M"):
        entry = {
            "SlurmParam": param,
            "memPerCPU": 7.6,
            "maxCPUNum": 32,
            "gpuHeadroom": 0,
            "cpuJobs": True,
        }
        places = S().places(main_settings(), entry)
        submit_pending(fake_slurm, places, 1000, 1, 0)
        S().set_timing(min_wait=0, look_interval=0)
        S().move([1000], [1])
        assert updates(fake_slurm) == [], param
    assert_matlab_history(fake_slurm, "placesWithOtherOptionsDoNotTrade")


@POSIX
def test_job_moved_by_something_else_is_left_to_it(fake_slurm):
    # (recorded at cpu_queue, found at gpu8: e.g. a user's script moved it; gpu9 has room)
    fake_slurm.set_nodes(
        [abc_node("n1", 24), gpu_node("g1", 128, 0), gpu_node("g2", 0, 0, "MIXED", "gpu9")]
    )
    entries = [
        {"SlurmParam": p, "memPerCPU": 7.6, "gpuHeadroom": 0, "cpuJobs": True}
        for p in (GPU, GPU.replace("gpu8", "gpu9"))
    ]
    places = S().places(main_settings(), entries)
    submit_pending(fake_slurm, places, 1000, 1, 0)
    fake_slurm.set_job_attributes("1000_1", "gpu8", 6, 6 * 21418, 0, "qos_normal")
    S().set_timing(min_wait=0, look_interval=0)
    S().move([1000], [1])
    # (then cpu_queue has room: still not moved)
    fake_slurm.set_nodes(
        [abc_node("n1", 0), gpu_node("g1", 128, 0), gpu_node("g2", 0, 0, "MIXED", "gpu9")]
    )
    S().set_timing(min_wait=0)
    S().move([1000], [1])
    assert updates(fake_slurm) == []
    assert_matlab_history(fake_slurm, "jobMovedBySomethingElseIsLeftToIt")


@POSIX
def test_failed_move_leaves_the_job(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    places = abc_first(0)
    submit_pending(fake_slurm, places, 1000, 1, 0)
    fake_slurm.fail_updates()
    S().set_timing(min_wait=0, look_interval=0)
    with contextlib.redirect_stdout(io.StringIO()) as out:
        S().move([1000], [1])
    assert "was not moved" in out.getvalue()
    assert fake_slurm.job_attributes("1000_1")["partition"] == "cpu_queue"
    assert_matlab_history(fake_slurm, "failedMoveLeavesTheJob")


@POSIX
def test_framework_submits_to_the_place_with_room(fake_slurm, tmp_path, virtual_clock):
    # a main place of GPU nodes (for CPU jobs), whose idle GPUs keep their share; cpu_queue has
    # room: the jobs go to cpu_queue with the CPUs of a cpu_queue config (4, not gpu8's 16),
    # and the time limit of its SlurmParam (not the framework's)
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    ins, outs, tasks = make_tasks(tmp_path, 2)
    conf = gpu_first_config()
    conf["SlurmAlternatives"][0]["SlurmParam"] = f"{ABC} -t 2:00:00"
    framework(ins, outs, tasks, tmp_path, write_config(tmp_path, conf))
    jobs = submissions(fake_slurm)
    assert len(jobs) == 2
    assert all(f"--cpus-per-task=4 {ABC} -t 2:00:00 " in x for x in jobs), jobs
    assert not any(" -t 24:00:00 " in x for x in jobs), jobs
    assert_matlab_history(fake_slurm, "frameworkSubmitsToThePlaceWithRoom", tmp_path)


@POSIX
def test_compiled_framework_submits_to_the_place_with_room(fake_slurm, tmp_path, virtual_clock):
    # (MATLAB's compiled mode: its keys and mccMode are taken and ignored, the jobs are those of
    # a framework call without them)
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    ins, outs, tasks = make_tasks(tmp_path, 2)
    conf = gpu_first_config()
    conf.update(
        MCRParam=str(tmp_path / "mcr"),
        MCCMasterStr=str(tmp_path / "run_mccMaster.sh"),
        MCRCacheRoot=str(tmp_path / "cache"),
        BashLaunchStr="echo",
    )
    framework(ins, outs, tasks, tmp_path, write_config(tmp_path, conf), mcc_mode=True)
    jobs = submissions(fake_slurm)
    assert len(jobs) == 2
    assert all(f"--cpus-per-task=4 {ABC} " in x for x in jobs), jobs
    assert all(" -t 24:00:00 " in x for x in jobs), jobs
    assert_matlab_history(fake_slurm, "compiledFrameworkSubmitsToThePlaceWithRoom", tmp_path)


@POSIX
def test_framework_moves_pending_jobs(fake_slurm, tmp_path, virtual_clock):
    # the job of batch 1 waits at cpu_queue (full, gpu8 too); then gpu8 has room (the second
    # look at the nodes): the wait loop moves it there, with its cpu_queue shape, and not again
    # once it fits there
    fake_slurm.set_nodes(
        [abc_node("n1", 24), gpu_node("g1", 128, 0)], [abc_node("n1", 24), gpu_node("g1", 0, 0)]
    )
    fake_slurm.set_job_states("1000_1", ["PD:Priority"] * 6 + ["done"])
    ins, outs, tasks = make_tasks(tmp_path, 1)
    S().set_timing(min_wait=0, look_interval=0, view_age=0)
    config_file = write_config(tmp_path, abc_first_config(0))
    framework(ins, outs, tasks, tmp_path, config_file, query_interval=0)
    jobs = submissions(fake_slurm)
    assert len(jobs) == 1
    assert f"--cpus-per-task=4 {ABC} " in jobs[0]
    assert updates(fake_slurm) == [
        "update JobId=1000_1 Partition=gpu8 QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    ]
    assert_matlab_history(fake_slurm, "frameworkMovesPendingJobs", tmp_path)


def single_job(
    places_conf: dict,
    folder: Path,
    *,
    array_id: int = 3,
    mem_allocate: float = 70,
    nice: float = 0,
    job_id: int = -1,
):
    """generic_single_job_submit_wrapper's steps that use the places, on the fake (J07's
    wrapper calls them so): sinfo, the places of the call's config at the job's first
    submission, the choice, sbatch and the record; returns (job id, its places)"""
    slurm = Slurm()
    check = check_slurm_cluster(str(folder), str(folder), slurm)
    cpus = max(1, min(places_conf["maxCPUNum"], math.ceil(mem_allocate / places_conf["memPerCPU"])))
    sp = alt.single_job_places(
        alternatives=places_conf["SlurmAlternatives"],
        slurm_param=places_conf["SlurmParam"],
        slurm_constraint="",
        mem_per_cpu=places_conf["memPerCPU"],
        max_cpu_num=places_conf["maxCPUNum"],
        min_cpu_num=1,
        gpu_headroom=places_conf.get("gpuHeadroom"),
        requested_cpus=1,
        mem_allocate=mem_allocate,
        main_cpus=cpus,
        job_time_limit=24,
        main_time_option=time_str(24, places_conf["SlurmParam"], TIME_FORMAT_SINGLE),
        nice=nice,
    )
    k, cpus, param, constraint, time_option = sp.choose()
    r = slurm.submit(
        array_id=array_id,
        out=check.job_log_fname,
        err=check.job_log_error_fname,
        cpus_per_task=cpus,
        wrap='date; echo Matlab command:  "disp(1)"',
        slurm_param=param,
        slurm_constraint=constraint,
        time_str=time_option,
    )
    assert r.submitted
    sp.submitted(r.job_id, array_id, k)
    return r.job_id, sp


def single_job_status(job_id: int, array_id: int) -> int | None:
    """check_single_job_status (and the wrapper's look at a job it submitted): pending jobs may
    move, ended ones are forgotten"""
    status = Slurm().job_status(job_id, array_id)
    if status == 0:
        alt.move([job_id], [array_id])
    elif status is not None and status < 0:
        alt.forget([job_id], [array_id])
    return status


@POSIX
def test_single_job_goes_to_the_place_with_room_and_moves(fake_slurm, tmp_path):
    # the pipeline's and decon's single jobs (main place gpu8 without GPU headroom): submitted
    # to cpu_queue (gpu8 full) with a cpu_queue shape, then moved back to gpu8 once it has room
    # and cpu_queue not (keeping the cpu_queue shape), by the wrapper's and
    # check_single_job_status's looks at pending jobs
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 128, 0)])
    conf = gpu_first_config()
    conf["gpuHeadroom"] = 0
    with contextlib.redirect_stdout(io.StringIO()):
        job_id, _ = single_job(conf, tmp_path)
    jobs = submissions(fake_slurm)
    assert len(jobs) == 1
    assert f"--cpus-per-task=4 {ABC} " in jobs[0]
    fake_slurm.set_job_states(f"{job_id}_3", ["PD:Priority"])
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    S().set_timing(min_wait=0, look_interval=0)
    with contextlib.redirect_stdout(io.StringIO()):
        assert single_job_status(job_id, 3) == 0
    assert updates(fake_slurm) == [
        f"update JobId={job_id}_3 Partition=gpu8 QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    ]
    # and check_single_job_status looks at it as well
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 128, 0)])
    S().set_timing(min_wait=0)
    with contextlib.redirect_stdout(io.StringIO()):
        single_job_status(job_id, 3)
    assert len(updates(fake_slurm)) == 2
    assert_matlab_history(fake_slurm, "singleJobGoesToThePlaceWithRoomAndMoves", tmp_path)


@POSIX
def test_config_without_alternatives_is_as_before(fake_slurm, tmp_path, virtual_clock):
    # an empty list: the same jobs as without the key, and no Slurm command but sinfo, sbatch
    # and squeue -j
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    ins, outs, tasks = make_tasks(tmp_path, 2)
    conf = abc_first_config(0)
    conf["SlurmAlternatives"] = []
    without = {k: v for k, v in conf.items() if k != "SlurmAlternatives"}
    configs = [write_config(tmp_path, without, "a.json"), write_config(tmp_path, conf, "b.json")]
    jobs = []
    for config_file in configs:
        fake_slurm.clear_calls()
        framework(ins, outs, tasks, tmp_path, config_file, uuid="u1")
        jobs.append([_CALL_FOLDER.sub("<CALL>", x) for x in fake_slurm.calls("sbatch")])
        cmds = fake_slurm.history()
        assert all(x.startswith(("sinfo", "sbatch", "squeue -j", "scancel")) for x in cmds), cmds
        assert not any("--test-only" in x for x in cmds)
    assert jobs[1] == jobs[0]
    assert_matlab_history(fake_slurm, "configWithoutAlternativesIsAsBefore", tmp_path)


def test_matlab_histories_are_compared_with_the_comment_set_aside():
    # MATLAB tags every submission of a framework call with --comment=pk5d_<uuid12>_<batch>,
    # the cluster loop with its own tag: the comparison sets both aside, and nothing else
    names = [
        "frameworkSubmitsToThePlaceWithRoom",
        "compiledFrameworkSubmitsToThePlaceWithRoom",
        "frameworkMovesPendingJobs",
        "configWithoutAlternativesIsAsBefore",
        "pendingJobOfACallWithNiceMoves",
    ]
    tag = re.compile(r"=pk5d_[0-9A-Za-z]{1,12}_(?=\d+ )")
    for name in names:
        subs = [x for x in expand(MATLAB_HISTORIES[name]) if x.startswith("sbatch --array=")]
        assert subs and all(re.search(rf" --comment{tag.pattern}", x) for x in subs), name
    matlab = expand(MATLAB_HISTORIES["frameworkSubmitsToThePlaceWithRoom"])
    ours = [tag.sub("=pk5d_call_20261112_000000_000_a1b2c3:", x) for x in matlab]
    assert ours != matlab and comment_aside(ours) == comment_aside(matlab)
    untagged = [re.sub(r" --comment=\S+", "", x) for x in matlab]
    assert untagged != matlab and comment_aside(untagged) == comment_aside(matlab)
    assert comment_aside([x.replace("--array=2", "--array=3") for x in ours]) != comment_aside(
        matlab
    )
    assert comment_aside([x.replace("-t 2:00:00", "-t 2:00:01") for x in ours]) != comment_aside(
        matlab
    )


# --- CoordinatorNiceTest/jobsAtAnotherPlaceKeepTheirNice --------------------------------------


@POSIX
def test_jobs_at_another_place_keep_their_nice(fake_slurm, tmp_path, virtual_clock):
    # with SlurmAlternatives a job goes to the first place with room (here the alternative
    # cpu_queue, as gpu8 is full): it keeps its nice there, as a single job and as a framework
    # job (a moved pending job keeps it, as Slurm keeps its nice)
    fake_slurm.set_nodes(
        [
            FakeSlurm.node_line("n1", "cpu_queue", 24, 0, 514048, 0, 0, 0, "MIXED", "c24"),
            FakeSlurm.node_line("g1", "gpu8", 128, 128, 1031600, 128 * 8000, 8, 0, "MIXED"),
        ]
    )
    conf = {
        "SlurmParam": GPU,
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "gpuHeadroom": 0,
        "SlurmAlternatives": [{"SlurmParam": ABC, "memPerCPU": 20.9, "maxCPUNum": 24}],
    }
    with contextlib.redirect_stdout(io.StringIO()):
        single_job(conf, tmp_path, nice=1000)
    ins, outs, tasks = make_tasks(tmp_path, 1)
    framework(ins, outs, tasks, tmp_path, write_config(tmp_path, conf), query_interval=0, nice=1000)
    jobs = submissions(fake_slurm)
    assert len(jobs) == 2
    assert all(ABC in x for x in jobs), "The jobs did not go to the place with room."
    assert all(" --nice=1000" in x for x in jobs), jobs


@POSIX
def test_a_pending_job_with_a_nice_moves_to_an_alternative_and_keeps_its_nice(
    fake_slurm, tmp_path, virtual_clock
):
    # --nice is the call's (the user's rule): a job of a call with a nice waits at the main
    # place, whose options hold the nice; when gpu8 has room, it moves there as a job without a
    # nice does (MATLAB's frameworkMovesPendingJobs), with the nice on each sbatch and an update
    # that leaves it (Slurm keeps a job's nice when it moves); MATLAB before its fix (oracle
    # PetaKit5D_1697257) didn't move it (the oracle case nice_job_moves_to_an_alternative)
    fake_slurm.set_nodes(
        [abc_node("n1", 24), gpu_node("g1", 128, 0)], [abc_node("n1", 24), gpu_node("g1", 0, 0)]
    )
    fake_slurm.set_job_states("1000_1", ["PD:Priority"] * 6 + ["done"])
    ins, outs, tasks = make_tasks(tmp_path, 1)
    S().set_timing(min_wait=0, look_interval=0, view_age=0)
    config_file = write_config(tmp_path, abc_first_config(0))
    framework(ins, outs, tasks, tmp_path, config_file, query_interval=0, nice=5)
    jobs = submissions(fake_slurm)
    assert len(jobs) == 1 and f"--cpus-per-task=4 {ABC} --nice=5 -t 24:00:00 " in jobs[0]
    assert updates(fake_slurm) == [
        "update JobId=1000_1 Partition=gpu8 QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    ]
    expected = [
        x.replace("{ABC} -t", "{ABC} --nice=5 -t").replace("{GPU_ABC} -t", "{GPU_ABC} --nice=5 -t")
        for x in MATLAB_HISTORIES["frameworkMovesPendingJobs"]
    ]
    assert comment_aside(history(fake_slurm, tmp_path)) == comment_aside(expand(expected))


# --- SlurmAlternativesNiceTest (its places: cpu, and gpu, which takes CPU jobs) ---------------

NICE_CPU = "-p cpu --qos normal -n1 --mem-per-cpu=21418M"
NICE_GPU = "-p gpu --qos normal -n1 --mem-per-cpu=7800M"


def nice_cpu_node(cpu_alloc: int) -> str:
    return FakeSlurm.node_line(
        "n1", "cpu", 24, cpu_alloc, 514048, cpu_alloc * 21418, 0, 0, "MIXED", "cpu24"
    )


def nice_gpu_node(cpu_alloc: int) -> str:
    return FakeSlurm.node_line(
        "g1", "gpu", 128, cpu_alloc, 1031600, cpu_alloc * 8000, 8, 0, "MIXED", "gpu8"
    )


def nice_conf(alt_param: str) -> dict:
    """a cluster config of cpu with gpu (alt_param) as its alternative"""
    alternative = {
        "SlurmParam": alt_param,
        "memPerCPU": 7.6,
        "maxCPUNum": 32,
        "gpuHeadroom": 0,
        "cpuJobs": True,
    }
    return {
        "SlurmParam": NICE_CPU,
        "memPerCPU": 20.9,
        "maxCPUNum": 24,
        "SlurmAlternatives": [alternative],
    }


def nice_call(folder: Path, config_file: str, nice: float) -> None:
    """a framework call of one task of 70 GB with nice (its job never runs: the call gives up
    after its try)"""
    ins, outs, tasks = make_tasks(folder, 1)
    framework(ins, outs, tasks, folder, config_file, query_interval=0, nice=nice)


@POSIX
def test_pending_job_of_a_call_with_nice_moves(fake_slurm, tmp_path, virtual_clock):
    # the job of batch 1 waits at cpu (full, gpu too); then gpu has room (the second look at the
    # nodes): the wait loop moves it there, and it keeps its nice
    fake_slurm.set_nodes(
        [nice_cpu_node(24), nice_gpu_node(128)], [nice_cpu_node(24), nice_gpu_node(0)]
    )
    fake_slurm.set_job_states("1000_1", ["PD:Priority"] * 6 + ["done"])
    S().set_timing(min_wait=0, look_interval=0, view_age=0)
    nice_call(tmp_path, write_config(tmp_path, nice_conf(NICE_GPU)), 1000)
    jobs = submissions(fake_slurm)
    assert len(jobs) == 1 and f"{NICE_CPU} --nice=1000" in jobs[0]
    assert updates(fake_slurm) == [
        "update JobId=1000_1 Partition=gpu QOS=normal NumCPUs=4 CPUsPerTask=4 MinMemoryCPU=21418"
    ], "The pending job of a call with nice did not move."
    assert_matlab_history(fake_slurm, "pendingJobOfACallWithNiceMoves", tmp_path)


@POSIX
def test_jobs_get_the_calls_nice_at_every_place(fake_slurm, tmp_path, virtual_clock):
    # a place's own --nice gave way to the call's only where it had none: here the jobs go to
    # gpu (cpu is full), whose SlurmParam has --nice=5, and get the call's 1000 there, as
    # framework jobs and as single jobs
    fake_slurm.set_nodes([nice_cpu_node(24), nice_gpu_node(0)])
    S().set_timing(view_age=0)
    conf = nice_conf(f"{NICE_GPU} --nice=5")
    nice_call(tmp_path, write_config(tmp_path, conf), 1000)
    with contextlib.redirect_stdout(io.StringIO()):
        single_job(conf, tmp_path, nice=1000)
    jobs = submissions(fake_slurm)
    assert len(jobs) == 2
    for x in jobs:
        assert "-p gpu " in x and re.findall(r"--nice=(\d+)", x) == ["1000"], x
    assert_matlab_history(fake_slurm, "jobsGetTheCallsNiceAtEveryPlace", tmp_path)


# --- the parts, and what differs from MATLAB -----------------------------------------------------


def test_options_parse_as_matlab_reads_them():
    o = alt.parse_options("--partition gpu8 --qos=q -Aacct --mem-per-cpu 7800 --gpus=a100:3")
    assert (o.partition, o.qos, o.account, o.mem_per_cpu, o.gpus) == ("gpu8", "q", "acct", 7800, 3)
    assert o.mem_token == "--mem-per-cpu=7800"
    o = alt.parse_options("--gres=gpu:2,shard:1 -n 2 --nice=5 -c 4 --exclusive")
    assert (o.gres, o.gpus, o.ntasks) == ("gpu:2,shard:1", 2, 2)
    # (unknown options and their values are other tokens, kept in order; --nice is apart, the
    # last one as given)
    assert o.other == ("-c", "4", "--exclusive") and o.nice == "--nice=5"
    assert alt.parse_options("-p a --nice=2 --nice -x n1").nice == "--nice"
    assert alt.parse_options("-p a --nicex").other == ("--nicex",)
    assert alt.parse_options("--gres=gpu:a100").gpus == 1
    assert alt.parse_options("--gres=gpu").gpus == 1
    assert alt.parse_options("--mem=1T --mem-per-cpu=512K").mem_token == "--mem-per-cpu=512K"
    o = alt.parse_options("--mem=1T --mem-per-cpu=512K")
    assert (o.mem, o.mem_per_cpu) == (1024**2, 0.5)
    assert alt.parse_options("--mem=abc").mem is None
    assert math.isnan(alt.parse_options("-n x").ntasks)
    # (an option at the end without its value)
    assert alt.parse_options("-p").partition == ""
    assert alt.parse_options("-p\tgpu8\t--qos  q").qos == "q"


def test_memory_options_are_replaced():
    w = alt.with_memory
    assert w(ABC, "--mem=1G") == "-p cpu_queue --qos qos_normal -n1 --mem=1G"
    assert w("-p a --mem 5G -n1", "--mem-per-cpu=1M") == "-p a --mem-per-cpu=1M -n1"
    assert w("-p a", "--mem=1G") == "-p a --mem=1G"
    assert w("--mem=1 -p a --mem-per-cpu=2", "--mem=3") == "--mem=3 -p a"
    assert w("", "--mem=3") == " --mem=3"


def test_host_lists_expand():
    e = alt.expand_hostlist
    assert e("nd0[01-03,07].x") == ["nd001.x", "nd002.x", "nd003.x", "nd007.x"]
    assert e("g[1-2]") == ["g1", "g2"]
    assert e("g1,g2.cluster") == ["g1", "g2.cluster"]
    assert e("n[08-10]") == ["n08", "n09", "n10"]
    assert e("n[3-1]") == []


def test_pending_job_ids_give_the_keys_of_their_tasks():
    keys = alt._pending_keys
    assert keys("12_[3]") == (["12_3"], 1)
    assert keys("12_3") == (["12_3"], 1)
    # (a job that is no array is never one of ours)
    assert keys("12") == ([], 1)
    assert keys("12_[1-4,7%2]") == (["12_1", "12_2", "12_3", "12_4", "12_7"], 5)
    k, n = keys("12_[0-1500]")
    assert n == 1000 and k[-1] == "12_999"


def test_move_fields_follow_slurms_rules():
    places = S().places(
        main_settings(),
        [
            {"SlurmParam": GPU, "memPerCPU": 7.6},
            {"SlurmParam": "-p gpu8 -n1 --mem-per-cpu=7800M", "memPerCPU": 7.6},
            {"SlurmParam": f"-A acct2 {ABC}", "memPerCPU": 20.9},
            {"SlurmParam": f"{ABC} -C c48 -x n1 -w n[2-3] -J pk", "memPerCPU": 20.9},
            {"SlurmParam": "-p cpu_queue --qos qos_normal -n1", "memPerCPU": 20.9},
        ],
    )
    cpus = [4, 10, 10, 4, 4, 4]
    main = alt._main_shape(places, cpus)

    def fields(a, b):
        cur = alt._place_shape(places, cpus, a, main)
        shape = alt._place_shape(places, cpus, b, cur)
        return alt.move_fields(places[a], places[b], cur, shape)

    assert fields(0, 1) == {
        "Partition": "gpu8",
        "QOS": "qos_normal",
        "NumCPUs": "4",
        "CPUsPerTask": "4",
        "MinMemoryCPU": "21418",
    }
    # a place without a QOS takes no jobs of another partition
    assert fields(0, 2) == {}
    # another account: the QOS and the account are set
    assert list(fields(0, 3)) == ["Partition", "QOS", "Account", "NumCPUs", "CPUsPerTask"] + [
        "MinMemoryCPU"
    ]
    # features, node lists and the job name are set, and cleared back ('wrap' for the name)
    assert fields(0, 4)["Features"] == "c48" and fields(0, 4)["ReqNodeList"] == "n[2-3]"
    back = fields(4, 0)
    assert (back["Features"], back["ExcNodeList"], back["ReqNodeList"], back["JobName"]) == (
        "",
        "",
        "",
        "wrap",
    )
    # a place without a memory option: its jobs get the memory of the place they move to, but
    # none comes back to it (the partition's default memory can't be set again)
    assert fields(5, 0)["MinMemoryCPU"] == "21418"
    assert fields(0, 5) == {}


def test_cluster_places_shape_each_place_as_a_config_of_it_alone():
    settings = {
        "slurm_param": GPU,
        "slurm_constraint": "",
        "mem_per_cpu": 7.6,
        "max_cpu_num": 32,
        "min_cpu_num": 1,
        "gpu_headroom": None,
        "mem_allocate": 70,
        "cpus_per_task": 16,
        "job_time_limit": 24,
        "nice": 0,
        "gpu_job": False,
        "slurm_alternatives": [
            {"SlurmParam": ABC, "memPerCPU": 20.9, "maxCPUNum": 24},
            {"SlurmParam": f"{ABC} -t 2:00:00", "memPerCPU": 20.9, "maxCPUNum": 24},
        ],
    }
    call = Call(["/in"], ["/out"], [{"bash": "true"}], settings)
    shape = JobShape(16, 1, 1, False, False)
    places = alt.cluster_places(call, shape, " -t 24:00:00 ", requested_cpus=1)
    assert [p.cpus_per_task for p in places.places] == [16, 4, 4]
    assert places.time_options == [" -t 24:00:00 ", " -t 24:00:00 ", ""]
    # (without the CPUs asked for, the main place's: never fewer CPUs than MATLAB gives)
    fallback = alt.cluster_places(call, shape, " -t 24:00:00 ")
    assert [p.cpus_per_task for p in fallback.places] == [16, 16, 16]
    # the in-job slots' time format
    slots = alt.cluster_places(call, JobShape(16, 2, 2, True, True), " -t 24:00:00 ", 1)
    assert slots.time_options[1] == " -t 24:00:00 "
    call.settings["job_time_limit"] = 1.5
    slots = alt.cluster_places(call, JobShape(16, 2, 2, True, True), " -t 01:30:00 ", 1)
    assert slots.time_options == [" -t 01:30:00 ", " -t 01:30:00 ", ""]
    assert time_str(1.5, ABC) == " -t 1:30:00 "
    # jobs of whole nodes: the main place only
    call.settings["whole_node_job"] = True
    assert alt.cluster_places(call, shape, "") is None
    call.settings["whole_node_job"] = False
    # one place: as without alternatives
    call.settings["slurm_alternatives"] = []
    assert alt.cluster_places(call, shape, "") is None


def test_cpus_of_a_try_at_each_place():
    p = [
        alt.Place("a", "", 20.9, 24, 1, None, True, alt.parse_options("-p a"), 4),
        alt.Place("b", "", 7.6, 32, 1, None, True, alt.parse_options("-p b"), 16),
    ]
    session = Session()
    # CPU jobs: CPUs x (try + 1) within maxCPUNum, GPU jobs keep theirs
    assert alt.ClusterPlaces(session, p, ["", ""], False, False).own_cpus(2) == [12, 32]
    assert alt.ClusterPlaces(session, p, ["", ""], False, True).own_cpus(2) == [4, 16]
    # in-job slots: the other way round
    assert alt.ClusterPlaces(session, p, ["", ""], True, False).own_cpus(2) == [4, 16]
    assert alt.ClusterPlaces(session, p, ["", ""], True, True).own_cpus(1) == [8, 32]


def test_the_calls_nice_goes_to_every_place_and_keeps_no_jobs_apart():
    # --nice is the call's, not a place's: every place's jobs get it in place of a --nice of
    # the place's own (with_nice), and the options keep --nice apart from the other tokens, so
    # the main place (whose SlurmParam holds the call's nice) and the alternatives trade jobs, as
    # do places with a --nice of their own; a move leaves the nice to Slurm, which keeps it
    # (MATLAB before its fix moved no job with a nice between the main place and an alternative)
    settings = {
        "slurm_param": f"{ABC} --nice=5",
        "slurm_constraint": "",
        "mem_per_cpu": 20.9,
        "max_cpu_num": 24,
        "min_cpu_num": 1,
        "gpu_headroom": None,
        "mem_allocate": 70,
        "job_time_limit": 24,
        "nice": 5,
        "gpu_job": False,
        "slurm_alternatives": [
            {"SlurmParam": GPU, "memPerCPU": 7.6, "maxCPUNum": 32},
            {"SlurmParam": GPU.replace("gpu8", "gpu9"), "memPerCPU": 7.6, "maxCPUNum": 32},
            {"SlurmParam": f"{ABC} --nice=2", "memPerCPU": 20.9},
        ],
    }
    call = Call(["/in"], ["/out"], [{"bash": "true"}], settings)
    places = alt.cluster_places(call, JobShape(4, 1, 1, False, False), " -t 24:00:00 ", 1).places
    assert [p.slurm_param for p in places] == [
        f"{ABC} --nice=5",
        f"{GPU} --nice=5",
        f"{GPU.replace('gpu8', 'gpu9')} --nice=5",
        f"{ABC} --nice=5",
    ]
    own = [4, 16, 16, 4]
    shape = alt._main_shape(places, own)
    to_gpu = alt._place_shape(places, own, 1, shape)
    assert alt.move_fields(places[0], places[1], shape, to_gpu) == {
        "Partition": "gpu8",
        "QOS": "qos_normal",
        "NumCPUs": "4",
        "CPUsPerTask": "4",
        "MinMemoryCPU": "21418",
    }
    back = alt._place_shape(places, own, 0, to_gpu)
    assert alt.move_fields(places[1], places[0], to_gpu, back)["Partition"] == "cpu_queue"
    assert alt.move_fields(places[1], places[2], to_gpu, to_gpu)["Partition"] == "gpu9"
    assert alt.move_fields(places[3], places[0], shape, shape)["Partition"] == "cpu_queue"
    # other options still keep places apart
    exclusive = replace(places[1], opts=alt.parse_options(f"{GPU} --exclusive --nice=5"))
    assert alt.move_fields(places[0], exclusive, shape, to_gpu) == {}
    # without a nice, each place's SlurmParam is as configured
    call.settings.update(slurm_param=ABC, nice=0)
    places = alt.cluster_places(call, JobShape(4, 1, 1, False, False), " -t 24:00:00 ", 1).places
    assert [p.slurm_param for p in places][::3] == [ABC, f"{ABC} --nice=2"]

    # single jobs: the call's nice in place of the chosen place's own, at the main place too
    class Chooses:
        def __init__(self, k: int, param: str) -> None:
            self.k, self.param = k, param

        def choose(self, places, own_cpus, time_options):
            return self.k, own_cpus[self.k], self.param

    for k, param in enumerate([ABC, GPU, GPU, f"{ABC} --nice=2"]):
        sp = alt.SingleJobPlaces(Chooses(k, param), places, own, [""] * 4, 5)
        assert sp.choose()[2] == f"{param.replace(' --nice=2', '')} --nice=5"
        assert alt.SingleJobPlaces(Chooses(k, param), places, own, [""] * 4, 0).choose()[2] == param
    # with_nice: every --nice token goes, with the whitespace before it, then the call's comes
    # (MATLAB's regexprep and sprintf; the oracle case with_nice has more)
    assert alt.with_nice("--nice=2 -p a\t--nice --comment=--nice=3", 2.5) == (
        " -p a --comment=--nice=3 --nice=3"
    )


def test_single_job_places_are_not_rounded():
    sp = alt.single_job_places(
        alternatives=[
            {"SlurmParam": ABC, "memPerCPU": 20.9, "maxCPUNum": 24},
            {"SlurmParam": "-p x -n1 --mem-per-cpu=4G", "maxCPUNum": 8, "minCPUNum": 2},
        ],
        slurm_param=GPU,
        slurm_constraint="",
        mem_per_cpu=7.6,
        max_cpu_num=32,
        min_cpu_num=1,
        gpu_headroom=None,
        requested_cpus=5,
        mem_allocate=70,
        main_cpus=10,
        job_time_limit=1.5,
        main_time_option=" -t 1:30:00 ",
        session=Session(),
    )
    # max(cpusPerTask, ceil(memAllocate / memPerCPU)) within minCPUNum and maxCPUNum
    assert sp.own_cpus == [10, 5, 8]
    assert sp.time_options == [" -t 1:30:00 ", " -t 1:30:00 ", " -t 1:30:00 "]
    no_mem = alt.single_job_places(
        alternatives=[{"SlurmParam": ABC, "memPerCPU": 20.9}],
        slurm_param=GPU,
        slurm_constraint="",
        mem_per_cpu=7.6,
        max_cpu_num=32,
        min_cpu_num=1,
        gpu_headroom=None,
        requested_cpus=1,
        mem_allocate=None,
        main_cpus=1,
        job_time_limit=None,
        main_time_option="",
        session=Session(),
    )
    assert no_mem.own_cpus == [1, 1] and no_mem.time_options == ["", ""]
    assert alt.single_job_places(alternatives=None, **_SINGLE_REST) is None


_SINGLE_REST = {
    "slurm_param": GPU,
    "slurm_constraint": "",
    "mem_per_cpu": 7.6,
    "max_cpu_num": 32,
    "min_cpu_num": 1,
    "gpu_headroom": None,
    "requested_cpus": 1,
    "mem_allocate": 70,
    "main_cpus": 10,
    "job_time_limit": 24,
    "main_time_option": "",
}


def test_entries_that_matlab_refuses():
    for bad in (["-p x"], [{"SlurmParam": 5}], "-p x", [{"SlurmParam": "-p x", "memPerCPU": "7"}]):
        with pytest.raises(AlternativesError) as err:
            S().places(main_settings(), bad)
        assert err.value.identifier == "slurm_alternatives:entry"
    # cpuJobs is true only for true or 1
    places = S().places(
        main_settings(),
        [
            {"SlurmParam": GPU, "memPerCPU": 7.6, "cpuJobs": v}
            for v in (True, 1, 1.0, [1], "true", 0, None, 2)
        ],
    )
    assert [p.cpu_jobs for p in places[1:]] == [True, True, True, True, False, False, False, False]


@POSIX
def test_an_entry_a_shell_cannot_split_is_not_used(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    entry = {"SlurmParam": f'{GPU} -J "pk', "memPerCPU": 7.6, "gpuHeadroom": 0, "cpuJobs": True}
    with pytest.warns(AlternativesWarning) as record:
        places = S().places(main_settings(), entry)
    assert identifiers(record) == ["slurm_alternatives:testOnly"]
    # (not as the place with room, nor as the place to wait in)
    assert S().choose(places, [4, 10], ["", ""]) == (0, 4, ABC)
    assert not any("--test-only" in x for x in fake_slurm.calls("sbatch"))
    gpu_main = {**main_settings(), "SlurmParam": GPU, "memPerCPU": 7.6, "gpuHeadroom": 0}
    with pytest.warns(AlternativesWarning):
        places = S().places(gpu_main, {"SlurmParam": f'{ABC} -J "x', "memPerCPU": 20.9})
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 128, 0)])
    assert S().choose(places, [10, 4], ["", ""])[0] == 0


@POSIX
def test_a_view_slurm_does_not_answer_sends_jobs_to_the_main_place(fake_slurm, capsys):
    # (MATLAB read a failed scontrol as no nodes, a failed squeue as no pending jobs: CPU jobs
    # could go where GPU jobs of the user wait)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    fake_slurm.set_job_states("900_1", ["PD:JobHeldUser"])
    fake_slurm.set_job_attributes("900_1", "gpu8", 4, 32000, 1, "qos_normal")
    places = abc_first(0)
    for cmd in ("squeue", "scontrol"):
        fake_slurm.fail(cmd)
        assert choose(places) == 0
        assert "Slurm did not answer" in capsys.readouterr().err
        fake_slurm.recover(cmd)
    assert not any("--test-only" in x for x in fake_slurm.calls("sbatch"))
    # no job moves on it either; once Slurm answers, the job moves (to cpu_queue: gpu8 has a
    # pending GPU job)
    submit_pending(fake_slurm, places, 1000, 1, 1)
    S().set_timing(min_wait=0, look_interval=0)
    fake_slurm.fail("squeue")
    S().move([1000], [1])
    assert updates(fake_slurm) == []
    fake_slurm.recover("squeue")
    with contextlib.redirect_stdout(io.StringIO()):
        S().move([1000], [1])
    assert updated_jobs(fake_slurm) == ["1000_1"]


@POSIX
def test_a_view_that_times_out_is_unknown(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    session = Session(slurm=Slurm(query_timeout=1))
    places = session.places(main_settings(), {"SlurmParam": GPU, "memPerCPU": 7.6, "cpuJobs": 1})
    session.set_timing(view_age=20)
    fake_slurm.hang("scontrol", 30)
    t0 = time.monotonic()
    with contextlib.redirect_stderr(io.StringIO()):
        assert session.choose(places, [4, 10], ["", ""])[0] == 0
        # (the unknown view is kept for viewAge: Slurm isn't asked again at once)
        assert session.choose(places, [4, 10], ["", ""])[0] == 0
    assert time.monotonic() - t0 < 15
    assert fake_slurm.call_count("scontrol") == 1


@POSIX
def test_test_only_errors_of_the_moment_leave_the_place(fake_slurm):
    import warnings

    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    places = abc_first(0)
    # a busy controller (MATLAB's words), and sbatch without an answer: the place stays
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        fake_slurm.fail("sbatch", "sbatch: error: Socket timed out on send/recv operation")
        assert choose(places) == 0
        fake_slurm.recover("sbatch")
        session = Session(slurm=Slurm(query_timeout=1))
        session.set_timing(view_age=0)
        fake_slurm.hang("sbatch", 30)
        assert session.choose(places, [4, 10], ["", ""])[0] == 0
        fake_slurm.recover("sbatch")
    assert record == []
    assert not S().disabled and not session.disabled
    with contextlib.redirect_stdout(io.StringIO()):
        assert choose(places) == 1


@POSIX
def test_foreign_pending_jobs_are_never_touched(fake_slurm, tmp_path, virtual_clock):
    """jobs the session didn't submit (another call's, another session's, the user's own
    script's) are never updated or cancelled: not when the caller names them, not through the
    framework, also when they wait for resources where an alternative has room"""
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0), gpu_node("g2", 0, 0)])
    foreign = {
        "2001_1": ("cpu_queue", 4, 4 * 21418, 0, "PD:Priority"),
        "2002_1": ("cpu_queue", 2, 2 * 21418, 0, "PD:Resources"),
        "2003_2": ("gpu8", 4, 4 * 8000, 0, "PD:Priority"),
        "2004_1": ("cpu_queue", 4, 4 * 21418, 0, "PD:QOSMaxCpuPerUserLimit"),
        "2005_3": ("cpu_queue", 1, 1000, 0, "PD:JobHeldUser"),
    }
    for key, (part, cpus, mem, gpus, state) in foreign.items():
        fake_slurm.set_job_states(key, [state])
        fake_slurm.set_job_attributes(key, part, cpus, mem, gpus, "qos_normal")
    before = {k: fake_slurm.job_attributes(k) for k in foreign}
    places = abc_first(0)
    submit_pending(fake_slurm, places, 1000, 1, 0)
    S().set_timing(min_wait=0, look_interval=0)
    ids = [1000] + [int(k.split("_")[0]) for k in foreign]
    arrays = [1] + [int(k.split("_")[1]) for k in foreign]
    with contextlib.redirect_stdout(io.StringIO()):
        S().move(ids, arrays)
        S().move(ids[1:], arrays[1:])
    assert updated_jobs(fake_slurm) == ["1000_1"]
    # through the framework with alternatives (its jobs submitted to cpu_queue while gpu8 is
    # full, then gpu8 has room): its own jobs only
    alt.reset()
    S().set_timing(view_age=0, min_wait=0, look_interval=0)
    fake_slurm.clear_calls()
    full = [abc_node("n1", 24), gpu_node("g1", 128, 0), gpu_node("g2", 128, 0)]
    fake_slurm.set_nodes(full, full, full, [abc_node("n1", 24), gpu_node("g1", 0, 0)])
    ins, outs, tasks = make_tasks(tmp_path, 3)
    for b in range(1, 4):
        fake_slurm.set_job_states(f"{999 + b}_{b}", ["PD:Priority"] * 3 + ["done"])
    framework(ins, outs, tasks, tmp_path, write_config(tmp_path, abc_first_config(0)))
    assert all(f"--cpus-per-task=4 {ABC} " in x for x in submissions(fake_slurm))
    assert updated_jobs(fake_slurm) == ["1000_1", "1001_2", "1002_3"]
    own = {"1000", "1001", "1002"}
    for line in fake_slurm.history():
        if line.startswith("scontrol update"):
            assert line.split()[2].split("=")[1].split("_")[0] in own, line
        elif line.startswith("scancel"):
            assert all(x.split("_")[0] in own for x in line.split()[1:]), line
    assert {k: fake_slurm.job_attributes(k) for k in foreign} == before


def test_cluster_places_bind_each_record_to_its_choice(monkeypatch):
    # (a choice and its job: by batch where the loop passes it, else the last choice, and a job
    # whose choice isn't known is not recorded: it never moves)
    p = [
        alt.Place("-p a", "", 20.9, 24, 1, None, True, alt.parse_options("-p a"), 4),
        alt.Place("-p b", "", 7.6, 32, 1, None, True, alt.parse_options("-p b"), 16),
    ]
    session = Session()
    picks = iter([1, 0, 1])
    monkeypatch.setattr(session, "choose", lambda places, own, t: (k := next(picks), own[k], "x"))
    places = alt.ClusterPlaces(session, p, ["", ""], False, False)
    places.choose(4, 0, 1)
    places.choose(4, 0, 2)
    places.submitted(1002, 2, 0, 2)
    places.submitted(1001, 1, 0, 1)
    assert {k: e.k for k, e in session.registry.items()} == {"1002_2": 0, "1001_1": 1}
    places.choose(4, 0)
    places.submitted(1003, 3, 0)
    places.submitted(1004, 4, 0)
    assert session.registry["1003_3"].k == 1 and "1004_4" not in session.registry
    places.forget(["1001_1", "1002_2", "x"])
    assert set(session.registry) == {"1003_3"}


@POSIX
def test_view_is_kept_for_its_age_and_jobs_looked_at_once_a_minute(fake_slurm):
    vc = VirtualClock(time.time())
    session = Session(clock=vc.clock)
    places = session.places(main_settings(), {"SlurmParam": GPU, "memPerCPU": 7.6})
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    for _ in range(5):
        session.choose(places, [4, 10], ["", ""])
        vc.sleep(3.9)
    assert fake_slurm.call_count("scontrol") == 1
    vc.sleep(1)
    session.choose(places, [4, 10], ["", ""])
    assert fake_slurm.call_count("scontrol") == 2
    # a job is looked at after 2 minutes in its place, then at most once a minute
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 128, 0)])
    session.submitted(1000, 1, places, [4, 10], 0, ["", ""])
    fake_slurm.set_job_states("1000_1", ["PD:Priority"])
    fake_slurm.set_job_attributes("1000_1", "cpu_queue", 4, 4 * 21418, 0, "qos_normal")
    looks, views = [], []
    with contextlib.redirect_stdout(io.StringIO()) as out:
        for _ in range(100):
            session.move([1000], [1])
            looks.append(out.getvalue().count("waiting jobs"))
            views.append(fake_slurm.call_count("scontrol"))
            vc.sleep(3)
    # (calls every 3.01 s: looks at 120.4 s, 180.6 s and 240.8 s; no Slurm call between)
    assert [i for i in range(1, 100) if looks[i] > looks[i - 1]] == [40, 60, 80]
    assert looks[0] == 0 and views[0] == 2
    assert [i for i in range(1, 100) if views[i] > views[i - 1]] == [40, 60, 80]


@POSIX
def test_a_moved_job_waits_2_minutes_in_its_new_place(fake_slurm):
    vc = VirtualClock(time.time())
    session = Session(clock=vc.clock)
    entry = {"SlurmParam": GPU, "memPerCPU": 7.6, "cpuJobs": 1, "gpuHeadroom": 0}
    places = session.places(main_settings(), entry)
    session.submitted(1000, 1, places, [4, 10], 0, ["", ""])
    fake_slurm.set_job_states("1000_1", ["PD:Priority"])
    fake_slurm.set_job_attributes("1000_1", "cpu_queue", 4, 4 * 21418, 0, "qos_normal")
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    with contextlib.redirect_stdout(io.StringIO()):
        vc.sleep(130)
        session.move([1000], [1])
        assert updated_jobs(fake_slurm) == ["1000_1"]
        # 70 s after the move (200 s after its submission) it stays, without a Slurm call
        fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 128, 0)])
        n = len(fake_slurm.history())
        vc.sleep(70)
        session.move([1000], [1])
        assert len(fake_slurm.history()) == n
        vc.sleep(60)
        session.move([1000], [1])
    assert updated_jobs(fake_slurm) == ["1000_1", "1000_1"]
    assert updates(fake_slurm)[1].startswith("update JobId=1000_1 Partition=cpu_queue")


@POSIX
def test_slurm_calls_with_places_at_most_one_view_per_20_s(fake_slurm, tmp_path, monkeypatch):
    # the benchmark gate: a framework call with places reads the nodes and the pending jobs at
    # most once per 20 s (MATLAB's view), however many jobs it submits and looks at
    vc = VirtualClock(time.time())
    monkeypatch.setattr(cluster, "_clock", vc.clock)
    monkeypatch.setattr(cluster, "_sleep", vc.sleep)
    alt.reset()
    S().clock = vc.clock
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 128, 0)])
    n = 30
    ins, outs, tasks = make_tasks(tmp_path, n)
    for b in range(1, n + 1):
        fake_slurm.set_job_states(f"{999 + b}_{b}", ["PD:Priority"] * 60 + ["R"], [outs[b - 1]])
    t0 = vc.now
    framework(ins, outs, tasks, tmp_path, write_config(tmp_path, abc_first_config(0)))
    elapsed = vc.now - t0
    views = fake_slurm.calls("scontrol").count("show node -o")
    pending = sum(x.startswith("--me -t PD") for x in fake_slurm.calls("squeue"))
    print(f"{n} jobs, {elapsed:.0f} s: {views} node views, {pending} pending-job queries")
    assert elapsed > 150 and len(submissions(fake_slurm)) == n
    assert 1 <= views <= 1 + elapsed / 20 and pending == views, (views, pending, elapsed)
    # (moves begin once jobs waited 2 minutes)
    assert views >= 2


@POSIX
def test_options_with_shell_characters_reach_slurm_as_written(fake_slurm):
    # (MATLAB ran its commands through a shell: a bare | of -C a|b made a pipe, and the place
    # was refused; quoted values reached Slurm without their quotes, as they do here)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    entry = {"SlurmParam": f"{GPU} -C c24|gpu8", "memPerCPU": 7.6, "gpuHeadroom": 0}
    entry["cpuJobs"] = True
    places = S().places(main_settings(), entry)
    assert places[1].opts.constraint == "c24|gpu8"
    with contextlib.redirect_stdout(io.StringIO()):
        assert choose(places) == 1
    assert fake_slurm.calls("sbatch")[-1].endswith(" -C c24|gpu8 --wrap=true")
    # quoted: the constraint keeps its quotes in the options (as MATLAB's), Slurm gets none
    entry["SlurmParam"] = f"{GPU} -C 'c24|gpu8'"
    places = S().places(main_settings(), entry)
    assert places[1].opts.constraint == "'c24|gpu8'"
    submit_pending(fake_slurm, places, 1000, 1, 0)
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    S().set_timing(min_wait=0, look_interval=0)
    with contextlib.redirect_stdout(io.StringIO()):
        S().move([1000], [1])
    # (no node has the feature 'c24|gpu8' with its quotes: the job waits where it is)
    assert updates(fake_slurm) == []
    alt.reset()
    S().set_timing(view_age=0, min_wait=0, look_interval=0)
    gpu_main = {**main_settings(), "SlurmParam": GPU, "memPerCPU": 7.6, "gpuHeadroom": 0}
    entry = {"SlurmParam": f"{ABC} -C 'c24|c48'", "memPerCPU": 20.9}
    places = S().places(gpu_main, entry)
    S().submitted(1001, 1, places, [10, 4], 0, ["", ""])
    fake_slurm.set_job_states("1001_1", ["PD:Resources"])
    fake_slurm.set_job_attributes("1001_1", "gpu8", 10, 78000, 0, "qos_normal")
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 128, 0)])
    with contextlib.redirect_stdout(io.StringIO()):
        S().move([1001], [1])
    assert updates(fake_slurm) == [
        "update JobId=1001_1 Partition=cpu_queue QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418 Features=c24|c48"
    ]


@POSIX
def test_room_counts_tasks_memory_spec_and_limits_of_any_qos(fake_slurm):
    # a place's -n2 doubles a job's CPUs and memory: 8 CPUs and 171344 MB don't fit the
    # 114048 MB left (with -n1 they would), so Slurm isn't asked
    node = FakeSlurm.node_line("n1", "cpu_queue", 24, 0, 514048, 400000, 0, 0, "MIXED", "c24")
    fake_slurm.set_nodes([node, gpu_node("g1", 128, 0)])
    main = {**main_settings(), "SlurmParam": ABC.replace("-n1", "-n2")}
    places = S().places(main, {"SlurmParam": GPU, "memPerCPU": 7.6, "cpuJobs": True})
    assert choose(places) == 0
    assert not any("--test-only" in x for x in fake_slurm.calls("sbatch"))
    # memory the node keeps for itself (MemSpecLimit) is no node memory: 1031600 - 2048 MB
    spec = gpu_node("g1", 0, 0) + " MemSpecLimit=2048"
    fake_slurm.set_nodes([abc_node("n1", 0), spec])
    places = S().places(main_settings(), {"SlurmParam": "-p gpu8 --qos qos_normal -n1"})
    assert places[1].mem_per_cpu == (1031600 - 2048) / 128 / 1024
    # a QOS or association limit that a pending job waits for closes a place without a QOS,
    # whatever the job's QOS
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    fake_slurm.set_job_states("900_1", ["PD:QOSMaxCpuPerUserLimit"])
    fake_slurm.set_job_attributes("900_1", "cpu_queue", 4, 4 * 21418, 0, "qos_other")
    main = {**main_settings(), "SlurmParam": "-p cpu_queue -n1 --mem-per-cpu=21418M"}
    entry = {"SlurmParam": GPU, "memPerCPU": 7.6, "gpuHeadroom": 0, "cpuJobs": True}
    with contextlib.redirect_stdout(io.StringIO()):
        assert choose(S().places(main, entry)) == 1
    # (and one of another QOS doesn't close a place with a QOS)
    with contextlib.redirect_stdout(io.StringIO()):
        assert choose(abc_first(0)) == 0


def test_pending_rows_are_read_per_node_and_task():
    v = alt._View(0.0)
    alt._read_pending(
        v,
        [
            ["12_[1]", "p", "q", "8", "cpu=8,mem=16G,node=2,billing=8", "Priority"],
            ["13_[2]", "p", "q", "6", "mem=600M", "Resources"],
            ["14", "p", "q", "2", "cpu=2", "Priority"],
            ["15_[1-3]", "p", "q", "1", "cpu=1", "Dependency"],
        ],
    )
    rows = {(p.key, p.cpus, p.mem) for p in v.pend}
    # (per node; NumCPUs without a cpu count; a job without array is a row without key; an
    # array a row per task)
    assert rows == {
        ("12_1", 4.0, 8192.0),
        ("13_2", 6.0, 600.0),
        ("", 2.0, 0.0),
        ("15_1", 1.0, 0.0),
        ("15_2", 1.0, 0.0),
        ("15_3", 1.0, 0.0),
    }
    assert v.pend_index["15_3"] == 5


@POSIX
def test_pending_jobs_without_array_start_first_too(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 16), gpu_node("g1", 0, 0)])
    fake_slurm.set_job_states("2001", ["PD:Priority"])
    fake_slurm.set_job_attributes("2001", "cpu_queue", 8, 8 * 21418)
    with contextlib.redirect_stdout(io.StringIO()):
        assert choose(abc_first(0)) == 1


@POSIX
def test_nodes_without_memory_give_no_memory_per_cpu(fake_slurm):
    # (a place's nodes that report no memory: the main place's memory per CPU)
    zero = FakeSlurm.node_line("g1", "gpu8", 128, 0, 0, 0, 8, 0, "MIXED", "gpu8")
    fake_slurm.set_nodes([abc_node("n1", 0), zero])
    places = S().places(main_settings(), {"SlurmParam": "-p gpu8 --qos qos_normal -n1"})
    assert places[1].mem_per_cpu == 20.9
    assert alt.expand_hostlist("n[8-10]") == ["n8", "n9", "n10"]


@POSIX
def test_an_unknown_view_asks_no_place_without_partition(fake_slurm):
    fake_slurm.set_nodes([abc_node("n1", 24), gpu_node("g1", 0, 0)])
    entry = {"SlurmParam": "--qos qos_normal -n1 --mem-per-cpu=8000M", "cpuJobs": True}
    places = S().places(main_settings(), entry)
    with contextlib.redirect_stdout(io.StringIO()):
        assert choose(places) == 1
    n = fake_slurm.call_count("sbatch")
    fake_slurm.fail("squeue")
    with contextlib.redirect_stderr(io.StringIO()):
        assert choose(places) == 0
    assert fake_slurm.call_count("sbatch") == n


def test_nan_is_matlabs():
    from petakit5d.core.dtypes import NAN

    def bits(x: float) -> bytes:
        import struct

        return struct.pack("<d", x)

    assert bits(alt._NAN) == bits(NAN)
    for x in (
        alt._str2double("x"),
        alt.memory_mb("1.2.3G"),
        alt._tres_count("cpu=1", "mem"),
        alt._mdiv(0.0, 0.0),
        alt.parse_options("-n x").ntasks,
    ):
        assert bits(x) == bits(NAN)
    assert alt.expand_hostlist("n[]") == ["nNaN"]
    assert alt._pending_keys("12_[3-x]") == (["12_NaN"], 1)


# --- one place behaves exactly as without alternatives -------------------------------------------

# (tasks, run_tasks parameters, job states with the tasks whose outputs a job makes when it
# first runs)
ONE_PLACE_SCENARIOS = {
    "jobs_pend_run_and_finish": (
        3,
        {"query_interval": 0},
        {"1000_1": (["PD", "R"], [1]), "1001_2": (["R"], [2]), "1002_3": (["PD", "PD", "R"], [3])},
    ),
    "retries_get_more_cpus": (
        1,
        {"cpus_per_task": 2, "job_time_limit": 1.5, "max_trial_num": 3},
        {"1000_1": (["R", "done"], []), "1001_1": (["R", "done"], []), "1002_1": (["R"], [1])},
    ),
    "master_takes_the_last_batch": (
        4,
        {"master_compute": True, "query_interval": 0, "max_trial_num": 3},
        {
            "1000_1": (["PD", "PD", "R"], [1]),
            "1001_2": (["PD", "R"], [2]),
            "1002_3": (["PD"], []),
        },
    ),
    "gpu_jobs": (
        2,
        {"gpu_job": True, "cpus_per_task": 4, "max_trial_num": 2},
        {
            "1000_1": (["R", "done"], [1]),
            "1001_2": (["R", "done"], []),
            "1002_2": (["PD", "R"], [2]),
        },
    ),
    "in_job_slots": (
        4,
        {"gnu_parallel": True, "para_job_num": 2, "max_trial_num": 2},
        {"1000_1": (["R", "done"], []), "1001_2": (["R"], [3, 4]), "1002_1": (["R"], [1, 2])},
    ),
    "nice_and_memory": (
        2,
        {"nice": 3, "mem_allocate": 70},
        {"1000_1": (["PD", "R"], [1]), "1001_2": (["R"], [2])},
    ),
    "whole_node_jobs": (
        4,
        {"whole_node_job": True, "cpus_per_task": 6},
        {"1000_1": (["R"], [1, 2, 3, 4])},
    ),
}


def one_place_run(fake: FakeSlurm, folder: Path, spec: tuple, conf: dict, name: str, **extra):
    """a framework call of a scenario with a config; its done flags and all its commands (the
    call folder's name as <CALL>)"""
    n, params, states = spec
    fake.clear_calls()
    (folder / "out").mkdir(parents=True, exist_ok=True)
    for f in (folder / "out").glob("out_*"):
        f.unlink()
    outs = [str(folder / "out" / f"out_{k}.txt") for k in range(1, n + 1)]
    for key, (st, files) in states.items():
        fake.set_job_states(key, st, [outs[k - 1] for k in files] or None)
    ins, outs, _ = make_tasks(folder, n)
    tasks = [Task(TOUCH, [o]) for o in outs]
    settings = {"uuid": "u1", "query_interval": 3, **params, **extra}
    flags = framework(ins, outs, tasks, folder, write_config(folder, conf, name), **settings)
    return flags, [_CALL_FOLDER.sub("<CALL>", x) for x in fake.history()]


@POSIX
@pytest.mark.parametrize("scenario", sorted(ONE_PLACE_SCENARIOS))
def test_one_place_gives_the_commands_of_no_alternatives(
    scenario, fake_slurm, tmp_path, virtual_clock
):
    """a config with one place (no key, [], {}, null, and run_tasks' parameter) gives the same
    Slurm commands in the same order, and the done flags; and a whole-node call ignores places"""
    spec = ONE_PLACE_SCENARIOS[scenario]
    fake_slurm.set_nodes([abc_node("n1", 0), gpu_node("g1", 0, 0)])
    base = {"SlurmParam": ABC, "memPerCPU": 20.9, "maxCPUNum": 24}
    variants = [dict(base, SlurmAlternatives=v) for v in ([], {}, None)]
    if spec[1].get("whole_node_job"):
        variants.append(abc_first_config(0))
    want = one_place_run(fake_slurm, tmp_path, spec, base, "base.json")
    assert any(x.startswith("sbatch") for x in want[1])
    assert not any(
        x.startswith(("scontrol", "squeue --me")) or "--test-only" in x for x in want[1]
    ), want[1]
    for i, conf in enumerate(variants):
        assert one_place_run(fake_slurm, tmp_path, spec, conf, f"v{i}.json") == want, conf
    for i, value in enumerate(([], {}, None)):
        got = one_place_run(
            fake_slurm, tmp_path, spec, base, f"p{i}.json", slurm_alternatives=value
        )
        assert got == want, value


# --- oracle: MATLAB's local functions on grids, and scenarios on a fake cluster ------------------

# MATLAB: the local functions of slurm_alternatives.m of the oracle tree as function files in a
# temporary folder on the path, around BODY
_LOCAL = """
src__ = fileread(which('slurm_alternatives'));
lines__ = regexp(src__, '\\n', 'split');
starts__ = find(startsWith(lines__, 'function '));
dir__ = tempname;
mkdir(dir__);
names__ = {};
for i__ = 2 : numel(starts__)
    last__ = numel(lines__);
    if i__ < numel(starts__)
        last__ = starts__(i__ + 1) - 1;
    end
    name__ = regexp(lines__{starts__(i__)}, '=\\s*(\\w+)\\s*\\(', 'tokens', 'once');
    names__{end + 1} = name__{1};
    fid__ = fopen(fullfile(dir__, [name__{1}, '.m']), 'w');
    fprintf(fid__, '%s\\n', lines__{starts__(i__) : last__});
    fclose(fid__);
end
addpath(dir__);
try
BODY
catch err__
    rmpath(dir__);
    clear(names__{:});
    rmdir(dir__, 's');
    rethrow(err__);
end
rmpath(dir__);
clear(names__{:});
rmdir(dir__, 's');
"""


def _local_code(body: str) -> str:
    return _LOCAL.replace("BODY", body)


OPTION_TEXTS = [
    "",
    ABC,
    "-A lab --partition=gpu4 -q high -n1 --mem=78G --gres=gpu:a100:2 -C c24 -x g3.cluster "
    "-J pk -t 1:00:00 --exclusive",
    "--partition gpu8 --qos=q --account=a --mem-per-cpu 7800 --gpus=a100:3",
    "-pgpu8 -qhigh -Aacct -Cc24 -xn1 -wn[1-3] -Jname -n2 -t2:00:00 -G2",
    "--gres=gpu",
    "--gres=gpu:a100",
    "--gres=gpu:2,shard:1",
    "--gres=shard:1,gpu:0",
    "--gres=gpu:a100:x",
    "--gpus-per-node=4",
    "--gpus=2 --gres=gpu:1",
    "--mem=1T --mem-per-cpu=512K",
    "--mem=10GB",
    "--mem=10gb --mem-per-cpu=.5G",
    "--mem=abc",
    "--mem-per-cpu=1.5G",
    "--mem-per-cpu=1.2.3G",
    "--ntasks=4 -n 2",
    "-n x",
    "--nice=5 -c 4 --exclusive --comment=a=b",
    "  -p   a    -q  b  ",
    "-p\tgpu8\t--qos  q",
    "-p",
    "--qos",
    "-",
    "-- -p x",
    "--constraint=c24&tv --exclude=g[1-2] --nodelist=n1,n2 --job-name=x --time=1-00:00:00",
    "-C 'a|b'",
    "--mem-per-cpu=7800M -x n1 -x n2",
    "-t 5 --time=10",
    "x --mem-per-cpu=21418M",
    "--partition=",
    '-J "my job" --account',
    "--mem 5G --mem-per-cpu",
]
# slurm_alternatives('with_nice', SlurmParam, nice): SlurmParams with and without a --nice of
# their own (MATLAB's regexprep: its whitespace, the ends of the text, tokens that only contain
# --nice), nices rounded half away from zero
WITH_NICE = [
    (ABC, 5),
    (f"{ABC} --nice=2", 5),
    (f"--nice=2 {ABC}", 5),
    ("-p cpu_queue --nice", 2.5),
    ("-p cpu_queue --nicex -n1", 1),
    ("-p cpu_queue\t--nice=3\t--nice=4 -n1", 7),
    ("", 3),
    ("-p cpu_queue  ", 1),
    ("--nice", 1),
    ("-p cpu_queue --nice=1 -t 1:00:00 --nice=2", 1000),
    ("-p cpu_queue --comment=--nice=3", 4),
    ("-p cpu_queue --nice=", 0.5),
    ("-p cpu_queue --nice=7x -n1", 2.4999),
    ("-p cpu_queue -n1 --nice=1\n", 3),
]
MEMORY_PAIRS = [
    (ABC, "--mem=1G"),
    ("-p a --mem 5G -n1", "--mem-per-cpu=1M"),
    ("-p a", "--mem=1G"),
    ("--mem=1 -p a --mem-per-cpu=2", "--mem=3"),
    ("", "--mem=3"),
    ("  -p a  --mem-per-cpu 7800M ", "--mem-per-cpu=21418M"),
    ("-p a --mem", "--mem=2G"),
    ("--mem-per-cpu=1 --mem-per-cpu=2 --mem 3 x", "--mem=4"),
    ('-J "x --mem 5G y"', "--mem=1"),
    ("-p a --mem-per-gpu=5G", "--mem=1"),
]
HOST_LISTS = [
    "n1",
    "a,b",
    "nd0[01-03,07].x",
    "g[1-2]",
    "g1,g2.cluster",
    "n[1-2],m[3-4]",
    "n[08-10]",
    "a[1]b",
    "[1-2]",
    "",
    "n[1-2]x[3-4]",
    "n[]",
    "n[3-1]",
    "n[1-3,,5]",
    "n[1-x]",
    "n[a-c]",
    "nd[0001-0003]x",
    "n[8-10]",
    "n[9-10,011-012]",
]
PENDING_IDS = [
    "12_[3]",
    "12_3",
    "12",
    "12_[1-4,7%2]",
    "12_[1-3]%4",
    "12_[]",
    "12_[5-3]",
    "12_[1,1]",
    "12_[0-1500]",
    "x_1",
    "12_[1--3]",
    "12_[,2]",
    "12_[3-x]",
    "12_[x-3]",
    "12_[2.5]",
]
MEMORY_TEXTS = [
    "21418M",
    "7800",
    "32G",
    "1T",
    "512K",
    "10GB",
    "10gb",
    "1.5G",
    "1.2.3G",
    "",
    "abc",
    "5MB",
    "5B",
    "G",
    ".5G",
    "0",
    "10 G",
    "1e3M",
]
GPU_COUNTS = ["1", "a100:3", "a100", "", "x2", "3:x", "007"]
TRES = [
    ("cpu=4,mem=10G,gres/gpu=1", "cpu"),
    ("cpu=4,mem=10G,gres/gpu=1", "gres/gpu"),
    ("cpu=4,mem=10G,gres/gpu=1", "node"),
    ("cpu=4", "gres/gpu"),
    ("cpu=4", "mem"),
    ("node=2,cpu=8", "node"),
    ("gres/gpu:a100=2,gres/gpu=2", "gres/gpu"),
    ("billing=4,cpu=x", "cpu"),
    ("", "cpu"),
    ("", "node"),
    ("xcpu=4,cpu=2", "cpu"),
]


def _opts_struct(o: alt.Options) -> dict:
    return {
        "partition": o.partition,
        "qos": o.qos,
        "account": o.account,
        "memPerCPU": o.mem_per_cpu,
        "mem": o.mem,
        "memToken": o.mem_token,
        "gres": o.gres,
        "gpus": float(o.gpus),
        "constraint": o.constraint,
        "exclude": o.exclude,
        "nodelist": o.nodelist,
        "jobName": o.job_name,
        "ntasks": float(o.ntasks),
        "time": o.time,
        "nice": o.nice,
        "other": list(o.other),
    }


def py_parse_options(texts):
    return [_opts_struct(alt.parse_options(t)) for t in texts]


def py_with_memory(params, tokens):
    return [alt.with_memory(p, t) for p, t in zip(params, tokens, strict=True)]


def py_with_nice(params, nices):
    return [alt.with_nice(p, n) for p, n in zip(params, nices, strict=True)]


def py_expand_hostlist(texts):
    return [alt.expand_hostlist(t) for t in texts]


def py_pending_keys(ids):
    out = [alt._pending_keys(x) for x in ids]
    return [",".join(k) for k, _ in out], np.array([float(n) for _, n in out]).reshape(1, -1)


def py_numbers(memory, gpus, tres, names):
    return (
        [alt.memory_mb(x) for x in memory],
        np.array([alt._gpu_count(x) for x in gpus]).reshape(1, -1),
        np.array([alt._tres_count(t, n) for t, n in zip(tres, names, strict=True)]).reshape(1, -1),
    )


GRID_MAIN = {
    "SlurmParam": ABC,
    "SlurmConstraint": "",
    "memPerCPU": 20.9,
    "maxCPUNum": 24,
    "minCPUNum": 1,
    "gpuHeadroom": None,
}
GRID_ALTS = [
    {"SlurmParam": GPU, "memPerCPU": 7.6, "maxCPUNum": 32},
    {"SlurmParam": "-p gpu8 --qos qos_high -n1 --mem-per-cpu=7800M", "memPerCPU": 7.6},
    {"SlurmParam": "-p gpu8 -n1 --mem-per-cpu=7800M", "memPerCPU": 7.6},
    {"SlurmParam": f"-A acct2 {ABC}", "memPerCPU": 20.9},
    {"SlurmParam": "-p cpu_queue --qos qos_normal -n1 --mem=500G", "memPerCPU": 20.9},
    {"SlurmParam": "-p cpu_queue --qos qos_normal -n1", "memPerCPU": 20.9},
    {"SlurmParam": f"{GPU} --gres=gpu:1", "memPerCPU": 7.6},
    {
        "SlurmParam": "-p gpu4 --qos qos_normal -n1 --mem-per-cpu=31000M --gres=gpu:a100:1 -x a2",
        "memPerCPU": 30,
    },
    {"SlurmParam": f"{ABC} -C c48 -w n[3-4] -J pk", "memPerCPU": 20.9},
    {"SlurmParam": f"{ABC} --exclusive", "memPerCPU": 20.9},
    {"SlurmParam": "-p cpu_queue --qos qos_normal -n2 --mem-per-cpu=21418M", "memPerCPU": 20.9},
    {"SlurmParam": f"{ABC} -t 2:00:00", "memPerCPU": 20.9},
    {"SlurmParam": "--qos qos_normal -n1 --mem-per-cpu=8000M", "memPerCPU": 7.8},
    {"SlurmParam": "-p big --qos qos_normal -n1 --mem-per-cpu=42000M -A acct2", "memPerCPU": 41},
    {"SlurmParam": f"{GPU} --gres=gpu:2", "memPerCPU": 7.6},
    {"SlurmParam": f"{GPU} --gpus=1", "memPerCPU": 7.6},
    {"SlurmParam": "-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M", "memPerCPU": 20.7},
    {"SlurmParam": "-p gpu8 --qos qos_normal -n2 --mem-per-cpu=7800M", "memPerCPU": 7.6},
    {"SlurmParam": "-p big --qos qos_normal -n1 --mem-per-cpu=21418.5M", "memPerCPU": 20.9},
    {"SlurmParam": "-p big --qos qos_normal -n1 --mem=1000.5", "memPerCPU": 30},
]
GRID_CPUS = [4, 10, 10, 10, 4, 4, 4, 10, 2, 4, 4, 4, 4, 10, 2, 10, 10, 5, 10, 4, 3]

_MOVE_GRID = """
clear slurm_alternatives
places = slurm_alternatives('places', in{1}, in{2});
cpus = in{3};
F = {};
S = {};
for a = 1 : numel(places)
    for b = 1 : numel(places)
        if a == b
            continue;
        end
        cur = place_shape(places, cpus, a, main_shape(places, cpus));
        shp = place_shape(places, cpus, b, cur);
        F{end + 1} = move_fields(places(a), places(b), cur, shp);
        S{end + 1} = shp;
    end
end
out = {F, [S{:}]};
"""


def py_move_grid(main, alts, cpus):
    places = Session().places(main, alts)
    fields, shapes = [], []
    for a in range(len(places)):
        for b in range(len(places)):
            if a == b:
                continue
            cur = alt._place_shape(places, cpus, a, alt._main_shape(places, cpus))
            shape = alt._place_shape(places, cpus, b, cur)
            f = alt.move_fields(places[a], places[b], cur, shape)
            fields.append(" ".join(f"{k}={v}" for k, v in f.items()))
            shapes.append(
                {
                    "cpus": float(shape.cpus),
                    "memToken": shape.mem_token,
                    "ratio": shape.ratio,
                    "param": shape.param,
                }
            )
    return fields, shapes


# MATLAB: the steps of a scenario (in{2}, JSON) against the fake Slurm in in{1}
_SCENARIO = """
fake__ = in{1};
if fake__(end) == '/'
    fake__ = fake__(1 : end - 1);
end
spec__ = jsondecode(in{2});
cmds__ = {'sinfo', 'sbatch', 'squeue', 'scancel', 'scontrol'};
for i__ = 1 : numel(cmds__)
    fileattrib(fullfile(fake__, cmds__{i__}), '+x');
end
path__ = getenv('PATH');
setenv('PATH', [fake__, ':', path__]);
clear slurm_alternatives
K = [];
C = [];
P = {};
W = [];
try
    [~, w__] = system('command -v sbatch');
    if ~strcmp(strtrim(w__), [fake__, '/sbatch'])
        error('j09:fake', 'the fake Slurm is not first on the PATH: %s', strtrim(w__));
    end
    steps__ = spec__.steps;
    if isstruct(steps__)
        steps__ = num2cell(steps__);
    end
    places = [];
    last__ = [];
    for s__ = 1 : numel(steps__)
        st__ = steps__{s__};
        switch st__.op
            case 'set'
                slurm_alternatives('set', 'viewAge', st__.viewAge, 'minWait', st__.minWait);
                slurm_alternatives('set', 'lookInterval', st__.lookInterval);
            case 'places'
                places = slurm_alternatives('places', st__.main, st__.alts);
                w__ = [[places.memPerCPU]; [places.maxCPUNum]; [places.minCPUNum]];
                W = [W, [w__; double([places.cpuJobs])]];
            case 'choose'
                times__ = cellstr(st__.times);
                [k__, c__, p__] = slurm_alternatives('choose', places, st__.cpus(:)', times__);
                K(end + 1) = k__;
                C(end + 1) = c__;
                P{end + 1} = p__;
                last__ = struct('k', k__, 'cpus', st__.cpus(:)', 'times', {times__});
            case 'submit'
                a__ = {st__.job, st__.array, places, last__.cpus, last__.k, last__.times};
                slurm_alternatives('submitted', a__{:});
                id__ = sprintf('%d_%d', st__.job, st__.array);
                writelines(cellstr(st__.states), fullfile(fake__, 'state', id__));
                t__ = regexp(P{end}, '--mem-per-cpu=(\\d+)M', 'tokens', 'once');
                mem__ = 1000 * C(end);
                if ~isempty(t__)
                    mem__ = str2double(t__{1}) * C(end);
                end
                t__ = regexp(P{end}, '--gres=gpu(?::\\w+)?:(\\d+)', 'tokens', 'once');
                gpus__ = 0;
                if ~isempty(t__)
                    gpus__ = str2double(t__{1});
                end
                o__ = places(last__.k).opts;
                line__ = sprintf('partition=%s qos=%s', o__.partition, o__.qos);
                line__ = sprintf('%s cpus=%d mem=%d gpus=%d', line__, C(end), mem__, gpus__);
                writelines(line__, fullfile(fake__, 'jobs', id__));
            case 'move'
                slurm_alternatives('move', st__.jobs(:), st__.arrays(:));
            case 'forget'
                slurm_alternatives('forget', st__.jobs(:), st__.arrays(:));
            case 'nodes'
                d__ = dir(fullfile(fake__, 'nodes_*'));
                for i__ = 1 : numel(d__)
                    delete(fullfile(fake__, d__(i__).name));
                end
                writelines(cellstr(st__.lines), fullfile(fake__, 'nodes_1'));
            case 'state'
                writelines(cellstr(st__.states), fullfile(fake__, 'state', st__.id));
            case 'attrs'
                line__ = sprintf('partition=%s qos=%s', st__.partition, st__.qos);
                line__ = sprintf('%s cpus=%d mem=%d', line__, st__.cpus, st__.mem);
                line__ = sprintf('%s gpus=%d', line__, st__.gpus);
                writelines(line__, fullfile(fake__, 'jobs', st__.id));
            case 'testonly'
                writelines(st__.answer, fullfile(fake__, 'testonly', st__.partition));
            case 'failupdates'
                fclose(fopen(fullfile(fake__, 'update_fail'), 'w'));
        end
    end
catch err__
    setenv('PATH', path__);
    rethrow(err__);
end
setenv('PATH', path__);
H = {};
fn__ = fullfile(fake__, 'history.log');
if isfile(fn__)
    H = regexp(fileread(fn__), '\\n', 'split');
    H = H(1 : end - 1);
end
d__ = dir(fullfile(fake__, 'jobs'));
names__ = sort({d__(~[d__.isdir]).name});
J = cell(1, numel(names__));
for i__ = 1 : numel(names__)
    J{i__} = [names__{i__}, ' ', strtrim(fileread(fullfile(fake__, 'jobs', names__{i__})))];
end
out = {K, C, strjoin(P, newline), W, strjoin(H, newline), strjoin(J, newline)};
"""


def _demand(param: str, cpus: int) -> tuple[int, int]:
    """a submitted job's memory (MB) and GPUs in a scenario (the fake's record)"""
    m = re.search(r"--mem-per-cpu=(\d+)M", param)
    mem = int(m.group(1)) * cpus if m else 1000 * cpus
    g = re.search(r"--gres=gpu(?::\w+)?:(\d+)", param)
    return mem, int(g.group(1)) if g else 0


def py_scenario(fake_dir: str, spec_text: str):
    """the steps of a scenario with a new Session against the fake in fake_dir: the choices
    (1-based places, CPUs, SlurmParams), the places (memPerCPU, maxCPUNum, minCPUNum, cpuJobs),
    the fake's commands and its jobs at the end (texts: lines joined by newlines)"""
    if os.name == "nt":
        pytest.skip("The fake Slurm needs a POSIX shell.")
    fake = FakeSlurm(fake_dir.rstrip("/"))
    spec = json.loads(spec_text)
    old = os.environ.get("PATH", "")
    os.environ["PATH"] = fake.path_with(old)
    ks, cs, ps, ws = [], [], [], []
    try:
        session = Session()
        places, last = None, None
        for st in spec["steps"]:
            op = st["op"]
            if op == "set":
                session.set_timing(
                    view_age=st["viewAge"], min_wait=st["minWait"], look_interval=st["lookInterval"]
                )
            elif op == "places":
                places = session.places(st["main"], st["alts"])
                ws += [
                    [p.mem_per_cpu, p.max_cpu_num, p.min_cpu_num, float(p.cpu_jobs)] for p in places
                ]
            elif op == "choose":
                k, c, p = session.choose(places, st["cpus"], st["times"])
                ks.append(k + 1)
                cs.append(c)
                ps.append(p)
                last = (k, list(st["cpus"]), list(st["times"]))
            elif op == "submit":
                k, cpus, times = last
                session.submitted(st["job"], st["array"], places, cpus, k, times)
                key = f"{st['job']}_{st['array']}"
                fake.set_job_states(key, st["states"])
                mem, gpus = _demand(ps[-1], cs[-1])
                o = places[k].opts
                fake.set_job_attributes(key, o.partition, cs[-1], mem, gpus, o.qos)
            elif op == "move":
                session.move(st["jobs"], st["arrays"])
            elif op == "forget":
                session.forget(st["jobs"], st["arrays"])
            elif op == "nodes":
                fake.set_nodes(st["lines"])
            elif op == "state":
                fake.set_job_states(st["id"], st["states"])
            elif op == "attrs":
                fake.set_job_attributes(
                    st["id"], st["partition"], st["cpus"], st["mem"], st["gpus"], st["qos"]
                )
            elif op == "testonly":
                fake.set_test_only(st["partition"], st["answer"])
            elif op == "failupdates":
                fake.fail_updates()
    finally:
        os.environ["PATH"] = old
    jobs = sorted(os.listdir(fake.folder / "jobs"))
    js = [f"{n} {(fake.folder / 'jobs' / n).read_text().strip()}" for n in jobs]
    w = np.array(ws, dtype=float).T if ws else np.zeros((0, 0))
    return ks, cs, "\n".join(ps), w, "\n".join(fake.history()), "\n".join(js)


def make_bare_fake(folder: Path) -> None:
    FakeSlurm(folder).install()


# the cluster of the scenarios: (name, partitions, CPUs, MB, GPUs, features)
SCENARIO_NODES = [
    ("n1", "cpu_queue", 24, 514048, 0, "c24"),
    ("n2", "cpu_queue", 24, 514048, 0, "c24"),
    ("n3", "cpu_queue,big", 48, 1031745, 0, "c48"),
    ("t1", "cpu_queue", 20, 515638, 4, "tv"),
    ("g1", "gpu8", 128, 1031600, 8, "gpu8"),
    ("g2", "gpu8", 128, 1031600, 8, "gpu8"),
    ("a1", "gpu4", 16, 514048, 4, "a100"),
    ("a2", "gpu4", 16, 514048, 4, "a100"),
    ("m1", "big", 48, 1031745, 0, "c48"),
]
STATES = ["MIXED"] * 9 + ["IDLE", "ALLOCATED", "MIXED+PLANNED", "IDLE+DRAIN", "DOWN+NOT_RESPONDING"]
# (SlurmParam, memPerCPU or None, maxCPUNum or None) of the places of CPU jobs and of GPU jobs
CPU_PLACES = [
    (ABC, 20.9, 24),
    ("-p cpu_queue --qos qos_high -n1 --mem-per-cpu=21418M", 20.9, 24),
    ("-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M -C c48", None, 48),
    # (quoted: MATLAB's shell made a bare | a pipe, see test_options_with_shell_characters_...)
    ("-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M -C 'c24|c48'", None, None),
    ("-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M -x n[1-2]", None, None),
    ("-p cpu_queue --qos qos_normal -n1 --mem-per-cpu=21418M -w n1,t1", 20.9, None),
    ("-p big --qos qos_normal -n1 --mem=200G", None, 48),
    ("-p big --qos qos_normal -n1", None, 48),
    (GPU, 7.6, 32),
    ("-p gpu8 -n1 --mem-per-cpu=7800M", 7.6, 32),
    ("-p gpu4 --qos qos_normal -n1 --mem-per-cpu=31000M", 30.0, 16),
    ("--qos qos_normal -n1 --mem-per-cpu=8000M", 7.8, None),
    (f"{ABC} -t 2:00:00", 20.9, 24),
    (f"-A acct2 {ABC}", 20.9, 24),
]
GPU_PLACES = [
    (f"{GPU} --gres=gpu:1", 7.6, 32),
    ("-p gpu4 --qos qos_normal -n1 --mem-per-cpu=31000M --gres=gpu:a100:1 -x a2", 30.0, 16),
    ("-p gpu4 --qos qos_normal -n1 --mem-per-cpu=31000M --gres=gpu:1", None, 16),
    ("-p gpu8 --qos qos_high -n1 --mem-per-cpu=7800M --gres=gpu:1", 7.6, 32),
    ("-p gpu8 --qos qos_normal -n1 --mem-per-cpu=7800M --gpus=1", 7.6, 32),
]
REASONS = ["Priority"] * 4 + ["Resources", "None", "JobHeldUser", "QOSMaxCpuPerUserLimit"]
FOREIGN_REASONS = REASONS + ["AssocGrpCpuLimit", "BeginTime", "Dependency", "Resources"]


def _nodes(rng: random.Random) -> list[str]:
    lines = []
    for name, parts, cpus, mem, gpus, feats in SCENARIO_NODES:
        cpu_alloc = rng.choice([0, cpus // 4, cpus // 2, cpus - 4, cpus, rng.randrange(cpus + 1)])
        cpu_alloc = max(0, min(cpus, cpu_alloc))
        mem_alloc = min(mem, cpu_alloc * rng.choice([4000, 8000, 21418, 30000]))
        gpu_alloc = rng.randrange(gpus + 1) if gpus else 0
        lines.append(
            FakeSlurm.node_line(
                name,
                parts,
                cpus,
                cpu_alloc,
                mem,
                mem_alloc,
                gpus,
                gpu_alloc,
                rng.choice(STATES),
                feats,
            )
        )
    return lines


def _entry(rng: random.Random, place: tuple) -> dict:
    param, mpc, maxc = place
    entry: dict = {"SlurmParam": param}
    if mpc is not None and rng.random() < 0.8:
        entry["memPerCPU"] = mpc
    if maxc is not None and rng.random() < 0.7:
        entry["maxCPUNum"] = maxc
    headroom = rng.choice(["-", "-", None, 0, [4, 32], 2])
    if headroom != "-":
        entry["gpuHeadroom"] = headroom
    jobs = rng.choice(["-", True, True, False, 1])
    if jobs != "-":
        entry["cpuJobs"] = jobs
    if rng.random() < 0.08:
        entry["maxCPU"] = 3
    if rng.random() < 0.1:
        entry = {k.lower() if k == "SlurmParam" else k: v for k, v in entry.items()}
    return entry


def scenario(seed: int) -> dict:
    """a scenario on the fake cluster (deterministic for a seed)"""
    rng = random.Random(seed)
    gpu_jobs = rng.random() < 0.25
    pool = GPU_PLACES if gpu_jobs else CPU_PLACES
    picks = rng.sample(range(len(pool)), rng.randint(2, min(4, len(pool))))
    main_param, main_mpc, main_max = pool[picks[0]]
    main = {
        "SlurmParam": main_param,
        "SlurmConstraint": rng.choice([""] * 8 + ["-C c24", "--constraint=c48"]),
        "memPerCPU": main_mpc or 20.9,
        "maxCPUNum": main_max or 24,
        "minCPUNum": rng.choice([1, 1, 2]),
        "gpuHeadroom": rng.choice([None, None, 0, [4, 32]]),
    }
    alts = [_entry(rng, pool[i]) for i in picks[1:]]
    steps: list[dict] = [
        {"op": "set", "viewAge": rng.choice([0, 1e6]), "minWait": 0, "lookInterval": 0}
    ]
    for part in ("cpu_queue", "gpu8", "gpu4", "big"):
        if rng.random() < 0.25:
            steps.append(
                {"op": "testonly", "partition": part, "answer": rng.choice(["later", "error"])}
            )
    steps.append({"op": "nodes", "lines": _nodes(rng)})
    for f in range(rng.randint(0, 5)):
        part = rng.choice(["cpu_queue", "gpu8", "gpu4", "big"])
        gpus = rng.choice([0, 1]) if part in ("gpu8", "gpu4") else 0
        cpus = rng.choice([1, 2, 4, 8, 16, 32])
        key = f"{2001 + f}" if rng.random() < 0.2 else f"{2001 + f}_{rng.randint(1, 9)}"
        steps.append({"op": "state", "id": key, "states": [f"PD:{rng.choice(FOREIGN_REASONS)}"]})
        steps.append(
            {
                "op": "attrs",
                "id": key,
                "partition": part,
                "cpus": cpus,
                "mem": cpus * rng.choice([4000, 8000, 21418]),
                "gpus": gpus,
                "qos": rng.choice(["qos_normal", "qos_high"]),
            }
        )
    steps.append({"op": "places", "main": main, "alts": alts})
    n_places = 1 + len(alts)
    mem = rng.choice([8, 20, 70, 100, 300])
    ours: list[tuple[int, int]] = []

    def choose_step() -> None:
        cpus = [
            max(1, min(32, math.ceil(mem / rng.choice([7.6, 20.9, 30])))) for _ in range(n_places)
        ]
        times = [rng.choice(["", "", " -t 24:00:00 "]) for _ in range(n_places)]
        steps.append({"op": "choose", "cpus": cpus, "times": times})

    for i in range(rng.randint(3, 6)):
        choose_step()
        if rng.random() < 0.7:
            job = (1001 + i, 1 + i % 9)
            ours.append(job)
            steps.append(
                {
                    "op": "submit",
                    "job": job[0],
                    "array": job[1],
                    "states": [f"PD:{rng.choice(REASONS)}"],
                }
            )
    if rng.random() < 0.6:
        steps.append({"op": "nodes", "lines": _nodes(rng)})
    if ours and rng.random() < 0.2:
        job = rng.choice(ours)
        steps.append(
            {
                "op": "attrs",
                "id": f"{job[0]}_{job[1]}",
                "partition": rng.choice(["gpu8", "big"]),
                "cpus": 4,
                "mem": 4000,
                "gpus": 0,
                "qos": "qos_normal",
            }
        )
    if rng.random() < 0.15:
        steps.append({"op": "failupdates"})
    moves = [j[0] for j in ours] + ([2001] if rng.random() < 0.3 else [])
    arrays = [j[1] for j in ours] + ([1] if len(moves) > len(ours) else [])
    if moves:
        steps.append({"op": "move", "jobs": moves, "arrays": arrays})
    for _ in range(2):
        choose_step()
    if moves:
        if rng.random() < 0.3:
            steps.append({"op": "forget", "jobs": moves[:1], "arrays": arrays[:1]})
        steps.append({"op": "move", "jobs": moves, "arrays": arrays})
    return {"steps": steps}


def _scenario_case(name: str, spec: dict, marks: tuple[str, ...] = (), note: str = "") -> Case:
    return Case(
        name,
        code=_SCENARIO,
        args=[InPath("fake"), json.dumps(spec, sort_keys=True)],
        nargout=6,
        python=py_scenario,
        inputs={"fake": Tree(make_bare_fake, "fake")},
        outputs=["K", "C", "P", "W", "H", "J"],
        files=False,
        warnings="ids",
        marks=marks,
        note=note,
    )


# a job of a call with a nice (in the main place's SlurmParam, as the framework passes it) waits
# at the main place, then gpu8 has room: the job moves there (MATLAB before its fix of the nice,
# oracle PetaKit5D_1697257 and older, kept it, as the main place's options held the nice and
# gpu8's didn't)
NICE_SCENARIO = {
    "steps": [
        {"op": "set", "viewAge": 0, "minWait": 0, "lookInterval": 0},
        {
            "op": "nodes",
            "lines": [
                FakeSlurm.node_line(
                    "n1", "cpu_queue", 24, 24, 514048, 514048, 0, 0, "MIXED", "c24"
                ),
                FakeSlurm.node_line(
                    "g1", "gpu8", 128, 128, 1031600, 1031600, 8, 0, "MIXED", "gpu8"
                ),
            ],
        },
        {
            "op": "places",
            "main": {**GRID_MAIN, "SlurmParam": f"{ABC} --nice=5"},
            "alts": [
                {
                    "SlurmParam": GPU,
                    "memPerCPU": 7.6,
                    "maxCPUNum": 32,
                    "gpuHeadroom": 0,
                    "cpuJobs": True,
                }
            ],
        },
        {"op": "choose", "cpus": [4, 10], "times": ["", ""]},
        {"op": "submit", "job": 1000, "array": 1, "states": ["PD:Resources"]},
        {
            "op": "nodes",
            "lines": [
                FakeSlurm.node_line(
                    "n1", "cpu_queue", 24, 24, 514048, 514048, 0, 0, "MIXED", "c24"
                ),
                FakeSlurm.node_line("g1", "gpu8", 128, 0, 1031600, 0, 8, 0, "MIXED", "gpu8"),
            ],
        },
        {"op": "move", "jobs": [1000], "arrays": [1]},
    ]
}


@POSIX
def test_the_nice_scenario_moves_the_job(tmp_path):
    # (the oracle case nice_job_moves_to_an_alternative, on the Python side alone)
    make_bare_fake(tmp_path / "fake")
    ks, _, ps, _, h, j = py_scenario(str(tmp_path / "fake"), json.dumps(NICE_SCENARIO))
    assert ks == [1] and ps == f"{ABC} --nice=5"
    assert h.splitlines()[-1] == (
        "scontrol update JobId=1000_1 Partition=gpu8 QOS=qos_normal NumCPUs=4 CPUsPerTask=4 "
        "MinMemoryCPU=21418"
    )
    assert j == "1000_1 partition=gpu8 qos=qos_normal cpus=4 mem=85672 gpus=0"


def _hand_scenarios() -> dict[str, dict]:
    """scenarios for paths the random ones may miss"""
    cpu_nodes = [
        FakeSlurm.node_line("n1", "cpu_queue", 24, 24, 514048, 514048, 0, 0, "MIXED", "c24"),
        FakeSlurm.node_line("g1", "gpu8", 128, 0, 1031600, 0, 8, 0, "MIXED", "gpu8"),
        FakeSlurm.node_line("m1", "big", 48, 0, 1031745, 0, 0, 0, "IDLE", "c48"),
    ]
    choose4 = {"op": "choose", "cpus": [4, 10, 10], "times": ["", "", ""]}
    main = dict(GRID_MAIN)
    out = {
        # a place without a partition: Slurm's check alone, one job per view; then refused
        "no_partition_one_job_per_view": [
            {"op": "set", "viewAge": 1e6, "minWait": 0, "lookInterval": 0},
            {"op": "nodes", "lines": cpu_nodes},
            {
                "op": "places",
                "main": main,
                "alts": [
                    {"SlurmParam": "--qos qos_normal -n1 --mem-per-cpu=8000M", "cpuJobs": True},
                    {"SlurmParam": "-p big --qos qos_normal -n1"},
                ],
            },
            choose4,
            choose4,
            {"op": "set", "viewAge": 0, "minWait": 0, "lookInterval": 0},
            {"op": "testonly", "partition": "cpu_queue", "answer": "error"},
            choose4,
        ],
        # moves: the GPU nodes' queue, a job moved by something else, a failing update
        "moves_of_every_kind": [
            {"op": "set", "viewAge": 0, "minWait": 0, "lookInterval": 0},
            {"op": "nodes", "lines": cpu_nodes},
            {"op": "state", "id": "2001_1", "states": ["PD:JobHeldUser"]},
            {
                "op": "attrs",
                "id": "2001_1",
                "partition": "gpu8",
                "cpus": 4,
                "mem": 32000,
                "gpus": 1,
                "qos": "qos_normal",
            },
            {
                "op": "places",
                "main": {**main, "SlurmParam": GPU, "memPerCPU": 7.6, "maxCPUNum": 32},
                "alts": [
                    {"SlurmParam": ABC, "memPerCPU": 20.9},
                    {"SlurmParam": "-p big --qos qos_normal -n1 --mem-per-cpu=21000M"},
                ],
            },
            {"op": "choose", "cpus": [10, 4, 4], "times": ["", "", ""]},
            {"op": "submit", "job": 1001, "array": 1, "states": ["PD:Resources"]},
            {"op": "choose", "cpus": [10, 4, 4], "times": ["", "", ""]},
            {"op": "submit", "job": 1002, "array": 2, "states": ["PD:Priority"]},
            {
                "op": "attrs",
                "id": "1002_2",
                "partition": "cpu_queue",
                "cpus": 4,
                "mem": 4000,
                "gpus": 0,
                "qos": "qos_normal",
            },
            {"op": "move", "jobs": [1001, 1002], "arrays": [1, 2]},
            {"op": "failupdates"},
            {
                "op": "nodes",
                "lines": [cpu_nodes[0], cpu_nodes[1].replace("CPUAlloc=0", "CPUAlloc=8")],
            },
            {"op": "move", "jobs": [1001], "arrays": [1]},
        ],
        # the memory per CPU of places from their nodes (mixed node types: the smallest)
        "memory_per_cpu_from_nodes": [
            {"op": "set", "viewAge": 1e6, "minWait": 0, "lookInterval": 0},
            {
                "op": "nodes",
                "lines": [
                    FakeSlurm.node_line(n, p, c, 0, m, 0, g, 0, "IDLE", f)
                    for n, p, c, m, g, f in SCENARIO_NODES
                ],
            },
            {
                "op": "places",
                "main": main,
                "alts": [
                    {"SlurmParam": "-p cpu_queue --qos qos_normal -n1 --mem=10G"},
                    {"SlurmParam": "-p cpu_queue --qos qos_normal -n1 -C c48"},
                    {"SlurmParam": "-p big --qos qos_normal -n1", "MEMPERCPU": []},
                    {"SlurmParam": "-p nowhere --qos qos_normal -n1"},
                    {"SlurmParam": "-p gpu4 --qos qos_normal -n1 -x a[1-2]"},
                ],
            },
        ],
    }
    return {k: {"steps": v} for k, v in out.items()}


ORACLE_CASES = [
    Case(
        "local_parse_options",
        code=_local_code(
            "R = cellfun(@(s) parse_options(s), in{1}, 'UniformOutput', false); out = {[R{:}]};"
        ),
        args=[OPTION_TEXTS],
        python=py_parse_options,
    ),
    Case(
        "local_with_memory",
        code=_local_code(
            "out = {cellfun(@(p, m) with_memory(p, m), in{1}, in{2}, 'UniformOutput', false)};"
        ),
        args=[[p for p, _ in MEMORY_PAIRS], [t for _, t in MEMORY_PAIRS]],
        python=py_with_memory,
    ),
    Case(
        "with_nice",
        code="out = {cellfun(@(p, n) slurm_alternatives('with_nice', p, n), in{1}(:)', "
        "num2cell(in{2}(:)'), 'UniformOutput', false)};",
        args=[[p for p, _ in WITH_NICE], [float(n) for _, n in WITH_NICE]],
        python=py_with_nice,
    ),
    Case(
        "local_expand_hostlist",
        code=_local_code(
            "out = {cellfun(@(s) expand_hostlist(s), in{1}, 'UniformOutput', false)};"
        ),
        args=[HOST_LISTS],
        python=py_expand_hostlist,
    ),
    Case(
        "local_pending_keys",
        code=_local_code(
            "[k__, n__] = cellfun(@(s) pending_keys(s), in{1}, 'UniformOutput', false); "
            "out = {cellfun(@(k) strjoin(k, ','), k__, 'UniformOutput', false), cell2mat(n__)};"
        ),
        args=[PENDING_IDS],
        nargout=2,
        python=py_pending_keys,
    ),
    Case(
        "local_memory_gpu_and_tres_counts",
        code=_local_code(
            "out = {cellfun(@(s) memory_mb(s), in{1}, 'UniformOutput', false), "
            "cellfun(@(s) gpu_count(s), in{2}), cellfun(@(t, n) tres_count(t, n), in{3}, in{4})};"
        ),
        args=[MEMORY_TEXTS, GPU_COUNTS, [t for t, _ in TRES], [n for _, n in TRES]],
        nargout=3,
        python=py_numbers,
    ),
    Case(
        "local_move_fields_and_shapes_of_place_pairs",
        code=_local_code(_MOVE_GRID),
        args=[GRID_MAIN, GRID_ALTS, GRID_CPUS],
        nargout=2,
        python=py_move_grid,
        outputs=["fields", "shapes"],
    ),
    *[_scenario_case(name, spec) for name, spec in sorted(_hand_scenarios().items())],
    *[_scenario_case(f"scenario_{seed:03d}", scenario(seed)) for seed in range(48)],
    _scenario_case("nice_job_moves_to_an_alternative", NICE_SCENARIO),
]


@oracle.parametrize(ORACLE_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_scenarios_are_deterministic():
    assert json.dumps(scenario(7)) == json.dumps(scenario(7))
    assert json.dumps(scenario(7)) != json.dumps(scenario(8))


# --- real cluster --------------------------------------------------------------------------------


class _Recorder(Slurm):
    """the real Slurm, with every command line kept"""

    def __init__(self) -> None:
        super().__init__()
        self.lines: list[list[str]] = []

    def run(self, argv, timeout=None):
        self.lines.append([str(a) for a in argv])
        return super().run(argv, timeout)


def _show(key: str) -> dict[str, str]:
    from petakit5d.jobs.slurm import parse_key_values

    r = subprocess.run(["scontrol", "show", "job", "-o", key], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return parse_key_values(r.stdout)


def _never_running(slurm: Slurm, cpus: int, param: str, constraint: str, how: str) -> str:
    """a job that never runs (how: --hold or --begin=now+1hour), its key"""
    argv = ["sbatch", how, "--array=1", "-o", "/dev/null", f"--cpus-per-task={cpus}"]
    argv += [*shlex.split(param), *shlex.split(constraint), "-t", "0:10:00", "--wrap=true"]
    r = slurm.sbatch(argv)
    assert r.submitted, r.result.output
    return f"{r.job_id}_1"


@pytest.mark.cluster
def test_real_cluster_places_view_and_moves_of_jobs_that_never_run(real_slurm):
    """On a real Slurm, with jobs that never run (held, or to begin in an hour): the view
    parses the real nodes and the user's pending jobs (ours among them, as squeue lists them),
    sbatch --test-only answers (and refuses an unknown QOS: that place is left out), a job
    submitted to the chosen place has its shape, move() leaves held and later-beginning jobs
    alone, and the fields of a move to another place (scontrol update of partition, QOS, CPUs,
    memory and node lists) are accepted and give the job that shape, and the held job keeps its
    nice (a call's, --nice=5) through every move. Every update and cancel names our jobs only;
    they are cancelled by their ids.

    The config (PETAKIT5D_TEST_SLURM_CONFIG, outside the repo) gives SlurmParam and
    SlurmAlternatives: places between which a pending job may move."""
    conf = dict(config.ClusterConfig.load(os.environ["PETAKIT5D_TEST_SLURM_CONFIG"]).entries)
    if not conf.get("SlurmAlternatives"):
        pytest.skip("the test config has no SlurmAlternatives")
    slurm = _Recorder()
    session = Session(slurm=slurm)
    main = {
        "SlurmParam": conf["SlurmParam"],
        "SlurmConstraint": conf.get("SlurmConstraint", ""),
        "memPerCPU": conf.get("memPerCPU", 20.9),
        "maxCPUNum": conf.get("maxCPUNum", 24),
        "minCPUNum": conf.get("minCPUNum", 1),
        "gpuHeadroom": conf.get("gpuHeadroom"),
    }
    places = session.places(main, conf["SlurmAlternatives"])
    t0 = time.perf_counter()
    view = session.current_view()
    t_view = time.perf_counter() - t0
    assert view.known and view.names and len(places) >= 2
    assert all(c > 0 for c in view.cpu_tot) and all(m > 0 for m in view.mem_tot)
    own = [1] * len(places)
    times = [" -t 0:10:00 "] * len(places)
    with contextlib.redirect_stdout(io.StringIO()):
        k, cpus, param = session.choose(places, own, times)
    keys: list[str] = []
    moved = 0
    try:
        nice = f"{param} --nice=5"
        held = _never_running(slurm, cpus, nice, places[k].slurm_constraint, "--hold")
        keys.append(held)
        later = _never_running(slurm, cpus, param, places[k].slurm_constraint, "--begin=now+1hour")
        keys.append(later)
        for key in keys:
            session.submitted(int(key.split("_")[0]), 1, places, own, k, times)
        job = _show(held)
        assert (job["JobState"], job["Reason"]) == ("PENDING", "JobHeldUser"), job
        assert job["NumCPUs"] == str(cpus) and job["Nice"] == "5", job
        if places[k].opts.partition:
            assert job["Partition"] == places[k].opts.partition
        # the view lists them as pending jobs of the user, with their shape and reason
        session.set_timing(view_age=0, min_wait=0, look_interval=0)
        view = session.current_view()
        rows = {key: view.pend[view.pend_index[key]] for key in keys}
        assert rows[held].reason == "JobHeldUser" and rows[held].cpus == cpus, rows
        # (BeginTime once the scheduler looked at it; None before: neither moves)
        assert rows[later].reason in ("BeginTime", "None"), rows
        assert rows[held].partition == job["Partition"]
        # held and later-beginning jobs stay where they are
        with contextlib.redirect_stdout(io.StringIO()):
            session.move([int(x.split("_")[0]) for x in keys], [1, 1])
        assert not any(x[:2] == ["scontrol", "update"] for x in slurm.lines)
        # the fields of a move to each other place, as move() applies them: Slurm takes them
        e = session.registry[held]
        for b in range(len(places)):
            if b == e.k:
                continue
            shape = alt._place_shape(places, own, b, e.shape)
            fields = alt.move_fields(places[e.k], places[b], e.shape, shape)
            if not fields:
                continue
            u = slurm.update_job(held, fields)
            assert u.ok, (fields, u.output)
            job = _show(held)
            assert (job["JobState"], job["Reason"]) == ("PENDING", "JobHeldUser"), job
            assert job["NumCPUs"] == fields["NumCPUs"] and job["Partition"] == fields["Partition"]
            if "MinMemoryCPU" in fields:
                assert alt.memory_mb(job["MinMemoryCPU"]) == float(fields["MinMemoryCPU"]), job
            if "QOS" in fields:
                assert job["QOS"] == fields["QOS"], job
            # (an update leaves the job's nice)
            assert job["Nice"] == "5", job
            e.k, e.shape = b, shape
            moved += 1
        assert moved >= 1, "no two places of the config can trade jobs"
        # Slurm refuses an unknown QOS: that place is left out (with MATLAB's warning)
        bad = dict(conf["SlurmAlternatives"][0])
        param = re.sub(r"(--qos[= ]|-q ?)\S+", "", bad["SlurmParam"])
        bad["SlurmParam"] = f"{param} --qos pk5d_no_such_qos"
        refused = session.places(main, [bad])
        with pytest.warns(AlternativesWarning) as record:
            shape = alt._place_shape(refused, [1, 1], 1, alt._main_shape(refused, [1, 1]))
            assert not session._slurm_says_now(view, refused[1], shape, "", False)
        assert identifiers(record) == ["slurm_alternatives:testOnly"]
        assert refused[1].slurm_param in session.disabled
    finally:
        for key in keys:
            subprocess.run(["scancel", key.split("_")[0]], capture_output=True)
    ids = [x.split("_")[0] for x in keys]
    deadline = time.time() + 60
    while time.time() < deadline:
        left = subprocess.run(["squeue", "-h", "-j", ",".join(ids)], capture_output=True, text=True)
        if not left.stdout.strip():
            break
        time.sleep(2)
    assert not left.stdout.strip(), left.stdout
    for x in slurm.lines:
        if x[:2] == ["scontrol", "update"]:
            assert x[2].split("=")[1].split("_")[0] in ids, x
        assert x[0] != "scancel", x
    print(
        f"real cluster: {len(view.names)} nodes, {len(view.pend)} pending rows of the user, "
        f"view read in {t_view:.3f} s; places {[p.name for p in places]}, jobs {ids}, "
        f"moves applied {moved}"
    )
