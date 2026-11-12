"""Claims of the task runners (jobs.claims): atomic symbolic links in the call folder, their
owners and whether they may still run their tasks (MATLAB's TaskRunnerTest, the claim parts)."""

import os
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from petakit5d.jobs import claims
from petakit5d.jobs.claims import Claims
from petakit5d.jobs.writers import node_name

pytestmark = pytest.mark.skipif(os.name == "nt", reason="task runners run on Linux and macOS")

SRC_ROOT = str(Path(claims.__file__).resolve().parents[2])


@pytest.fixture
def no_job(monkeypatch):
    for name in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURMD_NODENAME"):
        monkeypatch.delenv(name, raising=False)


def test_owner_names_job_node_process_and_slot(monkeypatch):
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "4438625")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "7")
    monkeypatch.setenv("SLURM_JOB_ID", "4438632")
    monkeypatch.setenv("SLURMD_NODENAME", "node12")
    assert claims.owner(3) == f"4438625_7:node12:{os.getpid()}:3"
    assert claims.parse_owner(claims.owner(3)) == (4438625, 7, "node12", os.getpid(), "3")
    assert claims.job_key_of(claims.owner(3)) == "4438625_7"


def test_owner_outside_slurm_has_job_zero(no_job):
    who = claims.owner()
    assert who.startswith("0_0:") and who.endswith(f":{os.getpid()}:1")
    assert claims.job_key_of(who) is None
    assert claims.parse_owner("not an owner") is None


def test_a_task_is_claimed_once(tmp_path, no_job):
    c = Claims(tmp_path)
    assert c.claim(0, "a") and not c.claim(0, "b")
    assert c.owner_of(0) == "a" and c.owner_of(1) is None
    assert c.names() == {"t1"}
    assert claims.claimed_set(c.names()) == {0}


def test_release_removes_claims_and_failed_marks_and_marks_resubmitted(tmp_path):
    c = Claims(tmp_path)
    for i in range(4):
        c.claim(i, "77_1:n1:5:1")
    c.mark_failed(1, "77_1:n1:5:1")
    c.release([1, 2])
    c.mark_resubmitted(3)
    assert c.names() == {"t1", "t4", "r3"}
    assert c.resubmitted(3) and not c.resubmitted(1)
    assert claims.resubmitted_set(c.names()) == {3}


def test_give_back_removes_only_its_own_claim(tmp_path):
    c = Claims(tmp_path)
    c.claim(0, "me")
    c.claim(1, "other")
    c.give_back(0, "me")
    c.give_back(1, "me")
    assert c.names() == {"t2"}


def test_claim_states_follow_their_owners(tmp_path, no_job):
    # MATLAB's claimStatesFollowTheirOwners: live while the owner's job is listed (running or
    # pending), or (no job) while its process runs on this node and is not the caller; failed
    # claims and claims of the caller itself are not
    c = Claims(tmp_path)
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        node = node_name()
        owners = [
            "77_1:n1:5:1",  # running job
            "78_1:n1:5:1",  # gone
            f"0_0:{node}:{sleeper.pid}:0",  # a live process here
            f"0_0:{node}:{os.getpid()}:0",  # this process
            "77_1:n1:5:1",  # marked failed below
            "79_1:n2:6:2",  # pending job
        ]
        for i, who in enumerate(owners):
            assert c.claim(i, who)
        c.mark_failed(4, "77_1:n1:5:1")
        asked = []

        def listed(keys):
            asked.append(set(keys))
            return {"77_1", "79_1"}

        claimed, live = c.states(list(range(7)), c.names(), listed)
        assert claimed == [True, True, True, True, True, True, False]
        # (a process outside Slurm is judged by /proc: Linux only, as the runners)
        here = sys.platform.startswith("linux")
        assert live == [True, False, here, False, False, True, False]
        # (one question to Slurm, with the owners' jobs)
        assert asked == [{"77_1", "78_1", "79_1"}]
    finally:
        sleeper.kill()
        sleeper.wait()
    # a process that ended
    claimed, live = c.states([2], c.names(), lambda keys: set())
    assert live == [False]


def test_unknown_slurm_answer_gives_unknown_liveness(tmp_path, no_job):
    c = Claims(tmp_path)
    c.claim(0, "77_1:n1:5:1")
    c.claim(1, "0_0:elsewhere:1:0")
    claimed, live = c.states([0, 1], c.names(), lambda keys: None)
    assert claimed == [True, True] and live == [None, False]


def test_an_owner_the_runners_did_not_write_is_live(tmp_path):
    c = Claims(tmp_path)
    c.claim(0, "someone")
    assert c.states([0], c.names(), lambda keys: set())[1] == [True]


def test_the_callers_own_job_is_judged_by_its_process(tmp_path, monkeypatch):
    # (the master in a Slurm job: its own claims are not live to it; another process of its
    # job on this node is)
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "500")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "1")
    monkeypatch.setenv("SLURMD_NODENAME", "node12")
    c = Claims(tmp_path)
    c.claim(0, claims.owner(0))
    c.claim(1, f"500_1:node12:{os.getppid()}:1")
    claimed, live = c.states([0, 1], c.names(), lambda keys: pytest.fail("no squeue"))
    assert live == [False, sys.platform.startswith("linux")]


def test_stealer_marks(tmp_path):
    c = Claims(tmp_path)
    c.mark_stealer("1000_3", "1000_3:n1:9:1")
    c.mark_stealer("1000_3", "1000_3:n1:9:1")
    assert c.stealers() == {"1000_3": "1000_3:n1:9:1"}


def test_claims_need_symbolic_links(tmp_path, monkeypatch):
    c = Claims(tmp_path)
    assert c.usable()
    assert not any(n.startswith(".probe") for n in c.names())

    def refuse(*args, **kwargs):
        raise PermissionError(1, "Operation not permitted")

    monkeypatch.setattr(os, "symlink", refuse)
    assert not Claims(tmp_path / "other").usable()


CLAIMER = """
import os
import sys
import time
sys.path.insert(0, sys.argv[1])
from petakit5d.jobs.claims import Claims
c = Claims(sys.argv[2])
n, me, out = int(sys.argv[3]), sys.argv[4], sys.argv[5]
while not os.path.exists(sys.argv[6]):
    time.sleep(0.005)
order = range(n) if int(me) % 2 else range(n - 1, -1, -1)
got = [i for i in order if c.claim(i, f"0_0:test:{me}:1")]
with open(out, "w") as fh:
    fh.write("\\n".join(map(str, got)))
"""


def test_claims_are_taken_once_under_concurrency(tmp_path):
    # MATLAB's claimsAreTakenOnceUnderConcurrency: 16 processes claim 400 tasks at once
    n, procs = 400, 16
    script = tmp_path / "claimer.py"
    script.write_text(textwrap.dedent(CLAIMER))
    outs = [tmp_path / f"got_{k}.txt" for k in range(procs)]
    go = tmp_path / "go"
    ps = [
        subprocess.Popen(
            [
                sys.executable,
                str(script),
                SRC_ROOT,
                str(tmp_path),
                str(n),
                str(k),
                str(outs[k]),
                str(go),
            ]
        )
        for k in range(procs)
    ]
    time.sleep(1.0)
    go.write_text("")
    for p in ps:
        assert p.wait(timeout=120) == 0
    got = []
    for o in outs:
        text = o.read_text().strip()
        got += [int(x) for x in text.split("\n")] if text else []
    assert sorted(got) == list(range(n)), "each task should be claimed exactly once"
    # (more than one claimer got tasks)
    assert sum(1 for o in outs if o.read_text().strip()) > 1
