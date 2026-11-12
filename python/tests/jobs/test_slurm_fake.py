"""The fake Slurm (the stubs of the MATLAB suite's FakeSlurmFixture) and the no-Slurm and
real-Slurm fixtures.

Ported from the MATLAB suite: NoSlurmFixtureTest/slurmJobOfTheTestsIsHiddenAndRestored. The
fake's behavior is the specification of FakeSlurmFixture's comment, tested here; its stubs equal
the MATLAB fixture's (hashes below, of the oracle tree's stubs with its default partition and
node name replaced), and recordings of MATLAB tests that use it (CheckBatchSlurmJobsStatusTest,
the cluster-mode tests of GenericComputingFrameworkTest, SlurmAlternativesTest) replay with
identical answers (data/fake_slurm_replay/, when present).
"""

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest

from petakit5d.testing import slurm_fake
from petakit5d.testing.slurm_fake import (  # noqa: F401  (fixtures)
    SLURM_JOB_VARIABLES,
    FakeSlurm,
    fake_slurm,
    hide_slurm,
    no_slurm,
    real_slurm,
    replay,
    replay_history,
)

pytestmark = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")

# sha256 of each MATLAB stub body (FakeSlurmFixture.m of oracle tree PetaKit5D_5df18b2, the lines
# after its 4-line header, with its default partition and node name replaced by cpu_queue and
# node12)
MATLAB_STUB_SHA256 = {
    "sbatch": "e5daf5e6aeaaaf2e845bff8f66a1ea04353b3934ff313b8a9bc233fcc1a6c8b3",
    "scancel": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "scontrol": "e84d74d8e6b536782716e4f32daabde1917bd03eb2f465a9edcbd7ee84448b39",
    "sinfo": "2da7328a20bf15a6d071ad9f39244a8d9043bf3fc910a59b70d23a90a9b9b39b",
    "squeue": "a4245e3d3259f92f3d5133d9d1fe23171909e99a6132ec2262e6746dd00189ff",
}
REPLAY_DATA = Path(__file__).parent / "data" / "fake_slurm_replay"


def sh(*argv, **kw):
    return subprocess.run(list(argv), capture_output=True, text=True, **kw)


def test_stubs_are_the_matlab_fixtures():
    for cmd, digest in MATLAB_STUB_SHA256.items():
        lines = slurm_fake._STUBS[cmd]
        body = "\n".join(lines) + "\n" if lines else ""
        body = body.replace("@PARTITION@", "cpu_queue").replace("@NODE@", "node12")
        assert hashlib.sha256(body.encode()).hexdigest() == digest, cmd
        text = slurm_fake.stub_text(cmd)
        assert text.startswith(f'#!/bin/sh\nD="$(dirname "$0")"\necho "$*" >> "$D/{cmd}.log"\n')
        assert text.endswith(body or "\n")


def test_sbatch_numbers_jobs_and_everything_is_logged(fake_slurm):  # noqa: F811
    assert sh("sinfo").stdout.startswith("PARTITION AVAIL")
    assert sh("sbatch", "--array=3", "--wrap=true").stdout == "Submitted batch job 1000\n"
    assert sh("sbatch", "--array=4", "--wrap=true").stdout == "Submitted batch job 1001\n"
    sh("scancel", "1000")
    assert fake_slurm.calls("sbatch") == ["--array=3 --wrap=true", "--array=4 --wrap=true"]
    assert fake_slurm.call_count("sbatch") == 2 and fake_slurm.call_count("squeue") == 0
    assert fake_slurm.history() == [
        "sinfo ",
        "sbatch --array=3 --wrap=true",
        "sbatch --array=4 --wrap=true",
        "scancel 1000",
    ]
    fake_slurm.clear_calls()
    assert fake_slurm.history() == [] and fake_slurm.calls("sbatch") == []
    assert sh("sbatch", "--array=1", "--wrap=true").stdout == "Submitted batch job 1000\n"


def test_states_are_used_one_per_squeue_call(fake_slurm):  # noqa: F811
    fake_slurm.set_job_states("123_7", ["PD", "PD", "R"])
    sts = [sh("squeue", "-j", "123_7").stdout.splitlines()[-1].split()[4] for _ in range(4)]
    assert sts == ["PD", "PD", "R", "R"]
    fake_slurm.set_job_states("124_1", ["R", "done"])
    assert "124_1" in sh("squeue", "-j", "124_1").stdout
    assert sh("squeue", "-j", "124_1").stdout.splitlines() == [
        "             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)"
    ]


def test_squeue_formats_reasons_and_listing(fake_slurm):  # noqa: F811
    fake_slurm.set_job_states("1_1", ["PD:None", "R"])
    fake_slurm.set_job_states("1_2", ["PD"])
    # without -j: all jobs in their current state, no state used
    assert sh("squeue", "-h", "-o", "%i %t %r").stdout.splitlines() == [
        "1_1 PD None",
        "1_2 PD Resources",
    ]
    assert sh("squeue", "-h", "-o", "%i %t %r").stdout.splitlines() == [
        "1_1 PD None",
        "1_2 PD Resources",
    ]
    assert sh("squeue", "-h", "-t", "PENDING", "-o", "%i %T").stdout.splitlines() == [
        "1_1 PENDING",
        "1_2 PENDING",
    ]
    assert sh("squeue", "-h", "-j", "1_1", "-o", "%T %r").stdout == "PENDING None\n"
    assert sh("squeue", "-h", "-j", "1_1", "-o", "%T %r").stdout == "RUNNING None\n"
    line = sh("squeue", "-j", "1_1").stdout.splitlines()[1]
    assert line == "             1_1       cpu_queue     wrap     user  R       0:01      1 node12"


def test_a_job_writes_its_outputs_when_first_running(fake_slurm, tmp_path):  # noqa: F811
    outs = [tmp_path / "a.txt", tmp_path / "b.txt"]
    fake_slurm.set_job_states("9_1", ["PD", "R"], outputs=outs)
    sh("squeue", "-j", "9_1")
    assert not any(o.exists() for o in outs)
    sh("squeue", "-j", "9_1")
    assert all(o.exists() for o in outs)


def test_rejected_submissions(fake_slurm):  # noqa: F811
    fake_slurm.reject_submissions(2)
    r = [sh("sbatch", "--wrap=true") for _ in range(3)]
    assert [x.returncode for x in r] == [1, 1, 0]
    assert "Socket timed out" in r[0].stderr and r[2].stdout == "Submitted batch job 1000\n"


def test_sbatch_records_the_jobs_resources(fake_slurm):  # noqa: F811
    sh(
        "sbatch",
        "--array=2",
        "-p",
        "gpu8",
        "--qos=q2",
        "--cpus-per-task=4",
        "--mem=32G",
        "--gres=gpu:2",
        "--wrap=x",
    )
    sh("sbatch", "--array=3", "-pfast", "--cpus-per-task=2", "--mem-per-cpu=1T", "--wrap=x")
    sh("sbatch", "--array=4", "--partition", "p3", "-q", "q3", "--wrap=x")
    assert fake_slurm.job_attributes("1000_2") == {
        "partition": "gpu8",
        "qos": "q2",
        "cpus": 4,
        "mem": 32768,
        "gpus": 2,
    }
    assert fake_slurm.job_attributes("1001_3") == {
        "partition": "fast",
        "qos": "normal",
        "cpus": 2,
        "mem": 2 * 1048576,
        "gpus": 0,
    }
    assert fake_slurm.job_attributes("1002_4") == {
        "partition": "p3",
        "qos": "q3",
        "cpus": 1,
        "mem": 1000,
        "gpus": 0,
    }
    assert fake_slurm.job_attributes("9_9") == {}


def test_pending_jobs_listed_with_long_format(fake_slurm):  # noqa: F811
    fake_slurm.set_job_attributes("50_1", "gpu8", 16, 124800, 1, "q2")
    fake_slurm.set_job_states("50_1", ["PD:Priority"])
    fake_slurm.set_job_states("51_2", ["R"])
    fake_slurm.set_job_states("52", ["PD"])
    fmt = "JobArrayID:|,Partition:|,QOS:|,NumCPUs:|,tres-alloc:|,Reason:|"
    lines = sh("squeue", "--me", "-t", "PD", "-h", "-O", fmt).stdout.splitlines()
    assert lines == [
        "50_[1]|gpu8|q2|16|cpu=16,mem=124800M,node=1,billing=16,gres/gpu=1|Priority|",
        "52|cpu_queue|normal|1|cpu=1,mem=1000M,node=1,billing=1|Resources|",
    ]
    assert sh("squeue", "-h", "-O", fmt, "-p", "gpu8").stdout.splitlines() == [lines[0]]
    assert sh("squeue", "-h", "-O", "JobArrayID:|,Reason:|").stdout.splitlines()[1] == "51_2|None|"


def test_scontrol_nodes_and_updates(fake_slurm):  # noqa: F811
    assert sh("scontrol", "show", "node", "-o").stdout == ""
    a = [FakeSlurm.node_line("n1", "cpu_queue", 24, 0, 1000, 0)]
    b = [FakeSlurm.node_line("g1", "gpu8", 128, 0, 1000, 0, 8, 2, "IDLE", "gpu_big")]
    fake_slurm.set_nodes(a, b)
    assert [sh("scontrol", "show", "node", "-o").stdout for _ in range(3)] == [
        a[0] + "\n",
        b[0] + "\n",
        b[0] + "\n",
    ]
    assert (
        "Gres=gpu:8" in b[0]
        and "AllocTRES=cpu=0,mem=0M,gres/gpu=2" in b[0]
        and "AvailableFeatures=gpu_big" in b[0]
    )
    # updates of pending jobs only
    assert sh("scontrol", "update", "JobId=7_1", "Partition=x").returncode == 1
    fake_slurm.set_job_attributes("7_1", "cpu_queue", 4, 4000)
    fake_slurm.set_job_states("7_1", ["PD"])
    assert (
        sh(
            "scontrol",
            "update",
            "JobId=7_1",
            "Partition=gpu8",
            "QOS=q2",
            "NumCPUs=2",
            "MinMemoryCPU=500",
        ).returncode
        == 0
    )
    assert fake_slurm.job_attributes("7_1") == {
        "partition": "gpu8",
        "qos": "q2",
        "cpus": 2,
        "mem": 1000,
        "gpus": 0,
    }
    assert sh("scontrol", "update", "JobId=7_1", "MinMemoryNode=3000").returncode == 0
    assert fake_slurm.job_attributes("7_1")["mem"] == 3000
    fake_slurm.fail_updates()
    r = sh("scontrol", "update", "JobId=7_1", "Partition=x")
    assert r.returncode == 1 and "Invalid qos" in r.stderr
    fake_slurm.set_job_states("7_1", ["R"])
    assert "no longer pending" in sh("scontrol", "update", "JobId=7_1", "QOS=y").stderr


def test_scontrol_runs_the_next_one_for_other_subcommands(fake_slurm, tmp_path, monkeypatch):  # noqa: F811
    other = tmp_path / "other"
    other.mkdir()
    (other / "scontrol").write_text('#!/bin/sh\necho "real $*"\n')
    (other / "scontrol").chmod(0o755)
    monkeypatch.setenv(
        "PATH", f"{fake_slurm.folder}{os.pathsep}{other}{os.pathsep}/usr/bin{os.pathsep}/bin"
    )
    assert sh("scontrol", "show", "config").stdout == "real show config\n"


def test_test_only_answers_per_partition(fake_slurm):  # noqa: F811
    fake_slurm.set_test_only("p2", "later")
    fake_slurm.set_test_only("p3", "error")
    now = sh("sbatch", "--test-only", "-p", "p1", "--cpus-per-task=3", "--wrap=true")
    assert (
        now.returncode == 0 and "to start at " in now.stderr and "using 3 processors" in now.stderr
    )
    assert (
        "2099-01-01T00:00:00" in sh("sbatch", "--test-only", "--partition=p2", "--wrap=true").stderr
    )
    err = sh("sbatch", "--test-only", "-p", "p3", "--wrap=true")
    assert err.returncode == 1 and "Invalid qos specification" in err.stderr
    assert not (fake_slurm.folder / "next_id").exists()


def test_log_files_are_made_only_in_the_temporary_folders(fake_slurm, tmp_path):  # noqa: F811
    fake_slurm.make_log_files()
    sh(
        "sbatch",
        "--array=5",
        "-o",
        f"{tmp_path}/job_%A_%a.out",
        "-e",
        f"{tmp_path}/job_%j.err",
        "--wrap=true",
    )
    assert (tmp_path / "job_1000_5.out").exists() and (tmp_path / "job_1000.err").exists()
    # outside the fake's parent folder (pytest's base temp folder): nothing
    outside = fake_slurm.folder.parent.parent / f"pk5d_fake_outside_{os.getpid()}.out"
    try:
        sh("sbatch", "-o", str(outside), "--wrap=true")
        assert not outside.exists()
    finally:
        outside.unlink(missing_ok=True)


def test_hang_and_fail_python_additions(fake_slurm):  # noqa: F811
    fake_slurm.fail("squeue", "boom")
    r = sh("squeue")
    assert (r.returncode, r.stderr) == (1, "boom\n") and fake_slurm.call_count("squeue") == 1
    fake_slurm.recover("squeue")
    fake_slurm.hang("sinfo", 30)
    with pytest.raises(subprocess.TimeoutExpired):
        sh("sinfo", timeout=0.5)
    fake_slurm.recover("sinfo")
    assert sh("sinfo").returncode == 0


def test_other_defaults_for_partition_and_node(tmp_path, monkeypatch):
    fake = FakeSlurm(tmp_path / "f", partition="batch", node="c7").install()
    monkeypatch.setenv("PATH", fake.path_with())
    fake.set_job_states("1_1", ["R"])
    assert sh("squeue", "-j", "1_1").stdout.splitlines()[1].split()[1:] == [
        "batch",
        "wrap",
        "user",
        "R",
        "0:01",
        "1",
        "c7",
    ]


# --- the no-Slurm and real-Slurm fixtures


def test_no_slurm_hides_the_commands(no_slurm):  # noqa: F811
    for cmd in ("sinfo", "sbatch", "squeue", "scancel"):
        r = sh(cmd)
        assert r.returncode == 1 and r.stdout == ""


def test_slurm_job_of_the_tests_is_hidden_and_restored(tmp_path):
    # NoSlurmFixtureTest/slurmJobOfTheTestsIsHiddenAndRestored
    with pytest.MonkeyPatch.context() as outer:
        outer.setenv("SLURM_JOB_ID", "123")
        outer.setenv("SLURM_MEM_PER_NODE", "16384")
        outer.setenv("SLURMD_NODENAME", "node12")
        path = os.environ["PATH"]
        with pytest.MonkeyPatch.context() as mp:
            hide_slurm(mp, tmp_path / "stubs")
            hidden = [
                os.environ.get(v, "")
                for v in ("SLURM_JOB_ID", "SLURM_MEM_PER_NODE", "SLURMD_NODENAME")
            ]
            child = sh("sh", "-c", 'printf "%s" "${SLURM_JOB_ID-unset}"').stdout
            assert os.environ[slurm_fake.STUBS_VARIABLE] == str(tmp_path / "stubs")
        assert hidden == ["", "", ""]
        assert child == "unset"
        names = ("SLURM_JOB_ID", "SLURM_MEM_PER_NODE", "SLURMD_NODENAME")
        assert [os.environ.get(v) for v in names] == ["123", "16384", "node12"]
        assert os.environ["PATH"] == path
    assert set(SLURM_JOB_VARIABLES) >= {"SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_MEM_PER_CPU"}


def test_real_slurm_takes_the_stubs_off_the_path(no_slurm, request):  # noqa: F811
    stubs = os.environ[slurm_fake.STUBS_VARIABLE]
    assert stubs in os.environ["PATH"].split(os.pathsep)
    request.getfixturevalue("real_slurm")
    assert stubs not in os.environ["PATH"].split(os.pathsep)


# --- the MATLAB recordings


def _recordings():
    if not REPLAY_DATA.is_dir():
        return []
    return sorted(REPLAY_DATA.glob("*.json"))


@pytest.mark.skipif(not _recordings(), reason="no recordings in tests/jobs/data/fake_slurm_replay")
@pytest.mark.parametrize("fn", _recordings(), ids=lambda p: p.stem)
def test_matlab_recordings_replay_with_identical_answers(fn, tmp_path):
    rec = json.loads(fn.read_text())
    for i, seg in enumerate(rec["segments"]):
        root = tmp_path / f"s{i}"
        root.mkdir()
        assert replay(seg["events"], root) == [], seg["test"]
        assert replay_history(root) == seg["history"], seg["test"]
