"""The Slurm command layer against the fake Slurm (and, marked cluster, read-only against a real
one).

Ported from the MATLAB suite: CheckBatchSlurmJobsStatusTest (mapsQueueStatesToJobs,
tasksSharingAJobGetItsStatus); check_slurm_job_status's and checkSlurmCluster's behavior as the
framework tests use it (SingleJobSubmitWrapperTest/logPathsStayWhenSinfoFails). New: a failed
or timed-out query is unknown (None), never "job gone"; the command lines.
"""

import json
import os
import shlex
import time

import pytest

from petakit5d.jobs import slurm as slurm_mod
from petakit5d.jobs.shapes import time_str
from petakit5d.jobs.slurm import (
    GONE,
    PENDING,
    RUNNING,
    Slurm,
    check_slurm_cluster,
    job_key,
    parse_key_values,
    split_options,
)
from petakit5d.testing.slurm_fake import DEFAULT_PARTITION, FakeSlurm, fake_slurm, no_slurm  # noqa: F401

posix_only = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")


# --- check_batch_slurm_jobs_status


@posix_only
def test_maps_queue_states_to_jobs(fake_slurm):  # noqa: F811
    # CheckBatchSlurmJobsStatusTest/mapsQueueStatesToJobs: 1 running, 0 pending, -1 any other
    # state or not in the queue
    fake_slurm.set_job_states("100_1", ["R"])
    fake_slurm.set_job_states("100_2", ["PD"])
    fake_slurm.set_job_states("100_3", ["CG"])
    fake_slurm.set_job_states("101_1", ["R"])
    # 100_4 isn't in the queue, and job id -1 was never submitted
    status = Slurm().batch_status([100, 100, 100, 100, 101, -1], [1, 2, 3, 4, 1, 1])
    assert status == [1, 0, -1, -1, 1, -1]


@posix_only
def test_tasks_sharing_a_job_get_its_status(fake_slurm):  # noqa: F811
    # CheckBatchSlurmJobsStatusTest/tasksSharingAJobGetItsStatus: one squeue call for all
    fake_slurm.set_job_states("200_1", ["R"])
    fake_slurm.set_job_states("200_2", ["PD"])
    assert Slurm().batch_status([200, 200, 200, 200], [1, 1, 2, 2]) == [1, 1, 0, 0]
    assert fake_slurm.call_count("squeue") == 1
    assert fake_slurm.calls("squeue") == ["-j 200_1,200_2"]


@posix_only
def test_each_job_gets_its_own_state_in_any_order(fake_slurm):  # noqa: F811
    # MATLAB gave [1 0 -1] here (the states in the order of the sorted job keys when no key
    # repeats: a MATLAB bug, oracle run r1); each job gets its own state
    fake_slurm.set_job_states("101_1", ["PD"])
    fake_slurm.set_job_states("100_1", ["R"])
    fake_slurm.set_job_states("102_1", ["done"])
    assert Slurm().batch_status([101, 100, 102], [1, 1, 1]) == [PENDING, RUNNING, GONE]


@posix_only
def test_jobs_without_array_ids(fake_slurm):  # noqa: F811
    fake_slurm.set_job_states("300", ["R"])
    fake_slurm.set_job_states("301", ["PD"])
    assert Slurm().batch_status([301, 300, 0]) == [0, 1, -1]
    assert fake_slurm.calls("squeue") == ["-j 300,301"]


@posix_only
def test_many_jobs_are_queried_in_batches(fake_slurm, monkeypatch):  # noqa: F811
    monkeypatch.setattr(slurm_mod, "QUERY_BATCH", 3)
    for a in range(1, 8):
        fake_slurm.set_job_states(f"400_{a}", ["R" if a % 2 else "PD"])
    assert Slurm().batch_status([400] * 7, list(range(1, 8))) == [1, 0, 1, 0, 1, 0, 1]
    assert fake_slurm.calls("squeue") == [
        "-j 400_1,400_2,400_3",
        "-j 400_4,400_5,400_6",
        "-j 400_7",
    ]


@posix_only
def test_a_failed_query_is_unknown_not_gone(fake_slurm):  # noqa: F811
    fake_slurm.set_job_states("100_1", ["R"])
    fake_slurm.fail("squeue")
    s = Slurm()
    assert s.batch_status([100, -1], [1, 1]) == [None, GONE]
    assert s.job_status(100, 1) is None
    assert s.job_state("100_1") is None
    assert s.live_jobs(["100_1"]) is None
    assert s.pending_jobs() is None
    assert s.pending_reasons("someone") is None
    fake_slurm.recover("squeue")
    assert s.batch_status([100], [1]) == [RUNNING]


@posix_only
def test_a_query_that_times_out_is_unknown(fake_slurm, monkeypatch):  # noqa: F811
    fake_slurm.set_job_states("100_1", ["R"])
    fake_slurm.hang("squeue", 600)
    s = Slurm(query_timeout=0.5)
    t0 = time.monotonic()
    assert s.batch_status([100], [1]) == [None]
    assert s.job_status(100, 1) is None
    assert time.monotonic() - t0 < 60
    # (the command gets SIGTERM at the timeout, so the call returns at once: not after the
    # wait it leaves a command that ignores SIGTERM, here 120 s, or the command's 600 s; the
    # bound keeps a busy machine's slowness apart from either)
    monkeypatch.setattr(slurm_mod, "TERMINATE_WAIT", 120.0)
    t0 = time.monotonic()
    r = s.run(["squeue"])
    assert time.monotonic() - t0 < 60
    assert r.timed_out and r.returncode is None and not r.answered


@posix_only
def test_squeue_lines_are_parsed_as_matlab_does(tmp_path, monkeypatch):
    # MATLAB drops the first line (the header), takes the last state column of a line (greedy:
    # a job name may hold " R "), and test-only counts only with exit 0
    stubs = tmp_path / "stubs"
    stubs.mkdir()
    (stubs / "squeue").write_text(
        "#!/bin/sh\n"
        'echo "             9_1       q      wrap     user  R       0:01      1 node12"\n'
        'echo "             9_2       q     a R b     user  PD       0:00      1 (Priority)"\n'
        'echo "             9_3       q      wrap     user  R       0:01      1 node12"\n'
    )
    (stubs / "sbatch").write_text(
        "#!/bin/sh\necho 'sbatch: Job 1 to start at 2026-01-01T00:00:00 a using 1 processors' >&2\n"
        "exit 1\n"
    )
    for f in stubs.iterdir():
        f.chmod(0o755)
    monkeypatch.setenv("PATH", f"{stubs}{os.pathsep}{os.environ['PATH']}")
    s = Slurm()
    assert s.batch_status([9, 9, 9], [1, 2, 3]) == [GONE, PENDING, RUNNING]
    assert s.test_only(cpus_per_task=1).start is None


def test_query_batch_size_is_matlabs():
    assert slurm_mod.QUERY_BATCH == 5000


@posix_only
def test_an_invalid_job_id_means_the_job_is_gone(fake_slurm):  # noqa: F811
    # (squeue -j <one unknown id> exits 1 with this message)
    fake_slurm.fail("squeue", "slurm_load_jobs error: Invalid job id specified")
    s = Slurm()
    assert s.batch_status([100, 101], [1, 1]) == [GONE, GONE]
    assert s.job_status(100, 1) == GONE
    assert s.job_state("100_1") == ""
    assert s.live_jobs(["100_1"]) == set()


@posix_only
def test_missing_slurm_commands_are_unknown(no_slurm, tmp_path):  # noqa: F811
    # (the no-Slurm stubs exit 1 without an answer)
    s = Slurm()
    assert s.available() is False
    assert s.batch_status([100], [1]) == [None]
    r = Slurm(env={"PATH": str(tmp_path)}).run(["squeue"])
    assert r.error and r.returncode is None and not r.answered


# --- check_slurm_job_status


@posix_only
def test_single_job_status(fake_slurm):  # noqa: F811
    fake_slurm.set_job_states("500_1", ["PD", "R", "CG", "done"])
    s = Slurm()
    assert [s.job_status(500, 1) for _ in range(4)] == [PENDING, RUNNING, PENDING, GONE]
    fake_slurm.set_job_states("501", ["R"])
    assert s.job_status(501) == RUNNING
    assert fake_slurm.calls("squeue")[-1] == "-j 501"


# --- the other squeue forms


@posix_only
def test_job_state_live_jobs_and_pending_reasons(fake_slurm):  # noqa: F811
    fake_slurm.set_job_states("600_1", ["R"])
    fake_slurm.set_job_states("600_2", ["PD:Priority"])
    fake_slurm.set_job_states("600_3", ["PD:None"])
    s = Slurm()
    assert s.job_state("600_1") == "RUNNING"
    assert s.job_state("600_9") == ""
    assert s.live_jobs(["600_1", "600_2", "600_9"]) == {"600_1", "600_2"}
    assert s.pending_reasons("someone") == {"600_2": "Priority", "600_3": "None"}
    assert fake_slurm.calls("squeue")[-3:] == [
        "-h -j 600_9 -o %T",
        "-h -r -o %i -j 600_1,600_2,600_9",
        "-r -h -u someone -t PENDING -o %i %r",
    ]


def _submit(s, array_id, cpus=1, param="", *, out="o", err="e", wrap="true"):
    """sbatch with the single-job wrapper's command line (the cluster loop adds --comment)"""
    return s.sbatch(
        ["sbatch", f"--array={array_id}", "-o", out, "-e", err, f"--cpus-per-task={cpus}",
         *split_options(param), f"--wrap={wrap}"]
    )  # fmt: skip


@posix_only
def test_pending_jobs_of_the_user(fake_slurm):  # noqa: F811
    s = Slurm()
    r = _submit(s, 3, 4, "-p gpu8 --qos q2 -n1 --mem-per-cpu=7800M --gres=gpu:1")
    fake_slurm.set_job_states(f"{r.job_id}_3", ["PD:Priority"])
    fake_slurm.set_job_states("700_1", ["R"])
    rows = s.pending_jobs()
    assert rows == [
        [
            f"{r.job_id}_[3]",
            "gpu8",
            "q2",
            "4",
            "cpu=4,mem=31200M,node=1,billing=4,gres/gpu=1",
            "Priority",
        ]
    ]
    assert (
        fake_slurm.calls("squeue")[-1]
        == "--me -t PD -h -O JobArrayID:|,Partition:|,QOS:|,NumCPUs:|,tres-alloc:|,Reason:|"
    )


# --- sbatch


@posix_only
def test_sbatch_issues_the_command_line_as_given(fake_slurm):  # noqa: F811
    s = Slurm()
    param = "-p cpu_queue --qos normal -n1 --mem-per-cpu=21418M"
    r = _submit(
        s,
        7,
        4,
        f"{param} --constraint=c1 {time_str(12, param)}",
        out="/logs/job_%A_%a.out",
        err="/logs/job_%A_%a.err",
        wrap='printf "Start time: "; python -m petakit5d.jobs.exec /calls/c1 7',
    )
    assert r.submitted and r.job_id == 1000
    assert fake_slurm.calls("sbatch") == [
        "--array=7 -o /logs/job_%A_%a.out -e /logs/job_%A_%a.err --cpus-per-task=4 -p cpu_queue "
        '--qos normal -n1 --mem-per-cpu=21418M --constraint=c1 -t 12:00:00 --wrap=printf "Start '
        'time: "; python -m petakit5d.jobs.exec /calls/c1 7'
    ]
    assert fake_slurm.job_attributes("1000_7") == {
        "partition": "cpu_queue", "qos": "normal", "cpus": 4, "mem": 4 * 21418, "gpus": 0,
    }  # fmt: skip
    r = _submit(s, 8)
    assert r.job_id == 1001
    assert fake_slurm.job_attributes("1001_8")["partition"] == DEFAULT_PARTITION


@posix_only
def test_a_refused_submission_fails_and_one_without_answer_is_unknown(fake_slurm):  # noqa: F811
    # (the submissions that answer get the default timeout: with 0.5 s, the fake's sbatch on a
    # busy machine sometimes answered too late, and a late answer is unknown, by design)
    s = Slurm()
    fake_slurm.reject_submissions(1)
    r = _submit(s, 1)
    assert r.state == "failed" and r.job_id is None and "Socket timed out" in r.result.stderr
    assert not r.result.timed_out
    assert _submit(s, 1).job_id == 1000
    fake_slurm.hang("sbatch", 60)
    t0 = time.monotonic()
    r = _submit(Slurm(submit_timeout=0.5), 2)
    assert r.state == "unknown" and r.job_id is None and r.result.timed_out
    assert time.monotonic() - t0 < 20


@posix_only
def test_test_only(fake_slurm):  # noqa: F811
    s = Slurm()
    fake_slurm.set_test_only("p2", "later")
    fake_slurm.set_test_only("p3", "error")
    now = s.test_only(cpus_per_task=2, slurm_param="-p p1 -n1", time_str=" -t 1:00:00 ")
    assert now.start is not None and abs((now.start - now.start.now()).total_seconds()) < 120
    later = s.test_only(cpus_per_task=2, slurm_param="-p p2")
    assert later.start.year == 2099 and not later.refused
    refused = s.test_only(cpus_per_task=2, slurm_param="-p p3")
    assert refused.start is None and refused.refused
    assert (
        fake_slurm.calls("sbatch")[0]
        == "--test-only --cpus-per-task=2 -p p1 -n1 -t 1:00:00 --wrap=true"
    )
    assert fake_slurm.call_count("sbatch") == 3 and not (fake_slurm.folder / "next_id").exists()


# --- scancel, scontrol


@posix_only
def test_cancel(fake_slurm):  # noqa: F811
    s = Slurm()
    assert s.cancel(["100_1", "101_2"]).ok
    assert s.cancel([]).ok and fake_slurm.call_count("scancel") == 1
    s.cancel([102, 103], background=True)
    deadline = time.monotonic() + 10
    while fake_slurm.call_count("scancel") < 2 and time.monotonic() < deadline:
        time.sleep(0.05)
    assert fake_slurm.calls("scancel") == ["100_1 101_2", "102 103"]


@posix_only
def test_show_nodes_and_update(fake_slurm):  # noqa: F811
    s = Slurm()
    assert s.show_nodes() == []
    line = FakeSlurm.node_line("n1", "cpu_queue", 24, 8, 514048, 100000)
    assert " Partitions=cpu_queue " in line and " RealMemory=514048 AllocMem=100000 " in line
    fake_slurm.set_nodes([line])
    assert s.show_nodes() == [line]
    kv = parse_key_values(line)
    assert (kv["NodeName"], kv["CPUTot"], kv["AvailableFeatures"], kv["AllocTRES"]) == (
        "n1",
        "24",
        "",
        "cpu=8,mem=100000M",
    )
    r = _submit(s, 1, 4, "-p cpu_queue")
    fake_slurm.set_job_states("1000_1", ["PD"])
    assert s.update_job(
        "1000_1", "Partition=gpu8 QOS=q2 NumCPUs=2 CPUsPerTask=2 MinMemoryCPU=1000"
    ).ok
    assert fake_slurm.job_attributes("1000_1") == {
        "partition": "gpu8",
        "qos": "q2",
        "cpus": 2,
        "mem": 2000,
        "gpus": 0,
    }
    fake_slurm.set_job_states("1000_1", ["R"])
    r = s.update_job("1000_1", {"Partition": "cpu_queue"})
    assert not r.ok and "no longer pending" in r.stderr
    assert fake_slurm.calls("scontrol")[-1] == "update JobId=1000_1 Partition=cpu_queue"


@posix_only
def test_queries_that_fail_are_unknown(fake_slurm, tmp_path, monkeypatch):  # noqa: F811
    # (the fake has no sacct, and its scontrol runs the next one on the PATH for other
    # subcommands: here failing ones, so that no real Slurm answers)
    failing = tmp_path / "failing"
    failing.mkdir()
    for cmd in ("sacct", "scontrol"):
        (failing / cmd).write_text(
            "#!/bin/sh\necho 'slurm_load_jobs error: Unable to contact slurm controller' >&2\n"
            "exit 1\n"
        )
        (failing / cmd).chmod(0o755)
    monkeypatch.setenv(
        "PATH", f"{fake_slurm.folder}{os.pathsep}{failing}{os.pathsep}{os.environ['PATH']}"
    )
    s = Slurm()
    assert s.accounting([1, 2]) is None
    assert s.accounting_state(1) is None
    assert s.show_config() is None
    assert s.show_job("1_1") is None


# --- checkSlurmCluster


@posix_only
def test_check_slurm_cluster(fake_slurm, tmp_path):  # noqa: F811
    logs = tmp_path / "logs"
    logs.mkdir()
    check = check_slurm_cluster(str(tmp_path), str(logs))
    assert check.parse_cluster
    assert (
        check.job_log_fname == f"{logs}/job_%A_%a.out"
        and check.job_log_error_fname == f"{logs}/job_%A_%a.err"
    )
    # the folder as an absolute path, '..' resolved, links kept (as MATLAB's dir().folder)
    (tmp_path / "real").mkdir()
    (tmp_path / "link").symlink_to(tmp_path / "real")
    assert (
        check_slurm_cluster(str(tmp_path), str(tmp_path / "link")).job_log_dir == f"{tmp_path}/link"
    )
    assert check_slurm_cluster(str(tmp_path), f"{tmp_path}/real/../logs").job_log_dir == str(logs)
    # only the cluster
    check = check_slurm_cluster(str(tmp_path))
    assert check.parse_cluster and check.job_log_fname == ""
    check = check_slurm_cluster(str(tmp_path), "")
    assert check.parse_cluster and check.job_log_fname == "" and check.job_log_dir == ""
    # a missing folder: <data>/job_logs, writable for all
    with pytest.warns(UserWarning, match=r"use .*/job_logs as job log directory"):
        check = check_slurm_cluster(str(tmp_path / "data"), str(tmp_path / "missing"))
    assert check.job_log_dir == f"{tmp_path}/data/job_logs"
    assert os.stat(check.job_log_dir).st_mode & 0o222 == 0o222


@posix_only
def test_log_paths_stay_when_sinfo_fails(fake_slurm, tmp_path):  # noqa: F811
    # SingleJobSubmitWrapperTest/logPathsStayWhenSinfoFails: no names when sinfo fails (the
    # caller keeps its own)
    fake_slurm.fail("sinfo")
    with pytest.warns(UserWarning, match="not exist"):
        check = check_slurm_cluster(str(tmp_path), str(tmp_path))
    assert (check.parse_cluster, check.job_log_fname, check.job_log_error_fname) == (False, "", "")
    # a sinfo that doesn't answer: Slurm is there (its queries are unknown when they fail)
    fake_slurm.recover("sinfo")
    fake_slurm.hang("sinfo", 60)
    check = check_slurm_cluster(str(tmp_path), str(tmp_path), Slurm(query_timeout=0.5))
    assert check.parse_cluster and check.job_log_fname.endswith("/job_%A_%a.out")


# --- helpers


def test_job_key_and_options():
    assert job_key(12) == "12" and job_key(12, 3) == "12_3" and job_key("12", "") == "12"
    assert split_options(" -p q  --qos 'a b' -n1 ") == ["-p", "q", "--qos", "a b", "-n1"]
    assert split_options("") == [] and split_options(None) == []
    assert parse_key_values("A=1 B=(null) C=x=y") == {"A": "1", "B": "", "C": "x=y"}


# --- a real Slurm (read-only queries and sbatch --test-only; nothing is submitted)


def _cluster_config() -> dict:
    with open(os.environ["PETAKIT5D_TEST_SLURM_CONFIG"]) as f:
        return json.load(f)


@pytest.mark.cluster
def test_real_slurm_answers_queries():
    s = Slurm()
    assert s.available() is True
    # job id 1 is long gone on any cluster that ran a few jobs
    assert s.batch_status([1], [1]) == [GONE]
    assert s.batch_status([1, 2], [1, 1]) == [GONE, GONE]
    assert s.job_status(1, 1) == GONE
    assert s.job_state("1_1") == ""
    assert s.live_jobs(["1_1"]) == set()
    assert isinstance(s.pending_jobs(), list)
    nodes = s.show_nodes()
    assert nodes and all("NodeName=" in x for x in nodes)
    assert "SLURM_VERSION" in (s.show_config() or "")
    assert s.accounting(["1"]) is not None


@pytest.mark.cluster
def test_real_slurm_test_only_for_the_test_config():
    param = _cluster_config().get("SlurmParam", "")
    r = Slurm().test_only(cpus_per_task=1, slurm_param=param, time_str=" -t 0:01:00 ")
    assert r.result.answered, r.result
    assert r.start is not None, r.result.output
    assert shlex.split(param) == list(r.result.argv[3 : 3 + len(shlex.split(param))])
