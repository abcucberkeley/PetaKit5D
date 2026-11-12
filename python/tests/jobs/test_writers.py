"""The writer records of a result (petakit5d.jobs.writers). Ported from
tests/framework/ResultWritersTest.m (11; the names keep the MATLAB names in snake case), and
new tests for what the port adds or changes: a failed squeue in the wait, the order of
records in whole seconds (as MATLAB's dir dates them), the context manager, the record format
MATLAB writes.

A task that makes a result in parts with jobs of its own records its process, and another
task of the same result (another uuid) waits while the first may still work on it. Owners
outside Slurm are live while their process runs on this node (Linux); owners in Slurm while
their job is in the queue (a fake Slurm); none once their record says they have ended.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

from petakit5d.jobs import files, writers
from petakit5d.testing.slurm_fake import FakeSlurm, fake_slurm, no_slurm  # noqa: F401

LINUX = sys.platform.startswith("linux")
POSIX = os.name != "nt"
linux_only = pytest.mark.skipif(
    not LINUX, reason="Writers outside Slurm are checked on Linux only."
)
posix_only = pytest.mark.skipif(not POSIX, reason="The fake Slurm needs a POSIX shell.")
SLURM_VARIABLES = ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID")


@pytest.fixture(autouse=True)
def forget_squeue_answers():
    # (squeue's answers are kept per job for a minute in the process, and the tests' fake
    # job ids repeat)
    writers.reset()
    yield
    writers.reset()


@pytest.fixture
def outside_slurm(monkeypatch):
    for v in SLURM_VARIABLES:
        monkeypatch.delenv(v, raising=False)


@pytest.fixture
def start_process():
    # processes that run during the test (killed at the end)
    procs = []

    def start() -> int:
        p = subprocess.Popen(["sleep", "120"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        procs.append(p)
        return p.pid

    yield start
    for p in procs:
        p.kill()
        p.wait()


@pytest.fixture
def wait_spy(monkeypatch):
    # wait_for_files replaced: it records its calls and returns at once (then runs action)
    calls = []
    actions = {}

    def spy(groups, wait_time):
        calls.append((groups, wait_time))
        if len(calls) in actions:
            actions[len(calls)]()
        return False

    monkeypatch.setattr(files, "wait_for_files", spy)
    spy.calls = calls
    spy.at_call = actions
    return spy


def make_result_folder(tmp_path: Path) -> tuple[str, str]:
    folder = tmp_path.as_posix() + "/"
    os.makedirs(folder + "zarr_flag/")
    return folder + "result.zarr", folder + "zarr_flag/"


def write_record(flag_root, result_fn, uuid, owner, age=0.0, extra=""):
    # the record of uuid's writer, as result_writers writes it, dated age seconds from now,
    # with extra lines (e.g. 'ended=1')
    fsn = files.fileparts(result_fn)[1]
    flag_path = f"{flag_root}{fsn}_{uuid}/"
    os.makedirs(flag_path, exist_ok=True)
    txt = f"result={result_fn}\nuuid={uuid}\nowner={owner}\n" + (f"{extra}\n" if extra else "")
    fn = flag_path + "writer.txt"
    with open(fn, "w", newline="") as f:
        f.write(txt)
    if age:
        st = os.stat(fn)
        os.utime(fn, ns=(st.st_atime_ns, st.st_mtime_ns + int(age * 1e9)))
    return flag_path


def node() -> str:
    return writers.node_name()


def process_start(pid: int) -> str:
    return writers.process_start(pid)


# ------------------------------------------------------------------------------------------
# ResultWritersTest
# ------------------------------------------------------------------------------------------


@linux_only
def test_live_writers_are_those_of_other_uuids(tmp_path, start_process, outside_slurm):
    result_fn, flag_root = make_result_folder(tmp_path)
    pid = start_process()
    n = node()
    write_record(flag_root, result_fn, "aaa", f"0_:{n}:{pid}")
    write_record(flag_root, result_fn, "bbb", f"0_:{n}:999999999")
    write_record(flag_root, result_fn, "ccc", f"0_:othernode:{pid}")
    write_record(flag_root, result_fn[:-5] + "_x.zarr", "ddd", f"0_:{n}:{pid}")
    write_record(flag_root, result_fn, "eee", f"0_:{n}:{pid}:1")
    fff = f"0_:{n}:{pid}:{process_start(pid)}"
    write_record(flag_root, result_fn, "fff", fff)
    write_record(flag_root, result_fn, "ggg", fff, 0, "ended=1")
    # live: aaa, and fff (its process, started at the time recorded); bbb's process ended,
    # ccc runs elsewhere (not checkable), ddd is of another result, eee's process id is used
    # by another process now (started at another time), ggg's writer has ended (its process
    # stays, e.g. an interactive session)
    assert sorted(writers.live(result_fn, "zzz")) == sorted([f"0_:{n}:{pid}", fff])
    # (and not one's own)
    assert writers.live(result_fn, "aaa") == [fff]


@posix_only
def test_owner_jobs_are_live_while_in_the_queue(tmp_path, fake_slurm):  # noqa: F811
    result_fn, flag_root = make_result_folder(tmp_path)
    fake_slurm.set_job_states("777_3", ["PD", "done"])
    fake_slurm.set_job_states("778", ["R"])
    write_record(flag_root, result_fn, "aaa", "777_3:node1:11")
    write_record(flag_root, result_fn, "bbb", "778_:node2:12")
    assert sorted(writers.live(result_fn, "zzz")) == ["777_3:node1:11", "778_:node2:12"]
    assert writers.live(result_fn, "zzz") == ["778_:node2:12"], "A job that ended still counted."
    calls = fake_slurm.calls("squeue")
    assert any("-j 778" in c for c in calls) and not any("-j 778_" in c for c in calls), (
        "A job that is not an array job was checked as one."
    )


@linux_only
def test_writers_of_this_processs_job_are_its_processes(
    tmp_path,
    fake_slurm,  # noqa: F811
    monkeypatch,
    start_process,
    wait_spy,
):
    # a record of another uuid owned by this process's own Slurm job (a call this process
    # made before and stopped, or a job started again with the same id): squeue lists the
    # job, so this process waited for itself. Now the process counts, as outside Slurm: this
    # one is not waited for, another of the job on this node is while it runs; squeue is not
    # asked about the job.
    result_fn, flag_root = make_result_folder(tmp_path)
    fake_slurm.set_job_states("4242", ["R"])
    monkeypatch.setenv("SLURM_JOB_ID", "4242")
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "")
    n = node()
    me = os.getpid()
    write_record(flag_root, result_fn, "aaa", f"4242_:{n}:{me}:{process_start(me)}", -60)
    os.makedirs(flag_root + "result_bbb/")
    writers.register(result_fn, "bbb", flag_root + "result_bbb/")
    assert not writers.wait(result_fn, "bbb"), "This process waited for itself."
    pid = start_process()
    other = f"4242_:{n}:{pid}:{process_start(pid)}"
    write_record(flag_root, result_fn, "ccc", other, -60)
    assert writers.live(result_fn, "bbb") == [other]
    assert not any("4242" in c for c in fake_slurm.calls("squeue")), (
        "squeue was asked about this process's job."
    )


@linux_only
def test_wait_ends_when_the_result_comes(tmp_path, start_process, wait_spy, outside_slurm, capsys):
    # a writer recorded before waits; the wait looks again each time wait_for_files returns
    # (a spy that makes the result at its second call)
    result_fn, flag_root = make_result_folder(tmp_path)
    pid = start_process()
    write_record(flag_root, result_fn, "aaa", f"0_:{node()}:{pid}", -60)
    flag_path = flag_root + "result_bbb/"
    os.makedirs(flag_path)
    writers.register(result_fn, "bbb", flag_path)
    wait_spy.at_call[2] = lambda: os.mkdir(result_fn)
    assert writers.wait(result_fn, "bbb")
    assert len(wait_spy.calls) == 2
    assert wait_spy.calls[0] == ([[result_fn]], writers.WAIT_INTERVAL)
    assert "wait for it" in capsys.readouterr().out


@linux_only
def test_writer_recorded_later_does_not_wait_for_an_earlier_one(
    tmp_path, start_process, wait_spy, outside_slurm
):
    # two writers: the later one waits for the earlier one, not the other way round (else both
    # would wait), and of two recorded at the same time, the one with the larger uuid
    result_fn, flag_root = make_result_folder(tmp_path)
    pid = start_process()
    os.makedirs(flag_root + "result_aaa/")
    writers.register(result_fn, "aaa", flag_root + "result_aaa/")
    write_record(flag_root, result_fn, "bbb", f"0_:{node()}:{pid}", 60)
    # aaa came first: it doesn't wait for bbb
    assert not writers.wait(result_fn, "aaa")
    # the same time: aaa < bbb, so aaa goes first
    aaa_fn = flag_root + "result_aaa/writer.txt"
    bbb_fn = flag_root + "result_bbb/writer.txt"
    st = os.stat(aaa_fn)
    os.utime(bbb_fn, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert not writers.wait(result_fn, "aaa")
    assert len(writers.live(result_fn, "aaa")) == 1
    assert wait_spy.calls == []


def test_release_marks_only_this_processs_record_ended(tmp_path):
    # a writer that ends but by publishing (an error, Ctrl+C) marks its record ended, once;
    # not a record of another process, and no record again once the flags are gone (the
    # result published)
    result_fn, flag_root = make_result_folder(tmp_path)
    flag_path = flag_root + "result_aaa/"
    os.makedirs(flag_path)
    writers.register(result_fn, "aaa", flag_path)
    # (of another uuid: not this one's)
    writers.release(result_fn, "zzz", flag_path)
    assert "ended=" not in Path(flag_path, "writer.txt").read_text(), (
        "A record of another uuid was marked."
    )
    writers.release(result_fn, "aaa", flag_path)
    writers.release(result_fn, "aaa", flag_path)
    txt = Path(flag_path, "writer.txt").read_bytes().decode()
    assert txt.endswith("\nended=1\n"), txt
    assert txt.count("ended=") == 1
    write_record(flag_root, result_fn, "bbb", "777_3:node1:11")
    before = Path(flag_root, "result_bbb", "writer.txt").read_bytes()
    writers.release(result_fn, "bbb", flag_root + "result_bbb/")
    assert Path(flag_root, "result_bbb", "writer.txt").read_bytes() == before, (
        "Another process's record was changed."
    )
    for f in Path(flag_path).iterdir():
        f.unlink()
    os.rmdir(flag_path)
    writers.release(result_fn, "aaa", flag_path)
    assert not os.path.isdir(flag_path), "The flags of a published result were made again."


@posix_only
def test_blocks_of_a_writer_whose_job_ended_are_orphans(tmp_path, fake_slurm, monkeypatch):  # noqa: F811
    # a block's writer whose Slurm job squeue no longer lists (e.g., cancelled), or whose
    # record says it has ended, is gone: nothing will use the block. A writer still in the
    # queue, one outside Slurm, one that is this process's own job, no record, and a failing
    # squeue: not orphans.
    result_fn, flag_root = make_result_folder(tmp_path)
    fake_slurm.set_job_states("778_2", ["R"])
    monkeypatch.setenv("SLURM_JOB_ID", "4242")
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "")
    cases = [
        ("777_3:node1:11", "", True),
        ("778_2:node1:12", "", False),
        ("0_:node1:13", "", False),
        ("4242_:node1:14", "", False),
        ("778_2:node1:12", "ended=1", True),
        ("0_:node1:13", "ended=1", True),
    ]
    for i, (owner, extra, want) in enumerate(cases, 1):
        uuid = f"u{i}"
        write_record(flag_root, result_fn, uuid, owner, 0, extra)
        flag_fn = f"{flag_root}result_{uuid}/blocks_1_5.mat"
        assert writers.orphaned(flag_fn) is want, f"{owner} {extra}"
    assert not writers.orphaned(flag_root + "result_none/blocks_1_5.mat")


@posix_only
def test_failing_squeue_does_not_make_orphans(tmp_path, no_slurm):  # noqa: F811
    # (the test Slurm stubs fail, as squeue does when the controller is busy)
    result_fn, flag_root = make_result_folder(tmp_path)
    write_record(flag_root, result_fn, "u1", "777_3:node1:11")
    assert not writers.orphaned(flag_root + "result_u1/blocks_1_5.mat")


@posix_only
def test_squeue_answers_are_kept_for_a_minute(tmp_path, fake_slurm, capsys):  # noqa: F811
    # a job that runs many blocks of one writer (task runners) asked squeue at each block;
    # squeue's answer about a job is kept in the process for a minute: 2 checks make 1 call.
    # An older answer (max_age) is asked again, and a record that says its writer has ended
    # counts at once (the record is read each time).
    result_fn, flag_root = make_result_folder(tmp_path)
    fake_slurm.set_job_states("778_2", ["R"])
    write_record(flag_root, result_fn, "u1", "778_2:node1:12")
    flag_fn = flag_root + "result_u1/blocks_1_5.mat"
    assert not writers.orphaned(flag_fn)
    assert not writers.orphaned(flag_root + "result_u1/blocks_6_10.mat")
    assert fake_slurm.call_count("squeue") == 1, "squeue was asked again within the minute."
    assert not writers.orphaned(flag_fn, 0)
    assert fake_slurm.call_count("squeue") == 2, "An older answer was used."
    with open(flag_root + "result_u1/writer.txt", "a", newline="") as f:
        f.write("ended=1\n")
    assert writers.orphaned(flag_fn), "A writer whose record says it has ended counted as running."
    assert fake_slurm.call_count("squeue") == 2
    assert "(778_2:node1:12) has ended, skip the block" in capsys.readouterr().out


@posix_only
def test_failed_squeue_is_kept_as_running(tmp_path, no_slurm, monkeypatch, capsys):  # noqa: F811
    # a failed squeue (the test Slurm stubs fail) counts as running, and is kept as that for
    # a minute; then squeue is asked again
    result_fn, flag_root = make_result_folder(tmp_path)
    write_record(flag_root, result_fn, "u1", "777_3:node1:11")
    flag_fn = flag_root + "result_u1/blocks_1_5.mat"
    assert not writers.orphaned(flag_fn)
    # (now squeue answers, without the job)
    slurm = FakeSlurm(tmp_path / "fake").install()
    monkeypatch.setenv("PATH", slurm.path_with())
    assert not writers.orphaned(flag_fn), "A failed squeue was not kept as running."
    assert slurm.call_count("squeue") == 0
    assert writers.orphaned(flag_fn, 0)
    assert "(job 777_3) has ended, skip the block" in capsys.readouterr().out


def test_no_wait_without_another_writer(tmp_path, wait_spy):
    result_fn, flag_root = make_result_folder(tmp_path)
    flag_path = flag_root + "result_aaa/"
    os.makedirs(flag_path)
    writers.register(result_fn, "aaa", flag_path)
    assert not writers.wait(result_fn, "aaa")
    txt = Path(flag_path, "writer.txt").read_bytes().decode()
    assert f"result={result_fn}\nuuid=aaa\nowner=" in txt
    assert wait_spy.calls == []


# ------------------------------------------------------------------------------------------
# new
# ------------------------------------------------------------------------------------------


@posix_only
def test_a_failed_squeue_keeps_the_wait_going(tmp_path, fake_slurm, wait_spy):  # noqa: F811
    # (MATLAB counted the job as gone and computed a second time; a failed query is unknown)
    result_fn, flag_root = make_result_folder(tmp_path)
    write_record(flag_root, result_fn, "aaa", "778_:node2:12", -60)
    write_record(flag_root, result_fn, "bbb", "779_4:node2:13", -60)
    os.makedirs(flag_root + "result_ccc/")
    writers.register(result_fn, "ccc", flag_root + "result_ccc/")
    fake_slurm.fail("squeue")
    assert sorted(writers.live(result_fn, "ccc")) == ["778_:node2:12", "779_4:node2:13"]
    wait_spy.at_call[2] = lambda: fake_slurm.recover("squeue")
    # (the jobs are not in the queue once squeue answers again)
    assert not writers.wait(result_fn, "ccc")
    assert len(wait_spy.calls) == 2


@posix_only
def test_without_squeue_slurm_owners_are_gone(tmp_path, monkeypatch):
    result_fn, flag_root = make_result_folder(tmp_path)
    write_record(flag_root, result_fn, "aaa", "778_:node2:12")
    write_record(flag_root, result_fn, "bbb", "779_4:node2:13")
    empty = tmp_path / "bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    assert writers.live(result_fn, "zzz") == []


@linux_only
def test_records_are_ordered_in_whole_seconds_then_by_uuid(
    tmp_path, start_process, wait_spy, outside_slurm
):
    # MATLAB's dir dates files in whole seconds: two records in the same second are ordered by
    # uuid, so MATLAB and Python writers of the same data agree on which one goes first
    result_fn, flag_root = make_result_folder(tmp_path)
    pid = start_process()
    os.makedirs(flag_root + "result_bbb/")
    writers.register(result_fn, "bbb", flag_root + "result_bbb/")
    write_record(flag_root, result_fn, "aaa", f"0_:{node()}:{pid}")
    t = (int(time.time()) - 10) * 1_000_000_000
    os.utime(flag_root + "result_bbb/writer.txt", ns=(t, t + 200_000_000))
    os.utime(flag_root + "result_aaa/writer.txt", ns=(t, t + 700_000_000))
    # aaa is 0.5 s later, but in the same second, and aaa < bbb: bbb waits for aaa
    wait_spy.at_call[1] = lambda: os.mkdir(result_fn)
    assert writers.wait(result_fn, "bbb")
    assert len(wait_spy.calls) == 1


def test_recorded_writer_marks_its_record_ended_on_an_error(tmp_path):
    result_fn, flag_root = make_result_folder(tmp_path)
    flag_path = flag_root + "result_aaa/"
    os.makedirs(flag_path)
    with pytest.raises(RuntimeError), writers.recorded_writer(result_fn, "aaa", flag_path):
        assert "ended=" not in Path(flag_path, "writer.txt").read_text()
        raise RuntimeError("task failed")
    assert Path(flag_path, "writer.txt").read_text().endswith("\nended=1\n")


def test_record_is_matlabs(tmp_path, monkeypatch):
    # the text MATLAB's result_writers writes, with \n line ends on every system
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "1000001")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "7")
    monkeypatch.setenv("SLURMD_NODENAME", "node12")
    result_fn, flag_root = make_result_folder(tmp_path)
    flag_path = flag_root + "result_aaa/"
    os.makedirs(flag_path)
    writers.register(result_fn.replace("/result", "//result") + "/", "aaa", flag_path)
    txt = Path(flag_path, "writer.txt").read_bytes().decode()
    start = writers.process_start(os.getpid())
    assert txt == f"result={result_fn}\nuuid=aaa\nowner=1000001_7:node12:{os.getpid()}:{start}\n"
    assert re.fullmatch(r"\d+", start) if LINUX else start == ""
    assert os.listdir(flag_path) == ["writer.txt"]
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "")
    monkeypatch.setenv("SLURM_JOB_ID", "99")
    assert writers.my_owner().startswith("99_:node12:")
    monkeypatch.delenv("SLURM_JOB_ID")
    assert writers.own_job() == (0, "")


@pytest.mark.parametrize(
    ("path", "want"),
    [
        ("/a//b///c.zarr/", "/a/b/c.zarr"),
        ("//server//share/r.zarr", "//server/share/r.zarr"),
        ("rel/r.zarr", "rel/r.zarr"),
    ],
)
def test_result_paths_are_normalized_as_matlab_does(path, want):
    assert writers.normalize_result_path(path) == want


@posix_only
def test_pending_jobs_are_live(tmp_path, fake_slurm):  # noqa: F811
    # (pending, as running: the job may still start and work on the result)
    result_fn, flag_root = make_result_folder(tmp_path)
    fake_slurm.set_job_states("779", ["PD"])
    fake_slurm.set_job_states("780_2", ["PD"])
    write_record(flag_root, result_fn, "aaa", "779_:node3:13")
    write_record(flag_root, result_fn, "bbb", "780_2:node3:14")
    assert sorted(writers.live(result_fn, "zzz")) == ["779_:node3:13", "780_2:node3:14"]


def test_refresh_listing_leaves_nothing_behind(tmp_path):
    writers.refresh_listing(tmp_path.as_posix())
    assert os.listdir(tmp_path) == []
    # (a folder that doesn't exist: nothing happens)
    writers.refresh_listing((tmp_path / "none").as_posix())
    assert os.listdir(tmp_path) == []


def test_no_flag_folder_no_writers(tmp_path):
    assert writers.live((tmp_path / "r.zarr").as_posix(), "aaa") == []
    assert not writers.orphaned((tmp_path / "zarr_flag" / "r_aaa" / "b.mat").as_posix())
