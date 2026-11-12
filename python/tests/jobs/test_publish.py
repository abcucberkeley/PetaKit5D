"""publish_result and keep_other_result, the publish-order spy, and the duplicate-writer race.

Ported from tests/framework/PublishResultTest.m (4; the names keep the MATLAB names in snake
case). The rest is new: files as results, the rename rules of each system (Windows emulated on
Linux and macOS), answers to repeated network requests, and races of 2 to 32 writers of one
result (petakit5d.testing.publish_spy; the two-machine race on a network file system runs on
a cluster with ``python -m petakit5d.testing.publish_spy node``).
"""

from __future__ import annotations

import errno
import os
import sys
from pathlib import Path

import pytest

from petakit5d.jobs import files, publish
from petakit5d.jobs.publish import (
    keep_other_result,
    publish_result,
    remove_path,
    rename_no_replace,
    temp_path,
)
from petakit5d.testing import publish_spy
from petakit5d.testing.publish_spy import PublishSpy, list_files, normalize

POSIX = os.name != "nt"
LINUX = sys.platform.startswith("linux")
posix_only = pytest.mark.skipif(not POSIX, reason="Windows renames as the emulation does.")


@pytest.fixture
def spy():
    with PublishSpy() as s:
        yield s


@pytest.fixture(autouse=True)
def no_pauses(monkeypatch):
    monkeypatch.setattr(publish, "REMOVE_PAUSE", 0.0)
    monkeypatch.setattr(publish, "RENAME_PAUSE", 0.0)


def make_result(tmp_path: Path, other_writer: bool) -> tuple[str, str, str]:
    # a temporary result (a folder with mine.txt) and its flags, and with other_writer the
    # result of another writer (other.txt)
    folder = tmp_path.as_posix()
    result_fn = f"{folder}/result.zarr"
    tmp_fn = f"{result_fn}_1234"
    flag_path = f"{folder}/zarr_flag/result_1234/"
    os.makedirs(tmp_fn)
    Path(tmp_fn, "mine.txt").touch()
    os.makedirs(flag_path)
    Path(flag_path, "blocks_1_1.mat").touch()
    if other_writer:
        os.makedirs(result_fn)
        Path(result_fn, "other.txt").touch()
    return tmp_fn, result_fn, flag_path


def listing(folder: str) -> list[str]:
    return sorted(os.listdir(folder))


# ------------------------------------------------------------------------------------------
# PublishResultTest
# ------------------------------------------------------------------------------------------


def test_result_is_moved_into_place(tmp_path):
    tmp_fn, result_fn, flag_path = make_result(tmp_path, False)
    assert publish_result(tmp_fn, result_fn, flag_path)
    assert listing(result_fn) == ["mine.txt"]
    assert not os.path.isdir(tmp_fn)
    assert not os.path.isdir(flag_path)


def test_result_of_another_writer_is_kept(tmp_path, capsys, spy):
    tmp_fn, result_fn, flag_path = make_result(tmp_path, True)
    assert not publish_result(tmp_fn, result_fn, flag_path)
    assert listing(result_fn) == ["other.txt"], "The other writer's result was replaced."
    assert not os.path.isdir(tmp_fn), "The temporary result was not removed."
    assert not os.path.isdir(flag_path), "The flags were not removed."
    assert "made by another job" in capsys.readouterr().out
    # (the result seen, the rename is not tried)
    assert not spy.was_moved(result_fn)


def test_result_of_a_writer_just_before_is_kept(tmp_path, monkeypatch, spy):
    # a result that comes between the check and the rename: the rename fails, nothing enters
    # the result (MATLAB's movefile put the temporary result into it, and it was removed
    # there). The stub makes the check miss the other result.
    tmp_fn, result_fn, flag_path = make_result(tmp_path, True)
    monkeypatch.setattr(files, "batch_file_exist", lambda paths, *a, **k: [False] * len(paths))
    assert not publish_result(tmp_fn, result_fn, flag_path)
    assert listing(result_fn) == ["other.txt"]
    assert not os.path.isdir(tmp_fn)
    assert not os.path.isdir(flag_path)
    assert [c.ok for c in spy.moves_to(result_fn)] == [False], "The rename was not tried."


def test_kept_result_removes_a_temporary_file(tmp_path):
    folder = tmp_path.as_posix() + "/"
    result_fn = folder + "result.tif"
    tmp_fn = folder + "result_1234.tif"
    Path(result_fn).touch()
    Path(tmp_fn).touch()
    assert keep_other_result(result_fn, tmp_fn)
    assert not os.path.isfile(tmp_fn)
    assert os.path.isfile(result_fn)
    assert not keep_other_result(folder + "none.tif", folder + "none_1234.tif")


# ------------------------------------------------------------------------------------------
# files, temporary names, guards
# ------------------------------------------------------------------------------------------


def write(path: str | Path, text: str) -> str:
    Path(path).write_text(text)
    return os.fspath(path)


def test_temp_path_is_the_final_name_and_a_uuid(tmp_path):
    final = (tmp_path / "a.zarr").as_posix()
    assert temp_path(final + "/", "1234") == final + "_1234"
    t1, t2 = temp_path(final), temp_path(final)
    assert t1 != t2 and t1.startswith(final + "_") and len(t1) == len(final) + 37


def test_a_file_is_published_and_another_writers_file_kept(tmp_path, capsys):
    result = (tmp_path / "result.tif").as_posix()
    tmp1 = write(temp_path(result, "u1"), "one")
    tmp2 = write(temp_path(result, "u2"), "two")
    assert publish_result(tmp1, result)
    assert not publish_result(tmp2, result)
    assert Path(result).read_text() == "one"
    assert listing(tmp_path.as_posix()) == ["result.tif"]
    assert "made by another job" in capsys.readouterr().out


def test_a_file_that_comes_before_the_rename_is_kept(tmp_path, monkeypatch):
    # (POSIX: rename would replace the file; the publish links instead, which fails)
    result = write(tmp_path / "result.tif", "other")
    tmp = write(temp_path(result, "u1"), "mine")
    monkeypatch.setattr(files, "batch_file_exist", lambda paths, *a, **k: [False] * len(paths))
    assert not publish_result(tmp, result)
    assert Path(result).read_text() == "other"
    assert not os.path.exists(tmp)


def test_no_temporary_result_publishes_nothing_and_removes_the_flags(tmp_path):
    result = (tmp_path / "result.zarr").as_posix()
    flag = tmp_path / "zarr_flag" / "result_u1"
    flag.mkdir(parents=True)
    assert not publish_result(result + "_u1", result, flag)
    assert not os.path.exists(result)
    assert not flag.exists()


def test_a_flag_path_that_is_a_file_is_left(tmp_path):
    # (MATLAB removed the flag path only if it was a folder)
    tmp, result, _ = make_result(tmp_path, False)
    flag = write(tmp_path / "flag.txt", "x")
    assert publish_result(tmp, result, flag)
    assert os.path.isfile(flag)


@pytest.mark.parametrize("bad", ["same", "parent"])
def test_the_temporary_result_is_never_the_result(tmp_path, bad):
    tmp, result, flag = make_result(tmp_path, True)
    t = result if bad == "same" else tmp_path.as_posix()
    with pytest.raises(ValueError, match="is the result"):
        publish_result(t, result, flag)
    with pytest.raises(ValueError, match="is the result"):
        keep_other_result(result, t + "/")
    assert listing(result) == ["other.txt"]


def test_rename_no_replace_onto_a_folder_with_files_keeps_both(tmp_path):
    tmp, result, _ = make_result(tmp_path, True)
    assert not rename_no_replace(tmp, result)
    assert listing(result) == ["other.txt"]
    assert listing(tmp) == ["mine.txt"]


def test_rename_no_replace_raises_other_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        rename_no_replace(tmp_path / "none", tmp_path / "result")
    src = write(tmp_path / "a.tif", "x")
    with pytest.raises(FileNotFoundError):
        rename_no_replace(src, tmp_path / "no_folder" / "a.tif")
    assert os.path.isfile(src)


def test_failed_removal_of_the_flags_does_not_fail_a_published_result(
    tmp_path, monkeypatch, capsys
):
    tmp, result, flag = make_result(tmp_path, False)

    def fail(path):
        raise OSError(errno.ENOTEMPTY, "Directory not empty", path)

    monkeypatch.setattr(publish, "_rmtree", fail)
    assert publish_result(tmp, result, flag)
    assert listing(result) == ["mine.txt"]
    assert os.path.isdir(flag)
    assert "Could not remove" in capsys.readouterr().err


def test_remove_path_takes_entries_another_process_removed(tmp_path, monkeypatch):
    folder = tmp_path / "t"
    (folder / "sub").mkdir(parents=True)
    write(folder / "sub" / "a", "x")
    real = os.unlink

    def unlink_twice(p, *a, **k):
        real(p, *a, **k)
        raise FileNotFoundError(errno.ENOENT, "gone", p)

    monkeypatch.setattr(os, "unlink", unlink_twice)
    assert remove_path(folder)
    assert not folder.exists()
    assert remove_path(folder)


# ------------------------------------------------------------------------------------------
# rename rules of the systems, and repeated requests on network file systems
# ------------------------------------------------------------------------------------------


@posix_only
@pytest.mark.parametrize("kind", ["folder", "file"])
def test_windows_rename_rules(tmp_path, monkeypatch, kind):
    publish_spy.emulate_windows(monkeypatch)
    result = (tmp_path / ("r.zarr" if kind == "folder" else "r.tif")).as_posix()
    tmps = []
    for u in ("u1", "u2"):
        t = temp_path(result, u)
        if kind == "folder":
            os.mkdir(t)
            write(f"{t}/{u}", u)
        else:
            write(t, u)
        tmps.append(t)
    with pytest.raises(FileExistsError):
        publish_spy.windows_rename(tmps[0], tmps[1])
    monkeypatch.setattr(files, "batch_file_exist", lambda paths, *a, **k: [False] * len(paths))
    assert publish_result(tmps[0], result)
    assert not publish_result(tmps[1], result)
    content = listing(result) if kind == "folder" else [Path(result).read_text()]
    assert content == ["u1"]
    assert listing(tmp_path.as_posix()) == [os.path.basename(result)]


@posix_only
def test_windows_rename_held_by_a_scanner_is_tried_again(tmp_path, monkeypatch):
    publish_spy.emulate_windows(monkeypatch)
    tmp, result, flag = make_result(tmp_path, False)
    calls = []

    def busy_twice(src, dst):
        calls.append(src)
        if len(calls) <= 2:
            raise PermissionError(errno.EACCES, "Access is denied", src)
        publish_spy.windows_rename(src, dst)

    monkeypatch.setattr(os, "rename", busy_twice)
    assert publish_result(tmp, result, flag)
    assert len(calls) == 3
    assert listing(result) == ["mine.txt"]


@posix_only
def test_windows_rename_that_stays_denied_raises(tmp_path, monkeypatch):
    publish_spy.emulate_windows(monkeypatch)
    tmp, result, flag = make_result(tmp_path, False)

    def denied(src, dst):
        raise PermissionError(errno.EACCES, "Access is denied", src)

    monkeypatch.setattr(os, "rename", denied)
    with pytest.raises(PermissionError):
        publish_result(tmp, result, flag)
    assert listing(tmp) == ["mine.txt"]


@posix_only
@pytest.mark.parametrize("err", [errno.ENOENT, errno.ENOTEMPTY, errno.EEXIST])
def test_a_folder_rename_answered_twice_counts_as_published(tmp_path, monkeypatch, err):
    # NFS: the rename was done, its reply lost, the request sent again and answered with an
    # error; the final name is the temporary result's folder, so this writer published
    tmp, result, flag = make_result(tmp_path, False)
    real = os.rename

    def done_then_error(src, dst):
        real(src, dst)
        raise OSError(err, os.strerror(err), src)

    monkeypatch.setattr(os, "rename", done_then_error)
    assert publish_result(tmp, result, flag)
    assert listing(result) == ["mine.txt"]


@posix_only
@pytest.mark.parametrize("err", [errno.EEXIST, errno.EIO])
def test_a_file_link_answered_twice_counts_as_published(tmp_path, monkeypatch, err):
    result = (tmp_path / "result.tif").as_posix()
    tmp = write(temp_path(result, "u1"), "mine")
    real = os.link

    def done_then_exists(src, dst, *a, **k):
        real(src, dst)
        raise OSError(err, os.strerror(err), dst)

    monkeypatch.setattr(os, "link", done_then_exists)
    assert publish_result(tmp, result)
    assert Path(result).read_text() == "mine"
    assert not os.path.exists(tmp)


@posix_only
@pytest.mark.parametrize("err", [errno.EPERM, errno.EOPNOTSUPP])
def test_without_hard_links_a_file_is_renamed(tmp_path, monkeypatch, err):
    result = (tmp_path / "result.tif").as_posix()
    tmp = write(temp_path(result, "u1"), "mine")

    def no_links(src, dst, *a, **k):
        raise OSError(err, os.strerror(err), src)

    monkeypatch.setattr(os, "link", no_links)
    assert publish_result(tmp, result)
    assert Path(result).read_text() == "mine"
    tmp2 = write(temp_path(result, "u2"), "two")
    monkeypatch.setattr(files, "batch_file_exist", lambda paths, *a, **k: [False] * len(paths))
    assert not publish_result(tmp2, result)
    assert Path(result).read_text() == "mine"


@posix_only
def test_a_file_link_to_another_file_system_raises(tmp_path, monkeypatch):
    result = (tmp_path / "result.tif").as_posix()
    tmp = write(temp_path(result, "u1"), "mine")

    def cross(src, dst, *a, **k):
        raise OSError(errno.EXDEV, os.strerror(errno.EXDEV), src)

    monkeypatch.setattr(os, "link", cross)
    with pytest.raises(OSError, match="cross-device|Invalid cross|EXDEV|link"):
        publish_result(tmp, result)
    assert os.path.isfile(tmp)


# ------------------------------------------------------------------------------------------
# the spy (MovefileSpyFixture)
# ------------------------------------------------------------------------------------------


def test_spy_records_moves_with_the_files_next_to_the_destination(tmp_path, spy):
    folder = tmp_path.as_posix()
    os.makedirs(f"{folder}/MIPs/z.zarr/0")
    write(f"{folder}/MIPs/z.zarr/.zarray", "{}")
    write(f"{folder}/MIPs/a_MIP_z.tif", "x")
    tmp, result, flag = make_result(tmp_path, False)
    assert publish_result(tmp, result, flag)
    calls = spy.moves_to(result + "//")
    assert len(calls) == 1 and calls[0].op == "rename" and calls[0].ok
    assert calls[0].src == normalize(tmp)
    # (the contents of zarr folders are left out, also of the temporary zarr)
    assert spy.files_when_moved(result) == [
        "MIPs",
        "MIPs/a_MIP_z.tif",
        "MIPs/z.zarr",
        "result.zarr_1234",
        "zarr_flag",
        "zarr_flag/result_1234",
        "zarr_flag/result_1234/blocks_1_1.mat",
    ]
    assert spy.was_moved(result)
    assert not spy.was_moved(f"{folder}/other.zarr")
    assert spy.files_when_moved(f"{folder}/other.zarr") == []


def test_spy_sees_the_publish_of_a_file_by_a_link(tmp_path, spy):
    result = (tmp_path / "result.tif").as_posix()
    tmp = write(temp_path(result, "u1"), "mine")
    assert publish_result(tmp, result)
    assert [c.op for c in spy.moves_to(result)] == (["link"] if POSIX else ["rename"])


def test_publish_order_secondary_outputs_come_first(tmp_path, spy):
    # the pattern every task follows: MIPs of the temporary result, then the one rename
    folder = tmp_path.as_posix()
    result = f"{folder}/r.zarr"
    tmp = temp_path(result, "u1")
    os.makedirs(tmp)
    os.makedirs(f"{folder}/MIPs")
    write(f"{folder}/MIPs/r_MIP_z.tif", "mip")
    publish_result(tmp, result)
    assert "MIPs/r_MIP_z.tif" in spy.files_when_moved(result)
    after = set(list_files(folder)) - set(spy.files_when_moved(result))
    assert after == {"r.zarr"}, "Something was written after the result was published."


def test_spy_is_removed_on_exit():
    real = os.rename
    with PublishSpy() as s:
        assert os.rename is not real
        assert os.rename.__wrapped__ is real
    assert os.rename is real
    assert s.calls == []


def test_normalize_and_list_files(tmp_path):
    assert normalize("a//b///c/") == "a/b/c"
    assert normalize("//server//share/") == "//server/share"
    assert normalize("C:\\data\\\\x\\") == "C:/data/x"
    assert normalize("/") == "/"
    assert list_files(tmp_path / "none") == []


# ------------------------------------------------------------------------------------------
# the race: 2 to 32 writers of one result on this machine
# ------------------------------------------------------------------------------------------


def assert_race(rr):
    assert rr.ok, "\n".join(rr.errors)
    assert len(rr.published) == 1


@pytest.mark.parametrize("n", [2, 4, 8, 16, 32])
def test_race_of_writers_of_a_folder(tmp_path, n):
    assert_race(publish_spy.race_round(tmp_path.as_posix(), n, "folder", prepare=1 + n / 8))


@pytest.mark.parametrize("n", [2, 8, 32])
def test_race_of_writers_of_a_file(tmp_path, n):
    assert_race(publish_spy.race_round(tmp_path.as_posix(), n, "file", prepare=1 + n / 8))


@posix_only
@pytest.mark.parametrize("kind", ["folder", "file"])
def test_race_with_windows_rename_rules(tmp_path, kind):
    rr = publish_spy.race_round(tmp_path.as_posix(), 8, kind, semantics="windows", prepare=2)
    assert_race(rr)


@pytest.mark.skipif(not LINUX, reason="Writers outside Slurm are checked on Linux only.")
@pytest.mark.parametrize("n", [2, 8])
def test_race_of_recorded_writers_only_the_first_computes(tmp_path, n, monkeypatch):
    for v in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID"):
        monkeypatch.delenv(v, raising=False)
    rr = publish_spy.race_round(tmp_path.as_posix(), n, "folder", mode="writers", prepare=2)
    assert_race(rr)
    assert len(rr.computed) == 1 and len(rr.waited) == n - 1
