import os
import stat
import sys

import pytest

from petakit5d.core import paths

POSIX = not sys.platform.startswith("win")


@pytest.mark.parametrize(
    "p,expected",
    [
        ("", ("", "", "")),
        ("a.b.c", ("", "a.b", ".c")),
        (".bashrc", ("", "", ".bashrc")),
        ("/a", ("/", "a", "")),
        ("/", ("/", "", "")),
        ("a/b/", ("a/b", "", "")),
        ("/a/b/c.tar.gz", ("/a/b", "c.tar", ".gz")),
        ("a.", ("", "a", ".")),
        ("x/..", ("x", ".", ".")),
    ],
)
@pytest.mark.skipif(not POSIX, reason="the POSIX rules of fileparts")
def test_fileparts(p, expected):
    assert paths.fileparts(p) == expected


@pytest.mark.parametrize(
    "p,expected",
    [
        ("C:\\x\\y.z", ("C:\\x", "y", ".z")),
        ("C:y.z", ("C:", "y", ".z")),
        ("C:\\y.z", ("C:\\", "y", ".z")),
        ("\\\\srv", ("\\\\srv", "", "")),
        ("a/b.c", ("a", "b", ".c")),
        ("\\a", ("\\", "a", "")),
    ],
)
def test_fileparts_windows_rules(p, expected):
    assert paths._fileparts_pc(p) == expected


def test_sort_key_is_utf16_code_units():
    names = ["b", "B", "a", "_x", "1", "\U0001f600", "\uffff", " "]
    assert sorted(names, key=paths.matlab_sort_key) == [
        " ",
        "1",
        "B",
        "_x",
        "a",
        "b",
        "\U0001f600",
        "\uffff",
    ]


def test_trailing_separators():
    assert paths.with_trailing_sep("/a/b") == "/a/b/" == paths.with_trailing_sep("/a/b//")
    assert paths.strip_trailing_sep("/a/b/") == "/a/b"
    assert paths.strip_trailing_sep("/") == "/"


def test_dir_simplify_and_recursive(tmp_path):
    (tmp_path / "sub" / "inner").mkdir(parents=True)
    (tmp_path / "sub" / "inner" / "deep.txt").write_text("x")
    (tmp_path / "B").write_text("x")
    (tmp_path / "a").write_text("x")
    root = str(tmp_path)
    names = [e.name for e in paths.matlab_dir(root)]
    assert names == [".", "..", "B", "a", "sub"]
    assert (
        paths.simplify_path(root + "/sub/") == os.path.join(root, "sub").replace("\\", "/")
        if POSIX
        else True
    )
    assert paths.simplify_path(os.path.join(root, "a")).endswith("a")
    assert paths.simplify_path(os.path.join(root, "missing")) == paths.filesep
    files = paths.dir_recursive(root)
    assert files == [
        root + paths.filesep + "B",
        root + paths.filesep + "a",
        paths.filesep.join([root, "sub", "inner", "deep.txt"]),
    ]
    assert paths.dir_recursive(root + paths.filesep, "folder") == [
        paths.filesep.join([root, "sub"]),
        paths.filesep.join([root, "sub", "inner"]),
    ]
    assert paths.dir_recursive(os.path.join(root, "nothing")) == []


@pytest.mark.skipif(not POSIX, reason="group permissions")
def test_mkdir_recursive_group_write(tmp_path):
    target = tmp_path / "x" / "y" / "z"
    old = os.umask(0o022)
    try:
        assert paths.mkdir_recursive(str(target), group_write=True) == 1
    finally:
        os.umask(old)
    assert target.is_dir()
    assert os.stat(target).st_mode & stat.S_IWGRP
    assert paths.mkdir_recursive(str(target)) == 1  # existing folders are fine


def test_mkdir_recursive_relative_path_starts_at_the_root(monkeypatch):
    # MATLAB's strjoin({'', 'rel'}, '/') makes 'rel/sub' into '/rel', '/rel/sub' (kept)
    made = []
    monkeypatch.setattr(paths.os, "mkdir", lambda p: made.append(p))
    monkeypatch.setattr(paths.os.path, "isdir", lambda p: False)
    monkeypatch.setattr(paths, "_IS_PC", False)
    paths.mkdir_recursive("rel_f02/sub")
    assert made == ["/rel_f02", "/rel_f02/sub"]
