"""The hash of a folder input (a case's key) walks its folders in the order of their '/'
paths on every system: Windows' '\\' sorts after digits and '/' before them, so the chunk
folders of a zarr with 10 or more chunks along an axis ('10' and '1\\2') were hashed in
another order there and the key changed."""

import os

import pytest

from petakit5d.testing import cases


def _tree(root):
    for a in ("0", "1", "2", "10", "11"):
        for b in ("0", "3", "12"):
            d = root / "t.zarr" / a / b
            d.mkdir(parents=True)
            (d / "0").write_bytes(f"{a}.{b}".encode())
    (root / "t.zarr" / ".zarray").write_text("{}\n", newline="\n")


def test_folders_sort_as_slash_paths_with_windows_separators(monkeypatch):
    names = ["1", "1/10", "1/10/3", "1/2", "1-a", "1/1", "10", "2"]
    want = sorted(f"C:/w/d/{n}" for n in names)
    entries = [(f"C:\\w\\d\\{n}".replace("/", "\\"), [], []) for n in names]
    monkeypatch.setattr(os, "sep", "\\")
    order = sorted(entries, key=cases._walk_order)
    monkeypatch.undo()
    assert [e[0].replace("\\", "/") for e in order] == want
    # (sorted as they are, the Windows paths come in another order)
    assert [e[0].replace("\\", "/") for e in sorted(entries)] != want


@pytest.mark.skipif(os.sep != "/", reason="the order the hash had on systems with '/' paths")
def test_keys_keep_their_order_where_paths_have_slashes(tmp_path):
    _tree(tmp_path)
    walk = list(os.walk(tmp_path))
    assert sorted(walk, key=cases._walk_order) == sorted(walk)


def test_the_hash_does_not_depend_on_the_walk_order(tmp_path, monkeypatch):
    _tree(tmp_path)
    want = cases.tree_hash(tmp_path)
    real = os.walk

    def reversed_walk(root, *a, **k):
        return list(reversed(list(real(root, *a, **k))))

    monkeypatch.setattr(cases.os, "walk", reversed_walk)
    assert cases.tree_hash(tmp_path) == want
