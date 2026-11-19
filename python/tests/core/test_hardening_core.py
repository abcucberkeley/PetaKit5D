"""Tests the mutation campaign asked for (E06, docs/hardening.md): core."""

import json
import os

import numpy as np
import pytest

from petakit5d.core import arrays, bbox, dtypes, paths, threads


def test_python_numbers_are_doubles():
    # MATLAB's 2^60 + 1 is a double: 2^60 (an int64 plus it is 2^60)
    got = dtypes.plus(np.array([0], dtype=np.int64), 2**60 + 1)
    assert got.dtype == np.int64 and int(got[0]) == 2**60
    assert dtypes.mround(5).dtype == np.float64
    assert dtypes.mround(True).dtype == np.bool_


def test_thread_counts_are_ints():
    n = threads.resolve_threads(1.0)
    assert n == 1 and type(n) is int


def test_uint32_saturates_and_rounds_half_away():
    got = dtypes.uint32(np.array([-3.0, 2.5, 3.5, 4294967295.4, 5e9]))
    assert got.dtype == np.uint32
    np.testing.assert_array_equal(got, [0, 3, 4, 4294967295, 4294967295])


def test_sub2ind_of_vectors_gives_a_vector_and_of_scalars_an_int():
    got = arrays.sub2ind((3, 4), [1, 2], [2, 3])
    assert isinstance(got, np.ndarray)
    np.testing.assert_array_equal(got, [4, 8])
    one = arrays.sub2ind((3, 4), 2, 3)
    assert isinstance(one, int) and one == 8


@pytest.mark.skipif(os.name == "nt", reason="group write permission is POSIX")
def test_mkdir_recursive_leaves_a_folder_another_process_made(tmp_path, monkeypatch):
    # another process makes the folder after the check: no chmod of a folder that isn't ours
    target = tmp_path / "a" / "b"
    real = os.makedirs
    chmods = []

    def racing(p, *a, **k):
        real(p, *a, **k)
        raise FileExistsError(p)

    monkeypatch.setattr(paths.os, "makedirs", racing)
    monkeypatch.setattr(paths.os, "chmod", lambda p, mode: chmods.append(p))
    assert paths.mkdir_recursive(str(target), group_write=True) == 1
    assert target.is_dir()
    assert chmods == []
    monkeypatch.setattr(paths.os, "makedirs", real)
    paths.mkdir_recursive(str(tmp_path / "c"), group_write=True)
    assert len(chmods) == 1 and os.path.samefile(chmods[0], tmp_path / "c")


def test_a_full_box_of_a_numpy_shape_is_json_data():
    # (boxes go into task records, which are JSON)
    b = bbox.full_box(np.array([4, 5, 6]))
    assert b == (0, 0, 0, 4, 5, 6) and json.loads(json.dumps(b)) == [0, 0, 0, 4, 5, 6]
    assert all(type(v) is int for v in bbox.full_box(np.array([7, 8], np.uint16)))
