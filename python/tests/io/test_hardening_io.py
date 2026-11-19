"""Tests the mutation campaign asked for (E06, docs/hardening.md): io."""

import json

import numpy as np
import pytest

from petakit5d.io import bioformats as bf
from petakit5d.jobs import framework


def test_bioformats_tasks_ask_for_three_times_the_largest_file(tmp_path, monkeypatch):
    d = tmp_path / "data"
    d.mkdir()
    (d / "one.nd2").write_bytes(bytes(3 * 2**20))
    (d / "two.czi").write_bytes(bytes(2**20))
    calls = []

    def spy(inputs, outputs, tasks, **kw):
        calls.append(kw["mem_allocate"])
        return [True] * len(outputs)

    monkeypatch.setattr(framework, "run_tasks", spy)
    bf.bioformats_to_tiff_or_zarr_dataset(str(d))
    assert calls == [3 * 2**20 / 2**30 * 3]


def test_a_short_block_index_leaves_out_leading_axes_and_one_outside_is_refused(tmp_path):
    from petakit5d.io.zarr import create_zarr, read_zarr, write_zarr_block

    fn = tmp_path / "b.zarr"
    create_zarr(fn, data_size=(4, 4, 4), block_size=(2, 2, 2))
    block = np.arange(1, 9, dtype=np.uint16).reshape(2, 2, 2)
    # (iy, ix): block (0, 1, 1), as MATLAB's trailing axes
    write_zarr_block(fn, (1, 1), block)
    expected = np.zeros((4, 4, 4), np.uint16)
    expected[:2, 2:, 2:] = block
    assert np.array_equal(read_zarr(fn), expected)
    # outside along one axis only: refused by its block index
    with pytest.raises(ValueError, match=r"block \[0, 0, 2\] is outside the array"):
        write_zarr_block(fn, (0, 0, 2), block)
    assert np.array_equal(read_zarr(fn), expected)


def test_the_zarr_v3_metadata_takes_numpy_sizes():
    from petakit5d.io.zarr import zarr_json_text

    shape, block, shard = (np.array(v, np.int64) for v in ([8, 8, 8], [4, 4, 4], [2, 2, 2]))
    meta = json.loads(zarr_json_text(shape, block, shard))
    assert meta["shape"] == [8, 8, 8]
    assert meta["chunk_grid"]["configuration"]["chunk_shape"] == [4, 4, 4]
