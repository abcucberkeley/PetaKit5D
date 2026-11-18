"""``io.zarr.read_zarr_matlab_layout``: zarr reads in MATLAB's layout ((z, x, y) C order) for the
kernels that run on it, from arrays MATLAB wrote (no transposing pass) and native ones (planes
turned in place), whole or a box, into a new array, an ``out`` array or a ``ReadBuffer``."""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.io import ReadBuffer, read_zarr
from petakit5d.io.zarr import read_zarr_matlab_layout
from petakit5d.testing import synth

SHAPE = (7, 19, 13)  # (z, y, x), no two sides equal


def _vol(dtype):
    a = np.arange(np.prod(SHAPE), dtype=np.float64).reshape(SHAPE) * 3 + 1
    return a.astype(dtype)


@pytest.fixture(params=["matlab", "native"])
def zarr_file(request, tmp_path):
    vol = _vol(np.uint16)
    p = tmp_path / f"{request.param}.zarr"
    synth.write_zarr(p, vol, (4, 8, 8), request.param)
    return str(p), vol


@pytest.mark.parametrize("box", [None, (1, 2, 3, 6, 17, 12), (0, 0, 0, 1, 19, 13)])
def test_bytes_are_read_zarr_with_the_last_two_axes_swapped(zarr_file, box):
    path, vol = zarr_file
    got = read_zarr_matlab_layout(path, input_bbox=box, num_threads=3)
    ref = vol if box is None else vol[box[0] : box[3], box[1] : box[4], box[2] : box[5]]
    assert got.flags.c_contiguous and got.ndim == 3
    np.testing.assert_array_equal(got, np.swapaxes(ref, 1, 2))
    want = read_zarr(path, input_bbox=box)
    np.testing.assert_array_equal(got, np.swapaxes(want.reshape(ref.shape), 1, 2))


def test_reads_into_out_and_a_read_buffer(zarr_file):
    path, vol = zarr_file
    out = np.full((SHAPE[0], SHAPE[2], SHAPE[1]), 7, np.uint16)
    assert read_zarr_matlab_layout(path, out=out) is out
    np.testing.assert_array_equal(out, np.swapaxes(vol, 1, 2))
    buf = ReadBuffer()
    a = read_zarr_matlab_layout(path, input_bbox=(0, 0, 0, 3, 19, 13), out=buf)
    b = read_zarr_matlab_layout(path, input_bbox=(2, 1, 1, 5, 9, 12), out=buf)
    # the second read reuses the first's memory
    assert np.shares_memory(a, b)
    np.testing.assert_array_equal(b, np.swapaxes(vol[2:5, 1:9, 1:12], 1, 2))


def test_wrong_out_and_missing_files_are_refused(zarr_file, tmp_path):
    path, _ = zarr_file
    with pytest.raises(ValueError):
        read_zarr_matlab_layout(path, out=np.zeros(SHAPE, np.uint16))  # (z, y, x), not (z, x, y)
    with pytest.raises(TypeError):
        read_zarr_matlab_layout(path, out=np.zeros((7, 13, 19), np.float32))
    with pytest.raises(FileNotFoundError):
        read_zarr_matlab_layout(str(tmp_path / "none.zarr"))


@pytest.mark.parametrize("dtype", [np.float32, np.uint8])
def test_other_classes(tmp_path, dtype):
    vol = _vol(dtype)
    for layout in ("matlab", "native"):
        p = tmp_path / f"{layout}.zarr"
        synth.write_zarr(p, vol, (4, 8, 8), layout)
        got = read_zarr_matlab_layout(str(p), input_bbox=(1, 1, 1, 7, 18, 12))
        assert got.dtype == dtype
        np.testing.assert_array_equal(got, np.swapaxes(vol[1:7, 1:18, 1:12], 1, 2))
