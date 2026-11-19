"""Tests the mutation campaign asked for (E06, docs/hardening.md): deskew_rotate."""

import numpy as np

from petakit5d.deskew_rotate import zarr as dz
from petakit5d.io import write_zarr
from petakit5d.jobs import framework

FSN = "scan_Iter_0000_CamA_ch0_stack0000_0000000msec_0000t"
KW = {
    "skew_angle": 32.45,
    "reverse": True,
    "batch_size": (1024, 16, 1024),
    "block_size": (64, 16, 64),
    "parse_cluster": False,
    "save_mip": False,
}


def test_large_file_blocks_that_all_finish_run_once(tmp_path, monkeypatch):
    fn = str(tmp_path / f"{FSN}.zarr")
    write_zarr(np.full((29, 37, 61), 100, np.uint16), fn, block_size=(8, 16, 16))
    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append(len(list(outputs)))
        return [True] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    dz.deskew_rotate_zarr(fn, 0.108, 0.3, **KW)
    assert len(calls) == 1 and calls[0] > 1
    assert (tmp_path / "DSR" / f"{FSN}.zarr").is_dir()
