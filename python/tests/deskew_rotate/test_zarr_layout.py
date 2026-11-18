"""Deskew/rotation of zarr files MATLAB wrote, read in the kernels' layout
(``io.zarr.read_zarr_matlab_layout``) without page turns: the large-file batches
(``deskew_rotate_zarr``) and the file path (``deskew_rotate_file``). An input MATLAB wrote and the
same data written natively give the same result, bit for bit (uint16 and single, batches that
don't divide the volume, flipZstack, tiff and zarr results, the cases that keep the old path).
The oracle cases in test_zarr.py and test_file.py check MATLAB-layout inputs against MATLAB."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from petakit5d.deskew_rotate import zarr as dz
from petakit5d.io import read_zarr
from petakit5d.testing import synth
from petakit5d.testing.synth import CounterRNG

RAW = (29, 37, 61)
FSN = "scan_Iter_0000_CamA_ch0_stack0000_0000000msec_0000t"


def _frame(kind: str) -> np.ndarray:
    rng = CounterRNG(3, f"dsrzarr-layout/{kind}")
    vol = rng.integers(0, 3000, RAW, dtype=np.uint16)
    return vol.astype(np.float32) if kind == "single" else vol


@pytest.mark.parametrize("kind", ["uint16", "single"])
@pytest.mark.parametrize("flip", [False, True])
def test_matlab_layout_input_gives_the_native_result(tmp_path, kind, flip):
    vol = _frame(kind)
    results = {}
    for layout in ("matlab", "native"):
        d = Path(tmp_path, layout)
        d.mkdir()
        fn = d / f"{FSN}.zarr"
        synth.write_zarr(fn, vol, (8, 16, 16), layout)
        dz.deskew_rotate_zarr(
            str(fn),
            0.108,
            0.3,
            skew_angle=32.45,
            reverse=True,
            flip_zstack=flip,
            save_16bit=kind == "uint16",
            batch_size=(1024, 16, 1024),
            block_size=(64, 16, 64),
            parse_cluster=False,
        )
        results[layout] = read_zarr(str(d / "DSR" / f"{FSN}.zarr"))
    a, b = results["matlab"], results["native"]
    assert a.dtype == b.dtype and a.shape == b.shape
    np.testing.assert_array_equal(a, b)
    assert np.count_nonzero(a) > 0


def test_large_file_batches_are_read_in_the_kernels_layout(tmp_path, monkeypatch):
    """no page turns: every batch comes from read_zarr_matlab_layout"""
    import petakit5d.io.zarr as iz
    from petakit5d.deskew_rotate import file as df

    calls = []
    real = iz.read_zarr_matlab_layout
    monkeypatch.setattr(
        iz, "read_zarr_matlab_layout", lambda *a, **k: calls.append(1) or real(*a, **k)
    )
    monkeypatch.setattr(df, "_turn_pages", lambda *a, **k: pytest.fail("a page turn"))
    fn = tmp_path / f"{FSN}.zarr"
    synth.write_zarr(fn, _frame("uint16"), (8, 16, 16), "matlab")
    dz.deskew_rotate_zarr(
        str(fn), 0.108, 0.3, reverse=True, batch_size=(1024, 16, 1024), parse_cluster=False
    )
    assert len(calls) == 3  # 37 rows in batches of 16


# the file path (XR_deskewRotateFrame) on a zarr input


def _dsr_file(folder: Path, vol, layout, **kw):
    from petakit5d.deskew_rotate import deskew_rotate_file

    folder.mkdir(parents=True)
    fn = folder / f"{FSN}.zarr"
    synth.write_zarr(fn, vol, (8, 16, 16), layout)
    deskew_rotate_file(str(fn), 0.108, 0.3, reverse=True, rotate=True, **kw)
    ext = ".zarr" if kw.get("save_zarr") else ".tif"
    from petakit5d.io import read_tiff

    out = folder / "DSR" / f"{FSN}{ext}"
    return read_zarr(str(out)) if ext == ".zarr" else read_tiff(str(out))


@pytest.mark.parametrize(
    ("kind", "kw"),
    [
        ("uint16", {"save_zarr": True}),
        ("uint16", {"save_zarr": False}),
        ("single", {"save_zarr": True, "save_16bit": False}),
        ("uint16", {"save_zarr": True, "interp_method": "cubic"}),
        ("uint16", {"save_zarr": True, "input_bbox": (0, 2, 1, 29, 30, 61)}),
        ("uint16", {"save_zarr": True, "flip_zstack": True}),
        ("uint16", {"save_zarr": False, "flip_zstack": True}),
        ("uint16", {"save_zarr": True, "input_axis_order": "zxy"}),
        ("uint16", {"save_zarr": True, "input_axis_order": "yxz"}),
    ],
)
def test_file_on_matlab_layout_zarr_gives_the_native_result(tmp_path, kind, kw):
    vol = _frame(kind)
    a = _dsr_file(tmp_path / "matlab", vol, "matlab", **kw)
    b = _dsr_file(tmp_path / "native", vol, "native", **kw)
    assert a.dtype == b.dtype and a.shape == b.shape
    np.testing.assert_array_equal(a, b)


def test_file_reads_a_matlab_layout_zarr_without_turning_it(tmp_path, monkeypatch):
    """the combined path on MATLAB's zarr: one read in the kernels' layout and no turn of the
    input (the result's turn for a tiff output stays)"""
    import petakit5d.io.zarr as iz
    from petakit5d.deskew_rotate import file as df

    reads, turns = [], []
    real_read, real_turn = iz.read_zarr_matlab_layout, df._turn_pages
    monkeypatch.setattr(
        iz, "read_zarr_matlab_layout", lambda *a, **k: reads.append(1) or real_read(*a, **k)
    )
    monkeypatch.setattr(df, "_turn_pages", lambda *a, **k: turns.append(1) or real_turn(*a, **k))
    _dsr_file(tmp_path / "m", _frame("uint16"), "matlab", save_zarr=True)
    assert (len(reads), len(turns)) == (1, 0)
    reads.clear()
    _dsr_file(tmp_path / "n", _frame("uint16"), "native", save_zarr=True)
    assert (len(reads), len(turns)) == (0, 1)
