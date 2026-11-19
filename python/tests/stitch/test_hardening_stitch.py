"""Tests the mutation campaign asked for (E06, docs/hardening.md): stitching."""

import os
import sqlite3
import stat
import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d.stitch import dataset as SD
from petakit5d.stitch import imagelist as IL
from petakit5d.stitch import imagelist_generate as ig
from petakit5d.testing.synth import cut_tiles, write_stitch_tiles

POSIX = pytest.mark.skipif(os.name == "nt", reason="POSIX file modes")

SLICE_HEADER = (
    "Filename,Channel,Is_CamB_,XFPGAPosition_um_,YFPGAPosition_um_,ZFPGAPosition_um_,"
    "XMCS2Position_um_,YMCS2Position_um_,ZMCS2Position_um_"
)


def _slice_lists(root: Path, names, mcs2, nslice=3, step=0.2):
    """SliceList csv files of tiles whose stage moves in x (as test_imagelist_generate's)"""
    for i, (fn, (x0, y, z)) in enumerate(zip(names, mcs2, strict=True)):
        lines = [SLICE_HEADER]
        for k in range(nslice):
            x = x0 + step * k
            lines.append(f"{fn},0,0,{x + 0.2:f},{y:f},{z:f},{x:f},{y:f},{z:f}")
        (root / f"SliceList_tile{i + 1}_part{i:04d}.csv").write_text(
            "\n".join(lines) + "\n", newline="\n"
        )


NAME = (
    "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_{abst:010d}msecAbs_"
    "000{X}_{y:03d}{Y}_000{Z}_0000{T}.tif"
)
MCS2 = [[490.05, -8052.12, -5473.28], [490.04, -7862.72, -5473.27]]


def test_encoder_tile_indices_are_read_ignoring_case(tmp_path):
    # (MATLAB's regexpi: _000X_001Y_000Z_0000T names the same tile as _000x_001y_000z_0000t)
    names = [NAME.format(abst=1000 * (i + 1), y=i, X="X", Y="Y", Z="Z", T="T") for i in range(2)]
    _slice_lists(tmp_path, names, MCS2)
    ig.generate_image_list(str(tmp_path), "encoder", dz=0.2, channel_patterns=("CamA_ch0",))
    assert (tmp_path / "ImageList_from_encoder.csv").is_file()


@POSIX
def test_the_encoders_step_file_is_group_writable(tmp_path):
    names = [NAME.format(abst=1000 * (i + 1), y=i, X="x", Y="y", Z="z", T="t") for i in range(2)]
    _slice_lists(tmp_path, names, MCS2)
    old = os.umask(0o022)
    try:
        ig.generate_image_list(str(tmp_path), "encoder", dz=0.2, channel_patterns=("CamA_ch0",))
    finally:
        os.umask(old)
    (txt,) = list(tmp_path.rglob("dz_actual_*.txt"))
    assert os.stat(txt).st_mode & stat.S_IWGRP


DB_COLUMNS = (
    ("ImageList_id", "INTEGER"),
    ("Filename", "TEXT"),
    ("Camera", "INTEGER"),
    ("Channel", "INTEGER"),
    ("Laser_nm", "INTEGER"),
    ("File_Stack_Index", "INTEGER"),
    ("Elapsed_ms", "REAL"),
    ("Absolute_ms", "INTEGER"),
    ("Tile_index_X", "INTEGER"),
    ("Tile_index_Y", "INTEGER"),
    ("Tile_index_Z", "INTEGER"),
    ("Tile_index_T", "INTEGER"),
    ("Frame_Type", "TEXT"),
) + tuple((f"Mean_{a}_{k}_Position_um", "REAL") for k in ("MCS2", "Target", "FPGA") for a in "XYZ")


def _db(path: Path, rows) -> None:
    con = sqlite3.connect(path)
    cols = ", ".join(f"{n} {t}" for n, t in DB_COLUMNS)
    con.execute(f"CREATE TABLE ImageListView_NoNull ({cols})")
    con.executemany(
        f"INSERT INTO ImageListView_NoNull VALUES ({', '.join('?' * len(DB_COLUMNS))})", rows
    )
    con.commit()
    con.close()


def test_an_empty_acquisition_database_is_skipped_with_a_warning(tmp_path):
    rows = []
    for x in (0, 1):
        pos = [100.0 * x + 0.125, -20.5, 3.0 + 0.1 * x, 101.0 * x, -21.0, 3.5, 0.0, 0, 0]
        rows.append((x + 1, f"acq_CamA_ch0_488nm_{x:03d}x.tif", 0, 0, 488, 0, 0.0, 5000 + x, x,
                     0, 0, 0, "", *pos))  # fmt: skip
    _db(tmp_path / "Scan_Iter_0000_slicelist.sqlite3", rows)
    _db(tmp_path / "Scan_Iter_0001_slicelist.sqlite3", [])
    with pytest.warns(Warning, match=r"Scan_Iter_0001_slicelist\.sqlite3 is empty, skip it!"):
        out = ig.generate_image_list_from_sqlite(str(tmp_path), channel_patterns=("CamA",))
    assert os.path.isfile(out)


def test_a_primary_channel_without_a_camera_needs_the_ch_format():
    fsn = ["Scan_Iter_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_000x_000y"]
    with pytest.raises(IL.ImageListError, match=r"format 'ch\[0-9\]'"):
        IL._primary_channel(fsn, [0], ["A"], "chan", False)


TILE = (16, 32, 36)
STITCH_KW = {
    "dsr": True, "resample_type": "isotropic", "xy_pixel_size": 0.108, "xcorr_shift": False,
    "channel_patterns": ("CamA_ch0",), "primary_ch": "CamA_ch0", "parse_cluster": False,
    "block_size": (16, 16, 16), "batch_size": (16, 32, 32),
}  # fmt: skip


def _stitch_data(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    gt = np.random.default_rng(9).integers(1000, 60000, size=(16, 60, 72), dtype=np.uint16)
    tiles = cut_tiles(gt, [(0, 0, 0), (0, 28, 0), (0, 0, 36), (0, 28, 36)], TILE)
    return str(data) + "/", write_stitch_tiles(str(data), tiles)[1]


def test_a_finished_volume_loses_its_stale_marker(tmp_path):
    data, csv = _stitch_data(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SD.stitch_dataset(data, csv, **STITCH_KW)
        (result,) = list(Path(data, "matlab_stitch").glob("*.zarr"))
        marker = result.with_suffix(".tmp")
        marker.write_text("")
        SD.stitch_dataset(data, csv, **STITCH_KW)
    assert result.is_dir()
    assert not marker.exists(), "the marker of a finished volume was left"


def test_distance_map_tasks_ask_for_ten_times_a_tile_in_single(tmp_path, monkeypatch):
    from petakit5d.jobs import framework

    data, csv = _stitch_data(tmp_path)
    real = framework.run_tasks
    mems = []

    def spy(inputs, outputs, tasks, **kw):
        if any("compute_tile_bwdist" in t.record()["func"] for t in tasks):
            mems.append(kw["mem_allocate"])
        return real(inputs, outputs, tasks, **kw)

    monkeypatch.setattr(framework, "run_tasks", spy)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SD.stitch_dataset(data, csv, **STITCH_KW)
    assert mems == [float(np.prod(TILE)) * 4 / 1024**3 * 10]
