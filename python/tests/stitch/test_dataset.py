"""Stitching of datasets (petakit5d.stitch.dataset) against MATLAB's XR_matlab_stitching_wrapper
(the oracle tree's), and StitchingTest's wrapper tests in intent.

Oracle cases ``wrapper.<scene>``: each side gets the same tiles and image list (and, where the
scene registers, the primary volumes' registration results ``xcorr/<name>.mat``, the shifts
given, since the registration is S04's), then runs the wrapper in local mode (MATLAB's
function, the port's through its MATLAB-name layer, with the same arguments). Compared: every
file the runs wrote (the parameters files, the stitched volumes and their MIPs, the stitch and
block information, the converted tiles, the distance maps), the zarr files by voxels and
metadata after the axis mapping. Scenes: a time series of two cameras registered once
('primaryFirst': every time point on the primary's canvas), 'primary' and 'all', a pinned stitch
information file, processed DSR zarr tiles, skewed space, tiles without camera names, several
folders (multiLoc), online stitching, streaming with the tiles there, local workers; and the
demo datasets (``demo.*``, with ``PETAKIT5D_DEMO_DATA``): demo_skewed_space_stitching and the
DSR stitching of demo_zarr_stitching.

The pytest functions port StitchingTest's wrapper methods (local mode, and cluster mode against
the fake Slurm: the jobs' sizes and nice, the master's last volume, the last tries, results,
tiles and stitch information other nodes just made, the waits of a second, the pipeline's
marker, streaming, local workers, secondaries) and test the wrapper's own logic (the function
string and its values, the names, the modes, online stitching, the parameters files, the tile
information file, MATLAB's quirks).
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pytest

from petakit5d.io import read_mat, read_zarr
from petakit5d.jobs import files as jfiles
from petakit5d.jobs import single
from petakit5d.stitch import dataset as SD
from petakit5d.stitch import volume as V
from petakit5d.testing import oracle, slurm_fake
from petakit5d.testing.cases import Case, In, InPath, Tree, Value
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG, cut_tiles, stitch_tile_name, write_stitch_tiles

fake_slurm = slurm_fake.fake_slurm  # (the fixture)

POSIX = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
LINUX = pytest.mark.skipif(
    not os.sys.platform.startswith("linux"), reason="Task runners run on Linux only."
)

PX = 0.108
NAME = "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs"
GRID_STARTS = [(0, 0, 0), (0, 36, 0), (0, 0, 48), (0, 36, 48)]  # (z0, y0, x0)
GRID_GT = (32, 100, 120)
GRID_TILE = (32, 64, 72)


def name_of(it: int = 0, cam: str = "A") -> str:
    return NAME.replace("Iter_0000", f"Iter_{it:04d}").replace("CamA", f"Cam{cam}")


def random_volume(shape, seed: int = 1, stream: str = "gt", lo: int = 1000, hi: int = 60000):
    """nonzero random uint16 values (0 means no data in stitching)"""
    return CounterRNG(seed, stream).integers(lo, hi, shape).astype(np.uint16)


def write_tiles(data, gt, starts=GRID_STARTS, tile=GRID_TILE, pixel=(PX, PX, PX), stage=None,
                tidx=None, **kw):  # fmt: skip
    """the tiles of cut_tiles and their image list (StitchingTest.writeTiles)"""
    tiles = cut_tiles(gt, starts, tile, pixel, stage, tidx)
    return write_stitch_tiles(str(data), tiles, **kw)[1]


def dp(path) -> str:
    """a data folder with a trailing separator, as MATLAB's tests and demos give them (the
    wrapper's paths then have MATLAB's '//')"""
    s = str(path)
    return s if s.endswith(("/", "\\")) else s + "/"


def stitch(data, csv, **kw):
    """the wrapper with StitchingTest's defaults for exact DSR stitching (later ones win)"""
    args = {
        "dsr": True,
        "resample_type": "isotropic",
        "xy_pixel_size": PX,
        "axis_order": "x,y,z",
        "xcorr_shift": False,
        "channel_patterns": ("CamA_ch0",),
        "primary_ch": "CamA_ch0",
        "parse_cluster": False,
        **kw,
    }
    SD.stitch_dataset(dp(data), str(csv), **args)


def result(data, name: str = NAME, folder: str = "matlab_stitch") -> np.ndarray:
    return read_zarr(f"{data}/{folder}/{name}.zarr")


def two_volumes(root: Path, n: int = 2):
    """two (or n) time points of the 2 x 2 tile layout, and the paths of their results"""
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    gt = random_volume(GRID_GT)
    csv = None
    for it in range(n):
        csv = write_tiles(data, gt, iter=it)
    res = [f"{data}/matlab_stitch/{name_of(it)}.zarr" for it in range(n)]
    return data, csv, res


def tiles_to_come(root: Path):
    """the raw tiles of the 2 x 2 layout, and their processed (DSR) tiles in stage/DSR/, to
    be moved to data/DSR/"""
    data, stage = root / "data", root / "stage"
    data.mkdir(parents=True)
    stage.mkdir()
    gt = random_volume(GRID_GT)
    csv = write_tiles(data, gt)
    write_tiles(stage, gt, tile_dir="DSR")
    (data / "DSR").mkdir()
    return data, csv, stage


class FrameSpy:
    """the frame's calls (keyword arguments), each running ``action`` if given"""

    def __init__(self, action=None):
        self.calls: list[dict] = []
        self.action = action

    def __call__(self, *a, **kw):
        self.calls.append(kw)
        if self.action is not None:
            self.action(kw)


@pytest.fixture
def frames(monkeypatch):
    spy = FrameSpy()
    monkeypatch.setattr(V, "stitch_volume", spy)
    return spy


def make_result(path: str):
    """an empty file in place of a volume's result (as the MATLAB tests' spies make it)"""

    def act(_kw=None):
        with open(path, "w"):
            pass

    return act


@pytest.fixture
def pauses(monkeypatch):
    """the waits' pauses, without waiting (the fake clock moves by each pause)"""
    now = [0.0]
    calls: list[float] = []

    def pause(t):
        calls.append(t)
        now[0] += t

    monkeypatch.setattr(jfiles, "_pause", pause)
    monkeypatch.setattr(jfiles, "_clock", lambda: now[0])
    return calls


def at_call(counter: list, n, fcn):
    """count the calls; at the n-th run fcn (n a list: at the n[i]-th, run fcn[i])"""
    ns = n if isinstance(n, (list, tuple)) else [n]
    fs = fcn if isinstance(fcn, (list, tuple)) else [fcn]

    def hook(*_a):
        counter[0] += 1
        for k, f in zip(ns, fs, strict=True):
            if counter[0] == k:
                f()

    return hook


def waits_of_a_second_problem(waits, n) -> str:
    """pktest.waitsOfASecondProblem: n waits, the first a second, none of 2 s or more"""
    if len(waits) != n:
        return f"There were {len(waits)} waits instead of {n}: {waits} s."
    if n > 0 and waits[0] != 1:
        return f"The first wait was not a second: {waits} s."
    if any(w <= 0 or w >= 2 for w in waits):
        return f"A wait was not between 0 and 2 s: {waits} s."
    return ""


def cluster(data, csv, logs: Path, **kw):
    logs.mkdir(exist_ok=True)
    kw.setdefault("parse_cluster", True)
    kw.setdefault("job_log_dir", str(logs))
    stitch(data, csv, **kw)


def cpus_of(fake) -> list[int]:
    return [int(re.search(r"--cpus-per-task=(\d+)", c).group(1)) for c in fake.calls("sbatch")]


@pytest.fixture
def submitted(monkeypatch):
    """the task record of each sbatch submission, in order"""
    out: list[dict] = []
    real = single.submit_single_job

    def spy(task, *a, **kw):
        r = real(task, *a, **kw)
        if r[2] > 0:
            out.append(task.record())
        return r

    monkeypatch.setattr(single, "submit_single_job", spy)
    return out


# ---------------------------------------------------------------------------------------------
# StitchingTest in intent: local mode
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["none", "feather", "max", "mean", "median"])
def test_exact_positions_reproduce_ground_truth(tmp_path, method):
    """(exactPositionsReproduceGroundTruth) the wrapper's result is the ground truth at the
    canvas convention: the first tile at voxel 0, two empty voxels at the end of each axis"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    stitch(tmp_path, csv, blend_method=method, batch_size=(16, 32, 32), block_size=(16, 16, 16))
    out = result(tmp_path)
    assert out.shape == tuple(s + 2 for s in GRID_GT)
    np.testing.assert_array_equal(out[:32, :100, :120], gt)
    assert not out[32:].any() and not out[:, 100:].any() and not out[:, :, 120:].any()


def test_tiles_2d(tmp_path):
    """(tiles2D) single page tiles"""
    gt = random_volume((100, 120))
    starts = [(0, 0), (36, 0), (0, 48), (36, 48)]
    csv = write_tiles(tmp_path, gt, starts, (64, 72))
    stitch(tmp_path, csv, blend_method="feather")
    out = result(tmp_path)
    assert out.shape[-2:] == (102, 122)
    np.testing.assert_array_equal(out.reshape(102, 122)[:100, :120], gt)


def test_io_scan_uses_stage_positions_as_is(tmp_path):
    """(ioScanUsesStagePositionsAsIs) IO scan tiles in x and z (z offsets in units of dz)"""
    gt = random_volume((40, 64, 104))
    starts = [(0, 0, 0), (0, 0, 40), (16, 0, 0), (16, 0, 40)]
    csv = write_tiles(tmp_path, gt, starts, (24, 64, 64), (0.5, PX, PX))
    stitch(tmp_path, csv, dsr=False, io_scan=True, dz=0.5, blend_method="none")
    out = result(tmp_path)
    assert out.shape == (42, 66, 106)
    np.testing.assert_array_equal(out[:40, :64, :104], gt)


@pytest.mark.parametrize("axis_order", ["x,y,z", "-x,y,z"])
def test_skewed_space_placement(tmp_path, axis_order):
    """(skewedSpacePlacement) a stage x step of n * dz is n raw frames in z; '-x' flips it"""
    n, dz = 10, 0.3
    gt = random_volume((50, 100, 64))
    z1, z2 = (0, n) if axis_order.startswith("-") else (n, 0)
    starts = [(z1, 0, 0), (z1, 36, 0), (z2, 0, 0), (z2, 36, 0)]
    stage = [[0, 0, 0], [0, 36 * PX, 0], [n * dz, 0, 0], [n * dz, 36 * PX, 0]]
    tidx = [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
    csv = write_tiles(tmp_path, gt, starts, (40, 64, 64), stage=stage, tidx=tidx)
    stitch(tmp_path, csv, dsr=False, dz=dz, skew_angle=32.45, reverse=True,
           axis_order=axis_order, blend_method="none")  # fmt: skip
    out = result(tmp_path)
    np.testing.assert_array_equal(out[:50, :100, :64], gt)


def test_processed_zarr_tiles(tmp_path):
    """(processedZarrTiles) DSR zarr tiles of the raw names in DSR/"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt, tile_dir="DSR", save_format="zarr")
    stitch(tmp_path, csv, zarr_file=True, processed_dir_str="DSR", blend_method="feather")
    np.testing.assert_array_equal(result(tmp_path)[:32, :100, :120], gt)


def test_volume_size_estimate_matches_the_stitched_volume(tmp_path):
    """(volumeSizeEstimateMatchesTheStitchedVolume, the wrapper's part) the wrapper sizes a
    volume's job by stitch_volume_size_estimate, which is the stitched size (S03 tests its
    layouts); here through the wrapper's own call"""
    from petakit5d.stitch.geometry import stitch_volume_size_estimate

    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    stitch(tmp_path, csv, out_bbox=(1, 4, 2, 30, 90, 100))
    out = result(tmp_path)
    tiles = cut_tiles(gt, GRID_STARTS, GRID_TILE)
    est, _, rng, _ = stitch_volume_size_estimate(
        np.array([t["xyz"] for t in tiles]), GRID_TILE,
        tile_idx=np.array([[*t["idx"], 1] for t in tiles], dtype=np.float64),
        dsr=True, resample_type="isotropic", xcorr_shift=False, out_bbox=(1, 4, 2, 30, 90, 100),
    )  # fmt: skip
    assert tuple(int(v) for v in est) == out.shape
    assert tuple(int(v) for v in rng[0]) == out.shape


def test_secondary_channel_uses_the_task_size_of_the_primary(tmp_path, monkeypatch):
    """(secondaryChannelUsesTheTaskSizeOfThePrimary) the secondary reads the primary's per-task
    block information, whose task size depends on the primary's local workers"""
    from petakit5d.core import memory

    monkeypatch.setattr(memory, "get_system_memory", lambda: 0.0)
    gt_a, gt_b = random_volume((32, 64, 104), 1), random_volume((32, 64, 104), 2)
    stage = [[0, 0, 0], [40 * PX, 0, 0]]
    tidx = [[0, 0, 0], [1, 0, 0]]
    starts = [(0, 0, 0), (0, 0, 40)]
    write_tiles(tmp_path, gt_a, starts, (32, 64, 64), stage=stage, tidx=tidx, cam="A")
    csv = write_tiles(tmp_path, gt_b, starts, (32, 64, 64), stage=stage, tidx=tidx, cam="B")
    opts = {"primary_ch": "CamA_ch0", "blend_method": "feather", "batch_size": (16, 32, 32),
            "block_size": (16, 16, 16)}  # fmt: skip
    stitch(tmp_path, csv, **opts, channel_patterns=("CamA_ch0",), local_worker_num=4)
    stitch(tmp_path, csv, **opts, channel_patterns=("CamA_ch0", "CamB_ch0"), local_worker_num=1)
    out_b = result(tmp_path, name_of(0, "B"))
    np.testing.assert_array_equal(out_b[:32, :64, :104], gt_b)


def test_secondary_with_a_short_tile_stitches_what_it_has(tmp_path):
    """(secondaryWithAShortTileStitchesWhatItHas) time point 1's tile 4 is shorter in z and
    its tile 2 one voxel off the primary's: it gets its own block information, with the
    primary's tile positions and shifts, and is its tiles stitched on their own"""
    data, alone = tmp_path / "data", tmp_path / "alone"
    data.mkdir()
    alone.mkdir()
    write_tiles(data, random_volume(GRID_GT, 1))
    tiles = cut_tiles(random_volume(GRID_GT, 2), GRID_STARTS, GRID_TILE)
    tiles[3]["im"] = np.ascontiguousarray(tiles[3]["im"][:20])
    _, csv_alone = write_stitch_tiles(str(alone), tiles, iter=1)
    tiles[1]["xyz"] = [tiles[1]["xyz"][0] + PX, *tiles[1]["xyz"][1:]]
    _, csv = write_stitch_tiles(str(data), tiles, iter=1)
    with pytest.warns(V.StitchVolumeWarning) as rec:
        stitch(data, csv)
    assert any(getattr(w.message, "identifier", "") == "XR_stitching_frame_zarr_dev_v1:tileSizes"
               for w in rec)  # fmt: skip
    stitch(alone, csv_alone)
    np.testing.assert_array_equal(result(data, name_of(1)), result(alone, name_of(1)))
    info = data / "matlab_stitch" / "stitchInfo"
    assert list(info.glob(f"{name_of(1)}_task_size_*")), "The secondary has no block info."
    assert list(info.glob(f"{NAME}_task_size_*"))


def test_secondary_of_the_primary_sizes_uses_its_block_info(tmp_path):
    """(secondaryOfThePrimarySizesUsesItsBlockInfo) the same tiles and size: no block
    information of its own, the same result"""
    data, csv, _ = two_volumes(tmp_path)
    stitch(data, csv)
    info = data / "matlab_stitch" / "stitchInfo"
    assert list(info.glob(f"{NAME}_task_size_*"))
    assert not list(info.glob(f"{name_of(1)}_task_size_*")), "The secondary made its own."
    np.testing.assert_array_equal(result(data, name_of(1)), result(data))


def test_secondary_of_another_stitched_size_gets_its_own_block_info(tmp_path):
    """(secondaryOfAnotherStitchedSizeGetsItsOwnBlockInfo) a later time point stitched with
    another out box makes its own block information: its result is its box of the volume"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    stitch(tmp_path, csv, out_bbox=(0, 0, 0, 32, 100, 120))
    write_tiles(tmp_path, gt, iter=1)
    with pytest.warns(V.StitchVolumeWarning) as rec:
        stitch(tmp_path, csv, out_bbox=(0, 10, 20, 32, 90, 100))
    assert any(
        getattr(w.message, "identifier", "") == "XR_stitching_frame_zarr_dev_v1:stitchedSize"
        for w in rec
    )
    np.testing.assert_array_equal(result(tmp_path, name_of(1)), gt[:, 10:90, 20:100])
    info = tmp_path / "matlab_stitch" / "stitchInfo"
    assert list(info.glob(f"{name_of(1)}_task_size_*")), "The secondary has no block info."


def test_a_time_series_is_registered_once_on_one_canvas(tmp_path):
    """('primaryFirst', acceptance) the first time point's primary is registered (here its
    registration result is given: S04), and every time point and channel is placed with its
    shifts on the same canvas"""
    gt = random_volume(GRID_GT)
    stage = [[t["xyz"][0] + (2 * PX if k == 2 else 0), *t["xyz"][1:]]
             for k, t in enumerate(cut_tiles(gt, GRID_STARTS, GRID_TILE))]  # fmt: skip
    for it in range(3):
        for cam in "AB":
            csv = write_tiles(tmp_path, gt, stage=stage, iter=it, cam=cam)
    _write_shift_result(tmp_path, gt, stage, [[0, 0, 0], [0, 0, 0], [-2, 0, 0], [0, 0, 0]])
    stitch(tmp_path, csv, xcorr_shift=True, channel_patterns=("CamA_ch0", "CamB_ch0"),
           blend_method="none")  # fmt: skip
    outs = {(it, cam): result(tmp_path, name_of(it, cam)) for it in range(3) for cam in "AB"}
    shapes = {o.shape for o in outs.values()}
    assert len(shapes) == 1, shapes
    for o in outs.values():
        np.testing.assert_array_equal(o, outs[(0, "A")])
    # registered once: only the primary's registration result and stitch information
    rt = tmp_path / "matlab_stitch"
    assert sorted(p.name for p in (rt / "xcorr").iterdir()) == [f"{NAME}.mat"]
    infos = [p.name for p in (rt / "stitchInfo").glob("*.mat") if "_block_info" not in p.name]
    assert infos == [f"{NAME}.mat"]


def _write_shift_result(data, gt, stage, shift_vox, name: str = NAME, folder="matlab_stitch",
                        **layout):  # fmt: skip
    """the registration result of the volume (MATLAB's xcorr/<name>.mat), shifts in voxels
    (x, y, z) per tile, as S08's cases write it"""
    from petakit5d.stitch import geometry as G

    tiles = cut_tiles(gt, GRID_STARTS, GRID_TILE, stage_xyz=stage)
    lay = G.prepare_layout(
        np.array([t["xyz"] for t in tiles]),
        np.array([t["im"].shape for t in tiles]),
        tile_idx=np.array([[*t["idx"], 1] for t in tiles], dtype=np.float64),
        axis_order="x,y,z", data_order="z,y,x", xy_pixel_size=PX, dz=0.5, dsr=True,
        resample_type="isotropic", **layout,
    )  # fmt: skip
    d = np.asarray(shift_vox, dtype=np.float64)[lay.sinds] * lay.space.pixel_sizes
    out = Path(data, folder, "xcorr")
    out.mkdir(parents=True, exist_ok=True)
    V.write_xcorr_result(str(out / f"{name}.mat"), lay.xyz + d, d)


def test_volumes_with_another_tile_count_are_skipped_in_primary_first(tmp_path, capsys,
                                                                       frames):  # fmt: skip
    """('primaryFirst') a volume of another number of tiles than the primary is skipped"""
    gt = random_volume(GRID_GT)
    write_tiles(tmp_path, gt)
    csv = write_tiles(tmp_path, gt, GRID_STARTS[:3], iter=1)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv)
    assert [Path(c["result_path"]).name for c in frames.calls] == [f"{NAME}.zarr"]
    assert "different from the number of tiles in the primary channel" in capsys.readouterr().out


def test_local_workers_give_the_same_result(tmp_path):
    """(localWorkersGiveTheSameResult, without the registration: S04) the frame's tasks on 2
    local workers give the result and MIP of one task at a time"""
    gt = random_volume(GRID_GT)
    outs = []
    for k in (2, 1):
        d = tmp_path / f"w{k}"
        d.mkdir()
        csv = write_tiles(d, gt)
        stitch(d, csv, blend_method="feather", batch_size=(16, 32, 32), block_size=(16, 16, 16),
               local_worker_num=k)  # fmt: skip
        outs.append((result(d), _read_mip(d)))
    np.testing.assert_array_equal(outs[0][0], outs[1][0])
    np.testing.assert_array_equal(outs[0][1], outs[1][1])


def _read_mip(d, name: str = NAME, folder: str = "matlab_stitch"):
    from petakit5d.io import read_tiff

    return read_tiff(f"{d}/{folder}/MIPs/{name}_MIP_z.tif")


def test_large_datasets_use_local_workers(tmp_path, frames, pauses, monkeypatch):
    """(largeDatasetsUseLocalWorkers) without local_worker_num, volumes with 32 or more tiles
    to register use local workers (one per 3 CPUs, up to 8), smaller ones none"""
    from petakit5d.core import threads

    monkeypatch.setattr(threads, "available_cpus", lambda: 12)
    for n_t in (4, 32):
        d = tmp_path / f"n{n_t}"
        d.mkdir()
        x, y = np.meshgrid(np.arange(4), np.arange(n_t // 4), indexing="ij")
        x, y = x.ravel(order="F"), y.ravel(order="F")
        tiles = [
            {
                "im": np.ones((2, 8, 8), np.uint16),
                "idx": [int(a), int(b), 0],
                "xyz": [a * 6 * PX, b * 6 * PX, 0.0],
            }
            for a, b in zip(x, y, strict=True)
        ]
        _, csv = write_stitch_tiles(str(d), tiles)
        stitch(d, csv, xcorr_shift=True, max_trial_num=1)
    assert len(frames.calls) == 2
    # (the function string names localWorkerNum only when it is over 1)
    assert "local_worker_num" not in frames.calls[0], "The small dataset used local workers."
    assert frames.calls[1]["local_worker_num"] == 4


def test_an_existing_result_is_kept(tmp_path, frames):
    """a volume whose result exists is done (no frame call)"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    os.makedirs(f"{tmp_path}/matlab_stitch/{NAME}.zarr")
    stitch(tmp_path, csv)
    assert frames.calls == []


def test_a_young_tmp_marker_is_another_process_volume(tmp_path, frames, monkeypatch):
    """(local mode) a volume whose .tmp marker is younger than unit_wait_time minutes per tile
    is left to the process that made it (a wait, no try); an older one is taken over"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    os.makedirs(f"{tmp_path}/matlab_stitch")
    marker = f"{tmp_path}/matlab_stitch/{NAME}.tmp"
    Path(marker).touch()
    waits = []

    def wait(groups, t):
        waits.append((groups, t))
        os.utime(marker, (1, 1))
        return False

    monkeypatch.setattr(jfiles, "wait_for_files", wait)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv)
    assert len(waits) == 1 and waits[0][1] == 5
    assert waits[0][0][0] == [f"{tmp_path}//matlab_stitch/{NAME}.zarr"]
    assert len(frames.calls) == 1 and not os.path.exists(marker)


def test_missing_processed_tiles_count_as_tries(tmp_path, frames, pauses):
    """(local mode) processed tiles that aren't there count as tries: max_trial_num rounds"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    stitch(tmp_path, csv, processed_dir_str="DSR", zarr_file=False, max_trial_num=2)
    assert frames.calls == []
    assert pauses == [5], "The tries were not rounds 5 s apart."


# ---------------------------------------------------------------------------------------------
# StitchingTest in intent: cluster mode against the fake Slurm
# ---------------------------------------------------------------------------------------------


@POSIX
def test_master_does_not_stitch_a_volume_whose_job_runs(tmp_path, fake_slurm, pauses, frames):
    """(masterDoesNotStitchAVolumeWhoseJobRuns) once the master's last volume is done, the
    last unfinished one is one whose job runs: the master waits for it"""
    data, csv, res = two_volumes(tmp_path)
    fake_slurm.set_job_states("1000_1", ["R", "done"], res[:1])
    frames.action = make_result(res[1])
    cluster(data, csv, tmp_path / "logs", master_compute=True, xcorr_mode="all")
    assert len(frames.calls) == 1, "The master stitched a volume whose job was running."
    assert "Scan_Iter_0001" in frames.calls[0]["result_path"]


@POSIX
def test_master_in_a_job_stitches_no_volume_that_does_not_fit_it(
    tmp_path, fake_slurm, pauses, frames, monkeypatch
):
    """(masterInAJobStitchesNoVolumeThatDoesNotFitIt) in a Slurm job without the memory for it,
    the master submits its last volume"""
    from petakit5d.core import memory

    data, csv, res = two_volumes(tmp_path)
    monkeypatch.setenv("SLURM_JOB_ID", "4242")
    monkeypatch.setattr(memory, "get_system_memory", lambda: 1e-6)
    fake_slurm.set_job_states("1000_1", ["R", "done"], res[:1])
    fake_slurm.set_job_states("1001_2", ["R", "done"], res[1:])
    cluster(data, csv, tmp_path / "logs", master_compute=True, xcorr_mode="all")
    assert frames.calls == [], "The master stitched a volume that does not fit its job."
    assert fake_slurm.call_count("sbatch") == 2


@POSIX
@pytest.mark.parametrize(
    "case",
    [
        ((16, 16, 16), "", 1, {}),
        ((16, 16, 16), "on", 2, {}),
        ((512, 512, 512), "", 8, {}),
        ((16, 16, 16), "", 8, {"xcorr_shift": True}),
        ((16, 16, 16), "", 1, {"xcorr_shift": True, "xy_max_offset": 2, "z_max_offset": 1}),
    ],
    ids=[
        "168_batches",
        "config_master_compute",
        "one_batch",
        "xcorr_could_be_small",
        "xcorr_small_offsets",
    ],  # fmt: skip
)
def test_cluster_jobs_ask_for_what_they_hold(tmp_path, fake_slurm, pauses, frames, monkeypatch,
                                             case):  # fmt: skip
    """(clusterJobsAskForWhatTheyHold) a volume of 30 batches or more asks for what it holds
    (one CPU), or with the config's masterCompute a block's job (2 CPUs); a smaller one (or one
    the registration could make so) is stitched in its job (cpus_per_task, 8); the wrapper's
    master_compute reaches the jobs, and they go first (nice)"""
    from petakit5d.jobs import coordinator

    monkeypatch.setattr(coordinator, "coordinator_nice", lambda slurm=None: 1000)
    batch, config, cpus, more = case
    cfg = ""
    if config:
        cfg = str(tmp_path / "cpu_config.json")
        Path(cfg).write_text('{"masterCompute": true}')
    data, csv, res = two_volumes(tmp_path)
    fake_slurm.set_job_states("1000_1", ["R", "done"], res[:1])
    fake_slurm.set_job_states("1001_2", ["R", "done"], res[1:])
    submitted, mems = [], []
    real = single.submit_single_job

    def spy(task, *a, **kw):
        r = real(task, *a, **kw)
        if r[2] > 0:
            submitted.append(task.record())
            mems.append(kw["mem_allocate"])
        return r

    monkeypatch.setattr(single, "submit_single_job", spy)
    cluster(data, csv, tmp_path / "logs", master_compute=False, xcorr_mode="all",
            batch_size=batch, block_size=(16, 16, 16), config_file=cfg, **more)  # fmt: skip
    lines = fake_slurm.calls("sbatch")
    assert len(lines) == 2
    assert cpus_of(fake_slurm) == [cpus, cpus], "the jobs' CPUs"
    # the memory: a coordinator's (the process and its z MIP, uint16), else all the tiles' in
    # single times 2 (MATLAB's 8 bytes per voxel)
    from petakit5d.jobs.coordinator import coordinator_resources
    from petakit5d.stitch.geometry import stitch_volume_size_estimate

    tiles = cut_tiles(random_volume(GRID_GT), GRID_STARTS, GRID_TILE)
    _, _, rng, _ = stitch_volume_size_estimate(
        np.array([t["xyz"] for t in tiles]), GRID_TILE,
        tile_idx=np.array([[*t["idx"], 1] for t in tiles], dtype=np.float64), dsr=True,
        resample_type="isotropic", xcorr_shift=more.get("xcorr_shift", False),
        xy_max_offset=more.get("xy_max_offset", 300), z_max_offset=more.get("z_max_offset", 50),
        batch_size=batch, block_size=(16, 16, 16),
    )  # fmt: skip
    nz, ny, nx = (float(v) for v in rng[1])
    mip = ny * nx * 2 / 2**30
    block = math.prod(min(b, v) for b, v in zip(batch, (nz, ny, nx), strict=True)) * 2 / 2**30 * 20
    if cpus == 8:
        want = 4 * math.prod(GRID_TILE) * 4 / 2**30 * 2
    else:
        want = coordinator_resources(mip, bool(config), 2, block)[1]
    assert mems == [pytest.approx(want, rel=1e-12)] * 2
    assert all("--nice=1000" in x for x in lines), "no nice"
    assert all(t["kwargs"]["master_compute"] is False for t in submitted), (
        "The wrapper's masterCompute did not reach the jobs."
    )


@POSIX
def test_job_sizes_use_the_first_tile_size(tmp_path, fake_slurm, pauses, monkeypatch):
    """(MATLAB's quirk) the jobs' memory comes from the first tile size measured, for every
    volume: a time point of larger tiles gets the first one's"""
    data = tmp_path / "data"
    data.mkdir()
    write_tiles(data, random_volume(GRID_GT), iter=0)
    csv = write_tiles(data, random_volume((40, 100, 120)), tile=(40, 64, 72), iter=1)
    res = [f"{data}/matlab_stitch/{name_of(it)}.zarr" for it in (0, 1)]
    fake_slurm.set_job_states("1000_1", ["R", "done"], res[:1])
    fake_slurm.set_job_states("1001_2", ["R", "done"], res[1:])
    mems = []
    real = single.submit_single_job

    def spy(task, *a, **kw):
        r = real(task, *a, **kw)
        if r[2] > 0:
            mems.append(kw["mem_allocate"])
        return r

    monkeypatch.setattr(single, "submit_single_job", spy)
    cluster(data, csv, tmp_path / "logs", master_compute=False, xcorr_mode="all")
    want = 4 * math.prod(GRID_TILE) * 4 / 2**30 * 2
    assert mems == [pytest.approx(want, rel=1e-12)] * 2


@POSIX
def test_master_waits_for_the_last_try_of_a_job(tmp_path, fake_slurm, pauses, frames):
    """(masterWaitsForTheLastTryOfAJob) with master_compute, once the last try of every
    unfinished volume was submitted, the wrapper waits for those jobs"""
    data, csv, res = two_volumes(tmp_path)
    fake_slurm.set_job_states("1000_1", ["PD", "R", "done"])
    frames.action = make_result(res[1])
    cluster(data, csv, tmp_path / "logs", master_compute=True, xcorr_mode="all", max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 1, "The last try was submitted again."
    # (the status checks; the call's end looks at its task folder's jobs once more)
    checks = [c for c in fake_slurm.history() if c.startswith("squeue -j") and "1000_1" in c]
    assert len(checks) == 3, "The job was not checked until it ended."


def _claim(logs: Path, i: int, job: str) -> None:
    for c in logs.glob(f"{single.POOL_PREFIX}*/claims"):
        try:
            os.symlink(f"{job}:othernode:5:1", c / f"t{i}")
        except FileExistsError:
            pass


@LINUX
@POSIX
def test_last_try_waits_for_a_task_runner_that_has_the_task(tmp_path, fake_slurm, monkeypatch,
                                                            frames):  # fmt: skip
    """(lastTryWaitsForATaskRunnerThatHasTheTask) the second volume's job runs the third
    volume's task (a task runner), so the third's own job ends at once: the wrapper counts it as
    running while the second's job runs, and does not submit it again"""
    data, csv, res = two_volumes(tmp_path, 3)
    logs = tmp_path / "logs"
    logs.mkdir()
    fake_slurm.set_job_states("1000_1", ["R", "done"], res[:1])
    fake_slurm.set_job_states("1001_2", ["R", "R", "R", "R", "done"], res[1:2])
    monkeypatch.setattr(jfiles, "_pause", lambda t: _claim(logs, 2, "1001_2"))
    cluster(data, csv, logs, master_compute=False, xcorr_mode="all", max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 3, "The third volume was submitted again."
    checks = [c for c in fake_slurm.history() if c.startswith("squeue") and "1001_2" in c]
    assert len(checks) >= 5, "The wrapper did not wait for the job that had the third volume."
    assert not list(logs.glob(f"{single.POOL_PREFIX}*")), "The pool is still there."


@POSIX
def test_cluster_mode_sees_results_other_nodes_just_made(tmp_path, fake_slurm, pauses,
                                                         monkeypatch):  # fmt: skip
    """(clusterModeSeesResultsOtherNodesJustMade) results are looked for with batch_file_exist
    (the stub reports them made, so no job is submitted)"""
    data, csv, _ = two_volumes(tmp_path, 1)
    monkeypatch.setattr(jfiles, "batch_file_exist", lambda p, *a, **k: [True] * len(p))
    cluster(data, csv, tmp_path / "logs", master_compute=False, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 0, "A job was submitted for a finished result."


def test_processed_tiles_other_nodes_just_made_are_seen(tmp_path, frames, pauses, monkeypatch):
    """(processedTilesOtherNodesJustMadeAreSeen) processed tiles are looked for with
    batch_file_exist: the stub reports them missing though they are on disk"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt, tile_dir="DSR", save_format="zarr")
    monkeypatch.setattr(jfiles, "batch_file_exist", lambda p, *a, **k: [False] * len(p))
    stitch(tmp_path, csv, zarr_file=True, processed_dir_str="DSR", max_trial_num=1)
    assert frames.calls == [], "The wrapper stitched tiles that batch_file_exist did not see."


@POSIX
def test_cluster_mode_waits_for_the_last_try(tmp_path, fake_slurm, pauses, frames):
    """(clusterModeWaitsForTheLastTry) the last try's job is checked until it ends, one check
    per round (none before submitting), without submitting again"""
    data, csv, _ = two_volumes(tmp_path, 1)
    fake_slurm.set_job_states("1000_1", ["PD", "R", "done"])
    cluster(data, csv, tmp_path / "logs", master_compute=False, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 1, "The last try was submitted again."
    assert sum(c.startswith("-j ") for c in fake_slurm.calls("squeue")) == 3


@POSIX
def test_cluster_mode_ends_soon_after_the_last_result(tmp_path, fake_slurm, pauses):
    """(clusterModeEndsSoonAfterTheLastResult) the results still to come are looked for each
    second; the result comes during the second second"""
    data, csv, res = two_volumes(tmp_path, 1)
    fake_slurm.set_job_states("1000_1", ["R"])
    hook = at_call([0], 2, lambda: os.makedirs(res[0]))
    real = jfiles._pause

    def pause(t):
        real(t)
        hook()

    jfiles._pause = pause
    try:
        cluster(data, csv, tmp_path / "logs", master_compute=False, max_trial_num=1)
    finally:
        jfiles._pause = real
    assert fake_slurm.call_count("sbatch") == 1
    assert not waits_of_a_second_problem(pauses, 2), waits_of_a_second_problem(pauses, 2)


@POSIX
def test_cluster_mode_waits_for_tiles_jobs_still_make(tmp_path, fake_slurm, pauses, frames):
    """(clusterModeWaitsForTilesJobsStillMake) while the pipeline's marker in their folder is
    recent, missing processed tiles are waited for without counting tries, looked for each
    second; they come at the third wait, and the master stitches the volume"""
    data, csv, stage = tiles_to_come(tmp_path)
    (data / "DSR" / SD.DESKEW_MARKER).touch()
    res = f"{data}/matlab_stitch/{NAME}.zarr"
    frames.action = make_result(res)
    real = jfiles._pause

    def move():
        for p in (stage / "DSR").iterdir():
            shutil.move(str(p), str(data / "DSR" / p.name))

    hook = at_call([0], 3, move)

    def pause(t):
        real(t)
        hook()

    jfiles._pause = pause
    try:
        cluster(data, csv, tmp_path / "logs", processed_dir_str="DSR", master_compute=True,
                max_trial_num=1)  # fmt: skip
    finally:
        jfiles._pause = real
    assert len(frames.calls) == 1, "The wrapper gave up tiles that jobs still made."
    assert not waits_of_a_second_problem(pauses, 3), waits_of_a_second_problem(pauses, 3)


@POSIX
def test_cluster_mode_looks_for_the_primary_stitch_info_each_second(tmp_path, fake_slurm,
                                                                    pauses):  # fmt: skip
    """(clusterModeLooksForThePrimaryStitchInfoEachSecond) time point 0: the primary is done,
    its stitch information comes at the second wait, and the secondary's job makes its result;
    time point 1 (primary channel) waits for its tiles until the marker goes (third wait)"""
    data = tmp_path / "data"
    data.mkdir()
    gt = random_volume((32, 64, 104))
    stage = [[0, 0, 0], [40 * PX, 0, 0]]
    tidx = [[0, 0, 0], [1, 0, 0]]
    t = ([(0, 0, 0), (0, 0, 40)], (32, 64, 64))
    write_tiles(data, gt, *t, stage=stage, tidx=tidx, cam="A")
    write_tiles(data, gt, *t, stage=stage, tidx=tidx, cam="B")
    csv = write_tiles(data, gt, *t, stage=stage, tidx=tidx, cam="A", iter=1)
    for cam in "AB":
        write_tiles(data, gt, *t, stage=stage, tidx=tidx, cam=cam, tile_dir="DSR",
                    image_list="ImageList_DSR.csv")  # fmt: skip
    marker = data / "DSR" / SD.DESKEW_MARKER
    marker.touch()
    sp = data / "matlab_stitch"
    (sp / "stitchInfo").mkdir(parents=True)
    Path(sp / f"{NAME}.zarr").touch()
    info = sp / "stitchInfo" / f"{NAME}.mat"
    res = str(sp / f"{name_of(0, 'B')}.zarr")
    # the secondary of time point 0 is the third volume
    fake_slurm.set_job_states("1000_3", ["R", "done"], [res])
    real = jfiles._pause
    hook = at_call([0], [2, 3], [info.touch, marker.unlink])

    def pause(t):
        real(t)
        hook()

    jfiles._pause = pause
    try:
        cluster(data, csv, tmp_path / "logs", channel_patterns=("CamA_ch0", "CamB_ch0"),
                processed_dir_str="DSR", master_compute=False, max_trial_num=1)  # fmt: skip
    finally:
        jfiles._pause = real
    assert fake_slurm.call_count("sbatch") == 1
    # waits of a second: two for the stitch information, five (a round) for the tiles
    assert not waits_of_a_second_problem(pauses, 7), waits_of_a_second_problem(pauses, 7)


@POSIX
@pytest.mark.parametrize("marker_age", [-1, 11], ids=["none", "old"])
def test_cluster_mode_gives_up_tiles_without_a_recent_marker(tmp_path, fake_slurm, pauses,
                                                             frames, marker_age):  # fmt: skip
    """(clusterModeGivesUpTilesWithoutARecentMarker) without a recent marker (none, or one
    left from a pipeline that ended) missing processed tiles count as tries"""
    data, csv, _ = tiles_to_come(tmp_path)
    if marker_age >= 0:
        m = data / "DSR" / SD.DESKEW_MARKER
        m.touch()
        old = time.time() - marker_age * 60
        os.utime(m, (old, old))
    cluster(data, csv, tmp_path / "logs", processed_dir_str="DSR", master_compute=True,
            max_trial_num=2)  # fmt: skip
    assert frames.calls == [], "The wrapper stitched missing tiles."
    assert len(pauses) == 1, "The tiles did not count as tries."


def test_streaming_sees_tiles_other_machines_just_made(tmp_path, pauses, monkeypatch):
    """(streamingSeesTilesOtherMachinesJustMade) streaming checks tiles with batch_file_exist:
    the stub reports every file made (the result too), while there are no raw tiles"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt, tile_dir="DSR")
    monkeypatch.setattr(jfiles, "batch_file_exist", lambda p, *a, **k: [True] * len(p))
    stitch(tmp_path, csv, processed_dir_str="DSR", streaming=True)
    assert pauses == [], "The wrapper waited for tiles batch_file_exist reported made."


# ---------------------------------------------------------------------------------------------
# the function string and its values
# ---------------------------------------------------------------------------------------------


def _values(**over):
    v = {k: "" for k, kind in SD.FRAME_FIELDS if kind == "s"}
    v.update({k: False for k, kind in SD.FRAME_FIELDS if kind == "l"})
    v.update({k: 0 for k, kind in SD.FRAME_FIELDS if kind == "n"})
    v.update(over)
    return v


def test_numbers_go_through_matlabs_conversions():
    s = SD.frame_task_string(
        _values(
            xy_pixel_size=0.1,
            dz=1 / 3,
            skew_angle=32.45,
            tile_offset=1.5,
            xy_max_offset=150,
            z_max_offset=150,
            xcorr_thresh=0.25,
            edge_artifacts=2,
        )  # fmt: skip
    )
    assert "'xyPixelSize',0.1000000000," in s
    assert "'dz',0.3333333333," in s
    assert "'tileOffset',1.500000e+00," in s  # (%d of a non-integer is %e)
    assert "'edgeArtifacts',2.0000000000," in s
    assert s.endswith("'mccMode',false,'configFile','')")


def test_a_vector_where_matlab_formats_a_scalar_raises():
    with pytest.raises(SD.StitchDatasetError, match="xy_max_offset"):
        SD.frame_task_string(_values(xy_max_offset=np.array([150.0, 100.0])))


def test_matrices_and_resample_factors_are_matlabs_text():
    assert SD._mat2str_commas(np.array([[1 / 3, 2.0], [3.0, 4.0]])) == ("[0.333333333333333,2;3,4]")
    assert SD._mat2str_commas(None) == "[]"
    assert SD._mat2str_commas(np.array([False, False, True])) == "[false,false,true]"
    assert SD._mat2str_commas(True) == "true"
    assert SD._num2str_commas(np.array([[2.0, 1.0, 1.0]]), "%.10d,") == (
        "0000000002,0000000001,0000000001,"
    )
    assert SD._num2str_commas(np.array([[1.5, 1.0]]), "%.10d,") == "1.5000000000e+00,0000000001,"
    assert SD._num2str_commas(np.array([[3.0, 5, 2, 90, np.inf, np.inf]]), "%d,") == (
        "3,5,2,90,Inf,Inf,"
    )
    assert SD._num2str_commas(None, "%d,") == ""
    with pytest.raises(SD.StitchDatasetError):
        SD._num2str_commas(np.ones((3, 2)), "%d,")


def test_the_frame_gets_the_values_of_matlabs_string(tmp_path):
    """the function string read back as MATLAB's parser does, bound to stitch_volume: Python
    values (z, y, x) with MATLAB's rounding"""
    v = _values(
        tile_fullpaths="{'/mnt/a/t1.tif','/mnt/a/t2.tif'}", xyz="[0,0,0;6.912,0,0]",
        axis_order="x,y,z", data_order="y,x,z", xy_pixel_size=0.1 / 3, dz=0.3,
        skew_angle=32.45, result_path="/mnt/a/matlab_stitch/r.zarr", stitch_info_dir="stitchInfo",
        pool_size="[]", batch_size="[64,32,16]", block_size="[16,16,8]", shard_size="[]",
        resample_type="isotropic", resample_factor="0000000002,0000000001,0000000001,",
        input_bbox="[]", tile_out_bbox="[]", blend_method="feather", xy_max_offset=300,
        z_max_offset=50, xcorr_downsample="[2,2,1]", xcorr_thresh=0.25, shift_method="grid",
        axis_weight="[1,0.1,10]", out_bbox="3,5,2,90,Inf,Inf,", tile_idx="[0,0,0,1;1,0,0,1]",
        flipped_tile="0,1,", process_fun_path="{''}", stitch_mip="[]", dist_bboxes="[]",
        save_16bit=True, is_primary_ch=True, uuid="ab12", save_mip=True, dsr=True,
    )  # fmt: skip
    task = SD.frame_task(SD.frame_task_string(v))
    kw = task.record()["kwargs"]
    assert task.func == SD.FRAME_FUNC
    assert kw["tile_fullpaths"] == ["/mnt/a/t1.tif", "/mnt/a/t2.tif"]
    assert kw["xy_pixel_size"] == 0.0333333333  # (%0.10f)
    assert kw["batch_size"] == [16, 64, 32] and kw["block_size"] == [8, 16, 16]
    assert kw["resample_factor"] == [1, 2, 1]
    assert kw["xcorr_downsample"] == [1, 2, 2] and kw["axis_weight"] == [10, 1, 0.1]
    assert kw["out_bbox"] == [1, 2, 4, None, 90, None]
    assert kw["data_order"] == "z,y,x"
    assert kw["flipped_tile"] == [0, 1]
    assert kw["tile_idx"] == [[0, 0, 0, 1], [1, 0, 0, 1]]
    assert kw["result_path"] == "/mnt/a/matlab_stitch/r.zarr"


def test_the_wrapper_passes_matlabs_values_to_the_frame(tmp_path, frames):
    """the call's values reach the frame through MATLAB's text (%0.10f, mat2str's 15 digits)"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    stitch(tmp_path, csv, xy_pixel_size=PX + 1e-12, dz=1 / 3, batch_size=(16, 32, 64),
           resample_factor=(1, 1, 1), axis_weight=(10, 1, 1 / 3), out_bbox=(1, 2, 3, 30, None,
           None), uuid="feed")  # fmt: skip
    kw = frames.calls[0]
    assert kw["xy_pixel_size"] == 0.108 and kw["dz"] == 0.3333333333
    assert tuple(kw["batch_size"]) == (16, 32, 64)
    assert tuple(kw["axis_weight"]) == (10, 1, 0.333333333333333)
    assert tuple(kw["out_bbox"]) == (1, 2, 3, 30, None, None)
    assert kw["uuid"] == "feed" and kw["is_primary_ch"] is True
    assert kw["parse_cluster"] is False and kw["stitch_info_dir"] == "stitchInfo"
    assert kw["result_path"] == f"{tmp_path}//matlab_stitch/{NAME}.zarr"


def test_axis_order_loses_its_spaces(tmp_path, frames):
    """MATLAB removes the spaces of axisOrder before its (never failing) check, and the frame
    gets it so"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, axis_order=" x, y ,z")
    assert frames.calls[0]["axis_order"] == "x,y,z"


def test_coordinates_keep_15_digits(tmp_path, frames):
    gt = random_volume(GRID_GT)
    stage = [[t["xyz"][0] + 1e-9 / 3, *t["xyz"][1:]] for t in cut_tiles(gt, GRID_STARTS,
                                                                        GRID_TILE)]  # fmt: skip
    csv = write_tiles(tmp_path, gt, stage=stage)
    stitch(tmp_path, csv)
    got = np.asarray(frames.calls[0]["coordinates"], dtype=np.float64)
    from petakit5d.core.text import through_mat2str

    want = through_mat2str(np.array([[float(f"{v:.15g}") for v in r] for r in stage]))
    # (the image list keeps 15 digits too; the rows in the image list's sorted order)
    assert sorted(map(tuple, got)) == sorted(map(tuple, want))


# ---------------------------------------------------------------------------------------------
# names, modes and files
# ---------------------------------------------------------------------------------------------


def test_names_with_and_without_cameras(tmp_path, frames):
    gt = random_volume(GRID_GT)
    tiles = cut_tiles(gt, GRID_STARTS, GRID_TILE)
    _, csv = write_stitch_tiles(str(tmp_path), tiles, cam="")
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, channel_patterns=("ch0",), primary_ch="ch0")
    assert Path(frames.calls[0]["result_path"]).name == (
        "Scan_Iter_0000_0000_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs.zarr"
    )


def test_garbled_name_when_the_msec_of_tiles_differ(tmp_path, frames):
    """MATLAB passes the unique msec values into the name's sprintf: with two, the format is
    used again from the start: the second msec value is the msecAbs field and the msecAbs value
    goes to the z text's %s, which makes it char(0); the .tmp marker of that name can't be made
    (MATLAB: fclose of fopen's -1). Kept."""
    from petakit5d.core.text import sprintf

    gt = random_volume(GRID_GT)
    _, csv = write_stitch_tiles(str(tmp_path), cut_tiles(gt, GRID_STARTS, GRID_TILE))
    # tile 3's msec field is 5
    old = stitch_tile_name([0, 1, 0])
    new = old.replace("0000000msec", "0000005msec")
    os.rename(f"{tmp_path}/{old}.tif", f"{tmp_path}/{new}.tif")
    Path(csv).write_text(Path(csv).read_text().replace(old, new))
    with pytest.raises(SD.StitchDatasetError) as err:
        stitch(tmp_path, csv)
    assert err.value.identifier == "MATLAB:FileIO:InvalidFid"
    assert frames.calls == []
    name = sprintf(SD._NAME_CAM, f"{tmp_path}//matlab_stitch", "", "0000_0000", "A", 0, 0, 488,
                   np.array([0.0, 5.0]), 0.0, "")  # fmt: skip
    # (the second pass of the format takes the empty z text and ends at its first %s)
    assert name.endswith("_0000000msec_0000000005msecAbs\x00/")


def test_primary_mode_registers_each_time_point(tmp_path, frames):
    """('primary') each time point's primary channel is a primary; the other channels wait for
    it and get its stitch information"""
    gt = random_volume(GRID_GT)
    for it in range(2):
        for cam in "AB":
            csv = write_tiles(tmp_path, gt, iter=it, cam=cam)

    def act(kw):
        os.makedirs(kw["result_path"])
        if kw["is_primary_ch"]:
            n = Path(kw["result_path"]).stem
            os.makedirs(f"{tmp_path}/matlab_stitch/stitchInfo", exist_ok=True)
            Path(f"{tmp_path}/matlab_stitch/stitchInfo/{n}.mat").touch()

    frames.action = act
    stitch(tmp_path, csv, xcorr_mode="primary", channel_patterns=("CamA_ch0", "CamB_ch0"))
    got = [(Path(c["result_path"]).stem, c["is_primary_ch"], Path(c["stitch_info_fullpath"]).stem)
           for c in frames.calls]  # fmt: skip
    # (time points outermost, then cameras)
    assert got == [
        (name_of(0, "A"), True, name_of(0, "A")),
        (name_of(0, "B"), False, name_of(0, "A")),
        (name_of(1, "A"), True, name_of(1, "A")),
        (name_of(1, "B"), False, name_of(1, "A")),
    ]


def test_all_mode_makes_every_volume_a_primary(tmp_path, frames):
    gt = random_volume(GRID_GT)
    for cam in "AB":
        csv = write_tiles(tmp_path, gt, cam=cam)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, xcorr_mode="all", channel_patterns=("CamA_ch0", "CamB_ch0"))
    assert [(c["is_primary_ch"], Path(c["stitch_info_fullpath"]).stem) for c in frames.calls] == [
        (True, name_of(0, "A")),
        (True, name_of(0, "B")),
    ]


def test_a_pinned_stitch_info_file_places_every_volume(tmp_path, frames):
    """(stitch_info_fullpath) every volume is a secondary of the given file ('stitchInfo')"""
    gt = random_volume(GRID_GT)
    for it in range(2):
        csv = write_tiles(tmp_path, gt, iter=it)
    pinned = tmp_path / "pinned.mat"
    pinned.touch()
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, stitch_info_fullpath=str(pinned))
    assert [(c["is_primary_ch"], c["stitch_info_fullpath"]) for c in frames.calls] == [
        (False, str(pinned)),
        (False, str(pinned)),
    ]


def test_a_secondary_without_its_time_points_primary_fails_in_primary_mode(tmp_path, frames):
    """('primary') a time point whose primary channel has no tiles: MATLAB's wrapper marks the
    first time point's volumes as missing (a quirk, kept) and fails on the secondary's
    primary information (index error)"""
    gt = random_volume(GRID_GT)
    write_tiles(tmp_path, gt, iter=0, cam="A")
    write_tiles(tmp_path, gt, iter=0, cam="B")
    csv = write_tiles(tmp_path, gt, iter=1, cam="B")

    def act(kw):
        os.makedirs(kw["result_path"])
        n = Path(kw["result_path"]).stem
        os.makedirs(f"{tmp_path}/matlab_stitch/stitchInfo", exist_ok=True)
        Path(f"{tmp_path}/matlab_stitch/stitchInfo/{n}.mat").touch()

    frames.action = act
    with pytest.raises(SD.StitchDatasetError) as err:
        stitch(tmp_path, csv, xcorr_mode="primary", channel_patterns=("CamA_ch0", "CamB_ch0"))
    assert err.value.identifier == "MATLAB:badsubscript"
    assert [Path(c["result_path"]).stem for c in frames.calls] == [name_of(0, "A"), name_of(0, "B")]


def test_online_stitch_names_by_tile_count_and_removes_fewer(tmp_path, frames):
    """(online_stitch) the result is named by its tile count; a result of fewer tiles goes"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    rt = tmp_path / "matlab_stitch"
    (rt / "MIPs").mkdir(parents=True)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, online_stitch=True, streaming=True)
    assert Path(frames.calls[0]["result_path"]).name == f"{NAME}_000z_004ntile.zarr"


def test_online_stitch_without_streaming_names_by_tile_count(tmp_path, frames):
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    rt = tmp_path / "matlab_stitch"
    old = rt / f"{NAME}_003ntile.zarr"
    old.mkdir(parents=True)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, online_stitch=True)
    assert Path(frames.calls[0]["result_path"]).name == f"{NAME}_004ntile.zarr"
    # (MATLAB's pattern cuts 32 characters of a name that has no z layer: it never matches
    # the results of fewer tiles, which stay)
    assert old.exists()


def test_online_stitch_removes_results_of_fewer_tiles(tmp_path, frames):
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    rt = tmp_path / "matlab_stitch"
    old = rt / f"{NAME}_000z_003ntile.zarr"
    old.mkdir(parents=True)
    (rt / "MIPs").mkdir()
    (rt / "MIPs" / f"{NAME}_000z_003ntile_MIP_z.tif").touch()
    keep = rt / f"{NAME}_000z_005ntile.zarr"
    keep.mkdir()
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, online_stitch=True, streaming=True)
    assert not old.exists() and not (rt / "MIPs" / f"{NAME}_000z_003ntile_MIP_z.tif").exists()
    assert keep.exists()


def test_streaming_waits_for_tiles_and_gives_up(tmp_path, frames, pauses):
    """(streaming) missing tiles are waited for in 5 s rounds; with online stitching the
    wrapper gives up after 5 rounds without a new tile"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    os.remove(f"{tmp_path}/{stitch_tile_name([1, 1, 0])}.tif")
    stitch(tmp_path, csv, streaming=True, online_stitch=True)
    assert frames.calls == []
    # each round counts twice (a missing tile, then the round), the counter passes 5 after 3
    assert pauses == [5, 5, 5]


def test_parameters_files_hold_matlabs_values(tmp_path, frames):
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, batch_size=(16, 32, 64), out_bbox=(1, 2, 3, 30, None, None))
    rt = tmp_path / "matlab_stitch"
    pr = json.loads((rt / "parameters.json").read_text())
    # (jsonencode writes Inf as null)
    assert pr["batchSize"] == [32, 64, 16] and pr["outBbox"] == [3, 4, 2, None, None, 30]
    assert pr["dataOrder"] == "y,x,z" and pr["xcorrDownsample"] == [2, 2, 1]
    assert pr["dataPaths"] == dp(tmp_path) and pr["uuid"] == ""
    mat = read_mat(str(rt / "parameters.mat"))["pr"]
    assert list(np.ravel(mat["batchSize"])) == [32, 64, 16]
    assert (rt / "stitchInfo").is_dir() and (rt / "tmp").is_dir()


@POSIX
def test_more_than_10_tiles_go_to_the_frame_in_a_file(tmp_path, fake_slurm, pauses, submitted):
    """(cluster mode) a volume of more than 10 tiles gets its tiles, coordinates, tile indices
    and flips in stitchInfo/<name>_tile_info.mat"""
    d = tmp_path / "data"
    d.mkdir()
    tiles = [{"im": np.ones((2, 8, 8), np.uint16) * (k + 1), "idx": [k, 0, 0],
              "xyz": [k * 6 * PX, 0.0, 0.0]} for k in range(12)]  # fmt: skip
    _, csv = write_stitch_tiles(str(d), tiles)
    res = f"{d}/matlab_stitch/{NAME}.zarr"
    fake_slurm.set_job_states("1000_1", ["R", "done"], [res])
    cluster(d, csv, tmp_path / "logs", master_compute=False, max_trial_num=1)
    kw = submitted[0]["kwargs"]
    fn = f"{d}//matlab_stitch/stitchInfo/{NAME}_tile_info.mat"
    assert kw["tile_info_fullpath"] == fn and kw["tile_fullpaths"] == []
    assert kw["coordinates"] is None
    a = read_mat(fn)
    paths = [str(p) for p in np.asarray(a["tile_fullpaths"], dtype=object).ravel()]
    assert len(paths) == 12 and paths[0].endswith(".tif")
    assert np.asarray(a["xyz"]).shape == (12, 3) and np.asarray(a["tileIdx"]).shape == (12, 4)
    assert np.size(a["flippedTile"]) == 0
    assert not list(Path(f"{d}/matlab_stitch/stitchInfo").glob("*_tile_info_*.mat"))


@POSIX
def test_the_primary_first_primary_is_stitched_by_the_master(tmp_path, fake_slurm, pauses,
                                                             frames, submitted):  # fmt: skip
    """(cluster mode, master_compute, registration) the registered primary of 'primaryFirst'
    is stitched here first, then the others are submitted but the last, which the master
    takes"""
    data, csv, res = two_volumes(tmp_path, 3)

    def act(kw):
        os.makedirs(kw["result_path"])
        os.makedirs(f"{data}/matlab_stitch/stitchInfo", exist_ok=True)
        Path(f"{data}/matlab_stitch/stitchInfo/{NAME}.mat").touch()

    frames.action = act
    fake_slurm.set_job_states("1000_2", ["R", "done"], res[1:2])
    cluster(data, csv, tmp_path / "logs", master_compute=True, xcorr_shift=True)
    assert [Path(c["result_path"]).stem for c in frames.calls] == [NAME, name_of(2)]
    assert fake_slurm.call_count("sbatch") == 1
    # (master_compute True is the frame's default: the function string doesn't name it)
    assert "master_compute" not in submitted[0]["kwargs"]
    assert "--array=2 " in fake_slurm.calls("sbatch")[0] + " "


@POSIX
def test_multi_loc_volumes_take_their_tile_count_as_their_number(tmp_path, fake_slurm, pauses,
                                                                 frames):  # fmt: skip
    """(multi_loc, MATLAB's quirk) the tile loop reuses the volume's number f: a volume's job
    array index is its number of tiles"""
    gt = random_volume(GRID_GT)
    folders = []
    for k in range(2):
        d = tmp_path / f"loc{k}"
        d.mkdir()
        tiles = cut_tiles(gt, GRID_STARTS, GRID_TILE)[2 * k : 2 * k + 2]
        _, csv = write_stitch_tiles(str(d), tiles, iter=0)
        for it in (1, 2):
            write_stitch_tiles(str(d), tiles, iter=it)
        folders.append((str(d), csv))
    res = [f"{folders[0][0]}/matlab_stitch/{name_of(it)}.zarr" for it in range(3)]
    for j, r in zip((1000, 1001, 1002), res, strict=True):
        fake_slurm.set_job_states(f"{j}_4", ["R", "done"], [r])
    SD.stitch_dataset(
        [dp(f[0]) for f in folders], [f[1] for f in folders], multi_loc=True, dsr=True,
        resample_type="isotropic", xy_pixel_size=PX, xcorr_shift=False,
        channel_patterns=("CamA_ch0",), primary_ch="CamA_ch0", parse_cluster=True,
        master_compute=False, xcorr_mode="all", job_log_dir=str(tmp_path / "logs"),
    )  # fmt: skip
    arrays = [re.search(r"--array=(\d+)", c).group(1) for c in fake_slurm.calls("sbatch")]
    assert arrays == ["4", "4", "4"]


def _settings(path: Path, interval: str) -> None:
    head = "X Stage Offset, Interval (um), # of Pixels for Excitation (0) :"
    # '\n' on every system (an oracle input: Windows' text mode would write '\r\n')
    path.write_text(f"Settings\n{head}\t0\t{interval}\t101\nEnd\n", newline="\n")


def test_setting_files_flip_tiles(tmp_path, frames):
    """(parse_setting_file) tiles whose setting files have a negative stage interval are
    flipped (a tile without one is not: the wrapper has no flip_zstack); the flags reach the
    frame through num2str's '%d,'"""
    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    names = [stitch_tile_name(t["idx"]) for t in cut_tiles(gt, GRID_STARTS, GRID_TILE)]
    for n, iv in zip(names[:3], ("0.5", "-0.5", "-0.3"), strict=True):
        _settings(tmp_path / f"{n}_Settings.txt", iv)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    with pytest.warns(Warning, match="No setting file found"):
        stitch(tmp_path, csv, parse_setting_file=True)
    # (in the image list's order: z, y, x of the tile indices)
    order = sorted(range(4), key=lambda k: cut_tiles(gt, GRID_STARTS, GRID_TILE)[k]["idx"][::-1])
    want = [[0, 1, 1, 0][k] for k in order]
    assert list(frames.calls[0]["flipped_tile"]) == want


def test_processed_tiles_coming(tmp_path):
    d = tmp_path / "DSR"
    d.mkdir()
    tiles = [str(d / "a.zarr")]
    assert not SD.processed_tiles_coming(tiles)
    m = d / SD.DESKEW_MARKER
    m.touch()
    assert SD.processed_tiles_coming(tiles)
    old = time.time() - 11 * 60
    os.utime(m, (old, old))
    assert not SD.processed_tiles_coming(tiles)


def test_process_functions_follow_the_first_tile_s_channel(tmp_path, frames):
    gt = random_volume(GRID_GT)
    for cam in "AB":
        csv = write_tiles(tmp_path, gt, cam=cam)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    stitch(tmp_path, csv, xcorr_mode="all", channel_patterns=("CamA_ch0", "CamB_ch0"),
           process_fun_path=("", ""))  # fmt: skip
    assert [c["process_fun_path"] for c in frames.calls] == [[""], [""]]


def test_paths_follow_matlab_with_windows_rules(monkeypatch):
    """(Windows rules forced) the wrapper's names are text joined with '/' as MATLAB's on any
    platform: a volume's name is the result's name after the last separator of either kind,
    and the pipeline's marker is looked for in each processed tile's folder"""
    from petakit5d.core import paths

    monkeypatch.setattr(paths, "_IS_PC", True)
    rt = "C:\\data\\run\\" + "/matlab_stitch"
    assert paths.fileparts(f"{rt}/{NAME}.zarr")[1] == NAME
    seen = []
    monkeypatch.setattr(SD, "_batch_exist", lambda p: seen.extend(p) or [False] * len(p))
    tiles = ["C:\\data\\run\\DSR\\a.zarr", "C:\\data\\run\\DSR\\b.zarr",
             "D:/other/DSR/c.zarr"]  # fmt: skip
    assert not SD.processed_tiles_coming(tiles)
    assert seen == ["C:\\data\\run\\DSR/deskew_in_progress.tmp",
                    "D:/other/DSR/deskew_in_progress.tmp"]  # fmt: skip


def test_feather_maps_under_the_windows_rule(tmp_path, monkeypatch):
    """(MATLAB's ispc branch forced) the wrapper's feather distance maps go into imdist/ itself,
    not imdist/<first tile>/, and the stitch information names them there; the stitched volume
    and its MIP are those of the other rule (the Linux references, which Windows skips)"""
    from petakit5d.stitch import distance

    gt = random_volume(GRID_GT)
    first = f"{NAME}_000x_000y_000z_0000t"
    got = []
    for pc in (False, True):
        monkeypatch.setattr(distance, "_is_pc", lambda pc=pc: pc)
        d = tmp_path / ("windows" if pc else "posix")
        d.mkdir()
        stitch(d, write_tiles(d, gt), blend_method="feather", batch_size=(16, 32, 32),
               block_size=(16, 16, 16))  # fmt: skip
        folder = d / "matlab_stitch" / "imdist"
        folder = folder if pc else folder / first
        si = V.read_stitch_info(f"{d}/matlab_stitch/stitchInfo/{NAME}.mat", ["imdistFullpaths"])
        assert [Path(p).parent for p in si["imdistFullpaths"]] == [folder]
        assert sorted(os.listdir(folder)) == [f"{first}_wr_10.zarr", "parameters.mat"]
        got.append((result(d), _read_mip(d)))
    for a, b in zip(*got, strict=True):
        np.testing.assert_array_equal(a, b)


def test_the_matlab_name_layer_takes_matlabs_values(tmp_path, frames):
    from petakit5d import matlab

    gt = random_volume(GRID_GT)
    csv = write_tiles(tmp_path, gt)
    frames.action = lambda kw: os.makedirs(kw["result_path"])
    matlab.XR_matlab_stitching_wrapper(
        dp(tmp_path), csv, "DSR", True, "resultDir", "matlab_stitch", "xyPixelSize", PX,
        "batchSize", [64, 32, 16], "xcorrShift", False, "ChannelPatterns", ["CamA_ch0"],
        "parseCluster", False, "outBbox", [1, 1, 1, 90, np.inf, 30],
    )  # fmt: skip
    kw = frames.calls[0]
    assert tuple(kw["batch_size"]) == (16, 64, 32)
    assert tuple(kw["out_bbox"]) == (0, 0, 0, 30, 90, None)


# ---------------------------------------------------------------------------------------------
# oracle cases: the wrapper in local mode, MATLAB's against the port's
# ---------------------------------------------------------------------------------------------

# the registration results a scene gives are written by each side (MATLAB's save, the port's
# write_xcorr_result) before its runs, as the frame's cases do (a .mat file holds its time)

# StitchingTest.stitch's defaults (exact DSR stitching), with small batches and blocks
_BASE = {
    "DSR": True,
    "resampleType": "isotropic",
    "xyPixelSize": PX,
    "axisOrder": "x,y,z",
    "xcorrShift": False,
    "channelPatterns": ["CamA_ch0"],
    "primaryCh": "CamA_ch0",
    "parseCluster": False,
    "batchSize": [32, 32, 16],
    "blockSize": [16, 16, 16],
}
_TWO_CAMS = {"channelPatterns": ["CamA_ch0", "CamB_ch0"]}
_SHIFT = [[0, 0, 0], [0, 0, 0], [-2, 1, 0], [0, 0, 0]]  # (x, y, z) voxels per tile


@dataclass
class Group:
    """the tiles of one volume, cut from a ground truth at the scene's starts"""

    it: int = 0
    cam: str = "A"
    folder: str = ""
    seed: int = 1
    keep: tuple | None = None  # the tiles kept (indices), else all
    tile_dir: str = ""
    fmt: str = "tif"
    intervals: tuple | None = None  # stage intervals of the tiles' setting files


@dataclass
class WScene:
    groups: list
    runs: list  # per run: MATLAB name-value pairs (with _BASE), its folders and image lists
    results: list  # result paths relative to the root
    gt: tuple = GRID_GT
    starts: list = field(default_factory=lambda: list(GRID_STARTS))
    tile: tuple = GRID_TILE
    pixel: tuple = (PX, PX, PX)
    stage: list | None = None
    tidx: list | None = None
    # registration results given: (result folder, group index, shifts in voxels (x, y, z))
    shifts: list = field(default_factory=list)


def _run(args=None, data=("",), lists=("ImageList.csv",), path_args=()):
    return {
        "args": {**_BASE, **(args or {})},
        "data": list(data),
        "lists": list(lists),
        "path_args": list(path_args),
    }


def _names(groups):
    return [f"matlab_stitch/{name_of(g.it, g.cam)}.zarr" for g in groups]


_TS = [Group(it, cam, seed=1 + (cam == "B")) for it in (0, 1) for cam in "AB"]
_SKEW_STAGE = [[0, 0, 0], [0, 36 * PX, 0], [10 * 0.3, 0, 0], [10 * 0.3, 36 * PX, 0]]
_NOCAM = "Scan_Iter_0000_0000_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs"

WSCENES: dict[str, WScene] = {
    # a time series of two cameras registered once: every volume on the primary's canvas
    "timeseries.feather.xcorr": WScene(
        _TS,
        [_run({**_TWO_CAMS, "xcorrShift": True, "BlendMethod": "feather"})],
        _names(_TS),
        shifts=[("matlab_stitch", 0, _SHIFT)],
    ),
    # 'primary': each time point's primary channel registered (its result given)
    "primary_mode.mean.xcorr": WScene(
        _TS,
        [_run({**_TWO_CAMS, "xcorrShift": True, "xcorrMode": "primary", "BlendMethod": "mean"})],
        _names(_TS),
        shifts=[
            ("matlab_stitch", 0, _SHIFT),
            ("matlab_stitch", 2, [[0, 0, 0], [1, 0, 0], [0, 0, 0], [0, -1, 0]]),
        ],
    ),
    "all_mode.none": WScene(
        _TS[:2], [_run({**_TWO_CAMS, "xcorrMode": "all", "BlendMethod": "none"})], _names(_TS[:2])
    ),
    # a pinned stitch information file (of a first run) for every volume of a second run
    "pinned_stitch_info.max": WScene(
        [Group(0, "A"), Group(0, "B", seed=2), Group(1, "A", seed=3)],
        [
            _run({"xcorrShift": True, "BlendMethod": "max"}),
            _run(
                {
                    **_TWO_CAMS,
                    "xcorrShift": True,
                    "BlendMethod": "max",
                    "resultDirName": "matlab_stitch_pinned",
                },
                path_args=("stitchInfoFullpath", f"matlab_stitch/stitchInfo/{NAME}.mat"),
            ),
        ],
        [f"matlab_stitch/{NAME}.zarr"]
        + [f"matlab_stitch_pinned/{name_of(i, c)}.zarr" for i, c in ((0, "A"), (0, "B"), (1, "A"))],
        shifts=[("matlab_stitch", 0, _SHIFT)],
    ),
    # processed DSR zarr tiles (MATLAB's layout) of the raw names
    "dsr_zarr.feather": WScene(
        [Group(tile_dir="DSR", fmt="zarr")],
        [_run({"zarrFile": True, "processedDirStr": "DSR", "BlendMethod": "feather"})],
        [f"matlab_stitch/{NAME}.zarr"],
    ),
    "skewed.none": WScene(
        [Group()],
        [
            _run(
                {
                    "DSR": False,
                    "dz": 0.3,
                    "skewAngle": 32.45,
                    "reverse": True,
                    "axisOrder": "-x,y,z",
                    "BlendMethod": "none",
                }
            )
        ],
        [f"matlab_stitch/{NAME}.zarr"],
        gt=(50, 100, 64),
        starts=[(0, 0, 0), (0, 36, 0), (10, 0, 0), (10, 36, 0)],
        tile=(40, 64, 64),
        stage=_SKEW_STAGE,
        tidx=[[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]],
    ),
    # tiles flipped by their setting files (negative stage intervals), in skewed space
    "settings.skewed_flipped.max": WScene(
        [Group(intervals=("0.3", "-0.3", "0.3", "-0.3"))],
        [
            _run(
                {
                    "DSR": False,
                    "dz": 0.3,
                    "skewAngle": 32.45,
                    "reverse": True,
                    "axisOrder": "-x,y,z",
                    "BlendMethod": "max",
                    "parseSettingFile": True,
                }
            )
        ],
        [f"matlab_stitch/{NAME}.zarr"],
        gt=(50, 100, 64),
        starts=[(0, 0, 0), (0, 36, 0), (10, 0, 0), (10, 36, 0)],
        tile=(40, 64, 64),
        stage=_SKEW_STAGE,
        tidx=[[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]],
    ),
    "no_cam.mean": WScene(
        [Group(cam="")],
        [_run({"channelPatterns": ["ch0"], "primaryCh": "ch0", "BlendMethod": "mean"})],
        [f"matlab_stitch/{_NOCAM}.zarr"],
    ),
    # tiles in two folders, an image list each
    "multi_loc.feather": WScene(
        [Group(folder="loc0/", keep=(0, 1)), Group(folder="loc1/", keep=(2, 3))],
        [
            _run(
                {"multiLoc": True, "BlendMethod": "feather"},
                data=("loc0/", "loc1/"),
                lists=("loc0/ImageList.csv", "loc1/ImageList.csv"),
            )
        ],
        [f"loc0/matlab_stitch/{NAME}.zarr"],
    ),
    "online.none": WScene(
        [Group()],
        [_run({"onlineStitch": True, "BlendMethod": "none"})],
        [f"matlab_stitch/{NAME}_004ntile.zarr"],
    ),
    "streaming_online.mean": WScene(
        [Group()],
        [_run({"streaming": True, "onlineStitch": True, "BlendMethod": "mean"})],
        [f"matlab_stitch/{NAME}_000z_004ntile.zarr"],
    ),
    "local_workers.feather": WScene(
        [Group()],
        [_run({"BlendMethod": "feather", "localWorkerNum": 2})],
        [f"matlab_stitch/{NAME}.zarr"],
    ),
    # 'primaryFirst': a time point of another tile count is skipped
    "primary_first.other_count": WScene(
        [Group(0), Group(1, keep=(0, 1, 2))],
        [_run({"BlendMethod": "none"})],
        [f"matlab_stitch/{NAME}.zarr"],
    ),
    # a resample factor through the function string's %.10d
    "resampled.feather": WScene(
        [Group()],
        [_run({"BlendMethod": "feather", "resampleType": "given", "resampleFactor": [2, 1, 1]})],
        [f"matlab_stitch/{NAME}.zarr"],
        gt=(32, 80, 108),
        starts=[(0, 0, 0), (0, 0, 36), (0, 30, 0), (0, 30, 36)],
        tile=(32, 50, 72),
        pixel=(PX, 2 * PX, PX),
    ),
}


def _scene_tiles(sc: WScene, g: Group) -> list[dict]:
    gt = random_volume(sc.gt, g.seed, "s09")
    tiles = cut_tiles(gt, sc.starts, sc.tile, sc.pixel, sc.stage, sc.tidx)
    return [t for k, t in enumerate(tiles) if g.keep is None or k in g.keep]


def _scene_layout(sc: WScene, g: Group, args: dict):
    """the frame's coordinates in stitching space (S03's), for a registration result"""
    from petakit5d.core import axes as ax
    from petakit5d.stitch import geometry as G

    tiles = _scene_tiles(sc, g)
    rf = args.get("resampleFactor")
    return G.prepare_layout(
        np.array([t["xyz"] for t in tiles]),
        np.array([t["im"].shape for t in tiles]),
        tile_idx=np.array([[*t["idx"], 1] for t in tiles], dtype=np.float64),
        axis_order=args.get("axisOrder", "x,y,z"),
        data_order="z,y,x",
        xy_pixel_size=args.get("xyPixelSize", PX),
        dz=args.get("dz", 0.5),
        skew_angle=args.get("skewAngle", 32.45),
        io_scan=args.get("IOScan", False),
        dsr=args.get("DSR", False),
        resample_type=args.get("resampleType", "xy_isotropic"),
        resample_factor=None if rf is None else ax.from_yxz(rf),
    )


def make_wroot(name: str):
    """the tiles and image lists of a scene"""

    def make(d):
        sc = WSCENES[name]
        for g in sc.groups:
            folder = os.path.join(os.fspath(d), g.folder)
            os.makedirs(folder, exist_ok=True)
            write_stitch_tiles(
                folder,
                _scene_tiles(sc, g),
                cam=g.cam,
                iter=g.it,
                tile_dir=g.tile_dir,
                save_format=g.fmt,
                zarr_layout="matlab",
            )
            for tl, iv in zip(_scene_tiles(sc, g), g.intervals or (), strict=False):
                n = stitch_tile_name(tl["idx"], g.cam, g.it)
                _settings(Path(folder, f"{n}_Settings.txt"), iv)

    return make


def wparams(name: str) -> dict:
    """the case's struct for MATLAB: per run its folders, image lists and name-value pairs
    (paths relative to the root), and the results to read"""
    sc = WSCENES[name]
    runs = []
    for r in sc.runs:
        flat = []
        for k, v in r["args"].items():
            flat += [k, v]
        runs.append(
            {"dataPaths": r["data"], "lists": r["lists"], "args": flat, "pathArgs": r["path_args"]}
        )
    xcorr = []
    for folder, gi, shift in sc.shifts:
        g = sc.groups[gi]
        lay = _scene_layout(sc, g, sc.runs[0]["args"])
        dx = np.asarray(shift, dtype=np.float64)[lay.sinds] * lay.space.pixel_sizes
        xcorr.append({"folder": folder, "name": name_of(g.it, g.cam), "xyz_shift": lay.xyz + dx,
                      "dxyz_shift": dx})  # fmt: skip
    return {"runs": runs, "results": list(sc.results), "xcorr": xcorr}


WRAPPER_CODE = r"""
p = in{1}; d = in{2};
for k = 1 : numel(p.xcorr)
	x = p.xcorr{k};
	mkdir([d x.folder '/xcorr']);
	xyz_shift = x.xyz_shift; dxyz_shift = x.dxyz_shift;
	save('-v7.3', [d x.folder '/xcorr/' x.name '.mat'], 'xyz_shift', 'dxyz_shift');
end
for k = 1 : numel(p.runs)
	q = p.runs{k};
	args = q.args;
	for j = 1 : 2 : numel(q.pathArgs)
		args = [args, {q.pathArgs{j}, [d q.pathArgs{j + 1}]}];
	end
	dps = cellfun(@(x) [d x], q.dataPaths, 'unif', 0);
	lists = cellfun(@(x) [d x], q.lists, 'unif', 0);
	if numel(dps) == 1
		dps = dps{1};
		lists = lists{1};
	end
	XR_matlab_stitching_wrapper(dps, lists, args{:});
end
out = cell(1, numel(p.results));
for k = 1 : numel(p.results)
	out{k} = readzarr([d p.results{k}]);
end
"""


def run_wrapper(name: str, d: str) -> list:
    """the Python side: the same runs through the MATLAB-name layer"""
    from petakit5d import matlab

    p = wparams(name)
    for x in p["xcorr"]:
        os.makedirs(f"{d}{x['folder']}/xcorr", exist_ok=True)
        V.write_xcorr_result(f"{d}{x['folder']}/xcorr/{x['name']}.mat", x["xyz_shift"],
                             x["dxyz_shift"])  # fmt: skip
    for q in p["runs"]:
        args = list(q["args"])
        pa = q["pathArgs"]
        for j in range(0, len(pa), 2):
            args += [pa[j], d + pa[j + 1]]
        dps = [d + x for x in q["dataPaths"]]
        lists = [d + x for x in q["lists"]]
        if len(dps) == 1:
            dps, lists = dps[0], lists[0]
        matlab.XR_matlab_stitching_wrapper(dps, lists, *args)
    return [read_zarr(d + r) for r in p["results"]]


def _wrapper_case(name: str) -> Case:
    n = len(WSCENES[name].results)
    return Case(
        f"wrapper.{name}",
        code=WRAPPER_CODE,
        args=[In("p"), InPath("root")],
        nargout=n,
        python=run_wrapper,
        py_args=[name, InPath("root")],
        out_roles=("vol",) * n,
        inputs={
            "p": Value(lambda name=name: wparams(name)),
            "root": Tree(make_wroot(name), "root"),
        },
        tree=tree(ignore_vars=("t", "t1", "ip")),
        marks=("slow",),
    )


CASES = [_wrapper_case(n) for n in WSCENES]

# MATLAB on Windows (ispc) writes the feather distance maps into imdist/ itself, not
# imdist/<first tile>/, and so does the port there: the references that hold such maps, made
# on Linux, hold the other folder and paths. The Windows layout is tested on every system by
# forcing that branch (test_feather_maps_under_the_windows_rule).
_WINDOWS_MAPS = (
    "the feather distance maps go into imdist/ on Windows (MATLAB's ispc branch); this "
    "reference is MATLAB's on Linux (imdist/<first tile>/): "
    "test_feather_maps_under_the_windows_rule covers it"
)


def _maps_in_tile_folder(case: Case) -> bool:
    """whether the case's MATLAB reference holds files in a folder of imdist/"""
    files = oracle.reference(case).files
    return files is not None and any(re.search(r"(^|/)imdist/[^/]+/", n) for n in files.names())


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    from petakit5d.stitch import distance

    if distance._is_pc() and _maps_in_tile_folder(case):
        pytest.skip(_WINDOWS_MAPS)
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# cluster mode on a real Slurm cluster: the same trees as local mode
# ---------------------------------------------------------------------------------------------


def _tree_files(root: Path) -> dict[str, Path]:
    """the results' files and zarr arrays (left out: job logs, flags, the parameters files,
    which name the mode, and the block information, whose task size is the mode's: all batches
    in one task here, 20 of 512^3 per cluster job)"""
    out = {}
    for d, dirs, files in os.walk(root):
        rel = Path(d).relative_to(root).as_posix()
        if rel.startswith(("job_logs", "logs")) or "zarr_flags" in rel or "_task_size_" in rel:
            dirs[:] = []
            continue
        files = [f for f in files if "_task_size_" not in f]
        for x in [x for x in dirs if x.endswith(".zarr")]:
            out[f"{rel}/{x}"] = Path(d, x)
        dirs[:] = [x for x in dirs if not x.endswith(".zarr")]
        for f in files:
            if not f.startswith("parameters."):
                out[f"{rel}/{f}"] = Path(d, f)
    return out


@pytest.mark.cluster
@pytest.mark.slow
def test_cluster_mode_gives_the_local_trees():
    """a time series of two cameras ('primaryFirst', registration results given) stitched in
    local mode and in cluster mode (one job per volume, the master computing none): the same
    volumes, MIPs, converted tiles, maps and stitch information"""
    import uuid

    from petakit5d.io import read_tiff

    scratch = os.environ.get("PETAKIT5D_TEST_SLURM_SCRATCH")
    if not scratch or not os.path.isdir(scratch):
        pytest.skip("PETAKIT5D_TEST_SLURM_SCRATCH names no folder")
    root = Path(scratch) / f"pk5d_stitch_test_{uuid.uuid4().hex[:12]}"
    name = "timeseries.feather.xcorr"
    try:
        trees = []
        for mode in ("local", "cluster"):
            d = root / mode
            d.mkdir(parents=True)
            make_wroot(name)(d)
            p = wparams(name)
            for x in p["xcorr"]:
                os.makedirs(f"{d}/{x['folder']}/xcorr", exist_ok=True)
                V.write_xcorr_result(
                    f"{d}/{x['folder']}/xcorr/{x['name']}.mat", x["xyz_shift"], x["dxyz_shift"]
                )
            logs = root / f"logs_{mode}"
            logs.mkdir()
            stitch(
                d, d / "ImageList.csv", channel_patterns=("CamA_ch0", "CamB_ch0"),
                xcorr_shift=True, blend_method="feather", batch_size=(16, 32, 32),
                block_size=(16, 16, 16), parse_cluster=mode == "cluster", master_compute=False,
                config_file=os.environ["PETAKIT5D_TEST_SLURM_CONFIG"], job_log_dir=str(logs),
            )  # fmt: skip
            trees.append(d)
        # the volumes were stitched by Slurm jobs (the master computes none): their logs say so
        logs = [f.read_text() for f in (root / "logs_cluster").rglob("job_*.out")]
        made = sum(t.count("Stitching result") + t.count("Stitch for tiles") for t in logs)
        assert made >= 4, logs
        a, b = _tree_files(trees[0]), _tree_files(trees[1])
        assert sorted(a) == sorted(b)
        for rel in a:
            if rel.endswith(".zarr"):
                np.testing.assert_array_equal(read_zarr(str(a[rel])), read_zarr(str(b[rel])))
            elif rel.endswith(".tif"):
                np.testing.assert_array_equal(read_tiff(str(a[rel])), read_tiff(str(b[rel])))
            elif rel.endswith(".json"):
                ta = a[rel].read_text().replace(str(trees[0]), "<R>")
                assert ta == b[rel].read_text().replace(str(trees[1]), "<R>"), rel
    finally:
        shutil.rmtree(root, ignore_errors=True)
