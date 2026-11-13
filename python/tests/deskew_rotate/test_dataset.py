"""Deskew/rotate of datasets (petakit5d.deskew_rotate.dataset, task D04).

Oracle cases against MATLAB: ``XR_deskew_rotate_data_wrapper`` in local mode on copies of
dataset folders (MATLAB runs the oracle tree, Python the port through the MATLAB-name layer;
the folders they leave are compared, tiffs byte for byte, zarr by voxels and metadata after the
layout mapping, parameters.mat by content, parameters.json as text): the combined path on a
dataset of two channels and two time points, channel patterns, two data paths, deskew only
(also with the combined default), the separate path, overwrite, zarr input and output, chunk
sizes, single output, input boxes, partial files, flat field per channel, no MIPs, resampling,
objective scans, and the large-file path (batches of rows that don't divide the result, input
boxes, no MIPs, single output, resampling, MIP masks, its errors).

Then DeskewRotateDataWrapperTest and DeskewRotateMemoryTest's wrapper methods in intent on the
port's (z, y, x) data: outputs equal the in-memory functions, the job shapes (a framework spy),
the publish order, the job of a large file going after its blocks (a fake Slurm).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.deskew_rotate import dataset as dd
from petakit5d.deskew_rotate import memory as dm
from petakit5d.io import create_zarr, read_tiff, read_zarr, write_tiff, write_zarr
from petakit5d.testing import oracle, slurm_fake, synth
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG

SKEW = 32.45
PX = 0.108
DZ = 0.3
RAW = (29, 37, 61)  # (z, y, x): not cubic; 37 rows give batches of 16, 16 and 5
FSN_A0 = "scan_Iter_0000_CamA_ch0_stack0000_0000000msec_0000t"
FSN_B0 = "scan_Iter_0000_CamB_ch0_stack0000_0000000msec_0000t"
FSN_A1 = "scan_Iter_0001_CamA_ch0_stack0000_0011000msec_0001t"
FSNS = (FSN_A0, FSN_B0, FSN_A1)
ERRORS = (dd.DeskewRotateDatasetError, RuntimeError, ValueError, OSError, IndexError)


def make_frame(shape=RAW, kind="uint16", seed=0):
    """beads on a background of 100 with noise, (z, y, x)"""
    nz, ny, nx = shape
    rng = CounterRNG(seed, f"dsrds/{kind}/{'x'.join(map(str, shape))}")
    pos = [
        (2 + (k * 7 + seed) % max(1, nz - 4), 3 + (k * 5) % max(1, ny - 6), 4 + (k * 11) % (nx - 8))
        for k in range(9)
    ]
    vol = synth.beads(shape, pos, 1.5, 1000.0, 100.0)
    noise = rng.integers(0, 50, shape, dtype=np.uint16)
    if kind == "single":
        return np.ascontiguousarray(vol + noise.astype(np.float32))
    v = np.floor(vol).astype(np.uint16) + noise
    if kind == "uint8":
        return np.ascontiguousarray((v // 8).astype(np.uint8))
    return np.ascontiguousarray(v)


# ------------------------------------------------------------------------------------------
# oracle cases
# ------------------------------------------------------------------------------------------


def dataset_tree(
    *, zarr=False, names=FSNS, parts=1, ff=False, kind="uint16", sub="", masks=False, chunks=None
):
    """a dataset folder: the files ``names`` (seed = their index), as tiff (``parts`` partial
    files each) or as zarr in MATLAB's layout; flat-field images ff_A.tif, ff_B.tif, bg.tif;
    MIP masks mask_{z,y,x}.tif of a box inside the volume"""

    def make(folder):
        d = Path(folder) / sub
        d.mkdir(parents=True, exist_ok=True)
        for k, name in enumerate(names):
            vol = make_frame(RAW, kind, k)
            if zarr:
                synth.write_zarr(d / f"{name}.zarr", vol, chunks or (16, 16, 32), layout="matlab")
                continue
            edges = np.linspace(0, RAW[0], parts + 1).astype(int)
            for p in range(parts):
                fn = name if p == 0 else f"{name}_part{p:04d}"
                synth.write_tiff(d / f"{fn}.tif", vol[edges[p] : edges[p + 1]])
        if ff:
            ny, nx = RAW[1], RAW[2]
            bg = np.full((ny, nx), 100, dtype=np.float32)
            ramp = np.linspace(500, 1000, nx, dtype=np.float64).astype(np.float32)
            synth.write_tiff(d / "ff_A.tif", np.ascontiguousarray(bg + ramp[np.newaxis]))
            synth.write_tiff(d / "ff_B.tif", np.ascontiguousarray(bg + ramp[np.newaxis, ::-1]))
            synth.write_tiff(d / "bg.tif", bg)
        if masks:
            m = np.zeros(RAW, dtype=np.uint8)
            m[8:22, 10:30, 15:45] = 1
            mdir = Path(folder) / "masks"
            mdir.mkdir(exist_ok=True)
            synth.write_tiff(mdir / "mask_z.tif", m.max(axis=0))
            synth.write_tiff(mdir / "mask_y.tif", np.ascontiguousarray(m.max(axis=1).T))
            synth.write_tiff(mdir / "mask_x.tif", np.ascontiguousarray(m.max(axis=2).T))

    return make


def two_paths_tree(folder):
    for k, t in enumerate(("tile1", "tile2")):
        d = Path(folder) / t
        d.mkdir(parents=True)
        synth.write_tiff(d / f"{FSN_A0}.tif", make_frame(RAW, "uint16", k))


def existing_tree(folder):
    """a dataset whose first file already has a DSR result (7s)"""
    dataset_tree()(folder)
    (Path(folder) / "DSR").mkdir()
    synth.write_tiff(Path(folder) / "DSR" / f"{FSN_A0}.tif", np.full((3, 4, 5), 7, np.uint16))


# MATLAB's cd to the repository root changes nothing in the work folder
_CODE = (
    "p = in{1}; if ischar(p), p = {p}; end; p = cellfun(@(x) [in{2}, x], p, 'unif', 0); "
    "a = in{3}; for k = 1 : numel(a), if ischar(a{k}) && startsWith(a{k}, '@d/'), "
    "a{k} = [in{2}, a{k}(4 : end)]; elseif iscell(a{k}), for j = 1 : numel(a{k}), "
    "if startsWith(a{k}{j}, '@d/'), a{k}{j} = [in{2}, a{k}{j}(4 : end)]; end, end, end, end; "
    "XR_deskew_rotate_data_wrapper(p, a{:}); out = {};"
)


def _expand(v, d):
    if isinstance(v, str) and v.startswith("@d/"):
        return d + v[3:]
    if isinstance(v, list) and all(isinstance(x, str) for x in v):
        return [_expand(x, d) for x in v]
    return v


def p_wrapper(sub, d, args):
    from petakit5d import matlab

    paths = [d + s for s in (sub if isinstance(sub, list) else [sub])]
    matlab.XR_deskew_rotate_data_wrapper(paths, *[_expand(a, d) for a in args])
    return []


BASE = [
    "xyPixelSize", PX, "dz", DZ, "skewAngle", SKEW, "reverse", True,
    "channelPatterns", ["CamA", "CamB"], "parseCluster", False, "masterCompute", True,
]  # fmt: skip
LARGE = ["zarrFile", True, "largeFile", True, "batchSize", [16, 1024, 1024],
         "blockSize", [16, 64, 64]]  # fmt: skip
# the writer records of a large file that failed (its owner is the process)
_TREE = tree(masks=[(r"owner=[^\n]*", "owner=<owner>")])


def ds_case(cid, make, *args, sub="d/", error=False):
    return Case(
        f"ds.{cid}",
        code=_CODE,
        nargout=0,
        args=[sub, InPath("d"), [*BASE, *args]],
        python=p_wrapper,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
        tree=_TREE,
    )


def _cases():
    t = dataset_tree(sub="d")
    z = dataset_tree(zarr=True, sub="d")
    z1 = dataset_tree(zarr=True, names=(FSN_A0,), sub="d")
    ffa = ["FFCorrection", True, "FFImagePaths", ["@d/d/ff_A.tif", "@d/d/ff_B.tif"],
           "backgroundPaths", ["@d/d/bg.tif", "@d/d/bg.tif"]]  # fmt: skip
    masks = ["maskFullpaths", ["@d/masks/mask_z.tif", "@d/masks/mask_y.tif", "@d/masks/mask_x.tif"]]
    return [
        ds_case("combined", t),
        ds_case("channel_b", t, "channelPatterns", ["CamB"]),
        ds_case("two_paths", two_paths_tree, sub=["tile1", "tile2"]),
        ds_case("deskew_only", t, "rotate", False, "DSRCombined", False),
        ds_case("deskew_only_combined_default", t, "rotate", False, "DSRCombined", True),
        ds_case("separate", t, "DSRCombined", False),
        ds_case("existing", existing_tree, sub=""),
        ds_case("existing_overwrite", existing_tree, "overwrite", True, sub=""),
        ds_case("save_zarr", t, "saveZarr", True, "blockSize", [16, 32, 16]),
        ds_case("zarr_input", z, "zarrFile", True),
        ds_case("single_out", t, "save16bit", False),
        ds_case("input_bbox", t, "inputBbox", [3, 5, 2, 30, 58, 27]),
        ds_case("partial_files", dataset_tree(parts=3, names=(FSN_A0,), sub="d")),
        ds_case("flatfield", dataset_tree(ff=True, sub="d"), *ffa),
        ds_case("no_mip", t, "saveMIP", False),
        ds_case("resample", t, "resampleFactor", [2, 2, 2]),
        ds_case("objective_scan", t, "objectiveScan", True, "dz", 0.2),
        ds_case("single_in", dataset_tree(kind="single", names=(FSN_A0,), sub="d")),
        # large files
        ds_case("large", z, *LARGE),
        ds_case(
            "large.uneven", z1, *LARGE[:4], "batchSize", [13, 1024, 1024], "blockSize", [8, 32, 32]
        ),
        ds_case("large.one_task", z1, *LARGE, "taskSize", 1),
        ds_case("large.no_mip", z1, *LARGE, "saveMIP", False),
        ds_case("large.single_out", z1, *LARGE, "save16bit", False),
        ds_case("large.input_bbox", z1, *LARGE, "inputBbox", [3, 5, 2, 30, 58, 27]),
        ds_case("large.forward", z1, *LARGE, "reverse", False),
        ds_case("large.resample", z1, *LARGE, "resampleFactor", [1.5, 1.5, 1.5]),
        ds_case(
            "large.uint8_in",
            dataset_tree(zarr=True, kind="uint8", names=(FSN_A0,), sub="d"),
            *LARGE,
        ),  # fmt: skip
        ds_case(
            "large.single_in",
            dataset_tree(zarr=True, kind="single", names=(FSN_A0,), sub="d"),
            *LARGE,
        ),  # fmt: skip
        ds_case(
            "large.masks",
            dataset_tree(zarr=True, names=(FSN_A0,), sub="d", masks=True),
            *LARGE,
            *masks,
        ),  # fmt: skip
        ds_case("large.objective", z1, *LARGE, "objectiveScan", True, error=True),
        ds_case("large.deskew_only", z1, *LARGE, "rotate", False, error=True),
        ds_case("large.resample_one", z1, *LARGE, "resampleFactor", 1.5, error=True),
    ]


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# DeskewRotateDataWrapperTest in intent
# ------------------------------------------------------------------------------------------

LOCAL = {
    "xy_pixel_size": PX,
    "dz": DZ,
    "skew_angle": SKEW,
    "reverse": True,
    "channel_patterns": ("CamA", "CamB"),
    "parse_cluster": False,
    "master_compute": True,
}
LARGE_KW = {"zarr_file": True, "large_file": True, "batch_size": (1024, 16, 1024),
            "block_size": (64, 16, 64)}  # fmt: skip


def run(paths, **kw):
    dd.deskew_rotate_dataset(paths, **{**LOCAL, **kw})


def expected_dsr(frame, dz=DZ):
    from petakit5d.deskew_rotate import deskew_rotate_volume

    return deskew_rotate_volume(frame, SKEW, dz, PX, reverse=True, save_16bit=True)


def write_dataset(folder, zarr=False, names=FSNS):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    frames = []
    for k, n in enumerate(names):
        f = make_frame(RAW, "uint16", k)
        frames.append(f)
        if zarr:
            write_zarr(f, str(folder / f"{n}.zarr"), block_size=(16, 16, 32))
        else:
            write_tiff(f, str(folder / f"{n}.tif"))
    return f"{folder}/", frames


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, b.dtype, a.shape, b.shape)
    np.testing.assert_array_equal(a, b)


def test_processes_every_matching_file(tmp_path):
    path, frames = write_dataset(tmp_path / "d")
    run(path)
    for name, f in zip(FSNS, frames, strict=True):
        dsr = read_tiff(f"{path}DSR/{name}.tif")
        same(dsr, expected_dsr(f))
        same(read_tiff(f"{path}DSR/MIPs/{name}_MIP_z.tif"), dsr.max(axis=0))
    assert os.path.isfile(f"{path}DSR/parameters.json")
    assert os.path.isfile(f"{path}DSR/parameters.mat")
    assert not list(Path(path, "DSR").glob("*.tmp")), "task flag files were left behind"
    pr = json.loads(Path(f"{path}DSR/parameters.json").read_text())
    assert pr["dz"] == DZ and pr["skewAngle"] == SKEW


def test_channel_patterns_select_files(tmp_path):
    path, _ = write_dataset(tmp_path / "d")
    run(path, channel_patterns=("CamB",))
    assert sorted(os.listdir(f"{path}DSR")) == sorted(
        ["MIPs", "parameters.json", "parameters.mat", f"{FSN_B0}.tif"]
    )


def test_multiple_data_paths_are_processed(tmp_path):
    frame = make_frame()
    paths = []
    for t in ("tile1", "tile2"):
        (tmp_path / t).mkdir()
        write_tiff(frame, str(tmp_path / t / f"{FSN_A0}.tif"))
        paths.append(str(tmp_path / t))
    run(paths)
    for p in paths:
        same(read_tiff(f"{p}/DSR/{FSN_A0}.tif"), expected_dsr(frame))


def _deskewed(frame):
    from petakit5d.deskew_rotate import deskew_volume
    from petakit5d.deskew_rotate.file import _to_class

    ds = deskew_volume(frame.astype(np.float32), SKEW, DZ, PX, True, crop=False)
    return _to_class(ds, "uint16", 1)


@pytest.mark.parametrize("combined", [False, True])
def test_deskew_only_writes_ds(tmp_path, combined):
    # DSRCombined=true (the default) used to send deskew only to the combined path
    path, frames = write_dataset(tmp_path / "d")
    run(path, rotate=False, dsr_combined=combined)
    same(read_tiff(f"{path}DS/{FSN_A0}.tif"), _deskewed(frames[0]))
    assert not os.path.isdir(f"{path}DSR")
    assert os.path.isfile(f"{path}DS/parameters.json")


def test_separate_deskew_rotate_writes_both(tmp_path):
    path, frames = write_dataset(tmp_path / "d")
    run(path, dsr_combined=False)
    assert os.path.isfile(f"{path}DS/{FSN_A0}.tif")
    assert read_tiff(f"{path}DSR/{FSN_A0}.tif").shape == expected_dsr(frames[0]).shape


def test_overwrite_removes_the_ds_folder_of_the_separate_path(tmp_path):
    path, _ = write_dataset(tmp_path / "d", names=(FSN_A0,))
    os.makedirs(f"{path}DS")
    Path(f"{path}DS/old.txt").write_text("x")
    run(path, dsr_combined=False, overwrite=True)
    assert not os.path.exists(f"{path}DS/old.txt"), "overwrite kept the DS folder"
    assert os.path.isfile(f"{path}DS/{FSN_A0}.tif")


def test_existing_outputs_are_skipped_unless_overwrite(tmp_path):
    path, frames = write_dataset(tmp_path / "d")
    run(path)
    fn = f"{path}DSR/{FSN_A0}.tif"
    write_tiff(frames[0] * 2, f"{path}{FSN_A0}.tif")
    run(path)
    same(read_tiff(fn), expected_dsr(frames[0]))
    run(path, overwrite=True)
    same(read_tiff(fn), expected_dsr(frames[0] * 2))


def test_save_zarr_matches_tiff(tmp_path):
    path, frames = write_dataset(tmp_path / "d")
    run(path, save_zarr=True)
    same(read_zarr(f"{path}DSR/{FSN_A0}.zarr"), expected_dsr(frames[0]))


def test_zarr_input_matches_tiff_input(tmp_path):
    path, frames = write_dataset(tmp_path / "d", zarr=True, names=(FSN_A0,))
    run(path, zarr_file=True)
    same(read_tiff(f"{path}DSR/{FSN_A0}.tif"), expected_dsr(frames[0]))


RAW16 = (29, 64, 61)  # 64 rows: batches of 16 (DeskewRotateDataWrapperTest's)


def test_large_file_matches_regular(tmp_path):
    # the large file path splits the output along y into batches (4 here) and must
    # reproduce the in-memory result exactly (dz 0.2 as in MATLAB's test: no skewed space
    # interpolation, whose uint16 tail rule depends on a batch's input rows)
    path = f"{tmp_path}/"
    frame = make_frame(RAW16)
    write_zarr(frame, f"{path}{FSN_A0}.zarr", block_size=(16, 16, 32))
    run(path, **LARGE_KW, dz=0.2)
    dsr = read_zarr(f"{path}DSR/{FSN_A0}.zarr")
    same(dsr, expected_dsr(frame, 0.2))
    same(read_tiff(f"{path}DSR/MIPs/{FSN_A0}_MIP_z.tif"), dsr.max(axis=0))
    # (MATLAB's orientation: the y MIP is nx x nz, the x MIP ny x nz)
    same(read_tiff(f"{path}DSR/MIPs/{FSN_A0}_MIP_y.tif"), dsr.max(axis=1).T)
    same(read_tiff(f"{path}DSR/MIPs/{FSN_A0}_MIP_x.tif"), dsr.max(axis=2).T)


@pytest.mark.parametrize("rows", [5, 13, 32, 37, 70])
def test_large_file_batches_change_only_the_uint16_tail(tmp_path, rows):
    # dz 0.3 interpolates in skewed space, whose uint16 sums take the first floor(n / 16) * 16
    # values of a y line pairwise and the rest left to right, n the rows of a batch's input
    # (MATLAB's mex): the rows near a batch's end differ by 1 count in a few voxels from the
    # in-memory result, as in MATLAB (the oracle cases compare them with MATLAB's)
    path = f"{tmp_path}/"
    frame = make_frame(RAW16)
    write_zarr(frame, f"{path}{FSN_A0}.zarr", block_size=(16, 16, 32))
    run(path, **{**LARGE_KW, "batch_size": (1024, rows, 1024), "block_size": (8, 8, 16)})
    got = read_zarr(f"{path}DSR/{FSN_A0}.zarr").astype(np.int32)
    want = expected_dsr(frame).astype(np.int32)
    diff = np.abs(got - want)
    assert diff.max() <= 1 and np.count_nonzero(diff) < diff.size * 1e-3
    if rows >= RAW16[1]:
        assert diff.max() == 0


def test_large_file_passes_job_settings_to_mips(tmp_path, monkeypatch):
    from petakit5d.tools import mip

    calls = []
    real = mip.mip_zarr
    monkeypatch.setattr(mip, "mip_zarr", lambda *a, **k: (calls.append(k), real(*a, **k)))
    path, _ = write_dataset(tmp_path / "d", zarr=True, names=(FSN_A0,))
    run(path, **LARGE_KW, job_log_dir=f"{path}logs")
    assert calls, "mip_zarr was not called"
    assert calls[0]["parse_cluster"] is False
    assert calls[0]["master_compute"] is True
    assert calls[0]["job_log_dir"] == f"{path}logs"


def test_large_file_makes_mips_before_publishing(tmp_path):
    from petakit5d.testing.publish_spy import PublishSpy

    path, _ = write_dataset(tmp_path / "d", zarr=True, names=(FSN_A0,))
    with PublishSpy() as spy:
        run(path, **LARGE_KW)
    fn = f"{path}DSR/{FSN_A0}.zarr"
    assert spy.was_moved(fn), "the result was not moved into place"
    when = spy.files_when_moved(fn)
    for lab in ("y", "x", "z"):
        assert any(w.endswith(f"MIPs/{FSN_A0}_MIP_{lab}.tif") for w in when), lab


def test_large_file_mips_are_merged_from_the_blocks(tmp_path, monkeypatch):
    from petakit5d.testing.publish_spy import PublishSpy
    from petakit5d.tools import mip

    def no_blocks(*a, **k):
        raise AssertionError("The MIPs were computed from the result.")

    monkeypatch.setattr(mip, "mip_block", no_blocks)
    path, frames = write_dataset(tmp_path / "d", zarr=True, names=(FSN_A0,))
    with PublishSpy() as spy:
        run(path, **LARGE_KW)
    flags = [c for c in spy.calls if "/zarr_flag/" in c.dst and c.dst.endswith(".mat")]
    assert flags, "the block flags were not moved into place"
    dsr = expected_dsr(frames[0])
    same(read_tiff(f"{path}DSR/MIPs/{FSN_A0}_MIP_z.tif"), dsr.max(axis=0))


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls (MATLAB's FrameworkSpyFixture); nothing runs"""
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        return [fake.done] * len(list(outputs))

    fake.done = True
    monkeypatch.setattr(framework, "run_tasks", fake)
    fake.calls = calls
    return fake


fake_slurm = slurm_fake.fake_slurm  # (the fixture)


@pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
def test_large_file_jobs_go_after_their_blocks(tmp_path, fake_slurm, framework_spy, monkeypatch):
    from petakit5d.jobs import coordinator

    create_zarr(str(tmp_path / f"{FSN_A0}.zarr"), data_size=(64, 512, 256), block_size=(64, 64, 64))
    monkeypatch.setattr(coordinator, "coordinator_nice", lambda slurm=None: 1000)
    (tmp_path / "logs").mkdir()
    for large in (True, False):
        run(
            str(tmp_path),
            zarr_file=True,
            large_file=large,
            parse_cluster=True,
            dsr_dir_name=f"DSR_{int(large)}",
            job_log_dir=str(tmp_path / "logs"),
        )
    assert framework_spy.calls[0]["nice"] == 1000
    assert framework_spy.calls[1]["nice"] == 0


def test_large_file_jobs_ask_for_what_they_hold(tmp_path, framework_spy):
    # the job of a large file asks for the process and the merge of the MIP slabs on one CPU,
    # and if it computes batches itself (master_compute, the config's if it sets it), for a
    # batch's estimate on the CPUs of a batch's job. The input is an 8 G-voxel header.
    create_zarr(
        str(tmp_path / f"{FSN_A0}.zarr"), data_size=(512, 8192, 2048), block_size=(256, 256, 256)
    )
    off = tmp_path / "cpu_config.json"
    off.write_text('{"masterCompute": false}', newline="\n")
    for i, (mc, cfg) in enumerate([(True, ""), (False, ""), (True, str(off))]):
        run(
            str(tmp_path),
            zarr_file=True,
            large_file=True,
            master_compute=mc,
            config_file=cfg,
            cpus_per_task=8,
            dsr_dir_name=f"DSR_{i}",
        )
    c = framework_spy.calls
    # a batch of 1024 rows across x and z, less than half of the input in single (32 GB)
    assert c[0]["cpus_per_task"] == 8
    assert 2 < c[0]["mem_allocate"] < 16
    for i in (1, 2):
        assert c[i]["cpus_per_task"] == 1, i
        assert c[i]["mem_allocate"] < 4, i
    # (the large-file tasks get MATLAB's task-string values)
    t = c[0]["tasks"][0]
    assert t.kwargs["cpus_per_task"] == 8 and list(t.kwargs["mask_fullpaths"]) == [""]


def test_flat_field_uses_each_channels_images(tmp_path):
    from petakit5d.flatfield import flatfield_correct_volume

    path, frames = write_dataset(tmp_path / "d")
    ny, nx = RAW[1:]
    bg = np.full((ny, nx), 100, dtype=np.float32)
    ff_a = bg + np.linspace(500, 1000, nx, dtype=np.float32)[np.newaxis]
    ff_b = bg + np.linspace(1000, 500, nx, dtype=np.float32)[np.newaxis]
    os.makedirs(f"{path}ff")
    for n, im in (("ffA", ff_a), ("ffB", ff_b), ("bg", bg)):
        write_tiff(np.ascontiguousarray(im), f"{path}ff/{n}.tif")
    run(
        path,
        ff_correction=True,
        lower_limit=0.4,
        ff_image_paths=(f"{path}ff/ffA.tif", f"{path}ff/ffB.tif"),
        background_paths=(f"{path}ff/bg.tif", f"{path}ff/bg.tif"),
    )
    for k, ff in ((0, ff_a), (1, ff_b)):
        corrected = flatfield_correct_volume(frames[k], ff, bg, lower_limit=0.4)
        same(read_tiff(f"{path}DSR/{FSNS[k]}.tif"), expected_dsr(corrected))


def test_flat_field_images_are_selected_as_matlab_indexes_them(tmp_path):
    # FFImagePaths{mask}: the first true place; none and a place past the end raise
    path, _ = write_dataset(tmp_path / "d", names=(FSN_B0,))
    with pytest.raises(dd.DeskewRotateDatasetError) as err:
        run(path, channel_patterns=("CamA", "CamB"), ff_correction=True, ff_image_paths=("",))
    assert err.value.identifier == "MATLAB:badsubscript"
    with pytest.raises(dd.DeskewRotateDatasetError) as err:
        dd._first_selected(("a",), [])
    assert err.value.identifier == "MATLAB:needMoreRhsOutputs"
    with pytest.raises(dd.DeskewRotateDatasetError) as err:
        dd._first_selected(("a",), [0, 1])
    assert err.value.identifier == "MATLAB:badsubscript"
    assert dd._first_selected(("a", "b", "c"), [1, 2]) == "b"


def test_partial_files_are_joined_along_z(tmp_path):
    frame = make_frame()
    (tmp_path / "d").mkdir()
    write_tiff(frame[:12], str(tmp_path / "d" / f"{FSN_A0}.tif"))
    write_tiff(frame[12:], str(tmp_path / "d" / f"{FSN_A0}_part0001.tif"))
    run(str(tmp_path / "d"))
    same(read_tiff(str(tmp_path / "d" / "DSR" / f"{FSN_A0}.tif")), expected_dsr(frame))
    assert not (tmp_path / "d" / "DSR" / f"{FSN_A0}_part0001.tif").exists()


def test_save_mip_false_writes_no_mip(tmp_path):
    path, _ = write_dataset(tmp_path / "d")
    run(path, save_mip=False)
    assert not os.path.isdir(f"{path}DSR/MIPs")


def test_block_size_is_used_for_zarr_output(tmp_path):
    from petakit5d.io import get_zarr_info

    path, _ = write_dataset(tmp_path / "d")
    run(path, save_zarr=True, block_size=(16, 16, 32))
    assert tuple(get_zarr_info(f"{path}DSR/{FSN_A0}.zarr").block_size) == (16, 16, 32)


def test_objective_scan_large_file_is_rejected(tmp_path):
    path, _ = write_dataset(tmp_path / "d", zarr=True, names=(FSN_A0,))
    with pytest.raises(dd.DeskewRotateDatasetError, match="objective scan is not supported"):
        run(path, **LARGE_KW, objective_scan=True)
    assert not os.path.isdir(f"{path}DSR"), "folders were made before the error"


def test_both_disabled_does_nothing(tmp_path, capsys):
    path, _ = write_dataset(tmp_path / "d")
    run(path, deskew=False, rotate=False)
    assert "Both deskew and rotation are disabled" in capsys.readouterr().out
    assert sorted(os.listdir(path)) == sorted(f"{n}.tif" for n in FSNS)


def test_objective_scan_keeps_the_default_folder_name(tmp_path, framework_spy):
    # MATLAB renames the folder to rotated/ only for 'DSR/', never for the default 'DSR'
    path, _ = write_dataset(tmp_path / "d", names=(FSN_A0,))
    run(path, objective_scan=True)
    run(path, objective_scan=True, dsr_dir_name="DSR/")
    assert framework_spy.calls[0]["outputs"][0].endswith(f"/DSR/{FSN_A0}.tif")
    assert os.path.normpath(framework_spy.calls[1]["outputs"][0]).endswith(
        os.path.normpath(f"/rotated/{FSN_A0}.tif")
    )


def test_task_values_take_matlabs_text_rounding(tmp_path, framework_spy):
    path, _ = write_dataset(tmp_path / "d", names=(FSN_A0,))
    run(path, const_offset=1 / 3, skew_angle=32.4512345678)
    t = framework_spy.calls[0]["tasks"][0]
    assert t.kwargs["const_offset"] == 0.3333333333
    assert t.kwargs["skew_angle"] == 32.4512345678
    path2, _ = write_dataset(tmp_path / "z", zarr=True, names=(FSN_A0,))
    run(path2, **LARGE_KW, skew_angle=32.4512345678)
    tz = framework_spy.calls[1]["tasks"][0]
    # (%d of a non-integer: %e, 7 significant digits)
    assert tz.kwargs["skew_angle"] == 32.45123


def test_tries_double_the_resources(tmp_path, framework_spy):
    path, _ = write_dataset(tmp_path / "d", names=(FSN_A0,))
    framework_spy.done = False
    with pytest.raises(dd.DeskewRotateDatasetError, match="1 / 1 failed"):
        run(path, cpus_per_task=2)
    c = framework_spy.calls
    assert [x["cpus_per_task"] for x in c] == [2, 4, 8]
    m = c[0]["mem_allocate"]
    assert [x["mem_allocate"] for x in c] == pytest.approx([m, 2 * m, 4 * m])
    assert all(x["max_trial_num"] == 2 for x in c)


# ------------------------------------------------------------------------------------------
# DeskewRotateMemoryTest's wrapper methods
# ------------------------------------------------------------------------------------------


def estimate(size, **kw):
    return dm.estimate_deskew_rotate_memory(
        size, PX, DZ, skew_angle=SKEW, block_size=(256,) * 3, **kw
    )[0]


def test_wrapper_requests_the_estimate(tmp_path, framework_spy):
    (tmp_path / "d").mkdir()
    write_tiff(np.zeros(RAW, np.uint16), str(tmp_path / "d" / f"{FSN_A0}.tif"))
    run(str(tmp_path / "d"))
    assert framework_spy.calls[0]["mem_allocate"] == pytest.approx(estimate(RAW), rel=1e-12)
    # options that change the code path change the estimate
    run(str(tmp_path / "d"), save_16bit=False, save_zarr=True)
    assert framework_spy.calls[1]["mem_allocate"] == pytest.approx(
        estimate(RAW, save_16bit=False, save_zarr=True), rel=1e-12
    )


def test_wrapper_counts_all_partial_files(tmp_path, framework_spy):
    frame = np.zeros(RAW, np.uint16)
    (tmp_path / "d").mkdir()
    write_tiff(frame[:12], str(tmp_path / "d" / f"{FSN_A0}.tif"))
    write_tiff(frame[12:], str(tmp_path / "d" / f"{FSN_A0}_part0001.tif"))
    run(str(tmp_path / "d"))
    assert framework_spy.calls[0]["mem_allocate"] == pytest.approx(
        estimate(RAW, partial_file=True), rel=1e-12
    )
    # (normalized: the separators of the given folder and of MATLAB's '/' join meet on Windows)
    assert [os.path.normpath(p) for p in framework_spy.calls[0]["tasks"][0].args[0]] == [
        os.path.normpath(f"{tmp_path}/d/{FSN_A0}.tif"),
        os.path.normpath(f"{tmp_path}/d/{FSN_A0}_part0001.tif"),
    ]


def test_wrapper_uses_the_largest_file(tmp_path, framework_spy):
    (tmp_path / "d").mkdir()
    write_zarr(np.zeros(RAW, np.uint16), str(tmp_path / "d" / f"{FSN_A0}.zarr"))
    large = (RAW[0] * 2, RAW[1], RAW[2])
    write_zarr(np.zeros(large, np.uint16), str(tmp_path / "d" / f"{FSN_B0}.zarr"))
    run(str(tmp_path / "d"), zarr_file=True)
    assert framework_spy.calls[0]["mem_allocate"] == pytest.approx(estimate(large), rel=1e-12)


def test_wrapper_counts_the_cropped_input(tmp_path, framework_spy):
    (tmp_path / "d").mkdir()
    write_tiff(np.zeros(RAW, np.uint16), str(tmp_path / "d" / f"{FSN_A0}.tif"))
    box = (0, 0, 0, RAW[0], 16, RAW[2])
    run(str(tmp_path / "d"), input_bbox=box)
    assert framework_spy.calls[0]["mem_allocate"] == pytest.approx(
        estimate((RAW[0], 16, RAW[2]), read_size=RAW), rel=1e-12
    )


def test_wrapper_counts_preprocessing_in_the_estimate(tmp_path, framework_spy):
    # flips, flat field and input axis orders hold another copy of the input
    (tmp_path / "d").mkdir()
    write_tiff(np.zeros(RAW, np.uint16), str(tmp_path / "d" / f"{FSN_A0}.tif"))
    run(str(tmp_path / "d"), input_axis_order="xyz")
    run(str(tmp_path / "d"), flip_zstack=True)
    for c in framework_spy.calls:
        assert c["mem_allocate"] == pytest.approx(estimate(RAW, preprocess=True), rel=1e-12)


def test_each_file_gets_its_own_flip(tmp_path, framework_spy, monkeypatch):
    # (the setting files' flips, per file, not the wrapper's flipZstack)
    from petakit5d.pipeline import filenames

    real = filenames.parse_dataset_filenames

    def flips(*a, **k):
        names = real(*a, **k)
        names.flip_zstack = np.array([k % 2 == 1 for k in range(len(names.fnames))])
        return names

    monkeypatch.setattr(filenames, "parse_dataset_filenames", flips)
    path, _ = write_dataset(tmp_path / "d")
    run(path, flip_zstack=False)
    got = [t.kwargs["flip_zstack"] for t in framework_spy.calls[0]["tasks"]]
    assert got == [False, True, False]


def test_controller_resources_hold_one_batch_of_rows():
    # the job of a large file: a batch of batch_size's y rows (and a border row on each side
    # of its input) across all of x and z, with the merge of the MIP slabs
    from petakit5d.jobs.coordinator import coordinator_resources
    from petakit5d.tools.mip import merge_memory

    size = (512, 8192, 2048)
    out = dm.estimate_deskew_rotate_memory(size, PX, DZ, skew_angle=SKEW)[2]
    for computes, mip in ((True, True), (False, True), (True, False)):
        got = dd.controller_resources(
            size, PX, DZ, SKEW, True, "linear", None, None, (1024, 512, 1024), (256,) * 3,
            mip, "uint16", computes, 8,
        )  # fmt: skip
        block = (min(256, out[0]), 256, 256)
        task = dm.estimate_deskew_rotate_memory(
            (512, 514, 2048), PX, DZ, skew_angle=SKEW, out_size=(out[0], 512, out[2]),
            large_file=True, block_size=block,
        )[0]  # fmt: skip
        merge = merge_memory(out, (out[0], 512, out[2]), block, "uint16") if mip else 0.0
        assert got == pytest.approx(coordinator_resources(merge, computes, 8, task), rel=1e-12)


# ------------------------------------------------------------------------------------------
# cluster mode on a real Slurm cluster: the same trees as local mode
# ------------------------------------------------------------------------------------------


def _tree_files(root):
    out = {}
    for dirpath, _, files in os.walk(root):
        for f in files:
            p = Path(dirpath, f)
            rel = str(p.relative_to(root))
            if rel.endswith((".mat", ".json")) or "/zarr_flag/" in f"/{rel}":
                continue
            out[rel] = p
    return out


@pytest.mark.cluster
@pytest.mark.slow
@pytest.mark.parametrize("large", [False, True])
def test_cluster_mode_gives_the_local_trees(large):
    import shutil
    import uuid

    scratch = os.environ.get("PETAKIT5D_TEST_SLURM_SCRATCH")
    if not scratch or not os.path.isdir(scratch):
        pytest.skip("PETAKIT5D_TEST_SLURM_SCRATCH names no folder")
    root = Path(scratch) / f"pk5d_dsr_test_{uuid.uuid4().hex[:12]}"
    try:
        kw = dict(LARGE_KW) if large else {}
        trees = []
        for mode in ("local", "cluster"):
            path, _ = write_dataset(root / mode, zarr=large)
            (root / f"logs_{mode}").mkdir()
            run(
                path,
                **kw,
                parse_cluster=mode == "cluster",
                master_compute=False,
                config_file=os.environ["PETAKIT5D_TEST_SLURM_CONFIG"],
                job_log_dir=str(root / f"logs_{mode}"),
            )
            trees.append(Path(path, "DSR"))
        # the files were made by Slurm jobs (the master computes none): their logs say so
        logs = [f.read_text() for f in (root / "logs_cluster").rglob("job_*.out")]
        made = sum(t.count(" done in ") for t in logs)
        assert made >= len(FSNS), logs
        if large:
            # (the block jobs, submitted by the job of the file)
            assert any("deskew_rotate_block" in t for t in logs), logs
        a, b = _tree_files(trees[0]), _tree_files(trees[1])
        assert sorted(a) == sorted(b)
        for rel in a:
            if rel.endswith(".tif"):
                same(read_tiff(str(a[rel])), read_tiff(str(b[rel])))
            elif ".zarr/" in rel and rel.endswith(".zarray"):
                z = str(a[rel].parent)
                same(read_zarr(z), read_zarr(str(b[rel].parent)))
    finally:
        shutil.rmtree(root, ignore_errors=True)
