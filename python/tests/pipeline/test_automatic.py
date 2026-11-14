"""The microscope pipeline (petakit5d.pipeline.automatic, task P02) against MATLAB's
XR_microscopeAutomaticProcessing (the oracle tree's), and AutomaticProcessingTest in intent.

Oracle cases ``pipeline.<scene>``: each side gets the same dataset folders (tiles cut from a raw
volume along y, their image list, flat-field and background images), then runs the pipeline
in local mode (MATLAB's function, the port's through its MATLAB-name layer, with the same
arguments). Compared: every file the runs wrote (the DS and DSR results and their MIPs, the
stitched volumes, their MIPs and stitch information, the stitching's parameters files, the
pixel-size markers), zarr by voxels and metadata after the axis mapping. Scenes: the parameter
sets of the demos (demo_microscope_pipeline_single_dataset: two cameras, two time points, flat
field, stitching; ..._nostitch; ..._multiple_dataset: two datasets with an image list each),
the separate deskew and rotation, deskew only, zarr results, partial files, an input box, a
blend method, resampling (a factor with rx != ry, which the stitching reads as [y x z]),
overwrite, and MATLAB's errors and warnings.

The pytest functions port AutomaticProcessingTest's methods (local mode, and cluster mode
against the fake Slurm: the waits, the last tries, task runners, the stitching job after the
deskew jobs, the deskew marker, results other nodes just made, streaming with files that come
during the run) and test the pipeline's own logic (its function strings and their values,
MATLAB's quirks). A real-Slurm test (``cluster``) checks that cluster mode makes the local
mode's trees.
"""

from __future__ import annotations

import os
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pytest

from petakit5d.io import read_tiff, read_zarr
from petakit5d.jobs import files as jfiles
from petakit5d.jobs import single
from petakit5d.pipeline import automatic as A
from petakit5d.stitch import distance as sd
from petakit5d.testing import oracle, slurm_fake, synth
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG, cut_tiles, stitch_tile_name, write_stitch_tiles

fake_slurm = slurm_fake.fake_slurm  # (the fixture)

POSIX = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
LINUX = pytest.mark.skipif(
    not os.sys.platform.startswith("linux"), reason="Task runners run on Linux only."
)

PX = 0.108
DZ = 0.3
SKEW = 32.45
RAW = (24, 56, 40)  # (z, y, x) of the raw volume: not cubic
TILE_Y = 36  # tiles along y overlapping by 16 rows: rows 0-35 and 20-55
STARTS = [(0, 0, 0), (0, 20, 0)]
STITCH = "Scan_Iter_{it:04d}_0000_Cam{cam}_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs"
ERRORS = (A.AutomaticProcessingError, ValueError, IndexError, RuntimeError, OSError)


def raw_volume(seed: int = 0, shape=RAW) -> np.ndarray:
    """beads on a background of 200 with noise, (z, y, x) uint16"""
    nz, ny, nx = shape
    pos = [(4 + (k * 5 + seed) % (nz - 8), 3 + (k * 13) % (ny - 6), 4 + (k * 7) % (nx - 8))
           for k in range(8)]  # fmt: skip
    vol = synth.beads(shape, pos, 1.5, 2000.0, 200.0)
    noise = CounterRNG(seed, f"p02/{nz}x{ny}x{nx}").integers(0, 40, shape, dtype=np.uint16)
    return np.ascontiguousarray(np.floor(vol).astype(np.uint16) + noise)


def tiles_of(vol: np.ndarray) -> list[dict]:
    return cut_tiles(vol, STARTS, (RAW[0], TILE_Y, RAW[2]), (PX * 2, PX, PX))


def ff_images(folder: Path, cams) -> None:
    """flat fields falling off to about half at the right edge, and a background of 100"""
    ny, nx = TILE_Y, RAW[2]
    ramp = 1100 - (np.arange(nx, dtype=np.int64) * 500) // (nx - 1)
    for c, cam in enumerate(cams):
        ff = np.broadcast_to((100 + ramp * (c + 1)).astype(np.uint16), (ny, nx))
        synth.write_tiff(folder / f"ff_{cam}.tif", np.ascontiguousarray(ff))
        synth.write_tiff(folder / f"bg_{cam}.tif", np.full((ny, nx), 100, np.uint16))


def dataset(cams=("A", "B"), iters=(0,), ff=False, parts=1, seed=0, sub="d"):
    """a dataset folder (``sub`` of the tree): per time point and camera the two tiles of a raw
    volume (partial files with parts > 1) and their image list ImageList.csv; flat-field images
    FF/ff_<cam>.tif and FF/bg_<cam>.tif"""

    def make(folder):
        d = Path(folder) / sub
        d.mkdir(parents=True, exist_ok=True)
        for it in iters:
            for k, cam in enumerate(cams):
                tiles = tiles_of(raw_volume(seed + 10 * it + k))
                if parts == 1:
                    write_stitch_tiles(str(d), tiles, cam=cam, iter=it)
                    continue
                for t in tiles:
                    name = stitch_tile_name(t["idx"], cam, it)
                    edges = np.linspace(0, RAW[0], parts + 1).astype(int)
                    for p in range(parts):
                        fn = name if p == 0 else f"{name}_part{p:04d}"
                        synth.write_tiff(d / f"{fn}.tif", t["im"][edges[p] : edges[p + 1]])
        if ff:
            (d / "FF").mkdir()
            ff_images(d / "FF", cams)

    return make


def two_datasets(folder):
    for k, sub in enumerate(("d1", "d2")):
        dataset(cams=("A",), seed=7 * k + 1, sub=sub)(folder)


def existing(folder):
    """a dataset whose first file has a DSR result (7s) already"""
    dataset(cams=("A",), sub="")(folder)
    (Path(folder) / "DSR").mkdir()
    name = stitch_tile_name((0, 0, 0), "A", 0)
    synth.write_tiff(Path(folder) / "DSR" / f"{name}.tif", np.full((3, 4, 5), 7, np.uint16))


# MATLAB's cd to the repository root changes nothing in the work folder
_CODE = (
    "p = in{1}; if ischar(p), p = {p}; end; p = cellfun(@(x) [in{2}, x], p, 'unif', 0); "
    "a = in{3}; for k = 1 : numel(a), if ischar(a{k}) && startsWith(a{k}, '@d/'), "
    "a{k} = [in{2}, a{k}(4 : end)]; elseif iscell(a{k}), for j = 1 : numel(a{k}), "
    "if startsWith(a{k}{j}, '@d/'), a{k}{j} = [in{2}, a{k}{j}(4 : end)]; end, end, end, end; "
    "XR_microscopeAutomaticProcessing(p, a{:}); out = {};"
)


def _expand(v, d):
    if isinstance(v, str) and v.startswith("@d/"):
        return d + v[3:]
    if isinstance(v, list) and all(isinstance(x, str) for x in v):
        return [_expand(x, d) for x in v]
    return v


def p_pipeline(sub, d, args):
    from petakit5d import matlab

    paths = [d + s for s in (sub if isinstance(sub, list) else [sub])]
    matlab.XR_microscopeAutomaticProcessing(paths, *[_expand(a, d) for a in args])
    return []


BASE = [
    "xyPixelSize", PX, "dz", DZ, "reverse", True, "streaming", False, "parseCluster", False,
    "cpusPerTask", 4,
]  # fmt: skip
FF = ["FFCorrection", True, "lowerLimit", 0.4, "constOffset", 1,
      "FFImagePaths", ["@d/d/FF/ff_B.tif", "@d/d/FF/ff_A.tif"],
      "backgroundPaths", ["@d/d/FF/bg_B.tif", "@d/d/FF/bg_A.tif"]]  # fmt: skip
STITCH_ON = ["stitch", True, "imageListFullpaths", "@d/d/ImageList.csv", "axisOrder", "-x,y,z"]
# (the stitching's .mat files hold MATLAB's inputParser object ip, which the port's don't: S08's
# listed difference)
_TREE = tree(ignore_vars=("t", "t1", "ip"))


def pl_case(cid, make, *args, sub="d/", error=False, marks=()):
    return Case(
        f"pipeline.{cid}",
        code=_CODE,
        nargout=0,
        args=[sub, InPath("d"), [*BASE, *args]],
        python=p_pipeline,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
        tree=_TREE,
        marks=marks,
    )


def _cases():
    demo = dataset(cams=("B", "A"), iters=(0, 1), ff=True)
    one = dataset(cams=("A",))
    chs = ["channelPatterns", ["CamB_ch0", "CamA_ch0"]]
    cha = ["channelPatterns", ["CamA_ch0"]]
    slow = ("slow",)
    return [
        # demo_microscope_pipeline_single_dataset's parameters
        pl_case(
            "single_dataset",
            demo,
            *chs,
            "save16bit",
            [True, True],
            "overwrite",
            False,
            "deskew",
            True,
            "rotate",
            True,
            *STITCH_ON,
            *FF,
            "minModifyTime",
            1,
            "maxModifyTime",
            10,
            "maxWaitLoopNum",
            10,
            "configFile",
            "",
            "mccMode",
            False,
            marks=slow,
        ),  # fmt: skip
        # demo_microscope_pipeline_single_dataset_nostitch's
        pl_case(
            "nostitch",
            demo,
            *chs,
            "save16bit",
            [True, True],
            "deskew",
            True,
            "rotate",
            True,
            "stitch",
            False,
            "DSRCombined",
            True,
            *FF,
        ),  # fmt: skip
        # demo_microscope_pipeline_multiple_dataset's (two datasets, an image list each)
        pl_case(
            "multiple_datasets",
            two_datasets,
            *cha,
            "save16bit",
            [True, True],
            "stitch",
            True,
            "DSRCombined",
            True,
            "imageListFullpaths",
            ["@d/d1/ImageList.csv", "@d/d2/ImageList.csv"],
            "axisOrder",
            "-x,y,z",
            "FFCorrection",
            False,
            sub=["d1/", "d2/"],
            marks=slow,
        ),  # fmt: skip
        # (AutomaticProcessingTest/imageListsStayWithTheirDatasets) datasets given out of order:
        # the image lists follow the sorted data paths
        pl_case(
            "datasets_out_of_order",
            two_datasets,
            *cha,
            "stitch",
            True,
            "DSRCombined",
            True,
            "imageListFullpaths",
            ["@d/d2/ImageList.csv", "@d/d1/ImageList.csv"],
            "axisOrder",
            "-x,y,z",
            sub=["d2/", "d1/"],
            marks=slow,
        ),  # fmt: skip
        pl_case("separate", one, *cha, "DSRCombined", False),
        pl_case("deskew_only", one, *cha, "rotate", False),
        pl_case("deskew_only_separate", one, *cha, "rotate", False, "DSRCombined", False),
        pl_case("save_zarr_stitch", one, *cha, "saveZarr", True, *STITCH_ON, marks=slow),
        pl_case("partial_files", dataset(cams=("A",), parts=2), *cha),
        pl_case("input_bbox", one, *cha, "inputBbox", [3, 5, 2, 30, 38, 20]),
        pl_case(
            "blend_feather",
            one,
            *cha,
            "blendMethod",
            "feather",
            "edgeArtifacts",
            0,
            *STITCH_ON,
            marks=slow,
        ),  # fmt: skip
        # a factor with rx != ry: the stitching gets it as [y x z]
        pl_case(
            "resample_rx_ry", one, *cha, "resampleFactor", [1.5, 1, 1], *STITCH_ON, marks=slow
        ),  # fmt: skip
        pl_case("resample_xy_isotropic", one, *cha, "resampleType", "xy_isotropic"),
        pl_case("existing", existing, *cha, sub=""),
        pl_case("existing_overwrite", existing, *cha, "overwrite", [False, True, False], sub=""),
        pl_case("flip_zstack", one, *cha, "flipZstack", True, "save16bit", [False, True]),
        # the image list missing: a warning, no stitching, its folder made all the same
        pl_case(
            "image_list_missing", one, *cha, "stitch", True, "imageListFullpaths", "@d/d/none.csv"
        ),  # fmt: skip
        # MATLAB's errors
        pl_case("ff_count", one, *cha, "FFCorrection", True, error=True),
        pl_case(
            "ff_missing",
            one,
            *cha,
            "FFCorrection",
            True,
            "FFImagePaths",
            ["@d/d/x.tif"],
            "backgroundPaths",
            ["@d/d/x.tif"],
            error=True,
        ),  # fmt: skip
        # one save16bit value for both steps (MATLAB's index error until e3fbccc)
        pl_case("save16bit_one_stitch", one, *cha, "save16bit", True, *STITCH_ON),
    ]


CASES = _cases()

# MATLAB on Windows (ispc) writes the stitching's feather distance maps into imdist/ itself, not
# imdist/<first tile>/, and so does the port there: the references that hold such maps, made on
# Linux, hold the other folder and paths. The Windows layout is tested on every system by
# forcing that branch (stitch/test_dataset.py::test_feather_maps_under_the_windows_rule).
_WINDOWS_MAPS = (
    "the feather distance maps go into imdist/ on Windows (MATLAB's ispc branch); this "
    "reference is MATLAB's on Linux (imdist/<first tile>/): "
    "stitch/test_dataset.py::test_feather_maps_under_the_windows_rule covers it"
)


def _maps_in_tile_folder(case: Case) -> bool:
    """whether the case's MATLAB reference holds files in a folder of imdist/"""
    files = oracle.reference(case).files
    return files is not None and any(re.search(r"(^|/)imdist/[^/]+/", n) for n in files.names())


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if sd._is_pc() and _maps_in_tile_folder(case):
        pytest.skip(_WINDOWS_MAPS)
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# AutomaticProcessingTest in intent
# ------------------------------------------------------------------------------------------

LOCAL = {
    "xy_pixel_size": PX,
    "dz": DZ,
    "reverse": True,
    "channel_patterns": ("CamA_ch0",),
    "streaming": False,
    "parse_cluster": False,
    "cpus_per_task": 4,
}
SMALL = (8, 16, 16)  # (z, y, x) of the tiles of the cluster tests


def dp(path) -> str:
    """a data folder with a trailing separator, as the MATLAB tests give them"""
    s = str(path)
    return s if s.endswith(("/", "\\")) else s + "/"


def tile_fn(data, ch: str = "CamA_ch0", yi: int = 0) -> str:
    return (
        f"{dp(data)}Scan_Iter_0000_0000_{ch}_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_"
        f"000x_{yi:03d}y_000z_0000t.tif"
    )


def dsr_of(fn: str) -> str:
    folder, name = os.path.split(fn)
    return f"{folder}/DSR/{name}"


def run(paths, **kw):
    A.automatic_processing(paths, **{**LOCAL, **kw})


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, b.dtype, a.shape, b.shape)
    np.testing.assert_array_equal(a, b)


def two_tiles(root: Path, raw=None):
    """two tiles of one volume, overlapping in y, and their image list (writeTwoTiles)"""
    data = root / "data"
    raw = raw_volume() if raw is None else raw
    _, csv = write_stitch_tiles(dp(data), tiles_of(raw))
    return dp(data), csv


def stitched(data, it: int = 0, cam: str = "A", folder: str = "matlab_stitch") -> str:
    return f"{dp(data)}{folder}/{STITCH.format(it=it, cam=cam)}.zarr"


def test_flat_field_and_deskew_rotate(tmp_path):
    """(flatFieldAndDeskewRotate) two channels with their flat-field and background images: the
    DSR of each equals the file's deskew/rotation with the flat field run directly"""
    from petakit5d.deskew_rotate import deskew_rotate_file

    data, ref = tmp_path / "data", tmp_path / "ref"
    for d in (data, ref):
        dataset(cams=("B", "A"), ff=True, sub="")(d)
    run(dp(data), channel_patterns=("CamB_ch0", "CamA_ch0"), ff_correction=True,
        const_offset=1, ff_image_paths=[f"{data}/FF/ff_B.tif", f"{data}/FF/ff_A.tif"],
        background_paths=[f"{data}/FF/bg_B.tif", f"{data}/FF/bg_A.tif"])  # fmt: skip
    for cam in ("A", "B"):
        fn = tile_fn(ref, f"Cam{cam}_ch0")
        deskew_rotate_file(
            fn, PX, DZ, skew_angle=SKEW, reverse=True, ff_correction=True,
            ff_image=f"{ref}/FF/ff_{cam}.tif", background_image=f"{ref}/FF/bg_{cam}.tif",
            const_offset=1, lower_limit=0.4, rotate=True,
        )  # fmt: skip
        got = read_tiff(dsr_of(tile_fn(data, f"Cam{cam}_ch0")))
        assert got.dtype == np.uint16
        same(got, read_tiff(dsr_of(fn)))
        name = os.path.basename(fn)[:-4]
        assert os.path.isfile(f"{data}/DSR/MIPs/{name}_MIP_z.tif")


def test_demo_parameter_sets_are_accepted(tmp_path):
    """(demoParameterSetIsAccepted) the demos' Save16bit (they once passed 4 values) is taken"""
    from petakit5d import matlab
    from petakit5d.params import literal
    from petakit5d.params.values import to_plain

    demos = Path(__file__).resolve().parents[3] / "demos"
    names = ("demo_microscope_pipeline_single_dataset",
             "demo_microscope_pipeline_single_dataset_nostitch",
             "demo_microscope_pipeline_multiple_dataset")  # fmt: skip
    if not all((demos / f"{n}.m").is_file() for n in names):
        pytest.skip("the MATLAB demos are not next to python/")
    data = tmp_path / "data"
    data.mkdir()
    synth.write_tiff(tile_fn(data), raw_volume())
    for n in names:
        m = re.search(r"(?m)^\s*Save16bit\s*=\s*([^;\r\n]+);", (demos / f"{n}.m").read_text())
        assert m, f"{n} no longer sets Save16bit."
        save16 = to_plain(literal.parse(m.group(1)))
        matlab.XR_microscopeAutomaticProcessing(
            dp(data), "xyPixelSize", PX, "dz", DZ, "Reverse", True, "ChannelPatterns",
            ["CamA_ch0"], "Save16bit", save16, "Overwrite", False, "Streaming", False,
            "Deskew", True, "Rotate", True, "Stitch", False, "DSRCombined", True,
            "minModifyTime", 1, "maxModifyTime", 10, "maxWaitLoopNum", 10, "parseCluster", False,
            "cpusPerTask", 4, "configFile", "", "mccMode", False,
        )  # fmt: skip
    assert os.path.isfile(dsr_of(tile_fn(data)))


# -- cluster mode against the fake Slurm ----------------------------------------------------


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


def small_tiles(data: Path, n: int) -> list[str]:
    data.mkdir(parents=True, exist_ok=True)
    fns = [tile_fn(data, yi=y) for y in range(n)]
    for fn in fns:
        synth.write_tiff(fn, np.zeros(SMALL, np.uint16))
    return fns


def cluster(paths, logs: Path, **kw):
    logs.mkdir(exist_ok=True)
    kw = {"parse_cluster": True, "job_log_dir": str(logs), **kw}
    run(paths, **kw)


def wait_pauses(pauses, n: int) -> None:
    """(verifyWaitsOfASecond) the pipeline went on within a second of the results: n waits of a
    second, the last of which saw them"""
    assert pauses == [1.0] * n, pauses


@POSIX
def test_cluster_mode_wait_grows_with_unfinished_files(tmp_path, fake_slurm, monkeypatch):
    """(clusterModeWaitGrowsWithUnfinishedFiles) 0.3 s per unfinished file, for their results"""
    fns = small_tiles(tmp_path / "data", 20)
    calls = []
    monkeypatch.setattr(jfiles, "wait_for_files", lambda g, t: calls.append((g, t)) or False)
    cluster(dp(tmp_path / "data"), tmp_path / "logs", max_trial_num=2)
    assert fake_slurm.call_count("sbatch") == 40
    assert [t for _, t in calls] == pytest.approx([6.0, 6.0])
    awaited = [re.sub("/+", "/", p) for p in calls[0][0][0]]
    assert sorted(awaited) == sorted(re.sub("/+", "/", dsr_of(f)) for f in fns)


@POSIX
def test_cluster_mode_ends_soon_after_the_last_results(tmp_path, fake_slurm, pauses, monkeypatch):
    """(clusterModeEndsSoonAfterTheLastResults) the results come during the second second"""
    fns = small_tiles(tmp_path / "data", 2)

    def pause(t):
        pauses.append(t)
        if len(pauses) == 2:
            for f in fns:
                Path(dsr_of(f)).touch()

    monkeypatch.setattr(jfiles, "_pause", pause)
    cluster(dp(tmp_path / "data"), tmp_path / "logs", max_trial_num=2)
    assert fake_slurm.call_count("sbatch") == 2, "A job was submitted again."
    wait_pauses(pauses, 2)


@POSIX
def test_cluster_mode_does_not_wait_for_steps_given_up(tmp_path, fake_slurm, pauses, monkeypatch):
    """(clusterModeDoesNotWaitForStepsGivenUp) tile 1's job is gone at the second round, tile 2's
    runs, and its result comes during the first second of the second wait"""
    fns = small_tiles(tmp_path / "data", 2)
    fake_slurm.set_job_states("1001_2", ["R"])

    def pause(t):
        pauses.append(t)
        if len(pauses) == 4:
            Path(dsr_of(fns[1])).touch()

    monkeypatch.setattr(jfiles, "_pause", pause)
    cluster(dp(tmp_path / "data"), tmp_path / "logs", max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 2
    wait_pauses(pauses, 4)


def _claim(logs: Path, i: int, job: str, kind: str | None = None) -> None:
    for c in logs.glob(f"{single.POOL_PREFIX}*/claims"):
        if kind is not None and kind not in (c.parent / "tasks.ndjson").read_text():
            continue
        try:
            os.symlink(f"{job}:othernode:5:1", c / f"t{i}")
        except FileExistsError:
            pass


@LINUX
@POSIX
def test_last_try_waits_for_a_task_runner_that_has_the_task(tmp_path, fake_slurm, monkeypatch):
    """(lastTryWaitsForATaskRunnerThatHasTheTask) tile 2's job runs tile 3's task too"""
    fns = small_tiles(tmp_path / "data", 3)
    logs = tmp_path / "logs"
    fake_slurm.set_job_states("1000_1", ["R", "done"], [dsr_of(fns[0])])
    fake_slurm.set_job_states("1001_2", ["R", "R", "R", "R", "R", "done"], [dsr_of(fns[1])])
    monkeypatch.setattr(jfiles, "_pause", lambda t: _claim(logs, 2, "1001_2"))
    cluster(dp(tmp_path / "data"), logs, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 3, "Tile 3 was submitted again."
    checks = [c for c in fake_slurm.history() if c.startswith("squeue") and "1001_2" in c]
    assert len(checks) >= 6, "The pipeline did not wait for the job that had tile 3."
    assert not list(logs.glob(f"{single.POOL_PREFIX}*")), "The pool is still there."


@LINUX
@POSIX
def test_stitching_last_try_waits_for_a_task_runner_that_has_the_task(
    tmp_path, fake_slurm, monkeypatch
):
    """(stitchingLastTryWaitsForATaskRunnerThatHasTheTask) the second dataset's stitching job
    runs the third's stitching task too"""
    paths, csvs = [], []
    for d in range(3):
        data = tmp_path / f"d{d}"
        tiles = [{"im": np.zeros(SMALL, np.uint16), "idx": (0, 0, 0), "xyz": (0, 0, 0)}]
        _, csv = write_stitch_tiles(dp(data), tiles)
        paths.append(dp(data))
        csvs.append(csv)
    logs = tmp_path / "logs"
    fsn = STITCH.format(it=0, cam="A")
    for d in range(3):
        fake_slurm.set_job_states(
            f"{1000 + 2 * d}_{d + 1}", ["R", "done"],
            [f"{paths[d]}DSR/{fsn}_000x_000y_000z_0000t.tif"],
        )  # fmt: skip
    fake_slurm.set_job_states("1001_1", ["R", "done"], [f"{paths[0]}matlab_stitch/{fsn}.zarr"])
    fake_slurm.set_job_states(
        "1003_2", ["R"] * 5 + ["done"], [f"{paths[1]}matlab_stitch/{fsn}.zarr"]
    )
    monkeypatch.setattr(jfiles, "_pause", lambda t: _claim(logs, 2, "1003_2", "stitch_dataset"))
    cluster(paths, logs, stitch=True, image_list_fullpaths=csvs, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 6, "A job was submitted again."
    checks = [c for c in fake_slurm.history() if c.startswith("squeue") and "1003_2" in c]
    assert len(checks) >= 6, "The pipeline did not wait for the job that had the third stitching."


@POSIX
def test_cluster_mode_sees_results_other_nodes_just_made(tmp_path, fake_slurm, monkeypatch, pauses):
    """(clusterModeSeesResultsOtherNodesJustMade) batch_file_exist sees them: no job"""
    small_tiles(tmp_path / "data", 1)
    monkeypatch.setattr(jfiles, "batch_file_exist", lambda paths, *a, **k: [True] * len(paths))
    cluster(dp(tmp_path / "data"), tmp_path / "logs", max_trial_num=2)
    assert fake_slurm.call_count("sbatch") == 0, "A deskew job was submitted for a finished result."


@pytest.mark.slow
def test_deskew_rotate_and_stitch_multiple_datasets(tmp_path):
    """(deskewRotateAndStitchMultipleDatasets) two datasets, an image list each: the stitched
    DSR tiles equal the DSR of the whole volume (with edge_artifacts 0; the deskew's uint16 sums
    depend on the y lines' lengths, so the tiles' and the volume's y sizes are multiples of 16,
    as in MATLAB's test)"""
    from petakit5d.deskew_rotate import deskew_rotate_volume

    paths, csvs, want = [], [], []
    for d in range(2):
        raw = raw_volume(seed=3 + d, shape=(40, 80, 64))
        want.append(deskew_rotate_volume(raw, SKEW, DZ, PX, reverse=True, save_16bit=True))
        tiles = cut_tiles(raw, [(0, 0, 0), (0, 32, 0)], (40, 48, 64), (PX * 2, PX, PX))
        data, csv = dp(tmp_path / f"ds{d}"), write_stitch_tiles(dp(tmp_path / f"ds{d}"), tiles)[1]
        paths.append(data)
        csvs.append(csv)
    run(paths, stitch=True, image_list_fullpaths=csvs, axis_order="-x,y,z", edge_artifacts=0)
    for d in range(2):
        fn = stitched(paths[d])
        assert os.path.isdir(fn), f"Stitched result {fn} is missing."
        out = read_zarr(fn)
        sz = want[d].shape
        assert out.shape == tuple(s + 2 for s in sz)
        same(out[: sz[0], : sz[1], : sz[2]], want[d])


@pytest.mark.slow
def test_local_stitching_finishes_without_waiting(tmp_path, monkeypatch):
    """(localStitchingFinishesWithoutWaiting) no wait of 30 s after the stitching"""
    data, csv = two_tiles(tmp_path)
    calls = []
    monkeypatch.setattr(jfiles, "_pause", lambda t: calls.append(t))
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z")
    assert os.path.isdir(stitched(data)), "The stitched result is missing."
    assert not any(t >= 30 for t in calls), f"The pipeline waited {calls} s."


def _stitch_args(csv, n):
    return {"stitch": True, "image_list_fullpaths": [csv], "max_trial_num": n}


@POSIX
def test_cluster_mode_waits_for_the_last_stitching_try(tmp_path, fake_slurm, monkeypatch):
    """(clusterModeWaitsForTheLastStitchingTry) checked until it ends, not submitted again"""
    data, csv = two_tiles(tmp_path)
    fake_slurm.set_job_states("1002_1", ["PD", "R", "R", "R", "done"])
    monkeypatch.setattr(jfiles, "_pause", lambda t: None)
    cluster(data, tmp_path / "logs", **_stitch_args(csv, 1))
    assert fake_slurm.call_count("sbatch") == 3, "A last try was submitted again."
    checks = [c for c in fake_slurm.history() if c.startswith("squeue") and "1002_1" in c]
    assert len(checks) >= 5, "The stitching job was not checked until it ended."


STITCH_FUNC = A.TASK_FUNCS["XR_matlab_stitching_wrapper"]


@POSIX
def test_cluster_mode_submits_stitching_after_the_deskew_jobs(
    tmp_path, fake_slurm, monkeypatch, submitted
):
    """(clusterModeSubmitsStitchingAfterTheDeskewJobs)"""
    data, csv = two_tiles(tmp_path)
    monkeypatch.setattr(jfiles, "_pause", lambda t: None)
    cluster(data, tmp_path / "logs", **_stitch_args(csv, 1))
    funcs = [r["func"] for r in submitted]
    assert len(funcs) == 3
    assert [i for i, f in enumerate(funcs) if f == STITCH_FUNC] == [2], (
        "The stitching job was not submitted after the deskew jobs."
    )
    # (the stitching job streams only when the pipeline does)
    assert submitted[2]["kwargs"]["streaming"] is False
    assert submitted[2]["kwargs"]["parse_cluster"] is True


@POSIX
def test_jobs_ask_for_matlabs_memory(tmp_path, fake_slurm, monkeypatch):
    """the combined deskew: size / 2^30 * 2 * (7 + 3 / prod(resample) + 8 * tries) GB; the
    separate one XR_estimateComputingMemory's; the stitching size / 2^30 * 2 * 20 (2 for MIP
    stitching); the stitching job's array index is the dataset's first file's"""
    from petakit5d.tools.memory_estimate import estimate_computing_memory

    data, csv = two_tiles(tmp_path)
    calls = []
    real = single.submit_single_job

    def spy(task, job_id, task_id, **kw):
        calls.append((task.func, task_id, kw["mem_allocate"], kw["cpus_per_task"]))
        return real(task, job_id, task_id, **kw)

    monkeypatch.setattr(single, "submit_single_job", spy)
    monkeypatch.setattr(jfiles, "_pause", lambda t: None)
    cluster(data, tmp_path / "logs", **_stitch_args(csv, 2), resample_factor=(2, 1, 1))
    gb = float(np.prod((RAW[0], TILE_Y, RAW[2]))) * 2 / 2**30
    deskew = [c for c in calls if c[0] != STITCH_FUNC]
    stitch = [c for c in calls if c[0] == STITCH_FUNC]
    assert [c[1] for c in deskew] == [1, 2, 1, 2]
    assert [c[2] for c in deskew] == pytest.approx([gb * 2 * 8.5] * 2 + [gb * 2 * 16.5] * 2)
    assert [c[1] for c in stitch] == [1, 1]
    assert [c[2] for c in stitch] == pytest.approx([gb * 2 * 20] * 2)
    assert all(c[3] == 4 for c in calls)
    calls.clear()
    cluster(data, tmp_path / "logs2", max_trial_num=1, dsr_combined=False, overwrite=True)
    est = float(estimate_computing_memory(tile_fn(data), ("deskew",))[0][0])
    assert [c[2] for c in calls] == pytest.approx([est, est])


@POSIX
def test_cluster_mode_submits_stitching_when_tiles_were_deskewed_before(
    tmp_path, fake_slurm, monkeypatch, submitted
):
    """(clusterModeSubmitsStitchingWhenTilesWereDeskewedBefore)"""
    data, csv = two_tiles(tmp_path)
    os.makedirs(f"{data}DSR")
    Path(dsr_of(tile_fn(data))).touch()
    monkeypatch.setattr(jfiles, "_pause", lambda t: None)
    cluster(data, tmp_path / "logs", **_stitch_args(csv, 1))
    funcs = [r["func"] for r in submitted]
    assert len(funcs) == 2 and funcs[1] == STITCH_FUNC


@POSIX
def test_cluster_mode_looks_for_the_stitched_result_each_second(
    tmp_path, fake_slurm, monkeypatch, pauses
):
    """(clusterModeLooksForTheStitchedResultEachSecond) its temporary folder comes during the
    second second (it doesn't count), the result during the third"""
    data, csv = two_tiles(tmp_path)
    os.makedirs(f"{data}DSR")
    for y in (0, 1):
        Path(dsr_of(tile_fn(data, yi=y))).touch()
    fn = stitched(data)

    def pause(t):
        pauses.append(t)
        if len(pauses) == 2:
            os.mkdir(fn + "_1234")
        if len(pauses) == 3:
            os.mkdir(fn)

    monkeypatch.setattr(jfiles, "_pause", pause)
    cluster(data, tmp_path / "logs", **_stitch_args(csv, 1))
    assert fake_slurm.call_count("sbatch") == 1
    wait_pauses(pauses, 3)


def _age_then_back(marker: str, ages: list[float], stop: bool = False) -> None:
    """the marker's age in s (NaN: none), then 5 minutes added to it (recordAndAgeFile)"""
    if os.path.exists(marker):
        ages.append(time.time() - os.path.getmtime(marker))
        old = time.time() - 300
        os.utime(marker, (old, old))
    else:
        ages.append(float("nan"))
    if stop:
        raise KeyboardInterrupt("Stopped by the test.")


@POSIX
def test_cluster_mode_tells_stitching_that_tiles_are_coming(tmp_path, fake_slurm, monkeypatch):
    """(clusterModeTellsStitchingThatTilesAreComing) the marker is kept (and written again once
    a minute old) while deskew jobs are pending or running, removed after and at the end"""
    data, csv = two_tiles(tmp_path)
    marker = f"{data}DSR/deskew_in_progress.tmp"
    fake_slurm.set_job_states("1002_1", ["R"] * 6 + ["done"])
    ages: list[float] = []
    monkeypatch.setattr(jfiles, "wait_for_files", lambda g, t: _age_then_back(marker, ages))
    cluster(data, tmp_path / "logs", **_stitch_args(csv, 2))
    assert len(ages) >= 4
    assert all(a < 60 for a in ages[:3]), f"The marker was missing or not written again: {ages}"
    assert all(np.isnan(ages[3:])), f"The marker was kept after the deskew tries: {ages}"
    assert not os.path.exists(marker), "The marker was left after the run."


@POSIX
def test_deskew_marker_is_removed_when_the_pipeline_stops(tmp_path, fake_slurm, monkeypatch):
    """(deskewMarkerIsRemovedWhenThePipelineStops) an error (or Ctrl-C) leaves no marker"""
    data, csv = two_tiles(tmp_path)
    marker = f"{data}DSR/deskew_in_progress.tmp"
    ages: list[float] = []
    monkeypatch.setattr(jfiles, "_pause", lambda t: _age_then_back(marker, ages, True))
    with pytest.raises(KeyboardInterrupt):
        cluster(data, tmp_path / "logs", **_stitch_args(csv, 2))
    assert ages and ages[0] < 60, "There was no marker while the deskew jobs were pending."
    assert not os.path.exists(marker), "The marker was left after the pipeline stopped."


# -- streaming --------------------------------------------------------------------------------


def _age(fn: str, minutes: float) -> None:
    t = time.time() - minutes * 60
    os.utime(fn, (t, t))


def test_streaming_leaves_out_recent_files_with_their_sizes(tmp_path):
    """(streamingLeavesOutRecentFilesWithTheirSizes) files changed in the last min_modify_time
    minutes are left for later, their sizes too"""
    from petakit5d.pipeline.filenames import parse_dataset_filenames

    data = tmp_path / "data"
    data.mkdir()
    fns = [tile_fn(data, yi=y) for y in range(3)]
    for fn, depth in zip(fns, (16, 20, 24), strict=True):
        synth.write_tiff(fn, np.zeros((depth, 16, 16), np.uint16))
    args = ([dp(data)], ["CamA_ch0"], False, False, False, [])
    every = parse_dataset_filenames(*args, False, 1, False)
    p = parse_dataset_filenames(*args, True, 1, False)
    assert p.fnames == [] and p.gfnames == [] and len(p.partialvols) == 0
    assert len(p.data_sizes) == 0, "The sizes of the files left out stayed in the list."
    _age(fns[2], 5)
    p = parse_dataset_filenames(*args, True, 1, False)
    name = os.path.basename(fns[2])
    assert p.fnames == [name]
    assert p.gfnames == [[name]], "The file got the group of another file."
    assert len(p.partialvols) == 1
    assert list(p.data_sizes) == [every.data_sizes[2]], "The file got the size of another file."


def test_streaming_processes_new_files(tmp_path, monkeypatch):
    """(streamingProcessesNewFiles) a file that arrives during the run (at the first wait) is
    processed"""
    data, stage = tmp_path / "data", tmp_path / "stage"
    data.mkdir()
    stage.mkdir()
    raw = raw_volume()
    synth.write_tiff(tile_fn(data), raw)
    staged, new = tile_fn(stage, yi=1), tile_fn(data, yi=1)
    synth.write_tiff(staged, raw)
    _age(tile_fn(data), 5)
    _age(staged, 5)
    calls = []

    def pause(t):
        calls.append(t)
        if os.path.exists(staged):
            shutil.move(staged, new)

    monkeypatch.setattr(jfiles, "_pause", pause)
    run(dp(data), streaming=True, max_modify_time=1, max_wait_loop_num=2)
    assert os.path.isfile(dsr_of(new)), "The new file was not processed."
    assert os.path.isfile(dsr_of(tile_fn(data)))
    assert calls == [30, 30], calls


@POSIX
@pytest.mark.slow
def test_local_run_submits_no_stitching_jobs(tmp_path, fake_slurm):
    """(localRunSubmitsNoStitchingJobs) the stitching gets parse_cluster False"""
    data, csv = two_tiles(tmp_path)
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z")
    assert fake_slurm.call_count("sbatch") == 0, "The local run submitted Slurm jobs."
    assert os.path.isdir(stitched(data)), "The stitched result is missing."


# -- the function strings, their values and MATLAB's quirks -----------------------------------


@pytest.fixture
def strings(monkeypatch):
    """the function strings the pipeline runs here, in order (and runs them)"""
    out: list[str] = []
    real = A._State._run_here

    def spy(self, func_str):
        out.append(func_str)
        return real(self, func_str)

    monkeypatch.setattr(A._State, "_run_here", spy)
    return out


def test_values_reach_the_steps_as_matlab_formats_them(tmp_path, strings):
    """%.20d of the pixel size, z step and angle, mat2str of the box and of the deskew's [x y z]
    resample factor, %.10d of the stitching's [y x z] (MATLAB's 0000000001 for 1); the stitching
    of a local run gets one try"""
    data, csv = two_tiles(tmp_path)
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z",
        resample_factor=(1, 1, 1.5), max_trial_num=3)  # fmt: skip
    deskew = [s for s in strings if s.startswith("XR_deskewRotateFrame(")]
    stitch = [s for s in strings if s.startswith("XR_matlab_stitching_wrapper(")]
    assert len(deskew) == 2 and stitch
    assert ",1.07999999999999998779e-01,2.99999999999999988898e-01,'skewAngle'," in deskew[0]
    assert "3.24500000000000028422e+01" in deskew[0]
    assert "'resampleFactor',[[1.5,1,1]]" in deskew[0] and "'inputBbox',[]" in deskew[0]
    s = stitch[0]
    assert "'resampleFactor',[0000000001,1.5000000000e+00,0000000001,]" in s
    assert "'maxTrialNum',1," in s
    assert "'ResultDir','matlab_stitch'" in s and "'ProcessedDirStr','DSR'" in s
    assert "'resampleType','isotropic'" in s and "'stitchMIP',[]" in s
    # read back: the deskew's [x y z] factor is (z, y, x) (1, 1, 1.5), and so is the stitching's
    # (AutomaticProcessingTest/stitchingGetsTheResampleFactorInItsAxisOrder; until e3fbccc the
    # stitching read the deskew's text as [y x z]: its y factor was the deskew's x factor)
    d = A.task_of(deskew[0]).kwargs
    assert d["xy_pixel_size"] == PX and d["dz"] == DZ and d["skew_angle"] == SKEW
    assert tuple(d["resample_factor"]) == (1, 1, 1.5)
    k = A.task_of(s).kwargs
    assert tuple(k["resample_factor"]) == (1, 1, 1.5)
    assert k["result_dir_name"] == "matlab_stitch" and k["processed_dir_str"] == "DSR"


def test_the_input_box_reaches_the_deskew_as_a_python_box(tmp_path, strings):
    data = tmp_path / "data"
    dataset(cams=("A",), sub="")(data)
    run(dp(data), input_bbox=(1, 2, 4, 20, 30, 38))
    assert "'inputBbox',[3,5,2,30,38,20]," in strings[0]
    assert tuple(A.task_of(strings[0]).kwargs["input_bbox"]) == (1, 2, 4, 20, 30, 38)


def test_image_lists_stay_with_their_datasets(tmp_path, strings, monkeypatch):
    """(AutomaticProcessingTest/imageListsStayWithTheirDatasets) the data paths are sorted
    (unique) and the image lists follow them (until e3fbccc they kept the given order, so
    datasets given out of order got each other's lists)"""
    paths, csvs = [], []
    for d in ("b", "a"):
        data, csv = two_tiles(tmp_path / d)
        paths.append(data)
        csvs.append(csv)
    run(paths, stitch=True, image_list_fullpaths=csvs, axis_order="-x,y,z")
    stitch = [s for s in strings if s.startswith("XR_matlab_stitching_wrapper(")]
    assert stitch
    for s in stitch:
        data, csv = A._split_call(s)[1][:2]
        assert csv == f"{data}ImageList.csv", f"Dataset {data} got another image list."


def test_a_file_flip_is_the_default_of_files_found_later(tmp_path, monkeypatch):
    """(MATLAB's quirk) flipZstack takes the last deskewed file's flip, and the streaming look
    for new files gets it"""
    from petakit5d.pipeline import filenames

    data = tmp_path / "data"
    data.mkdir()
    synth.write_tiff(tile_fn(data), np.zeros(SMALL, np.uint16))
    _age(tile_fn(data), 5)
    real = filenames.parse_dataset_filenames
    flips = []

    def spy(*a, **k):
        flips.append(a[3])
        r = real(*a, **k)
        r.flip_zstack[:] = True
        return r

    monkeypatch.setattr(filenames, "parse_dataset_filenames", spy)
    monkeypatch.setattr(jfiles, "_pause", lambda t: None)
    run(dp(data), streaming=True, max_modify_time=1, max_wait_loop_num=1)
    assert flips[0] is False and flips[1:] and all(f is True for f in flips[1:]), flips


def test_streaming_adds_new_files_named_as_another_datasets_files(tmp_path, monkeypatch):
    """(AutomaticProcessingTest/streamingAddsNewFilesNamedAsAnotherDatasetsFiles) a new file is
    one with no listed file of its name in its data folder (until e3fbccc a file was left out if
    any listed name of any dataset contained its name)"""
    d1, d2, stage = tmp_path / "d1", tmp_path / "d2", tmp_path / "stage"
    for d in (d1, d2, stage):
        d.mkdir()
    for fn in (tile_fn(d1), tile_fn(d2, yi=1), tile_fn(stage)):
        synth.write_tiff(fn, np.zeros(SMALL, np.uint16))
        _age(fn, 5)

    def pause(t):
        if os.path.exists(tile_fn(stage)):
            shutil.move(tile_fn(stage), tile_fn(d2))

    monkeypatch.setattr(jfiles, "_pause", pause)
    run([dp(d1), dp(d2)], streaming=True, max_modify_time=1, max_wait_loop_num=2)
    assert os.path.isfile(dsr_of(tile_fn(d1)))
    assert os.path.isfile(dsr_of(tile_fn(d2))), (
        "The new file of the second dataset was not processed."
    )


def test_the_axis_order_check_never_raises(tmp_path):
    """(MATLAB's quirk) ~regexpi(...) is false for a match and for none"""
    data = tmp_path / "data"
    dataset(cams=("A",), sub="")(data)
    with pytest.warns(UserWarning, match="does not exist"):
        run(dp(data), stitch=True, image_list_fullpaths=[f"{data}/none.csv"], axis_order="q")
    assert os.path.isdir(f"{data}/matlab_stitch"), "the stitching folder is made all the same"


def test_resample_settings():
    """XR_checkResampleSetting: 'given' completes one or two values; 'xy_isotropic' makes the
    rotated voxels isotropic; an empty 'given' fails"""
    assert A.check_resample_setting("given", [2], False, SKEW, PX, DZ)[0] == [2, 2, 2]
    assert A.check_resample_setting("given", [2, 3], False, SKEW, PX, DZ)[0] == [2, 2, 3]
    assert A.check_resample_setting("isotropic", None, False, SKEW, PX, DZ)[0] == [1, 1, 1]
    r, z_aniso = A.check_resample_setting("xy_isotropic", None, False, SKEW, PX, DZ)
    assert r[:2] == [1, 1] and r[2] == pytest.approx(1.17615334908839, rel=1e-14)
    assert z_aniso == pytest.approx(np.sin(np.deg2rad(SKEW)) * DZ / PX, rel=1e-15)
    assert A.check_resample_setting("xy_isotropic", None, True, SKEW, PX, DZ)[1] == DZ / PX
    with pytest.raises(A.AutomaticProcessingError, match="must not be empty"):
        A.check_resample_setting("given", None, False, SKEW, PX, DZ)


def test_ff_images_are_picked_by_contains(tmp_path):
    """a file that no channel pattern contains (a pattern that matched as an expression) has
    no flat field: MATLAB's error"""
    data = tmp_path / "data"
    dataset(cams=("A",), ff=True, sub="")(data)
    with pytest.raises(A.AutomaticProcessingError) as e:
        run(dp(data), channel_patterns=("CamA.*ch0",), ff_correction=True,
            ff_image_paths=[f"{data}/FF/ff_A.tif"],
            background_paths=[f"{data}/FF/bg_A.tif"])  # fmt: skip
    assert e.value.identifier == "MATLAB:needMoreRhsOutputs"


def test_young_tmp_marker_is_another_process_file(tmp_path, monkeypatch):
    """(local mode) a file whose .tmp marker is younger than unit_wait_time minutes is left;
    the round waits 30 s, then takes it over"""
    data = tmp_path / "data"
    small_tiles(data, 1)
    os.makedirs(f"{data}/DSR")
    marker = f"{data}/DSR/{os.path.basename(tile_fn(data))[:-4]}.tmp"
    Path(marker).touch()
    waits = []

    def wait(groups, t):
        waits.append(t)
        os.utime(marker, (1, 1))
        return False

    monkeypatch.setattr(jfiles, "wait_for_files", wait)
    run(dp(data))
    assert waits == [30.0]
    assert os.path.isfile(dsr_of(tile_fn(data))) and not os.path.exists(marker)


def test_a_file_deleted_during_the_run_is_finished(tmp_path, monkeypatch):
    data = tmp_path / "data"
    fns = small_tiles(data, 2)
    real = A._State._deskew

    def deskew(self, f, frame, *a):
        if f == 0:
            os.remove(fns[1])
        return real(self, f, frame, *a)

    monkeypatch.setattr(A._State, "_deskew", deskew)
    run(dp(data))
    assert os.path.isfile(dsr_of(fns[0])) and not os.path.exists(dsr_of(fns[1]))


@pytest.mark.cluster
@pytest.mark.slow
def test_cluster_mode_gives_the_local_trees():
    """the demo's parameters (two cameras, two time points, flat field, stitching) in local mode
    and in cluster mode (a job per file and one stitching job, the master computing none): the
    same DSR files, MIPs, stitched volumes and their MIPs"""
    import uuid

    scratch = os.environ.get("PETAKIT5D_TEST_SLURM_SCRATCH")
    if not scratch or not os.path.isdir(scratch):
        pytest.skip("PETAKIT5D_TEST_SLURM_SCRATCH names no folder")
    root = Path(scratch) / f"pk5d_pipeline_test_{uuid.uuid4().hex[:12]}"
    try:
        trees = []
        for mode in ("local", "cluster"):
            d = root / mode
            dataset(cams=("B", "A"), iters=(0, 1), ff=True, sub="")(d)
            logs = root / f"logs_{mode}"
            logs.mkdir(parents=True)
            run(dp(d), channel_patterns=("CamB_ch0", "CamA_ch0"), stitch=True,
                image_list_fullpaths=[f"{d}/ImageList.csv"], axis_order="-x,y,z",
                ff_correction=True, const_offset=1,
                ff_image_paths=[f"{d}/FF/ff_B.tif", f"{d}/FF/ff_A.tif"],
                background_paths=[f"{d}/FF/bg_B.tif", f"{d}/FF/bg_A.tif"],
                parse_cluster=mode == "cluster", master_compute=False,
                config_file=os.environ["PETAKIT5D_TEST_SLURM_CONFIG"],
                job_log_dir=str(logs))  # fmt: skip
            trees.append(d)
        logs = [f.read_text() for f in (root / "logs_cluster").rglob("job_*.out")]
        assert sum(t.count("Deskew, rotate and resample") for t in logs) >= 8, logs
        a, b = _results(trees[0]), _results(trees[1])
        # 8 DSR files and their MIPs, 4 stitched volumes and their MIPs
        assert sorted(a) == sorted(b) and len(a) == 24
        for rel in a:
            if rel.endswith(".zarr"):
                same(read_zarr(str(a[rel])), read_zarr(str(b[rel])))
            elif rel.endswith(".tif"):
                same(read_tiff(str(a[rel])), read_tiff(str(b[rel])))
    finally:
        shutil.rmtree(root, ignore_errors=True)


def _results(root: Path) -> dict[str, Path]:
    """the DSR files and MIPs, the stitched volumes and their MIPs"""
    out = {}
    for sub in ("DSR", "DSR/MIPs", "matlab_stitch", "matlab_stitch/MIPs"):
        for p in sorted((root / sub).glob("*")):
            if p.suffix in (".tif", ".zarr"):
                out[f"{sub}/{p.name}"] = p
    return out


def test_stitch_prefix_leaves_out_39_characters():
    name = (
        "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y"
    )
    name += "_000z_0000t"
    assert A._stitch_prefix(name) == name[: name.index("_0106060251")]
    assert A._stitch_prefix("short") == ""


def test_more_image_lists_than_datasets_turn_stitching_off(tmp_path):
    data, csv = two_tiles(tmp_path)
    with pytest.warns(UserWarning, match="does not match that of the data paths"):
        run(data, stitch=True, image_list_fullpaths=[csv, csv], axis_order="-x,y,z")
    assert not os.path.exists(stitched(data)) and os.path.isfile(dsr_of(tile_fn(data)))


def test_stitching_gets_the_second_save_16bit(tmp_path, strings):
    data, csv = two_tiles(tmp_path)
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z",
        save_16bit=(True, False), edge_artifacts=0)  # fmt: skip
    deskew = [x for x in strings if x.startswith("XR_deskewRotateFrame(")]
    stitch = [x for x in strings if x.startswith("XR_matlab_stitching_wrapper(")]
    assert A.task_of(deskew[0]).kwargs["save_16bit"] is True
    assert A.task_of(stitch[0]).kwargs["save_16bit"] is False


def test_scalar_save16bit_is_for_both_steps(tmp_path, strings):
    """(AutomaticProcessingTest/scalarSave16bitIsForBothSteps) one save16bit value is for both
    steps (the stitching string took its second, an index error, until e3fbccc)"""
    data, csv = two_tiles(tmp_path)
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z", save_16bit=False,
        edge_artifacts=0)  # fmt: skip
    stitch = [x for x in strings if x.startswith("XR_matlab_stitching_wrapper(")]
    assert stitch and A.task_of(stitch[0]).kwargs["save_16bit"] is False
    assert read_tiff(dsr_of(tile_fn(data))).dtype == np.float32


def test_local_stitching_without_result_is_tried_once(tmp_path, monkeypatch):
    """(AutomaticProcessingTest/localStitchingWithoutResultIsTriedOnce) locally a try counts once
    the dataset's files are deskewed, so that a stitching that makes no result ends (it was
    called again in every round until e3fbccc)"""
    data, csv = two_tiles(tmp_path)
    stitches = []
    real = A._State._run_here

    def spy(self, func_str):
        if func_str.startswith("XR_matlab_stitching_wrapper("):
            stitches.append(func_str)
            return None
        return real(self, func_str)

    monkeypatch.setattr(A._State, "_run_here", spy)
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z")
    # (after each file's deskew; the second one counts)
    assert len(stitches) == 2


def test_the_result_is_found_with_any_msec_field(tmp_path, strings):
    """(AutomaticProcessingTest/clusterModeFindsTheResultOfTilesWithOtherMsecFields, here
    locally) a stitched volume has its first tile's msec field, which another file of the
    volume may not have: the look for the result takes any (until e3fbccc each file looked for
    its own)"""
    data, csv = two_tiles(tmp_path)
    names = [os.path.basename(tile_fn(data, yi=k))[:-4] for k in (0, 1)]
    text = Path(csv).read_bytes()
    os.makedirs(f"{data}DSR")
    for k, n in enumerate(names):
        new = n.replace("_0000000msec_", f"_{5000 + 1000 * k:07d}msec_")
        os.rename(f"{data}{n}.tif", f"{data}{new}.tif")
        text = text.replace(n.encode(), new.encode())
        Path(f"{data}DSR/{new}.tif").touch()
    Path(csv).write_bytes(text)
    os.makedirs(stitched(data).replace("_0000000msec_", "_0005000msec_"))
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z")
    assert not [x for x in strings if x.startswith("XR_matlab_stitching_wrapper(")]


def test_a_stitched_result_removes_its_marker(tmp_path):
    """the stitched result's marker is MATLAB's <name>..tmp"""
    data, csv = two_tiles(tmp_path)
    os.makedirs(f"{data}matlab_stitch")
    fn = stitched(data)
    os.mkdir(fn)
    marker = fn[:-4] + ".tmp"
    assert marker.endswith("Abs..tmp")
    Path(marker).touch()
    run(data, stitch=True, image_list_fullpaths=[csv], axis_order="-x,y,z")
    assert not os.path.exists(marker)


def test_local_tries_end_a_file_without_result(tmp_path, monkeypatch):
    """a file whose deskew makes no result is tried max_trial_num times, 30 s apart (its own
    .tmp marker made older at each wait, so that it doesn't count as another process's)"""
    data = tmp_path / "data"
    fns = small_tiles(data, 1)
    marker = f"{data}/DSR/{os.path.basename(fns[0])[:-4]}.tmp"
    runs, waits = [], []

    def wait(groups, t):
        waits.append(t)
        os.utime(marker, (1, 1))

    monkeypatch.setattr(A._State, "_run_here", lambda self, s: runs.append(s))
    monkeypatch.setattr(jfiles, "wait_for_files", wait)
    run(dp(data), max_trial_num=3)
    assert len(runs) == 3 and waits == [30.0, 30.0]


@POSIX
def test_a_pending_last_try_is_waited_for(tmp_path, fake_slurm, monkeypatch):
    fns = small_tiles(tmp_path / "data", 1)
    fake_slurm.set_job_states("1000_1", ["PD", "PD", "R", "done"], [dsr_of(fns[0])])
    monkeypatch.setattr(jfiles, "_pause", lambda t: None)
    cluster(dp(tmp_path / "data"), tmp_path / "logs", max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 1
    checks = [c for c in fake_slurm.history() if c.startswith("squeue") and "1000_1" in c]
    assert len(checks) >= 3, "The pending job of the last try was not waited for."


def test_a_data_path_without_separator(tmp_path):
    data = tmp_path / "data"
    small_tiles(data, 1)
    run(str(data))
    assert os.path.isfile(dsr_of(tile_fn(data)))
