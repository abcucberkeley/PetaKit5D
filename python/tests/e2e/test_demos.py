"""The demo examples (python/examples/phase1) and the checks on the demo dataset
(tests/demo_data/DemoDatasetTest.m).

The examples' own logic (data paths, registration results, timings, the generic computing
framework's user function) is tested on small synthetic inputs. The DemoDatasetTest checks need
the demo dataset (``demodata``: PETAKIT5D_DEMO_DATA names the PetaKit5D_demo_cell_image_dataset
folder or its parent); the dataset is only read, outputs go to the test's folder. The demo-level
comparison with MATLAB's output trees is ``examples/phase1/compare_demo_trees.py``
(docs/parity/demos.md).
"""

from __future__ import annotations

import importlib
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / "examples" / "phase1"
CELL = "PetaKit5D_demo_cell_image_dataset"
FSN = (
    "Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_"
    "000x_003y_000z_0000t"
)
SKEW_ANGLE = 32.45
XY_PIXEL_SIZE = 0.108
DZ = 0.3


@pytest.fixture
def examples(monkeypatch):
    """the examples folder on sys.path; returns a module importer"""
    monkeypatch.syspath_prepend(str(EXAMPLES))
    return importlib.import_module


def example_modules() -> list[str]:
    return sorted(p.stem for p in EXAMPLES.glob("*.py") if not p.stem.startswith("_"))


# ---------------------------------------------------------------------------------------------
# the examples themselves
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", example_modules())
def test_example_imports_without_running(examples, name):
    """importing an example does no work: everything runs in main()"""
    before = set(sys.modules)
    mod = examples(name)
    assert callable(getattr(mod, "main", None)) or name == "crop_deskew_rotate_demo_function"
    heavy = {m for m in set(sys.modules) - before if m.split(".")[0] in ("numpy", "scipy")}
    assert not heavy, f"{name} imports {sorted(heavy)[:3]} at import time"


def test_paper_demos_are_the_matlab_list(examples):
    """run_all_demos_for_the_paper runs demo_run_all_demos_for_the_paper.m's demos, in order"""
    run_all = examples("run_all_demos_for_the_paper")
    names = [m for _, m in run_all.PAPER_DEMOS]
    assert names == [
        "demo_generic_computing_framework",
        "demo_fast_tiff_zarr_readers_writers",
        "demo_geometric_transformation",
        "demo_rl_deconvolution",
        "demo_zarr_stitching",
        "demo_large_scale_processing",
        "demo_useful_tools",
        "demo_phase_and_2photon_stitching",
        "demo_widefield_and_confocal_deconvolution",
    ]
    assert all(n in example_modules() for n in names)


def test_data_path_takes_the_folder_or_its_parent(examples, monkeypatch, tmp_path):
    common = examples("_common")
    (tmp_path / CELL).mkdir()
    for value in (str(tmp_path), str(tmp_path) + "/", str(tmp_path / CELL)):
        monkeypatch.setenv("PETAKIT5D_DEMO_DATA", value)
        # built as the MATLAB demos build it: [destPath, '/<dataset>/']
        assert common.data_path() == tmp_path.as_posix() + "//" + CELL + "/"
        assert os.path.isdir(common.data_path())
    monkeypatch.setenv("PETAKIT5D_DEMO_DATA", str(tmp_path / "nothing"))
    with pytest.raises(SystemExit, match="is not in"):
        common.data_path()


def test_background_path_is_found_and_written_as_the_demos_write_it(examples, tmp_path):
    common = examples("_common")
    folder = tmp_path / CELL / "FF" / "darks"
    folder.mkdir(parents=True)
    (folder / "AVG_DF400_CamA_10ms.tif").write_bytes(b"")
    data_path = tmp_path.as_posix() + "/" + CELL + "/"
    assert (
        common.background_path(data_path, "CamA") == data_path + "FF/darks/AVG_DF400_CamA_10ms.tif"
    )
    assert common.background_path(data_path, "CamA", "/FF/") == (
        data_path + "/FF/darks/AVG_DF400_CamA_10ms.tif"
    )
    with pytest.raises(SystemExit):
        common.background_path(data_path, "CamB")


@pytest.mark.parametrize(
    "name, calls",
    [
        ("demo_skewed_space_stitching", 1),
        ("demo_zarr_stitching", 1),
        ("demo_phase_and_2photon_stitching", 2),
    ],
)
def test_the_stitching_examples_register_their_tiles(name, calls):
    """every stitching step registers its tiles, as the MATLAB demos do (xcorrShift = true)"""
    import ast

    tree = ast.parse((EXAMPLES / f"{name}.py").read_text(encoding="utf-8"))
    values = [
        kw.value
        for c in ast.walk(tree)
        if isinstance(c, ast.Call) and getattr(c.func, "id", "") == "stitch_dataset"
        for kw in c.keywords
        if kw.arg == "xcorr_shift"
    ]
    assert [isinstance(v, ast.Constant) and v.value is True for v in values] == [True] * calls


def test_steps_are_timed_into_the_timings_file(examples, monkeypatch, tmp_path):
    common = examples("_common")
    log = tmp_path / "t.jsonl"
    monkeypatch.setenv("PETAKIT5D_DEMO_TIMINGS", str(log))
    steps = common.Steps("demo_x")
    with steps("one"):
        pass
    with steps("two"):
        pass
    steps.skip("three", "not here")
    rows = [json.loads(x) for x in log.read_text().splitlines()]
    assert [(r["demo"], r["step"]) for r in rows] == [("demo_x", "one"), ("demo_x", "two")]
    assert all(r["seconds"] >= 0 for r in rows)
    assert b"\r\n" not in log.read_bytes()


def test_in_memory_results_are_saved_only_when_asked(examples, monkeypatch, tmp_path):
    from petakit5d.io import read_tiff

    common = examples("_common")
    im = np.arange(2 * 3 * 4, dtype=np.uint16).reshape(2, 3, 4)
    monkeypatch.delenv("PETAKIT5D_DEMO_SAVE_DIR", raising=False)
    common.save_result("demo_x", "a", im)
    assert not any(tmp_path.iterdir())
    monkeypatch.setenv("PETAKIT5D_DEMO_SAVE_DIR", str(tmp_path))
    common.save_result("demo_x", "a", im)
    np.testing.assert_array_equal(read_tiff(str(tmp_path / "demo_x" / "a.tif")), im)


@pytest.mark.parametrize("nframe", [7, 10, 12, 3])
def test_replicate_frames_is_matlabs_repmat_cropped(examples, nframe):
    """repmat(im, 1, 1, ceil(nframe / nz)) cropped to nframe slices"""
    common = examples("_common")
    im = np.arange(5 * 3 * 4, dtype=np.uint16).reshape(5, 3, 4)
    expected = np.tile(im, (-(-nframe // 5), 1, 1))[:nframe]
    got = common.replicate_frames(im, nframe)
    assert got.flags.c_contiguous
    np.testing.assert_array_equal(got, expected)


def _skewed_tile(nz=24, ny=40, nx=30, seed=3) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(90, 2000, size=(nz, ny, nx)).astype(np.uint16)


def test_crop_deskew_rotate_demo_function_is_the_volume_call(examples, tmp_path):
    """the user function crops, deskews and rotates as deskew_rotate_volume on the crop, and
    publishes its result by one rename"""
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.io import read_tiff, write_tiff

    fn = examples("crop_deskew_rotate_demo_function").crop_deskew_rotate_demo_function
    im = _skewed_tile()
    src, dst = str(tmp_path / "in.tif"), str(tmp_path / "out.tif")
    write_tiff(im, src)
    box = (2, 5, 1, 22, 35, 29)
    fn(src, dst, SKEW_ANGLE, XY_PIXEL_SIZE, DZ, box)
    expected = deskew_rotate_volume(
        im[2:22, 5:35, 1:29].copy(), SKEW_ANGLE, DZ, XY_PIXEL_SIZE, reverse=True, save_16bit=True
    )
    np.testing.assert_array_equal(read_tiff(dst), expected)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["in.tif", "out.tif"]
    with pytest.raises(FileNotFoundError):
        fn(str(tmp_path / "none.tif"), dst, SKEW_ANGLE, XY_PIXEL_SIZE, DZ, box)


def test_generic_framework_example_runs_the_tasks(examples, monkeypatch, tmp_path):
    """the demo's main() on a small dataset: one task per tile, results in Cropped_DSR/, the
    example cluster configuration in the current folder"""
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.io import read_tiff, write_tiff
    from petakit5d.jobs.config import load_config

    data = tmp_path / "Downloads" / CELL
    data.mkdir(parents=True)
    tiles = {f"tile_{i}.tif": _skewed_tile(nz=12, ny=150, nx=20, seed=i) for i in range(2)}
    for name, im in tiles.items():
        write_tiff(im, str(data / name))
    monkeypatch.setenv("PETAKIT5D_DEMO_DATA", str(tmp_path / "Downloads"))
    monkeypatch.chdir(tmp_path)
    examples("demo_generic_computing_framework").main()
    for name, im in tiles.items():
        # MATLAB's box [101, 1, 1, 1700, 512, 500] limited to the image
        expected = deskew_rotate_volume(
            im[:, 100:150, :].copy(), SKEW_ANGLE, DZ, XY_PIXEL_SIZE, reverse=True, save_16bit=True
        )
        np.testing.assert_array_equal(read_tiff(str(data / "Cropped_DSR" / name)), expected)
    cfg = load_config(str(tmp_path / "test_slurm_cpu_config_file.json"))
    assert cfg is not None


def test_compare_demo_trees_rules(examples):
    cmp = examples("compare_demo_trees")
    c = CELL + "/"
    assert cmp.demo_of(c + "Cropped_DSR/a.tif") == "generic_computing_framework"
    assert cmp.demo_of(c + "DSR/a.zarr") == "geometric_transformation"
    assert cmp.demo_of(c + "DSR/a.tif") == "RL_deconvolution"
    assert cmp.demo_of(c + "matlab_stitch/matlab_decon_omw/a.zarr") == "large_scale_processing"
    assert cmp.demo_of(c + "matlab_decon_omw/DSR_for_useful_tools/a.zarr") == "useful_tools"
    assert cmp.demo_of("x/y.tif", "pipeline") == "pipeline"
    assert "legacy" in cmp.not_ported(c + "DS_separated/a.tif")
    assert cmp.not_ported(c + "matlab_stitch/xcorr/v/xcorr_tile_1.mat")
    assert not cmp.not_ported(c + "matlab_stitch/xcorr/v.mat")
    phase = cmp.OTHER + "/Phase/ff_corrected_1/a.tif"
    assert cmp.not_ported(phase) and not cmp.not_ported(phase, frozenset([phase]))


def test_compare_demo_trees_masks_the_roots(examples, tmp_path):
    from petakit5d.io import write_tiff

    cmp = examples("compare_demo_trees")
    im = _skewed_tile(nz=3, ny=5, nx=4)
    roots = {}
    for side, home in (("m", "~"), ("p", str(tmp_path / "p"))):
        d = tmp_path / side / "Downloads" / CELL / "Cropped"
        d.mkdir(parents=True)
        write_tiff(im, str(d / "a.tif"))
        (d / "parameters.json").write_text(
            json.dumps({"dataPaths": [f"{home}/Downloads//{CELL}/"]}), newline="\n"
        )
        (d / "empty_flag").mkdir()
        roots[side] = str(tmp_path / side / "Downloads")
    (tmp_path / "inputs" / CELL).mkdir(parents=True)
    args = [roots["m"], roots["p"], "--inputs", str(tmp_path / "inputs")]
    args += ["--python-home", str(tmp_path / "p"), "--workers", "1"]
    assert cmp.main(args) == 0
    write_tiff(im + 1, str(Path(roots["p"]) / CELL / "Cropped" / "a.tif"))
    assert cmp.main(args) == 1


# ---------------------------------------------------------------------------------------------
# DemoDatasetTest (tests/demo_data/DemoDatasetTest.m)
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def dataset() -> Path:
    """the demo cell dataset (PETAKIT5D_DEMO_DATA: the folder or its parent), or a skip"""
    value = os.environ.get("PETAKIT5D_DEMO_DATA", "")
    path = Path(value)
    if (path / CELL).is_dir():
        path = path / CELL
    if not value or not (path / (FSN + ".tif")).is_file():
        pytest.skip(f"demodata: {value!r} does not contain the demo dataset")
    return path


def copy_tile(dataset: Path, tmp_path: Path) -> str:
    shutil.copy(dataset / (FSN + ".tif"), tmp_path)
    return tmp_path.as_posix() + "/"


def gradient_energy(im: np.ndarray) -> float:
    """mean squared differences along y and x of the middle third of the z slices"""
    from petakit5d.core.dtypes import mround

    nz = im.shape[0]
    im = im[int(mround(nz / 3)) - 1 : int(mround(2 * nz / 3))].astype(np.float64)
    return float(np.mean(np.diff(im, axis=-2) ** 2) + np.mean(np.diff(im, axis=-1) ** 2))


@pytest.mark.demodata
@pytest.mark.slow
def test_cpp_tiff_reader_matches_libtiff(dataset):
    """cppTiffReaderMatchesLibtiff: cpp-tiff's reader against an independent reader (tifffile,
    which decodes LZW with imagecodecs: skipped without it), and a range read"""
    import importlib.util

    from petakit5d.core.dtypes import mround
    from petakit5d.io import get_image_size, read_tiff

    fn = str(dataset / (FSN + ".tif"))
    im = read_tiff(fn)
    sz = get_image_size(fn)
    assert im.shape == sz
    assert sz[1:] == (1800, 512)
    z0 = int(mround(sz[0] / 2)) - 1
    np.testing.assert_array_equal(read_tiff(fn, range=(z0, z0 + 10)), im[z0 : z0 + 10])
    if importlib.util.find_spec("imagecodecs") is None:
        pytest.skip("the independent reader needs tifffile and imagecodecs")
    import tifffile

    np.testing.assert_array_equal(im, tifffile.imread(fn))


@pytest.mark.demodata
@pytest.mark.slow
def test_writers_round_trip_real_data(dataset, tmp_path):
    """writersRoundTripRealData"""
    from petakit5d.io import read_tiff, read_zarr, write_tiff, write_zarr

    im = read_tiff(str(dataset / (FSN + ".tif")))
    write_tiff(im, str(tmp_path / "a.tif"))
    np.testing.assert_array_equal(read_tiff(str(tmp_path / "a.tif")), im)
    write_zarr(im, str(tmp_path / "a.zarr"), block_size=(256, 256, 256))
    np.testing.assert_array_equal(read_zarr(str(tmp_path / "a.zarr")), im)


@pytest.mark.demodata
@pytest.mark.slow
def test_combined_deskew_rotate_matches_separate(dataset, tmp_path):
    """combinedDeskewRotateMatchesSeparate: the comparison of demo_geometric_transformation"""
    from petakit5d.deskew_rotate import deskew_rotate_file
    from petakit5d.io import read_tiff

    out_path = copy_tile(dataset, tmp_path)
    fn = out_path + FSN + ".tif"
    deskew_rotate_file(
        fn,
        XY_PIXEL_SIZE,
        DZ,
        skew_angle=SKEW_ANGLE,
        reverse=True,
        dsr_combined=False,
        rotate=True,
        ds_dir_name="DS_separated/",
        dsr_dir_name="DSR_separated/",
    )
    deskew_rotate_file(fn, XY_PIXEL_SIZE, DZ, skew_angle=SKEW_ANGLE, reverse=True)
    separate = read_tiff(out_path + "DSR_separated/" + FSN + ".tif")
    combined = read_tiff(out_path + "DSR/" + FSN + ".tif")
    assert combined.shape == separate.shape
    # the separate path interpolates twice; on this data the difference is about 2%
    a, e = combined.astype(np.float64), separate.astype(np.float64)
    assert np.linalg.norm(a - e) / np.linalg.norm(e) <= 0.04


@pytest.mark.demodata
@pytest.mark.slow
def test_omw_deconvolution_on_cpu(dataset, tmp_path):
    """omwDeconvolutionOnCPU: OMW decon as in demo_RL_deconvolution, on the CPU"""
    from petakit5d.core.threads import available_cpus
    from petakit5d.decon import decon_dataset
    from petakit5d.io import read_tiff

    out_path = copy_tile(dataset, tmp_path)
    decon_dataset(
        [out_path],
        result_dir_name="matlab_decon_omw",
        xy_pixel_size=XY_PIXEL_SIZE,
        dz=DZ,
        channel_patterns=["CamB_ch0"],
        psf_fullpaths=[str(dataset / "PSF" / "488_2_c.tif")],
        dz_psf=0.3,
        rl_method="omw",
        wiener_alpha=0.005,
        otf_cum_thresh=0.87,
        skewed=True,
        background=100,
        decon_iter=2,
        save_16bit=True,
        parse_cluster=False,
        gpu_job=False,
        cpus_per_task=available_cpus(),
    )
    raw = read_tiff(out_path + FSN + ".tif")
    decon = read_tiff(out_path + "matlab_decon_omw/" + FSN + ".tif")
    assert decon.dtype == np.uint16
    assert decon.shape == raw.shape
    # deconvolution sharpens: more energy in high spatial frequencies
    assert gradient_energy(decon) > 2 * gradient_energy(np.maximum(raw.astype(np.float32) - 100, 0))
