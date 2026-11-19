"""The examples of the demos' Phase 2 steps (python/examples/phase2).

The examples' own logic is tested on small synthetic inputs: imports do no work, the run-all
runs each phase 1 example with its Phase 2 steps after it, the data tools steps call the tools
with the MATLAB demo's parameters (in Python's form), and the step-by-step comparison of output
trees applies the parity rows' tolerances and reports the numbers. The comparison of the demo
runs with MATLAB's is ``examples/phase2/compare_phase2_trees.py`` (docs/parity/phase2.md).
"""

from __future__ import annotations

import importlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

PHASE2 = Path(__file__).resolve().parents[2] / "examples" / "phase2"
CELL = "PetaKit5D_demo_cell_image_dataset"
OTHER = "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets"
TILE_B = (
    "Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_"
    "000x_003y_000z_0000t"
)


@pytest.fixture
def phase2(monkeypatch):
    """loads an example of examples/phase2 by path (as ``phase2_<name>``), with
    examples/phase1 first on the import path, as the examples have it"""
    monkeypatch.syspath_prepend(str(PHASE2))
    loader = importlib.import_module("_phase1")
    monkeypatch.syspath_prepend(loader.PHASE1)
    return loader.load


def example_modules() -> list[str]:
    return sorted(p.stem for p in PHASE2.glob("*.py") if not p.stem.startswith("_"))


def test_the_examples_are_these():
    assert example_modules() == [
        "compare_phase2_trees",
        "demo_useful_tools",
        "run_all_demos_for_the_paper",
    ]


@pytest.mark.parametrize("name", example_modules())
def test_example_imports_without_running(phase2, name):
    """importing an example does no work: everything runs in main()"""
    before = set(sys.modules)
    mod = phase2(name)
    assert callable(getattr(mod, "main", None))
    heavy = {m for m in set(sys.modules) - before if m.split(".")[0] in ("numpy", "scipy")}
    assert not heavy, f"{name} imports {sorted(heavy)[:3]} at import time"


def test_phase1_is_put_first_on_the_path(phase2, monkeypatch):
    import _phase1

    monkeypatch.setattr(sys, "path", [str(PHASE2), "elsewhere", _phase1.PHASE1])
    importlib.reload(_phase1)
    assert sys.path == [_phase1.PHASE1, str(PHASE2), "elsewhere"]
    assert Path(_phase1.PHASE1) == PHASE2.parent / "phase1"


def test_phase1_namesakes_stay_phase1s(phase2):
    """``import demo_useful_tools`` is phase 1's example; this folder's is loaded as
    phase2_demo_useful_tools"""
    p2 = phase2("demo_useful_tools")
    p1 = importlib.import_module("demo_useful_tools")
    assert Path(p1.__file__).parent.name == "phase1"
    assert Path(p2.__file__).parent.name == "phase2"
    assert p2 is sys.modules["phase2_demo_useful_tools"] is phase2("demo_useful_tools")


# ---------------------------------------------------------------------------------------------
# the run-all
# ---------------------------------------------------------------------------------------------


def test_run_all_is_the_papers_list_with_the_phase2_steps(phase2):
    run_all = phase2("run_all_demos_for_the_paper")
    names = [m for _, m in run_all.paper_demos()]
    assert names[6] == "demo_useful_tools" and len(names) == 9
    assert run_all.PHASE2_STEPS == {"demo_useful_tools": "demo_useful_tools"}
    assert all((PHASE2 / f"{m}.py").is_file() for m in run_all.PHASE2_STEPS.values())


def test_run_all_runs_the_phase2_steps_after_their_phase1_example(phase2, monkeypatch):
    run_all = phase2("run_all_demos_for_the_paper")
    calls = []
    for _, m in run_all.paper_demos():
        fake = type(sys)(m)
        fake.main = lambda m=m: calls.append(("phase1", m))
        monkeypatch.setitem(sys.modules, m, fake)
    monkeypatch.setattr(phase2("demo_useful_tools"), "main", lambda: calls.append(("phase2", "u")))
    run_all.main(["demo_phase_and_2photon_stitching", "demo_useful_tools"])
    assert calls == [
        ("phase1", "demo_useful_tools"),
        ("phase2", "u"),
        ("phase1", "demo_phase_and_2photon_stitching"),
    ]
    with pytest.raises(SystemExit, match="unknown demos"):
        run_all.main(["demo_nothing"])


# ---------------------------------------------------------------------------------------------
# the data tools demo's Phase 2 steps
# ---------------------------------------------------------------------------------------------

TOOLS = (
    "fft_spectrum_dataset",
    "psf_analysis_dataset",
    "fsc_analysis_dataset",
    "chromatic_shift_correction_dataset",
)


@pytest.fixture
def demo_data(monkeypatch, tmp_path):
    """a dataset folder with the phase 1 results these steps need (empty folders)"""
    cell = tmp_path / "Downloads" / CELL
    (cell / "matlab_decon_omw" / "DS_for_useful_tools").mkdir(parents=True)
    (cell / "matlab_decon_omw_debug" / f"{TILE_B}_debug").mkdir(parents=True)
    monkeypatch.setenv("PETAKIT5D_DEMO_DATA", str(tmp_path / "Downloads"))
    monkeypatch.delenv("PETAKIT5D_DEMO_TIMINGS", raising=False)
    return tmp_path.as_posix() + "/Downloads//" + CELL + "/"


@pytest.fixture
def spied_tools(monkeypatch):
    """the four tools replaced by recorders: [(name, args, kwargs)]"""
    import petakit5d.tools as tools

    calls = []
    for name in TOOLS:
        getattr(tools, name)
        monkeypatch.setattr(
            tools, name, lambda *a, _n=name, **k: calls.append((_n, a, k)), raising=True
        )
    return calls


def _sind_matlab(x: float) -> float:
    from petakit5d.matlabcompat.elementary import sind

    return float(sind(x))


def expected_calls(dp: str) -> list:
    dz = 0.3 * _sind_matlab(32.45)
    jobs = {"parse_cluster": False, "master_compute": True, "cpus_per_task": 24}
    jobs |= {"mcc_mode": False, "config_file": ""}
    return [
        (
            "fft_spectrum_dataset",
            ([dp + "matlab_decon_omw/DS_for_useful_tools/"],),
            {
                "result_dir_name": "FFT",
                "xy_pixel_size": 0.108,
                "dz": dz,
                "zarr_file": True,
                "out_size": (1001, 1001, 1001),
                "out_pixel_size": 0.1 * 0.488 / 1.33,
                "save_3d_stack": False,
                **jobs,
            },
        ),
        (
            "psf_analysis_dataset",
            ([dp + "PSF/"],),
            {
                "xy_pixel_size": 0.108,
                "dz": 0.3,
                "skew_angle": 32.45,
                "channel_patterns": ["488"],
                "channels": [488],
                "deskew": True,
                "objective_scan": False,
                "source_str": "Demo PSF",
                "rw_fn": [dp + "PSF/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif"],
                "visible": True,
                **jobs,
            },
        ),
        (
            "fsc_analysis_dataset",
            ([dp + f"matlab_decon_omw_debug/{TILE_B}_debug/"],),
            {
                "result_dir_name": "FSCs",
                "xy_pixel_size": 0.108,
                "dz": dz,
                "dr": 5,
                "dtheta": math.pi / 6,
                "res_thresh_method": "one-bit",
                "res_thresh": 0.2,
                "half_size": (201, 201, 201),
                "input_bbox": None,
                "channel_patterns": ["Iter"],
                "channels": [488],
                "suffix": "decon",
                "iter_interval": 2,
                **jobs,
            },
        ),
        (
            "chromatic_shift_correction_dataset",
            ([dp],),
            {
                # MATLAB's [0, 0, 0; 1, -5, -2] ([y x z] rows) as (z, y, x) rows
                "chromatic_offset": [(0, 0, 0), (-2, 1, -5)],
                "result_dir_name": "Chromatic_Shift_Corrected",
                "mode": "valid",
                "pad_value": 0,
                "new_origin": None,
                "channel_patterns": ["CamA", "CamB"],
                "psf_fullpaths": [],
                "max_offset": (20, 20, 20),
                "crop_length": (0, 0, 0),
                "zarr_file": False,
                "large_file": False,
                "save_zarr": False,
                "batch_size": (1024, 1024, 1024),
                "block_size": (256, 256, 256),
                "parse_cluster": True,
                "master_compute": True,
                "job_log_dir": "../job_logs",
                "cpus_per_task": 3,
                "config_file": "",
                "mcc_mode": False,
                "uuid": "",
                "debug": False,
            },
        ),
    ]


def test_useful_tools_steps_have_the_demos_parameters(phase2, demo_data, spied_tools, monkeypatch):
    """the four steps, in the demo's order, with MATLAB's values in Python's form"""
    import _common

    monkeypatch.setattr(_common, "figures_available", lambda: True)
    phase2("demo_useful_tools").main()
    assert spied_tools == expected_calls(demo_data)


def test_useful_tools_leaves_out_the_psf_analysis_without_matplotlib(
    phase2, demo_data, spied_tools, monkeypatch, capsys
):
    import _common

    monkeypatch.setattr(_common, "figures_available", lambda: False)
    phase2("demo_useful_tools").main()
    assert [c[0] for c in spied_tools] == [TOOLS[0], TOOLS[2], TOOLS[3]]
    out = capsys.readouterr().out
    assert "PSF analysis: skipped (matplotlib is not installed" in out
    assert "Imaris converter: skipped" in out


def test_useful_tools_needs_the_phase1_results(phase2, demo_data, spied_tools):
    import shutil

    shutil.rmtree(demo_data + "matlab_decon_omw_debug")
    with pytest.raises(SystemExit, match="run examples/phase1/demo_useful_tools.py first"):
        phase2("demo_useful_tools").main()
    assert spied_tools == []


def test_deskewed_z_step_is_matlabs(phase2):
    """0.3 * sind(32.45) as MATLAB computes it (the value in MATLAB's FFT parameters.mat)"""
    dz = phase2("demo_useful_tools").dz_deskewed()
    assert dz == 0.3 * _sind_matlab(32.45)
    assert dz.hex() == "0x1.49aa20621a5bfp-3"


def test_bioformats_conversion_is_matlabs_section_after_the_return(phase2, monkeypatch):
    import petakit5d.io.bioformats as bf

    calls = []
    monkeypatch.setattr(
        bf, "bioformats_to_tiff_or_zarr_dataset", lambda *a, **k: calls.append((a, k))
    )
    phase2("demo_useful_tools").bioformats_conversion()
    assert calls == [
        (
            (["/data/microscopy/sample_dataset_1/", "/data/microscopy/sample_dataset_2/"],),
            {
                "result_dir_name": "tiffs",
                "channel_patterns": [".nd2", ".czi"],
                "n_channels": 2,
                "data_formats": [".nd2", ".czi"],
                "save_zarr": False,
                "block_size": (256, 256, 256),
                "over_write": False,
                "uuid": "",
            },
        )
    ]


# ---------------------------------------------------------------------------------------------
# the comparison of the Phase 2 steps' outputs
# ---------------------------------------------------------------------------------------------


def test_compare_rules(phase2):
    cmp = phase2("compare_phase2_trees")
    c, o = CELL + "/", OTHER + "/"
    fft = c + "matlab_decon_omw/DS_for_useful_tools/FFT/"
    assert cmp.step_of(fft + "MIPs/a_MIP_z.tif") == ("useful_tools", "FFT analysis")
    assert cmp.step_of(c + "PSF/DS/PSFAnalysis/RW_info.mat")[1] == "PSF analysis"
    assert cmp.step_of(c + "PSF/RW_PSFs/x.tif") is None
    assert cmp.step_of(c + f"matlab_decon_omw_debug/{TILE_B}_debug/FSCs/Iter_0002.mat")[1] == (
        "FSC analysis"
    )
    assert cmp.step_of(c + f"matlab_decon_omw_debug/{TILE_B}_debug/Iter_0002.tif") is None
    assert cmp.step_of(c + "Chromatic_Shift_Corrected/a.tif")[1] == "chromatic shift correction"
    assert cmp.step_of(o + "Phase/ff_corrected_1/BaSiC_flatfield/a_flatfield.mat") == (
        "phase_and_2photon_stitching",
        "phase flat field",
    )
    assert cmp.step_of(o + "Phase/ff_corrected/matlab_stitch/a.zarr") is None
    assert cmp.step_of(c + "matlab_decon_omw/psfgen/x_figure.png")[1] == (
        "figures of the phase 1 steps"
    )
    assert cmp.step_of(c + "Cropped/a.tif") is None
    assert "Imaris" in cmp.pending(c + "matlab_decon_omw/DSR_for_useful_tools/imaris/a.ims")
    assert "F6" in cmp.pending(c + "x/figures/a.fig")
    assert not cmp.pending(fft + "parameters.mat")


def test_name_of_finds_the_output_of_a_problem(phase2):
    cmp = phase2("compare_phase2_trees")
    names = {"a/b.mat", "a/b.mat2", "a/c.tif"}
    assert cmp.name_of("files/a/b.mat", names) == "a/b.mat"
    assert cmp.name_of("files/a/b.mat.res_mu.points", names) == "a/b.mat"
    assert cmp.name_of("files/a/b.mat2{1}", names) == "a/b.mat2"
    assert cmp.name_of("files/a/c.tif (voxels)", names) == "a/c.tif"
    assert cmp.name_of("files/a/d.tif", names) is None


def test_array_stats(phase2):
    cmp = phase2("compare_phase2_trees")
    r = np.array([1.0, 2.0, 4.0, np.nan], np.float32)
    a = r.copy()
    assert cmp.array_stats(r, a)["differ"] == 0
    a[2] = np.nextafter(np.float32(4), np.float32(5))
    a[0] = 1.5
    s = cmp.array_stats(r, a)
    assert (s["differ"], s["max_abs"], s["max_rel"]) == (2, 0.5, 0.5)
    assert s["max_rel_to_max"] == 0.5 / 4 and s["max_ulp"] == 2**22
    assert s["rel_l2"] == pytest.approx(np.hypot(0.5, 2**-21) / np.sqrt(21))
    u = cmp.array_stats(np.array([3, 7], np.uint16), np.array([4, 7], np.uint16))
    assert (u["differ"], u["max_abs"], "max_ulp" in u) == (1, 1.0, False)
    assert "size" in cmp.array_stats(np.zeros(3), np.zeros(4))
    # percentiles of the relative differences of all values (0 where equal)
    r = np.ones(200)
    a = r.copy()
    a[:10] *= 1.5
    a[10] = 1.25
    s = cmp.array_stats(r, a)
    assert (s["rel_p50"], s["rel_p99"], s["max_rel"]) == (0.0, 0.5, 0.5)
    text = cmp.summary_text({"x.tif": [("voxels", s)], "y.mat": [("v", u)]})
    assert "12 of 202 values of 2 arrays differ" in text
    assert "largest relative 0.5" in text and "99th percentile relative 0.5" in text


def _res(points, res):
    return {
        "comp_angles": np.array([[0.0], [math.pi / 2]]),
        "points": np.asarray(points, np.float64).reshape(-1, 1),
        "thetas": np.array([[0.0], [math.pi / 6]]),
        "res": np.asarray(res, np.float64).reshape(-1, 1),
    }


def _fsc(v):
    return {
        "r": np.arange(3.0).reshape(-1, 1),
        "thetas": np.array([[0.0], [math.pi / 6]]),
        "fsc": np.asarray(v, np.float64),
        "npoints": np.ones((3, 2)),
    }


def _write_trees(root: Path, side: str, *, fft_mip, tile, flat, points, fsc, png=True, fig=False):
    from petakit5d.io import write_mat, write_tiff

    d = root / side / "Downloads"
    cell, other = d / CELL, d / OTHER
    fft = cell / "matlab_decon_omw" / "DS_for_useful_tools" / "FFT"
    (fft / "MIPs").mkdir(parents=True)
    write_tiff(fft_mip, str(fft / "MIPs" / "a_MIP_z.tif"))
    (fft / "a.zarr").write_bytes(b"")
    (cell / "Chromatic_Shift_Corrected").mkdir()
    write_tiff(tile, str(cell / "Chromatic_Shift_Corrected" / "a.tif"))
    fsc_dir = cell / "matlab_decon_omw_debug" / f"{TILE_B}_debug"
    (fsc_dir / "FSCs").mkdir(parents=True)
    write_mat(
        str(fsc_dir / "FSCs" / "Iter_0002.mat"),
        {"fsc_mu": _fsc(fsc), "res_mu": _res(points, [0.25, 0.3])},
    )
    (fsc_dir / "figures").mkdir()
    if png:
        (fsc_dir / "figures" / "fsc_res_vs_thetas_Iter_0002_decon.png").write_bytes(side.encode())
    if fig:
        (fsc_dir / "figures" / "fsc_x.fig").write_bytes(b"MATLAB 5.0 MAT-file")
        imaris = cell / "matlab_decon_omw" / "DSR_for_useful_tools" / "imaris"
        imaris.mkdir(parents=True)
        (imaris / "a.ims").write_bytes(b"hdf5")
    ff = other / "Phase" / "ff_corrected_1"
    (ff / "BaSiC_flatfield").mkdir(parents=True)
    write_tiff(flat[0], str(ff / "t.tif"))
    write_mat(
        str(ff / "BaSiC_flatfield" / "t_flatfield.mat"),
        {"minval": np.float32(flat[2]), "flatfield": flat[1]},
    )
    (cell / "Cropped").mkdir()
    write_tiff(tile + (side == "p"), str(cell / "Cropped" / "phase1_output.tif"))
    return str(d)


@pytest.fixture
def trees(tmp_path):
    """MATLAB's and Python's output trees of a few Phase 2 steps, and the folder of inputs"""
    rng = np.random.default_rng(7)
    mip = rng.random((1, 9, 11), dtype=np.float32) + np.float32(0.1)
    tile = rng.integers(100, 4000, size=(3, 5, 7)).astype(np.uint16)
    ff = rng.random((5, 7)) + 0.5
    flat = (tile, ff, 1234.5)
    fsc = rng.random((3, 2))
    (tmp_path / "inputs" / CELL).mkdir(parents=True)
    m = _write_trees(
        tmp_path, "m", fft_mip=mip, tile=tile, flat=flat, points=[0.3, 0.4], fsc=fsc, fig=True
    )
    # an input of the dataset, in the inputs folder and in both trees: not an output
    for root in (tmp_path / "inputs", tmp_path / "m" / "Downloads", tmp_path / "p" / "Downloads"):
        (root / CELL).mkdir(parents=True, exist_ok=True)
        (root / CELL / "raw.csv").write_bytes(b"x,y\n")
    return tmp_path, (mip, tile, flat, fsc), m


def _run(cmp, m, p, tmp_path, backend="exact", label="t", *more):
    out = tmp_path / f"out_{label}"
    args = [m, p, "--inputs", str(tmp_path / "inputs"), "--backend", backend, *more]
    args += ["--out", str(out), "--label", label, "--workers", "1"]
    rc = cmp.main(args)
    return (
        rc,
        json.loads((out / f"compare_{label}.json").read_text()),
        (out / f"compare_{label}.txt").read_text(),
    )


def test_compare_equal_trees(phase2, trees):
    cmp = phase2("compare_phase2_trees")
    tmp, (mip, tile, flat, fsc), m = trees
    p = _write_trees(tmp, "p", fft_mip=mip, tile=tile, flat=flat, points=[0.3, 0.4], fsc=fsc)
    rc, res, text = _run(cmp, m, p, tmp)
    assert rc == 0, text
    fft = res["steps"]["useful_tools: FFT analysis"]
    assert (fft["matlab"], fft["identical"], fft["differing"]) == (2, 2, [])
    fsc_step = res["steps"]["useful_tools: FSC analysis"]
    assert (fsc_step["identical"], fsc_step["figures"]) == (1, 1)
    # the Imaris output and MATLAB's .fig are listed, not compared; the phase 1 output is left
    # to the phase 1 tool
    assert sorted(len(v) for v in res["not_compared"].values()) == [1, 1]
    assert "outputs of phase 1 steps (examples/phase1/compare_demo_trees.py): 1" in text


def test_compare_finds_and_measures_differences(phase2, trees):
    cmp = phase2("compare_phase2_trees")
    tmp, (mip, tile, flat, fsc), m = trees
    tile2 = tile.copy()
    tile2[0, 0, 0] += 3
    p = _write_trees(
        tmp, "p", fft_mip=mip, tile=tile2, flat=flat, points=[0.3, 0.4], fsc=fsc, png=False
    )
    rc, res, text = _run(cmp, m, p, tmp)
    assert rc == 1
    chrom = res["steps"]["useful_tools: chromatic shift correction"]
    assert chrom["differing"] == [f"{CELL}/Chromatic_Shift_Corrected/a.tif"]
    ((key, stats),) = chrom["numbers"][chrom["differing"][0]]
    assert (key, stats["differ"], stats["max_abs"]) == ("voxels", 1, 3.0)
    assert "1 of 105 uint16 values differ, largest by 3" in text
    fsc_step = res["steps"]["useful_tools: FSC analysis"]
    assert [p["what"] for p in fsc_step["problems"]] == ["missing"]


def test_compare_applies_the_backends_tolerances(phase2, trees):
    """P9 (the spectra on the default backend), N3 (FSC radii within a step), N8 (BaSiCPy)"""
    cmp = phase2("compare_phase2_trees")
    tmp, (mip, tile, flat, fsc), m = trees
    mip2 = mip.copy()
    mip2[0, 3, 4] = np.nextafter(mip2[0, 3, 4], np.float32(2))
    ff2 = flat[1] * (1 + 0.01 * np.cos(np.arange(flat[1].size).reshape(flat[1].shape)))
    tiles2 = tile + tile // 200  # about 0.5% more
    flat2 = (tiles2, ff2, flat[2] + 0.004)
    fsc2 = fsc * (1 + 1e-15)
    p = _write_trees(tmp, "p", fft_mip=mip2, tile=tile, flat=flat2, points=[0.3001, 0.4], fsc=fsc2)
    rc, res, text = _run(cmp, m, p, tmp, backend="default", label="d")
    assert rc == 0, text
    assert res["steps"]["useful_tools: FFT analysis"]["differing"] != []
    phase = res["steps"]["phase_and_2photon_stitching: phase flat field"]
    assert len(phase["differing"]) == 2 and phase["problems"] == []
    # the exact backend allows no difference in the spectra or the FSC curves
    rc, res, text = _run(cmp, m, p, tmp, backend="exact", label="e")
    assert rc == 1
    assert res["steps"]["useful_tools: FFT analysis"]["problems"]
    assert res["steps"]["useful_tools: FSC analysis"]["problems"]
    assert res["steps"]["phase_and_2photon_stitching: phase flat field"]["problems"] == []


def test_compare_fsc_radii_beyond_a_step_differ(phase2, trees):
    cmp = phase2("compare_phase2_trees")
    tmp, (mip, tile, flat, fsc), m = trees
    p = _write_trees(tmp, "p", fft_mip=mip, tile=tile, flat=flat, points=[0.3002, 0.4], fsc=fsc)
    rc, res, text = _run(cmp, m, p, tmp, backend="default")
    assert rc == 1
    probs = res["steps"]["useful_tools: FSC analysis"]["problems"]
    assert len(probs) == 1 and "up to 2 steps of 1e-4" in " ".join(probs[0]["lines"])


@pytest.mark.parametrize("which", ["tile_rel", "flatfield_rel_l2"])
def test_compare_flat_field_beyond_n8_differs(phase2, trees, which):
    """N8: a tile value 20% off (relative 0.06), or the flat field 5% off everywhere (relative L2
    0.01), each beyond one limit only"""
    cmp = phase2("compare_phase2_trees")
    tmp, (mip, tile, flat, fsc), m = trees
    if which == "tile_rel":
        tile2 = flat[0].copy()
        k = np.unravel_index(np.argmin(tile2), tile2.shape)  # its smallest value: little L2
        tile2[k] = tile2[k] // 5 * 6
        flat2, name = (tile2, flat[1], flat[2]), "t.tif"
    else:
        flat2, name = (flat[0], flat[1] * 1.05, flat[2]), "t_flatfield.mat"
    p = _write_trees(tmp, "p", fft_mip=mip, tile=tile, flat=flat2, points=[0.3, 0.4], fsc=fsc)
    rc, res, text = _run(cmp, m, p, tmp, backend="default")
    assert rc == 1
    probs = res["steps"]["phase_and_2photon_stitching: phase flat field"]["problems"]
    assert [p["path"].split("/")[-1].split(" ")[0] for p in probs] == [name]


def test_compare_spectra_on_the_default_backend(phase2, trees):
    """P14 (P2's relative L2) end to end; P9 for the step on MATLAB's own decon results"""
    cmp = phase2("compare_phase2_trees")
    tmp, (mip, tile, flat, fsc), m = trees
    mip2 = mip.copy()
    mip2[0, 3, 4] *= np.float32(1.0001)
    p = _write_trees(tmp, "p", fft_mip=mip2, tile=tile, flat=flat, points=[0.3, 0.4], fsc=fsc)
    rc, res, text = _run(cmp, m, p, tmp, "default", "p14")
    assert rc == 0, text
    rc, res, text = _run(cmp, m, p, tmp, "default", "p9", "--matlab-inputs")
    assert rc == 1
    probs = res["steps"]["useful_tools: FFT analysis"]["problems"]
    assert "largest difference of |F|/max" in " ".join(probs[0]["lines"])
    assert "the steps on MATLAB's decon results" in text
    mip3 = mip * np.float32(1.0005)  # relative L2 5e-4
    p = _write_trees(tmp, "q", fft_mip=mip3, tile=tile, flat=flat, points=[0.3, 0.4], fsc=fsc)
    rc, res, text = _run(cmp, m, p, tmp, "default", "p14b")
    assert rc == 1
    assert "P14" in " ".join(res["steps"]["useful_tools: FFT analysis"]["problems"][0]["lines"])


def test_compare_fsc_curves_on_the_default_backend(phase2, trees):
    cmp = phase2("compare_phase2_trees")
    tmp, (mip, tile, flat, fsc), m = trees
    fsc2 = fsc * (1 + 1e-6)
    p = _write_trees(tmp, "p", fft_mip=mip, tile=tile, flat=flat, points=[0.3, 0.4], fsc=fsc2)
    assert _run(cmp, m, p, tmp, "default", "c14")[0] == 0
    rc, res, text = _run(cmp, m, p, tmp, "default", "c11", "--matlab-inputs")
    assert rc == 1
    assert "fsc_mu" in " ".join(res["steps"]["useful_tools: FSC analysis"]["problems"][0]["lines"])
