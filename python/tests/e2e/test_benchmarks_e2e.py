"""The end-to-end benchmark files (python/benchmarks/e2e, docs/benchmarks_e2e.md): the task-cost
spec, the task probe and its call folders for the job entry point, and the workflow driver in
local mode on small synthetic tiles (the cluster runs are the benchmarks themselves)."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

PY = Path(__file__).resolve().parents[2]
E2E = PY / "benchmarks" / "e2e"
EXAMPLES = PY / "examples" / "phase1"


def _module(name: str):
    spec = importlib.util.spec_from_file_location(f"_e2e_{name}", E2E / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _tile(nz: int, ny: int, nx: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(0, 4000, size=(nz, ny, nx), dtype=np.uint16)


def test_task_costs_spec_loads_and_names_existing_scripts():
    from petakit5d.testing.bench import load_spec

    spec = load_spec(str(E2E / "task_costs.json"))
    ids = [v["id"] for v in spec["variants"]]
    assert ids == [
        "ml_dsr",
        "mcc_dsr",
        "py_dsr",
        "pyexec_dsr",
        "ml_tiny",
        "mcc_tiny",
        "py_tiny",
        "pyexec_tiny",
    ]
    for v in spec["variants"]:
        if v.get("script"):
            assert (PY / v["script"]).is_file()
    # the compiled variants run only when a build is named
    assert {v["id"] for v in spec["variants"] if v.get("when")} == {"mcc_dsr", "mcc_tiny"}
    assert spec["gate"]["overhead_max"] == 1.0


def _env() -> dict[str, str]:
    env = dict(os.environ)
    paths = [p for p in sys.path if p]
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


def test_probe_call_folder_runs_the_tiny_crop_through_the_entry_point(tmp_path):
    """``task_probe.py call`` writes a call folder that ``python -m petakit5d.jobs.exec`` runs:
    the crop of the 32^3 box"""
    from petakit5d.io import read_tiff, write_tiff

    probe = _module("task_probe")
    vol = _tile(40, 48, 36, 1)
    src, dst = tmp_path / "x.tif", tmp_path / "out" / "tiny.tif"
    write_tiff(vol, str(src))
    probe.main(["call", str(tmp_path / "call"), "tiny", str(src), str(dst)])
    r = subprocess.run(
        [sys.executable, "-m", "petakit5d.jobs.exec", str(tmp_path / "call"), "1"],
        capture_output=True,
        text=True,
        env=_env(),
        cwd=str(tmp_path),
        timeout=300,
    )
    assert r.returncode == 0, r.stderr
    np.testing.assert_array_equal(read_tiff(str(dst)), vol[:32, :32, :32])


def test_probe_dsr_call_names_the_file_task(tmp_path):
    from petakit5d.jobs.tasks import read_call

    probe = _module("task_probe")
    probe.main(["call", str(tmp_path / "call"), "dsr", str(tmp_path / "tile.tif")])
    call = read_call(str(tmp_path / "call"))
    assert call.n == 1
    rec = call.records[0]
    assert rec["func"] == "petakit5d.deskew_rotate.file:deskew_rotate_file"
    assert rec["args"][1:] == [0.108, 0.3]
    assert rec["kwargs"]["reverse"] is True and rec["kwargs"]["rotate"] is True


def test_workflow_generic_in_local_mode(tmp_path):
    """the driver's generic workflow: the demo's task on every tile, as the demo runs it"""
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.io import read_tiff, write_tiff

    wf = _module("e2e_cluster_workflow")
    data = tmp_path / "data"
    data.mkdir()
    tiles = {f"tile_{i}.tif": _tile(12, 150, 20, i) for i in range(2)}
    for name, im in tiles.items():
        write_tiff(im, str(data / name))
    timing = tmp_path / "timing.txt"
    wf.main(["generic", str(data), "", str(timing), str(EXAMPLES)])
    for name, im in tiles.items():
        expected = deskew_rotate_volume(
            im[:, 100:150, :].copy(), 32.45, 0.3, 0.108, reverse=True, save_16bit=True
        )
        np.testing.assert_array_equal(read_tiff(str(data / "Cropped_DSR" / name)), expected)
    line = timing.read_text().split()
    assert line[0] == "generic" and float(line[1]) > 0


def test_workflow_dsr_in_local_mode(tmp_path):
    """the driver's dsr workflow: a zarr result per tile of the two channels"""
    from petakit5d.io import write_tiff

    wf = _module("e2e_cluster_workflow")
    data = tmp_path / "data"
    data.mkdir()
    names = ["s_CamA_ch0_t0.tif", "s_CamB_ch0_t0.tif", "other.tif"]
    for i, name in enumerate(names):
        write_tiff(_tile(12, 40, 20, i), str(data / name))
    wf.main(["dsr", str(data), "", str(tmp_path / "timing.txt")])
    done = sorted(p.name for p in (data / "DSR").glob("*.zarr"))
    assert done == ["s_CamA_ch0_t0.zarr", "s_CamB_ch0_t0.zarr"]


def test_workflow_unknown_name_is_refused(tmp_path):
    wf = _module("e2e_cluster_workflow")
    with pytest.raises(SystemExit):
        wf.main(["nothing", str(tmp_path), "", str(tmp_path / "t.txt")])
