"""The PSF sigma estimate (GU_estimateSigma3D) and the camera threshold
(XR_get_camera_standard_deviation_threshold, get_camera_type_from_setting_file), task C02b,
against MATLAB R2024b through the oracle harness, and their rules (no MATLAB).

The estimate's cases reuse C01's MATLAB references (``test_references.py``: the same code and
inputs, so the same keys). Its fits are the fitGaussian3D mex's nonlinear ones (parity D5), so
the sigmas carry D6.
"""

import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.detection import detect, run, sigma
from petakit5d.io import read_tiff
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.compare import EXACT, tol

_spec = importlib.util.spec_from_file_location(
    "detection_references_sigma", Path(__file__).with_name("test_references.py")
)
REF = importlib.util.module_from_spec(_spec)
REF.__file__ = str(_spec.origin) + ":helpers-sigma"
_spec.loader.exec_module(REF)

TOL_SIGMA = tol(reason="parity D6: the PSF sigma estimate (the mex's nonlinear fits, D5)", rel=1e-7)


def _py_estimate(rt):
    detect.clear_speedup_cache()  # the case starts a MATLAB session's state
    out = []
    for c in (1, 2):
        sxy, sz, xyz, _, _ = sigma.estimate_sigma3d(os.path.join(rt, f"psf_ch{c}.tif"))
        out += [sxy, sz, {k: xyz[k] + 1 for k in ("x", "y", "z")}]
    return out


def _source(rt: str) -> str:
    return os.path.join(rt, REF.EXPERIMENT, "ch1") + os.sep


def _py_camera_and_background(rt):
    src = _source(rt)
    ct = run.camera_type(src)
    th = run.camera_threshold(src)
    first = sorted(os.listdir(os.path.join(src, "DS")))[0]
    bi = detect.estimate_background(read_tiff(os.path.join(src, "DS", first)).astype(np.float64), 1)
    return [ct, th, bi[None, :]]


# designed condition folders for the camera search: {name: ({relative path: lines}, source)};
# each tree has a settings file at its root, so the upward search never leaves the tree
CAMERA_TREES = {
    "orca": ({"ex/A_Settings.txt": ["Model :\tC11440-22C", 'Camera type = "Orca4.0"']}, "ex/ch1/"),
    "andor": ({"ex/A_Settings.txt": ['Camera type = "Andor"']}, "ex/ch1/"),
    "both": ({"ex/A_Settings.txt": ['Camera Type = "orca and ANDOR"']}, "ex/ch1/"),
    "unknown": ({"ex/A_Settings.txt": ['Camera type = "Hamamatsu"']}, "ex/ch1/"),
    "two_lines": (
        {"ex/A_Settings.txt": ['Camera type = "Andor"', 'Camera type = "Orca4.0"']},
        "ex/ch1/",
    ),
    "indented": (
        {"ex/A_Settings.txt": ['   Camera type = "Andor"', "\tCamera type = Orca"]},
        "ex/ch1/",
    ),
    "two_up": ({}, "ex/ch1/"),
    "in_source": ({"ex/ch1/x_settings.TXT": ["Camera type = Andor"]}, "ex/ch1"),
    "two_files": (
        {"ex/A_Settings.txt": ["Camera type = Andor"], "ex/B_Settings.txt": ["Camera type = Orca"]},
        "ex/ch1/",
    ),
    "dropped_by_dir": ({"ex/#Settings.txt": ["Camera type = Orca"]}, "ex/ch1/"),
    "not_a_settings_file": ({"ex/A_Settings.txt.bak": ["Camera type = Orca"]}, "ex/ch1/"),
}


def _make_tree(files):
    def make(folder):
        root = Path(folder)
        (root / "ex" / "ch1").mkdir(parents=True, exist_ok=True)
        with open(root / "Z_Settings.txt", "w", newline="\n") as f:
            f.write("Camera type = Andor 2\n")
        for rel, lines in files.items():
            p = root / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            with open(p, "w", newline="\n") as f:
                f.write("".join(ln + "\n" for ln in lines))

    return make


def _py_camera(rt, rel):
    src = os.path.join(rt, *rel.strip("/").split("/")) + (os.sep if rel.endswith("/") else "")
    return [run.camera_type(src), run.camera_threshold(src)]


def _cases():
    cond = {"cond": Tree(REF.make_condition, "cond")}
    cs = [
        Case(
            "estimate_sigma",
            code=REF.CLEAR_CODE
            + "rt = in{1}; [a, b, c] = GU_estimateSigma3D([rt 'psf_ch1.tif'], []); "
            + "[d, e, f] = GU_estimateSigma3D([rt 'psf_ch2.tif'], []); out = {a, b, c, d, e, f};",
            args=[InPath("cond")],
            nargout=6,
            inputs=cond,
            python=_py_estimate,
            compare=TOL_SIGMA,
        ),
        Case(
            "camera_and_background",
            code=REF.CLEAR_CODE
            + REF.LOAD
            + "ct = get_camera_type_from_setting_file(data); "
            + "th = XR_get_camera_standard_deviation_threshold(data); "
            + "bi = XR_estimate_overall_background_information(double(readtiff("
            + "data.framePathsDS{1}{1})), 'randseed', 1); out = {ct, th, bi};",
            args=[InPath("cond")],
            nargout=3,
            inputs=cond,
            python=_py_camera_and_background,
            compare=EXACT,
        ),
    ]
    for name, (files, rel) in CAMERA_TREES.items():
        cs.append(
            Case(
                f"camera.{name}",
                code="d = struct('source', [in{1} in{2}]); "
                + "out = {get_camera_type_from_setting_file(d), "
                + "XR_get_camera_standard_deviation_threshold(d)};",
                args=[InPath("tree"), rel],
                nargout=2,
                inputs={"tree": Tree(_make_tree(files), "tree")},
                python=_py_camera,
                compare=EXACT,
            )
        )
    cs.append(
        Case(
            "camera.no_camera_line",
            code="d = struct('source', [in{1} 'ex/ch1/']); "
            "out = {get_camera_type_from_setting_file(d)};",
            args=[InPath("tree")],
            inputs={
                "tree": Tree(
                    _make_tree({"ex/A_Settings.txt": ["Model : Orca", "Cam type = Orca"]}), "tree"
                )
            },
            python=lambda rt: run.camera_type(os.path.join(rt, "ex", "ch1") + os.sep),
            expect_error="MATLAB:TooManyOutputsDueToMissingBraces",
            py_error=(run.CameraSettingsError,),
        )
    )
    return cs


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- rules (no MATLAB) ----------------------------------------------------------------------------


def test_a_bead_gives_its_sigmas():
    """PointDetectionTest's estimateSigmaFromBead in intent: a single bead's sigmas come back"""
    detect.clear_speedup_cache()
    vol = synth.beads((41, 48, 48), [(20.2, 23.7, 24.4)], (2.6, 1.4), 2000.0, 100.0).astype(
        np.float64
    )
    vol = vol + 2.0 * synth.CounterRNG(2026, "C02b.bead").normal(vol.shape)
    sxy, sz, xyz, ps, mask = sigma.estimate_sigma3d(vol)
    assert sxy == pytest.approx(1.4, rel=0.03)
    assert sz == pytest.approx(2.6, rel=0.03)
    assert [xyz["x"], xyz["y"], xyz["z"]] == pytest.approx([24.4, 23.7, 20.2], abs=0.05)


def test_the_bead_image_reads_zeros_as_nan():
    v = np.array([[[0.0, 1.0], [2.0, 0.0]]])
    p = sigma.psf_volume(v)
    assert np.isnan(p[0, 0, 0]) and np.isnan(p[0, 1, 1]) and p[0, 0, 1] == 1.0
    assert sigma.z_aniso_sigma(1.5, 3.0, 2.5) == (1.5, 1.2)


def test_the_camera_search_needs_a_settings_file(tmp_path):
    """MATLAB's search ends at the root's '' (dir lists nothing there): no settings file raises
    MATLAB:subsdeldimmismatch (measured); a missing folder on the way too"""
    (tmp_path / "ex" / "ch1").mkdir(parents=True)
    if any(
        p.is_dir() and any(n.lower().endswith("settings.txt") for n in os.listdir(p))
        for p in [tmp_path, *tmp_path.parents]
    ):
        pytest.skip("a settings file above the test's folder")
    with pytest.raises(run.CameraSettingsError) as e:
        run.camera_type(str(tmp_path / "ex" / "ch1") + os.sep)
    assert e.value.identifier == "MATLAB:subsdeldimmismatch"
    with pytest.raises(run.CameraSettingsError):
        run.camera_type(str(tmp_path / "missing") + os.sep)


def test_the_camera_threshold_options(tmp_path):
    (tmp_path / "ex" / "ch1").mkdir(parents=True)
    (tmp_path / "ex" / "A_Settings.txt").write_bytes(b'Camera type = "Orca4.0"\n')
    src = str(tmp_path / "ex" / "ch1") + os.sep
    assert run.camera_threshold(src) == 3.0
    assert run.camera_threshold(src, scmos=4.5) == 4.5
    (tmp_path / "ex" / "A_Settings.txt").write_bytes(b"Camera type = Hamamatsu\n")
    assert np.isnan(run.camera_threshold(src))
