"""The examples of the point-detection demos (python/examples/phase3), end to end.

The examples' own logic is tested on small inputs: imports do no work, the dataset's folders
are found as the demos name them, each demo's steps call the port with the MATLAB demo's
parameters (in Python's form), and the step-by-step comparison of output trees applies the
parity rows' tolerances and reports the numbers. The demos themselves run end to end on a
cut-down synthetic dataset as oracle cases: the demo dataset's layout (one condition, the two
bead images, experiments of two channels) with 16 small frames per channel, stage-scanned
spots that deskew into Gaussian spots, channel 2 shifted by 2 voxels in x. MATLAB runs the
demo's steps and Python the example's; the results and every file are compared with the
tolerances of ``compare_phase3_trees`` (the Amira files apart: W15). The multiple-experiments
demo runs up to its step 3: MATLAB's lifetime analysis needs more tracks than such data has
(docs/detection_scope.md, section 6); its step 4 is compared on the demo dataset. The
comparison of the demo runs on the full demo dataset is
``examples/phase3/compare_phase3_trees.py`` (docs/parity/phase3.md).
"""

from __future__ import annotations

import importlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, InPath, Tree

PHASE3 = Path(__file__).resolve().parents[2] / "examples" / "phase3"
DATASET = "PetaKit5D_point_detection_demo_dataset"
DEMOS = {
    "by_step": "demo_point_detection_3d_by_step",
    "integrated": "demo_point_detection_3d_integrated",
    "multiple": "demo_point_detection_3d_multiple_experiments",
}


def _example(name: str):
    """the example module ``name`` of examples/phase3 (the folder first on the import path,
    and examples/phase1 for its _common, as the examples have them)"""
    if str(PHASE3) not in sys.path:
        sys.path.insert(0, str(PHASE3))
    importlib.import_module("_phase3")
    return importlib.import_module(name)


@pytest.fixture
def phase3(monkeypatch):
    """examples/phase3 on the import path (and phase1 for _common); returns an importer"""
    monkeypatch.syspath_prepend(str(PHASE3))
    monkeypatch.delenv("PETAKIT5D_DEMO_TIMINGS", raising=False)
    monkeypatch.delenv("PETAKIT5D_DEMO_SAVE_DIR", raising=False)
    loader = importlib.import_module("_phase3")
    monkeypatch.syspath_prepend(loader.PHASE1)
    return importlib.import_module


def example_modules() -> list[str]:
    return sorted(p.stem for p in PHASE3.glob("*.py") if not p.stem.startswith("_"))


# ---------------------------------------------------------------------------------------------
# a cut-down dataset: the demo dataset's layout with small synthetic frames
# ---------------------------------------------------------------------------------------------

COND = "cond_synth"
EXPERIMENTS = ("Ex1_synth_z0p5", "Ex2_synth_z0p5")
# (z, y, x) of a raw frame: 24 planes of a stage scan. Track processing drops the tracks that
# come within ceil(4 sigma_xy) voxels (6 here) of a border of the deskewed frame, in z too
RAW = (24, 56, 48)
XSTEP = 4.1  # x shift per plane of the deskew (dz cos(31.5 deg) / 0.104 is 4.0992)
N_FRAMES = 16
FRAME_MS = 2000
BACKGROUND = 100.0
NOISE = 4.0
SIGMA_Z, SIGMA_XY = 1.3, 1.3  # the spots' sigmas in the deskewed frames (voxels)
CH2_DX = 2.0  # channel 2 sees the spots 2 voxels further in x, at 0.7 times the amplitude
PSF_SIGMA = {1: (3.2, 1.3), 2: (3.8, 1.5)}  # (z, yx) of the bead images
# per experiment: (first frame, last frame, (z, y, u) at the first frame, the velocity per
# frame, amplitude, frames missing); u is x in the raw frame at the spot's own plane (its x in
# the deskewed frame is u + XSTEP z: over its planes, 2.6 on either side, the spot spans 10.7
# voxels of raw x on either side, which stay in the frame). The spots keep 5 voxels or more
# apart, and a frame has 10 of them or fewer (with which the bad-point filter keeps every
# detection, as measured)
TRACKS = {
    1: (
        (0, 15, (11.5, 11.0, 13.0), (0.02, 0.15, 0.1), 300.0, ()),
        (0, 15, (11.5, 11.0, 24.0), (-0.02, 0.1, 0.12), 340.0, ()),
        (0, 3, (11.5, 11.0, 35.0), (0.0, -0.12, -0.1), 280.0, ()),
        (0, 15, (12.5, 20.0, 13.0), (0.03, -0.1, 0.15), 320.0, ()),
        (2, 15, (12.5, 20.0, 24.0), (0.0, 0.12, 0.1), 300.0, (9,)),
        (0, 13, (11.8, 29.0, 13.5), (-0.03, 0.1, 0.1), 360.0, ()),
        (4, 15, (11.8, 29.0, 24.5), (0.02, -0.1, 0.12), 290.0, ()),
        (0, 15, (11.8, 29.0, 35.0), (0.0, 0.12, -0.12), 330.0, ()),
        (0, 4, (12.8, 38.0, 13.0), (0.02, 0.15, 0.1), 300.0, ()),
        (5, 15, (12.8, 38.0, 24.0), (-0.02, -0.1, 0.1), 310.0, ()),
        (0, 15, (12.2, 47.0, 13.0), (0.0, -0.1, 0.12), 350.0, ()),
        (3, 12, (12.2, 47.0, 24.0), (0.03, 0.08, -0.1), 270.0, ()),
    ),
    2: (
        (0, 15, (11.8, 12.0, 16.0), (0.02, 0.12, 0.1), 320.0, ()),
        (0, 12, (12.4, 22.0, 29.0), (-0.02, -0.1, 0.1), 290.0, ()),
        (2, 15, (12.0, 31.0, 15.0), (0.0, 0.1, 0.12), 350.0, (8,)),
        (0, 15, (12.6, 40.0, 27.0), (0.03, 0.1, -0.1), 310.0, ()),
        (3, 15, (11.6, 47.0, 18.0), (0.0, -0.1, 0.1), 270.0, ()),
    ),
}
# the noise's stream (synth.CounterRNG): one whose weak noise detections (a few, amplitude 6)
# all have a residual sigma above the camera threshold; MATLAB's runner left out one below it
# that its detector returns (docs/parity/phase3.md, open item)
STREAM = "phase3h"


def spots(experiment: int, frame: int, channel: int) -> list[tuple[tuple[float, ...], float]]:
    """the spots of a frame: ((z, y, u), amplitude)"""
    out = []
    for first, last, p0, v, a, missing in TRACKS[experiment]:
        if first <= frame <= last and frame not in missing:
            t = frame - first
            p = [p0[i] + v[i] * t for i in range(3)]
            if channel == 2:
                p[2] += CH2_DX
            out.append((tuple(p), a if channel == 1 else 0.7 * a))
    return out


def _round_uint16(v: np.ndarray) -> np.ndarray:
    t = np.trunc(v)
    v = t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)
    return np.clip(v, 0, 65535).astype(np.uint16)


def raw_frame(experiment: int, frame: int, channel: int) -> np.ndarray:
    """a raw stage-scan frame (z, y, x), uint16: each spot sheared as the deskew un-shears it
    (plane z of the deskewed frame is the raw plane shifted by XSTEP z in x), noise from a
    bit-stable stream"""
    z = np.arange(RAW[0], dtype=np.float64).reshape(-1, 1, 1)
    y = np.arange(RAW[1], dtype=np.float64).reshape(1, -1, 1)
    x = np.arange(RAW[2], dtype=np.float64).reshape(1, 1, -1)
    vol = np.full(RAW, BACKGROUND)
    for (z0, y0, u0), amp in spots(experiment, frame, channel):
        dz, dy = z - z0, y - y0
        dx = x - u0 + XSTEP * dz
        e = -(dz * dz / (2 * SIGMA_Z * SIGMA_Z) + (dy * dy + dx * dx) / (2 * SIGMA_XY * SIGMA_XY))
        vol = vol + amp * synth.stable_exp(e)
    noise = synth.CounterRNG(2026, f"{STREAM}/ex{experiment}/ch{channel}/{frame}").normal(RAW)
    return _round_uint16(vol + NOISE * noise)


def psf_volume(channel: int) -> np.ndarray:
    vol = synth.beads((40, 48, 48), [(20.0, 24.0, 24.0)], PSF_SIGMA[channel], 1000.0, BACKGROUND)
    noise = synth.CounterRNG(2026, f"phase3/psf/ch{channel}").normal(vol.shape)
    return _round_uint16(vol + 5.0 * noise)


def frame_name(channel: int, frame: int) -> str:
    cam, wl = ("B", 560) if channel == 1 else ("A", 642)
    t = frame * FRAME_MS
    return (
        f"synth_Cam{cam}_ch{channel - 1}_stack{frame:04d}_{wl}nm_{t:07d}msec_"
        f"{1000000 + t:010d}msecAbs.tif"
    )


def make_dataset(folder, experiments: int = 1) -> None:
    """the dataset folder (``folder/PetaKit5D_point_detection_demo_dataset/``): one condition
    with the bead images and ``experiments`` experiments, each with its settings file (the
    camera, for the detection's residual threshold) and 16 frames of two channels"""
    cond = Path(folder) / DATASET / COND
    for c, name in ((1, "560totalPSF.tif"), (2, "642totalPSF.tif")):
        cond.mkdir(parents=True, exist_ok=True)
        synth.write_tiff(cond / name, psf_volume(c))
    for e in range(1, experiments + 1):
        ex = cond / EXPERIMENTS[e - 1]
        for c in (1, 2):
            (ex / f"ch{c}").mkdir(parents=True, exist_ok=True)
            for k in range(N_FRAMES):
                synth.write_tiff(ex / f"ch{c}" / frame_name(c, k), raw_frame(e, k, c))
        with open(ex / "synth_Settings.txt", "w", newline="\n") as f:
            f.write('Camera type = "Orca4.0"\n')


# ---------------------------------------------------------------------------------------------
# the demos end to end on the cut-down dataset: MATLAB's demo steps against the examples'
# ---------------------------------------------------------------------------------------------

# the functions of the demos' steps that keep state between calls (detection's t thresholds
# are cached per session); cleared so that each case starts as a new MATLAB session
CLEAR = (
    "clear XR_pointSourceDetection3D XR_pointSourceDetection3D_3 XR_fitGaussians3D "
    "XR_fitGaussian3D XR_fitGaussian3D_Ac fitGaussian3D estGaussianAmplitude3D XR_locmax3d "
    "filter_bad_points_in_detection XR_estimate_overall_background_information "
    "XR_get_camera_standard_deviation_threshold get_camera_type_from_setting_file "
    "XR_runDetection3D XR_runDetection3D_frame_parallel XR_runTracking3D XR_runTrackProcessing3D "
    "trackCloseGapsKalmanSparse linkFeaturesKalmanSparse costMatLinearMotionLink2 "
    "costMatLinearMotionCloseGaps2 kalmanGainLinearMotion kalmanInitLinearMotion "
    "getAveDispEllipseAll2 estimTrackTypeParamLM2 rotateTracks3D GU_estimateSigma3D "
    "XR_loadConditionData3D XR_correctXZoffsetData3D XZoffsetEstimation3D deskewData "
    "XR_runDetTrack3d; "
)
# the dataset's folders, as the examples find them; the demo's steps 1 and 2
LOAD = (
    "r = in{1}; ds = [r 'PetaKit5D_point_detection_demo_dataset']; d = dir(ds); "
    "d = d([d.isdir] & ~startsWith({d.name}, '.')); cond = [ds '/' d(1).name]; "
    "e = dir([cond '/Ex1*']); e = e([e.isdir]); rt = [cond '/' e(1).name]; "
)
STEPS_1_2 = (
    "data = XR_loadConditionData3D(rt, {'ch1', 'ch2'}, {'gfp', 'rfp'}); data0 = data; "
    "[sigmaXY_ch1, sigmaZ_ch1] = GU_estimateSigma3D([cond '/560totalPSF.tif'], []); "
    "[sigmaXY_ch2, sigmaZ_ch2] = GU_estimateSigma3D([cond '/642totalPSF.tif'], []); "
    "sigma_mat = [sigmaXY_ch1, sigmaZ_ch1 ./ data(1).zAniso; "
    "sigmaXY_ch2, sigmaZ_ch2 ./ data(1).zAniso]; "
)
OUT = "out = {rmfield(data0, 'date'), sigma_mat};"
BY_STEP = (
    CLEAR + LOAD + STEPS_1_2 + "aname = 'Analysis_update_3'; "
    "data = deskewData(data, 'Overwrite', true, 'SkewAngle', 31.5, 'aname', aname, "
    "'Rotate', false,'sCMOSCameraFlip', false, 'Crop', false); "
    "XR_correctXZoffsetData3D(data, 'Overwrite', true); "
    "apath = arrayfun(@(i) [i.source filesep aname filesep], data, 'unif', 0); "
    "XR_runDetection3D(data, 'Sigma', sigma_mat, 'Overwrite', true, 'WindowSize', [], "
    "'ResultsPath', apath,'Mode', 'xyzAc', 'FitMixtures', false, 'MaxMixtures', 1, "
    "'DetectionMethod', 'lowSNR', 'BackgroundCorrection', true); "
    "data = data(~([data.movieLength]==1)); "
    "XR_runTracking3D(data, loadTrackSettings('Radius', [3, 6], 'MaxGapLength', 2), "
    "'FileName', 'trackedFeatures.mat', 'Overwrite', true, 'ResultsPath', apath); "
    "XR_runTrackProcessing3D(data, 'Overwrite', true, 'TrackerOutput', 'trackedFeatures.mat', "
    "'FileName', 'ProcessedTracks.mat', 'Buffer', [3, 3], 'BufferAll', false, "
    "'FitMixtures', false, 'WindowSize', [], 'ResultsPath', apath); "
    "rotateTracks3D(data, 'Overwrite', true, 'Crop', false, 'ResultsPath', apath); " + OUT
)
INTEGRATED = (
    CLEAR
    + LOAD
    + STEPS_1_2
    + "XR_runDetTrack3d(data, 'Sigma', sigma_mat, 'DetectionMethod', 'lowSNR', 'aname', "
    "'Analysis_update_3', 'FitMixtures', false, 'BackgroundCorrection', true, "
    "'XZoffsetCorrection', true); " + OUT
)
MULTIPLE = (
    CLEAR
    + LOAD
    + "rt = cond; "
    + STEPS_1_2
    + "XR_runDetTrack3d(data, 'Sigma', sigma_mat, 'DetectionMethod', 'lowSNR', 'aname', "
    "'Analysis_update_3', 'FitMixtures', false, 'BackgroundCorrection', true, "
    "'zoffsetCorrection', true); " + OUT
)


def posix_paths(v):
    """the value with ``/`` for every ``\\`` in its texts: paths the port made with Windows'
    filesep (MATLAB's there) as the references, made on Linux, hold them"""
    if isinstance(v, str):
        return v.replace("\\", "/")
    if isinstance(v, dict):
        return {k: posix_paths(x) for k, x in v.items()}
    if isinstance(v, list | tuple):
        return type(v)(posix_paths(x) for x in v)
    if isinstance(v, np.ndarray) and v.dtype == object:
        out = np.empty(v.shape, dtype=object)
        for i, x in np.ndenumerate(v):
            out[i] = posix_paths(x)
        return out
    from petakit5d.io.matfile import StructArray

    if isinstance(v, StructArray):
        return StructArray(v.shape, v.fields, [posix_paths(e) for e in v])
    return v


class _DemoData:
    """the examples' environment for a run on ``root``: PETAKIT5D_DEMO_DATA, the in-memory
    results in ``saved/``, no timings file, the detection's caches cleared (a new session)"""

    def __init__(self, root: str):
        self.root = root
        self.env = {
            "PETAKIT5D_DEMO_DATA": root,
            "PETAKIT5D_DEMO_SAVE_DIR": os.path.join(os.getcwd(), "saved"),
            "PETAKIT5D_DEMO_TIMINGS": None,
        }

    def __enter__(self):
        from petakit5d.detection import detect

        self.old = {k: os.environ.get(k) for k in self.env}
        for k, v in self.env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        detect.clear_threshold_cache()
        detect.clear_speedup_cache()
        return self

    def __exit__(self, *exc):
        for k, v in self.old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def saved(self, demo: str) -> list:
        """the condition (without ``date``, which comes from the folders' digits) and the PSF
        sigmas, as the example saved them"""
        from petakit5d.io.matfile import StructArray, read_mat

        folder = os.path.join(self.env["PETAKIT5D_DEMO_SAVE_DIR"], DEMOS[demo])
        data = posix_paths(read_mat(os.path.join(folder, "load.mat"))["data"])
        if isinstance(data, dict):
            data = {k: v for k, v in data.items() if k != "date"}
        else:
            rows = [{k: v for k, v in e.items() if k != "date"} for e in data]
            data = StructArray(data.shape, [f for f in data.fields if f != "date"], rows)
        return [data, read_mat(os.path.join(folder, "sigma.mat"))["sigma_mat"]]


def _run_main(demo: str):
    def run(root):
        with _DemoData(root) as env:
            _example(DEMOS[demo]).main()
            return env.saved(demo)

    return run


def _run_multiple(root):
    """the multiple-experiments demo's steps 1 to 3 (the lifetime analysis needs more tracks
    than this dataset has: MATLAB's stops at an empty category, docs/detection_scope.md)"""
    with _DemoData(root) as env:
        p3 = _example("_phase3")
        ex = _example(DEMOS["multiple"])
        data = p3.load_data(p3.condition_path())
        p3.save_values(DEMOS["multiple"], "load", {"data": data})
        sigma = p3.psf_sigmas(data, p3.psf_paths())
        p3.save_values(DEMOS["multiple"], "sigma", {"sigma_mat": sigma})
        ex.detection_tracking(data, sigma)
        return env.saved("multiple")


def _case_tree():
    """the files compared: not the examples' in-memory results (``saved/``), nor the empty
    folder ``<source>Analysis/`` that MATLAB's XR_runTrackProcessing3D makes beside the given
    results folder (docs/parity/phase3.md); the Amira files print the fits' values (W15). On
    Windows not the frame-parallel runner's per-frame files (``Detections/``) either: the
    references are Linux runs, and MATLAB's XR_runDetTrack3d takes 'parfor' on a PC, whose
    XR_runDetection3D writes none, as run_detection_tracking does"""
    from petakit5d.testing.compare import ignore

    cmp = _example("compare_phase3_trees")
    per_frame = (r"/Detections(/|$)",) if os.name == "nt" else ()
    return cmp.tree_spec(
        ignore=(r"^saved/", r"/ch1/Analysis$", *per_frame),
        amira=ignore("W15 (approved): the Amira files print the fits' values (W11)"),
    )


def _cases() -> list[Case]:
    from petakit5d.testing.compare import EXACT

    cmp = _example("compare_phase3_trees")
    out = []
    for name, code, python, n in (
        ("demo.by_step", BY_STEP, _run_main("by_step"), 1),
        ("demo.integrated", INTEGRATED, _run_main("integrated"), 1),
        ("demo.multiple", MULTIPLE, _run_multiple, 2),
    ):
        out.append(
            Case(
                name,
                code=code,
                args=[InPath("data")],
                nargout=2,
                python=python,
                inputs={"data": Tree(lambda p, n=n: make_dataset(p, n), "data")},
                compare={0: EXACT, 1: cmp.tolerances()["*/sigma.mat"]},
                tree=_case_tree(),
                marks=("slow",),
            )
        )
    return out


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path, phase3):
    oracle.check(case, tmp_path)


def test_the_cut_down_dataset_deskews_into_round_spots():
    """the raw frames' spots, deskewed, are Gaussian spots at their places (z, y, u + XSTEP z),
    and channel 2's 2 voxels further in x"""
    from petakit5d.deskew_rotate.volume import deskew_volume

    for c in (1, 2):
        raw = raw_frame(1, 6, c).astype(np.float64)
        ds = deskew_volume(raw, 31.5, 0.5, 0.104, False, crop=False)
        assert ds.shape == (24, 56, 48 + 95)
        for (z0, y0, u0), _ in spots(1, 6, c):
            x0 = u0 + XSTEP * z0
            zc, yc, xc = synth.centroid(ds - BACKGROUND, (z0, y0, x0), 2)
            assert abs(zc - z0) < 0.35 and abs(yc - y0) < 0.35 and abs(xc - x0) < 0.35


def _results(root: Path) -> Path:
    return root / DATASET / COND / EXPERIMENTS[0] / "ch1" / "Analysis_update_3"


@pytest.mark.slow
def test_the_demos_run_on_the_cut_down_data(phase3, tmp_path, monkeypatch):
    """without MATLAB (the oracle cases need the golden store to decide): the step-by-step and
    the integrated examples on two copies of the cut-down dataset write MATLAB's files and the
    same results, every spot is detected where it is, and the tracks are the spots' tracks"""
    from petakit5d.detection import detect
    from petakit5d.io.matfile import read_mat
    from petakit5d.testing.compare import compare
    from petakit5d.testing.compare import read_mat as read_values

    roots = {}
    for demo in ("by_step", "integrated"):
        roots[demo] = tmp_path / demo
        make_dataset(roots[demo], 1)
        monkeypatch.setenv("PETAKIT5D_DEMO_DATA", str(roots[demo]))
        detect.clear_threshold_cache()
        detect.clear_speedup_cache()
        phase3(DEMOS[demo]).main()
    a = _results(roots["by_step"])
    files = sorted(p.relative_to(a).as_posix() for p in a.rglob("*") if p.is_file())
    amira = [f"Amira/Tracking/tracking_t_{k:04d}.am" for k in range(1, N_FRAMES + 1)]
    assert [f for f in files if f.startswith("Amira/")] == amira
    rest = [f for f in files if not f.startswith(("Amira/", "Masks/"))]
    assert rest == [
        "Detection3D.mat",
        "ProcessedTracks.mat",
        "RotatedTracks.mat",
        "dynRange.mat",
        "settings.mat",
        "trackedFeatures.mat",
    ]
    assert len([f for f in files if f.startswith("Masks/")]) == N_FRAMES
    # the integrated example: the same files and values (and the frame-parallel runner's)
    b = _results(roots["integrated"])
    for name in rest:
        rep = compare(read_values(str(a / name)), read_values(str(b / name)), path=name)
        assert rep.ok, rep.text(name)
    # its per-frame files: none on Windows, where MATLAB's XR_runDetTrack3d takes 'parfor'
    # (XR_runDetection3D) and so does run_detection_tracking
    per_frame = 0 if os.name == "nt" else N_FRAMES
    assert len(list((b / "Detections").glob("Detection_F_*.mat"))) == per_frame
    # every spot detected within 0.5 voxels (MATLAB's 1-based positions in the deskewed frame,
    # which the offset correction padded by channel 2's shift in x)
    fi = read_mat(str(a / "Detection3D.mat"))["frameInfo"]
    for k, f in enumerate(fi):
        got = np.stack([np.asarray(f[c], dtype=np.float64)[0] for c in ("z", "y", "x")], 1)
        for (z0, y0, u0), _ in spots(1, k, 1):
            want = np.array([z0, y0, u0 + XSTEP * z0 + CH2_DX]) + 1
            assert np.min(np.linalg.norm(got - want, axis=1)) < 0.5, (k, want)
    # each synthetic track of two frames or more is a track of the tracker, from its first frame
    tracks = read_mat(str(a / "trackedFeatures.mat"))["tracksFinal"]
    starts = {int(np.asarray(t["seqOfEvents"])[0, 0]) for t in tracks}
    assert {first + 1 for first, *_ in TRACKS[1]} <= starts


# ---------------------------------------------------------------------------------------------
# the examples themselves
# ---------------------------------------------------------------------------------------------


def test_the_examples_are_these():
    """one example per MATLAB demo (demos/demo_point_detection_3d_*.m), and the comparison"""
    assert example_modules() == sorted(["compare_phase3_trees", *DEMOS.values()])


@pytest.mark.parametrize("name", example_modules())
def test_example_imports_without_running(phase3, name):
    """importing an example does no work: everything runs in main()"""
    before = set(sys.modules)
    mod = phase3(name)
    assert callable(getattr(mod, "main", None))
    heavy = {m for m in set(sys.modules) - before if m.split(".")[0] in ("numpy", "scipy")}
    assert not heavy, f"{name} imports {sorted(heavy)[:3]} at import time"


def test_phase1_is_put_first_on_the_path(phase3, monkeypatch):
    p3 = phase3("_phase3")
    monkeypatch.setattr(sys, "path", [str(PHASE3), "elsewhere", p3.PHASE1])
    importlib.reload(p3)
    assert sys.path == [p3.PHASE1, str(PHASE3), "elsewhere"]
    assert Path(p3.PHASE1) == PHASE3.parent / "phase1"


def _layout(root: Path, experiments=("Ex1_a_z0p5", "Ex2_b_z0p5", "Ex10_c_z0p5")) -> Path:
    """a condition folder with experiments, a hidden folder and a second condition (the first
    by name is the demos')"""
    cond = root / DATASET / "some_condition"
    for e in experiments:
        (cond / e).mkdir(parents=True)
    (root / DATASET / ".hidden").mkdir()
    (root / DATASET / "z_other_condition").mkdir()
    return cond


@pytest.mark.parametrize("named", ["parent", "parent/", "dataset"])
def test_the_dataset_folders(phase3, monkeypatch, tmp_path, named):
    """PETAKIT5D_DEMO_DATA names the dataset's folder or its parent; the condition is its only
    folder, experiment n the folder Ex<n>... (not Ex1 for Ex10), the bead images in the
    condition folder; roots without a trailing separator, as the demos write them"""
    p3 = phase3("_phase3")
    cond = _layout(tmp_path)
    value = {"parent": str(tmp_path), "parent/": str(tmp_path) + "/", "dataset": str(cond.parent)}
    monkeypatch.setenv("PETAKIT5D_DEMO_DATA", value[named])
    root = tmp_path.as_posix() + "/" + DATASET
    assert p3.dataset_path() == root
    assert p3.condition_path() == root + "/some_condition"
    assert p3.experiment_path(1) == root + "/some_condition/Ex1_a_z0p5"
    assert p3.experiment_path(2) == root + "/some_condition/Ex2_b_z0p5"
    assert p3.experiment_path(10) == root + "/some_condition/Ex10_c_z0p5"
    assert p3.psf_paths() == [
        root + "/some_condition/560totalPSF.tif",
        root + "/some_condition/642totalPSF.tif",
    ]
    with pytest.raises(SystemExit, match="no experiment 3"):
        p3.experiment_path(3)


def test_a_missing_dataset_is_named(phase3, monkeypatch, tmp_path):
    p3 = phase3("_phase3")
    monkeypatch.setenv("PETAKIT5D_DEMO_DATA", str(tmp_path))
    with pytest.raises(SystemExit, match="is not in"):
        p3.dataset_path()
    (tmp_path / DATASET).mkdir()
    with pytest.raises(SystemExit, match="no condition folder"):
        p3.condition_path()


@pytest.fixture
def tiny(tmp_path, monkeypatch):
    """a dataset of two experiments of two frames (3 x 4 x 5), PETAKIT5D_DEMO_DATA set"""
    cond = tmp_path / DATASET / COND
    for e, n in (("Ex1_synth_z0p5", 2), ("Ex2_synth_z0p5", 2)):
        for c in (1, 2):
            d = cond / e / f"ch{c}"
            d.mkdir(parents=True)
            for k in range(n):
                synth.write_tiff(d / frame_name(c, k), np.full((3, 4, 5), 100 + k, np.uint16))
    monkeypatch.setenv("PETAKIT5D_DEMO_DATA", str(tmp_path))
    return tmp_path


def test_the_condition_is_loaded_with_the_demos_answers(phase3, tiny):
    """the channels ch1 and ch2 with the markers gfp and rfp, what the demos' comments say to
    enter at XR_loadConditionData3D's prompts"""
    p3 = phase3("_phase3")
    data = p3.load_data(p3.experiment_path(1))
    assert len(data) == 1
    m = data[0]
    assert m.markers == ["gfp", "rfp"]
    assert [os.path.basename(c.rstrip("/\\")) for c in m.channels] == ["ch1", "ch2"]
    assert m.source == m.channels[0] and m.movie_length == 2
    assert len(p3.load_data(p3.condition_path())) == 2


def test_the_psf_sigmas_take_sigma_z_over_the_z_anisotropy(phase3, tiny, monkeypatch):
    """[sigmaXY, sigmaZ ./ data(1).zAniso], a row per channel, from GU_estimateSigma3D"""
    import petakit5d.detection.sigma as sigma

    p3 = phase3("_phase3")
    data = p3.load_data(p3.experiment_path(1))
    seen = []
    est = {"a.tif": (1.25, 3.1), "b.tif": (1.5, 3.9)}
    monkeypatch.setattr(
        sigma, "estimate_sigma3d", lambda f: seen.append(f) or (*est[f], {}, {}, None)
    )
    s = p3.psf_sigmas(data, ["a.tif", "b.tif"])
    za = data[0].z_aniso
    assert seen == ["a.tif", "b.tif"]
    assert s.dtype == np.float64 and s.shape == (2, 2)
    assert s.tolist() == [[1.25, 3.1 / za], [1.5, 3.9 / za]]
    # dz 0.5 from the folder name: 0.5 sin(31.5 pi / 180) / 0.104
    assert za == pytest.approx(0.5 * np.sin(31.5 * np.pi / 180) / 0.104, rel=1e-15)


def test_in_memory_results_are_saved_as_mat_files_only_when_asked(phase3, tiny, monkeypatch):
    """as MATLAB's save writes them by default (version 7), a condition as its struct array"""
    from petakit5d.io.matfile import mat_version, read_mat

    p3 = phase3("_phase3")
    data = p3.load_data(p3.condition_path())
    out = tiny / "saved"
    p3.save_values("demo_x", "load", {"data": data, "v": np.eye(2)})
    assert not out.exists()
    monkeypatch.setenv("PETAKIT5D_DEMO_SAVE_DIR", str(out))
    p3.save_values("demo_x", "load", {"data": data, "v": np.eye(2)})
    assert mat_version(str(out / "demo_x" / "load.mat")) == "5"
    got = read_mat(str(out / "demo_x" / "load.mat"))
    assert [d["source"] for d in got["data"]] == [m.source for m in data]
    assert got["v"].tolist() == [[1.0, 0.0], [0.0, 1.0]]


def test_lists_of_structs_are_saved_as_struct_arrays(phase3):
    """MATLAB's cohorts are a 1 x n struct array; one element is a struct"""
    from petakit5d.io.matfile import StructArray

    p3 = phase3("_phase3")
    two = p3._struct_array([{"a": 1.0, "b": 2.0}, {"a": 3.0, "b": 4.0}])
    assert isinstance(two, StructArray) and two.shape == (1, 2) and two.fields == ("a", "b")
    assert p3._struct_array([{"a": 1.0}]) == {"a": 1.0}
    mixed = [{"a": 1.0}, {"b": 2.0}]
    assert p3._struct_array(mixed) is mixed
    assert p3._struct_array([]) == [] and p3._struct_array(5.0) == 5.0


# ---------------------------------------------------------------------------------------------
# each demo's steps call the port with the MATLAB demo's parameters
# ---------------------------------------------------------------------------------------------

STAGES = {
    "petakit5d.detection": (
        "deskew_data",
        "correct_xz_offset_data",
        "run_detection",
        "run_track_processing",
        "rotate_tracks",
        "run_detection_tracking",
    ),
    "petakit5d.detection.tracking": ("run_tracking", "track_settings"),
    "petakit5d.detection.lifetime": ("run_lifetime_analysis", "plot_intensity_cohorts"),
}


@pytest.fixture
def spied(monkeypatch, tiny):
    """the heavy stages replaced by recorders [(name, args, kwargs)]; the PSF sigmas fixed
    (1.25, 3.1) and (1.5, 3.9)"""
    import petakit5d.detection.sigma as sigma

    calls = []
    for mod, names in STAGES.items():
        m = importlib.import_module(mod)
        for n in names:
            monkeypatch.setattr(
                m, n, lambda *a, _n=n, **k: calls.append((_n, a, k)) or _n, raising=False
            )
    est = {"560totalPSF.tif": (1.25, 3.1), "642totalPSF.tif": (1.5, 3.9)}
    monkeypatch.setattr(
        sigma, "estimate_sigma3d", lambda f: (*est[os.path.basename(f)], {}, {}, None)
    )
    return calls


def _posix(path: str) -> str:
    """a path with slashes (the port joins with the platform's separator, as MATLAB)"""
    return path.replace("\\", "/")


def _sigma(data) -> list:
    za = data[0].z_aniso
    return [[1.25, 3.1 / za], [1.5, 3.9 / za]]


def _sources(arg) -> list[str]:
    rows = [arg] if isinstance(arg, dict) else list(arg)
    return [r["source"] if isinstance(r, dict) else r.source for r in rows]


def test_by_step_calls_the_stages_with_the_demos_parameters(phase3, spied):
    from petakit5d.core.paths import filesep

    ex = phase3(DEMOS["by_step"])
    ex.main()
    names = [c[0] for c in spied]
    assert names == [
        "deskew_data",
        "correct_xz_offset_data",
        "run_detection",
        "track_settings",
        "run_tracking",
        "run_track_processing",
        "rotate_tracks",
    ]
    calls = {c[0]: c for c in spied}
    data = calls["deskew_data"][1][0]
    assert len(data) == 1 and _posix(data[0].source).endswith("Ex1_synth_z0p5/ch1/")
    apath = [data[0].source + filesep + "Analysis_update_3" + filesep]
    assert calls["deskew_data"][2] == {
        "overwrite": True,
        "skew_angle": 31.5,
        "aname": "Analysis_update_3",
        "rotate": False,
        "s_cmos_camera_flip": False,
        "crop": False,
    }
    assert calls["correct_xz_offset_data"][1] == (data,)
    assert calls["correct_xz_offset_data"][2] == {"overwrite": True}
    det = calls["run_detection"]
    assert det[1] == (data,)
    assert det[2]["sigma"].tolist() == _sigma(data)
    assert {k: v for k, v in det[2].items() if k != "sigma"} == {
        "overwrite": True,
        "window_size": None,
        "results_path": apath,
        "mode": "xyzAc",
        "fit_mixtures": False,
        "max_mixtures": 1,
        "detection_method": "lowSNR",
        "background_correction": True,
    }
    assert calls["track_settings"][2] == {"radius": (3, 6), "max_gap_length": 2}
    trk = calls["run_tracking"]
    assert _sources(trk[1][0]) == [data[0].source] and trk[1][1] == "track_settings"
    assert trk[2] == {
        "file_name": "trackedFeatures.mat",
        "overwrite": True,
        "results_path": apath,
    }
    prc = calls["run_track_processing"]
    assert prc[2] == {
        "overwrite": True,
        "tracker_output": "trackedFeatures.mat",
        "file_name": "ProcessedTracks.mat",
        "buffer": (3, 3),
        "buffer_all": False,
        "fit_mixtures": False,
        "window_size": None,
        "results_path": apath,
    }
    assert calls["rotate_tracks"][2] == {"overwrite": True, "crop": False, "results_path": apath}


def test_by_step_results_folders_are_the_demos(phase3, tiny, monkeypatch):
    """[i.source filesep aname filesep]: the source ends with a separator, so two follow each
    other, with MATLAB's filesep on each platform"""
    import petakit5d.core.paths as paths

    p3, ex = phase3("_phase3"), phase3(DEMOS["by_step"])
    data = p3.load_data(p3.condition_path())
    for sep in ("/", "\\"):
        monkeypatch.setattr(paths, "filesep", sep)
        assert ex.results_paths(data) == [m.source + sep + "Analysis_update_3" + sep for m in data]


def test_by_step_tracks_each_movie_into_its_own_folder(phase3, tiny):
    """step 5 drops the one-frame movies (tracking needs two frames) with their folders: the
    demo drops only the movies, so the later movies would get the folders before theirs (the
    slip W3 of XR_runDetTrack3d)"""
    import dataclasses

    from petakit5d.detection.condition import Condition

    p3, ex = phase3("_phase3"), phase3(DEMOS["by_step"])
    two = p3.load_data(p3.condition_path())
    one = dataclasses.replace(two[0], movie_length=1)
    data = Condition([two[0], one, two[1]])
    movies, apath = ex.tracked_movies(data, ["a/", "b/", "c/"])
    assert _sources(movies) == [two[0].source, two[1].source]
    assert apath == ["a/", "c/"]
    movies, apath = ex.tracked_movies(Condition([two[0]]), ["a/"])
    assert isinstance(movies, dict) and apath == ["a/"]


def test_integrated_calls_the_wrapper_with_the_demos_parameters(phase3, spied):
    phase3(DEMOS["integrated"]).main()
    (call,) = spied
    assert call[0] == "run_detection_tracking"
    data = call[1][0]
    assert len(data) == 1 and _posix(data[0].source).endswith("Ex1_synth_z0p5/ch1/")
    assert call[2]["sigma"].tolist() == _sigma(data)
    assert {k: v for k, v in call[2].items() if k != "sigma"} == {
        "detection_method": "lowSNR",
        "aname": "Analysis_update_3",
        "fit_mixtures": False,
        "background_correction": True,
        "x_zoffset_correction": True,
    }


def test_multiple_experiments_calls_the_wrapper_without_the_offsets(phase3, spied):
    """all experiments of the condition; MATLAB's 'zoffsetCorrection' is no option of
    XR_runDetTrack3d (its parser drops it), so the offsets stay uncorrected"""
    phase3(DEMOS["multiple"]).main()
    call = spied[0]
    assert call[0] == "run_detection_tracking"
    data = call[1][0]
    assert [_posix(m.source).rstrip("/").split("/")[-2] for m in data] == [
        "Ex1_synth_z0p5",
        "Ex2_synth_z0p5",
    ]
    assert {k: v for k, v in call[2].items() if k != "sigma"} == {
        "detection_method": "lowSNR",
        "aname": "Analysis_update_3",
        "fit_mixtures": False,
        "background_correction": True,
    }


def test_multiple_experiments_runs_the_lifetime_analysis(phase3, spied, monkeypatch, tmp_path):
    """step 4: runLifetimeAnalysis3D(data, 'AnalysisPath', aname, 'overwrite', true) and
    plotIntensityCohorts(data, 'AnalysisPath', aname) on the condition's movies, after the
    tracks; lftRes and the cohorts (a struct array) saved"""
    import petakit5d.detection.lifetime as lifetime
    from petakit5d.io.matfile import read_mat

    cohorts = [
        {"bounds": np.array([[1.0, 2.0]]), "t": 2.0},
        {"bounds": np.zeros((0, 0)), "t": 3.0},
    ]

    def analysis(data, **kw):
        spied.append(("lft", data, kw))
        return {"v": 1.0}

    def cohort_numbers(data, **kw):
        spied.append(("coh", data, kw))
        return cohorts

    monkeypatch.setattr(lifetime, "run_lifetime_analysis", analysis)
    monkeypatch.setattr(lifetime, "plot_intensity_cohorts", cohort_numbers)
    monkeypatch.setenv("PETAKIT5D_DEMO_SAVE_DIR", str(tmp_path / "saved"))
    phase3(DEMOS["multiple"]).main()
    assert [c[0] for c in spied] == ["run_detection_tracking", "lft", "coh"]
    (det, lft, coh) = spied
    assert lft[2] == {"analysis_path": "Analysis_update_3", "overwrite": True}
    assert coh[2] == {"analysis_path": "Analysis_update_3"}
    assert lft[1] is det[1][0] and coh[1] is det[1][0]  # the loaded condition, as in the demo
    for c in (lft, coh):
        assert [_posix(s).rstrip("/").split("/")[-2] for s in _sources(c[1])] == [
            "Ex1_synth_z0p5",
            "Ex2_synth_z0p5",
        ]
    saved = tmp_path / "saved" / DEMOS["multiple"]
    assert read_mat(str(saved / "lifetime.mat"))["lftRes"]["v"] == 1.0
    back = read_mat(str(saved / "cohorts.mat"))["cohorts"]
    assert back.shape == (1, 2) and [e["t"] for e in back] == [2.0, 3.0]


def test_steps_are_timed_and_results_saved(phase3, spied, monkeypatch, tmp_path):
    from petakit5d.io.matfile import read_mat

    log = tmp_path / "t.jsonl"
    monkeypatch.setenv("PETAKIT5D_DEMO_TIMINGS", str(log))
    monkeypatch.setenv("PETAKIT5D_DEMO_SAVE_DIR", str(tmp_path / "saved"))
    phase3(DEMOS["by_step"]).main()
    rows = [json.loads(x) for x in log.read_text().splitlines()]
    assert [r["step"] for r in rows] == [
        "Step 1: load data",
        "Step 2: PSF sigmas",
        "Step 3: deskew",
        "Step 3: XZ offsets",
        "Step 4: detection",
        "Step 5.1: tracking",
        "Step 5.2: track processing",
        "Step 5.3: rotate tracks",
    ]
    assert {r["demo"] for r in rows} == {DEMOS["by_step"]}
    saved = tmp_path / "saved" / DEMOS["by_step"]
    assert sorted(p.name for p in saved.iterdir()) == ["load.mat", "sigma.mat"]
    assert list(read_mat(str(saved / "load.mat"))["data"]["markers"]) == ["gfp", "rfp"]


# ---------------------------------------------------------------------------------------------
# the comparison of the output trees
# ---------------------------------------------------------------------------------------------


def test_compare_rules(phase3):
    cmp = phase3("compare_phase3_trees")
    ex = f"{DATASET}/cond/Ex2_x_z0p5/"
    a = ex + "ch1/Analysis_update_3/"
    assert cmp.step_of("saved/demo_x/load.mat") == "load"
    assert cmp.step_of("saved/demo_x/sigma.mat") == "PSF sigmas"
    assert cmp.step_of(ex + "ch1/DS/f.tif") == "deskew"
    assert cmp.step_of(ex + "ch2/DS/f.tif") == "deskew"
    assert cmp.step_of(ex + "ch2/DS/xz_offset.mat") == "XZ offsets"
    assert cmp.step_of(a + "settings.mat") == "deskew"
    assert cmp.step_of(a + "dynRange.mat") == "deskew"
    for f in ("Detection3D.mat", "Masks/dmask_01.tif", "Detections/Detection_F_3.mat"):
        assert cmp.step_of(a + f) == "detection"
    assert cmp.step_of(a + "trackedFeatures.mat") == "tracking"
    assert cmp.step_of(a + "ProcessedTracks.mat") == "track processing"
    assert cmp.step_of(a + "Amira/Tracking/tracking_t_0001.am") == "track processing"
    assert cmp.step_of(a + "RotatedTracks.mat") == "track rotation"
    assert cmp.step_of(a + "lifetimeData.mat") == "lifetime analysis"
    assert cmp.step_of("saved/demo_x/cohorts.mat") == "intensity cohorts"
    assert cmp.step_of(ex + "ch1/Analysis/other.mat") is None
    assert cmp.step_of(ex + "ch1/raw.tif") is None


def test_compare_tolerances_are_the_rows(phase3):
    cmp = phase3("compare_phase3_trees")
    e2e, alone = cmp.tolerances(), cmp.tolerances(matlab_inputs=True)
    assert e2e["*/sigma.mat"].rel == 1e-7 and "D6" in e2e["*/sigma.mat"].reason
    assert alone["*/sigma.mat"] == e2e["*/sigma.mat"]
    # each step on MATLAB's inputs: its own rows
    assert alone["*/Detection3D.mat"].rel == 1e-9 and "D3" in alone["*/Detection3D.mat"].reason
    assert "*/trackedFeatures.mat" not in alone  # TK1-TK4: bit for bit
    assert alone["*/ProcessedTracks.mat"].rel == 1e-11
    assert alone["*/ProcessedTracks.mat"].abs == 1e-10
    # end to end: what the PSF sigmas feed (E1)
    for g in ("Detection3D.mat", "Detections/*.mat", "trackedFeatures.mat", "RotatedTracks.mat"):
        assert e2e["*/" + g].rel == 1e-8 and "E1" in e2e["*/" + g].reason
    spec = cmp.tree_spec()
    assert spec.spec_for(f"{DATASET}/c/Ex1/ch1/DS/f.tif") == spec.values
    assert spec.spec_for("x/ch1/A/RotatedTracks.mat") == e2e["*/RotatedTracks.mat"]


def test_lifetime_tolerances_are_the_rows(phase3):
    """the lifetime files: MATLAB's on MATLAB's tracks but the threshold (L4) and the cohorts'
    amplitudes (L7); end to end E1 for the others, and at least E1 for those"""
    from petakit5d.testing.compare import EXACT, Ignore

    cmp = phase3("compare_phase3_trees")
    specs, other = cmp.lifetime_specs(matlab_inputs=True)
    assert other == EXACT
    for k in ("lftRes.MaxIntensityThreshold", "MaxIntensityThreshold"):
        assert specs[k].rel == 1e-4 and "L4" in specs[k].reason
    for k in ("cohorts.A", "cohorts.Amin", "cohorts.Aplus"):
        assert (specs[k].rel, specs[k].abs) == (5e-3, 1e-9) and "L7" in specs[k].reason
    specs, other = cmp.lifetime_specs()
    assert (other.rel, other.abs) == (1e-8, 1e-8) and "E1" in other.reason
    assert specs["MaxIntensityThreshold"].rel == 1e-4 and specs["MaxIntensityThreshold"].abs == 1e-8
    assert (specs["cohorts.A"].rel, specs["cohorts.A"].abs) == (5e-3, 1e-8)
    spec = cmp.tree_spec()
    for f in ("lifetimeData.mat", "Parameters.mat"):  # not in the tree comparison: by fields
        assert isinstance(spec.spec_for(f"{DATASET}/c/Ex1/ch1/A/{f}"), Ignore)
    assert isinstance(spec.spec_for("saved/demo_x/cohorts.mat"), Ignore)


def _lifetime_values(tmp_path, name, side, **kw):
    from petakit5d.io import write_mat
    from petakit5d.testing.compare import read_mat

    path = tmp_path / f"{side}_{name}"
    write_mat(str(path), _lifetime_files(**kw)[name], version="7")
    return read_mat(str(path))


@pytest.mark.parametrize(
    ("name", "change", "alone", "e2e"),
    [
        ("lifetime.mat", {}, "equal", "equal"),
        ("lifetime.mat", {"threshold": 100.0001}, "within", "within"),  # L4
        ("lifetime.mat", {"threshold": 100.1}, "problem", "problem"),
        ("Parameters.mat", {"threshold": 100.0001}, "within", "within"),
        ("Parameters.mat", {"a2": 0.9 * (1 + 1e-9)}, "problem", "within"),  # exact; E1
        ("lifetimeData.mat", {"life": 4.0 * (1 + 1e-9)}, "problem", "within"),
        ("lifetimeData.mat", {"life": 6.0}, "problem", "problem"),
        ("cohorts.mat", {"cohort_a": 10.001}, "within", "within"),  # L7
        ("cohorts.mat", {"cohort_a": 10.1}, "problem", "problem"),
    ],
)
def test_lifetime_files_field_by_field(phase3, tmp_path, name, change, alone, e2e):
    cmp = phase3("compare_phase3_trees")
    mv = _lifetime_values(tmp_path, name, "m")
    pv = _lifetime_values(tmp_path, name, "p", **change)
    for matlab_inputs, want in ((True, alone), (False, e2e)):
        rep = cmp.lifetime_problems(mv, pv, name, *cmp.lifetime_specs(matlab_inputs))
        within = [n for n in rep.notes if cmp.WITHIN in n]
        got = "problem" if rep.problems else "within" if within else "equal"
        assert got == want, (matlab_inputs, rep.text("t"))
    if name == "cohorts.mat" and change:  # each element of the struct array by its index
        where = [p.path for p in rep.problems] + rep.notes
        assert any("files/cohorts.mat.cohorts(1).A" in w for w in where), where


def test_lifetime_files_have_matlabs_variables_and_fields(phase3, tmp_path):
    from petakit5d.io import write_mat
    from petakit5d.testing.compare import read_mat

    cmp = phase3("compare_phase3_trees")
    specs = cmp.lifetime_specs(True)
    mv = _lifetime_values(tmp_path, "Parameters.mat", "m")
    write_mat(str(tmp_path / "p.mat"), {"a": np.array([[1.0, 0.9]])}, version="7")
    rep = cmp.lifetime_problems(mv, read_mat(str(tmp_path / "p.mat")), "P.mat", *specs)
    assert [p.what for p in rep.problems] == ["names"]
    assert "MATLAB's only ['MaxIntensityThreshold']" in rep.problems[0].lines[0]
    mv = _lifetime_values(tmp_path, "lifetime.mat", "m")
    res = {"lftRes": {"nCCP": np.array([[5.0]]), "MaxIntensityThreshold": 100.0}}
    write_mat(str(tmp_path / "q.mat"), res, version="7")
    rep = cmp.lifetime_problems(mv, read_mat(str(tmp_path / "q.mat")), "L.mat", *specs)
    assert [p.what for p in rep.problems] == ["names"]  # MATLAB's order of the fields
    assert rep.problems[0].path == "files/L.mat.lftRes"


def test_amira_numbers(phase3):
    cmp = phase3("compare_phase3_trees")
    m = "# AmiraMesh\ndefine VERTEX 2\n@1\n1.5 2.25 -3\n@8\n120\n"
    assert cmp.amira_problems(m, m, 0, 0) == ([], {"n": 5, "differ": 0, "max_rel_to_1": 0.0})
    p = m.replace("2.25", "2.2500000000004")
    probs, st = cmp.amira_problems(m, p, 1e-11, 1e-10)
    assert probs == [] and st["differ"] == 1
    probs, _ = cmp.amira_problems(m, p, 1e-15, 1e-15)
    assert len(probs) == 1 and "token 7" in probs[0]
    probs, _ = cmp.amira_problems(m, m.replace("120", "121"), 0.1, 10)
    assert probs and "121" in probs[0]  # integers exact
    probs, _ = cmp.amira_problems(m, m.replace("VERTEX 2", "VERTEX 3"), 0.1, 10)
    assert probs
    assert cmp.amira_problems(m, m + "1\n", 0, 0)[0] == [
        "11 tokens in MATLAB's file, 12 in Python's"
    ]


def test_array_stats(phase3):
    cmp = phase3("compare_phase3_trees")
    r = np.linspace(1.0, 2.0, 200)
    a = r.copy()
    a[:10] *= 1 + 1e-9
    a[10] = np.nan
    s = cmp.array_stats(r, a)
    assert (s["differ"], s["nonfinite_differ"]) == (11, 1)
    assert s["max_rel"] == pytest.approx(1e-9, rel=1e-6)
    assert s["rel_p50"] == 0.0 and s["rel_p99"] == pytest.approx(1e-9, rel=1e-6)
    assert cmp.array_stats(r, r)["differ"] == 0
    assert "size" in cmp.array_stats(np.zeros(3), np.zeros(4))
    m = cmp._merge([s, cmp.array_stats(r, r)])
    assert (m["n"], m["differ"], m["arrays"]) == (400, 11, 2)


def _frame_info(n_per_frame, shift=0.0):
    rows = []
    for k, n in enumerate(n_per_frame):
        x = np.arange(1.0, n + 1.0).reshape(1, -1) + 10 * k
        rows.append(
            {
                "x": x + shift,
                "A": 100 + x,
                "hval_AD": np.zeros((1, n), dtype=bool),
                "dRange": np.array([[1.0, 2.0]]),
            }
        )
    from petakit5d.io.matfile import StructArray

    return StructArray((1, len(rows)), list(rows[0]), rows)


def _lifetime_files(threshold=100.0, life=4.0, cohort_a=10.0, a2=0.9) -> dict:
    """the lifetime analysis' files, small: {file: variables} (results folder, saved)"""
    from petakit5d.io.matfile import StructArray

    rows = [
        {
            "bounds": np.array([[1.0, 2.0]]),
            "t": np.array([[0.0, 2.0]]),
            "A": np.array([[cohort_a]]),
        },
        {"bounds": np.zeros((0, 0)), "t": np.array([[0.0, 2.0]]), "A": np.array([[cohort_a / 2]])},
    ]
    return {
        "lifetimeData.mat": {"lftData": {"lifetime_s": np.array([[2.0, life]]), "A": np.eye(2)}},
        "Parameters.mat": {"MaxIntensityThreshold": threshold, "a": np.array([[1.0, a2]])},
        "lifetime.mat": {"lftRes": {"MaxIntensityThreshold": threshold, "nCCP": np.array([[5.0]])}},
        "cohorts.mat": {"cohorts": StructArray((1, 2), ["bounds", "t", "A"], rows)},
    }


def _write_trees(
    root: Path, side: str, *, fi, ds, sigma, tracks_x=0.0, extra=False, am="1.5", lifetime=None
):
    from petakit5d.io import write_mat, write_tiff

    d = root / side
    ex = d / DATASET / "cond" / "Ex1_x_z0p5"
    (ex / "ch1" / "DS").mkdir(parents=True)
    write_tiff(ds, str(ex / "ch1" / "DS" / "f1.tif"))
    a = ex / "ch1" / "Analysis_update_3"
    (a / "Amira" / "Tracking").mkdir(parents=True)
    write_mat(str(a / "settings.mat"), {"dz": 0.5}, version="7")
    write_mat(str(a / "Detection3D.mat"), {"frameInfo": fi}, version="7.3")
    tf = {"tracksCoordAmpCG": np.array([[1.0, 2.0, 3.0]]) + tracks_x, "seqOfEvents": np.eye(2)}
    write_mat(str(a / "trackedFeatures.mat"), {"tracksFinal": tf}, version="7")
    (a / "Amira" / "Tracking" / "tracking_t_0001.am").write_text(
        f"# AmiraMesh\n@1\n{am} 2 3\n", newline="\n"
    )
    if extra:
        write_mat(str(ex / "ch1" / "Analysis_update_3" / "other.mat"), {"v": 1.0}, version="7")
    (ex / "ch1" / "raw.tif").write_bytes(b"input")
    saved = root / f"saved_{side}" / "demo_x"
    saved.mkdir(parents=True)
    write_mat(str(saved / "sigma.mat"), {"sigma_mat": sigma}, version="7")
    for name, values in (lifetime or {}).items():
        folder = a if name in ("lifetimeData.mat", "Parameters.mat") else saved
        write_mat(str(folder / name), values, version="7")
    return str(d), str(root / f"saved_{side}")


@pytest.fixture
def trees(tmp_path):
    """MATLAB's output tree of a few steps (and the folder of inputs)"""
    inputs = tmp_path / "in" / DATASET / "cond" / "Ex1_x_z0p5" / "ch1"
    inputs.mkdir(parents=True)
    (inputs / "raw.tif").write_bytes(b"input")
    ds = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
    sigma = np.array([[1.3, 1.25], [1.5, 1.5]])
    m = _write_trees(tmp_path, "m", fi=_frame_info([3, 2]), ds=ds, sigma=sigma)
    return tmp_path, m, ds, sigma


def _run(cmp, tmp_path, m, p, label="t", *more):
    out = tmp_path / f"out_{label}"
    args = [m[0], p[0], "--inputs", str(tmp_path / "in"), "--matlab-saved", m[1]]
    args += ["--python-saved", p[1], "--out", str(out), "--label", label, "--workers", "1", *more]
    rc = cmp.main(args)
    return (
        rc,
        json.loads((out / f"compare_{label}.json").read_text()),
        (out / f"compare_{label}.txt").read_text(),
    )


def test_compare_equal_trees(phase3, trees):
    cmp = phase3("compare_phase3_trees")
    tmp, m, ds, sigma = trees
    p = _write_trees(tmp, "p", fi=_frame_info([3, 2]), ds=ds, sigma=sigma)
    rc, res, text = _run(cmp, tmp, m, p)
    assert rc == 0, text
    assert list(res["steps"]) == [
        "PSF sigmas",
        "deskew",
        "detection",
        "tracking",
        "track processing",
    ]
    assert res["steps"]["deskew"]["identical"] == 2
    assert res["steps"]["detection"]["differing"] == []


def test_compare_within_the_rows(phase3, trees):
    """the sigmas within D6, the detections and tracks within E1 end to end, the Amira numbers
    within it: no problem, the numbers reported"""
    cmp = phase3("compare_phase3_trees")
    tmp, m, ds, sigma = trees
    p = _write_trees(
        tmp,
        "p",
        fi=_frame_info([3, 2], shift=5e-10),
        ds=ds,
        sigma=sigma * (1 + 5e-8),
        tracks_x=1e-12,
        am="1.5000000000001",
    )
    rc, res, text = _run(cmp, tmp, m, p, "w")
    assert rc == 0, text
    det = res["steps"]["detection"]
    assert len(det["differing"]) == 1 and det["problems"] == []
    assert "frameInfo.x: 5 of 5 float64 values differ" in text
    assert "detections: MATLAB 5, Python 5 in 2 frames" in text
    assert "Amira files: 1 of 1 differ, 1 of 3 numbers" in text
    # each step on MATLAB's inputs: the tracker's own rows allow no difference
    rc, res, text = _run(cmp, tmp, m, p, "a", "--matlab-inputs")
    assert rc == 1
    assert res["steps"]["tracking"]["problems"] and not res["steps"]["detection"]["problems"]


def test_compare_the_lifetime_files(phase3, tmp_path):
    """the lifetime step's files in the trees: compared field by field, by their rows"""
    cmp = phase3("compare_phase3_trees")
    inputs = tmp_path / "in" / DATASET / "cond" / "Ex1_x_z0p5" / "ch1"
    inputs.mkdir(parents=True)
    (inputs / "raw.tif").write_bytes(b"input")
    ds = np.zeros((2, 3, 4), dtype=np.float32)
    sigma = np.array([[1.3, 1.25], [1.5, 1.5]])
    trees = {}
    for side, kw in (("m", {}), ("p", {"threshold": 100.0001, "life": 4.0 * (1 + 1e-9)})):
        fi = _frame_info([3, 2])
        lifetime = _lifetime_files(**kw)
        trees[side] = _write_trees(tmp_path, side, fi=fi, ds=ds, sigma=sigma, lifetime=lifetime)
    rc, res, text = _run(cmp, tmp_path, trees["m"], trees["p"], "lt")
    assert rc == 0, text
    lt = res["steps"]["lifetime analysis"]
    assert (lt["matlab"], lt["python"], lt["identical"]) == (3, 3, 0)
    assert sorted(n.rsplit("/", 1)[-1] for n in lt["differing"]) == [
        "Parameters.mat",
        "lifetime.mat",
        "lifetimeData.mat",
    ]
    assert res["steps"]["intensity cohorts"]["identical"] == 1
    rc, res, text = _run(cmp, tmp_path, trees["m"], trees["p"], "lta", "--matlab-inputs")
    assert rc == 1  # on MATLAB's tracks the lifetimes are MATLAB's
    probs = res["steps"]["lifetime analysis"]["problems"]
    assert [p["path"].rsplit("/", 1)[-1] for p in probs] == ["lifetimeData.mat.lftData.lifetime_s"]


def test_compare_beyond_the_rows(phase3, trees):
    cmp = phase3("compare_phase3_trees")
    tmp, m, ds, sigma = trees
    ds2 = ds.copy()
    ds2[1, 2, 3] = np.nextafter(ds2[1, 2, 3], np.float32(100))
    p = _write_trees(
        tmp, "p", fi=_frame_info([3, 1], shift=0.0), ds=ds2, sigma=sigma * (1 + 1e-6), extra=True
    )
    rc, res, text = _run(cmp, tmp, m, p, "b")
    assert rc == 1
    assert res["steps"]["deskew"]["problems"]  # exact
    assert res["steps"]["PSF sigmas"]["problems"]  # beyond D6
    assert res["steps"]["detection"]["problems"]  # a detection less in frame 2
    assert "frames with other counts (MATLAB/Python): 1: 2 (2/1)" in text
    assert res["other"] == [
        ["Python", f"{DATASET}/cond/Ex1_x_z0p5/ch1/Analysis_update_3/other.mat"]
    ]
