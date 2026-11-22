"""The detection runners without a camera threshold (task C09; docs/parity/detection.md D13).

The runners give the lowSNR detector a residual threshold made from each frame's noise and the
camera threshold that ``XR_get_camera_standard_deviation_threshold`` finds in the settings file
above the condition (``run.residual_threshold``). The cases run both runners on one condition
per case: no settings file, a settings file without a camera line, one naming a camera the
search does not know, and a camera threshold of NaN or 0 (the threshold function stubbed in the
case's MATLAB code, the port's ``residual_sigma_threshold``); ``test_detect.py``'s
``runner_threshold`` cases have a known camera.

The condition's frames are those of ``test_detect.py``'s threshold condition (task C02d): three
frames whose noise lies in each range of the runners' rule, each with three bright spots and a
weak one whose residual sigma lies between the rule's level and another rule's.

MATLAB, measured (probe d13; docs/parity/detection.md D13):

- no settings file, or one without a camera line: the camera search stops with an error
  (``MATLAB:subsdeldimmismatch``, ``MATLAB:TooManyOutputsDueToMissingBraces``), which
  ``XR_runDetection3D`` catches and goes on with each frame's noise as the threshold; the
  frame-parallel runner stops with it before any frame (``matlab_frame_parallel.*``);
- an unknown camera: the threshold NaN; ``XR_runDetection3D`` takes the noise, the
  frame-parallel runner passes NaN and finds nothing;
- a threshold of NaN or 0 (stubbed): NaN as an unknown camera; 0 reaches the detector in both
  runners, which keeps every candidate.

The port's frame-parallel runner does what ``XR_runDetection3D`` does in all of them
(``camera.*.frame_parallel``: its files are MATLAB's ``XR_runDetection3D``'s), MATLAB's bug
fixed.
"""

from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.dtypes import NAN
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.compare import tol
from petakit5d.testing.compare import tree as tree_spec

TOL_DETECT = tol(
    reason="parity D3: the detector's fits (3-D filter order, covariance)", rel=1e-9, abs=1e-9
)

EXPERIMENT = "Ex1_cam_z0p5"
SIZE = (16, 56, 64)  # (z, y, x)
SIGMA = np.array([[1.3, 1.6], [1.3, 1.6]])  # [sigma_xy sigma_z] of both channels
BRIGHT = [((5.2, 12.4, 14.6), 220.0), ((9.6, 40.3, 50.2), 260.0), ((11.1, 20.7, 47.4), 190.0)]
# per frame: the noise, the weak spot's patch noise, its (z, y, x) and amplitude
FRAMES = (
    (1.6, 2.7, (8.3, 36.6, 22.4), 9.0),
    (3.5, 3.1, (8.4, 36.2, 22.7), 9.0),
    (5.5, 4.0, (7.6, 35.8, 23.1), 12.0),
)
FRAME_MS = 2000
# the settings file's text per camera case (None: no settings file above the condition); the
# cases use 'camera' (for the stubs), 'none', 'nocamera' and 'unknown', the measurement also the
# others (EMCCD: 58.5 to every frame's detector, no detections; an empty name as an unknown one)
SETTINGS = {
    "camera": 'Camera type = "Orca4.0"\n',  # sCMOS: the camera threshold 3
    "emccd": 'Camera type = "Andor iXon"\n',  # EMCCD: 65
    "none": None,
    "nocamera": "Laser = 560\nExposure = 40\n",
    "unknown": 'Camera type = "pco.edge 4.2"\n',
    "empty": 'Camera type = ""\n',
}


def frame_name(channel: int, frame: int) -> str:
    cam, wl = ("B", 560) if channel == 1 else ("A", 642)
    t = frame * FRAME_MS
    return (
        f"cam_Cam{cam}_ch{channel - 1}_stack{frame:04d}_{wl}nm_{t:07d}msec_"
        f"{1000000 + t:010d}msecAbs.tif"
    )


def _counts(v: np.ndarray) -> np.ndarray:
    t = np.trunc(v)
    v = t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)
    return np.clip(v, 0, 65535).astype(np.uint16)


def cam_frame(k: int) -> np.ndarray:
    """frame k (0-based) of the condition, (z, y, x) uint16 (test_detect.py's rt_frame)"""
    noise, pn, wpos, wamp = FRAMES[k]
    sig = synth.beads(SIZE, [p for p, _ in BRIGHT], (1.6, 1.3), [a for _, a in BRIGHT], 0.0)
    weak = synth.beads(SIZE, [wpos], (1.6, 1.3), [wamp], 0.0)
    sig = sig.astype(np.float64)
    std = np.full(SIZE, noise)
    z, y, x = (int(round(t)) for t in wpos)
    std[max(z - 5, 0) : z + 6, max(y - 5, 0) : y + 6, max(x - 5, 0) : x + 6] = pn
    std = np.sqrt(std**2 + sig)
    r = synth.CounterRNG(2026, "C02d.rt2." + ("low", "mid", "high")[k]).normal(SIZE)
    return _counts(100.0 + sig + weak.astype(np.float64) + std * r)


def make_condition(folder, camera: str) -> None:
    """the condition (two channels of the same frames, raw and DS) with the settings file of a
    camera case beside its channel folders"""
    ex = Path(folder) / EXPERIMENT
    for c in (1, 2):
        for sub in (ex / f"ch{c}", ex / f"ch{c}" / "DS"):
            sub.mkdir(parents=True, exist_ok=True)
        for k in range(len(FRAMES)):
            vol = cam_frame(k)
            synth.write_tiff(ex / f"ch{c}" / frame_name(c, k), vol)
            synth.write_tiff(ex / f"ch{c}" / "DS" / frame_name(c, k), vol)
    text = SETTINGS[camera]
    if text is not None:
        with open(ex / "cam_Settings.txt", "w", newline="\n") as f:
            f.write(text)


def make_case(camera: str):
    def make(folder) -> None:
        make_condition(folder, camera)

    make.__name__ = f"make_condition_{camera}"
    return make


# the functions that may keep state between calls (cleared per case, as test_references.py's)
CLEAR_CODE = (
    "clear XR_pointSourceDetection3D XR_pointSourceDetection3D_3 XR_fitGaussians3D "
    "XR_fitGaussian3D XR_fitGaussian3D_Ac fitGaussian3D estGaussianAmplitude3D XR_locmax3d "
    "filter_bad_points_in_detection XR_estimate_overall_background_information "
    "XR_get_camera_standard_deviation_threshold get_camera_type_from_setting_file "
    "XR_runDetection3D XR_runDetection3D_frame_parallel XR_loadConditionData3D "
    "XR_loadConditionData; "
)
LOAD = (
    "rt = in{1}; ex = [rt '" + EXPERIMENT + "' filesep]; "
    "data = XR_loadConditionData3D(ex, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "apath = {[data(1).source 'Analysis' filesep]}; "
)
# XR_runDetection3D as the step-by-step demo calls it, and both runners with the options
# XR_runDetTrack3d passes the frame-parallel one (the integrated demos')
RUN = (
    "XR_runDetection3D(data, 'Sigma', in{2}, 'Overwrite', true, 'WindowSize', [], "
    "'ResultsPath', apath, 'Mode', 'xyzAc', 'FitMixtures', false, 'MaxMixtures', 1, "
    "'DetectionMethod', 'lowSNR', 'BackgroundCorrection', true, 'WriteAmira', false); "
)
_DT_OPTIONS = (
    "'Sigma', in{2}, 'Overwrite', false, 'WindowSize', [], 'ResultsPath', apath, "
    "'Mode', 'xyzAc', 'FitMixtures', false, 'MaxMixtures', 2, 'DetectionMethod', 'lowSNR', "
    "'BackgroundCorrection', true"
)
RUN_AS_DT = "XR_runDetection3D(data, " + _DT_OPTIONS + "); "
FRAME_PARALLEL = "XR_runDetection3D_frame_parallel(data, " + _DT_OPTIONS + "); "
COUNTS = "s = load([apath{1} 'Detection3D.mat']); out = {arrayfun(@(e) size(e.x, 2), s.frameInfo)};"


def stubbed(value: str, call: str) -> str:
    """a call with XR_get_camera_standard_deviation_threshold stubbed to give value (MATLAB
    code), in a folder of its own first on the path, taken off again whatever the call does"""
    return (
        "stub = tempname; mkdir(stub); "
        "fid = fopen(fullfile(stub, 'XR_get_camera_standard_deviation_threshold.m'), 'w'); "
        "fprintf(fid, 'function [t] = XR_get_camera_standard_deviation_threshold(data, "
        "varargin)\\nt = " + value + ";\\nend\\n'); fclose(fid); addpath(stub); "
        "clear XR_get_camera_standard_deviation_threshold; "
        "try, " + call + "catch e, rmpath(stub); rethrow(e); end; rmpath(stub); "
        "clear XR_get_camera_standard_deviation_threshold; "
    )


def _data(rt):
    import os

    from petakit5d.detection import condition as cond

    data = cond.load_condition_data(
        os.path.join(rt, EXPERIMENT) + os.sep, ["ch1", "ch2"], ["gfp", "rfp"]
    )
    return data, [data[0].source + "Analysis" + os.sep]


def _counts_of(entries) -> np.ndarray:
    """the detections per frame of a runner's entries (MATLAB's size(frameInfo(k).x, 2))"""
    return np.array([float(np.asarray(e["x"]).size // 2) for e in entries])


def _py_run(threshold=None):
    def run_it(rt, sigma):
        from petakit5d.detection import detect, run

        detect.clear_threshold_cache()
        data, rp = _data(rt)
        (entries,) = run.run_detection(
            data,
            sigma=np.asarray(sigma),
            overwrite=True,
            window_size=None,
            results_path=rp,
            mode="xyzAc",
            fit_mixtures=False,
            max_mixtures=1,
            detection_method="lowSNR",
            background_correction=True,
            write_amira=False,
            residual_sigma_threshold=threshold,
        )
        return _counts_of(entries)

    return run_it


def _py_frame_parallel(threshold=None):
    def run_it(rt, sigma):
        from petakit5d.detection import detect, run

        detect.clear_threshold_cache()
        data, rp = _data(rt)
        (entries,) = run.run_detection_frame_parallel(
            data,
            sigma=np.asarray(sigma),
            overwrite=False,
            window_size=None,
            results_path=rp,
            mode="xyzAc",
            fit_mixtures=False,
            max_mixtures=2,
            detection_method="lowSNR",
            background_correction=True,
            residual_sigma_threshold=threshold,
        )
        return _counts_of(entries)

    return run_it


def _py_camera_search(rt, sigma):
    """MATLAB's frame-parallel runner up to its error: the camera search (``camera_threshold``,
    which raises MATLAB's error with its identifier)"""
    import os

    from petakit5d.detection import run

    if os.name == "nt":
        # on Windows the walk ends at the drive's root, which MATLAB lists: no settings file
        # there gives MATLAB:UndefinedFunction (deduced; the references are made on Linux)
        pytest.skip("MATLAB's error without a settings file is another one on Windows")
    data, _ = _data(rt)
    run.camera_threshold(data[0].source)


def _py_nan_threshold(rt, sigma):
    """MATLAB's frame-parallel runner with a NaN camera threshold, measured: each frame's
    detector gets the threshold NaN and [bg noise NaN] as its BackgroundInfo; the port's
    detector given that, per frame, its detections' count"""
    from petakit5d.detection import detect, run
    from petakit5d.io import read_tiff

    detect.clear_threshold_cache()
    data, _ = _data(rt)
    sig = run.runner_sigma(np.asarray(sigma))
    counts = []
    for k, path in enumerate(data[0].frame_paths_ds[0], start=1):
        vol = run.frame_volume(read_tiff(path))
        bi = detect.estimate_background(vol, k)
        ps, _ = detect.point_source_detection(
            vol,
            sig[0],
            background_correction=True,
            background_info=np.r_[bi, NAN],
            residual_sigma_threshold=NAN,
        )
        counts.append(float(np.asarray(ps["x"]).size))
    return np.array(counts)


def _cases():
    tree_mat = tree_spec(per_file={"*.mat": TOL_DETECT})
    # the frame-parallel runner's own per-frame files, which XR_runDetection3D does not write
    tree_fp_as_run = tree_spec(
        per_file={"*.mat": TOL_DETECT}, ignore=[r"Detections/Detection_F_\d+\.mat$"]
    )
    cs = []

    def add(name, camera, code, python, **kw):
        cs.append(
            Case(
                name,
                code=CLEAR_CODE + LOAD + code,
                args=[InPath("cond"), SIGMA],
                inputs={"cond": Tree(make_case(camera), "cond")},
                python=python,
                **kw,
            )
        )

    # (case, the condition's camera case, the MATLAB code around a call, the port's threshold)
    variants = (
        ("none", "none", lambda c: c, None),
        ("nocamera", "nocamera", lambda c: c, None),
        ("unknown", "unknown", lambda c: c, None),
        ("stub_nan", "camera", lambda c: stubbed("NaN", c), NAN),
        ("stub_zero", "camera", lambda c: stubbed("0", c), 0.0),
    )
    # each case twice: its files (the fits within D3: decided with the store) and its
    # detections per frame alone (decided by the manifest too, in CI). The frame_parallel calls
    # (the frame-parallel runner, and XR_runDetection3D with its options) write Amira files,
    # whose writer's warning (without identifier: third-party) has its text withheld
    for name, camera, wrap, thr in variants:
        # XR_runDetection3D: each frame's noise without a camera threshold, 0 as given
        run_code = wrap(RUN) + COUNTS
        # both runners pass 0 on (the frame-parallel runner's files, per-frame ones too); else
        # the port's frame-parallel runner gives XR_runDetection3D's files (D13)
        fp_code = wrap(FRAME_PARALLEL if name == "stub_zero" else RUN_AS_DT) + COUNTS
        fp_tree = tree_mat if name == "stub_zero" else tree_fp_as_run
        add(f"camera.{name}.run_detection", camera, run_code, _py_run(thr), tree=tree_mat)
        add(
            f"camera.{name}.frame_parallel",
            camera,
            fp_code,
            _py_frame_parallel(thr),
            tree=fp_tree,
            warning_text="withheld",
        )
        # (the same MATLAB calls; without the files their keys and references are their own)
        add(f"counts.{name}.run_detection", camera, run_code, _py_run(thr), files=False)
        add(
            f"counts.{name}.frame_parallel",
            camera,
            fp_code,
            _py_frame_parallel(thr),
            files=False,
            warning_text="withheld",
        )
    # MATLAB's frame-parallel runner itself (D13): it stops with the camera search's errors,
    # and with NaN as the threshold its detector finds nothing in any frame
    for name, ident in (
        ("none", "MATLAB:subsdeldimmismatch"),
        ("nocamera", "MATLAB:TooManyOutputsDueToMissingBraces"),
    ):
        add(
            f"matlab_frame_parallel.{name}",
            name,
            FRAME_PARALLEL + "out = {};",
            _py_camera_search,
            nargout=0,
            files=False,
            expect_error=ident,
            error_text="withheld",
        )
    for name, camera, wrap in (
        ("unknown", "unknown", lambda c: c),
        ("stub_nan", "camera", lambda c: stubbed("NaN", c)),
    ):
        add(
            f"matlab_frame_parallel.{name}",
            camera,
            wrap(FRAME_PARALLEL) + COUNTS,
            _py_nan_threshold,
            files=False,
            warning_text="withheld",  # the Amira writer's, as above
        )
    return cs


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_a_camera_threshold_of_0_reaches_the_detector():
    """a given threshold of 0 takes the camera threshold's place: 0 for every frame's noise, as
    MATLAB's runners pass on a camera threshold of 0 (probe d13; NaN: the noise)"""
    from petakit5d.detection import run

    for noise in (0.5, 1.6073164540857801, 3.5492398528252935, 5.5508490987101542, 90.0):
        assert run.residual_threshold(noise, 0.0) == 0.0
        assert run.residual_threshold(noise, NAN) == noise


@pytest.mark.parametrize("runner", ["run_detection", "run_detection_frame_parallel"])
def test_amira_detection_files_asked_for_give_the_ports_warning(tmp_path, runner):
    """write_amira: no Amira file of the detections is written (MATLAB's runners call an Amira
    writer that PetaKit5D does not include, warn without an identifier and leave the folder
    empty); the port leaves the same empty folder and warns once per movie, in its own words
    with its own identifier; nothing without write_amira"""
    import os
    import warnings

    from petakit5d.detection import detect, run

    make_condition(tmp_path, "camera")
    data, _ = _data(str(tmp_path))
    for write in (True, False):
        rp = [str(tmp_path / f"{write}_{i}") + os.sep for i in (1, 2)]
        detect.clear_threshold_cache()
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always")
            getattr(run, runner)(
                [data[0], data[0]],
                sigma=SIGMA,
                results_path=rp,
                detection_method="lowSNR",
                background_correction=True,
                write_amira=write,
            )
        ours = [w.message for w in wl if issubclass(w.category, run.AmiraNotWrittenWarning)]
        assert len(ours) == (2 if write else 0)
        for w, p in zip(ours, rp, strict=False):
            assert w.identifier == "petakit5d:detection:amiraNotWritten"
            assert os.path.join(p, "Amira", "Detection") in str(w)
        for p in rp:
            folder = os.path.join(p, "Amira", "Detection")
            assert os.path.isdir(folder) == write
            assert not write or not os.listdir(folder)
