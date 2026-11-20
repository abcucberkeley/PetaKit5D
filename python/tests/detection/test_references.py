"""Black-box references of the point-detection stages (task C01): MATLAB's outputs of each stage on
synthetic movies, recorded through the oracle harness, for the new implementations of C02 to C05
(docs/detection_scope.md lists the stages, their inputs and outputs, and these cases).

The cases are MATLAB-only (no Python side yet): the stage owners give them their Python calls
and comparisons (the case key is the MATLAB call and the inputs, so a case copied unchanged
keeps its reference). Each case clears the detection functions first, since some keep caches
in persistent variables across calls (detection's t threshold is cached per alpha).

Inputs are synthetic and bit-stable (``petakit5d.testing.synth``): Gaussian spots on a
background of 100 with Gaussian noise of standard deviation 5, moving along fixed tracks, in two
channels; a condition folder laid out as XR_loadConditionData3D reads it, with each frame
written twice, as the raw frame and as its 'DS' frame (no deskew: detection reads DS), and two
bead PSFs for the sigma estimate.
"""

import numpy as np

from petakit5d.core.arrays import to_matlab
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value

EXPERIMENT = "Ex1_synth_z0p5"
FRAME_MS = 2000
SIZE = (16, 56, 64)  # (z, y, x) of a frame
NOISE = 5.0
BACKGROUND = 100.0
# the functions the traced demos reach that may keep state between calls (cleared per case)
CLEAR = (
    "XR_pointSourceDetection3D XR_pointSourceDetection3D_3 XR_fitGaussians3D XR_fitGaussian3D "
    "XR_fitGaussian3D_Ac fitGaussian3D estGaussianAmplitude3D XR_locmax3d "
    "filter_bad_points_in_detection XR_estimate_overall_background_information "
    "XR_get_camera_standard_deviation_threshold get_camera_type_from_setting_file "
    "XR_runDetection3D XR_runDetection3D_frame_parallel XR_runTracking3D XR_runTrackProcessing3D "
    "trackCloseGapsKalmanSparse linkFeaturesKalmanSparse costMatLinearMotionLink2 "
    "costMatLinearMotionCloseGaps2 kalmanGainLinearMotion kalmanInitLinearMotion "
    "getAveDispEllipseAll2 estimTrackTypeParamLM2 rotateTracks3D GU_estimateSigma3D "
    "XR_loadConditionData3D XR_loadConditionData runLifetimeAnalysis3D"
)
CLEAR_CODE = "clear " + CLEAR + "; "

# tracks: (first frame, last frame, (z, y, x) at the first frame, velocity per frame, amplitude,
# frames missing). Channel 2 sees the same spots shifted by (0, 0.4, -0.3) at 0.7 times the
# amplitude.
TRACKS = (
    (0, 15, (5.0, 10.3, 12.6), (0.05, 0.31, 0.22), 260.0, ()),
    (0, 15, (8.2, 24.7, 40.4), (-0.04, 0.18, -0.35), 320.0, ()),
    (0, 15, (10.6, 40.1, 20.2), (0.0, -0.27, 0.41), 210.0, ()),
    (0, 15, (6.4, 46.5, 52.3), (0.06, -0.12, -0.2), 380.0, ()),
    (0, 9, (9.1, 16.2, 50.8), (0.02, 0.44, 0.1), 240.0, ()),
    (4, 15, (7.7, 30.4, 28.9), (-0.03, -0.2, 0.33), 300.0, ()),
    (3, 10, (11.3, 12.9, 30.1), (0.0, 0.25, 0.25), 280.0, ()),
    (5, 11, (4.6, 34.8, 8.7), (0.08, 0.1, 0.36), 230.0, ()),
    (4, 12, (8.9, 50.2, 36.4), (-0.05, -0.31, 0.05), 340.0, (8,)),
    (6, 13, (12.1, 22.6, 56.2), (0.0, 0.2, -0.28), 270.0, ()),
)
N_FRAMES = 16
SIGMA = {1: (1.6, 1.3), 2: (1.9, 1.5)}  # (z, yx) voxels of the spots in each channel
PSF_SIGMA = {1: (3.2, 1.3), 2: (3.8, 1.5)}  # (z, yx) of the bead PSFs (z before / zAniso)


def _spots(frame: int, channel: int) -> tuple[np.ndarray, np.ndarray]:
    pos, amp = [], []
    for first, last, p0, v, a, missing in TRACKS:
        if first <= frame <= last and frame not in missing:
            t = frame - first
            p = [p0[i] + v[i] * t for i in range(3)]
            if channel == 2:
                p = [p[0], p[1] + 0.4, p[2] - 0.3]
            pos.append(p)
            amp.append(a if channel == 1 else 0.7 * a)
    return np.array(pos, dtype=np.float64).reshape(-1, 3), np.array(amp, dtype=np.float64)


def _noisy(vol: np.ndarray, stream: str) -> np.ndarray:
    noise = synth.CounterRNG(2026, stream).normal(vol.shape)
    v = vol.astype(np.float64) + NOISE * noise
    t = np.trunc(v)
    v = t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)
    return np.clip(v, 0, 65535).astype(np.uint16)


def frame_volume(frame: int, channel: int) -> np.ndarray:
    """frame of a channel, (z, y, x) uint16"""
    pos, amp = _spots(frame, channel)
    vol = synth.beads(SIZE, pos, SIGMA[channel], amp, BACKGROUND)
    return _noisy(vol, f"frame{frame}/ch{channel}")


def psf_volume(channel: int) -> np.ndarray:
    vol = synth.beads((40, 48, 48), [(20.0, 24.0, 24.0)], PSF_SIGMA[channel], 1000.0, BACKGROUND)
    return _noisy(vol, f"psf/ch{channel}")


def frame_name(channel: int, frame: int) -> str:
    cam, wl = ("B", 560) if channel == 1 else ("A", 642)
    t = frame * FRAME_MS
    return (
        f"synth_Cam{cam}_ch{channel - 1}_stack{frame:04d}_{wl}nm_{t:07d}msec_"
        f"{1000000 + t:010d}msecAbs.tif"
    )


def make_condition(folder) -> None:
    """the synthetic condition: one experiment, two channels, raw and DS frames, the settings
    file's camera line, two bead PSFs"""
    from pathlib import Path

    root = Path(folder)
    ex = root / EXPERIMENT
    for c in (1, 2):
        for sub in (ex / f"ch{c}", ex / f"ch{c}" / "DS"):
            sub.mkdir(parents=True, exist_ok=True)
        for k in range(N_FRAMES):
            vol = frame_volume(k, c)
            synth.write_tiff(ex / f"ch{c}" / frame_name(c, k), vol)
            synth.write_tiff(ex / f"ch{c}" / "DS" / frame_name(c, k), vol)
        synth.write_tiff(root / f"psf_ch{c}.tif", psf_volume(c))
    with open(ex / "synth_Settings.txt", "w", newline="\n") as f:
        f.write('Camera type = "Orca4.0"\n')


def ds_frame(frame: int, channel: int) -> np.ndarray:
    """a DS frame as detection reads it: double, MATLAB's (y, x, z)"""
    return to_matlab(frame_volume(frame, channel).astype(np.float64))


# detection's options as the traced runs pass them (XR_runDetection3D_frame_parallel with
# 'lowSNR' and background correction, GU_estimateSigma3D, PointDetectionTest)
LOWSNR_OPTS = (
    "'Alpha', 0.05, 'Mask', [], 'RemoveRedundant', true, 'RedundancyRadius', 0.5, "
    "'RefineMaskLoG', false, 'WindowSize', [], 'Mode', 'xyzAc', 'FitMixtures', false, "
    "'MaxMixtures', 2, 'BackgroundCorrection', true, 'FilterByResidualSigma', true, "
    "'ResidualSigmaThresh', 3"
)
SIGMA_EST_OPTS = (
    "'Mode', 'xyzAsrc', 'Alpha', 0.05, 'Mask', [], 'RemoveRedundant', true, "
    "'RefineMaskLoG', false, 'WindowSize', [], 'FitGaussianMethod', 'original', "
    "'BackgroundCorrection', false"
)
SIGMA_DS = {1: [1.3, 1.6], 2: [1.5, 1.9]}  # [sigma_xy, sigma_z] for detection on the DS frames

LOAD = (
    "rt = in{1}; ex = [rt '" + EXPERIMENT + "' filesep]; "
    "data = XR_loadConditionData3D(ex, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
)
SIGMAS = (
    "[sxy1, sz1] = GU_estimateSigma3D([rt 'psf_ch1.tif'], []); "
    "[sxy2, sz2] = GU_estimateSigma3D([rt 'psf_ch2.tif'], []); "
    "sigma = [sxy1, sz1 ./ data(1).zAniso; sxy2, sz2 ./ data(1).zAniso]; "
)
APATH = "apath = {[data(1).source 'Analysis' filesep]}; "


def track_code(parallel: bool) -> str:
    """tracking, track processing and track rotation as the demos call them; ``parallel``:
    with XR_runDetTrack3d's 'ParallelMethod'"""
    pm = ", 'ParallelMethod', 'lock_distributed'" if parallel else ""
    return (
        "XR_runTracking3D(data, loadTrackSettings('Radius', [3, 6], 'MaxGapLength', 2), "
        "'FileName', 'trackedFeatures.mat', 'Overwrite', true, 'ResultsPath', apath" + pm + "); "
        "XR_runTrackProcessing3D(data, 'Overwrite', true, 'TrackerOutput', "
        "'trackedFeatures.mat', 'FileName', 'ProcessedTracks.mat', 'Buffer', [3, 3], "
        "'BufferAll', false, 'FitMixtures', false, 'WindowSize', [], 'ResultsPath', apath"
        + pm
        + "); rotateTracks3D(data, 'Overwrite', true, 'Crop', false, 'ResultsPath', apath); "
    )


def _cases():
    cond = {"cond": Tree(make_condition, "cond")}
    cs = [
        Case(
            "load_condition",
            code=CLEAR_CODE + LOAD + "out = {data};",
            args=[InPath("cond")],
            inputs=cond,
        ),
        Case(
            "estimate_sigma",
            code=CLEAR_CODE
            + "rt = in{1}; [a, b, c] = GU_estimateSigma3D([rt 'psf_ch1.tif'], []); "
            + "[d, e, f] = GU_estimateSigma3D([rt 'psf_ch2.tif'], []); out = {a, b, c, d, e, f};",
            args=[InPath("cond")],
            nargout=6,
            inputs=cond,
        ),
        Case(
            "camera_and_background",
            code=CLEAR_CODE
            + LOAD
            + "ct = get_camera_type_from_setting_file(data); "
            + "th = XR_get_camera_standard_deviation_threshold(data); "
            + "bi = XR_estimate_overall_background_information(double(readtiff("
            + "data.framePathsDS{1}{1})), 'randseed', 1); out = {ct, th, bi};",
            args=[InPath("cond")],
            nargout=3,
            inputs=cond,
        ),
        Case(
            "track_settings",
            code=CLEAR_CODE
            + "a = loadTrackSettings('Radius', [3, 6], 'MaxGapLength', 2); "
            + "b = loadTrackSettings(); out = {a, b};",
            nargout=2,
        ),
    ]
    for c in (1, 2):
        for k in (0, 8):
            cs.append(
                Case(
                    f"psd_lowsnr.ch{c}.frame{k}",
                    code=CLEAR_CODE
                    + "vol = in{1}; bi = XR_estimate_overall_background_information(vol, "
                    + "'randseed', 1); [ps, mask] = XR_pointSourceDetection3D_3(vol, in{2}, "
                    + LOWSNR_OPTS
                    + ", 'BackgroundInfo', [bi, 3]); out = {bi, ps, mask};",
                    args=[In("vol"), SIGMA_DS[c]],
                    nargout=3,
                    inputs={"vol": Value(lambda c=c, k=k: ds_frame(k, c))},
                )
            )
    cs.append(
        Case(
            "psd_threshold_cache",
            code=CLEAR_CODE
            + "v1 = in{1}; v2 = in{2}; b1 = XR_estimate_overall_background_information(v1, "
            + "'randseed', 1); b2 = XR_estimate_overall_background_information(v2, 'randseed', "
            + "1); a = XR_pointSourceDetection3D_3(v1, in{3}, "
            + LOWSNR_OPTS
            + ", 'BackgroundInfo', [b1, 3]); b = XR_pointSourceDetection3D_3(v2, in{4}, "
            + LOWSNR_OPTS
            + ", 'BackgroundInfo', [b2, 3]); "
            + CLEAR_CODE
            + "c = XR_pointSourceDetection3D_3(v2, in{4}, "
            + LOWSNR_OPTS
            + ", 'BackgroundInfo', [b2, 3]); out = {a, b, c, isequaln(b, c)};",
            args=[In("v1"), In("v2"), SIGMA_DS[1], SIGMA_DS[2]],
            nargout=4,
            inputs={"v1": Value(lambda: ds_frame(0, 1)), "v2": Value(lambda: ds_frame(0, 2))},
        )
    )
    for c in (1, 2):
        cs.append(
            Case(
                f"psd_speedup_form.ch{c}",
                code=CLEAR_CODE
                + "[ps, mask] = XR_pointSourceDetection3D(in{1}, in{2}, "
                + "'BackgroundCorrection', false); out = {ps, mask};",
                args=[In("vol"), SIGMA_DS[c]],
                nargout=2,
                inputs={"vol": Value(lambda c=c: ds_frame(4, c))},
            )
        )
        cs.append(
            Case(
                f"psd_sigma_estimate_form.ch{c}",
                code=CLEAR_CODE
                + "[ps, mask] = XR_pointSourceDetection3D(in{1}, in{2}, "
                + SIGMA_EST_OPTS
                + "); out = {ps, mask};",
                args=[In("psf"), [1.2, 2.5]],
                nargout=2,
                inputs={"psf": Value(lambda c=c: to_matlab(psf_volume(c).astype(np.float64)))},
            )
        )
    # the stages through the files they write, as the demos run them
    by_step = (
        CLEAR_CODE
        + LOAD
        + SIGMAS
        + APATH
        + "XR_runDetection3D(data, 'Sigma', sigma, 'Overwrite', true, 'WindowSize', [], "
        + "'ResultsPath', apath, 'Mode', 'xyzAc', 'FitMixtures', false, 'MaxMixtures', 1, "
        + "'DetectionMethod', 'lowSNR', 'BackgroundCorrection', true); "
        + track_code(False)
        + "out = {data, sigma};"
    )
    integrated = (
        CLEAR_CODE
        + LOAD
        + SIGMAS
        + APATH
        + "XR_runDetection3D_frame_parallel(data, 'Sigma', sigma, 'Overwrite', false, "
        + "'WindowSize', [], 'ResultsPath', apath, 'Mode', 'xyzAc', 'FitMixtures', false, "
        + "'MaxMixtures', 2, 'DetectionMethod', 'lowSNR', 'BackgroundCorrection', true); "
        + track_code(True)
        + "out = {data, sigma};"
    )
    speedup = (
        CLEAR_CODE
        + LOAD
        + APATH
        + "XR_runDetection3D(data, 'Sigma', in{2}, 'DetectionMethod', 'speedup', 'WriteAmira', "
        + "false, 'Overwrite', true, 'ResultsPath', apath); "
        + "XR_runTracking3D(data, loadTrackSettings('Radius', [3, 6], 'MaxGapLength', 2), "
        + "'FileName', 'trackedFeatures.mat', 'Overwrite', true, 'ResultsPath', apath); "
        + "out = {data};"
    )
    two_sigmas = np.array([[1.3, 1.6], [1.5, 1.9]])
    for name, code, args, nout in (
        ("chain.by_step", by_step, [InPath("cond")], 2),
        ("chain.integrated", integrated, [InPath("cond")], 2),
        ("chain.speedup_no_background", speedup, [InPath("cond"), two_sigmas], 1),
    ):
        cs.append(Case(name, code=code, args=args, nargout=nout, inputs=cond, marks=("slow",)))
    # the multiple-experiments demo: XR_runDetTrack3d as it calls it (with its unknown
    # 'zoffsetCorrection'), then the lifetime analysis and the intensity cohorts
    multiple = (
        CLEAR_CODE
        + "rt = in{1}; data = XR_loadConditionData3D(rt, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
        + "XR_runDetTrack3d(data, 'Sigma', in{2}, 'DetectionMethod', 'lowSNR', 'aname', "
        + "'Analysis', 'FitMixtures', false, 'BackgroundCorrection', true, "
        + "'zoffsetCorrection', true); "
        + "lftRes = runLifetimeAnalysis3D(data, 'AnalysisPath', 'Analysis', 'overwrite', true); "
        + "cohorts = plotIntensityCohorts(data, 'AnalysisPath', 'Analysis'); "
        + "out = {data, lftRes, cohorts};"
    )
    # on this small condition MATLAB's lifetime analysis stops with ecdf's error on an empty
    # category; the lifetime reference is the demo data's recorded run (docs/detection_scope.md)
    cs.append(
        Case(
            "multiple.lifetime_empty_category",
            code=multiple,
            args=[InPath("lcond"), two_sigmas],
            nargout=3,
            inputs={"lcond": Tree(make_lifetime_condition, "lcond")},
            expect_error="stats:ecdf:NotEnoughData",
            marks=("slow",),
        )
    )
    return cs


L_FRAMES = 30
L_SIZE = (14, 48, 56)
L_EXPERIMENTS = ("Ex1_synthL_z0p5", "Ex2_synthL_z0p5", "Ex3_synthL_z0p5")


def lifetime_tracks(experiment: int) -> list[tuple]:
    """22 tracks per experiment: lifetimes of 3 to 17 frames, starting after frame 4 and ending
    4 frames before the end (the buffers of track processing), positions inside the volume,
    slow motion; from a bit-stable stream"""
    u = synth.CounterRNG(2026, f"lifetime/ex{experiment}").uniform((22, 9))
    tracks = []
    for r in u:
        first = 4 + int(r[0] * (L_FRAMES - 12))
        last = min(L_FRAMES - 5, first + 2 + int(r[1] * 15))
        p0 = (
            2.5 + r[2] * (L_SIZE[0] - 5),
            4.0 + r[3] * (L_SIZE[1] - 8),
            4.0 + r[4] * (L_SIZE[2] - 8),
        )
        v = ((r[5] - 0.5) * 0.1, (r[6] - 0.5) * 0.6, (r[7] - 0.5) * 0.6)
        tracks.append((first, last, p0, v, 180.0 + 260.0 * r[8]))
    return tracks


def lifetime_frame(experiment: int, frame: int, channel: int) -> np.ndarray:
    pos, amp = [], []
    for first, last, p0, v, a in lifetime_tracks(experiment):
        if first <= frame <= last:
            t = frame - first
            p = [p0[i] + v[i] * t for i in range(3)]
            if channel == 2:
                p = [p[0], p[1] + 0.4, p[2] - 0.3]
            pos.append(p)
            amp.append(a if channel == 1 else 0.7 * a)
    pos = np.array(pos, dtype=np.float64).reshape(-1, 3)
    vol = synth.beads(L_SIZE, pos, SIGMA[channel], np.array(amp), BACKGROUND)
    return _noisy(vol, f"lifetime/ex{experiment}/frame{frame}/ch{channel}")


def make_lifetime_condition(folder) -> None:
    """three synthetic experiments of 30 frames, two channels, raw and DS frames"""
    from pathlib import Path

    root = Path(folder)
    for e, name in enumerate(L_EXPERIMENTS, start=1):
        ex = root / name
        for c in (1, 2):
            for sub in (ex / f"ch{c}", ex / f"ch{c}" / "DS"):
                sub.mkdir(parents=True, exist_ok=True)
            for k in range(L_FRAMES):
                vol = lifetime_frame(e, k, c)
                synth.write_tiff(ex / f"ch{c}" / frame_name(c, k), vol)
                synth.write_tiff(ex / f"ch{c}" / "DS" / frame_name(c, k), vol)
        with open(ex / "synthL_Settings.txt", "w", newline="\n") as f:
            f.write('Camera type = "Orca4.0"\n')


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)
