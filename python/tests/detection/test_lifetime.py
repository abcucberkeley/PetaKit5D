"""Tests of petakit5d.detection.lifetime (task C05) and of the MathWorks functions it needs, against
MATLAB R2024b (oracle cases): ``runLifetimeAnalysis3D`` and ``plotIntensityCohorts`` on synthetic
conditions, the lifetime analysis' functions as black boxes (``adtest1``, ``binterp``,
``fitGaussianModeToCDF``, ``scaleEDFs``, ``getVisitorIndex``), and matlabcompat's ``hist``,
``normcdf``, ``erf``, ``erfc`` and ``lsqnonlin``.

The synthetic conditions are folders as XR_loadConditionData3D reads them (experiments of 50 or
60 frames at 2 s in two channels, the frames tiny since the lifetime analysis does not read them),
each experiment with a ``ch1/Analysis/ProcessedTracks.mat`` written here: tracks with the fields
and shapes of MATLAB's track processing, drawn from a bit-stable stream (``testing.synth``) and
saved as MAT v6 with a fixed header, so the files are the same bytes everywhere. Their lifetimes
fill the analysis' six threshold cohorts; their amplitudes rise with the lifetime (``early``: from
the first frame, so the cohorts differ from 3 frames on; ``late``: from frame 6; ``flat``: not at
all), a few short tracks start bright (visitors), and the categories, the second channel's
significance and the gaps are mixed; the tracks are in random order (MATLAB sorts them).

Tolerances: docs/parity/detection.md, rows L1 to L7 (approved by the user on 2026-11-20). The
fits are MATLAB's lsqnonlin (trust-region-reflective) reproduced from its measured rules, not its
code: where the port's path through a fit is MATLAB's, the results agree to about 1e-9; where it
is not (some cases here, named by their tolerance), they differ by more than rounding.

The errors of adtest1 and binterp have no identifier and GPL messages: their cases withhold the
messages from the manifest and compare the (empty) identifier only (``error_text``, decided by
the user on 2026-11-20). ``lifetime.demo`` runs on MATLAB's tracks of the multiple-experiments
demo (``$PETAKIT5D_DEMO_DATA``: the demo dataset's frames for the condition loader, and
``PetaKit5D_point_detection_demo_tracks``, the demo run's ``ProcessedTracks.mat`` files).
"""

import math
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.detection import lifetime as lt
from petakit5d.io.matfile import StructArray, read_mat, write_mat
from petakit5d.matlabcompat import elementary as el
from petakit5d.matlabcompat import hist as mhist
from petakit5d.matlabcompat import optim, special
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, Data, In, InPath, Tree, Value
from petakit5d.testing.compare import ignore, tol

NAN = dtypes.NAN
INF = math.inf

# -- tolerances (docs/parity/detection.md, proposed) --------------------------------------------

TOL_LSQ = tol(reason="parity/detection.md L1: lsqnonlin's rounding", rel=1e-6, abs=1e-12)
TOL_SCALE = tol(reason="parity/detection.md L2: scaleEDFs' fits", rel=1e-6, abs=1e-9)
TOL_MODE = tol(reason="parity/detection.md L3: fitGaussianModeToCDF's path", rel=0.03)
TOL_THRESHOLD = tol(reason="parity/detection.md L4: the threshold's fit", rel=1e-4)
TOL_AD = tol(reason="parity/detection.md L5: adtest1's logarithms", rel=1e-13)
TOL_BINTERP = tol(reason="parity/detection.md L6: binterp", rel_max=2e-13, abs=1e-13)
TOL_SCALE_PATH = tol(
    reason="parity/detection.md L2: scaleEDFs' fits, another path", rel=1e-3, abs=5e-3
)
TOL_COHORT = tol(reason="parity/detection.md L7: scaled cohorts", rel=5e-3, abs=1e-9)
TOL_COHORT_ONE = tol(
    reason="parity/detection.md L7: one movie's cohorts, binterp",
    rel_max=2e-13,
    abs=1e-12,
)

# -- synthetic conditions -----------------------------------------------------------------------

FRAME_MS = 2000
FR = FRAME_MS / 1000
TINY = (2, 4, 4)  # (z, y, x) of the frames: the lifetime analysis never reads them
# frames of the tracks in each threshold cohort at 2 s (the last one up to the movie's length)
COHORT_FRAMES = ((3, 5), (6, 7), (8, 10), (11, 20), (21, 30), (31, 48))
CATEGORIES = (1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 5, 3, 4, 6, 7, 8)
TRACK_FIELDS = tuple(
    "t f x y z A c x_pstd y_pstd z_pstd A_pstd c_pstd sigma_r SE_sigma_r pval_Ar isPSF "
    "tracksFeatIndxCG gapVect gapStatus gapIdx seqOfEvents nSeg visibility lifetime_s start end "
    "startBuffer endBuffer MotionAnalysis RSS hval_Ar hval_AD catIdx hasMultiplePeaks "
    "significantMaster significantVsBackground significantSlave".split()
)  # MATLAB's track processing's order
MAT_HEADER = b"MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Thu Jan  1 00:00:00 2026"
CONDITIONS = {
    # name: (experiments, amplitude design, tracks per experiment, frames, all six cohorts)
    "three_early": (3, "early", 120, 60, True),
    "two_early": (2, "early", 150, 60, True),
    "one_late": (1, "late", 220, 50, True),
    "flat": (1, "flat", 200, 60, True),
    "no_long_tracks": (1, "early", 80, 60, False),
}


def _row(v) -> np.ndarray:
    return np.asarray(v, dtype=np.float64).reshape(1, -1)


def _buffer(frames, a1, a2, s1, s2) -> dict:
    """a start or end buffer of 3 frames (MATLAB's 2 x 3 per channel, t 1 x 3)"""
    A = np.vstack([a1, a2])
    S = np.vstack([s1, s2])
    two = np.ones((2, 3))
    t = _row((np.asarray(frames, dtype=np.float64) - 1) * FR)
    return {
        "t": t, "x": 20 * two, "y": 20 * two, "z": 5 * two, "A": A, "c": 100 * two,
        "A_pstd": 0.04 * np.abs(A) + 0.5, "c_pstd": 0.2 * two, "sigma_r": S, "SE_sigma_r": S / 20,
        "pval_Ar": 0.5 * two, "hval_AD": np.zeros((2, 3), dtype=bool),
    }  # fmt: skip


def _amplitude(design: str, c: int, L: int, u4: float, z: np.ndarray) -> np.ndarray:
    """the first channel's amplitudes of a track of L frames in threshold cohort c"""
    k = np.arange(L, dtype=np.float64)
    fall = np.minimum(1.0, (L - k) / 2.0)
    P = 10.0 + 6.0 * c + 4.0 * u4
    if design == "early":  # brighter with the lifetime from the first frame
        return P * np.minimum(1.0, (k + 1) / 3.0) * fall + 1.5 * z[:L]
    if design == "late":  # the same first 5 frames in every cohort
        return 8.0 + (P - 8.0) * np.clip((k - 4.0) / 3.0, 0.0, 1.0) * fall + 1.5 * z[:L]
    return (20.0 + 4.0 * u4) * np.minimum(1.0, (k + 1) / 3.0) * fall + 1.5 * z[:L]  # flat


def synth_tracks(experiment: int, design: str, n: int, frames: int, six: bool) -> list[dict]:
    """n tracks of one experiment (module docstring), in MATLAB's field order; ``six`` False
    leaves out the two longest threshold cohorts"""
    rng = synth.CounterRNG(2026, f"lifetime_c05/{design}/ex{experiment}")
    U = rng.uniform((n, 12))
    Z = rng.normal((n, 4, frames))
    ncoh = len(COHORT_FRAMES) if six else 4
    tracks = []
    for i in range(n):
        u, z = U[i], Z[i]
        cat = CATEGORIES[int(u[0] * len(CATEGORIES))]
        c = int(u[1] * ncoh)
        lo, hi = COHORT_FRAMES[c]
        hi = min(hi, frames - 10)
        L = lo + int(u[2] * (hi - lo + 1))
        if c == 0 and u[2] < 0.08:
            L = 2  # kept by loadTracks, dropped by the lifetime analysis' cutoff of 3 frames
        s = 4 + int(u[3] * (frames - 5 - L))  # the buffers inside the movie
        e = s + L - 1
        k = np.arange(L, dtype=np.float64)
        a1 = _amplitude(design, c, L, u[4], z[0])
        if cat == 1 and L <= 12 and u[8] < 0.1:
            a1[: min(3, L)] += 45.0  # a visitor: a bright start, shorter than 30 s
        sig2 = bool(u[5] < 0.55)
        a2 = 0.5 * a1 + z[1, :L] if sig2 else 1.5 + z[1, :L]
        A = np.vstack([a1, a2])
        S = np.vstack([4.0 + 0.3 * z[2, :L], 3.0 + 0.3 * z[2, :L][::-1]])
        gap = np.zeros((1, L), dtype=bool)
        if L >= 6 and u[9] < 0.25:
            gap[0, L // 2] = True
        f = np.arange(s, e + 1, dtype=np.float64)
        two = np.ones((2, L))
        zb = z[3]
        sb = _buffer([s - 3, s - 2, s - 1], 2 + 0.8 * zb[0:3], 1 + 0.8 * zb[3:6],
                     4 + 0.3 * zb[12:15], 3 + 0.3 * zb[15:18])  # fmt: skip
        eb = _buffer([e + 1, e + 2, e + 3], 2 + 0.8 * zb[6:9], 1 + 0.8 * zb[9:12],
                     4 + 0.3 * zb[18:21], 3 + 0.3 * zb[21:24])  # fmt: skip
        motion = {"MSD": _row(np.arange(1, 11) * 0.01), "MSDstd": _row(np.full(10, 0.001)),
                  "totalDisplacement": 0.5}  # fmt: skip
        tracks.append({
            "t": _row((f - 1) * FR), "f": _row(f),
            "x": np.vstack([10 + 30 * u[10] + 0.1 * k, 10.4 + 30 * u[10] + 0.1 * k]),
            "y": np.vstack([10 + 30 * u[11] - 0.1 * k, 9.7 + 30 * u[11] - 0.1 * k]),
            "z": np.vstack([5 + 0.05 * k, 5 + 0.05 * k]), "A": A, "c": 100 * two,
            "x_pstd": 0.1 * two, "y_pstd": 0.1 * two, "z_pstd": 0.1 * two,
            "A_pstd": 0.04 * np.abs(A) + 0.5, "c_pstd": 0.2 * two, "sigma_r": S,
            "SE_sigma_r": S / 20,
            "pval_Ar": np.vstack([1e-6 * two[0], 0.5 * two[1]]), "isPSF": two.copy(),
            "tracksFeatIndxCG": _row(k + 1 + i), "gapVect": gap, "gapStatus": np.zeros((0, 0)),
            "gapIdx": np.zeros((0, 0)), "seqOfEvents": np.array([[s, 1, 1, NAN], [e, 2, 1, NAN]]),
            "nSeg": 1.0, "visibility": 1.0, "lifetime_s": L * FR, "start": float(s),
            "end": float(e),
            "startBuffer": sb, "endBuffer": eb, "MotionAnalysis": motion, "RSS": 100 * two,
            "hval_Ar": np.vstack([two[0], 0 * two[1]]), "hval_AD": 0 * two, "catIdx": float(cat),
            "hasMultiplePeaks": False, "significantMaster": np.array([[True], [sig2]]),
            "significantVsBackground": np.vstack([two[0] > 0, np.full(L, sig2)]),
            "significantSlave": np.array([[bool(u[6] < 0.3)], [bool(u[7] < 0.5)]]),
        })  # fmt: skip
    return tracks


def write_fixed_mat(path, variables) -> None:
    """a MAT v6 file whose header has no time or platform (the same bytes everywhere)"""
    write_mat(path, variables, version="6")
    data = bytearray(Path(path).read_bytes())
    data[:116] = MAT_HEADER.ljust(116, b" ")
    Path(path).write_bytes(bytes(data))


def frame_name(channel: int, frame: int) -> str:
    cam, wl = ("B", 560) if channel == 1 else ("A", 642)
    t = frame * FRAME_MS
    return (
        f"synth_Cam{cam}_ch{channel - 1}_stack{frame:04d}_{wl}nm_{t:07d}msec_"
        f"{1000000 + t:010d}msecAbs.tif"
    )


def experiment_name(e: int) -> str:
    return f"Ex{e}_synthLT_z0p5"


def make_condition(folder, name: str) -> None:
    """a condition folder: per experiment two channels of tiny raw and DS frames, the settings
    file, and ch1/Analysis/ProcessedTracks.mat"""
    n_ex, design, n, frames, six = CONDITIONS[name]
    root = Path(folder)
    vol = np.full(TINY, 100, dtype=np.uint16)
    for e in range(1, n_ex + 1):
        ex = root / experiment_name(e)
        for c in (1, 2):
            for sub in (ex / f"ch{c}", ex / f"ch{c}" / "DS"):
                sub.mkdir(parents=True, exist_ok=True)
                for k in range(frames):
                    synth.write_tiff(sub / frame_name(c, k), vol)
        with open(ex / "synthLT_Settings.txt", "w", newline="\n") as fh:
            fh.write('Camera type = "Orca4.0"\n')
        (ex / "ch1" / "Analysis").mkdir()
        rows = synth_tracks(e, design, n, frames, six)
        tracks = StructArray((len(rows), 1), TRACK_FIELDS, rows)
        info = {"procFlag": np.ones((1, 2), dtype=bool)}
        write_fixed_mat(ex / "ch1" / "Analysis" / "ProcessedTracks.mat",
                        {"tracks": tracks, "processingInfo": info})  # fmt: skip


def condition_data(root) -> list[dict]:
    """the movies of a condition folder as the port's lifetime analysis takes them (MATLAB's
    data struct: source, framerate, movieLength, channels)"""
    root = Path(root)
    out = []
    for ex in sorted(p for p in root.iterdir() if p.is_dir() and p.name.startswith("Ex")):
        frames = len(list((ex / "ch1").glob("*.tif")))
        out.append({"source": f"{ex}/ch1/", "framerate": FR, "movieLength": frames,
                    "channels": [f"{ex}/ch1/", f"{ex}/ch2/"]})  # fmt: skip
    return out


# lftRes's fields in MATLAB's order, and how each is compared
LFT_FIELDS = tuple(
    "t cellArea lftHist_Ia lftHist_Ib lftHist_IIa nSamples_Ia MaxIntensityThreshold pctCCP pctCS "
    "pctVisit nCCP lftHistCCP lftHistCS lftHistCCP_counts lftHistCS_counts histWeight lftHistVisit "
    "slaveCombs lftHistSlaveAll lftHistSlaveCCP lftHistSlaveCS lftHistSlaveAll_counts "
    "lftHistSlaveCCP_counts lftHistSlaveCS_counts pctSlaveCCP pctSlaveCS startsPerFrameAll "
    "initDensityAll initDensityIa initDensityCCP meanLftHistCCP meanLftHistCS "
    "meanLftHistVisit".split()
)
COHORT_FIELDS = ("bounds", "t", "Amin", "Aplus", "A")
LIFETIME_OUT = (
    [f"lftRes.{k}" for k in LFT_FIELDS]
    + [f"cohorts{i}.{k}" for i in (1, 2) for k in COHORT_FIELDS]
    + [
        "lifetimeData",
        "Parameters.MaxIntensityThreshold",
        "Parameters.nFramesThreshold",
        "Parameters.a",
    ]
)
LIFETIME_SPEC = {
    "lftRes.MaxIntensityThreshold": TOL_THRESHOLD,
    "Parameters.MaxIntensityThreshold": TOL_THRESHOLD,
    **{f"cohorts{i}.{k}": TOL_COHORT for i in (1, 2) for k in ("Amin", "Aplus", "A")},
}
# the fields that follow from the CCP/CS classification by the threshold
CLASSIFIED = (
    "pctCCP pctCS nCCP lftHistCCP lftHistCS lftHistCCP_counts lftHistCS_counts lftHistSlaveCCP "
    "lftHistSlaveCS lftHistSlaveCCP_counts lftHistSlaveCS_counts pctSlaveCCP pctSlaveCS "
    "meanLftHistCCP meanLftHistCS".split()
)
TOL_OTHER_MINIMUM = tol(reason="parity/detection.md L4: another local minimum", rel=0.1)
OTHER_MINIMUM = ignore("parity/detection.md L4: classified by another threshold")
LIFETIME_SPEC_OTHER = {  # two movies: the threshold's fit ends in another local minimum
    **LIFETIME_SPEC,
    "lftRes.MaxIntensityThreshold": TOL_OTHER_MINIMUM,
    "Parameters.MaxIntensityThreshold": TOL_OTHER_MINIMUM,
    **{f"lftRes.{k}": OTHER_MINIMUM for k in CLASSIFIED},
}
LIFETIME_SPEC_ONE = {  # one movie: the median and quartiles of the tracks, no scaling
    **LIFETIME_SPEC,
    **{f"cohorts{i}.{k}": TOL_COHORT_ONE for i in (1, 2) for k in ("Amin", "Aplus", "A")},
}
LIFETIME_CODE = (
    "rt = in{1}; data = XR_loadConditionData3D(rt, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "lftRes = runLifetimeAnalysis3D(data, 'AnalysisPath', 'Analysis', 'overwrite', true); "
    "cohorts = plotIntensityCohorts(data, 'AnalysisPath', 'Analysis'); close all force; "
    "ld = cell(1, numel(data)); for i = 1 : numel(data), "
    "ld{i} = load([data(i).source 'Analysis' filesep 'lifetimeData.mat']); end; "
    "pa = load([data(1).source 'Analysis' filesep 'Parameters.mat']); "
    "out = {"
    + ", ".join(f"lftRes.{k}" for k in LFT_FIELDS)
    + ", "
    + ", ".join(f"cohorts({i}).{k}" for i in (1, 2) for k in COHORT_FIELDS)
    + ", ld, pa.MaxIntensityThreshold, pa.nFramesThreshold, pa.a};"
)


def run_analysis(data, analysis: str, **opts) -> list:
    """the port's lifetime analysis and intensity cohorts of a condition's movies, in the order of
    LIFETIME_OUT"""
    res = lt.run_lifetime_analysis(data, analysis_path=analysis, overwrite=True, **opts)
    coh = lt.plot_intensity_cohorts(data, analysis_path=analysis)
    ld = [read_mat(os.path.join(d["source"], analysis, "lifetimeData.mat")) for d in data]
    pa = read_mat(os.path.join(data[0]["source"], analysis, "Parameters.mat"))
    return (
        [res[k] for k in LFT_FIELDS]
        + [coh[i][k] for i in (0, 1) for k in COHORT_FIELDS]
        + [ld, pa["MaxIntensityThreshold"], pa["nFramesThreshold"], pa["a"]]
    )


def run_condition(root, **opts) -> list:
    """run_analysis of a synthetic condition folder"""
    return run_analysis(condition_data(root), "Analysis", **opts)


# 'Display' 'off'
DISPLAY_OFF_CODE = LIFETIME_CODE.replace(
    "'overwrite', true);", "'overwrite', true, 'Display', 'off');"
)
assert DISPLAY_OFF_CODE != LIFETIME_CODE
# the file names of the tracks and of the lifetime data (one movie)
RENAMED = ("tracks_c05b.mat", "lifetime_c05b.mat")
FILE_NAMES_OUT = [f"lftRes.{k}" for k in LFT_FIELDS] + [
    "lifetimeData",
    "Parameters.MaxIntensityThreshold",
    "Parameters.nFramesThreshold",
    "Parameters.a",
    "default_file_exists",
]
FILE_NAMES_CODE = (
    "rt = in{1}; data = XR_loadConditionData3D(rt, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "lftRes = runLifetimeAnalysis3D(data, 'AnalysisPath', 'Analysis', 'overwrite', true, "
    "'ProcessedTracks', in{2}, 'LifetimeData', in{3}); close all force; "
    "a = [data(1).source 'Analysis' filesep]; ld = load([a in{3}]); "
    "pa = load([a 'Parameters.mat']); out = {"
    + ", ".join(f"lftRes.{k}" for k in LFT_FIELDS)
    + ", ld, pa.MaxIntensityThreshold, pa.nFramesThreshold, pa.a, "
    "exist([a 'lifetimeData.mat'], 'file')};"
)


def make_renamed_condition(folder) -> None:
    """one_late's condition with its tracks under RENAMED[0]"""
    make_condition(folder, "one_late")
    for p in Path(folder).glob("Ex*/ch1/Analysis/ProcessedTracks.mat"):
        p.rename(p.with_name(RENAMED[0]))


def run_file_names(root, tracks: str, lifetime: str) -> list:
    data = condition_data(root)
    res = lt.run_lifetime_analysis(
        data,
        analysis_path="Analysis",
        overwrite=True,
        processed_tracks=tracks,
        lifetime_data=lifetime,
    )
    a = os.path.join(data[0]["source"], "Analysis")
    pa = read_mat(os.path.join(a, "Parameters.mat"))
    return [res[k] for k in LFT_FIELDS] + [
        read_mat(os.path.join(a, lifetime)),
        pa["MaxIntensityThreshold"],
        pa["nFramesThreshold"],
        pa["a"],
        2.0 if os.path.isfile(os.path.join(a, "lifetimeData.mat")) else 0.0,
    ]


# plotIntensityCohorts' LftDataName: the cohorts from the lifetime data file of that name (the
# renamed one the analysis wrote; lifetimeData.mat is not made)
COHORTS_FILE_OUT = [f"cohorts{i}.{k}" for i in (1, 2) for k in COHORT_FIELDS] + [
    "default_file_exists"
]
COHORTS_FILE_CODE = (
    "rt = in{1}; data = XR_loadConditionData3D(rt, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "runLifetimeAnalysis3D(data, 'AnalysisPath', 'Analysis', 'overwrite', true, "
    "'ProcessedTracks', in{2}, 'LifetimeData', in{3}); close all force; "
    "cohorts = plotIntensityCohorts(data, 'AnalysisPath', 'Analysis', 'LftDataName', in{3}); "
    "close all force; a = [data(1).source 'Analysis' filesep]; out = {"
    + ", ".join(f"cohorts({i}).{k}" for i in (1, 2) for k in COHORT_FIELDS)
    + ", exist([a 'lifetimeData.mat'], 'file')};"
)


def run_cohorts_file_name(root, tracks: str, lifetime: str) -> list:
    data = condition_data(root)
    lt.run_lifetime_analysis(
        data,
        analysis_path="Analysis",
        overwrite=True,
        processed_tracks=tracks,
        lifetime_data=lifetime,
    )
    coh = lt.plot_intensity_cohorts(data, analysis_path="Analysis", lft_data_name=lifetime)
    a = os.path.join(data[0]["source"], "Analysis")
    return [coh[i][k] for i in (0, 1) for k in COHORT_FIELDS] + [
        2.0 if os.path.isfile(os.path.join(a, "lifetimeData.mat")) else 0.0
    ]


# the multiple-experiments demo: its frames (for the condition loader: the frame rate and the
# movies' lengths) and the demo run's ProcessedTracks.mat, read and written in a copy
DEMO = "$PETAKIT5D_DEMO_DATA/PetaKit5D_point_detection_demo_dataset"
DEMO_TRACKS = "$PETAKIT5D_DEMO_DATA/PetaKit5D_point_detection_demo_tracks"
DEMO_ANALYSIS = "Analysis_update_3"
DEMO_CODE = (
    "r = in{1}; if r(end) ~= '/', r = [r '/']; end; f = dir(r); "
    "f = f([f.isdir] & ~startsWith({f.name}, '.')); "
    "data = XR_loadConditionData3D([r f(1).name '/'], {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "t = in{2}; if t(end) ~= '/', t = [t '/']; end; "
    "for i = 1 : numel(data), p = strsplit(data(i).source, '/'); e = [t p{end - 2} '/']; "
    "mkdir([e 'ch2']); data(i).source = [e 'ch1/']; data(i).channels{1} = [e 'ch1/']; "
    "data(i).channels{2} = [e 'ch2/']; end; "
    f"lftRes = runLifetimeAnalysis3D(data, 'AnalysisPath', '{DEMO_ANALYSIS}', 'Overwrite', true); "
    f"cohorts = plotIntensityCohorts(data, 'AnalysisPath', '{DEMO_ANALYSIS}'); close all force; "
    "ld = cell(1, numel(data)); for i = 1 : numel(data), "
    f"ld{{i}} = load([data(i).source '{DEMO_ANALYSIS}' filesep 'lifetimeData.mat']); end; "
    f"pa = load([data(1).source '{DEMO_ANALYSIS}' filesep 'Parameters.mat']); "
    "out = {"
    + ", ".join(f"lftRes.{k}" for k in LFT_FIELDS)
    + ", "
    + ", ".join(f"cohorts({i}).{k}" for i in (1, 2) for k in COHORT_FIELDS)
    + ", ld, pa.MaxIntensityThreshold, pa.nFramesThreshold, pa.a};"
)


def demo_data(frames, tracks) -> list[dict]:
    """the demo condition's movies as the port's condition loader gives them (MATLAB's struct),
    their sources and channels in the copy of the tracks"""
    from petakit5d.detection.condition import load_condition_data

    root = Path(frames)
    cond = sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith("."))[0]
    out = []
    for m in load_condition_data(str(cond) + os.sep, ["ch1", "ch2"], ["gfp", "rfp"]):
        ex = Path(tracks) / Path(m.source).parent.name
        (ex / "ch2").mkdir(exist_ok=True)
        d = m.to_matlab()
        d["source"] = str(ex / "ch1") + os.sep
        d["channels"] = [d["source"], str(ex / "ch2") + os.sep]
        out.append(d)
    return out


def run_demo(frames, tracks) -> list:
    return run_analysis(demo_data(frames, tracks), DEMO_ANALYSIS)


def _lifetime_cases() -> list:
    cs = []
    for name in ("three_early", "two_early", "one_late", "flat"):
        cs.append(
            Case(
                f"lifetime.{name}",
                code=LIFETIME_CODE,
                args=[InPath("cond")],
                nargout=len(LIFETIME_OUT),
                outputs=LIFETIME_OUT,
                python=run_condition,
                inputs={"cond": Tree(lambda p, n=name: make_condition(p, n), "cond")},
                compare={"three_early": LIFETIME_SPEC, "two_early": LIFETIME_SPEC_OTHER}.get(
                    name, LIFETIME_SPEC_ONE
                ),
                files=False,
                marks=("slow",),
            )
        )
    # a threshold cohort without tracks: ecdf's error, as in MATLAB (docs/detection_scope.md 6)
    cs.append(
        Case(
            "lifetime.no_long_tracks",
            code=LIFETIME_CODE,
            args=[InPath("cond")],
            nargout=len(LIFETIME_OUT),
            outputs=LIFETIME_OUT,
            python=run_condition,
            inputs={"cond": Tree(lambda p: make_condition(p, "no_long_tracks"), "cond")},
            expect_error="stats:ecdf:NotEnoughData",
            py_error=lt.LifetimeError,
            files=False,
            marks=("slow",),
        )
    )
    cs.append(
        Case(
            "lifetime.display_off",
            code=DISPLAY_OFF_CODE,
            args=[InPath("cond")],
            nargout=len(LIFETIME_OUT),
            outputs=LIFETIME_OUT,
            python=lambda root: run_condition(root, display="off"),
            inputs={"cond": Tree(lambda p: make_condition(p, "three_early"), "cond")},
            compare=LIFETIME_SPEC,
            files=False,
            marks=("slow",),
        )
    )
    cs.append(
        Case(
            "lifetime.file_names",
            code=FILE_NAMES_CODE,
            args=[InPath("cond"), *RENAMED],
            nargout=len(FILE_NAMES_OUT),
            outputs=FILE_NAMES_OUT,
            python=run_file_names,
            inputs={"cond": Tree(make_renamed_condition, "cond")},
            compare=LIFETIME_SPEC,
            files=False,
            marks=("slow",),
        )
    )
    cs.append(
        Case(
            "lifetime.cohorts_file_name",
            code=COHORTS_FILE_CODE,
            args=[InPath("cond"), *RENAMED],
            nargout=len(COHORTS_FILE_OUT),
            outputs=COHORTS_FILE_OUT,
            python=run_cohorts_file_name,
            inputs={"cond": Tree(make_renamed_condition, "cond")},
            compare={
                f"cohorts{i}.{k}": TOL_COHORT_ONE for i in (1, 2) for k in ("Amin", "Aplus", "A")
            },
            files=False,
            marks=("slow",),
        )
    )
    cs.append(
        Case(
            "lifetime.demo",
            code=DEMO_CODE,
            args=[InPath("frames"), InPath("tracks")],
            nargout=len(LIFETIME_OUT),
            outputs=LIFETIME_OUT,
            python=run_demo,
            inputs={"frames": Data(DEMO, copy=False), "tracks": Data(DEMO_TRACKS)},
            compare=LIFETIME_SPEC,
            files=False,
            marks=("demodata", "slow"),
        )
    )
    return cs


# -- the lifetime analysis' functions, as black boxes --------------------------------------------


def sample(n: int, stream: str, mu: float = 10.0, sigma: float = 2.0) -> np.ndarray:
    """n normal values (bit-stable)"""
    return mu + sigma * synth.CounterRNG(2026, f"lifetime_c05/{stream}").normal(n)


def track_maxima(name: str, channel: int) -> list[np.ndarray]:
    """the maxima of the category Ia tracks' amplitudes of each experiment of a synthetic condition
    (columns), as getLifetimeData gives them to scaleEDFs"""
    n_ex, design, n, frames, six = CONDITIONS[name]
    out = []
    for e in range(1, n_ex + 1):
        rows = synth_tracks(e, design, n, frames, six)
        rows = sorted(rows, key=lambda t: -t["lifetime_s"])
        out.append(
            np.array([t["A"][channel].max() for t in rows if t["catIdx"] == 1]).reshape(-1, 1)
        )
    return out


def lognormal(n: int, stream: str, scale: float = 1.0, shift: float = 0.0) -> np.ndarray:
    """n values scale * exp(N(3, 1)) + shift, as columns (the maxima scaleEDFs aligns)"""
    v = scale * el.exp(3.0 + synth.CounterRNG(2026, f"lifetime_c05/{stream}").normal(n)) + shift
    return v.reshape(-1, 1)


def visitor_input(variant: int):
    """getVisitorIndex's lftData: 12 tracks around 30 s, 10 short ones, amplitudes 10 to 15"""
    L = np.concatenate([np.linspace(28.0, 32.0, 12), FR * np.arange(3, 13)])
    if variant == 2:
        L[0] = 30.0
    if variant == 3:
        L[:12] = 30.0
    u = synth.CounterRNG(2026, f"lifetime_c05/visitors{variant}").uniform((L.size, 60, 2))
    A = 10.0 + 5.0 * u
    if variant == 4:
        A[:, 3:, :] = NAN
    if variant == 5:  # a short track exactly at the thresholds: not a visitor (">")
        from petakit5d.matlabcompat import stats as ms

        A[12, :3, 0] = ms.prctile(A[L >= 30, :3, 0], 95, 0).ravel()
    return L.reshape(-1, 1), A


VISITOR_CODE = (
    "L = in{1}; A = in{2}; n = numel(L); lf = struct(); lf.lifetime_s = L; "
    "lf.trackLengths = round(L / 2); lf.start = 5 * ones(n, 1); lf.catIdx = ones(n, 1); "
    "lf.index = (1 : n)'; lf.A = A; lf.A_pstd = A / 10; lf.sigma_r = A / 5; "
    "lf.SE_sigma_r = A / 50; "
    "lf.sbA = A(:, 1 : 3, :); lf.ebA = A(:, 1 : 3, :); lf.sbSigma_r = A(:, 1 : 3, :); "
    "lf.ebSigma_r = A(:, 1 : 3, :); lf.gapMat_Ia = false(n, 60); "
    "lf.significantMaster = true(n, 2); "
    "lf.significantSlave = false(n, 2); out = {getVisitorIndex(lf)};"
)


def _function_cases() -> list:
    cs = []
    for n in (5, 6, 10, 30, 100):
        x = sample(n, f"ad{n}", 4.0, 1.5)
        cs.append(
            Case(
                f"adtest1.estimated.n{n}",
                matlab="adtest1",
                args=[x],
                nargout=4,
                python=lambda x: lt.adtest1(x),
                compare={1: TOL_AD, 2: TOL_AD},
            )
        )
        for mu, sigma in ((4.0, 1.5), (5.0, 0.5)):
            cs.append(
                Case(
                    f"adtest1.known.n{n}.mu{mu:g}.sigma{sigma:g}",
                    matlab="adtest1",
                    args=[x, "mu", mu, "sigma", sigma],
                    nargout=4,
                    python=lambda x, _m, mu, _s, sigma: lt.adtest1(x, mu, sigma),
                    compare={1: TOL_AD, 2: TOL_AD},
                )
            )
    cs.append(
        Case(
            "adtest1.too_few",
            matlab="adtest1",
            args=[sample(4, "ad4")],
            nargout=4,
            python=lambda x: lt.adtest1(x),
            expect_error=True,
            py_error=lt.LifetimeError,
            error_text="withheld",  # a GPL message (decided by the user on 2026-11-20)
        )
    )
    # (one sample, and positions outside [1, n] with fewer than 4 samples, are left out: MATLAB's
    # result then depends on memory beyond its input; docs/parity/detection.md L6)
    for n in (2, 4, 7, 13, 30):
        f = sample(n, f"bi{n}", 0.0, 10.0)
        mirrored = np.array([0.0, 0.25, 0.5, 1.0, n - 0.5, float(n), n + 0.5, n + 1.0])
        for name, xi in (
            ("grid", np.linspace(1.0, n, 3 * n + 2)),
            ("mirrored", mirrored if n >= 4 else None),
            ("integers", np.arange(1.0, n + 1.0)),
        ):
            if xi is None:
                continue
            cs.append(
                Case(
                    f"binterp.n{n}.{name}",
                    matlab="binterp",
                    args=[f, xi],
                    python=lt.binterp,
                    compare=TOL_BINTERP,
                )
            )
        # positions beyond [0, n + 1]: MATLAB's error (no identifier; its message is withheld)
        xo = np.array([0.25, 0.5, 1.0, n - 0.5, float(n), n + 0.5, n + 3.25])
        cs.append(
            Case(
                f"binterp.n{n}.outside",
                matlab="binterp",
                args=[f, xo],
                python=lt.binterp,
                expect_error=True,
                py_error=lt.LifetimeError,
                error_text="withheld",
            )
        )
    cs.append(
        Case(
            "binterp.nan",
            matlab="binterp",
            args=[np.array([1.0, NAN, 3.0, 4.0]), np.array([1.0, 2.5])],
            python=lt.binterp,
            expect_error=True,
            py_error=lt.LifetimeError,
            error_text="withheld",
        )
    )
    for name, make in (
        ("normal_n200", lambda: sample(200, "m200", 10.0, 2.0).reshape(-1, 1)),
        ("normal_n1000", lambda: sample(1000, "m1000", 30.0, 6.0).reshape(-1, 1)),
        ("lognormal_n400", lambda: lognormal(400, "m400")),
    ):
        cs.append(
            Case(
                f"fitGaussianModeToCDF.{name}",
                matlab="fitGaussianModeToCDF",
                args=[In("x")],
                nargout=2,
                python=lt.fit_gaussian_mode_to_cdf,
                inputs={"x": Value(make)},
                compare=TOL_MODE,
            )
        )
    cs.append(
        Case(
            "fitGaussianModeToCDF.constant",
            matlab="fitGaussianModeToCDF",
            args=[np.full(4, 5.0)],
            nargout=2,
            python=lt.fit_gaussian_mode_to_cdf,
            expect_error="optimlib:commonMsgs:UsrObjUndefAtX0",
            py_error=lt.LifetimeError,
        )
    )
    for name, make in (
        ("one", lambda: [lognormal(300, "s1")]),
        ("two", lambda: [lognormal(300, "s1"), lognormal(250, "s2", 1.3)]),
        (
            "three",
            lambda: [
                lognormal(300, "s1"),
                lognormal(250, "s2", 1.3),
                lognormal(400, "s3", 0.8, 2.0),
            ],
        ),
        ("same", lambda: [lognormal(300, "s1"), lognormal(300, "s1")]),
        ("three_early_ch1", lambda: track_maxima("three_early", 0)),
        ("three_early_ch2", lambda: track_maxima("three_early", 1)),
        # the pooled 99th percentile beyond the reference's largest value (its EDF is 1 there)
        (
            "beyond",
            lambda: [lognormal(300, "s1"), lognormal(250, "s2", 1.1), lognormal(200, "s4", 10.0)],
        ),
    ):
        cs.append(
            Case(
                f"scaleEDFs.{name}",
                matlab="scaleEDFs",
                args=[In("s"), "Display", False],
                nargout=3,
                python=lambda s, _d, _v: lt.scale_edfs(s),
                inputs={"s": Value(make)},
                compare=TOL_SCALE if name in ("one", "two", "three", "same") else TOL_SCALE_PATH,
            )
        )
    for v in (1, 2, 3, 4, 5):
        cs.append(
            Case(
                f"getVisitorIndex.variant{v}",
                code=VISITOR_CODE,
                args=[In("L"), In("A")],
                python=lambda L, A: [lt.visitor_index({"lifetime_s": L, "A": A})],
                inputs={
                    "L": Value(lambda v=v: visitor_input(v)[0]),
                    "A": Value(lambda v=v: visitor_input(v)[1]),
                },
            )
        )
    return cs


# -- matlabcompat: hist, normcdf, erf, erfc, lsqnonlin -------------------------------------------

HIST_CENTERS = np.array([1.0, 2.0, 4.0, 7.0, 11.0])
HIST_Y = np.array(
    [0, 1, 1.5, 2, 3, 3.0000000000000004, 4, 5.5, 5.4999999999999991, 7, 9, 11, 12, NAN, INF, -INF]
)


def erf_grid() -> np.ndarray:
    """1/1024 steps in [-7, 7] (exact binary fractions), decades down to 1e-300, the edges"""
    dense = [k / 1024 for k in range(-7168, 7169)]
    small = [float(f"{m}e-{k}") for k in range(1, 301, 3) for m in (1, 3)]
    big = [8.0, 10.0, 20.0, 26.5, 27.25, 28.0, 30.0, 1e3, 1e300, INF]
    edges = [
        0.84375,
        1.25,
        2.857142857142857,
        6.0,
        28.0,
        5e-324,
        2.2250738585072014e-308,
        -0.0,
        NAN,
    ]
    v = dense + small + [-s for s in small] + big + [-b for b in big] + edges
    return np.array(v, dtype=np.float64)


def small_grid(g: np.ndarray, n: int = 200) -> np.ndarray:
    fin = g[np.isfinite(g)]
    idx = np.unique(np.linspace(0, fin.size - 1, n - 8).astype(np.int64))
    return np.concatenate([fin[idx], g[~np.isfinite(g)]])[:n]


def normcdf_grid() -> np.ndarray:
    return np.array([k / 100 for k in range(-4000, 6001)] + [-INF, INF, NAN], dtype=np.float64)


# lsqnonlin problems: (MATLAB function of (p, a, b), Python function, a, b, x0, lb, ub, tol)
def _decay_data():
    t = np.arange(20, dtype=np.float64) / 4
    y = 3.0 * el.exp(-0.7 * t) + 0.05 * np.cos(np.arange(20, dtype=np.float64))
    return t, y


def _cdf_data():
    x = np.sort(sample(150, "lsqcdf", 5.0, 1.2))
    f = (np.arange(1, 151, dtype=np.float64) - 0.5) / 150
    return x, f


def _decay(p, t, y):
    return p[0] * el.exp(-p[1] * t) - y


def _cdf(p, x, f):
    return (p[0] * (0.5 * (1 + special.erf((x - p[1]) / math.sqrt(2) / p[2]))) - f) / x.size


def _rosen(p, a, b):
    return np.array([10.0 * (p[1] - p[0] * p[0]), 1.0 - p[0]])


LSQ_PROBLEMS = {
    "decay": (
        "@(p, t, y) p(1) * exp(-p(2) * t) - y",
        _decay,
        _decay_data,
        [1.0, 0.1],
        [0.0, 0.0],
        [INF, 10.0],
        1e-8,
    ),
    "decay_bound": (
        "@(p, t, y) p(1) * exp(-p(2) * t) - y",
        _decay,
        _decay_data,
        [1.0, 0.1],
        [0.0, 0.0],
        [2.5, 10.0],
        1e-8,
    ),
    "cdf": (
        "@(p, x, f) (p(1) * (0.5 * (1 + erf((x - p(2)) / sqrt(2) / p(3)))) - f) / numel(x)",
        _cdf,
        _cdf_data,
        [0.5, 4.0, 2.0],
        [0.0, 0.0, 0.0],
        [1.0, INF, INF],
        1e-8,
    ),
    "cdf_tol6": (
        "@(p, x, f) (p(1) * (0.5 * (1 + erf((x - p(2)) / sqrt(2) / p(3)))) - f) / numel(x)",
        _cdf,
        _cdf_data,
        [0.5, 4.0, 2.0],
        [0.0, 0.0, 0.0],
        [1.0, INF, INF],
        1e-6,
    ),
    "rosenbrock_bound": (
        "@(p, a, b) [10 * (p(2) - p(1) ^ 2); 1 - p(1)]",
        _rosen,
        lambda: (0.0, 0.0),
        [-1.2, 1.0],
        [-2.0, -2.0],
        [0.5, 2.0],
        1e-10,
    ),
}


def _lsq_cases() -> list:
    cs = []
    for name, (mfun, pfun, data, x0, lb, ub, tl) in LSQ_PROBLEMS.items():
        code = (
            f"f = {mfun}; a = in{{1}}; b = in{{2}}; "
            "opts = optimset('Display', 'off', 'TolFun', in{6}, "
            "'TolX', in{6}, 'MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Jacobian', 'off'); "
            "[x, rn] = lsqnonlin(@(p) f(p, a, b), in{3}, in{4}, in{5}, opts); out = {x, rn};"
        )

        def py(a, b, x0, lb, ub, tl, pfun=pfun):
            r = optim.lsqnonlin(
                pfun, x0, lb, ub, (a, b), tol_fun=tl, tol_x=tl, max_fun_evals=10000, max_iter=10000
            )
            return r.x, r.resnorm

        cs.append(
            Case(
                f"lsqnonlin.{name}",
                code=code,
                args=[In("a"), In("b"), np.array(x0), np.array(lb), np.array(ub), tl],
                nargout=2,
                python=py,
                inputs={"a": Value(lambda d=data: d()[0]), "b": Value(lambda d=data: d()[1])},
                compare=TOL_LSQ,
            )
        )
    return cs


def _compat_cases() -> list:
    cs = [
        Case("hist.row", matlab="hist", args=[HIST_Y, HIST_CENTERS], python=mhist.hist),
        Case(
            "hist.column",
            matlab="hist",
            args=[HIST_Y.reshape(-1, 1), HIST_CENTERS],
            python=mhist.hist,
        ),
        Case(
            "hist.column_centers",
            matlab="hist",
            args=[HIST_Y, HIST_CENTERS.reshape(-1, 1)],
            python=mhist.hist,
        ),
        Case(
            "hist.int16",
            matlab="hist",
            args=[np.arange(13, dtype=np.int16), HIST_CENTERS],
            python=mhist.hist,
        ),
        Case(
            "hist.matrix",
            matlab="hist",
            args=[np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]]), HIST_CENTERS],
            python=lambda y, c: mhist.hist(y, c, matlab_layout=True),
        ),
        Case(
            "hist.lifetimes",
            matlab="hist",
            args=[FR * np.array([3, 3, 4, 9, 10, 25, 49, 50, 51, 60, 1]), FR * np.arange(3, 51)],
            python=mhist.hist,
        ),
        Case(
            "hist.empty",
            matlab="hist",
            args=[np.zeros((0, 1)), np.arange(1.0, 6.0)],
            python=mhist.hist,
        ),
    ]
    for fn, pf in (("erf", special.erf), ("erfc", special.erfc)):
        cs.append(
            Case(f"{fn}.grid", matlab=fn, args=[In("x")], python=pf, inputs={"x": Value(erf_grid)})
        )
        cs.append(Case(f"{fn}.small", matlab=fn, args=[small_grid(erf_grid())], python=pf))
    cs += [
        Case(
            "normcdf.grid",
            matlab="normcdf",
            args=[In("x"), 8.92968222709834, 0.7999357568764726],
            python=special.normcdf,
            inputs={"x": Value(normcdf_grid)},
        ),
        Case(
            "normcdf.small",
            matlab="normcdf",
            args=[small_grid(normcdf_grid()), 0.0, 1.0],
            python=special.normcdf,
        ),
    ]
    return cs + _lsq_cases()


CASES = _lifetime_cases() + _function_cases() + _compat_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- rules (no MATLAB) ----------------------------------------------------------------------------


def test_hist_counts_a_value_on_an_edge_in_the_lower_bin():
    assert mhist.hist(HIST_Y, HIST_CENTERS).tolist() == [4, 2, 4, 2, 3]
    assert mhist.hist([], [1.0, 2.0]).tolist() == [0, 0]
    with pytest.raises(ValueError):
        mhist.hist([1.0], [2.0, 1.0])


def test_normcdf_is_erfc_of_the_standardized_value():
    x = np.array([-3.0, 0.0, 0.5, 2.0, NAN])
    assert np.array_equal(
        special.normcdf(x, 1.0, 2.0),
        0.5 * special.erfc(-((x - 1.0) / 2.0) / math.sqrt(2)),
        equal_nan=True,
    )


def test_erf_edges():
    assert special.erf(np.array([INF, -INF, 0.0])).tolist() == [1.0, -1.0, 0.0]
    assert special.erfc(np.array([INF, -INF])).tolist() == [0.0, 2.0]
    assert np.isnan(special.erf(np.array([NAN])))[0]


def test_lsqnonlin_stays_in_the_bounds_and_raises_at_an_undefined_start():
    r = optim.lsqnonlin(
        _rosen, [-1.2, 1.0], [-2.0, -2.0], [0.5, 2.0], (0, 0), tol_fun=1e-10, tol_x=1e-10
    )
    assert r.x[0] <= 0.5 and abs(r.x[0] - 0.5) < 1e-6 and abs(r.x[1] - 0.25) < 1e-6
    with pytest.raises(optim.LsqnonlinError) as e:
        optim.lsqnonlin(lambda p: np.array([NAN, p[0]]), [1.0])
    assert e.value.identifier == "optimlib:commonMsgs:UsrObjUndefAtX0"


def test_lsqnonlin_steps_back_from_a_bound_and_differences_inward():
    # a trust-region step that would leave the box stops 0.95 of the way to the bound
    x, lb, ub = np.array([0.5, 0.5]), np.zeros(2), np.ones(2)
    g, H = np.array([-1.0, -1e-3]), np.array([[-1.0, 0.0], [0.0, 1.0]])
    s = optim._step(x, g, H, lb, ub, 10.0)[0]
    assert abs((x[0] + s[0]) - (0.5 + 0.95 * 0.5)) < 1e-15 and np.all(x + s < ub)
    # at the upper bound the forward difference goes down
    seen = []

    def lin(p):
        seen.append(float(p[0]))
        return np.array([p[0] - 10.0])

    optim.lsqnonlin(lin, [1.0], [-INF], [1.0], tol_fun=1e-12, tol_x=1e-12, max_iter=1)
    assert seen[1] == 1.0 - math.sqrt(np.finfo(float).eps)


def test_lsqnonlin_does_not_depend_on_the_blas():
    # the arithmetic is written out: the same bits from any thread count or array layout
    t, y = _decay_data()
    a = optim.lsqnonlin(
        _decay, [1.0, 0.1], [0.0, 0.0], [INF, 10.0], (t, y), tol_fun=1e-8, tol_x=1e-8
    )
    b = optim.lsqnonlin(
        _decay,
        [1.0, 0.1],
        [0.0, 0.0],
        [INF, 10.0],
        (np.ascontiguousarray(t[::-1][::-1]), y.copy()),
        tol_fun=1e-8,
        tol_x=1e-8,
    )
    assert np.array_equal(a.x, b.x) and a.func_count == b.func_count


def test_binterp_interpolates_its_samples_and_mirrors():
    f = sample(9, "bi_rule", 0.0, 5.0)
    assert np.allclose(lt.binterp(f, np.arange(1.0, 10.0)), f, rtol=0, atol=1e-12)
    xi = np.array([1.25, 1.75])
    assert np.allclose(lt.binterp(f, 2 - xi), lt.binterp(f, xi), rtol=0, atol=1e-12)  # mirror at 1
    with pytest.raises(lt.LifetimeError):
        lt.binterp(f, [0.5, 10.5])  # beyond n + 1


def test_adtest1_needs_five_values_and_flags_a_shifted_sample():
    with pytest.raises(lt.LifetimeError):
        lt.adtest1([1.0, 2.0, 3.0, 4.0])
    h, p, stat, cv = lt.adtest1(sample(50, "ad_rule", 0.0, 1.0) + 3.0, 0.0, 1.0)
    assert h and math.isnan(p) and stat > cv == 2.492


def test_the_lifetime_data_of_synthetic_tracks(tmp_path):
    tracks = synth_tracks(1, "early", 40, 60, True)
    path = tmp_path / "ProcessedTracks.mat"
    write_fixed_mat(path, {"tracks": StructArray((len(tracks), 1), TRACK_FIELDS, tracks)})
    tr = lt.load_tracks(str(path))
    d = lt.lifetime_data_from_tracks(tr, 60, 2)
    # sorted by lifetime, longest first (ties in file order), as MATLAB's loadTracks
    tracks = sorted(tracks, key=lambda t: -t["lifetime_s"])
    assert d["lifetime_s"].ravel().tolist() == [t["lifetime_s"] for t in tracks]
    ia = [t for t in tracks if t["catIdx"] == 1]
    assert list(d) == list(lt.LIFETIME_FIELDS)
    assert d["A"].shape == (len(ia), 60, 2) and d["sbA"].shape == (len(ia), 3, 2)
    assert d["significantMaster"].dtype == bool and d["significantMaster"].shape == (len(tracks), 2)
    first = ia[0]
    n = first["A"].shape[1]
    assert np.array_equal(d["A"][0, :n, :], first["A"].T) and np.isnan(d["A"][0, n:, :]).all()
    assert np.isnan(d["A"]).any() and d["A"][np.isnan(d["A"])].view(np.uint64)[0] == np.float64(
        NAN
    ).view(np.uint64)
    assert d["index"].ravel().tolist() == list(range(1, len(tracks) + 1))


def test_movies_of_different_frame_rates_are_refused(tmp_path):
    # MATLAB: MATLAB:InputParser:ArgumentFailedValidation (measured on three_early's movies)
    data = [
        {"source": f"{tmp_path}/m{i}/", "framerate": fr, "movieLength": 60, "channels": [1, 2]}
        for i, fr in enumerate((2.0, 3.0))
    ]
    with pytest.raises(lt.LifetimeError) as e:
        lt.run_lifetime_analysis(data)
    assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"


def test_the_interface_is_matlabs_and_unported_options_are_refused(tmp_path):
    from petakit5d import matlab as pm
    from petakit5d.params import MatlabArgumentError, find

    spec = find("runLifetimeAnalysis3D")
    assert spec.wrapper is lt.run_lifetime_analysis or spec.func is lt.run_lifetime_analysis
    assert [p.matlab for p in spec.params][:6] == [
        "data", "lb", "ub", "Display", "DisplayMode", "AnalysisPath"
    ]  # fmt: skip
    movie = {"source": f"{tmp_path}/m/ch1/", "framerate": 2.0, "movieLength": 60.0}
    # MATLAB's parser refuses these (tests/params/test_third_party_checks.py)
    for name, value in (("Overwrite", "yes"), ("AnalysisPath", 3.0), ("Display", "maybe")):
        with pytest.raises(MatlabArgumentError) as e:
            pm.runLifetimeAnalysis3D(movie, name, value)
        assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
    # options MATLAB takes that the port does not: before anything is read
    for kw in ({"cutoff_f": 4}, {"lb": (1, 11, 16, 21, 41, 60)}, {"print_": True},
               {"max_intensity_threshold": 120.0}, {"rescale": False}, {"type": "Ia"}):  # fmt: skip
        with pytest.raises(NotImplementedError, match="only the default is ported"):
            lt.run_lifetime_analysis([movie], **kw)
    with pytest.raises(NotImplementedError, match="Binning"):
        pm.runLifetimeAnalysis3D(movie, "Binning", 2.0)


def test_the_movies_come_as_structs_struct_arrays_or_conditions():
    class FakeMovie:  # a movie of detection.condition: its struct by to_matlab()
        def __init__(self, d):
            self.d = d

        def to_matlab(self):
            return dict(self.d)

    a, b = {"source": "/m1/ch1/", "framerate": 2.0}, {"source": "/m2/ch1/", "framerate": 2.0}
    assert lt._movies(a) == [a]
    assert lt._movies(StructArray((1, 2), ("source", "framerate"), [a, b])) == [a, b]
    assert lt._movies([FakeMovie(a), FakeMovie(b)]) == [a, b]
    assert lt._movies((a, b)) == [a, b]


def test_a_cell_mask_is_refused(tmp_path):
    # MATLAB would compute cell areas and initiation densities from it, which are not ported
    src = tmp_path / "m" / "ch1"
    (src / "Detection").mkdir(parents=True)
    (src / "Detection" / "cellmask.tif").write_bytes(b"")
    movie = {"source": f"{src}{os.sep}", "framerate": 2.0, "movieLength": 60.0}
    with pytest.raises(NotImplementedError, match="cell masks"):
        lt.run_lifetime_analysis(movie)


def test_the_source_takes_the_platforms_separator(monkeypatch):
    # MATLAB's [data.source AnalysisPath filesep ...]: the separator is the platform's
    from petakit5d.core import paths

    monkeypatch.setattr(paths, "filesep", "/")  # Linux and macOS
    assert lt._source({"source": "/data/m/ch1"}) == "/data/m/ch1/"
    assert lt._source({"source": "/data/m/ch1/"}) == "/data/m/ch1/"
    monkeypatch.setattr(paths, "filesep", "\\")  # Windows
    assert lt._source({"source": "C:\\data\\m\\ch1"}) == "C:\\data\\m\\ch1\\"
    assert lt._source({"source": "C:\\data\\m\\ch1\\"}) == "C:\\data\\m\\ch1\\"
    assert lt._source({"source": "C:/data/m/ch1/"}) == "C:/data/m/ch1/"


def test_a_category_without_tracks_is_ecdfs_error():
    with pytest.raises(lt.LifetimeError) as e:
        lt.fit_gaussian_mode_to_cdf(np.zeros(0))
    assert e.value.identifier == "stats:ecdf:NotEnoughData"
    with pytest.raises(lt.LifetimeError) as e:
        lt.fit_gaussian_mode_to_cdf(np.full(4, 5.0))
    assert e.value.identifier == "optimlib:commonMsgs:UsrObjUndefAtX0"


def test_the_reference_edf_is_one_past_its_largest_value(monkeypatch):
    # measured bit for bit on the arguments of MATLAB's fits (the demo's second channel)
    seen = []
    real = optim.lsqnonlin

    def spy(fun, x0, lb=None, ub=None, args=(), **kw):
        seen.append(args)
        return real(fun, x0, lb, ub, args, **kw)

    monkeypatch.setattr(optim, "lsqnonlin", spy)
    samples = [lognormal(300, "s1"), lognormal(250, "s2", 1.1), lognormal(200, "s4", 10.0)]
    a, c, ref = lt.scale_edfs(samples)
    x_ref = np.sort(samples[ref - 1].ravel())
    _, _, refE, grid = seen[0]
    beyond = grid > x_ref[-1]
    assert beyond.any() and np.all(refE[beyond] == 1.0) and not np.isnan(refE).any()


def test_visitors_are_the_short_tracks_above_the_long_tracks_percentile():
    L, A = visitor_input(1)
    vis = lt.visitor_index({"lifetime_s": L, "A": A})
    long_ = L.ravel() >= 30
    assert not vis[long_].any()
    A2 = A.copy()
    A2[12, :3, 0] = 1e3  # a short track with a bright start
    assert lt.visitor_index({"lifetime_s": L, "A": A2})[12]
