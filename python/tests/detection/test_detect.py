"""The 3-D point-source detection (task C02): the background estimate, the local maxima, the LoG
kernel, the Anderson-Darling test, the MathWorks functions detection reaches (histcounts,
findpeaks, medfilt1, rangesearch) and the lowSNR detector, against MATLAB R2024b through the
oracle harness, and the detection's rules on synthetic spots (MATLAB's PointDetectionTest in
intent: spots found within 0.25 voxel with amplitudes within 10 %, none in noise).

The detector cases reuse C01's MATLAB references (``test_references.py``: the same MATLAB code
and inputs, so the same keys) and add the Python side. The runner_threshold cases (task C02d)
run both runners on a designed condition whose frames tell the runners' residual threshold
apart from the camera threshold and from the noise, and ``XR_runDetection3D`` without a
settings file.
"""

import importlib.util
from pathlib import Path

import numpy as np

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.detection import detect, fit
from petakit5d.matlabcompat import neighbors, signal
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value
from petakit5d.testing.compare import EXACT, tol
from petakit5d.testing.compare import tree as tree_spec

_spec = importlib.util.spec_from_file_location(
    "detection_references", Path(__file__).with_name("test_references.py")
)
REF = importlib.util.module_from_spec(_spec)
# its cases register under this name, apart from the module pytest collects
REF.__file__ = str(_spec.origin) + ":helpers"
_spec.loader.exec_module(REF)

TOL_LOG = tol(reason="parity D1: fspecial3('log') arithmetic order", abs=1e-16, rel=2e-15)
TOL_DETECT = tol(
    reason="parity D3: the detector's fits (3-D filter order, covariance)", rel=1e-9, abs=1e-9
)
TOL_SPEEDUP_MEX = tol(
    reason="parity D7: the 'speedup' detector with the fitGaussian3D mex's fits (D5)",
    rel=1e-6,
    abs=1e-9,
)


def _val(a: np.ndarray) -> Value:
    return Value(lambda: a)


def _rng(stream: str) -> synth.CounterRNG:
    return synth.CounterRNG(2026, "C02." + stream)


def _nan_band(vol_m: np.ndarray) -> np.ndarray:
    """a DS-like NaN region (outside a slanted band in x over z) of a MATLAB (y, x, z) volume"""
    v = vol_m.copy()
    ny, nx, nz = v.shape
    x = np.arange(nx)[None, :, None]
    z = np.arange(nz)[None, None, :]
    v[np.broadcast_to((x < 2 * z + 3) | (x > 2 * z + nx - 25), v.shape)] = np.nan
    return v


def _bg_inputs():
    a = REF.ds_frame(0, 1)
    b = _nan_band(REF.ds_frame(3, 2))
    c = a[:30, :40, :8].copy()
    c[:5] = 0
    return {"a": a, "b": b, "c": c}


def _py_background(vol_m, _name, seed):
    return detect.estimate_background(from_matlab(vol_m), seed)[None, :]


def _py_locmax(vol_m, w):
    return to_matlab(detect.locmax3d(from_matlab(vol_m), w))


def _py_log(size, sigma):
    return to_matlab(detect.fspecial3_log(size, sigma))


def _ad_samples():
    g = _rng("ad").normal((6, 294))
    out = {}
    for i, s in enumerate((0.7, 0.9, 1.0, 1.1, 1.3, 1.0)):
        x = g[i] * s
        out[f"x{i}"] = np.sort(x) if i < 5 else x  # the last one out of order
    return out


def _py_adtest(x, mu, sigma, alpha, case):
    return float(fit.adtest(x, mu, sigma, alpha, int(case)))


def _signal_inputs():
    r = _rng("signal")
    sr = np.abs(r.normal((1, 800))[0]) * 4 + 3
    y = np.round(np.abs(r.normal((1, 120))[0]) * 20)
    y[40:44] = 55  # a flat top
    return {"sr": sr, "y": y}


def _py_histcounts(x, n):
    c, e = signal.histcounts(x, int(n))
    return [c[None, :], e[None, :]]


def _py_findpeaks(y, _name, m):
    pk, lc = signal.findpeaks(y, m)
    return [pk[None, :], lc[None, :]]


def _py_medfilt1(y):
    return signal.medfilt1(y)[None, :]


def _points():
    p = _rng("points").normal((60, 3)) * 2
    p[10] = p[11] + 0.25  # a close pair
    return {"P": p}


def _py_rangesearch(X, Y, r):
    idx, d = neighbors.rangesearch(X, Y, r)
    return [[(i + 1).astype(np.float64)[None, :] for i in idx], [v[None, :] for v in d]]


def _ps_matlab(ps: dict) -> dict:
    """detections (Python form) as MATLAB's struct: 1 x n rows, 1-based positions"""
    out = {}
    for k in (
        "x",
        "y",
        "z",
        "A",
        "s",
        "c",
        "x_pstd",
        "y_pstd",
        "z_pstd",
        "A_pstd",
        "s_pstd",
        "c_pstd",
        "x_init",
        "y_init",
        "z_init",
        "sigma_r",
        "SE_sigma_r",
        "RSS",
        "pval_Ar",
        "hval_Ar",
        "hval_AD",
        "isPSF",
    ):
        v = np.asarray(ps[k])
        if k in ("x", "y", "z", "x_init", "y_init", "z_init"):
            v = v + 1
        out[k] = v if v.ndim == 2 else v[None, :]
    return out


def _py_lowsnr(vol_m, sigma):
    detect.clear_threshold_cache()
    vol = from_matlab(vol_m)
    bi = detect.estimate_background(vol, 1)
    ps, mask = detect.point_source_detection(vol, sigma, background_info=np.r_[bi, 3.0])
    return [bi[None, :], _ps_matlab(ps), to_matlab(mask)]


def _py_speedup(vol_m, sigma):
    """XR_pointSourceDetection3D(vol, sigma, 'BackgroundCorrection', false): its defaults
    (mode 'xyzAc', 'CD_Newton', redundancy radius 0.25)"""
    detect.clear_speedup_cache()
    ps, mask = detect.point_source_detection_speedup(from_matlab(vol_m), sigma)
    return [_ps_matlab(ps), to_matlab(mask)]


def _py_speedup_sigma_form(psf_m, sigma):
    """the sigma estimate's call: mode 'xyzAsrc' with the fitGaussian3D mex"""
    detect.clear_speedup_cache()
    ps, mask = detect.point_source_detection_speedup(
        from_matlab(psf_m), sigma, mode="xyzAsrc", fit_method="original"
    )
    return [_ps_matlab(ps), to_matlab(mask)]


def _cases():
    bg = _bg_inputs()
    cs = []
    for name, seed in (("a", 1), ("b", 2), ("c", 7)):
        cs.append(
            Case(
                f"background.{name}",
                matlab="XR_estimate_overall_background_information",
                args=[In(name), "randseed", float(seed)],
                python=_py_background,
                inputs={name: _val(bg[name])},
            )
        )
    cs.append(
        Case(
            "locmax3d.nan",
            matlab="XR_locmax3d",
            args=[In("b"), [3.0, 3.0, 3.0]],
            python=_py_locmax,
            inputs={"b": _val(bg["b"])},
        )
    )
    for size, sigma in (([7.0, 7.0, 7.0], [1.3, 1.3, 1.6]), ([3.0, 3.0, 3.0], [1.5, 1.5, 1.9])):
        cs.append(
            Case(
                f"fspecial3_log.{int(size[0])}",
                matlab="fspecial3",
                args=["log", size, sigma],
                python=lambda _k, size, sigma: _py_log(size, sigma),
                compare=TOL_LOG,
            )
        )
    ad = _ad_samples()
    for k in ad:
        for s in (0.8, 1.0, 1.25):
            cs.append(
                Case(
                    f"adtest.{k}.s{s}",
                    matlab="adtest_mex",
                    args=[In(k), 0.0, s, 0.05, 2.0],
                    python=_py_adtest,
                    inputs={k: _val(ad[k])},
                )
            )
    si = _signal_inputs()
    for n in (60.0, 341.0, 17.0):
        cs.append(
            Case(
                f"histcounts.n{int(n)}",
                matlab="histcounts",
                args=[In("sr"), n],
                nargout=2,
                python=_py_histcounts,
                inputs={"sr": _val(si["sr"])},
            )
        )
    # one finite value (the bad-point filter's histogram of few residual sigmas: the limits
    # v -/+ 0.5 to multiples of 0.5, a tie toward v), NaN and infinities beside it, none finite
    inf = np.inf
    for name, x, n in (
        ("one", [2.814], 60.0),
        ("equal", [6.36, 6.36, 6.36], 17.0),
        ("tie", [0.75, 0.75], 13.0),
        ("with_nan", [0.25, NAN, 0.25, inf], 4.0),
        ("none_finite", [NAN, inf, -inf], 13.0),
        ("finite_and_not", [3.0, NAN, 5.5, -inf, 4.0], 13.0),
    ):
        cs.append(
            Case(
                f"histcounts.{name}.n{int(n)}",
                matlab="histcounts",
                args=[In("x"), n],
                nargout=2,
                python=_py_histcounts,
                inputs={"x": _val(np.array([x], dtype=np.float64))},
            )
        )
    cs.append(
        Case(
            "findpeaks.prominence5",
            matlab="findpeaks",
            args=[In("y"), "minpeakprominence", 5.0],
            nargout=2,
            python=_py_findpeaks,
            inputs={"y": _val(si["y"])},
        )
    )
    cs.append(
        Case(
            "medfilt1.order3",
            matlab="medfilt1",
            args=[In("y")],
            python=_py_medfilt1,
            inputs={"y": _val(si["y"])},
        )
    )
    pts = _points()
    cs.append(
        Case(
            "rangesearch.r1p5",
            matlab="rangesearch",
            args=[In("P"), In("P"), 1.5],
            nargout=2,
            python=_py_rangesearch,
            inputs={"P": _val(pts["P"])},
        )
    )
    for c in (1, 2):
        for k in (0, 8):
            cs.append(
                Case(
                    f"psd_lowsnr.ch{c}.frame{k}",
                    code=REF.CLEAR_CODE
                    + "vol = in{1}; bi = XR_estimate_overall_background_information(vol, "
                    + "'randseed', 1); [ps, mask] = XR_pointSourceDetection3D_3(vol, in{2}, "
                    + REF.LOWSNR_OPTS
                    + ", 'BackgroundInfo', [bi, 3]); out = {bi, ps, mask};",
                    args=[In("vol"), REF.SIGMA_DS[c]],
                    nargout=3,
                    inputs={"vol": Value(lambda c=c, k=k: REF.ds_frame(k, c))},
                    python=_py_lowsnr,
                    compare={0: EXACT, 1: TOL_DETECT, 2: EXACT},
                )
            )
    # a frame with two spots: the bad-point filter's scalar true for fewer than three points,
    # with which MATLAB's detector keeps both (task C02c)
    cs.append(
        Case(
            "psd_lowsnr.two_spots",
            code=REF.CLEAR_CODE
            + "vol = in{1}; bi = XR_estimate_overall_background_information(vol, "
            + "'randseed', 1); [ps, mask] = XR_pointSourceDetection3D_3(vol, in{2}, "
            + REF.LOWSNR_OPTS
            + ", 'BackgroundInfo', [bi, 3]); out = {bi, ps, mask};",
            args=[In("vol"), [1.3, 1.6]],
            nargout=3,
            inputs={"vol": Value(lambda: to_matlab(_two_spot_frame()))},
            python=_py_lowsnr,
            compare={0: EXACT, 1: TOL_DETECT, 2: EXACT},
        )
    )
    # C01's references of the 'speedup' detector (PointDetectionTest's form, the sigma
    # estimate's form): the same MATLAB code and inputs, so the same keys
    for c in (1, 2):
        cs.append(
            Case(
                f"psd_speedup_form.ch{c}",
                code=REF.CLEAR_CODE
                + "[ps, mask] = XR_pointSourceDetection3D(in{1}, in{2}, "
                + "'BackgroundCorrection', false); out = {ps, mask};",
                args=[In("vol"), REF.SIGMA_DS[c]],
                nargout=2,
                inputs={"vol": Value(lambda c=c: REF.ds_frame(4, c))},
                python=_py_speedup,
                compare={0: TOL_DETECT, 1: EXACT},
            )
        )
        cs.append(
            Case(
                f"psd_sigma_estimate_form.ch{c}",
                code=REF.CLEAR_CODE
                + "[ps, mask] = XR_pointSourceDetection3D(in{1}, in{2}, "
                + REF.SIGMA_EST_OPTS
                + "); out = {ps, mask};",
                args=[In("psf"), [1.2, 2.5]],
                nargout=2,
                inputs={"psf": Value(lambda c=c: to_matlab(REF.psf_volume(c).astype(np.float64)))},
                python=_py_speedup_sigma_form,
                compare={0: TOL_SPEEDUP_MEX, 1: EXACT},
            )
        )
    # XR_runDetTrack3d's detection stage (XR_runDetection3D_frame_parallel with the wrapper's
    # options: 'lowSNR' without the background correction) on C01's synthetic condition: every
    # file it writes (Detection3D.mat, the per-frame Detections/Detection_F_<k>.mat, the masks)
    cs.append(
        Case(
            "frame_parallel.wrapper_options",
            code=REF.CLEAR_CODE
            + REF.LOAD
            + "apath = {[data(1).source 'Analysis' filesep]}; "
            + "XR_runDetection3D_frame_parallel(data, 'Sigma', in{2}, 'Overwrite', false, "
            + "'WindowSize', [], 'ResultsPath', apath, 'Mode', 'xyzAc', 'FitMixtures', false, "
            + "'MaxMixtures', 2, 'DetectionMethod', 'lowSNR', 'BackgroundCorrection', false); "
            + "out = {};",
            args=[InPath("cond"), SIGMA_RUN],
            nargout=0,
            python=_py_frame_parallel,
            inputs={"cond": Tree(REF.make_condition, "cond")},
            tree=tree_spec(per_file={"*.mat": TOL_DETECT}),
            marks=("slow",),
        )
    )
    return cs


def _py_frame_parallel(rt, sigma):
    import os

    from petakit5d.detection import condition as cond
    from petakit5d.detection import run

    detect.clear_threshold_cache()  # the case starts a MATLAB session's state
    data = cond.load_condition_data(
        os.path.join(rt, REF.EXPERIMENT) + os.sep, ["ch1", "ch2"], ["gfp", "rfp"]
    )
    rp = [data[0].source + "Analysis" + os.sep]
    run.run_detection_frame_parallel(
        data,
        sigma=np.asarray(sigma),
        results_path=rp,
        overwrite=False,
        window_size=None,
        mode="xyzAc",
        fit_mixtures=False,
        max_mixtures=2,
        detection_method="lowSNR",
        background_correction=False,
    )


SIGMA_RUN = np.array([REF.SIGMA_DS[1], REF.SIGMA_DS[2]])

# -- the runners' residual threshold (task C02d) ---------------------------------------------------
# a condition of three frames whose noise lies in each range of the runners' rule (below 2/3 of
# the camera threshold 3, between 1 and 1.5 of it, above 1.5 of it), each with three bright
# spots (photon-like noise: residual sigmas above every level) and a weak one in a patch of
# other noise, whose residual sigma lies between the rule's level and another rule's: kept with
# 0.9 x 3 = 2.7 and not with 3 (frame 1), left out with the noise 3.55 and kept with 3 (frame 2),
# kept with 3 and not with the noise 5.55 (frame 3)
RT_EXPERIMENT = "Ex1_rt_z0p5"
RT_SIZE = (16, 56, 64)  # (z, y, x)
RT_SIGMA = np.array([[1.3, 1.6], [1.3, 1.6]])  # [sigma_xy sigma_z] of both channels
RT_BRIGHT = [((5.2, 12.4, 14.6), 220.0), ((9.6, 40.3, 50.2), 260.0), ((11.1, 20.7, 47.4), 190.0)]
# per frame: the noise, the weak spot's patch noise, its (z, y, x) and amplitude
RT_FRAMES = (
    (1.6, 2.7, (8.3, 36.6, 22.4), 9.0),
    (3.5, 3.1, (8.4, 36.2, 22.7), 9.0),
    (5.5, 4.0, (7.6, 35.8, 23.1), 12.0),
)


def _rt_counts(v: np.ndarray) -> np.ndarray:
    t = np.trunc(v)
    v = t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)
    return np.clip(v, 0, 65535).astype(np.uint16)


def rt_frame(k: int) -> np.ndarray:
    """frame k (0-based) of the threshold condition, (z, y, x) uint16"""
    noise, pn, wpos, wamp = RT_FRAMES[k]
    sig = synth.beads(
        RT_SIZE, [p for p, _ in RT_BRIGHT], (1.6, 1.3), [a for _, a in RT_BRIGHT], 0.0
    )
    weak = synth.beads(RT_SIZE, [wpos], (1.6, 1.3), [wamp], 0.0)
    sig = sig.astype(np.float64)
    std = np.full(RT_SIZE, noise)
    z, y, x = (int(round(t)) for t in wpos)
    std[max(z - 5, 0) : z + 6, max(y - 5, 0) : y + 6, max(x - 5, 0) : x + 6] = pn
    std = np.sqrt(std**2 + sig)
    r = synth.CounterRNG(2026, "C02d.rt2." + ("low", "mid", "high")[k]).normal(RT_SIZE)
    return _rt_counts(100.0 + sig + weak.astype(np.float64) + std * r)


def _make_rt(folder, settings: bool) -> None:
    root = Path(folder)
    ex = root / RT_EXPERIMENT
    for c in (1, 2):
        for sub in (ex / f"ch{c}", ex / f"ch{c}" / "DS"):
            sub.mkdir(parents=True, exist_ok=True)
        for k in range(len(RT_FRAMES)):
            vol = rt_frame(k)
            synth.write_tiff(ex / f"ch{c}" / REF.frame_name(c, k), vol)
            synth.write_tiff(ex / f"ch{c}" / "DS" / REF.frame_name(c, k), vol)
    if settings:
        with open(ex / "rt_Settings.txt", "w", newline="\n") as f:
            f.write('Camera type = "Orca4.0"\n')


def make_rt_condition(folder) -> None:
    """the threshold condition: two channels of the same frames, raw and DS, the settings file's
    camera line (camera threshold 3)"""
    _make_rt(folder, True)


def make_rt_condition_no_settings(folder) -> None:
    """the threshold condition without a settings file: no camera threshold (the noise)"""
    _make_rt(folder, False)


RT_LOAD = (
    "rt = in{1}; ex = [rt '" + RT_EXPERIMENT + "' filesep]; "
    "data = XR_loadConditionData3D(ex, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "apath = {[data(1).source 'Analysis' filesep]}; "
)


def _rt_data(rt):
    import os

    from petakit5d.detection import condition as cond

    data = cond.load_condition_data(
        os.path.join(rt, RT_EXPERIMENT) + os.sep, ["ch1", "ch2"], ["gfp", "rfp"]
    )
    return data, [data[0].source + "Analysis" + os.sep]


def _py_rt_run(rt, sigma):
    from petakit5d.detection import run

    detect.clear_threshold_cache()
    data, rp = _rt_data(rt)
    run.run_detection(
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
    )


def _py_rt_frame_parallel(rt, sigma):
    from petakit5d.detection import run

    detect.clear_threshold_cache()
    data, rp = _rt_data(rt)
    run.run_detection_frame_parallel(
        data,
        sigma=np.asarray(sigma),
        results_path=rp,
        overwrite=False,
        window_size=None,
        mode="xyzAc",
        fit_mixtures=False,
        max_mixtures=2,
        detection_method="lowSNR",
        background_correction=True,
    )


RT_RUN = (
    "XR_runDetection3D(data, 'Sigma', in{2}, 'Overwrite', true, 'WindowSize', [], "
    "'ResultsPath', apath, 'Mode', 'xyzAc', 'FitMixtures', false, 'MaxMixtures', 1, "
    "'DetectionMethod', 'lowSNR', 'BackgroundCorrection', true, 'WriteAmira', false); out = {};"
)
RT_FRAME_PARALLEL = (
    "XR_runDetection3D_frame_parallel(data, 'Sigma', in{2}, 'Overwrite', false, "
    "'WindowSize', [], 'ResultsPath', apath, 'Mode', 'xyzAc', 'FitMixtures', false, "
    "'MaxMixtures', 2, 'DetectionMethod', 'lowSNR', 'BackgroundCorrection', true); out = {};"
)


def _rt_cases():
    cs = []
    for name, code, py, make in (
        ("runner_threshold.run_detection", RT_RUN, _py_rt_run, make_rt_condition),
        (
            "runner_threshold.frame_parallel",
            RT_FRAME_PARALLEL,
            _py_rt_frame_parallel,
            make_rt_condition,
        ),
        # without a settings file the camera search stops with an error, which XR_runDetection3D
        # catches (the threshold is the noise); the frame-parallel runner stops (D13)
        ("runner_threshold.no_settings", RT_RUN, _py_rt_run, make_rt_condition_no_settings),
    ):
        cs.append(
            Case(
                name,
                code=REF.CLEAR_CODE + RT_LOAD + code,
                args=[InPath("cond"), RT_SIGMA],
                nargout=0,
                python=py,
                inputs={"cond": Tree(make, "cond")},
                tree=tree_spec(per_file={"*.mat": TOL_DETECT}),
            )
        )
    return cs


CASES = _cases() + _rt_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- rules (no MATLAB) ----------------------------------------------------------------------------


def _spots_volume(seed: str, spots, noise=5.0, bg=100.0, shape=(24, 64, 64), sigma=(1.3, 1.6)):
    """a (z, y, x) volume of Gaussian spots ((x, y, z), amplitude) on a background with noise"""
    Z, Y, X = np.meshgrid(*(np.arange(n, dtype=np.float64) for n in shape), indexing="ij")
    v = np.full(shape, bg)
    for (x, y, z), a in spots:
        v += a * np.exp(
            -((X - x) ** 2 + (Y - y) ** 2) / (2 * sigma[0] ** 2)
            - (Z - z) ** 2 / (2 * sigma[1] ** 2)
        )
    return v + noise * _rng(seed).normal(shape)


def _two_spot_frame() -> np.ndarray:
    """the oracle case's two-spot frame in whole counts, as a camera gives them (exp's last bits
    follow the CPU's vector unit, the rounded counts do not)"""
    v = _spots_volume("two", SPOTS[:2])
    t = np.trunc(v)
    return t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)


SPOTS = [
    ((12.3, 14.6, 8.2), 120.0),
    ((40.7, 20.1, 12.5), 90.0),
    ((25.2, 48.8, 15.9), 150.0),
    ((50.4, 45.3, 6.6), 110.0),
]


def test_spots_are_found_where_they_are():
    detect.clear_threshold_cache()
    v = _spots_volume("spots", SPOTS)
    bi = detect.estimate_background(v, 1)
    ps, _ = detect.point_source_detection(v, (1.3, 1.6), background_info=np.r_[bi, 3.0])
    P = np.stack([ps["x"], ps["y"], ps["z"]], 1)
    for (x, y, z), a in SPOTS:
        d = np.abs(P - np.array([x, y, z])).max(axis=1)
        i = int(np.argmin(d))
        assert d[i] < 0.25
        assert abs(ps["A"][i] - a) < 0.1 * a
    assert len(P) == len(SPOTS)


def test_two_spots_are_kept_by_the_bad_point_filter():
    """a frame with two detections: the bad-point filter's scalar true for fewer than three
    points keeps both (MATLAB's detector on frames of one and two spots, task C02c's probe fu1;
    the oracle case psd_lowsnr.two_spots)"""
    detect.clear_threshold_cache()
    v = _spots_volume("two", SPOTS[:2])
    bi = detect.estimate_background(v, 1)
    ps, _ = detect.point_source_detection(v, (1.3, 1.6), background_info=np.r_[bi, 3.0])
    P = np.stack([ps["x"], ps["y"], ps["z"]], 1)
    assert len(P) == 2
    for (x, y, z), _a in SPOTS[:2]:
        assert np.abs(P - np.array([x, y, z])).max(axis=1).min() < 0.25


def test_noise_gives_no_detection():
    detect.clear_threshold_cache()
    v = _spots_volume("noise", [])
    bi = detect.estimate_background(v, 1)
    ps, _ = detect.point_source_detection(v, (1.3, 1.6), background_info=np.r_[bi, 3.0])
    assert ps["x"].size == 0


def test_the_t_threshold_is_cached_per_alpha():
    detect.clear_threshold_cache()
    v = _spots_volume("cache", SPOTS[:2])
    a = detect.lowsnr_prefit(v, (1.3, 1.6)).t_threshold
    b = detect.lowsnr_prefit(v, (1.9, 2.5)).t_threshold
    detect.clear_threshold_cache()
    c = detect.lowsnr_prefit(v, (1.9, 2.5)).t_threshold
    assert a == b and b != c


def test_adtest_takes_its_sample_in_the_given_order():
    x = _rng("order").normal((1, 300))[0]
    assert not fit.adtest(np.sort(x), 0.0, 1.0)
    assert fit.adtest(x, 0.0, 1.0)


def test_fit_recovers_a_spot_in_a_window():
    s, sz = 1.3, 1.6
    Z, Y, X = np.meshgrid(np.arange(7.0), np.arange(7.0), np.arange(7.0), indexing="ij")
    w = 100 + 80 * np.exp(
        -((X - 3.3) ** 2 + (Y - 2.8) ** 2) / (2 * s * s) - (Z - 3.1) ** 2 / (2 * sz * sz)
    )
    R = fit.fit_cd([w], [[3, 3, 3, 50, s, sz, 90]])
    # MATLAB's stop rule (RSS change below 1e-7 of the first RSS, steps below 1e-4) ends early
    assert np.allclose(R.prm[0, :3], [3.3, 2.8, 3.1], atol=1e-5)
    assert np.allclose(R.prm[0, [3, 6]], [80, 100], rtol=1e-4)


def _locmax_ref(a, window):
    """XR_locmax3d's measured rule in plain NumPy: NaN as -1e8, strictly larger than every other
    voxel of the face-clipped window, the replaced value times 1 or 0"""
    from scipy import ndimage as ndi

    w = np.asarray(window, dtype=np.int64).ravel()
    b = np.where(np.isnan(a), -1e8, a)
    hz, hy, hx = int(w[2]) // 2, int(w[0]) // 2, int(w[1]) // 2
    d = ndi.maximum_filter(
        b, size=(2 * hz + 1, 2 * hy + 1, 2 * hx + 1), mode="constant", cval=-np.inf
    )
    cand = np.nonzero(b == d)
    unique = np.ones(cand[0].size, dtype=bool)
    nz, ny, nx = b.shape
    for dz in range(-hz, hz + 1):
        for dy in range(-hy, hy + 1):
            for dx in range(-hx, hx + 1):
                if dz == dy == dx == 0:
                    continue
                zz, yy, xx = cand[0] + dz, cand[1] + dy, cand[2] + dx
                ok = (zz >= 0) & (zz < nz) & (yy >= 0) & (yy < ny) & (xx >= 0) & (xx < nx)
                tie = np.zeros_like(ok)
                tie[ok] = b[zz[ok], yy[ok], xx[ok]] == b[cand][ok]
                unique &= ~tie
    m = np.zeros(b.shape)
    m[tuple(c[unique] for c in cand)] = 1.0
    return b * m


def test_locmax3d_follows_its_rule():
    a = np.round(_rng("locmax").normal((9, 14, 11)) * 3)  # many ties
    a[2, 3, 4] = np.nan
    a[5, :, 0] = np.nan
    a[0, 0, 0] = -0.0
    for w in ((3, 3, 3), (5, 5, 7), (3, 5, 1)):
        got, ref = detect.locmax3d(a, w), _locmax_ref(a, w)
        assert got.tobytes() == ref.tobytes()


def test_prefit_statistics_are_the_numpy_expressions():
    r = _rng("prefit")
    fg, fu, fu2 = (
        np.abs(r.normal((4, 6, 5))) * s + o for s, o in ((30, 200), (40, 3000), (500, 1e5))
    )
    fu2[0, 0, :2] = 0.0  # negative RSS clipped to 0
    gsum, g2sum, n, k = 54.88, 20.01, 343, 1.6448536269514722
    A, c, sres, T, df = detect._prefit_stats(fg, fu, fu2, gsum, g2sum, n, k, 1.0, 2)
    eA = (fg - gsum * fu / n) / (g2sum - gsum * gsum / n)
    ec = (fu - eA * gsum) / n
    with np.errstate(invalid="ignore"):
        rss = eA * eA * g2sum - 2 * eA * (fg - ec * gsum) + fu2 + n * ec * ec - 2 * ec * fu
        rss = np.where(rss < 0, 0.0, rss)
        es = np.sqrt(rss / (n - 1))
        C11 = float(np.linalg.inv(np.array([[g2sum, gsum], [gsum, float(n)]]))[0, 0])
        sA = np.sqrt(rss / (n - 3) * C11)
        SE = es / np.sqrt(2.0 * (n - 1)) * k
        eT = (eA - es * k) / np.sqrt((sA * sA + SE * SE) / n)
        u = sA * sA + SE * SE
        edf = (n - 1) * (u * u) / ((sA * sA) * (sA * sA) + (SE * SE) * (SE * SE))
    for got, ref in ((A, eA), (c, ec), (sres, es), (T, eT), (df, edf)):
        assert got.tobytes() == ref.tobytes()


def _movie(root, frames=3) -> dict:
    """the synthetic condition's movie as MATLAB's struct holds it (the fields run_detection
    reads; framePathsDS a column of paths per channel, as condition.Movie gives it)"""
    import os

    ex = Path(root) / REF.EXPERIMENT
    chans = [str(ex / "ch1") + os.sep, str(ex / "ch2") + os.sep]
    cols = []
    for c in (1, 2):
        paths = [str(p) for p in sorted((ex / f"ch{c}" / "DS").glob("*.tif"))[:frames]]
        col = np.empty((len(paths), 1), dtype=object)
        col[:, 0] = paths
        cols.append(col)
    return {"source": chans[0], "channels": chans, "framePathsDS": cols, "zAniso": 2.5}


def test_run_detection_writes_detections_and_masks(tmp_path):
    """PointDetectionTest's runDetectionWritesDetectionsAndMask in intent: the 'speedup' run
    writes Detection3D.mat with one entry per frame and one mask per frame, no Amira folder
    with WriteAmira false; the camera's threshold comes from the settings file"""
    import os

    from petakit5d.detection import run
    from petakit5d.io import matfile

    REF.make_condition(tmp_path)
    detect.clear_speedup_cache()
    res = tmp_path / "Analysis"
    (entries,) = run.run_detection(
        _movie(tmp_path),
        sigma=SIGMA_RUN,
        results_path=[str(res) + os.sep],
        detection_method="speedup",
        write_amira=False,
    )
    assert len(entries) == 3 and all(np.size(e["x"]) > 0 for e in entries)
    assert sorted(p.name for p in (res / "Masks").iterdir()) == [
        "dmask_1.tif",
        "dmask_2.tif",
        "dmask_3.tif",
    ]
    assert not (res / "Amira").exists()
    fi = matfile.read_mat(str(res / "Detection3D.mat"))["frameInfo"]
    assert len(fi) == 3


def test_run_detection_takes_a_condition(tmp_path):
    """XR_runDetection3D's condition, as the wrapper passes it (measured, probes rd1 and rd2):
    a 1 x n struct array, one results folder per movie; with Overwrite false (the default) a
    movie whose Detection3D.mat exists is skipped and its file kept; the default folder is
    [source 'Analysis' filesep] and the default method 'speedup'; a folder without the trailing
    separator gets 'Detection3D.mat' beside it (string concatenation) and its masks inside"""
    import os

    from petakit5d.detection import run
    from petakit5d.io import matfile
    from petakit5d.io.matfile import StructArray

    REF.make_condition(tmp_path)
    m = _movie(tmp_path, frames=1)
    data = StructArray((1, 2), list(m), [m, dict(m)])
    ra, rb = str(tmp_path / "RA") + os.sep, str(tmp_path / "RB")
    os.makedirs(ra)
    matfile.write_mat(ra + "Detection3D.mat", {"marker": 1.0})
    detect.clear_speedup_cache()
    out = run.run_detection(
        data, sigma=SIGMA_RUN, results_path=[ra, rb], detection_method="speedup", write_amira=False
    )
    assert out[0] is None and list(matfile.read_mat(ra + "Detection3D.mat")) == ["marker"]
    assert len(out[1]) == 1 and os.path.isfile(rb + "Detection3D.mat")
    assert os.listdir(rb) == ["Masks"] and os.listdir(os.path.join(rb, "Masks")) == ["dmask_1.tif"]
    # the defaults: the Analysis folder, 'speedup' (the same detections), skipped when there
    (got,) = run.run_detection(m, sigma=SIGMA_RUN, write_amira=False)
    default = m["source"] + "Analysis" + os.sep + "Detection3D.mat"
    assert os.path.isfile(default)
    assert np.array_equal(got[0]["x"], out[1][0]["x"], equal_nan=True)
    stamp = os.path.getmtime(default)
    assert run.run_detection(m, sigma=SIGMA_RUN, write_amira=False) == [None]
    assert os.path.getmtime(default) == stamp


def test_the_wrapper_runs_the_detection(tmp_path):
    """XR_runDetTrack3d's detection stage ('parfor') with run_detection: the condition as
    load_condition_data gives it and the wrapper's settings; Detection3D.mat in each movie's
    Analysis folder, a mask per frame"""
    import os

    from petakit5d.detection import condition as cond
    from petakit5d.detection import run, wrapper
    from petakit5d.io import matfile

    REF.make_condition(tmp_path)
    data = cond.load_condition_data(str(tmp_path), ["ch1", "ch2"], ["gfp", "rfp"])
    old = dict(wrapper.STAGES)
    try:
        wrapper.STAGES["run_detection"] = run.run_detection
        wrapper.run_detection_tracking(data, sigma=SIGMA_RUN, parallel_method="parfor", track=False)
    finally:
        wrapper.STAGES.clear()
        wrapper.STAGES.update(old)
    res = os.path.join(data[0].source, "Analysis")
    fi = matfile.read_mat(os.path.join(res, "Detection3D.mat"))["frameInfo"]
    assert len(fi) == REF.N_FRAMES and all(np.size(e["x"]) > 0 for e in fi)
    assert len(os.listdir(os.path.join(res, "Masks"))) == REF.N_FRAMES


def test_run_detection_sigma_rows_past_the_channels(tmp_path):
    """the wrapper's default sigma has three rows: for two channels MATLAB's frameInfo keeps two
    rows per field and the whole sigma in s (probe rd2: 16 bytes more per frame, no other
    change)"""
    import os

    from petakit5d.detection import run

    REF.make_condition(tmp_path)
    detect.clear_speedup_cache()
    s3 = np.vstack([SIGMA_RUN, [1.58, 1.608]])
    (entries,) = run.run_detection(
        _movie(tmp_path, frames=1),
        sigma=s3,
        results_path=[str(tmp_path / "R") + os.sep],
        write_amira=False,
    )
    e = entries[0]
    assert np.asarray(e["x"]).shape[0] == 2 and np.asarray(e["isPSF"]).shape[0] == 2
    assert np.array_equal(e["s"], s3) and len(e["dRange"]) == 2


def test_run_detection_refuses_what_is_not_ported(tmp_path):
    import os

    import pytest

    from petakit5d.detection import run

    REF.make_condition(tmp_path)
    m = _movie(tmp_path, frames=1)
    rp = [str(tmp_path / "R") + os.sep]
    for kw, err in (
        ({"window_size": [2, 2]}, NotImplementedError),
        ({"mode": "xyzAsc"}, NotImplementedError),
        ({"fit_mixtures": True}, NotImplementedError),
        ({"detection_method": "original"}, NotImplementedError),
        # MATLAB's parameters the port takes at their defaults only (the registration's)
        ({"remove_redundant": False}, NotImplementedError),
        ({"redundancy_radius": 0.25}, NotImplementedError),
        ({"master": 2}, NotImplementedError),
        ({"alpha": 0.01}, NotImplementedError),
        ({"file_name": "D.mat"}, NotImplementedError),
        ({"filter_by_residual_sigma": False}, NotImplementedError),
    ):
        with pytest.raises(err):
            run.run_detection(m, sigma=SIGMA_RUN, results_path=rp, **kw)
    # the defaults in MATLAB's forms pass the guard ([] for none, a logical, a 1 x 1 double)
    matlab_forms = (np.array([[True]]), np.array([[0.5]]), True, np.zeros((0, 0)), 0.05, [],
                    "Detection3D.mat", 1.0)  # fmt: skip
    run._defaults_only("run_detection", matlab_forms)
    with pytest.raises(TypeError):  # sigma and the folders by name only
        run.run_detection(m, SIGMA_RUN, rp)
    with pytest.raises(TypeError):  # MATLAB's 'ResultsPath' must be a cell
        run.run_detection(m, sigma=SIGMA_RUN, results_path=rp[0])
    with pytest.raises(ValueError):
        run.run_detection(m, sigma=None, results_path=rp)
    assert not os.path.exists(rp[0])


def _fp_entries(path):
    from petakit5d.io import matfile

    fi = matfile.read_mat(path, "frameInfo")["frameInfo"]
    return [fi] if isinstance(fi, dict) else list(fi)


def _same_entry(a, b) -> bool:
    return all(
        np.array_equal(np.asarray(a[k]), np.asarray(b[k]), equal_nan=True)
        for k in a
        if k != "dRange"
    ) and all(np.array_equal(x, y) for x, y in zip(a["dRange"], b["dRange"], strict=True))


def test_frame_parallel_writes_run_detection_s_files(tmp_path):
    """XR_runDetection3D_frame_parallel (measured, probes fp1-fp4): the same Detection3D.mat and
    masks as XR_runDetection3D, plus Detections/Detection_F_<k>.mat per frame (k not padded,
    one variable frameInfo_k, the frame's entry) and an empty Amira/Detection folder"""
    import os

    from petakit5d.detection import run
    from petakit5d.io import matfile, read_tiff

    REF.make_condition(tmp_path)
    m = _movie(tmp_path)
    ra, rb = str(tmp_path / "RA") + os.sep, str(tmp_path / "RB") + os.sep
    detect.clear_threshold_cache()
    (ea,) = run.run_detection(m, sigma=SIGMA_RUN, results_path=[ra], detection_method="lowSNR")
    detect.clear_threshold_cache()
    (eb,) = run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rb], detection_method="lowSNR"
    )
    fa, fb = _fp_entries(ra + "Detection3D.mat"), _fp_entries(rb + "Detection3D.mat")
    assert len(fa) == len(fb) == 3 and all(_same_entry(a, b) for a, b in zip(fa, fb, strict=True))
    assert sorted(os.listdir(rb + "Detections")) == [f"Detection_F_{k}.mat" for k in (1, 2, 3)]
    for k in (1, 2, 3):
        v = matfile.read_mat(rb + f"Detections/Detection_F_{k}.mat")
        assert list(v) == ["frameInfo_k"] and _same_entry(v["frameInfo_k"], fb[k - 1])
        name = run.dmask_name(k, 3)
        assert np.array_equal(read_tiff(ra + "Masks/" + name), read_tiff(rb + "Masks/" + name))
    assert os.listdir(rb + "Amira/Detection") == []
    assert len(eb) == 3 and np.array_equal(np.asarray(eb[0]["x"]), np.asarray(ea[0]["x"]))


def test_frame_parallel_reuses_frames_and_skips_movies(tmp_path):
    """MATLAB's job framework (probes fp1-fp3): a frame whose Detections file exists is not
    detected again (its entry goes into Detection3D.mat as it is), a movie whose Detection3D.mat
    exists is skipped with overwrite false; with overwrite true the port detects every frame
    again (MATLAB keeps both files then, ignoring the flag: parity D11)"""
    import os

    from petakit5d.detection import run
    from petakit5d.io import matfile

    REF.make_condition(tmp_path)
    m = _movie(tmp_path)
    rp = str(tmp_path / "R") + os.sep
    detect.clear_threshold_cache()
    (first,) = run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rp], detection_method="lowSNR"
    )
    # frame 2's file changed, Detection3D.mat removed: frame 2 is taken as it is
    f2 = rp + "Detections/Detection_F_2.mat"
    e2 = matfile.read_mat(f2)["frameInfo_k"]
    x = np.array(e2["x"], dtype=np.float64)
    x[0, 0] = 999.0
    e2["x"] = x
    matfile.write_mat(f2, {"frameInfo_k": e2}, version="7.3")
    os.remove(rp + "Detection3D.mat")
    (again,) = run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rp], detection_method="lowSNR"
    )
    assert np.asarray(_fp_entries(rp + "Detection3D.mat")[1]["x"])[0, 0] == 999.0
    assert np.asarray(again[1]["x"])[0, 0] == 999.0
    # Detection3D.mat exists: skipped
    stamp = os.path.getmtime(rp + "Detection3D.mat")
    assert run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rp], detection_method="lowSNR"
    ) == [None]
    assert os.path.getmtime(rp + "Detection3D.mat") == stamp
    # overwrite: every frame detected again
    detect.clear_threshold_cache()
    (fresh,) = run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rp], detection_method="lowSNR", overwrite=True
    )
    assert np.array_equal(np.asarray(fresh[1]["x"]), np.asarray(first[1]["x"]))
    assert np.asarray(_fp_entries(rp + "Detection3D.mat")[1]["x"])[0, 0] != 999.0


def test_frame_parallel_keeps_a_detection3d_another_caller_wrote(tmp_path, monkeypatch):
    """callers sharing a movie's frames: when Detection3D.mat appears while this caller's tasks
    end (another caller gathered the frames first), it is kept as that caller wrote it"""
    import os

    from petakit5d.detection import run
    from petakit5d.jobs import framework

    REF.make_condition(tmp_path)
    m = _movie(tmp_path)
    rp = str(tmp_path / "R") + os.sep
    orig = framework.run_tasks

    def tasks_then_another_caller(*a, **kw):
        out = orig(*a, **kw)
        with open(kw["final_out_fullpath"], "wb") as f:
            f.write(b"another caller's file")
        return out

    monkeypatch.setattr(framework, "run_tasks", tasks_then_another_caller)
    detect.clear_threshold_cache()
    (entries,) = run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rp], detection_method="lowSNR"
    )
    assert len(entries) == 3
    with open(rp + "Detection3D.mat", "rb") as f:
        assert f.read() == b"another caller's file"


def test_frame_parallel_takes_over_a_stale_marker(tmp_path):
    """a frame's Detection_F_<k>.tmp marker (the job framework's, as MATLAB's runner leaves it
    while a frame is detected) older than the framework's wait time is taken over: the frame is
    detected here and its marker removed (a fresh marker is left to its caller: the framework's
    rule, tests/jobs)"""
    import os
    import time

    from petakit5d.detection import run

    REF.make_condition(tmp_path)
    m = _movie(tmp_path)
    rp = str(tmp_path / "R") + os.sep
    os.makedirs(rp + "Detections")
    marker = rp + "Detections/Detection_F_2.tmp"
    open(marker, "w").close()
    old = time.time() - 600
    os.utime(marker, (old, old))
    detect.clear_threshold_cache()
    (out,) = run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rp], detection_method="lowSNR"
    )
    assert len(out) == 3 and os.path.exists(rp + "Detection3D.mat")
    assert sorted(os.listdir(rp + "Detections")) == [f"Detection_F_{k}.mat" for k in (1, 2, 3)]


def test_frame_parallel_empty_frame_gets_the_empty_entry(tmp_path, monkeypatch):
    """a frame without detections: the empty entry, as XR_runDetection3D gives it (MATLAB's
    frame-parallel runner repeats the session's previous frame's detections there: parity
    D12, probe fp3)"""
    import os

    from petakit5d.detection import run

    REF.make_condition(tmp_path)
    m = _movie(tmp_path)
    rp = str(tmp_path / "R") + os.sep
    orig = detect.point_source_detection
    calls = []

    def no_detections_in_frame_2(vol, sigma, **kw):
        calls.append(1)
        ps, mask = orig(vol, sigma, **kw)
        if len(calls) == 2:
            ps = {k: (np.asarray(v)[..., :0]) for k, v in ps.items()}
        return ps, mask

    monkeypatch.setattr(detect, "point_source_detection", no_detections_in_frame_2)
    detect.clear_threshold_cache()
    (entries,) = run.run_detection_frame_parallel(
        m, sigma=SIGMA_RUN, results_path=[rp], detection_method="lowSNR"
    )
    e = entries[1]
    assert np.asarray(e["x"]).shape == (0, 0) and np.asarray(e["isPSF"]).shape == (0, 0)
    assert np.array_equal(e["s"], SIGMA_RUN) and len(e["dRange"]) == 2
    assert np.asarray(entries[0]["x"]).size > 0 and np.asarray(entries[2]["x"]).size > 0


def test_frame_parallel_refuses_what_is_not_ported(tmp_path):
    import os

    import pytest

    from petakit5d.detection import run

    REF.make_condition(tmp_path)
    m = _movie(tmp_path, frames=1)
    rp = [str(tmp_path / "R") + os.sep]
    with pytest.raises(NotImplementedError, match="XR_pointSourceDetection3D_1"):
        run.run_detection_frame_parallel(m, sigma=SIGMA_RUN, results_path=rp)  # MATLAB's default
    for kw in ({"detection_method": "adaptive"}, {"detection_method": "original"},
               {"detection_method": "lowSNR", "window_size": [2, 2]},
               {"detection_method": "lowSNR", "fit_mixtures": True},
               {"detection_method": "lowSNR", "remove_secondary_redundant": False},
               {"detection_method": "lowSNR", "cell_mask": np.ones((2, 2))}):  # fmt: skip
        with pytest.raises(NotImplementedError):
            run.run_detection_frame_parallel(m, sigma=SIGMA_RUN, results_path=rp, **kw)
    assert not os.path.exists(rp[0])


def test_the_wrapper_runs_the_frame_parallel_detection(tmp_path):
    """XR_runDetTrack3d's default detection stage ('lock_distributed', not on Windows) with
    run_detection_frame_parallel: Detection3D.mat and a file per frame in each Analysis folder"""
    import os

    from petakit5d.detection import condition as cond
    from petakit5d.detection import run, wrapper
    from petakit5d.io import matfile

    REF.make_condition(tmp_path)
    data = cond.load_condition_data(str(tmp_path), ["ch1", "ch2"], ["gfp", "rfp"])
    old = dict(wrapper.STAGES)
    try:
        wrapper.STAGES["run_detection_frame_parallel"] = run.run_detection_frame_parallel
        wrapper.STAGES["run_detection"] = run.run_detection
        wrapper.run_detection_tracking(data, sigma=SIGMA_RUN, track=False)
    finally:
        wrapper.STAGES.clear()
        wrapper.STAGES.update(old)
    res = os.path.join(data[0].source, "Analysis")
    fi = matfile.read_mat(os.path.join(res, "Detection3D.mat"))["frameInfo"]
    assert len(fi) == REF.N_FRAMES
    if os.name != "nt":  # Windows: XR_runDetTrack3d's parfor (run_detection)
        assert len(os.listdir(os.path.join(res, "Detections"))) == REF.N_FRAMES


def test_fit_cd_counts_no_iteration_for_a_window_too_small():
    """a window of fewer than 7 voxels is not fitted: its stats are NaN and its iteration count
    0 on every platform (a NaN cast to int64 differs between x86-64 and ARM)"""
    w = np.full((3, 3, 3), np.nan)
    w[1, 1, :] = (100.0, 180.0, 100.0)
    w[1, 0, 1] = 120.0
    r = fit.fit_cd(
        [w, 100.0 + np.zeros((3, 3, 3))], np.array([[1, 1, 1, 50.0, 1.0, 1.0, 100.0]] * 2)
    )
    assert r.iterations.dtype == np.int64 and r.iterations[0] == 0 and np.isnan(r.rss[0])


def test_fit_is_nan_from_the_window_centre():
    """XR_fitGaussians3D's 'xyzAc' fit is NaN when it ends w or more from the window's centre,
    the rounded start (MATLAB's pattern on 2,315 traced fits from the other channel's
    positions): a start at 11.6 (centre 12, w = 2) with the spot at 9.8 is NaN though 1.8 from
    the start; with the spot at 13.8 it stands though 2.2 from the start, at 14.3 it is NaN"""
    z, y, x = np.mgrid[0:21, 0:21, 0:21].astype(np.float64)
    for spot, nan in ((9.8, True), (13.8, False), (14.3, True)):
        vol = 100.0 + 500.0 * np.exp(-((x - spot) ** 2 + (y - 10) ** 2 + (z - 10) ** 2) / 2.0)
        r = fit.fit_gaussians3d(
            vol,
            np.array([[11.6, 10.0, 10.0]]),
            np.array([500.0]),
            [1.0, 1.0],
            np.array([100.0]),
            "xyzAc",
        )
        assert bool(np.isnan(r["x"][0])) is nan, (spot, r["x"][0])
        if not nan:
            assert abs(r["x"][0] - spot) < 1e-3


def test_redundant_fits_keep_the_lowest_rss_of_each_ball():
    ps = {
        "x": np.array([1.0, 1.2, 1.4, 5.0]),
        "y": np.zeros(4),
        "z": np.zeros(4),
        "RSS": np.array([3.0, 1.0, 2.0, 9.0]),
    }
    # 0's ball holds 1 (a lower RSS), 2's ball holds 1 too: only 1 and the far 3 stay
    assert detect.remove_redundant(ps, 0.5).tolist() == [False, True, False, True]


# the residual thresholds MATLAB's runners passed to the lowSNR detector (task C02d, probes r8b
# and r8c: their detector calls traced on frames of a synthetic movie with designed noise):
# (camera threshold, the frame's background noise, the threshold)
RUNNER_THRESHOLDS = (
    (3.0, 1.6048004347761791, 2.7000000000000002),
    (3.0, 1.9420178870677436, 2.7000000000000002),
    (3.0, 2.0449257942204269, 3.0),
    (3.0, 2.9559334394370813, 3.0),
    (3.0, 3.0055279532877908, 3.0055279532877908),
    (3.0, 3.2117485606464795, 3.2117485606464795),
    (3.0, 4.4303870183389318, 4.4303870183389318),
    (3.0, 4.5190963828078026, 3.0),
    (3.0, 9.3965835665055764, 3.0),
    (65.0, 1.6048004347761791, 58.5),
    (65.0, 43.006900651041732, 58.5),
    (65.0, 44.509702496520148, 65.0),
    (65.0, 65.9991463408615, 65.9991463408615),
    (65.0, 96.00343051368958, 96.00343051368958),
    (65.0, 98.995456165427413, 65.0),
    (NAN, 1.5166356044723308, 1.5166356044723308),
    (NAN, 16.450795651558572, 16.450795651558572),
)


def test_the_runners_residual_threshold_follows_the_noise():
    """the runners' 'ResidualSigmaThresh' of a frame: 0.9 camera thresholds below 1/1.5 of it,
    the camera threshold up to it, the noise up to 1.5 camera thresholds, the camera threshold
    beyond; the noise without a camera threshold (XR_runDetection3D's)"""
    from petakit5d.detection import run

    for cam, noise, thr in RUNNER_THRESHOLDS:
        got = run.residual_threshold(noise, cam)
        assert got == thr, (cam, noise, got, thr)


def test_detect_frame_gives_the_detector_the_frames_threshold(monkeypatch):
    """detect_frame passes the frame's residual threshold to the detector and the camera
    threshold in 'BackgroundInfo' (MATLAB's runner, probe r8a: [bg noise cam] and the noise of
    a frame noisier than the camera threshold)"""
    from petakit5d.detection import run

    v = _spots_volume("runner_thr", SPOTS, noise=3.6).astype(np.uint16)
    seen = {}
    orig = detect.point_source_detection

    def spy(vol, sigma, **kw):
        seen.update(kw)
        return orig(vol, sigma, **kw)

    monkeypatch.setattr(detect, "point_source_detection", spy)
    detect.clear_threshold_cache()
    run.detect_frame([v], 4, np.array([[1.3, 1.6]]), 2.5, 0, 3.0)
    bi = detect.estimate_background(run.frame_volume(v), 4)
    assert 3.0 < bi[2] <= 4.5
    assert seen["residual_sigma_threshold"] == bi[2]
    assert np.array_equal(seen["background_info"], np.r_[bi, 3.0])


def test_the_runners_go_on_without_a_camera_threshold(tmp_path, monkeypatch):
    """the camera search's errors (no settings file, none with a camera line) give the runners
    no camera threshold, NaN (XR_runDetection3D, probe r8d), and the frames' noise as the
    threshold; the frame-parallel runner the same (MATLAB's stops there: parity D13)"""
    import os

    from petakit5d.detection import run

    (tmp_path / "ex" / "ch1").mkdir(parents=True)
    with open(tmp_path / "ex" / "A_Settings.txt", "w", newline="\n") as f:
        f.write("Laser = 560\n")
    assert np.isnan(run.runner_camera_threshold(str(tmp_path / "ex" / "ch1") + os.sep))

    def no_settings(_source):
        raise run.CameraSettingsError("MATLAB:subsdeldimmismatch", "no settings file")

    monkeypatch.setattr(run, "camera_threshold", no_settings)
    seen = []
    orig = detect.point_source_detection

    def spy(vol, sigma, **kw):
        seen.append((kw["residual_sigma_threshold"], kw["background_info"]))
        return orig(vol, sigma, **kw)

    monkeypatch.setattr(detect, "point_source_detection", spy)
    # the threshold condition's first frame (noise 1.6, its bright spots' residual sigmas about
    # 5): REF's frames are noisier than their spots' residuals, which the noise removes
    make_rt_condition_no_settings(tmp_path / "c")
    data, _ = _rt_data(str(tmp_path / "c"))
    for runner in (run.run_detection, run.run_detection_frame_parallel):
        detect.clear_threshold_cache()
        seen.clear()
        rp = str(tmp_path / runner.__name__) + os.sep
        (entries,) = runner(
            data, sigma=RT_SIGMA, results_path=[rp], detection_method="lowSNR", frames=[1]
        )
        thr, bi = seen[0]
        assert np.isnan(bi[3]) and thr == bi[2]
        assert np.asarray(entries[0]["x"]).size >= 3


def _rt_found(k, thr):
    """frame k of the threshold condition detected with the residual threshold thr: its
    background estimate, whether the weak spot is a detection, how many bright ones are"""
    from petakit5d.detection import run

    detect.clear_threshold_cache()
    V = run.frame_volume(rt_frame(k))
    bi = detect.estimate_background(V, k + 1)
    ps, _ = detect.point_source_detection(
        V, RT_SIGMA[0], background_info=np.r_[bi, 3.0], residual_sigma_threshold=thr
    )
    n = np.size(ps["x"])
    P = np.stack([ps["z"], ps["y"], ps["x"]], axis=1) if n else np.zeros((0, 3))

    def near(q):
        return bool(n) and float(np.min(np.abs(P - np.array(q)).max(axis=1))) < 1

    return bi, near(RT_FRAMES[k][2]), sum(near(q) for q, _ in RT_BRIGHT)


def test_the_threshold_condition_tells_the_rules_apart():
    """the runner_threshold cases' frames: the runners' rule and another rule decide the weak
    spot each way (frame 1: kept with 2.7, not with 3; frame 2: not with the noise 3.55, kept
    with 3; frame 3: kept with 3, not with the noise 5.55), the bright spots kept by both"""
    from petakit5d.detection import run

    for k, other, kept in ((0, 3.0, True), (1, 3.0, False), (2, None, True)):
        bi, _, _ = _rt_found(k, 3.0)
        thr = run.residual_threshold(bi[2], 3.0)
        alt = float(bi[2]) if other is None else other
        assert thr != alt
        _, weak, bright = _rt_found(k, thr)
        _, weak_alt, bright_alt = _rt_found(k, alt)
        assert weak == kept and weak_alt != kept and bright == bright_alt == 3, (k, thr, alt)
