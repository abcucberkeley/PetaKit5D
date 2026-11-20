"""The 3-D point-source detection (task C02): the background estimate, the local maxima, the LoG
kernel, the Anderson-Darling test, the MathWorks functions detection reaches (histcounts,
findpeaks, medfilt1, rangesearch) and the lowSNR detector, against MATLAB R2024b through the
oracle harness, and the detection's rules on synthetic spots (MATLAB's PointDetectionTest in
intent: spots found within 0.25 voxel with amplitudes within 10 %, none in noise).

The detector cases reuse C01's MATLAB references (``test_references.py``: the same MATLAB code
and inputs, so the same keys) and add the Python side.
"""

import importlib.util
from pathlib import Path

import numpy as np

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.detection import detect, fit
from petakit5d.matlabcompat import neighbors, signal
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.compare import EXACT, tol

_spec = importlib.util.spec_from_file_location(
    "detection_references", Path(__file__).with_name("test_references.py")
)
REF = importlib.util.module_from_spec(_spec)
# its cases register under this name, apart from the module pytest collects
REF.__file__ = str(_spec.origin) + ":helpers"
_spec.loader.exec_module(REF)

TOL_LOG = tol(reason="parity D2: fspecial3('log') arithmetic order", abs=1e-16, rel=2e-15)
TOL_DETECT = tol(
    reason="parity D4: the detector's fits (3-D filter order, covariance)", rel=1e-9, abs=1e-9
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
    return cs


CASES = _cases()


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


def test_redundant_fits_keep_the_lowest_rss_of_each_ball():
    ps = {
        "x": np.array([1.0, 1.2, 1.4, 5.0]),
        "y": np.zeros(4),
        "z": np.zeros(4),
        "RSS": np.array([3.0, 1.0, 2.0, 9.0]),
    }
    # 0's ball holds 1 (a lower RSS), 2's ball holds 1 too: only 1 and the far 3 stay
    assert detect.remove_redundant(ps, 0.5).tolist() == [False, True, False, True]
