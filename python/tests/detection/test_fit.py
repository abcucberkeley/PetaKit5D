"""The fitGaussian3D mex's fits (task C02b), estGaussianAmplitude3D and conv3fast, and
XR_fitGaussian3D's CD_Newton fit at the box (task C02d: the box's end at the last valid plane,
the stop test of a fit held there, the standard deviations of a negative entry on the
covariance's diagonal) and its stop test's absolute bound and iteration limit (task C08),
against MATLAB R2024b through the oracle harness, and their rules (no MATLAB).

The mex has no source; the port's fits are written from the method (Levenberg-Marquardt least
squares in Moré's trust-region form) and from the mex's measured behaviour. Its linear modes
land within a few ulps of MATLAB's (parity D4), its nonlinear ones within the noise of the
trust-region path's stopping point (D5). estGaussianAmplitude3D and conv3fast are bit for bit.

The stop test's cases (task C08) are designed weak spots at three scales (the window and the
start's A and c times 1/4, 1 and 4: every RSS of the fit scales exactly, so an absolute bound
on its change stops each scale elsewhere), designs at 1/64 and 1/256 of their values that take
more than 500 iterations (MATLAB's limit is 2,000) and, with the demo data (``demodata``),
windows of track processing's recorded calls on the multiple demo's frames, cut from the frames
as the port deskews them.
"""

import decimal
import functools
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.detection import detect, fit
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.compare import EXACT, tol

TOL_LINEAR = tol(
    reason="parity D4: the fitGaussian3D mex's linear fits (the solver's QR arithmetic)",
    rel_max=1e-13,
    rel=1e-13,
    abs=1e-12,
)
TOL_NONLINEAR = tol(
    reason="parity D5: the fitGaussian3D mex's nonlinear fits (the trust-region path's stop)",
    rel=1e-6,
    abs=1e-6,
)


def _rng(stream: str) -> synth.CounterRNG:
    return synth.CounterRNG(2026, "C02b." + stream)


def _gauss_window(shape, xyz, A, s, sz, c, noise, stream):
    """a MATLAB (y, x, z) window: a Gaussian on 0-based voxel coordinates plus noise"""
    ny, nx, nz = shape
    Y, X, Z = np.meshgrid(*(np.arange(n, dtype=np.float64) for n in shape), indexing="ij")
    x, y, z = xyz
    arg = -((X - x) * (X - x) + (Y - y) * (Y - y)) / (2 * s * s) - (Z - z) * (Z - z) / (2 * sz * sz)
    w = A * synth.stable_exp(arg) + c
    if noise:
        w = w + noise * _rng(stream).normal(shape)
    return w


def _windows():
    out = {
        "w977": _gauss_window((9, 9, 7), (4.3, 3.8, 3.2), 60.0, 1.6, 1.33, 100.0, 4.0, "w977"),
        "w777": _gauss_window((7, 7, 7), (3.1, 2.8, 3.3), 25.0, 1.38, 1.1, 100.0, 3.0, "w777"),
        "w7713": _gauss_window((7, 7, 13), (3.2, 3.4, 6.1), 900.0, 1.35, 2.8, 120.0, 5.0, "w7713"),
        "w5711": _gauss_window((5, 7, 11), (3.4, 2.2, 5.3), 300.0, 1.2, 2.1, 50.0, 2.0, "w5711"),
    }
    w = out["w977"].copy()
    w[0, :, :] = np.nan  # a border of voxels left out, as track processing's windows have
    w[4, 4, 1] = np.nan
    out["w977nan"] = w
    return out


def _ad_window(name):
    """a MATLAB (y, x, z) 7 x 7 x 7 window: 100 plus noise of standard deviation 3, and an
    outlier of X on about 3 % of the voxels; X chosen so that the residuals' A^2 of a mode 'c'
    fit lies near the critical value 2.308 (seed 1, X 10.5: 2.345; seed 0, X 11.75: 2.298)"""
    seed, X = {"ad_above": (1, 10.5), "ad_below": (0, 11.75)}[name]
    shape = (7, 7, 7)
    base = 100.0 + 3.0 * synth.CounterRNG(2026, f"C02c.ad{seed}").normal(shape)
    sel = synth.CounterRNG(2026, f"C02c.adsel{seed}").uniform(shape) < 0.03
    return base + X * sel


def _py_mex(win_m, prm, mode):
    p, ps, C, res = fit.fitgaussian3d(from_matlab(win_m), np.ravel(prm), mode)
    r = {
        "data": to_matlab(res["data"]),
        "hAD": np.array(res["hAD"]),
        "mean": res["mean"],
        "std": res["std"],
        "RSS": res["RSS"],
    }
    return [p[None, :], ps[None, :], C, r]


def _eg_inputs():
    r = _rng("eg")
    v = 100 + 5 * r.normal((20, 24, 16))
    v[9, 11, 7] += 300
    w = 50 + 3 * r.normal((17, 13, 11))
    wn = w.copy()
    wn[4, 0, 0] = np.nan
    wn[10, 6, 5] = np.nan
    return {"v": v, "w": w, "wn": wn}


def _py_eg(vol_m, sigma, *opts):
    window = opts[1] if len(opts) == 2 else None
    A, c = detect.est_gaussian_amplitude3d(from_matlab(vol_m), sigma, window)
    return [to_matlab(A), to_matlab(c)]


def _py_conv3fast(vol_m, h1, h2, h3):
    return to_matlab(detect.conv3fast(from_matlab(vol_m), h1, h2, h3))


# -- XR_fitGaussian3D's CD_Newton fit at the box (task C02d) -------------------------------------
TOL_CD = tol(
    reason="parity D3: XR_fitGaussian3D's CD_Newton fits (the sums' order)", rel=1e-9, abs=1e-9
)
S_CD, SZ_CD = 1.38, 1.10  # the detector's sigmas of the demo's first channel, rounded


def _gauss(shape, centre, s=S_CD, sz=SZ_CD):
    """a (z, y, x) Gaussian of peak 1 at centre (x, y, z), 0-based voxel coordinates"""
    nz, ny, nx = shape
    z, y, x = np.meshgrid(*(np.arange(n, dtype=np.float64) for n in shape), indexing="ij")
    arg = -((x - centre[0]) ** 2 + (y - centre[1]) ** 2) / (2 * s * s) - (z - centre[2]) ** 2 / (
        2 * sz * sz
    )
    return synth.stable_exp(arg)


def cd_window(name):
    """(z, y, x) windows of counts (100 plus a spot and normal noise of 2.5, rounded): a spot in
    the NaN planes at the window's end (``end_nan``: MATLAB's box ends at the last valid plane,
    z 4), one before the first valid plane (``start_nan``: the box starts at the window's first
    voxel all the same)"""
    shape = (7, 7, 7)
    centre, A, nan_z = {
        "end_nan": ((3.2, 2.9, 5.4), 60.0, (5, 6)),
        "start_nan": ((2.8, 3.3, 0.9), 60.0, (0, 1)),
    }[name]
    r = synth.CounterRNG(2026, "C02d.box." + name).normal(shape)
    w = np.round(100.0 + A * _gauss(shape, centre) + 2.5 * r)
    for k in nan_z:
        w[k] = np.nan
    return w


def searched_window(i):
    """(z, y, x) windows of a random search of designs (task C02d): one or two spots near or
    beyond the window's faces, normal noise, rounded; NaN beyond a slanted plane or from a z
    plane on. ``NEGC``: fits that end with a negative entry on the covariance's diagonal;
    ``HELD``: fits held on the box while A and c move (no NaN; the stop test's step there is the
    direction it was held from: two more iterations, A 1e-4 apart)"""
    shape = (7, 7, 7)
    r = synth.CounterRNG(2026, f"C02d.negc.{i}").uniform((16,))
    z, y, x = np.meshgrid(*(np.arange(n, dtype=np.float64) for n in shape), indexing="ij")
    c1 = (-2 + 10 * r[0], -2 + 10 * r[1], -2 + 10 * r[2])
    c2 = (-2 + 10 * r[3], -2 + 10 * r[4], -2 + 10 * r[5])
    A1, A2 = 10 + 120 * r[6], 120 * r[7] * (r[8] > 0.5)
    noise = 1 + 4 * r[9]
    n = synth.CounterRNG(2026, f"C02d.negc.n.{i}").normal(shape)
    w = np.round(100 + A1 * _gauss(shape, c1) + A2 * _gauss(shape, c2) + noise * n)
    kind = int(r[10] * 4)
    if kind == 1:
        w[(z + x * (0.5 + r[11])) < 1 + 4 * r[12]] = np.nan
    elif kind == 2:
        w[(z - x * (0.5 + r[11])) > 4 + 3 * r[12]] = np.nan
    elif kind == 3:
        w[int(r[11] * 3) + 4 :] = np.nan
    return w


NEGC = (139, 786, 834)
HELD = (14490, 19948)


def cd_start(w):
    """the start the detector would give such a window: its centre, the amplitude and the
    background from its values"""
    return np.array(
        [3.0, 3.0, 3.0, float(np.nanmax(w) - np.nanmedian(w)), S_CD, SZ_CD, float(np.nanmedian(w))]
    )


def _grid_m(shape_m):
    """X, Y, Z of a MATLAB (y, x, z) window: its 0-based voxel coordinates"""
    Y, X, Z = np.meshgrid(*(np.arange(n, dtype=np.float64) for n in shape_m), indexing="ij")
    return X, Y, Z


def _py_cd(win_m, X, Y, Z, p0, mode, n=10):
    """XR_fitGaussian3D(win, X, Y, Z, p0, 'xyzAc', 10) as XR_fitGaussians3D calls it (track
    processing calls it without the 10: the same results, task C08): the fit's parameters and the
    standard deviations of x, y, z, A and c"""
    R = fit.fit_cd([np.ascontiguousarray(from_matlab(win_m))], np.reshape(p0, (1, 7)))
    return [R.prm, R.pstd[:, [0, 1, 2, 3, 6]]]


def _cd_cases():
    cs = []
    wins = {f"cd.{k}": cd_window(k) for k in ("end_nan", "start_nan")}
    wins.update({f"cd.negc{i}": searched_window(i) for i in NEGC})
    wins.update({f"cd.held{i}": searched_window(i) for i in HELD})
    for name, w in wins.items():
        wm = to_matlab(w)
        X, Y, Z = _grid_m(wm.shape)
        cs.append(
            Case(
                f"xr_fit.{name.split('.', 1)[1]}",
                matlab="XR_fitGaussian3D",
                args=[In("win"), In("X"), In("Y"), In("Z"), cd_start(w), "xyzAc", 10.0],
                nargout=2,
                python=_py_cd,
                inputs={
                    "win": Value(lambda wm=wm: wm),
                    "X": Value(lambda X=X: X),
                    "Y": Value(lambda Y=Y: Y),
                    "Z": Value(lambda Z=Z: Z),
                },
                compare={0: TOL_CD, 1: TOL_CD},
            )
        )
    return cs


# -- XR_fitGaussian3D's stop test (task C08) ------------------------------------------------------
# MATLAB's fit stops when the RSS drops by less than 1e-7 of the first RSS with a step below
# 1e-4, or changes by less than 1e-6 in absolute value. The designs: weak spots whose fits that
# absolute bound stops where a bound relative to the first RSS would not, and the same windows
# scaled (the window and the start's A and c times a power of 2: every RSS scales exactly, so
# MATLAB's stop moves with the scale)
STOP_DESIGNS = (262, 355, 460)  # 1e-6 stops them at scale 1, an iteration before 1e-10 of rss0
STOP_SCALES = (0.25, 1.0, 4.0)
# MATLAB's iteration limit, 2,000: designs whose values times 1/64 or 1/256 take more than 500
# iterations (the 1e-6 bound no longer stops a slow creep at such values): 1,104 and 2,000 (the
# limit) at 1/64, 1,232 at 1/256 (probe d1: 22 of 29 such designs end at MATLAB's)
ITER_DESIGNS = ((3, 1 / 64), (76, 1 / 64), (572, 1 / 256))


def stop_window(i, scale=1.0):
    """(z, y, x) window of design i of C08's search and its start (the detector's: the window's
    centre, its peak over its median, the median), both scaled: a weak spot (A 1.5 to 8) on a
    background of 100 with normal noise of 1.5 to 4, rounded"""
    shape = (7, 7, 7)
    r = synth.CounterRNG(2026, f"C08.stop.{i}").uniform((6,))
    centre = (1.5 + 4 * r[0], 1.5 + 4 * r[1], 1.5 + 4 * r[2])
    A = 1.5 + 6.5 * r[3]
    noise = 1.5 + 2.5 * r[4]
    n = synth.CounterRNG(2026, f"C08.stop.n.{i}").normal(shape)
    w = np.round(100.0 + A * _gauss(shape, centre) + noise * n)
    p0 = np.array(
        [3.0, 3.0, 3.0, float(np.max(w) - np.median(w)), S_CD, SZ_CD, float(np.median(w))]
    )
    p0[[3, 6]] *= scale
    return w * scale, p0


# windows of track processing's XR_fitGaussian3D calls in the multiple demo (recorded from MATLAB
# by task C04c; cut from the deskewed frames without NaN voxels): (experiment, channel, frame,
# the window's first voxel (x, y, z), 1-based, its size (z, y, x), the start as recorded): of
# start buffers, gap frames and an end buffer, both channels. The absolute bound stops the first,
# third and sixth fits an iteration earlier than a bound of 1e-10 of the first RSS did, the
# others an iteration later
TP_WINDOWS = (
    (1, 1, 1, (154, 383, 17), (7, 7, 7),
     (2.8614374213620124, 3.4149086460757303, 2.863370398251746, -0.8602656859449797,
      1.3758153893715306, 1.1049523723652759, 102.11125466482194)),
    (1, 1, 12, (167, 323, 21), (7, 7, 7),
     (2.814540152204245, 3.443159838042334, 3.368981238942112, 4.46869154016484,
      1.3758153893715306, 1.1049523723652759, 121.76418163724342)),
    (2, 2, 6, (158, 162, 10), (7, 9, 9),
     (3.611760172064919, 4.258375287145782, 2.8022037952041856, -1.5554502840260338,
      1.6109238767286678, 1.32923089738702, 101.23864363402612)),
    (2, 1, 6, (208, 89, 8), (7, 7, 7),
     (3.089398218712205, 3.4527074292107187, 3.2721597571034433, 8.5076634341517,
      1.3758153893715306, 1.1049523723652759, 166.58607852207825)),
    (1, 2, 31, (244, 270, 13), (7, 9, 9),
     (3.8012601717815073, 3.7145234084187564, 2.7373922647664344, -11.263406859290317,
      1.6109238767286678, 1.32923089738702, 126.31350864622874)),
    (2, 2, 21, (277, 226, 14), (7, 9, 9),
     (3.6815542863997166, 3.5061265870119485, 3.182425810789976, 0.6693893528531146,
      1.6109238767286678, 1.32923089738702, 101.24961641907522)),
)  # fmt: skip
DEMO = "$PETAKIT5D_DEMO_DATA/PetaKit5D_point_detection_demo_dataset"


def demo_ds_frame(ex: int, ch: int, k: int) -> np.ndarray:
    """frame k of channel ch of experiment ex of the demo dataset deskewed as the demos' deskewData
    deskews it (single), as double: the frames track processing reads"""
    from petakit5d.deskew_rotate.volume import deskew_volume
    from petakit5d.detection.condition import load_condition_data
    from petakit5d.io import read_tiff

    root = Path(os.path.expandvars(DEMO))
    cond = sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith("."))[0]
    m = load_condition_data(str(cond) + os.sep, ["ch1", "ch2"], ["gfp", "rfp"])[ex - 1]
    raw = np.asarray(read_tiff(m.frame_paths[ch - 1][k - 1]), dtype=np.float64)
    ds = deskew_volume(raw, 31.5, float(m.dz), float(m.pixel_size), False, crop=False)
    return np.asarray(ds.astype(np.float32), dtype=np.float64)


@functools.cache
def tp_window(j) -> np.ndarray:
    """TP_WINDOWS[j]'s window (z, y, x)"""
    ex, ch, k, (x0, y0, z0), (nz, ny, nx), _ = TP_WINDOWS[j]
    v = demo_ds_frame(ex, ch, k)
    return np.ascontiguousarray(v[z0 - 1 : z0 - 1 + nz, y0 - 1 : y0 - 1 + ny, x0 - 1 : x0 - 1 + nx])


def _stop_cases():
    cs = []
    designs = [(i, f, f"stop{i}_x{f:g}") for i in STOP_DESIGNS for f in STOP_SCALES]
    designs += [(i, f, f"iters{i}_x{f:g}") for i, f in ITER_DESIGNS]
    for i, f, name in designs:
        w, p0 = stop_window(i, f)
        wm = to_matlab(w)
        X, Y, Z = _grid_m(wm.shape)
        cs.append(
            Case(
                f"xr_fit.{name}",
                matlab="XR_fitGaussian3D",
                args=[In("win"), In("X"), In("Y"), In("Z"), p0, "xyzAc", 10.0],
                nargout=2,
                python=_py_cd,
                inputs={
                    "win": Value(lambda wm=wm: wm),
                    "X": Value(lambda X=X: X),
                    "Y": Value(lambda Y=Y: Y),
                    "Z": Value(lambda Z=Z: Z),
                },
                compare={0: TOL_CD, 1: TOL_CD},
            )
        )
    for j, (ex, ch, k, _, size, p0) in enumerate(TP_WINDOWS):
        X, Y, Z = _grid_m((size[1], size[2], size[0]))
        cs.append(
            Case(
                f"xr_fit.demo_ex{ex}_ch{ch}_frame{k}_{j}",
                matlab="XR_fitGaussian3D",
                args=[In("win"), In("X"), In("Y"), In("Z"), np.array(p0), "xyzAc"],
                nargout=2,
                python=_py_cd,
                inputs={
                    "win": Value(lambda j=j: to_matlab(tp_window(j))),
                    "X": Value(lambda X=X: X),
                    "Y": Value(lambda Y=Y: Y),
                    "Z": Value(lambda Z=Z: Z),
                },
                compare={0: TOL_CD, 1: TOL_CD},
                marks=("demodata",),
            )
        )
    return cs


def _cases():
    W = _windows()
    cs = []
    linear = [
        ("w977", [4.0, 4.0, 3.0, 50.0, 1.6, 1.33, 95.0]),
        ("w777", [3.0, 3.0, 3.0, 0.0, 1.38, 1.1, 0.0]),
        ("w977nan", [3.6, 4.2, 2.9, 61.0, 1.6, 1.33, 101.0]),
        ("w5711", [3.0, 2.0, 5.0, 250.0, 1.2, 2.1, 40.0]),
    ]
    for name, p0 in linear:
        for mode in ("Ac", "A", "c"):
            cs.append(
                Case(
                    f"mex.{mode}.{name}",
                    matlab="fitGaussian3D",
                    args=[In(name), p0, mode],
                    nargout=4,
                    python=_py_mex,
                    inputs={name: Value(lambda name=name: W[name])},
                    compare={0: TOL_LINEAR, 1: TOL_LINEAR, 2: TOL_LINEAR, 3: TOL_LINEAR},
                )
            )
    nonlinear = [
        ("w977", [4.0, 4.0, 3.0, 50.0, 1.5, 1.2, 95.0]),
        ("w7713", [3.0, 3.0, 6.0, 800.0, 1.2, 2.5, 110.0]),
        ("w5711", [3.0, 2.0, 5.0, 250.0, 1.0, 2.0, 40.0]),
        ("w977nan", [4.0, 4.0, 3.0, 55.0, 1.6, 1.33, 98.0]),
    ]
    for name, p0 in nonlinear:
        for mode in ("xyzAsrc", "xyzAc"):
            cs.append(
                Case(
                    f"mex.{mode}.{name}",
                    matlab="fitGaussian3D",
                    args=[In(name), p0, mode],
                    nargout=4,
                    python=_py_mex,
                    inputs={name: Value(lambda name=name: W[name])},
                    compare={
                        0: TOL_NONLINEAR,
                        1: TOL_NONLINEAR,
                        2: TOL_NONLINEAR,
                        3: TOL_NONLINEAR,
                    },
                )
            )
    # the mex's Anderson-Darling decision near its critical value 2.308: noise windows with a few
    # outliers, fitted in mode 'c' (A 0), whose residuals' A^2 is 2.345 (just above) and 2.298
    # (just below)
    for name in ("ad_above", "ad_below"):
        cs.append(
            Case(
                f"mex.c.{name}",
                matlab="fitGaussian3D",
                args=[In(name), [3.0, 3.0, 3.0, 0.0, 1.3, 1.3, 100.0], "c"],
                nargout=4,
                python=_py_mex,
                inputs={name: Value(lambda name=name: _ad_window(name))},
                compare={0: TOL_LINEAR, 1: TOL_LINEAR, 2: TOL_LINEAR, 3: TOL_LINEAR},
            )
        )
    eg = _eg_inputs()
    for name, sigma, opts in (
        ("v", [1.3758, 1.1], ["WindowSize", None]),
        ("v", [1.6109, 1.3292], ["WindowSize", None]),
        ("v", [2.1, 0.8], ["WindowSize", None]),
        ("w", [1.3758, 1.1], ["WindowSize", [2.0, 3.0]]),
        ("wn", [1.3758, 1.1], ["WindowSize", None]),
    ):
        tag = f"{name}.s{sigma[0]}_{sigma[1]}" + ("" if opts[1] is None else ".w23")
        cs.append(
            Case(
                f"est_gaussian_amplitude3d.{tag}",
                matlab="estGaussianAmplitude3D",
                args=[In(name), sigma, *opts],
                nargout=2,
                python=_py_eg,
                inputs={name: Value(lambda name=name: eg[name])},
                compare=EXACT,
            )
        )
    for name, h1, h3 in (
        ("v", [1.0, 0.7, 0.3, 0.09], [1.0, 0.6, 0.2]),
        ("wn", [1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0]),
    ):
        cs.append(
            Case(
                f"conv3fast.{name}.{len(h1)}_{len(h3)}",
                matlab="conv3fast",
                args=[In(name), h1, h1, h3],
                python=_py_conv3fast,
                inputs={name: Value(lambda name=name: eg[name])},
                compare=EXACT,
            )
        )
    return cs


CASES = _cases() + _cd_cases() + _stop_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_the_cd_box_ends_at_the_last_valid_plane():
    """a spot in a window's NaN end planes: the fit held on the last valid plane (z 4), not the
    window's last (MATLAB's box, measured on its traced fits); a spot before the NaN first
    planes: the fit goes below the first valid plane (the box starts at the window's first
    voxel)"""
    w = cd_window("end_nan")
    R = fit.fit_cd([w], cd_start(w)[None, :])
    assert R.prm[0, 2] == 4.0
    w = cd_window("start_nan")
    R = fit.fit_cd([w], cd_start(w)[None, :])
    assert 0.0 < R.prm[0, 2] < 2.0


def test_the_cd_fits_negative_covariance_diagonal_gives_the_magnitude():
    """the designed windows' fits end with a negative entry on C's diagonal: their standard
    deviations are sqrt(|C_kk| RSS / (n - 6)), finite, as MATLAB's results have them"""
    for i in NEGC:
        w = searched_window(i)
        R = fit.fit_cd([w], cd_start(w)[None, :])
        d = np.diag(R.C[0])
        assert np.any(d < 0), i
        k = np.array([0, 1, 2, 3, 6])
        want = np.sqrt(np.abs(d) * R.rss[0] / (R.n[0] - 6))
        np.testing.assert_allclose(R.pstd[0, k], want, rtol=1e-12)


# -- rules (no MATLAB) ----------------------------------------------------------------------------


def test_cr_exp_is_correctly_rounded():
    decimal.getcontext().prec = 60
    xs = np.concatenate([_rng("exp").normal((1, 2000))[0] * 6, [-0.0, 0.5, -700.0, 700.0, 1e-20]])
    ext = fit._ext()
    for x in xs:
        assert ext.cr_exp(float(x)) == float(decimal.Decimal(float(x)).exp())


def test_the_ad_windows_straddle_the_critical_value():
    """the mex's Anderson-Darling decision is A^2 > 2.308 (measured on 33,686 traced calls): the
    designed windows lie at 2.345 (not normal) and 2.298 (normal), as MATLAB decides them
    (oracle cases mex.c.ad_above and ad_below)"""
    p0 = np.array([[3.0, 3.0, 3.0, 0.0, 1.3, 1.3, 100.0]])
    up = fit.fit_gaussian3d_mex([from_matlab(_ad_window("ad_above"))], p0, "c")
    lo = fit.fit_gaussian3d_mex([from_matlab(_ad_window("ad_below"))], p0, "c")
    assert 2.308 < up.a2[0] < 2.4 and bool(up.h_ad[0]) is True
    assert 2.25 < lo.a2[0] < 2.308 and bool(lo.h_ad[0]) is False


def test_mode_letters():
    assert fit.mex_free("cA").tolist() == fit.mex_free("Ac").tolist() == [3, 6]
    assert fit.mex_free("xyzAsrc").tolist() == list(range(7))
    assert fit.mex_free("ac").tolist() == [6]  # letters of 'xyzAsrc' only, as MATLAB's mex


def test_noiseless_gaussian_is_recovered_in_window_coordinates():
    w = from_matlab(_gauss_window((7, 9, 11), (4.3, 3.6, 5.2), 80.0, 1.3, 2.4, 20.0, 0, ""))
    p, ps, C, res = fit.fitgaussian3d(w, [4.3, 3.6, 3.2, 50, 1.5, 2.0, 10], "xyzAsrc")
    assert np.allclose(p, [4.3, 3.6, 5.2, 80, 1.3, 2.4, 20], rtol=1e-10, atol=1e-10)
    assert res["data"].shape == w.shape
    assert ps.shape == (7,) and C.shape == (7, 7)


def test_fixed_parameters_stay_and_nan_voxels_are_left_out():
    w = from_matlab(_windows()["w977nan"])
    p0 = [3.6, 4.2, 2.9, 61.0, 1.6, 1.33, 101.0]
    p, ps, C, res = fit.fitgaussian3d(w, p0, "Ac")
    assert p[[0, 1, 2, 4, 5]].tolist() == [p0[i] for i in (0, 1, 2, 4, 5)]
    assert np.array_equal(np.isnan(res["data"]), np.isnan(w))
    n = int(np.sum(~np.isnan(w)))
    assert ps[0] ** 2 / C[0, 0] == pytest.approx(res["RSS"] / (n - 3), rel=1e-12)


def test_too_few_voxels_raises_where_matlab_aborts():
    w = np.full((3, 3, 3), np.nan)
    w[1, 1, 1] = 5.0
    with pytest.raises(ValueError) as e:
        fit.fitgaussian3d(w, [1, 1, 1, 1, 1, 1, 0], "Ac")
    assert e.value.identifier == "petakit5d:detection:tooFewPoints"


def test_batched_fits_do_not_depend_on_the_thread_count():
    W = _windows()
    wins = [from_matlab(W[k]) for k in ("w977", "w777", "w977nan", "w5711")] * 3
    p0 = np.tile([3.5, 3.0, 3.0, 40.0, 1.4, 1.3, 90.0], (len(wins), 1))
    a = fit.fit_gaussian3d_mex(wins, p0, "xyzAc", num_threads=1)
    b = fit.fit_gaussian3d_mex(wins, p0, "xyzAc", num_threads=4)
    assert np.array_equal(a.prm, b.prm) and np.array_equal(a.C, b.C)


def _conv_ref(v, hy, hx, hz):
    """conv3fast by its measured rule, in NumPy: mirrored borders, centre-out symmetric sums,
    passes x, y, z"""

    def one(a, h, ax):
        h = np.asarray(h, dtype=np.float64)
        w = h.size - 1
        a = np.moveaxis(a, ax, 0)
        n = a.shape[0]
        idx = np.arange(-w, n + w)
        p2 = 2 * (n - 1)
        j = np.mod(idx, p2) if n > 1 else np.zeros_like(idx)
        j = np.where(j < n, j, p2 - j)
        p = a[j]
        out = h[0] * p[w : w + n]
        for k in range(1, w + 1):
            out = out + h[k] * (p[w - k : w - k + n] + p[w + k : w + k + n])
        return np.moveaxis(out, 0, ax)

    return one(one(one(v, hx, 2), hy, 1), hz, 0)


def test_conv3fast_follows_its_rule_on_small_axes():
    v = _rng("conv").normal((5, 3, 2)) * 10
    out = detect.conv3fast(v, [1.0, 0.5], [1.0, 0.25, 0.1], [1.0, 0.4, 0.2, 0.05])
    assert np.array_equal(out, _conv_ref(v, [1.0, 0.5], [1.0, 0.25, 0.1], [1.0, 0.4, 0.2, 0.05]))


def test_the_cd_stop_test_bounds_the_rss_change_absolutely():
    """the designs' fits stop when the RSS changes by less than 1e-6 in their last iteration,
    with a step above 1e-4 and a change of 1e-10 of the first RSS or more (neither the relative
    test nor a bound relative to the first RSS stops them there); scaled by 4 (every RSS 16 times
    larger) they take more iterations (task C08: MATLAB's ends at the three scales)"""
    for i in STOP_DESIGNS:
        w, p0 = stop_window(i)
        R = fit.fit_cd([w], p0[None, :])
        k = int(R.iterations[0])
        R1 = fit.fit_cd([w], p0[None, :], max_iter=k - 1)
        rss0 = fit.fit_cd([w], p0[None, :], max_iter=0).rss[0]
        d = abs(R1.rss[0] - R.rss[0])
        assert 1e-10 * rss0 <= d < 1e-6, i
        assert np.max(np.abs(R.prm[0, :3] - R1.prm[0, :3])) > 1e-4, i
        w4, p4 = stop_window(i, 4.0)
        assert int(fit.fit_cd([w4], p4[None, :]).iterations[0]) > k, i


def test_the_cd_fit_stops_at_2000_iterations():
    """at 1/64 and 1/256 of their values the designs ITER_DESIGNS take more than 500 iterations,
    the second the 2,000 of MATLAB's limit, so a limit of 500 leaves them elsewhere (task C08:
    MATLAB's ends, oracle cases xr_fit.iters*)"""
    for i, f in ITER_DESIGNS:
        w, p0 = stop_window(i, f)
        R = fit.fit_cd([w], p0[None, :])
        k = int(R.iterations[0])
        assert 500 < k <= 2000 and (k == 2000) == (i == 76), (i, k)
        R5 = fit.fit_cd([w], p0[None, :], max_iter=500)
        assert np.max(np.abs(R5.prm - R.prm)) > 1e-3, i


def test_the_cd_fit_held_on_the_box_goes_on_while_a_and_c_move():
    """the HELD windows: their fits end on the box (a coordinate 0 or 6) and take more than the
    iterations at which the position stopped moving (the stop test's step is the direction the
    fit is held from, not its move, 0 there)"""
    for i in HELD:
        w = searched_window(i)
        R = fit.fit_cd([w], cd_start(w)[None, :])
        assert np.isnan(w).sum() == 0
        assert np.any((R.prm[0, :3] == 0.0) | (R.prm[0, :3] == 6.0)), i
        # one iteration fewer leaves A short of its value
        R1 = fit.fit_cd([w], cd_start(w)[None, :], max_iter=int(R.iterations[0]) - 1)
        assert np.array_equal(R1.prm[0, :3], R.prm[0, :3]) and R1.prm[0, 3] != R.prm[0, 3], i
