"""The fitGaussian3D mex's fits (task C02b), estGaussianAmplitude3D and conv3fast, and
XR_fitGaussian3D's CD_Newton fit at the box (task C02d: the box's end at the last valid plane,
the stop test of a fit held there, the standard deviations of a negative entry on the
covariance's diagonal), against MATLAB R2024b through the oracle harness, and their rules (no
MATLAB).

The mex has no source; the port's fits are written from the method (Levenberg-Marquardt least
squares in Moré's trust-region form) and from the mex's measured behaviour. Its linear modes
land within a few ulps of MATLAB's (parity D4), its nonlinear ones within the noise of the
trust-region path's stopping point (D5). estGaussianAmplitude3D and conv3fast are bit for bit.
"""

import decimal

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


def _py_cd(win_m, X, Y, Z, p0, mode, n):
    """XR_fitGaussian3D(win, X, Y, Z, p0, 'xyzAc', 10) as XR_fitGaussians3D calls it: the fit's
    parameters and the standard deviations of x, y, z, A and c"""
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


CASES = _cases() + _cd_cases()


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
    assert e.value.identifier == "pk5d:fitGaussian3D:tooFewPoints"


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
