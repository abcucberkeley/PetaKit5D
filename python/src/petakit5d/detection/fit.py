"""The Gaussian fits of the 3-D point-source detection: ``XR_fitGaussian3D`` (method
``'CD_Newton'``, mode ``'xyzAc'``), ``XR_fitGaussian3D_Ac`` (amplitude and background at a fixed
position), ``XR_fitGaussians3D`` (the windows around candidates, the fits and their tests) and
``adtest_mex``'s Anderson-Darling test.

Written from the published method (Aguet et al. 2013: a Gaussian of fixed width plus a constant
fitted by least squares in a window, an amplitude test against the residuals' noise, and an
Anderson-Darling test of the residuals' normality) and from MATLAB R2024b's behaviour measured as
a black box (task C02; no GPL or MathWorks source). The fit kernel is C++
(``petakit5d._ext.detection``); ``docs/parity/detection.md`` lists what is not bit for bit.

The API for other stages (track processing reads amplitudes out with the same fits):

- ``fit_cd(windows, p0)``: batches of windows, each a ``(z, y, x)`` array (NaN: voxels left out),
  starting parameters ``[x y z A s sz c]`` in the window's 0-based coordinates; returns
  ``FitCD`` arrays (one row per window).
- ``fit_gaussians3d(vol, xyz, A, sigma, c, ...)``: MATLAB's ``XR_fitGaussians3D`` on a volume.

Positions are named: ``x`` along the array's last axis, ``y`` the middle one, ``z`` the first.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from petakit5d.core.threads import resolve_threads

__all__ = ["AD_CRITICAL", "FitCD", "adtest", "ad_statistic", "fit_cd", "fit_gaussians3d"]

# adtest_mex's critical values of A^2 (case 2: mean known, variance estimated), measured
AD_CRITICAL = {(2, 0.05): 2.308}

_EXT: dict = {}


def _ext():
    if "mod" not in _EXT:
        from petakit5d._ext import load

        _EXT["mod"] = load("detection")
    return _EXT["mod"]


def ad_statistic(x, mu: float = 0.0, sigma: float = 1.0) -> float:
    """the Anderson-Darling statistic ``A^2 = -n - sum((2i - 1)(log F(x_i) + log(1 - F(x_m))))/n``,
    ``m = n + 1 - i``,
    of ``x`` **in its given order** against N(mu, sigma^2), as ``adtest_mex`` computes it (it
    does not sort: a sample out of order gives a large statistic)"""
    v = np.ascontiguousarray(np.asarray(x, dtype=np.float64).ravel())
    return float(_ext().ad_stat(v, float(mu), float(sigma)))


def adtest(x, mu: float = 0.0, sigma: float = 1.0, alpha: float = 0.05, case: int = 2) -> bool:
    """``adtest_mex(x, mu, sigma, alpha, case)``: whether A^2 (in x's order) exceeds the critical
    value (measured for case 2 at alpha 0.05: 2.308)"""
    key = (int(case), float(alpha))
    if key not in AD_CRITICAL:
        raise NotImplementedError(
            f"adtest_mex's critical value for case {case}, alpha {alpha} was not measured"
        )
    return ad_statistic(x, mu, sigma) > AD_CRITICAL[key]


@dataclass
class FitCD:
    """results of ``fit_cd`` (k windows): ``prm`` (k, 7) ``[x y z A s sz c]`` in window
    coordinates, ``pstd`` (k, 7) their standard deviations (0 for the sigmas), ``C`` (k, 5, 5)
    the covariance in (x, y, z, A, c) without the residual factor, ``rss``, ``mean`` and
    ``std`` of the residuals, ``h_ad`` the Anderson-Darling test of the sorted residuals,
    ``a2`` its statistic, ``iterations``, ``n`` the voxels fitted"""

    prm: np.ndarray
    pstd: np.ndarray
    C: np.ndarray
    rss: np.ndarray
    mean: np.ndarray
    std: np.ndarray
    h_ad: np.ndarray
    a2: np.ndarray
    iterations: np.ndarray
    n: np.ndarray


def _pack(windows) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """windows ((z, y, x) arrays) as rows of MATLAB's column-major order, their [ny nx nz]"""
    k = len(windows)
    dims = np.empty((k, 3), dtype=np.int64)
    m = 1
    for i, w in enumerate(windows):
        nz, ny, nx = np.shape(w)
        dims[i] = (ny, nx, nz)
        m = max(m, ny * nx * nz)
    vals = np.full((k, m), np.nan)
    nvalid = np.empty(k, dtype=np.int64)
    for i, w in enumerate(windows):
        v = np.asarray(w, dtype=np.float64).transpose(0, 2, 1).ravel()
        vals[i, : v.size] = v
        nvalid[i] = int(np.count_nonzero(~np.isnan(v)))
    return vals, dims, nvalid


def fit_cd(windows, p0, *, max_iter: int = 500, num_threads: int | None = None) -> FitCD:
    """``XR_fitGaussian3D(win, X, Y, Z, p0, 'xyzAc', ...)`` with method ``'CD_Newton'`` on each
    window (a ``(z, y, x)`` array, NaN for voxels left out; its X, Y, Z the 0-based voxel
    coordinates): position, amplitude and background fitted, sigmas fixed. ``p0`` (k, 7)."""
    p = np.ascontiguousarray(np.asarray(p0, dtype=np.float64).reshape(-1, 7))
    vals, dims, nvalid = _pack(windows)
    if p.shape[0] != vals.shape[0]:
        raise ValueError("one row of p0 per window")
    k = p.shape[0]
    prm = np.empty((k, 7))
    pstd = np.empty((k, 7))
    C = np.empty((k, 25))
    st = np.empty((k, 6))
    _ext().fit_cd(
        vals,
        dims,
        p,
        prm,
        pstd,
        C,
        st,
        int(max_iter),
        AD_CRITICAL[(2, 0.05)],
        resolve_threads(num_threads),
    )
    return FitCD(
        prm,
        pstd,
        C.reshape(k, 5, 5),
        st[:, 0],
        st[:, 1],
        st[:, 2],
        st[:, 3] == 1,
        st[:, 4],
        st[:, 5].astype(np.int64),
        nvalid,
    )


# -- XR_fitGaussians3D -----------------------------------------------------------------------------


def _window_half_size(sigma, window) -> tuple[int, int]:
    if window is None or np.asarray(window).size == 0:
        s = np.asarray(sigma, dtype=np.float64).ravel()
        return int(math.ceil(2 * s[0])), int(math.ceil(2 * s[-1]))
    w = np.asarray(window, dtype=np.int64).ravel()
    return int(w[0]), int(w[-1])


def _mround(x):
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def fit_gaussians3d(
    vol,
    xyz,
    A,
    sigma,
    c,
    mode: str = "xyzAc",
    *,
    mask=None,
    alpha: float = 0.05,
    window=None,
    k_level: float | None = None,
    num_threads: int | None = None,
) -> dict:
    """``XR_fitGaussians3D(vol, X, A, sigma, c, mode, 'mask', mask, 'alpha', alpha,
    'WindowSize', [w wz], 'kLevel', k)``: a fit in the window around every position.

    ``xyz`` (k, 3) the positions ``[x y z]`` (0-based voxels; the window is centred on the
    rounded position and clipped to the volume), ``A`` and ``c`` the starting amplitudes and
    backgrounds, ``sigma = [s sz]``. With a mask, the window's voxels of other mask components
    (26-connected labels) than the position's are left out. Results (MATLAB's field names; positions
    0-based): ``x y z A c`` and their ``_pstd``, ``s`` (2, k), ``x_init ...``, ``sigma_r`` (the
    residuals' standard deviation), ``SE_sigma_r``, ``RSS``, ``pval_Ar`` (the one-sided Welch
    test of A against ``k_level`` sigma_r; MATLAB's default ``norminv(1 - alpha/2)``, the
    detector passes ``norminv(1 - alpha)``), ``hval_Ar``, ``hval_AD``. A fit whose position ends
    ``w`` or more from its start along an axis gives NaN everywhere. Mode ``'Ac'``: amplitude and
    background at the given positions (``fit_ac``)."""
    from petakit5d.detection import stats as _st
    from petakit5d.matlabcompat import label as _lb

    if mode not in ("xyzAc", "Ac"):
        raise NotImplementedError(f"mode {mode!r}: 'xyzAc' and 'Ac' are ported")
    v = np.asarray(vol, dtype=np.float64)
    P = np.asarray(xyz, dtype=np.float64).reshape(-1, 3)
    A = np.asarray(A, dtype=np.float64).ravel()
    c = np.asarray(c, dtype=np.float64).ravel()
    s, sz = (float(a) for a in np.asarray(sigma, dtype=np.float64).ravel()[:2])
    wxy, wz = _window_half_size((s, sz), window)
    if k_level is None:
        k_level = float(np.asarray(_st.norminv(1 - alpha / 2)).ravel()[0])
    labels = None
    if mask is not None:
        labels, _ = _lb.bwlabeln(np.asarray(mask, dtype=bool))
    nz, ny, nx = v.shape
    k = P.shape[0]
    wins = []
    p0 = np.empty((k, 7))
    starts = np.empty((k, 3), dtype=np.int64)  # (x, y, z) 0-based
    # positions in MATLAB's 1-based values, so that the window offsets and the results have its
    # bits (x - 1 and x + 1 are exact for the positions a volume has)
    P1 = P + 1
    for i in range(k):
        cx, cy, cz = (int(_mround(P1[i, j])) - 1 for j in range(3))
        x0, x1 = max(cx - wxy, 0), min(cx + wxy, nx - 1)
        y0, y1 = max(cy - wxy, 0), min(cy + wxy, ny - 1)
        z0, z1 = max(cz - wz, 0), min(cz + wz, nz - 1)
        w = v[z0 : z1 + 1, y0 : y1 + 1, x0 : x1 + 1].copy()
        if labels is not None:
            lab = labels[z0 : z1 + 1, y0 : y1 + 1, x0 : x1 + 1]
            own = labels[cz, cy, cx]
            w[(lab != 0) & (lab != own)] = np.nan
        wins.append(w)
        starts[i] = (x0, y0, z0)
        p0[i] = (P1[i, 0] - (x0 + 1), P1[i, 1] - (y0 + 1), P1[i, 2] - (z0 + 1), A[i], s, sz, c[i])
    R = fit_ac(wins, p0) if mode == "Ac" else fit_cd(wins, p0, num_threads=num_threads)
    return _results(R, P, p0, starts, wxy, wz, s, sz, k_level, alpha, fixed=mode == "Ac")


def fit_ac(windows, p0) -> FitCD:
    """``XR_fitGaussian3D_Ac``: amplitude and background at a fixed position by linear least
    squares (``A = (n sum(g d) - sum(g) sum(d)) / (n sum(g^2) - sum(g)^2)``, ``c = (sum(d) - A
    sum(g)) / n``, MATLAB's sums); ``C = inv(2 J'J)`` in (A, c), ``prmStd = sqrt(diag(C) RSS /
    (n - 3))``; the Anderson-Darling test on the residuals **unsorted**, as MATLAB passes them
    (so it nearly always rejects). Measured on 1,005 traced calls: A within 4e-12, the rest within
    1e-15."""
    from petakit5d.matlabcompat import reduce as _r

    p = np.asarray(p0, dtype=np.float64).reshape(-1, 7)
    k = p.shape[0]
    prm = p.copy()
    pstd = np.zeros((k, 7))
    C = np.full((k, 5, 5), np.nan)
    st = np.full((k, 6), np.nan)
    nvalid = np.zeros(k, dtype=np.int64)
    sm = lambda a: float(np.asarray(_r.sum(a)).ravel()[0])  # noqa: E731
    for i, w in enumerate(windows):
        wl = np.asarray(w, dtype=np.float64).transpose(0, 2, 1).ravel()  # MATLAB's order
        nz_, ny_, nx_ = np.shape(w)
        q = np.arange(wl.size)
        ok = ~np.isnan(wl)
        d = wl[ok]
        yf = (q % ny_)[ok].astype(np.float64)
        xf = ((q // ny_) % nx_)[ok].astype(np.float64)
        zf = (q // (ny_ * nx_))[ok].astype(np.float64)
        s, sz = p[i, 4], p[i, 5]
        from petakit5d.matlabcompat.elementary import exp as _mexp

        g = _mexp(
            -((xf - p[i, 0]) ** 2 + (yf - p[i, 1]) ** 2) / (2 * s * s)
            - (zf - p[i, 2]) ** 2 / (2 * sz * sz)
        )
        n = d.size
        nvalid[i] = n
        if n < 3:
            continue
        sg, sg2, sd, sgd = sm(g), sm(g * g), sm(d), sm(g * d)
        Ai = (n * sgd - sg * sd) / (n * sg2 - sg * sg)
        ci = (sd - Ai * sg) / n
        r = d - Ai * g - ci
        rss = float(np.sum(r * r))
        Ci = np.linalg.inv(2 * np.array([[sg2, sg], [sg, float(n)]]))
        sd_r = math.sqrt(rss / (n - 1))
        prm[i, 3], prm[i, 6] = Ai, ci
        pstd[i, 3] = math.sqrt(Ci[0, 0] * rss / (n - 3))
        pstd[i, 6] = math.sqrt(Ci[1, 1] * rss / (n - 3))
        C[i, 3, 3], C[i, 3, 4], C[i, 4, 3], C[i, 4, 4] = Ci[0, 0], Ci[0, 1], Ci[1, 0], Ci[1, 1]
        a2 = ad_statistic(r, 0.0, sd_r)
        st[i] = (rss, sm(r) / n, sd_r, float(a2 > AD_CRITICAL[(2, 0.05)]), a2, 0)
    return FitCD(
        prm,
        pstd,
        C,
        st[:, 0],
        st[:, 1],
        st[:, 2],
        st[:, 3] == 1,
        st[:, 4],
        np.zeros(k, dtype=np.int64),
        nvalid,
    )


def _results(
    R: FitCD,
    P: np.ndarray,
    p0: np.ndarray,
    starts: np.ndarray,
    wxy: int,
    wz: int,
    s: float,
    sz: float,
    k_level: float,
    alpha: float,
    fixed: bool = False,
) -> dict:
    """XR_fitGaussians3D's fields from the fits (``fixed``: the 'Ac' fits, whose positions stay)"""
    from petakit5d.detection import stats as _st

    k = P.shape[0]
    # MATLAB's x = prm + start (1-based), then 0-based (exact)
    x = (R.prm[:, 0] + (starts[:, 0] + 1)) - 1
    y = (R.prm[:, 1] + (starts[:, 1] + 1)) - 1
    z = (R.prm[:, 2] + (starts[:, 2] + 1)) - 1
    d = np.abs(R.prm[:, :3] - p0[:, :3])
    with np.errstate(invalid="ignore"):
        far = (d[:, 0] >= wxy) | (d[:, 1] >= wxy) | (d[:, 2] >= wz)
    n = R.n.astype(np.float64)
    sigma_r = R.std.copy()
    SE = sigma_r / np.sqrt(2 * (n - 1))
    A = R.prm[:, 3].copy()
    A_pstd = R.pstd[:, 3].copy()
    with np.errstate(invalid="ignore", divide="ignore"):
        T = (A - sigma_r * k_level) / np.sqrt((A_pstd**2 + (SE * k_level) ** 2) / n)
        df = (n - 1) * (A_pstd**2 + (SE * k_level) ** 2) ** 2 / (A_pstd**4 + (SE * k_level) ** 4)
    out = {
        "x": x,
        "y": y,
        "z": z,
        "A": A,
        "s": np.vstack([np.full(k, s), np.full(k, sz)]),
        "c": R.prm[:, 6].copy(),
        "x_pstd": R.pstd[:, 0].copy(),
        "y_pstd": R.pstd[:, 1].copy(),
        "z_pstd": R.pstd[:, 2].copy(),
        "A_pstd": A_pstd,
        "s_pstd": np.zeros((2, k)),
        "c_pstd": R.pstd[:, 6].copy(),
        "x_init": P[:, 0].copy(),
        "y_init": P[:, 1].copy(),
        "z_init": P[:, 2].copy(),
        "sigma_r": sigma_r,
        "SE_sigma_r": SE,
        "RSS": R.rss.copy(),
    }
    pval = np.asarray(_st.tcdf(-T, df), dtype=np.float64).ravel() if k else np.zeros(0)
    out["pval_Ar"] = pval
    out["hval_Ar"] = pval < alpha
    out["hval_AD"] = R.h_ad.copy()
    for key in (
        "x",
        "y",
        "z",
        "A",
        "c",
        "x_pstd",
        "y_pstd",
        "z_pstd",
        "A_pstd",
        "c_pstd",
        "sigma_r",
        "SE_sigma_r",
        "RSS",
        "pval_Ar",
    ):
        out[key] = np.where(far, np.nan, out[key])
    # the tests of a NaN fit are false (MATLAB: 136 of 136)
    out["hval_Ar"] = out["hval_Ar"] & ~far
    out["hval_AD"] = out["hval_AD"] & ~far
    return out
