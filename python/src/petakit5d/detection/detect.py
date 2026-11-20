"""3-D point-source detection as ``XR_runDetection3D`` runs it (method ``'lowSNR'``, the demos'):
the overall background estimate, the model-based candidates and their mask, the Gaussian fits
with their significance tests, the redundancy and bad-point filters, and the per-frame results
(``Detection3D.mat``'s ``frameInfo``).

A new implementation from the published method (Aguet et al. 2013, Dev. Cell 26: 279-291, the
model-based detection: a Gaussian of fixed width plus a constant background fitted in a window
around every voxel, a test of the amplitude against the residuals' noise, and a chi-square test
of the window against the constant model) and from MATLAB R2024b's behaviour, measured as a
black box (task C02: logging wrappers with full inputs and outputs, MATLAB's profiler and
designed inputs; no GPL or MathWorks source was read). ``docs/detection_scope.md`` lists the
stages; ``docs/parity/detection.md`` the measured differences.

Layout: volumes are ``(z, y, x)`` C arrays. Positions are named ``x``, ``y``, ``z`` in voxels,
0-based (MATLAB's minus 1); ``Detection3D.mat`` keeps MATLAB's 1-based values (converted at the
file, ``frame_info_to_matlab``). The voxel order of window sums is MATLAB's column-major order of
``(y, x, z)``: y fastest, then x, then z.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from petakit5d.matlabcompat import filters as _mf
from petakit5d.matlabcompat import reduce as _r
from petakit5d.matlabcompat import stats as _ms
from petakit5d.matlabcompat.elementary import exp as _mexp

__all__ = [
    "estimate_background",
    "fspecial3_log",
    "gaussian_kernel3",
    "lowsnr_prefit",
    "locmax3d",
    "window_half_size",
]


def _mround(x):
    """MATLAB's round (half away from zero)"""
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def _f(x) -> float:
    """a MATLAB-compat reduction's result as a float"""
    return float(np.asarray(x).ravel()[0])


def _mlin(a: np.ndarray) -> np.ndarray:
    """a (z, y, x) array as MATLAB's column-major vector of (y, x, z)"""
    return np.ascontiguousarray(np.asarray(a).transpose(0, 2, 1)).ravel()


# -- the overall background ----------------------------------------------------------------------


def estimate_background(
    vol, randseed=1, *, half_size=(10, 10, 2), n_blocks: int = 200
) -> np.ndarray:
    """``XR_estimate_overall_background_information(vol, 'randseed', s)``: the medians over 200
    random blocks of the blocks' mean, median and standard deviation, ``[mean, median, std]``.

    MATLAB's steps (measured bit for bit): NaN counts as 0; the blocks lie in the bounding box of
    the nonzero voxels. Each block's centre is drawn with MATLAB's twister seeded with ``s``
    (``rng(s)``; NumPy's legacy ``RandomState(s)`` gives the same doubles): ``randi`` between
    ``half_size`` and ``round(0.3 n)`` (either order) for x, y, z in turn, then ``rand(3, 1)``
    for x, y, z: an axis whose value is above 0.5 takes the centre ``n - 1`` if x's value is
    above 0.5, else ``n`` (1-based; the box's size ``n`` along the axis), instead of its draw. A
    centre on a NaN voxel is drawn again. The block is ``centre -/+ half_size`` (y, x, z) clipped
    to the box; its voxels that are not NaN count, in MATLAB's order (zeros count). ``half_size``
    is ``[y x z]``.
    """
    v = np.asarray(vol, dtype=np.float64)
    if v.ndim != 3:
        raise ValueError("estimate_background takes a 3-D volume (z, y, x)")
    nanm = np.isnan(v)
    v0 = np.where(nanm, 0.0, v)
    nz = np.nonzero(v0)
    if nz[0].size == 0:
        raise ValueError("the volume has no nonzero voxel")
    lo = [int(a.min()) for a in nz]
    hi = [int(a.max()) + 1 for a in nz]
    crop = v0[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]]
    ncrop = nanm[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]]
    # MATLAB's (y, x, z) sizes and half sizes
    n = np.array([crop.shape[1], crop.shape[2], crop.shape[0]], dtype=np.int64)
    hs = np.array(half_size, dtype=np.int64)
    rs = np.random.RandomState(int(randseed))
    M = np.empty((n_blocks, 3))
    for i in range(n_blocks):
        while True:
            c = np.zeros(3, dtype=np.int64)
            for ax in (1, 0, 2):  # x, y, z
                a_, b_ = int(hs[ax]), int(_mround(0.3 * n[ax]))
                lo_, hi_ = min(a_, b_), max(a_, b_)
                c[ax] = lo_ + int(math.floor(rs.random_sample() * (hi_ - lo_ + 1)))
            r = rs.random_sample(3)  # x, y, z
            back = 1 if r[0] > 0.5 else 0
            for ax, rv in ((0, r[1]), (1, r[0]), (2, r[2])):
                if rv > 0.5:
                    c[ax] = n[ax] - back
            if not ncrop[c[2] - 1, c[0] - 1, c[1] - 1]:
                break
        b0 = np.maximum(c - hs, 1) - 1
        b1 = np.minimum(c + hs, n)
        blk = crop[b0[2] : b1[2], b0[0] : b1[0], b0[1] : b1[1]]
        x = _mlin(blk)[~_mlin(ncrop[b0[2] : b1[2], b0[0] : b1[0], b0[1] : b1[1]])]
        M[i, 0] = _f(_r.mean(x))
        M[i, 1] = _f(_ms.median(x))
        M[i, 2] = _f(_ms.std(x))
    return np.array([_f(_ms.median(M[:, j])) for j in range(3)])


# -- kernels -------------------------------------------------------------------------------------


def window_half_size(sigma) -> tuple[int, int]:
    """the detector's default window half sizes ``(w_xy, w_z)`` for ``sigma = [s_xy, s_z]``:
    ``ceil(2 sigma)``"""
    s = np.asarray(sigma, dtype=np.float64).ravel()
    if s.size == 1:
        s = np.array([s[0], s[0]])
    return int(math.ceil(2 * s[0])), int(math.ceil(2 * s[1]))


def _grid(wxy: int, wz: int):
    """meshgrid(-wxy:wxy, -wxy:wxy, -wz:wz) in the (z, y, x) layout: X, Y, Z"""
    ax = np.arange(-wxy, wxy + 1, dtype=np.float64)
    az = np.arange(-wz, wz + 1, dtype=np.float64)
    Z, Y, X = np.meshgrid(az, ax, ax, indexing="ij")
    return X, Y, Z


def gaussian_kernel3(sigma, wxy: int, wz: int) -> np.ndarray:
    """the detector's unnormalised Gaussian ``exp(-(x^2 + y^2)/(2 s^2) - z^2/(2 sz^2))`` on
    ``[-w, w]`` (MATLAB's exp), (z, y, x)"""
    s, sz = float(sigma[0]), float(sigma[1])
    X, Y, Z = _grid(wxy, wz)
    return _mexp(-(X * X + Y * Y) / (2 * s * s) - Z * Z / (2 * sz * sz))


def fspecial3_log(size, sigma) -> np.ndarray:
    """MATLAB's ``fspecial3('log', size, sigma)`` for odd sizes ``[y x z]`` and ``sigma`` ``[sy
    sx sz]`` (or one value): the Laplacian of the Gaussian ``h_g = exp(-(x^2/(2 sx^2) + y^2/(2
    sy^2) + z^2/(2 sz^2)))``, ``h = h_g (x^2/sx^4 + y^2/sy^4 + z^2/sz^4 - (1/sx^2 + 1/sy^2 +
    1/sz^2)) / sum(h_g)``, minus its mean; (z, y, x)."""
    sz_ = np.asarray(size, dtype=np.int64).ravel()
    if sz_.size == 1:
        sz_ = np.repeat(sz_, 3)
    sg = np.asarray(sigma, dtype=np.float64).ravel()
    if sg.size == 1:
        sg = np.repeat(sg, 3)
    if np.any(sz_ % 2 == 0):
        raise ValueError("fspecial3_log takes odd sizes")
    hy, hx, hz = (int(v) // 2 for v in sz_)
    sy, sx, szz = (float(v) for v in sg)
    ay = np.arange(-hy, hy + 1, dtype=np.float64)
    ax_ = np.arange(-hx, hx + 1, dtype=np.float64)
    az = np.arange(-hz, hz + 1, dtype=np.float64)
    Z, Y, X = np.meshgrid(az, ay, ax_, indexing="ij")
    arg = X * X / (2 * sx * sx) + Y * Y / (2 * sy * sy) + Z * Z / (2 * szz * szz)
    hg = _mexp(-arg)
    tot = _f(_r.sum(_mlin(hg)))
    lap = X * X / sx**4 + Y * Y / sy**4 + Z * Z / szz**4 - (1 / sx**2 + 1 / sy**2 + 1 / szz**2)
    h = hg * lap / tot
    return h - _f(_r.sum(_mlin(h))) / h.size


# -- local maxima ---------------------------------------------------------------------------------


def locmax3d(vol, window) -> np.ndarray:
    """``XR_locmax3d(vol, [wy wx wz])`` as the detection calls it (odd sizes): with NaN replaced
    by -1e8, a voxel is a maximum when it is larger than every other voxel of its window (the
    window clipped at the volume's faces; a tie is no maximum); the result is the replaced
    volume times the maxima's mask (so -0 at NaN voxels and at negative non-maxima, as MATLAB)."""
    from scipy import ndimage as ndi

    a = np.asarray(vol, dtype=np.float64)
    w = np.asarray(window, dtype=np.int64).ravel()
    if w.size == 1:
        w = np.repeat(w, 3)
    if np.any(w % 2 == 0):
        raise ValueError("locmax3d takes odd window sizes")
    b = np.where(np.isnan(a), -1e8, a)
    hz, hy, hx = int(w[2]) // 2, int(w[0]) // 2, int(w[1]) // 2  # MATLAB's [y x z] in (z, y, x)
    d = ndi.maximum_filter(
        b, size=(2 * hz + 1, 2 * hy + 1, 2 * hx + 1), mode="constant", cval=-np.inf
    )
    cand = np.nonzero(b == d)
    unique = np.ones(cand[0].size, dtype=bool)
    val = b[cand]
    nz, ny, nx = b.shape
    for dz in range(-hz, hz + 1):
        for dy in range(-hy, hy + 1):
            for dx in range(-hx, hx + 1):
                if dz == 0 and dy == 0 and dx == 0:
                    continue
                zz, yy, xx = cand[0] + dz, cand[1] + dy, cand[2] + dx
                ok = (zz >= 0) & (zz < nz) & (yy >= 0) & (yy < ny) & (xx >= 0) & (xx < nx)
                tie = np.zeros_like(ok)
                tie[ok] = b[zz[ok], yy[ok], xx[ok]] == val[ok]
                unique &= ~tie
    is_max = np.zeros(b.shape, dtype=np.float64)
    is_max[tuple(c[unique] for c in cand)] = 1.0
    return b * is_max


# -- the model-based pre-fit ----------------------------------------------------------------------


@dataclass
class Prefit:
    """the voxel-wise estimates of ``XR_pointSourceDetection3D_3`` (all (z, y, x))"""

    vol: np.ndarray  # the background-corrected volume the fits see
    A: np.ndarray
    c: np.ndarray
    sres: np.ndarray  # the residuals' standard deviation
    T: np.ndarray  # the amplitude's t statistic
    df: np.ndarray
    flog: np.ndarray  # the volume filtered with the Laplacian of the Gaussian
    fu: np.ndarray  # the window sums of the volume
    fu2: np.ndarray  # and of its squares
    k_level: float
    t_threshold: float
    chi_threshold: float
    n: int  # voxels of a window
    wxy: int
    wz: int


_T_THRESHOLDS: dict[float, float] = {}


def t_threshold(alpha: float, df: float) -> float:
    """the mask's t threshold ``tinv(16 alpha, df)``, cached per ``alpha`` as MATLAB caches it
    for a session: the first call's value is kept whatever ``df`` later calls give, so a second
    channel (another sigma, another ``df``) uses the first channel's threshold."""
    from petakit5d.detection import stats as _st

    a = float(alpha)
    if a not in _T_THRESHOLDS:
        _T_THRESHOLDS[a] = _f(_st.tinv(16 * a, df))
    return _T_THRESHOLDS[a]


def clear_threshold_cache() -> None:
    """forget the cached t thresholds (a new MATLAB session)"""
    _T_THRESHOLDS.clear()


def lowsnr_prefit(
    vol,
    sigma,
    *,
    alpha: float = 0.05,
    window=None,
    background_correction: bool = True,
    num_threads: int | None = None,
) -> Prefit:
    """the voxel-wise model fits of the lowSNR detector: background correction, the filters, and
    the amplitude, background, residual and test statistics of a Gaussian plus constant fitted
    in the window around every voxel."""
    from petakit5d.detection import stats as _st

    v = np.asarray(vol, dtype=np.float64)
    if background_correction:
        med_all = _f(_ms.median(v.ravel(), omitnan=True))
        med_z = np.array([_f(_ms.median(v[k].ravel(), omitnan=True)) for k in range(v.shape[0])])
        v = v - med_z[:, None, None] + med_all
    if window is None or np.asarray(window).size == 0:
        wxy, wz = window_half_size(sigma)
    else:
        ww = np.asarray(window, dtype=np.int64).ravel()
        wxy, wz = int(ww[0]), int(ww[-1])
    s, sz = float(sigma[0]), float(sigma[1])
    g = gaussian_kernel3((s, sz), wxy, wz)
    u = np.ones_like(g)
    n = g.size
    fg = _mf.imfilter(v, g, "conv", "same", "symmetric", num_threads=num_threads)
    fu = _mf.imfilter(v, u, "conv", "same", "symmetric", num_threads=num_threads)
    fu2 = _mf.imfilter(v * v, u, "conv", "same", "symmetric", num_threads=num_threads)
    klog = fspecial3_log((2 * wxy + 1, 2 * wxy + 1, 2 * wz + 1), (s, s, sz))
    flog = _mf.imfilter(v, -klog, "conv", "same", "symmetric", num_threads=num_threads)
    gl = _mlin(g)
    gsum = _f(_r.sum(gl))
    g2sum = _f(_r.sum(gl * gl))
    A = (fg - gsum * fu / n) / (g2sum - gsum * gsum / n)
    c = (fu - A * gsum) / n
    with np.errstate(invalid="ignore"):
        RSS = A * A * g2sum - 2 * A * (fg - c * gsum) + fu2 + n * c * c - 2 * c * fu
        RSS = np.where(RSS < 0, 0.0, RSS)
        sres = np.sqrt(RSS / (n - 1))
        J = np.array([[g2sum, gsum], [gsum, float(n)]])
        C11 = float(np.linalg.inv(J)[0, 0])
        sigma_A = np.sqrt(RSS / (n - 3) * C11)
        k_level = _f(_st.norminv(1 - alpha))
        SE = sres / math.sqrt(2 * (n - 1)) * k_level
        scomb = np.sqrt((sigma_A * sigma_A + SE * SE) / n)
        T = (A - sres * k_level) / scomb
        df = (n - 1) * (sigma_A * sigma_A + SE * SE) ** 2 / (sigma_A**4 + SE**4)
    dfv = df[np.isfinite(df)]
    tthr = t_threshold(alpha, float(dfv.ravel()[0]) if dfv.size else float("nan"))
    chi = _f(_st.chi2inv(1 - 2 * alpha, n))
    return Prefit(v, A, c, sres, T, df, flog, fu, fu2, k_level, tthr, chi, n, wxy, wz)


def lowsnr_mask(P: Prefit) -> np.ndarray:
    """the candidates' mask: the amplitude test ``T > -t`` or the window's chi-square test
    against the constant model ``sum((v - c)^2) / sres^2 > chi2inv(1 - 2 alpha, n)``, without
    the voxels within ``w_xy - 1`` of a face (along z too)"""
    with np.errstate(invalid="ignore", divide="ignore"):
        S = P.fu2 - 2 * P.c * P.fu + P.n * P.c * P.c
        m = (P.T > -P.t_threshold) | (S / (P.sres * P.sres) > P.chi_threshold)
    # the border is w_xy - 1 along every axis (measured with w_xy = 3 and w_z = 3 and 4)
    bz = by = P.wxy - 1
    if bz > 0:
        m[:bz] = False
        m[-bz:] = False
    if by > 0:
        m[:, :by] = False
        m[:, -by:] = False
        m[:, :, :by] = False
        m[:, :, -by:] = False
    return m


def lowsnr_candidates(
    P: Prefit, mask: np.ndarray, *, residual_sigma_threshold: float | None = None
) -> np.ndarray:
    """the fits' starting voxels ``(k, 3)`` as ``(z, y, x)`` in MATLAB's linear order (as its
    ``unique`` sorts them): one per mask component (26-connected), the mean of its voxels'
    1-based subscripts rounded, and the local maxima of the LoG-filtered volume in the mask more
    than 2 voxels from their component's mean position; without those whose residual sigma is
    at most 0.95 times ``residual_sigma_threshold`` (``FilterByResidualSigma``). Exact on the
    demo's frame 1 (1,264 of 1,264 starts)."""
    from petakit5d.matlabcompat import label as _lb

    nz, ny, nx = mask.shape
    cc = _lb.bwconncomp(mask)
    out = np.empty((cc.num_objects, 3), dtype=np.int64)
    means = np.empty((cc.num_objects, 3))
    labels = np.zeros(mask.size, dtype=np.int64)
    for k, idx in enumerate(cc.pixel_idx_list):
        idx = np.asarray(idx, dtype=np.int64)  # C flat indices in MATLAB's linear order
        labels[idx] = k + 1
        for j, v in enumerate((idx // (ny * nx), (idx // nx) % ny, idx % nx)):  # z, y, x
            means[k, j] = _f(_r.mean((v + 1).astype(np.float64))) - 1
            out[k, j] = int(_mround(means[k, j] + 1)) - 1
    extra = _extra_maxima(P, mask, labels.reshape(mask.shape), means)
    allc = np.concatenate([out, extra]) if extra.size else out
    lin = allc[:, 1] + ny * (allc[:, 2] + nx * allc[:, 0])
    _, first = np.unique(lin, return_index=True)
    out = allc[first]  # MATLAB's unique: sorted by linear index
    if residual_sigma_threshold is not None:
        z, y, x = out.T
        out = out[P.sres[z, y, x] > 0.95 * residual_sigma_threshold]
    return out


def _extra_maxima(
    P: Prefit, mask: np.ndarray, labels: np.ndarray, means: np.ndarray, radius: float = 2.0
) -> np.ndarray:
    """the local maxima (3 x 3 x 3) of the LoG-filtered volume in the mask that lie more than
    ``radius`` voxels from the mean position of their component (the starts besides the
    centroids; measured: none of 414 excluded ones beyond 1.9996, none of 337 added ones within
    2.0169), ``(k, 3)`` (z, y, x) in MATLAB's linear order"""
    lm = locmax3d(P.flog, (3, 3, 3))
    pts = np.argwhere((lm != 0) & mask)
    if pts.size == 0:
        return np.zeros((0, 3), dtype=np.int64)
    comp = labels[tuple(pts.T)]
    d = np.sqrt(np.sum((pts - means[comp - 1]) ** 2, axis=1))
    pts = pts[d > radius]
    nz, ny, nx = mask.shape
    return pts[np.argsort(pts[:, 1] + ny * (pts[:, 2] + nx * pts[:, 0]), kind="stable")]


# -- after the fits -------------------------------------------------------------------------------

_FIELDS = (
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
)


def _take(ps: dict, sel) -> dict:
    """the entries ``sel`` (a boolean or index array over the fits) of every field"""
    out = {}
    for k, v in ps.items():
        a = np.asarray(v)
        out[k] = a[:, sel] if a.ndim == 2 else a[sel]
    return out


def remove_redundant(ps: dict, radius: float = 0.5) -> np.ndarray:
    """which fits survive the redundancy filter: in the ball of ``radius`` voxels around every
    fit (positions ``[x y z]``), all fits but the one with the lowest RSS are removed"""
    from scipy.spatial import cKDTree

    Q = np.stack([np.asarray(ps[k], dtype=np.float64) for k in ("x", "y", "z")], axis=1)
    rss = np.asarray(ps["RSS"], dtype=np.float64)
    keep = np.ones(Q.shape[0], dtype=bool)
    if Q.shape[0] == 0:
        return keep
    for ball in cKDTree(Q).query_ball_point(Q, radius):
        if len(ball) < 2:
            continue
        b = np.asarray(sorted(ball))
        best = b[int(np.argmin(rss[b]))]
        keep[b[b != best]] = False
    return keep


def postfit(
    ps: dict,
    vol: np.ndarray,
    *,
    remove_redundant_fits: bool = True,
    radius: float = 0.5,
    background_info=None,
    residual_sigma_threshold: float = 3.0,
) -> dict:
    """the detector's steps after the fits: the fits with NaN results dropped, the redundant
    ones (``remove_redundant``), those whose amplitude is not significant (``hval_Ar``), and,
    with ``background_info``, the bad points (``filter_bad_points``); ``isPSF`` added"""
    ps = _take(ps, ~np.isnan(np.asarray(ps["x"], dtype=np.float64)))
    if remove_redundant_fits:
        ps = _take(ps, remove_redundant(ps, radius))
    ps = _take(ps, np.asarray(ps["hval_Ar"], dtype=bool))
    if background_info is not None and np.asarray(ps["x"]).size > 0:
        keep = filter_bad_points(ps, residual_sigma_threshold, background_info, vol)
        ps = _take(ps, keep)
    ps["isPSF"] = ~np.asarray(ps["hval_AD"], dtype=bool)
    return ps


def filter_bad_points(ps: dict, threshold: float, background_info, vol) -> np.ndarray:
    """``filter_bad_points_in_detection(pos, sigma_r, thr, bgInfo, vol, A, c)``: which detections
    to keep (``detection.badpoints``)"""
    from petakit5d.detection import badpoints

    pos = np.stack([np.asarray(ps[k], dtype=np.float64) + 1 for k in ("x", "y", "z")], axis=1)
    keep = badpoints.filter_bad_points(
        pos,
        np.asarray(ps["sigma_r"], dtype=np.float64),
        float(threshold),
        np.asarray(background_info, dtype=np.float64).ravel(),
        vol,
        np.asarray(ps["A"], dtype=np.float64),
        np.asarray(ps["c"], dtype=np.float64),
    )
    return np.asarray(keep, dtype=bool).ravel()


def point_source_detection(
    vol,
    sigma,
    *,
    alpha: float = 0.05,
    background_correction: bool = True,
    filter_by_residual_sigma: bool = True,
    residual_sigma_threshold: float = 3.0,
    background_info=None,
    remove_redundant_fits: bool = True,
    redundancy_radius: float = 0.5,
    window=None,
    num_threads: int | None = None,
) -> tuple[dict, np.ndarray]:
    """``XR_pointSourceDetection3D_3(vol, [s sz], ...)`` (the lowSNR detector, mode 'xyzAc',
    no mixtures): the detections (a dict of MATLAB's fields, positions 0-based) and the mask."""
    from petakit5d.detection import fit as _fit

    P = lowsnr_prefit(
        vol,
        sigma,
        alpha=alpha,
        window=window,
        background_correction=background_correction,
        num_threads=num_threads,
    )
    mask = lowsnr_mask(P)
    thr = residual_sigma_threshold if filter_by_residual_sigma else None
    cand = lowsnr_candidates(P, mask, residual_sigma_threshold=thr)
    z, y, x = cand.T
    xyz = np.stack([x, y, z], axis=1).astype(np.float64)
    ps = _fit.fit_gaussians3d(
        P.vol,
        xyz,
        P.A[z, y, x],
        (float(sigma[0]), float(sigma[1])),
        P.c[z, y, x],
        "xyzAc",
        mask=mask,
        alpha=alpha,
        window=(P.wxy, P.wz),
        k_level=P.k_level,
        num_threads=num_threads,
    )
    ps = postfit(
        ps,
        P.vol,
        remove_redundant_fits=remove_redundant_fits,
        radius=redundancy_radius,
        background_info=background_info,
        residual_sigma_threshold=residual_sigma_threshold,
    )
    return ps, mask
