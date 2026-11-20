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
    "conv3fast",
    "est_gaussian_amplitude3d",
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


# -- estGaussianAmplitude3D -----------------------------------------------------------------------


def conv3fast(vol, hy, hx, hz, *, num_threads: int | None = None) -> np.ndarray:
    """MATLAB's ``conv3fast(vol, hy, hx, hz)`` as estGaussianAmplitude3D calls it (an image
    processing mex without source; measured bit for bit): per axis the symmetric kernel given by
    its half ``h[0..w]``, ``h0 v(i) + h1 (v(i-1) + v(i+1)) + ...`` from the centre out, the
    borders mirrored without repeating the edge voxel; the passes along x, then y, then z.
    ``vol`` (z, y, x) double; NaN propagates with the input's payload, as in MATLAB."""
    from petakit5d.core.threads import resolve_threads
    from petakit5d.detection.fit import _ext

    v = np.ascontiguousarray(np.asarray(vol, dtype=np.float64))
    if v.ndim != 3:
        raise ValueError("conv3fast takes a 3-D volume (z, y, x)")
    k = [np.ascontiguousarray(np.asarray(h, dtype=np.float64).ravel()) for h in (hy, hx, hz)]
    out = np.empty_like(v)
    work = np.empty_like(v)
    _ext().conv3fast(v, k[0], k[1], k[2], out, work, resolve_threads(num_threads))
    return out


def est_gaussian_amplitude3d(
    vol, sigma, window=None, *, num_threads: int | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """MATLAB's ``[A, c] = estGaussianAmplitude3D(vol, [s sz], 'WindowSize', [w wz])``: the
    amplitude and background of a Gaussian of sigmas ``(s, s, sz)`` plus a constant fitted by
    least squares in the window around every voxel (black box, measured bit for bit). The
    window's half sizes are ``ceil(2 sigma)`` unless given; the window sums are two
    ``conv3fast`` filterings, of the half kernels ``exp(-(0:w)^2 / (2 s^2))`` and of ones; the
    kernel's sums come from ``exp(-(x^2 + y^2) / (2 s^2)) exp(-z^2 / (2 sz^2))`` on the window,
    in MATLAB's order; ``A = (fg - sum(g) fu / n) / (sum(g^2) - sum(g)^2 / n)``,
    ``c = (fu - A sum(g)) / n``. ``vol`` (z, y, x)."""
    s, sz = (float(t) for t in np.asarray(sigma, dtype=np.float64).ravel()[:2])
    if window is None or np.asarray(window).size == 0:
        w, wz = int(math.ceil(2 * s)), int(math.ceil(2 * sz))
    else:
        ww = np.asarray(window, dtype=np.float64).ravel()
        w, wz = int(ww[0]), int(ww[-1])
    hx = _mexp(-(np.arange(w + 1, dtype=np.float64) ** 2) / (2 * s * s))
    hz = _mexp(-(np.arange(wz + 1, dtype=np.float64) ** 2) / (2 * sz * sz))
    fg = conv3fast(vol, hx, hx, hz, num_threads=num_threads)
    fu = conv3fast(vol, np.ones(w + 1), np.ones(w + 1), np.ones(wz + 1), num_threads=num_threads)
    a = np.arange(-w, w + 1, dtype=np.float64)
    az = np.arange(-wz, wz + 1, dtype=np.float64)
    Z, Y, X = np.meshgrid(az, a, a, indexing="ij")
    g = _mexp(-(X * X + Y * Y) / (2 * s * s)) * _mexp(-Z * Z / (2 * sz * sz))
    n = float(g.size)
    gl = _mlin(g)
    gsum = _f(_r.sum(gl))
    g2sum = _f(_r.sum(gl * gl))
    A = (fg - gsum * fu / n) / (g2sum - gsum * gsum / n)
    c = (fu - A * gsum) / n
    return A, c


# -- local maxima ---------------------------------------------------------------------------------


def locmax3d(vol, window) -> np.ndarray:
    """``XR_locmax3d(vol, [wy wx wz])`` as the detection calls it (odd sizes): with NaN replaced
    by -1e8, a voxel is a maximum when it is larger than every other voxel of its window (the
    window clipped at the volume's faces; a tie is no maximum); the result is the replaced
    volume times the maxima's mask (so -0 at NaN voxels and at negative non-maxima, as MATLAB).
    C++, threaded."""
    from petakit5d.core.threads import resolve_threads
    from petakit5d.detection.fit import _ext

    a = np.ascontiguousarray(np.asarray(vol, dtype=np.float64))
    if a.ndim != 3:
        raise ValueError("locmax3d takes a 3-D volume (z, y, x)")
    w = np.asarray(window, dtype=np.int64).ravel()
    if w.size == 1:
        w = np.repeat(w, 3)
    if np.any(w % 2 == 0):
        raise ValueError("locmax3d takes odd window sizes")
    out = np.empty_like(a)
    # MATLAB's [y x z] in (z, y, x)
    _ext().locmax3d(a, int(w[2]) // 2, int(w[0]) // 2, int(w[1]) // 2, out, resolve_threads(None))
    return out


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
        v = background_corrected(v)
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
    k_level = _f(_st.norminv(1 - alpha))
    A, c, sres, T, df = _prefit_stats(fg, fu, fu2, gsum, g2sum, n, k_level, 1.0, num_threads)
    dfv = df[np.isfinite(df)]
    tthr = t_threshold(alpha, float(dfv.ravel()[0]) if dfv.size else float("nan"))
    chi = _f(_st.chi2inv(1 - 2 * alpha, n))
    return Prefit(v, A, c, sres, T, df, flog, fu, fu2, k_level, tthr, chi, n, wxy, wz)


def _prefit_stats(fg, fu, fu2, gsum: float, g2sum: float, n: int, k_level: float, gfac: float, nt):
    """the voxel-wise model statistics (C++, threaded, each voxel's arithmetic in the order of
    ``A = (fg - gsum fu / n) / (g2sum - gsum^2 / n)``, ``c = (fu - A gsum) / n``,
    ``RSS = A A g2sum - 2 A (fg - c gsum) + fu2 + n c c - 2 c fu`` (at least 0),
    ``sres = sqrt(RSS / (n - 1))``, ``sigma_A = sqrt(RSS / (n - 3) C11)``,
    ``SE = sres / sqrt(2 (n - 1)) k``, ``T = (A - sres k) / sqrt((sigma_A^2 + SE^2) / n)``,
    ``df = (n - 1) (sigma_A^2 + SE^2)^2 / (sigma_A^4 + SE^4)``); ``gfac`` multiplies ``fg``
    first (the 'speedup' detector's normalized filter); returns A, c, sres, T, df"""
    from petakit5d.core.threads import resolve_threads
    from petakit5d.detection.fit import _ext

    J = np.array([[g2sum, gsum], [gsum, float(n)]])
    C11 = float(np.linalg.inv(J)[0, 0])
    shp = np.shape(fg)
    out = [np.empty(shp) for _ in range(5)]
    _ext().prefit_stats(
        np.ascontiguousarray(fg, dtype=np.float64),
        np.ascontiguousarray(fu, dtype=np.float64),
        np.ascontiguousarray(fu2, dtype=np.float64),
        float(gsum),
        float(g2sum),
        float(n),
        C11,
        float(k_level),
        math.sqrt(2 * (n - 1)),
        float(gfac),
        *out,
        resolve_threads(nt),
    )
    return tuple(out)


def lowsnr_mask(P: Prefit) -> np.ndarray:
    """the candidates' mask: the amplitude test ``T > -t`` or the window's chi-square test
    against the constant model ``sum((v - c)^2) / sres^2 > chi2inv(1 - 2 alpha, n)``, without
    the voxels within ``w_xy - 1`` of a face (along z too; measured with w_xy = 3 and w_z = 3
    and 4). ``S = fu2 - 2 c fu + n c c`` per voxel (C++, threaded)."""
    from petakit5d.core.threads import resolve_threads
    from petakit5d.detection.fit import _ext

    m = np.empty(P.T.shape, dtype=np.uint8)
    _ext().lowsnr_mask(
        P.T,
        P.c,
        P.fu,
        P.fu2,
        P.sres,
        float(P.n),
        float(P.t_threshold),
        float(P.chi_threshold),
        max(int(P.wxy) - 1, 0),
        m,
        resolve_threads(None),
    )
    return m.view(bool)


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
    K = int(cc.num_objects)
    labels = np.zeros(mask.size, dtype=np.int64)
    if K:
        # C flat indices; component k's 1-based subscripts sum exactly (integers), so MATLAB's
        # mean (its sum over the count) is the exact sum divided once by the count
        sizes = np.array([np.size(ix) for ix in cc.pixel_idx_list], dtype=np.int64)
        allix = np.concatenate([np.asarray(ix, dtype=np.int64).ravel() for ix in cc.pixel_idx_list])
        lab = np.repeat(np.arange(1, K + 1, dtype=np.int64), sizes)
        labels[allix] = lab
        cnt = np.bincount(lab, minlength=K + 1)[1:].astype(np.float64)
        subs = (allix // (ny * nx), (allix // nx) % ny, allix % nx)  # z, y, x
        means = np.stack(
            [
                np.bincount(lab, weights=(v + 1).astype(np.float64), minlength=K + 1)[1:] / cnt - 1
                for v in subs
            ],
            axis=1,
        )
    else:
        means = np.empty((0, 3))
    out = (_mround(means + 1) - 1).astype(np.int64)
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


# -- the 'speedup' detector (XR_pointSourceDetection3D) -----------------------------------------


def fspecial3_gauss(size, sigma) -> np.ndarray:
    """MATLAB's ``fspecial3('gaussian', size, sigma)`` for odd sizes ``[y x z]`` and sigma
    ``[sy sx sz]`` (or one value): ``exp(-(x^2/(2 sx^2) + y^2/(2 sy^2) + z^2/(2 sz^2)))`` over its
    sum (MATLAB's exp and sum; bit for bit on 12 traced calls); (z, y, x)"""
    sz_ = np.asarray(size, dtype=np.int64).ravel()
    if sz_.size == 1:
        sz_ = np.repeat(sz_, 3)
    sg = np.asarray(sigma, dtype=np.float64).ravel()
    if sg.size == 1:
        sg = np.repeat(sg, 3)
    if np.any(sz_ % 2 == 0):
        raise ValueError("fspecial3_gauss takes odd sizes")
    hy, hx, hz = (int(v) // 2 for v in sz_)
    sy, sx, szz = (float(v) for v in sg)
    Z, Y, X = np.meshgrid(
        np.arange(-hz, hz + 1.0), np.arange(-hy, hy + 1.0), np.arange(-hx, hx + 1.0), indexing="ij"
    )
    hg = _mexp(-(X * X / (2 * sx * sx) + Y * Y / (2 * sy * sy) + Z * Z / (2 * szz * szz)))
    return hg / _f(_r.sum(_mlin(hg)))


_T_THRESHOLDS_SPEEDUP: dict[float, float] = {}


def t_threshold_speedup(alpha: float, df: float) -> float:
    """the 'speedup' detector's mask threshold ``tinv(alpha, df)``, cached per ``alpha`` for the
    session as MATLAB caches it (its own cache, apart from the lowSNR detector's)"""
    from petakit5d.detection import stats as _st

    a = float(alpha)
    if a not in _T_THRESHOLDS_SPEEDUP:
        _T_THRESHOLDS_SPEEDUP[a] = _f(_st.tinv(a, df))
    return _T_THRESHOLDS_SPEEDUP[a]


def background_corrected(vol) -> np.ndarray:
    """the detectors' background correction: each z plane minus its median plus the volume's
    median (``nanmedian``: NaN left out)"""
    v = np.asarray(vol, dtype=np.float64)
    med_all = _f(_ms.median(v.ravel(), omitnan=True))
    med_z = np.array([_f(_ms.median(v[k].ravel(), omitnan=True)) for k in range(v.shape[0])])
    return v - med_z[:, None, None] + med_all


def speedup_prefit(
    vol,
    sigma,
    *,
    alpha: float = 0.05,
    window=None,
    background_correction: bool = False,
    num_threads: int | None = None,
) -> Prefit:
    """the 'speedup' detector's voxel-wise model fits: the Gaussian's window sums through the
    normalized ``fspecial3('gaussian')`` filter times the unnormalized kernel's sum ``G``,
    ``A = (fg G - G fu / n) / (sum(g^2) - G^2 / n)``, ``c = (fu - A G) / n``, the residuals'
    statistics and the amplitude's t statistic with ``k_level = norminv(1 - alpha/2)``; the
    t threshold ``tinv(alpha, mean(df))`` over the finite voxels, once per session"""
    from petakit5d.detection import stats as _st

    v = np.asarray(vol, dtype=np.float64)
    if background_correction:
        v = background_corrected(v)
    if window is None or np.asarray(window).size == 0:
        wxy, wz = window_half_size(sigma)
    else:
        ww = np.asarray(window, dtype=np.int64).ravel()
        wxy, wz = int(ww[0]), int(ww[-1])
    s, sz = float(sigma[0]), float(sigma[1])
    size = (2 * wxy + 1, 2 * wxy + 1, 2 * wz + 1)
    h = fspecial3_gauss(size, (s, s, sz))
    fg = _mf.imfilter(v, h, "conv", "same", "symmetric", num_threads=num_threads)
    u = np.ones_like(h)
    fu = _mf.imfilter(v, u, "conv", "same", "symmetric", num_threads=num_threads)
    fu2 = _mf.imfilter(v * v, u, "conv", "same", "symmetric", num_threads=num_threads)
    klog = fspecial3_log(size, (s, s, sz))
    flog = _mf.imfilter(v, -klog, "conv", "same", "symmetric", num_threads=num_threads)
    g = gaussian_kernel3((s, sz), wxy, wz)
    n = g.size
    gl = _mlin(g)
    gsum = _f(_r.sum(gl))
    g2sum = _f(_r.sum(gl * gl))
    k_level = _f(_st.norminv(1 - alpha / 2))
    A, c, sres, T, df = _prefit_stats(fg, fu, fu2, gsum, g2sum, n, k_level, gsum, num_threads)
    dfl = _mlin(df)
    dfl = dfl[np.isfinite(dfl)]
    tthr = t_threshold_speedup(alpha, _f(_r.mean(dfl)) if dfl.size else float("nan"))
    return Prefit(v, A, c, sres, T, df, flog, fu, fu2, k_level, tthr, float("nan"), n, wxy, wz)


def speedup_mask(P: Prefit) -> np.ndarray:
    """the 'speedup' detector's mask: ``T > -tinv(alpha, df)``, without the voxels within
    ``w_xy - 1`` of a face"""
    with np.errstate(invalid="ignore"):
        m = P.T > -P.t_threshold
    b = P.wxy - 1
    if b > 0:
        m[:b] = False
        m[-b:] = False
        m[:, :b] = False
        m[:, -b:] = False
        m[:, :, :b] = False
        m[:, :, -b:] = False
    return m


def point_source_detection_speedup(
    vol,
    sigma,
    *,
    mode: str = "xyzAc",
    alpha: float = 0.05,
    remove_redundant_fits: bool = True,
    redundancy_radius: float = 0.25,
    window=None,
    fit_method: str = "CD_Newton",
    background_correction: bool = False,
    num_threads: int | None = None,
) -> tuple[dict, np.ndarray]:
    """``XR_pointSourceDetection3D(vol, [s sz], ...)`` (the 'speedup' detector, also the PSF
    sigma estimate's; measured as a black box): the model fits in every voxel's window, the mask
    ``T > -tinv(alpha, df)``, the starts at the Laplacian-of-Gaussian maxima (windows
    ``2 ceil(sigma) + 1``) inside the mask, ``XR_fitGaussians3D`` (mode, method ``'CD_Newton'`` or
    ``'original'``: the fitGaussian3D mex), then the NaN fits, the redundant ones (radius) and
    those whose amplitude is not significant dropped; ``isPSF = ~hval_AD``. Positions 0-based."""
    from petakit5d.detection import fit as _fit

    P = speedup_prefit(
        vol,
        sigma,
        alpha=alpha,
        window=window,
        background_correction=background_correction,
        num_threads=num_threads,
    )
    mask = speedup_mask(P)
    s, sz = float(sigma[0]), float(sigma[1])
    lw = (2 * math.ceil(s) + 1, 2 * math.ceil(s) + 1, 2 * math.ceil(sz) + 1)
    lm = locmax3d(P.flog, lw)
    nz, ny, nx = mask.shape
    pts = np.argwhere((lm != 0) & mask)  # (z, y, x)
    pts = pts[np.argsort(pts[:, 1] + ny * (pts[:, 2] + nx * pts[:, 0]), kind="stable")]
    z, y, x = pts.T
    xyz = np.stack([x, y, z], axis=1).astype(np.float64)
    ps = _fit.fit_gaussians3d(
        P.vol,
        xyz,
        P.A[z, y, x],
        (s, sz),
        P.c[z, y, x],
        mode,
        mask=mask,
        alpha=alpha,
        window=(P.wxy, P.wz),
        method=fit_method,
        num_threads=num_threads,
    )
    ps = postfit(ps, P.vol, remove_redundant_fits=remove_redundant_fits, radius=redundancy_radius)
    return ps, mask


def clear_speedup_cache() -> None:
    """forget the 'speedup' detector's cached t thresholds (a new MATLAB session)"""
    _T_THRESHOLDS_SPEEDUP.clear()


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
