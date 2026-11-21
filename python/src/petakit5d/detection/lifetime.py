"""Lifetime analysis of processed tracks: MATLAB's ``runLifetimeAnalysis3D`` and
``plotIntensityCohorts`` (task C05), their numbers (the figures are left out).

A new implementation. MATLAB's lifetime analysis (``runLifetimeAnalysis3D``, ``getLifetimeData``,
``loadTracks``, ``scaleEDFs``, ``fitGaussianModeToCDF``, ``getVisitorIndex``, ``adtest1``,
``plotIntensityCohorts`` and the ``binterp`` mex) is GPL code that was not read: every rule below
was measured from MATLAB R2024b as a black box (the inputs and outputs of its calls and of the
MathWorks functions they call, recorded on the demo data, and probes with chosen inputs), with the
published method of Aguet et al. 2013 (Dev. Cell 26:279) as the specification of what is computed.

The data flow, as measured:

- **Per movie** (``get_lifetime_data``): the tracks of
  ``<source><AnalysisPath>/ProcessedTracks.mat``, sorted by lifetime (longest first, ties in file
  order), give ``lifetimeData.mat`` (written there unless it exists and ``overwrite`` is false):
  every track's ``lifetime_s``, ``trackLengths`` (end - start + 1), ``start``, ``catIdx``,
  ``index`` (1..N) and the significance flags; for the category Ia tracks (catIdx 1) their
  amplitudes, their uncertainties, residual noise and its error from frame 1 of the track (NaN
  after), the gap mask, and the start and end buffers. In memory the tracks of at least
  ``cutoff_f`` frames are kept (``*_all`` fields), and the Ia tracks among them; visitors
  (``visitor_index``) are moved to ``visitors``. With ``scale``, each channel's maximum
  amplitudes (of all Ia tracks of each movie) are aligned between movies (``scale_edfs``), and
  the amplitude fields of the kept tracks (not the visitors) multiplied by the movie's factor.
- **runLifetimeAnalysis3D**: lifetime histograms on ``t = (cutoff_f:N)*framerate``,
  ``N = movieLength - 2*buffer`` (buffer 5), weighted by ``N/(N - L + 1)`` (the observation bias of
  a track of L frames) and normalized to unit area; six lifetime cohorts of the kept Ia tracks
  ([1, 10], [11, 15], [16, 20], [21, 40], [41, 60], [61, 120] s); for 3 to 12 frames the mode of
  each cohort's maximum over the first frames is fitted (``fit_gaussian_mode_to_cdf``) and the
  cohorts' modes tested against the first's (``adtest1``); ``nFramesThreshold`` is one less than
  the first number of frames where they differ; the maximum intensity threshold is the 99th
  percentile of the Gaussian mode of the kept tracks' maxima over the first ``nFramesThreshold``
  frames; tracks above it are CCPs (clathrin-coated pits), below CSs.
- **plotIntensityCohorts**: the kept Ia tracks of at least 5 frames in lifetime cohorts
  ([10, 20), [20, 40), ..., [100, 120) s), each track with its buffers resampled (``binterp``) to
  the cohort's length, averaged per movie and then over movies, with the standard error over
  movies, for each combination of the other channel's significance; with one movie, the median
  of its tracks between their 25th and 75th percentiles.

The fits (``matlabcompat.optim.lsqnonlin``) follow MATLAB's path in most cases, and the results
then agree to about 1e-9; where MATLAB takes a step the port does not reproduce, they differ
more (the mode fit's non-smooth objective can end in another local minimum).
``docs/parity/detection.md`` (rows L1 to L7) lists what was measured: on the demo data every
field of ``lftRes`` but the threshold (9.9e-10) is MATLAB's bit for bit. Everything else is
MATLAB's arithmetic in MATLAB's order; MATLAB's quirks are kept and listed in that page. Cell
masks (MATLAB's cell areas and initiation densities) are not supported: those fields are NaN, as
in the demo's runs, and a mask file where ``MaskPath`` names it is refused.

``run_lifetime_analysis`` is registered with MATLAB's interface (names, order, defaults) and the
port's own checks of the options, made from MATLAB's accept and refuse behaviour; it takes
``AnalysisPath``, ``ProcessedTracks``, ``LifetimeData``, ``Overwrite`` and ``Display`` (measured:
the numbers and the two printed values are the same with 'off'), every other option at its
default only (``NotImplementedError`` otherwise). The figures are left out.
"""

from __future__ import annotations

import math
import os
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from petakit5d.core import dtypes
from petakit5d.params import Param, register

__all__ = [
    "LifetimeError",
    "adtest1",
    "binterp",
    "fit_gaussian_mode_to_cdf",
    "get_lifetime_data",
    "hist_weights",
    "lifetime_data_from_tracks",
    "load_tracks",
    "plot_intensity_cohorts",
    "read_processed_tracks",
    "run_lifetime_analysis",
    "scale_edfs",
    "visitor_index",
]

_SQRT2 = math.sqrt(2.0)
# the lifetime analysis' fixed settings, as measured (MATLAB's defaults)
CUTOFF_F = 3  # frames: the shortest tracks kept by runLifetimeAnalysis3D
COHORT_CUTOFF_F = 5  # plotIntensityCohorts' shortest tracks
BUFFER = 5  # frames: the histogram range is cutoff..movieLength - 2*BUFFER
LOADTRACKS_CUTOFF_F = 2  # loadTracks' own cutoff, as getLifetimeData calls it
FRAME_RANGE = tuple(range(3, 13))  # numbers of first frames tested for nFramesThreshold
COHORT_LB = (1.0, 11.0, 16.0, 21.0, 41.0, 61.0)  # s, the threshold cohorts (inclusive)
COHORT_UB = (10.0, 15.0, 20.0, 40.0, 60.0, 120.0)
INTENSITY_BOUNDS = (10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0)  # s, plotIntensityCohorts
VISITOR_LIFETIME = (
    30.0  # s: tracks this long give the visitors' thresholds; shorter ones are tested
)
TRACK_FIELDS = ("lifetime_s", "start", "end", "catIdx", "significantMaster", "significantSlave")
IA_FIELDS = ("A", "A_pstd", "sigma_r", "SE_sigma_r", "gapVect", "startBuffer", "endBuffer")
LIFETIME_FIELDS = (
    "lifetime_s",
    "trackLengths",
    "start",
    "catIdx",
    "index",
    "A",
    "A_pstd",
    "sigma_r",
    "SE_sigma_r",
    "sbA",
    "ebA",
    "sbSigma_r",
    "ebSigma_r",
    "gapMat_Ia",
    "significantMaster",
    "significantSlave",
)
_PER_TRACK = (
    "lifetime_s",
    "trackLengths",
    "start",
    "catIdx",
    "index",
    "significantMaster",
    "significantSlave",
)
_IA_ARRAYS = (
    "A",
    "A_pstd",
    "sigma_r",
    "SE_sigma_r",
    "sbA",
    "ebA",
    "sbSigma_r",
    "ebSigma_r",
    "gapMat_Ia",
)
_SCALED = ("A", "A_pstd", "sigma_r", "SE_sigma_r", "sbA", "ebA", "sbSigma_r", "ebSigma_r")


class LifetimeError(ValueError):
    """An error where MATLAB raises one; ``identifier`` is MATLAB's error id."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


# -- MathWorks functions as the lifetime analysis calls them ---------------------------------------


def _interp1(x, v, xq) -> np.ndarray:
    """MATLAB's ``interp1(x, v, xq)`` ('linear', NaN outside; ``matlabcompat.interp``)"""
    from petakit5d.matlabcompat.interp import interp1

    f64 = np.float64
    return interp1(np.asarray(x, dtype=f64), np.asarray(v, dtype=f64), np.asarray(xq, dtype=f64))


def _linspace(a: float, b: float, n: int) -> np.ndarray:
    """MATLAB's ``linspace(a, b, n)``: ``a + (0:n-1)*(b - a)/(n - 1)``, the last point b"""
    y = a + np.arange(n, dtype=np.float64) * (b - a) / (n - 1)
    y[-1] = b
    return y


def _colon(a: float, step: float, b: float) -> np.ndarray:
    """MATLAB's ``a:step:b`` for exact cases (the probability grid 0:0.001:1)"""
    n = int(math.floor((b - a) / step + 0.5)) + 1
    k = np.arange(n, dtype=np.float64)
    half = n // 2
    y = np.empty(n)
    y[:half] = a + k[:half] * step
    y[half:] = b - (n - 1 - k[half:]) * step
    return y


def _msum(v) -> float:
    from petakit5d.matlabcompat import reduce as mr

    return float(np.asarray(mr.sum(np.asarray(v, dtype=np.float64))).ravel()[0])


def _nanmax_rows(X: np.ndarray) -> np.ndarray:
    """MATLAB's ``max(X, [], 2)`` (and nanmax): NaN left out, the row's first NaN where it has
    no value"""
    from petakit5d.matlabcompat import reduce as mr

    if X.shape[0] == 0:
        return np.zeros(0)
    if X.shape[1] == 0:
        return np.full(X.shape[0], dtypes.NAN)
    return np.asarray(mr.max(X, 1)[0], dtype=np.float64)


def _ecdf(y):
    from petakit5d.detection import stats

    v = np.asarray(y, dtype=np.float64).ravel()
    if not np.any(~np.isnan(v)):
        raise LifetimeError(
            "stats:ecdf:NotEnoughData", "ecdf of no values (a category or cohort without tracks)"
        )
    return stats.ecdf(v)


def hist(y, centers) -> np.ndarray:
    """MATLAB's ``hist(y, centers)`` of a vector (``matlabcompat.hist``)"""
    from petakit5d.matlabcompat.hist import hist as _h

    return _h(y, centers)


# -- reading tracks --------------------------------------------------------------------------------


_LOGICAL = frozenset({"significantMaster", "significantSlave", "gapVect"})
BUFFER_FIELDS = ("A", "sigma_r")  # the buffers' fields the lifetime analysis reads


def _h5_array(oid, logical: bool) -> np.ndarray:
    """a MATLAB v7.3 dataset (by its low-level id) as an array of MATLAB's shape: an empty array
    (a uint64 vector of its dimensions with the attribute MATLAB_empty) as zeros of that size"""
    import h5py

    shape = oid.get_space().shape
    arr = np.empty(shape, dtype=oid.dtype)
    if arr.size:
        oid.read(h5py.h5s.ALL, h5py.h5s.ALL, arr)
    if arr.dtype == np.uint64 and arr.ndim == 1 and h5py.h5a.exists(oid, b"MATLAB_empty"):
        return np.zeros(tuple(int(v) for v in arr), dtype=bool if logical else np.float64)
    if logical:
        arr = arr.astype(bool)
    elif arr.dtype != np.float64:
        arr = arr.astype(np.float64)
    return arr.T if arr.ndim >= 2 else arr


def _h5_value(fid, ref, name: str):
    """one track's field of a v7.3 file by its object reference; a struct (the buffers) as a
    dict of its fields ``BUFFER_FIELDS``"""
    import h5py

    oid = h5py.h5r.dereference(ref, fid)
    if isinstance(oid, h5py.h5g.GroupID):
        out = {}
        for b in BUFFER_FIELDS:
            try:
                out[b] = _h5_array(h5py.h5d.open(oid, b.encode()), False)
            except KeyError:
                pass
        return out
    return _h5_array(oid, name in _LOGICAL)


def read_processed_tracks(path: str, fields: Sequence[str], ia_fields: Sequence[str] = ()) -> dict:
    """The fields of the ``tracks`` struct array of a ProcessedTracks.mat (v7.3 through h5py's
    low-level interface, v7 and v6 through ``petakit5d.io.matfile``), as lists in track order:
    ``fields`` for every track, ``ia_fields`` only for the category Ia tracks (catIdx 1; None for
    the others). Arrays have MATLAB's shapes (channels x frames); structs (the buffers) are
    dicts."""
    from petakit5d.io import matfile

    want = list(dict.fromkeys(list(fields) + ["catIdx"]))
    if matfile.mat_version(path) == "7.3":
        import h5py

        with h5py.File(path, "r") as f:
            g = f["tracks"]
            fid = f.id
            out = {}
            for k in want:
                refs = g[k][()].ravel(order="F")
                out[k] = [_h5_value(fid, r, k) for r in refs]
            cat = [float(v.ravel()[0]) for v in out["catIdx"]]
            for k in ia_fields:
                refs = g[k][()].ravel(order="F")
                out[k] = [
                    _h5_value(fid, r, k) if c == 1 else None for r, c in zip(refs, cat, strict=True)
                ]
        return out
    tr = matfile.read_mat(path, "tracks")["tracks"]
    els = [tr] if isinstance(tr, dict) else list(tr)
    out = {k: [e[k] for e in els] for k in want}
    cat = [float(np.asarray(v).ravel()[0]) for v in out["catIdx"]]
    for k in ia_fields:
        out[k] = [e[k] if c == 1 else None for e, c in zip(els, cat, strict=True)]
    return out


def _scalar(v) -> float:
    return float(np.asarray(v, dtype=np.float64).ravel()[0])


def load_tracks(path: str, *, cutoff_f: int = LOADTRACKS_CUTOFF_F) -> dict:
    """MATLAB's ``loadTracks`` as getLifetimeData calls it ('Category' 'all', 'Mask' false): the
    tracks of at least ``cutoff_f`` frames, sorted by lifetime, longest first (ties in file
    order; measured: track processing writes them in that order already)."""
    tr = read_processed_tracks(path, TRACK_FIELDS, IA_FIELDS)
    n = len(tr["catIdx"])
    lengths = [_scalar(tr["end"][i]) - _scalar(tr["start"][i]) + 1 for i in range(n)]
    keep = [i for i in range(n) if lengths[i] >= cutoff_f]
    life = np.array([_scalar(tr["lifetime_s"][i]) for i in keep], dtype=np.float64)
    order = [keep[j] for j in np.argsort(-life, kind="stable")]
    if order != list(range(n)):
        tr = {k: [v[i] for i in order] for k, v in tr.items()}
    return tr


def lifetime_data_from_tracks(tracks: dict, movie_length: int, n_channels: int) -> dict:
    """The variables of MATLAB's ``lifetimeData.mat`` (one movie; module docstring), in its order
    and with its classes and shapes (N x 1 doubles, N x nCh logicals, nIa x movieLength x nCh)."""
    n = len(tracks["catIdx"])
    lifetime = np.array([_scalar(v) for v in tracks["lifetime_s"]], dtype=np.float64)
    start = np.array([_scalar(v) for v in tracks["start"]], dtype=np.float64)
    end = np.array([_scalar(v) for v in tracks["end"]], dtype=np.float64)
    cat = np.array([_scalar(v) for v in tracks["catIdx"]], dtype=np.float64)

    def flags(name):
        out = np.zeros((n, n_channels), dtype=bool)
        for i, v in enumerate(tracks[name]):
            b = np.asarray(v, dtype=bool).ravel()
            out[i, : b.size] = b
        return out

    ia = np.flatnonzero(cat == 1)
    nia = ia.size
    nb = 3
    for i in ia:
        b = tracks["startBuffer"][i]
        if isinstance(b, dict) and "A" in b:
            nb = np.asarray(b["A"]).reshape(n_channels, -1).shape[1]
            break
    nan = dtypes.NAN
    arrays = {
        k: np.full((nia, movie_length, n_channels), nan)
        for k in ("A", "A_pstd", "sigma_r", "SE_sigma_r")
    }
    bufs = {
        k: np.full((nia, nb, n_channels), nan) for k in ("sbA", "ebA", "sbSigma_r", "ebSigma_r")
    }
    gaps = np.zeros((nia, movie_length), dtype=bool)
    for r, i in enumerate(ia):
        for k in ("A", "A_pstd", "sigma_r", "SE_sigma_r"):
            a = np.asarray(tracks[k][i], dtype=np.float64).reshape(n_channels, -1)
            arrays[k][r, : a.shape[1], :] = a.T
        g = np.asarray(tracks["gapVect"][i], dtype=bool).ravel()
        gaps[r, : g.size] = g
        for bk, src in (("startBuffer", ("sbA", "sbSigma_r")), ("endBuffer", ("ebA", "ebSigma_r"))):
            b = tracks[bk][i]
            if isinstance(b, dict) and "A" in b:
                bufs[src[0]][r] = np.asarray(b["A"], dtype=np.float64).reshape(n_channels, -1).T
                bufs[src[1]][r] = (
                    np.asarray(b["sigma_r"], dtype=np.float64).reshape(n_channels, -1).T
                )
    return {
        "lifetime_s": lifetime.reshape(-1, 1),
        "trackLengths": (end - start + 1).reshape(-1, 1),
        "start": start.reshape(-1, 1),
        "catIdx": cat.reshape(-1, 1),
        "index": np.arange(1, n + 1, dtype=np.float64).reshape(-1, 1),
        **arrays,
        **bufs,
        "gapMat_Ia": gaps,
        "significantMaster": flags("significantMaster"),
        "significantSlave": flags("significantSlave"),
    }


def _read_lifetime_file(path: str) -> dict:
    from petakit5d.io import matfile

    d = matfile.read_mat(path)
    out = {}
    for k in LIFETIME_FIELDS:
        v = d[k]
        a = np.asarray(v)
        if k in ("significantMaster", "significantSlave", "gapMat_Ia"):
            a = a.astype(bool)
        else:
            a = a.astype(np.float64)
        if k in ("lifetime_s", "trackLengths", "start", "catIdx", "index") and a.ndim < 2:
            a = a.reshape(-1, 1)
        out[k] = a
    return out


# -- visitors --------------------------------------------------------------------------------------


def visitor_index(lft: dict) -> np.ndarray:
    """MATLAB's ``getVisitorIndex(lftData)`` (measured): the tracks shorter than 30 s whose master
    channel amplitude in any of its first three frames exceeds the 95th percentile, frame by
    frame, of those frames over the tracks of 30 s and more; a logical column."""
    from petakit5d.matlabcompat import stats as ms

    L = np.asarray(lft["lifetime_s"], dtype=np.float64).ravel()
    A = np.asarray(lft["A"], dtype=np.float64)
    first = A[:, :3, 0]
    long_ = L >= VISITOR_LIFETIME
    thr = np.asarray(ms.prctile(first[long_], 95, 0, matlab_layout=True), dtype=np.float64).ravel()
    with np.errstate(invalid="ignore"):
        above = first > thr[None, :]
    return (L < VISITOR_LIFETIME) & np.any(above, axis=1)


# -- scaling between movies ------------------------------------------------------------------------


def _scale_cost(p, x, f, ref, grid):
    y = _interp1(x[1:], f[1:], grid[1:] / p[0])
    F = np.concatenate(([0.0], p[1] + (1 - p[1]) * y - ref[1:]))
    F[np.isnan(F)] = 0
    return F


def scale_edfs(samples: Sequence[Any]) -> tuple[np.ndarray, np.ndarray, int]:
    """MATLAB's ``[a, c, refIdx] = scaleEDFs(samples)`` (measured; 'Display' and the figure left
    out): the reference is the sample whose quantile function (on 0:0.001:1) is closest to the
    samples' median one; every other sample's empirical distribution, evaluated at ``x/a`` and
    mapped to ``c + (1 - c)*F``, is fitted to the reference's (1 above its largest value) on 1000
    points between the pooled 1st and 99th percentiles (``lsqnonlin`` from
    ``[median(ref)/median(sample), 0]``, bounds [0 -1] to [Inf 1], TolFun and TolX 1e-6). Returns
    a, c (one per sample; 1 and 0 for the reference) and the 1-based reference index."""
    from petakit5d.matlabcompat import optim
    from petakit5d.matlabcompat import reduce as mr
    from petakit5d.matlabcompat import stats as ms

    ss = [np.asarray(s, dtype=np.float64).ravel() for s in samples]
    ns = len(ss)
    if ns == 1:
        return np.array([1.0]), np.array([0.0]), 1
    edfs = [_ecdf(s) for s in ss]  # (f, x)
    pooled = np.concatenate(ss)
    lo = float(np.asarray(ms.prctile(pooled, 1)).ravel()[0])
    hi = float(np.asarray(ms.prctile(pooled, 99)).ravel()[0])
    pgrid = _colon(0.0, 0.001, 1.0)
    Q = np.vstack([_interp1(f, x, pgrid) for f, x in edfs])
    medQ = np.asarray(ms.median(Q, 0, matlab_layout=True), dtype=np.float64).reshape(1, -1)
    D = Q - medQ
    dev = np.asarray(mr.sum(D * D, 1, omitnan=True, matlab_layout=True), dtype=np.float64).ravel()
    ref = int(np.argmin(dev))
    grid = _linspace(lo, hi, 1000)
    fr_, xr_ = edfs[ref]
    refE = np.concatenate(([0.0], _interp1(xr_[1:], fr_[1:], grid[1:])))
    # outside the reference's values its EDF is 1 above (measured) and 0 below
    refE[(grid > xr_[-1]) & np.isnan(refE)] = 1.0
    refE[(grid < xr_[1]) & np.isnan(refE)] = 0.0
    a = np.ones(ns)
    c = np.zeros(ns)
    med_ref = float(np.asarray(ms.median(ss[ref])).ravel()[0])
    for i in range(ns):
        if i == ref:
            continue
        f, x = edfs[i]
        x0 = [med_ref / float(np.asarray(ms.median(ss[i])).ravel()[0]), 0.0]
        r = optim.lsqnonlin(
            _scale_cost,
            x0,
            [0.0, -1.0],
            [np.inf, 1.0],
            (x, f, refE, grid),
            tol_fun=1e-6,
            tol_x=1e-6,
            max_fun_evals=10000,
            max_iter=10000,
        )
        a[i], c[i] = r.x
    return a, c, ref + 1


# -- the Gaussian mode of a distribution -----------------------------------------------------------


def _mode_cost(p, x, f):
    from petakit5d.matlabcompat import special

    F = (p[0] * (0.5 * (1 + special.erf((x - p[1]) / _SQRT2 / p[2]))) - f) / x.size
    F[x > p[1]] = 0
    return F


def fit_gaussian_mode_to_cdf(samples) -> tuple[float, float]:
    """MATLAB's ``[mu, sigma] = fitGaussianModeToCDF(samples)`` (measured): ``A*Phi(x; mu, sigma)``
    fitted to the empirical distribution below its mode (``x <= mu``), from ``[0.5, mean, std/2]``
    with bounds [0 0 0] to [1 Inf Inf] (``lsqnonlin``, TolFun and TolX 1e-8). Constant samples
    (std 0) raise MATLAB's ``optimlib:commonMsgs:UsrObjUndefAtX0``."""
    from petakit5d.matlabcompat import optim
    from petakit5d.matlabcompat import reduce as mr
    from petakit5d.matlabcompat import stats as ms

    v = np.asarray(samples, dtype=np.float64).ravel()
    f, x = _ecdf(v)
    x0 = [
        0.5,
        float(np.asarray(mr.mean(v)).ravel()[0]),
        float(np.asarray(ms.std(v)).ravel()[0]) / 2,
    ]
    try:
        r = optim.lsqnonlin(
            _mode_cost,
            x0,
            [0.0, 0.0, 0.0],
            [1.0, np.inf, np.inf],
            (x, f),
            tol_fun=1e-8,
            tol_x=1e-8,
            max_fun_evals=10000,
            max_iter=10000,
        )
    except ValueError as e:
        ident = getattr(e, "identifier", "")
        raise LifetimeError(ident, str(e)) from e
    return float(r.x[1]), float(r.x[2])


# -- Anderson-Darling test -------------------------------------------------------------------------

_AD_CV = {"known": 2.492, "estimated": 0.752}  # 5 % critical values (Stephens 1974, cases 0 and 3)


def adtest1(
    x, mu: float | None = None, sigma: float | None = None
) -> tuple[bool, float, float, float]:
    """MATLAB's ``[H, p, adstat, cv] = adtest1(x, 'mu', mu, 'sigma', sigma)`` (measured, alpha
    0.05): the Anderson-Darling test of normality. With mu and sigma given, the statistic
    ``A^2 = -n - sum((2i - 1)(log F(x_i) + log(1 - F(x_{n+1-i}))))/n`` against 2.492 (p NaN);
    without them, mean and sample standard deviation, ``A^2 (1 + 0.75/n + 2.25/n^2)`` against 0.752
    and D'Agostino and Stephens' p-value. Fewer than 5 values are an error, as in MATLAB."""
    from petakit5d.matlabcompat import reduce as mr
    from petakit5d.matlabcompat import special
    from petakit5d.matlabcompat import stats as ms

    v = np.sort(np.asarray(x, dtype=np.float64).ravel(), kind="stable")
    n = v.size
    if n < 5:
        raise LifetimeError("", "adtest1 needs at least 5 values")
    known = mu is not None and sigma is not None
    if not known:
        mu = float(np.asarray(mr.mean(v)).ravel()[0])
        sigma = float(np.asarray(ms.std(v)).ravel()[0])
    F = 0.5 * special.erfc(-((v - mu) / sigma) / _SQRT2)
    i = np.arange(1, n + 1, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = (2 * i - 1) * (np.log(F) + np.log(1 - F[::-1]))
    s = _msum(terms)
    a2 = -n - s * (1 / n)  # (138 of 144 probed statistics bit for bit; MATLAB's log differs)
    if known:
        cv = _AD_CV["known"]
        p = float(dtypes.NAN)
        stat = a2
    else:
        cv = _AD_CV["estimated"]
        stat = a2 * (1 + 0.75 / n + 2.25 / (n * n))
        if stat < 0.2:
            p = 1 - math.exp(-13.436 + 101.14 * stat - 223.73 * stat * stat)
        elif stat < 0.34:
            p = 1 - math.exp(-8.318 + 42.796 * stat - 59.938 * stat * stat)
        elif stat < 0.6:
            p = math.exp(0.9177 - 4.279 * stat - 1.38 * stat * stat)
        else:
            p = math.exp(1.2937 - 5.709 * stat + 0.0186 * stat * stat)
    return bool(stat > cv), p, stat, cv


# -- B-spline interpolation (binterp) --------------------------------------------------------------

_Z = math.sqrt(3.0) - 2.0


def binterp(f, xi) -> np.ndarray:
    """MATLAB's ``binterp(f, xi)`` (CMEAnalysis' mex, no source) as measured: cubic B-spline
    interpolation of a 1-D signal with mirror-symmetric boundaries at 1-based positions xi in
    [0, n + 1] (outside 1..n mirrored), after Thevenaz, Blu and Unser, IEEE Trans. Med. Imaging
    19(7), 2000. Within ``docs/parity/detection.md`` L5 of MATLAB's. NaN samples and positions
    outside [0, n + 1] are errors, as in MATLAB (a single sample is not supported: MATLAB's
    result then depends on memory beyond its input)."""
    c = np.asarray(f, dtype=np.float64).ravel().copy()
    if np.isnan(c).any():
        raise LifetimeError("", "binterp: NaN samples")
    n = c.size
    xq = np.asarray(xi, dtype=np.float64).ravel()
    if not np.all((xq >= 0) & (xq <= n + 1)):
        raise LifetimeError(
            "",
            "binterp: the positions must lie in [0, n + 1] (the samples are at 1 to n)",
        )
    if n > 1:
        c *= 6.0  # the cubic B-spline's gain (a single sample is its own coefficient)
        z = _Z
        zn = z
        iz = 1.0 / z
        z2n = z ** (n - 1)
        s = c[0] + z2n * c[n - 1]
        z2n *= z2n * iz
        for k in range(1, n - 1):
            s += (zn + z2n) * c[k]
            zn *= z
            z2n *= iz
        c[0] = s / (1.0 - zn * zn)
        for k in range(1, n):
            c[k] += z * c[k - 1]
        c[n - 1] = (z / (z * z - 1.0)) * (z * c[n - 2] + c[n - 1])
        for k in range(n - 2, -1, -1):
            c[k] = z * (c[k + 1] - c[k])
    x = xq - 1.0
    fl = np.floor(x)
    t = x - fl
    w3 = (1.0 / 6.0) * t * t * t
    w0 = (1.0 / 6.0) + 0.5 * t * (t - 1.0) - w3
    w2 = t + w0 - 2.0 * w3
    w1 = 1.0 - w0 - w2 - w3
    i0 = fl.astype(np.int64) - 1
    p = 2 * (n - 1)

    def mirror(i):
        if n == 1:
            return np.zeros_like(i)
        i = np.abs(i) % p
        return np.where(i >= n, p - i, i)

    return (
        w0 * c[mirror(i0)]
        + w1 * c[mirror(i0 + 1)]
        + w2 * c[mirror(i0 + 2)]
        + w3 * c[mirror(i0 + 3)]
    )


# -- getLifetimeData -------------------------------------------------------------------------------


def _source(d: dict) -> str:
    """a movie's ``source`` with a trailing separator (the platform's, MATLAB's ``filesep``)"""
    from petakit5d.core import paths

    s = str(d["source"])
    return s if s.endswith(("/", paths.filesep)) else s + paths.filesep


def _movies(data) -> list[dict]:
    """the movies of ``data`` as dicts of MATLAB's fields: a struct (a dict), a struct array
    (``io.matfile.StructArray``), a condition (``detection.condition``'s movies) or a sequence of
    them"""
    if isinstance(data, Mapping):
        return [dict(data)]
    return [dict(d.to_matlab() if hasattr(d, "to_matlab") else d) for d in data]


def _n_channels(d: dict) -> int:
    ch = d.get("channels", [d["source"]])
    return len(ch) if isinstance(ch, list | tuple | np.ndarray) else 1


def _subset(lf: dict, rows: np.ndarray, ia_rows: np.ndarray, fields=None) -> dict:
    """the per-track fields at rows and the Ia arrays at ia_rows"""
    out = {}
    for k in fields or LIFETIME_FIELDS:
        a = lf[k]
        out[k] = a[ia_rows] if k in _IA_ARRAYS else a[rows]
    return out


def get_lifetime_data(
    data: Sequence[dict],
    *,
    overwrite: bool = False,
    return_valid_only: bool = False,
    exclude_visitors: bool = True,
    cutoff_f: int = CUTOFF_F,
    scale: bool = True,
    analysis_path: str = "Analysis",
    processed_tracks: str = "ProcessedTracks.mat",
    lifetime_data: str = "lifetimeData.mat",
) -> list[dict]:
    """MATLAB's ``getLifetimeData(data, ...)`` as the lifetime analysis calls it (module
    docstring): one dict per movie with MATLAB's fields and order; writes
    ``<source><analysis_path>/<lifetime_data>`` per movie unless it exists and ``overwrite`` is
    false (then it is read)."""
    from petakit5d.io import matfile

    movies = []
    raw_max = []
    for d in data:
        folder = os.path.join(_source(d), analysis_path)
        fn = os.path.join(folder, lifetime_data)
        nch = _n_channels(d)
        if overwrite or not os.path.exists(fn):
            tracks = load_tracks(os.path.join(folder, processed_tracks))
            lf = lifetime_data_from_tracks(tracks, int(d["movieLength"]), nch)
            matfile.write_mat(fn, lf, version="7")
        else:
            lf = _read_lifetime_file(fn)
        movies.append(lf)
    out = []
    for lf in movies:
        TL = lf["trackLengths"].ravel()
        cat = lf["catIdx"].ravel()
        keep_all = np.flatnonzero(TL >= cutoff_f)
        ia_all = np.flatnonzero(cat == 1)
        ia_keep = np.flatnonzero((cat == 1) & (TL >= cutoff_f))
        ia_rows = np.searchsorted(ia_all, ia_keep)
        raw_max.append([_nanmax_rows(lf["A"][:, :, ch]) for ch in range(lf["A"].shape[2])])
        e = _subset(lf, ia_keep, ia_rows)
        e["index"] = (np.searchsorted(keep_all, ia_keep) + 1).astype(np.float64).reshape(-1, 1)
        if not return_valid_only:
            for k in (
                "lifetime_s",
                "trackLengths",
                "start",
                "catIdx",
                "significantMaster",
                "significantSlave",
            ):
                e[k + "_all"] = lf[k][keep_all]
            e["index_all"] = np.arange(1, keep_all.size + 1, dtype=np.float64).reshape(-1, 1)
            # MATLAB's order: the _all fields after significantSlave, index_all after catIdx_all
            order = list(LIFETIME_FIELDS) + [
                "lifetime_s_all",
                "trackLengths_all",
                "start_all",
                "catIdx_all",
                "index_all",
                "significantMaster_all",
                "significantSlave_all",
            ]
            e = {k: e[k] for k in order}
        if exclude_visitors:
            vis = visitor_index(e)
            visitors = {k: e[k][vis] for k in (*_PER_TRACK, *_IA_ARRAYS)}
            visitors = {
                k: visitors[k]
                for k in (
                    "lifetime_s",
                    "trackLengths",
                    "start",
                    "catIdx",
                    "index",
                    "significantMaster",
                    "significantSlave",
                    "A",
                    "A_pstd",
                    "sigma_r",
                    "SE_sigma_r",
                    "sbA",
                    "ebA",
                    "sbSigma_r",
                    "ebSigma_r",
                    "gapMat_Ia",
                )
            }
            for k in (*_PER_TRACK, *_IA_ARRAYS):
                e[k] = e[k][~vis]
            e["visitors"] = [visitors]
        out.append(e)
    nch = max((lf["A"].shape[2] for lf in movies), default=1)
    if scale:
        factors = np.ones((len(out), nch))
        for ch in range(nch):
            a, _c, _ref = scale_edfs([rm[ch] for rm in raw_max])
            factors[:, ch] = a
        for i, e in enumerate(out):
            for k in _SCALED:  # the visitors stay unscaled (as MATLAB's)
                e[k] = e[k] * factors[i][None, None, :]
            e["a"] = factors[i].reshape(-1, 1)
    for e in out:
        e["maxA"] = (
            np.column_stack([_nanmax_rows(e["A"][:, :, ch]) for ch in range(e["A"].shape[2])])
            if e["A"].shape[0]
            else np.zeros((0, e["A"].shape[2]))
        )
    return out


# -- runLifetimeAnalysis3D -------------------------------------------------------------------------


def hist_weights(
    movie_length: int, cutoff_f: int = CUTOFF_F, buffer: int = BUFFER
) -> tuple[np.ndarray, np.ndarray]:
    """the histogram's frame counts ``L = cutoff_f:N`` and weights ``N/(N - L + 1)``, with
    ``N = movieLength - 2*buffer``"""
    N = int(movie_length) - 2 * buffer
    L = np.arange(cutoff_f, N + 1, dtype=np.float64)
    return L, N / (N - L + 1)


def _norm_hist(counts, w, framerate):
    cw = counts * w
    return cw / _msum(cw) / framerate


def _is_default(value: Any, default: Any) -> bool:
    """whether a value is the default: [] for None, text by text, numbers and logicals by value
    (element by element)"""
    if default is None:
        return value is None or np.size(value) == 0
    if isinstance(default, str):
        return isinstance(value, str) and value == default
    if isinstance(value, str):
        return False
    a = np.asarray(value, dtype=np.float64).reshape(-1)
    b = np.asarray(default, dtype=np.float64).reshape(-1)
    return a.shape == b.shape and bool(np.all(a == b))


_LB = (1, 11, 16, 21, 41, 61)
_UB = (10, 15, 20, 40, 60, 120)


@register(
    "runLifetimeAnalysis3D",
    # third-party code (its header names an outside author): its interface, with the port's own
    # checks (accepts=) made from MATLAB's behaviour (tests/params/test_third_party_checks.py:
    # an error of a check's own, not inputParser's, counts as taken, e.g. Display's cells of
    # other than one element); the options the port takes at their defaults only raise
    # NotImplementedError otherwise
    Param("data", None, "Movies (the structure XR_loadConditionData3D returns).", convert="same"),
    Param("lb", None, "Lower lifetime bounds of the threshold's cohorts, s.", convert="same"),
    Param("ub", None, "Upper lifetime bounds of the threshold's cohorts, s.", convert="same"),
    Param(
        "Display",
        None,
        "Figures, 'on' or 'off' (the port makes none).",
        accepts="in(off,on) | cell & !numel(1)",
    ),
    Param("DisplayMode", None, "Figures' mode (the default only).", accepts="in(print,screen)"),
    Param("AnalysisPath", None, "Analysis folder in each movie's source folder.", accepts="char"),
    Param("ProcessedTracks", None, "Track file in the analysis folder.", accepts="char"),
    Param("LifetimeData", None, "Lifetime data file in the analysis folder.", accepts="char"),
    Param("Type", None, "Not ported: the default only.", accepts="char"),
    Param("Cutoff_f", None, "Not ported: the default only (3 frames).", accepts="scalar"),
    Param("Print", None, "Figures to files (false only).", py="print_", accepts="logical"),
    Param("Buffer", None, "Not ported: the default only (5 frames)."),
    Param("MaxIntensityThreshold", None, "Not ported: [] only (the threshold is fitted)."),
    Param(
        "Overwrite",
        None,
        "Read the tracks again where lifetimeData.mat exists.",
        accepts="logical",
    ),
    Param("ClassificationSelector", None, "Not ported: the default only."),
    Param("ShowThresholdRange", None, "Not ported: false only.", accepts="logical"),
    Param("MaxP", None, "Not ported: the default only."),
    Param("YLim", None, "Figures' limits ([] only)."),
    Param("Rescale", None, "Not ported: true only (the movies are scaled).", accepts="logical"),
    Param("RemoveOutliers", None, "Not ported: true only.", accepts="logical"),
    Param(
        "ExcludeVisitors",
        None,
        "Not ported: true only (the visitors are left out).",
        accepts="logical",
    ),
    Param("FirstNFrames", None, "Not ported: [] only."),
    Param("PoolDatasets", None, "Not ported: true only.", convert="same", accepts="logical"),
    Param("ShowStatistics", None, "Not ported: true only.", accepts="logical"),
    Param("SelectIndex", None, "Not ported: [] only.", accepts="cell"),
    Param("SlaveNames", None, "Not ported: [] only."),
    Param(
        "Colormap",
        None,
        "Figures' colors (the default only).",
        matlab_default="jet(numel(data))",
    ),
    Param("CorrectObservationBias", None, "Not ported: true only.", accepts="logical"),
    Param("InitDensity", None, "Not ported: the default only.", accepts="in(mean,median)"),
    Param("AmplitudeCorrection", None, "Not ported: [] only."),
    Param(
        "MaskPath",
        None,
        "Cell mask in each movie's source folder (not ported: none may exist there).",
        matlab_default="['Detection' filesep 'cellmask.tif']",
        accepts="char",
    ),
    Param(
        "Binning",
        None,
        "Not ported: the default only.",
        accepts="numeric & scalar & positive-whole | cell | struct | string | handle",
    ),
    task="C05",
)
def run_lifetime_analysis(
    data,
    lb=_LB,
    ub=_UB,
    *,
    display: str = "on",
    display_mode: str = "screen",
    analysis_path: str = "Analysis",
    processed_tracks: str = "ProcessedTracks.mat",
    lifetime_data: str = "lifetimeData.mat",
    type: str = "all",
    cutoff_f: int = CUTOFF_F,
    print_: bool = False,
    buffer: int = BUFFER,
    max_intensity_threshold=None,
    overwrite: bool = False,
    classification_selector: str = "significantMaster",
    show_threshold_range: bool = False,
    max_p: int = 3,
    y_lim=None,
    rescale: bool = True,
    remove_outliers: bool = True,
    exclude_visitors: bool = True,
    first_n_frames=None,
    pool_datasets: bool = True,
    show_statistics: bool = True,
    select_index=None,
    slave_names=None,
    colormap=None,
    correct_observation_bias: bool = True,
    init_density: str = "mean",
    amplitude_correction=None,
    mask_path: str | None = None,
    binning: int = 1,
) -> dict:
    """MATLAB's ``lftRes = runLifetimeAnalysis3D(data, ...)`` (module docstring; the figures are
    left out). ``data``: the condition's movies (``Condition``, a struct array, dicts with MATLAB's
    fields: source, framerate, movieLength, channels). The options the port takes:
    ``analysis_path``, ``processed_tracks``, ``lifetime_data``, ``overwrite`` and ``display``
    (no figures either way); the others at their defaults only. Returns lftRes with MATLAB's
    fields and order; writes lifetimeData.mat per movie and Parameters.mat in the first movie's
    analysis folder."""
    from petakit5d.core import paths
    from petakit5d.io import matfile
    from petakit5d.matlabcompat import stats as ms

    for name, value, default in (
        ("lb", lb, _LB),
        ("ub", ub, _UB),
        ("DisplayMode", display_mode, "screen"),
        ("Type", type, "all"),
        ("Cutoff_f", cutoff_f, CUTOFF_F),
        ("Print", print_, False),
        ("Buffer", buffer, BUFFER),
        ("MaxIntensityThreshold", max_intensity_threshold, None),
        ("ClassificationSelector", classification_selector, "significantMaster"),
        ("ShowThresholdRange", show_threshold_range, False),
        ("MaxP", max_p, 3),
        ("YLim", y_lim, None),
        ("Rescale", rescale, True),
        ("RemoveOutliers", remove_outliers, True),
        ("ExcludeVisitors", exclude_visitors, True),
        ("FirstNFrames", first_n_frames, None),
        ("PoolDatasets", pool_datasets, True),
        ("ShowStatistics", show_statistics, True),
        ("SelectIndex", select_index, None),
        ("SlaveNames", slave_names, None),
        ("Colormap", colormap, None),
        ("CorrectObservationBias", correct_observation_bias, True),
        ("InitDensity", init_density, "mean"),
        ("AmplitudeCorrection", amplitude_correction, None),
        ("Binning", binning, 1),
    ):
        if not _is_default(value, default):
            raise NotImplementedError(
                f"run_lifetime_analysis: {name} {value!r}: only the default is ported"
            )
    data = _movies(data)
    if not data:
        raise LifetimeError("", "run_lifetime_analysis: no movies")
    mask = mask_path if mask_path is not None else "Detection" + paths.filesep + "cellmask.tif"
    for d in data:
        if os.path.isfile(_source(d) + mask):
            raise NotImplementedError(f"run_lifetime_analysis: cell masks are not ported ({mask})")
    if len({float(d["framerate"]) for d in data}) != 1:
        raise LifetimeError(
            "MATLAB:InputParser:ArgumentFailedValidation", "The movies' frame rates differ."
        )
    framerate = float(data[0]["framerate"])
    nf = int(data[0]["movieLength"])
    lft = get_lifetime_data(
        data,
        overwrite=overwrite,
        return_valid_only=False,
        exclude_visitors=True,
        cutoff_f=CUTOFF_F,
        scale=True,
        analysis_path=analysis_path,
        processed_tracks=processed_tracks,
        lifetime_data=lifetime_data,
    )
    Lf, w = hist_weights(nf)
    t = Lf * framerate
    cat_all = np.concatenate([e["catIdx_all"].ravel() for e in lft])
    L_all = np.concatenate([e["lifetime_s_all"].ravel() for e in lft])
    res: dict[str, Any] = {"t": t, "cellArea": dtypes.NAN}
    for name, ci in (("lftHist_Ia", 1), ("lftHist_Ib", 2), ("lftHist_IIa", 5)):
        res[name] = _norm_hist(hist(L_all[cat_all == ci], t), w, framerate)
    n_ia = int(np.sum(cat_all == 1))
    res["nSamples_Ia"] = float(n_ia)
    # kept Ia tracks, pooled over movies
    L = np.concatenate([e["lifetime_s"].ravel() for e in lft])
    A1 = np.concatenate([e["A"][:, :, 0] for e in lft], axis=0)
    sm = np.concatenate([e["significantMaster"] for e in lft], axis=0)
    vis_L = np.concatenate([e["visitors"][0]["lifetime_s"].ravel() for e in lft])
    # the threshold: cohorts of the first frames' maxima
    cohorts = [
        np.flatnonzero((L >= lo) & (L <= hi)) for lo, hi in zip(COHORT_LB, COHORT_UB, strict=True)
    ]
    nc = len(cohorts)
    hval = []
    for k in FRAME_RANGE:
        mus, sig = [], []
        for rows in cohorts:
            m_, s_ = fit_gaussian_mode_to_cdf(_nanmax_rows(A1[rows, :k]))
            mus.append(m_)
            sig.append(s_)
        h, *_ = adtest1(mus[1:], mu=mus[0], sigma=sig[0] / math.sqrt(nc))
        hval.append(h)
    first = next((i for i, h in enumerate(hval) if h), None)
    # (no frame count where they differ: the last count, not measured)
    n_frames = FRAME_RANGE[first] - 1 if first is not None else FRAME_RANGE[-1]
    mu, s = fit_gaussian_mode_to_cdf(_nanmax_rows(A1[:, :n_frames]))
    from petakit5d.detection import stats

    z99 = float(np.asarray(stats.norminv(0.99)).ravel()[0])
    thr = mu + s * z99
    p95 = float(np.asarray(ms.prctile(A1[:, 0], 95)).ravel()[0])
    # (MATLAB prints these two values, whatever 'Display' is; the port's words)
    print(f"Lifetime analysis: the CCP threshold is {thr:.2f} (from frames 1 to {n_frames})")
    print(f"Lifetime analysis: the tracks' first amplitudes have a 95th percentile of {p95:.2f}")
    res["MaxIntensityThreshold"] = thr
    maxA = _nanmax_rows(A1)
    ccp = maxA > thr
    cs = ~ccp
    n_tot = L.size + vis_L.size
    res["pctCCP"] = float(np.sum(ccp)) / n_tot
    res["pctCS"] = float(np.sum(cs)) / n_tot
    res["pctVisit"] = float(vis_L.size) / n_tot
    res["nCCP"] = float(np.sum(ccp))
    c_ccp = hist(L[ccp], t)
    c_cs = hist(L[cs], t)
    res["lftHistCCP"] = _norm_hist(c_ccp, w, framerate)
    res["lftHistCS"] = _norm_hist(c_cs, w, framerate)
    res["lftHistCCP_counts"] = c_ccp
    res["lftHistCS_counts"] = c_cs
    res["histWeight"] = w
    res["lftHistVisit"] = _norm_hist(hist(vis_L, t), w, framerate)
    nch = sm.shape[1]
    combos = _slave_combos(nch)
    res["slaveCombs"] = combos
    all_h, ccp_h, cs_h, all_c, ccp_c, cs_c, p_ccp, p_cs = [], [], [], [], [], [], [], []
    for row in combos:
        sel = np.all(sm[:, 1:] == row[None, :], axis=1)
        ca = hist(L[sel], t)
        cc = hist(L[sel & ccp], t)
        cx = hist(L[sel & cs], t)
        all_c.append(ca)
        ccp_c.append(cc)
        cs_c.append(cx)
        all_h.append(_norm_hist(ca, w, framerate))
        ccp_h.append(_norm_hist(cc, w, framerate))
        cs_h.append(cx * w / _msum(cx) / framerate)  # MATLAB normalizes these by the count (quirk)
        p_ccp.append(float(np.sum(sel & ccp)) / L.size)
        p_cs.append(float(np.sum(sel & cs)) / L.size)
    res["lftHistSlaveAll"] = all_h
    res["lftHistSlaveCCP"] = ccp_h
    res["lftHistSlaveCS"] = cs_h
    res["lftHistSlaveAll_counts"] = all_c
    res["lftHistSlaveCCP_counts"] = ccp_c
    res["lftHistSlaveCS_counts"] = cs_c
    res["pctSlaveCCP"] = np.array(p_ccp)
    res["pctSlaveCS"] = np.array(p_cs)
    starts = np.concatenate([e["start_all"].ravel() for e in lft])
    spf = hist(starts, np.arange(1, nf + 1, dtype=np.float64))
    res["startsPerFrameAll"] = spf[BUFFER : nf - 2]
    nan2 = np.full(2, dtypes.NAN)
    res["initDensityAll"] = nan2.copy()
    res["initDensityIa"] = nan2.copy()
    res["initDensityCCP"] = nan2.copy()
    res["meanLftHistCCP"] = res["lftHistCCP"].copy()
    res["meanLftHistCS"] = res["lftHistCS"].copy()
    res["meanLftHistVisit"] = res["lftHistVisit"].copy()
    params = {"MaxIntensityThreshold": thr, "nFramesThreshold": float(n_frames), "a": 1.0}
    matfile.write_mat(
        os.path.join(_source(data[0]), analysis_path, "Parameters.mat"), params, version="7"
    )
    return res


def _slave_combos(nch: int) -> np.ndarray:
    """the significance combinations of the other channels: rows of logicals, all true first
    (MATLAB's [1; 0] for two channels)"""
    k = nch - 1
    if k <= 0:
        return np.zeros((0, 0), dtype=bool)
    rows = [[(i >> (k - 1 - j)) & 1 == 0 for j in range(k)] for i in range(2**k)]
    return np.array(rows, dtype=bool).reshape(2**k, k)


def _quartiles(M: np.ndarray, nt: int):
    """the median, the 25th and the 75th percentiles of the rows of M per column (MATLAB's median
    and prctile); NaN rows without tracks"""
    from petakit5d.matlabcompat import stats as ms

    if not M.shape[0]:
        nan = np.full(nt, dtypes.NAN)
        return nan, nan.copy(), nan.copy()
    med = np.asarray(ms.median(M, 0), dtype=np.float64).reshape(-1)
    q25 = np.asarray(ms.prctile(M, 25, 0), dtype=np.float64).reshape(-1)
    q75 = np.asarray(ms.prctile(M, 75, 0), dtype=np.float64).reshape(-1)
    return med, q25, q75


# -- plotIntensityCohorts --------------------------------------------------------------------------


def plot_intensity_cohorts(
    data,
    *,
    analysis_path: str = "Analysis",
    lifetime_data: str = "lifetimeData.mat",
) -> list[dict]:
    """MATLAB's ``cohorts = plotIntensityCohorts(data, 'AnalysisPath', analysis_path)`` (module
    docstring; the figure left out): one dict per significance combination of the other channel
    (all true first), with ``bounds`` (in the first only), ``t``, ``Amin``, ``Aplus`` and ``A``
    (channels x cohorts object arrays of rows). Reads the lifetimeData.mat files (it writes them
    first where they are missing). With one movie, ``A`` is the median of the tracks and
    ``Amin``, ``Aplus`` their 25th and 75th percentiles (MATLAB's other branch, measured).
    ``data`` as ``run_lifetime_analysis`` takes it. Not registered: plotIntensityCohorts is
    third-party code outside the stored inputParser extract."""
    from petakit5d.matlabcompat import reduce as mr
    from petakit5d.matlabcompat import stats as ms

    data = _movies(data)
    framerate = float(data[0]["framerate"])
    lft = get_lifetime_data(
        data,
        overwrite=False,
        return_valid_only=True,
        exclude_visitors=True,
        cutoff_f=COHORT_CUTOFF_F,
        scale=True,
        analysis_path=analysis_path,
        lifetime_data=lifetime_data,
    )
    nch = lft[0]["A"].shape[2]
    b = lft[0]["sbA"].shape[1]
    bounds = np.array(INTENSITY_BOUNDS)
    ncoh = bounds.size - 1
    lens = [int(math.floor(((bounds[c] + bounds[c + 1]) / 2) / framerate)) for c in range(ncoh)]
    tvec = [np.arange(-b, L + b, dtype=np.float64) * framerate for L in lens]
    combos = _slave_combos(nch)
    if combos.shape[0] == 0:
        combos = np.ones((1, 0), dtype=bool)
    # interpolated tracks per movie, cohort and channel
    interp = []
    for e in lft:
        Lm = e["lifetime_s"].ravel()
        TL = e["trackLengths"].ravel().astype(np.int64)
        per = {}
        for c in range(ncoh):
            rows = np.flatnonzero((Lm >= bounds[c]) & (Lm < bounds[c + 1]))
            nt = lens[c] + 2 * b
            for ch in range(nch):
                M = np.empty((rows.size, nt))
                for r, i in enumerate(rows):
                    sig = np.concatenate(
                        [e["sbA"][i, :, ch], e["A"][i, : TL[i], ch], e["ebA"][i, :, ch]]
                    )
                    M[r] = binterp(sig, _linspace(1.0, float(sig.size), nt))
                per[(c, ch)] = (rows, M)
        interp.append(per)
    out = []
    for a_i, row in enumerate(combos):
        Am = np.empty((nch, ncoh), dtype=object)
        Amin = np.empty((nch, ncoh), dtype=object)
        Aplus = np.empty((nch, ncoh), dtype=object)
        for c in range(ncoh):
            nt = lens[c] + 2 * b
            for ch in range(nch):
                means = np.empty((len(lft), nt))
                for m, e in enumerate(lft):
                    rows, M = interp[m][(c, ch)]
                    sel = (
                        np.all(e["significantMaster"][rows][:, 1:] == row[None, :], axis=1)
                        if row.size
                        else np.ones(rows.size, dtype=bool)
                    )
                    sub = M[sel] if M.shape[0] else M
                    if len(lft) == 1:
                        # one movie: the median and the quartiles of its tracks (measured)
                        Am[ch, c], Amin[ch, c], Aplus[ch, c] = _quartiles(sub, nt)
                        break
                    means[m] = (
                        np.asarray(
                            mr.mean(sub, 0, omitnan=True, matlab_layout=True), dtype=np.float64
                        ).reshape(-1)
                        if sub.shape[0]
                        else np.full(nt, dtypes.NAN)
                    )
                if len(lft) == 1:
                    continue
                Am[ch, c] = np.asarray(
                    mr.mean(means, 0, omitnan=True, matlab_layout=True), dtype=np.float64
                ).reshape(-1)
                se = np.asarray(
                    ms.std(means, 0, 0, omitnan=True, matlab_layout=True), dtype=np.float64
                ).reshape(-1) / math.sqrt(len(lft))
                Amin[ch, c] = Am[ch, c] - se
                Aplus[ch, c] = Am[ch, c] + se
        out.append(
            {
                "bounds": bounds.copy() if a_i == 0 else np.zeros((0, 0)),
                "t": [tv.copy() for tv in tvec],
                "Amin": Amin,
                "Aplus": Aplus,
                "A": Am,
            }
        )
    return out
