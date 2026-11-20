"""Frame-to-frame linking with a Kalman filter for Brownian motion (u-track's linking step,
Jaqaman et al. 2008), as MATLAB's ``linkFeaturesKalmanSparse`` with ``costMatLinearMotionLink2``
and the ``kalman*LinearMotion`` functions computes it for ``linearMotion = 0``.

Written from the method and MATLAB's measured behaviour (the functions' inputs and outputs
recorded on the demo data and on probes; their sources were not read). Per frame pair:

- prediction: the state is not moved (Brownian motion); the position variance is the state's
  plus the process noise variance;
- search radius of a feature: ``s = sqrt(3 * q)`` with ``q`` its process noise variance, ``nn``
  the smallest nearest-neighbour distance of its track's last ``nn_window + 1`` frames, the
  multiple ``max(brown_std_mult, min(nn / s / 2, 100))``, the radius ``multiple * s`` clamped to
  ``[min_search_radius, max_search_radius]``; a link is allowed within it (inclusive);
- cost: the squared Euclidean distance (``pdist2``'s, squared); birth and death cost
  ``1.05`` times the largest cost of the links made so far (of this pass and the one before),
  or of the allowed costs when there is none; the lower right block the birth/death cost on the
  transposed links;
- update of a link ``i -> j``: ``k = p / (p + e^2 + eps)`` per axis (``p`` the predicted
  variance, ``e`` the observation's error), the state noise ``k (z - x)`` (kept on ``i``), the
  state the observation ``z``, the variance ``p - k p``, the process noise variance the variance
  (MATLAB's ``var``) of all the track's state noise of this pass, as a frames x 3 matrix in
  column order; a new feature takes the initial state (variance ``e^2`` and ``(max r / mult)^2``,
  noise ``(max r / mult)^2 / 3``) or, in the later passes, its state from the pass before.

Tracks of a pass are sorted by start frame, then end frame (descending), then the feature index
at the end, as MATLAB's.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["Frames", "KalmanInfo", "LinkResult", "link_pass", "reverse_kalman"]

EPS = float(np.finfo(np.float64).eps)
MAX_STD_MULT = 100.0
NAN = float("nan")


@dataclass
class Frames:
    """detections of a movie, concatenated: ``pos`` (N, 3) [x y z], ``err`` (N, 3), ``amp`` (N, 2)
    [A dA], ``nn`` (N,) nearest-neighbour distance in the frame, ``offsets`` (T + 1,)"""

    pos: np.ndarray
    err: np.ndarray
    amp: np.ndarray
    nn: np.ndarray
    offsets: np.ndarray

    @property
    def n_frames(self) -> int:
        return self.offsets.size - 1

    def count(self, t: int) -> int:
        return int(self.offsets[t + 1] - self.offsets[t])

    def sl(self, t: int) -> slice:
        return slice(int(self.offsets[t]), int(self.offsets[t + 1]))

    def reversed(self) -> Frames:
        counts = np.diff(self.offsets)[::-1]
        order = (
            np.concatenate(
                [
                    np.arange(self.offsets[t], self.offsets[t + 1])
                    for t in range(self.n_frames - 1, -1, -1)
                ]
            ).astype(np.int64)
            if self.n_frames
            else np.zeros(0, np.int64)
        )
        off = np.zeros(self.n_frames + 1, np.int64)
        np.cumsum(counts, out=off[1:])
        return Frames(self.pos[order], self.err[order], self.amp[order], self.nn[order], off)


@dataclass
class KalmanInfo:
    """the filter state of every detection, concatenated as ``Frames``: velocity (N, 3), position
    and velocity variances (N, 3) each, process noise variances of position and velocity (N,),
    state noise of the position (N, 3) and of the velocity (N, 3), scheme (N, 2)"""

    vel: np.ndarray
    pvar: np.ndarray
    vvar: np.ndarray
    qpos: np.ndarray
    qvel: np.ndarray
    sn: np.ndarray
    snv: np.ndarray
    scheme: np.ndarray

    @classmethod
    def empty(cls, n: int) -> KalmanInfo:
        z3 = np.zeros((n, 3))
        return cls(
            z3.copy(),
            z3.copy(),
            z3.copy(),
            np.zeros(n),
            np.zeros(n),
            z3.copy(),
            z3.copy(),
            np.zeros((n, 2)),
        )

    def take(self, order: np.ndarray) -> KalmanInfo:
        return KalmanInfo(
            *(
                getattr(self, f)[order]
                for f in ("vel", "pvar", "vvar", "qpos", "qvel", "sn", "snv", "scheme")
            )
        )


def reverse_kalman(info: KalmanInfo, frames: Frames) -> KalmanInfo:
    """``kalmanReverseLinearMotion``: frames in reverse order, velocities negated (0 becomes -0)"""
    rf = frames.reversed()
    order = (
        np.concatenate(
            [
                np.arange(frames.offsets[t], frames.offsets[t + 1])
                for t in range(frames.n_frames - 1, -1, -1)
            ]
        ).astype(np.int64)
        if frames.n_frames
        else np.zeros(0, np.int64)
    )
    out = info.take(order)
    out.vel = -out.vel
    del rf
    return out


@dataclass
class LinkResult:
    """one linking pass: ``tracks`` (n_tracks, T) 1-based feature indices (0: none), sorted as
    MATLAB's; ``nn_dist`` (n_tracks, T) NaN outside the tracks; ``kalman`` per detection;
    ``max_cost`` the largest link cost (of this pass and the ones before)"""

    tracks: np.ndarray
    nn_dist: np.ndarray
    kalman: KalmanInfo
    max_cost: float
    link_cost_all: list = field(default_factory=list)


def _mvar_rows(H: np.ndarray) -> np.ndarray:
    """MATLAB's ``var(v)`` of each row of ``H`` (rows are the vectors, all of one length)"""
    from petakit5d.matlabcompat.stats import var

    if H.shape[1] == 1:
        return np.zeros(H.shape[0])
    # MATLAB layout: columns are the vectors
    return np.asarray(var(np.ascontiguousarray(H.T), axis=0, matlab_layout=True)).reshape(-1)


def link_pass(
    frames: Frames,
    *,
    min_radius: float,
    max_radius: float,
    brown_std_mult: float,
    nn_window: int,
    prev: KalmanInfo | None = None,
    prev_max_cost: float = NAN,
) -> LinkResult:
    """one pass of frame-to-frame linking over the frames in their order"""
    from petakit5d.detection.tracking.lap import solve_assignment

    T = frames.n_frames
    N = int(frames.offsets[-1])
    off = frames.offsets
    kal = KalmanInfo.empty(N)
    v0 = (max_radius / brown_std_mult) ** 2  # initial velocity variance
    q0 = v0 / 3.0  # initial process noise variance
    # per detection: its track id, the track's nearest-neighbour distances of the last frames,
    # the track's history length (links so far)
    track = np.full(N, -1, np.int64)
    w = int(nn_window) + 1
    nnh = np.full((N, w), np.nan)
    hist_len = np.zeros(N, np.int64)
    # state noise history per track: (track, link, axis), filled as links are made
    hist = np.zeros((N, max(T - 1, 1), 3))
    n_tracks = 0
    starts: list[int] = []
    max_cost = prev_max_cost

    def init_new(idx: np.ndarray, t: int) -> None:
        nonlocal n_tracks
        if idx.size == 0:
            return
        g = idx + off[t]
        if prev is None:
            kal.vel[g] = 0.0
            kal.pvar[g] = frames.err[g] ** 2
            kal.vvar[g] = v0
            kal.qpos[g] = q0
            kal.qvel[g] = q0
        else:
            kal.vel[g] = prev.vel[g]
            kal.pvar[g] = prev.pvar[g]
            kal.vvar[g] = prev.vvar[g]
            kal.qpos[g] = prev.qpos[g]
            kal.qvel[g] = prev.qvel[g]
        ids = np.arange(n_tracks, n_tracks + idx.size)
        track[g] = ids
        starts.extend([t] * idx.size)
        n_tracks += idx.size
        nnh[g, :-1] = np.nan
        nnh[g, -1] = frames.nn[g]
        hist_len[g] = 0

    if T == 0:
        return LinkResult(np.zeros((0, 0)), np.zeros((0, 0)), kal, max_cost)
    init_new(np.arange(frames.count(0)), 0)
    for t in range(T - 1):
        s1, s2 = frames.sl(t), frames.sl(t + 1)
        n1, n2 = s1.stop - s1.start, s2.stop - s2.start
        if n2 == 0:
            continue
        if n1 == 0:
            init_new(np.arange(n2), t + 1)
            continue
        x1 = frames.pos[s1]
        x2 = frames.pos[s2]
        q = kal.qpos[s1]
        pp = kal.pvar[s1] + q[:, None]  # predicted position variance
        pv = kal.vvar[s1] + kal.qvel[s1][:, None]
        kstd = np.sqrt(3.0 * q)
        nn = np.fmin.reduce(nnh[s1], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = nn / kstd / 2.0
        ratio = np.where(np.isnan(ratio), MAX_STD_MULT, np.minimum(ratio, MAX_STD_MULT))
        mult = np.maximum(brown_std_mult, ratio)
        r = np.minimum(np.maximum(mult * kstd, min_radius), max_radius)
        r = np.where(np.isnan(r), min_radius, r)
        d = np.sqrt(
            ((x1[:, None, 0] - x2[None, :, 0]) ** 2 + (x1[:, None, 1] - x2[None, :, 1]) ** 2)
            + (x1[:, None, 2] - x2[None, :, 2]) ** 2
        )
        ok = d <= r[:, None]
        ii, jj = np.nonzero(ok)
        cost = d[ii, jj] ** 2
        if np.isnan(max_cost):
            bd = 1.05 * (max(float(cost.max()), 0.0) if cost.size else 0.0)
        else:
            bd = 1.05 * max_cost
        n = n1 + n2
        rows = np.concatenate([ii, np.arange(n1), n1 + np.arange(n2), n1 + jj])
        cols = np.concatenate([jj, n2 + np.arange(n1), np.arange(n2), n2 + ii])
        vals = np.concatenate([cost, np.full(n1 + n2 + ii.size, bd)])
        x = solve_assignment(rows, cols, vals, n)
        lk = np.nonzero(x[:n1] < n2)[0]
        lj = x[lk]
        if lk.size:
            c = d[lk, lj] ** 2
            cm = float(c.max())
            max_cost = cm if np.isnan(max_cost) else max(max_cost, cm)
        g1 = lk + off[t]
        g2 = lj + off[t + 1]
        # update of the linked features
        p = pp[lk]
        e2 = frames.err[g2] ** 2
        k = p / (p + (e2 + EPS))
        sn = k * (frames.pos[g2] - x1[lk])
        kal.sn[g1] = sn
        # the velocity rows of the gain are 0: their noise is -0 when every product is -0
        kal.snv[g1] = np.where(np.signbit(sn).all(axis=1), -0.0, 0.0)[:, None]
        kal.scheme[g1, 1] = 1.0
        kal.vel[g2] = 0.0
        kal.pvar[g2] = p - k * p
        kal.vvar[g2] = pv[lk]
        kal.scheme[g2, 0] = 1.0
        tid = track[g1]
        track[g2] = tid
        L = hist_len[g1] + 1
        hist[tid, L - 1] = sn
        hist_len[g2] = L
        nnh[g2, :-1] = nnh[g1, 1:]
        nnh[g2, -1] = frames.nn[g2]
        qn = np.empty(lk.size)
        for Lk in np.unique(L):
            sel = np.nonzero(L == Lk)[0]
            H = hist[tid[sel], :Lk].transpose(0, 2, 1).reshape(sel.size, 3 * Lk)
            qn[sel] = _mvar_rows(H)
        kal.qpos[g2] = qn
        kal.qvel[g2] = 0.0
        new = np.setdiff1d(np.arange(n2), lj, assume_unique=True)
        init_new(new, t + 1)
    # tracks matrix, sorted by start, end (descending), feature index at the end
    tracks = np.zeros((n_tracks, T), np.int64)
    nnd = np.full((n_tracks, T), np.nan)
    for t in range(T):
        s = frames.sl(t)
        tids = track[s]
        tracks[tids, t] = np.arange(1, s.stop - s.start + 1)
        nnd[tids, t] = frames.nn[s]
    start = np.asarray(starts, np.int64)
    present = tracks > 0
    end = T - 1 - np.argmax(present[:, ::-1], axis=1)
    last = tracks[np.arange(n_tracks), end]
    order = np.lexsort((last, -end, start))
    return LinkResult(tracks[order], nnd[order], kal, max_cost)
