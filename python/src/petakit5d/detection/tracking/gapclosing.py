"""Gap closing, merging and splitting of track segments (u-track's second step, Jaqaman et al.
2008), as MATLAB's ``costMatLinearMotionCloseGaps2`` with ``estimTrackTypeParamLM2`` and
``getAveDispEllipseAll2`` computes it for Brownian motion (``linearMotion = 0``).

Written from the method and MATLAB's measured behaviour (the functions' inputs and outputs on
the demo data and on synthetic track sets; their sources were not read):

- per segment: the mean step length ``md`` (segments of 2 or more frames), its noise std
  ``sqrt(q)`` at its end for segments of ``len_for_classify`` (5) or more frames, else 1;
  a 1 (also an estimate of exactly 1) is replaced by ``md / sqrt(2)``; segments of 1 frame use
  the 10th percentile of the estimated stds (none: NaN, and they close no gap);
- search radius at a segment's end (start) for a gap of ``g`` frames: with ``s`` the std times
  ``g ** 0.5`` and ``nn`` the smallest nearest-neighbour distance of the segment's last (first)
  ``nn_window + 1`` frames, ``max(brown_std_mult, min(nn / s / 2, 100)) * s * sqrt(3)``, at least
  ``min_search_radius`` and at most ``max_search_radius * g ** 0.5``;
- a gap from the end of ``i`` to the start of ``j`` is allowed when the distance is within
  ``j``'s start radius (``i`` shorter than 5 frames), or within ``i``'s end radius and, for ``j`` of
  5 frames or more, ``j``'s start radius; its cost ``d^2 / (m g^0.5)^2`` with ``m`` the mean of
  the two segments' ``md`` (a 1-frame segment's is the mean of all);
- ``i`` ending at ``t - 1`` merges into ``j`` (present at ``t - 1`` and ``t``) within ``i``'s end
  radius for one frame, if the amplitude of ``j`` after (mean of its frames ``t .. t + 4``) is
  larger than that of ``i`` (mean of its last 5 frames) and ``a / b + a / (b + o) > 2`` (``b``:
  ``j``'s mean of frames ``t - 5 .. t - 1``, ``o``: ``i``'s), at the cost
  ``d^2 rho / md_i^2`` (``rho = a / (b + o)``, ``rho^-2`` below 1); splits are the mirror image;
- the merge columns and split rows by time, then segment; birth and death cost 1.05 times the
  largest cost.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from petakit5d.detection.tracking.linking import MAX_STD_MULT, Frames, KalmanInfo

__all__ = ["GapProblem", "gap_cost_problem"]

AMP_WINDOW = 5


@dataclass
class GapProblem:
    """the gap-closing assignment problem: allowed entries of the upper left block (rows: segment
    ends then split rows; columns: segment starts then merge columns), the merge columns and split
    rows as (segment, time), the birth/death cost"""

    rows: np.ndarray
    cols: np.ndarray
    vals: np.ndarray
    n: int
    merge_seg: np.ndarray
    merge_time: np.ndarray
    split_seg: np.ndarray
    split_time: np.ndarray
    alt_cost: float


def _mean_rows(W: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """MATLAB's ``mean`` of the first ``counts[k]`` values of row ``k`` (sequential sums: at most
    ``AMP_WINDOW`` values)"""
    s = W[:, 0].copy()
    for c in range(1, W.shape[1]):
        m = counts > c
        s[m] = s[m] + W[m, c]
    return s / counts


def _window_means(A: np.ndarray, starts: np.ndarray, stops: np.ndarray) -> np.ndarray:
    """mean of ``A[k, starts[k]:stops[k]]`` (at most ``AMP_WINDOW`` values, MATLAB's order)"""
    cnt = stops - starts
    idx = starts[:, None] + np.arange(AMP_WINDOW)[None, :]
    idx = np.minimum(idx, A.shape[1] - 1)
    W = A[np.arange(A.shape[0])[:, None], idx]
    return _mean_rows(W, cnt)


def _amplitude_rule(aa: np.ndarray, bb: np.ndarray, oo: np.ndarray):
    """merges and splits: allowed where the joint amplitude ``aa`` is larger than the other
    segment's ``oo`` and ``aa / bb + aa / (bb + oo) > 2`` (``bb`` the joint segment's amplitude
    on the other side of the event); the amplitude factor of the allowed ones, ``rho`` or
    ``rho^-2`` below 1, ``rho = aa / (bb + oo)``"""
    keep = (aa > oo) & (aa / bb + aa / (bb + oo) > 2)
    rho = aa[keep] / (bb[keep] + oo[keep])
    return keep, np.where(rho >= 1, rho, rho**-2.0)


def _vector_mean(v: np.ndarray, omitnan: bool = False) -> float:
    from petakit5d.matlabcompat import reduce as mred

    return float(np.ravel(mred.mean(np.asarray(v, np.float64), omitnan=omitnan))[0])


def _mean_step(P: np.ndarray, ts: np.ndarray, te: np.ndarray) -> np.ndarray:
    """mean step length of each segment (NaN for one frame), MATLAB's ``mean`` of
    ``sqrt(sum(diff(xyz).^2, 2))``"""
    from petakit5d.matlabcompat import reduce as mred

    n = ts.size
    md = np.full(n, np.nan)
    L = te - ts + 1
    for Lk in np.unique(L):
        if Lk < 2:
            continue
        sel = np.nonzero(L == Lk)[0]
        idx = ts[sel, None] + np.arange(Lk)[None, :]
        X = P[sel[:, None], idx]  # (m, Lk, 3)
        dX = X[:, 1:] - X[:, :-1]
        st = np.sqrt((dX[..., 0] ** 2 + dX[..., 1] ** 2) + dX[..., 2] ** 2)  # (m, Lk - 1)
        md[sel] = np.asarray(
            mred.mean(np.ascontiguousarray(st.T), axis=0, matlab_layout=True)
        ).reshape(-1)
    return md


def _radii(sig, nn, W, bsm, min_r, max_r, t_reach, scaling):
    """(segments, W) search radii for gaps of 1..W frames"""
    out = np.empty((sig.size, W))
    for g in range(1, W + 1):
        tg = (
            g ** scaling[0] if g <= t_reach else t_reach ** scaling[0] * (g / t_reach) ** scaling[1]
        )
        st = sig * tg
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = nn / st / 2.0
        ratio = np.where(np.isnan(ratio), MAX_STD_MULT, np.minimum(ratio, MAX_STD_MULT))
        mult = np.maximum(bsm[g - 1], ratio)
        r = mult * st * np.sqrt(3.0)
        r = np.minimum(np.maximum(r, min_r), max_r * tg)
        # a std of 0 (a segment that does not move) or NaN gives no search region (MATLAB)
        r = np.where(np.isnan(sig) | (sig == 0), np.nan, r)
        out[:, g - 1] = r
    return out


def gap_cost_problem(
    tracks: np.ndarray,
    frames: Frames,
    kalman: KalmanInfo,
    nn_dist: np.ndarray,
    *,
    min_radius: float,
    max_radius: float,
    brown_std_mult,
    brown_scaling=(0.5, 0.01),
    time_reach_conf_b: float = 3,
    len_for_classify: int = 5,
    nn_window: int = 3,
    time_window: int = 3,
    merge_split: bool = True,
) -> GapProblem:
    """the gap-closing problem of the segments ``tracks`` (n, T) (1-based feature indices)"""
    from petakit5d.matlabcompat.stats import prctile

    n, T = tracks.shape
    present = tracks > 0
    ts = np.argmax(present, axis=1)
    te = T - 1 - np.argmax(present[:, ::-1], axis=1)
    L = te - ts + 1
    gidx = np.where(present, tracks - 1 + frames.offsets[:-1][None, :], -1)
    P = np.where(present[..., None], frames.pos[np.maximum(gidx, 0)], np.nan)  # (n, T, 3)
    A = np.where(present, frames.amp[np.maximum(gidx, 0), 0], np.nan)
    md = _mean_step(P, ts, te)
    md_all = _vector_mean(md, omitnan=True)
    det = L >= len_for_classify
    rows_all = np.arange(n)
    ns = np.ones(n)
    ge = gidx[rows_all, te]
    ns[det] = np.sqrt(kalman.qpos[ge[det]])
    est = ns.copy()
    repl = (ns == 1) & np.isfinite(md)
    ns[repl] = md[repl] / np.sqrt(2.0)
    good = det & (est != 1)
    undet = float(np.ravel(prctile(est[good], 10))[0]) if good.any() else np.nan
    sig = ns.copy()
    sig[~det & (ns == 1)] = undet
    w = int(nn_window) + 1
    nnE = np.empty(n)
    nnS = np.empty(n)
    for k in range(n):
        nnE[k] = np.fmin.reduce(nn_dist[k, max(te[k] - w + 1, ts[k]) : te[k] + 1])
        nnS[k] = np.fmin.reduce(nn_dist[k, ts[k] : min(ts[k] + w, te[k] + 1)])
    bsm = np.broadcast_to(np.asarray(brown_std_mult, float).ravel(), (time_window,))
    rE = _radii(
        sig, nnE, time_window, bsm, min_radius, max_radius, time_reach_conf_b, brown_scaling
    )
    rS = _radii(
        sig, nnS, time_window, bsm, min_radius, max_radius, time_reach_conf_b, brown_scaling
    )
    mdx = np.where(np.isnan(md), md_all, md)
    PE = P[rows_all, te]
    PS = P[rows_all, ts]
    ri, ci, vi = [], [], []
    # gaps
    for t in range(T - 1):
        ends = np.nonzero(te == t)[0]
        if ends.size == 0:
            continue
        for g in range(1, time_window + 1):
            if t + g >= T:
                break
            starts = np.nonzero(ts == t + g)[0]
            if starts.size == 0:
                continue
            a = PE[ends]
            b = PS[starts]
            d = np.sqrt(
                ((a[:, None, 0] - b[None, :, 0]) ** 2 + (a[:, None, 1] - b[None, :, 1]) ** 2)
                + (a[:, None, 2] - b[None, :, 2]) ** 2
            )
            okE = d <= rE[ends, g - 1][:, None]
            okS = d <= rS[starts, g - 1][None, :]
            dI = det[ends][:, None]
            dJ = det[starts][None, :]
            ok = np.where(dI, okE & np.where(dJ, okS, True), okS)
            p, q = np.nonzero(ok)
            if p.size == 0:
                continue
            m = (mdx[ends[p]] + mdx[starts[q]]) / 2.0
            c = d[p, q] ** 2 / (m * g**0.5) ** 2
            ri.append(ends[p])
            ci.append(starts[q])
            vi.append(c)
    merge_seg = np.zeros(0, np.int64)
    merge_time = np.zeros(0, np.int64)
    split_seg = np.zeros(0, np.int64)
    split_time = np.zeros(0, np.int64)
    if merge_split:
        mr, mc, mv, mt = [], [], [], []
        sr, sc, sv, st_ = [], [], [], []
        for t in range(1, T):
            mids = np.nonzero(present[:, t - 1] & present[:, t])[0]
            if mids.size == 0:
                continue
            # merges: segments ending at t - 1 into mids at t
            ends = np.nonzero(te == t - 1)[0]
            if ends.size:
                a = PE[ends]
                b = P[mids, t]
                d = np.sqrt(
                    ((a[:, None, 0] - b[None, :, 0]) ** 2 + (a[:, None, 1] - b[None, :, 1]) ** 2)
                    + (a[:, None, 2] - b[None, :, 2]) ** 2
                )
                p, q = np.nonzero(d <= rE[ends, 0][:, None])
                if p.size:
                    i_, j_ = ends[p], mids[q]
                    aa = _window_means(
                        A[j_], np.full(j_.size, t), np.minimum(t + AMP_WINDOW, te[j_] + 1)
                    )
                    bb = _window_means(
                        A[j_], np.maximum(t - AMP_WINDOW, ts[j_]), np.full(j_.size, t)
                    )
                    oo = _window_means(
                        A[i_], np.maximum(te[i_] - AMP_WINDOW + 1, ts[i_]), te[i_] + 1
                    )
                    keep, ac = _amplitude_rule(aa, bb, oo)
                    if keep.any():
                        dd = d[p[keep], q[keep]]
                        mr.append(i_[keep])
                        mc.append(j_[keep])
                        mt.append(np.full(int(keep.sum()), t))
                        mv.append(dd**2 * ac / mdx[i_[keep]] ** 2)
            # splits: segments starting at t from mids at t - 1
            starts = np.nonzero(ts == t)[0]
            if starts.size:
                a = P[mids, t - 1]
                b = PS[starts]
                d = np.sqrt(
                    ((a[:, None, 0] - b[None, :, 0]) ** 2 + (a[:, None, 1] - b[None, :, 1]) ** 2)
                    + (a[:, None, 2] - b[None, :, 2]) ** 2
                )
                q, p = np.nonzero(d <= rS[starts, 0][None, :])
                if p.size:
                    i_, j_ = mids[q], starts[p]
                    aa = _window_means(
                        A[i_], np.maximum(t - AMP_WINDOW, ts[i_]), np.full(i_.size, t)
                    )
                    bb = _window_means(
                        A[i_], np.full(i_.size, t), np.minimum(t + AMP_WINDOW, te[i_] + 1)
                    )
                    oo = _window_means(A[j_], ts[j_], np.minimum(ts[j_] + AMP_WINDOW, te[j_] + 1))
                    keep, ac = _amplitude_rule(aa, bb, oo)
                    if keep.any():
                        dd = d[q[keep], p[keep]]
                        sr.append(i_[keep])
                        sc.append(j_[keep])
                        st_.append(np.full(int(keep.sum()), t))
                        sv.append(dd**2 * ac / mdx[j_[keep]] ** 2)
        if mr:
            mr_, mc_, mt_, mv_ = (np.concatenate(x) for x in (mr, mc, mt, mv))
            key = np.unique(np.c_[mt_, mc_], axis=0)  # (time, segment) sorted
            merge_time, merge_seg = key[:, 0], key[:, 1]
            col = np.searchsorted(merge_time * (n + 1) + merge_seg, mt_ * (n + 1) + mc_)
            ri.append(mr_)
            ci.append(n + col)
            vi.append(mv_)
        if sr:
            sr_, sc_, st2, sv_ = (np.concatenate(x) for x in (sr, sc, st_, sv))
            key = np.unique(np.c_[st2, sr_], axis=0)
            split_time, split_seg = key[:, 0], key[:, 1]
            row = np.searchsorted(split_time * (n + 1) + split_seg, st2 * (n + 1) + sr_)
            ri.append(n + row)
            ci.append(sc_)
            vi.append(sv_)
    rows = np.concatenate(ri) if ri else np.zeros(0, np.int64)
    cols = np.concatenate(ci) if ci else np.zeros(0, np.int64)
    vals = np.concatenate(vi) if vi else np.zeros(0)
    alt = 1.05 * float(np.ravel(prctile(vals, 100))[0]) if vals.size else np.nan
    return GapProblem(rows, cols, vals, n, merge_seg, merge_time, split_seg, split_time, alt)
