"""The tracker as MATLAB's ``trackCloseGapsKalmanSparse`` runs it with the demos' settings: three
linking passes (forward, backward, forward; each later pass starts its new features from the pass
before, reversed), gap closing with merging and splitting on the last pass's segments, and the
compound tracks.

A compound track (one element of ``tracksFinal``) is a set of segments joined by gap closing,
merging and splitting; compound tracks are ordered by their first segment (segments are ordered
by start), their rows (segments joined by closed gaps) by first segment. Per compound track:

- ``tracksFeatIndxCG``: (rows, frames) feature indices (0: none) from its first frame to its
  last, double;
- ``tracksCoordAmpCG``: (rows, 8 * frames) ``[x y z A dx dy dz dA]`` per frame, NaN where none;
- ``seqOfEvents``: rows ``[frame, 1 (start) or 2 (end), row, partner]``: a start's partner is
  the row it splits from, an end's the row it merges into (the end's frame is then the frame of
  the merge, one after its last), NaN for births and deaths; the starts (by row), then the ends
  (by row), stably sorted by frame, a merge's end counted at its last frame.
"""

from __future__ import annotations

import numpy as np

from petakit5d.core.dtypes import NAN
from petakit5d.detection.tracking.gapclosing import gap_cost_problem
from petakit5d.detection.tracking.lap import solve_assignment
from petakit5d.detection.tracking.linking import Frames, KalmanInfo, link_pass, reverse_kalman

__all__ = ["frames_from_frame_info", "kalman_to_matlab", "track_frames"]

# the nearest-neighbour distance MATLAB gives the only detection of a frame (measured)
NN_ALONE = 1000.0


def _col(x, n: int, k: int) -> np.ndarray:
    a = np.asarray(x, dtype=np.float64)
    if a.size == 0:
        return np.zeros((n, k))
    return a.reshape(n, k) if a.ndim == 2 else a.reshape(-1, k)


def frames_from_frame_info(frame_info) -> Frames:
    """the detections of ``Detection3D.mat``'s ``frameInfo`` (a sequence of per-frame dicts with
    ``xCoord``, ``yCoord``, ``zCoord``, ``amp``: (n, 2) [value, error]) as ``Frames``, with each
    detection's distance to its nearest neighbour in the frame (MATLAB's ``knnsearch``; 1000
    for the only detection of a frame, as MATLAB gives it)"""
    from petakit5d.matlabcompat.neighbors import knnsearch

    pos, err, amp, nn, counts = [], [], [], [], []
    for f in frame_info:
        xc = np.asarray(f["xCoord"], dtype=np.float64)
        m = 0 if xc.size == 0 else xc.shape[0]
        xc = _col(xc, m, 2)
        yc = _col(f["yCoord"], m, 2)
        zc = _col(f["zCoord"], m, 2)
        am = _col(f["amp"], m, 2)
        P = np.c_[xc[:, 0], yc[:, 0], zc[:, 0]]
        pos.append(P)
        err.append(np.c_[xc[:, 1], yc[:, 1], zc[:, 1]])
        amp.append(am)
        if m >= 2:
            _, d = knnsearch(P, P, 2)
            nn.append(d[:, 1])
        else:
            nn.append(np.full(m, NN_ALONE))
        counts.append(m)
    off = np.zeros(len(counts) + 1, np.int64)
    np.cumsum(counts, out=off[1:])
    cat = lambda xs, k: np.concatenate(xs).reshape(-1, k) if xs else np.zeros((0, k))  # noqa: E731
    return Frames(
        cat(pos, 3), cat(err, 3), cat(amp, 2), np.concatenate(nn) if nn else np.zeros(0), off
    )


def kalman_to_matlab(kal: KalmanInfo, frames: Frames) -> list[dict]:
    """the filter state per frame as MATLAB's ``kalmanInfoLink``: ``stateVec`` (n, 6),
    ``stateCov`` and ``noiseVar`` (6, 6, n) diagonal, ``stateNoise`` (n, 6), ``scheme`` (n, 2)"""
    out = []
    for t in range(frames.n_frames):
        s = frames.sl(t)
        m = s.stop - s.start
        if m == 0:
            out.append(
                {
                    "stateVec": np.zeros((0, 6)),
                    "stateCov": np.zeros((6, 6, 0)),
                    "noiseVar": np.zeros((6, 6, 0)),
                    "stateNoise": np.zeros((0, 6)),
                    "scheme": np.zeros((0, 2)),
                }
            )
            continue
        sv = np.c_[frames.pos[s], kal.vel[s]]
        dcov = np.c_[kal.pvar[s], kal.vvar[s]]
        dnv = np.c_[np.repeat(kal.qpos[s][:, None], 3, 1), np.repeat(kal.qvel[s][:, None], 3, 1)]
        cov = np.zeros((6, 6, m))
        nv = np.zeros((6, 6, m))
        ar = np.arange(6)
        cov[ar, ar, :] = dcov.T
        nv[ar, ar, :] = dnv.T
        if m == 1:
            cov = cov[:, :, 0]
            nv = nv[:, :, 0]
        out.append(
            {
                "stateVec": sv,
                "stateCov": cov,
                "noiseVar": nv,
                "stateNoise": np.c_[kal.sn[s], kal.snv[s]],
                "scheme": kal.scheme[s].copy(),
            }
        )
    return out


def _compound(tracks: np.ndarray, frames: Frames, x: np.ndarray, gp) -> list[dict]:
    """the compound tracks from the segments and the gap-closing assignment ``x``"""
    n, T = tracks.shape
    nm = gp.merge_seg.size
    ns = gp.split_seg.size
    present = tracks > 0
    ts = np.argmax(present, axis=1)
    te = T - 1 - np.argmax(present[:, ::-1], axis=1)
    nxt = np.full(n, -1)  # gap link: segment continues as
    prv = np.full(n, -1)
    merge_into = np.full(n, -1)  # segment ends by merging into segment
    split_from = np.full(n, -1)  # segment starts by splitting from segment
    for r in range(n):
        c = x[r]
        if c < n:
            nxt[r] = c
            prv[c] = r
        elif c < n + nm:
            merge_into[r] = gp.merge_seg[c - n]
    for k in range(ns):
        c = x[n + k]
        if c < n:
            split_from[c] = gp.split_seg[k]
    # union of segments into compound tracks
    parent = np.arange(n)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    for a in range(n):
        if nxt[a] >= 0:
            union(a, nxt[a])
        if merge_into[a] >= 0:
            union(a, merge_into[a])
        if split_from[a] >= 0:
            union(a, split_from[a])
    roots = np.array([find(a) for a in range(n)])
    gidx = np.where(present, tracks - 1 + frames.offsets[:-1][None, :], -1)
    out = []
    order = np.unique(roots)  # the smallest segment of each compound track
    members_of = {}
    for a in range(n):
        members_of.setdefault(roots[a], []).append(a)
    for root in order:
        segs = members_of[root]
        heads = [s for s in segs if prv[s] < 0]
        heads.sort()
        row_of = {}
        chains = []
        for h in heads:
            chain = [h]
            while nxt[chain[-1]] >= 0:
                chain.append(nxt[chain[-1]])
            for s in chain:
                row_of[s] = len(chains)
            chains.append(chain)
        t0 = min(ts[s] for s in segs)
        t1 = max(te[s] for s in segs)
        nf = t1 - t0 + 1
        F = np.zeros((len(chains), nf))
        C = np.full((len(chains), 8 * nf), NAN)
        ev_s, ev_e = [], []
        for k, chain in enumerate(chains):
            for s in chain:
                for t in range(ts[s], te[s] + 1):
                    F[k, t - t0] = tracks[s, t]
                    g = gidx[s, t]
                    C[k, 8 * (t - t0) : 8 * (t - t0) + 8] = np.r_[
                        frames.pos[g], frames.amp[g, 0], frames.err[g], frames.amp[g, 1]
                    ]
            h, e = chain[0], chain[-1]
            sp_ = split_from[h]
            ev_s.append((ts[h] + 1, ts[h] + 1, 1, k + 1, NAN if sp_ < 0 else row_of[sp_] + 1))
            mi = merge_into[e]
            if mi < 0:
                ev_e.append((te[e] + 1, te[e] + 1, 2, k + 1, NAN))
            else:
                ev_e.append((te[e] + 1, te[e] + 2, 2, k + 1, row_of[mi] + 1))
        # starts, then ends, each by row, stably sorted by frame (a merge's end counted at the
        # segment's last frame, reported one frame later)
        ev = np.array(ev_s + ev_e, dtype=np.float64)
        ev = ev[np.argsort(ev[:, 0], kind="stable"), 1:]
        out.append({"tracksFeatIndxCG": F, "tracksCoordAmpCG": C, "seqOfEvents": ev})
    return out


def track_frames(frames: Frames, settings: dict) -> tuple[list[dict], KalmanInfo]:
    """``trackCloseGapsKalmanSparse(movieInfo, costMatrices, gapCloseParam, kalmanFunctions, 3,
    ...)``: the compound tracks (MATLAB's ``tracksFinal`` elements) and the last linking pass's
    filter state"""
    lp = settings["costMatrices"][0]["parameters"]
    gp = settings["costMatrices"][1]["parameters"]
    gc = settings["gapCloseParam"]
    if float(lp["linearMotion"]) != 0 or float(gp["linearMotion"]) != 0:
        raise NotImplementedError("tracking: only Brownian motion (linearMotion 0) is ported")
    kw = {
        "min_radius": float(lp["minSearchRadius"]),
        "max_radius": float(lp["maxSearchRadius"]),
        "brown_std_mult": float(lp["brownStdMult"]),
        "nn_window": int(lp["nnWindow"]),
    }
    r1 = link_pass(frames, **kw)
    rf = frames.reversed()
    r2 = link_pass(rf, prev=reverse_kalman(r1.kalman, frames), prev_max_cost=r1.max_cost, **kw)
    r3 = link_pass(frames, prev=reverse_kalman(r2.kalman, rf), prev_max_cost=r2.max_cost, **kw)
    tracks = r3.tracks
    T = frames.n_frames
    tw = int(gc["timeWindow"])
    prob = gap_cost_problem(
        tracks,
        frames,
        r3.kalman,
        r3.nn_dist,
        min_radius=float(gp["minSearchRadius"]),
        max_radius=float(gp["maxSearchRadius"]),
        brown_std_mult=np.asarray(gp["brownStdMult"], float).ravel(),
        brown_scaling=tuple(np.asarray(gp["brownScaling"], float).ravel()),
        time_reach_conf_b=float(gp["timeReachConfB"]),
        len_for_classify=int(gp["lenForClassify"]),
        nn_window=int(gp["nnWindow"]),
        time_window=tw,
        merge_split=bool(gc["mergeSplit"]),
    )
    n = tracks.shape[0]
    nm, ns = prob.merge_seg.size, prob.split_seg.size
    nr, nc = n + ns, n + nm
    if prob.vals.size and T > 1:
        bd = prob.alt_cost
        rows = np.concatenate([prob.rows, np.arange(nr), nr + np.arange(nc), nr + prob.cols])
        cols = np.concatenate([prob.cols, nc + np.arange(nr), np.arange(nc), nc + prob.rows])
        vals = np.concatenate([prob.vals, np.full(nr + nc + prob.vals.size, bd)])
        x = solve_assignment(rows, cols, vals, nr + nc)
    else:
        x = np.full(nr + nc, n + nm + 10**9)
    return _compound(tracks, frames, x, prob), r3.kalman
