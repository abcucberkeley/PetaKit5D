"""Track processing (``petakit5d.detection.processing``, ``detection.amira``; task C04b) against
MATLAB's ``XR_runTrackProcessing3D`` and ``rotateTracks3D`` on designed conditions, through the
files they write (``ProcessedTracks.mat``, ``RotatedTracks.mat``, the Amira files).

A designed condition is a folder as ``XR_loadConditionData3D`` reads it, with the tracker's
output and the detections written directly into ``ch1/Analysis`` (``trackedFeatures.mat``,
``Detection3D.mat``, ``Masks``), so that each track tests one rule of the processing: gaps of
one and two frames, with and without a spot in them (a "ghost": rendered, not detected), with
the p-values the gap test reads set (``pval_set``); spots in the buffer frames; an amplitude
that rises and falls (multiple peaks); tracks cut at the movie's ends; two-frame, single-frame
and border tracks; compound tracks with a split and a merge; a track that moves in large
steps. Frames are Gaussian spots on a background of 100 with noise of standard deviation 5,
bit-stable (``petakit5d.testing.synth``).
"""

from __future__ import annotations

import functools
import glob
import os

import numpy as np
import pytest

from petakit5d.core.dtypes import mround
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value
from petakit5d.testing.compare import ignore, tol
from petakit5d.testing.compare import tree as tree_spec

EXPERIMENT = "Ex1_probe_z0p5"
FRAME_MS = 2000
NF = 40
SIZE = (24, 90, 110)  # (z, y, x)
SIGMA = ((1.3, 1.27), (1.5, 1.51))  # per channel (sigma_xy, sigma_z) in voxels
Z_ANISO = 2.5120123303651387  # sin(31.5 deg) 0.5 / 0.104, as MATLAB's data.zAniso of z0p5
FI_FIELDS = (
    "x", "y", "z", "A", "s", "c", "x_pstd", "y_pstd", "z_pstd", "A_pstd", "c_pstd", "x_init",
    "y_init", "z_init", "sigma_r", "SE_sigma_r", "RSS", "pval_Ar", "hval_Ar", "hval_AD", "isPSF",
    "xCoord", "yCoord", "zCoord", "amp", "dRange",
)  # fmt: skip
CLEAR = (
    "clear XR_runTrackProcessing3D rotateTracks3D XR_amiraWriteTracks XR_loadConditionData3D "
    "XR_loadConditionData XR_fitGaussian3D fitGaussian3D estGaussianAmplitude3D; "
)


# ---------------------------------------------------------------------------------------------
# the designed condition
# ---------------------------------------------------------------------------------------------


def frame_name(channel: int, frame: int) -> str:
    cam, wl = ("B", 560) if channel == 1 else ("A", 642)
    t = (frame - 1) * FRAME_MS
    return (
        f"probe_Cam{cam}_ch{channel - 1}_stack{frame - 1:04d}_{wl}nm_{t:07d}msec_"
        f"{1000000 + t:010d}msecAbs.tif"
    )


def _jitter(name: str, n: int, s: float) -> np.ndarray:
    """(n, 3) bit-stable position offsets in [-0.3 s, 0.3 s) (x, y), [-0.2 s, 0.2 s) (z)"""
    u = synth.CounterRNG(2026, f"processing/jitter/{name}").uniform((n, 3)) - 0.5
    return u * np.array([0.6, 0.6, 0.4]) * s


def _track(cell, f0: int, amps, *, gaps=(), ghost=1.0, jit=0.0, name="", shift=None, pval_set=()):
    """a single-segment track at ``cell`` (x, y) from frame ``f0``, an amplitude per frame,
    gap points (0-based) holding a ghost at ``ghost`` times the amplitude (0: none)"""
    x, y = cell
    n = len(amps)
    J = _jitter(name or f"{x}/{y}/{f0}", n, jit) if jit else np.zeros((n, 3))
    pts, ghosts = [], []
    for k in range(n):
        dx = 0.0 if shift is None or k < shift[0] else shift[1]
        p = [x + J[k, 0] + dx, y + J[k, 1], 12.3 + J[k, 2], float(amps[k])]
        if k in gaps:
            pts.append(None)
            if ghost > 0:
                ghosts.append([f0 + k, p[0], p[1], p[2], p[3] * ghost])
        else:
            pts.append(p)
    return {
        "segments": [{"start": f0, "pts": pts}],
        "pval_set": [list(v) for v in pval_set],
    }, ghosts


def _spec() -> dict:
    cells = [(x, y) for y in range(12, 82, 12) for x in range(12, 102, 12)]
    tracks, ghosts = [], []

    def add(tr_g):
        tracks.append(tr_g[0])
        ghosts.extend(tr_g[1])

    # the pool: plain tracks, amplitudes 150 to 425 (the lifetime cohorts), small steps
    for i in range(12):
        add(_track(cells[i], 6 + i % 5, [150.0 + 25.0 * i] * (8 + i), jit=1.0, name=f"pool{i}"))
    ramp = [280.0 + 3.0 * k for k in range(14)]
    pv = ((3, 0, 0.5), (3, 1, 0.5))  # the gap columns 7-8 reach column 4 of both channels
    c = iter(cells[12:])
    add(_track(next(c), 6, ramp, gaps=(6, 7), pval_set=pv, jit=1.0, name="gap2_pass"))
    add(_track(next(c), 6, ramp, gaps=(6, 7), jit=1.0, name="gap2_fail"))
    add(_track(next(c), 6, ramp, gaps=(6,), pval_set=((3, 0, 0.5),), jit=1.0, name="gap1_pass"))
    add(_track(next(c), 8, [280.0, 283.0, 286.0], gaps=(1,), pval_set=((0, 1, 0.5),), name="n3"))
    add(_track(next(c), 8, [28.0, 28.3, 28.6], gaps=(1,), pval_set=((0, 1, 0.5),), name="n3low"))
    add(
        _track(
            next(c), 8, [280.0] * 4, gaps=(1, 2), ghost=0.0, pval_set=((0, 1, 0.5),), name="dense4"
        )
    )
    # buffers: a spot in the nearest start frame, a brighter one in the farthest start frame,
    # a dimmer one in the farthest end frame
    tr = _track(next(c), 8, [300.0] * 12, name="buf_near")
    add(tr)
    ghosts.append([7, *tr[0]["segments"][0]["pts"][0][:3], 300.0])
    tr = _track(next(c), 8, [300.0] * 12, name="buf_far_bright")
    add(tr)
    ghosts.append([5, *tr[0]["segments"][0]["pts"][0][:3], 360.0])
    tr = _track(next(c), 8, [300.0] * 12, name="buf_far_dim")
    add(tr)
    ghosts.append([22, *tr[0]["segments"][0]["pts"][-1][:3], 240.0])
    # multiple peaks (40 s > 30 s), cut tracks, a two-frame and a single-frame track, a border
    add(_track(next(c), 5, [200.0 + 150.0 * ((k // 4) % 2 == 0) for k in range(20)], name="peaks"))
    add(_track(next(c), 1, [250.0] * 10, name="cut_start"))
    add(_track(next(c), 31, [250.0] * 10, name="cut_end"))
    add(_track(next(c), 1, [220.0] * NF, name="whole"))
    add(_track(next(c), 12, [260.0, 262.0], name="two"))
    add(_track(next(c), 14, [260.0], name="one"))
    tr = _track((5.4, 40.0), 10, [300.0] * 9, name="border")
    add(tr)
    # large steps: five or more above the 95th percentile of the median steps
    add(_track(next(c), 6, ramp, jit=6.0, name="steps"))
    # compound: a branch splits off at frame 12 and dies at 20; a branch born at 8 merges at 18
    x, y = next(c)
    a = [[x, y, 12.3, 300.0] for _ in range(20)]
    b = [[x + 3.2, y + 0.4, 12.3, 260.0] for _ in range(9)]
    tracks.append(
        {
            "segments": [{"start": 6, "pts": a}, {"start": 12, "pts": b}],
            "events": [[6, 1, 1, np.nan], [12, 1, 2, 1], [20, 2, 2, np.nan], [25, 2, 1, np.nan]],
        }
    )
    x, y = next(c)
    a = [[x, y, 12.3, 300.0] for _ in range(20)]
    b = [[x - 3.1, y + 0.5, 12.3, 250.0] for _ in range(10)]
    tracks.append(
        {
            "segments": [{"start": 6, "pts": a}, {"start": 8, "pts": b}],
            "events": [[6, 1, 1, np.nan], [8, 1, 2, np.nan], [18, 2, 2, 1], [25, 2, 1, np.nan]],
        }
    )
    return {"tracks": tracks, "ghosts": ghosts}


def _detections(spec: dict) -> list[list[list[float]]]:
    det = [[] for _ in range(NF)]
    for tr in spec["tracks"]:
        tr["_idx"] = []
        for seg in tr["segments"]:
            idx = []
            for j, p in enumerate(seg["pts"]):
                if p is None:
                    idx.append(0)
                    continue
                f = seg["start"] - 1 + j
                det[f].append(list(p))
                idx.append(len(det[f]))
            tr["_idx"].append(idx)
    return det


def _render(spec: dict, det, channel: int, f: int) -> np.ndarray:
    s = SIGMA[channel - 1]
    D = np.asarray(det[f], dtype=np.float64).reshape(-1, 4)
    G = [g[1:] for g in spec["ghosts"] if int(g[0]) == f + 1]
    if G:
        D = np.vstack([D, np.asarray(G, dtype=np.float64).reshape(-1, 4)])
    sh = (0.0, 0.0, 0.0) if channel == 1 else (0.0, 0.0, 0.0)
    pos = np.stack([D[:, 2] - 1 + sh[2], D[:, 1] - 1 + sh[1], D[:, 0] - 1 + sh[0]], axis=1)
    amp = D[:, 3] * (1.0 if channel == 1 else 0.7)
    vol = synth.beads(SIZE, pos, (s[1], s[0]), amp, 100.0)
    noise = synth.CounterRNG(2026, f"processing/f{f}/ch{channel}").normal(vol.shape)
    v = vol.astype(np.float64) + 5.0 * noise
    t = np.trunc(v)
    v = t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)
    return np.clip(v, 1, 65535).astype(np.uint16)


def _mask(det, f: int) -> np.ndarray:
    nz, ny, nx = SIZE
    s = SIGMA[0]
    m = np.zeros(SIZE, dtype=bool)
    zz, yy, xx = np.meshgrid(
        np.arange(nz) + 1.0, np.arange(ny) + 1.0, np.arange(nx) + 1.0, indexing="ij"
    )
    for x, y, z, _a in det[f]:
        m |= ((xx - x) ** 2 + (yy - y) ** 2) / (2 * s[0]) ** 2 + (zz - z) ** 2 / (
            2 * s[1]
        ) ** 2 <= 1
    return m.astype(np.uint8) * 255


def _frame_info(spec: dict, det, vols):
    from petakit5d.io.matfile import StructArray

    rows = []
    for f, D in enumerate(det):
        n = len(D)
        D = np.asarray(D, dtype=np.float64).reshape(n, 4)
        e = {
            "x": np.vstack([D[:, 0], D[:, 0]]),
            "y": np.vstack([D[:, 1], D[:, 1]]),
            "z": np.vstack([D[:, 2], D[:, 2]]),
            "A": np.vstack([D[:, 3], D[:, 3] * 0.7]),
            "s": np.array([list(SIGMA[0]), list(SIGMA[1])], dtype=np.float64),
            "c": np.full((2, n), 100.0),
            "x_pstd": np.full((2, n), 0.05),
            "y_pstd": np.full((2, n), 0.05),
            "z_pstd": np.full((2, n), 0.08),
            "A_pstd": np.full((2, n), 4.0),
            "c_pstd": np.full((2, n), 0.5),
            "x_init": mround(D[:, 0]).reshape(1, n),
            "y_init": mround(D[:, 1]).reshape(1, n),
            "z_init": mround(D[:, 2]).reshape(1, n),
            "sigma_r": np.full((2, n), 5.0),
            "SE_sigma_r": np.full((2, n), 0.25),
            "RSS": np.full((2, n), 8550.0),
            "pval_Ar": np.zeros((2, n)),
            "hval_Ar": np.ones((2, n), dtype=bool),
            "hval_AD": np.zeros((2, n), dtype=bool),
            "isPSF": np.ones((2, n), dtype=bool),
        }
        e["xCoord"] = np.stack([D[:, 0], e["x_pstd"][0]], axis=1)
        e["yCoord"] = np.stack([D[:, 1], e["y_pstd"][0]], axis=1)
        e["zCoord"] = np.stack([D[:, 2], e["z_pstd"][0]], axis=1) * Z_ANISO
        e["amp"] = np.stack([D[:, 3], e["A_pstd"][0]], axis=1)
        e["dRange"] = [
            np.array([[float(v.min()), float(v.max())]]) for v in (vols[0][f], vols[1][f])
        ]
        if n == 0:
            for k in FI_FIELDS:
                if k not in ("s", "dRange"):
                    e[k] = np.zeros((0, 0))
        rows.append({k: e[k] for k in FI_FIELDS})
    for tr in spec["tracks"]:
        seg = tr["segments"][0]
        for j, ch, v in tr.get("pval_set", []):
            f = seg["start"] - 1 + int(j)
            rows[f]["pval_Ar"][int(ch), tr["_idx"][0][int(j)] - 1] = float(v)
    return StructArray((1, len(rows)), FI_FIELDS, rows)


def _tracks_final(spec: dict, det):
    from petakit5d.io.matfile import StructArray

    out = []
    for tr in spec["tracks"]:
        segs = tr["segments"]
        f0 = min(s["start"] for s in segs)
        f1 = max(s["start"] + len(s["pts"]) - 1 for s in segs)
        n = f1 - f0 + 1
        fi = np.zeros((len(segs), n))
        ca = np.full((len(segs), 8 * n), np.nan)
        for k, (seg, idx) in enumerate(zip(segs, tr["_idx"], strict=True)):
            for j, d in enumerate(idx):
                f = seg["start"] + j
                fi[k, f - f0] = d
                if d:
                    x, y, z, a = det[f - 1][d - 1]
                    ca[k, 8 * (f - f0) : 8 * (f - f0) + 8] = [
                        x,
                        y,
                        z * Z_ANISO,
                        a,
                        0.05,
                        0.05,
                        0.08 * Z_ANISO,
                        4.0,
                    ]
        if "events" in tr:
            ev = np.asarray(tr["events"], dtype=np.float64).reshape(-1, 4)
        else:
            s, e = segs[0]["start"], segs[0]["start"] + len(segs[0]["pts"]) - 1
            ev = np.array([[s, 1, 1, np.nan], [e, 2, 1, np.nan]], dtype=np.float64)
        out.append({"tracksFeatIndxCG": fi, "tracksCoordAmpCG": ca, "seqOfEvents": ev})
    return StructArray((len(out), 1), ("tracksFeatIndxCG", "tracksCoordAmpCG", "seqOfEvents"), out)


def _compound_spec() -> dict:
    """compound tracks as the tracker writes them, one structure each (a main segment of 8
    frames from frame 8, a branch beside it): branches that split off and die, are born and
    merge, end or start with the main, continue a main that ends at the split, merge into a
    segment born 1 to 3 frames before (the joins), loops, two branches, branches with gaps, and
    structures at the movie's first and last frames; and the pool of plain tracks"""
    cells = [(x, y) for y in range(12, 82, 12) for x in range(12, 102, 12)]
    tracks = []
    for i in range(12):
        tracks.append(
            _track(cells[i], 6 + i % 5, [150.0 + 25.0 * i] * (8 + i), jit=1.0, name=f"cpool{i}")[0]
        )
    nan = np.nan
    c = iter(cells[12:])

    def comp(segs, ev):
        """segs: (start, frames, dx, dy, amplitude, gap points); ev: [frame, 1|2, segment,
        partner or NaN]; the tracker's event order: starts, then ends, stably by frame, a
        merge's end at the merging segment's last frame"""
        x, y = next(c)
        out = []
        for f0, n, dx, dy, a, gaps in segs:
            pts = [None if k in gaps else [x + dx, y + dy, 12.3, a] for k in range(n)]
            out.append({"start": f0, "pts": pts})
        starts = sorted([e for e in ev if e[1] == 1], key=lambda e: e[2])
        ends = sorted([e for e in ev if e[1] == 2], key=lambda e: e[2])

        def key(e):
            return e[0] - (1 if (e[1] == 2 and not np.isnan(e[3])) else 0)

        tracks.append({"segments": out, "events": sorted(starts + ends, key=key)})

    M = (8, 8, 0.0, 0.0, 300.0, ())  # the main segment, frames 8-15
    for L in (1, 2, 3, 4):  # splits off at 10 and dies after L frames
        comp(
            [M, (10, L, 3.0, 2.0, 150.0, ())],
            [[8, 1, 1, nan], [10, 1, 2, 1], [9 + L, 2, 2, nan], [15, 2, 1, nan]],
        )
    for L in (1, 2, 3, 4):  # born at 10, merges into the main after L frames
        comp(
            [M, (10, L, 3.0, 2.0, 150.0, ())],
            [[8, 1, 1, nan], [10, 1, 2, nan], [10 + L, 2, 2, 1], [15, 2, 1, nan]],
        )
    for L in (2, 3):  # splits off and ends with the main
        comp(
            [M, (16 - L, L, 3.0, 2.0, 150.0, ())],
            [[8, 1, 1, nan], [16 - L, 1, 2, 1], [15, 2, 2, nan], [15, 2, 1, nan]],
        )
    for L in (2, 3):  # born with the main, merges after L frames
        comp(
            [M, (8, L, 3.0, 2.0, 150.0, ())],
            [[8, 1, 1, nan], [8, 1, 2, nan], [8 + L, 2, 2, 1], [15, 2, 1, nan]],
        )
    for L in (1, 2, 3):  # the main ends at the split, the branch continues L frames
        comp(
            [(8, 5, 0.0, 0.0, 300.0, ()), (12, L, 3.0, 2.0, 300.0, ())],
            [[8, 1, 1, nan], [12, 1, 2, 1], [12, 2, 1, nan], [11 + L, 2, 2, nan]],
        )
    # L frames merge at 12 into a segment that starts at 12 (L = 1 is left out: the branch goes,
    # the single track's span is longer than its columns, and MATLAB's XR_amiraWriteTracks stops
    # with MATLAB:badsubscript in a parfor, after some frames, W18)
    for L in (2, 3):
        comp(
            [(12 - L, L, 3.0, 2.0, 300.0, ()), (12, 4, 0.0, 0.0, 300.0, ())],
            [[12 - L, 1, 1, nan], [12, 1, 2, nan], [12, 2, 1, 2], [15, 2, 2, nan]],
        )
    for L in (1, 2, 3):  # 4 frames merge at 12 into a segment born L frames before (joins)
        comp(
            [(8, 4, 3.0, 2.0, 300.0, ()), (12 - L, 4 + L, 0.0, 0.0, 300.0, ())],
            [[8, 1, 1, nan], [12 - L, 1, 2, nan], [12, 2, 1, 2], [15, 2, 2, nan]],
        )
    for L in (1, 2, 3, 4):  # loops: split at 10, merge back after L frames
        comp(
            [M, (10, L, 3.0, 2.0, 150.0, ())],
            [[8, 1, 1, nan], [10, 1, 2, 1], [10 + L, 2, 2, 1], [15, 2, 1, nan]],
        )
    comp(
        [M, (9, 2, 3.0, 2.0, 150.0, ()), (12, 1, -3.0, 2.0, 150.0, ())],
        [
            [8, 1, 1, nan],
            [9, 1, 2, 1],
            [10, 2, 2, nan],
            [12, 1, 3, 1],
            [12, 2, 3, nan],
            [15, 2, 1, nan],
        ],
    )
    comp(
        [M, (9, 2, 3.0, 2.0, 150.0, ()), (12, 2, -3.0, 2.0, 150.0, ())],
        [
            [8, 1, 1, nan],
            [9, 1, 2, 1],
            [10, 2, 2, nan],
            [12, 1, 3, nan],
            [14, 2, 3, 1],
            [15, 2, 1, nan],
        ],
    )
    comp(
        [M, (10, 3, 3.0, 2.0, 150.0, (1,))],
        [[8, 1, 1, nan], [10, 1, 2, 1], [12, 2, 2, nan], [15, 2, 1, nan]],
    )
    comp(
        [(8, 8, 0.0, 0.0, 300.0, (3,)), (10, 2, 3.0, 2.0, 150.0, ())],
        [[8, 1, 1, nan], [10, 1, 2, 1], [11, 2, 2, nan], [15, 2, 1, nan]],
    )
    comp(
        [(1, 8, 0.0, 0.0, 300.0, ()), (1, 2, 3.0, 2.0, 150.0, ())],
        [[1, 1, 1, nan], [1, 1, 2, nan], [3, 2, 2, 1], [8, 2, 1, nan]],
    )
    comp(
        [(NF - 7, 8, 0.0, 0.0, 300.0, ()), (NF - 1, 2, 3.0, 2.0, 150.0, ())],
        [[NF - 7, 1, 1, nan], [NF - 1, 1, 2, 1], [NF, 2, 2, nan], [NF, 2, 1, nan]],
    )
    # a long main with a merge, a split and gaps (vis 1, both kinds of gaps)
    comp(
        [
            (6, 24, 0.0, 0.0, 300.0, (5, 12, 13)),
            (10, 6, 3.0, 2.0, 200.0, ()),
            (20, 5, -3.0, 2.0, 200.0, ()),
        ],
        [
            [6, 1, 1, nan],
            [10, 1, 2, nan],
            [16, 2, 2, 1],
            [20, 1, 3, 1],
            [24, 2, 3, nan],
            [29, 2, 1, nan],
        ],
    )
    # a join where the choice at the overlap frame (11) needs the interpolation: the merging
    # segment's detection there is the one near the middle of its frame 10 and the receiver's
    # frame 12, the receiver's own is near the merging segment's frame 10
    x, y = next(c)
    a = [[x, y, 12.3, 300.0]] * 3 + [[x + 2.2, y, 12.3, 300.0]]
    b = [[x + 0.5, y + 0.3, 12.3, 280.0]] + [[x + 4.0, y, 12.3, 280.0]] * 5
    tracks.append(
        {
            "segments": [{"start": 8, "pts": a}, {"start": 11, "pts": b}],
            "events": [[8, 1, 1, nan], [11, 1, 2, nan], [12, 2, 1, 2], [16, 2, 2, nan]],
        }
    )
    # the same join, the merging segment with a branch that split off it at 9: after the join
    # the branch's partner is the joined segment
    x, y = next(c)
    a = [[x, y, 12.3, 300.0]] * 3 + [[x + 2.2, y, 12.3, 300.0]]
    br = [[x - 3.0, y + 2.0, 12.3, 200.0]] * 6
    b = [[x + 0.5, y + 0.3, 12.3, 280.0]] + [[x + 4.0, y, 12.3, 280.0]] * 5
    tracks.append(
        {
            "segments": [{"start": 8, "pts": a}, {"start": 9, "pts": br}, {"start": 11, "pts": b}],
            "events": [
                [8, 1, 1, nan],
                [9, 1, 2, 1],
                [11, 1, 3, nan],
                [12, 2, 1, 3],
                [14, 2, 2, nan],
                [16, 2, 3, nan],
            ],
        }
    )
    # a compound track over the whole movie (visibility 3): a branch splits off at 20 for 5 frames
    comp(
        [(1, NF, 0.0, 0.0, 300.0, ()), (20, 5, 3.0, 2.0, 200.0, ())],
        [[1, 1, 1, nan], [20, 1, 2, 1], [24, 2, 2, nan], [NF, 2, 1, nan]],
    )
    return {"tracks": tracks, "ghosts": []}


SPECS = {"designed": _spec, "compound": _compound_spec}


@functools.lru_cache(maxsize=2)
def _designed(kind: str = "designed"):
    """the spec, its detections per frame and the rendered frames (both channels)"""
    spec = SPECS[kind]()
    det = _detections(spec)
    vols = {c: [_render(spec, det, c, f) for f in range(NF)] for c in (1, 2)}
    return spec, det, vols


def make_condition(folder, kind: str = "designed") -> None:
    """a designed condition's frames (raw and DS) and the detection masks in ch1/Analysis; the
    detections and the tracker's tracks are the case's values, saved by each side"""
    from pathlib import Path

    from petakit5d.detection.run import dmask_name

    spec, det, vols = _designed(kind)
    ex = Path(folder) / EXPERIMENT
    ana = ex / "ch1" / "Analysis"
    (ana / "Masks").mkdir(parents=True, exist_ok=True)
    for c in (1, 2):
        (ex / f"ch{c}" / "DS").mkdir(parents=True, exist_ok=True)
        for f in range(NF):
            synth.write_tiff(ex / f"ch{c}" / frame_name(c, f + 1), vols[c][f])
            synth.write_tiff(ex / f"ch{c}" / "DS" / frame_name(c, f + 1), vols[c][f])
    with open(ex / "probe_Settings.txt", "w", newline="\n") as fh:
        fh.write('Camera type = "Orca4.0"\n')
    for f in range(NF):
        synth.write_tiff(ana / "Masks" / dmask_name(f + 1, NF), _mask(det, f))


def frame_info(kind: str = "designed"):
    spec, det, vols = _designed(kind)
    return _frame_info(spec, det, [vols[1], vols[2]])


def tracks_final(kind: str = "designed"):
    spec, det, _ = _designed(kind)
    return _tracks_final(spec, det)


def _mstruct(sa):
    """a StructArray as MATLAB's struct array (the case's value)"""
    from petakit5d.testing.compare import MStruct, to_mvalue

    return MStruct(
        sa.shape, list(sa.fields), [{k: to_mvalue(el[k]) for k in sa.fields} for el in sa.ravel()]
    )


def condition_data(cond: str) -> dict:
    """the movie's ``data`` as XR_loadConditionData3D gives it, for the fields processing reads"""
    ex = os.path.join(str(cond), EXPERIMENT)
    ds = [sorted(glob.glob(os.path.join(ex, f"ch{c}", "DS", "*.tif"))) for c in (1, 2)]
    return {
        "source": os.path.join(ex, "ch1") + os.sep,
        "framerate": FRAME_MS / 1000.0,
        "movieLength": float(NF),
        "framePathsDS": [np.array(p, dtype=object).reshape(-1, 1) for p in ds],
        "zAniso": Z_ANISO,
        "angle": 31.5,
        "objectiveScan": False,
        "reversed": False,
    }


def _py_run(cond, kind: str = "designed"):
    from petakit5d.detection.processing import rotate_tracks, run_track_processing
    from petakit5d.io.matfile import write_mat

    data = condition_data(cond)
    rp = data["source"] + "Analysis" + os.sep
    write_mat(rp + "Detection3D.mat", {"frameInfo": frame_info(kind)}, version="7.3")
    write_mat(rp + "trackedFeatures.mat", {"tracksFinal": tracks_final(kind)}, version="7")
    run_track_processing(data, overwrite=True, results_path=[rp])
    rotate_tracks(data, overwrite=True, results_path=[rp])


# the gap and buffer fits' values and what follows from them (positions, amplitudes, their
# deviations and p-values, MotionAnalysis, the rotated positions) carry C02b's fits' last bits
# (docs/parity/detection.md W11, W12); every other value is equal at any tolerance, and classes,
# sizes, fields and NaN are compared as they are
FIT_REL = 1e-11
FIT_ABS = 1e-10
TRACKS_TOL = tol(
    reason="docs/parity/detection.md W11, W12: the gap and buffer fits' last bits and their "
    "p-values",
    rel=FIT_REL,
    abs=FIT_ABS,
)


MATLAB_RUN = (
    CLEAR
    + "rt = in{1}; ex = [rt '"
    + EXPERIMENT
    + "' filesep]; ana = [ex 'ch1' filesep 'Analysis' "
    "filesep]; frameInfo = in{2}; save([ana 'Detection3D.mat'], 'frameInfo', '-v7.3'); "
    "tracksFinal = in{3}; save([ana 'trackedFeatures.mat'], 'tracksFinal'); "
    "data = XR_loadConditionData3D(ex, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "apath = {[data(1).source 'Analysis' filesep]}; "
    "XR_runTrackProcessing3D(data, 'Overwrite', true, 'TrackerOutput', 'trackedFeatures.mat', "
    "'FileName', 'ProcessedTracks.mat', 'Buffer', [3, 3], 'BufferAll', false, 'FitMixtures', "
    "false, 'WindowSize', [], 'ResultsPath', apath); "
    "rotateTracks3D(data, 'Overwrite', true, 'Crop', false, 'ResultsPath', apath); out = {};"
)


def _case(kind: str) -> Case:
    return Case(
        f"run.{kind}",
        code=MATLAB_RUN,
        args=[InPath("cond"), In("fi"), In("tf")],
        nargout=0,
        python=_py_run,
        py_args=[InPath("cond"), kind],
        inputs={
            "cond": Tree(functools.partial(make_condition, kind=kind), "cond"),
            "fi": Value(lambda: _mstruct(frame_info(kind))),
            "tf": Value(lambda: _mstruct(tracks_final(kind))),
        },
        tree=tree_spec(
            ignore=(r"Detection3D\.mat$", r"trackedFeatures\.mat$"),
            per_file={
                "*Tracks.mat": TRACKS_TOL,
                # text: the amplitudes of gap frames carry the fits' last bits (W15); the writer
                # itself is checked byte for byte on MATLAB's tracks (test_amira_on_matlab_tracks)
                "*.am": ignore("docs/parity/detection.md W15: gap frames' values carry W11"),
            },
        ),
        marks=("slow",),
    )


CASES = [_case("designed"), _case("compound")]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_amira_on_matlab_tracks(case):
    """the Amira writer on MATLAB's processed tracks of a case gives MATLAB's Amira files byte
    for byte (needs the golden store: MATLAB's files)"""
    from petakit5d.detection import amira
    from petakit5d.io.matfile import read_mat

    store = oracle.Store.from_env()
    if store is None or not store.has(case.key()):
        pytest.skip(f"needs MATLAB's files of {case.id} (the golden store)")
    ana = store.case_dir(case.key()) / "files" / "cond" / EXPERIMENT / "ch1" / "Analysis"
    pt = read_mat(ana / "ProcessedTracks.mat", "tracks")["tracks"]
    pt = [pt] if isinstance(pt, dict) else list(pt)
    keep = [t for t in pt if int(t["catIdx"]) in (1, 2, 3, 4)]
    for k in range(1, NF + 1):
        want = (ana / "Amira" / "Tracking" / f"tracking_t_{k:04d}.am").read_bytes().decode()
        assert amira.frame_text(keep, k, [1.0, 1.0, -Z_ANISO], SIZE[0] * Z_ANISO) == want, k


def test_amira_write_tracks(tmp_path):
    from petakit5d.detection import amira

    tr = {
        "x": np.array([2.0, 3.0]),
        "y": np.array([5.0, 6.0]),
        "z": np.array([1.0, 1.0]),
        "A": np.array([11.0, 12.0]),
        "start": 1.0,
        "end": 2.0,
        "lifetime_s": 5.6,
        "catIdx": 1.0,
    }
    out = amira.write_tracks(str(tmp_path / "Amira" / "tracking"), [tr], 3, 2.0, 4)
    assert [os.path.basename(p) for p in out] == [f"tracking_t_000{k}.am" for k in (1, 2, 3)]
    txt = open(out[1]).read()
    # FlipZ: z' = (z - 1) * -zAniso + slices * zAniso; one-row vectors index as matrices do
    assert "\n@1\n1 4 8\n2 5 8\n" in txt and "\n@8\n10\n11\n" in txt
    out = amira.write_tracks(str(tmp_path / "flat" / "tracking"), [tr], 2, 2.0, 4, flip_z=False)
    txt = open(out[1]).read()
    assert "\n@1\n1 4 0\n2 5 0\n" in txt  # z' = (z - 1) * zAniso, no offset


# ---------------------------------------------------------------------------------------------
# the measured rules on small inputs (no MATLAB)
# ---------------------------------------------------------------------------------------------


def test_existing_results_are_kept(tmp_path):
    """without Overwrite a movie whose result exists is left alone (nothing else is read)"""
    from petakit5d.detection.processing import rotate_tracks, run_track_processing

    rp = str(tmp_path) + os.sep
    names = ("ProcessedTracks.mat", "RotatedTracks.mat")
    for f in names:
        (tmp_path / f).write_bytes(b"kept")
    data = {"source": str(tmp_path), "movieLength": 1.0, "framerate": 1.0, "zAniso": 1.0}
    run_track_processing(data, results_path=[rp])
    rotate_tracks(data, results_path=[rp])
    assert all((tmp_path / f).read_bytes() == b"kept" for f in names)
    with pytest.raises(NotImplementedError):
        run_track_processing(data, results_path=[rp], buffer=(5, 5))
    with pytest.raises(NotImplementedError):
        rotate_tracks(data, results_path=[rp], crop=True)


def test_the_wrappers_stages(tmp_path):
    """XR_runDetTrack3d's wrapper finds track processing and the rotation here and calls them
    with its keywords (its defaults; ``parfor`` as on Windows)"""
    from petakit5d.detection import processing, wrapper

    assert wrapper._stage("run_track_processing") is processing.run_track_processing
    assert wrapper._stage("rotate_tracks") is processing.rotate_tracks
    rp = str(tmp_path) + os.sep
    for f in ("ProcessedTracks.mat", "RotatedTracks.mat"):
        (tmp_path / f).write_bytes(b"kept")
    data = {"source": str(tmp_path), "movieLength": 1.0, "framerate": 1.0, "zAniso": 1.0}
    for method in ("lock_distributed", "parfor"):
        wrapper._call(
            "run_track_processing",
            data,
            overwrite=False,
            tracker_output="trackedFeatures.mat",
            file_name="ProcessedTracks.mat",
            buffer=(3, 3),
            buffer_all=False,
            fit_mixtures=False,
            window_size=None,
            results_path=[rp],
            parallel_method=method,
        )
    wrapper._call("rotate_tracks", data, overwrite=False, crop=False, results_path=[rp])
    assert (tmp_path / "ProcessedTracks.mat").read_bytes() == b"kept"


def _same(u, v) -> bool:
    """equal as read back: the same keys in the same order, classes, sizes and bits"""
    if isinstance(v, dict):
        return isinstance(u, dict) and list(u) == list(v) and all(_same(u[k], v[k]) for k in v)
    a, b = np.asarray(u), np.asarray(v)
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


@pytest.mark.parametrize("n", [0, 1, 2, 3])
def test_rotated_tracks_file(tmp_path, n):
    """RotatedTracks.mat: ``rotated_tracks`` of the processed tracks, ``processingInfo`` left
    out; two or more tracks are rotated in a copy of the file (a struct array of references),
    one track (a scalar struct) and none the plain way, and both read back alike"""
    from petakit5d.detection import processing as P
    from petakit5d.detection.processing import rotate_tracks, rotated_tracks
    from petakit5d.io.matfile import StructArray, read_mat, write_mat

    frame = str(tmp_path / "f1.tif")
    synth.write_tiff(frame, np.zeros((9, 20, 31), dtype=np.uint16))
    rng = np.random.default_rng(7)

    def buf(k):
        return {
            "t": np.arange(1.0, k + 1).reshape(1, -1),
            **{a: rng.random((2, k)) * 30 for a in "xyzA"},
        }

    rows = []
    for i in range(n):
        tr = {"t": np.arange(1.0, 6).reshape(1, -1), **{a: rng.random((2, 5)) * 30 for a in "xyzA"}}
        tr["x"][1, 2] = np.nan  # a channel without a value
        tr["startBuffer"] = buf(3)
        tr["endBuffer"] = np.zeros((0, 0)) if i == 1 else buf(3)
        tr["catIdx"] = float(i + 1)
        rows.append(tr)
    rp = str(tmp_path) + os.sep
    st = (
        StructArray((1, n), tuple(rows[0]), rows)
        if n
        else StructArray.empty(P.TRACK_FIELDS, (1, 0))
    )
    write_mat(rp + "ProcessedTracks.mat", {"tracks": st, "processingInfo": {"procFlag": 1.0}})
    data = {
        "source": rp,
        "framePathsDS": [np.array([frame], dtype=object).reshape(-1, 1)],
        "angle": 31.5,
        "zAniso": 2.5,
        "objectiveScan": False,
        "reversed": False,
    }
    assert P._has_tracks_struct(rp + "ProcessedTracks.mat") == (n > 1)  # the path taken
    rotate_tracks(data, results_path=[rp])
    out = read_mat(rp + "RotatedTracks.mat")
    assert list(out) == ["tracks"]
    got = out["tracks"]
    got = [got] if isinstance(got, dict) else list(got)
    pt = read_mat(rp + "ProcessedTracks.mat", "tracks")["tracks"]
    pt = [pt] if isinstance(pt, dict) else list(pt)
    want = rotated_tracks(pt, (9, 20, 31), 31.5, 2.5)
    assert len(got) == n
    assert all(_same(g, w) for g, w in zip(got, want, strict=True))
    if n:
        assert "Y" in got[0]["endBuffer"] and "Y" not in got[0]["startBuffer"]
    else:  # no tracks: an empty struct array of the processing's fields
        assert out["tracks"].shape == (1, 0) and out["tracks"].fields == P.TRACK_FIELDS


def test_cohort_thresholds_of_a_probe():
    """the lifetime cohorts of a designed probe (12 plain tracks of 8 to 19 frames, amplitudes 150
    to 425, frames of 2 s, 40 frames) against MATLAB's thresholds recorded from that probe (the
    outputs of its interpln calls)"""
    from petakit5d.detection import processing as P

    tracks = [
        {"lifetime_s": 2.0 * (8 + i), "A": np.full((2, 8 + i), 150.0 + 25.0 * i)} for i in range(12)
    ]
    lo, hi = P.cohort_thresholds(tracks, list(range(12)), 2.0, 40)
    want_lo = [425.0] * 2 + [150.0] * 10 + [151.25, 176.25] + [200.0 + 25.0 * k for k in range(10)]
    want_lo += [425.0] * 16
    want_hi = [425.0] * 2 + [150.0 + 25.0 * k for k in range(10)] + [398.75, 423.75]
    want_hi += [425.0] * 26
    assert lo.tolist() == want_lo
    assert hi.tolist() == want_hi


def _tr(gaps: str, nseg: int = 1) -> dict:
    gv = np.array([c == "G" for c in gaps]).reshape(1, -1)
    f = np.array([np.nan if c == "|" else 1.0 for c in gaps]).reshape(1, -1)
    return {"gapVect": gv, "f": f, "nSeg": float(nseg)}


@pytest.mark.parametrize(
    "gaps, status",
    [
        (".G.", [4]),
        (".GG.", [5]),
        (".GG..", [5]),
        ("..GG..", [4]),
        (".G.GG...", [4, 5]),
        ("...G..GG....G.GG.", [4, 4, 4, 5]),
        ("..GG.|.G..", [5, 4]),
        ("..GG..|.G..", [4, 4]),
        (".GG|..G..", [5, 4]),
    ],
)
def test_gap_status(gaps, status):
    from petakit5d.detection import processing as P

    tr = _tr(gaps)
    runs = P.gap_runs(tr)
    assert P.gap_status(tr, runs).tolist() == status


def test_gap_runs_one_based():
    from petakit5d.detection import processing as P

    runs = P.gap_runs(_tr(".GG...G."))
    assert [r.tolist() for r in runs] == [[[2.0, 3.0]], [[7.0]]]


def test_frame_paths():
    """``framePathsDS`` in the forms the inputs arrive in (loadmat's object arrays, the
    registered wrapper's lists and string arrays): file names per channel, column-major"""
    from petakit5d.detection.processing import _channel_paths, _paths

    assert _paths("a.tif") == ["a.tif"]
    assert _paths(np.array(["a.tif", "b.tif"])) == ["a.tif", "b.tif"]
    assert _paths(np.array([list("a.tif")])) == ["a.tif"]  # a char row
    assert _paths([["a", "b"], "c"]) == ["a", "b", "c"]
    assert _paths(np.array([["a", "b"], ["c", "d"]], dtype=object)) == ["a", "c", "b", "d"]
    col = [np.array(["a1", "a2"], dtype=object).reshape(-1, 1), np.array(["b1"], dtype=object)]
    assert _channel_paths(col) == [["a1", "a2"], ["b1"]]
    assert _channel_paths(np.array(["a1", "a2"], dtype=object)) == [["a1", "a2"]]
    assert _channel_paths(["a1", "a2"]) == [["a1", "a2"]]
    cell = np.empty((1, 2), dtype=object)  # loadmat's 1 x nCh cell of columns
    cell[0, 0], cell[0, 1] = col
    assert _channel_paths(cell) == [["a1", "a2"], ["b1"]]


@pytest.mark.parametrize(
    ("sep", "source", "want"),
    [
        ("/", "/data/ex1/ch1/", "/data/ex1/ch1/Analysis/"),
        ("\\", "D:\\data\\ex1\\ch1\\", "D:\\data\\ex1\\ch1\\Analysis\\"),
    ],
)
def test_default_results_folder(monkeypatch, tmp_path, sep, source, want):
    """the default ResultsPath is MATLAB's [source 'Analysis' filesep], with each platform's
    filesep (both forms, on any platform); without ResultsPath both functions use it"""
    from petakit5d.core import paths
    from petakit5d.detection.processing import (
        _default_results,
        rotate_tracks,
        run_track_processing,
    )

    monkeypatch.setattr(paths, "filesep", sep)
    assert _default_results([{"source": source}, {"source": source}]) == [want, want]
    monkeypatch.setattr(paths, "filesep", os.sep)
    ana = tmp_path / "Analysis"
    ana.mkdir()
    for f in ("ProcessedTracks.mat", "RotatedTracks.mat"):
        (ana / f).write_bytes(b"kept")
    data = {"source": str(tmp_path) + os.sep, "movieLength": 1.0, "framerate": 1.0}
    run_track_processing(data)
    rotate_tracks(data)
    assert (ana / "ProcessedTracks.mat").read_bytes() == b"kept"


def test_border_and_peaks_rules():
    """a track goes when any position (rounded) of any channel is within the margin of a
    border; multiple peaks are read from channel 1's amplitudes, for single visibility-1 tracks
    over 30 s"""
    from petakit5d.detection.processing import has_multiple_peaks, inside_volume

    size = (90, 110, 24)  # (y, x, z)
    x = np.array([[9.0, 8.0, 7.0, 7.0, 7.0], [20.0, 20.0, 20.0, 20.0, np.nan]])
    tr = {"x": x, "y": np.full((2, 5), 40.0), "z": np.full((2, 5), 12.0)}
    assert inside_volume(tr, 6, size)
    tr["x"][0, 3] = 6.5  # rounds to 7 (MATLAB's round): 7 - 6 >= 1
    assert inside_volume(tr, 6, size)
    tr["x"][0, 3] = 6.49  # one position at 6: the track goes
    assert not inside_volume(tr, 6, size)
    tr["x"][0, 3] = 7.0
    tr["x"][1, 3] = 104.6  # 105 + 6 > 110, in channel 2
    assert not inside_volume(tr, 6, size)

    peaks = [200.0 + 150.0 * ((k // 4) % 2 == 0) for k in range(20)]
    flat = [300.0] * 20

    def t(a1, a2, lifetime=40.0, vis=1.0, nseg=1.0):
        return {"A": np.array([a1, a2]), "nSeg": nseg, "visibility": vis, "lifetime_s": lifetime}

    assert has_multiple_peaks(t(peaks, flat))
    assert not has_multiple_peaks(t(flat, peaks))
    assert not has_multiple_peaks(t(peaks, flat, lifetime=30.0))
    assert not has_multiple_peaks(t(peaks, flat, vis=2.0))
    assert not has_multiple_peaks(t(peaks, flat, nseg=2.0))


def test_join_positions():
    """a join's positions: detection k (1-based) of a frame, and the interpolation at the
    overlap frame (the ends themselves at the ends)"""
    from petakit5d.detection.processing import _interp_xy, _xy

    fi = [{"x": [[1.0, 2.0]], "y": [[3.0, 4.0]]}, {"x": [[5.0, 6.0, 7.0]], "y": [[8.0, 9.0, 10.0]]}]
    assert _xy(fi, 2, 1).tolist() == [5.0, 8.0]
    assert _xy(fi, 2, 3).tolist() == [7.0, 10.0]
    assert _xy(fi, 1, 2).tolist() == [2.0, 4.0]
    v0, v1 = [1.0, 2.0], [3.0, 6.0]
    assert _interp_xy(10, v0, 12, v1, 11).tolist() == [2.0, 4.0]
    assert _interp_xy(11, v0, 12, v1, 11).tolist() == v0
    assert _interp_xy(10, v0, 11, v1, 11).tolist() == v1


def test_interpln():
    from petakit5d.detection import processing as P

    x = [4.0, 6.0, 12.0]
    y = [150.0, 175.0, 425.0]
    out = P.interpln(x, y, [0.0, 2.0, 4.0, 5.0, 6.0, 8.0, 12.0, 14.0]).ravel()
    assert out.tolist() == [425.0, 425.0, 150.0, 150.0, 175.0, 175.0, 425.0, 425.0]


def test_locmax1d():
    from petakit5d.detection import processing as P

    x = np.array([5, 4, 3, 2, 3, 4, 9, 4, 3, 2, 1, 0], dtype=float)
    assert P.locmax1d(x, 7).ravel().tolist() == [7.0]
    assert P.locmax1d(np.arange(10.0)[::-1], 7).size == 0


def test_linear_indexing_quirk():
    from petakit5d.detection import processing as P

    M = np.array([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]])
    assert P._linear(M, [1, 2, 3, 4]).tolist() == [1.0, 10.0, 2.0, 20.0]


def test_dense_gaps():
    from petakit5d.detection import processing as P

    assert P._dense_gaps(np.array([0, 1, 1, 0], bool))
    assert not P._dense_gaps(np.array([0, 1, 1, 0, 0], bool))
    assert P._dense_gaps(np.array([0, 1, 1, 0, 1, 1, 0, 0], bool))
    assert not P._dense_gaps(np.array([0, 1, 1, 0, 1, 0, 0, 0], bool))


def test_rotation_matrix_and_points():
    from petakit5d.detection import processing as P

    T = P.rotation_matrix((16, 48, 64), 31.5, Z_ANISO)
    assert T[1] == [0.0, 1.0, 0.0, 0.0]
    x, y, z = P._forward(T, np.array([14.0]), np.array([14.0]), np.array([9.7]))
    assert y.tolist() == [14.0]
    # the volume's centre goes to the rotated box's centre
    cx, cz = 32.5, 8.5
    xc, _, zc = P._forward(T, np.array([cx]), np.array([1.0]), np.array([cz]))
    assert abs(xc[0] - 38.5) < 1e-12 and abs(zc[0] - 34.5) < 1e-12
    tracks = [
        {
            "x": np.array([[14.0, np.nan]]),
            "y": np.array([[3.0, 4.0]]),
            "z": np.array([[9.7, 2.0]]),
            "catIdx": 1.0,
        }
    ]
    out = P.rotated_tracks(tracks, (16, 48, 64), 31.5, Z_ANISO)
    assert out[0]["y"] is tracks[0]["y"] and np.isnan(out[0]["x"][0, 1])
    assert tracks[0]["x"][0, 0] == 14.0  # the input is not changed


def test_amira_frame_text():
    from petakit5d.detection import amira

    tr = {
        "x": np.array([[2.0, 3.0, 4.0], [2.5, 3.5, 4.5]]),
        "y": np.array([[5.0, 6.0, 7.0], [5.0, 6.0, 7.0]]),
        "z": np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]),
        "A": np.array([[11.0, 12.0, 13.0], [1.0, 1.0, 1.0]]),
        "start": 2.0,
        "end": 4.0,
        "lifetime_s": 8.4,
        "catIdx": 1.0,
    }
    assert "define EDGE 0 \n" in amira.frame_text([tr], 2, [1, 1, 1], 0.0)
    txt = amira.frame_text([tr], 3, [1, 1, -2.0], 10.0)
    assert txt.startswith(
        "# AmiraMesh 3D ASCII 2.0\n\ndefine VERTEX 2\ndefine EDGE 1 \ndefine POINT 2\n"
    )
    # the edge's end and its points index the channel x frame matrices linearly (MATLAB's)
    assert "\n@1\n1 4 10\n1.5 4 10\n" in txt
    assert "\n@4\n1 4 10\n1.5 4 10\n" in txt
    assert "\n@8\n10\n11\n" in txt and txt.endswith("\n@9\n1\n")


def test_matlab_track_classes():
    from petakit5d.detection import processing as P

    tr = {k: np.zeros((2, 3)) for k in P.TRACK_FIELDS}
    tr.update(
        t=np.zeros((1, 3)),
        f=np.zeros((1, 3)),
        gapVect=np.zeros((1, 3), bool),
        gapIdx=[np.array([2.0])],
    )
    tr.update(
        nSeg=1.0,
        visibility=2.0,
        lifetime_s=8.4,
        start=1.0,
        end=3.0,
        catIdx=3.0,
        hasMultiplePeaks=False,
    )
    tr.update(startBuffer=np.zeros((0, 0)), endBuffer=np.zeros((0, 0)))
    tr.update(
        MotionAnalysis={
            "totalDisplacement": 1.0,
            "MSD": np.zeros((1, 10)),
            "MSDstd": np.zeros((1, 10)),
        }
    )
    tr.update(significantMaster=np.ones((2, 1), bool), significantSlave=np.zeros((2, 1), bool))
    tr.update(significantVsBackground=np.ones((2, 3), bool))
    m = P.matlab_track(tr)
    assert list(m) == list(P.TRACK_FIELDS)
    assert m["gapVect"].dtype == bool and m["hasMultiplePeaks"].dtype == bool
    assert m["gapIdx"][0].shape == (1, 1) and isinstance(m["catIdx"], float)
