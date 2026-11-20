"""``XR_runTrackProcessing3D`` and ``rotateTracks3D``: the tracker's tracks
(``trackedFeatures.mat``) turned into MATLAB's processed tracks (``ProcessedTracks.mat``): the
detections' values along each track, short compound-track branches decoupled, tracks too near
the volume's faces dropped, gaps filled by Gaussian fits at interpolated positions, the buffer
frames before and after each track read out with fits at its first and last positions, the
categories and the per-track statistics; one Amira file per frame (``detection.amira``); and
the tracks rotated into the deskewed-rotated frame (``RotatedTracks.mat``).

``run_track_processing`` (``XR_runTrackProcessing3D``) and ``rotate_tracks`` (``rotateTracks3D``)
take a condition's movies and their results folders, and are the stages ``detection.wrapper``'s
``XR_runDetTrack3d`` calls; ``process_movie`` and ``rotated_tracks`` do the same in memory.

New implementations from MATLAB's measured behaviour (task C04b): both files' headers name an
outside author, so their code was not read; the behaviour comes from recordings of the
functions they call (inputs, outputs and caller line numbers), the profiler's executed line
numbers and counts, designed probes and the files they write. ``docs/parity/detection.md``
(rows W11 on) lists what is not bit for bit.

The tracks keep the file's layout and values (``docs/conventions.md``, exceptions 4 and 5): the
fields are channel x frame matrices, compound tracks have their segments side by side with a
NaN column between them, positions and frame numbers are MATLAB's 1-based ones, and the
rotation's matrix is MATLAB's row-vector ``[x y z 1] * T``. Volumes are read as ``(z, y, x)``.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

import numpy as np

from petakit5d.core.dtypes import NAN, mround
from petakit5d.matlabcompat.interp import interp1
from petakit5d.params import Param, register

__all__ = [
    "FIELDS17",
    "TRACK_FIELDS",
    "convert_tracks",
    "process_movie",
    "rotate_tracks",
    "rotated_tracks",
    "run_track_processing",
    "statistics",
    "write_processed_tracks",
]

# the per-frame fields read from the detections, in the processed tracks' order
FIELDS17 = (
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
    "pval_Ar",
    "isPSF",
    "RSS",
    "hval_Ar",
    "hval_AD",
)
# the processed track's fields in MATLAB's order
TRACK_FIELDS = (
    "t",
    "f",
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
    "pval_Ar",
    "isPSF",
    "tracksFeatIndxCG",
    "gapVect",
    "gapStatus",
    "gapIdx",
    "seqOfEvents",
    "nSeg",
    "visibility",
    "lifetime_s",
    "start",
    "end",
    "startBuffer",
    "endBuffer",
    "MotionAnalysis",
    "RSS",
    "hval_Ar",
    "hval_AD",
    "catIdx",
    "hasMultiplePeaks",
    "significantMaster",
    "significantVsBackground",
    "significantSlave",
)


def _mat2(a) -> np.ndarray:
    """a MATLAB value as a 2-D float array (a scalar as 1 x 1)"""
    v = np.asarray(a, dtype=np.float64)
    if v.ndim == 0:
        return v.reshape(1, 1)
    if v.ndim == 1:
        return v.reshape(1, -1)
    return v


# -- the tracker's tracks -------------------------------------------------------------------------


@dataclass
class _Raw:
    """one track of the tracker: feature indices (segments x frames, 0: none), the events
    ``[frame, 1 start | 2 end, segment, partner or NaN]`` and the first frame (1-based)"""

    F: np.ndarray
    E: np.ndarray
    start: int
    order: int = 0
    joined: list = field(default_factory=list)


def raw_tracks(tracks_final) -> list[_Raw]:
    """the tracker's ``tracksFinal`` (a sequence of dicts) as raw tracks"""
    out = []
    for i, t in enumerate(tracks_final):
        F = _mat2(t["tracksFeatIndxCG"]).astype(np.int64)
        E = _mat2(t["seqOfEvents"]).reshape(-1, 4).copy()
        if F.shape[0] > 1:
            E[np.isnan(E)] = NAN  # MATLAB writes a compound track's events anew (measured)
        out.append(_Raw(F, E, int(np.min(E[:, 0])), i))
    return out


def _events_of(E: np.ndarray, seg: int) -> tuple[np.ndarray, np.ndarray]:
    """the start and end event rows of a segment (1-based number)"""
    s = E[(E[:, 2] == seg) & (E[:, 1] == 1)]
    e = E[(E[:, 2] == seg) & (E[:, 1] == 2)]
    return s[0], e[0]


def _span(E: np.ndarray, seg: int) -> tuple[int, int]:
    """a segment's first and last frame from its events (a merge's end is reported one frame
    after its last detection)"""
    s, e = _events_of(E, seg)
    last = int(e[0]) - (0 if np.isnan(e[3]) else 1)
    return int(s[0]), last


def _xy(frame_info, frame: int, k: int) -> np.ndarray:
    fi = frame_info[frame - 1]
    return np.array([float(_mat2(fi["x"])[0, k - 1]), float(_mat2(fi["y"])[0, k - 1])])


def _interp_xy(x0: float, v0, x1: float, v1, xq: float) -> np.ndarray:
    """``interp1([x0 x1], [v0; v1], xq)`` of the positions (x and y, each a column)"""
    g = np.array([x0, x1], dtype=np.float64)
    v = np.array([v0, v1], dtype=np.float64)
    return np.array([float(interp1(g, v[:, j], float(xq))) for j in range(v.shape[1])])


def _join_late_receivers(t: _Raw, frame_info) -> None:
    """a merge into a segment that started one frame before the merge becomes one segment: the
    merging segment's frames, at the overlap frame the detection nearer (in x and y) to the track
    interpolated across it, then the receiver's frames (main:204-255, measured)"""
    while True:
        E = t.E
        nseg = t.F.shape[0]
        hit = None
        for a in range(1, nseg + 1):
            sa, ea = _events_of(E, a)
            if np.isnan(ea[3]):
                continue
            b = int(ea[3])
            sb, _ = _events_of(E, b)
            if int(ea[0]) - int(sb[0]) == 1:
                hit = (a, b, int(ea[0]))
                break
        if hit is None:
            return
        a, b, fm = hit
        o = fm - 1  # the overlap frame, both segments detected
        c = o - t.start
        da, db = int(t.F[a - 1, c]), int(t.F[b - 1, c])
        prev = o - 1
        if prev >= t.start and t.F[a - 1, c - 1] > 0:
            p0 = _xy(frame_info, prev, int(t.F[a - 1, c - 1]))
            x0 = prev
        else:
            p0 = _xy(frame_info, o, da) if da > 0 else None
            x0 = o
        nxt = int(t.F[b - 1, c + 1]) if c + 1 < t.F.shape[1] else 0
        keep = da
        if p0 is not None and nxt > 0 and da > 0 and db > 0:
            q = _interp_xy(x0, p0, o + 1, _xy(frame_info, o + 1, nxt), o)
            d_a = math.hypot(*(_xy(frame_info, o, da) - q))
            d_b = math.hypot(*(_xy(frame_info, o, db) - q))
            keep = db if d_b < d_a else da
        elif da == 0:
            keep = db
        # the joined segment takes the receiver's row and number; the merging segment's row
        # goes and the segments after it move up one (measured)
        row = t.F[a - 1].copy()
        row[c] = keep
        row[c + 1 :] = t.F[b - 1, c + 1 :]
        F = t.F.copy()
        F[b - 1] = row
        F = np.delete(F, a - 1, axis=0)
        # events: the receiver's start row takes the merging segment's start frame (its
        # partner stays); the merging segment's start and merge rows go; the merging segment
        # becomes the receiver everywhere (a receiver split from it then names itself); the
        # rows in frame order (stable)
        E = E.copy()
        a_start = _events_of(E, a)[0][0]
        E[(E[:, 2] == b) & (E[:, 1] == 1), 0] = a_start
        E = E[~(E[:, 2] == a)].copy()
        E[E[:, 3] == a, 3] = b
        E[E[:, 2] > a, 2] -= 1
        E[E[:, 3] > a, 3] -= 1
        E = E[np.argsort(E[:, 0], kind="stable")]
        # the segments in the order of their starts (stable), as the tracker orders them
        n = F.shape[0]
        starts = [_events_of(E, k)[0][0] for k in range(1, n + 1)]
        order = sorted(range(n), key=lambda i: starts[i])
        if order != list(range(n)):
            new = {old + 1: i + 1 for i, old in enumerate(order)}
            F = F[order]
            E[:, 2] = [new[int(v)] for v in E[:, 2]]
            E[:, 3] = [new[int(v)] if not np.isnan(v) else v for v in E[:, 3]]
        t.F, t.E = F, E
        t.joined.append((a, b, fm))


def _remove_short_segments(t: _Raw) -> None:
    """short branches go (measured): a segment of one frame attached to another; a segment of
    two or three frames that splits off and dies before the track's last frame, is born after
    the track's first frame and merges, or splits off and merges back into the same segment (a
    loop); a segment between two different segments stays. The other events keep their
    segment numbers."""
    E = t.E
    nseg = t.F.shape[0]
    first = int(np.min(E[:, 0]))
    last = t.start + t.F.shape[1] - 1
    drop = []
    for k in range(1, nseg + 1):
        s, e = _events_of(E, k)
        ps = int(s[3]) if not np.isnan(s[3]) else 0
        pe = int(e[3]) if not np.isnan(e[3]) else 0
        if not (ps or pe):
            continue
        f0, f1 = _span(E, k)
        n = f1 - f0 + 1
        if n == 1:
            drop.append(k)
        elif n <= 3:
            if ps and pe:
                if ps == pe:
                    drop.append(k)
            elif ps and f1 < last:
                drop.append(k)
            elif pe and f0 > first:
                drop.append(k)
    if not drop:
        return
    keep = [k for k in range(1, nseg + 1) if k not in drop]
    t.F = t.F[[k - 1 for k in keep]]
    t.E = E[~np.isin(E[:, 2], drop)].copy()
    t.seg_ids = keep


# -- conversion -----------------------------------------------------------------------------------


def convert_tracks(tracks_final, frame_info, framerate: float) -> list[dict]:
    """the tracker's tracks as MATLAB's processed tracks before the fits: single-frame tracks
    dropped, compound tracks decoupled, the detections' 17 fields per frame (z as the
    detections have it; ``reverse_z`` turns it at the end), frames, times, gaps, events, span
    and lifetime"""
    nch = _mat2(frame_info[0]["x"]).shape[0] if len(frame_info) else 1
    out = []
    for t in raw_tracks(tracks_final):
        if t.F.shape[1] == 1:
            continue
        t.seg_ids = list(range(1, t.F.shape[0] + 1))
        if t.F.shape[0] > 1:
            _join_late_receivers(t, frame_info)
            t.seg_ids = list(range(1, t.F.shape[0] + 1))
            _remove_short_segments(t)
            t.E = t.E[np.argsort(t.E[:, 0], kind="stable")]
        out.append(_convert_one(t, frame_info, framerate, nch))
    return out


def _convert_one(t: _Raw, frame_info, framerate: float, nch: int) -> dict:
    ncol = t.F.shape[1]
    cols_f, cols_gap, vals = [], [], {k: [] for k in FIELDS17}
    for r, seg in enumerate(t.seg_ids):
        if r:
            cols_f.append(np.nan)
            cols_gap.append(False)
            for k in FIELDS17:
                vals[k].append(np.full(nch, NAN))
        f0, f1 = _span(t.E, seg)
        for f in range(f0, f1 + 1):
            d = int(t.F[r, f - t.start])
            cols_f.append(float(f))
            cols_gap.append(d == 0)
            fi = frame_info[f - 1]
            for k in FIELDS17:
                if d == 0:
                    vals[k].append(np.full(nch, NAN))
                else:
                    vals[k].append(_mat2(fi[k])[:, d - 1].astype(np.float64))
    f = np.array(cols_f, dtype=np.float64).reshape(1, -1)
    tr = {
        "t": np.where(np.isnan(f), NAN, (f - 1) * framerate),
        "f": np.where(np.isnan(f), NAN, f),
    }
    for k in FIELDS17:
        tr[k] = np.stack(vals[k], axis=1)
    tr["tracksFeatIndxCG"] = t.F.astype(np.float64)
    tr["gapVect"] = np.array(cols_gap, dtype=bool).reshape(1, -1)
    tr["seqOfEvents"] = t.E
    tr["nSeg"] = float(t.F.shape[0])
    tr["start"] = float(t.start)
    tr["end"] = float(t.start + ncol - 1)
    tr["lifetime_s"] = ncol * framerate
    tr["_order"] = t.order
    return tr


# -- the border filter and the classes -----------------------------------------------------------


def border_margin(sigma) -> int:
    """the border filter's margin: ``ceil(4 sigma_xy)`` of the first channel (measured)"""
    s = _mat2(sigma)
    return int(math.ceil(4 * float(s[0, 0])))


def inside_volume(tr: dict, w: int, size_yxz) -> bool:
    """whether every position of every channel keeps ``round(p) - w >= 1`` and
    ``round(p) + w <= n`` in x, y and z (``size_yxz``: the DS frame's MATLAB size)"""
    ny, nx, nz = (int(v) for v in size_yxz)
    for k, n in (("x", nx), ("y", ny), ("z", nz)):
        v = tr[k]
        v = v[~np.isnan(v)]
        if v.size == 0:
            continue
        r = mround(v)
        if np.any(r - w < 1) or np.any(r + w > n):
            return False
    return True


def visibility(tr: dict, movie_length: int, buffer=(3, 3)) -> float:
    """3: the whole movie; 1: room for both buffers; 2: otherwise (measured)"""
    s, e = int(tr["start"]), int(tr["end"])
    if s == 1 and e == movie_length:
        return 3.0
    if s > buffer[0] and e <= movie_length - buffer[1]:
        return 1.0
    return 2.0


def reverse_z(tr: dict, n_slices: int) -> None:
    """MATLAB's processed tracks hold z reversed: ``n_slices + 1 - z`` (the track's frames and
    its buffers; measured bit for bit), done after the fits"""
    tr["z"] = (n_slices + 1) - tr["z"]
    for b in ("startBuffer", "endBuffer"):
        if isinstance(tr.get(b), dict):
            tr[b]["z"] = (n_slices + 1) - tr[b]["z"]


def sort_by_lifetime(tracks: list[dict]) -> list[dict]:
    """MATLAB's stable sort by lifetime, longest first"""
    order = sorted(range(len(tracks)), key=lambda i: -tracks[i]["lifetime_s"])
    return [tracks[i] for i in order]


# -- the fits at fixed positions (gaps and buffers) ----------------------------------------------

# the fields a fit gives a gap or buffer frame, in MATLAB's order
FIT_FIELDS = (
    "x",
    "y",
    "z",
    "A",
    "c",
    "A_pstd",
    "c_pstd",
    "sigma_r",
    "SE_sigma_r",
    "pval_Ar",
    "hval_AD",
)


@dataclass
class _Point:
    """a position to read out in one frame and channel: the track, the column (a gap) or the
    buffer index (-1: a gap), the kind ('gap', 'start', 'end') and the position (1-based,
    x, y, and z in the frame's orientation)"""

    track: int
    col: int
    kind: str
    x: float
    y: float
    z: float


def window_half_sizes(sigma) -> tuple[int, int]:
    """the read-out windows' half sizes: ``ceil(2 sigma)`` in xy and z (measured: 7 x 7 x 7 and
    9 x 9 x 7 windows for the demo's two channels)"""
    s = np.asarray(sigma, dtype=np.float64).ravel()
    return int(math.ceil(2 * s[0])), int(math.ceil(2 * s[1]))


def gap_positions(tr: dict) -> dict:
    """per gap column (index into the concatenated vectors), per channel the interpolated
    position (x, y, z): interp1 between the detected frames around the gap run (main:508-510)"""
    f = tr["f"].ravel()
    gv = tr["gapVect"].ravel()
    nch = tr["x"].shape[0]
    out = {}
    n = f.size
    i = 0
    while i < n:
        if not gv[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and gv[j + 1]:
            j += 1
        a, b = i - 1, j + 1  # the detected columns around the run
        cols = np.arange(i, j + 1)
        for ch in range(nch):
            pos = []
            for k in ("x", "y", "z"):
                v = tr[k][ch]
                pos.append(interp1(f[[a, b]], v[[a, b]], f[cols]))
            for m, c in enumerate(cols):
                out.setdefault(int(c), {})[ch] = (pos[0][m], pos[1][m], pos[2][m])
        i = j + 1
    return out


def mask_labels(mask) -> np.ndarray:
    """the detection mask's 26-connected components (MATLAB's bwconncomp and labelmatrix; only
    equal or different labels matter here)"""
    from scipy import ndimage

    lab, _ = ndimage.label(np.asarray(mask) > 0, structure=np.ones((3, 3, 3), dtype=bool))
    return lab


def read_windows(vol, labels, points, sigma, A_map, c_map):
    """per point the window ``[round(p) - w, round(p) + w]`` of the frame (z, y, x), NaN on the
    voxels of the detection mask's components other than the one at the rounded position, and
    the starting parameters ``[x y z A s sz c]`` in the window's 0-based coordinates, A and c
    from the amplitude maps at the rounded position"""
    wxy, wz = window_half_sizes(sigma)
    s, sz = (float(v) for v in np.asarray(sigma, dtype=np.float64).ravel()[:2])
    wins, p0 = [], np.empty((len(points), 7))
    corner = np.empty((len(points), 3), dtype=np.int64)
    for i, q in enumerate(points):
        rx, ry, rz = (int(mround(v)) for v in (q.x, q.y, q.z))
        x0, y0, z0 = rx - wxy, ry - wxy, rz - wz  # 1-based first voxel
        w = vol[z0 - 1 : rz + wz, y0 - 1 : ry + wxy, x0 - 1 : rx + wxy].copy()
        if labels is not None:
            lab = labels[z0 - 1 : rz + wz, y0 - 1 : ry + wxy, x0 - 1 : rx + wxy]
            own = labels[rz - 1, ry - 1, rx - 1]
            w[(lab > 0) & (lab != own)] = NAN
        wins.append(w)
        corner[i] = (x0, y0, z0)
        p0[i] = (
            q.x - x0,
            q.y - y0,
            q.z - z0,
            A_map[rz - 1, ry - 1, rx - 1],
            s,
            sz,
            c_map[rz - 1, ry - 1, rx - 1],
        )
    return wins, p0, corner


K_LEVEL = 1.959963984540054  # norminv(0.975, 0, 1), as MATLAB returns it (main:132)


def fit_windows(wins, p0, corner, sigma, *, num_threads: int | None = None) -> dict:
    """interpTrack's read-out of each window: the free fit (``XR_fitGaussian3D``, mode
    'xyzAc'), kept when its position is within 2 voxels of the window's centre in x, y and z,
    else the fit of amplitude and background at the given position (the ``fitGaussian3D`` mex,
    mode 'Ac'); then the fields ``FIT_FIELDS`` (positions 1-based, z as the frame has it) and
    the amplitude's one-sided Welch test against kLevel residual sigmas (measured)"""
    from petakit5d.detection import fit as _fit
    from petakit5d.detection import stats as _st

    k = len(wins)
    out = {f: np.full(k, NAN) for f in FIT_FIELDS}
    if k == 0:
        return out
    wxy, wz = window_half_sizes(sigma)
    R = _fit.fit_cd(wins, p0, num_threads=num_threads)
    ctr = np.array([wxy, wxy, wz], dtype=np.float64)
    kept = np.all(np.abs(R.prm[:, :3] - ctr) <= 2, axis=1)
    prm = R.prm.copy()
    A_pstd = R.pstd[:, 3].copy()
    c_pstd = R.pstd[:, 6].copy()
    sd = R.std.copy()
    had = R.h_ad.copy()
    n = R.n.astype(np.float64)
    ref = np.nonzero(~kept)[0]
    if ref.size:
        M = _fit.fit_gaussian3d_mex([wins[i] for i in ref], p0[ref], "Ac", num_threads=num_threads)
        prm[ref] = M.prm
        A_pstd[ref] = M.prm_std[:, 0]
        c_pstd[ref] = M.prm_std[:, 1]
        sd[ref] = M.std
        had[ref] = M.h_ad
        n[ref] = M.n
    out["x"] = corner[:, 0] + prm[:, 0]
    out["y"] = corner[:, 1] + prm[:, 1]
    out["z"] = corner[:, 2] + prm[:, 2]
    out["A"] = prm[:, 3]
    out["c"] = prm[:, 6]
    out["A_pstd"] = A_pstd
    out["c_pstd"] = c_pstd
    out["sigma_r"] = sd
    SE = sd / np.sqrt(2 * (n - 1))
    out["SE_sigma_r"] = SE
    with np.errstate(invalid="ignore", divide="ignore"):
        q = SE * K_LEVEL
        T = (prm[:, 3] - sd * K_LEVEL) / np.sqrt((A_pstd * A_pstd + q**2) / n)
        df = (n - 1) * (A_pstd * A_pstd + q**2) ** 2 / (A_pstd**4 + q**4)
    out["pval_Ar"] = np.asarray(_st.tcdf(-T, df), dtype=np.float64).ravel()
    out["hval_AD"] = had.astype(np.float64)
    out["_kept"] = kept
    return out


BUFFER_FIELDS = (
    "t",
    "x",
    "y",
    "z",
    "A",
    "c",
    "A_pstd",
    "c_pstd",
    "sigma_r",
    "SE_sigma_r",
    "pval_Ar",
)


def _n_frames(tr: dict, framerate: float) -> int:
    return int(round(float(tr["lifetime_s"]) / framerate))


def init_buffers(tracks: list[dict], framerate: float, buffer=(3, 3)) -> None:
    """visibility-1 tracks get NaN buffers of ``buffer`` frames before and after (fields
    ``BUFFER_FIELDS``, t the frames' times); the others ``[]`` (measured)"""
    for tr in tracks:
        if int(tr["visibility"]) != 1:
            tr["startBuffer"] = np.zeros((0, 0))
            tr["endBuffer"] = np.zeros((0, 0))
            continue
        nch = tr["x"].shape[0]
        s, e = int(tr["start"]), int(tr["end"])
        for key, frames in (
            ("startBuffer", np.arange(s - buffer[0], s)),
            ("endBuffer", np.arange(e + 1, e + 1 + buffer[1])),
        ):
            b = {"t": ((frames - 1) * framerate).reshape(1, -1).astype(np.float64)}
            for f in BUFFER_FIELDS[1:]:
                b[f] = np.full((nch, frames.size), NAN)
            tr[key] = b


def frame_points(
    tracks: list[dict], gaps: list[dict], k: int, ch: int, framerate: float, buffer=(3, 3)
):
    """the positions to read out in frame k (1-based) and channel ch, in MATLAB's order: the
    gap frames (track by track), then the start buffers, then the end buffers of the
    visibility-1 tracks of three frames or more (tracks in the tracker's order)"""
    pts = []
    for j, tr in enumerate(tracks):
        f = tr["f"].ravel()
        for c in np.nonzero(tr["gapVect"].ravel() & (f == k))[0]:
            x, y, z = gaps[j][int(c)][ch]
            pts.append(_Point(j, int(c), "gap", x, y, z))
    for kind in ("start", "end"):
        for j, tr in enumerate(tracks):
            if int(tr["visibility"]) != 1 or _n_frames(tr, framerate) < 3:
                continue
            s, e = int(tr["start"]), int(tr["end"])
            if kind == "start" and s - buffer[0] <= k <= s - 1:
                col, idx = 0, k - (s - buffer[0])
            elif kind == "end" and e + 1 <= k <= e + buffer[1]:
                col = int(np.nonzero(~np.isnan(tr["f"].ravel()))[0][-1])
                idx = k - (e + 1)
            else:
                continue
            pts.append(_Point(j, idx, kind, tr["x"][ch, col], tr["y"][ch, col], tr["z"][ch, col]))
    return pts


def read_out(
    tracks: list[dict],
    read_frame,
    read_mask,
    sigma,
    movie_length: int,
    framerate: float,
    buffer=(3, 3),
    *,
    num_threads: int | None = None,
) -> None:
    """the gap and buffer read-out (main:539-608): per frame the detection mask's components,
    per channel the frame, the amplitude maps (``estGaussianAmplitude3D``) and the fits at the
    frame's positions; ``tracks`` in the tracker's order. ``read_frame(ch, k)`` gives the DS
    frame (z, y, x) as double, ``read_mask(k)`` the detection mask."""
    from petakit5d.detection import detect as _d

    sig = _mat2(sigma)
    nch = sig.shape[0]
    gaps = [gap_positions(tr) for tr in tracks]
    for k in range(1, movie_length + 1):
        pts_ch = [frame_points(tracks, gaps, k, ch, framerate, buffer) for ch in range(nch)]
        if not any(pts_ch):
            continue
        labels = mask_labels(read_mask(k))
        for ch in range(nch):
            pts = pts_ch[ch]
            if not pts:
                continue
            vol = read_frame(ch, k)
            A_map, c_map = _d.est_gaussian_amplitude3d(vol, sig[ch], num_threads=num_threads)
            wins, p0, corner = read_windows(vol, labels, pts, sig[ch], A_map, c_map)
            res = fit_windows(wins, p0, corner, sig[ch], num_threads=num_threads)
            for m, q in enumerate(pts):
                tr = tracks[q.track]
                if q.kind == "gap":
                    for f in FIT_FIELDS:
                        tr[f][ch, q.col] = res[f][m]
                else:
                    b = tr["startBuffer" if q.kind == "start" else "endBuffer"]
                    if "hval_AD" not in b:
                        # the field takes the class of its first value: double from the free
                        # fit, logical from the refit at the given position (measured on all
                        # 3,980 demo buffers)
                        dt = np.float64 if res["_kept"][m] else np.bool_
                        b["hval_AD"] = np.zeros(
                            (nch, buffer[0 if q.kind == "start" else 1]), dtype=dt
                        )
                    for f in FIT_FIELDS[:-1]:
                        b[f][ch, q.col] = res[f][m]
                    b["hval_AD"][ch, q.col] = res["hval_AD"][m]


# -- the statistics -------------------------------------------------------------------------------


def locmax1d(x, w: int) -> np.ndarray:
    """1-based indices of the local maxima of ``x``: elements equal to the maximum of their
    window of ``w`` (odd), the window wholly inside the vector (measured on 494 calls)"""
    x = np.asarray(x, dtype=np.float64).ravel()
    h = int(w) // 2
    out = [i + 1 for i in range(h, x.size - h) if x[i] == np.max(x[i - h : i + h + 1])]
    return np.array(out, dtype=np.float64).reshape(-1, 1)


def has_multiple_peaks(tr: dict) -> bool:
    """single-segment, visibility-1 tracks longer than 30 s: whether the autocorrelation of
    channel 1's amplitude (``conv(A, fliplr(A))`` from lag 0) has a local maximum (window 7;
    measured: 11 frames of 2.8 s and 16 of 2 s are tested, 10 and 15 not)"""
    if int(tr["nSeg"]) != 1 or int(tr["visibility"]) != 1 or not float(tr["lifetime_s"]) > 30:
        return False
    a = tr["A"][0]
    ac = np.convolve(a, a[::-1])
    return locmax1d(ac[a.size - 1 :], 7).size > 0


def gap_runs(tr: dict) -> list[np.ndarray]:
    """the gaps: runs of consecutive gap columns, as 1-based indices into the concatenated
    vectors (MATLAB's gapIdx)"""
    gv = tr["gapVect"].ravel()
    runs = []
    i = 0
    while i < gv.size:
        if gv[i]:
            j = i
            while j + 1 < gv.size and gv[j + 1]:
                j += 1
            runs.append(np.arange(i + 1, j + 2, dtype=np.float64).reshape(1, -1))
            i = j + 1
        else:
            i += 1
    return runs


def gap_status(tr: dict, runs: list[np.ndarray]) -> np.ndarray:
    """per gap 4 (valid) or 5 (invalid): a gap longer than the run of detections before or
    after it (up to the next gap, the track's end or a compound track's separator) is invalid
    (measured on all 6,085 demo gaps)"""
    gv = tr["gapVect"].ravel()
    f = tr["f"].ravel()
    det = ~gv & ~np.isnan(f)
    n = det.size
    st = np.empty(len(runs))
    for k, r in enumerate(runs):
        idx = r.ravel().astype(np.int64) - 1
        a, nb = idx[0] - 1, 0
        while a >= 0 and det[a]:
            nb, a = nb + 1, a - 1
        b, na = idx[-1] + 1, 0
        while b < n and det[b]:
            na, b = na + 1, b + 1
        st[k] = 5.0 if idx.size > min(nb, na) else 4.0
    return st


def _buffers_bad(tr: dict, alpha: float = 0.05) -> bool:
    """channel 1's buffers: significant signal in the two frames nearest the track, or a buffer
    frame brighter than the track's maximum (measured: all but 2 of the demo's 2,171 single
    visibility-1 tracks, docs/parity/detection.md)"""
    b, e = tr["startBuffer"], tr["endBuffer"]
    sp, ep = b["pval_Ar"][0], e["pval_Ar"][0]
    with np.errstate(invalid="ignore"):
        m = np.nanmax(tr["A"][0])
        return bool(
            np.any(sp[1:] < alpha)
            or np.any(ep[:2] < alpha)
            or np.any(b["A"][0] > m)
            or np.any(e["A"][0] > m)
        )


def step_sizes(tr: dict) -> np.ndarray:
    """channel 1's frame-to-frame 3-D steps (every column: gap frames at their fitted
    positions, NaN across a compound track's separators; the input of main:962)"""
    x, y, z = tr["x"][0], tr["y"][0], tr["z"][0]
    dx, dy, dz = np.diff(x), np.diff(y), np.diff(z)
    return np.sqrt(dx * dx + dy * dy + dz * dz)


def interpln(x, y, xi) -> np.ndarray:
    """the lifetime thresholds' interpolation (measured on its recorded calls): at each ``xi``
    the value at the last ``x`` at or before it (``x`` ascending), and ``y``'s last value outside
    ``[x(1), x(end)]``; a column"""
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    xi = np.asarray(xi, dtype=np.float64).ravel()
    out = np.full(xi.size, y[-1] if y.size else NAN)
    if x.size:
        inside = (xi >= x[0]) & (xi <= x[-1])
        out[inside] = y[np.searchsorted(x, xi[inside], side="right") - 1]
    return out.reshape(-1, 1)


def cohort_thresholds(tracks: list[dict], pool, framerate: float, movie_length: int):
    """main:744-775: per lifetime cohort k = 1..movie_length (frames) the 5th and 95th
    percentiles of channel 1's maximum amplitude over the pool's tracks whose lifetime is within
    ``round(10 / framerate)`` frames of k (measured: 4 at 2.8 s, 5 at 2 s), then ``interpln``
    from the cohorts that have tracks (at ``(k - 1) framerate``) to every cohort. Returns
    ``(lo, hi)``, index k - 1 for a track of k frames."""
    from petakit5d.matlabcompat import stats as _s

    w = int(mround(10.0 / framerate))
    life = np.array([_n_frames(tracks[i], framerate) for i in pool], dtype=np.int64)
    amax = np.array([float(np.nanmax(tracks[i]["A"][0])) for i in pool])
    pr = np.full((movie_length, 2), NAN)
    for k in range(1, movie_length + 1):
        v = amax[np.abs(life - k) <= w]
        if v.size:
            pr[k - 1] = np.asarray(_s.prctile(v, [5, 95]), dtype=np.float64)
    t = np.arange(movie_length, dtype=np.float64) * framerate
    ok = ~np.isnan(pr[:, 0])
    return interpln(t[ok], pr[ok, 0], t).ravel(), interpln(t[ok], pr[ok, 1], t).ravel()


def _linear(M, idx) -> np.ndarray:
    """MATLAB's ``M(idx)`` of a channel x frame matrix: linear, column-major 1-based indices, as
    the gap test indexes the per-frame fields (measured; with two channels the indices reach
    channel 2 and other frames)"""
    return np.asarray(M, dtype=np.float64).T.ravel()[np.asarray(idx, dtype=np.int64) - 1]


def _dense_gaps(gv: np.ndarray) -> bool:
    """more than half the frames gaps: four or more in some 7 frames (``conv(gapVect,
    ones(1, 7), 'valid')``, main:799-800), or half or more of a shorter track (measured)"""
    n = gv.size
    if n >= 7:
        return bool(np.convolve(gv.astype(np.float64), np.ones(7), "valid").max() >= 4)
    return bool(gv.sum() >= n / 2)


def statistics(tracks: list[dict], framerate: float, movie_length: int, master: int = 0) -> dict:
    """the per-track statistics and categories (main:609-1021, measured; ``tracks`` in their
    final order): gaps and their status, significance, multiple peaks, then the categories:

    - compound tracks: 5 (6 with an invalid gap), 7, 8 for visibility 1, 2, 3; single tracks of
      visibility 2 and 3: 3 and 4;
    - single tracks of visibility 1: 2 with signal in their buffers or multiple peaks; a track
      with gaps is 1 only when every ``pval_Ar`` its gap indices reach (linearly: ``_linear``)
      is 0.05 or more, channel 1's maximum amplitude reaches the lower threshold of its lifetime
      cohort (``cohort_thresholds``, from the gap-free tracks of category 1) and its gaps are not
      dense (``_dense_gaps``).

    Returns the cohorts' thresholds ``{'lo', 'hi'}`` (per lifetime, index frames - 1)."""
    single_v1 = []
    for i, tr in enumerate(tracks):
        nch, ncol = tr["x"].shape
        runs = gap_runs(tr)
        st = gap_status(tr, runs)
        tr["_runs"] = runs
        if runs:
            tr["gapIdx"] = runs
            tr["gapStatus"] = st.reshape(1, -1)
        elif int(tr["nSeg"]) > 1:  # compound tracks: a 1 x 0 cell and a 1 x 0 double
            tr["gapIdx"] = np.empty((1, 0), dtype=object)
            tr["gapStatus"] = np.zeros((1, 0))
        else:
            tr["gapIdx"] = np.zeros((0, 0))
            tr["gapStatus"] = np.zeros((0, 0))
        tr["significantVsBackground"] = np.ones((nch, ncol), dtype=bool)
        sig = np.zeros((nch, 1), dtype=bool)
        for ch in range(nch):
            sig[ch, 0] = True if ch == master else bool(np.nansum(tr["hval_Ar"][ch]) >= 2)
        tr["significantMaster"] = sig
        tr["significantSlave"] = np.zeros((nch, 1), dtype=bool)
        tr["hasMultiplePeaks"] = has_multiple_peaks(tr)
        vis = int(tr["visibility"])
        if int(tr["nSeg"]) > 1:
            cat = {1: 6.0 if np.any(st == 5) else 5.0, 2: 7.0, 3: 8.0}[vis]
        elif vis == 1:
            if _n_frames(tr, framerate) < 3:
                cat = 1.0  # two frames: NaN buffers, no test
            elif tr["hasMultiplePeaks"] or _buffers_bad(tr):
                cat = 2.0
            else:
                cat = 2.0 if runs else 1.0
                if runs:
                    single_v1.append(i)
        else:
            cat = {2: 3.0, 3: 4.0}[vis]
        tr["catIdx"] = cat
    pool = [
        i
        for i, tr in enumerate(tracks)
        if tr["catIdx"] == 1.0 and int(tr["nSeg"]) == 1 and not tr["_runs"]
    ]
    lo, hi = cohort_thresholds(tracks, pool, framerate, movie_length)
    for i in single_v1:
        tr = tracks[i]
        idx = np.concatenate([r.ravel() for r in tr["_runs"]])
        n = min(_n_frames(tr, framerate), movie_length)
        if (
            np.all(_linear(tr["pval_Ar"], idx) >= 0.05)
            and float(np.nanmax(tr["A"][0])) >= lo[n - 1]
            and not _dense_gaps(tr["gapVect"].ravel())
        ):
            tr["catIdx"] = 1.0
    # main:955-969: the 95th percentile of the tracks' median steps; a track of category 1
    # with five or more steps above it goes to 2 (measured: P8, the demo's two)
    from petakit5d.matlabcompat import stats as _s

    steps = [step_sizes(tr) for tr in tracks]
    med = np.array([float(_s.median(s, omitnan=True)) if s.size else NAN for s in steps])
    p95 = float(np.asarray(_s.prctile(med, 95)).ravel()[0]) if len(tracks) else NAN
    for tr, s in zip(tracks, steps, strict=True):
        if tr["catIdx"] == 1.0 and int(np.sum(s > p95)) >= 5:
            tr["catIdx"] = 2.0
    for tr in tracks:
        tr.pop("_runs", None)
    return {"lo": lo, "hi": hi, "step_p95": p95}


def motion_analysis(tr: dict) -> dict | np.ndarray:
    """single-segment tracks: the distance from the first to the last position and the mean
    squared displacement with its standard deviation at lags 1 to 10 (NaN beyond the track),
    channel 1's positions; compound tracks: []"""
    from petakit5d.matlabcompat import reduce as _r
    from petakit5d.matlabcompat import stats as _s

    if int(tr["nSeg"]) != 1:
        return np.zeros((0, 0))
    x, y, z = tr["x"][0], tr["y"][0], tr["z"][0]
    td = math.sqrt((x[-1] - x[0]) ** 2 + (y[-1] - y[0]) ** 2 + (z[-1] - z[0]) ** 2)
    msd = np.full((1, 10), NAN)
    sd = np.full((1, 10), NAN)
    for lag in range(1, min(10, x.size - 1) + 1):
        d2 = (x[lag:] - x[:-lag]) ** 2 + (y[lag:] - y[:-lag]) ** 2 + (z[lag:] - z[:-lag]) ** 2
        msd[0, lag - 1] = float(np.asarray(_r.mean(d2)).ravel()[0])
        sd[0, lag - 1] = float(np.asarray(_s.std(d2)).ravel()[0])
    return {"totalDisplacement": td, "MSD": msd, "MSDstd": sd}


# -- rotateTracks3D -------------------------------------------------------------------------------


def rotation_matrix(size_zyx, angle: float, z_aniso: float) -> list[list[float]]:
    """rotateTracks3D's affine matrix (MATLAB's row-vector convention, ``[x y z 1] * T``) for a
    DS volume of ``(nz, ny, nx)``: the volume's centre to the origin, z scaled by ``z_aniso``,
    the rotation by ``angle`` degrees in the x-z plane, then the centre of the rotated box
    ``round([nx cos + nz zAniso sin, nx sin + nz zAniso cos])`` (measured on 7 geometries, bit
    for bit; MATLAB's cosine and sine, its 4 x 4 product)."""
    from petakit5d.deskew_rotate.geometry import _chain, _cos, _eye, _round, _sin, _translation

    nz, _ny, nx = (int(v) for v in size_zyx)
    th = float(angle) * math.pi / 180
    c, s = _cos(th), _sin(th)
    S = float(z_aniso)
    wx = _round(nx * c + nz * S * s)
    wz = _round(nx * s + nz * S * c)
    t1 = _translation([-(nx + 1) / 2, 0.0, -(nz + 1) / 2])
    sc = _eye()
    sc[2][2] = S
    r = [[c, 0.0, s, 0.0], [0.0, 1.0, 0.0, 0.0], [-s, 0.0, c, 0.0], [0.0, 0.0, 0.0, 1.0]]
    t2 = _translation([(wx + 1) / 2, 0.0, (wz + 1) / 2])
    return _chain(t1, sc, r, t2)


def _forward(T, x, y, z):
    """``transformPointsForward`` of affine3d: each output coordinate ``K1 x + K2 y + K3 z + K4``
    summed in that order (imlincomb's, measured)"""
    out = []
    for j in range(3):
        out.append(((T[0][j] * x + T[1][j] * y) + T[2][j] * z) + T[3][j])
    return out


def rotated_tracks(tracks: list[dict], size_zyx, angle: float, z_aniso: float) -> list[dict]:
    """``rotateTracks3D`` (``Crop`` false): copies of the processed tracks with x and z (every
    channel, the buffers too) rotated into the deskewed-rotated frame; y and every other field
    unchanged, but for a field ``Y`` MATLAB adds to the end buffer, the rotated y (measured)"""
    T = rotation_matrix(size_zyx, angle, z_aniso)
    out = []
    for tr in tracks:
        t = dict(tr)
        srcs = [(t, "")]
        for b in ("startBuffer", "endBuffer"):
            if isinstance(t.get(b), dict):
                t[b] = dict(t[b])
                srcs.append((t[b], b))
        for d, b in srcs:
            x = np.asarray(d["x"], dtype=np.float64)
            y = np.asarray(d["y"], dtype=np.float64)
            z = np.asarray(d["z"], dtype=np.float64)
            xr, yr, zr = _forward(T, x, y, z)
            d["x"], d["z"] = xr, zr
            if b == "endBuffer":
                d["Y"] = yr
        out.append(t)
    return out


# -- the function ---------------------------------------------------------------------------------


def process_movie(
    tracks_final,
    frame_info,
    read_frame,
    read_mask,
    *,
    framerate: float,
    movie_length: int,
    size_yxz,
    n_slices: int,
    buffer=(3, 3),
    num_threads: int | None = None,
) -> tuple[list[dict], dict]:
    """one movie's track processing in memory: the tracker's ``tracksFinal`` and the
    detections' ``frameInfo`` (sequences of dicts), ``read_frame(ch, k)`` the DS frame (z, y, x)
    of channel ch (0-based) and frame k (1-based) as double, ``read_mask(k)`` the detection mask
    of frame k. ``size_yxz`` the DS frames' MATLAB size, ``n_slices`` their slice count. Returns
    the processed tracks in MATLAB's order (dicts of ``TRACK_FIELDS``) and the cohorts'
    thresholds (``statistics``)."""
    frame_info = list(frame_info)
    tracks = convert_tracks(tracks_final, frame_info, framerate)
    sigma = _mat2(frame_info[0]["s"])
    w = border_margin(sigma)
    tracks = [t for t in tracks if inside_volume(t, w, size_yxz)]
    for t in tracks:
        t["visibility"] = visibility(t, movie_length, buffer)
    init_buffers(tracks, framerate, buffer)
    read_out(
        tracks,
        read_frame,
        read_mask,
        sigma,
        movie_length,
        framerate,
        buffer,
        num_threads=num_threads,
    )
    tracks = sort_by_lifetime(tracks)
    for t in tracks:
        reverse_z(t, n_slices)
    info = statistics(tracks, framerate, movie_length)
    for t in tracks:
        t["MotionAnalysis"] = motion_analysis(t)
    return tracks, info


def matlab_track(tr: dict) -> dict:
    """a processed track as MATLAB's struct element holds it: the fields of ``TRACK_FIELDS``
    in that order with MATLAB's classes (doubles, logical gapVect, hasMultiplePeaks and
    significance, cells of gap indices, structs of buffers and MotionAnalysis, ``[]``)"""
    out = {}
    for k in TRACK_FIELDS:
        v = tr[k]
        if k in ("nSeg", "visibility", "lifetime_s", "start", "end", "catIdx"):
            out[k] = float(v)
        elif k == "hasMultiplePeaks":
            out[k] = np.bool_(v)
        elif k in ("gapVect", "significantMaster", "significantVsBackground", "significantSlave"):
            out[k] = np.asarray(v, dtype=bool)
        elif k == "gapIdx":
            out[k] = (
                [np.asarray(r, dtype=np.float64).reshape(1, -1) for r in v]
                if isinstance(v, list)
                else v
            )
        elif isinstance(v, dict):
            out[k] = {kk: (np.asarray(vv) if not np.isscalar(vv) else vv) for kk, vv in v.items()}
        else:
            out[k] = np.asarray(v, dtype=np.float64)
    return out


def write_processed_tracks(path: str, tracks: list[dict], proc_flag=(True, True)) -> None:
    """``ProcessedTracks.mat`` as MATLAB saves it (v7.3): ``tracks`` (1 x N struct array) and
    ``processingInfo.procFlag`` (logical 1 x 2)"""
    from petakit5d.io.matfile import StructArray, write_mat

    rows = [matlab_track(t) for t in tracks]
    st = (
        StructArray((1, len(rows)), TRACK_FIELDS, rows)
        if rows
        else StructArray.empty(TRACK_FIELDS, (1, 0))
    )
    info = {"procFlag": np.array([[bool(proc_flag[0]), bool(proc_flag[1])]])}
    write_mat(path, {"tracks": st, "processingInfo": info}, version="7.3")


def _paths(cell) -> list[str]:
    """a MATLAB cell (lists, object arrays, char, nested) of file names as a flat list, in
    MATLAB's (column-major) order"""
    if isinstance(cell, str):
        return [cell]
    if isinstance(cell, np.ndarray):
        if cell.dtype.kind == "U":
            return (
                ["".join(cell.ravel())]
                if cell.ndim and cell.shape[0] == 1
                else [str(v) for v in cell.ravel()]
            )
        items = list(cell.ravel(order="F"))
    else:
        items = list(cell)
    out = []
    for v in items:
        out.extend(_paths(v))
    return out


def _channel_paths(v) -> list[list[str]]:
    """``framePathsDS`` (a 1 x nCh cell of cells of file names) as one list of names per channel"""
    if isinstance(v, np.ndarray) and v.dtype == object:
        items = list(v.ravel(order="F"))
    elif isinstance(v, (list, tuple)):
        items = list(v)
    else:
        items = [v]
    if items and all(isinstance(x, str) for x in items):
        return [list(items)]
    return [_paths(x) for x in items]


def _movies(data) -> list[dict]:
    if isinstance(data, dict):
        return [data]
    return [dict(d) for d in data]


def _default_results(movies) -> list[str]:
    """MATLAB's default ``ResultsPath``, each movie's ``[source 'Analysis' filesep]`` with the
    platform's ``filesep`` (``core.paths``; a movie's ``source`` ends with one)"""
    from petakit5d.core import paths

    return [str(m["source"]) + "Analysis" + paths.filesep for m in movies]


_DEFAULTS_ONLY = ("fit_mixtures", "preprocess", "postprocess", "window_size", "buffer_all")


@register(
    "XR_runTrackProcessing3D",
    Param(
        "data",
        "@isstruct",
        "Movies (the structure XR_loadConditionData3D returns).",
        convert="same",
    ),
    Param(
        "Buffer",
        "@(x) numel(x)==2",
        "Buffer frames read out before and after each track.",
        convert="same",
    ),
    Param("Overwrite", "@islogical", "Process movies whose result exists again."),
    Param("FitMixtures", "@islogical", "Fit mixtures (not ported: false only)."),
    Param("TrackerOutput", "@ischar", "Tracker result file in the results folder."),
    Param("FileName", "@ischar", "Result file name."),
    Param("DetectionFile", "@ischar", "Detection file name in the results folder."),
    Param(
        "Frames",
        "@(x) numel(unique(diff(x)))==1",
        "Frames (0-based; default all; all only).",
        matlab_default="arrayfun(@(x) 1:x.movieLength, data, 'unif', 0)",
        convert="index1",
    ),
    Param("Preprocess", "@islogical", "Preprocessing (true only)."),
    Param("Postprocess", "@islogical", "Postprocessing (true only)."),
    Param("WindowSize", None, "Fit window size (default only).", convert="same"),
    Param("BufferAll", "@islogical", "Buffers for every track (false only)."),
    Param(
        "ResultsPath",
        "@iscell",
        "Results folder of each movie (default <source>/Analysis/).",
        matlab_default="arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0)",
    ),
    Param("ParallelMethod", "@ischar", "Accepted; movies are processed in this process."),
    Param("WriteAmira", "@islogical", "Write the Amira files (<results>/Amira/Tracking)."),
    task="C04",
)
def run_track_processing(
    data,
    *,
    buffer=(3, 3),
    overwrite: bool = False,
    fit_mixtures: bool = False,
    tracker_output: str = "trackedFeatures.mat",
    file_name: str = "ProcessedTracks.mat",
    detection_file: str = "Detection3D.mat",
    frames=None,
    preprocess: bool = True,
    postprocess: bool = True,
    window_size=None,
    buffer_all: bool = False,
    results_path=None,
    parallel_method: str = "lock_distributed",
    write_amira: bool = True,
) -> None:
    """``XR_runTrackProcessing3D``: each movie's tracks (``<results>/<tracker_output>``) processed
    with its detections, frames and masks into ``<results>/<file_name>`` and, with
    ``write_amira``, ``<results>/Amira/Tracking/tracking_t_<k>.am`` (the tracks of categories 1
    to 4). The defaults of the options that change the processing are the only ones ported."""
    from petakit5d.detection import amira
    from petakit5d.detection.run import dmask_name
    from petakit5d.io import read_tiff
    from petakit5d.io.matfile import read_mat

    del parallel_method
    for nm, v, want in zip(
        _DEFAULTS_ONLY,
        (fit_mixtures, preprocess, postprocess, window_size, buffer_all),
        (False, True, True, None, False),
        strict=True,
    ):
        empty = want is None and (v is None or np.asarray(v).size == 0)
        if not (empty or (want is not None and bool(v) == want)):
            raise NotImplementedError(
                f"XR_runTrackProcessing3D: {nm}={v!r} is not ported (only the default)"
            )
    if [int(v) for v in np.asarray(buffer).ravel()] != [3, 3]:
        raise NotImplementedError("XR_runTrackProcessing3D: Buffer other than [3 3] is not ported")
    movies = _movies(data)
    if results_path is None:
        results_path = _default_results(movies)
    elif isinstance(results_path, str):
        results_path = [results_path]
    for k, m in enumerate(movies):
        rp = str(results_path[k])
        out = os.path.join(rp, str(file_name))
        if os.path.exists(out) and not overwrite:
            continue
        n = int(np.asarray(m["movieLength"]).ravel()[0])
        if frames is not None and [int(v) for v in np.asarray(frames).ravel()] != list(range(n)):
            raise NotImplementedError(
                "XR_runTrackProcessing3D: Frames other than all frames is not ported"
            )
        framerate = float(np.asarray(m["framerate"]).ravel()[0])
        tf = read_mat(os.path.join(rp, str(tracker_output)), "tracksFinal")["tracksFinal"]
        tf = [tf] if isinstance(tf, dict) else list(tf)
        fi = read_mat(os.path.join(rp, str(detection_file)), "frameInfo")["frameInfo"]
        fi = [fi] if isinstance(fi, dict) else list(fi)
        ds = _channel_paths(m["framePathsDS"])
        first = read_tiff(ds[0][0])
        nz, ny, nx = first.shape

        def read_frame(ch, f, ds=ds):
            return np.asarray(read_tiff(ds[ch][f - 1]), dtype=np.float64)

        def read_mask(f, rp=rp, n=n):
            return read_tiff(os.path.join(rp, "Masks", dmask_name(f, n)))

        tracks, _ = process_movie(
            tf,
            fi,
            read_frame,
            read_mask,
            framerate=framerate,
            movie_length=n,
            size_yxz=(ny, nx, nz),
            n_slices=nz,
        )
        tmp = out + ".tmp_" + str(os.getpid())
        write_processed_tracks(tmp, tracks, (preprocess, postprocess))
        os.replace(tmp, out)
        if write_amira:
            keep = [t for t in tracks if int(t["catIdx"]) in (1, 2, 3, 4)]
            za = float(np.asarray(m["zAniso"]).ravel()[0])
            amira.write_tracks(
                os.path.join(rp, "Amira", "Tracking", "tracking"), keep, n, za, nz, flip_z=True
            )


def _rotate_v73(src: str, dst: str, T) -> None:
    """``RotatedTracks.mat`` from a v7.3 ``ProcessedTracks.mat`` written by MATLAB or the port
    (a struct array of references): the file copied and only what the rotation changes written
    into it (x and z of the tracks and of the buffers, the end buffers' new field ``Y``,
    ``processingInfo`` removed). The same values as ``rotated_tracks`` gives, read back, and as
    MATLAB's ``RotatedTracks.mat`` from the same file; reading and writing every object of the
    struct array (about 200,000 on the demo) takes 12 times as long."""
    import shutil

    import h5py

    from petakit5d.io import matfile as _mf

    shutil.copyfile(src, dst)
    vt = h5py.vlen_dtype(np.dtype("S1"))
    with h5py.File(dst, "r+", locking=False) as f:
        if "processingInfo" in f:
            del f["processingInfo"]
        g = f["tracks"]
        w = _mf._H5Writer(f, True)

        def rotate(obj, add_y: bool) -> None:
            dx, dy, dz = obj["x"], obj["y"], obj["z"]
            x, y, z = (np.asarray(d[()], dtype=np.float64).T for d in (dx, dy, dz))
            xr, yr, zr = _forward(T, x, y, z)
            dx[...] = np.ascontiguousarray(np.asarray(xr).T)
            dz[...] = np.ascontiguousarray(np.asarray(zr).T)
            if add_y:
                names = [*_mf._H5Reader._fields(obj), "Y"]
                w.write(obj, "Y", _mf._to_node(np.asarray(yr), "Y"))
                fields = np.empty(len(names), dtype=vt)
                for i, k in enumerate(names):
                    fields[i] = np.frombuffer(k.encode("ascii"), dtype="S1")
                del obj.attrs["MATLAB_fields"]
                obj.attrs.create("MATLAB_fields", fields, dtype=vt)

        class _Track:
            """one element's x, y, z through the struct array's references"""

            def __init__(self, i):
                self.i = i

            def __getitem__(self, k):
                return f[refs[k][self.i]]

        refs = {k: g[k][()].ravel() for k in ("x", "y", "z", "startBuffer", "endBuffer")}
        for i in range(refs["x"].size):
            rotate(_Track(i), False)
            for b in ("startBuffer", "endBuffer"):
                obj = f[refs[b][i]]
                if isinstance(obj, h5py.Group):  # a struct; [] is a dataset
                    rotate(obj, b == "endBuffer")


def rotate_tracks(data, *, overwrite: bool = False, crop: bool = False, results_path=None) -> None:
    """``rotateTracks3D(data, 'Crop', false, ...)``: ``<results>/RotatedTracks.mat`` (v7.3,
    ``tracks``) from ``<results>/ProcessedTracks.mat``, the geometry from the first DS frame of
    channel 1 and the movie's ``angle`` and ``zAniso`` (stage-scan, not reversed: measured). A
    v7.3 file of two or more tracks is rotated in a copy (``_rotate_v73``), any other file
    through ``rotated_tracks``."""
    from petakit5d.io import read_tiff
    from petakit5d.io.matfile import StructArray, mat_version, read_mat, write_mat

    if bool(np.asarray(crop).ravel()[0]):
        raise NotImplementedError("rotateTracks3D: 'Crop' true is not ported")
    movies = _movies(data)
    if results_path is None:
        results_path = _default_results(movies)
    elif isinstance(results_path, str):
        results_path = [results_path]
    for k, m in enumerate(movies):
        if bool(np.asarray(m.get("objectiveScan", False)).ravel()[0]) or bool(
            np.asarray(m.get("reversed", False)).ravel()[0]
        ):
            raise NotImplementedError(
                "rotateTracks3D: objective scans and reversed scans are not ported"
            )
        rp = str(results_path[k])
        out = os.path.join(rp, "RotatedTracks.mat")
        if os.path.exists(out) and not overwrite:
            continue
        src = os.path.join(rp, "ProcessedTracks.mat")
        ds = _channel_paths(m["framePathsDS"])[0]
        size = read_tiff(ds[0]).shape
        angle = float(np.asarray(m["angle"]).ravel()[0])
        za = float(np.asarray(m["zAniso"]).ravel()[0])
        tmp = out + ".tmp_" + str(os.getpid())
        if mat_version(src) == "7.3" and _has_tracks_struct(src):
            _rotate_v73(src, tmp, rotation_matrix(size, angle, za))
        else:
            pt = read_mat(src, "tracks")["tracks"]
            pt = [pt] if isinstance(pt, dict) else list(pt)
            rows = rotated_tracks(pt, size, angle, za)
            fields = tuple(pt[0].keys()) if pt else TRACK_FIELDS
            st = (
                StructArray((1, len(rows)), fields, rows)
                if rows
                else StructArray.empty(fields, (1, 0))
            )
            write_mat(tmp, {"tracks": st}, version="7.3")
        os.replace(tmp, out)


def _has_tracks_struct(path: str) -> bool:
    """whether a v7.3 file's 'tracks' is a struct array of two or more tracks, its fields
    arrays of references (the in-place rotation's case); one track (a scalar struct, its fields
    stored in place) and none go the plain way"""
    import h5py

    with h5py.File(path, "r", locking=False) as f:
        g = f.get("tracks")
        if not isinstance(g, h5py.Group):
            return False
        keys = ("x", "y", "z", "startBuffer", "endBuffer")
        return all(
            isinstance(g.get(k), h5py.Dataset) and h5py.check_dtype(ref=g[k].dtype) is not None
            for k in keys
        )
