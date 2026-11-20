"""Tracking (``petakit5d.detection.tracking``, task C03) against MATLAB's tracker on synthetic
detections: MATLAB's ``trackCloseGapsKalmanSparse`` and ``XR_runTracking3D`` with the demos'
settings, its mex functions (``createDistanceMatrix``, ``KDTreeBallQuery``, ``lap``) and
``loadTrackSettings``, compared bit for bit; and PointDetectionTest's tracking check in intent.

The detections are synthetic and bit-stable (``petakit5d.testing.synth.CounterRNG``): spots in
Brownian motion with births, deaths and missing frames, a dense field of them, spots that merge
and split (their amplitudes adding up), the 3-frame movie of PointDetectionTest, and frames with
no or one detection and repeated positions.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Work
from petakit5d.testing.compare import MArray, MStruct

TRACK = "petakit5d.detection.tracking"
SETTINGS = "loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2)"


# ---------------------------------------------------------------------------------------------
# synthetic detections: per frame (x, y, z) positions, (dx, dy, dz) errors, amplitude and error
# ---------------------------------------------------------------------------------------------


def _frame(P, E, A, dA) -> dict:
    P = np.asarray(P, float).reshape(-1, 3)
    E = np.asarray(E, float).reshape(-1, 3)
    A = np.asarray(A, float).reshape(-1)
    dA = np.asarray(dA, float).reshape(-1)
    if P.shape[0] == 0:  # MATLAB's [] for a frame without detections
        e = np.zeros((0, 0))
        return {"xCoord": e, "yCoord": e.copy(), "zCoord": e.copy(), "amp": e.copy()}
    return {
        "xCoord": np.c_[P[:, 0], E[:, 0]],
        "yCoord": np.c_[P[:, 1], E[:, 1]],
        "zCoord": np.c_[P[:, 2], E[:, 2]],
        "amp": np.c_[A, dA],
    }


def _spots_movie(
    seed: str, n: int, T: int, box, step: float, gap_p: float, life: bool
) -> list[dict]:
    """``n`` spots in Brownian motion in ``box`` (x, y, z), steps of std ``step``; each spot lives
    a random part of the movie (``life``) and misses a frame with probability ``gap_p``"""
    from petakit5d.testing.synth import CounterRNG

    rng = CounterRNG(2026, f"tracking/{seed}")
    box = np.asarray(box, float)
    p0 = rng.uniform((n, 3)) * box
    steps = rng.normal((T, n, 3)) * step
    pos = p0[None] + np.cumsum(steps, axis=0)
    if life:
        u = rng.uniform((n, 2))
        t0 = np.floor(u[:, 0] * (T // 2)).astype(int)
        t1 = np.minimum(T - 1, t0 + 2 + np.floor(u[:, 1] * T).astype(int))
    else:
        t0 = np.zeros(n, int)
        t1 = np.full(n, T - 1)
    miss = rng.uniform((T, n)) < gap_p
    amp0 = 60.0 + 90.0 * rng.uniform(n)
    ampn = rng.normal((T, n)) * 4.0
    err = 0.03 + 0.15 * rng.uniform((T, n, 3))
    frames = []
    for t in range(T):
        on = (t >= t0) & (t <= t1) & ~miss[t]
        k = np.nonzero(on)[0]
        frames.append(_frame(pos[t, k], err[t, k], amp0[k] + ampn[t, k], 1.0 + 0.5 * err[t, k, 0]))
    return frames


def movie_brownian() -> list[dict]:
    return _spots_movie("brownian", 40, 25, (40.0, 40.0, 20.0), 0.6, 0.06, True)


def movie_dense() -> list[dict]:
    return _spots_movie("dense", 300, 20, (50.0, 50.0, 25.0), 0.8, 0.04, True)


def movie_merge_split() -> list[dict]:
    """pairs of spots that come together and merge, and spots that split in two (amplitudes add
    up), among spots that only move; 24 frames"""
    from petakit5d.testing.synth import CounterRNG

    rng = CounterRNG(2026, "tracking/mergesplit")
    T = 24
    frames_p = [[] for _ in range(T)]
    frames_a = [[] for _ in range(T)]

    def add(t, p, a):
        frames_p[t].append(p)
        frames_a[t].append(a)

    jit = rng.normal((T, 64, 3)) * 0.12
    k = 0
    for e in range(6):
        c = np.array([12.0 + 22.0 * (e % 3), 12.0 + 24.0 * (e // 3), 8.0])
        a1, a2 = 80.0 + 10.0 * e, 60.0 + 7.0 * e
        tm = 8 + e
        for t in range(T):
            if e < 3:  # merge at tm: two spots 1.6 apart come together
                if t < tm:
                    sep = 0.8 * max(0.25, (tm - t) / 4)
                    add(t, c + [sep, 0, 0] + jit[t, k], a1)
                    add(t, c - [sep, 0, 0] + jit[t, k + 1], a2)
                else:
                    add(t, c + jit[t, k], a1 + a2)
            else:  # split at tm: one spot becomes two that move apart
                if t < tm:
                    add(t, c + jit[t, k], a1 + a2)
                else:
                    sep = 0.8 * min(2.0, 0.25 + (t - tm) / 4)
                    add(t, c + [0, sep, 0] + jit[t, k], a1)
                    add(t, c - [0, sep, 0] + jit[t, k + 1], a2)
        k += 2
    for s in range(10):
        c = np.array([6.0 + 6.5 * s, 62.0, 4.0 + s])
        for t in range(T):
            if not (s == 3 and t in (10, 11)) and not (s == 7 and t == 15):
                add(t, c + [0.3 * t, 0, 0] + jit[t, 12 + s], 70.0 + 3.0 * s)
    frames = []
    for t in range(T):
        P = np.array(frames_p[t])
        frames.append(_frame(P, np.full(P.shape, 0.05), np.array(frames_a[t]), np.ones(len(P))))
    return frames


def movie_pdt(velocity=(0.6, 1.0, 0.4), z_aniso: float = 1.5) -> list[dict]:
    """PointDetectionTest's movie as the tracker sees its detections: 5 spots ([y x z] at
    frame 1), 3 frames, moving by ``velocity`` [y x z] voxels per frame; z scaled by zAniso"""
    from petakit5d.testing.synth import CounterRNG

    spots = np.array([[16.3, 20.7, 18.2], [16, 48, 32], [48, 24, 40], [40, 40, 20], [32, 32, 48]])
    rng = CounterRNG(2026, "tracking/pdt")
    noise = rng.normal((3, 5, 3)) * 0.03
    frames = []
    for t in range(3):
        yxz = spots + t * np.asarray(velocity) + noise[t]
        P = np.c_[yxz[:, 1], yxz[:, 0], yxz[:, 2] * z_aniso]
        E = np.tile([0.04, 0.04, 0.06 * z_aniso], (5, 1))
        frames.append(_frame(P, E, np.full(5, 200.0) + noise[t, :, 0] * 50, np.full(5, 3.0)))
    return frames


def movie_edges() -> list[dict]:
    """frames with no detection and with one, repeated positions (a zero-cost link)"""
    e = np.full((2, 3), 0.05)
    return [
        _frame([[1, 1, 1], [10, 10, 10]], e, [100, 90], [1, 1]),
        _frame([[1.2, 1, 1]], e[:1], [100], [1]),
        _frame(np.zeros((0, 3)), np.zeros((0, 3)), [], []),
        _frame([[1.5, 1.1, 1], [10, 10.2, 10]], e, [100, 90], [1, 1]),
        _frame([[1.5, 1.1, 1], [10, 10.2, 10]], e, [100, 90], [1, 1]),
        _frame(
            [[1.6, 1.0, 1], [10.1, 10.1, 10], [20, 20, 20]],
            np.full((3, 3), 0.05),
            [100, 90, 80],
            [1, 1, 1],
        ),
    ]


def movie_radius() -> list[dict]:
    """the search radius rules, one scenario per 300 units of x, 11 frames: a link at exactly the
    radius; the nearest neighbour of the 4th frame back; half the neighbour distance; the cap of
    the std multiple; a gap of 2 frames beyond the largest radius of 1 frame; the std of a short
    segment from its mean step"""
    T = 11
    pts: list[list] = [[] for _ in range(T)]

    def put(t, x, y=0.0, z=0.0):
        pts[t].append((x, y, z))

    # S1: frame 0, then a resting spot exactly 6 (the largest radius) away from frame 1 on
    for t in range(6):
        put(t, 300.0 if t == 0 else 306.0)
    # S2: a slow spot with a neighbour 2.5 away in frames 0-3 only; 4.5 away at frame 7
    for t in range(7):
        put(t, 600.0 + 0.3 * t)
    for t in range(4):
        put(t, 600.0 + 0.3 * t, 2.5)
    for t in range(7, 11):
        put(t, 600.0 + 0.3 * 6 + 4.5 + 0.3 * (t - 7))
    # S3: a slow spot with a neighbour 10 away throughout; 4.2 away at frame 7
    for t in range(11):
        put(t, 900.0 + 0.3 * min(t, 6) + (4.2 + 0.3 * (t - 7) if t >= 7 else 0.0), 0.0)
        put(t, 900.0 + 0.3 * t, 10.0)
    # S4: a spot that hardly moves (tiny noise variance); 4.0 away at frame 7
    for t in range(11):
        put(t, 1200.0 + 0.001 * min(t, 6) + (4.0 if t >= 7 else 0.0))
    # S5: a spot in large steps (frames 0-6), then a short segment 7 away at frames 8-9
    for t in range(7):
        put(t, 1500.0 + 2.0 * (t % 2))
    for t in (8, 9):
        put(t, 1500.0 + 7.0)
    # S6: a short segment (frames 0-2), a short one 3.4 on (frames 3-5), a resting neighbour
    for t in range(3):
        put(t, 1800.0 + 0.5 * t)
    for t in range(3, 6):
        put(t, 1800.0 + 1.0 + 3.4 + (t - 3))
    for t in range(6):
        put(t, 1800.0 + 2.7, 1.8)
    frames = []
    for t in range(T):
        P = np.array(pts[t])
        frames.append(_frame(P, np.full(P.shape, 0.05), np.full(len(P), 100.0), np.ones(len(P))))
    return frames


GAP_RULES_D = 4.0


def movie_gap_rules() -> list[dict]:
    """the rules of merging and gap closing, one scenario at each of x 100, 300, 600 and 900, 14
    frames, and 100 spots in Brownian motion of different step sizes far away (the 10th
    percentile of their noise stds is the std of 1-frame segments): a merge that the means of 5
    frames of amplitude allow (4 would not); a segment that can merge into two, where rho^-2
    (rho < 1) decides; a gap to a segment that does not move (std 0: no search region); a gap
    between two 1-frame segments within the radius of that percentile. (MATLAB loses a spot at
    x = 0 from its tracks, so none is there.)"""
    from petakit5d.testing.synth import CounterRNG

    T = 14
    pts: list[list] = [[] for _ in range(T)]

    def put(t, x, y, z, a):
        pts[t].append((x, y, z, a))

    # 1: i (amplitude 50) ends 1 away from j at frame 6; j's amplitude 80, then 70 for 4 frames
    # and 350
    for t in range(T):
        put(t, 100.0, 0.02 * (-1) ** t, 0.0, 80.0 if t < 7 else (350.0 if t == 11 else 70.0))
    for t in range(7):
        put(t, 96.0 + 0.5 * t, 0.0, 0.0, 50.0)
    # 2: i (amplitude 60) ends between j1 (1.0 away; amplitude 65, then 100: rho 0.8) and j2
    # (1.1 away; 65, then 150: rho 1.2)
    for t in range(T):
        put(t, 300.0 + 0.01 * (-1) ** t, 1.0, 0.0, 65.0 if t < 7 else 100.0)
        put(t, 300.0 - 0.01 * (-1) ** t, -1.1, 0.0, 65.0 if t < 7 else 150.0)
    for t in range(7):
        put(t, 297.0 + 0.5 * t, 0.0, 0.0, 60.0)
    # 3: a segment of 3 frames, then 2 frames later and 1 away one of 3 frames that does not move
    for t in range(3):
        put(t, 600.0 + 0.5 * t, 0.0, 0.0, 100.0)
    for t in range(4, 7):
        put(t, 602.0, 0.0, 0.0, 100.0)
    # 4: 1-frame detections at frames 5 and 7 (amplitude 200), GAP_RULES_D apart, the second 2
    # away from a resting spot
    for t in range(T):
        put(t, 900.0, 0.0, 0.02 * (-1) ** t, 100.0)
    put(5, 900.0 - GAP_RULES_D, 2.0, 0.0, 200.0)
    put(7, 900.0, 2.0, 0.0, 200.0)
    # the spots far away, steps of std 0.45 to 1.44
    rng = CounterRNG(2026, "tracking/gap_rules")
    k = np.arange(100)
    steps = rng.normal((T, 100, 3)) * (0.45 + 0.01 * k)[None, :, None]
    p0 = np.c_[2000.0 + 60.0 * (k % 10), 60.0 * (k // 10), np.zeros(100)]
    pos = p0[None] + np.cumsum(steps, axis=0)
    for t in range(T):
        for j in range(100):
            put(t, *pos[t, j], 100.0)
    frames = []
    for t in range(T):
        Q = np.array(pts[t])
        frames.append(_frame(Q[:, :3], np.full((len(Q), 3), 0.05), Q[:, 3], np.ones(len(Q))))
    return frames


MOVIES = {
    "radius": movie_radius,
    "gap_rules": movie_gap_rules,
    "brownian": movie_brownian,
    "dense": movie_dense,
    "merge_split": movie_merge_split,
    "pdt": movie_pdt,
    "edges": movie_edges,
}


def _mframes(frames: list[dict]) -> MStruct:
    """MATLAB's frameInfo (1 x T struct array) of the frames"""
    fields = ["xCoord", "yCoord", "zCoord", "amp"]
    items = []
    for f in frames:
        el = {}
        for k in fields:
            a = np.asarray(f[k], float)
            el[k] = MArray("double", a.reshape(0, 0) if a.size == 0 else a)
        items.append(el)
    return MStruct((1, len(frames)), fields, items)


# ---------------------------------------------------------------------------------------------
# Python sides
# ---------------------------------------------------------------------------------------------


def _track(frames):
    from petakit5d.detection.tracking import track_detections

    return track_detections(frames)


def _run_tracking_case(work, frames):
    """the Python side of XR_runTracking3D on a results folder holding Detection3D.mat"""
    import os

    from petakit5d.detection.tracking import run_tracking
    from petakit5d.io.matfile import StructArray, write_mat

    rp = os.path.join(work, "Analysis") + os.sep
    os.makedirs(rp, exist_ok=True)
    fi = StructArray.from_list([dict(f) for f in frames])
    write_mat(os.path.join(rp, "Detection3D.mat"), {"frameInfo": fi})
    data = {"source": work, "movieLength": float(len(frames))}
    run_tracking(data, file_name="trackedFeatures.mat", overwrite=True, results_path=[rp])
    return None


def _ball_query(X, Y, r):
    from petakit5d.detection.tracking import kd_tree_ball_query

    idx, dist = kd_tree_ball_query(X, Y, r)
    return [(i + 1).astype(np.float64) for i in idx], dist


def _lap(cc, marker, extended):
    from petakit5d.detection.tracking import solve_assignment

    C = np.asarray(cc, float)
    r, c = np.nonzero(C != marker)
    x = solve_assignment(r, c, C[r, c], C.shape[0])
    y = np.empty_like(x)
    y[x] = np.arange(x.size)
    return (x + 1).astype(np.int32).reshape(-1, 1), (y + 1).astype(np.int32).reshape(-1, 1)


def _settings(radius=None, gap=None):
    from petakit5d.detection.tracking import track_settings

    kw = {}
    if radius is not None:
        kw["radius"] = radius
    if gap is not None:
        kw["max_gap_length"] = gap
    return track_settings(**kw)


# ---------------------------------------------------------------------------------------------
# inputs of the mex cases
# ---------------------------------------------------------------------------------------------


def _points(seed: str, n: int, scale: float, d: int = 3) -> np.ndarray:
    from petakit5d.testing.synth import CounterRNG

    return CounterRNG(2026, f"tracking/points/{seed}").uniform((n, d)) * scale


def _aug(seed: str, n1: int, n2: int, integer: bool) -> np.ndarray:
    """an augmented linking problem as the tracker builds it (-5: no link)"""
    from petakit5d.testing.synth import CounterRNG

    rng = CounterRNG(2026, f"tracking/lap/{seed}")
    u = rng.uniform((n1, n2))
    C = np.floor(1 + 4 * u) if integer else 36 * u
    C[rng.uniform((n1, n2)) < 0.5] = -5
    bd = 5.0 if integer else 37.8
    A = np.full((n1 + n2, n1 + n2), -5.0)
    A[:n1, :n2] = C
    A[np.arange(n1), n2 + np.arange(n1)] = bd
    A[n1 + np.arange(n2), np.arange(n2)] = bd
    L = C.T.copy()
    L[L != -5] = bd
    A[n1:, n2:] = L
    return A


def _cases():
    cs = []
    for name, mk in MOVIES.items():
        cs.append(
            Case(
                f"track.{name}",
                code="s = "
                + SETTINGS
                + "; [tf, kl] = trackCloseGapsKalmanSparse(in{1}, s.costMatrices, s.gapCloseParam, "
                "s.kalmanFunctions, 3, 0, 0); out = {tf, kl};",
                args=[In("frames")],
                nargout=2,
                python=_track,
                py_args=[mk()],
                inputs={"frames": Value(lambda mk=mk: _mframes(mk()))},
            )
        )
    for name in ("brownian", "edges"):
        mk = MOVIES[name]
        cs.append(
            Case(
                f"run_tracking.{name}",
                code="rt = in{1}; frameInfo = in{2}; mkdir([rt 'Analysis']); "
                "save([rt 'Analysis' filesep 'Detection3D.mat'], 'frameInfo', '-v7.3'); "
                "data = struct('source', rt, 'channels', {{rt}}, 'movieLength', numel(frameInfo)); "
                "XR_runTracking3D(data, " + SETTINGS + ", 'FileName', 'trackedFeatures.mat', "
                "'Overwrite', true, 'ResultsPath', {[rt 'Analysis' filesep]}); out = {};",
                args=[Work(), In("frames")],
                nargout=0,
                python=_run_tracking_case,
                py_args=[Work(), mk()],
                inputs={"frames": Value(lambda mk=mk: _mframes(mk()))},
            )
        )
    X = _points("x", 400, 20.0)
    Y = _points("y", 30, 20.0)
    cs += [
        Case(
            "create_distance_matrix.random",
            matlab="createDistanceMatrix",
            args=[In("a"), In("b")],
            python="petakit5d.detection.tracking:create_distance_matrix",
            inputs={
                "a": Value(lambda: _points("a", 37, 30.0)),
                "b": Value(lambda: _points("b", 23, 30.0)),
            },
        ),
        Case(
            "create_distance_matrix.empty",
            matlab="createDistanceMatrix",
            args=[np.zeros((0, 3)), np.array([[1.0, 2.0, 3.0]])],
            python="petakit5d.detection.tracking:create_distance_matrix",
        ),
        Case(
            "kd_tree_ball_query.random",
            matlab="KDTreeBallQuery",
            args=[In("x"), In("y"), 1.7 * np.ones((30, 1))],
            nargout=2,
            python=_ball_query,
            inputs={"x": Value(lambda: X), "y": Value(lambda: Y)},
        ),
        Case(
            "kd_tree_ball_query.boundary_and_empty",
            matlab="KDTreeBallQuery",
            args=[
                np.array(
                    [
                        [0.0, 0, 0],
                        [1, 0, 0],
                        [1 - 2.0**-52, 0, 0],
                        [1 + 2.0**-52, 0, 0],
                        [0.5, 0.25, 0],
                    ]
                ),
                np.array([[0.0, 0, 0], [9.0, 9, 9]]),
                1.0,
            ],
            nargout=2,
            python=_ball_query,
        ),
        Case(
            "kd_tree_ball_query.self_2d",
            matlab="KDTreeBallQuery",
            args=[In("p"), In("p"), 0.6],
            nargout=2,
            python=_ball_query,
            inputs={"p": Value(lambda: _points("p2", 40, 4.0, 2))},
        ),
    ]
    for name, n1, n2, integer in (
        ("random", 40, 37, False),
        ("ties", 9, 8, True),
        ("ties_wide", 30, 24, True),
    ):
        cs.append(
            Case(
                f"lap.{name}",
                matlab="lap",
                args=[In("cc"), -5.0, 0.0],
                nargout=2,
                python=_lap,
                inputs={
                    "cc": Value(
                        lambda n1=n1, n2=n2, integer=integer, name=name: _aug(name, n1, n2, integer)
                    )
                },
            )
        )
    zero = np.array([[0.0, 1, 5, -5], [-5, 2, -5, 5], [5, -5, 5, -5], [-5, 5, 5, 5]])
    cs.append(Case("lap.zero_cost", matlab="lap", args=[zero, -5.0, 0.0], nargout=2, python=_lap))
    for name, args, kw in (
        ("demo", "'Radius', [3 6], 'MaxGapLength', 2", {"radius": (3, 6), "gap": 2}),
        ("default", "", {}),
        ("radius_2_8_gap_4", "'Radius', [2 8], 'MaxGapLength', 4", {"radius": (2, 8), "gap": 4}),
        ("gap_1", "'MaxGapLength', 1", {"gap": 1}),
    ):
        cs.append(
            Case(
                f"track_settings.{name}",
                code=f"out = {{loadTrackSettings({args})}};",
                python=lambda kw=kw: _settings(**kw),
            )
        )
    return cs


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# tests of our own
# ---------------------------------------------------------------------------------------------


def test_tracking_links_spots_across_frames():
    """PointDetectionTest.trackingLinksSpotsAcrossFrames in intent: each of the 5 spots forms one
    track over the 3 frames, its mean step the velocity ([y x z], z scaled by zAniso) within 0.1"""
    from petakit5d.detection.tracking import track_detections

    velocity = np.array([0.6, 1.0, 0.4])
    z_aniso = 1.5
    tracks, _ = track_detections(movie_pdt(velocity, z_aniso))
    assert len(tracks) == 5
    for t in tracks:
        coords = np.asarray(t["tracksCoordAmpCG"]).reshape(-1, 8)
        assert coords.shape == (3, 8)
        step = np.mean(np.diff(coords[:, [1, 0, 2]], axis=0), axis=0)
        np.testing.assert_allclose(step, velocity * [1, 1, z_aniso], atol=0.1)


def test_frames_without_detections_and_zero_cost_links():
    from petakit5d.detection.tracking import track_detections

    tracks, kal = track_detections(movie_edges())
    # both spots of frame 1 are joined across the empty frame 3 to frames 4 to 6, their repeated
    # positions of frames 4 and 5 linked at cost 0; the spot of frame 6 alone
    rows = [np.asarray(t["tracksFeatIndxCG"]) for t in tracks]
    assert [r.shape for r in rows] == [(1, 6), (1, 6), (1, 1)]
    assert rows[0].tolist() == [[1, 1, 0, 1, 1, 1]] and rows[1].tolist() == [[2, 0, 0, 2, 2, 2]]
    assert kal[2]["stateVec"].shape == (0, 6) and kal[2]["stateCov"].shape == (6, 6, 0)
    assert kal[1]["stateCov"].shape == (6, 6)  # one detection: a 6 x 6 matrix, as MATLAB's


def test_saved_file_reads_back(tmp_path):
    """trackedFeatures.mat: the variables MATLAB saves, the tracks a column struct array with
    double indices, one filter state per frame"""
    from petakit5d.detection.tracking import (
        track_detections,
        track_settings,
        write_tracked_features,
    )
    from petakit5d.io.matfile import read_mat

    frames = movie_dense()
    tracks, kal = track_detections(frames)
    fn = str(tmp_path / "trackedFeatures.mat")
    write_tracked_features(fn, track_settings((3, 6), 2), tracks, kal)
    m = read_mat(fn)
    assert list(m) == [
        "costMatrices",
        "gapCloseParam",
        "kalmanFunctions",
        "tracksFinal",
        "kalmanInfoLink",
    ]
    assert m["tracksFinal"].shape == (len(tracks), 1)
    assert m["kalmanInfoLink"].shape == (1, len(frames))
    for a, b in zip(m["tracksFinal"], tracks, strict=True):
        assert np.asarray(a["tracksFeatIndxCG"]).dtype == np.float64
        assert np.array_equal(a["seqOfEvents"], b["seqOfEvents"], equal_nan=True)


def test_nearest_neighbour_of_a_lone_detection():
    """MATLAB gives the only detection of a frame the nearest-neighbour distance 1000 (measured)"""
    from petakit5d.detection.tracking.tracks import frames_from_frame_info

    e = np.full((2, 3), 0.05)
    f = frames_from_frame_info(
        [_frame([[1, 2, 3]], e[:1], [100], [1]), _frame([[1, 2, 3], [1, 2, 7]], e, [1, 2], [1, 1])]
    )
    assert f.nn.tolist() == [1000.0, 4.0, 4.0]


def test_events_starts_before_ends_and_merges_at_their_last_frame():
    """seqOfEvents of a compound track in MATLAB's order (as on the demo's tracks): the events
    by frame, a merge's end counted at the merging segment's last frame (so it comes before the
    events of the frame it is reported at), at one frame the starts before the ends"""
    from types import SimpleNamespace

    from petakit5d.detection.tracking.linking import Frames
    from petakit5d.detection.tracking.tracks import _compound

    T, n = 10, 4
    tracks = np.zeros((n, T), dtype=np.int64)
    tracks[0, 0:10] = 1  # row 1
    tracks[1, 3:6] = 2  # splits from row 1 at frame 4 (1-based), merges into row 3 at frame 7
    tracks[2, 4:9] = 3
    tracks[3, 6] = 4  # splits from row 1 at frame 7, ends there
    off = np.arange(0, n * T + 1, n)
    frames = Frames(
        np.zeros((n * T, 3)), np.zeros((n * T, 3)), np.ones((n * T, 2)), np.zeros(n * T), off
    )
    gp = SimpleNamespace(
        merge_seg=np.array([2]),
        merge_time=np.array([6]),
        split_seg=np.array([0, 0]),
        split_time=np.array([3, 6]),
    )
    x = np.full(2 * n + 3, 10**6)
    x[1] = n  # segment 2 ends in the merge column
    x[n] = 1  # the first split starts segment 2
    x[n + 1] = 3  # the second split starts segment 4
    (ct,) = _compound(tracks, frames, x, gp)
    assert ct["tracksFeatIndxCG"].shape == (4, 10)
    assert np.array_equal(
        ct["seqOfEvents"],
        np.array(
            [
                [1, 1, 1, np.nan],
                [4, 1, 2, 1],
                [5, 1, 3, np.nan],
                [7, 2, 2, 3],
                [7, 1, 4, 1],
                [7, 2, 4, np.nan],
                [9, 2, 3, np.nan],
                [10, 2, 1, np.nan],
            ]
        ),
        equal_nan=True,
    )


def test_results_are_repeatable():
    from petakit5d.detection.tracking import track_detections

    a, _ = track_detections(movie_brownian())
    track_detections(movie_dense())
    b, _ = track_detections(movie_brownian())
    assert len(a) == len(b)
    for x, y in zip(a, b, strict=True):
        assert np.array_equal(x["tracksFeatIndxCG"], y["tracksFeatIndxCG"])
        assert np.array_equal(x["seqOfEvents"], y["seqOfEvents"], equal_nan=True)


def test_linear_motion_is_not_ported():
    from petakit5d.detection.tracking import track_detections, track_settings

    s = track_settings()
    s["costMatrices"][0]["parameters"]["linearMotion"] = 1.0
    with pytest.raises(NotImplementedError):
        track_detections(movie_pdt(), s)
