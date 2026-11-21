"""The lowSNR detector's bad-point filter (``filter_bad_points_in_detection``, task C02c): its
second valley search against MATLAB R2024b through the oracle harness, on synthetic points whose
residual sigmas form a designed histogram (a background peak, a valley at 0.2 of it, a small bump,
the main peak), and the search's measured rules on designed histograms (no MATLAB).

The rules were measured with designed histograms in MATLAB's filter as a black box (mocked
histogram and peak list, every point given two neighbours so that the keep mask shows the
removal threshold; no source read): docs/parity/detection.md, the bad-point filter.
"""

import numpy as np
import pytest

from petakit5d.detection import badpoints as BP
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, Value

# the designed histogram's counts per bin (no data of any frame): the sigmas are placed evenly
# inside the bins of the binning the filter makes of them (iterated with the port's binning,
# which is MATLAB's on every recorded call)
HEAD = [90, 26, 23, 17, 19, 25, 17, 20, 40, 46, 55, 70, 80, 84, 74, 66, 58, 51, 45, 39, 34, 30,
        27, 24, 21, 19, 17, 15, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 3, 2, 1, 2, 1, 0, 1, 0, 1,
        0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1]  # fmt: skip
TAIL = 14


def _place(counts, e):
    sr = []
    for b, n in enumerate(counts):
        if n:
            a, w = e[b], e[b + 1] - e[b]
            sr.append(a + (np.arange(n) + 0.5) / n * w)
    sr.append(e[len(counts)] + 0.5 + np.arange(TAIL) * 0.37)
    return np.concatenate(sr)


def designed_sigma_r(head=HEAD) -> np.ndarray:
    e = 2.81 + np.arange(61) * 0.4
    for _ in range(8):
        hd = (list(head) + [0] * 80)[: len(e) - 1]
        sr = _place(hd, e)
        _, _, seq = BP.stage12(sr)
        _, cnt, e2, _ = seq[-1]
        if len(e2) == len(e) and np.allclose(e2, e) and np.array_equal(cnt, hd):
            break
        e = np.asarray(e2, dtype=np.float64)
    return sr


# the design with bin 4 at 19 (0.21 of the background peak): the second search's A262
HEAD19 = HEAD[:3] + [19] + HEAD[4:]


def _inputs(head=HEAD) -> dict:
    sr = designed_sigma_r(head)
    n = sr.size
    box = np.array([120.0, 110.0, 30.0])
    P = synth.CounterRNG(2026, "C02c.badpoints.pos").uniform((n, 3)) * box
    k = np.arange(1, n + 1)
    return {
        "P": P + 1.0,
        "sr": sr.reshape(1, -1),
        "vol": synth.CounterRNG(2026, "C02c.badpoints.vol").uniform((9, 8, 7)) * 10 + 100,
        "A": (50.0 + 10.0 * np.mod(k, 7)).reshape(1, -1),
        "c": (100.0 + np.mod(k, 5)).reshape(1, -1),
    }


def _py_filter(P, sr, thr, bg, vol, A, c):
    return BP.filter_bad_points(P, np.ravel(sr), thr, np.ravel(bg), vol, np.ravel(A), np.ravel(c))


def _args(bg):
    return [In("P"), In("sr"), 3.0, bg, In("vol"), In("A"), In("c")]


def _cases():
    cs = []
    # camera thresholds 3 (B with the second valley's threshold T2 = 4.86 below 2 x 3); 2.4, 2.1
    # and 2 (T2 at or above 2 cam, the first bin's centre 2.795 below 1.5 cam and bin 4 below
    # 0.2 of the background peak: B); 1.85 (1.5 cam below 2.795: the first valley's A)
    vals = {k: Value(lambda v=v: v) for k, v in _inputs().items()}
    for cam in (3.0, 2.4, 2.1, 2.0, 1.85):
        cs.append(
            Case(
                f"second_valley.cam{cam:g}",
                matlab="filter_bad_points_in_detection",
                args=_args([100.0, 100.0, 2.6, cam]),
                python=_py_filter,
                inputs=vals,
            )
        )
    # bin 4 at 0.21 of the peak (camera threshold 2): A262 with min(noise, 1.5 cam), 2.6 with
    # noise 2.6, 3 with noise 3.2 (the noise rule's 2.88 below it)
    vals19 = {k: Value(lambda v=v: v) for k, v in _inputs(HEAD19).items()}
    for nz in (2.6, 3.2):
        cs.append(
            Case(
                f"second_valley.a262.noise{nz:g}",
                matlab="filter_bad_points_in_detection",
                args=_args([100.0, 100.0, nz, 2.0]),
                python=_py_filter,
                inputs=vals19,
            )
        )
    cs += _few_cases()
    return cs


def _few_inputs(sr, seed: str, box=30.0) -> dict:
    """inputs of the filter for the sigmas ``sr`` (NaN allowed) at uniform positions"""
    sr = np.asarray(sr, dtype=np.float64)
    n = sr.size
    k = np.arange(1, n + 1)
    return {
        "P": synth.CounterRNG(2026, "C02c.badpoints.few." + seed).uniform((n, 3)) * box + 1.0,
        "sr": sr.reshape(1, -1),
        "vol": synth.CounterRNG(2026, "C02c.badpoints.vol").uniform((9, 8, 7)) * 10 + 100,
        "A": (50.0 + 10.0 * np.mod(k, 7)).reshape(1, -1),
        "c": (100.0 + np.mod(k, 5)).reshape(1, -1),
    }


def _few_sigmas() -> dict:
    """the few-point cases' sigmas: one and two points (MATLAB's scalar true), 40 within a ratio
    of 2 (the short path), with NaN, 8 beyond it (the long path's early return), and the
    designed histogram's sigmas with NaN among them (the long path)"""
    r = synth.CounterRNG(2026, "C02c.badpoints.few.sr")
    s40 = 3.0 + 2.5 * r.uniform((40,))
    s40n = s40.copy()
    s40n[[3, 11, 25, 31]] = np.nan
    big = designed_sigma_r()
    bign = big.copy()
    bign[[0, 100, 400, 401, 700]] = np.nan
    return {
        "n1": [3.4],
        "n2": [np.nan, 3.3],
        "short40": s40,
        "short40_nan": s40n,
        "early8": 3.0 + 13.0 * (np.arange(8) + 0.5) / 8,
        "long_nan": bign,
    }


def _few_cases():
    cs = []
    for name, sr in _few_sigmas().items():
        vals = {k: Value(lambda v=v: v) for k, v in _few_inputs(sr, name).items()}
        bgs = [[100.0, 100.0, 2.6, 3.0]]
        if name.startswith("short40"):
            bgs.append([100.0, 100.0, 4.4, 3.0])  # the noise rule at min(0.9 noise, thr 4)
        for bg in bgs:
            cs.append(
                Case(
                    f"few.{name}.noise{bg[2]:g}",
                    matlab="filter_bad_points_in_detection",
                    args=[
                        In("P"),
                        In("sr"),
                        4.0 if name.startswith("short") else 3.0,
                        bg,
                        In("vol"),
                        In("A"),
                        In("c"),
                    ],  # fmt: skip
                    python=_py_filter,
                    inputs=vals,
                )
            )
    return cs


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- rules (no MATLAB) ----------------------------------------------------------------------------

E = 2.8 + 0.445 * np.arange(55)  # 54 bins of width 0.445 from 2.8


def _hist(head):
    """a 54-bin histogram: ``head`` from bin 1, then the design's main peak and tail"""
    c = np.zeros(54)
    h = list(head) + HEAD[len(head) :]
    c[:54] = h[:54]
    return c


def test_the_search_starts_on_the_peaks_share():
    """lines 221-222: the three bins from max(bp - 1, 1) above 0.12 of all counts, or the first
    two of them above 0.09 (strict: a share of exactly 0.12 does not start it)"""
    rest = [51.0] * 16 + [64.0]  # 880 in 17 bins
    c = np.array([60.0, 30.0, 30.0, *rest])  # 120 of 1000: 0.12; bins 1-2: 0.09
    assert BP.second_trigger(c, 1) is False
    c[2] += 1  # 121 of 1001
    assert BP.second_trigger(c, 1) is True
    d = np.array([70.0, 21.0, 0.0, *([53.0] * 16 + [61.0])])  # bins 1-2: 91 of 1000
    assert BP.second_trigger(d, 1) is True
    g = np.array([0.0] * 7 + [60.0, 70.0, 70.0] + [40.0] * 20)  # around bin 9: 200 of 1000
    assert BP.second_trigger(g, 9) is True
    assert BP.second_trigger(g, 1) is False


def test_the_levels_are_the_peak_over_4_to_2():
    """lines 226-229: the levels cbp / 4, / 3.5, / 3, / 2.5, / 2 in turn; the first bin after the
    peak and before the next listed peak strictly below the level; T2 = its upper edge (B with it
    for camera threshold 3 when found at one of the first four levels)"""
    for b4, level in ((57, 2), (58, 3), (66, 3), (67, 4), (79, 4), (80, 5), (99, 5)):
        c = _hist([200, 150, 120, b4, b4 + 3, 48])
        info = {}
        sv = BP.second_valley(c, E, 1, 6, 3.0, info=info)
        assert info["second"] == (level, 4), b4
        assert sv == (("B", E[4]) if level <= 4 else None), b4
    c = _hist([200, 150, 120, 100, 103, 48])  # 100 is not below 200 / 2: not found
    info = {}
    assert BP.second_valley(c, E, 1, 6, 3.0, info=info) is None and info["second"] is None
    c = _hist([80, 60, 45, 20, 19, 24])  # 20 is not below 80 / 4: bin 5 is the first
    info = {}
    assert BP.second_valley(c, E, 1, 6, 3.0, info=info) == ("B", E[5])
    assert info["second"] == (1, 5)


def test_the_camera_threshold_decides_the_path():
    """lines 238-252: B with T2 = e(5) = 4.58 when T2 < 2 cam (cam not 5); otherwise, for cam
    below 5 and the first bin's centre 3.0225 below 1.5 cam, B when bin 4 (17) is below 0.2 of
    the peak (90); else the first valley's decision"""
    c = _hist([90, 26, 23, 17, 19, 25])
    want = {2.0: None, 2.01: None, 2.02: ("B", E[4]), 2.1: ("B", E[4]), 2.25: ("B", E[4]),
            2.3: ("B", E[4]), 3.0: ("B", E[4]), 4.99: ("B", E[4]), 5.0: None,
            5.01: ("B", E[4]), 6.0: ("B", E[4])}  # fmt: skip
    for cam, sv in want.items():
        for nz in (1.5, 2.6, 3.5):  # the noise does not enter
            assert BP.second_valley(c, E, 1, 6, cam, nz) == sv, (cam, nz)


def test_the_count_at_the_second_valley_decides_b_or_a262():
    """lines 248 and 258-262: with bin 4 at 0.2 of the peak or above (18 of 90), A262 with
    min(noise, 1.5 cam), when found at one of the first three levels; at level 4 or 5 the first
    valley's decision"""
    c = _hist([90, 26, 23, 18, 19, 25])
    assert BP.second_valley(c, E, 1, 6, 2.1, 2.6) == ("A262", 2.6)
    assert BP.second_valley(c, E, 1, 6, 2.1, 3.5) == ("A262", 1.5 * 2.1)
    assert BP.second_valley(c, E, 1, 6, 5.0, 2.6) is None  # camera threshold 5
    for head, level, sv in (([90, 30, 27, 25, 27, 32], 2, ("A262", 2.6)),
                            ([90, 35, 31, 29, 30, 36], 3, ("A262", 2.6)),
                            ([90, 40, 37, 35, 36, 42], 4, None),
                            ([90, 50, 47, 44, 46, 52], 5, None)):  # fmt: skip
        info = {}
        assert BP.second_valley(_hist(head), E, 1, 6, 2.1, 2.6, info) == sv, level
        assert info["second"] == (level, 4)
    # T2 = e(14) = 8.585 at or above 2 cam for cam 4, bin 13 (21) above 0.2 of 90: A262
    far = _hist([90, 85, 80, 76, 72, 68, 64, 60, 56, 52, 48, 44, 21, 23, 31])
    assert BP.second_valley(far, E, 1, 15, 4.0, 2.6) == ("A262", 2.6)
    assert BP.second_valley(far, E, 1, 15, 6.0, 2.6) == ("B", E[13])  # cam above 5: 8.585 < 12
    # cam above 5 with T2 = e(26) = 13.925 at or above 2 cam: the first valley's decision (a
    # recorded call with cam 6 took neither block)
    farther = _hist([*np.round(np.linspace(90, 24, 24)), 21, 23, 31])
    assert BP.second_trigger(farther, 1)
    assert BP.second_valley(farther, E, 1, 28, 4.5, 2.6) == ("A262", 2.6)
    assert BP.second_valley(farther, E, 1, 28, 5.5, 2.6) is None


def test_the_first_bins_centre_bounds_the_camera_threshold():
    """lines 247 and 257 take the first bin's centre (3.0225), not the background peak's: with
    the peak at bin 3 (centre 3.9125) and cam 2.1, B with T2 = e(7)"""
    c = _hist([10, 40, 90, 26, 23, 17, 19, 25])
    assert BP.second_trigger(c, 3)
    assert BP.second_valley(c, E, 3, 8, 2.1, 2.6) == ("B", E[6])
    assert BP.second_valley(c, E, 3, 8, 2.0, 2.6) is None


def test_the_second_valley_decides_in_the_decision():
    """the decision on the A paths: with the second valley found, B with T2 (also from a valley
    0.32 of the peak, whose first decision is A with thr), or A262"""
    sr = np.linspace(2.8, 20.0, 500)
    for head in ([90, 26, 23, 17, 19, 25], [90, 42, 32, 29, 31, 36]):
        cnt = _hist(head)
        st = BP.stage3(sr, cnt.size, cnt, E, [6.0, 14.0, 39.0])
        path, t, info = BP.decide(st, 3.0, 2.65, 3.0)
        assert (path, t) == ("B", E[4]), head
        assert info["second"][1] == 4
    cnt = _hist([90, 26, 23, 18, 19, 25])
    st = BP.stage3(sr, cnt.size, cnt, E, [6.0, 14.0, 39.0])
    assert BP.decide(st, 3.0, 2.65, 2.1)[:2] == ("A262", 2.65)


def test_the_first_valleys_b_skips_camera_threshold_5():
    """a valley below 0.15 of the peak gives B with T_c below 2 cam, but not with cam 5 (A with
    T_c: two designed histograms)"""
    sr = np.linspace(2.8, 20.0, 500)
    cnt = _hist([90, 26, 23, 10, 19, 25])
    st = BP.stage3(sr, cnt.size, cnt, E, [6.0, 14.0, 39.0])
    assert BP.decide(st, 3.0, 2.65, 4.9)[0] == "B"
    path, tc, _ = BP.decide(st, 3.0, 2.65, 5.0)
    assert path == "A317" and tc == pytest.approx(E[0] + 3.5 * 0.445)


def test_the_a262_removal():
    """A262 (probes s11, s12, every point given k neighbours): counts 2-3 below its threshold T,
    the noise rule (counts 2..U) below min(T_noise, T); U by B's rule"""
    sr = np.linspace(2.9, 4.5, 400)
    for cnt_v, nz, cam, want in ((3, 3.6, 2.5, 3.6), (4, 3.6, 2.5, 0.9 * 3.6),
                                 (5, 3.6, 2.5, 0.9 * 3.6), (2, 3.6, 2.1, 1.5 * 2.1),
                                 (4, 3.6, 2.1, 1.5 * 2.1)):  # fmt: skip
        cnt = np.full(sr.size, cnt_v)
        t = min(nz, 1.5 * cam)
        U = BP.median_bound(sr, cnt, "A262", t)
        keep = BP.removal(sr, cnt, 3.0, nz, cam, "A262", t, U=U)
        assert np.array_equal(keep, ~(sr < want)), (cnt_v, nz, cam)
    # U: B's rule (the candidates' median count 8 or more), not A's
    cnt = np.full(sr.size, 6)
    assert BP.median_bound(sr, cnt, "A262", 4.0) == 8
    assert BP.median_bound(sr, cnt, "A317", 4.0) == 6


def test_b_takes_the_median_when_most_candidates_have_many_neighbours():
    """B's bound (lines 351-357, C02's rule): taken when more than 100 candidates (sr < T) and
    more than half of them above 8 neighbours; U = min(their median, 10)"""
    sr = np.linspace(2.9, 4.5, 400)[1:]
    k = np.flatnonzero(sr < 4.0)  # 274 candidates
    for above, want in ((138, 9), (137, 8)):  # 137 of 274 above 8: half, not more
        cnt = np.full(sr.size, 3)
        cnt[k[: k.size - above]] = 8
        cnt[k[k.size - above :]] = 9
        assert BP.median_bound(sr, cnt, "B", 4.0) == want, above
    # 100 candidates (all with 12 neighbours) do not take it, 101 do
    assert BP.median_bound(sr[:100], np.full(100, 12), "B", 4.0) == 8
    assert BP.median_bound(sr[:101], np.full(101, 12), "B", 4.0) == 10


def test_the_designed_input_reaches_the_second_valley():
    """the oracle cases' inputs: the filter's histogram is the design; the second search changes
    the decision with camera threshold 3, 2.1 and 2 (B at T2) but not with 1.85, and gives A262
    with bin 4 at 19"""
    I = _inputs()  # noqa: E741
    out = {}
    for cam in (3.0, 2.1, 2.0, 1.85):
        out[cam] = BP.filter_bad_points_info(I["P"], I["sr"].ravel(), 3.0, [100, 100, 2.6, cam])
    assert out[3.0][1]["path"] == "B" and out[3.0][1]["decide"]["second"] == (1, 4)
    for cam in (2.1, 2.0):
        assert out[cam][1]["path"] == "B" and out[cam][1]["tc"] == out[3.0][1]["tc"]
    assert out[1.85][1]["path"] == "A317"
    assert out[3.0][0].sum() < out[1.85][0].sum()
    J = _inputs(HEAD19)
    for nz, t in ((2.6, 2.6), (3.2, 3.0)):
        keep, info = BP.filter_bad_points_info(J["P"], J["sr"].ravel(), 3.0, [100, 100, nz, 2.0])
        assert info["path"] == "A262" and info["tc"] == t, nz
        assert info["stage3"]["cnt"][3] == 19


@pytest.mark.parametrize("cam", [3.0, 2.1])
def test_the_filter_is_deterministic(cam):
    I = _inputs()  # noqa: E741
    a = BP.filter_bad_points(I["P"], I["sr"].ravel(), 3.0, [100, 100, 2.6, cam])
    b = BP.filter_bad_points(I["P"], I["sr"].ravel(), 3.0, [100, 100, 2.6, cam])
    assert np.array_equal(a, b)


# -- few points (task C02c's follow-up, probes fu1 and fu3 to fu3c) --------------------------------


def test_fewer_than_three_points_give_one_true():
    """MATLAB's filter returns a scalar true for one and two points, NaN among them (no MathWorks
    call is made), with which the detector keeps them"""
    P = np.array([[5.0, 5.0, 5.0], [9.0, 5.0, 5.0]])
    for sr in ([3.4], [np.nan, 3.3], [0.5, 40.0]):
        keep = BP.filter_bad_points(P[: len(sr)], np.array(sr), 3.0, [100, 100, 2.6, 3.0])
        assert keep.shape == (1,) and keep[0]


def test_the_short_path_is_taken_up_to_the_ratio_2():
    """prctile(sr, 95) <= 2 prctile(sr, 5) skips the histogram analysis (on 5 points the two
    percentiles are the extremes; the ratio 2 itself takes it, 2.0000000003 does not)"""
    assert BP.short_path(np.array([3.0, 4.0, 6.0, 5.0, 5.5]))
    assert not BP.short_path(np.array([3.0, 4.0, 6.000000001, 5.0, 5.5]))
    assert BP.short_path(np.full(5, 4.2))
    assert BP.short_path(np.array([2.0, np.nan, 4.0, 3.0]))
    assert not BP.short_path(designed_sigma_r())


def test_the_short_path_removes_as_a_with_thr_and_nan_below_every_threshold():
    """with the neighbour counts given (MATLAB's lists mocked in the probes): 1 neighbour removed;
    2 and 3 below thr removed; 2 to 8 below min(the noise threshold, thr) removed; a NaN sigma
    removed with 1 to 8 neighbours and kept with more"""
    sr = np.array(
        [2.0, 2.5, 2.9, 3.1, 3.5, 4.0, 2.2, 2.7, 3.3, 3.8, np.nan, np.nan, np.nan, np.nan]
    )
    cnt = np.array([1, 1, 2, 2, 2, 2, 3, 3, 5, 5, 1, 3, 8, 9])
    P = np.zeros((sr.size, 3))
    keep, info = BP.filter_bad_points_info(
        P, sr, 3.0, [100, 100, 2.65, 3.0], oracle={"counts": cnt}
    )
    assert info["path"] == "A318" and info["stage3"] is None
    assert keep.astype(int).tolist() == [0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1]
    keep, _ = BP.filter_bad_points_info(P, sr, 3.4, [100, 100, 4.0, 3.0], oracle={"counts": cnt})
    assert keep[8] == 0 and keep[9] == 1  # 3.3 below min(0.9 x 4, 3.4) with 5 neighbours
    # three neighbours: below min(thr, camera threshold, 2 noise) only (probe fu10: kept from 2.6
    # below thr 4 with camera 3.5 and noise 1); two: below thr
    s3 = np.array([2.6, 2.75, 3.2, 3.95, 2.7, 3.3])
    c3 = np.array([3, 3, 3, 3, 2, 2])
    keep, _ = BP.filter_bad_points_info(np.zeros((6, 3)), s3, 4.0, [100, 100, 1.0, 3.5],
                                        oracle={"counts": c3})  # fmt: skip
    assert keep.astype(int).tolist() == [1, 1, 1, 1, 0, 0]


def test_few_sigmas_below_the_99th_percentile_keep_every_point():
    """the long path with fewer than 10 sigmas below their 99th percentile returns every point
    (lines 58-60 to 115; 5 to 10 points kept whatever their neighbours, 11 went on)"""
    for n in (5, 10):
        sr = 3.0 + 13.0 * (np.arange(n) + 0.5) / n
        keep, info = BP.filter_bad_points_info(np.zeros((n, 3)), sr, 3.0, [100, 100, 2.6, 3.0])
        assert info["path"] == "early" and keep.size == n and keep.all()
    sr = 3.0 + 13.0 * (np.arange(11) + 0.5) / 11
    P = np.arange(33, dtype=np.float64).reshape(11, 3)
    _, info = BP.filter_bad_points_info(P, sr, 3.0, [100, 100, 2.6, 3.0])
    assert info["path"] != "early"


def test_nan_sigmas_on_the_long_path():
    """NaN left out of the percentiles, the histogram and the bin count (MATLAB's max and min
    leave the NaN of std out), and never removed for its sigma (probe fu3c: frame 28's call with
    five NaN, each point given two neighbours, kept them)"""
    sr = designed_sigma_r()
    srn = sr.copy()
    srn[[0, 100, 400, 401, 700]] = np.nan
    cnt = np.full(sr.size, 2)
    P = np.zeros((sr.size, 3))
    keep, info = BP.filter_bad_points_info(
        P, srn, 3.0, [100, 100, 2.6, 3.0], oracle={"counts": cnt}
    )
    assert info["path"] == "B" and info["nb0"] == 60
    assert keep[[0, 100, 400, 401, 700]].all()
    # MATLAB's min leaves NaN out: the range of the bin count is that of the finite sigmas (a
    # long tail: the Freedman-Diaconis count above 60)
    r = synth.CounterRNG(2026, "C02c.badpoints.tail")
    tail = np.r_[3.0 + np.abs(r.normal((1000,))) * 0.6, np.linspace(8.0, 30.0, 40)]
    nb_tail = BP.stage12(tail)[1]
    tail[[5, 50]] = np.nan
    assert nb_tail > 60 and BP.stage12(tail)[1] > 60


def test_line_160s_bound():
    """C at line 160 when the background peak's bin centre is at or above 2 x max(camera
    threshold, noise) while the noise is at most 2 camera thresholds (designed histograms: 2
    noise to 0.02 at noises 2.4 to 3.5, 2 cam at cameras 5 and 6; the demo's camera 3: 6), 2.5 x
    noise above that, never with a camera threshold of 1.25 or less below the noise"""
    assert BP.early_c_bound(2.65, 3.0) == 6.0
    assert BP.early_c_bound(2.65, 1.5) == 5.3
    assert BP.early_c_bound(3.5, 2.0) == 7.0
    assert BP.early_c_bound(2.65, 5.0) == 10.0
    assert BP.early_c_bound(3.0, 1.5) == 6.0  # noise at 2 cam: the first form (5.95 to 6.05)
    assert BP.early_c_bound(3.5, 1.5) == 8.75
    assert BP.early_c_bound(6.0, 2.5) == 15.0
    assert BP.early_c_bound(1.2, 1.25) == 2.5
    assert BP.early_c_bound(2.2, 1.25) == np.inf
    assert BP.early_c_bound(2.65, 1.0) == np.inf


def test_line_160s_bound_in_the_decision():
    """the designed histogram with its first bin's centre moved: C from the bound on, the first
    valley's analysis below it (camera threshold 1.5, noise 2.65: 5.3; C02's rule took C from 3)"""
    c = _hist(HEAD[:6])
    for centre, want in ((3.1, False), (5.2, False), (5.4, True)):
        e = E - E[0] + centre - 0.445 / 2
        path, _, info = BP.decide({"cnt": c, "e": e, "lst": [1, 6, 14, 39]}, 3.0, 2.65, 1.5)
        assert (path == "C") == want and abs(info["sig_bp"] - centre) < 1e-9, centre
