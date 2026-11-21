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
