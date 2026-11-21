"""``filter_bad_points_in_detection(pos, sigma_r, thr, bgInfo, vol, A, c)``: which detections the
lowSNR detector keeps after its fits, learned as a black box (task C02; no GPL source read): the
residuals' sigma histogram (coarsened while ragged) gives a background peak and a valley between
it and the next peak; points with too few neighbours and a sigma below the valley, the residual
threshold or the noise threshold are removed. Measured: the keep mask of MATLAB on the demo's
120 frames with the background correction (both channels; channel 1 frame 28 takes the second
valley search, tasks C02 and C02c), on 15 of 15 traced calls, on 248 of 285 replayed
variants and on 118 designed calls of few points (task C02c). Inputs that matter: the positions
(through the neighbour counts), sigma_r, ``thr`` and ``bgInfo(3:4)`` (noise, camera threshold).

Stages (as measured): 0. fewer than three points: one true (MATLAB's result, with which the
detector keeps them); ``prctile(sr, 95) <= 2 prctile(sr, 5)``: straight to stage 5 and the
removal of the path A with thr (``short_path``); fewer than 10 sigmas below their 99th
percentile: every point kept (task C02c); 1. ``x`` = sigma_r below its 99th percentile,
``max(60, Freedman-Diaconis, ceil(10 (p99 - min)/std))`` bins (NaN left out), peaks of
prominence 5; 2. ``nb <- ceil(0.9 nb)`` while the histogram is ragged; 3. the background bin
search with optional smoothing; 4. the path (C: thr only, A: the valley or thr, B: the valley and
0.8 of it), with the second valley search (``second_valley``: B, or A262, A with min(noise, 1.5
camera thresholds)) on the A paths; 5. neighbour counts (radius mean + 3 std of the distance to
the second-nearest other point); 6. the removal rules.

Open (docs/parity/detection.md, the bad-point filter): the smoothing choice of line 104 (conv or
medfilt1), wrong on 3 of the demo's 60 frames without the background correction, and the B
path's threshold on one more; line 160's bound in the cases listed there.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.spatial import cKDTree


class _MW:
    """MATLAB's functions in the forms this filter calls them (petakit5d.matlabcompat)"""

    @staticmethod
    def histcounts(x, nbins):
        from petakit5d.matlabcompat import signal

        return signal.histcounts(x, nbins)

    @staticmethod
    def findpeaks(y, m):
        from petakit5d.matlabcompat import signal

        return signal.findpeaks(y, m)

    @staticmethod
    def medfilt1(y):
        from petakit5d.matlabcompat import signal

        return signal.medfilt1(y)

    @staticmethod
    def conv3_same(y, h):
        """``conv(y, h, 'same')`` of a 3-tap h, MATLAB's chain of fused multiply-adds:
        ``fma(h3, y(i-1), fma(h2, y(i), h1 y(i+1)))``, zeros beyond the ends"""
        from fractions import Fraction

        def fma(a, b, c):
            return float(Fraction(a) * Fraction(b) + Fraction(c))

        v = np.asarray(y, dtype=np.float64).ravel()
        w = [float(t) for t in np.asarray(h, dtype=np.float64).ravel()]
        p = np.concatenate([[0.0], v, [0.0]])
        return np.array(
            [fma(w[2], p[i], fma(w[1], p[i + 1], w[0] * p[i + 2])) for i in range(v.size)]
        )

    @staticmethod
    def prctile(x, p):
        from petakit5d.matlabcompat import stats

        return float(
            np.asarray(stats.prctile(np.asarray(x, dtype=np.float64).ravel(), p)).ravel()[0]
        )

    @staticmethod
    def median(x):
        from petakit5d.matlabcompat import stats

        return float(np.asarray(stats.median(np.asarray(x, dtype=np.float64).ravel())).ravel()[0])

    @staticmethod
    def std(x):
        from petakit5d.matlabcompat import stats

        return float(np.asarray(stats.std(np.asarray(x, dtype=np.float64).ravel())).ravel()[0])

    @staticmethod
    def mean(x):
        from petakit5d.matlabcompat import reduce

        v = np.asarray(x, dtype=np.float64).ravel()
        return float(np.asarray(reduce.sum(v)).ravel()[0]) / v.size


M = _MW


# ---------------------------------------------------------------- stage 1-2 (exact on all recorded)


def coarsen_cond(cnt, locs, nb):
    """lines 64 / 69: more than 20 peaks, or at least 4 gaps of at most 4 bins between the peaks
    below nb/4 (nb: the bin count asked of histcounts), or at least 3 pairs of equal neighbouring
    bins (counts(i) == counts(i+1)) whose count is above the 25th percentile of the counts"""
    locs = np.asarray(locs, float)
    if locs.size > 20:
        return True
    low = locs[locs < nb / 4]
    if int(np.sum(np.diff(low) <= 4)) >= 4:
        return True
    cnt = np.asarray(cnt, float)
    p25 = M.prctile(cnt, 25)
    z = np.flatnonzero(np.diff(cnt) == 0)
    return int(np.sum(cnt[z] > p25)) >= 3


def stage12(sr):
    """(x, nb0, [(nb, counts, edges, locs), ...]) the histogram sequence, the last one final"""
    sr = np.asarray(sr, float).ravel()
    p99 = M.prctile(sr, 99)
    x = sr[sr < p99]
    iqr = M.prctile(sr, 75) - M.prctile(sr, 25)
    rng = p99 - float(np.nanmin(sr))
    # MATLAB's max and min leave NaN out (a NaN sigma makes std NaN; prctile leaves it out)
    terms = [rng / (2 * iqr * sr.size ** (-1 / 3)), 10 * rng / M.std(sr)]
    terms = [t for t in terms if not math.isnan(t)]
    if any(math.isinf(t) for t in terms):
        raise NotImplementedError("the bad-point filter's bin count is infinite (no spread)")
    nb0 = max([60] + [math.ceil(t) for t in terms])
    seq = []
    nb = nb0
    cnt, e = M.histcounts(x, nb)
    _, locs = M.findpeaks(cnt, 5)
    seq.append((nb, cnt, e, locs))
    if coarsen_cond(cnt, locs, nb):
        while True:
            nb = math.ceil(0.9 * nb)
            cnt, e = M.histcounts(x, nb)
            _, locs = M.findpeaks(cnt, 5)
            seq.append((nb, cnt, e, locs))
            if not coarsen_cond(cnt, locs, nb):
                break
    return x, nb0, seq


# ---------------------------------------------------------------- stage 3


def bg_bin(cnt, e, pq):
    """the bin with the largest count among bins 1..round((pq - e(1)) / w) (first on ties)"""
    w = e[1] - e[0]
    k = int(round((pq - e[0]) / w))
    k = max(1, min(k, cnt.size))
    return int(np.argmax(cnt[:k])) + 1, k


def nbr3(b, nb):
    a = np.array([b - 1, b, b + 1])
    return a[(a >= 1) & (a <= nb)]


def smooth_rule(cnt, lst, pre, locs, nb):
    """(smooth, conv): lines 94, 99 and 104.
    Line 94 (exact on 1710 recorded evaluations): branch '95' when the entries prepended at lines
    85/89 plus the found peaks below nb/5 number at least 4, else branch '96'.
    Line 99, branch 96 (exact on 1154): smoothing when two neighbouring entries of the list (in list
    order, prepended entries first) among those below nb/5 are at most 3 bins apart (difference
    <= 3, negative differences included). Branch 95 (PROVISIONAL, 545 of 556 recorded): at least
    two differences <= 3 between the sorted entries below nb/5.
    Line 104 (PROVISIONAL, right on the 8 real-run frames that smooth): conv [0.1 0.8 0.1] when
    bin 1 is above the median of the list's counts, else medfilt1."""
    low = [p for p in lst if p < nb / 5]
    n94 = len(pre) + sum(1 for p in locs if p < nb / 5)
    if n94 >= 4:
        smooth = int(np.sum(np.diff(np.array(sorted(low), float)) <= 3)) >= 2
    else:
        smooth = bool(np.any(np.diff(np.array(low, float)) <= 3))
    if not smooth:
        return False, False
    med = M.median(np.asarray(cnt, float)[np.array(lst) - 1])
    return True, bool(cnt[0] > med)


def stage3(sr, nb, cnt, e, locs, smooth=None, conv=None):
    sr = np.asarray(sr, float).ravel()
    locs = [int(v) for v in locs]
    p30 = M.prctile(sr, 30)
    b, k = bg_bin(cnt, e, p30)
    lst = list(locs)
    pre = []
    if not np.any(np.isin(nbr3(b, nb), lst)):
        lst = [b] + lst
        pre.append(b)
    if cnt[0] > cnt[1] and 1 not in lst:
        lst = [1] + lst
        pre.append(1)
    if smooth is None:
        smooth, conv = smooth_rule(cnt, lst, pre, locs, nb)
    cur = np.asarray(cnt, float).copy()
    if smooth:
        cur = M.conv3_same(cnt, [0.1, 0.8, 0.1]) if conv else M.medfilt1(cnt)
        _, l2 = M.findpeaks(cur, 5)
        lst = [int(v) for v in l2]
    for q in (30.0, 20.0, 10.0):
        bq, _ = bg_bin(cur, e, M.prctile(sr, q))
        if not np.any(np.isin(nbr3(bq, nb), lst)):
            nxt = cur[bq] if bq < cur.size else -np.inf
            if cur[bq - 1] >= nxt:
                lst = [bq] + lst
    if cur[0] > cur[1] and 1 not in lst:
        lst = [1] + lst
    return {
        "cnt": cur,
        "raw": np.asarray(cnt, float),
        "e": e,
        "lst": lst,
        "b30": b,
        "k30": k,
        "smooth": bool(smooth),
        "conv": bool(conv),
    }


# ---------------------------------------------------------------- stage 4


def valley(cnt, e, lst):
    """background peak bp (bin 1 when above bin 2, else the first listed peak), the next listed
    peak after it, the lowest bin between them (first on ties) and its centre T_c"""
    w = e[1] - e[0]
    pl = sorted({int(p) for p in lst})
    bp = 1 if cnt[0] > cnt[1] else (pl[0] if pl else 1)
    later = [p for p in pl if p > bp]
    nxt = later[0] if later else cnt.size
    kv = bp + int(np.argmin(cnt[bp - 1 : nxt]))
    return bp, kv, nxt, e[0] + (kv - 0.5) * w


def early_c_bound(noise, cam):
    """line 160's bound on the background peak's bin centre (C at or above it; task C02c:
    frame 28's histogram with its edges moved, 29 groups of camera threshold, noise and thr, and
    5475 recorded calls): 2 x max(camera threshold, noise) while the noise is at most 2 camera
    thresholds (to 0.02 at noises 2.4 to 3.5 and cameras 1.5 to 6), 2.5 x noise above that (line
    17's case: between 8.02 and 9.02 at noise 3.5, camera 1.5; 13.02 and 15.02 at 6, 2.5; above
    7.5 at 3.01, 1.5); with a camera threshold of 1.25 or less and the noise above it, no C (A up
    to the largest centres designed: 5.6 at camera 1.25 and noise 2.2, 33 at camera 1; where
    between 1.25 and 1.5 that ends is open); thr does not enter (2, 3 and 4 alike)"""
    if cam < 1.5 and noise > cam:
        return math.inf
    if noise > 2 * cam:
        return 2.5 * noise
    return 2 * max(cam, noise)


def decide(st, thr, noise, cam):
    """the removal path ('A317': A with T_c, 'A318': A with thr, 'B', 'C', and the second
    search's 'A262': A with its threshold) and its threshold T_c. Measured on synthetic
    histograms (probes e8/e9) and checked on the recorded calls:
      C when the list has one entry, or the background peak is low (cbp <= 0.34 x prctile 85 of
      the counts: the exact form is open), or its bin centre is at or above line 160's bound
      (``early_c_bound``: 2 x max(camera threshold, noise) for the demo's), or the valley is
      shallow (ckv >= 0.6 cbp, line 196); (a close, taller next peak also gives C, lines
      178-179: not implemented, no real-run frame has it);
      B when the valley is below 0.15 of the background peak and T_c < 2 x camera threshold
      (camera threshold not 5: two designed histograms with 5 took A with T_c, task C02c);
      otherwise (the A paths) the second valley search's decision when it makes one
      (``second_valley``, lines 221-267: B or A262);
      A with T_c when the valley is below 0.3 of the background peak and
      max(T_noise, thr) < T_c < max(2.5 noise, 2 camera threshold), else A with thr.
    Not covered: the 2-entry list branch (148-149), the camera threshold >= 5 and noise > 2 x
    camera branches."""
    cur, lst, e = st["cnt"], st["lst"], st["e"]
    if len(lst) < 2:
        return "C", None, {}
    bp, kv, nxt, tc = valley(cur, e, lst)
    w = e[1] - e[0]
    cbp, ckv = cur[bp - 1], cur[kv - 1]
    sig = e[0] + (bp - 0.5) * w
    p85 = M.prctile(cur, 85) if len(lst) > 2 else np.nan
    tn = max(0.9 * noise, min(noise, cam))
    info = {
        "bp": bp,
        "kv": kv,
        "nxt": nxt,
        "tc": tc,
        "cbp": cbp,
        "ckv": ckv,
        "sig_bp": sig,
        "p85": p85,
    }
    if len(lst) > 2 and not cbp > 0.34 * p85:
        return "C", tc, info
    if not sig < early_c_bound(noise, cam):
        return "C", tc, info
    if ckv >= 0.6 * cbp:
        return "C", tc, info
    if ckv < 0.15 * cbp and cam != 5 and tc < 2 * cam:
        return "B", tc, info
    sv = second_valley(cur, e, bp, nxt, cam, noise, info)
    if sv is not None:
        return sv[0], sv[1], info
    if ckv < 0.3 * cbp and tc > max(tn, thr) and tc < max(2.5 * noise, 2 * cam):
        return "A317", tc, info
    return "A318", tc, info


# the second valley search's levels: the background peak's count over these
SECOND_DIVISORS = (4.0, 3.5, 3.0, 2.5, 2.0)


def _msum(v) -> float:
    from petakit5d.matlabcompat import reduce

    v = np.asarray(v, dtype=np.float64).ravel()
    return float(np.asarray(reduce.sum(v)).ravel()[0]) if v.size else 0.0


def second_trigger(cur, bp: int) -> bool:
    """lines 221-222: the second valley search runs when the bins around the background peak
    hold a large share of all counts: the three bins from max(bp - 1, 1) more than 0.12 of
    them (line 221), or the first two of those more than 0.09 (line 222; 786 of 786 designed
    and real calls)"""
    cur = np.asarray(cur, dtype=np.float64)
    total = _msum(cur)
    if not total > 0:
        return False
    lo = max(bp - 1, 1) - 1
    return _msum(cur[lo : lo + 3]) / total > 0.12 or _msum(cur[lo : lo + 2]) / total > 0.09


def second_valley(
    cur, e, bp: int, nxt: int, cam: float, noise: float = 0.0, info: dict | None = None
):
    """the second valley search (lines 221-267; reached on the A paths, after the first valley's
    B test): (path, threshold) of the decision it makes, or None to keep the first valley's.
    Measured with about 620 designed histograms and the demo's calls (MATLAB's line numbers and
    keep masks only, no source read):

    - when ``second_trigger``: for the levels cbp / d, d = 4, 3.5, 3, 2.5, 2 in turn (lines
      226-229), the first bin j after the background peak and before the next listed peak
      whose count is below the level (strict); T2 = its upper edge e(j + 1) (line 232);
    - ('B', T2) when j was found at one of the first four levels and T2 < 2 cam (cam, the
      camera threshold, not 5; lines 238-243; for cam above 5 the bound is open: 8.585 at 5.5
      gave B);
    - otherwise, for cam below 5 and the first bin's centre below 1.5 cam (lines 247 and 257;
      the background peak was bin 1 in the designs, a recorded call with it at bin 9 or 10 took
      the lines with the first bin's centre below 1.5 cam, and one with cam 6 did not):
      ('B', T2) when the count of bin j is below cbp / 5 (lines 248-252), else, with j found at
      one of the first three levels (line 258), ('A262', min(noise, 1.5 cam)) (line 262: A's
      removal with that threshold, the noise rule's threshold at most it, B's median bound);
    - otherwise the first valley's decision stands (with cam 5 MATLAB enters both blocks and
      keeps it: lines 249-250, 259-260)."""
    cur = np.asarray(cur, dtype=np.float64)
    if not second_trigger(cur, bp):
        return None
    cbp = cur[bp - 1]
    found = None
    for i, d in enumerate(SECOND_DIVISORS, start=1):
        lv = cbp / d
        for j in range(bp + 1, nxt):  # 1-based bins bp+1 .. nxt-1
            if cur[j - 1] < lv:
                found = (i, j)
                break
        if found:
            break
    if info is not None:
        info["second"] = found
    if found is None:
        return None
    level, j = found
    t2 = float(e[j])  # e(j + 1), 1-based
    if level <= 4 and cam != 5 and t2 < 2 * cam:
        return "B", t2
    s1 = e[0] + (e[1] - e[0]) / 2  # the first bin's centre
    if cam >= 5 or not s1 < 1.5 * cam:
        return None
    if cur[j - 1] < cbp / 5:
        return "B", t2
    if level <= 3:
        return "A262", min(noise, 1.5 * cam)
    return None


# ---------------------------------------------------------------- stage 5-6


def neighbour_counts(pos):
    """(R, counts, d3): R = mean + 3 std (MATLAB's mean, std) of the third column of
    knnsearch(pos, pos, 'K', 3); counts of points within R, the point itself included"""
    P = np.asarray(pos, float)
    tree = cKDTree(P)
    d, _ = tree.query(P, k=3)
    d3 = d[:, 2]
    R = M.mean(d3) + 3 * M.std(d3)
    cnt = np.array([len(v) for v in tree.query_ball_point(P, R)])
    return R, cnt, d3


def noise_threshold(noise, cam):
    return max(0.9 * noise, min(noise, cam))


def removal(sr, cnt, thr, noise, cam, path, tc, U=8, short=False):
    """the keep mask of the removal stage (``short``: the short path's A with thr). With noise >
    2 x camera threshold (line 17; PROVISIONAL, not in the real-run data): A with thr removes
    counts 2..10 below thr and neither A with thr nor C applies the noise rule"""
    sr = np.asarray(sr, float).ravel()
    tn = noise_threshold(noise, cam)
    flag17 = cam < 5 and noise > 2 * cam
    rem = cnt == 1
    if path == "C":
        rem |= (cnt == 2) & (sr < thr)
    elif path == "B":
        rem |= (cnt >= 2) & (cnt <= 4) & (sr < tc)
        rem |= (cnt >= 5) & (cnt <= U) & (sr < 0.8 * tc)
    elif short and not flag17:
        # the short path's A with thr: 2 neighbours below thr, 3 below min(thr, cam, 2 noise)
        # (task C02c, probes sp1, fu9 and fu10: three neighbours kept below thr 4 with cameras 3
        # to 5.5 and noise 1 or 2.6, removed below thr 3 with camera 3 and noises 1.5 to 4; the
        # long path removes them below thr, camera 2 included: C02's recorded variants)
        rem |= (cnt == 2) & (sr < thr)
        rem |= (cnt == 3) & (sr < min(thr, cam, 2 * noise))
    else:
        t1 = tc if path in ("A317", "A262") else thr
        rem |= (cnt >= 2) & (cnt <= (10 if flag17 and path == "A318" else 3)) & (sr < t1)
    if not (flag17 and path in ("A318", "C")):
        if path in ("A318", "C"):
            tn = min(tn, thr)
        elif path == "A262":  # the second search's A: the noise rule at most its threshold
            tn = min(tn, tc)
        rem |= (cnt >= 2) & (cnt <= U) & (sr < tn)
    return ~rem


def median_bound(sr, cnt, path, tc):
    """U, the upper neighbour count of the B class-2 and noise rules: min(median of the
    candidates' counts, 10) when MATLAB takes that median, else 8. Candidates: sr < T_c.
    B and A262 (line 352; A262 runs B's lines 351-357; exact on the 72 recorded B calls and 18
    of the 19 B frames of the demo without the background correction): the median is taken
    when there are more than 100 candidates and more than half of them have a count above 8.
    Open: that frame (103 candidates, median count 8, 51 above 8) took it with U 9 or 10.
    A (line 359, PROVISIONAL): taken when there are more than 100 candidates (right on every
    real-run frame; the rule is open for low-noise and A-with-thr calls). C: never."""
    if path == "C" or tc is None:
        return 8
    sr = np.asarray(sr, float).ravel()
    cand = sr < tc
    nc = int(np.sum(cand))
    if path in ("B", "A262"):  # A262: the lines of B's median (351-357)
        take = nc > 100 and int(np.sum(cnt[cand] > 8)) > nc / 2
    else:
        take = nc > 100
    if not take:
        return 8
    return min(M.median(cnt[cand]), 10)


def short_path(sr) -> bool:
    """lines 34-39: the histogram analysis (lines 40-279) is skipped and the points are removed
    as on the path A with thr (lines 273 and 318) when ``prctile(sr, 95) <= 2 prctile(sr, 5)``
    (task C02c, probes fu3, fu3b and fu3c: taken on 78 calls with ratios up to 2 itself, not on
    47 from 2.0000000003, whatever the spread)"""
    p95, p5 = M.prctile(sr, 95), M.prctile(sr, 5)
    return p95 <= 2 * p5


# the long path returns every point when fewer sigmas than this lie below their 99th percentile
FEW_BELOW_P99 = 10

# the information of a call that returns before the analysis
_NO_ANALYSIS = {
    "tc": None,
    "U": None,
    "nb0": None,
    "nsteps": 0,
    "nb": None,
    "stage3": None,
    "decide": {},
    "counts": None,
}


def filter_bad_points_info(pos, sigma_r, thr, bginfo, vol=None, A=None, c=None, *, oracle=None):
    sr = np.asarray(sigma_r, float).ravel()
    thr = float(np.asarray(thr).ravel()[0])
    bg = np.asarray(bginfo, float).ravel()
    noise, cam = float(bg[2]), float(bg[3])
    oracle = oracle or {}
    if sr.size < 3:
        # MATLAB returns a scalar true for fewer than three points (task C02c, probes fu1 and
        # fu3: one and two points, NaN among them; no MathWorks call is made), with which the
        # detector keeps them all
        return np.ones(1, dtype=bool), {**_NO_ANALYSIS, "path": "few"}
    short = short_path(sr)
    if short:
        nb0, seq, st, dinfo = None, [], None, {}
        nb, path, tc = None, "A318", None
    else:
        x = sr[sr < M.prctile(sr, 99)]
        if x.size < FEW_BELOW_P99:
            # lines 58-60 -> 113-115: every point kept, no neighbour count made (task C02c,
            # probes fu3b and fu3c: 5 to 10 points with ratios from 2 kept all, 11 went on)
            return np.ones(sr.size, dtype=bool), {**_NO_ANALYSIS, "path": "early"}
        x, nb0, seq = stage12(sr)
        nb, cnt, e, locs = seq[-1]
        st = stage3(sr, nb, cnt, e, locs, smooth=oracle.get("smooth"), conv=oracle.get("conv"))
        path, tc, dinfo = decide(st, thr, noise, cam)
    if "path" in oracle:
        path = oracle["path"]
    ncnt = oracle["counts"] if "counts" in oracle else neighbour_counts(pos)[1]
    U = oracle.get("U", median_bound(sr, ncnt, path, tc))
    # the short path counts a NaN sigma below every threshold (lines 294-310: removed with 1 to
    # U neighbours, kept with more; probe fu3c); the long path never removes it by its sigma
    keep = removal(np.where(np.isnan(sr), -np.inf, sr) if short else sr, ncnt, thr, noise, cam,
                   path, tc, U=U, short=short)  # fmt: skip
    return keep, {
        "nb0": nb0,
        "nsteps": len(seq),
        "nb": nb,
        "path": path,
        "tc": tc,
        "U": U,
        "stage3": st,
        "decide": dinfo,
        "counts": ncnt,
    }


def filter_bad_points(pos, sigma_r, thr, bginfo, vol=None, A=None, c=None):
    """the keep mask (bool, one per point); see the module docstring for what is provisional"""
    return filter_bad_points_info(pos, sigma_r, thr, bginfo, vol, A, c)[0]
