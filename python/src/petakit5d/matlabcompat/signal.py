"""MATLAB's ``histcounts`` (with a bin count), ``findpeaks`` (with ``'MinPeakProminence'``) and
``medfilt1`` (order 3) in the forms the point detection calls them (R2024b), written from
MATLAB's documentation and black-box measurements of their outputs (task C02; no MathWorks
source was read): bit for bit on every recorded call of the detection's bad-point filter (136
histograms, 149 peak searches).
"""

from __future__ import annotations

import math

import numpy as np

from petakit5d.matlabcompat import elementary as _el

__all__ = ["findpeaks", "histcounts", "histcounts_edges", "medfilt1"]


def _pow10(k: int) -> float:
    """MATLAB's ``10^k`` (fdlibm's pow: ``10^-5`` is 9.999999999999999e-06, not the nearest
    double)"""
    return float(_el.power(np.array(10.0), np.array(float(k))))


def histcounts_edges(x, nbins: int) -> np.ndarray:
    """the bin edges of ``histcounts(x, nbins)`` for finite data with more than one value: the
    left edge is ``min(x)`` rounded down to a multiple of the raw width ``(max - min)/nbins``
    truncated to one significant digit; the width is the smallest multiple of the largest power
    of ten not above ``ul - ll`` that is at least ``ll``, with ``ll = (max - left)/nbins`` and
    ``ul = (max - left)/(nbins - 1)``; the last edge is at least ``max(x)``."""
    v = np.asarray(x, dtype=np.float64).ravel()
    nbins = int(nbins)
    if nbins < 1:
        raise ValueError("nbins must be a positive integer")
    xmin, xmax = float(v.min()), float(v.max())
    if not (math.isfinite(xmin) and math.isfinite(xmax)) or xmax == xmin:
        raise NotImplementedError("histcounts of constant or non-finite data is not ported")
    raw = (xmax - xmin) / nbins
    p = _pow10(math.floor(math.log10(raw)))
    u = p * math.floor(raw / p)
    left = min(u * math.floor(xmin / u), xmin)
    if nbins > 1:
        ll = (xmax - left) / nbins
        ul = (xmax - left) / (nbins - 1)
        p10 = _pow10(math.floor(math.log10(ul - ll)))
        w = p10 * math.ceil(ll / p10)
    else:
        w = u
    e = left + np.arange(nbins + 1, dtype=np.float64) * w
    e[-1] = max(e[-1], xmax)
    return e


def histcounts(x, nbins: int) -> tuple[np.ndarray, np.ndarray]:
    """``[N, edges] = histcounts(x, nbins)``: counts (double) of the bins ``[e_k, e_k+1)``, the
    last one closed, and the edges (``histcounts_edges``)"""
    v = np.asarray(x, dtype=np.float64).ravel()
    e = histcounts_edges(v, nbins)
    idx = np.clip(np.searchsorted(e, v, side="right") - 1, 0, int(nbins) - 1)
    return np.bincount(idx, minlength=int(nbins)).astype(np.float64), e


def findpeaks(y, min_prominence: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """``[pks, locs] = findpeaks(y, 'MinPeakProminence', m)`` of a finite vector: a peak is a
    sample above its left neighbour and above the next sample that differs from it (a flat top
    counts once, at its first sample; the ends are not peaks). Its prominence: from the peak,
    left and right until a sample higher than the peak (equal ones do not stop) or the end; the
    larger of the two intervals' minima is the reference. Peaks with prominence at least ``m``
    are kept; ``locs`` are 1-based (double)."""
    v = np.asarray(y, dtype=np.float64).ravel()
    n = v.size
    locs = []
    i = 1
    while i < n - 1:
        if v[i] > v[i - 1]:
            j = i
            while j + 1 < n and v[j + 1] == v[i]:
                j += 1
            if j + 1 < n and v[j + 1] < v[i]:
                locs.append(i)
            i = j + 1
        else:
            i += 1
    keep = []
    for i in locs:
        h = v[i]
        a, lmin = i - 1, h
        while a >= 0 and v[a] <= h:
            lmin = min(lmin, v[a])
            a -= 1
        b, rmin = i + 1, h
        while b < n and v[b] <= h:
            rmin = min(rmin, v[b])
            b += 1
        if h - max(lmin, rmin) >= min_prominence:
            keep.append(i)
    k = np.asarray(keep, dtype=np.int64)
    return v[k], k.astype(np.float64) + 1


def medfilt1(y) -> np.ndarray:
    """``medfilt1(y)``: the median of 3, zeros beyond the ends"""
    v = np.asarray(y, dtype=np.float64).ravel()
    p = np.concatenate([[0.0], v, [0.0]])
    return np.median(np.stack([p[:-2], p[1:-1], p[2:]]), axis=0)
