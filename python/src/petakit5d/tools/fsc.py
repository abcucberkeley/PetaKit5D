"""Fourier shell correlation (MATLAB's ``XR_FSC``, ``XR_FSC_resolution``,
``XR_one_image_FSC_analysis``, ``XR_one_image_FSC_analysis_frame``, ``XR_FSC_analysis_wrapper``
and ``XR_two_image_FSC_analysis``): the resolution of an image per angle in the x-z plane, from
the correlation of two volumes' Fourier transforms over shells of radius and sectors of angle.

``fourier_shell_correlation(im1, im2)`` (``XR_FSC``), as MATLAB computes it:

- each volume zero-padded (``floor`` before, ``ceil`` after) and cropped about
  ``round((size + 1) / 2)`` to ``half_size``, multiplied by the cube root of the product of three
  Hamming windows (``(h_y h_x) h_z``, fdlibm's ``pow``), Fourier transformed and shifted;
- the spectrum split into shells by radius (``round((r + dr / 2) / dr)``) and sectors by the
  angle about the y axis (``round(az / dtheta) + 1``, the sector at pi merged into the one at 0,
  everything beyond ``(N - 1) / 2 + 1`` in one outer shell), with the cones along z (within
  pi/36) and x (within pi/48) left out under ``skip_cone_region``; the grid follows MATLAB's
  ``meshgrid``, ``cart2sph`` and ``permute(azimuth, [3, 1, 2])``;
- per cell ``real(sum(F1 conj(F2))) / sqrt(sum(|F1|^2) sum(|F2|^2))``, summed in MATLAB's
  linear order (``accumarray``), MATLAB's ``abs`` of complex values (msun's ``hypot``); empty
  cells and the last remaining one dropped; ``fsc['fsc']`` (radius x angle) with 1 where the
  first radius of an angle is missing; ``fsc['r']`` the radii in units of the Nyquist frequency,
  ``fsc['thetas']`` the angles, ``fsc['npoints']`` the points per cell.

With the exact FFT backend (the FFTW add-on) the curves are MATLAB's bit for bit; with the
permissive one (pocketfft, the default) within decision record 0006's double tolerance.

``fsc_resolution(fsc)`` (``XR_FSC_resolution``): per angle the first radius on a grid of 1e-4
(from the first radius to 1, then to sqrt(2)) where the smoothed curve falls below the threshold
(``'fixed'``) or below the smoothed one-bit or half-bit curve of the points per cell, and the
resolution ``2 psz / radius`` with ``psz = xy + (dz - xy) |cos(theta)|``, put on a full circle
(0 for x, pi/2 for z) by the angles' sectors (those the FSC lacks, inside the cones, are NaN; an
odd number of sectors in the half circle is refused, ``XR_FSC_resolution:oddSectorNumber``). A
curve that never falls below it is fitted with a quadratic; failing that the radius is sqrt(2).
Images without signal (all zero) have no FSC (``XR_FSC:noSignal``).

**The smoothing is not MathWorks's code.** MATLAB smooths with the Curve Fitting Toolbox's
``fit(r, fsc, 'smoothingSpline')``, which its documentation says is based on ``csaps``. Here it is
the cubic smoothing spline of the published method (Reinsch 1967; de Boor, *A Practical Guide
to Splines*, ch. XIV, whose SMOOTH csaps implements; the algorithm as Green and Silverman 1994,
section 2.3, give it): minimizing ``p sum (y - f)^2 + (1 - p) int f''^2`` with csaps's documented
default smoothing parameter (``p trace(A) = (1 - p) trace(B)`` for the system matrix
``p A + (1 - p) B``: ``default_smoothing``); evaluated as a piecewise cubic, past the last point
by the last piece. The quadratic is the least-squares fit, solved in exact rational arithmetic.
The radii found are MATLAB's or a grid step or two off (the last bits of the spline differ);
the tolerance on the resolution values is approved in ``docs/parity/tools.md`` (N3). The arithmetic
is plain IEEE operations in a fixed order, so the results are the same on every platform.

``one_image_fsc_analysis(fn)`` (``XR_one_image_FSC_analysis``): the image (a tiff file or a
``(z, y, x)`` array) clipped at a percentile (``clip_per``), cropped to ``bbox``, permuted for
``res_axis`` (``'yz'`` swaps x and y; ``'xy'`` makes MATLAB's ``permute(im, [3, 2, 1])``, with the
pixel sizes of the permuted axes and the FSC's plane of x and y), resized to isotropic voxels in
double (``imresize3_fast``, linear, only the region kept) and cropped or padded to twice
``half_size``; the FSC of four pairs of its interleaved sub-lattices (every other voxel along each
axis) and the resolution of their mean. ``two_image_fsc_analysis`` (``XR_two_image_FSC_analysis``):
two images resized (``imresize3``, cubic), cropped about their centre or padded to twice ``n``, and
the mean FSC of their eight same sub-lattices.

``fsc_analysis_file`` (``XR_one_image_FSC_analysis_frame``) writes the analysis of an image (of
each box of ``bbox``) to a v7.3 MAT-file under a temporary name and renames it.
``fsc_analysis_dataset`` (``XR_FSC_analysis_wrapper``) runs it for the tiff files of datasets
through the job framework (with ``multi_regions``, regions about each file's own centre), saves
``fsc_res_info_c<c>.mat`` per channel (the mean resolutions leave out the NaN sectors) and, with
matplotlib (the optional ``figures`` dependency), the figures as PNG (MATLAB's ``.fig`` files are
not written). ``fsc['npoints']`` and ``fsc['fsc']`` hold the angles the FSC found
(``fsc['thetas']``), column by column.

Kept from MATLAB:

- the crop to twice ``half_size`` starts at ``max(1, floor((size - 2 h) / 2))``, one voxel before
  the centred box (the two-image analysis's at the centred box);
- ``XR_FSC_resolution`` without ``dtheta`` takes the interval of the first two angles (a single
  angle fails);
- ``XR_FSC_analysis_wrapper`` filters by channel pattern on the full paths.

Not kept: MATLAB keeps the last image read by file name in a persistent variable; here a file
is read on each call.
"""

from __future__ import annotations

import math
import os
import threading
from collections.abc import Sequence
from fractions import Fraction
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "FSCError",
    "fourier_shell_correlation",
    "fsc_analysis_dataset",
    "fsc_analysis_file",
    "fsc_resolution",
    "hamming",
    "one_image_fsc_analysis",
    "poly2_fit",
    "smoothing_spline",
    "two_image_fsc_analysis",
]

_GRAIN = 1 << 18  # elements per piece of the threaded element-wise passes


class FSCError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


class FSCWarning(UserWarning):
    """a warning MATLAB gives, with its identifier"""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# MATLAB's pieces: colon, the Hamming window, abs of complex values
# ---------------------------------------------------------------------------------------------


def _colon(base: float, inc: float, limit: float) -> np.ndarray:
    """MATLAB's ``base : inc : limit`` of doubles (as ``params.literal``: the values symmetric
    about the midpoint, the last one the limit within 2 eps)"""
    from petakit5d.core.dtypes import mround

    base, inc, limit = float(base), float(inc), float(limit)
    if inc == 0 or (base < limit and inc < 0) or (limit < base and inc > 0):
        return np.zeros(0)
    eps = float(np.finfo(np.float64).eps)
    tol = 2.0 * eps * max(abs(base), abs(limit))
    sig = 1.0 if inc > 0 else -1.0
    if base == math.floor(base) and inc == math.floor(inc):
        n = int(math.floor((limit - base) / inc))
    else:
        n = int(mround((limit - base) / inc))
        if sig * (base + n * inc - limit) > tol:
            n -= 1
    last = base + n * inc
    if sig * (last - limit) >= -tol:
        last = limit
    k = np.arange(n // 2 + 1, dtype=np.float64)
    out = np.empty(n + 1)
    out[: n // 2 + 1] = base + k * inc
    out[n - np.arange(n // 2 + 1)] = last - k * inc
    if n % 2 == 0:
        out[n // 2] = (base + last) / 2
    return out


def hamming(n: int) -> np.ndarray:
    """MATLAB's ``hamming(n)`` (symmetric; a column, as a 1-D array): ``0.54 - 0.46 cos(2 pi x)``
    with ``x = (0 : m - 1) / (n - 1)`` on the first ``m = ceil(n / 2)`` points (fdlibm's cosine),
    mirrored; ``[1]`` for ``n = 1``"""
    from petakit5d.matlabcompat.elementary import cos

    n = int(n)
    if n <= 0:
        return np.zeros(0)
    if n == 1:
        return np.ones(1)
    m = (n + 1) // 2
    x = np.arange(m, dtype=np.float64) / (n - 1)
    w = 0.54 - 0.46 * cos(2 * np.pi * x)
    return np.concatenate([w, w[: n - m][::-1]])


def _threads(num_threads: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(num_threads)


def _pieces(fn, n: int, nt: int, grain: int = _GRAIN) -> None:
    """``fn(lo, hi)`` over [0, n) in pieces of ``grain``, on ``nt`` threads (NumPy releases the
    GIL); the result does not depend on the pieces"""
    starts = list(range(0, n, grain))
    if nt <= 1 or len(starts) <= 1:
        for lo in starts:
            fn(lo, min(n, lo + grain))
        return
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(min(nt, len(starts))) as ex:
        for f in [ex.submit(fn, lo, min(n, lo + grain)) for lo in starts]:
            f.result()


_HI_MASK = np.uint64(0xFFFFFFFF00000000)


def _abs_squared(re: np.ndarray, im: np.ndarray) -> np.ndarray:
    """``abs(complex(re, im)) .^ 2`` as MATLAB computes it: msun's ``hypot`` (as
    ``matlabcompat.elementary.hypot``; its medium-range path written out here, the rest through
    it), squared"""
    from petakit5d.matlabcompat.elementary import hypot

    xa, ya = np.abs(re), np.abs(im)
    hx = (xa.view(np.uint64) >> np.uint64(32)).astype(np.int64)
    hy = (ya.view(np.uint64) >> np.uint64(32)).astype(np.int64)
    sw = hy > hx
    a = np.where(sw, ya, xa)
    b = np.where(sw, xa, ya)
    ha = np.where(sw, hy, hx)
    hb = np.where(sw, hx, hy)
    special = ((ha - hb) > 0x3C00000) | (ha > 0x5F300000) | (hb < 0x23D00000)
    with np.errstate(invalid="ignore", over="ignore"):
        w = a - b
        t1 = (a.view(np.uint64) & _HI_MASK).view(np.float64)
        t2 = a - t1
        far = np.sqrt(t1 * t1 - (b * (-b) - t2 * (a + t1)))
        y1 = (b.view(np.uint64) & _HI_MASK).view(np.float64)
        y2 = b - y1
        t1 = ((ha + 0x00100000).astype(np.uint64) << np.uint64(32)).view(np.float64)
        t2 = (a + a) - t1
        near = np.sqrt(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)))
        h = np.where(w > b, far, near)
    if special.any():
        h[special] = hypot(re[special], im[special])
    with np.errstate(over="ignore"):
        return h * h


def _hi_words(x: np.ndarray) -> np.ndarray:
    return (x.view(np.uint64) >> np.uint64(32)).astype(np.int64)


def _from_words(hi: np.ndarray, lo) -> np.ndarray:
    h = (np.asarray(hi, dtype=np.int64) & 0xFFFFFFFF).astype(np.uint64) << np.uint64(32)
    return (h | np.asarray(lo, dtype=np.uint64)).view(np.float64)


def _lo0(x: np.ndarray) -> np.ndarray:
    return (x.view(np.uint64) & _HI_MASK).view(np.float64)


def _pow_positive(x: np.ndarray, y: float) -> np.ndarray:
    """fdlibm's ``pow(x, y)`` (``matlabcompat.elementary.power``, MATLAB's ``.^``) for positive
    normal finite ``x`` and an exponent with ``2^-1022 < |x^y| < 2^1023``, element-wise in NumPy:
    e_pow.c's path for finite arguments, one IEEE operation per C operation (the scalar port of
    elementary.py made element-wise; the other cases through it)"""
    from petakit5d.matlabcompat import elementary as E

    x = np.ascontiguousarray(x, dtype=np.float64)
    y = float(y)
    hx = _hi_words(x)
    lx = x.view(np.uint64) & np.uint64(0xFFFFFFFF)
    iy = (_hi_words(np.asarray([y]))[0]) & 0x7FFFFFFF
    ok = (hx >= 0x00100000) & (hx < 0x7FF00000) & ~((hx == 0x3FF00000) & (lx == 0))
    out = np.empty(x.shape)
    if (~ok).any():
        out[~ok] = E.power(x[~ok], y)
    if iy > 0x41E00000 or not ok.any() or not math.isfinite(y) or y in (0.0, 1.0, -1.0, 2.0, 0.5):
        if ok.any():
            out[ok] = E.power(x[ok], y)
        return out
    xs, ix, lo = x[ok], hx[ok], lx[ok]
    n = (ix >> 20) - 0x3FF
    j = ix & 0x000FFFFF
    ix = j | 0x3FF00000
    k = ((j > 0x3988E) & (j < 0xBB67A)).astype(np.int64)
    up = j >= 0xBB67A
    n = n + up
    ix = ix - up * 0x00100000
    ax = _from_words(ix, lo)
    bp = np.where(k == 1, E._BP[1], E._BP[0])
    dp_h = np.where(k == 1, E._DP_H[1], E._DP_H[0])
    dp_l = np.where(k == 1, E._DP_L[1], E._DP_L[0])
    L = E._L
    u = ax - bp
    v = 1.0 / (ax + bp)
    ss = u * v
    s_h = _lo0(ss)
    t_h = _from_words(((ix >> 1) | 0x20000000) + 0x00080000 + (k << 18), 0)
    t_l = ax - (t_h - bp)
    s_l = v * ((u - s_h * t_h) - s_h * t_l)
    s2 = ss * ss
    r = s2 * s2 * (L[0] + s2 * (L[1] + s2 * (L[2] + s2 * (L[3] + s2 * (L[4] + s2 * L[5])))))
    r = r + s_l * (s_h + ss)
    s2 = s_h * s_h
    t_h = _lo0(3.0 + s2 + r)
    t_l = r - ((t_h - 3.0) - s2)
    u = s_h * t_h
    v = s_l * t_h + t_l * ss
    p_h = _lo0(u + v)
    p_l = v - (p_h - u)
    z_h = E._CP_H * p_h
    z_l = E._CP_L * p_h + p_l * E._CP + dp_l
    t = n.astype(np.float64)
    t1 = _lo0(((z_h + z_l) + dp_h) + t)
    t2 = z_l - (((t1 - t) - dp_h) - z_h)
    y1 = float(_lo0(np.asarray([y]))[0])
    p_l = (y - y1) * t1 + y * t2
    p_h = y1 * t1
    z = p_l + p_h
    jz = _hi_words(z)
    jz = np.where(jz >= 0x80000000, jz - 0x100000000, jz)  # C's signed high word
    i = jz & 0x7FFFFFFF
    edge = (jz >= 0x40900000) | (i >= 0x4090CC00)
    kk = (i >> 20) - 0x3FF
    big = i > 0x3FE00000
    nn = np.where(big, jz + (0x00100000 >> np.clip(kk + 1, 0, 31)), 0)
    k2 = ((nn & 0x7FFFFFFF) >> 20) - 0x3FF
    tt = _from_words(nn & ~(0x000FFFFF >> np.clip(k2, 0, 31)), 0)
    nn = np.where(big, ((nn & 0x000FFFFF) | 0x00100000) >> np.clip(20 - k2, 0, 63), 0)
    nn = np.where(jz < 0, -nn, nn)
    p_h = np.where(big, p_h - tt, p_h)
    t = _lo0(p_l + p_h)
    u = t * E._LG2_H
    v = (p_l - (t - p_h)) * E._LG2 + t * E._LG2_L
    z = u + v
    w = v - (z - u)
    t = z * z
    P = E._PP
    t1 = z - t * (P[0] + t * (P[1] + t * (P[2] + t * (P[3] + t * P[4]))))
    r = (z * t1) / (t1 - 2.0) - (w + z * w)
    z = 1.0 - (r - z)
    j2 = _hi_words(z) + (nn << 20)
    edge |= (j2 >> 20) <= 0
    res = _from_words(j2, z.view(np.uint64) & np.uint64(0xFFFFFFFF))
    if edge.any():
        res[edge] = E.power(xs[edge], y)
    out[ok] = res
    return out


# ---------------------------------------------------------------------------------------------
# the window and the shells (kept between calls, as MATLAB keeps them in persistent variables)
# ---------------------------------------------------------------------------------------------

_LOCK = threading.Lock()
_WINDOW: dict[int, np.ndarray] = {}
_SHELLS: dict[tuple, tuple] = {}


def _window(n: int) -> np.ndarray:
    """``(hwd .* hwd' .* reshape(hwd, 1, 1, n)) .^ (1 / 3)`` as a (z, y, x) array: ``(h_y h_x)
    h_z`` raised by fdlibm's pow; computed on the first half of each axis and mirrored (the
    window is symmetric, so the products are)"""
    with _LOCK:
        w = _WINDOW.get(n)
    if w is not None:
        return w
    h = hamming(n)
    if np.array_equal(h, h[::-1]):
        m = (n + 1) // 2
        hh = h[:m]
        p = (hh[None, :, None] * hh[None, None, :]) * hh[:, None, None]
        idx = np.minimum(np.arange(n), n - 1 - np.arange(n))
        w = np.ascontiguousarray(_pow_positive(p, 1.0 / 3.0)[np.ix_(idx, idx, idx)])
    else:
        p = (h[None, :, None] * h[None, None, :]) * h[:, None, None]
        w = _pow_positive(p, 1.0 / 3.0)
    w.setflags(write=False)
    with _LOCK:
        _WINDOW.clear()
        _WINDOW[n] = w
    return w


def _grid_values(n: int) -> np.ndarray:
    """MATLAB's ``-(n - 1) / 2 : (n - 1) / 2``"""
    return _colon(-(n - 1) / 2, 1.0, (n - 1) / 2)


def _shells(n: int, dr: float, dtheta: float, skip_cone_region: bool, nt: int) -> tuple:
    """(the voxels' positions in MATLAB's linear order sorted by cell, stably; where each cell
    starts in them; the points per cell; max_r; max_e) for a cube of size n"""
    key = (n, float(dr), float(dtheta), bool(skip_cone_region))
    with _LOCK:
        got = _SHELLS.get(key)
    if got is not None:
        return got
    from petakit5d.core.dtypes import mround
    from petakit5d.matlabcompat.elementary import atan2, hypot

    c = _grid_values(n)
    # MATLAB's order (y fastest, then x, then z) is the C order of a (z, x, y) array; with
    # X(i, j, k) = c(j), Y = c(i), Z = c(k): r(i, j, k) = hypot(hypot(c(j), c(i)), c(k)) and,
    # after permute(azimuth, [3, 1, 2]), azimuth(i, j, k) = atan2(c(j), c(k))
    hxy = hypot(c[:, None], c[None, :])  # [x j, y i]
    az = atan2(c[None, :], c[:, None])  # [z k, x j]
    az = np.pi * (az < 0) + az
    rr = np.empty((n, n, n))

    def fill_r(lo, hi):
        k0, k1 = lo // (n * n), -(-hi // (n * n))
        for k in range(k0, k1):
            rr[k] = hypot(hxy, c[k])

    _pieces(fill_r, n * n * n, nt)
    shells_r = mround((rr + dr / 2) / dr)
    rr = None
    shells_a = mround(az / dtheta) + 1
    shells_a[shells_a == shells_a.max()] = shells_a.min()
    outer = mround(((n - 1) / 2 + 1) / dr) + 1
    shells_r[(shells_r - 1) * dr > (n - 1) / 2 + 1] = outer
    max_a = shells_a.max()
    sa = np.broadcast_to(shells_a[:, :, None], (n, n, n)).copy()
    sa[shells_r == shells_r.max()] = max_a + 1
    max_r = int(shells_r.max())
    max_e = int(sa.max())
    shells = (shells_r - 1) + (sa - 1) * max_r
    del sa, shells_r
    if skip_cone_region:
        cone = (az < (np.pi / 36)) | (az > np.pi - (np.pi / 36))
        cone |= np.abs(az - np.pi / 2) < np.pi / 48
        shells[np.broadcast_to(cone[:, :, None], (n, n, n))] = max_r * max_e - 1
    need = shells != max_r * max_e - 1  # the voxels not in the last cell, (z, x, y)
    shells = shells.astype(np.intp).ravel()
    counts = np.bincount(shells, minlength=max_r * max_e)
    # the voxels of each cell, in MATLAB's linear order (a stable sort), and where each starts
    order = np.argsort(shells, kind="stable")
    starts = np.concatenate([[0], np.cumsum(counts)])
    got = (order, starts, counts.astype(np.float64), max_r, max_e, need)
    with _LOCK:
        _SHELLS.clear()
        _SHELLS[key] = got
    return got


def _cell_sums(
    values: np.ndarray, order: np.ndarray, starts: np.ndarray, nt: int, last: bool = True
) -> np.ndarray:
    """``accumarray(shells(:), values(:))``: each cell's sum from +0 in MATLAB's linear order
    (``values`` in that order), as ``np.bincount`` adds; the values gathered by cell and each
    cell's run added in sequence (``np.add.accumulate``), on threads. ``last`` False: the last
    cell is left out (0)."""
    nc = starts.size - 1
    v = np.empty(int(starts[nc if last else nc - 1]))

    def gather(lo, hi):
        np.take(values, order[lo:hi], out=v[lo:hi])

    _pieces(gather, v.size, nt)
    out = np.zeros(nc)
    nc = nc if last else nc - 1

    def cells(lo, hi):
        for c in range(lo, hi):
            a, b = int(starts[c]), int(starts[c + 1])
            if b > a:
                # 0 + v1 + v2 + ...: + 0.0 turns a sum of -0s into +0, as bincount's
                out[c] = np.add.accumulate(v[a:b])[-1] + 0.0

    _pieces(cells, nc, nt, grain=max(1, nc // (4 * max(1, nt))))
    return out


def _frc(f1: np.ndarray, f2: np.ndarray, cells: tuple, nt: int, fast: bool = True) -> np.ndarray:
    """``real(accumarray(shells, F1 .* conj(F2))) ./ sqrt(accumarray(shells, |F1|^2) .*
    accumarray(shells, |F2|^2))`` of two FFTs (z, y, x), fftshifted, per cell.

    The last cell (the outer shell and the cones) is dropped by XR_FSC unless its value is NaN,
    so only that matters: it is NaN when F1 or F2 is zero all over it (0 / 0). ``fast`` computes
    it so, without the cell's sums, unless the largest value there is so small or large that
    the product of the sums could underflow or overflow, or is not finite."""
    order, starts, _, _, _, need = cells
    if fast:
        cross, p1, p2, peaks = _spectrum_values(f1, f2, nt, need)
        top = peaks.max(axis=0) if peaks.size else np.zeros(2)
        sure = bool(np.all(np.isfinite(peaks))) and all(
            v == 0 or 2.0**-250 <= v < 2.0**400 for v in top
        )
        if sure:
            sums = [_cell_sums(v, order, starts, nt, last=False) for v in (cross, p1, p2)]
            empty = 0.0 if (top == 0).any() else 1.0
            for s in sums:
                s[-1] = empty  # 0 / 0 is NaN; 1 / 1 is a value
            with np.errstate(invalid="ignore", divide="ignore"):
                return sums[0] / np.sqrt(sums[1] * sums[2])
    cross, p1, p2 = _spectrum_values(f1, f2, nt)
    sc, s1, s2 = (_cell_sums(v, order, starts, nt) for v in (cross, p1, p2))
    with np.errstate(invalid="ignore", divide="ignore"):
        return sc / np.sqrt(s1 * s2)


def _times(a: np.ndarray, w: np.ndarray, nt: int) -> np.ndarray:
    """``a .* w`` of two double volumes (``a`` may be a strided view), plane by plane"""
    out = np.empty(w.shape)

    def planes(lo, hi):
        np.multiply(a[lo:hi], w[lo:hi], out=out[lo:hi])

    _pieces(planes, w.shape[0], nt, grain=1)
    return out


def clear_caches() -> None:
    """forget the kept windows and shells (MATLAB's ``clear XR_FSC``)"""
    with _LOCK:
        _WINDOW.clear()
        _SHELLS.clear()


# ---------------------------------------------------------------------------------------------
# XR_FSC
# ---------------------------------------------------------------------------------------------


def _half3(v: Any, name: str) -> tuple[int, ...]:
    a = np.asarray(v, dtype=np.float64).ravel()
    if a.size != 3:
        raise FSCError(f"{name}: three values (z, y, x)")
    return tuple(float(x) for x in a)


def _fit_size(im: np.ndarray, half: tuple) -> np.ndarray:
    """XR_FSC's pad (floor before, ceil after) and crop about round((size + 1) / 2)"""
    from petakit5d.core.dtypes import mround

    sz = im.shape
    if all(s == h for s, h in zip(sz, half, strict=True)):
        return im
    pre = [max(0, math.floor((h - s) / 2)) for s, h in zip(sz, half, strict=True)]
    post = [max(0, math.ceil((h - s) / 2)) for s, h in zip(sz, half, strict=True)]
    im = np.pad(im, list(zip(pre, post, strict=True)))
    sl = []
    for s, h in zip(im.shape, half, strict=True):
        cc = float(mround((s + 1) / 2))
        hf = (h - 1) / 2
        lo, hi = cc - hf, cc + hf
        if lo != math.floor(lo) or lo < 1 or hi > s:
            raise FSCError(
                "Array indices must be positive integers or logical values.",
                "MATLAB:badsubscript",
            )
        sl.append(slice(int(lo) - 1, int(hi)))
    return im[tuple(sl)]


def _as_double(im: Any, name: str) -> np.ndarray:
    a = np.asarray(im)
    if a.dtype != np.float64:
        raise NotImplementedError(
            f"{name}: only double volumes are ported (MATLAB's XR_FSC of other classes is not)"
        )
    if a.ndim != 3:
        raise FSCError(f"{name}: a 3D volume (z, y, x)")
    return a


def _spectrum_values(
    f1: np.ndarray, f2: np.ndarray, nt: int, need: np.ndarray | None = None
) -> tuple[np.ndarray, ...]:
    """``real(F1 .* conj(F2))``, ``abs(F1) .^ 2`` and ``abs(F2) .^ 2`` of the fftshifted spectra
    of two (z, y, x) FFTs, flat in MATLAB's linear order (the C order of (z, x, y)); plane by
    plane. The FFTs of real volumes are Hermitian (MATLAB's fill, bit for bit), so the values at
    -k are those at k (|conj(c)| = |c|, and the products' real parts are equal): the plane of
    -z is gathered from that of z."""
    nz, ny, nx = f1.shape
    src = [(np.arange(m) - m // 2) % m for m in (nz, ny, nx)]  # shifted index -> unshifted
    mir = [(2 * (m // 2) - np.arange(m)) % m for m in (nz, ny, nx)]  # -k, shifted
    out = [np.empty((nz, nx, ny)) for _ in range(3)]
    todo = [q for q in range(nz) if q <= mir[0][q]]
    gat = np.ix_(src[1], src[2])
    back = np.ix_(mir[2], mir[1])
    # with need (the voxels outside the last cell, (z, x, y)): |F|^2 only where a plane or its
    # mirror needs it, and per plane the largest |part| of F1 and F2 in the last cell
    peaks = np.zeros((len(todo), 2))

    def rest_flags(t, sel, parts):
        for col, (r, i) in enumerate(parts):
            a = np.maximum(np.abs(r[sel]), np.abs(i[sel]))
            peaks[t, col] = a.max() if a.size else 0.0

    def planes(lo, hi):
        cross, p1, p2 = out
        for t in range(lo, hi):
            q = todo[t]
            m = mir[0][q]
            g1 = f1[src[0][q]][gat].T
            g2 = f2[src[0][q]][gat].T
            r1, i1 = np.ascontiguousarray(g1.real), np.ascontiguousarray(g1.imag)
            r2, i2 = np.ascontiguousarray(g2.real), np.ascontiguousarray(g2.imag)
            # real(F1 .* conj(F2)): r1 r2 - i1 (-i2), no fused multiply-add
            cross[q] = r1 * r2 + i1 * i2
            if need is None:
                p1[q] = _abs_squared(r1, i1)
                p2[q] = _abs_squared(r2, i2)
            else:
                nq, nm = need[q], need[m][back]
                mk = nq | nm
                for o, r, i in ((p1, r1, i1), (p2, r2, i2)):
                    v = np.zeros((nx, ny))
                    v[mk] = _abs_squared(r[mk], i[mk])
                    o[q] = v
                rest_flags(t, ~nq | ~nm, ((r1, i1), (r2, i2)))
            if m != q:
                for o in out:
                    o[m] = o[q][back]

    _pieces(planes, len(todo), nt, grain=1)
    vals = tuple(o.reshape(-1) for o in out)
    return vals if need is None else (*vals, peaks)


@register(
    "XR_FSC",
    Param("im1", "@isnumeric", "First volume (z, y, x), double.", convert="volume", py="im1"),
    Param("im2", "@isnumeric", "Second volume (z, y, x), double.", convert="volume", py="im2"),
    Param("xyPixelSize", "@isnumeric", "Pixel size in x and y (um)."),
    Param("dz", "@isnumeric", "Voxel size in z (um)."),
    Param("dr", "@isnumeric", "Radial interval (voxels of the spectrum)."),
    Param("dtheta", "@isnumeric", "Angular interval (rad).", matlab_default="pi / 12"),
    Param("resThreshMethod", "@ischar", "'fixed', 'one-bit' or 'half-bit'."),
    Param("resThresh", "@isnumeric", "Threshold of 'fixed'."),
    Param("skipConeRegion", "@islogical", "Leave out the cones along z and x."),
    Param(
        "halfSize", "@isnumeric", "Size (z, y, x) of the volumes (all equal).", convert="size_yxz"
    ),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    task="B05",
)
def fourier_shell_correlation(
    im1: np.ndarray,
    im2: np.ndarray,
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.108,
    dr: float = 1,
    dtheta: float = math.pi / 12,
    res_thresh_method: str = "fixed",
    res_thresh: float = 0.2,
    skip_cone_region: bool = True,
    half_size: Sequence[float] = (501, 501, 501),
    debug: bool = False,
    num_threads: int | None = None,
) -> tuple[dict, dict]:
    """The FSC of two double volumes (z, y, x) per radius and angle, and its resolution
    (``XR_FSC``; see the module): ``(fsc, res)``, dicts of MATLAB's struct fields (vectors as
    1-D arrays, ``fsc['fsc']`` and ``fsc['npoints']`` radius x angle)."""
    from petakit5d.matlabcompat.fft import fftn

    nt = _threads(num_threads)
    half = _half3(half_size, "half_size")
    a = _fit_size(_as_double(im1, "im1"), half)
    b = _fit_size(_as_double(im2, "im2"), half)
    n = int(half[1])  # MATLAB's halfSize(1): y
    w = _window(n)
    if a.shape != w.shape or b.shape != w.shape:
        raise FSCError(
            "Arrays have incompatible sizes for this operation.", "MATLAB:sizeDimensionsMustMatch"
        )
    cells = _shells(n, float(dr), float(dtheta), bool(skip_cone_region), nt)
    npoints, max_r, max_e = cells[2:5]
    num_xi = max_r * max_e
    f1 = fftn(_times(a, w, nt), threads=nt)
    f2 = fftn(_times(b, w, nt), threads=nt)
    frc = _frc(f1, f2, cells, nt)
    del f1, f2
    keep = ~np.isnan(frc)
    if not keep.any():
        # (every cell is 0 / 0: no FSC)
        raise FSCError(
            "The FSC is undefined: the images have no signal (e.g., they are all zero).",
            "XR_FSC:noSignal",
        )
    xi = np.arange(1, num_xi + 1)[keep]
    xi, n_points, frc = xi[:-1], npoints[keep][:-1], frc[keep][:-1]
    xri = (xi - 1) % max_r + 1
    xei = (xi - 1) // max_r + 1
    unique_e = np.unique(xei)
    npoints_mat = np.zeros((max_r, max_e))
    frc_mat = np.zeros((max_r, max_e))
    for i, e in enumerate(unique_e):
        inds = xei == e
        npoints_mat[xri[inds] - 1, i] = n_points[inds]
        frc_mat[xri[inds] - 1, i] = frc[inds]
    npoints_mat = npoints_mat[:-1, :-1]
    frc_mat = frc_mat[:-1, :-1].copy()
    if frc_mat.shape[0] and frc_mat.shape[1]:
        row = frc_mat[0]
        row[row == 0] = 1
    fsc = {
        "r": (np.arange(max_r - 1, dtype=np.float64) * dr) / ((n - 1) / 2),
        "thetas": (unique_e.astype(np.float64) - 1) * dtheta,
        "fsc": frc_mat,
        "npoints": np.ascontiguousarray(npoints_mat),
    }
    res = fsc_resolution(
        fsc,
        xy_pixel_size=xy_pixel_size,
        dz=dz,
        res_thresh_method=res_thresh_method,
        res_thresh=res_thresh,
        dtheta=dtheta,
        debug=debug,
    )
    return fsc, res


# ---------------------------------------------------------------------------------------------
# smoothing and the resolution (XR_FSC_resolution)
# ---------------------------------------------------------------------------------------------


class _Pieces:
    """a piecewise polynomial in local coordinates (MATLAB's ppform): ``coefs[i]`` (highest power
    first) on ``[breaks[i], breaks[i + 1])``, the end pieces extended beyond the breaks"""

    def __init__(self, breaks: np.ndarray, coefs: np.ndarray):
        self.breaks = np.asarray(breaks, dtype=np.float64)
        self.coefs = np.asarray(coefs, dtype=np.float64)

    def __call__(self, x) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        nb = self.breaks.size
        if nb < 2:
            return np.full(x.shape, self.coefs[0, -1])
        i = np.clip(np.searchsorted(self.breaks, x, side="right") - 1, 0, nb - 2)
        t = x - self.breaks[i]
        c = self.coefs[i]
        v = c[..., 0]
        for k in range(1, self.coefs.shape[1]):
            v = v * t + c[..., k]
        return v


def default_smoothing(x: Sequence[float]) -> float:
    """The documented default smoothing parameter of ``csaps`` (and so of ``fit``'s
    'smoothingSpline', which the documentation says is based on csaps): the ``p`` that makes
    ``p trace(A)`` equal ``(1 - p) trace(B)``, where ``p A + (1 - p) B`` is the coefficient matrix
    of the spline's linear system (de Boor's SMOOTH, unit weights). Here that matrix is
    ``p R + (1 - p) Q'Q`` (the system below times ``p``), so ``A = R`` and ``B = Q'Q``:
    ``p = trace(Q'Q) / (trace(R) + trace(Q'Q))``. 1 for two points."""
    xs = [float(v) for v in x]
    h = [xs[i + 1] - xs[i] for i in range(len(xs) - 1)]
    ta = tb = 0.0
    for j in range(len(h) - 1):
        ta += (h[j] + h[j + 1]) / 3.0
        q0, q2 = 1.0 / h[j], 1.0 / h[j + 1]
        q1 = -q0 - q2
        tb += q0 * q0 + q1 * q1 + q2 * q2
    if ta + tb == 0:
        return 1.0
    return tb / (ta + tb)


def smoothing_spline(x: Sequence[float], y: Sequence[float], p: float | None = None) -> _Pieces:
    """The cubic smoothing spline of the points ``(x, y)`` (``x`` increasing), minimizing
    ``p sum (y_i - f(x_i))^2 + (1 - p) int f''(t)^2 dt``; ``p`` None: ``default_smoothing(x)``,
    csaps's documented default. Reinsch's algorithm (Green and Silverman 1994, 2.3.3): ``(R + a Q'Q)
    gamma = Q'y``, ``g = y - a Q gamma`` with ``a = (1 - p) / p``, solved by a banded LDL'
    factorization in Python floats; returned as a callable piecewise cubic."""
    xs = [float(v) for v in x]
    ys = [float(v) for v in y]
    n = len(xs)
    if n != len(ys) or n < 2:
        raise FSCError("Not enough data for a smoothing spline (two points or more).")
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    if any(not v > 0 for v in h):
        raise FSCError("The data sites must be distinct and increasing.")
    if p is None:
        p = default_smoothing(xs)
    p = float(p)
    if not 0 < p <= 1:
        raise FSCError("The smoothing parameter must be in (0, 1].")
    al = (1.0 - p) / p
    m = n - 2
    gam = [0.0] * n
    g = list(ys)
    if m > 0:
        q0 = [1.0 / h[j] for j in range(m)]
        q1 = [-1.0 / h[j] - 1.0 / h[j + 1] for j in range(m)]
        q2 = [1.0 / h[j + 1] for j in range(m)]
        d0 = [
            (h[j] + h[j + 1]) / 3.0 + al * (q0[j] * q0[j] + q1[j] * q1[j] + q2[j] * q2[j])
            for j in range(m)
        ]
        d1 = [h[j + 1] / 6.0 + al * (q1[j] * q0[j + 1] + q2[j] * q1[j + 1]) for j in range(m - 1)]
        d2 = [al * (q2[j] * q0[j + 2]) for j in range(m - 2)]
        rhs = [(ys[j + 2] - ys[j + 1]) / h[j + 1] - (ys[j + 1] - ys[j]) / h[j] for j in range(m)]
        # LDL' of the symmetric pentadiagonal matrix (d0, d1, d2)
        dd = [0.0] * m
        l1 = [0.0] * m
        l2 = [0.0] * m
        for j in range(m):
            v = d0[j]
            if j >= 1:
                v -= l1[j - 1] * l1[j - 1] * dd[j - 1]
            if j >= 2:
                v -= l2[j - 2] * l2[j - 2] * dd[j - 2]
            dd[j] = v
            if j + 1 < m:
                u = d1[j]
                if j >= 1:
                    u -= l1[j - 1] * l2[j - 1] * dd[j - 1]
                l1[j] = u / v
            if j + 2 < m:
                l2[j] = d2[j] / v
        z = [0.0] * m
        for j in range(m):
            v = rhs[j]
            if j >= 1:
                v -= l1[j - 1] * z[j - 1]
            if j >= 2:
                v -= l2[j - 2] * z[j - 2]
            z[j] = v
        c = [z[j] / dd[j] for j in range(m)]
        for j in range(m - 1, -1, -1):
            v = c[j]
            if j + 1 < m:
                v -= l1[j] * c[j + 1]
            if j + 2 < m:
                v -= l2[j] * c[j + 2]
            c[j] = v
        for j in range(m):
            gam[j + 1] = c[j]
        for i in range(n):
            qg = 0.0
            if i < m:
                qg += q0[i] * c[i]
            if 1 <= i <= m:
                qg += q1[i - 1] * c[i - 1]
            if 2 <= i <= m + 1:
                qg += q2[i - 2] * c[i - 2]
            g[i] = ys[i] - al * qg
    coefs = np.empty((n - 1, 4))
    for i in range(n - 1):
        hi = h[i]
        coefs[i, 0] = (gam[i + 1] - gam[i]) / (6.0 * hi)
        coefs[i, 1] = gam[i] / 2.0
        coefs[i, 2] = (g[i + 1] - g[i]) / hi - hi * (2.0 * gam[i] + gam[i + 1]) / 6.0
        coefs[i, 3] = g[i]
    return _Pieces(np.asarray(xs), coefs)


def poly2_fit(x: Sequence[float], y: Sequence[float]):
    """The least-squares quadratic ``p1 x^2 + p2 x + p3`` of the points, from the normal
    equations in exact rational arithmetic, rounded to doubles; returned as a callable
    (``(p1 x^2 + p2 x) + p3`` in double)"""
    fx = [Fraction(float(v)) for v in x]
    fy = [Fraction(float(v)) for v in y]
    s = [sum(v**k for v in fx) for k in range(5)]
    t = [sum((v**k) * w for v, w in zip(fx, fy, strict=True)) for k in range(3)]
    a = [[s[4], s[3], s[2]], [s[3], s[2], s[1]], [s[2], s[1], s[0]]]
    rhs = [t[2], t[1], t[0]]

    def det(mm):
        return (
            mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
            - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
            + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0])
        )

    d = det(a)
    if d == 0:
        raise FSCError("A quadratic needs three distinct points.")
    coef = []
    for k in range(3):
        mk = [[rhs[i] if j == k else a[i][j] for j in range(3)] for i in range(3)]
        coef.append(float(det(mk) / d))
    p1, p2, p3 = coef

    def f(xv):
        xv = np.asarray(xv, dtype=np.float64)
        return (p1 * (xv * xv) + p2 * xv) + p3

    f.coefficients = (p1, p2, p3)  # type: ignore[attr-defined]
    return f


def _first(mask: np.ndarray) -> int | None:
    hit = np.flatnonzero(mask)
    return int(hit[0]) if hit.size else None


@register(
    "XR_FSC_resolution",
    Param("fsc", "@isstruct", "The FSC (fields r, thetas, fsc, npoints)."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in x and y (um)."),
    Param("dz", "@isnumeric", "Voxel size in z (um)."),
    Param("resThreshMethod", "@ischar", "'fixed', 'one-bit' or 'half-bit'."),
    Param("resThresh", "@isnumeric", "Threshold of 'fixed'."),
    Param(
        "dtheta",
        "@isnumeric",
        "Interval of the angles (rad; None: that of the first two angles).",
    ),
    Param("debug", "@islogical", "Not used."),
    task="B05",
)
def fsc_resolution(
    fsc: dict,
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.108,
    res_thresh_method: str = "fixed",
    res_thresh: float = 0.2,
    dtheta: float | None = None,
    debug: bool = False,
) -> dict:
    """The resolution per angle of an FSC (``XR_FSC_resolution``; see the module): a dict of
    ``comp_angles`` (the FSC's angles), ``points`` (the radius found per angle), ``thetas``
    (``0 : dtheta : 2 pi``) and ``res`` (the resolution on that circle: 0 for x, pi/2 for z;
    each angle's sector at its place, NaN for the sectors the FSC lacks). The smoothing spline
    is not MATLAB's (the module)."""
    from petakit5d.matlabcompat.elementary import cos

    frc_mat = np.asarray(fsc["fsc"], dtype=np.float64)
    frc_mat = frc_mat.reshape(frc_mat.shape[0], -1) if frc_mat.ndim != 2 else frc_mat
    r = np.asarray(fsc["r"], dtype=np.float64).ravel()
    thetas = np.asarray(fsc["thetas"], dtype=np.float64).ravel()
    npoints = np.asarray(fsc["npoints"], dtype=np.float64)
    nth = thetas.size
    sqrt2 = math.sqrt(2)
    if res_thresh_method == "fixed":
        thr = float(res_thresh)
        points = np.zeros(nth)
        for i in range(nth):
            frc_i = frc_mat[:, i]
            f = smoothing_spline(r, frc_i)
            rs = _colon(r[0], 0.0001, 1)
            ind = _first(f(rs) < thr)
            if ind is None:
                rs = _colon(r[0], 0.0001, sqrt2)
                ind = _first(f(rs) < thr)
                if ind is None:
                    f = poly2_fit(r, frc_i)
                    ind = _first(f(rs) < thr)
            points[i] = sqrt2 if ind is None else rs[ind]
    elif res_thresh_method in ("one-bit", "half-bit"):
        points = np.zeros(nth)
        for i in range(nth):
            frc_i = frc_mat[:, i]
            f = smoothing_spline(r, frc_i)
            sn = np.sqrt(npoints[:, i])
            if res_thresh_method == "one-bit":
                tp = ((sn + 2 * sqrt2) + 2) / (3 * sn + 2 * sqrt2)
                lim = 1 / 3
            else:
                tp = (sn * 0.2071 + 1.9102) / (sn * 1.2071 + 0.9102)
                lim = 0.2071 / 1.2071
            ts = smoothing_spline(np.append(r[1:], 2.0), np.append(tp[1:], lim))
            rs = _colon(r[0], 0.0001, 1)
            ind = _first(f(rs) - ts(rs) <= 0)
            if ind is None:
                rs = _colon(r[0], 0.0001, sqrt2)
                f = poly2_fit(r, frc_i)
                ind = _first(f(rs) - ts(rs) <= 0)
            points[i] = sqrt2 if ind is None else rs[ind]
    else:
        raise FSCError("Unrecognized function or variable 'points'.", "MATLAB:UndefinedFunction")

    psz = xy_pixel_size + (dz - xy_pixel_size) * np.abs(cos(thetas))
    res_mat = (2 * psz) / points
    if dtheta is None or np.size(dtheta) == 0:
        if nth < 2:
            raise FSCError(
                "Index exceeds the number of array elements. Index must not exceed 0.",
                "MATLAB:badsubscript",
            )
        dth = thetas[1] - thetas[0]
    else:
        dth = float(np.asarray(dtheta, dtype=np.float64).ravel()[0])
    res_theta = _colon(0, dth, 2 * math.pi)
    # the full circle (0 for x, pi/2 for z): the angle phi takes the sector of thetas at
    # phi + pi/2 (mod pi), each angle's sector at its place (round(theta / dtheta)), the ones
    # the FSC lacks NaN; an odd number of sectors has none at pi/2
    from petakit5d.core.dtypes import NAN, mround

    n_sector = int(mround(math.pi / dth))
    if n_sector % 2 != 0:
        raise FSCError(
            f"The half circle has {n_sector} sectors (pi / dtheta): an even number is needed, "
            "so that x (pi / 2) is one of them.",
            "XR_FSC_resolution:oddSectorNumber",
        )
    idx = mround(thetas / dth).astype(np.int64)
    res_s = np.full(max(n_sector, int(idx.max(initial=-1)) + 1), NAN)
    res_s[n_sector:] = 0  # (MATLAB grows the array with zeros)
    res_s[idx] = res_mat
    cind = n_sector // 2 + 1
    res_out = np.concatenate([res_s[cind - 1 :], res_s, res_s[:cind]])
    return {"comp_angles": thetas, "points": points, "thetas": res_theta, "res": res_out}


# ---------------------------------------------------------------------------------------------
# one image: XR_one_image_FSC_analysis
# ---------------------------------------------------------------------------------------------


def _read_image(fn: str) -> np.ndarray:
    from petakit5d.io import read_tiff

    return read_tiff(fn)


def _path_or_volume_in(v: Any) -> Any:
    from petakit5d.core.arrays import from_matlab

    return v if isinstance(v, (str, os.PathLike)) else from_matlab(np.asarray(v))


def _path_or_volume_out(v: Any) -> Any:
    from petakit5d.core.arrays import to_matlab

    return v if isinstance(v, (str, os.PathLike)) else to_matlab(np.asarray(v))


def _clip(im: np.ndarray, clip_per: float) -> np.ndarray:
    """``im .* (im <= v) + v .* (im > v)`` with ``v = prctile(im(im > 0), clip_per)``, in double"""
    from petakit5d.matlabcompat.stats import prctile

    im = np.asarray(im).astype(np.float64)
    v = float(np.asarray(prctile(im[im > 0], float(clip_per))).ravel()[0])
    with np.errstate(invalid="ignore"):
        return im * (im <= v) + v * (im > v)


def _isotropic_region(
    im: np.ndarray, xy: float, dz: float, half: tuple, scale: tuple[float, float, float]
) -> np.ndarray:
    """the image resized to isotropic voxels (linear, in double), only the region kept, then
    padded to twice ``half``; ``scale``: the pixel sizes of its axes (z, y, x)"""
    from petakit5d.core.dtypes import mround
    from petakit5d.matlabcompat.resize import imresize3_fast

    sz = im.shape  # (z, y, x)
    psz = min(xy, dz)
    if dz != xy:
        rsz = tuple(int(mround((s * f) / psz)) for s, f in zip(sz, scale, strict=True))
    else:
        rsz = tuple(sz)
    s = [1, 1, 1]
    t = list(rsz)
    if any(rv > hv * 2 for rv, hv in zip(rsz, half, strict=True)):
        for k in range(3):
            hs = math.floor((rsz[k] - half[k] * 2) / 2)
            s[k] = max(1, hs)
            t[k] = int(min(rsz[k], s[k] + half[k] * 2 - 1))
    if dz != xy:
        box = (s[0] - 1, s[1] - 1, s[2] - 1, t[0], t[1], t[2])
        im = imresize3_fast(im, rsz, "linear", box, "double")
    else:
        im = np.asarray(im[s[0] - 1 : t[0], s[1] - 1 : t[1], s[2] - 1 : t[2]]).astype(np.float64)
    sz = im.shape
    if any(sv < hv * 2 for sv, hv in zip(sz, half, strict=True)):
        lsz = [math.floor((half[k] * 2 - sz[k]) / 2) for k in range(3)]
        rpad = [int(half[k] * 2 - sz[k] - lsz[k]) for k in range(3)]
        im = np.pad(im, [(int(lv), rv) for lv, rv in zip(lsz, rpad, strict=True)])
    return np.ascontiguousarray(im)


# sub-lattice pairs: MATLAB's (y, x, z) parities of im(:, :, :, a, b, c)
_PAIRS = (
    ((1, 1, 1), (2, 2, 2)),
    ((2, 1, 1), (1, 2, 2)),
    ((1, 2, 1), (2, 1, 2)),
    ((1, 1, 2), (2, 2, 1)),
)


def _sub_lattice(im: np.ndarray, parity: tuple[int, int, int]) -> np.ndarray:
    """``im(a : 2 : end, b : 2 : end, c : 2 : end)`` for MATLAB's (y, x, z) parities, a view"""
    a, b, c = parity
    return im[c - 1 :: 2, a - 1 :: 2, b - 1 :: 2]


@register(
    "XR_one_image_FSC_analysis",
    Param(
        "fn",
        "@(x) ischar(x) || isnumeric(x)",
        "Tiff file, or the image (z, y, x).",
        convert=(_path_or_volume_in, _path_or_volume_out),
    ),
    Param("xyPixelSize", "@isnumeric", "Pixel size in x and y (um)."),
    Param("dz", "@isnumeric", "Voxel size in z (um)."),
    Param("dr", "@isnumeric", "Radial interval (voxels of the spectrum)."),
    Param("dtheta", "@isnumeric", "Angular interval (rad).", matlab_default="pi / 12"),
    Param("resThreshMethod", "@ischar", "'fixed', 'one-bit' or 'half-bit'."),
    Param("resThresh", "@isnumeric", "Threshold of 'fixed'."),
    Param(
        "halfSize",
        "@isnumeric",
        "Half size (z, y, x) of the region (all equal).",
        convert="size_yxz",
    ),
    Param(
        "bbox",
        "@isnumeric",
        "Box (z0, y0, x0, z1, y1, x1) of the region (None: all).",
        convert="bbox",
    ),
    Param("resAxis", "@ischar", "'xz', 'yz' (x and y swapped) or 'xy'.", convert="same"),
    Param("skipConeRegion", "@islogical", "Leave out the cones along z and x."),
    Param("clipPer", "@isnumeric", "Clip above this percentile of the positive voxels (None: no)."),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    task="B05",
)
def one_image_fsc_analysis(
    fn: str | np.ndarray,
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.108,
    dr: float = 1,
    dtheta: float = math.pi / 12,
    res_thresh_method: str = "fixed",
    res_thresh: float = 0.2,
    half_size: Sequence[float] = (501, 501, 501),
    bbox: Sequence[int] | None = None,
    res_axis: str = "xz",
    skip_cone_region: bool = True,
    clip_per: float | None = None,
    debug: bool = False,
    num_threads: int | None = None,
) -> tuple[dict, dict, list, list]:
    """The FSC of an image's interleaved sub-lattices (``XR_one_image_FSC_analysis``; see the
    module): ``(fsc_mu, res_mu, fsc_cell, res_cell)``, the mean FSC of the four pairs and its
    resolution, and each pair's FSC and resolution."""
    half = _half3(half_size, "half_size")
    im = _read_image(os.fspath(fn)) if isinstance(fn, (str, os.PathLike)) else np.asarray(fn)
    if im.ndim == 2:
        im = im[None]
    if clip_per is not None and np.size(clip_per):
        im = _clip(im, float(np.asarray(clip_per).ravel()[0]))
    if bbox is not None and len(bbox):
        z0, y0, x0, z1, y1, x1 = (int(v) for v in bbox)
        if z1 > im.shape[0] or y1 > im.shape[1] or x1 > im.shape[2] or min(z0, y0, x0) < 0:
            raise FSCError("Index in position 1 exceeds array bounds.", "MATLAB:badsubscript")
        im = im[z0:z1, y0:y1, x0:x1]
    xy, dz_ = float(xy_pixel_size), float(dz)
    # the pixel sizes of the axes (z, y, x) after the permutation, and the FSC's dz (of its
    # plane's second axis: for 'xy', y)
    scale, fsc_dz = (dz_, xy, xy), dz
    if res_axis == "yz":
        im = np.swapaxes(im, 1, 2)
    elif res_axis == "xy":
        im = np.transpose(im, (1, 0, 2))
        scale, fsc_dz = (xy, dz_, xy), xy_pixel_size
    im = _isotropic_region(np.ascontiguousarray(im), xy, dz_, half, scale)
    kw = {
        "xy_pixel_size": xy_pixel_size,
        "dz": fsc_dz,
        "dr": dr,
        "dtheta": dtheta,
        "half_size": half_size,
        "skip_cone_region": skip_cone_region,
        "num_threads": num_threads,
    }
    fsc_cell, res_cell = [], []
    for pa, pb in _PAIRS:
        f, r = fourier_shell_correlation(_sub_lattice(im, pa), _sub_lattice(im, pb), **kw)
        fsc_cell.append(f)
        res_cell.append(r)
    fsc_mu = dict(fsc_cell[0])
    fsc_mu["fsc"] = (
        ((fsc_cell[0]["fsc"] + fsc_cell[1]["fsc"]) + fsc_cell[2]["fsc"]) + fsc_cell[3]["fsc"]
    ) / 4
    res_mu = fsc_resolution(
        fsc_mu,
        xy_pixel_size=xy_pixel_size,
        dz=fsc_dz,
        res_thresh_method=res_thresh_method,
        res_thresh=res_thresh,
        dtheta=dtheta,
        debug=debug,
    )
    return fsc_mu, res_mu, fsc_cell, res_cell


# ---------------------------------------------------------------------------------------------
# a file: XR_one_image_FSC_analysis_frame
# ---------------------------------------------------------------------------------------------


def _col(v: Any) -> np.ndarray:
    a = np.asarray(v, dtype=np.float64).ravel()
    return a.reshape(-1, 1)


def _fsc_struct(f: dict) -> dict:
    """an FSC dict in MATLAB's form (columns)"""
    return {
        "r": _col(f["r"]),
        "thetas": _col(f["thetas"]),
        "fsc": np.asarray(f["fsc"], dtype=np.float64),
        "npoints": np.asarray(f["npoints"], dtype=np.float64),
    }


def _res_struct(r: dict) -> dict:
    return {k: _col(r[k]) for k in ("comp_angles", "points", "thetas", "res")}


def _column_cell(items: list) -> np.ndarray:
    out = np.empty((len(items), 1), dtype=object)
    for i, v in enumerate(items):
        out[i, 0] = v
    return out


def _boxes(bbox: Any) -> list[tuple[int, ...]]:
    if bbox is None:
        return []
    a = np.asarray(bbox, dtype=np.int64)
    if a.size == 0:
        return []
    a = a.reshape(-1, 6)
    return [tuple(int(v) for v in row) for row in a]


@register(
    "XR_one_image_FSC_analysis_frame",
    Param("fn", "@(x) ischar(x) || isnumeric(x)", "Tiff file of the image."),
    Param("fnout", "@ischar", "Result file (.mat)."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in x and y (um)."),
    Param("dz", "@isnumeric", "Voxel size in z (um)."),
    Param("dr", "@isnumeric", "Radial interval (voxels of the spectrum)."),
    Param("dtheta", "@isnumeric", "Angular interval (rad).", matlab_default="pi / 12"),
    Param("resThreshMethod", "@ischar", "'fixed', 'one-bit' or 'half-bit'."),
    Param("resThresh", "@isnumeric", "Threshold of 'fixed'."),
    Param(
        "halfSize",
        "@isnumeric",
        "Half size (z, y, x) of the regions (all equal).",
        convert="size_yxz",
    ),
    Param(
        "bbox",
        "@isnumeric",
        "Boxes (z0, y0, x0, z1, y1, x1), one region each (None: all).",
        convert="bboxes",
    ),
    Param("resAxis", "@ischar", "'xz', 'yz' (x and y swapped) or 'xy'.", convert="same"),
    Param("skipConeRegion", "@islogical", "Leave out the cones along z and x."),
    Param("clipPer", "@isnumeric", "Clip above this percentile of the positive voxels (None: no)."),
    Param("debug", "@islogical", "Not used."),
    task="B05",
)
def fsc_analysis_file(
    fn: str,
    fnout: str,
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.108,
    dr: float = 1,
    dtheta: float = math.pi / 12,
    res_thresh_method: str = "fixed",
    res_thresh: float = 0.2,
    half_size: Sequence[float] = (501, 501, 501),
    bbox: Sequence[Sequence[int]] | None = None,
    res_axis: str = "xz",
    skip_cone_region: bool = True,
    clip_per: float | None = None,
    debug: bool = False,
):
    """The FSC analysis of a tiff file (of each box of ``bbox``) saved as ``fnout``
    (``XR_one_image_FSC_analysis_frame``): v7.3 with ``fsc_mu``, ``res_mu``, ``fsc_cell``,
    ``res_cell`` (n x 1 cells for several boxes), ``centers`` (MATLAB's ``[y x z]``: the image
    size, or the boxes' centres) and ``bbox`` (MATLAB's boxes), written under a temporary name
    and renamed. Nothing if ``fnout`` exists. Returns what MATLAB returns (None then)."""
    from petakit5d.core.arrays import size_to_matlab
    from petakit5d.core.bbox import box_to_matlab
    from petakit5d.core.dtypes import mround
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import write_mat
    from petakit5d.io.imagesize import get_image_size
    from petakit5d.jobs.files import exist

    if exist(fnout):
        print(f"Result file {fnout} already exists!", flush=True)
        return None
    boxes = _boxes(bbox)
    mboxes = np.asarray([box_to_matlab(b) for b in boxes], dtype=np.float64).reshape(-1, 6)
    if boxes:
        centers = mround((mboxes[:, 3:6] + mboxes[:, 0:3]) / 2)
    else:
        centers = np.asarray(size_to_matlab(get_image_size(fn)), dtype=np.float64)[None, :]
    im = _read_image(fn) if isinstance(fn, (str, os.PathLike)) else np.asarray(fn)
    out = []
    for b in boxes or [None]:
        out.append(
            one_image_fsc_analysis(
                im,
                xy_pixel_size=xy_pixel_size,
                dz=dz,
                dr=dr,
                dtheta=dtheta,
                half_size=half_size,
                res_thresh_method=res_thresh_method,
                res_thresh=res_thresh,
                res_axis=res_axis,
                skip_cone_region=skip_cone_region,
                bbox=b,
                clip_per=clip_per,
            )
        )
    del im
    if len(out) == 1:
        fsc_mu, res_mu, fsc_cell, res_cell = out[0]
        v = {
            "fsc_mu": _fsc_struct(fsc_mu),
            "res_mu": _res_struct(res_mu),
            "fsc_cell": [_fsc_struct(f) for f in fsc_cell],
            "res_cell": [_res_struct(r) for r in res_cell],
        }
        result = out[0]
    else:
        v = {
            "fsc_mu": _column_cell([_fsc_struct(o[0]) for o in out]),
            "res_mu": _column_cell([_res_struct(o[1]) for o in out]),
            "fsc_cell": _column_cell([[_fsc_struct(f) for f in o[2]] for o in out]),
            "res_cell": _column_cell([[_res_struct(r) for r in o[3]] for o in out]),
        }
        result = tuple([o[k] for o in out] for k in range(4))
    v["centers"] = centers
    v["bbox"] = mboxes if boxes else np.zeros((0, 0))
    tmp = f"{fnout[:-4]}_{get_uuid()}.mat"
    write_mat(tmp, v)
    os.replace(tmp, fnout)
    return result


# ---------------------------------------------------------------------------------------------
# datasets: XR_FSC_analysis_wrapper
# ---------------------------------------------------------------------------------------------


def _grid_in(v: Any) -> Any:
    """MATLAB's region centres (rows of 1-based [y x z]) as (z, y, x) 0-based rows"""
    if v is None:
        return None
    a = np.asarray(v, dtype=np.float64)
    if a.size == 0:
        return None
    a = a.reshape(-1, 3)
    return tuple((int(r[2]) - 1, int(r[0]) - 1, int(r[1]) - 1) for r in a)


def _grid_out(v: Any) -> Any:
    if v is None:
        return None
    rows = [(r[1] + 1, r[2] + 1, r[0] + 1) for r in v]
    return np.asarray(rows, dtype=np.float64).reshape(-1, 3)


def _interval_in(v: Any) -> Any:
    from petakit5d.core.axes import from_yxz

    a = np.asarray(v, dtype=np.float64).ravel()
    return tuple(float(x) for x in from_yxz(tuple(a))) if a.size == 3 else tuple(a)


def _interval_out(v: Any) -> Any:
    from petakit5d.core.axes import to_yxz

    return [float(x) for x in to_yxz(tuple(v))]


def _mkdir_warn(path: str) -> None:
    """MATLAB's ``mkdir(path)`` without outputs: a warning when the folder exists"""
    import warnings

    if os.path.isdir(path):
        warnings.warn(
            FSCWarning("Directory already exists.", "MATLAB:MKDIR:DirectoryExists"), stacklevel=3
        )
        return
    os.makedirs(path)


def _ycs(c: float, d: float, h: float, n: float, cf: float) -> np.ndarray:
    """one axis of MATLAB's region centres: ``[c : -d : h + 1](end : -1 : 2)`` then
    ``c : d : n - h + 1``; ``cf`` (the file's own centre) for d = -1"""
    if d == -1:
        return np.asarray([cf], dtype=np.float64)
    down = _colon(c, -d, h + 1)
    up = _colon(c, d, n - h + 1)
    return np.concatenate([down[::-1][: max(0, down.size - 1)], up])


def _clamp_boxes(b: np.ndarray, size: np.ndarray) -> np.ndarray:
    if np.any((b[:, 0:3] < 1) & (b[:, 3:6] > size)):
        b = b.copy()
        b[:, 0:3] = b[:, 0:3] * (b[:, 0:3] >= 1) + (b[:, 0:3] < 1)
        b[:, 3:6] = b[:, 3:6] * (b[:, 3:6] <= size) + (b[:, 3:6] > size) * size
    return b


def _region_boxes(
    im_sizes: np.ndarray, interval: list, grid: np.ndarray | None, half: list
) -> list[np.ndarray]:
    """MATLAB's boxes ([ymin xmin zmin ymax xmax zmax], 1-based) of multi_regions per file;
    sizes, the interval, the grid and the half size in MATLAB's [y x z]"""
    from petakit5d.core.dtypes import mround

    nf = im_sizes.shape[0]
    hv = np.asarray(half, dtype=np.float64)
    out: list[np.ndarray] = [np.zeros((0, 6))] * nf
    if interval is not None and len(interval) and grid is None:
        cs = mround((im_sizes + 1) / 2)
        d = [float(v) for v in interval]
        for f in range(nf):
            # (each file's centre and the half size of each axis)
            axes = []
            for k in range(3):
                c = float(cs[f, k])
                axes.append(_ycs(c, d[k], float(hv[k]), float(im_sizes[f, k]), c))
            ys, xs, zs = (np.asarray(a) for a in axes)
            Y, X, Z = np.meshgrid(ys, xs, zs, indexing="ij")
            ctr = np.stack([Y.ravel(order="F"), X.ravel(order="F"), Z.ravel(order="F")], axis=1)
            b = np.concatenate([ctr - hv, ctr + hv - 1], axis=1)
            out[f] = _clamp_boxes(b, im_sizes[f])
    if grid is not None:
        for f in range(nf):
            b = np.concatenate([grid - hv, grid + hv - 1], axis=1)
            out[f] = _clamp_boxes(b, im_sizes[f])
    return out


@register(
    "XR_FSC_analysis_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in x and y (um)."),
    Param("dz", "@isnumeric", "Scan interval (um)."),
    Param("dr", "@isnumeric", "Radial interval (voxels of the spectrum)."),
    Param("dtheta", "@isnumeric", "Angular interval (rad).", matlab_default="pi / 12"),
    Param("resThreshMethod", "@ischar", "'fixed', 'half-bit' or 'one-bit'."),
    Param("resThresh", "@isnumeric", "Threshold of 'fixed'."),
    Param(
        "halfSize",
        "@isnumeric",
        "Half size (z, y, x) of the regions (all equal).",
        convert="size_yxz",
    ),
    Param(
        "inputBbox",
        "@isnumeric",
        "Box (z0, y0, x0, z1, y1, x1) of the region (None: centre).",
        convert="bbox",
    ),
    Param("resAxis", "@ischar", "'xz' or 'yz' (x and y swapped).", convert="same"),
    Param("skipConeRegion", "@islogical", "Leave out the cones along z and x."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files."),
    Param("channels", "@isnumeric", "Wavelengths of the channels (figure names).", convert="same"),
    Param("multiRegions", "@islogical", "Several regions per image."),
    Param(
        "regionInterval",
        "@isnumeric",
        "Region interval (z, y, x) of multi_regions (-1: the centre only).",
        convert=(_interval_in, _interval_out),
    ),
    Param(
        "regionGrid",
        "@isnumeric",
        "Region centres (z, y, x), 0-based (None: from region_interval).",
        convert=(_grid_in, _grid_out),
    ),
    Param("clipPer", "@isnumeric", "Clip above this percentile (None: no)."),
    Param("suffix", "@ischar", "Suffix of the figure names."),
    Param("iterInterval", "@isnumeric", "Iteration interval of the figures."),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("cpusPerTask", "@isscalar", "CPUs of a task's job."),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only)."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="B05",
)
def fsc_analysis_dataset(
    data_paths: str | Sequence[str],
    *,
    result_dir_name: str = "FSCs",
    xy_pixel_size: float = 0.108,
    dz: float = 0.1,
    dr: float = 10,
    dtheta: float = math.pi / 12,
    res_thresh_method: str = "fixed",
    res_thresh: float = 0.2,
    half_size: Sequence[float] = (251, 251, 251),
    input_bbox: Sequence[int] | None = None,
    res_axis: str = "xz",
    skip_cone_region: bool = True,
    channel_patterns: Sequence[str] = ("tif",),
    channels: Sequence[float] = (488, 560),
    multi_regions: bool = False,
    region_interval: Sequence[float] = (-1, 50, 50),
    region_grid: Sequence[Sequence[int]] | None = None,
    clip_per: float | None = None,
    suffix: str = "decon",
    iter_interval: float = 5,
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 4,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """FSC analysis of the tiff files of datasets (``XR_FSC_analysis_wrapper``):
    ``fsc_analysis_file`` per file into ``<data path>/<result_dir_name>/<name>.mat`` through the
    job framework, then per channel ``fsc_res_info_c<c>.mat`` (``iter``, ``thetas``,
    ``res_mu_mat``, ``res_mat``, ``fsns_c``) and, with matplotlib, the figures in
    ``<data path>/figures/``."""
    from petakit5d.core.arrays import size_to_matlab
    from petakit5d.core.axes import to_yxz
    from petakit5d.core.bbox import box_to_matlab, boxes_from_matlab
    from petakit5d.core.paths import matlab_dir
    from petakit5d.core.text import through_mat2str, through_text
    from petakit5d.io.imagesize import get_image_size
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import matches_channel_pattern

    paths = [data_paths] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    if os.name == "nt":
        paths = [p.replace("\\", "/") for p in paths]
    exps = [p + "/" for p in paths]
    print(exps, flush=True)
    half_m = [float(v) for v in to_yxz(tuple(_half3(half_size, "half_size")))]

    inputs: list[str] = []
    outputs: list[str] = []
    sizes: list[list[float]] = []
    for d in exps:
        if not os.path.isdir(d):
            continue
        names = [e.name for e in matlab_dir(d + "*tif")]
        full = [d + n for n in names]
        if multi_regions:
            sizes.extend([float(v) for v in size_to_matlab(get_image_size(f))] for f in full)
        fsc_path = d + result_dir_name + "/"
        _mkdir_warn(fsc_path)
        inputs.extend(full)
        outputs.extend(fsc_path + n[:-4] + ".mat" for n in names)

    flag = [False] * len(inputs)
    for p in channel_patterns:
        flag = [a or b for a, b in zip(flag, matches_channel_pattern(inputs, p), strict=True)]
    keep = [i for i, f in enumerate(flag) if f]
    inputs = [inputs[i] for i in keep]
    outputs = [outputs[i] for i in keep]
    nf = len(inputs)

    in_box = None if input_bbox is None or not len(input_bbox) else box_to_matlab(input_bbox)
    bboxes: list[Any] = [in_box] * nf
    if multi_regions:
        im_sizes = np.asarray([sizes[i] for i in keep], dtype=np.float64).reshape(-1, 3)
        grid = None if region_grid is None else _grid_out(region_grid)
        interval = None if region_interval is None else _interval_out(region_interval)
        bboxes = list(_region_boxes(im_sizes, interval, grid, half_m))

    # the task's values as MATLAB's func strings hand them over (%0.20f, mat2str)
    xy_t = float(through_text(xy_pixel_size, "%0.20f"))
    dz_t = float(through_text(dz, "%0.20f"))
    dr_t = float(through_text(dr, "%0.20f"))
    dth_t = float(through_text(dtheta, "%0.20f"))
    thr_t = float(through_text(res_thresh, "%0.20f"))
    half_t = through_mat2str(np.asarray(half_m, dtype=np.float64))
    half_py = (float(half_t[2]), float(half_t[0]), float(half_t[1]))
    clip_t = (
        None
        if clip_per is None or not np.size(clip_per)
        else float(through_text(float(np.asarray(clip_per).ravel()[0]), "%0.20f"))
    )
    tasks = []
    for f in range(nf):
        b = bboxes[f]
        if b is None or np.size(b) == 0:
            box_py = None
        else:
            bm = np.asarray(through_mat2str(np.asarray(b, dtype=np.float64)), dtype=np.float64)
            box_py = [list(map(int, r)) for r in boxes_from_matlab(bm.reshape(-1, 6))]
        tasks.append(
            Task.of(
                fsc_analysis_file,
                inputs[f],
                outputs[f],
                xy_pixel_size=xy_t,
                dz=dz_t,
                dr=dr_t,
                dtheta=dth_t,
                half_size=list(half_py),
                res_thresh_method=res_thresh_method,
                res_thresh=thr_t,
                res_axis=res_axis,
                skip_cone_region=bool(skip_cone_region),
                bbox=box_py,
                clip_per=clip_t,
            )
        )
    mem_allocate = float(np.prod(half_m)) * 4 * 200 / 1024**3
    framework.run_tasks(
        inputs,
        outputs,
        tasks,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        parse_cluster=parse_cluster,
        master_compute=master_compute,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )
    _summaries(
        exps,
        result_dir_name,
        channel_patterns,
        channels,
        dtheta,
        iter_interval,
        suffix,
        multi_regions,
    )


def _summaries(
    exps, result_dir_name, channel_patterns, channels, dtheta, iter_interval, suffix, multi_regions
) -> None:
    """the per-channel summary files and figures of XR_FSC_analysis_wrapper"""
    from petakit5d.core.paths import matlab_dir
    from petakit5d.io import read_mat, write_mat
    from petakit5d.matlabcompat import reduce as R
    from petakit5d.pipeline.filenames import matches_channel_pattern
    from petakit5d.tools import figures

    it = float(iter_interval)
    ntheta = int(round(2 * math.pi / dtheta)) + 1
    draw = figures.available()
    chans = np.asarray(channels, dtype=np.float64).ravel()
    for d in exps:
        names = [e.name for e in matlab_dir(d + "*tif")]
        full = [d + n for n in names]
        fig_path = d + "figures/"
        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)
        fsc_path = d + result_dir_name + "/"
        for c, pat in enumerate(channel_patterns):
            hit = matches_channel_pattern(full, pat)
            fs = [n for n, h in zip(names, hit, strict=True) if h]
            if not fs:
                continue
            n_reg = 1
            if multi_regions:
                first = read_mat(fsc_path + fs[0][:-4] + ".mat")
                n_reg = int(np.asarray(first["bbox"]).reshape(-1, 6).shape[0])
            res_mat = np.zeros((len(fs), ntheta, n_reg))
            for i, n in enumerate(fs):
                a = read_mat(fsc_path + n[:-4] + ".mat")
                if not multi_regions:
                    rv = np.asarray(a["res_mu"]["res"], dtype=np.float64).ravel()
                    if rv.size != ntheta:
                        raise FSCError(
                            "Unable to perform assignment because the size of the left side is "
                            f"1-by-{ntheta} and the size of the right side is {rv.size}-by-1.",
                            "MATLAB:subsassigndimmismatch",
                        )
                    res_mat[i, :, 0] = rv
                    thetas = a["res_mu"]["thetas"]
                else:
                    rm = a["res_mu"]
                    # (one region's result holds its struct itself, not a cell of them)
                    cells = (
                        [rm]
                        if isinstance(rm, dict)
                        else list(np.asarray(rm, dtype=object).ravel(order="F"))
                    )
                    for j in range(n_reg):
                        res_mat[i, :, j] = np.asarray(cells[j]["res"], dtype=np.float64).ravel()
                    thetas = cells[0]["thetas"]
            mm = res_mat if multi_regions else res_mat[:, :, 0]
            # (the mean of the measured angles: the sectors inside the cones are NaN)
            res_mu_mat = np.asarray(
                R.mean(
                    np.ascontiguousarray(mm[:, :-1, ...]),
                    axis=1,
                    keepdims=True,
                    matlab_layout=True,
                    omitnan=bool(np.isnan(mm).any()),
                )
            )
            write_mat(
                f"{fsc_path}/fsc_res_info_c{c + 1}.mat",
                {
                    "iter": it,
                    "thetas": np.asarray(thetas, dtype=np.float64),
                    "res_mu_mat": res_mu_mat,
                    "res_mat": mm,
                    "fsns_c": _column_cell(fs),
                },
            )
            if draw:
                ch = float(chans[c]) if c < chans.size else math.nan
                th = np.asarray(thetas, dtype=np.float64).ravel()
                _draw_figures(fig_path, fs, mm, th, it, ch, suffix, multi_regions)


def _fill_linear(a: np.ndarray) -> np.ndarray:
    """``fillmissing(a, 'linear', 2)`` of rows (NaN filled by the line through the neighbouring
    values, extrapolated at the ends); for the figures only"""
    out = np.array(a, dtype=np.float64)
    x = np.arange(out.shape[1], dtype=np.float64)
    for row in out.reshape(-1, out.shape[1]) if out.ndim == 2 else []:
        ok = ~np.isnan(row)
        if ok.sum() >= 2 and not ok.all():
            xs, ys = x[ok], row[ok]
            v = np.interp(x, xs, ys)
            lo, hi = x < xs[0], x > xs[-1]
            v[lo] = ys[0] + (x[lo] - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
            v[hi] = ys[-1] + (x[hi] - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
            row[:] = v
    return out


def _save_png(fig, path: str) -> None:
    from petakit5d.core.ids import get_uuid

    root, ext = os.path.splitext(path)
    tmp = f"{root}_{get_uuid()}{ext}"
    fig.savefig(tmp, format="png", facecolor="w")
    os.replace(tmp, path)


def _draw_figures(fig_path, names, res_mat, thetas, it, channel, suffix, multi) -> None:
    """XR_FSC_analysis_wrapper's figures (matplotlib): the mean resolution against the
    iteration (``iter : iter : n iter``) with its minimum, again with y in [0.2, 0.5], and the
    resolution against the angle of each file (polar, r in [0, 2])"""
    from matplotlib.figure import Figure

    from petakit5d.core.text import sprintf

    r = np.array(res_mat, dtype=np.float64)
    r[np.isinf(r)] = np.nan
    if multi:
        r = np.stack([_fill_linear(r[:, :, j]) for j in range(r.shape[2])], axis=2)
        per = r[:, :-1, :].mean(axis=1)
        mu = per.mean(axis=1)
    else:
        r = _fill_linear(r)
        per = None
        mu = r[:, :-1].mean(axis=1)
    n = mu.size
    xs = it * np.arange(1, n + 1)
    fig = Figure(figsize=(12, 8), dpi=100, facecolor="w")
    ax = fig.add_subplot(1, 1, 1)
    if per is not None:
        ax.plot(xs, per, "--")
    ax.plot(xs, mu)
    mind = int(np.nanargmin(mu)) if np.isfinite(mu).any() else 0
    ax.plot([xs[mind]], [mu[mind]], "o")
    top = np.nanmax(mu) if np.isfinite(mu).any() else 0.0
    ty = min(mu[mind] * 1.3, mu[mind] * 0.9 + top * 0.1)
    if per is None:
        ax.plot(it * np.arange(0, n + 1), np.full(n + 1, mu[mind] * 1.01), "--")
    ax.text(xs[mind], ty, sprintf("Iter = %d", xs[mind]))
    if per is not None:
        sd = per.std(axis=1, ddof=1) if per.shape[1] > 1 else np.zeros(n)
        ms = mind + int(np.nanargmin(np.abs(mu[mind:] - (mu[mind] + sd[mind]))))
        ty2 = min(mu[ms] * 1.3, mu[ms] * 0.9 + top * 0.1)
        ax.text(xs[ms], ty2, sprintf("Iter = %d (mu+std)", xs[ms]))
    ax.grid(True)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Resolution (um)")
    base = sprintf(
        "%s/fsc_average_resolution_vs_decon_iterations_ch_%d_%s", fig_path, channel, suffix
    )
    _save_png(fig, base + ".png")
    ax.set_ylim(0.2, 0.5)
    ax.set_yticks(_colon(0.2, 0.025, 0.5))
    _save_png(fig, base + "_res_0p2_0p5.png")
    for i, name in enumerate(names):
        fig = Figure(figsize=(12, 8), dpi=100, facecolor="w")
        ax = fig.add_subplot(1, 1, 1, projection="polar")
        ax.plot(thetas, r[i].mean(axis=1) if multi else r[i])
        ax.grid(True)
        ax.set_rlim(0, 2)
        ax.set_title(sprintf("Iter %d", (i + 1) * it))
        _save_png(fig, sprintf("%s/fsc_res_vs_thetas_%s_%s.png", fig_path, name[:-4], suffix))


# ---------------------------------------------------------------------------------------------
# two images: XR_two_image_FSC_analysis
# ---------------------------------------------------------------------------------------------


def _isotropic_cubic(im: np.ndarray, xy: float, dz: float, size: tuple) -> np.ndarray:
    """resample_crop_pad_image_isotropic: the image (double) resized to isotropic voxels
    (MATLAB's ``imresize3``, cubic), cropped about its centre (``floor((size - n) / 2)``) or
    padded (``floor`` before) to ``size``"""
    from petakit5d.core.dtypes import mround
    from petakit5d.matlabcompat.resize import imresize3

    psz = min(xy, dz)
    if dz != xy:
        im = imresize3(
            im,
            tuple(int(mround((s * f) / psz)) for s, f in zip(im.shape, (dz, xy, xy), strict=True)),
        )
    sz = im.shape
    if any(s > n for s, n in zip(sz, size, strict=True)):
        s0 = [max(1, math.floor((s - n) / 2)) for s, n in zip(sz, size, strict=True)]
        t = [int(min(s, a + n - 1)) for s, a, n in zip(sz, s0, size, strict=True)]
        im = im[s0[0] - 1 : t[0], s0[1] - 1 : t[1], s0[2] - 1 : t[2]]
        sz = im.shape
    if any(s < n for s, n in zip(sz, size, strict=True)):
        pre = [math.floor((n - s) / 2) for s, n in zip(sz, size, strict=True)]
        post = [int(n - s - a) for s, n, a in zip(sz, size, pre, strict=True)]
        im = np.pad(im, list(zip(pre, post, strict=True)))
    return np.ascontiguousarray(im)


@register(
    "XR_two_image_FSC_analysis",
    Param(
        "im1",
        "@(x) ischar(x) || isnumeric(x)",
        "Tiff file, or the image (z, y, x).",
        convert=(_path_or_volume_in, _path_or_volume_out),
        py="im1",
    ),
    Param(
        "im2",
        "@(x) ischar(x) || isnumeric(x)",
        "Tiff file, or the image (z, y, x).",
        convert=(_path_or_volume_in, _path_or_volume_out),
        py="im2",
    ),
    Param("xyPixelSize", "@isnumeric", "Pixel size in x and y (um)."),
    Param("dz", "@isnumeric", "Voxel size in z (um)."),
    Param("dr", "@isnumeric", "Radial interval (voxels of the spectrum)."),
    Param("dtheta", "@isnumeric", "Angular interval (rad).", matlab_default="pi / 12"),
    Param("resThreshMethod", "@ischar", "'fixed', 'one-bit' or 'half-bit'."),
    Param("resThresh", "@isnumeric", "Threshold of 'fixed'."),
    Param("N", "@isnumeric", "Half size (z, y, x).", convert="size_yxz"),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    task="B05",
)
def two_image_fsc_analysis(
    im1: str | np.ndarray,
    im2: str | np.ndarray,
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.108,
    dr: float = 1,
    dtheta: float = math.pi / 12,
    res_thresh_method: str = "one-bit",
    res_thresh: float = 0.2,
    n: Sequence[float] = (501, 501, 501),
    debug: bool = False,
    num_threads: int | None = None,
) -> tuple[dict, dict]:
    """The FSC of two images (``XR_two_image_FSC_analysis``; tiff files or ``(z, y, x)``
    arrays): each resized to isotropic voxels in double (MATLAB's ``imresize3``, cubic), cropped
    about its centre or padded to twice ``n``, the FSC (``fourier_shell_correlation``, its
    defaults) of each of the eight sub-lattices of every other voxel of the first image with the
    same of the second, and the resolution of their mean: ``(fsc, res)`` (``fsc`` the last
    pair's, with the mean curves)."""
    from petakit5d.matlabcompat import reduce as R

    half = _half3(n, "n")
    vols = []
    for im in (im1, im2):
        a = _read_image(os.fspath(im)) if isinstance(im, (str, os.PathLike)) else np.asarray(im)
        if a.ndim == 2:
            a = a[None]
        a = np.asarray(a, dtype=np.float64)
        vols.append(
            _isotropic_cubic(a, float(xy_pixel_size), float(dz), tuple(2 * h for h in half))
        )
    kw = {
        "xy_pixel_size": xy_pixel_size,
        "dz": dz,
        "dr": dr,
        "dtheta": dtheta,
        "half_size": n,
        "num_threads": num_threads,
    }
    fsc_cell = []
    fsc_i = None
    for i in range(8):
        # MATLAB's [y, x, z] = ind2sub([2, 2, 2], i)
        par = (i % 2 + 1, (i // 2) % 2 + 1, i // 4 + 1)
        fsc_i, _ = fourier_shell_correlation(
            _sub_lattice(vols[0], par), _sub_lattice(vols[1], par), **kw
        )
        fsc_cell.append(fsc_i)
    shapes = {np.shape(f["fsc"]) for f in fsc_cell}
    if len(shapes) != 1:
        raise FSCError(
            "Dimensions of arrays being concatenated are not consistent.",
            "MATLAB:catenate:dimensionMismatch",
        )
    stack = np.stack([np.asarray(f["fsc"], dtype=np.float64) for f in fsc_cell], axis=2)
    fsc_mu = dict(fsc_i)
    fsc_mu["fsc"] = np.ascontiguousarray(
        np.asarray(R.mean(stack, axis=2, matlab_layout=True)).reshape(stack.shape[:2])
    )
    res_mu = fsc_resolution(
        fsc_mu,
        xy_pixel_size=xy_pixel_size,
        dz=dz,
        res_thresh_method=res_thresh_method,
        res_thresh=res_thresh,
        dtheta=dtheta,
        debug=debug,
    )
    return fsc_mu, res_mu
