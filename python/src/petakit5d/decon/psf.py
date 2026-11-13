"""PSFs for deconvolution: cleaning, centering and resampling (``psf_gen_new``), the crops
RLdecon and the decon wrapper make, the cached reader of PSF files (``decon_readtiff_cached``)
and the rotation of PSF files (``XR_rotate_PSF``).

Arrays are the port's: a PSF is a ``(z, y, x)`` volume or a ``(y, x)`` image (one z plane, as
MATLAB's 2D PSF), C-contiguous. Results equal MATLAB's bit for bit (psf_gen_new.m,
decon_readtiff_cached.m, XR_rotate_PSF.m and the PSF lines of RLdecon.m and
XR_decon_data_wrapper.m at the oracle commit):

- **psf_gen** (``generate_psf`` with the FFT backend and thread count): the PSF in single; the
  background step of the method, MATLAB's for volumes with a positive voxel in their first or
  last five z slices and for images:

  - ``'masked'`` (the default): the 3x3x3 median (``medfilt3``, omw_backprojector_mex's for
    volumes; ``medfilt2`` with zero padding for images); the threshold per x and z from the
    first and last ten y rows of the median, ``3 max(sqrt(|m - 100|)) + mean(m)`` over them
    (``3 max(sqrt(m)) + mean(m)`` for images, complex as MATLAB's ``sqrt`` where a value is
    negative); the voxels above it with their 26-connected (8 for images) components of at
    least 300 voxels (50); the component of the median's first maximum (all voxels outside
    such components where the maximum is in none, as ``L == L(peakInd)`` with label 0); its
    closing with ``strel('sphere', 3)`` (``strel('disk', 3)``); the PSF times the mask minus
    the mean of its first and last ten y rows per x, negative values set to 0.
  - ``'median'``: the PSF minus ``medFactor`` times the median of the positive values of the
    first and last five z slices (of the whole image), negative values set to 0 (``max``:
    NaN to 0, -0 kept), the component of the median's peak (``medfilt3`` of the double data,
    symmetric padding; ``medfilt2``) above 0 as above (30 voxels for images), closed, times the
    PSF.
  - another name: nothing (MATLAB's switch has no case for it).

  Then the PSF is padded with zero slices where its nonzero slices reach further than half its
  length from the peak (the first maximum in MATLAB's order), shifted circularly so that the
  peak is at ``size // 2`` (0-based), and, for ``dz_data / dz_psf > 1``, resampled along z by
  truncating its spectrum: ``fftn`` in single, ``new_nz = uint16(round(nz / ratio))``, the first
  ``h = new_nz / 2`` (uint16 division: rounded half away from zero) and the last ``h`` planes of
  the spectrum in double, the second over the first where they overlap (odd ``new_nz``), and
  ``real(ifftn)`` in double: a resampled PSF is double, as MATLAB's. Negative values are set to
  0 at the end. The FFTs run on ``backend`` (``petakit5d.fft``): the exact one gives MATLAB's
  bits, the permissive one (the default) the last bits of its library (decision record 0006).
- **crop_psf_support** (RLdecon and XR_decon_data_wrapper, methods other than OMW): the box of
  the nonzero planes along each axis, one more plane on each side, as wide on both sides of the
  center as the nearer side allows. **crop_psf_to_data** (RLdecon): a PSF larger than the data
  along an axis is cut to the data's length about its center, with MATLAB's warning.
- **read_tiff_cached** (``decon_readtiff_cached``): the last four files read are kept with their
  modification time and size; a file is read again once either changes. Cached arrays are
  read-only (MATLAB's arrays are values); copy one to change it.
- **rotate_psf** (``XR_rotate_PSF``): the PSF file rotated (``deskew_rotate.rotate_volume``),
  its zero voxels set to the median of the positive values below their ``HighPrctile``
  percentile, written as ``<folder>/Rotated/<name>.tif`` (uint16 or single) with the
  parameters.

Nothing here depends on the thread count.
"""

from __future__ import annotations

import os
import stat
import threading
import warnings
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "PSFError",
    "clear_read_cache",
    "crop_psf_support",
    "crop_psf_to_data",
    "generate_psf",
    "psf_gen",
    "read_tiff_cached",
    "rotate_psf",
]

_BAD_SUBSCRIPT = "MATLAB:badsubscript"
_INVALID_INDEX = "MATLAB:badsubscript"  # also for indices that are not positive
_AREA = {("masked", 3): 300, ("median", 3): 300, ("masked", 2): 50, ("median", 2): 30}


class PSFError(IndexError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# small MATLAB rules, and passes over planes on threads
# ---------------------------------------------------------------------------------------------


def _par(fn: Callable[[int, int], None], n: int, threads: int) -> None:
    """``fn(start, stop)`` on chunks of ``range(n)`` on up to ``threads`` threads (NumPy's
    element-wise loops and the kernels release the interpreter lock)"""
    k = min(threads, n)
    if k <= 1:
        fn(0, n)
        return
    edges = [n * i // k for i in range(k + 1)]
    with ThreadPoolExecutor(k) as ex:
        for f in [ex.submit(fn, edges[i], edges[i + 1]) for i in range(k)]:
            f.result()


def _single(x: Any, threads: int) -> np.ndarray:
    """MATLAB's single(x), C-contiguous (no copy of single data)"""
    a = np.asarray(x)
    if a.dtype == np.float32:
        return np.ascontiguousarray(a)
    if a.dtype.kind in "biu" and a.dtype.itemsize < 8 or a.dtype == np.float64:
        from petakit5d.kernels import ops

        return ops.convert(np.ascontiguousarray(a), "single", num_threads=threads)
    from petakit5d.core import dtypes

    return np.ascontiguousarray(dtypes.single(a))


def _edges(n: int, k: int, position: int) -> np.ndarray:
    """``[1 : k, end - k + 1 : end]`` along an axis of length n as 0-based indices, with MATLAB's
    errors of the indexing (``position``: the index's position in MATLAB's subscript)"""
    idx = np.r_[np.arange(1, k + 1), np.arange(n - k + 1, n + 1)]
    if idx.min() < 1:
        raise PSFError(
            _INVALID_INDEX,
            f"Index in position {position} is invalid. Array indices must be positive integers "
            "or logical values.",
        )
    if idx.max() > n:
        raise PSFError(
            _BAD_SUBSCRIPT,
            f"Index in position {position} exceeds array bounds. Index must not exceed {n}.",
        )
    return idx - 1


class _Planes:
    """per z plane of a (z, y, x) volume: the largest value (NaN ignored) and whether a value is
    nonzero (NaN counts), filled by ``take`` on the planes of a pass"""

    def __init__(self, nz: int, dtype):
        self.max = np.empty(nz, dtype=dtype)
        self.nonzero = np.empty(nz, dtype=bool)

    def take(self, vol: np.ndarray, s: int, e: int, nonzero: bool = True) -> None:
        part = vol[s:e].reshape(e - s, -1)
        self.max[s:e] = np.fmax.reduce(part, axis=1)
        if nonzero:
            self.nonzero[s:e] = [np.count_nonzero(p) > 0 for p in part]

    def first_max(self, vol: np.ndarray) -> tuple[int, int, int]:
        """the (z, y, x) index of ``max(vol(:))``: the first maximum in MATLAB's order (z, then
        x, then y), NaN ignored unless every value is NaN"""
        v = np.fmax.reduce(self.max)
        if np.isnan(v):
            from petakit5d.matlabcompat import reduce

            return reduce.max(vol)[1]
        z = int(np.argmax(self.max == v))
        k = int(np.argmax(np.ascontiguousarray((vol[z] == v).T)))
        x, y = divmod(k, vol.shape[1])
        return z, y, x


def _place(src: np.ndarray, nz2: int, pre: int, shift: Sequence[int], threads: int):
    """``circshift(padarray(src, [0 0 pre], 'pre') padded to nz2 planes, shift)``: src's planes at
    ``pre`` in an array of nz2 planes (the others +0), shifted circularly by ``shift``
    (z, y, x), plane by plane on threads"""
    nz, ny, nx = src.shape
    out = np.empty((nz2, ny, nx), dtype=src.dtype)
    d0 = (pre + shift[0]) % nz2  # the plane src[0] goes to
    covered = np.zeros(nz2, dtype=bool)
    covered[(d0 + np.arange(nz)) % nz2] = True
    out[~covered] = 0
    ry, rx = shift[1] % ny, shift[2] % nx

    def copy(s: int, e: int) -> None:
        for i in range(s, e):
            p, q = src[i], out[(d0 + i) % nz2]
            q[ry:, rx:] = p[: ny - ry, : nx - rx]
            q[ry:, :rx] = p[: ny - ry, nx - rx :]
            q[:ry, rx:] = p[ny - ry :, : nx - rx]
            q[:ry, :rx] = p[ny - ry :, nx - rx :]

    _par(copy, nz, threads)
    return out


# ---------------------------------------------------------------------------------------------
# the background step
# ---------------------------------------------------------------------------------------------


def _mask_and_peak(med: np.ndarray, a, threads: int) -> tuple[np.ndarray, tuple[int, ...]]:
    """``med - a > 0`` (``med > a``, the same for every value with gradual underflow) and the
    first maximum of med, in one pass on threads"""
    mask = np.empty(med.shape, dtype=bool)
    planes = _Planes(med.shape[0], med.dtype)
    full = np.ndim(a) == 3

    def run(s: int, e: int) -> None:
        np.greater(med[s:e], a[s:e] if full else a, out=mask[s:e])
        planes.take(med, s, e, nonzero=False)

    _par(run, med.shape[0], threads)
    return mask, planes.first_max(med)


def _peak_component(mask: np.ndarray, peak: tuple[int, ...], area: int, threads: int):
    """``BW = bwareaopen(mask, area, conn); L = bwlabeln(BW); BW = L == L(peakInd)`` with the
    maximal connectivity (26, 8 for an image): the component of the median's first maximum if
    it has at least ``area`` voxels, else (L(peakInd) is 0) every voxel outside the components
    that large; with the planes, rows and columns it reaches"""
    from petakit5d.matlabcompat import label

    lab, n = label._label(mask, 26 if mask.shape[0] > 1 else 8, "any", threads)
    lab = lab.reshape(mask.shape)
    lp = int(lab[peak])
    bw = np.empty(mask.shape, dtype=bool)
    nz = mask.shape[0]
    count = np.zeros(nz, dtype=np.int64)
    if lp:

        def run(s: int, e: int) -> None:
            np.equal(lab[s:e], lp, out=bw[s:e])
            count[s:e] = [np.count_nonzero(p) for p in bw[s:e]]

        _par(run, nz, threads)
        if int(count.sum()) >= area:
            return bw
    sizes = np.bincount(lab.reshape(-1), minlength=n + 1)
    keep = sizes >= area
    keep[0] = False

    def rest(s: int, e: int) -> None:
        np.logical_not(keep[lab[s:e]], out=bw[s:e])

    _par(rest, nz, threads)
    return bw


def _close(bw: np.ndarray, threads: int) -> np.ndarray:
    """``imclose(bw, strel('sphere', 3))`` (``strel('disk', 3)`` for an image): the closing of
    the mask extended with false, computed on the mask's bounding box widened by the
    neighborhood (the closing is inside the dilation, so it is false beyond that box)"""
    from petakit5d.matlabcompat import morphology as M

    vol = bw.shape[0] > 1
    se = M.strel("sphere", 3) if vol else M.strel("disk", 3)
    img = bw if vol else bw[0]
    # the neighborhood's reach from its center floor((size + 1) / 2) (1-based) on either side
    reach = [max((s + 1) // 2 - 1, s - (s + 1) // 2) for s in se.neighborhood.shape]
    nd = img.ndim
    hits = [np.zeros(n, dtype=bool) for n in img.shape]
    hits_parts: list = []

    def run(s: int, e: int) -> None:
        part = img[s:e]
        hits_parts.append(
            (s, [np.any(part, axis=tuple(a for a in range(nd) if a != ax)) for ax in range(nd)])
        )

    _par(run, img.shape[0], threads)
    for s, h in hits_parts:
        hits[0][s : s + h[0].size] |= h[0]
        for ax in range(1, nd):
            hits[ax] |= h[ax]
    box = []
    for ax in range(nd):
        hit = np.flatnonzero(hits[ax])
        if hit.size == 0:
            return bw.copy()
        box.append(
            (max(int(hit[0]) - reach[ax], 0), min(int(hit[-1]) + 1 + reach[ax], img.shape[ax]))
        )
    sl = tuple(slice(a, b) for a, b in box)
    out = np.zeros(img.shape, dtype=bool)
    out[sl] = M.imclose(np.ascontiguousarray(img[sl]), se, threads=threads)
    return out.reshape(bw.shape)


def _masked_threshold(med_rows: np.ndarray, image: bool, threads: int) -> np.ndarray:
    """``max(sqrt(abs(m - 100)), [], 1) * 3 + mean(m, 1)`` of the rows (``sqrt(m)`` for an image)
    as (z, 1, x) single; for an image with a negative value the real part of MATLAB's complex
    result (the largest magnitude, on a tie the larger phase angle: the negative value, whose
    square root is imaginary)"""
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat import reduce

    if image:
        v = med_rows
        neg = v < 0
        if neg.any():
            mag = np.sqrt(np.abs(v))
            best = reduce.max(mag, axis=1, keepdims=True)[0]
            top = mag == best
            imag = np.any(top & neg, axis=1, keepdims=True)
            root = np.where(imag, np.float32(0), best)
        else:
            root = reduce.max(np.sqrt(v), axis=1, keepdims=True)[0]
    else:
        v = np.abs(dtypes.minus(med_rows, 100.0))
        root = reduce.max(np.sqrt(v), axis=1, keepdims=True)[0]
    mean = reduce.mean(med_rows, axis=1, keepdims=True, num_threads=threads)
    return dtypes.plus(dtypes.times(root, 3.0), mean)


def _background_masked(raw: np.ndarray, out: np.ndarray, planes: _Planes, threads: int) -> None:
    """the 'masked' step (see the module) of raw (nz, ny, nx) single into out, with the planes'
    statistics of the result"""
    from petakit5d.matlabcompat import morphology as M
    from petakit5d.matlabcompat import reduce

    image = raw.shape[0] == 1
    if image:
        med = M.medfilt2(raw[0], (3, 3), threads=threads)[None]
    else:
        med = M.medfilt3(raw, (3, 3, 3), threads=threads)
    rows = _edges(raw.shape[1], 10, 1)
    a = _masked_threshold(np.ascontiguousarray(med[:, rows, :]), image, threads)
    mask, peak = _mask_and_peak(med, a, threads)
    del med
    bw = _close(_peak_component(mask, peak, _AREA["masked", 2 if image else 3], threads), threads)
    m = reduce.mean(np.ascontiguousarray(raw[:, rows, :]), axis=(0, 1), num_threads=threads)

    def run(s: int, e: int) -> None:
        o = out[s:e]
        np.multiply(raw[s:e], bw[s:e], out=o, dtype=np.float32)
        np.subtract(o, m, out=o)
        np.copyto(o, np.float32(0), where=o < 0)
        planes.take(out, s, e)

    _par(run, raw.shape[0], threads)


def _positive_box(p: np.ndarray, threads: int) -> list[tuple[int, int]] | None:
    """the bounding box (per axis (start, stop)) of the positive values of a (z, y, x) volume,
    None without one"""
    nd = p.ndim
    hits = [np.zeros(n, dtype=bool) for n in p.shape]
    parts: list = []

    def run(s: int, e: int) -> None:
        pos = p[s:e] > 0
        parts.append(
            (s, [np.any(pos, axis=tuple(a for a in range(nd) if a != ax)) for ax in range(nd)])
        )

    _par(run, p.shape[0], threads)
    for s, h in parts:
        hits[0][s : s + h[0].size] |= h[0]
        for ax in range(1, nd):
            hits[ax] |= h[ax]
    box = []
    for h in hits:
        i = np.flatnonzero(h)
        if i.size == 0:
            return None
        box.append((int(i[0]), int(i[-1]) + 1))
    return box


def _median_mask_and_peak(p: np.ndarray, threads: int) -> tuple[np.ndarray, tuple[int, ...]]:
    """``psf_med = medfilt3(double(p), [3, 3, 3])`` (symmetric padding) of a volume without
    negative values or NaN, as the 'median' method uses it: ``psf_med > 0`` and the first
    maximum of psf_med. A median is positive exactly where 14 or more of the 27 values of its
    window are, elsewhere 0; so the windows are counted on the box of the positive values
    (widened by 1: the windows beyond it hold none) and the medians are taken there only, as
    the 14th smallest of their windows' values. Where that box is large, medfilt3 itself."""
    from itertools import product

    from petakit5d.kernels import ops
    from petakit5d.matlabcompat import morphology as M

    shape = p.shape
    box = _positive_box(p, threads)
    if box is None:  # every median is 0: the first one is the first maximum
        return np.zeros(shape, dtype=bool), (0, 0, 0)
    box = [(max(a - 1, 0), min(b + 1, n)) for (a, b), n in zip(box, shape, strict=True)]
    if 8 * int(np.prod([b - a for a, b in box])) > p.size:
        med = M.medfilt3(ops.convert(p, "double", num_threads=threads), (3, 3, 3), threads=threads)
        return _mask_and_peak(med, 0, threads)
    sl = tuple(slice(a, b) for a, b in box)
    c = (p[sl] > 0).astype(np.uint8)
    # the windows' values beyond the box: the image's edge repeated (symmetric padding of one),
    # else no positive value
    for ax, ((a, b), n) in enumerate(zip(box, shape, strict=True)):
        first = np.take(c, [0], axis=ax) if a == 0 else np.zeros_like(np.take(c, [0], axis=ax))
        last = np.take(c, [-1], axis=ax) if b == n else np.zeros_like(first)
        c = np.concatenate([first, c, last], axis=ax)
    for ax in range(3):  # the 3 x 3 x 3 window sums, one axis after the other
        n = c.shape[ax] - 2
        part = [slice(None)] * 3
        terms = []
        for o in range(3):
            part[ax] = slice(o, o + n)
            terms.append(c[tuple(part)])
        c = terms[0] + terms[1] + terms[2]
    inner = c >= 14
    mask = np.zeros(shape, dtype=bool)
    mask[sl] = inner
    zi, yi, xi = np.nonzero(inner)
    if zi.size == 0:
        return mask, (0, 0, 0)
    z, y, x = zi + box[0][0], yi + box[1][0], xi + box[2][0]
    vals = np.empty((z.size, 27), dtype=np.float64)
    for k, (dz, dy, dx) in enumerate(product((-1, 0, 1), repeat=3)):
        vals[:, k] = p[
            np.clip(z + dz, 0, shape[0] - 1),
            np.clip(y + dy, 0, shape[1] - 1),
            np.clip(x + dx, 0, shape[2] - 1),
        ]
    med = np.partition(vals, 13, axis=1)[:, 13]
    top = np.flatnonzero(med == med.max())
    lin = y[top] + shape[1] * (x[top] + shape[2] * z[top])  # MATLAB's linear index
    k = top[int(np.argmin(lin))]
    return mask, (int(z[k]), int(y[k]), int(x[k]))


def _background_median(
    raw: np.ndarray, sel: np.ndarray, med_factor: Any, out: np.ndarray, planes, threads: int
) -> None:
    """the 'median' step (see the module) of raw (nz, ny, nx) single into out; sel the voxels
    whose positive values give the median"""
    from petakit5d.core import dtypes
    from petakit5d.kernels import ops
    from petakit5d.matlabcompat import morphology as M
    from petakit5d.matlabcompat import stats

    image = raw.shape[0] == 1
    # medFactor * median(...): single, computed in double (the median of nothing is NaN)
    level = np.float32(dtypes.times(med_factor, stats.median(sel[sel > 0])))

    def sub(s: int, e: int) -> None:
        np.subtract(raw[s:e], level, out=out[s:e])

    _par(sub, raw.shape[0], threads)
    ops.clamp(out, lo=0.0, out=out, num_threads=threads)  # max(psf_raw, 0)
    if image:
        mask, peak = _mask_and_peak(M.medfilt2(out[0], (3, 3), threads=threads)[None], 0, threads)
    else:
        mask, peak = _median_mask_and_peak(out, threads)
    bw = _close(_peak_component(mask, peak, _AREA["median", 2 if image else 3], threads), threads)

    def run(s: int, e: int) -> None:
        o = out[s:e]
        np.multiply(o, bw[s:e], out=o, dtype=np.float32)
        planes.take(out, s, e)

    _par(run, raw.shape[0], threads)


# ---------------------------------------------------------------------------------------------
# resampling along z
# ---------------------------------------------------------------------------------------------


def _u16(x: float) -> int:
    """MATLAB's uint16 of a double: rounded half away from zero, saturated"""
    from petakit5d.core import dtypes

    return int(dtypes.uint16(float(x)))


def _mround(x: float) -> float:
    from petakit5d.core import dtypes

    return float(dtypes.mround(x))


def _truncate_z(psf: np.ndarray, ratio: float, backend: str | None, threads: int) -> np.ndarray:
    """psf_gen_new's Fourier truncation along z with MATLAB's uint16 index arithmetic:
    ``real(ifftn(trunc))`` in double of the single spectrum's first and last ``new_nz / 2``
    planes (the second over the first where they overlap)"""
    from petakit5d.matlabcompat import fft as F

    nz, ny, nx = psf.shape
    spec = F.fftn(psf, backend=backend, threads=threads)
    new_nz = _u16(_mround(nz / ratio))
    h = _u16(new_nz / 2)
    lo = _u16(_u16(new_nz - h) + 1)  # new_nz - new_nz / 2 + 1, in uint16
    src = _u16(_u16(nz - h) + 1)  # nz - new_nz / 2 + 1, in uint16
    n_dst, n_src = max(new_nz - lo + 1, 0), max(_u16(nz) - src + 1, 0)
    if n_dst != n_src:
        raise PSFError(
            "MATLAB:subsassigndimmismatch",
            f"Unable to perform assignment because the size of the left side is {ny}-by-{nx}-by-"
            f"{n_dst} and the size of the right side is {ny}-by-{nx}-by-{n_src}.",
        )
    # the plane of the spectrum each plane of trunc gets: 1 : h, then lo : new_nz over them
    take = np.full(new_nz, -1, dtype=np.int64)
    take[:h] = np.arange(h)
    take[lo - 1 : new_nz] = np.arange(src - 1, src - 1 + n_src)
    trunc = np.empty((new_nz, ny, nx), dtype=np.complex128)

    def fill(s: int, e: int) -> None:
        for d in range(s, e):
            if take[d] < 0:
                trunc[d] = 0
            else:
                trunc[d] = spec[take[d]]

    _par(fill, new_nz, threads)
    spec = None  # freed before the inverse transform
    out = F.ifftn(trunc, backend=backend, threads=threads)
    return np.ascontiguousarray(out.real) if np.iscomplexobj(out) else out


# ---------------------------------------------------------------------------------------------
# psf_gen_new
# ---------------------------------------------------------------------------------------------


def generate_psf(
    psf,
    dz_psf: float,
    dz_data: float,
    med_factor: float = 1.5,
    method: str = "masked",
    *,
    backend: str | None = None,
    threads: int | None = None,
) -> np.ndarray:
    """``psf_gen_new(psf, dz_psf, dz_data, medFactor, PSFGenMethod)`` (see the module): the
    cleaned PSF, centered at ``size // 2`` and resampled to ``dz_data`` along z. ``psf``: a
    (z, y, x) volume or a (y, x) image of any real class; the result has its number of
    dimensions (a resampled image with no planes left is ``(0, ny, nx)``), single, or double when
    resampled. ``backend``: the FFT backend of the resampling (``petakit5d.fft``; None: the
    default); ``threads``: None for all available CPUs."""
    from petakit5d.core.threads import resolve_threads
    from petakit5d.kernels import ops

    t = resolve_threads(threads)
    raw = _single(psf, t)
    if raw.ndim not in (2, 3):
        raise ValueError("the PSF must be a (z, y, x) volume or a (y, x) image")
    image_in = raw.ndim == 2
    vol = raw[None] if image_in else raw
    nz, ny, nx = vol.shape
    m = method.lower()
    step = m in ("masked", "median")
    if step and nz > 1:
        edge = vol[_edges(nz, 5, 3)]
        step = bool(np.any(edge > 0))
    else:
        edge = vol
    planes = _Planes(nz, np.float32)
    if step:
        cur = np.empty(vol.shape, dtype=np.float32)
        if m == "median":
            _background_median(vol, edge, med_factor, cur, planes, t)
        else:
            _background_masked(vol, cur, planes, t)
        zst = np.flatnonzero(planes.nonzero)
    else:
        # the PSF as it is: MATLAB's sums of the planes (values may be negative)
        from petakit5d.matlabcompat import reduce

        cur = vol
        _par(lambda s, e: planes.take(cur, s, e, nonzero=False), nz, t)
        zst = np.flatnonzero(reduce.sum(cur, axis=(1, 2), num_threads=t) != 0)
    pz, py, px = planes.first_max(cur)
    if zst.size == 0:
        raise PSFError(
            _BAD_SUBSCRIPT, "Index exceeds the number of array elements. Index must not exceed 0."
        )
    zs, zt = int(zst[0]), int(zst[-1])
    post = pre = 0
    if pz - zs > nz / 2:
        post = pz - zs - nz // 2 + 3  # round((nz - 1) / 2) is nz // 2
    nz2 = nz + post
    if zt - pz > nz2 / 2:
        pre = zt - pz - nz2 // 2 + 3
    nz2 += pre
    shift = (nz2 // 2 - (pz + pre), ny // 2 - py, nx // 2 - px)
    out = _place(cur, nz2, pre, shift, t)
    cur = None

    ratio = dz_data / dz_psf
    if ratio > 1:
        out = _truncate_z(out, ratio, backend, t)
    ops.clamp(out, lo=0.0, out=out, num_threads=t)  # max(psf, 0)
    if image_in and out.shape[0] == 1:
        out = out[0]
    return out


def _psf_to_matlab(v):
    from petakit5d.core.arrays import to_matlab

    return to_matlab(v)


@register(
    "psf_gen_new",
    Param("psf", None, "Raw PSF, a (z, y, x) volume or a (y, x) image.", convert="volume"),
    Param("dz_psf", None, "z step of the PSF."),
    Param("dz_data", None, "z step of the data the PSF is resampled to."),
    Param("medFactor", None, "Factor of the median subtracted by the 'median' method."),
    Param("PSFGenMethod", None, "Background removal: 'masked' or 'median'."),
    task="R01",
    matlab_out=_psf_to_matlab,
)
def psf_gen(psf, dz_psf, dz_data, med_factor=1.5, psf_gen_method="masked") -> np.ndarray:
    """Clean, center and resample a raw PSF (psf_gen_new): the background removed, the peak at
    the FFT center ``size // 2``, z resampled to the data's step (see the module)."""
    return generate_psf(psf, dz_psf, dz_data, med_factor, psf_gen_method)


# ---------------------------------------------------------------------------------------------
# the PSF crops of RLdecon and XR_decon_data_wrapper
# ---------------------------------------------------------------------------------------------


def _vol3(psf: Any) -> np.ndarray:
    a = np.asarray(psf)
    if a.ndim not in (2, 3):
        raise ValueError("the PSF must be a (z, y, x) volume or a (y, x) image")
    return a[None] if a.ndim == 2 else a


def _support(a: np.ndarray, axis: int) -> np.ndarray:
    """``find(squeeze(sum(psf, other dims)))`` along an axis of the (z, y, x) view: the planes
    with a nonzero sum (MATLAB's sums where a value is negative; else any nonzero value)"""
    other = tuple(i for i in range(3) if i != axis)
    if a.dtype.kind in "fiu" and not (a.dtype.kind != "u" and np.any(a < 0)):
        return np.flatnonzero(np.any(a != 0, axis=other))
    from petakit5d.matlabcompat import reduce

    return np.flatnonzero(reduce.sum(a, axis=other) != 0)


def crop_psf_support(psf) -> np.ndarray:
    """RLdecon's crop of a generated PSF for methods other than OMW (also
    XR_decon_data_wrapper's): along each axis the planes from one before the first to one after
    the last with a nonzero sum, as many cut on both sides as the nearer side allows (the center
    stays). A PSF without nonzero planes raises MATLAB's error. (z, y, x) or (y, x) in and out."""
    a = _vol3(psf)
    lo, hi = [], []
    for ax in range(3):
        s = _support(a, ax)
        if s.size == 0:
            raise PSFError(
                _BAD_SUBSCRIPT,
                "Index exceeds the number of array elements. Index must not exceed 0.",
            )
        c = max(min(int(s[0]), a.shape[ax] - 1 - int(s[-1])) - 1, 0)
        lo.append(c)
        hi.append(a.shape[ax] - c)
    out = np.ascontiguousarray(a[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]])
    return out[0] if np.asarray(psf).ndim == 2 else out


def _msize(shape: Sequence[int]) -> str:
    """mat2str of MATLAB's size [y x z] of a (z, y, x) shape"""
    z, y, x = shape
    return f"[{y} {x} {z}]"


def crop_psf_to_data(psf, data_shape: Sequence[int]) -> np.ndarray:
    """RLdecon's crop of a PSF larger than the data along an axis: the data's length about the
    PSF's center (``floor((p - n) / 2)`` cut before), with MATLAB's warning (identifier '').
    ``data_shape``: (z, y, x) or (y, x)."""
    a = _vol3(psf)
    d = tuple(int(v) for v in data_shape)
    d = (1,) * (3 - len(d)) + d
    if not any(p > n for p, n in zip(a.shape, d, strict=True)):
        return np.asarray(psf)
    w = UserWarning(
        f"The psf size {_msize(a.shape)} is larger than the data size {_msize(d)}, crop it!"
    )
    w.identifier = ""
    warnings.warn(w, stacklevel=2)
    s = [max(0, (p - n) // 2) for p, n in zip(a.shape, d, strict=True)]
    t = [b + min(p, n) for b, p, n in zip(s, a.shape, d, strict=True)]
    out = np.ascontiguousarray(a[s[0] : t[0], s[1] : t[1], s[2] : t[2]])
    return out[0] if np.asarray(psf).ndim == 2 else out


# ---------------------------------------------------------------------------------------------
# decon_readtiff_cached
# ---------------------------------------------------------------------------------------------

_CACHE_SIZE = 4
_CACHE: list[tuple[str, tuple[int, int], np.ndarray]] = []  # the most recent first
_CACHE_LOCK = threading.Lock()


def clear_read_cache() -> None:
    """Forget the files read_tiff_cached keeps (MATLAB's ``clear decon_readtiff_cached``)."""
    with _CACHE_LOCK:
        _CACHE.clear()


def read_tiff_cached(path) -> np.ndarray:
    """``decon_readtiff_cached(fn)``: the tiff file as ``petakit5d.io.read_tiff`` reads it, kept
    for the next calls with the same path (the last four files read): a file whose modification
    time or size changed is read again; a path that is not a file is read without the cache (and
    raises). The array is read-only."""
    from petakit5d import io

    fn = os.fspath(path)
    try:
        st = os.stat(fn)
    except OSError:
        st = None
    if st is None or not stat.S_ISREG(st.st_mode):
        return io.read_tiff(fn)
    info = (st.st_mtime_ns, st.st_size)
    with _CACHE_LOCK:
        for f, i, im in _CACHE:
            if f == fn and i == info:
                return im
    im = io.read_tiff(fn)
    im.setflags(write=False)
    with _CACHE_LOCK:
        _CACHE[:] = [(fn, info, im)] + [e for e in _CACHE if e[0] != fn][: _CACHE_SIZE - 1]
    return im


# ---------------------------------------------------------------------------------------------
# XR_rotate_PSF
# ---------------------------------------------------------------------------------------------


@register(
    "XR_rotate_PSF",
    Param("PSFfile", None, "The PSF tiff file.", py="psf_file"),
    Param("Overwrite", "@islogical", "Rotate again where the result exists."),
    Param("objectiveScan", "@islogical", "Objective scan (zAniso dz / px)."),
    Param("xyPixelSize", "@isscalar", "Pixel size."),
    Param("dz", "@isscalar", "Scan step."),
    Param("SkewAngle", "@isscalar", "Skew angle in degrees."),
    Param("HighPrctile", "@isscalar", "Percentile below which positive values give the fill."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("Save16bit", "@islogical", "Save uint16 (else single).", convert="same"),
    Param("uuid", "@isstr", "Id in the temporary name (empty: a new one)."),
    task="R01",
)
def rotate_psf(
    psf_file,
    *,
    overwrite=False,
    objective_scan=False,
    xy_pixel_size=0.108,
    dz=0.1,
    skew_angle=32.45,
    high_prctile=99,
    reverse=False,
    save_16bit=True,
    uuid="",
) -> None:
    """Rotate a PSF file, as for rotated sample-scan data (``XR_rotate_PSF``): the PSF as double,
    rotated by ``rotate_volume`` (``zAniso = sin(angle) dz / px``, ``dz / px`` for objective
    scans), its zero voxels set to the median of its positive values below their
    ``high_prctile`` percentile, written to ``<folder>/Rotated/<name>.tif`` (uint16, or single
    without ``save_16bit``; first as ``<name>_<uuid>.tif``, then renamed) with parameters.mat
    (v7.3) and parameters.json of the call there. An existing result is kept unless
    ``overwrite``."""
    values = dict(locals())
    from petakit5d.core import dtypes
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.deskew_rotate import rotate_volume
    from petakit5d.deskew_rotate.geometry import _sin
    from petakit5d.io import read_tiff, write_mat, write_tiff
    from petakit5d.matlabcompat import stats
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    uid = uuid if uuid else get_uuid()
    fn = os.fspath(psf_file)
    # MATLAB's fileparts: the folder without its last separator, the name without extension
    folder, base = os.path.split(fn)
    fsname = os.path.splitext(base)[0]
    rot_dir = folder + "/Rotated/"
    os.makedirs(rot_dir, exist_ok=True)
    final = f"{rot_dir}{fsname}.tif"
    if os.path.exists(final) and not overwrite:
        return
    if objective_scan:
        z_aniso = dz / xy_pixel_size
    else:
        theta = skew_angle * np.pi / 180
        z_aniso = _sin(theta) * dz / xy_pixel_size
    im = dtypes.double(read_tiff(fn))
    rt = rotate_volume(im, skew_angle, z_aniso, reverse, crop=True, objective_scan=objective_scan)
    pos = rt[rt > 0]
    hi = stats.prctile(pos, high_prctile)[0]
    med = stats.median(rt[(rt > 0) & (rt < hi)])
    rt = np.where(rt == 0, med, rt)
    tmp = f"{rot_dir}{fsname}_{uid}.tif"
    out = dtypes.uint16(rt) if save_16bit else dtypes.single(rt)
    write_tiff(np.ascontiguousarray(out), tmp)
    os.replace(tmp, final)
    pr = {k: to_plain(v) for k, v in matlab_results(rotate_psf, values).items()}
    write_mat(rot_dir + "/parameters.mat", {"pr": pr})
    write_json_file(pr, rot_dir + "/parameters.json")
