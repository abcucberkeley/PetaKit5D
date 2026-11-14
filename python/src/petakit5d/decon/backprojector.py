"""The OTF masked Wiener (OMW) back projector (``omw_backprojector_generation``) and its files in
the psfgen folder (``omw_backprojector_files``).

Arrays are the port's: a PSF is a ``(z, y, x)`` volume or a ``(y, x)`` image (one z plane, as
MATLAB's 2D PSF; a ``(1, y, x)`` volume is that image too), C-contiguous. Results equal MATLAB's
bit for bit (omw_backprojector_generation.m with omw_backprojector_mex at the oracle commit) on
the exact FFT backend:

- **omw_backprojector** (``generate_backprojector`` with the FFT backend and thread count):
  the PSF in single, divided by its fixed-order double sum (``decon_sum_double``, the quotient
  in double rounded to single); its OTF (``decon_psf2otf``), ``|OTF|`` (MATLAB's abs of complex
  single: msun's hypotf in single, measured) and its ``fftshift``. The threshold: the |OTF|
  values sorted descending in double, the first whose running sum (in order) exceeds
  ``OTFCumThresh`` times their fixed-order sum, over the DC value (single), at least
  ``minIntThrsh``. The mask: the shifted |OTF| above DC times the threshold, opened,
  area-opened (``OTFAreaThresh``, 26-connected components), closed, opened and area-opened
  again (``strel('sphere', 2)``).
- **Skewed PSFs** (``skewed`` None: detected as MATLAB detects it, from the 8-connected
  components of the mask's projection along y and where its z profile peaks): the
  ``bwlabeln`` labels, the 3 largest components where there are more than three, labels 3, 2
  and 1 stacked along z; otherwise the largest component. Then each z, y and x slice with
  voxels (in that order) is replaced by its ``bwconvhull`` (the mex's hulls), and the edge is
  the mask without its erosion by ``strel('sphere', 1)``.
- **The window** (``hannWinBounds = [l, u]``): ``win(x) = (x < l) + cos(pi (x - l) / 2 /
  (u - l))^2 (l <= x < u)`` of each voxel's radius over ``rn``, a combination of the radii of
  the 2 edge points nearest to its direction (``cart2sph``, ``sph2cart`` and ``knnsearch`` as
  MATLAB computes them). The weights are kept as MATLAB indexes them: the distances of the
  voxels numbered by the edge points' indices in MATLAB's linear order, not the voxel's own
  (speed/reports/omw.md). Where the mex's bounds decide a voxel (0 < l < u < Inf, at least two
  edge points, positive weights and radii) the window is 1 or 0 without its nearest points, and
  the decon kernels compute the others (fdlibm's and msun's functions in C, the nearest edge
  points on a grid of cells; ``_window_fast_numpy`` is the same in NumPy); otherwise every
  voxel gets MATLAB's full computation. A skewed window and mask are folded back by adding
  their z thirds.
- The back projector: ``fftshift(real(ifftn(ifftshift(window) .* conj(OTF) ./ (|OTF|.^2 +
  alpha))))`` with MATLAB's mixed-class arithmetic (measured: ``conj`` gives ``0 - b``, the
  double window times a complex single is ``times_complex``).

The outputs (``BackProjector``) are MATLAB's four: the back projector (single), the masked
Wiener filter (complex single), the shifted |OTF| (single) and the mask (bool; double counts
for skewed PSFs). Nothing depends on the thread count. Known differences:

- **The permissive FFT backend** (pocketfft, the default; decision record 0006): the back
  projector and the masked Wiener filter within a relative L2 error of 1e-5 and 1e-5 of the
  largest value (accepted; measured up to 1.9e-6: the Wiener division amplifies the OTF's last
  bits where |OTF| is small), the |OTF| within 0006, the masks identical (in every case
  measured; a voxel at the threshold could move).
- **NaN signs:** with degenerate bounds (``l == u``) the window and the back projector are NaN;
  the inverse FFT propagates NaN signs differently in MATLAB's FFTW and the exact backend, so
  only the NaN positions are compared there.
- ``times_complex`` differs from MATLAB where both parts of a value are non-finite (one of them
  infinite) and the window is +-0 (only a non-finite OTF has such values).

**omw_backprojector_files**:
``<psfgen>/<name>_back_projector_alpha_<%0.6f>_otf_cum_thresh_<%0.6f>.tif`` is read
(``read_tiff_cached``) if it exists, else made, with ``<name>_OTF.tif`` (the shifted |OTF|) and
``<name>_OTF_mask_otf_cum_thresh_<%0.6f>.tif`` (uint8), each written under a name with a uuid
and renamed. The back projector is published last (MATLAB writes it first), so that the other
files exist once it does, and a back projector another writer published meanwhile is kept. The
figure of the mask's outline on the OTF (``<...>_figure.png``) is drawn with matplotlib when it
is installed; the data files don't depend on it.
"""

from __future__ import annotations

import math
import os
import threading
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any, NamedTuple

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "BackProjector",
    "BackProjectorError",
    "backprojector_file_names",
    "generate_backprojector",
    "omw_backprojector",
    "omw_backprojector_files",
    "window_function",
]

_EMPTY_PSF = "A valid nonzero PSF needs to be provided!"
_EMPTY_MASK = "The OTF mask is empty, check the PSF and OTF-related parameters!"
_PI = 3.141592653589793
_CHUNK = 1 << 16  # elements per piece of a threaded element-wise pass


class BackProjectorError(ValueError):
    """an error MATLAB raises, with its identifier ('' for ``error(message)``)"""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


class BackProjector(NamedTuple):
    """omw_backprojector_generation's outputs, in the port's arrays: ``b_omw`` the back projector
    (single), ``otf_bp_omw`` the masked Wiener filter (complex single; single where the OTF is
    real), ``abs_otf_c`` the fftshifted |OTF| (single), ``otf_mask`` the OTF mask (bool; for
    skewed PSFs its z thirds added, double)."""

    b_omw: np.ndarray
    otf_bp_omw: np.ndarray
    abs_otf_c: np.ndarray
    otf_mask: np.ndarray


# ---------------------------------------------------------------------------------------------
# threads
# ---------------------------------------------------------------------------------------------


def _threads(threads: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(threads)


_POOLS: dict[int, ThreadPoolExecutor] = {}
_POOLS_LOCK = threading.Lock()


def _pool(k: int) -> ThreadPoolExecutor:
    """a pool of k threads, kept for the next passes (a pool per thread count)"""
    with _POOLS_LOCK:
        ex = _POOLS.get(k)
        if ex is None:
            ex = _POOLS[k] = ThreadPoolExecutor(k, thread_name_prefix="pk5d-bp")
        return ex


def _par(fn: Callable[[int, int], None], n: int, nt: int, chunk: int = _CHUNK) -> None:
    """``fn(start, stop)`` on pieces of ``range(n)`` of about ``chunk`` elements, on up to ``nt``
    threads (NumPy's loops release the interpreter lock)"""
    pieces = max(1, -(-n // chunk))
    k = min(nt, pieces)
    edges = [n * i // pieces for i in range(pieces + 1)]
    if k <= 1 or threading.current_thread().name.startswith("pk5d-bp"):
        for i in range(pieces):
            fn(edges[i], edges[i + 1])
        return
    ex = _pool(k)
    for f in [ex.submit(fn, edges[i], edges[i + 1]) for i in range(pieces)]:
        f.result()


# ---------------------------------------------------------------------------------------------
# the OTF and its threshold
# ---------------------------------------------------------------------------------------------


def _single(x: Any, nt: int) -> np.ndarray:
    """MATLAB's single(x), C-contiguous"""
    a = np.asarray(x)
    if a.dtype == np.float32:
        return np.ascontiguousarray(a)
    if a.dtype.kind in "biu" and a.dtype.itemsize < 8 or a.dtype == np.float64:
        from petakit5d.kernels import ops

        return ops.convert(np.ascontiguousarray(a), "single", num_threads=nt)
    from petakit5d.core import dtypes

    return np.ascontiguousarray(dtypes.single(a))


def _divide_single(a: np.ndarray, s: float, nt: int) -> np.ndarray:
    """``a ./ s`` of single a and a double scalar s: the quotient in double, rounded to single"""
    from petakit5d.core import dtypes

    out = np.empty(a.shape, np.float32)
    fa, fo = a.reshape(-1), out.reshape(-1)

    def run(i0: int, i1: int) -> None:
        q = fa[i0:i1].astype(np.float64)
        np.divide(q, s, out=q)
        fo[i0:i1] = dtypes.x86_nan(q, fa[i0:i1])

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        _par(run, fa.size, nt)
    return out


_HI = np.int32(-4096)  # 0xfffff000: the high 12 bits of a single's significand and its exponent


def _f32(w: np.ndarray) -> np.ndarray:
    return w.view(np.float32)


def _hypotf_medium(ha: np.ndarray, hb: np.ndarray) -> np.ndarray:
    """msun's hypotf for medium a >= b given as their bits (signs cleared): sqrt(a^2 + b^2)
    from the products of their high parts, in single"""
    a, b = _f32(ha), _f32(hb)
    w = a - b
    t1 = _f32(ha & _HI)
    far = np.sqrt(t1 * t1 - (b * (-b) - (a - t1) * (a + t1)))
    y1 = _f32(hb & _HI)
    t1 = _f32((ha + 0x00800000) & _HI)
    near = np.sqrt(t1 * y1 - (w * (-w) - (t1 * (b - y1) + ((a + a) - t1) * b)))
    return np.where(w > b, far, near)


def _hypotf_rest(ha: np.ndarray, hb: np.ndarray) -> np.ndarray:
    """msun's hypotf where a / b > 2^30, a is Inf or NaN, or a or b is scaled first"""
    a, b = _f32(ha).copy(), _f32(hb).copy()
    direct = (ha - hb) > 0x0F000000
    special = ~direct & (ha >= 0x7F800000)
    big = ~direct & ~special & (ha > 0x58800000)
    small = ~direct & ~special & (hb < 0x26800000)
    zero = small & (hb == 0)
    sub = small & ~zero & (hb <= 0x007FFFFF)
    tiny = small & ~zero & ~sub
    k = np.zeros(ha.shape, np.int32)
    ha, hb = ha.copy(), hb.copy()
    # a > 2^50: both scaled by 2^-68; b < 2^-50: by 2^68, or for a subnormal b by 2^126 (the
    # high words of a and b are not read again there)
    ha[big] -= 0x22000000
    hb[big] -= 0x22000000
    k[big] = 68
    ha[tiny] += 0x22000000
    hb[tiny] += 0x22000000
    k[tiny] = -68
    a = np.where(big | tiny, _f32(ha), a)
    b = np.where(big | tiny, _f32(hb), b)
    a[sub] *= np.float32(2.0**126)
    b[sub] *= np.float32(2.0**126)
    k[sub] = -126
    w = a - b
    t1 = _f32(ha & _HI)
    far = np.sqrt(t1 * t1 - (b * (-b) - (a - t1) * (a + t1)))
    y1 = _f32(hb & _HI)
    t1 = _f32((ha + 0x00800000) & _HI)
    near = np.sqrt(t1 * y1 - (w * (-w) - (t1 * (b - y1) + ((a + a) - t1) * b)))
    h = np.where(w > b, far, near)
    h = np.where(k != 0, _f32(np.int32(0x3F800000) + (k << 23)) * h, h)
    h = np.where(direct, _f32(ha) + _f32(hb), h)
    h = np.where(zero, _f32(ha), h)
    # MATLAB's NaN for a NaN part (also next to Inf, where hypotf gives Inf), else Inf
    nan = special & (ha > 0x7F800000)
    h = np.where(special & ~nan, np.float32(np.inf), h)
    return np.where(nan, _f32(np.int32(0x7FC00000)), h)


def _ext():
    from petakit5d._ext import load

    return load("decon")


def _abs(otf: np.ndarray, nt: int) -> np.ndarray:
    """MATLAB's ``abs`` of the single OTF: for complex values msun's ``hypotf`` (e_hypotf.c in
    single precision, with MATLAB's NaN for a NaN part; MATLAB's bits, measured), in the decon
    kernels (``_abs_numpy`` writes it out in NumPy); for real values the magnitude"""
    if not np.iscomplexobj(otf):
        return np.abs(otf)
    c = np.ascontiguousarray(otf, dtype=np.complex64).reshape(-1).view(np.float32)
    return _ext().omw_abs(c, nt).reshape(otf.shape)


def _abs_numpy(otf: np.ndarray, nt: int) -> np.ndarray:
    """``_abs`` of complex single values in NumPy (the reference of the kernel)"""
    out = np.empty(otf.shape, np.float32)
    bits = np.ascontiguousarray(otf).reshape(-1).view(np.int32).reshape(-1, 2)
    fo = out.reshape(-1)

    def run(i0: int, i1: int) -> None:
        hx = bits[i0:i1, 0] & 0x7FFFFFFF
        hy = bits[i0:i1, 1] & 0x7FFFFFFF
        ha, hb = np.maximum(hx, hy), np.minimum(hx, hy)
        h = _hypotf_medium(ha, hb)
        rest = np.flatnonzero(((ha - hb) > 0x0F000000) | (ha > 0x58800000) | (hb < 0x26800000))
        if rest.size:
            h[rest] = _hypotf_rest(ha[rest], hb[rest])
        fo[i0:i1] = h

    with np.errstate(all="ignore"):
        _par(run, fo.size, nt)
    return out


def _sorted_descending(a: np.ndarray, nt: int) -> np.ndarray:
    """``sort(double(a(:)), 'descend')``: the mex's radix sort, MATLAB's sort where a has NaN
    (first when descending)"""
    from petakit5d.kernels import decon as kd

    v = kd.omw_sortdesc(np.ascontiguousarray(a, np.float32), threads=nt)
    if v is not None:
        return v
    from petakit5d.matlabcompat import reduce

    return np.asarray(reduce.sort(a.astype(np.float64).reshape(-1), descend=True))


def _first_above(v: np.ndarray, t: float) -> int:
    """``find(cumsum(v) > t, 1, 'first')`` (0-based; -1 if none): the running sum in order, in
    pieces that carry the sum so far"""
    acc = 0.0
    for i0 in range(0, v.size, 1 << 20):
        c = v[i0 : i0 + (1 << 20)].copy()
        c[0] = acc + c[0]
        np.cumsum(c, out=c)
        hit = np.flatnonzero(c > t)
        if hit.size:
            return i0 + int(hit[0])
        acc = float(c[-1])
    return -1


def _otf_threshold(abs_otf, dc, otf_thresh, min_int_thrsh, nt, verbose) -> np.float32:
    """the OTF threshold over DC: the |OTF| value where the running sum of the sorted values
    first exceeds ``otf_thresh`` of their sum, over DC (single), at least ``min_int_thrsh``"""
    from petakit5d.core import text
    from petakit5d.kernels import decon as kd

    vals = _sorted_descending(abs_otf, nt)
    total = kd.sum_double(vals, threads=nt)
    t = _first_above(vals, total * float(otf_thresh))
    if t < 0:
        # dc_thresh is empty: the comparison with the 3D |OTF| fails
        if verbose:
            print("OTF threshold: ", flush=True)
        raise BackProjectorError(
            "MATLAB:sizeDimensionsMustMatch", "Arrays have incompatible sizes for this operation."
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.float32(float(vals[t]) / float(dc))
    m = np.float32(min_int_thrsh)
    if np.isnan(q) or (not np.isnan(m) and m > q):
        q = m
    if verbose:
        print(text.sprintf("OTF threshold: %0.2d\n", float(q)), end="", flush=True)
    return q


# ---------------------------------------------------------------------------------------------
# the mask
# ---------------------------------------------------------------------------------------------


def _area_open(bw: np.ndarray, p: Any, nt: int) -> np.ndarray:
    """``bwareaopen(bw, p)``: the mex's components for 3D masks (MATLAB's ndims 3), bwareaopen
    for images"""
    if bw.ndim == 3:
        from petakit5d.kernels import decon as kd

        pv = float(p)
        if pv >= 0 and pv == math.floor(pv) and pv <= 9e15:
            return kd.omw_components(bw, "minsize", pv, threads=nt)
    from petakit5d.matlabcompat.label import bwareaopen

    return bwareaopen(bw, p, threads=nt)


def _largest_components(bw: np.ndarray, k: int, nt: int) -> np.ndarray:
    """the k largest 26-connected components (ties to the first in bwconncomp's order)"""
    if bw.ndim == 3:
        from petakit5d.kernels import decon as kd

        return kd.omw_components(bw, "largest", k, threads=nt)
    from petakit5d.matlabcompat.label import bwconncomp

    cc = bwconncomp(bw, 26, threads=nt)
    sizes = np.array([p.size for p in cc.pixel_idx_list], dtype=np.int64)
    keep = np.argsort(-sizes, kind="stable")[: max(int(k), 0)]
    out = np.zeros(bw.size, dtype=bool)
    for i in keep:
        out[cc.pixel_idx_list[i]] = True
    return out.reshape(bw.shape)


def _otf_mask(abs_c: np.ndarray, thr: np.float32, area: Any, nt: int) -> np.ndarray:
    from petakit5d.kernels import ops
    from petakit5d.matlabcompat.morphology import imclose, imopen, strel

    se = strel("sphere", 2)
    m = ops.compare(abs_c, "gt", float(thr), num_threads=nt)
    m = imopen(m, se, threads=nt)
    m = _area_open(m, area, nt)
    m = imclose(m, se, threads=nt)
    m = imopen(m, se, threads=nt)
    m = _area_open(m, area, nt)
    return m


def _is_skewed(m: np.ndarray, nt: int) -> bool:
    """MATLAB's detection: the 8-connected components of ``squeeze(sum(mask, 1)) > 0``; with
    more than one, skewed where every entry of its sum along the first dimension is positive and
    the first maximum is more than 1 from the middle"""
    from petakit5d.matlabcompat.label import bwconncomp

    proj = m.any(axis=-2)  # 3D: (nz, nx); image: (nx,)
    if m.ndim == 2:
        img = proj.reshape(1, -1)  # MATLAB's 1 x nx row
    elif proj.shape[1] == 1:
        img = proj.reshape(-1, 1)  # squeeze of 1 x 1 x nz: an nz x 1 column
    else:
        img = np.ascontiguousarray(proj.T)  # nx x nz
    if bwconncomp(img, 8, threads=nt).num_objects <= 1:
        return False
    line = img.sum(axis=0)
    pind = int(np.argmax(line)) + 1
    return bool(np.all(line > 0)) and abs(pind - (line.size + 1) / 2) > 1


def _three_lobes(m: np.ndarray, nt: int) -> np.ndarray:
    """the skewed mask: the bwlabeln labels (the 3 largest components where there are more than
    three labels besides the background), labels 3, 2 and 1 stacked along z"""
    from petakit5d.matlabcompat.label import bwlabeln

    lab, n = bwlabeln(m, threads=nt)
    nl = n + int(not m.all())
    if nl > 4:
        m = _largest_components(m, 3, nt)
        lab, n = bwlabeln(m, threads=nt)
    elif nl == 2:
        lab = (lab + 1) * m
    lab3 = lab.reshape((1,) + lab.shape) if lab.ndim == 2 else lab
    k = lab3.shape[0]
    out = np.empty((3 * k,) + lab3.shape[1:], dtype=bool)

    def run(i0: int, i1: int) -> None:
        for j, label in enumerate((3.0, 2.0, 1.0)):
            np.equal(lab3[i0:i1], label, out=out[j * k + i0 : j * k + i1])

    _par(run, k, nt, chunk=max(1, (1 << 16) // max(1, lab3[0].size)))
    return out


def _slice_image(m3: np.ndarray, axis: int, k: int) -> tuple[np.ndarray, Callable]:
    """MATLAB's squeezed slice of the (nz, ny, nx) mask (MATLAB's ny x nx x nz) at k along axis
    (0: ``A(:, :, z)``, 1: ``squeeze(A(y, :, :))``, 2: ``squeeze(A(:, x, :))``) as a (rows,
    columns) image, and the function that writes an image of that shape back"""
    nz, ny, nx = m3.shape
    if axis == 0:
        img = m3[k]
    elif axis == 1:
        img = (
            m3[0, k, :].reshape(1, nx)
            if nz == 1
            else (m3[:, k, 0].reshape(nz, 1) if nx == 1 else m3[:, k, :].T)
        )
    else:
        img = (
            m3[0, :, k].reshape(ny, 1)
            if nz == 1
            else (m3[:, 0, k].reshape(nz, 1) if ny == 1 else m3[:, :, k].T)
        )
    shape = img.shape

    def put(h: np.ndarray) -> None:
        h = h.reshape(shape)
        if axis == 0:
            m3[k] = h
        elif axis == 1:
            if nz == 1:
                m3[0, k, :] = h[0]
            elif nx == 1:
                m3[:, k, 0] = h[:, 0]
            else:
                m3[:, k, :] = h.T
        elif nz == 1:
            m3[0, :, k] = h[:, 0]
        elif ny == 1:
            m3[:, 0, k] = h[:, 0]
        else:
            m3[:, :, k] = h.T

    return np.ascontiguousarray(img), put


def _hulls(m: np.ndarray, nt: int) -> np.ndarray:
    """each z, then y, then x slice with voxels replaced by its ``bwconvhull(slice, 'union')``
    (the mex's hulls where every size is above 1, MATLAB's loops of bwconvhull otherwise)"""
    m3 = m if m.ndim == 3 else m.reshape((1,) + m.shape)
    zi = np.flatnonzero(m3.any(axis=(1, 2)))
    yi = np.flatnonzero(m3.any(axis=(0, 2)))
    xi = np.flatnonzero(m3.any(axis=(0, 1)))
    if min(m3.shape) > 1:
        from petakit5d.kernels import decon as kd

        return kd.omw_convhull(m3, zi, yi, xi, threads=nt)
    from petakit5d.matlabcompat.morphology import bwconvhull

    out = np.array(m3, order="C", copy=True)
    for axis, inds in ((0, zi), (1, yi), (2, xi)):
        for k in inds:
            img, put = _slice_image(out, axis, int(k))
            put(bwconvhull(img, "union"))
    return out.reshape(m.shape)


# ---------------------------------------------------------------------------------------------
# the window
# ---------------------------------------------------------------------------------------------


def window_function(x, lo: float, hi: float, *, threads: int = 1) -> np.ndarray:
    """``(x < l) + cos(pi * (x - l) / 2 / (u - l)) .^ 2 .* (x >= l & x < u)`` (double), with
    MATLAB's cos, and NaN as MATLAB has it"""
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat import elementary

    x = np.asarray(x, dtype=np.float64)
    lo, hi = float(lo), float(hi)
    below = (x < lo).astype(np.float64)
    taper = (x >= lo) & (x < hi)
    with np.errstate(all="ignore"):
        d = dtypes.x86_nan(x - lo, x)
        t = dtypes.x86_nan((_PI * d) / 2.0, d)
        t = dtypes.x86_nan(t / (hi - lo), t)
    # elsewhere the product is +0, unless the cosine is NaN (x or the argument not finite)
    need = np.flatnonzero(taper | ~np.isfinite(t))
    out = below
    if need.size:
        tn, c = t.reshape(-1)[need], np.empty(need.size)

        def run(i0: int, i1: int) -> None:
            c[i0:i1] = elementary.cos(tn[i0:i1])

        _par(run, need.size, threads, chunk=1 << 15)
        out = below.copy()
        fo, fb, ft = out.reshape(-1), below.reshape(-1), taper.reshape(-1)
        fo[need] = (c * c) * ft[need] + fb[need]
    return out


class _Directions:
    """``[az, el, r] = cart2sph(x - c(2), y - c(1), z - c(3)); [Xn, Yn, Zn] = sph2cart(az, el,
    1)`` of the voxels of a grid of MATLAB size psz with center c = (psz + 1) / 2, as MATLAB
    computes them elementwise (matlabcompat.neighbors). hypot(x, y), az and its cosine and sine
    depend on the voxel's (y, x) only, el and r on hypot(x, y) and z: each is computed once per
    distinct argument (the same functions of the same values, so the same bits)."""

    def __init__(self, psz: Sequence[int]):
        from petakit5d.matlabcompat.elementary import atan2, fdlibm_cos, fdlibm_sin, hypot

        ny, nx, nz = (int(v) for v in psz)
        self.psz = (ny, nx, nz)
        # (nx, ny) arrays: the flat index x * ny + y is MATLAB's linear index in a plane
        dx = (np.arange(1, nx + 1, dtype=np.float64) - (nx + 1) / 2).reshape(-1, 1)
        dy = (np.arange(1, ny + 1, dtype=np.float64) - (ny + 1) / 2).reshape(1, -1)
        dx, dy = np.broadcast_arrays(dx, dy)
        h = hypot(dx, dy)
        az = atan2(dy, dx)
        self.caz = np.ascontiguousarray(fdlibm_cos(az)).reshape(-1)
        self.saz = np.ascontiguousarray(fdlibm_sin(az)).reshape(-1)
        self.h, inv = np.unique(h.reshape(-1), return_inverse=True)
        self.hidx = inv.reshape(-1).astype(np.int64)  # by MATLAB's (y, x): y + ny * x
        self.dz = np.arange(1, nz + 1, dtype=np.float64) - (nz + 1) / 2

    def __call__(self, t: np.ndarray, nt: int = 1):
        """(Yn, Xn, Zn, r) of the voxels at the 0-based MATLAB linear indices t"""
        from petakit5d.matlabcompat.elementary import atan2, fdlibm_cos, fdlibm_sin, hypot

        ny, nx, nz = self.psz
        nh = self.h.size
        yx = t % (ny * nx)
        z = t // (ny * nx)
        # dz and -dz: the same r and cos(el), el and sin(el) of the other sign (atan2 is odd in
        # y for x >= 0, sin odd and cos even, in fdlibm), so the tables hold z >= center
        neg = 2 * z < nz - 1
        zk = np.where(neg, nz - 1 - z, z) - (nz - 1) // 2
        key = zk * nh + self.hidx[yx]
        used = np.zeros((nz - (nz - 1) // 2) * nh, dtype=bool)
        used[key] = True
        u = np.flatnonzero(used)
        r_u, ce_u, se_u = (np.empty(u.size) for _ in range(3))
        dz = self.dz[(nz - 1) // 2 :]

        def run(i0: int, i1: int) -> None:
            hv, zv = self.h[u[i0:i1] % nh], dz[u[i0:i1] // nh]
            el = atan2(zv, hv)
            r_u[i0:i1], ce_u[i0:i1], se_u[i0:i1] = hypot(hv, zv), fdlibm_cos(el), fdlibm_sin(el)

        _par(run, u.size, nt, chunk=1 << 15)
        pos = np.zeros(used.size, dtype=np.int64)
        pos[u] = np.arange(u.size)
        p = pos[key]
        ce = ce_u[p]
        se = se_u[p]
        np.negative(se, out=se, where=neg)
        # sph2cart with radius 1: z = 1 .* sin(el), x = (1 .* cos(el)) .* cos(az), ...
        return ce * self.saz[yx], ce * self.caz[yx], se, r_u[p]


def _directions(dx: np.ndarray, dy: np.ndarray, dz: np.ndarray, nt: int):
    """``[az, el, r] = cart2sph(dx, dy, dz); [xn, yn, zn] = sph2cart(az, el, 1)`` elementwise:
    (yn, xn, zn, r)"""
    from petakit5d.matlabcompat.neighbors import cart2sph, sph2cart

    az, el, r = cart2sph(dx, dy, dz)
    x, y, z = sph2cart(az, el, np.ones(dx.size))
    return y, x, z, r


def _offsets(t: np.ndarray, psz: Sequence[int]):
    """the offsets (dx, dy, dz) from the center (psz + 1) / 2 of the voxels at the 0-based MATLAB
    linear indices t, in MATLAB's 1-based coordinates (exact)"""
    ny, nx, nz = psz
    y = t % ny
    x = (t // ny) % nx
    z = t // (ny * nx)
    return (
        (x + 1).astype(np.float64) - (nx + 1) / 2,
        (y + 1).astype(np.float64) - (ny + 1) / 2,
        (z + 1).astype(np.float64) - (nz + 1) / 2,
    )


def _knn2(points: np.ndarray, queries: np.ndarray, nt: int):
    """``[idx, d] = knnsearch(points, queries, 'K', 2)`` as matlabcompat.neighbors.knnsearch
    gives it (squared distances summed over the columns in order, ties to the lower index,
    0-based), with the k-d tree's candidates on ``nt`` threads (knnsearch's too)"""
    from petakit5d.matlabcompat.neighbors import knnsearch

    n, nq = points.shape[0], queries.shape[0]
    if n < 64 or nq * n <= (1 << 22) or not np.all(np.isfinite(points)):
        return knnsearch(points, queries, 2, num_threads=nt)
    from scipy.spatial import cKDTree

    k = 2
    tree = cKDTree(points)
    idx = np.empty((nq, k), np.int64)
    dk = np.empty((nq, k))

    def rank(rows: np.ndarray, m: int) -> np.ndarray:
        """the rows whose 2 nearest points are among the tree's m nearest (every point within
        the 2nd distance and a margin), ranked exactly; returns the other rows"""
        dist, cand = tree.query(queries[rows], k=m, workers=nt)
        complete = (dist[:, m - 1] > dist[:, k - 1] * (1 + 1e-9) + 1e-300) | (m >= n)

        def run(i0: int, i1: int) -> None:
            c, q = cand[i0:i1], queries[rows[i0:i1]]
            df = q[:, None, :] - points[c]
            df *= df
            # summed over the columns in order (0 + the first is the first)
            d2 = df[:, :, 0] + df[:, :, 1]
            for j in range(2, points.shape[1]):
                d2 += df[:, :, j]
            ok = complete[i0:i1]
            done = rows[i0:i1]
            # in the tree's order the exact squared distances mostly increase strictly: the
            # first two are the nearest; other rows are sorted by squared distance, then index
            first = (d2[:, 0] < d2[:, 1]) & (d2[:, 1] < d2[:, 2:].min(axis=1))
            idx[done[ok & first]] = c[ok & first, :k]
            dk[done[ok & first]] = np.sqrt(d2[ok & first, :k])
            other = np.flatnonzero(ok & ~first)
            if other.size:
                co, do = c[other], d2[other]
                pos = np.lexsort((co, do), axis=-1)[:, :k]
                idx[done[other]] = np.take_along_axis(co, pos, axis=1)
                dk[done[other]] = np.sqrt(np.take_along_axis(do, pos, axis=1))

        _par(run, rows.size, nt, chunk=1 << 15)
        return rows[~complete]

    rest = rank(np.arange(nq), 3)
    if rest.size:
        rest = rank(rest, min(64, n))
    for r in rest:
        radius = tree.query(queries[r], k=k)[0][k - 1] * (1 + 1e-9) + 1e-300
        ball = np.asarray(sorted(tree.query_ball_point(queries[r], radius)), np.int64)
        i, dd = knnsearch(points[ball], queries[r], k, num_threads=nt)
        idx[r] = ball[i[0]]
        dk[r] = dd[0]
    return idx, dk


def _window_fast(edge: np.ndarray, psz, lo: float, hi: float, nt: int) -> np.ndarray | None:
    """the window with the mex's bounds (omw_window_with_mex), as a (z, y, x) float64 array, or
    None where MATLAB computes it in full (fewer than two edge points, a zero weight or radius):
    the decon kernels' window values (``_window_fast_numpy`` is the same in NumPy)"""
    from petakit5d.core.arrays import matlab_order

    e_t = np.flatnonzero(matlab_order(edge))  # find(OTF_mask_edge), 0-based
    if e_t.size < 2:
        return None
    w = _ext().omw_window_values(list(psz), e_t.astype(np.int64), lo, hi, nt)
    return None if w is None else w.reshape(edge.shape)


def _window_fast_numpy(edge: np.ndarray, psz, lo: float, hi: float, nt: int):
    """``_window_fast`` in NumPy, on the kernels' bounds (the reference of the kernel)"""
    from petakit5d.core.arrays import matlab_order, zyx_to_matlab_index
    from petakit5d.kernels import decon as kd

    e_t = np.flatnonzero(matlab_order(edge))
    ne = e_t.size
    if ne < 2:
        return None
    directions = _Directions(psz)
    eyn, exn, ezn, rn = directions(e_t, nt)
    shape = edge.shape
    mask, u_zyx = kd.omw_window(shape, np.stack([ezn, eyn, exn], axis=1), rn, (lo, hi), threads=nt)
    u_t = np.asarray(zyx_to_matlab_index(u_zyx, shape), np.int64)
    # the first nE voxels (their distances weigh the radii) and the voxels to compute
    q_t = np.concatenate([np.arange(ne, dtype=np.int64), u_t])
    qyn, qxn, qzn, rs = directions(q_t, nt)
    pts = np.stack([eyn, exn, ezn], axis=1)
    inds, d = _knn2(pts, np.stack([qyn, qxn, qzn], axis=1), nt)
    dmat = d[:ne]
    if not (dmat[:, 1].min() > 0 and rn.min() > 0):
        return None
    inds = inds[ne:]
    d1 = dmat[inds[:, 0], 0]
    d2 = dmat[inds[:, 1], 1]
    r = (rn[inds[:, 0]] * d2 + rn[inds[:, 1]] * d1) / (d1 + d2)
    mask[u_zyx] = window_function(rs[ne:] / r, lo, hi, threads=nt)
    return mask


def _window_full(edge: np.ndarray, psz, lo: float, hi: float, nt: int) -> np.ndarray:
    """the window of MATLAB's code without the mex: every voxel's nearest edge points"""
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import matlab_order

    yn, xn, zn, rs = _Directions(psz)(np.arange(edge.size, dtype=np.int64), nt)
    sel = matlab_order(edge).reshape(-1)
    pts = np.stack([yn[sel], xn[sel], zn[sel]], axis=1)
    inds, d = _knn2(pts, np.stack([yn, xn, zn], axis=1), nt)
    if inds.shape[1] < 2:
        raise BackProjectorError("MATLAB:badsubscript", "Index in position 2 exceeds array bounds.")
    rn = rs[sel]
    d1 = d[inds[:, 0], 0]
    d2 = d[inds[:, 1], 1]
    with np.errstate(all="ignore"):
        a = dtypes.x86_nan(rn[inds[:, 0]] * d2, d2)
        b = dtypes.x86_nan(rn[inds[:, 1]] * d1, d1)
        num = dtypes.x86_nan(a + b, a, b)
        den = d1 + d2
        r = dtypes.x86_nan(num / den, num, den)
        x = dtypes.x86_nan(rs / r, r)
    w = window_function(x, lo, hi, threads=nt)
    ny, nx, nz = psz
    # MATLAB's linear order is the C order of the (z, x, y) view
    return np.ascontiguousarray(w.reshape(nz, nx, ny).swapaxes(1, 2)).reshape(edge.shape)


def _window(mask: np.ndarray, bounds, nt: int) -> np.ndarray:
    from petakit5d.matlabcompat.morphology import imerode, strel

    edge = mask & ~imerode(mask, strel("sphere", 1), threads=nt)
    m3 = edge if edge.ndim == 3 else edge.reshape((1,) + edge.shape)
    psz = (m3.shape[1], m3.shape[2], m3.shape[0])
    lo, hi = float(bounds[0]), float(bounds[1])
    w = None
    if lo > 0 and lo < hi and hi < math.inf:
        w = _window_fast(m3, psz, lo, hi, nt)
    if w is None:
        w = _window_full(m3, psz, lo, hi, nt)
    return w.reshape(mask.shape)


def _fold(a: np.ndarray, nt: int = 1) -> np.ndarray:
    """``sum(reshape(a, ny, nx, nz / 3, 3), 4)``: the z thirds added in order (double; a mask
    gives its counts), plane by plane on threads"""
    k = a.shape[0] // 3
    p = np.ascontiguousarray(a).reshape(3, k, -1)
    out = np.empty((k, p.shape[2]))

    def run(i0: int, i1: int) -> None:
        if a.dtype == np.bool_:
            out[i0:i1] = p[0, i0:i1].view(np.uint8) + p[1, i0:i1].view(np.uint8) + p[2, i0:i1]
        else:
            np.add(p[0, i0:i1], p[1, i0:i1], out=out[i0:i1])
            out[i0:i1] += p[2, i0:i1]

    _par(run, k, nt, chunk=max(1, (1 << 16) // max(1, p.shape[2])))
    return out.reshape((k,) + a.shape[1:])


# ---------------------------------------------------------------------------------------------
# the Wiener filter and the back projector
# ---------------------------------------------------------------------------------------------


def _masked_wiener(otf, abs_otf, om, alpha: float, nt: int):
    """``om .* (conj(OTF) ./ (|OTF|.^2 + alpha))`` as MATLAB computes it (measured): ``.^ 2`` a
    single product, ``+ alpha`` in double rounded to single; ``conj`` (``conj_imag``); each part
    divided by the real denominator in single; the double window times it as ``times_complex``
    (a real OTF: the product in double, rounded). ``om`` None: the Wiener filter alone."""
    from petakit5d.core import dtypes

    cplx = np.iscomplexobj(otf)
    out = np.empty(otf.shape, np.complex64 if cplx else np.float32)
    fa = abs_otf.reshape(-1)
    fm = None if om is None else om.reshape(-1)
    if cplx:
        src = np.ascontiguousarray(otf).reshape(-1).view(np.float32).reshape(-1, 2)
        dst = out.reshape(-1).view(np.float32).reshape(-1, 2)
    else:
        src = np.ascontiguousarray(otf).reshape(-1, 1)
        dst = out.reshape(-1, 1)

    def run(i0: int, i1: int) -> None:
        a = fa[i0:i1]
        den = np.multiply(a, a).astype(np.float64)
        den += alpha
        den32 = dtypes.x86_nan(den, a).astype(np.float32)
        re = src[i0:i1, 0]
        wr = dtypes.x86_nan(np.divide(re, den32), re, den32)
        if not cplx:
            if fm is None:
                dst[i0:i1, 0] = wr
            else:
                m, w64 = fm[i0:i1], wr.astype(np.float64)
                dst[i0:i1, 0] = dtypes.x86_nan(m * w64, m, w64)
            return
        im = conj_imag(src[i0:i1, 1])
        wi = dtypes.x86_nan(np.divide(im, den32), im, den32)
        if fm is None:
            dst[i0:i1, 0], dst[i0:i1, 1] = wr, wi
        else:
            dst[i0:i1, 0], dst[i0:i1, 1] = times_complex(fm[i0:i1], wr, wi)

    with np.errstate(all="ignore"):
        _par(run, fa.size, nt)
    return out


def conj_imag(b: np.ndarray) -> np.ndarray:
    """the imaginary parts of MATLAB's ``conj`` of single values: ``0 - b`` (+0 for both zeros,
    measured)"""
    return np.subtract(np.float32(0), b)


def times_complex(m: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """MATLAB's ``m .* x`` of double reals m and complex singles x = a + bi (measured on every
    combination of zeros, Inf and NaN): the complex product ``(m + 0i)(a + bi)`` in double,
    rounded to single, with the products by the zero imaginary part signed zeros (also for an
    infinite part: ``m a - copysign(0, b)``, ``m b + copysign(0, a)``); ``b == 0`` gives
    ``(m a, +0)``; with ``a == 0`` the real part is +0 where m or b is not finite; otherwise a
    NaN part makes the other part that NaN. Not matched: both parts non-finite (one infinite) and
    m = +-0, where MATLAB gives (NaN, 0) or (0, NaN) by the signs (only for a non-finite OTF)"""
    from petakit5d.core import dtypes

    a64, b64 = a.astype(np.float64), b.astype(np.float64)
    with np.errstate(all="ignore"):
        ma = dtypes.x86_nan(m * a64, m, a64)
        mb = dtypes.x86_nan(m * b64, m, b64)
        re = ma - np.copysign(0.0, b64)
        im = mb + np.copysign(0.0, a64)
    if not (np.isfinite(re).all() and np.isfinite(im).all() and np.all(b64 != 0)):
        bz = b64 == 0
        az = (a64 == 0) & ~bz
        # a NaN part makes the other part that NaN (a and b nonzero)
        both = ~az & ~bz
        im_nan, re_nan = np.isnan(im), np.isnan(re)
        re, im = np.where(both & im_nan, im, re), np.where(both & re_nan & ~im_nan, re, im)
        re = np.where(az & ~(np.isfinite(m) & np.isfinite(b64)), 0.0, re)
        re = np.where(bz, ma, re)
        im = np.where(bz, 0.0, im)
    return re.astype(np.float32), im.astype(np.float32)


def _back_projector(otf_bp_omw: np.ndarray, backend: str | None, nt: int) -> np.ndarray:
    """``fftshift(real(ifftn(OTF_bp_omw)))`` (single)"""
    from petakit5d.matlabcompat import fft as F

    b = F.ifftn(otf_bp_omw, backend=backend, threads=nt)
    return _shift(b.real if np.iscomplexobj(b) else b, False, nt)


def _shift(a: np.ndarray, inverse: bool, nt: int) -> np.ndarray:
    """``fftshift(a)`` (``ifftshift`` with inverse) of a volume or an image, C-contiguous, plane
    by plane on threads"""
    a3 = a.reshape((1,) + a.shape) if a.ndim == 2 else a
    s = [(-(n // 2) if inverse else n // 2) % n for n in a3.shape]
    out = np.empty(a3.shape, dtype=a.dtype)
    nz, ny, nx = a3.shape
    sy, sx = s[1], s[2]

    def run(i0: int, i1: int) -> None:
        for z in range(i0, i1):
            src, dst = a3[(z - s[0]) % nz], out[z]
            dst[sy:, sx:] = src[: ny - sy, : nx - sx]
            dst[sy:, :sx] = src[: ny - sy, nx - sx :]
            dst[:sy, sx:] = src[ny - sy :, : nx - sx]
            dst[:sy, :sx] = src[ny - sy :, nx - sx :]

    _par(run, nz, nt, chunk=max(1, (1 << 16) // max(1, ny * nx)))
    return out.reshape(a.shape)


# ---------------------------------------------------------------------------------------------
# omw_backprojector_generation
# ---------------------------------------------------------------------------------------------


def _as_psf(psf: Any) -> np.ndarray:
    a = np.asarray(psf)
    if a.ndim == 3 and a.shape[0] == 1:
        a = a.reshape(a.shape[1:])
    if a.ndim not in (2, 3):
        raise ValueError("the PSF must be a (z, y, x) volume or a (y, x) image")
    return a


def generate_backprojector(
    psf,
    alpha,
    skewed: bool | None = None,
    otf_cum_thresh=0.9,
    otf_area_thresh=100,
    hann_win_bounds: Sequence[float] = (0.8, 1.0),
    min_int_thrsh=1e-3,
    *,
    backend: str | None = None,
    threads: int | None = None,
    verbose: bool = True,
) -> BackProjector:
    """omw_backprojector_generation(psf, alpha, skewed, 'OTFCumThresh', otf_cum_thresh,
    'OTFAreaThresh', otf_area_thresh, 'hannWinBounds', hann_win_bounds, 'minIntThrsh',
    min_int_thrsh) on the FFT ``backend`` (None: the default one) with ``threads`` threads (None:
    the CPUs available), printing the threshold as MATLAB does unless ``verbose`` is False (see
    the module)."""
    from petakit5d.decon.otf import psf2otf
    from petakit5d.kernels import decon as kd

    nt = _threads(threads)
    shape_in = np.shape(psf)
    a = _as_psf(psf)
    if a.size == 0 or not np.any(a != 0):
        raise BackProjectorError("", _EMPTY_PSF)
    bounds = np.asarray(hann_win_bounds, dtype=np.float64).reshape(-1)
    if bounds.size != 2:
        raise BackProjectorError(
            "MATLAB:InputParser:ArgumentFailedValidation",
            "The value of 'hannWinBounds' is invalid.",
        )

    p = _single(a, nt)
    p = _divide_single(p, kd.sum_double(p, threads=nt), nt)
    otf = psf2otf(p, backend=backend, threads=nt)
    abs_otf = _abs(otf, nt)
    abs_c = _shift(abs_otf, False, nt)
    dc = abs_otf.reshape(-1)[0]
    q = _otf_threshold(abs_otf, dc, otf_cum_thresh, min_int_thrsh, nt, verbose)

    mask = _otf_mask(abs_c, np.float32(dc * q), otf_area_thresh, nt)
    if not mask.any():
        raise BackProjectorError("", _EMPTY_MASK)
    if skewed is None:
        skewed = _is_skewed(mask, nt)
    if skewed:
        mask = _three_lobes(mask, nt)
    else:
        mask = _largest_components(mask, 1, nt)
    mask = _hulls(mask, nt)

    win = _window(mask, bounds, nt)
    if skewed:
        win = _fold(win, nt)
        mask = _fold(mask, nt)
        if a.ndim == 2:
            win, mask = win.reshape(a.shape), mask.reshape(a.shape)
    om = _shift(win, True, nt)
    otf_bp_omw = _masked_wiener(otf, abs_otf, om, float(alpha), nt)
    b = _back_projector(otf_bp_omw, backend, nt)
    out_shape = shape_in if len(shape_in) == 3 else a.shape
    return BackProjector(
        b.reshape(out_shape),
        otf_bp_omw.reshape(out_shape),
        abs_c.reshape(out_shape),
        mask.reshape(out_shape),
    )


def _bp_to_matlab(r: BackProjector):
    from petakit5d.core.arrays import to_matlab

    return tuple(to_matlab(v) for v in r)


@register(
    "omw_backprojector_generation",
    Param("psf", None, "PSF, a (z, y, x) volume or a (y, x) image.", convert="volume"),
    Param("alpha", None, "Wiener parameter."),
    Param("skewed", None, "PSF in the skewed space; None (MATLAB's []): detected."),
    Param("OTFCumThresh", "@isnumeric", "Accumulated OTF energy quantile for the OTF mask."),
    Param("OTFAreaThresh", "@isnumeric", "Smallest component of the OTF mask, in voxels."),
    Param(
        "hannWinBounds",
        "@(x) isnumeric(x) && numel(x) == 2",
        "Apodization range of the relative distance to the mask's edge.",
    ),
    Param("minIntThrsh", "@(x) isnumeric(x)", "Lower bound of the OTF threshold over DC."),
    task="R02",
    matlab_out=_bp_to_matlab,
)
def omw_backprojector(
    psf,
    alpha,
    skewed,
    *,
    otf_cum_thresh=0.9,
    otf_area_thresh=100,
    hann_win_bounds=(0.8, 1),
    min_int_thrsh=0.001,
) -> BackProjector:
    """The OTF masked Wiener back projector of a PSF (omw_backprojector_generation): the back
    projector, the masked Wiener filter, the shifted |OTF| and the OTF mask (see the module)."""
    return generate_backprojector(
        psf, alpha, skewed, otf_cum_thresh, otf_area_thresh, hann_win_bounds, min_int_thrsh
    )


# ---------------------------------------------------------------------------------------------
# omw_backprojector_files
# ---------------------------------------------------------------------------------------------


def backprojector_file_names(
    psfgen_path: str, psf_fsn: str, wiener_alpha: float, otf_cum_thresh: float
) -> dict[str, str]:
    """the psfgen files of omw_backprojector_files, as MATLAB names them (``sprintf`` with
    ``%0.6f``; the folder as given, a trailing ``/`` kept): 'bp', 'figure', 'otf', 'mask'"""
    from petakit5d.core import text

    stem = f"{psfgen_path}/{psf_fsn}"
    bp = text.sprintf(
        "%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f", stem, wiener_alpha, otf_cum_thresh
    )
    return {
        "bp": bp + ".tif",
        "figure": bp + "_figure.png",
        "otf": f"{stem}_OTF.tif",
        "mask": text.sprintf("%s_OTF_mask_otf_cum_thresh_%0.6f.tif", stem, otf_cum_thresh),
    }


def _tmp_name(final: str, uuid: str) -> str:
    """MATLAB's temporary name: the uuid before the extension"""
    root, ext = os.path.splitext(final)
    return f"{root}_{uuid}{ext}"


def _figure(abs_otf_c: np.ndarray, otf_mask: np.ndarray, path: str) -> bool:
    """visualize_OTF_and_mask_outline: log10 of the |OTF| (normalized, -3 to 0, 'hot') of the
    central xy, yz and xz planes with the outline of the mask, as a PNG; False without
    matplotlib"""
    try:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    from petakit5d.matlabcompat.morphology import imerode, strel

    a = abs_otf_c if abs_otf_c.ndim == 3 else abs_otf_c.reshape((1,) + abs_otf_c.shape)
    m = np.asarray(otf_mask) != 0
    m = m if m.ndim == 3 else m.reshape((1,) + m.shape)
    nz, ny, nx = a.shape
    cz, cy, cx = (math.ceil((n + 1) / 2) - 1 for n in (nz, ny, nx))
    # MATLAB's planes: xy (y, x), yz (y, z) and xz (x, z), each displayed transposed
    planes = [(a[cz], m[cz]), (a[:, :, cx].T, m[:, :, cx].T), (a[:, cy, :].T, m[:, cy, :].T)]
    if nz == 1:
        planes = planes[:1]
    disk = strel("disk", 1)
    labels = [("k_x", "k_y"), ("k_x", "k_z"), ("k_y", "k_z")]
    fig, axes = plt.subplots(1, len(planes), figsize=(19.2, 10.8) if nz > 1 else (8, 8))
    axes = np.atleast_1d(axes)
    with np.errstate(divide="ignore", invalid="ignore"):
        for ax, (img, msk), (xl, yl) in zip(axes, planes, labels, strict=False):
            v = np.log10(img / img.max()).T
            im = ax.imshow(v, cmap="hot", vmin=-3, vmax=0, origin="lower", interpolation="nearest")
            outline = msk & ~imerode(np.ascontiguousarray(msk), disk)
            r, c = np.nonzero(outline)
            ax.scatter(r, c, s=1)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
    fig.colorbar(im, ax=axes.tolist(), label="log$_{10}$(OTF Intensity)")
    fig.savefig(path, dpi=100, facecolor="w")
    plt.close(fig)
    return True


def omw_backprojector_files(
    psf,
    psfgen_path: str,
    psf_fsn: str,
    wiener_alpha: float,
    otf_cum_thresh: float,
    skewed: bool | None,
    hann_win_bounds: Sequence[float],
    psf_fn: str,
    save_figure: bool = True,
    *,
    backend: str | None = None,
    threads: int | None = None,
) -> np.ndarray:
    """omw_backprojector_files: the back projector of ``psf`` from the psfgen folder, or made
    (``generate_backprojector`` with the default OTFAreaThresh and minIntThrsh) and saved there
    with the OTF, the OTF mask and, with ``save_figure`` and matplotlib, the figure (see the
    module). ``psf_fn`` names the PSF in the messages. Returns the back projector (a cached
    read-only array when it was read)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.decon.psf import read_tiff_cached
    from petakit5d.io import write_tiff
    from petakit5d.jobs.publish import remove_path, rename_no_replace

    names = backprojector_file_names(psfgen_path, psf_fsn, wiener_alpha, otf_cum_thresh)
    if os.path.exists(names["bp"]):
        print(f"Load existing OMW back projector {names['bp']} for {psf_fn} ...", flush=True)
        return read_tiff_cached(names["bp"])
    print(f"OMW back projector generation for {psf_fn} ...", flush=True)
    uuid = get_uuid()
    r = generate_backprojector(
        psf,
        wiener_alpha,
        skewed,
        otf_cum_thresh,
        hann_win_bounds=hann_win_bounds,
        backend=backend,
        threads=threads,
    )
    if save_figure:
        _figure(r.abs_otf_c, r.otf_mask, names["figure"])
    nt = _threads(threads)
    for key, data in (("otf", r.abs_otf_c), ("mask", _uint8(r.otf_mask))):
        tmp = _tmp_name(names[key], uuid)
        write_tiff(data, tmp, num_threads=nt)
        os.replace(tmp, names[key])
    tmp = _tmp_name(names["bp"], uuid)
    write_tiff(r.b_omw, tmp, num_threads=nt)
    if not rename_no_replace(tmp, names["bp"]):
        # another writer published the same back projector: keep it
        remove_path(tmp)
    return r.b_omw


def _uint8(m: np.ndarray) -> np.ndarray:
    """``uint8(OTF_mask)``: logical as 0 and 1, the folded counts as they are"""
    if m.dtype == np.bool_:
        return m.view(np.uint8)
    from petakit5d.core import dtypes

    return np.asarray(dtypes.uint8(m))
