"""PSF analysis (MATLAB's ``XR_psf_analysis_wrapper``, ``XR_psf_analysis_plot`` and
``Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update``) and PSF detection and cropping
(``XR_psf_detection_and_cropping``, ``XR_psf_detection_and_analysis_wrapper``).

**The analysis of a PSF file** (``psf_otf_analysis``, MATLAB's
``Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update``): the PSF in double, its background
(``bg_factor`` times the median of its positive values) subtracted and clipped at 0, padded in
y and x to a square of at least ``psf_subpix``'s y and x (and in z to its z), its maximum (the
first in MATLAB's order) moved to the center ``(n + 1) / 2`` (rounded). Then, as MATLAB
computes them on its simulation grid (1001^3 voxels of 0.1 excitation wavelengths in the
medium, zero outside the data):

- the PSF resized (``imresize3``, linear) to the grid's pixel size, normalized by its maximum
  and centered on the grid, then moved so that the maximum of its Gaussian blur (sigma 2) is at
  the grid's center; its xy slice through the maximum, its xz slice through the maximum and the
  yz slice at the brightest of 11 x planes around the brightest column sum of the central xy
  plane, each 101 x 101 around the grid's center;
- the OTF: ``abs(fftshift(fftn(...)))`` of the ``psf_subpix`` center of the PSF, normalized,
  interpolated (``interpn``, linear) to the grid's frequency spacing, normalized and centered;
  its central xy, xz and yz planes (1001 x 1001) and the line cuts of the central 401 x 401
  parts of the xz and yz planes: along kx (ky) and kz through the center and along kz through
  the "fattest" column ``round(201 + 401 NA / 1.33 / 4)``, normalized by the lateral line's
  maximum.

The 2D outputs keep MATLAB's orientation (rows, columns): xy is ``(y, x)``, xz ``(z, x)``, yz
``(z, y)``; the line cuts are 1-D. As in MATLAB the 1001^3 grids are never made: the data
regions are kept with their offsets, and only the part of the resized PSF that lands on the
grid is resized (the whole PSF when the part's maximum may not be the whole one's).

**Figures** (matplotlib, the optional ``figures`` dependency, ``petakit5d.tools.figures``'s
conventions): the analysis figure of MATLAB's 2 x 5 panels (the PSF slices, gamma-adjusted
slices, the OTF planes and the line cuts; gids ``pk5d:psf_xy``, ``pk5d:psf_xz``,
``pk5d:psf_yz``, ``pk5d:psf_xy_gamma`` ..., ``pk5d:otf_xy``, ``pk5d:otf_xz``,
``pk5d:otf_yz``, ``pk5d:linecuts_xz``, ``pk5d:linecuts_yz``) and the line cuts against the
Richards & Wolf PSF's (``pk5d:linecuts_rw``). Pixels and fonts are matplotlib's.

**Files of a PSF** (``psf_analysis_plot``): ``<name>_infos.mat`` (v7.3, the 11 outputs),
``comp_<name>.png`` and ``.pdf`` (the analysis figure) and the figure of the line cuts,
``wT_<name>.pdf`` and then ``wT_<name>.png``, the task's result, written under a temporary
name and renamed (it needs matplotlib). The dataset function deskews the PSFs (the pipeline,
``petakit5d.pipeline``), analyzes the Richards & Wolf PSFs of the channels (``RW_info.mat`` in
each result folder ``<data>/DS/PSFAnalysis/``) and runs one task per PSF through the job
framework.

**Detection** (``psf_detection_and_cropping``): for each y the Radon transform of the
``(x, z)`` plane at the skewed angle ``atand(tand(a) * xy / (dz * sind(a)))``, its regional
maxima (8-connected) that rise above their 11 x 11 neighborhood (peak minus median above
2 (0.95 local std + 0.05 global std)), each peak's (x, z) the brightest voxel within 10 pixels
of its line; isolated peaks (no other within ``dist_thresh``) away from the border, cropped to
``crop_size`` (``<prefix>_crop_<y>_<x>_<z>.tif`` and its z MIP in ``MIPs/``, saturated ones
(above 60000) left out) and ``<prefix>.mat`` (``pcoords``, 1-based ``[y x z]`` rows, and
``pr``), saved under a temporary name and renamed: the task's result.

Kept from MATLAB: the analysis wrapper's ``objective_scan`` only sets the z pixel size (the
deskew does not get it); the detection wrapper runs its tasks with the framework's default
``parse_cluster`` (True) and the analysis without its ``xy_pixel_size``; a channel other than
488, 560 and 642 nm keeps the wavelengths of the channel before it (an error if it is the
first); a file that matches a channel pattern only as a regular expression gets no task call
(the framework retries it and gives up).

Not ported: a PSF without a positive value after the background is removed (MATLAB divides by
a zero maximum and plots NaN), and OTF regions that would not fit on the grid
(``PsfAnalysisError``).
"""

from __future__ import annotations

import math
import os
import re
from collections.abc import Sequence
from typing import Any, NamedTuple

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "PsfAnalysis",
    "PsfAnalysisError",
    "psf_analysis_dataset",
    "psf_analysis_plot",
    "psf_detection_and_analysis_dataset",
    "psf_detection_and_cropping",
    "psf_otf_analysis",
]

# MATLAB's tick labels of the line cut plots, [-1:0.2:1]
_LINECUT_LABELS = ["-1", "-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1"]
# the simulation grid of MATLAB's plots
_SIM_PIXSIZE = 0.1  # pixel size in excitation wavelengths in the medium
_SIM_HALF_FOV = 50  # half width of the field of view
_MAX_CALC_K = 0.25 / 0.1  # OTF plotted to +/- MaxCalcK
_GRID = int(2 * _SIM_HALF_FOV / _SIM_PIXSIZE + 1)  # 1001
_GAUSS_SIGMA = 2
_WAVELENGTHS = {488: (488, 515), 560: (560, 605), 642: (642, 680)}
_SATURATION = 60000

_OUTPUTS = (
    "xy_exp_PSF",
    "xz_exp_PSF",
    "yz_exp_PSF",
    "xy_exp_OTF",
    "xz_exp_OTF",
    "yz_exp_OTF",
    "xOTF_linecut",
    "yOTF_linecut",
    "zOTF_linecut",
    "zOTF_bowtie_linecut",
    "zOTF_bowtie_linecut_yz",
)


class PsfAnalysisError(ValueError):
    """an error MATLAB raises (``identifier``: its id), or a case the port does not take"""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


class PsfAnalysis(NamedTuple):
    """The outputs of ``psf_otf_analysis`` in MATLAB's order (2D slices in MATLAB's rows and
    columns, line cuts 1-D), and the figure (None without matplotlib or when not drawn)."""

    xy_exp_psf: np.ndarray
    xz_exp_psf: np.ndarray
    yz_exp_psf: np.ndarray
    xy_exp_otf: np.ndarray
    xz_exp_otf: np.ndarray
    yz_exp_otf: np.ndarray
    x_otf_linecut: np.ndarray
    y_otf_linecut: np.ndarray
    z_otf_linecut: np.ndarray
    z_otf_bowtie_linecut: np.ndarray
    z_otf_bowtie_linecut_yz: np.ndarray
    figure: Any = None

    def outputs(self) -> tuple[np.ndarray, ...]:
        """the 11 numeric outputs"""
        return tuple(self)[:11]


# ---------------------------------------------------------------------------------------------
# MATLAB's arithmetic used here
# ---------------------------------------------------------------------------------------------


def _round(v: float) -> int:
    """MATLAB's round of a double, as an int"""
    from petakit5d.core.dtypes import mround

    return int(mround(np.float64(v)))


def _colon(a: float, d: float, b: float) -> np.ndarray:
    """MATLAB's ``a:d:b`` for ``d > 0`` (measured): ``n`` steps, the count rounded and checked
    against ``b`` with a tolerance of ``2 eps max(|a|, |b|)``, the end ``c = a + n d`` (``b`` when
    within the tolerance), the first half ``a + k d`` and the second ``c - (n - k) d``, the
    middle of an even ``n`` ``(a + c) / 2``"""
    a, d, b = float(a), float(d), float(b)
    if not (d > 0) or b < a:
        return np.zeros(0)
    tol = 2.0 * 2.0**-52 * max(abs(a), abs(b))
    n = int(round((b - a) / d))
    if a + n * d - b > tol:
        n -= 1
    c = a + n * d
    if c - b > -tol:
        c = b
    k = np.arange(n // 2 + 1, dtype=np.float64)
    v = np.empty(n + 1)
    v[: k.size] = a + k * d
    v[n - k.astype(np.int64)] = c - k * d
    if n % 2 == 0:
        v[n // 2] = (a + c) / 2
    return v


def _interp_axis(v: np.ndarray, axis: int, q: np.ndarray, nt: int | None) -> np.ndarray:
    """``interpn``'s linear step along one axis (0, 1, 2) of the 1-based grid ``1..n`` at the
    points ``q``: ``(1 - t) .* v(i) + t .* v(i + 1)`` with ``i = floor(q)`` (``n - 1`` at the
    end)"""
    n = v.shape[axis]
    i = np.clip(np.floor(q).astype(np.int64), 1, max(n - 1, 1))
    t = q - i
    shape = [1] * 3
    shape[axis] = -1
    t = t.reshape(shape)
    u = 1 - t
    i0, i1 = i - 1, np.minimum(i, n - 1)
    oshape = list(v.shape)
    oshape[axis] = q.size
    out = np.empty(oshape)
    c = 1 if axis == 0 else 0  # the axis the threads split

    def run(a: int, b: int) -> None:
        sl = [slice(None)] * 3
        sl[c] = slice(a, b)
        vc = v[tuple(sl)]
        lo = np.take(vc, i0, axis=axis)
        hi = np.take(vc, i1, axis=axis)
        np.multiply(lo, u, out=lo)
        np.multiply(hi, t, out=hi)
        np.add(lo, hi, out=out[tuple(sl)])

    _par(run, v.shape[c], nt)
    return out


def _interpn(v: np.ndarray, qy: np.ndarray, qx: np.ndarray, qz: np.ndarray, nt: int | None = None
             ) -> np.ndarray:  # fmt: skip
    """MATLAB's ``interpn(x, y, z, V, xq, yq, zq)`` (linear) of a ``(z, y, x)`` volume on its
    ``ndgrid`` of ``1..n``, queried on the ``ndgrid`` of the vectors (MATLAB's dimensions y, x,
    z; measured: the steps along y, then x, then z)"""
    out = _interp_axis(v, 1, qy, nt)
    out = _interp_axis(out, 2, qx, nt)
    return _interp_axis(out, 0, qz, nt)


def _tand(x: float) -> float:
    """``tand`` of a double scalar: reduced by the nearest multiple of 90 degrees, ``tan`` of the
    rest in radians, ``-1 / t`` for an odd multiple (MATLAB's within 1 ulp; 81% of 20,016
    probed angles bit for bit)"""
    from petakit5d.matlabcompat.elementary import tan

    n = _round(x / 90)
    r = x - n * 90
    t = float(np.asarray(tan(np.float64(r * (math.pi / 180)))))
    if n % 2:
        return -1 / t if t != 0 else math.inf
    return t


def _atand(x: float) -> float:
    """``atand`` of a double scalar: ``180 / pi`` times fdlibm's ``atan`` (MATLAB's within 1 ulp;
    81% of 20,009 probed values bit for bit)"""
    from petakit5d.matlabcompat.elementary import atan

    return (180 / math.pi) * float(np.asarray(atan(np.float64(x))))


def _sind(x: float) -> float:
    from petakit5d.matlabcompat.elementary import sind

    return float(np.asarray(sind(np.float64(x))))


def _cosd(x: float) -> float:
    from petakit5d.matlabcompat.elementary import cosd

    return float(np.asarray(cosd(np.float64(x))))


def _first_max(a: np.ndarray) -> tuple[float, tuple[int, ...]]:
    """``[m, i] = max(a(:))``: the value and the index (the first in MATLAB's order)"""
    from petakit5d.matlabcompat.reduce import max as mmax

    m, idx = mmax(a)
    return float(m), tuple(int(i) for i in np.atleast_1d(np.asarray(idx)).reshape(-1))


# ---------------------------------------------------------------------------------------------
# the data on the simulation grid: a box at an offset, zero elsewhere
# ---------------------------------------------------------------------------------------------


def _par(fn, n: int, nt: int | None) -> list:
    """``fn(a, b)`` over pieces ``a:b`` of ``range(n)`` in threads (NumPy releases the
    interpreter lock in its loops); the results in order"""
    from petakit5d.core.threads import resolve_threads

    k = max(1, min(resolve_threads(nt), n))
    cuts = [n * i // k for i in range(k + 1)]
    pieces = [(cuts[i], cuts[i + 1]) for i in range(k) if cuts[i + 1] > cuts[i]]
    if len(pieces) == 1:
        return [fn(*pieces[0])]
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(len(pieces)) as ex:
        return list(ex.map(lambda ab: fn(*ab), pieces))


def _divide(v: np.ndarray, m: float, nt: int | None) -> None:
    """``v ./ m`` in place"""
    _par(lambda a, b: np.divide(v[a:b], m, out=v[a:b]), v.shape[0], nt)


def _max_places(v: np.ndarray, nt: int | None) -> tuple[float, np.ndarray]:
    """the maximum of ``v`` (no NaN) and the places (rows (z, y, x)) that hold it"""
    m = max(float(x) for x in _par(lambda a, b: v[a:b].max(), v.shape[0], nt))
    pl = _par(lambda a, b: np.argwhere(v[a:b] == m) + (a, 0, 0), v.shape[0], nt)
    return m, np.concatenate(pl, axis=0)


def _first_in_matlab_order(places: np.ndarray) -> tuple[int, int, int]:
    """the first of the places (rows (z, y, x)) in MATLAB's linear order (y fastest, then x,
    then z)"""
    i = np.lexsort((places[:, 1], places[:, 2], places[:, 0]))[0]
    return tuple(int(x) for x in places[i])


class _Grid:
    """MATLAB's 1001^3 simulation grid, zero except for the box ``v`` at ``off``: along the
    axes in ``rot`` the box spans the grid and the grid is turned by ``rot`` (``np.roll``, kept
    lazily); ``box_slice``, ``box_max`` and ``box_circshift``"""

    def __init__(self, v: np.ndarray, off: Sequence[int], nt: int | None):
        self.v = v
        self.off = [int(o) for o in off]
        self.rot = [0, 0, 0]
        self.nt = nt
        self._max: tuple[float, np.ndarray] | None = None

    def max(self) -> tuple[float, tuple[int, int, int]]:
        """``box_max``: the maximum and its first place on the grid"""
        if self._max is None:
            self._max = _max_places(self.v, self.nt)
        m, pl = self._max
        if not m > 0:
            raise PsfAnalysisError("a PSF without a positive value is not taken (MATLAB: NaN)")
        g = pl + np.asarray(self.off)
        for d in range(3):
            if self.rot[d]:
                g[:, d] = (pl[:, d] + self.rot[d]) % _GRID
        return m, _first_in_matlab_order(g)

    def slice(self, axis: int, k: int) -> np.ndarray:
        """plane ``k`` (0-based) along ``axis``: ``box_slice``"""
        if not 0 <= k < _GRID:
            dim = {0: 3, 1: 1, 2: 2}[axis]
            raise PsfAnalysisError(
                f"Slice {k + 1} along dimension {dim} is outside the {_GRID} planes."
            )
        v = self.v
        other = [d for d in range(3) if d != axis]
        out = np.zeros((_GRID, _GRID), dtype=v.dtype)
        kk = (k - self.rot[axis]) % _GRID if self.rot[axis] else k - self.off[axis]
        if 0 <= kk < v.shape[axis]:
            sl = np.take(v, kk, axis=axis)
            for j, d in enumerate(other):
                if self.rot[d]:
                    sl = np.roll(sl, self.rot[d], axis=j)
            r0, c0 = self.off[other[0]], self.off[other[1]]
            out[r0 : r0 + sl.shape[0], c0 : c0 + sl.shape[1]] = sl
        return out

    def circshift(self, shift: Sequence[int]) -> None:
        """``box_circshift``: the box moved by ``shift``, or the grid turned along the axes the
        box would wrap around (the box extended over the grid along them first)"""
        for d in range(3):
            n = self.v.shape[d]
            sh = int(shift[d])
            if sh == 0:
                continue
            if not self.rot[d] and self.off[d] + sh >= 0 and self.off[d] + sh + n <= _GRID:
                self.off[d] += sh
                continue
            if n < _GRID:
                pad = [(0, 0)] * 3
                pad[d] = (self.off[d], _GRID - self.off[d] - n)
                self.v = np.pad(self.v, pad)
                self._max = None
            self.off[d] = 0
            self.rot[d] = (self.rot[d] + sh) % _GRID


def _filtered_max(
    v: np.ndarray, off: Sequence[int], lo: Sequence[int], hi: Sequence[int], nt: int | None
) -> tuple[float, tuple[int, int, int]]:
    """``filtered_max``: the maximum of ``imgaussfilt3(grid, 2, 'FilterDomain', 'spatial')``
    over the voxels ``lo .. hi`` (inclusive) and its first place, from the region of whole y
    columns and a margin of the filter's radius in z and x (MATLAB's region: its values are the
    whole grid's)"""
    from petakit5d.matlabcompat.filters import imgaussfilt3

    r = math.ceil(2 * _GAUSS_SIGMA)
    lo = [max(int(x), 0) for x in lo]
    hi = [min(int(x), _GRID - 1) for x in hi]
    s = [max(lo[0] - r, 0), 0, max(lo[2] - r, 0)]
    t = [min(hi[0] + r, _GRID - 1), _GRID - 1, min(hi[2] + r, _GRID - 1)]
    reg = np.zeros([t[d] - s[d] + 1 for d in range(3)], dtype=v.dtype)
    bs = [max(s[d], off[d]) for d in range(3)]
    bt = [min(t[d], off[d] + v.shape[d] - 1) for d in range(3)]
    if all(bt[d] >= bs[d] for d in range(3)):
        reg[tuple(slice(bs[d] - s[d], bt[d] - s[d] + 1) for d in range(3))] = v[
            tuple(slice(bs[d] - off[d], bt[d] - off[d] + 1) for d in range(3))
        ]
    reg = imgaussfilt3(reg, _GAUSS_SIGMA, filter_domain="spatial", num_threads=nt)
    sub = reg[tuple(slice(lo[d] - s[d], hi[d] - s[d] + 1) for d in range(3))]
    m, idx = _first_max(np.ascontiguousarray(sub))
    return m, (idx[0] + lo[0], idx[1] + lo[1], idx[2] + lo[2])


def _gaussian_peak(g: _Grid) -> tuple[int, int, int]:
    """``box_gaussian_peak``: the place of the maximum of the grid's Gaussian blur (sigma 2):
    within the filter's radius of the voxels of the box that reach the blur's value at the
    box's maximum"""
    v, off, nt = g.v, g.off, g.nt
    r = math.ceil(2 * _GAUSS_SIGMA)
    _, p = g.max()
    m0, _ = _filtered_max(v, off, p, p, nt)
    thr = m0 / (1 + 1e-12)
    # the box of the voxels >= thr, from the maxima of the box's planes and lines
    zy = np.concatenate(_par(lambda a, b: v[a:b].max(axis=2), v.shape[0], nt), axis=0)
    zx = np.concatenate(_par(lambda a, b: v[a:b].max(axis=1), v.shape[0], nt), axis=0)
    zs, ys = np.nonzero(zy >= thr)
    xs = np.nonzero(zx >= thr)[1]
    lo = [int(zs.min()) + off[0] - r, int(ys.min()) + off[1] - r, int(xs.min()) + off[2] - r]
    hi = [int(zs.max()) + off[0] + r, int(ys.max()) + off[1] + r, int(xs.max()) + off[2] + r]
    m, place = _filtered_max(v, off, lo, hi, nt)
    if not m > 0:
        raise PsfAnalysisError("a PSF without a positive value is not taken (MATLAB: NaN)")
    return place


def _resize_part(p: np.ndarray, size: Sequence[int], box: Sequence[int], nt: int | None):
    """``resize_array_part``: the region ``box`` (z0, y0, x0, z1, y1, x1) of
    ``imresize3(p, size, 'linear')`` for ``p >= 0`` and the maximum of the whole resize: only
    the region is resized when its maximum exceeds the input voxels the rest of the resize
    uses (its values are their weighted means)"""
    from petakit5d.matlabcompat.resize import contributions, imresize3_fast

    size = tuple(int(s) for s in size)
    if tuple(box) != (0, 0, 0, *size):
        outer_max = -math.inf
        for d in range(3):
            _, ind = contributions(p.shape[d], size[d], size[d] / p.shape[d], "triangle", 2, True)
            rows = np.r_[0 : box[d], box[d + 3] : size[d]]
            if rows.size:
                used = np.unique(ind[rows])
                outer_max = max(outer_max, float(np.take(p, used, axis=d).max()))
        v = imresize3_fast(p, size, "linear", out_bbox=box, num_threads=nt)
        m = float(v.max())
        if m > outer_max * (1 + 1e-12):
            return v, m
    v = imresize3_fast(p, size, "linear", num_threads=nt)
    m = float(v.max())
    z0, y0, x0, z1, y1, x1 = box
    return np.ascontiguousarray(v[z0:z1, y0:y1, x0:x1]), m


def _grid_box(size: Sequence[int]) -> tuple[tuple[int, ...], list[int]]:
    """the part of a volume of ``size`` (z, y, x) on the grid (its center) and its offset there
    (MATLAB's xyi, zi: the y offset for x too)"""
    nz, ny, _ = size
    xyi = max(_round(_SIM_HALF_FOV / _SIM_PIXSIZE - ny / 2), 1)
    zi = max(_round(_SIM_HALF_FOV / _SIM_PIXSIZE - nz / 2), 1)
    lo, hi = [], []
    for n in size:
        c = _round((n + 1) / 2)
        k = min(n, _GRID)
        lh = math.ceil((k - 1) / 2)
        rh = k - lh - 1
        lo.append(c - lh - 1)
        hi.append(c + rh)
    return (*lo, *hi), [zi - 1, xyi - 1, xyi - 1]


def _check_on_grid(v: np.ndarray, off: Sequence[int]) -> None:
    if any(off[d] + v.shape[d] > _GRID for d in range(3)):
        raise PsfAnalysisError("an OTF region that does not fit on the grid is not taken")


# ---------------------------------------------------------------------------------------------
# the analysis
# ---------------------------------------------------------------------------------------------


def _prepare(
    psf: np.ndarray, psf_subpix: Sequence[int], bg_factor: float
) -> tuple[np.ndarray, tuple[int, ...]]:
    """the PSF in double with its background removed, padded and its maximum at the center,
    and the maximum's place before the move"""
    from petakit5d.matlabcompat.reduce import maximum
    from petakit5d.matlabcompat.stats import median

    a = np.asarray(psf, dtype=np.float64)
    if a.ndim == 2:
        a = a[None]
    a = np.ascontiguousarray(a)
    pos = a[a > 0]
    backgnd = float(bg_factor) * _scalar(median(pos))
    a = np.asarray(maximum(a - backgnd, 0.0), dtype=np.float64)
    nz, ny, nx = a.shape
    sz_, sy, sx = (int(v) for v in psf_subpix)
    size_xy = max(sy, sx, ny, nx)
    pad = [(0, 0), (0, 0), (0, 0)]
    for d, n in ((1, ny), (2, nx)):
        h = math.ceil((size_xy - n - 1) / 2)
        pad[d] = (h, size_xy - n - h)
    if nz < sz_:
        h = math.ceil((sz_ - nz - 1) / 2)
        pad[0] = (h, abs(sz_ - nz) - h)
    if any(p != (0, 0) for p in pad):
        a = np.pad(a, pad)
    _, peak = _first_max(a)
    shift = [_round((a.shape[d] + 1) / 2 - (peak[d] + 1)) for d in range(3)]
    return np.ascontiguousarray(np.roll(a, shift, axis=(0, 1, 2))), peak


def _crop101(s: np.ndarray) -> np.ndarray:
    h = [_round((n + 1) / 2) - 1 for n in s.shape]
    return np.ascontiguousarray(s[h[0] - 50 : h[0] + 51, h[1] - 50 : h[1] + 51])


def _otf(
    p: np.ndarray,
    psf_subpix: Sequence[int],
    xypixsize: float,
    zpixsize: float,
    index: float,
    exc_lambda: float,
    nt: int | None,
) -> tuple[np.ndarray, list[int]]:
    """the OTF of the center of the PSF on the grid: the region and its offset"""
    from petakit5d.matlabcompat.elementary import hypot
    from petakit5d.matlabcompat.fft import fftn

    sub = [int(v) for v in psf_subpix]
    start = [max(_round((p.shape[d] - sub[d]) / 2), 1) - 1 for d in range(3)]
    subset = np.ascontiguousarray(p[tuple(slice(start[d], start[d] + sub[d]) for d in range(3))])
    f = fftn(subset, threads=nt)
    otf = hypot(f.real, f.imag) if np.iscomplexobj(f) else np.abs(f)
    otf = np.fft.fftshift(np.asarray(otf, dtype=np.float64))
    otf = otf / float(otf.max())
    cz, cy, cx = otf.shape
    xy_half_fov = cy * xypixsize / 2
    z_half_fov = cz * zpixsize / 2
    sim_half_fovnm = _SIM_HALF_FOV * (exc_lambda / index)
    sim_otf = 1 / sim_half_fovnm
    magxy = (1 / xy_half_fov) / sim_otf
    magz = (1 / z_half_fov) / sim_otf
    d = [_round(magz * cz), _round(magxy * cy), _round(magxy * cx)]
    qz = _colon(1, cz / d[0], cz)
    qy = _colon(1, cy / d[1], cy)
    qx = _colon(1, cx / d[2], cx)
    o3 = _interpn(otf, qy, qx, qz, nt)
    o3 = o3 / float(o3.max())
    box, off = _grid_box(o3.shape)
    z0, y0, x0, z1, y1, x1 = box
    v = np.ascontiguousarray(o3[z0:z1, y0:y1, x0:x1])
    _check_on_grid(v, off)
    return v, off


def _linecuts(plane: np.ndarray, na_det: float, rounded: bool = True):
    """the central part of an OTF plane and its line cuts (lateral, axial, axial at the
    fattest column), normalized by the lateral line's maximum, and the fattest column
    (1-based)"""
    n = plane.shape[0]
    lo = (1 - 1 / _MAX_CALC_K) * (n - 1) / 2 + 1
    hi = (1 + 1 / _MAX_CALC_K) * (n - 1) / 2 + 1
    if rounded:
        lo, hi = _round(lo), _round(hi)
    part = np.ascontiguousarray(plane[int(lo) - 1 : int(hi), int(lo) - 1 : int(hi)])
    a = part.shape[0]
    mid = int((a + 1) / 2)
    fattest = _round(mid + a * na_det / 1.33 / 4)
    if not 1 <= fattest <= part.shape[1]:
        raise PsfAnalysisError(
            f"Index in position 2 exceeds array bounds. Index must not exceed {part.shape[1]}.",
            "MATLAB:badsubscript",
        )
    lateral = part[mid - 1, :]
    m = float(lateral.max())
    return part, lateral / m, part[:, mid - 1] / m, part[:, fattest - 1] / m, fattest


def _read_psf(filenm: str) -> np.ndarray:
    from petakit5d.io import read_tiff

    return read_tiff(filenm)


def psf_otf_analysis(
    filenm: str,
    source_descrip: str,
    xypixsize: float,
    zpixsize: float,
    na_det: float,
    index: float,
    exc_lambda: float,
    det_lambda: float,
    psf_subpix: Sequence[int] = (101, 128, 128),
    gamma: float = 0.5,
    bg_factor: float = 1.5,
    visible: bool = True,
    *,
    plot: bool = True,
    num_threads: int | None = None,
) -> PsfAnalysis:
    """MATLAB's ``Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update`` (see the module): the PSF
    and OTF slices and line cuts of the PSF tiff ``filenm`` (pixel sizes in nm, wavelengths in
    nm, ``psf_subpix`` the OTF's region ``(z, y, x)``) and, with ``plot`` (Python only) and
    matplotlib, the analysis figure (a pyplot figure when ``visible``). ``num_threads`` (Python
    only): the threads of the kernels; the results do not depend on it."""
    del det_lambda  # (MATLAB's argument, not used)
    xypixsize, zpixsize, index = float(xypixsize), float(zpixsize), float(index)
    exc_lambda, na_det = float(exc_lambda), float(na_det)
    psf_subpix = [int(v) for v in psf_subpix]
    p, raw_peak = _prepare(_read_psf(filenm), psf_subpix, bg_factor)

    # the PSF resized to the grid's pixel size, the part on the grid
    sim_xzpixsize = _SIM_PIXSIZE * (exc_lambda / index)
    magxy = xypixsize / sim_xzpixsize
    magz = zpixsize / sim_xzpixsize
    nz, ny, nx = p.shape
    b = (_round(magz * nz), _round(magxy * ny), _round(magxy * nx))
    box, off = _grid_box(b)
    nt = num_threads
    v, vmax = _resize_part(p, b, box, nt)
    if not vmax > 0:
        raise PsfAnalysisError("a PSF without a positive value is not taken (MATLAB: NaN)")
    _divide(v, vmax, nt)
    g = _Grid(v, off, nt)

    # centered on the maximum of its Gaussian blur
    c = (_GRID + 1) / 2
    peak = _gaussian_peak(g)
    if any(abs(c - (peak[d] + 1)) > 0.5 for d in range(3)):
        g.circshift([_round(c - (peak[d] + 1)) for d in range(3)])
        _, peak = g.max()

    # the brightest column sum of the central xy plane, then the brightest of 11 x planes
    from petakit5d.matlabcompat.reduce import sum as msum

    mid = (_GRID + 1) // 2 - 1
    xy = g.slice(0, mid)
    xy = xy / float(xy.max())
    _, (peakx,) = _first_max(np.asarray(msum(xy, axis=-2)).reshape(-1))
    maxval = [float(g.slice(2, peakx + k).max()) for k in range(-5, 6)]
    _, (best,) = _first_max(np.asarray(maxval))
    yz_psf = _crop101(g.slice(2, peakx + best - 5))
    _, (gz, gy, _) = g.max()
    xy_psf = _crop101(g.slice(0, gz))
    xz_psf = _crop101(g.slice(1, gy))
    del g, v

    # the OTF
    go = _Grid(*_otf(p, psf_subpix, xypixsize, zpixsize, index, exc_lambda, nt), nt)
    yz_otf = go.slice(2, mid)
    xy_otf = go.slice(0, mid)
    xz_otf = go.slice(1, mid)
    del go
    yz_part, y_cut, _, bowtie_yz, f_yz = _linecuts(yz_otf, na_det)
    xz_part, x_cut, z_cut, bowtie, f_xz = _linecuts(xz_otf, na_det)

    res = PsfAnalysis(
        xy_psf, xz_psf, yz_psf, xy_otf, xz_otf, yz_otf, x_cut, y_cut, z_cut, bowtie, bowtie_yz
    )
    if plot:
        from petakit5d.tools import figures

        if figures.available():
            xy_part = _linecuts(xy_otf, na_det, rounded=False)[0]
            fig = _analysis_figure(
                res,
                {"xy": xy_part, "xz": xz_part, "yz": yz_part},
                {"xz": f_xz, "yz": f_yz},
                peak,
                raw_peak,
                str(source_descrip),
                float(gamma),
                bool(visible),
            )
            res = res._replace(figure=fig)
    return res


# ---------------------------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------------------------


def _hot256():
    """MATLAB's ``hot(256)`` as a matplotlib colormap (``image`` of ``256 .* x`` indexes it
    directly: row ``fix(256 x)``, clamped to 1..256)"""
    from matplotlib.colors import ListedColormap

    n = 256
    k = np.arange(1, n + 1, dtype=np.float64)
    n1 = math.floor(3 / 8 * n)
    n2 = math.floor(3 / 8 * n)
    r = np.minimum(k / n1, 1)
    g = np.clip((k - n1) / n2, 0, 1)
    b = np.clip((k - n1 - n2) / (n - n1 - n2), 0, 1)
    return ListedColormap(np.stack([r, g, b], axis=1), name="hot256")


def _direct_index(x: np.ndarray) -> np.ndarray:
    """the colormap rows (0-based) of ``image(256 .* x)``"""
    with np.errstate(invalid="ignore"):
        k = np.fix(256 * np.asarray(x, dtype=np.float64))
    return np.clip(np.nan_to_num(k, nan=1), 1, 256) - 1


def _ticks(ax, n: int, labels_x, labels_y) -> None:
    t = 1 + (n - 1) / 4 * np.arange(5)
    ax.set_xticks(t)
    ax.set_xticklabels(labels_x)
    ax.set_yticks(t)
    ax.set_yticklabels(labels_y)


def _num(v: float) -> str:
    from petakit5d.core.text import num2str

    return num2str(float(v))


def _image_panel(fig, pos: int, gid: str, data: np.ndarray, xl: str, yl: str, title: str,
                 psf: bool, cmap) -> Any:  # fmt: skip
    ax = fig.add_subplot(2, 5, pos)
    ax.set_gid(gid)
    n = data.shape[0]
    ax.imshow(
        _direct_index(data),
        cmap=cmap,
        vmin=0,
        vmax=255,
        origin="upper",
        extent=(0.5, data.shape[1] + 0.5, data.shape[0] + 0.5, 0.5),
        interpolation="nearest",
    )
    ax.set_xlim(1, n)
    ax.set_ylim(data.shape[1], 1)
    ax.set_aspect("equal")
    if psf:
        lab = [_num((n - 1) / 2 * _SIM_PIXSIZE * f) for f in (-1, -0.5, 0, 0.5, 1)]
        _ticks(ax, n, lab, lab)
    else:
        lab = ["-1", "-0.5", "0", "0.5", "1"]
        _ticks(ax, n, lab, lab)
    ax.set_xlabel(xl, fontsize=14)
    ax.set_ylabel(yl, fontsize=14)
    ax.set_title(title, fontsize=12)
    return ax


def _linecut_panel(ax, cuts, styles, names, colors) -> None:
    n = cuts[0].size
    with np.errstate(divide="ignore", invalid="ignore"):
        for c, s in zip(cuts, styles, strict=True):
            ax.plot(np.arange(1, n + 1), np.log10(c), s, linewidth=2)
    ax.set_xlim(1, n)
    ax.set_ylim(-3, 0)
    ax.set_box_aspect(1)
    ax.grid(True)
    ax.set_xticks(1 + (n - 1) / 10 * np.arange(11))
    ax.set_xticklabels(_LINECUT_LABELS)
    ax.set_yticks([-3, -2, -1, 0])
    ax.set_yticklabels(["0.001", "0.01", "0.1", "1"])
    ax.set_xlabel(r"k / (4$\pi$/$\lambda$)", fontsize=14)
    ax.set_ylabel("OTF Strength", fontsize=14)
    for k, (t, col) in enumerate(zip(names, colors, strict=True)):
        ax.text(0.6 * n, -0.15 * (k + 1), t, color=col, fontsize=14)


def _analysis_figure(res: PsfAnalysis, parts, fattest, peak, raw_peak, source: str,
                     gamma: float, visible: bool):  # fmt: skip
    """MATLAB's figure of the analysis: 2 x 5 panels (2500 x 1100 pixels)"""
    from petakit5d.tools.figures import _figure

    fig = _figure(2500, 1100, visible)
    cmap = _hot256()
    g = _num_gamma(gamma)
    lam = r"$\lambda_{exc}$/n"
    k4 = r"4$\pi$n/$\lambda_{exc}$"
    gx = (_GRID + 1) / 2
    pz, py, px = peak
    panels = [
        (1, "pk5d:psf_xy", res.xy_exp_psf, f"x / ({lam})", f"y / ({lam})",
         f"Overall PSF From {source}", True, (px + 1 - (gx - 51), py + 1 - (gx - 51))),
        (2, "pk5d:psf_xz", res.xz_exp_psf, f"x / ({lam})", f"z / ({lam})",
         f"Overall PSF From {source}", True, (px + 1 - (gx - 51), px + 1 - (gx - 51))),
        (4, "pk5d:psf_yz", res.yz_exp_psf, f"y / ({lam})", f"z / ({lam})",
         f"Overall PSF From {source}", True, (py + 1 - (gx - 51), pz + 1 - (gx - 51))),
        (7, "pk5d:psf_xz_gamma", res.xz_exp_psf**gamma, f"x / ({lam})", f"z / ({lam})",
         f"Overall PSF From {source}, gamma = {g}", True, None),
        (9, "pk5d:psf_yz_gamma", res.yz_exp_psf**gamma, f"y / ({lam})", f"z / ({lam})",
         f"Overall PSF From {source}, gamma = {g}", True, None),
        (3, "pk5d:otf_xz", parts["xz"] ** gamma, f"$k_x$ / ({k4})", f"$k_z$ / ({k4})",
         f"Overall OTF From {source}, gamma = {g}", False, None),
        (5, "pk5d:otf_yz", parts["yz"] ** gamma, f"$k_y$ / ({k4})", f"$k_z$ / ({k4})",
         f"Overall OTF From {source}, gamma = {g}", False, None),
        (6, "pk5d:otf_xy", parts["xy"] ** gamma, f"$k_x$ / ({k4})", f"$k_y$ / ({k4})",
         f"Overall OTF From {source}, gamma = {g}", False, None),
    ]  # fmt: skip
    for pos, gid, data, xl, yl, title, psf, star in panels:
        ax = _image_panel(fig, pos, gid, data, xl, yl, title, psf, cmap)
        n = data.shape[0]
        if star is not None:
            ax.plot((n + 1) / 2, (n + 1) / 2, "o")
            ax.plot(star[0], star[1], "*")
        if not psf and gid != "pk5d:otf_xy":
            m = (n + 1) / 2
            ax.plot([1, n], [m, m], color=(0, 0, 1), linewidth=1)
            ax.plot([m, m], [1, n], color=(1, 0, 0), linewidth=1)
            f = fattest[gid[-2:]]
            ax.plot([f, f], [1, n], color=(0, 1, 0), linewidth=1)
    for pos, gid, cuts, lat in (
        (8, "pk5d:linecuts_xz", (res.z_otf_linecut, res.x_otf_linecut, res.z_otf_bowtie_linecut),
         "kx"),
        (10, "pk5d:linecuts_yz", (_axial_yz(parts["yz"], res), res.y_otf_linecut,
                                  res.z_otf_bowtie_linecut_yz), "ky"),
    ):  # fmt: skip
        ax = fig.add_subplot(2, 5, pos)
        ax.set_gid(gid)
        _linecut_panel(
            ax,
            cuts,
            ("r", "b", "g"),
            (f"OTF along {lat}", "OTF along kz", "Bowtie OTF along kz"),
            ((0, 0, 1), (1, 0, 0), (0, 0.75, 0)),
        )
        ax.set_title(f"Overall OTF linecuts From {source}", fontsize=12)
    fig.text(
        0.5,
        0.02,
        "Image y = x galvo dither direction & y-stage scan direction\n"
        "In deskewed data image x = LLS(y) propogation direction\n"
        f"Peak coordinate (xyz): ({raw_peak[2] + 1}, {raw_peak[1] + 1}, {raw_peak[0] + 1})",
        ha="center",
        fontsize=14,
    )
    return fig


def _axial_yz(part: np.ndarray, res: PsfAnalysis) -> np.ndarray:
    """the axial line cut of the yz OTF plane (plotted, not returned by MATLAB)"""
    mid = (part.shape[0] + 1) // 2 - 1
    lateral = part[mid, :]
    return part[:, mid] / float(lateral.max())


def _num_gamma(gamma: float) -> str:
    from petakit5d.core.text import num2str

    return num2str(float(gamma), "%1.2f")


def _rw_figure(res: PsfAnalysis, rw: Sequence[np.ndarray], source: str, visible: bool):
    """MATLAB's figure of the line cuts against the Richards & Wolf PSF's (600 x 600)"""
    from petakit5d.tools.figures import _figure

    fig = _figure(600, 600, visible)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_gid("pk5d:linecuts_rw")
    n = res.z_otf_linecut.size
    x = np.arange(1, n + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        for c, s in (
            (res.z_otf_linecut, "r"),
            (res.x_otf_linecut, "b"),
            (res.z_otf_bowtie_linecut, "g"),
        ):
            ax.plot(x, np.log10(c), s, linewidth=2)
        for c, s in ((rw[4], "--r"), (rw[2], "--b"), (rw[5], "--g")):
            c = np.asarray(c, dtype=np.float64).reshape(-1)
            ax.plot(np.arange(1, c.size + 1), np.log10(c), s, linewidth=2)
    m = np.asarray(rw[4]).size
    ax.set_xlim(1, m)
    ax.set_ylim(-3, 0)
    ax.set_box_aspect(1)
    ax.grid(True)
    ax.set_xticks(1 + (m - 1) / 10 * np.arange(11))
    ax.set_xticklabels(_LINECUT_LABELS)
    ax.set_yticks([-3, -2, -1, 0])
    ax.set_yticklabels(["0.001", "0.01", "0.1", "1"])
    ax.set_xlabel(r"k / (4$\pi$/$\lambda$)", fontsize=14)
    ax.set_ylabel("OTF Strength", fontsize=14)
    a = res.xz_exp_otf.shape[1]
    for k, (t, col) in enumerate(
        (("OTF along kx", (0, 0, 1)), ("OTF along kz", (1, 0, 0)),
         ("Bowtie OTF along kz", (0, 0.75, 0)))
    ):  # fmt: skip
        ax.text(0.6 * a, -0.15 * (k + 1), t, color=col, fontsize=14)
    ax.set_title(f"Overall OTF linecuts From {source}", fontsize=14)
    ax.legend(
        [
            "OTF along kz",
            "OTF along kx",
            "Bowtie OTF along kx",
            "RW OTF along kz",
            "RW OTF along kx",
            "RW Bowtie OTF along kz",
        ]
    )
    return fig


# ---------------------------------------------------------------------------------------------
# the task of a PSF file
# ---------------------------------------------------------------------------------------------


def _rw_entry(rw_info: Any, ch_ind: int) -> list[np.ndarray]:
    """``RW_info{ch_ind}`` of a loaded RW_info.mat (a cell column)"""
    cells = rw_info
    if isinstance(cells, np.ndarray):
        cells = list(cells.reshape(-1, order="F"))
    elif not isinstance(cells, list):
        cells = [cells]
    if not 0 <= ch_ind < len(cells):
        raise PsfAnalysisError(
            f"Index exceeds the number of array elements. Index must not exceed {len(cells)}.",
            "MATLAB:badsubscript",
        )
    e = cells[ch_ind]
    if isinstance(e, np.ndarray):
        e = list(e.reshape(-1, order="F"))
    return [np.asarray(v, dtype=np.float64) for v in e]


def _savefig(fig, path: str, fmt: str) -> None:
    fig.savefig(path, format=fmt, facecolor=fig.get_facecolor(), dpi=100)


def _close(fig) -> None:
    """MATLAB's ``close``: a pyplot figure is closed, another one dropped"""
    import sys

    plt = sys.modules.get("matplotlib.pyplot")
    if plt is not None and fig is not None and plt.fignum_exists(getattr(fig, "number", -1)):
        plt.close(fig)


@register(
    "XR_psf_analysis_plot",
    Param("frameFullname", "@ischar", "The PSF file (tiff)."),
    Param("figureFullname", "@ischar", "The figure of its line cuts (png), the task's result."),
    Param("RW_info_Fullname", "@ischar", "RW_info.mat of the Richards & Wolf PSFs.",
          py="rw_info_fullname"),
    Param("ch_ind", "@isscalar", "Index of the channel in RW_info.", convert="index1"),
    Param("source_descrip", "@ischar", "Description of the source shown in the figures."),
    Param("xypixsize", "@isscalar", "xy pixel size (nm).", convert="same"),
    Param("zpixsize", "@isscalar", "z pixel size (nm).", convert="same"),
    Param("NAdet", "@isscalar", "Numerical aperture of the detection objective.", py="na_det"),
    Param("index", "@isscalar", "Refractive index of the medium."),
    Param("exc_lambda", "@isscalar", "Excitation wavelength (nm)."),
    Param("det_lambda", "@isscalar", "Detection wavelength (nm)."),
    Param("PSFsubpix", "@isvector", "Region (z, y, x) of the PSF for the OTF.",
          convert="size_yxz", py="psf_subpix"),
    Param("gamma", "@isscalar", "Contrast factor (0 < gamma < 1) of the secondary plots."),
    Param("bgFactor", "@isscalar", "Background factor (times the median of positive values)."),
    Param("visible", "@islogical", "Show the figures."),
    task="B06",
)  # fmt: skip
def psf_analysis_plot(
    frame_fullname: str,
    figure_fullname: str,
    rw_info_fullname: str,
    ch_ind: int,
    source_descrip: str,
    xypixsize: float,
    zpixsize: float,
    na_det: float,
    index: float,
    exc_lambda: float,
    det_lambda: float,
    psf_subpix: Sequence[int],
    gamma: float,
    bg_factor: float,
    *,
    visible: bool = True,
) -> None:
    """The analysis of one PSF file and its files (see the module): ``<name>_infos.mat``,
    ``comp_<name>.png`` / ``.pdf`` and the line cuts against the Richards & Wolf PSF of channel
    ``ch_ind`` (0-based) as ``figure_fullname`` (png, written last) and its pdf. Nothing if the
    figure exists. Needs matplotlib (the ``figures`` extra)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.io import read_mat, write_mat
    from petakit5d.jobs.files import exist
    from petakit5d.tools import figures

    if exist(figure_fullname):
        print("The figure already exists, skip it!", flush=True)
        return
    if not figures.available():
        raise ImportError(
            "psf_analysis_plot writes figures: it needs matplotlib (pip install "
            "'petakit5d[figures]')"
        )
    rw = _rw_entry(read_mat(rw_info_fullname)["RW_info"], int(ch_ind))
    res = psf_otf_analysis(
        frame_fullname,
        source_descrip,
        xypixsize,
        zpixsize,
        na_det,
        index,
        exc_lambda,
        det_lambda,
        psf_subpix,
        gamma,
        bg_factor,
        visible,
    )
    fsname = fileparts(frame_fullname)[1]
    result_dir = fileparts(figure_fullname)[0]
    uuid = get_uuid()
    tmp = f"{result_dir}/{fsname}_info_{uuid}.mat"
    write_mat(tmp, dict(zip(_OUTPUTS, res.outputs(), strict=True)))
    os.replace(tmp, f"{result_dir}/{fsname}_infos.mat")

    fig = res.figure
    _savefig(fig, os.path.join(result_dir, f"comp_{fsname}.png"), "png")
    _savefig(fig, os.path.join(result_dir, f"comp_{fsname}.pdf"), "pdf")
    _close(fig)

    fig = _rw_figure(res, rw, str(source_descrip), bool(visible))
    root = figure_fullname[:-4]
    _savefig(fig, f"{root}.pdf", "pdf")
    # the PNG is the result the job framework checks: written last, under a temporary name
    tmp = f"{root}_{uuid}.png"
    _savefig(fig, tmp, "png")
    os.replace(tmp, figure_fullname)
    if not visible:
        _close(fig)


# ---------------------------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------------------------


def _paths(data_paths: str | Sequence[str]) -> list[str]:
    paths = [data_paths] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    if os.name == "nt":
        paths = [p.replace("\\", "/") for p in paths]
    return paths


def _wavelengths(channel: float, prev: tuple[float, float] | None) -> tuple[float, float]:
    """MATLAB's switch of the channel: the excitation and detection wavelengths (the previous
    channel's for another channel; an error for the first)"""
    w = _WAVELENGTHS.get(int(channel)) if float(channel) == int(channel) else None
    if w is not None:
        return float(w[0]), float(w[1])
    if prev is None:
        raise PsfAnalysisError(
            "Unrecognized function or variable 'exc_lambda'.", "MATLAB:UndefinedFunction"
        )
    return prev


def _noop() -> None:
    """the task of a file without a channel (MATLAB's empty function string)"""


def _txt(v: float) -> float:
    """a value after MATLAB's ``%.20f`` in a task's function string"""
    from petakit5d.core.text import through_text

    return float(through_text(float(v), "%.20f"))


@register(
    "XR_psf_analysis_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders with the PSFs."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Scan interval in um."),
    Param("skewAngle", "@isnumeric", "Skew angle (in degree) of the stage."),
    Param("deskew", "@islogical", "Deskew the PSFs first."),
    Param("flipZstack", "@islogical", "Flip z stacks."),
    Param("objectiveScan", "@islogical", "Objective scan (sets the z pixel size only)."),
    Param("zStageScan", "@islogical", "Z stage scan (orthogonal to objective scan)."),
    Param("channelPatterns", "@iscell", "Channel identifiers for included channels."),
    Param("channels", "@isnumeric", "Wavelength (nm) of each channel."),
    Param("save16bit", "@islogical", "Save 16bit result for the deskew.", convert="same"),
    Param("bgFactor", "@isnumeric", "PSF background factor (times the median)."),
    Param("RWFn", "@iscell", "Richards & Wolf widefield PSF of each channel.", py="rw_fn"),
    Param("sourceStr", "@ischar", "Source of the PSF shown in the titles of the plots."),
    Param("visible", "@islogical", "Show the figures."),
    Param("parseCluster", "@islogical", "Use the Slurm cluster for the processing."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("cpusPerTask", "@isscalar", "CPUs per task of the deskew jobs."),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only)."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="B06",
)
def psf_analysis_dataset(
    data_paths: str | Sequence[str],
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.1,
    skew_angle: float = 32.45,
    deskew: bool = True,
    flip_zstack: bool = False,
    objective_scan: bool = False,
    z_stage_scan: bool = False,
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0"),
    channels: Sequence[float] = (488, 560),
    save_16bit: bool = True,
    bg_factor: float = 1.5,
    rw_fn: Sequence[str] = ("", ""),
    source_str: str = "test",
    visible: bool = True,
    parse_cluster: bool = True,
    master_compute: bool = False,
    cpus_per_task: float = 8,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """PSF analysis of the isolated PSFs of datasets (``XR_psf_analysis_wrapper``; see the
    module): the PSFs deskewed into ``DS/`` (with ``deskew``), the Richards & Wolf PSFs of the
    channels analyzed into ``RW_info.mat``, then one ``psf_analysis_plot`` task per PSF that
    matches a channel pattern, into ``<data>/DS/PSFAnalysis/``."""
    from petakit5d.core.paths import matlab_dir
    from petakit5d.core.text import through_mat2str
    from petakit5d.io import get_image_size, write_mat
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import contains_text, matches_channel_pattern

    paths = _paths(data_paths)
    exps = [f"{p}/" for p in paths]
    print(exps, flush=True)
    channel_patterns = [str(c) for c in channel_patterns]
    channels = [float(c) for c in np.atleast_1d(np.asarray(channels, dtype=np.float64))]
    rw_fn = [str(f) for f in rw_fn]

    if deskew:
        from petakit5d.pipeline.automatic import automatic_processing

        automatic_processing(
            exps,
            xy_pixel_size=xy_pixel_size,
            dz=dz,
            skew_angle=skew_angle,
            reverse=True,
            z_stage_scan=z_stage_scan,
            channel_patterns=tuple(channel_patterns),
            save_16bit=save_16bit,
            overwrite=False,
            streaming=False,
            cpus_per_task=cpus_per_task,
            parse_cluster=parse_cluster,
            master_compute=master_compute,
            mcc_mode=mcc_mode,
            config_file=config_file,
            deskew=True,
            rotate=False,
            dsr_combined=False,
            parse_setting_file=False,
            flip_zstack=flip_zstack,
            ff_correction=False,
        )
        exps = [f"{p}/DS/" for p in paths]
    print(exps, flush=True)

    na_det, index, gamma = 1.0, 1.33, 0.5
    xypixsize = xy_pixel_size * 1000
    if objective_scan:
        zpixsize = dz * 1000
    elif z_stage_scan:
        zpixsize = dz * _cosd(skew_angle) * 1000
    else:
        zpixsize = dz * _sind(skew_angle) * 1000
    psf_subpix = (_round(100 * 100 / zpixsize) + 1, 128, 128)

    # the Richards & Wolf PSFs of the channels
    rw_info = np.empty((len(channel_patterns), 1), dtype=object)
    wl = None
    for c in range(len(channel_patterns)):
        wl = _wavelengths(channels[c], wl)
        r = psf_otf_analysis(
            rw_fn[c], source_str, xypixsize, 0.1 * 1000, na_det, index, wl[0], wl[1],
            (101, 128, 128), gamma, 0, visible, plot=False,
        )  # fmt: skip
        rw_info[c, 0] = [r.xz_exp_psf, r.xz_exp_otf, r.x_otf_linecut, r.y_otf_linecut,
                         r.z_otf_linecut, r.z_otf_bowtie_linecut]  # fmt: skip

    rw_files = []
    for rtd in exps:
        result_dir = f"{rtd}PSFAnalysis{os.sep}"
        os.makedirs(result_dir, exist_ok=True)
        rw_files.append(f"{result_dir}/RW_info.mat")
        write_mat(rw_files[-1], {"RW_info": rw_info})

    sub = through_mat2str(np.asarray([psf_subpix[1], psf_subpix[2], psf_subpix[0]], np.float64))
    sub_py = [int(sub[2]), int(sub[0]), int(sub[1])]

    frames, figs, tasks = [], [], []
    for d, rtd in enumerate(exps):
        result_dir = f"{rtd}PSFAnalysis{os.sep}"
        fsn = [e.name for e in matlab_dir(f"{rtd}*.tif")]
        fn = [rtd + x for x in fsn]
        keep = [False] * len(fn)
        for p in channel_patterns:
            keep = [a or b for a, b in zip(keep, matches_channel_pattern(fn, p), strict=True)]
        fn = [x for x, k in zip(fn, keep, strict=True) if k]
        fsn = [x for x, k in zip(fsn, keep, strict=True) if k]
        for f, name in zip(fn, fsn, strict=True):
            frames.append(f)
            figs.append(f"{result_dir}wT_{name[:-4]}.png")
            hits = [i for i, p in enumerate(channel_patterns) if contains_text([f], p)[0]]
            if not hits:
                tasks.append(Task.of(_noop))
                continue
            if len(hits) > 1:
                raise PsfAnalysisError(
                    "SWITCH expression must be a scalar or a character vector.",
                    "MATLAB:switch:caseExpressionWrongType",
                )
            exc, det = _wavelengths(channels[hits[0]], wl)
            tasks.append(
                Task.of(
                    psf_analysis_plot,
                    f,
                    figs[-1],
                    rw_files[d],
                    hits[0],
                    str(source_str),
                    _txt(xypixsize),
                    _txt(zpixsize),
                    _txt(na_det),
                    _txt(index),
                    _txt(exc),
                    _txt(det),
                    sub_py,
                    _txt(gamma),
                    _txt(bg_factor),
                    visible=bool(visible),
                )
            )
    if not frames:
        raise PsfAnalysisError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    mem = float(np.prod(get_image_size(frames[0]), dtype=np.float64)) * 4 / 1024**3 * 100
    done = [False]
    for i in range(1, 4):
        if all(done):
            break
        done = framework.run_tasks(
            frames,
            figs,
            tasks,
            max_trial_num=2,
            parse_cluster=parse_cluster,
            master_compute=master_compute,
            cpus_per_task=4,
            mem_allocate=mem * i,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )


# ---------------------------------------------------------------------------------------------
# PSF detection
# ---------------------------------------------------------------------------------------------


def _radon_geometry(m: int, n: int) -> tuple[int, int, int]:
    """the center (1-based row and column) and the last radial bin of radon's grid"""
    yc = (m + 1) // 2
    xc = (n + 1) // 2
    rlast = math.ceil(math.sqrt((m - (m - 1) // 2 - 1) ** 2 + (n - (n - 1) // 2 - 1) ** 2)) + 1
    return yc, xc, rlast


def _radon_xp(m: int, n: int) -> np.ndarray:
    """radon's radial coordinates ``xp`` of an ``m x n`` image"""
    _, _, rlast = _radon_geometry(m, n)
    return np.arange(-rlast, rlast + 1, dtype=np.float64)


# subpixel offsets (x, y) in the order they are added
_SUBPIXELS = ((-0.25, 0.25), (-0.25, -0.25), (0.25, -0.25), (0.25, 0.25))


def _radon_weights(m: int, n: int, theta: float) -> tuple[np.ndarray, np.ndarray, int]:
    """the radial bins of radon's grid that each pixel of an ``m x n`` image reaches at
    ``theta`` and its weights there: the first bin (n, m) and the weights (n, m, k), the
    subpixels' shares added in their order; the number of bins ``2 rlast + 1``"""
    from petakit5d.matlabcompat.elementary import fdlibm_cos, fdlibm_sin

    yc, xc, rlast = _radon_geometry(m, n)
    t = np.array([theta * math.pi / 180])
    c = float(fdlibm_cos(t)[0])
    s = float(fdlibm_sin(t)[0])
    xs = (np.arange(1, n + 1, dtype=np.float64) - xc) * c
    ys = (yc - np.arange(1, m + 1, dtype=np.float64)) * s
    # r - rfirst of each subpixel (column, row)
    rs = [xs[:, None] + a * c + ys[None, :] + b * s + rlast for a, b in _SUBPIXELS]
    r1 = [r.astype(np.int64) for r in rs]
    b0 = np.minimum.reduce(r1)
    k = int(max(int((q - b0).max()) for q in r1)) + 2
    w = np.zeros((n, m, k))
    for r, q in zip(rs, r1, strict=True):
        d = r - q
        j = (q - b0)[..., None] == np.arange(k)
        w += np.where(j, (1 - d)[..., None], 0)
        w += np.where(np.roll(j, 1, axis=-1), d[..., None], 0)
    return b0, w, 2 * rlast + 1


def _radon(planes: np.ndarray, theta: float, nt: int | None = None) -> np.ndarray:
    """MATLAB's ``radon(I, theta)`` of each ``m x n`` image of ``planes`` (k, m, n) for one
    angle (degrees), from its documented method: each pixel is split into four subpixels whose
    values (a quarter of the pixel's) are projected onto the radial axis, ``x cos + y sin`` from
    the center pixel ``floor((size + 1) / 2)`` (y up), and split linearly between the two
    nearest bins; ``(k, 2 rlast + 1)``. The sums run in a fixed order (column-major pixels);
    their last bits are not MATLAB's (``docs/parity/tools.md``)."""
    k, m, n = planes.shape
    b0, w, nb = _radon_weights(m, n, theta)
    kw = w.shape[-1]
    bins = (b0[..., None] + np.arange(kw)).reshape(-1)
    w = w.reshape(n * m, kw)
    out = np.empty((k, nb))
    per = max(1, (1 << 21) // bins.size)

    def run(a: int, b: int) -> None:
        for p0 in range(a, b, per):
            p1 = min(b, p0 + per)
            pix = np.transpose(planes[p0:p1], (0, 2, 1)).reshape(p1 - p0, n * m, 1) * 0.25
            idx = bins[None, :] + (np.arange(p1 - p0) * nb)[:, None]
            out[p0:p1] = np.bincount(
                idx.reshape(-1), (pix * w[None]).reshape(-1), (p1 - p0) * nb
            ).reshape(-1, nb)

    _par(run, k, nt)
    return out


def _regional_max(a: np.ndarray) -> np.ndarray:
    """MATLAB's ``imregionalmax`` of a 2D image (8-connected; B08's port)"""
    from petakit5d.tools.puncta_removal import imregionalmax

    return np.asarray(imregionalmax(a), dtype=bool)


def _scalar(x) -> float:
    return float(np.asarray(x).reshape(-1)[0])


def _column_major(a: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(a.T).reshape(-1)


def _detect(im: np.ndarray, xy_pixel_size: float, dz: float, skew_angle: float) -> np.ndarray:
    """the peaks of the PSFs of ``im`` (z, y, x): rows ``[y x z]``, 1-based (MATLAB's
    ``pcoords`` before the isolation and border tests)"""
    from petakit5d.matlabcompat.stats import median, std

    nz, ny, nx = im.shape
    theta = float(skew_angle)
    theta_1 = _atand(_tand(theta) * xy_pixel_size / (dz * _sind(theta)))
    # the Radon transform of each (x, z) plane (MATLAB's squeeze(im(y, :, :))): ny x nbins
    planes = np.transpose(im, (1, 2, 0))
    im_radon = _radon(planes, theta_1)
    bw = _regional_max(im_radon)
    px, py = np.nonzero(bw.T)  # find(BW): column-major
    im_std = _scalar(std(_column_major(im_radon)))
    h, w = im_radon.shape
    # the local test: peak - median > 2 (0.95 std + 0.05 global std) of the 11 x 11 window
    # (fewer at the border). The median of a whole window is its 61st value (MATLAB's exactly);
    # NumPy's std (within a few ulps of MATLAB's) decides where the two sides are apart by more
    # than 1e-9 of their size, MATLAB's (matlabcompat) the rest and the border windows
    ys, xs = py.astype(np.int64), px.astype(np.int64)
    peak = im_radon[ys, xs]
    lhs = np.empty(ys.size)
    rhs = np.full(ys.size, np.inf)
    inner = (ys >= 5) & (ys + 6 <= h) & (xs >= 5) & (xs + 6 <= w)
    if inner.any():
        d = np.arange(-5, 6)
        # the windows column-major: x outer, y inner
        win = im_radon[
            ys[inner, None, None] + d[None, None, :], xs[inner, None, None] + d[None, :, None]
        ]
        win = win.reshape(-1, 121)
        lhs[inner] = peak[inner] - np.partition(win, 60, axis=1)[:, 60]
        rhs[inner] = 2 * (np.std(win, axis=1, ddof=1) * 0.95 + im_std * 0.05)
    with np.errstate(invalid="ignore"):
        exact = ~inner | ~(np.abs(lhs - rhs) > 1e-9 * (np.abs(lhs) + np.abs(rhs)))
    for i in np.flatnonzero(exact):
        yi, xi = int(ys[i]), int(xs[i])
        vec = _column_major(
            im_radon[max(0, yi - 5) : min(h, yi + 6), max(0, xi - 5) : min(w, xi + 6)]
        )
        lhs[i] = peak[i] - _scalar(median(vec))
        rhs[i] = 2 * (_scalar(std(vec)) * 0.95 + im_std * 0.05)
    kept = [(int(ys[i]), int(xs[i])) for i in np.flatnonzero(lhs > rhs)]

    # the brightest voxel of each peak's (x, z) plane within 10 pixels of its line
    xp = _radon_xp(nx, nz)
    k = _tand(theta_1 + 90)
    zg = (nx - np.arange(nx, dtype=np.float64))[:, None] - (nx + 1) // 2  # flip(Z, 1) - rc(1)
    xg = np.arange(1, nz + 1, dtype=np.float64)[None, :] - (nz + 1) // 2
    s1, c1 = _sind(theta_1), _cosd(theta_1)
    pcoords = np.zeros((len(kept), 3))
    for i, (yi, xi) in enumerate(kept):
        img = im[:, yi, :].T.astype(np.float64)  # (x, z)
        c = xp[xi] * s1 - xp[xi] * c1 * k
        line = np.abs(xg * k + c - zg) < 10
        _, (px_, pz_) = _first_max(np.ascontiguousarray(line * img))
        pcoords[i] = (yi + 1, px_ + 1, pz_ + 1)
    return pcoords


@register(
    "XR_psf_detection_and_cropping",
    Param("fn", None, "The calibration image (tiff), or its parts in order (along z)."),
    Param("resultDir", None, "Result folder (with a trailing separator)."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Scan interval in um."),
    Param("skewAngle", "@isnumeric", "Skew angle (in degree) of the stage."),
    Param("cropSize", "@isnumeric", "Size (z, y, x) of the crops.", convert="size_yxz"),
    Param(
        "distThresh",
        "@isnumeric",
        "Distance (z, y, x) within which another PSF makes a PSF not isolated.",
        convert="size_yxz",
    ),
    Param("prefix", "@ischar", "Prefix of the result files."),
    task="B06",
)
def psf_detection_and_cropping(
    fn: str | Sequence[str],
    result_dir: str,
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.1,
    skew_angle: float = 32.45,
    crop_size: Sequence[float] = (201, 256, 128),
    dist_thresh: Sequence[float] = (201, 256, 128),
    prefix: str = "test_",
) -> None:
    """Detection of the PSFs of a calibration image and crops of the isolated ones
    (``XR_psf_detection_and_cropping``; see the module): ``<result_dir><prefix>_crop_*.tif``,
    their z MIPs in ``MIPs/`` and ``<result_dir><prefix>.mat``, written last. Nothing if that
    exists."""
    import time

    from petakit5d.core.axes import to_yxz
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import read_tiff, write_mat, write_tiff
    from petakit5d.jobs.files import exist
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {
        k: to_plain(v)
        for k, v in matlab_results(psf_detection_and_cropping, dict(locals())).items()
    }
    t0 = time.perf_counter()
    info_fn = f"{result_dir}{prefix}.mat"
    if exist(info_fn):
        return
    fns = [fn] if isinstance(fn, (str, os.PathLike)) else list(fn)
    fns = [os.fspath(f) for f in fns]
    print(f"Start PSF detection for {fns[0]} ...", flush=True)
    ims = [read_tiff(f) for f in fns]
    ims = [a[None] if a.ndim == 2 else a for a in ims]
    im = np.concatenate(ims, axis=0) if len(ims) > 1 else ims[0]
    if im.dtype.kind == "f":
        im = np.where(np.isnan(im), im.dtype.type(0), im)
    im = np.ascontiguousarray(im)
    nz, ny, nx = im.shape

    pcoords = _detect(im, float(xy_pixel_size), float(dz), float(skew_angle))

    # isolated peaks away from the border (MATLAB's [y x z])
    dist = np.asarray(to_yxz(tuple(float(v) for v in dist_thresh)), dtype=np.float64)
    crop = np.asarray(to_yxz(tuple(float(v) for v in crop_size)), dtype=np.float64)
    keep = [int(np.sum(np.all(np.abs(pc - pcoords) <= dist, axis=1))) == 1 for pc in pcoords]
    pcoords = pcoords[np.asarray(keep, dtype=bool)]
    border = (crop - 1) / 2
    msize = np.array([ny, nx, nz], dtype=np.float64)
    pcoords = pcoords[np.all((pcoords >= border) & (msize - pcoords >= border), axis=1)]

    hsize = np.array([_round(v) for v in (crop - 1) / 2], dtype=np.int64)
    os.makedirs(result_dir, exist_ok=True)
    mip_dir = f"{result_dir}MIPs/"
    os.makedirs(mip_dir, exist_ok=True)
    included = np.ones(pcoords.shape[0], dtype=bool)
    for i, pc in enumerate(pcoords.astype(np.int64)):
        lo = pc - hsize
        hi = pc + hsize
        if np.any(lo < 1):
            raise PsfAnalysisError(
                "Index in position 1 is invalid. Array indices must be positive integers or "
                "logical values.",
                "MATLAB:index:expected_one_positive_integer",
            )
        im_i = np.ascontiguousarray(im[lo[2] - 1 : hi[2], lo[0] - 1 : hi[0], lo[1] - 1 : hi[1]])
        if float(im_i.max()) > _SATURATION:
            included[i] = False
            continue
        y, x, z = (int(v) for v in pc)
        name = f"{prefix}_crop_{y}_{x}_{z}"
        write_tiff(im_i, f"{result_dir}/{name}.tif")
        write_tiff(np.ascontiguousarray(im_i.max(axis=0)), f"{mip_dir}/{name}_MIP_z.tif")
    pcoords = pcoords[included]

    # saved under a temporary name and renamed: the job framework counts the task as done
    # once info_fn exists
    tmp = f"{info_fn[:-4]}_{get_uuid()}.mat"
    write_mat(tmp, {"pcoords": np.ascontiguousarray(pcoords.reshape(-1, 3)), "pr": pr})
    os.replace(tmp, info_fn)
    print("Done!", flush=True)
    print(f"Elapsed time is {time.perf_counter() - t0:.6f} seconds.", flush=True)


# MATLAB's regexp '^(?!.*_part[0-9]*.tif).*$': not a part file
_PART_RX = re.compile(r"^(?!.*_part[0-9]*.tif).*$")


@register(
    "XR_psf_detection_and_analysis_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders with calibration images."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Scan interval in um."),
    Param("skewAngle", "@isnumeric", "Skew angle (in degree) of the stage."),
    Param("cropSize", "@isnumeric", "Size (z, y, x) to crop isolated PSFs.", convert="size_yxz"),
    Param("flipZstack", "@islogical", "Flip z stacks."),
    Param(
        "distThresh",
        "@isnumeric",
        "Distance (z, y, x) between PSF centers below which PSFs are not isolated.",
        convert="size_yxz",
    ),
    Param("channelPatterns", "@iscell", "Channel identifiers for included channels."),
    Param("channels", "@isnumeric", "Wavelength (nm) of each channel."),
    Param("RWFn", "@iscell", "Richards & Wolf widefield PSF of each channel.", py="rw_fn"),
    Param("sourceStr", "@ischar", "Source of the PSF shown in the titles of the plots."),
    Param("visible", "@islogical", "Show the figures."),
    Param("parseCluster", "@islogical", "Use the Slurm cluster for the analysis."),
    Param("masterCompute", "@islogical", "This process computes tasks of the analysis too."),
    Param("cpusPerTask", "@isscalar", "CPUs per task."),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only)."),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="B06",
)
def psf_detection_and_analysis_dataset(
    data_paths: str | Sequence[str],
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.1,
    skew_angle: float = 32.45,
    crop_size: Sequence[float] = (201, 256, 128),
    flip_zstack: bool = False,
    dist_thresh: Sequence[float] = (201, 256, 128),
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0"),
    channels: Sequence[float] = (488, 560),
    rw_fn: Sequence[str] = ("", ""),
    source_str: str = "test",
    visible: bool = True,
    parse_cluster: bool = True,
    master_compute: bool = False,
    cpus_per_task: float = 8,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """PSF detection, cropping and analysis of calibration images
    (``XR_psf_detection_and_analysis_wrapper``; see the module): one
    ``psf_detection_and_cropping`` task per image (its parts ``_partNNN`` together) into
    ``<data>/Cropped/``, then ``psf_analysis_dataset`` of the crops (deskewed)."""
    from petakit5d.core.axes import to_yxz
    from petakit5d.core.paths import matlab_dir
    from petakit5d.core.text import through_text
    from petakit5d.io import get_image_size
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import contains_text

    paths = _paths(data_paths)
    print(paths, flush=True)
    prefixes: list[str] = []
    fns: list[str] = []
    for p in paths:
        names = [e.name for e in matlab_dir(f"{p}/*.tif")]
        prefixes.extend(sorted({n[:-4] for n in names if _PART_RX.match(n)}))
        fns.extend(f"{p}{n}" for n in names)
    keep = [False] * len(fns)
    for c in channel_patterns:
        keep = [a or b for a, b in zip(keep, contains_text(fns, str(c)), strict=True)]
    fns = [f for f, k in zip(fns, keep, strict=True) if k]

    def ten(v):
        return [float(through_text(x, "%.10f")) for x in to_yxz(tuple(float(x) for x in v))]

    crop_m, dist_m = ten(crop_size), ten(dist_thresh)
    frames, outs, tasks = [], [], []
    for prefix in prefixes:
        fns_i = [f for f in fns if prefix in f]
        if not fns_i:
            continue
        result = f"{os.path.dirname(fns_i[0])}/Cropped/"
        os.makedirs(result, exist_ok=True)
        frames.append(fns_i[0])
        outs.append(f"{result}{prefix}.mat")
        tasks.append(
            Task.of(
                psf_detection_and_cropping,
                fns_i,
                result,
                xy_pixel_size=_txt(xy_pixel_size),
                dz=_txt(dz),
                skew_angle=_txt(skew_angle),
                crop_size=[crop_m[2], crop_m[0], crop_m[1]],
                dist_thresh=[dist_m[2], dist_m[0], dist_m[1]],
                prefix=prefix,
            )
        )
    if not frames:
        raise PsfAnalysisError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    mem = float(np.prod(get_image_size(frames[0]), dtype=np.float64)) * 4 / 1024**3 * 100
    # (MATLAB's call: the framework's default parseCluster)
    opts = {
        "master_compute": True,
        "cpus_per_task": cpus_per_task,
        "mem_allocate": mem,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    if not all(framework.run_tasks(frames, outs, tasks, **opts)):
        framework.run_tasks(frames, outs, tasks, **opts)

    psf_analysis_dataset(
        [f"{p}/Cropped/" for p in paths],
        dz=dz,
        skew_angle=skew_angle,
        channel_patterns=channel_patterns,
        channels=channels,
        deskew=True,
        flip_zstack=flip_zstack,
        objective_scan=False,
        z_stage_scan=False,
        source_str=source_str,
        rw_fn=rw_fn,
        visible=visible,
        parse_cluster=parse_cluster,
        master_compute=master_compute,
        cpus_per_task=cpus_per_task,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )
