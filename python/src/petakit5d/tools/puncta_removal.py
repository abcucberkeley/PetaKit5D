"""Puncta removal for expansion microscopy (MATLAB's ``microscopeDataProcessing/puncta_removal``):
the zarr workflow (``XR_ExM_PunctaRemoval_zarr``), its block tasks (the updated method, the
interpolated background, the background estimation, and the first method) and the removal of
bright puncta by point detection (``XR_ExM_PunctaRemovalPointDetection``).

Arrays are ``(z, y, x)`` C order; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive;
sizes, batch and block sizes and window sizes are ``(z, y, x)``; batch numbers are 0-based. The
MATLAB-name layer converts MATLAB's values. Files keep MATLAB's names and contents:
``<data>/Puncta_Removed/<name>.zarr`` (native zarr) with its MIPs (``mip_zarr``),
``parameters.mat`` and ``parameters.json`` of the call, the block flags
``zarr_flag/<name>_<uuid>/blocks_<i>_<j>.mat``; the background estimation writes
``<data>/Puncta_background_estimation/<name>_batch_<y>_<x>_<z>_background_values.zarr`` (batch
size in MATLAB's ``y, x, z`` order, as MATLAB names it and the interpolated background reads
it) from the per-task results ``blocks_<i>_<j>.mat`` (``bg_mat``, ``batchInds``).

MATLAB's helper functions are written from their method and from MATLAB's measured behaviour
(oracle cases and probes; none of their code was read):

- ``threshold_otsu`` (``thresholdOtsu``): the values scaled to [0, 1] by the smallest positive
  value and the largest (``(x - min(x(x > 0))) / (max - min)``, below 0 counted as 0), put in
  256 bins (``round(255 * v)``), Otsu's threshold as the mean of the bins that maximize the
  between-class variance, scaled back (``min + level * (max - min)``, ``level = bin / 255``);
  one distinct value gives that value.
- ``locmax3d`` (``XR_locmax3d``): the values larger than every other value in the window
  (clipped at the borders) centred on them, else 0.
- ``imregionalmax`` (2-D, 8-connected plateaus with no higher neighbour).
- ``fillmissing2_linear`` (``fillmissing2(A, 'linear', MissingLocations=M)``): linear
  interpolation in MATLAB's Delaunay triangulation of the known pixels (points ``(row,
  column)``; where four or more are on one circle, the triangulation of the symbolic
  perturbation that lifts points in lexicographic order, so that the largest point of a
  cocircular polygon is cut off first), NaN outside their convex hull, NaN pixels left out of
  the known ones; fewer than three known pixels, or all on a line, leave the image as it is.
  Computed in double, rounded to the class.
- ``fillmissing2_nearest`` (``fillmissing2(A, 'nearest')``): each NaN pixel takes the value of
  the nearest pixel that isn't NaN (Euclidean distance; of equal ones, the last in column-major
  order); fewer than three known pixels, or all on a line, leave the image as it is.
- ``imresize3_linear`` (``imresize3(A, size, 'linear')`` for enlarging): per dimension, in the
  order of increasing scale, triangle weights normalized, in double, rounded to the class after
  each dimension; ``imresize_bilinear`` (``imresize(A, size, 'bilinear')`` for enlarging) is
  the same of an image.
- regionprops3's ``PrincipalAxisLength(:, 1) > 35``: the largest eigenvalue of the
  coordinates' covariance (normalized by the count) plus 1/12, decided in exact arithmetic.

MATLAB quirks kept (each has a test):

- the watershed of each punctum's region in the point detection doesn't change the result
  (its mask is computed after the voxels to remove were taken), so it is not computed;
- the threshold from the MIP mask (``T_mask``) is computed only when ``base_threshold`` is
  given, as only then it is used (MATLAB computes it always; nothing else depends on it);
- ``XR_ExM_PunctaRemoval_block`` (the first method, not used by the zarr workflow): an empty
  ``max_threshold`` doesn't bound ``T`` (``T > []`` is false).

parameters.json is written with ``write_json_file``, as MATLAB writes it from 5ecc0b7 on.
"""

from __future__ import annotations

import math
import os
import re
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "PunctaRemovalError",
    "fillmissing2_linear",
    "fillmissing2_nearest",
    "imregionalmax",
    "imresize3_linear",
    "imresize_bilinear",
    "locmax3d",
    "puncta_removal_background_estimation_block",
    "puncta_removal_block",
    "puncta_removal_interpolated_background_block",
    "puncta_removal_point_detection",
    "puncta_removal_updated_block",
    "puncta_removal_zarr",
    "threshold_otsu",
]


class PunctaRemovalError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# MATLAB's arithmetic
# ---------------------------------------------------------------------------------------------


def _f64(x: Any) -> np.float64:
    return np.float64(x)


def _times(a: Any, b: Any) -> Any:
    """``a .* b`` with MATLAB's classes (single with double: in double, rounded)"""
    from petakit5d.core.dtypes import times

    return times(a, b)


def _plus(a: Any, b: Any) -> Any:
    from petakit5d.core.dtypes import plus

    return plus(a, b)


def _scalar(x: Any) -> Any:
    """a 0-d result as a NumPy scalar of its class"""
    a = np.asarray(x)
    return a.reshape(()).dtype.type(a.reshape(-1)[0]) if a.size == 1 else a


def _concat_max(a: float, b: Any) -> Any:
    """``max([a, b])`` with ``a`` double and ``b`` of a class: MATLAB's concatenation converts
    ``a`` to ``b``'s class (single or an integer class: rounded, saturated)"""
    from petakit5d.core.dtypes import cast
    from petakit5d.matlabcompat import reduce as mred

    b = _scalar(b)
    cls = np.asarray(b).dtype
    if cls == np.float64:
        v = np.array([a, b], dtype=np.float64)
    else:
        from petakit5d.core.arrays import matlab_class

        v = np.array([_scalar(cast(np.float64(a), matlab_class(np.asarray(b)))), b], dtype=cls)
    return _scalar(mred.max(v)[0])


def _mmax(a: Any, b: Any) -> Any:
    """MATLAB's ``max(a, b)`` of two scalars (NaN ignored; the class of mixed operands: single
    or the integer class, the double one converted to it)"""
    from petakit5d.core.dtypes import cast, result_class

    cls = result_class(np.asarray(a), np.asarray(b))
    x = _scalar(cast(np.asarray(a), cls))
    y = _scalar(cast(np.asarray(b), cls))
    if np.isnan(x):
        return y
    if np.isnan(y):
        return x
    return y if y > x else x


def _mmin(a: Any, b: Any) -> Any:
    """MATLAB's ``min(a, b)`` of two scalars"""
    from petakit5d.core.dtypes import cast, result_class

    cls = result_class(np.asarray(a), np.asarray(b))
    x = _scalar(cast(np.asarray(a), cls))
    y = _scalar(cast(np.asarray(b), cls))
    if np.isnan(x):
        return y
    if np.isnan(y):
        return x
    return y if y < x else x


def _min_array(a: np.ndarray, bound: float) -> np.ndarray:
    """``min(a, bound)`` of an array and a double scalar (the array's class)"""
    if math.isinf(bound) and bound > 0:
        return a
    from petakit5d.core.arrays import matlab_class
    from petakit5d.core.dtypes import cast

    b = _scalar(cast(np.float64(bound), matlab_class(a)))
    out = np.minimum(a, b)
    if a.dtype.kind == "f":
        nan = np.isnan(a)
        if nan.any():
            out[nan] = b
    return out


def _vol3(a: np.ndarray) -> np.ndarray:
    return a[np.newaxis] if a.ndim == 2 else a


def _ints(v: Any) -> list[int]:
    out = []
    for x in np.asarray(v).ravel():
        f = float(x)
        if not math.isfinite(f) or f != math.floor(f):
            raise ValueError(f"expected integers, got {v!r}")
        out.append(int(f))
    return out


def _size3(v: Any, name: str) -> list[int]:
    s = _ints(v)
    if len(s) == 1:
        s = s * 3
    if len(s) != 3:
        raise ValueError(f"{name} has 3 values (z, y, x), got {len(s)}")
    return s


def _nums(v: Any) -> list[float]:
    return [float(x) for x in np.atleast_1d(np.asarray(v, dtype=np.float64)).ravel()]


def _empty(v: Any) -> bool:
    return v is None or np.asarray(v).size == 0


def _toc(t0: float) -> None:
    print(f"Elapsed time is {time.perf_counter() - t0:.6f} seconds.", flush=True)


# ---------------------------------------------------------------------------------------------
# thresholdOtsu, XR_locmax3d, imregionalmax
# ---------------------------------------------------------------------------------------------


def _otsu_bins(t: np.ndarray, lo: float, r: float) -> np.ndarray:
    """the bins of values (float64, changed in place): ``round(255 * (t - lo) / r)`` in [0, 255]"""
    t -= lo
    t /= r
    t *= 255
    t += 0.5
    np.floor(t, out=t)
    np.clip(t, 0, 255, out=t)
    return t.astype(np.intp)


def _otsu_counts(v: np.ndarray) -> tuple[np.ndarray, float, float, int]:
    """the 256 bin counts of thresholdOtsu, the smallest positive value, the largest, the count"""
    if v.dtype.kind == "f" and np.isnan(v).any():
        v = v[~np.isnan(v)]
    pos = v[v > 0]
    if pos.size == 0:
        raise PunctaRemovalError("thresholdOtsu takes positive values", "MATLAB:innerdim")
    lo = float(pos.min())
    hi = float(v.max())
    counts = np.zeros(256, dtype=np.float64)
    if hi == lo:
        return counts, lo, hi, v.size
    step = 1 << 20
    for k in range(0, v.size, step):
        b = _otsu_bins(v[k : k + step].astype(np.float64), lo, hi - lo)
        counts += np.bincount(b, minlength=256)
    return counts, lo, hi, v.size


def _chunks(n: int, nt: int) -> list[slice]:
    k = max(1, min(4 * nt, n // (1 << 18)))
    step = -(-n // k)
    return [slice(i, min(n, i + step)) for i in range(0, n, step)]


def _bincount(a: np.ndarray, minlength: int, num_threads: int | None) -> np.ndarray:
    """``np.bincount`` of a flat non-negative integer array in threads (NumPy's releases the
    GIL); integer counts, so the result doesn't depend on the split"""
    from petakit5d.core.threads import resolve_threads

    nt = resolve_threads(num_threads)
    if nt == 1 or a.size < 1 << 20:
        return np.bincount(a, minlength=minlength)
    from concurrent.futures import ThreadPoolExecutor

    parts = _chunks(a.size, nt)
    with ThreadPoolExecutor(nt) as ex:
        res = list(ex.map(lambda sl: np.bincount(a[sl], minlength=minlength), parts))
    m = max(r.size for r in res)
    out = np.zeros(m, dtype=np.int64)
    for r in res:
        out[: r.size] += r
    return out


def _take(table: np.ndarray, idx: np.ndarray, num_threads: int | None) -> np.ndarray:
    """``table[idx]`` of a flat index array in threads"""
    from petakit5d.core.threads import resolve_threads

    nt = resolve_threads(num_threads)
    out = np.empty(idx.shape, dtype=table.dtype)
    if nt == 1 or idx.size < 1 << 20:
        np.take(table, idx, out=out)
        return out
    from concurrent.futures import ThreadPoolExecutor

    fi, fo = idx.reshape(-1), out.reshape(-1)
    with ThreadPoolExecutor(nt) as ex:
        list(ex.map(lambda sl: np.take(table, fi[sl], out=fo[sl]), _chunks(fi.size, nt)))
    return out


def _area_open(bw: np.ndarray, p: float, num_threads: int | None) -> np.ndarray:
    """``bwareaopen(bw, p, 26)`` of a volume: its 26-connected components of at least ``p``
    voxels (matlabcompat.label's labels; the sizes and the mask in threads)"""
    from petakit5d.matlabcompat import label

    if not (p >= 0 and p == int(p)):
        return label.bwareaopen(bw, p, 26, threads=num_threads)
    lab, n = label._label(bw, 26, "any", num_threads)
    sizes = _bincount(lab.reshape(-1), n + 1, num_threads)
    keep = sizes >= p
    keep[0] = False
    return _take(keep, lab, num_threads).reshape(bw.shape)


def _value_counts(v: np.ndarray, num_threads: int | None = 1) -> tuple[np.ndarray, np.ndarray]:
    """the distinct values of 8- or 16-bit integers (ascending, int64) and their counts"""
    base = int(np.iinfo(v.dtype).min) if v.dtype.kind == "i" else 0
    flat = (v.ravel().astype(np.int64) - base) if base else v.ravel()
    raw = _bincount(flat, 1, num_threads)
    idx = np.nonzero(raw)[0]
    return idx + base, raw[idx]


def _otsu_counts_values(vals: np.ndarray, cnt: np.ndarray) -> tuple[np.ndarray, float, float, int]:
    """``_otsu_counts`` of integer values given as distinct values and their counts: the bins
    of the values present (the same bins)"""
    pos = vals[vals > 0]
    if pos.size == 0:
        raise PunctaRemovalError("thresholdOtsu takes positive values", "MATLAB:innerdim")
    lo = float(pos[0])
    hi = float(vals[-1])
    n = int(cnt.sum())
    counts = np.zeros(256, dtype=np.float64)
    if hi == lo:
        return counts, lo, hi, n
    b = _otsu_bins(vals.astype(np.float64), lo, hi - lo)
    counts += np.bincount(b, weights=cnt.astype(np.float64), minlength=256)
    return counts, lo, hi, n


def _otsu_from_counts(counts: np.ndarray, lo: float, hi: float, n: int) -> float:
    if hi == lo:
        return lo
    p = counts / float(n)
    omega = np.cumsum(p)
    mu = np.cumsum(p * np.arange(1, 257, dtype=np.float64))
    with np.errstate(divide="ignore", invalid="ignore"):
        sb = (mu[-1] * omega - mu) ** 2 / (omega * (1 - omega))
    good = ~np.isnan(sb)
    if not good.any():
        return lo
    m = sb[good].max()
    idx = float(np.mean(np.nonzero(sb == m)[0] + 1.0))
    level = (idx - 1) / 255
    return level * (hi - lo) + lo


def threshold_otsu(x: Any) -> float:
    """MATLAB's ``thresholdOtsu`` (module docstring): Otsu's threshold of the values (any order,
    any class) as a double. No values raise ``MATLAB:innerdim``, as MATLAB does."""
    v = np.asarray(x).ravel()
    if v.size == 0:
        raise PunctaRemovalError("Inner matrix dimensions must agree.", "MATLAB:innerdim")
    if v.dtype.kind in "uib" and v.dtype.itemsize <= 2:
        return _otsu_from_counts(*_otsu_counts_values(*_value_counts(v)))
    return _otsu_from_counts(*_otsu_counts(v))


def locmax3d(vol: Any, window: Any = 5) -> np.ndarray:
    """MATLAB's ``XR_locmax3d(vol, [w w w])`` (``w`` odd, the same along each axis): ``vol``
    where a value is larger than all others in the ``w^3`` window centred on it (clipped at the
    borders), else 0."""
    from scipy import ndimage as ndi

    a = np.asarray(vol)
    w = _ints(window)
    if len(set(w)) != 1 or w[0] % 2 != 1:
        raise ValueError("locmax3d takes an odd window, the same along each axis")
    w0 = w[0]
    h = w0 // 2
    out = np.zeros(a.shape, dtype=a.dtype)
    if a.size == 0:
        return out
    m = ndi.maximum_filter(a, size=w0, mode="nearest")
    cand = np.nonzero((a == m) & ~np.isnan(a) if a.dtype.kind == "f" else a == m)
    if cand[0].size == 0:
        return out
    # a candidate is the maximum; it is the only one if no other voxel of its window equals it
    val = a[cand]
    unique = np.ones(val.size, dtype=bool)
    shape = a.shape
    for dz in range(-h, h + 1):
        for dy in range(-h, h + 1):
            for dx in range(-h, h + 1):
                if dz == 0 and dy == 0 and dx == 0:
                    continue
                zz, yy, xx = cand[0] + dz, cand[1] + dy, cand[2] + dx
                ok = (
                    (zz >= 0)
                    & (zz < shape[0])
                    & (yy >= 0)
                    & (yy < shape[1])
                    & (xx >= 0)
                    & (xx < shape[2])
                    & unique
                )
                if not ok.any():
                    continue
                k = np.nonzero(ok)[0]
                unique[k[a[zz[k], yy[k], xx[k]] == val[k]]] = False
    keep = tuple(c[unique] for c in cand)
    out[keep] = a[keep]
    return out


def imregionalmax(im: Any) -> np.ndarray:
    """MATLAB's ``imregionalmax`` of an image ``(y, x)`` (8-connected): the plateaus of equal
    values without a higher neighbour. NaN raises ``MATLAB:imregionalmax:expectedNonNaN``."""
    from scipy import ndimage as ndi

    a = np.asarray(im)
    if a.dtype.kind == "f" and np.isnan(a).any():
        raise PunctaRemovalError(
            "Expected input number 1, I, to be non-NaN.", "MATLAB:imregionalmax:expectedNonNaN"
        )
    if a.size == 0:
        return np.zeros(a.shape, dtype=bool)
    m = ndi.maximum_filter(a, size=3, mode="nearest")
    cand = a == m
    # a plateau with a voxel that has a higher neighbour is not a maximum: a candidate next to
    # an equal non-candidate is in such a plateau (adjacent candidates have equal values)
    ny, nx = a.shape
    bad = np.zeros(a.shape, dtype=bool)
    non = ~cand
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0:
                continue
            ys = slice(max(0, dy), ny + min(0, dy))
            yd = slice(max(0, -dy), ny + min(0, -dy))
            xs = slice(max(0, dx), nx + min(0, dx))
            xd = slice(max(0, -dx), nx + min(0, -dx))
            bad[yd, xd] |= cand[yd, xd] & non[ys, xs] & (a[ys, xs] == a[yd, xd])
    lab, n = ndi.label(cand, structure=np.ones((3, 3), dtype=bool))
    if n == 0:
        return cand
    worse = np.zeros(n + 1, dtype=bool)
    worse[np.unique(lab[bad])] = True
    worse[0] = True
    return ~worse[lab]


# ---------------------------------------------------------------------------------------------
# fillmissing2(A, 'linear', MissingLocations=M)
# ---------------------------------------------------------------------------------------------


def _orient(a: tuple[int, int], b: tuple[int, int], c: tuple[int, int]) -> int:
    v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (v > 0) - (v < 0)


def _cocircular(a, b, c, d) -> bool:
    """d on the circle through a, b, c (exact, integers)"""
    adx, ady = a[0] - d[0], a[1] - d[1]
    bdx, bdy = b[0] - d[0], b[1] - d[1]
    cdx, cdy = c[0] - d[0], c[1] - d[1]
    ad = adx * adx + ady * ady
    bd = bdx * bdx + bdy * bdy
    cd = cdx * cdx + cdy * cdy
    return (
        adx * (bdy * cd - bd * cdy) - ady * (bdx * cd - bd * cdx) + ad * (bdx * cdy - bdy * cdx)
        == 0
    )


def _cell_triangles(verts: list[tuple[int, int]]) -> list[tuple]:
    """the triangles of a convex polygon of cocircular points under the symbolic perturbation:
    the lexicographically largest vertex is cut off as an ear, repeatedly"""
    cy = sum(p[0] for p in verts) / len(verts)
    cx = sum(p[1] for p in verts) / len(verts)
    poly = sorted(verts, key=lambda p: math.atan2(p[1] - cx, p[0] - cy))
    tris = []
    while len(poly) > 3:
        k = max(range(len(poly)), key=lambda i: poly[i])
        n = len(poly)
        tris.append((poly[k - 1], poly[k], poly[(k + 1) % n]))
        del poly[k]
    tris.append(tuple(poly))
    return tris


def _interp(tris, q, a64: np.ndarray) -> float | None:
    for t in tris:
        A, B, C = t
        if _orient(A, B, C) < 0:
            B, C = C, B
        dA = (B[0] - q[0]) * (C[1] - q[1]) - (B[1] - q[1]) * (C[0] - q[0])
        dB = (C[0] - q[0]) * (A[1] - q[1]) - (C[1] - q[1]) * (A[0] - q[0])
        dC = (A[0] - q[0]) * (B[1] - q[1]) - (A[1] - q[1]) * (B[0] - q[0])
        if dA >= 0 and dB >= 0 and dC >= 0:
            D = dA + dB + dC
            return (dA / D) * a64[A] + (dB / D) * a64[B] + (dC / D) * a64[C]
    return None


def fillmissing2_linear(a: Any, missing: Any) -> np.ndarray:
    """MATLAB's ``fillmissing2(A, 'linear', MissingLocations=M)`` of an image ``(y, x)``
    (module docstring): the pixels of ``M`` filled by linear interpolation of the others."""
    from petakit5d.core.dtypes import NAN

    a = np.asarray(a)
    m = np.asarray(missing, dtype=bool)
    out = a.copy()
    if not m.any():
        return out
    known = ~m
    if a.dtype.kind == "f":
        known &= ~np.isnan(a)
    rr, cc = np.nonzero(known)
    if rr.size < 3:
        return out
    from scipy.spatial import Delaunay, QhullError

    pts = np.column_stack([rr, cc]).astype(np.float64)
    try:
        tri = Delaunay(pts)
    except QhullError:
        return out
    simp = tri.simplices
    nbrs = tri.neighbors
    qr, qc = np.nonzero(m)
    qs = tri.find_simplex(np.column_stack([qr, qc]).astype(np.float64))
    a64 = a.astype(np.float64)
    vals = np.full(qr.size, NAN, dtype=np.float64)
    cells: dict[int, list[tuple]] = {}

    def vert(j: int) -> tuple[int, int]:
        return (int(rr[j]), int(cc[j]))

    for i in range(qr.size):
        s = int(qs[i])
        if s < 0:
            continue
        tris = cells.get(s)
        if tris is None:
            # s's cell: the triangles reached through edges whose far vertex is on s's circle
            t0 = [vert(j) for j in simp[s]]
            seen = {s}
            stack = [s]
            vs: set[int] = set()
            while stack:
                t = stack.pop()
                vs.update(int(j) for j in simp[t])
                for _e, nb in enumerate(nbrs[t]):
                    if nb < 0 or nb in seen:
                        continue
                    far = [int(j) for j in simp[nb] if j not in simp[t]][0]
                    if _cocircular(t0[0], t0[1], t0[2], vert(far)):
                        seen.add(int(nb))
                        stack.append(int(nb))
            tris = [tuple(t0)] if len(vs) == 3 else _cell_triangles([vert(j) for j in vs])
            for t in seen:
                cells[t] = tris
        v = _interp(tris, (int(qr[i]), int(qc[i])), a64)
        if v is not None:
            vals[i] = v
    if a.dtype.kind == "f":
        out[qr, qc] = vals.astype(a.dtype)
    else:
        from petakit5d.core.arrays import matlab_class
        from petakit5d.core.dtypes import cast

        out[qr, qc] = cast(vals, matlab_class(a))
    return out


def fillmissing2_nearest(a: Any) -> np.ndarray:
    """MATLAB's ``fillmissing2(A, 'nearest')`` of an image ``(y, x)`` (module docstring): each
    NaN pixel takes the value of the nearest pixel that isn't NaN (Euclidean distance; of equal
    ones, the last in MATLAB's column-major order); with fewer than three such pixels, or all on
    a line, the image is kept."""
    a = np.asarray(a)
    out = a.copy()
    if a.dtype.kind != "f":
        return out
    miss = np.isnan(a)
    # the known pixels in column-major order, reversed, so that argmin's first minimum is the
    # last of equal ones
    kc, kr = (v[::-1].astype(np.int64) for v in np.nonzero(~miss.T))
    if kr.size < 3 or not miss.any():
        return out
    dr, dc = kr - kr[0], kc - kc[0]
    j = int(np.flatnonzero((dr != 0) | (dc != 0))[0])
    if not np.any(dr[j] * dc - dc[j] * dr):
        # (all on a line: no triangle)
        return out
    qr, qc = np.nonzero(miss)
    vals = a[kr, kc]
    step = max(1, 2**22 // kr.size)
    for s in range(0, qr.size, step):
        r = qr[s : s + step].astype(np.int64)[:, None]
        c = qc[s : s + step].astype(np.int64)[:, None]
        d2 = (r - kr[None, :]) ** 2 + (c - kc[None, :]) ** 2
        out[qr[s : s + step], qc[s : s + step]] = vals[np.argmin(d2, axis=1)]
    return out


# ---------------------------------------------------------------------------------------------
# imresize3(A, size, 'linear') (enlarging)
# ---------------------------------------------------------------------------------------------


def _linear_weights(n_in: int, n_out: int) -> tuple[np.ndarray, np.ndarray]:
    s = n_out / n_in
    x = np.arange(1, n_out + 1, dtype=np.float64)
    u = x / s + 0.5 * (1 - 1 / s)
    left = np.floor(u - 1)
    idx = left[:, None] + np.arange(3, dtype=np.float64)[None, :]
    w = np.maximum(0.0, 1 - np.abs(u[:, None] - idx))
    w = w / np.sum(w, axis=1, keepdims=True)
    ind = (np.clip(idx, 1, n_in) - 1).astype(np.intp)
    return w, ind


def imresize3_linear(a: Any, size: Sequence[int]) -> np.ndarray:
    """MATLAB's ``imresize3(A, [ny nx nz], 'linear')`` of a volume ``(z, y, x)`` to ``size``
    ``(z, y, x)``, no smaller than ``A`` along any axis (module docstring)."""
    v = np.asarray(a)
    out_size = _size3(size, "size")
    if v.ndim != 3:
        raise ValueError("imresize3_linear takes a volume (z, y, x)")
    if any(o < n for o, n in zip(out_size, v.shape, strict=True)):
        raise ValueError("imresize3_linear only enlarges")
    # MATLAB's dimensions (y, x, z) in the order of increasing scale (stable)
    dims = [(-2, 1), (-1, 2), (-3, 0)]
    scales = [out_size[k] / v.shape[k] for _, k in dims]
    order = sorted(range(3), key=lambda i: scales[i])
    cls = v.dtype
    b = v
    for i in order:
        ax, k = dims[i]
        w, ind = _linear_weights(b.shape[ax], out_size[k])
        src = np.moveaxis(b, ax, 0)
        acc = np.zeros((out_size[k],) + src.shape[1:], dtype=np.float64)
        shp = (-1,) + (1,) * (src.ndim - 1)
        for j in range(w.shape[1]):
            acc += w[:, j].reshape(shp) * src[ind[:, j]].astype(np.float64)
        if cls.kind == "f":
            r = acc.astype(cls)
        else:
            from petakit5d.core.arrays import matlab_class
            from petakit5d.core.dtypes import cast

            r = cast(acc, matlab_class(v))
        b = np.ascontiguousarray(np.moveaxis(r, 0, ax))
    return b


def imresize_bilinear(a: Any, size: Sequence[int]) -> np.ndarray:
    """MATLAB's ``imresize(A, [ny nx], 'bilinear')`` of an image ``(y, x)`` to ``size``
    ``(y, x)``, no smaller than ``A`` along either axis: ``imresize3_linear`` of the image as
    one slice"""
    v = np.asarray(a)
    if v.ndim != 2:
        raise ValueError("imresize_bilinear takes an image (y, x)")
    ny, nx = (int(s) for s in size)
    return imresize3_linear(v[np.newaxis], (1, ny, nx))[0]


# ---------------------------------------------------------------------------------------------
# point detection
# ---------------------------------------------------------------------------------------------


def _bbox_window(c: Sequence[int], h: Sequence[int], shape: Sequence[int]) -> tuple[slice, ...]:
    return tuple(slice(max(0, c[d] - h[d]), min(shape[d], c[d] + h[d] + 1)) for d in range(len(c)))


def _matlab_values(a: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """``a(:)`` (or ``a(mask)``) in MATLAB's linear order"""
    from petakit5d.core.arrays import matlab_order

    v = matlab_order(a)
    if mask is None:
        return np.ascontiguousarray(v).reshape(-1)
    return v[matlab_order(mask)]


def _gt_scalar(a: np.ndarray, t: Any) -> np.ndarray:
    """MATLAB's ``a > t`` of an array and a scalar: compared exactly"""
    t = _scalar(t)
    if np.asarray(t).dtype == a.dtype:
        return a > t
    if a.dtype.kind in "ui" and math.isfinite(float(t)):
        # (integers above t are those above floor(t))
        info = np.iinfo(a.dtype)
        f = math.floor(float(t))
        if f < info.min:
            return np.ones(a.shape, dtype=bool)
        if f >= info.max:
            return np.zeros(a.shape, dtype=bool)
        return a > a.dtype.type(f)
    return a > np.float64(t)


def _component_at(bw: np.ndarray, p: tuple[int, ...]) -> np.ndarray | None:
    """the 26-connected component of ``bw`` that holds ``p`` (None if ``bw[p]`` is false)"""
    from scipy import ndimage as ndi

    if not bw[p]:
        return None
    lab, _ = ndi.label(bw, structure=np.ones((3,) * bw.ndim, dtype=bool))
    return lab == lab[p]


def _large_components(bw: np.ndarray, det_vol_thrsh: float) -> np.ndarray:
    """the voxels of the 26-connected components of ``bw`` with more than ``2 *
    det_vol_thrsh`` voxels or a principal axis longer than 35 (regionprops3's
    ``PrincipalAxisLength(:, 1)``: ``2 * sqrt(5 * (lambda + 1/12))``, lambda the largest
    eigenvalue of the coordinates' covariance normalized by the count), decided exactly"""
    from fractions import Fraction

    from scipy import ndimage as ndi

    lab, n = ndi.label(bw, structure=np.ones((3, 3, 3), dtype=bool))
    if n == 0:
        return np.zeros(bw.shape, dtype=bool)
    idx = np.nonzero(lab)
    lv = lab[idx]
    cnt = np.bincount(lv, minlength=n + 1)
    coords = [c.astype(np.int64) for c in idx]
    s1 = [np.bincount(lv, weights=c.astype(np.float64), minlength=n + 1) for c in coords]
    s2 = {}
    for i in range(3):
        for j in range(i, 3):
            s2[i, j] = np.bincount(
                lv, weights=(coords[i] * coords[j]).astype(np.float64), minlength=n + 1
            )
    remove = np.zeros(n + 1, dtype=bool)
    lim = 2 * det_vol_thrsh
    # PAL > 35 <=> lambda > 1225/20 - 1/12 = 367/6
    c0 = Fraction(367, 6)
    for k in range(1, n + 1):
        nk = int(cnt[k])
        if nk > lim:
            remove[k] = True
            continue
        # n^2 * covariance (exact integers): n * sum(xi xj) - sum(xi) sum(xj)
        S = [int(s1[i][k]) for i in range(3)]
        M = [[0] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(i, 3):
                M[i][j] = M[j][i] = nk * int(s2[i, j][k]) - S[i] * S[j]
        # lambda_max(M / n^2) > c0 <=> c0 * n^2 * I - M is not positive semidefinite
        c = c0 * nk * nk
        A = [[(c if i == j else 0) - M[i][j] for j in range(3)] for i in range(3)]
        minors = [A[0][0], A[1][1], A[2][2]]
        minors += [
            A[0][0] * A[1][1] - A[0][1] * A[1][0],
            A[0][0] * A[2][2] - A[0][2] * A[2][0],
            A[1][1] * A[2][2] - A[1][2] * A[2][1],
        ]
        det = (
            A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
            - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
            + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
        )
        minors.append(det)
        if any(mv < 0 for mv in minors):
            remove[k] = True
    return remove[lab]


def _detect_2d(vol_sm: np.ndarray, hbbox: Sequence[int], sigma_thrsh: float):
    from petakit5d.matlabcompat import stats

    if vol_sm.dtype.kind == "f" and np.isnan(vol_sm).any():
        from petakit5d.matlabcompat import reduce as mred

        mip, zinds = mred.max(vol_sm, axis=-3)
    else:
        mip = vol_sm.max(axis=0)
        zinds = None
    bw = imregionalmax(mip)
    xs, ys = np.nonzero(bw.T)
    if zinds is None:
        # the first z of each peak's maximum (MATLAB's index of max along z)
        zs = np.argmax(vol_sm[:, ys, xs], axis=0)
    else:
        zs = zinds[ys, xs]
    kept = np.zeros(ys.size, dtype=bool)
    im_std = _scalar(stats.std(mip))
    ref = _times(im_std, 0.05)
    ny, nx = mip.shape
    hy, hx = int(hbbox[1]), int(hbbox[2])
    for i in range(ys.size):
        y, x = int(ys[i]), int(xs[i])
        win = mip[max(0, y - hy) : min(ny, y + hy + 1), max(0, x - hx) : min(nx, x + hx + 1)]
        peak = mip[y, x]
        bg = _scalar(stats.median(win))
        sd = _scalar(stats.std(win))
        lhs = _scalar(peak - bg)
        rhs = _scalar(_times(sigma_thrsh, _plus(_times(sd, 0.95), ref)))
        if lhs > rhs:
            kept[i] = True
    return zs, ys, xs, kept


def _detect_3d(vol_sm: np.ndarray, hbbox: Sequence[int], sigma_thrsh: float):
    from petakit5d.core.arrays import matlab_order
    from petakit5d.matlabcompat import stats

    lm = locmax3d(vol_sm, 5)
    # find(BW > 0) in MATLAB's order: (z, x, y) view, y fastest
    zs, xs, ys = np.nonzero(matlab_order(lm) > 0)
    kept = np.zeros(zs.size, dtype=bool)
    im_std = _scalar(stats.std(_matlab_values(vol_sm, vol_sm > 0)))
    ref = _times(im_std, 0.05)
    for i in range(zs.size):
        p = (int(zs[i]), int(ys[i]), int(xs[i]))
        win = vol_sm[_bbox_window(p, hbbox, vol_sm.shape)]
        vals = _matlab_values(win, win > 0)
        peak = vol_sm[p]
        bg = _scalar(stats.median(vals))
        sd = _scalar(stats.std(vals))
        lhs = _scalar(peak - bg)
        rhs = _scalar(_times(sigma_thrsh, _plus(_times(sd, 0.95), ref)))
        if lhs > rhs:
            kept[i] = True
    return zs, ys, xs, kept


def _fill_slices(vol: np.ndarray, bw_s: np.ndarray, border: int = 3) -> None:
    """the pixels of ``bw_s`` filled slice by slice, each 8-connected component's pixels in a
    window of its box with ``border`` around it; later components' fills win (MATLAB's
    assignment), all from the slice before filling"""
    from scipy import ndimage as ndi

    ny, nx = vol.shape[1:]
    st = np.ones((3, 3), dtype=bool)
    for z in range(vol.shape[0]):
        bz = bw_s[z]
        if not bz.any():
            continue
        # components numbered in MATLAB's order (column-major): label the transpose
        lab, n = ndi.label(bz.T, structure=st)
        if n == 0:
            continue
        objs = ndi.find_objects(lab)
        plane = vol[z]
        fills = []
        for sl in objs:
            xs_, ys_ = sl
            y0, y1 = max(0, ys_.start - border), min(ny, ys_.stop + border)
            x0, x1 = max(0, xs_.start - border), min(nx, xs_.stop + border)
            win = plane[y0:y1, x0:x1]
            m = bz[y0:y1, x0:x1]
            f = fillmissing2_linear(win, m)
            # pixels outside the convex hull of the known ones are left NaN: they take the
            # value of the nearest pixel
            if np.isnan(f[m]).any():
                f = fillmissing2_nearest(f)
            fills.append((y0, x0, m, f[m]))
        out = plane.copy()
        for y0, x0, m, v in fills:
            yy, xx = np.nonzero(m)
            out[yy + y0, xx + x0] = v
        vol[z] = out


def _point_detection(
    vol: np.ndarray,
    sigma: float,
    local_win_size: Sequence[int],
    sigma_thrsh: float,
    int_thrsh: float,
    init_detect: str,
    det_vol_thrsh: float,
    num_threads: int | None,
) -> np.ndarray:
    """XR_ExM_PunctaRemovalPointDetection on ``vol`` (changed in place and returned)"""
    from petakit5d.core.dtypes import mround
    from petakit5d.matlabcompat import stats
    from petakit5d.matlabcompat.filters import imgaussfilt3

    hbbox = [int(mround((float(w) - 1) / 2)) for w in local_win_size]
    vol_sm = imgaussfilt3(vol, sigma, num_threads=num_threads)
    if init_detect == "2d":
        zs, ys, xs, kept = _detect_2d(vol_sm, hbbox, sigma_thrsh)
    elif init_detect == "3d":
        zs, ys, xs, kept = _detect_3d(vol_sm, hbbox, sigma_thrsh)
    else:
        # (MATLAB: kept_flag undefined)
        raise PunctaRemovalError(
            "Unrecognized function or variable 'kept_flag'.", "MATLAB:UndefinedFunction"
        )
    if not kept.any():
        return vol

    sz = vol.shape
    rm = np.zeros(sz, dtype=bool)
    any_rm = False
    for i in np.nonzero(kept)[0]:
        t0 = time.perf_counter()
        p = (int(zs[i]), int(ys[i]), int(xs[i]))
        peak = _mmax(vol_sm[p], vol[p])
        if float(peak) < float(int_thrsh):
            continue
        sl = _bbox_window(p, hbbox, sz)
        pl = tuple(p[d] - sl[d].start for d in range(3))
        im_i = vol_sm[sl]
        t_i = threshold_otsu(im_i[im_i > 0])
        med = _scalar(stats.median(im_i))
        t_s = _scalar(_plus(min(t_i, float(int_thrsh)) * 0.2, _times(med, 0.8)))
        bw = _component_at(_gt_scalar(im_i, t_s), pl)
        if bw is None:
            continue
        while np.count_nonzero(bw) > det_vol_thrsh:
            t_s = _scalar(_times(t_s, 1.1))
            bw = _component_at(_gt_scalar(im_i, t_s), pl)
            if bw is None:
                break
        if bw is None:
            continue
        # (MATLAB's watershed and its mask come after the voxels are taken: not computed)
        rm[sl] |= bw
        any_rm = True
        _toc(t0)
    if not any_rm:
        return vol
    del vol_sm

    rm[_large_components(rm, det_vol_thrsh)] = False
    _fill_slices(vol, rm)
    return vol


@register(
    "XR_ExM_PunctaRemovalPointDetection",
    Param("vol", "@isnumeric", "Volume (z, y, x).", convert="volume"),
    Param("Sigma", "@isnumeric", "Sigma of the 3D Gaussian filter."),
    Param("OTSUMaxPer", "@isnumeric", "Not used."),
    Param(
        "localWinSize",
        "@isvector",
        "Window (z, y, x) around a peak for its local statistics and its segmentation.",
        convert="size_yxz",
    ),
    Param("SigmaThrsh", "@isvector", "A peak is kept above its background by this many sigmas."),
    Param("intThrsh", "@isvector", "Peaks below this intensity are not removed."),
    Param("initDetect", "@ischar", "Initial detection: '2d' (z MIP) or '3d' (the volume)."),
    Param("detVolThrsh", "@isnumeric", "Largest volume of a punctum's region."),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads of the Gaussian filter (None: all).", python_only=True),
    task="B08",
)
def puncta_removal_point_detection(
    vol: np.ndarray,
    *,
    sigma: float = 2.5,
    otsu_max_per: float = 99.9,
    local_win_size: Sequence[int] = (41, 41, 41),
    sigma_thrsh: float = 4,
    int_thrsh: float = 5000,
    init_detect: str = "2d",
    det_vol_thrsh: float = 5000,
    debug: bool = False,
    num_threads: int | None = None,
) -> np.ndarray:
    """``XR_ExM_PunctaRemovalPointDetection``: bright puncta of ``vol`` (a ``(z, y, x)`` single
    or double volume) found as peaks of its Gaussian filter (local maxima of the z MIP, or of
    the volume with ``'3d'``) above their windows' median by ``sigma_thrsh`` sigmas, each
    segmented around its peak (a local Otsu threshold raised until the region has at most
    ``det_vol_thrsh`` voxels), large regions dropped, and the rest filled slice by slice by
    linear interpolation of the pixels around them. Returns a new volume."""
    del otsu_max_per, debug
    v = np.array(vol, copy=True)
    if v.ndim == 2:
        v = v[np.newaxis]
    out = _point_detection(
        v,
        float(sigma),
        _size3(local_win_size, "local_win_size"),
        float(np.asarray(sigma_thrsh).ravel()[0]),
        float(np.asarray(int_thrsh).ravel()[0]),
        str(init_detect),
        float(np.asarray(det_vol_thrsh).ravel()[0]),
        num_threads,
    )
    return out.reshape(np.shape(vol))


# ---------------------------------------------------------------------------------------------
# block tasks
# ---------------------------------------------------------------------------------------------


def _check_flag(flag_fullname: str, inds: list[int], overwrite: bool) -> bool:
    """MATLAB's checks of a block task's flag: True to go on"""
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise PunctaRemovalError(
            f"The block directory {flag_path} does not exist, skip the processing!"
        )
    if os.path.isfile(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            print(
                f"The batch files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!",
                flush=True,
            )
            return False
    return True


def _boxes(v: Any) -> np.ndarray:
    return np.asarray(v, dtype=np.int64).reshape(-1, 6)


def _region_shape(ob: Sequence[int]) -> tuple[int, int, int]:
    return (ob[3] - ob[0], ob[4] - ob[1], ob[5] - ob[2])


def _read_batch(path: str, box: Sequence[int], num_threads: int | None) -> np.ndarray:
    from petakit5d.io import read_zarr

    return _vol3(read_zarr(path, input_bbox=list(box), num_threads=num_threads))


def _write_region(a: np.ndarray, path: str, box: Sequence[int], num_threads: int | None) -> None:
    from petakit5d.io import write_zarr

    write_zarr(a, path, bbox=list(box), create=False, num_threads=num_threads)


def _out_dtype(path: str) -> np.dtype:
    from petakit5d.core.arrays import numpy_dtype
    from petakit5d.io import get_zarr_info

    return numpy_dtype(get_zarr_info(path).dtype)


def _batch_threshold(
    im_g: np.ndarray,
    sigma: float,
    otsu_max_per: float,
    min_threshold: float,
    max_threshold: Any,
    base_threshold: Any,
    num_threads: int | None,
) -> Any:
    """the threshold of a batch's Gaussian filter (the updated method and the background
    estimation); None where MATLAB writes an empty batch (no values for Otsu)"""
    from petakit5d.matlabcompat import morphology as morph
    from petakit5d.matlabcompat import stats
    from petakit5d.matlabcompat.filters import imgaussfilt

    max_bound = math.inf if _empty(max_threshold) else float(np.asarray(max_threshold).ravel()[0])
    base_bound = (
        -math.inf if _empty(base_threshold) else float(np.asarray(base_threshold).ravel()[0])
    )

    mip = im_g.max(axis=0) if im_g.dtype.kind != "f" or not np.isnan(im_g).any() else None
    if mip is None:
        from petakit5d.matlabcompat import reduce as mred

        mip = mred.max(im_g, axis=-3)[0]
    pos = mip[mip > 0]
    p_all = _scalar(stats.prctile(mip, otsu_max_per))
    if im_g.dtype.kind in "ui" and im_g.dtype.itemsize <= 2:
        # (the values in (0, p): the same multiset as im_g(im_g > 0 & im_g < p))
        vals, cnt = _value_counts(im_g, num_threads)
        keep = (vals > 0) & (vals < np.float64(p_all))
        if not keep.any():
            keep = (vals > 0) & (vals < np.float64(_scalar(stats.prctile(pos, otsu_max_per))))
        if not keep.any():
            return None
        t = _otsu_from_counts(*_otsu_counts_values(vals[keep], cnt[keep]))
    else:
        sel = (im_g > 0) & (im_g < p_all)
        if not sel.any():
            sel = (im_g > 0) & (im_g < _scalar(stats.prctile(pos, otsu_max_per)))
        if not sel.any():
            return None
        t = threshold_otsu(im_g[sel])
        del sel
    if _empty(base_threshold):
        if t < min_threshold:
            t = min_threshold
        return t

    # the threshold from the MIP's background mask (used with a base threshold only)
    t_min = _scalar(stats.prctile(pos, 1))
    t_med = _scalar(stats.median(pos))
    mip_bg = imgaussfilt(_min_array(mip, max_bound), sigma * 20, num_threads=num_threads)
    mask = mip > mip_bg
    beta = 1.1
    while np.count_nonzero(mask) / mask.size > 0.9:
        mask = mip > _times(mip_bg, beta)
        beta = beta * 1.1
    disk = morph.strel("disk", 5)
    lo = _mmax(base_bound, t_min)
    hi = _mmax(t_med, max_bound)
    mask = ~morph.imdilate(mask, disk, threads=num_threads) & (
        _gt_scalar(mip, lo) & _lt_scalar(mip, hi)
    )
    mask = morph.imopen(mask, disk, threads=num_threads)
    mask = morph.imerode(mask, disk, threads=num_threads)
    vals = mip[mask]
    t_mask = _scalar(stats.mode(vals))
    if t_mask < _scalar(stats.prctile(vals, 25)):
        t_mask = _scalar(stats.median(vals))
    t_mask = _mmax(base_bound, _mmin(t_mask, max_bound))
    if t > max_bound:
        t = float(np.asarray(max_threshold).ravel()[0])
    base = float(np.asarray(base_threshold).ravel()[0])
    return _scalar(_plus(_times(0.9, _concat_max(base, t_mask)), 0.1 * t))


def _lt_scalar(a: np.ndarray, t: Any) -> np.ndarray:
    t = _scalar(t)
    if np.asarray(t).dtype == a.dtype:
        return a < t
    return a < np.float64(t)


def _masked_single(batch: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """``single(batch) .* keep``"""
    out = batch.astype(np.float32)
    out *= keep
    return out


def _finish_batch(
    out: np.ndarray,
    nz_mask: np.ndarray | None,
    offset: float,
    dtype: np.dtype,
    local_box: Sequence[int],
) -> np.ndarray:
    """offset on the input's nonzero voxels, the output class, the batch's part of the region"""
    from petakit5d.core.arrays import matlab_class
    from petakit5d.core.dtypes import cast

    if offset != 0:
        for z in range(out.shape[0]):
            out[z] = (out[z].astype(np.float64) + nz_mask[z] * float(offset)).astype(np.float32)
    lb = local_box
    part = out[lb[0] : lb[3], lb[1] : lb[4], lb[2] : lb[5]]
    if out.dtype != dtype:
        return np.ascontiguousarray(cast(part, matlab_class(np.empty(0, dtype))))
    return np.ascontiguousarray(part)


def _remove_points(
    out: np.ndarray,
    sigma: float,
    local_win_size: Sequence[int],
    sigma_thrsh: Any,
    int_thrsh: Any,
    init_detect: str,
    det_vol_thrsh: Any,
    num_threads: int | None,
) -> np.ndarray:
    if not out.any():
        return out
    st = _nums(sigma_thrsh)
    it = _nums(int_thrsh)
    dv = _nums(det_vol_thrsh)
    for j in range(len(st)):
        out = _point_detection(
            out, sigma, local_win_size, st[j], it[j], init_detect, dv[j], num_threads
        )
    return out


_BLOCK_COMMON = [
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("Sigma", "@isnumeric", "Sigma of the 3D Gaussian filter."),
    Param("OTSUMaxPer", "@isnumeric", "Percentile of the MIP below which values enter Otsu."),
]


@register(
    "XR_ExM_PunctaRemoval_updated_block",
    Param("batchInds", "@isnumeric", "Batch numbers.", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("outFullpath", "@(x) ischar(x)", "Result zarr."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file of the task."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch in the input.",
        convert="bboxes",
        py="batch_bboxes",
    ),
    Param(
        "RegionBBoxes",
        "@isnumeric",
        "Box of each batch in the result.",
        convert="bboxes",
        py="region_bboxes",
    ),
    Param(
        "localBBoxes",
        "@isnumeric",
        "Box of each result region in its batch.",
        convert="bboxes",
        py="local_bboxes",
    ),
    *_BLOCK_COMMON,
    Param("MinThreshold", "@isnumeric", "Smallest threshold."),
    Param("MaxThreshold", "@isnumeric", "Largest threshold (empty: none)."),
    Param(
        "BaseThreshold", "@isnumeric", "Base threshold weighed with the local ones (empty: none)."
    ),
    Param("volThrsh", "@isvector", "Objects smaller than this are removed."),
    Param("offset", "@isnumeric", "Added to the input's nonzero voxels.", convert="same"),
    Param(
        "localWinSize", "@isvector", "Window (z, y, x) of the point detection.", convert="size_yxz"
    ),
    Param("SigmaThrsh", "@isnumeric", "Sigma thresholds of the point detection (one pass each)."),
    Param("intThrsh", "@isnumeric", "Intensity thresholds of the point detection."),
    Param("initDetect", "@ischar", "Initial detection: '2d' or '3d'."),
    Param("detVolThrsh", "@isnumeric", "Largest punctum volumes of the point detection."),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all).", python_only=True),
    task="B08",
)
def puncta_removal_updated_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    out_fullpath: str,
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    region_bboxes: Sequence[Sequence[int]],
    local_bboxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    sigma: float = 2.5,
    otsu_max_per: float = 99.9,
    min_threshold: float = 0,
    max_threshold: float | None = None,
    base_threshold: float | None = None,
    vol_thrsh: float = 1000,
    offset: float = 0,
    local_win_size: Sequence[int] = (23, 15, 15),
    sigma_thrsh: Any = 4,
    int_thrsh: Any = 10000,
    init_detect: str = "2d",
    det_vol_thrsh: Any = 5000,
    debug: bool = False,
    num_threads: int | None = None,
) -> list[bool] | None:
    """``XR_ExM_PunctaRemoval_updated_block``: each batch (``batch_bboxes`` in the input zarr)
    thresholded on its Gaussian filter (Otsu's threshold of the values below the MIP's
    ``otsu_max_per`` percentile, bounded by ``min_threshold``, or weighed with
    ``base_threshold`` and the MIP background mask's threshold), small objects removed,
    bright puncta removed by point detection, ``offset`` added, cast to the result's class and
    written into its region (``region_bboxes``, the batch's part ``local_bboxes``); then the
    flag. Nothing is done if the flag exists (unless ``overwrite``)."""
    from petakit5d.matlabcompat.filters import imgaussfilt3

    del debug
    inds = _ints(batch_inds)
    if not _check_flag(flag_fullname, inds, overwrite):
        return None
    if not os.path.isdir(zarr_fullpath):
        raise PunctaRemovalError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    if not os.path.isdir(out_fullpath):
        raise PunctaRemovalError(f"The output zarr file {out_fullpath} doesnot exist!")
    dtype = _out_dtype(out_fullpath)
    ib_all, ob_all, lb_all = _boxes(batch_bboxes), _boxes(region_bboxes), _boxes(local_bboxes)
    win = _size3(local_win_size, "local_win_size")
    done = [False] * len(inds)
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        ib, ob, lb = ([int(v) for v in r[i]] for r in (ib_all, ob_all, lb_all))
        batch = _read_batch(zarr_fullpath, ib, num_threads)
        nz_mask = batch != 0 if offset != 0 else None
        if not batch.any():
            _write_region(np.zeros(_region_shape(ob), dtype), out_fullpath, ob, num_threads)
            done[i] = True
            continue
        im_g = imgaussfilt3(batch, sigma, num_threads=num_threads)
        t = _batch_threshold(
            im_g, sigma, otsu_max_per, min_threshold, max_threshold, base_threshold, num_threads
        )
        if t is None:
            _write_region(np.zeros(_region_shape(ob), dtype), out_fullpath, ob, num_threads)
            done[i] = True
            continue
        keep = _area_open(_gt_scalar(im_g, t), float(vol_thrsh), num_threads)
        del im_g
        out = _masked_single(batch, keep)
        del keep, batch
        out = _remove_points(
            out, sigma, win, sigma_thrsh, int_thrsh, init_detect, det_vol_thrsh, num_threads
        )
        res = _finish_batch(out, nz_mask, offset, dtype, lb)
        del out
        _write_region(res, out_fullpath, ob, num_threads)
        done[i] = True
        _toc(t0)
    if all(done):
        open(flag_fullname, "w").close()
    return done


@register(
    "XR_ExM_PunctaRemoval_interpolated_background_block",
    Param("batchInds", "@isnumeric", "Batch numbers.", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("outFullpath", "@(x) ischar(x)", "Result zarr."),
    Param("bgFullpath", "@(x) ischar(x)", "Background zarr of the background estimation."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file of the task."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch in the input.",
        convert="bboxes",
        py="batch_bboxes",
    ),
    Param(
        "RegionBBoxes",
        "@isnumeric",
        "Box of each batch in the result.",
        convert="bboxes",
        py="region_bboxes",
    ),
    Param(
        "localBBoxes",
        "@isnumeric",
        "Box of each result region in its batch.",
        convert="bboxes",
        py="local_bboxes",
    ),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("bgFactor", "@isnumeric", "Factor of the background."),
    Param("Sigma", "@isnumeric", "Sigma of the 3D Gaussian filter."),
    Param("OTSUMaxPer", "@isnumeric", "Passed to the point detection (not used there)."),
    Param("volThrsh", "@isvector", "Objects smaller than this are removed."),
    Param("offset", "@isnumeric", "Added to the input's nonzero voxels.", convert="same"),
    Param(
        "localWinSize", "@isvector", "Window (z, y, x) of the point detection.", convert="size_yxz"
    ),
    Param("SigmaThrsh", "@isnumeric", "Sigma thresholds of the point detection (one pass each)."),
    Param("intThrsh", "@isnumeric", "Intensity thresholds of the point detection."),
    Param("initDetect", "@ischar", "Initial detection: '2d' or '3d'."),
    Param("detVolThrsh", "@isnumeric", "Largest punctum volumes of the point detection."),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all).", python_only=True),
    task="B08",
)
def puncta_removal_interpolated_background_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    out_fullpath: str,
    bg_fullpath: str,
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    region_bboxes: Sequence[Sequence[int]],
    local_bboxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    bg_factor: float = 1.0,
    sigma: float = 2.5,
    otsu_max_per: float = 99.9,
    vol_thrsh: float = 1000,
    offset: float = 0,
    local_win_size: Sequence[int] = (23, 15, 15),
    sigma_thrsh: Any = 4,
    int_thrsh: Any = 10000,
    init_detect: str = "2d",
    det_vol_thrsh: Any = 5000,
    debug: bool = False,
    num_threads: int | None = None,
) -> list[bool] | None:
    """``XR_ExM_PunctaRemoval_interpolated_background_block``: as the updated method, with the
    threshold of each voxel from the estimated background (``bg_fullpath``, one value per
    background batch, enlarged linearly to the batch) times ``bg_factor``."""
    from petakit5d.matlabcompat.filters import imgaussfilt3

    del debug, otsu_max_per
    inds = _ints(batch_inds)
    if not _check_flag(flag_fullname, inds, overwrite):
        return None
    if not os.path.isdir(zarr_fullpath):
        raise PunctaRemovalError(f"The input zarr file {zarr_fullpath} does not exist!")
    if not os.path.isdir(bg_fullpath):
        raise PunctaRemovalError(f"The background zarr file {bg_fullpath} does not exist!")
    if not os.path.isdir(out_fullpath):
        raise PunctaRemovalError(f"The output zarr file {out_fullpath} does not exist!")
    from petakit5d.io import get_image_size

    dtype = _out_dtype(out_fullpath)
    tok = re.search(r"batch_(\d+)_(\d+)_(\d+)_background_values", bg_fullpath)
    if tok is None:
        raise PunctaRemovalError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    # MATLAB's name holds the batch size in (y, x, z)
    bgb = [int(tok.group(3)), int(tok.group(1)), int(tok.group(2))]
    bg_size = [int(v) for v in get_image_size(bg_fullpath)[:3]]
    ib_all, ob_all, lb_all = _boxes(batch_bboxes), _boxes(region_bboxes), _boxes(local_bboxes)
    win = _size3(local_win_size, "local_win_size")
    done = [False] * len(inds)
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        ib, ob, lb = ([int(v) for v in r[i]] for r in (ib_all, ob_all, lb_all))
        batch = _read_batch(zarr_fullpath, ib, num_threads)
        nz_mask = batch != 0 if offset != 0 else None
        if not batch.any():
            _write_region(np.zeros(_region_shape(ob), dtype), out_fullpath, ob, num_threads)
            done[i] = True
            continue
        im_g = imgaussfilt3(batch, sigma, num_threads=num_threads)
        n = [ib[q + 3] - ib[q] for q in range(3)]
        b0 = [ib[q] // bgb[q] for q in range(3)]
        b1 = [min(bg_size[q], b0[q] + math.ceil(n[q] / bgb[q])) for q in range(3)]
        bg = _read_batch(bg_fullpath, b0 + b1, num_threads)
        if bg.shape[0] == 1:
            # one background value deep in z (MATLAB's array is 2-D): the plane's
            # interpolation at every z
            im_bg = np.broadcast_to(imresize_bilinear(bg[0], n[1:]), n)
        else:
            im_bg = imresize3_linear(bg, n)
        thr = _times(im_bg, float(bg_factor))
        del im_bg
        above = np.empty(im_g.shape, dtype=bool)
        for z in range(im_g.shape[0]):
            above[z] = im_g[z].astype(np.float64) > thr[z].astype(np.float64)
        del im_g, thr
        keep = _area_open(above, float(vol_thrsh), num_threads)
        del above
        out = _masked_single(batch, keep)
        del keep, batch
        out = _remove_points(
            out, sigma, win, sigma_thrsh, int_thrsh, init_detect, det_vol_thrsh, num_threads
        )
        res = _finish_batch(out, nz_mask, offset, dtype, lb)
        del out
        _write_region(res, out_fullpath, ob, num_threads)
        done[i] = True
        _toc(t0)
    if all(done):
        open(flag_fullname, "w").close()
    return done


@register(
    "XR_ExM_PunctaRemoval_background_estimation_block",
    Param("batchInds", "@isnumeric", "Batch numbers.", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("resultFullname", "@(x) ischar(x)", "Result .mat file of the task."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch in the input.",
        convert="bboxes",
        py="batch_bboxes",
    ),
    *_BLOCK_COMMON,
    Param("MinThreshold", "@isnumeric", "Smallest threshold."),
    Param("MaxThreshold", "@isnumeric", "Largest threshold (empty: none)."),
    Param(
        "BaseThreshold", "@isnumeric", "Base threshold weighed with the local ones (empty: none)."
    ),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all).", python_only=True),
    task="B08",
)
def puncta_removal_background_estimation_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    result_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    sigma: float = 2.5,
    otsu_max_per: float = 99.9,
    min_threshold: float = 0,
    max_threshold: float | None = None,
    base_threshold: float | None = None,
    uuid: str = "",
    debug: bool = False,
    num_threads: int | None = None,
) -> list[bool] | None:
    """``XR_ExM_PunctaRemoval_background_estimation_block``: the threshold of each batch (as
    the updated method finds it; 0 for an empty batch) saved as ``bg_mat`` (a double column,
    with ``batchInds`` 1-based) in ``result_fullname`` (written to a temporary name, then
    renamed). Nothing is done if the result exists (unless ``overwrite``)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.io import write_mat
    from petakit5d.matlabcompat.filters import imgaussfilt3

    del debug
    if not uuid:
        uuid = get_uuid()
    inds = _ints(batch_inds)
    if os.path.isfile(result_fullname):
        if overwrite:
            os.remove(result_fullname)
        else:
            print(
                f"The batch files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!",
                flush=True,
            )
            return None
    if not os.path.isdir(zarr_fullpath):
        raise PunctaRemovalError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    ib_all = _boxes(batch_bboxes)
    done = [False] * len(inds)
    bg = np.zeros(len(inds), dtype=np.float64)
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        ib = [int(v) for v in ib_all[i]]
        batch = _read_batch(zarr_fullpath, ib, num_threads)
        if not batch.any():
            done[i] = True
            continue
        im_g = imgaussfilt3(batch, sigma, num_threads=num_threads)
        del batch
        t = _batch_threshold(
            im_g, sigma, otsu_max_per, min_threshold, max_threshold, base_threshold, num_threads
        )
        bg[i] = 0.0 if t is None else float(t)
        done[i] = True
        _toc(t0)
    if all(done):
        path, fsn, ext = fileparts(result_fullname)
        tmp = f"{path}/{fsn}_{uuid}{ext}"
        write_mat(
            tmp,
            {
                "bg_mat": bg.reshape(-1, 1),
                "batchInds": (np.asarray(inds, dtype=np.float64) + 1).reshape(1, -1),
            },
        )
        os.replace(tmp, result_fullname)
    return done


@register(
    "XR_ExM_PunctaRemoval_block",
    Param("blockInds", "@isnumeric", "Batch numbers.", convert="index1", py="batch_inds"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param("outFullpath", "@(x) ischar(x)", "Result zarr."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file of the task."),
    Param(
        "BatchBBoxes",
        "@isnumeric",
        "Box of each batch in the input.",
        convert="bboxes",
        py="batch_bboxes",
    ),
    Param(
        "RegionBBoxes",
        "@isnumeric",
        "Box of each batch in the result.",
        convert="bboxes",
        py="region_bboxes",
    ),
    Param(
        "localBBoxes",
        "@isnumeric",
        "Box of each result region in its batch.",
        convert="bboxes",
        py="local_bboxes",
    ),
    *_BLOCK_COMMON,
    Param("MinThreshold", "@isnumeric", "Smallest threshold."),
    Param("MaxThreshold", "@isnumeric", "Largest threshold (empty: none)."),
    Param(
        "BaseThreshold", "@isnumeric", "Base threshold weighed with the local one (empty: none)."
    ),
    Param("volThrsh", "@isvector", "Objects smaller than this are removed."),
    Param("uuid", "@ischar", "Not used."),
    Param("debug", "@islogical", "Not used."),
    Param("num_threads", None, "Threads (None: all).", python_only=True),
    task="B08",
)
def puncta_removal_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    out_fullpath: str,
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    region_bboxes: Sequence[Sequence[int]],
    local_bboxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    sigma: float = 2.5,
    otsu_max_per: float = 99.9,
    min_threshold: float = 0,
    max_threshold: float | None = None,
    base_threshold: float | None = None,
    vol_thrsh: float = 1000,
    uuid: str = "",
    debug: bool = False,
    num_threads: int | None = None,
) -> list[bool] | None:
    """``XR_ExM_PunctaRemoval_block`` (the first method; the zarr workflow doesn't use it): each
    batch thresholded on its Gaussian filter at Otsu's threshold of its values below their
    ``otsu_max_per`` percentile (bounded by ``min_threshold``, or ``0.75 * base_threshold +
    0.25 * T``), small objects removed, written into its region; then the flag."""
    from petakit5d.matlabcompat import label, stats
    from petakit5d.matlabcompat.filters import imgaussfilt3

    del uuid, debug
    inds = _ints(batch_inds)
    if not _check_flag(flag_fullname, inds, overwrite):
        return None
    if not os.path.isdir(zarr_fullpath):
        raise PunctaRemovalError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    if not os.path.isdir(out_fullpath):
        raise PunctaRemovalError(f"The output zarr file {out_fullpath} doesnot exist!")
    dtype = _out_dtype(out_fullpath)
    ib_all, ob_all, lb_all = _boxes(batch_bboxes), _boxes(region_bboxes), _boxes(local_bboxes)
    done = [False] * len(inds)
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        ib, ob, lb = ([int(v) for v in r[i]] for r in (ib_all, ob_all, lb_all))
        batch = _read_batch(zarr_fullpath, ib, num_threads)
        if not batch.any():
            _write_region(np.zeros(_region_shape(ob), dtype), out_fullpath, ob, num_threads)
            done[i] = True
            continue
        im_g = imgaussfilt3(batch, sigma, num_threads=num_threads)
        p = _scalar(stats.prctile(batch, otsu_max_per))
        t = threshold_otsu(batch[(batch > 0) & (batch < p)])
        if not _empty(base_threshold):
            if not _empty(max_threshold) and t > float(np.asarray(max_threshold).ravel()[0]):
                t = float(np.asarray(max_threshold).ravel()[0])
            t = 0.75 * float(np.asarray(base_threshold).ravel()[0]) + 0.25 * t
        elif t < min_threshold:
            t = min_threshold
        keep = label.bwareaopen(_gt_scalar(im_g, t), vol_thrsh, 26, threads=num_threads)
        del im_g
        out = _masked_single(batch, keep)
        del keep, batch
        res = _finish_batch(out, None, 0, dtype, lb)
        _write_region(res, out_fullpath, ob, num_threads)
        done[i] = True
        _toc(t0)
    if all(done):
        open(flag_fullname, "w").close()
    return done


# ---------------------------------------------------------------------------------------------
# the zarr workflow
# ---------------------------------------------------------------------------------------------


def _yxz_name(v: Sequence[int]) -> str:
    """a (z, y, x) size in MATLAB's file names: ``<y>_<x>_<z>``"""
    return f"{v[1]}_{v[2]}_{v[0]}"


def _save_parameters(out_path: str, values: dict[str, Any]) -> None:
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(puncta_removal_zarr, values).items()}
    write_mat(out_path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, out_path + "/parameters.json")


@register(
    "XR_ExM_PunctaRemoval_zarr",
    Param("zarrFullpath", "@(x) ischar(x) || iscell(x)", "Input zarr."),
    Param(
        "mode",
        "@(x) ischar(x) || iscell(x)",
        "'updated', 'bg_interpolated' (thresholds from bg_fullpath) or 'bg_estimation'.",
    ),
    Param(
        "bgFullpath", "@(x) ischar(x)", "Background zarr of 'bg_estimation' ('bg_interpolated')."
    ),
    Param("bgFactor", "@(x) isnumeric(x)", "Factor of the background ('bg_interpolated')."),
    Param("Sigma", "@isnumeric", "Sigma of the 3D Gaussian filter."),
    Param("OTSUMaxPer", "@isnumeric", "Percentile of the MIP below which values enter Otsu."),
    Param("MinThreshold", "@isnumeric", "Smallest threshold."),
    Param("MaxThreshold", "@isnumeric", "Largest threshold (empty: none)."),
    Param(
        "BaseThreshold", "@isnumeric", "Base threshold weighed with the local ones (empty: none)."
    ),
    Param("volThrsh", "@isnumeric", "Objects smaller than this are removed."),
    Param("offset", "@isnumeric", "Added to the input's nonzero voxels.", convert="same"),
    Param(
        "localWinSize", "@isvector", "Window (z, y, x) of the point detection.", convert="size_yxz"
    ),
    Param("SigmaThrsh", "@isnumeric", "Sigma thresholds of the point detection (one pass each)."),
    Param("intThrsh", "@isnumeric", "Intensity thresholds of the point detection."),
    Param("initDetect", "@ischar", "Initial detection: '2d' (z MIP) or '3d' (the volume)."),
    Param("detVolThrsh", "@isnumeric", "Largest punctum volumes of the point detection."),
    Param(
        "bbox",
        "@(x) isempty(x) || isvector(x)",
        "Box (z0, y0, x0, z1, y1, x1) to process.",
        convert="bbox",
    ),
    Param("batchSize", "@isvector", "Batch size (z, y, x).", convert="size_yxz"),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of the result.", convert="size_yxz"),
    Param("Save16bit", "@islogical", "Save the result as uint16 (else single)."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param("parseCluster", "@islogical", "Run the batches as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes batches too."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="B08",
)
def puncta_removal_zarr(
    zarr_fullpath: str,
    *,
    mode: str = "updated",
    bg_fullpath: str = "",
    bg_factor: float = 1.0,
    sigma: float = 2.5,
    otsu_max_per: float = 99.9,
    min_threshold: float = 0,
    max_threshold: float | None = None,
    base_threshold: float | None = None,
    vol_thrsh: float = 1000,
    offset: float = 0,
    local_win_size: Sequence[int] = (23, 15, 15),
    sigma_thrsh: Any = 4,
    int_thrsh: Any = 10000,
    init_detect: str = "2d",
    det_vol_thrsh: Any = 5000,
    bbox: Sequence[int] | None = None,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    save_16bit: bool = False,
    uuid: str = "",
    debug: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 1,
    mcc_mode: bool = False,
    config_file: str = "",
) -> list[bool] | None:
    """Puncta removal of a large zarr file in batches (``XR_ExM_PunctaRemoval_zarr``): the
    batches (``batch_size`` with borders of 100, whole chunks of ``block_size``) as block tasks
    of ``mode`` through the job framework into ``<data>/Puncta_Removed/<name>.zarr_<uuid>``,
    its MIPs, then the result moved into place; with ``'bg_estimation'`` the thresholds of the
    batches into ``<data>/Puncta_background_estimation/<name>_batch_<y>_<x>_<z>
    _background_values.zarr`` instead. Nothing is done if the result exists. Returns the tasks'
    done flags (None if nothing was done)."""
    values = dict(locals())
    from petakit5d.core.blocks import zarr_chunk_coordinates
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive, strip_trailing_sep
    from petakit5d.core.text import through_mat2str, through_text
    from petakit5d.io import create_zarr, get_image_size, get_zarr_info, read_mat, write_zarr
    from petakit5d.jobs import framework
    from petakit5d.jobs.publish import publish_result, remove_path
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.tools.mip import mip_zarr

    # (parameters.mat/.json hold the call's uuid, empty when it is made here, as in MATLAB)
    if not uuid:
        uuid = get_uuid()
    bg_estimate = mode == "bg_estimation"
    zarr_fullpath = os.fspath(zarr_fullpath)
    if not os.path.isdir(zarr_fullpath):
        raise PunctaRemovalError(f"The input zarr file {zarr_fullpath} does not exist!")
    batch = _size3(batch_size, "batch_size")
    block = _size3(block_size, "block_size")
    data_path, fsname, _ext = fileparts(strip_trailing_sep(zarr_fullpath))
    if bg_estimate:
        out_rt = f"{data_path}/Puncta_background_estimation/"
        os.makedirs(out_rt, exist_ok=True)
        out_path = f"{out_rt}/{fsname}_batch_{_yxz_name(batch)}/"
        out_fullpath = f"{out_rt}/{fsname}_batch_{_yxz_name(batch)}_background_values.zarr"
    else:
        out_path = f"{data_path}/Puncta_Removed/"
        out_fullpath = f"{out_path}/{fsname}.zarr"
        print(f"Start Large-file puncta removal for {fsname}...", flush=True)
    if os.path.isdir(out_fullpath):
        print("The output result exists, skip it!", flush=True)
        return None
    tmp_path = f"{out_fullpath}_{uuid}"
    if not os.path.isdir(out_path):
        os.makedirs(out_path, exist_ok=True)
    _save_parameters(out_path, values)
    print(f"reading {fsname}...", flush=True)
    dtype = "uint16" if save_16bit else "single"
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(data_path).parse_cluster)
    t_start = time.perf_counter()

    im_size = [int(v) for v in get_image_size(zarr_fullpath)[:3]]
    border = [100, 100, 100]
    if bg_estimate:
        from petakit5d.core.dtypes import mround

        border = [min(int(mround(b / 2)), c) for b, c in zip(batch, border, strict=True)]
    batch_boxes, region_boxes, local_boxes = zarr_chunk_coordinates(
        im_size,
        batch_size=batch,
        block_size=block,
        same_batch_size=False,
        border_size=border,
        bbox=None if _empty(bbox) else list(bbox),
    )
    if _empty(bbox):
        out_size = im_size
    else:
        bb = _ints(bbox)
        out_size = [bb[q + 3] - bb[q] for q in range(3)]

    fresh = not bg_estimate
    if not bg_estimate and os.path.isdir(tmp_path):
        info = get_zarr_info(tmp_path)
        if list(info.size[:3]) == out_size and list(info.block_size[:3]) == block:
            fresh = False
        else:
            remove_path(tmp_path)
    if fresh:
        n_chunks = 1
        for o, b in zip(out_size, block, strict=True):
            n_chunks *= math.ceil(o / b)
        create_zarr(
            tmp_path,
            data_size=out_size,
            block_size=block,
            dtype=dtype,
            dim_separator="/" if n_chunks > 10000 else ".",
        )
    flag_path = ""
    if not bg_estimate:
        flag_path = f"{out_path}/zarr_flag/{fsname}_{uuid}/"
        if fresh and os.path.isdir(flag_path):
            remove_path(flag_path)
        mkdir_recursive(flag_path)

    task_size = 5 * math.ceil(math.prod(1024 / b for b in batch))
    num_batch = int(batch_boxes.shape[0])
    num_tasks = math.ceil(num_batch / task_size)
    if num_tasks > 20000:
        task_size = min(100000, math.ceil(num_batch / 20000))
        num_tasks = math.ceil(num_batch / task_size)

    # the values as the task strings pass them (%0.20d exact, %d of a fraction as %e, mat2str)
    off = float(through_text(offset, "%d"))
    mt = through_mat2str(min_threshold)
    mx = None if _empty(max_threshold) else through_mat2str(max_threshold)
    bs = None if _empty(base_threshold) else through_mat2str(base_threshold)
    win = [int(v) for v in through_mat2str(np.asarray(_size3(local_win_size, "local_win_size")))]
    st = through_mat2str(np.asarray(sigma_thrsh, dtype=np.float64))
    it = through_mat2str(np.asarray(int_thrsh, dtype=np.float64))
    dv = through_mat2str(np.asarray(det_vol_thrsh, dtype=np.float64))
    tasks, outputs = [], []
    for t in range(num_tasks):
        inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
        bb_i = batch_boxes[inds].tolist()
        if bg_estimate:
            out_t = f"{out_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
            task = Task.of(
                puncta_removal_background_estimation_block,
                inds,
                zarr_fullpath,
                out_t,
                bb_i,
                overwrite=False,
                sigma=float(sigma),
                otsu_max_per=float(otsu_max_per),
                min_threshold=mt,
                max_threshold=mx,
                base_threshold=bs,
                uuid=uuid,
                debug=bool(debug),
            )
        else:
            out_t = f"{flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
            common = {
                "overwrite": False,
                "sigma": float(sigma),
                "otsu_max_per": float(otsu_max_per),
                "vol_thrsh": float(vol_thrsh),
                "offset": off,
                "local_win_size": win,
                "sigma_thrsh": st,
                "int_thrsh": it,
                "init_detect": init_detect,
                "det_vol_thrsh": dv,
                "debug": bool(debug),
            }
            rb_i, lb_i = region_boxes[inds].tolist(), local_boxes[inds].tolist()
            if mode == "bg_interpolated":
                task = Task.of(
                    puncta_removal_interpolated_background_block,
                    inds,
                    zarr_fullpath,
                    tmp_path,
                    bg_fullpath,
                    out_t,
                    bb_i,
                    rb_i,
                    lb_i,
                    bg_factor=float(bg_factor),
                    **common,
                )
            elif mode == "updated":
                task = Task.of(
                    puncta_removal_updated_block,
                    inds,
                    zarr_fullpath,
                    tmp_path,
                    out_t,
                    bb_i,
                    rb_i,
                    lb_i,
                    min_threshold=mt,
                    max_threshold=mx,
                    base_threshold=bs,
                    **common,
                )
            else:
                # (MATLAB leaves the task string empty: the framework fails the task)
                raise PunctaRemovalError(f"Unknown puncta removal mode {mode!r}")
        tasks.append(task)
        outputs.append(out_t)

    inputs = [zarr_fullpath] * num_tasks
    mem = math.prod(b + c * 2 for b, c in zip(batch, border, strict=True)) * 4 / 1024**3 * 5
    done = framework.run_tasks(
        inputs,
        outputs,
        tasks,
        final_out_fullpath=out_fullpath,
        parse_cluster=parse_cluster,
        master_compute=master_compute,
        task_batch_num=1,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )
    if not all(done):
        raise PunctaRemovalError("The computing for some batchs is missing!")

    if bg_estimate:
        print("Collect local background values...", flush=True)
        bg_size = [math.ceil(o / b) for o, b in zip(out_size, batch, strict=True)]
        bg = np.zeros(bg_size, dtype=np.float32)
        for t in range(num_tasks):
            inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
            a = read_mat(f"{out_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat")
            vals = np.asarray(a["bg_mat"], dtype=np.float64).reshape(-1).astype(np.float32)
            rb = region_boxes[inds]
            for k in range(len(inds)):
                g = tuple(int(rb[k][q]) // batch[q] for q in range(3))
                bg[g] = vals[k]
        bs_bg = [math.ceil(1024 / b) for b in batch]
        create_zarr(
            tmp_path, data_size=bg_size, block_size=bs_bg, dtype="single", dim_separator="."
        )
        write_zarr(bg, tmp_path, create=False)
        print("Done!", flush=True)
        publish_result(tmp_path, out_fullpath)
    else:
        # the MIPs before the result is moved into place (the framework counts the task as
        # done once it exists); MATLAB makes them from the result after the move
        if os.path.isdir(tmp_path) and not os.path.isdir(out_fullpath):
            mip_zarr(
                tmp_path, parse_cluster=parse_cluster, mcc_mode=mcc_mode, config_file=config_file
            )
        publish_result(tmp_path, out_fullpath, flag_path)
    print(f"Elapsed time is {time.perf_counter() - t_start:.6f} seconds.", flush=True)
    return done
