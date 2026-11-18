"""The stitching registration's pieces of one tile pair in memory: the 2D normalized
cross-correlation, the crop of a long template to its brightest part and the choice of the
template's slices.

==================================  ===========================================================
Python                              MATLAB
==================================  ===========================================================
``normxcorr2``                      normxcorr2 (Image Processing Toolbox): a new implementation
                                    from its documentation and J. P. Lewis, "Fast Normalized
                                    Cross-Correlation" (Vision Interface 1995)
``normxcorr2_max_shift``            normxcorr2_max_shift (on ``kernels.xcorr``'s window and
                                    first maximum, with ``normxcorr2``)
``normxcorr3_max_shift``            normxcorr3_max_shift (``kernels.xcorr``, the 3D NCC)
``crop_subregion_by_intensity``     crop_subregion_by_intensity
``template_indices``                cross_correlation_registration_3d/_2d: the template's rows,
                                    columns and slices with enough data
``cap_region_inds``                 cross_correlation_registration_3d's cap_region_2_inds
==================================  ===========================================================

Arrays are the port's: volumes ``(z, y, x)``, images ``(y, x)``, C order; boxes
``(z0, y0, x0, z1, y1, x1)``, 0-based and end-exclusive. Indices along an axis are 0-based.

``normxcorr2`` follows Lewis's method, which normxcorr2's documentation cites: the correlation of
the image with the template less its mean by FFTs, the windows' sums and sums of squares from
running sums, then the correlation coefficient. Everything is in double; the template and the image
are each taken less their means (with the image's zero padding accounted for), so offsets in the
data don't cost precision, and each window's sum of squares is added from its image part and its
padding (terms of one sign). Windows whose values are all equal (the padding's zeros included) give
0. The coefficients are not clipped to [-1, 1].

MATLAB's own function chooses its FFTs and sums by itself (its results change in the last bits
with its thread count), so the port's coefficients are not MATLAB's bits: they agree to about
1e-13 where the data have contrast and are closer to the exact NCC where they don't (MATLAB's error
reaches 1e-7 on 16-bit data at an offset of 2^15 with a contrast of 12, and there it sets its true
peak of 1 to 0); ``docs/parity/stitch.md`` (proposed). Nearly flat windows (variance far below
the data's rounding) give MATLAB rounding noise and the port their NCC. The input checks and errors
are MATLAB's (from probes of the function).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

__all__ = [
    "RegistrationError",
    "cap_region_inds",
    "crop_subregion_by_intensity",
    "normxcorr2",
    "normxcorr2_max_shift",
    "normxcorr3_max_shift",
    "template_indices",
]

# the template's margins cut by cap_region_2_inds: [0.1, 0.25] of its length
_CAP_BOUNDS = (0.1, 0.25)
_NUMERIC = (
    "logical",
    "double",
    "single",
    "uint8",
    "uint16",
    "uint32",
    "uint64",
    "int8",
    "int16",
    "int32",
    "int64",
)


class RegistrationError(ValueError):
    """An error MATLAB raises; ``identifier`` is its MATLAB error identifier."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _bad_index(n: int) -> RegistrationError:
    return RegistrationError(
        f"Index exceeds the number of array elements. Index must not exceed {n}.",
        "MATLAB:badsubscript",
    )


# ---------------------------------------------------------------------------------------------
# normxcorr2
# ---------------------------------------------------------------------------------------------


def _checked(x, k: int, name: str) -> np.ndarray:
    """validateattributes(x, {numeric classes, logical}, {'real', 'nonsparse', '2d',
    'finite'}) as normxcorr2 checks its inputs"""
    import scipy.sparse as sp

    if sp.issparse(x):
        raise RegistrationError(
            f"Expected input number {k}, {name}, to be nonsparse.",
            "MATLAB:normxcorr2:expectedNonsparse",
        )
    a = np.asarray(x)
    if a.dtype.kind not in "biufc" or a.dtype.itemsize > 16:
        raise RegistrationError(
            f"Expected input number {k}, {name}, to be one of these types:\n\n"
            + ", ".join(_NUMERIC)
            + f"\n\nInstead its type was {type(x).__name__}.",
            "MATLAB:normxcorr2:invalidType",
        )
    if a.dtype.kind == "c":
        raise RegistrationError(
            f"Expected input number {k}, {name}, to be real.", "MATLAB:normxcorr2:expectedReal"
        )
    if a.ndim > 2:
        raise RegistrationError(
            f"Expected input number {k}, {name}, to be two-dimensional.",
            "MATLAB:normxcorr2:expected2D",
        )
    if a.ndim < 2:
        a = a.reshape(1, -1)  # a vector is a row
    if a.dtype.kind == "f" and not np.isfinite(a).all():
        raise RegistrationError(
            f"Expected input number {k}, {name}, to be finite.",
            "MATLAB:normxcorr2:expectedFinite",
        )
    return a


def _fast_length(n: int) -> int:
    """the smallest 2,3,5,7-smooth length >= n (the FFTs' size; any size >= the full output
    gives the same linear correlation)"""
    best = 1 << max(0, (n - 1).bit_length())
    p7 = 1
    while p7 < best:
        p5 = p7
        while p5 < best:
            p3 = p5
            while p3 < best:
                v = p3
                while v < n:
                    v *= 2
                best = min(best, v)
                p3 *= 3
            p5 *= 5
        p7 *= 7
    return best


def _window_ranges(n_win: int, n_img: int) -> tuple[np.ndarray, np.ndarray]:
    """for each full-output position k along an axis (0 .. n_img + n_win - 2), the image's
    first and last index (inclusive) under the window"""
    k = np.arange(n_img + n_win - 1)
    return np.maximum(0, k - n_win + 1), np.minimum(k, n_img - 1)


def _box_sums(x: np.ndarray, lo: Sequence[np.ndarray], hi: Sequence[np.ndarray]) -> np.ndarray:
    """sums of x over the boxes lo[d][k] .. hi[d][k] (inclusive) along each axis d, for every
    combination of the k: running sums along y, then x (in double, in index order)"""
    s = x
    for d in range(2):
        p = np.cumsum(s, axis=d, dtype=np.float64)
        pad = [(0, 0), (0, 0)]
        pad[d] = (1, 0)
        p = np.pad(p, pad)
        s = np.take(p, hi[d] + 1, axis=d) - np.take(p, lo[d], axis=d)
    return s


def _window_sums(x: np.ndarray, m: int, n: int) -> np.ndarray:
    """the sums of the zero-padded image x over every m x n window of the full output: running
    sums along y, then x, in index order (the values of ``_box_sums`` of the image's windows;
    the padding adds exact zeros)"""
    M, N = x.shape
    p = np.zeros((M + 2 * m - 1, N), dtype=np.float64)
    np.cumsum(x, axis=0, dtype=np.float64, out=p[m : m + M])
    p[m + M :] = p[m + M - 1]
    s = p[m:] - p[:-m]
    q = np.zeros((s.shape[0], N + 2 * n - 1), dtype=np.float64)
    np.cumsum(s, axis=1, dtype=np.float64, out=q[:, n : n + N])
    q[:, n + N :] = q[:, n + N - 1 : n + N]
    return q[:, n:] - q[:, :-n]


def _sum(x: np.ndarray) -> float:
    from petakit5d.matlabcompat import reduce

    return float(np.asarray(reduce.sum(x)).reshape(-1)[0])


def _mean(x: np.ndarray) -> float:
    return _sum(x) / x.size


def normxcorr2(template, image, *, threads: int | None = None) -> np.ndarray:
    """``C = normxcorr2(template, A)``: the normalized cross-correlation of the image ``A`` with the
    template at every position where they overlap (``(ny + my - 1, nx + mx - 1)``, double; the
    template's last element on the image's ``(i, j)`` at ``C[i, j]``), the image zero-padded.

    Inputs: 2-D, real, finite arrays of any numeric class or bool (a vector is a row), a template
    of at least 2 elements that are not all the same, an image at least as large along each axis.
    Raises ``RegistrationError`` with MATLAB's identifiers otherwise.
    """
    T = _checked(template, 1, "T")
    A = _checked(image, 2, "A")
    if T.size < 2:
        raise RegistrationError(
            "TEMPLATE must contain at least 2 elements.", "images:normxcorr2:invalidTemplate"
        )
    if A.shape[0] < T.shape[0] or A.shape[1] < T.shape[1]:
        raise RegistrationError(
            "A must be the same size or larger than TEMPLATE.",
            "images:normxcorr2:invalidSizeForA",
        )
    T = T.astype(np.float64)
    A = A.astype(np.float64)
    if np.all(T == T.flat[0]):
        raise RegistrationError(
            "The values of TEMPLATE cannot all be the same.",
            "images:normxcorr2:sameElementsInTemplate",
        )
    from petakit5d import fft as pfft

    m, n = T.shape
    M, N = A.shape
    nt = m * n
    out_shape = (M + m - 1, N + n - 1)
    # the template and the image less their means; the image's padding is -am then (sums in a
    # fixed order: the same bits on every platform)
    T0 = T - _mean(T)
    am = _mean(A)
    A0 = A - am
    # the correlation of A0 (zero-padded) with T0, by FFTs of smooth sizes
    fshape = (_fast_length(out_shape[0]), _fast_length(out_shape[1]))
    hs = pfft.HalfSpectrum(fshape, np.float64, backend="permissive", threads=threads)
    a_pad = np.zeros(fshape)
    a_pad[:M, :N] = A0
    t_pad = np.zeros(fshape)
    t_pad[:m, :n] = T0[::-1, ::-1]
    s = hs.forward(a_pad)
    hs.multiply(s, hs.forward(t_pad))
    corr = hs.inverse(s)[: out_shape[0], : out_shape[1]]
    # window sums of A0 and A0^2 over the image, and of T0 over its part on the image
    lo_a, hi_a = zip(_window_ranges(m, M), _window_ranges(n, N), strict=True)
    a2 = A0 * A0
    s_a = _window_sums(A0, m, n)
    s_a2 = _window_sums(a2, m, n)
    lo_t = [
        lo - k + w - 1
        for lo, k, w in zip(
            lo_a, (np.arange(out_shape[0]), np.arange(out_shape[1])), (m, n), strict=True
        )
    ]
    hi_t = [
        hi - k + w - 1
        for hi, k, w in zip(
            hi_a, (np.arange(out_shape[0]), np.arange(out_shape[1])), (m, n), strict=True
        )
    ]
    s_t = _box_sums(T0, lo_t, hi_t)
    cnt = np.outer(hi_a[0] - lo_a[0] + 1, hi_a[1] - lo_a[1] + 1)
    n_pad = nt - cnt
    # the window's sum of A - am (its padding -am) and its mean
    s_b = s_a - n_pad * am
    mean_w = s_b / nt
    sum_t0 = _sum(T0)
    # the numerator of Lewis's coefficient: sum((a - mean_w) .* t0) over the window
    num = corr - am * (sum_t0 - s_t) - mean_w * sum_t0
    # sum((a - mean_w)^2) over the window: its image part and its padding, each a sum of terms
    # of one sign (no cancellation where the padding dominates)
    var_w = (s_a2 - 2 * mean_w * s_a + cnt * (mean_w * mean_w)) + n_pad * ((am + mean_w) ** 2)
    den_t = _sum(T0 * T0) - sum_t0 * sum_t0 / nt
    # windows whose values are all equal (the padding's zeros included) give 0. A flat window's
    # computed variance is the rounding of the running sums, below about 10 eps times this bound
    # of their magnitude; where any window is below 1e4 eps times it, the flat ones are found
    # exactly (max == min)
    valid = var_w > 0
    bound = (_sum(a2) + nt * (am * am)) * (A.size / nt + 1)
    low = valid & (var_w <= 1e4 * np.finfo(np.float64).eps * bound)
    if low.any():
        valid &= ~_flat_windows(A, m, n)
    C = np.zeros(out_shape)
    den = np.sqrt(np.maximum(var_w, 0.0) * den_t)
    np.divide(num, den, out=C, where=valid)
    return C


def _flat_windows(A: np.ndarray, m: int, n: int) -> np.ndarray:
    """whether every window of the full output holds one value (max == min over the zero-padded
    image)"""
    from scipy import ndimage

    M, N = A.shape
    P = np.zeros((M + 2 * (m - 1), N + 2 * (n - 1)))
    P[m - 1 : m - 1 + M, n - 1 : n - 1 + N] = A
    out = []
    for fn in (ndimage.maximum_filter1d, ndimage.minimum_filter1d):
        v = fn(P, m, axis=0, mode="constant")[m // 2 : m // 2 + M + m - 1]
        v = fn(v, n, axis=1, mode="constant")[:, n // 2 : n // 2 + N + n - 1]
        out.append(v)
    return out[0] == out[1]


def normxcorr2_max_shift(T, A, max_shifts, *, threads: int | None = None):
    """``[max_off, max_corr, C] = normxcorr2_max_shift(T, A, maxShifts)`` for images ``(y, x)``:
    ``normxcorr2`` in the window of the shifts within ``max_shifts`` (``(y, x)`` or ``(z, y, x)``
    values, or rows of lower and upper bounds), its first maximum in MATLAB's order and the
    template's offset there, ``(0, y, x)``."""
    from petakit5d.kernels import xcorr as kx

    def fn(t, a):
        return normxcorr2(t, a, threads=threads)

    return kx.normxcorr2_max_shift(T, A, max_shifts, normxcorr2=fn)


def normxcorr3_max_shift(T, A, max_shifts, fast_fft: bool = False, **kwargs):
    """``[max_off, max_corr, C] = normxcorr3_max_shift(T, A, maxShifts, fastFFT)``
    (``kernels.xcorr.normxcorr3_max_shift``)."""
    from petakit5d.kernels import xcorr as kx

    return kx.normxcorr3_max_shift(T, A, max_shifts, fast_fft, **kwargs)


# ---------------------------------------------------------------------------------------------
# crop_subregion_by_intensity
# ---------------------------------------------------------------------------------------------

# MATLAB dimension (1 y, 2 x, 3 z) -> axis of a (z, y, x) volume
_AXIS = {1: 1, 2: 2, 3: 0}


def _msize(a: np.ndarray) -> list[int]:
    """MATLAB's size(a) of a (z, y, x) volume or (y, x) image (trailing ones dropped)"""
    if a.ndim == 2:
        return [a.shape[0], a.shape[1]]
    s = [a.shape[1], a.shape[2], a.shape[0]]
    return s if s[2] != 1 else s[:2]


def _zero_runs(v: np.ndarray, min_len: float) -> np.ndarray:
    """1-based positions of the runs of zeros of v of at least min_len elements (bwconncomp of
    v == 0 on a vector), in increasing order"""
    z = np.concatenate(([False], v == 0, [False]))
    d = np.diff(z.astype(np.int8))
    starts = np.flatnonzero(d == 1)
    ends = np.flatnonzero(d == -1)
    keep = (ends - starts) >= min_len
    if not keep.any():
        return np.zeros(0, dtype=np.int64)
    return np.concatenate(
        [np.arange(s, e) + 1 for s, e in zip(starts[keep], ends[keep], strict=True)]
    )


def crop_subregion_by_intensity(
    region,
    dim_num_thrsh: float = 2000,
    blank_num_trsh: float | None = None,
    *,
    num_threads: int | None = None,
) -> tuple[np.ndarray, tuple[int, ...]]:
    """``[region, crop_bbox] = crop_subregion_by_intensity(region, dimNumThrsh, blankNumTrsh)``:
    the part of a long region with the most signal.

    Each axis longer than ``dim_num_thrsh`` (longest first) is cut to ``dim_num_thrsh`` voxels
    around its first position with the most voxels above ``mean + 3 std`` of the region (from its
    start or to its end when that position is near them), within the runs of slices without such
    voxels of at least ``blank_num_trsh`` around it (``None``: none). The mean and the standard
    deviation (of the region as single) are MATLAB's at one thread. Returns the crop and its box
    ``(z0, y0, x0, z1, y1, x1)`` (``(y0, x0, y1, x1)`` for an image); a region within the limit is
    returned as it is with its whole box.
    """
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat import reduce, stats

    a = np.asarray(region)
    if a.ndim not in (2, 3):
        raise ValueError("region must be a (z, y, x) volume or a (y, x) image")
    msz = _msize(a)
    crop = [1, 1, 1, *msz]
    if all(s <= dim_num_thrsh for s in msz):
        return a, _box_of(crop, a)
    # the axes to crop: longest first (a stable sort, as MATLAB's), over the limit
    order = sorted(range(len(msz)), key=lambda d: -msz[d])
    adj = [d + 1 for d in order if msz[d] > dim_num_thrsh]
    mu = reduce.mean(a, num_threads=num_threads)
    sd = stats.std(a.astype(np.float32), num_threads=num_threads)
    # mu + 3 * std in single; the comparison is exact (NumPy compares uint16 or single data with
    # a single value in single, other classes in double)
    thr = np.float32(dtypes.plus(mu, dtypes.times(3.0, sd)))
    bw = a > thr
    axis_of = _AXIS if a.ndim == 3 else {1: 0, 2: 1}
    for ax in adj:
        others = tuple(i for i in range(a.ndim) if i != axis_of[ax])
        fz = np.count_nonzero(bw, axis=others)
        pind = int(np.argmax(fz)) + 1
        s, t = 1, fz.size
        if blank_num_trsh is not None:
            pix = _zero_runs(fz, blank_num_trsh)
            before = pix[pix < pind]
            after = pix[pix > pind]
            s = int(before.max()) if before.size else 1
            t = int(after.min()) if after.size else fz.size
        thrsh = float(dim_num_thrsh)
        if t - s + 1 <= thrsh:
            sa, ta = s, t
        elif pind - s + 1 <= thrsh / 2:
            sa = s
            ta = min(t, sa + thrsh - 1)
        elif t - pind + 1 <= thrsh / 2:
            ta = t
            sa = max(s, ta - thrsh + 1)
        else:
            sa = pind - math.ceil((thrsh - 1) / 2) + 1
            ta = sa + thrsh - 1
        crop[ax - 1] = sa
        crop[ax + 2] = ta
    if len(crop) < 6:
        # crop3d_mex takes 6 values; its MATLAB fallback indexes the 6th
        raise _bad_index(len(crop))
    return a[_slices(crop)].copy(), _box_of(crop, a)


def _slices(crop: Sequence[float]) -> tuple[slice, slice, slice]:
    y0, x0, z0, y1, x1, z1 = (int(v) for v in crop)
    return slice(z0 - 1, z1), slice(y0 - 1, y1), slice(x0 - 1, x1)


def _box_of(crop: Sequence[float], a: np.ndarray) -> tuple[int, ...]:
    if a.ndim == 2:
        return (int(crop[0]) - 1, int(crop[1]) - 1, int(crop[3]), int(crop[4]))
    if len(crop) == 5:
        return (0, int(crop[0]) - 1, int(crop[1]) - 1, 1, int(crop[3]), int(crop[4]))
    y0, x0, z0, y1, x1, z1 = (int(v) for v in crop)
    return (z0 - 1, y0 - 1, x0 - 1, z1, y1, x1)


# ---------------------------------------------------------------------------------------------
# the template's rows, columns and slices
# ---------------------------------------------------------------------------------------------


def cap_region_inds(inds, sz: float, max_off: float, bounds: Sequence[float] = _CAP_BOUNDS):
    """``cap_region_2_inds``: the 1-based positions of ``inds`` (a mask, or positions); above 5
    of them, only those inside ``(max(sz b1, min(sz b2, max_off / 2)), min(sz (1 - b1), max(sz
    (1 - b2), sz - max_off / 2)))``."""
    v = np.asarray(inds)
    pos = np.flatnonzero(v) + 1 if v.dtype == np.bool_ else v.astype(np.int64)
    if pos.size > 5:
        b1, b2 = float(bounds[0]), float(bounds[1])
        lo = max(sz * b1, min(sz * b2, max_off / 2))
        hi = min(sz * (1 - b1), max(sz * (1 - b2), sz - max_off / 2))
        pos = pos[(pos > lo) & (pos < hi)]
    return pos


def template_indices(nonzero: np.ndarray, fraction: float, two_d: bool = False):
    """The template's slices (dimension 3), columns (2) and rows (1) with data in more than
    ``1 / fraction`` of their area, each counted within the ones chosen before, as masks along
    MATLAB's dimensions 1, 2 and 3 (``(y, x, z)``); ``two_d``: the first slice only."""
    nz3, n1, n2 = nonzero.shape
    if two_d:
        z = np.zeros(nz3, dtype=bool)
        z[0] = True
        zc = 1
    else:
        z = np.count_nonzero(nonzero, axis=(1, 2)) > n1 * n2 / fraction
        zc = int(np.count_nonzero(z))
    sub = nonzero[z]
    x = np.count_nonzero(sub, axis=(0, 1)) > n1 * zc / fraction
    y = np.count_nonzero(sub[:, :, x], axis=(0, 2)) > int(np.count_nonzero(x)) * zc / fraction
    return y, x, z
