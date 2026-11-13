"""The edge erosion of the decon masks (``decon_mask_edge_erosion.m``).

``decon_mask_edge_erosion(mask, edge_erosion)`` returns the mask with its edges eroded, as the
decon functions apply it to the mask of nonzero input voxels before they multiply the result by
it (``EdgeErosion``). The erosion is taken as a double (MATLAB's ``double(EdgeErosion)``: single
and integer values give the same masks as the double).

- A full mask (every voxel true) gets its first and last ``edge_erosion`` rows, columns and (3D)
  slices set false directly: MATLAB's ``unique([1 : E, end - E + 1 : end])`` kept within the
  mask, so a dimension shorter than the erosion is all edge; MATLAB's index error where the
  indices are not integers (a fractional erosion).
- Other masks are eroded by the decomposition of ``strel('disk', edge_erosion - 1)`` in each z
  slice and by a line of ``2 * edge_erosion - 1`` voxels along z, with the first and last rows,
  columns and slices set false first (``kernels.decon.mask_edge_erosion``, the mex's kernel, for
  erosions the mex takes: an integer of at least 1; MATLAB's ``imerode`` code for others, with
  strel's errors).

Masks are ``(z, y, x)`` volumes or ``(y, x)`` images; a volume with one z slice is MATLAB's 2D
mask (its size has no third dimension). Numeric masks count their positive values. The input is
never changed.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

__all__ = ["MatlabError", "decon_mask_edge_erosion", "disk_decomposition"]

# MATLAB's error of an index that is not a positive integer
_BAD_INDEX = "MATLAB:badsubscript"


class MatlabError(ValueError):
    """An error MATLAB raises, with its identifier."""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


_DISK: dict[int, tuple[list[np.ndarray], tuple[int, int], tuple[int, int]]] = {}


def disk_decomposition(r: int) -> tuple[list[np.ndarray], tuple[int, int], tuple[int, int]]:
    """The neighborhoods ((rows y, columns x) bool arrays) of the decomposition of
    ``strel('disk', r)`` that imerode applies in turn, and the padding (before, after) of the
    image for them, each (rows, columns): morphop's ``getpadsize`` with before and after swapped
    for erosion when the decomposition has several non-orthogonal parts, else none."""
    if r not in _DISK:
        from petakit5d.matlabcompat.morphology import strel

        se = strel("disk", r)
        seq = se.decompose()
        nhoods = [np.ascontiguousarray(s.nhood) for s in seq]
        pad_ul, pad_lr = (0, 0), (0, 0)
        if len(seq) > 1 and not se.is_decomposition_orthogonal():
            before, after = se.getpadsize()
            pad_ul, pad_lr = (int(after[0]), int(after[1])), (int(before[0]), int(before[1]))
        _DISK[r] = (nhoods, pad_ul, pad_lr)
    nhoods, pad_ul, pad_lr = _DISK[r]
    return list(nhoods), pad_ul, pad_lr


def _logical(mask: Any) -> tuple[np.ndarray, bool]:
    """the mask as bool (numeric: ``mask > 0``) and whether it is a new array"""
    a = np.asarray(mask)
    if a.dtype == np.bool_:
        return a, False
    if np.iscomplexobj(a):
        a = a.real
    if a.dtype.kind not in "iuf":
        raise TypeError(f"the mask must be logical or numeric, not {a.dtype}")
    from petakit5d.kernels import ops

    return ops.compare(a, "gt", 0.0), True


def _scalar(edge_erosion: Any) -> float:
    """MATLAB's ``double(EdgeErosion)`` of a real numeric scalar"""
    e = np.asarray(edge_erosion)
    if e.size != 1 or e.dtype.kind not in "biuf":
        raise TypeError("EdgeErosion must be a real numeric scalar")
    return float(e.reshape(-1)[0])


def _colon(a: float, b: float) -> list[float]:
    """MATLAB's ``a : b`` of doubles (a NaN end gives NaN)"""
    if math.isnan(a) or math.isnan(b):
        return [math.nan]
    if a > b:
        return []
    if math.isinf(b - a):
        raise MatlabError("MATLAB:pmaxsize", "Maximum variable size allowed is exceeded.")
    return [a + k for k in range(int(math.floor(b - a)) + 1)]


def _edge_rows(n: int, e: float, position: int) -> np.ndarray:
    """0-based indices of MATLAB's ``unique(min(max([1 : e, end - e + 1 : end], 1), end))``
    along a dimension of n (min and max ignore NaN); MATLAB's error where they are not integers"""
    v = set()
    for x in _colon(1.0, e) + _colon(n - e + 1, float(n)):
        v.add(min(1.0 if math.isnan(x) else max(x, 1.0), float(n)))
    if any(x != math.floor(x) for x in v):
        raise MatlabError(
            _BAD_INDEX,
            f"Index in position {position} is invalid. Array indices must be positive integers "
            "or logical values.",
        )
    return np.array(sorted(int(x) - 1 for x in v), dtype=np.int64)


def _full_mask(a: np.ndarray, e: float, ndims: int) -> np.ndarray:
    """a full mask's first and last e rows, columns (and slices) false"""
    rows = _edge_rows(a.shape[-2], e, 1)
    cols = _edge_rows(a.shape[-1], e, 2)
    planes = _edge_rows(a.shape[0], e, 3) if ndims == 3 else None
    out = a.copy()
    out[..., rows, :] = False
    out[..., :, cols] = False
    if planes is not None:
        out[planes] = False
    return out


def _disk_radius(e: float) -> int:
    """strel('disk', e - 1)'s check of its radius (MATLAB's validateattributes: an integer, NaN
    and Inf not, then nonnegative)"""
    r = e - 1
    if not math.isfinite(r) or r != math.floor(r):
        raise MatlabError(
            "MATLAB:strel:expectedInteger", "Expected input number 2, R, to be integer-valued."
        )
    if r < 0:
        raise MatlabError(
            "MATLAB:strel:expectedNonnegative", "Expected input number 2, R, to be nonnegative."
        )
    return int(r)


def _imerode_code(a: np.ndarray, e: float, ndims: int, threads: int | None) -> np.ndarray:
    """decon_mask_edge_erosion.m's imerode code (erosions the mex doesn't take)"""
    from petakit5d.matlabcompat.morphology import Strel, imerode, strel

    out = a.copy()
    out[..., 0, :] = False
    out[..., -1, :] = False
    out[..., :, 0] = False
    out[..., :, -1] = False
    if ndims == 3:
        out[0] = False
        out[-1] = False
    r = _disk_radius(e)
    out = imerode(out, strel("disk", r), threads=threads)
    if ndims == 3:
        out = imerode(out, Strel(np.ones((2 * r + 1, 1, 1), dtype=bool)), threads=threads)
    return out


def _mex_takes(e: float) -> bool:
    """whether the mex erodes with e (an integer of at least 1)"""
    return e >= 1 and e == math.floor(e) and math.isfinite(e)


def decon_mask_edge_erosion(mask: Any, edge_erosion: Any, *, threads: int | None = None):
    """The mask with its edges eroded by ``edge_erosion`` (decon_mask_edge_erosion.m): see the
    module. ``mask``: a bool or numeric (z, y, x) volume or (y, x) image; returns a new bool
    array of its shape."""
    a, _ = _logical(mask)
    if a.ndim not in (2, 3):
        raise ValueError("decon_mask_edge_erosion: a (z, y, x) volume or a (y, x) image")
    e = _scalar(edge_erosion)
    ndims = 3 if a.ndim == 3 and a.shape[0] > 1 else 2
    if bool(np.all(a)):
        return _full_mask(a, e, ndims)
    if _mex_takes(e):
        from petakit5d.kernels.decon import mask_edge_erosion

        nhoods, pad_ul, pad_lr = disk_decomposition(int(e) - 1)
        rz = int(e) - 1 if ndims == 3 else 0
        return mask_edge_erosion(a, nhoods, pad_ul, pad_lr, rz, True, threads=threads)
    return _imerode_code(a, e, ndims, threads)
