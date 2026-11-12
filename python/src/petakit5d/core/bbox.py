"""MATLAB bounding boxes: ``[ymin, xmin, zmin, ymax, xmax, zmax]``, 1-based and inclusive.

Parameters, configs and the oracle pass boxes with MATLAB's values; code converts them to
slices here and nowhere else. A box of a 2D image has ``zmin == zmax == 1``.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = [
    "BBox",
    "bbox_intersect",
    "bbox_size",
    "bbox_to_slices",
    "check_bbox",
    "clamp_bbox",
    "crop",
    "full_bbox",
    "pad_widths",
    "slices_to_bbox",
]

BBox = tuple[int, int, int, int, int, int]


def _as_int(v, what: str) -> int:
    f = float(v)
    if f != np.floor(f) or not np.isfinite(f):
        raise ValueError(f"{what} must be integers, got {v!r}")
    return int(f)


def check_bbox(bbox, size: Sequence[int] | None = None, *, allow_outside: bool = False) -> BBox:
    """A validated box as a tuple of 6 ints.

    Values must be integral and ``min <= max`` on each axis. With ``size`` given (y, x, z;
    a 2-element size means nz = 1), the box must lie inside ``[1, size]`` unless
    ``allow_outside`` (the crop wrappers' pad option).
    """
    b = np.asarray(bbox).ravel()
    if b.size != 6:
        raise ValueError(
            f"a bounding box has 6 values [ymin, xmin, zmin, ymax, xmax, zmax], got {b.size}"
        )
    v = tuple(_as_int(x, "bounding box values") for x in b)
    if any(v[i] > v[i + 3] for i in range(3)):
        raise ValueError(f"bounding box {list(v)} has a start after its end")
    if size is not None and not allow_outside:
        sz = _size3(size)
        if any(v[i] < 1 or v[i + 3] > sz[i] for i in range(3)):
            raise ValueError(f"bounding box {list(v)} is outside the volume of size {list(sz)}")
    return v  # type: ignore[return-value]


def _size3(size: Sequence[int]) -> tuple[int, int, int]:
    s = [int(x) for x in size]
    if len(s) == 2:
        s.append(1)
    if len(s) != 3:
        raise ValueError(f"a volume size has 2 or 3 values (y, x, z), got {len(s)}")
    return s[0], s[1], s[2]


def full_bbox(size: Sequence[int]) -> BBox:
    """The box of a whole volume: ``[1, 1, 1, ny, nx, nz]``."""
    ny, nx, nz = _size3(size)
    return (1, 1, 1, ny, nx, nz)


def bbox_size(bbox) -> tuple[int, int, int]:
    """``bbox(4:6) - bbox(1:3) + 1``."""
    b = check_bbox(bbox)
    return (b[3] - b[0] + 1, b[4] - b[1] + 1, b[5] - b[2] + 1)


def bbox_to_slices(bbox) -> tuple[slice, slice, slice]:
    """0-based half-open slices ``(y, x, z)`` of a MATLAB box: ``a(ys:yt, xs:xt, zs:zt)``."""
    b = check_bbox(bbox)
    return (slice(b[0] - 1, b[3]), slice(b[1] - 1, b[4]), slice(b[2] - 1, b[5]))


def slices_to_bbox(slices: Sequence[slice], size: Sequence[int] | None = None) -> BBox:
    """The MATLAB box of 0-based slices with step 1; open ends need ``size``."""
    sl = list(slices)
    if len(sl) == 2:
        sl.append(slice(0, 1))
    if len(sl) != 3:
        raise ValueError("expected 2 or 3 slices (y, x, z)")
    sz = _size3(size) if size is not None else (None, None, None)
    out_s, out_t = [], []
    for s, n in zip(sl, sz, strict=True):
        if s.step not in (None, 1):
            raise ValueError("slices with a step have no bounding box")
        if n is not None:
            start, stop, _ = s.indices(n)
        else:
            if s.start is None and s.stop is None:
                raise ValueError("an open slice needs the volume size")
            start = 0 if s.start is None else s.start
            stop = s.stop
            if stop is None or start < 0 or stop < 0:
                raise ValueError("an open or negative slice needs the volume size")
        if stop <= start:
            raise ValueError("an empty slice has no bounding box")
        out_s.append(start + 1)
        out_t.append(stop)
    return (out_s[0], out_s[1], out_s[2], out_t[0], out_t[1], out_t[2])


def crop(a: np.ndarray, bbox) -> np.ndarray:
    """``a(ys:yt, xs:xt, zs:zt)`` as a view; a 2D array takes a box with z range 1..1."""
    sl = bbox_to_slices(bbox)
    if a.ndim == 2:
        b = check_bbox(bbox)
        if b[2] != 1 or b[5] != 1:
            raise ValueError("a 2D image takes a box with zmin = zmax = 1")
    check_bbox(bbox, a.shape if a.ndim in (2, 3) else None)
    return a[sl[0], sl[1]] if a.ndim == 2 else a[sl]


def clamp_bbox(bbox, size: Sequence[int]) -> BBox:
    """The box limited to ``[1, size]`` on each axis (raises if nothing is left)."""
    b = np.asarray(bbox).ravel()
    if b.size != 6:
        raise ValueError("a bounding box has 6 values")
    sz = _size3(size)
    v = [_as_int(x, "bounding box values") for x in b]
    s = [max(1, v[i]) for i in range(3)]
    t = [min(sz[i], v[i + 3]) for i in range(3)]
    if any(s[i] > t[i] for i in range(3)):
        raise ValueError(f"bounding box {v} does not overlap the volume of size {list(sz)}")
    return (s[0], s[1], s[2], t[0], t[1], t[2])


def bbox_intersect(a, b) -> BBox | None:
    """The intersection of two boxes, or None if they don't overlap."""
    x = check_bbox(a)
    y = check_bbox(b)
    s = [max(x[i], y[i]) for i in range(3)]
    t = [min(x[i + 3], y[i + 3]) for i in range(3)]
    if any(s[i] > t[i] for i in range(3)):
        return None
    return (s[0], s[1], s[2], t[0], t[1], t[2])


def pad_widths(bbox, size: Sequence[int]) -> tuple[tuple[int, int], ...]:
    """Pre and post padding per axis for a box that reaches outside ``[1, size]``.

    ``((pre_y, post_y), (pre_x, post_x), (pre_z, post_z))``, as ``np.pad`` takes them: the
    crop of ``clamp_bbox(bbox, size)`` padded by these has the size of ``bbox``.
    """
    b = check_bbox(bbox, allow_outside=True)
    sz = _size3(size)
    return tuple((max(0, 1 - b[i]), max(0, b[i + 3] - sz[i])) for i in range(3))
