"""Boxes and ranges of the Python API: 0-based, end-exclusive, in ``(z, y, x)`` order.

A box is ``(z0, y0, x0, z1, y1, x1)`` and selects ``vol[z0:z1, y0:y1, x0:x1]``; a 2D box is
``(y0, x0, y1, x1)``. A range is ``(start, stop)``. An open end is ``None`` (MATLAB's ``Inf``).

MATLAB's box ``[ymin xmin zmin ymax xmax zmax]`` (1-based, inclusive) is
``(zmin-1, ymin-1, xmin-1, zmax, ymax, xmax)``; its 2D box ``[ymin xmin ymax xmax]`` is
``(ymin-1, xmin-1, ymax, xmax)``; its range ``[first last]`` is ``(first-1, last)``. Only the
MATLAB-name layer and files in MATLAB's formats hold MATLAB's values; convert with
``box_from_matlab`` / ``box_to_matlab`` and ``range_from_matlab`` / ``range_to_matlab``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

__all__ = [
    "Box",
    "box_from_matlab",
    "box_shape",
    "box_slices",
    "box_to_matlab",
    "boxes_from_matlab",
    "boxes_to_matlab",
    "check_box",
    "clamp_box",
    "crop",
    "full_box",
    "intersect_boxes",
    "pad_widths",
    "range_from_matlab",
    "range_to_matlab",
]

Box = tuple  # (z0, y0, x0, z1, y1, x1) or (y0, x0, y1, x1); ends may be None

# MATLAB box position of each Python box position: 3D (z0 y0 x0 z1 y1 x1) <- [y x z y x z]
_FROM_M3 = (2, 0, 1, 5, 3, 4)
_TO_M3 = (1, 2, 0, 4, 5, 3)


def _int(v, what: str) -> int:
    f = float(v)
    if not math.isfinite(f) or f != math.floor(f):
        raise ValueError(f"{what} must be integers, got {v!r}")
    return int(f)


def _end_from_matlab(v):
    f = float(v)
    if f == math.inf:
        return None
    return _int(f, "box ends")


def _end_to_matlab(v):
    return math.inf if v is None else int(v)


def box_from_matlab(b) -> Box:
    """MATLAB's ``[ymin xmin zmin ymax xmax zmax]`` (or 2D ``[ymin xmin ymax xmax]``) as our box
    (``Inf`` ends become ``None``)."""
    v = list(np.asarray(b, dtype=np.float64).ravel())
    if len(v) == 6:
        m = v
        return tuple(_int(m[i], "box starts") - 1 for i in (2, 0, 1)) + tuple(
            _end_from_matlab(m[i]) for i in (5, 3, 4)
        )
    if len(v) == 4:
        return (
            _int(v[0], "box starts") - 1,
            _int(v[1], "box starts") - 1,
            _end_from_matlab(v[2]),
            _end_from_matlab(v[3]),
        )
    raise ValueError(
        f"a MATLAB box has 6 values [ymin xmin zmin ymax xmax zmax] (or 4 in 2D), got {len(v)}"
    )


def box_to_matlab(box) -> tuple:
    """Our box as MATLAB's ``[ymin xmin zmin ymax xmax zmax]`` (2D: ``[ymin xmin ymax xmax]``);
    ``None`` ends become ``Inf``."""
    b = _box_values(box)
    if len(b) == 6:
        z0, y0, x0, z1, y1, x1 = b
        return (y0 + 1, x0 + 1, z0 + 1, _end_to_matlab(y1), _end_to_matlab(x1), _end_to_matlab(z1))
    y0, x0, y1, x1 = b
    return (y0 + 1, x0 + 1, _end_to_matlab(y1), _end_to_matlab(x1))


def boxes_from_matlab(rows) -> np.ndarray:
    """Rows of MATLAB boxes (an ``(n, 6)`` array, finite) as an ``(n, 6)`` int64 array of ours."""
    m = np.asarray(rows)
    if m.ndim != 2 or m.shape[1] not in (4, 6):
        raise ValueError(f"expected an (n, 6) or (n, 4) array of MATLAB boxes, got shape {m.shape}")
    if not np.all(np.isfinite(m)) or np.any(m != np.floor(m)):
        raise ValueError("MATLAB boxes must be integral (open ends: use box_from_matlab)")
    m = m.astype(np.int64)
    if m.shape[1] == 6:
        out = m[:, list(_FROM_M3)].copy()
        out[:, :3] -= 1
    else:
        out = m.copy()
        out[:, :2] -= 1
    return out


def boxes_to_matlab(rows) -> np.ndarray:
    """An ``(n, 6)`` (or ``(n, 4)``) array of our boxes as MATLAB's rows (int64)."""
    b = np.asarray(rows).astype(np.int64)
    if b.ndim != 2 or b.shape[1] not in (4, 6):
        raise ValueError(f"expected an (n, 6) or (n, 4) array of boxes, got shape {b.shape}")
    if b.shape[1] == 6:
        out = b[:, list(_TO_M3)].copy()
        out[:, :3] += 1
    else:
        out = b.copy()
        out[:, :2] += 1
    return out


def range_from_matlab(r) -> tuple:
    """MATLAB's ``[first last]`` (1-based, inclusive; ``last`` may be Inf) as ``(start, stop)``."""
    v = list(np.asarray(r, dtype=np.float64).ravel())
    if len(v) != 2:
        raise ValueError(f"a MATLAB range has 2 values [first last], got {len(v)}")
    return (_int(v[0], "range starts") - 1, _end_from_matlab(v[1]))


def range_to_matlab(r) -> tuple:
    """``(start, stop)`` as MATLAB's ``[first last]`` (``None`` stop becomes ``Inf``)."""
    if len(r) != 2:
        raise ValueError(f"a range has 2 values (start, stop), got {len(r)}")
    start, stop = r
    return (
        _int(start, "range starts") + 1,
        _end_to_matlab(None if stop is None else _int(stop, "range stops")),
    )


def _box_values(box) -> tuple:
    if box is None:
        raise ValueError("no box")
    v = list(box.tolist() if isinstance(box, np.ndarray) else box)
    if len(v) not in (4, 6):
        raise ValueError(
            f"a box has 6 values (z0, y0, x0, z1, y1, x1) or 4 (y0, x0, y1, x1), got {len(v)}"
        )
    h = len(v) // 2
    starts = [_int(s, "box starts") for s in v[:h]]
    ends = [None if e is None else _int(e, "box ends") for e in v[h:]]
    return tuple(starts) + tuple(ends)


def _shape_like(shape: Sequence[int] | None, h: int):
    if shape is None:
        return None
    s = [int(v) for v in shape]
    if len(s) < h:
        raise ValueError(f"a {h}D box needs a shape of at least {h} dimensions, got {s}")
    return s[-h:]


def check_box(box, shape: Sequence[int] | None = None, *, allow_outside: bool = False) -> Box:
    """A validated box (a tuple of ints, ``None`` for open ends): ``start <= stop`` per axis and,
    with ``shape`` (the last 2 or 3 dimensions are used), inside ``[0, shape]`` unless
    ``allow_outside`` (the crop wrappers' pad option)."""
    b = _box_values(box)
    h = len(b) // 2
    for i in range(h):
        if b[h + i] is not None and b[i] > b[h + i]:
            raise ValueError(f"box {list(b)} has a start after its stop")
    s = _shape_like(shape, h)
    if s is not None and not allow_outside:
        for i in range(h):
            stop = s[i] if b[h + i] is None else b[h + i]
            if b[i] < 0 or stop > s[i]:
                raise ValueError(f"box {list(b)} is outside an array of shape {s}")
    return b


def full_box(shape: Sequence[int]) -> Box:
    """The box of a whole volume ``(0, 0, 0, nz, ny, nx)`` (an image's: ``(0, 0, ny, nx)``)."""
    s = [int(v) for v in shape]
    if len(s) not in (2, 3):
        raise ValueError(f"expected a shape of 2 or 3 dimensions, got {s}")
    return tuple([0] * len(s) + s)


def box_shape(box, shape: Sequence[int] | None = None) -> tuple[int, ...]:
    """``(z1 - z0, y1 - y0, x1 - x0)``; open ends need ``shape``."""
    b = _box_values(box)
    h = len(b) // 2
    s = _shape_like(shape, h)
    out = []
    for i in range(h):
        stop = b[h + i]
        if stop is None:
            if s is None:
                raise ValueError("a box with an open end needs the shape")
            stop = s[i]
        out.append(stop - b[i])
    return tuple(out)


def box_slices(box) -> tuple[slice, ...]:
    """The slices ``(slice(z0, z1), slice(y0, y1), slice(x0, x1))`` of a box (open ends kept)."""
    b = _box_values(box)
    h = len(b) // 2
    return tuple(slice(b[i], b[h + i]) for i in range(h))


def crop(a: np.ndarray, box) -> np.ndarray:
    """``a[z0:z1, y0:y1, x0:x1]`` (a view; leading axes kept), checked to be inside ``a``."""
    b = check_box(box, a.shape)
    return a[(Ellipsis, *box_slices(b))]


def clamp_box(box, shape: Sequence[int]) -> Box:
    """The box limited to ``[0, shape]`` (open ends closed); ``ValueError`` if nothing is left."""
    b = _box_values(box)
    h = len(b) // 2
    s = _shape_like(shape, h)
    starts = [max(0, b[i]) for i in range(h)]
    stops = [s[i] if b[h + i] is None else min(s[i], b[h + i]) for i in range(h)]
    if any(starts[i] >= stops[i] for i in range(h)):
        raise ValueError(f"box {list(b)} does not overlap an array of shape {s}")
    return tuple(starts + stops)


def intersect_boxes(a, b) -> Box | None:
    """The intersection of two boxes of the same dimension, or None if it is empty."""
    x, y = _box_values(a), _box_values(b)
    if len(x) != len(y):
        raise ValueError("boxes of different dimensions")
    h = len(x) // 2
    starts = [max(x[i], y[i]) for i in range(h)]
    stops = []
    for i in range(h):
        e = [v for v in (x[h + i], y[h + i]) if v is not None]
        stops.append(min(e) if e else None)
    if any(stops[i] is not None and starts[i] >= stops[i] for i in range(h)):
        return None
    return tuple(starts + stops)


def pad_widths(box, shape: Sequence[int]) -> tuple[tuple[int, int], ...]:
    """Pre and post padding per axis for a box that reaches outside ``[0, shape]``, as ``np.pad``
    takes them: ``np.pad(crop(a, clamp_box(box, a.shape)), pad_widths(box, a.shape))`` has the
    box's shape."""
    b = _box_values(box)
    h = len(b) // 2
    s = _shape_like(shape, h)
    return tuple(
        (max(0, -b[i]), 0 if b[h + i] is None else max(0, b[h + i] - s[i])) for i in range(h)
    )
