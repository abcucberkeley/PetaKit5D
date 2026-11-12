"""MATLAB's structuring elements, binary morphology, convex hulls and median filters, exactly
(Image Processing Toolbox R2024b, measured with oracle cases).

Arrays are NumPy's: an image is ``(y, x)``, a volume ``(z, y, x)``, C-contiguous; logical
arrays are ``bool``. A neighborhood follows the same rule (MATLAB's ``rows x cols x planes`` is
``(planes, rows, cols)``); per-axis values (sizes, pads, offsets) are ``(z, y, x)`` or ``(y, x)``.

- ``strel``: structuring elements as ``strel.m`` builds them, with their neighborhoods and
  decompositions: ``strel('disk', r)`` is a sequence of periodic lines for r >= 3 (N = 4 by
  default, plus a horizontal and a vertical line when the lines fall short of r), not a
  Euclidean disk; ``strel('sphere', r)`` is ``(x/r)^2 + (y/r)^2 + (z/r)^2 <= 1`` in double.
- ``imerode`` and ``imdilate`` of logical images: erosion is the AND over the neighborhood with
  pixels outside the image counting as true, dilation the OR over the reflected neighborhood
  with pixels outside false, about the center ``floor((size + 1) / 2)`` (MATLAB's, 1-based). A
  decomposed strel is applied as ``morphop`` applies it: element by element, on the image
  padded (``getpadsize``, with the outside value) when the decomposition is not orthogonal,
  then cropped. ``morphop_fast``'s paths (Halide, 2D full neighborhoods) use the whole
  neighborhood; both are followed, and they agree for the strels ``strel`` makes.
- ``imopen`` is ``imdilate(imerode(A, se), se)``; ``imclose`` pads the image with false by
  ``ceil(size(nhood) / 2)``, dilates, erodes and crops, so it is the closing of the image
  extended with false (MATLAB measured on random neighborhoods, Halide and morphop paths).
- ``bwconvhull``: regionprops' ``ConvexImage`` (the hull of the pixel-edge midpoints filled
  with poly2mask's rules) in MATLAB's rows and columns; ``bwconvhull_slices`` for the slices of
  a volume along an axis.
- ``medfilt2`` (zero padding by default) and ``medfilt3`` (symmetric padding by default): the
  median of each window, exact for data without NaN and -0. With NaN or both signs of zero in a
  window MATLAB's result depends on its code path and on the positions in the window; here a
  window with NaN gives its first NaN in MATLAB's order and a zero result is +0 (proposed
  exception, see the module tests).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.core.arrays import from_matlab, size_to_matlab, to_matlab

__all__ = [
    "Strel",
    "bwconvhull",
    "bwconvhull_slices",
    "imclose",
    "imdilate",
    "imerode",
    "imopen",
    "medfilt2",
    "medfilt3",
    "strel",
]


class MatlabError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# structuring elements (built in MATLAB's form, (rows, cols[, planes]) Fortran, as strel.m does)
# ---------------------------------------------------------------------------------------------


def _mshape(shape: Sequence[int]) -> tuple[int, ...]:
    """MATLAB's size: at least 2 dimensions, no trailing singletons beyond the second"""
    s = list(shape) + [1] * (2 - len(shape))
    while len(s) > 2 and s[-1] == 1:
        s.pop()
    return tuple(int(v) for v in s)


def _as_mnhood(nhood: Any) -> np.ndarray:
    """a neighborhood in MATLAB's form (bool, MATLAB's size; dropping trailing singletons
    keeps the element order)"""
    a = np.asarray(nhood)
    if a.ndim == 0:
        a = a.reshape(1, 1)
    elif a.ndim == 1:
        a = a.reshape(1, -1)
    a = a != 0
    return a.reshape(_mshape(a.shape))


def _py_nhood(nhood: Any) -> np.ndarray:
    """a neighborhood given in Python's form as MATLAB's form"""
    a = np.asarray(nhood)
    if a.ndim == 0:
        a = a.reshape(1, 1)
    elif a.ndim == 1:
        a = a.reshape(1, -1)
    if a.ndim > 3:
        raise NotImplementedError("neighborhoods of more than 3 dimensions are not supported")
    return _as_mnhood(to_matlab(a))


def _mstrel(mnhood: Any, decomposition: Sequence[Strel] = (), type: str = "arbitrary") -> Strel:
    """a Strel from a neighborhood in MATLAB's form"""
    se = Strel.__new__(Strel)
    se._m = _as_mnhood(mnhood)
    se._m.flags.writeable = False
    se.decomposition = tuple(decomposition)
    se.type = type
    se._minkowski_exact = None
    se._py = None
    return se


def _py_vector(v: np.ndarray) -> tuple[int, ...]:
    """a per-dimension vector of MATLAB's (rows, cols[, planes]) as Python's ((planes,) rows,
    cols)"""
    v = [int(x) for x in v]
    return tuple(v[2:3] + v[:2]) if len(v) >= 3 else tuple(v)


class Strel:
    """A flat structuring element as MATLAB's ``strel`` object holds it.

    ``nhood``: the neighborhood (bool, C-contiguous, ``(y, x)`` or ``(z, y, x)``);
    ``decomposition``: the sequence of strels it is applied as (empty: the neighborhood
    itself); ``type``: the shape name. Build them with ``strel(...)``; ``Strel(nhood)`` is an
    element without decomposition.
    """

    __slots__ = ("_m", "_minkowski_exact", "_py", "decomposition", "type")

    def __init__(self, nhood: Any, decomposition: Sequence[Strel] = (), type: str = "arbitrary"):
        m = _mstrel(_py_nhood(nhood), decomposition, type)
        self._m = m._m
        self.decomposition = m.decomposition
        self.type = type
        self._minkowski_exact = None
        self._py = None

    def __repr__(self) -> str:
        return (
            f"Strel(type={self.type!r}, size={self.nhood.shape}, neighbors={int(self._m.sum())},"
            f" decomposition={len(self.decomposition)})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Strel):
            return NotImplemented
        return (
            self._m.shape == other._m.shape
            and bool(np.array_equal(self._m, other._m))
            and self.decomposition == other.decomposition
        )

    __hash__ = None  # compared by value

    @property
    def nhood(self) -> np.ndarray:
        """the neighborhood (``se.Neighborhood``) in Python's form"""
        if self._py is None:
            self._py = from_matlab(self._m)
            self._py.flags.writeable = False
        return self._py

    @property
    def neighborhood(self) -> np.ndarray:
        return self.nhood

    def getnhood(self) -> np.ndarray:
        return self.nhood

    @property
    def matlab_nhood(self) -> np.ndarray:
        """the neighborhood in MATLAB's form (rows, cols[, planes]), Fortran order"""
        return self._m

    @property
    def dimensionality(self) -> int:
        """MATLAB's ``se.Dimensionality`` (``ndims`` of the neighborhood)"""
        return self._m.ndim

    def decompose(self) -> list[Strel]:
        """``decompose(se)`` / ``getsequence(se)``: the strels applied in turn"""
        return list(self.decomposition) if self.decomposition else [self]

    def getneighbors(self) -> np.ndarray:
        """``getneighbors(se)``: the offsets of the neighbors from the center, one row per
        neighbor in MATLAB's order (column-major in MATLAB's form), columns ``(dy, dx)`` or
        ``(dz, dy, dx)`` (int64)"""
        o = _neighbors(self._m)
        return o[:, [2, 0, 1]] if o.shape[1] == 3 else o

    def is_decomposition_orthogonal(self) -> bool:
        """``isdecompositionorthogonal(se)``"""
        return _is_orthogonal(self.decompose())

    def getpadsize(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        """``[pad_ul, pad_lr] = getpadsize(se)`` as ``(before, after)``, each ``(y, x)`` or
        ``(z, y, x)``"""
        ul, lr = _padsize(self.decompose())
        return _py_vector(ul), _py_vector(lr)

    def reflect(self) -> Strel:
        """``reflect(se)``: the neighborhood reflected through its center (even sizes grow by
        one at the end so that the center stays a pixel)"""
        nh = self._m[tuple(slice(None, None, -1) for _ in range(self._m.ndim))]
        nh = np.pad(nh, [(0, 1 - n % 2) for n in nh.shape])
        return _mstrel(nh, tuple(s.reflect() for s in self.decomposition), self.type)

    def minkowski_exact(self) -> bool:
        """whether the full dilation of the decomposition is the neighborhood (then the
        decomposition and the neighborhood give the same erosions and dilations)"""
        if self._minkowski_exact is None:
            if not self.decomposition:
                self._minkowski_exact = True
            else:
                full = _full_dilate(np.ones((1, 1), dtype=bool), self.decomposition)
                self._minkowski_exact = full.shape == self._m.shape and bool(
                    np.array_equal(full, self._m)
                )
        return self._minkowski_exact


def _neighbors(nhood: np.ndarray) -> np.ndarray:
    idx = np.flatnonzero(nhood.ravel(order="F"))
    subs = np.unravel_index(idx, nhood.shape, order="F")
    center = [(n + 1) // 2 - 1 for n in nhood.shape]
    if len(idx) == 0:
        return np.zeros((0, nhood.ndim), dtype=np.int64)
    return np.stack([s - c for s, c in zip(subs, center, strict=True)], axis=1).astype(np.int64)


def _is_orthogonal(seq: Sequence[Strel]) -> bool:
    nd = max([2] + [s._m.ndim for s in seq])
    P = np.ones((len(seq), nd), dtype=np.int64)
    for i, s in enumerate(seq):
        P[i, : s._m.ndim] = s._m.shape
    P[P == 0] = 1
    return bool(np.any(np.sum(P != 1, axis=0) == 1))


def _padsize(seq: Sequence[Strel]) -> tuple[np.ndarray, np.ndarray]:
    offs = [_neighbors(s._m) for s in seq]
    if not offs:
        return np.zeros(2, dtype=np.int64), np.zeros(2, dtype=np.int64)
    nd = max(o.shape[1] for o in offs)
    pad_ul = np.zeros(nd, dtype=np.int64)
    pad_lr = np.zeros(nd, dtype=np.int64)
    for o in offs:
        if o.shape[0] == 0:
            continue
        o = np.pad(o, ((0, 0), (0, nd - o.shape[1])))
        pad_ul += np.maximum(0, -o.min(axis=0))
        pad_lr += np.maximum(0, o.max(axis=0))
    return pad_ul, pad_lr


def _ref_step(a: np.ndarray, offsets: np.ndarray, erode: bool) -> np.ndarray:
    """one erosion (AND, outside true) or dilation (OR of the reflected offsets, outside
    false) of a small bool array by offsets (rows of at most a.ndim columns)"""
    nd = a.ndim
    off = np.zeros((offsets.shape[0], nd), dtype=np.int64)
    off[:, : offsets.shape[1]] = offsets
    if not erode:
        off = -off
    out = np.full(a.shape, erode, dtype=bool)
    if off.shape[0] == 0:
        return out
    g = np.abs(off).max(axis=0)
    p = np.pad(a, [(int(v), int(v)) for v in g], constant_values=erode)
    for o in off:
        sl = tuple(slice(int(g[d] + o[d]), int(g[d] + o[d] + a.shape[d])) for d in range(nd))
        if erode:
            out &= p[sl]
        else:
            out |= p[sl]
    return out


def _full_dilate(a: np.ndarray, seq: Sequence[Strel]) -> np.ndarray:
    """``imdilate(a, seq, 'full')`` of a small bool array: padded by getpadsize with false,
    dilated by each strel in turn"""
    pad_ul, pad_lr = _padsize(seq)
    nd = max(a.ndim, len(pad_ul))
    a = a.reshape(a.shape + (1,) * (nd - a.ndim))
    pu = list(pad_ul) + [0] * (nd - len(pad_ul))
    pl = list(pad_lr) + [0] * (nd - len(pad_lr))
    b = np.pad(a, [(int(u), int(v)) for u, v in zip(pu, pl, strict=True)])
    for s in seq:
        b = _ref_step(b, _neighbors(s._m), erode=False)
    return b.reshape(_mshape(b.shape))


def _mround(x: float) -> float:
    """MATLAB's round of a double scalar: half away from zero, exactly (sign of zero kept)"""
    t = float(math.trunc(x))
    if abs(x - t) >= 0.5:
        t += math.copysign(1.0, x)
    return math.copysign(t, x) if t == 0 else t


def _intline(x1: float, x2: float, y1: float, y2: float) -> tuple[np.ndarray, np.ndarray]:
    """``iptui.intline``: the integer points of the segment (x1, y1) - (x2, y2)"""
    dx, dy = abs(x2 - x1), abs(y2 - y1)
    if dx == 0 and dy == 0:
        return np.array([x1], dtype=np.float64), np.array([y1], dtype=np.float64)
    flip = False
    if dx >= dy:
        if x1 > x2:
            x1, x2, y1, y2 = x2, x1, y2, y1
            flip = True
        m = (y2 - y1) / (x2 - x1)
        x = np.arange(x1, x2 + 0.5, dtype=np.float64)
        y = _round_array(y1 + m * (x - x1))
    else:
        if y1 > y2:
            x1, x2, y1, y2 = x2, x1, y2, y1
            flip = True
        m = (x2 - x1) / (y2 - y1)
        y = np.arange(y1, y2 + 0.5, dtype=np.float64)
        x = _round_array(x1 + m * (y - y1))
    if flip:
        x, y = x[::-1], y[::-1]
    return x, y


def _round_array(a: np.ndarray) -> np.ndarray:
    from petakit5d.core.dtypes import mround

    return np.asarray(mround(np.asarray(a, dtype=np.float64)))


def _line_nhood(r: np.ndarray, c: np.ndarray) -> np.ndarray:
    r = r.astype(np.int64)
    c = c.astype(np.int64)
    mr, mc = int(np.abs(r).max()), int(np.abs(c).max())
    nh = np.zeros((2 * mr + 1, 2 * mc + 1), dtype=bool)
    nh[r + mr, c + mc] = True
    return nh


def _strel_line(length: float, deg: float) -> Strel:
    if length >= 1:
        theta = math.fmod(deg, 180.0)
        if theta < 0:
            theta += 180.0
        theta = theta * math.pi / 180
        x = _mround((length - 1) / 2 * math.cos(theta))
        y = -_mround((length - 1) / 2 * math.sin(theta))
        c, r = _intline(-x, x, -y, y)
        return _mstrel(_line_nhood(r, c), (), "line")
    return _mstrel(np.zeros((0, 0), dtype=bool), (), "line")


def _strel_periodicline(p: int, v: Sequence[int]) -> Strel:
    k = np.arange(-p, p + 1, dtype=np.int64)
    return _mstrel(_line_nhood(k * int(v[0]), k * int(v[1])), (), "periodicline")


def _strel_arbitrary(mnhood: Any, type: str = "arbitrary") -> Strel:
    """MakeArbitraryStrel of a neighborhood in MATLAB's form: an all-true neighborhood is
    decomposed into lines along each dimension when that is at least twice cheaper"""
    nh = _as_mnhood(mnhood)
    dec: list[Strel] = []
    if nh.size and nh.all():
        size = nh.shape
        if math.prod(size) / sum(size) >= 2:
            for k in range(nh.ndim):
                sk = [1] * nh.ndim
                sk[k] = size[k]
                dec.append(_mstrel(np.ones(sk, dtype=bool), (), "arbitrary"))
    return _mstrel(nh, dec, type)


def _strel_disk(r: int, n: int = 4) -> Strel:
    if r < 3:
        n = 0
    if n == 0:
        xx, yy = np.meshgrid(np.arange(-r, r + 1), np.arange(-r, r + 1))
        return _mstrel(xx**2 + yy**2 <= r**2, (), "disk")
    v = {
        4: [(1, 0), (1, 1), (0, 1), (-1, 1)],
        6: [(1, 0), (1, 2), (2, 1), (0, 1), (-1, 2), (-2, 1)],
        8: [(1, 0), (2, 1), (1, 1), (1, 2), (0, 1), (-1, 2), (-1, 1), (-2, 1)],
    }[n]
    theta = math.pi / (2 * n)
    k = 2 * r / (1 / math.tan(theta) + 1 / math.sin(theta))
    dec = [_strel_periodicline(math.floor(k / math.sqrt(a * a + b * b)), (a, b)) for a, b in v]
    nhood = _full_dilate(np.ones((1, 1), dtype=bool), dec)
    rd = np.nonzero(nhood)[0]
    rd = rd + 1 - (nhood.shape[0] + 1) // 2
    radial_difference = r - int(rd.max())
    length = 2 * (radial_difference - 1) + 1
    if length >= 3:
        lines = [_strel_line(length, 0), _strel_line(length, 90)]
        dec += lines
        nhood = _full_dilate(nhood, lines)
    return _mstrel(nhood, dec, "disk")


def _strel_sphere(r: int) -> Strel:
    v = np.arange(-r, r + 1, dtype=np.float64)
    y, x, z = np.meshgrid(v, v, v, indexing="ij")
    with np.errstate(invalid="ignore", divide="ignore"):
        nh = ((x / r) ** 2 + (y / r) ** 2) + (z / r) ** 2 <= 1
    return _mstrel(nh, (), "sphere")


def _strel_diamond(m: int) -> Strel:
    rr, cc = np.meshgrid(np.arange(-m, m + 1), np.arange(-m, m + 1))
    nh = (np.abs(rr) + np.abs(cc)) <= m
    dec: list[Strel] = []
    if m > 2:
        n = math.floor(math.log2(m))
        dec.append(_mstrel(np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])))
        sizes = [2 ** (k + 1) + 1 for k in range(n)]
        q = m - 2**n
        if q > 0:
            sizes.append(2 * q + 1)
        for p in sizes:
            mid = (p + 1) // 2 - 1
            e = np.zeros((p, p), dtype=bool)
            e[0, mid] = e[p - 1, mid] = e[mid, 0] = e[mid, p - 1] = True
            dec.append(_mstrel(e))
    return _mstrel(nh, dec, "diamond")


def _strel_octagon(m: int) -> Strel:
    k = m // 3
    rr, cc = np.meshgrid(np.arange(-m, m + 1), np.arange(-m, m + 1))
    nh = np.abs(rr) + np.abs(cc) <= m + k
    dec: list[Strel] = []
    if k >= 1:
        a = np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]], dtype=bool)
        four = [a, a.T, np.eye(3, dtype=bool), np.rot90(np.eye(3, dtype=bool))]
        dec = [_mstrel(e) for _ in range(k) for e in four]
    return _mstrel(nh, dec, "octagon")


def _strel_pair(mn: Sequence[int]) -> Strel:
    size = [abs(int(v)) * 2 + 1 for v in mn]
    nh = np.zeros(size, dtype=bool)
    c = [(s + 1) // 2 - 1 for s in size]
    nh[c[0], c[1]] = True
    nh[c[0] + int(mn[0]), c[1] + int(mn[1])] = True
    return _mstrel(nh, (), "pair")


_TYPES = (
    "arbitrary",
    "square",
    "diamond",
    "rectangle",
    "octagon",
    "line",
    "disk",
    "sphere",
    "cube",
    "cuboid",
)


def _strel_type(name: str) -> str:
    low = name.lower()
    for full in ("ball", "pair", "periodicline"):
        if low and full.startswith(low):
            return full
    hits = [t for t in _TYPES if low and t.startswith(low)]
    if len(hits) == 1 or low in _TYPES:
        return low if low in _TYPES else hits[0]
    raise MatlabError(
        "MATLAB:strel:unrecognizedStringChoice", f"Expected STREL_TYPE to match one of {_TYPES}"
    )


def _nonneg_int(v: Any, name: str) -> int:
    f = float(v)
    if not (f >= 0 and f == math.floor(f) and math.isfinite(f)):
        raise MatlabError(f"MATLAB:strel:expected{name}", f"{name} must be a nonnegative integer")
    return int(f)


def strel(shape: Any, *params: Any) -> Strel:
    """MATLAB's ``strel`` with Python's per-axis order: ``strel(nhood)`` (a ``(y, x)`` or
    ``(z, y, x)`` 0/1 array), ``strel('disk', r[, n])``, ``strel('sphere', r)``,
    ``strel('square', w)``, ``strel('rectangle', (ny, nx))``, ``strel('line', len, deg)``,
    ``strel('periodicline', p, (dy, dx))``, ``strel('diamond', r)``, ``strel('octagon', r)``,
    ``strel('pair', (dy, dx))``, ``strel('cube', w)``, ``strel('cuboid', (nz, ny, nx))``,
    ``strel('arbitrary', nhood)``. Nonflat strels (``'ball'``, heights) are not supported."""
    if not isinstance(shape, str):
        if params:
            raise NotImplementedError("nonflat structuring elements (heights) are not supported")
        return _strel_arbitrary(_py_nhood(shape))
    t = _strel_type(shape)
    if t == "arbitrary":
        if len(params) != 1:
            raise NotImplementedError("strel('arbitrary', nhood) takes one neighborhood (flat)")
        return _strel_arbitrary(_py_nhood(params[0]))
    if t == "ball":
        raise NotImplementedError("strel('ball', ...) is nonflat; not supported")
    if t == "square":
        w = _nonneg_int(params[0], "Size")
        return _strel_arbitrary(np.ones((w, w), dtype=bool), "square")
    if t == "rectangle":
        m, n = (_nonneg_int(v, "Size") for v in params[0])
        return _strel_arbitrary(np.ones((m, n), dtype=bool), "rectangle")
    if t == "cube":
        w = _nonneg_int(params[0], "Width")
        return _strel_arbitrary(np.ones((w, w, w), dtype=bool), "cube")
    if t == "cuboid":
        p, m, n = (_nonneg_int(v, "Size") for v in params[0])
        return _strel_arbitrary(np.ones((m, n, p), dtype=bool), "cuboid")
    if t == "disk":
        r = _nonneg_int(params[0], "R")
        n = 4 if len(params) < 2 else int(params[1])
        if n not in (0, 4, 6, 8):
            raise MatlabError("images:strel:invalidN", "N must be 0, 4, 6, or 8")
        return _strel_disk(r, n)
    if t == "sphere":
        return _strel_sphere(_nonneg_int(params[0], "R"))
    if t == "diamond":
        return _strel_diamond(_nonneg_int(params[0], "Size"))
    if t == "octagon":
        m = _nonneg_int(params[0], "Size")
        if m % 3:
            raise MatlabError("images:strel:notMultipleOf3", "M must be a multiple of 3")
        return _strel_octagon(m)
    if t == "line":
        return _strel_line(float(params[0]), float(params[1]))
    if t == "periodicline":
        return _strel_periodicline(_nonneg_int(params[0], "P"), [int(v) for v in params[1]])
    if t == "pair":
        return _strel_pair([int(v) for v in params[0]])
    raise MatlabError("images:strel:unknownStrelType", f"unknown strel type {shape!r}")


# ---------------------------------------------------------------------------------------------
# erosion and dilation
# ---------------------------------------------------------------------------------------------


def _strels(se: Any, func: str) -> list[Strel]:
    """images.internal.strelcheck: a strel, a sequence of strels, or a 0/1 neighborhood (in
    Python's form)"""
    if isinstance(se, Strel):
        return [se]
    if isinstance(se, (list, tuple)) and se and all(isinstance(s, Strel) for s in se):
        return list(se)
    a = np.asarray(se)
    if a.dtype.kind not in "biuf":
        raise MatlabError("images:strelcheck:invalidStrelType", f"{func}: SE must be a strel")
    if a.size and not np.all((a == 0) | (a == 1)):
        raise MatlabError(
            "images:strelcheck:invalidStrelValues", f"{func}: SE must contain only 0s and 1s"
        )
    return [_strel_arbitrary(_py_nhood(a))]


def _logical(a: Any, func: str) -> np.ndarray:
    a = np.asarray(a)
    if a.dtype != np.bool_:
        raise NotImplementedError(f"{func}: only logical (bool) images are supported")
    if a.ndim not in (2, 3):
        raise NotImplementedError(f"{func}: images must be (y, x) or (z, y, x)")
    return np.ascontiguousarray(a)


def _msize(a: np.ndarray) -> tuple[int, ...]:
    """MATLAB's size of a Python (y, x) or (z, y, x) array"""
    return size_to_matlab(a.shape)


def _fast_nhood(a: np.ndarray, ses: list[Strel], full_2d: bool) -> np.ndarray | None:
    """the neighborhood morphop_fast applies directly (the Halide kernels, or the builtin of 2D
    full neighborhoods when full_2d), or None when morphop applies the decomposition"""
    if len(ses) != 1 or a.size == 0:
        return None
    nh = ses[0]._m
    nd = len(_msize(a))
    if nh.size == 0 or nd > 3 or nd != nh.ndim:
        return None
    nn = int(nh.sum())
    if nd == 2 and nn == nh.size:
        return nh if full_2d else None
    if nn < 600 and all(s <= 15 for s in nh.shape):
        return nh
    return None


def _offsets_zyx(mnhood: np.ndarray) -> np.ndarray:
    """the neighbors of a neighborhood in MATLAB's form as rows (dz, dy, dx)"""
    o = _neighbors(mnhood)
    if o.shape[1] > 3:
        raise NotImplementedError("structuring elements of more than 3 dimensions")
    o = np.pad(o, ((0, 0), (0, 3 - o.shape[1])))
    return np.ascontiguousarray(o[:, [2, 0, 1]], dtype=np.int64)


def _pad_zyx(p: Sequence[int]) -> np.ndarray:
    """a MATLAB pad vector (rows, cols[, planes]) as (z, y, x)"""
    if len(p) > 3:
        raise NotImplementedError("structuring elements of more than 3 dimensions")
    v = [int(x) for x in p] + [0] * (3 - len(p))
    return np.array([v[2], v[0], v[1]], dtype=np.int64)


def _kernel_steps(
    a: np.ndarray,
    mnhoods: Sequence[np.ndarray],
    erode: bool,
    pad: tuple[Sequence[int], Sequence[int]] | None,
    threads: int | None,
) -> np.ndarray:
    """the kernel: a (y, x) or (z, y, x) bool image, neighborhoods in MATLAB's form, pads as
    MATLAB's vectors (before, after)"""
    from petakit5d import _ext
    from petakit5d.core.threads import resolve_threads

    morph = _ext.load("morph")
    offs = [_offsets_zyx(nh) for nh in mnhoods]
    offsets = np.concatenate(offs, axis=0) if offs else np.zeros((0, 3), dtype=np.int64)
    offsets = np.ascontiguousarray(offsets, dtype=np.int64)
    counts = np.array([o.shape[0] for o in offs], dtype=np.int64)
    if pad is None:
        pb = pa = np.zeros(3, dtype=np.int64)
    else:
        pb, pa = _pad_zyx(pad[0]), _pad_zyx(pad[1])
    v = np.ascontiguousarray(a.reshape((1,) * (3 - a.ndim) + a.shape))
    out = np.empty(v.shape, dtype=bool)
    morph.binary_morph(v, out, offsets, counts, erode, pb, pa, resolve_threads(threads))
    return out.reshape(a.shape)


def _apply_nhood(a: np.ndarray, se: Strel, erode: bool, threads: int | None) -> np.ndarray:
    """erosion or dilation by the whole neighborhood (outside true or false), computed with the
    decomposition (padded, which gives the same result) when that is exact and cheaper"""
    seq = list(se.decomposition)
    if (
        len(seq) > 1
        and se.minkowski_exact()
        and sum(int(s._m.sum()) for s in seq) < int(se._m.sum())
    ):
        pad_ul, pad_lr = _padsize(seq)
        pad = (pad_lr, pad_ul) if erode else (pad_ul, pad_lr)
        return _kernel_steps(a, [s._m for s in seq], erode, pad, threads)
    return _kernel_steps(a, [se._m], erode, None, threads)


def _morphop(a: np.ndarray, ses: list[Strel], erode: bool, threads: int | None) -> np.ndarray:
    """images.internal.morphop: the sequence of the decomposed strels, on the image padded by
    getpadsize (swapped for erosion) when the sequence is not orthogonal"""
    seq = [s for se in ses for s in se.decompose()]
    pad = None
    if len(seq) > 1 and not _is_orthogonal(seq):
        pad_ul, pad_lr = _padsize(seq)
        pad = (pad_lr, pad_ul) if erode else (pad_ul, pad_lr)
    return _kernel_steps(a, [s._m for s in seq], erode, pad, threads)


def _erode_dilate(a: Any, se: Any, erode: bool, func: str, threads: int | None) -> np.ndarray:
    a = _logical(a, func)
    ses = _strels(se, func)
    if a.size == 0:
        return a.copy()
    if _fast_nhood(a, ses, full_2d=True) is not None:
        return _apply_nhood(a, ses[0], erode, threads)
    return _morphop(a, ses, erode, threads)


def imerode(a: Any, se: Any, *, threads: int | None = None) -> np.ndarray:
    """MATLAB's ``imerode`` of a logical image ((y, x) or (z, y, x), bool): the AND over the
    neighborhood, pixels outside the image true. ``se``: a ``Strel``, a list of them, or a 0/1
    neighborhood array in Python's form. ``threads``: None for all available CPUs."""
    return _erode_dilate(a, se, True, "imerode", threads)


def imdilate(a: Any, se: Any, *, threads: int | None = None) -> np.ndarray:
    """MATLAB's ``imdilate`` of a logical image: the OR over the reflected neighborhood, pixels
    outside the image false."""
    return _erode_dilate(a, se, False, "imdilate", threads)


def _single_strel(se: Any, func: str) -> list[Strel]:
    ses = _strels(se, func)
    if len(ses) != 1:
        raise MatlabError(f"images:{func}:nonscalarStrel", f"{func}: SE must be a single strel")
    return ses


def imopen(a: Any, se: Any, *, threads: int | None = None) -> np.ndarray:
    """MATLAB's ``imopen`` of a logical image: ``imdilate(imerode(a, se), se)``."""
    a = _logical(a, "imopen")
    ses = _single_strel(se, "imopen")
    if a.size == 0:
        # MATLAB packs a 2D image for a 2D strel (bwpack), and the packed erosion of an image
        # without rows needs its row count: imopen errors there (kept)
        if len(_msize(a)) == 2 and _msize(a)[0] == 0 and ses[0]._m.ndim == 2:
            raise MatlabError(
                "images:morphop:missingPackedM", "M must be provided for packed erosion."
            )
        return a.copy()
    if _fast_nhood(a, ses, full_2d=False) is not None:
        e = _apply_nhood(a, ses[0], True, threads)
        return _apply_nhood(e, ses[0], False, threads)
    return imdilate(imerode(a, ses, threads=threads), ses, threads=threads)


def imclose(a: Any, se: Any, *, threads: int | None = None) -> np.ndarray:
    """MATLAB's ``imclose`` of a logical image: the image padded with false by
    ``ceil(size(nhood) / 2)`` (MATLAB's dimensions, as many as both have), dilated, eroded and
    cropped: the closing of the image extended with false (MATLAB's Halide path gives the
    same)."""
    a = _logical(a, "imclose")
    ses = _single_strel(se, "imclose")
    if a.size == 0:
        return a.copy()
    nh = ses[0]._m
    mpad = [(n + 1) // 2 for n in nh.shape][: min(nh.ndim, len(_msize(a)))]
    mpad += [0] * (3 - len(mpad))
    # MATLAB's (rows, cols, planes) as Python's axes
    pad = [mpad[2], mpad[0], mpad[1]] if a.ndim == 3 else mpad[:2]
    ap = np.pad(a, [(p, p) for p in pad])
    b = imerode(imdilate(ap, ses, threads=threads), ses, threads=threads)
    return np.ascontiguousarray(
        b[tuple(slice(p, p + n) for p, n in zip(pad, a.shape, strict=True))]
    )


# ---------------------------------------------------------------------------------------------
# convex hull images
# ---------------------------------------------------------------------------------------------


def bwconvhull_slices(
    vol: Any, axis: int, indices: Sequence[int], *, threads: int | None = None
) -> np.ndarray:
    """The volume with its slices along ``axis`` at ``indices`` (0-based) replaced by their
    ``bwconvhull(slice, 'union')``, as MATLAB computes the hull of the slice it squeezes out:
    axis 0 (``V(:, :, k)``): rows y, columns x; axis 1 (``squeeze(V(k, :, :))``): rows x,
    columns z; axis 2 (``squeeze(V(:, k, :))``): rows y, columns z. ``vol``: a (z, y, x) bool
    volume with every size above 1 (MATLAB's squeeze changes the slices otherwise). The slices
    of one call are independent; calls for several axes depend on their order."""
    from petakit5d import _ext
    from petakit5d.core.threads import resolve_threads

    v = np.asarray(vol)
    if v.dtype != np.bool_ or v.ndim != 3:
        raise ValueError("bwconvhull_slices: vol must be a (z, y, x) bool volume")
    if min(v.shape) < 2:
        raise ValueError("bwconvhull_slices: every size of vol must be at least 2")
    if axis not in (0, 1, 2):
        raise ValueError("bwconvhull_slices: axis must be 0, 1 or 2")
    out = np.array(v, order="C", copy=True)
    inds = np.ascontiguousarray(np.asarray(indices, dtype=np.int64).reshape(-1))
    _ext.load("morph").bwconvhull_slices(out, int(axis), inds, resolve_threads(threads))
    return out


def bwconvhull(bw: Any, method: str = "union", conn: Any = 8) -> np.ndarray:
    """MATLAB's ``bwconvhull(BW, method, conn)`` of a 2D image (nonzero is true): ``'union'``
    the convex hull image of all true pixels, ``'objects'`` the union of the hull images of
    the components (connectivity 4 or 8)."""
    from petakit5d import _ext

    a = np.asarray(bw)
    if a.ndim != 2:
        raise MatlabError("MATLAB:bwconvhull:expected2D", "bwconvhull: BW must be 2-D")
    a = np.ascontiguousarray(a if a.dtype == np.bool_ else a != 0)
    m = method.lower()
    if m and "union".startswith(m):
        m = "union"
    elif m and "objects".startswith(m):
        m = "objects"
    else:
        raise MatlabError(
            "MATLAB:bwconvhull:unrecognizedStringChoice", "METHOD must be 'union' or 'objects'"
        )
    morph = _ext.load("morph")
    if a.size == 0:
        return np.zeros(a.shape, dtype=bool)
    if m == "union":
        out = np.ascontiguousarray(a.reshape((1,) + a.shape))
        morph.bwconvhull_slices(out, 0, np.zeros(1, dtype=np.int64), 1)
        return out.reshape(a.shape)
    from petakit5d.matlabcompat.label import bwconncomp

    c = np.asarray(conn)
    if c.size != 1 or float(c.reshape(-1)[0]) not in (4, 8):
        if not (c.shape == (3, 3) and np.array_equal(c, c[::-1, ::-1]) and c[1, 1]):
            raise MatlabError(
                "images:bwconvhull:invalidConnectivity", "bwconvhull: invalid connectivity"
            )
    out = np.zeros(a.shape, dtype=bool)
    cc = bwconncomp(a, conn)
    for idx in cc.pixel_idx_list:
        ys, xs = np.unravel_index(idx, a.shape)
        y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
        sub = np.zeros((1, y1 - y0, x1 - x0), dtype=bool)
        sub[0, ys - y0, xs - x0] = True
        morph.bwconvhull_slices(sub, 0, np.zeros(1, dtype=np.int64), 1)
        out[y0:y1, x0:x1] |= sub[0]
    return out


# ---------------------------------------------------------------------------------------------
# median filters
# ---------------------------------------------------------------------------------------------

_ORDER_CLASSES = {
    np.dtype(np.bool_): "bool",
    np.dtype(np.uint8): "uint8",
    np.dtype(np.int8): "int8",
    np.dtype(np.uint16): "uint16",
    np.dtype(np.int16): "int16",
    np.dtype(np.uint32): "uint32",
    np.dtype(np.int32): "int32",
    np.dtype(np.float32): "single",
    np.dtype(np.float64): "double",
}


def _pad_mode(a: np.ndarray, widths: list[tuple[int, int]], padopt: str) -> np.ndarray:
    if padopt == "zeros":
        return np.pad(a, widths, mode="constant", constant_values=0)
    if padopt == "ones":
        return np.pad(a, widths, mode="constant", constant_values=1)
    if padopt == "symmetric":
        return np.pad(a, widths, mode="symmetric")
    if padopt == "replicate":
        return np.pad(a, widths, mode="edge")
    raise ValueError(padopt)


def _order_filter(
    a: np.ndarray, size: tuple[int, int, int], k: int, padopt: str, threads: int | None
) -> np.ndarray:
    """the k-th smallest (0-based) of each window of a (z, y, x) volume, the window placed as
    ordfilt2 places it (center floor((size + 1) / 2), MATLAB's 1-based), padded with padopt"""
    from petakit5d import _ext
    from petakit5d.core.threads import resolve_threads

    morph = _ext.load("morph")
    widths = [((s + 1) // 2 - 1, s - (s + 1) // 2) for s in size]
    p = np.ascontiguousarray(_pad_mode(a, widths, padopt))
    out = np.empty(a.shape, dtype=a.dtype)
    nt = resolve_threads(threads)
    if a.dtype == np.float32 and size == (3, 3, 3) and k == 13:
        # the selection network reads 16 floats past the end of the last line
        buf = np.empty(p.size + 16, dtype=np.float32)
        buf[: p.size] = p.reshape(-1)
        buf[p.size :] = 0
        pv = buf[: p.size].reshape(p.shape)
        if morph.median27_single(pv, out, nt):
            return out
    getattr(morph, "order_filter_" + _ORDER_CLASSES[a.dtype])(p, out, *size, k, nt)
    return out


def _check_class(a: np.ndarray, func: str) -> None:
    if a.dtype not in _ORDER_CLASSES:
        raise NotImplementedError(f"{func}: class {a.dtype} is not supported")


def _imlincomb_half(b: np.ndarray, b2: np.ndarray) -> np.ndarray:
    """imlincomb(0.5, b, 0.5, b2): computed in double, then in b's class (integers rounded half
    away from zero and saturated)"""
    from petakit5d.core.dtypes import cast

    r = 0.5 * b.astype(np.float64) + 0.5 * b2.astype(np.float64)
    if b.dtype == np.float64:
        return r
    if b.dtype == np.float32:
        return r.astype(np.float32)
    return np.asarray(cast(r, _ORDER_CLASSES[b.dtype]))


def _padopt(padopt: str, options: Sequence[str], func: str) -> str:
    p = padopt.lower()
    hits = [o for o in options if p and o.startswith(p)]
    if len(hits) != 1:
        raise MatlabError(
            f"MATLAB:{func}:unrecognizedStringChoice", f"{func}: PADOPT must be one of {options}"
        )
    return hits[0]


def medfilt2(
    a: Any, size: Sequence[int] = (3, 3), padopt: str = "zeros", *, threads: int | None = None
) -> np.ndarray:
    """MATLAB's ``medfilt2(A, [m n], padopt)`` of a (y, x) image: the median of each m x n
    window (``size = (m, n)``), padded with ``'zeros'`` (default), ``'symmetric'`` or
    ``'indexed'`` (ones for double, else zeros); for even m * n the mean of the two middle
    values in double (OR for logical), as MATLAB's imlincomb."""
    a = np.asarray(a)
    if a.ndim != 2:
        raise MatlabError("MATLAB:medfilt2:expected2D", "medfilt2: A must be 2-D")
    _check_class(a, "medfilt2")
    p = _padopt(padopt, ("indexed", "zeros", "symmetric"), "medfilt2")
    if p == "indexed":
        p = "ones" if a.dtype == np.float64 else "zeros"
    m, n = (int(v) for v in size)
    if m < 1 or n < 1:
        raise MatlabError("MATLAB:medfilt2:expectedPositive", "medfilt2: [M N] must be positive")
    if a.size == 0:
        return a.copy()
    v = np.ascontiguousarray(a.reshape((1,) + a.shape))
    if (m * n) % 2:
        return _order_filter(v, (1, m, n), (m * n + 1) // 2 - 1, p, threads).reshape(a.shape)
    b = _order_filter(v, (1, m, n), m * n // 2 - 1, p, threads).reshape(a.shape)
    b2 = _order_filter(v, (1, m, n), m * n // 2, p, threads).reshape(a.shape)
    return b | b2 if a.dtype == np.bool_ else _imlincomb_half(b, b2)


def medfilt3(
    a: Any,
    size: Sequence[int] = (3, 3, 3),
    padopt: str = "symmetric",
    *,
    threads: int | None = None,
) -> np.ndarray:
    """MATLAB's ``medfilt3(A, [m n p], padopt)`` with Python's per-axis order: the median of each
    window of ``size = (p, m, n)`` (odd sizes) of a (z, y, x) volume, padded with
    ``'symmetric'`` (default), ``'zeros'`` or ``'replicate'``. A 2D image (MATLAB's matrix, also
    (1, y, x)) goes to medfilt2 with the window ``(m, n)`` and padopt, as MATLAB does."""
    a = np.asarray(a)
    if a.ndim not in (2, 3):
        raise MatlabError("MATLAB:medfilt3:expected3D", "medfilt3: A must be 2-D or 3-D")
    _check_class(a, "medfilt3")
    sz = [int(v) for v in size]
    if len(sz) != 3 or any(v < 1 or v % 2 == 0 for v in sz):
        raise MatlabError("MATLAB:medfilt3:expectedOdd", "medfilt3: [m n p] must be 3 odd sizes")
    p = _padopt(padopt, ("replicate", "zeros", "symmetric"), "medfilt3")
    if a.size == 0:
        return a.copy()
    if a.ndim == 2 or a.shape[0] == 1:
        img = a if a.ndim == 2 else a[0]
        return medfilt2(img, (sz[1], sz[2]), p, threads=threads).reshape(a.shape)
    v = np.ascontiguousarray(a)
    if v.dtype == np.float32 and sz == [3, 3, 3] and p == "symmetric":
        # psf_gen_new's call: the selection network with the padding made in the kernel
        from petakit5d import _ext
        from petakit5d.core.threads import resolve_threads

        out = np.empty(v.shape, dtype=np.float32)
        if _ext.load("morph").median27_symmetric_single(v, out, resolve_threads(threads)):
            return out
    return _order_filter(v, (sz[0], sz[1], sz[2]), (sz[0] * sz[1] * sz[2] + 1) // 2 - 1, p, threads)
