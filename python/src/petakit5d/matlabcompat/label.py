"""MATLAB's connected components, exactly (Image Processing Toolbox R2024b, measured with oracle
cases).

Images are NumPy's ``(y, x)`` or ``(z, y, x)``, C-contiguous; a connectivity array follows the
array rule (MATLAB's 3 x 3 x 3 ``rows x cols x planes`` is ``(planes, rows, cols)``; a 3 x 3
array is ``(y, x)``). Label numbers and component lists follow MATLAB's linear order (y
fastest, then x, then z), so they equal MATLAB's:

- ``bwconncomp`` and ``bwlabel`` number the components by their first voxel in MATLAB's order;
  a ``ConnComp``'s pixel lists hold C flat indices (``np.ravel`` order) sorted in MATLAB's
  order.
- ``bwlabeln`` numbers them differently (measured on MATLAB): voxels in MATLAB's order; a
  voxel takes the provisional label of its first labeled earlier neighbor, in the column-major
  order of MATLAB's connectivity matrix, and the sets of its other labeled earlier neighbors
  are united with its set in that order, by size in provisional labels (on a tie the voxel's
  set keeps its root); a voxel without labeled earlier neighbors gets a new provisional label;
  the final labels go to the roots in increasing provisional order.
- ``bwareaopen`` keeps the components of at least ``p`` voxels; ``imfill(bw, 'holes')`` fills
  the background not connected to the border (connectivity 4 or 6 by default).
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np

from petakit5d.core.arrays import size_to_matlab, to_matlab, zyx_to_matlab_index

__all__ = [
    "ConnComp",
    "bwareaopen",
    "bwconncomp",
    "bwlabel",
    "bwlabeln",
    "conndef",
    "connectivity_matrix",
    "imfill",
    "labelmatrix",
]


class MatlabError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _msize(a: np.ndarray) -> tuple[int, ...]:
    """MATLAB's size of a Python (y, x) or (z, y, x) array"""
    return size_to_matlab(a.shape)


def conndef(num_dims: int, type: str) -> np.ndarray:
    """MATLAB's ``conndef(num_dims, 'minimal' | 'maximal')`` (double, 3 x ... x 3; symmetric, so
    the same in both forms)"""
    t = type.lower()
    if t and "minimal".startswith(t):
        conn = np.zeros((3,) * num_dims)
        c = (1,) * num_dims
        conn[c] = 1
        for k in range(num_dims):
            for d in (-1, 1):
                idx = list(c)
                idx[k] += d
                conn[tuple(idx)] = 1
        return conn
    if t and "maximal".startswith(t):
        return np.ones((3,) * num_dims)
    raise MatlabError("MATLAB:conndef:unrecognizedStringChoice", "type must be minimal or maximal")


def _check_conn(conn: Any, func: str) -> Any:
    """iptcheckconn (in either form: its checks are symmetric)"""
    c = np.asarray(conn)
    # a Python int stands for MATLAB's double; NumPy integer classes are refused as in MATLAB
    if c.dtype.kind not in "bf" and not isinstance(conn, int):
        raise MatlabError(f"MATLAB:{func}:invalidType", f"{func}: CONN must be double or logical")
    if c.size == 1:
        if float(c.reshape(-1)[0]) not in (1, 4, 6, 8, 18, 26):
            raise MatlabError("images:validate:badScalarConn", f"{func}: invalid scalar CONN")
        return conn
    if any(s != 3 for s in c.shape):
        raise MatlabError("images:validate:badConnSize", f"{func}: CONN must be 3 x 3 x ... x 3")
    flat = c.reshape(-1)
    if np.any((flat != 0) & (flat != 1)):
        raise MatlabError("images:validate:badConnValue", f"{func}: CONN must contain 0s and 1s")
    if flat[(flat.size - 1) // 2] == 0:
        raise MatlabError("images:validate:badConnCenter", f"{func}: CONN must have its center")
    if not np.array_equal(flat, flat[::-1]):
        raise MatlabError("images:validate:nonsymmetricConn", f"{func}: CONN must be symmetric")
    return conn


def connectivity_matrix(conn: Any) -> np.ndarray:
    """``images.internal.getBinaryConnectivityMatrix`` in Python's form: a bool (3, 3) or
    (3, 3, 3) array (MATLAB's 3 x 3 x 3 matrix as (planes, rows, cols))"""
    c = np.asarray(conn)
    if c.size == 1:
        v = float(c.reshape(-1)[0])
        if v == 1:
            return np.ones((1, 1), dtype=bool)
        if v == 4:
            return conndef(2, "minimal") != 0
        if v == 6:
            return conndef(3, "minimal") != 0
        if v == 8:
            return np.ones((3, 3), dtype=bool)
        if v == 18:
            m = np.ones((3, 3, 3), dtype=bool)
            m[::2, ::2, ::2] = False
            return m
        if v == 26:
            return np.ones((3, 3, 3), dtype=bool)
        raise MatlabError("images:validate:badScalarConn", "invalid scalar connectivity")
    return np.ascontiguousarray(c != 0)


def _offsets(connb: np.ndarray) -> np.ndarray:
    """the neighbors of a Python-form connectivity before its center in C order, rows
    (dz, dy, dx)"""
    if connb.ndim > 3:
        raise NotImplementedError("connectivities of more than 3 dimensions are not supported")
    c3 = connb.reshape((1,) * (3 - connb.ndim) + connb.shape)
    subs = np.argwhere(c3)
    d = subs - (np.array(c3.shape) - 1) // 2
    lin = 9 * d[:, 0] + 3 * d[:, 1] + d[:, 2]
    return np.ascontiguousarray(d[lin < 0], dtype=np.int64)


def _matlab_offsets(connb: np.ndarray) -> np.ndarray:
    """the neighbors of the connectivity before its center in MATLAB's order, rows (dy, dx, dz)
    in the column-major order of MATLAB's matrix (bwlabeln's neighbor order)"""
    if connb.ndim > 3:
        raise NotImplementedError("connectivities of more than 3 dimensions are not supported")
    m = to_matlab(connb)  # (rows, cols[, planes]), Fortran
    m = m.reshape(m.shape + (1,) * (3 - m.ndim), order="F")
    idx = np.flatnonzero(m.ravel(order="F"))
    subs = np.stack(np.unravel_index(idx, m.shape, order="F"), axis=1)
    d = subs - (np.array(m.shape) - 1) // 2
    lin = d[:, 0] + 3 * d[:, 1] + 9 * d[:, 2]
    return np.ascontiguousarray(d[lin < 0], dtype=np.int64)


def _as_bw(bw: Any, func: str) -> np.ndarray:
    a = np.asarray(bw)
    if a.dtype.kind == "c":
        raise MatlabError(f"MATLAB:{func}:expectedReal", f"{func}: BW must be real")
    if a.dtype != np.bool_:
        a = a != 0
    if a.ndim > 3:
        raise NotImplementedError(f"{func}: images of more than 3 dimensions are not supported")
    if a.ndim < 2:
        a = a.reshape((1, -1) if a.ndim == 1 else (1, 1))
    return np.ascontiguousarray(a)


def _vol(a: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(a.reshape((1,) * (3 - a.ndim) + a.shape))


def _label(a: np.ndarray, conn: Any, kind: str, threads: int | None):
    """(labels int32 or int64 as a (z, y, x) volume, number of components); kind 'first'
    (MATLAB's first-voxel order) or 'any' (C order, for sizes and holes)"""
    from petakit5d import _ext
    from petakit5d.core.threads import resolve_threads

    morph = _ext.load("morph")
    off = _offsets(connectivity_matrix(conn))
    v = _vol(a)
    dt, suffix = (np.int32, "int32") if v.size < 2**31 - 1 else (np.int64, "int64")
    lab = np.empty(v.shape, dtype=dt)
    n = getattr(morph, f"label_{kind}_{suffix}")(v, lab, off, resolve_threads(threads))
    return lab, int(n)


def _default_conn(a: np.ndarray) -> int:
    return 8 if len(_msize(a)) == 2 else 26


def bwlabeln(bw: Any, conn: Any = None, *, threads: int | None = None) -> tuple[np.ndarray, int]:
    """MATLAB's ``[L, num] = bwlabeln(BW, conn)``: L double, numbered as bwlabeln numbers (see
    the module docstring). ``conn``: 4, 8, 6, 18, 26 or a connectivity array; default the
    maximal connectivity (8 in 2D, 26 in 3D)."""
    from petakit5d import _ext
    from petakit5d.core.threads import resolve_threads

    a = _as_bw(bw, "bwlabeln")
    conn = _default_conn(a) if conn is None else _check_conn(conn, "bwlabeln")
    off = _matlab_offsets(connectivity_matrix(conn))
    v = _vol(a)
    out = np.empty(v.shape, dtype=np.float64)
    morph = _ext.load("morph")
    f = morph.bwlabeln_int32 if v.size < 2**31 - 1 else morph.bwlabeln_int64
    n = f(v, out, off, resolve_threads(threads))
    return out.reshape(a.shape), int(n)


def bwlabel(bw: Any, n: Any = 8, *, threads: int | None = None) -> tuple[np.ndarray, int]:
    """MATLAB's ``[L, num] = bwlabel(BW, n)`` of a (y, x) image (n = 4 or 8): L double,
    components in the order of their first pixel in MATLAB's order"""
    a = _as_bw(bw, "bwlabel")
    if a.ndim != 2:
        raise MatlabError("MATLAB:bwlabel:expected2D", "bwlabel: BW must be 2-D")
    if float(n) not in (4, 8):
        raise MatlabError("images:bwlabel:badConn", "bwlabel: N must be 4 or 8")
    lab, num = _label(a, int(n), "first", threads)
    return lab.reshape(a.shape).astype(np.float64), num


def _matlab_index(idx: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    """C flat indices of an array of ``shape`` as MATLAB's 0-based linear indices"""
    sub = np.unravel_index(np.asarray(idx, dtype=np.int64), shape)
    return np.asarray(zyx_to_matlab_index(sub, shape), dtype=np.int64)


@dataclasses.dataclass
class ConnComp:
    """MATLAB's ``bwconncomp`` structure: ``image_size`` the image's shape, ``pixel_idx_list``
    per component its C flat indices (``np.ravel`` order) sorted in MATLAB's linear order,
    components in the order of their first voxel in that order"""

    connectivity: Any
    image_size: tuple[int, ...]
    num_objects: int
    pixel_idx_list: list[np.ndarray]

    def to_matlab(self) -> dict:
        """the structure as MATLAB holds it (MATLAB's size, 1-based linear indices as double
        columns, a connectivity array in MATLAB's form)"""
        c = np.asarray(self.connectivity)
        conn = float(c.reshape(-1)[0]) if c.size == 1 else to_matlab(c).astype(np.float64)
        msize = size_to_matlab(self.image_size)
        return {
            "Connectivity": conn,
            "ImageSize": [float(s) for s in msize],
            "NumObjects": float(self.num_objects),
            "PixelIdxList": [
                (_matlab_index(p, self.image_size) + 1).astype(np.float64)
                for p in self.pixel_idx_list
            ],
        }


def bwconncomp(bw: Any, conn: Any = None, *, threads: int | None = None) -> ConnComp:
    """MATLAB's ``CC = bwconncomp(BW, conn)``; default connectivity 8 (2D) or 26 (3D)"""
    a = _as_bw(bw, "bwconncomp")
    if conn is None:
        conn = _default_conn(a)
    else:
        _check_conn(conn, "bwconncomp")
        c = np.asarray(conn)
        if c.shape == (3, 3) and np.array_equal(c != 0, conndef(2, "minimal") != 0):
            conn = 4
        elif c.shape == (3, 3) and np.all(c != 0):
            conn = 8
    lab, n = _label(a, conn, "first", threads)
    return ConnComp(conn, tuple(a.shape), n, _pixel_lists(lab, n))


def _pixel_lists(lab: np.ndarray, n: int) -> list[np.ndarray]:
    from petakit5d import _ext

    morph = _ext.load("morph")
    group = morph.group_by_label_int32 if lab.dtype == np.int32 else morph.group_by_label_int64
    idx, starts = group(lab, n)
    return [idx[starts[k] : starts[k + 1]] for k in range(n)]


def labelmatrix(cc: ConnComp) -> np.ndarray:
    """MATLAB's ``labelmatrix(CC)``: uint8, uint16, uint32 or double, the smallest class that
    holds NumObjects"""
    n = cc.num_objects
    if n <= 255:
        dt = np.uint8
    elif n <= 65535:
        dt = np.uint16
    elif n <= 2**32 - 1:
        dt = np.uint32
    else:
        dt = np.float64
    out = np.zeros(int(np.prod(cc.image_size)), dtype=dt)
    for k, p in enumerate(cc.pixel_idx_list):
        out[p] = k + 1
    return out.reshape(cc.image_size)


def bwareaopen(bw: Any, p: Any, conn: Any = None, *, threads: int | None = None) -> np.ndarray:
    """MATLAB's ``bwareaopen(BW, P, conn)``: the components with at least P voxels (default
    connectivity: the maximal one, 8 in 2D, 26 in 3D)"""
    a = _as_bw(bw, "bwareaopen")
    pv = float(p)
    if not (pv >= 0 and pv == int(pv)):
        raise MatlabError(
            "MATLAB:bwareaopen:expectedNonnegative", "bwareaopen: P must be a nonnegative integer"
        )
    conn = _default_conn(a) if conn is None else _check_conn(conn, "bwareaopen")
    lab, n = _label(a, conn, "any", threads)
    sizes = np.bincount(lab.reshape(-1), minlength=n + 1)
    keep = sizes >= pv
    keep[0] = False
    return keep[lab].reshape(a.shape)


def imfill(bw: Any, option: Any = "holes", conn: Any = None, *, threads: int | None = None):
    """MATLAB's ``imfill(BW, 'holes')`` and ``imfill(BW, conn, 'holes')`` of a logical image:
    the background pixels not connected (connectivity ``conn``, default 4 in 2D, 6 in 3D) to
    the image border become true. The option matches as MATLAB's validatestring does
    (``'hole'`` is ``'holes'``)."""
    if not isinstance(option, str):
        option, conn = conn, option  # imfill(BW, conn, 'holes')
    if not isinstance(option, str) or not option or not "holes".startswith(option.lower()):
        raise NotImplementedError("imfill: only the 'holes' option is supported")
    a = np.asarray(bw)
    if a.dtype != np.bool_:
        raise NotImplementedError("imfill: only logical (bool) images are supported")
    a = _as_bw(a, "imfill")
    conn = conndef(len(_msize(a)), "minimal") if conn is None else _check_conn(conn, "imfill")
    lab, n = _label(~a, conn, "any", threads)
    lab = lab.reshape(a.shape)
    border = np.zeros(n + 1, dtype=bool)
    # MATLAB pads (and so borders) its own dimensions: a (1, y, x) volume is a matrix there
    for d in range(a.ndim - len(_msize(a)), a.ndim):
        for end in (0, a.shape[d] - 1):
            sl = [slice(None)] * a.ndim
            sl[d] = end
            border[lab[tuple(sl)].reshape(-1)] = True
    border[0] = True
    return a | ~border[lab]
