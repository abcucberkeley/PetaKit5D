"""MATLAB's ``imwarp`` for 3D affine transforms (linear, nearest, cubic), value for value.

``imwarp(V, affine3d(T), method, 'FillValues', fill, 'OutputView', R)`` as MATLAB R2024b
computes it, for the output views PetaKit5D uses: ``imref3d(sz, 1, 1, 1)`` (``output_size``) and
the view of a box of it (``output_box``).

Arrays are ``(z, y, x)``. The transform keeps MATLAB's convention (the conventions' exception
4): ``T`` is ``affine3d``'s 4 x 4 matrix, ``[x y z 1] * T`` maps a point of the input to the
output, in world coordinates whose voxel centers are MATLAB's 1-based subscripts (x the
column, y the row, z the plane).

What MATLAB does, and this module with it:

- The composite matrix ``tComp = I2W / T * W2I`` maps output subscripts to input ones (``I2W``
  of the output view, ``W2I`` of the input's default view, the identity), its last column set
  to ``[0 0 0 1]``. MATLAB's ``a / T`` is a triangular solve when ``T`` is triangular, else an LU
  of ``T`` (partial pivoting, multipliers by division) and the transposed solve, each in
  LAPACK's dot-product order; computed here in Python floats in the same order.
- ``linear`` and ``nearest`` run MATLAB's ``warp3d`` built-in, ``cubic`` runs ``interp3``'s
  cubic convolution (griddedInterpolant); the arithmetic per output voxel is in
  ``python/cpp/warp/warp.cpp``. Values never depend on the thread count.
- The class of the result is the input's: ``uint8``, ``int16``, ``uint16``, ``single`` or
  ``double`` (other classes take other paths in MATLAB and are refused).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = ["ImwarpError", "composite_matrix", "imwarp", "mrdivide"]

_CLASSES = {
    np.dtype(np.uint8): "uint8",
    np.dtype(np.int16): "int16",
    np.dtype(np.uint16): "uint16",
    np.dtype(np.float32): "single",
    np.dtype(np.float64): "double",
}
_METHODS = ("nearest", "linear", "cubic", "bilinear", "bicubic")
_METHOD_CODE = {"linear": 0, "nearest": 1, "cubic": 2}


class ImwarpError(ValueError):
    """An error where MATLAB raises one; ``identifier`` is MATLAB's error id."""

    def __init__(self, message: str, identifier: str):
        super().__init__(message)
        self.identifier = identifier


def _method(method: str) -> str:
    """validatestring over imwarp's methods (a unique case-insensitive prefix); 'bilinear' and
    'bicubic' are 'linear' and 'cubic'"""
    if not isinstance(method, str) or not method:
        raise ImwarpError("method must be a nonempty string", "MATLAB:imwarp:invalidType")
    m = method.lower()
    hits = [c for c in _METHODS if c.startswith(m)]
    if m in _METHODS:
        hits = [m]
    if len(hits) != 1:
        raise ImwarpError(
            f"Expected method to match one of these values: {', '.join(_METHODS)}; "
            f"the input, '{method}', "
            + ("matches several of them." if hits else "did not match any of them."),
            "MATLAB:imwarp:unrecognizedStringChoice"
            if not hits
            else "MATLAB:imwarp:ambiguousStringChoice",
        )
    return {"bilinear": "linear", "bicubic": "cubic"}.get(hits[0], hits[0])


def _matrix(tform) -> list[list[float]]:
    """affine3d's checks of T, in its order: size, finite, not singular, last column"""
    a = np.asarray(tform)
    if a.dtype == np.float32:
        raise TypeError(
            "a single precision T is not supported (MATLAB computes its matrix in single)"
        )
    if a.shape != (4, 4):
        raise ImwarpError("Expected T to be of size 4x4.", "MATLAB:affine3d:set:T:incorrectSize")
    a = a.astype(np.float64)
    if not np.all(np.isfinite(a)):
        raise ImwarpError("Expected T to be finite.", "MATLAB:affine3d:set:T:expectedFinite")
    m = a.tolist()
    if _det(m) == 0:
        raise ImwarpError(
            "The transformation matrix is singular.", "images:geotrans:singularTransformationMatrix"
        )
    if [r[3] for r in m] != [0.0, 0.0, 0.0, 1.0]:
        raise ImwarpError(
            "The last column of T must be [0 0 0 1]'.", "images:geotrans:invalidAffineMatrix"
        )
    return m


def _lu(A: list[list[float]]) -> tuple[list[list[float]], list[int]]:
    """LU with partial pivoting (the first largest magnitude), multipliers by division"""
    A = [list(r) for r in A]
    piv = []
    for j in range(4):
        p = j
        for i in range(j + 1, 4):
            if abs(A[i][j]) > abs(A[p][j]):
                p = i
        piv.append(p)
        if p != j:
            A[j], A[p] = A[p], A[j]
        if A[j][j] == 0:
            continue
        for i in range(j + 1, 4):
            A[i][j] = A[i][j] / A[j][j]
        for k in range(j + 1, 4):
            for i in range(j + 1, 4):
                A[i][k] = A[i][k] - A[i][j] * A[j][k]
    return A, piv


def _det(A: list[list[float]]) -> float:
    LU, piv = _lu(A)
    d = 1.0
    for j in range(4):
        d = d * LU[j][j]
        if piv[j] != j:
            d = -d
    return d


def _solve_lu_t(LU: list[list[float]], piv: list[int], b: list[float]) -> list[float]:
    """T' x = b with P T = L U: U' then L' in dot-product form, then the swaps in reverse"""
    x = list(b)
    for i in range(4):
        t = x[i]
        for k in range(i):
            t = t - LU[k][i] * x[k]
        x[i] = t / LU[i][i]
    for i in range(3, -1, -1):
        t = x[i]
        for k in range(i + 1, 4):
            t = t - LU[k][i] * x[k]
        x[i] = t
    for j in range(3, -1, -1):
        p = piv[j]
        if p != j:
            x[j], x[p] = x[p], x[j]
    return x


def _solve_lower_t(L: list[list[float]], b: list[float]) -> list[float]:
    """L' x = b for a lower triangular L (backward, dot-product form)"""
    x = list(b)
    for i in range(3, -1, -1):
        t = x[i]
        for k in range(i + 1, 4):
            t = t - L[k][i] * x[k]
        x[i] = t / L[i][i]
    return x


def _solve_upper_t(U: list[list[float]], b: list[float]) -> list[float]:
    """U' x = b for an upper triangular U (forward, dot-product form)"""
    x = list(b)
    for i in range(4):
        t = x[i]
        for k in range(i):
            t = t - U[k][i] * x[k]
        x[i] = t / U[i][i]
    return x


def mrdivide(a, T) -> list[list[float]]:
    """MATLAB's ``a / T`` for 4 x 4 double matrices, bit for bit: the rows ``x`` of ``x * T = a``.

    As MATLAB R2024b computes it: a triangular solve when ``T`` is lower or upper triangular,
    else an LU of ``T`` (partial pivoting, the first largest magnitude, multipliers by division)
    and the transposed solve, each in LAPACK's dot-product order. ``a`` and ``T``: 4 x 4 array-
    likes (rows); the result is a list of 4 rows of Python floats.
    """
    a = [[float(v) for v in row] for row in a]
    T = [[float(v) for v in row] for row in T]
    if len(a) != 4 or len(T) != 4 or any(len(r) != 4 for r in a + T):
        raise ValueError("mrdivide: 4 x 4 matrices")
    if all(T[i][j] == 0 for i in range(4) for j in range(i + 1, 4)):
        return [_solve_lower_t(T, row) for row in a]
    if all(T[i][j] == 0 for i in range(4) for j in range(i)):
        return [_solve_upper_t(T, row) for row in a]
    LU, piv = _lu(T)
    return [_solve_lu_t(LU, piv, row) for row in a]


_mrdivide = mrdivide  # the earlier private name (deskew_rotate.geometry imports it)


def _view(
    output_size: Sequence[int] | None, output_box: Sequence[int] | None
) -> tuple[tuple[int, int, int], list[float]]:
    """(nz, ny, nx) of the output and the translation row of the output view's
    intrinsic-to-world matrix ([x0, y0, z0]: the box's first voxel, 0 for the full view)"""
    if (output_size is None) == (output_box is None):
        raise ValueError("give exactly one of output_size and output_box")
    if output_size is not None:
        size = tuple(int(s) for s in output_size)
        if len(size) != 3 or any(s < 1 for s in size):
            raise ValueError(f"output_size: three positive sizes (nz, ny, nx), got {output_size}")
        return size, [0.0, 0.0, 0.0]
    b = [int(v) for v in output_box]
    if len(b) != 6 or any(b[i + 3] <= b[i] for i in range(3)):
        raise ValueError(
            f"output_box: (z0, y0, x0, z1, y1, x1) with z1 > z0, ..., got {output_box}"
        )
    size = (b[3] - b[0], b[4] - b[1], b[5] - b[2])
    return size, [float(b[2]), float(b[1]), float(b[0])]


def composite_matrix(
    tform, output_size: Sequence[int] | None = None, output_box: Sequence[int] | None = None
) -> np.ndarray:
    """imwarp's composite matrix tComp (4 x 4, MATLAB's convention): output subscripts to input
    subscripts, as MATLAB computes it for the output view."""
    T = _matrix(tform)
    _, shift = _view(output_size, output_box)
    a = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [*shift, 1.0]]
    X = mrdivide(a, T)
    # times the input's world-to-intrinsic matrix (the identity), as a matrix product
    eye = [[1.0 if r == c else 0.0 for c in range(4)] for r in range(4)]
    t = []
    for r in range(4):
        row = []
        for c in range(4):
            s = 0.0
            for m in range(4):
                s = s + X[r][m] * eye[m][c]
            row.append(s)
        t.append(row)
    for r in range(4):
        t[r][3] = 1.0 if r == 3 else 0.0
    return np.array(t, dtype=np.float64)


def imwarp(
    vol: np.ndarray,
    tform,
    method: str = "linear",
    *,
    output_size: Sequence[int] | None = None,
    output_box: Sequence[int] | None = None,
    fill_value: float = 0,
    num_threads: int | None = None,
) -> np.ndarray:
    """MATLAB's imwarp(vol, affine3d(T), method, 'FillValues', fill, 'OutputView', R).

    vol: (z, y, x) array of class uint8, int16, uint16, single or double (C order; another
    layout is copied). tform: affine3d's T (4 x 4, MATLAB's convention). method: 'linear',
    'nearest' or 'cubic' ('bilinear', 'bicubic', unique prefixes). The output view, exactly
    one of: output_size (nz, ny, nx), MATLAB's
    ``imref3d([ny nx nz], 1, 1, 1)``; output_box (z0, y0, x0, z1, y1, x1), 0-based and
    end-exclusive, that part of the full view (MATLAB's ``imref3d`` with world limits
    ``[x0 + 0.5, x1 + 0.5]``, ...). fill_value: for voxels mapped outside the input, cast like
    vol. num_threads: None or 0 for all CPUs this process may use.

    Returns the (nz, ny, nx) result, of vol's class.
    """
    from petakit5d._ext import load
    from petakit5d.core import arrays, dtypes, threads

    m = _method(method)
    a = np.asarray(vol)
    if a.dtype not in _CLASSES:
        raise TypeError(
            f"imwarp: class {a.dtype} is not supported (uint8, int16, uint16, single, double)"
        )
    # a 2D MATLAB image (one plane) does not agree with a 3D transform
    if a.ndim != 3 or a.shape[0] == 1:
        raise ImwarpError(
            "The input image A must be a 3D array for a 3D geometric transformation.",
            "images:imwarp:tformDoesNotAgreeWithSizeOfA",
        )
    a = arrays.as_volume(a)
    size, _ = _view(output_size, output_box)
    tcomp = composite_matrix(tform, output_size, output_box)
    if m == "cubic" and min(a.shape) < 3:
        if min(a.shape) == 1:
            raise ImwarpError(
                "Interpolation requires at least two sample points in each dimension.",
                "MATLAB:griddedInterpolant:DegenerateGridErrId",
            )
        raise NotImplementedError(
            "imwarp cubic with 2 samples along an axis: MATLAB falls back to griddedInterpolant's "
            "linear method, which this port does not reproduce"
        )
    fill = float(dtypes.cast(np.float64(fill_value), _CLASSES[a.dtype]))
    out = np.empty(size, dtype=a.dtype)
    load("warp").warp_affine(
        a, tcomp, out, _METHOD_CODE[m], fill, threads.resolve_threads(num_threads)
    )
    return out
