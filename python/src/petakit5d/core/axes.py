"""Per-axis values and axis names: ``(z, y, x)`` in the Python API.

Every per-axis vector of the Python API follows ``(z, y, x)``: block, batch and chunk sizes, pool
sizes, borders and paddings, resample factors, voxel sizes given as a vector, shifts in voxels,
per-axis flags. MATLAB's are ``[y x z]`` (most parameters) or ``[x y z]`` (some, by name):
``from_yxz``, ``to_yxz``, ``from_xyz``, ``to_xyz`` convert. A scalar applies to all axes and passes
through unchanged; a 2-value vector is 2D (``[y x]`` and ``(y, x)``, ``[x y]`` and ``(y, x)``).

Exceptions (by name, unchanged): ``xy_pixel_size``, ``dz``, ``skew_angle``, ``xy_max_offset``,
``z_max_offset``; stage coordinates and ``axis_order`` strings (``'x,y,z'``, ``'-x,y,z'``).

Axes: MATLAB's dimension 1 (y) is NumPy axis -2, dimension 2 (x) axis -1, dimension 3 (z) axis
-3, dimension k > 3 axis -k. Axis strings that name an array's axes (``data_order``,
``input_axis_order``, ``output_axis_order``) name them in the array's order, so MATLAB's
``'y,x,z'`` and ``'yxz'`` are our ``'z,y,x'`` and ``'zyx'``.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "axis_permutation",
    "axis_string_from_matlab",
    "axis_string_to_matlab",
    "axis_to_dim",
    "dim_to_axis",
    "from_xyz",
    "from_yxz",
    "to_xyz",
    "to_yxz",
]


def _vec(v):
    """A per-axis value as (list, kind, was_array): kind 'scalar', or the vector's length."""
    if np.ndim(v) == 0:
        return None
    a = np.asarray(v)
    if a.ndim > 1 and sum(n != 1 for n in a.shape) > 1:
        raise ValueError(f"a per-axis value is a scalar or a vector, got shape {a.shape}")
    return a.reshape(-1)


def _back(orig, vals: np.ndarray):
    if isinstance(orig, np.ndarray):
        return vals
    if isinstance(orig, list):
        return vals.tolist()
    return tuple(vals.tolist())


def _reorder(v, order3, order2):
    flat = _vec(v)
    if flat is None:
        return v
    if flat.size == 3:
        out = flat[list(order3)]
    elif flat.size == 2:
        out = flat[list(order2)]
    else:
        raise ValueError(f"a per-axis vector has 2 or 3 values, got {flat.size}")
    return _back(v, out)


def from_yxz(v):
    """MATLAB's ``[y x z]`` (or ``[y x]``) as ``(z, y, x)`` (or ``(y, x)``); scalars pass."""
    return _reorder(v, (2, 0, 1), (0, 1))


def to_yxz(v):
    """``(z, y, x)`` (or ``(y, x)``) as MATLAB's ``[y x z]`` (or ``[y x]``); scalars pass."""
    return _reorder(v, (1, 2, 0), (0, 1))


def from_xyz(v):
    """``[x y z]`` (or ``[x y]``) as ``(z, y, x)`` (or ``(y, x)``); scalars pass."""
    return _reorder(v, (2, 1, 0), (1, 0))


def to_xyz(v):
    """``(z, y, x)`` (or ``(y, x)``) as ``[x y z]`` (or ``[x y]``); scalars pass."""
    return _reorder(v, (2, 1, 0), (1, 0))


def dim_to_axis(dim: int, ndim: int | None = None) -> int:
    """The NumPy axis of MATLAB's dimension ``dim`` (1: y, 2: x, 3: z, k > 3: the k-th): a negative
    axis, or with ``ndim`` the non-negative one of an array of that many dimensions."""
    dim = int(dim)
    if dim < 1:
        raise ValueError(f"MATLAB dimensions start at 1, got {dim}")
    ax = {1: -2, 2: -1}.get(dim, -dim)
    if ndim is None:
        return ax
    if -ax > ndim:
        raise ValueError(f"dimension {dim} does not exist in an array of {ndim} dimensions")
    return ax + ndim


def axis_to_dim(axis: int, ndim: int) -> int:
    """MATLAB's dimension of NumPy axis ``axis`` of an array of ``ndim`` dimensions."""
    axis, ndim = int(axis), int(ndim)
    if not -ndim <= axis < ndim:
        raise ValueError(f"axis {axis} out of range for {ndim} dimensions")
    neg = axis - ndim if axis >= 0 else axis
    return {-2: 1, -1: 2}.get(neg, -neg)


def _letters(s: str) -> tuple[list[str], str]:
    sep = "," if "," in s else ""
    parts = s.split(",") if sep else list(s)
    return [p.strip() for p in parts], sep


def axis_string_from_matlab(s: str) -> str:
    """MATLAB's axis string of an array (``'y,x,z'``, ``'yxz'``, ``'x,y,z'``) as ours (``'z,y,x'``,
    ``'zyx'``, ``'z,x,y'``): the names of MATLAB's dimensions 3, 1, 2 (separators kept)."""
    p, sep = _letters(s)
    if len(p) != 3:
        raise ValueError(f"an axis string names 3 axes, got {s!r}")
    return sep.join([p[2], p[0], p[1]])


def axis_string_to_matlab(s: str) -> str:
    """Our axis string (``'z,y,x'``, ``'zyx'``) as MATLAB's (``'y,x,z'``, ``'yxz'``)."""
    p, sep = _letters(s)
    if len(p) != 3:
        raise ValueError(f"an axis string names 3 axes, got {s!r}")
    return sep.join([p[1], p[2], p[0]])


def axis_permutation(input_axis_order: str, output_axis_order: str = "zyx") -> tuple[int, ...]:
    """The ``np.transpose`` axes that turn an array whose axes are named ``input_axis_order`` into
    one whose axes are named ``output_axis_order`` (our strings: ``'zyx'`` is the convention).

    ``np.transpose(vol, axis_permutation('zyx', 'xzy'))`` stores ``(z, y, x)`` data as
    ``(x, z, y)``; ``axis_permutation('xzy', 'zyx')`` brings it back.
    """
    src = "".join(_letters(input_axis_order.lower())[0])
    dst = "".join(_letters(output_axis_order.lower())[0])
    if sorted(src) != ["x", "y", "z"] or sorted(dst) != ["x", "y", "z"]:
        raise ValueError(
            "axis orders must name x, y and z once each, got "
            f"{input_axis_order!r} and {output_axis_order!r}"
        )
    return tuple(src.index(c) for c in dst)
