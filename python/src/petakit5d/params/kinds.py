"""Conversion kinds: how a parameter's MATLAB value becomes its Python value and back.

The snake_case API takes Python values (CONVENTION_BRIEF.md: boxes ``(z0, y0, x0, z1, y1, x1)``
0-based and end-exclusive, sizes and per-axis vectors ``(z, y, x)``, NumPy axes, ``'zyx'`` axis
strings, ``(z, y, x)`` arrays). The MATLAB-name layer, the command line, configs and the saved
parameter files hold MATLAB's values. Each parameter names its kind (``Param(convert=...)``):

======================  =====================================  ==================================
kind                    MATLAB                                 Python
======================  =====================================  ==================================
``same``                checked: no conversion
``bbox``                ``[ymin xmin zmin ymax xmax zmax]``    ``(z0, y0, x0, z1, y1, x1)``
                        (2D ``[ymin xmin ymax xmax]``)          (2D ``(y0, x0, y1, x1)``); Inf: None
``bboxes``              rows of boxes (``n x 6`` or ``n x 4``)  tuple of boxes
``size_yxz``            ``[y x z]``, also 6 or 9 values         ``(z, y, x)`` per group of three;
                        (groups of three); 2D ``[y x]``         2D ``(y, x)``; scalars pass
``per_axis_xyz``        ``[x y z]`` (2D ``[x y]``)              ``(z, y, x)`` (2D ``(y, x)``)
``resample_yxz``        ``[y x z]``, or ``[xy z]``              ``(z, y, x)``, or ``(z, xy)``
``resample_xyz``        ``[x y z]``, or ``[xy z]``              ``(z, y, x)``, or ``(z, xy)``
``flags_yxz``           ``[y x z]`` flags (MIP ``[0 0 1]``)     ``(z, y, x)`` (``(1, 0, 0)``)
``index1``              1-based indices                         0-based indices
``point1_yxz``          1-based ``[y x z]`` coordinates         0-based ``(z, y, x)``
``range1``              ``[first last]``                        ``(start, stop)``; Inf: None
``dim``                 dimension numbers (1 y, 2 x, 3 z)       NumPy axes (-2, -1, -3)
``axis_string``         ``'yxz'`` (names the array's axes)      ``'zyx'``
``data_order``          ``'y,x,z'``                             ``'z,y,x'``
``volume``              ``(y, x, z)`` Fortran image array       ``(z, y, x)`` C array
======================  =====================================  ==================================

Empty values (``[]``, ``''``) are None (or ``""``) on both sides. Anything else is a pair of
functions ``convert=(from_matlab, to_matlab)``. NumPy and the core helpers are imported only when
a conversion runs, so importing the registry stays light.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

__all__ = ["KINDS", "Kind", "get_kind", "looks_converted"]


class Kind:
    """A named pair of conversions (both take and give plain Python values)."""

    __slots__ = ("name", "from_matlab", "to_matlab", "doc")

    def __init__(self, name: str, from_matlab: Callable, to_matlab: Callable, doc: str = ""):
        self.name = name
        self.from_matlab = from_matlab
        self.to_matlab = to_matlab
        self.doc = doc

    def __repr__(self) -> str:
        return f"<conversion kind {self.name}>"


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    if type(v).__module__ == "numpy":
        return v.size == 0
    try:
        return len(v) == 0
    except TypeError:
        return False


def _empty_out(v: Any) -> Any:
    return "" if isinstance(v, str) else None


def _skip_empty(f: Callable) -> Callable:
    def g(v: Any) -> Any:
        if _empty(v):
            return _empty_out(v)
        return f(v)

    return g


def _flat(v: Any) -> list:
    import numpy as np

    a = np.asarray(v)
    if a.ndim > 1 and sum(n != 1 for n in a.shape) > 1:
        raise ValueError(f"a per-axis value is a scalar or a vector, got shape {a.shape}")
    return a.reshape(-1).tolist()


def _is_scalar(v: Any) -> bool:
    import numpy as np

    return np.ndim(v) == 0


def _same_type(orig: Any, vals: list) -> Any:
    if type(orig).__module__ == "numpy":
        import numpy as np

        return np.asarray(vals, dtype=orig.dtype)  # the class stays (uint16 stays uint16)
    return tuple(vals)


def _triplets(order: tuple[int, int, int], order2: tuple[int, int] | None) -> Callable:
    def f(v: Any) -> Any:
        if _is_scalar(v):
            return v
        x = _flat(v)
        if len(x) == 2 and order2 is not None:
            return _same_type(v, [x[i] for i in order2])
        if len(x) % 3 or not x:
            raise ValueError(f"expected 2, 3, 6 or 9 per-axis values, got {len(x)}")
        out = []
        for k in range(0, len(x), 3):
            out.extend(x[k + i] for i in order)
        return _same_type(v, out)

    return _skip_empty(f)


# [y x z] <-> (z, y, x); [x y z] <-> (z, y, x)
_YXZ_IN, _YXZ_OUT = (2, 0, 1), (1, 2, 0)
_XYZ = (2, 1, 0)


def _box_in(v: Any) -> Any:
    from petakit5d.core.bbox import box_from_matlab

    return box_from_matlab(v)


def _box_out(v: Any) -> Any:
    from petakit5d.core.bbox import box_to_matlab

    return box_to_matlab(v)


def _rows(v: Any) -> list[list]:
    import numpy as np

    a = np.asarray(v, dtype=object)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.ndim != 2:
        raise ValueError(f"expected rows of boxes, got shape {a.shape}")
    return [list(r) for r in a.tolist()]


def _boxes_in(v: Any) -> Any:
    return tuple(_box_in(r) for r in _rows(v))


def _boxes_out(v: Any) -> Any:
    return tuple(_box_out(r) for r in _rows(v))


def _map_numbers(f: Callable[[Any], Any]) -> Callable:
    def g(v: Any) -> Any:
        if _is_scalar(v):
            return f(v)
        if type(v).__module__ == "numpy":
            import numpy as np

            return np.vectorize(f, otypes=[np.int64])(v)
        return tuple(g(e) for e in v)

    return _skip_empty(g)


def _int(v: Any) -> int:
    f = float(v)
    if f != int(f):
        raise ValueError(f"an index must be an integer, got {v!r}")
    return int(f)


def _point_in(v: Any) -> Any:
    return _map_numbers(lambda e: _int(e) - 1)(_triplets(_YXZ_IN, (0, 1))(v))


def _point_out(v: Any) -> Any:
    return _triplets(_YXZ_OUT, (0, 1))(_map_numbers(lambda e: _int(e) + 1)(v))


def _range_in(v: Any) -> Any:
    from petakit5d.core.bbox import range_from_matlab

    return range_from_matlab(v)


def _range_out(v: Any) -> Any:
    from petakit5d.core.bbox import range_to_matlab

    return range_to_matlab(v)


def _dim_in(v: Any) -> int:
    from petakit5d.core.axes import dim_to_axis

    return dim_to_axis(_int(v))


def _dim_out(v: Any) -> int:
    from petakit5d.core.axes import axis_to_dim

    a = _int(v)
    # a non-negative axis is read for a 3D volume
    return axis_to_dim(a, max(3, -a))


def _axis_in(v: Any) -> Any:
    from petakit5d.core.axes import axis_string_from_matlab

    return axis_string_from_matlab(v)


def _axis_out(v: Any) -> Any:
    from petakit5d.core.axes import axis_string_to_matlab

    return axis_string_to_matlab(v)


def _vol_in(v: Any) -> Any:
    from petakit5d.core.arrays import from_matlab

    return from_matlab(v)


def _vol_out(v: Any) -> Any:
    from petakit5d.core.arrays import to_matlab

    return to_matlab(v)


def _resample(order: tuple[int, int, int]) -> Callable:
    # a 2-value factor is [xy z] (MATLAB completes it to [xy xy z]): (z, xy) in Python
    return _triplets(order, (1, 0))


def _same(v: Any) -> Any:
    return v


KINDS: dict[str, Kind] = {
    k.name: k
    for k in (
        Kind("same", _same, _same, "checked: no conversion"),
        Kind("bbox", _skip_empty(_box_in), _skip_empty(_box_out), "a box"),
        Kind("bboxes", _skip_empty(_boxes_in), _skip_empty(_boxes_out), "rows of boxes"),
        Kind("size_yxz", _triplets(_YXZ_IN, (0, 1)), _triplets(_YXZ_OUT, (0, 1)), "[y x z]"),
        Kind("per_axis_xyz", _triplets(_XYZ, (1, 0)), _triplets(_XYZ, (1, 0)), "[x y z]"),
        Kind("resample_yxz", _resample(_YXZ_IN), _resample(_YXZ_OUT), "[y x z] or [xy z]"),
        Kind("resample_xyz", _resample(_XYZ), _resample(_XYZ), "[x y z] or [xy z]"),
        Kind("flags_yxz", _triplets(_YXZ_IN, (0, 1)), _triplets(_YXZ_OUT, (0, 1)), "flags"),
        Kind(
            "index1",
            _map_numbers(lambda e: _int(e) - 1),
            _map_numbers(lambda e: _int(e) + 1),
            "1-based index",
        ),  # fmt: skip
        Kind("point1_yxz", _point_in, _point_out, "1-based [y x z] point"),
        Kind("range1", _skip_empty(_range_in), _skip_empty(_range_out), "[first last]"),
        Kind("dim", _map_numbers(_dim_in), _map_numbers(_dim_out), "dimension number"),
        Kind("axis_string", _skip_empty(_axis_in), _skip_empty(_axis_out), "'yxz' axis names"),
        Kind("data_order", _skip_empty(_axis_in), _skip_empty(_axis_out), "'y,x,z' data order"),
        Kind("volume", _skip_empty(_vol_in), _skip_empty(_vol_out), "image array"),
    )
}


def get_kind(convert: str | tuple[Callable, Callable] | None) -> Kind | None:
    """The Kind of ``Param(convert=...)``: a name, a ``(from_matlab, to_matlab)`` pair, or None."""
    if convert is None:
        return None
    if isinstance(convert, str):
        try:
            return KINDS[convert]
        except KeyError:
            raise ValueError(
                f"unknown conversion kind {convert!r}; one of {sorted(KINDS)}"
            ) from None
    if isinstance(convert, tuple) and len(convert) == 2 and all(callable(f) for f in convert):
        name = getattr(convert[0], "__name__", "custom")
        return Kind(f"custom:{name}", convert[0], convert[1], "own conversion")
    raise ValueError("convert is a kind name or a (from_matlab, to_matlab) pair")


# MATLAB parameter names that look like boxes, sizes, per-axis vectors, indices, axes or image
# arrays: a registered one needs a kind (convert="same" when it is checked and needs none)
_LOOKS = (
    r"bbox|bboxes|box$|size$|sizes$|order$|resample|downsample|axisweight|^axis$|mipaxis"
    r"|maxzproj|stitchmip|offset$|coords?$|coord$|^range$|inds?$|indices$|idx$|^dims?$|pool"
    r"|border|padsize|padding|^frame$|^vol$|^data$|^img$|^im$|^i$|^psf$|^psf_[fb]$|^background$|rawdata|lsimage"
)


# scalars that name their axes keep their names and values (exception 1)
_NAMED_SCALARS = r"pixelsize$|^dz|skewangle|^(xy|z)maxoffset$"


def looks_converted(matlab_name: str) -> bool:
    """Whether a MATLAB parameter name suggests a value that may need a conversion kind."""
    import re

    if re.search(_NAMED_SCALARS, matlab_name, re.I):
        return False
    return bool(re.search(_LOOKS, matlab_name, re.I))
