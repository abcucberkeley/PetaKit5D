"""MATLAB's ``imresize`` and ``imresize3`` (R2024b, Linux x86-64) bit for bit on the port's
(z, y, x) C arrays, and the lab's ``imresize3_fast`` and ``imresize3_average``.

The built-ins are written from MathWorks's documentation, the published method (separable
resampling with Keys' cubic kernel, a = -1/2, the triangle and box kernels, the kernel widened by
1/scale when shrinking with antialiasing) and black-box probes of the built-ins and of the
internals ``imresize3_fast`` calls (``contributions``, ``imresizemex``, the kernels), on inputs
and outputs only (task N01; no toolbox source was read). What was measured:

- **Output size**: the ceiling of a size given, or ``ceil(scale * size)`` for a scale. With NaN
  in the size (MATLAB keeps the aspect ratio: one NaN for imresize, two for imresize3) every
  dimension takes the scale of the given one and its size ``ceil(scale * size)``. The weights
  use ``out / in`` for a size and the given scale for a scale (``imresize(A, 0.5)`` of 11 rows
  makes 6 rows with the weights of 0.5, not of 6/11).
- **Passes**: one per dimension (MATLAB's y, x, z: axes -2, -1, -3), in ascending order of
  scale, y first among equal ones (a stable sort).
- **Weights** (MATLAB's ``contributions``): output ``x = 1..m`` sits at
  ``u = x/scale + 0.5*(1 - 1/scale)``; shrinking with antialiasing the kernel is
  ``h(x) = scale*k(scale*x)`` and its width ``w/scale``; the taps are
  ``floor(u - w/2) + (0 : ceil(w) + 1)``, their weights ``h(u - tap)`` divided by the row's sum
  (added from the left); taps past the ends are mirrored (``[1..n, n..1]``, periodic); columns
  with no nonzero weight are dropped. The cubic kernel is
  ``(1.5a^3 - 2.5a^2 + 1)(a <= 1) + (-0.5a^3 + 2.5a^2 - 4a + 2)(1 < a <= 2)`` of ``a = |x|``,
  ``a^2`` as ``a*a`` and ``a^3`` MATLAB's power (fdlibm's pow, which is not ``a*a*a``).
- **A pass**: if every output has a single nonzero weight (nearest; box and nearest when
  enlarging; a scale of 1; linear and cubic at odd integer shrink factors without
  antialiasing), the copy of that input value (``-0`` and NaN payloads kept). Otherwise each
  output is a sum in double from +0 over every tap in order, a weight of 0 included (Inf there
  gives NaN), then rounded to the class: single to nearest, integers ``floor(s + 0.5)``
  (``ceil(s - 0.5)`` below 0) saturated, so 0.49999999999999994 gives 1. Integer data is
  rounded after every pass; single data stays single between passes. A NaN sum keeps the
  first NaN in tap order (the sum's first operand keeps its NaN): a NaN input's payload,
  quieted, or MATLAB's NaN for Inf * 0 and Inf - Inf, on every processor.
- ``imresize`` resizes MATLAB's dimensions 1 and 2 (the last two axes) of an array of any number
  of dimensions; ``imresize3`` takes a 3-D array only (one plane is a MATLAB error,
  ``MATLAB:badsubscript``, as MATLAB indexes its third dimension).
- Complex data: the real and imaginary parts separately. Classes: uint8, int8, uint16, int16,
  uint32, int32, single, double.

Not ported: the lanczos kernels, logical data (imresize3's rule for it was not found by the
probes), imresize's colormap forms and GPU arrays (``NotImplementedError``).

The passes run in C++ (``petakit5d._ext.resize``): each output value is one thread's sum in a
fixed order, so results do not depend on the thread count. ``imresize3_fast`` (the lab's
``imresize3_fast.m``) uses the same weights and passes, without the passes of dimensions that
keep their size and, with ``out_bbox``, only on the input region the output region needs.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from petakit5d.core import dtypes
from petakit5d.core.threads import resolve_threads
from petakit5d.params import Param, register

__all__ = [
    "MatlabError",
    "contributions",
    "cubic",
    "imresize",
    "imresize3",
    "imresize3_average",
    "imresize3_fast",
    "kernel_values",
    "resize_axis",
]

_ext_state: dict = {"mod": None}

# kernel name and width of each method
_KERNELS = {"box": ("box", 1.0), "triangle": ("triangle", 2.0), "cubic": ("cubic", 4.0)}
_METHODS3 = {
    "nearest": "box",
    "linear": "triangle",
    "cubic": "cubic",
    "box": "box",
    "triangle": "triangle",
    "lanczos2": None,
    "lanczos3": None,
}
_METHODS2 = {
    "nearest": "box",
    "bilinear": "triangle",
    "bicubic": "cubic",
    "box": "box",
    "triangle": "triangle",
    "cubic": "cubic",
    "lanczos2": None,
    "lanczos3": None,
}
_RESIZE_CLASSES = ("uint8", "int8", "uint16", "int16", "uint32", "int32", "single", "double")


class MatlabError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _ext():
    if _ext_state["mod"] is None:
        from petakit5d._ext import load

        _ext_state["mod"] = load("resize")
    return _ext_state["mod"]


# ---------------------------------------------------------------------------------------------
# kernels and weights
# ---------------------------------------------------------------------------------------------


def cubic(x) -> np.ndarray:
    """MATLAB's cubic kernel (Keys, a = -1/2) of ``x`` (double), with fdlibm's pow for |x|^3"""
    x = np.ascontiguousarray(x, dtype=np.float64)
    out = np.empty_like(x)
    _ext().cubic(x.reshape(-1), out.reshape(-1))
    return out


def kernel_values(kernel: str, x) -> np.ndarray:
    """The kernel ``'cubic'``, ``'triangle'`` or ``'box'`` at ``x`` (double), as MATLAB's
    resize kernels compute it"""
    x = np.asarray(x, dtype=np.float64)
    if kernel == "cubic":
        return cubic(x)
    if kernel == "triangle":
        return (x + 1) * ((-1 <= x) & (x < 0)) + (1 - x) * ((0 <= x) & (x <= 1))
    if kernel == "box":
        return ((-0.5 <= x) & (x < 0.5)).astype(np.float64)
    raise ValueError(f"unknown kernel {kernel!r}")


def contributions(
    in_length: int,
    out_length: int,
    scale: float,
    kernel: str = "cubic",
    kernel_width: float | None = None,
    antialiasing: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """The weights and input indices of each output of one pass, as MATLAB's
    ``contributions(in_length, out_length, scale, kernel, kernel_width, antialiasing)`` (module
    docstring): ``(weights, indices)``, both ``out_length x P``, indices 0-based."""
    n = int(in_length)
    m = int(out_length)
    kw = float(_KERNELS[kernel][1] if kernel_width is None else kernel_width)
    scale = float(scale)
    widen = scale < 1 and bool(antialiasing)
    if widen:
        kw = kw / scale
    x = np.arange(1, m + 1, dtype=np.float64)
    u = x / scale + 0.5 * (1 - 1 / scale)
    left = np.floor(u - kw / 2)
    p = int(math.ceil(kw)) + 2
    ind = left[:, None] + np.arange(p, dtype=np.float64)[None, :]
    d = u[:, None] - ind
    w = scale * kernel_values(kernel, scale * d) if widen else kernel_values(kernel, d)
    s = np.zeros(m)
    for k in range(p):
        s = s + w[:, k]
    w = w / s[:, None]
    aux = np.concatenate([np.arange(n), np.arange(n - 1, -1, -1)]).astype(np.int64)
    ind0 = aux[np.mod(ind.astype(np.int64) - 1, 2 * n)]
    keep = np.any(w != 0, axis=0)
    return np.ascontiguousarray(w[:, keep]), np.ascontiguousarray(ind0[:, keep])


def _one_hot(w: np.ndarray) -> bool:
    return bool(w.size) and bool(np.all(np.count_nonzero(w, axis=1) == 1))


# ---------------------------------------------------------------------------------------------
# one pass
# ---------------------------------------------------------------------------------------------


def _class(a: np.ndarray) -> str:
    if a.dtype == np.bool_:
        raise NotImplementedError(
            "resizing logical data is not ported (MATLAB's rule not measured)"
        )
    try:
        cls = dtypes.matlab_class(a.dtype)
    except (TypeError, KeyError, ValueError):
        cls = str(a.dtype)
    if cls not in _RESIZE_CLASSES:
        raise TypeError(
            f"resizing {cls} data is not supported (MATLAB's classes: {_RESIZE_CLASSES})"
        )
    return cls


def _layout(shape: tuple[int, ...], axis: int) -> tuple[int, int, int]:
    outer = int(np.prod(shape[:axis], dtype=np.int64))
    inner = int(np.prod(shape[axis + 1 :], dtype=np.int64))
    return outer, int(shape[axis]), inner


def _target(a: np.ndarray, axis: int, n_out: int, out: np.ndarray | None) -> np.ndarray:
    shape = a.shape[:axis] + (n_out,) + a.shape[axis + 1 :]
    if out is None:
        return np.empty(shape, dtype=a.dtype)
    if out.shape != shape or out.dtype != a.dtype or not out.flags.c_contiguous:
        raise ValueError(f"out must be a C-contiguous {a.dtype} array of shape {shape}")
    if np.shares_memory(out, a):
        raise ValueError("out must not overlap the input")
    return out


def _gather(
    a: np.ndarray, axis: int, idx: np.ndarray, nt: int, out: np.ndarray | None
) -> np.ndarray:
    n_out = int(idx.size)
    dst = _target(a, axis, n_out, out)
    outer, n_in, inner = _layout(a.shape, axis)
    size = a.dtype.itemsize
    if size not in (1, 2, 4, 8):
        raise TypeError(f"no index copy for {a.dtype}")
    view = np.dtype(f"u{size}")
    fn = getattr(_ext(), f"gather_{size}")
    fn(
        a.view(view),
        dst.view(view),
        outer,
        n_in,
        inner,
        n_out,
        np.ascontiguousarray(idx, np.int64),
        nt,
    )
    return dst


def _weighted(
    a: np.ndarray, axis: int, w: np.ndarray, ind: np.ndarray, nt: int, out: np.ndarray | None
) -> np.ndarray:
    n_out = int(w.shape[0])
    dst = _target(a, axis, n_out, out)
    outer, n_in, inner = _layout(a.shape, axis)
    fn = getattr(_ext(), f"resize_{_class(a)}")
    fn(a, dst, outer, n_in, inner, n_out, w, ind, nt)
    return dst


def _pass(
    a: np.ndarray,
    axis: int,
    w: np.ndarray,
    ind: np.ndarray,
    nt: int,
    out: np.ndarray | None = None,
    *,
    copy_one_hot: bool = True,
) -> np.ndarray:
    """One pass along ``axis`` of a C array: the index copy where the weights are one-hot (and
    ``copy_one_hot``), else the weighted sums of MATLAB's resize"""
    axis = axis % a.ndim
    if np.iscomplexobj(a):
        re = _pass(np.ascontiguousarray(a.real), axis, w, ind, nt, copy_one_hot=copy_one_hot)
        im = _pass(np.ascontiguousarray(a.imag), axis, w, ind, nt, copy_one_hot=copy_one_hot)
        dst = _target(a, axis, int(w.shape[0]), out)
        dst.real = re
        dst.imag = im
        return dst
    _class(a)
    if copy_one_hot and _one_hot(w):
        idx = ind[np.arange(w.shape[0]), np.argmax(w != 0, axis=1)]
        if idx.size == a.shape[axis] and np.array_equal(idx, np.arange(idx.size)):
            if out is None:
                return a
            np.copyto(_target(a, axis, int(idx.size), out), a)
            return out
        return _gather(a, axis, idx, nt, out)
    return _weighted(a, axis, w, ind, nt, out)


def resize_axis(
    a,
    axis: int,
    weights,
    indices,
    *,
    num_threads: int | None = None,
    out: np.ndarray | None = None,
) -> np.ndarray:
    """MATLAB's resize pass (``imresizemex``) along one axis of a C array: output ``i`` is the
    sum in double of ``weights[i, k] * a[indices[i, k]]`` over ``k`` in order from +0, rounded
    to the class (module docstring); ``indices`` are 0-based. ``out`` takes the result."""
    a = np.ascontiguousarray(a)
    w = np.ascontiguousarray(weights, dtype=np.float64)
    ind = np.ascontiguousarray(indices, dtype=np.int64)
    if w.ndim != 2 or w.shape != ind.shape:
        raise ValueError("weights and indices must be out_length x P")
    return _pass(a, axis, w, ind, resolve_threads(num_threads), out, copy_one_hot=False)


# ---------------------------------------------------------------------------------------------
# arguments
# ---------------------------------------------------------------------------------------------


def _match(value, options: dict, fn: str) -> str:
    """MATLAB's validatestring: case-insensitive, a unique prefix or an exact name"""
    if isinstance(value, str) and value:
        v = value.lower()
        exact = [o for o in options if o == v]
        hits = exact or [o for o in options if o.startswith(v)]
        if len(hits) == 1:
            if options[hits[0]] is None:
                raise NotImplementedError(f"{fn} method {hits[0]!r} is not ported")
            return hits[0]
    if fn == "imresize3":
        raise MatlabError(
            "images:imresize3:oddNumberArgs",
            "Function IMRESIZE3 expected an even number of parameter/value arguments.",
        )
    raise MatlabError(
        "MATLAB:images:imresize:oddNumberArgs",
        "Function IMRESIZE expected an even number of parameter/value arguments.",
    )


def _is_nan(v) -> bool:
    return v is None or (isinstance(v, float | np.floating) and math.isnan(v))


def _numbers(v) -> list:
    if np.isscalar(v) or v is None:
        return [v]
    return list(np.asarray(v, dtype=object).reshape(-1))


def _sizes(
    shape: Sequence[int],
    scale_or_size,
    scale,
    output_size,
    nd: int,
    fn: str,
) -> tuple[list[int], list[float]]:
    """(output sizes, scales) per axis of ``shape`` (Python order) from MATLAB's forms"""
    if scale_or_size is not None:
        vals = _numbers(scale_or_size)
        if len(vals) == 1 and output_size is None:
            scale = vals[0]
        else:
            output_size = scale_or_size
    if output_size is not None:
        vals = _numbers(output_size)
        if len(vals) != nd:
            if fn == "imresize3":
                raise MatlabError(
                    "images:imresize3:badScaleOrSize", "Invalid scale or size input argument."
                )
            if len(vals) == 3:
                raise MatlabError(
                    "MATLAB:images:imresize:invalidOutputSize",
                    "When the second argument is a 1-by-3 vector, IMRESIZE interprets it as an "
                    "M-by-3 colormap.",
                )
            raise MatlabError(
                "MATLAB:images:imresize:badScaleOrSize", "Invalid scale or size input argument."
            )
        nan = [_is_nan(v) for v in vals]
        if any(nan):
            given = [i for i in range(nd) if not nan[i]]
            if fn == "imresize3" and len(given) != 1:
                raise MatlabError(
                    "images:imresize3:invalidOutputSize",
                    "Output size with exactly one NaN is forbidden.  Either all three elements of "
                    "output size can be numeric, or exactly one element can be numeric while the "
                    "other two elements are NaN.",
                )
            if fn == "imresize" and not given:
                raise MatlabError(
                    "MATLAB:images:imresize:allNaN",
                    "Specified output size cannot contain two NaNs.",
                )
            g = given[0]
            _positive([vals[g]], fn)
            og = int(math.ceil(float(vals[g])))
            s = og / shape[g]
            return [og if i == g else int(math.ceil(s * shape[i])) for i in range(nd)], [s] * nd
        _positive(vals, fn)
        out = [int(math.ceil(float(v))) for v in vals]
        return out, [out[i] / shape[i] for i in range(nd)]
    if scale is None:
        raise MatlabError(f"MATLAB:{fn}:notEnoughInputs", "Not enough input arguments.")
    vals = _numbers(scale)
    if len(vals) == 1:
        vals = vals * nd
    if len(vals) != nd:
        raise MatlabError(f"images:{fn}:badScaleOrSize", "Invalid scale or size input argument.")
    _positive(vals, fn)
    sc = [float(v) for v in vals]
    return [int(math.ceil(sc[i] * shape[i])) for i in range(nd)], sc


def _positive(vals, fn: str) -> None:
    for v in vals:
        if _is_nan(v) or not float(v) > 0:
            ident = (
                "MATLAB:imresize3:expectedPositive"
                if fn == "imresize3"
                else ("MATLAB:resizeParseInputs:expectedPositive")
            )
            raise MatlabError(ident, "The output size or scale must be positive.")


def _run(
    a: np.ndarray,
    passes: list[tuple[int, int, float]],
    kernel: str,
    aa: bool,
    nt: int,
    out: np.ndarray | None,
) -> np.ndarray:
    """the passes (axis, output length, scale) in order"""
    b = a
    last = len(passes) - 1
    for j, (axis, n_out, s) in enumerate(passes):
        w, ind = contributions(b.shape[axis], n_out, s, kernel, None, aa)
        b = _pass(b, axis, w, ind, nt, out if j == last else None)
    # every pass an identity: a copy (into out, the last pass did that)
    return a.copy() if b is a else b


# ---------------------------------------------------------------------------------------------
# imresize3, imresize
# ---------------------------------------------------------------------------------------------


def imresize3(
    v,
    scale_or_size=None,
    method: str = "cubic",
    *,
    scale=None,
    output_size=None,
    antialiasing: bool | None = None,
    num_threads: int | None = None,
    out: np.ndarray | None = None,
) -> np.ndarray:
    """MATLAB's ``imresize3(V, scale_or_size, method)`` of a ``(z, y, x)`` volume.

    ``scale_or_size``: a number is the scale of every axis, three values the output size
    ``(nz, ny, nx)`` (rounded up; two of them None or NaN: the third one's scale for all).
    ``scale`` (a number or ``(z, y, x)``) and ``output_size`` are MATLAB's ``'Scale'`` and
    ``'OutputSize'``. ``method``: ``'cubic'`` (default), ``'linear'``, ``'nearest'``, ``'box'``,
    ``'triangle'`` (unique prefixes, any case). ``antialiasing``: when shrinking (default True,
    False for nearest). ``out`` takes the result (C-contiguous, the class of ``v``).
    """
    a = np.asarray(v)
    while a.ndim > 3 and a.shape[0] == 1:
        a = a[0]
    if a.size == 0:
        raise MatlabError(
            "MATLAB:imresize3:expectedNonempty", "Expected input number 1, V, to be nonempty."
        )
    if a.ndim > 3:
        raise MatlabError(
            "images:imresize3:incorrectDimensions",
            "Input image does not have the correct number of dimensions.",
        )
    if a.ndim < 3 or a.shape[0] == 1:
        raise MatlabError(
            "MATLAB:badsubscript",
            "Index exceeds the number of array elements. Index must not exceed 2.",
        )
    a = np.ascontiguousarray(a)
    _class(a.real if np.iscomplexobj(a) else a)
    meth = _match(method, _METHODS3, "imresize3")
    kernel = _METHODS3[meth]
    aa = (meth != "nearest") if antialiasing is None else bool(antialiasing)
    osz, sc = _sizes(a.shape, scale_or_size, scale, output_size, 3, "imresize3")
    # MATLAB's dimensions y, x, z (axes 1, 2, 0) in ascending order of scale
    dims = sorted((1, 2, 0), key=lambda ax: sc[ax])
    return _run(
        a, [(ax, osz[ax], sc[ax]) for ax in dims], kernel, aa, resolve_threads(num_threads), out
    )


def imresize(
    a,
    scale_or_size=None,
    method: str = "bicubic",
    *,
    scale=None,
    output_size=None,
    antialiasing: bool | None = None,
    num_threads: int | None = None,
    out: np.ndarray | None = None,
) -> np.ndarray:
    """MATLAB's ``imresize(A, scale_or_size, method)`` of an image ``(y, x)`` or a stack
    ``(..., y, x)`` (every plane alike).

    ``scale_or_size``: a number is the scale, two values the output size ``(ny, nx)`` (rounded
    up; one None or NaN: the other one's scale for both). ``scale`` (a number or ``(y, x)``) and
    ``output_size`` are MATLAB's ``'Scale'`` and ``'OutputSize'``. ``method``: ``'bicubic'``
    (default), ``'bilinear'``, ``'nearest'``, ``'box'``, ``'triangle'``, ``'cubic'`` (unique
    prefixes, any case). ``antialiasing``: when shrinking (default True, False for nearest).
    """
    arr = np.asarray(a)
    if arr.size == 0:
        raise MatlabError(
            "MATLAB:resizeParseInputs:expectedNonempty",
            "Expected input number 1, A, to be nonempty.",
        )
    if arr.ndim < 2:
        arr = arr.reshape(1, -1)
    arr = np.ascontiguousarray(arr)
    _class(arr.real if np.iscomplexobj(arr) else arr)
    meth = _match(method, _METHODS2, "imresize")
    kernel = _METHODS2[meth]
    aa = (meth != "nearest") if antialiasing is None else bool(antialiasing)
    osz, sc = _sizes(arr.shape[-2:], scale_or_size, scale, output_size, 2, "imresize")
    nd = arr.ndim
    dims = sorted((0, 1), key=lambda i: sc[i])
    passes = [(nd - 2 + i, osz[i], sc[i]) for i in dims]
    return _run(arr, passes, kernel, aa, resolve_threads(num_threads), out)


# ---------------------------------------------------------------------------------------------
# the lab's imresize3_fast and imresize3_average
# ---------------------------------------------------------------------------------------------

_FAST_KERNELS = {"linear": "triangle", "cubic": "cubic"}


@register(
    "imresize3_fast",
    Param("im", None, "Volume (z, y, x).", convert="volume"),
    Param("outSize", None, "Output size (z, y, x), rounded up.", convert="size_yxz"),
    Param("method", None, "'linear' or 'cubic' (other methods: imresize3's result)."),
    Param(
        "outBbox",
        None,
        "Output region (z0, y0, x0, z1, y1, x1) to compute (None: all).",
        convert="bbox",
    ),
    Param("outClass", None, "MATLAB class to resize in (None: the class of im)."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    task="N01",
)
def imresize3_fast(
    im,
    out_size: Sequence[float],
    method: str = "cubic",
    out_bbox: Sequence[int] | None = None,
    out_class: str | None = None,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """PetaKit5D's ``imresize3_fast(im, outSize, method, outBbox, outClass)``: imresize3 for
    'linear' and 'cubic' with its weights and passes, skipping the axes that keep their size.

    ``out_size`` is ``(nz, ny, nx)`` (rounded up). ``out_bbox`` (``(z0, y0, x0, z1, y1, x1)``)
    computes only that region of the output, from the input region it needs (the values of the
    whole output's region). ``out_class``: the MATLAB class to resize in (default: the class of
    ``im``); only the input region used is converted. Other methods crop imresize3's result.
    """
    from petakit5d.core.bbox import refuse_nan_end

    refuse_nan_end(out_bbox)  # (MATLAB indexes with it)
    a = np.asarray(im)
    cls = dtypes.matlab_class(a.dtype) if out_class is None else str(out_class)
    nt = resolve_threads(num_threads)
    if not isinstance(method, str) or method.lower() not in _FAST_KERNELS:
        b = a if dtypes.matlab_class(a.dtype) == cls else np.asarray(dtypes.cast(a, cls))
        r = imresize3(b, out_size, method, num_threads=nt)
        if out_bbox is not None:
            z0, y0, x0, z1, y1, x1 = (int(v) for v in out_bbox)
            r = np.ascontiguousarray(r[z0:z1, y0:y1, x0:x1])
        return r
    kernel = _FAST_KERNELS[method.lower()]
    v = a.reshape((1,) * (3 - a.ndim) + a.shape) if a.ndim < 3 else a
    sz = v.shape
    osz = [int(math.ceil(float(s))) for s in _numbers(out_size)]
    if len(osz) != 3:
        raise ValueError("out_size must have three values (nz, ny, nx)")
    sc = [osz[i] / sz[i] for i in range(3)]
    # MATLAB's sort of [sy sx sz]: y, x, z in ascending order of scale (stable)
    order = sorted((1, 2, 0), key=lambda ax: sc[ax])
    weights: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    region = []
    for ax in range(3):
        lo, hi = (0, osz[ax]) if out_bbox is None else (int(out_bbox[ax]), int(out_bbox[ax + 3]))
        if osz[ax] == sz[ax]:
            region.append((lo, hi))
            continue
        w, ind = contributions(sz[ax], osz[ax], sc[ax], kernel, None, True)
        if out_bbox is not None:
            w, ind = w[lo:hi], ind[lo:hi]
            s = int(ind.min())
            region.append((s, int(ind.max()) + 1))
            ind = ind - s
        else:
            region.append((0, sz[ax]))
        weights[ax] = (np.ascontiguousarray(w), np.ascontiguousarray(ind))
    sub = v[tuple(slice(lo, hi) for lo, hi in region)]
    if dtypes.matlab_class(sub.dtype) != cls:
        sub = np.asarray(dtypes.cast(sub, cls))
    b = np.ascontiguousarray(sub)
    for ax in order:
        if ax in weights:
            w, ind = weights[ax]
            b = _pass(b, ax, w, ind, nt, copy_one_hot=False)
    if np.shares_memory(b, a):
        b = b.copy()
    return b.reshape(b.shape[3 - a.ndim :]) if a.ndim < 3 else b


def _squeeze_matlab(m: np.ndarray) -> np.ndarray:
    """MATLAB's squeeze of an array in MATLAB's layout (at least two dimensions remain)"""
    if m.ndim <= 2:
        return m
    shape = [s for s in m.shape if s != 1]
    # a vector is a column, a single element 1 x 1
    shape = shape + [1] * (2 - len(shape)) if len(shape) < 2 else shape
    return m.reshape(shape, order="F")


@register(
    "imresize3_average",
    Param("im", None, "Volume (z, y, x).", convert="volume"),
    Param("resampleFactor", None, "Block size (z, y, x) to average.", convert="size_yxz"),
    task="N01",
)
def imresize3_average(im, resample_factor: Sequence[int]) -> np.ndarray:
    """PetaKit5D's ``imresize3_average(im, resampleFactor)``: the single mean of each block of
    ``resample_factor`` ``(z, y, x)`` voxels; sizes that the factor does not divide are padded
    with zeros at the end and the edge blocks divided by their voxels inside the image. The
    result has MATLAB's squeezed shape (a factor that leaves one plane along an axis drops it).
    A NaN result has the NaN MATLAB's sum gives it (``reduce.sum``): for blocks of up to 43
    values along y the first NaN of the block in MATLAB's order (y fastest, then x, z)."""
    from petakit5d.core.arrays import from_matlab, to_matlab
    from petakit5d.matlabcompat import reduce

    v = np.asarray(im)
    v = v.reshape((1,) * (3 - v.ndim) + v.shape) if v.ndim < 3 else v
    rz, ry, rx = (int(r) for r in resample_factor)
    rs = (ry, rx, rz)
    nz, ny, nx = v.shape
    sz = (ny, nx, nz)
    s = np.asarray(dtypes.single(v))
    if any(sz[i] % rs[i] for i in range(3)):
        # zero-padded at the end; the edge blocks divided by their voxels inside the image
        sz1 = [(-(-sz[i] // rs[i])) * rs[i] for i in range(3)]
        p = np.zeros((sz1[2], sz1[0], sz1[1]), dtype=np.float32)
        p[:nz, :ny, :nx] = s
        c = np.zeros(p.shape, dtype=np.bool_)
        c[:nz, :ny, :nx] = True
        shp6 = (rs[0], sz1[0] // rs[0], rs[1], sz1[1] // rs[1], rs[2], sz1[2] // rs[2])
        p6 = to_matlab(p).reshape(shp6, order="F")
        c6 = to_matlab(c).reshape(shp6, order="F")
        num = reduce.sum(p6, axis=(0, 2, 4), keepdims=True, matlab_layout=True)
        cnt = reduce.sum(c6, axis=(0, 2, 4), keepdims=True, matlab_layout=True)
        r = dtypes.rdivide(_squeeze_matlab(np.asarray(num)), _squeeze_matlab(np.asarray(cnt)))
    else:
        shp6 = (rs[0], sz[0] // rs[0], rs[1], sz[1] // rs[1], rs[2], sz[2] // rs[2])
        s6 = to_matlab(s).reshape(shp6, order="F")
        r = reduce.mean(s6, axis=(0, 2, 4), keepdims=True, matlab_layout=True)
        r = _squeeze_matlab(np.asarray(r))
    return from_matlab(np.asarray(r))
