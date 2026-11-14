"""MATLAB's Gaussian and linear filters on the port's (z, y, x) C arrays (R2024b, Linux x86-64):
``fspecial('gaussian')``, ``imgaussfilt``, ``imgaussfilt3``, ``imfilter``, ``convn``, ``conv2`` and
``edgetaper``, written from their documentation and from black-box measurements of their results
(task N04; no toolbox source was read). Each output value is MATLAB's arithmetic: the kernels'
formulas, the padding, the class of the computation, and the order and fusing of the
multiply-adds, which MATLAB's code chooses by the pass and by the place of the output:

- ``fspecial('gaussian', hsize, sigma)``: ``exp(-(x.*x + y.*y)/(2*sigma*sigma))`` on the grid
  centred on the kernel, values below ``eps*max`` set to 0, divided by MATLAB's ``sum``.
- ``convn``, ``conv2``: one fused multiply-add chain per output over the kernel in MATLAB's
  column-major order (``conv2``'s order; ``convn`` with a 3-D kernel visits z slowest, x
  descending, y fastest); single data with a single kernel; a kernel larger than the image in
  every dimension swaps the roles and adds plainly. ``conv2(u, v, A)`` filters the columns with
  ``u``, then the rows with ``v``.
- ``imfilter`` with a vector kernel (and so every pass of ``imgaussfilt3``): along MATLAB's rows
  (x) and planes (z) a fused chain in the correlation's order; along the columns (y, MATLAB's
  contiguous dimension) fused chains whose tap order MATLAB's vector code sets by the row: the
  rows in vectors of 8 (double) or 16 (single) outputs visit the taps from the top in chunks of
  4 (8), the two halves of a chunk interleaved. With zero padding these vectors run while the
  kernel's lower reach stays inside the column; the rows after them take the top chunk
  ascending, and the rows whose window passes the column's end the top chunk descending. With
  the other paddings the padded column goes in blocks of 8 KB of padded rows, each block after
  the first starting m - 1 rows back. Each block's whole vectors take the body order (in single
  with a three-tap remainder, the last vector of a block of whole vectors adds the remainder as
  a chain of its own). Its other rows take each chunk as an ascending chain, the chains added
  from the top, then the remainder's chain; a one-tap remainder is fused last into the sum.
  Along x, rows and planes read their data first and the padding after it (along x with
  replicate or constant padding: the padded taps at an end as one product with their weights'
  sum). Arrays with less than a tenth nonzero
  values (kernels up to 15 taps) take a plain path: each product rounded, added in the
  correlation's order.
- ``imgaussfilt``: ``padarray`` by the kernel's half size, then for kernels of 49 values or more
  ``conv2(gy, gx, A, 'valid')`` in double (single and integer data too), smaller kernels a 2-D
  kernel through ``imfilter`` in the data's class (single for integers); the frequency domain
  for large images and kernels (``FilterDomain`` 'auto', MATLAB's choice with the Intel IPP
  library, the default).
- ``imgaussfilt3``: three ``imfilter`` passes with the 1-D kernels, x, then y, then z, in single
  (double for double data), cast to the class at the end.
- ``edgetaper``: the weights from the circular autocorrelation of the PSF's projections over
  ``n - 1`` points, closed periodically; the blur is the circular convolution with the PSF.

What is not bit for bit MATLAB's (measured by the oracle cases; tolerances proposed for
``docs/parity/matlabcompat.md``): the column pass with zero padding, the rows of its last
vectors (double 7 taps or more, single 11 or more; 'full' from 5 taps), imgaussfilt's small 2-D
kernels and imfilter's general 2-D and 3-D kernels (MATLAB's tap order there not measured),
imgaussfilt's separable path for kernels of 21 or more taps (conv2's long kernels), the
frequency domain path and edgetaper (rounding level: half spectra of ``petakit5d.fft``).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from petakit5d.core import dtypes
from petakit5d.core.threads import resolve_threads

__all__ = [
    "MatlabError",
    "conv2",
    "convn",
    "edgetaper",
    "fspecial",
    "gaussian_kernel_1d",
    "imfilter",
    "imgaussfilt",
    "imgaussfilt3",
]

_EPS = float(np.finfo(np.float64).eps)
_ext_state: dict = {"mod": None}

# imfilter's switch to the plain path: fewer than a tenth nonzero values, kernels up to this
# many values (measured: 15 taps switch, 17 do not)
_SPARSE_FRACTION = 10
_SPARSE_MAX_TAPS = 15
# imgaussfilt: kernels of this many values or more are filtered separably (conv2)
_SEPARABLE_FROM = 49


class MatlabError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _ext():
    if _ext_state["mod"] is None:
        from petakit5d._ext import load

        _ext_state["mod"] = load("filters")
    return _ext_state["mod"]


# ---------------------------------------------------------------------------------------------
# kernels
# ---------------------------------------------------------------------------------------------


def _exp(x: np.ndarray) -> np.ndarray:
    from petakit5d.matlabcompat import elementary

    return elementary.exp(np.asarray(x, dtype=np.float64))


def _msum(v: np.ndarray) -> float:
    from petakit5d.matlabcompat import reduce

    s = reduce.sum(np.asarray(v, dtype=np.float64), axis=None, matlab_layout=True)
    return float(np.asarray(s).reshape(-1)[0])


def gaussian_kernel_1d(sigma: float, size: int) -> np.ndarray:
    """The normalized 1-D Gaussian of imgaussfilt and imgaussfilt3: ``exp(-x.^2/(2 sigma^2))`` on
    ``-(size-1)/2 : (size-1)/2``, divided by MATLAB's sum (double)."""
    r = (int(size) - 1) / 2
    x = np.arange(-r, r + 0.5, 1.0)
    g = _exp(-(x * x) / (2 * sigma * sigma))
    return g / _msum(g)


def _gaussian_2d(sy: float, sx: float, ky: int, kx: int) -> np.ndarray:
    """The 2-D Gaussian of imgaussfilt's small kernels, (ky, kx): exp of the sum of both terms,
    divided by MATLAB's sum over the kernel."""
    ry, rx = (ky - 1) / 2, (kx - 1) / 2
    X, Y = np.meshgrid(np.arange(-rx, rx + 0.5, 1.0), np.arange(-ry, ry + 0.5, 1.0))
    h = _exp(-(X * X / (2 * sx * sx) + Y * Y / (2 * sy * sy)))
    return h / _msum(h)


def fspecial(kind: str, hsize=(3, 3), sigma: float = 0.5) -> np.ndarray:
    """``fspecial('gaussian', hsize, sigma)`` (the only type the port's callers use): an
    ``(rows, cols)`` double kernel; a scalar ``hsize`` is square."""
    if not isinstance(kind, str) or kind.lower() != "gaussian":
        raise MatlabError("images:fspecial:unsupported", f"fspecial type {kind!r} is not ported")
    hs = np.atleast_1d(np.asarray(hsize, dtype=np.float64)).ravel()
    if hs.size == 1:
        hs = np.array([hs[0], hs[0]])
    if hs.size != 2 or np.any(hs < 1) or np.any(hs != np.floor(hs)):
        raise MatlabError(
            "images:fspecial:wrongSizeHSize",
            "HSIZE must be a positive integer or a 2-element vector",
        )
    sigma = float(sigma)
    siz = (hs - 1) / 2
    x, y = np.meshgrid(np.arange(-siz[1], siz[1] + 0.5, 1.0), np.arange(-siz[0], siz[0] + 0.5, 1.0))
    h = _exp(-(x * x + y * y) / (2 * sigma * sigma))
    h[h < _EPS * h.max()] = 0
    s = _msum(h)
    if s != 0:
        h = h / s
    return h


# ---------------------------------------------------------------------------------------------
# arguments
# ---------------------------------------------------------------------------------------------

_PAD_CODES = {"replicate": 1, "symmetric": 2, "circular": 3}
_FLOAT = {"single": np.float32, "double": np.float64}


def _match(value: str, options: Sequence[str]) -> str | None:
    """MATLAB's validatestring: a case-insensitive unique prefix (an exact name wins)."""
    v = value.lower()
    hits = [o for o in options if o.lower().startswith(v)]
    exact = [o for o in options if o.lower() == v]
    if exact:
        return exact[0]
    return hits[0] if len(hits) == 1 else None


def _padding(p) -> tuple[int, float]:
    """(mode code, constant): 0 a constant value, 1 replicate, 2 symmetric, 3 circular."""
    if isinstance(p, str):
        m = _match(p, ("replicate", "symmetric", "circular"))
        if m is None:
            raise MatlabError("images:imgaussfilt:invalidPadding", f"unknown padding {p!r}")
        return _PAD_CODES[m], 0.0
    v = np.asarray(p)
    if v.size != 1 or v.dtype.kind not in "biuf":
        raise MatlabError("images:imgaussfilt:invalidPadding", "Padding must be a name or a number")
    return 0, float(v.ravel()[0])


def _pad_value(value: float, cls: str) -> float:
    """a constant padding value in the class of the array, as padarray converts it"""
    if cls in _FLOAT or value == 0:
        return value
    if cls == "logical":
        return float(value != 0)
    return float(np.asarray(dtypes.cast(value, cls)))


def _compute_class(cls: str) -> type:
    """the class imgaussfilt3 and imgaussfilt's small kernels compute in: single, double for
    double and the 32-bit integers"""
    return np.float64 if cls in ("double", "int32", "uint32") else np.float32


def _class_of(a: np.ndarray) -> str:
    return dtypes.matlab_class(a.dtype)


def _out_cast(y: np.ndarray, cls: str, nt: int | None = None) -> np.ndarray:
    """the result in the class: MATLAB's cast (round half away from zero, saturate, NaN to 0)"""
    if cls in _FLOAT:
        return y.astype(_FLOAT[cls], copy=False)
    if cls == "logical":
        return np.asarray(dtypes.cast(y, cls))
    from petakit5d.kernels import ops

    return ops.convert(np.ascontiguousarray(y), cls, mode="cast", num_threads=nt)


# ---------------------------------------------------------------------------------------------
# tap programs (python/cpp/filters): one row of ints per program; >= 0 a fused tap, <= -10 a
# plain tap (-op - 10), -1 the end of a chain (the chains are added in order), -3 the end of a
# chain after which the taps continue from the sum so far, -2 padding
# ---------------------------------------------------------------------------------------------


def _programs(rows: list[list[int]]) -> tuple[np.ndarray, np.ndarray]:
    """unique programs and the row map"""
    uniq: dict[tuple, int] = {}
    rowmap = np.empty(len(rows), dtype=np.int32)
    for i, r in enumerate(rows):
        rowmap[i] = uniq.setdefault(tuple(r), len(uniq))
    plen = max((len(k) for k in uniq), default=1)
    prog = np.full((max(len(uniq), 1), plen), -2, dtype=np.int32)
    for k, j in uniq.items():
        prog[j, : len(k)] = k
    return prog, rowmap


def _body_order(m: int, half: int) -> list[int]:
    """the column pass's vector body: taps from the top in chunks of 2*half, each chunk's upper
    and lower halves interleaved (descending in each)"""
    order: list[int] = []
    a = m - 1
    while a >= 0:
        top = [t for t in range(a, a - half, -1) if t >= 0]
        bot = [t for t in range(a - half, a - 2 * half, -1) if t >= 0]
        for k in range(max(len(top), len(bot))):
            if k < len(top):
                order.append(top[k])
            if k < len(bot):
                order.append(bot[k])
        a -= 2 * half
    return order


def _column_rows(n_out: int, m: int, pre: int, double: bool, zero_pad: bool) -> list[list[int]]:
    """the tap order of every output row of imfilter's column pass (MATLAB's y). With zero
    padding the vector body ends before the vectors whose windows pass the column's end. A padded
    column is taken in blocks of 8 KB of padded rows (1,024 double, 2,048 single), each block
    after the first starting m - 1 rows back, and the body of each block covers its whole
    vectors."""
    vec, half = (8, 2) if double else (16, 4)
    body = _body_order(m, half)
    c = min(m, 2 * half)
    if not zero_pad:
        tail = _padded_tail(m, 2 * half)
        rows: list[list[int]] = []
        step = (1024 if double else 2048) - (m - 1)
        if step <= 0:
            step = n_out
        while len(rows) < n_out:
            cnt = min(step, n_out - len(rows))
            nbody = vec * (cnt // vec)
            block = [body] * nbody + [tail] * (cnt - nbody)
            if not double and nbody == cnt and cnt and m % (2 * half) == 3:
                # single, a block of whole vectors and a remainder chunk of three taps: the
                # last vector adds the remainder as an ascending chain of its own to the whole
                # chunks (measured for 3, 11 and 19 taps, columns of 32 to 5,000 rows)
                nfull = m - 3
                last = body[:nfull] + ([-1] if nfull else []) + [0, 1, 2]
                block[cnt - vec :] = [last] * vec
            rows.extend(block)
        return rows
    post = m - 1 - pre
    nbody = vec * max(0, (n_out - post) // vec)
    tail = sorted(body[:c]) + body[c:]
    border = sorted(body[:c], reverse=True) + body[c:]
    rows = []
    for i in range(n_out):
        if i < nbody:
            rows.append(body)
        elif n_out - 1 - i >= post:
            rows.append(tail)
        else:
            rows.append(border)
    return rows


def _padded_tail(m: int, chunk: int) -> list[int]:
    """the column pass's rows after the vector body of an explicitly padded column: the chunks
    of taps from the top, each an ascending chain, added from the top chunk down, then the
    remainder's chain; a remainder of one tap is fused last into the sum of the chunks (measured
    for 3 to 31 taps)"""
    if m <= chunk:
        return list(range(m))
    full = []
    a = m
    while a - chunk >= 0:
        full.append(list(range(a - chunk, a)))
        a -= chunk
    rem = list(range(a))
    out: list[int] = []
    for k, ch in enumerate(full):
        if k:
            out.append(-1)
        out.extend(ch)
    if len(rem) == 1:
        if len(full) > 1:
            out.append(-3)
        out.extend(rem)
    elif rem:
        out.append(-1)
        out.extend(rem)
    return out


def _line_pass(
    x: np.ndarray,
    w: np.ndarray,
    axis: int,
    rows: list[list[int]],
    pre: int,
    n_out: int,
    mode: int,
    value: float,
    nt: int,
    span: int | None = None,
) -> np.ndarray:
    """a 1-D pass along ``axis`` of a C array with per-row tap programs (``w`` may hold the edge
    taps' weights after the ``span`` taps of the kernel)"""
    shape = x.shape
    ax = axis % x.ndim
    outer = int(np.prod(shape[:ax], dtype=np.int64))
    inner = int(np.prod(shape[ax + 1 :], dtype=np.int64))
    out_shape = shape[:ax] + (n_out,) + shape[ax + 1 :]
    out = np.empty(out_shape, dtype=x.dtype)
    prog, rowmap = _programs(rows)
    fn = _ext().line_program_f64 if x.dtype == np.float64 else _ext().line_program_f32
    fn(
        np.ascontiguousarray(x),
        out,
        outer,
        shape[ax],
        inner,
        n_out,
        np.ascontiguousarray(w, dtype=x.dtype),
        prog,
        rowmap,
        pre,
        len(w) if span is None else span,
        mode,
        value,
        nt,
    )
    return out


_EDGE_LOW = 1 << 21
_EDGE_HIGH = 1 << 22


def _flat_rows(
    n_in: int, n_out: int, m: int, pre: int, mode: int, value: float, w: np.ndarray, rows_axis: bool
) -> tuple[list[list[int]], np.ndarray]:
    """the tap order of imfilter's row (x) and plane (z) passes: a fused chain in the
    correlation's order; at the borders of an explicitly padded array the taps that read data
    first, then the padded ones; along the rows (x) with replicate or constant padding the padded
    taps at each end as one product, last, of the edge value or the constant with their weights'
    sum (measured)."""
    asc = list(range(m))
    if mode == 0 and value == 0:
        return [asc] * n_out, w
    rows: list[list[int]] = []
    extra: list = []
    for i in range(n_out):
        low = [t for t in range(m) if i + t - pre < 0]
        high = [t for t in range(m) if i + t - pre >= n_in]
        if not low and not high:
            rows.append(asc)
            continue
        valid = [t for t in range(m) if 0 <= i + t - pre < n_in]
        merge_low = merge_high = rows_axis and mode in (0, 1)
        if not (merge_low and low) and not (merge_high and high):
            rows.append(valid + low + high)
            continue
        r = list(valid) + ([] if merge_low else low)
        for side, taps, merge in ((_EDGE_LOW, low, merge_low), (_EDGE_HIGH, high, merge_high)):
            if taps and merge:
                ws = w.dtype.type(0)
                for t in taps:
                    ws = w.dtype.type(ws + w[t])
                extra.append(ws)
                r.append(-(side + m + len(extra) - 1))
        rows.append(r)
    if extra:
        w = np.concatenate([w, np.array(extra, dtype=w.dtype)])
    return rows, w


def _sparse(x: np.ndarray, numel_h: int) -> bool:
    """imfilter's plain path: fewer than a tenth nonzero values (counted until the tenth is
    reached)"""
    if numel_h > _SPARSE_MAX_TAPS:
        return False
    flat = x.reshape(-1)
    need = -(-flat.size // _SPARSE_FRACTION)  # nonzeros that make the array dense
    step = max(1 << 20, need)
    found = 0
    for k in range(0, flat.size, step):
        found += int(np.count_nonzero(flat[k : k + step]))
        if found >= need:
            return False
    return True


# ---------------------------------------------------------------------------------------------
# imfilter
# ---------------------------------------------------------------------------------------------


def _imfilter_options(options) -> tuple[int, float, str, str]:
    mode, value, shape, corr = 0, 0.0, "same", "corr"
    for o in options:
        if isinstance(o, str):
            m = _match(o, ("same", "full", "corr", "conv", "replicate", "symmetric", "circular"))
            if m is None:
                raise MatlabError("images:imfilter:invalidInputStr", f"unknown option {o!r}")
            if m in ("same", "full"):
                shape = m
            elif m in ("corr", "conv"):
                corr = m
            else:
                mode, value = _PAD_CODES[m], 0.0
        else:
            mode, value = 0, float(np.asarray(o).ravel()[0])
    return mode, value, shape, corr


def _vector_axis(h: np.ndarray) -> int | None:
    """the axis of a kernel with one non-singleton axis (None: a general kernel)"""
    big = [k for k, n in enumerate(h.shape) if n > 1]
    if len(big) == 1:
        return big[0] - h.ndim
    if not big:
        return -1
    return None


def _imfilter_float(
    a: np.ndarray, h: np.ndarray, mode: int, value: float, shape: str, nt: int
) -> np.ndarray:
    """correlation of a float array (its class) with a double kernel (Python layout), MATLAB's
    padding and output size"""
    dt = a.dtype
    nd = max(a.ndim, h.ndim)
    a = a.reshape((1,) * (nd - a.ndim) + a.shape)
    h = np.asarray(h, dtype=np.float64).reshape((1,) * (nd - h.ndim) + h.shape)
    sparse = _sparse(a, h.size)
    ax = _vector_axis(h)
    if ax is not None:
        m = h.shape[ax]
        w = h.reshape(-1).astype(dt)
        n = a.shape[ax]
        pre = (m - 1) // 2 if shape == "same" else m - 1
        n_out = n if shape == "same" else n + m - 1
        if sparse:
            rows = [[-(t + 10) for t in range(m)]] * n_out
        elif ax == -2:
            rows = _column_rows(n_out, m, pre, dt == np.float64, mode == 0 and value == 0)
        else:
            rows, w = _flat_rows(n, n_out, m, pre, mode, value, w, ax == -1)
        return _line_pass(a, w, ax, rows, pre, n_out, mode, value, nt, span=m)
    # a general kernel: one chain per output (MATLAB's own order not measured for every size)
    return _taps_correlation(a, h, mode, value, shape, not sparse, nt)


def _pad_array(a: np.ndarray, pre, post, mode: int, value: float) -> np.ndarray:
    widths = list(zip(pre, post, strict=True))
    if mode == 0:
        return np.pad(a, widths, mode="constant", constant_values=a.dtype.type(value))
    return np.pad(a, widths, mode={1: "edge", 2: "symmetric", 3: "wrap"}[mode])


def _taps_correlation(a, h, mode, value, shape, fused, nt) -> np.ndarray:
    """correlation with a general kernel through the tap kernel: the chain over the taps in
    MATLAB's column-major order (fused: x descending, as conv2's order for the flipped kernel)"""
    dt = a.dtype
    a3 = a.reshape((1,) * (3 - a.ndim) + a.shape) if a.ndim < 3 else a
    h3 = h.reshape((1,) * (3 - h.ndim) + h.shape) if h.ndim < 3 else h
    if a3.ndim != 3:
        raise MatlabError(
            "images:imfilter:dims", "imfilter of arrays of more than 3 dimensions is not ported"
        )
    m = h3.shape
    if shape == "same":
        pre = [(k - 1) // 2 for k in m]
        post = [k - 1 - p for k, p in zip(m, pre, strict=True)]
    else:
        pre = post = [k - 1 for k in m]
    P = _pad_array(a3, pre, post, mode, value)
    out_shape = tuple(s - k + 1 for s, k in zip(P.shape, m, strict=True))
    taps = []
    for z in range(m[0]):
        for x in range(m[2]):
            for y in range(m[1]):
                taps.append((z, y, x))
    if fused:
        # conv2 order of the flipped kernel: z, x descending, y ascending over its taps
        taps.sort(key=lambda t: (-t[0], -t[2], t[1]))
    else:
        taps.sort(key=lambda t: (t[0], t[2], t[1]))
    offs = np.array(taps, dtype=np.int64).reshape(-1, 3)
    w = np.array([h3[t] for t in taps], dtype=dt)
    out = np.empty(out_shape, dtype=dt)
    fn = _ext().taps_nd_f64 if dt == np.float64 else _ext().taps_nd_f32
    fn(np.ascontiguousarray(P), out, offs, w, bool(fused), nt)
    return out.reshape(out_shape[3 - a.ndim :] if a.ndim < 3 else out_shape)


def imfilter(a, h, *options, num_threads: int | None = None) -> np.ndarray:
    """``imfilter(A, h, options...)``: options as MATLAB's (a padding name or value, 'same' or
    'full', 'corr' or 'conv'); ``h`` a double kernel in the Python layout (a column kernel is
    ``(m, 1)``, a row ``(1, m)``, a z kernel ``(m, 1, 1)``). The result has A's class (integers
    computed in double, rounded and saturated)."""
    a = np.asarray(a)
    h = np.asarray(h)
    if h.dtype != np.float64:
        raise MatlabError("MATLAB:conv2:NonDoubleFilter", "The filter must be double")
    mode, value, shape, corr = _imfilter_options(options)
    if corr == "conv":
        h = h[tuple(slice(None, None, -1) for _ in range(h.ndim))]
        # MATLAB's center of the flipped kernel for even sizes
    cls = _class_of(a)
    nt = resolve_threads(num_threads)
    if a.size == 0:
        return a.copy()
    value = _pad_value(value, cls)
    if cls == "logical":
        # computed in double; true where the value rounds to a nonzero integer
        y = _imfilter_float(a.astype(np.float64), h, mode, value, shape, nt)
        return np.asarray(dtypes.mround(y)) != 0
    dt = np.float32 if cls == "single" else np.float64
    y = _imfilter_float(a.astype(dt, copy=False), h, mode, value, shape, nt)
    return _out_cast(y, cls, nt)


# ---------------------------------------------------------------------------------------------
# convn and conv2
# ---------------------------------------------------------------------------------------------


def _conv_common(a: np.ndarray, b: np.ndarray, shape: str, nt: int, order: str) -> np.ndarray:
    """full / same / valid convolution of (z, y, x) arrays in MATLAB's tap order"""
    shape = _match(shape, ("full", "same", "valid")) if isinstance(shape, str) else None
    if shape is None:
        raise MatlabError("MATLAB:conv2:unknownShapeParameter", "Unknown shape parameter")
    dt = np.float32 if (a.dtype == np.float32 or b.dtype == np.float32) else np.float64
    nd = max(a.ndim, b.ndim, 3)
    A = np.asarray(a, dtype=dt).reshape((1,) * (nd - a.ndim) + a.shape)
    Bk = np.asarray(b, dtype=dt).reshape((1,) * (nd - b.ndim) + b.shape)
    if nd > 3:
        raise MatlabError(
            "MATLAB:convn:dims", "convn of arrays of more than 3 dimensions is not ported"
        )
    na, nb = A.shape, Bk.shape
    full = (
        tuple(x + y - 1 for x, y in zip(na, nb, strict=True))
        if A.size and Bk.size
        else tuple(max(x, y) for x, y in zip(na, nb, strict=True))
    )
    swap = all(y > x for x, y in zip(na, nb, strict=True) if x > 1 or y > 1) and A.size < Bk.size
    if swap:
        A, Bk = Bk, A
        na, nb = A.shape, Bk.shape
    # the full convolution as a correlation of A padded by nb - 1 with the flipped kernel
    P = np.pad(A, [(k - 1, k - 1) for k in nb])
    taps = [(z, y, x) for z in range(nb[0]) for y in range(nb[1]) for x in range(nb[2])]
    # tap (kz, ky, kx) of the kernel reads P at (o + nb - 1 - k): MATLAB's order over k
    if swap:
        taps.sort(key=lambda t: (t[0], t[2], t[1]))
    elif order == "3d":
        taps.sort(key=lambda t: (t[0], -t[2], t[1]))
    else:
        taps.sort(key=lambda t: (t[0], t[2], t[1]))
    offs = np.array(
        [(nb[0] - 1 - z, nb[1] - 1 - y, nb[2] - 1 - x) for z, y, x in taps], dtype=np.int64
    )
    w = np.array([Bk[t] for t in taps], dtype=dt) if taps else np.zeros(0, dtype=dt)
    out = np.zeros(full, dtype=dt)
    if A.size and Bk.size:
        fn = _ext().taps_nd_f64 if dt == np.float64 else _ext().taps_nd_f32
        fn(np.ascontiguousarray(P), out, offs.reshape(-1, 3), w, not swap, nt)
    na0 = np.asarray(a).shape
    na0 = (1,) * (nd - len(na0)) + na0
    nb0 = (1,) * (nd - np.ndim(b)) + np.shape(b)
    if shape == "same":
        st = [int(math.ceil((k - 1) / 2)) for k in nb0]
        sl = tuple(slice(s, s + n) for s, n in zip(st, na0, strict=True))
        out = out[sl]
    elif shape == "valid":
        sz = [max(n - k + 1, 0) if k else n for n, k in zip(na0, nb0, strict=True)]
        sl = tuple(slice(k - 1, k - 1 + s) for k, s in zip(nb0, sz, strict=True))
        out = out[sl]
        if any(s == 0 for s in sz):
            out = np.zeros(sz, dtype=dt)
    keep = max(np.ndim(a), np.ndim(b))
    return np.ascontiguousarray(out.reshape(out.shape[nd - keep :]))


def convn(a, b, shape: str = "full", *, num_threads: int | None = None) -> np.ndarray:
    """``convn(A, B, shape)`` ('full', 'same', 'valid') of (z, y, x) arrays; single if either is
    single."""
    a = np.asarray(a)
    b = np.asarray(b)
    order = (
        "3d"
        if (np.ndim(b) == 3 and np.shape(b)[0] > 1 and (np.shape(b)[1] > 1 or np.shape(b)[2] > 1))
        else "2d"
    )
    return _conv_common(a, b, shape, resolve_threads(num_threads), order)


def conv2(*args, num_threads: int | None = None) -> np.ndarray:
    """``conv2(A, B, shape)`` or ``conv2(u, v, A, shape)`` (columns with ``u``, then rows with
    ``v``) on 2-D arrays."""
    nt = resolve_threads(num_threads)
    shape = "full"
    if args and isinstance(args[-1], str):
        shape = args[-1]
        args = args[:-1]
    if len(args) == 2:
        a, b = (np.asarray(v) for v in args)
        if a.ndim > 2 or b.ndim > 2:
            raise MatlabError("MATLAB:conv2:inputsMustBe2D", "N-D arrays are not supported")
        return _conv_common(a, b, shape, nt, "2d")
    if len(args) == 3:
        u, v, a = (np.asarray(x) for x in args)
        u = u.reshape(-1, 1)
        v = v.reshape(1, -1)
        t = _conv_common(a, u, shape, nt, "2d")
        return _conv_common(t, v, shape, nt, "2d")
    raise MatlabError("MATLAB:conv2:inputs", "conv2 takes (A, B[, shape]) or (u, v, A[, shape])")


# ---------------------------------------------------------------------------------------------
# imgaussfilt and imgaussfilt3
# ---------------------------------------------------------------------------------------------

_GAUSS_CLASSES = ("uint8", "uint16", "uint32", "int8", "int16", "int32", "single", "double")


def _check_gauss_input(a: np.ndarray, fname: str) -> str:
    if np.iscomplexobj(a):
        raise MatlabError(f"MATLAB:{fname}:expectedReal", "Expected A to be real.")
    cls = _class_of(a)
    if cls not in _GAUSS_CLASSES:
        raise MatlabError(
            f"MATLAB:{fname}:invalidType",
            "Expected A to be one of these types:\n\nuint8, uint16, uint32, int8, int16, int32, "
            f"single, double\n\nInstead its type was {cls}.",
        )
    return cls


def _sigmas(sigma, nd: int, fname: str) -> np.ndarray:
    s = np.atleast_1d(np.asarray(sigma, dtype=np.float64)).ravel()
    if s.size == 1:
        s = np.repeat(s, nd)
    if s.size != nd:
        msg = (
            "Sigma can be a scalar or 3-element vector."
            if nd == 3
            else "Sigma can be a scalar or 2-element vector."
        )
        raise MatlabError(f"images:imgaussfilt:invalidLength{'3' if nd == 3 else ''}", msg)
    if np.any(~(s > 0)) or np.any(~np.isfinite(s)):
        raise MatlabError(f"MATLAB:{fname}:expectedPositive", "Expected Sigma to be positive.")
    return s


def _filter_sizes(filter_size, s: np.ndarray, nd: int) -> np.ndarray:
    if filter_size is None:
        return np.array([2 * int(math.ceil(2 * q)) + 1 for q in s])
    f = np.atleast_1d(np.asarray(filter_size, dtype=np.float64)).ravel()
    if f.size == 1:
        f = np.repeat(f, nd)
    tag = "Three" if nd == 3 else "Two"
    if f.size != nd:
        raise MatlabError(
            f"MATLAB:validate{tag}DFilterSize:incorrectSize", "filterSize has the wrong size."
        )
    if np.any(f < 1) or np.any(f != np.floor(f)):
        raise MatlabError(
            f"MATLAB:validate{tag}DFilterSize:expectedPositive",
            "Expected filterSize to be positive.",
        )
    if np.any(f % 2 == 0):
        raise MatlabError(
            f"MATLAB:validate{tag}DFilterSize:expectedOdd", "Expected filterSize to be odd."
        )
    return f.astype(int)


def _domain(d: str) -> str:
    m = _match(d, ("auto", "spatial", "frequency")) if isinstance(d, str) else None
    if m is None:
        raise MatlabError("images:imgaussfilt:invalidFilterDomain", f"unknown FilterDomain {d!r}")
    return m


def _auto_frequency_2d(shape: tuple[int, int], k: tuple[int, int], cls: str) -> bool:
    """imgaussfilt's 'auto' choice with the IPP library (measured on images of 1.2 to 9 million
    pixels: the frequency domain above 2^21 pixels, single images from 121-point kernels,
    double and integer images from 201)."""
    n = shape[0] * shape[1]
    kk = max(k)
    if n <= 2**21:
        return False
    return kk >= (121 if cls == "single" else 201)


def _spatial_2d(
    a: np.ndarray, s: np.ndarray, f: np.ndarray, mode: int, value: float, cls: str, nt: int
):
    ky, kx = int(f[0]), int(f[1])
    if ky * kx >= _SEPARABLE_FROM:
        gy = gaussian_kernel_1d(s[0], ky)
        gx = gaussian_kernel_1d(s[1], kx)
        ry, rx = (ky - 1) // 2, (kx - 1) // 2
        ny, nx = a.shape[-2], a.shape[-1]
        dy, dx = [list(range(ky - 1, -1, -1))] * ny, [list(range(kx - 1, -1, -1))] * nx
        if mode == 0 and value != 0:
            # padarray first: the padded columns are filtered along y too
            P = _pad_array(a.astype(np.float64), (ry, rx), (ry, rx), mode, value)
            t = _line_pass(P, gy, -2, dy, 0, ny, 0, 0.0, nt)
            y = _line_pass(t, gx, -1, dx, 0, nx, 0, 0.0, nt)
        else:
            # the same values with each pass padding its own axis
            a64 = a.astype(np.float64, copy=False)
            t = _line_pass(a64, gy, -2, dy, ry, ny, mode, value, nt)
            y = _line_pass(t, gx, -1, dx, rx, nx, mode, value, nt)
        return _out_cast(y, cls, nt)
    h = _gaussian_2d(s[0], s[1], ky, kx)
    dt = _compute_class(cls)
    y = _imfilter_float(a.astype(dt, copy=False), h, mode, value, "same", nt)
    return _out_cast(y, cls, nt)


def _frequency(
    a: np.ndarray, kernels: list[np.ndarray], mode: int, value: float, cls: str, nt: int
):
    """the frequency domain: the array padded by the kernels' half sizes, multiplied by the
    transform of the kernel (the kernel's centre at the origin), transformed back, cropped; half
    spectra of the CPU FFT backend (petakit5d.fft; MATLAB's bits not reproduced: within
    rounding)"""
    from petakit5d import fft

    dt = np.float64 if cls == "double" else np.float32
    rs = [(len(g) - 1) // 2 for g in kernels]
    P = _pad_array(a.astype(dt), rs, rs, mode, value)
    h = kernels[0]
    for g in kernels[1:]:
        h = np.multiply.outer(h, g)
    hs = fft.HalfSpectrum(P.shape, dt, threads=nt)
    Y = hs.convolve(P, hs.otf(h.astype(dt)))
    y = Y[tuple(slice(r, r + n) for r, n in zip(rs, a.shape, strict=True))]
    return _out_cast(np.ascontiguousarray(y), cls, nt)


def imgaussfilt(
    a,
    sigma=0.5,
    *,
    filter_size=None,
    padding="replicate",
    filter_domain: str = "auto",
    num_threads: int | None = None,
) -> np.ndarray:
    """``imgaussfilt(A, sigma, 'FilterSize', f, 'Padding', p, 'FilterDomain', d)``: a 2-D Gaussian
    filter of an image ``(y, x)`` (a stack ``(..., y, x)`` plane by plane); ``sigma`` and
    ``filter_size`` a scalar or ``(y, x)``."""
    a = np.asarray(a)
    cls = _check_gauss_input(a, "imgaussfilt")
    s = _sigmas(sigma, 2, "imgaussfilt")
    f = _filter_sizes(filter_size, s, 2)
    mode, value = _padding(padding)
    dom = _domain(filter_domain)
    nt = resolve_threads(num_threads)
    if a.size == 0:
        return a.copy()
    if a.ndim < 2:
        a2 = a.reshape(1, -1)
    else:
        a2 = a
    lead = a2.shape[:-2]
    planes = a2.reshape((-1,) + a2.shape[-2:])
    use_freq = dom == "frequency" or (
        dom == "auto" and _auto_frequency_2d(planes.shape[1:], tuple(f), cls)
    )
    out = []
    for k in range(planes.shape[0]):
        p = planes[k]
        if use_freq:
            g = [gaussian_kernel_1d(s[0], f[0]), gaussian_kernel_1d(s[1], f[1])]
            out.append(_frequency(p, g, mode, value, cls, nt))
        else:
            out.append(_spatial_2d(p, s, f, mode, value, cls, nt))
    y = np.stack(out).reshape(lead + a2.shape[-2:])
    return y.reshape(a.shape)


def imgaussfilt3(
    a,
    sigma=0.5,
    *,
    filter_size=None,
    padding="replicate",
    filter_domain: str = "auto",
    num_threads: int | None = None,
) -> np.ndarray:
    """``imgaussfilt3(A, sigma, ...)``: a 3-D Gaussian filter of a volume ``(z, y, x)`` (an image
    is one plane); ``sigma`` and ``filter_size`` a scalar or ``(z, y, x)``."""
    a = np.asarray(a)
    cls = _check_gauss_input(a, "imgaussfilt3")
    s = _sigmas(sigma, 3, "imgaussfilt3")
    f = _filter_sizes(filter_size, s, 3)
    mode, value = _padding(padding)
    dom = _domain(filter_domain)
    nt = resolve_threads(num_threads)
    if a.size == 0:
        return a.copy()
    v = a.reshape((1,) * (3 - a.ndim) + a.shape) if a.ndim < 3 else a
    if v.ndim > 3:
        raise MatlabError(
            "images:imgaussfilt3:dims",
            "imgaussfilt3 of arrays of more than 3 dimensions is not ported",
        )
    # Python sigma (z, y, x)
    g = [gaussian_kernel_1d(s[d], f[d]) for d in range(3)]
    if dom == "frequency":
        y = _frequency(v, g, mode, value, cls, nt)
        return y.reshape(a.shape)
    dt = _compute_class(cls)
    x = v.astype(dt, copy=False)
    # MATLAB filters x (its dimension 2), then y, then z, each an imfilter pass
    for ax, w in ((-1, g[2]), (-2, g[1]), (-3, g[0])):
        if len(w) == 1 and w[0] == 1.0:
            continue
        h = w.reshape([-1 if d == ax % 3 else 1 for d in range(3)])
        x = _imfilter_float(x, h, mode, value, "same", nt)
    return _out_cast(x, cls, nt).reshape(a.shape)


# ---------------------------------------------------------------------------------------------
# edgetaper
# ---------------------------------------------------------------------------------------------


def edgetaper(i, psf, *, num_threads: int | None = None) -> np.ndarray:
    """``edgetaper(I, PSF)``: ``alpha .* I + (1 - alpha) .* blurred`` with the circular blur by the
    PSF and ``alpha`` the product of ``1 - beta`` along each axis, ``beta`` the normalized
    circular autocorrelation of the PSF's projection over ``n - 1`` points, closed with its
    first value (rounding level against MATLAB)."""
    img = np.asarray(i)
    p = np.asarray(psf, dtype=np.float64)
    cls = _class_of(img)
    nt = resolve_threads(num_threads)
    if p.ndim > img.ndim or any(
        k > 1 and k > n / 2 for k, n in zip(p.shape[::-1], img.shape[::-1], strict=False)
    ):
        raise MatlabError(
            "images:edgetaper:psfMustBeSmallerThanHalfImage",
            "EDGETAPER: PSF size must be smaller than half of the image size in any nonsingleton "
            "dimension.",
        )
    pp = p.reshape((1,) * (img.ndim - p.ndim) + p.shape)
    dt = np.float64
    I = img.astype(dt)
    alpha = np.ones(I.shape, dtype=dt)
    for ax in range(I.ndim):
        if pp.shape[ax] == 1:
            continue
        other = tuple(k for k in range(I.ndim) if k != ax)
        proj = pp.sum(axis=other) if other else pp.copy()
        n = I.shape[ax]
        P = np.fft.fft(proj, n - 1)
        b = np.abs(np.fft.ifft(np.abs(P) ** 2))
        b = b / b.max()
        b = np.concatenate([b, b[:1]])
        shp = [1] * I.ndim
        shp[ax] = n
        alpha = alpha * (1 - b).reshape(shp)
    from petakit5d import fft

    hs = fft.HalfSpectrum(I.shape, dt, threads=nt)
    blurred = hs.convolve(I, hs.otf(pp / pp.sum()))
    J = alpha * I + (1 - alpha) * blurred
    return _out_cast(J, cls)
