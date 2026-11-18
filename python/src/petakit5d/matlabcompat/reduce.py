"""MATLAB's sum, mean, cumsum, max and min with their index, sort and unique, exactly: the order of
additions, the classes and the tie, NaN and signed-zero rules of R2024b (measured, N05;
``docs/decisions/0005-sums-and-math-library.md``).

Arrays are the port's ``(z, y, x)`` C arrays (an image ``(y, x)``; a 1-D array is a row); ``axis``
is a NumPy axis or a tuple of them, None for all elements (MATLAB's ``x(:)`` or ``'all'``).
``matlab_layout=True`` takes an array in MATLAB's own form instead (``P[i-1, j-1] = M(i, j)``:
matrices of parameters, such as a list of tile sizes, and kernels that keep MATLAB's layout); its
axis k is MATLAB's dimension k + 1. Results come back in the layout of the input, reduced axes
dropped unless ``keepdims``.

MATLAB adds in one of two ways, single and double data alike (single data in single; double,
integer and logical data and outtype 'double' in double):

- **the vector rule**, for elements that are contiguous in MATLAB's memory (``x(:)``, ``'all'``,
  dimension 1, and reduced dimensions that only singletons precede): in index order up to 43
  elements; up to 88,999 in four lanes (element t in lane mod(t, 4), each lane in index order,
  added as ((l1 + l2) + l3) + l4); from 89,000 in blocks of 4,096 elements, each in four lanes
  (a last block of up to 43 elements in index order: X05), the block sums in index order. From
  89,000 elements MATLAB splits the blocks among its threads, so there its own result depends on
  its thread count; the result here is MATLAB's at one thread. ``MANY_VECTORS_FROM`` (17) or more
  vectors in an array of ``MANY_NUMEL_FROM`` (39,600) or more elements are summed in 4,096-blocks
  at any length: the same values up to 4,096 elements, other values above (measured, X05).
- **the slice rule**, for dimensions that a non-singleton one precedes (``sum(x, 2)``,
  ``sum(x, 3)``): out(i) = x(i, 1) + x(i, 2) + ... in index order, at any thread count; in chunks
  of 1,024 (each chunk in order, the chunk sums in order) from 2,048 elements along the
  dimensions, and from 1,025 in arrays of more than ``SLICE_CHUNK_NUMEL`` elements when the
  dimensions before have ``SLICE_CHUNK_INNER`` or more elements or more than one slice follows
  (measured on 120 shapes; the inner bound lies between 14 and 16).

Reduced dimensions next to each other count as one (``[1 2]``: the vector rule on ``x(:, :, k)``;
``[2 3]``: the slice rule over the ``nx * nz`` slices); runs of them are reduced one after the
other, in increasing order (``[1 3]``: the vector rule along y, then the slice rule along z).
``mean`` divides the sum by the count in double and rounds to the class (``single(double(s) /
n)``; single division differs above 2^24 elements). A sum of -0s is -0, an empty sum +0;
``omitnan`` skips NaN, and a NaN that starts a sum (a lane, a block, a chunk) counts as +0.
Every addition is one IEEE operation in the class, made by NumPy's element-wise ``add`` or a
cumulative sum (never a NumPy reduction, whose order and SIMD paths differ), so results do not
depend on the platform, the thread count or NumPy's SIMD dispatch.

Each NaN result has the NaN MATLAB gives that output (measured, task X05). An x86 addition
returns its first NaN operand, quieted, and Inf - Inf MATLAB's NaN (sign bit set), so the
payload follows MATLAB's operand order: the running sums are the first operand of their
additions, except the four-lane sum's last mod(L, 4) terms (and those of the last 4,096-block of
single data summed in double), the block sums (single and double sums of their own class) and
the slices of single data summed in double; in a block of single data summed in double l2 is
the first operand of l1 + l2. A sum of one value is a copy (a signaling NaN stays signaling);
with 'omitnan' a NaN result is MATLAB's NaN. The fast paths leave the payload to the compiler
and the CPU, so the NaN outputs alone are summed again in that order (``_exact_nans``), on every
platform; results without NaN cost one ``isnan`` pass over the output. From 89,000 elements
MATLAB's payloads, like its values, depend on its thread count: they are MATLAB's at one thread.
"""

from __future__ import annotations

import builtins
from concurrent.futures import ThreadPoolExecutor
from typing import NamedTuple

import numpy as np

from petakit5d.core import dtypes
from petakit5d.core.arrays import numpy_dtype
from petakit5d.core.axes import axis_to_dim
from petakit5d.core.threads import resolve_threads

__all__ = [
    "BLOCK",
    "BLOCKS_FROM",
    "MANY_NUMEL_FROM",
    "MANY_VECTORS_FROM",
    "SEQUENTIAL_UPTO",
    "SLICE_CHUNK",
    "SLICE_CHUNK_FROM",
    "SLICE_CHUNK_INNER",
    "SLICE_CHUNK_NUMEL",
    "MatlabError",
    "cumsum",
    "default_axis",
    "max",
    "maximum",
    "mean",
    "min",
    "minimum",
    "sort",
    "sum",
    "unique",
]

SEQUENTIAL_UPTO = 43
BLOCKS_FROM = 89000
BLOCK = 4096
SLICE_CHUNK = 1024
SLICE_CHUNK_FROM = 2048
SLICE_CHUNK_NUMEL = 20000
SLICE_CHUNK_INNER = 14

_SEGMENT_BLOCKS = 256  # blocks of 4,096 per segment of a long vector (1 Mi elements)
MANY_VECTORS_FROM = 17
MANY_NUMEL_FROM = 39600
_FEW_VECTORS = MANY_VECTORS_FROM - 1  # up to this many vectors, each is summed as one run
_QUIET = {4: np.uint32(0x00400000), 8: np.uint64(0x0008000000000000)}
_UINT = {4: np.uint32, 8: np.uint64}


class MatlabError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, identifier: str, message: str):
        super().__init__(f"{message} ({identifier})")
        self.identifier = identifier


# -- layout ------------------------------------------------------------------------------------


def _array(x) -> np.ndarray:
    a = np.asarray(x)
    if a.dtype.kind not in "biuf":
        raise TypeError(f"real numeric or logical array required, got {a.dtype}")
    return a


def _mview(a: np.ndarray, matlab_layout: bool) -> np.ndarray:
    """MATLAB's array as a view indexed as MATLAB indexes it (``m[i, j, k]`` is x(i+1, j+1, k+1)),
    at least 2-D (a 1-D array is a row)"""
    if a.ndim == 0:
        return a.reshape(1, 1)
    if a.ndim == 1:
        return a.reshape(1, -1)
    if matlab_layout:
        return a
    n = a.ndim
    return np.transpose(a, (n - 2, n - 1, *range(n - 3, -1, -1)))


def _py_axes(dims, ndim: int, matlab_layout: bool) -> tuple[int, ...]:
    """the input's axes of dimensions of MATLAB's view"""
    if ndim <= 1:
        return tuple(0 for d in dims if d == 1) if ndim == 1 else ()
    if matlab_layout:
        return tuple(dims)
    perm = (ndim - 2, ndim - 1, *range(ndim - 3, -1, -1))
    return tuple(perm[d] for d in dims)


def _mdims(axis, ndim: int, matlab_layout: bool) -> tuple[int, ...]:
    """the 0-based dimensions of MATLAB's view to reduce over"""
    vnd = builtins.max(ndim, 2)
    if axis is None:
        return tuple(range(vnd))
    axes = (axis,) if isinstance(axis, int | np.integer) else tuple(axis)
    nax = builtins.max(ndim, 1)
    out = set()
    for ax in axes:
        ax = int(ax)
        if not -nax <= ax < nax:
            raise ValueError(f"axis {ax} out of range for {ndim} dimensions")
        if ndim <= 1:
            out.add(1)
        elif matlab_layout:
            out.add(ax % ndim)
        else:
            out.add(axis_to_dim(ax, ndim) - 1)
    return tuple(sorted(out))


def _finish(rv: np.ndarray, a: np.ndarray, dims, matlab_layout: bool, keepdims: bool):
    """a result in MATLAB's view (reduced dims of size 1) in the input's layout"""
    if a.ndim == 0:
        return rv.reshape(())[()]
    if a.ndim == 1:
        r = rv.reshape(rv.shape[1:])
    elif matlab_layout:
        r = rv
    else:
        n = rv.ndim
        r = np.transpose(rv, (*range(n - 1, 1, -1), 0, 1))
    if not keepdims:
        axes = set(_py_axes(dims, a.ndim, matlab_layout))
        r = r.reshape([s for i, s in enumerate(r.shape) if i not in axes])
    r = np.ascontiguousarray(r)
    return r[()] if r.ndim == 0 else r


def default_axis(x, *, matlab_layout: bool = False) -> int:
    """The axis of MATLAB's default dimension: the first non-singleton one in MATLAB's order
    (dimension 1 when all are singletons), as ``sum(x)``, ``median(x)`` or ``max(x)`` take it."""
    a = np.asarray(x)
    if a.ndim == 0:
        return 0
    m = _mview(a, matlab_layout)
    d = next((k for k, s in enumerate(m.shape) if s != 1), 0)
    if a.ndim == 1:
        return 0
    return _py_axes((d,), a.ndim, matlab_layout)[0]


# -- classes and NaN ---------------------------------------------------------------------------


def _sum_class(a: np.ndarray, outtype: str) -> np.dtype:
    """the class the sum adds in: single data in single, other data in double ('default');
    'double'; 'native' (the input's class; integers saturate)"""
    if outtype not in ("default", "double", "native"):
        raise ValueError(f"outtype must be 'default', 'double' or 'native', got {outtype!r}")
    if a.dtype.kind in "iu" and a.dtype.itemsize == 8:
        raise TypeError("64-bit integers are not supported")
    if outtype == "double":
        return np.dtype(np.float64)
    if outtype == "native":
        if a.dtype.kind == "b":
            raise TypeError("'native' sums of logical data are not supported")
        return a.dtype
    return np.dtype(np.float32) if a.dtype == np.float32 else np.dtype(np.float64)


def _quiet(v: np.ndarray) -> np.ndarray:
    """v (float) with its NaNs quieted, as an x86 addition returns a NaN operand"""
    v = np.asarray(v)
    u = _UINT[v.dtype.itemsize]
    return np.where(np.isnan(v), (v.view(u) | _QUIET[v.dtype.itemsize]).view(v.dtype), v)


def _xadd(a: np.ndarray, b: np.ndarray, dt) -> np.ndarray:
    """a + b in class dt with x86's NaN rule on every CPU: a's NaN, else b's (quieted), else
    MATLAB's NaN (Inf - Inf)"""
    with np.errstate(invalid="ignore", over="ignore"):
        r = np.add(a, b, dtype=dt)
    nan = np.isnan(r)
    if nan.any():
        na, nb = np.isnan(a), np.isnan(b)
        r = np.where(na, _quiet(a), np.where(nb, _quiet(b), np.where(nan, dt.type(dtypes.NAN), r)))
    return r


def _fold_first(W: np.ndarray, dt) -> np.ndarray:
    """((-0 + w1) + w2) + ... of each row of W (k, L), L >= 1, the running sum the first operand
    of each addition: a NaN result has the NaN of the first addition that made one (the term's,
    quieted, or MATLAB's NaN for Inf - Inf)"""
    with np.errstate(invalid="ignore", over="ignore"):
        c = np.cumsum(W, axis=1, dtype=dt)
    s = np.array(c[:, -1], copy=True)
    rows = np.flatnonzero(np.isnan(s))
    if rows.size:
        k = np.argmax(np.isnan(c[rows]), axis=1)
        x = W[rows, k]
        s[rows] = np.where(np.isnan(x), _quiet(x), dt.type(dtypes.NAN))
    return s


def _fold_last(W: np.ndarray, dt) -> np.ndarray:
    """the same additions with each term the first operand: a NaN result has the last NaN term's
    NaN (quieted), or MATLAB's NaN (Inf - Inf) without one"""
    with np.errstate(invalid="ignore", over="ignore"):
        s = np.array(np.cumsum(W, axis=1, dtype=dt)[:, -1], copy=True)
    rows = np.flatnonzero(np.isnan(s))
    if rows.size:
        wn = np.isnan(W[rows])
        last = W.shape[1] - 1 - np.argmax(wn[:, ::-1], axis=1)
        x = W[rows, last]
        s[rows] = np.where(wn.any(axis=1), _quiet(x), dt.type(dtypes.NAN))
    return s


def _lanes_exact(W: np.ndarray, dt, tail_first: bool, swap: bool = False) -> np.ndarray:
    """the four-lane sum of the rows of W (k, L), L >= 4, in MATLAB's operand order: the lanes'
    running sums first, but the last mod(L, 4) terms first with ``tail_first``; the lanes added
    as ((l1 + l2) + l3) + l4, the running sum first (``swap``: l2 first in l1 + l2)"""
    L = W.shape[1]
    main = L - L % 4
    acc = [_fold_first(W[:, j:main:4], dt) for j in range(4)]
    for t in range(main, L):
        j = t - main
        acc[j] = _xadd(W[:, t], acc[j], dt) if tail_first else _xadd(acc[j], W[:, t], dt)
    a01 = _xadd(acc[1], acc[0], dt) if swap else _xadd(acc[0], acc[1], dt)
    return _xadd(_xadd(a01, acc[2], dt), acc[3], dt)


class _Order(NamedTuple):
    """MATLAB's operand order where NaNs meet, for one (data class, sum class)"""

    tail: bool  # the four-lane sum's last mod(L, 4) terms are the first operand
    block_tail: bool  # those of the last 4,096-block are
    blocks_term_first: bool  # each block sum is the first operand of its addition
    block_swap: bool  # in a block, l2 is the first operand of l1 + l2
    slice_term_first: bool  # in the slice rule, each slice is the first operand


# measured in R2024b at one thread (task X05's probes; from 89,000 elements MATLAB's own
# payloads depend on its thread count, as its values do)
_NAN_ORDER = {
    (np.dtype(np.float32), np.dtype(np.float32)): _Order(True, False, True, False, False),
    (np.dtype(np.float32), np.dtype(np.float64)): _Order(True, True, False, True, True),
    (np.dtype(np.float64), np.dtype(np.float64)): _Order(True, False, True, False, False),
}


def _vector_exact(W: np.ndarray, dt, o: _Order, many: bool = False) -> np.ndarray:
    """the vector rule over the rows of W (k, L) in MATLAB's operand order; ``many``: the order
    of MATLAB's sums of many vectors (``MANY_VECTORS_FROM``), its 4,096-blocks at any length"""
    k, L = W.shape
    if L == 1:
        return W[:, 0].copy()  # (a copy: a NaN not quieted)
    if L <= SEQUENTIAL_UPTO:
        return _fold_first(W, dt)
    if L < BLOCKS_FROM and not many:
        return _lanes_exact(W, dt, o.tail)
    nf = L // BLOCK
    full = W[:, : nf * BLOCK].reshape(k * nf, BLOCK)
    bs = _lanes_exact(full, dt, False, o.block_swap).reshape(k, nf)
    if L % BLOCK:
        W1 = W[:, nf * BLOCK :]
        if W1.shape[1] > SEQUENTIAL_UPTO:
            last = _lanes_exact(W1, dt, o.block_tail, o.block_swap)
        else:
            last = _fold_first(W1, dt)
        bs = np.concatenate([bs, last[:, None]], axis=1)
    return _fold_last(bs, dt) if o.blocks_term_first else _fold_first(bs, dt)


def _slice_exact(W: np.ndarray, dt, chunked: bool, term_first: bool) -> np.ndarray:
    """the slice rule over the rows of W (k, R) in MATLAB's operand order: each sum from its first
    slice (a copy), then the running sum the first operand of each addition (``term_first``: the
    slice); with ``chunked`` the chunk sums added in order, the running total first"""
    fold = _fold_last if term_first else _fold_first
    R = W.shape[1]
    if R == 1:
        return W[:, 0].astype(dt)
    if not chunked or R <= SLICE_CHUNK:
        return fold(W, dt)
    parts = []
    for c0 in range(0, R, SLICE_CHUNK):
        P = W[:, c0 : c0 + SLICE_CHUNK]
        parts.append(P[:, 0].astype(dt) if P.shape[1] == 1 else fold(P, dt))
    return _fold_first(np.stack(parts, axis=1), dt)


def _exact_nans(m: np.ndarray, dims, dt, nan: np.ndarray) -> np.ndarray:
    """the NaN outputs (``nan``, MATLAB's view) of the sum over dims of float view m, summed
    again in MATLAB's operand order, which decides their payloads: each output's values gathered
    alone, the dims reduced group by group by the rule ``_sum_core`` takes"""
    kept = [d for d in range(m.ndim) if d not in dims]
    idx = np.nonzero(nan)
    mm = np.moveaxis(m, kept, list(range(len(kept))))
    cur = mm[tuple(idx[d] for d in kept)] if kept else mm[None]
    cur = cur.astype(dt, copy=False)
    order = _NAN_ORDER.get((m.dtype, dt), _NAN_ORDER[(np.dtype(np.float64),) * 2])
    shape = list(m.shape)
    for d0, d1 in _groups(dims):
        ng = d1 - d0 + 1
        red = int(np.prod(cur.shape[1 : 1 + ng], dtype=np.int64))
        inner = int(np.prod(shape[:d0], dtype=np.int64))
        outer = int(np.prod(shape[d1 + 1 :], dtype=np.int64))
        # the group's elements last, in MATLAB's order
        t = np.moveaxis(cur, list(range(1, 1 + ng)), list(range(-1, -ng - 1, -1)))
        rest = t.shape[: t.ndim - ng]
        W = np.ascontiguousarray(t).reshape(-1, red)
        if inner == 1:
            numel = int(np.prod(shape, dtype=np.int64))
            many = outer >= MANY_VECTORS_FROM and numel >= MANY_NUMEL_FROM
            s = _vector_exact(W, dt, order, many)
        else:
            chunked = _slice_chunked(red, inner, outer, m.size)
            s = _slice_exact(W, dt, chunked, order.slice_term_first)
        cur = s.reshape(rest)
        shape[d0 : d1 + 1] = [1] * ng
        # the next dims add the sums so far, of the sum's class
        order = _NAN_ORDER.get((dt, dt), order)
    return cur.reshape(-1)


def _nan_fix(r: np.ndarray, m: np.ndarray, dims, dt, omit: bool) -> np.ndarray:
    """r (the sum over dims of view m, MATLAB's view) with each NaN output as MATLAB gives it:
    with 'omitnan' MATLAB's NaN (only Inf - Inf makes one), else the NaN of the first addition
    that made one in MATLAB's operand order"""
    if dt.kind != "f":
        return r
    nan = np.isnan(r)
    if not nan.any():
        return r
    r = np.array(r, dtype=dt, copy=True)
    if omit or m.dtype.kind != "f":
        r[nan] = dt.type(dtypes.NAN)
    else:
        r[nan] = _exact_nans(m, dims, dt, nan)
    return r


# -- the vector rule ---------------------------------------------------------------------------


def _combine4(l0, l1, l2, l3, dt):
    return np.add(np.add(np.add(l0, l1, dtype=dt), l2, dtype=dt), l3, dtype=dt)


def _pad(w: np.ndarray, mult: int, dt) -> np.ndarray:
    p = -w.shape[-1] % mult
    if not p:
        return w
    return np.concatenate([w, np.full(w.shape[:-1] + (p,), -0.0, dtype=dt)], axis=-1)


def _omit_starts(w: np.ndarray, period: int, nstart: int) -> None:
    """in place: NaN as -0, and as +0 where a lane starts (the first ``nstart`` elements of every
    ``period``)"""
    nan = np.isnan(w)
    if not nan.any():
        return
    w[nan] = -0.0
    n = w.shape[-1]
    st = (np.arange(0, n, period)[:, None] + np.arange(nstart)[None, :]).ravel()
    st = st[st < n]
    s = w[..., st]
    s[nan[..., st]] = 0.0
    w[..., st] = s


def _rows_rule(w: np.ndarray, dt, omit: bool, nt: int) -> np.ndarray:
    """the vector rule over the rows of w (k, L), a contiguous copy in class dt (changed)"""
    k, L = w.shape
    if L == 0:
        return np.zeros(k, dtype=dt)
    if L <= SEQUENTIAL_UPTO:
        if omit:
            _omit_starts(w, L, 1)
        return np.cumsum(w, axis=1, dtype=dt)[:, -1]
    if L < BLOCKS_FROM:
        if omit:
            _omit_starts(w, L, 4)
        r = np.cumsum(_pad(w, 4, dt).reshape(k, -1, 4), axis=1, dtype=dt)[:, -1, :]
        return _combine4(r[:, 0], r[:, 1], r[:, 2], r[:, 3], dt)
    nb = -(-L // BLOCK)
    bs = np.empty((k, nb), dtype=dt)
    jobs = [(b0, builtins.min(b0 + _SEGMENT_BLOCKS, nb)) for b0 in range(0, nb, _SEGMENT_BLOCKS)]

    def run(job):
        b0, b1 = job
        seg = w[:, b0 * BLOCK : builtins.min(b1 * BLOCK, L)]
        if omit:
            seg = seg.copy()
            _omit_starts(seg, BLOCK, 4)
        bs[:, b0:b1] = _block_sums(_pad(seg, BLOCK, dt), dt)

    _parallel(run, jobs, nt)
    if L - (nb - 1) * BLOCK <= SEQUENTIAL_UPTO:
        bs[:, -1] = _seq_rows(w[:, (nb - 1) * BLOCK :], dt, omit)
    return np.cumsum(bs, axis=1, dtype=dt)[:, -1]


def _short_last_block(n: int) -> bool:
    """whether a sum of n values in 4,096-blocks ends with a block of 4 to 43 values, which
    MATLAB adds in index order (``_seq_rows``; four lanes give other bits)"""
    return n >= BLOCKS_FROM and 4 <= n % BLOCK <= SEQUENTIAL_UPTO


def _seq_rows(w: np.ndarray, dt, omit: bool) -> np.ndarray:
    """the sums of the rows of w (k, m) in index order (the last block of a block sum)"""
    w = np.array(w, dtype=dt, copy=True)
    if omit:
        _omit_starts(w, w.shape[1], 1)
    return np.cumsum(w, axis=1, dtype=dt)[:, -1]


def _block_sums(w: np.ndarray, dt) -> np.ndarray:
    """the sums of the 4,096-blocks of the rows of w (k, nb * 4096), each in four lanes"""
    k = w.shape[0]
    lanes = np.cumsum(w.reshape(k, -1, BLOCK // 4, 4), axis=2, dtype=dt)[:, :, -1, :]
    return _combine4(lanes[..., 0], lanes[..., 1], lanes[..., 2], lanes[..., 3], dt)


def _parallel(fn, jobs, nt: int) -> None:
    if nt <= 1 or len(jobs) <= 1:
        for j in jobs:
            fn(j)
        return
    with ThreadPoolExecutor(max_workers=builtins.min(nt, len(jobs))) as ex:
        list(ex.map(fn, jobs))


def _fortran_fill(G: np.ndarray, t0: int, t1: int, out: np.ndarray) -> None:
    """out[:] = the elements t0 .. t1-1 of G in Fortran order (MATLAB's linear order)"""
    if t1 <= t0:
        return
    if G.ndim <= 1:
        out[:] = G.reshape(-1)[t0:t1]
        return
    if G.ndim == 2:
        g0 = G.shape[0]
        c0, c1 = t0 // g0, (t1 - 1) // g0
        blk = np.ascontiguousarray(G[:, c0 : c1 + 1].T).reshape(-1)
        out[:] = blk[t0 - c0 * g0 : t1 - c0 * g0]
        return
    P = int(np.prod(G.shape[:-1], dtype=np.int64))
    pos = 0
    for s in range(t0 // P, (t1 - 1) // P + 1):
        lo, hi = builtins.max(t0 - s * P, 0), builtins.min(t1 - s * P, P)
        _fortran_fill(G[..., s], lo, hi, out[pos : pos + hi - lo])
        pos += hi - lo


def _fortran_flat(G: np.ndarray, dt) -> np.ndarray:
    """the elements of G in Fortran order, contiguous, in class dt"""
    if G.flags.f_contiguous:
        return np.array(G.reshape(-1, order="F"), dtype=dt, copy=True)
    out = np.empty(G.size, dtype=G.dtype)
    _fortran_fill(G, 0, G.size, out)
    return out.astype(dt, copy=False)


def _long_rule(G: np.ndarray, dt, omit: bool, nt: int):
    """the vector rule (blocks of 4,096) over all elements of G in Fortran order, read in
    segments of 256 blocks (in threads; the block sums are added in order afterwards)"""
    L = G.size
    nb = -(-L // BLOCK)
    bs = np.empty(nb, dtype=dt)
    jobs = [(b0, builtins.min(b0 + _SEGMENT_BLOCKS, nb)) for b0 in range(0, nb, _SEGMENT_BLOCKS)]

    def run(job):
        b0, b1 = job
        t0, t1 = b0 * BLOCK, builtins.min(b1 * BLOCK, L)
        seg = np.empty(t1 - t0, dtype=G.dtype)
        _fortran_fill(G, t0, t1, seg)
        w = seg.astype(dt, copy=False).reshape(1, -1)
        if omit:
            w = w.copy() if w.base is seg else w
            _omit_starts(w, BLOCK, 4)
        bs[b0:b1] = _block_sums(_pad(w, BLOCK, dt), dt)[0]

    _parallel(run, jobs, nt)
    m = L - (nb - 1) * BLOCK
    if m <= SEQUENTIAL_UPTO:
        seg = np.empty(m, dtype=G.dtype)
        _fortran_fill(G, L - m, L, seg)
        bs[-1] = _seq_rows(seg.reshape(1, -1), dt, omit)[0]
    return np.cumsum(bs, dtype=dt)[-1]


def _add_slice(acc: np.ndarray, v: np.ndarray, dt, omit: bool) -> None:
    if omit:
        v = np.where(np.isnan(v), dt.type(-0.0), v)
    np.add(acc, v, out=acc, dtype=dt, casting="unsafe")


def _start(v: np.ndarray, dt, omit: bool) -> np.ndarray:
    acc = np.array(v, dtype=dt, copy=True)
    if omit:
        acc[np.isnan(acc)] = 0.0
    return acc


def _many_rule(V: np.ndarray, dt, omit: bool) -> np.ndarray:
    """the vector rule over the last axis of V (..., L), L < 89,000, for many vectors at once:
    one vectorized addition per element of the vectors"""
    L = V.shape[-1]
    if L <= SEQUENTIAL_UPTO:
        acc = _start(V[..., 0], dt, omit)
        for j in range(1, L):
            _add_slice(acc, V[..., j], dt, omit)
        return acc
    lanes = []
    for k in range(4):
        acc = _start(V[..., k], dt, omit)
        for j in range(k + 4, L, 4):
            _add_slice(acc, V[..., j], dt, omit)
        lanes.append(acc)
    return _combine4(*lanes, dt)


def _many_blocks(V: np.ndarray, dt, omit: bool) -> np.ndarray:
    """the block rule over the last axis of V (..., L) for many vectors: blocks of 4,096, each
    by ``_many_rule`` (a last block of up to 43 values in index order), the block sums added in
    order (MATLAB sums 17 or more vectors this way at any length)"""
    total = None
    for b0 in range(0, V.shape[-1], BLOCK):
        s = _many_rule(V[..., b0 : b0 + BLOCK], dt, omit)
        total = s if total is None else np.add(total, s, dtype=dt)
    return total


def _vector_rule(m: np.ndarray, d0: int, d1: int, dt, omit: bool, nt: int) -> np.ndarray:
    """the vector rule over dims d0..d1 of view m (the dims before d0 are singletons); the
    result with those dims of size 1"""
    S = m.shape
    ng = d1 - d0 + 1
    red = int(np.prod(S[d0 : d1 + 1], dtype=np.int64))
    oshape = S[d1 + 1 :]
    outer = int(np.prod(oshape, dtype=np.int64))
    out_shape = S[:d0] + (1,) * ng + oshape
    if red == 0 or outer == 0:
        return np.zeros(out_shape, dtype=dt)
    g = m.reshape(S[d0:])  # the leading singletons dropped (a view)
    if red >= BLOCKS_FROM:
        res = np.empty(outer, dtype=dt)
        for i, idx in enumerate(np.ndindex(*oshape[::-1])):
            res[i] = _long_rule(g[(Ellipsis, *idx[::-1])], dt, omit, nt)
        return res.reshape(oshape[::-1]).transpose().reshape(out_shape)
    if outer <= _FEW_VECTORS:
        w = np.empty((outer, red), dtype=dt)
        for i, idx in enumerate(np.ndindex(*oshape[::-1])):
            w[i] = _fortran_flat(g[(Ellipsis, *idx[::-1])], dt)
        res = _rows_rule(w, dt, omit, nt)
        return res.reshape(oshape[::-1]).transpose().reshape(out_shape)
    # many vectors: the group's elements last, in Fortran order (a view for one dimension);
    # MATLAB sums 17 or more in 4,096-blocks at any length (MANY_VECTORS_FROM)
    t = np.moveaxis(g, list(range(ng)), list(range(-1, -ng - 1, -1)))
    V = t.reshape(t.shape[: t.ndim - ng] + (red,))
    if red > BLOCK:
        return _many_blocks(V, dt, omit).reshape(out_shape)
    return _many_rule(V, dt, omit).reshape(out_shape)


# -- the slice rule ----------------------------------------------------------------------------


def _slice_chunked(red: int, inner: int, outer: int, numel: int) -> bool:
    if red >= SLICE_CHUNK_FROM:
        return True
    return (
        red > SLICE_CHUNK
        and numel > SLICE_CHUNK_NUMEL
        and (inner >= SLICE_CHUNK_INNER or outer > 1)
    )


def _slice_rule(m: np.ndarray, d0: int, d1: int, dt, omit: bool, chunked: bool) -> np.ndarray:
    """out(i, o) = the sum of x(i, r, o) over the dims d0..d1 of view m, slice by slice in
    MATLAB's order (in chunks of 1,024 slices when chunked); the dims kept with size 1"""
    S = m.shape
    gshape = S[d0 : d1 + 1]
    red = int(np.prod(gshape, dtype=np.int64))
    out_shape = S[:d0] + (1,) * len(gshape) + S[d1 + 1 :]
    if red == 0:
        return np.zeros(out_shape, dtype=dt)
    pre = (slice(None),) * d0
    chunk = SLICE_CHUNK if chunked else red
    total = part = None
    for r, gi in enumerate(np.ndindex(*gshape[::-1])):
        v = m[pre + gi[::-1]]
        if r % chunk == 0:
            if part is not None:
                total = part if total is None else np.add(total, part, dtype=dt)
            part = _start(v, dt, omit)
        else:
            _add_slice(part, v, dt, omit)
    total = part if total is None else np.add(total, part, dtype=dt)
    return total.reshape(out_shape)


# -- sum and mean ------------------------------------------------------------------------------


def _groups(dims) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for d in dims:
        if out and d == out[-1][1] + 1:
            out[-1][1] = d
        else:
            out.append([d, d])
    return [(a, b) for a, b in out]


def _native_int_sum(m: np.ndarray, dims, dt) -> np.ndarray:
    """'native' sums of integers: saturating at every addition, in index order"""
    info = np.iinfo(dt)
    r = m.astype(np.int64)
    for d in sorted(dims):
        if r.shape[d] == 0:
            r = np.zeros(r.shape[:d] + (1,) + r.shape[d + 1 :], np.int64)
            continue
        acc = np.take(r, 0, axis=d)
        for j in range(1, r.shape[d]):
            acc = np.clip(acc + np.take(r, j, axis=d), info.min, info.max)
        r = np.expand_dims(acc, d)
    return r.astype(dt)


def _from_view(v: np.ndarray, ndim: int, matlab_layout: bool) -> np.ndarray:
    """an array of MATLAB's view (reduced dims of size 1) in the layout of an input of ndim
    dimensions (the inverse of ``_mview``)"""
    if matlab_layout or ndim < 2:
        return v
    perm = (ndim - 2, ndim - 1, *range(ndim - 3, -1, -1))
    return np.transpose(v, np.argsort(perm))


def _kernel_sum(a, m, dims, dt, nt: int, matlab_layout: bool, center=None):
    """the sum by K08's kernels (``kernels.ops``) where they cover it, in MATLAB's view with the
    reduced dims of size 1, else None: all elements (the vector rule), and one dimension of the
    slice rule on a C- or Fortran-contiguous array; single sums of single data, double sums of
    any data. ``center`` (MATLAB's view, of a's class): the sum of ``abs(x - center).^2`` (var's).
    The same additions in the same order as the NumPy code below, which stays the reference."""
    from petakit5d.kernels import ops

    if not ops.extension_available() or (a.dtype.kind in "iu" and a.dtype.itemsize == 8):
        return None
    single = dt == np.float32
    if (single and a.dtype != np.float32) or (center is not None and dt != a.dtype):
        return None
    if len(dims) == m.ndim:
        if _short_last_block(m.size):
            return None  # (the kernels add a short last block in four lanes)
        c = None if center is None else float(np.asarray(center).reshape(-1)[0])
        s = ops.sum_ordered(a, single=single, matlab_layout=matlab_layout, center=c, num_threads=nt)
        return None if s is None else np.full((1,) * m.ndim, s, dtype=dt)
    if len(dims) != 1 or a.ndim < 2:
        return None
    d = dims[0]
    inner = int(np.prod(m.shape[:d], dtype=np.int64))
    if inner == 1:
        return None  # the vector rule on many vectors (NumPy)
    red = m.shape[d]
    outer = int(np.prod(m.shape[d + 1 :], dtype=np.int64))
    chunked = _slice_chunked(red, inner, outer, m.size)
    ax = _py_axes((d,), a.ndim, matlab_layout)[0]
    if a.flags.c_contiguous:
        A, axA = a, ax
    elif a.flags.f_contiguous:
        A, axA = a.T, a.ndim - 1 - ax
    else:
        return None
    O = int(np.prod(A.shape[:axA], dtype=np.int64))
    I = int(np.prod(A.shape[axA + 1 :], dtype=np.int64))  # noqa: E741
    cen = None
    if center is not None:
        cA = _from_view(np.asarray(center), a.ndim, matlab_layout)
        cA = cA if A is a else cA.T
        cen = np.ascontiguousarray(cA, dtype=a.dtype).reshape(O, I)
    out = ops.sum_slices(
        A.reshape(O, A.shape[axA], I), single=single, chunked=chunked, center=cen, num_threads=nt
    )
    if out is None:
        return None
    rA = out.reshape(A.shape[:axA] + (1,) + A.shape[axA + 1 :])
    return _mview(rA if A is a else rA.T, matlab_layout)


def _sum_core(
    a: np.ndarray, m: np.ndarray, dims, dt, omit: bool, nt: int, matlab_layout: bool
) -> np.ndarray:
    """the sum over dims of view m of array a, in MATLAB's view with the dims of size 1"""
    omit = bool(omit) and m.dtype.kind == "f"
    if dt.kind in "iu":
        return _native_int_sum(m, dims, dt)
    full = len(dims) == m.ndim
    if (
        dt == np.float64
        and not omit
        and full
        and m.size >= BLOCKS_FROM
        and not _short_last_block(m.size)
    ):
        # sum(x(:), 'double') in MATLAB's one-thread order: decon_sum_double's kernel (K08)
        from petakit5d.kernels import ops

        src = a if a.dtype.kind in "biuf" else a.astype(np.float64)
        s = ops.sum_double(src, matlab_layout=matlab_layout, num_threads=nt)
        # (a NaN: decon_sum_double's payload, which isn't sum's for single data)
        return _nan_fix(np.full((1,) * m.ndim, s), m, dims, dt, False)
    if not omit and dt.kind == "f":
        r = _kernel_sum(a, m, dims, dt, nt, matlab_layout)
        if r is not None:
            return _nan_fix(r, m, dims, dt, False)
    numel = m.size
    r = m
    for d0, d1 in _groups(dims):
        inner = int(np.prod(r.shape[:d0], dtype=np.int64))
        if inner == 1:
            r = _vector_rule(r, d0, d1, dt, omit, nt)
        else:
            red = int(np.prod(r.shape[d0 : d1 + 1], dtype=np.int64))
            outer = int(np.prod(r.shape[d1 + 1 :], dtype=np.int64))
            r = _slice_rule(r, d0, d1, dt, omit, _slice_chunked(red, inner, outer, numel))
    return _nan_fix(np.asarray(r, dtype=dt), m, dims, dt, omit)


def sum(  # noqa: A001
    x,
    axis=None,
    *,
    outtype: str = "default",
    omitnan: bool = False,
    keepdims: bool = False,
    matlab_layout: bool = False,
    num_threads: int | None = None,
):
    """MATLAB's ``sum(x, dim, outtype, nanflag)``, added in MATLAB's order at one thread (module
    docstring), whatever ``num_threads`` (threads only read in parallel).

    ``outtype``: 'default' (single for single data, double otherwise), 'double' or 'native'
    (integers saturate at every addition); ``omitnan``: MATLAB's 'omitnan'.
    """
    a = _array(x)
    dt = _sum_class(a, outtype)
    dims = _mdims(axis, a.ndim, matlab_layout)
    m = _mview(a, matlab_layout)
    r = _sum_core(a, m, dims, dt, omitnan, resolve_threads(num_threads), matlab_layout)
    return _finish(r, a, dims, matlab_layout, keepdims)


def _count(m: np.ndarray, dims, omit: bool):
    if omit and m.dtype.kind == "f":
        return np.count_nonzero(~np.isnan(m), axis=tuple(dims), keepdims=True).astype(np.float64)
    return float(np.prod([m.shape[d] for d in dims], dtype=np.int64))


def _mean_view(a, m, dims, outtype, omit, nt, matlab_layout):
    """the mean in MATLAB's view, in double, and the class to round it to"""
    if outtype not in ("default", "double", "native"):
        raise ValueError(f"outtype must be 'default', 'double' or 'native', got {outtype!r}")
    sum_out = "double" if a.dtype.kind in "biu" else outtype
    s = _sum_core(a, m, dims, _sum_class(a, sum_out), omit, nt, matlab_layout)
    n = _count(m, dims, omit)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = s.astype(np.float64) / n
    if np.any(np.asarray(n) == 0):
        q = np.where(np.asarray(n) == 0, dtypes.NAN, q)
    if a.dtype.kind in "iu" and outtype == "native":
        return q, a.dtype.name
    if a.dtype.kind in "biu" or outtype == "double":
        return q, "double"
    return q, a.dtype.name


def mean(
    x,
    axis=None,
    *,
    outtype: str = "default",
    omitnan: bool = False,
    keepdims: bool = False,
    matlab_layout: bool = False,
    num_threads: int | None = None,
):
    """MATLAB's ``mean(x, dim, outtype, nanflag)``: ``sum``'s sum divided by the count in double,
    rounded to the class ('default': single for single data, double otherwise; 'double';
    'native': the input's class, integers rounded half away from zero and saturated). The mean of
    nothing is NaN."""
    a = _array(x)
    dims = _mdims(axis, a.ndim, matlab_layout)
    m = _mview(a, matlab_layout)
    q, cls = _mean_view(a, m, dims, outtype, omitnan, resolve_threads(num_threads), matlab_layout)
    r = np.asarray(dtypes.cast(q, cls)) if cls in dtypes.INT_CLASSES else q.astype(numpy_dtype(cls))
    return _finish(r, a, dims, matlab_layout, keepdims)


# -- cumsum ------------------------------------------------------------------------------------


def cumsum(x, axis: int, *, omitnan: bool = False, matlab_layout: bool = False) -> np.ndarray:
    """MATLAB's ``cumsum(x, dim)``: each element the previous result plus itself, in index order
    along ``axis`` (the same computation in either layout), in the class (single for single data,
    double for double and logical; integers saturate at every step). ``omitnan``: NaN adds
    nothing (+0 at the start). A NaN has MATLAB's payload (``_cumsum_nans``), which depends on
    whether ``axis`` is MATLAB's dimension 1: ``matlab_layout`` as in ``sum``."""
    a = _array(x)
    if a.ndim == 0:
        a = a.reshape(1)
    ax = int(axis) % a.ndim
    if a.dtype.kind == "b":
        return np.cumsum(a, axis=ax, dtype=np.float64)
    if a.dtype.kind in "iu":
        if a.dtype.itemsize == 8:
            raise TypeError("64-bit integers are not supported")
        info = np.iinfo(a.dtype)
        w = np.moveaxis(a, ax, 0).astype(np.int64)
        out = np.empty(w.shape, dtype=a.dtype)
        acc = None
        for j in range(w.shape[0]):
            acc = w[j] if acc is None else np.clip(acc + w[j], info.min, info.max)
            out[j] = acc
        return np.ascontiguousarray(np.moveaxis(out, 0, ax))
    w = a
    if omitnan:
        nan = np.isnan(a)
        if nan.any():
            w = np.where(nan, a.dtype.type(-0.0), a)
            lead = np.moveaxis(w, ax, 0)
            first = np.moveaxis(nan, ax, 0)[0]
            lead[0, ...] = np.where(first, a.dtype.type(0.0), lead[0, ...])
    with np.errstate(invalid="ignore", over="ignore"):
        r = np.cumsum(w, axis=ax, dtype=a.dtype)
    # (a NaN stays NaN to the end of the axis: the last slice tells)
    if r.shape[ax] and np.isnan(np.take(r, -1, axis=ax)).any():
        d = _mdims(ax, a.ndim, matlab_layout)[0]
        inner = int(np.prod(_mview(a, matlab_layout).shape[:d], dtype=np.int64))
        _cumsum_nans(r, a, ax, omitnan, inner > 1)
    return r


def _cumsum_nans(r: np.ndarray, a: np.ndarray, ax: int, omit: bool, term_first: bool) -> None:
    """in place: each NaN of cumsum's r as MATLAB gives it (measured), on the lines along the
    axis that end in NaN: the first element a copy of x(1) (a NaN not quieted); along MATLAB's
    dimension 1 (or a row) the running sum is the first operand of each addition, so a NaN keeps
    the first NaN that came, a term's (quieted) or MATLAB's NaN for Inf - Inf; along a later
    dimension (``term_first``) each term is, so a NaN is the last NaN term so far, or MATLAB's
    NaN without one; MATLAB's NaN with 'omitnan'"""
    if r.ndim == 1:
        _cumsum_nans(r[None], a[None], 1, omit, term_first)  # (views: written in place)
        return
    rv_all, av_all = np.moveaxis(r, ax, -1), np.moveaxis(a, ax, -1)
    lines = np.nonzero(np.isnan(rv_all[..., -1]))
    rv, av = rv_all[lines], av_all[lines]  # (lines, L) copies
    nv = np.isnan(rv)
    if omit:
        rv_all[lines] = np.where(nv, a.dtype.type(dtypes.NAN), rv)
        return
    if term_first:
        last = np.where(np.isnan(av), np.arange(av.shape[-1]), -1)
        k = np.maximum.accumulate(last, axis=-1)
        x = np.take_along_axis(av, np.maximum(k, 0), axis=-1)
        pay = np.where(k >= 0, _quiet(x), a.dtype.type(dtypes.NAN))
    else:
        k = np.argmax(nv, axis=-1)[..., None]
        x = np.take_along_axis(av, k, axis=-1)
        pay = np.where(np.isnan(x), _quiet(x), a.dtype.type(dtypes.NAN))
    rv = np.where(nv, pay, rv)
    rv[..., 0] = av[..., 0]  # (a copy)
    rv_all[lines] = rv


# -- max and min -------------------------------------------------------------------------------


def _first_true_matlab(mask: np.ndarray, matlab_layout: bool) -> int:
    """the 0-based MATLAB linear index of the first True of mask (-1 if none)"""
    return _nth_true_matlab(mask, 0, matlab_layout)


def _nth_true_matlab(mask: np.ndarray, j: int, matlab_layout: bool) -> int:
    """the 0-based MATLAB linear index of the j-th (0-based) True of mask (-1 if none), the
    planes of a (..., y, x) array scanned in MATLAB's order (counts per plane, then per column)"""
    if mask.ndim <= 1 or matlab_layout:
        v = mask.reshape(-1, order="F") if mask.ndim > 1 else mask.reshape(-1)
        k = np.flatnonzero(v)
        return int(k[j]) if j < k.size else -1
    ny, nx = mask.shape[-2:]
    planes = mask.reshape(-1, ny, nx)
    for p in range(planes.shape[0]):
        pl = planes[p]
        cols = np.count_nonzero(pl, axis=0)
        c = int(cols.sum())
        if j >= c:
            j -= c
            continue
        cum = np.cumsum(cols)
        xi = int(np.searchsorted(cum, j, side="right"))
        j -= int(cum[xi - 1]) if xi else 0
        yi = int(np.flatnonzero(pl[:, xi])[j])
        return yi + ny * (xi + nx * p)
    return -1


def _index_tuple(k: int, a: np.ndarray, matlab_layout: bool) -> tuple:
    """the index into a of its element k in MATLAB's linear order"""
    if a.ndim == 0:
        return ()
    if a.ndim == 1:
        return (k,)
    if matlab_layout:
        return tuple(int(i) for i in np.unravel_index(k, a.shape, order="F"))
    ny, nx = a.shape[-2:]
    y, r = k % ny, k // ny
    xi, p = r % nx, r // nx
    lead = np.unravel_index(p, a.shape[:-2]) if a.ndim > 2 else ()
    return (*(int(i) for i in lead), int(y), int(xi))


def _extreme_all(a: np.ndarray, omitnan: bool, greater: bool, matlab_layout: bool):
    red = np.fmax if greater else np.fmin
    if a.dtype.kind == "f":
        nan = np.isnan(a)
        if nan.any() and (not omitnan or nan.all()):
            k = _first_true_matlab(nan, matlab_layout)
        else:
            v = red.reduce(a, axis=None)
            k = _first_true_matlab(a == v, matlab_layout)
    else:
        v = np.max(a) if greater else np.min(a)
        k = _first_true_matlab(a == v, matlab_layout)
    idx = _index_tuple(k, a, matlab_layout)
    return a[idx], idx


def _extreme_axis(w: np.ndarray, omitnan: bool, greater: bool) -> np.ndarray:
    """the index along the last axis: the first of the equal extremes, NaN ignored unless all are
    NaN (then the first NaN); with omitnan False the first NaN"""
    if w.dtype.kind != "f":
        v = np.max(w, axis=-1, keepdims=True) if greater else np.min(w, axis=-1, keepdims=True)
        return np.argmax(w == v, axis=-1)
    red = np.fmax if greater else np.fmin
    v = red.reduce(w, axis=-1, keepdims=True)
    nan = np.isnan(w)
    k = np.argmax(w == v, axis=-1)
    use_nan = np.isnan(v[..., 0]) | (nan.any(axis=-1) & (not omitnan))
    return np.where(use_nan, np.argmax(nan, axis=-1), k)


def _extreme(x, axis, omitnan, keepdims, matlab_layout, greater):
    a = _array(x)
    if a.size == 0:
        raise ValueError("max and min of an empty array: MATLAB returns []")
    if axis is None:
        return _extreme_all(a, omitnan, greater, matlab_layout)
    if not isinstance(axis, int | np.integer):
        raise TypeError("max and min take one axis, or None for all elements")
    ax = int(axis) % builtins.max(a.ndim, 1)
    w = np.moveaxis(a, ax, -1) if a.ndim else a.reshape(1)
    k = _extreme_axis(w, omitnan, greater)
    val = np.take_along_axis(w, k[..., None], axis=-1)[..., 0]
    if keepdims:
        val, k = np.expand_dims(val, ax), np.expand_dims(k, ax)
    return (val[()] if val.ndim == 0 else val), (k[()] if k.ndim == 0 else k)


def max(  # noqa: A001
    x, axis=None, *, omitnan: bool = True, keepdims: bool = False, matlab_layout: bool = False
):
    """MATLAB's ``[m, i] = max(x, [], dim)``: the largest value and its 0-based index along
    ``axis``; for ``axis=None`` the index tuple into ``x`` of the first maximum in MATLAB's
    linear order (``max(x(:))``). Equal values: the first (``max([0 -0])`` is +0, ``max([-0 0])``
    -0). NaN is ignored unless all are NaN (MATLAB's default; ``omitnan=False``: 'includenan')."""
    return _extreme(x, axis, omitnan, keepdims, matlab_layout, True)


def min(  # noqa: A001
    x, axis=None, *, omitnan: bool = True, keepdims: bool = False, matlab_layout: bool = False
):
    """MATLAB's ``[m, i] = min(x, [], dim)`` (as ``max``)."""
    return _extreme(x, axis, omitnan, keepdims, matlab_layout, False)


def _pairwise(a, b, greater: bool, omitnan: bool):
    A, B = np.asarray(a), np.asarray(b)
    cls = dtypes.result_class(A, B)
    dt = numpy_dtype(cls)
    if dt.kind == "f":
        Af, Bf = A.astype(dt, copy=False), B.astype(dt, copy=False)
    else:
        Af = A if A.dtype == dt else np.asarray(dtypes.cast(A, cls))
        Bf = B if B.dtype == dt else np.asarray(dtypes.cast(B, cls))
    Af, Bf = np.broadcast_arrays(Af, Bf)
    pick_b = (Bf > Af) if greater else (Bf < Af)
    if dt.kind == "f":
        na, nb = np.isnan(Af), np.isnan(Bf)
        pick_b = (pick_b | (na & ~nb)) if omitnan else ((pick_b & ~na) | (nb & ~na))
    r = np.where(pick_b, Bf, Af).astype(dt)
    return r[()] if r.ndim == 0 else r


def maximum(a, b, *, omitnan: bool = True):
    """MATLAB's ``max(a, b)``: element-wise in MATLAB's class of ``a`` and ``b`` (an integer and a
    double scalar: the integer class, the double converted first); on a tie ``a``
    (``max(0, -0)`` is 0, ``max(-0, 0)`` -0); NaN gives the other value (``omitnan=False``: NaN)."""
    return _pairwise(a, b, True, omitnan)


def minimum(a, b, *, omitnan: bool = True):
    """MATLAB's ``min(a, b)`` (as ``maximum``)."""
    return _pairwise(a, b, False, omitnan)


# -- sort and unique ---------------------------------------------------------------------------


def sort(x, axis: int = -1, *, descend: bool = False):
    """MATLAB's ``[s, i] = sort(x, dim, 'ascend' | 'descend')`` along ``axis``: stable (equal
    values, -0 and +0 among them, keep their order), NaN last ascending and first descending;
    0-based indices."""
    a = _array(x)
    if a.ndim == 0:
        return a.copy(), np.zeros((), np.intp)
    ax = int(axis) % a.ndim
    w = np.moveaxis(a, ax, -1)
    if descend:
        key = -w if w.dtype.kind == "f" else (~w if w.dtype.kind == "b" else np.invert(w))
        idx = np.argsort(key, axis=-1, kind="stable")
        if a.dtype.kind == "f":
            nan = np.take_along_axis(np.isnan(w), idx, axis=-1)
            idx = np.take_along_axis(idx, np.argsort(~nan, axis=-1, kind="stable"), axis=-1)
    else:
        idx = np.argsort(w, axis=-1, kind="stable")
    s = np.take_along_axis(w, idx, axis=-1)
    return np.ascontiguousarray(np.moveaxis(s, -1, ax)), np.ascontiguousarray(
        np.moveaxis(idx, -1, ax)
    )


def unique(x, *, stable: bool = False, occurrence: str = "first", matlab_layout: bool = False):
    """MATLAB's ``[u, ia, ic] = unique(x)`` of all elements in MATLAB's order (``x(:)``): sorted
    (``stable``: in order of first occurrence), each NaN a value of its own (last), -0 and +0
    one value (the first in sorted order); ``ia`` the 0-based MATLAB linear indices of the 'first'
    (or 'last') occurrences, ``ic`` the index in ``u`` of each element. 1-D results."""
    if occurrence not in ("first", "last"):
        raise ValueError("occurrence must be 'first' or 'last'")
    a = _array(x)
    v = _fortran_flat(_mview(a, matlab_layout), a.dtype) if a.ndim >= 2 else a.reshape(-1)
    n = v.size
    if n == 0:
        e = np.zeros(0, np.intp)
        return v.copy(), e, e.copy()
    s, order = sort(v)
    new = np.ones(n, bool)
    new[1:] = s[1:] != s[:-1]
    if s.dtype.kind == "f":
        new[1:] |= np.isnan(s[1:])
    grp = np.cumsum(new) - 1
    starts = np.flatnonzero(new)
    ends = np.append(starts[1:], n) - 1
    ia = order[starts] if occurrence == "first" else order[ends]
    ic = np.empty(n, np.intp)
    ic[order] = grp
    u = s[starts]
    if stable:
        perm = np.argsort(order[starts], kind="stable")
        rank = np.empty_like(perm)
        rank[perm] = np.arange(perm.size)
        return u[perm], ia[perm], rank[ic]
    return u, ia, ic
