"""The element-wise primitives and the fixed-order sum (petakit5d.kernels.ops, K08).

Oracle cases compare every primitive with MATLAB R2024b (every class pair of the conversions,
the decon expressions with scale, offset, clamp and mask, comparisons with double scalars, and
sum(x(:), 'double') at one thread); in each case the Python side also checks that the
compiled kernel equals the NumPy reference bit for bit. The other tests check thread-count
independence, layouts, ``out``, errors and the fall-back without the extension.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.kernels import ops
from petakit5d.kernels.ops import SUM_BLOCK
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Vol
from petakit5d.testing.synth import CounterRNG

SOURCES = (
    "logical",
    "int8",
    "uint8",
    "int16",
    "uint16",
    "int32",
    "uint32",
    "int64",
    "uint64",
    "single",
    "double",
)
TARGETS = SOURCES[1:]
INT_LIMITS = {
    "int8": (-128, 127),
    "uint8": (0, 255),
    "int16": (-32768, 32767),
    "uint16": (0, 65535),
    "int32": (-(2**31), 2**31 - 1),
    "uint32": (0, 2**32 - 1),
    "int64": (-(2**63), 2**63 - 1),
    "uint64": (0, 2**64 - 1),
}
NP = {
    "logical": np.bool_,
    "single": np.float32,
    "double": np.float64,
    **{c: np.dtype(c).type for c in INT_LIMITS},
}


# -- inputs ------------------------------------------------------------------------------------


def _bits(words, dtype) -> np.ndarray:
    return np.array(words, dtype=np.uint64 if dtype == np.float64 else np.uint32).view(dtype)


def _edge_double_list() -> np.ndarray:
    v = [0.0, -0.0, np.inf, -np.inf, 0.49999999999999994, -0.49999999999999994]
    v += [float(np.float32(0.49999997)), -float(np.float32(0.49999997))]
    for k in range(-6, 7):
        h = k + 0.5
        v += [h, np.nextafter(h, np.inf), np.nextafter(h, -np.inf), k, k + 1e-9, k - 1e-9]
    limits = sorted({b for lo_hi in INT_LIMITS.values() for b in lo_hi})
    for b in limits + [128, -129, 256, 32768, -32769, 65536, 2**31, -(2**31) - 1, 2**32]:
        b = float(b)
        for d in (0.0, 0.5, -0.5, 0.49999, -0.49999, 1.0, -1.0):
            v.append(b + d)
        v += [np.nextafter(b, np.inf), np.nextafter(b, -np.inf)]
        v += [np.nextafter(b + 0.5, np.inf), np.nextafter(b - 0.5, -np.inf)]
    v += [2.0**53, 2.0**53 + 2, 2.0**63, -(2.0**63), 2.0**64, 1e19, 1e20, -1e20]
    v += [3.4028234663852886e38, 3.5e38, -3.5e38, 1e-40, 1e-46, 5e-324, 2.2250738585072014e-308]
    v += [1 + 2.0**-24, 1 + 3 * 2.0**-24, 16777217.0, 16777219.0]
    return np.array(v, dtype=np.float64)


def _edge_doubles() -> np.ndarray:
    nans = _bits(
        [0xFFF8000000000000, 0x7FF8000000000000, 0x7FF4000000000001, 0xFFF0000000000123],
        np.float64,
    )
    rng = CounterRNG(8, "ops/doubles")
    scale = 2.0 ** rng.integers(-40, 70, 600).reshape(-1).astype(np.float64)
    rnd = rng.normal(600).reshape(-1) * scale
    return np.concatenate([_edge_double_list(), nans, rnd])


def _edge_singles() -> np.ndarray:
    with np.errstate(over="ignore"):
        a = _edge_double_list().astype(np.float32)
    half = np.float32(0.5)
    extra = [np.nextafter(half, np.float32(0)), -np.nextafter(half, np.float32(0))]
    for k in range(-6, 7):
        h = np.float32(k) + half
        extra += [np.nextafter(h, np.float32(np.inf)), np.nextafter(h, np.float32(-np.inf))]
    nans = _bits([0xFFC00000, 0x7FC00000, 0x7FA00001, 0xFF800123], np.float32)
    rng = CounterRNG(8, "ops/singles")
    scale = 2.0 ** rng.integers(-30, 40, 600).reshape(-1).astype(np.float64)
    rnd = (rng.normal(600).reshape(-1) * scale).astype(np.float32)
    return np.concatenate([a, np.array(extra, np.float32), nans, rnd])


def _edge_ints(cls: str) -> np.ndarray:
    lo, hi = INT_LIMITS[cls]
    v = [0, 1, lo, hi, lo + 1, hi - 1]
    for c in INT_LIMITS.values():
        for b in c + (c[0] - 1, c[1] + 1):
            if lo <= b <= hi:
                v.append(b)
    big = [2**24 + 1, 2**53 + 1, 2**53 + 2**29 + 1, 2**62 + 2**38 + 1, 2**63 - 2**39 - 1]
    big += [2**63 + 2**39 + 1, 2**64 - 2**39 - 1, 16777217, 2**31 - 1, 2**32 - 129]
    for b in big + [-b for b in big]:
        if lo <= b <= hi:
            v.append(b)
    rng = CounterRNG(8, f"ops/{cls}")
    if np.dtype(cls).itemsize == 8:
        r = np.ascontiguousarray(rng.uint32(1200).reshape(-1)).view(np.dtype(cls))
    else:
        r = rng.integers(lo, hi, 600, dtype=np.int64).reshape(-1)
    return np.array(v + [int(x) for x in r], dtype=np.dtype(cls))


def _values(cls: str) -> np.ndarray:
    """edge values and random values of a class, as a Fortran-ordered (n, 3, 2) array"""
    if cls == "double":
        v = _edge_doubles()
    elif cls == "single":
        v = _edge_singles()
    elif cls == "logical":
        v = CounterRNG(8, "ops/logical").uint32(301).reshape(-1) % 3 == 0
    else:
        v = _edge_ints(cls)
    n = -(-v.size // 6) * 6
    p = np.zeros(n, dtype=v.dtype)
    p[: v.size] = v
    return np.asfortranarray(p.reshape(n // 6, 3, 2, order="F"))


def _same(a, b) -> bool:
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    return a.tobytes(order="A") == b.tobytes(order="A") if a.size else True


def _bitwise(a, b) -> bool:
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    return np.array_equal(
        np.ascontiguousarray(a).view(np.uint8), np.ascontiguousarray(b).view(np.uint8)
    )


def _both(kernel, reference):
    """the kernel's result, after checking that it equals the reference bit for bit"""

    def run(*args, **kwargs):
        got, ref = kernel(*args, **kwargs), reference(*args, **kwargs)
        g = got if isinstance(got, tuple) else (got,)
        r = ref if isinstance(ref, tuple) else (ref,)
        for i, (x, y) in enumerate(zip(g, r, strict=True)):
            if not _bitwise(x, y):
                raise AssertionError(f"output {i + 1}: kernel and NumPy reference differ")
        return got

    return run


# -- conversions -------------------------------------------------------------------------------


def _convert_all(mode):
    def k(x):
        return tuple(ops.convert(x, c, mode=mode) for c in TARGETS)

    def r(x):
        return tuple(ops.convert_ref(x, c, mode=mode) for c in TARGETS)

    return _both(k, r)


_CLASS_CODE = "x = in{1}; out = {" + ", ".join(f"{c}(x)" for c in TARGETS) + "};"
_CAST_CODE = "x = in{1}; out = {" + ", ".join(f"cast(x, '{c}')" for c in TARGETS) + "};"

CONVERT_CASES = [
    Case(
        f"convert_{mode}_{src}",
        code=_CLASS_CODE if mode == "class" else _CAST_CODE,
        args=[Vol(In("x"))],
        nargout=len(TARGETS),
        outputs=list(TARGETS),
        out_roles=("vol",) * len(TARGETS),
        python=_convert_all(mode),
        inputs={"x": Value(lambda src=src: _values(src))},
    )
    for mode in ("class", "cast")
    for src in SOURCES
]


# -- scale, offset, clamp, mask ----------------------------------------------------------------


def _decon_volume(dtype) -> np.ndarray:
    """a (y, x, z) volume like decon data (values from -200 to 2800, every 7th 0), with the
    edge values of doubles at every third voxel from the second on"""
    rng = CounterRNG(9, f"ops/decon/{np.dtype(dtype).name}")
    # the values in MATLAB's linear order (1-D draws: the same bytes as before the convention)
    v = rng.uniform(37 * 29 * 11) * 3000.0 - 200.0
    v[::7] = 0.0
    e = _edge_doubles()
    k = min(e.size, v[1::3].size)
    v[1 : 1 + 3 * k : 3] = e[:k]
    with np.errstate(over="ignore", invalid="ignore"):
        return np.asfortranarray(v.astype(dtype).reshape(37, 29, 11, order="F"))


def _int_volume(cls) -> np.ndarray:
    v = _values(cls).reshape(-1, order="F")
    rng = CounterRNG(9, f"ops/intvol/{cls}")
    lo, hi = INT_LIMITS[cls]
    r = rng.integers(max(lo, -3000), min(hi, 70000), 37 * 29 * 11 - v.size, dtype=np.int64)
    a = np.concatenate([v, r.reshape(-1).astype(np.dtype(cls))])
    return np.asfortranarray(a.reshape(37, 29, 11, order="F"))


def _so(**kw):
    def k(x, *args):
        return ops.scale_offset(x, *args, **kw)

    def r(x, *args):
        return ops.scale_offset_ref(x, *args, **kw)

    return _both(k, r)


def _so_mask(**kw):
    def k(x, s, o, m):
        y = ops.scale_offset(x, s, o, mask=m)
        return y, ops.scale_offset(x, s, o, mask=m, out_class="uint16", **kw)

    def r(x, s, o, m):
        y = ops.scale_offset_ref(x, s, o, mask=m)
        return y, ops.scale_offset_ref(x, s, o, mask=m, out_class="uint16", **kw)

    return _both(k, r)


def _mask_of(x) -> np.ndarray:
    return np.asfortranarray(ops.compare_ref(x, "ne", 0.0))


AFFINE_CASES = (
    [
        # decon_lucy_function: J_2 = (J_2 * scaleFactor + deconOffset) .* I_mask; uint16(J_2)
        Case(
            f"decon_scale_offset_mask_{cls}",
            code="y = (in{1} * in{2} + in{3}) .* in{4}; out = {y, uint16(y)};",
            args=[Vol(In("x")), 1.7, -3.25, Vol(In("m"))],
            nargout=2,
            out_roles=("vol", "vol"),
            python=_so_mask(),
            inputs={
                "x": Value(lambda cls=cls: _decon_volume(NP[cls])),
                "m": Value(lambda cls=cls: _mask_of(_decon_volume(NP[cls]))),
            },
        )
        for cls in ("single", "double")
    ]
    + [
        # decon_lucy_function: J_2 = J_2 * scaleFactor; uint16(J_2) (and cast)
        Case(
            f"decon_scale_uint16_{cls}_{mode}",
            code="y = in{1} * in{2}; out = {"
            + ("uint16(y)" if mode == "class" else "cast(y, 'uint16')")
            + "};",
            args=[Vol(In("x")), 0.37],
            out_roles=("vol",),
            python=_so(out_class="uint16", mode=mode),
            inputs={"x": Value(lambda cls=cls: _decon_volume(NP[cls]))},
        )
        for cls in ("single", "double")
        for mode in ("class", "cast")
    ]
    + [
        # decon_lucy_function / RLdecon: max(I - Background, 0) (offset -Background, lower bound 0)
        Case(
            f"background_max0_{cls}",
            code="out = {max(in{1} - in{2}, 0)};",
            args=[Vol(In("x")), 100.0],
            out_roles=("vol",),
            python=_so(lo=0.0),
            py_args=[Vol(In("x")), None, -100.0],
            inputs={"x": Value(lambda cls=cls: _decon_volume(NP[cls]))},
        )
        for cls in ("single", "double")
    ]
    + [
        # integer data with double scalars: rounded half away and saturated at every step
        Case(
            f"int_scale_offset_{cls}",
            code="x = in{1}; out = {x * in{2}, x * in{2} + in{3}, x + in{3}, x * -1.5, x * NaN,"
            " x * 3.3e9, min(x * in{2}, 1000), max(x - 100.5, 7)};",
            args=[Vol(In("x")), 1.5, -100.5],
            nargout=8,
            out_roles=("vol",) * 8,
            python=_both(
                lambda x, s, o, f=ops.scale_offset: (
                    f(x, s),
                    f(x, s, o),
                    f(x, None, o),
                    f(x, -1.5),
                    f(x, dtypes.NAN),
                    f(x, 3.3e9),
                    f(x, s, hi=1000.0),
                    f(x, None, -100.5, lo=7.0),
                ),
                lambda x, s, o, f=ops.scale_offset_ref: (
                    f(x, s),
                    f(x, s, o),
                    f(x, None, o),
                    f(x, -1.5),
                    f(x, dtypes.NAN),
                    f(x, 3.3e9),
                    f(x, s, hi=1000.0),
                    f(x, None, -100.5, lo=7.0),
                ),
            ),
            inputs={"x": Value(lambda cls=cls: _int_volume(cls))},
        )
        for cls in ("int8", "uint8", "int16", "uint16", "int32", "uint32")
    ]
    + [
        # single and double with scalars that round, overflow and make NaN
        Case(
            f"float_scale_offset_{cls}",
            code="x = in{1}; out = {x * 0.1, x * 0.1 + 1e-3, x + 1e30, x * 1e300, x * NaN, x * -0,"
            " (x * 3 - 2) .* (x ~= 0), int16(x * 0.01 + 0.5), uint8(x * 0.01)};",
            args=[Vol(In("x"))],
            nargout=9,
            out_roles=("vol",) * 9,
            python=_both(
                lambda x, f=ops.scale_offset, m=ops.compare: (
                    f(x, 0.1),
                    f(x, 0.1, 1e-3),
                    f(x, None, 1e30),
                    f(x, 1e300),
                    f(x, dtypes.NAN),
                    f(x, -0.0),
                    f(x, 3.0, -2.0, mask=m(x, "ne", 0.0)),
                    f(x, 0.01, 0.5, out_class="int16"),
                    f(x, 0.01, out_class="uint8"),
                ),
                lambda x, f=ops.scale_offset_ref, m=ops.compare_ref: (
                    f(x, 0.1),
                    f(x, 0.1, 1e-3),
                    f(x, None, 1e30),
                    f(x, 1e300),
                    f(x, dtypes.NAN),
                    f(x, -0.0),
                    f(x, 3.0, -2.0, mask=m(x, "ne", 0.0)),
                    f(x, 0.01, 0.5, out_class="int16"),
                    f(x, 0.01, out_class="uint8"),
                ),
            ),
            inputs={"x": Value(lambda cls=cls: _values(cls))},
        )
        for cls in ("single", "double")
    ]
)


def _clamp_tuple(f):
    return lambda x: (
        f(x, 0.0, 1.0),
        f(x, None, 1.0),
        f(x, 0.0),
        f(x, -0.0),
        f(x, None, -0.0),
        f(x, dtypes.NAN, 2.0),
        f(x, -3.0, dtypes.NAN),
        f(x, 1.0, 0.0),
        f(x, -1e40, 1e40),
        f(x, 0.1, 0.7),
    )


CLAMP_CASES = [
    Case(
        f"clamp_{cls}",
        code="x = in{1}; out = {max(min(x, 1), 0), min(x, 1), max(x, 0), max(x, -0), min(x, -0),"
        " max(min(x, 2), NaN), max(min(x, NaN), -3), max(min(x, 0), 1), max(min(x, 1e40), -1e40),"
        " max(min(x, 0.7), 0.1)};",
        args=[Vol(In("x"))],
        nargout=10,
        out_roles=("vol",) * 10,
        python=_both(_clamp_tuple(ops.clamp), _clamp_tuple(ops.clamp_ref)),
        inputs={"x": Value(lambda cls=cls: _values(cls))},
    )
    for cls in ("single", "double", "uint16", "int16", "uint8", "int32", "uint32", "int8")
]


# -- comparisons -------------------------------------------------------------------------------

_CMP_OPS = ("eq", "ne", "lt", "le", "gt", "ge")
_CMP_CODE = "x = in{1}; v = in{2}; out = {x == v, x ~= v, x < v, x <= v, x > v, x >= v};"


def _cmp_all(f):
    return lambda x, v: tuple(f(x, op, v) for op in _CMP_OPS)


_CMP_VALUES = {
    "float": (0.0, -0.0, 0.1, 0.5, float("nan"), float("inf"), 16777217.0, 3.4028235677973366e38),
    "int": (0.0, 0.5, -0.5, 255.0, 65535.4, -1e10, 1e10, float("nan"), float("-inf")),
    "int64": (
        0.0,
        -0.5,
        2.0**53,
        2.0**53 + 2,
        2.0**63,
        -(2.0**63),
        2.0**64,
        9.2e18,
        -9.2e18,
        1.5e19,
        float("inf"),
    ),
    "logical": (0.0, 1.0, 0.5, -0.0, float("nan")),
}


def _cmp_kind(cls: str) -> str:
    if cls in ("single", "double"):
        return "float"
    if cls in ("int64", "uint64"):
        return "int64"
    return "logical" if cls == "logical" else "int"


def _vid(v: float) -> str:
    return (
        repr(v)
        .replace("-", "m")
        .replace(".", "p")
        .replace("+", "")
        .replace("inf", "Inf")
        .replace("nan", "NaN")
    )


COMPARE_CASES = [
    Case(
        f"compare_{cls}_{_vid(v)}",
        code=_CMP_CODE,
        args=[Vol(In("x")), v],
        nargout=6,
        outputs=list(_CMP_OPS),
        out_roles=("vol",) * 6,
        python=_both(_cmp_all(ops.compare), _cmp_all(ops.compare_ref)),
        inputs={"x": Value(lambda cls=cls: _values(cls))},
    )
    for cls in SOURCES
    for v in _CMP_VALUES[_cmp_kind(cls)]
]


# -- sums --------------------------------------------------------------------------------------

_ONE_THREAD = (
    "n0 = maxNumCompThreads(1); s = sum(in{1}(:), 'double'); maxNumCompThreads(n0); out = {s};"
)

# decon_sum_double's MATLAB fall-back (block_sum, parked decon change), for 89000+ elements
_BLOCK_SUM = """x = in{1}; B = 4096; nPiece = 256; x = x(:); n = numel(x); nb = ceil(n / B);
bs = zeros(nb, 1);
for b = 1 : nPiece : nb
    e = min(b + nPiece - 1, nb);
    v = x((b - 1) * B + 1 : min(e * B, n));
    v(end + 1 : (e - b + 1) * B) = -0;
    L = sum(reshape(v, 4, B / 4, []), 2, 'double');
    bs(b : e) = ((L(1, :) + L(2, :)) + L(3, :)) + L(4, :);
end
s = bs(1);
for b = 2 : nb
    s = s + bs(b);
end
out = {s};"""


def _sum_data(n: int, cls: str, kind: str = "wide") -> np.ndarray:
    """n values whose sum depends on the order (magnitudes from 2^-30 to 2^30)"""
    rng = CounterRNG(10, f"ops/sum/{cls}/{kind}/{n}")
    if kind == "negzero":
        v = np.full(n, -0.0)
    elif cls in ("single", "double"):
        v = rng.normal(n).reshape(-1) * 2.0 ** rng.integers(-30, 30, n).reshape(-1)
        if kind == "special" and n >= 8:
            v[n // 3] = np.inf
            v[n // 2] = dtypes.NAN if n % 2 else -np.inf
        if kind in NAN_PROBES:
            # distinct payloads where each order of operands gives another result
            v = v.astype(NP[cls])
            for j, a in enumerate(NAN_PROBES[kind][1]):
                v[a] = _NAN_PAYLOADS[cls + "8"][j]
            return np.asfortranarray(v)
        if kind == "nans":
            v = v.astype(NP[cls])
            pays = _NAN_PAYLOADS[cls]
            for j, at in enumerate(sorted({n // 2, n // 3 + 1, n - 1, n // 4 + 2})):
                v[at] = pays[j % len(pays)]
            return np.asfortranarray(v)
    elif cls == "logical":
        return np.asfortranarray(rng.uint32(n).reshape(-1) % 2 == 0)
    else:
        lo, hi = INT_LIMITS[cls]
        v = rng.integers(lo, hi, n, dtype=np.int64).reshape(-1)
        return np.asfortranarray(v.astype(np.dtype(cls)))
    return np.asfortranarray(v.astype(NP[cls]))


_NAN_PAYLOADS = {
    "double": _bits([0x7FF8000000000001, 0xFFF8000000000002, 0x7FF0000000000003], np.float64),
    "single": _bits([0x7FC00001, 0xFFC00002, 0x7F800003], np.float32),
    "double8": _bits([0x7FF8000000000000 + k for k in range(1, 9)], np.float64),
    "single8": _bits([0x7FC00000 + k for k in range(1, 9)], np.float32),
}


_B = SUM_BLOCK
# name: (n, positions of the payloads 1, 2, ...): in lanes 1 and 2 of the four-lane sum and of a
# block, in a lane's tail, in lanes 0 and 3, one per block
NAN_PROBES = {
    "nanseq": (12, (5, 9, 2, 6)),
    "nanlanes": (60, (5, 9, 2, 6)),
    "nanlanes03": (60, (4, 7)),
    "nantail": (50, (1, 49)),
    "nanblocklanes": (100003, tuple(3 * _B + a for a in (5, 9, 2, 6))),
    "nanblocktail": (100003, (24 * _B + 1, 24 * _B + 1697)),
    "nanblocks": (100003, (2 * _B + 7, 5 * _B)),
    "nanblocks3": (100003, (_B + 3, 3 * _B + 2, 7 * _B + 1)),
    "nantailp60": (63, (56, 60)),
    "nantailp56": (63, (52, 56)),
    "nantailp48": (63, (44, 48)),
    "nantailp32": (63, (28, 32)),
    "nanbtail1656": (24 * _B + 1660, (24 * _B + 1652, 24 * _B + 1656)),
    "nanbtail1648": (24 * _B + 1660, (24 * _B + 1644, 24 * _B + 1648)),
    "nanbtail1632": (24 * _B + 1660, (24 * _B + 1628, 24 * _B + 1632)),
}


def _sum_vol(shape, cls) -> np.ndarray:
    n = int(np.prod(shape))
    return np.asfortranarray(_sum_data(n, cls).reshape(shape, order="F"))


def _sum(x):
    return ops.sum_double(x)


# MATLAB sizes of volumes and images whose 4096-blocks and lanes cross columns and planes:
# ny above 4096, planes of exactly 4096, ny of 3, a column, odd sizes in each regime
SUM_VOLUMES = (
    (5000, 3, 7),
    (37, 29, 101),
    (64, 64, 25),
    (3, 40000),
    (89001, 1),
    (211, 223),
    (7, 5),
    (13, 11, 3),
)
# NaN payloads at MATLAB linear indices of a volume: lanes, blocks and the last block's tail
NAN_VOLUMES = {
    (37, 29, 101): (5, 9, 2 * _B + 7, 5 * _B, 26 * _B + 3, 108371, 108372),
    (211, 223): (1, 5, 47049, 47050, 47052),
}


def _nan_vol(shape, cls) -> np.ndarray:
    n = int(np.prod(shape))
    v = _sum_data(n, cls).reshape(-1)
    for j, a in enumerate(NAN_VOLUMES[shape]):
        v[a] = _NAN_PAYLOADS[cls + "8"][j]
    return np.asfortranarray(v.reshape(shape, order="F"))


def _regime_data(n: int) -> np.ndarray:
    v = np.zeros(n)
    v[0], v[1], v[5] = 2.0**53, 1.0, 1.0
    return v


def _sid(shape) -> str:
    return "x".join(str(v) for v in shape)


SUM_SIZES = (0, 1, 2, 3, 5, 43, 44, 45, 47, 4095, 4096, 4097, 88999, 89000, 89001, 102403)
SUM_CASES = (
    [
        Case(
            f"sum_{cls}_{n}",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda n=n, cls=cls: _sum_data(n, cls))},
        )
        for cls in ("double", "single")
        for n in SUM_SIZES
    ]
    + [
        Case(
            f"sum_{cls}_{n}",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda n=n, cls=cls: _sum_data(n, cls))},
        )
        for cls in ("uint16", "int32", "uint32", "int8", "logical")
        for n in (5, 1001, 89001)
    ]
    + [
        Case(
            f"sum_{kind}_{cls}_{n}",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda n=n, cls=cls, kind=kind: _sum_data(n, cls, kind))},
        )
        for kind in ("negzero", "special", "nans")
        for cls in ("double", "single")
        for n in (3, 50, 100003)
    ]
    + [
        Case(
            f"sum_{kind}_{cls}_{n}",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda n=n, cls=cls, kind=kind: _sum_data(n, cls, kind))},
        )
        for kind, (n, _) in NAN_PROBES.items()
        for cls in ("double", "single")
    ]
    + [
        # values where index order and four lanes differ: 2^53 + 1 + 1 is 2^53 in index order,
        # 2^53 + 2 with the ones in one lane (MATLAB switches at 44 elements)
        Case(
            f"sum_regime_double_{n}",
            code=_ONE_THREAD,
            args=[In("x")],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda n=n: _regime_data(n))},
        )
        for n in (43, 44)
    ]
    + [
        Case(
            "sum_volume_single_1000003",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda: _sum_vol((101, 97, 103), "single"))},
        ),
        Case(
            "sum_volume_double_odd",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda: _sum_vol((63, 41, 37), "double"))},
        ),
    ]
    + [
        Case(
            f"sum_mlorder_{cls}_{_sid(shape)}",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda shape=shape, cls=cls: _sum_vol(shape, cls))},
        )
        for cls in ("double", "single")
        for shape in SUM_VOLUMES
    ]
    + [
        Case(
            f"sum_mlorder_nans_{cls}_{_sid(shape)}",
            code=_ONE_THREAD,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda shape=shape, cls=cls: _nan_vol(shape, cls))},
        )
        for cls in ("double", "single")
        for shape in NAN_VOLUMES
    ]
    + [
        # below 89000 elements MATLAB's sum is the same at any thread count (this session's: 4+)
        Case(
            f"sum_threads_{cls}_{n}",
            code="out = {sum(in{1}(:), 'double')};",
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda n=n, cls=cls: _sum_data(n, cls))},
        )
        for cls in ("double", "single")
        for n in (44, 50001, 88999)
    ]
    + [
        # decon_sum_double's MATLAB fall-back gives the same bits
        Case(
            f"sum_block_sum_{cls}_{n}",
            code=_BLOCK_SUM,
            args=[Vol(In("x"))],
            python=_both(_sum, ops.sum_double_ref),
            inputs={"x": Value(lambda n=n, cls=cls: _sum_data(n, cls))},
        )
        for cls in ("double", "single")
        for n in (89000, 1048577)
    ]
)

CASES = CONVERT_CASES + AFFINE_CASES + CLAMP_CASES + COMPARE_CASES + SUM_CASES


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- kernel against reference, threads, layouts ------------------------------------------------

needs_ext = pytest.mark.skipif(not ops.extension_available(), reason="extension not built")


def _vol(make, *args) -> np.ndarray:
    """an input generator's MATLAB array as the port's (z, y, x) C array"""
    return from_matlab(make(*args))


@needs_ext
@pytest.mark.parametrize("src", SOURCES)
@pytest.mark.parametrize("mode", ["class", "cast"])
def test_convert_equals_reference_for_every_pair(src, mode):
    x = _vol(_values, src)
    for c in TARGETS:
        assert _bitwise(ops.convert(x, c, mode=mode), ops.convert_ref(x, c, mode=mode)), (src, c)


def _ext():
    from petakit5d import _ext

    return _ext.load("core")


@needs_ext
@pytest.mark.parametrize("nt", [1, 2, 3, 5, 8])
def test_results_do_not_depend_on_threads(nt):
    m = _ext()
    rng = CounterRNG(11, "ops/threads")
    x = (rng.normal((17, 129, 131)).reshape(17, 129, 131) * 1000).astype(np.float32)
    x = np.ascontiguousarray(x)
    n = x.size
    flat = x.reshape(-1)
    ref_u16 = ops.convert(x, "uint16", num_threads=1).reshape(-1)
    out = np.empty(n, np.uint16)
    m.ops_convert(flat, out, False, nt)
    assert _bitwise(out, ref_u16)
    ref_a = ops.scale_offset(x, 1.7, -3.0, lo=0.0, out_class="uint16", num_threads=1)
    out = np.empty(n, np.uint16)
    m.ops_affine(flat, out, True, 1.7, True, -3.0, False, 0.0, True, 0.0, None, False, nt)
    assert _bitwise(out, ref_a.reshape(-1))
    ref_c = ops.compare(x, "gt", 0.5, num_threads=1).reshape(-1)
    out = np.empty(n, np.bool_)
    m.ops_compare(flat, out, 4, 0.5, nt)
    assert _bitwise(out, ref_c)
    for size in (89000, 4096 * 37 + 5, n):
        s1 = m.ops_sum_double(flat[:size], 1)
        assert _bitwise(np.float64(m.ops_sum_double(flat[:size], nt)), np.float64(s1))
        assert _bitwise(np.float64(s1), np.float64(ops.sum_double_ref(flat[:size])))
    # MATLAB's order of the (z, y, x) volume, read strided
    for vol in (x, x[:6], x[:, :40, :]):
        v = np.ascontiguousarray(vol)
        s1 = m.ops_sum_double_yx(v, 1)
        assert _bitwise(np.float64(m.ops_sum_double_yx(v, nt)), np.float64(s1))
        assert _bitwise(np.float64(s1), np.float64(ops.sum_double_ref(v)))


@needs_ext
@pytest.mark.parametrize(
    "shape",
    [(3, 1000, 37), (2, 4100, 11), (5, 3, 7000), (1, 89001, 2), (7, 64, 64), (30, 31, 97)],
)
def test_strided_sum_tiles_and_thread_ranges(shape):
    # tiles of columns, blocks across columns and planes, thread ranges that split tiles
    m = _ext()
    n = int(np.prod(shape))
    v = from_matlab(_sum_data(n, "single").reshape(shape[1], shape[2], shape[0], order="F"))
    want = np.float64(ops.sum_double_ref(v))
    for nt in (1, 2, 3, 7):
        assert _bitwise(np.float64(m.ops_sum_double_yx(v, nt)), want), nt


def test_sum_order_matters_and_is_matlabs():
    # the data make the order visible: MATLAB's order, memory order, pairwise and sequential
    # sums differ
    v = _vol(_sum_vol, (61, 47, 37), "double")
    s = ops.sum_double(v)
    assert s != ops.sum_double(v.reshape(-1))
    assert s != float(np.cumsum(v.reshape(-1))[-1]) or s != float(np.sum(v))
    assert _bitwise(np.float64(s), np.float64(ops.sum_double_ref(v)))


def test_sum_signed_zero_and_empty():
    assert np.signbit(ops.sum_double(np.full(5, -0.0)))
    assert np.signbit(ops.sum_double(np.full(100003, -0.0, np.float32)))
    assert np.signbit(ops.sum_double(np.full((3, 211, 223), -0.0, np.float32)))
    assert not np.signbit(ops.sum_double(np.zeros(0)))
    assert ops.sum_double(np.zeros((0, 3))) == 0.0
    assert ops.sum_double(np.zeros((4, 0, 3))) == 0.0


@pytest.mark.parametrize("shape", [(61, 47, 37), (5000, 3, 7), (3, 40000), (211, 223), (7, 5)])
def test_sum_of_matlabs_array_equals_the_volumes(shape):
    m = _sum_vol(shape, "single")  # MATLAB's array, Fortran order
    v = from_matlab(m)
    s = ops.sum_double(v)
    assert _bitwise(np.float64(ops.sum_double(m, matlab_layout=True)), np.float64(s))
    assert _bitwise(np.float64(ops.sum_double(to_matlab(v), matlab_layout=True)), np.float64(s))
    c = np.ascontiguousarray(m)  # MATLAB's array held in C order
    assert _bitwise(np.float64(ops.sum_double(c, matlab_layout=True)), np.float64(s))
    assert _bitwise(np.float64(ops.sum_double_ref(m, matlab_layout=True)), np.float64(s))
    # a non-contiguous view is read in MATLAB's order too
    w = np.zeros((v.shape[0], v.shape[1] + 2, v.shape[-1]) if v.ndim == 3 else v.shape, v.dtype)
    if v.ndim == 3:
        w[:, 1:-1, :] = v
        assert _bitwise(np.float64(ops.sum_double(w[:, 1:-1, :])), np.float64(s))


def test_sum_nan_payloads_follow_matlabs_order_in_volumes():
    v = _vol(_nan_vol, (37, 29, 101), "double")
    got = np.float64(ops.sum_double(v)).view(np.uint64)
    assert got == np.float64(ops.sum_double_ref(v)).view(np.uint64)
    m = to_matlab(v)
    assert got == np.float64(ops.sum_double(m, matlab_layout=True)).view(np.uint64)


def test_layouts_c_fortran_and_strided():
    x = _vol(_decon_volume, np.float32)
    f = to_matlab(x)  # a MATLAB-layout Fortran array, as a kernel may hold one
    for arr in (x, f, x[:, ::2, :], x[0], x[:, :, 3]):
        y = ops.convert(arr, "uint16")
        assert _bitwise(np.ascontiguousarray(y), ops.convert_ref(arr, "uint16"))
        assert y.flags.f_contiguous or y.flags.c_contiguous
    assert ops.convert(x, "uint16").flags.c_contiguous
    assert ops.convert(x[:, ::2, :], "uint16").flags.c_contiguous
    assert ops.convert(f, "uint16").flags.f_contiguous
    m = x != 0
    for mk in (m, np.asfortranarray(m), np.repeat(m, 2, axis=-1)[..., ::2]):
        y = ops.scale_offset(x, 2.0, mask=mk)
        assert _bitwise(y, ops.scale_offset_ref(x, 2.0, mask=m))


def test_zero_dim_and_empty():
    for a in (np.float32(0.49999997), np.zeros((0, 4), np.float32), np.array(2.5)):
        for mode in ("class", "cast"):
            assert _bitwise(ops.convert(a, "uint16", mode=mode), ops.convert_ref(a, "uint16"))
    assert ops.scale_offset(np.zeros((3, 0)), 2.0).shape == (3, 0)
    assert ops.compare(np.zeros(0, np.int64), "eq", 1.0).shape == (0,)


def test_class_rounding_of_a_single_scalar_is_exact_unless_scalar_fast():
    h = np.full(1, np.nextafter(np.float32(0.5), np.float32(0)))
    assert ops.convert(h, "uint16")[0] == 0
    assert ops.convert(h, "uint16", scalar_fast=True)[0] == 1
    assert ops.convert(np.repeat(h, 2), "uint16")[0] == 1
    assert ops.convert(np.repeat(h, 2), "uint16", mode="cast")[0] == 0
    assert ops.convert(-np.repeat(h, 2), "int16")[0] == -1


def test_out_and_in_place():
    x = _vol(_decon_volume, np.float32)
    out = np.empty(x.shape, np.uint16)
    r = ops.convert(x, "uint16", out=out)
    assert r is out and _bitwise(out, ops.convert_ref(x, "uint16"))
    y = x.copy()
    ops.scale_offset(y, 2.0, 1.0, out=y)
    assert _bitwise(y, ops.scale_offset_ref(x, 2.0, 1.0))
    buf = np.zeros(x.size + 1, np.float32)
    a = buf[:-1].reshape(x.shape)
    with pytest.raises(ValueError, match="overlaps"):
        ops.convert(a, "single", out=buf[1:].reshape(x.shape))
    with pytest.raises(ValueError, match="must be a uint8"):
        ops.convert(x, "uint8", out=out)
    with pytest.raises(ValueError, match="memory order"):
        ops.convert(x, "uint16", out=np.asfortranarray(np.empty(x.shape, np.uint16)))


def test_errors():
    x = _vol(_values, "uint16")
    with pytest.raises(TypeError, match="MATLAB:mixedClasses"):
        ops.scale_offset(x, 2.0, mask=np.ones(x.shape, bool))
    with pytest.raises(TypeError, match="not supported"):
        ops.scale_offset(_vol(_values, "int64"), 2.0)
    with pytest.raises(TypeError, match="64-bit"):
        ops.sum_double(_vol(_values, "uint64"))
    with pytest.raises(ValueError, match="unsupported target"):
        ops.convert(x, "logical")
    with pytest.raises(ValueError, match="mode"):
        ops.convert(x, "single", mode="round")
    with pytest.raises(ValueError, match="op"):
        ops.compare(x, "lte", 1.0)
    with pytest.raises(ValueError, match="mask"):
        ops.scale_offset(_vol(_values, "single"), 2.0, mask=np.ones(3, bool))
    with pytest.raises(TypeError):
        ops.convert(np.array(["a"]), "single")


@needs_ext
def test_bindings_refuse_conversions():
    m = _ext()
    x = np.zeros(10, np.float32)
    with pytest.raises(TypeError):
        m.ops_convert(x[::2], np.empty(5, np.uint16), True, 1)
    with pytest.raises(ValueError, match="sizes"):
        m.ops_convert(x, np.empty(5, np.uint16), True, 1)
    with pytest.raises(ValueError, match="num_threads"):
        m.ops_sum_double(x, 0)
    with pytest.raises(TypeError):
        m.ops_sum_double_yx(to_matlab(np.zeros((2, 3, 4), np.float32)), 1)
    with pytest.raises(ValueError, match="logical"):
        m.ops_compare(x, np.empty(10, np.uint8), 0, 0.0, 1)


def test_fallback_without_extension(monkeypatch):
    x = _vol(_decon_volume, np.float32)
    want = (
        ops.convert(x, "uint16"),
        ops.scale_offset(x, 0.5, 2.0, lo=0.0, mask=x != 0, out_class="uint8"),
        ops.compare(x, "le", 0.0),
        ops.sum_double(x),
        ops.sum_double(to_matlab(x), matlab_layout=True),
    )
    monkeypatch.setattr(ops, "_core", lambda: None)
    got = (
        ops.convert(x, "uint16"),
        ops.scale_offset(x, 0.5, 2.0, lo=0.0, mask=x != 0, out_class="uint8"),
        ops.compare(x, "le", 0.0),
        ops.sum_double(x),
        ops.sum_double(to_matlab(x), matlab_layout=True),
    )
    for a, b in zip(want, got, strict=True):
        assert _bitwise(np.asarray(a), np.asarray(b))
    out = np.empty(x.shape, np.uint16)
    assert ops.convert(x, "uint16", out=out) is out


def test_matlab_semantics_spot_checks():
    s = np.array([-0.0, 0.0, dtypes.NAN, 2.0], np.float32)
    r = ops.clamp(s, 0.0, 1.0)
    assert np.signbit(r[0]) and not np.signbit(r[1]) and r[2] == 1.0 and r[3] == 1.0
    assert ops.clamp(s, None, 1.0)[2] == 1.0 and ops.clamp(s, 0.0)[2] == 0.0
    u = np.array([1, 2, 3, 65535], np.uint16)
    assert list(ops.scale_offset(u, 1.5)) == [2, 3, 5, 65535]
    assert list(ops.scale_offset(u, None, -1.5)) == [0, 1, 2, 65534]
    assert list(ops.scale_offset(u, dtypes.NAN)) == [0, 0, 0, 0]
    i = np.array([2**53 + 1, -(2**63)], np.int64)
    assert list(ops.compare(i, "gt", 2.0**53)) == [True, False]
    assert list(ops.compare(i, "eq", -(2.0**63))) == [False, True]
    assert _same(ops.sum_double(np.arange(5.0)), 10.0)
