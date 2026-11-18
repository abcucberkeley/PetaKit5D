"""Tests of petakit5d.matlabcompat.reduce: MATLAB's sum, mean, cumsum, max and min with their
index, sort and unique, against MATLAB R2024b (oracle cases) and against a scalar model of the
order of additions (unit tests).

The inputs are integers of up to 24 (single) or 53 (double) bits times powers of two, so that
every order of additions gives its own result: a case passes only with MATLAB's order. MATLAB's
sums of 89,000 elements or more depend on its thread count; those cases run MATLAB at one
thread (``maxNumCompThreads(1)``), whose result the port gives at any thread count. No MATLAB
test covers these built-ins alone (DeconAccelTest, R03's, covers decon_sum_double).
"""

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.matlabcompat import reduce as R
from petakit5d.matlabcompat import stats as S
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, Dim, In, Value, Vol

_C1 = np.uint64(0xBF58476D1CE4E5B9)
_C2 = np.uint64(0x94D049BB133111EB)
_GOLD = np.uint64(0x9E3779B97F4A7C15)
NZ = -0.0


def _bits64(seed: int, n: int) -> np.ndarray:
    """splitmix64 outputs"""
    with np.errstate(over="ignore"):
        z = np.uint64(seed) + np.arange(1, n + 1, dtype=np.uint64) * _GOLD
        z = (z ^ (z >> np.uint64(30))) * _C1
        z = (z ^ (z >> np.uint64(27))) * _C2
        return z ^ (z >> np.uint64(31))


def wide(seed: int, shape, cls: str = "single") -> np.ndarray:
    """values m * 2^e of both signs, m of 24 (single) or 53 (double) bits: exact in the class"""
    n = int(np.prod(shape))
    b = _bits64(seed, n)
    if cls == "single":
        m = (b >> np.uint64(40)).astype(np.float64)
        e = ((b >> np.uint64(8)) % np.uint64(21)).astype(np.int64) - 10
    else:
        m = (b >> np.uint64(11)).astype(np.float64)
        e = ((b >> np.uint64(3)) % np.uint64(41)).astype(np.int64) - 20
    s = np.where(b & np.uint64(1), -1.0, 1.0)
    v = np.ldexp(s * m, e)
    return v.astype(np.float32 if cls == "single" else np.float64).reshape(shape, order="F")


def ints(seed: int, shape, lo: int, hi: int, cls: str = "double") -> np.ndarray:
    n = int(np.prod(shape))
    v = (_bits64(seed, n) % np.uint64(hi - lo + 1)).astype(np.int64) + lo
    return v.astype(np.dtype(cls)).reshape(shape, order="F")


def _one_thread(expr: str, nout: int = 1) -> str:
    outs = ", ".join(f"o{i}" for i in range(nout))
    return (
        f"t__ = maxNumCompThreads(1); [{outs}] = {expr}; maxNumCompThreads(t__); out = {{{outs}}};"
    )


# -- Python sides ------------------------------------------------------------------------------


def py_sum(x, dim=None, *rest):
    kw = {}
    for r in rest:
        if r in ("double", "native", "default"):
            kw["outtype"] = r
        elif r == "omitnan":
            kw["omitnan"] = True
    return R.sum(x, dim, keepdims=dim is not None, **kw)


def py_mean(x, dim=None, *rest):
    kw = {}
    for r in rest:
        if r in ("double", "native", "default"):
            kw["outtype"] = r
        elif r == "omitnan":
            kw["omitnan"] = True
    return R.mean(x, dim, keepdims=dim is not None, **kw)


def py_cumsum(x, dim, *rest):
    return R.cumsum(x, dim, omitnan="omitnan" in rest)


def py_max(x, dim=None, nanflag="omitnan"):
    m, i = R.max(x, dim, omitnan=nanflag == "omitnan", keepdims=dim is not None)
    return m, (np.asarray(i, dtype=np.float64) if dim is not None else i)


def py_min(x, dim=None, nanflag="omitnan"):
    m, i = R.min(x, dim, omitnan=nanflag == "omitnan", keepdims=dim is not None)
    return m, (np.asarray(i, dtype=np.float64) if dim is not None else i)


def py_sort(x, dim, mode="ascend"):
    s, i = R.sort(x, dim, descend=mode == "descend")
    return s, i.astype(np.float64)


def py_unique(x, *opts):
    u, ia, ic = R.unique(
        x, stable="stable" in opts, occurrence="last" if "last" in opts else "first"
    )
    return u, ia.astype(np.float64) + 1, ic.astype(np.float64) + 1


# -- cases -------------------------------------------------------------------------------------


def _vector_cases():
    cs = []
    seed = 100
    for cls in ("single", "double"):
        for n in (1, 2, 7, 43, 44, 45, 47, 100, 1023, 4095, 4096, 4097, 8193, 50000, 88999):
            seed += 1
            cs.append(
                Case(
                    f"sum_row_{cls}_{n}",
                    matlab="sum",
                    args=[In("x")],
                    python=py_sum,
                    inputs={"x": Value(lambda s=seed, n=n, c=cls: wide(s, (1, n), c).ravel())},
                )
            )
            if n in (43, 44, 4097, 88999):
                cs.append(
                    Case(
                        f"mean_row_{cls}_{n}",
                        matlab="mean",
                        args=[In("x")],
                        python=py_mean,
                        inputs={"x": Value(lambda s=seed, n=n, c=cls: wide(s, (1, n), c).ravel())},
                    )
                )
        for n in (89000, 89001, 100003, 409601):
            seed += 1
            make = Value(lambda s=seed, n=n, c=cls: wide(s, (1, n), c).ravel())
            cs.append(
                Case(
                    f"sum_row_{cls}_{n}_one_thread",
                    code=_one_thread("sum(in{1})"),
                    args=[In("x")],
                    python=py_sum,
                    inputs={"x": make},
                )
            )
            cs.append(
                Case(
                    f"mean_row_{cls}_{n}_one_thread",
                    code=_one_thread("mean(in{1})"),
                    args=[In("x")],
                    python=py_mean,
                    inputs={"x": make},
                )
            )
    for n in (44, 4097, 100003):
        seed += 1
        cs.append(
            Case(
                f"sum_row_single_{n}_double_out",
                code=_one_thread("sum(in{1}, 'double')"),
                args=[In("x")],
                python=lambda x: R.sum(x, outtype="double"),
                inputs={"x": Value(lambda s=seed, n=n: wide(s, (1, n), "single").ravel())},
            )
        )
        cs.append(
            Case(
                f"mean_row_single_{n}_double_out",
                code=_one_thread("mean(in{1}, 'double')"),
                args=[In("x")],
                python=lambda x: R.mean(x, outtype="double"),
                inputs={"x": Value(lambda s=seed, n=n: wide(s, (1, n), "single").ravel())},
            )
        )
    # a mean of more than 2^24 single values: the division in double
    cs.append(
        Case(
            "mean_single_over_2pow24",
            code=_one_thread("mean(in{1})"),
            args=[In("x")],
            python=py_mean,
            inputs={
                "x": Value(
                    lambda: (ints(7, (1, 2**24 + 3), 0, 999).ravel() * 1.0).astype(np.float32)
                )
            },
            marks=("slow",),
        )
    )
    return cs


_VOL_SHAPES = (
    ((44, 30, 70), "single"),
    ((130, 80, 9), "double"),
    ((6, 2100, 5), "single"),
    ((7, 300, 50), "double"),
    ((20, 1100, 3), "single"),
    ((3, 1, 9000), "single"),
    ((11, 1500, 2), "double"),  # chunked slices: under 14 rows, but two pages
)
_DIMS = (1, 2, 3, (1, 2), (1, 3), (2, 3), "all")


def _dim_args(d):
    if d == "all":
        return ["all"], [None]
    if isinstance(d, tuple):
        from petakit5d.core.axes import dim_to_axis

        return [list(d)], [tuple(dim_to_axis(k) for k in d)]
    return [Dim(d)], [Dim(d)]


def _dname(d):
    return (
        "all"
        if d == "all"
        else ("d" + "".join(str(k) for k in d) if isinstance(d, tuple) else f"d{d}")
    )


def _volume_cases():
    cs = []
    seed = 300
    for shape, cls in _VOL_SHAPES:
        seed += 1
        tag = "x".join(map(str, shape))
        make = Value(lambda s=seed, sh=shape, c=cls: wide(s, sh, c))
        for d in _DIMS:
            margs, pargs = _dim_args(d)
            for fn, py in (("sum", py_sum), ("mean", py_mean)):
                if fn == "mean" and d not in (2, 3, (1, 3), "all"):
                    continue
                cs.append(
                    Case(
                        f"{fn}_{cls}_{tag}_{_dname(d)}",
                        code=_one_thread(f"{fn}(in{{1}}, in{{2}})"),
                        args=[Vol(In("x")), *margs],
                        py_args=[Vol(In("x")), *pargs],
                        python=py,
                        out_roles=("vol",),
                        inputs={"x": make},
                    )
                )
    # omitnan means over tiles (dimension 4), as the blending's nanmean
    cs.append(
        Case(
            "mean_tiles_omitnan_single",
            code="out = {mean(in{1}, 4, 'omitnan')};",
            args=[Vol(In("t"))],
            python=lambda t: R.mean(t, 0, omitnan=True, keepdims=True),
            out_roles=("vol",),
            inputs={"t": Value(_tiles_with_nan)},
        )
    )
    cs.append(
        Case(
            "sum_tiles_omitnan_single",
            code="out = {sum(in{1}, 4, 'omitnan')};",
            args=[Vol(In("t"))],
            python=lambda t: R.sum(t, 0, omitnan=True, keepdims=True),
            out_roles=("vol",),
            inputs={"t": Value(_tiles_with_nan)},
        )
    )
    cs.append(
        Case(
            "mean_psf_rows_d13_single",
            code="x = in{1}; out = {mean(x([1:10, end-9:end], :, :), [1, 3])};",
            args=[Vol(In("x"))],
            python=lambda x: R.mean(
                np.concatenate([x[:, :10], x[:, -10:]], axis=1), (-3, -2), keepdims=True
            ),
            out_roles=("vol",),
            inputs={"x": Value(lambda: wide(331, (64, 50, 101), "single"))},
        )
    )
    return cs


def _tiles_with_nan():
    t = wide(341, (9, 7, 5, 4), "single")
    b = _bits64(342, t.size).reshape(t.shape, order="F")
    t[(b % np.uint64(3)) == 0] = np.float32(dtypes.NAN)
    t[..., 1] = np.float32(dtypes.NAN)
    return t


def _int_cases():
    cs = []
    u = Value(lambda: ints(401, (1, 3000), 0, 65535, "uint16").ravel())
    for opt in ((), ("native",), ("double",)):
        cs.append(
            Case(
                f"sum_uint16_{'_'.join(opt) or 'default'}",
                matlab="sum",
                args=[In("u"), *opt],
                py_args=[In("u"), None, *opt],
                python=py_sum,
                inputs={"u": u},
            )
        )
    cs.append(
        Case(
            "sum_int8_native_saturates",
            matlab="sum",
            args=[In("i"), "native"],
            py_args=[In("i"), None, "native"],
            python=py_sum,
            inputs={"i": Value(lambda: np.array([100, 100, -50, -127, -100, 90], np.int8))},
        )
    )
    cs.append(
        Case(
            "sum_logical_d1_d2",
            code="out = {sum(in{1}, 1), sum(in{1}, 2)};",
            args=[Vol(In("b"))],
            nargout=2,
            python=lambda b: (R.sum(b, -2, keepdims=True), R.sum(b, -1, keepdims=True)),
            out_roles=("vol", "vol"),
            inputs={"b": Value(lambda: ints(402, (30, 50, 4), 0, 1, "uint8") > 0)},
        )
    )
    for cls, vals in (
        ("uint8", [255, 255, 254]),
        ("int8", [-3, -2]),
        ("uint16", [65535, 65535, 1, 2]),
    ):
        cs.append(
            Case(
                f"mean_{cls}_native",
                matlab="mean",
                args=[In("v"), "native"],
                py_args=[In("v"), None, "native"],
                python=py_mean,
                inputs={"v": Value(lambda c=cls, v=vals: np.array(v, dtype=c))},
            )
        )
        cs.append(
            Case(
                f"mean_{cls}_default",
                matlab="mean",
                args=[In("v")],
                python=py_mean,
                inputs={"v": Value(lambda c=cls, v=vals: np.array(v, dtype=c))},
            )
        )
    return cs


def _cumsum_cases():
    cs = []
    x = Value(lambda: wide(501, (20, 9, 11), "single"))
    for d in (1, 2, 3):
        cs.append(
            Case(
                f"cumsum_single_d{d}",
                matlab="cumsum",
                args=[Vol(In("x")), Dim(d)],
                python=py_cumsum,
                out_roles=("vol",),
                inputs={"x": x},
            )
        )
    cs.append(
        Case(
            "cumsum_uint8_saturates",
            matlab="cumsum",
            args=[In("v"), 2],
            py_args=[In("v"), -1],
            python=py_cumsum,
            inputs={"v": Value(lambda: np.array([200, 100, 50], np.uint8))},
        )
    )
    cs.append(
        Case(
            "cumsum_int8_saturates",
            matlab="cumsum",
            args=[In("v"), 2],
            py_args=[In("v"), -1],
            python=py_cumsum,
            inputs={"v": Value(lambda: np.array([100, 100, -100, -100, -100], np.int8))},
        )
    )
    cs.append(
        Case(
            "cumsum_omitnan",
            matlab="cumsum",
            args=[In("v"), 2, "omitnan"],
            py_args=[In("v"), -1, "omitnan"],
            python=py_cumsum,
            inputs={"v": Value(lambda: np.array([dtypes.NAN, 1.0, dtypes.NAN, 2, NZ]))},
        )
    )
    return cs


def _special_vectors():
    nan = dtypes.NAN
    return {
        "signed_zeros": [0.0, NZ, NZ, 0.0, NZ],
        "zeros_neg_first": [NZ, 0.0, 0.0],
        "nan_first": [nan, 1.0, 2.0, 2.0, 1.0],
        "all_nan": [nan, nan, nan],
        "inf": [np.inf, -np.inf, 1.0, np.inf, -np.inf],
        "ties": [3.0, 1.0, 3.0, 1.0, 2.0, 1.0],
        "nan_ties": [2.0, nan, nan, 2.0, 1.0, nan, 1.0],
        "one": [5.0],
    }


def _index_cases():
    cs = []
    for name, v in _special_vectors().items():
        for cls in ("double", "single"):
            make = Value(lambda v=v, c=cls: np.array(v, dtype=np.float64).astype(c))
            for fn, py in (("max", py_max), ("min", py_min)):
                cs.append(
                    Case(
                        f"{fn}_{name}_{cls}",
                        code=f"[m, i] = {fn}(in{{1}}); out = {{m, i}};",
                        args=[In("v")],
                        nargout=2,
                        python=lambda v, py=py: (lambda m, i: (m, np.float64(i[0] + 1)))(*py(v)),
                        inputs={"v": make},
                    )
                )
                cs.append(
                    Case(
                        f"{fn}_{name}_{cls}_includenan",
                        code=f"[m, i] = {fn}(in{{1}}, [], 'includenan'); out = {{m, i}};",
                        args=[In("v")],
                        nargout=2,
                        python=lambda v, py=py: (lambda m, i: (m, np.float64(i[0] + 1)))(
                            *py(v, None, "includenan")
                        ),
                        inputs={"v": make},
                    )
                )
            cs.append(
                Case(
                    f"sort_{name}_{cls}",
                    code=(
                        "[s, i] = sort(in{1}); [sd, id] = sort(in{1}, 'descend');"
                        " out = {s, i, sd, id};"
                    ),
                    args=[In("v")],
                    nargout=4,
                    python=lambda v: (*py_sort(v, -1), *py_sort(v, -1, "descend")),
                    out_roles=(None, "index1", None, "index1"),
                    inputs={"v": make},
                )
            )
            cs.append(
                Case(
                    f"unique_{name}_{cls}",
                    code=(
                        "[u, ia, ic] = unique(in{1}); [us, ias, ics] = unique(in{1}, 'stable');"
                        " [ul, ial] = unique(in{1}, 'last'); out = {u, ia, ic, us, ias, ics, ial};"
                    ),
                    args=[In("v")],
                    nargout=7,
                    python=lambda v: (
                        *py_unique(v),
                        *py_unique(v, "stable"),
                        py_unique(v, "last")[1],
                    ),
                    inputs={"v": make},
                )
            )
    # volumes: the first maximum in MATLAB's linear order, and along each dimension
    vol = Value(_tie_volume)
    for fn, py in (("max", py_max), ("min", py_min)):
        cs.append(
            Case(
                f"{fn}_volume_all_linear",
                code=f"[m, i] = {fn}(in{{1}}, [], 'all', 'linear'); out = {{m, i}};",
                args=[Vol(In("t"))],
                nargout=2,
                python=lambda t, py=py: py(t),
                out_roles=(None, [("linindex", "t"), float]),
                inputs={"t": vol},
            )
        )
        for d in (1, 2, 3):
            cs.append(
                Case(
                    f"{fn}_volume_d{d}",
                    code=f"[m, i] = {fn}(in{{1}}, [], in{{2}}); out = {{m, i}};",
                    args=[Vol(In("t")), Dim(d)],
                    nargout=2,
                    python=py,
                    out_roles=("vol", ["index1", "vol"]),
                    inputs={"t": vol},
                )
            )
    cs.append(
        Case(
            "sort_volume_d3_descend",
            code="[s, i] = sort(in{1}, 3, 'descend'); out = {s, i};",
            args=[Vol(In("t"))],
            nargout=2,
            python=lambda t: py_sort(t, -3, "descend"),
            out_roles=("vol", ["index1", "vol"]),
            inputs={"t": vol},
        )
    )
    return cs


def _tie_volume():
    v = ints(601, (7, 9, 6), -3, 3, "double")
    b = _bits64(602, v.size).reshape(v.shape, order="F")
    v[(b % np.uint64(11)) == 0] = dtypes.NAN
    v[(v == 0) & ((b % np.uint64(2)) == 0)] = NZ
    return v


def _pair_cases():
    cs = []
    nan = dtypes.NAN
    pairs = {
        "signed_zero": ([0.0, NZ, NZ, 0.0], [NZ, 0.0, NZ, 0.0]),
        "nan": ([nan, 1.0, nan, 2.0], [1.0, nan, nan, 3.0]),
    }
    for name, (a, b) in pairs.items():
        for fn, py in (("max", R.maximum), ("min", R.minimum)):
            cs.append(
                Case(
                    f"{fn}_pair_{name}",
                    matlab=fn,
                    args=[In("a"), In("b")],
                    python=py,
                    inputs={
                        "a": Value(lambda a=a: np.array(a)),
                        "b": Value(lambda b=b: np.array(b)),
                    },
                )
            )
    cs.append(
        Case(
            "max_pair_uint16_double_scalar",
            matlab="max",
            args=[In("a"), 200.7],
            python=R.maximum,
            inputs={"a": Value(lambda: np.array([1, 300, 200, 201], np.uint16))},
        )
    )
    cs.append(
        Case(
            "min_pair_single_double_scalar",
            matlab="min",
            args=[In("a"), 1.25],
            python=R.minimum,
            inputs={"a": Value(lambda: np.array([1.0, 2.0, dtypes.NAN], np.float32))},
        )
    )
    return cs


def _zero_cases():
    cs = []
    z = {
        "sum_negzeros": ("sum", [NZ, NZ, NZ]),
        "sum_negzeros_50": ("sum", [NZ] * 50),
        "sum_empty": ("sum", np.zeros((0, 3))),
        "sum_omitnan_lead": ("sum", [dtypes.NAN, NZ]),
        "sum_omitnan_mid": ("sum", [NZ, dtypes.NAN, NZ]),
        "sum_omitnan_all": ("sum", [dtypes.NAN, dtypes.NAN]),
        "mean_negzeros": ("mean", [NZ, NZ]),
        "mean_empty": ("mean", np.zeros((1, 0))),
        "mean_omitnan_all": ("mean", [dtypes.NAN, dtypes.NAN]),
    }
    cs.append(
        Case(
            "sum_empty_columns",
            matlab="sum",
            args=[In("v")],
            python=lambda v: R.sum(v, 0, keepdims=True, matlab_layout=True),
            inputs={"v": Value(lambda: np.zeros((0, 3)))},
        )
    )
    del z["sum_empty"]
    for cid, (fn, v) in z.items():
        omit = "omitnan" in cid
        args = [In("v"), "omitnan"] if omit else [In("v")]
        py = py_sum if fn == "sum" else py_mean
        cs.append(
            Case(
                cid,
                matlab=fn,
                args=args,
                py_args=[In("v"), None, "omitnan"] if omit else [In("v")],
                python=py,
                inputs={"v": Value(lambda v=v: np.array(v, dtype=np.float64))},
            )
        )
    cs.append(
        Case(
            "sum_slices_negzeros_chunked",
            matlab="sum",
            args=[In("v"), 2],
            py_args=[In("v"), 1],
            python=lambda v, ax: R.sum(v, ax, keepdims=True, matlab_layout=True),
            inputs={"v": Value(lambda: np.full((3, 2100), NZ))},
        )
    )
    return cs


_QB = {
    "single": (np.uint32, 0x7F800000, 0x00400000, 31),
    "double": (np.uint64, 0x7FF0 << 48, 1 << 51, 63),
}


def nan_wide(seed: int, shape, cls: str = "single", every: int = 9, infs: int = 0, snan=False):
    """``wide`` values (MATLAB's form) with NaNs at about one in ``every`` places, each with its
    own payload: its 1-based MATLAB index, the sign bit set at odd ones, quiet (with ``snan``
    every third signaling); ``infs`` pairs of Inf and -Inf at other places (Inf - Inf)"""
    v = np.array(wide(seed, shape, cls), order="F")
    flat = v.reshape(-1, order="F")  # a view
    n = flat.size
    b = _bits64(seed + 7919, n)
    pos = np.flatnonzero(b % np.uint64(every) == 0)
    ut, exp, quiet, sb = _QB[cls]
    t = pos.astype(np.uint64) + 1
    bits = np.uint64(exp) | t | ((t & np.uint64(1)) << np.uint64(sb))
    sig = (t % np.uint64(3) == 0) if snan else np.zeros(t.size, bool)
    bits = np.where(sig, bits, bits | np.uint64(quiet))
    flat[pos] = bits.astype(ut).view(flat.dtype)
    if infs:
        rest = np.flatnonzero(b % np.uint64(every) != 0)
        k = rest[(_bits64(seed + 104729, rest.size) % np.uint64(rest.size)).astype(np.int64)]
        k = np.unique(k)[: 2 * infs]
        flat[k[0::2]] = np.inf
        flat[k[1::2]] = -np.inf
    return v


def _payload_cases():
    """each output's NaN as MATLAB gives it (X05): two or more payloads per output, in the vector
    rule's regimes (index order, four lanes and their tails, 4,096-blocks at one thread), the
    slice rule and its chunks, merged and successive dimensions, single and double, 'double'
    sums of single data, mean, 'omitnan' with Inf - Inf, cumsum, var, signaling NaNs"""
    cs = []
    seed = 900
    for cls in ("single", "double"):
        outs = (("sum", ()), ("mean", ()))
        if cls == "single":
            outs += (("sum", ("double",)),)
        for n, every, infs in ((5, 2, 0), (43, 6, 1), (47, 5, 0), (50, 4, 1), (51, 7, 0),
                               (4097, 300, 2), (50001, 4000, 0), (89003, 9000, 0),
                               (90114, 20000, 1), (90122, 3, 0), (90155, 5000, 0)):  # fmt: skip
            seed += 1
            make = Value(
                lambda s=seed, n=n, c=cls, e=every, i=infs: nan_wide(s, (1, n), c, e, i).ravel()
            )
            for fn, opt in outs:
                tag = f"{fn}_payloads_row_{cls}_{n}" + ("_double_out" if opt else "")
                args = ", " + ", ".join(f"'{o}'" for o in opt) if opt else ""
                cs.append(
                    Case(
                        tag,
                        code=_one_thread(f"{fn}(in{{1}}{args})"),
                        args=[In("x")],
                        py_args=[In("x"), None, *opt],
                        python=py_sum if fn == "sum" else py_mean,
                        inputs={"x": make},
                    )
                )
    shapes = (
        ((44, 30, 70), "single", 9),
        ((130, 80, 9), "double", 11),
        ((6, 2100, 5), "single", 300),
        ((11, 1500, 2), "double", 200),
        ((7, 300, 50), "single", 40),
    )
    for shape, cls, every in shapes:
        seed += 1
        tag = "x".join(map(str, shape))
        make = Value(lambda s=seed, sh=shape, c=cls, e=every: nan_wide(s, sh, c, e, 3, True))
        for d in _DIMS:
            margs, pargs = _dim_args(d)
            outs = [("sum", ()), ("mean", ())]
            if cls == "single":
                outs.append(("sum", ("double",)))
            for fn, opt in outs:
                o = "".join(f", '{x}'" for x in opt)
                cs.append(
                    Case(
                        f"{fn}_payloads_{cls}_{tag}_{_dname(d)}" + ("_double_out" if opt else ""),
                        code=_one_thread(f"{fn}(in{{1}}, in{{2}}{o})"),
                        args=[Vol(In("x")), *margs],
                        py_args=[Vol(In("x")), *pargs, *opt],
                        python=py_sum if fn == "sum" else py_mean,
                        out_roles=("vol",),
                        inputs={"x": make},
                    )
                )
    # 4,096-blocks in a merged group, and in a group after another
    big = (
        ((300, 300, 2), "single", (1, 2)),
        ((300, 300, 2), "double", (1, 2)),
        ((300, 300, 2), "single", "all"),
        ((9, 3, 12000), "single", (1, 3)),
    )
    for shape, cls, d in big:
        seed += 1
        tag = "x".join(map(str, shape))
        make = Value(lambda s=seed, sh=shape, c=cls: nan_wide(s, sh, c, 9000, 2))
        margs, pargs = _dim_args(d)
        opts = ((), ("double",)) if cls == "single" else ((),)
        for opt in opts:
            o = "".join(f", '{x}'" for x in opt)
            cs.append(
                Case(
                    f"sum_payloads_{cls}_{tag}_{_dname(d)}" + ("_double_out" if opt else ""),
                    code=_one_thread(f"sum(in{{1}}, in{{2}}{o})"),
                    args=[Vol(In("x")), *margs],
                    py_args=[Vol(In("x")), *pargs, *opt],
                    python=py_sum,
                    out_roles=("vol",),
                    inputs={"x": make},
                )
            )
    # many vectors (17 or more in 39,600 or more elements): 4,096-blocks at any length, values
    # and payloads; a last block of up to 43 values in index order
    many = (
        ((5000, 17, 2), "single", 1, 0),
        ((5000, 16), "single", 1, 0),
        ((4101, 20), "double", 1, 0),
        ((8193, 18), "single", 1, 0),
        ((50, 800), "single", 9, 0),
        ((51, 800), "double", 9, 0),
        ((4140, 20), "single", 37, 0),
        ((4140, 20), "double", 37, 1),
        ((90, 21, 23), "single", 11, 0),
    )
    for shape, cls, every, infs in many:
        seed += 1
        tag = "x".join(map(str, shape))
        if every == 1:
            make = Value(lambda s=seed, sh=shape, c=cls: wide(s, sh, c))
        else:
            make = Value(lambda s=seed, sh=shape, c=cls, e=every, i=infs: nan_wide(s, sh, c, e, i))
        for opt in ((), ("double",)) if cls == "single" else ((),):
            o = "".join(f", '{x}'" for x in opt)
            for d in (1, (1, 2)) if len(shape) == 3 and every > 1 else (1,):
                margs, pargs = _dim_args(d)
                cs.append(
                    Case(
                        f"sum_many_{cls}_{tag}_{_dname(d)}" + ("_double_out" if opt else ""),
                        code=_one_thread(f"sum(in{{1}}, in{{2}}{o})"),
                        args=[Vol(In("x")), *margs],
                        py_args=[Vol(In("x")), *pargs, *opt],
                        python=py_sum,
                        out_roles=("vol",),
                        inputs={"x": make},
                    )
                )
    # 'omitnan': only Inf - Inf makes a NaN, MATLAB's
    for cls in ("single", "double"):
        seed += 1
        make = Value(lambda s=seed, c=cls: nan_wide(s, (12, 40, 30), c, 5, 40))
        for d in (1, 2, 3, "all"):
            margs, pargs = _dim_args(d)
            for fn in ("sum", "mean"):
                cs.append(
                    Case(
                        f"{fn}_payloads_omitnan_{cls}_{_dname(d)}",
                        code=f"out = {{{fn}(in{{1}}, in{{2}}, 'omitnan')}};",
                        args=[Vol(In("x")), *margs],
                        py_args=[Vol(In("x")), *pargs, "omitnan"],
                        python=py_sum if fn == "sum" else py_mean,
                        out_roles=("vol",),
                        inputs={"x": make},
                    )
                )
    # cumsum: the running sum first along dimension 1 and rows, each term first after it
    for cls in ("single", "double"):
        seed += 1
        make = Value(lambda s=seed, c=cls: nan_wide(s, (20, 9, 11), c, 7, 4, True))
        for d in (1, 2, 3):
            cs.append(
                Case(
                    f"cumsum_payloads_{cls}_d{d}",
                    matlab="cumsum",
                    args=[Vol(In("x")), Dim(d)],
                    python=py_cumsum,
                    out_roles=("vol",),
                    inputs={"x": make},
                )
            )
        seed += 1
        cs.append(
            Case(
                f"cumsum_payloads_row_{cls}",
                matlab="cumsum",
                args=[In("v"), 2],
                py_args=[In("v"), -1],
                python=py_cumsum,
                inputs={
                    "v": Value(lambda s=seed, c=cls: nan_wide(s, (1, 60), c, 4, 1, True).ravel())
                },
            )
        )
    # var: abs(x - mean).^2 summed (abs clears the sign of a NaN)
    for cls in ("single", "double"):
        seed += 1
        make = Value(lambda s=seed, c=cls: nan_wide(s, (44, 30, 7), c, 13, 2))
        for d in (1, 2, 3):
            cs.append(
                Case(
                    f"var_payloads_{cls}_d{d}",
                    code="out = {var(in{1}, 0, in{2})};",
                    args=[Vol(In("x")), Dim(d)],
                    python=lambda x, ax: S.var(x, 0, ax, keepdims=True),
                    out_roles=("vol",),
                    inputs={"x": make},
                )
            )
    # one value: a copy (a signaling NaN stays signaling)
    for cls in ("single", "double"):
        seed += 1
        make = Value(lambda s=seed, c=cls: nan_wide(s, (6, 5), c, 2, 0, True))
        for d in (1, 2, 3):
            cs.append(
                Case(
                    f"sum_payloads_singleton_{cls}_d{d}",
                    code="x = in{1}; out = {sum(x(1, :), in{2}), sum(x(:, 1), in{2})};",
                    args=[In("x"), d],
                    nargout=2,
                    python=_singleton_sums,
                    inputs={"x": make},
                )
            )
    # a last 4,096-block of 4 to 43 values: in index order (values), few and many vectors
    for cls in ("single", "double"):
        for n, k in ((90122, 1), (90155, 1), (8235, 17), (90155, 17)):
            seed += 1
            make = Value(lambda s=seed, n=n, k=k, c=cls: last_block_wide(s, n, k, c))
            outs = ((), ("double",)) if cls == "single" else ((),)
            for opt in outs:
                o = "".join(f", '{x}'" for x in opt)
                cs.append(
                    Case(
                        f"sum_short_last_block_{cls}_{n}x{k}" + ("_double_out" if opt else ""),
                        code=_one_thread(f"sum(in{{1}}, 1{o})"),
                        args=[In("x")],
                        py_args=[In("x"), *opt],
                        python=lambda x, *o: R.sum(
                            x,
                            0,
                            outtype=o[0] if o else "default",
                            keepdims=True,
                            matlab_layout=True,
                        ),  # fmt: skip
                        inputs={"x": make},
                    )
                )
    return cs


def last_block_wide(seed: int, n: int, k: int, cls: str) -> np.ndarray:
    """(n, k) in MATLAB's form: small values, and in the last 4,096-block values m * 2^e of a
    wide range (e in [-60, 60]), whose sums differ between index order and four lanes"""
    x = np.array(wide(seed, (n, k), cls), order="F")
    x *= x.dtype.type(2.0**-40)
    m = n % R.BLOCK
    b = _bits64(seed + 31, m * k).reshape(m, k)
    e = (b % np.uint64(121)).astype(np.int64) - 60
    v = np.ldexp(((b >> np.uint64(40)) % np.uint64(2**20)).astype(np.float64) - 2**19, e)
    x[n - m :, :] = v.astype(x.dtype)
    return x


def _singleton_sums(x, d):
    """sum(x(1, :), d) and sum(x(:, 1), d) of a matrix in MATLAB's form"""
    return tuple(
        R.sum(v[..., None], d - 1, keepdims=True, matlab_layout=True)[..., 0]
        for v in (x[:1, :], x[:, :1])
    )


CASES = (
    _vector_cases()
    + _volume_cases()
    + _int_cases()
    + _cumsum_cases()
    + _index_cases()
    + _pair_cases()
    + _zero_cases()
    + _payload_cases()
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- a scalar model of MATLAB's order (unit tests) ---------------------------------------------


def _model_vector(v: np.ndarray, dt) -> float:
    """the vector rule, one addition at a time in Python"""
    T = dt.type
    n = v.size
    if n == 0:
        return T(0)
    if n <= R.SEQUENTIAL_UPTO:
        s = T(v[0])
        for x in v[1:]:
            s = T(s + T(x))
        return s

    def lanes(w):
        acc = [T(w[k]) if k < w.size else T(NZ) for k in range(4)]
        for t in range(4, w.size):
            acc[t % 4] = T(acc[t % 4] + T(w[t]))
        return T(T(T(acc[0] + acc[1]) + acc[2]) + acc[3])

    if n < R.BLOCKS_FROM:
        return lanes(v)
    s = None
    for i in range(0, n, R.BLOCK):
        w = v[i : i + R.BLOCK]
        b = lanes(w) if w.size > R.SEQUENTIAL_UPTO else _model_vector(w, dt)  # a short last block
        s = b if s is None else T(s + b)
    return s


@pytest.mark.parametrize(
    "n", [1, 2, 43, 44, 45, 46, 47, 1000, 4097, 88999, 89000, 89001, 90122, 90155, 95000]
)
@pytest.mark.parametrize("cls", ["single", "double"])
def test_vector_rule_equals_the_scalar_model(n, cls):
    v = wide(n + (7 if cls == "single" else 11), (n,), cls)
    dt = v.dtype
    with np.errstate(all="ignore"):
        ref = _model_vector(v, dt)
    assert np.asarray(R.sum(v)).tobytes() == np.asarray(ref, dtype=dt).tobytes()


def _model_slices(x2: np.ndarray, chunked: bool, dt) -> np.ndarray:
    """out(i) = x(i, 1) + x(i, 2) + ..., in chunks of 1,024 when chunked"""
    T = dt.type
    m, k = x2.shape
    out = np.empty(m, dt)
    for i in range(m):
        total = None
        step = R.SLICE_CHUNK if chunked else k
        for c0 in range(0, k, step):
            p = T(x2[i, c0])
            for j in range(c0 + 1, min(c0 + step, k)):
                p = T(p + T(x2[i, j]))
            total = p if total is None else T(total + p)
        out[i] = total
    return out


@pytest.mark.parametrize(
    ("shape", "chunked"),
    [
        ((6, 1000), False),
        ((6, 2047), False),
        ((6, 2048), True),
        ((16, 1250), False),
        ((16, 1260), True),
        ((12, 2000), False),
        ((19, 1100), True),
        ((2, 3000), True),
    ],
)
def test_slice_rule_and_its_chunks(shape, chunked):
    x = wide(int(np.prod(shape)), shape, "single")
    r = R.sum(x, 1, matlab_layout=True)
    with np.errstate(all="ignore"):
        ref = _model_slices(x, chunked, x.dtype)
    assert r.tobytes() == ref.tobytes()
    assert R._slice_chunked(shape[1], shape[0], 1, x.size) is chunked


@pytest.mark.parametrize("cls", ["single", "double"])
@pytest.mark.parametrize("axis", [None, 0, 1, 2, (0, 2), (1, 2), (0, 1)])
def test_layouts_agree(cls, axis):
    """the port's (z, y, x) array and MATLAB's own layout give the same sums"""
    m = wide(91 + (axis is None), (44, 30, 70), cls)
    from petakit5d.core.arrays import from_matlab, to_matlab
    from petakit5d.core.axes import dim_to_axis

    p = from_matlab(m)
    if axis is None:
        pax = None
    else:
        dims = (axis,) if isinstance(axis, int) else axis
        pax = tuple(dim_to_axis(d + 1) for d in dims)
    a = R.sum(m, axis, keepdims=True, matlab_layout=True)
    b = R.sum(p, pax, keepdims=True)
    assert np.asarray(a).tobytes() == np.ascontiguousarray(to_matlab(np.asarray(b))).tobytes() or (
        np.ndim(a) == 0 and np.asarray(a).tobytes() == np.asarray(b).tobytes()
    )


@pytest.mark.parametrize("fn", ["sum", "mean"])
def test_results_do_not_depend_on_the_thread_count(fn):
    x = wide(5, (4096 * 300 + 17,), "single")
    f = getattr(R, fn)
    r1 = f(x, num_threads=1)
    for nt in (2, 3, 8):
        assert np.asarray(f(x, num_threads=nt)).tobytes() == np.asarray(r1).tobytes()


def test_mean_divides_in_double():
    """above 2^24 elements, single(double(s)/n) differs from s/single(n)"""
    x = np.full(2**24 + 1, 3.0, np.float32)
    x[0] = np.float32(7e6)
    s = R.sum(x)
    assert R.mean(x) == np.float32(np.float64(s) / x.size)


def test_nan_results_are_matlabs_nan():
    v = np.array([1.0, np.inf, -np.inf, 2.0])
    r = R.sum(v)
    assert np.asarray(r).tobytes() == np.float64(dtypes.NAN).tobytes()
    q = np.array([1.0, np.float64(np.nan)])  # a positive NaN keeps its payload
    assert np.asarray(R.sum(q)).view(np.uint64) == np.float64(np.nan).view(np.uint64)


def _winner_vector(pos, L: int, o, many: bool = False) -> int | None:
    """the 0-based place (of ``pos``, the NaNs among zeros) whose NaN a vector-rule sum of L
    values keeps, from the operand order alone (``reduce._NAN_ORDER``'s fields); ``many``: in
    MATLAB's order for many vectors (4,096-blocks at any length)"""

    def first(ps):
        return min(ps) if ps else None

    def lanes(ps, m, tail_first, swap):
        main = m - m % 4
        win = []
        for j in range(4):
            w = first([p for p in ps if p < main and p % 4 == j])
            t = [p for p in ps if p >= main and p - main == j]
            if t and (tail_first or w is None):
                w = t[0]
            win.append(w)
        order = (1, 0, 2, 3) if swap else (0, 1, 2, 3)
        return next((win[j] for j in order if win[j] is not None), None)

    if L == 1:
        return first(pos)
    if L <= R.SEQUENTIAL_UPTO:
        return first(pos)
    if L < R.BLOCKS_FROM and not many:
        return lanes(pos, L, o.tail, False)
    wins = []
    for b0 in range(0, L, R.BLOCK):
        m = min(R.BLOCK, L - b0)
        ps = [p - b0 for p in pos if b0 <= p < b0 + m]
        last = b0 + m == L and L % R.BLOCK
        if m <= R.SEQUENTIAL_UPTO:
            w = first(ps)
        else:
            w = lanes(ps, m, o.block_tail if last else False, o.block_swap)
        wins.append(None if w is None else w + b0)
    wins = [w for w in wins if w is not None]
    if not wins:
        return None
    return wins[-1] if o.blocks_term_first else wins[0]


def _winner_slices(pos, n: int, chunked: bool, term_first: bool) -> int | None:
    """the same for the slice rule along n slices"""
    size = R.SLICE_CHUNK if chunked and n > R.SLICE_CHUNK else n
    wins = []
    for c0 in range(0, n, size):
        ps = sorted(p for p in pos if c0 <= p < c0 + size)
        if ps:
            wins.append(ps[-1] if term_first and len(ps) else ps[0])
    return wins[0] if wins else None


def _with_nans(shape, places, cls):
    """zeros with NaNs at MATLAB linear places, each with its place + 1 as payload"""
    ut, exp, quiet, _ = _QB[cls]
    x = np.zeros(int(np.prod(shape)), np.float32 if cls == "single" else np.float64)
    for p in places:
        x[p] = ut(exp | quiet | (p + 1)).view(x.dtype)
    return x.reshape(shape, order="F")


def _bits(r) -> int:
    r = np.asarray(r).reshape(1)
    return int(r.view(np.uint32 if r.dtype == np.float32 else np.uint64)[0])


def _payload_id(r, cls, outtype):
    b = _bits(r)
    if np.asarray(r).dtype == np.float32:
        return b & 0x3FFFFF
    return (b >> 29) & 0x3FFFFF if cls == "single" else b & 0xFFFFFFFFFFFF


_COMBOS = (("single", "default"), ("double", "default"), ("single", "double"))


@pytest.mark.parametrize("cls,outtype", _COMBOS)
@pytest.mark.parametrize("n", [5, 43, 44, 50, 51, 4097, 89000, 89001, 89003, 90114])
def test_vector_rule_payloads_follow_the_operand_order(cls, outtype, n):
    """the NaN a sum keeps, against a model of the operand order alone (zeros elsewhere), in
    rows, in columns of either layout and at 1 and 3 threads"""
    o = R._NAN_ORDER[
        (
            np.dtype(np.float32 if cls == "single" else np.float64),
            np.dtype(np.float32 if cls == "single" and outtype == "default" else np.float64),
        )
    ]
    b = _bits64(n * 7 + len(cls) + len(outtype), 120)
    special = [0, 1, 2, 3, n - 1, n - 2, n - 3, n - n % 4, max(n - 4096, 0), min(4096, n - 1)]
    for t in range(20):
        k = 2 + int(b[t] % np.uint64(3))
        ps = sorted({int(v % np.uint64(n)) for v in b[20 + 4 * t : 20 + 4 * t + k]})
        if t < 5:
            ps = sorted({special[(t + j) % len(special)] for j in range(3)} | set(ps[:1]))
        want = _winner_vector(ps, n, o) + 1
        row = _with_nans((1, n), ps, cls)
        col = np.ascontiguousarray(np.concatenate([row.T, np.zeros((n, 1), row.dtype)], 1))
        got = [
            R.sum(row.ravel(), outtype=outtype, num_threads=1),
            R.sum(row, 1, outtype=outtype, matlab_layout=True, num_threads=3)[0],
            R.sum(col, 0, outtype=outtype, matlab_layout=True, num_threads=3)[0],
            R.sum(np.asfortranarray(col), 0, outtype=outtype, matlab_layout=True)[0],
            R.sum(col[None], -2, outtype=outtype, num_threads=1)[0, 0],
        ]
        for g in got:
            assert _payload_id(g, cls, outtype) == want, (ps, want)


@pytest.mark.parametrize("cls,outtype", _COMBOS)
@pytest.mark.parametrize(
    "n,k", [(50, 792), (50, 791), (2475, 16), (2475, 17), (45, 880), (4101, 17), (4140, 20)]
)
def test_many_vector_payloads_follow_matlabs_order_for_many_vectors(cls, outtype, n, k):
    """17 or more columns in 39,600 or more elements: MATLAB's 4,096-blocks at any length
    (``MANY_VECTORS_FROM``, ``MANY_NUMEL_FROM``), in either layout"""
    o = R._NAN_ORDER[
        (
            np.dtype(np.float32 if cls == "single" else np.float64),
            np.dtype(np.float32 if cls == "single" and outtype == "default" else np.float64),
        )
    ]
    many = k >= R.MANY_VECTORS_FROM and n * k >= R.MANY_NUMEL_FROM
    b = _bits64(n * 3 + k, 60)
    tails = [n - 1, n - 2, n - 3, 4096, 4097, 1, 2, 5]
    for t in range(15):
        ps = sorted({tails[t % 8], tails[(t * 3 + 1) % 8], int(b[t] % np.uint64(n))} - {n, n + 1})
        ps = [p for p in ps if p < n]
        want = _winner_vector(ps, n, o, many) + 1
        x = np.zeros((n, k), np.float32 if cls == "single" else np.float64, order="F")
        x[:, 0] = _with_nans((n,), ps, cls)
        for got in (
            R.sum(x, 0, outtype=outtype, matlab_layout=True)[0],
            R.sum(np.ascontiguousarray(x)[None], -2, outtype=outtype)[0, 0],
        ):
            assert _payload_id(got, cls, outtype) == want, (ps, want)


@pytest.mark.parametrize("cls", ["single", "double"])
@pytest.mark.parametrize("n,k", [(90155, 1), (90122, 2), (90155, 2), (8235, 17)])
def test_a_short_last_block_is_summed_in_index_order(cls, n, k):
    """a last 4,096-block of 4 to 43 values: in index order (four lanes give other bits here),
    on the kernel paths too (all elements, a column) and in either layout"""
    x = last_block_wide(5 + n + k, n, k, cls)
    dt = x.dtype
    T = dt.type

    def lanes4(w):
        acc = [T(NZ)] * 4
        for t, v in enumerate(w):
            acc[t % 4] = T(acc[t % 4] + v)
        return T(T(T(acc[0] + acc[1]) + acc[2]) + acc[3])

    def blocks(v, last):
        s = None
        for i in range(0, v.size, R.BLOCK):
            w = v[i : i + R.BLOCK]
            b = lanes4(w) if w.size == R.BLOCK else last(w)
            s = b if s is None else T(s + b)
        return s

    with np.errstate(all="ignore"):
        want = blocks(x[:, 0], lambda w: _model_vector(w, dt))  # (in index order)
        lanes_last = blocks(x[:, 0], lanes4)
    if k == 1:
        assert _bits(R.sum(x)) == _bits(want) != _bits(lanes_last)
        assert _bits(R.sum(x, outtype="double")) == _bits(R.sum(x.astype(np.float64)))
    got = R.sum(x, 0, matlab_layout=True)
    assert _bits(got[0]) == _bits(want)
    assert _bits(R.sum(np.ascontiguousarray(x)[None], -2)[0, 0]) == _bits(want)


@pytest.mark.parametrize("cls", ["single", "double"])
def test_many_vectors_of_more_than_4096_values_are_summed_in_blocks(cls):
    """values: 17 columns of 5,000 add as two blocks (4,096 + 904), 16 in four lanes"""
    for k in (16, 17):
        x = wide(77 + k, (5000, k), cls)
        got = R.sum(x, 0, matlab_layout=True)
        dt = x.dtype.type
        for c in (0, k - 1):
            v = x[:, c]
            if k >= R.MANY_VECTORS_FROM:
                want = dt(_model_vector(v[:4096], x.dtype) + _model_vector(v[4096:], x.dtype))
            else:
                want = _model_vector(v, x.dtype)
            assert _bits(got[c]) == _bits(want)


@pytest.mark.parametrize("cls,outtype", _COMBOS)
@pytest.mark.parametrize("n", [50, 1024, 1025, 2100, 3000])
def test_slice_rule_payloads_follow_the_operand_order(cls, outtype, n):
    o = R._NAN_ORDER[
        (
            np.dtype(np.float32 if cls == "single" else np.float64),
            np.dtype(np.float32 if cls == "single" and outtype == "default" else np.float64),
        )
    ]
    b = _bits64(n * 13 + len(cls) + len(outtype), 100)
    for t in range(20):
        ps = sorted({int(v % np.uint64(n)) for v in b[4 * t : 4 * t + 2 + t % 3]})
        chunked = R._slice_chunked(n, 2, 1, 2 * n)
        want = _winner_slices(ps, n, chunked, o.slice_term_first) + 1
        x = _with_nans((2, n), [2 * p for p in ps], cls)
        for got in (
            R.sum(x, 1, outtype=outtype, matlab_layout=True)[0],
            R.sum(np.ascontiguousarray(x), -1, outtype=outtype, num_threads=3)[0],
        ):
            assert _payload_id(got, cls, outtype) // 2 == want - 1, (ps, want)


def test_payloads_of_merged_and_successive_dims():
    """[1 2] is one vector of ny * nx values per plane; [1 3] reduces y first, then z, of the
    sums so far (of the sum's class: 'double' sums of single data take double's order there)"""
    for cls, outtype in _COMBOS:
        x = _with_nans((7, 9, 5), [3, 10, 70, 71, 200, 214], cls)  # MATLAB's (y, x, z)
        dt = np.dtype(np.float32 if cls == "single" and outtype == "default" else np.float64)
        o = R._NAN_ORDER[(x.dtype, dt)]
        r = R.sum(x, (0, 1), outtype=outtype, keepdims=True, matlab_layout=True)
        for z in range(5):
            ps = [p - 63 * z for p in (3, 10, 70, 71, 200, 214) if 63 * z <= p < 63 * (z + 1)]
            if ps:
                assert (
                    _payload_id(r[0, 0, z], cls, outtype) == _winner_vector(ps, 63, o) + 63 * z + 1
                )
        r = R.sum(x, (0, 2), outtype=outtype, keepdims=True, matlab_layout=True)
        # (y, x, z) of the places: 3 (3, 0, 0); 10 (3, 1, 0), 70 (0, 1, 1), 71 (1, 1, 1),
        # 200 (4, 1, 3); 214 (4, 3, 3): x = 1 keeps z = 0's NaN, the running sum first
        assert _payload_id(r[0, 0, 0], cls, outtype) == 4
        assert _payload_id(r[0, 1, 0], cls, outtype) == 11
        assert _payload_id(r[0, 3, 0], cls, outtype) == 215


def test_payloads_are_x86s_on_every_cpu():
    """a signaling NaN meets a quiet one: x86 keeps the first operand, quieted (ARM would keep
    the signaling one); Inf - Inf is MATLAB's NaN (ARM's default NaN is positive)"""
    for cls in ("single", "double"):
        ut, exp, quiet, _ = _QB[cls]
        dt = np.float32 if cls == "single" else np.float64
        q = ut(exp | quiet | 5).view(dt)
        s = ut(exp | 9).view(dt)
        v = np.array([q, s, 1.0], dtype=dt)
        assert _bits(R.sum(v)) == exp | quiet | 5
        assert _bits(R.sum(v[::-1])) == exp | quiet | 9
        v = np.array([np.inf, -np.inf, 2.0], dtype=dt)
        for r in (R.sum(v), R.sum(v, omitnan=True), R.cumsum(v, 0)[1], R.mean(v)):
            assert np.asarray(r).astype(np.float64).tobytes() == np.float64(dtypes.NAN).tobytes()


def test_one_value_is_a_copy():
    for cls in ("single", "double"):
        ut, exp, _, _ = _QB[cls]
        s = ut(exp | 9).view(np.float32 if cls == "single" else np.float64)
        x = np.full((1, 3), s)
        assert np.asarray(R.sum(x[:, :1])).view(ut) == ut(exp | 9)  # sum(x(1))
        assert np.all(np.asarray(R.sum(x, 0, keepdims=True)).view(ut) == ut(exp | 9))
        assert np.asarray(R.cumsum(x, 1)).view(ut)[0, 0] == ut(exp | 9)
        assert np.asarray(R.cumsum(x, 1)).view(ut)[0, 1] != ut(exp | 9)  # added: quieted


def test_cumsum_payloads():
    """along MATLAB's dimension 1 the first NaN stays; along a later one the last NaN term"""
    x = _with_nans((3, 4), [1, 2, 4, 7, 10], "double")  # (1-based) 2, 3, 5, 8, 11
    ids = np.vectorize(lambda r: int(np.float64(r).view(np.uint64)) & 0xFFFF)
    d1 = R.cumsum(x, 0, matlab_layout=True)
    assert ids(d1[:, 0]).tolist() == [0, 2, 2] and ids(d1[:, 1]).tolist() == [0, 5, 5]
    d2 = R.cumsum(x, 1, matlab_layout=True)
    assert ids(d2[1]).tolist() == [2, 5, 8, 11] and ids(d2[2]).tolist() == [3, 3, 3, 3]
    p = np.ascontiguousarray(x)  # an image (y, x): axis -1 is MATLAB's dimension 2
    assert ids(R.cumsum(p, -1)[1]).tolist() == [2, 5, 8, 11]
    assert ids(R.cumsum(p, -2)[:, 1]).tolist() == [0, 5, 5]


def test_omitnan_lead_counts_as_plus_zero():
    assert np.signbit(R.sum(np.array([np.nan, NZ]), omitnan=True)) == False  # noqa: E712
    assert np.signbit(R.sum(np.array([NZ, np.nan, NZ]), omitnan=True))
    assert np.signbit(R.cumsum(np.array([np.nan, NZ]), 0, omitnan=True)[0]) == False  # noqa: E712


def test_signed_zero_ties_keep_the_first():
    m, i = R.max(np.array([0.0, NZ]))
    assert not np.signbit(m) and i == (0,)
    m, i = R.max(np.array([NZ, 0.0]))
    assert np.signbit(m) and i == (0,)
    assert np.signbit(R.maximum(NZ, 0.0)) and not np.signbit(R.maximum(0.0, NZ))


def test_first_maximum_in_matlabs_order():
    """(z, y, x) = (0, 1, 0) and (0, 0, 1) hold the maximum: MATLAB's first is y fastest"""
    v = np.zeros((2, 3, 4))
    v[0, 1, 0] = 5
    v[0, 0, 1] = 5
    assert R.max(v)[1] == (0, 1, 0)
    v[0, 1, 0] = 4
    assert R.max(v)[1] == (0, 0, 1)


def test_default_axis():
    assert R.default_axis(np.zeros((1, 3)), matlab_layout=True) == 1
    assert R.default_axis(np.zeros((4, 3)), matlab_layout=True) == 0
    assert R.default_axis(np.zeros((5, 1, 3))) == 2  # MATLAB 1 x 3 x 5: dimension 2, here x
    assert R.default_axis(np.zeros((5, 2, 3))) == 1


def test_errors():
    with pytest.raises(ValueError):
        R.sum(np.zeros(3), outtype="single")
    with pytest.raises(TypeError):
        R.sum(np.zeros(3, np.int64))
    with pytest.raises(TypeError):
        R.sum(np.array(["a"]))
    with pytest.raises(ValueError):
        R.max(np.zeros(0))
    with pytest.raises(ValueError):
        R.unique(np.zeros(3), occurrence="middle")


def test_native_integer_sums_saturate_in_order():
    assert R.sum(np.array([100, 100, -50], np.int8), outtype="native") == 77
    assert R.sum(np.array([-50, 100, 100], np.int8), outtype="native") == 127
