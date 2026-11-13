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


CASES = (
    _vector_cases()
    + _volume_cases()
    + _int_cases()
    + _cumsum_cases()
    + _index_cases()
    + _pair_cases()
    + _zero_cases()
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
        b = lanes(v[i : i + R.BLOCK])
        s = b if s is None else T(s + b)
    return s


@pytest.mark.parametrize("n", [1, 2, 43, 44, 45, 46, 47, 1000, 4097, 88999, 89000, 89001, 95000])
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
