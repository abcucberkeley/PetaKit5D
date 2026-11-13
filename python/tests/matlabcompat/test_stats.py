"""Tests of petakit5d.matlabcompat.stats: MATLAB's var, std, median, prctile and mode against
MATLAB R2024b (oracle cases) and their rules (unit tests).

Inputs as in test_reduce (exact integers times powers of two, so the order of the sums shows);
sums of 89,000 elements or more run MATLAB at one thread. Signed zeros at a median or mode of
equal values: MATLAB's choice of -0 or +0 depends on its selection code for larger arrays (decision
record 0005); the cases here hold either no such tie or a small one, where MATLAB takes the
stable order. No MATLAB test covers these built-ins alone.
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


def py_var(x, w=0, dim=None, *rest):
    return S.var(x, w, dim, omitnan="omitnan" in rest, keepdims=dim is not None)


def py_std(x, w=0, dim=None, *rest):
    return S.std(x, w, dim, omitnan="omitnan" in rest, keepdims=dim is not None)


def py_median(x, dim=None, *rest):
    return S.median(x, dim, omitnan="omitnan" in rest, keepdims=dim is not None)


def py_prctile(x, p, dim=None):
    return S.prctile(x, p, dim)


def _var_cases():
    cs = []
    seed = 700
    for cls in ("single", "double"):
        for n in (1, 2, 3, 43, 44, 1000, 50000):
            seed += 1
            make = Value(lambda s=seed, n=n, c=cls: wide(s, (1, n), c).ravel())
            for fn, py in (("var", py_var), ("std", py_std)):
                for w in (0, 1):
                    cs.append(
                        Case(
                            f"{fn}_{cls}_{n}_w{w}",
                            matlab=fn,
                            args=[In("x"), w],
                            python=py,
                            inputs={"x": make},
                        )
                    )
        for n in (89000, 100003):
            seed += 1
            make = Value(lambda s=seed, n=n, c=cls: wide(s, (1, n), c).ravel())
            cs.append(
                Case(
                    f"std_{cls}_{n}_one_thread",
                    code=_one_thread("std(in{1})"),
                    args=[In("x")],
                    python=py_std,
                    inputs={"x": make},
                )
            )
        vol = Value(lambda s=seed, c=cls: wide(s + 50, (20, 9, 11), c))
        for d in (1, 2, 3):
            cs.append(
                Case(
                    f"std_{cls}_vol_d{d}",
                    matlab="std",
                    args=[Vol(In("v")), 0, Dim(d)],
                    python=py_std,
                    out_roles=("vol",),
                    inputs={"v": vol},
                )
            )
        cs.append(
            Case(
                f"std_{cls}_vol_d12",
                matlab="std",
                args=[Vol(In("v")), 0, [1, 2]],
                py_args=[Vol(In("v")), 0, (-2, -1)],
                python=py_std,
                out_roles=("vol",),
                inputs={"v": vol},
            )
        )
        cs.append(
            Case(
                f"std_{cls}_vol_all",
                matlab="std",
                args=[Vol(In("v")), 0, "all"],
                py_args=[Vol(In("v")), 0, None],
                python=py_std,
                inputs={"v": vol},
            )
        )
    nan = dtypes.NAN
    for name, v in (("nan", [2.0, nan, 4.0, 7.0]), ("negzero", [NZ]), ("equal", [1.5, 1.5])):
        for omit in ((), ("omitnan",)):
            cs.append(
                Case(
                    f"var_{name}{'_omitnan' if omit else ''}",
                    matlab="var",
                    args=[In("v"), 0, 2, *omit] if omit else [In("v")],
                    py_args=[In("v"), 0, None, *omit] if omit else [In("v")],
                    python=py_var,
                    inputs={"v": Value(lambda v=v: np.array(v))},
                )
            )
    cs.append(
        Case(
            "var_empty",
            matlab="var",
            args=[In("v")],
            python=py_var,
            inputs={"v": Value(lambda: np.zeros((1, 0)))},
        )
    )
    cs.append(
        Case(
            "std_uint16_errors",
            matlab="std",
            args=[In("v")],
            python=py_std,
            expect_error="MATLAB:var:integerClass",
            py_error=R.MatlabError,
            inputs={"v": Value(lambda: np.array([1, 2, 3], np.uint16))},
        )
    )
    # crop_subregion_by_intensity: mean and std of a single region, at one thread
    cs.append(
        Case(
            "mean_std_region_single_one_thread",
            code=_one_thread("deal(mean(in{1}(:)), std(single(in{1}(:))))", 2),
            args=[Vol(In("r"))],
            nargout=2,
            python=lambda r: (R.mean(r), S.std(r)),
            inputs={"r": Value(lambda: wide(771, (64, 2100, 3), "single"))},
        )
    )
    return cs


def _median_vectors():
    nan = dtypes.NAN
    return {
        "odd": [5.0, 1.0, 3.0],
        "even": [4.0, 1.0, 3.0, 2.0],
        "inf_pair": [np.inf, -np.inf],
        "inf_same": [np.inf, np.inf],
        "neg_inf_same": [-np.inf, -np.inf],
        "big_opposite": [-1.7976931348623157e308, 1.7976931348623157e308],
        "big_same": [1e308, 1.5e308],
        "zeros_even": [NZ, 0.0],
        "zeros_odd": [0.0, NZ, 0.0],
        "zeros_odd2": [NZ, 0.0, NZ],
        "nan": [1.0, nan, 3.0],
        "one": [7.0],
        "empty": [],
    }


def _median_cases():
    cs = []
    for name, v in _median_vectors().items():
        for cls in ("double", "single"):
            make = Value(lambda v=v, c=cls: np.array(v, dtype=np.float64).astype(c).reshape(1, -1))
            cs.append(
                Case(
                    f"median_{name}_{cls}",
                    matlab="median",
                    args=[In("v")],
                    python=lambda v: S.median(v.ravel()),
                    inputs={"v": make},
                )
            )
            cs.append(
                Case(
                    f"median_{name}_{cls}_omitnan",
                    matlab="median",
                    args=[In("v"), "omitnan"],
                    python=lambda v, flag: S.median(v.ravel(), omitnan=flag == "omitnan"),
                    inputs={"v": make},
                )
            )
            cs.append(
                Case(
                    f"prctile_{name}_{cls}",
                    matlab="prctile",
                    args=[In("v"), [0, 1, 12.5, 25, 33.3, 50, 66.7, 75, 99, 100]],
                    python=lambda v, p: S.prctile(v.ravel(), p),
                    inputs={"v": make},
                )
            )
    seed = 800
    for cls in ("double", "single"):
        for n in (2, 5, 6, 21, 22, 999, 1000, 100000):
            seed += 1
            make = Value(lambda s=seed, n=n, c=cls: wide(s, (1, n), c).ravel())
            cs.append(
                Case(
                    f"median_rand_{cls}_{n}",
                    matlab="median",
                    args=[In("x")],
                    python=py_median,
                    inputs={"x": make},
                )
            )
            cs.append(
                Case(
                    f"prctile_rand_{cls}_{n}",
                    matlab="prctile",
                    args=[In("x"), [0.5, 2.5, 10, 33.3, 50, 66.7, 90, 97.5, 99, 99.5]],
                    python=py_prctile,
                    inputs={"x": make},
                )
            )
        ties = Value(lambda s=seed, c=cls: ints(s, (1, 41), -6, 6, c).ravel())
        cs.append(
            Case(
                f"median_ties_{cls}",
                matlab="median",
                args=[In("x")],
                python=py_median,
                inputs={"x": ties},
            )
        )
        cs.append(
            Case(
                f"prctile_ties_{cls}",
                matlab="prctile",
                args=[In("x"), [1, 10, 25, 37, 50, 63, 75, 90, 99]],
                python=py_prctile,
                inputs={"x": ties},
            )
        )
        vol = Value(lambda s=seed, c=cls: wide(s + 9, (8, 7, 6), c))
        for d in (1, 2, 3):
            cs.append(
                Case(
                    f"median_{cls}_vol_d{d}",
                    matlab="median",
                    args=[Vol(In("v")), Dim(d)],
                    python=py_median,
                    out_roles=("vol",),
                    inputs={"v": vol},
                )
            )
            cs.append(
                Case(
                    f"prctile_{cls}_vol_d{d}",
                    matlab="prctile",
                    args=[Vol(In("v")), [1, 50, 99.5], Dim(d)],
                    python=py_prctile,
                    out_roles=("vol",),
                    inputs={"v": vol},
                )
            )
        cs.append(
            Case(
                f"median_{cls}_vol_all",
                matlab="median",
                args=[Vol(In("v")), "all"],
                py_args=[Vol(In("v")), None],
                python=py_median,
                inputs={"v": vol},
            )
        )
        cs.append(
            Case(
                f"median_{cls}_tiles_omitnan",
                code="out = {median(in{1}, 4, 'omitnan')};",
                args=[Vol(In("t"))],
                python=lambda t: S.median(t, 0, omitnan=True, keepdims=True),
                out_roles=("vol",),
                inputs={"t": Value(lambda c=cls: _tiles(c))},
            )
        )
    for cls, vals in (
        ("uint16", [1, 2]),
        ("uint16", [65535, 65535, 1, 2]),
        ("int8", [-3, -2]),
        ("int16", [-7, 4]),
        ("int16", [-32768, 32767]),
        ("uint8", [255, 254, 255, 255]),
    ):
        tag = f"{cls}_{'_'.join(str(abs(v)) for v in vals)}"
        make = Value(lambda c=cls, v=vals: np.array(v, dtype=c))
        cs.append(
            Case(
                f"median_{tag}",
                matlab="median",
                args=[In("v")],
                python=py_median,
                inputs={"v": make},
            )
        )
        cs.append(
            Case(
                f"prctile_{tag}",
                matlab="prctile",
                args=[In("v"), [10, 37, 50, 81]],
                python=py_prctile,
                inputs={"v": make},
            )
        )
    # psf_gen_new: median of the positive values of single data; XR_rotate_PSF's prctile
    cs.append(
        Case(
            "median_positive_single_psf",
            code="x = in{1}; out = {median(single(x(x > 0)))};",
            args=[Vol(In("x"))],
            python=lambda x: S.median(_positive(x)),
            inputs={"x": Value(lambda: wide(881, (40, 40, 10), "single"))},
        )
    )
    cs.append(
        Case(
            "prctile_high_positive_single",
            code="x = in{1}; out = {prctile(x(x > 0), 99.99)};",
            args=[Vol(In("x"))],
            python=lambda x: S.prctile(_positive(x), 99.99),
            inputs={"x": Value(lambda: wide(882, (40, 40, 10), "single"))},
        )
    )
    # stitching: median(imSizes) of a list of tile sizes, a single tile's size giving a scalar
    for name, sizes in (
        ("tiles", [[512, 1024, 201], [512, 1024, 199], [512, 1024, 201], [510, 1024, 200]]),
        ("one_tile", [[512, 1024, 201]]),
    ):
        cs.append(
            Case(
                f"median_imsizes_{name}",
                matlab="median",
                args=[In("s")],
                python=lambda s: S.median(
                    s, R.default_axis(s, matlab_layout=True), matlab_layout=True, keepdims=True
                ),
                inputs={"s": Value(lambda z=sizes: np.array(z, dtype=np.float64))},
            )
        )
    return cs


def _positive(x):
    """x(x > 0) of a (z, y, x) array, in MATLAB's order"""
    from petakit5d.core.arrays import matlab_order

    v = matlab_order(x).reshape(-1)
    return v[v > 0]


def _tiles(cls):
    t = wide(891, (6, 5, 4, 4), cls)
    b = _bits64(892, t.size).reshape(t.shape, order="F")
    t[(b % np.uint64(4)) == 0] = dtypes.NAN
    t[..., 2] = dtypes.NAN
    return t


def _mode_cases():
    cs = []
    nan = dtypes.NAN
    vecs = {
        "ties": [3.0, 1.0, 3.0, 1.0, 2.0],
        "nan": [nan, 2.0, nan, 2.0],
        "all_nan": [nan, nan],
        "signed_zeros": [0.0, NZ],
        "signed_zeros2": [NZ, 0.0],
        "unique": [5.0, 6.0, 7.0],
        "empty": [],
    }
    for name, v in vecs.items():
        make = Value(lambda v=v: np.array(v, dtype=np.float64).reshape(1, -1))
        cs.append(
            Case(
                f"mode_{name}",
                code="[m, f] = mode(in{1}); out = {mode(in{1}), m, f};",
                args=[In("v")],
                nargout=3,
                python=lambda v: (S.mode(v.ravel()), *S.mode(v.ravel(), with_count=True)),
                inputs={"v": make},
            )
        )
        cs.append(
            Case(
                f"mode3_{name}",
                code="[m, f] = mode(in{1}); [m3, f3, c3] = mode(in{1}); out = {m3, f3};",
                args=[In("v")],
                nargout=2,
                python=lambda v: S.mode(v.ravel(), with_count=True, first_of_run=True),
                inputs={"v": make},
            )
        )
    sizes = Value(
        lambda: np.array(
            [[512, 1024, 201], [512, 1024, 199], [512, 1024, 201], [510, 1024, 199]],
            dtype=np.float64,
        )
    )
    cs.append(
        Case(
            "mode_imsizes_d1",
            matlab="mode",
            args=[In("s"), 1],
            python=lambda s, d: S.mode(s, 0, keepdims=True, matlab_layout=True),
            inputs={"s": sizes},
        )
    )
    cs.append(
        Case(
            "mode_groups",
            matlab="mode",
            args=[In("g")],
            python=lambda g: S.mode(g),
            inputs={"g": Value(lambda: ints(901, (1, 60), 1, 5, "double").ravel())},
        )
    )
    return cs


CASES = _var_cases() + _median_cases() + _mode_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- unit tests --------------------------------------------------------------------------------


def test_var_uses_abs_so_nan_loses_its_sign():
    v = S.var(np.array([2.0, dtypes.NAN]))
    assert np.isnan(v) and not np.signbit(v)


def test_median_meanof_rules():
    assert S.median(np.array([np.inf, np.inf])) == np.inf
    assert np.isnan(S.median(np.array([-np.inf, np.inf])))
    assert S.median(np.array([-1.7976931348623157e308, 1.7976931348623157e308])) == 0
    assert S.median(np.array([1e308, 1.5e308])) == 1.25e308
    assert not np.signbit(S.median(np.array([NZ, NZ])))  # -0 + (-0 - -0)/2 is +0
    assert S.median(np.array([-3, -2], np.int8)) == -3


def test_median_odd_takes_the_stable_middle():
    assert np.signbit(S.median(np.array([0.0, NZ, 0.0])))
    assert not np.signbit(S.median(np.array([NZ, 0.0, NZ])))
    # whole-array selection (partition) and the sorted path agree
    x = ints(5, (9, 7, 11), -3, 3, "double")
    x[x == 0] = np.where(_bits64(6, int(np.sum(x == 0))) % np.uint64(2) == 0, NZ, 0.0)
    a = S.median(x)
    b = S.median(R._fortran_flat(R._mview(x, False), x.dtype), 0)
    assert np.asarray(a).tobytes() == np.asarray(b).tobytes()


@pytest.mark.parametrize("n", [1, 2, 3, 10, 101])
def test_prctile_endpoints(n):
    x = wide(n, (n,), "double")
    s = np.sort(x)
    p = S.prctile(x, [0, 100])
    assert p[0] == s[0] and p[1] == s[-1]


def test_prctile_axis_equals_vector():
    v = wide(3, (8, 7, 6), "single")
    p = [1, 37.5, 50, 99.5]
    r = S.prctile(v, p, 0)
    for i in range(v.shape[1]):
        for j in range(v.shape[2]):
            assert r[:, i, j].tobytes() == S.prctile(v[:, i, j].copy(), p).tobytes()


def test_mode_returns_last_or_first_of_run():
    assert np.signbit(S.mode(np.array([0.0, NZ])))
    assert not np.signbit(S.mode(np.array([0.0, NZ]), first_of_run=True))
    m, f = S.mode(np.array([3.0, 1, 3, 1, 2]), with_count=True)
    assert m == 1 and f == 2


def test_errors():
    with pytest.raises(R.MatlabError):
        S.var(np.array([1, 2], np.int32))
    with pytest.raises(ValueError):
        S.prctile(np.array([1.0]), 101)
    with pytest.raises(TypeError):
        S.median(np.zeros((2, 2)), (0, 1))
