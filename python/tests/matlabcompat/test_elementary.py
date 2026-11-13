"""Tests of petakit5d.matlabcompat.elementary against MATLAB R2024b: sin, cos, tan, cot, exp,
sqrt, hypot, atan2 and power of doubles on dense grids and random values, sind and cosd,
PetaKit5D's fastPower, tukeywin, padarray and circshift.

Inputs are made with IEEE operations only (k/64, k*pi/1024, integers times powers of two), so
they are the same on every platform; the results are compared bit for bit, sind and cosd (S03's
rule) on a dense grid, S03's 1,979 angles and 27,000 random angles too. The Hann window of
tukeywin (r >= 1) differs in the last bit for some n above 11: only n up to 11 are cases. No
MATLAB test covers these functions alone.
"""

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.matlabcompat import elementary as E
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, PerAxis, Value, Vol

_C1 = np.uint64(0xBF58476D1CE4E5B9)
_C2 = np.uint64(0x94D049BB133111EB)
_GOLD = np.uint64(0x9E3779B97F4A7C15)
NZ = -0.0


def _bits64(seed: int, n: int) -> np.ndarray:
    with np.errstate(over="ignore"):
        z = np.uint64(seed) + np.arange(1, n + 1, dtype=np.uint64) * _GOLD
        z = (z ^ (z >> np.uint64(30))) * _C1
        z = (z ^ (z >> np.uint64(27))) * _C2
        return z ^ (z >> np.uint64(31))


def rand53(seed: int, n: int, scale_exp: int) -> np.ndarray:
    """53-bit integers of both signs times 2^scale_exp"""
    b = _bits64(seed, n)
    m = (b >> np.uint64(11)).astype(np.float64)
    return np.ldexp(np.where(b & np.uint64(1), -m, m), scale_exp)


def dense() -> np.ndarray:
    k = np.arange(-4096, 4097, dtype=np.float64)
    j = np.arange(-20000, 20001, dtype=np.float64)
    special = [
        0.0,
        NZ,
        np.inf,
        -np.inf,
        dtypes.NAN,
        2.2250738585072014e-308,
        1e-300,
        5e-324,
        np.pi / 2,
        np.pi,
        2 * np.pi,
        1e5,
        -1e5,
    ]
    return np.concatenate(
        [
            k / 64,
            j * np.pi / 1024,
            rand53(11, 20000, -46),
            rand53(12, 20000, -60),
            rand53(13, 5000, -40),
            special,
        ]
    )


def _fn_cases():
    cs = []
    x = Value(dense)
    for fn in ("sin", "cos", "tan", "cot", "exp", "sqrt"):
        args = [In("x")]
        code = None
        if fn == "sqrt":
            code = "out = {sqrt(abs(in{1}))};"
        if fn == "exp":
            code = "x = in{1}; out = {exp(x(abs(x) < 800))};"
        cs.append(
            Case(
                f"{fn}_double_dense",
                matlab=None if code else fn,
                code=code,
                args=args,
                python={
                    "sqrt": lambda v: E.sqrt(np.abs(v)),
                    "exp": lambda v: E.exp(v[np.abs(v) < 800]),
                }.get(fn, getattr(E, fn)),
                inputs={"x": x},
            )
        )
    cs.append(
        Case(
            "exp_double_grid",
            matlab="exp",
            args=[In("x")],
            python=E.exp,
            inputs={
                "x": Value(
                    lambda: np.concatenate(
                        [np.arange(-8000, 8001) / 16.0, [-745.2, 709.78, -708.4, 1.0, 0.5]]
                    )
                )
            },
        )
    )
    for fn in ("hypot", "atan2"):
        cs.append(
            Case(
                f"{fn}_double_dense",
                code=f"x = in{{1}}; out = {{{fn}(x, fliplr(x))}};",
                args=[In("x")],
                python=lambda v, f=E.hypot if fn == "hypot" else E.atan2: f(v, v[::-1]),
                inputs={"x": x},
            )
        )
    return cs


_ANGLES = [
    -720.0,
    -540,
    -450,
    -360,
    -315,
    -270,
    -225,
    -180,
    -135,
    -90,
    -45,
    -30,
    0,
    30,
    45,
    60,
    90,
    120,
    135,
    150,
    180,
    210,
    270,
    300,
    360,
    405,
    450,
    540,
    720,
    1e15,
    1e16,
    NZ,
    0.5,
    89.99999999999999,
    90.00000000000001,
    179.99999999999997,
    1e300,
    np.inf,
    dtypes.NAN,
    45.000000000000007,
    32.45,
    31.5,
    30.1,
    -32.45,
    32.4,
    31.8,
    30.0,
    33.0,
    45.0,
]


def s03_angles() -> np.ndarray:
    """the 1,979 angles of S03's elem.trig case (stitching geometry)"""
    a = [k / 100 for k in range(-36000, 36001, 137)] + [k / 4 for k in range(-720, 721)]
    a += [32.45, 31.5, 30.0, 45.0, 60.0, 90.0, -30.0, 32.8, 0.0, 359.95, 0.01, 1e-9]
    return np.array(a, dtype=np.float64)


def wide_angles() -> np.ndarray:
    """random angles of every size, multiples of 7.5 degrees and their neighbours, +-45 + k ulp,
    tiny and signed-zero angles"""
    k = np.arange(-96, 97, dtype=np.float64) * 7.5
    near = np.concatenate([np.nextafter(k, np.inf), np.nextafter(k, -np.inf)])
    q = np.array([45.0, 135.0, 225.0, 315.0, -45.0, -135.0])
    ulps = np.concatenate([q + i * np.spacing(q) for i in range(-3, 4)])
    tiny = [
        NZ,
        0.0,
        5e-324,
        -5e-324,
        1e-300,
        -1e-300,
        2.2250738585072014e-308,
        359.99999999999994,
        -720.0000000000001,
        1e9 + 0.5,
        -3e12 - 0.25,
    ]
    return np.concatenate(
        [
            rand53(31, 20000, -40),
            rand53(32, 5000, -26),
            rand53(33, 2000, -50),
            k,
            near,
            ulps,
            np.array(tiny),
        ]
    )


def _deg_cases():
    cs = []
    # angles above 2^56 where x - fix(x/360)*360 and the exact remainder differ
    huge = [
        8.78662275217687e16,
        1.9541282940265766e17,
        4.999914160116846e17,
        1.0039373046294822e17,
        1.4845544263400305e18,
        9.053390975729444e17,
    ]
    for fn in ("sind", "cosd"):
        cs.append(
            Case(
                f"{fn}_huge_reduction",
                matlab=fn,
                args=[In("a")],
                python=getattr(E, fn),
                inputs={"a": Value(lambda h=huge: np.array(h))},
            )
        )
    for fn in ("sind", "cosd"):
        # (cosd's set lacks 30.1, as when cosd was 1 ulp off there: the case's key stays)
        angles = [v for v in _ANGLES if fn == "sind" or v != 30.1]
        cs.append(
            Case(
                f"{fn}_special",
                matlab=fn,
                args=[In("a")],
                python=getattr(E, fn),
                inputs={"a": Value(lambda a=angles: np.array(a))},
            )
        )
        cs.append(
            Case(
                f"{fn}_dense_grid",
                matlab=fn,
                args=[In("a")],
                python=getattr(E, fn),
                inputs={"a": Value(lambda: np.arange(-3600, 3601) / 10.0)},
            )
        )
        for name, make in (("s03_angles", s03_angles), ("wide_angles", wide_angles)):
            cs.append(
                Case(
                    f"{fn}_{name}",
                    matlab=fn,
                    args=[In("a")],
                    python=getattr(E, fn),
                    inputs={"a": Value(make)},
                )
            )
    return cs


def _power_cases():
    cs = []
    xp = Value(lambda: np.abs(rand53(21, 3000, -50)))
    for p in (-15.0, -2.5, -1.0, 0.0, 0.1, 1 / 3, 0.5, 1.0, 2.0, 2.5, 3.0, 7.0, 10.0, 15.0):
        cs.append(
            Case(
                f"power_double_{str(p).replace('.', 'p').replace('-', 'm')}",
                code="out = {in{1} .^ in{2}};",
                args=[In("x"), p],
                python=E.power,
                inputs={"x": xp},
            )
        )
    # distance_weight_single_axis: dfactor .^ (1 : n)
    cs.append(
        Case(
            "power_dfactor_series",
            code="out = {0.99 .^ (1 : 1500)};",
            args=[],
            python=lambda: E.power(0.99, np.arange(1, 1501, dtype=np.float64)),
        )
    )
    cs.append(
        Case(
            "power_root_degree",
            code="out = {(1e-50) .^ (1 ./ (1 : 20)), 10 .^ (1/10)};",
            args=[],
            nargout=2,
            python=lambda: (
                E.power(1e-50, 1 / np.arange(1, 21, dtype=np.float64)),
                E.power(10.0, 1 / 10),
            ),
        )
    )
    base = Value(
        lambda: ((np.arange(1, 4001, dtype=np.float32) + np.float32(1)) / np.float32(10)).astype(
            np.float32
        )
    )
    for w in (0, 1, 2, 3, 5, 10, 13, -1, -2):
        cs.append(
            Case(
                f"fastpower_single_w{w}".replace("-", "m"),
                matlab="fastPower",
                args=[In("b"), float(w)],
                python=E.fast_power,
                inputs={"b": base},
            )
        )
    cs.append(
        Case(
            "fastpower_tukey_window",
            code="w = tukeywin(110, 0.5); out = {fastPower(w(5 : 104), 10)};",
            args=[],
            python=lambda: E.fast_power(E.tukeywin(110, 0.5)[4:104], 10),
        )
    )
    return cs


def _window_cases():
    cs = []
    for r in (0.5, 0.05, 0.1, 0.25, 1 / 3, 0.75, 0.9, 0.0, -0.2):
        for n in (1, 2, 3, 10, 11, 44, 99, 100, 277):
            cs.append(
                Case(
                    f"tukeywin_{n}_{str(round(r, 4)).replace('.', 'p').replace('-', 'm')}",
                    matlab="tukeywin",
                    args=[float(n), r],
                    python=lambda n, r: E.tukeywin(int(n), r),
                )
            )
    for n in range(1, 12):
        cs.append(
            Case(
                f"tukeywin_hann_{n}",
                matlab="tukeywin",
                args=[float(n), 1.0],
                python=lambda n, r: E.tukeywin(int(n), r),
            )
        )
    # compute_tile_bwdist's window along z
    for nz in (5, 67, 201, 500):
        cs.append(
            Case(
                f"tukeywin_tile_window_{nz}",
                code=(
                    "sz3 = in{1}; nw = max(round(sz3 * 1.1), sz3 + 2); w = tukeywin(nw, 0.5);"
                    " ws = max(round(sz3 * 0.05), 2); out = {w(ws : ws + sz3 - 1)};"
                ),
                args=[float(nz)],
                python=lambda s: _tile_window(int(s)),
            )
        )
    return cs


def _tile_window(sz3: int) -> np.ndarray:
    from petakit5d.core import dtypes as D

    nw = max(int(D.mround(sz3 * 1.1)), sz3 + 2)
    w = E.tukeywin(nw, 0.5)
    ws = max(int(D.mround(sz3 * 0.05)), 2)
    return w[ws - 1 : ws - 1 + sz3]


def _pad_cases():
    cs = []
    a = Value(lambda: np.arange(1, 61, dtype=np.float64).reshape(3, 4, 5, order="F"))
    u = Value(lambda: np.arange(1, 61, dtype=np.uint16).reshape(3, 4, 5, order="F"))
    specs = [
        ("const_post", [1, 2, 1], [0, "post"], {"value": 0, "direction": "post"}),
        ("const_both_7", [1, 1, 2], [7], {"value": 7}),
        (
            "replicate_pre",
            [1, 2, 3],
            ["replicate", "pre"],
            {"method": "replicate", "direction": "pre"},
        ),
        ("symmetric_both", [4, 5, 6], ["symmetric", "both"], {"method": "symmetric"}),
        ("circular_both", [3, 9, 11], ["circular", "both"], {"method": "circular"}),
        (
            "symmetric_post",
            [5, 1, 0],
            ["symmetric", "post"],
            {"method": "symmetric", "direction": "post"},
        ),
    ]
    for name, pad, margs, kw in specs:
        cs.append(
            Case(
                f"padarray_{name}",
                matlab="padarray",
                args=[Vol(In("a")), PerAxis(pad, "yxz"), *margs],
                py_args=[Vol(In("a")), PerAxis(pad, "yxz")],
                py_kwargs=kw,
                python=E.padarray,
                out_roles=("vol",),
                inputs={"a": a},
            )
        )
    cs.append(
        Case(
            "padarray_uint16_value_rounded",
            matlab="padarray",
            args=[Vol(In("u")), PerAxis([1, 1, 1], "yxz"), 3.6, "both"],
            py_args=[Vol(In("u")), PerAxis([1, 1, 1], "yxz")],
            py_kwargs={"value": 3.6},
            python=E.padarray,
            out_roles=("vol",),
            inputs={"u": u},
        )
    )
    for name, k in (("mixed", [1, -1, 2]), ("big", [-7, 9, -13])):
        cs.append(
            Case(
                f"circshift_{name}",
                matlab="circshift",
                args=[Vol(In("a")), PerAxis(k, "yxz")],
                python=lambda v, k: E.circshift(v, tuple(int(s) for s in k), (0, 1, 2)),
                out_roles=("vol",),
                inputs={"a": a},
            )
        )
    return cs


CASES = _fn_cases() + _deg_cases() + _power_cases() + _window_cases() + _pad_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- unit tests --------------------------------------------------------------------------------


def test_results_do_not_use_the_c_library():
    """the fdlibm ports give fixed bits (values MATLAB returns; glibc's differ for these)"""
    assert E.exp(1.0) == 2.718281828459045
    assert E.sind(30.0) == 0.5 and E.cosd(60.0) == 0.5
    assert E.sind(180.0) == 0 and not np.signbit(E.sind(180.0))
    assert np.signbit(E.sind(-180.0))
    assert E.cosd(90.0) == 0 and not np.signbit(E.cosd(-90.0))
    assert E.power(2.0, 0.5) == np.sqrt(2.0)


def test_single_data_rounds_the_double_result():
    x = np.array([0.5, 1.5, 2.5], np.float32)
    assert E.sin(x).dtype == np.float32
    assert E.sin(x).tobytes() == E.sin(x.astype(np.float64)).astype(np.float32).tobytes()


@pytest.mark.parametrize("w", [0, 1, 2, 3, 10, -1, -3])
def test_fast_power_is_binary_powering(w):
    b = np.array([1.1, 0.9, 2.0], np.float32)
    r = E.fast_power(b, w)
    if w == 0:
        assert r == 1.0 and isinstance(r, float)
        return
    e, base = abs(w), (np.float32(1) / b if w < 0 else b)
    ref, sq = None, base.copy()
    while e:
        if e & 1:
            ref = sq.copy() if ref is None else (ref * sq).astype(np.float32)
        sq = (sq * sq).astype(np.float32)
        e >>= 1
    assert np.asarray(r).tobytes() == ref.tobytes()


def test_tukeywin_is_symmetric_where_matlab_is():
    for n in (10, 11, 44, 99):
        w = E.tukeywin(n, 0.5)
        assert w[0] == 0 and w[-1] == 0 and np.all(w <= 1)


def test_padarray_and_circshift():
    a = np.arange(24.0).reshape(2, 3, 4)
    assert E.padarray(a, (1, 0, 2), method="replicate", direction="post").shape == (3, 3, 6)
    assert E.circshift(a, 1, 2)[0, 0, 1] == a[0, 0, 0]
    with pytest.raises(ValueError):
        E.padarray(a, (1, 2), method="mirror")
    with pytest.raises(ValueError):
        E.fast_power(a, 1.5)


def test_sind_and_cosd_rules():
    assert np.signbit(E.sind(NZ)) and np.signbit(E.sind(-180.0)) and not np.signbit(E.sind(180.0))
    assert E.sind(60.0) == E.cosd(30.0) == float.fromhex("0x1.bb67ae8584caap-1")
    assert E.cosd(90.0) == 0 and not np.signbit(E.cosd(-270.0))
    # at +-45 the cosine kernel, for both
    assert E.sind(45.0) == E.cosd(45.0) == E.kernel_cos(np.pi / 4, 0.0)
    a = wide_angles()
    assert E.sind(-a).tobytes() == (-E.sind(a)).tobytes()  # odd
    assert E.cosd(-a).tobytes() == E.cosd(a).tobytes()  # even
    nan_s, nan_c = E.sind(np.inf), E.cosd(np.inf)
    assert np.signbit(nan_s) and not np.signbit(nan_c)


def test_kernels_and_raw_functions():
    x = np.array([0.1, -0.7, 0.785])
    assert E.kernel_sin(x, 0.0).tobytes() == E.fdlibm_sin(x).tobytes()
    assert E.kernel_cos(x, 0.0).tobytes() == E.fdlibm_cos(x).tobytes()
    assert E.kernel_tan(x, 0.0, 1).tobytes() == E.tan(x).tobytes()
    n, y0, y1 = E.rem_pio2(np.array([10.0, -3.0]))
    assert n.tolist() == [6, -2] and abs(y0[0] - (10 - 6 * np.pi / 2)) < 1e-15
    # fdlibm's sine and cosine of Inf are MATLAB's NaN on every CPU; MATLAB's cos gives +NaN
    r = E.fdlibm_cos(np.array([np.inf, -np.inf]))
    assert all(np.signbit(r)) and not any(np.signbit(E.cos(np.array([np.inf, -np.inf]))))
    assert E.kernel_sin(0.3, 0.0).shape == () and E.rem_pio2(5.0)[0].shape == ()
