"""Tests of petakit5d.detection.stats: MATLAB's norminv, tcdf, tinv, chi2inv, fcdf, ksdensity,
ecdf, mad and kstest2 against MATLAB R2024b (oracle cases), on dense grids with their tails, edge
values, NaN and Inf, and in the forms the point-detection code calls them (task C01's traces).

Grids are made by exact arithmetic or from decimal literals, so they are the same bytes on every
platform. The large grids are compared with the golden store; each function also has a small
case (under 2 KB) that the manifest keeps whole, so CI compares it too. Tolerances: the rows of
docs/detection_scope.md, statistics (proposed). No MATLAB test covers these functions.
"""

import math

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.detection import stats
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.compare import tol

NAN = dtypes.NAN
INF = math.inf
EDGES = [0.0, -0.0, 1.0, -0.5, 1.5, -INF, INF, NAN]


def _f(*v) -> np.ndarray:
    return np.array(v, dtype=np.float64)


def _decades(lo: int, hi: int, mantissas=(1, 2, 5)) -> list[float]:
    """m * 10^-k for k in [lo, hi], from decimal literals"""
    return [float(f"{m}e-{k}") for k in range(hi, lo - 1, -1) for m in mantissas]


def p_core() -> np.ndarray:
    """probabilities in [1e-6, 1 - 1e-6]: 1e-4 steps, decades, around 1/2, the levels the
    detection code uses, and the edges"""
    mid = [k / 10000 for k in range(1, 10000)]
    low = _decades(4, 6)
    high = [1 - v for v in low]
    half = [0.5 + k * 2.0**-53 for k in range(-4, 5)]
    levels = [0.05, 0.8, 0.9, 0.95, 0.975, 0.025, 0.2, 0.1, 0.001, 0.999]
    return _f(*(mid + low + high + half + levels + EDGES))


def p_tail() -> np.ndarray:
    """the far tails: realmin to 1e-7, and 1 - 1e-7 to 1 - 2^-53"""
    low = [2.2250738585072014e-308] + _decades(7, 300, (1,))[::7] + _decades(7, 20)
    high = [1 - float(f"1e-{k}") for k in range(7, 17)] + [1 - 2.0**-52, 1 - 2.0**-53]
    return _f(*(low + high))


def p_subnormal() -> np.ndarray:
    return _f(5e-324, 1e-320, 1e-315, 1e-310)


def t_core() -> np.ndarray:
    """t values: 1/100 steps in [-40, 40], the edges"""
    return _f(*([k / 100 for k in range(-4000, 4001)] + [-INF, INF, NAN, -0.0]))


def t_tail() -> np.ndarray:
    big = [50.0, 75.75, 100.0, 125.5, 1e3, 1e4, 1e6, 1e10, 1e100, 1e300]
    return _f(*(big + [-b for b in big]))


def f_grid() -> np.ndarray:
    """F values: 1/200 steps in [0, 20], tails to 1e300, Inf, below 0, NaN"""
    dense = [k / 200 for k in range(4001)]
    big = [25.5, 50.0, 100.0, 1e3, 1e6, 1e10, 1e100, 1e300, INF]
    return _f(*(dense + big + [-0.0, -1.0, -INF, NAN, 5e-324, 1e-300, 1e-10]))


def small(grid: np.ndarray, n: int = 200) -> np.ndarray:
    """an evenly spread subset of a grid, edges kept (for the cases CI compares)"""
    fin = grid[np.isfinite(grid)]
    idx = np.unique(np.linspace(0, fin.size - 1, n - 8).astype(np.int64))
    return np.concatenate([fin[idx], grid[~np.isfinite(grid)]])[:n]


T_NU = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7.25, 10, 30, 100, 343.5, 555.75, 593.125, 1000, 1e4]
T_NU_LARGE = [1e6, 1e7, 1.5e7]  # MATLAB's own large-nu form, up to 1.3e-7 off
T_NU_NORMAL = [3e7, 1e8, 1e10, INF]  # the normal cdf
T_NU_BAD = [0, -1, NAN]
TINV_NU = [0.5, 1, 1.5, 2, 3, 4, 5, 7.25, 10, 30, 100, 343.5, 555.75, 593.125, 1000, 1e4]
TINV_NU += [1e6, INF, 0, -1, NAN]
CHI_NU = [0.5, 1, 2, 3, 4.5, 10, 27, 64.5, 125, 343, 729, 1000, 1331, 1e4, 1e6, 0, -1, NAN, INF]
F_NU = [(1, 1), (1, 10), (2, 5), (3, 1000), (5, 2), (10, 10), (30.5, 7.25), (100, 100)]
F_NU += [(343, 686), (1e4, 1e4), (INF, 10), (10, INF), (0, 5), (-1, 5), (NAN, 5)]
F_NU_LARGE = [(1e6, 1)]  # MATLAB up to 7.9e-9 off


def _tag(v) -> str:
    return repr(float(v)).replace("-", "m").replace(".", "p").replace("+", "")


def welch_t() -> tuple[np.ndarray, np.ndarray]:
    """t statistics and Welch degrees of freedom in the form detection and track processing
    pass to tcdf (vectors of one size, large negative t, non-integer df)"""
    k = np.arange(1, 1201, dtype=np.float64)
    return -6000.0 + k * 5.25, 300.0 + k * 0.4375


# Tolerances, proposed from the measurements of task C01 (docs/detection_scope.md, statistics;
# the rows S1 to S8). MATLAB's special functions are closed; where the differences are large,
# MATLAB is the one off the exact value.
TOL_NORMINV = tol(reason="detection_scope.md statistics S1 (proposed)", ulp=4, rel_max=2e-15)
TOL_NORMINV_TAIL = tol(reason="detection_scope.md statistics S1 (proposed)", rel=1e-12)
TOL_SUBNORMAL = tol(reason="detection_scope.md statistics S2 (proposed)", rel=2e-4)
TOL_T = tol(reason="detection_scope.md statistics S3 (proposed)", rel=1e-11, abs=1e-300)
TOL_T_LARGE = tol(reason="detection_scope.md statistics S4 (proposed)", rel=2e-7, abs=1e-300)
TOL_T_NORMAL = tol(reason="detection_scope.md statistics S3 (proposed)", rel=1e-13, abs=1e-300)
TOL_TINV = tol(reason="detection_scope.md statistics S5 (proposed)", rel=2e-10, abs=1e-14)
TOL_CHI = tol(reason="detection_scope.md statistics S6 (proposed)", rel=1e-11)
TOL_F = tol(reason="detection_scope.md statistics S7 (proposed)", rel=1e-11, abs=1e-300)
TOL_F_LARGE = tol(reason="detection_scope.md statistics S8 (proposed)", rel=1e-8, abs=1e-300)
# ksdensity: the density only (grid and bandwidth exact); not called by the traced code
TOL_KSD = tol(reason="detection_scope.md statistics S9 (proposed)", ulp=16)
TOL_KSD_LARGE = tol(reason="detection_scope.md statistics S10 (proposed)", rel=0.02)


def _norminv_cases():
    cs = []
    for mu, sigma in [(0.0, 1.0), (3.0, 2.0), (-1.5, 0.25), (0.0, 0.0), (0.0, -1.0)]:
        cs.append(
            Case(
                f"norminv.core.mu{_tag(mu)}.sigma{_tag(sigma)}",
                matlab="norminv",
                args=[In("p"), mu, sigma],
                python=stats.norminv,
                inputs={"p": Value(p_core)},
                compare=TOL_NORMINV,
            )
        )
    for mu, sigma in [(NAN, 1.0), (INF, 1.0), (0.0, INF), (0.0, NAN)]:
        cs.append(
            Case(
                f"norminv.edge.mu{_tag(mu)}.sigma{_tag(sigma)}",
                matlab="norminv",
                args=[_f(0.0, 0.25, 0.5, 0.975, 1.0, NAN), mu, sigma],
                python=stats.norminv,
                compare=TOL_NORMINV,
            )
        )
    cs += [
        Case(
            "norminv.one_argument",
            matlab="norminv",
            args=[In("p")],
            python=stats.norminv,
            inputs={"p": Value(p_core)},
            compare=TOL_NORMINV,
        ),
        Case(
            "norminv.tail",
            matlab="norminv",
            args=[In("p")],
            python=stats.norminv,
            inputs={"p": Value(p_tail)},
            compare=TOL_NORMINV_TAIL,
        ),
        Case(
            "norminv.subnormal",
            matlab="norminv",
            args=[p_subnormal()],
            python=stats.norminv,
            compare=TOL_SUBNORMAL,
        ),
        Case(
            "norminv.detection_level",
            matlab="norminv",
            args=[0.975, 0, 1],
            python=stats.norminv,
            compare=TOL_NORMINV,
        ),
        Case(
            "norminv.small",
            matlab="norminv",
            args=[small(p_core()), 0.5, 1.5],
            python=stats.norminv,
            compare=TOL_NORMINV,
        ),
    ]
    return cs


def _t_cases():
    cs = []
    groups = [
        (T_NU + T_NU_BAD, TOL_T),
        (T_NU_LARGE, TOL_T_LARGE),
        (T_NU_NORMAL, TOL_T_NORMAL),
    ]
    for nus, spec in groups:
        for nu in nus:
            cs.append(
                Case(
                    f"tcdf.core.nu{_tag(nu)}",
                    matlab="tcdf",
                    args=[In("x"), nu],
                    python=stats.tcdf,
                    inputs={"x": Value(t_core)},
                    compare=spec,
                )
            )
            cs.append(
                Case(
                    f"tcdf.tail.nu{_tag(nu)}",
                    matlab="tcdf",
                    args=[t_tail(), nu],
                    python=stats.tcdf,
                    compare=spec,
                )
            )
    for nu in (1, 3, 593.125, INF):
        cs.append(
            Case(
                f"tcdf.upper.nu{_tag(nu)}",
                matlab="tcdf",
                args=[In("x"), nu, "upper"],
                python=lambda x, nu, _u: stats.tcdf(x, nu, upper=True),
                inputs={"x": Value(t_core)},
                compare=TOL_T_NORMAL if math.isinf(nu) else TOL_T,
            )
        )
    cs += [
        Case(
            "tcdf.welch_vectors",
            matlab="tcdf",
            args=[In("x"), In("nu")],
            python=stats.tcdf,
            inputs={"x": Value(lambda: welch_t()[0]), "nu": Value(lambda: welch_t()[1])},
            compare=TOL_T,
        ),
        Case(
            "tcdf.small",
            matlab="tcdf",
            args=[small(t_core()), 555.75],
            python=stats.tcdf,
            compare=TOL_T,
        ),
        Case(
            "tcdf.small_cauchy",
            matlab="tcdf",
            args=[small(t_core()), 1],
            python=stats.tcdf,
        ),
    ]
    for nu in TINV_NU:
        cs.append(
            Case(
                f"tinv.core.nu{_tag(nu)}",
                matlab="tinv",
                args=[In("p"), nu],
                python=stats.tinv,
                inputs={"p": Value(p_core)},
                compare=TOL_NORMINV if math.isinf(nu) else TOL_TINV,
            )
        )
    cs += [
        Case(
            "tinv.detection_levels",
            matlab="tinv",
            args=[_f(0.05, 0.8, 0.95, 0.2), _f(593.125, 365.25, 812.5, 299.75)],
            python=stats.tinv,
            compare=TOL_TINV,
        ),
        Case(
            "tinv.small",
            matlab="tinv",
            args=[small(p_core()), 343.5],
            python=stats.tinv,
            compare=TOL_TINV,
        ),
        Case(
            "tinv.small_cauchy",
            matlab="tinv",
            args=[small(p_core()), 1],
            python=stats.tinv,
        ),
    ]
    return cs


def _chi_f_cases():
    cs = []
    for nu in CHI_NU:
        cs.append(
            Case(
                f"chi2inv.core.nu{_tag(nu)}",
                matlab="chi2inv",
                args=[In("p"), nu],
                python=stats.chi2inv,
                inputs={"p": Value(p_core)},
                compare=TOL_CHI,
            )
        )
    cs += [
        Case(
            "chi2inv.window_sizes",
            matlab="chi2inv",
            args=[0.9, _f(27, 125, 147, 245, 343, 405, 567, 729, 1331)],
            python=stats.chi2inv,
            compare=TOL_CHI,
        ),
        Case(
            "chi2inv.small",
            matlab="chi2inv",
            args=[small(p_core()), 343],
            python=stats.chi2inv,
            compare=TOL_CHI,
        ),
    ]
    for (nu1, nu2), spec in [(v, TOL_F) for v in F_NU] + [(v, TOL_F_LARGE) for v in F_NU_LARGE]:
        cs.append(
            Case(
                f"fcdf.grid.nu{_tag(nu1)}.{_tag(nu2)}",
                matlab="fcdf",
                args=[In("x"), nu1, nu2],
                python=stats.fcdf,
                inputs={"x": Value(f_grid)},
                compare=spec,
            )
        )
    for nu1, nu2 in [(2, 5), (343, 686)]:
        cs.append(
            Case(
                f"fcdf.upper.nu{_tag(nu1)}.{_tag(nu2)}",
                matlab="fcdf",
                args=[In("x"), nu1, nu2, "upper"],
                python=lambda x, a, b, _u: stats.fcdf(x, a, b, upper=True),
                inputs={"x": Value(f_grid)},
                compare=TOL_F,
            )
        )
    cs.append(
        Case(
            "fcdf.small",
            matlab="fcdf",
            args=[small(f_grid()), 30.5, 7.25],
            python=stats.fcdf,
            compare=TOL_F,
        )
    )
    return cs


def sample(n: int, stream: str, mu: float = 10.0, sigma: float = 2.0) -> np.ndarray:
    """n normal values (bit-stable), as a MATLAB row"""
    from petakit5d.testing import synth

    return mu + sigma * synth.CounterRNG(2026, f"stats/{stream}").normal(n)


def sample_ties(n: int, stream: str) -> np.ndarray:
    """values on a 0.1 grid (ties), with two NaN"""
    v = np.floor(sample(n, stream) * 10.0 + 0.5) / 10.0
    v[[3, n // 2]] = NAN
    return v


def _sample_cases():
    cs = []
    for name, make in [
        ("n100", lambda: sample(100, "a")),
        ("n37_ties_nan", lambda: sample_ties(37, "b")),
        ("n1000", lambda: sample(1000, "c", 0.0, 1.0)),
    ]:
        for flag in (0, 1):
            cs.append(
                Case(
                    f"mad.flag{flag}.{name}",
                    matlab="mad",
                    args=[In("x"), flag],
                    python=lambda x, f: stats.mad(x, f),
                    inputs={"x": Value(make)},
                )
            )
    cs.append(
        Case(
            "mad.columns",
            matlab="mad",
            args=[In("x"), 1],
            python=lambda x, f: stats.mad(x, f, matlab_layout=True),
            inputs={"x": Value(lambda: sample(60, "d").reshape(4, 15, order="F"))},
        )
    )
    for name, make in [
        ("n60_ties_nan", lambda: sample_ties(60, "e")),
        ("n500", lambda: sample(500, "f")),
    ]:
        cs.append(
            Case(
                f"ecdf.{name}",
                matlab="ecdf",
                args=[In("y")],
                nargout=2,
                python=stats.ecdf,
                inputs={"y": Value(make)},
            )
        )
    for name, make, spec in [
        ("n100", lambda: sample(100, "g"), TOL_KSD),
        ("n500", lambda: sample(500, "h", 3.0, 0.5), TOL_KSD_LARGE),
    ]:
        cs.append(
            Case(
                f"ksdensity.default.{name}",
                matlab="ksdensity",
                args=[In("x")],
                nargout=3,
                python=stats.ksdensity,
                inputs={"x": Value(make)},
                compare={0: spec},
            )
        )
    cs.append(
        Case(
            "ksdensity.at_points",
            matlab="ksdensity",
            args=[In("x"), _f(*[k / 4 for k in range(16, 57)])],
            nargout=3,
            python=stats.ksdensity,
            inputs={"x": Value(lambda: sample(200, "i"))},
            compare={0: TOL_KSD},
        )
    )
    for tail in ("unequal", "larger", "smaller"):
        cs.append(
            Case(
                f"kstest2.{tail}",
                matlab="kstest2",
                args=[In("a"), In("b"), "Tail", tail],
                nargout=3,
                python=lambda a, b, _t, t: stats.kstest2(a, b, tail=t),
                inputs={
                    "a": Value(lambda: sample(50, "j")),
                    "b": Value(lambda: sample(80, "k", 10.5)),
                },
            )
        )
    cs.append(
        Case(
            "kstest2.ties_nan_alpha",
            matlab="kstest2",
            args=[In("a"), In("b"), "Alpha", 0.01],
            nargout=3,
            python=lambda a, b, _a, al: stats.kstest2(a, b, alpha=al),
            inputs={
                "a": Value(lambda: sample_ties(40, "l")),
                "b": Value(lambda: sample_ties(70, "m")),
            },
        )
    )
    return cs


CASES = _norminv_cases() + _t_cases() + _chi_f_cases() + _sample_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- rules (no MATLAB) --------------------------------------------------------------------------


def test_scalar_arguments_give_a_scalar_and_single_gives_single():
    assert np.ndim(stats.norminv(0.5)) == 0
    assert stats.tinv(np.float32(0.9), 10).dtype == np.float32
    assert stats.tcdf(np.zeros((2, 3)), 5).shape == (2, 3)


def test_sizes_must_match_unless_scalar():
    with pytest.raises(ValueError) as e:
        stats.tcdf(np.zeros(3), np.ones(4))
    assert e.value.identifier == "stats:tcdf:InputSizeMismatch"


def test_edges_of_the_quantile_functions():
    assert stats.norminv(0.0) == -INF and stats.norminv(1.0) == INF
    assert np.isnan(stats.norminv(0.5, 0, 0))
    assert stats.tinv(0.0, 7) == -INF and stats.tinv(1.0, 7) == INF
    assert np.isnan(stats.tinv(0.5, 0)) and np.isnan(stats.tinv(1.5, 3))
    assert stats.chi2inv(0.0, 4) == 0.0 and stats.chi2inv(1.0, 4) == INF
    assert stats.chi2inv(0.5, 0) == 0.0 and np.isnan(stats.chi2inv(1.0, 0))
    assert np.isnan(stats.chi2inv(0.5, INF))


def test_the_exact_cauchy_forms_are_symmetric():
    x = np.array([0.5, 2.0, 11.66, 1e10])
    assert np.array_equal(stats.tcdf(-x, 1), stats.tcdf(x, 1, upper=True))
    assert stats.tcdf(-1e300, 1) > 0  # atan2(1, 1e300)/pi, not 0.5 + atan(x)/pi
