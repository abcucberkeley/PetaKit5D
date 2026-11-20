"""MATLAB's corr, Pearson, of two columns (matlabcompat.correlation; task C04a): oracle cases and
the forms MATLAB's own arithmetic was probed in. MATLAB computes x'*y / norm(x) / norm(y) on the
centered columns with its BLAS (whose order depends on the CPU and, from 4,096 values, on its
thread count); the port sums in MATLAB's order: within W1 of docs/parity/detection.md."""

import numpy as np
import pytest

from petakit5d.matlabcompat.correlation import CorrError, corr
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, Value, mcell
from petakit5d.testing.compare import tol

W1 = tol(reason="detection.md W1", abs=4e-15)
LOOP = (
    "c = in{1}; r = zeros(numel(c), 1); "
    "for i = 1 : numel(c), a = c{i}; r(i) = corr(a{1}, a{2}); end; out = {r};"
)


def small_pairs() -> list:
    """small integer columns (3 to 9 values), means not always exact"""
    u = synth.CounterRNG(2026, "corr/small").uniform((400, 18))
    out = []
    for k, r in enumerate(u):
        n = 3 + k % 7
        x = np.floor(r[:n] * 41) - 20
        y = np.floor(r[9 : 9 + n] * 41) - 20
        out.append((x.reshape(-1, 1), y.reshape(-1, 1)))
    return out


def formula_pair(n: int) -> tuple[np.ndarray, np.ndarray]:
    k = np.arange(1, n + 1, dtype=np.float64)
    x = np.mod(k * 7919, 10007) / 7 + np.mod(k * 104729, 4099) / 13
    y = 0.5 * x + np.mod(k * 15485863, 30011) / 11 - np.mod(k * 49979687, 1021) / 3
    return x.reshape(-1, 1), y.reshape(-1, 1)


def _loop(pairs):
    return np.array([corr(a, b) for a, b in pairs]).reshape(-1, 1)


def _case(cid, pairs, compare=W1):
    return Case(
        cid,
        code=LOOP,
        args=[In("c")],
        python=_loop,
        py_args=[pairs],
        inputs={"c": Value(lambda: mcell(*[mcell(a, b) for a, b in pairs]))},
        compare=compare,
    )


EDGES = [
    ([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]),  # 1 - eps/2, as MATLAB gives it
    ([1.0, 2.0, 3.0], [3.0, 2.0, 1.0]),
    ([1.0, 2.0], [2.0, 1.0]),
    ([1e300, -1e300, 1.0], [1.0, 2.0, 3.0]),  # no overflow of the norm
    ([1e-200, 2e-200, 3e-200], [3e-200, 1e-200, 2e-200]),  # no underflow
    ([0.1, 0.2, 0.3, 0.4], [0.3, 0.1, 0.4, 0.2]),
    ([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]),  # NaN
    ([1.0, np.nan, 3.0], [1.0, 2.0, 3.0]),  # NaN
    ([5.0], [7.0]),  # NaN
]
CASES = [
    _case("corr.small_integers", small_pairs()),
    _case(
        "corr.formula", [formula_pair(n) for n in (2, 3, 5, 10, 100, 1000, 4095, 4096, 4097, 65536)]
    ),
    _case(
        "corr.edges",
        [(np.reshape(a, (-1, 1)), np.reshape(b, (-1, 1))) for a, b in EDGES],
        compare=None,
    ),  # fmt: skip
    Case(
        "corr.single",
        code="out = {corr(single(in{1}), single(in{2})), corr(single(in{1}), in{2})};",
        args=[np.reshape([1.0, 2.0, 4.0, 7.0], (-1, 1)), np.reshape([1.0, 3.0, 2.0, 5.0], (-1, 1))],
        python=lambda x, y: (
            corr(x.astype(np.float32), y.astype(np.float32)),
            corr(x.astype(np.float32), y),
        ),  # fmt: skip
        nargout=2,
        compare=tol(reason="detection.md W1", abs=2e-7),
    ),
    Case(
        "corr.error.empty",
        matlab="corr",
        args=[np.zeros((0, 1)), np.zeros((0, 1))],
        python=corr,
        expect_error="stats:corr:TooFewInputs",
        py_error=(CorrError,),
    ),
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_probed_values_exactly():
    """the edge values MATLAB gave in the black-box probes, bit for bit"""
    want = [0.9999999999999999, -0.9999999999999999, -0.9999999999999999, -0.49999999999999994,
            0.0, 0.0]  # fmt: skip
    for (a, b), w in zip(EDGES[:6], want, strict=True):
        assert corr(np.array(a), np.array(b)) == w
    from petakit5d.core.dtypes import NAN

    for a, b in EDGES[6:]:  # MATLAB's NaN for 0/0, the data's NaN (here MATLAB's) propagated
        r = corr(np.where(np.isnan(a), NAN, a), np.array(b))
        assert np.float64(r).view(np.uint64) == 0xFFF8000000000000


def test_shapes_and_errors():
    with pytest.raises(NotImplementedError):
        corr(np.ones((4, 2)), np.ones((4, 1)))
    with pytest.raises(CorrError) as e:
        corr(np.ones(3), np.ones(4))
    assert e.value.identifier == "stats:corr:InputSizeMismatch"
