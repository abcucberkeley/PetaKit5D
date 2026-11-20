"""MATLAB's interp1 with the linear method (matlabcompat.interp; task C04a): oracle cases in the
forms track processing and the lifetime analysis call it, and its errors."""

import numpy as np
import pytest

from petakit5d.matlabcompat.interp import InterpError, interp1
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, Value, mcell, mtyped

LOOP = (
    "c = in{1}; out1 = cell(numel(c), 1); "
    "for i = 1 : numel(c), a = c{i}; out1{i} = interp1(a{:}); end; out = {out1};"
)


def _u(stream: str, shape) -> np.ndarray:
    return synth.CounterRNG(2026, f"interp1/{stream}").uniform(shape)


def two_point_calls() -> list:
    """the track processing form: two knots (integer frames), one or two queries between"""
    u = _u("two", (300, 6))
    calls = []
    for r in u:
        x0 = float(np.floor(r[0] * 50))
        x = [x0, x0 + 1 + float(np.floor(r[1] * 3))]
        v = [r[2] * 300, r[3] * 300]
        q = [x[0] + r[4] * (x[1] - x[0])] if r[5] < 0.5 else [x[0] + 1, x[1] - 1 + r[4] * 0.5]
        calls.append([x, v, q])
    return calls


def grid_calls() -> list:
    """uneven grids, queries inside, on knots, outside and NaN, every extrapolation form"""
    calls = []
    for k in range(40):
        n = 2 + k % 23
        u = _u(f"grid{k}", (3, n + 30))
        x = np.cumsum(0.1 + 3 * u[0, :n]) - 7
        v = (u[1, :n] - 0.5) * 100
        q = np.concatenate(
            (
                [x[0] - 1, x[0], x[-1], x[-1] + 0.5],
                x[0] + u[2, : n + 20] * (x[-1] - x[0]),
                [x[n // 2]],
            )
        )
        for extra in ([], ["linear"], ["linear", "extrap"], ["linear", 0.0]):
            calls.append([x.reshape(-1, 1), v.reshape(-1, 1), q.reshape(1, -1), *extra])
        calls.append([x[::-1].copy(), v[::-1].copy(), q])
    return calls


def lifetime_calls() -> list:
    """the lifetime analysis' form: 910 increasing points, 999 queries on [0, 120]"""
    u = _u("lifetime", (2, 910))
    x = (np.sort(u[0]) * 120).reshape(-1, 1)
    v = np.cumsum(u[1]).reshape(-1, 1)
    v = v / v[-1]
    q = (np.arange(999, dtype=np.float64) * 120 / 998).reshape(1, -1)
    return [[x, v, q, "linear"], [x, v, q]]


def knot_calls() -> list:
    """queries on the knots and between them next to Inf, NaN (MATLAB's) and signed zeros: the
    interval a query on a knot takes, and the products with 0"""
    from petakit5d.core.dtypes import NAN

    x = [1.0, 2.0, 3.0]
    q = [1.0, 1.5, 2.0, 2.5, 3.0]
    vals = (
        [1.0, 2.0, np.inf], [np.inf, 2.0, 3.0], [1.0, NAN, 3.0], [NAN, 2.0, 3.0], [1.0, 2.0, NAN],
        [-np.inf, 2.0, np.inf], [np.inf, np.inf, 3.0], [1.0, -np.inf, -np.inf],
        [-3.0, -0.0, 5.0], [3.0, -0.0, -5.0], [-3.0, 0.0, -5.0], [3.0, 0.0, 5.0],
    )  # fmt: skip
    calls = [[x, list(v), q] for v in vals]
    calls += [[x, list(v), [0.5, 3.5], "linear", "extrap"] for v in vals[:2]]
    calls.append([x[::-1], [np.inf, 2.0, 1.0], q])
    return calls


def zero_knot_calls() -> list:
    """signed zeros on the knots when a value is not finite: a term of weight 0 left out or
    added as 0, and only where its value is not finite or wherever the values hold one"""
    from petakit5d.core.dtypes import NAN

    return [
        [[1.0, 2.0, 3.0, 4.0], [np.inf, 3.0, -0.0, 5.0], [3.0, 1.0, 1.5]],
        [[1.0, 2.0, 3.0], [NAN, -0.0, 5.0], [2.0]],
        [[1.0, 2.0, 3.0], [1.0, -0.0, np.inf], [2.0]],
        [[1.0, 2.0, 3.0], [-0.0, NAN, 3.0], [1.0]],
        [[1.0, 2.0], [np.inf, -0.0], [2.0, 1.0]],
        [[1.0, 2.0, 3.0], [-0.0, 5.0, np.inf], [1.0, 2.0]],
        [[1.0, 2.0], [np.inf, -np.inf], [1.5, 1.0, 2.0]],
    ]


def _loop(calls):
    return [interp1(*a) for a in calls]


def _case(cid, calls):
    return Case(
        cid,
        code=LOOP,
        args=[In("c")],
        python=_loop,
        py_args=[calls],
        inputs={"c": Value(lambda: mcell(*[mcell(*a) for a in calls]))},
    )


CASES = [
    _case("interp1.two_points", two_point_calls()),
    _case("interp1.grids", grid_calls()),
    _case("interp1.lifetime", lifetime_calls()),
    _case("interp1.knots_nonfinite", knot_calls()),
    _case("interp1.knots_signed_zeros", zero_knot_calls()),
    Case(
        "interp1.equal_knots",
        code="out = {interp1(in{1}, in{2}, in{3}), interp1(in{1}, in{4}, in{3}), "
        "interp1(in{1}, in{5}, in{3})};",
        args=[[1.0, 2.0], [0.1, 0.1], [1.3, 1.5, 1.71], [-0.0, 0.0], [0.0, -0.0]],
        python=lambda x, a, q, b, c: (interp1(x, a, q), interp1(x, b, q), interp1(x, c, q)),
        nargout=3,
    ),
    Case(
        "interp1.single",
        code="out = {interp1(single(in{1}), single(in{2}), single(in{3})), "
        "interp1(single(in{1}), in{2}, in{3}), interp1(in{1}, in{2}, single(in{3}))};",
        args=[[1.0, 2.0, 3.0], [0.1, 0.3, 0.7], [1.37, 2.21, 2.9]],
        python=lambda x, v, q: (
            interp1(np.float32(x), np.float32(v), np.float32(q)),
            interp1(np.float32(x), np.asarray(v), np.asarray(q)),
            interp1(np.asarray(x), np.asarray(v), np.float32(q)),
        ),
        nargout=3,
    ),
    Case(
        "interp1.shapes",
        code="out = {interp1(in{1}, in{2}, in{3}), interp1(in{1}', in{2}', in{4}), "
        "interp1(in{1}, in{2}, in{5}), interp1(in{1}, in{2}, zeros(1, 0))};",
        args=[
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 4.0],
            np.array([[1.5], [2.5]]),
            [1.5, 2.5],
            np.array([[1.5, 2.5], [2.0, 2.25]]),
        ],  # fmt: skip
        python=lambda x, v, q1, q2, q3: (
            interp1(x, v, q1),
            interp1(np.reshape(x, (-1, 1)), np.reshape(v, (-1, 1)), q2),
            interp1(x, v, q3),
            interp1(x, v, np.zeros((1, 0))),
        ),
        nargout=4,
    ),
]
GI = "MATLAB:griddedInterpolant:"
for cid, args, ident in (
    ("duplicate_points", [[1.0, 1.0, 2.0], [1.0, 2.0, 3.0], 1.5], GI + "NonUniqueCompVecsPtsErrId"),
    ("nan_point", [[1.0, np.nan, 3.0], [1.0, 2.0, 3.0], 2.0], GI + "NonFiniteInputPtsErrId"),
    ("inf_point", [[1.0, 2.0, np.inf], [1.0, 2.0, 3.0], 2.5], GI + "NonFiniteInputPtsErrId"),
    ("one_point", [1.0, 5.0, 1.0], GI + "DegenerateGridErrId"),
    ("integer_grid", [mtyped([1, 2, 3], "int32"), [1.0, 2.0, 4.0], 2.5], "MATLAB:interp1:Xnumeric"),
    ("integer_values", [[1.0, 2.0, 3.0], mtyped([1, 2, 4], "int16"), 2.5],
     "MATLAB:interp1:NonFloatValues"),
):  # fmt: skip
    py = [np.asarray(a.data).reshape(-1) if hasattr(a, "mclass") else a for a in args]
    CASES.append(
        Case(
            f"interp1.error.{cid}",
            matlab="interp1",
            args=args,
            py_args=py,
            python=interp1,
            expect_error=ident,
            py_error=(InterpError,),
        )
    )


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_formula_on_the_traced_values():
    """XR_runTrackProcessing3D's calls recorded by C01's trace (17 digits)"""
    rec = [
        ([25, 27], [141.92105280916775, 144.962873691113], [26], [143.44196325014036]),
        ([45, 48], [302.48543430611107, 303.37470520716096], [46, 47],
         [302.78185793979441, 303.07828157347768]),
        ([45, 48], [13.680338402948045, 16.25462585886811], [46, 47],
         [14.538434221588068, 15.39653004022809]),
    ]  # fmt: skip
    for x, v, q, r in rec:
        assert interp1(np.array(x, float), np.array(v), np.array(q, float)).tolist() == r


def test_unported():
    with pytest.raises(NotImplementedError):
        interp1([1.0, 2.0], [1.0, 2.0], 1.5, "nearest")
    with pytest.raises(NotImplementedError):
        interp1([1.0, 2.0], np.float32([1.0, 2.0]), 1.5)
