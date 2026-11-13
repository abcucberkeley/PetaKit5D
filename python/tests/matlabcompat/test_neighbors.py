"""Tests of petakit5d.matlabcompat.neighbors: MATLAB's knnsearch (its default k-d tree and its
exhaustive search) and cart2sph / sph2cart.

The knnsearch cases are exact: indices and distances, with many exact ties (directions of
integer offsets, as the OMW back projector's edge points; duplicated points; integer grids),
K = 1 to 4, 1-D points, rows with NaN and K above the number of points. The cart2sph and
sph2cart cases are exact too (MATLAB's atan2, sin and cos are fdlibm's, which the module
reproduces): grid offsets, random values, bands of magnitudes from 2^-42 to 2^42, signed zeros,
infinities and subnormals, angles up to 1.6e6 and near multiples of pi/2. No MATLAB test covers
these functions alone (DeconBackProjectorTest covers the back projector; the decon task's).
"""

import numpy as np
import pytest

from petakit5d.matlabcompat import elementary as E
from petakit5d.matlabcompat import neighbors as nbrs
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.compare import exact, tol

_C1 = np.uint64(0xBF58476D1CE4E5B9)
_C2 = np.uint64(0x94D049BB133111EB)
_GOLD = np.uint64(0x9E3779B97F4A7C15)


def _u(seed: int, n: int) -> np.ndarray:
    """n uniforms in [0, 1) of splitmix64 started at seed"""
    with np.errstate(over="ignore"):
        z = np.uint64(seed) + np.arange(1, n + 1, dtype=np.uint64) * _GOLD
        z = (z ^ (z >> np.uint64(30))) * _C1
        z = (z ^ (z >> np.uint64(27))) * _C2
        z = z ^ (z >> np.uint64(31))
    return (z >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)


def directions(psz, seed, frac=7):
    """unit vectors of the offsets of a grid from its center (y, x, z), as the window's edge
    points and queries, with sqrt only: (edge points, queries)"""
    c = (np.array(psz, dtype=np.float64) + 1) / 2
    Y, X, Z = np.meshgrid(*(np.arange(1, n + 1, dtype=np.float64) for n in psz), indexing="ij")
    P = np.column_stack(
        [Y.ravel(order="F") - c[0], X.ravel(order="F") - c[1], Z.ravel(order="F") - c[2]]
    )
    r = np.sqrt((P[:, 0] * P[:, 0] + P[:, 1] * P[:, 1]) + P[:, 2] * P[:, 2])
    r[r == 0] = 1.0
    U = P / r[:, None]
    e = np.unique(np.floor(_u(seed, max(2, len(U) // frac)) * len(U)).astype(np.int64))
    return U[e], U


def _knn_case(cid, make, k, nsmethod=None):
    opt = f", 'NSMethod', '{nsmethod}'" if nsmethod else ""
    return Case(
        cid,
        code=f"[i, d] = knnsearch(in{{1}}, in{{2}}, 'K', {k}{opt}); out = {{i, d}};",
        args=[In("X"), In("Y")],
        nargout=2,
        python=lambda X, Y: _py_knn(X, Y, k),
        inputs={"X": Value(lambda: make()[0]), "Y": Value(lambda: make()[1])},
    )


def _py_knn(X, Y, k):
    i, d = nbrs.knnsearch(X, Y, k)
    return [(i + 1).astype(np.float64), d]


def _knn_cases():
    cs = []
    s = 300
    for psz in ((5, 5, 5), (6, 6, 4), (7, 9, 5), (12, 12, 12), (16, 14, 9), (21, 21, 11)):
        for k in (1, 2):
            s += 1
            cs.append(
                _knn_case(
                    f"knn_dirs_{'x'.join(map(str, psz))}_k{k}",
                    lambda psz=psz, s=s: directions(psz, s),
                    k,
                )
            )
    cs.append(
        _knn_case(
            "knn_dirs_12x12x12_exhaustive", lambda: directions((12, 12, 12), 303), 2, "exhaustive"
        )
    )
    for n, nq, k in (
        (10, 30, 2),
        (50, 200, 2),
        (200, 500, 3),
        (40, 100, 1),
        (40, 100, 4),
        (3, 9, 4),
        (1, 5, 2),
    ):
        s += 1

        def make(s=s, n=n, nq=nq):
            X = np.floor(_u(s, 3 * n) * 4).reshape(n, 3)
            Q = np.floor(_u(s + 1, 3 * nq) * 4).reshape(nq, 3) + 0.5 * (
                _u(s + 2, 3 * nq) < 0.3
            ).reshape(nq, 3)
            return X, Q

        cs.append(_knn_case(f"knn_intgrid_{n}_{nq}_k{k}", make, k))
    for n, nq in ((30, 100), (300, 1000), (2000, 3000)):
        s += 1

        def make(s=s, n=n, nq=nq):
            X = (_u(s, 3 * n) * 2 - 1).reshape(n, 3)
            X = np.vstack([X, X[:5]])
            Q = np.vstack([(_u(s + 1, 3 * nq) * 2 - 1).reshape(nq, 3), X[:10]])
            return X, Q

        cs.append(_knn_case(f"knn_rand_{n}_{nq}", make, 2))
    for n in (5, 30, 200):
        s += 1

        def make(s=s, n=n):
            X = np.unique(np.floor(_u(s, n) * 3 * n) + 1)[:, None]
            Q = (np.floor(_u(s + 1, 4 * n) * 3 * n) + 1)[:, None]
            return X, Q

        cs.append(_knn_case(f"knn_oned_{n}", make, 1))
    nanx = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [0, 0], [0.5, 0.5], [np.nan, 1], [9, 9]])
    nanq = np.array([[0.5, 0.5], [0, 0], [2, 2], [np.nan, 0]])
    cs.append(_knn_case("knn_nan_rows", lambda: (nanx, nanq), 4))
    cs.append(_knn_case("knn_k_above_n", lambda: (nanx[:6], nanq[:3]), 8))
    cs.append(
        _knn_case(
            "knn_rows_are_points", lambda: (np.array([1.0, 2.0, 3.0]), np.array([2.0, 2.0, 2.0])), 1
        )
    )
    cs.append(
        _knn_case(
            "knn_inf_point",
            lambda: (
                np.array([[0.0, 0.0], [np.inf, 0.0], [1.0, 1.0]]),
                np.array([[0.2, 0.1], [5.0, 5.0]]),
            ),
            3,
        )
    )
    for cid, X, Y, k in (
        ("dims", np.zeros((3, 2)), np.zeros((4, 3)), 1),
        ("k0", np.zeros((3, 2)), np.zeros((4, 2)), 0),
        ("kfrac", np.zeros((3, 2)), np.zeros((4, 2)), 1.5),
    ):
        cs.append(
            Case(
                f"knn_err_{cid}",
                code=f"[i, d] = knnsearch(in{{1}}, in{{2}}, 'K', {k});",
                args=[X, Y],
                python=lambda X, Y, k=k: nbrs.knnsearch(X, Y, k),
                expect_error=True,
                py_error=nbrs.KnnError,
            )
        )
    return cs


# --- cart2sph / sph2cart -----------------------------------------------------------------------

_SPH_CODE = """
[az, el, r] = cart2sph(in{1}(:, 1), in{1}(:, 2), in{1}(:, 3));
[x, y, z] = sph2cart(az, el, r);
[x1, y1, z1] = sph2cart(az, el, ones(size(r)));
out = {az, el, r, x, y, z, x1, y1, z1};
"""


def _py_sph(V):
    az, el, r = nbrs.cart2sph(V[:, 0], V[:, 1], V[:, 2])
    x, y, z = nbrs.sph2cart(az, el, r)
    x1, y1, z1 = nbrs.sph2cart(az, el, np.ones(r.shape))
    return [v.reshape(-1, 1) for v in (az, el, r, x, y, z, x1, y1, z1)]


def _sph_values(kind, part):
    if kind == "grid":
        g = np.arange(-6.0, 6.5, 0.5)
        Y, X, Z = np.meshgrid(g, g, g[::3], indexing="ij")
        V = np.column_stack([X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F")])
    elif kind == "rand":
        V = ((_u(4000, 3 * 2400) * 2 - 1) * 100).reshape(-1, 3)
    elif kind == "wide":
        # one band of magnitudes per case: 2^(b - 2) to 2^(b + 2)
        b = (-40, -20, -8, 8, 20, 40)[part]
        e = (np.floor(_u(4001 + part, 900) * 4) + b - 2).astype(np.int64)
        return np.ldexp(_u(4100 + part, 900) * 2 - 1, e).reshape(-1, 3)
    else:
        from petakit5d.core.dtypes import NAN

        v = np.array([0.0, -0.0, 1.0, -1.0, np.inf, -np.inf, 1e-310, -1e-310, 1e308, NAN])
        Y, X, Z = np.meshgrid(v, v, v, indexing="ij")
        V = np.column_stack([X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F")])
        nan = np.isnan(V).any(axis=1)
        if kind == "nan":
            # 200 points: inline in the manifest, so CI decides the comparison of NaN-ness
            return V[nan][part * 200 : (part + 1) * 200]
        V = V[~nan]
    # 300 points: each output is 2,400 bytes, kept in the manifest by its hash
    return V[part * 300 : (part + 1) * 300]


# MATLAB's NaN payloads where NaN meets NaN in sph2cart's products are not reproduced (hypot's
# positive NaN and atan2's MATLAB NaN are): NaN where MATLAB has NaN, other values exact
_NAN_ONLY = tol(reason="proposed (N06): NaN payloads of sph2cart with NaN inputs", abs=0.0)


def _sph_cases():
    cs = []
    for kind, parts in (("grid", 19), ("rand", 8), ("wide", 6), ("special", 3), ("nan", 2)):
        for part in range(parts):
            cs.append(
                Case(
                    f"sph_{kind}_{part}",
                    code=_SPH_CODE,
                    args=[In("V")],
                    nargout=9,
                    python=_py_sph,
                    inputs={"V": Value(lambda k=kind, p=part: _sph_values(k, p))},
                    compare=_NAN_ONLY if kind == "nan" else exact(),
                )
            )
    # sph2cart of any angles up to 2^20 pi/2 (fdlibm's medium argument reduction), multiples of
    # pi/2 (its third reduction step) and huge radii
    for part in range(8):
        cs.append(
            Case(
                f"sph2cart_angles_{part}",
                code="[x, y, z] = sph2cart(in{1}(:, 1), in{1}(:, 2), in{1}(:, 3));"
                " out = {x, y, z};",
                args=[In("A")],
                nargout=3,
                python=lambda A: [
                    v.reshape(-1, 1) for v in nbrs.sph2cart(A[:, 0], A[:, 1], A[:, 2])
                ],
                inputs={"A": Value(lambda p=part: _angles(p))},
            )
        )
    return cs


def _angles(part):
    n = 300
    lim = (4.0, 50.0, 1e3, 1e5, 1.6e6, 4.0, 1e3, 1e5)[part]
    a = (_u(5000 + part, 2 * n) * 2 - 1) * lim
    if part >= 5:
        # near multiples of pi/2, where the reduction needs its later steps
        k = (
            np.floor(_u(5100 + part, 2 * n) * (2 * lim / 1.5707963267948966))
            - lim / 1.5707963267948966
        )
        a = k * 1.5707963267948966 + (_u(5200 + part, 2 * n) - 0.5) * 2.0**-30
    r = (_u(5300 + part, n) + 0.5) * np.ldexp(1.0, (part * 37) % 120 - 60)
    return np.column_stack([a[:n], a[n:], r])


# --- the elementary functions alone, where fdlibm's branches meet ------------------------------


def _elem_values(kind, part):
    """300 arguments per case: (y, x) for atan2 and hypot, x for sin and cos"""
    u = _u(6000 + 31 * part + len(kind), 600)
    if kind == "atan2":
        # |y/x| next to atan's branch limits (7/16, 11/16, 19/16, 39/16), tiny and huge
        lims = (0.4375, 0.6875, 1.1875, 2.4375, 2.0**-29, 2.0**66, 1.0, 3.0)
        r = lims[part] * (1 + (u[:300] - 0.5) * 0.08)
        sx = np.where(u[300:] < 0.5, -1.0, 1.0)
        sy = np.where((u[300:] * 4) % 1 < 0.5, -1.0, 1.0)
        x = sx * (1 + np.floor(u[300:] * 7))
        return np.column_stack([sy * r * np.abs(x), x])
    if kind == "hypot":
        return (u.reshape(300, 2) * 2 - 1) * np.ldexp(1.0, part * 30 - 60)
    # sin and cos: the kernels' limits (2^-27, 0.3, 0.78125, pi/4), the first reduction step
    # (pi/2 to 3pi/4), the second and third (exact products k pi/2, the table of n < 32)
    if part == 0:
        x = np.ldexp(u[:300] + 0.5, -28 + (np.floor(u[300:] * 3)).astype(np.int64))
    elif part == 1:
        x = 0.3 + (u[:300] - 0.5) * 0.02
    elif part == 2:
        x = 0.78125 + u[:300] * 0.0042
    elif part == 3:
        x = 1.5707963267948966 + (u[:300] - 0.5) * 1.6
    elif part == 4:
        x = (np.floor(u[:300] * 31) + 1) * 1.5707963267948966
    elif part == 5:
        x = (np.floor(u[:300] * 100000) + 32) * 1.5707963267948966
    else:
        x = (np.floor(u[:300] * 1000000) + 1) * 1.5707963267948966 + (u[300:] - 0.5) * 2.0**-20
    return np.where(u[300:] < 0.3, -x, x)[:, None]


def _nan_pairs():
    from petakit5d.core.dtypes import NAN

    pos = float(np.array([0x7FF8000000000000], dtype=np.uint64).view(np.float64)[0])
    v = np.array([NAN, pos, 0.0, -1.5, 2.0, np.inf, -np.inf])
    Y, X = np.meshgrid(v, v, indexing="ij")
    V = np.column_stack([Y.ravel(order="F"), X.ravel(order="F")])
    return V[np.isnan(V).any(axis=1)]


def _elem_cases():
    cs = []
    for fn, nparts, code, py in (
        (
            "atan2",
            8,
            "out = {atan2(in{1}(:, 1), in{1}(:, 2))};",
            lambda V: E.atan2(V[:, 0], V[:, 1]).reshape(-1, 1),
        ),
        (
            "hypot",
            4,
            "out = {hypot(in{1}(:, 1), in{1}(:, 2))};",
            lambda V: E.hypot(V[:, 0], V[:, 1]).reshape(-1, 1),
        ),
        ("sin", 7, "out = {sin(in{1})};", lambda V: E.sin(V[:, 0]).reshape(-1, 1)),
        ("cos", 7, "out = {cos(in{1})};", lambda V: E.cos(V[:, 0]).reshape(-1, 1)),
    ):
        kind = "sincos" if fn in ("sin", "cos") else fn
        if fn in ("atan2", "hypot"):
            # NaN of both signs against finite values, Inf and NaN: MATLAB's NaN of atan2,
            # hypot's positive NaN
            cs.append(
                Case(
                    f"elem_{fn}_nan",
                    code=code,
                    args=[In("V")],
                    python=py,
                    inputs={"V": Value(_nan_pairs)},
                )
            )
        for part in range(nparts):
            cs.append(
                Case(
                    f"elem_{fn}_{part}",
                    code=code,
                    args=[In("V")],
                    python=py,
                    inputs={"V": Value(lambda k=kind, p=part: _elem_values(k, p))},
                )
            )
    return cs


CASES = _knn_cases() + _sph_cases() + _elem_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- unit tests --------------------------------------------------------------------------------


@pytest.mark.parametrize("seed", range(3))
def test_tree_search_equals_exhaustive(seed, monkeypatch):
    """the k-d tree only proposes candidates: with many duplicated points the result is the
    exhaustive one"""
    E, Q = directions((15, 13, 11), 50 + seed, frac=3)
    E = np.vstack([E, E[::7]])  # more exact duplicates
    for k in (1, 2, 3):
        ie, de = nbrs.knnsearch(E, Q, k)
        monkeypatch.setattr(nbrs, "_EXHAUSTIVE_MAX", 0)
        it, dt = nbrs.knnsearch(E, Q, k)
        monkeypatch.undo()
        assert np.array_equal(ie, it)
        assert de.tobytes() == dt.tobytes()


def test_num_threads_are_the_tree_queries_and_keep_the_result(monkeypatch):
    """num_threads (Python only; None: the CPUs available) go to the k-d tree's queries; the
    result is the exhaustive search's at any count"""
    from scipy import spatial

    from petakit5d.core.threads import available_cpus

    E, Q = directions((15, 13, 11), 53, frac=3)
    E = np.vstack([E, E[::7]])
    want_i, want_d = nbrs.knnsearch(E, Q, 2)  # exhaustive
    workers = []

    class Spy(spatial.cKDTree):
        def query(self, *args, **kwargs):
            workers.append(kwargs.get("workers"))
            return super().query(*args, **kwargs)

    monkeypatch.setattr(spatial, "cKDTree", Spy)
    monkeypatch.setattr(nbrs, "_EXHAUSTIVE_MAX", 0)
    avail = available_cpus()
    for nt, expect in ((1, 1), (3, min(3, avail)), (None, avail)):
        i, d = nbrs.knnsearch(E, Q, 2, num_threads=nt)
        assert np.array_equal(i, want_i) and d.tobytes() == want_d.tobytes()
        assert workers[-1] == expect
    with pytest.raises(ValueError, match="thread count"):
        nbrs.knnsearch(E, Q, 2, num_threads=-1)


def test_ties_go_to_the_lowest_index():
    X = np.array([[1.0, 0, 0]] * 12)
    X[4] = [0, 1, 0]
    i, d = nbrs.knnsearch(X, np.array([[1.0, 0, 0], [0.9, 0.1, 0]]), 3)
    assert i.tolist() == [[0, 1, 2], [0, 1, 2]]


def test_order_is_by_squared_distance():
    """two squared distances with the same square root: the smaller squared distance first"""
    q = np.array([[0.0]])
    a = 1.0
    b = np.nextafter(1.0, 2.0)  # b^2 and a^2 differ, sqrt of the squares may not
    X = np.array([[b], [a]])
    i, d = nbrs.knnsearch(X, q, 2)
    assert i.tolist() == [[1, 0]]


def test_nan_handling():
    from petakit5d.core.dtypes import NAN

    X = np.array([[np.nan, 0.0], [0.0, 0.0], [1.0, 1.0]])
    i, d = nbrs.knnsearch(X, np.array([[0.9, 0.9], [np.nan, 1.0]]), 2)
    assert i.tolist() == [[2, 1], [1, 2]]
    assert d[1].tobytes() == np.full(2, NAN).tobytes()


def test_one_dimensional_points():
    i, d = nbrs.knnsearch(np.array([[1.0], [2.0], [3.0]]), np.array([[2.5], [0.0]]))
    assert i[:, 0].tolist() == [1, 0]
    assert d[:, 0].tolist() == [0.5, 1.0]


def test_hypot_is_the_reference_platforms():
    """msun's hypot is the reference platform's C library's (glibc 2.28) and MATLAB's; glibc
    2.35 and later round the first three the other way (fdlibm 5.3's own version the next two:
    it does not read the scaled arguments again); then both scalings, and two subnormals of the
    same high word, larger second, which are swapped after the scaling"""
    pairs = [
        (1.7799444501991235e-05, 0.4917817083544318, 0.4917817086765465),
        (69806544.69088814, 69803484.08146675, 98719198.0903033),
        (1934.5564907201772, 967.2782453604309, 2162.899909782032),
        (2.186954278671333e-307, 1.4686997726e-314, 2.1869542786713378e-307),
        (1.0040762481921225e-306, 6.2987061748e-313, 1.0040762481923201e-306),
        (3e200, 4e200, 4.9999999999999995e200),
        (1.2345e300, 6.789e299, 1.4088631800142978e300),
        (4.564467125769093e180, 2.6971851197726455e180, 5.30180796639691e180),
        (3e-140, 4e-140, 5e-140),
        (1e-140, 9.999999990686774e-141, 1.4142135617145506e-140),
        (2.086006742350502e-308, 2.086008864346292e-308, 2.9500605267117925e-308),
        (1.582266161533255e-309, 1.58227677151221e-309, 2.237669767324993e-309),
        (2.1219957924e-314, 4.2439915814e-314, 4.7449268366e-314),
    ]
    x, y, want = (np.array(v) for v in zip(*pairs, strict=True))
    assert E.hypot(x, y).tobytes() == want.tobytes()
    assert E.hypot(-y, x).tobytes() == want.tobytes()
    assert E.hypot(3.0, 4.0) == 5.0 and E.hypot(np.ones((2, 3)), 1.0).shape == (2, 3)


def test_cart2sph_and_sph2cart_call_no_c_library_function(monkeypatch):
    """hypot, atan2, sin and cos are written out, so that the results are MATLAB's on any C
    library and CPU (NumPy's ufuncs call the C library's functions or SIMD versions of them;
    angles up to 2^20 pi/2)"""
    V = np.vstack([_sph_values("rand", 0), _sph_values("special", 0), _sph_values("nan", 0)])
    A = np.vstack([_angles(p) for p in range(8)])
    want = _py_sph(V) + [v.copy() for v in nbrs.sph2cart(A[:, 0], A[:, 1], A[:, 2])]

    def refuse(*args, **kwargs):
        raise AssertionError("a C library function was called")

    for name in ("hypot", "arctan2", "arctan", "sin", "cos", "tan", "exp", "log"):
        monkeypatch.setattr(np, name, refuse)
    got = _py_sph(V) + list(nbrs.sph2cart(A[:, 0], A[:, 1], A[:, 2]))
    assert all(g.tobytes() == w.tobytes() for g, w in zip(got, want, strict=True))


def test_cart2sph_round_trip_and_hypot():
    az, el, r = nbrs.cart2sph(np.array([3.0, 0.0]), np.array([4.0, 0.0]), np.array([12.0, 0.0]))
    assert r.tolist() == [13.0, 0.0]
    x, y, z = nbrs.sph2cart(az, el, r)
    assert np.allclose([x[0], y[0], z[0]], [3, 4, 12], rtol=1e-15, atol=0)
    assert np.array_equal(az, E.atan2(np.array([4.0, 0.0]), np.array([3.0, 0.0])))
