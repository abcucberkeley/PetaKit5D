"""Tests of petakit5d.matlabcompat.graph: MATLAB's graph, minspantree, dfsearch and neighbors,
with MATLAB's tie rules.

The oracle cases build graphs with many equal weights (as stitch_local_assignment makes them:
every pair below the correlation threshold weighs -0.001), parallel edges, self-loops,
infinite and NaN weights and isolated nodes, and compare the edges, the trees of both methods
and types, the predecessors, the depth-first orders and the neighbors exactly; and they run
stitch_local_assignment itself against its steps written here with this module (its port is
the stitching task's). No MATLAB test covers these functions alone (StitchingTest uses
stitch_local_assignment through the stitching wrapper; the stitching task's).
"""

import numpy as np
import pytest

from petakit5d.matlabcompat import graph as G
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value

_C1 = np.uint64(0xBF58476D1CE4E5B9)
_C2 = np.uint64(0x94D049BB133111EB)
_GOLD = np.uint64(0x9E3779B97F4A7C15)


class Rand:
    """splitmix64 streams (as in test_lsq): each draw uses its own sub-seed"""

    def __init__(self, seed: int):
        self.seed = seed
        self.k = 0

    def u64(self, n: int) -> np.ndarray:
        self.k += 1
        s = np.uint64((self.seed * 1000003 + self.k * 7919) & 0xFFFFFFFFFFFFFFFF)
        with np.errstate(over="ignore"):
            z = s + np.arange(1, n + 1, dtype=np.uint64) * _GOLD
            z = (z ^ (z >> np.uint64(30))) * _C1
            z = (z ^ (z >> np.uint64(27))) * _C2
            return z ^ (z >> np.uint64(31))

    def rand(self, n: int) -> np.ndarray:
        return (self.u64(n) >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)

    def randi(self, lo: int, hi: int, n: int) -> np.ndarray:
        return lo + np.floor(self.rand(n) * (hi - lo + 1))


def grid_problem(seed, gx, gy, gz=1, *, filt=0.05, thr=0.25, corners=False):
    """pairs of neighboring tiles (and corner neighbors) with relative shifts and correlations,
    ``filt`` of them below the threshold"""
    r = Rand(seed)
    idx = np.array(
        [(x, y, z) for z in range(gz) for x in range(gx) for y in range(gy)], dtype=np.float64
    ).reshape(-1, 3)
    nF = len(idx)
    pairs = []
    for i in range(nF):
        dd = np.abs(idx[i + 1 :] - idx[i])
        ok = dd.sum(axis=1) == 1
        if corners:
            ok |= (dd.max(axis=1) == 1) & (dd[:, 2] == 0)
        pairs.extend((i, i + 1 + j) for j in np.nonzero(ok)[0])
    pr = np.array(pairs, dtype=np.float64).reshape(-1, 2) + 1
    nP = len(pr)
    rel = r.randi(-30, 30, 3 * nP).reshape(nP, 3) + r.randi(0, 7, 3 * nP).reshape(nP, 3) / 8.0
    xc = 0.3 + 0.7 * r.rand(nP)
    xc[r.rand(nP) < filt] = 0.1
    return {
        "nF": float(nF),
        "mx": np.column_stack([pr, xc]),
        "ash": np.column_stack([pr, rel]),
        "thr": float(thr),
    }


def random_graph(seed, n, kind):
    """(s, t, w) 1-based, as MATLAB takes them"""
    r = Rand(seed)
    m = int(n * (0.6 if kind == "forest" else 1.5)) + 1
    s = r.randi(1, n, m)
    t = r.randi(1, n, m)
    if kind not in ("selfloop", "multi"):
        keep = s != t
        s, t = s[keep], t[keep]
    if kind not in ("multi", "selfloop"):
        key = np.minimum(s, t) * 100000 + np.maximum(s, t)
        _, first = np.unique(key, return_index=True)
        first = np.sort(first)
        s, t = s[first], t[first]
    if s.size == 0:
        s, t = np.array([1.0]), np.array([2.0])
    if kind == "cont":
        w = -(0.3 + 0.7 * r.rand(s.size))
    elif kind == "allsame":
        w = -0.001 * np.ones(s.size)
    elif kind == "inf":
        w = np.array([np.inf, -np.inf, 1.0, -1.0])[r.randi(0, 3, s.size).astype(np.int64)]
    elif kind == "nan":
        w = np.array([np.nan, 1.0, 2.0])[r.randi(0, 2, s.size).astype(np.int64)]
    else:
        w = -np.array([0.001, 0.5, 0.7])[r.randi(0, 2, s.size).astype(np.int64)]
    return s, t, w


_GRAPH_CODE = """
s = in{1}; t = in{2}; w = in{3}; nn = in{4};
if isempty(nn), Gr = graph(s, t, w); else, Gr = graph(s, t, w, nn); end
[T, p] = minspantree(Gr, 'Type', 'forest');
[T2, p2] = minspantree(Gr);
[T3, p3] = minspantree(Gr, 'Type', 'forest', 'Method', 'sparse');
[T4, p4] = minspantree(Gr, 'Method', 'sparse');
r = 1 + mod(in{5}, numnodes(Gr));
[T5, p5] = minspantree(Gr, 'Type', 'forest', 'Root', r);
nb = cell(1, numnodes(Gr));
for k = 1 : numnodes(Gr)
    nb{k} = neighbors(Gr, k);
end
if numnodes(Gr) > 0
    v = dfsearch(T, 1); vG = dfsearch(Gr, 1); vr = dfsearch(Gr, r);
else
    v = []; vG = []; vr = [];
end
out = {Gr.Edges.EndNodes, Gr.Edges.Weight, T.Edges.EndNodes, T.Edges.Weight, p, ...
       T2.Edges.EndNodes, p2, T3.Edges.EndNodes, p3, T4.Edges.EndNodes, p4, ...
       T5.Edges.EndNodes, p5, v, vG, vr, nb};
"""


def _ends(T):
    return (T.edges + 1).astype(np.float64).reshape(-1, 2)


def _col(v):
    return (np.asarray(v, dtype=np.float64) + 1).reshape(-1, 1)


def py_graph(s, t, w, nn, root):
    g = G.graph(
        np.asarray(s).reshape(-1) - 1,
        np.asarray(t).reshape(-1) - 1,
        np.asarray(w).reshape(-1),
        None if nn is None else int(nn),
    )
    out = [_ends(g), g.weights.reshape(-1, 1)]
    r0 = int(root) % g.numnodes
    T, p = g.minspantree(type="forest")
    out += [_ends(T), T.weights.reshape(-1, 1), (p + 1).reshape(1, -1)]
    for kw in (
        {},
        {"type": "forest", "method": "sparse"},
        {"method": "sparse"},
        {"type": "forest", "root": r0},
    ):
        Tk, pk = g.minspantree(**kw)
        out += [_ends(Tk), (pk + 1).reshape(1, -1)]
    if g.numnodes:
        out += [_col(T.dfsearch(0)), _col(g.dfsearch(0)), _col(g.dfsearch(r0))]
    else:
        out += [None, None, None]
    out.append([_col(g.neighbors(k)) for k in range(g.numnodes)])
    return out


def _graph_case(cid, make, nn=None, root=0):
    s, t, w = make()
    nan = bool(np.any(np.isnan(w)))
    return Case(
        cid,
        code=_GRAPH_CODE,
        args=[In("s"), In("t"), In("w"), nn, float(root)],
        nargout=17,
        python=py_graph,
        inputs={
            "s": Value(lambda: make()[0][:, None]),
            "t": Value(lambda: make()[1][:, None]),
            "w": Value(lambda: make()[2][:, None]),
        },
        # minspantree refuses NaN weights
        expect_error="MATLAB:graphfun:graphbuiltin:NaNWeight" if nan else False,
        py_error=G.GraphError,
    )


def _graph_cases():
    cs = []
    s = 7000
    for kind in ("cont", "ties3", "allsame", "multi", "forest", "selfloop", "inf", "nan"):
        for n in (2, 3, 5, 8, 13, 40, 120):
            s += 1
            root = s % n
            cs.append(
                _graph_case(
                    f"graph_{kind}_{n}_s{s}",
                    lambda s=s, n=n, k=kind: random_graph(s, n, k),
                    root=root,
                )
            )
    # more nodes than the edges name; node 1 isolated
    cs.append(
        _graph_case(
            "graph_numnodes",
            lambda: (np.array([2.0, 3.0]), np.array([3.0, 4.0]), np.array([1.0, 1.0])),
            nn=6.0,
            root=3,
        )
    )
    cs.append(
        _graph_case(
            "graph_scalar_t",
            lambda: (np.array([1.0, 2.0, 4.0]), np.array([3.0]), np.array([-1.0, -1.0, -2.0])),
        )
    )
    cs.append(
        _graph_case("graph_one_edge", lambda: (np.array([1.0]), np.array([2.0]), np.array([0.5])))
    )
    cs.append(
        _graph_case(
            "graph_selfloop_only", lambda: (np.array([1.0]), np.array([1.0]), np.array([0.5]))
        )
    )
    return cs


_ERRORS = [
    ("err_zero_node", "graph([0 1], [1 2])"),
    ("err_frac_node", "graph(1.5, 2)"),
    ("err_t_zero", "graph([1 2], [0 1])"),
    ("err_sizes", "graph([1 2 3], [2 3])"),
    ("err_weight_size", "graph([1 2], [2 3], [1 1 1])"),
    ("err_numnodes", "graph([1 2], [2 3], [1 1], 2)"),
    ("err_dfs_node", "dfsearch(graph([1 2], [2 3]), 4)"),
    ("err_neighbors_node", "neighbors(graph([1 2], [2 3]), 0)"),
    ("err_root", "minspantree(graph([1 2], [2 3]), 'Root', 5)"),
    ("err_type", "minspantree(graph([1 2], [2 3]), 'Type', 'x')"),
    ("err_method", "minspantree(graph([1 2], [2 3]), 'Method', 'x')"),
]


def _py_error(expr):
    calls = {
        "err_zero_node": lambda: G.graph([-1, 0], [0, 1]),
        "err_frac_node": lambda: G.graph([0.5], [1]),
        "err_t_zero": lambda: G.graph([0, 1], [-1, 0]),
        "err_sizes": lambda: G.graph([0, 1, 2], [1, 2]),
        "err_weight_size": lambda: G.graph([0, 1], [1, 2], [1, 1, 1]),
        "err_numnodes": lambda: G.graph([0, 1], [1, 2], [1, 1], 2),
        "err_dfs_node": lambda: G.graph([0, 1], [1, 2]).dfsearch(3),
        "err_neighbors_node": lambda: G.graph([0, 1], [1, 2]).neighbors(-1),
        "err_root": lambda: G.graph([0, 1], [1, 2]).minspantree(root=4),
        "err_type": lambda: G.graph([0, 1], [1, 2]).minspantree(type="x"),
        "err_method": lambda: G.graph([0, 1], [1, 2]).minspantree(method="x"),
    }
    return calls[expr]


def _error_cases():
    return [
        Case(
            f"graph_{cid}",
            code=f"out = {{{expr}}};",
            python=lambda c=cid: _py_error(c)(),
            expect_error=True,
            py_error=G.GraphError,
        )
        for cid, expr in _ERRORS
    ]


# --- stitch_local_assignment: MST and DFS, each tile placed from its predecessor ---------------


def local_assignment(nF, mx, ash, xcorr_thresh):
    """stitch_local_assignment's steps, written with this module (0-based inside)"""
    mx = mx.copy()
    mx[mx[:, 2] < xcorr_thresh, 2] = 0.001
    g = G.graph(mx[:, 0].astype(np.int64) - 1, mx[:, 1].astype(np.int64) - 1, -mx[:, 2])
    T, _ = g.minspantree(type="forest")
    v = T.dfsearch(0)
    d = np.zeros((int(nF), 3))
    seen = np.zeros(int(nF), bool)
    for i, n in enumerate(v):
        if i == 0:
            seen[n] = True
            continue
        nb = T.neighbors(n)
        pr = int(nb[seen[nb]][0])
        if n > pr:
            s_, t_, sign = pr, n, 1.0
        else:
            s_, t_, sign = n, pr, -1.0
        ind = np.nonzero((ash[:, 0] == s_ + 1) & (ash[:, 1] == t_ + 1))[0][0]
        d[n] = d[pr] + sign * ash[ind, 2:5]
        seen[n] = True
    return d


def _local_case(cid, seed, gx, gy, gz, kw):
    def prob():
        return grid_problem(seed, gx, gy, gz, **kw)

    return Case(
        cid,
        code="p = in{1}; ov = false(p.nF); d = stitch_local_assignment(p.nF, p.mx, p.ash, ov, "
        "p.thr); out = {d};",
        args=[In("p")],
        python=lambda p: local_assignment(p["nF"], p["mx"], p["ash"], p["thr"]),
        inputs={"p": Value(prob)},
    )


def _local_cases():
    cs = []
    s = 8000
    for gx, gy, gz in [
        (2, 1, 1),
        (2, 2, 1),
        (3, 3, 1),
        (4, 3, 1),
        (5, 5, 1),
        (8, 6, 1),
        (10, 10, 1),
        (3, 3, 2),
        (5, 4, 3),
        (20, 20, 1),
    ]:
        for name, kw in (
            ("conn", {"filt": 0.0}),
            ("ties", {"filt": 0.5}),
            ("allbelow", {"filt": 1.0}),
            ("corners", {"filt": 0.3, "corners": True}),
            ("thr", {"thr": 0.7, "filt": 0.0}),
        ):
            s += 1
            cs.append(_local_case(f"local_{name}_{gx}x{gy}x{gz}_s{s}", s, gx, gy, gz, kw))
    return cs


CASES = _graph_cases() + _error_cases() + _local_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- unit tests --------------------------------------------------------------------------------


def test_edges_are_sorted_and_parallel_edges_keep_their_order():
    g = G.graph([2, 0, 1, 1, 0], [0, 1, 0, 1, 2], [5.0, 4.0, 3.0, 2.0, 1.0])
    assert g.edges.tolist() == [[0, 1], [0, 1], [0, 2], [0, 2], [1, 1]]
    assert g.weights.tolist() == [4.0, 3.0, 5.0, 1.0, 2.0]
    assert g.numnodes == 3


def test_prim_ties_go_to_the_lowest_node_and_the_first_edge():
    # a square with equal weights: from 0, nodes 1 and 3 tie (1 first), then 2 is reached from
    # 1 first, and 3's edge from 0 is kept (not strictly larger from 2)
    g = G.graph([0, 1, 2, 0], [1, 2, 3, 3], [1.0, 1.0, 1.0, 1.0])
    T, pred = g.minspantree()
    assert T.edges.tolist() == [[0, 1], [0, 3], [1, 2]]
    assert pred.tolist() == [-1, 0, 1, 0]


def test_infinite_weights_reach_nodes():
    g = G.graph([0, 1, 2, 0], [1, 2, 3, 3], [np.inf, -np.inf, 1.0, np.inf])
    T, pred = g.minspantree(type="forest")
    assert T.edges.tolist() == [[0, 1], [1, 2], [2, 3]]
    assert pred.tolist() == [-1, 0, 1, 2]


def test_forest_restarts_at_the_lowest_node():
    g = G.graph([1, 2], [2, 3], [1.0, 1.0], num_nodes=6)
    T, pred = g.minspantree(type="forest")
    assert T.edges.tolist() == [[1, 2], [2, 3]]
    assert pred.tolist() == [-1, -1, 1, 2, -1, -1]
    _, p2 = g.minspantree(root=1)
    assert np.isnan(p2[[0, 4, 5]]).all() and p2[[1, 2, 3]].tolist() == [-1, 1, 2]
    from petakit5d.core.dtypes import NAN

    assert p2[0].tobytes() == np.float64(NAN).tobytes()


def test_dfsearch_order_and_depth():
    g = G.graph([0, 0, 1, 1, 2], [1, 2, 3, 4, 5])
    assert g.dfsearch(0).tolist() == [0, 1, 3, 4, 2, 5]
    n = 20000
    path = G.graph(np.arange(n - 1), np.arange(1, n))
    assert path.dfsearch(0).size == n  # no recursion limit


def test_neighbors_once_and_self_loops():
    g = G.graph([0, 0, 1, 2, 3], [1, 1, 2, 2, 4], [5, 4, 3, 2, 1])
    assert g.neighbors(0).tolist() == [1]
    assert g.neighbors(2).tolist() == [1, 2]
    assert g.neighbors(1).tolist() == [0, 2]


def test_kruskal_is_stable():
    g = G.graph([0, 1, 0], [1, 2, 2], [1.0, 1.0, 1.0])
    T, pred = g.minspantree(method="sparse")
    assert T.edges.tolist() == [[0, 1], [0, 2]]
    assert pred.tolist() == [-1, 0, 0]


def test_errors_carry_matlab_ids():
    with pytest.raises(G.GraphError) as e:
        G.graph([0, 1], [1, 2], [1, 1], 2)
    assert e.value.identifier == "MATLAB:graphfun:graph:InvalidNumNodes"


def test_spanning_forest_of_1600_tiles_is_fast():
    import time

    p = grid_problem(5, 40, 40, filt=0.5)
    t0 = time.perf_counter()
    local_assignment(p["nF"], p["mx"], p["ash"], p["thr"])
    assert time.perf_counter() - t0 < 5.0
