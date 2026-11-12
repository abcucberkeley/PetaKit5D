"""MATLAB's undirected ``graph`` with ``minspantree``, ``dfsearch`` and ``neighbors``, as the
stitching shift assignment (``stitch_local_assignment``) uses them, with MATLAB's tie rules.

Node ids are 0-based (MATLAB's minus one). The rules, measured against MATLAB R2024b:

- ``graph(s, t, w)``: an edge per pair (a scalar ``s``, ``t`` or ``w`` is repeated), stored with
  its ends ordered (smaller first) and the edges sorted by their ends; parallel edges (a
  multigraph) keep their input order. The graph has ``max(s, t) + 1`` nodes unless ``num_nodes``
  is given. Weights may be infinite or NaN (but ``minspantree`` refuses NaN, as MATLAB).
- ``minspantree`` with the default method ``'dense'`` is Prim's algorithm: the next node is the
  one with the smallest connecting weight, the lowest node id among equal weights; the first edge
  that reaches a node sets its weight (also an infinite one), and it changes only for a strictly
  smaller weight (so among equal weights the edge from the node added first wins, and among
  parallel edges the first). The tree grows from
  ``root``; with ``type='forest'`` it then restarts at the lowest node not yet reached. Self-loops
  are never in a tree. ``method='sparse'`` is Kruskal's algorithm on the weights sorted stably.
- ``pred``: the parent of each node in the tree (float64), -1 for roots, NaN (MATLAB's) for
  nodes the tree doesn't reach (``type='tree'``).
- ``dfsearch(G, s)``: the nodes reachable from ``s`` in the order a depth-first search
  discovers them, visiting the neighbors of a node in ascending order.
- ``neighbors(G, n)``: the nodes joined to ``n`` by an edge, ascending, each once.
"""

from __future__ import annotations

import numpy as np

__all__ = ["Graph", "graph", "minspantree", "dfsearch", "neighbors", "GraphError"]


class GraphError(ValueError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


class Graph:
    """An undirected weighted graph as MATLAB stores it: ``edges`` (m x 2, int64, smaller end
    first, sorted), ``weights`` (m,), ``numnodes``."""

    def __init__(self, s, t, weights=None, num_nodes=None):
        s = _node_ids(s, "s", "MATLAB:graphfun:graphbuiltin:InvalidSRC")
        t = _node_ids(t, "t", "MATLAB:graphfun:graphbuiltin:InvalidTARG")
        if s.size == 1 and t.size != 1:
            s = np.full(t.size, s[0])
        elif t.size == 1 and s.size != 1:
            t = np.full(s.size, t[0])
        if s.size != t.size:
            raise GraphError(
                "s and t must have the same number of elements",
                "MATLAB:graphfun:graphbuiltin:EqualNumel",
            )
        if weights is None:
            w = np.ones(s.size)
        else:
            w = np.asarray(weights, dtype=np.float64).reshape(-1)
            if w.size == 1 and s.size != 1:
                w = np.full(s.size, w[0])
            if w.size != s.size:
                raise GraphError(
                    "weights must have one value per edge",
                    "MATLAB:graphfun:graph:InvalidSizeWeight",
                )
        nn = int(max(s.max(initial=-1), t.max(initial=-1))) + 1
        if num_nodes is not None:
            if int(num_nodes) < nn:
                raise GraphError(
                    "num_nodes is smaller than the largest node id",
                    "MATLAB:graphfun:graph:InvalidNumNodes",
                )
            nn = int(num_nodes)
        a = np.minimum(s, t)
        b = np.maximum(s, t)
        order = np.lexsort((np.arange(a.size), b, a))
        self.edges = np.column_stack([a[order], b[order]]).astype(np.int64)
        self.weights = w[order]
        self.numnodes = nn
        self._adj = None

    @property
    def numedges(self) -> int:
        return int(self.edges.shape[0])

    def _adjacency(self):
        """per node: neighbor ids and edge indices, by neighbor then edge (self-loops left out)"""
        if self._adj is None:
            a, b = self.edges[:, 0], self.edges[:, 1]
            keep = a != b
            e = np.nonzero(keep)[0]
            u = np.concatenate([a[keep], b[keep]])
            v = np.concatenate([b[keep], a[keep]])
            ee = np.concatenate([e, e])
            order = np.lexsort((ee, v, u))
            u, v, ee = u[order], v[order], ee[order]
            ptr = np.searchsorted(u, np.arange(self.numnodes + 1))
            self._adj = (ptr, v, ee)
        return self._adj

    def neighbors(self, node) -> np.ndarray:
        n = _check_node(self, node)
        ptr, v, _ = self._adjacency()
        nb = np.unique(v[ptr[n] : ptr[n + 1]])
        if np.any((self.edges[:, 0] == n) & (self.edges[:, 1] == n)):
            nb = np.union1d(nb, [n])
        return nb.astype(np.int64)

    def minspantree(self, *, type="tree", root=0, method="dense"):
        """(tree as a Graph with the same number of nodes, pred)"""
        if type not in ("tree", "forest"):
            raise GraphError(
                "type must be 'tree' or 'forest'", "MATLAB:graphfun:minspantree:ParseType"
            )
        if method not in ("dense", "sparse"):
            raise GraphError(
                "method must be 'dense' or 'sparse'", "MATLAB:graphfun:minspantree:ParseMethod"
            )
        restart = type == "forest"
        if self.numnodes == 0:
            return _subgraph(self, np.zeros(0, np.int64)), np.zeros(0)
        root = _check_node(self, root)
        if np.any(np.isnan(self.weights)):
            raise GraphError(
                "Graph edge weights contain NaN.", "MATLAB:graphfun:graphbuiltin:NaNWeight"
            )
        if method == "dense":
            chosen, pred = _prim(self, root, restart)
        else:
            chosen, pred = _kruskal(self, root, restart)
        return _subgraph(self, chosen), pred

    def dfsearch(self, s) -> np.ndarray:
        s = _check_node(self, s)
        ptr, v, _ = self._adjacency()
        seen = np.zeros(self.numnodes, bool)
        out = [s]
        seen[s] = True
        # an explicit stack of (node, position in its neighbor list): the recursive order
        stack = [(s, ptr[s])]
        while stack:
            u, k = stack[-1]
            end = ptr[u + 1]
            while k < end and seen[v[k]]:
                k += 1
            if k == end:
                stack.pop()
                continue
            stack[-1] = (u, k + 1)
            w = int(v[k])
            seen[w] = True
            out.append(w)
            stack.append((w, ptr[w]))
        return np.asarray(out, dtype=np.int64)


def graph(s, t, weights=None, num_nodes=None) -> Graph:
    """MATLAB's ``graph(s, t, weights)`` with 0-based node ids."""
    return Graph(s, t, weights, num_nodes)


def minspantree(G: Graph, *, type="tree", root=0, method="dense"):
    """``[T, pred] = minspantree(G, 'Type', type, 'Root', root, 'Method', method)``."""
    return G.minspantree(type=type, root=root, method=method)


def dfsearch(G: Graph, s) -> np.ndarray:
    """``v = dfsearch(G, s)``: node ids in discovery order."""
    return G.dfsearch(s)


def neighbors(G: Graph, node) -> np.ndarray:
    """``neighbors(G, node)``."""
    return G.neighbors(node)


def _node_ids(v, name, ident):
    a = np.asarray(v).reshape(-1)
    if a.size and (a.dtype.kind not in "iuf" or np.any(a != np.floor(a)) or np.any(a < 0)):
        raise GraphError(
            f"{name} must hold node ids (nonnegative integers)",
            ident,
        )
    return a.astype(np.int64)


def _check_node(G, n):
    n = int(n)
    if n < 0 or n >= G.numnodes:
        raise GraphError("node id out of range", "MATLAB:graphfun:graph:InvalidNodeID")
    return n


def _subgraph(G, chosen):
    T = Graph.__new__(Graph)
    chosen = np.sort(np.asarray(chosen, dtype=np.int64))
    T.edges = G.edges[chosen].copy()
    T.weights = G.weights[chosen].copy()
    T.numnodes = G.numnodes
    T._adj = None
    return T


def _prim(G, root, restart):
    """Prim's algorithm on an array of keys (MATLAB's 'dense' method)"""
    from petakit5d.core.dtypes import NAN

    n = G.numnodes
    ptr, v, ee = G._adjacency()
    w = G.weights
    key = np.full(n, np.inf)
    pedge = np.full(n, -1, np.int64)
    pred = np.full(n, NAN)
    intree = np.zeros(n, bool)
    cand = np.zeros(n, bool)  # reached, not in the tree
    chosen = []
    for r in [root] + [k for k in range(n) if k != root]:
        if intree[r]:
            continue
        if not restart and r != root:
            break
        pred[r] = -1.0
        u = r
        while True:
            intree[u] = True
            cand[u] = False
            if pedge[u] >= 0:
                chosen.append(pedge[u])
            for k in range(ptr[u], ptr[u + 1]):
                t = v[k]
                if intree[t]:
                    continue
                we = w[ee[k]]
                if not cand[t] or we < key[t]:
                    key[t] = we
                    pedge[t] = ee[k]
                    pred[t] = u
                    cand[t] = True
            idx = np.nonzero(cand)[0]
            if idx.size == 0:
                break
            kv = key[idx]
            ok = ~np.isnan(kv)
            # the smallest key, the lowest node among equal keys
            u = int(idx[ok][np.argmin(kv[ok])]) if np.any(ok) else int(idx[0])
    if restart:
        pred[np.isnan(pred)] = -1.0
    return np.asarray(chosen, np.int64), pred


def _kruskal(G, root, restart):
    """Kruskal's algorithm on the weights sorted stably (MATLAB's 'sparse' method); pred from a
    breadth-first search of the tree from root (and, for a forest, from each next lowest node)"""
    from petakit5d.core.dtypes import NAN

    n = G.numnodes
    a, b = G.edges[:, 0], G.edges[:, 1]
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    chosen = []
    for e in np.argsort(G.weights, kind="stable"):
        if a[e] == b[e]:
            continue
        ra, rb = find(a[e]), find(b[e])
        if ra != rb:
            parent[ra] = rb
            chosen.append(e)
    chosen = np.asarray(chosen, np.int64)
    T = _subgraph(G, chosen)
    ptr, v, _ = T._adjacency()
    pred = np.full(n, NAN)
    seen = np.zeros(n, bool)
    for r in [root] + [k for k in range(n) if k != root]:
        if seen[r]:
            continue
        if not restart and r != root:
            break
        seen[r] = True
        pred[r] = -1.0
        queue = [r]
        head = 0
        while head < len(queue):
            u = queue[head]
            head += 1
            for k in range(ptr[u], ptr[u + 1]):
                t = int(v[k])
                if not seen[t]:
                    seen[t] = True
                    pred[t] = u
                    queue.append(t)
    if not restart:
        # the tree of root only
        keep = seen[T.edges[:, 0]]
        chosen = np.sort(chosen)[keep]
    return chosen, pred
