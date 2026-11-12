"""Tests of petakit5d.matlabcompat.label against MATLAB R2024b, bit for bit: bwlabeln's and
bwlabel's numbers, bwconncomp's structure (component order and pixel lists), labelmatrix,
bwareaopen and imfill 'holes', with connectivity 4, 8, 6, 18, 26 and connectivity arrays.

Oracle cases call MATLAB's built-ins themselves with arguments in MATLAB's form; the Python side
converts by hand (``from_matlab``, ``to_matlab``, MATLAB's linear indices) until the oracle's
roles (task F14) are in the tree.

MATLAB tests whose intent is covered (their classes belong to other tasks):
DeconBackProjectorTest/componentsMatchBwareaopenAndBwconncomp and
largestComponentTiesGoToTheFirst (bwconncomp's order and bwareaopen: ``bwconncomp.*``,
``bwareaopen.*``).
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.matlabcompat import label as L
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.synth import CounterRNG


def uniform(seed: int, shape, stream: str = "n03-label") -> np.ndarray:
    return CounterRNG(seed, stream).uniform(shape)


def mask(seed: int, shape, p: float) -> np.ndarray:
    return np.asfortranarray(uniform(seed, shape) > p)


def conn_array(seed: int, three_d: bool) -> np.ndarray:
    """a random symmetric connectivity with its center (MATLAB's form, double)"""
    c = uniform(seed, (3, 3, 3) if three_d else (3, 3), "n03-conn") > 0.5
    c = c | c[::-1, ::-1, ::-1] if three_d else c | c[::-1, ::-1]
    c[(1,) * c.ndim] = True
    return np.asfortranarray(c.astype(np.float64))


def py_conn(conn):
    """a MATLAB connectivity (scalar or matrix) in Python's form"""
    c = np.asarray(conn)
    return conn if c.size == 1 else from_matlab(c)


# ---------------------------------------------------------------------------------------------
# structured masks (MATLAB's form)
# ---------------------------------------------------------------------------------------------


def checker(shape) -> np.ndarray:
    idx = np.indices(shape).sum(axis=0)
    return np.asfortranarray(idx % 2 == 0)


def rings(n: int) -> np.ndarray:
    y, x = np.indices((n, n)) - (n - 1) / 2
    r = np.sqrt(y * y + x * x)
    return np.asfortranarray((np.floor(r) % 3) == 0)


def spiral(n: int) -> np.ndarray:
    a = np.zeros((n, n), bool)
    t, b, lft, rgt = 0, n - 1, 0, n - 1
    while t <= b and lft <= rgt:
        a[t, lft : rgt + 1] = True
        a[t : b + 1, rgt] = True
        a[b, lft : rgt + 1] = True
        a[t + 2 : b + 1, lft] = True
        t, b, lft, rgt = t + 2, b - 2, lft + 2, rgt - 2
        if lft <= rgt and t <= b:
            a[t - 1, lft] = True
    return np.asfortranarray(a)


def diagonals(shape) -> np.ndarray:
    idx = np.indices(shape)
    return np.asfortranarray(
        (idx[0] + 2 * idx[1] + (3 * idx[2] if len(shape) == 3 else 0)) % 5 == 0
    )


STRUCTURED = {
    "checker2": lambda: checker((15, 17)),
    "checker3": lambda: checker((7, 8, 6)),
    "rings": lambda: rings(41),
    "spiral": lambda: spiral(37),
    "diag2": lambda: diagonals((30, 33)),
    "diag3": lambda: diagonals((14, 15, 12)),
    "single": lambda: np.asfortranarray(np.eye(1, 9, 4, dtype=bool).reshape(3, 3)),
    "alltrue3": lambda: np.ones((5, 6, 4), bool, order="F"),
    "allfalse2": lambda: np.zeros((6, 7), bool, order="F"),
    "row": lambda: mask(11, (1, 50), 0.4),
    "col": lambda: mask(12, (50, 1), 0.4),
    "zcol": lambda: mask(13, (1, 1, 40), 0.4),
    "slab": lambda: mask(14, (1, 20, 9), 0.4),
    "empty": lambda: np.zeros((0, 4), bool, order="F"),
}

# ---------------------------------------------------------------------------------------------
# Python sides
# ---------------------------------------------------------------------------------------------


def py_bwlabeln(a, conn=None):
    lab, n = L.bwlabeln(from_matlab(np.asarray(a)), None if conn is None else py_conn(conn))
    return to_matlab(lab), n


def py_bwlabel(a, n=8):
    lab, num = L.bwlabel(from_matlab(np.asarray(a)), n)
    return to_matlab(lab), num


def py_bwconncomp(a, conn=None):
    return L.bwconncomp(
        from_matlab(np.asarray(a)), None if conn is None else py_conn(conn)
    ).to_matlab()


def py_labelmatrix(a, conn):
    return to_matlab(L.labelmatrix(L.bwconncomp(from_matlab(np.asarray(a)), py_conn(conn))))


def py_bwareaopen(a, p, conn=None):
    c = None if conn is None else py_conn(conn)
    return to_matlab(L.bwareaopen(from_matlab(np.asarray(a)), p, c))


def py_imfill(a, *args):
    args = [py_conn(x) if not isinstance(x, str) else x for x in args]
    return to_matlab(L.imfill(from_matlab(np.asarray(a)), *args))


# ---------------------------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------------------------

CONNS_2D = (4, 8)
CONNS_3D = (6, 18, 26)
LABEL_CASES: list[Case] = []


def label_case(cid, make, conn=None, func="bwlabeln", marks=()):
    pyf = {"bwlabeln": py_bwlabeln, "bwlabel": py_bwlabel}[func]
    return Case(
        cid,
        matlab=func,
        args=[In("A")] + ([] if conn is None else [conn]),
        nargout=2,
        outputs=("L", "num"),
        python=pyf,
        inputs={"A": Value(make)},
        marks=marks,
    )


for k in range(30):
    shape = (5 + (7 * k) % 50, 4 + (11 * k) % 60)
    p = 0.2 + 0.6 * ((k * 37) % 100) / 100
    for conn in CONNS_2D:
        LABEL_CASES.append(
            label_case(f"bwlabeln.2d.c{conn}.{k}", lambda k=k, s=shape, p=p: mask(k, s, p), conn)
        )
        LABEL_CASES.append(
            label_case(
                f"bwlabel.c{conn}.{k}", lambda k=k, s=shape, p=p: mask(k, s, p), conn, "bwlabel"
            )
        )
    LABEL_CASES.append(
        label_case(f"bwlabeln.2d.default.{k}", lambda k=k, s=shape, p=p: mask(100 + k, s, p))
    )
for k in range(30):
    shape = (3 + (5 * k) % 30, 2 + (7 * k) % 33, 2 + (3 * k) % 25)
    p = 0.25 + 0.5 * ((k * 41) % 100) / 100
    for conn in CONNS_3D:
        LABEL_CASES.append(
            label_case(
                f"bwlabeln.3d.c{conn}.{k}", lambda k=k, s=shape, p=p: mask(200 + k, s, p), conn
            )
        )
    LABEL_CASES.append(
        label_case(f"bwlabeln.3d.default.{k}", lambda k=k, s=shape, p=p: mask(300 + k, s, p))
    )
for k in range(16):
    three = k % 2 == 0
    shape = (19, 23, 11) if three else (31, 29)
    LABEL_CASES.append(
        Case(
            f"bwlabeln.connarray.{k}",
            matlab="bwlabeln",
            args=[In("A"), In("conn")],
            nargout=2,
            python=py_bwlabeln,
            inputs={
                "A": Value(lambda k=k, s=shape: mask(400 + k, s, 0.55)),
                "conn": Value(lambda k=k, t=three: conn_array(500 + k, t)),
            },
        )
    )
# 2D connectivity on 3D images and 3D connectivity on 2D images
for conn, shape in [
    (8, (20, 21, 6)),
    (4, (20, 21, 6)),
    (26, (40, 41)),
    (6, (40, 41)),
    (18, (33, 30)),
]:
    LABEL_CASES.append(
        label_case(
            f"bwlabeln.mixed.c{conn}.{len(shape)}d",
            lambda s=shape: mask(600 + sum(s), s, 0.5),
            conn,
        )
    )
for name, make in STRUCTURED.items():
    for conn in (4, 8, 6, 26):
        LABEL_CASES.append(label_case(f"bwlabeln.{name}.c{conn}", make, conn))
LABEL_CASES += [
    label_case("bwlabeln.large.c26", lambda: mask(700, (200, 210, 100), 0.62), 26),
    label_case("bwlabeln.large.c6", lambda: mask(701, (180, 190, 90), 0.55), 6),
    label_case("bwlabeln.large2d.c8", lambda: mask(702, (1024, 1100), 0.45), 8),
]

CC_CASES: list[Case] = []
for k in range(6):
    for conn in (None, 4, 8, 6, 26):
        three = conn in (6, 26) or (conn is None and k % 2)
        shape = (17 + k, 13 + 2 * k, 5 + k % 4) if three else (25 + 3 * k, 22 + k)
        CC_CASES.append(
            Case(
                f"bwconncomp.{'3d' if three else '2d'}.c{conn or 'default'}.{k}",
                matlab="bwconncomp",
                args=[In("A")] + ([] if conn is None else [conn]),
                python=py_bwconncomp,
                inputs={"A": Value(lambda k=k, s=shape: mask(800 + k, s, 0.5 + 0.03 * k))},
            )
        )
CC_CASES += [
    Case(
        f"bwconncomp.connarray.{k}",
        matlab="bwconncomp",
        args=[In("A"), In("conn")],
        python=py_bwconncomp,
        inputs={
            "A": Value(lambda k=k: mask(900 + k, (14, 15, 9) if k % 2 == 0 else (21, 24), 0.5)),
            "conn": Value(lambda k=k: conn_array(950 + k, k % 2 == 0)),
        },
    )
    for k in range(6)
]
CC_CASES += [
    Case(
        f"bwconncomp.{name}",
        matlab="bwconncomp",
        args=[In("A"), conn],
        python=py_bwconncomp,
        inputs={"A": Value(make)},
    )
    for name, make, conn in [
        (
            "cross4",
            lambda: np.asfortranarray(np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], float)),
            8,
        ),
        ("allfalse", lambda: np.zeros((5, 6, 3), bool, order="F"), 26),
        ("empty", lambda: np.zeros((0, 3), bool, order="F"), 8),
        ("checker3", lambda: checker((6, 7, 5)), 18),
    ]
]
CC_CASES += [
    Case(
        "bwconncomp.matrix4",
        matlab="bwconncomp",
        args=[In("A"), In("conn")],
        python=py_bwconncomp,
        inputs={
            "A": Value(lambda: mask(990, (20, 22), 0.5)),
            "conn": Value(
                lambda: np.asfortranarray(np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], float))
            ),
        },
    ),
    Case(
        "labelmatrix.2d",
        code="out = {labelmatrix(bwconncomp(in{1}, in{2}))};",
        args=[In("A"), 8],
        python=py_labelmatrix,
        inputs={"A": Value(lambda: mask(991, (30, 31), 0.5))},
    ),
    Case(
        "labelmatrix.3d.many",
        code="out = {labelmatrix(bwconncomp(in{1}, in{2}))};",
        args=[In("A"), 6],
        python=py_labelmatrix,
        inputs={"A": Value(lambda: mask(992, (40, 41, 20), 0.75))},
    ),
]

AREA_CASES: list[Case] = []
for k in range(10):
    three = k % 2 == 1
    shape = (25, 30, 12) if three else (48, 47)  # over 2 KB: hashes in the manifest
    for p in (0, 1, 3, 8, 30):
        AREA_CASES.append(
            Case(
                f"bwareaopen.{'3d' if three else '2d'}.p{p}.{k}",
                matlab="bwareaopen",
                args=[In("A"), p],
                python=py_bwareaopen,
                inputs={"A": Value(lambda k=k, s=shape: mask(1000 + k, s, 0.55 + 0.02 * k))},
            )
        )
    for conn in (6, 18) if three else (4,):
        AREA_CASES.append(
            Case(
                f"bwareaopen.{'3d' if three else '2d'}.p5.c{conn}.{k}",
                matlab="bwareaopen",
                args=[In("A"), 5, conn],
                python=py_bwareaopen,
                inputs={"A": Value(lambda k=k, s=shape: mask(1100 + k, s, 0.55))},
            )
        )
AREA_CASES.append(
    # a component of exactly P voxels is kept (DeconBackProjectorTest's threshold case)
    Case(
        "bwareaopen.exact_size",
        matlab="bwareaopen",
        args=[In("A"), 4, 26],
        python=py_bwareaopen,
        inputs={"A": Value(lambda: exact_size_mask())},
    )
)


def exact_size_mask() -> np.ndarray:
    a = np.zeros((10, 10, 10), bool, order="F")
    a[1:3, 1:3, 1] = True  # 4 voxels: kept
    a[6:9, 6, 6] = True  # 3 voxels: removed
    return a


FILL_CASES: list[Case] = []
for k in range(8):
    three = k % 2 == 1
    shape = (21, 23, 11) if three else (35, 37)
    FILL_CASES.append(
        Case(
            f"imfill.holes.{'3d' if three else '2d'}.{k}",
            matlab="imfill",
            args=[In("A"), "holes"],
            python=py_imfill,
            inputs={"A": Value(lambda k=k, s=shape: mask(1200 + k, s, 0.35 + 0.05 * k))},
        )
    )
    for conn in (6, 26) if three else (4, 8):
        FILL_CASES.append(
            Case(
                f"imfill.holes.c{conn}.{'3d' if three else '2d'}.{k}",
                matlab="imfill",
                args=[In("A"), conn, "holes"],
                python=py_imfill,
                inputs={"A": Value(lambda k=k, s=shape: mask(1300 + k, s, 0.4))},
            )
        )
FILL_CASES += [
    Case(
        f"imfill.{name}",
        matlab="imfill",
        args=[In("A"), opt],
        python=py_imfill,
        inputs={"A": Value(make)},
    )
    for name, make, opt in [
        ("rings_hole", lambda: rings(31), "hole"),
        ("spiral", lambda: spiral(29), "holes"),
        ("slab", lambda: mask(1401, (1, 20, 9), 0.3), "holes"),
        ("zcol", lambda: mask(1402, (1, 1, 30), 0.3), "holes"),
        ("onez", lambda: np.asfortranarray(rings(21)[:, :, None]), "holes"),
    ]
]

ERROR_CASES = [
    Case(
        f"error.{func}.{name}",
        matlab=func,
        args=[In("A"), *extra],
        nargout=nout,
        python=pyf,
        inputs={"A": Value(lambda: mask(1500, (6, 7), 0.5))},
        expect_error=True,
        py_error=L.MatlabError,
    )
    for func, pyf, nout, name, extra in [
        ("bwlabeln", py_bwlabeln, 2, "conn5", [5]),
        ("bwlabeln", py_bwlabeln, 2, "conn2x2", [np.ones((2, 2))]),
        ("bwconncomp", py_bwconncomp, 1, "connvalue", [np.full((3, 3), 2.0)]),
        (
            "bwconncomp",
            py_bwconncomp,
            1,
            "nocenter",
            [np.array([[1.0, 0, 1], [0, 0, 0], [1, 0, 1]])],
        ),
        (
            "bwareaopen",
            py_bwareaopen,
            1,
            "nonsymmetric",
            [3, np.array([[1.0, 1, 0], [0, 1, 0], [0, 0, 0]])],
        ),
    ]
]

ALL_CASES = LABEL_CASES + CC_CASES + AREA_CASES + FILL_CASES + ERROR_CASES


@oracle.parametrize(ALL_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# tests without MATLAB
# ---------------------------------------------------------------------------------------------


def first_order(lab: np.ndarray) -> np.ndarray:
    """labels renumbered by first voxel in MATLAB's order (the (z, x, y) view)"""
    v = np.swapaxes(lab.reshape((1,) * (3 - lab.ndim) + lab.shape), -1, -2).reshape(-1)
    u, first = np.unique(v[v > 0], return_index=True)
    order = u[np.argsort(first)]
    perm = np.zeros(int(lab.max()) + 1, dtype=np.int64)
    perm[order] = np.arange(1, len(order) + 1)
    return perm[lab.astype(np.int64)]


@pytest.mark.parametrize("conn", [4, 8, 6, 18, 26])
def test_bwconncomp_order_and_components(conn):
    from scipy import ndimage

    rng = np.random.default_rng(conn)
    a = rng.random((9, 21, 23)) > 0.55 if conn in (6, 18, 26) else rng.random((40, 45)) > 0.5
    cc = L.bwconncomp(a, conn)
    lab = L.labelmatrix(cc).astype(np.int64)
    s = L.connectivity_matrix(conn)
    ref, n = ndimage.label(a, structure=s)
    assert cc.num_objects == n
    assert np.array_equal(lab, first_order(ref))
    for p in cc.pixel_idx_list:
        mi = L._matlab_index(p, a.shape)
        assert np.all(np.diff(mi) > 0)


def test_labels_do_not_depend_on_the_thread_count():
    rng = np.random.default_rng(7)
    a = rng.random((30, 40, 50)) > 0.6
    ref = L.bwconncomp(a, 26, threads=1)
    for t in (2, 3, 7):
        cc = L.bwconncomp(a, 26, threads=t)
        assert cc.num_objects == ref.num_objects
        assert all(
            np.array_equal(x, y) for x, y in zip(cc.pixel_idx_list, ref.pixel_idx_list, strict=True)
        )
        assert np.array_equal(L.bwlabeln(a, 6, threads=t)[0], L.bwlabeln(a, 6, threads=1)[0])


def test_bwlabeln_labels_are_bwconncomps_up_to_numbering():
    rng = np.random.default_rng(8)
    for conn in (4, 8, 6, 26):
        a = rng.random((12, 30, 31)) > 0.5
        lab, n = L.bwlabeln(a, conn)
        ref = L.labelmatrix(L.bwconncomp(a, conn)).astype(np.int64)
        assert n == ref.max()
        pairs = np.unique(np.stack([lab.reshape(-1), ref.reshape(-1)]), axis=1)
        assert pairs.shape[1] == n + 1  # a bijection (background included)


def test_imfill_and_bwareaopen_against_scipy():
    from scipy import ndimage

    rng = np.random.default_rng(9)
    a = rng.random((50, 55)) > 0.45
    assert np.array_equal(L.imfill(a, "holes"), ndimage.binary_fill_holes(a))
    lab, n = ndimage.label(a, structure=np.ones((3, 3)))
    sizes = np.bincount(lab.reshape(-1))
    keep = sizes >= 6
    keep[0] = False
    assert np.array_equal(L.bwareaopen(a, 6), keep[lab])


def test_errors_carry_matlab_identifiers():
    with pytest.raises(L.MatlabError) as e:
        L.bwlabeln(np.ones((3, 3), bool), 5)
    assert e.value.identifier == "images:validate:badScalarConn"
