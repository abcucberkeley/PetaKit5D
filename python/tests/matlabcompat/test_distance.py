"""Tests of petakit5d.matlabcompat.distance.bwdist against MATLAB R2024b, bit for bit (single
distances of 2D images up to 2048 x 2048 and of volumes; empty, full and numeric inputs; long
distances).

MATLAB tests whose intent is covered (their classes belong to other tasks):
DeconWrapperTest/largeFileErosionMaskMatchesDistanceMap (``bwdist(~bw) > e - 1``:
``bwdist.erosion_mask.*``); the per-slice distance maps of compute_tile_bwdist
(``bwdist.slices.*``).
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.matlabcompat import distance as D
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.synth import CounterRNG


def mask(seed: int, shape, p: float) -> np.ndarray:
    return np.asfortranarray(CounterRNG(seed, "n03-dist").uniform(shape) > p)


def py_bwdist(a):
    return to_matlab(D.bwdist(from_matlab(np.asarray(a))))


def edge_mask(shape, k: int) -> np.ndarray:
    """skewed-data-like masks: true inside a band, with a few holes"""
    a = np.zeros(shape, bool, order="F")
    y, x, z = np.indices(shape)
    a[(x > 0.5 * z + 2) & (x < 0.5 * z + shape[1] * 0.6)] = True
    h = CounterRNG(k, "n03-dist-holes").integers(0, a.size - 1, 5)
    a.reshape(-1, order="F")[h] = False
    return a


def py_erosion_mask(a, e):
    """RLdecon_large_in_memory: the border set false, then bwdist(~mask) > EdgeErosion - 1"""
    b = from_matlab(np.asarray(a)).copy()
    b[:, [0, -1], :] = False
    b[:, :, [0, -1]] = False
    b[[0, -1], :, :] = False
    return to_matlab(D.bwdist(~b) > e - 1)


EROSION_CODE = """
b = in{1};
b([1, end], :, :) = false;
b(:, [1, end], :) = false;
b(:, :, [1, end]) = false;
out = {bwdist(~b) > in{2} - 1};
"""

CASES: list[Case] = []
for k, (shape, p) in enumerate(
    [
        ((37, 41), 0.97),
        ((64, 64), 0.5),
        ((100, 120), 0.999),
        ((1, 50), 0.9),
        ((50, 1), 0.9),
        ((513, 600), 0.9995),
        ((256, 256), 0.99),
        ((7, 5), 0.0),
        ((300, 301), 0.99999),
    ]
):
    CASES.append(
        Case(
            f"bwdist.2d.{k}",
            matlab="bwdist",
            args=[In("A")],
            python=py_bwdist,
            inputs={"A": Value(lambda k=k, s=shape, p=p: mask(k, s, p))},
        )
    )
for n, p in [(1024, 0.9999), (2048, 0.99995), (2048, 0.9)]:
    CASES.append(
        Case(
            f"bwdist.2d.n{n}.p{str(p)[2:]}",
            matlab="bwdist",
            args=[In("A")],
            python=py_bwdist,
            inputs={"A": Value(lambda n=n, p=p: mask(n + int(p * 1e5), (n, n), p))},
        )
    )
for k, (shape, p) in enumerate(
    [((30, 35, 20), 0.995), ((64, 70, 40), 0.9999), ((10, 11, 12), 0.5)]
):
    CASES.append(
        Case(
            f"bwdist.3d.{k}",
            matlab="bwdist",
            args=[In("A")],
            python=py_bwdist,
            inputs={"A": Value(lambda k=k, s=shape, p=p: mask(100 + k, s, p))},
        )
    )
CASES += [
    Case(
        f"bwdist.{name}",
        matlab="bwdist",
        args=[In("A")],
        python=py_bwdist,
        inputs={"A": Value(make)},
    )
    for name, make in [
        ("empty_image", lambda: np.zeros((12, 13), bool, order="F")),
        ("full_image", lambda: np.ones((12, 13, 3), bool, order="F")),
        ("zero_size", lambda: np.zeros((0, 4), bool, order="F")),
        ("long_row", lambda: long_row()),
        ("corners3", lambda: corners((40, 50, 30))),
        ("uint16", lambda: np.asfortranarray((mask(201, (40, 45), 0.95) * 7).astype(np.uint16))),
        ("double", lambda: np.asfortranarray(CounterRNG(202, "n03-dist").uniform((30, 31)) - 0.99)),
    ]
]


def long_row() -> np.ndarray:
    a = np.zeros((2, 30000), bool, order="F")
    a[0, 0] = True
    return a


def corners(shape) -> np.ndarray:
    a = np.zeros(shape, bool, order="F")
    a[0, 0, 0] = True
    a[-1, -1, -1] = True
    return a


CASES += [
    Case(
        f"bwdist.erosion_mask.e{e}.{k}",
        code=EROSION_CODE,
        args=[In("A"), e],
        python=py_erosion_mask,
        inputs={"A": Value(lambda k=k, s=s: edge_mask(s, k))},
    )
    for k, s in enumerate([(40, 50, 30), (64, 80, 33)])
    for e in (1, 2, 5, 10)
]


def py_slices(a):
    """compute_tile_bwdist: the 2D distance map of each z slice"""
    v = from_matlab(np.asarray(a))
    return to_matlab(np.stack([D.bwdist(v[z]) for z in range(v.shape[0])]))


CASES.append(
    Case(
        "bwdist.slices.0",
        code=(
            "b = in{1}; d = zeros(size(b), 'single');"
            " for z = 1 : size(b, 3), d(:, :, z) = bwdist(b(:, :, z)); end; out = {d};"
        ),
        args=[In("A")],
        python=py_slices,
        inputs={"A": Value(lambda: mask(300, (60, 70, 9), 0.98))},
    )
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def brute(a: np.ndarray) -> np.ndarray:
    pts = np.argwhere(a)
    q = np.argwhere(np.ones(a.shape, bool))
    if len(pts) == 0:
        return np.full(a.shape, np.inf, np.float32)
    d2 = ((q[:, None, :] - pts[None, :, :]) ** 2).sum(-1).min(1)
    return np.sqrt(d2.astype(np.float32)).reshape(a.shape)


@pytest.mark.parametrize("shape", [(17, 23), (9, 10, 11), (1, 40), (5, 1, 7)])
def test_bwdist_equals_brute_force(shape):
    rng = np.random.default_rng(len(shape))
    for p in (0.5, 0.95, 1.0):
        a = rng.random(shape) > p
        assert np.array_equal(D.bwdist(a), brute(a))


def test_bwdist_does_not_depend_on_the_thread_count():
    a = np.random.default_rng(3).random((40, 300, 310)) > 0.999
    ref = D.bwdist(a, threads=1)
    for t in (2, 5):
        assert np.array_equal(D.bwdist(a, threads=t), ref)
