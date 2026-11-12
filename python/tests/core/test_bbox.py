import math

import numpy as np
import pytest

from petakit5d.core import bbox as bb


def test_box_from_and_to_matlab():
    # MATLAB [ymin xmin zmin ymax xmax zmax] = [2 3 4 5 6 7] selects a(2:5, 3:6, 4:7)
    box = bb.box_from_matlab([2, 3, 4, 5, 6, 7])
    assert box == (3, 1, 2, 7, 5, 6)
    assert bb.box_to_matlab(box) == (2, 3, 4, 5, 6, 7)
    assert bb.box_from_matlab([2, 3, 5, 6]) == (1, 2, 5, 6)  # 2D [ymin xmin ymax xmax]
    assert bb.box_from_matlab([1, 1, 1, math.inf, 10, math.inf]) == (0, 0, 0, None, None, 10)
    assert bb.box_to_matlab((0, 0, 0, None, None, 10)) == (1, 1, 1, math.inf, 10, math.inf)
    with pytest.raises(ValueError):
        bb.box_from_matlab([1, 2, 3])
    with pytest.raises(ValueError):
        bb.box_from_matlab([1.5, 1, 1, 2, 2, 2])


def test_the_box_selects_what_matlab_selects():
    from petakit5d.core.arrays import from_matlab

    m = np.arange(6 * 7 * 8).reshape((6, 7, 8), order="F")  # MATLAB's (y, x, z)
    mbox = [2, 3, 4, 5, 6, 7]
    sel_matlab = m[1:5, 2:6, 3:7]  # a(2:5, 3:6, 4:7)
    p = from_matlab(m)
    sel = bb.crop(p, bb.box_from_matlab(mbox))
    np.testing.assert_array_equal(sel, from_matlab(sel_matlab))
    assert np.shares_memory(sel, p)


@pytest.mark.parametrize("dim", [2, 3])
def test_conversions_are_exact_inverses_on_generated_values(dim):
    rng = np.random.default_rng(dim)
    for _ in range(10000):
        starts = rng.integers(-50, 1000, size=dim)
        stops = starts + rng.integers(0, 1000, size=dim)
        ends = [None if rng.random() < 0.1 else int(s) for s in stops]
        box = tuple(int(s) for s in starts) + tuple(ends)
        m = bb.box_to_matlab(box)
        assert bb.box_from_matlab(m) == box
        assert bb.box_to_matlab(bb.box_from_matlab(m)) == m
        r = (int(starts[0]), ends[0])
        assert bb.range_from_matlab(bb.range_to_matlab(r)) == r
    rows = rng.integers(0, 500, size=(10000, 6))
    rows[:, 3:] += rows[:, :3] + 1
    mrows = bb.boxes_to_matlab(rows)
    np.testing.assert_array_equal(bb.boxes_from_matlab(mrows), rows)
    for r, mr in zip(rows[:200], mrows[:200], strict=True):
        assert tuple(mr) == bb.box_to_matlab(tuple(r))


def test_ranges():
    assert bb.range_from_matlab([3, 7]) == (2, 7)
    assert bb.range_from_matlab([1, math.inf]) == (0, None)
    assert bb.range_to_matlab((2, 7)) == (3, 7)
    assert bb.range_to_matlab((0, None)) == (1, math.inf)


def test_slices_shape_full():
    assert bb.box_slices((1, 2, 3, 4, 5, 6)) == (slice(1, 4), slice(2, 5), slice(3, 6))
    assert bb.box_slices((1, 2, None, 5)) == (slice(1, None), slice(2, 5))
    assert bb.box_shape((1, 2, 3, 4, 5, 6)) == (3, 3, 3)
    assert bb.box_shape((0, 0, 0, None, 4, None), (7, 8, 9)) == (7, 4, 9)
    with pytest.raises(ValueError):
        bb.box_shape((0, 0, None, 1))
    assert bb.full_box((5, 6, 7)) == (0, 0, 0, 5, 6, 7)
    assert bb.full_box((6, 7)) == (0, 0, 6, 7)


def test_crop_2d_and_leading_axes():
    im = np.arange(30).reshape(5, 6)
    assert bb.crop(im, (1, 2, 3, 5)).shape == (2, 3)
    tiles = np.zeros((4, 5, 6, 7))
    assert bb.crop(tiles, (1, 1, 1, 3, 4, 5)).shape == (4, 2, 3, 4)
    with pytest.raises(ValueError, match="outside"):
        bb.crop(im, (0, 0, 6, 1))


def test_check_box():
    assert bb.check_box(np.array([1, 2, 3, 4, 5, 6])) == (1, 2, 3, 4, 5, 6)
    for bad in ([1, 2, 3], [1.5, 1, 1, 2, 2, 2], [3, 1, 1, 2, 2, 2]):
        with pytest.raises(ValueError):
            bb.check_box(bad)
    assert bb.check_box((0, 0, 0, 0, 0, 0)) == (0, 0, 0, 0, 0, 0)  # an empty box is a box
    with pytest.raises(ValueError, match="outside"):
        bb.check_box((-1, 0, 0, 2, 2, 2), (5, 5, 5))
    assert bb.check_box((-1, 0, 0, 9, 2, 2), (5, 5, 5), allow_outside=True) == (-1, 0, 0, 9, 2, 2)


def test_clamp_intersect_and_pad():
    assert bb.clamp_box((-3, 1, 0, 9, 4, None), (5, 6, 7)) == (0, 1, 0, 5, 4, 7)
    with pytest.raises(ValueError, match="overlap"):
        bb.clamp_box((6, 0, 0, 9, 2, 2), (5, 6, 7))
    assert bb.intersect_boxes((0, 0, 0, 4, 4, 4), (2, 1, 3, 9, 9, 9)) == (2, 1, 3, 4, 4, 4)
    assert bb.intersect_boxes((0, 0, 0, 2, 2, 2), (2, 2, 2, 4, 4, 4)) is None
    assert bb.intersect_boxes((0, 0, None, 5), (1, 1, 3, None)) == (1, 1, 3, 5)
    box = (-2, 1, 0, 7, 6, 9)
    pads = bb.pad_widths(box, (5, 6, 7))
    assert pads == ((2, 2), (0, 0), (0, 2))
    a = np.ones((5, 6, 7))
    padded = np.pad(bb.crop(a, bb.clamp_box(box, a.shape)), pads)
    assert padded.shape == bb.box_shape(box)
