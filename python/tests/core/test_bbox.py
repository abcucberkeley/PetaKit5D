import numpy as np
import pytest

from petakit5d.core import bbox as bb


def test_bbox_to_slices_is_matlab_indexing():
    a = np.asfortranarray(np.arange(5 * 6 * 7).reshape(5, 6, 7))
    box = [2, 3, 4, 4, 6, 7]  # MATLAB a(2:4, 3:6, 4:7)
    sl = bb.bbox_to_slices(box)
    assert sl == (slice(1, 4), slice(2, 6), slice(3, 7))
    c = bb.crop(a, box)
    assert c.shape == bb.bbox_size(box) == (3, 4, 4)
    assert c[0, 0, 0] == a[1, 2, 3]
    assert np.shares_memory(c, a)
    assert bb.slices_to_bbox(sl) == tuple(box)


def test_full_bbox_and_2d():
    assert bb.full_bbox((5, 6, 7)) == (1, 1, 1, 5, 6, 7)
    assert bb.full_bbox((5, 6)) == (1, 1, 1, 5, 6, 1)
    im = np.zeros((5, 6))
    assert bb.crop(im, [1, 2, 1, 3, 4, 1]).shape == (3, 3)
    with pytest.raises(ValueError, match="2D"):
        bb.crop(im, [1, 2, 1, 3, 4, 2])
    assert bb.slices_to_bbox([slice(0, 2), slice(1, 3)]) == (1, 2, 1, 2, 3, 1)


def test_check_bbox():
    assert bb.check_bbox(np.array([1.0, 2, 3, 4, 5, 6])) == (1, 2, 3, 4, 5, 6)
    for bad in ([1, 2, 3], [1.5, 1, 1, 2, 2, 2], [3, 1, 1, 2, 2, 2]):
        with pytest.raises(ValueError):
            bb.check_bbox(bad)
    with pytest.raises(ValueError, match="outside"):
        bb.check_bbox([0, 1, 1, 2, 2, 2], (5, 5, 5))
    assert bb.check_bbox([0, 1, 1, 9, 2, 2], (5, 5, 5), allow_outside=True) == (0, 1, 1, 9, 2, 2)


def test_clamp_intersect_and_pad():
    assert bb.clamp_bbox([-3, 2, 1, 9, 4, 10], (5, 6, 7)) == (1, 2, 1, 5, 4, 7)
    with pytest.raises(ValueError, match="overlap"):
        bb.clamp_bbox([7, 1, 1, 9, 2, 2], (5, 6, 7))
    assert bb.bbox_intersect([1, 1, 1, 4, 4, 4], [3, 2, 4, 9, 9, 9]) == (3, 2, 4, 4, 4, 4)
    assert bb.bbox_intersect([1, 1, 1, 2, 2, 2], [3, 3, 3, 4, 4, 4]) is None
    box = [-1, 2, 1, 7, 6, 9]
    pads = bb.pad_widths(box, (5, 6, 7))
    assert pads == ((2, 2), (0, 0), (0, 2))
    a = np.ones((5, 6, 7))
    padded = np.pad(bb.crop(a, bb.clamp_bbox(box, a.shape)), pads)
    assert padded.shape == bb.bbox_size(box)


def test_slices_to_bbox_errors():
    with pytest.raises(ValueError):
        bb.slices_to_bbox([slice(0, 2, 2), slice(0, 1), slice(0, 1)])
    with pytest.raises(ValueError):
        bb.slices_to_bbox([slice(None), slice(0, 1), slice(0, 1)])
    assert bb.slices_to_bbox([slice(None), slice(1, None), slice(0, 1)], (4, 5, 6)) == (
        1,
        2,
        1,
        4,
        5,
        1,
    )
    with pytest.raises(ValueError, match="empty"):
        bb.slices_to_bbox([slice(2, 2), slice(0, 1), slice(0, 1)])
