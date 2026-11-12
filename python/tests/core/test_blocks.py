import numpy as np
import pytest

from petakit5d.core import bbox as bb
from petakit5d.core import blocks


def _cover(boxes, shape):
    count = np.zeros(shape, dtype=np.int32)
    for b in boxes:
        count[b[0] : b[3], b[1] : b[4], b[2] : b[5]] += 1
    return count


@pytest.mark.parametrize("same", [True, False])
def test_zarr_chunks_tile_the_volume_once(same):
    vol = (41, 130, 77)  # (nz, ny, nx)
    b, r, loc = blocks.zarr_chunk_coordinates(
        vol,
        batch_size=(20, 50, 50),
        block_size=(8, 16, 16),
        border_size=(3, 5, 0),
        same_batch_size=same,
    )
    assert np.all(_cover(r, vol) == 1)
    # MATLAB's batch numbering: y fastest, then x, then z
    assert r[1, 1] > r[0, 1] and r[1, 0] == r[0, 0] and r[1, 2] == r[0, 2]
    assert r[0, 4] == 64  # 50 rounded up to whole blocks of 16
    for bi, ri, li in zip(b, r, loc, strict=True):
        assert np.all(bi[:3] + li[:3] == ri[:3])
        assert np.all(li[3:] - li[:3] == ri[3:] - ri[:3])
        assert np.all(bi[:3] >= 0) and np.all(bi[3:] <= vol)


def test_zarr_chunks_match_the_matlab_form():
    rng = np.random.default_rng(0)
    for _ in range(200):
        vol = tuple(int(v) for v in rng.integers(1, 300, 3))
        args = {
            "batch_size": tuple(int(v) for v in rng.integers(1, 120, 3)),
            "block_size": tuple(int(v) for v in rng.integers(1, 40, 3)),
            "border_size": tuple(int(v) for v in rng.integers(0, 9, 3)),
            "same_batch_size": bool(rng.integers(0, 2)),
        }
        p = blocks.zarr_chunk_coordinates(vol, **args)
        m = blocks.zarr_chunk_coordinates_matlab(
            [vol[1], vol[2], vol[0]],
            **{k: ([v[1], v[2], v[0]] if isinstance(v, tuple) else v) for k, v in args.items()},
        )
        for pp, mm in zip(p, m, strict=True):
            np.testing.assert_array_equal(bb.boxes_to_matlab(pp), mm)


def test_zarr_chunks_with_a_box():
    b, r, loc = blocks.zarr_chunk_coordinates(
        (100, 100, 100), bbox=(30, 10, 20, 40, 60, 70), border_size=0
    )
    np.testing.assert_array_equal(b[0], [30, 10, 20, 40, 60, 70])
    np.testing.assert_array_equal(r[0], [0, 0, 0, 10, 50, 50])  # relative to the box
    np.testing.assert_array_equal(loc[0], [0, 0, 0, 10, 50, 50])
    b2, _, _ = blocks.zarr_chunk_coordinates(
        (100, 100, 100), bbox=(30, 10, 20, None, 60, None), border_size=0
    )
    np.testing.assert_array_equal(b2[0], [30, 10, 20, 100, 60, 100])


def test_zarr_chunks_checks():
    with pytest.raises(ValueError):
        blocks.zarr_chunk_coordinates((10, 10))
    with pytest.raises(ValueError):
        blocks.zarr_chunk_coordinates((10, 10, 10), batch_size=(0, 1, 1))
    with pytest.raises(ValueError):
        blocks.zarr_chunk_coordinates((10, 10, 10), border_size=(1, 2))


def test_sub_volume_coordinates_cover_with_overlap():
    vol = (900, 2000, 1800)
    boxes = blocks.sub_volume_coordinates(vol, chunk_size=(500, 1000, 1000), overlap_size=50)
    assert boxes.shape[1] == 6
    # MATLAB's meshgrid order: x fastest, then y, then z
    assert boxes[1, 2] != boxes[0, 2] and boxes[1, 1] == boxes[0, 1]
    assert boxes[:, :3].min() == 0 and tuple(boxes[:, 3:].max(axis=0)) == vol
    with pytest.raises(ValueError, match="not valid"):
        blocks.sub_volume_coordinates((100, 1, 100))


def test_bbox_to_blocks():
    s = blocks.bbox_to_blocks(
        (2, 99, 199, 60, 700, 900), (0, 0, 0), (32, 256, 256), (128, 1024, 1024), (4, 8, 8)
    )
    assert s.b_sub_size == (4, 4, 4)
    np.testing.assert_array_equal(
        s.block_sub[:4], [[0, 0, 0], [0, 1, 0], [0, 2, 0], [0, 0, 1]]
    )  # y fastest
    np.testing.assert_array_equal(s.block_ind[:4], [0, 1, 2, 4])
    np.testing.assert_array_equal(s.b_coords[0], [2, 99, 199, 36, 264, 264])
    np.testing.assert_array_equal(s.w_coords[1], [2, 248, 199, 36, 520, 264])
    s0 = blocks.bbox_to_blocks(
        (2, 99, 199, 60, 700, 900), (0, 0, 0), (32, 256, 256), (128, 1024, 1024)
    )
    count = np.zeros((128, 1024, 1024), dtype=np.int8)
    for w in s0.w_coords:
        count[w[0] : w[3], w[1] : w[4], w[2] : w[5]] += 1
    assert count[2:60, 99:700, 199:900].min() == 1 and count.sum() == 58 * 601 * 701
    np.testing.assert_array_equal(s0.bbox_coords, s0.w_coords - np.tile(np.array([2, 99, 199]), 2))
    with pytest.raises(IndexError):
        blocks.bbox_to_blocks(
            (1, 2000, 1, 2, 2100, 2), (0, 0, 0), (32, 256, 256), (128, 1024, 1024)
        )


@pytest.mark.parametrize(
    "n,up,down,odd",
    [
        (1031, 1050, 1024, 1050),
        (1, 1, 1, 1),
        (11, 12, 10, 12),
        (49, 50, 48, 49),
        (7, 8, 6, 7),
        (14, 14, 14, 14),
    ],
)
def test_find_good_factor_number(n, up, down, odd):
    # 1029 = 3 * 7^3 is odd with 7 as its largest factor: skipped unless allowed
    assert blocks.find_good_factor_number(n, 1) == up
    assert blocks.find_good_factor_number(n, -1) == down
    assert blocks.find_good_factor_number(n, 1, True) == odd


def test_prime_factors():
    assert blocks.prime_factors(0) == [0] and blocks.prime_factors(1) == [1]
    assert blocks.prime_factors(360) == [2, 2, 2, 3, 3, 5]
    with pytest.raises(ValueError):
        blocks.prime_factors(-4)
