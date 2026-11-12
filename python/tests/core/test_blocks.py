import numpy as np
import pytest

from petakit5d.core import blocks


def _cover(region, size):
    """Each voxel of the box [1, size] in exactly one region."""
    count = np.zeros(size, dtype=np.int32)
    for r in region:
        count[r[0] - 1 : r[3], r[1] - 1 : r[4], r[2] - 1 : r[5]] += 1
    return count


@pytest.mark.parametrize("same", [True, False])
def test_zarr_chunks_tile_the_volume_once(same):
    vol = (130, 77, 41)
    b, r, l = blocks.zarr_chunk_coordinates(
        vol,
        batch_size=(50, 50, 20),
        block_size=(16, 16, 8),
        border_size=(5, 0, 3),
        same_batch_size=same,
    )
    assert np.all(_cover(r, vol) == 1)
    # column-major batch order: y fastest
    assert r[1, 0] > r[0, 0] and r[1, 1] == r[0, 1]
    # whole blocks per batch: 50 -> 64 (4 blocks of 16)
    assert r[0, 3] == 64
    # local boxes map each output region into its batch
    for bi, ri, li in zip(b, r, l, strict=True):
        assert np.all(bi[:3] + li[:3] - 1 == ri[:3])
        assert np.all(li[3:] - li[:3] == ri[3:] - ri[:3])
        assert np.all(bi[:3] >= 1) and np.all(bi[3:] <= vol)
    if same:
        sizes = {tuple(x) for x in (np.minimum(b[:, 3:], vol) - b[:, :3] + 1)}
        assert len(sizes) > 1  # the borders are cut at the volume edges
        inner = b[(b[:, 0] > 5) & (b[:, 3] < vol[0]) & (b[:, 2] > 3) & (b[:, 5] < vol[2])]
        if inner.size:
            assert len({tuple(x) for x in inner[:, 3:] - inner[:, :3]}) == 1


def test_zarr_chunks_with_bbox_and_sizes_smaller_than_blocks():
    b, r, l = blocks.zarr_chunk_coordinates(
        (100, 100, 100),
        batch_size=(1024, 1024, 1024),
        block_size=(256, 256, 256),
        bbox=(11, 21, 31, 60, 70, 40),
        border_size=0,
    )
    assert b.shape == (1, 6)
    np.testing.assert_array_equal(b[0], [11, 21, 31, 60, 70, 40])
    np.testing.assert_array_equal(r[0], [1, 1, 1, 50, 50, 10])  # relative to the box
    np.testing.assert_array_equal(l[0], [1, 1, 1, 50, 50, 10])


def test_zarr_chunks_checks():
    with pytest.raises(ValueError):
        blocks.zarr_chunk_coordinates((10, 10))
    with pytest.raises(ValueError):
        blocks.zarr_chunk_coordinates((10, 10, 10), batch_size=(0, 1, 1))
    with pytest.raises(ValueError):
        blocks.zarr_chunk_coordinates((10, 10, 10), border_size=(1, 2))


def test_sub_volume_coordinates_cover_with_overlap():
    vol = (2000, 1800, 900)
    xmin, xmax, ymin, ymax, zmin, zmax, nn = blocks.sub_volume_coordinates(
        vol, chunk_size=(1000, 1000, 500), overlap_size=50
    )
    assert nn == len(xmin) == len(zmax)
    # x varies fastest (meshgrid order)
    assert xmin[1] != xmin[0] and ymin[1] == ymin[0]
    assert ymin.min() == 1 and ymax.max() == vol[0] and zmax.max() == vol[2]
    for lo, hi in ((ymin, ymax), (xmin, xmax), (zmin, zmax)):
        starts = sorted(set(lo.tolist()))
        ends = sorted(set(hi.tolist()))
        assert all(s <= e + 1 for s, e in zip(starts[1:], ends[:-1], strict=True))
    for c in (ymax - ymin + 1, xmax - xmin + 1):
        assert all(max(blocks.prime_factors(int(v))) <= 7 or v == c.max() for v in c)


def test_sub_volume_coordinates_errors():
    with pytest.raises(ValueError, match="not valid"):
        blocks.sub_volume_coordinates((1, 100, 100))


def test_bbox_to_blocks():
    s = blocks.bbox_to_blocks(
        [100, 200, 3, 700, 900, 60], [1, 1, 1], [256, 256, 32], [1024, 1024, 128], [8, 8, 4]
    )
    assert s.b_sub_size == (4, 4, 4)
    np.testing.assert_array_equal(s.block_sub[:4], [[1, 1, 1], [2, 1, 1], [3, 1, 1], [1, 2, 1]])
    np.testing.assert_array_equal(s.block_ind[:4], [1, 2, 3, 5])
    np.testing.assert_array_equal(s.b_coords[0], [100, 200, 3, 264, 264, 36])
    np.testing.assert_array_equal(s.w_coords[1], [249, 200, 3, 520, 264, 36])
    # without borders the world coordinates tile the box exactly
    s0 = blocks.bbox_to_blocks(
        [100, 200, 3, 700, 900, 60], [1, 1, 1], [256, 256, 32], [1024, 1024, 128]
    )
    count = np.zeros((1024, 1024, 128), dtype=np.int8)
    for w in s0.w_coords:
        count[w[0] - 1 : w[3], w[1] - 1 : w[4], w[2] - 1 : w[5]] += 1
    assert count[99:700, 199:900, 2:60].min() == 1 and count.sum() == 601 * 701 * 58
    np.testing.assert_array_equal(s0.bbox_coords, s0.w_coords - np.tile(np.array([99, 199, 2]), 2))


def test_bbox_to_blocks_box_outside_raises():
    with pytest.raises(IndexError):
        blocks.bbox_to_blocks(
            [2000, 1, 1, 2100, 2, 2], [1, 1, 1], [256, 256, 32], [1024, 1024, 128]
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
    assert blocks.prime_factors(1031) == [1031]
    with pytest.raises(ValueError):
        blocks.prime_factors(-4)
