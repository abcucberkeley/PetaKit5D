import numpy as np
import pytest

from petakit5d.core import axes
from petakit5d.core.arrays import from_matlab


def test_per_axis_vectors():
    assert axes.from_yxz([10, 20, 30]) == [30, 10, 20]  # [y x z] -> (z, y, x)
    assert axes.to_yxz((30, 10, 20)) == (10, 20, 30)
    assert axes.from_xyz((1, 2, 3)) == (3, 2, 1)
    assert axes.to_xyz((3, 2, 1)) == (1, 2, 3)
    assert axes.from_yxz(5) == 5 and axes.to_xyz(2.5) == 2.5
    assert axes.from_yxz((4, 5)) == (4, 5)  # 2D [y x] is (y, x)
    assert axes.from_xyz((4, 5)) == (5, 4)
    np.testing.assert_array_equal(axes.from_yxz(np.array([1, 2, 3])), [3, 1, 2])
    with pytest.raises(ValueError):
        axes.from_yxz([1, 2, 3, 4])


def test_per_axis_vectors_follow_the_array_mapping():
    m = np.zeros((4, 5, 6), order="F")  # MATLAB size [4 5 6] = [y x z]
    assert from_matlab(m).shape == tuple(axes.from_yxz([4, 5, 6]))


def test_inverses_on_generated_values():
    rng = np.random.default_rng(1)
    for _ in range(10000):
        n = int(rng.integers(2, 4))
        v = tuple(int(x) for x in rng.integers(-100, 100, size=n))
        assert axes.to_yxz(axes.from_yxz(v)) == v == axes.from_yxz(axes.to_yxz(v))
        assert axes.to_xyz(axes.from_xyz(v)) == v == axes.from_xyz(axes.to_xyz(v))


@pytest.mark.parametrize("ndim", [2, 3, 4, 5])
def test_dims_and_axes(ndim):
    assert axes.dim_to_axis(1) == -2 and axes.dim_to_axis(2) == -1 and axes.dim_to_axis(3) == -3
    m = np.zeros(tuple(range(2, 2 + ndim)), order="F")
    p = from_matlab(m)
    for dim in range(1, ndim + 1):
        ax = axes.dim_to_axis(dim, ndim)
        assert p.shape[ax] == m.shape[dim - 1]  # size(a, dim) is the axis' length
        assert axes.axis_to_dim(ax, ndim) == dim
        assert axes.axis_to_dim(ax - ndim, ndim) == dim
    with pytest.raises(ValueError):
        axes.dim_to_axis(ndim + 1, ndim)
    with pytest.raises(ValueError):
        axes.dim_to_axis(0)


def test_a_reduction_over_a_matlab_dim_is_over_its_axis():
    rng = np.random.default_rng(2)
    m = rng.random((3, 4, 5))
    p = from_matlab(m)
    for dim in (1, 2, 3):
        mx = m.max(axis=dim - 1, keepdims=True)  # MATLAB max(a, [], dim)
        np.testing.assert_array_equal(
            from_matlab(mx), p.max(axis=axes.dim_to_axis(dim, 3), keepdims=True)
        )


def test_axis_strings():
    assert axes.axis_string_from_matlab("y,x,z") == "z,y,x"
    assert axes.axis_string_from_matlab("yxz") == "zyx"
    assert axes.axis_string_from_matlab("x,y,z") == "z,x,y"
    for s in ("y,x,z", "x,y,z", "zxy", "xzy"):
        assert axes.axis_string_to_matlab(axes.axis_string_from_matlab(s)) == s
    with pytest.raises(ValueError):
        axes.axis_string_from_matlab("xy")


def test_axis_permutation_matches_matlab_axis_order_mapping():
    from petakit5d.core.arrays import axis_order_mapping, mpermute

    m = np.arange(2 * 3 * 4).reshape((2, 3, 4), order="F")  # MATLAB's yxz data
    p = from_matlab(m)
    for morder in ("yxz", "xyz", "zxy", "yzx", "zyx", "xzy"):
        stored_m = mpermute(m, axis_order_mapping(morder, "yxz"))  # MATLAB stores it in morder
        porder = axes.axis_string_from_matlab(morder)
        stored_p = np.transpose(p, axes.axis_permutation("zyx", porder))
        np.testing.assert_array_equal(from_matlab(stored_m), stored_p)
        back = np.transpose(stored_p, axes.axis_permutation(porder, "zyx"))
        np.testing.assert_array_equal(back, p)
    with pytest.raises(ValueError):
        axes.axis_permutation("zyy", "zyx")
