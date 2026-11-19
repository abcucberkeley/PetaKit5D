"""Tests the mutation campaign asked for (E06, docs/hardening.md): kernels."""

import numpy as np
import pytest

from petakit5d.kernels import stitch as K


def test_feather_blending_into_a_4d_block_is_out_of_bound():
    # MATLAB's code (double tiles): a block of more than 3 dimensions is refused
    # (the box fits the block's first three dimensions)
    rng = np.random.default_rng(2)
    f = rng.uniform(0, 100, size=(2, 2, 4, 5))
    d = rng.uniform(0, 1, size=(2, 2, 4, 5)).astype(np.float32)
    block = np.zeros((2, 4, 6, 7))
    with pytest.raises(K.StitchKernelError, match="out of bound") as e:
        K.feather_blending_3d(f, d, block, (0, 0, 0, 2, 4, 5))
    assert e.value.identifier == "feather_blending:inputError"
    # (a 3-D block takes it)
    out, mex = K.feather_blending_3d(f, d, np.zeros((2, 4, 6)), (0, 0, 0, 2, 4, 5))
    assert not mex and out.shape == (2, 4, 6)


@pytest.mark.parametrize("size_t", [(0, 2, 2), (2, 1.5, 2)])
def test_integral_image_of_matlab_classes_needs_whole_window_sizes(size_t):
    a = np.arange(4 * 6 * 5, dtype=np.float64).reshape(4, 6, 5)
    with pytest.raises(K.StitchKernelError, match="size_t must be 3 integers") as e:
        K.integral_image_3d(a, size_t)
    assert e.value.identifier == "integral_image_3d:inputError"
