"""Tests the mutation campaign asked for (E06, docs/hardening.md): matlabcompat."""

import math

import numpy as np
import pytest

from petakit5d.matlabcompat import label, morphology, stats, warp


def test_median_of_a_scalar_along_an_axis():
    assert stats.median(np.float64(5.0), axis=0) == 5.0
    assert stats.median(np.float32(2.5), axis=-1) == np.float32(2.5)


@pytest.mark.parametrize("method", [5, "", None])
def test_imwarp_method_must_be_a_nonempty_string(method):
    with pytest.raises(warp.ImwarpError) as e:
        warp._method(method)
    assert e.value.identifier == "MATLAB:imwarp:invalidType"


@pytest.mark.parametrize(
    ("r", "attribute"),
    [(-1, "Nonnegative"), (1.5, "Integer"), (math.inf, "Integer"), (math.nan, "Integer")],
)
def test_strel_sizes_must_be_nonnegative_integers(r, attribute):
    # validateattributes' identifiers (probed in MATLAB; the oracle cases strel.error.*): an
    # infinite size was Python's OverflowError
    with pytest.raises(morphology.MatlabError, match="input number 2, R,") as e:
        morphology.strel("disk", r)
    assert e.value.identifier == f"MATLAB:strel:expected{attribute}"


def test_conndef_type_must_be_minimal_or_maximal():
    with pytest.raises(label.MatlabError) as e:
        label.conndef(3, "medium")
    assert e.value.identifier == "MATLAB:conndef:unrecognizedStringChoice"
    assert label.conndef(2, "max").sum() == 9 and label.conndef(2, "min").sum() == 5


@pytest.mark.parametrize("size", [(0, 3), (3, -1)])
def test_medfilt2_window_must_be_positive(size):
    with pytest.raises(morphology.MatlabError) as e:
        morphology.medfilt2(np.ones((5, 6)), size)
    assert e.value.identifier == "MATLAB:medfilt2:expectedPositive"


@pytest.mark.parametrize("dtype", [np.complex128, np.complex64])
def test_ifftn_is_real_only_for_exactly_conjugate_symmetric_input(dtype):
    # (the extension checks a row that is its own mirror up to its middle: the Nyquist element of
    # a 12 x 12 array lies past the first third of such a row along either axis)
    from petakit5d.matlabcompat import fft

    rng = np.random.default_rng(12)
    z = rng.standard_normal((12, 12)) + 1j * rng.standard_normal((12, 12))
    # (z(k) + conj(z(-k))) / 2 is conjugate symmetric bit for bit
    h = ((z + np.conj(np.roll(z[::-1, ::-1], 1, axis=(0, 1)))) / 2).astype(dtype)
    assert not np.iscomplexobj(fft.ifftn(h))
    h[6, 6] += 0.5j
    assert np.iscomplexobj(fft.ifftn(h))


def test_fdlibm_pow_kernel_of_negative_bases_with_whole_exponents():
    # (fdlibm's test of an odd or even exponent, in both of its branches: exponents under and
    # over 2^21)
    from petakit5d._ext import available, load
    from petakit5d.matlabcompat import elementary

    if not available("resize"):
        pytest.skip("the resize extension is not built")
    x = np.array([-2.0, -3.0, -1.5, -0.75, -7.0])
    for y in (3.0, 5.0, 6.0, 7.0, 10.0, 2.0**21 + 1, 2.0**22):
        out = np.empty_like(x)
        load("resize").fdlibm_pow(x, y, out)
        ref = np.asarray(elementary.power(x, y))
        assert np.array_equal(out.view(np.uint64), ref.view(np.uint64)), y


@pytest.mark.parametrize("shape", [(2, 4, 4), (4, 2, 4), (4, 4, 2)])
def test_cubic_warp_kernel_refuses_fewer_than_3_samples_along_an_axis(shape):
    from petakit5d._ext import available, load

    if not available("warp"):
        pytest.skip("the warp extension is not built")
    vol, out = np.zeros(shape, np.float32), np.zeros((3, 3, 3), np.float32)
    with pytest.raises(ValueError, match="cubic needs at least 3 samples"):
        load("warp").warp_affine(vol, np.eye(4), out, 2, 0.0, 1)


def test_only_rows_of_opposite_pairs_are_difference_rows():
    import scipy.sparse as sp

    from petakit5d.matlabcompat import lsq

    ok = lsq._difference_rows(sp.csr_matrix(np.array([[1.0, -1.0, 0.0], [0.0, 0.0, 0.0]])))
    assert ok is not None
    ends, vals = ok
    np.testing.assert_array_equal(ends, [[1, 0], [-1, -1]])
    np.testing.assert_array_equal(vals, [1.0, 0.0])
    bad = sp.csr_matrix(np.array([[1.0, -1.0, 0.0], [0.0, 1.0, 2.0]]))
    assert lsq._difference_rows(bad) is None
