import subprocess
import sys

import numpy as np
import pytest

from petakit5d.core import arrays


def test_import_core_loads_nothing_heavy():
    code = (
        "import sys, petakit5d.core, petakit5d.core.paths, petakit5d.core.threads; "
        "print('numpy' in sys.modules)"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    ).stdout.strip()
    assert out == "False"


def test_lazy_names_resolve():
    import petakit5d.core as core

    assert core.cast is __import__("petakit5d.core.dtypes", fromlist=["cast"]).cast
    assert core.bbox_to_slices(core.full_bbox((3, 4, 5))) == (slice(0, 3), slice(0, 4), slice(0, 5))
    with pytest.raises(AttributeError):
        core.not_a_name  # noqa: B018


@pytest.mark.parametrize(
    "shape,size",
    [
        ((3, 4, 1), (3, 4)),
        ((3, 1, 1, 1), (3, 1)),
        ((1, 1, 5), (1, 1, 5)),
        ((5,), (1, 5)),
        ((), (1, 1)),
        ((0, 3), (0, 3)),
        ((2, 3, 4, 1, 1), (2, 3, 4)),
    ],
)
def test_msize_drops_trailing_singletons(shape, size):
    a = np.zeros(shape)
    assert arrays.msize(a) == size
    assert arrays.mndims(a) == len(size)
    assert arrays.msize(a, 3) == (size[2] if len(size) > 2 else 1)
    assert arrays.mnumel(a) == int(np.prod(size))


@pytest.mark.parametrize(
    "shape,squeezed",
    [
        ((1, 1, 5), (5, 1)),
        ((1, 5), (1, 5)),
        ((1, 3, 1, 4), (3, 4)),
        ((3, 4, 1), (3, 4)),
        ((1, 1, 1, 1, 7), (7, 1)),
        ((2, 1, 3), (2, 3)),
    ],
)
def test_msqueeze_as_matlab(shape, squeezed):
    a = np.asfortranarray(np.arange(np.prod(shape), dtype=float).reshape(shape, order="F"))
    s = arrays.msqueeze(a)
    assert s.shape == squeezed
    # column-major element order is kept
    np.testing.assert_array_equal(s.ravel(order="F"), a.ravel(order="F"))


def test_as_volume_views_and_copies():
    v = np.zeros((4, 5, 3), dtype=np.uint16, order="F")
    assert arrays.as_volume(v) is v
    c = np.zeros((4, 5, 3), dtype=np.uint16)  # C order
    f = arrays.as_volume(c)
    assert f.flags.f_contiguous and f.shape == (4, 5, 3)
    im = np.zeros((4, 5), dtype=np.float32, order="F")
    v2 = arrays.as_volume(im)
    assert v2.shape == (4, 5, 1) and np.shares_memory(v2, im)
    with pytest.raises(ValueError):
        arrays.as_volume(np.zeros(5))


def test_check_volume():
    v = arrays.fzeros((4, 5, 3), "uint16")
    assert arrays.check_volume(v, dtypes=("uint16", "single")) is v
    with pytest.raises(ValueError, match="Fortran"):
        arrays.check_volume(np.zeros((4, 5, 3), dtype=np.uint16))
    with pytest.raises(TypeError, match="class"):
        arrays.check_volume(v, dtypes=("single",))
    with pytest.raises(ValueError, match="dimensions"):
        arrays.check_volume(np.zeros(4), ndim=(3,))


def test_tiff_layout_round_trip():
    yxz = np.asfortranarray(np.arange(24).reshape(2, 3, 4))
    zyx = arrays.yxz_to_zyx(yxz)
    assert zyx.shape == (4, 2, 3) and zyx.flags.c_contiguous
    assert zyx[3, 1, 2] == yxz[1, 2, 3]
    back = arrays.zyx_to_yxz(zyx)
    assert back.flags.f_contiguous
    np.testing.assert_array_equal(back, yxz)


def test_linear_indices_are_matlab_minus_one():
    a = np.asfortranarray(np.arange(60).reshape(3, 4, 5))
    # MATLAB a(2, 3, 4) is linear index sub2ind([3 4 5], 2, 3, 4) = 2 + 3*2 + 12*3 = 44
    assert arrays.sub2ind((3, 4, 5), 2, 3, 4) == 44
    assert a.ravel(order="F")[44 - 1] == a[1, 2, 3]
    assert arrays.ind2sub((3, 4, 5), 44) == (2, 3, 4)
    s = arrays.ind2sub((3, 4, 5), np.array([1, 60]))
    np.testing.assert_array_equal(s[0], [1, 3])
    # fewer outputs fold the trailing dimensions
    assert arrays.ind2sub((3, 4, 5), 44, nout=2) == (2, 15)
    with pytest.raises(IndexError):
        arrays.sub2ind((3, 4, 5), 4, 1, 1)
    with pytest.raises(IndexError):
        arrays.ind2sub((3, 4, 5), 61)


def test_mpermute():
    a = np.asfortranarray(np.arange(24).reshape(2, 3, 4))
    p = arrays.mpermute(a, (3, 1, 2))
    assert p.shape == (4, 2, 3) and p.flags.f_contiguous
    assert p[3, 1, 2] == a[1, 2, 3]
    with pytest.raises(ValueError):
        arrays.mpermute(a, (1, 2))


@pytest.mark.parametrize("order", ["yxz", "xyz", "zxy", "yzx", "zyx", "xzy"])
def test_axis_order_mapping_stores_yxz_data_in_the_given_order(order):
    # AxisOrderMappingTest/storesYXZDataInTheGivenOrder
    vol = np.arange(1, 2 * 3 * 4 + 1).reshape((2, 3, 4), order="F")
    stored = arrays.mpermute(vol, arrays.axis_order_mapping(order, "yxz"))
    assert stored.shape == tuple(vol.shape["yxz".index(c)] for c in order)


def test_axis_order_mapping_checks():
    assert arrays.axis_order_mapping("XYZ") == (2, 1, 3)  # the input is lowered
    with pytest.raises(ValueError, match="Output"):
        arrays.axis_order_mapping("xyz", "XYZ")
    with pytest.raises(ValueError, match="Input"):
        arrays.axis_order_mapping("xy")


def test_class_names():
    assert arrays.numpy_dtype("single") == np.float32
    assert arrays.matlab_class(np.zeros(2, np.uint16)) == "uint16"
    assert arrays.matlab_class(1.5) == "double"
    assert arrays.matlab_class(True) == "logical"
    assert arrays.matlab_class(np.float32(1)) == "single"
    with pytest.raises(TypeError):
        arrays.numpy_dtype("complex128")


def test_data_type_to_byte_number():
    # MasterComputeFitTest uses it for the memory of a dataset
    assert [
        arrays.data_type_to_byte_number(t) for t in ("uint8", "uint16", "single", "double")
    ] == [1, 2, 4, 8]
    with pytest.raises(ValueError, match="Unsupported data type"):
        arrays.data_type_to_byte_number("int32")


def test_fortran_allocators():
    for f in (arrays.fzeros, arrays.fones, arrays.fempty):
        a = f((3, 4, 2), "uint16")
        assert a.flags.f_contiguous and a.dtype == np.uint16
    assert arrays.ffull((2, 2), 7, "single")[1, 1] == 7
