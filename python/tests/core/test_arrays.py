import subprocess
import sys

import numpy as np
import pytest

from petakit5d.core import arrays

CLASSES = (
    "double",
    "single",
    "int8",
    "uint8",
    "int16",
    "uint16",
    "int32",
    "uint32",
    "int64",
    "uint64",
    "logical",
)


def _rand(shape, cls, rng):
    a = rng.integers(0, 2**63 - 1, size=shape, dtype=np.int64)
    return (a & 1).astype(bool) if cls == "logical" else a.astype(arrays.numpy_dtype(cls))


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
    assert core.box_slices(core.full_box((3, 4, 5))) == (slice(0, 3), slice(0, 4), slice(0, 5))
    assert core.from_yxz((1, 2, 3)) == (3, 1, 2)
    with pytest.raises(AttributeError):
        core.not_a_name  # noqa: B018
    for gone in ("fzeros", "fempty", "bbox_to_slices", "full_bbox", "zyx_to_yxz"):
        assert not hasattr(core, gone)


@pytest.mark.parametrize("ndim", [2, 3, 4, 5])
@pytest.mark.parametrize("cls", CLASSES)
def test_round_trip_and_layout(ndim, cls):
    rng = np.random.default_rng(ndim)
    shape = tuple(int(v) for v in rng.integers(1, 6, size=ndim))
    m = np.asfortranarray(_rand(shape, cls, rng))
    p = arrays.from_matlab(m)
    assert p.flags.c_contiguous and p.dtype == m.dtype
    assert p.shape == arrays.size_from_matlab(shape)
    back = arrays.to_matlab(p)
    assert back.flags.f_contiguous and back.shape == shape
    np.testing.assert_array_equal(back, m)
    # element (z, y, x) of ours is MATLAB's (y, x, z): a(y, x, z, d4, ...) = p[..., d4, z, y, x]
    idx = tuple(int(rng.integers(0, n)) for n in shape)
    assert p[tuple(reversed(idx[2:])) + idx[:2]] == m[idx]


def test_3d_is_not_the_same_bytes():
    m = np.asfortranarray(np.arange(24, dtype=np.uint16).reshape((2, 3, 4), order="F"))
    p = arrays.from_matlab(m)
    assert p.shape == (4, 2, 3)
    np.testing.assert_array_equal(p, np.ascontiguousarray(np.transpose(m, (2, 0, 1))))
    assert p.tobytes() != m.tobytes(order="F")  # every z plane is transposed
    assert arrays.from_matlab(np.zeros((4, 5))).shape == (4, 5)  # an image keeps its shape
    assert arrays.from_matlab(np.arange(3)).shape == (3,)  # a row stays


@pytest.mark.parametrize("seed", range(20))
def test_matlab_order_iterates_as_matlab_colon(seed):
    rng = np.random.default_rng(seed)
    shape = tuple(int(v) for v in rng.integers(1, 7, size=int(rng.integers(2, 6))))
    m = np.asfortranarray(rng.random(shape))
    p = arrays.from_matlab(m)
    np.testing.assert_array_equal(arrays.matlab_order(p).ravel(), m.ravel(order="F"))
    # the linear index helpers agree with it
    k = int(rng.integers(0, m.size))
    idx = arrays.matlab_index_to_zyx(k, p.shape)
    assert p[idx] == m.ravel(order="F")[k]
    assert arrays.zyx_to_matlab_index(idx, p.shape) == k


def test_linear_index_formula():
    nz, ny, nx = 3, 4, 5
    z, y, x = 2, 1, 3
    assert arrays.zyx_to_matlab_index((z, y, x), (nz, ny, nx)) == y + ny * (x + nx * z)
    zz, yy, xx = np.meshgrid(range(nz), range(ny), range(nx), indexing="ij")
    k = arrays.zyx_to_matlab_index((zz, yy, xx), (nz, ny, nx))
    np.testing.assert_array_equal(np.sort(k.ravel()), np.arange(nz * ny * nx))
    with pytest.raises(IndexError):
        arrays.zyx_to_matlab_index((3, 0, 0), (nz, ny, nx))


def test_threaded_copy_is_bit_exact(monkeypatch):
    monkeypatch.setattr(arrays, "_THREADED_BYTES", 1024)
    rng = np.random.default_rng(3)
    m = np.asfortranarray(rng.random((37, 23, 19)).astype(np.float32))
    p = arrays.from_matlab(m, num_threads=4)
    np.testing.assert_array_equal(p, np.transpose(m, (2, 0, 1)))
    np.testing.assert_array_equal(arrays.to_matlab(p, num_threads=4), m)


@pytest.mark.parametrize(
    "shape,msize",
    [
        ((3, 4, 5), (4, 5, 3)),
        ((1, 4, 5), (4, 5)),
        ((4, 5), (4, 5)),
        ((2, 3, 4, 5), (4, 5, 3, 2)),
        ((1, 1, 4, 5), (4, 5)),
        ((2, 1, 4, 5), (4, 5, 1, 2)),
        ((7,), (1, 7)),
        ((), (1, 1)),
    ],
)
def test_sizes(shape, msize):
    assert arrays.size_to_matlab(shape) == msize


def test_size_from_matlab():
    assert arrays.size_from_matlab((4, 5, 3)) == (3, 4, 5)
    assert arrays.size_from_matlab((4, 5)) == (4, 5)
    assert arrays.size_from_matlab((4, 5), 3) == (1, 4, 5)
    assert arrays.size_from_matlab((4, 5, 3, 2)) == (2, 3, 4, 5)
    assert arrays.size_from_matlab((4, 5, 3, 1), 3) == (3, 4, 5)
    with pytest.raises(ValueError):
        arrays.size_from_matlab((4, 5, 3, 2), 3)


def test_as_volume_and_check_volume():
    v = np.zeros((3, 4, 5), dtype=np.uint16)
    assert arrays.as_volume(v) is v
    im = np.zeros((4, 5), dtype=np.float32)
    v2 = arrays.as_volume(im)
    assert v2.shape == (1, 4, 5) and np.shares_memory(v2, im)
    f = arrays.as_volume(np.asfortranarray(v))
    assert f.flags.c_contiguous
    with pytest.raises(ValueError):
        arrays.as_volume(np.zeros(5))
    assert arrays.check_volume(v, dtypes=("uint16", "single")) is v
    with pytest.raises(ValueError, match="from_matlab"):
        arrays.check_volume(np.asfortranarray(v))
    with pytest.raises(TypeError, match="class"):
        arrays.check_volume(v, dtypes=("single",))
    with pytest.raises(ValueError, match="dimensions"):
        arrays.check_volume(np.zeros(4), ndim=(3,))


# MATLAB's rules on MATLAB-form values (parameters, .mat contents)


@pytest.mark.parametrize(
    "shape,size",
    [
        ((3, 4, 1), (3, 4)),
        ((3, 1, 1, 1), (3, 1)),
        ((1, 1, 5), (1, 1, 5)),
        ((5,), (1, 5)),
        ((), (1, 1)),
    ],
)
def test_msize_drops_trailing_singletons(shape, size):
    a = np.zeros(shape)
    assert arrays.msize(a) == size
    assert arrays.mndims(a) == len(size)
    assert arrays.mnumel(a) == int(np.prod(size))


@pytest.mark.parametrize(
    "shape,squeezed",
    [((1, 1, 5), (5, 1)), ((1, 5), (1, 5)), ((1, 3, 1, 4), (3, 4)), ((2, 1, 3), (2, 3))],
)
def test_msqueeze_as_matlab(shape, squeezed):
    a = np.arange(np.prod(shape), dtype=float).reshape(shape, order="F")
    s = arrays.msqueeze(a)
    assert s.shape == squeezed
    np.testing.assert_array_equal(s.ravel(order="F"), a.ravel(order="F"))


def test_sub2ind_ind2sub_are_matlabs():
    assert arrays.sub2ind((3, 4, 5), 2, 3, 4) == 44
    assert arrays.ind2sub((3, 4, 5), 44) == (2, 3, 4)
    assert arrays.ind2sub((3, 4, 5), 44, nout=2) == (2, 15)
    with pytest.raises(IndexError):
        arrays.sub2ind((3, 4, 5), 4, 1, 1)


def test_mpermute_and_axis_order_mapping_on_matlab_values():
    a = np.arange(24).reshape((2, 3, 4), order="F")
    p = arrays.mpermute(a, (3, 1, 2))
    assert p.shape == (4, 2, 3)
    np.testing.assert_array_equal(p, arrays.from_matlab(a))  # permute(a, [3 1 2]) is from_matlab
    for order in ("yxz", "xyz", "zxy", "yzx", "zyx", "xzy"):
        stored = arrays.mpermute(a, arrays.axis_order_mapping(order, "yxz"))
        assert stored.shape == tuple(a.shape["yxz".index(c)] for c in order)
    with pytest.raises(ValueError, match="Output"):
        arrays.axis_order_mapping("xyz", "XYZ")


def test_class_names():
    assert arrays.numpy_dtype("single") == np.float32
    assert arrays.matlab_class(np.zeros(2, np.uint16)) == "uint16"
    assert arrays.matlab_class(1.5) == "double" and arrays.matlab_class(True) == "logical"
    with pytest.raises(TypeError):
        arrays.numpy_dtype("complex128")
    assert [
        arrays.data_type_to_byte_number(t) for t in ("uint8", "uint16", "single", "double")
    ] == [1, 2, 4, 8]
