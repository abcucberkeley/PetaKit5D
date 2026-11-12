"""petakit5d.io sizes, types and zarr metadata (task IO2).

Ports the intent of tests/io/ImageInfoTest.m (test names kept), plus the tiff header
reader and getImageBoundingBox.
"""

import struct
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "bindings"))
import iobind_cases as ic  # noqa: E402

pytest.importorskip("cpptiff")
pytest.importorskip("cppzarr")

from petakit5d.io import (  # noqa: E402
    create_zarr,
    get_image_bounding_box,
    get_image_dtype,
    get_image_size,
    get_image_size_batch,
    get_zarr_info,
    write_tiff,
    write_zarr,
)
from petakit5d.io import imagesize as S  # noqa: E402

DTYPES = ["uint8", "uint16", "single", "double"]


@pytest.mark.parametrize("dtype", DTYPES)
def test_tiffSizeAndType(tmp_path, dtype):
    nt = S.numpy_dtype(dtype)
    write_tiff(np.zeros((23, 31, 7), nt, order="F"), tmp_path / "im.tif")
    assert get_image_size(tmp_path / "im.tif") == (23, 31, 7)
    assert get_image_dtype(tmp_path / "im.tif") == dtype
    write_tiff(np.zeros((23, 31), nt), tmp_path / "im2d.tif")
    assert get_image_size(tmp_path / "im2d.tif") == (23, 31, 1)


@pytest.mark.parametrize("dtype", DTYPES)
def test_zarrSizeAndType(tmp_path, dtype):
    fn = tmp_path / "im.zarr"
    write_zarr(np.zeros((23, 31, 7), S.numpy_dtype(dtype), order="F"), fn, block_size=(16, 16, 4))
    assert get_image_size(fn) == (23, 31, 7)
    assert get_image_dtype(fn) == dtype
    info = get_zarr_info(fn)
    assert info.size == (23, 31, 7) and info.block_size == (16, 16, 4)
    assert info.dtype == dtype and info.order == "F" and info.zarr_format == 2
    assert info.fill_value == 0 and info.filters is None
    assert info.compressor == {
        "blocksize": 0,
        "clevel": 1,
        "cname": "zstd",
        "id": "blosc",
        "shuffle": 1,
    }


def test_zarrFolderWithoutExtension(tmp_path):
    fn = tmp_path / "im_zarr_folder"
    create_zarr(fn, data_size=(23, 31, 7), block_size=(16, 16, 4), dtype="single")
    assert get_image_size(fn) == (23, 31, 7)
    assert get_image_dtype(fn, True) == "single"
    with pytest.raises(ValueError):
        get_image_dtype(fn)


def test_batchSizesOfMixedFiles(tmp_path):
    fns = [tmp_path / "a.tif", tmp_path / "b.tif", tmp_path / "c.zarr"]
    write_tiff(np.zeros((10, 20, 3), np.uint16, order="F"), fns[0])
    write_tiff(np.zeros((11, 21, 4), np.uint16, order="F"), fns[1])
    write_zarr(np.zeros((12, 22, 5), np.uint16, order="F"), fns[2])
    out = get_image_size_batch(fns, False)
    assert out.dtype == np.int64
    assert out.tolist() == [[10, 20, 3], [11, 21, 4], [12, 22, 5]]


def test_batchSizesOfManyTiffs(tmp_path):
    fns = [tmp_path / f"im_{i:03d}.tif" for i in range(1, 102)]
    sz = [[8, 9, 2 + i % 3] for i in range(1, 102)]
    for fn, s in zip(fns, sz, strict=False):
        write_tiff(np.zeros(s, np.uint16, order="F"), fn)
    assert get_image_size_batch(fns).tolist() == sz
    assert get_image_size_batch(fns, False).tolist() == sz


def test_batchOfDifferentDimensionsErrors(tmp_path):
    write_tiff(np.zeros((4, 4, 2), np.uint16, order="F"), tmp_path / "a.tif")
    create_zarr(tmp_path / "b.zarr", data_size=(4, 4), expand_2d_dim=False)
    with pytest.raises(ValueError):
        get_image_size_batch([tmp_path / "a.tif", tmp_path / "b.zarr"])


def test_getZarrInfoOfMissingFileErrors(tmp_path):
    with pytest.raises(FileNotFoundError):
        get_zarr_info(tmp_path / "missing.zarr")
    with pytest.raises(FileNotFoundError):
        get_zarr_info(tmp_path / "missing")


def test_unknownFormatErrors(tmp_path):
    fn = tmp_path / "im.png"
    fn.write_bytes(b"\x89PNG\r\n\x1a\n")
    with pytest.raises(ValueError):
        get_image_dtype(fn)
    with pytest.raises(ValueError):
        get_image_size(fn)


def test_queryingMetadataKeepsWarningState(tmp_path):
    # MATLAB's queries turned all warnings off; Python's change no warning filters
    write_tiff(np.zeros((4, 4, 2), np.uint16, order="F"), tmp_path / "im.tif")
    write_zarr(np.zeros((4, 4, 2), np.uint16, order="F"), tmp_path / "im.zarr")
    before = list(warnings.filters)
    get_image_size(tmp_path / "im.tif")
    get_zarr_info(tmp_path / "im.zarr")
    assert warnings.filters == before


def test_unknownZarrTypeWarns(tmp_path):
    import cppzarr

    fn = tmp_path / "i.zarr"
    cppzarr.write_zarr(
        str(fn), np.ones((4, 4, 4), np.int16), order="F", chunks=[4, 4, 4], dimension_separator="."
    )
    with pytest.warns(UserWarning):
        assert get_image_dtype(fn) == ""
    assert get_zarr_info.__name__ and get_image_size(fn) == (4, 4, 4)


@pytest.mark.parametrize("dtype", ic.TIFF_DTYPES)
def test_tiffTypeOfEveryClass(tmp_path, dtype):
    fn = tmp_path / "im.tif"
    write_tiff(ic.make_data(dtype, (5, 6, 2), "random", 1), fn)
    assert get_image_dtype(fn) == ic.MATLAB_CLASS[dtype]


def _classic_tiff(fn, endian, bits, sample_format, desc=None, spp=1, compression=1):
    # one uncompressed image, 3 x 2, in classic TIFF of either byte order
    e = "<" if endian == "II" else ">"
    data = b"\0" * (6 * bits // 8 * spp)
    n = 10 + (desc is not None)
    after = 8 + 2 + 12 * n + 4
    data_off = after + (len(desc) if desc is not None else 0)
    entries = [
        (256, 3, 1, 2),
        (257, 3, 1, 3),
        (258, 3, 1, bits),
        (259, 3, 1, compression),
        (262, 3, 1, 2 if spp == 3 else 1),
        (273, 4, 1, data_off),
        (277, 3, 1, spp),
        (278, 3, 1, 3),
        (279, 4, 1, len(data)),
        (339, 3, 1, sample_format),
    ]
    if desc is not None:
        entries.append((270, 2, len(desc), after))
    entries.sort()
    with open(fn, "wb") as f:
        f.write(endian.encode() + struct.pack(e + "HI", 42, 8))
        f.write(struct.pack(e + "H", len(entries)))
        for tag, typ, cnt, val in entries:
            f.write(struct.pack(e + "HHI", tag, typ, cnt))
            f.write(struct.pack(e + "HH", val, 0) if typ == 3 else struct.pack(e + "I", val))
        f.write(struct.pack(e + "I", 0))
        if desc is not None:
            f.write(desc)
        f.write(data)


@pytest.mark.parametrize("endian", ["II", "MM"])
def test_tiffHeaderReader(tmp_path, endian):
    fn = tmp_path / "h.tif"
    _classic_tiff(fn, endian, 32, 3, b"ImageJ=1.54f\nimages=12\n\0")
    h = S._tiff_header(str(fn))
    assert (h.width, h.length, h.bits, h.sample_format, h.compression) == (2, 3, 32, 3, 1)
    assert S._imagej_images(h.description) == 12
    assert get_image_dtype(fn) == "single"
    # one uncompressed directory with ImageJ's images= is a stack of that many slices
    assert get_image_size(fn) == (3, 2, 12)
    # ... only when uncompressed (getImageSizeMex)
    _classic_tiff(fn, endian, 32, 3, b"ImageJ=1.54f\nimages=12\n\0", compression=5)
    assert get_image_size(fn) == (3, 2, 1)
    _classic_tiff(fn, endian, 16, 2)
    assert get_image_dtype(fn) == "int16" and get_image_size(fn) == (3, 2, 1)
    _classic_tiff(fn, endian, 8, 1, spp=3)
    assert get_image_size(fn) == (3, 2, 1, 3)


def test_imagejImagesParsing():
    f = S._imagej_images
    assert f(b"ImageJ=1\nimages=  42abc") == 42
    assert f(b"ImageJ=1\nimages=x") == 0
    assert f(b"images=5") == 0  # not an ImageJ description
    assert f(None) == 0


def test_bigTiffHeaderReader(tmp_path):
    fn = tmp_path / "b.tif"
    write_tiff(np.zeros((1031, 7, 3), np.int32, order="F"), fn)
    assert fn.read_bytes()[2:4] == b"\x2b\x00"  # BigTIFF
    h = S._tiff_header(str(fn))
    assert (h.width, h.length, h.bits, h.sample_format, h.compression) == (7, 1031, 32, 2, 5)


# --- getImageBoundingBox --------------------------------------------------------------------


def _bbox_reference(a):
    # find(sum(I, other axes) > 0, 1, 'first' / 'last'), in exact arithmetic
    a3 = a.reshape(a.shape + (1,) * (3 - a.ndim)).astype(object)
    nd = len(ic.matlab_shape(a.shape))
    lo, hi = [], []
    for ax in range(3 if nd == 3 else 2):
        s = np.moveaxis(a3, ax, 0).reshape(a3.shape[ax], -1).sum(axis=1)
        idx = [i for i, v in enumerate(s) if v > 0]
        if idx:
            lo.append(idx[0] + 1)
            hi.append(idx[-1] + 1)
    return tuple(lo + hi)


@pytest.mark.parametrize(
    "dtype", ["bool", "uint8", "uint16", "int16", "int32", "float32", "float64"]
)
def test_imageBoundingBox(dtype):
    a = np.zeros((20, 30, 10), dtype, order="F")
    a[3:9, 5:21, 2:7] = 1
    assert get_image_bounding_box(a) == (4, 6, 3, 9, 21, 7)
    assert get_image_bounding_box(a[:, :, 4]) == (4, 6, 9, 21)
    assert get_image_bounding_box(a[:, :, 4:5]) == (4, 6, 9, 21)  # MATLAB: a 2D image
    assert get_image_bounding_box(np.zeros((4, 4, 4), dtype)) == ()
    rng = np.random.default_rng(3)
    if np.dtype(dtype).kind in "iu":
        lo = -3 if np.dtype(dtype).kind == "i" else 0
        b = rng.integers(lo, 4, (13, 9, 6)).astype(dtype)
        b[(rng.random(b.shape) < 0.8)] = 0
        assert get_image_bounding_box(b) == _bbox_reference(b)


def test_imageBoundingBoxNanAndNegative():
    a = np.zeros((5, 6, 3), np.float32)
    a[1, 2, 1] = 2.0
    a[3, 4, 2] = np.nan
    # the slices with a NaN sum to NaN, which is not > 0
    assert get_image_bounding_box(a) == (2, 3, 2, 2, 3, 2)
    a[1, 4, 2] = np.nan
    a[2, 3, 1] = 1.0
    assert get_image_bounding_box(a) == (3, 3, 2, 3, 4, 2)
    a[2, 3, 1] = 0.0
    a[1, 2, 1] = -1.0
    assert get_image_bounding_box(a) == ()
