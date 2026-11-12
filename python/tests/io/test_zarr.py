"""petakit5d.io zarr reading, writing and creation (task IO2).

Ports the intent of tests/io/ZarrIOTest.m (test names kept), plus the .zarray text of
createzarr, MATLAB's conversions in region writes (gap Z3) and the options of writezarr.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "bindings"))
import iobind_cases as ic  # noqa: E402

cppzarr = pytest.importorskip("cppzarr")

from petakit5d.io import (  # noqa: E402
    create_zarr,
    get_image_size,
    read_zarr,
    write_zarr,
    write_zarr_block,
)
from petakit5d.io import zarr as Z  # noqa: E402

IM = (37, 41, 13)

# .zarray as createzarr/writezarr (cpp-zarr mex, PetaKit5D at 5df18b2) write it for a
# 37 x 41 x 13 uint16 array with blockSize [16, 16, 5] (IO1)
MATLAB_ZARRAY = (
    '{\n    "chunks": [\n        16,\n        16,\n        5\n    ],\n    "compressor": {\n'
    '        "blocksize": 0,\n        "clevel": 1,\n        "cname": "zstd",\n'
    '        "id": "blosc",\n        "shuffle": 1\n    },\n    "dtype": "<u2",\n'
    '    "fill_value": 0,\n    "filters": null,\n    "order": "F",\n    "shape": [\n'
    '        37,\n        41,\n        13\n    ],\n    "zarr_format": 2\n}\n'
)


def _rand(dtype, shape, seed=1):
    return ic.make_data(dtype, shape, "random", seed)


def _chunks(fn):
    out = []
    for root, _, files in os.walk(fn):
        out += [
            os.path.relpath(os.path.join(root, f), fn).replace(os.sep, "/")
            for f in files
            if not f.startswith(".")
        ]
    return sorted(out)


def _meta(fn):
    return json.loads((Path(fn) / ".zarray").read_text())


# --- ZarrIOTest ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "block",
    [(16, 16, 5), (37, 41, 13), (64, 64, 64)],
    ids=["nonDivisible", "singleChunk", "larger"],
)
@pytest.mark.parametrize("dtype", ic.ZARR_WRAPPER_DTYPES)
def test_roundTripPreservesTypeAndValues(tmp_path, dtype, block):
    im = _rand(dtype, IM, 2)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, block_size=block)
    out = read_zarr(fn)
    assert out.flags.f_contiguous
    assert ic.bits_equal(out, im)


@pytest.mark.parametrize("compressor", ic.ZARR_COMPRESSORS)
def test_compressorsRoundTrip(tmp_path, compressor):
    im = _rand("uint16", IM, 3)
    fn = tmp_path / "im.zarr"
    create_zarr(
        fn, data_size=im.shape, block_size=(16, 16, 5), dtype="uint16", compressor=compressor
    )
    write_zarr(im, fn)
    assert ic.bits_equal(read_zarr(fn), im)
    assert _meta(fn)["compressor"]["id" if compressor == "gzip" else "cname"] == compressor


@pytest.mark.parametrize("chunks", [(16, 8, 5), (8, 16, 13)], ids=["zSmaller", "zLarger"])
@pytest.mark.parametrize("order", ["F", "C"])
def test_chunkOrderRoundTrip(tmp_path, order, chunks):
    im = _rand("uint16", IM, 4)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=im.shape, block_size=chunks, dtype="uint16", order=order)
    write_zarr(im, fn)
    out = read_zarr(fn)
    assert out.flags.f_contiguous and ic.bits_equal(out, im)
    assert _meta(fn)["order"] == order


def test_defaultMetadata(tmp_path):
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("uint16", IM, 5), fn, block_size=(16, 16, 5))
    assert (fn / ".zarray").read_text() == MATLAB_ZARRAY
    assert (fn / "2.2.2").is_file() and not (fn / "3.0.0").exists()
    assert sorted(os.listdir(fn))[0] == ".zarray" and len(os.listdir(fn)) == 1 + 27


def test_blockSizeIsCappedAtImageSize(tmp_path):
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("uint16", IM, 6), fn, block_size=(64, 64, 64))
    assert _meta(fn)["chunks"] == list(IM)


def test_nestedDimensionSeparator(tmp_path):
    im = _rand("uint16", IM, 7)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, block_size=(16, 16, 5), dim_separator="/")
    assert ic.bits_equal(read_zarr(fn), im)
    assert (fn / "2" / "2" / "2").is_file()
    assert _meta(fn)["dimension_separator"] == "/"


def test_shardedRoundTrip(tmp_path):
    # createzarr's v2 .zarray with a sharding_indexed codec; cpp-zarr writes regions into it
    im = _rand("uint16", IM, 8)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, block_size=(32, 32, 10), shard_size=(16, 16, 5))
    assert ic.bits_equal(read_zarr(fn), im)
    meta = _meta(fn)
    assert meta["zarr_format"] == 2 and meta["chunks"] == [32, 32, 10]
    assert meta["codecs"][0]["name"] == "sharding_indexed"
    assert meta["codecs"][0]["configuration"]["chunk_shape"] == [16, 16, 5]
    assert len(_chunks(fn)) == 8


def test_image2dIsStoredWithSingletonZ(tmp_path):
    im = _rand("uint16", (21, 17), 9)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn)
    out = read_zarr(fn)
    assert out.shape == (21, 17) and ic.bits_equal(out, im)
    assert get_image_size(fn) == (21, 17, 1)
    assert _meta(fn)["chunks"] == [21, 17, 1]


def test_readBoundingBox(tmp_path):
    im = _rand("uint16", IM, 10)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, block_size=(16, 16, 5))
    for b in ([3, 5, 2, 30, 33, 11], [17, 17, 6, 17, 17, 6], [1, 1, 7, 37, 41, 7], [1, 1, 1, *IM]):
        want = ic.matlab_shape(im[b[0] - 1 : b[3], b[1] - 1 : b[4], b[2] - 1 : b[5]].shape)
        out = read_zarr(fn, input_bbox=b)
        assert out.shape == want, b
        assert ic.bits_equal_m(out, im[b[0] - 1 : b[3], b[1] - 1 : b[4], b[2] - 1 : b[5]]), b


@pytest.mark.parametrize(
    "bbox", [[1, 1, 1, 38, 41, 13], [0, 1, 1, 5, 5, 5], [5, 5, 5, 4, 5, 5], [1, 1, 1, 5, 5]]
)
def test_readBoundingBoxPastTheEndErrors(tmp_path, bbox):
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("uint16", IM, 11), fn)
    with pytest.raises(ValueError):
        read_zarr(fn, input_bbox=bbox)


def test_readMissingFileErrors(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_zarr(tmp_path / "missing.zarr")


def test_partialWritesUpdateOnlyTheRegion(tmp_path):
    im = _rand("uint16", IM, 12)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=IM, block_size=(16, 16, 5), dtype="uint16")
    expected = np.zeros(IM, np.uint16, order="F")
    assert ic.bits_equal(read_zarr(fn), expected)
    b = [5, 7, 2, 30, 20, 9]
    write_zarr(im[4:30, 6:20, 1:9], fn, bbox=b, create=False)
    expected[4:30, 6:20, 1:9] = im[4:30, 6:20, 1:9]
    assert ic.bits_equal(read_zarr(fn), expected)
    b = [10, 1, 1, 37, 15, 13]
    region = _rand("uint16", (28, 15, 13), 13)
    write_zarr(region, fn, bbox=b, create=False)
    expected[9:37, 0:15, 0:13] = region
    assert ic.bits_equal(read_zarr(fn), expected)
    assert all(
        len(n.split(".")) == 3 and all(p.isdigit() for p in n.split(".")) for n in _chunks(fn)
    )
    assert (fn / ".zarray").read_text() == MATLAB_ZARRAY


def test_overwriteReplacesExistingArray(tmp_path):
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("float32", (20, 20, 20), 14), fn)
    im = _rand("uint16", IM, 15)
    write_zarr(im, fn, overwrite=True)
    assert ic.bits_equal(read_zarr(fn), im)
    assert _meta(fn)["dtype"] == "<u2"


def test_writeIntoExistingArrayKeepsItsMetadata(tmp_path):
    # without overwrite, the data goes into the array there, converted to its type
    fn = tmp_path / "im.zarr"
    create_zarr(
        fn, data_size=IM, block_size=(8, 16, 13), dtype="single", order="C", compressor="gzip"
    )
    text = (fn / ".zarray").read_text()
    im = _rand("uint16", IM, 16)
    write_zarr(im, fn, block_size=(16, 16, 5))
    assert (fn / ".zarray").read_text() == text
    out = read_zarr(fn)
    assert out.dtype == np.float32 and np.array_equal(out, im.astype(np.float32))


def test_sparseWriteSkipsEmptyChunks(tmp_path):
    im = np.zeros((32, 32, 10), np.uint16, order="F")
    im[:16, :16, :5] = 7
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, block_size=(16, 16, 5), sparse_data=True)
    assert ic.bits_equal(read_zarr(fn), im)
    assert _chunks(fn) == ["0.0.0"]


@pytest.mark.xfail(
    "sparse" not in Z._write_options(),
    strict=True,
    reason="layout difference Z1: cpp-zarr 1.7.3 never writes all-zero chunks; "
    "writezarr (sparseData false) writes them",
)
def test_allZeroChunksAreWritten(tmp_path):
    im = np.zeros((32, 32, 10), np.uint16, order="F")
    im[:16, :16, :5] = 7
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, block_size=(16, 16, 5))
    assert ic.bits_equal(read_zarr(fn), im)
    assert len(_chunks(fn)) == 8


@pytest.mark.parametrize("dtype", ["int16", "uint32", "int8", "int64", "bool"])
def test_createzarrRejectsUnsupportedType(tmp_path, dtype):
    with pytest.raises(ValueError):
        create_zarr(tmp_path / "im.zarr", data_size=(4, 4, 4), dtype=dtype)
    with pytest.raises((ValueError, TypeError)):
        write_zarr(np.zeros((4, 4, 4), dtype), tmp_path / "w.zarr")
    assert not (tmp_path / "w.zarr").exists()


@pytest.mark.skipif("uuid" in Z._write_options(), reason="cpp-zarr takes uuid: tested below")
def test_uuidFalseIsAcceptedWithoutEffect(tmp_path):
    # gap Z4: until cpp-zarr takes uuid, chunks are always written to a temporary name
    fn = tmp_path / "im.zarr"
    data = _rand("uint16", (32, 32, 8), 17)
    write_zarr(data, fn, block_size=(16, 16, 8), uuid=False)
    assert ic.bits_equal(read_zarr(fn), data)


@pytest.mark.skipif("uuid" not in Z._write_options(), reason="cpp-zarr has no uuid (gap Z4)")
def test_uuidFalseWritesChunksInPlace(tmp_path):
    data = ic.make_data("uint16", (32, 32, 8), "smooth", 18)
    for i, (opts, in_place) in enumerate(
        [({}, False), ({"uuid": True}, False), ({"uuid": False}, True)]
    ):
        fn = tmp_path / f"im{i}.zarr"
        write_zarr(ic.make_data("uint16", (32, 32, 8), "smooth", 19), fn, block_size=(16, 16, 8))
        link = tmp_path / f"link{i}"
        os.link(fn / "0.0.0", link)
        write_zarr(data, fn, **opts)
        assert ic.bits_equal(read_zarr(fn), data)
        assert (link.read_bytes() == (fn / "0.0.0").read_bytes()) == in_place


# --- createzarr's .zarray -------------------------------------------------------------------


@pytest.mark.parametrize("sep", [".", "/"])
@pytest.mark.parametrize("order", ["F", "C"])
@pytest.mark.parametrize("compressor", ic.ZARR_COMPRESSORS)
@pytest.mark.parametrize("dtype", ic.ZARR_WRAPPER_DTYPES)
def test_zarrayIsTheMexText(tmp_path, dtype, compressor, order, sep):
    # cpp-zarr writes the mex's text (IO1, checked against the mex for every compressor,
    # order and separator); create_zarr writes the same
    a = _rand(dtype, (9, 7, 5), 20)
    ref = tmp_path / "ref.zarr"
    cppzarr.write_zarr(
        str(ref),
        a,
        order=order,
        chunks=[4, 7, 2],
        dimension_separator=sep,
        cname=compressor,
        clevel=1,
    )
    fn = tmp_path / "im.zarr"
    create_zarr(
        fn,
        data_size=a.shape,
        block_size=(4, 7, 2),
        dtype=ic.MATLAB_CLASS[dtype],
        order=order,
        compressor=compressor,
        dim_separator=sep,
    )
    assert (fn / ".zarray").read_text() == (ref / ".zarray").read_text()
    assert os.listdir(fn) == [".zarray"]


def test_createzarrOptions(tmp_path):
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=(21, 17), block_size=(500, 500, 500), dtype="float32", clevel=5)
    m = _meta(fn)
    assert m["shape"] == [21, 17, 1] and m["chunks"] == [21, 17, 1] and m["dtype"] == "<f4"
    assert m["compressor"]["clevel"] == 5
    create_zarr(fn, data_size=(21, 17), block_size=(8, 8), expand_2d_dim=False)
    assert _meta(fn)["shape"] == [21, 17] and _meta(fn)["chunks"] == [8, 8]
    # the folder is replaced, not added to
    write_zarr(_rand("uint16", (21, 17), 21), fn)
    create_zarr(fn, data_size=(4, 4, 4), compressor="none")
    assert os.listdir(fn) == [".zarray"] and _meta(fn)["compressor"] is None
    with pytest.raises(ValueError):
        create_zarr(fn, compressor="zfp")
    with pytest.raises(ValueError):
        create_zarr(fn, dim_separator="_")


def test_shardSizeNotDividingTheBlockWarns(tmp_path):
    fn = tmp_path / "im.zarr"
    with pytest.warns(UserWarning):
        create_zarr(fn, data_size=IM, block_size=(32, 32, 10), shard_size=(16, 12, 5))
    assert _meta(fn)["codecs"][0]["configuration"]["chunk_shape"] == [16, 32, 5]


def test_zarrSubSize(tmp_path):
    im = _rand("uint16", IM, 22)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, block_size=(8, 8, 4), zarr_sub_size=(2, 2, 2))
    assert _meta(fn)["subfolders"] == [2, 2, 2]
    assert ic.bits_equal(read_zarr(fn), im)
    # more subfolders than chunks: none
    fn2 = tmp_path / "im2.zarr"
    write_zarr(im, fn2, block_size=(16, 16, 5), zarr_sub_size=(10, 10, 10))
    assert "subfolders" not in _meta(fn2)


# --- region writes, conversions -------------------------------------------------------------


def test_regionWriteConvertsAsMatlab(tmp_path):
    # gap Z3: cpp-zarr converts with astype; the mex like uint16() (values measured, IO1)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=(4, 4, 1), block_size=(4, 4, 1), dtype="uint16")
    v = ic.convert_values((4, 4, 1))
    with np.errstate(invalid="raise"):
        write_zarr(v, fn, bbox=[1, 1, 1, 4, 4, 1], create=False)
    assert np.array_equal(read_zarr(fn), ic.matlab_uint16(v)[:, :, 0])


def test_castLikeMatlab():
    c = Z.cast_like_matlab
    s = np.array(
        [0.5, 1.5, 2.5, -0.5, -1.5, 0.49999997, 65535.5, np.nan, np.inf, -np.inf, 1e9], np.float32
    )
    assert c(s, np.uint16).tolist() == [1, 2, 3, 0, 0, 0, 65535, 0, 65535, 0, 65535]
    assert c(s, np.int16).tolist() == [1, 2, 3, -1, -2, 0, 32767, 0, 32767, -32768, 32767]
    assert c(np.array([0.49999999999999994, 2.5, -2.5]), np.int8).tolist() == [0, 3, -3]
    assert c(np.array([300, -5, 70000], np.int32), np.uint8).tolist() == [255, 0, 255]
    assert c(np.array([2**63 + 5], np.uint64), np.int64).tolist() == [2**63 - 1]
    assert c(np.array([-1], np.int64), np.uint64).tolist() == [0]
    assert c(np.array([1e19, -1e19, 2.0**62]), np.int64).tolist() == [2**63 - 1, -(2**63), 2**62]
    assert c(np.array([1e20, 3.5]), np.uint64).tolist() == [2**64 - 1, 4]
    assert c(np.array([1 + 2**-30]), np.float32).dtype == np.float32
    assert c(np.array([True, False]), np.uint16).tolist() == [1, 0]


def test_regionWritesValidateTheBox(tmp_path):
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=IM, block_size=(16, 16, 5))
    with pytest.raises(ValueError):
        write_zarr(np.zeros((5, 5, 5), np.uint16), fn, bbox=[1, 1, 1, 5, 5, 4], create=False)
    with pytest.raises(ValueError):
        write_zarr(np.zeros((5, 5, 5), np.uint16), fn, bbox=[34, 1, 1, 38, 5, 5], create=False)
    with pytest.raises(FileNotFoundError):
        write_zarr(np.zeros((5, 5, 5), np.uint16), tmp_path / "no.zarr", create=False)


def test_2dWithoutExpandedZ(tmp_path):
    # writezarr with expand2dDim false: createzarr still stores [ny, nx, 1]; the 4-value
    # bbox is fitted to the 3 axes
    im = _rand("uint16", (21, 17), 23)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn, expand_2d_dim=False, block_size=(8, 8, 8))
    assert ic.bits_equal_m(read_zarr(fn), im)


def test_writeZarrBlock(tmp_path):
    im = _rand("uint16", IM, 24)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=IM, block_size=(16, 16, 5))
    expected = np.zeros(IM, np.uint16, order="F")
    # an inner block, and an edge block given with data past the edge (trimmed)
    write_zarr_block(fn, [1, 2, 1], im[:16, 16:32, :5])
    expected[:16, 16:32, :5] = im[:16, 16:32, :5]
    edge = _rand("uint16", (16, 16, 5), 25)
    write_zarr_block(fn, [3, 3, 3], edge)
    expected[32:, 32:, 10:] = edge[:5, :9, :3]
    assert ic.bits_equal(read_zarr(fn), expected)


def test_layoutsAndThreadsGiveTheSameChunks(tmp_path):
    im = ic.make_data("uint16", (300, 280, 40), "smooth", 26)
    ref = tmp_path / "ref.zarr"
    cppzarr.write_zarr(str(ref), im, order="F", chunks=[128, 128, 16], dimension_separator=".")
    variants = {
        "F1": (im, 1),
        "F7": (im, 7),
        "C": (np.ascontiguousarray(im), 3),
        "strided": (np.asfortranarray(np.repeat(im, 2, axis=1))[:, ::2], 4),
    }
    for name, (a, t) in variants.items():
        fn = tmp_path / f"{name}.zarr"
        write_zarr(a, fn, block_size=(128, 128, 16), num_threads=t)
        assert _chunks(fn) == _chunks(ref), name
        assert all((fn / c).read_bytes() == (ref / c).read_bytes() for c in _chunks(ref)), name
        assert (fn / ".zarray").read_text() == (ref / ".zarray").read_text()


def test_readsCOrderAndV3ArraysAsFortran(tmp_path):
    im = _rand("float64", IM, 27)
    fn = tmp_path / "c.zarr"
    cppzarr.write_zarr(str(fn), im, order="C", chunks=[16, 16, 5], dimension_separator=".")
    out = read_zarr(fn, input_bbox=[2, 3, 4, 20, 30, 10])
    assert out.flags.f_contiguous and ic.bits_equal(out, im[1:20, 2:30, 3:10])
    v3 = tmp_path / "v3.zarr"
    cppzarr.write_zarr(str(v3), im, chunks=[16, 16, 5], zarr_format=3)
    out = read_zarr(v3)
    assert out.flags.f_contiguous and ic.bits_equal(out, im)


def test_noTemporaryFilesAndGroupWrite(tmp_path):
    fn = tmp_path / "im.zarr"
    old = os.umask(0o022)
    try:
        write_zarr(_rand("uint16", IM, 28), fn, block_size=(16, 16, 5))
        if os.name != "nt":
            assert os.stat(fn).st_mode & 0o020
        fn2 = tmp_path / "im2.zarr"
        create_zarr(fn2, data_size=IM, group_write=False)
        if os.name != "nt":
            assert not os.stat(fn2).st_mode & 0o020
    finally:
        os.umask(old)
    assert [n for n in os.listdir(fn) if n.startswith(".")] == [".zarray"]
