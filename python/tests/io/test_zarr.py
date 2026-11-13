"""petakit5d.io zarr reading, writing and creation (tasks IO2, IO3).

Ports the intent of tests/io/ZarrIOTest.m (test names kept) on the port's native arrays
([nz, ny, nx], order C), plus the reading of MATLAB's arrays ([ny, nx, nz], order F), the
.zarray text, MATLAB's conversions in region writes (gap Z3), the options of writezarr, the
plane transpose of petakit5d._ext.io and the MATLAB oracle cases.
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

from petakit5d.core.arrays import from_matlab, to_matlab  # noqa: E402
from petakit5d.core.bbox import box_from_matlab, box_slices  # noqa: E402
from petakit5d.core.dtypes import cast  # noqa: E402
from petakit5d.io import (  # noqa: E402
    create_zarr,
    get_image_size,
    get_zarr_info,
    read_zarr,
    write_zarr,
    write_zarr_block,
)
from petakit5d.io import zarr as Z  # noqa: E402
from petakit5d.testing import oracle  # noqa: E402
from petakit5d.testing.cases import Case, In, Value, Work  # noqa: E402

IM = (13, 37, 41)  # (z, y, x), never cubic

# .zarray of a native 13 x 37 x 41 uint16 array with blockSize (5, 16, 16)
NATIVE_ZARRAY = (
    '{\n    "chunks": [\n        5,\n        16,\n        16\n    ],\n    "compressor": {\n'
    '        "blocksize": 0,\n        "clevel": 1,\n        "cname": "zstd",\n'
    '        "id": "blosc",\n        "shuffle": 1\n    },\n    "dtype": "<u2",\n'
    '    "fill_value": 0,\n    "filters": null,\n    "order": "C",\n    "shape": [\n'
    '        13,\n        37,\n        41\n    ],\n    "zarr_format": 2\n}\n'
)


def _rand(dtype, shape, seed=1, pattern="random"):
    return np.ascontiguousarray(ic.make_data(dtype, shape, pattern, seed))


def _chunks(fn):
    out = []
    for root, _, files in os.walk(fn):
        out += [
            os.path.relpath(os.path.join(root, f), fn).replace(os.sep, "/")
            for f in files
            if not f.startswith(".") and f != "zarr.json"
        ]
    return sorted(out)


def _meta(fn):
    return json.loads((Path(fn) / ".zarray").read_text())


def _matlab_zarr(fn, m, chunks_yxz, sep="."):
    """an array in MATLAB's layout, as writezarr writes it: m is MATLAB's (y, x, z) array;
    cpp-zarr's order-F write gives the mex's .zarray and chunk bytes (IO1, IO2)"""
    m = m if m.ndim == 3 else m[:, :, np.newaxis]
    cppzarr.write_zarr(
        str(fn),
        np.asfortranarray(m),
        order="F",
        chunks=[min(a, b) for a, b in zip(m.shape, chunks_yxz, strict=True)],
        dimension_separator=sep,
    )


# --- ZarrIOTest ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "block", [(5, 16, 16), IM, (64, 64, 64)], ids=["nonDivisible", "singleChunk", "larger"]
)
@pytest.mark.parametrize("dtype", ic.ZARR_WRAPPER_DTYPES)
def test_roundTripPreservesTypeAndValues(tmp_path, dtype, block):
    vol = _rand(dtype, IM, 2)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=block)
    out = read_zarr(fn)
    assert out.flags.c_contiguous
    assert ic.bits_equal(out, vol)


@pytest.mark.parametrize("compressor", ic.ZARR_COMPRESSORS)
def test_compressorsRoundTrip(tmp_path, compressor):
    vol = _rand("uint16", IM, 3)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=vol.shape, block_size=(5, 16, 16), compressor=compressor)
    write_zarr(vol, fn)
    assert ic.bits_equal(read_zarr(fn), vol)
    assert _meta(fn)["compressor"]["id" if compressor == "gzip" else "cname"] == compressor


@pytest.mark.parametrize("chunks", [(5, 16, 8), (13, 8, 16)], ids=["zSmaller", "zLarger"])
def test_chunkOrderRoundTrip(tmp_path, chunks):
    # (native arrays are C order; MATLAB's F-order files are read in the tests below)
    vol = _rand("uint16", IM, 4)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=vol.shape, block_size=chunks)
    write_zarr(vol, fn)
    out = read_zarr(fn)
    assert out.flags.c_contiguous and ic.bits_equal(out, vol)
    assert _meta(fn)["order"] == "C" and _meta(fn)["chunks"] == list(chunks)


def test_defaultMetadata(tmp_path):
    # a native array: shape and chunks in (z, y, x), order C, chunk files z.y.x
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("uint16", IM, 5), fn, block_size=(5, 16, 16))
    assert (fn / ".zarray").read_text() == NATIVE_ZARRAY
    fn2 = tmp_path / "im2.zarr"
    write_zarr(_rand("uint16", IM, 5), fn2, block_size=(5, 16, 8))  # 3 x 3 x 6 chunks
    assert (fn2 / "2.2.5").is_file() and not (fn2 / "5.2.2").exists()
    assert len(_chunks(fn2)) == 3 * 3 * 6


def test_defaultBlockSize(tmp_path):
    # MATLAB's blockSize [500, 500, 500] ([y x z]) through axes.from_yxz, capped at the size
    assert Z._BLOCK == (500, 500, 500)
    fn = tmp_path / "im.zarr"
    vol = _rand("uint8", (3, 600, 700), 6)
    write_zarr(vol, fn)
    assert _meta(fn)["chunks"] == [3, 500, 500] and _meta(fn)["shape"] == [3, 600, 700]


def test_zarrPythonReadsTheSameArray(tmp_path):
    zarr = pytest.importorskip("zarr")
    vol = _rand("float32", IM, 7)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(5, 16, 8))
    z = zarr.open_array(str(fn), mode="r")
    assert z.shape == IM and z.chunks == (5, 16, 8)
    assert ic.bits_equal(np.asarray(z[...]), vol)


def test_blockSizeIsCappedAtImageSize(tmp_path):
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("uint16", IM, 6), fn, block_size=(64, 64, 64))
    assert _meta(fn)["chunks"] == list(IM)


def test_nestedDimensionSeparator(tmp_path):
    vol = _rand("uint16", IM, 7)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(5, 16, 8), dim_separator="/")
    assert ic.bits_equal(read_zarr(fn), vol)
    assert (fn / "2" / "2" / "5").is_file()
    assert _meta(fn)["dimension_separator"] == "/"


def _zjson(fn):
    return json.loads((Path(fn) / "zarr.json").read_text())


def test_shardedRoundTrip(tmp_path):
    # shardSize: a standard Zarr v3 array with sharding (the user's decision), in (z, y, x):
    # shards of blockSize holding chunks of shardSize
    vol = _rand("uint16", IM, 8)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(10, 32, 32), shard_size=(5, 16, 16))
    assert ic.bits_equal(read_zarr(fn), vol)
    assert not (fn / ".zarray").exists()
    meta = _zjson(fn)
    assert meta["zarr_format"] == 3 and meta["shape"] == list(IM)
    assert meta["chunk_grid"]["configuration"]["chunk_shape"] == [10, 32, 32]
    assert meta["codecs"][0]["name"] == "sharding_indexed"
    assert meta["codecs"][0]["configuration"]["chunk_shape"] == [5, 16, 16]
    assert len(_chunks(fn)) == 8 and (fn / "c.1.1.1").is_file()
    info = get_zarr_info(fn)
    assert info.zarr_format == 3 and info.block_size == (10, 32, 32) and info.dtype == "uint16"
    assert get_image_size(fn) == IM
    # a region write into it
    region = _rand("uint16", (3, 20, 10), 9)
    write_zarr(region, fn, bbox=(4, 10, 30, 7, 30, 40), create=False)
    vol[4:7, 10:30, 30:40] = region
    assert ic.bits_equal(read_zarr(fn), vol)


@pytest.mark.parametrize("sep", [".", "/"])
@pytest.mark.parametrize("compressor", ["zstd", "lz4", "gzip"])
@pytest.mark.parametrize("dtype", ic.ZARR_WRAPPER_DTYPES)
def test_shardedArrayIsCppZarrsArray(tmp_path, dtype, compressor, sep):
    # the zarr.json text and the shard bytes of cpp-zarr's own v3 write
    vol = _rand(dtype, IM, 10, "smooth")
    ref = tmp_path / "ref.zarr"
    cppzarr.write_zarr(
        str(ref),
        vol,
        chunks=[5, 16, 16],
        shards=[10, 32, 32],
        zarr_format=3,
        order="C",
        dimension_separator=sep,
        cname=compressor,
        clevel=1,
    )
    fn = tmp_path / "im.zarr"
    create_zarr(
        fn,
        data_size=IM,
        block_size=(10, 32, 32),
        shard_size=(5, 16, 16),
        dtype=ic.MATLAB_CLASS[dtype],
        compressor=compressor,
        dim_separator=sep,
    )
    assert (fn / "zarr.json").read_text() == (ref / "zarr.json").read_text()
    write_zarr(vol, fn)
    assert _chunks(fn) == _chunks(ref)
    assert all((fn / c).read_bytes() == (ref / c).read_bytes() for c in _chunks(ref))
    zarr = pytest.importorskip("zarr")
    assert ic.bits_equal(np.asarray(zarr.open_array(str(fn), mode="r")[...]), vol)


def test_readsMatlabsV2ShardedArrays(tmp_path):
    # MATLAB's createzarr with shardSize: a v2 .zarray with a sharding_indexed codec, written
    # as by IO1 (cpp-zarr writes regions into it as the mex does); read as (z, y, x)
    m = ic.make_data("uint16", (37, 41, 13), "random", 11)
    fn = tmp_path / "m.zarr"
    fn.mkdir()
    meta = {
        "chunks": [32, 32, 10],
        "codecs": [
            {
                "configuration": {
                    "chunk_shape": [16, 16, 5],
                    "codecs": [
                        {
                            "configuration": {
                                "blocksize": 0,
                                "clevel": 1,
                                "cname": "zstd",
                                "shuffle": "shuffle",
                                "typesize": 2,
                            },
                            "name": "blosc",
                        }
                    ],
                },
                "name": "sharding_indexed",
            }
        ],
        "compressor": {"blocksize": 0, "clevel": 1, "cname": "zstd", "id": "blosc", "shuffle": 1},
        "dtype": "<u2",
        "fill_value": 0,
        "filters": None,
        "order": "F",
        "shape": [37, 41, 13],
        "zarr_format": 2,
    }
    (fn / ".zarray").write_text(json.dumps(meta, indent=4, sort_keys=True) + "\n")
    cppzarr.write_zarr(str(fn), np.asfortranarray(m), [0, 0, 0], [37, 41, 13])
    assert ic.bits_equal(read_zarr(fn), from_matlab(m))
    assert ic.bits_equal(
        read_zarr(fn, input_bbox=(1, 2, 3, 9, 30, 40)), from_matlab(m)[1:9, 2:30, 3:40]
    )


def test_image2dIsStoredWithSingletonZ(tmp_path):
    im = _rand("uint16", (21, 17), 9)
    fn = tmp_path / "im.zarr"
    write_zarr(im, fn)
    out = read_zarr(fn)
    assert out.shape == (21, 17) and out.flags.c_contiguous and ic.bits_equal(out, im)
    assert get_image_size(fn) == (1, 21, 17)
    assert _meta(fn)["shape"] == [1, 21, 17] and _meta(fn)["chunks"] == [1, 21, 17]
    # a 2D box takes the image's slice
    assert ic.bits_equal(read_zarr(fn, input_bbox=(2, 3, 20, 15)), im[2:20, 3:15])


def test_readBoundingBox(tmp_path):
    vol = _rand("uint16", IM, 10)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(5, 16, 16))
    # chunk-crossing, at the borders, a single voxel, a single slice, the whole volume, open
    boxes = [
        (1, 2, 4, 11, 30, 33),
        (0, 0, 0, 13, 37, 41),
        (12, 36, 40, 13, 37, 41),
        (5, 16, 16, 6, 17, 17),
        (6, 0, 0, 7, 37, 41),
        (3, 20, 25, None, None, None),
        (4, 0, 0, 13, 1, 41),
    ]
    for b in boxes:
        want = vol[box_slices(b)]
        out = read_zarr(fn, input_bbox=b)
        assert out.flags.c_contiguous, b
        assert ic.bits_equal(out, want[0] if want.shape[0] == 1 else want), b


@pytest.mark.parametrize(
    "bbox",
    [
        (0, 0, 0, 14, 37, 41),
        (-1, 0, 0, 5, 5, 5),
        (5, 5, 5, 4, 6, 6),
        (5, 5, 5, 5, 6, 6),
        (0, 0, 5, 5),
        (0, 0, 0, 5, 5),
    ],
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
    vol = _rand("uint16", IM, 12)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=IM, block_size=(5, 16, 16))
    expected = np.zeros(IM, np.uint16)
    assert ic.bits_equal(read_zarr(fn), expected)
    b = (1, 4, 6, 9, 30, 20)
    write_zarr(vol[box_slices(b)], fn, bbox=b, create=False)
    expected[box_slices(b)] = vol[box_slices(b)]
    assert ic.bits_equal(read_zarr(fn), expected)
    b = (0, 9, 0, 13, 37, 15)
    region = _rand("uint16", (13, 28, 15), 13)
    write_zarr(region, fn, bbox=b, create=False)
    expected[box_slices(b)] = region
    assert ic.bits_equal(read_zarr(fn), expected)
    assert all(
        len(n.split(".")) == 3 and all(p.isdigit() for p in n.split(".")) for n in _chunks(fn)
    )
    assert (fn / ".zarray").read_text() == NATIVE_ZARRAY


def test_overwriteReplacesExistingArray(tmp_path):
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("float32", (20, 21, 22), 14), fn)
    vol = _rand("uint16", IM, 15)
    write_zarr(vol, fn, overwrite=True)
    assert ic.bits_equal(read_zarr(fn), vol)
    assert _meta(fn)["dtype"] == "<u2"


def test_writeIntoExistingArrayKeepsItsMetadata(tmp_path):
    # without overwrite, the data goes into the array there, converted to its type
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=IM, block_size=(13, 8, 16), dtype="single", compressor="gzip")
    text = (fn / ".zarray").read_text()
    vol = _rand("uint16", IM, 16)
    write_zarr(vol, fn, block_size=(5, 16, 16))
    assert (fn / ".zarray").read_text() == text
    out = read_zarr(fn)
    assert out.dtype == np.float32 and np.array_equal(out, vol.astype(np.float32))


def test_sparseWriteSkipsEmptyChunks(tmp_path):
    vol = np.zeros((10, 32, 32), np.uint16)
    vol[:5, :16, :16] = 7
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(5, 16, 16), sparse_data=True)
    assert ic.bits_equal(read_zarr(fn), vol)
    assert _chunks(fn) == ["0.0.0"]


@pytest.mark.xfail(
    "sparse" not in Z._write_options(),
    strict=True,
    reason="layout difference Z1: cpp-zarr 1.7.3 never writes all-zero chunks; "
    "writezarr (sparseData false) writes them",
)
def test_allZeroChunksAreWritten(tmp_path):
    vol = np.zeros((10, 32, 32), np.uint16)
    vol[:5, :16, :16] = 7
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(5, 16, 16))
    assert ic.bits_equal(read_zarr(fn), vol)
    assert len(_chunks(fn)) == 8


@pytest.mark.parametrize("dtype", ["int16", "uint32", "int8", "int64", "bool"])
def test_createzarrRejectsUnsupportedType(tmp_path, dtype):
    with pytest.raises(ValueError):
        create_zarr(tmp_path / "im.zarr", data_size=(4, 5, 6), dtype=dtype)
    with pytest.raises((ValueError, TypeError)):
        write_zarr(np.zeros((4, 5, 6), dtype), tmp_path / "w.zarr")
    assert not (tmp_path / "w.zarr").exists()


@pytest.mark.skipif("uuid" in Z._write_options(), reason="cpp-zarr takes uuid: tested below")
def test_uuidFalseIsAcceptedWithoutEffect(tmp_path):
    # gap Z4: until cpp-zarr takes uuid, chunks are always written to a temporary name
    fn = tmp_path / "im.zarr"
    data = _rand("uint16", (8, 32, 33), 17)
    write_zarr(data, fn, block_size=(8, 16, 16), uuid=False)
    assert ic.bits_equal(read_zarr(fn), data)


@pytest.mark.skipif("uuid" not in Z._write_options(), reason="cpp-zarr has no uuid (gap Z4)")
def test_uuidFalseWritesChunksInPlace(tmp_path):
    data = _rand("uint16", (8, 32, 33), 18, "smooth")
    for i, (opts, in_place) in enumerate(
        [({}, False), ({"uuid": True}, False), ({"uuid": False}, True)]
    ):
        fn = tmp_path / f"im{i}.zarr"
        write_zarr(_rand("uint16", (8, 32, 33), 19, "smooth"), fn, block_size=(8, 16, 16))
        link = tmp_path / f"link{i}"
        os.link(fn / "0.0.0", link)
        write_zarr(data, fn, **opts)
        assert ic.bits_equal(read_zarr(fn), data)
        assert (link.read_bytes() == (fn / "0.0.0").read_bytes()) == in_place


# --- createzarr's .zarray -------------------------------------------------------------------


@pytest.mark.parametrize("sep", [".", "/"])
@pytest.mark.parametrize("compressor", ic.ZARR_COMPRESSORS)
@pytest.mark.parametrize("dtype", ic.ZARR_WRAPPER_DTYPES)
def test_zarrayIsTheMexText(tmp_path, dtype, compressor, sep):
    # cpp-zarr writes the mex's text (IO1); create_zarr writes the same for native arrays
    a = _rand(dtype, (5, 9, 7), 20)
    ref = tmp_path / "ref.zarr"
    cppzarr.write_zarr(
        str(ref),
        a,
        order="C",
        chunks=[2, 4, 7],
        dimension_separator=sep,
        cname=compressor,
        clevel=1,
    )
    fn = tmp_path / "im.zarr"
    create_zarr(
        fn,
        data_size=a.shape,
        block_size=(2, 4, 7),
        dtype=ic.MATLAB_CLASS[dtype],
        compressor=compressor,
        dim_separator=sep,
    )
    assert (fn / ".zarray").read_text() == (ref / ".zarray").read_text()
    assert os.listdir(fn) == [".zarray"]


def test_createzarrOptions(tmp_path):
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=(21, 17), block_size=(500, 500, 500), dtype="float32", clevel=5)
    m = _meta(fn)
    assert m["shape"] == [1, 21, 17] and m["chunks"] == [1, 21, 17] and m["dtype"] == "<f4"
    assert m["compressor"]["clevel"] == 5 and m["order"] == "C"
    create_zarr(fn, data_size=(21, 17), block_size=(8, 9), expand_2d_dim=False)
    assert _meta(fn)["shape"] == [21, 17] and _meta(fn)["chunks"] == [8, 9]
    write_zarr(_rand("uint16", (21, 17), 21), fn)  # a 2D array written as it is
    create_zarr(fn, data_size=(4, 5, 6), compressor="none")
    assert os.listdir(fn) == [".zarray"] and _meta(fn)["compressor"] is None
    with pytest.raises(ValueError):
        create_zarr(fn, compressor="zfp")
    with pytest.raises(ValueError):
        create_zarr(fn, dim_separator="_")


def test_shardSizeNotDividingTheBlockWarns(tmp_path):
    fn = tmp_path / "im.zarr"
    with pytest.warns(UserWarning):
        create_zarr(fn, data_size=IM, block_size=(10, 32, 32), shard_size=(5, 12, 16))
    assert _zjson(fn)["codecs"][0]["configuration"]["chunk_shape"] == [5, 32, 16]
    with pytest.raises(ValueError):
        create_zarr(
            fn, data_size=IM, block_size=(4, 8, 8), shard_size=(2, 4, 4), zarr_sub_size=(2, 2, 2)
        )


def test_zarrSubSize(tmp_path):
    vol = _rand("uint16", IM, 22)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(4, 8, 8), zarr_sub_size=(2, 2, 3))
    assert _meta(fn)["subfolders"] == [2, 2, 3]
    assert ic.bits_equal(read_zarr(fn), vol)
    fn2 = tmp_path / "im2.zarr"
    write_zarr(vol, fn2, block_size=(5, 16, 16), zarr_sub_size=(10, 10, 10))
    assert "subfolders" not in _meta(fn2)


# --- region writes, conversions -------------------------------------------------------------


def test_regionWriteConvertsAsMatlab(tmp_path):
    # gap Z3: cpp-zarr converts with astype; the mex like MATLAB's cast (values measured, IO1)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=(1, 4, 4), block_size=(1, 4, 4), dtype="uint16")
    v = np.ascontiguousarray(ic.convert_values((4, 4)))
    with np.errstate(invalid="raise"):
        write_zarr(v, fn, bbox=(0, 0, 0, 1, 4, 4), create=False)
    assert np.array_equal(read_zarr(fn), ic.matlab_uint16(v))
    assert np.array_equal(read_zarr(fn), cast(v, "uint16"))


def test_regionWritesValidateTheBox(tmp_path):
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=IM, block_size=(5, 16, 16))
    with pytest.raises(ValueError):
        write_zarr(np.zeros((5, 5, 5), np.uint16), fn, bbox=(0, 0, 0, 4, 5, 5), create=False)
    with pytest.raises(ValueError):
        write_zarr(np.zeros((5, 5, 5), np.uint16), fn, bbox=(0, 33, 0, 5, 38, 5), create=False)
    with pytest.raises(FileNotFoundError):
        write_zarr(np.zeros((5, 5, 5), np.uint16), tmp_path / "no.zarr", create=False)
    with pytest.raises(ValueError):
        write_zarr(np.zeros((3, 4, 5), np.uint16, order="F"), tmp_path / "f.zarr")


def test_writeZarrBlock(tmp_path):
    vol = _rand("uint16", IM, 24)
    fn = tmp_path / "im.zarr"
    create_zarr(fn, data_size=IM, block_size=(5, 16, 16))
    expected = np.zeros(IM, np.uint16)
    # an inner block, and an edge block given with data past the edge (trimmed)
    write_zarr_block(fn, (0, 1, 0), vol[:5, 16:32, :16])
    expected[:5, 16:32, :16] = vol[:5, 16:32, :16]
    edge = _rand("uint16", (5, 16, 16), 25)
    write_zarr_block(fn, (2, 2, 2), edge)
    expected[10:, 32:, 32:] = edge[:3, :5, :9]
    assert ic.bits_equal(read_zarr(fn), expected)


def test_layoutsAndThreadsGiveTheSameChunks(tmp_path):
    vol = _rand("uint16", (40, 300, 280), 26, "smooth")
    ref = tmp_path / "ref.zarr"
    cppzarr.write_zarr(str(ref), vol, order="C", chunks=[16, 128, 128], dimension_separator=".")
    variants = {
        "C1": (vol, 1),
        "C7": (vol, 7),
        "strided": (np.repeat(vol, 2, axis=2)[:, :, ::2], 4),
    }
    for name, (a, t) in variants.items():
        fn = tmp_path / f"{name}.zarr"
        write_zarr(a, fn, block_size=(16, 128, 128), num_threads=t)
        assert _chunks(fn) == _chunks(ref), name
        assert all((fn / c).read_bytes() == (ref / c).read_bytes() for c in _chunks(ref)), name
        assert (fn / ".zarray").read_text() == (ref / ".zarray").read_text()


def test_readsV3Arrays(tmp_path):
    vol = _rand("float64", IM, 27)
    v3 = tmp_path / "v3.zarr"
    cppzarr.write_zarr(str(v3), vol, chunks=[5, 16, 16], zarr_format=3)
    out = read_zarr(v3, input_bbox=(1, 2, 3, 10, 20, 30))
    assert out.flags.c_contiguous and ic.bits_equal(out, vol[1:10, 2:20, 3:30])


def test_noTemporaryFilesAndGroupWrite(tmp_path):
    fn = tmp_path / "im.zarr"
    old = os.umask(0o022)
    try:
        write_zarr(_rand("uint16", IM, 28), fn, block_size=(5, 16, 16))
        if os.name != "nt":
            assert os.stat(fn).st_mode & 0o020
        fn2 = tmp_path / "im2.zarr"
        create_zarr(fn2, data_size=IM, group_write=False)
        if os.name != "nt":
            assert not os.stat(fn2).st_mode & 0o020
    finally:
        os.umask(old)
    assert [n for n in os.listdir(fn) if n.startswith(".")] == [".zarray"]


# --- arrays in MATLAB's layout ([ny, nx, nz], order F) ----------------------------------------


@pytest.mark.parametrize("dtype", ic.ZARR_WRAPPER_DTYPES + ic.ZARR_OTHER_DTYPES)
def test_readsMatlabLayoutAsFromMatlab(tmp_path, dtype):
    # every class, chunks not dividing the size; whole and by boxes at the borders
    m = ic.make_data(dtype, (37, 41, 13), "random", ic.seed_of(dtype))  # MATLAB's (y, x, z)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (16, 16, 5))
    want = from_matlab(m)
    out = read_zarr(fn)
    assert out.flags.c_contiguous and ic.bits_equal(out, want)
    for mb in (
        [3, 5, 2, 30, 33, 11],
        [1, 1, 1, 37, 41, 13],
        [37, 41, 13, 37, 41, 13],
        [17, 1, 6, 17, 41, 6],
        [1, 30, 1, 37, 41, 13],
    ):
        b = box_from_matlab(mb)
        sub = want[box_slices(b)]
        got = read_zarr(fn, input_bbox=b)
        assert got.flags.c_contiguous, mb
        assert ic.bits_equal(got, sub[0] if sub.shape[0] == 1 else sub), mb


def test_readsMatlabLayoutImage(tmp_path):
    m = ic.make_data("uint16", (21, 17), "random", 30)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (500, 500, 500))  # [21, 17, 1], as writezarr stores 2D data
    out = read_zarr(fn)
    assert out.shape == (21, 17) and ic.bits_equal(out, m)
    assert get_image_size(fn) == (1, 21, 17)
    assert ic.bits_equal(read_zarr(fn, input_bbox=(2, 3, 20, 15)), m[2:20, 3:15])
    # MATLAB's 2D array stored [ny, nx] (no z)
    fn2 = tmp_path / "m2.zarr"
    cppzarr.write_zarr(
        str(fn2), np.asfortranarray(m), order="F", chunks=[8, 8], dimension_separator="."
    )
    assert ic.bits_equal(read_zarr(fn2), m) and get_image_size(fn2) == (21, 17)


def test_matlabLayoutInfo(tmp_path):
    m = ic.make_data("float32", (37, 41, 13), "random", 31)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (16, 8, 5))
    info = get_zarr_info(fn)
    assert info.axes == "yxz" and info.order == "F"
    assert info.size == (13, 37, 41) and info.block_size == (5, 16, 8)
    assert get_image_size(fn) == (13, 37, 41)
    native = tmp_path / "n.zarr"
    write_zarr(from_matlab(m), native, block_size=(5, 16, 8))
    ninfo = get_zarr_info(native)
    assert ninfo.axes == "zyx" and ninfo.size == info.size and ninfo.block_size == (5, 16, 8)


def test_axesOverridesTheDetection(tmp_path):
    vol = _rand("uint16", (6, 7, 9), 32)
    fn = tmp_path / "n.zarr"
    write_zarr(vol, fn)
    # the native array read as if it were MATLAB's [ny, nx, nz]: from_matlab of that array
    assert ic.bits_equal(read_zarr(fn, axes="yxz"), from_matlab(vol))
    assert ic.bits_equal(read_zarr(fn, axes="zyx"), vol)
    m = ic.make_data("uint16", (7, 9, 6), "random", 33)
    fm = tmp_path / "m.zarr"
    _matlab_zarr(fm, m, (4, 4, 4))
    # an F-order array taken as (z, y, x) as it is
    out = read_zarr(fm, axes="zyx")
    assert out.shape == (7, 9, 6) and out.flags.c_contiguous and ic.bits_equal(out, m)
    with pytest.raises(ValueError):
        read_zarr(fm, axes="xyz")


def test_regionWriteIntoMatlabLayoutRaises(tmp_path):
    m = ic.make_data("uint16", (37, 41, 13), "random", 34)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (16, 16, 5))
    before = {c: (fn / c).read_bytes() for c in _chunks(fn)}
    with pytest.raises(ValueError, match="MATLAB's layout"):
        write_zarr(np.zeros((2, 3, 4), np.uint16), fn, bbox=(0, 0, 0, 2, 3, 4), create=False)
    with pytest.raises(ValueError, match="MATLAB's layout"):
        write_zarr(from_matlab(m), fn)
    with pytest.raises(ValueError, match="MATLAB's layout"):
        write_zarr_block(fn, (0, 0, 0), np.zeros((5, 16, 16), np.uint16))
    assert {c: (fn / c).read_bytes() for c in _chunks(fn)} == before
    # overwrite replaces it with a native array
    write_zarr(from_matlab(m), fn, overwrite=True)
    assert get_zarr_info(fn).axes == "zyx" and ic.bits_equal(read_zarr(fn), from_matlab(m))


def test_matlabLayoutReadWithoutTheExtension(tmp_path, monkeypatch):
    m = ic.make_data("uint16", (37, 41, 13), "random", 35)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (16, 16, 5))
    monkeypatch.setattr(Z, "_inplace_kernel", lambda: None)
    assert ic.bits_equal(read_zarr(fn), from_matlab(m))
    assert ic.bits_equal(
        read_zarr(fn, input_bbox=(1, 2, 3, 9, 30, 40)), from_matlab(m)[1:9, 2:30, 3:40]
    )


def test_matlabLayoutReadIsOneVolume(tmp_path):
    # the planes are turned around in cpp-zarr's buffer: no second volume
    import tracemalloc

    from petakit5d import _ext

    if not _ext.available("io"):
        pytest.skip("petakit5d._ext.io is not built")
    m = ic.make_data("uint16", (300, 200, 20), "random", 36)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (128, 128, 8))
    tracemalloc.start()
    out = read_zarr(fn)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert ic.bits_equal(out, from_matlab(m))
    assert peak <= 1.05 * m.nbytes


def test_readThreadsGiveTheSameResult(tmp_path):
    m = ic.make_data("uint16", (300, 257, 9), "random", 37)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (128, 128, 4))
    for t in (1, 2, 5):
        assert ic.bits_equal(read_zarr(fn, num_threads=t), from_matlab(m))


# --- reading into a given array or a reused buffer (out=) ------------------------------------


@pytest.mark.parametrize("alloc_empty", [False, True], ids=["ownOutput", "npEmptyOutput"])
def test_readIntoOut(tmp_path, monkeypatch, alloc_empty):

    monkeypatch.setattr(Z, "_ALLOC_EMPTY", alloc_empty)
    vol = _rand("uint16", IM, 40)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(5, 16, 16))
    assert ic.bits_equal(read_zarr(fn), vol) and read_zarr(fn).flags.c_contiguous
    o = np.empty(IM, np.uint16)
    assert read_zarr(fn, out=o) is o and ic.bits_equal(o, vol)
    b = (1, 2, 3, 9, 30, 40)
    o = np.zeros((8, 28, 37), np.uint16)
    assert read_zarr(fn, input_bbox=b, out=o) is o and ic.bits_equal(o, vol[box_slices(b)])
    # one slice: out (ny, nx) or (1, ny, nx)
    for shape in [(37, 41), (1, 37, 41)]:
        o = np.empty(shape, np.uint16)
        assert read_zarr(fn, input_bbox=(6, 0, 0, 7, 37, 41), out=o) is o
        assert ic.bits_equal(o.reshape(37, 41), vol[6])


def test_readBufferIsReusedOverBlocks(tmp_path):
    # edge blocks read into the first bytes of the buffer give a fresh read's bits
    from petakit5d.io import ReadBuffer

    vol = _rand("float32", IM, 41)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(5, 16, 16))
    buf = ReadBuffer()
    first = None
    for z0 in range(0, 13, 5):
        for y0 in range(0, 37, 16):
            for x0 in range(0, 41, 16):
                b = (z0, y0, x0, min(z0 + 5, 13), min(y0 + 16, 37), min(x0 + 16, 41))
                got = read_zarr(fn, input_bbox=b, out=buf)
                first = got.base if first is None else first
                assert got.flags.c_contiguous and got.base is first  # one buffer
                want = read_zarr(fn, input_bbox=b)
                assert ic.bits_equal(got, want), b
    # a larger read or another class makes a new buffer
    assert ic.bits_equal(read_zarr(fn, out=buf), vol)
    write_zarr(vol.astype(np.float64), tmp_path / "d.zarr")
    assert ic.bits_equal(read_zarr(tmp_path / "d.zarr", out=buf), vol.astype(np.float64))


def test_readIntoOutChecksIt(tmp_path):
    fn = tmp_path / "im.zarr"
    write_zarr(_rand("uint16", IM, 42), fn)
    bad = {
        "dtype": np.empty(IM, np.uint8),
        "shape": np.empty((13, 37, 40), np.uint16),
        "permuted": np.empty((41, 37, 13), np.uint16),  # as many elements, other axes
        "fortran": np.empty(IM, np.uint16, order="F"),
        "strided": np.empty((13, 37, 82), np.uint16)[:, :, ::2],
        "readonly": np.empty(IM, np.uint16),
        "list": [0] * 10,
    }
    bad["readonly"].flags.writeable = False
    for name, o in bad.items():
        with pytest.raises((TypeError, ValueError)):
            read_zarr(fn, out=o)
        assert name


@pytest.mark.parametrize("kernel", [True, False], ids=["inPlace", "noExtension"])
def test_readMatlabLayoutIntoOut(tmp_path, monkeypatch, kernel):
    from petakit5d.io import ReadBuffer

    if not kernel:
        monkeypatch.setattr(Z, "_inplace_kernel", lambda: None)
    m = ic.make_data("uint16", (37, 41, 13), "random", 43)
    fn = tmp_path / "m.zarr"
    _matlab_zarr(fn, m, (16, 16, 5))
    want = from_matlab(m)
    o = np.empty(IM, np.uint16)
    assert read_zarr(fn, out=o) is o and ic.bits_equal(o, want)
    buf = ReadBuffer()
    for b in [
        (0, 0, 0, 13, 37, 41),
        (1, 2, 3, 9, 30, 40),
        (12, 36, 40, 13, 37, 41),
        (5, 0, 0, 6, 37, 41),
    ]:
        sub = want[box_slices(b)]
        got = read_zarr(fn, input_bbox=b, out=buf)
        assert ic.bits_equal(got, sub[0] if sub.shape[0] == 1 else sub), b
    # MATLAB's 2D array [ny, nx]
    m2 = ic.make_data("uint16", (21, 17), "random", 44)
    fn2 = tmp_path / "m2.zarr"
    cppzarr.write_zarr(
        str(fn2), np.asfortranarray(m2), order="F", chunks=[8, 8], dimension_separator="."
    )
    o2 = np.empty((21, 17), np.uint16)
    assert read_zarr(fn2, out=o2) is o2 and ic.bits_equal(o2, m2)


def test_readIntoOutMakesNoAllocation(tmp_path):
    import tracemalloc

    from petakit5d.io import ReadBuffer

    vol = _rand("uint16", (20, 300, 200), 45)
    fn = tmp_path / "im.zarr"
    write_zarr(vol, fn, block_size=(8, 128, 128))
    buf = ReadBuffer()
    read_zarr(fn, out=buf)
    tracemalloc.start()
    for b in [(0, 0, 0, 10, 300, 200), (10, 0, 0, 20, 300, 200), (0, 0, 0, 20, 300, 200)]:
        read_zarr(fn, input_bbox=b, out=buf)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    # (a new output of the smallest box would be half the volume; small temporaries of the
    # packages stay far below)
    assert peak < 0.2 * vol.nbytes


# --- the plane transpose (petakit5d._ext.io) --------------------------------------------------


def _ext_io():
    from petakit5d import _ext

    if not _ext.available("io"):
        pytest.skip("petakit5d._ext.io is not built")
    return _ext.load("io")


@pytest.mark.parametrize("mode", [0, 1, 2])
@pytest.mark.parametrize("num_threads", [1, 3, 8])
@pytest.mark.parametrize(
    "dtype", ["uint8", "int16", "uint16", "float32", "float64", "uint64", "bool"]
)
def test_transposeKernelIsExact(dtype, num_threads, mode):
    # pages of rows (x fastest) into pages of columns, in place; sizes not divisible by the
    # tiles (64) or the SIMD blocks (2, 4, 8); fewer and more pages than threads; the modes:
    # chosen, whole pages per thread, pages shared by the threads
    k = _ext_io().transpose_pages_inplace
    for nx, ny, nz in [
        (1, 1, 1),
        (7, 5, 3),
        (64, 64, 2),
        (130, 67, 3),
        (9, 200, 1),
        (1031, 3, 2),
        (33, 17, 20),
        (300, 129, 9),
    ]:
        vol = _rand(dtype, (nz, ny, nx), nx * ny)
        buf = vol.reshape(-1).copy()
        k(buf, nx, ny, nz, num_threads, mode)
        # each plane transposed: (z, x, y) in C order
        assert ic.bits_equal(buf.reshape(nz, nx, ny), np.ascontiguousarray(vol.transpose(0, 2, 1)))


def test_inplaceKernelOfLargePages():
    # pages too large for a scratch page per thread: the threads share each page
    k = _ext_io().transpose_pages_inplace
    vol = _rand("uint16", (5, 3001, 2900), 9)  # 17 MB pages, 87 MB in all
    buf = vol.reshape(-1).copy()
    k(buf, 2900, 3001, 5, 8)
    assert ic.bits_equal(buf.reshape(5, 2900, 3001), np.ascontiguousarray(vol.transpose(0, 2, 1)))


def test_transposeKernelChecksItsArguments():
    k = _ext_io().transpose_pages_inplace
    with pytest.raises((TypeError, ValueError)):
        k(np.zeros(40, np.uint16), 4, 5, 3, 1)  # not nx * ny * nz elements
    with pytest.raises((TypeError, ValueError)):
        k(np.zeros(40, np.uint16)[::2], 4, 5, 1, 1)  # not contiguous, not copied
    with pytest.raises((TypeError, ValueError)):
        k(np.zeros(40, np.uint16), 4, 5, 2, 0)
    with pytest.raises((TypeError, ValueError)):
        k(np.zeros(40, np.uint16), 4, 5, 2, 1, 3)  # no such mode


# --- MATLAB oracle -------------------------------------------------------------------------------
# MATLAB writes with writezarr and reads with readzarr (boxes at the borders, chunk remainders);
# Python reads MATLAB's layout (cpp-zarr's order-F write is the mex's file, IO1) and writes and
# reads its native array: both must give from_matlab of MATLAB's result. Boxes, sizes and
# outputs are converted by hand with the core helpers until F14's roles do it.


def _zcase_py(m, block_yxz, bbox, work):
    from petakit5d.core.axes import from_yxz

    fm, fn = work + "m.zarr", work + "n.zarr"
    _matlab_zarr(fm, m, block_yxz)
    b = box_from_matlab(bbox) if len(bbox) else None
    got_m = read_zarr(fm, input_bbox=b)
    write_zarr(from_matlab(m), fn, block_size=from_yxz(block_yxz))
    got_n = read_zarr(fn, input_bbox=b)
    if not ic.bits_equal(got_m, got_n):
        raise AssertionError("the MATLAB-layout and the native read differ")
    return to_matlab(got_m)


def _zregion_py(size_yxz, cls, region, bbox, work):
    # createzarr, then a region write of another class (converted as MATLAB converts)
    from petakit5d.core.axes import from_yxz

    fn = work + "n.zarr"
    create_zarr(fn, data_size=from_yxz(size_yxz), block_size=from_yxz((16, 16, 5)), dtype=cls)
    write_zarr(from_matlab(region), fn, bbox=box_from_matlab(bbox), create=False)
    return to_matlab(read_zarr(fn))


def _special_single(shape):
    s = np.array(
        [
            0.5,
            1.5,
            2.5,
            -0.5,
            -1.5,
            0.49999997,
            65535.5,
            np.nan,
            np.inf,
            -np.inf,
            1e9,
            255.5,
            254.49998,
            -0.0,
            3.4999998,
            100.5,
        ],
        np.float32,
    )
    n = int(np.prod(shape))
    v = np.resize(s, n) + np.repeat(np.arange((n + 15) // 16), 16)[:n].astype(np.float32) * 3.25
    v[:16] = s
    return v.reshape(shape, order="F")


_ZBOXES = [
    [],
    [3, 5, 2, 30, 33, 11],
    [1, 1, 13, 37, 41, 13],
    [37, 41, 1, 37, 41, 13],
    [17, 17, 6, 17, 17, 6],
    [33, 1, 11, 37, 16, 13],
]
_ZCASES = []
for _dt in ic.ZARR_WRAPPER_DTYPES:
    for _k, _bb in enumerate(_ZBOXES):
        _ZCASES.append(
            Case(
                f"zarr_{_dt}_box{_k}",
                code=(
                    "fn = fullfile(work, 'm.zarr'); writezarr(in{1}, fn, 'blockSize', in{2}); "
                    "if isempty(in{3}), out{1} = readzarr(fn); "
                    "else, out{1} = readzarr(fn, 'inputBbox', in{3}); end"
                ),
                args=[In("m"), [16, 16, 5], _bb],
                python=_zcase_py,
                py_args=[In("m"), (16, 16, 5), _bb, Work()],
                inputs={
                    "m": Value(
                        lambda dt=_dt: ic.make_data(dt, (37, 41, 13), "random", ic.seed_of(dt))
                    )
                },
                files=False,
            )
        )
_ZCASES.append(
    Case(
        "zarr_image2d",
        code="fn = fullfile(work, 'm.zarr'); writezarr(in{1}, fn); out{1} = readzarr(fn);",
        args=[In("m")],
        python=lambda m, work: _zcase_py(m, (500, 500, 500), [], work),
        py_args=[In("m"), Work()],
        inputs={"m": Value(lambda: ic.make_data("uint16", (21, 17), "random", 38))},
        files=False,
    )
)
for _name, _cls, _make in [
    ("dbl_u16", "uint16", lambda: ic.convert_values((26, 14, 8))),
    ("sgl_u16", "uint16", lambda: _special_single((26, 14, 8))),
    ("sgl_u8", "uint8", lambda: _special_single((26, 14, 8))),
    ("dbl_sgl", "single", lambda: ic.convert_values((26, 14, 8))),
]:
    _ZCASES.append(
        Case(
            f"zarr_convert_{_name}",
            code=(
                "fn = fullfile(work, 'n.zarr'); createzarr(fn, 'dataSize', in{1}, 'blockSize', "
                "[16, 16, 5], 'dtype', in{2}); writezarr(in{3}, fn, 'bbox', in{4}, "
                "'create', false); out{1} = readzarr(fn);"
            ),
            args=[[37, 41, 13], _cls, In("r"), [5, 7, 2, 30, 20, 9]],
            python=_zregion_py,
            py_args=[(37, 41, 13), _cls, In("r"), [5, 7, 2, 30, 20, 9], Work()],
            inputs={"r": Value(_make)},
            files=False,
        )
    )


@oracle.parametrize(_ZCASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)
