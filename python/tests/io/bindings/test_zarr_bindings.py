"""cpp-zarr's Python package as PetaKit5D's zarr reader and writer.

Ports the intent of tests/io/ZarrIOTest.m (test names kept) for the binding, plus
the layout and gap checks of task IO1. PetaKit5D's arrays: Zarr v2, order F,
shape [ny, nx, nz], blosc zstd level 1, '.' separator, all chunks written.
"""

import os

import numpy as np
import pytest

import iobind_cases as ic

cppzarr = pytest.importorskip('cppzarr')

IM = (37, 41, 13)

# .zarray as createzarr/writezarr (cpp-zarr mex, PetaKit5D at 5df18b2) write it for a
# 37 x 41 x 13 uint16 array with blockSize [16, 16, 5]
MATLAB_ZARRAY = (
    '{\n    "chunks": [\n        16,\n        16,\n        5\n    ],\n    "compressor": {\n'
    '        "blocksize": 0,\n        "clevel": 1,\n        "cname": "zstd",\n        "id": "blosc",\n'
    '        "shuffle": 1\n    },\n    "dtype": "<u2",\n    "fill_value": 0,\n    "filters": null,\n'
    '    "order": "F",\n    "shape": [\n        37,\n        41,\n        13\n    ],\n    "zarr_format": 2\n}\n')


def _write(fn, a, chunks=(16, 16, 5), **kw):
    """write_zarr with PetaKit5D's settings; chunks capped at the size, as writezarr does."""
    a = a.reshape(a.shape + (1,) * (3 - a.ndim), order='F')
    kw = dict(dict(order='F', dimension_separator='.', cname='zstd', clevel=1), **kw)
    cppzarr.write_zarr(str(fn), a, chunks=[min(s, c) for s, c in zip(a.shape, chunks)], **kw)


def _chunks(fn):
    out = []
    for root, _, files in os.walk(fn):
        out += [os.path.relpath(os.path.join(root, f), fn).replace(os.sep, '/') for f in files if not f.startswith('.')]
    return sorted(out)


@pytest.mark.parametrize('chunks', [(16, 16, 5), (37, 41, 13), (64, 64, 64)], ids=['nonDivisible', 'singleChunk', 'larger'])
@pytest.mark.parametrize('dtype', ic.ZARR_WRAPPER_DTYPES + ic.ZARR_OTHER_DTYPES)
def test_roundTripPreservesTypeAndValues(tmp_path, dtype, chunks):
    im = ic.make_data(dtype, IM, 'random', 2)
    fn = tmp_path / 'im.zarr'
    _write(fn, im, chunks)
    out = cppzarr.read_zarr(str(fn))
    assert out.flags.f_contiguous
    assert ic.bits_equal(out, im)


@pytest.mark.parametrize('compressor', ic.ZARR_COMPRESSORS)
def test_compressorsRoundTrip(tmp_path, compressor):
    import json
    im = ic.make_data('uint16', IM, 'random', 3)
    fn = tmp_path / 'im.zarr'
    _write(fn, im, cname=compressor)
    assert ic.bits_equal(cppzarr.read_zarr(str(fn)), im)
    meta = json.load(open(fn / '.zarray'))
    assert meta['compressor']['id' if compressor == 'gzip' else 'cname'] == compressor


@pytest.mark.parametrize('chunks', [(16, 8, 5), (8, 16, 13)], ids=['zSmaller', 'zLarger'])
@pytest.mark.parametrize('order', ['F', 'C'])
def test_chunkOrderRoundTrip(tmp_path, order, chunks):
    import json
    im = ic.make_data('uint16', IM, 'random', 4)
    fn = tmp_path / 'im.zarr'
    _write(fn, im, chunks, order=order)
    assert ic.bits_equal(cppzarr.read_zarr(str(fn), order='F'), im)
    assert json.load(open(fn / '.zarray'))['order'] == order
    # without order, the array comes back in the file's order
    assert cppzarr.read_zarr(str(fn)).flags['C_CONTIGUOUS' if order == 'C' else 'F_CONTIGUOUS']


def test_defaultMetadata(tmp_path):
    # the .zarray text is MATLAB's, byte for byte; chunk files are y.x.z, 0-based
    fn = tmp_path / 'im.zarr'
    _write(fn, ic.make_data('uint16', IM, 'random', 5))
    assert (fn / '.zarray').read_text() == MATLAB_ZARRAY
    assert (fn / '2.2.2').is_file() and not (fn / '3.0.0').exists()


def test_packageDefaultsAreNotPetaKits(tmp_path):
    # write_zarr's own defaults: 256 chunks and '/' (PetaKit5D passes its settings)
    import json
    fn = tmp_path / 'im.zarr'
    cppzarr.write_zarr(str(fn), ic.make_data('uint16', (300, 20, 3), 'random', 6))
    meta = json.load(open(fn / '.zarray'))
    assert meta['chunks'] == [256, 256, 256] and meta.get('dimension_separator') == '/'


def test_chunksLargerThanTheArrayAreKept(tmp_path):
    # MATLAB's writezarr caps blockSize at the size (blockSizeIsCappedAtImageSize);
    # the package writes the chunks it is given, so the caller caps them
    import json
    fn = tmp_path / 'im.zarr'
    cppzarr.write_zarr(str(fn), ic.make_data('uint16', IM, 'random', 7), chunks=[64, 64, 64], order='F', dimension_separator='.')
    assert json.load(open(fn / '.zarray'))['chunks'] == [64, 64, 64]


def test_blockSizeIsCappedAtImageSize(tmp_path):
    import json
    fn = tmp_path / 'im.zarr'
    _write(fn, ic.make_data('uint16', IM, 'random', 7), (64, 64, 64))
    assert json.load(open(fn / '.zarray'))['chunks'] == list(IM)


def test_matlabUint16Helper():
    # the reference for region-write conversion: MATLAB R2024b's uint16() of these doubles
    # (and the mex's conversion in writezarr, which gave the same values)
    v = ic.convert_values((16,))
    assert ic.matlab_uint16(v).tolist() == [3, 4, 0, 0, 65535, 0, 65535, 0, 65535, 65535, 0, 2, 101, 65535, 0, 1]


def test_matlabShapeHelper():
    assert ic.matlab_shape((5,)) == (5, 1)
    assert ic.matlab_shape((3, 1, 1)) == (3, 1)
    assert ic.matlab_shape((1, 1, 1)) == (1, 1)
    assert ic.matlab_shape((2, 3, 1, 4)) == (2, 3, 1, 4)


def test_nestedDimensionSeparator(tmp_path):
    im = ic.make_data('uint16', IM, 'random', 8)
    fn = tmp_path / 'im.zarr'
    _write(fn, im, dimension_separator='/')
    assert ic.bits_equal(cppzarr.read_zarr(str(fn)), im)
    assert (fn / '2' / '2' / '2').is_file()


def test_shardedRoundTrip(tmp_path):
    # the package shards only Zarr v3 arrays (see test_gap_v2ShardedLayout)
    im = ic.make_data('uint16', IM, 'random', 9)
    fn = tmp_path / 'im.zarr'
    cppzarr.write_zarr(str(fn), im, chunks=[16, 16, 5], shards=[32, 32, 10], zarr_format=3, order='F')
    assert ic.bits_equal(cppzarr.read_zarr(str(fn)), im)


def test_image2dIsStoredWithSingletonZ(tmp_path):
    # the package keeps a 2D array 2D; PetaKit5D's [ny, nx, 1] needs the caller to add the axis
    import json
    im = ic.make_data('uint16', (21, 17), 'random', 10)
    fn2, fn3 = tmp_path / '2d.zarr', tmp_path / '3d.zarr'
    cppzarr.write_zarr(str(fn2), im, chunks=[21, 17], order='F', dimension_separator='.')
    _write(fn3, im)
    assert json.load(open(fn2 / '.zarray'))['shape'] == [21, 17]
    assert json.load(open(fn3 / '.zarray'))['shape'] == [21, 17, 1]
    assert ic.bits_equal(cppzarr.read_zarr(str(fn3)), im[:, :, None])


def test_readBoundingBox(tmp_path):
    # MATLAB's bbox [ymin, xmin, zmin, ymax, xmax, zmax] (1-based, inclusive) is
    # start [ymin-1, xmin-1, zmin-1], end [ymax, xmax, zmax]
    im = ic.make_data('uint16', IM, 'random', 11)
    fn = tmp_path / 'im.zarr'
    _write(fn, im)
    for b in ([3, 5, 2, 30, 33, 11], [17, 17, 6, 17, 17, 6], [1, 1, 7, 37, 41, 7], [1, 1, 1, *IM]):
        out = cppzarr.read_zarr(str(fn), [b[0] - 1, b[1] - 1, b[2] - 1], b[3:])
        assert ic.bits_equal(out, im[b[0] - 1:b[3], b[1] - 1:b[4], b[2] - 1:b[5]]), b


def test_readIntoOut(tmp_path):
    im = ic.make_data('uint16', IM, 'random', 12)
    fn = tmp_path / 'im.zarr'
    _write(fn, im)
    out = np.empty((10, 11, 4), np.uint16, order='F')
    assert cppzarr.read_zarr(str(fn), [2, 3, 4], [12, 14, 8], out=out) is out
    assert np.array_equal(out, im[2:12, 3:14, 4:8])


def test_readBoundingBoxPastTheEndErrors(tmp_path):
    fn = tmp_path / 'im.zarr'
    _write(fn, ic.make_data('uint16', IM, 'random', 13))
    with pytest.raises(Exception):
        cppzarr.read_zarr(str(fn), [0, 0, 0], [IM[0] + 1, IM[1], IM[2]])


def test_readMissingFileErrors(tmp_path):
    with pytest.raises(Exception):
        cppzarr.read_zarr(str(tmp_path / 'missing.zarr'))


def test_partialWritesUpdateOnlyTheRegion(tmp_path):
    im = ic.make_data('uint16', IM, 'random', 14)
    fn = tmp_path / 'im.zarr'
    _write(fn, np.zeros(IM, np.uint16, order='F'))
    expected = np.zeros(IM, np.uint16, order='F')
    assert ic.bits_equal(cppzarr.read_zarr(str(fn)), expected)
    for b, region in (([5, 7, 2, 30, 20, 9], None), ([10, 1, 1, 37, 15, 13], ic.make_data('uint16', (28, 15, 13), 'random', 15))):
        if region is None:
            region = im[b[0] - 1:b[3], b[1] - 1:b[4], b[2] - 1:b[5]]
        cppzarr.write_zarr(str(fn), region, [b[0] - 1, b[1] - 1, b[2] - 1], b[3:])
        expected[b[0] - 1:b[3], b[1] - 1:b[4], b[2] - 1:b[5]] = region
        assert ic.bits_equal(cppzarr.read_zarr(str(fn)), expected), b
    # no temporary files left; names are y.x.z
    assert all(len(n.split('.')) == 3 and all(p.isdigit() for p in n.split('.')) for n in _chunks(fn))


def test_regionWriteIntoMatlabStyleMetadataOnly(tmp_path):
    # the stand-in for createzarr: the .zarray text written by the caller, no chunks
    fn = tmp_path / 'im.zarr'
    fn.mkdir()
    (fn / '.zarray').write_text(MATLAB_ZARRAY)
    assert ic.bits_equal(cppzarr.read_zarr(str(fn)), np.zeros(IM, np.uint16, order='F'))
    region = ic.make_data('uint16', (20, 20, 6), 'random', 16)
    cppzarr.write_zarr(str(fn), region, [3, 4, 5], [23, 24, 11])
    assert np.array_equal(cppzarr.read_zarr(str(fn), [3, 4, 5], [23, 24, 11]), region)
    assert (fn / '.zarray').read_text() == MATLAB_ZARRAY


def test_overwriteReplacesExistingArray(tmp_path):
    fn = tmp_path / 'im.zarr'
    _write(fn, ic.make_data('float32', (20, 20, 20), 'random', 17))
    im = ic.make_data('uint16', IM, 'random', 18)
    _write(fn, im)
    assert ic.bits_equal(cppzarr.read_zarr(str(fn)), im)


def test_sparseWriteSkipsEmptyChunks(tmp_path):
    im = ic.make_data('uint16', (32, 32, 10), 'blocks', 19)
    fn = tmp_path / 'im.zarr'
    _write(fn, im)
    assert ic.bits_equal(cppzarr.read_zarr(str(fn)), im)
    assert _chunks(fn) == ['0.0.0']


def test_zeroRegionRemovesChunkFiles(tmp_path):
    # an all-zero region over written chunks removes their files (reads as the fill value)
    fn = tmp_path / 'im.zarr'
    _write(fn, ic.make_data('uint16', (32, 32, 10), 'smooth', 20))
    assert len(_chunks(fn)) == 8
    cppzarr.write_zarr(str(fn), np.zeros((16, 16, 5), np.uint16), [0, 0, 0], [16, 16, 5])
    assert '0.0.0' not in _chunks(fn) and len(_chunks(fn)) == 7


@pytest.mark.parametrize('layout', ['C', 'F', 'strided'])
def test_writeLayoutsGiveTheSameChunks(tmp_path, layout):
    im = ic.make_data('uint16', IM, 'smooth', 21)
    a = {'C': np.ascontiguousarray(im), 'F': im, 'strided': np.asfortranarray(np.repeat(im, 2, axis=1))[:, ::2]}[layout]
    ref, fn = tmp_path / 'ref.zarr', tmp_path / 'im.zarr'
    _write(ref, im)
    _write(fn, a)
    assert _chunks(ref) == _chunks(fn)
    assert all((ref / c).read_bytes() == (fn / c).read_bytes() for c in _chunks(ref))


def test_threadsGiveTheSameChunks(tmp_path):
    im = ic.make_data('uint16', (300, 280, 40), 'smooth', 22)
    files = []
    for t in (1, 3, 8):
        fn = tmp_path / f'{t}.zarr'
        _write(fn, im, (128, 128, 16), num_threads=t)
        files.append({c: (fn / c).read_bytes() for c in _chunks(fn)})
        assert ic.bits_equal(cppzarr.read_zarr(str(fn), num_threads=t), im)
    assert files[0] == files[1] == files[2]


def test_readBatch(tmp_path):
    im = ic.make_data('uint16', IM, 'random', 23)
    fn = tmp_path / 'im.zarr'
    _write(fn, im)
    out = cppzarr.read_zarr_batch([(str(fn), [0, 0, 0], [10, 10, 10]), (str(fn), None, None)])
    assert np.array_equal(out[0], im[:10, :10, :10]) and np.array_equal(out[1], im)


# --- gaps (each is in the IO1 report with a reproducer) ----------------------------------------

@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='write_zarr always skips all-zero chunks (sparse is not exposed); '
                                       "PetaKit5D's writezarr writes them (sparseData false)")
def test_gap_allZeroChunksCanBeWritten(tmp_path):
    fn = tmp_path / 'im.zarr'
    _write(fn, np.zeros((32, 32, 10), np.uint16, order='F'))
    assert len(_chunks(fn)) == 8


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='no metadata-only creation (createZarrFile) in the Python package')
def test_gap_metadataOnlyCreate():
    assert any(hasattr(cppzarr, n) for n in ('create_zarr', 'create_zarr_file', 'create'))


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason="region writes convert with NumPy's astype (truncate, wrap); the mex "
                                       'converts as MATLAB does (round half away from zero, saturate, NaN to 0)')
def test_gap_regionWriteConvertsAsMatlab(tmp_path):
    fn = tmp_path / 'im.zarr'
    _write(fn, np.zeros((4, 4, 1), np.uint16))
    v = ic.convert_values((4, 4, 1))
    with np.errstate(invalid='ignore'):
        cppzarr.write_zarr(str(fn), v, [0, 0, 0], [4, 4, 1])
    assert np.array_equal(cppzarr.read_zarr(str(fn)), ic.matlab_uint16(v))


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason="uuid=False (chunks written in place, the mex's 'uuid' option) is not exposed")
def test_gap_uuidFalseWritesChunksInPlace(tmp_path):
    fn = tmp_path / 'im.zarr'
    _write(fn, ic.make_data('uint16', (16, 16, 8), 'smooth', 24), (16, 16, 8))
    link = tmp_path / 'link'
    os.link(fn / '0.0.0', link)
    data = ic.make_data('uint16', (16, 16, 8), 'smooth', 25)
    try:
        cppzarr.write_zarr(str(fn), data, order='F', chunks=[16, 16, 8], dimension_separator='.', uuid=False)
    except TypeError:
        pytest.fail('no uuid option')
    assert link.read_bytes() == (fn / '0.0.0').read_bytes()


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason="createzarr with shardSize writes a Zarr v2 .zarray with a 'codecs' "
                                       'sharding_indexed entry; the package shards only Zarr v3 arrays')
def test_gap_v2ShardedLayout(tmp_path):
    fn = tmp_path / 'im.zarr'
    cppzarr.write_zarr(str(fn), ic.make_data('uint16', IM, 'random', 26), chunks=[16, 16, 5], shards=[32, 32, 10],
                       zarr_format=2, order='F', dimension_separator='.')
    assert (fn / '.zarray').is_file()
