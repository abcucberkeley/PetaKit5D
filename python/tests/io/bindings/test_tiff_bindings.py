"""cpp-tiff's Python package as PetaKit5D's tiff reader and writer.

Ports the intent of tests/io/TiffIOTest.m (test names kept) for the binding,
plus the layout and gap checks of task IO1. Volumes are (y, x, z) Fortran
order; cpp-tiff takes and returns (z, y, x) C order.
"""

import os
import struct
import tracemalloc

import numpy as np
import pytest

import iobind_cases as ic

cpptiff = pytest.importorskip('cpptiff')

SHAPES = {'image2d': (37, 53), 'stack': (37, 53, 11), 'singleRow': (1, 64, 5), 'singleColumn': (64, 1, 5)}


def _rand(dtype, shape, seed=1):
    return ic.make_data(dtype, shape, 'random', seed)


def _tags(fn):
    """{tag: (type, count, value or offset)} of the first IFD of a little-endian BigTIFF."""
    with open(fn, 'rb') as f:
        head = f.read(16)
        assert head[:4] == b'II+\x00', 'not a little-endian BigTIFF'
        off = struct.unpack('<Q', head[8:16])[0]
        f.seek(off)
        n = struct.unpack('<Q', f.read(8))[0]
        out = {}
        for _ in range(n):
            tag, typ, cnt = struct.unpack('<HHQ', f.read(12))
            raw = f.read(8)
            val = struct.unpack('<H', raw[:2])[0] if typ == 3 and cnt == 1 else struct.unpack('<Q', raw)[0]
            out[tag] = (typ, cnt, val)
    return out


@pytest.mark.parametrize('dtype', ic.TIFF_DTYPES)
def test_roundTripPreservesTypeAndValues(tmp_path, dtype):
    im = _rand(dtype, (37, 53, 11))
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im)
    out = ic.tiff_read_yxz(fn)
    assert out.dtype == im.dtype and out.flags.f_contiguous
    assert ic.bits_equal(out, im)


@pytest.mark.parametrize('name', SHAPES)
def test_roundTripPreservesShape(tmp_path, name):
    shape = SHAPES[name]
    im = _rand('uint16', shape)
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im)
    assert ic.bits_equal_m(ic.tiff_read_yxz(fn), im)
    nz = shape[2] if len(shape) == 3 else 1
    # get_image_shape is (z, y, x), also for one page
    assert tuple(cpptiff.get_image_shape(str(fn))) == (nz, shape[0], shape[1])


@pytest.mark.parametrize('compression', ic.TIFF_COMPRESSIONS)
@pytest.mark.parametrize('dtype', ['uint8', 'uint16', 'float32', 'float64', 'int64'])
def test_directWriteMatchesPublicWrite(tmp_path, dtype, compression):
    # the extension's transpose=True entry takes the Fortran (y, x, z) memory as
    # the mex does, and writes the same bytes as write_tiff of the (z, y, x) copy
    for shape in ((37, 53, 11), (37, 53), (1031, 517, 2)):
        im = _rand(dtype, shape)
        a, b = tmp_path / 'public.tif', tmp_path / 'direct.tif'
        ic.tiff_write_yxz(a, im, compression)
        ic.tiff_write_yxz_direct(b, im, compression)
        assert a.read_bytes() == b.read_bytes(), shape


def test_inputLayoutsGiveTheSameFile(tmp_path):
    im = _rand('uint16', (40, 30, 6))
    zyx = np.ascontiguousarray(im.transpose(2, 0, 1))
    variants = {'C': zyx, 'F': np.asfortranarray(zyx), 'yxzview': im.transpose(2, 0, 1),
                'strided': np.ascontiguousarray(np.repeat(zyx, 2, axis=2))[:, :, ::2]}
    ref = None
    for name, arr in variants.items():
        fn = tmp_path / f'{name}.tif'
        cpptiff.write_tiff(str(fn), arr)
        ref = ref or fn.read_bytes()
        assert fn.read_bytes() == ref, name


def test_fileLayoutIsBigTiffLzwWithPetaKitStrips(tmp_path):
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, _rand('uint16', (1031, 53, 3)))
    t = _tags(fn)
    # width, length, bits, compression, photometric, strip offsets, samples, rows per strip,
    # strip byte counts, planar configuration, sample format: no other tags
    assert sorted(t) == [256, 257, 258, 259, 262, 273, 277, 278, 279, 284, 339]
    assert t[259][2] == 5 and t[278][2] == 512 and t[273][1] == 3


@pytest.mark.parametrize('compression,tag', [('lzw', 5), ('zstd', 50000), ('none', 1)])
def test_compressionOfWriteModes(tmp_path, compression, tag):
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, _rand('uint16', (37, 53, 4)), compression)
    assert _tags(fn)[259][2] == tag
    fn2 = tmp_path / 'default.tif'
    cpptiff.write_tiff(str(fn2), _rand('uint16', (4, 37, 53)))
    assert _tags(fn2)[259][2] == 5


@pytest.mark.parametrize('dtype', ic.TIFF_DTYPES)
def test_zstdRoundTrip(tmp_path, dtype):
    im = _rand(dtype, (37, 53, 11))
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im, 'zstd')
    assert ic.bits_equal(ic.tiff_read_yxz(fn), im)
    assert ic.bits_equal(ic.tiff_read_yxz(fn, [2, 7]), im[:, :, 2:7])


def test_rangeReadsInclusiveSlices(tmp_path):
    # MATLAB's range [first, last] (1-based, inclusive) is z_range [first - 1, last]
    im = _rand('uint16', (37, 53, 11))
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im)
    assert ic.bits_equal(ic.tiff_read_yxz(fn, [2, 7]), im[:, :, 2:7])
    assert ic.bits_equal(ic.tiff_read_yxz(fn, [0, 11]), im)
    # one slice comes back 2D: (y, x), as MATLAB's size of a single slice
    assert ic.tiff_read_yxz(fn, [4, 5]).shape == (37, 53)
    assert ic.bits_equal(ic.tiff_read_yxz(fn, 4), im[:, :, 4])
    with pytest.raises(Exception):
        cpptiff.read_tiff(str(fn), [5, 5])
    with pytest.raises(Exception):
        cpptiff.read_tiff(str(fn), [3, 2])


def test_rangePastTheEndErrors(tmp_path):
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, _rand('uint16', (5, 7, 6)))
    for zr in ([0, 7], [6, 7], 6):
        with pytest.raises(Exception):
            cpptiff.read_tiff(str(fn), zr)


def test_moreThan65535Pages(tmp_path):
    im = _rand('uint8', (4, 4, 70000))
    fn = tmp_path / 'pages.tif'
    ic.tiff_write_yxz(fn, im)
    assert tuple(cpptiff.get_image_shape(str(fn))) == (70000, 4, 4)
    assert ic.bits_equal(ic.tiff_read_yxz(fn), im)
    assert ic.bits_equal(ic.tiff_read_yxz(fn, [69989, 70000]), im[:, :, 69989:70000])


def test_overwritesExistingFile(tmp_path):
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, _rand('uint16', (20, 30, 6)))
    im = _rand('float32', (10, 12, 3))
    ic.tiff_write_yxz(fn, im)
    assert ic.bits_equal(ic.tiff_read_yxz(fn), im)


def test_createsMissingOutputFolders(tmp_path):
    fn = tmp_path / 'a' / 'b' / 'c' / 'im.tif'
    im = _rand('uint16', (20, 30, 6))
    ic.tiff_write_yxz(fn, im)
    assert ic.bits_equal(ic.tiff_read_yxz(fn), im)


def test_readIntoOut(tmp_path):
    im = _rand('uint16', (20, 30, 6))
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im)
    out = np.empty((6, 20, 30), np.uint16)
    assert cpptiff.read_tiff(str(fn), out=out) is out
    assert np.array_equal(out, im.transpose(2, 0, 1))


@pytest.mark.parametrize('threads', [2, 3, 8, 20])
def test_threadedCopyToYxzEqualsOneThread(tmp_path, threads):
    im = _rand('float32', (33, 21, 13))
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im)
    out = ic.tiff_read_yxz(fn, convert_threads=threads)
    assert out.flags.f_contiguous and ic.bits_equal(out, im)


def test_threadsGiveTheSameFile(tmp_path):
    im = ic.make_data('uint16', (1031, 517, 5), 'smooth', 3)
    files = []
    for t in (1, 2, 7):
        fn = tmp_path / f'{t}.tif'
        ic.tiff_write_yxz(fn, im, num_threads=t)
        files.append(fn.read_bytes())
        assert ic.bits_equal(ic.tiff_read_yxz(fn, num_threads=t), im)
    assert files[0] == files[1] == files[2]


def test_unsupportedInputIsNotSilentlyDropped(tmp_path):
    # (4D data is color, (z, y, x, c) with 3 or 4 samples: 5 samples is not)
    for name, arr in {'4-D': np.zeros((2, 2, 4, 5), np.uint16), 'bool': np.zeros((2, 4, 4), bool),
                      'complex': np.zeros((2, 4, 4), np.complex64)}.items():
        fn = tmp_path / f'{name}.tif'
        with pytest.raises(Exception):
            cpptiff.write_tiff(str(fn), arr)


def _write_imagej_stack(im, fn, single_ifd):
    """uint16 stack as ImageJ writes it (port of TiffIOTest.writeImageJStack): big-endian,
    'ImageJ=...\\nimages=N' in the first IFD, all slices contiguous, one IFD per slice
    after the data (or only the first IFD, as ImageJ does for stacks over 4 GB)."""
    ny, nx, nz = im.shape
    n_ifd = 1 if single_ifd else nz
    desc = f'ImageJ=1.54f\nimages={nz}\nslices={nz}\nloop=false\n'.encode() + b'\0'
    n_entries = [11] + [10] * (n_ifd - 1)
    ifd_bytes = [2 + n * 12 + 4 for n in n_entries]
    desc_off = 8 + ifd_bytes[0]
    data_off = desc_off + len(desc)
    slice_bytes = ny * nx * 2
    ifd_offs = [8] + [data_off + nz * slice_bytes + s for s in np.concatenate([[0], np.cumsum(ifd_bytes[1:-1])]).astype(int)]
    pixels = np.ascontiguousarray(im.transpose(2, 0, 1)).astype('>u2').tobytes()
    with open(fn, 'wb') as f:
        f.write(b'MM' + struct.pack('>HI', 42, 8))
        for z in range(n_ifd):
            if z == 1:
                f.write(desc + pixels)
            entries = [(254, 4, 1, 0), (256, 4, 1, nx), (257, 4, 1, ny), (258, 3, 1, 16), (259, 3, 1, 1), (262, 3, 1, 1),
                       (270, 2, len(desc), desc_off), (273, 4, 1, data_off + z * slice_bytes), (277, 3, 1, 1),
                       (278, 4, 1, ny), (279, 4, 1, slice_bytes)]
            if z > 0:
                del entries[6]
            f.write(struct.pack('>H', len(entries)))
            for tag, typ, cnt, val in entries:
                f.write(struct.pack('>HHI', tag, typ, cnt))
                f.write(struct.pack('>HH', val, 0) if typ == 3 else struct.pack('>I', val))
            f.write(struct.pack('>I', ifd_offs[z + 1] if z < n_ifd - 1 else 0))
        if n_ifd == 1:
            f.write(desc + pixels)


@pytest.mark.parametrize('single_ifd', [False, True], ids=['readsImageJStack', 'readsSingleIfdImageJStack'])
def test_readsImageJStack(tmp_path, single_ifd):
    im = _rand('uint16', (20, 30, 9))
    fn = tmp_path / 'imagej.tif'
    _write_imagej_stack(im, fn, single_ifd)
    assert ic.bits_equal(ic.tiff_read_yxz(fn), im)
    if not single_ifd:
        # (single-IFD stacks: test_gap_imageShapeOfSingleIfdImageJStack)
        assert tuple(cpptiff.get_image_shape(str(fn))) == (9, 20, 30)
    # rangeReadOfImageJStack: the slices of the range, not the first ones
    assert ic.bits_equal(ic.tiff_read_yxz(fn, [2, 6]), im[:, :, 2:6])


@pytest.mark.parametrize('compression', ['none', 'lzw', 'packbits', 'deflate'])
def test_readsTiffWrittenByLibtiffWriters(tmp_path, compression):
    # readsTiffWrittenByMatlab: multi-page classic TIFF from another libtiff-style
    # writer (tifffile here; MATLAB's imwrite files are read in the oracle matrix)
    tifffile = pytest.importorskip('tifffile')
    if compression != 'none':
        pytest.importorskip('imagecodecs')
    im = _rand('uint16', (37, 53, 11))
    fn = tmp_path / f'{compression}.tif'
    tifffile.imwrite(str(fn), np.ascontiguousarray(im.transpose(2, 0, 1)), compression=None if compression == 'none' else compression)
    assert ic.bits_equal(ic.tiff_read_yxz(fn), im)
    assert ic.bits_equal(ic.tiff_read_yxz(fn, [3, 9]), im[:, :, 3:9])


def test_zstdIsReadByOtherReaders(tmp_path):
    tifffile = pytest.importorskip('tifffile')
    pytest.importorskip('imagecodecs')
    im = _rand('uint16', (37, 53, 11))
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im, 'zstd')
    assert np.array_equal(tifffile.imread(str(fn)), im.transpose(2, 0, 1))


# --- gaps (each is in the IO1 report with a reproducer) ----------------------------------------

@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='read_tiff cannot fill a (y, x, z) Fortran array as the mex does; '
                                       'out must be (z, y, x) C order, so every read costs a transposing copy')
def test_gap_readIntoYxzFortranOut(tmp_path):
    im = _rand('uint16', (20, 30, 6))
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, im)
    out = np.empty((20, 30, 6), np.uint16, order='F')
    cpptiff.read_tiff(str(fn), out=out)
    assert np.array_equal(out, im)


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='write_tiff copies a (y, x, z) Fortran volume to (z, y, x) C order '
                                       'before writing; the extension could take it directly (transpose=True)')
def test_gap_publicWriteOfYxzWithoutCopy(tmp_path):
    im = _rand('uint16', (512, 256, 64))
    tracemalloc.start()
    try:
        ic.tiff_write_yxz(tmp_path / 'im.tif', im)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert peak < im.nbytes / 4


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='float16 is written as 16-bit IEEE floats but read back as uint16')
def test_gap_float16RoundTrip(tmp_path):
    im = np.arange(60, dtype=np.float16).reshape(3, 4, 5)
    fn = tmp_path / 'f16.tif'
    cpptiff.write_tiff(str(fn), im)
    assert cpptiff.read_tiff(str(fn)).dtype == np.float16


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='an empty array is written without an error, and reading the file '
                                       'returns uninitialized memory')
def test_gap_emptyArrayIsRejected(tmp_path):
    with pytest.raises(Exception):
        cpptiff.write_tiff(str(tmp_path / 'e.tif'), np.zeros((0, 4, 5), np.uint16))


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason="an unknown or upper-case compression ('jpeg', 'ZSTD') silently writes LZW")
def test_gap_unknownCompressionIsRejected(tmp_path):
    with pytest.raises(Exception):
        cpptiff.write_tiff(str(tmp_path / 'q.tif'), np.zeros((2, 4, 5), np.uint16), compression='jpeg')


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='big-endian input is written without conversion (wrong values)')
def test_gap_bigEndianInput(tmp_path):
    im = np.arange(40, dtype=np.uint16).reshape(2, 4, 5)
    fn = tmp_path / 'be.tif'
    cpptiff.write_tiff(str(fn), im.astype('>u2'))
    assert np.array_equal(cpptiff.read_tiff(str(fn)), im)


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='a z_range past the last page fails with "the file is damaged or uses an '
                                       'unsupported layout" and prints retries for directory 0 on stderr')
def test_gap_rangePastTheEndMessage(tmp_path):
    fn = tmp_path / 'im.tif'
    ic.tiff_write_yxz(fn, _rand('uint16', (5, 7, 6)))
    # (the message, not the file name, which holds the test's name)
    with pytest.raises(Exception, match='(?i)out of bounds|z_range|past the last|has only'):
        cpptiff.read_tiff(str(fn), [0, 7])


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='get_image_shape reports 1 page for an ImageJ stack with one IFD (over 4 GB), '
                                       'which read_tiff reads in full and getImageSizeMex sizes right')
def test_gap_imageShapeOfSingleIfdImageJStack(tmp_path):
    fn = tmp_path / 'imagej.tif'
    _write_imagej_stack(_rand('uint16', (20, 30, 9)), fn, True)
    assert tuple(cpptiff.get_image_shape(str(fn))) == (9, 20, 30)


@pytest.mark.gap
@pytest.mark.xfail(strict=True, reason='read_tiff, write_tiff and get_image_shape take str only (os.PathLike fails '
                                       'in the extension with a TypeError)')
def test_gap_pathLikeNames(tmp_path):
    fn = tmp_path / 'im.tif'
    cpptiff.write_tiff(fn, np.zeros((2, 4, 5), np.uint16))
    cpptiff.read_tiff(fn)
