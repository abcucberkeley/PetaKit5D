"""petakit5d.io tiff reading and writing (task IO2).

Ports the intent of tests/io/TiffIOTest.m (test names kept), plus the checks of the layout
change, the writer entries and the input checks cpp-tiff 1.7.3 lacks (T3-T9, IO1 report).
"""

import os
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "bindings"))
import iobind_cases as ic  # noqa: E402

cpptiff = pytest.importorskip("cpptiff")

from petakit5d.io import get_image_dtype, get_image_size, read_tiff, write_tiff  # noqa: E402
from petakit5d.io import tiff as T  # noqa: E402
from petakit5d.io.imagesize import _tiff_header  # noqa: E402

DTYPES = ic.TIFF_DTYPES
SHAPES = {"image2d": (37, 53), "stack": (37, 53, 11), "singleRow": (1, 64, 5)}


def _rand(dtype, shape, seed=1):
    return ic.make_data(dtype, shape, "random", seed)


def _compression_tag(fn):
    return _tiff_header(str(fn)).compression


def _write_imagej_stack(im, fn, single_ifd):
    """uint16 stack as ImageJ writes it (port of TiffIOTest.writeImageJStack): big-endian,
    'ImageJ=...\\nimages=N' in the first IFD, slices contiguous, one IFD per slice after the
    data (or only the first IFD, as ImageJ does for stacks over 4 GB)."""
    ny, nx, nz = im.shape
    n_ifd = 1 if single_ifd else nz
    desc = f"ImageJ=1.54f\nimages={nz}\nslices={nz}\nloop=false\n".encode() + b"\0"
    ifd_bytes = [2 + n * 12 + 4 for n in [11] + [10] * (n_ifd - 1)]
    desc_off = 8 + ifd_bytes[0]
    data_off = desc_off + len(desc)
    slice_bytes = ny * nx * 2
    ifd_offs = [8] + [
        data_off + nz * slice_bytes + int(s)
        for s in np.concatenate([[0], np.cumsum(ifd_bytes[1:-1])])
    ]
    pixels = np.ascontiguousarray(im.transpose(2, 0, 1)).astype(">u2").tobytes()
    with open(fn, "wb") as f:
        f.write(b"MM" + struct.pack(">HI", 42, 8))
        for z in range(n_ifd):
            if z == 1:
                f.write(desc + pixels)
            entries = [
                (254, 4, 1, 0),
                (256, 4, 1, nx),
                (257, 4, 1, ny),
                (258, 3, 1, 16),
                (259, 3, 1, 1),
                (262, 3, 1, 1),
                (270, 2, len(desc), desc_off),
                (273, 4, 1, data_off + z * slice_bytes),
                (277, 3, 1, 1),
                (278, 4, 1, ny),
                (279, 4, 1, slice_bytes),
            ]
            if z > 0:
                del entries[6]
            f.write(struct.pack(">H", len(entries)))
            for tag, typ, cnt, val in entries:
                f.write(struct.pack(">HHI", tag, typ, cnt))
                f.write(struct.pack(">HH", val, 0) if typ == 3 else struct.pack(">I", val))
            f.write(struct.pack(">I", ifd_offs[z + 1] if z < n_ifd - 1 else 0))
        if n_ifd == 1:
            f.write(desc + pixels)


# --- TiffIOTest ------------------------------------------------------------------------------


@pytest.mark.parametrize("dtype", DTYPES)
def test_roundTripPreservesTypeAndValues(tmp_path, dtype):
    im = _rand(dtype, (37, 53, 11))
    fn = tmp_path / "im.tif"
    write_tiff(im, fn)
    out = read_tiff(fn)
    assert out.flags.f_contiguous
    assert ic.bits_equal(out, im)
    assert get_image_dtype(fn) == ic.MATLAB_CLASS[dtype]


@pytest.mark.parametrize("name", SHAPES)
def test_roundTripPreservesShape(tmp_path, name):
    im = _rand("uint16", SHAPES[name])
    fn = tmp_path / "im.tif"
    write_tiff(im, fn)
    out = read_tiff(fn)
    assert out.shape == im.shape and out.flags.f_contiguous
    assert ic.bits_equal(out, im)
    assert get_image_size(fn) == (SHAPES[name] + (1,))[:3]


@pytest.mark.parametrize(
    "mode",
    [
        {"compression": "lzw", "mode": "parallel"},
        {"compression": "none", "mode": "parallel"},
        {"mode": "libtiff"},
    ],
    ids=["parallelLzw", "parallelNone", "libtiff"],
)
def test_writeModesAreReadableByAllReaders(tmp_path, mode):
    im = _rand("uint16", (37, 53, 11))
    fn = tmp_path / "im.tif"
    write_tiff(im, fn, **mode)
    assert ic.bits_equal(read_tiff(fn), im)
    tifffile = pytest.importorskip("tifffile")
    if mode.get("compression") != "none":
        pytest.importorskip("imagecodecs")
    assert np.array_equal(tifffile.imread(str(fn)), im.transpose(2, 0, 1))


@pytest.mark.parametrize("compression", ["none", "lzw", "packbits", "deflate"])
def test_readsTiffWrittenByMatlab(tmp_path, compression):
    # classic multi-page TIFF of another libtiff-style writer (MATLAB's imwrite files are
    # read in the oracle cases)
    tifffile = pytest.importorskip("tifffile")
    if compression != "none":
        pytest.importorskip("imagecodecs")
    im = _rand("uint16", (37, 53, 11))
    fn = tmp_path / f"{compression}.tif"
    tifffile.imwrite(
        str(fn),
        np.ascontiguousarray(im.transpose(2, 0, 1)),
        compression=None if compression == "none" else compression,
    )
    assert ic.bits_equal(read_tiff(fn), im)
    assert ic.bits_equal(read_tiff(fn, range=[4, 9]), im[:, :, 3:9])


def test_rangeReadsInclusiveSlices(tmp_path):
    im = _rand("uint16", (37, 53, 11))
    fn = tmp_path / "im.tif"
    write_tiff(im, fn)
    assert ic.bits_equal(read_tiff(fn, range=[3, 7]), im[:, :, 2:7])
    one = read_tiff(fn, range=[5, 5])
    assert one.shape == (37, 53) and ic.bits_equal(one, im[:, :, 4])
    assert ic.bits_equal(read_tiff(fn, range=[1, 11]), im)
    assert ic.bits_equal(read_tiff(fn, range=np.array([3.0, 7.0])), im[:, :, 2:7])
    # readtiff_parallel takes a list of slices
    assert ic.bits_equal(T.read_tiff_parallel(fn, [3, 4, 5, 9, 2]), im[:, :, [2, 3, 4, 8, 1]])


@pytest.mark.parametrize("rng", [[0, 3], [3, 12], [12, 12], [5, 4], [2], [1.5, 3]])
def test_rangeOutOfBoundsErrors(tmp_path, rng):
    # the mex's "Range is out of bounds" (cpp-tiff 1.7.3 reports a damaged file, T7)
    fn = tmp_path / "im.tif"
    write_tiff(_rand("uint16", (5, 6, 11)), fn)
    with pytest.raises(ValueError):
        read_tiff(fn, range=rng)


@pytest.mark.slow
def test_moreThan65535Pages(tmp_path):
    im = _rand("uint8", (4, 4, 70000))
    fn = tmp_path / "pages.tif"
    write_tiff(im, fn)
    assert get_image_size(fn) == (4, 4, 70000)
    assert ic.bits_equal(read_tiff(fn), im)
    assert ic.bits_equal(read_tiff(fn, range=[69990, 70000]), im[:, :, 69989:])


def test_overwritesExistingFile(tmp_path):
    fn = tmp_path / "im.tif"
    write_tiff(_rand("uint16", (20, 30, 6)), fn)
    im = _rand("float32", (10, 12, 3), 2)
    write_tiff(im, fn)
    assert ic.bits_equal(read_tiff(fn), im)


def test_createsMissingOutputFolders(tmp_path):
    fn = tmp_path / "a" / "b" / "c" / "im.tif"
    im = _rand("uint16", (20, 30, 6))
    write_tiff(im, fn)
    assert ic.bits_equal(read_tiff(fn), im)


@pytest.mark.parametrize(
    "single_ifd", [False, True], ids=["readsImageJStack", "readsSingleIfdImageJStack"]
)
def test_readsImageJStack(tmp_path, single_ifd):
    im = _rand("uint16", (20, 30, 9))
    fn = tmp_path / "imagej.tif"
    _write_imagej_stack(im, fn, single_ifd)
    assert ic.bits_equal(read_tiff(fn), im)
    # getImageSizeMex counts the slices of a single-IFD stack (cpp-tiff's shape: 1, T8)
    assert get_image_size(fn) == (20, 30, 9)
    # rangeReadOfImageJStack: the slices of the range
    assert ic.bits_equal(read_tiff(fn, range=[3, 6]), im[:, :, 2:6])
    with pytest.raises(ValueError):
        read_tiff(fn, range=[3, 10])


@pytest.mark.parametrize(
    "bad", ["4-D", "empty", "zeroY", "zeroZ", "bool", "float16", "complex", "1-D"]
)
def test_unsupportedInputIsNotSilentlyDropped(tmp_path, bad):
    # cpp-tiff writes some of these without an error (T3, T4: a read past the buffer)
    arr = {
        "4-D": np.zeros((4, 4, 2, 2), np.uint16),
        "empty": np.zeros((0, 0), np.uint16),
        "zeroY": np.zeros((0, 4, 3), np.uint16),
        "zeroZ": np.zeros((4, 4, 0), np.uint16),
        "bool": np.zeros((4, 4), bool),
        "float16": np.zeros((4, 4), np.float16),
        "complex": np.zeros((4, 4), np.complex64),
        "1-D": np.zeros(5, np.uint16),
    }[bad]
    fn = tmp_path / "x.tif"
    with pytest.raises((ValueError, TypeError)):
        write_tiff(arr, fn)
    assert not fn.exists()


@pytest.mark.parametrize("dtype", DTYPES)
def test_zstdRoundTrip(tmp_path, dtype):
    im = _rand(dtype, (37, 53, 11))
    fn = tmp_path / "im.tif"
    write_tiff(im, fn, compression="zstd")
    assert _compression_tag(fn) == 50000
    assert ic.bits_equal(read_tiff(fn), im)
    assert ic.bits_equal(read_tiff(fn, range=[3, 7]), im[:, :, 2:7])


def test_zstdIsReadByProcessingReaders(tmp_path):
    # (MPageTiffAdapter regions: not ported, blockedImage is replaced by region reads)
    im = _rand("uint16", (37, 53, 11))
    fn = tmp_path / "im.tif"
    write_tiff(im, fn, compression="zstd")
    assert get_image_size(fn) == (37, 53, 11)
    assert get_image_dtype(fn) == "uint16"
    tifffile = pytest.importorskip("tifffile")
    pytest.importorskip("imagecodecs")
    assert np.array_equal(tifffile.imread(str(fn)), im.transpose(2, 0, 1))


def test_compressionOfWriteModes(tmp_path):
    im = _rand("uint16", (37, 53, 4))
    write_tiff(im, tmp_path / "default.tif")
    write_tiff(im, tmp_path / "libtiff.tif", compression="zstd", mode="libtiff")
    write_tiff(im, tmp_path / "none.tif", compression="none")
    assert _compression_tag(tmp_path / "default.tif") == 5
    assert _compression_tag(tmp_path / "libtiff.tif") == 5
    assert _compression_tag(tmp_path / "none.tif") == 1
    assert ic.bits_equal(read_tiff(tmp_path / "libtiff.tif"), im)


# --- the file is the mex's ------------------------------------------------------------------


@pytest.mark.parametrize("compression", ["lzw", "zstd", "none"])
@pytest.mark.parametrize("dtype", ["uint8", "uint16", "float32", "float64", "int64"])
def test_fileIsCpptiffsFile(tmp_path, dtype, compression):
    # byte for byte the file of cpp-tiff's public write_tiff (which IO1 found equal to the
    # mex's), for F, C and strided input
    im = _rand(dtype, (1031, 37, 3), 3)
    ref = tmp_path / "ref.tif"
    cpptiff.write_tiff(
        str(ref), np.ascontiguousarray(im.transpose(2, 0, 1)), compression=compression
    )
    layouts = {
        "F": im,
        "C": np.ascontiguousarray(im),
        "strided": np.asfortranarray(np.repeat(im, 2, axis=1))[:, ::2],
    }
    for name, a in layouts.items():
        fn = tmp_path / f"{name}.tif"
        write_tiff(a, fn, compression=compression)
        assert fn.read_bytes() == ref.read_bytes(), name


def test_2dFileIsCpptiffsFile(tmp_path):
    im = _rand("uint16", (37, 53), 4)
    ref = tmp_path / "ref.tif"
    cpptiff.write_tiff(str(ref), np.ascontiguousarray(im))
    for name, a in {"F": im, "C": np.ascontiguousarray(im)}.items():
        fn = tmp_path / f"{name}.tif"
        write_tiff(a, fn)
        assert fn.read_bytes() == ref.read_bytes(), name


def test_writerEntryIsDetected():
    # the direct entry while cpp-tiff has no order='F'; the copying public write otherwise
    assert T._writer() in ("order", "direct")
    assert T._reader() in ("order", "transpose")


def test_publicWriteFallbackWritesTheSameFile(tmp_path, monkeypatch):
    im = _rand("uint16", (40, 30, 5), 5)
    write_tiff(im, tmp_path / "a.tif")
    monkeypatch.setattr(T, "_writer", lambda: "copy")
    write_tiff(im, tmp_path / "b.tif")
    write_tiff(im[:, :, 0], tmp_path / "c.tif")
    monkeypatch.undo()
    write_tiff(im[:, :, 0], tmp_path / "d.tif")
    assert (tmp_path / "a.tif").read_bytes() == (tmp_path / "b.tif").read_bytes()
    assert (tmp_path / "c.tif").read_bytes() == (tmp_path / "d.tif").read_bytes()


def test_bigEndianAndPathLikeInput(tmp_path):
    # cpp-tiff writes big-endian data unswapped (T6) and takes str names only (T9)
    im = _rand("uint16", (9, 7, 3), 6)
    write_tiff(im.astype(">u2"), tmp_path / "be.tif")
    assert ic.bits_equal(read_tiff(tmp_path / "be.tif"), im)
    write_tiff(im, str(tmp_path / "s.tif"))
    assert ic.bits_equal(read_tiff(str(tmp_path / "s.tif")), im)


@pytest.mark.parametrize("compression", ["LZW", "NONE", "Zstd", "jpeg"])
def test_onlyMatlabsCompressionNames(tmp_path, compression):
    # mustBeMember is case-sensitive; cpp-tiff writes LZW for any other name (T5)
    with pytest.raises(ValueError):
        write_tiff(_rand("uint16", (4, 4)), tmp_path / "x.tif", compression=compression)


def test_missingFileErrors(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_tiff(tmp_path / "missing.tif")


def test_groupWrite(tmp_path):
    if os.name == "nt":
        pytest.skip("no group permissions on Windows")
    fn = tmp_path / "im.tif"
    old = os.umask(0o022)
    try:
        write_tiff(_rand("uint16", (4, 4)), fn, group_write=False)
        assert not os.stat(fn).st_mode & 0o020
        write_tiff(_rand("uint16", (4, 4)), fn)
        assert os.stat(fn).st_mode & 0o020
    finally:
        os.umask(old)


# --- the layout change -----------------------------------------------------------------------


def _kernels():
    from petakit5d import _ext

    out = [("numpy", T._transpose_numpy)]
    if _ext.available("io"):
        ext = _ext.load("io")
        out.append(("ext", ext.transpose_pages))

        for mode in (0, 1, 2):  # chosen, whole pages per thread, shared pages

            def inplace(src, dst, num_threads, mode=mode):
                nx, ny, nz = src.shape
                buf = src.ravel(order="F").copy()
                ext.transpose_pages_inplace(buf, nx, ny, nz, num_threads, mode)
                dst[...] = buf.reshape((ny, nx, nz), order="F")

            out.append((f"inplace{mode}", inplace))
    return out


@pytest.mark.parametrize("num_threads", [1, 3, 8])
@pytest.mark.parametrize(
    "dtype", ["uint8", "int16", "uint16", "float32", "float64", "uint64", "bool"]
)
def test_transposeKernelIsExact(dtype, num_threads):
    # sizes not divisible by the tiles (64) or the SIMD blocks (2, 4, 8); fewer and more
    # pages than threads (the in-place kernel's two ways of splitting the work)
    for shape in [
        (1, 1, 1),
        (7, 5, 3),
        (64, 64, 2),
        (130, 67, 3),
        (9, 200, 1),
        (1031, 3, 2),
        (33, 17, 20),
        (300, 129, 9),
    ]:
        nx, ny, nz = shape
        a = ic.make_data(dtype, (ny, nx, nz), "random", nx * ny)
        src = np.ascontiguousarray(a.transpose(2, 0, 1)).transpose(2, 1, 0)  # (x, y, z) F
        assert src.flags.f_contiguous
        for name, k in _kernels():
            dst = np.empty((ny, nx, nz), a.dtype, order="F")
            k(src, dst, num_threads)
            assert ic.bits_equal(dst, a), (name, shape)


def test_inplaceKernelOfLargePages():
    # pages too large for a scratch page per thread: the threads share each page
    from petakit5d import _ext

    if not _ext.available("io"):
        pytest.skip("petakit5d._ext.io is not built")
    a = ic.make_data("uint16", (2900, 3001, 5), "random", 9)  # 17 MB pages, 87 MB in all
    buf = np.ascontiguousarray(a.transpose(2, 0, 1)).reshape(-1)
    _ext.load("io").transpose_pages_inplace(buf, 3001, 2900, 5, 8)
    assert ic.bits_equal(buf.reshape(a.shape, order="F"), a)


def test_transposeKernelChecksItsArguments():
    from petakit5d import _ext

    if not _ext.available("io"):
        pytest.skip("petakit5d._ext.io is not built")
    k = _ext.load("io").transpose_pages
    src = np.zeros((4, 5, 2), np.uint16, order="F")
    with pytest.raises((TypeError, ValueError)):
        k(src, np.zeros((4, 5, 2), np.uint16, order="F"), 1)  # not (ny, nx, nz)
    with pytest.raises((TypeError, ValueError)):
        k(src, np.zeros((5, 4, 2), np.uint8, order="F"), 1)  # another type
    with pytest.raises((TypeError, ValueError)):
        k(src, np.zeros((5, 4, 2), np.uint16, order="C"), 1)  # C order, not converted
    with pytest.raises((TypeError, ValueError)):
        k(src, np.zeros((5, 4, 2), np.uint16, order="F"), 0)
    ip = _ext.load("io").transpose_pages_inplace
    with pytest.raises((TypeError, ValueError)):
        ip(np.zeros(40, np.uint16), 4, 5, 3, 1)  # not nx * ny * nz elements
    with pytest.raises((TypeError, ValueError)):
        ip(np.zeros(40, np.uint16)[::2], 4, 5, 1, 1)  # not contiguous, not copied
    with pytest.raises((TypeError, ValueError)):
        ip(np.zeros(40, np.uint16), 4, 5, 2, 1, 3)  # no such mode


def test_readIsOneLayoutChange(tmp_path):
    # the result is cpp-tiff's buffer, turned around in place: one volume in memory, as the
    # mex (peak memory of whole tiles: benchmarks)
    import tracemalloc

    im = _rand("uint16", (300, 200, 20), 7)
    fn = tmp_path / "im.tif"
    write_tiff(im, fn, compression="none")
    tracemalloc.start()
    out = read_tiff(fn)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert out.flags.f_contiguous
    assert ic.bits_equal(out, im)
    from petakit5d import _ext

    limit = 1.05 if T._reader() == "order" or _ext.available("io") else 2.05
    assert peak <= im.nbytes * limit


def test_readThreadsGiveTheSameResult(tmp_path):
    im = _rand("uint16", (300, 257, 9), 8)
    fn = tmp_path / "im.tif"
    write_tiff(im, fn)
    for t in (1, 2, 5):
        assert ic.bits_equal(read_tiff(fn, num_threads=t), im)
