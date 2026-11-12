"""petakit5d.io tiff reading and writing (tasks IO2, IO3).

Ports the intent of tests/io/TiffIOTest.m (test names kept) on the port's (z, y, x) arrays,
plus the file identity with MATLAB's files, the input checks cpp-tiff 1.7.3 lacks (T3-T9, IO1
report) and the MATLAB oracle cases.
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

from petakit5d.core.arrays import from_matlab, to_matlab  # noqa: E402
from petakit5d.core.bbox import range_from_matlab  # noqa: E402
from petakit5d.io import get_image_dtype, get_image_size, read_tiff, write_tiff  # noqa: E402
from petakit5d.io import tiff as T  # noqa: E402
from petakit5d.io.imagesize import _tiff_header  # noqa: E402
from petakit5d.testing import oracle  # noqa: E402
from petakit5d.testing.cases import Case, In, Value, Work  # noqa: E402

DTYPES = ic.TIFF_DTYPES
# (z, y, x) shapes, never cubic
SHAPES = {"image2d": (37, 53), "stack": (11, 37, 53), "singleRow": (5, 1, 64)}


def _rand(dtype, shape, seed=1, pattern="random"):
    """a C-ordered array of the shape (values over the type's range, specials for floats)"""
    return np.ascontiguousarray(ic.make_data(dtype, shape, pattern, seed))


def _compression_tag(fn):
    return _tiff_header(str(fn)).compression


def _write_imagej_stack(vol, fn, single_ifd):
    """uint16 (z, y, x) stack as ImageJ writes it (port of TiffIOTest.writeImageJStack):
    big-endian, 'ImageJ=...\\nimages=N' in the first IFD, slices contiguous, one IFD per slice
    after the data (or only the first IFD, as ImageJ does for stacks over 4 GB)."""
    nz, ny, nx = vol.shape
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
    pixels = vol.astype(">u2").tobytes()
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


def _mex_bytes(tmp_path, m, compression="lzw"):
    """the file the mex writes for MATLAB's array m ((y, x, z) Fortran): IO1 found cpp-tiff's
    transpose=True entry, given m's memory, byte-identical to parallelWriteTiff"""
    from cpptiff.cpptiff import pybind11_write_tiff

    m = np.asfortranarray(m)
    fn = tmp_path / "mex.tif"
    pybind11_write_tiff(
        str(fn), m.transpose(1, 0, 2) if m.ndim == 3 else m.T, True, compression, 0, 1
    )
    return fn.read_bytes()


# --- TiffIOTest ------------------------------------------------------------------------------


@pytest.mark.parametrize("dtype", DTYPES)
def test_roundTripPreservesTypeAndValues(tmp_path, dtype):
    vol = _rand(dtype, (11, 37, 53))
    fn = tmp_path / "im.tif"
    write_tiff(vol, fn)
    out = read_tiff(fn)
    assert out.flags.c_contiguous
    assert ic.bits_equal(out, vol)
    assert get_image_dtype(fn) == ic.MATLAB_CLASS[dtype]


@pytest.mark.parametrize("name", SHAPES)
def test_roundTripPreservesShape(tmp_path, name):
    vol = _rand("uint16", SHAPES[name])
    fn = tmp_path / "im.tif"
    write_tiff(vol, fn)
    out = read_tiff(fn)
    assert out.shape == vol.shape and out.flags.c_contiguous
    assert ic.bits_equal(out, vol)
    assert get_image_size(fn) == ((1,) + SHAPES[name])[-3:]


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
    vol = _rand("uint16", (11, 37, 53))
    fn = tmp_path / "im.tif"
    write_tiff(vol, fn, **mode)
    assert ic.bits_equal(read_tiff(fn), vol)
    tifffile = pytest.importorskip("tifffile")
    if mode.get("compression") != "none":
        pytest.importorskip("imagecodecs")
    assert np.array_equal(tifffile.imread(str(fn)), vol)


@pytest.mark.parametrize("compression", ["none", "lzw", "packbits", "deflate"])
def test_readsTiffWrittenByMatlab(tmp_path, compression):
    # classic multi-page TIFF of another libtiff-style writer (MATLAB's own files: oracle cases)
    tifffile = pytest.importorskip("tifffile")
    if compression != "none":
        pytest.importorskip("imagecodecs")
    vol = _rand("uint16", (11, 37, 53))
    fn = tmp_path / f"{compression}.tif"
    tifffile.imwrite(str(fn), vol, compression=None if compression == "none" else compression)
    assert ic.bits_equal(read_tiff(fn), vol)
    assert ic.bits_equal(read_tiff(fn, range=(3, 9)), vol[3:9])


def test_rangeReadsInclusiveSlices(tmp_path):
    # ranges are (start, stop), 0-based, end-exclusive: MATLAB's [3, 7] is (2, 7)
    vol = _rand("uint16", (11, 37, 53))
    fn = tmp_path / "im.tif"
    write_tiff(vol, fn)
    assert ic.bits_equal(read_tiff(fn, range=range_from_matlab([3, 7])), vol[2:7])
    one = read_tiff(fn, range=(4, 5))
    assert one.shape == (37, 53) and ic.bits_equal(one, vol[4])
    assert ic.bits_equal(read_tiff(fn, range=(0, 11)), vol)
    assert ic.bits_equal(read_tiff(fn, range=(6, None)), vol[6:])
    assert ic.bits_equal(read_tiff(fn, range=np.array([2.0, 7.0])), vol[2:7])
    # readtiff_parallel takes a list of slices (0-based)
    assert ic.bits_equal(T.read_tiff_parallel(fn, [2, 3, 4, 8, 1]), vol[[2, 3, 4, 8, 1]])
    assert ic.bits_equal(T.read_tiff_parallel(fn, [5]), vol[5])


@pytest.mark.parametrize("rng", [(-1, 3), (2, 12), (11, 12), (5, 5), (5, 4), (2,), (1.5, 3)])
def test_rangeOutOfBoundsErrors(tmp_path, rng):
    # the mex's "Range is out of bounds" (cpp-tiff 1.7.3 reports a damaged file, T7)
    fn = tmp_path / "im.tif"
    write_tiff(_rand("uint16", (11, 5, 6)), fn)
    with pytest.raises(ValueError):
        read_tiff(fn, range=rng)


@pytest.mark.slow
def test_moreThan65535Pages(tmp_path):
    vol = _rand("uint8", (70000, 4, 5))
    fn = tmp_path / "pages.tif"
    write_tiff(vol, fn)
    assert get_image_size(fn) == (70000, 4, 5)
    assert ic.bits_equal(read_tiff(fn), vol)
    assert ic.bits_equal(read_tiff(fn, range=(69989, 70000)), vol[69989:])


def test_overwritesExistingFile(tmp_path):
    fn = tmp_path / "im.tif"
    write_tiff(_rand("uint16", (6, 20, 30)), fn)
    vol = _rand("float32", (3, 10, 12), 2)
    write_tiff(vol, fn)
    assert ic.bits_equal(read_tiff(fn), vol)


def test_createsMissingOutputFolders(tmp_path):
    fn = tmp_path / "a" / "b" / "c" / "im.tif"
    vol = _rand("uint16", (6, 20, 30))
    write_tiff(vol, fn)
    assert ic.bits_equal(read_tiff(fn), vol)


@pytest.mark.parametrize(
    "single_ifd", [False, True], ids=["readsImageJStack", "readsSingleIfdImageJStack"]
)
def test_readsImageJStack(tmp_path, single_ifd):
    vol = _rand("uint16", (9, 20, 30))
    fn = tmp_path / "imagej.tif"
    _write_imagej_stack(vol, fn, single_ifd)
    assert ic.bits_equal(read_tiff(fn), vol)
    # getImageSizeMex counts the slices of a single-IFD stack (cpp-tiff's shape: 1, T8)
    assert get_image_size(fn) == (9, 20, 30)
    # rangeReadOfImageJStack: the slices of the range
    assert ic.bits_equal(read_tiff(fn, range=(2, 6)), vol[2:6])
    with pytest.raises(ValueError):
        read_tiff(fn, range=(2, 10))


@pytest.mark.parametrize(
    "bad", ["4-D", "empty", "zeroY", "zeroZ", "bool", "float16", "complex", "1-D", "fortran"]
)
def test_unsupportedInputIsNotSilentlyDropped(tmp_path, bad):
    # cpp-tiff writes some of these without an error (T3, T4: a read past the buffer); a
    # Fortran-ordered volume is likely MATLAB's layout, not (z, y, x)
    arr = {
        "4-D": np.zeros((2, 2, 4, 4), np.uint16),
        "empty": np.zeros((0, 0), np.uint16),
        "zeroY": np.zeros((3, 0, 4), np.uint16),
        "zeroZ": np.zeros((0, 4, 4), np.uint16),
        "bool": np.zeros((4, 4), bool),
        "float16": np.zeros((4, 4), np.float16),
        "complex": np.zeros((4, 4), np.complex64),
        "1-D": np.zeros(5, np.uint16),
        "fortran": np.zeros((3, 4, 5), np.uint16, order="F"),
    }[bad]
    fn = tmp_path / "x.tif"
    with pytest.raises((ValueError, TypeError)):
        write_tiff(arr, fn)
    assert not fn.exists()


@pytest.mark.parametrize("dtype", DTYPES)
def test_zstdRoundTrip(tmp_path, dtype):
    vol = _rand(dtype, (11, 37, 53))
    fn = tmp_path / "im.tif"
    write_tiff(vol, fn, compression="zstd")
    assert _compression_tag(fn) == 50000
    assert ic.bits_equal(read_tiff(fn), vol)
    assert ic.bits_equal(read_tiff(fn, range=(2, 7)), vol[2:7])


def test_zstdIsReadByProcessingReaders(tmp_path):
    # (MPageTiffAdapter regions: not ported, blockedImage is replaced by region reads)
    vol = _rand("uint16", (11, 37, 53))
    fn = tmp_path / "im.tif"
    write_tiff(vol, fn, compression="zstd")
    assert get_image_size(fn) == (11, 37, 53)
    assert get_image_dtype(fn) == "uint16"
    tifffile = pytest.importorskip("tifffile")
    pytest.importorskip("imagecodecs")
    assert np.array_equal(tifffile.imread(str(fn)), vol)


def test_compressionOfWriteModes(tmp_path):
    vol = _rand("uint16", (4, 37, 53))
    write_tiff(vol, tmp_path / "default.tif")
    write_tiff(vol, tmp_path / "libtiff.tif", compression="zstd", mode="libtiff")
    write_tiff(vol, tmp_path / "none.tif", compression="none")
    assert _compression_tag(tmp_path / "default.tif") == 5
    assert _compression_tag(tmp_path / "libtiff.tif") == 5
    assert _compression_tag(tmp_path / "none.tif") == 1
    assert ic.bits_equal(read_tiff(tmp_path / "libtiff.tif"), vol)


# --- the file is MATLAB's ----------------------------------------------------------------------


@pytest.mark.parametrize("compression", ["lzw", "zstd", "none"])
@pytest.mark.parametrize("dtype", ["uint8", "uint16", "float32", "float64", "int64"])
def test_fileIsTheMexsFile(tmp_path, dtype, compression):
    # write_tiff(from_matlab(m)) writes the bytes the mex writes for MATLAB's array m; for C
    # input and strided views alike
    m = ic.make_data(dtype, (1031, 37, 3), "random", 3)  # MATLAB's (y, x, z)
    ref = _mex_bytes(tmp_path, m, compression)
    vol = from_matlab(m)
    layouts = {"C": vol, "strided": np.repeat(vol, 2, axis=2)[:, :, ::2]}
    for name, a in layouts.items():
        fn = tmp_path / f"{name}.tif"
        write_tiff(a, fn, compression=compression)
        assert fn.read_bytes() == ref, name
        assert ic.bits_equal(to_matlab(read_tiff(fn)), m)


def test_imageFileIsTheMexsFile(tmp_path):
    m = ic.make_data("uint16", (37, 53), "random", 4)
    write_tiff(from_matlab(m), tmp_path / "py.tif")
    assert (tmp_path / "py.tif").read_bytes() == _mex_bytes(tmp_path, m)


def test_bigEndianAndPathLikeInput(tmp_path):
    # cpp-tiff writes big-endian data unswapped (T6) and takes str names only (T9)
    vol = _rand("uint16", (3, 9, 7), 6)
    write_tiff(vol.astype(">u2"), tmp_path / "be.tif")
    assert ic.bits_equal(read_tiff(tmp_path / "be.tif"), vol)
    write_tiff(vol, str(tmp_path / "s.tif"))
    assert ic.bits_equal(read_tiff(str(tmp_path / "s.tif")), vol)


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


def test_readAndWriteMakeNoCopy(tmp_path):
    # read_tiff returns cpp-tiff's array, write_tiff passes the C-ordered memory: at most one
    # volume in memory either way
    import tracemalloc

    vol = _rand("uint16", (20, 300, 200), 7)
    fn = tmp_path / "im.tif"
    tracemalloc.start()
    write_tiff(vol, fn, compression="none")
    peak_w = tracemalloc.get_traced_memory()[1]
    tracemalloc.reset_peak()
    out = read_tiff(fn)
    peak_r = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert out.flags.c_contiguous and ic.bits_equal(out, vol)
    assert peak_w <= 0.05 * vol.nbytes
    assert peak_r <= 1.05 * vol.nbytes


def test_readThreadsGiveTheSameResult(tmp_path):
    vol = _rand("uint16", (9, 300, 257), 8)
    fn = tmp_path / "im.tif"
    write_tiff(vol, fn)
    for t in (1, 2, 5):
        assert ic.bits_equal(read_tiff(fn, num_threads=t), vol)


# --- MATLAB oracle: files MATLAB writes are read as from_matlab of its arrays -----------------


def _roundtrip(m, comp, rng, work):
    # Python's side of a case: MATLAB's array as ours, written, read back (range given in
    # MATLAB's form until F14's roles convert it), returned in MATLAB's form
    fn = work + "a.tif"
    write_tiff(from_matlab(m), fn, compression=comp)
    r = None if rng is None else range_from_matlab(rng)
    return to_matlab(read_tiff(fn, range=r))


_TIFF_CASES = []
for _dt in DTYPES:
    for _sname, _shape in {"s3d": (37, 53, 11), "s2d": (37, 53), "row": (1, 64, 5)}.items():
        for _comp in ("lzw", "zstd"):
            _rng = [2, _shape[2] - 1] if len(_shape) == 3 else None
            _TIFF_CASES.append(
                Case(
                    f"tiff_{_dt}_{_sname}_{_comp}",
                    code=(
                        "fn = fullfile(work, 'a.tif'); writetiff(in{1}, fn, 'Compression', in{2}); "
                        "if isempty(in{3}), out{1} = readtiff(fn); "
                        "else, out{1} = readtiff(fn, 'range', in{3}); end"
                    ),
                    args=[In("m"), _comp, [] if _rng is None else _rng, Work()],
                    python=_roundtrip,
                    py_args=[In("m"), _comp, _rng, Work()],
                    inputs={
                        "m": Value(
                            lambda dt=_dt, sh=_shape: ic.make_data(dt, sh, "random", ic.seed_of(dt))
                        )
                    },
                    files=True,
                )
            )


@oracle.parametrize(_TIFF_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)
