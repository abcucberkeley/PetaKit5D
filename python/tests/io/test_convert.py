"""Tiff and zarr conversion (petakit5d.io.convert, task T01).

Oracle cases against MATLAB: tiffToZarr (uint8, uint16, single with NaN payloads and -0, double;
chunk sizes that don't divide the image; the default chunks; flip, inputBbox, tileOutBbox and
their combinations; a user function; reading by z blocks with and without boxes; partial files;
images with and without expand2dDim; frames in memory; existing results; another compressor;
the '/' separator above 10,000 chunks; an unsupported class; resampling, which waits for T04),
zarrToTiff (both layouts' classes, an image, a user function, a missing folder, an existing
result), XR_tiffToZarr_wrapper (channel patterns, boxes, reading by blocks, flipped tiles,
partial files, boxes per tile location, files given by path or folder) and
XR_zarrToTiff_wrapper (channel patterns as text and as expressions). Inputs are tiffs and zarr
arrays in MATLAB's layout, read by both sides; the zarr results are compared by data and
metadata after the layout mapping, tiffs byte for byte, parameter files by content.

Then TiffZarrConversionTest (tests/tools) in intent on the port's (z, y, x) data and native
zarr arrays, the job settings of the wrappers (framework spy), the publish order, and the
user-function and resampling hooks.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.core.dtypes import NAN, cast
from petakit5d.io import (
    convert,
    get_zarr_info,
    read_tiff,
    read_zarr,
    write_tiff,
    write_zarr,
)
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, Work
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG

# ---------------------------------------------------------------------------------------------
# data: (z, y, x) volumes from counter-based random numbers (bit-stable everywhere)
# ---------------------------------------------------------------------------------------------

NAN32 = np.float32(NAN)
NAN32_B = np.array([0x7FC00001], dtype=np.uint32).view(np.float32)[0]
SHAPE = (13, 41, 37)  # MATLAB 41 x 37 x 13
CHUNKS = [16, 16, 4]  # MATLAB's [y x z]: remainders on every axis
BOX = [5, 7, 3, 30, 33, 10]  # MATLAB's [ymin xmin zmin ymax xmax zmax]
OUT_BOX = [3, 4, 2, 20, 30, 7]


def volume(shape, cls, name, special=False):
    """a deterministic (z, y, x) volume: integers over the class's range (uint16 below 60000),
    floats normal * 100; with special: NaN of two payloads and -0"""
    rng = CounterRNG(0, f"t01/{name}")
    shape = tuple(int(s) for s in shape)
    if cls == "uint8":
        return rng.integers(0, 255, shape, dtype=np.uint8)
    if cls == "uint16":
        return rng.integers(0, 60000, shape, dtype=np.uint16)
    if cls == "int16":
        return (rng.integers(0, 60000, shape, dtype=np.uint16).astype(np.int32) - 30000).astype(
            np.int16
        )
    a = rng.normal(shape) * 100.0
    a = np.ascontiguousarray(a.astype(np.float32) if cls == "single" else a)
    if special:
        u = rng.uniform((3, *shape))
        a[u[0] < 0.05] = a.dtype.type(NAN)
        a[u[1] < 0.05] = -0.0
        if a.dtype == np.float32:
            a[u[2] < 0.05] = NAN32_B
    return a


def tif_tree(files):
    """a folder of tiff files: {name: (shape, class, special)}"""

    def make(folder):
        for name, (shape, cls, special) in files.items():
            a = volume(shape, cls, name, special)
            synth.write_tiff(Path(folder) / name, a[0] if len(shape) == 3 and shape[0] == 1 else a)

    return make


def zarr_tree(files):
    """a folder of zarr arrays in MATLAB's layout: {name: (shape, class, chunks (z, y, x))}"""

    def make(folder):
        for name, (shape, cls, chunks) in files.items():
            a = volume(shape, cls, name, cls == "single")
            synth.write_zarr(Path(folder) / name, a, block_size=chunks, layout="matlab")

    return make


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

# MATLAB's user functions and the Python functions the port's cases give for them
_FCNS = {
    "@(x) single(x) / 3": lambda x: (x.astype(np.float64) / 3.0).astype(np.float32),
    "@(x) double(x) * 1.7 - 20": lambda x: x.astype(np.float64) * 1.7 - 20.0,
    "@(x) x(:, end : -1 : 1, :)": lambda x: np.ascontiguousarray(x[..., ::-1]),
    "@(x) single(x) * 2": lambda x: x.astype(np.float32) * np.float32(2),
}


def _with_functions(args):
    """the MATLAB arguments with each usrFcn code string as its Python function"""
    out = list(args)
    for i in range(len(out) - 1):
        if isinstance(out[i], str) and out[i].lower() == "usrfcn" and out[i + 1] in _FCNS:
            out[i + 1] = _FCNS[out[i + 1]]
    return out


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*_with_functions(args))
    return []


def py_t2z(*args):
    return call_matlab("tiffToZarr", *args)


def py_z2t(*args):
    return call_matlab("zarrToTiff", *args)


def py_t2z_wrapper(*args):
    if _named(args, "resampleFactor"):
        needs_resampling()
    return call_matlab("XR_tiffToZarr_wrapper", *args)


def py_z2t_wrapper(*args):
    return call_matlab("XR_zarrToTiff_wrapper", *args)


def _named(args, name):
    return any(isinstance(a, str) and a == name for a in args)


def needs_resampling():
    """skips the Python side of a case until T04's imresize3_fast is there (the MATLAB
    reference is made already)"""
    try:
        from petakit5d.tools.resample import imresize3_fast  # noqa: F401
    except ImportError:
        pytest.skip("resampling arrives with T04 (petakit5d.tools.resample.imresize3_fast)")


def _py_t2z_file(d, name, *args):
    if _named(args, "resampleFactor"):
        needs_resampling()
    return py_t2z(d + name, *args)


def _py_t2z_parts(d, names, *args):
    return py_t2z([d + n for n in names], *args)


def _py_t2z_frame(*args):
    return py_t2z("", *args)


def t2z_case(cid, files, src, *args, out="out.zarr", frame=None, error=False):
    """tiffToZarr(<d>/src, <work>/out, [], args...), src a name or a list of names (partial
    files); with frame (src None): tiffToZarr('', <work>/out, frame, args...)"""
    inputs = {"d": Tree(tif_tree(files), "d")}
    if frame is not None:
        inputs["im"] = Value(frame)
        code = "tiffToZarr('', in{:}); out = {};"
        cargs, python = [Work(out), In("im"), *args], _py_t2z_frame
    elif isinstance(src, list):
        code = "d = in{1}; tiffToZarr(cellfun(@(x) [d, x], in{2}, 'unif', 0), in{3 : end});"
        code += " out = {};"
        cargs, python = [InPath("d"), src, Work(out), None, *args], _py_t2z_parts
    else:
        code = "tiffToZarr([in{1}, in{2}], in{3 : end}); out = {};"
        cargs, python = [InPath("d"), src, Work(out), None, *args], _py_t2z_file
    return Case(
        cid,
        code=code,
        nargout=0,
        args=cargs,
        python=python,
        inputs=inputs,
        expect_error=error,
        py_error=(convert.ConvertError, ValueError, OSError),
    )


def _t2z_cases():
    one = {"im.tif": (SHAPE, "uint16", False)}
    bs = ["blockSize", CHUNKS]
    cases = []
    for cls in ("uint8", "uint16", "single", "double"):
        files = {"im.tif": (SHAPE, cls, cls in ("single", "double"))}
        cases.append(t2z_case(f"t2z.{cls}", files, "im.tif", *bs))
    cases += [
        t2z_case("t2z.default_chunks", one, "im.tif"),
        t2z_case("t2z.flip", one, "im.tif", *bs, "flipZstack", True),
        t2z_case("t2z.input_bbox", one, "im.tif", *bs, "inputBbox", BOX),
        t2z_case("t2z.input_bbox_one_slice", one, "im.tif", *bs, "inputBbox", [5, 7, 3, 30, 33, 3]),
        t2z_case("t2z.tile_out_bbox", one, "im.tif", *bs, "tileOutBbox", OUT_BOX),
        t2z_case(
            "t2z.flip_both_boxes",
            one,
            "im.tif",
            *bs,
            "flipZstack",
            True,
            "inputBbox",
            BOX,
            "tileOutBbox",
            [2, 3, 1, 20, 25, 6],
        ),
        t2z_case(
            "t2z.usrfcn_single",
            one,
            "im.tif",
            *bs,
            "flipZstack",
            True,
            "tileOutBbox",
            OUT_BOX,
            "usrFcn",
            "@(x) single(x) / 3",
        ),
        t2z_case(
            "t2z.usrfcn_uint8_saturates",
            {"im.tif": (SHAPE, "uint8", False)},
            "im.tif",
            *bs,
            "usrFcn",
            "@(x) double(x) * 1.7 - 20",
        ),
        t2z_case("t2z.regions", one, "im.tif", *bs, "readWholeTiff", False),
        t2z_case(
            "t2z.regions_input_bbox", one, "im.tif", *bs, "readWholeTiff", False, "inputBbox", BOX
        ),
        t2z_case(
            "t2z.regions_tile_out_bbox",
            one,
            "im.tif",
            *bs,
            "readWholeTiff",
            False,
            "tileOutBbox",
            OUT_BOX,
        ),
        t2z_case(
            "t2z.regions_both_boxes",
            one,
            "im.tif",
            *bs,
            "readWholeTiff",
            False,
            "inputBbox",
            BOX,
            "tileOutBbox",
            [2, 3, 1, 20, 25, 6],
        ),
        t2z_case(
            "t2z.regions_single_one_chunk_layer",
            {"im.tif": (SHAPE, "single", True)},
            "im.tif",
            "blockSize",
            [16, 16, 32],
            "readWholeTiff",
            False,
        ),
        t2z_case(
            "t2z.partial_files",
            {"p1.tif": ((5, 41, 37), "uint16", False), "p2.tif": ((8, 41, 37), "uint16", False)},
            ["p1.tif", "p2.tif"],
            *bs,
        ),
        t2z_case(
            "t2z.partial_files_input_bbox",
            {"p1.tif": ((5, 41, 37), "uint16", False), "p2.tif": ((8, 41, 37), "uint16", False)},
            ["p1.tif", "p2.tif"],
            *bs,
            "inputBbox",
            BOX,
        ),
        t2z_case("t2z.image", {"im.tif": ((1, 41, 37), "uint16", False)}, "im.tif", *bs),
        t2z_case(
            "t2z.image_not_expanded",
            {"im.tif": ((1, 41, 37), "uint16", False)},
            "im.tif",
            *bs,
            "expand2dDim",
            False,
        ),
        t2z_case(
            "t2z.image_wide_not_expanded",
            {"im.tif": ((1, 23, 61), "single", True)},
            "im.tif",
            "expand2dDim",
            False,
        ),
        t2z_case(
            "t2z.frame_flip",
            {},
            None,
            *bs,
            "flipZstack",
            True,
            frame=lambda: to_matlab(volume(SHAPE, "uint16", "frame")),
        ),
        t2z_case(
            "t2z.frame_single_tile_out_bbox",
            {},
            None,
            *bs,
            "tileOutBbox",
            OUT_BOX,
            frame=lambda: to_matlab(volume(SHAPE, "single", "frame_single", True)),
        ),
        t2z_case(
            "t2z.frame_image_not_expanded",
            {},
            None,
            *bs,
            "expand2dDim",
            False,
            frame=lambda: to_matlab(volume((1, 29, 33), "uint16", "frame_2d")),
        ),
        t2z_case("t2z.compressor_lz4", one, "im.tif", *bs, "compressor", "lz4", out="lz4.zarr"),
        t2z_case(
            "t2z.slash_separator",
            {"im.tif": ((40, 25, 23), "uint8", False)},
            "im.tif",
            "blockSize",
            [1, 1, 2],
        ),
        t2z_case(
            "t2z.int16_unsupported",
            {"im.tif": (SHAPE, "int16", False)},
            "im.tif",
            *bs,
            error=True,
        ),
        t2z_case(
            "t2z.resample",
            one,
            "im.tif",
            *bs,
            "resampleFactor",
            [1, 1, 2],
        ),
    ]

    # an existing result: kept, or replaced with Overwrite
    def existing(folder):
        tif_tree(one)(folder)
        synth.write_zarr(
            Path(folder) / "old.zarr", volume((3, 5, 7), "uint16", "old"), (3, 5, 7), "matlab"
        )

    for cid, extra in (
        ("t2z.existing_kept", []),
        ("t2z.existing_overwritten", ["Overwrite", True]),
    ):
        cases.append(
            Case(
                cid,
                code="tiffToZarr([in{1}, 'im.tif'], [in{1}, 'old.zarr'], in{2 : end}); out = {};",
                nargout=0,
                args=[InPath("d"), None, *bs, *extra],
                python=lambda d, *a: py_t2z(d + "im.tif", d + "old.zarr", *a),
                inputs={"d": Tree(existing, "d")},
            )
        )
    return cases


# MATLAB's zarrToTiff with a user function fails on getIORegion of the applied image's InMemory
# adapter and falls back to blockedImage's write with MPageTiffAdapter (MATLAB's Tiff: classic
# TIFF, LZW, a Software tag); the port writes the mex's BigTIFF of the same voxels, a proposed
# listed difference like writetiff's 'none' through saveastiff (decision 0004)
FCN_TIFF = tree(tiff_bytes=False)


def z2t_case(cid, files, src, *args, out="out.tif"):
    fcn = _named(args, "usrFcn")
    return Case(
        cid,
        code="zarrToTiff([in{1}, in{2}], in{3 : end}); out = {};",
        nargout=0,
        args=[InPath("d"), src, Work(out), *args],
        python=lambda d, n, *a: py_z2t(d + n, *a),
        inputs={"d": Tree(zarr_tree(files), "d")},
        tree=FCN_TIFF if fcn else tree(),
        note="MATLAB writes classic TIFF with a user function" if fcn else "",
    )


def _z2t_cases():
    cases = [
        z2t_case(f"z2t.{cls}", {"im.zarr": (SHAPE, cls, (4, 16, 16))}, "im.zarr")
        for cls in ("uint8", "uint16", "single")
    ]
    cases += [
        z2t_case("z2t.image", {"im.zarr": ((1, 41, 37), "uint16", (1, 16, 16))}, "im.zarr"),
        z2t_case(
            "z2t.usrfcn_flip_x",
            {"im.zarr": (SHAPE, "uint16", (4, 16, 16))},
            "im.zarr",
            "usrFcn",
            "@(x) x(:, end : -1 : 1, :)",
        ),
        z2t_case(
            "z2t.usrfcn_to_single",
            {"im.zarr": (SHAPE, "uint16", (4, 16, 16))},
            "im.zarr",
            "usrFcn",
            "@(x) single(x) * 2",
        ),
        z2t_case(
            "z2t.missing_folder",
            {"im.zarr": (SHAPE, "uint16", (4, 16, 16))},
            "im.zarr",
            out="a/b/out.tif",
        ),
    ]

    def existing(folder):
        zarr_tree({"im.zarr": (SHAPE, "uint16", (4, 16, 16))})(folder)
        synth.write_tiff(Path(folder) / "old.tif", volume((3, 5, 7), "uint16", "old"))

    for cid, extra in (
        ("z2t.existing_kept", []),
        ("z2t.existing_overwritten", ["Overwrite", True]),
    ):
        cases.append(
            Case(
                cid,
                code="zarrToTiff([in{1}, 'im.zarr'], [in{1}, 'old.tif'], in{2 : end}); out = {};",
                nargout=0,
                args=[InPath("d"), *extra],
                python=lambda d, *a: py_z2t(d + "im.zarr", d + "old.tif", *a),
                inputs={"d": Tree(existing, "d")},
            )
        )
    return cases


DATASET = ("scan_CamA_ch0_t0000", "scan_CamB_ch0_t0000", "scan_CamA_ch1_t0000")
LOCAL = ["parseCluster", False, "masterCompute", True]


def dataset_tree(fmt, shape=(12, 40, 36), extra=()):
    """TiffZarrConversionTest's writeDataset: three files, the third matching no channel
    pattern of the tests; extra: more tiff names"""

    def make(folder):
        for k, n in enumerate((*DATASET, *extra)):
            a = volume(shape, "uint16", f"ds/{fmt}/{k}")
            if fmt == "zarr":
                synth.write_zarr(Path(folder) / f"{n}.zarr", a, (4, 16, 16), "matlab")
            else:
                synth.write_tiff(Path(folder) / f"{n}.tif", a)

    return make


def t2zw_case(cid, make, *args):
    return Case(
        cid,
        code="XR_tiffToZarr_wrapper(in{:}); out = {};",
        nargout=0,
        args=[InPath("d"), *args, *LOCAL],
        python=py_t2z_wrapper,
        inputs={"d": Tree(make, "d")},
    )


def _t2z_wrapper_cases():
    chans = ["channelPatterns", ["CamA_ch0", "CamB_ch0"]]
    bs = ["blockSize", CHUNKS]
    cases = [
        t2zw_case("t2zw.default", dataset_tree("tif"), *chans, *bs),
        t2zw_case("t2zw.all_tifs_default_chunks", dataset_tree("tif")),
        t2zw_case(
            "t2zw.input_bbox",
            dataset_tree("tif"),
            "channelPatterns",
            ["CamA_ch0"],
            "inputBbox",
            BOX,
            "resultDirName",
            "tiff_to_zarr",
        ),
        t2zw_case(
            "t2zw.regions",
            dataset_tree("tif"),
            "channelPatterns",
            ["CamA_ch0"],
            "readWholeTiff",
            False,
            "inputBbox",
            [1, 1, 1, 40, 36, 12],
            *bs,
        ),
        t2zw_case(
            "t2zw.flipped_tiles_tile_out_bbox",
            dataset_tree("tif"),
            *chans,
            *bs,
            "flippedTile",
            [True, False],
            "tileOutBbox",
            OUT_BOX,
        ),
        t2zw_case(
            "t2zw.partial_files",
            dataset_tree("tif", extra=("scan_CamA_ch0_t0000_part001",)),
            "channelPatterns",
            ["CamA_ch0"],
            "partialFile",
            True,
            *bs,
        ),
        t2zw_case(
            "t2zw.boxes_per_location",
            dataset_tree("tif"),
            *chans,
            *bs,
            "locIds",
            [1, 2],
            "inputBbox",
            np.array([BOX, [1, 1, 1, 40, 36, 12]], dtype=np.float64),
            "tileOutBbox",
            [2, 2, 2, 20, 20, 6],
        ),
        t2zw_case("t2zw.resample", dataset_tree("tif"), *chans, *bs, "resampleFactor", [1, 2]),
    ]

    def paths_code(sub):
        return (
            "d = in{1}; XR_tiffToZarr_wrapper('', 'tiffFullpaths', "
            + sub
            + ", in{2 : end}); out = {};"
        )

    cases += [
        Case(
            "t2zw.tiff_fullpaths",
            code=paths_code("{[d, 'scan_CamA_ch0_t0000.tif'], [d, 'scan_CamB_ch0_t0000.tif']}"),
            nargout=0,
            args=[InPath("d"), "channelPatterns", ["CamA_ch0", "CamB_ch0"], *bs, *LOCAL],
            python=lambda d, *a: py_t2z_wrapper(
                "",
                "tiffFullpaths",
                [d + "scan_CamA_ch0_t0000.tif", d + "scan_CamB_ch0_t0000.tif"],
                *a,
            ),
            inputs={"d": Tree(dataset_tree("tif"), "d")},
        ),
        Case(
            "t2zw.tiff_folder",
            code=paths_code("d(1 : end - 1)"),
            nargout=0,
            args=[InPath("d"), "channelPatterns", ["Cam"], *bs, *LOCAL],
            python=lambda d, *a: py_t2z_wrapper("", "tiffFullpaths", d[:-1], *a),
            inputs={"d": Tree(dataset_tree("tif"), "d")},
        ),
    ]
    return cases


def _z2t_wrapper_cases():
    def case(cid, *args):
        return Case(
            cid,
            code="XR_zarrToTiff_wrapper(in{:}); out = {};",
            nargout=0,
            args=[InPath("d"), *args, *LOCAL],
            python=py_z2t_wrapper,
            inputs={"d": Tree(dataset_tree("zarr"), "d")},
        )

    return [
        case("z2tw.default", "channelPatterns", ["CamA_ch0", "CamB_ch0"]),
        case("z2tw.expression", "channelPatterns", ["Cam[AB]_ch0"], "resultDirName", "tifs"),
    ]


CASES = _t2z_cases() + _z2t_cases() + _t2z_wrapper_cases() + _z2t_wrapper_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# TiffZarrConversionTest in intent, on (z, y, x) data and native zarr arrays
# ---------------------------------------------------------------------------------------------

IM_SIZE = (12, 40, 36)  # MATLAB's ImSize [40 36 12]
CHANNELS = ("CamA_ch0", "CamB_ch0")


def write_dataset(folder, fmt):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ims = []
    for k, n in enumerate(DATASET):
        a = volume(IM_SIZE, "uint16", f"py/{fmt}/{k}")
        ims.append(a)
        if fmt == "zarr":
            write_zarr(a, str(folder / f"{n}.zarr"), block_size=(4, 16, 16))
        else:
            write_tiff(a, str(folder / f"{n}.tif"))
    return str(folder) + "/", ims, list(DATASET)


def assert_same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, a.shape, b.dtype, b.shape)
    assert a.tobytes() == b.tobytes()


def test_tiff_to_zarr_wrapper(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "tif")
    convert.tiff_to_zarr_dataset(
        path, channel_patterns=CHANNELS, block_size=(4, 16, 16), parse_cluster=False
    )
    for i in range(2):
        fn = f"{path}zarr/{names[i]}.zarr"
        assert_same(read_zarr(fn), ims[i])
        info = get_zarr_info(fn)
        assert info.block_size == (4, 16, 16)
        assert info.compressor["cname"] == "zstd"
        assert info.order == "C"
    assert not os.path.isdir(f"{path}zarr/{names[2]}.zarr"), "a file of no channel was converted"
    for f in ("parameters.mat", "parameters.json"):
        assert os.path.isfile(f"{path}zarr/{f}")


def test_tiff_to_zarr_wrapper_with_input_bbox(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "tif")
    b = (2, 4, 6, 10, 30, 33)
    convert.tiff_to_zarr_dataset(
        path,
        channel_patterns=CHANNELS[:1],
        input_bbox=b,
        result_dir_name="tiff_to_zarr",
        parse_cluster=False,
    )
    out = read_zarr(f"{path}tiff_to_zarr/{names[0]}.zarr")
    assert_same(out, ims[0][2:10, 4:30, 6:33])


def test_tiff_to_zarr_wrapper_without_reading_whole_tiff(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "tif")
    convert.tiff_to_zarr_dataset(
        path,
        channel_patterns=CHANNELS[:1],
        read_whole_tiff=False,
        input_bbox=(0, 0, 0, *IM_SIZE),
        block_size=(4, 16, 16),
        parse_cluster=False,
    )
    assert_same(read_zarr(f"{path}zarr/{names[0]}.zarr"), ims[0])


def test_zarr_to_tiff_wrapper(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "zarr")
    convert.zarr_to_tiff_dataset(path, channel_patterns=CHANNELS, parse_cluster=False)
    for i in range(2):
        assert_same(read_tiff(f"{path}tiffs/{names[i]}.tif"), ims[i])
    assert not os.path.exists(f"{path}tiffs/{names[2]}.tif"), "a file of no channel was converted"


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls (MATLAB's FrameworkSpyFixture); nothing runs"""
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        return [True] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    return calls


def test_zarr_to_tiff_wrapper_passes_job_settings(tmp_path, framework_spy):
    # each task converts one file: the memory per task is for one file
    path, _, _ = write_dataset(tmp_path / "d", "zarr")
    convert.zarr_to_tiff_dataset(
        path, channel_patterns=CHANNELS, parse_cluster=False, master_compute=False, cpus_per_task=3
    )
    c = framework_spy[0]
    assert c["parse_cluster"] is False
    assert c["master_compute"] is False
    assert c["cpus_per_task"] == 3
    assert c["max_trial_num"] == 2
    assert c["mem_allocate"] == pytest.approx(np.prod(IM_SIZE) * 4 / 1024**3 * 2.5, rel=1e-12)


def test_tiff_to_zarr_wrapper_job_settings(tmp_path, framework_spy):
    # memory of the last file's group times 4 (MATLAB: a cell of functions is never empty), or
    # 2.75 for the flat-field function; twice that without big_data
    path, _, _ = write_dataset(tmp_path / "d", "tif")
    fn = tmp_path / "ff.txt"
    fn.write_text("petakit5d.flatfield:flatfield_correct_file\n", newline="\n")
    gb = np.prod(IM_SIZE) * 2 / 1024**3
    convert.tiff_to_zarr_dataset(path, channel_patterns=CHANNELS, parse_cluster=False)
    convert.tiff_to_zarr_dataset(
        path, channel_patterns=CHANNELS, parse_cluster=False, big_data=False, cpus_per_task=2
    )
    convert.tiff_to_zarr_dataset(
        path, channel_patterns=CHANNELS, parse_cluster=False, process_fun_path=(str(fn), str(fn))
    )
    assert framework_spy[0]["mem_allocate"] == pytest.approx(gb * 4, rel=1e-12)
    assert framework_spy[1]["mem_allocate"] == pytest.approx(gb * 8, rel=1e-12)
    assert framework_spy[1]["cpus_per_task"] == 2
    assert framework_spy[2]["mem_allocate"] == pytest.approx(gb * 2.75, rel=1e-12)
    assert framework_spy[0]["max_trial_num"] == 2
    t = framework_spy[2]["tasks"][0]
    assert t.func == "petakit5d.io.convert:tiff_to_zarr"
    assert t.kwargs["usr_fcn"] == "petakit5d.flatfield:flatfield_correct_file"
    assert framework_spy[0]["outputs"][0].endswith("/zarr/scan_CamA_ch0_t0000.zarr/")


def test_tiff_to_zarr_wrapper_tries_again_with_more(tmp_path, monkeypatch):
    from petakit5d.jobs import framework

    calls = []
    monkeypatch.setattr(
        framework, "run_tasks", lambda i, o, t, **kw: calls.append(kw) or [False] * len(o)
    )
    path, _, _ = write_dataset(tmp_path / "d", "tif")
    convert.tiff_to_zarr_dataset(path, channel_patterns=CHANNELS, parse_cluster=False)
    assert len(calls) == 2
    assert calls[1]["cpus_per_task"] == 2
    assert calls[1]["mem_allocate"] == pytest.approx(2 * calls[0]["mem_allocate"])


def test_wrapper_user_functions_follow_channels_and_keywords(tmp_path, framework_spy):
    # rows of process_fun_path are channels, columns alternatives: the first whose keywords one
    # is in the file name (.mat files with usrFun and keywords; .txt files: the first line, with
    # the keyword '', which every name contains), else the first
    from petakit5d.io import write_mat

    path, _, _ = write_dataset(tmp_path / "d", "tif")
    a = tmp_path / "a.txt"
    a.write_text("pkg_a:f\n", newline="\n")
    b = tmp_path / "b.mat"
    write_mat(str(b), {"usrFun": "pkg_b:g", "keywords": np.array(["CamB"], dtype=object)})
    c = tmp_path / "c.mat"
    write_mat(str(c), {"usrFun": "pkg_c:h", "keywords": np.array(["nowhere"], dtype=object)})
    convert.tiff_to_zarr_dataset(
        path,
        channel_patterns=("CamA", "CamB"),
        parse_cluster=False,
        process_fun_path=((str(a), ""), (str(c), str(b))),
    )
    got = [t.kwargs["usr_fcn"] for t in framework_spy[0]["tasks"]]
    # files: scan_CamA_ch0, scan_CamA_ch1 (dir order), scan_CamB_ch0
    names = [os.path.basename(t.args[0][0]) for t in framework_spy[0]["tasks"]]
    want = {"scan_CamA_ch0_t0000.tif": "pkg_a:f", "scan_CamA_ch1_t0000.tif": "pkg_a:f"}
    want["scan_CamB_ch0_t0000.tif"] = "pkg_b:g"
    assert dict(zip(names, got, strict=True)) == want


def test_wrapper_boxes_per_location_share_a_single_box(tmp_path, framework_spy):
    # one box per tile location for one setting, one box (or none) for the other: every
    # location gets that one box
    path, _, _ = write_dataset(tmp_path / "d", "tif")
    rows = ((0, 0, 0, 12, 40, 36), (1, 1, 1, 10, 30, 30))
    for k, (inb, outb) in enumerate(((rows[1], rows), (rows, rows[1]), (None, rows))):
        convert.tiff_to_zarr_dataset(
            path,
            channel_patterns=CHANNELS,
            parse_cluster=False,
            loc_ids=(1, 2),
            input_bbox=inb,
            tile_out_bbox=outb,
            result_dir_name=f"z{k}",
        )
        kw = [t.kwargs for t in framework_spy[k]["tasks"]]
        one = "input_bbox" if k != 1 else "tile_out_bbox"
        per = "tile_out_bbox" if k != 1 else "input_bbox"
        want = None if k == 2 else list(rows[1])
        assert [w[one] for w in kw] == [want, want], k
        assert [w[per] for w in kw] == [list(r) for r in rows], k


def test_wrapper_passes_resample_factors_as_matlab_text(tmp_path, framework_spy):
    # mat2str keeps 15 significant digits
    path, _, _ = write_dataset(tmp_path / "d", "tif")
    convert.tiff_to_zarr_dataset(
        path, channel_patterns=CHANNELS, parse_cluster=False, resample_factor=(1 / 3, 1.0, 1.0)
    )
    rs = framework_spy[0]["tasks"][0].kwargs["resample_factor"]
    assert rs == [0.333333333333333, 1.0, 1.0]


def test_single_file_round_trip(tmp_path):
    im = volume(IM_SIZE, "uint16", "round_trip")
    write_tiff(im, str(tmp_path / "im.tif"))
    convert.tiff_to_zarr(
        str(tmp_path / "im.tif"), str(tmp_path / "im.zarr"), block_size=(4, 16, 16)
    )
    assert_same(read_zarr(str(tmp_path / "im.zarr")), im)
    assert get_zarr_info(str(tmp_path / "im.zarr")).block_size == (4, 16, 16)
    convert.zarr_to_tiff(str(tmp_path / "im.zarr"), str(tmp_path / "back.tif"))
    assert_same(read_tiff(str(tmp_path / "back.tif")), im)


def test_tiff_to_zarr_from_memory(tmp_path):
    im = volume(IM_SIZE, "single", "memory", True)
    convert.tiff_to_zarr("", str(tmp_path / "im.zarr"), im, block_size=(4, 16, 16))
    assert_same(read_zarr(str(tmp_path / "im.zarr")), im)


def test_tiff_to_zarr_processing_options(tmp_path):
    # flip in z, then the crop (tile_out_bbox) and the user function, in that order; the result
    # has the input's class (MATLAB's rounding)
    im = volume(IM_SIZE, "uint16", "options")
    write_tiff(im, str(tmp_path / "im.tif"))
    b = (1, 2, 3, 11, 38, 30)
    convert.tiff_to_zarr(
        str(tmp_path / "im.tif"),
        str(tmp_path / "flip.zarr"),
        block_size=(4, 16, 16),
        flip_zstack=True,
        tile_out_bbox=b,
        usr_fcn=lambda x: (x.astype(np.float64) / 3.0).astype(np.float32),
    )
    want = im[::-1][1:11, 2:38, 3:30].astype(np.float64) / 3.0
    assert_same(read_zarr(str(tmp_path / "flip.zarr")), cast(want.astype(np.float32), "uint16"))


def test_tiff_to_zarr_resampling_needs_t04(tmp_path):
    im = volume(IM_SIZE, "uint16", "rs")
    write_tiff(im, str(tmp_path / "im.tif"))
    try:
        from petakit5d.tools.resample import imresize3_fast  # noqa: F401
    except ImportError:
        with pytest.raises(NotImplementedError, match="T04") as err:
            convert.tiff_to_zarr(
                str(tmp_path / "im.tif"), str(tmp_path / "rs.zarr"), resample_factor=(2, 1)
            )
        assert err.value.identifier == "petakit5d:convert:resampleNotPorted"
        assert not list(tmp_path.glob("rs.zarr*")), "a result was left"
        return
    convert.tiff_to_zarr(
        str(tmp_path / "im.tif"), str(tmp_path / "rs.zarr"), resample_factor=(2, 1)
    )
    assert read_zarr(str(tmp_path / "rs.zarr")).shape == (6, 40, 36)


def test_resample_factors_complete_as_matlab():
    rf = convert._resample_factors
    assert rf(None) is None
    assert rf((1, 1, 1)) is None
    assert rf(2) == [2.0, 2.0, 2.0]
    assert rf((2, 0.5)) == [2.0, 0.5, 0.5]  # (z, xy): MATLAB's [xy z] -> [xy xy z]
    assert rf((2, 3, 4)) == [2.0, 3.0, 4.0]


def test_resampling_calls_the_t04_interface(tmp_path, monkeypatch):
    # out size round(size ./ rs) in (z, y, x), before the user function
    import sys
    import types

    seen = []

    def imresize3_fast(vol, out_size):
        seen.append((vol.shape, out_size))
        return np.zeros(out_size, vol.dtype)

    mod = types.ModuleType("petakit5d.tools.resample")
    mod.imresize3_fast = imresize3_fast
    monkeypatch.setitem(sys.modules, "petakit5d.tools.resample", mod)
    im = volume(IM_SIZE, "uint16", "rs2")
    convert.tiff_to_zarr(
        "", str(tmp_path / "rs.zarr"), im, resample_factor=(2.5, 3), usr_fcn=lambda x: x + 1
    )
    assert seen == [((12, 40, 36), (5, 13, 12))]
    assert_same(read_zarr(str(tmp_path / "rs.zarr")), np.ones((5, 13, 12), np.uint16))


def test_tiff_to_zarr_reads_regions(tmp_path):
    im = volume(IM_SIZE, "uint16", "regions")
    write_tiff(im, str(tmp_path / "im.tif"))
    src = str(tmp_path / "im.tif")
    convert.tiff_to_zarr(
        src, str(tmp_path / "whole.zarr"), block_size=(4, 16, 16), read_whole_tiff=False
    )
    assert_same(read_zarr(str(tmp_path / "whole.zarr")), im)
    b = (1, 2, 3, 11, 38, 30)
    convert.tiff_to_zarr(
        src,
        str(tmp_path / "crop.zarr"),
        block_size=(4, 16, 16),
        read_whole_tiff=False,
        tile_out_bbox=b,
    )
    assert_same(read_zarr(str(tmp_path / "crop.zarr")), im[1:11, 2:38, 3:30])


def test_regions_are_read_by_chunk_layers(tmp_path, monkeypatch):
    # one read of the block's slices per chunk layer, as blockedImage read them
    from petakit5d import io

    im = volume(IM_SIZE, "uint16", "layers")
    write_tiff(im, str(tmp_path / "im.tif"))
    ranges = []
    real = io.read_tiff

    def spy(path, rng=None, **kw):
        ranges.append(rng)
        return real(path, rng, **kw)

    monkeypatch.setattr(io, "read_tiff", spy)
    convert.tiff_to_zarr(
        str(tmp_path / "im.tif"),
        str(tmp_path / "r.zarr"),
        block_size=(5, 16, 16),
        read_whole_tiff=False,
        input_bbox=(1, 0, 0, 12, 40, 36),
    )
    assert ranges == [(1, 6), (6, 11), (11, 12)]
    assert_same(read_zarr(str(tmp_path / "r.zarr")), im[1:])


def test_tiff_to_zarr_partial_files_and_frames(tmp_path):
    im = volume(IM_SIZE, "uint16", "parts")
    write_tiff(im[:5], str(tmp_path / "part1.tif"))
    write_tiff(im[5:], str(tmp_path / "part2.tif"))
    convert.tiff_to_zarr(
        [str(tmp_path / "part1.tif"), str(tmp_path / "part2.tif")],
        str(tmp_path / "parts.zarr"),
        block_size=(4, 16, 16),
    )
    assert_same(read_zarr(str(tmp_path / "parts.zarr")), im)
    convert.tiff_to_zarr(
        "", str(tmp_path / "frame.zarr"), im, block_size=(4, 16, 16), flip_zstack=True
    )
    assert_same(read_zarr(str(tmp_path / "frame.zarr")), im[::-1])
    write_tiff(im[0], str(tmp_path / "im2d.tif"))
    convert.tiff_to_zarr(
        str(tmp_path / "im2d.tif"), str(tmp_path / "im2d.zarr"), block_size=(4, 16, 16)
    )
    assert_same(read_zarr(str(tmp_path / "im2d.zarr")), im[0])
    assert get_zarr_info(str(tmp_path / "im2d.zarr")).size == (1, 40, 36)


def test_tiff_to_zarr_writes_private_result_in_place(tmp_path, monkeypatch):
    # with a full uuid made by the call no other process writes the temporary result, so its
    # chunks may be written in place (uuid false); a given uuid can be shared by copies of a
    # task, which keep cpp-zarr's temporary chunk files
    from petakit5d import io

    write_tiff(volume(IM_SIZE, "uint16", "private"), str(tmp_path / "im.tif"))
    calls = []
    real = io.write_zarr

    def spy(data, path, **kw):
        calls.append((path, kw.get("uuid", True)))
        return real(data, path, **kw)

    monkeypatch.setattr(io, "write_zarr", spy)
    src = str(tmp_path / "im.tif")
    given = "01234567-89ab-cdef-0123-456789abcdef"
    convert.tiff_to_zarr(src, str(tmp_path / "whole.zarr"), block_size=(4, 16, 16))
    convert.tiff_to_zarr(
        src, str(tmp_path / "regions.zarr"), block_size=(4, 16, 16), read_whole_tiff=False
    )
    convert.tiff_to_zarr(src, str(tmp_path / "given.zarr"), block_size=(4, 16, 16), uuid=given)
    assert len(calls) == 5  # the whole image, three z blocks of 4, the given uuid
    from petakit5d.core.ids import get_uuid

    private = len(get_uuid()) >= 32
    for path, flag in calls[:4]:
        assert flag is (not private), path
    assert calls[4] == (str(tmp_path / "given.zarr") + "_" + given, True)


def test_results_are_moved_into_place(tmp_path):
    # a task is done once its result exists: it is written under <name>_<uuid> and renamed
    from petakit5d.testing.publish_spy import PublishSpy

    im = volume(IM_SIZE, "uint16", "moved")
    write_tiff(im, str(tmp_path / "im.tif"))
    with PublishSpy() as spy:
        convert.tiff_to_zarr(str(tmp_path / "im.tif"), str(tmp_path / "im.zarr"))
        convert.zarr_to_tiff(str(tmp_path / "im.zarr"), str(tmp_path / "back.tif"))
    assert spy.was_moved(tmp_path / "im.zarr")
    assert spy.was_moved(tmp_path / "back.tif")
    assert sorted(os.listdir(tmp_path)) == ["back.tif", "im.tif", "im.zarr"]


def test_a_result_another_writer_published_is_kept(tmp_path, monkeypatch):
    # also with overwrite: a published result is never deleted (MATLAB replaced it)
    from petakit5d import io

    im = volume(IM_SIZE, "uint16", "other")
    write_tiff(im, str(tmp_path / "im.tif"))
    other = np.full((2, 3, 4), 7, np.uint16)
    final = str(tmp_path / "im.zarr")
    real = io.write_zarr

    def racing(data, path, **kw):
        real(data, path, **kw)
        if not os.path.isdir(final):
            real(other, final)

    monkeypatch.setattr(io, "write_zarr", racing)
    convert.tiff_to_zarr(str(tmp_path / "im.tif"), final, overwrite=True)
    assert_same(read_zarr(final), other)
    assert sorted(os.listdir(tmp_path)) == ["im.tif", "im.zarr"], "the temporary result was left"


def test_existing_results_are_kept_or_replaced(tmp_path, capsys):
    im = volume(IM_SIZE, "uint16", "exists")
    write_tiff(im, str(tmp_path / "im.tif"))
    out = str(tmp_path / "out.zarr")
    write_zarr(np.zeros((2, 2, 2), np.uint16), out)
    convert.tiff_to_zarr(str(tmp_path / "im.tif"), out + "/")
    assert "already exists, skip it!" in capsys.readouterr().out
    assert read_zarr(out).shape == (2, 2, 2)
    convert.tiff_to_zarr(str(tmp_path / "im.tif"), out, overwrite=True)
    assert_same(read_zarr(out), im)
    tif = str(tmp_path / "back.tif")
    write_tiff(np.zeros((2, 2), np.uint16), tif)
    convert.zarr_to_tiff(out, tif)
    assert read_tiff(tif).shape == (2, 2)
    convert.zarr_to_tiff(out, tif, overwrite=True)
    assert_same(read_tiff(tif), im)


def test_matlab_layout_zarr_converts_to_the_same_tiff(tmp_path):
    im = volume(IM_SIZE, "single", "layout", True)
    synth.write_zarr(tmp_path / "m.zarr", im, (4, 16, 16), "matlab")
    write_zarr(im, str(tmp_path / "n.zarr"), block_size=(4, 16, 16))
    convert.zarr_to_tiff(str(tmp_path / "m.zarr"), str(tmp_path / "m.tif"))
    convert.zarr_to_tiff(str(tmp_path / "n.zarr"), str(tmp_path / "n.tif"))
    assert (tmp_path / "m.tif").read_bytes() == (tmp_path / "n.tif").read_bytes()


def test_image_not_expanded_gets_matlabs_chunks(tmp_path):
    # an image without expand_2d_dim is stored [ny, nx] with chunks keeping prod(block_size):
    # along the shorter side whole, the other round(prod / that side)
    im = volume((1, 23, 61), "uint16", "chunks2d")[0]
    convert.tiff_to_zarr("", str(tmp_path / "a.zarr"), im, expand_2d_dim=False)
    info = get_zarr_info(str(tmp_path / "a.zarr"))
    assert info.size == (23, 61)
    assert info.block_size == (23, 61)
    convert.tiff_to_zarr(
        "", str(tmp_path / "b.zarr"), im, expand_2d_dim=False, block_size=(1, 8, 8)
    )
    # [8 8] keeps 64 voxels: ny = 23 whole, round(64 / 23) = 3 along x
    assert get_zarr_info(str(tmp_path / "b.zarr")).block_size == (23, 3)
    assert_same(read_zarr(str(tmp_path / "b.zarr")), im)


def test_crop_boxes_outside_the_image_fail_as_matlab(tmp_path):
    im = volume(IM_SIZE, "uint16", "outside")
    with pytest.raises(convert.ConvertError) as err:
        convert.tiff_to_zarr("", str(tmp_path / "a.zarr"), im, tile_out_bbox=(0, 0, 0, 13, 40, 36))
    assert err.value.identifier == "crop:inputError"
    write_tiff(im, str(tmp_path / "im.tif"))
    with pytest.raises(convert.ConvertError):
        convert.tiff_to_zarr(
            str(tmp_path / "im.tif"),
            str(tmp_path / "b.zarr"),
            read_whole_tiff=False,
            tile_out_bbox=(0, 0, 0, 12, 41, 36),
        )


# ---------------------------------------------------------------------------------------------
# user functions
# ---------------------------------------------------------------------------------------------


def add_seven(vol, extra=0):
    return vol + np.uint16(7 + extra)


def test_user_function_forms():
    uf = convert.user_function
    assert uf("") is None and uf(None) is None
    g = uf("numpy:negative")
    assert g is np.negative
    rec = uf('{"func": "numpy:multiply", "args": [3]}')
    assert rec(np.array([2])).tolist() == [6]
    rec2 = uf({"func": "numpy:clip", "kwargs": {"a_min": 0, "a_max": 1}})
    assert rec2(np.array([-2, 5])).tolist() == [0, 1]
    assert uf(np.sqrt) is np.sqrt
    with pytest.raises(ValueError, match="package.module:function"):
        uf("not_an_import_path")


def test_matlab_code_needs_the_tile_function_parser(monkeypatch):
    import sys

    monkeypatch.setitem(sys.modules, "petakit5d.stitch.tile_functions", None)
    with pytest.raises(NotImplementedError) as err:
        convert.user_function("@(x) x + 1")
    assert err.value.identifier == "petakit5d:convert:matlabCodeFunction"


def test_user_function_by_import_path_through_a_task(tmp_path):
    # the dataset wrapper's tasks carry the function as text (import path or JSON record)
    path, ims, names = write_dataset(tmp_path / "d", "tif")
    fn = tmp_path / "fn.txt"
    fn.write_text('{"func": "numpy:add", "args": [5]}\n', newline="\n")
    convert.tiff_to_zarr_dataset(
        path, channel_patterns=CHANNELS[:1], parse_cluster=False, process_fun_path=str(fn)
    )
    assert_same(read_zarr(f"{path}zarr/{names[0]}.zarr"), ims[0] + np.uint16(5))


def test_matlab_layer_takes_a_python_callable(tmp_path):
    from petakit5d import matlab

    im = volume(IM_SIZE, "uint16", "callable")
    matlab.tiffToZarr("", str(tmp_path / "a.zarr"), to_matlab(im), "usrFcn", add_seven)
    assert_same(read_zarr(str(tmp_path / "a.zarr")), im + np.uint16(7))


def test_zarr_to_tiff_wrapper_user_function_by_import_path(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "zarr")
    convert.zarr_to_tiff_dataset(
        path, channel_patterns=CHANNELS[:1], parse_cluster=False, usr_fcn="numpy:negative"
    )
    assert_same(read_tiff(f"{path}tiffs/{names[0]}.tif"), np.negative(ims[0]))


def test_wrapper_parameters_files_hold_matlabs_values(tmp_path):
    path, _, _ = write_dataset(tmp_path / "d", "tif")
    convert.tiff_to_zarr_dataset(
        path,
        channel_patterns=CHANNELS,
        parse_cluster=False,
        block_size=(4, 16, 8),
        input_bbox=(1, 2, 3, 10, 30, 33),
    )
    pr = json.loads((tmp_path / "d" / "zarr" / "parameters.json").read_text())
    assert pr["blockSize"] == [16, 8, 4]
    assert pr["inputBbox"] == [3, 4, 2, 30, 33, 10]
