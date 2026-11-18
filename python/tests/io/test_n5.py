"""N5 and zarr conversion and the BigStitcher XML (petakit5d.io.n5, task B09).

MATLAB has no tests of these functions; the oracle cases are new:

- zarrToN5 of MATLAB-layout zarr arrays (uint16 and single, chunks and batches that don't
  divide the image, a box, all-zero chunks, defaults, uint8's error): attributes.json's text,
  the chunk files, their headers and the data read back by MATLAB's N5Adapter. MATLAB's side
  needs zarr-python 2 and numcodecs in the active environment (``pkpy_n5_setup``); the Python
  side needs numcodecs (the n5 extra) for blosc, and these cases skip without it.
- N5ToZarr of N5 datasets in MATLAB's layout (raw chunks, so the inputs don't depend on a codec
  library's version): uint16 and single, flipEmptyValue, a box, edge chunks cut as n5-java
  writes them, a folder with a .zgroup (its s0), uint8's error; the zarr results compared by
  data and metadata after the layout mapping.
- constructXmlFile with writestruct, removeTextTags and cleanXmlFile (the demo's steps): tiff
  tiles, MATLAB-layout zarr tiles (the demo), native zarr tiles with dataOrder 'zyx', mapped
  file names without cameras, a channel filtered out (MATLAB's warning), the defaults, and the
  errors (no tiles, an xcorr file with a dataOrder MATLAB doesn't map): the XML text equal.

Then the N5 layer itself (round trips of every codec, region writes, cut and missing chunks, a
cross-check with zarr-python 2 where it is installed), the converters' round trips, tasks and
flags, and the XML pieces.
"""

from __future__ import annotations

import json
import os
import struct
from pathlib import Path

import numpy as np
import pytest

from petakit5d.io import n5
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, InPath, PerAxis, Tree
from petakit5d.testing.synth import CounterRNG

try:
    import numcodecs  # noqa: F401

    HAS_NUMCODECS = True
except ImportError:
    HAS_NUMCODECS = False

needs_numcodecs = pytest.mark.skipif(not HAS_NUMCODECS, reason="numcodecs (the n5 extra) missing")
HELPERS = ("tests/io/matlab_n5",)
RAW = {"type": "raw"}


def volume(shape, cls, name):
    rng = CounterRNG(0, f"b09/n5/{name}")
    shape = tuple(int(s) for s in shape)
    if cls == "uint8":
        return rng.integers(0, 255, shape, dtype=np.uint8)
    if cls == "uint16":
        return rng.integers(0, 60000, shape, dtype=np.uint16)
    return np.ascontiguousarray((rng.normal(shape) * 100.0).astype(np.float32))


def write_raw_n5(folder, vol, block, cut_edges=False, axes="yxz"):
    """an N5 dataset of ``vol`` with raw chunks, in MATLAB's layout (``axes`` 'yxz', the
    oracle inputs) or standard N5 ('zyx'); ``cut_edges``: edge chunks stored only as far as
    the data goes (as n5-java writes them)"""
    arr = n5.create_n5(folder, vol.shape, block, vol.dtype, compression=RAW, axes=axes)
    n5.write_n5(vol, arr, num_threads=1)
    if cut_edges:
        vshape, vblock = arr.view_shape, arr.view_block
        for root, _, files in os.walk(folder):
            for f in files:
                if f == "attributes.json":
                    continue
                p = Path(root) / f
                rel = p.relative_to(folder).parts
                vidx = [int(i) for i in reversed(rel)]
                raw = p.read_bytes()
                full = np.frombuffer(raw[16:], dtype=arr.dtype.newbyteorder(">")).reshape(vblock)
                ext = [min(vblock[d], vshape[d] - vidx[d] * vblock[d]) for d in range(3)]
                if ext == list(vblock):
                    continue
                part = np.ascontiguousarray(full[: ext[0], : ext[1], : ext[2]])
                p.write_bytes(struct.pack(">HH3I", 0, 3, *reversed(ext)) + part.tobytes())


def n5_summary(path, axes="zyx"):
    """pkpy_n5_summary in Python: attributes.json's text, the chunk files, their headers, the
    data (z, y, x). A standard N5 (the port's) is given in MATLAB's layout: dimensions, block
    sizes, chunk names and header dimensions permuted from [x, y, z] to [z, x, y], so that it is
    compared with MATLAB's file after the permutation"""
    path = os.fspath(path)
    with open(os.path.join(path, "attributes.json"), encoding="utf-8", newline="") as f:
        txt = f.read()
    rel = []
    for root, _, files in os.walk(path):
        for f in files:
            r = os.path.relpath(os.path.join(root, f), path).replace(os.sep, "/")
            if r != "attributes.json":
                rel.append(r)
    rel.sort()
    heads = np.zeros((len(rel), 16), dtype=np.uint8)
    for i, r in enumerate(rel):
        with open(os.path.join(path, r), "rb") as f:
            h = f.read(16)
        heads[i, : len(h)] = np.frombuffer(h, dtype=np.uint8)
    data = n5.read_n5(path, axes=axes)
    if axes == "zyx":
        perm = (2, 0, 1)  # [x y z] -> [z x y]
        meta = json.loads(txt)
        for k in ("dimensions", "blockSize"):
            meta[k] = [meta[k][i] for i in perm]
        txt = n5._attributes_text(meta)
        rel = sorted("/".join(r.split("/")[i] for i in perm) for r in rel)
        heads = np.zeros((len(rel), 16), dtype=np.uint8)
        for i, r in enumerate(rel):
            m = r.split("/")
            orig = "/".join(m[j] for j in (1, 2, 0))
            with open(os.path.join(path, orig), "rb") as f:
                h = bytearray(f.read(16))
            d = struct.unpack(">3I", h[4:16])
            h[4:16] = struct.pack(">3I", *(d[j] for j in perm))
            heads[i, : len(h)] = np.frombuffer(bytes(h), dtype=np.uint8)
    return [txt, rel, heads, data]


# ---------------------------------------------------------------------------------------------
# oracle cases: zarrToN5, N5ToZarr
# ---------------------------------------------------------------------------------------------


def zarr_tree(shape, cls, name, chunks=(4, 4, 4), zeros=None):
    def make(folder):
        a = volume(shape, cls, name)
        if zeros is not None:
            a[zeros] = 0
        synth.write_zarr(Path(folder), a, chunks, "matlab")

    return make


def n5_tree(shape, cls, name, block, cut_edges=False, flips=False, zgroup=False):
    def make(folder):
        a = volume(shape, cls, name)
        if flips:
            a.reshape(-1)[::7] = 65535
        target = Path(folder)
        if zgroup:
            target.mkdir(parents=True, exist_ok=True)
            (target / ".zgroup").write_bytes(b'{\n    "zarr_format": 2\n}')
            target = target / "s0"
        write_raw_n5(target, a, block, cut_edges)

    return make


def py_zarr_to_n5(path, *nv):
    from petakit5d import matlab

    z = path.rstrip("/")
    matlab.zarrToN5(z, *nv, "parseCluster", False)
    return n5_summary(os.path.dirname(z) + "/n5/t.n5")


def py_n5_to_zarr(path, *nv):
    from petakit5d import matlab

    # (the inputs are MATLAB's layout, which the file can't tell)
    matlab.N5ToZarr(path.rstrip("/"), *nv, "parseCluster", False, "axes", "yxz")
    return []


Z2N5_CODE = (
    "pkpy_n5_setup(); z = in{1}(1 : end - 1); zarrToN5(z, in{2 : end}, 'parseCluster', false); "
    "out = pkpy_n5_summary([fileparts(z), '/n5/t.n5']);"
)
N5Z_CODE = (
    "pkpy_n5_setup(); N5ToZarr(in{1}(1 : end - 1), in{2 : end}, 'parseCluster', false); out = {};"
)
NEEDS_BLOSC = set()


def z2n5_case(cid, make, *nv, error=False):
    NEEDS_BLOSC.add(cid)
    return Case(
        cid,
        code=Z2N5_CODE,
        nargout=4,
        args=[InPath("z"), *nv],
        python=py_zarr_to_n5,
        inputs={"z": Tree(make, "t.zarr")},
        out_roles=(None, None, None, "vol"),
        files=False,
        helpers=HELPERS,
        expect_error=error,
        py_error=(n5.N5Error,),
    )


def n5z_case(cid, make, *nv, error=False):
    return Case(
        cid,
        code=N5Z_CODE,
        nargout=0,
        args=[InPath("n"), *nv],
        python=py_n5_to_zarr,
        inputs={"n": Tree(make, "in.n5")},
        helpers=HELPERS,
        expect_error=error,
        py_error=(n5.N5Error,),
    )


N5_CASES = [
    z2n5_case(
        "z2n5_uint16_odd",
        zarr_tree((9, 13, 11), "uint16", "z1", (3, 4, 5)),
        "blockSize", [4, 5, 3], "batchSize", [8, 5, 6],
    ),
    z2n5_case(
        "z2n5_single_bbox",
        zarr_tree((7, 12, 10), "single", "z2"),
        "blockSize", [4, 4, 4], "batchSize", [4, 4, 4], "bbox", [2, 3, 2, 11, 9, 6],
    ),
    z2n5_case(
        "z2n5_zero_chunks",
        zarr_tree((8, 10, 9), "uint16", "z3", zeros=np.s_[:4]),
        "blockSize", [5, 5, 4], "batchSize", [5, 5, 4],
    ),
    z2n5_case("z2n5_defaults", zarr_tree((6, 7, 5), "uint16", "z4")),
    z2n5_case(
        "z2n5_uint8_error", zarr_tree((4, 5, 6), "uint8", "z5"), error="MATLAB:UndefinedFunction"
    ),
    n5z_case(
        "n5z_raw_uint16",
        n5_tree((9, 13, 11), "uint16", "n1", (3, 4, 5)),
        "blockSize", [5, 4, 3], "batchSize", [10, 8, 6],
    ),
    n5z_case(
        "n5z_flip_single_bbox",
        n5_tree((7, 12, 10), "single", "n2", (4, 4, 4), flips=True),
        "blockSize", [4, 4, 4], "batchSize", [4, 4, 4], "bbox", [2, 3, 2, 11, 9, 6],
        "flipEmptyValue", True,
    ),
    n5z_case(
        "n5z_flip_uint16",
        n5_tree((6, 9, 8), "uint16", "n3", (4, 4, 4), flips=True),
        "flipEmptyValue", True,
    ),
    n5z_case(
        "n5z_cut_edges",
        n5_tree((9, 13, 11), "uint16", "n4", (4, 5, 3), cut_edges=True),
        "blockSize", [6, 6, 6],
    ),
    n5z_case("n5z_zgroup_s0", n5_tree((5, 6, 7), "uint16", "n5", (4, 4, 4), zgroup=True)),
    n5z_case(
        "n5z_uint8_error",
        n5_tree((4, 5, 6), "uint8", "n6", (4, 4, 4)),
        error="MATLAB:nonExistentField",
    ),
]  # fmt: skip


# ---------------------------------------------------------------------------------------------
# oracle cases: constructXmlFile, writestruct, removeTextTags, cleanXmlFile
# ---------------------------------------------------------------------------------------------

TILE_SHAPE = (4, 6, 5)  # MATLAB's [6 5 4]


def tile_name(cam, k, x, y, ch=0):
    c = f"Cam{cam}_" if cam else ""
    return (
        f"Scan_Iter_0000_0000_{c}ch{ch}_CAM1_stack0000_488nm_0000000msec_"
        f"{1000 + k:010d}msecAbs_{x:03d}x_{y:03d}y_000z_0000t"
    )


# (camera, x index, y index, StageX, StageY, StageZ) of the demo-like 2 x 2 grid and a CamA tile
GRID = [
    ("B", 0, 0, "10.5", "-3", "1.5"),
    ("B", 0, 1, "10.5", "27.125", "3"),
    ("B", 1, 0, "30.75", "-3", "4.5"),
    ("B", 1, 1, "30.75", "27.125", "6"),
    ("A", 0, 0, "10.5", "-3", "1.5"),
]
MAT_HEADER = b"MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Thu Jan  1 00:00:00 2026"
SHIFTS = np.array([[1.5, -2, 3], [0, 0, 0], [-1, 2.25, 0.5], [4, 5, 6]])


def xml_tree(grid=GRID, kind="tiff", mapped=False, xcorr=True, tiles=True):
    def make(folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        head = "Filepath,Filename,StageX_um_,StageY_um_,StageZ_um_"
        rows = [head + (",mappedFilename" if mapped else "")]
        for k, (cam, xi, yi, x, y, z) in enumerate(grid):
            name = tile_name(cam, k, xi, yi)
            fname = f"tile_{k}" if mapped else name
            vol = np.full(TILE_SHAPE, k + 1, dtype=np.uint16)
            if tiles:
                if kind == "tiff":
                    synth.write_tiff(folder / f"{fname}.tif", vol)
                else:
                    layout = "matlab" if kind == "mzarr" else "native"
                    synth.write_zarr(folder / f"{fname}.zarr", vol, TILE_SHAPE, layout)
            row = f"/mnt/projects/run1/{fname}.tif,{fname}.tif,{x},{y},{z}"
            rows.append(row + (f",{name}.tif" if mapped else ""))
        with open(folder / "ImageList.csv", "w", newline="\n") as f:
            f.write("\n".join(rows) + "\n")
        if xcorr:
            from petakit5d.io import write_mat

            # (MAT v6 with a header without time or platform: the same input everywhere)
            write_mat(folder / "xcorr.mat", {"dxyz_shift": SHIFTS}, version="6")
            data = bytearray((folder / "xcorr.mat").read_bytes())
            data[:116] = MAT_HEADER.ljust(116, b" ")
            (folder / "xcorr.mat").write_bytes(bytes(data))

    return make


XML_STEPS = (
    "writestruct(s, 'out.xml', 'StructNodeName', 'SpimData'); removeTextTags('out.xml'); "
    "cleanXmlFile('out.xml'); out = {vs, fileread('out.xml'), t.Filename'};"
)
XML_CODE = (
    "d = in{1}; xc = ''; if in{7}, xc = [d, 'xcorr.mat']; end; "
    "[s, t, vs] = constructXmlFile(d, [d, 'out.n5'], [d, 'ImageList.csv'], in{2 : 6}, xc); "
    + XML_STEPS
)
XML_DEFAULTS_CODE = (
    "d = in{1}; [s, t, vs] = constructXmlFile(d, [d, 'out.n5'], [d, 'ImageList.csv'], in{2}); "
    + XML_STEPS
)


def _xml_steps(s, t, vs):
    n5.write_struct_xml(s, "out.xml", "SpimData")
    n5.remove_text_tags("out.xml")
    n5.clean_xml_file("out.xml")
    with open("out.xml", encoding="utf-8", newline="") as f:
        text = f.read()
    return [np.array(vs), text, list(t["Filename"])]


def py_xml(d, voxel, patterns, zarr_file, axis_order, data_order, xcorr):
    s, t, vs = n5.construct_xml_file(
        d,
        d + "out.n5",
        d + "ImageList.csv",
        voxel,
        patterns,
        zarr_file,
        axis_order,
        data_order,
        d + "xcorr.mat" if xcorr else "",
    )
    return _xml_steps(s, t, vs)


def py_xml_defaults(d, voxel):
    # the defaults but data_order: MATLAB's 'yxz' describes its own N5 layout, the port's
    # default 'zyx' the standard N5 it writes (docs/decisions/0009-n5-layout.md)
    s, t, vs = n5.construct_xml_file(d, d + "out.n5", d + "ImageList.csv", voxel, data_order="yxz")
    return _xml_steps(s, t, vs)


def xml_case(cid, make, *args, error=False, code=XML_CODE, python=py_xml):
    return Case(
        cid,
        code=code,
        nargout=3,
        args=[InPath("d"), *args],
        python=python,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=(n5.N5Error, ValueError),
    )


VOXEL = PerAxis([0.108, 0.2, 0.3], "xyz")
XML_CASES = [
    xml_case("xml_demo_tiff", xml_tree(), VOXEL, ["CamB_ch0"], False, "-xyz", "yxz", True),
    xml_case(
        "xml_demo_zarr", xml_tree(kind="mzarr"), VOXEL, ["CamB_ch0"], True, "-xyz", "yxz", True
    ),
    xml_case("xml_axis_order_tiff", xml_tree(), VOXEL, ["CamB"], False, "y,-x,z", "yxz", True),
    xml_case(
        "xml_native_zarr_zyx",
        xml_tree(kind="nzarr"),
        VOXEL,
        ["CamB_ch0"],
        True,
        "x,-y,z",
        "zyx",
        True,
    ),
    xml_case(
        "xml_mapped_nocam",
        xml_tree([(None, *g[1:]) for g in GRID[:4]], mapped=True, xcorr=False),
        PerAxis([0.1, 0.1, 0.25], "xyz"),
        ["ch0"],
        False,
        "xyz",
        "yxz",
        False,
    ),
    xml_case(
        "xml_defaults",
        xml_tree([("A", *g[1:]) for g in GRID[:4]], xcorr=False),
        VOXEL,
        code=XML_DEFAULTS_CODE,
        python=py_xml_defaults,
    ),
    xml_case(
        "xml_no_tiles_error",
        xml_tree(tiles=False),
        VOXEL,
        ["CamB_ch0"],
        False,
        "xyz",
        "yxz",
        False,
        error=True,
    ),
    xml_case(
        "xml_commas_xcorr_error",
        xml_tree(),
        VOXEL,
        ["CamB_ch0"],
        False,
        "xyz",
        "y,x,z",
        True,
        error="MATLAB:UndefinedFunction",
    ),
]

CASES = N5_CASES + XML_CASES


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.id in NEEDS_BLOSC and not HAS_NUMCODECS:
        pytest.skip("numcodecs (the n5 extra) is missing: the N5 results are blosc-compressed")
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# the N5 layer
# ---------------------------------------------------------------------------------------------

# MATLAB's attributes.json for [7 9 5] uint16 data in [4 4 3] blocks (zarr-python 2's text)
MATLAB_ATTRIBUTES = """{
    "blockSize": [
        3,
        4,
        4
    ],
    "compression": {
        "blocksize": 0,
        "clevel": 1.0,
        "cname": "zstd",
        "shuffle": 1,
        "type": "blosc"
    },
    "dataType": "uint16",
    "dimensions": [
        5,
        9,
        7
    ],
    "n5": "2.0.0"
}"""


def test_attributes_are_zarr_python_2s_in_both_layouts(tmp_path):
    arr = n5.create_n5(tmp_path / "a.n5", (5, 7, 9), (3, 4, 4), "uint16", axes="yxz")
    assert (tmp_path / "a.n5" / "attributes.json").read_text() == MATLAB_ATTRIBUTES
    assert (arr.dims, arr.block, arr.shape, arr.block_size) == (
        (5, 9, 7),
        (3, 4, 4),
        (5, 7, 9),
        (3, 4, 4),
    )
    assert n5.open_n5(tmp_path / "a.n5", "yxz") == arr
    assert arr.matlab_class == "uint16"
    # standard N5: [nx ny nz]
    std = n5.create_n5(tmp_path / "b.n5", (5, 7, 9), (3, 4, 4), "uint16")
    meta = json.loads((tmp_path / "b.n5" / "attributes.json").read_text())
    assert (meta["dimensions"], meta["blockSize"]) == ([9, 7, 5], [4, 4, 3])
    assert (tmp_path / "b.n5" / "attributes.json").read_text() == MATLAB_ATTRIBUTES.replace(
        "3,\n        4,\n        4", "4,\n        4,\n        3"
    ).replace("5,\n        9,\n        7", "9,\n        7,\n        5")
    assert (std.shape, std.block_size, std.axes) == ((5, 7, 9), (3, 4, 4), "zyx")
    assert n5.open_n5(tmp_path / "b.n5") == std
    with pytest.raises(ValueError, match="axes"):
        n5.open_n5(tmp_path / "b.n5", "xyz")


def test_standard_chunks_are_whole_big_endian_x_fastest(tmp_path):
    a = volume((5, 7, 9), "uint16", "layout")
    arr = n5.create_n5(tmp_path / "a.n5", a.shape, (3, 4, 4), a.dtype, compression=RAW)
    n5.write_n5(a, arr)
    names = sorted(
        os.path.relpath(os.path.join(r, f), tmp_path / "a.n5").replace(os.sep, "/")
        for r, _, fs in os.walk(tmp_path / "a.n5")
        for f in fs
        if f != "attributes.json"
    )
    # dimensions [x y z]: 3 x 2 x 2 chunks named ix/iy/iz
    assert names == [f"{i}/{j}/{k}" for i in range(3) for j in range(2) for k in range(2)]
    raw = (tmp_path / "a.n5" / "2" / "1" / "1").read_bytes()
    assert struct.unpack(">HH3I", raw[:16]) == (0, 3, 4, 4, 3)
    blk = np.frombuffer(raw[16:], dtype=">u2").reshape(3, 4, 4)  # (z, y, x) of the chunk
    # the chunk of x 8:12, y 4:8, z 3:6 (padded with zeros past the data)
    ref = np.zeros((3, 4, 4), dtype=np.uint16)
    ref[:2, :3, :1] = a[3:5, 4:7, 8:9]
    np.testing.assert_array_equal(blk, ref)


def test_matlab_layout_chunks_are_whole_big_endian_in_its_order(tmp_path):
    a = volume((5, 7, 9), "uint16", "layout")
    arr = n5.create_n5(tmp_path / "a.n5", a.shape, (3, 4, 4), a.dtype, compression=RAW, axes="yxz")
    n5.write_n5(a, arr)
    names = sorted(
        os.path.relpath(os.path.join(r, f), tmp_path / "a.n5").replace(os.sep, "/")
        for r, _, fs in os.walk(tmp_path / "a.n5")
        for f in fs
        if f != "attributes.json"
    )
    # dimensions [z x y]: 2 x 3 x 2 chunks named iz/ix/iy
    assert names == [f"{i}/{j}/{k}" for i in range(2) for j in range(3) for k in range(2)]
    raw = (tmp_path / "a.n5" / "1" / "2" / "1").read_bytes()
    assert struct.unpack(">HH3I", raw[:16]) == (0, 3, 3, 4, 4)
    blk = np.frombuffer(raw[16:], dtype=">u2").reshape(4, 4, 3)  # (y, x, z) of the chunk
    # the chunk of z 3:6, x 8:12, y 4:8 (padded with zeros past the data)
    ref = np.zeros((4, 4, 3), dtype=np.uint16)
    ref[:3, :1, :2] = a[3:5, 4:7, 8:9].transpose(1, 2, 0)
    np.testing.assert_array_equal(blk, ref)


CODECS = [
    {"type": "raw"},
    {"type": "gzip", "level": -1},
    {"type": "gzip", "level": 3, "useZlib": True},
    {"type": "bzip2", "blockSize": 9},
    {"type": "xz", "preset": 6},
    pytest.param(dict(n5.MATLAB_COMPRESSION), marks=needs_numcodecs),
    pytest.param(
        {"type": "blosc", "cname": "lz4", "clevel": 5, "shuffle": 2}, marks=needs_numcodecs
    ),
    pytest.param(
        {"type": "zstd", "level": 3},
        marks=pytest.mark.skipif(not HAS_NUMCODECS and n5._zstd_module() is None, reason="no zstd"),
    ),
]


@pytest.mark.parametrize("comp", CODECS, ids=lambda c: c["type"] + str(len(c)))
@pytest.mark.parametrize("cls", ["uint16", "single", "uint8"])
def test_round_trips_of_every_codec(tmp_path, comp, cls):
    a = volume((7, 10, 9), cls, f"codec/{cls}")
    arr = n5.create_n5(tmp_path / "a.n5", a.shape, (3, 4, 5), a.dtype, compression=comp)
    n5.write_n5(a, arr, num_threads=3)
    b = n5.read_n5(tmp_path / "a.n5")
    assert b.dtype == a.dtype and b.flags.c_contiguous
    np.testing.assert_array_equal(b.view(np.uint8), a.view(np.uint8))
    box = (1, 2, 3, 6, 9, 8)
    np.testing.assert_array_equal(n5.read_n5(arr, bbox=box), a[1:6, 2:9, 3:8])


def test_region_writes_keep_the_rest_of_their_chunks(tmp_path):
    a = volume((6, 9, 7), "uint16", "regions")
    arr = n5.create_n5(tmp_path / "a.n5", a.shape, (4, 4, 4), a.dtype, compression=RAW)
    n5.write_n5(a, arr)
    b = volume((3, 5, 2), "uint16", "regions/b")
    n5.write_n5(b, arr, start=(2, 3, 4))
    ref = a.copy()
    ref[2:5, 3:8, 4:6] = b
    np.testing.assert_array_equal(n5.read_n5(arr), ref)
    # a region into chunks that have no files yet: the rest of them is 0
    arr2 = n5.create_n5(tmp_path / "b.n5", a.shape, (4, 4, 4), a.dtype, compression=RAW)
    n5.write_n5(b, arr2, start=(1, 1, 1))
    ref = np.zeros_like(a)
    ref[1:4, 1:6, 1:3] = b
    np.testing.assert_array_equal(n5.read_n5(arr2), ref)


def test_cut_and_missing_chunks_read_as_zero_past_their_data(tmp_path):
    a = volume((9, 13, 11), "single", "cut")
    write_raw_n5(tmp_path / "a.n5", a, (4, 5, 3), cut_edges=True)
    raw = (tmp_path / "a.n5" / "2" / "3" / "2").read_bytes()
    assert struct.unpack(">HH3I", raw[:16]) == (0, 3, 1, 2, 3)
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "a.n5", axes="yxz"), a)
    os.remove(tmp_path / "a.n5" / "0" / "0" / "0")
    ref = a.copy()
    ref[:4, :5, :3] = 0
    # a chunk stored smaller than the data it covers (the dataset grew): 0 past its data
    one = np.array([[[123.5]]], dtype=">f4")
    (tmp_path / "a.n5" / "0" / "0" / "1").write_bytes(
        struct.pack(">HH3I", 0, 3, 1, 1, 1) + one.tobytes()
    )
    ref[:4, 5:10, :3] = 0
    ref[0, 5, 0] = 123.5
    out = np.full(a.shape, 7, dtype=a.dtype)
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "a.n5", out=out, axes="yxz"), ref)
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "a.n5", axes="yxz"), ref)
    # standard N5 with cut edges
    write_raw_n5(tmp_path / "s.n5", a, (4, 5, 3), cut_edges=True, axes="zyx")
    raw = (tmp_path / "s.n5" / "3" / "2" / "2").read_bytes()
    assert struct.unpack(">HH3I", raw[:16]) == (0, 3, 2, 3, 1)
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "s.n5"), a)


def test_images_and_casts(tmp_path):
    a = volume((10, 9), "uint16", "image")
    arr = n5.create_n5(tmp_path / "a.n5", a.shape, (4, 4), "uint16", compression=RAW)
    assert arr.dims == (9, 10)
    n5.write_n5(a, arr)
    np.testing.assert_array_equal(n5.read_n5(arr), a)
    arr3 = n5.create_n5(tmp_path / "b.n5", (1, 4, 3), (1, 4, 4), "single", compression=RAW)
    n5.write_n5(np.array([[1.4, -2.5, 70000.0]] * 4), arr3)  # double data: cast to single
    np.testing.assert_array_equal(n5.read_n5(arr3)[0, 0], np.float32([1.4, -2.5, 70000.0]))


def test_errors(tmp_path):
    with pytest.raises(n5.N5Error, match="Unable to open"):
        n5.open_n5(tmp_path / "none.n5")
    (tmp_path / "g").mkdir()
    (tmp_path / "g" / "attributes.json").write_text('{"n5": "2.0.0"}')
    with pytest.raises(n5.N5Error, match="not an N5 dataset"):
        n5.open_n5(tmp_path / "g")
    (tmp_path / "z").mkdir()
    (tmp_path / "z" / ".zgroup").write_text("{}")
    with pytest.raises(n5.N5Error, match="Unable to open"):
        n5.open_n5(tmp_path / "z")
    arr = n5.create_n5(tmp_path / "a.n5", (2, 3, 4), (2, 2, 2), "uint16", compression=RAW)
    with pytest.raises(n5.N5Error, match="out of the bound"):
        n5.read_n5(arr, bbox=(0, 0, 0, 3, 3, 4))
    bad = n5.create_n5(
        tmp_path / "c.n5", (2, 3, 4), (2, 2, 2), "uint16", compression={"type": "lz4"}
    )
    with pytest.raises(n5.N5Error, match="not supported"):
        n5.write_n5(np.zeros((2, 3, 4), np.uint16), bad)


def test_legacy_attributes(tmp_path):
    a = volume((3, 4, 5), "uint16", "legacy")
    write_raw_n5(tmp_path / "a.n5", a, (2, 2, 2), axes="zyx")
    meta = json.loads((tmp_path / "a.n5" / "attributes.json").read_text())
    del meta["compression"]
    meta["compressionType"] = "raw"
    (tmp_path / "a.n5" / "attributes.json").write_text(json.dumps(meta))
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "a.n5"), a)


def test_zarr_python_2_reads_and_writes_the_same_datasets(tmp_path):
    zarr = pytest.importorskip("zarr")
    if not hasattr(zarr, "N5Store"):
        pytest.skip("zarr-python without N5Store (version 3)")
    import warnings

    a = volume((7, 10, 9), "uint16", "zarr2")
    arr = n5.create_n5(tmp_path / "a.n5", a.shape, (3, 4, 5), a.dtype)
    n5.write_n5(a, arr)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        # an independent standard reader: the image itself, (z, y, x)
        z = zarr.open(zarr.N5Store(str(tmp_path / "a.n5")), mode="r")
        assert z.shape == a.shape
        np.testing.assert_array_equal(z[...], a)
        w = zarr.open(
            store=zarr.N5Store(str(tmp_path / "b.n5")),
            mode="w",
            shape=a.shape,
            chunks=(3, 4, 5),
            dtype="u2",
            compressor=numcodecs.Blosc(cname="zstd", clevel=1.0),
            order="C",
        )
        w[...] = a
        # MATLAB's N5Adapter: its (y, x, z) array, read back with axes 'yxz'
        m = zarr.open(
            store=zarr.N5Store(str(tmp_path / "m.n5")),
            mode="w",
            shape=(10, 9, 7),
            chunks=(4, 5, 3),
            dtype="u2",
            compressor=numcodecs.Blosc(cname="zstd", clevel=1.0),
            order="C",
        )
        m[...] = a.transpose(1, 2, 0)
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "m.n5", axes="yxz"), a)
    assert (tmp_path / "b.n5" / "attributes.json").read_bytes() == (
        tmp_path / "a.n5" / "attributes.json"
    ).read_bytes()
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "b.n5"), a)
    for r, _, fs in os.walk(tmp_path / "b.n5"):
        for f in fs:
            p = Path(r) / f
            q = tmp_path / "a.n5" / p.relative_to(tmp_path / "b.n5")
            assert q.is_file(), p
            assert p.read_bytes()[:16] == q.read_bytes()[:16]


# ---------------------------------------------------------------------------------------------
# the converters
# ---------------------------------------------------------------------------------------------


@needs_numcodecs
@pytest.mark.parametrize("layout", ["native", "matlab"])
@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_zarr_n5_zarr_round_trips(tmp_path, layout, cls):
    from petakit5d.io import get_zarr_info, read_zarr

    a = volume((9, 13, 11), cls, f"rt/{cls}")
    synth.write_zarr(tmp_path / "t.zarr", a, (4, 4, 4), layout)
    n5.zarr_to_n5(
        str(tmp_path / "t.zarr"), block_size=(3, 4, 5), batch_size=(6, 8, 5), parse_cluster=False
    )
    arr = n5.open_n5(tmp_path / "n5" / "t.n5")
    assert (arr.shape, arr.block_size, arr.compression) == (
        a.shape,
        (3, 4, 5),
        n5.MATLAB_COMPRESSION,
    )
    assert not (tmp_path / "n5" / "zarr_flag").exists() or not any(
        (tmp_path / "n5" / "zarr_flag").iterdir()
    )
    n5.n5_to_zarr(str(tmp_path / "n5" / "t.n5"), block_size=(5, 5, 5), parse_cluster=False)
    out = str(tmp_path / "n5" / "zarr" / "t.zarr")
    info = get_zarr_info(out)
    assert (info.order, info.size, info.block_size) == ("C", a.shape, (5, 5, 5))
    np.testing.assert_array_equal(read_zarr(out), a)
    # boxes
    n5.zarr_to_n5(
        str(tmp_path / "t.zarr"),
        result_dir_str="box/",
        bbox=(1, 2, 3, 8, 12, 9),
        block_size=(4, 4, 4),
        parse_cluster=False,
    )
    np.testing.assert_array_equal(n5.read_n5(tmp_path / "box" / "t.n5"), a[1:8, 2:12, 3:9])


def test_n5_to_zarr_flips_empty_values(tmp_path):
    from petakit5d.io import read_zarr

    a = volume((5, 6, 7), "uint16", "flip")
    a.reshape(-1)[::5] = 65535
    write_raw_n5(tmp_path / "a.n5", a, (2, 3, 4), axes="zyx")
    n5.n5_to_zarr(str(tmp_path / "a.n5"), flip_empty_value=True, parse_cluster=False)
    write_raw_n5(tmp_path / "m" / "a.n5", a, (2, 3, 4))
    n5.n5_to_zarr(str(tmp_path / "m" / "a.n5"), axes="yxz", parse_cluster=False)
    np.testing.assert_array_equal(read_zarr(str(tmp_path / "m" / "zarr" / "a.zarr")), a)
    ref = a.copy()
    ref[ref == 65535] = 0
    np.testing.assert_array_equal(read_zarr(str(tmp_path / "zarr" / "a.zarr")), ref)


def test_existing_results_are_kept(tmp_path, capsys):
    a = volume((3, 4, 5), "uint16", "kept")
    write_raw_n5(tmp_path / "a.n5", a, (2, 2, 2))
    (tmp_path / "zarr" / "a.zarr").mkdir(parents=True)
    n5.n5_to_zarr(str(tmp_path / "a.n5"), parse_cluster=False)
    assert "The output result exists, skip it!" in capsys.readouterr().out
    assert not any((tmp_path / "zarr" / "a.zarr").iterdir())
    synth.write_zarr(tmp_path / "t.zarr", a, (2, 2, 2), "matlab")
    (tmp_path / "n5" / "t.n5").mkdir(parents=True)
    n5.zarr_to_n5(str(tmp_path / "t.zarr"), parse_cluster=False)
    assert not any((tmp_path / "n5" / "t.n5").iterdir())


def test_matlabs_types_only(tmp_path):
    a = volume((3, 4, 5), "uint8", "types")
    write_raw_n5(tmp_path / "a.n5", a, (2, 2, 2))
    with pytest.raises(n5.N5Error) as e:
        n5.n5_to_zarr(str(tmp_path / "a.n5"), parse_cluster=False)
    assert e.value.identifier == "MATLAB:nonExistentField"
    synth.write_zarr(tmp_path / "t.zarr", a, (2, 2, 2), "matlab")
    with pytest.raises(n5.N5Error) as e:
        n5.zarr_to_n5(str(tmp_path / "t.zarr"), parse_cluster=False)
    assert e.value.identifier == "MATLAB:UndefinedFunction"
    with pytest.raises(n5.N5Error, match="does not exist"):
        n5.zarr_to_n5(str(tmp_path / "none.zarr"), parse_cluster=False)


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls; the tasks don't run, the outputs count as done"""
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        return [True] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    return calls


def test_tasks_flags_and_memory(tmp_path, framework_spy):
    a = volume((10, 10, 10), "uint16", "tasks")
    write_raw_n5(tmp_path / "a.n5", a, (2, 2, 2))
    n5.n5_to_zarr(str(tmp_path / "a.n5"), block_size=(2, 2, 2), batch_size=(2, 2, 2), uuid="u1")
    (call,) = framework_spy
    # 125 batches, 10 per task (N5ToZarr's minimum), flags named by MATLAB's batch numbers
    assert len(call["tasks"]) == 13
    flag = f"{tmp_path}/zarr//zarr_flag/a_u1/"
    assert call["outputs"][0] == flag + "/blocks_1_10.mat"
    assert call["outputs"][-1] == flag + "/blocks_121_125.mat"
    assert call["mem_allocate"] == pytest.approx(8 * 4 / 2**30 * 2.5)
    rec = call["tasks"][1].record()
    assert rec["args"][0] == list(range(10, 20))
    # batch 11 of 5 x 5 x 5 (y fastest, then x, then z): y 0, x 2, z 0
    assert rec["args"][4][0] == [0, 0, 4, 2, 2, 6]
    framework_spy.clear()
    synth.write_zarr(tmp_path / "t.zarr", a, (2, 2, 2), "matlab")
    n5.zarr_to_n5(str(tmp_path / "t.zarr"), block_size=(2, 2, 2), batch_size=(2, 2, 2))
    (call,) = framework_spy
    assert len(call["tasks"]) == 25  # 5 per task


def test_batches_are_numbered_y_fastest(tmp_path):
    _, block, batch, inb, outb = n5._batches((6, 4, 4), (1, 0, 0, 5, 4, 4), (2, 2, 2), (2, 2, 2))
    assert (block, batch) == ([2, 2, 2], [2, 2, 2])
    assert outb[:3].tolist() == [[0, 0, 0, 2, 2, 2], [0, 2, 0, 2, 4, 2], [0, 0, 2, 2, 2, 4]]
    assert (inb - outb)[:, :3].tolist() == [[1, 0, 0]] * len(outb)


def test_block_tasks_skip_done_flags_and_need_their_folder(tmp_path, capsys):
    a = volume((3, 4, 5), "uint16", "blocks")
    write_raw_n5(tmp_path / "a.n5", a, (2, 2, 2))
    flag = tmp_path / "flags" / "blocks_1_1.mat"
    args = ([0], str(tmp_path / "a.n5"), str(tmp_path / "o.zarr"), str(flag))
    with pytest.raises(n5.N5Error, match="does not exist, skip the processing"):
        n5.n5_to_zarr_block(*args, [[0, 0, 0, 3, 4, 5]], [[0, 0, 0, 3, 4, 5]])
    flag.parent.mkdir()
    flag.write_bytes(b"")
    assert n5.n5_to_zarr_block(*args, [[0, 0, 0, 3, 4, 5]], [[0, 0, 0, 3, 4, 5]]) == []
    assert "The batch files (1 - 1) already exist, skip them!" in capsys.readouterr().out


def test_matlab_names_convert_matlabs_values(tmp_path):
    pytest.importorskip("numcodecs")
    from petakit5d import matlab

    a = volume((5, 6, 7), "uint16", "names")
    synth.write_zarr(tmp_path / "t.zarr", a, (2, 2, 2), "matlab")
    # MATLAB's [y x z] block size and its 1-based inclusive box
    matlab.zarrToN5(
        str(tmp_path / "t.zarr"),
        "BlockSize",
        [2, 3, 4],
        "bbox",
        [1, 2, 1, 6, 7, 5],
        "parseCluster",
        False,
    )
    arr = n5.open_n5(tmp_path / "n5" / "t.n5")
    assert (arr.block_size, arr.dims) == ((4, 2, 3), (6, 6, 5))
    np.testing.assert_array_equal(n5.read_n5(arr), a[:, :, 1:])


# ---------------------------------------------------------------------------------------------
# the XML pieces
# ---------------------------------------------------------------------------------------------


def test_write_struct_xml(tmp_path):
    s = {
        "versionAttribute": 0.2,
        "A": {"typeAttribute": "x&y", "Text": "a<b"},
        "B": [{"id": 0, "name": 1.0}, {"id": 2, "name": None}],
        "C": None,
        "D": "",
        "E": {"onlyAttribute": 1},
    }
    n5.write_struct_xml(s, tmp_path / "s.xml", "Root")
    assert (tmp_path / "s.xml").read_bytes() == (
        b'<?xml version="1.0" encoding="UTF-8"?>\n<Root version="0.2">\n'
        b'    <A type="x&amp;y">\n        <Text>a&lt;b</Text>\n    </A>\n'
        b"    <B>\n        <id>0</id>\n        <name>1</name>\n    </B>\n"
        b'    <B>\n        <id>2</id>\n    </B>\n    <D></D>\n    <E only="1"/>\n</Root>\n'
    )


def test_remove_text_tags_and_clean_xml_file(tmp_path):
    f = tmp_path / "a.xml"
    f.write_bytes(
        b'<r>\n  <BasePath type="relative">\n    <Text>.</Text>\n  </BasePath>\n'
        b"  <x><Text>keep</Text> tail</x>\n  <Text>open\n  close</Text>\n"
        b'  <n5 type="relative">\n\t<Text> a\\nb.n5\\t </Text>\n  </n5>\n</r>\n'
    )
    n5.remove_text_tags(f)
    assert f.read_bytes() == (
        b'<r>\n  <BasePath type="relative">\n    .\n  </BasePath>\n'
        b"  <x><Text>keep</Text> tail</x>\n  <Text>open\n  close</Text>\n"
        b'  <n5 type="relative">\n\t a\\nb.n5\\t \n  </n5>\n</r>\n'
    )
    n5.clean_xml_file(f)
    # white space around the texts goes, and MATLAB's '\n' and '\t' (backslash pairs)
    assert f.read_bytes() == (
        b'<r>\n  <BasePath type="relative">.</BasePath>\n'
        b"  <x><Text>keep</Text> tail</x>\n  <Text>open\n  close</Text>\n"
        b'  <n5 type="relative">ab.n5</n5>\n</r>\n'
    )


def test_xml_of_n5_in_xyz_order_and_zyx_shifts(tmp_path):
    d = str(tmp_path / "d") + "/"
    xml_tree(kind="nzarr")(d)
    with pytest.warns(UserWarning, match="Some files in the image list file do not exist"):
        s, t, vs = n5.construct_xml_file(
            d, d + "out.n5", d + "ImageList.csv", (0.3, 0.2, 0.108), ["CamB"], True, "xyz", "zyx",
            d + "xcorr.mat",
        )  # fmt: skip
    assert vs == (0.108, 0.2, 0.3)
    setups = s["SequenceDescription"]["ViewSetups"]["ViewSetup"]
    assert [v["size"] for v in setups] == ["5 6 4"] * 4  # [nx ny nz]
    assert t.height == 4
    reg = s["ViewRegistrations"]["ViewRegistration"][0]["ViewTransform"]
    # dxyz_shift(:, [3 1 2]) along [x y z]
    assert reg[2]["affine"] == ("1.0 0.0 0.0 27.777778 0.0 1.0 0.0 7.500000 0.0 0.0 1.0 -6.666667")


@needs_numcodecs
def test_xml_describes_the_n5_zarr_to_n5_writes(tmp_path):
    d = str(tmp_path / "d") + "/"
    xml_tree(kind="nzarr", xcorr=False)(d)
    with pytest.warns(UserWarning, match="Some files in the image list file do not exist"):
        s, t, vs = n5.construct_xml_file(
            d, d + "out.n5", d + "ImageList.csv", (0.3, 0.2, 0.108), ["CamB"], True
        )
    assert vs == (0.108, 0.2, 0.3)  # [x y z]: the default data_order 'zyx'
    setups = s["SequenceDescription"]["ViewSetups"]["ViewSetup"]
    for fn, setup in zip(t["Filename"], setups, strict=True):
        n5.zarr_to_n5(d + fn, block_size=(2, 4, 4), parse_cluster=False)
        arr = n5.open_n5(d + "n5/" + fn.replace(".zarr", ".n5"))
        assert setup["size"] == " ".join(str(v) for v in arr.dims)
