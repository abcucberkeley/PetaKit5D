"""Maximum intensity projections (petakit5d.tools.mip, task T02).

Oracle cases against MATLAB: saveMIP_tiff, saveMIP_zarr, saveMIP_frame, XR_MIP_wrapper (tiff,
zarr and large zarr datasets), XR_MIP_zarr (slabs only, pool sizes, a box, existing slabs, the
merge of the blocks' MIPs), MIP_zarr_batch_size, MIP_zarr_merge_memory, MIP_zarr_block_mips,
MIP_zarr_save_block_mips and MIP_zarr_merge_block_mips, for uint8, uint16 and single (NaN with
two payloads, -0 and 0, values a cast to uint16 rounds and saturates), sizes not divisible by
the batch, chunk or pool sizes, each axis choice, thin and 2D images. The MIP tiffs are
compared byte for byte (MATLAB's orientation), the slab zarrs by data and metadata after the
layout mapping, the flags and parameter files by content. Inputs are MATLAB's: tiffs, zarr
arrays in MATLAB's layout and flags in MATLAB's format, read by both sides.

Then MIPTest and BlockMIPTest (tests/tools) in intent on the port's (z, y, x) data and native
zarr arrays, the publish order of the large-file path (spy), the job shapes (framework spy),
the nice of large-file jobs (fake Slurm), and results that don't depend on threads or on the
zarr layout.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.dtypes import NAN, cast
from petakit5d.io import (
    create_zarr,
    read_mat,
    read_tiff,
    read_zarr,
    write_mat,
    write_tiff,
    write_zarr,
)
from petakit5d.testing import oracle, slurm_fake, synth
from petakit5d.testing.cases import (
    Box,
    Case,
    Dim,
    File,
    In,
    InPath,
    PerAxis,
    Tree,
    Value,
    Vol,
    Work,
)
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import mip

# ---------------------------------------------------------------------------------------------
# data: (z, y, x) volumes made from counter-based random numbers (bit-stable everywhere)
# ---------------------------------------------------------------------------------------------

NAN32 = np.float32(NAN)
# a second NaN payload (MATLAB: typecast(uint32(hex2dec('7FC00001')), 'single'))
NAN32_B = np.array([0x7FC00001], dtype=np.uint32).view(np.float32)[0]
NEG0 = np.float32(-0.0)


def volume(shape, cls, name, special=False):
    """a deterministic (z, y, x) volume: integers in [0, max], single normal * 100; with
    special (single): NaN of two payloads, -0, a line of NaN along x and z, lines of -0 and 0
    (BlockMIPTest's makeImage)"""
    rng = CounterRNG(0, f"t02/{name}")
    shape = tuple(int(s) for s in shape)
    if cls == "uint8":
        return rng.integers(0, 255, shape, dtype=np.uint8)
    if cls == "uint16":
        return rng.integers(0, 60000, shape, dtype=np.uint16)
    a = np.ascontiguousarray((rng.normal(shape) * 100.0).astype(np.float32))
    if special:
        u = rng.uniform((3, *shape))
        a[u[0] < 0.05] = NAN32
        a[u[1] < 0.05] = NEG0
        a[u[2] < 0.05] = NAN32_B
        if shape[1] > 1:
            a[:, 1, :] = NAN32
        if shape[2] > 2:
            a[:, :, 2] = NEG0
            a[1::3, :, 2] = 0
        if shape[1] > 3 and shape[2] > 4:
            a[:, 3, 4] = NEG0
            a[-1, 3, 4] = 0
    return a


def to_cast(shape, name):
    """single values a cast to uint16 rounds (x.5), saturates (below 0, above 65535) and maps
    NaN to 0, so a MIP cast differs from a cast's MIP nowhere"""
    rng = CounterRNG(0, f"t02/cast/{name}")
    a = (rng.uniform(shape) * 80000.0 - 5000.0).astype(np.float32)
    a[rng.uniform(shape) < 0.1] = np.float32(65535.5)
    a[rng.uniform(shape) < 0.05] = NAN32
    flat = a.reshape(-1)
    flat[::7] = np.floor(flat[::7]) + np.float32(0.5)
    return a


def ref_max(a, axis):
    """MATLAB's max along an axis, folded in order: NaN ignored (a line of NaN keeps its first
    NaN), of equal values the first"""
    a = np.moveaxis(np.asarray(a), axis, 0)
    best = a[0].copy()
    for s in a[1:]:
        if best.dtype.kind == "f":
            take = (s > best) | (np.isnan(best) & ~np.isnan(s))
        else:
            take = s > best
        best[take] = s[take]
    return np.expand_dims(best, axis)


def same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype, (a.dtype, b.dtype)
    assert a.shape == b.shape, (a.shape, b.shape)
    if a.dtype.kind == "f":
        assert a.view(f"u{a.itemsize}").tobytes() == b.view(f"u{b.itemsize}").tobytes()
    else:
        assert np.array_equal(a, b)


def tiff_image(m):
    """the MIP tiff (MATLAB's squeeze, MATLAB's orientation) of a projection (axis kept)"""
    m = to_matlab(m)
    s = list(m.shape)
    while len(s) > 2 and s[-1] == 1:  # MATLAB's size: no trailing singletons
        s.pop()
    if len(s) > 2:  # squeeze
        s = [d for d in s if d != 1]
        s += [1] * (2 - len(s))
    return np.ascontiguousarray(m.reshape(s, order="F"))


AXES = {"y": -2, "x": -1, "z": -3}

# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

M = "petakit5d.matlab"


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def py_save_mip_tiff(*args):
    return call_matlab("saveMIP_tiff", *args)


def py_save_mip_zarr(d, *args):
    return call_matlab("saveMIP_zarr", d + "im.zarr", *args)


def py_wrapper(*args):
    return call_matlab("XR_MIP_wrapper", *args)


def py_mip_zarr(d, *args):
    return call_matlab("XR_MIP_zarr", d + "im.zarr", *args)


def py_mip_zarr_flags(d, names, *args):
    return call_matlab(
        "XR_MIP_zarr", d + "im.zarr", *args, "blockFlagFullpaths", [d + n for n in names]
    )


def tif_file(shape, cls, name, special=False, maker=None):
    def make(path):
        synth.write_tiff(path, maker(shape, name) if maker else volume(shape, cls, name, special))

    return make


def tif_2d(shape, cls, name):
    def make(path):
        synth.write_tiff(path, volume((1, *shape), cls, name)[0])

    return make


def zarr_tree(shape, cls, name, chunks, special=False, maker=None, fname="im.zarr"):
    """a folder with one zarr array in MATLAB's layout (chunks (z, y, x))"""

    def make(folder):
        a = maker(shape, name) if maker else volume(shape, cls, name, special)
        synth.write_zarr(Path(folder) / fname, a, block_size=chunks, layout="matlab")

    return make


DATASET = ("scan_CamA_ch0_t0000", "scan_CamB_ch0_t0000", "scan_CamA_ch1_t0000")


def dataset_tree(fmt, cls="uint16", shape=(12, 40, 36), chunks=(4, 16, 16), special=False):
    """MIPTest's writeDataset: three files, the third matching no channel pattern of the tests"""

    def make(folder):
        for k, n in enumerate(DATASET):
            a = volume(shape, cls, f"ds/{fmt}/{cls}/{k}", special)
            if fmt == "zarr":
                synth.write_zarr(Path(folder) / f"{n}.zarr", a, block_size=chunks, layout="matlab")
            else:
                synth.write_tiff(Path(folder) / f"{n}.tif", a)

    return make


def tif_case(cid, make, *args):
    return Case(
        cid,
        matlab="saveMIP_tiff",
        nargout=0,
        args=[InPath("im"), Work("out_MIP_z.tif"), *args],
        python=py_save_mip_tiff,
        inputs={"im": File(make, "im.tif")},
    )


def zarr_case(cid, make, *args):
    return Case(
        cid,
        code="saveMIP_zarr([in{1}, 'im.zarr'], in{2 : end}); out = {};",
        nargout=0,
        args=[InPath("d"), Work("out_MIP_z.tif"), *args],
        python=py_save_mip_zarr,
        inputs={"d": Tree(make, "d")},
    )


def wrapper_case(cid, make, *args):
    return Case(
        cid,
        code="XR_MIP_wrapper(in{:}); out = {};",
        nargout=0,
        args=[InPath("d"), *args, "parseCluster", False, "masterCompute", True],
        python=py_wrapper,
        inputs={"d": Tree(make, "d")},
    )


def mip_zarr_case(cid, make, *args, error=False):
    return Case(
        cid,
        code="XR_MIP_zarr([in{1}, 'im.zarr'], in{2 : end}); out = {};",
        nargout=0,
        args=[InPath("d"), *args, "parseCluster", False],
        python=py_mip_zarr,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=(mip.MipError, ValueError, OSError, RuntimeError),
    )


def _tif_cases():
    s = (23, 37, 29)  # MATLAB 37 x 29 x 23
    return [
        tif_case("tif.uint16.all", tif_file(s, "uint16", "t1"), "axis", [1, 1, 1]),
        tif_case("tif.uint8.all", tif_file(s, "uint8", "t2"), "axis", [1, 1, 1]),
        tif_case(
            "tif.single.all",
            tif_file(s, "single", "t3", True),
            "axis",
            [1, 1, 1],
            "dtype",
            "single",
        ),
        tif_case(
            "tif.single_to_uint16.all",
            tif_file(s, "single", "t4", maker=to_cast),
            "axis",
            [1, 1, 1],
        ),
        tif_case("tif.uint16.y_and_z", tif_file(s, "uint16", "t5"), "axis", [1, 0, 1]),
        tif_case("tif.uint16.x", tif_file(s, "uint16", "t6"), "axis", [0, 1, 0]),
        tif_case("tif.uint16.default_z", tif_file(s, "uint16", "t7")),
        tif_case(
            "tif.uint16.bbox",
            tif_file(s, "uint16", "t8"),
            "axis",
            [1, 1, 1],
            "inputBbox",
            [5, 7, 3, 30, 25, 20],
        ),
        tif_case(
            "tif.single.bbox_one_slice",
            tif_file(s, "single", "t9", True),
            "axis",
            [1, 1, 1],
            "inputBbox",
            [2, 3, 7, 30, 25, 7],
            "dtype",
            "single",
        ),
        tif_case("tif.uint16.thin_x", tif_file((23, 37, 1), "uint16", "t10"), "axis", [1, 1, 1]),
        tif_case("tif.uint16.thin_y", tif_file((23, 1, 29), "uint16", "t11"), "axis", [1, 1, 1]),
        tif_case("tif.uint16.line_z", tif_file((23, 1, 1), "uint16", "t12"), "axis", [1, 1, 1]),
        tif_case("tif.uint16.2d", tif_2d((37, 29), "uint16", "t13"), "axis", [1, 1, 1]),
        tif_case(
            "tif.single.2d", tif_2d((37, 29), "single", "t14"), "axis", [1, 1, 1], "dtype", "single"
        ),
    ]


def _zarr_cases():
    big = (12, 128, 128)  # 64 chunks of 16 x 16 per layer of 4: 3 slabs
    return [
        zarr_case(
            "zarr.uint16.slabs", zarr_tree(big, "uint16", "z1", (4, 16, 16)), "axis", [1, 1, 1]
        ),
        zarr_case(
            "zarr.single.slabs",
            zarr_tree(big, "single", "z2", (4, 16, 16), True),
            "axis",
            [1, 1, 1],
        ),
        zarr_case(
            "zarr.single.slabs_bbox",
            zarr_tree(big, "single", "z3", (4, 16, 16), True),
            "axis",
            [1, 1, 1],
            "inputBbox",
            [5, 7, 3, 120, 110, 10],
        ),
        zarr_case("zarr.uint16.slabs_default_z", zarr_tree(big, "uint16", "z4", (4, 16, 16))),
        zarr_case(
            "zarr.uint16.uneven_slabs",
            zarr_tree((11, 64, 64), "uint16", "z5", (2, 16, 16)),
            "axis",
            [1, 1, 1],
        ),
        zarr_case(
            "zarr.uint8.one_read",
            zarr_tree((12, 40, 36), "uint8", "z6", (4, 16, 16)),
            "axis",
            [1, 1, 1],
        ),
        zarr_case(
            "zarr.uint16.as_single",
            zarr_tree(big, "uint16", "z7", (4, 16, 16)),
            "axis",
            [1, 1, 1],
            "dtype",
            "single",
        ),
        zarr_case(
            "zarr.single.to_uint16",
            zarr_tree((12, 40, 36), "single", "z8", (4, 16, 16), maker=to_cast),
            "axis",
            [0, 1, 1],
            "dtype",
            "uint16",
        ),
    ]


def _frame_cases():
    cases = []
    for cid, shape, cls, special, args in (
        ("frame.single_to_uint16.all", (11, 17, 13), "single", False, ["axis", [1, 1, 1]]),
        ("frame.single.all", (11, 17, 13), "single", True, ["axis", [1, 1, 1], "dtype", "single"]),
        ("frame.uint16.default_z", (11, 17, 13), "uint16", False, []),
        ("frame.uint8.x", (11, 17, 13), "uint8", False, ["axis", [0, 1, 0]]),
        ("frame.uint16.none", (11, 17, 13), "uint16", False, ["axis", [0, 0, 0]]),
    ):
        maker = to_cast if cid.startswith("frame.single_to") else None
        cases.append(
            Case(
                cid,
                matlab="saveMIP_frame",
                nargout=0,
                args=[Vol(In("im")), Work("sub/out_MIP_z.tif"), *args],
                python=_py_frame,
                inputs={
                    "im": Value(
                        lambda s=shape, c=cls, n=cid, sp=special, mk=maker: to_matlab(
                            mk(s, n) if mk else volume(s, c, n, sp)
                        )
                    )
                },
            )
        )
    return cases


def _py_frame(im, path, *args):
    """saveMIP_frame's name-value pairs, MATLAB's [y x z] axis as (z, y, x)"""
    from petakit5d.core.axes import from_yxz

    kw = dict(zip(args[0::2], args[1::2], strict=True))
    if "axis" in kw:
        kw["axis"] = tuple(from_yxz(tuple(kw["axis"])))
    mip.save_mip_frame(im, path, **kw)
    return []


def _wrapper_cases():
    chans = ["channelPatterns", ["CamA_ch0", "CamB_ch0"]]
    one = ["channelPatterns", ["CamA_ch0"]]
    cases = [
        wrapper_case("wrapper.tif.uint16.all", dataset_tree("tif"), "axis", [1, 1, 1], *chans),
        wrapper_case("wrapper.tif.uint16.default_z", dataset_tree("tif"), *one),
        wrapper_case(
            "wrapper.tif.single.save16bit_false",
            dataset_tree("tif", "single", special=True),
            "axis",
            [1, 1, 1],
            "save16bit",
            False,
            *one,
        ),
        wrapper_case(
            "wrapper.tif.single.to_uint16",
            dataset_tree("tif", "single", special=True),
            "axis",
            [1, 1, 1],
            *one,
        ),
        wrapper_case(
            "wrapper.tif.uint16.bbox",
            dataset_tree("tif"),
            "axis",
            [1, 1, 1],
            "inputBbox",
            [5, 7, 3, 30, 33, 10],
            *one,
        ),
        wrapper_case(
            "wrapper.zarr.uint16.all",
            dataset_tree("zarr"),
            "axis",
            [1, 1, 1],
            "zarrFile",
            True,
            *chans,
        ),
        wrapper_case(
            "wrapper.large.uint16.all",
            dataset_tree("zarr"),
            "axis",
            [1, 1, 1],
            "zarrFile",
            True,
            "largeFile",
            True,
            "batchSize",
            [16, 16, 4],
            *one,
        ),
    ]
    # the large-file path: each axis choice, sizes not divisible by the batch size, classes
    for cls, special in (("uint8", False), ("uint16", False), ("single", True)):
        for ax in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0]):
            name = "".join(c for c, f in zip("yxz", ax, strict=True) if f)
            cases.append(
                wrapper_case(
                    f"wrapper.large.{cls}.{name}",
                    dataset_tree(
                        "zarr", cls, shape=(23, 37, 29), chunks=(5, 7, 6), special=special
                    ),
                    "axis",
                    ax,
                    "zarrFile",
                    True,
                    "largeFile",
                    True,
                    "batchSize",
                    [16, 12, 10],
                    *one,
                )
            )
    return cases


def _mip_zarr_cases():
    s = (23, 37, 29)
    u16 = zarr_tree(s, "uint16", "mz1", (5, 7, 6))
    sgl = zarr_tree(s, "single", "mz2", (5, 7, 6), True)
    cases = [
        mip_zarr_case("mip_zarr.uint16.default", u16, "batchSize", [16, 12, 10]),
        mip_zarr_case("mip_zarr.single.default", sgl, "batchSize", [16, 12, 10]),
        mip_zarr_case(
            "mip_zarr.uint16.slabs_only", u16, "batchSize", [16, 12, 10], "mipSlab", True
        ),
        mip_zarr_case(
            "mip_zarr.uint16.one_batch", u16, "axis", [0, 1, 1], "batchSize", [64, 64, 64]
        ),
        mip_zarr_case(
            "mip_zarr.uint16.pool3", u16, "batchSize", [16, 12, 10], "poolSize", [4, 3, 2]
        ),
        mip_zarr_case(
            "mip_zarr.single.pool3", sgl, "batchSize", [16, 12, 10], "poolSize", [4, 3, 2]
        ),
        mip_zarr_case(
            "mip_zarr.uint16.pool6",
            u16,
            "batchSize",
            [16, 12, 10],
            "poolSize",
            [4, 6, 4, 2, 2, 1],
        ),
        mip_zarr_case(
            "mip_zarr.uint16.pool6_not_multiple",
            u16,
            "batchSize",
            [16, 12, 10],
            "poolSize",
            [4, 3, 4, 2, 2, 1],
        ),
        mip_zarr_case(
            "mip_zarr.uint16.pool9_no_resampling",
            u16,
            "batchSize",
            [16, 12, 10],
            "poolSize",
            [4, 6, 4, 2, 2, 1, 1, 1, 1],
        ),
        mip_zarr_case(
            "mip_zarr.uint16.bbox",
            u16,
            "batchSize",
            [16, 12, 10],
            "inputBbox",
            [3, 2, 4, 30, 25, 20],
        ),
        # boxes that end before the image's end along one axis or all, with batch remainders
        # (batches 21 x 12 x 10 here: MATLAB's last batches end at the box's end since 11c7d0b)
        mip_zarr_case(
            "mip_zarr.uint16.bbox_end_y",
            u16,
            "batchSize",
            [16, 12, 10],
            "inputBbox",
            [1, 1, 1, 30, 29, 23],
        ),
        mip_zarr_case(
            "mip_zarr.uint16.bbox_end_x",
            u16,
            "batchSize",
            [16, 12, 10],
            "inputBbox",
            [5, 1, 1, 37, 20, 23],
        ),
        mip_zarr_case(
            "mip_zarr.uint16.bbox_end_z",
            u16,
            "batchSize",
            [16, 12, 10],
            "inputBbox",
            [1, 2, 4, 37, 29, 17],
        ),
        mip_zarr_case(
            "mip_zarr.single.bbox_end_all",
            sgl,
            "batchSize",
            [16, 12, 10],
            "inputBbox",
            [2, 3, 2, 25, 17, 13],
        ),
        mip_zarr_case(
            "mip_zarr.uint16.bbox_slabs_only",
            u16,
            "batchSize",
            [16, 12, 10],
            "inputBbox",
            [3, 2, 4, 30, 25, 20],
            "mipSlab",
            True,
        ),
        mip_zarr_case(
            "mip_zarr.uint16.bbox_to_the_end",
            u16,
            "batchSize",
            [16, 12, 10],
            "inputBbox",
            [3, 2, 4, 37, 29, 23],
        ),
        mip_zarr_case(
            "mip_zarr.uint16.existing_slabs",
            existing_slabs_tree(),
            "axis",
            [1, 1, 1],
            "batchSize",
            [16, 12, 10],
        ),
    ]
    for cls, special in (("uint16", False), ("single", True)):
        cases.append(
            Case(
                f"mip_zarr.{cls}.merge_flags",
                code=(
                    "fns = cellfun(@(x) [in{1}, x], in{2}, 'unif', 0); "
                    "XR_MIP_zarr([in{1}, 'im.zarr'], in{3 : end}, 'parseCluster', false, "
                    "'blockFlagFullpaths', fns); out = {};"
                ),
                nargout=0,
                args=[InPath("d"), FLAG_NAMES, "batchSize", [8, 6, 5]],
                python=py_mip_zarr_flags,
                py_args=[InPath("d"), FLAG_NAMES, "batchSize", [8, 6, 5], "parseCluster", False],
                inputs={"d": Tree(flags_tree(cls, special, with_zarr=True), "d")},
            )
        )
    cases.append(
        Case(
            "mip_zarr.uint16.merge_flag_without_mips",
            code=(
                "fns = cellfun(@(x) [in{1}, x], in{2}, 'unif', 0); "
                "XR_MIP_zarr([in{1}, 'im.zarr'], in{3 : end}, 'parseCluster', false, "
                "'blockFlagFullpaths', fns); out = {};"
            ),
            nargout=0,
            args=[InPath("d"), FLAG_NAMES, "batchSize", [8, 6, 5]],
            python=py_mip_zarr_flags,
            py_args=[InPath("d"), FLAG_NAMES, "batchSize", [8, 6, 5], "parseCluster", False],
            inputs={"d": Tree(flags_tree("uint16", False, True, wrong=True, drop_last=True), "d")},
        )
    )
    return cases


# BlockMIPTest's image, MIP batches and blocks not aligned with them (MATLAB [37 29 23],
# [8 6 5], [7 9 4]), and its flags (4 blocks a task, blocks in reverse order)
FLAG_SHAPE = (23, 37, 29)
FLAG_MIP_BATCH = (5, 8, 6)
FLAG_BLOCK = (4, 7, 9)


def _flag_blocks():
    nz, ny, nx = FLAG_SHAPE
    bz, by, bx = FLAG_BLOCK
    boxes = []
    for iz in range(-(-nz // bz)):
        for ix in range(-(-nx // bx)):
            for iy in range(-(-ny // by)):
                s = (iz * bz, iy * by, ix * bx)
                boxes.append((*s, min(s[0] + bz, nz), min(s[1] + by, ny), min(s[2] + bx, nx)))
    return boxes[::-1]


def _flag_names():
    n = len(_flag_blocks())
    return [f"flags/blocks_{s + 1}_{min(s + 4, n)}.mat" for s in range(0, n, 4)]


FLAG_NAMES = _flag_names()
_MAT_HEADER = b"MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Thu Jan  1 00:00:00 2026"


def fix_mat_header(path):
    """a MAT v5 file's header text without the platform and time it was written on (inputs
    the same everywhere)"""
    data = bytearray(Path(path).read_bytes())
    assert data[:10] == b"MATLAB 5.0"
    data[:116] = _MAT_HEADER.ljust(116, b" ")
    Path(path).write_bytes(bytes(data))


def write_flags(folder, im, cls, wrong=False, drop_last=False):
    """the block tasks' flags of im (z, y, x) as MATLAB's MIP_zarr_save_block_mips saves them"""
    os.makedirs(folder, exist_ok=True)
    boxes = _flag_blocks()
    src = im + 1 if wrong else im
    names = []
    for s in range(0, len(boxes), 4):
        bb = boxes[s : s + 4]
        mips = [
            mip.block_mips(src[b[0] : b[3], b[1] : b[4], b[2] : b[5]], b, FLAG_MIP_BATCH, cls)
            for b in bb
        ]
        fn = Path(folder) / f"blocks_{s + 1}_{min(s + 4, len(boxes))}.mat"
        mip.save_block_mips(str(fn), bb, mips, FLAG_MIP_BATCH, 1)
        fix_mat_header(fn)
        names.append(fn)
    if drop_last:
        write_mat(names[-1], {"t": 1.0}, version="6")
        fix_mat_header(names[-1])
    return names


def flags_tree(cls, special, with_zarr=False, wrong=False, drop_last=False):
    def make(folder):
        im = volume(FLAG_SHAPE, cls, f"flags/{cls}", special)
        write_flags(Path(folder) / "flags", im, cls, wrong=wrong, drop_last=drop_last)
        if with_zarr:
            synth.write_zarr(Path(folder) / "im.zarr", im, block_size=(5, 4, 3), layout="matlab")

    return make


def existing_slabs_tree():
    """a zarr with slab zarrs of an earlier run (MATLAB's layout) and its z MIP"""

    def make(folder):
        im = volume((23, 37, 29), "uint16", "mz/slabs")
        synth.write_zarr(Path(folder) / "im.zarr", im, block_size=(5, 7, 6), layout="matlab")
        os.makedirs(Path(folder) / "MIPs")
        for lab, ax in AXES.items():
            # an earlier run's slab: here a constant, so the MIP shows where it came from
            shape = list(im.shape)
            shape[ax] = 2
            slab = np.full(shape, 7 + ord(lab), dtype=np.uint16)
            synth.write_zarr(
                Path(folder) / "MIPs" / f"im_{lab}.zarr",
                slab,
                block_size=(5, 7, 6),
                layout="matlab",
            )

    return make


def _batch_size_cases():
    cases = []
    for k, (ins, blk, cls, req, slab) in enumerate(
        (
            ([5000, 3000, 900], [256, 256, 256], "uint16", None, None),
            ([5000, 3000, 900], [256, 300, 100], "uint16", None, None),
            ([5000, 3000, 900], [256, 256, 256], "uint16", [100, 100, 100], None),
            ([100000, 100000, 11000], [256, 256, 256], "uint16", None, None),
            ([200000, 100000, 10000], [256, 256, 256], "uint16", None, None),
            ([200000, 100000, 10000], [256, 256, 256], "uint16", None, True),
            ([150000, 150000, 30000], [512, 256, 128], "single", None, None),
            ([90000, 90000, 90000], [100, 128, 300], "uint8", [1000, 700, 2000], None),
            ([37, 29, 23], [7, 9, 4], "single", [8, 6, 5], None),
            ([60000, 60000, 60000], [256, 256, 256], "uint16", [512, 512, 512], None),
        )
    ):
        args = [PerAxis(ins, "yxz"), PerAxis(blk, "yxz"), cls]
        kw = {}
        if req is not None:
            args += ["batchSize", PerAxis(req, "yxz")]
            kw["batch_size"] = PerAxis(req, "yxz")
        if slab is not None:
            args += ["mipSlab", slab]
            kw["mip_slab"] = slab
        cases.append(
            Case(
                f"batch_size.{k}",
                matlab="MIP_zarr_batch_size",
                args=args,
                python=mip.mip_zarr_batch_size,
                py_args=args[:3],
                py_kwargs=kw,
                out_roles=("yxz",),
                py_out=lambda r: [np.asarray(r, dtype=np.float64)],
            )
        )
    for k, (ims, reg, blk, cls) in enumerate(
        (
            ([1800, 22651, 270], [1024, 22651, 270], [256, 256, 256], "uint16"),
            ([37, 29, 23], [7, 9, 4], [5, 4, 3], "single"),
            ([200000, 100000, 10000], [2048, 2048, 2048], [256, 256, 256], "uint16"),
            ([5000, 3000, 900], [1000, 4000, 200], [256, 300, 100], "uint8"),
        )
    ):
        args = [PerAxis(ims, "yxz"), PerAxis(reg, "yxz"), PerAxis(blk, "yxz"), cls]
        cases.append(
            Case(
                f"merge_memory.{k}",
                matlab="MIP_zarr_merge_memory",
                args=args,
                python=mip.merge_memory,
            )
        )
    return cases


def _py_block_mips(block, box, mbs, cls):
    z, y, x = mip.block_mips(block, box, mbs, cls)
    return [y, x, z]  # (a list: one output, a cell)


def _block_mips_cases():
    cases = []
    specs = (
        ("uint16", (13, 20, 15), [3, 2, 4, 22, 16, 16], [8, 6, 5], None, False),
        ("single", (13, 20, 15), [3, 2, 4, 22, 16, 16], [8, 6, 5], None, True),
        ("single_to_uint16", (9, 12, 10), [1, 1, 1, 12, 10, 9], [5, 4, 4], "uint16", False),
        ("uint8", (13, 20, 15), [9, 7, 6, 28, 21, 18], [4, 3, 3], None, False),
        ("uint16.one_batch", (13, 20, 15), [1, 1, 1, 20, 15, 13], [64, 64, 64], None, False),
        ("single.thin", (1, 20, 15), [3, 2, 4, 22, 16, 4], [8, 6, 5], None, True),
    )
    for name, shape, box, mbs, out_cls, special in specs:
        cls = name.split(".")[0].split("_to_")[0]
        maker = to_cast if "_to_" in name else None
        cases.append(
            Case(
                f"block_mips.{name}",
                matlab="MIP_zarr_block_mips",
                args=[Vol(In("b")), Box(box), PerAxis(mbs, "yxz"), out_cls or cls],
                python=_py_block_mips,
                out_roles=(("each", "vol"),),
                inputs={
                    "b": Value(
                        lambda s=shape, c=cls, n=name, sp=special, mk=maker: to_matlab(
                            mk(s, n) if mk else volume(s, c, f"bm/{n}", sp)
                        )
                    )
                },
            )
        )
    return cases


def _py_save_flag(path, boxes, cell, mbs, t):
    """MIP_zarr_save_block_mips with MATLAB's boxes, cell (rows: blocks; y, x, z) and batch"""
    from petakit5d.core.axes import from_yxz
    from petakit5d.core.bbox import boxes_from_matlab

    rows = boxes_from_matlab(np.asarray(boxes, dtype=np.float64).reshape(-1, 6))
    mips = [
        [
            from_matlab(np.asarray(r[2])),
            from_matlab(np.asarray(r[0])),
            from_matlab(np.asarray(r[1])),
        ]
        for r in cell
    ]
    mip.save_block_mips(
        path, [tuple(int(v) for v in r) for r in rows], mips, tuple(from_yxz(tuple(mbs))), t
    )
    return []


def _save_flag_cases():
    im = volume((4, 6, 5), "uint16", "flag/one")
    box = (0, 0, 0, 4, 6, 5)
    z, y, x = mip.block_mips(im, box, (8, 8, 8), "uint16")
    cell1 = [[to_matlab(y), to_matlab(x), to_matlab(z)]]
    ims = volume((9, 11, 7), "single", "flag/two", True)
    b1, b2 = (0, 0, 0, 9, 5, 7), (0, 5, 0, 9, 11, 7)
    m1 = mip.block_mips(ims[:, :5], b1, (4, 4, 3), "single")
    m2 = mip.block_mips(ims[:, 5:], b2, (4, 4, 3), "single")
    cell2 = [[to_matlab(m[1]), to_matlab(m[2]), to_matlab(m[0])] for m in (m1, m2)]
    from petakit5d.core.bbox import box_to_matlab

    def mcell_rows(rows):
        from petakit5d.testing.cases import MCell, to_mvalue

        items = []
        for c in range(3):
            for r in rows:
                items.append(to_mvalue(r[c]))
        return MCell((len(rows), 3), items)

    return [
        Case(
            "save_flag.one_block",
            code="MIP_zarr_save_block_mips(in{:}); out = {};",
            nargout=0,
            args=[
                Work("blocks_1_1.mat"),
                list(box_to_matlab(box)),
                mcell_rows(cell1),
                [8, 8, 8],
                1.5,
            ],
            python=_py_save_flag,
            py_args=[Work("blocks_1_1.mat"), list(box_to_matlab(box)), cell1, [8, 8, 8], 1.5],
            tree=tree(ignore_vars=()),
        ),
        Case(
            "save_flag.two_blocks",
            code="MIP_zarr_save_block_mips(in{:}); out = {};",
            nargout=0,
            args=[
                Work("blocks_1_2.mat"),
                np.array([box_to_matlab(b1), box_to_matlab(b2)], dtype=np.float64),
                mcell_rows(cell2),
                [4, 3, 4],
                None,
            ],
            python=_py_save_flag,
            py_args=[
                Work("blocks_1_2.mat"),
                np.array([box_to_matlab(b1), box_to_matlab(b2)], dtype=np.float64),
                cell2,
                [4, 3, 4],
                None,
            ],
        ),
    ]


def _py_merge(d, names, axis, im_size, batch, cls, *max_mem):
    slab = mip.merge_block_mips([d + n for n in names], axis, im_size, batch, cls, *max_mem)
    return None if slab is None else to_matlab(slab)


MERGE_CODE = (
    "fns = cellfun(@(x) [in{1}, x], in{2}, 'unif', 0); "
    "out = {MIP_zarr_merge_block_mips(fns, in{3 : end})};"
)


def _merge_cases():
    cases = []
    ims = [37, 29, 23]
    for cls, special in (("uint16", False), ("single", True)):
        for dim in (1, 2, 3):
            cases.append(
                Case(
                    f"merge.{cls}.dim{dim}",
                    code=MERGE_CODE,
                    args=[
                        InPath("d"),
                        FLAG_NAMES,
                        Dim(dim),
                        PerAxis(ims, "yxz"),
                        PerAxis([8, 6, 5], "yxz"),
                        cls,
                    ],
                    python=_py_merge,
                    inputs={"d": Tree(flags_tree(cls, special), "d")},
                )
            )
    tr = Tree(flags_tree("uint16", False), "d")
    for cid, args in (
        ("merge.other_batch", [Dim(1), PerAxis(ims, "yxz"), PerAxis([16, 6, 5], "yxz"), "uint16"]),
        ("merge.other_class", [Dim(1), PerAxis(ims, "yxz"), PerAxis([8, 6, 5], "yxz"), "single"]),
        (
            "merge.other_size",
            [Dim(2), PerAxis([37, 29, 24], "yxz"), PerAxis([8, 6, 5], "yxz"), "uint16"],
        ),
        (
            "merge.over_memory",
            [Dim(1), PerAxis(ims, "yxz"), PerAxis([8, 6, 5], "yxz"), "uint16", 2.1e-5],
        ),
        (
            "merge.under_memory",
            [Dim(1), PerAxis(ims, "yxz"), PerAxis([8, 6, 5], "yxz"), "uint16", 2.2e-5],
        ),
        (
            "merge.unknown_memory",
            [Dim(1), PerAxis(ims, "yxz"), PerAxis([8, 6, 5], "yxz"), "uint16", NAN],
        ),
    ):
        cases.append(
            Case(
                cid,
                code=MERGE_CODE,
                args=[InPath("d"), FLAG_NAMES, *args],
                python=_py_merge,
                inputs={"d": tr},
            )
        )
    cases.append(
        Case(
            "merge.missing_flag",
            code=MERGE_CODE,
            args=[
                InPath("d"),
                FLAG_NAMES[1:],
                Dim(3),
                PerAxis(ims, "yxz"),
                PerAxis([8, 6, 5], "yxz"),
                "uint16",
            ],
            python=_py_merge,
            inputs={"d": tr},
        )
    )
    cases.append(
        Case(
            "merge.flag_with_only_t",
            code=MERGE_CODE,
            args=[
                InPath("d"),
                FLAG_NAMES,
                Dim(3),
                PerAxis(ims, "yxz"),
                PerAxis([8, 6, 5], "yxz"),
                "uint16",
            ],
            python=_py_merge,
            inputs={"d": Tree(flags_tree("uint16", False, drop_last=True), "d")},
        )
    )
    return cases


CASES = (
    _tif_cases()
    + _zarr_cases()
    + _frame_cases()
    + _wrapper_cases()
    + _mip_zarr_cases()
    + _batch_size_cases()
    + _block_mips_cases()
    + _save_flag_cases()
    + _merge_cases()
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# MIPTest in intent, on (z, y, x) data and native zarr arrays
# ---------------------------------------------------------------------------------------------

IM_SIZE = (12, 40, 36)  # MATLAB's ImSize [40 36 12]
CHANNELS = ("CamA_ch0", "CamB_ch0")


def write_dataset(folder, fmt, cls="uint16", shape=IM_SIZE, chunks=(4, 16, 16), special=False):
    """MIPTest's writeDataset with the port's files (native zarr arrays)"""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ims = []
    for k, n in enumerate(DATASET):
        a = volume(shape, cls, f"py/{fmt}/{cls}/{k}", special)
        ims.append(a)
        if fmt == "zarr":
            write_zarr(a, str(folder / f"{n}.zarr"), block_size=chunks)
        else:
            write_tiff(a, str(folder / f"{n}.tif"))
    return str(folder) + "/", ims, list(DATASET)


def verify_mips(base, im, dtype):
    for lab, ax in AXES.items():
        got = read_tiff(f"{base}_MIP_{lab}.tif")
        same_bits(got, tiff_image(cast(ref_max(im, ax), dtype)))


def test_all_axes_of_tiff(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "tif")
    mip.mip_dataset(path, axis=(1, 1, 1), channel_patterns=CHANNELS, parse_cluster=False)
    for i in range(2):
        verify_mips(f"{path}MIPs/{names[i]}", ims[i], "uint16")
    assert not os.path.exists(f"{path}MIPs/{names[2]}_MIP_z.tif"), "a file of no channel pattern"


def test_default_is_z_only(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "tif")
    mip.mip_dataset(path, channel_patterns=CHANNELS[:1], parse_cluster=False)
    base = f"{path}MIPs/{names[0]}"
    same_bits(read_tiff(base + "_MIP_z.tif"), ref_max(ims[0], -3)[0])
    assert not os.path.exists(base + "_MIP_x.tif")
    assert not os.path.exists(base + "_MIP_y.tif")


def test_single_output_without_save16bit(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "tif", "single", special=True)
    mip.mip_dataset(
        path, axis=(1, 1, 1), channel_patterns=CHANNELS[:1], save_16bit=False, parse_cluster=False
    )
    verify_mips(f"{path}MIPs/{names[0]}", ims[0], "single")


def test_all_axes_of_zarr(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "zarr")
    mip.mip_dataset(
        path, axis=(1, 1, 1), channel_patterns=CHANNELS, zarr_file=True, parse_cluster=False
    )
    for i in range(2):
        verify_mips(f"{path}MIPs/{names[i]}", ims[i], "uint16")


def test_all_axes_of_large_zarr(tmp_path):
    # batch-wise MIPs of a zarr file, merged across batches
    path, ims, names = write_dataset(tmp_path / "d", "zarr")
    mip.mip_dataset(
        path,
        axis=(1, 1, 1),
        channel_patterns=CHANNELS[:1],
        zarr_file=True,
        large_file=True,
        batch_size=(4, 16, 16),
        parse_cluster=False,
    )
    verify_mips(f"{path}MIPs/{names[0]}", ims[0], "uint16")
    for lab in AXES:
        assert os.path.isdir(f"{path}MIPs/{names[0]}_{lab}.zarr"), f"no {lab} slab"
    assert os.listdir(f"{path}MIPs/zarr_flag") == [], "the flag folder was left behind"


def test_input_bbox_projects_the_cropped_region(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d", "tif")
    b = (2, 4, 6, 10, 30, 33)  # MATLAB's [5, 7, 3, 30, 33, 10]
    mip.mip_dataset(
        path, axis=(1, 1, 1), input_bbox=b, channel_patterns=CHANNELS[:1], parse_cluster=False
    )
    verify_mips(f"{path}MIPs/{names[0]}", ims[0][b[0] : b[3], b[1] : b[4], b[2] : b[5]], "uint16")


def test_save_mip_tiff_writes_requested_axes(tmp_path):
    im = volume(IM_SIZE, "uint16", "req")
    write_tiff(im, str(tmp_path / "im.tif"))
    # MATLAB's axis [1, 0, 1]: y and z
    mip.save_mip_tiff(str(tmp_path / "im.tif"), str(tmp_path / "out_MIP_z.tif"), axis=(1, 1, 0))
    same_bits(read_tiff(str(tmp_path / "out_MIP_y.tif")), tiff_image(ref_max(im, -2)))
    same_bits(read_tiff(str(tmp_path / "out_MIP_z.tif")), ref_max(im, -3)[0])
    assert not os.path.exists(tmp_path / "out_MIP_x.tif")


def test_mip_tiffs_keep_matlabs_orientation(tmp_path):
    # the y MIP is MATLAB's nx x nz image, the x MIP ny x nz: transposes of NumPy's
    im = volume((5, 7, 3), "uint16", "orient")
    write_tiff(im, str(tmp_path / "im.tif"))
    mip.save_mip_tiff(str(tmp_path / "im.tif"), str(tmp_path / "o_MIP_z.tif"), axis=(1, 1, 1))
    y = read_tiff(str(tmp_path / "o_MIP_y.tif"))
    x = read_tiff(str(tmp_path / "o_MIP_x.tif"))
    assert y.shape == (3, 5)
    assert np.array_equal(y, im.max(axis=1).T)
    assert x.shape == (7, 5)
    assert np.array_equal(x, im.max(axis=2).T)
    assert np.array_equal(read_tiff(str(tmp_path / "o_MIP_z.tif")), im.max(axis=0))


@pytest.mark.parametrize(
    ("shape", "want"),
    [
        ((5, 1, 3), {"y": (3, 5), "x": (5, 1), "z": (1, 3)}),
        ((5, 4, 1), {"y": (5, 1), "x": (4, 5), "z": (4, 1)}),
        ((5, 1, 1), {"y": (5, 1), "x": (5, 1), "z": (1, 1)}),
        ((1, 4, 3), {"y": (1, 3), "x": (4, 1), "z": (4, 3)}),
    ],
)
def test_thin_mips_are_squeezed_as_matlab_squeezes(shape, want):
    im = volume(shape, "uint16", f"thin/{shape}")
    for lab, ax in AXES.items():
        got = mip._matlab_image(mip.axis_max(im, ax))
        assert got.shape == want[lab], lab
        same_bits(got, tiff_image(ref_max(im, ax)))


@pytest.mark.parametrize("threads", [1, 5])
def test_large_mip_images_are_turned_in_threads(threads, monkeypatch):
    # y and x MIPs from 4 MB are turned to MATLAB's orientation in row blocks in threads
    from petakit5d.core import threads as core_threads

    monkeypatch.setattr(core_threads, "available_cpus", lambda: threads)
    rng = CounterRNG(0, "t02/turn")
    for shape in ((1100, 1, 2101), (1037, 2203, 1)):
        m = rng.integers(0, 65535, shape, dtype=np.uint16)
        got = mip._matlab_image(m)
        assert got.flags.c_contiguous
        same_bits(got, tiff_image(m))


def test_slab_another_writer_made_is_kept(tmp_path):
    # a slab that another writer (of the same result) had put in place is kept, and the MIP
    # of its axis is made from it
    path, ims, names = write_dataset(tmp_path / "d", "zarr")
    os.makedirs(f"{path}MIPs")
    slab_fn = f"{path}MIPs/{names[0]}_y.zarr"
    other = np.full((IM_SIZE[0], 1, IM_SIZE[2]), 7, dtype=np.uint16)
    write_zarr(other, slab_fn)
    mip.mip_zarr(f"{path}{names[0]}.zarr", axis=(1, 1, 1), parse_cluster=False)
    same_bits(read_zarr(slab_fn), other)
    same_bits(
        read_tiff(f"{path}MIPs/{names[0]}_MIP_y.tif"),
        np.full((IM_SIZE[2], IM_SIZE[0]), 7, np.uint16),
    )
    same_bits(read_tiff(f"{path}MIPs/{names[0]}_MIP_z.tif"), ref_max(ims[0], -3)[0])
    left = [n for n in os.listdir(f"{path}MIPs") if n.startswith(f"{names[0]}_y_")]
    assert left == [], "the temporary slab was not removed"
    assert not [n for n in os.listdir(slab_fn) if n.endswith(".zarr")], "moved into the other's"


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls (MATLAB's FrameworkSpyFixture); nothing runs"""
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        return [False] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    return calls


fake_slurm = slurm_fake.fake_slurm  # (the fixture)


@pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
def test_large_zarr_jobs_go_after_their_batches(tmp_path, fake_slurm, framework_spy, monkeypatch):
    # the job of a large file waits for its MIP batch jobs: it asks to go after them (nice,
    # coordinator_nice); the jobs of other files don't
    from petakit5d.jobs import coordinator

    create_zarr(
        str(tmp_path / "scan_CamA_ch0_t0000.zarr"),
        data_size=(64, 512, 512),
        block_size=(64, 64, 64),
    )
    monkeypatch.setattr(coordinator, "coordinator_nice", lambda slurm=None: 1000)
    (tmp_path / "logs").mkdir()
    for large in (True, False):
        mip.mip_dataset(
            str(tmp_path),
            axis=(1, 1, 1),
            channel_patterns=("CamA_ch0",),
            zarr_file=True,
            large_file=large,
            parse_cluster=True,
            result_dir_name=f"MIPs_{int(large)}",
            job_log_dir=str(tmp_path / "logs"),
        )
    assert framework_spy[0]["nice"] == 1000
    assert framework_spy[1]["nice"] == 0
    # a file is done once the MIP of its last flagged axis exists
    assert framework_spy[0]["outputs"][0].endswith("scan_CamA_ch0_t0000_MIP_z.tif")


def test_large_zarr_jobs_ask_for_what_they_hold(tmp_path, framework_spy):
    # the job of a large file (mip_zarr) runs its MIP batches as jobs; it asks for the process
    # with one CPU, and if it computes batches itself (master_compute, the config's if it sets
    # it), for cpus_per_task and a batch's estimate (2048^3 uint16: 16 GiB times the port's
    # BATCH_MEMORY_FACTOR 1.2, 19.2 GB; MATLAB asks for 4 times, 64 GB). The input is a 34
    # G-voxel header. The process alone is PROCESS_MEMORY (1 GB in Python, MATLAB's 2 GB).
    from petakit5d.jobs.coordinator import PROCESS_MEMORY

    create_zarr(
        str(tmp_path / "scan_CamA_ch0_t0000.zarr"),
        data_size=(2048, 4096, 4096),
        block_size=(256, 256, 256),
    )
    off = tmp_path / "cpu_config.json"
    off.write_text('{"masterCompute": false}', newline="\n")
    for i, (mc, cfg) in enumerate([(True, ""), (False, ""), (True, str(off))]):
        mip.mip_dataset(
            str(tmp_path),
            axis=(1, 1, 1),
            channel_patterns=("CamA_ch0",),
            zarr_file=True,
            large_file=True,
            parse_cluster=False,
            master_compute=mc,
            config_file=cfg,
            result_dir_name=f"MIPs_{i}",
        )
    assert framework_spy[0]["cpus_per_task"] == 3
    assert framework_spy[0]["mem_allocate"] == pytest.approx(16 * 1.2, abs=1e-9)
    for i in (1, 2):
        assert framework_spy[i]["cpus_per_task"] == 1, i
        assert framework_spy[i]["mem_allocate"] == PROCESS_MEMORY, i
    assert all("master_compute" not in t.kwargs for t in framework_spy[0]["tasks"])
    assert all(t.kwargs["master_compute"] is False for t in framework_spy[1]["tasks"])


def test_batch_jobs_ask_for_1_2_times_their_batch(tmp_path, framework_spy):
    # the port's batch task holds one batch and its MIPs (peak 1.07 times a batch, measured):
    # it asks for 1.2 times the batch (MATLAB: 4 times, 2.5 in compiled mode), doubled at each
    # of the three tries; mcc_mode doesn't change it
    assert mip.BATCH_MEMORY_FACTOR == 1.2
    fn = str(tmp_path / "im.zarr")
    create_zarr(fn, data_size=(64, 512, 512), block_size=(64, 64, 64))
    for k, mcc in enumerate((False, True)):
        with pytest.raises(mip.MipError, match="not finished"):
            mip.mip_zarr(
                fn, batch_size=(64, 256, 256), parse_cluster=False, mcc_mode=mcc, uuid=f"u{k}"
            )
    gb = 64 * 256 * 256 * 2 / 2**30
    got = [c["mem_allocate"] for c in framework_spy]
    assert got == pytest.approx([gb * 1.2 * 2**i for i in range(3)] * 2, rel=1e-12)


def test_mips_are_moved_into_place(tmp_path):
    # a task counts as done once its MIP exists, so a MIP written in place could be seen, and
    # its job cancelled, while half written
    from petakit5d.testing.publish_spy import PublishSpy

    im = volume(IM_SIZE, "uint16", "moved")
    write_tiff(im, str(tmp_path / "im.tif"))
    write_zarr(im, str(tmp_path / "im.zarr"), block_size=(4, 16, 16))
    with PublishSpy() as spy:
        mip.save_mip_tiff(str(tmp_path / "im.tif"), str(tmp_path / "tif_MIP_z.tif"), axis=(1, 1, 1))
        mip.save_mip_zarr(
            str(tmp_path / "im.zarr"), str(tmp_path / "zarr_MIP_z.tif"), axis=(1, 1, 1)
        )
    for f in ("tif", "zarr"):
        verify_mips(str(tmp_path / f), im, "uint16")
        for lab in AXES:
            fn = tmp_path / f"{f}_MIP_{lab}.tif"
            assert spy.was_moved(fn), f"{fn} was written in place"
    assert len(list(tmp_path.glob("*.tif"))) == 7, "temporary files were left behind"


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_zarr_is_read_in_slabs_of_whole_chunks(tmp_path, cls, monkeypatch):
    # save_mip_zarr reads zarr files in z slabs of whole chunk layers with at least 64 chunks
    # (3 slabs of one layer here, 2 slabs of two layers in the box) and keeps a running max
    # along z. For single, lines of NaN and of -0 and 0, where max takes the first (0 first in
    # the second and third slab here).
    from petakit5d import io

    im = volume((12, 128, 128), cls, f"slabs/{cls}")
    if cls == "single":
        im[:, 2, :] = NAN32
        im[:, :, 4] = NEG0
        im[[1, 8], :, 4] = 0
        im[:, 6, 8] = NEG0
        im[[4, 11], 6, 8] = 0
    fn = str(tmp_path / "im.zarr")
    write_zarr(im, fn, block_size=(4, 16, 16))
    reads = []
    real = io.read_zarr

    def counting(path, *, input_bbox=None, **kw):
        reads.append(None if input_bbox is None else tuple(input_bbox))
        return real(path, input_bbox=input_bbox, **kw)

    monkeypatch.setattr(io, "read_zarr", counting)
    mip.save_mip_zarr(fn, str(tmp_path / "all_MIP_z.tif"), axis=(1, 1, 1))
    assert reads == [(0, 0, 0, 4, 128, 128), (4, 0, 0, 8, 128, 128), (8, 0, 0, 12, 128, 128)]
    mip.save_mip_zarr(fn, str(tmp_path / "z_MIP_z.tif"))
    b = (2, 4, 6, 10, 120, 110)  # MATLAB's [5, 7, 3, 120, 110, 10]
    reads.clear()
    mip.save_mip_zarr(fn, str(tmp_path / "bbox_MIP_z.tif"), axis=(1, 1, 1), input_bbox=b)
    assert reads == [(2, 4, 6, 8, 120, 110), (8, 4, 6, 10, 120, 110)]
    # a box from the last voxel of a chunk still touches 8 x 8 chunks a layer: one layer a slab
    reads.clear()
    mip.save_mip_zarr(fn, str(tmp_path / "edge_MIP_z.tif"), input_bbox=(0, 15, 15, 12, 128, 128))
    assert reads == [(z, 15, 15, z + 4, 128, 128) for z in (0, 4, 8)]
    crop = im[b[0] : b[3], b[1] : b[4], b[2] : b[5]]
    for lab, ax in AXES.items():
        same_bits(read_tiff(str(tmp_path / f"all_MIP_{lab}.tif")), tiff_image(ref_max(im, ax)))
        same_bits(read_tiff(str(tmp_path / f"bbox_MIP_{lab}.tif")), tiff_image(ref_max(crop, ax)))
    same_bits(read_tiff(str(tmp_path / "z_MIP_z.tif")), ref_max(im, -3)[0])
    assert not os.path.exists(tmp_path / "z_MIP_y.tif")


def test_zarr_mips_fall_back_to_blocks_when_a_read_fails(tmp_path, monkeypatch):
    # saveMIP_zarr_axis: if the slab reads fail, each MIP is made from blocks across its axis
    im = volume((12, 40, 36), "single", "fallback", True)
    fn = str(tmp_path / "im.zarr")
    write_zarr(im, fn, block_size=(4, 3, 2))

    def failing(*a, **k):
        raise MemoryError("out of memory")

    monkeypatch.setattr(mip, "_slab_mips", failing)
    mip.save_mip_zarr(fn, str(tmp_path / "f_MIP_z.tif"), axis=(1, 1, 1))
    for lab, ax in AXES.items():
        same_bits(read_tiff(str(tmp_path / f"f_MIP_{lab}.tif")), tiff_image(ref_max(im, ax)))
    with pytest.raises(mip.MipError, match="input BBox"):
        mip.save_mip_zarr(
            fn, str(tmp_path / "g_MIP_z.tif"), axis=(1, 1, 1), input_bbox=(0, 0, 0, 2, 2, 2)
        )


def test_save_mip_frame_writes_matlabs_unsqueezed_mips(tmp_path):
    im = volume((5, 7, 3), "single", "frame", True)
    fn = str(tmp_path / "sub" / "f_MIP_z.tif")
    mip.save_mip_frame(im, fn, axis=(1, 1, 1), dtype="single")
    for lab, ax in AXES.items():
        want = ref_max(im, ax)
        same_bits(read_tiff(str(tmp_path / "sub" / f"f_MIP_{lab}.tif")).reshape(want.shape), want)
    mip.save_mip_frame(im, fn, axis=(1, 1, 1))  # exists: skipped
    assert sorted(os.listdir(tmp_path / "sub")) == [f"f_MIP_{lab}.tif" for lab in "xyz"]


def test_project_3d_image_to_2d():
    im = volume((5, 7, 3), "uint16", "proj")
    assert np.array_equal(mip.project_3d_image_to_2d(im, "MIP_xy"), im.max(axis=0))
    assert np.array_equal(mip.project_3d_image_to_2d(im, "mip_xz"), im.max(axis=1))
    assert np.array_equal(mip.project_3d_image_to_2d(im, "mip_yz"), im.max(axis=2))
    assert np.array_equal(mip.project_3d_image_to_2d(im, "central_xy"), im[2])
    assert np.array_equal(mip.project_3d_image_to_2d(im, "central_xz"), im[:, 3])
    assert np.array_equal(mip.project_3d_image_to_2d(im, "central_yz"), im[:, :, 1])
    assert mip.project_3d_image_to_2d(im, "other") is im
    with pytest.raises(NotImplementedError):
        mip.project_3d_image_to_2d(im, "mean_xy")


# ---------------------------------------------------------------------------------------------
# BlockMIPTest in intent
# ---------------------------------------------------------------------------------------------


def make_image(cls):
    return volume(FLAG_SHAPE, cls, f"block/{cls}", special=cls == "single")


def block_boxes(shape, block):
    """the blocks of a grid in MATLAB's order (y fastest), reversed"""
    nz, ny, nx = shape
    boxes = []
    for iz in range(-(-nz // block[0])):
        for ix in range(-(-nx // block[2])):
            for iy in range(-(-ny // block[1])):
                s = (iz * block[0], iy * block[1], ix * block[2])
                e = (min(s[0] + block[0], nz), min(s[1] + block[1], ny), min(s[2] + block[2], nx))
                boxes.append((*s, *e))
    return boxes[::-1]


def write_py_flags(folder, im, block=FLAG_BLOCK, mbs=FLAG_MIP_BATCH, task_size=4):
    """the flags of tasks that wrote im in blocks, task_size blocks each, blocks in reverse
    order (BlockMIPTest's writeFlags)"""
    boxes = block_boxes(im.shape, block)
    cls = {np.dtype(np.uint16): "uint16", np.dtype(np.float32): "single"}[im.dtype]
    os.makedirs(folder, exist_ok=True)
    fns = []
    for s in range(0, len(boxes), task_size):
        bb = boxes[s : s + task_size]
        mips = [mip.block_mips(im[b[0] : b[3], b[1] : b[4], b[2] : b[5]], b, mbs, cls) for b in bb]
        fn = f"{folder}/blocks_{s + 1}_{s + len(bb)}.mat"
        mip.save_block_mips(fn, bb, mips, mbs, 1)
        fns.append(fn)
    return fns


def range_ref(im, axis, starts):
    ends = [*starts[1:], im.shape[axis]]
    parts = [
        ref_max(np.take(im, range(s, e), axis=axis), axis)
        for s, e in zip(starts, ends, strict=True)
    ]
    return np.concatenate(parts, axis=axis)


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_merged_slabs_match_the_image(tmp_path, cls):
    # blocks not aligned with the MIP batches, written by tasks in an order that isn't their
    # order in the image; for single, a region of NaN and of -0 and 0
    im = make_image(cls)
    fns = write_py_flags(str(tmp_path / "f"), im)
    for p in range(3):
        slab = mip.merge_block_mips(fns, p, im.shape, FLAG_MIP_BATCH, cls)
        same_bits(slab, range_ref(im, p, list(range(0, im.shape[p], FLAG_MIP_BATCH[p]))))


def test_block_mips_without_the_kernel_match(monkeypatch):
    # classes (MATLAB: platforms) without range_max_3d take the max of each part
    im = volume((13, 20, 15), "uint16", "nomex") // 2  # (fits int16)
    box = (3, 2, 1, 16, 22, 16)
    expected = mip.block_mips(im, box, (5, 8, 6), "uint16")
    monkeypatch.setattr(mip, "_RANGE_DTYPES", frozenset())
    got = mip.block_mips(im, box, (5, 8, 6), "uint16")
    for p in range(3):
        same_bits(got[p], expected[p])
    other = mip.block_mips(im.astype(np.int16), box, (5, 8, 6), "uint16")
    for p in range(3):
        same_bits(other[p], expected[p])


def test_block_mips_are_cast_as_the_zarr_stores_them(tmp_path):
    # single blocks written to a uint16 zarr (large-file decon with save16bit) are rounded and
    # saturated as written
    rng = CounterRNG(0, "t02/castblock")
    block = (rng.normal((9, 12, 10)) * 40000 + 30000).astype(np.float32)
    block[2, 0:5, 1] = [0.5, 1.5, -3, NAN32, 70000]
    block[:, :, 3] = NAN32
    fn = str(tmp_path / "out.zarr")
    create_zarr(fn, data_size=(9, 12, 10), block_size=(4, 4, 4), dtype="uint16")
    write_zarr(block, fn, bbox=(0, 0, 0, 9, 12, 10), create=False)
    stored = read_zarr(fn)
    mips = mip.block_mips(block, (0, 0, 0, 9, 12, 10), (4, 5, 4), "uint16")
    for p, bs in enumerate((4, 5, 4)):
        same_bits(mips[p], range_ref(stored, p, list(range(0, stored.shape[p], bs))))


def test_flag_is_moved_into_place(tmp_path):
    # the flag marks the task as done, so it must not be seen half written
    from petakit5d.testing.publish_spy import PublishSpy

    im = volume((4, 6, 5), "uint16", "flagmove")
    mips = mip.block_mips(im, (0, 0, 0, 4, 6, 5), (8, 8, 8), "uint16")
    fn = str(tmp_path / "blocks_1_1.mat")
    with PublishSpy() as spy:
        mip.save_block_mips(fn, [(0, 0, 0, 4, 6, 5)], [mips], (8, 8, 8), 1.5)
    assert spy.was_moved(fn), "the flag was written in place"
    assert [p.name for p in tmp_path.glob("*.mat")] == ["blocks_1_1.mat"]
    a = read_mat(fn)
    assert float(a["t"]) == 1.5
    z = a["zMIPs"][0] if isinstance(a["zMIPs"], list) else a["zMIPs"].ravel()[0]
    same_bits(from_matlab(np.asarray(z)), ref_max(im, -3)[0])
    assert np.array_equal(np.asarray(a["blockBboxes"]).ravel(), [1, 1, 1, 6, 5, 4])
    assert np.asarray(a["mipBatchSize"]).dtype == np.float64


def test_merge_needs_mips_of_all_blocks(tmp_path):
    # flags of tasks without MIPs, MIPs for other batches, or blocks that don't cover the image
    # give no slab, so that the caller reads the image
    im = make_image("uint16")
    sz = im.shape
    fns = write_py_flags(str(tmp_path / "f"), im)

    def merge(f=fns, bs=FLAG_MIP_BATCH, cls="uint16"):
        return mip.merge_block_mips(f, 1, sz, bs, cls)

    assert merge() is not None
    assert merge(bs=(5, 16, 6)) is None, "other batch size"
    # MIPs saved for other batches along z only: their y slab would fit
    a = read_mat(fns[0])
    a["mipBatchSize"] = np.array([[8.0, 6.0, 6.0]])
    write_mat(fns[0], a, version="6")
    assert merge() is None, "flag for other batches"
    a["mipBatchSize"] = np.array([[8.0, 6.0, 5.0]])
    write_mat(fns[0], a, version="6")
    # the last flag without the MIP of its last block
    a = read_mat(fns[-1])
    y = a["yMIPs"]
    items = list(y) if isinstance(y, list) else list(np.asarray(y, dtype=object).ravel())
    cut = np.empty((len(items) - 1, 1), dtype=object)
    for k, m in enumerate(items[:-1]):
        cut[k, 0] = m
    a["yMIPs"] = cut
    write_mat(fns[-1], a, version="6")
    assert merge() is None, "flag with a MIP missing"
    a["yMIPs"] = y
    write_mat(fns[-1], a, version="6")
    assert merge() is not None
    assert merge(cls="single") is None, "other data type"
    assert merge(f=fns[1:]) is None, "missing flag"
    assert merge(f=[*fns, fns[0]]) is None, "repeated flag"
    write_mat(fns[1], {"t": 1.0}, version="6")
    assert merge() is None, "flag with only t"
    Path(fns[1]).write_bytes(b"")
    assert merge() is None, "empty flag"


def test_merge_needs_blocks_that_cover_the_image_once(tmp_path):
    # two overlapping blocks and a gap of the same size add up to the image size
    im = volume((6, 10, 8), "uint16", "cover")
    boxes = [(0, 0, 0, 6, 5, 8), (0, 0, 0, 6, 5, 8)]
    mips = [mip.block_mips(im[:, :5], b, (4, 4, 4), "uint16") for b in boxes]
    fn = str(tmp_path / "blocks_1_2.mat")
    mip.save_block_mips(fn, boxes, mips, (4, 4, 4))
    for p in range(3):
        assert mip.merge_block_mips([fn], p, im.shape, (4, 4, 4), "uint16") is None, p


def test_merge_keeps_to_the_memory_limit(tmp_path):
    # the y slab and its coverage count (5 x 29 x 23 x 3 bytes) and the blocks' MIPs (blocks
    # overlap 10 batches along y in all: 10 x 29 x 23 x 2 bytes) are 23345 bytes, 2.174e-5 GiB
    im = make_image("uint16")
    fns = write_py_flags(str(tmp_path / "f"), im)
    args = (fns, -2, im.shape, FLAG_MIP_BATCH, "uint16")
    assert mip.merge_block_mips(*args, 2.1e-5) is None, "over the limit"
    assert mip.merge_block_mips(*args, 2.2e-5) is not None, "under the limit"
    assert mip.merge_block_mips(*args, NAN) is not None, "unknown memory"


def test_batch_size_is_rounded_to_blocks_and_capped():
    bs = mip.mip_zarr_batch_size
    assert bs((900, 5000, 3000), (256, 256, 256), "uint16") == (900, 2048, 2048)
    assert bs((900, 5000, 3000), (100, 256, 300), "uint16") == (900, 2048, 2100)
    want = (256, 256, 256)
    assert bs((900, 5000, 3000), (256, 256, 256), "uint16", batch_size=(100, 100, 100)) == want


def test_batch_size_grows_for_slabs_over_100gb():
    # the batch grows by a block along the axis of the largest slab until no slab is over 100
    # GiB (3 steps here), or the batch is (21 steps, to a 104 GiB batch)
    bs = mip.mip_zarr_batch_size
    assert bs((11000, 100000, 100000), (256, 256, 256), "uint16") == (2304, 2304, 2304)
    assert bs((10000, 200000, 100000), (256, 256, 256), "uint16") == (4352, 3584, 3584)
    want = (2048, 2048, 2048)
    assert bs((10000, 200000, 100000), (256, 256, 256), "uint16", mip_slab=True) == want


def zarr_copies(folder, cls, native=True):
    """the same zarr in two folders (BlockMIPTest's writeZarrCopies)"""
    im = make_image(cls)
    fns = []
    for i in range(2):
        sub = Path(folder) / f"c{i}"
        sub.mkdir(parents=True)
        fn = str(sub / "im.zarr")
        if native:
            write_zarr(im, fn, block_size=(5, 4, 3))
        else:
            synth.write_zarr(fn, im, block_size=(5, 4, 3), layout="matlab")
        fns.append(fn)
    return fns, im


def verify_same_mips(fn1, fn2):
    p1 = os.path.dirname(fn1) + "/MIPs/"
    p2 = os.path.dirname(fn2) + "/MIPs/"
    for lab in AXES:
        s1, s2 = f"{p1}im_{lab}.zarr", f"{p2}im_{lab}.zarr"
        assert os.path.isdir(s2), f"no {lab} slab"
        z1 = json.loads(Path(s1, ".zarray").read_text())
        assert json.loads(Path(s2, ".zarray").read_text()) == z1, lab
        same_bits(read_zarr(s2), read_zarr(s1))
        same_bits(read_tiff(f"{p2}im_MIP_{lab}.tif"), read_tiff(f"{p1}im_MIP_{lab}.tif"))
    assert os.listdir(p2 + "zarr_flag") == [], "the MIP flag folder was left behind"


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_xr_mip_zarr_merges_the_block_mips(tmp_path, cls, monkeypatch):
    # the same slab zarrs (data, chunks, data type) and MIPs as from the zarr, without job
    # framework tasks
    from petakit5d.jobs import framework

    fns, im = zarr_copies(tmp_path, cls)
    flags = write_py_flags(str(tmp_path / "f"), im)
    mip.mip_zarr(fns[0], batch_size=FLAG_MIP_BATCH, parse_cluster=False)
    calls = []
    monkeypatch.setattr(framework, "run_tasks", lambda *a, **k: calls.append(1))
    mip.mip_zarr(fns[1], batch_size=FLAG_MIP_BATCH, parse_cluster=False, block_flag_fullpaths=flags)
    assert calls == [], "the MIPs were computed by reading the zarr"
    verify_same_mips(fns[0], fns[1])
    for lab, ax in AXES.items():
        got = read_tiff(os.path.dirname(fns[0]) + f"/MIPs/im_MIP_{lab}.tif")
        same_bits(got, tiff_image(ref_max(im, ax)))


def test_xr_mip_zarr_reads_the_zarr_if_a_flag_has_no_mips(tmp_path):
    # all or nothing: the other flags' MIPs are wrong here, so a partial merge would show
    fns, im = zarr_copies(tmp_path, "uint16")
    flags = write_py_flags(str(tmp_path / "f"), im + 1)
    write_mat(flags[-1], {"t": 1.0}, version="6")
    mip.mip_zarr(fns[0], batch_size=FLAG_MIP_BATCH, parse_cluster=False)
    mip.mip_zarr(fns[1], batch_size=FLAG_MIP_BATCH, parse_cluster=False, block_flag_fullpaths=flags)
    verify_same_mips(fns[0], fns[1])


def test_xr_mip_zarr_reads_the_zarr_without_memory_for_the_merge(tmp_path, monkeypatch):
    # the flags' MIPs are wrong here, so a merge would show
    from petakit5d.core import memory

    fns, im = zarr_copies(tmp_path, "uint16")
    flags = write_py_flags(str(tmp_path / "f"), im + 1)
    mip.mip_zarr(fns[0], batch_size=FLAG_MIP_BATCH, parse_cluster=False)
    monkeypatch.setattr(memory, "get_system_memory", lambda: 1e-6)
    mip.mip_zarr(fns[1], batch_size=FLAG_MIP_BATCH, parse_cluster=False, block_flag_fullpaths=flags)
    verify_same_mips(fns[0], fns[1])


def test_xr_mip_zarr_resumes_from_existing_slabs(tmp_path):
    # slab zarrs from an earlier run are used for missing MIPs (the flags' MIPs are wrong here)
    fns, im = zarr_copies(tmp_path, "uint16")
    flags = write_py_flags(str(tmp_path / "f"), im + 1)
    mip.mip_zarr(fns[0], batch_size=FLAG_MIP_BATCH, parse_cluster=False)
    mip_fn = os.path.dirname(fns[0]) + "/MIPs/im_MIP_x.tif"
    expected = read_tiff(mip_fn)
    os.remove(mip_fn)
    mip.mip_zarr(fns[0], batch_size=FLAG_MIP_BATCH, parse_cluster=False, block_flag_fullpaths=flags)
    same_bits(read_tiff(mip_fn), expected)


# ---------------------------------------------------------------------------------------------
# the port's own checks
# ---------------------------------------------------------------------------------------------


def test_large_file_path_publishes_slabs_and_mips_before_its_flag(tmp_path):
    # the wrapper counts a large file as done once its last flagged MIP exists: the slabs and
    # the other MIPs are in place before that MIP is moved there
    from petakit5d.testing.publish_spy import PublishSpy

    path, ims, names = write_dataset(tmp_path / "d", "zarr")
    with PublishSpy() as spy:
        mip.mip_dataset(
            path,
            axis=(1, 1, 1),
            channel_patterns=CHANNELS[:1],
            zarr_file=True,
            large_file=True,
            batch_size=(4, 16, 16),
            parse_cluster=False,
        )
    flag = f"{path}MIPs/{names[0]}_MIP_z.tif"
    seen = spy.files_when_moved(flag)
    for lab in "yx":
        assert f"{names[0]}_MIP_{lab}.tif" in seen, lab
    for lab in "yxz":
        assert f"{names[0]}_{lab}.zarr" in seen, lab
    mips = os.path.normpath(f"{path}MIPs")
    order = [c.dst for c in spy.calls if c.ok and os.path.normpath(c.dst).startswith(mips)]
    slabs = [i for i, d in enumerate(order) if d.endswith(".zarr")]
    tifs = [i for i, d in enumerate(order) if d.endswith("tif")]
    assert slabs and tifs and max(slabs) < min(tifs), order


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_native_and_matlab_layout_inputs_give_the_same_files(tmp_path, cls):
    # existing data in MATLAB's layout: the same MIPs and slabs as from a native array
    fn_n, _ = zarr_copies(tmp_path / "n", cls)
    fn_m, _ = zarr_copies(tmp_path / "m", cls, native=False)
    for fns in (fn_n, fn_m):
        mip.mip_zarr(fns[0], batch_size=(5, 8, 6), parse_cluster=False)
        out = os.path.dirname(fns[1]) + "/MIPs/s_MIP_z.tif"
        os.makedirs(os.path.dirname(out))
        mip.save_mip_zarr(fns[1], out, axis=(1, 1, 1))
    verify_same_mips(fn_n[0], fn_m[0])
    for lab in AXES:
        a = read_tiff(os.path.dirname(fn_n[1]) + f"/MIPs/s_MIP_{lab}.tif")
        same_bits(a, read_tiff(os.path.dirname(fn_m[1]) + f"/MIPs/s_MIP_{lab}.tif"))


@pytest.mark.parametrize("threads", [1, 3])
def test_results_do_not_depend_on_threads(tmp_path, threads, monkeypatch):
    from petakit5d.kernels import mip as kmip

    monkeypatch.setattr(kmip, "_available", threads)
    im = make_image("single")
    fns = write_py_flags(str(tmp_path / "f"), im)
    for p in range(3):
        same_bits(
            mip.merge_block_mips(fns, p, im.shape, FLAG_MIP_BATCH, "single"),
            range_ref(im, p, list(range(0, im.shape[p], FLAG_MIP_BATCH[p]))),
        )
        same_bits(mip.axis_max(im, p - 3), ref_max(im, p - 3))


def test_merge_memory_counts_the_slab_and_the_blocks_mips():
    # the slabs of a 18000 x 40000 x 30000 uint16 image whose block tasks wrote columns 1024
    # deep along z
    shape, region = (18000, 40000, 30000), (1024, 40000, 30000)
    mem = mip.merge_memory(shape, region, (256, 256, 256), "uint16")
    mbs = mip.mip_zarr_batch_size(shape, (256, 256, 256), "uint16")
    want = 0.0
    total = float(18000 * 40000 * 30000)
    for n, r, b in zip(shape, region, mbs, strict=True):
        n_mip, n_reg = math.ceil(n / b), math.ceil(n / r)
        plane = total / n
        want = max(want, (plane * n_mip * 3 + plane * (n_reg + n_mip - 1) * 2) / 2**30)
    assert mem == pytest.approx(want, rel=1e-15)


def test_mip_block_skips_done_batches_and_checks_its_folders(tmp_path, capsys):
    im = volume((6, 10, 8), "uint16", "blk")
    fn = str(tmp_path / "im.zarr")
    write_zarr(im, fn, block_size=(3, 5, 4))
    box = [(0, 0, 0, 6, 10, 8)]
    slabs = ["a", "b", "c"]
    with pytest.raises(mip.MipError, match="block directory"):
        mip.mip_block([0], fn, slabs, str(tmp_path / "no" / "f.mat"), box, (0, 0, 0), (6, 10, 8))
    flag = tmp_path / "blocks_1_1.mat"
    flag.write_bytes(b"")
    assert mip.mip_block([0], fn, slabs, str(flag), box, (0, 0, 0), (6, 10, 8)) is None
    assert "already exist" in capsys.readouterr().out
    with pytest.raises(mip.MipError, match="doesnot exist"):
        mip.mip_block(
            [0],
            str(tmp_path / "x.zarr"),
            slabs,
            str(flag),
            box,
            (0, 0, 0),
            (6, 10, 8),
            overwrite=True,
        )


def test_resampling_before_pooling_is_not_ported(tmp_path):
    im = volume((6, 10, 8), "uint16", "rs")
    fn = str(tmp_path / "im.zarr")
    write_zarr(im, fn, block_size=(3, 5, 4))
    with pytest.raises(NotImplementedError, match="imresize3"):
        mip.mip_block(
            [0],
            fn,
            ["a", "b", "c"],
            str(tmp_path / "f.mat"),
            [(0, 0, 0, 6, 10, 8)],
            (0, 0, 0),
            (2, 2, 2, 1, 1, 1, 2, 2, 2),
        )


def test_matlab_layer_converts_axis_and_box(tmp_path):
    # the MATLAB-name functions take [y x z] flags and MATLAB's 1-based boxes
    from petakit5d import matlab

    im = volume(IM_SIZE, "uint16", "layer")
    write_tiff(im, str(tmp_path / "im.tif"))
    matlab.saveMIP_tiff(
        str(tmp_path / "im.tif"),
        str(tmp_path / "a_MIP_z.tif"),
        "axis",
        [1, 0, 0],
        "inputBbox",
        [5, 7, 3, 30, 33, 10],
    )
    crop = im[2:10, 4:30, 6:33]
    same_bits(read_tiff(str(tmp_path / "a_MIP_y.tif")), tiff_image(ref_max(crop, -2)))
    assert not os.path.exists(tmp_path / "a_MIP_z.tif")


@pytest.mark.parametrize(
    "box",
    [(1, 1, 2, 20, 30, 27), (0, 0, 0, 23, 30, 29), (0, 0, 0, 23, 37, 25), (3, 0, 0, 19, 37, 29)],
)
def test_xr_mip_zarr_with_a_box_ending_before_the_image_projects_the_box(tmp_path, box):
    # the last batches end at the box's end (MATLAB since 11c7d0b; before, at the image's end,
    # past the slab: zarr:inputError): slabs of the box's size and the MIPs of the box
    fns, im = zarr_copies(tmp_path, "uint16")
    mip.mip_zarr(fns[0], batch_size=(5, 8, 6), input_bbox=box, parse_cluster=False)
    crop = im[box[0] : box[3], box[1] : box[4], box[2] : box[5]]
    out = os.path.dirname(fns[0]) + "/MIPs/"
    for lab, ax in AXES.items():
        same_bits(read_tiff(f"{out}im_MIP_{lab}.tif"), tiff_image(ref_max(crop, ax)))
        slab = read_zarr(f"{out}im_{lab}.zarr")
        n = crop.shape[ax]
        k = {"z": 5, "y": 8, "x": 6}[lab]
        assert slab.shape[ax] == -(-n // k), lab
        same_bits(slab, range_ref(crop, ax % 3, list(range(0, n, k))))


def test_xr_mip_zarr_with_a_box_reads_the_zarr_not_the_flags(tmp_path):
    # the blocks' MIPs are of the whole image: with a box the slabs come from the zarr (the
    # flags' MIPs are wrong here, so a merge would show)
    fns, im = zarr_copies(tmp_path, "uint16")
    flags = write_py_flags(str(tmp_path / "f"), im + 1)
    box = (3, 5, 4, FLAG_SHAPE[0], FLAG_SHAPE[1], FLAG_SHAPE[2])
    mip.mip_zarr(fns[0], batch_size=FLAG_MIP_BATCH, input_bbox=box, parse_cluster=False)
    mip.mip_zarr(
        fns[1],
        batch_size=FLAG_MIP_BATCH,
        input_bbox=box,
        parse_cluster=False,
        block_flag_fullpaths=flags,
    )
    verify_same_mips(fns[0], fns[1])
    crop = im[box[0] :, box[1] :, box[2] :]
    got = read_tiff(os.path.dirname(fns[1]) + "/MIPs/im_MIP_z.tif")
    same_bits(got, tiff_image(ref_max(crop, -3)))
