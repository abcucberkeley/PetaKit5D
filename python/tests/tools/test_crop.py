"""Cropping (petakit5d.tools.crop, task T03).

Oracle cases against MATLAB: XR_crop_dataset (tiff, zarr and large zarr datasets, several
channels and datasets, boxes touching and leaving the image, padding, single slices, the
'center' and 'moving' crop types, existing results, tiff saved as zarr), XR_crop_frame (pad
values, MATLAB's errors), XR_crop_zarr (the start before the image without pad, the box past
the end, no MIPs, sizes not divisible by the batch and chunk sizes) and trimBorder, for uint8,
uint16 and single (NaN of two payloads, -0 and 0). Output tiffs are compared byte for byte, zarr
results and MIP slabs by data and metadata after the layout mapping, parameters.mat by
content. Inputs are MATLAB's: tiffs and zarr arrays in MATLAB's layout, read by both sides.

Then CropDatasetTest (tests/crop) in intent on the port's (z, y, x) data and native zarr
arrays, CropMexTest.trimBorderModes, the publish order and the writer records of the large-file
path (spies, a fake Slurm), the job shapes (a framework spy), and results that don't depend on
the zarr layout or on threads.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.io import create_zarr, read_mat, read_tiff, read_zarr, write_tiff, write_zarr
from petakit5d.testing import oracle, slurm_fake, synth
from petakit5d.testing.cases import Case, File, In, InPath, PerAxis, Tree, Value, Vol, Work
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import crop

# ---------------------------------------------------------------------------------------------
# data: (z, y, x) volumes made from counter-based random numbers (bit-stable everywhere)
# ---------------------------------------------------------------------------------------------

NAN32 = np.float32(NAN)
# a second NaN payload (MATLAB: typecast(uint32(hex2dec('7FC00001')), 'single'))
NAN32_B = np.array([0x7FC00001], dtype=np.uint32).view(np.float32)[0]


def volume(shape, cls, name, special=False):
    """a deterministic (z, y, x) volume: integers in [0, max), single normal * 100; with
    special (single): NaN of two payloads, -0 and 0"""
    rng = CounterRNG(0, f"t03/{name}")
    shape = tuple(int(s) for s in shape)
    if cls == "uint8":
        return rng.integers(0, 255, shape, dtype=np.uint8)
    if cls == "uint16":
        return rng.integers(0, 60000, shape, dtype=np.uint16)
    a = np.ascontiguousarray((rng.normal(shape) * 100.0).astype(np.float32))
    if special:
        u = rng.uniform((3, *shape))
        a[u[0] < 0.05] = NAN32
        a[u[1] < 0.05] = np.float32(-0.0)
        a[u[2] < 0.05] = NAN32_B
    return a


def same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype, (a.dtype, b.dtype)
    assert a.shape == b.shape, (a.shape, b.shape)
    if a.dtype.kind == "f":
        assert a.view(f"u{a.itemsize}").tobytes() == b.view(f"u{b.itemsize}").tobytes()
    else:
        assert np.array_equal(a, b)


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

IM_SIZE = (12, 40, 36)  # CropDatasetTest's ImSize [40 36 12]
BBOX = [5, 7, 3, 30, 33, 10]  # MATLAB's [ymin xmin zmin ymax xmax zmax]
DATASET = ("scan_CamA_ch0_t0000", "scan_CamB_ch0_t0000", "scan_CamA_ch1_t0000")
CHANS = ["channelPatterns", ["CamA_ch0", "CamB_ch0"]]
ONE = ["channelPatterns", ["CamA_ch0"]]
LOCAL = ["parseCluster", False, "masterCompute", True]
SMALL = ["batchSize", [16, 16, 4], "blockSize", [8, 8, 4]]
ERRORS = (crop.CropError, ValueError, OSError, RuntimeError, IndexError)


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def py_dataset(*args):
    return call_matlab("XR_crop_dataset", *args)


def py_frame(*args):
    return call_matlab("XR_crop_frame", *args)


def py_zarr(d, *args):
    return call_matlab("XR_crop_zarr", d + "im.zarr", *args)


def dataset_tree(fmt, cls="uint16", shape=IM_SIZE, chunks=(4, 16, 16), special=False, tag=""):
    """CropDatasetTest's writeDataset: three files, the third matching no channel pattern of
    the tests (zarr arrays in MATLAB's layout)"""

    def make(folder):
        for k, n in enumerate(DATASET):
            a = volume(shape, cls, f"ds/{fmt}/{cls}/{tag}{k}", special)
            if fmt == "zarr":
                synth.write_zarr(Path(folder) / f"{n}.zarr", a, block_size=chunks, layout="matlab")
            else:
                synth.write_tiff(Path(folder) / f"{n}.tif", a)

    return make


def time_points_tree(nt, tag):
    """movingCrop's dataset: nt time points of one channel"""

    def make(folder):
        for t in range(nt):
            a = volume(IM_SIZE, "uint16", f"tp/{tag}/{t}")
            synth.write_tiff(Path(folder) / f"scan_CamA_ch0_t{t:04d}.tif", a)

    return make


def existing_tree():
    """existingOutputIsNotOverwritten: a result of the first file is there"""

    def make(folder):
        dataset_tree("tif", tag="ex")(folder)
        (Path(folder) / "Cropped").mkdir()
        sentinel = np.arange(16, dtype=np.uint16).reshape(1, 4, 4)
        synth.write_tiff(Path(folder) / "Cropped" / f"{DATASET[0]}.tif", sentinel)

    return make


def zarr_tree(shape, cls, name, chunks, special=False):
    """a folder with one zarr array im.zarr in MATLAB's layout (chunks (z, y, x))"""

    def make(folder):
        a = volume(shape, cls, name, special)
        synth.write_zarr(Path(folder) / "im.zarr", a, block_size=chunks, layout="matlab")

    return make


def tif_file(shape, cls, name, special=False):
    def make(path):
        synth.write_tiff(path, volume(shape, cls, name, special))

    return make


def dataset_case(cid, inputs, args, error=False):
    return Case(
        cid,
        code="XR_crop_dataset(in{:}); out = {};",
        nargout=0,
        args=args,
        python=py_dataset,
        inputs=inputs,
        expect_error=error,
        py_error=ERRORS,
    )


def ds(cid, make, *args, error=False):
    return dataset_case(cid, {"d": Tree(make, "d")}, [InPath("d"), *args, *LOCAL], error)


def frame_case(cid, make, name, out, *args, error=False):
    return Case(
        cid,
        matlab="XR_crop_frame",
        nargout=0,
        args=[InPath("im"), Work(out), *args, *LOCAL],
        python=py_frame,
        inputs={"im": File(make, name)},
        expect_error=error,
        py_error=ERRORS,
    )


def zarr_case(cid, make, *args, error=False):
    return Case(
        cid,
        code="XR_crop_zarr([in{1}, 'im.zarr'], in{2 : end}); out = {};",
        nargout=0,
        args=[InPath("d"), *args, "parseCluster", False],
        python=py_zarr,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        # (the writer record of a crop that failed: its owner is the job and the process)
        tree=tree(masks=[(r"owner=[^\n]*", "owner=<owner>")]),
        py_error=ERRORS,
    )


def _dataset_cases():
    pad_box = [-2, -1, 0, 45, 40, 15]
    cases = [
        ds("ds.tif.uint16.fixed", dataset_tree("tif"), BBOX, *CHANS),
        ds("ds.tif.uint8.fixed", dataset_tree("tif", "uint8"), BBOX, *ONE),
        ds("ds.tif.single.fixed", dataset_tree("tif", "single", special=True), BBOX, *ONE),
        ds(
            "ds.zarr.uint16.zarr_out",
            dataset_tree("zarr"),
            BBOX,
            *CHANS,
            "zarrFile",
            True,
            "saveZarr",
            True,
        ),
        ds(
            "ds.zarr.uint16.tif_out",
            dataset_tree("zarr"),
            BBOX,
            *CHANS,
            "zarrFile",
            True,
            "saveZarr",
            False,
        ),
        ds(
            "ds.zarr.single.zarr_out_block",
            dataset_tree("zarr", "single", special=True),
            BBOX,
            *ONE,
            "zarrFile",
            True,
            "saveZarr",
            True,
            "blockSize",
            [8, 16, 4],
        ),
        ds(
            "ds.large.uint16",
            dataset_tree("zarr"),
            BBOX,
            *CHANS,
            "zarrFile",
            True,
            "largeFile",
            True,
            *SMALL,
        ),
        ds(
            "ds.large.uint8",
            dataset_tree("zarr", "uint8"),
            BBOX,
            *ONE,
            "zarrFile",
            True,
            "largeFile",
            True,
            *SMALL,
        ),
        ds(
            "ds.large.single",
            dataset_tree("zarr", "single", special=True),
            BBOX,
            *ONE,
            "zarrFile",
            True,
            "largeFile",
            True,
            *SMALL,
        ),
        ds(
            "ds.large.uint16.uneven",
            dataset_tree("zarr", shape=(23, 37, 29), chunks=(5, 7, 6), tag="un"),
            [2, 3, 2, 37, 29, 23],
            *ONE,
            "zarrFile",
            True,
            "largeFile",
            True,
            "batchSize",
            [16, 12, 10],
            "blockSize",
            [8, 6, 5],
        ),
        ds(
            "ds.large.uint16.whole_image_one_batch",
            dataset_tree("zarr"),
            [1, 1, 1, 40, 36, 12],
            *ONE,
            "zarrFile",
            True,
            "largeFile",
            True,
        ),
        ds("ds.tif.uint16.pad", dataset_tree("tif"), pad_box, *ONE, "pad", True),
        ds(
            "ds.tif.single.pad",
            dataset_tree("tif", "single", special=True),
            pad_box,
            *ONE,
            "pad",
            True,
        ),
        ds(
            "ds.zarr.uint16.pad",
            dataset_tree("zarr"),
            pad_box,
            *ONE,
            "pad",
            True,
            "zarrFile",
            True,
            "saveZarr",
            True,
        ),
        ds(
            "ds.large.uint16.pad",
            dataset_tree("zarr"),
            pad_box,
            *ONE,
            "pad",
            True,
            "zarrFile",
            True,
            "largeFile",
            True,
            *SMALL,
        ),
        ds(
            "ds.large.single.pad",
            dataset_tree("zarr", "single", special=True),
            pad_box,
            *ONE,
            "pad",
            True,
            "zarrFile",
            True,
            "largeFile",
            True,
            *SMALL,
        ),
        ds(
            "ds.tif.uint16.pad_end_only",
            dataset_tree("tif"),
            [30, 20, 8, 60, 50, 20],
            *ONE,
            "pad",
            True,
        ),
        ds("ds.tif.uint16.clip", dataset_tree("tif"), [30, 20, 8, 60, 50, 20], *ONE),
        ds(
            "ds.zarr.uint16.clip",
            dataset_tree("zarr"),
            [30, 20, 8, 60, 50, 20],
            *ONE,
            "zarrFile",
            True,
            "saveZarr",
            True,
        ),
        ds("ds.tif.uint16.touching_borders", dataset_tree("tif"), [1, 1, 1, 40, 36, 12], *ONE),
        ds(
            "ds.zarr.uint16.single_slice",
            dataset_tree("zarr"),
            [5, 7, 6, 30, 33, 6],
            *ONE,
            "zarrFile",
            True,
            "saveZarr",
            True,
            "resultDirName",
            "Cropped_small",
        ),
        ds(
            "ds.large.uint16.single_slice",
            dataset_tree("zarr"),
            [5, 7, 6, 30, 33, 6],
            *ONE,
            "zarrFile",
            True,
            "largeFile",
            True,
            *SMALL,
            "resultDirName",
            "Cropped_large",
        ),
        ds("ds.tif.uint16.single_slice", dataset_tree("tif"), [5, 7, 6, 30, 33, 6], *ONE),
        ds(
            "ds.tif.uint16.center",
            dataset_tree("tif"),
            [1, 1, 1, 10, 9, 4],
            *ONE,
            "cropType",
            "center",
        ),
        ds(
            "ds.tif.uint16.center_odd",
            dataset_tree("tif"),
            [3, 3, 2, 13, 11, 6],
            *ONE,
            "cropType",
            "center",
        ),
        ds(
            "ds.tif.uint16.moving",
            time_points_tree(3, "mv"),
            [3, 3, 2, 12, 12, 6],
            *ONE,
            "cropType",
            "moving",
            "lastStartCoords",
            [21, 13, 6],
        ),
        ds(
            "ds.tif.uint16.moving_halves",
            time_points_tree(4, "mh"),
            [3, 3, 2, 12, 12, 6],
            *ONE,
            "cropType",
            "moving",
            "lastStartCoords",
            [6, 4, 3],
        ),
        ds(
            "ds.tif.uint16.moving_out_of_image",
            time_points_tree(2, "mo"),
            [3, 3, 2, 12, 12, 6],
            *ONE,
            "cropType",
            "moving",
            "lastStartCoords",
            [32, 13, 6],
            error=True,
        ),
        ds("ds.tif.uint16.existing", existing_tree(), BBOX, *CHANS),
        ds("ds.tif.uint16.zarr_out", dataset_tree("tif"), BBOX, *ONE, "saveZarr", True),
        ds("ds.tif.uint16.result_dir", dataset_tree("tif"), BBOX, *ONE, "resultDirName", "c2"),
    ]
    cases.append(
        dataset_case(
            "ds.two_datasets.tif",
            {
                "d1": Tree(dataset_tree("tif", tag="a"), "d1"),
                "d2": Tree(dataset_tree("tif", tag="b"), "d2"),
            },
            [[InPath("d1"), InPath("d2")], BBOX, *CHANS, *LOCAL],
        )
    )
    cases.append(
        dataset_case(
            "ds.two_datasets.large",
            {
                "d1": Tree(dataset_tree("zarr", tag="a"), "d1"),
                "d2": Tree(dataset_tree("zarr", tag="b"), "d2"),
            },
            [
                [InPath("d1"), InPath("d2")],
                BBOX,
                *ONE,
                "zarrFile",
                True,
                "largeFile",
                True,
                *SMALL,
                *LOCAL,
            ],
        )
    )
    # (MATLAB's ch_tps(f) of the last dataset, numbered over all files: an index error)
    cases.append(
        dataset_case(
            "ds.two_datasets.moving_quirk",
            {
                "d1": Tree(time_points_tree(2, "q1"), "d1"),
                "d2": Tree(time_points_tree(2, "q2"), "d2"),
            },
            [
                [InPath("d1"), InPath("d2")],
                [3, 3, 2, 12, 12, 6],
                *ONE,
                "cropType",
                "moving",
                "lastStartCoords",
                [6, 4, 3],
                *LOCAL,
            ],
            error=True,
        )
    )
    return cases


def _frame_cases():
    s = (11, 17, 13)  # MATLAB 17 x 13 x 11
    u16 = tif_file(s, "uint16", "f1")
    sgl = tif_file(s, "single", "f2", True)
    pad = [-1, -2, 0, 20, 15, 13]
    return [
        frame_case(
            "frame.tif.uint16.pad_value", u16, "im.tif", "o.tif", pad, "pad", True, "padValue", 7.5
        ),
        frame_case(
            "frame.tif.uint16.pad_value_saturate",
            u16,
            "im.tif",
            "o.tif",
            pad,
            "pad",
            True,
            "padValue",
            70000,
        ),
        frame_case(
            "frame.tif.uint16.pad_value_negative",
            u16,
            "im.tif",
            "o.tif",
            pad,
            "pad",
            True,
            "padValue",
            -3,
        ),
        frame_case(
            "frame.tif.single.pad_value",
            sgl,
            "im.tif",
            "o.tif",
            pad,
            "pad",
            True,
            "padValue",
            -3.25,
        ),
        frame_case(
            "frame.tif.single.pad_nan",
            sgl,
            "im.tif",
            "o.tif",
            pad,
            "pad",
            True,
            "padValue",
            float("nan"),
        ),
        frame_case(
            "frame.tif.uint16.zarr_out_pad",
            u16,
            "im.tif",
            "o.zarr",
            pad,
            "pad",
            True,
            "saveZarr",
            True,
            "blockSize",
            [8, 8, 4],
        ),
        frame_case(
            "frame.tif.uint16.overwrite",
            u16,
            "im.tif",
            "o.tif",
            [2, 2, 2, 9, 9, 9],
            "overwrite",
            True,
        ),
        frame_case(
            "frame.tif.uint16.y_before_image",
            u16,
            "im.tif",
            "o.tif",
            [0, 2, 2, 9, 9, 9],
            error=True,
        ),
        frame_case(
            "frame.tif.uint16.z_before_image",
            u16,
            "im.tif",
            "o.tif",
            [2, 2, 0, 9, 9, 9],
            error=True,
        ),
        frame_case("frame.tif.uint16.missing_input", u16, "im.tif", "o.tif", [2, 2, 2, 9, 9, 9]),
    ]


def _zarr_cases():
    s = (23, 37, 29)
    u16 = zarr_tree(s, "uint16", "z1", (5, 7, 6))
    sgl = zarr_tree(s, "single", "z2", (5, 7, 6), True)
    small = ["batchSize", [16, 12, 10], "blockSize", [8, 6, 5]]
    return [
        zarr_case("zarr.uint16.default", u16, Work("crop.zarr"), [2, 3, 2, 37, 29, 23], *small),
        zarr_case("zarr.single.default", sgl, Work("crop.zarr"), [4, 1, 3, 30, 29, 20], *small),
        zarr_case(
            "zarr.uint16.no_mips",
            u16,
            Work("crop.zarr"),
            [2, 3, 2, 37, 29, 23],
            *small,
            "saveMIP",
            False,
        ),
        # (MATLAB pads a box that starts before the image without pad)
        zarr_case(
            "zarr.uint16.start_before_no_pad",
            u16,
            Work("crop.zarr"),
            [-3, 0, -1, 20, 20, 10],
            *small,
        ),
        zarr_case(
            "zarr.uint16.past_end_no_pad",
            u16,
            Work("crop.zarr"),
            [2, 3, 2, 38, 29, 23],
            *small,
            error=True,
        ),
        zarr_case(
            "zarr.uint16.past_end_pad",
            u16,
            Work("crop.zarr"),
            [30, 25, 20, 45, 34, 26],
            *small,
            "pad",
            True,
        ),
        zarr_case(
            "zarr.uint16.outside_image_pad",
            u16,
            Work("crop.zarr"),
            [41, 2, 2, 50, 9, 9],
            *small,
            "pad",
            True,
        ),
    ]


def _trim_cases():
    s = (9, 14, 11)  # MATLAB 14 x 11 x 9
    im = {"im": Value(lambda: to_matlab(volume(s, "uint16", "trim")))}
    sgl = {"im": Value(lambda: to_matlab(volume(s, "single", "trim_s", True)))}
    cases = []
    for cid, args, inputs, err in (
        ("trim.uint16.both", [PerAxis([2, 3, 1], "yxz"), "both"], im, False),
        ("trim.uint16.pre", [PerAxis([2, 3, 1], "yxz"), "pre"], im, False),
        ("trim.uint16.post", [PerAxis([2, 3, 1], "yxz"), "post"], im, False),
        ("trim.uint16.default_both", [PerAxis([1, 0, 2], "yxz")], im, False),
        ("trim.single.both", [PerAxis([2, 3, 1], "yxz"), "both"], sgl, False),
        ("trim.uint16.zero", [PerAxis([0, 0, 0], "yxz"), "both"], im, False),
        ("trim.uint16.larger_than_image", [PerAxis([8, 1, 1], "yxz"), "both"], im, False),
        ("trim.uint16.post_scalar", [2, "post"], im, False),
        ("trim.uint16.both_scalar", [2, "both"], im, True),
        ("trim.uint16.negative", [PerAxis([-1, 0, 0], "yxz"), "pre"], im, True),
        ("trim.uint16.negative_post", [PerAxis([0, -1, 0], "yxz"), "post"], im, True),
        ("trim.uint16.unknown_method", [PerAxis([1, 1, 1], "yxz"), "middle"], im, True),
    ):
        cases.append(
            Case(
                cid,
                matlab="trimBorder",
                args=[Vol(In("im")), *args],
                python=crop.trim_border,
                out_roles=None if err else ("vol",),
                inputs=inputs,
                expect_error=err,
                py_error=ERRORS,
            )
        )
    return cases


CASES = _dataset_cases() + _frame_cases() + _zarr_cases() + _trim_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# CropDatasetTest in intent, on (z, y, x) data and native zarr arrays
# ---------------------------------------------------------------------------------------------

PY_BOX = (2, 4, 6, 10, 30, 33)  # BBOX as a Python box
CHANNELS = ("CamA_ch0", "CamB_ch0")


def write_dataset(folder, fmt, cls="uint16", shape=IM_SIZE, chunks=(4, 16, 16)):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ims = []
    for k, n in enumerate(DATASET):
        a = volume(shape, cls, f"py/{fmt}/{cls}/{k}")
        ims.append(a)
        if fmt == "zarr":
            write_zarr(a, str(folder / f"{n}.zarr"), block_size=chunks)
        else:
            write_tiff(a, str(folder / f"{n}.tif"))
    return str(folder) + "/", ims


def cut(im, b):
    return im[b[0] : b[3], b[1] : b[4], b[2] : b[5]]


def pad_cut(im, b):
    """the box with zeros outside the image (CropDatasetTest's padCrop)"""
    out = np.zeros(tuple(b[q + 3] - b[q] for q in range(3)), dtype=im.dtype)
    s = [max(b[q], 0) for q in range(3)]
    t = [min(b[q + 3], im.shape[q]) for q in range(3)]
    o = [s[q] - b[q] for q in range(3)]
    out[o[0] : o[0] + t[0] - s[0], o[1] : o[1] + t[1] - s[1], o[2] : o[2] + t[2] - s[2]] = im[
        s[0] : t[0], s[1] : t[1], s[2] : t[2]
    ]
    return out


LOCAL_KW = {"parse_cluster": False, "master_compute": True}
LARGE_KW = {
    "zarr_file": True,
    "large_file": True,
    "batch_size": (4, 16, 16),
    "block_size": (4, 8, 8),
}


def test_tiff_to_tiff(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "tif")
    crop.crop_dataset(path, PY_BOX, channel_patterns=CHANNELS, **LOCAL_KW)
    for i in range(2):
        same_bits(read_tiff(f"{path}Cropped/{DATASET[i]}.tif"), cut(ims[i], PY_BOX))
    assert not os.path.exists(f"{path}Cropped/{DATASET[2]}.tif"), "a file of no channel"
    assert os.path.isfile(f"{path}Cropped/parameters.mat")
    pr = read_mat(f"{path}Cropped/parameters.mat")["pr"]
    # (MATLAB's values: the box 1-based [y x z])
    assert [float(v) for v in np.asarray(pr["inputBbox"]).ravel()] == BBOX


def test_zarr_to_zarr(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    crop.crop_dataset(
        path, PY_BOX, channel_patterns=CHANNELS, zarr_file=True, save_zarr=True, **LOCAL_KW
    )
    for i in range(2):
        same_bits(read_zarr(f"{path}Cropped/{DATASET[i]}.zarr"), cut(ims[i], PY_BOX))


def test_zarr_to_tiff(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    crop.crop_dataset(path, PY_BOX, channel_patterns=CHANNELS, zarr_file=True, **LOCAL_KW)
    for i in range(2):
        same_bits(read_tiff(f"{path}Cropped/{DATASET[i]}.tif"), cut(ims[i], PY_BOX))


def test_large_zarr_matches_in_memory_crop(tmp_path):
    # small batches so the crop is assembled from several blocks
    path, ims = write_dataset(tmp_path / "d", "zarr")
    crop.crop_dataset(path, PY_BOX, channel_patterns=CHANNELS, **LARGE_KW, **LOCAL_KW)
    for i in range(2):
        same_bits(read_zarr(f"{path}Cropped/{DATASET[i]}.zarr"), cut(ims[i], PY_BOX))


def test_large_zarr_makes_mips_before_publishing(tmp_path):
    # the job framework counts a crop as done, and may cancel its job, once the result exists
    from petakit5d.testing.publish_spy import PublishSpy
    from petakit5d.tools.mip import axis_max

    path, ims = write_dataset(tmp_path / "d", "zarr")
    with PublishSpy() as spy:
        crop.crop_dataset(path, PY_BOX, channel_patterns=CHANNELS, **LARGE_KW, **LOCAL_KW)
    for i in range(2):
        fn = f"{path}Cropped/{DATASET[i]}.zarr"
        assert spy.was_moved(fn), f"{DATASET[i]} was not moved into place"
        when = spy.files_when_moved(fn)
        for lab in ("y", "x", "z"):
            assert any(w.endswith(f"MIPs/{DATASET[i]}_MIP_{lab}.tif") for w in when), lab
        mip = read_tiff(f"{path}Cropped/MIPs/{DATASET[i]}_MIP_z.tif")
        same_bits(mip, axis_max(cut(ims[i], PY_BOX), -3)[0])


def test_large_zarr_mips_are_merged_from_the_blocks(tmp_path, monkeypatch):
    # the crop tasks save the MIPs of their blocks in their flags, and mip_zarr merges them
    # into the MIP slabs instead of reading the result again in mip_block tasks
    from petakit5d.testing.publish_spy import PublishSpy
    from petakit5d.tools import mip

    def no_blocks(*a, **k):
        raise AssertionError("The MIPs were computed from the result.")

    monkeypatch.setattr(mip, "mip_block", no_blocks)
    path, ims = write_dataset(tmp_path / "d", "zarr")
    with PublishSpy() as spy:
        crop.crop_dataset(path, PY_BOX, channel_patterns=CHANNELS[:1], **LARGE_KW, **LOCAL_KW)
    flags = [c for c in spy.calls if "/zarr_flag/" in c.dst and c.dst.endswith(".mat")]
    assert flags, "The block flags were not moved into place."
    sub = cut(ims[0], PY_BOX)
    # (MATLAB's orientation: the y MIP is nx x nz, the x MIP ny x nz)
    want = {"y": sub.max(axis=1).T, "x": sub.max(axis=2).T, "z": sub.max(axis=0)}
    for lab, m in want.items():
        same_bits(read_tiff(f"{path}Cropped/MIPs/{DATASET[0]}_MIP_{lab}.tif"), m)


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls (MATLAB's FrameworkSpyFixture); nothing runs. A test may set
    ``spy.action`` (called at each call) and ``spy.done`` (the done flags returned)."""
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        if fake.action is not None:
            fake.action()
        return [fake.done] * len(list(outputs))

    fake.action = None
    fake.done = False
    monkeypatch.setattr(framework, "run_tasks", fake)
    fake.calls = calls
    return fake


def test_large_zarr_blocks_get_the_mip_option_only_for_mips(tmp_path, framework_spy):
    path, _ = write_dataset(tmp_path / "d", "zarr")
    framework_spy.done = True
    for save in (False, True):
        crop.crop_zarr(
            f"{path}{DATASET[0]}.zarr",
            f"{path}crop_{int(save)}.zarr",
            PY_BOX,
            batch_size=(4, 16, 16),
            block_size=(4, 8, 8),
            save_mip=save,
            parse_cluster=False,
        )
    c0, c1 = framework_spy.calls[0], framework_spy.calls[1]
    assert all("mip_batch_size" not in t.kwargs for t in c0["tasks"])
    assert all(t.kwargs["mip_batch_size"] == [8, 26, 27] for t in c1["tasks"])
    # 8 batches, 5 per task: MATLAB's flag names, 1-based
    names = [os.path.basename(o) for o in c1["outputs"]]
    assert names == ["blocks_1_5.mat", "blocks_6_8.mat"]
    assert [t.args[0] for t in c1["tasks"]] == [[0, 1, 2, 3, 4], [5, 6, 7]]


def test_large_zarr_keeps_a_result_another_writer_published(tmp_path, framework_spy):
    # two job sets on the same data: this one keeps the other's crop and removes its own
    # (incomplete here: the spy reports the blocks done while the other crop comes)
    path, _ = write_dataset(tmp_path / "d", "zarr")
    crop_fn = f"{path}crop.zarr"
    other = np.full((5, 5, 5), 7, dtype=np.uint16)
    framework_spy.done = True
    framework_spy.action = lambda: write_zarr(other, crop_fn)
    crop.crop_zarr(
        f"{path}{DATASET[0]}.zarr",
        crop_fn,
        PY_BOX,
        batch_size=(4, 16, 16),
        block_size=(4, 8, 8),
        parse_cluster=False,
    )
    same_bits(read_zarr(crop_fn), other)
    assert not list(Path(path).glob("crop.zarr_*")), "the temporary crop was not removed"
    assert not list(Path(path, "zarr_flag").glob("crop_*")), "the flags were not removed"
    assert not list(Path(path).glob("MIPs/*.tif")), "MIPs were made from the incomplete crop"
    # (and the blocks stop once the crop is there)
    assert framework_spy.calls[0]["final_out_fullpath"] == crop_fn


def test_large_zarr_keeps_a_result_published_during_its_mips(tmp_path, monkeypatch):
    from petakit5d.tools import mip

    path, _ = write_dataset(tmp_path / "d", "zarr")
    crop_fn = f"{path}crop.zarr"
    other = np.full((5, 5, 5), 7, dtype=np.uint16)
    monkeypatch.setattr(mip, "mip_zarr", lambda *a, **k: write_zarr(other, crop_fn))
    crop.crop_zarr(
        f"{path}{DATASET[0]}.zarr",
        crop_fn,
        PY_BOX,
        batch_size=(4, 16, 16),
        block_size=(4, 8, 8),
        parse_cluster=False,
    )
    same_bits(read_zarr(crop_fn), other)
    assert not list(Path(path).glob("crop.zarr_*")), "the temporary crop was not removed"
    assert not [n for n in os.listdir(crop_fn) if "zarr_" in n], "moved into the other's"


@pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="Writers outside Slurm are checked on Linux only."
)
def test_large_zarr_waits_for_a_writer_of_another_call(tmp_path, framework_spy, monkeypatch):
    # a crop of another call (another uuid) that may still make the result: this one waits
    # for it (a process on this node, whose crop comes during the wait); no block runs
    import subprocess

    from petakit5d.jobs import files, writers

    for v in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID"):
        monkeypatch.delenv(v, raising=False)
    writers.reset()
    path, _ = write_dataset(tmp_path / "d", "zarr")
    crop_fn = f"{path}crop.zarr"
    proc = subprocess.Popen(["sleep", "120"])
    try:
        flag = f"{path}zarr_flag/crop_otheruuid/"
        os.makedirs(flag)
        owner = f"0_:{writers.node_name()}:{proc.pid}:{writers.process_start(proc.pid)}"
        with open(flag + "writer.txt", "w", newline="") as f:
            f.write(f"result={crop_fn}\nuuid=otheruuid\nowner={owner}\n")
        st = os.stat(flag + "writer.txt")
        os.utime(flag + "writer.txt", ns=(st.st_atime_ns, st.st_mtime_ns - 60 * 10**9))
        other = np.full((5, 5, 5), 7, dtype=np.uint16)

        def spy(groups, wait_time):
            write_zarr(other, crop_fn)
            return False

        monkeypatch.setattr(files, "wait_for_files", spy)
        crop.crop_zarr(
            f"{path}{DATASET[0]}.zarr",
            crop_fn,
            PY_BOX,
            batch_size=(4, 16, 16),
            block_size=(4, 8, 8),
            parse_cluster=False,
        )
    finally:
        proc.kill()
        proc.wait()
    assert framework_spy.calls == [], "the crop was made a second time"
    same_bits(read_zarr(crop_fn), other)


fake_slurm = slurm_fake.fake_slurm  # (the fixture)


@pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
def test_large_zarr_blocks_of_an_ended_task_are_skipped(tmp_path, fake_slurm):
    # the block jobs of a large crop whose job was cancelled ran on, for a result nothing
    # would publish; a block skips its work when its task's job has ended
    from petakit5d.jobs import writers

    writers.reset()
    path, _ = write_dataset(tmp_path / "d", "zarr")
    crop_fn = f"{path}crop.zarr"
    tmp_fn = f"{crop_fn}_u1"
    create_zarr(tmp_fn, data_size=(8, 26, 27), block_size=(4, 8, 8), dtype="uint16")
    flag = f"{path}zarr_flag/crop_u1/"
    os.makedirs(flag)
    with open(flag + "writer.txt", "w", newline="") as f:
        f.write(f"result={crop_fn}\nuuid=u1\nowner=777_3:othernode:5\n")
    box = (PY_BOX[0], PY_BOX[1], PY_BOX[2], PY_BOX[0] + 4, PY_BOX[1] + 8, PY_BOX[2] + 8)
    out = crop.crop_block(
        [0],
        f"{path}{DATASET[0]}.zarr",
        tmp_fn,
        flag + "blocks_1_1.mat",
        [box],
        [(0, 0, 0, 4, 8, 8)],
        uuid="u1",
    )
    assert out is None
    assert not os.path.isfile(flag + "blocks_1_1.mat"), "the block of an ended task was cropped"
    assert int(np.count_nonzero(read_zarr(tmp_fn))) == 0
    writers.reset()


def test_large_zarr_task_that_fails_marks_its_record_ended(tmp_path, framework_spy):
    path, _ = write_dataset(tmp_path / "d", "zarr")
    crop_fn = f"{path}crop.zarr"
    with pytest.raises(crop.CropError, match="not finished"):
        crop.crop_zarr(
            f"{path}{DATASET[0]}.zarr",
            crop_fn,
            PY_BOX,
            batch_size=(4, 16, 16),
            block_size=(4, 8, 8),
            parse_cluster=False,
        )
    records = list(Path(path, "zarr_flag").glob("crop_*/writer.txt"))
    assert len(records) == 1, "the crop did not record itself"
    assert "\nended=1" in records[0].read_text(), "its record was not marked ended"
    # (two tries, the second with twice the memory)
    mems = [c["mem_allocate"] for c in framework_spy.calls]
    assert len(mems) == 2 and mems[1] == pytest.approx(2 * mems[0])


@pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
def test_large_zarr_jobs_go_after_their_blocks(tmp_path, fake_slurm, framework_spy, monkeypatch):
    from petakit5d.jobs import coordinator

    create_zarr(
        str(tmp_path / "scan_CamA_ch0_t0000.zarr"),
        data_size=(64, 512, 512),
        block_size=(64, 64, 64),
    )
    monkeypatch.setattr(coordinator, "coordinator_nice", lambda slurm=None: 1000)
    (tmp_path / "logs").mkdir()
    for large in (True, False):
        crop.crop_dataset(
            str(tmp_path),
            (0, 0, 0, 64, 512, 512),
            channel_patterns=("CamA_ch0",),
            zarr_file=True,
            large_file=large,
            parse_cluster=True,
            result_dir_name=f"Cropped_{int(large)}",
            job_log_dir=str(tmp_path / "logs"),
        )
    assert framework_spy.calls[0]["nice"] == 1000
    assert framework_spy.calls[1]["nice"] == 0


def test_large_zarr_jobs_ask_for_what_they_hold(tmp_path, framework_spy):
    # the job of a large crop (crop_zarr) runs its blocks as jobs (1 CPU each); it asks for
    # the process and the merge of the MIP slabs with one CPU, and if it crops blocks itself
    # (master_compute, the config's if it sets it), for cpus_per_task (2) and a block's
    # estimate (a 512 x 1024 x 1024 uint16 batch, 1 GiB, times BLOCK_MEMORY_FACTOR; MATLAB:
    # 8 bytes per voxel, 4 GB). The input is an 8 G-voxel header.
    create_zarr(
        str(tmp_path / "scan_CamA_ch0_t0000.zarr"),
        data_size=(512, 4096, 4096),
        block_size=(256, 256, 256),
    )
    off = tmp_path / "cpu_config.json"
    off.write_text('{"masterCompute": false}')
    for i, (mc, cfg) in enumerate([(True, ""), (False, ""), (True, str(off))]):
        crop.crop_dataset(
            str(tmp_path),
            (0, 0, 0, 512, 4096, 4096),
            channel_patterns=("CamA_ch0",),
            zarr_file=True,
            large_file=True,
            parse_cluster=False,
            master_compute=mc,
            config_file=cfg,
            result_dir_name=f"Cropped_{i}",
        )
    c = framework_spy.calls
    assert c[0]["cpus_per_task"] == 2
    assert c[0]["mem_allocate"] == pytest.approx(crop.BLOCK_MEMORY_FACTOR, abs=1e-9)
    for i in (1, 2):
        assert c[i]["cpus_per_task"] == 1, i
        assert c[i]["mem_allocate"] < crop.BLOCK_MEMORY_FACTOR, i
    assert all("master_compute" not in t.kwargs for t in c[0]["tasks"])
    assert all(t.kwargs["master_compute"] is False for t in c[1]["tasks"])


def test_block_jobs_ask_for_the_measured_factor_of_their_batch(tmp_path, framework_spy):
    # MATLAB asks for 8 bytes per voxel of a batch whatever the class; the port for
    # BLOCK_MEMORY_FACTOR times the batch in its class, doubled at the second try
    for cls, byte in (("uint8", 1), ("uint16", 2), ("single", 4)):
        fn = str(tmp_path / f"{cls}.zarr")
        create_zarr(fn, data_size=(64, 512, 512), block_size=(64, 64, 64), dtype=cls)
        with pytest.raises(crop.CropError):
            crop.crop_zarr(
                fn,
                str(tmp_path / f"c_{cls}.zarr"),
                (0, 0, 0, 64, 512, 512),
                batch_size=(64, 256, 256),
                parse_cluster=False,
            )
        gb = 64 * 256 * 256 * byte / 2**30 * crop.BLOCK_MEMORY_FACTOR
        got = [c["mem_allocate"] for c in framework_spy.calls[-2:]]
        assert got == pytest.approx([gb, 2 * gb], rel=1e-12)


def test_large_zarr_passes_job_settings_to_mips(tmp_path, monkeypatch):
    from petakit5d.tools import mip

    calls = []
    monkeypatch.setattr(mip, "mip_zarr", lambda *a, **k: calls.append(k))
    path, _ = write_dataset(tmp_path / "d", "zarr")
    crop.crop_zarr(
        f"{path}{DATASET[0]}.zarr",
        f"{path}crop.zarr",
        PY_BOX,
        batch_size=(4, 16, 16),
        block_size=(4, 8, 8),
        parse_cluster=False,
        master_compute=False,
        job_log_dir=f"{path}logs",
    )
    assert calls, "mip_zarr was not called"
    assert calls[0]["parse_cluster"] is False
    assert calls[0]["master_compute"] is False
    assert calls[0]["job_log_dir"] == f"{path}logs"
    assert len(calls[0]["block_flag_fullpaths"]) > 0


PAD_BOX = (-1, -3, -2, 15, 45, 40)  # MATLAB's [-2, -1, 0, 45, 40, 15]


def test_pad_past_image_edge(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "tif")
    crop.crop_dataset(path, PAD_BOX, channel_patterns=CHANNELS, pad=True, **LOCAL_KW)
    same_bits(read_tiff(f"{path}Cropped/{DATASET[0]}.tif"), pad_cut(ims[0], PAD_BOX))


def test_pad_past_image_edge_large_zarr(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    crop.crop_dataset(path, PAD_BOX, channel_patterns=CHANNELS, pad=True, **LARGE_KW, **LOCAL_KW)
    same_bits(read_zarr(f"{path}Cropped/{DATASET[0]}.zarr"), pad_cut(ims[0], PAD_BOX))


def test_clip_past_image_edge_without_pad(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "tif")
    crop.crop_dataset(path, (7, 29, 19, 20, 60, 50), channel_patterns=CHANNELS, **LOCAL_KW)
    same_bits(read_tiff(f"{path}Cropped/{DATASET[0]}.tif"), ims[0][7:, 29:, 19:])


def test_single_slice_zarr_crop(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    box = (5, 4, 6, 6, 30, 33)
    crop.crop_dataset(
        path,
        box,
        channel_patterns=CHANNELS[:1],
        zarr_file=True,
        save_zarr=True,
        result_dir_name="Cropped_small",
        **LOCAL_KW,
    )
    crop.crop_dataset(
        path,
        box,
        channel_patterns=CHANNELS[:1],
        **LARGE_KW,
        result_dir_name="Cropped_large",
        **LOCAL_KW,
    )
    for r in ("Cropped_small", "Cropped_large"):
        got = read_zarr(f"{path}{r}/{DATASET[0]}.zarr")
        same_bits(got.reshape(cut(ims[0], box).shape), cut(ims[0], box))


def test_center_crop(tmp_path):
    from petakit5d.core.dtypes import mround

    path, ims = write_dataset(tmp_path / "d", "tif")
    size = (4, 10, 9)
    crop.crop_dataset(
        path, (0, 0, 0, *size), crop_type="center", channel_patterns=CHANNELS[:1], **LOCAL_KW
    )
    # start = round((imSize + 1) / 2) - floor(bboxSize / 2) (1-based)
    s = [int(mround((n + 1) / 2)) - e // 2 - 1 for n, e in zip(IM_SIZE, size, strict=True)]
    box = (*s, *(a + e for a, e in zip(s, size, strict=True)))
    same_bits(read_tiff(f"{path}Cropped/{DATASET[0]}.tif"), cut(ims[0], box))


def test_moving_crop(tmp_path):
    # the box moves linearly from its start to last_start_coords over a channel's time points
    from petakit5d.core.dtypes import mround

    path = str(tmp_path) + "/"
    nt = 3
    ims = [volume(IM_SIZE, "uint16", f"mv{t}") for t in range(nt)]
    for t in range(nt):
        write_tiff(ims[t], f"{path}scan_CamA_ch0_t{t:04d}.tif")
    box = (1, 2, 2, 6, 12, 12)
    last = (5, 20, 12)
    crop.crop_dataset(
        path,
        box,
        crop_type="moving",
        last_start_coords=last,
        channel_patterns=CHANNELS[:1],
        **LOCAL_KW,
    )
    size = [box[q + 3] - box[q] for q in range(3)]
    for t in range(nt):
        s = [box[q] + int(mround((last[q] - box[q]) / (nt - 1) * t)) for q in range(3)]
        got = read_tiff(f"{path}Cropped/scan_CamA_ch0_t{t:04d}.tif")
        same_bits(got, cut(ims[t], (*s, *(a + e for a, e in zip(s, size, strict=True)))))


def test_moving_crop_counts_time_points_per_channel(tmp_path):
    # two channels of 2 time points each, the files sorted CamA t0, t1, CamB t0, t1: each
    # file's time point is its number in its channel (MATLAB's cumsum of the channel flags)
    path = str(tmp_path) + "/"
    names = [f"scan_Cam{c}_ch0_t{t:04d}" for c in "AB" for t in range(2)]
    ims = {n: volume(IM_SIZE, "uint16", f"mc/{n}") for n in names}
    for n, a in ims.items():
        write_tiff(a, f"{path}{n}.tif")
    box = (1, 2, 2, 6, 12, 12)
    last = (3, 10, 6)
    crop.crop_dataset(
        path,
        box,
        crop_type="moving",
        last_start_coords=last,
        channel_patterns=CHANNELS,
        **LOCAL_KW,
    )
    for n in names:
        s = box[:3] if n.endswith("t0000") else last
        got = read_tiff(f"{path}Cropped/{n}.tif")
        same_bits(got, cut(ims[n], (*s, *(a + 5 if q == 0 else a + 10 for q, a in enumerate(s)))))


def test_center_crop_of_odd_sizes(tmp_path):
    # MATLAB's center round((imSize + 1) / 2) of odd sizes: (11, 39, 35) -> 6, 20, 18 (1-based)
    path = str(tmp_path) + "/"
    im = volume((11, 39, 35), "uint16", "odd")
    write_tiff(im, f"{path}scan_CamA_ch0_t0000.tif")
    crop.crop_dataset(
        path, (0, 0, 0, 3, 8, 5), crop_type="center", channel_patterns=CHANNELS[:1], **LOCAL_KW
    )
    # start = center - floor(size / 2): 6 - 1, 20 - 4, 18 - 2 (1-based), 0-based one less
    same_bits(read_tiff(f"{path}Cropped/scan_CamA_ch0_t0000.tif"), im[4:7, 15:23, 15:20])


def test_moving_crop_of_two_datasets_fails_as_matlab(tmp_path):
    # MATLAB numbers the files over all datasets but takes the time points of the last one
    for d in ("a", "b"):
        for t in range(2):
            (tmp_path / d).mkdir(exist_ok=True)
            write_tiff(
                volume(IM_SIZE, "uint16", f"q{d}{t}"),
                str(tmp_path / d / f"scan_CamA_ch0_t{t:04d}.tif"),
            )
    with pytest.raises(crop.CropError) as err:
        crop.crop_dataset(
            [str(tmp_path / "a"), str(tmp_path / "b")],
            (1, 2, 2, 6, 12, 12),
            crop_type="moving",
            last_start_coords=(2, 3, 3),
            channel_patterns=CHANNELS[:1],
            **LOCAL_KW,
        )
    assert err.value.identifier == "MATLAB:badsubscript"


def test_existing_output_is_not_overwritten(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "tif")
    os.makedirs(f"{path}Cropped")
    sentinel = np.arange(16, dtype=np.uint16).reshape(4, 4)
    write_tiff(sentinel, f"{path}Cropped/{DATASET[0]}.tif")
    crop.crop_dataset(path, PY_BOX, channel_patterns=CHANNELS, **LOCAL_KW)
    same_bits(read_tiff(f"{path}Cropped/{DATASET[0]}.tif"), sentinel)
    same_bits(read_tiff(f"{path}Cropped/{DATASET[1]}.tif"), cut(ims[1], PY_BOX))


def test_tiff_input_saved_as_zarr(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "tif")
    crop.crop_dataset(path, PY_BOX, channel_patterns=CHANNELS[:1], save_zarr=True, **LOCAL_KW)
    fn = f"{path}Cropped/{DATASET[0]}.zarr"
    assert os.path.isdir(fn), "the zarr result was not saved as <name>.zarr"
    same_bits(read_zarr(fn), cut(ims[0], PY_BOX))


# ---------------------------------------------------------------------------------------------
# CropMexTest.trimBorderModes, and trim_border's rules
# ---------------------------------------------------------------------------------------------


def test_trim_border_modes():
    im = volume((9, 14, 11), "uint16", "trim_modes")
    bs = (1, 2, 3)  # MATLAB's [2, 3, 1]
    same_bits(crop.trim_border(im, bs, "both"), im[1:-1, 2:-2, 3:-3])
    same_bits(crop.trim_border(im, bs, "pre"), im[1:, 2:, 3:])
    same_bits(crop.trim_border(im, bs, "post"), im[:-1, :-2, :-3])
    same_bits(crop.trim_border(im, bs), im[1:-1, 2:-2, 3:-3])
    assert crop.trim_border(im, (1, 2, 3)).flags.c_contiguous


def test_trim_border_edge_cases():
    im = volume((9, 14, 11), "uint16", "trim_edges")
    assert crop.trim_border(im, (0, 0, 0)) is im
    assert crop.trim_border(im, (5, 1, 1)).shape == (0, 12, 9)
    assert crop.trim_border(im, (12, 0, 0), "post").shape == (0, 14, 11)
    same_bits(crop.trim_border(im, 2, "post"), im[:-2, :-2, :-2])
    for args, ident in (
        ((2, "both"), "MATLAB:badsubscript"),
        (((-1, 0, 0), "pre"), "MATLAB:badsubscript"),
        (((0, -1, 0), "post"), "MATLAB:badsubscript"),
        (((1, 1, 1), "middle"), "MATLAB:UndefinedFunction"),
    ):
        with pytest.raises(crop.CropError) as err:
            crop.trim_border(im, *args)
        assert err.value.identifier == ident, args


# ---------------------------------------------------------------------------------------------
# the port's own rules
# ---------------------------------------------------------------------------------------------


def test_matlab_layer_converts_boxes_and_sizes(tmp_path):
    from petakit5d import matlab

    path, ims = write_dataset(tmp_path / "d", "zarr")
    matlab.XR_crop_dataset(
        path,
        BBOX,
        "channelPatterns",
        ["CamA_ch0"],
        "zarrFile",
        True,
        "largeFile",
        True,
        "batchSize",
        [16, 16, 4],
        "blockSize",
        [8, 8, 4],
        "parseCluster",
        False,
    )
    fn = f"{path}Cropped/{DATASET[0]}.zarr"
    same_bits(read_zarr(fn), cut(ims[0], PY_BOX))
    from petakit5d.io import get_zarr_info

    assert tuple(get_zarr_info(fn).block_size) == (4, 8, 8)


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_native_and_matlab_layout_inputs_give_the_same_crops(tmp_path, cls):
    im = volume((23, 37, 29), cls, "layouts", True)
    box = (-2, 3, 1, 20, 40, 25)
    out = []
    for layout in ("native", "matlab"):
        d = tmp_path / layout
        d.mkdir()
        synth.write_zarr(d / "im.zarr", im, block_size=(5, 7, 6), layout=layout)
        crop.crop_zarr(
            str(d / "im.zarr"),
            str(d / "c.zarr"),
            box,
            pad=True,
            batch_size=(10, 16, 12),
            block_size=(5, 8, 6),
            parse_cluster=False,
        )
        out.append(read_zarr(str(d / "c.zarr")))
        for lab in ("y", "x", "z"):
            out.append(read_tiff(str(d / f"MIPs/c_MIP_{lab}.tif")))
    for a, b in zip(out[:4], out[4:], strict=True):
        same_bits(a, b)
    same_bits(out[0], pad_cut(im, box))


@pytest.mark.parametrize("threads", [1, 3])
def test_results_do_not_depend_on_threads(tmp_path, threads, monkeypatch):
    from petakit5d.core import threads as th

    monkeypatch.setattr(th, "_AVAILABLE", threads)
    im = volume((23, 37, 29), "single", "threads", True)
    synth.write_zarr(tmp_path / "im.zarr", im, block_size=(5, 7, 6))
    crop.crop_zarr(
        str(tmp_path / "im.zarr"),
        str(tmp_path / "c.zarr"),
        (1, 2, 3, 20, 30, 25),
        batch_size=(10, 16, 12),
        block_size=(5, 8, 6),
        parse_cluster=False,
    )
    same_bits(read_zarr(str(tmp_path / "c.zarr")), im[1:20, 2:30, 3:25])


def test_crop_zarr_refuses_a_box_past_the_end_without_pad(tmp_path):
    path, _ = write_dataset(tmp_path / "d", "zarr")
    with pytest.raises(crop.CropError, match="out of bound"):
        crop.crop_zarr(
            f"{path}{DATASET[0]}.zarr", f"{path}c.zarr", (0, 0, 0, 13, 40, 36), parse_cluster=False
        )
    # (and pads a box that starts before the image without pad, as its blocks always pad)
    crop.crop_zarr(
        f"{path}{DATASET[0]}.zarr",
        f"{path}c.zarr",
        (-1, -2, -3, 10, 30, 30),
        batch_size=(4, 16, 16),
        block_size=(4, 8, 8),
        parse_cluster=False,
    )
    im = read_zarr(f"{path}{DATASET[0]}.zarr")
    same_bits(read_zarr(f"{path}c.zarr"), pad_cut(im, (-1, -2, -3, 10, 30, 30)))


def test_crop_zarr_skips_an_existing_result(tmp_path, framework_spy):
    path, _ = write_dataset(tmp_path / "d", "zarr")
    os.makedirs(f"{path}c.zarr")
    assert (
        crop.crop_zarr(f"{path}{DATASET[0]}.zarr", f"{path}c.zarr", PY_BOX, parse_cluster=False)
        is None
    )
    assert framework_spy.calls == []


def test_many_chunks_use_the_nested_separator(tmp_path, framework_spy):
    # MATLAB: '/' when the result has more than 10,000 chunks
    import json

    fn = str(tmp_path / "im.zarr")
    create_zarr(fn, data_size=(1, 101, 101), block_size=(1, 1, 1))
    framework_spy.done = True
    for k, box in enumerate(((0, 0, 0, 1, 101, 101), (0, 0, 0, 1, 100, 100))):
        crop.crop_zarr(
            fn,
            str(tmp_path / f"c{k}.zarr"),
            box,
            block_size=(1, 1, 1),
            save_mip=False,
            parse_cluster=False,
            uuid=f"u{k}",
        )
        meta = json.loads(Path(tmp_path / f"c{k}.zarr" / ".zarray").read_text())
        assert meta.get("dimension_separator", ".") == ("/" if k == 0 else "."), k


def test_crop_block_skips_done_batches_and_checks_its_folders(tmp_path, capsys):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    src = f"{path}{DATASET[0]}.zarr"
    tmp_fn = f"{path}c.zarr_u2"
    create_zarr(tmp_fn, data_size=(4, 8, 8), block_size=(4, 8, 8), dtype="uint16")
    box = [(0, 0, 0, 4, 8, 8)]
    with pytest.raises(crop.CropError, match="does not exist"):
        crop.crop_block([0], src, tmp_fn, f"{path}nofolder/blocks_1_1.mat", box, box)
    os.makedirs(f"{path}zarr_flag/c_u2")
    flag = f"{path}zarr_flag/c_u2/blocks_1_1.mat"
    with pytest.raises(crop.CropError, match="doesnot exist"):
        crop.crop_block([0], f"{path}none.zarr", tmp_fn, flag, box, box)
    assert crop.crop_block([0], src, tmp_fn, flag, box, box) == [True]
    assert os.path.getsize(flag) == 0
    same_bits(read_zarr(tmp_fn), ims[0][:4, :8, :8])
    assert crop.crop_block([0], src, tmp_fn, flag, box, box) is None
    assert "already exist, skip them" in capsys.readouterr().out
    # overwrite: the batch again, the flag with the MIPs of the blocks
    assert crop.crop_block(
        [0], src, tmp_fn, flag, box, box, overwrite=True, mip_batch_size=(4, 8, 8)
    ) == [True]
    flags = read_mat(flag)
    assert set(flags) >= {"blockBboxes", "mipBatchSize", "zMIPs"}


def test_crop_block_fills_batches_outside_the_image_with_zeros(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    src = f"{path}{DATASET[0]}.zarr"
    tmp_fn = f"{path}c.zarr_u3"
    create_zarr(tmp_fn, data_size=(4, 16, 8), block_size=(4, 8, 8), dtype="uint16")
    os.makedirs(f"{path}zarr_flag/c_u3")
    # one batch inside, one partly outside, one wholly outside the image
    inb = [(0, 0, 0, 4, 8, 8), (10, 36, 32, 14, 44, 40)]
    outb = [(0, 0, 0, 4, 8, 8), (0, 8, 0, 4, 16, 8)]
    crop.crop_block([0, 1], src, tmp_fn, f"{path}zarr_flag/c_u3/blocks_1_2.mat", inb, outb)
    got = read_zarr(tmp_fn)
    same_bits(got[:, :8], ims[0][:4, :8, :8])
    same_bits(got[:, 8:], pad_cut(ims[0], inb[1]))
    tmp2 = f"{path}c.zarr_u4"
    create_zarr(tmp2, data_size=(4, 8, 8), block_size=(4, 8, 8), dtype="uint16")
    os.makedirs(f"{path}zarr_flag/c_u4")
    crop.crop_block(
        [0],
        src,
        tmp2,
        f"{path}zarr_flag/c_u4/blocks_1_1.mat",
        [(20, 50, 50, 24, 58, 58)],
        [(0, 0, 0, 4, 8, 8)],
    )
    assert int(np.count_nonzero(read_zarr(tmp2))) == 0
    # batches that touch the image only at its edges: starting at its end, ending at its start
    for k, ib in enumerate(((12, 0, 0, 16, 8, 8), (-4, 0, 0, 0, 8, 8), (0, -8, 0, 4, 0, 8))):
        tmp3 = f"{path}c.zarr_e{k}"
        create_zarr(tmp3, data_size=(4, 8, 8), block_size=(4, 8, 8), dtype="uint16")
        os.makedirs(f"{path}zarr_flag/c_e{k}")
        done = crop.crop_block(
            [0], src, tmp3, f"{path}zarr_flag/c_e{k}/blocks_1_1.mat", [ib], [(0, 0, 0, 4, 8, 8)]
        )
        assert done == [True], ib
        assert int(np.count_nonzero(read_zarr(tmp3))) == 0, ib


def test_crop_file_pads_with_its_class(tmp_path):
    im = volume((5, 6, 7), "uint16", "padval")
    write_tiff(im, str(tmp_path / "im.tif"))
    box = (-1, 0, 0, 5, 6, 8)
    for k, (v, want) in enumerate(((7.5, 8), (70000, 65535), (-3, 0), (2.4999, 2))):
        out = str(tmp_path / f"o{k}.tif")
        crop.crop_file(str(tmp_path / "im.tif"), out, box, pad=True, pad_value=v)
        got = read_tiff(out)
        assert got.shape == (6, 6, 8)
        assert np.all(got[0] == want) and np.all(got[:, :, 7] == want), (v, want)
        same_bits(got[1:, :, :7], im)
    # single keeps the pad value's bits: -0 stays -0 (cast(-0, 'single'))
    sgl = volume((5, 6, 7), "single", "padval_s")
    write_tiff(sgl, str(tmp_path / "s.tif"))
    for k, v in enumerate((-0.0, 0.0)):
        out = str(tmp_path / f"s{k}.tif")
        crop.crop_file(str(tmp_path / "s.tif"), out, box, pad=True, pad_value=v)
        got = read_tiff(out)
        assert np.all(np.signbit(got[0]) == np.signbit(v)) and np.all(got[0] == 0), v
        same_bits(got[1:, :, :7], sgl)


def test_crop_file_skips_missing_inputs_and_existing_results(tmp_path, capsys):
    im = volume((5, 6, 7), "uint16", "skip")
    write_tiff(im, str(tmp_path / "im.tif"))
    with pytest.warns(UserWarning, match="does not exist"):
        crop.crop_file(str(tmp_path / "none.tif"), str(tmp_path / "o.tif"), (0, 0, 0, 2, 2, 2))
    assert not os.path.exists(tmp_path / "o.tif")
    crop.crop_file(str(tmp_path / "im.tif"), str(tmp_path / "o.tif"), (0, 0, 0, 2, 2, 2))
    crop.crop_file(str(tmp_path / "im.tif"), str(tmp_path / "o.tif"), (0, 0, 0, 3, 3, 3))
    assert "already exist, skip it" in capsys.readouterr().out
    assert read_tiff(str(tmp_path / "o.tif")).shape == (2, 2, 2)
    crop.crop_file(
        str(tmp_path / "im.tif"), str(tmp_path / "o.tif"), (0, 0, 0, 3, 3, 3), overwrite=True
    )
    assert read_tiff(str(tmp_path / "o.tif")).shape == (3, 3, 3)


def test_crop_file_overwrite_replaces_a_zarr_result(tmp_path):
    # MATLAB's movefile moved the new folder into the old one; the port replaces it
    im = volume((5, 6, 7), "uint16", "ow")
    write_tiff(im, str(tmp_path / "im.tif"))
    out = str(tmp_path / "o.zarr")
    crop.crop_file(str(tmp_path / "im.tif"), out, (0, 0, 0, 2, 2, 2), save_zarr=True)
    crop.crop_file(
        str(tmp_path / "im.tif"), out, (0, 0, 0, 3, 4, 5), save_zarr=True, overwrite=True
    )
    same_bits(read_zarr(out), im[:3, :4, :5])
    assert [n for n in os.listdir(out) if n.endswith(".zarr")] == []


def test_crop_file_errors_for_a_box_before_the_image_without_pad(tmp_path):
    im = volume((5, 6, 7), "uint16", "err")
    write_tiff(im, str(tmp_path / "im.tif"))
    for box in ((0, -1, 0, 3, 3, 3), (-1, 0, 0, 3, 3, 3)):
        with pytest.raises(ERRORS):
            crop.crop_file(str(tmp_path / "im.tif"), str(tmp_path / "o.tif"), box)
    assert not os.path.exists(tmp_path / "o.tif")


def test_results_are_moved_into_place(tmp_path):
    from petakit5d.testing.publish_spy import PublishSpy

    im = volume((5, 6, 7), "uint16", "moved")
    write_tiff(im, str(tmp_path / "im.tif"))
    with PublishSpy() as spy:
        crop.crop_file(str(tmp_path / "im.tif"), str(tmp_path / "o.tif"), (0, 0, 0, 2, 2, 2))
        crop.crop_file(
            str(tmp_path / "im.tif"), str(tmp_path / "o.zarr"), (0, 0, 0, 2, 2, 2), save_zarr=True
        )
    assert spy.was_moved(tmp_path / "o.tif")
    assert spy.was_moved(tmp_path / "o.zarr")
    assert sorted(os.listdir(tmp_path)) == ["im.tif", "o.tif", "o.zarr"]


def test_block_memory_is_the_factor_of_the_batch():
    assert crop.block_memory((1024, 1024, 1024), "uint16") == pytest.approx(
        2 * crop.BLOCK_MEMORY_FACTOR
    )
    assert crop.block_memory((512, 512, 512), "single") == pytest.approx(
        0.5 * crop.BLOCK_MEMORY_FACTOR
    )
