"""Resampling (petakit5d.tools.resample, task T04).

Oracle cases against MATLAB: XR_resampleFrame (factors up and down per axis, non-integer
factors, one number, [xy z], each method, uint8, uint16 and single with NaN of two payloads,
-0 and 0, save16bit off, boxes, one plane, images with [y x] and [y x z] factors, zarr input
and output, unknown and capitalized methods, an existing result, constant images),
XR_resampleSingleZarr (each method, borders that are and are not multiples of the factor,
non-integer factors and factors that lose digits in the block tasks' text, upsampling, sizes
not divisible by the batch and chunk sizes, boxes, every factor 1; each case also gives the
whole image resized at once, so the cases pin where MATLAB's batches equal it),
XR_resample_dataset (tiff, zarr and large zarr datasets, the factors' 10 decimals, mean and
max, a box, no matching file), resampleZarrBlock's checks and XR_checkResampleSetting. Files
are compared byte for byte (tiff), by data and metadata after the layout mapping (zarr) and by
content (parameters.mat and .json).

Then ResampleTest (tests/tools) in intent on the port's (z, y, x) data and native zarr arrays,
the batches' placement, the factors' text at the task boundaries (a framework spy) and the
results another writer published. ResampleTest's imresize3_fast methods are N01's oracle cases
(tests/matlabcompat/test_resize.py, resample_test.* and fast.*).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.core.dtypes import NAN
from petakit5d.io import read_tiff, read_zarr, write_tiff, write_zarr
from petakit5d.matlabcompat import resize as R
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, File, InPath, Tree, Work
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import resample as rs

# ---------------------------------------------------------------------------------------------
# data: (z, y, x) volumes from counter-based random numbers (bit-stable everywhere)
# ---------------------------------------------------------------------------------------------

NAN32 = np.float32(NAN)
# a second NaN payload (MATLAB: typecast(uint32(hex2dec('7FC00001')), 'single'))
NAN32_B = np.array([0x7FC00001], dtype=np.uint32).view(np.float32)[0]
IM_SIZE = (12, 40, 36)  # ResampleTest's ImSize [40 36 12]
DATASET = ("scan_CamA_ch0_t0000", "scan_CamB_ch0_t0000", "scan_CamA_ch1_t0000")
CHANS = ["channelPatterns", ["CamA_ch0", "CamB_ch0"]]
ONE = ["channelPatterns", ["CamA_ch0"]]
LOCAL = ["parseCluster", False, "masterCompute", True]
SMALL = ["batchSize", [16, 16, 4], "blockSize", [8, 8, 2]]
ERRORS = (rs.ResampleError, ValueError, OSError, RuntimeError, IndexError)


def volume(shape, cls, name, special=False):
    """a deterministic (z, y, x) volume: integers in [0, max), single normal * 100 + 300; with
    special (single): NaN of two payloads, -0 and 0 (``special='one'``: MATLAB's NaN only)"""
    rng = CounterRNG(0, f"t04/{name}")
    shape = tuple(int(s) for s in shape)
    if cls == "uint8":
        return rng.integers(0, 255, shape, dtype=np.uint8)
    if cls == "uint16":
        return rng.integers(0, 60000, shape, dtype=np.uint16)
    a = np.ascontiguousarray((rng.normal(shape) * 100.0 + 300.0).astype(np.float32))
    if special:
        u = rng.uniform((3, *shape))
        a[u[0] < 0.03] = NAN32
        a[u[1] < 0.05] = np.float32(-0.0)
        a[u[2] < 0.03] = NAN32 if special == "one" else NAN32_B
    return a


def same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype, (a.dtype, b.dtype)
    assert a.shape == b.shape, (a.shape, b.shape)
    if a.dtype.kind == "f":
        assert a.view(f"u{a.itemsize}").tobytes() == b.view(f"u{b.itemsize}").tobytes()
    else:
        assert np.array_equal(a, b)


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def tif_file(shape, cls, name, special=False):
    def make(path):
        synth.write_tiff(path, volume(shape, cls, name, special))

    return make


def const_file(path):
    synth.write_tiff(path, np.full(IM_SIZE, 1000, dtype=np.uint16))


def zarr_file(shape, cls, name, chunks=(4, 16, 16), special=False):
    """a zarr array in MATLAB's layout (chunks (z, y, x))"""

    def make(path):
        synth.write_zarr(
            path, volume(shape, cls, name, special), block_size=chunks, layout="matlab"
        )

    return make


def zarr_tree(shape, cls, name, chunks=(4, 16, 16), special=False):
    """a folder with im.zarr"""

    def make(folder):
        zarr_file(shape, cls, name, chunks, special)(Path(folder) / "im.zarr")

    return make


def dataset_tree(fmt, cls="uint16", tag=""):
    """ResampleTest's writeDataset: three files, the third matching no channel pattern"""

    def make(folder):
        for k, n in enumerate(DATASET):
            a = volume(IM_SIZE, cls, f"ds/{fmt}/{cls}/{tag}{k}")
            if fmt == "zarr":
                synth.write_zarr(
                    Path(folder) / f"{n}.zarr", a, block_size=(4, 16, 16), layout="matlab"
                )
            else:
                synth.write_tiff(Path(folder) / f"{n}.tif", a)

    return make


def existing_tree(folder):
    synth.write_tiff(Path(folder) / "im.tif", volume(IM_SIZE, "uint16", "ex"))
    synth.write_tiff(Path(folder) / "o.tif", np.arange(16, dtype=np.uint16).reshape(1, 4, 4))


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------


def frame_case(cid, make, name, out, rf, *args, error=False):
    if name.endswith(".zarr"):
        # (a zarr input is a folder: d/im.zarr)
        return Case(
            cid,
            code="XR_resampleFrame([in{1}, 'im.zarr'], in{2 : end}); out = {};",
            nargout=0,
            args=[InPath("d"), Work(out), rf, *args],
            python=lambda d, *a: call_matlab("XR_resampleFrame", d + "im.zarr", *a),
            inputs={"d": Tree(make, "d")},
            expect_error=error,
            py_error=ERRORS,
        )
    return Case(
        cid,
        matlab="XR_resampleFrame",
        nargout=0,
        args=[InPath("im"), Work(out), rf, *args],
        python=lambda *a: call_matlab("XR_resampleFrame", *a),
        inputs={"im": File(make, name)},
        expect_error=error,
        py_error=ERRORS,
    )


def _frame_cases():
    u16 = tif_file(IM_SIZE, "uint16", "f/u16")
    u8 = tif_file(IM_SIZE, "uint8", "f/u8")
    sgl = tif_file(IM_SIZE, "single", "f/sgl", True)
    sgl1 = tif_file(IM_SIZE, "single", "f/sgl1", "one")
    im2 = tif_file((1, 40, 36), "uint16", "f/2d")
    factors = {
        "down2": [2, 2, 2],
        "z_only": [1, 1, 2],
        "up2xy": [0.5, 0.5, 1],
        "non_integer": [1.5, 1.5, 1.5],
        "down3": [3, 3, 3],
        "mixed": [0.7, 1.6, 3],
        "scalar": 2,
        "xy_z": [2, 3],
    }
    cases = []
    for fid, rf in factors.items():
        for m in ("linear", "cubic", "nearest"):
            cases.append(
                frame_case(f"frame.uint16.{fid}.{m}", u16, "im.tif", "o.tif", rf, "interpMethod", m)
            )
    for fid, rf in (
        ("down2", [2, 2, 2]),
        ("z_only", [1, 1, 2]),
        ("down3", [3, 3, 3]),
        ("blocks", [3, 5, 5]),
        ("mixed", [2, 3, 2]),
        ("scalar", 4),
    ):
        for m in ("max", "mean"):
            cases.append(
                frame_case(f"frame.uint16.{fid}.{m}", u16, "im.tif", "o.tif", rf, "interpMethod", m)
            )
    no16 = ["save16bit", False]
    cases += [
        frame_case(
            "frame.uint16.non_integer.max",
            u16,
            "im.tif",
            "o.tif",
            [1.5, 2.5, 2],
            "interpMethod",
            "max",
        ),
        frame_case(
            "frame.uint16.non_integer.mean",
            u16,
            "im.tif",
            "o.tif",
            [1.5, 1.5, 1.5],
            "interpMethod",
            "mean",
            error=True,
        ),
        # (sizes 39 x 36 x 12 divide by 1.5: reshape's error; 41 pads by 1: reshape's too)
        frame_case(
            "frame.uint16.non_integer.mean.divisible",
            tif_file((12, 39, 36), "uint16", "f/u16_39"),
            "im.tif",
            "o.tif",
            [1.5, 1.5, 1.5],
            "interpMethod",
            "mean",
            error=True,
        ),
        frame_case(
            "frame.uint16.non_integer.mean.integer_padding",
            tif_file((12, 41, 36), "uint16", "f/u16_41"),
            "im.tif",
            "o.tif",
            [1.5, 1.5, 1.5],
            "interpMethod",
            "mean",
            error=True,
        ),
        frame_case("frame.uint8.down2.linear", u8, "im.tif", "o.tif", [2, 2, 2]),
        frame_case(
            "frame.uint8.mixed.cubic", u8, "im.tif", "o.tif", [0.7, 1.6, 3], "interpMethod", "cubic"
        ),
        frame_case("frame.single.down2.linear", sgl, "im.tif", "o.tif", [2, 2, 2], *no16),
        frame_case("frame.single.down2.linear.one_nan", sgl1, "im.tif", "o.tif", [2, 2, 2], *no16),
        frame_case(
            "frame.single.non_integer.cubic.one_nan",
            sgl1,
            "im.tif",
            "o.tif",
            [1.5, 1.3, 2.5],
            "interpMethod",
            "cubic",
            *no16,
        ),
        frame_case(
            "frame.single.blocks.mean.one_nan",
            sgl1,
            "im.tif",
            "o.tif",
            [3, 5, 5],
            "interpMethod",
            "mean",
            *no16,
        ),
        frame_case(
            "frame.zarr.single.cubic.one_nan",
            zarr_tree(IM_SIZE, "single", "f/zs1", special="one"),
            "im.zarr",
            "o.zarr",
            [1.5, 2, 0.5],
            "zarrFile",
            True,
            "saveZarr",
            True,
            "interpMethod",
            "cubic",
            *no16,
        ),
        frame_case(
            "frame.single.non_integer.cubic",
            sgl,
            "im.tif",
            "o.tif",
            [1.5, 1.3, 2.5],
            "interpMethod",
            "cubic",
            *no16,
        ),
        frame_case(
            "frame.single.up.nearest",
            sgl,
            "im.tif",
            "o.tif",
            [0.5, 0.7, 1],
            "interpMethod",
            "nearest",
            *no16,
        ),
        frame_case(
            "frame.single.blocks.max",
            sgl,
            "im.tif",
            "o.tif",
            [3, 5, 5],
            "interpMethod",
            "max",
            *no16,
        ),
        frame_case(
            "frame.single.blocks.mean",
            sgl,
            "im.tif",
            "o.tif",
            [3, 5, 5],
            "interpMethod",
            "mean",
            *no16,
        ),
        frame_case("frame.single.down2.to_uint16", sgl, "im.tif", "o.tif", [2, 2, 2]),
        frame_case("frame.uint16.down2.single_out", u16, "im.tif", "o.tif", [2, 2, 2], *no16),
        frame_case(
            "frame.uint16.bbox",
            u16,
            "im.tif",
            "o.tif",
            [2, 2, 2],
            "inputBbox",
            [5, 7, 3, 30, 33, 10],
        ),
        frame_case(
            "frame.uint16.bbox.cubic",
            u16,
            "im.tif",
            "o.tif",
            [1.5, 0.8, 2],
            "inputBbox",
            [2, 3, 1, 37, 35, 11],
            "interpMethod",
            "cubic",
        ),
        # (one plane: the image path)
        frame_case(
            "frame.uint16.bbox.one_plane",
            u16,
            "im.tif",
            "o.tif",
            [2, 3, 2],
            "inputBbox",
            [5, 7, 6, 30, 33, 6],
        ),
        frame_case(
            "frame.zarr.zarr_out",
            zarr_tree(IM_SIZE, "uint16", "f/z"),
            "im.zarr",
            "o.zarr",
            [2, 2, 2],
            "zarrFile",
            True,
            "saveZarr",
            True,
        ),
        frame_case(
            "frame.zarr.tif_out.bbox",
            zarr_tree(IM_SIZE, "uint16", "f/z"),
            "im.zarr",
            "o.tif",
            [2, 2, 2],
            "zarrFile",
            True,
            "inputBbox",
            [5, 7, 3, 30, 33, 10],
        ),
        frame_case(
            "frame.zarr.single.cubic",
            zarr_tree(IM_SIZE, "single", "f/zs", special=True),
            "im.zarr",
            "o.zarr",
            [1.5, 2, 0.5],
            "zarrFile",
            True,
            "saveZarr",
            True,
            "interpMethod",
            "cubic",
            *no16,
        ),
        frame_case(
            "frame.tif.zarr_out.block",
            u16,
            "im.tif",
            "o.zarr",
            [2, 2, 2],
            "saveZarr",
            True,
            "blockSize",
            [8, 8, 4],
        ),
        frame_case("frame.image.down2", im2, "im.tif", "o.tif", [2, 2]),
        frame_case("frame.image.yxz_factor", im2, "im.tif", "o.tif", [1, 2, 4]),
        frame_case("frame.image.scalar.cubic", im2, "im.tif", "o.tif", 2, "interpMethod", "cubic"),
        frame_case(
            "frame.image.y_x.nearest", im2, "im.tif", "o.tif", [2, 3], "interpMethod", "nearest"
        ),
        frame_case("frame.image.up.linear", im2, "im.tif", "o.tif", [0.5, 0.7]),
        frame_case(
            "frame.image.max", im2, "im.tif", "o.tif", [2, 2], "interpMethod", "max", error=True
        ),
        frame_case(
            "frame.image.mean", im2, "im.tif", "o.tif", [2, 2], "interpMethod", "mean", error=True
        ),
        frame_case(
            "frame.uint16.unknown_method", u16, "im.tif", "o.tif", [2, 2, 2], "interpMethod", "box"
        ),
        frame_case(
            "frame.uint16.capitalized_method",
            u16,
            "im.tif",
            "o.tif",
            [2, 2, 2],
            "interpMethod",
            "Linear",
        ),
        # (the 10 decimals of the dataset wrapper give 23 x columns here, the full factor 22)
        frame_case("frame.uint16.full_precision", u16, "im.tif", "o.tif", [2, 1.60000000004, 2]),
        Case(
            "frame.existing_result",
            code="XR_resampleFrame([in{1}, 'im.tif'], [in{1}, 'o.tif'], [2, 2, 2]); out = {};",
            nargout=0,
            args=[InPath("d")],
            python=lambda d: call_matlab("XR_resampleFrame", d + "im.tif", d + "o.tif", [2, 2, 2]),
            inputs={"d": Tree(existing_tree, "d")},
        ),
    ]
    for m in ("linear", "cubic", "nearest", "max", "mean"):
        cases.append(
            frame_case(
                f"frame.constant.{m}", const_file, "im.tif", "o.tif", [2, 3, 2], "interpMethod", m
            )
        )
    return cases


# XR_resampleSingleZarr, then the result and the whole image resized at once
_ZARR_CODE = """d = in{1}; rf = in{2}; m = in{3};
XR_resampleSingleZarr([d, 'im.zarr'], [d, 'rs.zarr'], rf, 'interpMethod', m, in{4 : end});
im = readzarr([d, 'im.zarr']);
r = rf;
if isscalar(r)
    r = [r, r, r];
elseif numel(r) == 2
    r = [r(1), r(1), r(2)];
end
switch m
    case {'linear', 'cubic', 'nearest'}
        w = imresize3_fast(im, round(size(im) ./ r), m);
    case 'max'
        w = max_pooling_3d(im, r);
    otherwise
        w = imresize3_average(im, r);
end
out = {readzarr([d, 'rs.zarr']), w};"""


def _matlab_rf(rf):
    r = [float(v) for v in np.asarray(rf, dtype=np.float64).ravel()]
    if len(r) == 1:
        return r * 3
    if len(r) == 2:
        return [r[0], r[0], r[1]]
    return r


def whole_image(im, rf, method):
    """the whole image resized at once as the cases' MATLAB code does (rf: MATLAB's values)"""
    from petakit5d.kernels.mip import max_pooling_3d_matlab

    r = _matlab_rf(rf)
    rz = (r[2], r[0], r[1])
    if method in ("linear", "cubic", "nearest"):
        size = [float(v) for v in dtypes.mround(np.asarray(im.shape, np.float64) / np.asarray(rz))]
        return R.imresize3_fast(im, size, method)
    if method == "max":
        return max_pooling_3d_matlab(im, np.asarray(r))
    return R.imresize3_average(im, tuple(int(v) for v in rz))


def py_zarr(d, rf, method, *rest):
    from petakit5d import matlab

    matlab.XR_resampleSingleZarr(d + "im.zarr", d + "rs.zarr", rf, "interpMethod", method, *rest)
    return [read_zarr(d + "rs.zarr"), whole_image(read_zarr(d + "im.zarr"), rf, method)]


def zarr_case(cid, make, rf, method, *args, error=False, files=None):
    return Case(
        cid,
        code=_ZARR_CODE,
        nargout=2,
        args=[InPath("d"), rf, method, *args, "parseCluster", False],
        python=py_zarr,
        out_roles=None if error else ("vol", "vol"),
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
        files=files,
    )


def _zarr_cases():
    u16 = zarr_tree(IM_SIZE, "uint16", "z/u16")
    cases = []
    for m in ("linear", "cubic", "nearest", "max", "mean"):
        cases += [
            zarr_case(
                f"zarr.uint16.down2.{m}.border4", u16, [2, 2, 2], m, *SMALL, "borderSize", [4, 4, 4]
            ),
            zarr_case(f"zarr.uint16.down2.{m}.default_border", u16, [2, 2, 2], m, *SMALL),
        ]
    uneven = zarr_tree((23, 37, 29), "uint16", "z/un", (5, 7, 6))
    small_un = ["batchSize", [16, 12, 10], "blockSize", [8, 6, 5]]
    cases += [
        zarr_case(
            "zarr.uint16.down2.linear.no_border",
            u16,
            [2, 2, 2],
            "linear",
            *SMALL,
            "borderSize",
            [0, 0, 0],
        ),
        zarr_case(
            "zarr.uint16.down2.cubic.border1",
            u16,
            [2, 2, 2],
            "cubic",
            *SMALL,
            "borderSize",
            [1, 1, 1],
        ),
        zarr_case("zarr.uint16.down4.cubic", u16, [4, 4, 2], "cubic", *SMALL),
        zarr_case("zarr.uint16.down3.linear", u16, [3, 3, 3], "linear", *SMALL),
        zarr_case("zarr.uint16.non_integer.linear", u16, [1.5, 1.5, 1.5], "linear", *SMALL),
        zarr_case("zarr.uint16.non_integer.cubic", u16, [1.5, 1.3, 2.5], "cubic", *SMALL),
        zarr_case("zarr.uint16.digits.linear", u16, [1.2345678, 1.2345678, 2], "linear", *SMALL),
        zarr_case("zarr.uint16.up.linear", u16, [0.5, 0.5, 1], "linear", *SMALL),
        zarr_case("zarr.uint16.up.cubic", u16, [0.5, 0.7, 0.5], "cubic", *SMALL),
        zarr_case("zarr.uint16.mixed.cubic", u16, [0.7, 1.6, 3], "cubic", *SMALL),
        zarr_case("zarr.uint16.ones.linear", u16, [1, 1, 1], "linear", *SMALL),
        # (every factor 1: the batches are copied, whatever the method)
        zarr_case("zarr.uint16.ones.capitalized_method", u16, [1, 1, 1], "Linear", *SMALL),
        # (a batch size that is not whole chunks: rounded up to them)
        zarr_case(
            "zarr.uint16.batch_not_chunks.linear",
            u16,
            [1.5, 1.5, 1.5],
            "linear",
            "batchSize",
            [12, 12, 3],
            "blockSize",
            [8, 8, 2],
        ),
        zarr_case("zarr.uint16.scalar.linear", u16, 2, "linear", *SMALL),
        zarr_case("zarr.uint16.xy_z.cubic", u16, [2, 3], "cubic", *SMALL),
        zarr_case("zarr.uint16.one_batch.linear", u16, [2, 2, 2], "linear"),
        zarr_case("zarr.uint16.blocks.max", u16, [3, 5, 5], "max", *SMALL),
        zarr_case("zarr.uint16.blocks.mean", u16, [3, 5, 5], "mean", *SMALL),
        zarr_case(
            "zarr.uint16.bbox.linear",
            u16,
            [2, 2, 2],
            "linear",
            *SMALL,
            "inputBbox",
            [5, 7, 3, 30, 33, 10],
        ),
        zarr_case(
            "zarr.uint8.down2.linear",
            zarr_tree(IM_SIZE, "uint8", "z/u8"),
            [2, 2, 2],
            "linear",
            *SMALL,
        ),
        zarr_case(
            "zarr.single.down2.cubic",
            zarr_tree(IM_SIZE, "single", "z/sgl", special=True),
            [2, 2, 2],
            "cubic",
            *SMALL,
        ),
        zarr_case(
            "zarr.single.down2.mean",
            zarr_tree(IM_SIZE, "single", "z/sgl", special=True),
            [2, 2, 2],
            "mean",
            *SMALL,
        ),
        zarr_case(
            "zarr.single.down2.cubic.one_nan",
            zarr_tree(IM_SIZE, "single", "z/sgl1", special="one"),
            [2, 2, 2],
            "cubic",
            *SMALL,
        ),
        zarr_case(
            "zarr.single.down2.mean.one_nan",
            zarr_tree(IM_SIZE, "single", "z/sgl1", special="one"),
            [2, 2, 2],
            "mean",
            *SMALL,
        ),
        zarr_case("zarr.uint16.uneven.linear", uneven, [2, 3, 2], "linear", *small_un),
        zarr_case("zarr.uint16.uneven.cubic", uneven, [1.5, 2, 1.3], "cubic", *small_un),
        zarr_case("zarr.uint16.uneven.mean", uneven, [2, 3, 2], "mean", *small_un),
        # (the temporary results MATLAB and the port leave are not compared: their names don't
        # end with .zarr)
        zarr_case(
            "zarr.uint16.capitalized_method",
            u16,
            [2, 2, 2],
            "Linear",
            *SMALL,
            error=True,
            files=False,
        ),
    ]
    cases.append(
        Case(
            "zarr.missing_input",
            code=(
                "XR_resampleSingleZarr([in{1}, 'no.zarr'], [in{1}, 'rs.zarr'], [2, 2, 2], "
                "'parseCluster', false); out = {};"
            ),
            nargout=0,
            args=[InPath("d")],
            python=lambda d: call_matlab(
                "XR_resampleSingleZarr",
                d + "no.zarr",
                d + "rs.zarr",
                [2, 2, 2],
                "parseCluster",
                False,
            ),
            inputs={"d": Tree(u16, "d")},
            expect_error=True,
            py_error=ERRORS,
        )
    )
    return cases


def ds(cid, make, rf, *args, error=False):
    return Case(
        cid,
        code="XR_resample_dataset(in{:}); out = {};",
        nargout=0,
        args=[InPath("d"), rf, *args, *LOCAL],
        python=lambda *a: call_matlab("XR_resample_dataset", *a),
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
    )


def py_tif_names(*args):
    call_matlab("XR_resample_dataset", *args)
    return read_zarr(args[0] + "resampled/scan_CamA_ch0_t0000.tif")


def _dataset_cases():
    tif, zarr = dataset_tree("tif"), dataset_tree("zarr")
    large = ["zarrFile", True, "largeFile", True, "saveZarr", True, *SMALL]
    return [
        ds("ds.tif.down2", tif, [2, 2, 2], *CHANS),
        ds("ds.zarr.down2", zarr, [2, 2, 2], *CHANS, "zarrFile", True, "saveZarr", True),
        ds("ds.large.border4", zarr, [2, 2, 2], *ONE, *large, "borderSize", [4, 4, 4]),
        ds("ds.large.default_border", zarr, [2, 2, 2], *ONE, *large),
        ds("ds.large.digits", zarr, [1.23456789012, 2, 2], *ONE, *large),
        ds("ds.large.mean", zarr, [2, 3, 2], *ONE, *large, "interpMethod", "mean"),
        # (large files without saveZarr: zarr results named .tif, compared as arrays)
        Case(
            "ds.large.tif_names",
            code=(
                "XR_resample_dataset(in{:}); "
                "out = {readzarr([in{1}, 'resampled/scan_CamA_ch0_t0000.tif'])};"
            ),
            args=[
                InPath("d"),
                [2, 2, 2],
                *ONE,
                "zarrFile",
                True,
                "largeFile",
                True,
                *SMALL,
                *LOCAL,
            ],
            python=py_tif_names,
            out_roles=("vol",),
            inputs={"d": Tree(zarr, "d")},
            tree=tree(ignore=[r"scan_CamA_ch0_t0000\.tif"]),
        ),
        ds("ds.tif.ten_decimals", tif, [2, 1.60000000004, 2], *ONE),
        ds(
            "ds.tif.cubic.zarr_out",
            tif,
            [1.5, 2, 2],
            *ONE,
            "interpMethod",
            "cubic",
            "saveZarr",
            True,
            "blockSize",
            [8, 8, 4],
        ),
        ds("ds.tif.max", tif, [2, 3, 2], *ONE, "interpMethod", "max"),
        ds("ds.tif.mean.single", tif, [2, 3, 2], *ONE, "interpMethod", "mean", "save16bit", False),
        ds(
            "ds.tif.bbox",
            tif,
            [2, 2, 2],
            *ONE,
            "inputBbox",
            [5, 7, 3, 30, 33, 10],
            "resultDirName",
            "rs_box",
        ),
        ds("ds.tif.xy_z", dataset_tree("tif", "uint8"), [2, 3], *ONE),
        ds("ds.tif.no_match", tif, [2, 2, 2], "channelPatterns", ["CamC_ch9"], error=True),
    ]


def _block_cases():
    def tree_with_flag(flag_dir, flag):
        def make(folder):
            zarr_tree(IM_SIZE, "uint16", "b")(folder)
            zarr_file((6, 20, 18), "uint16", "b_out", (2, 8, 8))(Path(folder) / "rs.zarr")
            if flag_dir:
                (Path(folder) / "flags").mkdir()
                if flag:
                    (Path(folder) / "flags" / "blocks_1_2.mat").write_bytes(b"")

        return make

    code = (
        "resampleZarrBlock([1, 2], [in{1}, 'im.zarr'], [in{1}, 'rs.zarr'], "
        "[in{1}, 'flags/blocks_1_2.mat'], [2, 2, 2], 'batchSize', [16, 16, 4], "
        "'blockSize', [8, 8, 2]); out = {};"
    )

    def py(d):
        return call_matlab(
            "resampleZarrBlock",
            [1, 2],
            d + "im.zarr",
            d + "rs.zarr",
            d + "flags/blocks_1_2.mat",
            [2, 2, 2],
            "batchSize",
            [16, 16, 4],
            "blockSize",
            [8, 8, 2],
        )

    return [
        Case(
            "block.flag_exists",
            code=code,
            nargout=0,
            args=[InPath("d")],
            python=py,
            inputs={"d": Tree(tree_with_flag(True, True), "d")},
        ),
        Case(
            "block.no_flag_folder",
            code=code,
            nargout=0,
            args=[InPath("d")],
            python=py,
            inputs={"d": Tree(tree_with_flag(False, False), "d")},
            expect_error=True,
            py_error=ERRORS,
        ),
    ]


def _setting_cases():
    a, px, dz = 32.45, 0.108, 0.3
    cases = []
    for cid, args, err in (
        ("given.scalar", ["given", 2, False, a, px, dz], False),
        ("given.two", ["given", [2, 3], False, a, px, dz], False),
        ("given.three", ["given", [1, 2, 3], True, a, px, dz], False),
        ("given.empty", ["given", None, False, a, px, dz], True),
        ("isotropic", ["isotropic", None, False, a, px, dz], False),
        ("xy_isotropic", ["xy_isotropic", None, False, a, px, dz], False),
        ("xy_isotropic.objective_scan", ["xy_isotropic", None, True, a, px, 0.5], False),
        ("xy_isotropic.other_angle", ["xy_isotropic", [9], False, 31.5, 0.1, 0.2], False),
        ("unknown_type", ["other", [4, 5], False, a, px, dz], False),
    ):
        cases.append(
            Case(
                f"setting.{cid}",
                matlab="XR_checkResampleSetting",
                nargout=2,
                args=args,
                python=rs.check_resample_setting,
                expect_error=err,
                py_error=ERRORS,
            )
        )
    return cases


CASES = _frame_cases() + _zarr_cases() + _dataset_cases() + _block_cases() + _setting_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    # (the cases whose inputs hold two NaN payloads check that the first NaN in MATLAB's order
    # is kept where they meet, as in matlabcompat.resize)
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# ResampleTest in intent, on (z, y, x) data and native zarr arrays
# ---------------------------------------------------------------------------------------------

PY_FACTORS = {
    "down2": (2, 2, 2),
    "z_only": (2, 1, 1),
    "up2xy": (1, 0.5, 0.5),
    "non_integer": (1.5, 1.5, 1.5),
    "down3": (3, 3, 3),
}


def write_image(folder, fmt="tif", cls="uint16", name="im", shape=IM_SIZE):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    im = volume(shape, cls, f"py/{name}")
    if fmt == "zarr":
        fn = folder / f"{name}.zarr"
        write_zarr(im, str(fn), block_size=(4, 16, 16))
    else:
        fn = folder / f"{name}.tif"
        write_tiff(im, str(fn))
    return im, str(fn)


def write_dataset(folder, fmt):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ims = []
    for k, n in enumerate(DATASET):
        im = volume(IM_SIZE, "uint16", f"pyds/{fmt}/{k}")
        ims.append(im)
        if fmt == "zarr":
            write_zarr(im, str(folder / f"{n}.zarr"), block_size=(4, 16, 16))
        else:
            write_tiff(im, str(folder / f"{n}.tif"))
    return str(folder) + "/", ims


def round_size(shape, rf):
    return tuple(int(v) for v in dtypes.mround(np.asarray(shape, np.float64) / np.asarray(rf)))


def linear_reference(im, rf):
    """MATLAB's uint16(imresize3(single(im), round(size ./ rf), 'Method', 'linear'))"""
    return np.asarray(
        dtypes.uint16(R.imresize3(dtypes.single(im), round_size(im.shape, rf), "linear"))
    )


@pytest.mark.parametrize("fid", list(PY_FACTORS))
def test_linear_matches_imresize3(tmp_path, fid):
    im, fn = write_image(tmp_path)
    rs.resample_file(fn, str(tmp_path / "out.tif"), PY_FACTORS[fid])
    same_bits(read_tiff(str(tmp_path / "out.tif")), linear_reference(im, PY_FACTORS[fid]))


@pytest.mark.parametrize("method", ["linear", "cubic", "nearest", "max", "mean"])
def test_constant_image_stays_constant(tmp_path, method):
    write_tiff(np.full(IM_SIZE, 1000, dtype=np.uint16), str(tmp_path / "im.tif"))
    rf = (2, 2, 3)
    rs.resample_file(str(tmp_path / "im.tif"), str(tmp_path / "out.tif"), rf, interp_method=method)
    out = read_tiff(str(tmp_path / "out.tif"))
    if method in ("max", "mean"):
        size = tuple(-(-n // r) for n, r in zip(IM_SIZE, rf, strict=True))
    else:
        size = round_size(IM_SIZE, rf)
    assert out.shape == size
    assert np.unique(out).tolist() == [1000]


def block_stats(im, rf):
    """ResampleTest's blockStats: max and mean of each block of the voxels inside"""
    out = tuple(-(-n // r) for n, r in zip(im.shape, rf, strict=True))
    bmax = np.zeros(out, np.float32)
    bmean = np.zeros(out, np.float32)
    for k in range(out[0]):
        for i in range(out[1]):
            for j in range(out[2]):
                b = im[
                    k * rf[0] : (k + 1) * rf[0],
                    i * rf[1] : (i + 1) * rf[1],
                    j * rf[2] : (j + 1) * rf[2],
                ]
                bmax[k, i, j] = b.max()
                bmean[k, i, j] = b.astype(np.float64).mean()
    return bmax, bmean


def test_max_and_mean_match_block_reference(tmp_path):
    # non-divisible size: edge blocks use only the voxels inside
    im, fn = write_image(tmp_path)
    rf = (5, 3, 5)
    rs.resample_file(fn, str(tmp_path / "max.tif"), rf, interp_method="max", save_16bit=False)
    rs.resample_file(fn, str(tmp_path / "mean.tif"), rf, interp_method="mean", save_16bit=False)
    bmax, bmean = block_stats(im.astype(np.float32), rf)
    same_bits(read_tiff(str(tmp_path / "max.tif")), bmax)
    np.testing.assert_allclose(read_tiff(str(tmp_path / "mean.tif")), bmean, rtol=1e-6)


def test_image_2d(tmp_path):
    im = volume((40, 36), "uint16", "py/2d")
    write_tiff(im, str(tmp_path / "im.tif"))
    # (an image's two values are MATLAB's [y x]: here (z, xy) = (2, 2))
    rs.resample_file(str(tmp_path / "im.tif"), str(tmp_path / "out.tif"), (2, 2))
    expected = np.asarray(dtypes.uint16(R.imresize(dtypes.single(im), (20, 18), "bilinear")))
    same_bits(read_tiff(str(tmp_path / "out.tif")), expected)


def test_image_2d_takes_matlab_y_x_values(tmp_path):
    # MATLAB's [y x] = [2 3] reaches the port as (3, 2): y by 2, x by 3
    im = volume((40, 36), "uint16", "py/2dyx")
    write_tiff(im, str(tmp_path / "im.tif"))
    call_matlab("XR_resampleFrame", str(tmp_path / "im.tif"), str(tmp_path / "out.tif"), [2, 3])
    assert read_tiff(str(tmp_path / "out.tif")).shape == (20, 12)
    rs.resample_file(str(tmp_path / "im.tif"), str(tmp_path / "out2.tif"), (3, 2))
    same_bits(read_tiff(str(tmp_path / "out2.tif")), read_tiff(str(tmp_path / "out.tif")))


def test_input_bbox_resamples_the_cropped_region(tmp_path):
    im, fn = write_image(tmp_path)
    b = (2, 4, 6, 10, 30, 33)
    rs.resample_file(fn, str(tmp_path / "out.tif"), (2, 2, 2), input_bbox=b)
    region = im[b[0] : b[3], b[1] : b[4], b[2] : b[5]]
    same_bits(read_tiff(str(tmp_path / "out.tif")), linear_reference(region, (2, 2, 2)))


def test_single_output_without_save_16bit(tmp_path):
    im, fn = write_image(tmp_path)
    rs.resample_file(fn, str(tmp_path / "out.tif"), (2, 2, 2), save_16bit=False)
    expected = R.imresize3(dtypes.single(im), round_size(IM_SIZE, (2, 2, 2)), "linear")
    same_bits(read_tiff(str(tmp_path / "out.tif")), expected)


def test_zarr_input_and_output(tmp_path):
    im, fn = write_image(tmp_path, "zarr")
    out = str(tmp_path / "out.zarr")
    rs.resample_file(fn, out, (2, 2, 2), zarr_file=True, save_zarr=True)
    same_bits(read_zarr(out), linear_reference(im, (2, 2, 2)))


def test_image_2d_with_yxz_factor(tmp_path):
    # a 2-D image takes the y and x factors of a (z, y, x) factor
    write_tiff(volume((40, 36), "uint16", "py/2dz"), str(tmp_path / "im.tif"))
    rs.resample_file(str(tmp_path / "im.tif"), str(tmp_path / "out.tif"), (4, 1, 2))
    assert read_tiff(str(tmp_path / "out.tif")).shape == (40, 18)


def test_dataset_wrapper_tiff(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "tif")
    rs.resample_dataset(
        path, (2, 2, 2), channel_patterns=("CamA_ch0", "CamB_ch0"), parse_cluster=False
    )
    result = Path(path) / "resampled"
    for k in range(2):
        same_bits(read_tiff(str(result / f"{DATASET[k]}.tif")), linear_reference(ims[k], (2, 2, 2)))
    assert not (result / f"{DATASET[2]}.tif").exists(), (
        "a file that matches no pattern was resampled"
    )
    import json

    params = json.loads((result / "parameters.json").read_text())
    assert params["resampleFactor"] == [2, 2, 2]


def test_dataset_wrapper_zarr(tmp_path):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    rs.resample_dataset(
        path,
        (2, 2, 2),
        channel_patterns=("CamA_ch0", "CamB_ch0"),
        zarr_file=True,
        save_zarr=True,
        parse_cluster=False,
    )
    for k in range(2):
        out = read_zarr(f"{path}resampled/{DATASET[k]}.zarr")
        same_bits(out, linear_reference(ims[k], (2, 2, 2)))


def _large(tmp_path, border):
    path, ims = write_dataset(tmp_path / "d", "zarr")
    kw = {} if border is None else {"border_size": border}
    rs.resample_dataset(
        path,
        (2, 2, 2),
        channel_patterns=("CamA_ch0",),
        zarr_file=True,
        large_file=True,
        save_zarr=True,
        batch_size=(4, 16, 16),
        block_size=(2, 8, 8),
        parse_cluster=False,
        **kw,
    )
    return read_zarr(f"{path}resampled/{DATASET[0]}.zarr"), ims[0]


@pytest.mark.parametrize("border", [(4, 4, 4), None], ids=["border4", "default_border"])
def test_large_zarr_matches_in_memory(tmp_path, border):
    # batch-wise resampling with border buffers vs the whole image (largeZarrMatchesInMemory,
    # largeZarrWithDefaultBorderSize): the batches equal the whole image resized in uint16,
    # and differ from the single resize of the file path by rounding only (+/-1)
    out, im = _large(tmp_path, border)
    same_bits(out, R.imresize3_fast(im, round_size(IM_SIZE, (2, 2, 2)), "linear"))
    ref = linear_reference(im, (2, 2, 2))
    assert np.abs(out.astype(np.int64) - ref).max() <= 1


@pytest.fixture
def framework_spy(monkeypatch):
    """runs the tasks as run_tasks does locally, after ``spy.action`` (MATLAB's
    FrameworkSpyFixture); records the calls"""
    from petakit5d.jobs import framework

    real = framework.run_tasks
    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        if fake.action is not None:
            fake.action()
        if not fake.run:
            return [True] * len(list(outputs))
        return real(inputs, outputs, tasks, **kw)

    fake.action = None
    fake.run = True
    fake.calls = calls
    monkeypatch.setattr(framework, "run_tasks", fake)
    return fake


def test_large_zarr_keeps_result_of_another_copy(tmp_path, framework_spy):
    # two copies of a task share the temporary result: the other moves it into place while
    # this one runs
    _, fn = write_image(tmp_path, "zarr")
    out = str(tmp_path / "resampled.zarr")
    framework_spy.action = lambda: os.replace(out + "_test", out)
    framework_spy.run = False
    rs.resample_zarr(
        fn,
        out,
        (2, 2, 2),
        batch_size=(4, 16, 16),
        block_size=(2, 8, 8),
        parse_cluster=False,
        uuid="test",
    )
    assert os.path.isdir(out), "the result of the other copy was deleted"


def test_large_zarr_keeps_result_of_another_run(tmp_path, framework_spy):
    # a copy from another run (its own uuid) finds the result published: it keeps it
    im, fn = write_image(tmp_path, "zarr")
    out = str(tmp_path / "resampled.zarr")
    published = np.ascontiguousarray(im[::2, ::2, ::2])
    framework_spy.action = lambda: write_zarr(published, out)
    rs.resample_zarr(
        fn,
        out,
        (2, 2, 2),
        batch_size=(4, 16, 16),
        block_size=(2, 8, 8),
        parse_cluster=False,
        uuid="test",
    )
    same_bits(read_zarr(out), published)
    assert not os.path.isdir(out + "_test"), "the temporary result was left behind"


# ---------------------------------------------------------------------------------------------
# the batches and the factors at the task boundaries
# ---------------------------------------------------------------------------------------------


def test_batch_place_follows_matlab_arithmetic():
    # z batch 2 (0-based 1) of 4 output planes, factor 2, border 6: input planes 9..16 (1-based)
    # plus the border, the region from batch plane 4 (1-based) on
    p = rs.batch_place((1, 0, 0), (4, 16, 16), (2, 2, 2), (6, 6, 6), (12, 40, 36), (6, 20, 18))
    assert p["in_start_orig"] == [9, 1, 1]
    assert p["in_start"] == [3, 1, 1] and p["in_end"] == [12, 38, 36]
    assert p["out_start"] == [5, 1, 1] and p["out_end"] == [6, 16, 16]
    assert p["ba_start"] == [4, 1, 1] and p["ba_end"] == [5, 16, 16]
    # with a box, the end is cut at the image's end before the box's start is added
    q = rs.batch_place(
        (1, 0, 0), (4, 16, 16), (2, 2, 2), (0, 0, 0), (12, 40, 36), (4, 13, 14), (2, 4, 6)
    )
    assert q["in_start_orig"] == [11, 5, 7] and q["in_end"] == [12, 36, 36]


def test_block_order_is_matlab_ind2sub():
    # y fastest, then x, then z
    n = (2, 3, 4)
    subs = [rs._batch_sub(b, n) for b in range(24)]
    assert subs[:5] == [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 0, 1), (0, 1, 1)]
    assert subs[12] == (1, 0, 0)


def test_factors_reach_the_tasks_as_matlab_text(tmp_path, framework_spy):
    framework_spy.run = False
    path, _ = write_dataset(tmp_path / "d", "zarr")
    rs.resample_dataset(
        path,
        (2, 2.000000000049, 1 / 3),
        channel_patterns=("CamA_ch0",),
        zarr_file=True,
        parse_cluster=False,
    )
    # mat2str_comma(resampleFactor, 10)
    assert framework_spy.calls[-1]["tasks"][0].args[2] == [2.0, 2.0, 0.3333333333]
    fn = f"{path}{DATASET[0]}.zarr"
    rs.resample_zarr(
        fn,
        str(tmp_path / "rs.zarr"),
        (1.2345678, 2, 1 / 3),
        batch_size=(4, 16, 16),
        block_size=(2, 8, 8),
        parse_cluster=False,
    )
    # num2str(resampleFactor, '%d,'): %e of a fraction
    t = framework_spy.calls[-1]["tasks"][0]
    assert t.args[4] == [1.234568, 2.0, 0.3333333]
    assert t.kwargs["block_size"] == [2, 8, 8] and t.kwargs["border_size"] == [5.0, 5.0, 5.0]


def test_matlab_layer_converts_factors(tmp_path, framework_spy):
    framework_spy.run = False
    _, fn = write_image(tmp_path, "zarr")
    # [y x z] = [3 2 4], and [xy z] = [3 4]
    call_matlab(
        "XR_resampleSingleZarr", fn, str(tmp_path / "a.zarr"), [3, 2, 4], "parseCluster", False
    )
    assert framework_spy.calls[-1]["tasks"][0].args[4] == [4.0, 3.0, 2.0]
    call_matlab(
        "XR_resampleSingleZarr", fn, str(tmp_path / "b.zarr"), [3, 4], "parseCluster", False
    )
    assert framework_spy.calls[-1]["tasks"][0].args[4] == [4.0, 3.0, 3.0]


def test_unknown_method_saves_the_image_unresampled(tmp_path):
    im, fn = write_image(tmp_path)
    rs.resample_file(
        fn, str(tmp_path / "out.tif"), (2, 2, 2), interp_method="Linear", save_16bit=False
    )
    same_bits(read_tiff(str(tmp_path / "out.tif")), im.astype(np.float32))


def test_block_task_fails_for_a_capitalized_method(tmp_path):
    _, fn = write_image(tmp_path, "zarr")
    # (locally the framework raises the task's error)
    with pytest.raises(rs.ResampleError, match="out_batch|tasks failed"):
        rs.resample_zarr(
            fn,
            str(tmp_path / "rs.zarr"),
            (2, 2, 2),
            batch_size=(4, 16, 16),
            block_size=(2, 8, 8),
            interp_method="Linear",
            parse_cluster=False,
        )


def test_block_skips_done_batches_and_checks_its_folders(tmp_path, capsys):
    _, fn = write_image(tmp_path, "zarr")
    with pytest.raises(rs.ResampleError, match="does not exist"):
        rs.resample_zarr_block(
            [0], fn, fn, str(tmp_path / "no" / "f.mat"), (2, 2, 2), block_size=(2, 8, 8)
        )
    (tmp_path / "flags").mkdir()
    flag = tmp_path / "flags" / "blocks_1_1.mat"
    flag.write_bytes(b"")
    assert rs.resample_zarr_block([0], fn, fn, str(flag), (2, 2, 2), block_size=(2, 8, 8)) is None
    assert "already exist, skip them" in capsys.readouterr().out


@pytest.mark.parametrize("layout", ["native", "matlab"])
def test_native_and_matlab_layout_inputs_give_the_same_result(tmp_path, layout):
    im = volume(IM_SIZE, "uint16", "py/layout")
    fn = str(tmp_path / "im.zarr")
    synth.write_zarr(fn, im, block_size=(4, 16, 16), layout=layout)
    out = str(tmp_path / "rs.zarr")
    rs.resample_zarr(
        fn, out, (2, 1.5, 2), batch_size=(4, 16, 16), block_size=(2, 8, 8), parse_cluster=False
    )
    ref = str(tmp_path / "ref.zarr")
    ref_in = str(tmp_path / "ref_in.zarr")
    synth.write_zarr(ref_in, im, block_size=(4, 16, 16), layout="native")
    rs.resample_zarr(
        ref_in, ref, (2, 1.5, 2), batch_size=(4, 16, 16), block_size=(2, 8, 8), parse_cluster=False
    )
    same_bits(read_zarr(out), read_zarr(ref))


def test_imresize3_functions_are_matlabcompat_resize():
    assert rs.imresize3_fast is R.imresize3_fast
    assert rs.imresize3_average is R.imresize3_average


# (z, y, x) factors and borders; whether MATLAB's batches equal the whole image resized at once
# (the zarr.* oracle cases give both)
BORDER_CASES = [
    ("linear", (2, 2, 2), (5, 5, 5), True),
    ("cubic", (2, 2, 2), (4, 4, 4), True),
    ("cubic", (2, 4, 4), (5, 5, 5), True),
    ("linear", (3, 3, 3), (5, 5, 5), True),
    ("linear", (1, 0.5, 0.5), (5, 5, 5), True),
    ("nearest", (2, 2, 2), (0, 0, 0), True),
    ("linear", (2, 2, 2), (0, 0, 0), False),
    ("cubic", (2, 2, 2), (1, 1, 1), False),
    ("linear", (1.5, 1.5, 1.5), (5, 5, 5), False),
    ("cubic", (0.5, 0.5, 0.7), (5, 5, 5), False),
]


@pytest.mark.parametrize("method, rf, border, equal", BORDER_CASES)
def test_batches_equal_the_whole_image_where_matlab_s_do(tmp_path, method, rf, border, equal):
    im, fn = write_image(tmp_path, "zarr")
    out = str(tmp_path / "rs.zarr")
    rs.resample_zarr(
        fn,
        out,
        rf,
        batch_size=(4, 16, 16),
        block_size=(2, 8, 8),
        border_size=border,
        interp_method=method,
        parse_cluster=False,
    )
    whole = R.imresize3_fast(im, round_size(IM_SIZE, rf), method)
    assert np.array_equal(read_zarr(out), whole) == equal
