"""Channel unmixing (petakit5d.tools.unmix, task B02).

Oracle cases against MATLAB: unmix_channels_sum (the mex path) for every result class over
uint8, uint16, single and double channels and their mixes (single and double with NaN, Inf, -0,
huge values and values that round at .5), backgrounds and offsets (0, NaN, single), factors
that single can't represent, 3D, 2D, vectors, empty arrays and sizes around the mex's blocks of
8192; partial sums (no class) and sums added to a partial sum; the MATLAB code of the other
classes (int16 and logical channels, int32 and int8 results). Then XR_unmix_channels_frame
(tiff and zarr, linear and gaussian with isotropic and anisotropic sigmas, flat-field
correction, backgrounds, offsets, factor 0, 2 to 4 channels at once and one at a time, uint8,
uint16 and single channels, uint16 and single results, 2D frames, all factors 0),
XR_unmix_channels_zarr (linear and gaussian blocks, borders, batches that don't divide the
image, single and uint8 channels, 3 channels) and XR_unmix_channels_data_wrapper (the demo's
settings, the second channel, flat fields, zarr in and out, large files, gaussian sigmas
through the task strings, factors that the task strings round, one channel per folder, groups
without every channel, MATLAB's errors). Output tiffs are compared byte for byte, zarr results
by data and metadata after the layout mapping. Inputs are MATLAB's: tiffs and zarr arrays in
MATLAB's layout, read by both sides.

Then UnmixChannelsTest (tests/analysis) in intent on the port's (z, y, x) data and native zarr
arrays, the kernel against its NumPy reference, the text rounding of the task strings, the job
shapes (a framework spy), and results that don't depend on the threads or the path.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.io import read_tiff, read_zarr, write_tiff, write_zarr
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, Vol, Work, mcell
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import unmix

# ---------------------------------------------------------------------------------------------
# data: (z, y, x) arrays from counter-based random numbers (bit-stable everywhere)
# ---------------------------------------------------------------------------------------------

F32 = np.float32
NAN32 = F32(NAN)
EDGES = [NAN, np.inf, -np.inf, -0.0, 0.0, 0.49999997, 0.5, 65535.5, 254.5, 255.5, 1e-40, 3e38]
EDGES += [-3e38, 1e300]
DT = {"uint8": np.uint8, "uint16": np.uint16, "single": np.float32, "double": np.float64}
DT.update({"int16": np.int16, "int8": np.int8, "int32": np.int32, "logical": np.bool_})


def channel(shape, cls, name):
    """UnmixChannelsTest's randomChannel: uint8 and uint16 over their range, single and double
    (rand - 0.3) * 70000 with NaN, Inf, -0, huge values and values that round at .5"""
    rng = CounterRNG(2, f"b02/{name}/{cls}/{'x'.join(map(str, shape))}")
    shape = tuple(int(s) for s in shape)
    if cls == "uint8":
        return rng.integers(0, 255, shape, dtype=np.uint8)
    if cls == "uint16":
        return rng.integers(0, 65535, shape, dtype=np.uint16)
    if cls == "int16":
        return rng.integers(-100, 3000, shape, dtype=np.int16)
    if cls == "logical":
        return rng.integers(0, 1, shape, dtype=np.uint8).astype(bool)
    x = (rng.uniform(shape) - 0.3) * 70000.0
    flat = x.reshape(-1)
    n = min(flat.size, 3 * len(EDGES))
    if n:
        pos = np.argsort(rng.spawn("pos").uniform(flat.size), kind="stable")[:n]
        flat[pos] = [EDGES[k % len(EDGES)] for k in range(n)]
    with np.errstate(over="ignore"):
        return x.reshape(shape).astype(DT[cls])


def texture(shape, int_range, name):
    """pktest.texture as uint16 (rounded as MATLAB's uint16 rounds positive values)"""
    t = synth.texture(shape, 2.0, int_range, CounterRNG(5, f"b02/tex/{name}"))
    return np.floor(t.astype(np.float64) + 0.5).astype(np.uint16)


def same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype, (a.dtype, b.dtype)
    assert a.shape == b.shape, (a.shape, b.shape)
    assert a.tobytes() == b.tobytes()


# ---------------------------------------------------------------------------------------------
# oracle cases: unmix_channels_sum
# ---------------------------------------------------------------------------------------------

SIZES = [(23, 19, 11), (17, 13), (1, 9), (8193, 1), (8191, 2), (0, 4), (1, 1), (700, 600, 3)]
CLASS_SETS = [
    ("uint16", "uint16"),
    ("uint8", "single", "double"),
    ("single",),
    ("double", "uint16"),
    ("single", "uint8"),
]
FACTORS = [[1, -0.1], [0.1, -2.345678901234, 1 / 3], 1e3, np.array([-0.25, 1], F32), [1, 0.5]]
BACKGROUNDS = [[100, 99.7], [0, 99.7, -1 / 3], None, np.array([99.5, 0], F32), [NAN, 1 / 3]]
OFFSETS = [1, -3.3, None, F32(1.7), 0, NAN, 65535]


def py_shape(msize):
    """the Python shape of MATLAB's size [ny nx nz] or [ny nx]"""
    return tuple(msize) if len(msize) == 2 else (msize[2], msize[0], msize[1])


def ch_value(msize, cls, name):
    return Value(lambda: to_matlab(channel(py_shape(msize), cls, name)))


def sum_case(cid, ims, factors, bgs, offset, dtype, im0=None):
    inputs = {f"c{k}": v for k, v in enumerate(ims)}
    args = [[Vol(In(f"c{k}")) for k in range(len(ims))], factors, bgs, offset, dtype]
    if im0 is not None:
        inputs["im0"] = im0
        args.append(Vol(In("im0")))
    py_args = None
    if not ims:
        # (MATLAB's empty cell; the Python call gets an empty list)
        py_args = [[], *args[1:]]
        args[0] = mcell()
    return Case(
        cid,
        matlab="unmix_channels_sum",
        args=args,
        py_args=py_args,
        python=unmix.unmix_channels_sum,
        out_roles=("vol",),
        inputs=inputs,
    )


def _sum_cases():
    cases = []
    for i, msize in enumerate(SIZES):
        for k, classes in enumerate(CLASS_SETS):
            ims = [ch_value(msize, c, f"sum{i}/{k}/{j}") for j, c in enumerate(classes)]
            # (MATLAB's mod(i + k, numel(offsets)) + 1, with 1-based i and k)
            o = OFFSETS[(i + k + 2) % len(OFFSETS)]
            for out in ("uint8", "uint16", "single", "double"):
                sz = "x".join(map(str, msize))
                cid = f"sum.{out}.{sz}.{'_'.join(classes)}"
                cases.append(sum_case(cid, ims, FACTORS[k], BACKGROUNDS[k], o, out))
    # factors 0 and Inf, a background of 0 for every channel, a zero offset
    ims = [ch_value((23, 19, 11), c, f"inf/{c}") for c in ("uint16", "single")]
    cases.append(sum_case("sum.uint16.factor_inf", ims, [np.inf, 0], [0, 0], 0, "uint16"))
    cases.append(sum_case("sum.single.factor_inf", ims, [np.inf, 0], [0, 0], 0, "single"))
    cases.append(sum_case("sum.single.factor_nan", ims, [NAN, 1], [3, NAN], 2, "single"))
    return cases


def _partial_cases():
    msize = (23, 19, 11)
    classes = ("uint16", "single", "uint8", "double")
    ims = [ch_value(msize, c, f"partial/{c}") for c in classes]
    factors = [1, -0.1, 1 / 3, -2.345678901234]
    bgs = [100, 99.7, 0, 1 / 3]
    im0 = Value(lambda: to_matlab(channel(py_shape(msize), "single", "partial/im0")))
    cases = [
        sum_case("partial.all", ims, factors, bgs, None, ""),
        sum_case("partial.into_im0", ims[:2], factors[:2], bgs[:2], None, "", im0),
        sum_case("partial.one_channel_into_im0", ims[2:3], factors[2:3], bgs[2:3], None, "", im0),
    ]
    for out in ("uint8", "uint16", "single", "double"):
        cases.append(sum_case(f"partial.finish.{out}", [], [], [], 1.7, out, im0))
        cases.append(sum_case(f"partial.finish_no_offset.{out}", [], [], [], None, out, im0))
    cases.append(sum_case("partial.finish_nan_offset", [], [], [], NAN, "single", im0))
    return cases


def _other_class_cases():
    msize = (9, 7, 3)
    x = ch_value(msize, "int16", "other/x")
    y = ch_value(msize, "uint16", "other/y")
    lg = ch_value(msize, "logical", "other/l")
    return [
        sum_case("other.int16_channel", [x, y], [1, -0.1], [10, 5], 2, "uint16"),
        sum_case("other.int32_result", [y, y], [1, -0.1], None, 2, "int32"),
        sum_case("other.int8_result", [y, x], [0.01, -0.1], None, -3.5, "int8"),
        sum_case("other.logical_channel", [lg, y], [100, 0.5], [0.25, 1], None, "single"),
        sum_case(
            "other.single_factor_offset", [y, x], np.array([1, -0.1], F32), None, F32(2.5), "uint16"
        ),
    ]


# ---------------------------------------------------------------------------------------------
# oracle cases: frames, zarr files and datasets (MATLAB's layout on disk)
# ---------------------------------------------------------------------------------------------

FSN = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs"
    "_000x_000y_000z_0000t"
)
SHAPE = (24, 48, 40)  # UnmixChannelsTest's [48 40 24]
LOCAL = ["parseCluster", False, "masterCompute", True]
ERRORS = (unmix.UnmixError, ValueError, OSError, IndexError, RuntimeError)


def channels_abc(shape=SHAPE, tag=""):
    """twoChannelData's CamA and CamB (CamB bleeds into CamA), and frameMatchesTheOriginalCode's
    CamC"""
    b = texture(shape, (100, 3000), f"{tag}B")
    a_t = synth.texture(shape, 2.0, (100, 1500), CounterRNG(5, f"b02/tex/{tag}A"))
    a = a_t.astype(np.float64) + 0.1 * (b.astype(np.float32).astype(np.float64) - 100)
    a = np.floor(np.maximum(a, 0) + 0.5).astype(np.uint16)
    c = texture(shape, (0, 900), f"{tag}C")
    return {"CamA": a, "CamB": b, "CamC": c}


def ff_images(shape=SHAPE[1:]):
    """frameMatchesTheOriginalCode's flat fields (300 + 700 x / 40 + 100 c) and backgrounds
    (90 + 5 c)"""
    x = np.arange(1, shape[1] + 1, dtype=np.float64)[np.newaxis, :].repeat(shape[0], 0)
    ffs, bks = [], []
    for c in range(1, 4):
        ffs.append(np.floor(300 + 700 * x / 40 + 100 * c + 0.5).astype(np.uint16))
        bks.append(np.full(shape, 90 + 5 * c, dtype=np.uint16))
    return ffs, bks


def frame_tree(cams=("CamA", "CamB", "CamC"), cls="uint16", fmt="tif", shape=SHAPE, tag=""):
    """a folder with one file per channel (zarr in MATLAB's layout) and the flat fields"""

    def make(folder):
        chans = channels_abc(shape, tag)
        ffs, bks = ff_images(shape[-2:])
        for c, cam in enumerate(cams):
            a = chans[cam]
            if cls == "single":
                a = (a.astype(F32) + F32(0.37)).astype(F32)
            elif cls == "uint8":
                a = (a // 16).astype(np.uint8)
            name = FSN.format(cam)
            if fmt == "zarr":
                synth.write_zarr(Path(folder) / f"{name}.zarr", a, (8, 16, 16), layout="matlab")
            else:
                synth.write_tiff(Path(folder) / f"{name}.tif", a)
            synth.write_tiff(Path(folder) / f"ff_{c + 1}.tif", ffs[c])
            synth.write_tiff(Path(folder) / f"bk_{c + 1}.tif", bks[c])

    return make


def files(cams, fmt="tif"):
    return [Work(f"d/{FSN.format(c)}.{fmt}") for c in cams]


FF = [Work(f"d/ff_{c}.tif") for c in (1, 2, 3)]
BK = [Work(f"d/bk_{c}.tif") for c in (1, 2, 3)]


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def py_frame(*args):
    return call_matlab("XR_unmix_channels_frame", *args)


def py_zarr(*args):
    return call_matlab("XR_unmix_channels_zarr", *args)


def py_dataset(*args):
    return call_matlab("XR_unmix_channels_data_wrapper", *args)


def frame_case(cid, make, cams, factors, *args, fmt="tif", error=False):
    return Case(
        cid,
        matlab="XR_unmix_channels_frame",
        nargout=0,
        args=[files(cams, fmt), factors, "uuid", "t", *args],
        python=py_frame,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
    )


SIG3 = np.array([[0, 0, 0], [1, 1, 1], [0.7, 0.7, 2]], dtype=np.float64)


def _frame_cases():
    abc = ("CamA", "CamB", "CamC")
    t3 = frame_tree()
    ff = ["FFImagePaths", FF, "backgroundPaths", BK]
    bg3 = ["constBackground", [100, 99.5, 98]]
    return [
        # frameMatchesTheOriginalCode's settings
        frame_case("frame.linear_bg_offset", t3, abc, [1, -0.1, 0], *bg3, "constOffset", 1),
        frame_case(
            "frame.gaussian_bg_offset_single",
            t3,
            abc,
            [1, -0.1, 0.37],
            "mode",
            "gaussian",
            "unmixSigmas",
            SIG3,
            *bg3,
            "constOffset",
            0.5,
            "save16bit",
            False,
        ),
        frame_case(
            "frame.ff_negative_offset",
            t3,
            abc,
            [-0.25, 1, -1 / 3],
            "FFCorrection",
            True,
            *ff,
            "constOffset",
            -0.3,
        ),
        frame_case(
            "frame.gaussian_ff_bg_single",
            t3,
            abc,
            [1, 0, -0.1234567],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0, 0, 0], [1, 1, 1], [0, 0, 0]], dtype=np.float64),
            "FFCorrection",
            True,
            *ff,
            "constBackground",
            [3, 2, 1],
            "save16bit",
            False,
        ),
        frame_case(
            "frame.linear_bg_negative_offset_single",
            t3,
            abc,
            [1, -0.1, 0.37],
            *bg3,
            "constOffset",
            -0.3,
            "save16bit",
            False,
        ),
        frame_case("frame.linear_three_one_at_a_time", t3, abc, [1, -0.1, 0.37]),
        # single channels (backgrounds subtracted in place), uint16 and single results
        frame_case(
            "frame.single_channels_uint16",
            frame_tree(cls="single", tag="s"),
            abc,
            [1, -0.1, 0.37],
            "constBackground",
            [100, 99.7, 98],
            "constOffset",
            1,
        ),
        frame_case(
            "frame.single_channels_single",
            frame_tree(cls="single", tag="s"),
            abc,
            [1, -0.1, 0.37],
            "constBackground",
            [100, 99.7, 98],
            "constOffset",
            1,
            "save16bit",
            False,
        ),
        # two uint16 channels at once; the second channel unmixed; uint8 channels at once
        frame_case(
            "frame.two_at_once",
            frame_tree(("CamA", "CamB")),
            ("CamA", "CamB"),
            [1, -0.1],
            "constBackground",
            [100, 100],
            "constOffset",
            1,
        ),
        frame_case(
            "frame.second_channel",
            frame_tree(("CamA", "CamB")),
            ("CamA", "CamB"),
            [-0.25, 1],
            "channelInd",
            2,
            "save16bit",
            False,
        ),
        frame_case(
            "frame.uint8_channels",
            frame_tree(cls="uint8", tag="u8"),
            abc,
            [1, -0.5, 0.25],
            "constBackground",
            [3, 2.5, 1],
        ),
        # four channels (one at a time), all factors 0, one channel
        frame_case(
            "frame.all_factors_zero",
            frame_tree(("CamA", "CamB")),
            ("CamA", "CamB"),
            [0, 0],
            "constOffset",
            2.5,
        ),
        frame_case(
            "frame.one_channel_gaussian_aniso",
            frame_tree(("CamA",)),
            ("CamA",),
            [0.5],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[1.5, 0.8, 2.5]]),
        ),
        frame_case(
            "frame.gaussian_scalar_sigmas",
            frame_tree(("CamA", "CamB")),
            ("CamA", "CamB"),
            [1, -0.2],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0.0], [1.2]]),
        ),
        # zarr in and out, 2D frames
        frame_case(
            "frame.zarr_in_out",
            frame_tree(("CamA", "CamB"), fmt="zarr", tag="z"),
            ("CamA", "CamB"),
            [1, -0.1],
            "constBackground",
            [100, 100],
            "zarrFile",
            True,
            "saveZarr",
            True,
            "blockSize",
            [16, 8, 8],
            fmt="zarr",
        ),
        frame_case(
            "frame.tiff_to_zarr",
            frame_tree(("CamA", "CamB"), tag="tz"),
            ("CamA", "CamB"),
            [1, -0.1],
            "saveZarr",
            True,
        ),
        frame_case(
            "frame.two_d",
            frame_tree(("CamA", "CamB"), shape=(1, 48, 40), tag="2d"),
            ("CamA", "CamB"),
            [1, -0.1],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0, 0, 0], [1, 1, 1]], dtype=np.float64),
        ),
        frame_case(
            "frame.existing_result",
            _existing_frame_tree(),
            ("CamA", "CamB"),
            [1, -0.1],
        ),
    ]


def _existing_frame_tree():
    def make(folder):
        frame_tree(("CamA", "CamB"), tag="ex")(folder)
        (Path(folder) / "Unmixed").mkdir()
        sentinel = np.arange(16, dtype=np.uint16).reshape(1, 4, 4)
        synth.write_tiff(Path(folder) / "Unmixed" / f"{FSN.format('CamA')}.tif", sentinel)

    return make


def zarr_case(cid, make, cams, factors, *args, error=False):
    return Case(
        cid,
        matlab="XR_unmix_channels_zarr",
        nargout=0,
        args=[files(cams, "zarr"), factors, "uuid", "t", *args, "parseCluster", False],
        python=py_zarr,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
    )


def _zarr_cases():
    ab = ("CamA", "CamB")
    t = frame_tree(ab, fmt="zarr", tag="zt")
    small = ["batchSize", [32, 24, 16], "blockSize", [16, 8, 8]]
    sig = ["unmixSigmas", np.array([[0, 0, 0], [1, 1, 1]], dtype=np.float64)]
    return [
        zarr_case("zarr.linear", t, ab, [1, -0.1], *small),
        zarr_case("zarr.linear_default_batches", t, ab, [1, -0.1]),
        zarr_case(
            "zarr.linear_one_batch_per_chunk",
            t,
            ab,
            [1, -1 / 3],
            "batchSize",
            [32, 32, 16],
            "blockSize",
            [16, 16, 16],
        ),
        zarr_case(
            "zarr.gaussian_border",
            t,
            ab,
            [1, -0.1],
            "mode",
            "gaussian",
            *sig,
            *small,
            "borderSize",
            [4, 4, 4],
        ),
        zarr_case("zarr.gaussian_no_border", t, ab, [1, -0.1], "mode", "gaussian", *sig, *small),
        zarr_case(
            "zarr.gaussian_aniso_border",
            t,
            ab,
            [0.8, -0.3],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0.5, 0.5, 1.0], [1.0, 1.5, 2.0]]),
            *small,
            "borderSize",
            [5, 3, 4],
        ),
        zarr_case(
            "zarr.gaussian_unfiltered_at_once",
            t,
            ab,
            [1, -0.1],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0, 0, 0], [0, 0, 0]], dtype=np.float64),
            *small,
            "borderSize",
            [2, 2, 2],
        ),
        zarr_case(
            "zarr.single_channels",
            frame_tree(ab, "single", "zarr", tag="zs"),
            ab,
            [1, -0.1],
            *small,
        ),
        zarr_case(
            "zarr.single_channels_gaussian",
            frame_tree(ab, "single", "zarr", tag="zs"),
            ab,
            [1, -0.1],
            "mode",
            "gaussian",
            *sig,
            *small,
            "borderSize",
            [4, 4, 4],
        ),
        zarr_case(
            "zarr.uint8_channels", frame_tree(ab, "uint8", "zarr", tag="z8"), ab, [1, -0.5], *small
        ),
        zarr_case(
            "zarr.three_channels_second",
            frame_tree(("CamA", "CamB", "CamC"), fmt="zarr", tag="z3"),
            ("CamA", "CamB", "CamC"),
            [-0.25, 1, 0],
            "channelInd",
            2,
            *small,
        ),
        zarr_case("zarr.save16bit_false_ignored", t, ab, [1, -0.1], "save16bit", False, *small),
    ]


def dataset_tree(cams=("CamA", "CamB"), fmt="tif", cls="uint16", tag="ds", extra=()):
    def make(folder):
        frame_tree(cams, cls, fmt, tag=tag)(folder)
        for t, name in enumerate(extra, 1):
            a = texture(SHAPE, (100, 1000), f"{tag}/extra{t}")
            synth.write_tiff(Path(folder) / f"{name}.tif", a)

    return make


def two_folders_tree(one_channel_each=False, incomplete=False):
    """two dataset folders d/a and d/b: both channels in each, or CamA in a and CamB in b; with
    incomplete, b has two more CamA time points"""

    def make(folder):
        a, b = Path(folder) / "a", Path(folder) / "b"
        a.mkdir()
        b.mkdir()
        if one_channel_each:
            frame_tree(("CamA",), tag="fa")(a)
            frame_tree(("CamB",), tag="fb")(b)
            return
        frame_tree(("CamA", "CamB"), tag="fa")(a)
        frame_tree(("CamA", "CamB"), tag="fb")(b)
        if incomplete:
            for t in (1, 2):
                name = FSN.format("CamA").replace("Iter_0000_", f"Iter_{t:04d}_")
                synth.write_tiff(b / f"{name}.tif", np.zeros(SHAPE, np.uint16))

    return make


def dataset_case(cid, make, *args, paths=None, error=False):
    return Case(
        cid,
        code="XR_unmix_channels_data_wrapper(in{:}); out = {};",
        nargout=0,
        args=[InPath("d") if paths is None else paths, *args, "uuid", "t", *LOCAL],
        python=py_dataset,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
    )


def _dataset_cases():
    pats = ["channelPatterns", ["CamA", "CamB"]]
    t = dataset_tree()
    return [
        # the demo's settings
        dataset_case(
            "ds.demo_linear",
            t,
            "unmixFactors",
            [1, -0.1],
            "mode",
            "linear",
            *pats,
            "channelInd",
            1,
            "constBackground",
            [100, 100],
            "constOffset",
            1,
            "save16bit",
            True,
        ),
        dataset_case(
            "ds.second_single",
            t,
            "unmixFactors",
            [-0.25, 1],
            "channelInd",
            2,
            *pats,
            "save16bit",
            False,
        ),
        dataset_case(
            "ds.ff",
            t,
            "unmixFactors",
            [1, -0.1],
            *pats,
            "FFCorrection",
            True,
            "lowerLimit",
            0.4,
            "FFImagePaths",
            FF[:2],
            "backgroundPaths",
            BK[:2],
        ),
        dataset_case(
            "ds.ff_lower_limit_text",
            t,
            "unmixFactors",
            [1, -0.1],
            *pats,
            "FFCorrection",
            True,
            "lowerLimit",
            1 / 3,
            "FFImagePaths",
            FF[:2],
            "backgroundPaths",
            BK[:2],
            "save16bit",
            False,
        ),
        dataset_case(
            "ds.zarr_in_out",
            dataset_tree(fmt="zarr", tag="dz"),
            "unmixFactors",
            [1, -0.1],
            *pats,
            "constBackground",
            100,
            "zarrFile",
            True,
            "saveZarr",
            True,
        ),
        dataset_case(
            "ds.large_file",
            dataset_tree(fmt="zarr", tag="dl"),
            "unmixFactors",
            [1, -0.1],
            *pats,
            "constBackground",
            100,
            "constOffset",
            5,
            "zarrFile",
            True,
            "saveZarr",
            True,
            "largeFile",
            True,
            "batchSize",
            [32, 24, 16],
            "blockSize",
            [16, 8, 8],
        ),
        dataset_case(
            "ds.large_file_gaussian",
            dataset_tree(fmt="zarr", tag="dl"),
            "unmixFactors",
            [1, -0.1],
            *pats,
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0, 0, 0], [1, 1, 1]], dtype=np.float64),
            "zarrFile",
            True,
            "saveZarr",
            True,
            "largeFile",
            True,
            "batchSize",
            [32, 24, 16],
            "blockSize",
            [16, 8, 8],
            "borderSize",
            [4, 4, 4],
        ),
        dataset_case(
            "ds.gaussian",
            t,
            "unmixFactors",
            [1, -0.1],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0, 0, 0], [1, 1, 1]], dtype=np.float64),
            *pats,
            "save16bit",
            False,
        ),
        # factors and sigmas the task strings round (10 decimals, 10 digits), the backgrounds
        # and the offset (15 digits)
        dataset_case(
            "ds.text_rounding",
            t,
            "unmixFactors",
            [1.000000000049, -1 / 3],
            "mode",
            "gaussian",
            "unmixSigmas",
            np.array([[0, 0, 0], [2 / 3, 1.23456789012345, 4 / 3]]),
            *pats,
            "constBackground",
            [100 / 3, 99.12345678901234],
            "constOffset",
            1 / 7,
            "save16bit",
            False,
        ),
        dataset_case(
            "ds.three_channels",
            dataset_tree(("CamA", "CamB", "CamC"), tag="d3"),
            "unmixFactors",
            [1, -0.1, -0.05],
            "channelPatterns",
            ["CamA", "CamB", "CamC"],
            "constBackground",
            [100, 100, 50],
        ),
        dataset_case(
            "ds.one_channel_per_folder",
            two_folders_tree(one_channel_each=True),
            "unmixFactors",
            [1, -0.1],
            *pats,
            paths=[Work("d/a"), Work("d/b")],
        ),
        dataset_case(
            "ds.two_folders_incomplete",
            two_folders_tree(incomplete=True),
            "unmixFactors",
            [1, -0.1],
            *pats,
            paths=[Work("d/a"), Work("d/b")],
        ),
        dataset_case("ds.wrong_factor_count", t, "unmixFactors", [1, -0.1, 0.5], *pats, error=True),
        dataset_case(
            "ds.wrong_background_count",
            t,
            "unmixFactors",
            [1, -0.1],
            *pats,
            "constBackground",
            [1, 2, 3],
            error=True,
        ),
        dataset_case(
            "ds.missing_channel",
            dataset_tree(("CamA",), tag="dm"),
            "unmixFactors",
            [1, -0.1],
            *pats,
            error=True,
        ),
    ]


def py_gauss(a, sigma):
    from petakit5d.core.axes import from_yxz
    from petakit5d.matlabcompat.filters import imgaussfilt3

    return imgaussfilt3(a, tuple(from_yxz(sigma)))


def _gauss_cases():
    """imgaussfilt3 at the sigmas of the unmixing cases (single data, MATLAB's [sy sx sz])"""
    vol = Value(lambda: to_matlab(synth.texture(SHAPE, 2.0, (100, 3000), CounterRNG(5, "b02/g"))))
    sigmas = {
        "text": [0.6666666667, 1.23456789, 1.333333333],
        "x": [1, 0.6666666667, 1],
        "z": [1, 1, 1.23456789],
    }
    return [
        Case(
            f"gauss.{k}",
            matlab="imgaussfilt3",
            args=[Vol(In("v")), v],
            python=py_gauss,
            out_roles=("vol",),
            inputs={"v": vol},
        )
        for k, v in sigmas.items()
    ]


CASES = (
    _gauss_cases()
    + _sum_cases()
    + _partial_cases()
    + _other_class_cases()
    + _frame_cases()
    + _zarr_cases()
    + _dataset_cases()
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# UnmixChannelsTest in intent, on (z, y, x) data and native zarr arrays
# ---------------------------------------------------------------------------------------------

LOCAL_KW = {"parse_cluster": False, "master_compute": True}


def kernel_cases():
    """mexMatchesTheMatlabCode's cases in Python's form"""
    out = []
    for i, msize in enumerate(SIZES):
        for k, classes in enumerate(CLASS_SETS):
            ims = [channel(py_shape(msize), c, f"k{i}/{k}/{j}") for j, c in enumerate(classes)]
            o = OFFSETS[(i + k + 2) % len(OFFSETS)]
            out.append((ims, FACTORS[k], BACKGROUNDS[k], o))
    return out


def ref(ims, factors, bgs, offset, dtype, im=None):
    nc = len(ims)
    f = np.asarray(factors, dtype=np.float64).reshape(-1)[:nc]
    b = np.zeros(nc) if bgs is None else np.asarray(bgs, dtype=np.float64).reshape(-1)[:nc]
    o = 0.0 if offset is None else float(offset)
    return unmix.unmix_channels_sum_ref(ims, f, b, o, dtype, im)


@pytest.mark.parametrize("out_class", ["uint8", "uint16", "single", "double"])
def test_mex_matches_the_matlab_code(out_class):
    for ims, f, b, o in kernel_cases():
        got = unmix.unmix_channels_sum(ims, f, b, o, out_class)
        same_bits(got, ref(ims, f, b, o, out_class))


def test_partial_sums_give_the_same_result():
    shape = (11, 23, 19)
    ims = [channel(shape, c, f"ps/{c}") for c in ("uint16", "single", "uint8", "double")]
    factors = [1, -0.1, 1 / 3, -2.345678901234]
    bgs = [100, 99.7, 0, 1 / 3]
    for out in ("uint16", "single", "uint8", "double"):
        expected = unmix.unmix_channels_sum(ims, factors, bgs, 1.7, out)
        im = None
        for c in range(4):
            im = unmix.unmix_channels_sum(
                ims[c : c + 1], factors[c : c + 1], bgs[c : c + 1], None, "", im
            )
        assert im.dtype == F32
        same_bits(unmix.unmix_channels_sum([], [], None, 1.7, out, im), expected)
        # in place: the partial sum is the output
        acc = None
        for c in range(4):
            acc = unmix.unmix_channels_sum(
                ims[c : c + 1], factors[c : c + 1], bgs[c : c + 1], None, "", acc, out=acc
            )
        same_bits(acc, im)


def test_sum_rounds_as_matlab():
    x = np.array([3, 0.49999997, 0.5, 1.5, 2.5, 65534.5, 65535.5, np.inf, NAN, -5, 1e-3], F32)
    same_bits(
        unmix.unmix_channels_sum([x], 1, None, None, "uint16"),
        np.array([3, 0, 1, 2, 3, 65535, 65535, 65535, 0, 0, 0], np.uint16),
    )
    same_bits(
        unmix.unmix_channels_sum([x], 1, None, None, "uint8"),
        np.array([3, 0, 1, 2, 3, 255, 255, 255, 0, 0, 0], np.uint8),
    )
    x = (CounterRNG(3, "b02/round").uniform(1000) * 4000).astype(F32)
    out = unmix.unmix_channels_sum([x], 0.1, None, None, "single")
    same_bits(out, (x.astype(np.float64) * 0.1).astype(F32))
    assert not np.array_equal(out, x * F32(0.1)), "The product was computed in single."
    out = unmix.unmix_channels_sum([x, x], [1, -0.1], None, None, "single")
    p = (x.astype(np.float64) * -0.1).astype(F32)
    same_bits(out, (x.astype(np.float64) + p.astype(np.float64)).astype(F32))
    x = (CounterRNG(3, "b02/round2").uniform(1000) * 1e-3).astype(F32)
    out = unmix.unmix_channels_sum([x], 1, None, 0.1, "single")
    same_bits(out, (x.astype(np.float64) + 0.1).astype(F32))
    assert not np.array_equal(out, x + F32(0.1)), "The offset was added in single."
    # max(0, -0) and max(0, NaN) are +0
    out = unmix.unmix_channels_sum([np.array([-0.0, NAN], F32)], 1, None, None, "single")
    assert out.view(np.uint32).tolist() == [0, 0]


def test_mex_rejects_invalid_input():
    ims = [np.zeros((4, 5), np.uint16), np.zeros((4, 5), np.uint16)]
    with pytest.raises(unmix.UnmixError) as e:
        unmix.unmix_channels_sum([], [], None, 0, "uint16")
    assert e.value.identifier == "unmix_channels:inputError"
    with pytest.raises(unmix.UnmixError) as e:
        unmix.unmix_channels_sum([ims[0], np.zeros((5, 4), np.uint16)], [1, 2], None, 0, "uint16")
    assert e.value.identifier == "unmix_channels:inputError"
    with pytest.raises(unmix.UnmixError) as e:
        unmix.unmix_channels_sum(
            [ims[0], np.zeros((2, 4, 5), np.uint16)], [1, 2], None, 0, "uint16"
        )
    assert e.value.identifier == "unmix_channels:inputError"
    with pytest.raises(unmix.UnmixError) as e:
        unmix.unmix_channels_sum(ims, [1], None, 0, "uint16")
    assert e.value.identifier == "MATLAB:badsubscript"
    with pytest.raises(unmix.UnmixError) as e:
        unmix.unmix_channels_sum(ims, [1, 2], [0], 0, "uint16")
    assert e.value.identifier == "MATLAB:badsubscript"
    with pytest.raises(unmix.UnmixError):
        unmix.unmix_channels_sum(ims, [1, 2], None, [0, 1], "uint16")
    with pytest.raises(unmix.UnmixError) as e:
        unmix.unmix_channels_sum([ims[0].astype(np.complex64)], [1], None, 0, "uint16")
    assert e.value.identifier == "unmix_channels:dataTypeError"
    # partial sums: single, of the channels' size
    with pytest.raises(unmix.UnmixError) as e:
        unmix.unmix_channels_sum(ims, [1, 2], None, 0, "", np.zeros((4, 5)))
    assert e.value.identifier == "unmix_channels:dataTypeError"
    with pytest.raises(unmix.UnmixError):
        unmix.unmix_channels_sum(ims, [1, 2], None, 0, "", np.zeros((5, 4), F32))
    # the kernel's own checks: no copies of other layouts or classes
    if unmix.extension_available():
        ext = unmix._ext()
        f = np.ones(1)
        out = np.empty(20, np.uint16)
        with pytest.raises(TypeError):
            ext.unmix_sum(
                [np.zeros((4, 5), np.uint16, order="F").T.reshape(-1)[::2]],
                f,
                f,
                0.0,
                True,
                None,
                out,
                1,
            )
        with pytest.raises(ValueError):
            ext.unmix_sum([np.zeros(20, np.int16)], f, f, 0.0, True, None, out, 1)
        with pytest.raises(ValueError):
            ext.unmix_sum([np.zeros(20, np.uint16)], f, f, 0.0, False, None, out, 1)
        with pytest.raises(TypeError):
            ext.unmix_sum([np.zeros(20, np.uint16)], np.ones(1, F32), f, 0.0, True, None, out, 1)


def test_other_classes_use_the_matlab_code():
    rng = CounterRNG(4, "b02/other")
    x = rng.integers(-100, 3000, (3, 9, 7), dtype=np.int16)
    y = rng.integers(1, 3000, (3, 9, 7), dtype=np.uint16)
    out1 = unmix.unmix_channels_sum([x, y], [1, -0.1], [10, 5], 2, "uint16")
    same_bits(out1, unmix.unmix_channels_sum([x.astype(F32), y], [1, -0.1], [10, 5], 2, "uint16"))
    out2 = unmix.unmix_channels_sum([y, y], [1, -0.1], None, 2, "int32")
    from petakit5d.core.dtypes import cast

    same_bits(out2, cast(unmix.unmix_channels_sum([y, y], [1, -0.1], None, 2, "single"), "int32"))


def two_channel_data(folder, zarr_file=False, tag=""):
    """twoChannelData: CamA and CamB of (24, 48, 40), CamB bleeding into CamA"""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    chans = channels_abc(SHAPE, tag)
    for cam in ("CamA", "CamB"):
        fn = folder / FSN.format(cam)
        if zarr_file:
            write_zarr(chans[cam], f"{fn}.zarr", block_size=(16, 16, 16))
        else:
            write_tiff(chans[cam], f"{fn}.tif")
    return str(folder) + "/", chans["CamA"], chans["CamB"]


def s(a):
    return np.asarray(a).astype(F32)


def mul(a, f):
    """single * double, rounded to single once"""
    return (s(a).astype(np.float64) * f).astype(F32)


def sub(a, b):
    return (s(a).astype(np.float64) - b).astype(F32)


def cast16(a):
    return unmix.unmix_channels_sum_ref([s(a)], [1], [0], 0, "uint16")


def test_linear_unmixing_as_in_demo(tmp_path):
    path, a, b = two_channel_data(tmp_path / "d")
    unmix.unmix_channels_dataset(
        path,
        unmix_factors=(1, -0.1),
        mode="linear",
        channel_patterns=("CamA", "CamB"),
        channel_ind=0,
        const_background=(100, 100),
        const_offset=1,
        save_16bit=True,
        **LOCAL_KW,
    )
    out = read_tiff(f"{path}Unmixed/{FSN.format('CamA')}.tif")
    v = sub(a, 100) + mul(sub(b, 100), -0.1)
    expected = cast16((np.maximum(v, 0).astype(np.float64) + 1).astype(F32))
    same_bits(out, expected)
    # only the unmixed channel is written
    assert len(list(Path(path, "Unmixed").glob("*.tif"))) == 1


def test_unmix_second_channel_as_single(tmp_path):
    path, a, b = two_channel_data(tmp_path / "d")
    unmix.unmix_channels_dataset(
        path,
        unmix_factors=(-0.25, 1),
        channel_ind=1,
        channel_patterns=("CamA", "CamB"),
        save_16bit=False,
        **LOCAL_KW,
    )
    out = read_tiff(f"{path}Unmixed/{FSN.format('CamB')}.tif")
    assert out.dtype == F32
    same_bits(out, np.maximum(mul(a, -0.25) + mul(b, 1), F32(0)))


def original_frame_unmix(
    raws,
    factors,
    ffs=(),
    bks=(),
    mode="linear",
    sigmas=None,
    ff=False,
    bg=None,
    offset=None,
    save16=True,
):
    """the per-channel code of XR_unmix_channels_frame before unmix_channels_sum"""
    from petakit5d.flatfield import flatfield_correct_volume
    from petakit5d.matlabcompat.filters import imgaussfilt3

    im = np.zeros(raws[0].shape, F32)
    for c, r in enumerate(raws):
        if factors[c] == 0:
            continue
        im_c = s(r)
        if ff:
            im_c = flatfield_correct_volume(im_c, ffs[c], bks[c], lower_limit=0.4, const_offset=0)
        if bg is not None:
            im_c = sub(im_c, bg[c])
        if mode == "gaussian" and not np.all(np.asarray(sigmas[c]) == 0):
            im_c = imgaussfilt3(im_c, sigmas[c])
        im = im + mul(im_c, factors[c])
    im = np.maximum(im, F32(0))
    im[im == 0] = 0  # (+0)
    if offset is not None and offset != 0:
        im = (im.astype(np.float64) + offset).astype(F32)
    return cast16(im) if save16 else im


def test_flat_field_correction_before_unmixing(tmp_path, recwarn):
    path, a, b = two_channel_data(tmp_path / "d")
    ff_path = Path(path) / "FF"
    ff_path.mkdir()
    x = np.arange(1, 41, dtype=np.float64)[np.newaxis].repeat(48, 0)
    bk = np.full((48, 40), 100, np.uint16)
    ffs = [np.floor(v + 0.5).astype(np.uint16) for v in (300 + 800 * x / 40, 1100 - 500 * x / 40)]
    fns = [str(ff_path / f"ff_Cam{c}.tif") for c in "AB"]
    bns = [str(ff_path / f"bk_Cam{c}.tif") for c in "AB"]
    for c in range(2):
        write_tiff(ffs[c], fns[c])
        write_tiff(bk, bns[c])
    unmix.unmix_channels_dataset(
        path,
        unmix_factors=(1, -0.1),
        channel_patterns=("CamA", "CamB"),
        ff_correction=True,
        lower_limit=0.4,
        ff_image_paths=fns,
        background_paths=bns,
        **LOCAL_KW,
    )
    out = read_tiff(f"{path}Unmixed/{FSN.format('CamA')}.tif")
    same_bits(out, original_frame_unmix([a, b], [1, -0.1], ffs, [bk, bk], ff=True))


def test_mem_allocate_is_in_gigabytes(tmp_path, framework_spy):
    path, _, _ = two_channel_data(tmp_path / "d")
    unmix.unmix_channels_dataset(
        path, unmix_factors=(1, -0.1), channel_patterns=("CamA", "CamB"), **LOCAL_KW
    )
    assert framework_spy.calls[0]["mem_allocate"] == pytest.approx(
        48 * 40 * 24 * 4 * 2.5 / 1024**3, rel=1e-12
    )


def test_images_without_all_channels_are_left_out(tmp_path, framework_spy):
    p1, _, _ = two_channel_data(tmp_path / "d1")
    p2, _, _ = two_channel_data(tmp_path / "d2", tag="2")
    for t in (1, 2):
        fn = FSN.format("CamA").replace("Iter_0000_", f"Iter_{t:04d}_")
        write_tiff(np.zeros(SHAPE, np.uint16), f"{p2}{fn}.tif")
    with pytest.warns(UserWarning, match="do not have all channels"):
        unmix.unmix_channels_dataset(
            [p1, p2], unmix_factors=(1, -0.1), channel_patterns=("CamA", "CamB"), **LOCAL_KW
        )
    assert len(framework_spy.calls) == 1
    call = framework_spy.calls[0]

    def unify(p):
        return os.path.normpath(p).replace("\\", "/")

    fn_a = FSN.format("CamA") + ".tif"
    assert [unify(p) for p in call["inputs"]] == [unify(p1 + fn_a), unify(p2 + fn_a)]
    assert len(call["tasks"]) == 2
    assert unify(call["tasks"][1].args[0][1]) == unify(f"{p2}{FSN.format('CamB')}.tif"), (
        "Not unmixed with CamB."
    )


def test_zarr_input_and_output(tmp_path):
    path, a, b = two_channel_data(tmp_path / "d", zarr_file=True)
    unmix.unmix_channels_dataset(
        path,
        unmix_factors=(1, -0.1),
        channel_patterns=("CamA", "CamB"),
        const_background=100,
        zarr_file=True,
        save_zarr=True,
        **LOCAL_KW,
    )
    out = read_zarr(f"{path}Unmixed/{FSN.format('CamA')}.zarr")
    same_bits(out, cast16(np.maximum(sub(a, 100) + mul(sub(b, 100), -0.1), F32(0))))


def zarr_fns(path):
    return [f"{path}{FSN.format(c)}.zarr" for c in ("CamA", "CamB")]


def test_zarr_function_matches_frame_function(tmp_path):
    path, a, b = two_channel_data(tmp_path / "d", zarr_file=True)
    unmix.unmix_channels_zarr(
        zarr_fns(path),
        (1, -0.1),
        result_dir_name="Unmixed",
        mode="linear",
        channel_ind=0,
        batch_size=(16, 32, 32),
        block_size=(16, 16, 16),
        **LOCAL_KW,
    )
    out = read_zarr(f"{path}Unmixed/{FSN.format('CamA')}.zarr")
    unmix.unmix_channels_file(
        zarr_fns(path), (1, -0.1), zarr_file=True, save_zarr=True, result_dir_name="Frame"
    )
    same_bits(out, read_zarr(f"{path}Frame/{FSN.format('CamA')}.zarr"))
    same_bits(out, cast16(np.maximum(mul(a, 1) + mul(b, -0.1), F32(0))))


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


def test_zarr_function_keeps_result_of_another_run(tmp_path, framework_spy):
    path, a, _ = two_channel_data(tmp_path / "d", zarr_file=True)
    out_fn = f"{path}Unmixed/{FSN.format('CamA')}.zarr"
    framework_spy.done = True
    framework_spy.action = lambda: write_zarr(a, out_fn)
    unmix.unmix_channels_zarr(
        zarr_fns(path),
        (1, -0.1),
        batch_size=(16, 32, 32),
        block_size=(16, 16, 16),
        parse_cluster=False,
        uuid="test",
    )
    same_bits(read_zarr(out_fn), a)
    assert not list(Path(path, "Unmixed").glob("*_test.zarr")), (
        "The temporary result was left behind."
    )
    assert os.path.normpath(framework_spy.calls[0]["final_out_fullpath"]) == os.path.normpath(
        out_fn
    )


def test_gaussian_unmixing_blurs_the_other_channel(tmp_path):
    from petakit5d.matlabcompat.filters import imgaussfilt3

    path, a, b = two_channel_data(tmp_path / "d", zarr_file=True)
    unmix.unmix_channels_zarr(
        zarr_fns(path), (1, -0.1), mode="gaussian", unmix_sigmas=(0, (1, 1, 1)), **LOCAL_KW
    )
    # the output keeps the input class (uint16)
    out = read_zarr(f"{path}Unmixed/{FSN.format('CamA')}.zarr")
    expected = cast16(np.maximum(mul(a, 1) + mul(imgaussfilt3(s(b), 1.0), -0.1), F32(0)))
    same_bits(out, expected)


def test_gaussian_unmixing_through_data_wrapper(tmp_path):
    path, a, b = two_channel_data(tmp_path / "d")
    unmix.unmix_channels_dataset(
        path,
        unmix_factors=(1, -0.1),
        mode="gaussian",
        unmix_sigmas=((0, 0, 0), (1, 1, 1)),
        channel_patterns=("CamA", "CamB"),
        save_16bit=False,
        **LOCAL_KW,
    )
    out = read_tiff(f"{path}Unmixed/{FSN.format('CamA')}.tif")
    expected = original_frame_unmix(
        [a, b], [1, -0.1], mode="gaussian", sigmas=[0, 1.0], save16=False
    )
    same_bits(out, expected)


FRAME_SETTINGS = [
    ([1, -0.1, 0], {"const_background": (100, 99.5, 98), "const_offset": 1, "save_16bit": True}),
    (
        [1, -0.1, 0.37],
        {
            "mode": "gaussian",
            "unmix_sigmas": (0, 1, (2, 0.7, 0.7)),
            "const_background": (100, 99.5, 98),
            "const_offset": 0.5,
            "save_16bit": False,
        },
    ),
    ([-0.25, 1, -1 / 3], {"ff_correction": True, "const_offset": -0.3, "save_16bit": True}),
    (
        [1, 0, -0.1234567],
        {
            "mode": "gaussian",
            "unmix_sigmas": (0, 1, 0),
            "ff_correction": True,
            "const_background": (3, 2, 1),
            "save_16bit": False,
        },
    ),
    (
        [1, -0.1, 0.37],
        {"const_background": (100, 99.5, 98), "const_offset": -0.3, "save_16bit": False},
    ),
    ([1, -0.1, 0.37], {"save_16bit": True}),
]


def _frame_folder(tmp_path):
    path = tmp_path / "d"
    path.mkdir()
    chans = channels_abc(SHAPE, "")
    raws = [chans[c] for c in ("CamA", "CamB", "CamC")]
    fns = []
    for c, cam in enumerate(("CamA", "CamB", "CamC")):
        fn = str(path / f"{FSN.format(cam)}.tif")
        write_tiff(raws[c], fn)
        fns.append(fn)
    ffs, bks = ff_images()
    ffn = [str(path / f"ff_{c}.tif") for c in (1, 2, 3)]
    bkn = [str(path / f"bk_{c}.tif") for c in (1, 2, 3)]
    for c in range(3):
        write_tiff(ffs[c], ffn[c])
        write_tiff(bks[c], bkn[c])
    return path, raws, fns, ffs, bks, ffn, bkn


def test_frame_matches_the_original_code(tmp_path):
    path, raws, fns, ffs, bks, ffn, bkn = _frame_folder(tmp_path)
    for i, (factors, kw) in enumerate(FRAME_SETTINGS):
        unmix.unmix_channels_file(
            fns,
            factors,
            result_dir_name=f"Unmixed{i}",
            ff_image_paths=ffn,
            background_paths=bkn,
            uuid=f"t{i}",
            **kw,
        )
        out = read_tiff(f"{path}/Unmixed{i}/{FSN.format('CamA')}.tif")
        sig = kw.get("unmix_sigmas")
        expected = original_frame_unmix(
            raws,
            factors,
            ffs,
            bks,
            kw.get("mode", "linear"),
            None
            if sig is None
            else [np.asarray(v, float) if np.ndim(v) else float(v) for v in sig],
            kw.get("ff_correction", False),
            kw.get("const_background"),
            kw.get("const_offset"),
            kw.get("save_16bit", True),
        )
        same_bits(out, expected)
    # single channels, their backgrounds subtracted in place
    spath = path / "single"
    spath.mkdir()
    sraws = [(s(r) + F32(0.37)).astype(F32) for r in raws]
    sfns = [str(spath / f"{FSN.format(c)}.tif") for c in ("CamA", "CamB", "CamC")]
    for c in range(3):
        write_tiff(sraws[c], sfns[c])
    for save16 in (True, False):
        unmix.unmix_channels_file(
            sfns,
            [1, -0.1, 0.37],
            result_dir_name=f"UnmixedSingle{int(save16)}",
            const_background=(100, 99.7, 98),
            const_offset=1,
            save_16bit=save16,
            uuid="s",
        )
        out = read_tiff(f"{spath}/UnmixedSingle{int(save16)}/{FSN.format('CamA')}.tif")
        same_bits(
            out,
            original_frame_unmix(
                sraws, [1, -0.1, 0.37], bg=(100, 99.7, 98), offset=1, save16=save16
            ),
        )


def test_zarr_blocks_match_the_original_code(tmp_path):
    from petakit5d.core.blocks import zarr_chunk_coordinates
    from petakit5d.matlabcompat.filters import imgaussfilt3

    path, a, b = two_channel_data(tmp_path / "d", zarr_file=True)
    factors = (1, -0.1)
    small = {"batch_size": (16, 32, 24), "block_size": (8, 16, 8)}
    unmix.unmix_channels_zarr(
        zarr_fns(path), factors, result_dir_name="Linear", **small, **LOCAL_KW
    )
    unmix.unmix_channels_zarr(
        zarr_fns(path),
        factors,
        result_dir_name="Gaussian",
        mode="gaussian",
        unmix_sigmas=(0, 1),
        border_size=(4, 4, 4),
        **small,
        **LOCAL_KW,
    )
    out = read_zarr(f"{path}Linear/{FSN.format('CamA')}.zarr")
    same_bits(out, cast16(np.maximum(mul(a, 1) + mul(b, -0.1), F32(0))))

    # the gaussian blocks filter each batch with its border
    out = read_zarr(f"{path}Gaussian/{FSN.format('CamA')}.zarr")
    bb, rb, lb = zarr_chunk_coordinates(
        SHAPE,
        batch_size=(16, 32, 24),
        block_size=(8, 16, 8),
        same_batch_size=False,
        border_size=(4, 4, 4),
    )
    expected = np.zeros(SHAPE, np.uint16)
    for i in range(bb.shape[0]):
        q, r, loc = bb[i], rb[i], lb[i]
        sl = (slice(q[0], q[3]), slice(q[1], q[4]), slice(q[2], q[5]))
        batch = mul(a[sl], 1) + mul(imgaussfilt3(s(b[sl]), 1.0), -0.1)
        batch = cast16(np.maximum(batch, F32(0)))
        expected[r[0] : r[3], r[1] : r[4], r[2] : r[5]] = batch[
            loc[0] : loc[3], loc[1] : loc[4], loc[2] : loc[5]
        ]
    same_bits(out, expected)

    # single channels and result
    sfns = [f"{path}{FSN.format(c)}_single.zarr" for c in ("CamA", "CamB")]
    sa, sb = (s(a) + F32(0.37)).astype(F32), (s(b) - F32(0.37)).astype(F32)
    write_zarr(sa, sfns[0], block_size=(16, 16, 16))
    write_zarr(sb, sfns[1], block_size=(16, 16, 16))
    unmix.unmix_channels_zarr(sfns, factors, result_dir_name="LinearSingle", **small, **LOCAL_KW)
    out = read_zarr(f"{path}LinearSingle/{FSN.format('CamA')}_single.zarr")
    same_bits(out, np.maximum(mul(sa, 1) + mul(sb, -0.1), F32(0)))


# ---------------------------------------------------------------------------------------------
# the port's own
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("threads", [1, 3, 7])
def test_results_do_not_depend_on_threads(threads):
    ims = [channel((3, 700, 600), c, f"th/{c}") for c in ("uint16", "single", "double")]
    a = unmix.unmix_channels_sum(
        ims, [1, -0.1, 0.3], [100, 99.7, 0], 1.7, "uint16", num_threads=threads
    )
    same_bits(a, ref(ims, [1, -0.1, 0.3], [100, 99.7, 0], 1.7, "uint16"))
    p = unmix.unmix_channels_sum(ims, [1, -0.1, 0.3], [100, 99.7, 0], None, "", num_threads=threads)
    same_bits(p, ref(ims, [1, -0.1, 0.3], [100, 99.7, 0], None, ""))


def test_subtract_is_single_minus_double():
    x = channel((5, 17, 13), "single", "subtract")
    y = x.copy()
    unmix._subtract(y, 99.7, 3)
    with np.errstate(all="ignore"):
        same_bits(y, (x.astype(np.float64) - 99.7).astype(F32))


@pytest.mark.parametrize("n_channels", [2, 3, 4])
def test_at_once_and_one_at_a_time_give_the_same_frames(tmp_path, n_channels, monkeypatch):
    path = tmp_path / "d"
    path.mkdir()
    cams = [f"Cam{c}" for c in "ABCD"[:n_channels]]
    fns = []
    for k, cam in enumerate(cams):
        fn = str(path / f"{FSN.format(cam)}.tif")
        write_tiff(texture(SHAPE, (100, 3000), f"oat{k}"), fn)
        fns.append(fn)
    factors = [1, -0.1, 0.2, -0.05][:n_channels]
    unmix.unmix_channels_file(
        fns, factors, const_background=[100] * n_channels, const_offset=1, result_dir_name="A"
    )
    # forced one at a time
    monkeypatch.setattr(unmix, "_bytes", lambda cls: 8)
    unmix.unmix_channels_file(
        fns, factors, const_background=[100] * n_channels, const_offset=1, result_dir_name="B"
    )
    same_bits(
        read_tiff(f"{path}/A/{FSN.format('CamA')}.tif"),
        read_tiff(f"{path}/B/{FSN.format('CamA')}.tif"),
    )


def test_frame_skips_channels_with_factor_zero(tmp_path, monkeypatch):
    path, raws, fns, *_ = _frame_folder(tmp_path)
    read = []
    orig = unmix._read
    monkeypatch.setattr(unmix, "_read", lambda p, z: read.append(p) or orig(p, z))
    unmix.unmix_channels_file(fns, [1, 0, -0.1])
    assert read == [fns[0], fns[2]]


def test_frame_skips_an_existing_result(tmp_path, capsys):
    path, raws, fns, *_ = _frame_folder(tmp_path)
    unmix.unmix_channels_file(fns, [1, -0.1, 0])
    capsys.readouterr()
    unmix.unmix_channels_file(fns, [1, -0.1, 0])
    assert "The output results exist, skip it!" in capsys.readouterr().out


def test_dataset_tasks_get_matlabs_text_values(tmp_path, framework_spy):
    path, _, _ = two_channel_data(tmp_path / "d")
    unmix.unmix_channels_dataset(
        path,
        unmix_factors=(1.000000000049, -1 / 3),
        mode="gaussian",
        unmix_sigmas=(0, (1.23456789012345, 1 / 3, 2 / 3)),
        channel_patterns=("CamA", "CamB"),
        const_background=(100 / 3, 99.12345678901234),
        const_offset=1 / 7,
        lower_limit=1 / 3,
        **LOCAL_KW,
    )
    t = framework_spy.calls[0]["tasks"][0]
    assert t.args[1] == [1.0, -0.3333333333]
    assert t.kwargs["unmix_sigmas"] == [[0.0], [1.234567890, 0.3333333333, 0.6666666667]]
    assert t.kwargs["const_background"] == [33.3333333333333, 99.1234567890123]
    assert t.kwargs["const_offset"] == 0.142857142857143
    assert t.kwargs["lower_limit"] == 1 / 3  # %.20d keeps a double
    assert t.kwargs["channel_ind"] == 0


def test_zarr_tasks_get_matlabs_text_values(tmp_path, framework_spy):
    path, _, _ = two_channel_data(tmp_path / "d", zarr_file=True)
    framework_spy.done = True
    unmix.unmix_channels_zarr(
        zarr_fns(path), (1, -1 / 3), mode="gaussian", unmix_sigmas=(0, 1 / 3), parse_cluster=False
    )
    t = framework_spy.calls[0]["tasks"][0]
    assert t.args[3] == [1.0, -0.333333333333333]
    assert t.args[4] == [[0.0], [0.333333333333333]]


def test_large_files_get_no_backgrounds_offsets_or_flat_fields(tmp_path, framework_spy):
    # MATLAB quirk: the zarr tasks of the dataset have no such parameters
    path, _, _ = two_channel_data(tmp_path / "d", zarr_file=True)
    unmix.unmix_channels_dataset(
        path,
        unmix_factors=(1, -0.1),
        channel_patterns=("CamA", "CamB"),
        const_background=100,
        const_offset=5,
        zarr_file=True,
        large_file=True,
        **LOCAL_KW,
    )
    call = framework_spy.calls[0]
    t = call["tasks"][0]
    assert t.func.endswith(":unmix_channels_zarr")
    assert not {"const_background", "const_offset", "ff_correction"} & set(t.kwargs)
    assert call["mem_allocate"] == pytest.approx(1024**3 / 1024**3 * 4 * 2.5)


def test_large_files_without_save_zarr_wait_for_a_tiff(tmp_path, framework_spy):
    # MATLAB quirk: the dataset's outputs follow save_zarr, but the zarr tasks write .zarr, so
    # without save_zarr the framework waits for a tiff that never comes (it runs each task
    # max_trial_num times, minutes in local mode)
    path, _, _ = two_channel_data(tmp_path / "d", zarr_file=True)
    for save_zarr in (False, True):
        unmix.unmix_channels_dataset(
            path,
            unmix_factors=(1, -0.1),
            channel_patterns=("CamA", "CamB"),
            zarr_file=True,
            large_file=True,
            save_zarr=save_zarr,
            **LOCAL_KW,
        )
    assert framework_spy.calls[0]["outputs"][0].endswith(f"{FSN.format('CamA')}.tif")
    assert framework_spy.calls[1]["outputs"][0].endswith(f"{FSN.format('CamA')}.zarr")


def test_extra_factors_and_backgrounds_are_not_used():
    # (MATLAB takes the first numel(ims) values)
    x = channel((3, 9, 7), "uint16", "extra")
    same_bits(
        unmix.unmix_channels_sum([x], [0.5, 9], [5, 7, 11], 1, "single"),
        unmix.unmix_channels_sum([x], [0.5], [5], 1, "single"),
    )


def test_a_zero_sigma_of_a_filtered_channel_is_matlabs_error(tmp_path):
    from petakit5d.matlabcompat.filters import MatlabError

    path, raws, fns, *_ = _frame_folder(tmp_path)
    with pytest.raises(MatlabError) as e:
        unmix.unmix_channels_file(
            fns, [1, -0.1, 0], mode="gaussian", unmix_sigmas=(0, (1, 0, 1), 0)
        )
    assert e.value.identifier == "MATLAB:imgaussfilt3:expectedPositive"


def test_zarr_result_has_the_class_of_the_unmixed_channel(tmp_path):
    path = tmp_path / "d"
    path.mkdir()
    chans = channels_abc(SHAPE, "cls")
    fns = [str(path / f"{FSN.format(c)}.zarr") for c in ("CamA", "CamB")]
    write_zarr(s(chans["CamA"]), fns[0], block_size=(16, 16, 16))
    write_zarr(chans["CamB"], fns[1], block_size=(16, 16, 16))
    unmix.unmix_channels_zarr(fns, (1, -0.1), **LOCAL_KW)
    assert read_zarr(f"{path}/Unmixed/{FSN.format('CamA')}.zarr").dtype == F32
    unmix.unmix_channels_zarr(fns, (-0.1, 1), channel_ind=1, **LOCAL_KW)
    assert read_zarr(f"{path}/Unmixed/{FSN.format('CamB')}.zarr").dtype == np.uint16


def test_zarr_result_keeps_the_channel_class(tmp_path):
    # MATLAB quirk: save_16bit is not used by the zarr function
    path, _, _ = two_channel_data(tmp_path / "d", zarr_file=True)
    unmix.unmix_channels_zarr(zarr_fns(path), (1, -0.1), save_16bit=False, **LOCAL_KW)
    assert read_zarr(f"{path}Unmixed/{FSN.format('CamA')}.zarr").dtype == np.uint16


def test_blocks_check_their_folders_and_skip_done_flags(tmp_path, capsys):
    path, _, _ = two_channel_data(tmp_path / "d", zarr_file=True)
    fns = zarr_fns(path)
    with pytest.raises(unmix.UnmixError, match="block directory"):
        unmix.unmix_channels_block(
            [0],
            fns,
            str(tmp_path / "o.zarr"),
            (1, -0.1),
            str(tmp_path / "no" / "f.mat"),
            [[0, 0, 0, 8, 8, 8]],
        )
    flag = tmp_path / "flags" / "blocks_1_1.mat"
    flag.parent.mkdir()
    flag.write_bytes(b"")
    assert (
        unmix.unmix_channels_block(
            [0], fns, str(tmp_path / "o.zarr"), (1, -0.1), str(flag), [[0, 0, 0, 8, 8, 8]]
        )
        is None
    )
    assert "already exist, skip them" in capsys.readouterr().out
    flag.unlink()
    with pytest.raises(unmix.UnmixError, match="output zarr file"):
        unmix.unmix_channels_block(
            [0], fns, str(tmp_path / "o.zarr"), (1, -0.1), str(flag), [[0, 0, 0, 8, 8, 8]]
        )
    with pytest.raises(unmix.UnmixError, match="number of unmix factors"):
        unmix.unmix_channels_block(
            [0], fns, str(tmp_path / "o.zarr"), (1,), str(flag), [[0, 0, 0, 8, 8, 8]]
        )


def test_linear_blocks_read_channels_with_factor_zero(tmp_path, monkeypatch):
    # MATLAB quirk: the blocks include channels with factor 0 (the frames skip them)
    from petakit5d.io import create_zarr

    path, a, _ = two_channel_data(tmp_path / "d", zarr_file=True)
    out = str(tmp_path / "o.zarr")
    create_zarr(out, data_size=SHAPE, block_size=(8, 16, 16), dtype="uint16")
    flag = tmp_path / "f" / "blocks_1_1.mat"
    flag.parent.mkdir()
    reads = []
    orig = unmix._Buffers.read
    monkeypatch.setattr(
        unmix._Buffers, "read", lambda self, c, p, box: reads.append(p) or orig(self, c, p, box)
    )
    unmix.unmix_channels_block([0], zarr_fns(path), out, (1, 0), str(flag), [[0, 0, 0, 24, 48, 40]])
    assert len(reads) == 2 and flag.exists()
    same_bits(read_zarr(out), a)


def test_matlab_layer_converts_sigmas_and_indices(tmp_path):
    from petakit5d import matlab
    from petakit5d.params import REGISTRY

    spec = REGISTRY.find("XR_unmix_channels_zarr")
    p = {q.matlab: q for q in spec.params}["unmixSigmas"].conv
    assert p.from_matlab(np.array([[0, 0, 0], [1, 2, 3.0]])) == ((0.0, 0.0, 0.0), (3.0, 1.0, 2.0))
    assert p.from_matlab((0.5, 1.0)) == (0.5, 1.0)
    assert p.to_matlab(((0.0, 0.0, 0.0), (3.0, 1.0, 2.0))) == ((0.0, 0.0, 0.0), (1.0, 2.0, 3.0))
    assert p.to_matlab((0.5, (3.0, 1.0, 2.0))) == ((0.5, 0.5, 0.5), (1.0, 2.0, 3.0))
    path, a, b = two_channel_data(tmp_path / "d", zarr_file=True)
    matlab.XR_unmix_channels_zarr(
        zarr_fns(path),
        [-0.25, 1],
        "channelInd",
        2,
        "batchSize",
        [32, 24, 16],
        "blockSize",
        [16, 8, 8],
        "parseCluster",
        False,
    )
    out = read_zarr(f"{path}Unmixed/{FSN.format('CamB')}.zarr")
    same_bits(out, cast16(np.maximum(mul(a, -0.25) + mul(b, 1), F32(0))))


def test_instruction_sets_give_the_same_bits():
    if not unmix.extension_available():
        pytest.skip("the extension is not built")
    ext = unmix._ext()
    if ext.best_isa() < 1:
        pytest.skip("no AVX2 on this CPU")
    ims = [channel((3, 70, 61), c, f"isa/{c}") for c in ("uint16", "single", "double", "uint8")]
    f = np.array([1, -0.1, 0.3, np.inf])
    b = np.array([100, 99.7, 0, -1 / 3])
    for out_cls, fin in (
        ("uint16", True),
        ("uint8", True),
        ("single", True),
        ("double", True),
        ("single", False),
    ):
        outs = []
        for isa in (0, 1):
            o = np.empty(ims[0].size, DT[out_cls])
            ext.unmix_sum([x.reshape(-1) for x in ims], f, b, -2.5, fin, None, o, 3, True, isa)
            outs.append(o)
        same_bits(outs[0], outs[1])
    x = ims[1].reshape(-1)
    r = [np.empty_like(x), np.empty_like(x)]
    for isa in (0, 1):
        ext.subtract(x, np.inf, r[isa], 2, True, isa)
    same_bits(r[0], r[1])


def test_the_tools_package_exports_the_unmixing():
    """petakit5d.tools exports the unmixing's functions and error lazily, and its table names
    each registered function it exports with its MATLAB name"""
    import re

    import petakit5d.tools as tools
    from petakit5d.params import REGISTRY

    for name in ("unmix_channels_sum", "unmix_channels_file", "unmix_channels_block",
                 "unmix_channels_gaussian_block", "unmix_channels_zarr", "unmix_channels_dataset",
                 "UnmixError"):  # fmt: skip
        assert name in tools.__all__ and name in dir(tools), name
        assert getattr(tools, name) is getattr(unmix, name)
    for mname, (where, kind) in REGISTRY.index().items():
        if kind != "register" or not where.startswith("petakit5d.tools."):
            continue
        pname = REGISTRY.find(mname).python_name.rsplit(".", 1)[-1]
        if pname in tools.__all__:
            assert re.search(rf"^``{pname}``\s+{mname}\b", tools.__doc__, re.M), (pname, mname)
