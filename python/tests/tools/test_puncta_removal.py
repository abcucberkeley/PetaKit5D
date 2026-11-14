"""Puncta removal for expansion microscopy (petakit5d.tools.puncta_removal, task B08).

Oracle cases against MATLAB: the helpers as MATLAB has them (thresholdOtsu, XR_locmax3d,
imregionalmax, fillmissing2 'linear' with MissingLocations, imresize3 'linear', regionprops3's
volume and principal axis filter), XR_ExM_PunctaRemovalPointDetection ('2d' and '3d', puncta
at the borders and corners, no puncta, several sigma thresholds) and XR_ExM_PunctaRemoval_zarr
on zarr inputs in MATLAB's layout (the updated method with default, scalar and base
thresholds, single and uint16 data, uint16 results, an offset, a box, batches with borders,
'3d' detection, several passes, all-zero data, the background estimation and the interpolated
background). Results are compared by voxels and metadata (zarr after the layout mapping),
parameters.mat and parameters.json by content and text, and the MIPs byte for byte.

The pixels the point detection fills (fillmissing2) may differ by one unit in the last place
of single (measured 26 of 10,647 filled pixels: MATLAB's barycentric weights are not exact, so
values at a rounding tie may round the other way); the cases that fill say so.

Then PunctaRemovalTest (tests/analysis) in intent on the port's (z, y, x) data and native zarr
arrays, and unit tests of the helpers, the quirks and the thread independence.
"""

from __future__ import annotations

import itertools
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.io import read_zarr, write_zarr
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, PerAxis, Tree, Value, Vol
from petakit5d.testing.compare import tol, tree
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import puncta_removal as pr

# ---------------------------------------------------------------------------------------------
# data: (z, y, x) volumes from bit-stable random numbers and Gaussian spots
# ---------------------------------------------------------------------------------------------

# MATLAB's barycentric weights are not exact (0.49999999999999994, -2^-55 where the exact
# weights are 1/2 and 0): its fills differ from the exact interpolation by about 1e-16 of the
# window's values, so a filled pixel may round to the next single, and an exact 0 between
# zeros may be +-3e-14
FILL_REASON = "puncta_removal: fillmissing2 fills: 1 ulp of single, or 1e-12 of the largest value"
FILL_TOL = tol(reason=FILL_REASON, ulp=1, rel_max=1e-12, max_frac=0.02)
FILL_TOL_U16 = tol(reason=FILL_REASON + " (uint16: 1)", abs=1, max_frac=0.01)


def expanded_volume(
    shape=(64, 64, 64),
    puncta=((19, 19, 19), (29, 39, 43), (44, 29, 49)),
    amplitude=20000.0,
    structure=((4, 35), (44, 60), (4, 30)),
    name="vol",
    cls="single",
    background=100.0,
):
    """PunctaRemovalTest's volume in (z, y, x): a dark background, a stained block with noise
    and bright puncta (Gaussian spots, sigma 1.5); 0-based positions and ranges"""
    vol = synth.beads(shape, np.asarray(puncta, dtype=np.float64), 1.5, amplitude, 0.0)
    v = vol.astype(np.float64) + background
    if structure is not None:
        sl = tuple(slice(a, b) for a, b in structure)
        n = CounterRNG(0, f"b08/{name}").normal(v[sl].shape)
        v[sl] = v[sl] + 4900 + 200 * n
    if cls == "uint16":
        return np.clip(np.floor(v + 0.5), 0, 65535).astype(np.uint16)
    return v.astype(np.float32)


def detection_volume(shape, name, puncta, amp=8000.0, noise=40.0, base=300.0, zero_frac=0.0):
    """a single volume for the point detection: a smooth background, noise and spots"""
    rng = CounterRNG(0, f"b08/pd/{name}")
    z = np.arange(shape[0], dtype=np.float64)[:, None, None]
    y = np.arange(shape[1], dtype=np.float64)[None, :, None]
    x = np.arange(shape[2], dtype=np.float64)[None, None, :]
    v = base + 2.0 * y + 1.5 * x + 3.0 * z + noise * rng.normal(shape)
    v = v + synth.beads(shape, np.asarray(puncta, dtype=np.float64), 1.5, amp, 0.0).astype(
        np.float64
    )
    if zero_frac:
        u = rng.spawn("zero").uniform(shape)
        v[u < zero_frac] = 0
    return np.maximum(v, 0).astype(np.float32)


def zarr_tree(vol, block=(32, 32, 32)):
    """a folder with expanded.zarr in MATLAB's layout (chunks (z, y, x))"""

    def make(folder):
        synth.write_zarr(Path(folder) / "expanded.zarr", vol, block_size=block, layout="matlab")

    return make


# ---------------------------------------------------------------------------------------------
# oracle cases: helpers
# ---------------------------------------------------------------------------------------------


def otsu_values():
    rng = CounterRNG(0, "b08/otsu")
    return {
        "single": (rng.uniform(4000) * 1000 + 0.5).astype(np.float32),
        "bimodal": np.concatenate(
            [rng.spawn("a").uniform(500) * 10 + 1, rng.spawn("b").uniform(700) * 10 + 60]
        ).astype(np.float32),
        "uint16": (rng.spawn("u").integers(1, 4000, 3000)).astype(np.uint16),
        "double": rng.spawn("d").uniform(2500) * 1e4 + 1e-3,
        "ties": np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype=np.float32),
        "three": np.array([1, 2, 3], dtype=np.float32),
        "const": np.array([5, 5, 5], dtype=np.float32),
        "one": np.array([7], dtype=np.float32),
        "grid": (1 + np.arange(256) / 255).astype(np.float32),
        "tiny": np.array([1e-6, 2e-6, 1], dtype=np.float32),
    }


def helper_volumes():
    rng = CounterRNG(0, "b08/helpers")
    out = {
        "rand": rng.uniform((11, 12, 13)).astype(np.float32),
        "ints": np.floor(rng.spawn("i").uniform((13, 15, 14)) * 4).astype(np.float32),
        "smooth": synth.gaussian_filter(rng.spawn("s").uniform((19, 20, 21)), 1.0).astype(
            np.float32
        ),
    }
    v = np.zeros((9, 10, 11), dtype=np.float32)
    v[4, 4, 4] = v[4, 4, 5] = 3
    v[1, 1, 1] = 1
    v[7, 8, 9] = 2
    v[0, 0, 0] = 5
    out["plateau"] = v
    return out


def regmax_images():
    rng = CounterRNG(0, "b08/regmax")
    m = np.zeros((10, 12), dtype=np.float32)
    m[2:4, 2:5] = 2
    m[3, 5] = 2
    m[7, 9] = 1
    m[0, :] = 0.5
    return {
        "rand": rng.uniform((20, 23)).astype(np.float32),
        "plateaus": m,
        "ints": np.floor(rng.spawn("i").uniform((31, 30)) * 4).astype(np.float32),
        "flat": np.ones((6, 7), dtype=np.float32),
    }


def fill_windows():
    """(image, missing) pairs: holes inside, at edges and corners, several holes, rings"""
    rng = CounterRNG(0, "b08/fill")
    out = {}
    shapes = [(8, 9), (12, 13), (9, 10), (17, 23), (25, 19), (30, 28)]
    for i, sh in enumerate(shapes):
        r = rng.spawn(i)
        a = (r.uniform(sh) * 1000 + 50).astype(np.float32)
        m = r.spawn("m").uniform(sh) > 0.8
        m = np.asarray(m)
        # grow the holes by one pixel along y and x
        g = m.copy()
        g[1:] |= m[:-1]
        g[:, 1:] |= m[:, :-1]
        if i % 2:
            g[: 1 + i % 3, :] |= r.spawn("e").uniform(sh[1]) > 0.5
        out[f"w{i}"] = (a, g)
    a = (rng.spawn("c").uniform((8, 8)) * 100).astype(np.float32)
    m = np.zeros((8, 8), dtype=bool)
    m[0:3, 0:3] = True
    out["corner"] = (a, m)
    m = np.zeros((10, 11), dtype=bool)
    m[2:8, 3:8] = True
    m[4:6, 4:7] = False
    out["ring"] = (a.repeat(2, 0)[:10].repeat(2, 1)[:, :11], m)
    a2 = out["w0"][0].copy()
    a2[3, 2] = np.float32(np.nan)
    out["nan_known"] = (a2, out["w0"][1])
    m = np.ones((5, 5), dtype=bool)
    m[2, :] = False
    out["collinear"] = (out["w0"][0][:5, :5].copy(), m)
    return out


def resize_inputs():
    rng = CounterRNG(0, "b08/resize")
    return {
        "a": ((rng.uniform((2, 2, 3)) * 100).astype(np.float32), (9, 17, 25)),
        "b": ((rng.spawn("b").uniform((2, 3, 1)) * 100).astype(np.float32), (11, 30, 4)),
        "c": ((rng.spawn("c").uniform((4, 4, 4)) * 100).astype(np.float32), (64, 64, 64)),
        "d": ((rng.spawn("d").uniform((2, 2, 2)) * 100).astype(np.float32), (3, 224, 225)),
        "same": ((rng.spawn("e").uniform((5, 3, 4)) * 100).astype(np.float32), (5, 3, 4)),
    }


def components_volume():
    bw = np.zeros((29, 30, 31), dtype=bool)
    bw[4:7, 2:8, 3:12] = True
    bw[2:25, 19:22, 19:21] = True
    z = np.arange(29)[:, None, None]
    y = np.arange(30)[None, :, None]
    x = np.arange(31)[None, None, :]
    bw |= ((x - 19) / 6) ** 2 + ((y - 9) / 3) ** 2 + ((z - 17) / 4) ** 2 <= 1
    bw[4, 14, 4] = True
    bw[9:12, 24:27, 2:5] = True
    bw[10, 25, 3] = False
    # a long thin line (principal axis over 35)
    for t in np.linspace(0, 1, 60):
        bw[int(2 + 25 * t), int(27 - 2 * t), int(1 + 29 * t)] = True
    return bw


def _ot(name):
    return Case(
        f"otsu.{name}",
        matlab="thresholdOtsu",
        args=[In("x")],
        python=pr.threshold_otsu,
        inputs={"x": Value(lambda: otsu_values()[name].reshape(-1, 1))},
    )


def _locmax(name):
    return Case(
        f"locmax3d.{name}",
        matlab="XR_locmax3d",
        args=[Vol(In("v")), [5, 5, 5]],
        python=pr.locmax3d,
        out_roles=("vol",),
        inputs={"v": Value(lambda: to_matlab(helper_volumes()[name]))},
    )


def _regmax(name):
    return Case(
        f"imregionalmax.{name}",
        matlab="imregionalmax",
        args=[Vol(In("m"))],
        python=pr.imregionalmax,
        out_roles=("vol",),
        inputs={"m": Value(lambda: regmax_images()[name])},
    )


def _fill(name):
    return Case(
        f"fillmissing2.{name}",
        code="out = {fillmissing2(in{1}, 'linear', MissingLocations=in{2})};",
        args=[Vol(In("a")), Vol(In("m"))],
        python=pr.fillmissing2_linear,
        out_roles=("vol",),
        compare=FILL_TOL,
        inputs={
            "a": Value(lambda: fill_windows()[name][0]),
            "m": Value(lambda: fill_windows()[name][1]),
        },
    )


def _resize(name):
    a, size = resize_inputs()[name]
    msize = [size[1], size[2], size[0]]
    return Case(
        f"imresize3.{name}",
        code="out = {imresize3(in{1}, in{2}, 'linear')};",
        args=[Vol(In("a")), PerAxis(msize, "yxz")],
        python=pr.imresize3_linear,
        out_roles=("vol",),
        inputs={"a": Value(lambda: to_matlab(resize_inputs()[name][0]))},
    )


def py_large_components(bw, det_vol_thrsh):
    return pr._large_components(np.asarray(bw, dtype=bool), float(det_vol_thrsh))


_LARGE = (
    "s = regionprops3(in{1}, 'voxelIdxList', 'volume', 'PrincipalAxisLength'); "
    "s = s(s.Volume > in{2} * 2 | s.PrincipalAxisLength(:, 1) > 35, :); "
    "B = false(size(in{1})); B(cat(1, cat(1, s.VoxelIdxList{:}))) = true; out = {B};"
)


def _helper_cases():
    cases = [_ot(n) for n in otsu_values()]
    cases += [_locmax(n) for n in helper_volumes()]
    cases += [_regmax(n) for n in regmax_images()]
    cases += [_fill(n) for n in fill_windows()]
    cases += [_resize(n) for n in resize_inputs()]
    for d in (5, 20, 200):
        cases.append(
            Case(
                f"regionprops3_filter.{d}",
                code=_LARGE,
                args=[Vol(In("bw")), d],
                python=py_large_components,
                out_roles=("vol",),
                inputs={"bw": Value(lambda: to_matlab(components_volume()))},
            )
        )
    return cases


# ---------------------------------------------------------------------------------------------
# oracle cases: point detection
# ---------------------------------------------------------------------------------------------

PD_PUNCTA = ((10, 12, 14), (20, 30, 9), (5, 40, 33), (25, 8, 40), (14, 22, 22))


def pd_inputs():
    return {
        "basic": detection_volume((31, 46, 44), "basic", PD_PUNCTA),
        "edges": detection_volume(
            (24, 37, 41), "edges", ((0, 0, 0), (12, 0, 20), (23, 36, 40), (11, 18, 0), (6, 36, 5))
        ),
        "sparse": detection_volume((26, 33, 35), "sparse", PD_PUNCTA[:3], zero_frac=0.6),
        "none": detection_volume((20, 21, 22), "none", (), noise=5.0),
        "dense": detection_volume(
            (28, 40, 38),
            "dense",
            [(z, y, x) for z in (6, 20) for y in (6, 18, 30) for x in (7, 19, 31)],
            amp=12000.0,
        ),
    }


def py_point_detection(v, *args):
    kw = dict(zip(args[::2], args[1::2], strict=True))
    names = {
        "Sigma": "sigma",
        "localWinSize": "local_win_size",
        "SigmaThrsh": "sigma_thrsh",
        "intThrsh": "int_thrsh",
        "initDetect": "init_detect",
        "detVolThrsh": "det_vol_thrsh",
    }
    return pr.puncta_removal_point_detection(v, **{names[k]: w for k, w in kw.items()})


def _pd(cid, name, *args, compare=FILL_TOL, error=False):
    return Case(
        f"point_detection.{cid}",
        matlab="XR_ExM_PunctaRemovalPointDetection",
        args=[Vol(In("v")), *args],
        python=py_point_detection,
        out_roles=("vol",),
        compare=compare,
        inputs={"v": Value(lambda: to_matlab(pd_inputs()[name]))},
        expect_error=error,
        py_error=(pr.PunctaRemovalError, ValueError, TypeError),
    )


WIN = ["localWinSize", PerAxis([15, 15, 23], "yxz")]
WIN25 = ["localWinSize", PerAxis([25, 25, 25], "yxz")]


def _pd_cases():
    return [
        _pd("basic.2d", "basic", *WIN25, "intThrsh", 3000),
        _pd("basic.defaults", "basic"),
        _pd("basic.3d", "basic", *WIN, "intThrsh", 3000, "initDetect", "3d"),
        _pd("basic.smallvol", "basic", *WIN25, "intThrsh", 3000, "detVolThrsh", 40),
        _pd("basic.sigma", "basic", *WIN, "intThrsh", 2000, "SigmaThrsh", 2, "Sigma", 1.5),
        _pd("edges.2d", "edges", *WIN, "intThrsh", 1000, "SigmaThrsh", 2),
        _pd("edges.3d", "edges", *WIN, "intThrsh", 3000, "initDetect", "3d"),
        _pd("sparse.2d", "sparse", *WIN, "intThrsh", 1000, "SigmaThrsh", 2),
        _pd("sparse.3d", "sparse", *WIN25, "intThrsh", 500, "SigmaThrsh", 2, "initDetect", "3d"),
        _pd("none.2d", "none", *WIN, "intThrsh", 3000, "SigmaThrsh", 2),
        _pd("dense.2d", "dense", "localWinSize", PerAxis([9, 9, 11], "yxz"), "SigmaThrsh", 2),
        _pd(
            "dense.oddwin",
            "dense",
            "localWinSize",
            PerAxis([8, 10, 12], "yxz"),
            "intThrsh",
            2000,
            "SigmaThrsh",
            2,
            "detVolThrsh",
            20000,
        ),
    ]


# ---------------------------------------------------------------------------------------------
# oracle cases: the zarr workflow
# ---------------------------------------------------------------------------------------------

LOCAL = ["parseCluster", False, "masterCompute", True]
SMALL = ["batchSize", [64, 64, 64], "blockSize", [32, 32, 32]]
ERRORS = (pr.PunctaRemovalError, ValueError, OSError, RuntimeError, IndexError)


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def py_zarr(d, *args):
    return call_matlab("XR_ExM_PunctaRemoval_zarr", d + "expanded.zarr", *args)


def py_bg_then_interp(d, est_args, interp_args):
    call_matlab(
        "XR_ExM_PunctaRemoval_zarr", d + "expanded.zarr", "mode", "bg_estimation", *est_args
    )
    bg = f"{d}Puncta_background_estimation/expanded_batch_32_32_32_background_values.zarr"
    call_matlab(
        "XR_ExM_PunctaRemoval_zarr",
        d + "expanded.zarr",
        "mode",
        "bg_interpolated",
        "bgFullpath",
        bg,
        *interp_args,
    )
    return []


def workflow_inputs():
    return {
        "test": expanded_volume(name="test"),
        "odd": expanded_volume(
            shape=(50, 70, 61),
            puncta=((10, 15, 12), (40, 60, 50), (25, 35, 30), (45, 5, 55)),
            structure=((5, 45), (20, 60), (10, 50)),
            name="odd",
        ),
        "u16": expanded_volume(
            shape=(48, 66, 57),
            puncta=((12, 50, 40), (30, 20, 10), (40, 40, 45)),
            structure=((4, 40), (10, 56), (20, 50)),
            name="u16",
            cls="uint16",
        ),
        "zeros": np.zeros((40, 36, 44), dtype=np.float32),
        # longer than a batch and its borders along y
        "tall": expanded_volume(
            shape=(36, 150, 40),
            puncta=((10, 20, 12), (25, 130, 30), (18, 75, 20)),
            structure=((4, 32), (30, 120), (5, 35)),
            name="tall",
        ),
        "holes": holes_volume(),
    }


def holes_volume():
    """an evenly stained volume with dark columns at 6% of (y, x) and two puncta: its z MIP
    lies above its blur nearly everywhere (the background mask loop of the base threshold)"""
    shape = (30, 64, 70)
    u = CounterRNG(0, "b08/holes").uniform(shape[1:])
    v = np.full(shape, 5000.0)
    v[:, u < 0.06] = 100.0
    v = v + synth.beads(shape, np.array([[10.0, 20, 30], [20, 40, 50]]), 1.5, 20000.0, 0.0)
    return v.astype(np.float32)


# the thresholds of batches whose y size is a multiple of 16: MATLAB's imgaussfilt3 of single
# data differs there from matlabcompat.filters in the last 16 rows by 1 ulp (reported to N04),
# which moves Otsu's threshold by about 1e-7 of its value
BG_TOL = tol(
    reason="puncta_removal: background thresholds of batches of 16k rows (imgaussfilt3, N04)",
    rel=1e-6,
)
BG_FILES = {"*blocks_*.mat": BG_TOL, "*_background_values.zarr": BG_TOL}


def _wf(cid, name, *args, values=None, block=(32, 32, 32), note="", per_file=None, folder="d"):
    spec = tree(values=values if values is not None else FILL_TOL, per_file=per_file)
    return Case(
        f"zarr.{cid}",
        code="XR_ExM_PunctaRemoval_zarr([in{1}, 'expanded.zarr'], in{2 : end}); out = {};",
        nargout=0,
        args=[InPath("d"), *args, *LOCAL],
        python=py_zarr,
        inputs={"d": Tree(zarr_tree(workflow_inputs()[name], block), folder)},
        tree=spec,
        py_error=ERRORS,
        note=note,
    )


def _wf_cases():
    cases = [
        _wf("default", "test", *SMALL),
        _wf(
            "scalar_thresholds",
            "test",
            *SMALL,
            "MinThreshold",
            0,
            "MaxThreshold",
            3000,
            "BaseThreshold",
            1000,
        ),
        _wf(
            "base_only",
            "odd",
            "batchSize",
            [32, 32, 32],
            "blockSize",
            [16, 16, 16],
            "BaseThreshold",
            800,
        ),
        _wf(
            "odd_batches",
            "odd",
            "batchSize",
            [32, 32, 32],
            "blockSize",
            [16, 16, 16],
            "intThrsh",
            5000,
            "offset",
            50,
        ),
        _wf(
            "uint16_save16",
            "u16",
            "batchSize",
            [40, 40, 40],
            "blockSize",
            [20, 20, 20],
            "Save16bit",
            True,
            "intThrsh",
            8000,
            values=FILL_TOL_U16,
        ),
        _wf("uint16_single", "u16", *SMALL, "MaxThreshold", 2500),
        _wf(
            "bbox",
            "odd",
            "batchSize",
            [32, 32, 32],
            "blockSize",
            [16, 16, 16],
            "bbox",
            [5, 3, 2, 60, 55, 44],
        ),
        _wf("detect3d", "test", *SMALL, "initDetect", "3d"),
        _wf(
            "passes",
            "test",
            *SMALL,
            "SigmaThrsh",
            [4, 3],
            "intThrsh",
            [10000, 4000],
            "detVolThrsh",
            [5000, 2000],
        ),
        _wf("zeros", "zeros", *SMALL),
        _wf(
            "bg_estimation",
            "odd",
            "mode",
            "bg_estimation",
            "batchSize",
            [32, 32, 32],
            "blockSize",
            [16, 16, 16],
            per_file=BG_FILES,
        ),
        _wf("tall", "tall", "batchSize", [32, 32, 32], "blockSize", [16, 16, 16]),
        _wf(
            "holes_base",
            "holes",
            "batchSize",
            [64, 64, 64],
            "blockSize",
            [32, 32, 32],
            "BaseThreshold",
            1000,
        ),
        # PunctaRemovalTest/parametersFileHoldsTheParameters: a % in the folder name
        _wf(
            "percent_folder",
            "test",
            *SMALL,
            folder="data_100%",
        ),
    ]
    est = ["batchSize", [32, 32, 32], "blockSize", [16, 16, 16], *LOCAL]
    interp = ["batchSize", [32, 32, 32], "blockSize", [16, 16, 16], "intThrsh", 5000, *LOCAL]
    for cid, name in (("zarr.bg_interpolated", "odd"), ("zarr.bg_interpolated_tall", "tall")):
        cases.append(_bg_interp_case(cid, name, est, interp))
    return cases


def _bg_interp_case(cid, name, est, interp):
    return Case(
        cid,
        code=(
            "d = in{1}; XR_ExM_PunctaRemoval_zarr([d, 'expanded.zarr'], 'mode', "
            "'bg_estimation', in{2}{:}); XR_ExM_PunctaRemoval_zarr([d, 'expanded.zarr'], "
            "'mode', 'bg_interpolated', 'bgFullpath', [d, "
            "'Puncta_background_estimation/expanded_batch_32_32_32_background_values.zarr'], "
            "in{3}{:}); out = {};"
        ),
        nargout=0,
        args=[InPath("d"), est, interp],
        python=py_bg_then_interp,
        inputs={"d": Tree(zarr_tree(workflow_inputs()[name]), "d")},
        tree=tree(values=FILL_TOL, per_file=BG_FILES),
        py_error=ERRORS,
    )


CASES = _helper_cases() + _pd_cases() + _wf_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# PunctaRemovalTest
# ---------------------------------------------------------------------------------------------

TEST_PUNCTA = ((19, 19, 19), (29, 39, 43), (44, 29, 49))  # MATLAB's [y x z] 1-based, as (z, y, x)
TEST_INNER = (slice(7, 32), slice(47, 57), slice(7, 27))  # the structure's interior


def _punctum_volume(tmp_path, folder=""):
    data = tmp_path / folder if folder else tmp_path
    data.mkdir(parents=True, exist_ok=True)
    vol = expanded_volume(name="test")
    fn = str(data / "expanded.zarr")
    write_zarr(vol, fn, block_size=(32, 32, 32))
    return fn, vol


def _verify_removed(fn, vol):
    out_fn = os.path.join(os.path.dirname(fn), "Puncta_Removed", "expanded.zarr")
    assert os.path.isdir(out_fn), "No puncta removal output."
    out = read_zarr(out_fn)
    assert out.shape == vol.shape
    for p in TEST_PUNCTA:
        assert out[p] < 0.1 * 20000, f"Punctum at {p} was not removed."
    a = float(np.mean(out[TEST_INNER], dtype=np.float64))
    b = float(np.mean(vol[TEST_INNER], dtype=np.float64))
    assert abs(a - b) <= 0.01 * abs(b), "The stained structure was not kept."


SMALL_PY = {"batch_size": (64, 64, 64), "block_size": (32, 32, 32), "parse_cluster": False}


def test_removes_puncta_with_default_thresholds(tmp_path):
    fn, vol = _punctum_volume(tmp_path)
    pr.puncta_removal_zarr(fn, **SMALL_PY)
    _verify_removed(fn, vol)


def test_removes_puncta_with_scalar_thresholds(tmp_path):
    fn, vol = _punctum_volume(tmp_path)
    pr.puncta_removal_zarr(fn, min_threshold=0, max_threshold=3000, base_threshold=1000, **SMALL_PY)
    _verify_removed(fn, vol)


def test_parameters_file_holds_the_parameters(tmp_path):
    fn, vol = _punctum_volume(tmp_path, "data_100%")
    pr.puncta_removal_zarr(fn, **SMALL_PY)
    _verify_removed(fn, vol)
    with open(os.path.join(os.path.dirname(fn), "Puncta_Removed", "parameters.json")) as f:
        p = json.load(f)
    assert p["zarrFullpath"] == fn
    # MATLAB's values: [y x z] sizes
    assert p["batchSize"] == [64, 64, 64]
    assert p["localWinSize"] == [15, 15, 23]


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _sign(v):
    return (v > 0) - (v < 0)


def _brute_triangles(pts):
    """the Delaunay triangles of points under the symbolic perturbation (lexicographic order:
    a cocircular point larger than the triangle's is outside; the largest of the triangle
    lifted), by testing every triple"""

    def orient(a, b, c):
        return _sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))

    def incircle(a, b, c, d):
        m = [[p[0] - d[0], p[1] - d[1], (p[0] - d[0]) ** 2 + (p[1] - d[1]) ** 2] for p in (a, b, c)]
        return _sign(
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        )

    def perturbed(a, b, c, d):
        s = incircle(a, b, c, d)
        if s:
            return s
        srt = sorted([a, b, c, d])
        for p in (srt[3], srt[2]):
            if p == d:
                return -1
            o = orient(a, b, d) if p == c else orient(a, d, c) if p == b else orient(d, b, c)
            if o:
                return o
        return -1

    tris = set()
    for a, b, c in itertools.combinations(pts, 3):
        o = orient(a, b, c)
        if o == 0:
            continue
        if o < 0:
            b, c = c, b
        if all(perturbed(a, b, c, d) < 0 for d in pts if d not in (a, b, c)):
            tris.add(frozenset((a, b, c)))
    return tris


@pytest.mark.parametrize("r2", [25, 50, 65, 85, 325])
def test_cocircular_cells_are_cut_as_the_perturbation_does(r2):
    pts = sorted((x, y) for x in range(-20, 21) for y in range(-20, 21) if x * x + y * y == r2)
    got = {frozenset(t) for t in pr._cell_triangles(pts)}
    assert got == _brute_triangles(pts)


def test_fill_matches_the_brute_force_triangulation():
    rng = CounterRNG(0, "b08/brute")
    a = rng.uniform((9, 10)) * 100
    m = np.zeros((9, 10), dtype=bool)
    m[3:6, 3:7] = True
    m[6, 5] = True
    f = pr.fillmissing2_linear(a, m)
    known = [(int(r), int(c)) for r, c in zip(*np.nonzero(~m), strict=True)]
    tris = _brute_triangles(known)
    for q in zip(*np.nonzero(m), strict=True):
        q = (int(q[0]), int(q[1]))
        vals = [pr._interp([tuple(t)], q, a) for t in tris]
        vals = [v for v in vals if v is not None]
        assert vals, q
        assert f[q] == vals[0]


def test_fill_is_nan_outside_the_known_pixels_hull():
    a = np.arange(64, dtype=np.float32).reshape(8, 8)
    m = np.zeros((8, 8), dtype=bool)
    m[0:3, 0:3] = True
    f = pr.fillmissing2_linear(a, m)
    assert np.isnan(f[0, 0]) and np.isnan(f[0, 1]) and np.isnan(f[1, 0])
    # MATLAB's NaN
    assert f.view(np.uint32)[0, 0] == 0xFFC00000
    assert not np.isnan(f[2, 2])


def test_fill_leaves_nan_known_pixels_out_and_keeps_them():
    a = np.arange(1, 50, dtype=np.float32).reshape(7, 7)
    a[2, 3] = np.nan
    m = np.zeros((7, 7), dtype=bool)
    m[3, 3] = True
    f = pr.fillmissing2_linear(a, m)
    assert np.isnan(f[2, 3])
    assert np.isfinite(f[3, 3])


def test_fill_without_a_triangle_leaves_the_image():
    a = np.arange(25, dtype=np.float32).reshape(5, 5)
    m = np.ones((5, 5), dtype=bool)
    m[2, :] = False
    np.testing.assert_array_equal(pr.fillmissing2_linear(a, m), a)
    m[:] = True
    np.testing.assert_array_equal(pr.fillmissing2_linear(a, m), a)


def test_fill_of_a_plane_is_exact():
    y, x = np.mgrid[0:12, 0:13]
    a = (3.0 * y + 2.0 * x + 1.0).astype(np.float32)
    m = np.zeros(a.shape, dtype=bool)
    m[3:9, 4:10] = True
    np.testing.assert_array_equal(pr.fillmissing2_linear(a, m), a)


def test_otsu_scales_from_the_smallest_positive_value():
    # (MATLAB: the threshold of [0; 1] is 1, of [1; 2; 3] 1 + 2 * 63.5/255)
    assert pr.threshold_otsu(np.array([0, 1], np.float32)) == 1.0
    assert pr.threshold_otsu(np.array([1, 2, 3], np.float32)) == 63.5 / 255 * 2 + 1
    assert pr.threshold_otsu(np.array([5, 5], np.float32)) == 5.0
    with pytest.raises(pr.PunctaRemovalError) as e:
        pr.threshold_otsu(np.zeros(0, np.float32))
    assert e.value.identifier == "MATLAB:innerdim"


def test_otsu_doesnt_depend_on_the_order():
    v = otsu_values()["bimodal"]
    assert pr.threshold_otsu(v) == pr.threshold_otsu(v[::-1]) == pr.threshold_otsu(np.sort(v))


def test_locmax3d_needs_a_unique_maximum():
    v = np.zeros((7, 8, 9), dtype=np.float32)
    v[3, 3, 3] = v[3, 3, 4] = 2
    v[5, 6, 7] = 1
    out = pr.locmax3d(v, 5)
    assert out[3, 3, 3] == 0 and out[3, 3, 4] == 0
    assert out[5, 6, 7] == 1
    v[5, 6, 3] = 1  # (another 1 two voxels away in x: still the maximum? no, a tie)
    assert pr.locmax3d(v, 5)[5, 6, 7] == 1  # (5,6,3) is 4 voxels away: outside the window
    v[5, 6, 5] = 1
    assert pr.locmax3d(v, 5)[5, 6, 7] == 0


def test_regionalmax_plateau_next_to_a_higher_value():
    m = np.zeros((5, 6), dtype=np.float32)
    m[1, 1:4] = 2
    m[1, 4] = 3
    bw = pr.imregionalmax(m)
    assert not bw[1, 1:4].any() and bw[1, 4]
    with pytest.raises(pr.PunctaRemovalError) as e:
        pr.imregionalmax(np.array([[1.0, np.nan]], np.float32))
    assert e.value.identifier == "MATLAB:imregionalmax:expectedNonNaN"


def test_resize_keeps_a_constant_and_enlarges_in_order_of_scale():
    a = np.full((2, 3, 4), 7.25, dtype=np.float32)
    np.testing.assert_array_equal(
        pr.imresize3_linear(a, (5, 9, 4)), np.full((5, 9, 4), 7.25, np.float32)
    )
    with pytest.raises(ValueError):
        pr.imresize3_linear(a, (1, 9, 9))


def test_principal_axis_threshold_is_exact():
    # a straight line of n voxels along x: variance (n^2 - 1)/12, PAL = 2 sqrt(5 (var + 1/12))
    # = sqrt(5/3) n: over 35 from n = 28 (27.11 for 35 exactly)
    for n, removed in ((27, False), (28, True)):
        bw = np.zeros((3, 3, 40), dtype=bool)
        bw[1, 1, :n] = True
        assert pr._large_components(bw, 1e9)[1, 1, 0] == removed
    bw = np.zeros((3, 3, 40), dtype=bool)
    bw[1, 1, :10] = True
    assert pr._large_components(bw, 4)[1, 1, 0]  # 10 > 2 * 4 voxels
    assert not pr._large_components(bw, 5)[1, 1, 0]


# ---------------------------------------------------------------------------------------------
# quirks, the first method, threads
# ---------------------------------------------------------------------------------------------


def test_removed_regions_are_the_thresholded_components():
    # MATLAB's watershed mask comes after the voxels to remove are taken: the filled pixels are
    # those of the thresholded component around each peak, not of its watershed region
    v = detection_volume((20, 30, 31), "one", ((10, 15, 15),), amp=12000.0)
    out = pr.puncta_removal_point_detection(v, local_win_size=(25, 25, 25), int_thrsh=3000)
    from petakit5d.matlabcompat import stats
    from petakit5d.matlabcompat.filters import imgaussfilt3

    sm = imgaussfilt3(v, 2.5)
    p = np.unravel_index(int(np.argmax(sm)), sm.shape)
    sl = pr._bbox_window(p, (12, 12, 12), sm.shape)
    w = sm[sl]
    t = pr._plus(
        min(pr.threshold_otsu(w[w > 0]), 3000.0) * 0.2, pr._times(stats.median(w)[()], 0.8)
    )
    comp = np.zeros(v.shape, dtype=bool)
    comp[sl] = pr._component_at(w > t, tuple(p[d] - sl[d].start for d in range(3)))
    changed = out != v
    assert changed.any() and not (changed & ~comp).any()


def test_corner_fill_gives_nan_then_regionalmax_raises():
    # a punctum at a corner: its fill window has no known pixels around the corner
    v = detection_volume((16, 20, 21), "corner", ((8, 0, 0),), amp=20000.0, noise=0.0)
    out = pr.puncta_removal_point_detection(v, local_win_size=(23, 15, 15), int_thrsh=3000)
    assert np.isnan(out).any()
    with pytest.raises(pr.PunctaRemovalError) as e:
        pr.puncta_removal_point_detection(out, local_win_size=(23, 15, 15), int_thrsh=3000)
    assert e.value.identifier == "MATLAB:imregionalmax:expectedNonNaN"


def test_interpolated_background_one_value_deep_in_z_raises(tmp_path):
    vol = expanded_volume(
        shape=(20, 40, 44), puncta=(), structure=((2, 18), (5, 30), (5, 40)), name="flat"
    )
    fn = str(tmp_path / "expanded.zarr")
    write_zarr(vol, fn, block_size=(16, 16, 16))
    kw = {"batch_size": (32, 32, 32), "block_size": (16, 16, 16), "parse_cluster": False}
    pr.puncta_removal_zarr(fn, mode="bg_estimation", **kw)
    bg = str(
        tmp_path / "Puncta_background_estimation" / "expanded_batch_32_32_32_background_values.zarr"
    )
    assert read_zarr(bg).ndim == 2  # (one z value)
    with pytest.raises(pr.PunctaRemovalError) as e:
        pr.puncta_removal_zarr(fn, mode="bg_interpolated", bg_fullpath=bg, **kw)
    assert e.value.identifier == "MATLAB:badsubscript"


def test_mask_threshold_only_with_a_base_threshold(monkeypatch):
    calls = []
    from petakit5d.matlabcompat import filters

    orig = filters.imgaussfilt

    def spy(*a, **k):
        calls.append(1)
        return orig(*a, **k)

    monkeypatch.setattr(filters, "imgaussfilt", spy)
    im = filters.imgaussfilt3(expanded_volume(name="test"), 2.5)
    t0 = pr._batch_threshold(im, 2.5, 99.9, 0, None, None, None)
    assert not calls
    t1 = pr._batch_threshold(im, 2.5, 99.9, 0, None, 1000, None)
    assert calls and t1 != t0
    assert isinstance(t0, float)


def test_first_method_block(tmp_path):
    from petakit5d.core.blocks import zarr_chunk_coordinates
    from petakit5d.io import create_zarr
    from petakit5d.matlabcompat import label, stats
    from petakit5d.matlabcompat.filters import imgaussfilt3

    vol = expanded_volume(
        shape=(40, 48, 44), puncta=(), structure=((5, 30), (10, 40), (5, 30)), name="first"
    )
    fn = str(tmp_path / "in.zarr")
    write_zarr(vol, fn, block_size=(16, 16, 16))
    out = str(tmp_path / "out.zarr")
    create_zarr(out, data_size=vol.shape, block_size=(16, 16, 16), dtype="single")
    (tmp_path / "flags").mkdir()
    bb, rb, lb = zarr_chunk_coordinates(
        vol.shape,
        batch_size=(64, 64, 64),
        block_size=(16, 16, 16),
        same_batch_size=False,
        border_size=0,
    )
    flag = str(tmp_path / "flags" / "b.mat")
    pr.puncta_removal_block(
        list(range(len(bb))), fn, out, flag, bb.tolist(), rb.tolist(), lb.tolist()
    )
    assert os.path.isfile(flag)
    p = stats.prctile(vol, 99.9)[0]
    t = pr.threshold_otsu(vol[(vol > 0) & (vol < p)])
    keep = label.bwareaopen(imgaussfilt3(vol, 2.5) > np.float64(t), 1000, 26)
    np.testing.assert_array_equal(read_zarr(out), vol * keep)
    # an existing flag: skipped
    assert (
        pr.puncta_removal_block([0], fn, out, flag, bb.tolist(), rb.tolist(), lb.tolist()) is None
    )


def test_point_detection_doesnt_change_its_input():
    v = pd_inputs()["basic"]
    c = v.copy()
    pr.puncta_removal_point_detection(v, local_win_size=(23, 15, 15), int_thrsh=3000)
    np.testing.assert_array_equal(v, c)


@pytest.mark.parametrize("init", ["2d", "3d"])
def test_point_detection_doesnt_depend_on_threads(init):
    v = pd_inputs()["basic"]
    kw = {"local_win_size": (23, 15, 15), "int_thrsh": 3000, "init_detect": init}
    a = pr.puncta_removal_point_detection(v, num_threads=1, **kw)
    b = pr.puncta_removal_point_detection(v, num_threads=4, **kw)
    assert a.tobytes() == b.tobytes()


def test_workflow_doesnt_depend_on_the_batches(tmp_path):
    # batches cover the image with borders of 100 (larger than these volumes): one batch or
    # many give the same result
    vol = workflow_inputs()["odd"]
    res = []
    for k, b in enumerate(((64, 96, 96), (32, 32, 32))):
        d = tmp_path / f"r{k}"
        d.mkdir()
        fn = str(d / "expanded.zarr")
        write_zarr(vol, fn, block_size=(16, 16, 16))
        pr.puncta_removal_zarr(
            fn, batch_size=b, block_size=(16, 16, 16), parse_cluster=False, int_thrsh=5000
        )
        res.append(read_zarr(str(d / "Puncta_Removed" / "expanded.zarr")))
    np.testing.assert_array_equal(res[0], res[1])


def test_workflow_skips_an_existing_result(tmp_path, capsys):
    fn, _ = _punctum_volume(tmp_path)
    os.makedirs(tmp_path / "Puncta_Removed" / "expanded.zarr")
    assert pr.puncta_removal_zarr(fn, **SMALL_PY) is None
    assert "exists, skip" in capsys.readouterr().out


def test_workflow_makes_mips_before_publishing(tmp_path, monkeypatch):
    from petakit5d.tools import mip

    seen = []
    orig = mip.mip_zarr

    def spy(path, **kw):
        seen.append((path, os.path.isdir(str(tmp_path / "Puncta_Removed" / "expanded.zarr"))))
        return orig(path, **kw)

    monkeypatch.setattr(mip, "mip_zarr", spy)
    fn, _ = _punctum_volume(tmp_path)
    pr.puncta_removal_zarr(fn, **SMALL_PY)
    assert seen and not seen[0][1] and "_" in os.path.basename(seen[0][0])
    assert os.path.isfile(tmp_path / "Puncta_Removed" / "MIPs" / "expanded_MIP_z.tif")
    # the flags are removed after publishing
    assert os.listdir(tmp_path / "Puncta_Removed" / "zarr_flag") == []


def test_offset_goes_through_matlab_text():
    from petakit5d.core.text import through_text

    # %d of a fraction is %e with 6 digits
    assert through_text(0.123456789, "%d") == 0.1234568
    assert math.isclose(through_text(50, "%d"), 50)


@pytest.mark.parametrize("dtype", [np.uint16, np.uint8, np.int16])
def test_otsu_of_small_integers_counts_as_doubles(dtype):
    # (the counts of each integer value, then their bins: the bins of the values as doubles)
    rng = CounterRNG(0, f"b08/otsu_int/{np.dtype(dtype).name}")
    hi = 200 if dtype == np.uint8 else 3000
    v = rng.integers(0, hi, 5000).astype(dtype)
    v[:50] = 0
    if dtype == np.int16:
        v[50:80] = -7
    assert pr.threshold_otsu(v) == pr.threshold_otsu(v.astype(np.float64))
    w = np.array([0, 0, 10, 20, 30], dtype)
    assert pr.threshold_otsu(w) == pr.threshold_otsu(w.astype(np.float32))


@pytest.mark.parametrize("threads", [1, 4])
def test_area_open_equals_bwareaopen(threads):
    from petakit5d.matlabcompat import label

    u = CounterRNG(0, "b08/areaopen").uniform((70, 130, 125))
    bw = u > 0.7
    for p in (0, 1, 5, 40):
        np.testing.assert_array_equal(
            pr._area_open(bw, p, threads), label.bwareaopen(bw, p, 26, threads=threads)
        )


def test_integer_threshold_comparison_is_exact():
    a = np.array([0, 1, 2, 3, 65535], dtype=np.uint16)
    for t in (-1.5, 0.0, 0.5, 1.0, 2.999999999, 1e9, 65534.5):
        np.testing.assert_array_equal(pr._gt_scalar(a, t), a.astype(np.float64) > t)
    np.testing.assert_array_equal(pr._gt_scalar(a, np.inf), np.zeros(5, bool))
