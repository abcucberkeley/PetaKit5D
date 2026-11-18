"""Fourier shell correlation (petakit5d.tools.fsc, task B05).

MATLAB's FSCAnalysisTest (5 tests) is ported below: noisierIterationHasWorseResolution,
existingResultsAreNotRecomputed and resultIsMovedIntoPlace as tests of the dataset function
(and oracle cases of the wrapper), oneImageMatchesResizingTheWholeImage (its 8 image cases) as a
test against resizing the whole image and as oracle cases, shellCacheFollowsParameters as a test
of the kept shells and as oracle cases.

Oracle cases against MATLAB:

- hamming(n) (the Signal Processing Toolbox's window, MATLAB's bits);
- XR_FSC of double volumes: odd and even sizes, sizes other than halfSize (padded and cropped),
  dr, dtheta, skipConeRegion, the three threshold methods, and its errors (an odd number of
  angles, all-zero volumes, a halfSize that is not cubic);
- XR_FSC_resolution of FSC structs MATLAB made;
- XR_one_image_FSC_analysis on FSCAnalysisTest's eight images (resized up in z and in x and y,
  padded, isotropic, an even half size, bbox with 'yz', 'xy' with clipPer) and the one-bit and
  half-bit methods;
- XR_one_image_FSC_analysis_frame with one box and several, and with an existing result;
- XR_FSC_analysis_wrapper on decon-like iterations (the MATLAB test's), a rerun, the multi-region
  modes, and XR_two_image_FSC_analysis's failure.

Each FFT case runs on both FFT backends: ``.exact`` (the FFTW add-on; skipped without it or
where MATLAB's FFT bits don't hold) gives MATLAB's FSC curves bit for bit, ``.permissive``
within decision record 0006's double tolerance. The resolutions come from the published
smoothing spline with csaps's documented default, not MathWorks's code (the module): every
radius measured is MATLAB's, and they are compared with the tolerance proposed in
docs/parity/tools.md (one grid step).
"""

from __future__ import annotations

import math
import os

import numpy as np
import pytest

from petakit5d import fft as F
from petakit5d.core import dtypes
from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, File, In, InPath, Tree, Value, Work, mcell
from petakit5d.testing.compare import (
    Custom,
    MArray,
    MCell,
    MStruct,
    ignore,
    to_mvalue,
    tol,
    tree,
)
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import figures
from petakit5d.tools import fsc as FS

HAS_MPL = figures.available()
BACKENDS = ("exact", "permissive")
PI = math.pi

# ---------------------------------------------------------------------------------------------
# comparisons
# ---------------------------------------------------------------------------------------------

PERMISSIVE = tol(
    reason="0006: FFT-based outputs with the permissive backend, double precision (max <= 1e-13 "
    "of the largest value, relative L2 <= 1e-14)",
    rel_max=1e-13,
    rel_l2=1e-14,
)

# the resolution (docs/parity/tools.md, proposed): the port's smoothing spline is the published
# method with csaps's documented default, not MathWorks's code, so its last bits differ; every
# radius found on MATLAB's grid of 1e-4 is MATLAB's (all measured are), allowing RES_STEPS step
# where the smoothed curve meets the threshold within those last bits; the resolutions within
# RES_REL (what one step at the smallest radii gives)
STEP = 1e-4
RES_STEPS = 1
RES_REL = 0.01


def _data(v) -> np.ndarray:
    if not isinstance(v, MArray) or v.data is None:
        raise TypeError("not an array with data")
    return np.asarray(v.data, dtype=np.float64).ravel(order="F")


def _structs(v, path: str) -> list[tuple[str, dict]]:
    if isinstance(v, MCell):
        out = []
        for k, item in enumerate(v.items or []):
            out += _structs(item, f"{path}{{{k + 1}}}")
        return out
    if isinstance(v, MStruct):
        return [(path, v.items[0])]
    raise TypeError(f"{path}: not a struct or cell")


def _bits_equal(a: np.ndarray, b: np.ndarray) -> bool:
    return a.shape == b.shape and np.array_equal(a.view(np.uint64), b.view(np.uint64))


def _res_problems(ref, act, steps: float, rel: float) -> list[str]:
    probs = []
    try:
        rs, acts = _structs(ref, ""), _structs(act, "")
    except TypeError as e:
        return [str(e)]
    if len(rs) != len(acts):
        return [f"{len(rs)} structs in MATLAB's, {len(acts)} in Python's"]
    for (path, r), (_, a) in zip(rs, acts, strict=True):
        if sorted(r) != sorted(a):
            probs.append(f"{path}: fields {sorted(r)} against {sorted(a)}")
            continue
        for f in ("comp_angles", "thetas"):
            if not _bits_equal(_data(r[f]), _data(a[f])):
                probs.append(f"{path}.{f} differs")
        pr, pa = _data(r["points"]), _data(a["points"])
        if pr.shape != pa.shape:
            probs.append(f"{path}.points: {pr.size} values against {pa.size}")
            continue
        d = np.round(np.abs(pr - pa) / STEP)
        if np.nanmax(d, initial=0) > steps:
            probs.append(f"{path}.points: up to {np.nanmax(d):.0f} steps of 1e-4 (limit {steps:g})")
        probs += _rel_problems(_data(r["res"]), _data(a["res"]), rel, f"{path}.res")
    return probs


def _rel_problems(r: np.ndarray, a: np.ndarray, rel: float, path: str) -> list[str]:
    if r.shape != a.shape:
        return [f"{path}: {r.size} values against {a.size}"]
    fin = np.isfinite(r) & np.isfinite(a)
    if not np.array_equal(np.isfinite(r), np.isfinite(a)) or not np.array_equal(
        np.isnan(r), np.isnan(a)
    ):
        return [f"{path}: Inf or NaN at other places"]
    if not np.array_equal(r[~fin], a[~fin]):
        return [f"{path}: other infinities"]
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.abs(r[fin] - a[fin]) / np.abs(r[fin])
    if d.size and np.nanmax(d) > rel:
        return [f"{path}: relative difference up to {np.nanmax(d):.3g} (limit {rel:g})"]
    return []


RES = Custom(
    lambda r, a: _res_problems(r, a, RES_STEPS, RES_REL),
    f"resolution: the radii within {RES_STEPS} step of 1e-4 (proposed, docs/parity/tools.md)",
)
RES_VALUES = Custom(
    lambda r, a: _rel_problems(_data(r), _data(a), RES_REL, "res"),
    f"resolution values within {RES_REL:g} relative (proposed, docs/parity/tools.md)",
)


def _compare(be: str, fsc_out: tuple[int, ...], res_out: tuple[int, ...], other=None) -> dict:
    out: dict = {}
    for k in fsc_out:
        if be == "permissive":
            out[k] = PERMISSIVE
    for k in res_out:
        out[k] = RES
    out.update(other or {})
    return out


# ---------------------------------------------------------------------------------------------
# data (exact arithmetic: the same bits everywhere); MATLAB's form
# ---------------------------------------------------------------------------------------------


def normal_pair(n: int, seed: str):
    """FSCAnalysisTest's a = randn(n, n, n), b = a + randn(n, n, n) (CounterRNG)"""
    rng = CounterRNG(0, f"b05/{seed}")
    a = rng.normal((n, n, n))
    b = a + rng.spawn("b").normal((n, n, n))
    return to_matlab(a), to_matlab(b)


def fsc_image(sz_zyx, cls: str, name: str) -> np.ndarray:
    """``cast(1000 + 300 * texture(sz, 1.5, [0, 1]) + 30 * randn(sz), cls)`` (MATLAB's form)"""
    rng = CounterRNG(0, f"b05/img/{name}")
    t = synth.texture(sz_zyx, 1.5, (0.0, 1.0), rng=rng).astype(np.float64)
    v = 1000.0 + 300.0 * t + 30.0 * rng.spawn("noise").normal(sz_zyx)
    return to_matlab(np.asarray(dtypes.cast(v, cls)))


def iterations(folder, n: int = 64, shape=None) -> None:
    """noisierIterationHasWorseResolution's data: smooth structure with little noise, the same
    with much noise, and noise alone, as Iter_0002/4/6.tif (single)"""
    rng = CounterRNG(0, "b05/iters")
    sh = shape or (n, n, n)
    b = synth.gaussian_filter(rng.normal(sh), 2.0)
    s = 1000.0 + 300.0 * b / np.sqrt(np.mean((b - b.mean()) ** 2) * b.size / (b.size - 1))
    noise = rng.spawn("noise").normal(sh)
    for name, v in (
        ("Iter_0002.tif", s + 20 * noise),
        ("Iter_0004.tif", s + 600 * noise),
        ("Iter_0006.tif", 1000.0 + 600 * noise),
    ):
        synth.write_tiff(os.path.join(folder, name), v.astype(np.float32))


def smooth_iteration(folder, n: int = 48, name: str = "Iter_0002.tif") -> None:
    """existingResultsAreNotRecomputed's image: 1000 + 300 * a smoothed normal volume"""
    rng = CounterRNG(0, "b05/smooth")
    v = 1000.0 + 300.0 * synth.gaussian_filter(rng.normal((n, n, n)), 1.0)
    synth.write_tiff(os.path.join(folder, name), v.astype(np.float32))


# ---------------------------------------------------------------------------------------------
# the Python side
# ---------------------------------------------------------------------------------------------


def _fsc_m(f: dict) -> dict:
    return FS._fsc_struct(f)


def _res_m(r: dict) -> dict:
    return FS._res_struct(r)


def py_hamming(n):
    return FS.hamming(int(n)).reshape(-1, 1)


def _with_backend(be, fn, *args):
    if be is None:
        return fn(*args)
    with F.using(be):
        return fn(*args)


def py_xr_fsc(*args, be=None):
    from petakit5d import matlab

    f, r = _with_backend(be, matlab.XR_FSC, *args)
    return [_fsc_m(f), _res_m(r)]


def py_xr_fsc_resolution(*args):
    from petakit5d import matlab

    return _res_m(matlab.XR_FSC_resolution(*args))


def py_one_image(*args, be=None):
    from petakit5d import matlab

    f, r, fc, rc = _with_backend(be, matlab.XR_one_image_FSC_analysis, *args)
    return [_fsc_m(f), _res_m(r), mcell(*[_fsc_m(x) for x in fc]), mcell(*[_res_m(x) for x in rc])]


def py_matlab_call(name, *args, be=None):
    from petakit5d import matlab

    _with_backend(be, getattr(matlab, name), *args)
    return []


def _load(path):
    from petakit5d.io import read_mat

    return read_mat(path)


def _mat_cell(v):
    """a cell read back from a MAT-file (a list, or an object array) as MATLAB's cell"""
    if isinstance(v, list):
        return mcell(*v)
    a = np.asarray(v, dtype=object)
    items = [mcell(*x) if isinstance(x, list) else to_mvalue(x) for x in a.ravel(order="F")]
    return MCell(a.shape, items)


def py_frame_loaded(fn, fnout, *args, be=None):
    from petakit5d import matlab

    _with_backend(be, matlab.XR_one_image_FSC_analysis_frame, fn, fnout, *args)
    a = _load(fnout)
    out = [a["fsc_mu"], a["res_mu"], a["fsc_cell"], a["res_cell"], a["centers"], a["bbox"]]
    return [
        _mat_cell(v) if isinstance(v, (list, np.ndarray)) and k < 4 else v
        for k, v in enumerate(out)
    ]


def py_wrapper_loaded(path, *args, be=None):
    from petakit5d import matlab

    _with_backend(be, matlab.XR_FSC_analysis_wrapper, [path], *args)
    a = _load(path + "FSCs/Iter_0002.mat")
    info = _load(path + "FSCs/fsc_res_info_c1.mat")
    return [
        a["fsc_mu"],
        a["res_mu"],
        _mat_cell(a["fsc_cell"]),
        _mat_cell(a["res_cell"]),
        a["centers"],
        a["bbox"],
        info["iter"],
        info["thetas"],
        info["res_mu_mat"],
        info["res_mat"],
        _mat_cell(info["fsns_c"]),
    ]


def py_wrapper_twice(path, *args):
    py_matlab_call("XR_FSC_analysis_wrapper", [path], *args)
    py_matlab_call("XR_FSC_analysis_wrapper", [path], *args)
    return []


def py_wrapper_regions(path, *args, be=None):
    from petakit5d import matlab

    _with_backend(be, matlab.XR_FSC_analysis_wrapper, [path], *args)
    info = _load(path + "FSCs/fsc_res_info_c1.mat")
    a = _load(path + "FSCs/Iter_0002.mat")
    return [info["res_mat"], info["res_mu_mat"], a["centers"], a["bbox"]]


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

HAMMING_NS = (
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    15,
    16,
    17,
    20,
    21,
    24,
    25,
    31,
    32,
    50,
    51,
    64,
    99,
    100,
    101,
    128,
    151,
    200,
    201,
    250,
    251,
    256,
    400,
    401,
    500,
    501,
    512,
)
CASES: list[Case] = [
    Case(f"hamming.{n}", matlab="hamming", args=[n], python=py_hamming) for n in HAMMING_NS
]


def fsc_cases(cid, make, nv, *, error=False, backends=BACKENDS, **kw):
    """XR_FSC(a, b, nv{:}) on each backend"""
    out = []
    for be in (None,) if error else backends:
        out.append(
            Case(
                f"{cid}.{be}" if be else cid,
                matlab="XR_FSC",
                args=[In("a"), In("b"), *nv],
                nargout=2,
                python=py_xr_fsc,
                py_kwargs={"be": be},
                inputs={"a": Value(lambda m=make: m()[0]), "b": Value(lambda m=make: m()[1])},
                compare={} if error else _compare(be, (0,), (1,)),
                expect_error=error,
                **({"py_error": (FS.FSCError,)} if error else {}),
                **kw,
            )
        )
    return out


H21 = ["halfSize", [21, 21, 21]]
CASES += fsc_cases("fsc.cache.dr2", lambda: normal_pair(21, "a"), ["dr", 2, *H21])
CASES += fsc_cases("fsc.cache.dr3", lambda: normal_pair(21, "a"), ["dr", 3, *H21])
CASES += fsc_cases(
    "fsc.cache.dtheta4", lambda: normal_pair(21, "a"), ["dr", 2, "dtheta", PI / 4, *H21]
)
CASES += fsc_cases(
    "fsc.cache.no_cone", lambda: normal_pair(21, "a"), ["dr", 2, "skipConeRegion", False, *H21]
)
CASES += fsc_cases(
    "fsc.cache.size25", lambda: normal_pair(25, "c"), ["dr", 2, "halfSize", [25, 25, 25]]
)
CASES += fsc_cases("fsc.defaults21", lambda: normal_pair(21, "d"), H21)
CASES += fsc_cases(
    "fsc.even20", lambda: normal_pair(20, "e"), ["dr", 2, "dtheta", PI / 6, "halfSize", [20] * 3]
)
CASES += fsc_cases(
    "fsc.pad_crop",
    lambda: (normal_pair(23, "f")[0][1:20, :, 2:21], normal_pair(23, "f")[1][1:20, :, 2:21]),
    ["dr", 2, "dtheta", PI / 6, *H21],
)
CASES += fsc_cases(
    "fsc.one_bit",
    lambda: normal_pair(25, "g"),
    ["dr", 2, "dtheta", PI / 6, "resThreshMethod", "one-bit", "halfSize", [25] * 3],
)
CASES += fsc_cases(
    "fsc.half_bit",
    lambda: normal_pair(25, "g"),
    ["dr", 2, "dtheta", PI / 6, "resThreshMethod", "half-bit", "halfSize", [25] * 3],
)
CASES += fsc_cases(
    "fsc.fixed_0p5",
    lambda: normal_pair(25, "g"),
    [
        "dr",
        1,
        "dtheta",
        PI / 6,
        "resThresh",
        0.5,
        "xyPixelSize",
        0.1,
        "dz",
        0.3,
        "halfSize",
        [25] * 3,
    ],
)
CASES += fsc_cases(
    "fsc.odd_angles",
    lambda: normal_pair(21, "a"),
    ["dr", 2, "dtheta", PI / 5, *H21],
    warnings="ids",
)
# pi/24: the sectors at 0 and pi/2 lie inside the cones, so two angles are missing (npoints
# filled by the angles' positions)
CASES += fsc_cases(
    "fsc.missing_sectors",
    lambda: normal_pair(21, "a"),
    ["dr", 2, "dtheta", PI / 24, *H21],
    warnings="ids",
)
CASES += fsc_cases(
    "fsc.err.zero",
    lambda: (np.zeros((21, 21, 21)), np.zeros((21, 21, 21))),
    ["dr", 2, *H21],
    error=True,
)
CASES += fsc_cases(
    "fsc.err.not_cubic",
    lambda: (normal_pair(21, "a")[0][:, :, :15], normal_pair(21, "a")[1][:, :, :15]),
    ["dr", 2, "halfSize", [21, 21, 15]],
    error=True,
)

# FSCAnalysisTest's image cases: {MATLAB size [y x z], class, xy, dz, halfSize, bbox, resAxis,
# clipPer}
FSC_IMAGES = {
    "zUp": ((60, 57, 23), "uint16", 0.108, 0.3, (21, 21, 21), None, "xz", None),
    "xyUp": ((30, 28, 50), "uint16", 0.108, 0.05, (21, 21, 21), None, "xz", None),
    "padZ": ((60, 57, 9), "uint16", 0.108, 0.3, (21, 21, 21), None, "xz", None),
    "padAll": ((20, 22, 9), "uint16", 0.108, 0.3, (21, 21, 21), None, "xz", None),
    "isotropic": ((50, 47, 45), "single", 0.1, 0.1, (21, 21, 21), None, "xz", None),
    "evenHalfSize": ((60, 57, 23), "uint16", 0.108, 0.3, (20, 20, 20), None, "xz", None),
    "yzBbox": ((60, 57, 23), "uint16", 0.108, 0.2, (15, 15, 15), (5, 3, 2, 55, 50, 21), "yz", None),
    "xyClip": ((40, 37, 35), "double", 0.108, 0.15, (17, 17, 17), None, "xy", 99),
}


def _image_of(name):
    msz, cls = FSC_IMAGES[name][:2]
    return fsc_image((msz[2], msz[0], msz[1]), cls, name)


def one_image_cases(name, extra=(), cid=None, dtheta=PI / 6):
    msz, cls, xy, dz, h, bbox, ax, cp = FSC_IMAGES[name]
    nv = [
        "xyPixelSize",
        xy,
        "dz",
        dz,
        "dr",
        2,
        "dtheta",
        dtheta,
        "halfSize",
        list(h),
        "bbox",
        list(bbox) if bbox else None,
        "resAxis",
        ax,
        "clipPer",
        cp,
        *extra,
    ]
    out = []
    for be in BACKENDS:
        out.append(
            Case(
                f"{cid or 'one_image.' + name}.{be}",
                matlab="XR_one_image_FSC_analysis",
                args=[In("im"), *nv],
                nargout=4,
                python=py_one_image,
                py_kwargs={"be": be},
                inputs={"im": Value(lambda n=name: _image_of(n))},
                compare=_compare(be, (0, 2), (1, 3)),
                warnings="ids" if dtheta != PI / 6 else "record",
            )
        )
    return out


for _name in FSC_IMAGES:
    CASES += one_image_cases(_name)
CASES += one_image_cases("zUp", ["resThreshMethod", "one-bit"], "one_image.zUp_one_bit")
CASES += one_image_cases("zUp", ["resThreshMethod", "half-bit"], "one_image.zUp_half_bit")
CASES += one_image_cases("xyUp", ["skipConeRegion", False], "one_image.xyUp_no_cone")
CASES += one_image_cases("zUp", (), "one_image.zUp_odd_angles", PI / 5)


def _tif(name):
    def make(path):
        synth.write_tiff(path, from_matlab(_image_of(name)))

    return make


FRAME_CODE = (
    "XR_one_image_FSC_analysis_frame(in{1}, in{2}, in{3 : end}); a = load(in{2}); "
    "out = {a.fsc_mu, a.res_mu, a.fsc_cell, a.res_cell, a.centers, a.bbox};"
)


def frame_cases(cid, image, nv):
    out = []
    for be in BACKENDS:
        out.append(
            Case(
                f"{cid}.{be}",
                code=FRAME_CODE,
                args=[InPath("im"), Work("r.mat"), *nv],
                nargout=6,
                python=py_frame_loaded,
                py_kwargs={"be": be},
                inputs={"im": File(_tif(image), "im.tif")},
                tree=tree(per_file={"*.mat": ignore("compared as loaded outputs")}),
                compare=_compare(be, (0, 2), (1, 3)),
            )
        )
    return out


FRAME_NV = ["xyPixelSize", 0.108, "dz", 0.3, "dr", 2, "dtheta", PI / 6, "halfSize", [15, 15, 15]]
CASES += frame_cases("frame.whole", "zUp", FRAME_NV)
CASES += frame_cases("frame.one_box", "zUp", [*FRAME_NV, "bbox", [5, 3, 2, 55, 50, 21]])
CASES += frame_cases(
    "frame.two_boxes",
    "zUp",
    [*FRAME_NV, "bbox", np.array([[5, 3, 2, 55, 50, 21], [1, 1, 1, 40, 40, 23]], float)],
)
CASES.append(
    Case(
        "frame.exists_skips",
        matlab="XR_one_image_FSC_analysis_frame",
        nargout=0,
        args=[InPath("im"), InPath("old"), *FRAME_NV],
        python=lambda *a: py_matlab_call("XR_one_image_FSC_analysis_frame", *a),
        inputs={
            "im": File(_tif("zUp"), "im.tif"),
            "old": File(lambda p: open(p, "wb").close(), "old.mat"),
        },
    )
)

WRAPPER_CODE = (
    "XR_FSC_analysis_wrapper({in{1}}, in{2 : end}); a = load([in{1}, 'FSCs/Iter_0002.mat']); "
    "info = load([in{1}, 'FSCs/fsc_res_info_c1.mat']); "
    "out = {a.fsc_mu, a.res_mu, a.fsc_cell, a.res_cell, a.centers, a.bbox, info.iter, "
    "info.thetas, info.res_mu_mat, info.res_mat, info.fsns_c};"
)
WRAPPER_TREE = tree(
    # the figures (MATLAB's PNG and .fig; Python's PNG only with matplotlib: else an empty folder)
    ignore=[r".*figures(/.*)?$"],
    per_file={"*.mat": ignore("compared as loaded outputs")},
)
NOISIER_NV = [
    "resultDirName",
    "FSCs",
    "xyPixelSize",
    0.1,
    "dz",
    0.1,
    "dr",
    2,
    "dtheta",
    PI / 6,
    "resThreshMethod",
    "fixed",
    "resThresh",
    0.2,
    "halfSize",
    [32, 32, 32],
    "channelPatterns",
    ["Iter"],
    "channels",
    488,
    "iterInterval",
    2,
    "parseCluster",
    False,
    "masterCompute",
    True,
]


def wrapper_cases(
    cid,
    make,
    nv,
    code=WRAPPER_CODE,
    nargout=11,
    python=py_wrapper_loaded,
    res_values=(1, 3),
    res_arrays=(8, 9),
    fsc_out=(0, 2),
):
    out = []
    for be in BACKENDS:
        other = dict.fromkeys(res_arrays, RES_VALUES)
        out.append(
            Case(
                f"{cid}.{be}",
                code=code,
                args=[InPath("d"), *nv],
                nargout=nargout,
                python=python,
                py_kwargs={"be": be},
                inputs={"d": Tree(make, "d")},
                tree=WRAPPER_TREE,
                compare=_compare(be, fsc_out, res_values, other),
            )
        )
    return out


CASES += wrapper_cases("wrapper.noisier", iterations, NOISIER_NV)
SMALL_NV = [
    "xyPixelSize",
    0.1,
    "dz",
    0.1,
    "dr",
    2,
    "halfSize",
    [24, 24, 24],
    "channelPatterns",
    ["Iter"],
    "channels",
    488,
    "parseCluster",
    False,
    "masterCompute",
    True,
]
CASES.append(
    Case(
        "wrapper.rerun_warns",
        code="XR_FSC_analysis_wrapper({in{1}}, in{2 : end}); "
        "XR_FSC_analysis_wrapper({in{1}}, in{2 : end}); out = {};",
        args=[InPath("d"), *SMALL_NV],
        nargout=0,
        python=py_wrapper_twice,
        inputs={"d": Tree(smooth_iteration, "d")},
        tree=WRAPPER_TREE,
        warnings="ids",
    )
)
REGION_CODE = (
    "XR_FSC_analysis_wrapper({in{1}}, in{2 : end}); "
    "info = load([in{1}, 'FSCs/fsc_res_info_c1.mat']); a = load([in{1}, 'FSCs/Iter_0002.mat']); "
    "out = {info.res_mat, info.res_mu_mat, a.centers, a.bbox};"
)
REGION_NV = [
    "xyPixelSize",
    0.1,
    "dz",
    0.1,
    "dr",
    2,
    "dtheta",
    PI / 6,
    "halfSize",
    [12, 12, 12],
    "channelPatterns",
    ["Iter"],
    "channels",
    488,
    "multiRegions",
    True,
    "parseCluster",
    False,
    "masterCompute",
    True,
]
CASES += wrapper_cases(
    "wrapper.regions.grid",
    lambda d: iterations(d, shape=(40, 48, 44)),
    [*REGION_NV, "regionGrid", np.array([[24, 22, 20], [13, 13, 13], [30, 30, 25]], float)],
    code=REGION_CODE,
    nargout=4,
    python=py_wrapper_regions,
    res_values=(),
    res_arrays=(0, 1),
    fsc_out=(),
)
CASES += wrapper_cases(
    "wrapper.regions.interval",
    lambda d: iterations(d, shape=(40, 48, 44)),
    [*REGION_NV, "regionInterval", [10, 12, -1]],
    code=REGION_CODE,
    nargout=4,
    python=py_wrapper_regions,
    res_values=(),
    res_arrays=(0, 1),
    fsc_out=(),
)
CASES.append(
    Case(
        "wrapper.regions.single_region_fails",
        code=REGION_CODE,
        args=[InPath("d"), *REGION_NV, "regionGrid", [24, 22, 20]],
        nargout=4,
        python=py_wrapper_regions,
        inputs={"d": Tree(lambda d: iterations(d, shape=(40, 48, 44)), "d")},
        tree=WRAPPER_TREE,
        expect_error=True,
        py_error=(FS.FSCError,),
    )
)
CASES.append(
    Case(
        "two_image.fails",
        matlab="XR_two_image_FSC_analysis",
        args=[In("a"), In("b"), "dz", 0.2, "N", [10, 10, 10]],
        nargout=2,
        python=lambda *a: py_matlab_call("XR_two_image_FSC_analysis", *a),
        inputs={"a": Value(lambda: _image_of("padAll")), "b": Value(lambda: _image_of("padAll"))},
        expect_error="MATLAB:InputParser:UnmatchedParameter",
        py_error=(FS.FSCError,),
        error_message=True,
    )
)


def _check(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    _check(case, tmp_path)


_KW = {
    "xyPixelSize": "xy_pixel_size",
    "dz": "dz",
    "resThreshMethod": "res_thresh_method",
    "resThresh": "res_thresh",
}
_NAMES = {"xyPixelSize", "dz", "dr", "dtheta", "halfSize", "resultDirName", "resThreshMethod"}


def _method_kw(case) -> dict:
    args = list(case.args)
    k0 = next(k for k, a in enumerate(args) if isinstance(a, str) and a in _NAMES)
    nv = args[k0:]
    return {_KW[n]: v for n, v in zip(nv[::2], nv[1::2], strict=True) if n in _KW}


def _pairs_of(case) -> list[tuple[int, int, dict]]:
    """(FSC output, resolution output, the resolution's options) of a case; the pairs' own
    resolutions (fsc_cell, res_cell) are XR_FSC's: the pixel sizes, 'fixed' at 0.2"""
    kw = _method_kw(case)
    if case.matlab == "XR_FSC":
        return [(0, 1, kw)]
    sizes = {k: v for k, v in kw.items() if k in ("xy_pixel_size", "dz")}
    return [(0, 1, kw), (2, 3, sizes)]


def _fsc_of(fm: dict) -> dict:
    fsc = {k: np.asarray(fm[k].data, dtype=np.float64) for k in ("fsc", "npoints")}
    fsc.update({k: _data(fm[k]) for k in ("r", "thetas")})
    return fsc


SPLINE_CASES = [
    c
    for c in CASES
    if c.id.endswith(".permissive")
    and c.id.startswith(("fsc.", "one_image.", "frame.", "wrapper.noisier"))
]


@pytest.mark.parametrize("case", SPLINE_CASES, ids=lambda c: c.id[: -len(".permissive")])
def test_resolution_of_matlabs_fsc(case):
    # XR_FSC_resolution of MATLAB's own FSC curves (no FFT difference): only the smoothing
    # spline differs, and the radii are MATLAB's (docs/parity/tools.md, proposed)
    ref = oracle.reference(case)
    checked, probs = 0, []
    for i_f, i_r, kw in _pairs_of(case):
        for (path, fm), (_, rm) in zip(
            _structs(ref.outputs[i_f], ""), _structs(ref.outputs[i_r], ""), strict=True
        ):
            if any(v.data is None for v in (*fm.values(), *rm.values())):
                continue
            mine = to_mvalue(_res_m(FS.fsc_resolution(_fsc_of(fm), **kw)))
            ref_res = MStruct((1, 1), sorted(rm), [rm])
            found = _res_problems(ref_res, mine, RES_STEPS, RES_REL)
            probs += [f"out{i_r + 1}{path}{p}" for p in found]
            checked += 1
    if not checked:
        pytest.skip("MATLAB's curves are not in the manifest (hash only)")
    assert not probs, probs


# ---------------------------------------------------------------------------------------------
# FSCAnalysisTest, ported
# ---------------------------------------------------------------------------------------------


def _wrapper_kw(**kw):
    base = {
        "xy_pixel_size": 0.1,
        "dz": 0.1,
        "dr": 2,
        "channel_patterns": ("Iter",),
        "channels": (488,),
        "parse_cluster": False,
        "master_compute": True,
    }
    base.update(kw)
    return base


def test_noisier_iteration_has_worse_resolution(tmp_path):
    d = str(tmp_path) + "/"
    iterations(d)
    FS.fsc_analysis_dataset(
        [d],
        **_wrapper_kw(dtheta=PI / 6, half_size=(32, 32, 32), iter_interval=2),
    )
    from petakit5d.io import read_mat

    for name in ("Iter_0002", "Iter_0004", "Iter_0006"):
        a = read_mat(f"{d}FSCs/{name}.mat")
        assert {"fsc_mu", "res_mu", "fsc_cell", "res_cell"} <= set(a)
        assert len(a["fsc_cell"]) == 4
    info = read_mat(f"{d}FSCs/fsc_res_info_c1.mat")
    ntheta = round(2 * PI / (PI / 6)) + 1
    res = np.array(info["res_mat"], dtype=np.float64)
    assert res.shape == (3, ntheta)
    res[~np.isfinite(res)] = np.nan
    mean_res = np.nanmean(res, axis=1)
    assert mean_res[0] < mean_res[1], "More noise should give a worse resolution."
    assert mean_res[1] < mean_res[2], "Pure noise should give the worst resolution."
    # no better than the voxel size allows (Nyquist on the diagonal)
    assert np.all(mean_res >= 2 * 0.1 / math.sqrt(2) - 1e-6)
    if HAS_MPL:
        figs = [f for f in os.listdir(d + "figures") if f.endswith(".png")]
        assert len(figs) >= 4


def test_existing_results_are_not_recomputed(tmp_path):
    d = str(tmp_path) + "/"
    smooth_iteration(d)
    kw = _wrapper_kw(half_size=(24, 24, 24))
    FS.fsc_analysis_dataset([d], **kw)
    fn = d + "FSCs/Iter_0002.mat"
    first = os.stat(fn).st_mtime_ns
    with pytest.warns(FS.FSCWarning, match="Directory already exists"):
        FS.fsc_analysis_dataset([d], **kw)
    assert os.stat(fn).st_mtime_ns == first


def test_result_is_moved_into_place(tmp_path, monkeypatch):
    # the job framework counts a task as done once its result exists: the result is written
    # under a temporary name and renamed
    d = str(tmp_path) + "/"
    smooth_iteration(d)
    moved = []
    real = os.replace

    def spy(src, dst, *a, **k):
        moved.append((os.fspath(src), os.fspath(dst)))
        return real(src, dst, *a, **k)

    monkeypatch.setattr(os, "replace", spy)
    FS.fsc_analysis_dataset([d], **_wrapper_kw(half_size=(24, 24, 24)))
    target = os.path.normpath(d + "FSCs/Iter_0002.mat")
    assert any(os.path.normpath(dst) == target for _, dst in moved), "saved in place"
    left = [f for f in os.listdir(d + "FSCs") if f.startswith("Iter_0002_")]
    assert not left, "a temporary file was left behind"


def _reference_fsc(im, xy, dz, h, bbox, res_axis, clip_per):
    """XR_one_image_FSC_analysis before only the kept region of the resize was computed: resize
    the whole image, crop, pad (FSCAnalysisTest's referenceFSC, (z, y, x))"""
    from petakit5d.core.dtypes import mround
    from petakit5d.matlabcompat.resize import imresize3

    im = np.asarray(im).astype(np.float64)
    if clip_per is not None:
        im = FS._clip(im, clip_per)
    if bbox is not None:
        z0, y0, x0, z1, y1, x1 = bbox
        im = im[z0:z1, y0:y1, x0:x1]
    if res_axis == "yz":
        im = np.swapaxes(im, 1, 2)
    elif res_axis == "xy":
        im = np.transpose(im, (1, 0, 2))
    im = np.ascontiguousarray(im)
    if dz != xy:
        f = (dz, xy, xy)
        osz = [int(mround(s * k / min(xy, dz))) for s, k in zip(im.shape, f, strict=True)]
        im = imresize3(im, osz, "linear")
    sz = im.shape
    if any(s > 2 * hv for s, hv in zip(sz, h, strict=True)):
        s0 = [max(1, math.floor((s - 2 * hv) / 2)) for s, hv in zip(sz, h, strict=True)]
        t0 = [min(s, a + 2 * hv - 1) for s, a, hv in zip(sz, s0, h, strict=True)]
        im = im[s0[0] - 1 : t0[0], s0[1] - 1 : t0[1], s0[2] - 1 : t0[2]]
        sz = im.shape
    if any(s < 2 * hv for s, hv in zip(sz, h, strict=True)):
        lo = [math.floor((2 * hv - s) / 2) for s, hv in zip(sz, h, strict=True)]
        im = np.pad(im, [(a, 2 * hv - s - a) for a, s, hv in zip(lo, sz, h, strict=True)])
    im = np.ascontiguousarray(im)
    kw = {"xy_pixel_size": xy, "dz": dz, "dr": 2, "dtheta": PI / 6, "half_size": h}
    fscs = [
        FS.fourier_shell_correlation(FS._sub_lattice(im, a), FS._sub_lattice(im, b), **kw)[0]
        for a, b in FS._PAIRS
    ]
    mu = dict(fscs[0])
    mu["fsc"] = (((fscs[0]["fsc"] + fscs[1]["fsc"]) + fscs[2]["fsc"]) + fscs[3]["fsc"]) / 4
    return mu, FS.fsc_resolution(mu, xy_pixel_size=xy, dz=dz)


def _assert_tree_equal(a, b):
    if isinstance(a, dict):
        assert sorted(a) == sorted(b)
        for k in a:
            _assert_tree_equal(a[k], b[k])
    else:
        x, y = np.asarray(a), np.asarray(b)
        assert x.shape == y.shape and x.dtype == y.dtype
        assert np.array_equal(x.view(np.uint8), y.view(np.uint8))


@pytest.mark.parametrize("name", list(FSC_IMAGES))
def test_one_image_matches_resizing_the_whole_image(name, monkeypatch):
    msz, cls, xy, dz, h, bbox, res_axis, clip_per = FSC_IMAGES[name]
    im = from_matlab(_image_of(name))
    half = (h[2], h[0], h[1])
    box = None
    if bbox:
        from petakit5d.core.bbox import box_from_matlab

        box = box_from_matlab(bbox)
    fsc_ref, res_ref = _reference_fsc(im, xy, dz, half, box, res_axis, clip_per)
    # the resize must not fall back to imresize3 of the whole image
    from petakit5d.matlabcompat import resize

    def no_whole(*a, **k):
        raise AssertionError("imresize3 of the whole image")

    monkeypatch.setattr(resize, "imresize3", no_whole)
    fsc_mu, res_mu, _, _ = FS.one_image_fsc_analysis(
        im,
        xy_pixel_size=xy,
        dz=dz,
        dr=2,
        dtheta=PI / 6,
        half_size=half,
        bbox=box,
        res_axis=res_axis,
        clip_per=clip_per,
    )
    _assert_tree_equal(fsc_mu, fsc_ref)
    _assert_tree_equal(res_mu, res_ref)


def test_shell_cache_follows_parameters():
    # the window and the frequency shells are kept between calls: the same call twice and calls
    # with another dr, dtheta, skip_cone_region or half_size give the results of fresh calls
    a, b = (from_matlab(v) for v in normal_pair(21, "a"))
    c, d = (from_matlab(v) for v in normal_pair(25, "c"))
    h = {"half_size": (21, 21, 21)}
    calls = [
        (a, b, {"dr": 2, **h}),
        (a, b, {"dr": 2, **h}),
        (a, b, {"dr": 3, **h}),
        (a, b, {"dr": 2, "dtheta": PI / 4, **h}),
        (a, b, {"dr": 2, "skip_cone_region": False, **h}),
        (c, d, {"dr": 2, "half_size": (25, 25, 25)}),
        (a, b, {"dr": 2, **h}),
    ]
    FS.clear_caches()
    kept = [FS.fourier_shell_correlation(x, y, **kw) for x, y, kw in calls]
    for (x, y, kw), got in zip(calls, kept, strict=True):
        FS.clear_caches()
        fresh = FS.fourier_shell_correlation(x, y, **kw)
        _assert_tree_equal(got[0], fresh[0])
        _assert_tree_equal(got[1], fresh[1])


# ---------------------------------------------------------------------------------------------
# the pieces
# ---------------------------------------------------------------------------------------------


def test_colon_matches_the_literal_parser():
    from petakit5d.params.literal import parse

    for a, d, b in (
        (0, 1e-4, 1),
        (0, 1e-4, math.sqrt(2)),
        (0.04, 1e-4, 1),
        (0, PI / 6, 2 * PI),
        (0, PI / 4, 2 * PI),
        (-9.5, 1, 9.5),
        (24, -10, 13),
        (5, 2, 4),
        (0.121, 0.1, 18.721),
        (0.118, 0.3, 36.718),
        (0.892, 0.1, 15.292),
        (0.606, 0.0001, 0.621),
    ):
        got = FS._colon(a, d, b)
        ref = parse(f"{a!r} : {d!r} : {b!r}")
        want = np.asarray(ref.data if hasattr(ref, "data") else ref, dtype=np.float64).ravel()
        assert np.array_equal(got.view(np.uint64), want.view(np.uint64)), (a, d, b)


def test_abs_squared_equals_matlab_hypot():
    from petakit5d.matlabcompat.elementary import hypot

    rng = CounterRNG(0, "b05/hypot")
    u = rng.uniform((4, 20000))
    re = (u[0] - 0.5) * np.exp2(np.floor(u[1] * 200 - 100))
    im = (u[2] - 0.5) * np.exp2(np.floor(u[3] * 200 - 100))
    re[:50] = 0.0
    im[50:100] = -0.0
    re[100:110] = np.inf
    im[110:120] = np.nan
    re[120:130] = 5e-310
    im[130:140] = 1e300
    re[140:150] = im[140:150] * (1 + 2**-40)
    re[150:160] = 2 * im[150:160]  # a - b == b: msun's choice between its two formulas
    h = hypot(re, im)
    with np.errstate(over="ignore"):
        want = h * h
    got = FS._abs_squared(re, im)
    assert np.array_equal(np.isnan(got), np.isnan(want))
    ok = ~np.isnan(want)
    assert np.array_equal(got[ok].view(np.uint64), want[ok].view(np.uint64))


def test_smoothing_spline_matches_scipy():
    scipy_interp = pytest.importorskip("scipy.interpolate")
    if not hasattr(scipy_interp, "make_smoothing_spline"):
        pytest.skip("SciPy without make_smoothing_spline")
    rng = CounterRNG(0, "b05/spline")
    x = np.cumsum(0.05 + rng.uniform(30) * 0.1)
    y = np.cos(3 * x) + 0.1 * rng.normal(30)
    for p in (0.5, 0.9, 0.999, 1.0):
        f = FS.smoothing_spline(x, y, p)
        s = scipy_interp.make_smoothing_spline(x, y, lam=(1 - p) / p)
        xs = np.linspace(x[0], x[-1], 501)
        assert np.allclose(f(xs), s(xs), rtol=0, atol=1e-9), p


def test_default_smoothing_balances_the_traces():
    # p trace(R) = (1 - p) trace(Q'Q); for spacing h: trace(R) = m 2h/3, trace(Q'Q) = 6m/h^2,
    # so p = 1 / (1 + h^3 / 9), in the documented favorable range near 1 / (1 + h^3 / 6)
    for h, n in ((0.08, 13), (0.5, 7), (2.0, 30)):
        x = np.arange(n) * h
        assert math.isclose(FS.default_smoothing(x), 1 / (1 + h**3 / 9), rel_tol=1e-12)
    x = np.array([0.0, 0.1, 0.3, 0.35, 1.0, 2.0])
    p = FS.default_smoothing(x)
    hs = np.diff(x)
    ta = sum((hs[j] + hs[j + 1]) / 3 for j in range(4))
    q = [(1 / hs[j], 1 / hs[j] + 1 / hs[j + 1], 1 / hs[j + 1]) for j in range(4)]
    tb = sum(a * a + b * b + c * c for a, b, c in q)
    assert math.isclose(p * ta, (1 - p) * tb, rel_tol=1e-12)
    assert FS.default_smoothing([0.0, 1.0]) == 1.0


def test_smoothing_spline_limits():
    x = np.linspace(0, 1, 12)
    y = 0.3 + 2.0 * x
    # a line is its own smoothing spline, for any p; p = 1 interpolates
    for p in (0.1, 0.7, None):
        assert np.allclose(FS.smoothing_spline(x, y, p)(x), y, atol=1e-12)
    rng = CounterRNG(0, "b05/interp")
    z = rng.normal(12)
    assert np.allclose(FS.smoothing_spline(x, z, 1.0)(x), z, atol=1e-12)
    # beyond the last point: the last piece continued
    f = FS.smoothing_spline(x, z, 0.99)
    c = f.coefs[-1]
    t = 1.25 - x[-2]
    assert np.isclose(f(1.25), ((c[0] * t + c[1]) * t + c[2]) * t + c[3])


def test_poly2_fit_is_exact_least_squares():
    x = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    y = 1.5 - 0.5 * x + 2.0 * x * x
    f = FS.poly2_fit(x, y)
    assert f.coefficients == (2.0, -0.5, 1.5)
    rng = CounterRNG(0, "b05/poly")
    y2 = y + rng.normal(5)
    p = np.polyfit(x, y2, 2)
    f = FS.poly2_fit(x, y2)
    assert np.allclose(f.coefficients, p, atol=1e-12)


def test_pow_positive_equals_fdlibm_pow():
    from petakit5d.matlabcompat.elementary import power

    rng = CounterRNG(0, "b05/pow")
    u = rng.uniform((3, 4000))
    x = np.concatenate(
        [u[0], u[1] * 1e-3, np.exp2(np.floor(u[2] * 400 - 200)) * (1 + u[0]), [1.0, 0.5, 2.0]]
    )
    for y in (1 / 3, -1 / 3, 0.75, 3.5, -2.25):
        with np.errstate(over="ignore", under="ignore"):
            want = power(x, y)
            got = FS._pow_positive(x, y)
        assert np.array_equal(got.view(np.uint64), want.view(np.uint64)), y


def test_cell_sums_add_as_bincount():
    rng = CounterRNG(0, "b05/sums")
    cells = rng.integers(0, 40, 100000)
    v = rng.normal(100000) * np.exp2(np.floor(rng.uniform(100000) * 60 - 30))
    v[cells == 7] = -0.0
    order = np.argsort(cells, kind="stable")
    starts = np.concatenate([[0], np.cumsum(np.bincount(cells, minlength=41))])
    want = np.bincount(cells, weights=v, minlength=41)
    for nt in (1, 3):
        got = FS._cell_sums(v, order, starts, nt)
        assert np.array_equal(got.view(np.uint64), want.view(np.uint64))


@pytest.mark.parametrize("n", [15, 16])
@pytest.mark.parametrize("skip", [True, False])
@pytest.mark.parametrize("scale", ["normal", "zero", "one_zero", "tiny", "huge"])
def test_frc_without_the_last_cells_sums(n, skip, scale):
    # the last cell's value only decides whether it is NaN: the fast path equals the sums
    from petakit5d.matlabcompat.fft import fftn

    rng = CounterRNG(0, f"b05/frc/{n}")
    a, b = rng.normal((n, n, n)), rng.spawn("b").normal((n, n, n))
    if scale == "zero":
        a, b = a * 0, b * 0
    elif scale == "one_zero":
        b = b * 0
    elif scale == "tiny":
        a, b = a * 1e-200, b * 1e-200
    elif scale == "huge":
        a = a * 1e150
    cells = FS._shells(n, 2.0, PI / 6, skip, 2)
    f1, f2 = fftn(a), fftn(b)
    with np.errstate(over="ignore", invalid="ignore"):
        fast = FS._frc(f1, f2, cells, 2)
        full = FS._frc(f1, f2, cells, 2, fast=False)
    # the same cells are NaN; the others (but the last, whose value XR_FSC drops) are equal
    assert np.array_equal(np.isnan(fast), np.isnan(full))
    ok = ~np.isnan(full[:-1])
    assert np.array_equal(fast[:-1][ok].view(np.uint64), full[:-1][ok].view(np.uint64))
    if scale in ("zero", "one_zero"):
        assert np.isnan(full[-1])


@pytest.mark.parametrize("shape", [(7, 7, 7), (8, 8, 8), (6, 7, 5)])
def test_spectrum_values_use_the_hermitian_symmetry(shape):
    # the plane of -z gathered from that of z: the values of the whole shifted spectrum
    from petakit5d.matlabcompat.fft import fftn

    rng = CounterRNG(0, f"b05/spec/{shape}")
    f1, f2 = fftn(rng.normal(shape)), fftn(rng.spawn("b").normal(shape))
    got = FS._spectrum_values(f1, f2, 2)

    def order(v):
        return np.ascontiguousarray(np.swapaxes(np.fft.fftshift(v), 1, 2)).ravel()

    r1, i1, r2, i2 = (order(v) for v in (f1.real, f1.imag, f2.real, f2.imag))
    want = (r1 * r2 + i1 * i2, FS._abs_squared(r1, i1), FS._abs_squared(r2, i2))
    for g, w in zip(got, want, strict=True):
        assert np.array_equal(g.view(np.uint64), w.view(np.uint64))


def test_window_is_the_cube_root_of_the_hamming_products():
    from petakit5d.matlabcompat.elementary import power

    for n in (5, 20, 21):
        h = FS.hamming(n)
        want = power((h[None, :, None] * h[None, None, :]) * h[:, None, None], 1 / 3)
        assert np.array_equal(FS._window(n), want)


def test_odd_number_of_angles_rounds_the_indices():
    # [res(2.5 : 3), res, res(1 : 2.5)]: MATLAB warns twice and rounds the colon's values
    f = {
        "r": np.linspace(0, 1, 11),
        "thetas": np.array([0, 1, 2]) * PI / 5,
        "fsc": np.tile(np.linspace(1, 0, 11)[:, None], (1, 3)),
        "npoints": np.ones((11, 3)),
    }
    with pytest.warns(FS.FSCWarning, match="Integer operands") as rec:
        res = FS.fsc_resolution(f, xy_pixel_size=0.1, dz=0.3)
    assert [w.message.identifier for w in rec] == ["MATLAB:colon:nonIntegerIndex"] * 2
    r = (2 * (0.1 + 0.2 * np.abs(np.cos(f["thetas"])))) / res["points"]
    assert np.allclose(res["res"], np.concatenate([r[2:], r, r[:2]]))
    f["thetas"] = f["thetas"][:1]
    f["fsc"] = f["fsc"][:, :1]
    with pytest.raises(FS.FSCError) as e:
        FS.fsc_resolution(f)
    assert e.value.identifier == "MATLAB:badsubscript"
    f["thetas"] = np.array([0, 1]) * PI / 2
    f["fsc"] = np.tile(np.linspace(1, 0, 11)[:, None], (1, 2))
    f["fsc"] = f["fsc"][:, :2]
    res = FS.fsc_resolution(f, xy_pixel_size=0.1, dz=0.3)
    assert res["res"].size == 5
    # the curve 1 - r crosses 0.2 at 0.8 (the spline of a line is the line)
    assert np.allclose(res["points"], 0.8, atol=1e-4)


def test_unknown_method_fails():
    f = {
        "r": np.linspace(0, 1, 5),
        "thetas": np.array([0.0, 1.0]),
        "fsc": np.ones((5, 2)),
        "npoints": np.ones((5, 2)),
    }
    with pytest.raises(FS.FSCError) as e:
        FS.fsc_resolution(f, res_thresh_method="two-bit")
    assert e.value.identifier == "MATLAB:UndefinedFunction"


def test_results_do_not_depend_on_threads():
    a, b = (from_matlab(v) for v in normal_pair(25, "g"))
    kw = {"dr": 2, "dtheta": PI / 6, "half_size": (25, 25, 25)}
    one = FS.fourier_shell_correlation(a, b, num_threads=1, **kw)
    FS.clear_caches()
    four = FS.fourier_shell_correlation(a, b, num_threads=4, **kw)
    _assert_tree_equal(one[0], four[0])
    _assert_tree_equal(one[1], four[1])


def test_only_double_volumes():
    a = np.zeros((5, 5, 5), np.float32)
    with pytest.raises(NotImplementedError):
        FS.fourier_shell_correlation(a, a, half_size=(5, 5, 5))
