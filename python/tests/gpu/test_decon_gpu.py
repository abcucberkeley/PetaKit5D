"""The deconvolution on the GPU (petakit5d.decon.rl_gpu) against MATLAB's gpuArray code paths of
decon_lucy_function and decon_lucy_omw_function (useGPU true).

Ports tests/decon/DeconGPUTest.m (gpuMatchesCPU, gpuEdgeErosionUsesTheCPUMask,
gpuLoopMatchesThePreviousFormulation), the decon-level methods of DeconRFFTGPUTest.m
(gpuMexMatchesGpuMatlabLoop, gpuMexMatchesForOddAnd2DSizes, gpuMexAppliesTheOptions,
withoutTheGpuMexTheGpuLoopRuns) and DeconOverlapTest.m/gpuOverlapFunctionRunsWhileTheGPUComputes.

Oracle cases (GPU packs: MATLAB R2024b on the same GPU model; its cuFFT is 11.0.8, the port's pin):

- ``<gpu>.<method>.<path>...`` (gpu: A100, B200, TitanV; method rl or omw; path ``loop`` = useMex
  false, ``dflt`` = useMex [] (MATLAB's default, the loop), ``mex`` = useMex true):
  DeconRFFTGPUTest's data, sizes (odd, 2D, a single row) and options, DeconGPUTest's classes
  (uint16, single, double), sizes (smooth, prime factors) and option sets, backgrounds MATLAB's
  gpuArray code rounds, iterations 0 to 30, data with zero borders (EdgeErosion), the debug
  outputs of the GPU loop. Cases of another GPU model are skipped.
- ``<gpu>.sum.*``: MATLAB's gpuArray ``sum`` (single, and ``'double'``) of vectors of 1 to 10^7
  elements (``rl_gpu.matlab_gpu_sum``).
- ``any.*``: errors of MATLAB's GPU code paths, on any GPU.

Comparisons: bit for bit (MATLAB's gpuArray sums included: rl_gpu.matlab_gpu_sum). The native
half spectra (``layout='native'``, Python only) are within K09's accepted tolerance of the mex.

Tests without the ``gpu`` mark check what needs no GPU: the path choice, its errors, the boxes,
the memory sizes, the sum's block plan, the hand-off from petakit5d.decon.rl, the CPU fallbacks.
"""

from __future__ import annotations

import functools
import math
import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.decon import rl
from petakit5d.testing import oracle
from petakit5d.testing.cases import Box, Case, In, Value, Vol, Work
from petakit5d.testing.synth import CounterRNG, beads, decon_psf


def _rl_gpu():
    """petakit5d.decon.rl_gpu, imported when a test runs (not when the tests are collected: a test
    of petakit5d.decon.rl stands in for the module through sys.modules, which an imported module's
    attribute of its package would shadow)"""
    import importlib

    return importlib.import_module("petakit5d.decon.rl_gpu")


# the GPU models of the oracle cases: case id prefix -> part of MATLAB's (and CUDA's) GPU name
MODELS = {"A100": "A100", "B200": "B200", "TitanV": "TITAN V"}


def V(make):
    return Value(make)


# ---------------------------------------------------------------------------------------------
# inputs, in MATLAB's form (R03's and K09's generators: DeconRFFTTest's blurred beads, the PSFs)
# ---------------------------------------------------------------------------------------------

BLUR_SIGMA = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)


def size_name(s) -> str:
    return "x".join(str(v) for v in s)


def _py_size(msize) -> tuple[int, ...]:
    sz = tuple(int(v) for v in msize)
    return (sz[2], sz[0], sz[1]) if len(sz) == 3 else sz


def blurred_beads(msize, seed: int, cls: str = "uint16") -> np.ndarray:
    """beads blurred by the test PSF plus a little noise on a background of 100 (uint16); single:
    plus 0.37; double: the uint16 values; uint8: MATLAB's uint8(im / 8); MATLAB's form"""
    py = _py_size(msize)
    rng = CounterRNG(seed, "beads")
    n = max(3, round(math.prod(py) / 4000))
    pos = rng.uniform((n, len(py))) * (np.array(py, dtype=np.float64) - 1)
    im = beads(py, pos, BLUR_SIGMA if len(py) == 3 else BLUR_SIGMA[1], 1000.0, 100.0)
    im = im + np.float32(3) * rng.uniform(py, np.float32)
    u = np.floor(im + np.float32(0.5)).astype(np.uint16)  # im >= 0
    if cls == "uint16":
        out = u
    elif cls == "uint8":
        out = np.floor(u.astype(np.float64) / 8 + 0.5).astype(np.uint8)
    elif cls == "double":
        out = u.astype(np.float64)
    else:
        out = (u.astype(np.float32) + np.float32(0.37)).astype(np.float32)
    return to_matlab(out)


def zero_border(msize, seed: int) -> np.ndarray:
    """DeconGPUTest/gpuEdgeErosionUsesTheCPUMask: zero triangles along x-z, as skewed data has;
    MATLAB's form"""
    a = from_matlab(blurred_beads(msize, seed))
    nz, _, nx = a.shape
    z = np.arange(1, nz + 1, dtype=np.float64).reshape(-1, 1, 1)
    x = np.arange(1, nx + 1, dtype=np.float64).reshape(1, 1, -1)
    keep = (x > 0.5 * z + 2) & (x < 0.5 * z + 30)
    return to_matlab(np.where(keep, a, np.uint16(0)))


def below_background(msize, seed: int) -> np.ndarray:
    """blurred beads with a slab at 60, below the background of 100; MATLAB's form"""
    a = blurred_beads(msize, seed, "uint16")
    a[: a.shape[0] // 3, :, : a.shape[2] // 2] = 60
    return a


def make_psf(msize, which: str = "f") -> np.ndarray:
    """the PSF (or a narrower one as the back projector; 'neg': one with negative side lobes, as
    OMW back projectors have) cropped to the data as RLdecon crops, a plane for 2D data; MATLAB's
    form"""
    sigma = (2.0, 1.5) if which == "f" else (1.2, 1.0)
    psf = decon_psf((31, 31, 31), sigma)
    if which == "neg":
        wide = decon_psf((31, 31, 31), (3.0, 2.5)) * np.float32(0.25)
        psf = (psf * np.float32(2.0) - wide).astype(np.float32)
    py = _py_size(msize)
    if len(py) == 2:
        psf = psf[15]
    s = [max(0, (p - n) // 2) for p, n in zip(psf.shape, py, strict=True)]
    t = [a + min(p, n) for a, p, n in zip(s, psf.shape, py, strict=True)]
    psf = psf[tuple(slice(a, b) for a, b in zip(s, t, strict=True))]
    return to_matlab(np.ascontiguousarray(psf))


def formula_vector(n: int, a: float) -> np.ndarray:
    """single(mod((1 : n)' * a, 1) - 0.5) (R04's sum probes), a MATLAB column"""
    k = np.arange(1, n + 1, dtype=np.float64)
    return (np.mod(k * a, 1.0) - 0.5).astype(np.float32).reshape(-1, 1)


# ---------------------------------------------------------------------------------------------
# the Python side and the cases
# ---------------------------------------------------------------------------------------------


def _iterations_out(r):
    return float(r.iterations)  # (0 without iterations, as MATLAB since M07)


def py_gpu(im, *args, method, nargout=1, **kw):
    if method == "omw":
        psf, psfb, numit = args
        r = _rl_gpu().decon_lucy_omw(im, psf, psfb, numit, use_gpu=True, **kw)
        out = (r.decon, r.err_mat)
    else:
        psf, numit = args
        r = _rl_gpu().decon_lucy(im, psf, numit, use_gpu=True, **kw)
        out = (r.decon, r.err_mat, _iterations_out(r))
    return out[0] if nargout == 1 else out[:nargout]


NAMES = {
    "Background": "background",
    "save16bit": "save_16bit",
    "dampFactor": "damp_factor",
    "scaleFactor": "scale_factor",
    "deconOffset": "decon_offset",
    "EdgeErosion": "edge_erosion",
    "bbox": "bbox",
    "debug": "debug",
    "debug_folder": "debug_folder",
    "saveStep": "save_step",
    "useMex": "use_mex",
}


def _voxels(msize) -> int:
    return math.prod(int(v) for v in msize)


def gpu_case(
    cid: str,
    method: str,
    inputs: dict,
    numit,
    opts: dict,
    *,
    msize,
    device: str,
    nargout: int = 1,
    files: bool | None = None,
    expect_error: bool | str = False,
    marks=(),
) -> Case:
    """a decon call on the GPU with MATLAB's options (name: value)"""
    vols = [Vol(In("im")), Vol(In("psf"))] + ([Vol(In("psfb"))] if method == "omw" else [])
    args = [*vols, numit, "useGPU", True]
    py_kwargs = {"method": method, "nargout": nargout}
    for k, v in opts.items():
        args += [k, Box(v) if k == "bbox" else v]
        py_kwargs[NAMES[k]] = Box(v) if k == "bbox" else v
    return Case(
        cid,
        matlab="decon_lucy_function" if method == "simplified" else "decon_lucy_omw_function",
        args=args,
        nargout=nargout,
        python=py_gpu,
        py_args=[*vols, numit],
        py_kwargs=py_kwargs,
        inputs=inputs,
        out_roles=("vol", None, None)[:nargout],
        files=files,
        expect_error=expect_error,
        py_error=(rl.DeconError, ValueError, TypeError, MemoryError) if expect_error else Exception,
        device=device,
        marks=marks,
    )


def _inputs(msize, seed, cls="uint16", data=None):
    return {
        "im": V(
            data if data is not None else (lambda m=msize, s=seed, c=cls: blurred_beads(m, s, c))
        ),
        "psf": V(lambda m=msize: make_psf(m, "f")),
        "psfb": V(lambda m=msize: make_psf(m, "b")),
    }


METHODS = {"rl": ("simplified", 10), "omw": ("omw", 2)}
PATHS = {"loop": {"useMex": False}, "dflt": {}, "mex": {"useMex": True}}
RFFT_OPTIONS = {
    "damp": {"dampFactor": 1.2},
    "erosion": {"EdgeErosion": 3},
    "bbox": {"bbox": [3, 4, 5, 40, 41, 30]},
    "all": {
        "EdgeErosion": 2,
        "scaleFactor": 1.5,
        "deconOffset": 2,
        "bbox": [2, 2, 2, 44, 40, 33],
        "save16bit": True,
    },
    # options single can't hold exactly (single with double: in double, rounded)
    "inexact": {"dampFactor": 1.37, "scaleFactor": 1.37, "deconOffset": 0.3},
}
DGT_OPTIONS = {
    "u16": {"save16bit": True},
    "single": {"save16bit": False},
    "dampErosionOffsetScale": {
        "save16bit": False,
        "dampFactor": 1.2,
        "EdgeErosion": 4,
        "deconOffset": 3,
        "scaleFactor": 1.5,
    },
}
GPU_SIZES = {"smooth": [48, 48, 40], "oddPrimes": [61, 53, 37]}


def model_cases(model: str) -> list[Case]:
    dev = f"gpu:{MODELS[model]}"
    out: list[Case] = []

    def add(cid, method, inputs, numit, opts, msize, **kw):
        out.append(
            gpu_case(f"{model}.{cid}", method, inputs, numit, opts, msize=msize, device=dev, **kw)
        )

    ms = [48, 46, 36]
    # DeconRFFTGPUTest/gpuMexMatchesGpuMatlabLoop: each path, uint16 and single data
    for m, (method, n) in METHODS.items():
        for cls in ("uint16", "single"):
            for p, po in PATHS.items():
                add(
                    f"{m}.{p}.{cls}.{size_name(ms)}",
                    method,
                    _inputs(ms, 1, cls),
                    n,
                    {"Background": 100, "save16bit": False, **po},
                    ms,
                )
    # gpuMexMatchesForOddAnd2DSizes (and prime factors, a larger volume)
    for seed, sz in enumerate(
        ([47, 45, 33], [64, 50], [1, 64, 33], [61, 53, 37], [181, 101, 17], [30, 28, 26]), 2
    ):
        for m, (method, n) in METHODS.items():
            for p in ("loop", "mex"):
                add(
                    f"{m}.{p}.size.{size_name(sz)}",
                    method,
                    _inputs(sz, seed),
                    n,
                    {"Background": 100, "save16bit": False, **PATHS[p]},
                    sz,
                )
    # gpuMexAppliesTheOptions
    for on, o in RFFT_OPTIONS.items():
        for m, (method, n) in METHODS.items():
            for p in ("loop", "mex"):
                add(
                    f"{m}.{p}.opt.{on}",
                    method,
                    _inputs(ms, 8),
                    n,
                    {"Background": 100, "save16bit": False, **o, **PATHS[p]},
                    ms,
                )
    # DeconGPUTest/gpuLoopMatchesThePreviousFormulation: classes, sizes and option sets (bbox: the
    # test's [3, 4, 2, gpuSize - [4, 3, 2]])
    for gn, gs in GPU_SIZES.items():
        opts = dict(
            DGT_OPTIONS,
            bbox={"save16bit": True, "bbox": [3, 4, 2, gs[0] - 4, gs[1] - 3, gs[2] - 2]},
        )
        for cls in ("uint16", "single", "double"):
            for on, o in opts.items():
                for m, (method, _) in METHODS.items():
                    n = 4 if method == "simplified" else 2
                    add(
                        f"{m}.dgt.{gn}.{cls}.{on}",
                        method,
                        _inputs(gs, 9, cls),
                        n,
                        {"Background": 100, **o},
                        gs,
                    )
    # backgrounds (uint16: MATLAB's uint16(Background) on the GPU), iterations
    for m, (method, _) in METHODS.items():
        for bg, cls in (
            (100.5, "uint16"),
            (99.4, "uint16"),
            (0.49999999999999994, "uint16"),
            (NAN, "uint16"),
            (-3.5, "uint16"),
            (70000.0, "uint16"),
            (100.3, "single"),
            (-5.0, "single"),
        ):
            bgn = "nan" if bg != bg else f"{bg:g}".replace("-", "m").replace(".", "p")
            for it in (0, 1):
                add(
                    f"{m}.loop.bg{bgn}.{cls}.iter{it}",
                    method,
                    _inputs(ms, 10, cls),
                    it,
                    {"Background": bg, "save16bit": False},
                    ms,
                )
        for it in (2, 3, 30) if method == "simplified" else (3, 10):
            for p in ("loop", "mex"):
                add(
                    f"{m}.{p}.iter{it}",
                    method,
                    _inputs(ms, 11),
                    it,
                    {"Background": 100.3, "save16bit": False, **PATHS[p]},
                    ms,
                    nargout=3 if method == "simplified" else 2,
                )
    # DeconGPUTest/gpuEdgeErosionUsesTheCPUMask: data with zero borders
    zb = [48, 48, 40]
    for m, (method, n) in METHODS.items():
        n = 5 if method == "simplified" else 2
        for e in (0, 4):
            inp = {**_inputs(zb, 12), "im": V(lambda: zero_border(zb, 12))}
            add(
                f"{m}.loop.zeroborder.e{e}",
                method,
                inp,
                n,
                {"Background": 100, "save16bit": False, **({"EdgeErosion": e} if e else {})},
                zb,
            )
    # a back projector with negative lobes on data below the background: J is +0 there and its
    # products with negative back projections -0, which max(J, 0) makes +0
    for p in ("loop", "mex"):
        inp = {
            **_inputs(ms, 14),
            "im": V(lambda: below_background(ms, 14)),
            "psfb": V(lambda: make_psf(ms, "neg")),
        }
        add(
            f"omw.{p}.negbp", "omw", inp, 3, {"Background": 100, "save16bit": False, **PATHS[p]}, ms
        )
    # debug output of MATLAB's GPU loop (err_mat, the saved estimates)
    for m, (method, _) in METHODS.items():
        for step, it in ((1, 5), (2, 6)):
            add(
                f"{m}.loop.debug.step{step}",
                method,
                _inputs([40, 44, 36], 13),
                it,
                {
                    "Background": 100,
                    "save16bit": False,
                    "debug": True,
                    "debug_folder": Work("dbg"),
                    "saveStep": step,
                },
                [40, 44, 36],
                nargout=3 if method == "simplified" else 2,
                files=True,
            )
    # MATLAB's gpuArray sum (single and double) of column vectors
    for n in (1, 7, 64, 128, 129, 256, 1000, 4095, 4096, 8192, 8193, 79488, 1000003, 10000000):
        for which in ("single", "double"):
            code = (
                "out = {gather(sum(gpuArray(in{1})))};"
                if which == "single"
                else "out = {gather(sum(gpuArray(in{1}), 'double'))};"
            )
            out.append(
                Case(
                    f"{model}.sum.{which}.n{n}",
                    code=code,
                    args=[In("x")],
                    python=py_sum,
                    py_kwargs={"double": which == "double"},
                    inputs={"x": V(lambda n=n: formula_vector(n, 0.6180339887498949))},
                    device=dev,
                )
            )
    return out


def py_sum(x, *, double):
    import cupy as cp

    return _rl_gpu().matlab_gpu_sum(cp.asarray(np.ascontiguousarray(x)), double=double).get()[0]


CASES: list[Case] = []
for _model in MODELS:
    CASES += model_cases(_model)

# errors of MATLAB's GPU code paths (the same on any GPU)
_ems = [20, 22, 12]
_ERR = {
    "bbox.outside": (
        {"bbox": [1, 1, 1, 21, 22, 12]},
        "uint16",
        "parallel:array:SubsrefSubscriptExceedsDim",
    ),
    "bbox.zero": (
        {"bbox": [0, 1, 1, 20, 22, 12]},
        "uint16",
        "parallel:array:SubsrefInvalidSubscript",
    ),
    "uint8": ({}, "uint8", "MATLAB:mixedClasses"),
    "mex.double": ({"useMex": True}, "double", True),
    "mex.debug": ({"useMex": True, "debug": True, "debug_folder": Work("dbg")}, "uint16", True),
}
for _en, (_o, _cls, _id) in _ERR.items():
    for _m, (_method, _) in METHODS.items():
        CASES.append(
            gpu_case(
                f"any.error.{_m}.{_en}",
                _method,
                _inputs(_ems, 20, _cls),
                2,
                {"Background": 100, **_o},
                msize=_ems,
                device="gpu",
                expect_error=_id,
            )
        )
for _m, (_method, _) in METHODS.items():
    for _p in ("loop", "mex"):
        CASES.append(
            gpu_case(
                f"any.error.{_m}.{_p}.psflarger",
                _method,
                {
                    **_inputs([5, 5, 4], 21),
                    "psf": V(lambda: make_psf([31, 31, 31], "f")),
                    "psfb": V(lambda: make_psf([31, 31, 31], "b")),
                },
                2,
                {"Background": 100, **PATHS[_p]},
                msize=[5, 5, 4],
                device="gpu",
                expect_error="decon_lucy_rfft_gpu_mex:outSizeIsSmallerThanPsfSize"
                if _p == "mex"
                else True,
            )
        )
# a reversed box gives an empty result on the GPU; a 2D result takes a box whose z range is 1 : 1
for _m, (_method, _) in METHODS.items():
    CASES.append(
        gpu_case(
            f"any.{_m}.bbox.reversed",
            _method,
            _inputs(_ems, 22),
            2,
            {"Background": 100, "bbox": [5, 1, 1, 3, 22, 12]},
            msize=_ems,
            device="gpu",
        )
    )
    CASES.append(
        gpu_case(
            f"any.{_m}.bbox.2d",
            _method,
            _inputs([30, 25], 23),
            2,
            {"Background": 100, "bbox": [2, 3, 1, 20, 21, 1]},
            msize=[30, 25],
            device="gpu",
        )
    )
    CASES.append(
        gpu_case(
            f"any.error.{_m}.bbox.2d",
            _method,
            _inputs([30, 25], 23),
            2,
            {"Background": 100, "bbox": [2, 3, 1, 20, 21, 2]},
            msize=[30, 25],
            device="gpu",
            expect_error="parallel:array:SubsrefSubscriptExceedsDim",
        )
    )


@functools.cache
def gpu_name() -> str:
    import cupy as cp

    name = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)["name"]
    return name.decode() if isinstance(name, bytes) else str(name)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.device.startswith("gpu:"):
        model = case.device[4:]
        if model.lower() not in gpu_name().lower():
            pytest.skip(f"the reference is MATLAB's on a {model}; this GPU is a {gpu_name()}")
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# without a GPU
# ---------------------------------------------------------------------------------------------


def _vol(shape, dtype=np.uint16):
    return np.ones(shape, dtype=dtype)


@pytest.mark.parametrize(
    "use_mex,dtype,debug,bg,expected",
    [
        (None, np.uint16, False, 100.0, (False, True)),  # MATLAB's GPU default is the loop
        (None, np.float32, False, 100.0, (False, True)),
        (None, np.float64, False, 100.0, (False, False)),
        (False, np.uint16, False, 100.0, (False, True)),
        (True, np.uint16, False, 100.0, (True, True)),
        (True, np.float32, False, None, (True, True)),
        (None, np.uint16, True, 100.0, (False, False)),  # no debug in the GPU mex
        (None, np.uint16, False, np.ones((4, 5, 6)), (False, False)),
    ],
)
def test_the_gpu_path_is_chosen_as_matlab_chooses_it(use_mex, dtype, debug, bg, expected):
    # DeconRFFTGPUTest/gpuMexMatchesGpuMatlabLoop (the default is the MATLAB loop) and
    # withoutTheGpuMexTheGpuLoopRuns (decon_use_rfft_mex with a GPU)
    g = _rl_gpu()
    assert g.use_gpu_kernels(use_mex, "simplified", _vol((4, 5, 6), dtype), debug, bg) == expected
    assert g.use_gpu_kernels(use_mex, "omw", _vol((4, 5, 6), dtype), debug, bg) == expected


def test_a_single_row_runs_in_the_gpu_mex():
    # the GPU mex takes data the CPU mex doesn't (one row along y)
    assert _rl_gpu().use_gpu_kernels(True, "omw", _vol((33, 1, 64)), False, 100.0) == (True, True)
    assert rl.use_rfft_kernels(None, "omw", _vol((33, 1, 64)), False, 100.0)[1] is False


@pytest.mark.parametrize(
    "dtype,debug,cls",
    [(np.float64, False, "double"), (np.uint8, False, "uint8"), (np.uint16, True, "uint16")],
)
def test_the_gpu_mex_refuses_what_matlab_refuses(dtype, debug, cls):
    with pytest.raises(rl.DeconError) as e:
        _rl_gpu().use_gpu_kernels(True, "simplified", _vol((12, 20, 22), dtype), debug, 100.0)
    assert e.value.identifier == ""
    assert str(e.value) == (
        f"decon_lucy_rfft_gpu_mex does not support the data (class {cls}, size [20 22 12]) or the "
        "options (debug on the GPU or a non-scalar background)."
    )
    with pytest.raises(ValueError):
        _rl_gpu().use_gpu_kernels(None, "other", _vol((2, 2)), False, 0)


def _g(shape):
    return _rl_gpu()._geom(shape)


def test_boxes_are_indexed_as_gpuarrays_index():
    g = _g((12, 20, 22))
    assert _rl_gpu()._box(None, g) == (0, 0, 0, 12, 20, 22)
    assert _rl_gpu()._box((1, 2, 3, 10, 18, 20), g) == (1, 2, 3, 9, 16, 17)
    # a reversed range is empty (MATLAB's 5 : 3), also with a start before the first element
    assert _rl_gpu()._box((0, 4, 0, 12, 3, 22), g) == (0, 4, 0, 12, 0, 22)
    assert _rl_gpu()._box((0, -1, 0, 12, -1, 22), g)[4] == 0
    with pytest.raises(rl.DeconError) as e:
        _rl_gpu()._box((0, 0, 0, 12, 21, 22), g)
    assert e.value.identifier == "parallel:array:SubsrefSubscriptExceedsDim"
    with pytest.raises(rl.DeconError) as e:
        _rl_gpu()._box((0, -1, 0, 12, 20, 22), g)
    assert e.value.identifier == "parallel:array:SubsrefInvalidSubscript"
    # an open end is the end of the array
    assert _rl_gpu()._box((2, 0, 0, None, None, None), g) == (2, 0, 0, 10, 20, 22)
    # a 2D result: a 4-value box, or a 6-value one whose z range is 0 .. 1
    g2 = _g((30, 25))
    assert _rl_gpu()._box((2, 3, 20, 21), g2) == (0, 2, 3, 1, 18, 18)
    assert _rl_gpu()._box((0, 2, 3, 1, 20, 21), g2) == (0, 2, 3, 1, 18, 18)
    with pytest.raises(rl.DeconError) as e:
        _rl_gpu()._box((0, 2, 3, 2, 20, 21), g2)
    assert e.value.identifier == "parallel:array:SubsrefSubscriptExceedsDim"
    # MATLAB checks y, then x, then z
    with pytest.raises(rl.DeconError) as e:
        _rl_gpu()._box((-1, 0, 0, 13, 21, 22), g)
    assert e.value.identifier == "parallel:array:SubsrefSubscriptExceedsDim"


def test_the_geometry_is_matlabs_layout():
    assert _g((5, 7, 9)) == ((5, 7, 9), (5, 7, 9), (5, 9, 7))
    assert _g((7, 9)) == ((7, 9), (1, 7, 9), (9, 7))
    assert _g((5, 7, 9)).n == 315


def test_memory_need_counts_the_buffers_of_each_path():
    sz = (500, 1800, 512)  # the demo tile
    n = math.prod(sz)
    # MATLAB's GPU loops take 38 (OMW) and 42 (RL) bytes per voxel for uint16 data (the sizing
    # report); the port's loops the same buffers, with the plan's work area one spectrum
    omw = _rl_gpu().memory_need("omw", sz, np.uint16)
    rlv = _rl_gpu().memory_need("simplified", sz, np.uint16)
    assert omw.total == 38 * n and omw.otfs == 16 * n and omw.plan == 8 * n
    assert rlv.total == 42 * n and rlv.otfs == 8 * n
    assert omw.per_voxel == 38.0
    # single data take 2 bytes per voxel more; the mask one; a known work area replaces the estimate
    assert _rl_gpu().memory_need("omw", sz, np.float32).total == 40 * n
    assert _rl_gpu().memory_need("omw", sz, np.uint16, mask=True).total == 39 * n
    assert _rl_gpu().memory_need("omw", sz, np.uint16, plan_bytes=0).total == 30 * n
    # the half spectra: the GPU mex's 16 bytes per voxel of an OMW call (K09's measurement) plus
    # the data and the half OTFs
    half = _rl_gpu().memory_need("omw", sz, np.uint16, use_mex=True)
    assert 25.9 < half.per_voxel < 26.1
    with pytest.raises(ValueError):
        _rl_gpu().memory_need("other", sz)


@pytest.mark.parametrize("n", [1, 128, 129, 8192, 8193, 79488, 461373440])
def test_the_sum_plan_is_matlabs(n):
    # MATLAB's gpuArray sum (R04's probes): one chunk up to 8192 elements, else chunks of 1408 k
    # elements, at most SMs * ceil(1024 / SMs) of them (A100 1080, B200 1036, Titan V 1040)
    for sms, gmax in ((80, 1040), (108, 1080), (148, 1036)):
        g, c = _rl_gpu().sum_plan(n, sms)
        if n <= 8192:
            assert (g, c) == (1, n)
            continue
        assert c % 1408 == 0 and g == -(-n // c) and g <= gmax
        assert c == 1408 or -(-n // (c - 1408)) > gmax


@pytest.mark.parametrize(
    "n,k",
    [
        (1048577, 1),
        (2097152, 2),
        (3000000, 2),
        (4194304, 3),
        (6000000, 4),
        (8388608, 6),
        (12000000, 8),
        (16777216, 12),
        (33554432, 23),
    ],
)
def test_the_chunk_lengths_measured_on_an_a100(n, k):
    # the chunk length MATLAB uses on an A100 (108 multiprocessors), found per n by R04's probes
    assert _rl_gpu().sum_plan(n, 108) == (-(-n // (1408 * k)), 1408 * k)


def test_importing_the_module_imports_no_cupy():
    code = (
        "import sys; import petakit5d.decon.rl_gpu as m; "
        "assert 'cupy' not in sys.modules, 'cupy imported'; print(m.__name__)"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=os.environ)
    assert r.returncode == 0, r.stderr


def test_cuda_source_is_device_code_without_atomics():
    src = _rl_gpu()._SOURCE
    assert "atomic" not in src and "#include" not in src
    for name in _rl_gpu()._NAMES:
        assert name.split("<")[0] in src, name
    # the products and sums that MATLAB rounds one by one are written as rounding intrinsics
    assert "__fmul_rn" in src and "__fadd_rn" in src and "__fdiv_rn" in src and "__dmul_rn" in src


def _no_gpu(monkeypatch, reason="no GPU here"):
    from petakit5d.gpu import device

    def fail(*a, **k):
        raise device.GPUUnavailableError(reason)

    monkeypatch.setattr(device, "select_device", fail)
    monkeypatch.setattr(device, "gpu_count", lambda refresh=False: device.GPUCount(0, 0, reason))


def _small():
    im = from_matlab(blurred_beads([20, 22, 12], 30))
    psf = from_matlab(make_psf([20, 22, 12], "f"))
    return im, psf


def test_without_a_usable_gpu_auto_runs_on_the_cpu_with_a_warning(monkeypatch):
    # DeconFunctionTest/defaultFallsBackToTheCPUWithAWarning for the GPU module called directly
    from petakit5d.gpu.device import GPUFallbackWarning

    _no_gpu(monkeypatch)
    im, psf = _small()
    with pytest.warns(GPUFallbackWarning):
        r = _rl_gpu().decon_lucy(im, psf, 3, background=100.0, use_gpu=None)
    assert np.array_equal(r.decon, rl.lucy(im, psf, 3, background=100.0, use_gpu=False).decon)
    with pytest.warns(GPUFallbackWarning):
        r = _rl_gpu().decon_lucy_omw(im, psf, psf, 2, background=100.0, use_gpu=None)
    assert np.array_equal(
        r.decon, rl.lucy_omw(im, psf, psf, 2, background=100.0, use_gpu=False).decon
    )


def test_without_a_usable_gpu_required_raises(monkeypatch):
    from petakit5d.gpu.device import NoUsableGPUError

    _no_gpu(monkeypatch)
    im, psf = _small()
    for call in (
        lambda: _rl_gpu().decon_lucy(im, psf, 3, use_gpu=True),
        lambda: _rl_gpu().decon_lucy_omw(im, psf, psf, 2, use_gpu=True),
    ):
        with pytest.raises(NoUsableGPUError) as e:
            call()
        assert e.value.identifier == "decon:noUsableGPU"


def test_use_gpu_false_runs_on_the_cpu(monkeypatch):
    _no_gpu(monkeypatch, "must not be asked")
    im, psf = _small()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        r = _rl_gpu().decon_lucy(im, psf, 3, background=100.0, use_gpu=False)
    assert np.array_equal(r.decon, rl.lucy(im, psf, 3, background=100.0).decon)


def test_the_cpu_module_hands_a_gpu_call_to_this_module(monkeypatch):
    # petakit5d.decon.rl: a selected GPU gives the call to rl_gpu with the registered keywords
    calls = []
    monkeypatch.setattr(rl, "_gpu_device", lambda use_gpu: object())

    def fake(name):
        def f(**kw):
            calls.append((name, kw))
            return name

        return f

    monkeypatch.setattr(_rl_gpu(), "decon_lucy", fake("decon_lucy"))
    monkeypatch.setattr(_rl_gpu(), "decon_lucy_omw", fake("decon_lucy_omw"))
    im, psf = _small()
    rl.clear_otf_cache()
    assert rl.decon_lucy(im, psf, 3, background=100.0, use_gpu=True) == "decon_lucy"
    assert rl.decon_lucy_omw(im, psf, psf, 2, use_gpu=None, edge_erosion=2) == "decon_lucy_omw"
    # data only the GPU mex takes (one row along y) with use_mex=True reach the GPU's rules
    row = from_matlab(blurred_beads([1, 64, 33], 31))
    prow = from_matlab(make_psf([1, 64, 33], "f"))
    assert rl.decon_lucy_omw(row, prow, prow, 2, use_gpu=True, use_mex=True) == "decon_lucy_omw"
    assert calls[-1][1]["use_mex"] is True
    names = set(calls[0][1])
    assert {
        "i",
        "psf",
        "numit",
        "background",
        "use_gpu",
        "save_16bit",
        "bbox",
        "use_mex",
        "overlap_fcn",
    } <= names
    assert calls[0][1]["background"] == 100.0 and calls[1][1]["edge_erosion"] == 2


def test_the_gpu_functions_take_the_registered_keywords():
    import inspect

    for cpu, gpu in (
        (rl.decon_lucy, _rl_gpu().decon_lucy),
        (rl.decon_lucy_omw, _rl_gpu().decon_lucy_omw),
    ):
        c = inspect.signature(cpu).parameters
        g = inspect.signature(gpu).parameters
        assert list(c) == list(g)[: len(c)]
        assert all(
            c[k].default == g[k].default for k in c if c[k].default is not inspect.Parameter.empty
        )
        assert set(g) - set(c) == {"layout", "keep_plans"}


# ---------------------------------------------------------------------------------------------
# on a GPU
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def gpu():
    cp = pytest.importorskip("cupy")
    from petakit5d.gpu import device, fft

    n, reason, _ = device.usable_gpu_count()
    if n == 0:
        pytest.skip(f"no usable GPU: {reason}")
    device.select_device()
    _rl_gpu().clear_otf_cache()
    yield cp
    _rl_gpu().clear_otf_cache()
    fft.clear_plan_cache()
    cp.get_default_memory_pool().free_all_blocks()


def _data(msize=(48, 46, 36), seed=3, cls="uint16"):
    return (
        from_matlab(blurred_beads(list(msize), seed, cls)),
        from_matlab(make_psf(list(msize), "f")),
        from_matlab(make_psf(list(msize), "b")),
    )


def _noise_data(msize, seed):
    """random uint16 data on a background of 100 (cheap to make at any size) and the PSFs"""
    py = _py_size(msize)
    im = CounterRNG(seed, "noise").integers(100, 1100, py, np.uint16)
    return im, from_matlab(make_psf(list(msize), "f")), from_matlab(make_psf(list(msize), "b"))


def _rel(a, b) -> float:
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def _run(method, im, psf, psfb, n, **kw):
    if method == "omw":
        return _rl_gpu().decon_lucy_omw(im, psf, psfb, n, **kw)
    return _rl_gpu().decon_lucy(im, psf, n, **kw)


def _cpu(method, im, psf, psfb, n, **kw):
    if method == "omw":
        return rl.lucy_omw(im, psf, psfb, n, use_gpu=False, **kw)
    return rl.lucy(im, psf, n, use_gpu=False, **kw)


@pytest.mark.gpu
@pytest.mark.parametrize("method", ["simplified", "omw"])
@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_gpu_matches_cpu(gpu, method, cls):
    # DeconGPUTest/gpuMatchesCPU: within MATLAB's own CPU-to-GPU difference, on both CPU FFT
    # backends (relative L2 1e-4 as the MATLAB test; measured: about 1e-6)
    from petakit5d import fft as P

    im, psf, psfb = _data(cls=cls)
    n = 5 if method == "simplified" else 2
    g = _run(method, im, psf, psfb, n, background=100.0, save_16bit=False, use_gpu=True).decon
    assert g.dtype == np.float32 and g.shape == im.shape
    for be in P.available():
        c = _cpu(method, im, psf, psfb, n, background=100.0, save_16bit=False, backend=be).decon
        assert _rel(g, c) < 1e-4, be
    u = _run(method, im, psf, psfb, n, background=100.0, use_gpu=True).decon
    cu = _cpu(method, im, psf, psfb, n, background=100.0).decon
    d = np.abs(u.astype(np.int64) - cu.astype(np.int64))
    assert d.max() <= 1 and (d > 0).mean() < 1e-3


@pytest.mark.gpu
@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_gpu_edge_erosion_uses_the_eroded_mask(gpu, method):
    # DeconGPUTest/gpuEdgeErosionUsesTheCPUMask: the result with erosion is the result without it
    # times the mask eroded on the CPU
    from petakit5d.decon.erosion import decon_mask_edge_erosion

    im = from_matlab(zero_border([48, 48, 40], 12))
    _, psf, psfb = _data([48, 48, 40], 12)
    mask = decon_mask_edge_erosion(im > 0, 4)
    assert mask.any() and not np.array_equal(mask, im > 0)
    n = 5 if method == "simplified" else 2
    kw = {"background": 100.0, "save_16bit": False, "use_gpu": True}
    plain = _run(method, im, psf, psfb, n, **kw).decon
    eroded = _run(method, im, psf, psfb, n, edge_erosion=4, **kw).decon
    assert np.array_equal(eroded, plain * mask.astype(np.float32))


@pytest.mark.gpu
@pytest.mark.parametrize("method", ["simplified", "omw"])
@pytest.mark.parametrize("use_mex", [None, True])
def test_the_cached_otfs_give_the_same_bits(gpu, method, use_mex):
    # DeconGPUTest/gpuLoopMatchesThePreviousFormulation's second call: the cached OTFs
    im, psf, psfb = _data((48, 48, 40), 9)
    kw = {"background": 100.0, "use_gpu": True, "use_mex": use_mex, "save_16bit": False}
    a = _run(method, im, psf, psfb, 4, **kw).decon
    b = _run(method, im, psf, psfb, 4, **kw).decon
    _rl_gpu().clear_otf_cache()
    c = _run(method, im, psf, psfb, 4, **kw).decon
    assert a.tobytes() == b.tobytes() == c.tobytes()
    # a PSF with other bits (-0 for +0) is a new OTF
    p2 = psf.copy()
    p2[0, 0, 0] = np.float32(-0.0) if p2[0, 0, 0] == 0 else p2[0, 0, 0]
    assert _run(method, im, p2, psfb, 4, **kw).decon.shape == a.shape


@pytest.mark.gpu
@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_gpu_overlap_function_runs_while_the_gpu_computes(gpu, method):
    # DeconOverlapTest/gpuOverlapFunctionRunsWhileTheGPUComputes: the same result, the function's
    # output returned; the function starts before the GPU has finished the iterations
    im, psf, psfb = _noise_data((256, 256, 96), 4)
    n = 10 if method == "simplified" else 4
    kw = {"background": 100.0, "use_gpu": True}
    expected = _run(method, im, psf, psfb, n, **kw)
    state = {}

    def overlap():
        # the queued iterations are still running on the stream
        state["busy"] = not gpu.cuda.get_current_stream().done
        return "read"

    r = _run(method, im, psf, psfb, n, overlap_fcn=overlap, **kw)
    assert np.array_equal(r.decon, expected.decon) and r.overlap_out == "read"
    assert state["busy"], "the overlap function ran after the GPU had finished"
    assert _run(method, im, psf, psfb, 1, **kw).overlap_out is None


@pytest.mark.gpu
@pytest.mark.parametrize("method", ["simplified", "omw"])
@pytest.mark.parametrize("use_mex", [False, True])
def test_peak_memory_is_at_most_matlabs(gpu, method, use_mex):
    # the GPU memory of a call: at most memory_need, which is MATLAB's 38 (OMW) and 42 (RL) bytes
    # per voxel of uint16 data for the loop at FFT-friendly sizes
    from petakit5d.gpu import fft, memory

    im, psf, psfb = _noise_data((128, 160, 64), 5)
    n = 3 if method == "simplified" else 2
    _rl_gpu().clear_otf_cache()
    fft.clear_plan_cache()
    memory.release_cached()
    with memory.track_peak() as p:
        _run(method, im, psf, psfb, n, background=100.0, use_gpu=True, use_mex=use_mex)
    need = _rl_gpu().memory_need(method, im.shape, im.dtype, use_mex=use_mex)
    vox = im.size
    assert p.peak_above_start <= need.total + (1 << 20), (p.peak_above_start / vox, need.per_voxel)
    if not use_mex:
        assert p.peak_above_start / vox <= (38 if method == "omw" else 42) + 0.5


@pytest.mark.gpu
def test_out_of_memory_is_matlabs_error_and_leaves_no_arrays(gpu, monkeypatch):
    # the identifiers XR_RLdeconFrame3D catches to deconvolve on the CPU instead, and the GPU's
    # memory free again for that
    from petakit5d.gpu import memory
    from petakit5d.kernels.decon_gpu import GPUOutOfMemoryError

    im, psf, psfb = _data((48, 46, 36), 6)
    pool = gpu.get_default_memory_pool()
    _rl_gpu().clear_otf_cache()
    memory.release_cached(plans=True)
    before = pool.used_bytes()
    big = _rl_gpu().MemoryNeed(1 << 60, 0, 0, 0.0)
    monkeypatch.setattr(_rl_gpu(), "memory_need", lambda *a, **k: big)
    for method in ("simplified", "omw"):
        with pytest.raises(GPUOutOfMemoryError) as e:
            _run(method, im, psf, psfb, 3, background=100.0, use_gpu=True)
        assert e.value.identifier == "parallel:gpu:array:OOM"
    _rl_gpu().clear_otf_cache()
    memory.release_cached(plans=True)
    assert pool.used_bytes() <= before
    # an allocation that fails inside the iterations: the same error, nothing left behind
    monkeypatch.setattr(
        _rl_gpu(), "memory_need", lambda *a, **k: _rl_gpu().MemoryNeed(0, 0, 0, 0.0)
    )
    real = _rl_gpu()._single

    def fail(a):
        raise gpu.cuda.memory.OutOfMemoryError(1 << 40, 0)

    monkeypatch.setattr(_rl_gpu(), "_single", fail)
    with pytest.raises(GPUOutOfMemoryError) as e:
        _run("simplified", im, psf, psfb, 3, background=100.0, use_gpu=True)
    assert e.value.identifier == "parallel:gpu:array:OOM"
    monkeypatch.setattr(_rl_gpu(), "_single", real)
    _rl_gpu().clear_otf_cache()
    memory.release_cached(plans=True)
    assert pool.used_bytes() <= before


@pytest.mark.gpu
@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_the_native_half_spectra_are_within_k09s_tolerance(gpu, method):
    # layout='native' (Python only): the mex's steps on half spectra of the (z, y, x) arrays,
    # within the accepted native-layout tolerance of the mex's own result
    im, psf, psfb = _data((48, 46, 36), 7)
    n = 10 if method == "simplified" else 2
    kw = {"background": 100.0, "use_gpu": True, "use_mex": True, "save_16bit": False}
    a = _run(method, im, psf, psfb, n, **kw).decon
    b = _run(method, im, psf, psfb, n, layout="native", **kw).decon
    d = np.abs(a.astype(np.float64) - b)
    assert _rel(b, a) <= 1e-5 and d.max() <= 1e-5 * np.abs(a).max()


@pytest.mark.gpu
def test_large_transfers_go_through_the_pinned_buffers_intact(gpu):
    # data and results larger than a staging buffer, with a short last piece
    n = 3 * _rl_gpu()._CHUNK // 2 + 4099
    a = CounterRNG(5, "transfer").integers(0, 65535, (n,), np.uint16)
    d = _rl_gpu()._to_gpu(a)
    assert d.dtype == np.uint16 and d.shape == a.shape
    assert np.array_equal(d.get(), a)
    assert np.array_equal(_rl_gpu()._to_host(d), a)
    s = np.arange(7, dtype=np.float32)
    assert np.array_equal(_rl_gpu()._to_host(_rl_gpu()._to_gpu(s)), s)


@pytest.mark.gpu
def test_keep_plans_false_frees_the_plans(gpu):
    from petakit5d.gpu import fft

    im, psf, psfb = _data((48, 46, 36), 8)
    fft.clear_plan_cache()
    _rl_gpu().decon_lucy_omw(im, psf, psfb, 2, background=100.0, use_gpu=True, keep_plans=False)
    assert fft.cached_plan_bytes() == 0
    _rl_gpu().decon_lucy_omw(im, psf, psfb, 2, background=100.0, use_gpu=True)
    assert fft.cached_plan_bytes() > 0


def np_block(x: np.ndarray) -> np.generic:
    """one block of 256 threads as MATLAB's gpuArray sum (CUB's block reduction), in NumPy"""
    n = len(x)
    it = -(-n // 256)
    pad = np.zeros(it * 256, x.dtype)
    pad[:n] = x
    r = pad.reshape(it, 256)
    acc = r[0].copy()
    for row in r[1:]:
        acc = acc + row
    w = acc.reshape(8, 32).copy()
    while w.shape[1] > 1:
        w = w[:, 0::2] + w[:, 1::2]
    t = w[0, 0]
    for k in range(1, 8):
        if 32 * k < min(n, 256):
            t = x.dtype.type(t + w[k, 0])
    return t


def np_matlab_sum(x: np.ndarray, sms: int) -> np.generic:
    """MATLAB's gpuArray sum in NumPy (float32 or float64 x: the accumulator's class)"""
    n = len(x)
    if n <= 128:
        v = np.zeros(128, x.dtype)
        v[:n] = x
        h = 64
        while h:
            v[:h] = v[:h] + v[h : 2 * h]
            h //= 2
        return v[0]
    if n <= 8192:
        return np_block(x)
    g, c = _rl_gpu().sum_plan(n, sms)
    parts = np.array([np_block(x[b * c : min(n, (b + 1) * c)]) for b in range(g)], x.dtype)
    return np_matlab_sum(parts, sms)


@pytest.mark.gpu
@pytest.mark.parametrize("n", [100, 200, 8192, 8193, 50000, 1048583, 3000003])
def test_the_sum_follows_matlabs_order(gpu, n):
    # the gpuArray sum's order (R04's probes; MATLAB's results are the *.sum.* oracle cases): a
    # halving tree up to 128 elements, one CUB block up to 8192, chunks of 1408 k beyond
    sms = gpu.cuda.Device().attributes["MultiProcessorCount"]
    x = formula_vector(n, 0.7548776662466927).reshape(-1)
    s = _rl_gpu().matlab_gpu_sum(gpu.asarray(x)).get()[0]
    assert s.view(np.uint32) == np_matlab_sum(x, sms).view(np.uint32)
    d = _rl_gpu().matlab_gpu_sum(gpu.asarray(x), double=True).get()[0]
    assert d.view(np.uint64) == np_matlab_sum(x.astype(np.float64), sms).view(np.uint64)


@pytest.mark.gpu
def test_the_sum_of_nothing_is_zero(gpu):
    assert _rl_gpu().matlab_gpu_sum(gpu.asarray(np.zeros(0, np.float32))).get()[0] == 0
    with pytest.raises(ValueError):
        _rl_gpu().matlab_gpu_sum(gpu.asarray(np.zeros(3)))


@pytest.mark.gpu
def test_results_do_not_change_from_run_to_run(gpu):
    im, psf, psfb = _data((61, 53, 37), 9)
    a = _rl_gpu().decon_lucy(im, psf, 10, background=100.0, use_gpu=True, save_16bit=False).decon
    b = _rl_gpu().decon_lucy(im, psf, 10, background=100.0, use_gpu=True, save_16bit=False).decon
    assert a.tobytes() == b.tobytes()


@pytest.mark.gpu
def test_a_volume_of_one_slice_is_a_2d_image(gpu):
    im, psf, psfb = _data((64, 50), 10)
    for method in ("simplified", "omw"):
        a = _run(method, im, psf, psfb, 5, background=100.0, use_gpu=True, save_16bit=False).decon
        b = _run(
            method,
            im[None],
            psf[None],
            psfb[None],
            5,
            background=100.0,
            use_gpu=True,
            save_16bit=False,
        ).decon
        assert b.shape == (1, 64, 50) and a.tobytes() == b[0].tobytes()


@pytest.mark.gpu
def test_the_registered_functions_run_on_the_gpu(gpu):
    # petakit5d.decon.rl's registered functions with use_gpu None hand the call here
    im, psf, psfb = _data((48, 46, 36), 11)
    a = rl.decon_lucy_omw(im, psf, psfb, 2, background=100.0, save_16bit=False).decon
    b = (
        _rl_gpu()
        .decon_lucy_omw(im, psf, psfb, 2, background=100.0, save_16bit=False, use_gpu=True)
        .decon
    )
    assert a.tobytes() == b.tobytes()
