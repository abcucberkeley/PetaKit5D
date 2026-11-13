"""The deconvolution kernels on the GPU (petakit5d.kernels.decon_gpu) against the CUDA mex files.

Ports the intent of tests/decon/DeconRFFTGPUTest.m at the kernel level (gpuHalfOTFIsTheFirstHalfOf
TheOTF, gpuBackgroundIsSubtractedAsInTheGpuLoop, outOfMemoryLeavesNoGpuArraysBehind, the kernel
parts of gpuMexMatchesGpuMatlabLoop and gpuMexMatchesForOddAnd2DSizes) and the gpuArray methods of
DeconMaskErosionTest.m (gpuArrayMasksEqualTheCPUMasks, fullGpuArrayMasksThinnerThanTheErosion).

Oracle cases (GPU packs: MATLAB R2024b runs the oracle tree's decon_lucy_rfft_gpu_mex and
decon_mask_edge_erosion_gpu_mex; its cuFFT is 11.0.8, the port's pin):

- ``<gpu>.*`` (gpu: A100, B200, TitanV): the mex's 'otf', 'omw' and 'rlstep' on the same GPU model
  (cuFFT's kernels, and so its bits, depend on the GPU). Each case runs three Python sides:
  ``.zyx`` and ``.int`` in MATLAB's layout ((z, y, x) volumes, and RL's internal arrays), bit for
  bit; ``.native`` with the native half spectra, within the proposed tolerance (NATIVE_TOL; the
  counts are in the K09 report). Cases of another GPU model are skipped.
- ``any.*``: results that don't depend on the GPU (the erosion is logical; the mex's argument
  errors), on any GPU.

Case ids: ``<gpu>.otf.<size>``, ``<gpu>.omw2.uint16.<size>``, ``<gpu>.omw2.single.<size>``,
``A100.omw0.u16.bg100p5.48x46x36`` (iterations, data, background, size), ``<gpu>.rlstep.<size>``,
``<gpu>.rl5.uint16.<size>``; ``any.erosion.<mask>.<size>``, ``any.erosion.full.*``,
``any.erosion.mex.*``, ``any.erosion.error.*``, ``any.rfft.error.*``.

The tests without the ``gpu`` mark check what needs no GPU: shapes and plans, the mex's argument
checks and background rounding, the neighborhood offsets, the CUDA sources.
"""

from __future__ import annotations

import dataclasses
import functools
import gc
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.kernels import decon_gpu as G
from petakit5d.kernels.decon import DeconKernelError
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, PerAxis, Value, Vol, mcell
from petakit5d.testing.compare import tol
from petakit5d.testing.synth import CounterRNG, beads, decon_psf

SRC_ROOT = str(Path(G.__file__).resolve().parents[2])

# the GPU models of the oracle cases: case id prefix -> part of MATLAB's (and CUDA's) GPU name
MODELS = {"A100": "A100", "B200": "B200", "TitanV": "TITAN V"}

NATIVE_TOL = tol(
    reason="proposed (K09): the native (z, y, x) half spectra on the GPU, cuFFT's other plans",
    rel_max=1e-5,
    rel_l2=1e-5,
)
NATIVE_OTF = tol(
    reason="proposed (K09): native half OTFs on the GPU, cuFFT's other plans",
    rel_max=1e-6,
    rel_l2=1e-6,
)


def V(make):
    return Value(make)


# ---------------------------------------------------------------------------------------------
# inputs, in MATLAB's form (as K05's: DeconRFFTTest's blurred beads and PSFs)
# ---------------------------------------------------------------------------------------------

BLUR_SIGMA = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)


def blurred_beads(msize, seed: int, cls: str = "uint16") -> np.ndarray:
    """DeconRFFTGPUTest.blurredBeads: beads blurred by the test PSF plus a little noise, on a
    background of 100; MATLAB's form"""
    sz = tuple(int(v) for v in msize)
    py = (sz[2], sz[0], sz[1]) if len(sz) == 3 else sz
    rng = CounterRNG(seed, "beads")
    n = max(3, round(math.prod(sz) / 4000))
    pos = rng.uniform((n, len(py))) * (np.array(py, dtype=np.float64) - 1)
    im = beads(py, pos, BLUR_SIGMA if len(py) == 3 else BLUR_SIGMA[1], 1000.0, 100.0)
    im = im + np.float32(3) * rng.uniform(py, np.float32)
    if cls == "uint16":
        out = np.floor(im + np.float32(0.5)).astype(np.uint16)  # im >= 0
    else:
        out = (im + np.float32(0.37)).astype(np.float32)
    return to_matlab(out)


def make_psf(msize, which: str = "f") -> np.ndarray:
    """the PSF (or a back projector) cropped to the data as RLdecon crops, scaled by a power of
    two to about unit sum; 'neg' has negative side lobes, as OMW back projectors have; MATLAB's
    form"""
    sigma = (2.0, 1.5) if which == "f" else (1.2, 1.0)
    psf = decon_psf((31, 31, 31), sigma)
    psf = np.ascontiguousarray(psf * np.float32(2.0**-19))
    if which == "neg":
        wide = decon_psf((31, 31, 31), (3.0, 2.5)) * np.float32(2.0**-21)
        psf = np.ascontiguousarray(psf * np.float32(2.0) - wide).astype(np.float32)
    if len(msize) == 2:
        psf = psf[15]
        py = tuple(msize)
    else:
        py = (msize[2], msize[0], msize[1])
    s = [max(0, (p - n) // 2) for p, n in zip(psf.shape, py, strict=True)]
    t = [a + min(p, n) for a, p, n in zip(s, psf.shape, py, strict=True)]
    psf = psf[tuple(slice(a, b) for a, b in zip(s, t, strict=True))]
    return to_matlab(np.ascontiguousarray(psf))


def below_background(msize, seed: int) -> np.ndarray:
    """blurred beads with a slab below the background (J is +0 there, and products with it -0
    where the back projection is negative: the update makes them +0)"""
    a = blurred_beads(msize, seed, "uint16")
    if a.ndim == 3:
        a[: a.shape[0] // 3, :, : a.shape[2] // 2] = 60
    else:
        a[: a.shape[0] // 3, :] = 60
    return a


def size_name(s) -> str:
    return "x".join(str(v) for v in s)


# ---------------------------------------------------------------------------------------------
# the Python sides (MATLAB's arguments mapped by the roles; results on the host)
# ---------------------------------------------------------------------------------------------


def native_half_to_matlab(h: np.ndarray, shape) -> np.ndarray:
    """MATLAB's half OTF H (halved along y: (ny // 2 + 1, nx, nz), or (ny // 2 + 1, nx)) of a
    native half spectrum (halved along x), through the full spectrum's Hermitian symmetry"""
    vol = h if len(shape) == 3 else h.reshape((1, *h.shape))
    nz, ny, nx = (1, *shape) if len(shape) == 2 else shape
    nxh = vol.shape[-1]
    full = np.empty((nz, ny, nx), dtype=vol.dtype)
    full[..., :nxh] = vol
    kz, ky, kx = (-np.arange(nz)) % nz, (-np.arange(ny)) % ny, nx - np.arange(nxh, nx)
    full[..., nxh:] = np.conj(vol[kz][:, ky][:, :, kx])
    m = to_matlab(full)[: ny // 2 + 1]  # (ky, kx, kz)
    return m if len(shape) == 3 else m[..., 0]


def py_otf(psf, out_size, *, layout="matlab", internal=False):
    shape = tuple(int(v) for v in out_size)
    h = G.LucyRFFTGPU(shape, layout=layout).otf(psf).get()
    return h if layout == "matlab" else native_half_to_matlab(h, shape)


def py_omw(psf, raw, psfb, n_iter, bg, *, layout="matlab", internal=False):
    k = G.LucyRFFTGPU(raw.shape, layout=layout)
    hf, hb = k.otf(psf), k.otf(psfb)
    if internal:
        return k.from_internal(k.omw(k.to_internal(raw), hf, hb, n_iter, bg, internal=True)).get()
    return k.omw(raw, hf, hb, n_iter, bg).get()


def py_rlstep(psf, y, i, *, layout="matlab", internal=False):
    k = G.LucyRFFTGPU(y.shape, layout=layout)
    h = k.otf(psf)
    if internal:
        return k.from_internal(k.rlstep(k.to_internal(y), k.to_internal(i), h, internal=True)).get()
    return k.rlstep(y, i, h).get()


def py_rl_loop(psf, i, n, *, layout="matlab", internal=False):
    """n RL steps in a row, J <- rlstep(J, I, H) from J = single(I)"""
    import cupy as cp

    k = G.LucyRFFTGPU(i.shape, layout=layout)
    h = k.otf(psf)
    if internal:
        ii = k.to_internal(i)
        j = k.to_internal(i.astype(np.float32))
        for _ in range(int(n)):
            j = k.rlstep(j, ii, h, internal=True)
        return k.from_internal(j).get()
    ii = cp.asarray(i)
    j = ii.astype(cp.float32)
    for _ in range(int(n)):
        j = k.rlstep(j, ii, h)
    return j.get()


# MATLAB's code of the composed cases (in{} as listed in args)
OMW_CODE = (
    "sz = size(in{2}); H = decon_lucy_rfft_gpu_mex('otf', gpuArray(in{1}), sz); "
    "Hb = decon_lucy_rfft_gpu_mex('otf', gpuArray(in{3}), sz); "
    "out = {gather(decon_lucy_rfft_gpu_mex('omw', gpuArray(in{2}), H, Hb, in{4}, in{5}))};"
)
RLSTEP_CODE = (
    "H = decon_lucy_rfft_gpu_mex('otf', gpuArray(in{1}), size(in{2})); "
    "out = {gather(decon_lucy_rfft_gpu_mex('rlstep', gpuArray(in{2}), gpuArray(in{3}), H))};"
)
RL_LOOP_CODE = (
    "H = decon_lucy_rfft_gpu_mex('otf', gpuArray(in{1}), size(in{2})); I = gpuArray(in{2}); "
    "J = single(I); for k = 1 : in{3}, J = decon_lucy_rfft_gpu_mex('rlstep', J, I, H); end; "
    "out = {gather(J)};"
)


def three_sides(case: Case, spec_native=NATIVE_TOL) -> list[Case]:
    """the case's Python sides: MATLAB's layout on (z, y, x) volumes and on internal arrays (bit
    for bit), and the native layout (within the tolerance); one MATLAB reference"""
    out = []
    for name, kw, spec in (
        ("zyx", {"layout": "matlab"}, case.compare),
        ("int", {"layout": "matlab", "internal": True}, case.compare),
        ("native", {"layout": "native"}, spec_native),
    ):
        k = dict(case.py_kwargs or {})
        k.update(kw)
        out.append(dataclasses.replace(case, id=f"{case.id}.{name}", py_kwargs=k, compare=spec))
    return out


# MATLAB sizes [ny nx nz]: DeconRFFTGPUTest's (odd and 2D, a single row), prime and smooth sizes
RFFT_SIZES = [
    [48, 46, 36],
    [47, 45, 33],
    [64, 50],
    [1, 64, 33],
    [61, 53, 37],
    [128, 96, 40],
    [33, 1, 20],
    [181, 101, 17],
    [30, 28, 26],  # the PSF cropped to even sizes along every axis (its shift by floor(size / 2))
]


def model_cases(model: str) -> list[Case]:
    dev = f"gpu:{MODELS[model]}"
    cases: list[Case] = []

    def add(c: Case, spec_native=NATIVE_TOL):
        cases.extend(three_sides(dataclasses.replace(c, device=dev), spec_native))

    for seed, ms in enumerate(RFFT_SIZES, 1):
        sn = size_name(ms)
        inp = {
            "psf": V(lambda m=ms: make_psf(m, "f")),
            "psfb": V(lambda m=ms: make_psf(m, "b")),
            "u16": V(lambda m=ms, s=seed: blurred_beads(m, s, "uint16")),
            "f32": V(lambda m=ms, s=seed: blurred_beads(m, s, "single")),
        }
        add(
            Case(
                f"{model}.otf.{sn}",
                matlab="decon_lucy_rfft_gpu_mex",
                args=["otf", Vol(In("psf")), [float(v) for v in ms]],
                python=py_otf,
                py_args=[Vol(In("psf")), PerAxis([float(v) for v in ms], "yxz")],
                inputs={"psf": inp["psf"]},
            ),
            NATIVE_OTF,
        )
        add(
            Case(
                f"{model}.omw2.uint16.{sn}",
                code=OMW_CODE,
                args=[Vol(In("psf")), Vol(In("u16")), Vol(In("psfb")), 2, 100],
                python=py_omw,
                inputs={k: inp[k] for k in ("psf", "psfb", "u16")},
                out_roles=("vol",),
            )
        )
        add(
            Case(
                f"{model}.rlstep.{sn}",
                code=RLSTEP_CODE,
                args=[Vol(In("psf")), Vol(In("y")), Vol(In("i"))],
                python=py_rlstep,
                inputs={
                    "psf": inp["psf"],
                    "y": V(
                        lambda m=ms, s=seed: np.maximum(
                            blurred_beads(m, s, "single") - np.float32(99.5), np.float32(0)
                        )
                    ),
                    "i": V(lambda m=ms, s=seed: blurred_beads(m, s, "single") - np.float32(100)),
                },
                out_roles=("vol",),
            )
        )
        if seed <= 4:
            add(
                Case(
                    f"{model}.omw2.single.{sn}",
                    code=OMW_CODE,
                    args=[Vol(In("psf")), Vol(In("f32")), Vol(In("psfb")), 2, 100.3],
                    python=py_omw,
                    inputs={k: inp[k] for k in ("psf", "psfb", "f32")},
                    out_roles=("vol",),
                )
            )
            add(
                Case(
                    f"{model}.rlstep.uint16.{sn}",
                    code=RLSTEP_CODE,
                    args=[Vol(In("psf")), Vol(In("y")), Vol(In("u"))],
                    python=py_rlstep,
                    inputs={
                        "psf": inp["psf"],
                        "y": V(lambda m=ms, s=seed: blurred_beads(m, s, "single")),
                        "u": V(lambda m=ms, s=seed: blurred_beads(m, s + 100, "uint16")),
                    },
                    out_roles=("vol",),
                )
            )
    # iterations 0, 1, 5; backgrounds the mex rounds or saturates; data below the background
    ms = [48, 46, 36]
    sn = size_name(ms)
    inp = {
        "psf": V(lambda: make_psf(ms, "f")),
        "psfb": V(lambda: make_psf(ms, "b")),
        "neg": V(lambda: make_psf(ms, "neg")),
        "u16": V(lambda: blurred_beads(ms, 77, "uint16")),
        "f32": V(lambda: blurred_beads(ms, 77, "single")),
        "low": V(lambda: below_background(ms, 90)),
    }
    for n_iter, bg, data, bp in [
        (0, 100.5, "u16", "psfb"),
        (0, 99.4, "u16", "psfb"),
        (0, 100.3, "f32", "psfb"),
        (0, NAN, "u16", "psfb"),
        (0, -3.5, "u16", "psfb"),
        (0, 70000.0, "u16", "psfb"),
        (0, 0.49999999999999994, "u16", "psfb"),
        (1, 100, "u16", "psfb"),
        (5, 100, "u16", "psfb"),
        (5, 100.3, "f32", "psfb"),
        (2, -5.0, "f32", "psfb"),
        # a back projector with negative side lobes on data below the background: -0 products
        (3, 100, "low", "neg"),
    ]:
        bgn = "nan" if bg != bg else f"{bg:g}".replace("-", "m").replace(".", "p")
        add(
            Case(
                f"{model}.omw{n_iter}.{data}.bg{bgn}.{sn}",
                code=OMW_CODE,
                args=[Vol(In("psf")), Vol(In(data)), Vol(In(bp)), n_iter, bg],
                python=py_omw,
                inputs={k: inp[k] for k in ("psf", bp, data)},
                out_roles=("vol",),
            )
        )
    for ms in ([48, 46, 36], [64, 50]):
        add(
            Case(
                f"{model}.rl5.uint16.{size_name(ms)}",
                code=RL_LOOP_CODE,
                args=[Vol(In("psf")), Vol(In("i")), 5],
                python=py_rl_loop,
                inputs={
                    "psf": V(lambda m=ms: make_psf(m, "f")),
                    "i": V(lambda m=ms: np.maximum(blurred_beads(m, 78), 100) - np.uint16(100)),
                },
                out_roles=("vol",),
            )
        )
    return cases


CASES: list[Case] = []
for _model in MODELS:
    CASES += model_cases(_model)


# ---------------------------------------------------------------------------------------------
# the erosion on gpuArray masks (decon_mask_edge_erosion.m and decon_mask_edge_erosion_gpu_mex)
# ---------------------------------------------------------------------------------------------

MASK_SIZES = {
    "odd3D": [37, 41, 23],
    "small": [5, 5, 5],
    "twoSlices": [17, 19, 2],
    "oneRow": [1, 30, 20],
    "oneColumn": [30, 1, 20],
    "oneSlice": [30, 20, 1],
    "image2D": [45, 47],
    "row2D": [1, 30],
}


def make_mask(shape: str, msize, seed: int = 0) -> np.ndarray:
    """DeconMaskErosionTest.makeMask in MATLAB's form: random, skewed data's zero triangles along
    x-z, a thin plane, one true voxel, blobs"""
    sz3 = list(msize) + [1] * (3 - len(msize))
    rng = CounterRNG(seed, f"mask/{shape}/{size_name(msize)}")
    py = tuple(int(v) for v in ((sz3[2], sz3[0], sz3[1]) if len(msize) == 3 else msize))
    if shape == "random":
        a = rng.uniform(py) > 0.3
    elif shape == "skewed":
        z = np.arange(1, sz3[2] + 1, dtype=np.float64).reshape(-1, 1, 1)
        x = np.arange(1, sz3[1] + 1, dtype=np.float64).reshape(1, 1, -1)
        a = (x > 0.6 * z + sz3[1] * 0.05) & (x < 0.6 * z + sz3[1] * 0.5)
        a = np.broadcast_to(a, (sz3[2], sz3[0], sz3[1])).reshape(py)
    elif shape == "single":
        mo = np.zeros(math.prod(sz3), bool)
        mo[math.ceil(mo.size / 2) - 1] = True
        return np.reshape(mo, msize, order="F")
    else:  # blobs: imdilate(rand > 0.99, ones(3, 3))
        from petakit5d.matlabcompat.morphology import imdilate

        a = imdilate(rng.uniform(py) > 0.99, np.ones((3, 3), bool))
    return to_matlab(np.ascontiguousarray(a))


def gpu_erosion(mask, edge_erosion):
    return G.decon_mask_edge_erosion_gpu(mask, edge_erosion).get()


def gpu_erosion_all(mask, erosions):
    return [gpu_erosion(mask, float(e)) for e in erosions]


ERODE_CODE = (
    "A = gpuArray(in{1}); out = {arrayfun(@(E) gather(decon_mask_edge_erosion(A, E)), in{2}, "
    "'unif', 0)};"
)
ERODE_ONE = "out = {gather(decon_mask_edge_erosion(gpuArray(in{1}), in{2}))};"
EROSIONS = [float(e) for e in range(1, 21)] + [36.0, 41.0]

for _shape in ["random", "skewed", "single", "blobs"]:
    for _sname, _ms in MASK_SIZES.items():
        if make_mask(_shape, _ms).all():
            continue
        CASES.append(
            Case(
                f"any.erosion.{_shape}.{_sname}",
                code=ERODE_CODE,
                args=[Vol(In("a")), EROSIONS],
                python=gpu_erosion_all,
                inputs={"a": V(lambda s=_shape, m=_ms: make_mask(s, m))},
                out_roles=(("each", "vol"),),
                device="gpu",
            )
        )


def skewed_with_holes():
    """DeconMaskErosionTest/gpuArrayMasksEqualTheCPUMasks: a skewed mask with sparse holes"""
    a = make_mask("skewed", [61, 57, 33])
    return a & (to_matlab(CounterRNG(11, "holes").uniform((33, 61, 57))) > 0.001)


CASES.append(
    Case(
        "any.erosion.skewedHoles",
        code=ERODE_CODE,
        args=[Vol(In("a")), [1.0, 2.0, 3.0, 4.0, 9.0, 20.0, 36.0]],
        python=gpu_erosion_all,
        inputs={"a": V(skewed_with_holes)},
        out_roles=(("each", "vol"),),
        device="gpu",
    )
)
CASES.append(
    Case(
        "any.erosion.numeric",
        code=ERODE_ONE,
        args=[Vol(In("a")), 7],
        python=gpu_erosion,
        inputs={
            "a": V(
                lambda: (
                    make_mask("skewed", [40, 50, 30]).astype(np.uint16)
                    * to_matlab(CounterRNG(3).integers(1, 1000, (30, 40, 50), np.uint16))
                ).astype(np.uint16)
            )
        },
        out_roles=("vol",),
        device="gpu",
    )
)
# full masks get their edges set directly, also when thinner than 2 * EdgeErosion, and all edge
# when thinner than EdgeErosion (M07, 834604d; before, MATLAB's GPU code errored there as the
# CPU's did); a fractional erosion errors, NaN takes index 1
# (fullGpuArrayMasksThinnerThanTheErosion)
for _ms, _es in (
    ([3, 20, 20], (1.0, 2.0, 3.0)),
    ([20, 3, 20], (1.0, 2.0, 3.0)),
    ([20, 20, 3], (1.0, 2.0, 3.0, 2.5, NAN)),
    ([1, 30, 20], (1.0, 2.0, 3.0)),
    ([3, 20], (1.0, 2.0, 3.0)),
    ([30, 1], (1.0, 2.0, 3.0)),
):
    for _e in _es:
        _bad = _e == _e and _e != math.floor(_e)
        _en = "nan" if _e != _e else f"{_e:g}".replace(".", "p")
        CASES.append(
            Case(
                f"any.erosion.full.{size_name(_ms)}.e{_en}",
                code=ERODE_ONE,
                args=[Vol(In("a")), _e],
                python=gpu_erosion,
                inputs={"a": V(lambda m=_ms: np.ones(m, dtype=bool))},
                out_roles=("vol",),
                device="gpu",
                expect_error=True if _bad else None,
                py_error=Exception if _bad else None,
            )
        )
# erosion values the mex doesn't take: MATLAB's imerode code (and strel's errors); a single 3 is
# taken as the double 3 since M07 (834604d), so the mex erodes with it (the case keeps its name)
for _e, _err in ((2.5, True), (np.float32(3.0), False), (0.0, True), (-1.0, True)):
    CASES.append(
        Case(
            f"any.erosion.nonmex.{np.asarray(_e).dtype.name}.{float(_e):g}".replace("-", "m"),
            code=ERODE_ONE,
            args=[Vol(In("a")), _e],
            python=gpu_erosion,
            inputs={"a": V(lambda: make_mask("random", [12, 14, 6]))},
            out_roles=("vol",),
            device="gpu",
            expect_error=True if _err else None,
            py_error=Exception if _err else None,
        )
    )


def random_nhoods(trial: int):
    """DeconMaskErosionTest/mexEqualsImerodeForOtherNeighborhoods: 1 to 3 random (asymmetric,
    even-sized) neighborhoods, a random mask, and morphop's padding for them"""
    from petakit5d.matlabcompat.morphology import Strel

    rng = CounterRNG(trial, "nhoods")
    ns = int(rng.integers(1, 3, 1)[0])
    nh = []
    for _ in range(ns):
        sh = tuple(int(v) for v in rng.integers(1, 6, 2))
        e = rng.uniform(sh) > 0.5
        e.reshape(-1)[int(rng.integers(0, e.size - 1, 1)[0])] = True
        nh.append(e)
    msz = [int(v) for v in rng.integers(1, 30, 2)] + [int(rng.integers(2, 4, 1)[0])]
    a = rng.uniform((msz[2], msz[0], msz[1])) > 0.2
    pu = pl = (0, 0)
    if ns > 1:
        se = Strel(np.ones((1, 1), bool), [Strel(x) for x in nh])
        if not se.is_decomposition_orthogonal():
            before, after = se.getpadsize()
            pu, pl = tuple(after[:2]), tuple(before[:2])
    return to_matlab(np.ascontiguousarray(a)), nh, [float(v) for v in pu], [float(v) for v in pl]


def gpu_erosion_mex(a, nh, pu, pl, rz, clear):
    nh = list(getattr(nh, "items", nh))  # (an empty MATLAB cell, mcell())
    return G.mask_edge_erosion(a, nh, pu, pl, rz, clear).get()


EROSION_MEX = "decon_mask_edge_erosion_gpu_mex"
for _t in range(16):
    _a, _nh, _pu, _pl = random_nhoods(_t)
    CASES.append(
        Case(
            f"any.erosion.mex.random.t{_t}",
            matlab=EROSION_MEX,
            args=[Vol(In("a")), _nh, _pu, _pl, float(_t % 3), _t % 2 == 1],
            python=gpu_erosion_mex,
            inputs={"a": V(lambda t=_t: random_nhoods(t)[0])},
            out_roles=("vol",),
            device="gpu",
        )
    )
for _name, _args, _make in [
    (
        "zline",
        [[np.array([[True]])], [0.0, 0.0], [0.0, 0.0], 4.0, False],
        lambda: to_matlab(CounterRNG(5).uniform((40, 20, 30)) > 0.1),
    ),
    (
        "zline_long",
        [[np.array([[True]])], [0.0, 0.0], [0.0, 0.0], 30.0, False],
        lambda: to_matlab(CounterRNG(9).uniform((12, 20, 30)) > 0.02),
    ),
    (
        "numeric_nhood",
        [[np.array([[0.0, 2.0, 1.0]])], [0.0, 0.0], [0.0, 0.0], 0.0, False],
        lambda: to_matlab(CounterRNG(6).uniform((5, 10, 12)) > 0.5),
    ),
    (
        "empty",
        [[np.array([[True]])], [0.0, 0.0], [0.0, 0.0], 1.0, True],
        lambda: np.zeros((0, 3, 2), bool),
    ),
    (
        "no_nhoods",
        [mcell(), [0.0, 0.0], [0.0, 0.0], 2.0, True],
        lambda: to_matlab(CounterRNG(10).uniform((7, 9, 11)) > 0.2),
    ),
    (
        "clear_border_pads",
        [[np.ones((3, 5), bool), np.ones((4, 2), bool)], [2.0, 1.0], [0.0, 3.0], 2.0, True],
        lambda: to_matlab(CounterRNG(7).uniform((9, 21, 18)) > 0.15),
    ),
    (
        "wide_slices",
        [[np.ones((1, 7), bool), np.ones((7, 1), bool)], [0.0, 0.0], [0.0, 0.0], 1.0, True],
        lambda: to_matlab(CounterRNG(12).uniform((3, 70, 300)) > 0.05),
    ),
]:
    CASES.append(
        Case(
            f"any.erosion.mex.{_name}",
            matlab=EROSION_MEX,
            args=[Vol(In("a")), *_args],
            python=gpu_erosion_mex,
            inputs={"a": V(_make)},
            out_roles=("vol",),
            device="gpu",
        )
    )


def erosion_error(which: str):
    a = CounterRNG(8).uniform((5, 10, 12)) > 0.5
    t = [np.array([[True]])]
    z = [0.0, 0.0]
    if which == "double":
        return G.mask_edge_erosion(a.astype(np.float64), t, z, z, 1, True)
    if which == "4d":
        return G.mask_edge_erosion(np.ones((2, 3, 3, 3), bool), t, z, z, 1, True)
    if which == "not_cell":
        return G.mask_edge_erosion(a, np.array([[True]]), z, z, 1, True)
    if which == "false_nhood":
        return G.mask_edge_erosion(a, [np.zeros((3, 3), bool)], z, z, 1, True)
    if which == "3d_nhood":
        return G.mask_edge_erosion(a, [np.ones((2, 2, 2), bool)], z, z, 1, True)
    if which == "neg_pad":
        return G.mask_edge_erosion(a, t, [0.0, -1.0], z, 1, True)
    if which == "int_pad":
        return G.mask_edge_erosion(a, t, np.array([0, 0], np.int32), z, 1, True)
    if which == "scalar_pad":
        return G.mask_edge_erosion(a, t, z, 0, 1, True)
    if which == "neg_rz":
        return G.mask_edge_erosion(a, t, z, z, -1, True)
    if which == "frac_rz":
        return G.mask_edge_erosion(a, t, z, z, 1.5, True)
    if which == "nan_nhood":
        return G.mask_edge_erosion(a, [np.array([[1.0, NAN]])], z, z, 0, True)
    if which == "too_many":
        return G.mask_edge_erosion(a, [np.ones((65, 65), bool)], z, z, 0, True)
    raise AssertionError(which)


_EROSION_ID = "decon_mask_edge_erosion_gpu_mex:inputError"
for _w, _code in [
    ("double", "decon_mask_edge_erosion_gpu_mex(double(A), {true}, [0, 0], [0, 0], 1, true)"),
    ("4d", "decon_mask_edge_erosion_gpu_mex(true(3, 3, 3, 2), {true}, [0, 0], [0, 0], 1, true)"),
    ("not_cell", "decon_mask_edge_erosion_gpu_mex(A, true, [0, 0], [0, 0], 1, true)"),
    ("false_nhood", "decon_mask_edge_erosion_gpu_mex(A, {false(3)}, [0, 0], [0, 0], 1, true)"),
    ("3d_nhood", "decon_mask_edge_erosion_gpu_mex(A, {true(2, 2, 2)}, [0, 0], [0, 0], 1, true)"),
    ("neg_pad", "decon_mask_edge_erosion_gpu_mex(A, {true}, [0, -1], [0, 0], 1, true)"),
    ("int_pad", "decon_mask_edge_erosion_gpu_mex(A, {true}, int32([0, 0]), [0, 0], 1, true)"),
    ("scalar_pad", "decon_mask_edge_erosion_gpu_mex(A, {true}, [0, 0], 0, 1, true)"),
    ("neg_rz", "decon_mask_edge_erosion_gpu_mex(A, {true}, [0, 0], [0, 0], -1, true)"),
    ("frac_rz", "decon_mask_edge_erosion_gpu_mex(A, {true}, [0, 0], [0, 0], 1.5, true)"),
    ("nan_nhood", "decon_mask_edge_erosion_gpu_mex(A, {[1, NaN]}, [0, 0], [0, 0], 0, true)"),
    ("too_many", "decon_mask_edge_erosion_gpu_mex(A, {true(65)}, [0, 0], [0, 0], 0, true)"),
]:
    CASES.append(
        Case(
            f"any.erosion.error.{_w}",
            code=f"A = gpuArray(in{{1}}); {_code};",
            args=[Vol(In("a"))],
            python=erosion_error,
            py_args=[],
            py_kwargs={"which": _w},
            inputs={"a": V(lambda: to_matlab(CounterRNG(8).uniform((5, 10, 12)) > 0.5))},
            # (the mex's call of logical() raises MATLAB's error, not the mex's)
            expect_error="MATLAB:nologicalnan" if _w == "nan_nhood" else _EROSION_ID,
            py_error=DeconKernelError,
            device="gpu",
        )
    )


# the RL mex's argument errors
def rfft_error(im, y, *, which: str):
    k = G.LucyRFFTGPU(im.shape)
    h = k.otf(np.full((3, 3, 3), np.float32(1 / 27)))
    if which == "short_otf":
        return k.omw(im, h[..., :-1], h, 2, 100)
    if which == "double":
        return k.omw(im.astype(np.float64), h, h, 2, 100)
    if which == "real_otf":
        return k.omw(im, h.T.real.copy().T, h, 2, 100)  # (no cupy.asfortranarray: cuBLAS)
    if which == "host_otf":
        return k.omw(im, h.get(order="F"), h, 2, 100)
    if which == "rl_double":
        return k.rlstep(y.astype(np.float64), im, h)
    if which == "rl_sizes":
        return k.rlstep(y, y[:-1], h)
    if which == "complex_data":
        return k.omw(im.astype(np.complex64), h, h, 2, 100)
    if which == "bg_vector":
        return k.omw(im, h, h, 2, np.array([100.0, 1.0]))
    if which == "niter_vector":
        return k.omw(im, h, h, np.array([2.0, 3.0]), 100)
    if which == "psf_larger":
        return G.otf(np.ones((5, 5, 5), np.float32), (10, 4, 12))
    if which == "psf_double":
        return G.otf(np.ones((3, 3, 3)), (10, 16, 12))
    raise AssertionError(which)


def _rfft_err_inputs():
    im = to_matlab(np.ones((10, 16, 12), np.uint16) * np.uint16(150))
    return {"im": V(lambda: im), "y": V(lambda: im.astype(np.float32))}


_RFFT_ID = "decon_lucy_rfft_gpu_mex:inputError"
_ERR_H = "H = decon_lucy_rfft_gpu_mex('otf', gpuArray(single(ones(3, 3, 3) / 27)), [16, 12, 10]); "
for _name, _code in [
    ("short_otf", "decon_lucy_rfft_gpu_mex('omw', in{1}, H(:, :, 1 : end - 1), H, 2, 100)"),
    ("double", "decon_lucy_rfft_gpu_mex('omw', double(in{1}), H, H, 2, 100)"),
    ("real_otf", "decon_lucy_rfft_gpu_mex('omw', in{1}, real(H), H, 2, 100)"),
    ("host_otf", "decon_lucy_rfft_gpu_mex('omw', in{1}, gather(H), H, 2, 100)"),
    ("rl_double", "decon_lucy_rfft_gpu_mex('rlstep', double(in{2}), in{1}, H)"),
    ("rl_sizes", "decon_lucy_rfft_gpu_mex('rlstep', in{2}, in{2}(:, :, 1 : end - 1), H)"),
    ("complex_data", "decon_lucy_rfft_gpu_mex('omw', complex(single(in{1})), H, H, 2, 100)"),
    ("bg_vector", "decon_lucy_rfft_gpu_mex('omw', in{1}, H, H, 2, [100, 1])"),
    ("niter_vector", "decon_lucy_rfft_gpu_mex('omw', in{1}, H, H, [2, 3], 100)"),
    ("psf_larger", "decon_lucy_rfft_gpu_mex('otf', single(ones(5, 5, 5)), [4, 12, 10])"),
    ("psf_double", "decon_lucy_rfft_gpu_mex('otf', ones(3, 3, 3), [16, 12, 10])"),
]:
    CASES.append(
        Case(
            f"any.rfft.error.{_name}",
            code=_ERR_H + _code + ";",
            args=[Vol(In("im")), Vol(In("y"))],
            python=rfft_error,
            py_kwargs={"which": _name},
            inputs=_rfft_err_inputs(),
            expect_error="decon_lucy_rfft_gpu_mex:outSizeIsSmallerThanPsfSize"
            if _name == "psf_larger"
            else _RFFT_ID,
            py_error=DeconKernelError,
            device="gpu",
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
# without a GPU: shapes, plans, the mex's checks, the sources
# ---------------------------------------------------------------------------------------------


def mex_fft_dims(msize):
    """the mex's fft_dims of MATLAB's size [ny nx (nz)]: (nz if > 1, nx if > 1, ny)"""
    d = list(msize) + [1] * (3 - len(msize))
    return tuple([v for v in (d[2], d[1]) if v > 1] + [d[0]])


@pytest.mark.parametrize("msize", [*RFFT_SIZES, [5, 1, 1], [1, 1, 7], [1, 9], [9, 1], [1, 1]])
def test_geometry_has_the_mexs_plans_and_matlabs_sizes(msize):
    shape = from_matlab(np.zeros(msize, dtype=np.uint8)).shape
    g = G.geometry(shape)
    ny, nx = msize[0], msize[1]
    nz = msize[2] if len(msize) == 3 else 1
    assert g.zyx == (nz, ny, nx)
    assert g.plan_dims == mex_fft_dims(msize)
    assert g.fft_shape == ((nz, nx, ny) if len(msize) == 3 else (nx, ny))
    assert g.half_shape == (*g.fft_shape[:-1], ny // 2 + 1)
    # MATLAB's H and internal arrays, as the mex sizes them
    assert g.spectrum_shape == ((ny // 2 + 1, nx, nz) if len(msize) == 3 else (ny // 2 + 1, nx))
    assert g.internal_shape == tuple(msize)
    assert math.prod(g.half_shape) == g.nh and g.n == ny * nx * nz
    assert g.psf_strides() == (ny * nx, 1, ny)
    n = G.geometry(shape, "native")
    assert n.fft_shape == shape and n.internal_shape == shape
    assert n.half_shape == n.spectrum_shape == (*shape[:-1], nx // 2 + 1)
    assert n.psf_strides() == (ny * nx, nx, 1)


def test_geometry_rejects_what_the_mex_rejects():
    for bad in [(5,), (2, 3, 4, 5), (0, 3, 4), (3, 2**31, 2)]:
        with pytest.raises(DeconKernelError) as e:
            G.geometry(bad)
        assert e.value.identifier == "decon_lucy_rfft_gpu_mex:inputError"
    with pytest.raises(ValueError):
        G.geometry((4, 5, 6), "fortran")


@pytest.mark.parametrize(
    ("bg", "expected"),
    [
        (100.5, 101),
        (99.4, 99),
        (99.5, 100),
        (0.0, 0),
        (-3.5, 0),
        (NAN, 0),
        (65534.5, 65535),
        (70000.0, 65535),
        (float("inf"), 65535),
        (float("-inf"), 0),
        # floor(b + 0.5) in double: 1, where MATLAB's uint16 gives 0 (the mex's rounding)
        (0.49999999999999994, 1),
        (2.5, 3),
    ],
)
def test_background_is_rounded_as_the_mex_rounds(bg, expected):
    assert G.background_uint16(bg) == expected


def _host_args(shape=(6, 8, 10)):
    raw = np.full(shape, 150, np.uint16)
    return raw, G.LucyRFFTGPU(shape)


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda k, r: k.omw(r, None, None, "2", 100), "nIter must be a real numeric scalar."),
        (lambda k, r: k.omw(r, None, None, [2, 3], 100), "nIter must be a real numeric scalar."),
        (lambda k, r: k.omw(r, None, None, 3e9, 100), "nIter is too large."),
        (
            lambda k, r: k.omw(r, None, None, 2, [1, 2]),
            "background must be a real numeric scalar.",
        ),
        (lambda k, r: k.omw(r, None, None, 2, 1j), "background must be a real numeric scalar."),
        (
            lambda k, r: k.omw(r.astype(np.float64), None, None, 2, 100),
            "The data must be a real, non-empty uint16 or single array.",
        ),
        (
            lambda k, r: k.omw(r[:0], None, None, 2, 100),
            "The data must be a real, non-empty uint16 or single array.",
        ),
        (
            lambda k, r: k.omw(r, np.zeros((4, 8, 6), np.complex64), None, 2, 100),
            "The OTFs must be complex single half OTF gpuArrays.",
        ),
        (
            lambda k, r: k.rlstep(r, r, None),
            "Y must be a real single and I a real uint16 or single array.",
        ),
        (
            lambda k, r: k.rlstep(r.astype(np.float32), r[:-1], None),
            "Y and I must have the same size.",
        ),
        (
            lambda k, r: k.otf(np.ones((3, 3), np.float64)),
            "The PSF must be a real, non-empty single array.",
        ),
        (lambda k, r: k.otf(np.ones((7, 3, 3), np.float32)), "outSizeIsSmallerThanPsfSize"),
        (
            lambda k, r: k.otf(np.ones((1, 1, 1, 1), np.float32)),
            "Only 2D and 3D data with dimensions below 2^31 are supported.",
        ),
    ],
)
def test_the_mexs_argument_errors_need_no_gpu(call, message):
    raw, k = _host_args()
    with pytest.raises(DeconKernelError) as e:
        call(k, raw)
    assert str(e.value) == message
    assert e.value.identifier.startswith("decon_lucy_rfft_gpu_mex:")


def test_one_over_n_is_the_single_of_the_double_quotient():
    # the mex's invN = (float)(1.0 / (double)n): at N = 101 x 181 x 1031 (not a single) it differs
    # from a division in single
    n = 101 * 181 * 1031
    k = G.LucyRFFTGPU((101, 181, 1031))
    assert k._inv_n == np.float32(1.0 / n)
    assert k._inv_n != np.float32(1) / np.float32(n)


def test_counts_are_matlab_loop_counts():
    assert G._count(2.7, "n") == 2
    assert G._count(0.5, "n") == 0 and G._count(-4, "n") == 0 and G._count(NAN, "n") == 0
    assert G._count(np.uint8(3), "n") == 3 and G._count(True, "n") == 1


def test_nhood_offsets_are_in_matlabs_column_major_order():
    nh = [np.array([[0, 1, 0], [1, 1, 1]], bool), np.array([[2.0], [0.0], [-1.0], [0.5]])]
    start, dr, dc = G.nhood_offsets(nh)
    # first: rows 2, columns 3, center (1, 2) 1-based; second: 4 x 1, center (2, 1)
    assert start == [0, 4, 7]
    assert list(zip(dr, dc, strict=True))[:4] == [(1, -1), (0, 0), (1, 0), (1, 1)]
    assert list(zip(dr, dc, strict=True))[4:] == [(-1, 0), (1, 0), (2, 0)]
    with pytest.raises(DeconKernelError) as e:
        G.nhood_offsets([np.array([[1.0, NAN]])])
    assert e.value.identifier == "MATLAB:nologicalnan"


@pytest.mark.parametrize(
    "which",
    [
        "double",
        "4d",
        "not_cell",
        "false_nhood",
        "3d_nhood",
        "neg_pad",
        "int_pad",
        "scalar_pad",
        "neg_rz",
        "frac_rz",
        "too_many",
    ],
)
def test_the_erosion_mexs_argument_errors_need_no_gpu(which):
    with pytest.raises(DeconKernelError) as e:
        erosion_error(which)
    assert e.value.identifier == _EROSION_ID


def test_pads_must_be_doubles_as_matlab_classes_them():
    assert G._pad_arg([2, 3], "p") == (2, 3)
    assert G._pad_arg((0.0, 1.0), "p") == (0, 1)
    assert G._pad_arg(np.array([1.0, 4.0]), "p") == (1, 4)
    for bad in (np.array([1, 2]), [True, False], [1.5, 2], [1, 2, 3]):
        with pytest.raises(DeconKernelError):
            G._pad_arg(bad, "p")


def test_cuda_sources_are_device_code_with_the_kernels_python_calls():
    from petakit5d.gpu import kernels

    # (the files the package compiles: next to the extension modules, or in the source tree)
    rl = kernels.cuda_source("decon", "rl_rfft.cu").read_text(encoding="utf-8")
    er = kernels.cuda_source("decon", "edge_erosion.cu").read_text(encoding="utf-8")
    for src in (rl, er):
        assert not re.search(r"#include\s*[<\"](cufft|cuda_runtime|mex|stdint|float|math)", src)
        assert not re.search(r"\batomic[A-Z]\w*\(", src)  # element-wise and per-column work only
    for name in ("multiply_kernel", "update_kernel", "scatter_psf_kernel", "update_out_t"):
        assert re.search(rf'extern "C" __global__ void {name}\(', rl), name
    assert re.search(r'extern "C" __global__ void update_y_t\(', rl)
    for expr in G._RL_NAMES:
        assert re.search(rf"__global__ void {expr.split('<')[0]}\(", rl), expr
    for name in ("load_kernel", "erode_kernel", "z_kernel"):
        assert re.search(rf'extern "C" __global__ void {name}\(', er), name
    # the mex's rounding intrinsics stay (they pin the bits whatever the compiler options)
    assert rl.count("__fmul_rn") >= 6 and "__fdiv_rn" in rl


def test_importing_the_module_imports_no_cupy():
    code = (
        "import sys\nimport petakit5d.kernels.decon_gpu as G\n"
        "G.geometry((5, 6, 7)); G.background_uint16(3.5)\n"
        "print('cupy' in sys.modules)"
    )
    env = dict(os.environ, PYTHONPATH=SRC_ROOT + os.pathsep + os.environ.get("PYTHONPATH", ""))
    r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "False"


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
    fft.clear_plan_cache()
    yield cp
    fft.clear_plan_cache()
    cp.get_default_memory_pool().free_all_blocks()


def bits(a) -> np.ndarray:
    a = np.ascontiguousarray(a.get() if hasattr(a, "get") else a)
    return a.view(np.uint32 if a.dtype == np.float32 else np.uint64)


def data_and_psfs(msize, seed=3):
    raw = from_matlab(blurred_beads(msize, seed, "uint16"))
    psf = from_matlab(make_psf(msize, "f"))
    psfb = from_matlab(make_psf(msize, "b"))
    return raw, psf, psfb


def matlab_gpu_loop(raw, psf, psfb, n_iter, bg, method):
    """MATLAB's own GPU loops of decon_lucy_omw_function and decon_lucy_function (useMex false) on
    full spectra (gpuArray fftn/ifftn: petakit5d.gpu.fft.MatlabFFTPlan), with CuPy's element-wise
    operations (which can flush subnormals and contract products: a reference within rounding)"""
    import cupy as cp

    from petakit5d.gpu import fft

    shape = raw.shape
    plan = fft.MatlabFFTPlan(shape)

    def full_otf(p):
        pad = cp.zeros(shape, cp.float32)
        pad[tuple(slice(0, n) for n in p.shape)] = cp.asarray(p)
        pad = cp.roll(pad, tuple(-(n // 2) for n in p.shape), axis=tuple(range(len(shape))))
        return plan.forward(pad)

    i = cp.asarray(raw).astype(cp.float32)
    i = cp.maximum(i - cp.float32(G.background_uint16(bg)), 0)
    eps = cp.float32(2.0**-52)
    hf = full_otf(psf)
    if method == "omw":
        hb = full_otf(psfb)
        j = i.copy()
        for _ in range(n_iter):
            cx = cp.maximum(plan.inverse(plan.forward(j) * hf, real=True), eps)
            j = cp.maximum(plan.inverse(plan.forward(i / cx) * hb, real=True) * j, 0)
        return j
    y = i.copy()
    for _ in range(n_iter):
        cx = cp.maximum(plan.inverse(plan.forward(y) * hf, real=True), eps)
        y = cp.maximum(y * plan.inverse(cp.conj(hf) * plan.forward(i / cx), real=True), 0)
    return y


def rel_l2(a, b) -> float:
    a = np.asarray(a.get() if hasattr(a, "get") else a, np.float64)
    b = np.asarray(b.get() if hasattr(b, "get") else b, np.float64)
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


@pytest.mark.gpu
@pytest.mark.parametrize("layout", G.LAYOUTS)
@pytest.mark.parametrize("msize", [[48, 46, 36], [47, 45, 33], [64, 50], [1, 64, 33]])
def test_kernels_match_matlabs_gpu_loops(gpu, layout, msize):
    # DeconRFFTGPUTest/gpuMexMatchesGpuMatlabLoop and gpuMexMatchesForOddAnd2DSizes at the kernel
    # level: the half-spectrum steps against MATLAB's full-spectrum GPU loops, relative L2 3e-5
    raw, psf, psfb = data_and_psfs(msize)
    k = G.LucyRFFTGPU(raw.shape, layout=layout)
    hf, hb = k.otf(psf), k.otf(psfb)
    j = k.omw(raw, hf, hb, 2, 100)
    assert rel_l2(j, matlab_gpu_loop(raw, psf, psfb, 2, 100, "omw")) < 3e-5
    i = gpu.maximum(gpu.asarray(raw).astype(gpu.float32) - 100, 0)
    y = i.copy()
    for _ in range(4):
        y = k.rlstep(y, i, hf)
    assert rel_l2(y, matlab_gpu_loop(raw, psf, psfb, 4, 100, "rl")) < 3e-5


@pytest.mark.gpu
@pytest.mark.parametrize("layout", G.LAYOUTS)
def test_kernels_match_the_cpu_kernels(gpu, layout):
    # against K05's CPU kernels (the default FFT backend): GPU and CPU differ by rounding only
    from petakit5d.kernels import decon as K

    raw, psf, psfb = data_and_psfs([61, 53, 37])
    k = G.LucyRFFTGPU(raw.shape, layout=layout)
    c = K.LucyRFFT(raw.shape)
    j = k.omw(raw, k.otf(psf), k.otf(psfb), 3, 100)
    jc = c.omw(raw, c.otf(psf), c.otf(psfb), 3, 100)
    assert rel_l2(j, jc) < 3e-5
    y = np.maximum(raw.astype(np.float32) - 99.5, np.float32(0))
    i = (raw.astype(np.float32) - np.float32(100)).astype(np.float32)
    assert rel_l2(k.rlstep(y, i, k.otf(psf)), c.rlstep(y, i, c.otf(psf))) < 3e-5


@pytest.mark.gpu
@pytest.mark.parametrize("msize", [[48, 46, 36], [47, 45, 33], [64, 50], [1, 64, 33], [33, 1, 20]])
def test_internal_arrays_and_volumes_give_the_same_bits(gpu, msize):
    raw, psf, psfb = data_and_psfs(msize)
    k = G.LucyRFFTGPU(raw.shape)
    hf, hb = k.otf(psf), k.otf(psfb)
    ri = k.to_internal(raw)
    assert ri.shape == k.internal_shape and ri.flags.f_contiguous and ri.dtype == gpu.uint16
    assert np.array_equal(ri.get(), to_matlab(raw))  # MATLAB's array
    assert np.array_equal(k.from_internal(ri).get(), raw)
    for n in (0, 1, 3):
        a = k.omw(raw, hf, hb, n, 100.5)
        b = k.from_internal(k.omw(ri, hf, hb, n, 100.5, internal=True))
        assert np.array_equal(bits(a), bits(b)), n
    y = gpu.asarray(np.maximum(raw.astype(np.float32) - 99.5, np.float32(0)))
    i16 = gpu.asarray(np.maximum(raw, 100) - np.uint16(100))
    a = k.rlstep(y, i16, hf)
    b = k.from_internal(k.rlstep(k.to_internal(y), k.to_internal(i16), hf, internal=True))
    assert np.array_equal(bits(a), bits(b))
    # uint16 I equals single I of the same values (the kernels read uint16 as it is)
    assert np.array_equal(bits(a), bits(k.rlstep(y, i16.astype(gpu.float32), hf)))


@pytest.mark.gpu
@pytest.mark.parametrize("layout", G.LAYOUTS)
def test_background_is_subtracted_as_in_the_gpu_loop(gpu, layout):
    # DeconRFFTGPUTest/gpuBackgroundIsSubtractedAsInTheGpuLoop: uint16 data, the saturating
    # subtraction of uint16(background); single data, max(I - background, 0) in double
    raw, psf, _ = data_and_psfs([48, 46, 36])
    k = G.LucyRFFTGPU(raw.shape, layout=layout)
    h = k.otf(psf)
    for bg in (100.5, 99.4):
        out = k.omw(raw, h, h, 0, bg).get()
        b = np.uint16(G.background_uint16(bg))
        assert np.array_equal(out, np.where(raw > b, raw - b, 0).astype(np.float32))
    ims = (raw.astype(np.float32) + np.float32(0.37)).astype(np.float32)
    out = k.omw(ims, h, h, 0, 100.3).get()
    exp = np.maximum((ims.astype(np.float64) - 100.3).astype(np.float32), np.float32(0))
    assert np.array_equal(bits(out), bits(exp))


@pytest.mark.gpu
@pytest.mark.parametrize("msize", [[48, 46, 36], [47, 45, 33], [64, 50], [1, 64, 33]])
def test_half_otf_is_the_first_half_of_the_otf(gpu, msize):
    # DeconRFFTGPUTest/gpuHalfOTFIsTheFirstHalfOfTheOTF: the rows up to the middle of MATLAB's
    # first dimension of decon_psf2otf's OTF (here NumPy's FFT in double), within 1e-6
    psf = from_matlab(make_psf(msize, "f"))
    shape = from_matlab(np.zeros(msize, np.uint8)).shape
    pad = np.zeros(shape)
    pad[tuple(slice(0, n) for n in psf.shape)] = psf
    full = np.fft.fftn(
        np.roll(pad, tuple(-(n // 2) for n in psf.shape), axis=(0, 1, 2)[: pad.ndim])
    )
    h = G.LucyRFFTGPU(shape).otf(psf)
    assert h.shape == G.geometry(shape).spectrum_shape and h.flags.f_contiguous
    assert np.abs(h.get() - to_matlab(full)[: msize[0] // 2 + 1]).max() < 1e-6
    hn = G.LucyRFFTGPU(shape, layout="native").otf(psf).get()
    assert np.abs(hn - full[..., : shape[-1] // 2 + 1]).max() < 1e-6


@pytest.mark.gpu
def test_results_do_not_change_from_run_to_run(gpu):
    raw, psf, psfb = data_and_psfs([61, 53, 37])
    for layout in G.LAYOUTS:
        k = G.LucyRFFTGPU(raw.shape, layout=layout, keep_plans=False)
        first = bits(k.omw(raw, k.otf(psf), k.otf(psfb), 3, 100))
        for _ in range(3):
            assert np.array_equal(first, bits(k.omw(raw, k.otf(psf), k.otf(psfb), 3, 100)))


@pytest.mark.gpu
def test_kernels_run_on_the_current_stream(gpu):
    raw, psf, psfb = data_and_psfs([48, 46, 36])
    k = G.LucyRFFTGPU(raw.shape)
    ref = bits(k.omw(raw, k.otf(psf), k.otf(psfb), 2, 100))
    s = gpu.cuda.Stream(non_blocking=True)
    with s:
        d = gpu.asarray(raw)
        out = k.omw(d, k.otf(psf), k.otf(psfb), 2, 100)
        e = s.record()
    e.synchronize()
    assert np.array_equal(ref, bits(out))


@pytest.mark.gpu
def test_ptx_compiled_by_the_driver_gives_the_same_bits(gpu, monkeypatch):
    raw, psf, psfb = data_and_psfs([47, 45, 33])
    k = G.LucyRFFTGPU(raw.shape)
    ref = bits(k.omw(raw, k.otf(psf), k.otf(psfb), 2, 100))
    monkeypatch.setenv("PETAKIT5D_GPU_TARGET", "ptx")
    assert np.array_equal(ref, bits(k.omw(raw, k.otf(psf), k.otf(psfb), 2, 100)))
    assert G._module("rl_rfft.cu").info.jit


def _expected_bytes(g: G.Geometry, *, ws: int, kind: str) -> int:
    """the mex's buffers of a call (each rounded to the pool's 512 bytes): 'omw' J, R, the half
    spectrum and the work area; 'rlstep' R, the half spectrum and the work area"""

    def b(n):
        return -(-n // 512) * 512

    real, half = b(4 * g.n), b(8 * g.nh)
    return (2 * real if kind == "omw" else real) + half + b(ws)


@pytest.mark.gpu
@pytest.mark.parametrize("layout", G.LAYOUTS)
@pytest.mark.parametrize("internal", [False, True])
def test_peak_memory_is_the_mexs(gpu, layout, internal):
    # at most two complex half volumes (the spectrum and the work area) besides J and R; uint16
    # input is read as it is (no single copy), the (z, y, x) result takes the spectrum's memory
    from petakit5d.gpu import fft
    from petakit5d.gpu.memory import track_peak

    raw, psf, psfb = data_and_psfs([128, 96, 40])
    k = G.LucyRFFTGPU(raw.shape, layout=layout, keep_plans=False)
    hf, hb = k.otf(psf), k.otf(psfb)
    d = k.to_internal(raw) if internal else gpu.asarray(raw)
    ws = fft.RealFFTPlan(k.geometry.fft_shape).work_bytes
    gpu.cuda.Device().synchronize()
    with track_peak() as p:
        j = k.omw(d, hf, hb, 3, 100, internal=internal)
        gpu.cuda.Device().synchronize()
    assert p.peak_above_start <= _expected_bytes(k.geometry, ws=ws, kind="omw")
    y = j if internal else gpu.asarray(np.maximum(raw.astype(np.float32) - 99.5, np.float32(0)))
    with track_peak() as p:
        k.rlstep(y, d, hf, internal=internal)
        gpu.cuda.Device().synchronize()
    assert p.peak_above_start <= _expected_bytes(k.geometry, ws=ws, kind="rlstep")


@pytest.mark.gpu
def test_plans_are_cached_or_made_per_call(gpu):
    from petakit5d.gpu import fft

    raw, psf, _ = data_and_psfs([48, 46, 36])
    shape = G.geometry(raw.shape).fft_shape
    k = G.LucyRFFTGPU(raw.shape, keep_plans=False)
    k.omw(raw, k.otf(psf), k.otf(psf), 1, 100)
    assert fft.cached_plan_bytes() == 0
    k = G.LucyRFFTGPU(raw.shape)
    h = k.otf(psf)
    assert any(key[1:3] == ("real", shape) for key in fft._plans)
    a = bits(k.omw(raw, h, h, 1, 100))
    b = bits(G.LucyRFFTGPU(raw.shape, keep_plans=False).omw(raw, h, h, 1, 100))
    assert np.array_equal(a, b)


@pytest.mark.gpu
def test_timing_gives_the_mexs_four_parts(gpu):
    raw, psf, _ = data_and_psfs([48, 46, 36])
    k = G.LucyRFFTGPU(raw.shape, timing=True, keep_plans=False)
    h = k.otf(psf)
    k.omw(raw, h, h, 2, 100)
    plan_t, r2c, c2r, elem = k.last_times
    assert min(plan_t, r2c, c2r, elem) > 0
    k.timing = False
    k.omw(raw, h, h, 2, 100)
    assert k.last_times == (0.0, 0.0, 0.0, 0.0)


@pytest.mark.gpu
def test_out_of_memory_leaves_no_gpu_arrays_behind(gpu):
    # DeconRFFTGPUTest/outOfMemoryLeavesNoGpuArraysBehind: an allocation that fails inside a call
    # frees what the call allocated before, with MATLAB's identifier
    from petakit5d.gpu.memory import available_bytes

    shape = (128, 256, 256)
    k = G.LucyRFFTGPU(shape, keep_plans=False)
    im = gpu.asarray(CounterRNG(4).integers(100, 1100, shape, np.uint16))
    h = k.otf(np.full((3, 3, 3), np.float32(1 / 27)))
    pool = gpu.get_default_memory_pool()
    pool.free_all_blocks()
    leave = int(2.5 * 4 * math.prod(shape))
    ballast = []
    try:
        while available_bytes() - leave > 2**30:
            ballast.append(gpu.empty(2**30, gpu.uint8))
        ballast.append(gpu.empty(max(1, available_bytes() - leave), gpu.uint8))
        gc.collect()  # (arrays of earlier tests that wait for the collector)
        before = pool.used_bytes()
        with pytest.raises(G.GPUOutOfMemoryError) as e:
            k.omw(im, h, h, 2, 100)
        assert e.value.identifier == "parallel:gpu:array:OOM"
        assert pool.used_bytes() == before, "GPU memory was left allocated"
    finally:
        ballast.clear()
        pool.free_all_blocks()


@pytest.mark.gpu
def test_plan_memory_errors_have_matlabs_identifiers(gpu, monkeypatch):
    from petakit5d.gpu import fft

    k = G.LucyRFFTGPU((8, 8, 8), keep_plans=False)
    for code, ident in (
        (2, "parallel:gpu:array:OOMForOperation"),
        (5, "parallel:gpu:array:FFTPlanInternalError"),
    ):

        def fail(*a, c=code, **kw):
            raise fft.CuFFTMemoryError(c, "cufftMakePlanMany")

        monkeypatch.setattr(fft, "RealFFTPlan", fail)
        with pytest.raises(G.GPUOutOfMemoryError) as e:
            k.otf(np.ones((3, 3, 3), np.float32))
        assert e.value.identifier == ident


@pytest.mark.gpu
def test_gpu_erosion_equals_the_cpu_kernel(gpu):
    # DeconMaskErosionTest/gpuArrayMasksEqualTheCPUMasks: the GPU mex and the CPU mex (K05's
    # kernel) give the same masks, for the disk decompositions and random neighborhoods
    from petakit5d.decon import erosion as E
    from petakit5d.kernels import decon as K

    for name, ms in MASK_SIZES.items():
        for shape in ("random", "skewed", "blobs"):
            a = from_matlab(make_mask(shape, ms))
            for e in (1, 2, 3, 5, 9, 20):
                nh, pu, pl = E.disk_decomposition(e - 1)
                rz = e - 1 if a.ndim == 3 and a.shape[0] > 1 else 0
                ref = K.mask_edge_erosion(a, nh, pu, pl, rz, True)
                out = G.mask_edge_erosion(gpu.asarray(a), nh, pu, pl, rz, True)
                assert out.dtype == gpu.bool_ and out.shape == a.shape
                assert np.array_equal(out.get(), ref), (name, shape, e)
                assert np.array_equal(gpu_erosion(a, e), E.decon_mask_edge_erosion(a, e))
    for t in range(16):
        am, nh, pu, pl = random_nhoods(t)
        a = from_matlab(am)
        for rz, clear in ((0, False), (2, True)):
            ref = K.mask_edge_erosion(a, nh, pu, pl, rz, clear)
            assert np.array_equal(G.mask_edge_erosion(a, nh, pu, pl, rz, clear).get(), ref), t


@pytest.mark.gpu
def test_gpu_erosion_wrapper_equals_the_cpu_wrapper(gpu):
    # full masks (also thinner than the erosion), numeric masks, erosions the mex doesn't take
    from petakit5d.decon import erosion as E

    for shape in ((3, 20, 20), (20, 3, 20), (20, 20, 3), (20, 1, 30), (3, 20), (30, 1)):
        for e in (1, 2, 3, 2.5, NAN):
            a = np.ones(shape, bool)
            try:
                ref = E.decon_mask_edge_erosion(a, e)
            except E.MatlabError as err:
                # beyond a full mask's size, MATLAB's gpuArray assignment errors with its own id
                assert err.identifier == "MATLAB:badsubscript"
                with pytest.raises(E.MatlabError) as g_err:
                    gpu_erosion(a, e)
                assert g_err.value.identifier == "parallel:array:SubsasgnInvalidSubscript"
                continue
            assert np.array_equal(gpu_erosion(a, e), ref), (shape, e)
    a = from_matlab(make_mask("random", [12, 14, 6])).astype(np.float32) * np.float32(3) - 1
    for e in (0, 4, 2.5):
        try:
            ref = E.decon_mask_edge_erosion(a, e)
        except E.MatlabError as err:
            with pytest.raises(E.MatlabError) as g_err:
                gpu_erosion(a, e)
            assert g_err.value.identifier == err.identifier
            continue
        assert np.array_equal(gpu_erosion(a, e), ref), e
