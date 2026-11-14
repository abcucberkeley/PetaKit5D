"""Deconvolution of large files (petakit5d.decon.large, decon.block) against MATLAB's
RLdecon_large_in_memory, RLdecon_large_in_place and RLdecon_for_zarr_block (the oracle tree's,
on the CPU and on an A100).

Oracle cases: volumes of blurred beads (MATLAB 61 x 55 x 37, background 100) cut into 2 x 2 x 2
batches with odd remainders, a Gaussian PSF of 31^3 (and a noisy one for OMW); the files each
side writes are compared (results, MIPs and MIP slabs, the generated PSF, masks): the in-memory
decon with erosion (data without zeros, and with zeros: the distance map), a saved mask, a mask
file, zarr input and output, single results, OMW, 2D masks that skip batches, psfGen false,
flipZstack, the options, the default background, and Deskew and Rotate (the batches are only
deconvolved: Deskew sets the border PSF's z step); the
in-place decon (through MATLAB's local job framework) with the same options, a deskewed PSF
border, an existing result and options that change on their way through MATLAB's task strings;
the block task with useGPU on a node without a GPU and on the CPU. Each runs on both CPU FFT
backends (``.exact``: MATLAB's bits; ``.permissive``: decision record 0006).

The pytest functions port DeconWrapperTest's large-file methods (and DeconFunctionTest's
gpuBlockTaskWithoutUsableGPUErrors) in intent, through ``decon_file`` (XR_RLdeconFrame3D) and
the large-file functions, with spies for the framework, the MIPs and the moves.
"""

from __future__ import annotations

import math
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.decon import large as L
from petakit5d.decon import rl
from petakit5d.testing import oracle, slurm_fake
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.compare import tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, beads, decon_psf, write_tiff, write_zarr

BACKENDS = ("exact", "permissive")
XY, DZ, BG = 0.108, 0.2, 100
IM = (37, 61, 55)  # (z, y, x): MATLAB's 61 x 55 x 37
BLUR = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)
FSN = (
    "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_"
    "000x_000y_000z_0000t"
)
# MATLAB's [y x z]: 2 x 2 x 2 batches with odd remainders (61 = 32 + 29, 55 = 30 + 25 or
# 32 + 23, 37 = 20 + 17 or 24 + 13)
MEM_BATCH = [32, 30, 20]
PLACE_BATCH = [32, 32, 24]
PLACE_BLOCK = [16, 16, 8]

PERMISSIVE_U16 = tol(
    reason="0006: uint16 deconvolution with the permissive FFT backend (at most 1 count in at "
    "most 0.1% of the voxels)",
    abs=1,
    max_frac=1e-3,
)
PERMISSIVE_TOL = tol(
    reason="0006: deconvolution with the permissive FFT backend (relative L2 <= 1e-4)", rel_l2=1e-4
)
PERMISSIVE_BP = tol(
    reason="accepted (R02, 2026-11-13): the back projector and its masked Wiener filter with the "
    "permissive backend: max <= 1e-5 of the largest value, relative L2 <= 1e-5 (masks identical)",
    rel_max=1e-5,
    rel_l2=1e-5,
)
PERMISSIVE_FFT = tol(
    reason="0006: FFT-based outputs with the permissive backend, single precision (max <= 1e-5 "
    "of the largest value, relative L2 <= 1e-6)",
    rel_max=1e-5,
    rel_l2=1e-6,
)


# ---------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------


def frame(seed: int, shape=IM, n=14) -> np.ndarray:
    """blurred beads on a background of 100 (uint16, (z, y, x))"""
    rng = CounterRNG(seed, "large_frame")
    nz, ny, nx = shape
    u = rng.uniform((n, 3))
    pos = np.column_stack([2 + (nz - 4) * u[:, 0], 4 + (ny - 8) * u[:, 1], 4 + (nx - 8) * u[:, 2]])
    im = beads(shape, pos, BLUR, 1000.0, float(BG))
    return np.floor(im + np.float32(0.5)).astype(np.uint16)


def zero_corner(a: np.ndarray) -> np.ndarray:
    """a slanted corner without data (DeconWrapperTest's y + 2 x + 3 z < 60, 1-based)"""
    nz, ny, nx = a.shape
    z = np.arange(1, nz + 1, dtype=np.float64).reshape(-1, 1, 1)
    y = np.arange(1, ny + 1, dtype=np.float64).reshape(1, -1, 1)
    x = np.arange(1, nx + 1, dtype=np.float64).reshape(1, 1, -1)
    return np.where(y + 2 * x + 3 * z < 60, np.uint16(0), a)


def psf_raw() -> np.ndarray:
    return np.floor(decon_psf((31, 31, 31), (2.0, 1.5)) + np.float32(0.5)).astype(np.uint16)


def psf_noisy() -> np.ndarray:
    """the PSF with 5% deterministic noise in its nonzero slices (no exactly tied |OTF| values
    at the OMW mask's threshold; see R05's test_frame)"""
    p = decon_psf((31, 31, 31), (2.0, 1.5))
    u = CounterRNG(21, "psf_noise").uniform(p.shape, np.float32)
    p = p * (np.float32(1) + np.float32(0.05) * u)
    return np.floor(p + np.float32(0.5)).astype(np.uint16)


def write_masks(folder, shape, y_end, layouts=("matlab", "matlab", "matlab")) -> list[str]:
    """the 2D masks mask_{xy,xz,yz}.zarr of a volume of shape (z, y, x): MATLAB's ny x nx,
    nx x nz and ny x nz images (one-slice volumes, as MATLAB's createzarr stores images), the xy
    mask zero from y_end on (0-based)"""
    nz, ny, nx = shape
    xy = np.ones((1, ny, nx), np.uint8)
    xy[:, y_end:] = 0
    masks = [xy, np.ones((1, nx, nz), np.uint8), np.ones((1, ny, nz), np.uint8)]
    out = []
    for name, m, lay in zip(("xy", "xz", "yz"), masks, layouts, strict=True):
        fn = os.path.join(folder, f"mask_{name}.zarr")
        write_zarr(fn, m, (1, 32, 32), lay)
        out.append(fn)
    return out


def make_root(*, data=None, zarr=False, noisy=False, masks=False, mask_file=False, existing=False):
    """data/<FSN>.tif (or .zarr in MATLAB's layout), psf/PSF.tif; masks: 2D masks
    masks/mask_{xy,xz,yz}.zarr (MATLAB's layout; the xy mask zero from y 32 on, 0-based);
    mask_file: masks/eroded.tif (uint8: 0 at the faces, 2 in one corner, 1 elsewhere);
    existing: a result decon/<FSN>.zarr there"""

    def make(folder):
        folder = os.fspath(folder)
        os.makedirs(os.path.join(folder, "data"))
        im = frame(1) if data is None else data(1)
        if zarr:
            write_zarr(os.path.join(folder, "data", FSN + ".zarr"), im, (8, 16, 16), "matlab")
        else:
            write_tiff(os.path.join(folder, "data", FSN + ".tif"), im)
        os.makedirs(os.path.join(folder, "psf"))
        write_tiff(os.path.join(folder, "psf", "PSF.tif"), psf_noisy() if noisy else psf_raw())
        if masks or mask_file:
            os.makedirs(os.path.join(folder, "masks"))
        if masks:
            write_masks(os.path.join(folder, "masks"), IM, 32)
        if mask_file:
            mf = np.ones(IM, np.uint8)
            mf[:2] = 0
            mf[-2:] = 0
            mf[:, :3] = 0
            mf[:, :, -4:] = 0
            mf[10:20, 30:40, 5:15] = 2
            write_tiff(os.path.join(folder, "masks", "eroded.tif"), mf)
        if existing:
            os.makedirs(os.path.join(folder, "decon"))
            write_zarr(
                os.path.join(folder, "decon", FSN + ".zarr"), frame(9), (8, 16, 16), "matlab"
            )

    return make


# ---------------------------------------------------------------------------------------------
# the cases
# ---------------------------------------------------------------------------------------------

# option values with '<R>/' are paths below the case's root folder
_ROOTED = (
    "a = in(7 : end); for i = 1 : numel(a), if ischar(a{i}) || iscell(a{i}), "
    "a{i} = strrep(a{i}, '<R>/', d); end; end; "
)
_MEM_CODE = (
    "d = in{1}; " + _ROOTED + "RLdecon_large_in_memory([d, in{2}], [d, in{3}], [d, in{4}], "
    "in{5}, in{6}, a{:}); out = {};"
)
_PLACE_CODE = (
    "d = in{1}; " + _ROOTED + "RLdecon_large_in_place([d, in{2}], in{5}, in{6}, [d, in{4}], "
    "[d, in{3}], a{:}); out = {};"
)
_BLOCK_CODE = (
    "d = in{1}; " + _ROOTED + "outFn = [d, 'decon/out.zarr']; mkdir([d, 'decon/flags']); "
    "createzarr(outFn, dataSize=[61, 55, 37], blockSize=[16, 16, 8], dtype='uint16'); "
    "RLdecon_for_zarr_block(1, [d, in{2}], [d, in{3}], outFn, [d, in{4}], "
    "[1, 1, 1, 61, 55, 37], [1, 1, 1, 61, 55, 37], in{5}, in{6}, a{:}); out = {};"
)


def _rooted(opts, d: str) -> list:
    out = []
    for v in opts:
        if isinstance(v, str):
            v = v.replace("<R>/", d)
        elif isinstance(v, (list, tuple)) and all(isinstance(x, str) for x in v):
            v = [x.replace("<R>/", d) for x in v]
        out.append(v)
    return out


def p_mem(d, frame_fn, psf_fn, out_fn, xy, dz, *opts, backend):
    from petakit5d import matlab

    rl.clear_otf_cache()
    with P.using(backend):
        matlab.RLdecon_large_in_memory(
            d + frame_fn, d + psf_fn, d + out_fn, xy, dz, *_rooted(opts, d)
        )
    return []


def p_place(d, frame_fn, psf_fn, decon_dir, xy, dz, *opts, backend):
    from petakit5d import matlab

    rl.clear_otf_cache()
    with P.using(backend):
        matlab.RLdecon_large_in_place(
            d + frame_fn, xy, dz, d + decon_dir, d + psf_fn, *_rooted(opts, d)
        )
    return []


def p_block(d, frame_fn, psf_fn, flag_fn, xy, dz, *opts, backend):
    from petakit5d import matlab
    from petakit5d.io import create_zarr

    out_fn = d + "decon/out.zarr"
    os.makedirs(d + "decon/flags", exist_ok=True)
    create_zarr(out_fn, data_size=IM, block_size=(8, 16, 16), dtype="uint16")
    rl.clear_otf_cache()
    box = [1, 1, 1, 61, 55, 37]
    with P.using(backend):
        matlab.RLdecon_for_zarr_block(
            1, d + frame_fn, d + psf_fn, out_fn, d + flag_fn, box, box, xy, dz, *_rooted(opts, d)
        )
    return []


# (R02: MATLAB also saves a figure of the back projector, which the port does not draw)
FIGURE = r".*_figure\.png$"


def large_tree(single=False):
    return tree(
        ignore=[FIGURE],
        per_file={
            "*_back_projector_*": PERMISSIVE_BP,
            "*_OTF.tif": PERMISSIVE_FFT,
            "*/psfgen/*": PERMISSIVE_FFT,
            "*": PERMISSIVE_TOL if single else PERMISSIVE_U16,
        },
        tiff_bytes=False,
    )


MEM_BASE = [
    "Background",
    BG,
    "dzPSF",
    DZ,
    "DeconIter",
    2,
    "RLMethod",
    "simplified",
    "EdgeErosion",
    0,
    "batchSize",
    MEM_BATCH,
    "useGPU",
    False,
]
PLACE_BASE = [
    "Background",
    BG,
    "dzPSF",
    DZ,
    "DeconIter",
    2,
    "RLMethod",
    "simplified",
    "batchSize",
    PLACE_BATCH,
    "blockSize",
    PLACE_BLOCK,
    "parseCluster",
    False,
    "jobLogDir",
    "<R>/logs",
]


def _without(base, name):
    """the base options without a name-value pair"""
    i = base.index(name)
    return base[:i] + base[i + 2 :]


def _cases(
    cid,
    kind,
    opts,
    *,
    root=None,
    base=None,
    zarr_input=False,
    single=False,
    expect_error=False,
    backends=BACKENDS,
    device="cpu",
):
    if kind == "mem":
        code, python, defaults = _MEM_CODE, p_mem, MEM_BASE
        frame_fn = "data/" + FSN + (".zarr" if zarr_input else ".tif")
        out = [
            frame_fn,
            "psf/PSF.tif",
            "decon/" + FSN + (".zarr" if "saveZarr" in opts else ".tif"),
        ]
    else:
        code, python, defaults = _PLACE_CODE, p_place, PLACE_BASE
        out = ["data/" + FSN + ".zarr", "psf/PSF.tif", "decon"]
    res = []
    for be in backends:
        res.append(
            Case(
                f"{cid}.{be}",
                device=device,
                code=code,
                args=[InPath("root"), *out, XY, DZ, *(defaults if base is None else base), *opts],
                nargout=0,
                python=python,
                py_kwargs={"backend": be},
                inputs={"root": Tree(root or MEM_ROOT, "root")},
                expect_error=expect_error,
                py_error=Exception,
                tree=large_tree(single) if be == "permissive" else tree(ignore=[FIGURE]),
            )
        )
    return res


MEM_ROOT = make_root()
ZARR_ROOT = make_root(zarr=True)
ZERO_ROOT = make_root(data=lambda s: zero_corner(frame(s)))
NOISY_ROOT = make_root(noisy=True)
NOISY_ZARR_ROOT = make_root(zarr=True, noisy=True)
MASK_ROOT = make_root(masks=True)
MASK_ZARR_ROOT = make_root(zarr=True, masks=True)
MASKFILE_ROOT = make_root(mask_file=True)
EXISTING_ROOT = make_root(zarr=True, existing=True)
_OMW = ["RLMethod", "omw", "skewed", False, "wienerAlpha", 0.005, "OTFCumThresh", 0.9]

CASES: list[Case] = []
# --- in memory: DeconWrapperTest's largeFile* methods on 2 x 2 x 2 batches
CASES += _cases("mem.simplified", "mem", [])
CASES += _cases("mem.erosion", "mem", ["EdgeErosion", 3])
CASES += _cases(
    "mem.erosion_zeros", "mem", ["EdgeErosion", 3, "SaveMaskfile", True], root=ZERO_ROOT
)
CASES += _cases("mem.erosion_saved", "mem", ["EdgeErosion", 3, "SaveMaskfile", True])
CASES += _cases(
    "mem.erode_maskfile",
    "mem",
    ["EdgeErosion", 3, "ErodeMaskfile", "<R>/masks/eroded.tif"],
    root=MASKFILE_ROOT,
)
CASES += _cases("mem.zarr", "mem", ["saveZarr", True, "blockSize", [16, 16, 8]])
CASES += _cases("mem.zarr_default_chunks", "mem", ["saveZarr", True])
CASES += _cases("mem.zarr_input", "mem", ["EdgeErosion", 2], root=ZARR_ROOT, zarr_input=True)
CASES += _cases("mem.single", "mem", ["save16bit", False], single=True)
CASES += _cases("mem.omw_noisy", "mem", _OMW, root=NOISY_ROOT)
CASES += _cases(
    "mem.masks",
    "mem",
    [
        "maskFullpaths",
        ["<R>/masks/mask_xy.zarr", "<R>/masks/mask_xz.zarr", "<R>/masks/mask_yz.zarr"],
    ],
    root=MASK_ROOT,
)
CASES += _cases("mem.psfgen_false", "mem", ["psfGen", False])
CASES += _cases("mem.flip", "mem", ["flipZstack", True, "EdgeErosion", 2])
CASES += _cases(
    "mem.options", "mem", ["dampFactor", 1.2, "scaleFactor", 2, "deconOffset", 5, "fixIter", False]
)
CASES += _cases("mem.background_default", "mem", ["Background", []])
CASES += _cases("mem.batch_one", "mem", ["batchSize", [64, 64, 64], "EdgeErosion", 3])
# the batches are only deconvolved (until c24801c each batch's RLdecon skipped with Deskew, and
# its empty result couldn't be placed, and failed on its undefined file name with Rotate)
CASES += _cases("mem.deskew", "mem", ["Deskew", True])
CASES += _cases("mem.rotate", "mem", ["Rotate", True])
# --- in place: through MATLAB's local job framework (tasks of 5 batches, flags with the
# blocks' MIPs, the MIPs merged from them)
CASES += _cases("place.simplified", "place", [], root=ZARR_ROOT)
CASES += _cases("place.erosion", "place", ["EdgeErosion", 3], root=ZARR_ROOT)
CASES += _cases(
    "place.masks",
    "place",
    [
        "maskFullpaths",
        ["<R>/masks/mask_xy.zarr", "<R>/masks/mask_xz.zarr", "<R>/masks/mask_yz.zarr"],
    ],
    root=MASK_ZARR_ROOT,
)
CASES += _cases("place.omw_noisy", "place", _OMW, root=NOISY_ZARR_ROOT)
CASES += _cases("place.single", "place", ["save16bit", False], root=ZARR_ROOT, single=True)
CASES += _cases("place.deskew", "place", ["Deskew", True, "EdgeErosion", 2], root=ZARR_ROOT)
# (dz sind(-32.45) / dzPSF < -1: no resampling of the border PSF, where |angle| would)
CASES += _cases(
    "place.deskew_negative_ratio",
    "place",
    ["Deskew", True, "dzPSF", 0.05],
    base=_without(PLACE_BASE, "dzPSF"),
    root=ZARR_ROOT,
)
CASES += _cases("place.psfgen_false", "place", ["psfGen", False], root=ZARR_ROOT)
CASES += _cases("place.flip_ignored", "place", ["flipZstack", True], root=ZARR_ROOT)
# (the blocks get 1.234568 and 2.5 through MATLAB's '%d' of the task string)
CASES += _cases(
    "place.options_text",
    "place",
    ["dampFactor", 1.2345678901, "scaleFactor", 2.5, "deconOffset", 3],
    root=ZARR_ROOT,
)
CASES += _cases("place.existing", "place", [], root=EXISTING_ROOT)
CASES += _cases("place.one_batch", "place", ["batchSize", [64, 64, 64]], root=ZARR_ROOT)
# --- the block task: a GPU job's block without a usable GPU fails; on the CPU it runs
_BLOCK_BASE = ["Background", BG, "dzPSF", DZ, "DeconIter", 2, "psfGen", False]
for _be in BACKENDS:
    CASES.append(
        Case(
            f"block.cpu.{_be}",
            code=_BLOCK_CODE,
            args=[
                InPath("root"),
                "data/" + FSN + ".zarr",
                "psf/PSF.tif",
                "decon/flags/blocks_1_1.mat",
                XY,
                DZ,
                *_BLOCK_BASE,
                "useGPU",
                False,
                "mipBatchSize",
                [32, 32, 24],
            ],
            nargout=0,
            python=p_block,
            py_kwargs={"backend": _be},
            inputs={"root": Tree(ZARR_ROOT, "root")},
            tree=large_tree() if _be == "permissive" else tree(ignore=[FIGURE]),
        )
    )
CASES.append(
    Case(
        "error.block.gpu_without_gpu",
        code=_BLOCK_CODE,
        args=[
            InPath("root"),
            "data/" + FSN + ".zarr",
            "psf/PSF.tif",
            "decon/flags/blocks_1_1.mat",
            XY,
            DZ,
            *_BLOCK_BASE,
            "useGPU",
            True,
        ],
        nargout=0,
        python=p_block,
        py_kwargs={"backend": "exact"},
        inputs={"root": Tree(ZARR_ROOT, "root")},
        expect_error="decon:noUsableGPU",
        py_error=Exception,
        tree=tree(ignore=[FIGURE]),
    )
)
# --- GPU (A100): the in-memory decon with useGPU and an in-place GPU job (MATLAB's gpuArray
# bits through R04's GPU path)
CASES += _cases(
    "gpu.mem.simplified",
    "mem",
    ["useGPU", True],
    base=_without(MEM_BASE, "useGPU"),
    backends=("exact",),
    device="gpu:A100",
)
CASES += _cases(
    "gpu.place.simplified",
    "place",
    ["GPUJob", True],
    root=ZARR_ROOT,
    backends=("exact",),
    device="gpu:A100",
)


def _gpu_name() -> str:
    from petakit5d.gpu.device import select_device

    return select_device().name


def _check(case, tmp_path):
    if case.device.startswith("gpu:"):
        name = _gpu_name()
        if case.device[4:].lower() not in name.lower():
            pytest.skip(f"the reference is MATLAB's on a {case.device[4:]}; this GPU is a {name}")
    if ".exact" in case.id or case.id.startswith("error.block"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        oracle.check(case, tmp_path)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    _check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# helpers of the tests below
# ---------------------------------------------------------------------------------------------

SMALL = (32, 48, 48)  # (z, y, x): DeconWrapperTest's 48 x 48 x 32
SMALL_FSN = FSN


def same_bits(a, b) -> None:
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, a.shape, b.dtype, b.shape)
    assert np.array_equal(
        np.ascontiguousarray(a).view(np.uint8), np.ascontiguousarray(b).view(np.uint8)
    )


def close(a, b, rel=0.01, abs_=2.0) -> None:
    """DeconWrapperTest's verifyImagesClose"""
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    assert np.all(np.abs(a - b) <= np.maximum(abs_, rel * np.abs(b))), np.max(np.abs(a - b))


def dataset(tmp_path, *, zarr=False, data=None, n=1) -> tuple[str, list[str], str]:
    """DeconWrapperTest's createDataset: n time points of 48 x 48 x 32 beads (tiff, or zarr in
    MATLAB's layout) and a PSF; (data folder, file paths, PSF path)"""
    d = tmp_path / "data"
    d.mkdir(parents=True)
    fns = []
    for t in range(n):
        im = frame(10 + t, SMALL, 6) if data is None else data(t)
        name = SMALL_FSN.replace("_0000t", f"_{t:04d}t")
        if zarr:
            fn = str(d / (name + ".zarr"))
            write_zarr(fn, im, (16, 16, 16), "matlab")
        else:
            fn = str(d / (name + ".tif"))
            write_tiff(fn, im)
        fns.append(fn)
    psf = str(tmp_path / "PSF.tif")
    write_tiff(psf, psf_raw())
    return str(d) + "/", fns, psf


def run_frame(fn, decon_path, psf, **kw):
    """decon_file (XR_RLdeconFrame3D) as the wrapper runs it: the CPU, no cluster"""
    from petakit5d.decon.frame import decon_file

    base = {
        "psf_file": psf,
        "dz_psf": DZ,
        "background": BG,
        "rl_method": "simplified",
        "edge_erosion": 0,
        "parse_cluster": False,
        "gpu_job": False,
    }
    base.update(kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        decon_file(fn, XY, DZ, decon_path, **base)


def read(p):
    from petakit5d.io import read_tiff, read_zarr

    p = os.fspath(p)
    return read_zarr(p) if p.rstrip("/").endswith(".zarr") else read_tiff(p)


def verify_large_file_mips(result: str) -> None:
    """pktest.verifyLargeFileMIPs: the three MIPs (MATLAB's orientation) and the slabs' MIPs are
    those of the result"""
    out = read(result)
    folder, name = os.path.split(result.rstrip("/"))
    name = name[: -len(".zarr")]
    want = {"y": out.max(axis=1).T, "x": out.max(axis=2).T, "z": out.max(axis=0)}
    for lab, m in want.items():
        same_bits(read(f"{folder}/MIPs/{name}_MIP_{lab}.tif"), m)
        slab = np.asarray(read(f"{folder}/MIPs/{name}_{lab}.zarr"))
        slab = slab.reshape((-1,) * (3 - slab.ndim) + slab.shape) if slab.ndim < 3 else slab
        axis = {"y": 1, "x": 2, "z": 0}[lab]
        same_bits(np.asarray(slab).max(axis=axis), out.max(axis=axis))


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


@pytest.fixture
def mip_spy(monkeypatch):
    """records mip_zarr calls (FunctionSpyFixture('XR_MIP_zarr')) and runs them; ``spy.action``
    runs first"""
    from petakit5d.tools import mip

    real = mip.mip_zarr
    calls = []

    def spy(*a, **k):
        calls.append((a, k))
        if spy.action is not None:
            return spy.action()
        return real(*a, **k)

    spy.action = None
    spy.calls = calls
    monkeypatch.setattr(mip, "mip_zarr", spy)
    return spy


@pytest.fixture
def no_mip_blocks(monkeypatch):
    """mip_block tasks fail at once (they would read the result again)"""
    from petakit5d.tools import mip

    def no_blocks(*a, **k):
        raise AssertionError("The MIPs were computed from the result.")

    monkeypatch.setattr(mip, "mip_block", no_blocks)


IN_PLACE = {
    "rl_method": "simplified",
    "decon_iter": 2,
    "large_file": True,
    "large_method": "inplace",
    "save_zarr": True,
    "batch_size": (64, 64, 64),
    "block_size": (16, 16, 16),
}


# ---------------------------------------------------------------------------------------------
# DeconWrapperTest in intent: the in-memory decon
# ---------------------------------------------------------------------------------------------


def test_large_file_single_batch_matches_regular(tmp_path):
    # with one batch the in-memory decon runs the same decon on the same data, including the
    # edge erosion
    d, fns, psf = dataset(tmp_path)
    kw = {"decon_iter": 5, "edge_erosion": 3}
    run_frame(fns[0], d + "regular", psf, **kw)
    run_frame(fns[0], d + "large", psf, large_file=True, batch_size=(64, 64, 64), **kw)
    large = read(d + "large/" + SMALL_FSN + ".tif")
    same_bits(large, read(d + "regular/" + SMALL_FSN + ".tif"))
    same_bits(read(d + "large/MIPs/" + SMALL_FSN + "_MIP_z.tif"), large.max(axis=0))


def test_large_file_multiple_batches_match_regular(tmp_path):
    # batches are deconvolved with a border of half the PSF; with 2 iterations the stitched
    # result matches the whole-volume decon away from the image border
    d, fns, psf = dataset(tmp_path)
    run_frame(fns[0], d + "regular", psf, decon_iter=2)
    run_frame(fns[0], d + "large", psf, decon_iter=2, large_file=True, batch_size=(16, 24, 24))
    regular = read(d + "regular/" + SMALL_FSN + ".tif")
    large = read(d + "large/" + SMALL_FSN + ".tif")
    b = 8
    close(large[b:-b, b:-b, b:-b], regular[b:-b, b:-b, b:-b])
    assert not np.array_equal(large, regular), "the batches were not deconvolved apart"


def test_large_file_makes_mip_before_publishing(tmp_path):
    # the wrapper counts a file as done once its result exists
    from petakit5d.testing.publish_spy import PublishSpy

    d, fns, psf = dataset(tmp_path)
    with PublishSpy() as spy:
        run_frame(fns[0], d + "large", psf, decon_iter=2, large_file=True, batch_size=(64,) * 3)
    fn = d + "large/" + SMALL_FSN + ".tif"
    assert spy.was_moved(fn), "The result was not moved into place."
    published = spy.files_when_moved(fn)
    assert f"MIPs/{SMALL_FSN}_MIP_z.tif" in published, "The MIP was written after the move."
    from petakit5d.testing.publish_spy import list_files

    added = set(list_files(d + "large")) - set(published) - {SMALL_FSN + ".tif"}
    assert not added, f"Files were written after the result was moved into place: {added}"


def test_large_file_erosion_mask_matches_distance_map(tmp_path):
    # the eroded mask of the whole volume: without zeros the volume without EdgeErosion voxels
    # at each face (no distance map), with zeros the distance map
    from petakit5d.matlabcompat.distance import bwdist

    d, fns, psf = dataset(tmp_path)
    im = read(fns[0])
    zeros_fn = fns[0].replace("_0000t", "_0001t")
    write_tiff(zeros_fn, zero_corner(im))
    e = 3
    decon = str(tmp_path / "large") + "/"
    for fn in (fns[0], zeros_fn):
        name = os.path.basename(fn)[:-4]
        L.decon_large_in_memory(
            fn,
            psf,
            decon + name + ".tif",
            XY,
            DZ,
            background=BG,
            dz_psf=DZ,
            decon_iter=1,
            edge_erosion=e,
            save_maskfile=True,
            batch_size=(64, 64, 64),
            use_gpu=False,
        )
        bw = read(fn) > 0
        bw[:, [0, -1], :] = False
        bw[:, :, [0, -1]] = False
        bw[[0, -1], :, :] = False
        expected = bwdist(~bw) > e - 1
        mask = read(f"{decon}Masks/{name}_eroded.tif") > 0
        np.testing.assert_array_equal(mask, expected, name)
        # the result is zero outside the mask
        assert not np.any(read(decon + name + ".tif")[~expected])


def test_eroded_mask_without_zeros_skips_the_distance_map(monkeypatch):
    from petakit5d.matlabcompat import distance

    def no_bwdist(*a, **k):
        raise AssertionError("bwdist was used for data without zeros")

    im = frame(3, (12, 20, 18), 3)
    monkeypatch.setattr(distance, "bwdist", no_bwdist)
    m = L.eroded_volume_mask(im, 3)
    want = np.zeros(im.shape, bool)
    want[3:-3, 3:-3, 3:-3] = True
    np.testing.assert_array_equal(m, want)


def test_eroded_mask_of_an_image_is_empty():
    # MATLAB's 2D mask: its only plane is its first and last, so nothing is kept
    im = frame(3, (1, 20, 18), 3)
    assert not np.any(L.eroded_volume_mask(im, 2))
    assert not np.any(L.eroded_volume_mask(im[0], 2))
    # (bwdist of nothing is 0 everywhere: an erosion below 1 keeps it all)
    assert np.all(L.eroded_volume_mask(im, 0.4))


def test_large_file_in_memory_uses_block_size(tmp_path):
    # the in-memory zarr result had writezarr's default chunks (500^3) instead of blockSize
    from petakit5d.io import get_zarr_info

    d, fns, psf = dataset(tmp_path, zarr=True)
    run_frame(
        fns[0],
        d + "large",
        psf,
        decon_iter=1,
        save_zarr=True,
        large_file=True,
        large_method="inmemory",
        batch_size=(64, 64, 64),
        block_size=(16, 16, 16),
    )
    info = get_zarr_info(d + "large/" + SMALL_FSN + ".zarr")
    assert tuple(info.block_size) == (16, 16, 16)


def test_large_file_with_erode_by_first_time_point_completes(tmp_path):
    # (d83468d) the large-file path writes no mask of the first time point (the other files
    # don't wait for one) and erodes each file itself
    d, fns, psf = dataset(tmp_path, n=2)
    e = 3
    for fn in fns:
        run_frame(
            fn,
            "",
            psf,
            decon_iter=2,
            edge_erosion=e,
            save_maskfile=True,
            large_file=True,
            batch_size=(64, 64, 64),
        )
    assert not os.path.exists(d + "matlab_decon/Masks")
    shell = np.ones(SMALL, bool)
    shell[e:-e, e:-e, e:-e] = False
    for fn in fns:
        out = read(d + "matlab_decon/" + os.path.basename(fn))
        assert not np.any(out[shell]), fn
        assert np.any(out[~shell]), fn


def test_large_file_in_memory_batches_are_only_deconvolved(tmp_path, monkeypatch):
    """(MATLAB: largeFileInMemoryBatchesAreOnlyDeconvolved) Deskew and Rotate reached the RLdecon
    of each batch: with Deskew each batch returned nothing to place
    (MATLAB:subsassigndimmismatch), with Rotate each failed on its undefined file name. The
    batches are only deconvolved, as the in-place blocks: Deskew only changes the PSF that sets
    the batches' borders, so with one batch the results are equal"""
    from petakit5d.decon import file as F

    d, fns, psf = dataset(tmp_path)
    kw = {"background": BG, "dz_psf": DZ, "decon_iter": 2, "edge_erosion": 0, "use_gpu": False,
          "rl_method": "simplified"}  # fmt: skip
    seen = []
    real = F.rl_decon_file

    def spy(*a, **k):
        seen.append((k["deskew"], k["rotate"]))
        return real(*a, **k)

    monkeypatch.setattr(F, "rl_decon_file", spy)
    runs = {
        "plain": {"batch_size": (16, 24, 24)},
        "rotate": {"batch_size": (16, 24, 24), "rotate": True},
        "one": {"batch_size": (64, 64, 64)},
        "deskew": {"batch_size": (64, 64, 64), "deskew": True},
    }
    out = {}
    for name, opts in runs.items():
        fn = f"{d}large/{name}/{os.path.basename(fns[0])}"
        L.decon_large_in_memory(fns[0], psf, fn, XY, DZ, **kw, **opts)
        out[name] = read(fn)
    same_bits(out["rotate"], out["plain"])
    same_bits(out["deskew"], out["one"])
    assert out["plain"].shape == SMALL and out["plain"].max() > 0, "not deconvolved"
    assert seen and set(seen) == {(False, False)}
    for name in ("rotate", "deskew"):
        assert not os.path.isdir(f"{d}large/{name}/DSR"), f"{name}: a batch was rotated"


def test_in_memory_masks_skip_empty_batches(tmp_path, monkeypatch):
    from petakit5d.decon import file as F

    d, fns, psf = dataset(tmp_path)
    nz, ny, nx = SMALL
    xy = np.ones((ny, nx), np.uint8)
    xy[24:] = 0
    masks = [str(tmp_path / f"m_{k}.zarr") for k in ("xy", "xz", "yz")]
    write_zarr(masks[0], xy, (16, 16))
    write_zarr(masks[1], np.ones((nx, nz), np.uint8), (16, 16))
    write_zarr(masks[2], np.ones((ny, nz), np.uint8), (16, 16), "matlab")
    real = F.rl_decon_file
    calls = []

    def spy(*a, **k):
        calls.append(k["decon_bbox"])
        return real(*a, **k)

    monkeypatch.setattr(F, "rl_decon_file", spy)
    out_fn = d + "x/out.tif"
    L.decon_large_in_memory(
        fns[0],
        psf,
        out_fn,
        XY,
        DZ,
        background=BG,
        dz_psf=DZ,
        decon_iter=2,
        edge_erosion=0,
        batch_size=(16, 24, 24),
        mask_fullpaths=masks,
        use_gpu=False,
    )
    out = read(out_fn)
    assert not np.any(out[:, 24:]), "The batches outside the masks were deconvolved."
    assert np.any(out[:, :24])
    assert len(calls) == 4  # of 8 batches


def test_a_region_with_some_mask_zeros_is_deconvolved(tmp_path):
    nz, ny, nx = SMALL
    masks = write_masks(str(tmp_path), SMALL, ny)
    from petakit5d.io import write_zarr as wz

    xz = np.ones((1, nx, nz), np.uint8)
    xz[0, 5, 7] = 0
    wz(xz, masks[1], block_size=(1, 16, 16), overwrite=True)
    assert not L.masks_empty(masks, (0, 0, 0, nz, ny, nx))
    xz[...] = 0
    xz[0, 40, 20] = 1
    wz(xz, masks[1], block_size=(1, 16, 16), overwrite=True)
    assert L.masks_empty(masks, (0, 0, 0, 16, 24, 24))
    assert not L.masks_empty(masks, (16, 0, 24, 32, 24, 48))


def test_in_memory_batches_of_another_class_are_cast_to_the_result(tmp_path, monkeypatch):
    # (RLMethod 'original' gives single batches): cast as MATLAB's cast, rounded half away
    d, fns, psf = dataset(tmp_path)

    def fake(rawdata, out_fn, psf_fn, xy, dz, dz_psf, box, kw):
        shape = tuple(box[q + 3] - box[q] for q in range(3))
        return np.full(shape, 2.5, np.float32)

    monkeypatch.setattr(L, "_rl_batch", fake)
    out_fn = d + "x/out.tif"
    L.decon_large_in_memory(
        fns[0], psf, out_fn, XY, DZ, background=BG, edge_erosion=0, batch_size=(16, 24, 24)
    )
    same_bits(read(out_fn), np.full(SMALL, 3, np.uint16))


def test_in_place_background_defaults_to_100(tmp_path, framework_spy, mip_spy):
    d, fns, psf = dataset(tmp_path, zarr=True)
    framework_spy.done = True
    mip_spy.action = lambda: None
    L.decon_large_in_place(
        fns[0], XY, DZ, d + "large", psf, parse_cluster=False, batch_size=(64, 64, 64)
    )
    assert framework_spy.calls[0]["tasks"][0].kwargs["background"] == 100


# ---------------------------------------------------------------------------------------------
# DeconWrapperTest in intent: the in-place decon
# ---------------------------------------------------------------------------------------------


def test_large_file_in_place_matches_regular(tmp_path):
    # zarr in, zarr out (demo_large_scale_processing); one batch
    d, fns, psf = dataset(tmp_path, zarr=True)
    run_frame(fns[0], d + "regular", psf, decon_iter=3, save_zarr=True)
    run_frame(fns[0], d + "large", psf, **{**IN_PLACE, "decon_iter": 3})
    large = read(d + "large/" + SMALL_FSN + ".zarr")
    same_bits(large, read(d + "regular/" + SMALL_FSN + ".zarr"))
    same_bits(read(d + "large/MIPs/" + SMALL_FSN + "_MIP_z.tif"), large.max(axis=0))


def test_large_file_in_place_multiple_batches_match_regular(tmp_path):
    # the demo_large_scale_processing setup: zarr chunks deconvolved batch by batch, with edge
    # erosion
    d, fns, psf = dataset(tmp_path, zarr=True)
    e = 3
    run_frame(fns[0], d + "regular", psf, decon_iter=2, edge_erosion=e, save_zarr=True)
    run_frame(
        fns[0], d + "large", psf, **{**IN_PLACE, "edge_erosion": e, "batch_size": (16, 32, 32)}
    )
    regular = read(d + "regular/" + SMALL_FSN + ".zarr")
    large = read(d + "large/" + SMALL_FSN + ".zarr")
    shell = np.ones(SMALL, bool)
    shell[e:-e, e:-e, e:-e] = False
    assert not np.any(large[shell]), "The image border is not eroded."
    b = 8
    close(large[b:-b, b:-b, b:-b], regular[b:-b, b:-b, b:-b])


def test_large_file_in_place_makes_mips_before_publishing(tmp_path):
    from petakit5d.testing.publish_spy import PublishSpy, list_files

    d, fns, psf = dataset(tmp_path, zarr=True)
    with PublishSpy() as spy:
        run_frame(fns[0], d + "large", psf, **IN_PLACE)
    fn = d + "large/" + SMALL_FSN + ".zarr"
    assert spy.was_moved(fn), "The result was not moved into place."
    published = spy.files_when_moved(fn)
    for lab in ("y", "x", "z"):
        assert f"MIPs/{SMALL_FSN}_MIP_{lab}.tif" in published, lab
    added = {
        f
        for f in set(list_files(d + "large")) - set(published)
        if not f.startswith(SMALL_FSN + ".zarr")
    }
    assert not added, f"Files were written after the result was moved into place: {added}"


def test_large_file_in_place_mips_are_merged_from_the_blocks(tmp_path, no_mip_blocks):
    # the decon tasks save the MIPs of their blocks in their flags, and mip_zarr merges them
    from petakit5d.testing.publish_spy import PublishSpy

    d, fns, psf = dataset(tmp_path, zarr=True)
    with PublishSpy() as spy:
        run_frame(fns[0], d + "large", psf, **{**IN_PLACE, "batch_size": (16, 32, 32)})
    flags = [c for c in spy.calls if "/zarr_flag/" in c.dst and c.dst.endswith(".mat")]
    assert flags, "The block flags were not moved into place."
    verify_large_file_mips(d + "large/" + SMALL_FSN + ".zarr")


def test_large_file_in_place_mips_include_masked_blocks(tmp_path, no_mip_blocks):
    # blocks outside the 2D masks are written as zeros without decon, and their MIPs are saved
    d, fns, psf = dataset(tmp_path, zarr=True)
    nz, ny, nx = SMALL
    xy = np.ones((ny, nx), np.uint8)
    xy[32:] = 0
    masks = [str(tmp_path / f"mask_{k}.zarr") for k in ("xy", "xz", "yz")]
    write_zarr(masks[0], xy, (16, 16), "matlab")
    write_zarr(masks[1], np.ones((nx, nz), np.uint8), (16, 16), "matlab")
    write_zarr(masks[2], np.ones((ny, nz), np.uint8), (16, 16), "matlab")
    run_frame(
        fns[0], d + "large", psf, **{**IN_PLACE, "batch_size": (16, 32, 32)}, mask_fullpaths=masks
    )
    result = d + "large/" + SMALL_FSN + ".zarr"
    out = read(result)
    assert not np.any(out[:, 32:]), "The blocks outside the masks were deconvolved."
    assert np.any(out[:, :32])
    verify_large_file_mips(result)


def test_large_file_in_place_passes_job_settings_to_mips(tmp_path, mip_spy):
    # parse_cluster, master_compute and job_log_dir reach mip_zarr
    d, fns, psf = dataset(tmp_path, zarr=True)
    run_frame(fns[0], d + "large", psf, **IN_PLACE, job_log_dir=str(tmp_path / "logs"))
    assert mip_spy.calls, "mip_zarr was not called."
    k = mip_spy.calls[0][1]
    assert k["parse_cluster"] is False
    assert k["master_compute"] is True
    assert k["job_log_dir"] == str(tmp_path / "logs")
    assert len(k["block_flag_fullpaths"]) == 1


def test_large_file_in_place_gpu_job_mips_keep_master_compute(tmp_path, framework_spy, mip_spy):
    # GPU jobs turn master_compute off for the decon blocks only; the MIPs keep the caller's
    d, fns, psf = dataset(tmp_path, zarr=True)
    framework_spy.done = True
    mip_spy.action = lambda: None
    run_frame(fns[0], d + "large", psf, **IN_PLACE, gpu_job=True)
    call = framework_spy.calls[0]
    assert call["master_compute"] is False and call["gpu_job"] is True
    assert call["cpus_per_task"] == 4
    assert all(t.kwargs["use_gpu"] is True for t in call["tasks"])
    assert mip_spy.calls and mip_spy.calls[0][1]["master_compute"] is True


def test_large_file_in_place_keeps_a_result_another_writer_published(tmp_path, framework_spy):
    # the other result comes while the blocks run (the spy reports them done): this one keeps
    # it and removes its own incomplete one, its flags, and makes no MIPs
    from petakit5d.io import write_zarr as wz

    d, fns, psf = dataset(tmp_path, zarr=True)
    result = d + "large/" + SMALL_FSN + ".zarr"
    other = np.full(SMALL, 7, np.uint16)
    framework_spy.done = True
    framework_spy.action = lambda: wz(other, result, block_size=(16, 16, 16))
    run_frame(fns[0], d + "large", psf, **IN_PLACE)
    same_bits(read(result), other)
    assert not list(Path(d, "large").glob(SMALL_FSN + ".zarr_*")), "a temporary result was left"
    assert not list(Path(d, "large", "zarr_flag").glob(SMALL_FSN + "_*")), "flags were left"
    assert not list(Path(d, "large").glob("MIPs/*.tif")), "MIPs of the incomplete result"
    assert framework_spy.calls[0]["final_out_fullpath"] == result


def test_large_file_in_place_keeps_a_result_published_during_its_mips(tmp_path, mip_spy):
    from petakit5d.io import write_zarr as wz

    d, fns, psf = dataset(tmp_path, zarr=True)
    result = d + "large/" + SMALL_FSN + ".zarr"
    other = np.full(SMALL, 7, np.uint16)
    mip_spy.action = lambda: wz(other, result, block_size=(16, 16, 16))
    run_frame(fns[0], d + "large", psf, **{**IN_PLACE, "batch_size": (16, 32, 32)})
    same_bits(read(result), other)
    assert not list(Path(d, "large").glob(SMALL_FSN + ".zarr_*")), "a temporary result was left"
    assert not [n for n in os.listdir(result) if ".zarr_" in n], "moved into the other's"


@pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="Writers outside Slurm are checked on Linux only."
)
def test_large_file_in_place_waits_for_a_controller_of_another_call(
    tmp_path, framework_spy, monkeypatch
):
    # a controller of another call (another uuid) that may still make the volume: this one
    # waits for it (a process on this node, whose result comes during the wait); no block runs
    import subprocess

    from petakit5d.io import write_zarr as wz
    from petakit5d.jobs import files, writers

    for v in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID"):
        monkeypatch.delenv(v, raising=False)
    writers.reset()
    d, fns, psf = dataset(tmp_path, zarr=True)
    result = d + "large/" + SMALL_FSN + ".zarr"
    proc = subprocess.Popen(["sleep", "120"])
    try:
        flag = f"{d}large/zarr_flag/{SMALL_FSN}_otheruuid/"
        os.makedirs(flag)
        owner = f"0_:{writers.node_name()}:{proc.pid}:{writers.process_start(proc.pid)}"
        with open(flag + "writer.txt", "w", newline="") as f:
            f.write(f"result={result}\nuuid=otheruuid\nowner={owner}\n")
        st = os.stat(flag + "writer.txt")
        os.utime(flag + "writer.txt", ns=(st.st_atime_ns, st.st_mtime_ns - 60 * 10**9))
        other = np.full(SMALL, 7, np.uint16)

        def spy(groups, wait_time):
            wz(other, result, block_size=(16, 16, 16))
            return False

        monkeypatch.setattr(files, "wait_for_files", spy)
        run_frame(fns[0], d + "large", psf, **IN_PLACE)
    finally:
        proc.kill()
        proc.wait()
    assert framework_spy.calls == [], "The volume was deconvolved a second time."
    same_bits(read(result), other)
    assert not list(Path(d, "large").glob(SMALL_FSN + ".zarr_*")), "a temporary result was left"


fake_slurm = slurm_fake.fake_slurm  # (the fixture)


@pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
def test_large_file_in_place_blocks_of_an_ended_controller_are_skipped(
    tmp_path, framework_spy, fake_slurm
):
    # the block jobs of a controller whose job was cancelled ran on, for a result nothing would
    # publish; a block skips its work once its controller's job has ended (its writer record; a
    # fake Slurm that doesn't list the job). The block task is the one the controller made.
    from petakit5d.jobs import writers
    from petakit5d.jobs.tasks import run_record

    d, fns, psf = dataset(tmp_path, zarr=True)
    with pytest.raises(L.DeconLargeError, match="missing"):
        run_frame(fns[0], d + "large", psf, **IN_PLACE)
    call = framework_spy.calls[0]
    flag_fn = call["outputs"][0]
    record_fn = os.path.join(os.path.dirname(flag_fn), "writer.txt")
    assert os.path.isfile(record_fn), "The controller did not record itself."
    txt = Path(record_fn).read_text()
    assert "\nended=1" in txt, "The controller that failed did not mark its record ended."
    # a controller whose job was cancelled can't say so: its job is gone
    import re

    txt = re.sub(r"ended=[^\n]*\n", "", re.sub(r"owner=[^\n]*", "owner=777_3:othernode:5", txt))
    with open(record_fn, "w", newline="") as f:
        f.write(txt)
    writers.reset()
    assert run_record(call["tasks"][0].record()) is None
    assert not os.path.isfile(flag_fn), "The block of an ended controller was deconvolved."
    writers.reset()


def test_large_file_in_place_blocks_log_to_the_callers_job_log_dir(tmp_path, framework_spy):
    d, fns, psf = dataset(tmp_path, zarr=True)
    logs = str(tmp_path / "caller_logs")
    os.makedirs(logs)
    with pytest.raises(L.DeconLargeError):
        run_frame(fns[0], d + "large", psf, **IN_PLACE, job_log_dir=logs)
    assert framework_spy.calls, "The controller ran no block framework."
    assert framework_spy.calls[0]["job_log_dir"] == logs
    # (two tries, as MATLAB's)
    assert len(framework_spy.calls) == 2


@pytest.mark.parametrize(
    "gpu_job, master_compute, config, cpus, want",
    [
        (False, True, "", 4, 4),
        (False, True, "off", 4, 1),
        (False, False, "", 4, 1),
        (True, True, "", 2, 1),
    ],
)
def test_large_file_in_place_controllers_ask_for_what_they_hold(
    tmp_path, gpu_job, master_compute, config, cpus, want
):
    # the job of a volume runs its block jobs: it asks for the process and its MIP merge with
    # one CPU, and for the CPUs of a block's job if it deconvolves blocks itself (a 34 G-voxel
    # volume; MATLAB's [4096 4096 2048] with 1024^3 batches and 256^3 chunks)
    from petakit5d.jobs.coordinator import MERGE_MEMORY_FACTOR, PROCESS_MEMORY
    from petakit5d.tools.mip import merge_memory

    cfg = ""
    if config:
        cfg = str(tmp_path / "cpu_config.json")
        Path(cfg).write_text('{"masterCompute": false}')
    size, batch, block = (2048, 4096, 4096), (1024,) * 3, (256,) * 3
    c, mem = L.in_place_controller_resources(
        size,
        batch,
        block,
        gpu_job=gpu_job,
        master_compute=master_compute,
        config_file=cfg,
        cpus_per_task=cpus,
    )
    assert c == want
    merge = merge_memory(size, batch, block, "uint16")
    assert mem == max(PROCESS_MEMORY + merge, MERGE_MEMORY_FACTOR * merge)
    assert mem < 4


# ---------------------------------------------------------------------------------------------
# the batches, the tasks and their values
# ---------------------------------------------------------------------------------------------


def test_border_size_is_half_the_psf_and_covers_the_erosion():
    # round((size + 10) / 2), half away from zero
    assert L.border_size((31, 21, 22)) == (21, 16, 16)
    assert L.border_size((5, 7, 9), 12) == (13, 13, 13)
    assert L.border_size((51, 7, 9), 12) == (31, 13, 13)
    with pytest.raises(ValueError):
        L.border_size((7, 9))


def test_in_place_batches_are_sorted_by_size_then_by_matlab_box():
    from petakit5d.core.blocks import zarr_chunk_coordinates

    # equal sizes: MATLAB's y-fastest order becomes y start first, then x start
    b, r = L.in_place_batches((8, 32, 32), (8, 16, 16), (8, 8, 8), (0, 0, 0))
    b0, r0, _ = zarr_chunk_coordinates(
        (8, 32, 32), batch_size=(8, 16, 16), block_size=(8, 8, 8), border_size=0
    )
    same_bits(b, b0[[0, 2, 1, 3]])
    same_bits(r, r0[[0, 2, 1, 3]])
    # sizes first: the same batch size, cut by the border at the volume's ends (y: 20, 24, 20
    # voxels); the batches of equal size keep their order (MATLAB's sortrows is stable)
    b, r = L.in_place_batches((8, 40, 8), (8, 16, 8), (8, 8, 8), (0, 4, 0))
    assert b[:, [1, 4]].tolist() == [[0, 20], [20, 40], [12, 36]]
    assert r[:, [1, 4]].tolist() == [[0, 16], [32, 40], [16, 32]]
    # y start before z start
    b, _ = L.in_place_batches((16, 32, 8), (8, 16, 8), (8, 8, 8), (0, 0, 0))
    assert b[:, [0, 1]].tolist() == [[0, 0], [8, 0], [0, 16], [8, 16]]


def test_in_place_batches_keep_their_regions():
    from petakit5d.core.blocks import zarr_chunk_coordinates

    size, batch, block, border = (37, 61, 55), (24, 32, 32), (8, 16, 16), (5, 6, 7)
    b, r = L.in_place_batches(size, batch, block, border)
    b0, r0, _ = zarr_chunk_coordinates(
        size, batch_size=batch, block_size=block, same_batch_size=True, border_size=border
    )
    assert b.shape == b0.shape == (8, 6)
    pairs0 = sorted(zip(map(tuple, b0.tolist()), map(tuple, r0.tolist()), strict=True))
    pairs = sorted(zip(map(tuple, b.tolist()), map(tuple, r.tolist()), strict=True))
    assert pairs == pairs0
    vols = [int(np.prod(row[3:] - row[:3])) for row in b]
    assert vols == sorted(vols)


def test_in_memory_batches_are_matlab_order_without_same_size():
    b, r = L.in_memory_batches(IM, (20, 32, 30), (10, 10, 10))
    assert b.shape == (8, 6)
    # y fastest, then x, then z
    starts = [tuple(row[[0, 2, 1]]) for row in r]
    assert starts == sorted(starts)
    # the last region along y ends at the image
    assert r[1][4] == IM[1] and r[1][1] == 32


def test_in_place_task_sizes():
    assert L.in_place_task_size(8, False, False, "simplified") == 5
    assert L.in_place_task_size(8, True, False, "simplified") == 20
    assert L.in_place_task_size(1001, False, False, "simplified") == 5
    assert L.in_place_task_size(1001, False, True, "simplified") == 100
    assert L.in_place_task_size(1001, False, False, "omw") == 100
    assert L.in_place_task_size(1001, True, False, "simplified") == 100
    assert L.in_place_task_size(750001, False, True, "simplified") == 150


def test_block_task_values_keep_matlabs_text_rounding():
    v = L.block_task_values(
        pixel_size=0.108,
        dz=0.2,
        skew_angle=-32.45,
        background=100,
        dz_psf=0.2,
        decon_iter=2,
        wiener_alpha=1e-7 / 3,
        otf_cum_thresh=0.9,
        hann_win_bounds=(0.8, 1.0000000000000002),
        damp_factor=1.2345678901,
        scale_factor=2.5,
        decon_offset=3,
        edge_erosion=0,
    )
    # %0.20d keeps doubles; %d of a fraction is %e with 6 decimals; %0.20f 20 decimals
    assert v["xy_pixel_size"] == 0.108 and v["dz"] == 0.2 and v["skew_angle"] == -32.45
    w = L.block_task_values(**{**_values_kw(), "pixel_size": 0.1083333333333, "dz": 1 / 3})
    assert w["xy_pixel_size"] == 0.1083333333333 and w["dz"] == 1 / 3
    assert v["damp_factor"] == 1.234568
    assert v["scale_factor"] == 2.5
    assert v["wiener_alpha"] == float("3.333333333333e-08")
    assert v["hann_win_bounds"] == [0.8, 1]
    assert v["decon_iter"] == 2 and isinstance(v["decon_iter"], int)
    assert L.block_task_values(**{**_values_kw(), "scale_factor": None})["scale_factor"] is None


def _values_kw():
    return {
        "pixel_size": 0.108,
        "dz": 0.2,
        "skew_angle": 32.45,
        "background": 100,
        "dz_psf": 0.2,
        "decon_iter": 2,
        "wiener_alpha": 0.005,
        "otf_cum_thresh": 0.9,
        "hann_win_bounds": (0.8, 1),
        "damp_factor": 1,
        "scale_factor": 1,
        "decon_offset": 0,
        "edge_erosion": 0,
    }


def test_mask_boxes_are_matlabs_projections():
    region = (2, 3, 4, 12, 13, 14)
    assert L.mask_boxes(region) == [
        (0, 3, 4, 1, 13, 14),
        (0, 4, 2, 1, 14, 12),
        (0, 3, 2, 1, 13, 12),
    ]


def test_in_place_tasks_take_five_batches_with_matlabs_flag_names(tmp_path, framework_spy, mip_spy):
    d, fns, psf = dataset(tmp_path, zarr=True)
    framework_spy.done = True
    mip_spy.action = lambda: None
    # 2 x 3 x 3 = 18 batches
    run_frame(fns[0], d + "large", psf, **{**IN_PLACE, "batch_size": (16, 16, 16)})
    call = framework_spy.calls[0]
    names = [os.path.basename(o) for o in call["outputs"]]
    assert names == ["blocks_1_5.mat", "blocks_6_10.mat", "blocks_11_15.mat", "blocks_16_18.mat"]
    assert [t.args[0] for t in call["tasks"]][-1] == [15, 16, 17]
    t0 = call["tasks"][0]
    assert t0.func == "petakit5d.decon.block:decon_zarr_block"
    assert t0.kwargs["mip_batch_size"] == [32, 48, 48]
    assert t0.kwargs["fix_iter"] is False and t0.kwargs["background"] == 100


def test_in_place_skips_an_existing_result(tmp_path, framework_spy):
    d, fns, psf = dataset(tmp_path, zarr=True)
    os.makedirs(d + "large/" + SMALL_FSN + ".zarr")
    L.decon_large_in_place(fns[0], XY, DZ, d + "large", psf, parse_cluster=False)
    assert framework_spy.calls == []
    assert not os.path.exists(d + "large/zarr_flag")


# ---------------------------------------------------------------------------------------------
# the block task (DeconFunctionTest/gpuBlockTaskWithoutUsableGPUErrors and more)
# ---------------------------------------------------------------------------------------------


def _block_setup(tmp_path):
    from petakit5d.io import create_zarr

    d, fns, psf = dataset(tmp_path, zarr=True)
    out = str(tmp_path / "out.zarr")
    create_zarr(out, data_size=SMALL, block_size=(16, 16, 16), dtype="uint16")
    os.makedirs(tmp_path / "flags")
    flag = str(tmp_path / "flags" / "blocks_1_1.mat")
    box = (0, 0, 0, *SMALL)
    return fns[0], psf, out, flag, box


def _no_usable_gpu() -> bool:
    try:
        from petakit5d.gpu.device import usable_gpu_count

        return usable_gpu_count()[0] == 0
    except Exception:  # noqa: BLE001 (no GPU support at all)
        return True


@pytest.mark.skipif(not _no_usable_gpu(), reason="needs a node without a usable GPU")
def test_gpu_block_task_without_usable_gpu_errors(tmp_path):
    # a large-file GPU job's block fails (so that the framework retries it and the call
    # reports it) instead of running on the CPU; with use_gpu False it runs
    from petakit5d.decon.block import decon_zarr_block

    fn, psf, out, flag, box = _block_setup(tmp_path)
    args = ([0], fn, psf, out, flag, [box], [box], XY, DZ)
    kw = {"background": BG, "decon_iter": 2, "psf_gen": False, "rl_method": "simplified"}
    rl.clear_otf_cache()
    with pytest.raises(Exception) as e:
        decon_zarr_block(*args, **kw, use_gpu=True)
    assert getattr(e.value, "identifier", "") == "decon:noUsableGPU"
    assert not os.path.isfile(flag)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        decon_zarr_block(*args, **kw, use_gpu=False)
    assert os.path.isfile(flag)


def test_block_skips_an_existing_flag_and_overwrite_redoes_it(tmp_path):
    from petakit5d.decon.block import decon_zarr_block
    from petakit5d.io import read_mat

    fn, psf, out, flag, box = _block_setup(tmp_path)
    args = ([0], fn, psf, out, flag, [box], [box], XY, DZ)
    kw = {"background": BG, "decon_iter": 1, "psf_gen": False}
    Path(flag).write_bytes(b"")
    assert decon_zarr_block(*args, **kw) is None
    assert not np.any(read(out))
    assert decon_zarr_block(*args, **kw, overwrite=True) == [True]
    assert np.any(read(out))
    assert sorted(read_mat(flag)) == ["t"]


def test_block_needs_its_flag_folder_and_zarrs(tmp_path):
    from petakit5d.decon.block import decon_zarr_block

    fn, psf, out, flag, box = _block_setup(tmp_path)
    kw = {"background": BG, "decon_iter": 1}
    with pytest.raises(L.DeconLargeError, match="block directory"):
        decon_zarr_block(
            [0], fn, psf, out, str(tmp_path / "nope" / "b.mat"), [box], [box], XY, DZ, **kw
        )
    with pytest.raises(L.DeconLargeError, match="input zarr"):
        decon_zarr_block([0], fn + "x", psf, out, flag, [box], [box], XY, DZ, **kw)
    with pytest.raises(L.DeconLargeError, match="output zarr"):
        decon_zarr_block([0], fn, psf, out + "x", flag, [box], [box], XY, DZ, **kw)


def test_block_saves_the_mips_of_its_blocks(tmp_path):
    from petakit5d.decon.block import decon_zarr_block
    from petakit5d.io import read_mat

    fn, psf, out, flag, box = _block_setup(tmp_path)
    boxes = [(0, 0, 0, 16, 48, 48), (16, 0, 0, 32, 48, 48)]
    decon_zarr_block(
        [0, 1],
        fn,
        psf,
        out,
        flag,
        boxes,
        boxes,
        XY,
        DZ,
        background=BG,
        decon_iter=1,
        psf_gen=False,
        mip_batch_size=(16, 32, 32),
    )
    v = read_mat(flag)
    assert sorted(v) == ["blockBboxes", "mipBatchSize", "t", "xMIPs", "yMIPs", "zMIPs"]
    np.testing.assert_array_equal(v["mipBatchSize"], [[32, 32, 16]])
    np.testing.assert_array_equal(v["blockBboxes"], [[1, 1, 1, 48, 48, 16], [1, 1, 17, 48, 48, 32]])
    assert np.asarray(v["zMIPs"], dtype=object).size == 2


def test_block_stops_every_fifth_batch_once_its_flag_exists(tmp_path, monkeypatch):
    from petakit5d.decon import block as B
    from petakit5d.decon import large

    fn, psf, out, flag, box = _block_setup(tmp_path)
    boxes = [(k * 4, 0, 0, k * 4 + 4, 48, 48) for k in range(7)]
    n = []

    def fake(rawdata, *a, **k):
        n.append(1)
        if len(n) == 5:
            Path(flag).write_bytes(b"")
        return np.zeros((4, 48, 48), np.uint16)

    monkeypatch.setattr(large, "_rl_batch", fake)
    done = B.decon_zarr_block(
        list(range(7)), fn, psf, out, flag, boxes, boxes, XY, DZ, background=BG
    )
    assert len(n) == 5 and done == [True] * 5 + [False] * 2
    # (a GPU task doesn't look)
    os.remove(flag)
    n.clear()
    monkeypatch.setattr(large, "_rl_batch", fake)
    done = B.decon_zarr_block(
        list(range(7)), fn, psf, out, flag, boxes, boxes, XY, DZ, background=BG, use_gpu=True
    )
    assert len(n) == 7 and all(done)
