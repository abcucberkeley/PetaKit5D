"""Deconvolution of one file or an array (petakit5d.decon.file.rl_decon_file) against MATLAB's
RLdecon (the oracle tree's, on the CPU).

Oracle cases (inputs in MATLAB's form; the Python side calls ``petakit5d.matlab.RLdecon``, the
MATLAB-name layer, so the registry's conversions are part of each case):

- file.*: RLdecon of a frame of blurred beads (DeconFrameTest's frame, 48 x 44 x 32, a
  background of 100) and its PSF in a folder; the files written are compared (names, classes,
  voxels, info files): the deskew and rotation after the decon (separate with linear and cubic
  interpolation, combined, with Reverse and Resample, objective scan, rotation without deskew,
  deskew alone), PSF generation (written, read back, a resampled PSF, the raw PSF written as the
  psfgen file without psfGen), zarr input and output, single output, all MIP axes, deconBbox,
  EdgeErosion and an eroded-mask file, dampFactor, scaleFactor and deconOffset, debug output,
  flipZstack, no iterations, all-zero data, the skip rules, save3Dstack, OMW, a PSF larger than
  the data, no output file;
- mem.*: RLdecon of data in memory (DeconFunctionTest's and DeconOverlapTest's calls): the
  result, the overlap function's output; mem.files.*: the files of data in memory, named after
  the output file (its MIP, the deskew and rotation), and none without an output file (the PSF
  and the OMW back projector made without their files);
- error.*: MATLAB's errors ('cudagen', a method name RLdecon doesn't have, the deskew and
  rotation of data in memory without an output file, ...).

Each case runs on both CPU FFT backends: ``.exact`` (MATLAB's bits; skipped without the add-on)
and ``.permissive`` (decision record 0006: uint16 voxels at most 1 count off in at most 0.1% of
them; single outputs relative L2 <= 1e-4; R02's tolerance for back projector files).

The pytest functions below port DeconFrameTest's RLdecon methods and DeconFunctionTest's
in-memory RLdecon methods in intent.
"""

from __future__ import annotations

import math
import os
import warnings

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.decon import file as F
from petakit5d.decon import rl
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, InPath, Tree, Value
from petakit5d.testing.compare import tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, beads, decon_psf, write_tiff, write_zarr_matlab

BACKENDS = ("exact", "permissive")
XY, DZ, BG = 0.108, 0.2, 100
IM = (32, 48, 44)  # (z, y, x): MATLAB's 48 x 44 x 32
BLUR = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)
FSN = (
    "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_"
    "000x_000y_000z_0000t"
)

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
FIGURE = r".*_figure\.png$"


# ---------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------


def frame(seed: int, shape=IM) -> np.ndarray:
    """DeconFrameTest.writeFrame: six blurred beads on a background of 100 (uint16, (z, y, x))"""
    rng = CounterRNG(seed, "frame")
    nz, ny, nx = shape
    u = rng.uniform((6, 3))
    pos = np.column_stack(
        [2 + (nz - 4) * u[:, 0], 6 + (ny - 12) * u[:, 1], 6 + (nx - 12) * u[:, 2]]
    )
    im = beads(shape, pos, BLUR, 1000.0, float(BG))
    return np.floor(im + np.float32(0.5)).astype(np.uint16)


def psf_raw() -> np.ndarray:
    """pktest.deconPSF([31, 31, 31], [1.5, 2]) in uint16: zero in its first and last 5 slices, so
    psf_gen_new keeps it as it is"""
    return np.floor(decon_psf((31, 31, 31), (2.0, 1.5)) + np.float32(0.5)).astype(np.uint16)


def psf_noisy() -> np.ndarray:
    """the PSF with 5% deterministic noise in its nonzero slices: no exactly tied |OTF| values
    (a symmetric PSF's ties at the OMW mask's threshold fall either way with the permissive FFT
    backend's last bits)"""
    p = decon_psf((31, 31, 31), (2.0, 1.5))
    u = CounterRNG(21, "psf_noise").uniform(p.shape, np.float32)
    p = p * (np.float32(1) + np.float32(0.05) * u)
    return np.floor(p + np.float32(0.5)).astype(np.uint16)


def zero_border(a: np.ndarray) -> np.ndarray:
    """zero triangles along x-z, as skewed data has (an erosion mask that isn't full)"""
    nz, _, nx = a.shape
    z = np.arange(nz, dtype=np.float64).reshape(-1, 1, 1)
    x = np.arange(nx, dtype=np.float64).reshape(1, 1, -1)
    keep = (x > 0.5 * z + 2) & (x < 0.5 * z + nx - 14)
    return np.where(keep, a, np.uint16(0))


def make_root(
    *,
    data=None,
    shape=IM,
    zarr_input=False,
    psf=True,
    noisy=False,
    psfgen=None,
    existing=(),
    mask=False,
):
    """a folder with data/<FSN>.tif (or .zarr in MATLAB's layout) and psf/PSF_CamA.tif; psfgen:
    a PSF to put in out/psfgen/PSF_CamA_simplified.tif; existing: result files to put there
    first (relative paths: frames of other beads); mask: an eroded-mask zarr in MATLAB's layout"""

    def make(folder):
        folder = os.fspath(folder)
        os.makedirs(os.path.join(folder, "data"))
        im = frame(1, shape) if data is None else data()
        if zarr_input:
            write_zarr_matlab(
                os.path.join(folder, "data", FSN + ".zarr"), to_matlab(im), (16, 16, 16)
            )
        else:
            write_tiff(os.path.join(folder, "data", FSN + ".tif"), im)
        if psf:
            os.makedirs(os.path.join(folder, "psf"))
            write_tiff(
                os.path.join(folder, "psf", "PSF_CamA.tif"), psf_noisy() if noisy else psf_raw()
            )
        if psfgen is not None:
            os.makedirs(os.path.join(folder, "out", "psfgen"))
            write_tiff(os.path.join(folder, "out", "psfgen", "PSF_CamA_simplified.tif"), psfgen())
        for rel in existing:
            p = os.path.join(folder, *rel.split("/"))
            os.makedirs(os.path.dirname(p), exist_ok=True)
            write_tiff(p, frame(5, shape))
        if mask:
            # all but a corner block: the mask, not an erosion of the data's edges, shapes the
            # result
            m = np.ones(shape, np.uint8)
            m[: shape[0] // 3, : shape[1] // 3, : shape[2] // 3] = 0
            write_zarr_matlab(os.path.join(folder, "mask_eroded.zarr"), to_matlab(m), (16, 16, 16))

    return make


# ---------------------------------------------------------------------------------------------
# file cases
# ---------------------------------------------------------------------------------------------

_FILE_CODE = (
    "d = in{1}; fn = [d, 'data/', in{2}]; o = in{3}; if ~isempty(o), o = [d, o]; end; "
    "p = in{4}; if ~isempty(p), p = [d, p]; end; a = in(8 : end); "
    "i = find(strcmp(a, 'ErodeMaskfile')); if ~isempty(i), a{i + 1} = [d, a{i + 1}]; end; "
    "RLdecon(fn, o, p, in{5}, in{6}, in{7}, a{:}); out = {};"
)


def p_file(d, name, out, psf, xy, dz, dz_psf, *opts, backend):
    from petakit5d import matlab

    a = list(opts)
    if "ErodeMaskfile" in a:
        i = a.index("ErodeMaskfile")
        a[i + 1] = d + a[i + 1]
    with P.using(backend):
        matlab.RLdecon(
            d + "data/" + name, d + out if out else "", d + psf if psf else "", xy, dz, dz_psf, *a
        )
    return []


BASE = [
    "Background",
    BG,
    "DeconIter",
    2,
    "RLMethod",
    "simplified",
    "psfGen",
    False,
    "useGPU",
    False,
]


def file_tree(single: bool = False):
    per_file = {
        "*_back_projector_*": PERMISSIVE_BP,
        "*_OTF.tif": PERMISSIVE_FFT,
        "*/psfgen/*": PERMISSIVE_FFT,
        "*_info.mat": PERMISSIVE_TOL,
        "*": PERMISSIVE_TOL if single else PERMISSIVE_U16,
    }
    return tree(ignore=[FIGURE], per_file=per_file, tiff_bytes=False)


def file_cases(
    cid,
    opts,
    *,
    root=None,
    name=FSN + ".tif",
    out="out/" + FSN + ".tif",
    psf="psf/PSF_CamA.tif",
    dz=DZ,
    base=BASE,
    single=False,
    expect_error=False,
    backends=BACKENDS,
):
    out_list = []
    args = [InPath("root"), name, out, psf, XY, dz, DZ, *base, *opts]
    for be in backends:
        out_list.append(
            Case(
                f"{cid}.{be}",
                code=_FILE_CODE,
                args=args,
                nargout=0,
                python=p_file,
                py_kwargs={"backend": be},
                inputs={"root": Tree(root or make_root(), "root")},
                expect_error=expect_error,
                py_error=Exception,
                **(
                    {"tree": file_tree(single)}
                    if be == "permissive"
                    else {"tree": tree(ignore=[FIGURE])}
                ),
            )
        )
    return out_list


CASES: list[Case] = []
# The rotations after the decon run on the exact backend only: rotateFrame3D's Crop sizes the
# output from where the data is exactly zero, which the permissive backend's last bits can move
ROTATED = ("exact",)
# DeconFrameTest/deskewAndRotationAfterTheDeconRunSeparately
for _interp in ("linear", "cubic"):
    CASES += file_cases(
        f"file.separate.{_interp}",
        ["Deskew", True, "Rotate", True, "DSRCombined", False, "interpMethod", _interp],
        backends=ROTATED,
    )
CASES += file_cases(
    "file.separate.reverse_crop",
    ["Deskew", True, "Rotate", True, "DSRCombined", False, "Reverse", True, "Crop", True],
    backends=ROTATED,
)
CASES += file_cases("file.deskew_separate", ["Deskew", True, "DSRCombined", False])
CASES += file_cases("file.deskew_only", ["Deskew", True])  # MATLAB writes no DS file
# DeconFrameTest/objectiveScanFramesAreRotatedAfterTheDecon
CASES += file_cases(
    "file.objective_scan", ["objectiveScan", True, "Rotate", True], backends=ROTATED
)
CASES += file_cases(
    "file.objective_scan_cubic",
    ["objectiveScan", True, "Rotate", True, "interpMethod", "cubic", "Reverse", True],
    backends=ROTATED,
)
# DeconFrameTest/rotationWithoutDeskewIsTheCombinedDeskewAndRotation (through RLdecon)
CASES += file_cases("file.rotate_only", ["Rotate", True])
CASES += file_cases("file.combined", ["Deskew", True, "Rotate", True, "Reverse", True])
CASES += file_cases(
    "file.combined_resample", ["Deskew", True, "Rotate", True, "Resample", [1, 1, 1.5]]
)
CASES += file_cases(
    "file.combined_single", ["Deskew", True, "Rotate", True, "save16bit", False], single=True
)
CASES += file_cases(
    "file.combined_zarr",
    ["Deskew", True, "Rotate", True, "saveZarr", True],
    out="out/" + FSN + ".zarr",
)
CASES += file_cases("file.save3d_dsr_only", ["Rotate", True, "save3Dstack", [False, True, True]])
CASES += file_cases("file.save3d_scalar_false", ["save3Dstack", False])
# PSF handling
_PSFGEN = ["Background", BG, "DeconIter", 2, "RLMethod", "simplified", "useGPU", False]
CASES += file_cases("file.psfgen", [], base=_PSFGEN)
CASES += file_cases(
    "file.psfgen_cached",
    [],
    base=_PSFGEN,
    root=make_root(psfgen=lambda: (psf_raw()[2:-2, 3:-3, 1:-1] // 2).astype(np.float32)),
)
CASES += file_cases("file.psfgen_false_writes_raw", [])  # (BASE has psfGen false)
CASES += file_cases("file.psfgen_resampled", [], base=_PSFGEN, dz=0.5)
CASES += file_cases(
    "file.psf_larger",
    ["Deskew", True, "Rotate", True],
    base=_PSFGEN,
    root=make_root(shape=(20, 26, 24)),
)
# OMW: the symmetric PSF only with the exact backend (see psf_noisy), the noisy one on both
_OMW = ["Background", BG, "DeconIter", 2, "RLMethod", "omw", "useGPU", False]
CASES += file_cases("file.omw", ["skewed", False], base=_OMW, backends=("exact",))
CASES += file_cases(
    "file.omw_psfgen_false", ["skewed", False, "psfGen", False], base=_OMW, backends=("exact",)
)
CASES += file_cases("file.omw_noisy", ["skewed", False], base=_OMW, root=make_root(noisy=True))
CASES += file_cases(
    "file.omw_noisy_psfgen_false",
    ["skewed", False, "psfGen", False],
    base=_OMW,
    root=make_root(noisy=True),
)
# data and decon options
CASES += file_cases("file.single", ["save16bit", False], single=True)
CASES += file_cases("file.zarr", ["saveZarr", True], out="out/" + FSN + ".zarr")
CASES += file_cases(
    "file.zarr_input",
    ["blockSize", [16, 32, 8], "saveZarr", True],
    root=make_root(zarr_input=True),
    name=FSN + ".zarr",
    out="out/" + FSN + ".zarr",
)
CASES += file_cases("file.mip_all", ["mipAxis", [1, 1, 1]])
CASES += file_cases("file.mip_none", ["mipAxis", [0, 0, 0]])
CASES += file_cases("file.bbox", ["deconBbox", [3, 5, 2, 40, 38, 30]])
CASES += file_cases(
    "file.erosion", ["EdgeErosion", 3], root=make_root(data=lambda: zero_border(frame(1)))
)
CASES += file_cases(
    "file.erode_maskfile",
    ["EdgeErosion", 3, "ErodeMaskfile", "mask_eroded.zarr"],
    root=make_root(mask=True, data=lambda: zero_border(frame(1))),
)
CASES += file_cases(
    "file.damp_scale_offset", ["dampFactor", 1.2, "scaleFactor", 1.5, "deconOffset", 3]
)
CASES += file_cases("file.debug", ["debug", True, "saveStep", 1])
CASES += file_cases("file.flip", ["flipZstack", True])
CASES += file_cases(
    "file.zero_iter",
    ["Deskew", True, "Rotate", True],
    base=[
        "Background",
        BG,
        "DeconIter",
        0,
        "RLMethod",
        "simplified",
        "psfGen",
        False,
        "useGPU",
        False,
    ],
)
CASES += file_cases("file.zero_data", [], root=make_root(data=lambda: np.zeros(IM, np.uint16)))
CASES += file_cases(
    "file.original_zero_iter",
    [],
    base=[
        "Background",
        BG,
        "DeconIter",
        0,
        "RLMethod",
        "original",
        "psfGen",
        False,
        "useGPU",
        False,
    ],
    single=True,
)
# the skip rules
_OUT = "out/" + FSN + ".tif"
CASES += file_cases("file.exists", [], root=make_root(existing=(_OUT,)))
CASES += file_cases(
    "file.exists_rotate",
    ["Rotate", True],
    root=make_root(existing=(_OUT,)),
)
CASES += file_cases(
    "file.exists_all",
    ["Rotate", True],
    root=make_root(existing=(_OUT, "out/DSR/" + FSN + ".tif")),
)
CASES += file_cases(
    "file.exists_dsr",
    ["Rotate", True, "save3Dstack", [False, True, True]],
    root=make_root(existing=("out/DSR/" + FSN + ".tif",)),
)
CASES += file_cases("file.no_output", [], out="")
# errors
CASES += file_cases(
    "error.cudagen",
    [],
    base=[
        "Background",
        BG,
        "DeconIter",
        2,
        "RLMethod",
        "cudagen",
        "psfGen",
        False,
        "useGPU",
        False,
    ],
    expect_error=True,
)
CASES += file_cases(
    "error.unknown_method",
    [],
    base=[
        "Background",
        BG,
        "DeconIter",
        2,
        "RLMethod",
        "Simplified",
        "psfGen",
        False,
        "useGPU",
        False,
    ],
    expect_error=True,
)
CASES += file_cases(
    "error.original_no_background",
    [],
    base=["DeconIter", 0, "RLMethod", "original", "psfGen", False, "useGPU", False],
    expect_error=True,
)


# ---------------------------------------------------------------------------------------------
# in-memory cases
# ---------------------------------------------------------------------------------------------

_MEM_CODE = (
    "d = in{1}; o = in{2}; if ~isempty(o), o = [d, o]; end; f = in{3}; "
    "if f, a = [in(5 : end), {'overlapFcn', @() 'next'}]; else, a = in(5 : end); end; "
    "[r, s, t, v] = RLdecon('', o, [d, 'psf/PSF_CamA.tif'], 0.108, 0.2, 0.2, 'rawdata', in{4}, "
    "a{:}); "
    "out = {r, v};"
)


def p_mem(d, o, overlap, im, *opts, backend):
    from petakit5d import matlab

    a = list(opts)
    if overlap:
        a += ["overlapFcn", lambda: "next"]
    with P.using(backend):
        r = matlab.RLdecon(
            "", d + o if o else "", d + "psf/PSF_CamA.tif", 0.108, 0.2, 0.2, "rawdata", im, *a
        )
    return [r[0], r[3]]


_MEM = ["psfGen", False, "useGPU", False, "Background", BG]


def mem_cases(
    cid,
    opts,
    *,
    out="",
    overlap=False,
    data=None,
    single=True,
    expect_error=False,
    noisy=False,
    backends=BACKENDS,
):
    res = []
    make = data or (lambda: to_matlab(frame(2)))
    for be in backends:
        spec = None
        if be == "permissive" and not expect_error:
            spec = {0: PERMISSIVE_TOL if single else PERMISSIVE_U16}
        res.append(
            Case(
                f"{cid}.{be}",
                code=_MEM_CODE,
                args=[InPath("root"), out, overlap, In("im"), *_MEM, *opts],
                nargout=2,
                python=p_mem,
                py_kwargs={"backend": be},
                inputs={"root": Tree(make_root(noisy=noisy), "root"), "im": Value(make)},
                expect_error=expect_error,
                py_error=Exception,
                files=False,
                **({"compare": spec} if spec else {}),
            )
        )
    return res


# DeconFunctionTest/rldeconInMemoryMatchesDirectFunction (both input classes)
CASES += mem_cases(
    "mem.simplified.uint16", ["RLMethod", "simplified", "DeconIter", 5, "save16bit", False]
)
CASES += mem_cases(
    "mem.simplified.single",
    ["RLMethod", "simplified", "DeconIter", 5, "save16bit", False],
    data=lambda: to_matlab(frame(2).astype(np.float32)),
)
# DeconFunctionTest/rldeconWithZeroIterationsReturnsTheInput
CASES += mem_cases("mem.zero_iter", ["RLMethod", "simplified", "DeconIter", 0], single=False)
# DeconOverlapTest/rldeconReturnsTheOverlapOutput
for _m in ("omw", "simplified"):
    _ov = [
        "RLMethod",
        _m,
        "DeconIter",
        2,
        "skewed",
        False,
        "save3Dstack",
        [False, False, False],
        "mipAxis",
        [0, 0, 0],
    ]
    be = ("exact",) if _m == "omw" else BACKENDS
    CASES += mem_cases(
        f"mem.overlap.{_m}", _ov, out="out/result.tif", overlap=True, single=False, backends=be
    )
    CASES += mem_cases(f"mem.no_overlap.{_m}", _ov, out="out/result.tif", single=False, backends=be)
    if _m == "omw":
        CASES += mem_cases(
            "mem.overlap.omw_noisy",
            _ov,
            out="out/result.tif",
            overlap=True,
            single=False,
            noisy=True,
        )
# data in memory with an output file: the results are named after it (they took the undefined
# input file's name, fsname, until c24801c); without one, OMW makes its back projector in memory
CASES += mem_cases(
    "mem.deskew_rotate",
    ["RLMethod", "simplified", "DeconIter", 2, "Deskew", True, "Rotate", True],
    out="out/result.tif",
    single=False,
)
CASES += mem_cases(
    "mem.output_mip", ["RLMethod", "simplified", "DeconIter", 2], out="out/result.tif", single=False
)
# (the symmetric PSF: exact only, see psf_noisy; mem.files.no_output.omw has the noisy one)
CASES += mem_cases(
    "mem.omw_no_output",
    ["RLMethod", "omw", "DeconIter", 2, "skewed", False],
    single=False,
    backends=("exact",),
)


def mem_file_cases(cid, opts, *, out="", noisy=False, expect_error=False, backends=BACKENDS):
    """RLdecon of data in memory with the files written compared (DeconFunctionTest's
    rldeconInMemoryWithAnOutputFileNamesItsFilesAfterIt and
    rldeconInMemoryWithoutAnOutputFileWritesNothing)"""
    res = []
    for be in backends:
        spec = None
        if be == "permissive" and not expect_error:
            spec = {0: PERMISSIVE_U16}
        res.append(
            Case(
                f"{cid}.{be}",
                code=_MEM_CODE,
                args=[
                    InPath("root"),
                    out,
                    False,
                    In("im"),
                    "useGPU",
                    False,
                    "Background",
                    BG,
                    *opts,
                ],  # fmt: skip
                nargout=2,
                python=p_mem,
                py_kwargs={"backend": be},
                inputs={
                    "root": Tree(make_root(noisy=noisy), "root"),
                    "im": Value(lambda: to_matlab(frame(2))),
                },
                expect_error=expect_error,
                py_error=Exception,
                **({"compare": spec} if spec else {}),
                **(
                    {"tree": file_tree()} if be == "permissive" else {"tree": tree(ignore=[FIGURE])}
                ),
            )
        )
    return res


_SIMPLE = ["RLMethod", "simplified", "DeconIter", 2]
CASES += mem_file_cases("mem.files.output", [*_SIMPLE, "psfGen", False], out="out/result.tif")
CASES += mem_file_cases("mem.files.output_psfgen", [*_SIMPLE, "psfGen", True], out="out/result.tif")
CASES += mem_file_cases(
    "mem.files.deskew_rotate",
    [*_SIMPLE, "psfGen", False, "Deskew", True, "Rotate", True],
    out="out/result.tif",
    backends=ROTATED,
)
CASES += mem_file_cases("mem.files.no_output.simplified", [*_SIMPLE, "psfGen", True])
CASES += mem_file_cases(
    "mem.files.no_output.omw",
    ["RLMethod", "omw", "DeconIter", 2, "psfGen", True, "skewed", False],
    noisy=True,
)
# without an output file there is nowhere to write: MATLAB's undefined outFsn first
CASES += mem_file_cases(
    "error.mem_no_output.deskew_rotate",
    [*_SIMPLE, "psfGen", False, "Deskew", True, "Rotate", True],
    expect_error="MATLAB:UndefinedFunction",
    backends=("permissive",),
)


def _check(case, tmp_path):
    if case.id.endswith(".exact"):
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
# DeconFrameTest and DeconFunctionTest in intent
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def root(tmp_path):
    make_root()(tmp_path)
    return tmp_path


def _norm(psf):
    from petakit5d.kernels import decon as K

    p = np.asarray(psf, np.float32)
    return (p / np.float64(K.sum_double(p))).astype(np.float32)


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_rldecon_in_memory_matches_direct_function(root, cls):
    """the background is subtracted once, also for float input; RLdecon normalizes the PSF (its
    double sum) before passing it on"""
    im = frame(2) if cls == "uint16" else frame(2).astype(np.float32)
    r = F.rl_decon_file(
        "",
        "",
        str(root / "psf" / "PSF_CamA.tif"),
        0.108,
        0.2,
        0.2,
        rawdata=im,
        rl_method="simplified",
        psf_gen=False,
        use_gpu=False,
        background=BG,
        decon_iter=5,
        save_16bit=False,
    )
    rl.clear_otf_cache()
    exp = rl.decon_lucy(im, _norm(psf_raw()), 5, background=BG, use_gpu=False, save_16bit=False)
    assert r.decon.dtype == np.float32
    np.testing.assert_array_equal(r.decon, exp.decon)


def test_rldecon_with_zero_iterations_returns_the_input(root):
    im = frame(2)
    r = F.rl_decon_file(
        "",
        "",
        str(root / "psf" / "PSF_CamA.tif"),
        0.108,
        0.2,
        0.2,
        rawdata=im,
        rl_method="simplified",
        psf_gen=False,
        use_gpu=False,
        background=BG,
        decon_iter=0,
    )
    np.testing.assert_array_equal(r.decon, im)
    assert r.decon.dtype == np.uint16


def test_rldecon_returns_the_overlap_output(root):
    """as the large-file batches call it (data in memory, no result file); simplified here (OMW
    needs the back projector, covered by the mem.overlap.omw cases)"""
    out = str(root / "out" / "result.tif")
    os.makedirs(os.path.dirname(out))
    kw = {
        "rawdata": frame(2),
        "decon_iter": 2,
        "background": BG,
        "use_gpu": False,
        "psf_gen": False,
        "rl_method": "simplified",
        "save_3d_stack": (False, False, False),
        "mip_axis": (0, 0, 0),
    }
    psf = str(root / "psf" / "PSF_CamA.tif")
    a = F.rl_decon_file("", out, psf, 0.108, 0.2, 0.2, **kw, overlap_fcn=lambda: "next")
    b = F.rl_decon_file("", out, psf, 0.108, 0.2, 0.2, **kw)
    assert a.overlap_out == "next" and b.overlap_out is None
    np.testing.assert_array_equal(a.decon, b.decon)
    assert not os.path.exists(out)


def _run(root, name, *, out_name=None, **kw):
    fn = str(root / "data" / (FSN + ".tif"))
    out_dir = root / name
    os.makedirs(out_dir, exist_ok=True)
    out = str(out_dir / (out_name or FSN + ".tif"))
    base = {
        "background": BG,
        "decon_iter": 2,
        "rl_method": "simplified",
        "psf_gen": False,
        "use_gpu": False,
    }
    base.update(kw)
    r = F.rl_decon_file(fn, out, str(root / "psf" / "PSF_CamA.tif"), XY, DZ, DZ, **base)
    return r, str(out_dir)


@pytest.mark.parametrize("interp", ["linear", "cubic"])
def test_deskew_and_rotation_after_the_decon_run_separately(root, interp):
    """DeconFrameTest: the separate deskew and rotation (DSRCombined false) pass the interpolation
    on (M07's 8146dc5)"""
    from petakit5d.deskew_rotate import deskew_volume, rotate_volume
    from petakit5d.io import read_tiff
    from petakit5d.matlabcompat.elementary import sind

    r, d = _run(root, interp, deskew=True, rotate=True, dsr_combined=False, interp_method=interp)
    decon = read_tiff(os.path.join(d, FSN + ".tif"))
    ds = deskew_volume(decon.astype(np.float32), 32.45, DZ, XY, False, crop=False, interp=interp)
    np.testing.assert_array_equal(read_tiff(os.path.join(d, "DS", FSN + ".tif")), _u16(ds))
    z_aniso = float(sind(32.45)) * DZ / XY
    dsr = rotate_volume(ds, 32.45, z_aniso, False, crop=True, interp=interp)
    np.testing.assert_array_equal(read_tiff(os.path.join(d, "DSR", FSN + ".tif")), _u16(dsr))
    assert os.path.isfile(os.path.join(d, "DS", "MIPs", FSN + "_MIP_z.tif"))
    assert os.path.isfile(os.path.join(d, "DSR", "MIPs", FSN + "_MIP_z.tif"))
    assert r.ds is None and r.dsr.dtype == np.uint16


def test_cubic_interpolation_is_used(root):
    a, _ = _run(root, "lin", deskew=True, dsr_combined=False, save_3d_stack=(False, True, False))
    b, _ = _run(
        root,
        "cub",
        deskew=True,
        dsr_combined=False,
        interp_method="cubic",
        save_3d_stack=(False, True, False),
    )
    assert not np.array_equal(a.ds, b.ds)


def _u16(a):
    from petakit5d.kernels import ops

    return ops.convert(a, "uint16")


def test_objective_scan_frames_are_rotated_after_the_decon(root):
    from petakit5d.deskew_rotate import rotate_volume
    from petakit5d.io import read_tiff

    _, d = _run(root, "decon", objective_scan=True, rotate=True)
    decon = read_tiff(os.path.join(d, FSN + ".tif"))
    exp = _u16(
        rotate_volume(
            decon.astype(np.float32), 32.45, DZ / XY, False, crop=True, objective_scan=True
        )
    )
    np.testing.assert_array_equal(read_tiff(os.path.join(d, "Rotated", FSN + ".tif")), exp)
    mip = read_tiff(os.path.join(d, "Rotated", "MIPs", FSN + "_MIP_z.tif"))
    np.testing.assert_array_equal(mip, exp.max(axis=0))


def test_rotation_without_deskew_is_the_combined_deskew_and_rotation(root):
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.io import read_tiff

    _, d = _run(root, "decon", rotate=True, reverse=True)
    decon = read_tiff(os.path.join(d, FSN + ".tif"))
    exp = deskew_rotate_volume(decon, 32.45, DZ, XY, reverse=True, crop=True, save_16bit=True)
    np.testing.assert_array_equal(read_tiff(os.path.join(d, "DSR", FSN + ".tif")), _u16(exp))


def test_deskew_alone_under_dsr_combined_writes_nothing(root):
    """MATLAB's behaviour (M07's report, not fixed): Deskew without Rotate and with the default
    DSRCombined makes no DS result"""
    r, d = _run(root, "decon", deskew=True)
    assert r.ds is None and r.dsr is None
    assert not os.path.exists(os.path.join(d, "DS"))
    assert os.path.isfile(os.path.join(d, FSN + ".tif"))


def test_the_separate_path_ignores_resample(root):
    """MATLAB's behaviour (M07's report, not fixed): Resample reaches only the combined path"""
    a, _ = _run(root, "a", deskew=True, rotate=True, dsr_combined=False)
    b, _ = _run(root, "b", deskew=True, rotate=True, dsr_combined=False, resample=(2, 1, 1))
    np.testing.assert_array_equal(a.dsr, b.dsr)


def _files(folder) -> list[str]:
    """the files under folder, relative, with forward slashes"""
    return sorted(
        os.path.relpath(os.path.join(d, f), folder).replace(os.sep, "/")
        for d, _, fs in os.walk(folder)
        for f in fs
    )


def test_rldecon_in_memory_with_an_output_file_names_its_files_after_it(root):
    """(MATLAB: rldeconInMemoryWithAnOutputFileNamesItsFilesAfterIt) data in memory with an output
    file: its MIP took the input file's name (fsname undefined) and failed; the result and its
    MIP are named after the output file and equal those of the decon of the file"""
    from petakit5d.io import read_tiff, write_tiff

    im = frame(2)
    fn = str(root / "in.tif")
    write_tiff(im, fn)
    psf = str(root / "psf" / "PSF_CamA.tif")
    kw = {"rl_method": "simplified", "psf_gen": False, "use_gpu": False, "background": BG,
          "decon_iter": 2}  # fmt: skip
    F.rl_decon_file(fn, str(root / "file" / "in.tif"), psf, XY, DZ, DZ, **kw)
    os.makedirs(root / "memory")
    F.rl_decon_file("", str(root / "memory" / "out.tif"), psf, XY, DZ, DZ, rawdata=im, **kw)
    np.testing.assert_array_equal(
        read_tiff(str(root / "memory" / "out.tif")), read_tiff(str(root / "file" / "in.tif"))
    )
    np.testing.assert_array_equal(
        read_tiff(str(root / "memory" / "MIPs" / "out_MIP_z.tif")),
        read_tiff(str(root / "file" / "MIPs" / "in_MIP_z.tif")),
    )
    assert _files(root / "memory") == [
        "MIPs/out_MIP_z.tif",
        "decon_info/out_actual_iterations_2.txt",
        "decon_info/out_info.mat",
        "out.tif",
    ]


def test_rldecon_in_memory_deskews_and_rotates_into_the_output_folder(root):
    """data in memory with an output file and Deskew and Rotate: the DSR and its MIP are named
    after the output file, as the decon's MIP"""
    from petakit5d.io import read_tiff

    out = root / "out"
    os.makedirs(out)
    r = F.rl_decon_file(
        "", str(out / "r.tif"), str(root / "psf" / "PSF_CamA.tif"), XY, DZ, DZ, rawdata=frame(2),
        rl_method="simplified", psf_gen=False, use_gpu=False, background=BG, decon_iter=2,
        deskew=True, rotate=True,
    )  # fmt: skip
    np.testing.assert_array_equal(read_tiff(str(out / "DSR" / "r.tif")), r.dsr)
    np.testing.assert_array_equal(
        read_tiff(str(out / "DSR" / "MIPs" / "r_MIP_z.tif")), r.dsr.max(axis=0)
    )
    assert os.path.isfile(out / "MIPs" / "r_MIP_z.tif")


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_rldecon_in_memory_without_an_output_file_writes_nothing(root, tmp_path, method,
                                                                  monkeypatch):  # fmt: skip
    """(MATLAB: rldeconInMemoryWithoutAnOutputFileWritesNothing) without an input or an output
    file there is no folder: OMW (the back projector's files) and psfGen (the generated PSF's
    file) failed on the undefined folder. Both are made in memory, nothing is written, and the
    result equals the decon of the file"""
    from petakit5d.io import read_tiff, write_tiff

    im = frame(2)
    fn = str(root / "in.tif")
    write_tiff(im, fn)
    psf = str(root / "psf" / "PSF_CamA.tif")
    kw = {"rl_method": method, "psf_gen": True, "use_gpu": False, "background": BG,
          "decon_iter": 2}  # fmt: skip
    if method == "omw":
        kw["skewed"] = False
    F.rl_decon_file(fn, str(root / "file" / "in.tif"), psf, XY, DZ, DZ, **kw)
    here = tmp_path / "here"
    here.mkdir()
    monkeypatch.chdir(here)
    rl.clear_otf_cache()
    r = F.rl_decon_file("", "", psf, XY, DZ, DZ, rawdata=im, **kw)
    np.testing.assert_array_equal(r.decon, read_tiff(str(root / "file" / "in.tif")))
    assert _files(here) == [], "Files were written without an output file."
    if method == "omw":
        assert any("_back_projector_" in f for f in _files(root / "file"))


def test_rldecon_in_memory_without_an_output_file_cannot_deskew(tmp_path, root, monkeypatch):
    """MATLAB's behaviour kept: the deskew and rotation of data in memory without an output file
    fail on the undefined output name (outFsn) before anything is written"""
    monkeypatch.chdir(tmp_path)
    before = _files(tmp_path)
    with pytest.raises(rl.DeconError) as e:
        F.rl_decon_file(
            "", "", str(root / "psf" / "PSF_CamA.tif"), XY, DZ, DZ, rawdata=frame(2),
            rl_method="simplified", psf_gen=False, use_gpu=False, background=BG, decon_iter=1,
            deskew=True, rotate=True,
        )  # fmt: skip
    assert e.value.identifier == "MATLAB:UndefinedFunction"
    assert "'outFsn'" in str(e.value)
    assert _files(tmp_path) == before


def test_rldecon_in_memory_debug_folder_is_named_after_the_output(root):
    out = root / "out"
    os.makedirs(out)
    F.rl_decon_file(
        "", str(out / "r.tif"), str(root / "psf" / "PSF_CamA.tif"), XY, DZ, DZ, rawdata=frame(2),
        rl_method="simplified", psf_gen=False, use_gpu=False, background=BG, decon_iter=2,
        debug=True, save_step=1,
    )  # fmt: skip
    assert os.path.isdir(out / "r_debug") and os.listdir(out / "r_debug")


@pytest.mark.parametrize("name", ["Simplified", "OMW", "rl", "", 3])
def test_rldecon_rejects_an_unknown_method(root, name):
    """(MATLAB: rldeconRejectsAnUnknownMethod) RLdecon's method names are case-sensitive: another
    name left the result empty, and the decon failed writing its MIP; it fails with the
    arguments, before anything is written"""
    from petakit5d.params.spec import MatlabArgumentError

    with pytest.raises(MatlabArgumentError) as e:
        F.rl_decon_file(
            str(root / "data" / (FSN + ".tif")), str(root / "out" / "in.tif"),
            str(root / "psf" / "PSF_CamA.tif"), XY, DZ, DZ, rl_method=name, psf_gen=False,
            use_gpu=False, background=BG, decon_iter=2,
        )  # fmt: skip
    assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
    assert "'RLMethod'" in str(e.value)
    assert not os.path.exists(root / "out"), "The decon of an unknown method wrote files."


def test_methods_that_are_not_ported(root):
    psf = str(root / "psf" / "PSF_CamA.tif")
    kw = {
        "rawdata": frame(2),
        "psf_gen": False,
        "use_gpu": False,
        "background": BG,
        "decon_iter": 1,
    }
    with pytest.raises(NotImplementedError, match="B11"):
        F.rl_decon_file("", "", psf, XY, DZ, DZ, rl_method="original", **kw)
    with pytest.raises(rl.DeconError) as e:
        F.rl_decon_file("", "", psf, XY, DZ, DZ, rl_method="cudagen", **kw)
    assert e.value.identifier == "MATLAB:UndefinedFunction"
    with pytest.raises(NotImplementedError, match="edgetaper"):
        F.rl_decon_file("", "", psf, XY, DZ, DZ, rl_method="simplified", n_tapering=3, **kw)


def test_original_without_iterations_subtracts_the_background(root):
    r = F.rl_decon_file(
        "",
        "",
        str(root / "psf" / "PSF_CamA.tif"),
        XY,
        DZ,
        DZ,
        rawdata=frame(2),
        rl_method="original",
        psf_gen=False,
        use_gpu=False,
        background=BG,
        decon_iter=0,
    )
    exp = np.maximum(frame(2).astype(np.float32) - np.float32(BG), np.float32(0))
    np.testing.assert_array_equal(r.decon, exp)


def test_a_mat_psf_is_read(root, tmp_path):
    from petakit5d.io.matfile import write_mat

    p = str(tmp_path / "psf.mat")
    write_mat(p, {"psf": to_matlab(psf_raw().astype(np.float32))})
    kw = {
        "rawdata": frame(2),
        "psf_gen": False,
        "use_gpu": False,
        "background": BG,
        "decon_iter": 2,
        "rl_method": "simplified",
        "save_16bit": False,
    }
    a = F.rl_decon_file("", "", p, XY, DZ, DZ, **kw)
    b = F.rl_decon_file("", "", str(root / "psf" / "PSF_CamA.tif"), XY, DZ, DZ, **kw)
    np.testing.assert_array_equal(a.decon, b.decon)


def test_existing_results_are_skipped_and_kept(root):
    from petakit5d.io import read_tiff

    r, d = _run(root, "decon")
    before = read_tiff(os.path.join(d, FSN + ".tif"))
    r2, _ = _run(root, "decon", decon_iter=5)
    assert r2.decon is None
    np.testing.assert_array_equal(read_tiff(os.path.join(d, FSN + ".tif")), before)


def test_an_existing_tiff_result_is_written_again_and_a_zarr_kept(root, monkeypatch):
    """the decon result exists (DSR to make): it is read; MATLAB moves the result over it again
    (a tiff: the same voxels), a zarr result is kept"""
    from petakit5d.io import read_tiff, read_zarr

    moves = []
    real = F._move
    monkeypatch.setattr(F, "_move", lambda a, b: moves.append(b) or real(a, b))
    _, d = _run(root, "decon")
    before = read_tiff(os.path.join(d, FSN + ".tif"))
    moves.clear()
    _run(root, "decon", rotate=True)
    assert [os.path.basename(os.path.dirname(p)) for p in moves] == ["DSR", "decon"]
    np.testing.assert_array_equal(read_tiff(os.path.join(d, FSN + ".tif")), before)
    assert not [f for f in os.listdir(d) if f.endswith(".tif") and f != FSN + ".tif"]
    _, z = _run(root, "zarr", save_zarr=True, out_name=FSN + ".zarr")
    before = read_zarr(os.path.join(z, FSN + ".zarr"))
    moves.clear()
    _run(root, "zarr", save_zarr=True, out_name=FSN + ".zarr", rotate=True)
    assert [os.path.basename(os.path.dirname(p)) for p in moves] == ["DSR"]
    np.testing.assert_array_equal(read_zarr(os.path.join(z, FSN + ".zarr")), before)
    assert sorted(os.listdir(z)) == ["DSR", "MIPs", FSN + ".zarr", "decon_info", "psfgen"]


def test_an_eroded_mask_file_replaces_the_erosion(tmp_path):
    """EdgeErosion with an existing mask file: no erosion in the decon (data with zero borders,
    whose erosion would show), the result times the mask"""
    from petakit5d.io import read_tiff

    make_root(mask=True, data=lambda: zero_border(frame(1)))(tmp_path)
    root = tmp_path
    mask = str(tmp_path / "mask_eroded.zarr")
    _, d = _run(root, "a", edge_erosion=3, erode_maskfile=mask)
    b, _ = _run(root, "b", edge_erosion=0)
    m = np.ones(IM, bool)
    m[: IM[0] // 3, : IM[1] // 3, : IM[2] // 3] = False
    np.testing.assert_array_equal(read_tiff(os.path.join(d, FSN + ".tif")), np.where(m, b.decon, 0))


def test_rawdata_given_with_a_file_is_used_and_not_changed(root):
    """XR_RLdeconFrame3D's read-ahead hands RLdecon the data of the file: flipped in a copy"""
    im = frame(3)
    keep = im.copy()
    r, _ = _run(root, "a", rawdata=im, flip_zstack=True, save_3d_stack=(False, False, False))
    np.testing.assert_array_equal(im, keep)
    b, _ = _run(root, "b", rawdata=im[::-1].copy(), save_3d_stack=(False, False, False))
    np.testing.assert_array_equal(r.decon, b.decon)


def test_from_matlab_roundtrip_of_the_frame():
    """the inputs are MATLAB's 48 x 44 x 32 (never cubic)"""
    m = to_matlab(frame(1))
    assert m.shape == (48, 44, 32)
    np.testing.assert_array_equal(from_matlab(m), frame(1))
