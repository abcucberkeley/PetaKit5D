"""Deconvolution of datasets (petakit5d.decon.dataset.decon_dataset) against MATLAB's
XR_decon_data_wrapper (the oracle tree's, on the CPU, in local mode).

Oracle cases: XR_decon_data_wrapper on folders of frames of blurred beads (DeconWrapperTest's
light-sheet tiles, 48 x 44 x 32) with one PSF per camera, through
``petakit5d.matlab.XR_decon_data_wrapper``; the files written are compared (results, MIPs,
generated PSFs, back projectors, masks, info files, parameters.mat and parameters.json): RL and
OMW, tiff and zarr, uint16 and single; two channels with their own PSFs, backgrounds and OMW
parameters; the first time point's erosion; two datasets with their own iterations; PSFs
resampled to the data's z step; psfGen off; deconRotate; overwrite and existing results; the
widefield demo's settings; parameters with long decimals (MATLAB's rounding in the function
string); and the error for a missing PSF. Each runs on both CPU FFT backends (``.exact``:
MATLAB's bits; ``.permissive``: decision record 0006).

The pytest functions port DeconWrapperTest's methods in intent: the outputs against the decon
functions, and the cluster mode against the fake Slurm (jobs.single's single jobs, task runners,
waits, GPU jobs of several files, the OMW back projector made once).
"""

from __future__ import annotations

import json
import math
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.decon import dataset as D
from petakit5d.decon import rl
from petakit5d.jobs import config, single
from petakit5d.jobs import files as jfiles
from petakit5d.testing import oracle, slurm_fake
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.compare import tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.slurm_fake import SLURM_JOB_VARIABLES, hide_slurm
from petakit5d.testing.synth import CounterRNG, beads, decon_psf, write_tiff, write_zarr_matlab

fake_slurm = slurm_fake.fake_slurm  # (the fixture)

BACKENDS = ("exact", "permissive")
XY, DZ, BG = 0.108, 0.2, 100
IM = (32, 48, 44)  # (z, y, x): MATLAB's 48 x 44 x 32
BLUR = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)
POSIX = pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")


def fsn(cam: str, t: int, x: int = 0) -> str:
    return (
        f"Scan_Iter_0000_0000_Cam{cam}_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_"
        f"{x:03d}x_000y_000z_{t:04d}t"
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
    reason="0006: FFT-based outputs with the permissive FFT backend, single precision (max <= "
    "1e-5 of the largest value, relative L2 <= 1e-6)",
    rel_max=1e-5,
    rel_l2=1e-6,
)
FIGURE = r".*_figure\.png$"


# ---------------------------------------------------------------------------------------------
# data (DeconWrapperTest.createDataset, deterministic)
# ---------------------------------------------------------------------------------------------


def frame(seed: int, shape=IM) -> np.ndarray:
    """six blurred beads on a background of 100 (uint16, (z, y, x))"""
    rng = CounterRNG(seed, "wrapper_frame")
    nz, ny, nx = shape
    u = rng.uniform((6, 3))
    pos = np.column_stack(
        [2 + (nz - 4) * u[:, 0], 6 + (ny - 12) * u[:, 1], 6 + (nx - 12) * u[:, 2]]
    )
    im = beads(shape, pos, BLUR, 1000.0, float(BG))
    return np.floor(im + np.float32(0.5)).astype(np.uint16)


def psf_raw(sigma=(2.0, 1.5), noisy=False, seed=21) -> np.ndarray:
    """pktest.deconPSF([31 31 31], sigma); noisy: 5% deterministic noise, so that no |OTF|
    values tie at the OMW mask's threshold (a symmetric PSF's ties fall either way with the
    permissive backend's last bits)"""
    p = decon_psf((31, 31, 31), sigma)
    if noisy:
        u = CounterRNG(seed, "psf_noise").uniform(p.shape, np.float32)
        p = p * (np.float32(1) + np.float32(0.05) * u)
    return np.floor(p + np.float32(0.5)).astype(np.uint16)


def make_root(
    *,
    cams=("A",),
    n_t=1,
    sigmas=None,
    zarr=False,
    noisy=False,
    datasets=("data",),
    names=None,
    psf_names=None,
    zero_cols=False,
    existing=(),
):
    """<dataset>/<frames> (.tif, or MATLAB-layout .zarr with 16^3 chunks) and psf/PSF_Cam<c>.tif;
    names: plain file names instead of the tile names; zero_cols: the first time point's first
    10 columns zero; existing: result files there first (relative paths)"""

    def make(folder):
        folder = os.fspath(folder)
        os.makedirs(os.path.join(folder, "psf"))
        for c, cam in enumerate(cams):
            sig = (sigmas or {}).get(cam, (2.0, 1.5))
            pname = (psf_names or {}).get(cam, f"PSF_Cam{cam}")
            write_tiff(os.path.join(folder, "psf", pname + ".tif"), psf_raw(sig, noisy, 21 + c))
        seed = 1
        for ds in datasets:
            os.makedirs(os.path.join(folder, ds))
            for c, cam in enumerate(cams):
                for t in range(n_t):
                    name = names[c] if names else fsn(cam, t)
                    im = frame(seed)
                    seed += 1
                    if zero_cols and t == 0:
                        im[:, :, :10] = 0
                    if zarr:
                        write_zarr_matlab(
                            os.path.join(folder, ds, name + ".zarr"),
                            im.transpose(1, 2, 0),
                            (16, 16, 16),
                        )
                    else:
                        write_tiff(os.path.join(folder, ds, name + ".tif"), im)
        for rel in existing:
            p = os.path.join(folder, *rel.split("/"))
            os.makedirs(os.path.dirname(p), exist_ok=True)
            write_tiff(p, frame(99))

    return make


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

# in{1}: the root; in{2}: the dataset folder(s) under it; then the parameters, with the PSFs
# (and masks) relative to the root
_CODE = (
    "d = in{1}; p = in{2}; if iscell(p), p = cellfun(@(x) [d, x], p, 'UniformOutput', false); "
    "else, p = [d, p]; end; a = in(3 : end); "
    "for i = find(cellfun(@(x) ischar(x) && strcmpi(x, 'psfFullpaths'), a)), "
    "a{i + 1} = cellfun(@(x) [d, x], a{i + 1}, 'UniformOutput', false); end; "
    "XR_decon_data_wrapper(p, a{:}); out = {};"
)


def p_wrapper(d, data, *opts, backend):
    from petakit5d import matlab

    p = [d + x for x in data] if isinstance(data, (list, tuple)) else d + data
    a = list(opts)
    for i, v in enumerate(opts):
        if isinstance(v, str) and v.lower() == "psffullpaths":
            a[i + 1] = [d + x for x in a[i + 1]]
    rl.clear_otf_cache()
    with P.using(backend):
        matlab.XR_decon_data_wrapper(p, *a)
    return []


BASE = [
    "channelPatterns",
    ["CamA_ch0"],
    "psfFullpaths",
    ["psf/PSF_CamA.tif"],
    "xyPixelSize",
    XY,
    "dz",
    DZ,
    "dzPSF",
    DZ,
    "background",
    BG,
    "edgeErosion",
    0,
    "save16bit",
    True,
    "parseCluster",
    False,
    "GPUJob",
    False,
    "deconIter",
    3,
]


def wrapper_tree(single=False):
    return tree(
        ignore=[FIGURE],
        per_file={
            "*_back_projector_*": PERMISSIVE_BP,
            "*_OTF.tif": PERMISSIVE_FFT,
            "*/psfgen/*": PERMISSIVE_FFT,
            "*_info.mat": PERMISSIVE_TOL,
            "*": PERMISSIVE_TOL if single else PERMISSIVE_U16,
        },
        tiff_bytes=False,
    )


def wrapper_cases(
    cid,
    opts,
    *,
    root=None,
    data="data/",
    base=BASE,
    single=False,
    expect_error=False,
    backends=BACKENDS,
    device="cpu",
):
    res = []
    for be in backends:
        res.append(
            Case(
                f"{cid}.{be}",
                device=device,
                code=_CODE,
                args=[InPath("root"), data, *base, *opts],
                nargout=0,
                python=p_wrapper,
                py_kwargs={"backend": be},
                inputs={"root": Tree(root or make_root(), "root")},
                expect_error=expect_error,
                py_error=Exception,
                tree=wrapper_tree(single) if be == "permissive" else tree(ignore=[FIGURE]),
            )
        )
    return res


_OMW = ["RLMethod", "omw", "skewed", False, "deconIter", 2]
_NOISY = make_root(noisy=True)
CASES: list[Case] = []
# simplifiedWritesOutputsThatMatchTheFunction, save16bitFalseWritesSingle,
# zarrOutputMatchesTiffOutput (RL and OMW, tiff and zarr, uint16 and single)
CASES += wrapper_cases("rl.tif.uint16", [], root=make_root(n_t=2))
CASES += wrapper_cases("rl.tif.single", ["save16bit", False], single=True)
CASES += wrapper_cases(
    "rl.zarr.uint16", ["zarrFile", True, "saveZarr", True], root=make_root(zarr=True)
)
CASES += wrapper_cases(
    "rl.zarr.single",
    ["zarrFile", True, "saveZarr", True, "save16bit", False],
    root=make_root(zarr=True),
    single=True,
)
CASES += wrapper_cases("rl.tif_to_zarr", ["saveZarr", True, "blockSize", [16, 32, 8]])
# omwWritesBackProjectorFilesAndMatchesTheFunction (the symmetric PSF on the exact backend only)
CASES += wrapper_cases("omw.tif.uint16", _OMW, backends=("exact",))
CASES += wrapper_cases("omw.tif.uint16_noisy", _OMW, root=_NOISY)
CASES += wrapper_cases("omw.tif.single", [*_OMW, "save16bit", False], root=_NOISY, single=True)
CASES += wrapper_cases(
    "omw.zarr.uint16",
    [*_OMW, "zarrFile", True, "saveZarr", True],
    root=make_root(zarr=True, noisy=True),
)
CASES += wrapper_cases(
    "omw.zarr.single",
    [*_OMW, "zarrFile", True, "saveZarr", True, "save16bit", False],
    root=make_root(zarr=True, noisy=True),
    single=True,
)
# eachChannelUsesItsOwnPSFAndBackground, omwParametersArePerChannel
_TWO = make_root(cams=("A", "B"), sigmas={"B": (3.0, 2.5)}, noisy=True)
_TWO_BASE = [
    "channelPatterns",
    ["CamA_ch0", "CamB_ch0"],
    "psfFullpaths",
    ["psf/PSF_CamA.tif", "psf/PSF_CamB.tif"],
    *BASE[4:],
]
CASES += wrapper_cases("channels.rl", ["background", [100, 80]], root=_TWO, base=_TWO_BASE)
CASES += wrapper_cases(
    "channels.omw",
    [*_OMW, "wienerAlpha", [0.05, 0.005], "OTFCumThresh", [0.88, 0.9], "scaleFactor", [1, 2]],
    root=_TWO,
    base=_TWO_BASE,
)
# erodeByFirstTimePointAppliesItsMaskToAllFiles
CASES += wrapper_cases(
    "erode_by_ftp", ["edgeErosion", 2, "erodeByFTP", True], root=make_root(n_t=2, zero_cols=True)
)
CASES += wrapper_cases(
    "erosion_own", ["edgeErosion", 2, "erodeByFTP", False], root=make_root(n_t=2)
)
# two datasets, iterations per dataset
CASES += wrapper_cases(
    "datasets.iterations",
    ["deconIter", [2, 3]],
    root=make_root(datasets=("d1", "d2")),
    data=["d1/", "d2/"],
)
# the PSF resampled to the data's z step (dzPSF 0.1, dz 0.2), and psfGen off
CASES += wrapper_cases("psf.resampled", ["dzPSF", 0.1])
CASES += wrapper_cases("psf.omw_resampled", [*_OMW, "dzPSF", 0.1], root=_NOISY)
CASES += wrapper_cases("psf.gen_false", ["psfGen", False])
# deconRotateWritesTheDeskewedRotatedResult
CASES += wrapper_cases("decon_rotate", ["deconRotate", True, "deconIter", 2])
# overwrite (the folder removed first) and existing results (kept)
_EXIST = make_root(existing=("data/matlab_decon/" + fsn("A", 0) + ".tif",))
CASES += wrapper_cases("existing_kept", [], root=_EXIST)
CASES += wrapper_cases("overwrite", ["overwrite", True], root=_EXIST)
CASES += wrapper_cases("result_dir", ["resultDirName", "decon_omw", *_OMW], root=_NOISY)
# options passed on: damp factor, scale, offset, debug
CASES += wrapper_cases(
    "options",
    ["dampFactor", 1.2, "scaleFactor", 2, "deconOffset", 5, "debug", True, "saveStep", 2],
)
CASES += wrapper_cases("default_background", ["background", []])
# parameters with long decimals: their values in the files are MATLAB's after its
# function string (%.10f, %d of a non-integer, %.20f, mat2str)
_LONG = [
    "xyPixelSize",
    0.10812345678901234,
    "dz",
    0.20000000001234,
    "dzPSF",
    0.1234567890123456,
    "background",
    99.123456789,
    "dampFactor",
    1.2345678901,
    "scaleFactor",
    1.0123456789,
    "deconOffset",
    0.123456789,
    "skewAngle",
    31.123456789,
    "deconRotate",
    True,
]
CASES += wrapper_cases("long_decimals.rl", _LONG)
CASES += wrapper_cases(
    "long_decimals.omw",
    [
        *_OMW,
        *_LONG[:-2],
        "wienerAlpha",
        0.0051234567890123456789,
        "OTFCumThresh",
        0.8812345678901234567,
        "hannWinBounds",
        [0.81234567890123456, 1.0987654321098765],
    ],
    root=_NOISY,
)
# widefieldDemoSettingsWithPlainFilenames (Reverse and fixIter are not the wrapper's: ignored)
_WF = make_root(cams=("A",), names=["WF_raw"], psf_names={"A": "WF_PSF"}, noisy=True)
CASES += wrapper_cases(
    "widefield_demo",
    [
        "resultDirName",
        "matlab_decon_omw",
        "xyPixelSize",
        0.157,
        "dz",
        0.1,
        "Reverse",
        True,
        "ChannelPatterns",
        ["WF_raw"],
        "PSFFullpaths",
        ["psf/WF_PSF.tif"],
        "dzPSF",
        0.1,
        "RLmethod",
        "omw",
        "wienerAlpha",
        0.02,
        "OTFCumThresh",
        0.65,
        "skewed",
        True,
        "hannWinBounds",
        [0.8, 1.2],
        "DeconIter",
        3,
        "fixIter",
        True,
        "EdgeErosion",
        0,
        "Save16bit",
        True,
        "cpusPerTask",
        4,
    ],
    root=_WF,
)
# GPU jobs in local mode on an A100 (R04's GPU path gives MATLAB's gpuArray bits there)
CASES += wrapper_cases("gpu.rl", ["GPUJob", True], root=make_root(n_t=2), device="gpu:A100")
CASES += wrapper_cases("gpu.omw", [*_OMW, "GPUJob", True], root=_NOISY, device="gpu:A100")
# errors
CASES += wrapper_cases("error.missing_psf", ["psfFullpaths", ["psf/none.tif"]], expect_error=True)


def _gpu_name() -> str:
    from petakit5d.gpu.device import select_device

    return select_device().name


def _check(case, tmp_path):
    if case.device.startswith("gpu:"):
        name = _gpu_name()
        if case.device[4:].lower() not in name.lower():
            pytest.skip(f"the reference is MATLAB's on a {case.device[4:]}; this GPU is a {name}")
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
# the function string: MATLAB's rounding in transit
# ---------------------------------------------------------------------------------------------

# MATLAB R2024b's function strings (probe p1_funcstr at the oracle commit) of three value sets
_PROBE_BASE = {
    "frame": "/w/data/a.tif",
    "decon_path": "/w/data/matlab_decon/",
    "psf": "/w/psf/p.tif",
    "erode_mask": "",
    "save_maskfile": False,
    "rotate": False,
    "rl_method": "omw",
    "save_step": 5,
    "psf_gen": True,
    "save_zarr": False,
    "parse_cluster": False,
    "parse_parfor": False,
    "gpu_job": False,
    "save_16bit": True,
    "large_file": False,
    "large_method": "inmemory",
    "batch_size": "[1024,1024,1024]",
    "block_size": "[256,256,256]",
    "mask_fullpaths": "{''}",
    "uuid": "u",
    "cpus_per_task": 2,
    "mcc_mode": False,
    "config_file": "",
    "gpu_config_file": "",
    "debug": False,
    "flip_zstack": False,
}
_TAIL = (
    "'saveStep',5,'psfGen',true,'saveZarr',false,'parseCluster',false,'parseParfor',false,"
    "'GPUJob',false,'save16bit',true,'largeFile',false,'largeMethod','inmemory','batchSize',"
    "[1024,1024,1024],'blockSize',[256,256,256],"
)
_END = (
    "'maskFullpaths',{''},'uuid','u','cpusPerTask',2,'mccMode',false,'configFile','',"
    "'GPUConfigFile','')"
)
_HEAD = "XR_RLdeconFrame3D('/w/data/a.tif',"
_MID = "'/w/data/matlab_decon/','PSFfile','/w/psf/p.tif','dzPSF',"
PROBE_STRINGS = [
    (
        {
            "xy_pixel_size": 0.108,
            "dz": 0.3,
            "dz_psf": 0.3,
            "background": 100,
            "skew_angle": 32.45,
            "edge_erosion": 0,
            "decon_iter": 15,
            "wiener_alpha": 0.005,
            "otf_cum_thresh": 0.9,
            "hann_win_bounds": (0.8, 1),
            "skewed": None,
            "damp_factor": 1,
            "scale_factor": 1,
            "decon_offset": 0,
        },
        _HEAD + "0.1080000000,0.3000000000," + _MID + "0.3000000000,'background',[100],"
        "'skewAngle',3.245000e+01,'flipZstack',false,'edgeErosion',0,'ErodeMaskfile','',"
        "'SaveMaskfile',false,'Rotate',false,'deconIter',15,'RLMethod','omw','wienerAlpha',"
        "0.00500000000000000010,'OTFCumThresh',0.90000000000000002220,'hannWinBounds',[0.8,1],"
        "'skewed',[],'debug',false,"
        + _TAIL
        + "'dampFactor',1,'scaleFactor',[1],'deconOffset',0,"
        + _END,
    ),
    (
        {
            "xy_pixel_size": 0.10812345678901234,
            "dz": 0.20000000001234,
            "dz_psf": 0.123456789012345,
            "background": 99.123456789,
            "skew_angle": -32.4512345678,
            "edge_erosion": 2.5,
            "decon_iter": 3,
            "wiener_alpha": 0.00123456789012345678,
            "otf_cum_thresh": 0.87654321098765432,
            "hann_win_bounds": (0.81234567890123456, 1.2),
            "skewed": True,
            "damp_factor": 1.2345678901,
            "scale_factor": 2.5,
            "decon_offset": 0.123456789,
        },
        _HEAD + "0.1081234568,0.2000000000," + _MID + "0.1234567890,'background',"
        "[9.912346e+01],'skewAngle',-3.245123e+01,'flipZstack',false,'edgeErosion',2.500000e+00,"
        "'ErodeMaskfile','','SaveMaskfile',false,'Rotate',false,'deconIter',3,'RLMethod','omw',"
        "'wienerAlpha',0.00123456789012345671,'OTFCumThresh',0.87654321098765430875,"
        "'hannWinBounds',[0.812345678901235,1.2],'skewed',[true],'debug',false,"
        + _TAIL
        + "'dampFactor',1.234568e+00,'scaleFactor',[2.500000e+00],'deconOffset',1.234568e-01,"
        + _END,
    ),
    (
        {
            "xy_pixel_size": 1e-7,
            "dz": 12345.678901234567,
            "dz_psf": 1e20,
            "background": 1e30,
            "skew_angle": 1e300,
            "edge_erosion": 1e-300,
            "decon_iter": 2**60,
            "wiener_alpha": 1e-25,
            "otf_cum_thresh": 1 / 3,
            "hann_win_bounds": (1 / 3, 2 / 3),
            "skewed": False,
            "damp_factor": -0.5,
            "scale_factor": math.inf,
            "decon_offset": math.nan,
        },
        _HEAD + "0.0000001000,12345.6789012346," + _MID + "100000000000000000000.0000000000,"
        "'background',[1.000000e+30],'skewAngle',1.000000e+300,'flipZstack',false,'edgeErosion',"
        "1.000000e-300,'ErodeMaskfile','','SaveMaskfile',false,'Rotate',false,'deconIter',"
        "1152921504606846976,'RLMethod','omw','wienerAlpha',0.00000000000000000000,"
        "'OTFCumThresh',0.33333333333333331483,'hannWinBounds',[0.333333333333333,"
        "0.666666666666667],'skewed',[false],'debug',false," + _TAIL + "'dampFactor',"
        "-5.000000e-01,'scaleFactor',[Inf],'deconOffset',NaN," + _END,
    ),
]


@pytest.mark.parametrize("i", range(len(PROBE_STRINGS)))
def test_function_string_is_matlabs(i):
    values, expected = PROBE_STRINGS[i]
    v = dict(_PROBE_BASE, **values)
    v["hann_win_bounds"] = D._hann_text(v["hann_win_bounds"])
    assert D.frame_task_string(v) == expected


def test_long_decimals_reach_the_file_rounded_as_in_matlab():
    """the values decon_file gets are those of MATLAB's evaluation of the string"""
    kw = D.frame_task(PROBE_STRINGS[1][1]).kwargs
    assert kw["frame_fullpaths"] == "/w/data/a.tif"
    assert kw["xy_pixel_size"] == 0.1081234568 and kw["dz"] == 0.2
    assert kw["dz_psf"] == 0.123456789
    assert kw["background"] == 99.12346 and kw["skew_angle"] == -32.45123
    assert kw["edge_erosion"] == 2.5 and kw["decon_iter"] == 3
    assert kw["wiener_alpha"] == 0.00123456789012345671
    assert kw["otf_cum_thresh"] == 0.87654321098765430875
    assert kw["hann_win_bounds"] == [0.812345678901235, 1.2]
    assert kw["damp_factor"] == 1.234568 and kw["scale_factor"] == 2.5
    assert kw["decon_offset"] == 0.1234568
    assert kw["skewed"] is True and kw["mask_fullpaths"] == [""]
    # batchSize and blockSize come back as (z, y, x)
    assert kw["batch_size"] == [1024, 1024, 1024] and kw["block_size"] == [256, 256, 256]
    kw = D.frame_task(PROBE_STRINGS[2][1]).kwargs
    # (an alpha below 1e-20 reaches the file as 0)
    assert kw["wiener_alpha"] == 0 and kw["skewed"] is False and kw["decon_offset"] != 0
    assert D.frame_task(PROBE_STRINGS[0][1]).kwargs["skewed"] is None


def test_sizes_reach_the_file_in_their_axes():
    v = dict(_PROBE_BASE, **PROBE_STRINGS[0][0])
    v["hann_win_bounds"] = D._hann_text(v["hann_win_bounds"])
    v["batch_size"] = D._mat2str_commas(D._size_matlab((16, 64, 32)))
    v["block_size"] = D._mat2str_commas(D._size_matlab((8, 32, 16)))
    s = D.frame_task_string(v)
    assert "'batchSize',[64,32,16],'blockSize',[32,16,8]" in s
    kw = D.frame_task(s).kwargs
    assert kw["batch_size"] == [16, 64, 32] and kw["block_size"] == [8, 32, 16]


def test_a_vector_where_matlab_formats_one_number_raises():
    v = dict(_PROBE_BASE, **PROBE_STRINGS[0][0])
    v["hann_win_bounds"] = D._hann_text(v["hann_win_bounds"])
    v["xy_pixel_size"] = (0.108, 0.5)
    with pytest.raises(D.DeconDatasetError, match="shift"):
        D.frame_task_string(v)


def test_several_files_in_one_string():
    s = "XR_RLdeconFrame3D({'/w/a.tif','/w/b.tif'}" + PROBE_STRINGS[0][1][len(_HEAD) - 1 :]
    assert D.frame_task(s).kwargs["frame_fullpaths"] == ["/w/a.tif", "/w/b.tif"]


# ---------------------------------------------------------------------------------------------
# DeconWrapperTest in intent: local mode
# ---------------------------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def fresh_session(tmp_path_factory, monkeypatch):
    """Slurm hidden unless a test asks for the fake; a new session of single jobs"""
    hide_slurm(monkeypatch, tmp_path_factory.mktemp("no_slurm"))
    for name in (*SLURM_JOB_VARIABLES, "SLURM_CPUS_PER_TASK", "PETAKIT5D_RUNNER_STATE"):
        monkeypatch.delenv(name, raising=False)
    single.reset()
    config.clear_cache()
    rl.clear_otf_cache()
    yield
    single.reset()
    config.clear_cache()


class DS:
    """a dataset made by make_root in a folder"""

    def __init__(self, root: Path, cams=("A",), n_t=1, **kw):
        self.root = root
        make_root(cams=cams, n_t=n_t, **kw)(root)
        self.data = str(root / "data") + "/"
        self.cams = cams
        self.channels = [f"Cam{c}_ch0" for c in cams]
        self.psfs = [str(root / "psf" / f"PSF_Cam{c}.tif") for c in cams]
        self.fsns = [fsn(c, t) for c in cams for t in range(n_t)]
        ext = ".zarr" if kw.get("zarr") else ".tif"
        self.raws = [self.data + f + ext for f in self.fsns]

    def out(self, name, folder="matlab_decon", ext=".tif"):
        return f"{self.data}{folder}/{name}{ext}"


def run_wrapper(d: DS, **kw):
    args = {
        "channel_patterns": d.channels,
        "psf_fullpaths": d.psfs,
        "xy_pixel_size": XY,
        "dz": DZ,
        "dz_psf": DZ,
        "background": BG,
        "edge_erosion": 0,
        "save_16bit": True,
        "parse_cluster": False,
        "gpu_job": False,
    }
    args.update(kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D.decon_dataset([d.data], **args)


def _read(p):
    from petakit5d.io import read_tiff, read_zarr

    p = str(p)
    return read_zarr(p) if p.rstrip("/").endswith(".zarr") else read_tiff(p)


def _norm(psf):
    from petakit5d.kernels import decon as K

    p = np.asarray(psf, np.float32)
    return (p / np.float64(K.sum_double(p))).astype(np.float32)


def test_simplified_writes_outputs_that_match_the_function(tmp_path):
    d = DS(tmp_path)
    run_wrapper(d, rl_method="simplified", decon_iter=5)
    out = _read(d.out(d.fsns[0]))
    assert out.dtype == np.uint16 and out.shape == IM
    dp = d.data + "matlab_decon/"
    np.testing.assert_array_equal(_read(dp + f"MIPs/{d.fsns[0]}_MIP_z.tif"), out.max(axis=0))
    psf = _read(dp + "psfgen/PSF_CamA_simplified.tif")
    assert os.path.isfile(dp + "parameters.mat")
    params = json.loads(Path(dp + "parameters.json").read_text())
    assert params["RLMethod"] == "simplified"
    rl.clear_otf_cache()
    exp = rl.decon_lucy(_read(d.raws[0]), _norm(psf), 5, background=BG, use_gpu=False).decon
    np.testing.assert_array_equal(out, exp)


def test_omw_writes_back_projector_files_and_matches_the_function(tmp_path):
    d = DS(tmp_path)
    run_wrapper(
        d, rl_method="omw", decon_iter=2, wiener_alpha=0.005, otf_cum_thresh=0.9, skewed=False
    )
    pg = d.data + "matlab_decon/psfgen/PSF_CamA"
    bp = pg + "_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif"
    for f in (pg + "_omw.tif", bp, pg + "_OTF.tif", pg + "_OTF_mask_otf_cum_thresh_0.900000.tif"):
        assert os.path.isfile(f), f
    rl.clear_otf_cache()
    exp = rl.decon_lucy_omw(
        _read(d.raws[0]), _norm(_read(pg + "_omw.tif")), _read(bp), 2, background=BG, use_gpu=False
    ).decon
    np.testing.assert_array_equal(_read(d.out(d.fsns[0])), exp)


def test_decon_rotate_writes_the_deskewed_rotated_result(tmp_path):
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.kernels import ops

    d = DS(tmp_path)
    run_wrapper(d, decon_iter=2, decon_rotate=True)
    decon = _read(d.out(d.fsns[0]))
    exp = ops.convert(
        deskew_rotate_volume(decon, 32.45, DZ, XY, reverse=True, crop=True, save_16bit=True),
        "uint16",
    )
    dp = d.data + "matlab_decon/DSR/"
    np.testing.assert_array_equal(_read(dp + d.fsns[0] + ".tif"), exp)
    np.testing.assert_array_equal(_read(dp + f"MIPs/{d.fsns[0]}_MIP_z.tif"), exp.max(axis=0))


def test_save16bit_false_writes_single(tmp_path):
    d = DS(tmp_path)
    run_wrapper(d, decon_iter=3, save_16bit=False)
    out = _read(d.out(d.fsns[0]))
    assert out.dtype == np.float32
    psf = _read(d.data + "matlab_decon/psfgen/PSF_CamA_simplified.tif")
    rl.clear_otf_cache()
    exp = rl.decon_lucy(
        _read(d.raws[0]), _norm(psf), 3, background=BG, use_gpu=False, save_16bit=False
    ).decon
    np.testing.assert_array_equal(out, exp)


def test_zarr_output_matches_tiff_output(tmp_path):
    d = DS(tmp_path)
    run_wrapper(d, decon_iter=3, result_dir_name="decon_tif")
    run_wrapper(d, decon_iter=3, result_dir_name="decon_zarr", save_zarr=True)
    np.testing.assert_array_equal(
        _read(d.out(d.fsns[0], "decon_zarr", ".zarr")), _read(d.out(d.fsns[0], "decon_tif"))
    )


def test_each_channel_uses_its_own_psf_and_background(tmp_path):
    d = DS(tmp_path, cams=("A", "B"), sigmas={"B": (3.0, 2.5)})
    bg = [100, 80]
    run_wrapper(d, decon_iter=3, background=bg)
    for c, cam in enumerate(d.cams):
        psf = _read(d.data + f"matlab_decon/psfgen/PSF_Cam{cam}_simplified.tif")
        rl.clear_otf_cache()
        exp = rl.decon_lucy(_read(d.raws[c]), _norm(psf), 3, background=bg[c], use_gpu=False).decon
        np.testing.assert_array_equal(_read(d.out(d.fsns[c])), exp, err_msg=f"channel {cam}")


def test_omw_parameters_are_per_channel(tmp_path):
    d = DS(tmp_path, cams=("A", "B"))
    run_wrapper(
        d, rl_method="omw", decon_iter=1, wiener_alpha=[0.05, 0.005], otf_cum_thresh=[0.88, 0.9],
        skewed=False,
    )  # fmt: skip
    pg = d.data + "matlab_decon/psfgen/"
    assert os.path.isfile(pg + "PSF_CamA_back_projector_alpha_0.050000_otf_cum_thresh_0.880000.tif")
    assert os.path.isfile(pg + "PSF_CamB_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif")


def test_erode_by_first_time_point_applies_its_mask_to_all_files(tmp_path):
    from petakit5d.decon.erosion import decon_mask_edge_erosion

    d = DS(tmp_path, n_t=2, zero_cols=True)
    e = 2
    run_wrapper(d, decon_iter=3, edge_erosion=e, erode_by_ftp=True)
    mask = _read(d.data + f"matlab_decon/Masks/{d.fsns[0]}_eroded.zarr") > 0
    np.testing.assert_array_equal(mask, decon_mask_edge_erosion(_read(d.raws[0]) > 0, e))
    for t in range(2):
        out = _read(d.out(d.fsns[t]))
        assert np.count_nonzero(out[~mask]) == 0, f"time point {t}"
        assert np.count_nonzero(out[mask]) > 0


def test_a_mask_of_another_size_is_not_used(tmp_path):
    """(the first time point is narrower) MATLAB's warning, and the file's own erosion"""
    d = DS(tmp_path, n_t=2)
    write_tiff(d.raws[0], _read(d.raws[0])[:, :, :-4])
    with pytest.warns(
        UserWarning,
        match=r"image size \[48 44 32\] does not match the defined mask size \[48 40 32\]",
    ):
        D.decon_dataset(
            [d.data], channel_patterns=d.channels, psf_fullpaths=d.psfs, xy_pixel_size=XY, dz=DZ,
            dz_psf=DZ, background=BG, edge_erosion=2, parse_cluster=False, decon_iter=2,
        )  # fmt: skip
    assert os.path.isfile(d.out(d.fsns[1]))


def test_a_finished_first_time_point_without_its_mask_is_deconvolved_again(tmp_path):
    d = DS(tmp_path, n_t=2)
    run_wrapper(d, decon_iter=2)
    first = d.out(d.fsns[0])
    os.utime(first, (1, 1))
    run_wrapper(d, decon_iter=2, edge_erosion=2)
    assert os.stat(first).st_mtime > 1, "The first time point was not deconvolved again."
    assert os.path.isdir(d.data + f"matlab_decon/Masks/{d.fsns[0]}_eroded.zarr")


def test_widefield_demo_settings_with_plain_filenames(tmp_path):
    from petakit5d import matlab
    from petakit5d.testing.synth import centroid

    data = tmp_path / "Widefield"
    data.mkdir()
    pos = (16, 24, 22)
    im = beads(IM, np.array([pos], float), BLUR, 1000.0, float(BG))
    write_tiff(str(data / "WF_raw.tif"), np.floor(im + np.float32(0.5)).astype(np.uint16))
    write_tiff(str(data / "WF_PSF.tif"), psf_raw())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        matlab.XR_decon_data_wrapper(
            [str(data) + "/"], "resultDirName", "matlab_decon_omw", "xyPixelSize", 0.157, "dz",
            0.1, "Reverse", True, "ChannelPatterns", ["WF_raw"], "PSFFullpaths",
            [str(data / "WF_PSF.tif")], "dzPSF", 0.1, "parseSettingFile", False, "RLmethod",
            "omw", "wienerAlpha", 0.02, "OTFCumThresh", 0.65, "skewed", True, "hannWinBounds",
            [0.8, 1.2], "Background", 100, "DeconIter", 3, "fixIter", True, "EdgeErosion", 0,
            "Save16bit", True, "zarrFile", False, "saveZarr", False, "parseCluster", False,
            "largeFile", False, "GPUJob", False, "debug", False, "cpusPerTask", 4, "ConfigFile",
            "", "GPUConfigFile", "", "mccMode", False,
        )  # fmt: skip
    out = _read(data / "matlab_decon_omw" / "WF_raw.tif")
    assert out.shape == IM
    assert out.max() > 1.5 * 1000, "The bead was not sharpened."
    np.testing.assert_allclose(centroid(out.astype(np.float32), pos, 6), pos, atol=0.1)


def test_deconvolved_data_can_be_deskewed_and_rotated(tmp_path):
    """demo_RL_deconvolution step 5 (the dataset wrapper is D04's; here its file function)"""
    from petakit5d.deskew_rotate import deskew_rotate_file

    d = DS(tmp_path)
    run_wrapper(d, decon_iter=3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for f in (d.raws[0], d.out(d.fsns[0])):
            deskew_rotate_file(f, XY, DZ, reverse=True)
    raw = _read(d.data + f"DSR/{d.fsns[0]}.tif")
    decon = _read(d.data + f"matlab_decon/DSR/{d.fsns[0]}.tif")
    assert decon.shape == raw.shape and decon.max() > raw.max()


def test_parameters_are_saved_as_given(tmp_path):
    from petakit5d.io import read_mat

    d = DS(tmp_path)
    run_wrapper(d, decon_iter=2, batch_size=(16, 64, 32), skewed=None)
    dp = d.data + "matlab_decon/"
    pr = json.loads(Path(dp + "parameters.json").read_text())
    assert pr["deconIter"] == 2 and pr["batchSize"] == [64, 32, 16]
    assert pr["dataPaths"] == [d.data] and pr["skewed"] == [] and pr["uuid"] == ""
    assert pr["channelPatterns"] == d.channels and pr["psfFullpaths"] == d.psfs
    mat = read_mat(dp + "parameters.mat")["pr"]
    assert list(np.ravel(mat["batchSize"])) == [64, 32, 16]


def test_missing_psf_raises(tmp_path):
    d = DS(tmp_path)
    with pytest.raises(D.DeconDatasetError, match="PSF file .*none.tif does not exist!"):
        run_wrapper(d, psf_fullpaths=[str(tmp_path / "none.tif")])


def test_a_file_of_two_channels_raises(tmp_path):
    """MATLAB's function string would hold both channels' values"""
    d = DS(tmp_path)
    with pytest.raises(D.DeconDatasetError, match="2 channel patterns"):
        run_wrapper(d, channel_patterns=["CamA_ch0", "CamA"], psf_fullpaths=d.psfs * 2)


def test_a_young_tmp_marker_is_another_process_file(tmp_path, monkeypatch):
    """(local mode) a file whose .tmp marker is younger than unit_wait_time is left to the
    process that made it; an older one is taken over"""
    d = DS(tmp_path, n_t=2)
    os.makedirs(d.data + "matlab_decon")
    marker = d.data + f"matlab_decon/{d.fsns[0]}.tmp"
    Path(marker).touch()
    waits = []

    def wait(groups, t):
        waits.append(t)
        os.utime(marker, (1, 1))
        return False

    monkeypatch.setattr(jfiles, "wait_for_files", wait)
    run_wrapper(d, decon_iter=2)
    assert waits == [30.0]
    assert os.path.isfile(d.out(d.fsns[0])) and not os.path.exists(marker)


def test_zarr_results_have_matlabs_marker_name(tmp_path, monkeypatch):
    d = DS(tmp_path)
    seen = []
    real = D._touch
    monkeypatch.setattr(D, "_touch", lambda p: (seen.append(p), real(p)))
    run_wrapper(d, decon_iter=1, save_zarr=True)
    assert seen and seen[0].endswith(f"/{d.fsns[0]}..tmp")
    assert not os.path.exists(seen[0])


@pytest.mark.gpu
def test_wrapper_gpu_job_matches_cpu(tmp_path):
    """(DeconGPUTest/wrapperGPUJobMatchesCPU) uint16 results within rounding (1 count)"""
    d = DS(tmp_path, noisy=True)
    for gpu in (False, True):
        run_wrapper(
            d, rl_method="omw", decon_iter=2, skewed=False, gpu_job=gpu,
            result_dir_name=f"decon_gpu_{int(gpu)}",
        )  # fmt: skip
    a = _read(d.out(d.fsns[0], "decon_gpu_1")).astype(np.float64)
    b = _read(d.out(d.fsns[0], "decon_gpu_0")).astype(np.float64)
    assert np.all(np.abs(a - b) <= np.maximum(1, 1e-3 * np.abs(b)))


def test_original_method_is_phase_two(tmp_path):
    d = DS(tmp_path)
    with pytest.raises(NotImplementedError):
        run_wrapper(d, rl_method="original", decon_iter=2)


# ---------------------------------------------------------------------------------------------
# DeconWrapperTest in intent: cluster mode against the fake Slurm
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def pauses(monkeypatch):
    """the waits' pauses, without waiting (the fake clock moves by each pause)"""
    now = [0.0]
    calls: list[float] = []

    def pause(t):
        calls.append(t)
        now[0] += t

    monkeypatch.setattr(jfiles, "_pause", pause)
    monkeypatch.setattr(jfiles, "_clock", lambda: now[0])
    return calls


def gpu_config(folder: Path, **over) -> str:
    """a GPU job config like the production one for 16-core nodes with 4 GPUs (made-up names)"""
    cfg = {
        "SlurmParam": "-p gpu_test -n1 --mem-per-cpu=32128M --gres=gpu:1",
        "MemPerCPU": 31.375,
        "maxCPUNum": 16,
        "minCPUNum": 4,
    }
    cfg.update(over)
    n = len(list(folder.glob("gpu_config_*.json")))
    fn = folder / f"gpu_config_{n}.json"
    fn.write_text(json.dumps(cfg))
    return str(fn)


def cluster(d: DS, logs: Path | None = None, **kw):
    logs = logs or d.root / "logs"
    logs.mkdir(exist_ok=True)
    kw.setdefault("parse_cluster", True)
    kw.setdefault("job_log_dir", str(logs))
    run_wrapper(d, **kw)


@pytest.fixture
def submitted(monkeypatch):
    """the task record of each sbatch submission, in order (the call's task folders go with it)"""
    out: list[dict] = []
    real = single.submit_single_job

    def spy(task, *a, **kw):
        r = real(task, *a, **kw)
        if r[2] > 0:
            out.append(task.record())
        return r

    monkeypatch.setattr(single, "submit_single_job", spy)
    return out


def job_files(records) -> list[list[str]]:
    res = []
    for t in records:
        f = t["kwargs"]["frame_fullpaths"]
        f = [f] if isinstance(f, str) else f
        res.append([os.path.splitext(os.path.basename(x))[0] for x in f])
    return res


def cpus_of(fake) -> list[int]:
    return [int(re.search(r"--cpus-per-task=(\d+)", c).group(1)) for c in fake.calls("sbatch")]


@POSIX
def test_cluster_mode_checks_few_jobs_soon(tmp_path, fake_slurm, monkeypatch):
    d = DS(tmp_path)
    waits = []
    monkeypatch.setattr(jfiles, "wait_for_files", lambda g, t: waits.append(t) or False)
    cluster(d, max_trial_num=2)
    assert fake_slurm.call_count("sbatch") == 2
    assert waits == [3.0, 3.0]


@POSIX
def test_cluster_mode_ends_soon_after_the_last_result(
    tmp_path, fake_slurm, pauses, monkeypatch, submitted
):
    d = DS(tmp_path)
    out = d.out(d.fsns[0])
    n = [0]
    real = jfiles.batch_file_exist

    def look(paths, *a, **kw):
        n[0] += 1
        if n[0] == 2:
            Path(out).touch()
        return real(paths, *a, **kw)

    monkeypatch.setattr(jfiles, "batch_file_exist", look)
    cluster(d, max_trial_num=2)
    assert fake_slurm.call_count("sbatch") == 1, "The job was submitted again."
    assert pauses == [1.0, 1.0], "The wrapper did not go on within a second of the result."


@POSIX
def test_cluster_mode_submits_files_soon_after_their_mask(
    tmp_path, fake_slurm, pauses, monkeypatch
):
    from petakit5d.io import write_zarr

    d = DS(tmp_path, n_t=2)
    mask = d.data + f"matlab_decon/Masks/{d.fsns[0]}_eroded.zarr"
    os.makedirs(os.path.dirname(mask))
    fake_slurm.set_job_states("1000_1", ["R", "done"], [d.out(d.fsns[0])])
    calls = [0]

    def pause(t):
        calls[0] += 1
        pauses.append(t)
        if calls[0] == 2:
            write_zarr(np.ones(IM, np.uint8), mask)

    monkeypatch.setattr(jfiles, "_pause", pause)
    cluster(d, edge_erosion=2, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 2
    # two seconds until the mask, then a whole wait for the second file's job
    assert pauses[:2] == [1.0, 1.0] and sum(pauses[2:]) == pytest.approx(3.0)


@POSIX
def test_cluster_mode_sees_masks_other_nodes_just_made(
    tmp_path, fake_slurm, pauses, monkeypatch, submitted
):
    """the mask check is the wait's (batch_file_exist): here only it sees the mask"""
    from petakit5d import io

    d = DS(tmp_path, n_t=2)
    real = jfiles.batch_file_exist
    monkeypatch.setattr(
        jfiles,
        "batch_file_exist",
        lambda paths, *a, **kw: [
            str(p).endswith("_eroded.zarr") or r
            for p, r in zip(paths, real(paths, *a, **kw), strict=True)
        ],
    )
    sz = io.get_image_size
    monkeypatch.setattr(io, "get_image_size", lambda p: IM if p.endswith("_eroded.zarr") else sz(p))
    cluster(d, edge_erosion=2, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 2
    assert sum(pauses) == pytest.approx(3.0), (
        "The second file was not submitted in the first round."
    )


@POSIX
def test_cluster_mode_waits_for_the_last_try(tmp_path, fake_slurm, monkeypatch):
    d = DS(tmp_path)
    fake_slurm.set_job_states("1000_1", ["PD", "R", "R", "done"])
    waits = []
    monkeypatch.setattr(jfiles, "wait_for_files", lambda g, t: waits.append(t) or False)
    cluster(d, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 1, "The last try was submitted again."
    # (status checks; the call's end looks at the jobs of its task folders once more)
    assert sum(c.startswith("-j ") for c in fake_slurm.calls("squeue")) == 4
    assert len(waits) == 4


@POSIX
def test_cluster_mode_waits_for_the_last_try_of_a_job_of_several_files(
    tmp_path, fake_slurm, pauses, submitted
):
    d = DS(tmp_path, n_t=2)
    fake_slurm.set_job_states("1000_1", ["PD", "R", "R", "done"])
    cluster(d, max_trial_num=1, task_batch_num=2)
    assert fake_slurm.call_count("sbatch") == 1, "The last try was submitted again."
    checks = [c for c in fake_slurm.history() if c.startswith("squeue -j") and "1000_1" in c]
    assert len(checks) == 4, "The job was not checked until it ended."


@POSIX
def test_cluster_mode_wait_grows_with_unfinished_files(tmp_path, fake_slurm, monkeypatch):
    d = DS(tmp_path, n_t=20)
    waits = []
    monkeypatch.setattr(jfiles, "wait_for_files", lambda g, t: waits.append(t) or False)
    cluster(d, max_trial_num=2)
    assert fake_slurm.call_count("sbatch") == 40
    assert waits == pytest.approx([6.0, 6.0])


def _claim(logs: Path, i: int, job: str) -> None:
    for c in logs.glob(f"{single.POOL_PREFIX}*/claims"):
        try:
            os.symlink(f"{job}:othernode:5:1", c / f"t{i}")
        except FileExistsError:
            pass


@pytest.mark.skipif(
    not os.sys.platform.startswith("linux"), reason="Task runners run on Linux only."
)
@POSIX
def test_last_try_waits_for_a_task_runner_that_has_the_task(
    tmp_path, fake_slurm, monkeypatch, submitted
):
    d = DS(tmp_path, n_t=3)
    logs = tmp_path / "logs"
    logs.mkdir()
    outs = [d.out(f) for f in d.fsns]
    fake_slurm.set_job_states("1000_1", ["R", "done"], outs[:1])
    fake_slurm.set_job_states("1001_2", ["R", "R", "R", "R", "done"], outs[1:2])
    monkeypatch.setattr(jfiles, "_pause", lambda t: _claim(logs, 2, "1001_2"))
    cluster(d, logs, max_trial_num=1)
    assert fake_slurm.call_count("sbatch") == 3, "File 3 was submitted again."
    checks = [c for c in fake_slurm.history() if c.startswith("squeue") and "1001_2" in c]
    assert len(checks) >= 5, "The wrapper did not wait for the job that had file 3."
    assert not list(logs.glob(f"{single.POOL_PREFIX}*")), "The pool is still there."


@pytest.mark.skipif(
    not os.sys.platform.startswith("linux"), reason="Task runners run on Linux only."
)
@POSIX
def test_files_of_a_job_wait_for_a_task_runner_that_has_them(
    tmp_path, fake_slurm, monkeypatch, submitted
):
    d = DS(tmp_path, n_t=6)
    logs = tmp_path / "logs"
    logs.mkdir()
    outs = [d.out(f) for f in d.fsns]
    fake_slurm.set_job_states("1000_1", ["R", "done"], outs[0:2])
    fake_slurm.set_job_states("1001_3", ["R", "R", "R", "R", "done"], outs[2:4])
    monkeypatch.setattr(jfiles, "_pause", lambda t: _claim(logs, 2, "1001_3"))
    cluster(
        d, logs, gpu_job=True, task_batch_num=2, max_trial_num=2,
        gpu_config_file=gpu_config(tmp_path),
    )  # fmt: skip
    cmds = fake_slurm.history()
    submits = [i for i, c in enumerate(cmds) if c.startswith("sbatch")]
    assert len(submits) == 5
    running = [i for i, c in enumerate(cmds) if c.startswith("squeue") and "1001_3" in c]
    assert len(running) >= 5
    assert submits[3] > running[4], "Files were submitted again while the job that had them ran."
    assert single.POOL_PREFIX not in cmds[submits[0]], "The first job of its kind has task runners."
    assert all(single.POOL_PREFIX in cmds[i] for i in submits[1:3]), (
        "The pool jobs have no runners."
    )
    assert not any(single.POOL_PREFIX in cmds[i] for i in submits[3:]), "Resubmitted with runners."


@POSIX
@pytest.mark.parametrize(
    "case", [({"minCPUNum": 1}, True, 4), ({"minCPUNum": 8}, True, 8), ({}, False, 24)]
)
def test_gpu_jobs_take_at_most_four_cpus(tmp_path, fake_slurm, case):
    over, gpu, cpus = case
    d = DS(tmp_path, n_t=2)
    cluster(
        d, gpu_job=gpu, gpu_config_file=gpu_config(tmp_path, **over), cpus_per_task=24,
        max_trial_num=1,
    )  # fmt: skip
    assert cpus_of(fake_slurm) == [cpus, cpus]


@POSIX
def test_gpu_jobs_take_several_files_with_the_same_settings(tmp_path, fake_slurm, submitted):
    """and the jobs give the results of one file at a time (run here on the CPU)"""
    from petakit5d.jobs.tasks import run_record

    d = DS(tmp_path, cams=("A", "B"), n_t=3, sigmas={"B": (3.0, 2.5)})
    cluster(
        d, decon_iter=3, gpu_job=True, result_dir_name="decon_jobs",
        gpu_config_file=gpu_config(tmp_path), task_batch_num=2, max_trial_num=1,
    )  # fmt: skip
    assert job_files(submitted) == [d.fsns[0:2], d.fsns[2:3], d.fsns[3:5], d.fsns[5:6]]
    tasks = submitted
    assert all(t["kwargs"]["gpu_job"] is True for t in tasks)
    for t in tasks:
        t["kwargs"]["gpu_job"] = False
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run_record(t)
    run_wrapper(d, decon_iter=3, result_dir_name="decon_local")
    for f in d.fsns:
        np.testing.assert_array_equal(
            _read(d.out(f, "decon_jobs")), _read(d.out(f, "decon_local")), err_msg=f
        )


@POSIX
def test_running_job_is_not_submitted_again(tmp_path, fake_slurm, pauses, submitted):
    d = DS(tmp_path, n_t=2)
    fake_slurm.set_job_states("1000_1", ["R", "R", "done"])
    cluster(
        d, gpu_job=True, task_batch_num=2, max_trial_num=2, gpu_config_file=gpu_config(tmp_path)
    )
    assert job_files(submitted) == [d.fsns, d.fsns[0:1], d.fsns[1:2]]
    assert sum(c.startswith("-j ") and "1000_1" in c for c in fake_slurm.calls("squeue")) == 3


@POSIX
def test_retried_files_get_a_job_each(tmp_path, fake_slurm, pauses, submitted):
    d = DS(tmp_path, n_t=2)
    cluster(
        d, gpu_job=True, task_batch_num=2, max_trial_num=2, gpu_config_file=gpu_config(tmp_path)
    )
    assert job_files(submitted) == [d.fsns, d.fsns[0:1], d.fsns[1:2]]


@POSIX
def test_files_of_a_job_that_just_ended_are_looked_at_again(
    tmp_path, fake_slurm, pauses, monkeypatch, submitted
):
    """(jobs of several files) the round's listing missed the results (a folder listed before
    another node's job made them); the job is gone at the round's query, and the files are
    listed again before they are submitted again"""
    d = DS(tmp_path, n_t=2)
    fake_slurm.set_job_states("1000_1", ["R", "done"], [d.out(f) for f in d.fsns])
    real = D._listing
    calls = [0]

    def listing(folder):
        calls[0] += 1
        return set() if calls[0] <= 3 else real(folder)

    monkeypatch.setattr(D, "_listing", listing)
    # (and a look by name misses them too, as on NFS)
    exist = D._exist
    monkeypatch.setattr(D, "_exist", lambda p: "/matlab_decon/" not in p and exist(p))
    cluster(
        d, gpu_job=True, task_batch_num=2, max_trial_num=2, gpu_config_file=gpu_config(tmp_path)
    )
    assert job_files(submitted) == [d.fsns], "Finished files were submitted again."


@POSIX
def test_the_job_log_names_come_from_the_last_dataset(tmp_path, fake_slurm, pauses):
    """(MATLAB's checkSlurmCluster after the loop) a job log folder that doesn't exist is
    replaced by the last dataset's job_logs"""
    d = DS(tmp_path / "a")
    d2 = DS(tmp_path / "b")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D.decon_dataset(
            [d.data, d2.data], channel_patterns=d.channels, psf_fullpaths=d.psfs,
            xy_pixel_size=XY, dz=DZ, dz_psf=DZ, background=BG, parse_cluster=True,
            max_trial_num=1, job_log_dir=str(tmp_path / "none"),
        )  # fmt: skip
    cmds = fake_slurm.calls("sbatch")
    assert len(cmds) == 2
    assert all(f" -o {tmp_path}/b/data/job_logs/" in c for c in cmds)


def test_each_dataset_erodes_by_its_own_first_time_point(tmp_path):
    from petakit5d.decon.erosion import decon_mask_edge_erosion

    make_root(datasets=("d1", "d2"), n_t=2, zero_cols=True)(tmp_path)
    data = [str(tmp_path / x) + "/" for x in ("d1", "d2")]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D.decon_dataset(
            data, channel_patterns=["CamA_ch0"], psf_fullpaths=[str(tmp_path / "psf/PSF_CamA.tif")],
            xy_pixel_size=XY, dz=DZ, dz_psf=DZ, background=BG, parse_cluster=False, decon_iter=2,
            edge_erosion=2,
        )  # fmt: skip
    for p in data:
        mask = _read(p + f"matlab_decon/Masks/{fsn('A', 0)}_eroded.zarr") > 0
        raw = _read(p + fsn("A", 0) + ".tif")
        np.testing.assert_array_equal(mask, decon_mask_edge_erosion(raw > 0, 2))
        for t in range(2):
            assert np.count_nonzero(_read(p + f"matlab_decon/{fsn('A', t)}.tif")[~mask]) == 0


@POSIX
def test_files_of_different_sizes_go_to_different_jobs(tmp_path, fake_slurm, pauses, submitted):
    d = DS(tmp_path, n_t=3)
    write_tiff(d.raws[1], _read(d.raws[1])[:-4])
    cluster(
        d, gpu_job=True, task_batch_num=3, max_trial_num=1, gpu_config_file=gpu_config(tmp_path)
    )
    assert job_files(submitted) == [[d.fsns[0], d.fsns[2]], [d.fsns[1]]]


@POSIX
def test_gpu_jobs_take_up_to_eight_files_while_there_are_sixteen_jobs(
    tmp_path, fake_slurm, pauses, submitted
):
    d = DS(tmp_path, n_t=40)
    cluster(d, gpu_job=True, max_trial_num=1, gpu_config_file=gpu_config(tmp_path))
    assert [len(f) for f in job_files(submitted)] == [2] * 20


@POSIX
def test_omw_back_projector_is_made_once_before_the_jobs(tmp_path, fake_slurm, pauses):
    from petakit5d.decon.backprojector import generate_backprojector

    d = DS(tmp_path / "a", n_t=2)
    d2 = DS(tmp_path / "b")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        logs = tmp_path / "logs"
        logs.mkdir()
        D.decon_dataset(
            [d.data, d2.data], channel_patterns=d.channels, psf_fullpaths=d.psfs, xy_pixel_size=XY,
            dz=DZ, dz_psf=DZ, background=BG, rl_method="omw", decon_iter=2, wiener_alpha=0.005,
            otf_cum_thresh=0.9, skewed=False, parse_cluster=True, gpu_job=True, max_trial_num=1,
            gpu_config_file=gpu_config(tmp_path), job_log_dir=str(logs),
        )  # fmt: skip
    for data in (d.data, d2.data):
        pg = data + "matlab_decon/psfgen/PSF_CamA"
        bp = pg + "_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif"
        assert os.path.isfile(bp), "The wrapper did not make the back projector."
        r = generate_backprojector(
            _read(pg + "_omw.tif"), 0.005, False, 0.9, hann_win_bounds=(0.8, 1)
        )
        np.testing.assert_array_equal(_read(bp), r.b_omw)
        assert os.path.isfile(pg + "_OTF.tif")
        assert os.path.isfile(pg + "_OTF_mask_otf_cum_thresh_0.900000.tif")
        assert not os.path.exists(bp[:-4] + "_figure.png")


@POSIX
def test_the_premade_back_projector_gets_the_jobs_values(tmp_path, fake_slurm, pauses, monkeypatch):
    """the master makes the back projector with the values the jobs get from the function
    string (%.20f of alpha and threshold, mat2str of the bounds): a small alpha is rounded"""
    from petakit5d.core.text import through_text
    from petakit5d.decon import backprojector

    seen = []
    monkeypatch.setattr(backprojector, "omw_backprojector_files", lambda *a, **kw: seen.append(a))
    d = DS(tmp_path)
    alpha, thresh = 1.2345678901234567e-5, 0.8812345678901234567
    cluster(
        d, rl_method="omw", wiener_alpha=alpha, otf_cum_thresh=thresh, skewed=False,
        hann_win_bounds=(0.81234567890123456, 1.0), gpu_job=True, max_trial_num=1,
        gpu_config_file=gpu_config(tmp_path),
    )  # fmt: skip
    (a,) = seen
    assert a[3] == float(through_text(alpha, "%.20f")) != alpha
    assert a[4] == float(through_text(thresh, "%.20f"))
    assert a[6] == [0.812345678901235, 1] and a[8] is False


@POSIX
def test_a_job_of_several_files_gets_all_its_files(tmp_path, fake_slurm, pauses, submitted):
    """(MATLAB: compiledJobGetsAllFilesOfItsBatch) the task record holds the batch's files and
    the job's settings"""
    d = DS(tmp_path, n_t=2)
    cluster(
        d, gpu_job=True, mcc_mode=True, task_batch_num=2, max_trial_num=1,
        gpu_config_file=gpu_config(tmp_path),
    )  # fmt: skip
    (t,) = submitted
    assert t["func"] == D.FRAME_FUNC
    assert t["kwargs"]["frame_fullpaths"] == d.raws
    assert t["kwargs"]["gpu_job"] is True and t["kwargs"]["mcc_mode"] is True


@POSIX
def test_cluster_jobs_get_the_values_of_the_local_mode(
    tmp_path, fake_slurm, pauses, monkeypatch, submitted
):
    """the same function string, so the same values, in both modes"""
    from petakit5d.jobs import tasks

    d = DS(tmp_path)
    seen = []
    real = tasks.run_record
    monkeypatch.setattr(tasks, "run_record", lambda rec: seen.append(rec) or real(rec))
    run_wrapper(d, decon_iter=1, xy_pixel_size=0.10812345678901234, background=99.123456789)
    monkeypatch.setattr(tasks, "run_record", real)
    cluster(
        d, decon_iter=1, xy_pixel_size=0.10812345678901234, background=99.123456789,
        result_dir_name="decon_c", max_trial_num=1,
    )  # fmt: skip
    (job,) = submitted
    local = seen[0]["kwargs"]
    assert job["kwargs"]["xy_pixel_size"] == local["xy_pixel_size"] == 0.1081234568
    assert job["kwargs"]["background"] == local["background"] == 99.12346
    diff = {k for k in local if local[k] != job["kwargs"][k]}
    assert diff == {"decon_path", "parse_cluster", "uuid"}


@POSIX
def test_cluster_jobs_ask_for_the_memory_of_the_file(tmp_path, fake_slurm, pauses):
    """10 x 4 bytes per voxel: here a config with 1 GB per CPU and a 2 G-voxel header"""
    from petakit5d.io import create_zarr

    d = DS(tmp_path, zarr=True)
    cfg = tmp_path / "cpu.json"
    cfg.write_text(json.dumps({"MemPerCPU": 1, "maxCPUNum": 128}))
    import shutil

    shutil.rmtree(d.raws[0])
    create_zarr(d.raws[0], data_size=(512, 1024, 1024), block_size=(256, 256, 256), dtype="uint16")
    cluster(d, zarr_file=True, config_file=str(cfg), max_trial_num=1)
    assert cpus_of(fake_slurm) == [math.ceil(512 * 1024 * 1024 * 4 / 2**30 * 10)]


@POSIX
def test_cluster_mode_without_slurm_runs_here(tmp_path):
    d = DS(tmp_path)
    with pytest.warns(UserWarning, match="slurm-based computing cluster is not exist"):
        D.decon_dataset(
            [d.data], channel_patterns=d.channels, psf_fullpaths=d.psfs, xy_pixel_size=XY, dz=DZ,
            dz_psf=DZ, background=BG, decon_iter=1, parse_cluster=True,
        )  # fmt: skip
    assert os.path.isfile(d.out(d.fsns[0]))


# ---------------------------------------------------------------------------------------------
# in-place large files (the controllers; R07's functions)
# ---------------------------------------------------------------------------------------------


def _need_r07():
    from petakit5d import decon

    if not hasattr(
        decon, "in_place_controller_resources"
    ) and "in_place_controller_resources" not in getattr(decon, "_EXPORTS", {}):
        pytest.skip("R07's large-file functions are not in this tree")


@POSIX
@pytest.mark.parametrize(
    "case",
    [
        (False, True, None, 4, 4),
        (False, True, "off", 4, 1),
        (False, False, None, 4, 1),
        (True, True, None, 2, 1),
    ],
)
def test_in_place_controllers_ask_for_what_they_hold(
    tmp_path, fake_slurm, case, monkeypatch, submitted
):
    """(MATLAB: largeFileInPlaceControllersAskForWhatTheyHold, the sbatch line) a 34 G-voxel
    header; the controller's CPUs, its masterCompute and its nice"""
    _need_r07()
    import shutil

    from petakit5d.io import create_zarr
    from petakit5d.jobs import coordinator

    gpu, mc, cfg, cpus_task, cpus = case
    d = DS(tmp_path, zarr=True)
    shutil.rmtree(d.raws[0])
    create_zarr(d.raws[0], data_size=(2048, 4096, 4096), block_size=(256, 256, 256), dtype="uint16")
    if cfg:
        cfg = tmp_path / "cpu_config.json"
        cfg.write_text('{"masterCompute": false}')
    monkeypatch.setattr(coordinator, "coordinator_nice", lambda slurm=None: 1000)
    monkeypatch.setattr(jfiles, "wait_for_files", lambda g, t: False)
    cluster(
        d, zarr_file=True, save_zarr=True, large_file=True, large_method="inplace",
        gpu_job=gpu, master_compute=mc, config_file=str(cfg or ""), cpus_per_task=cpus_task,
        max_trial_num=1,
    )  # fmt: skip
    lines = fake_slurm.calls("sbatch")
    assert len(lines) == 1
    assert cpus_of(fake_slurm) == [cpus]
    (t,) = submitted
    assert t["kwargs"].get("master_compute", True) is mc
    assert t["kwargs"]["job_log_dir"] == str(tmp_path / "logs")
    assert " --nice=1000 " in lines[0], "The controller did not go before its blocks."


@POSIX
def test_cluster_mode_waits_for_a_controller_of_another_call(
    tmp_path, fake_slurm, monkeypatch, submitted
):
    """a writer record of another uuid whose job is running: no controller is submitted"""
    _need_r07()
    from petakit5d.jobs import writers

    d = DS(tmp_path, zarr=True)
    result = d.out(d.fsns[0], ext=".zarr")
    fake_slurm.set_job_states("777_1", ["R"])
    monkeypatch.setattr(
        writers, "live", lambda r, u: ["777_1:othernode:5"] if not os.path.isdir(result) else []
    )
    n = [0]

    def wait(groups, t):
        n[0] += 1
        if n[0] == 2:
            os.makedirs(result)
        return False

    monkeypatch.setattr(jfiles, "wait_for_files", wait)
    cluster(
        d, zarr_file=True, save_zarr=True, large_file=True, large_method="inplace",
        batch_size=(64, 64, 64), block_size=(16, 16, 16), max_trial_num=1,
    )  # fmt: skip
    assert fake_slurm.call_count("sbatch") == 0
