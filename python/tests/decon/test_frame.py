"""Deconvolution of the files of a job (petakit5d.decon.frame.decon_file) against MATLAB's
XR_RLdeconFrame3D (the oracle tree's, on the CPU).

Oracle cases: XR_RLdeconFrame3D on a folder with frames of blurred beads (DeconFrameTest's,
48 x 44 x 32) and a PSF, through ``petakit5d.matlab.XR_RLdeconFrame3D``; the files written are
compared (names, classes, voxels, info files and masks) for the option matrix of DeconFrameTest:
the default result folder and another one, existing results (kept, or replaced with Overwrite),
psfGen with a generated PSF of an earlier file, the eroded mask of the first time point (made,
replaced) and the erosion in the decon, flipZstack, an all-zero frame, several files, a missing
file, OMW, the rotation without deskew (the wrapper's deconRotate), deskew and rotation, zarr
and single results, the default background, the method's case, a GPU job without a usable GPU,
and MATLAB's errors. Each runs on both CPU FFT backends (``.exact``: MATLAB's bits;
``.permissive``: decision record 0006).

The pytest functions port DeconFrameTest's methods in intent: the GPU job's read-ahead, its CPU
rerun for files that don't fit on the GPU, and the large-file switch run with stand-ins (a GPU
with room, an RLdecon that runs on the CPU or runs out of GPU memory, the large-file functions).
"""

from __future__ import annotations

import math
import os
import warnings

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.decon import file as F
from petakit5d.decon import frame as FR
from petakit5d.decon import rl
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, InPath, Tree, Work
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
FSN2 = FSN.replace("_000x_", "_001x_")
FSN3 = FSN.replace("_000x_", "_002x_")

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
    nz, _, nx = a.shape
    z = np.arange(nz, dtype=np.float64).reshape(-1, 1, 1)
    x = np.arange(nx, dtype=np.float64).reshape(1, 1, -1)
    keep = (x > 0.5 * z + 2) & (x < 0.5 * z + nx - 14)
    return np.where(keep, a, np.uint16(0))


def make_root(
    *, frames=((FSN, 1),), data=None, existing=(), psfgen=False, old_mask=False, noisy=False
):
    """data/<frames>.tif, psf/PSF_CamA.tif; existing: result files there first (relative paths);
    psfgen: a generated PSF of an earlier file in data/matlab_decon/psfgen; old_mask: an earlier
    eroded mask (MATLAB's layout) in data/matlab_decon/Masks"""

    def make(folder):
        folder = os.fspath(folder)
        os.makedirs(os.path.join(folder, "data"))
        for name, seed in frames:
            im = frame(seed) if data is None else data(seed)
            write_tiff(os.path.join(folder, "data", name + ".tif"), im)
        os.makedirs(os.path.join(folder, "psf"))
        write_tiff(os.path.join(folder, "psf", "PSF_CamA.tif"), psf_noisy() if noisy else psf_raw())
        for rel in existing:
            p = os.path.join(folder, *rel.split("/"))
            os.makedirs(os.path.dirname(p), exist_ok=True)
            write_tiff(p, frame(9))
        if psfgen:
            p = os.path.join(folder, "data", "matlab_decon", "psfgen")
            os.makedirs(p)
            write_tiff(
                os.path.join(p, "PSF_CamA_simplified.tif"),
                (psf_raw()[2:-2, 3:-3, 1:-1] // 2).astype(np.float32),
            )
        if old_mask:
            p = os.path.join(folder, "data", "matlab_decon", "Masks")
            os.makedirs(p)
            m = np.zeros(IM, np.uint8)
            m[4:-4, 4:-4, 4:-4] = 1
            write_zarr_matlab(
                os.path.join(p, FSN + "_eroded.zarr"), m.transpose(1, 2, 0), (16, 16, 16)
            )

    return make


_CODE = (
    "d = in{1}; f = in{2}; if iscell(f), f = cellfun(@(x) [d, 'data/', x], f, "
    "'UniformOutput', false); else, f = [d, 'data/', f]; end; a = in(6 : end); "
    "for i = find(strcmp(a, 'PSFfile')), if ~isempty(a{i + 1}), a{i + 1} = [d, a{i + 1}]; end; "
    "end; "
    "XR_RLdeconFrame3D(f, in{3}, in{4}, in{5}, a{:}); out = {};"
)


def p_frame(d, files, xy, dz, decon_path, *opts, backend):
    from petakit5d import matlab

    if isinstance(files, (list, tuple)):
        f = [d + "data/" + x for x in files]
    else:
        f = d + "data/" + files
    a = list(opts)
    for i, v in enumerate(opts):
        if v == "PSFfile" and a[i + 1]:
            a[i + 1] = d + a[i + 1]
    rl.clear_otf_cache()
    with P.using(backend):
        matlab.XR_RLdeconFrame3D(f, xy, dz, decon_path, *a)
    return []


BASE = [
    "PSFfile",
    "psf/PSF_CamA.tif",
    "dzPSF",
    DZ,
    "Background",
    BG,
    "RLMethod",
    "simplified",
    "EdgeErosion",
    0,
    "parseCluster",
    False,
    "GPUJob",
    False,
    "DeconIter",
    3,
]


def frame_tree(single=False):
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


def frame_cases(
    cid,
    opts,
    *,
    root=None,
    files=FSN + ".tif",
    decon_path="",
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
                args=[InPath("root"), files, XY, DZ, decon_path, *base, *opts],
                nargout=0,
                python=p_frame,
                py_kwargs={"backend": be},
                inputs={"root": Tree(root or make_root(), "root")},
                expect_error=expect_error,
                py_error=Exception,
                tree=frame_tree(single) if be == "permissive" else tree(ignore=[FIGURE]),
            )
        )
    return res


_DEC = "data/matlab_decon/"
CASES: list[Case] = []
# DeconFrameTest/writesResultsToMatlabDeconByDefault
CASES += frame_cases("frame.default", [])
CASES += frame_cases("frame.decon_path", [], decon_path=Work("decon/"))
# existingResultIsSkipped, overwriteRecomputesTheResult
CASES += frame_cases("frame.existing", [], root=make_root(existing=(_DEC + FSN + ".tif",)))
CASES += frame_cases(
    "frame.overwrite",
    ["Overwrite", True, "Rotate", True],
    root=make_root(existing=(_DEC + FSN + ".tif", _DEC + "DSR/" + FSN + ".tif")),
)
# psfGenFalseUsesTheRawPSF (a generated PSF of an earlier file in psfgen/)
CASES += frame_cases("frame.psfgen_false", ["psfGen", False], root=make_root(psfgen=True))
CASES += frame_cases("frame.psfgen_cached", [], root=make_root(psfgen=True))
# savedErosionMaskCanBeReplaced, and the erosion in the decon
_ZB = make_root(data=lambda s: zero_border(frame(s)))
CASES += frame_cases("frame.mask", ["EdgeErosion", 2, "SaveMaskfile", True], root=_ZB)
CASES += frame_cases(
    "frame.mask_replaced",
    ["EdgeErosion", 2, "SaveMaskfile", True],
    root=make_root(data=lambda s: zero_border(frame(s)), old_mask=True),
)
CASES += frame_cases("frame.erosion", ["EdgeErosion", 3], root=_ZB)
CASES += frame_cases(
    "frame.erode_maskfile_two",
    ["EdgeErosion", 2, "SaveMaskfile", True],
    root=make_root(frames=((FSN, 1), (FSN2, 2)), data=lambda s: zero_border(frame(s))),
    files=[FSN + ".tif", FSN2 + ".tif"],
)
# flipZstackFlipsTheResult, allZeroFrameIsDeconvolved, simplifiedInfoFilesGoToTheDeconInfoFolder
CASES += frame_cases("frame.flip", ["flipZstack", True])
CASES += frame_cases("frame.zero", [], root=make_root(data=lambda s: np.zeros(IM, np.uint16)))
# several files, a missing one
_THREE = make_root(frames=((FSN, 1), (FSN2, 2), (FSN3, 3)))
CASES += frame_cases(
    "frame.three",
    ["DeconIter", 2],
    root=_THREE,
    files=[FSN + ".tif", FSN2 + ".tif", FSN3 + ".tif"],
)
CASES += frame_cases("frame.missing", [], files=["missing.tif", FSN + ".tif"])
# the methods and the options passed on
# OMW: the symmetric PSF only with the exact backend (see psf_noisy), the noisy one on both
_OMW = ["RLMethod", "omw", "skewed", False, "DeconIter", 2]
CASES += frame_cases("frame.omw", _OMW, backends=("exact",))
CASES += frame_cases("frame.omw_noisy", _OMW, root=make_root(noisy=True))
CASES += frame_cases(
    "frame.omw_noisy_two",
    _OMW,
    root=make_root(frames=((FSN, 1), (FSN2, 2)), noisy=True),
    files=[FSN + ".tif", FSN2 + ".tif"],
)
CASES += frame_cases("frame.method_case", ["RLMethod", "Simplified"])
CASES += frame_cases("frame.background_default", ["Background", []])
CASES += frame_cases("frame.single", ["save16bit", False], single=True)
CASES += frame_cases("frame.zarr", ["saveZarr", True, "blockSize", [16, 32, 8]])
CASES += frame_cases("frame.mip_all", ["mipAxis", [1, 1, 1]])
CASES += frame_cases(
    "frame.options",
    ["dampFactor", 1.2, "scaleFactor", 2, "deconOffset", 5, "fixIter", True, "errThresh", 0.01],
)
CASES += frame_cases("frame.debug", ["debug", True, "saveStep", 2])
# rotationWithoutDeskewIsTheCombinedDeskewAndRotation, and deskew with rotation
CASES += frame_cases("frame.rotate", ["Rotate", True, "DeconIter", 2])
CASES += frame_cases("frame.deskew_rotate", ["Deskew", True, "Rotate", True, "Reverse", False])
CASES += frame_cases("frame.save3d_dsr_only", ["Rotate", True, "save3Dstack", [False, True, True]])
CASES += frame_cases(
    "frame.dsr_exists",
    ["Rotate", True, "save3Dstack", [False, True, True]],
    root=make_root(existing=(_DEC + "DSR/" + FSN + ".tif",)),
)
# GPU jobs of two files on an A100 (the next file read during the decon; R04's GPU path gives
# MATLAB's gpuArray bits there)
CASES += frame_cases(
    "gpu.frame.simplified_two",
    ["GPUJob", True, "DeconIter", 3],
    root=make_root(frames=((FSN, 1), (FSN2, 2))),
    files=[FSN + ".tif", FSN2 + ".tif"],
    device="gpu:A100",
)
CASES += frame_cases(
    "gpu.frame.omw_noisy_two",
    ["GPUJob", True, *_OMW],
    root=make_root(frames=((FSN, 1), (FSN2, 2)), noisy=True),
    files=[FSN + ".tif", FSN2 + ".tif"],
    device="gpu:A100",
)
# errors
CASES += frame_cases("error.no_psf", ["PSFfile", ""], base=BASE[2:], expect_error=True)
CASES += frame_cases(
    "error.nothing_to_save",
    ["save3Dstack", [False, False, False], "mipAxis", [0, 0, 0]],
    expect_error=True,
)
CASES += frame_cases(
    "error.gpu_job_without_gpu", ["GPUJob", True], expect_error="decon:noUsableGPU"
)


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
# DeconFrameTest in intent
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def root(tmp_path):
    make_root()(tmp_path)
    return tmp_path


def _fn(root, name=FSN):
    return str(root / "data" / (name + ".tif"))


def run(fn, decon_path, root, **kw):
    base = {
        "psf_file": str(root / "psf" / "PSF_CamA.tif"),
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
        FR.decon_file(fn, XY, DZ, decon_path, **base)


def _read(p):
    from petakit5d.io import read_tiff

    return read_tiff(str(p))


def _norm(psf):
    from petakit5d.kernels import decon as K

    p = np.asarray(psf, np.float32)
    return (p / np.float64(K.sum_double(p))).astype(np.float32)


def _expected(fn, psf, n=3):
    rl.clear_otf_cache()
    return rl.decon_lucy(_read(fn), _norm(psf), n, background=BG, use_gpu=False).decon


def test_writes_results_to_matlab_decon_by_default(root):
    fn = _fn(root)
    run(fn, "", root, decon_iter=3)
    d = root / "data" / "matlab_decon"
    out = _read(d / (FSN + ".tif"))
    assert out.dtype == np.uint16
    np.testing.assert_array_equal(_read(d / "MIPs" / (FSN + "_MIP_z.tif")), out.max(axis=0))
    psf = _read(d / "psfgen" / "PSF_CamA_simplified.tif")
    np.testing.assert_array_equal(out, _expected(fn, psf))


def test_existing_result_is_skipped(root):
    fn = _fn(root)
    run(fn, str(root / "decon"), root, decon_iter=3)
    before = _read(root / "decon" / (FSN + ".tif"))
    write_tiff(fn, frame(7))
    run(fn, str(root / "decon"), root, decon_iter=3)
    np.testing.assert_array_equal(_read(root / "decon" / (FSN + ".tif")), before)


def test_psf_gen_false_uses_the_raw_psf(tmp_path):
    raw = (psf_raw().astype(np.int32) + BG).astype(np.uint16)  # background in every slice
    make_root()(tmp_path)
    write_tiff(tmp_path / "psf" / "PSF_CamA.tif", raw)
    fn1 = _fn(tmp_path)
    fn2 = str(tmp_path / "data" / ("second_" + FSN + ".tif"))
    write_tiff(fn2, _read(fn1))
    d = tmp_path / "decon"
    run(fn1, str(d), tmp_path, decon_iter=3, psf_gen=True)
    run(fn2, str(d), tmp_path, decon_iter=3, psf_gen=False)
    out2 = _read(d / ("second_" + FSN + ".tif"))
    np.testing.assert_array_equal(out2, _expected(fn2, raw))
    assert not np.array_equal(out2, _read(d / (FSN + ".tif")))


def test_saved_erosion_mask_can_be_replaced(tmp_path):
    from petakit5d.decon.erosion import decon_mask_edge_erosion
    from petakit5d.io import read_zarr

    make_root(data=lambda s: zero_border(frame(s)))(tmp_path)
    fn = _fn(tmp_path)
    d = tmp_path / "decon"
    kw = {"decon_iter": 2, "edge_erosion": 2, "save_maskfile": True}
    run(fn, str(d), tmp_path, **kw)
    mask = d / "Masks" / (FSN + "_eroded.zarr")
    assert mask.is_dir()
    os.remove(d / (FSN + ".tif"))
    run(fn, str(d), tmp_path, **kw)
    assert (d / (FSN + ".tif")).is_file()
    np.testing.assert_array_equal(
        read_zarr(str(mask)) > 0, decon_mask_edge_erosion(_read(fn) > 0, 2)
    )
    assert sorted(os.listdir(d / "Masks")) == [FSN + "_eroded.zarr"]


def test_overwrite_recomputes_the_result(root):
    fn = _fn(root)
    d = root / "decon"
    run(fn, str(d), root, decon_iter=3)
    write_tiff(fn, frame(7))
    run(fn, str(d), root, decon_iter=3, overwrite=True)
    psf = _read(d / "psfgen" / "PSF_CamA_simplified.tif")
    np.testing.assert_array_equal(_read(d / (FSN + ".tif")), _expected(fn, psf))


def test_overwrite_removes_the_dsr_result_too(root):
    """a zarr DSR result that is there would be kept, so Overwrite removes it with the decon"""
    from petakit5d.io import read_zarr

    fn = _fn(root)
    d = root / "decon"
    kw = {"decon_iter": 2, "rotate": True, "save_zarr": True}
    run(fn, str(d), root, **kw)
    dsr = d / "DSR" / (FSN + ".zarr")
    before = read_zarr(str(dsr))
    write_tiff(fn, frame(7))
    run(fn, str(d), root, overwrite=True, **kw)
    assert not np.array_equal(read_zarr(str(dsr)), before)


def test_flip_zstack_flips_the_result(root):
    fn = _fn(root)
    run(fn, str(root / "noflip"), root, decon_iter=3, flip_zstack=False)
    run(fn, str(root / "flip"), root, decon_iter=3, flip_zstack=True)
    a = _read(root / "flip" / (FSN + ".tif")).astype(np.float64)
    b = _read(root / "noflip" / (FSN + ".tif"))[::-1].astype(np.float64)
    # the PSF is symmetric in z, so the decon commutes with the flip
    assert np.linalg.norm(a - b) <= 1e-3 * np.linalg.norm(b)


def test_all_zero_frame_is_deconvolved(tmp_path):
    make_root(data=lambda s: np.zeros(IM, np.uint16))(tmp_path)
    run(_fn(tmp_path), str(tmp_path / "decon"), tmp_path, decon_iter=3)
    np.testing.assert_array_equal(
        _read(tmp_path / "decon" / (FSN + ".tif")), np.zeros(IM, np.uint16)
    )
    assert not (tmp_path / "decon" / "decon_info").exists()


def test_simplified_info_files_go_to_the_decon_info_folder(root):
    from petakit5d.io.matfile import read_mat

    d = root / "decon"
    run(_fn(root), str(d), root, decon_iter=3)
    info = read_mat(str(d / "decon_info" / (FSN + "_info.mat")))
    assert sorted(info) == ["debug", "errThresh", "err_mat", "fixIter", "nIter"]
    assert info["err_mat"].shape == (3, 4)
    txt = (d / "decon_info" / (FSN + "_actual_iterations_3.txt")).read_text()
    assert txt == "   3.0000000e+00\n"
    assert not (d / (FSN + "_info.mat")).exists()


class _Calls:
    def __init__(self):
        self.use_gpu = []


def _fake(monkeypatch, err_id=None, run_cpu=False):
    """an RLdecon that records use_gpu, raises err_id on the GPU (or runs the real one on the
    CPU)"""
    calls = _Calls()
    real = F.rl_decon_file

    def fake(*args, **kw):
        calls.use_gpu.append(kw.get("use_gpu"))
        if run_cpu:
            kw["use_gpu"] = False
            return real(*args, **kw)
        if kw.get("use_gpu") and err_id is not None:
            from petakit5d.kernels.decon_gpu import GPUOutOfMemoryError

            if err_id.startswith("parallel:gpu:array"):
                raise GPUOutOfMemoryError(err_id, "No GPU for this test.")
            raise rl.DeconError(err_id, "No GPU for this test.")
        from petakit5d.io import write_tiff as w

        w(np.zeros((4, 4, 4), np.uint16), args[1])
        return F.RLDeconResult(None, None, None, None)

    monkeypatch.setattr(F, "rl_decon_file", fake)
    return calls


@pytest.mark.parametrize("err_id", FR.GPU_CAPACITY_IDS)
def test_file_that_does_not_fit_on_the_gpu_is_deconvolved_on_the_cpu(root, monkeypatch, err_id):
    calls = _fake(monkeypatch, err_id)
    run(_fn(root), str(root / "decon"), root, decon_iter=3, gpu_job=True, overwrite=True)
    assert calls.use_gpu == [True, False]
    assert (root / "decon" / (FSN + ".tif")).is_file()


def test_other_gpu_job_errors_are_not_retried_on_the_cpu(root, monkeypatch):
    calls = _fake(monkeypatch, "MATLAB:test:otherError")
    with pytest.raises(rl.DeconError) as e:
        run(_fn(root), str(root / "decon"), root, gpu_job=True)
    assert e.value.identifier == "MATLAB:test:otherError"
    assert calls.use_gpu == [True]


def _spy(monkeypatch):
    """a log of reads ('read <file>') and of renames onto results ('move <destination>')"""
    from petakit5d import io

    events = []
    read, rename = io.read_tiff, F._move

    def read_spy(path, *a, **kw):
        events.append(f"read {os.path.normpath(os.fspath(path))}")
        return read(path, *a, **kw)

    def rename_spy(src, dst):
        events.append(f"move {os.path.normpath(os.fspath(dst))}")
        return rename(src, dst)

    monkeypatch.setattr(io, "read_tiff", read_spy)
    monkeypatch.setattr(F, "_move", rename_spy)
    return events


def _three(tmp_path):
    make_root(frames=((FSN, 1), (FSN2, 2), (FSN3, 3)))(tmp_path)
    return [_fn(tmp_path, n) for n in (FSN, FSN2, FSN3)]


def test_gpu_job_reads_the_next_file_during_the_decon(tmp_path, monkeypatch):
    """each file is read once, the next one before the current result is published, with the
    results of one file per job (a GPU with room; RLdecon on the CPU)"""
    monkeypatch.setattr(FR, "gpu_free_bytes", lambda: 10**15)
    calls = _fake(monkeypatch, run_cpu=True)
    fns = _three(tmp_path)
    kw = {"decon_iter": 2, "rl_method": "simplified", "gpu_job": True}
    for f in fns:
        run(f, str(tmp_path / "single"), tmp_path, **kw)
    events = _spy(monkeypatch)
    run(fns, str(tmp_path / "multi"), tmp_path, **kw)
    for f in fns:
        name = os.path.basename(f)
        np.testing.assert_array_equal(
            _read(tmp_path / "multi" / name), _read(tmp_path / "single" / name)
        )
        assert events.count(f"read {os.path.normpath(f)}") == 1
    for i in range(2):
        nxt = events.index(f"read {os.path.normpath(fns[i + 1])}")
        cur = events.index(
            f"move {os.path.normpath(str(tmp_path / 'multi' / os.path.basename(fns[i])))}"
        )
        assert nxt < cur
    assert all(calls.use_gpu)


def test_gpu_job_without_room_on_the_gpu_reads_the_files_in_turn(tmp_path, monkeypatch):
    monkeypatch.setattr(FR, "gpu_free_bytes", lambda: 0)
    _fake(monkeypatch, run_cpu=True)
    fns = _three(tmp_path)[:2]
    events = _spy(monkeypatch)
    run(fns, str(tmp_path / "decon"), tmp_path, decon_iter=2, gpu_job=True)
    r2 = f"read {os.path.normpath(fns[1])}"
    assert events.count(r2) == 1
    first = os.path.normpath(str(tmp_path / "decon" / (FSN + ".tif")))
    assert events.index(r2) > events.index(f"move {first}")


def test_gpu_job_does_not_read_a_finished_next_file(tmp_path, monkeypatch):
    monkeypatch.setattr(FR, "gpu_free_bytes", lambda: 10**15)
    _fake(monkeypatch, run_cpu=True)
    fns = _three(tmp_path)[:2]
    d = str(tmp_path / "decon")
    run(fns[1], d, tmp_path, decon_iter=2, gpu_job=True)
    events = _spy(monkeypatch)
    run(fns, d, tmp_path, decon_iter=2, gpu_job=True)
    assert events.count(f"read {os.path.normpath(fns[1])}") == 0
    assert (tmp_path / "decon" / (FSN + ".tif")).is_file()


def test_gpu_has_room_is_twice_the_decon(monkeypatch):
    monkeypatch.setattr(FR, "gpu_free_bytes", lambda: 2 * 42 * 1000)
    assert FR.gpu_has_room((10, 10, 10), "simplified")
    assert not FR.gpu_has_room((10, 10, 11), "simplified")
    monkeypatch.setattr(FR, "gpu_free_bytes", lambda: 2 * 38 * 1000)
    assert FR.gpu_has_room((10, 10, 10), "omw")
    monkeypatch.setattr(FR, "gpu_free_bytes", lambda: None)
    assert not FR.gpu_has_room((1, 1, 1), "omw")


def test_without_a_gpu_there_is_no_room():
    assert FR.gpu_free_bytes() is None or isinstance(FR.gpu_free_bytes(), int)


def test_rotation_without_deskew_is_the_combined_deskew_and_rotation(root):
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.kernels import ops

    run(_fn(root), "", root, decon_iter=2, rotate=True)
    d = root / "data" / "matlab_decon"
    decon = _read(d / (FSN + ".tif"))
    exp = ops.convert(
        deskew_rotate_volume(decon, 32.45, DZ, XY, reverse=True, crop=True, save_16bit=True),
        "uint16",
    )
    np.testing.assert_array_equal(_read(d / "DSR" / (FSN + ".tif")), exp)
    np.testing.assert_array_equal(_read(d / "DSR" / "MIPs" / (FSN + "_MIP_z.tif")), exp.max(axis=0))


def test_large_files_go_to_the_large_file_functions(root, monkeypatch):
    """a file above CPUMaxMem, and every file after it, goes to R07's in-memory function (and the
    call returns after it, as MATLAB's)"""
    from petakit5d import decon

    seen = []
    monkeypatch.setattr(decon, "decon_large_in_memory", lambda **kw: seen.append(kw), raising=False)
    run(
        [_fn(root), _fn(root)], str(root / "decon"), root, cpu_max_mem=1e-6, batch_size=(16, 32, 32)
    )
    assert len(seen) == 1
    kw = seen[0]
    assert kw["frame_fullpaths"] == _fn(root)
    assert kw["decon_fullpath"].endswith(FSN + ".tif")
    assert kw["batch_size"] == (16, 32, 32) and kw["background"] == BG
    assert not (root / "decon" / (FSN + ".tif")).exists()


def test_without_the_large_file_port_a_large_file_raises(root, monkeypatch):
    from petakit5d import decon

    monkeypatch.delattr(decon, "decon_large_in_memory", raising=False)
    with pytest.raises(NotImplementedError, match="R07"):
        run(_fn(root), str(root / "decon"), root, large_file=True)


def test_matlab_errors(root):
    with pytest.raises(FR.DeconFileError):
        run(_fn(root), "", root, psf_file="")
    with pytest.raises(FR.DeconFileError, match="needs to be saved"):
        run(_fn(root), "", root, save_3d_stack=(False, False, False), mip_axis=(0, 0, 0))


def test_a_missing_file_is_skipped_with_a_warning(root):
    with pytest.warns(UserWarning, match="does not exist"):
        FR.decon_file(
            [str(root / "nope.tif"), _fn(root)],
            XY,
            DZ,
            str(root / "decon"),
            psf_file=str(root / "psf" / "PSF_CamA.tif"),
            dz_psf=DZ,
            background=BG,
            edge_erosion=0,
            parse_cluster=False,
            decon_iter=1,
        )
    assert (root / "decon" / (FSN + ".tif")).is_file()
