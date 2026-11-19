"""Camera background estimation and the phase-image flat field (petakit5d.flatfield).

Oracle cases (MATLAB: XR_flatfield_background_estimation, XR_flatfield_background_wrapper and
XR_phase_image_flat_field_correction): the dark stacks of FlatfieldBackgroundEstimationTest
(noise that grows with z, bright slices), frames concatenated along z, uint16, single and zarr
frames, more than 200 slices (the percentiles), the options, a stack whose slices are all
excluded, a result in the current folder (M10's fix), MATLAB's error for a char list; the
wrapper with two cameras, two time points and channels, names without a camera, MIPs, partial
volume files, zarr, two datasets, the 'cam_ch' and 'came_ch_iter' modes and the modes MATLAB
fails on. Inputs are made in Python's form and handed to MATLAB as files.

The phase flat field: the pre-filter (``imgaussfilt`` sigma 100, frequency domain: within the
filters' tolerance) and the BaSiC flat field, which BaSiCPy estimates as another implementation
of the method (the tolerance of docs/parity/flatfield.md, approved). Its cases need BaSiCPy.
Then FlatfieldBackgroundEstimationTest's and PhaseFlatFieldCorrectionTest's methods.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.core.arrays import to_matlab
from petakit5d.flatfield import basic as fb
from petakit5d.flatfield import estimation as fe
from petakit5d.io import read_mat, read_tiff
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, File, In, InPath, Tree, Value, Work
from petakit5d.testing.compare import tol, tree
from petakit5d.testing.synth import CounterRNG

FSN = (
    "Scan_Iter_{it:04d}_0000_{cam}ch{ch}_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_"
    "000x_000y_000z_0000t"
)
DARK = 100

# ---------------------------------------------------------------------------------------------
# inputs (Python's form: stacks (z, y, x))
# ---------------------------------------------------------------------------------------------


def _rng(name):
    return CounterRNG(7, f"flatfield_b07/{name}")


def pattern(shape=(48, 40), kind="flat"):
    """a (y, x) dark pattern: flat at 100, or 100 + 0..20 varying over the image"""
    if kind == "flat":
        return np.full(shape, float(DARK))
    y, x = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing="ij")
    return (DARK + (3 * y + 5 * x) % 21).astype(np.float64)


def dark_stack(name, nz=40, shape=(48, 40), kind="flat", lo=2.0, hi=6.0, cls="uint16"):
    """FlatfieldBackgroundEstimationTest's darkStack: nz slices of the pattern plus noise whose
    sigma grows from lo to hi with z"""
    sig = lo + (hi - lo) * np.arange(nz, dtype=np.float64) / max(nz - 1, 1)
    v = pattern(shape, kind)[None] + sig[:, None, None] * _rng(name).normal((nz, *shape))
    return np.asarray(dtypes.cast(v, cls))


def bright_stack(name, nz=40):
    """darkStack plus 2000 in a square of the first 6 slices (brightSlicesAreExcluded)"""
    st = np.asarray(dark_stack(name, nz), dtype=np.float32)
    st[:6, 9:30, 9:30] += 2000
    return np.asarray(dtypes.cast(st, "uint16"))


def spiked_stack(name, nz, shape, lo, hi, slices, npix, amp):
    """darkStack with npix pixels (every 97th) raised by amp in the given slices: few bright
    pixels, so the slices stay below the median noise and reach the removal loop"""
    st = np.asarray(dark_stack(name, nz, shape, lo=lo, hi=hi), dtype=np.float64)
    flat = st.reshape(nz, -1)
    for z in slices:
        flat[z, ::97][:npix] += amp
    return np.asarray(dtypes.cast(st, "uint16"))


def spotted_stack(name, nz, shape, lo, hi, spots):
    """darkStack with spots {z: [(npix, amp), ...]}: npix pixels raised by amp in slice z, at
    places that differ from slice to slice and group to group (so no pixel's sigma over z is
    raised much)"""
    st = np.asarray(dark_stack(name, nz, shape, lo=lo, hi=hi), dtype=np.float64)
    flat = st.reshape(nz, -1)
    for z, groups in spots.items():
        for g, (npix, amp) in enumerate(groups):
            flat[z, (7 * z + 31 * g) % 97 :: 97][:npix] += amp
    return np.asarray(dtypes.cast(st, "uint16"))


def tif_file(make):
    def write(path):
        synth.write_tiff(path, make())

    return write


def tif_tree(files):
    """a folder of tiff files: {name: maker}"""

    def write(folder):
        for name, make in files.items():
            synth.write_tiff(Path(folder) / name, make())

    return write


def zarr_tree(files):
    def write(folder):
        for name, make in files.items():
            synth.write_zarr(Path(folder) / name, make(), block_size=(8, 16, 16), layout="matlab")

    return write


def fsn(cam="CamA_", ch=0, it=0, ext=".tif"):
    return FSN.format(it=it, cam=cam, ch=ch) + ext


# ---------------------------------------------------------------------------------------------
# the Python sides (the MATLAB-name layer, as MATLAB is called)
# ---------------------------------------------------------------------------------------------


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def py_estimation(frames, result, *nv):
    return call_matlab("XR_flatfield_background_estimation", list(frames), result, *nv)


def py_estimation_mips(frames, result, *nv):
    os.makedirs(os.path.join(os.path.dirname(result) or ".", "MIPs"), exist_ok=True)
    return py_estimation(frames, result, *nv)


def py_wrapper(*args):
    return call_matlab("XR_flatfield_background_wrapper", *args)


def py_prefilter(im, const):
    return fb.phase_prefilter(im, const)


def py_phase(d, const):
    names = sorted(n for n in os.listdir(d) if n.endswith(".tif"))
    out = os.path.join(os.getcwd(), "ff") + os.sep
    os.makedirs(out, exist_ok=True)
    call_matlab("XR_phase_image_flat_field_correction", [d + n for n in names], out, const)
    a = read_mat(f"{out}BaSiC_flatfield/{names[0][:-4]}_flatfield.mat")
    tiles = np.stack([read_tiff(out + n) for n in names])
    return [a["flatfield"], a["minval"], to_matlab(tiles)]


EST = "XR_flatfield_background_estimation"
MIPS = "mkdir([work, 'MIPs']); "


def est_case(cid, make, *nv, mips=False, name="im.tif", **kw):
    code = None
    if mips:
        code = MIPS + f"{EST}(in{{1}}, in{{2}}, in{{3 : end}}); out = {{}};"
    return Case(
        cid,
        matlab=None if mips else EST,
        code=code,
        nargout=0,
        args=[[InPath("im")], Work("bg.mat"), *nv],
        python=py_estimation_mips if mips else py_estimation,
        inputs={"im": File(make, name)},
        **kw,
    )


def est_tree_case(cid, make, names, *nv):
    files = ", ".join(f"[in{{1}}, '{n}']" for n in names)
    return Case(
        cid,
        code=f"{EST}({{{files}}}, in{{2 : end}}); out = {{}};",
        nargout=0,
        args=[InPath("d"), Work("bg.mat"), *nv],
        python=lambda d, r, *a: py_estimation([d + n for n in names], r, *a),
        inputs={"d": Tree(make, "d")},
    )


def py_wrapper_twice(*args):
    py_wrapper(*args)
    return py_wrapper(*args)


def wrapper_case(cid, make, *nv, error=False, folders=("d",), twice=False):
    args = [InPath(folders[0])] if len(folders) == 1 else [[InPath(f) for f in folders]]
    code = "XR_flatfield_background_wrapper(in{:}); " * (2 if twice else 1) + "out = {};"
    return Case(
        cid,
        code=code,
        nargout=0,
        args=[*args, *nv, "parseCluster", False, "masterCompute", True],
        python=py_wrapper_twice if twice else py_wrapper,
        inputs={f: Tree(m, f) for f, m in zip(folders, make, strict=True)},
        expect_error=error,
        py_error=(fe.MatlabError,),
        tree=tree(ignore_vars=()),
    )


def _cases_estimation():
    wavy = lambda: dark_stack("wavy", kind="wavy")  # noqa: E731
    return [
        est_case("est.dark_flat", tif_file(lambda: dark_stack("flat"))),
        est_case("est.dark_wavy", tif_file(wavy)),
        est_case("est.bright_slices", tif_file(lambda: bright_stack("bright"))),
        est_case("est.single_frames", tif_file(lambda: dark_stack("single", cls="single"))),
        est_case("est.mip", tif_file(wavy), "saveMIP", True, mips=True),
        est_case(
            "est.options",
            tif_file(lambda: bright_stack("opts")),
            "sigmaThresh",
            4,
            "sigmaFactor",
            2.5,
            "probThresh",
            0.01,
        ),
        est_case(
            "est.percentiles",
            tif_file(lambda: dark_stack("many", nz=230, shape=(20, 28), kind="wavy")),
        ),
        # the removal loop: one quiet slice with 0.12% bright pixels among ~220 kept (removed: a
        # step that removes under 0.1% stops the loop only from 1,000 slices)
        est_case(
            "est.one_spiked_slice",
            tif_file(lambda: spiked_stack("spk1", 440, (64, 64), 2.0, 6.0, [3], 5, 70)),
        ),
        # noisy slices (sigma above sigmaThresh) with 0.5% bright pixels are removed
        est_case(
            "est.noisy_spiked_slices",
            tif_file(
                lambda: spiked_stack("spk2", 60, (48, 40), 6.0, 20.0, range(0, 60, 3), 10, 60)
            ),
        ),
        # one quiet slice of ~220 with 6 bright pixels and one far out: a step that removes
        # one slice (0.45%) goes on
        est_case(
            "est.single_removal",
            tif_file(lambda: spotted_stack("sr", 440, (64, 64), 2.0, 2.6, {0: [(6, 11), (1, 50)]})),
        ),
        # the same with the far pixel at 15 sigma: below the range limit (20 sigma), kept
        est_case(
            "est.range_limit",
            tif_file(lambda: spotted_stack("rl", 440, (64, 64), 2.0, 2.6, {0: [(6, 11), (1, 30)]})),
        ),
        # slices with noise between sigmaThresh and the median and 0.16% bright pixels go
        est_case(
            "est.noise_limit",
            tif_file(
                lambda: spotted_stack(
                    "nl", 60, (48, 40), 6.0, 20.0, {z: [(3, 120)] for z in range(2, 7)}
                )
            ),
        ),
        # a median noise above sigmaThresh * 5: that is the first threshold
        est_case("est.very_noisy", tif_file(lambda: dark_stack("vn", lo=30.0, hi=42.0))),
        est_case("est.prob_thresh_zero", tif_file(lambda: bright_stack("pz")), "probThresh", 0),
        # more than 200 slices kept: the percentiles
        est_case(
            "est.percentiles_440",
            tif_file(lambda: dark_stack("many440", nz=440, shape=(20, 28), lo=3.0, hi=4.0)),
        ),
        est_tree_case(
            "est.planes",
            tif_tree(
                {f"p{k:02d}.tif": lambda k=k: dark_stack(f"pl{k}", nz=1)[0] for k in range(30)}
            ),
            [f"p{k:02d}.tif" for k in range(30)],
        ),
        est_case(
            "est.constant_stack",
            tif_file(lambda: np.full((12, 20, 16), 100, np.uint16)),
        ),
        est_case(
            "est.quiet_noise",
            tif_file(lambda: dark_stack("quiet", lo=0.5, hi=12.0, shape=(33, 47))),
        ),
        est_tree_case(
            "est.two_frames",
            tif_tree(
                {"a.tif": lambda: dark_stack("p0", nz=25), "b.tif": lambda: bright_stack("p1")}
            ),
            ["a.tif", "b.tif"],
        ),
        Case(
            "est.zarr",
            code=f"{EST}({{[in{{1}}, 'im.zarr']}}, in{{2 : end}}); out = {{}};",
            nargout=0,
            args=[InPath("d"), Work("bg.mat"), "zarrFile", True],
            python=lambda d, r, *a: py_estimation([d + "im.zarr"], r, *a),
            inputs={"d": Tree(zarr_tree({"im.zarr": lambda: dark_stack("zarr")}), "d")},
        ),
        # M10's b4eade6: a name without a folder is in the current folder (the work folder)
        Case(
            "est.current_folder",
            code=MIPS + f"{EST}({{'im.tif'}}, 'bg.mat', 'saveMIP', true); out = {{}};",
            nargout=0,
            python=lambda: py_estimation_mips(["im.tif"], "bg.mat", "saveMIP", True),
            inputs={"im": File(tif_file(wavy), "im.tif")},
            files=True,
        ),
        Case(
            "est.char_list_error",
            code=f"{EST}([in{{1}}], in{{2}}); out = {{}};",
            nargout=0,
            args=[InPath("im"), Work("bg.mat")],
            python=lambda f, r: call_matlab(EST, f, r),
            inputs={"im": File(tif_file(wavy), "im.tif")},
            expect_error="MATLAB:cellRefFromNonCell",
            py_error=(fe.MatlabError,),
        ),
    ]


def two_cams():
    return tif_tree(
        {
            fsn("CamA_"): lambda: dark_stack("camA", kind="wavy"),
            fsn("CamB_"): lambda: dark_stack("camB", kind="wavy"),
        }
    )


def time_points():
    """two time points of two channels of camera A, and camera B's channel 0"""
    files = {}
    for it in range(2):
        for ch in range(2):
            files[fsn("CamA_", ch, it)] = lambda it=it, ch=ch: (
                bright_stack(f"tp{it}{ch}") if (it, ch) == (1, 0) else dark_stack(f"tp{it}{ch}")
            )
        files[fsn("CamB_", 0, it)] = lambda it=it: dark_stack(f"tpB{it}", nz=30)
    return tif_tree(files)


def _cases_wrapper():
    chans = ("channelPatterns", ["CamA_ch0", "CamB_ch0"])
    tp_chans = ("channelPatterns", ["CamA_ch0", "CamA_ch1", "CamB_ch0"])
    no_cam = tif_tree(
        {
            fsn("", 0): lambda: dark_stack("nc0"),
            fsn("", 1): lambda: dark_stack("nc1", nz=30),
        }
    )
    parts = tif_tree(
        {
            fsn("CamA_"): lambda: dark_stack("pv0", nz=24),
            fsn("CamA_")[:-4] + "_part0001.tif": lambda: dark_stack("pv1", nz=20),
        }
    )
    return [
        wrapper_case("wrapper.two_cams", [two_cams()], *chans),
        wrapper_case(
            "wrapper.bright",
            [tif_tree({fsn("CamA_"): lambda: bright_stack("wb")})],
            "channelPatterns",
            ["CamA_ch0"],
        ),
        wrapper_case(
            "wrapper.mcc_false",
            [tif_tree({fsn("CamA_"): lambda: dark_stack("mcc")})],
            "channelPatterns",
            ["CamA_ch0"],
            "mccMode",
            False,
        ),
        wrapper_case("wrapper.time_points_mips", [time_points()], *tp_chans, "saveMIP", True),
        # the second call finds the tasks' results and the backgrounds, and keeps them
        wrapper_case("wrapper.twice", [two_cams()], *chans, twice=True),
        # a prefix before Scan_Iter_ names the backgrounds; a task without slices is left out
        wrapper_case(
            "wrapper.prefix_and_empty_task",
            [
                tif_tree(
                    {
                        "pre_" + fsn("CamA_", 0, 0): lambda: dark_stack("pre0"),
                        "pre_" + fsn("CamA_", 0, 1): lambda: np.full((12, 48, 40), 99, np.uint16),
                    }
                )
            ],
            "channelPatterns",
            ["CamA_ch0"],
        ),
        wrapper_case("wrapper.no_cam", [no_cam], "channelPatterns", ["ch0", "ch1"]),
        wrapper_case("wrapper.partial_volumes", [parts], "channelPatterns", ["CamA_ch0"]),
        wrapper_case(
            "wrapper.zarr",
            [
                zarr_tree(
                    {
                        fsn("CamA_", ext=".zarr"): lambda: dark_stack("wz0"),
                        fsn("CamB_", ext=".zarr"): lambda: bright_stack("wz1"),
                    }
                )
            ],
            *chans,
            "zarrFile",
            True,
        ),
        wrapper_case(
            "wrapper.two_datasets", [two_cams(), time_points()], *chans, folders=("d1", "d2")
        ),
        wrapper_case("wrapper.cam_ch", [time_points()], *tp_chans, "estimationMode", "cam_ch"),
        wrapper_case(
            "wrapper.came_ch_iter",
            [time_points()],
            *tp_chans,
            "estimationMode",
            "came_ch_iter",
        ),
        wrapper_case(
            "wrapper.came_ch_iter_first_tp",
            [time_points()],
            *tp_chans,
            "estimationMode",
            "came_ch_iter",
            "onlyFirstTP",
            True,
        ),
        wrapper_case(
            "wrapper.cam_ch_two_datasets",
            [two_cams(), time_points()],
            *chans,
            "estimationMode",
            "cam_ch",
            folders=("d1", "d2"),
        ),
        wrapper_case(
            "wrapper.cam_ch_iter_stack",
            [two_cams()],
            *chans,
            "estimationMode",
            "cam_ch_iter_stack",
            error="MATLAB:UndefinedFunction",
        ),
        wrapper_case(
            "wrapper.unknown_mode",
            [two_cams()],
            *chans,
            "estimationMode",
            "per_tile",
            error="MATLAB:UndefinedFunction",
        ),
    ]


# ---------------------------------------------------------------------------------------------
# the phase flat field
# ---------------------------------------------------------------------------------------------


def texture(shape, name, lo=500, hi=3000):
    """a smooth random texture in [lo, hi]: uniform noise, box-filtered (integer sums)"""
    ny, nx = shape
    r = 4
    n = _rng(name).integers(0, 1000, (ny + 2 * r, nx + 2 * r))
    c = np.cumsum(np.cumsum(n, axis=0), axis=1)
    c = np.pad(c, ((1, 0), (1, 0)))
    k = 2 * r + 1
    s = c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]
    s = s.astype(np.float64)
    return lo + (hi - lo) * (s - s.min()) / max(s.max() - s.min(), 1.0)


def shading(shape):
    """PhaseFlatFieldCorrectionTest's vignetting: 1 at the center, 0.6 at the corners"""
    ny, nx = shape
    y = (2 * np.arange(ny) - (ny - 1)) / (ny - 1)
    x = (2 * np.arange(nx) - (nx - 1)) / (nx - 1)
    return np.asarray(dtypes.cast(1 - 0.2 * (x[None, :] ** 2 + y[:, None] ** 2), "single"))


def phase_tile(shape, k, name):
    v = texture(shape, f"{name}/{k}") * shading(shape).astype(np.float64)
    return np.asarray(dtypes.cast(v, "uint16"))


def phase_tree(shape, n, name):
    def write(folder):
        for k in range(n):
            f = (
                "Scan_Iter_0000_CamB_ch0_CAM1_stack0000_642nm_0000000msec_0000006528774msecAbs_"
                f"000x_{k:03d}y_000z_0000t.tif"
            )
            synth.write_tiff(Path(folder) / f, phase_tile(shape, k, name))

    return write


# the filters' frequency-domain row of docs/parity/matlabcompat.md; minval is the smallest
# value of the pre-filtered tiles, whose rounding is relative to the largest
PREFILTER_TOL = tol(reason="matlabcompat: imgaussfilt frequency domain", rel_max=1e-6)
MINVAL_TOL = tol(reason="matlabcompat: imgaussfilt frequency domain", abs=0.01)
# docs/parity/flatfield.md (approved 2026-11-18, N8): BaSiCPy's flat field against MATLAB's
# BaSiC, and the tiles divided by it (uint16: plus their rounding): 12 tiles or more (the phase
# demo's 20 measured 5.6% at most, 0.64% in L2); 8 tiles with strong vignetting differ more
BASIC = "flatfield: BaSiCPy flat field (approved, N8)"
FEW_TILES = {"rel": 0.12, "rel_l2": 0.03}
TILES = {"rel": 0.06, "rel_l2": 0.01}


def basic_tols(t):
    return {
        "flatfield": tol(reason=BASIC, **t),
        "minval": MINVAL_TOL,
        "tiles": tol(reason=BASIC, abs=1, **t),
    }


PHASE_CODE = (
    "d = in{1}; s = dir([d, '*.tif']); names = sort({s.name}); "
    "fns = cellfun(@(x) [d, x], names(:), 'unif', 0); outPath = [work, 'ff/']; mkdir(outPath); "
    "XR_phase_image_flat_field_correction(fns, outPath, in{2}); "
    "a = load(sprintf('%sBaSiC_flatfield/%s_flatfield.mat', outPath, names{1}(1 : end - 4))); "
    "t = cellfun(@(x) readtiff([outPath, x]), names, 'unif', 0); "
    "out = {a.flatfield, a.minval, cat(3, t{:})};"
)


def phase_case(cid, shape, n, const, t=TILES):
    return Case(
        cid,
        code=PHASE_CODE,
        nargout=3,
        outputs=("flatfield", "minval", "tiles"),
        args=[InPath("d"), const],
        python=py_phase,
        inputs={"d": Tree(phase_tree(shape, n, cid), "d")},
        files=False,
        compare=basic_tols(t),
    )


def _prefilter_case(cid, shape, const):
    return Case(
        cid,
        code="imf = single(in{1}); out = {imf - imgaussfilt(imf, 100) + in{2}};",
        args=[In("im"), const],
        python=py_prefilter,
        inputs={"im": Value(lambda: phase_tile(shape, 0, cid))},
        compare=PREFILTER_TOL,
    )


PREFILTER_CASES = [
    _prefilter_case("phase.prefilter_128", (128, 128), 1100),
    _prefilter_case("phase.prefilter_96x160", (96, 160), 6000),
]
PHASE_CASES = [
    phase_case("phase.tiles_128x8", (128, 128), 8, 1100, FEW_TILES),
    phase_case("phase.tiles_96x160x12", (96, 160), 12, 6000),
    phase_case("phase.tiles_200x256x20", (200, 256), 20, 6000),
]
CASES = _cases_estimation() + _cases_wrapper()


@oracle.parametrize(CASES + PREFILTER_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


@oracle.parametrize(PHASE_CASES)
def test_matlab_basic(case, tmp_path):
    if not fb.basic_available():
        pytest.skip("needs BaSiCPy (the basic extra)")
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# FlatfieldBackgroundEstimationTest
# ---------------------------------------------------------------------------------------------


def _write(path, vol):
    synth.write_tiff(path, vol)


def _bg(path, cam="A"):
    return read_tiff(f"{path}/matlab_flat_field_estimation/Scan_Cam{cam}_ch0_CAM1_background.tif")


FSN0 = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_"
    "000z_0000t.tif"
)


def test_estimates_dark_level_per_camera(tmp_path):
    pat = pattern(kind="wavy")
    for cam in ("CamA", "CamB"):
        _write(tmp_path / FSN0.format(cam), dark_stack(f"t/{cam}", kind="wavy"))
    fe.flatfield_background_dataset(
        str(tmp_path), channel_patterns=("CamA_ch0", "CamB_ch0"), parse_cluster=False
    )
    for cam in "AB":
        bg = _bg(tmp_path, cam)
        assert bg.shape == (48, 40)
        err = bg.astype(np.float64) - pat
        assert abs(err.mean()) < 0.3
        assert np.abs(err).max() < 5
    a = read_mat(
        f"{tmp_path}/matlab_flat_field_estimation/{FSN0.format('CamA')[:-4]}.mat", ["numSlice"]
    )
    assert 0 < float(a["numSlice"]) < 40, "Noisy or bright slices should be excluded."


def test_bright_slices_are_excluded(tmp_path):
    _write(tmp_path / FSN0.format("CamA"), bright_stack("t/bright"))
    fe.flatfield_background_dataset(
        str(tmp_path), channel_patterns=("CamA_ch0",), parse_cluster=False
    )
    err = _bg(tmp_path).astype(np.float64) - DARK
    assert abs(err.mean()) < 0.3
    assert np.abs(err).max() < 5


def test_accepts_mcc_mode_false(tmp_path):
    from petakit5d import matlab

    _write(tmp_path / FSN0.format("CamA"), dark_stack("t/mcc"))
    matlab.XR_flatfield_background_wrapper(
        str(tmp_path), "channelPatterns", ["CamA_ch0"], "parseCluster", False, "mccMode", False
    )
    assert os.path.isfile(
        f"{tmp_path}/matlab_flat_field_estimation/Scan_CamA_ch0_CAM1_background.tif"
    )


def test_estimate_of_a_file_in_the_current_folder_stays_there(tmp_path, monkeypatch):
    fn = FSN0.format("CamA")
    _write(tmp_path / fn, dark_stack("t/cur"))
    (tmp_path / "full" / "MIPs").mkdir(parents=True)
    (tmp_path / "MIPs").mkdir()
    fe.flatfield_background_estimation(
        [str(tmp_path / fn)], str(tmp_path / "full" / "bg.mat"), save_mip=True
    )
    monkeypatch.chdir(tmp_path)
    fe.flatfield_background_estimation([fn], "bg.mat", save_mip=True)
    a, b = read_mat(tmp_path / "bg.mat"), read_mat(tmp_path / "full" / "bg.mat")
    assert a.keys() == b.keys()
    for k in a:
        np.testing.assert_array_equal(np.asarray(a[k]), np.asarray(b[k]))
    np.testing.assert_array_equal(
        read_tiff(tmp_path / "MIPs" / "bg_MIP_z.tif"),
        read_tiff(tmp_path / "full" / "MIPs" / "bg_MIP_z.tif"),
    )


# ---------------------------------------------------------------------------------------------
# PhaseFlatFieldCorrectionTest
# ---------------------------------------------------------------------------------------------


@pytest.mark.skipif(not fb.basic_available(), reason="needs BaSiCPy (the basic extra)")
def test_corrects_tiles_with_estimated_flat_field(tmp_path):
    from petakit5d.matlabcompat.filters import imgaussfilt

    shape, const = (128, 128), 1100
    out = str(tmp_path / "ff_corrected") + os.sep
    os.makedirs(out)
    fns, raw = [], []
    for i in range(8):
        raw.append(phase_tile(shape, i, "t/phase"))
        fns.append(str(tmp_path / f"tile_{i:03d}.tif"))
        _write(fns[-1], raw[-1])
    fb.phase_image_flatfield_correction(fns, out, const)
    a = read_mat(f"{out}BaSiC_flatfield/tile_000_flatfield.mat")
    ff = np.asarray(a["flatfield"])
    assert ff.shape == shape
    r = np.corrcoef(ff.astype(np.float64).ravel(), shading(shape).astype(np.float64).ravel())
    assert r[0, 1] > 0.8, "The estimated flat field does not follow the shading."
    for i in range(8):
        got = read_tiff(f"{out}tile_{i:03d}.tif")
        imf = raw[i].astype(np.float32)
        pre = imf - imgaussfilt(imf, 100) + np.float32(const)
        expected = np.asarray(dtypes.cast(pre / ff, "uint16"))
        assert got.dtype == np.uint16
        assert np.abs(got.astype(np.int64) - expected.astype(np.int64)).max() <= 1


@pytest.mark.skipif(not fb.basic_available(), reason="needs BaSiCPy (the basic extra)")
def test_phase_tiles_are_the_same_on_one_and_on_several_threads(tmp_path, monkeypatch):
    """the tiles are pre-filtered and divided in threads: the files are those of one thread"""
    import petakit5d.core.threads as th

    fns = []
    for i in range(6):
        fns.append(str(tmp_path / f"tile_{i:03d}.tif"))
        _write(fns[-1], phase_tile((96, 160), i, "t/threads"))
    outs = {}
    for n in (1, 4):
        monkeypatch.setattr(th, "resolve_threads", lambda k=None, n=n: n)
        out = str(tmp_path / f"out{n}") + os.sep
        os.makedirs(out)
        fb.phase_image_flatfield_correction(fns, out, 6000)
        outs[n] = out
    for i in range(6):
        np.testing.assert_array_equal(
            read_tiff(f"{outs[1]}tile_{i:03d}.tif"), read_tiff(f"{outs[4]}tile_{i:03d}.tif")
        )
    a = read_mat(f"{outs[1]}BaSiC_flatfield/tile_000_flatfield.mat")
    b = read_mat(f"{outs[4]}BaSiC_flatfield/tile_000_flatfield.mat")
    np.testing.assert_array_equal(np.asarray(a["flatfield"]), np.asarray(b["flatfield"]))
    assert np.asarray(a["minval"]).tobytes() == np.asarray(b["minval"]).tobytes()


# ---------------------------------------------------------------------------------------------
# more of the port's own
# ---------------------------------------------------------------------------------------------


def test_an_existing_result_is_kept(tmp_path, capsys):
    f = tmp_path / "im.tif"
    _write(f, dark_stack("t/keep"))
    res = tmp_path / "bg.mat"
    res.write_bytes(b"old")
    fe.flatfield_background_estimation([str(f)], str(res))
    assert res.read_bytes() == b"old"
    assert "already exists, skip it!" in capsys.readouterr().out
    fe.flatfield_background_estimation([str(f)], str(res), overwrite=True)
    assert float(read_mat(res, ["numSlice"])["numSlice"]) > 0


def test_the_estimate_does_not_depend_on_the_thread_count():
    st = bright_stack("t/threads", nz=60)
    a = fe.estimate_background(st, num_threads=1)
    b = fe.estimate_background(st, num_threads=4)
    for k in ("mu", "sigma", "med"):
        assert a[k].tobytes() == b[k].tobytes()
    assert a["numSlice"] == b["numSlice"]


def test_estimate_keeps_the_quieter_half_without_bright_slices():
    st = bright_stack("t/drop", nz=40)
    r = fe.estimate_background(st)
    # the 6 bright slices are noisier than the median slice: gone with the noisier half
    assert r["numSlice"] == 20
    assert abs(float(r["mu"].astype(np.float64).mean()) - DARK) < 0.3
    assert r["imp"] is None
    assert r["mu"].dtype == np.float32 and r["mu"].shape == (48, 40)
    np.testing.assert_array_equal(r["pmat"], np.arange(0, 101, 5, dtype=np.float64))


@pytest.mark.parametrize("nz", [7, 8])
def test_median_by_selection_is_matlabs_median(nz):
    from petakit5d.matlabcompat import stats

    v = np.asarray(_rng(f"med{nz}").integers(-3, 3, (nz, 9, 11)), dtype=np.float32)
    v[v == -3] = -0.0  # signed zeros: MATLAB's stable sort picks the sign
    v[1, 2, 3] = dtypes.NAN
    v[:, 4, 5] = np.inf
    v[: nz // 2, 6, 7] = -np.inf
    for nt in (1, 3):
        got = fe._median_z(v, nt)
        want = np.asarray(stats.median(v, axis=0), dtype=np.float32)
        assert got.view(np.uint32).tobytes() == want.view(np.uint32).tobytes()


def test_an_all_excluded_stack_has_no_slices():
    r = fe.estimate_background(np.full((6, 8, 10), 7, np.uint16))
    assert r["numSlice"] == 0
    assert np.isnan(r["mu"]).all() and np.isnan(r["med"]).all()


def test_cell_order_is_matlabs_linear_order():
    # X = cell(1, 2, 2); X{1, 2} = 'a'; X{2, 2} = 'b'; X{2, 1} = 'c'  -> c, a, b
    assert fe._cell_order({(1, 2): "a", (2, 2): "b", (2, 1): "c"}) == ["c", "a", "b"]
    # X = cell(1, 2, 1, 2); X{1, 2, 2} = 1; X{2, 1, 1} = 2  -> 2, 1
    assert fe._cell_order({(1, 2, 2): 1, (2, 1, 1): 2}) == [2, 1]


def test_phase_without_basicpy_raises(monkeypatch, tmp_path):
    import builtins

    real = builtins.__import__

    def no_basicpy(name, *a, **k):
        if name == "basicpy" or name.startswith("basicpy."):
            raise ImportError("no basicpy")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_basicpy)
    with pytest.raises(ImportError, match="BaSiCPy"):
        fb.basic_flatfield(np.ones((3, 16, 16), np.float32))


def test_phase_skips_when_all_outputs_exist(tmp_path, capsys):
    f = tmp_path / "t.tif"
    _write(f, phase_tile((32, 32), 0, "t/skip"))
    out = str(tmp_path / "o") + os.sep
    os.makedirs(out)
    Path(out + "t.tif").write_bytes(b"x")
    fb.phase_image_flatfield_correction([str(f)], out, 1100)
    assert "All output files exist, skip it!" in capsys.readouterr().out


def test_phase_divides_and_rounds_as_matlab(tmp_path, monkeypatch):
    # a given flat field (no BaSiCPy): single(pre) ./ flatfield, then uint16 (rounded, saturated)
    ff = np.full((32, 40), 1.3) + np.arange(40) / 400.0
    monkeypatch.setattr(fb, "basic_flatfield", lambda st: ff)
    fns = []
    for k in range(3):
        fns.append(str(tmp_path / f"t{k}.tif"))
        _write(fns[-1], phase_tile((32, 40), k, "t/round"))
    out = str(tmp_path / "o") + os.sep
    os.makedirs(out)
    fb.phase_image_flatfield_correction(fns, out, 1100)
    pre = [fb.phase_prefilter(read_tiff(f), 1100) for f in fns]
    for k in range(3):
        want = np.asarray(dtypes.cast(dtypes.rdivide(pre[k], ff), "uint16"))
        np.testing.assert_array_equal(read_tiff(f"{out}t{k}.tif"), want)
    a = read_mat(f"{out}BaSiC_flatfield/t0_flatfield.mat")
    assert a["flatfield"].dtype == np.float64
    assert a["minval"] == np.min(np.stack(pre))


def test_phase_char_list_is_matlabs_error():
    with pytest.raises(fb.MatlabError) as e:
        fb.phase_image_flatfield_correction("a.tif", "out/", 1100)
    assert e.value.identifier == "MATLAB:cellRefFromNonCell"
