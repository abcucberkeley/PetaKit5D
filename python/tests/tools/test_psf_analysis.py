"""PSF analysis, detection and cropping (petakit5d.tools.psf_analysis, task B06).

MATLAB's tests (tests/analysis/PSFAnalysisTest.m) and where they went:

- plotDataMatchWholeGrid (4 cases: centered, paddedZ, largerThanGrid, maxOutsideGrid): the
  PSF and OTF slices of Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update, which MATLAB computes from
  the data regions of its 1001^3 grids, against a whole-grid reference written in the test.
  Here the same four PSFs are oracle cases (``otf.<case>``) of the function itself: the port
  computes from data regions too, and MATLAB's results are the reference. More cases: a
  Richards & Wolf-like single PSF without background, an off-center PSF with odd sizes and
  other wavelengths, NA, gamma and background factor, a PSF longer in z than the OTF's region,
  and the fattest column outside the OTF (MATLAB's index error).
- analysisWritesFiguresAndInfo: the oracle case ``wrapper.deskew`` (XR_psf_analysis_wrapper on
  the test's PSF with deskew and a Richards & Wolf PSF: the deskewed PSF, RW_info.mat and the
  info file compared, the four figure files present on both sides) and
  ``test_analysis_writes_figures_and_info`` (the files, the deskewed size, and the PNG renamed
  into place last, after the PDF).
- detectionInfoIsMovedIntoPlace: the oracle case ``detect.two_beads`` (the test's image) and
  ``test_detection_info_is_moved_into_place`` (the info file renamed into place).

Each analysis case runs on both FFT backends: ``.exact`` (bit for bit; skipped without the
add-on or where FFTW's twiddle factors differ from the reference platform's) and
``.permissive`` (the PSF slices bit for bit, the OTF's within ``OTF_TOL``,
``docs/parity/tools.md``). Detection uses no FFT; its Radon transform is not MATLAB's in the
last bits (``docs/parity/tools.md``), so its cases compare what it decides: the peaks and the
crops, exactly.

Inputs are made with exact arithmetic (``petakit5d.testing.synth``); MATLAB reads the files.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.core.axes import from_yxz
from petakit5d.core.dtypes import uint16
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, File, InPath, Tree, Work
from petakit5d.testing.compare import ignore, tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import figures
from petakit5d.tools import psf_analysis as PA

HAS_MPL = figures.available()
needs_mpl = pytest.mark.skipif(not HAS_MPL, reason="matplotlib (the figures extra) is missing")
BACKENDS = ("exact", "permissive")
OTF_TOL = tol(
    reason="0006: double results of the permissive FFT backend (max <= 1e-13 of the largest "
    "value, relative L2 <= 1e-14); the OTF slices and line cuts, docs/parity/tools.md",
    rel_max=1e-13,
    rel_l2=1e-14,
)

# ---------------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------------


def bead_volume(msize, centers, sigma, amp=2000.0, bg=100.0, noise=0.0, seed="b06", second=()):
    """pktest.beads in the port's form: MATLAB's size [y x z], 1-based centers [y x z] and sigma
    [yx z]; plus ``noise`` times uniform noise; single (z, y, x)"""
    sz = from_yxz(tuple(int(v) for v in msize))
    pos = np.atleast_2d(np.asarray(centers, dtype=np.float64))
    pos = np.stack([pos[:, 2] - 1, pos[:, 0] - 1, pos[:, 1] - 1], axis=1)
    im = synth.beads(sz, pos, (sigma[1], sigma[0]), amp, bg)
    for c2, s2, a2 in second:
        p2 = np.array([[c2[2] - 1, c2[0] - 1, c2[1] - 1]], dtype=np.float64)
        im = im + synth.beads(sz, p2, (s2[1], s2[0]), a2, 0.0)
    if noise:
        im = im + np.float32(noise) * CounterRNG(0, seed).uniform(sz, np.float32)
    return im


def psf_file(msize, centers, sigma, cls="uint16", **kw):
    def make(path):
        im = bead_volume(msize, centers, sigma, **kw)
        synth.write_tiff(path, uint16(im) if cls == "uint16" else im.astype(np.float32))

    return make


# ---------------------------------------------------------------------------------------------
# oracle cases: Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update
# ---------------------------------------------------------------------------------------------


def py_analysis(fn, src, xy, z, na, idx, exc, det, sub, gamma, bg, visible, *, backend):
    """the MATLAB arguments (PSFsubpix [y x z]) to psf_otf_analysis, on one FFT backend"""
    with P.using(backend):
        r = PA.psf_otf_analysis(
            fn, src, xy, z, na, idx, exc, det, from_yxz(tuple(sub)), gamma, bg, visible, plot=False
        )
    return list(r.outputs())


def otf_cases(cid, make, zpix, sub, *, cls="uint16", exc=488, det=515, na=1.0, gamma=0.5, bg=1.5,
              marks=(), error=False):  # fmt: skip
    out = []
    for be in BACKENDS:
        compare = {} if be == "exact" else dict.fromkeys(range(3, 11), OTF_TOL)
        out.append(
            Case(
                f"otf.{cid}.{be}",
                matlab="Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update",
                nargout=11,
                args=[InPath("psf"), "test", 108, zpix, na, 1.33, exc, det, sub, gamma, bg, False],
                python=py_analysis,
                py_kwargs={"backend": be},
                inputs={"psf": File(make, "psf.tif")},
                files=False,
                expect_error="MATLAB:badsubscript" if error else False,
                py_error=(PA.PsfAnalysisError,),
                marks=marks,
                **({"compare": compare} if compare else {}),
            )
        )
    if not error:
        # the PSF slices alone (no FFT): bit for bit on every backend, in CI too
        out.append(
            Case(
                f"otf.{cid}.psf",
                matlab="Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update",
                nargout=3,
                args=[InPath("psf"), "test", 108, zpix, na, 1.33, exc, det, sub, gamma, bg, False],
                python=lambda *a: py_analysis(*a, backend="permissive")[:3],
                inputs={"psf": File(make, "psf.tif")},
                files=False,
                marks=marks,
            )
        )
    return out


def _otf_cases():
    slow = ("slow",)
    noisy = {"noise": 5.0}
    return [
        # PSFAnalysisTest's four PSFs
        *otf_cases("centered", psf_file([64, 64, 101], [32, 32, 51], [1.6, 3], seed="c1", **noisy),
                   160, [128, 128, 63]),
        *otf_cases("paddedZ", psf_file([64, 60, 21], [30, 33, 11], [1.6, 2], seed="c2", **noisy),
                   300, [128, 128, 35]),
        *otf_cases("largerThanGrid",
                   psf_file([380, 370, 31], [185, 190, 15], [2, 3], seed="c3", **noisy),
                   160, [128, 128, 63], marks=slow),
        *otf_cases("maxOutsideGrid",
                   psf_file([360, 360, 31], [180, 180, 16], [0.5, 0.8], seed="c4",
                            second=(([180, 355, 16], [4, 3], 1990.0),), **noisy),
                   160, [128, 128, 63], marks=slow),
        # the Richards & Wolf PSFs' form: single, no background, bgFactor 0
        *otf_cases("rw", psf_file([128, 128, 101], [64, 64, 51], [2, 6], cls="single", amp=1.0,
                                  bg=0.0), 100, [128, 128, 101], bg=0),
        *otf_cases("offcenter", psf_file([57, 63, 41], [20.3, 40.6, 15.2], [1.5, 2.5], seed="c6",
                                         **noisy),
                   200, [128, 128, 51], exc=560, det=605, na=1.1, gamma=0.7, bg=1.2),
        *otf_cases("tallZ", psf_file([64, 64, 151], [33, 31, 80], [1.6, 4], seed="c7", **noisy),
                   100, [128, 128, 101]),
        # the fattest column outside the OTF's central part: MATLAB's index error
        *otf_cases("fattestOutside", psf_file([64, 64, 101], [32, 32, 51], [1.6, 3], seed="c1",
                                              **noisy), 160, [128, 128, 63], na=3.0, error=True),
    ]  # fmt: skip


# ---------------------------------------------------------------------------------------------
# oracle cases: XR_psf_analysis_wrapper
# ---------------------------------------------------------------------------------------------

FIGS = {"*.png": ignore("figures: matplotlib's pixels"), "*.pdf": ignore("figures: matplotlib's")}
WRAPPER_TREE = tree(per_file=FIGS, ignore_vars=("t", "t1", "ip"))


def psf_dataset(names, deskewed=False):
    """PSF/ with the PSFs ``names`` (the test's bead, one per name, shifted a little) and
    PSF/RW_PSFs/rw<k>.tif (the test's Richards & Wolf PSF, two widths)"""

    def make(folder):
        root = os.path.join(folder, "PSF")
        os.makedirs(os.path.join(root, "RW_PSFs"))
        for k, name in enumerate(names):
            c = [32 + 0.3 * k, 32 - 0.2 * k, 51]
            im = bead_volume([64, 64, 101], c, [1.6, 3] if not deskewed else [1.5, 2.5])
            synth.write_tiff(os.path.join(root, name), uint16(im))
        for k in range(2):
            rw = bead_volume([128, 128, 101], [64, 64, 51], [2 + k * 0.4, 6], amp=1.0, bg=0.0)
            synth.write_tiff(os.path.join(root, "RW_PSFs", f"rw{k}.tif"), rw)

    return make


def py_wrapper(d, *args):
    from petakit5d import matlab

    a = list(args)
    rw = a[a.index("RWFn") + 1]
    a[a.index("RWFn") + 1] = [d + r for r in rw]
    matlab.XR_psf_analysis_wrapper([d + "PSF/"], *a)
    return []


def py_wrapper_mpl(d, *args, backend):
    pytest.importorskip("matplotlib")
    with P.using(backend):
        return py_wrapper(d, *args)


_WRAPPER = (
    "a = in(2 : end); i = find(strcmp(a, 'RWFn')); a{i + 1} = strcat(in{1}, a{i + 1}); "
    "XR_psf_analysis_wrapper({[in{1}, 'PSF/']}, a{:}); out = {};"
)
_RUN = ["sourceStr", "test", "visible", False, "parseCluster", False, "masterCompute", True]


def wrapper_case(cid, make, *args, error=False):
    return [
        Case(
            f"wrapper.{cid}.{be}",
            code=_WRAPPER,
            nargout=0,
            args=[InPath("d"), *args, *_RUN],
            python=py_wrapper_mpl,
            py_kwargs={"backend": be},
            inputs={"d": Tree(make, "d")},
            tree=WRAPPER_TREE if be == "exact" else WRAPPER_TREE_PERMISSIVE,
            expect_error=error,
            py_error=(PA.PsfAnalysisError,),
            marks=("slow",),
        )
        for be in BACKENDS
    ]


def py_plot(*args):
    from petakit5d import matlab

    matlab.XR_psf_analysis_plot(*args)
    return []


def plot_exists(folder):
    """a PSF and an existing figure: the task does nothing"""
    psf_file([64, 64, 101], [32, 32, 51], [1.6, 3])(os.path.join(folder, "p.tif"))
    os.makedirs(os.path.join(folder, "PSFAnalysis"))
    with open(os.path.join(folder, "PSFAnalysis", "wT_p.png"), "wb") as f:
        f.write(b"not a png")


# the info files hold the OTF's slices and line cuts
WRAPPER_TREE_PERMISSIVE = tree(values=OTF_TOL, per_file=FIGS, ignore_vars=("t", "t1", "ip"))


def _wrapper_cases():
    return [
        # PSFAnalysisTest/analysisWritesFiguresAndInfo
        *wrapper_case(
            "deskew",
            psf_dataset(["488_CamA_ch0.tif"]),
            "xyPixelSize", 0.108, "dz", 0.3, "skewAngle", 32.45, "channelPatterns", ["CamA_ch0"],
            "channels", 488, "deskew", True, "RWFn", ["PSF/RW_PSFs/rw0.tif"],
        ),
        # two channels, PSFs already deskewed (objective scan: z pixel size dz), a file of no
        # channel
        *wrapper_case(
            "two_channels",
            psf_dataset(["a_CamA_ch0.tif", "b_CamB_ch0.tif", "c_CamC_ch0.tif"], deskewed=True),
            "dz", 0.2, "channelPatterns", ["CamA_ch0", "CamB_ch0"], "channels", [560, 488],
            "deskew", False, "objectiveScan", True, "bgFactor", 1.2,
            "RWFn", ["PSF/RW_PSFs/rw1.tif", "PSF/RW_PSFs/rw0.tif"],
        ),
        # a first channel MATLAB's switch doesn't know: exc_lambda is undefined
        *wrapper_case(
            "unknown_channel",
            psf_dataset(["a_CamA_ch0.tif"], deskewed=True),
            "channelPatterns", ["CamA_ch0"], "channels", 405, "deskew", False,
            "RWFn", ["PSF/RW_PSFs/rw0.tif"], error="MATLAB:UndefinedFunction",
        ),
    ]  # fmt: skip


PLOT_CASES = [
    # the figure exists: nothing is done
    Case(
        "plot.exists_skips",
        code=(
            "XR_psf_analysis_plot([in{1}, 'p.tif'], [in{1}, 'PSFAnalysis/wT_p.png'], "
            "[in{1}, 'none.mat'], in{2 : end}); out = {};"
        ),
        nargout=0,
        args=[
            InPath("d"),
            1,
            "t",
            108,
            160,
            1,
            1.33,
            488,
            515,
            [128, 128, 63],
            0.5,
            1.5,
            "visible",
            False,
        ],  # fmt: skip
        python=lambda d, *a: py_plot(f"{d}p.tif", f"{d}PSFAnalysis/wT_p.png", f"{d}none.mat", *a),
        inputs={"d": Tree(plot_exists, "d")},
        tree=WRAPPER_TREE,
    ),
]


# ---------------------------------------------------------------------------------------------
# oracle cases: XR_psf_detection_and_cropping
# ---------------------------------------------------------------------------------------------

DETECT_TREE = tree(ignore_vars=("t", "t1"))


def py_detect(fn, result_dir, *args):
    from petakit5d import matlab

    matlab.XR_psf_detection_and_cropping(fn, result_dir, *args)
    return []


def beads_image(msize, centers, sigma, amps, noise, seed):
    def make(path):
        im = bead_volume(msize, centers, sigma, amp=np.asarray(amps, np.float64), noise=noise,
                         seed=seed)  # fmt: skip
        synth.write_tiff(path, uint16(im))

    return make


def beads_parts(msize, centers, sigma, amps, noise, seed, split):
    def make(folder):
        im = uint16(bead_volume(msize, centers, sigma, amp=np.asarray(amps, np.float64),
                                noise=noise, seed=seed))  # fmt: skip
        synth.write_tiff(os.path.join(folder, "im_part000.tif"), im[:split])
        synth.write_tiff(os.path.join(folder, "im_part001.tif"), im[split:])

    return make


MANY = {
    "msize": [160, 140, 120],
    # isolated ones, a close pair, one at the border, a saturated one, a dim one
    "centers": [[40, 40, 30], [110, 95, 80], [70, 30, 95], [72, 36, 100], [5, 70, 60],
                [130, 110, 30], [100, 40, 40]],
    "sigma": [1.6, 3],
    "amps": [2000, 1500, 1800, 1700, 2100, 70000, 300],
    "noise": 30.0,
    "seed": "many",
}  # fmt: skip


def _detect_cases():
    two = beads_image([128, 128, 101], [[40, 40, 30], [90, 80, 70]], [1.6, 3], [2000, 2000], 0.0,
                      "two")  # fmt: skip
    many = beads_image(**MANY)
    nv = ["cropSize", [40, 36, 31], "distThresh", [30, 30, 40], "dz", 0.2, "prefix", "many"]
    return [
        # PSFAnalysisTest/detectionInfoIsMovedIntoPlace
        Case(
            "detect.two_beads",
            matlab="XR_psf_detection_and_cropping",
            nargout=0,
            args=[
                InPath("im"),
                Work("Cropped/"),
                "cropSize",
                [32, 32, 21],
                "distThresh",
                [32, 32, 21],
                "prefix",
                "beads",
            ],  # fmt: skip
            python=py_detect,
            inputs={"im": File(two, "beads.tif")},
            tree=DETECT_TREE,
        ),
        Case(
            "detect.many",
            matlab="XR_psf_detection_and_cropping",
            nargout=0,
            args=[InPath("im"), Work("Cropped/"), *nv],
            python=py_detect,
            inputs={"im": File(many, "many.tif")},
            tree=DETECT_TREE,
        ),
        # the defaults (crops of 256 x 128 x 201: none fits) and another angle
        Case(
            "detect.defaults",
            matlab="XR_psf_detection_and_cropping",
            nargout=0,
            args=[InPath("im"), Work("Cropped/"), "skewAngle", 45, "xyPixelSize", 0.1],
            python=py_detect,
            inputs={"im": File(many, "many.tif")},
            tree=DETECT_TREE,
        ),
        # an image in two parts along z
        Case(
            "detect.parts",
            code=(
                "XR_psf_detection_and_cropping({[in{1}, 'im_part000.tif'], "
                "[in{1}, 'im_part001.tif']}, in{2 : end}); out = {};"
            ),
            nargout=0,
            args=[InPath("d"), Work("Cropped/"), *nv],
            python=lambda d, *a: py_detect([d + "im_part000.tif", d + "im_part001.tif"], *a),
            inputs={"d": Tree(beads_parts(**MANY, split=47), "d")},
            tree=DETECT_TREE,
        ),
        # the info file exists: nothing is done
        Case(
            "detect.exists_skips",
            code="XR_psf_detection_and_cropping([in{1}, 'beads.tif'], [in{1}, 'C/'], in{2 : end}); "
            "out = {};",
            nargout=0,
            args=[InPath("d"), "prefix", "beads"],
            python=lambda d, *a: py_detect(f"{d}beads.tif", f"{d}C/", *a),
            inputs={"d": Tree(detect_exists, "d")},
            tree=DETECT_TREE,
        ),
    ]


def detect_exists(folder):
    beads_image([128, 128, 101], [[40, 40, 30]], [1.6, 3], [2000], 0.0, "ex")(
        os.path.join(folder, "beads.tif")
    )
    os.makedirs(os.path.join(folder, "C"))
    with open(os.path.join(folder, "C", "beads.mat"), "wb") as f:
        f.write(b"not a mat file")


CASES = [*_otf_cases(), *_wrapper_cases(), *PLOT_CASES, *_detect_cases()]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# PSFAnalysisTest, ported
# ---------------------------------------------------------------------------------------------


@needs_mpl
@pytest.mark.slow
def test_analysis_writes_figures_and_info(tmp_path, monkeypatch):
    """the deskewed PSF, the figures and info files; the PNG (the task's result) is renamed into
    place last, after the PDF was written"""
    from petakit5d.io import read_tiff

    psf_dataset(["488_CamA_ch0.tif"])(str(tmp_path))
    root = f"{tmp_path}/PSF/"
    moves = {}
    real = os.replace

    def spy(src, dst):
        moves[os.fspath(dst)] = set(os.listdir(os.path.dirname(os.fspath(dst))))
        real(src, dst)

    monkeypatch.setattr(PA.os, "replace", spy)
    PA.psf_analysis_dataset(
        root,
        xy_pixel_size=0.108,
        dz=0.3,
        skew_angle=32.45,
        channel_patterns=("CamA_ch0",),
        channels=(488,),
        deskew=True,
        rw_fn=(f"{root}RW_PSFs/rw0.tif",),
        source_str="test",
        visible=False,
        parse_cluster=False,
        master_compute=True,
    )
    ds = read_tiff(f"{root}DS/488_CamA_ch0.tif")
    assert ds.shape[0] == 101 and ds.shape[1] == 64 and ds.shape[2] > 64
    res = f"{root}DS/PSFAnalysis/"
    fsn = "488_CamA_ch0"
    for f in (f"wT_{fsn}.png", f"wT_{fsn}.pdf", f"comp_{fsn}.png", f"comp_{fsn}.pdf",
              f"{fsn}_infos.mat", "RW_info.mat"):  # fmt: skip
        assert os.path.isfile(res + f), f
    png = os.path.normpath(res + f"wT_{fsn}.png")
    hit = [k for k in moves if os.path.normpath(k) == png]
    assert hit, "the PNG was not renamed into place"
    assert f"wT_{fsn}.pdf" in moves[hit[0]], "the PDF was not written before the PNG"


def test_detection_info_is_moved_into_place(tmp_path, monkeypatch):
    """the detection's info file (the task's result) is saved under another name and renamed"""
    from petakit5d.io import read_mat

    fn = str(tmp_path / "beads.tif")
    beads_image([128, 128, 101], [[40, 40, 30], [90, 80, 70]], [1.6, 3], [2000, 2000], 0.0,
                "two")(fn)  # fmt: skip
    result = f"{tmp_path}/Cropped/"
    moved = []
    real = os.replace

    def spy(src, dst):
        moved.append((os.fspath(src), os.fspath(dst)))
        real(src, dst)

    monkeypatch.setattr(PA.os, "replace", spy)
    PA.psf_detection_and_cropping(fn, result, crop_size=(21, 32, 32), dist_thresh=(21, 32, 32),
                                  prefix="beads")  # fmt: skip
    assert [d for _, d in moved] == [f"{result}beads.mat"]
    assert moved[0][0] != moved[0][1]
    a = read_mat(f"{result}beads.mat")
    assert "pcoords" in a
    pc = np.asarray(a["pcoords"])
    np.testing.assert_array_equal(pc[np.argsort(pc[:, 0])], [[40, 40, 30], [90, 80, 70]])


# ---------------------------------------------------------------------------------------------
# the port's own rules
# ---------------------------------------------------------------------------------------------


def test_colon_matches_the_symmetric_rule():
    """MATLAB's a:d:b (probed): the count with its tolerance, the symmetric halves, the middle"""
    v = PA._colon(1, 0.1, 1.3)
    assert v.size == 4 and v[-1] == 1.3
    v = PA._colon(1, 101 / 303, 101)
    assert v.size == 301 and v[-1] == 101
    v = PA._colon(1, 0.7, 19.9)
    assert v[0] == 1 and v[-1] == 1 + 27 * 0.7
    assert PA._colon(1, 0.5, 0.5).size == 0


def test_interpn_is_linear_along_y_then_x_then_z():
    rng = CounterRNG(1, "interp")
    v = rng.uniform((5, 6, 7))
    qy = PA._colon(1, 6 / 11, 6)
    qx = PA._colon(1, 7 / 9, 7)
    qz = PA._colon(1, 5 / 8, 5)
    out = PA._interpn(v, qy, qx, qz)
    assert out.shape == (qz.size, qy.size, qx.size)
    # the steps one after the other, as written
    ref = v
    for ax, q in ((1, qy), (2, qx), (0, qz)):
        n = ref.shape[ax]
        i = np.clip(np.floor(q).astype(int), 1, n - 1)
        t = (q - i).reshape([-1 if d == ax else 1 for d in range(3)])
        ref = (1 - t) * np.take(ref, i - 1, axis=ax) + t * np.take(ref, i, axis=ax)
    np.testing.assert_array_equal(out, ref)
    np.testing.assert_array_equal(PA._interpn(v, qy, qx, qz, 3), out)
    # at the grid points the values themselves
    np.testing.assert_array_equal(PA._interpn(v, *(np.arange(1.0, n + 1) for n in (6, 7, 5))), v)


def test_grid_circshift_equals_the_whole_grid(monkeypatch):
    """the box moved, or the grid turned (lazily), equals np.roll of the whole grid, in its
    slices along each axis and its first maximum (small grid)"""
    monkeypatch.setattr(PA, "_GRID", 11)
    rng = CounterRNG(2, "shift")
    for v_shape, off, shifts in (
        ((3, 4, 11), (2, 3, 0), [(1, -2, 3)]),
        ((3, 4, 5), (8, 0, 6), [(2, -1, 1), (-3, 4, 0)]),
        ((11, 11, 11), (0, 0, 0), [(-4, 5, 1), (2, 2, -7)]),
    ):
        v = rng.uniform(v_shape)
        full = np.zeros((11, 11, 11))
        full[tuple(slice(o, o + n) for o, n in zip(off, v_shape, strict=True))] = v
        g = PA._Grid(v.copy(), off, 2)
        for sh in shifts:
            g.circshift(sh)
            full = np.roll(full, sh, axis=(0, 1, 2))
            for axis in range(3):
                for k in range(11):
                    np.testing.assert_array_equal(g.slice(axis, k), np.take(full, k, axis=axis))
            m, place = g.max()
            assert m == full.max()
            assert place == PA._first_in_matlab_order(np.argwhere(full == full.max()))


def test_first_maximum_is_first_in_matlabs_order():
    """ties: the smallest z, then x, then y (MATLAB's linear order, y fastest)"""
    places = np.array([[1, 0, 5], [1, 3, 2], [1, 2, 2], [2, 0, 0]])
    assert PA._first_in_matlab_order(places) == (1, 2, 2)
    v = np.zeros((3, 4, 5))
    v[1, 3, 2] = v[1, 2, 2] = v[1, 0, 5 - 1] = 1.0
    g = PA._Grid(v, (0, 0, 0), 1)
    assert g.max() == (1.0, (1, 2, 2))


def _detected(tmp_path, centers, amps, crop, dist):
    from petakit5d.io import read_mat

    fn = str(tmp_path / "im.tif")
    beads_image([128, 128, 101], centers, [1.6, 3], amps, 0.0, "iso")(fn)
    PA.psf_detection_and_cropping(
        fn, f"{tmp_path}/C/", crop_size=from_yxz(crop), dist_thresh=from_yxz(dist), prefix="b"
    )
    pc = np.asarray(read_mat(f"{tmp_path}/C/b.mat")["pcoords"]).reshape(-1, 3)
    return pc[np.argsort(pc[:, 0])] if pc.size else pc


def test_detection_keeps_isolated_psfs_only(tmp_path):
    """two PSFs within the distance in y, x and z exclude each other; an isolated one stays"""
    pc = _detected(tmp_path, [[40, 40, 30], [40, 70, 30], [95, 95, 75]], [2000] * 3,
                   [21, 21, 15], [32, 50, 21])  # fmt: skip
    np.testing.assert_array_equal(pc, [[95, 95, 75]])


def test_saturated_psfs_are_left_out(tmp_path):
    """a crop whose maximum is above 60000 is not written and its peak not kept"""
    pc = _detected(tmp_path, [[40, 40, 30], [95, 95, 75]], [40000, 80000], [21, 21, 15],
                   [32, 32, 21])  # fmt: skip
    np.testing.assert_array_equal(pc, [[40, 40, 30]])
    assert os.path.isfile(f"{tmp_path}/C/b_crop_40_40_30.tif")
    assert not os.path.isfile(f"{tmp_path}/C/b_crop_95_95_75.tif")


def test_local_threshold_is_the_plain_rule():
    """the peaks the local test keeps are those of MATLAB's rule written plainly (matlabcompat's
    median and std of every window): the NumPy pre-decision changes nothing"""
    from petakit5d.matlabcompat.stats import median, std

    im = uint16(bead_volume(**{k: v for k, v in MANY.items() if k != "amps"},
                            amp=np.asarray(MANY["amps"], np.float64)))  # fmt: skip
    nz, ny, nx = im.shape
    theta = PA._atand(PA._tand(32.45) * 0.108 / (0.2 * PA._sind(32.45)))
    r = PA._radon(np.transpose(im, (1, 2, 0)), theta)
    px, py = np.nonzero(PA._regional_max(r).T)
    s_all = float(np.asarray(std(PA._column_major(r))).reshape(-1)[0])
    h, w = r.shape
    kept = []
    for yi, xi in zip(py.tolist(), px.tolist(), strict=True):
        v = PA._column_major(r[max(0, yi - 5) : min(h, yi + 6), max(0, xi - 5) : min(w, xi + 6)])
        bg = float(np.asarray(median(v)).reshape(-1)[0])
        sd = float(np.asarray(std(v)).reshape(-1)[0])
        if r[yi, xi] - bg > 2 * (sd * 0.95 + s_all * 0.05):
            kept.append(yi + 1)
    got = PA._detect(im, 0.108, 0.2, 32.45)
    assert 0 < len(kept) < py.size
    np.testing.assert_array_equal(got[:, 0], kept)


def test_regional_max_takes_plateaus_and_borders():
    a = np.zeros((6, 7))
    a[1:3, 1:3] = 5
    a[5, 6] = 2
    a[3, 5] = 4
    a[3, 6] = 4
    bw = PA._regional_max(a)
    assert bw[1:3, 1:3].all() and bw[5, 6] and bw[3, 5] and bw[3, 6]
    assert bw.sum() == 7
    assert PA._regional_max(-np.ones((4, 5))).all()


def test_radon_projects_a_pixel_to_its_bins():
    """a single pixel's four quarters land in the bins around x cos + y sin, summing to it"""
    im = np.zeros((1, 17, 9))
    im[0, 0, 0] = 1.0
    r = PA._radon(im, 51.99173)[0]
    assert r.size == 2 * 10 + 1
    assert abs(r.sum() - 1) < 1e-15
    np.testing.assert_array_equal(np.flatnonzero(r), [13, 14, 15])


def test_results_do_not_depend_on_the_thread_count(tmp_path):
    fn = str(tmp_path / "psf.tif")
    psf_file([64, 64, 101], [32, 32, 51], [1.6, 3], noise=5.0, seed="c1")(fn)
    args = (fn, "t", 108, 160, 1.0, 1.33, 488, 515, (63, 128, 128))
    ref = PA.psf_otf_analysis(*args, plot=False, num_threads=1).outputs()
    got = PA.psf_otf_analysis(*args, plot=False, num_threads=3).outputs()
    for a, b in zip(ref, got, strict=True):
        np.testing.assert_array_equal(a, b)


@needs_mpl
def test_analysis_figure_has_matlabs_panels(tmp_path):
    fn = str(tmp_path / "psf.tif")
    psf_file([64, 64, 101], [32, 32, 51], [1.6, 3], noise=5.0, seed="c1")(fn)
    r = PA.psf_otf_analysis(fn, "demo", 108, 160, 1.0, 1.33, 488, 515, (63, 128, 128), 0.5, 1.5,
                            False)  # fmt: skip
    gids = {ax.get_gid() for ax in r.figure.axes}
    assert {"pk5d:psf_xy", "pk5d:psf_xz", "pk5d:psf_yz", "pk5d:otf_xy", "pk5d:otf_xz",
            "pk5d:otf_yz", "pk5d:linecuts_xz", "pk5d:linecuts_yz"} <= gids  # fmt: skip
    assert len(r.figure.axes) == 10
    w, h = r.figure.get_size_inches() * r.figure.dpi
    assert (round(w), round(h)) == (2500, 1100)


def test_unknown_first_channel_is_matlabs_error():
    with pytest.raises(PA.PsfAnalysisError) as e:
        PA._wavelengths(405, None)
    assert e.value.identifier == "MATLAB:UndefinedFunction"
    assert PA._wavelengths(405, (560.0, 605.0)) == (560.0, 605.0)
    assert PA._wavelengths(642, None) == (642.0, 680.0)
