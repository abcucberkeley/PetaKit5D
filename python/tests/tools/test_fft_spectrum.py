"""FFT power spectra (petakit5d.tools.fft_spectrum, task B04).

Oracle cases against MATLAB: XR_fftSpectrumComputingFrame (FFTSpectrumTest's cases: the
resized image cropped, padded in y, padded in every axis, the 3D stack, background
subtraction, cubic interpolation, an even size equal to the resized one; and an image larger
than the region kept, so the read box has non-integer ends, zarr input, nearest interpolation,
shrinking, the output pixel size left empty, a non-integer background, single input, one
number as the size, the errors of even sizes and of one z plane, an existing result) and
XR_fftSpectrumComputingWrapper (FFTSpectrumTest's cosine and anisotropic data, the 3D stack,
channel patterns, a zarr dataset, values through the task text, no matching file). Each case
runs on both CPU FFT backends (``.exact``: MATLAB's bits; ``.permissive``: decision record
0006, the normalized spectrum within its tolerance). Probes of MATLAB's single ``log10`` and
``.^ 0.5`` pin the arithmetic of the outputs.

Then FFTSpectrumTest in intent (the cosine's peaks, the log spectrum, the anisotropic input,
the region-only resize against the whole spectrum, the even-size error) and unit tests of the
read box, the orientation of the images and the task text.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.matlabcompat import elementary as E
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value
from petakit5d.testing.compare import Custom, MArray, tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import fft_spectrum as fs

BACKENDS = ("exact", "permissive")
ERRORS = (fs.FFTSpectrumError, ValueError, OSError, RuntimeError, IndexError)
# FFTSpectrumTest's sizes (MATLAB's [y x z]) as (z, y, x)
CROP, PAD_Y, EVEN = (12, 31, 30), (12, 28, 30), (12, 30, 30)
N_COS, K_COS = 65, 8


# ---------------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------------


def textured(shape, name, cls="uint16"):
    """FFTSpectrumTest's image: 1000 + 300 texture + 30 noise (uint16), deterministic"""
    rng = CounterRNG(4, f"b04/{name}")
    t = synth.texture(shape, 1.5, (0.0, 1.0), rng=rng.spawn("t")).astype(np.float64)
    v = 1000.0 + 300.0 * t + 30.0 * rng.spawn("n").normal(shape)
    if cls == "single":
        return np.ascontiguousarray(v.astype(np.float32))
    return np.floor(v + 0.5).astype(np.uint16)


def cosine():
    """FFTSpectrumTest.cosineImage: 1000 + 500 cos(2 pi K x / N) along x (single)"""
    x = np.arange(N_COS, dtype=np.float64)
    c = E.cos(2 * np.pi * K_COS * x / N_COS)
    row = (1000.0 + 500.0 * c).astype(np.float32)
    return np.ascontiguousarray(np.broadcast_to(row, (N_COS, N_COS, N_COS)))


def frame_tree(vol, fmt="tif", existing=False):
    def make(folder):
        folder = Path(folder)
        if fmt == "zarr":
            synth.write_zarr(folder / "im.zarr", vol(), block_size=(8, 16, 16), layout="matlab")
        else:
            synth.write_tiff(folder / "im.tif", vol())
        if existing:
            (folder / "FFT").mkdir()
            synth.write_tiff(folder / "FFT" / "im.tif", np.arange(4, dtype=np.uint16).reshape(2, 2))

    return make


def dataset_tree(files, fmt="tif"):
    """files: {name: volume maker}"""

    def make(folder):
        for name, vol in files.items():
            p = Path(folder) / f"{name}.{fmt}"
            if fmt == "zarr":
                synth.write_zarr(p, vol(), block_size=(8, 16, 16), layout="matlab")
            else:
                synth.write_tiff(p, vol())

    return make


def probe_values():
    """single values in (0, 1]: every exponent, dense near 1, subnormals, 1 and 0.5"""
    rng = CounterRNG(4, "b04/probe")
    a = rng.uint32(1 << 19) % np.uint32(0x3F800000) + np.uint32(1)
    b = np.uint32(0x3F000000) + rng.uint32(1 << 18) % np.uint32(0x00800000)
    c = np.uint32(0x3F7F0000) + np.arange(1 << 16, dtype=np.uint32) % np.uint32(0x00010000)
    d = np.array([0x3F800000, 0x3F000000, 1, 0x007FFFFF, 0x00800000], dtype=np.uint32)
    return np.concatenate([a, b, c, d]).view(np.float32)


# ---------------------------------------------------------------------------------------------
# comparisons with the permissive backend
# ---------------------------------------------------------------------------------------------


def _power_close(power):
    """0006's single-precision tolerance on the normalized spectrum |F| / max: largest
    difference <= 1e-5 (of the maximum, 1), relative L2 <= 1e-6; ``power`` maps a file's
    values to it (gamma 0.5: squared; log10: 10^)"""

    def fn(ref, act) -> list[str]:
        if not isinstance(ref, MArray) or not isinstance(act, MArray):
            return ["not an array"]
        if ref.mclass != act.mclass or tuple(ref.size) != tuple(act.size):
            return [f"{ref.mclass} {ref.size} against {act.mclass} {act.size}"]
        r = power(np.asarray(ref.data, np.float64))
        a = power(np.asarray(act.data, np.float64))
        d = float(np.max(np.abs(r - a))) if r.size else 0.0
        l2 = float(np.linalg.norm(r - a) / max(np.linalg.norm(r), 1e-300)) if r.size else 0.0
        out = []
        if d > 1e-5:
            out.append(f"largest difference of |F|/max {d:g} > 1e-5")
        if l2 > 1e-6:
            out.append(f"relative L2 of |F|/max {l2:g} > 1e-6")
        return out

    return fn


GAMMA_CLOSE = Custom(
    _power_close(np.square),
    "proposed (B04): MIPs and slices with the permissive FFT backend, |F|/max (their square) "
    "within 0006's single tolerance",
)
STACK_CLOSE = Custom(
    _power_close(lambda v: np.power(10.0, v)),
    "proposed (B04): the log10 spectrum with the permissive FFT backend, |F|/max (10^) within "
    "0006's single tolerance",
)


# MATLAB's single log10 is not correctly rounded (1 or 2 ulp off at 14% of the probe's values,
# as B10 found); the port's is (the probe case pins both)
LOG10_ULP = tol(reason="proposed (B04, as B10): log10 of single within 2 ulp of MATLAB's", ulp=2)


def spectrum_tree(backend, stack=False):
    if backend == "exact":
        if stack:
            return tree(per_file={"*/FFT/*.tif": LOG10_ULP}, tiff_bytes=False)
        return tree()
    per_file = {"*_MIP_*.tif": GAMMA_CLOSE, "*_center_*.tif": GAMMA_CLOSE}
    if stack:
        per_file["*/FFT/*.tif"] = STACK_CLOSE
    return tree(per_file=per_file, tiff_bytes=False)


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

_FRAME_CODE = (
    "f = [in{1}, in{2}]; o = [in{1}, 'FFT/', in{3}]; "
    "if ~exist([in{1}, 'FFT'], 'dir'), mkdir([in{1}, 'FFT']); end; "
    "XR_fftSpectrumComputingFrame(f, o, in{4 : end}); out = {};"
)


def p_frame(d, name, result, *opts, backend):
    from petakit5d import matlab

    os.makedirs(d + "FFT", exist_ok=True)
    with P.using(backend):
        matlab.XR_fftSpectrumComputingFrame(d + name, d + "FFT/" + result, *opts)
    return []


def frame_cases(cid, make, *opts, name="im.tif", result="im.tif", error=False, stack=False):
    return [
        Case(
            f"{cid}.{be}",
            code=_FRAME_CODE,
            args=[InPath("d"), name, result, *opts],
            nargout=0,
            python=p_frame,
            py_kwargs={"backend": be},
            inputs={"d": Tree(make, "d")},
            expect_error=error,
            py_error=ERRORS,
            tree=spectrum_tree(be, stack),
        )
        for be in BACKENDS
    ]


# FFTSpectrumTest.matchesWholeSpectrumComputation: 0.108, 0.3, outPixelSize 0.05
BASE = ["xyPixelSize", 0.108, "dz", 0.3, "outPixelSize", 0.05]


def _frame_cases():
    crop = frame_tree(lambda: textured(CROP, "crop"))
    pad_y = frame_tree(lambda: textured(PAD_Y, "pad_y"))
    even = frame_tree(lambda: textured(EVEN, "even"))
    big = frame_tree(lambda: textured((40, 60, 50), "big"))
    n61, n81 = [61, 61, 61], [81, 81, 81]
    cases = []
    cases += frame_cases("frame.crop", crop, *BASE, "outSize", n61)
    cases += frame_cases("frame.pad_y", pad_y, *BASE, "outSize", n61)
    cases += frame_cases("frame.pad_all", pad_y, *BASE, "outSize", n81)
    cases += frame_cases(
        "frame.stack", crop, *BASE, "outSize", n61, "save3DStack", True, stack=True
    )
    cases += frame_cases("frame.background", crop, *BASE, "outSize", n61, "background", 1200)
    cases += frame_cases("frame.cubic", crop, *BASE, "outSize", n61, "interpMethod", "cubic")
    cases += frame_cases("frame.even_equal", even, *BASE, "outSize", [65, 65, 72])
    # the read box's ends: ceil(61 ./ [2.16 2.16 6]) .* 1.1 = [31.9 31.9 12.1]
    cases += frame_cases(
        "frame.region", big, *BASE, "outSize", n61, "save3DStack", True, stack=True
    )
    cases += frame_cases(
        "frame.region.zarr",
        frame_tree(lambda: textured((40, 60, 50), "big"), "zarr"),
        *BASE,
        "outSize",
        n61,
        "zarrFile",
        True,
        name="im.zarr",
    )
    cases += frame_cases("frame.nearest", crop, *BASE, "outSize", n61, "interpMethod", "nearest")
    cases += frame_cases(
        "frame.shrink",
        big,
        "xyPixelSize",
        0.1,
        "dz",
        0.2,
        "outPixelSize",
        0.25,
        "outSize",
        [15, 21, 19],
    )
    cases += frame_cases(
        "frame.default_pixel", crop, "xyPixelSize", 0.108, "dz", 0.3, "outSize", [41, 41, 41]
    )
    cases += frame_cases(
        "frame.background_fraction", crop, *BASE, "outSize", n61, "background", 1100.37
    )
    cases += frame_cases(
        "frame.single",
        frame_tree(lambda: textured(CROP, "single", "single")),
        *BASE,
        "outSize",
        n61,
        "save3DStack",
        True,
        stack=True,
    )
    # one number: no error where the resized image has that size in every axis
    cube = frame_tree(lambda: textured((10, 30, 30), "cube"))
    cases += frame_cases("frame.scalar_size", cube, "xyPixelSize", 0.1, "dz", 0.3, "outSize", 30)
    cases += frame_cases(
        "frame.scalar_size.error",
        cube,
        "xyPixelSize",
        0.1,
        "dz",
        0.3,
        "outSize",
        31,
        error=True,
    )
    # FFTSpectrumTest.evenSizeThatNeedsCropOrPadErrors
    cases += frame_cases(
        "frame.even_error",
        frame_tree(lambda: textured((10, 20, 20), "even_error")),
        "xyPixelSize",
        0.1,
        "dz",
        0.2,
        "outPixelSize",
        0.1,
        "outSize",
        [33, 32, 33],
        error=True,
    )
    cases += frame_cases(
        "frame.one_plane_error",
        frame_tree(lambda: textured((1, 30, 30), "plane")),
        *BASE,
        "outSize",
        n61,
        error=True,
    )
    cases += frame_cases(
        "frame.existing",
        frame_tree(lambda: textured(CROP, "crop"), existing=True),
        *BASE,
        "outSize",
        n61,
    )
    return cases


_WRAP_CODE = (
    "here = pwd; try, XR_fftSpectrumComputingWrapper({in{1}}, in{2 : end}); "
    "catch ME, cd(here); rethrow(ME); end; cd(here); out = {};"
)
LOCAL = ["parseCluster", False, "masterCompute", True]


def p_wrapper(d, *opts, backend):
    from petakit5d import matlab

    with P.using(backend):
        matlab.XR_fftSpectrumComputingWrapper([d], *opts)
    return []


def wrapper_cases(cid, make, *opts, error=False, stack=False):
    return [
        Case(
            f"{cid}.{be}",
            code=_WRAP_CODE,
            args=[InPath("d"), *opts, *LOCAL],
            nargout=0,
            python=p_wrapper,
            py_kwargs={"backend": be},
            inputs={"d": Tree(make, "d")},
            expect_error=error,
            py_error=ERRORS,
            tree=spectrum_tree(be, stack),
        )
        for be in BACKENDS
    ]


SCAN = (
    "Scan_Iter_0000_0000_Cam{}_ch{}_CAM1_stack0000_488nm_0000000msec_0000001234msecAbs_"
    "000x_000y_000z_0000t"
)


def _wrapper_cases():
    cos = dataset_tree({"cosine": cosine})
    n = [N_COS] * 3
    one = ["xyPixelSize", 0.1, "dz", 0.1]
    cases = []
    # FFTSpectrumTest.spectrumOfCosineHasExpectedPeaks, save3DStackWritesLogSpectrum
    cases += wrapper_cases(
        "wrapper.cosine",
        cos,
        "resultDirName",
        "FFT",
        *one,
        "outSize",
        n,
        "outPixelSize",
        0.1,
        "save3DStack",
        False,
    )
    cases += wrapper_cases(
        "wrapper.cosine.stack",
        cos,
        *one,
        "outSize",
        n,
        "outPixelSize",
        0.1,
        "save3DStack",
        True,
        stack=True,
    )
    # resamplesAnisotropicVoxels (outPixelSize empty: min(xyPixelSize, dz))
    aniso = dataset_tree({"aniso": lambda: textured((40, 80, 80), "aniso", "single")})
    cases += wrapper_cases("wrapper.aniso", aniso, "xyPixelSize", 0.1, "dz", 0.2, "outSize", n)
    chans = dataset_tree(
        {
            SCAN.format("A", 0): lambda: textured(CROP, "ca0"),
            SCAN.format("B", 0): lambda: textured(CROP, "cb0"),
            SCAN.format("A", 1): lambda: textured(CROP, "ca1"),
        }
    )
    cases += wrapper_cases(
        "wrapper.channels",
        chans,
        *BASE,
        "outSize",
        [41, 41, 41],
        "channelPatterns",
        ["CamA_ch0", "CamB"],
        "resultDirName",
        "spectra",
    )
    cases += wrapper_cases(
        "wrapper.channel_char",
        chans,
        *BASE,
        "outSize",
        [41, 41, 41],
        "channelPatterns",
        "CamA",
        "background",
        1000.123456789,
    )
    cases += wrapper_cases(
        "wrapper.zarr",
        dataset_tree(
            {"z1": lambda: textured(CROP, "z1"), "z2": lambda: textured(PAD_Y, "z2")}, "zarr"
        ),
        *BASE,
        "outSize",
        [41, 41, 41],
        "zarrFile",
        True,
    )
    cases += wrapper_cases(
        "wrapper.digits",
        dataset_tree({"im": lambda: textured(CROP, "digits")}),
        "xyPixelSize",
        0.1083333333333333,
        "dz",
        0.29999999999999993,
        "outPixelSize",
        0.0512345678901234,
        "outSize",
        [41, 43, 45],
        "interpMethod",
        "cubic",
    )
    cases += wrapper_cases(
        "wrapper.no_file",
        dataset_tree({"im": lambda: textured(CROP, "crop")}),
        *BASE,
        "channelPatterns",
        ["CamC"],
        error=True,
    )
    return cases


def p_probe(x):
    return [fs._log10_single(x, 1), np.sqrt(x)]


# MATLAB's log10 and .^ 0.5 of single data (the outputs' arithmetic)
PROBES = [
    Case(
        "probe.log10_sqrt_single",
        code="x = in{1}; out = {log10(x), x .^ 0.5};",
        args=[In("x")],
        nargout=2,
        python=p_probe,
        inputs={"x": Value(probe_values)},
        compare={0: LOG10_ULP},
    )
]

CASES = PROBES + _frame_cases() + _wrapper_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# FFTSpectrumTest in intent, on (z, y, x) data
# ---------------------------------------------------------------------------------------------

C = (N_COS - 1) // 2  # the DC voxel (0-based)


def run_dataset(folder, **kw):
    fs.fft_spectrum_dataset(str(folder), parse_cluster=False, master_compute=True, **kw)


def write_cosine(folder):
    synth.write_tiff(Path(folder) / "cosine.tif", cosine())


def test_spectrum_of_cosine_has_expected_peaks(tmp_path):
    from petakit5d.io import read_tiff

    write_cosine(tmp_path)
    run_dataset(
        tmp_path,
        result_dir_name="FFT",
        xy_pixel_size=0.1,
        dz=0.1,
        out_size=(N_COS,) * 3,
        out_pixel_size=0.1,
        save_3d_stack=False,
    )
    fft = tmp_path / "FFT"
    # without save3DStack, the stack is an empty placeholder
    assert (fft / "cosine.tif").stat().st_size == 0
    xy = read_tiff(fft / "central_slices" / "cosine_center_xy.tif").astype(np.float64)
    assert xy.shape == (N_COS, N_COS)
    # gamma 0.5 of the spectrum normalized to DC: the cosine's peaks at +-K are sqrt(0.25)
    assert abs(xy[C, C] - 1) <= 1e-5
    assert abs(xy[C, C + K_COS] - 0.5) <= 1e-4
    assert abs(xy[C, C - K_COS] - 0.5) <= 1e-4
    others = np.ones(xy.shape, bool)
    others[C, [C - K_COS, C, C + K_COS]] = False
    assert xy[others].max() < 0.01
    # point symmetric about DC
    np.testing.assert_allclose(xy, xy[::-1, ::-1], rtol=1e-5, atol=1e-7)
    for axis in "yxz":
        mip = read_tiff(fft / "MIPs" / f"cosine_MIP_{axis}.tif")
        assert mip.shape == (N_COS, N_COS)
        assert abs(float(mip.max()) - 1) <= 1e-5
    for plane in ("yz", "xz"):
        assert (fft / "central_slices" / f"cosine_center_{plane}.tif").is_file()


def test_save_3d_stack_writes_log_spectrum(tmp_path):
    from petakit5d.io import read_tiff

    write_cosine(tmp_path)
    run_dataset(
        tmp_path,
        xy_pixel_size=0.1,
        dz=0.1,
        out_size=(N_COS,) * 3,
        out_pixel_size=0.1,
        save_3d_stack=True,
    )
    spectrum = read_tiff(tmp_path / "FFT" / "cosine.tif")
    assert spectrum.shape == (N_COS,) * 3
    assert spectrum.dtype == np.float32
    assert abs(float(spectrum[C, C, C])) <= 1e-5
    assert abs(float(spectrum[C, C, C + K_COS]) - np.log10(0.25)) <= 1e-4


def test_resamples_anisotropic_voxels(tmp_path):
    from petakit5d.io import read_tiff

    vol = synth.texture((40, 80, 80), 2, rng=CounterRNG(4, "b04/aniso_py"))
    synth.write_tiff(tmp_path / "aniso.tif", vol)
    run_dataset(tmp_path, xy_pixel_size=0.1, dz=0.2, out_size=(N_COS,) * 3)
    xz = read_tiff(tmp_path / "FFT" / "central_slices" / "aniso_center_xz.tif")
    assert xz.shape == (N_COS, N_COS)


def reference_spectrum(im, xy, dz, px, n, method, background):
    """FFTSpectrumTest.referenceSpectrum: the whole image resized (imresize3), padded or
    cropped to n (z, y, x) and the shifted |FFT| normalized; (MIPs y, x, z; slices xy, yz, xz;
    the log10 stack) as fft_spectrum_file writes them"""
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat.fft import fftn
    from petakit5d.matlabcompat.resize import imresize3

    a = im.astype(np.float32)
    if background:
        a = dtypes.minus(a, float(background))
        a = a * (a > 0)
    scale = np.array([dz, xy, xy]) / px
    rsz = dtypes.mround(np.array(a.shape, np.float64) * scale)
    a = imresize3(a, output_size=[int(v) for v in rsz], method=method)
    n = np.asarray(n)
    sz = np.array(a.shape)
    if np.any(sz != n):
        if np.any(sz < n):
            pre = np.maximum(0, (n - sz) // 2)
            post = np.maximum(0, -(-(n - sz) // 2))
            a = np.pad(a, list(zip(pre, post, strict=True)))
        c = dtypes.mround((np.array(a.shape) + 1) / 2).astype(int) - 1
        h = (n - 1) // 2
        a = a[tuple(slice(ci - hi, ci + hi + 1) for ci, hi in zip(c, h, strict=True))]
    s = np.fft.fftshift(fs._abs(fftn(a), 1))
    q = s / s.max()
    stack = fs._log10_single(q, 1)
    g = np.sqrt(q)
    nz, ny, nx = g.shape
    mips = [g.max(axis=1).T, g.max(axis=2).T, g.max(axis=0)]
    # (MATLAB's round((n + 1) / 2), 0-based: n // 2, where fftshift puts DC)
    slices = [g[nz // 2], g[:, :, nx // 2].T, g[:, ny // 2, :].T]
    return mips, slices, stack


SPECTRUM_CASES = {
    "crop": (CROP, (61, 61, 61), {}),
    "padY": (PAD_Y, (61, 61, 61), {}),
    "padAll": (PAD_Y, (81, 81, 81), {}),
    "stack": (CROP, (61, 61, 61), {"save_3d_stack": True}),
    "background": (CROP, (61, 61, 61), {"background": 1200}),
    "cubic": (CROP, (61, 61, 61), {"interp_method": "cubic"}),
    "evenEqual": (EVEN, (72, 65, 65), {}),
}


@pytest.mark.parametrize("cid", list(SPECTRUM_CASES))
def test_matches_whole_spectrum_computation(cid, tmp_path, monkeypatch):
    from petakit5d.io import read_tiff
    from petakit5d.matlabcompat import resize as R

    sz, n, opts = SPECTRUM_CASES[cid]
    im = textured(sz, f"py/{cid}")
    synth.write_tiff(tmp_path / "im.tif", im)
    mips, slices, stack = reference_spectrum(
        im, 0.108, 0.3, 0.05, n, opts.get("interp_method", "linear"), opts.get("background", 0)
    )

    # the resize must not fall back to imresize3 of the whole image
    def no_imresize3(*a, **k):
        raise AssertionError("imresize3 called")

    monkeypatch.setattr(R, "imresize3", no_imresize3)
    (tmp_path / "FFT").mkdir()
    fs.fft_spectrum_file(
        str(tmp_path / "im.tif"),
        str(tmp_path / "FFT" / "im.tif"),
        xy_pixel_size=0.108,
        dz=0.3,
        out_pixel_size=0.05,
        out_size=n,
        **opts,
    )
    for k, (axis, plane) in enumerate(zip("yxz", ("xy", "yz", "xz"), strict=True)):
        got = read_tiff(tmp_path / "FFT" / "MIPs" / f"im_MIP_{axis}.tif")
        np.testing.assert_array_equal(got, mips[k], err_msg=f"MIP along {axis}")
        got = read_tiff(tmp_path / "FFT" / "central_slices" / f"im_center_{plane}.tif")
        np.testing.assert_array_equal(got, slices[k], err_msg=f"central {plane} slice")
    if opts.get("save_3d_stack"):
        np.testing.assert_array_equal(read_tiff(tmp_path / "FFT" / "im.tif"), stack)


def test_even_size_that_needs_crop_or_pad_errors(tmp_path):
    synth.write_tiff(tmp_path / "im.tif", textured((10, 20, 20), "py/even"))
    (tmp_path / "FFT").mkdir()
    with pytest.raises(fs.FFTSpectrumError, match="three odd numbers"):
        fs.fft_spectrum_file(
            str(tmp_path / "im.tif"),
            str(tmp_path / "FFT" / "im.tif"),
            xy_pixel_size=0.1,
            dz=0.2,
            out_pixel_size=0.1,
            out_size=(33, 33, 32),
        )
    assert not (tmp_path / "FFT" / "im.tif").exists()


# ---------------------------------------------------------------------------------------------
# unit tests
# ---------------------------------------------------------------------------------------------


def test_the_read_box_truncates_non_integer_ends():
    # keep [y x z] = min(size, ceil(61 ./ [2.16 2.16 6]) .* 1.1) = [31.9 31.9 12.1]
    assert fs.spectrum_region((40, 60, 50), 0.108, 0.3, 0.05, 61) == (13, 14, 9, 25, 45, 40)
    # a small image is read whole
    assert fs.spectrum_region(CROP, 0.108, 0.3, 0.05, (61, 61, 61)) == (0, 0, 0, 12, 31, 30)
    # without an output pixel size: min(xy, dz)
    box = fs.spectrum_region((40, 60, 50), 0.1, 0.3, None, 21)
    assert box == fs.spectrum_region((40, 60, 50), 0.1, 0.3, 0.1, 21)


def test_log10_single_is_correctly_rounded():
    from decimal import Context, Decimal

    x = probe_values()[::997]
    got = fs._log10_single(x, 2)
    ctx = Context(prec=40)
    for xi, gi in zip(x.tolist(), got.tolist(), strict=True):
        v = ctx.log10(Decimal(xi))
        g = np.float32(gi)
        # no single nearer to log10(x) on either side
        for nb in (np.nextafter(g, np.float32(-np.inf)), np.nextafter(g, np.float32(np.inf))):
            assert abs(Decimal(float(g)) - v) <= abs(Decimal(float(nb)) - v)
    assert fs._log10_single(np.array([1, 0], np.float32), 1).tolist() == [0.0, -np.inf]


def test_log10_ties_are_decided_in_40_digits():
    from decimal import Context, Decimal

    x = 0.7
    r = np.float32(np.log10(x))
    for nb in (np.nextafter(r, np.float32(np.inf)), np.nextafter(r, np.float32(-np.inf))):
        mid = (float(r) + float(nb)) / 2
        v = Context(prec=40).log10(Decimal(x))
        want = min(r, nb) if v < Decimal(mid) else max(r, nb)
        assert fs._log10_decided(x, mid, r, nb) == want


def test_images_have_matlabs_orientation():
    a = np.arange(3 * 5 * 7, dtype=np.float32).reshape(3, 5, 7)
    # squeeze(max(im, [], 1)) of MATLAB's ny x nx x nz is nx x nz
    assert fs._image(a.max(axis=1, keepdims=True), -2).shape == (7, 3)
    assert fs._image(a.max(axis=2, keepdims=True), -1).shape == (5, 3)
    assert fs._image(a.max(axis=0, keepdims=True), -3).shape == (5, 7)
    np.testing.assert_array_equal(fs._image(a[:, 0:1, :], -2), a[:, 0, :].T)


def test_an_existing_result_is_kept(tmp_path):
    synth.write_tiff(tmp_path / "im.tif", textured(CROP, "py/existing"))
    (tmp_path / "FFT").mkdir()
    (tmp_path / "FFT" / "im.tif").write_bytes(b"x")
    fs.fft_spectrum_file(
        str(tmp_path / "im.tif"), str(tmp_path / "FFT" / "im.tif"), out_size=(21, 21, 21)
    )
    assert (tmp_path / "FFT" / "im.tif").read_bytes() == b"x"
    assert not (tmp_path / "FFT" / "MIPs").exists()


def test_one_plane_errors(tmp_path):
    synth.write_tiff(tmp_path / "im.tif", textured((1, 30, 30), "py/plane"))
    (tmp_path / "FFT").mkdir()
    with pytest.raises(fs.FFTSpectrumError) as err:
        fs.fft_spectrum_file(
            str(tmp_path / "im.tif"), str(tmp_path / "FFT" / "im.tif"), out_size=(21,) * 3
        )
    assert err.value.identifier == "MATLAB:sizeDimensionsMustMatch"


def test_the_tasks_get_matlabs_text(tmp_path, monkeypatch, capsys):
    from petakit5d.jobs import framework

    seen = {}

    def fake(inputs, outputs, tasks, **kw):
        seen.update(inputs=inputs, outputs=outputs, tasks=tasks, kw=kw)
        return [True] * len(tasks)

    monkeypatch.setattr(framework, "run_tasks", fake)
    synth.write_tiff(tmp_path / "a_CamA.tif", textured(CROP, "py/t1"))
    synth.write_tiff(tmp_path / "b_CamB.tif", textured(CROP, "py/t2"))
    fs.fft_spectrum_dataset(
        str(tmp_path),
        xy_pixel_size=0.1083333333333333,
        dz=0.3,
        channel_patterns="CamB",
        background=1000.123456789,
        out_size=(41, 43, 45),
    )
    assert "Out pixel size: 1.083333e-01 um" in capsys.readouterr().out
    # the folder as the code gives it (with '/' on Windows too), joined with '/'
    posix = tmp_path.as_posix()
    assert [p.replace("\\", "/") for p in seen["inputs"]] == [f"{posix}/b_CamB.tif"]
    assert [p.replace("\\", "/") for p in seen["outputs"]] == [f"{posix}/FFT/b_CamB.tif"]
    kw = seen["tasks"][0].kwargs
    # %0.20d of a non-integer is %e with 20 digits (exact); %d of the background %e (7 digits)
    assert kw["xy_pixel_size"] == 0.1083333333333333
    assert kw["out_pixel_size"] == 0.1083333333333333
    assert kw["background"] == 1000.123
    assert kw["out_size"] == [41.0, 43.0, 45.0]
    assert seen["kw"]["mem_allocate"] == pytest.approx(45 * 41 * 43 * 4 / 1024**3 * 15)
    assert (tmp_path / "FFT" / "parameters.mat").is_file()


def test_no_matching_file_errors(tmp_path):
    synth.write_tiff(tmp_path / "a.tif", textured(CROP, "py/none"))
    with pytest.raises(fs.FFTSpectrumError) as err:
        fs.fft_spectrum_dataset(str(tmp_path), channel_patterns=["CamC"], parse_cluster=False)
    assert err.value.identifier == "MATLAB:badsubscript"
