"""Chromatic shift estimation and correction (petakit5d.tools.chromatic, task B03).

Oracle cases against MATLAB: XR_chromatic_shift_correction_data_wrapper on
ChromaticShiftCorrectionTest's data (two channels cut from one texture with the demo's offsets
[0 0 0; 1 -5 -2]) in the 'valid', 'same' and 'full' modes (and in upper case), a 'same' origin, a
pad value, zarr results, zarr input, large zarr files cropped in batches, two datasets, three
channels; offsets measured from bead PSFs (ChromaticShiftCorrectionTest's beads, and with a crop
and bounds that differ along y and x), PSFs with zarr data; MATLAB's errors (another mode, no
offsets, too few and too many offset rows, a pad value of large zarr files, channel and PSF
counts).
XR_chromatic_shift_estimation directly: three PSFs, a crop, one number as the bounds, the result
file, too few PSFs. Cases with PSFs run on both CPU FFT backends (``.exact``: MATLAB's bits;
``.permissive``: the peaks of the NCC within 1e-4, decision record 0006, the shifts exact).
Results are compared byte for byte (tiff), by data and metadata (zarr) and by content (.mat).

Then ChromaticShiftCorrectionTest in intent on the port's (z, y, x) data, the regions of the
modes against MATLAB's formulas, the bounds in (z, y, x), the MATLAB-name layer's conversions and
the job shape (a framework spy).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.core.dtypes import mround
from petakit5d.io import read_mat, read_tiff, read_zarr
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, InPath, Tree
from petakit5d.testing.compare import tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import chromatic as ch

BACKENDS = ("exact", "permissive")
ERRORS = (ch.ChromaticError, ValueError, OSError, RuntimeError, IndexError)
# ChromaticShiftCorrectionTest: ImSize [48 40 24] (MATLAB's [y x z]) as (z, y, x)
IM = (24, 48, 40)
# the demo's offsets, MATLAB's rows [y x z] and the port's (z, y, x)
OFFSET_M = np.array([[0.0, 0.0, 0.0], [1.0, -5.0, -2.0]])
OFFSET = OFFSET_M[:, [2, 0, 1]]
FSN = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_"
    "000x_000y_000z_0000t"
)
RESULT = "Chromatic_Shift_Corrected"
LOCAL = ["parseCluster", False, "masterCompute", True]
LOCAL_KW = {"parse_cluster": False, "master_compute": True}
PSF_NAMES = ("PSF/psf_CamA.tif", "PSF/psf_CamB.tif", "PSF/psf_CamC.tif")


# ---------------------------------------------------------------------------------------------
# data: (z, y, x) arrays from counter-based random numbers (bit-stable everywhere)
# ---------------------------------------------------------------------------------------------


def scene_of(size, name):
    """pktest.texture(size, 1.5, [500 5000]) as uint16 (rounded as MATLAB's uint16 rounds)"""
    t = synth.texture(size, 1.5, (500.0, 5000.0), CounterRNG(6, f"b03/{name}"))
    return np.asarray(mround(t.astype(np.float64)), dtype=np.uint16)


def shifted_channels(offsets, im=IM, name="scene"):
    """ChromaticShiftCorrectionTest.shiftedChannels: channel c voxel i is scene voxel
    i + offsets[c] + base; returns the channels, the scene and base"""
    o = np.asarray(offsets, dtype=np.int64)
    base = np.maximum(0, -o.min(axis=0))
    scene = scene_of(tuple(np.asarray(im) + o.max(axis=0) + base), name)
    chans = []
    for c in range(o.shape[0]):
        s = o[c] + base
        chans.append(
            np.ascontiguousarray(
                scene[s[0] : s[0] + im[0], s[1] : s[1] + im[1], s[2] : s[2] + im[2]]
            )
        )
    return chans, scene, base


def bead_psf(offset, size=(24, 40, 40), bead=(11, 19, 19)):
    """ChromaticShiftCorrectionTest's PSF: a bead at ``bead - offset`` (0-based (z, y, x)),
    sigma 1.5 (y, x) and 2 (z), 2000 over 100, as uint16"""
    pos = np.asarray(bead, dtype=np.float64) - np.asarray(offset, dtype=np.float64)
    v = synth.beads(size, pos.reshape(1, 3), (2.0, 1.5), 2000.0, 100.0)
    return np.asarray(mround(v.astype(np.float64)), dtype=np.uint16)


def dataset(offsets=OFFSET, cams=("CamA", "CamB"), fmt="tif", psfs=False, name="scene"):
    """a folder maker: the channels as <FSN>.tif (or MATLAB-layout zarr), bead PSFs in PSF/"""

    def make(folder):
        folder = Path(folder)
        chans, _, _ = shifted_channels(offsets, name=name)
        for c, cam in enumerate(cams):
            p = folder / (FSN.format(cam) + "." + fmt)
            if fmt == "zarr":
                synth.write_zarr(p, chans[c], block_size=(8, 16, 16), layout="matlab")
            else:
                synth.write_tiff(p, chans[c])
        if psfs:
            (folder / "PSF").mkdir()
            for c in range(len(cams)):
                synth.write_tiff(folder / PSF_NAMES[c], bead_psf(np.asarray(offsets)[c]))

    return make


def two_datasets(make_a, make_b):
    def make(folder):
        for sub, mk in (("a", make_a), ("b", make_b)):
            (Path(folder) / sub).mkdir()
            mk(Path(folder) / sub)

    return make


def psfs_only(offsets, n):
    def make(folder):
        for c in range(n):
            synth.write_tiff(Path(folder) / f"psf_{c}.tif", bead_psf(np.asarray(offsets)[c]))

    return make


# ---------------------------------------------------------------------------------------------
# oracle cases: XR_chromatic_shift_correction_data_wrapper
# ---------------------------------------------------------------------------------------------

# in{1}: the data folder; in{2}: its subfolders (empty: the folder); in{3}: PSF names in it
_WRAP_CODE = (
    "d = in{1}; s = in{2}; p = in{3}; o = in(4 : end); "
    "if isempty(s), ds = {d}; else, ds = cellfun(@(x) [d, x], s, 'UniformOutput', false); end; "
    "if ~isempty(p), p = cellfun(@(x) [d, x], p, 'UniformOutput', false); "
    "o = [o, {'psfFullpaths', p}]; end; "
    "here = pwd; try, XR_chromatic_shift_correction_data_wrapper(ds, o{:}); "
    "catch ME, cd(here); rethrow(ME); end; cd(here); out = {};"
)


def p_wrapper(d, subs, psf, *opts, backend="permissive"):
    from petakit5d import matlab

    ds = [d + s for s in subs] if subs else [d]
    if psf:
        opts = (*opts, "psfFullpaths", [d + x for x in psf])
    with P.using(backend):
        matlab.XR_chromatic_shift_correction_data_wrapper(ds, *opts)
    return []


# the NCC's peak with the permissive backend (0006: the NCC to 5e-3 in K06's cases; these
# peaks are within 1e-6), the shifts exact
PEAK_CLOSE = tol(
    reason="proposed (B03, as K06): the NCC peaks with the permissive FFT backend", abs=1e-4
)


def chroma_tree(backend):
    if backend == "permissive":
        return tree(per_file={"*.mat": PEAK_CLOSE})
    return tree()


# large zarr files: cropped in batches by XR_crop_zarr
LARGE_ZARR = ["zarrFile", True, "largeFile", True, "batchSize", [16, 16, 8], "blockSize", [8, 8, 8]]


def wrapper_case(cid, make, *opts, subs=(), psf=(), error=False, backend=None):
    kw = {} if backend is None else {"backend": backend}
    return Case(
        cid if backend is None else f"{cid}.{backend}",
        code=_WRAP_CODE,
        args=[InPath("d"), list(subs), list(psf), *opts, *LOCAL],
        nargout=0,
        python=p_wrapper,
        py_kwargs=kw,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
        tree=chroma_tree(backend or "exact"),
    )


def psf_cases(cid, make, *opts, psf=PSF_NAMES[:2], error=False):
    return [wrapper_case(cid, make, *opts, psf=psf, error=error, backend=b) for b in BACKENDS]


def _wrapper_cases():
    data = dataset()
    off = ["chromaticOffset", OFFSET_M]
    cases = [
        # ChromaticShiftCorrectionTest.validModeAlignsChannels
        wrapper_case(
            "wrapper.valid", data, "mode", "valid", *off, "channelPatterns", ["CamA", "CamB"]
        ),
        # sameModeKeepsImageSize, fullModeContainsBothChannels
        wrapper_case("wrapper.same", data, "mode", "same", *off),
        wrapper_case("wrapper.full", data, "mode", "full", *off),
        # zarrOutput
        wrapper_case("wrapper.save_zarr", data, *off, "saveZarr", True, "blockSize", [16, 16, 8]),
        wrapper_case("wrapper.same.origin", data, "mode", "same", *off, "newOrigin", [3, -2, 2]),
        wrapper_case("wrapper.full.pad", data, "mode", "full", *off, "padValue", 7),
        wrapper_case("wrapper.full.pad_fraction", data, "mode", "full", *off, "padValue", 7.6),
        # %d gives 7.500000e+00: the tasks pad with 8, not 7
        wrapper_case("wrapper.full.pad_text", data, "mode", "full", *off, "padValue", 7.4999999),
        # (another mode is refused; it cropped the channels shifted against their offsets
        # before 45b2bb1)
        wrapper_case("wrapper.other_mode.error", data, "mode", "shift", *off, error=True),
        wrapper_case("wrapper.upper_mode", data, "mode", "Valid", *off),
        wrapper_case("wrapper.zarr_in", dataset(fmt="zarr"), *off, "zarrFile", True),
        wrapper_case(
            "wrapper.zarr_in.same_zarr_out",
            dataset(fmt="zarr"),
            "mode",
            "same",
            *off,
            "zarrFile",
            True,
            "saveZarr",
            True,
        ),
        wrapper_case(
            "wrapper.two_datasets",
            two_datasets(dataset(), dataset(name="scene2")),
            "mode",
            "full",
            *off,
            "resultDirName",
            "aligned",
            subs=("a/", "b"),
        ),
        wrapper_case(
            "wrapper.three_channels",
            dataset(np.array([[0, 0, 0], [-2, 1, -5], [1, -3, 2]]), ("CamA", "CamB", "CamC")),
            "mode",
            "valid",
            "chromaticOffset",
            np.array([[0.0, 0, 0], [1, -5, -2], [-3, 2, 1]]),
            "channelPatterns",
            ["CamA", "CamB", "CamC"],
        ),
        wrapper_case("wrapper.large_file", data, *off, "largeFile", True),
        # large zarr files are cropped in batches into zarr (only sized the jobs before 45b2bb1)
        *[
            wrapper_case(
                f"wrapper.large_zarr.{m}",
                dataset(fmt="zarr"),
                "mode",
                m,
                *off,
                *LARGE_ZARR,
            )
            for m in ("valid", "full")
        ],
        wrapper_case(
            "wrapper.large_zarr.pad_value.error",
            dataset(fmt="zarr"),
            "mode",
            "full",
            *off,
            *LARGE_ZARR,
            "padValue",
            7,
            error=True,
        ),
        # MATLAB's errors (a row per channel: ...:noOffsets, ...:offsetSize)
        wrapper_case("wrapper.no_offset.error", data, error=True),
        wrapper_case("wrapper.one_row.error", data, "chromaticOffset", [1, 2, 3], error=True),
        wrapper_case(
            "wrapper.extra_row.error",
            data,
            "chromaticOffset",
            np.concatenate([OFFSET_M, [[0.0, 9.0, 0.0]]]),
            error=True,
        ),
        wrapper_case(
            "wrapper.channels.error",
            data,
            *off,
            "channelPatterns",
            ["CamA", "CamB", "CamC"],
            error=True,
        ),
    ]
    psf_data = dataset(psfs=True)
    # offsetsEstimatedFromPSFs
    cases += psf_cases("wrapper.psf", psf_data, "maxOffset", [8, 8, 6])
    # bounds of 3 along x and 8 along y: x = -5 is out of reach (MATLAB's swapped bounds found it
    # before 45b2bb1)
    cases += psf_cases(
        "wrapper.psf.bounds",
        psf_data,
        "mode",
        "full",
        "maxOffset",
        [8, 3, 6],
        "cropLength",
        [2, 3, 1],
    )
    cases += psf_cases("wrapper.psf.default_bounds", psf_data, "mode", "same")
    # an origin places the output by the offsets less their minimum
    cases += psf_cases(
        "wrapper.psf.origin", psf_data, "mode", "same", "newOrigin", [2, 3, 1], "maxOffset", 8
    )
    cases += psf_cases(
        "wrapper.psf.other_mode.error", psf_data, "mode", "shift", "maxOffset", 8, error=True
    )
    # the PSFs are read as tiff, the zarr data as zarr (zarrFile stayed false for the data after
    # the measurement, and its tasks failed, before 45b2bb1)
    cases += psf_cases(
        "wrapper.psf.zarr_data",
        dataset(fmt="zarr", psfs=True),
        "maxOffset",
        [8, 8, 6],
        "zarrFile",
        True,
    )
    cases += psf_cases("wrapper.psf.count.error", psf_data, psf=PSF_NAMES[:1], error=True)
    return cases


# ---------------------------------------------------------------------------------------------
# oracle cases: XR_chromatic_shift_estimation
# ---------------------------------------------------------------------------------------------

_EST_CODE = (
    "d = in{1}; p = cellfun(@(x) [d, x], in{2}, 'UniformOutput', false); "
    "r = in{3}; if ~isempty(r), r = [d, r]; end; "
    "out = {XR_chromatic_shift_estimation(p, 'resultFullpath', r, in{4 : end})};"
)


def p_estimation(d, names, result, *opts, backend):
    from petakit5d import matlab

    r = d + result if result else ""
    with P.using(backend):
        return matlab.XR_chromatic_shift_estimation(
            [d + x for x in names], "resultFullpath", r, *opts
        )


THREE = np.array([[0, 0, 0], [-2, 1, -5], [3, -2, 4]])


def est_cases(cid, make, names, result, *opts, error=False):
    return [
        Case(
            f"{cid}.{b}",
            code=_EST_CODE,
            args=[InPath("d"), list(names), result, *opts],
            nargout=1,
            python=p_estimation,
            py_kwargs={"backend": b},
            inputs={"d": Tree(make, "d")},
            expect_error=error,
            py_error=ERRORS,
            tree=chroma_tree(b),
            **({"compare": PEAK_CLOSE} if b == "permissive" else {}),
        )
        for b in BACKENDS
    ]


def _estimation_cases():
    names = [f"psf_{c}.tif" for c in range(3)]
    three = psfs_only(THREE, 3)
    cases = []
    cases += est_cases("est.three", three, names, "m.mat", "maxOffset", [6, 6, 6])
    cases += est_cases(
        "est.crop", three, names, "", "maxOffset", [7, 7, 6], "cropLength", [3, 2, 2]
    )
    cases += est_cases("est.scalar_bounds", three, names[:2], "m.mat", "maxOffset", 6)
    cases += est_cases("est.default", three, names[:2], "")
    cases += est_cases("est.one.error", three, names[:1], "", error=True)
    return cases


CASES = _wrapper_cases() + _estimation_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# ChromaticShiftCorrectionTest in intent, on (z, y, x) data
# ---------------------------------------------------------------------------------------------


def write_dataset(folder, offsets=OFFSET, fmt="tif", psfs=False):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    dataset(offsets, fmt=fmt, psfs=psfs)(folder)
    return str(folder) + "/"


def run(path, mode="valid", **kw):
    kw.setdefault("chromatic_offset", OFFSET)
    ch.chromatic_shift_correction_dataset(
        path, mode=mode, channel_patterns=("CamA", "CamB"), **LOCAL_KW, **kw
    )


def outputs(path, ext="tif"):
    rd = read_zarr if ext == "zarr" else read_tiff
    return [rd(f"{path}{RESULT}/{FSN.format(c)}.{ext}") for c in ("CamA", "CamB")]


def test_valid_mode_aligns_channels(tmp_path):
    path = write_dataset(tmp_path / "d")
    run(path, "valid")
    a, b = outputs(path)
    _, scene, base = shifted_channels(OFFSET)
    s = OFFSET.max(axis=0).astype(int)
    t = (np.asarray(IM) + OFFSET).min(axis=0).astype(int)
    assert a.shape == tuple(t - s)
    np.testing.assert_array_equal(b, a)
    s, t = s + base, t + base
    np.testing.assert_array_equal(a, scene[s[0] : t[0], s[1] : t[1], s[2] : t[2]])


def test_offsets_estimated_from_psfs(tmp_path):
    path = write_dataset(tmp_path / "d", psfs=True)
    run(
        path,
        "valid",
        chromatic_offset=None,
        psf_fullpaths=[path + n for n in PSF_NAMES[:2]],
        max_offset=(6, 8, 8),
    )
    a, b = outputs(path)
    np.testing.assert_array_equal(b, a)
    m = read_mat(f"{path}{RESULT}/measurement/chromatic_shift_measurement_result.mat")
    np.testing.assert_array_equal(m["relative_shift_yxz_mat"], OFFSET_M)
    assert m["max_xcorr_mat"].shape == (2, 1) and m["max_xcorr_mat"][1, 0] > 0.5


def test_mem_allocate_is_in_gigabytes(tmp_path, monkeypatch):
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs_, tasks, **kw):
        calls.append({"tasks": list(tasks), **kw})
        return [True] * len(list(outputs_))

    monkeypatch.setattr(framework, "run_tasks", fake)
    path = write_dataset(tmp_path / "d")
    run(path, "valid", job_log_dir="logs")
    assert len(calls) == 1
    assert calls[0]["mem_allocate"] == pytest.approx(48 * 40 * 24 * 4 * 2.5 / 1024**3, rel=1e-12)
    assert calls[0]["job_log_dir"] == "logs"
    assert len(calls[0]["tasks"]) == 2
    calls.clear()
    # (a large tiff file is cropped in memory: the image's memory)
    run(path, "valid", large_file=True, batch_size=(64, 128, 256))
    assert calls[0]["mem_allocate"] == pytest.approx(48 * 40 * 24 * 4 * 2.5 / 1024**3, rel=1e-12)
    assert "large_file" not in calls[0]["tasks"][0].kwargs


def test_failed_tasks_are_retried_with_twice_the_memory(tmp_path, monkeypatch):
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs_, tasks, **kw):
        calls.append(kw["mem_allocate"])
        return [False] * len(list(outputs_))

    monkeypatch.setattr(framework, "run_tasks", fake)
    path = write_dataset(tmp_path / "d")
    with pytest.warns(UserWarning, match="Some files are not processed!"):
        run(path, "valid")
    assert calls == [pytest.approx(48 * 40 * 24 * 4 * 2.5 / 1024**3), calls[0] * 2]


def test_zarr_output(tmp_path):
    path = write_dataset(tmp_path / "d")
    run(path, "valid", save_zarr=True)
    a, b = outputs(path, "zarr")
    np.testing.assert_array_equal(b, a)
    run(path, "valid", result_dir_name="tif")
    np.testing.assert_array_equal(a, read_tiff(f"{path}tif/{FSN.format('CamA')}.tif"))


def test_same_mode_keeps_image_size(tmp_path):
    path = write_dataset(tmp_path / "d")
    run(path, "same")
    a, b = outputs(path)
    assert a.shape == IM and b.shape == IM
    both = (a > 0) & (b > 0)
    assert both.sum() > 0
    np.testing.assert_array_equal(b[both], a[both])


def test_full_mode_contains_both_channels(tmp_path):
    path = write_dataset(tmp_path / "d")
    run(path, "full")
    a, b = outputs(path)
    union = tuple(((np.asarray(IM) + OFFSET).max(axis=0) - OFFSET.min(axis=0)).astype(int))
    assert a.shape == union and b.shape == union
    both = (a > 0) & (b > 0)
    np.testing.assert_array_equal(b[both], a[both])
    assert np.count_nonzero(a) == np.prod(IM)


# ---------------------------------------------------------------------------------------------
# unit tests
# ---------------------------------------------------------------------------------------------


def matlab_boxes(im_size_m, off_m, mode, origin_m=None):
    """MATLAB's formulas of XR_chromatic_shift_correction_data_wrapper (1-based, [y x z])"""
    n = np.asarray(im_size_m, dtype=np.float64)
    o = np.asarray(off_m, dtype=np.float64)
    b = np.concatenate([1 + o, n + o], axis=1)
    if mode == "valid":
        ov = np.concatenate([b[:, :3].max(0), b[:, 3:].min(0)])
        return ov - np.concatenate([b[:, :3], b[:, :3]], axis=1) + 1
    if mode == "same":
        if origin_m is None:
            ov = np.concatenate([b[:, :3].max(0), b[:, 3:].min(0)])
            origin_m = ov[:3] + np.floor((ov[3:] - ov[:3] + 1 - n) / 2)
        return np.concatenate([origin_m - b[:, :3] + 1, origin_m - b[:, :3] + n], axis=1)
    if mode == "full":
        origin = b[:, :3].min(0)
        full = b[:, 3:].max(0) - origin + 1
        return np.concatenate([origin - b[:, :3] + 1, origin - b[:, :3] + full], axis=1)
    return b


@pytest.mark.parametrize("mode", ["valid", "same", "full", "other"])
@pytest.mark.parametrize(
    "off_m",
    [OFFSET_M, np.array([[2.0, 0, -1], [-3, 4, 5], [0, -7, 1]]), np.zeros((2, 3))],
    ids=["demo", "three", "zero"],
)
def test_boxes_follow_matlab(mode, off_m):
    from petakit5d.core.bbox import box_from_matlab

    im_m = (48, 41, 23)
    boxes, pad = ch.correction_boxes((23, 48, 41), off_m[:, [2, 0, 1]], mode)
    want = matlab_boxes(im_m, off_m, mode)
    assert pad == (mode != "valid")
    for r in range(off_m.shape[0]):
        assert tuple(boxes[r]) == tuple(box_from_matlab(want[r]))


def test_same_mode_origin():
    from petakit5d.core.bbox import box_from_matlab

    boxes, _ = ch.correction_boxes((23, 48, 41), OFFSET, "same", (2, 3, -2))
    want = matlab_boxes((48, 41, 23), OFFSET_M, "same", np.array([4.0, -1, 3]))
    assert [tuple(b) for b in boxes] == [tuple(box_from_matlab(w)) for w in want]


def test_same_mode_keeps_unshifted_data_in_place():
    # the centered origin: unshifted channels are not moved (tests/README.md)
    boxes, _ = ch.correction_boxes((23, 48, 41), np.zeros((2, 3)), "same")
    assert [tuple(b) for b in boxes] == [(0, 0, 0, 23, 48, 41)] * 2


def test_box_errors():
    with pytest.raises(ch.ChromaticError) as e:
        ch.correction_boxes((4, 5, 6), None)
    assert e.value.identifier == "MATLAB:sizeDimensionsMustMatch"
    with pytest.raises(ValueError):
        ch.correction_boxes((4, 5, 6), [[0, 0, 0.5]])


def test_psf_bounds_are_in_zyx(tmp_path):
    # (ChromaticShiftCorrectionTest/psfBoundsAreInYXZ) a bead 5 voxels off in y and 1 in x is
    # found within bounds of 6 in y and 2 in x as within 6 in all axes (MATLAB bounded y by the
    # x values and x by the y values before 45b2bb1)
    paths = []
    for c, off in enumerate(((0, 0, 0), (0, -5, -1))):
        p = tmp_path / f"psf_{c}.tif"
        synth.write_tiff(p, bead_psf(off))
        paths.append(str(p))
    wide = ch.chromatic_shift_estimation(paths, max_offset=(6, 6, 6))
    np.testing.assert_array_equal(np.abs(wide[1]), [0, 5, 1])
    np.testing.assert_array_equal(ch.chromatic_shift_estimation(paths, max_offset=(2, 6, 2)), wide)
    # x's own bound of 4 finds x = -5 of the demo's offsets
    paths = []
    for c in range(2):
        p = tmp_path / f"demo_{c}.tif"
        synth.write_tiff(p, bead_psf(OFFSET[c]))
        paths.append(str(p))
    np.testing.assert_array_equal(
        ch.chromatic_shift_estimation(paths, max_offset=(6, 4, 8))[1, 2], -5
    )


def test_estimation_crop(tmp_path):
    paths = []
    for c in range(2):
        p = tmp_path / f"psf_{c}.tif"
        synth.write_tiff(p, bead_psf(OFFSET[c]))
        paths.append(str(p))
    s = ch.chromatic_shift_estimation(paths, max_offset=(6, 8, 8), crop_length=(2, 3, 1))
    np.testing.assert_array_equal(s, OFFSET)


def test_estimation_errors(tmp_path):
    with pytest.raises(ch.ChromaticError, match="At least two"):
        ch.chromatic_shift_estimation([str(tmp_path / "a.tif")])
    # (a char is one path: its characters were taken as the files)
    with pytest.raises(ch.ChromaticError, match="At least two PSF files are needed!"):
        ch.chromatic_shift_estimation("ab.tif")


def test_offsets_have_a_row_per_channel(tmp_path):
    # (ChromaticShiftCorrectionTest/offsetsHaveARowPerChannel)
    path = write_dataset(tmp_path / "d")
    ids = []
    for off in (None, OFFSET[1:], np.concatenate([OFFSET, [[0, 9, 0]]])):
        with pytest.raises(ch.ChromaticError) as e:
            run(path, "valid", chromatic_offset=off)
        ids.append(e.value.identifier.split(":")[-1])
    assert ids == ["noOffsets", "offsetSize", "offsetSize"]


def test_mode_is_one_of_three_in_any_case(tmp_path):
    # (ChromaticShiftCorrectionTest/modeIsOneOfThreeInAnyCase)
    path = write_dataset(tmp_path / "d")
    with pytest.raises(ch.ChromaticError) as e:
        run(path, "none")
    assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
    assert not os.path.exists(f"{path}{RESULT}")
    run(path, "Valid", result_dir_name="upper")
    run(path, "valid")
    fn = f"{FSN.format('CamB')}.tif"
    np.testing.assert_array_equal(read_tiff(f"{path}upper/{fn}"), read_tiff(f"{path}{RESULT}/{fn}"))


def test_large_zarr_files_are_cropped_in_batches(tmp_path):
    # (ChromaticShiftCorrectionTest/largeZarrFilesAreCroppedInBatches) into zarr results equal to
    # the regular path's, with the MIPs of the batches; their pad value is 0
    path = write_dataset(tmp_path / "d", fmt="zarr")
    kw = {"zarr_file": True, "batch_size": (8, 16, 16), "block_size": (8, 8, 8)}
    for mode in ("valid", "full"):
        run(path, mode, save_zarr=True, result_dir_name=f"regular_{mode}", **kw)
        run(path, mode, large_file=True, result_dir_name=f"large_{mode}", **kw)
        for cam in ("CamA", "CamB"):
            fn = f"{FSN.format(cam)}.zarr"
            np.testing.assert_array_equal(
                read_zarr(f"{path}large_{mode}/{fn}"), read_zarr(f"{path}regular_{mode}/{fn}")
            )
            assert os.path.isfile(f"{path}large_{mode}/MIPs/{FSN.format(cam)}_MIP_z.tif")
    with pytest.raises(ch.ChromaticError) as e:
        run(path, "full", large_file=True, pad_value=7, result_dir_name="large_pad", **kw)
    assert e.value.identifier == "XR_chromatic_shift_correction_data_wrapper:padValue"


def test_zarr_data_with_offsets_from_psfs(tmp_path):
    # (ChromaticShiftCorrectionTest/zarrDataWithOffsetsFromPSFs) the PSFs are read as tiff, the
    # data as zarr
    path = write_dataset(tmp_path / "d", fmt="zarr", psfs=True)
    run(
        path,
        "valid",
        chromatic_offset=None,
        psf_fullpaths=[path + n for n in PSF_NAMES[:2]],
        max_offset=(6, 8, 8),
        zarr_file=True,
    )
    a, b = outputs(path)
    np.testing.assert_array_equal(b, a)


def test_matlab_layer_converts_offsets_and_origin(tmp_path, monkeypatch):
    from petakit5d import matlab

    seen = {}
    real = ch.correction_boxes

    def spy(n, o, m, origin):
        seen.update(o=o, origin=origin)
        return real(n, o, m, origin)

    monkeypatch.setattr(ch, "correction_boxes", spy)
    path = write_dataset(tmp_path / "d")
    matlab.XR_chromatic_shift_correction_data_wrapper(
        [path], "chromaticOffset", OFFSET_M, "mode", "same", "newOrigin", [4, -1, 3], *LOCAL
    )
    np.testing.assert_array_equal(np.asarray(seen["o"]), OFFSET)
    assert tuple(seen["origin"]) == (2, 3, -2)


def test_matlab_layer_returns_yxz_shifts(tmp_path):
    from petakit5d import matlab

    paths = []
    for c in range(2):
        p = tmp_path / f"psf_{c}.tif"
        synth.write_tiff(p, bead_psf(OFFSET[c]))
        paths.append(str(p))
    r = matlab.XR_chromatic_shift_estimation(paths, "maxOffset", [8, 8, 6])
    np.testing.assert_array_equal(np.asarray(r), OFFSET_M)


def test_results_do_not_depend_on_threads(tmp_path):
    paths = []
    for c in range(2):
        p = tmp_path / f"psf_{c}.tif"
        synth.write_tiff(p, bead_psf(OFFSET[c]))
        paths.append(str(p))
    r1 = tmp_path / "r1.mat"
    r4 = tmp_path / "r4.mat"
    ch.chromatic_shift_estimation(paths, result_fullpath=str(r1), num_threads=1)
    ch.chromatic_shift_estimation(paths, result_fullpath=str(r4), num_threads=4)
    a, b = read_mat(r1), read_mat(r4)
    for k in ("relative_shift_yxz_mat", "max_xcorr_mat"):
        assert np.asarray(a[k]).tobytes() == np.asarray(b[k]).tobytes()
    assert not any(n.startswith("r1_") for n in os.listdir(tmp_path))
