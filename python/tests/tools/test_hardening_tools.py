"""Tests the mutation campaign asked for (E06, docs/hardening.md): tools."""

import numpy as np
import pytest

from petakit5d.io.tiff import read_tiff, write_tiff
from petakit5d.jobs import framework
from petakit5d.matlabcompat import elementary as E
from petakit5d.testing import synth
from petakit5d.tools import fft_spectrum as FS
from petakit5d.tools import fsc, psf_analysis, puncta_removal, rsync, unmix


def test_file_checks_on_the_cluster_raise_when_a_check_job_failed(tmp_path, monkeypatch):
    monkeypatch.setattr(framework, "run_tasks", lambda i, o, t, **kw: [False] * len(o))
    with pytest.raises(RuntimeError, match="not every check job finished"):
        rsync.batch_file_exist_cluster([str(tmp_path / "a")], str(tmp_path), parse_cluster=False)


def test_numpy_values_keep_their_form_through_the_text_rounding():
    r = np.array([1.5, 2.0])
    assert unmix._like(np.array([1.4, 2.1]), r) is r
    s = np.array(3.0)
    assert unmix._like(np.float32(2), s) is s
    assert unmix._like((1, 2), [4.0, 5.0]) == (4.0, 5.0)


@pytest.mark.parametrize("y", [0.37, 2.5, -1.3, 3.0, 1 / 3])
def test_fsc_power_is_matlabs_power_on_its_fast_path(y, monkeypatch):
    x = np.linspace(0.25, 50.0, 103)
    want = E.power(x, y)
    slow = []
    real = E.power

    def counted(a, b):
        slow.append(np.size(a))
        return real(a, b)

    monkeypatch.setattr(E, "power", counted)
    np.testing.assert_array_equal(fsc._pow_positive(x, y).view(np.uint64), want.view(np.uint64))
    assert sum(slow) == 0, "normal values went to the scalar power"


class _Stop(Exception):
    pass


def test_psf_detection_tasks_ask_for_100_times_the_frame_in_single(tmp_path, monkeypatch):
    im = np.zeros((21, 40, 36), dtype=np.uint16)
    im[10, 20, 18] = 1000
    write_tiff(im, str(tmp_path / "beads.tif"))
    mems = []

    def spy(inputs, outputs, tasks, **kw):
        mems.append(kw["mem_allocate"])
        raise _Stop

    monkeypatch.setattr(framework, "run_tasks", spy)
    with pytest.raises(_Stop):
        psf_analysis.psf_detection_and_analysis_dataset(
            str(tmp_path) + "/", crop_size=(9, 9, 9), dist_thresh=(9, 9, 9),
            channel_patterns=("beads",), channels=(488,),
        )  # fmt: skip
    assert mems == [21 * 40 * 36 * 4 / 1024**3 * 100]


def _two_frames(folder):
    rng = np.random.default_rng(5)
    for n in ("a_CamA", "b_CamA"):
        write_tiff(
            rng.integers(0, 2000, size=(8, 20, 18), dtype=np.uint16), str(folder / f"{n}.tif")
        )


def test_fft_spectra_of_tasks_not_all_done_are_tried_again(tmp_path, monkeypatch):
    _two_frames(tmp_path)
    calls = []

    def spy(inputs, outputs, tasks, **kw):
        calls.append(kw["mem_allocate"])
        return [len(calls) > 1, True]

    monkeypatch.setattr(framework, "run_tasks", spy)
    done = FS.fft_spectrum_dataset(str(tmp_path) + "/", out_size=(15, 15, 15), parse_cluster=False)
    assert done == [True, True]
    assert len(calls) == 2 and calls[1] == 2 * calls[0]


def test_fft_spectrum_takes_a_numpy_integer_background(tmp_path):
    _two_frames(tmp_path)
    fn = str(tmp_path / "a_CamA.tif")
    outs = []
    for k, bg in enumerate((np.int64(100), 100.0)):
        res = str(tmp_path / f"fft_{k}.tif")
        FS.fft_spectrum_file(fn, res, out_size=(15, 15, 15), background=bg, save_3d_stack=True)
        outs.append(read_tiff(res))
    np.testing.assert_array_equal(outs[0], outs[1])


def test_3d_detection_keeps_puncta_above_the_local_noise_only():
    rng = np.random.default_rng(12)
    vol = (100 + 10 * rng.standard_normal((20, 24, 22))).astype(np.float32)
    vol += synth.beads(vol.shape, np.array([[10.0, 12.0, 11.0]]), 1.5, 2000.0, 0.0)
    zs, ys, xs, kept = puncta_removal._detect_3d(vol, (3, 3, 3), 3.0)
    assert kept.size > 5, "noise gives local maxima"
    strong = [i for i in range(zs.size) if (zs[i], ys[i], xs[i]) == (10, 12, 11)]
    assert strong and kept[strong[0]]
    assert kept.sum() < kept.size / 2, "noise maxima were kept"


def test_mip_mask_warnings_and_errors_keep_their_identifiers():
    from petakit5d.tools.mip_masks import MipMaskError, MipMaskWarning

    assert MipMaskWarning("m", "a:b").identifier == "a:b" and MipMaskWarning("m").identifier == ""
    assert MipMaskError("m", "a:b").identifier == "a:b"
