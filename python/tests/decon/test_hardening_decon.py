"""Tests the GPU round of the mutation campaign asked for (docs/hardening.md): decon."""

import importlib

import numpy as np
import pytest

from petakit5d.testing.synth import CounterRNG, decon_psf


def _rl_gpu():
    # imported when a test runs, as in tests/gpu/test_decon_gpu.py
    return importlib.import_module("petakit5d.decon.rl_gpu")


@pytest.fixture
def gpu():
    cp = pytest.importorskip("cupy")
    from petakit5d.gpu import device, fft

    n, reason, _ = device.usable_gpu_count()
    if n == 0:
        pytest.skip(f"no usable GPU: {reason}")
    device.select_device()
    _rl_gpu().clear_otf_cache()
    yield cp
    _rl_gpu().clear_otf_cache()
    fft.clear_plan_cache()
    cp.get_default_memory_pool().free_all_blocks()


def _data(dtype=np.uint16):
    im = CounterRNG(7, "hardening").integers(100, 1100, (20, 24, 22), np.uint16)
    return im.astype(dtype), decon_psf((15, 15, 15), (2.0, 1.5))


@pytest.mark.gpu
def test_without_kept_plans_the_otfs_plan_is_closed_at_once(gpu, monkeypatch):
    # keep_plans=False: the FFT plan the OTFs are made with (and its GPU work area) is closed
    # when they are made, not left to the garbage collector
    from petakit5d.gpu import fft

    made = []

    class Recorded(fft.ComplexFFTPlan):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            made.append(self)

    monkeypatch.setattr(fft, "ComplexFFTPlan", Recorded)
    im, psf = _data()
    _rl_gpu().decon_lucy(
        im, psf, 2, background=100.0, use_gpu=True, use_mex=False, keep_plans=False
    )
    assert made and all(p.closed for p in made)


@pytest.mark.gpu
def test_the_frames_gpu_check_reads_the_free_memory(gpu):
    # the GPU's free memory as a number (MATLAB's gpuDevice().AvailableMemory): None, which sends
    # a frame's decon to the CPU, only without a usable GPU
    from petakit5d.decon import frame
    from petakit5d.gpu import device

    free = frame.gpu_free_bytes()
    assert type(free) is int and 0 < free <= device.select_device().total_bytes
    assert frame.gpu_has_room((20, 24, 22), "omw")


@pytest.mark.gpu
def test_the_debug_rows_of_double_data_count_every_iteration(gpu, tmp_path):
    # double data take MATLAB's first two iterations in double; the debug output books them as
    # iterations 1 and 2 (err_mat's first column, Iter_0001.tif and Iter_0002.tif)
    im, psf = _data(np.float64)
    res = _rl_gpu().decon_lucy(
        im,
        psf,
        3,
        background=100.0,
        use_gpu=True,
        use_mex=False,
        debug=True,
        debug_folder=str(tmp_path),
        save_step=1,
    )
    assert res.err_mat[:, 0].tolist() == [1.0, 2.0, 3.0]
    assert sorted(p.name for p in tmp_path.glob("Iter_*.tif")) == [
        "Iter_0001.tif",
        "Iter_0002.tif",
        "Iter_0003.tif",
    ]
