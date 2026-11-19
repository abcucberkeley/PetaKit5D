"""Tests the GPU round of the mutation campaign asked for (docs/hardening.md): gpu."""

import ctypes
import sys

import pytest

from petakit5d.gpu import device, fft, kernels


@pytest.fixture
def gpu():
    cp = pytest.importorskip("cupy")
    n, reason, _ = device.usable_gpu_count()
    if n == 0:
        pytest.skip(f"no usable GPU: {reason}")
    device.select_device()
    fft.clear_plan_cache()
    yield cp
    fft.clear_plan_cache()
    cp.get_default_memory_pool().free_all_blocks()


def test_launch_dimensions_take_one_to_three_numbers():
    # a grid or block of no dimension or of more than three is refused before a launch
    assert kernels._dim3(5) == (5, 1, 1)
    assert kernels._dim3((2, 3)) == (2, 3, 1)
    assert kernels._dim3([4, 3, 2]) == (4, 3, 2)
    for bad in ((), (1, 2, 3, 4)):
        with pytest.raises(ValueError, match="1 to 3 dimensions"):
            kernels._dim3(bad)


@pytest.mark.gpu
def test_a_plan_used_again_is_the_last_the_cache_drops(gpu, monkeypatch):
    # the cache keeps the plans used last: a hit makes a plan the newest
    monkeypatch.setattr(fft, "_plan_cache_size", 2)
    a = fft.get_plan((4, 6, 10))
    b = fft.get_plan((4, 6, 12))
    assert fft.get_plan((4, 6, 10)) is a
    fft.get_plan((4, 6, 14))  # drops the plan used least recently: b
    assert fft.get_plan((4, 6, 10)) is a
    assert fft.get_plan((4, 6, 12)) is not b


@pytest.mark.gpu
def test_cufft_found_nowhere_is_loaded_by_the_platforms_library_name(gpu, monkeypatch):
    # without PETAKIT5D_CUFFT_LIBRARY and without the wheel's cuFFT (cuda-pathfinder), the
    # library is loaded by the name of CUDA 12's cuFFT on the platform
    pathfinder = pytest.importorskip("cuda.pathfinder")
    tried = []

    def cdll(path, *args, **kwargs):
        tried.append(path)
        raise OSError("not loaded by this test")

    def not_found(name):
        raise RuntimeError(f"{name} not found")

    monkeypatch.setattr(fft, "_lib", None)
    monkeypatch.delenv("PETAKIT5D_CUFFT_LIBRARY", raising=False)
    monkeypatch.setattr(pathfinder, "load_nvidia_dynamic_lib", not_found)
    monkeypatch.setattr(ctypes, "CDLL", cdll)
    with pytest.raises(OSError, match="not loaded by this test"):
        fft._library()
    assert tried == ["cufft64_11.dll" if sys.platform == "win32" else "libcufft.so.11"]
