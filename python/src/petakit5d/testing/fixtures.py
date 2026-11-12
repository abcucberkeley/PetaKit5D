"""Shared pytest fixtures (the helpers of the MATLAB suite's pktest.TestCase).

Use them as a plugin (``pytest_plugins = ["petakit5d.testing.fixtures"]`` in the root
conftest) or import the ones a module needs (``from petakit5d.testing.fixtures import
pk_temp_folder``). Nothing here imports NumPy or the package's areas at import time.

    pk_temp_folder   a new folder as a str with '/' separators and a trailing '/', as
                     createTempFolder (the demos' and wrappers' path style on every system)
    rng              CounterRNG(0): the same numbers in every test, as rng(0) in pktest
    thread_count     parametrized over THREAD_COUNTS; with limit_threads for OpenMP children
    limit_threads    sets OMP_NUM_THREADS (and the BLAS variables) for processes a test starts
    restore_cwd      restores the current folder after the test
    require_gpu      skips the test without a usable GPU (same rule as the `gpu` marker)
    demo_data        the demo data folder (PETAKIT5D_DEMO_DATA), or skips
    write_volume     writes a (z, y, x) array as .tif or .zarr (native or MATLAB's layout) into
                     the test's folder, returns the path
    zarr_layout      parametrized over 'native' and 'matlab' (the two zarr layouts)
    paired_input     writes one volume for both sides: tiff shared, zarr in each side's layout
"""

from __future__ import annotations

import functools
import os
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pytest

# thread counts for tests that results don't depend on the thread count: one, an odd count
# that divides no block size, and more than most CI runners have
THREAD_COUNTS = (1, 3, 8)
_THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "PETAKIT5D_NUM_THREADS",
)


def available_cpus() -> int:
    """the CPUs this process may run on (its affinity where the system has one)"""
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def temp_folder_str(path: str | os.PathLike) -> str:
    """a folder path with '/' separators and a trailing '/'"""
    return os.fspath(path).replace("\\", "/").rstrip("/") + "/"


@functools.cache
def gpu_skip_reason() -> str | None:
    """why GPU tests can't run here, or None (CuPy is imported only if GPUs are visible)"""
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None and visible.strip() in ("", "-1", "NoDevFiles"):
        return "CUDA_VISIBLE_DEVICES hides all GPUs"
    try:
        import cupy
    except Exception as err:  # ImportError, or a CUDA library that fails to load
        return f"CuPy not importable ({type(err).__name__}: {err})"
    try:
        count = cupy.cuda.runtime.getDeviceCount()
    except Exception as err:
        return f"no CUDA device ({type(err).__name__}: {err})"
    return None if count > 0 else "no CUDA device"


@pytest.fixture
def pk_temp_folder(tmp_path: Path) -> str:
    """a new, empty folder of the test: '/' separated, with a trailing '/'"""
    d = tmp_path / "pk"
    d.mkdir()
    return temp_folder_str(d)


@pytest.fixture
def rng():
    """CounterRNG(0), the same stream in every test"""
    from petakit5d.testing.synth import CounterRNG

    return CounterRNG(0)


@pytest.fixture
def limit_threads(monkeypatch: pytest.MonkeyPatch) -> Callable[[int], int]:
    """limit_threads(n) sets the thread variables for processes the test starts; returns n"""

    def set_threads(n: int) -> int:
        n = int(n)
        if n < 1:
            raise ValueError("thread count must be >= 1")
        for v in _THREAD_VARIABLES:
            monkeypatch.setenv(v, str(n))
        return n

    return set_threads


@pytest.fixture(params=THREAD_COUNTS, ids=lambda n: f"threads{n}")
def thread_count(request: pytest.FixtureRequest, limit_threads: Callable[[int], int]) -> int:
    """each of THREAD_COUNTS in turn (also set for child processes)"""
    return limit_threads(request.param)


@pytest.fixture
def restore_cwd() -> Iterator[Path]:
    """the current folder at the start; restored after the test"""
    cwd = Path.cwd()
    try:
        yield cwd
    finally:
        os.chdir(cwd)


@pytest.fixture
def require_gpu() -> None:
    reason = gpu_skip_reason()
    if reason:
        pytest.skip(f"gpu: {reason}")


@pytest.fixture
def demo_data() -> Path:
    """the folder with the demo datasets (PETAKIT5D_DEMO_DATA)"""
    value = os.environ.get("PETAKIT5D_DEMO_DATA")
    if not value or not os.path.isdir(value):
        pytest.skip("demodata: PETAKIT5D_DEMO_DATA is not set or not a folder")
    return Path(value)


@pytest.fixture
def write_volume(pk_temp_folder: str) -> Callable[..., str]:
    """write_volume(name, vol, fmt='tif', block_size=(32, 32, 32), layout='native') writes a
    (z, y, x) volume into the test's folder (name without extension; a zarr native, or in
    MATLAB's layout with layout='matlab'; block_size in (z, y, x)) and returns the path as a
    '/' separated str"""
    from petakit5d.testing import synth

    def write(
        name: str,
        vol: Any,
        fmt: str = "tif",
        block_size: Any = (32, 32, 32),
        layout: str = "native",
    ) -> str:
        path = f"{pk_temp_folder}{name}.{fmt}"
        if fmt == "tif":
            synth.write_tiff(path, vol)
        elif fmt == "zarr":
            synth.write_zarr(path, vol, block_size, layout)
        else:
            raise ValueError("fmt: 'tif' or 'zarr'")
        return path

    return write


@pytest.fixture(params=("native", "matlab"), ids=lambda s: f"zarr_{s}")
def zarr_layout(request: pytest.FixtureRequest) -> str:
    """each zarr layout in turn, 'native' (the port's) and 'matlab' (MATLAB's writezarr): for
    the tests that read both"""
    return request.param


@pytest.fixture
def paired_input(pk_temp_folder: str) -> Callable[..., Any]:
    """paired_input(name, vol, fmt='zarr', block_size=(32, 32, 32)) writes a (z, y, x) volume
    for both sides (synth.write_paired_inputs: a tiff shared; a zarr in MATLAB's layout for
    MATLAB and native for the port) and returns its PairedInput(matlab=..., python=...)"""
    from petakit5d.testing import synth

    def write(name: str, vol: Any, fmt: str = "zarr", block_size: Any = (32, 32, 32)) -> Any:
        return synth.write_paired_inputs(pk_temp_folder, name, vol, fmt, block_size)

    return write
