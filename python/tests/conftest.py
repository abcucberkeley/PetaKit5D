"""Shared test set-up.

Markers (declared in pyproject.toml):
  gpu       skipped without a usable GPU (CuPy, a device, CUDA_VISIBLE_DEVICES not empty)
  cluster   skipped unless PETAKIT5D_TEST_SLURM_CONFIG names a cluster config file
  demodata  skipped unless PETAKIT5D_DEMO_DATA names the demo data folder
  slow      runs by default; deselect with -m "not slow"
  oracle    runs by default; the cases compare with the reference store
            (PETAKIT5D_ORACLE_DIR, PETAKIT5D_GOLDEN_DIR) or with the repo manifests
Site values come only from the environment: PETAKIT5D_ORACLE_DIR, PETAKIT5D_GOLDEN_DIR,
PETAKIT5D_DEMO_DATA, PETAKIT5D_TEST_SLURM_CONFIG, PETAKIT5D_TEST_SLURM_SCRATCH,
PETAKIT5D_TEST_FFTW.
"""

import functools
import os
import sys

import pytest

SITE_VARIABLES = (
    "PETAKIT5D_ORACLE_DIR",
    "PETAKIT5D_GOLDEN_DIR",
    "PETAKIT5D_DEMO_DATA",
    "PETAKIT5D_TEST_SLURM_CONFIG",
    "PETAKIT5D_TEST_SLURM_SCRATCH",
    "PETAKIT5D_TEST_FFTW",
    "PETAKIT5D_EXT_DIR",
)


@functools.cache
def gpu_skip_reason() -> str | None:
    """Why GPU tests can't run here, or None. CuPy is imported only if devices are visible."""
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


def _path_skip_reason(variable: str) -> str | None:
    value = os.environ.get(variable)
    if not value:
        return f"{variable} is not set"
    if not os.path.exists(value):
        return f"{variable} does not exist"
    return None


def pytest_collection_modifyitems(config, items):
    reasons = {}
    for item in items:
        for marker, check in (
            ("cluster", lambda: _path_skip_reason("PETAKIT5D_TEST_SLURM_CONFIG")),
            ("demodata", lambda: _path_skip_reason("PETAKIT5D_DEMO_DATA")),
            ("gpu", gpu_skip_reason),
        ):
            if item.get_closest_marker(marker) is None:
                continue
            if marker not in reasons:
                reasons[marker] = check()
            if reasons[marker]:
                item.add_marker(pytest.mark.skip(reason=f"{marker}: {reasons[marker]}"))


def pytest_report_header(config):
    import petakit5d

    lines = [
        f"petakit5d {petakit5d.__version__} from {os.path.dirname(petakit5d.__file__)}",
        f"python {sys.version.split()[0]}, OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}, "
        f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES')}",
        "site variables set: "
        + (", ".join(v for v in SITE_VARIABLES if os.environ.get(v)) or "none"),
    ]
    return lines
