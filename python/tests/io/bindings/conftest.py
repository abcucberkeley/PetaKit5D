"""Checks of the published cpp-tiff and cpp-zarr Python packages (task IO1)."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import iobind_cases as ic  # noqa: E402


def pytest_configure(config):
    # (the shared markers are in pyproject.toml)
    config.addinivalue_line(
        "markers", "gap: a known gap in cpp-tiff or cpp-zarr, listed for the maintainers"
    )


def pytest_report_header(config):
    v = ic.versions()
    return f"io bindings: cpp-tiff {v['cpp-tiff']}, cpp-zarr {v['cpp-zarr']}, numpy {v['numpy']}"


@pytest.fixture(scope="session")
def binding_versions():
    return ic.versions()
