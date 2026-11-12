"""Fixtures of the parameter registry tests."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
DATA = HERE / "data"


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def sample_crop():
    """The worked example: XR_crop_dataset declared in a registry of its own."""
    return _load("sample_crop")


@pytest.fixture(scope="session")
def extract():
    """The stored inputParser extract of the oracle tree."""
    from petakit5d.params.check import read_extract

    return read_extract(DATA / "inputparser_extract.json")


@pytest.fixture(scope="session")
def extract_json():
    return json.loads((DATA / "inputparser_extract.json").read_text())
