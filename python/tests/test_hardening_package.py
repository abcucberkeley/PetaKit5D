"""Tests the mutation campaign asked for (E06, docs/hardening.md): the package's lazy imports."""

import importlib

import pytest

import petakit5d


def test_a_subpackage_whose_dependency_is_missing_raises_that(monkeypatch):
    real = importlib.import_module

    def fake(name, *a, **k):
        if name == "petakit5d.hardening_probe":
            raise ModuleNotFoundError("No module named 'missing_dep'", name="missing_dep")
        return real(name, *a, **k)

    monkeypatch.setattr(importlib, "import_module", fake)
    with pytest.raises(ModuleNotFoundError, match="missing_dep"):
        petakit5d.__getattr__("hardening_probe")
    with pytest.raises(AttributeError, match="no attribute 'no_such_subpackage'"):
        petakit5d.__getattr__("no_such_subpackage")
