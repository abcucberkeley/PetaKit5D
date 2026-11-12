"""The oracle harness finds the repo's MATLAB side, and the python/ folder that cases' helper
folders are relative to, also when petakit5d is installed (CI)."""

from pathlib import Path

import pytest

from petakit5d.testing import oracle


def test_harness_dir_is_next_to_the_tests():
    want = Path(__file__).resolve().parents[1] / "oracle_matlab"
    assert oracle.harness_dir().resolve() == want
    assert (oracle.harness_dir() / "ORACLE_COMMIT").is_file()


def test_python_root_is_the_tests_python_folder(monkeypatch, tmp_path):
    # also when petakit5d is installed: its own folder is then in site-packages
    from petakit5d.testing import cases

    want = Path(__file__).resolve().parents[2]
    assert oracle.python_root() == want and cases._python_root() == want
    site = tmp_path / "site-packages" / "petakit5d" / "testing"
    monkeypatch.setattr(cases, "__file__", str(site / "cases.py"))
    monkeypatch.setattr(oracle, "__file__", str(site / "oracle.py"))
    assert oracle.python_root() == want and cases._python_root() == want


def test_a_missing_helper_folder_is_an_error():
    # not the hash of an empty folder, which gave every golden case another key in CI
    from petakit5d.testing.cases import Case

    c = Case("missing_helpers", matlab="disp", args=["x"], helpers=["tests/no_such_folder"])
    with pytest.raises(FileNotFoundError, match="helper folder tests/no_such_folder"):
        c.key()
