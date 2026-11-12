"""The core oracle cases against MATLAB's references (oracle_f02.json, by hash).

Stand-in until F04's harness lands: the manifest was made by running ``oracle_cases`` in MATLAB
R2024b at the oracle commit named in it. Cases that need a file tree run only on a
case-sensitive POSIX file system (the tree has names that differ only in case, and links).
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import _f02_oracle as _oracle  # noqa: E402
import f02_oracle_cases as oc  # noqa: E402

MANIFEST = Path(__file__).with_name("oracle_f02.json")
_REFS = json.loads(MANIFEST.read_text(encoding="utf-8"))["cases"] if MANIFEST.exists() else {}

_TREE_CASES = ("dir.", "dir_recursive.", "simplifyPath")
_POSIX_ONLY = ("fileparts.",)


def _case_sensitive_posix(tmp) -> bool:
    if sys.platform.startswith("win"):
        return False
    (tmp / "A").write_text("")
    ok = not (tmp / "a").exists()
    try:
        os.symlink("A", tmp / "l")
    except OSError:
        ok = False
    return ok


@pytest.mark.oracle
@pytest.mark.parametrize("c", oc.CASES, ids=[c.id for c in oc.CASES])
def test_matches_matlab(c, tmp_path):
    ref = _REFS.get(c.id)
    if ref is None:
        pytest.fail(f"no MATLAB reference for {c.id} in {MANIFEST.name}")
    if (
        c.id.startswith(_TREE_CASES)
        and c.id != "dir_recursive.missing"
        and not _case_sensitive_posix(tmp_path)
    ):
        pytest.skip("needs a case-sensitive POSIX file system with symbolic links")
    if c.id.startswith(_POSIX_ONLY) and sys.platform.startswith("win"):
        pytest.skip("MATLAB's fileparts follows other rules on Windows")
    inputs = {k: mk() for k, mk in c.inputs.items()}
    if ref["kind"] == "error":
        with pytest.raises(Exception):  # noqa: B017 (MATLAB errors; the Python error type is the port's)
            c.python(**inputs)
        return
    got = _oracle.entry_of_python(c.kind, c.python(**inputs))
    assert _oracle.same(ref, got), _oracle.describe(ref, got)


def test_every_case_has_a_reference():
    assert MANIFEST.exists(), "oracle_f02.json is missing"
    missing = [c.id for c in oc.CASES if c.id not in _REFS]
    assert not missing, missing
