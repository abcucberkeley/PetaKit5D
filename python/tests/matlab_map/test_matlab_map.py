"""The map of the MATLAB test suite (tests/matlab_map/*.json) and its checker.

test_map_is_consistent is the CI check: every MATLAB test is in the map once, the counts
match suite.json, every fate points at something that exists, and every test of a finished
class has a fate.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path

import pytest

PYTHON_ROOT = Path(__file__).resolve().parents[2]
MAP_DIR = PYTHON_ROOT / "tests" / "matlab_map"


def _tool():
    spec = importlib.util.spec_from_file_location(
        "check_test_map", PYTHON_ROOT / "tools" / "check_test_map.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ctm = _tool()


def test_map_is_consistent():
    errs, counts = ctm.check(MAP_DIR, PYTHON_ROOT)
    assert not errs, "\n".join(errs[:50])
    suite = json.loads((MAP_DIR / "suite.json").read_text())
    assert counts == {"n_all": suite["n_all"], "n_default": suite["n_default"]}


def test_map_is_of_the_repos_oracle_commit():
    f = PYTHON_ROOT / "tests" / "oracle_matlab" / "ORACLE_COMMIT"
    if not f.is_file():
        pytest.skip("no tests/oracle_matlab/ORACLE_COMMIT (oracle harness not in this tree)")
    commit = f.read_text().split()[0]
    mapped = json.loads((MAP_DIR / "suite.json").read_text())["oracle_commit"]
    assert commit.startswith(mapped) or mapped.startswith(commit), (
        f"the map is of oracle {mapped}, the repo's oracle is {commit}: "
        "run check_test_map.py update"
    )


def test_report_and_docs_render():
    suite, folders = ctm.load_map(MAP_DIR)
    text = "\n".join(ctm.report_lines(suite, folders))
    assert f"{suite['n_default']} in a default run" in text
    assert "| total |" in text
    assert "## Classes" in ctm.docs_text(suite, folders)


# ------------------------------------------------------------------------------------------
# the checker on small maps
# ------------------------------------------------------------------------------------------


def _small_map(tmp_path: Path) -> Path:
    d = tmp_path / "map"
    d.mkdir(parents=True)
    suite = {
        "oracle_commit": "abc1234",
        "matlab": "R2024b",
        "n_all": 4,
        "n_default": 3,
        "folders": {"crop": 4},
    }
    (d / "suite.json").write_text(json.dumps(suite))
    folder = {
        "folder": "crop",
        "classes": {
            "CropMexTest": {
                "file": "tests/crop/CropMexTest.m",
                "pins": "crop",
                "category": "data",
                "owner": "K01",
                "phase": 1,
                "finished": False,
                "fate": [],
                "methods": {
                    "crop3dMatchesIndexing": {
                        "pins": "p",
                        "tags": ["Unit"],
                        "cases": ["dtype=uint8", "dtype=single"],
                        "fate": [],
                    },
                    "rejectsBbox": {"pins": "p", "tags": ["Unit"], "cases": None, "fate": []},
                    "negativeValues": {
                        "pins": "p",
                        "tags": ["KnownIssue", "Unit"],
                        "cases": None,
                        "fate": [],
                        "known_issue": True,
                    },
                },
            }
        },
    }
    ctm.write_folder(d, "crop", folder)
    return d


def _edit(d: Path, fn) -> None:
    data = json.loads((d / "crop.json").read_text())
    fn(data["classes"]["CropMexTest"])
    ctm.write_folder(d, "crop", data)


def test_checker_accepts_a_pending_unfinished_class(tmp_path):
    d = _small_map(tmp_path)
    assert ctm.check(d, PYTHON_ROOT)[0] == []
    errs, _ = ctm.check(d, PYTHON_ROOT, strict=True)
    assert len(errs) == 4 and all("no fate" in e for e in errs)


def test_checker_wants_fates_in_finished_classes(tmp_path):
    d = _small_map(tmp_path)
    _edit(d, lambda c: c.update(finished=True))
    errs, _ = ctm.check(d, PYTHON_ROOT)
    assert len(errs) == 4
    me = "tests/matlab_map/test_matlab_map.py"

    def fill(c):
        c["methods"]["crop3dMatchesIndexing"]["fate"] = [f"pytest:{me}::test_map_is_consistent"]
        c["methods"]["crop3dMatchesIndexing"]["case_fates"] = {
            "dtype=single": ["na:single is checked elsewhere"]
        }
        c["methods"]["rejectsBbox"]["fate"] = [
            "oracle:tests/matlab_map/test_golden.py::golden.decon_omw.*"
        ]
        c["fate"] = ["na:the rest of the class is not applicable"]

    _edit(d, fill)
    assert ctm.check(d, PYTHON_ROOT)[0] == []


@pytest.mark.parametrize(
    "fate, message",
    [
        ("pytest:tests/nothing_here.py::test_x", "no file"),
        ("pytest:tests/matlab_map/test_matlab_map.py::test_does_not_exist", "no test"),
        ("pytest:tests/matlab_map/test_matlab_map.py", "needs <path>::<name>"),
        ("oracle:tests/matlab_map/test_golden.py::golden.no_such_scenario.*", "no case"),
        ("na:", "give the reason"),
        ("maybe:later", "bad fate kind"),
    ],
)
def test_checker_rejects_bad_fates(tmp_path, fate, message):
    d = _small_map(tmp_path)
    _edit(d, lambda c: c["methods"]["rejectsBbox"].update(fate=[fate]))
    errs, _ = ctm.check(d, PYTHON_ROOT)
    assert len(errs) == 1 and message in errs[0], errs


def test_pytest_fates_may_name_parameters_and_classes(tmp_path):
    d = _small_map(tmp_path)
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_x.py").write_text(
        "class TestCrop:\n    def test_a(self):\n        pass\n"
    )
    _edit(
        d,
        lambda c: c["methods"]["rejectsBbox"].update(
            fate=["pytest:tests/test_x.py::TestCrop::test_a[uint8]"]
        ),
    )
    assert ctm.check(d, tmp_path)[0] == []


def test_checker_counts_tests(tmp_path):
    d = _small_map(tmp_path)
    _edit(d, lambda c: c["methods"].pop("rejectsBbox"))
    errs, _ = ctm.check(d, PYTHON_ROOT)
    assert any("the map has 3 tests, suite.json 4" in e for e in errs)
    d = _small_map(tmp_path / "b")
    _edit(d, lambda c: c["methods"]["negativeValues"].pop("known_issue"))
    errs, _ = ctm.check(d, PYTHON_ROOT)
    assert any("without KnownIssue" in e for e in errs)


def test_checker_wants_category_and_owner(tmp_path):
    d = _small_map(tmp_path)
    _edit(d, lambda c: c.update(category="unclassified", owner="unassigned"))
    errs, _ = ctm.check(d, PYTHON_ROOT)
    assert any("unclassified" in e for e in errs) and any("no owner" in e for e in errs)


def test_checker_finds_duplicates_and_unknown_case_fates(tmp_path):
    d = _small_map(tmp_path)

    def dup(c):
        m = c["methods"]["crop3dMatchesIndexing"]
        m["cases"] = ["dtype=uint8", "dtype=uint8"]
        m["case_fates"] = {"dtype=double": ["na:no such case here"]}

    _edit(d, dup)
    errs, _ = ctm.check(d, PYTHON_ROOT)
    assert any("listed twice" in e for e in errs)
    assert any("unknown case (dtype=double)" in e for e in errs)


def _dump(tests: list[tuple[str, str, list[str]]]) -> dict:
    out = []
    for cls_meth, tags, params in tests:
        cls, meth = cls_meth.split("/")
        for p in params or [None]:
            out.append(
                {
                    "name": cls_meth + (f"({p})" if p else ""),
                    "procedure": meth,
                    "class": cls,
                    "folder": "crop",
                    "tags": tags,
                    "params": [],
                }
            )
    n_ki = sum(len(p or [None]) for _, t, p in tests if "KnownIssue" in t)
    return {
        "matlab": "R2024b",
        "computer": "GLNXA64",
        "n_all": len(out),
        "n_default": len(out) - n_ki,
        "tests": out,
    }


def test_update_keeps_fates_and_adds_new_tests(tmp_path):
    d = _small_map(tmp_path)
    _edit(
        d,
        lambda c: (
            c["methods"]["rejectsBbox"].update(fate=["na:kept fate"]),
            c.update(finished=True),
        ),
    )
    tests_dir = tmp_path / "tests"
    (tests_dir / "crop").mkdir(parents=True)
    source = [
        "classdef CropMexTest < pktest.TestCase",
        "    % crop mex functions",
        "    methods (Test)",
        "        function rejectsBbox(testCase)",
        "            % invalid boxes are rejected",
        "        end",
        "        function newCheck(testCase)",
        "            % a new check",
        "        end",
        "    end",
        "end",
    ]
    (tests_dir / "crop" / "CropMexTest.m").write_text("\n".join(source) + "\n")
    dump = _dump(
        [
            (
                "CropMexTest/crop3dMatchesIndexing",
                ["Unit"],
                ["dtype=uint8", "dtype=single", "dtype=double"],
            ),
            ("CropMexTest/rejectsBbox", ["Unit"], []),
            ("CropMexTest/newCheck", ["Unit"], []),
            ("NewClassTest/works", ["Unit"], []),
        ]
    )
    (tmp_path / "dump.json").write_text(json.dumps(dump))
    notes = ctm.update(tmp_path / "dump.json", tests_dir, "def5678", d)
    assert "new class crop/NewClassTest (unclassified, no owner)" in notes
    assert "new method CropMexTest/newCheck" in notes
    assert "removed crop/CropMexTest/negativeValues" in notes
    assert any("crop3dMatchesIndexing: cases changed" in n for n in notes)
    suite, folders = ctm.load_map(d)
    assert suite["oracle_commit"] == "def5678" and suite["n_all"] == 6
    cls = folders["crop"]["classes"]["CropMexTest"]
    assert cls["owner"] == "K01" and cls["finished"] is False  # a new method reopens the class
    assert cls["methods"]["rejectsBbox"]["fate"] == ["na:kept fate"]
    assert cls["methods"]["rejectsBbox"]["pins"] == "invalid boxes are rejected"
    assert cls["methods"]["newCheck"]["fate"] == []
    errs, _ = ctm.check(d, PYTHON_ROOT)
    assert any("NewClassTest: unclassified" in e for e in errs)


def test_update_from_the_repos_map_is_a_fixed_point(tmp_path):
    # the map rebuilt from its own inventory (without the MATLAB files) is the same map
    d = tmp_path / "map"
    shutil.copytree(MAP_DIR, d, ignore=shutil.ignore_patterns("*.py", "matlab", "__pycache__"))
    suite, folders = ctm.load_map(d)
    tests = []
    for folder, cls, meth, p, _cr, mr in ctm.iter_cases(folders):
        tests.append(
            {
                "name": ctm.case_name(cls, meth, p),
                "procedure": meth,
                "class": cls,
                "folder": folder,
                "tags": mr["tags"],
                "params": [],
            }
        )
    dump = {
        "matlab": suite["matlab"],
        "computer": suite.get("computer", ""),
        "n_all": suite["n_all"],
        "n_default": suite["n_default"],
        "tests": tests,
    }
    (tmp_path / "dump.json").write_text(json.dumps(dump))
    notes = ctm.update(tmp_path / "dump.json", tmp_path / "no_tests_dir", suite["oracle_commit"], d)
    assert notes == []
    _, again = ctm.load_map(d)
    assert again == folders
