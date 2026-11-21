"""Third-party MATLAB code: the repository keeps its interface, not its validator texts.

The stored inputParser extract withholds the validators of every function whose file the
provenance rules call third-party (its header names an outside author) and keeps their names,
kinds, order and defaults; check_registry compares those only; the mapping table and the spec
skeleton show no validator of such code (decided on 2026-11-20).
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

from petakit5d.params import Param, Registry
from petakit5d.params.check import check_spec, extract_text, read_extract, write_extract
from petakit5d.params.mapping import _planned_section, _spec_section
from petakit5d.params.matlab_source import MatlabFunction, MatlabParam
from petakit5d.params.provenance import WITHHELD, ProvenanceError, interface_only, third_party

PYROOT = Path(__file__).resolve().parents[2]
EXTRACT = PYROOT / "tests" / "params" / "data" / "inputparser_extract.json"
TABLE = PYROOT / "docs" / "MATLAB_MAPPING.md"

# the extract's third-party functions, by the author lines of each file's header
THIRD_PARTY_FUNCTIONS = {
    "deskewFrame3D",  # an outside-author folder (llsmtools)
    "rotateFrame3D",  # the same
    "XR_loadConditionData3D",  # the header names an outside author (cmeAnalysis)
    "XR_runDetection3D",  # the same
    "XR_runTrackProcessing3D",  # the same
    "XR_runTracking3D",  # the same
    "runLifetimeAnalysis3D",  # the same
}

LLSM = "an outside-author folder (llsmtools)"


def _tool(name: str):
    spec = importlib.util.spec_from_file_location(name, PYROOT / "tools" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _fake(file: str, parser: dict | None = None, keep_unmatched: bool = False) -> MatlabFunction:
    return MatlabFunction(
        "XR_fake", file, ["data", "varargin"],
        [MatlabParam("required", "data", None, "@(x) isstruct(x) && fake_check_a(x)"),
         MatlabParam("optional", "mode", "'w'", None),
         MatlabParam("parameter", "Overwrite", "false", "@(x) islogical(x) && fake_check_b(x)")],
        keep_unmatched=keep_unmatched, parser=parser, notes=["a note"],
    )  # fmt: skip


# ----------------------------------------------------------------------------------------------
# the stored extract


def test_the_extract_keeps_only_the_interface_of_third_party_code(extract_json):
    """Every function the provenance rules call third-party has no validator text in the
    stored extract (nor in an mcc parser), and is marked; no other function is."""
    found = set()
    for name, e in extract_json["functions"].items():
        why = third_party(e["file"])
        if why is None:
            assert "third_party" not in e, name
            assert not any("withheld" in p for p in e["params"]), name
            continue
        found.add(name)
        assert e["third_party"] == why, name
        assert all(p["validator"] is None for p in e["params"]), name
        assert any(p.get("withheld") for p in e["params"]), name
        assert all(v is None for v in ((e["parser"] or {}).get("validators") or {}).values())
    assert found == THIRD_PARTY_FUNCTIONS


def test_the_interface_of_third_party_code_is_kept(extract):
    """Names, kinds, order and defaults stay (XR_runTracking3D as the example)."""
    fn = extract["XR_runTracking3D"]
    results = "arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0)"
    assert [(p.kind, p.name, p.default) for p in fn.params] == [
        ("required", "data", None),
        ("optional", "settings", "[]"),
        ("parameter", "Overwrite", "false"),
        ("parameter", "FileName", "'trackedFeatures'"),
        ("parameter", "Frames", "[]"),
        ("parameter", "DetectionFile", "'Detection3D.mat'"),
        ("parameter", "ResultsPath", results),
        ("parameter", "ParallelMethod", "'lock_distributed'"),
    ]
    assert all(p.withheld and p.validator is None for p in fn.params)


def test_the_stored_extract_is_in_the_makers_form(extract_json):
    """The file is what write_extract writes for its entries, so no line was edited by hand."""
    funcs = {k: MatlabFunction.from_json(v).to_json() for k, v in extract_json["functions"].items()}
    assert extract_text(funcs, extract_json["oracle_commit"]) == EXTRACT.read_text()


# ----------------------------------------------------------------------------------------------
# whose code


@pytest.mark.parametrize(
    ("file", "third"),
    [
        ("pointDetection/llsmtools/transforms/deskewFrame3D.m", True),
        ("pointDetection\\llsmtools\\transforms\\rotateFrame3D.m", True),
        ("pointDetection/XR_runDetection3D.m", True),
        ("pointDetection/runLifetimeAnalysis3D.m", True),
        ("pointDetection/rotateTracks3D.m", True),
        ("pointDetection/deskewData.m", True),
        ("imageProcessing/filters/x.m", True),
        ("third_parties/u-track3D/tracking/x.m", True),
        ("third_parties/sim3D/perdecomp_3D.m", True),
        ("third_parties/sim3D/imtranslate_function.m", True),
        ("third_parties/sim3D/simRecon.m", False),
        ("pointDetection/XR_runDetTrack3d.m", False),
        ("pointDetection/src/XR_amiraWriteTracks.m", False),
        ("microscopeDataProcessing/deskew_rotate/XR_deskewRotateFrame.m", False),
        ("utils/getImageSize.m", False),
    ],
)
def test_whose_code_it_is(file, third):
    assert (third_party(file) is not None) is third


@pytest.mark.parametrize(
    "file",
    [
        "pointDetection/XR_runDetection3D_frame_parallel.m",
        "pointDetection/src/XR_new_fake.m",
        "third_parties/normxcorr3/normxcorr3.m",
        "third_parties/sim3D_old/x.m",
    ],
)
def test_a_file_of_a_folder_with_outside_code_needs_a_decision(file):
    with pytest.raises(ProvenanceError, match="author lines of its header"):
        third_party(file)


def test_interface_only_withholds_the_validators_and_keeps_the_rest():
    parser = {
        "conversions": {"Overwrite": "str2num"},
        "defaults": {"data": None, "mode": "'w'", "Overwrite": "false"},
        "validators": {"data": "@(x) fake_check_c(x)", "mode": None, "Overwrite": "@islogical"},
    }
    fn = _fake("pointDetection/llsmtools/XR_fake.m", parser, keep_unmatched=True)
    got = interface_only(fn)
    assert [(p.kind, p.name, p.default) for p in got.params] == [
        (p.kind, p.name, p.default) for p in fn.params
    ]
    assert [(p.validator, p.withheld) for p in got.params] == [
        (None, True), (None, False), (None, True)
    ]  # fmt: skip
    assert got.parser == {**parser, "validators": {"data": None, "mode": None, "Overwrite": None}}
    assert (got.third_party, got.args, got.keep_unmatched, got.notes) == (
        LLSM, fn.args, True, ["a note"]
    )  # fmt: skip
    assert interface_only(got) == got
    assert fn.params[0].validator.startswith("@(x) isstruct(x)")  # the input is not changed
    lab = _fake("microscopeDataProcessing/XR_fake.m", parser)
    assert interface_only(lab) is lab


def test_the_marks_round_trip_and_the_labs_entries_keep_their_form():
    fn = interface_only(_fake("pointDetection/llsmtools/XR_fake.m"))
    d = fn.to_json()
    assert MatlabFunction.from_json(d) == fn
    assert d["third_party"] == LLSM
    assert [p.get("withheld") for p in d["params"]] == [True, None, True]
    lab = _fake("microscopeDataProcessing/XR_fake.m").to_json()
    assert "third_party" not in lab and not any("withheld" in p for p in lab["params"])
    assert MatlabFunction.from_json(lab) == _fake("microscopeDataProcessing/XR_fake.m")


# ----------------------------------------------------------------------------------------------
# check_registry


def _spec(validator: str | None, *, mode_kw: bool = False, overwrite: bool = False):
    reg = Registry()
    params = (
        Param("data", validator, convert="same"),
        Param("mode", None),
        Param("Overwrite", validator),
    )
    if mode_kw:

        @reg.register("XR_fake", *params)
        def fake_kw(data, *, mode="w", overwrite=overwrite):
            return data

    else:

        @reg.register("XR_fake", *params)
        def fake(data, mode="w", *, overwrite=overwrite):
            return data

    return reg.find("XR_fake")


def test_check_spec_compares_the_interface_of_third_party_code():
    outside = interface_only(_fake("pointDetection/llsmtools/XR_fake.m"))
    assert check_spec(_spec(None), outside) == []
    assert check_spec(_spec("@ischar"), outside) == []  # validators are not compared
    problems = check_spec(_spec(None, overwrite=True), outside)
    assert len(problems) == 1 and "default True, MATLAB has false" in problems[0], problems
    problems = check_spec(_spec(None, mode_kw=True), outside)
    assert len(problems) == 1 and "parameter here, optional in MATLAB" in problems[0], problems
    # the lab's code: validators are compared
    problems = check_spec(_spec(None), _fake("microscopeDataProcessing/XR_fake.m"))
    assert len(problems) == 2 and all("validator None, MATLAB has" in p for p in problems)


def test_a_withheld_validator_counts_for_text_by_position():
    """inputParser reads text at an optional argument without a validator as a parameter name
    (check_spec reports it where the function has no name-value parameters); an optional
    argument whose validator is withheld has one."""
    reg = Registry()

    @reg.register("XR_mode_fake", Param("x", None), Param("mode", None))
    def mode_fake(x, mode="w"):
        return mode

    def fn(withheld):
        return MatlabFunction(
            "XR_mode_fake", "pointDetection/llsmtools/XR_mode_fake.m", ["x", "mode"],
            [MatlabParam("required", "x", None, None),
             MatlabParam("optional", "mode", "'w'", None, withheld)],
        )  # fmt: skip

    spec = reg.find("XR_mode_fake")
    assert check_spec(spec, fn(True)) == []
    problems = check_spec(spec, fn(False))
    assert len(problems) == 1 and "takes it by position" in problems[0]


OUTSIDE_FUNC = """\
function [] = XR_outside_fake(data, varargin)
% XR_outside_fake: a stand-in for third-party code
%   Overwrite : true|false (default: false). Overwrite the results.
ip = inputParser;
ip.addRequired('data', @(x) isstruct(x) && fake_check_a(x));
ip.addOptional('mode', 'w');
ip.addParameter('Overwrite', false, @(x) islogical(x) && fake_check_b(x));
ip.parse(data, varargin{:});
end
"""


def _tree(root: Path) -> Path:
    files = {
        "pointDetection/llsmtools/XR_outside_fake.m": "XR_outside_fake",
        "mcc/parsers/XR_outside_fake_parser.m": "XR_outside_fake_parser",
        "pointDetection/XR_runDetTrack3d.m": "XR_runDetTrack3d",
        "microscopeDataProcessing/XR_lab_fake.m": "XR_lab_fake",
    }
    for rel, name in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(OUTSIDE_FUNC.replace("XR_outside_fake", name), newline="\n")
    return root


def test_write_extract_keeps_only_the_interface_of_third_party_code(tmp_path):
    root = _tree(tmp_path / "tree")
    out = tmp_path / "extract.json"
    write_extract(root, {"XR_outside_fake", "XR_runDetTrack3d", "XR_lab_fake"}, out)
    # the two lab files' validators, nothing of the outside function or its mcc parser
    assert out.read_text().count("fake_check_") == 4
    got = read_extract(out)
    o = got["XR_outside_fake"]
    assert o.third_party == LLSM
    assert [(p.kind, p.name, p.default, p.validator, p.withheld) for p in o.params] == [
        ("required", "data", None, None, True),
        ("optional", "mode", "'w'", None, False),
        ("parameter", "Overwrite", "false", None, True),
    ]
    assert o.parser["validators"] == {"data": None, "mode": None, "Overwrite": None}
    assert o.parser["defaults"] == {"data": None, "mode": "'w'", "Overwrite": "false"}
    for lab in ("XR_runDetTrack3d", "XR_lab_fake"):
        assert got[lab].third_party is None
        assert got[lab].params[2].validator == "@(x) islogical(x) && fake_check_b(x)"


def test_write_extract_refuses_a_file_without_a_decision(tmp_path):
    p = tmp_path / "tree" / "pointDetection" / "XR_new_fake.m"
    p.parent.mkdir(parents=True)
    p.write_text(OUTSIDE_FUNC.replace("XR_outside_fake", "XR_new_fake"), newline="\n")
    with pytest.raises(ProvenanceError, match="XR_new_fake.m"):
        write_extract(tmp_path / "tree", {"XR_new_fake"}, tmp_path / "extract.json")
    assert not (tmp_path / "extract.json").exists()


def test_the_skeleton_of_third_party_code_takes_nothing_from_its_file(tmp_path, capsys):
    root = _tree(tmp_path)
    argv = ["skeleton", "XR_outside_fake", "XR_lab_fake", "--matlab-root", str(root)]
    assert _tool("check_registry").main(argv) == 0
    outside, lab = capsys.readouterr().out.split("@register(")[1:]
    assert "fake_check" not in outside and "Overwrite the results" not in outside
    assert "stand-in" not in outside and "Param('Overwrite', None, '')" in outside
    assert "fake_check_b" in lab and "Overwrite the results" in lab


def test_the_tree_check_compares_the_interface_of_third_party_code(tmp_path, monkeypatch, capsys):
    import petakit5d.params as params

    root = _tree(tmp_path)
    reg = Registry()

    @reg.register(
        "XR_outside_fake",
        Param("data", None, convert="same"),
        Param("mode", None),
        Param("Overwrite", None),
    )
    def outside_fake(data, mode="w", *, overwrite=False):
        return data

    monkeypatch.setattr(params, "REGISTRY", reg)
    assert _tool("check_registry").main(["--matlab-root", str(root)]) == 0
    out = capsys.readouterr().out
    assert "1 in step" in out and "fake_check" not in out


# ----------------------------------------------------------------------------------------------
# the mapping table


def test_the_table_marks_withheld_validators_of_registered_and_planned_functions():
    outside = interface_only(_fake("pointDetection/llsmtools/XR_fake.m"))
    lab = _fake("microscopeDataProcessing/XR_fake.m")
    reg = Registry()

    @reg.register(
        "XR_fake",
        Param("data", "@(x) fake_spec_check(x)", convert="same"),
        Param("mode", None),
        Param("Overwrite", "@islogical"),
        Param("numThreads", "@isnumeric", python_only=True),
    )
    def fake(data, mode="w", *, overwrite=False, num_threads=None):
        return data

    spec = reg.find("XR_fake")
    sec = "\n".join(_spec_section(spec, outside))
    assert sec.count(WITHHELD) == 2 and "fake_spec_check" not in sec and "`@islogical`" not in sec
    assert "`@isnumeric`" in sec  # a Python-only parameter's own validator
    assert f"Third-party MATLAB code: {LLSM}. Its validators are not shown." in sec
    sec = "\n".join(_spec_section(spec, lab))
    assert WITHHELD not in sec and "fake_spec_check" in sec and "Third-party" not in sec
    sec = "\n".join(_planned_section("XR_fake", "petakit5d.x.fake", outside))
    assert sec.count(WITHHELD) == 2 and "fake_check" not in sec and "Third-party MATLAB" in sec
    assert "fake_check_b" in "\n".join(_planned_section("XR_fake", "petakit5d.x.fake", lab))


def _sections(text: str) -> dict[str, str]:
    heading = re.compile(r"(?ms)^### `(\w+)`.*?(?=^### |^## |\Z)")
    return {m.group(1): m.group(0) for m in heading.finditer(text)}


def _validator_cells(section: str) -> list[tuple[str, str]]:
    """(validator cell, row) of each parameter row of a section."""
    lines = section.splitlines()
    head = next(r for r in lines if r.startswith("| MATLAB |"))
    col = [c.strip() for c in head.strip("|").split("|")].index("Validator")
    out = []
    for r in lines:
        if r.startswith("| `"):
            cells = re.split(r"(?<!\\)\|", r.strip()[1:-1])
            out.append((cells[col].strip(), r))
    return out


def test_the_table_shows_no_validator_of_third_party_code(extract):
    """docs/MATLAB_MAPPING.md: in each third-party function's section the note, and the
    marker in place of every withheld validator (registered or planned); nowhere else."""
    sections = _sections(TABLE.read_text())
    for name, sec in sections.items():
        fn = extract.get(name)
        if name not in THIRD_PARTY_FUNCTIONS:
            assert WITHHELD not in sec and "Third-party MATLAB code" not in sec, name
            continue
        assert f"Third-party MATLAB code: {fn.third_party}." in sec, name
        withheld = {p.name for p in fn.params if p.withheld}
        cells = _validator_cells(sec)
        assert len(cells) >= len(fn.params), name
        for cell, row in cells:
            matlab = row.split("`")[1]
            if "Python only" in row:
                continue
            assert cell == (WITHHELD if matlab in withheld else ""), (name, matlab)
    assert THIRD_PARTY_FUNCTIONS <= set(sections)
