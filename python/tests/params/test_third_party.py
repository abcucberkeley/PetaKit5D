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
    "XR_runDetection3D_frame_parallel",  # the same
    "XR_runTrackProcessing3D",  # the same
    "XR_runTracking3D",  # the same
    "runLifetimeAnalysis3D",  # the same
    "plotIntensityCohorts",  # the same (CMEAnalysis' copy, the one that runs)
    "loadTrackSettings",  # the same (CMEAnalysis' copy, the one that runs)
    "deskewData",  # its header names an outside author
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
        # (loadTrackSettings validates none of its options: nothing to withhold)
        assert any(p.get("withheld") for p in e["params"]) or name == "loadTrackSettings", name
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
        ("pointDetection/XR_runDetection3D_frame_parallel.m", True),
        ("pointDetection/runLifetimeAnalysis3D.m", True),
        ("pointDetection/rotateTracks3D.m", True),
        ("pointDetection/deskewData.m", True),
        ("pointDetection/CMEAnalysis_104_XR/software/plotIntensityCohorts.m", True),
        ("pointDetection/llsmtools/cmeAnalysis3D/plotIntensityCohorts.m", True),
        ("pointDetection/CMEAnalysis_104_XR/software/loadTrackSettings.m", True),
        ("pointDetection/llsmtools/cmeAnalysis3D/deskewData.m", True),
        ("pointDetection/XR_correctXZoffsetData3D.m", False),
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
        "pointDetection/XR_new_fake_frame_parallel.m",
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
    # validators are not compared, and a validator text there is refused
    problems = check_spec(_spec("@ischar"), outside)
    assert len(problems) == 2 and all("instead of a validator" in p for p in problems), problems
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
            if matlab not in withheld:
                assert cell == "", (name, matlab)
            else:  # the marker, with the port's own check where the spec has one
                assert cell == WITHHELD or cell.startswith(f"{WITHHELD} the port's: `"), (
                    name,
                    matlab,
                )
    assert THIRD_PARTY_FUNCTIONS <= set(sections)


# ----------------------------------------------------------------------------------------------
# plotIntensityCohorts, registered by its interface


def test_a_files_own_decision_wins_over_its_folders():
    from petakit5d.params.provenance import THIRD_PARTY

    cme = "pointDetection/CMEAnalysis_104_XR/"
    own = cme + "software/plotIntensityCohorts.m"
    assert third_party(own) == THIRD_PARTY[own] != THIRD_PARTY[cme]
    assert third_party(cme + "software/loadTracks.m") == THIRD_PARTY[cme]


def test_plot_intensity_cohorts_is_registered_by_its_interface(extract):
    """MATLAB's names, kinds, order and defaults (check_registry's rules), no validator text,
    and a check of the port's own for each withheld validator but data's (MATLAB took every
    probe of data at its parse: tests/params/test_third_party_checks.py)"""
    from petakit5d.params import find
    from petakit5d.params.check import check_signature_defaults

    spec = find("plotIntensityCohorts")
    fn = extract["plotIntensityCohorts"]
    assert spec.python_name == "petakit5d.detection.lifetime.plot_intensity_cohorts"
    assert check_spec(spec, fn) == [] and check_signature_defaults(spec) == []
    assert [(p.kind, p.matlab) for p in spec.params] == [(m.kind, m.name) for m in fn.params]
    assert all(p.validator is None for p in spec.params)
    withheld = {m.name for m in fn.params if m.withheld}
    assert {p.matlab for p in spec.params if p.accepts} == withheld - {"data"}


# a value other than the default that the port's check takes, per parameter taken at its
# default only
NOT_DEFAULT = {
    "ch": (1.0,),
    "overwrite": True,
    "cohort_bounds_s": (10.0, 20.0),
    "show_variation": False,
    "fill_mode": "sem",
    "front_layer": True,
    "show_background": True,
    "rescale": False,
    "rescaling_reference": "MED",
    "scale_slave_channel": False,
    "scaling_factor": (1.0, 2.0),
    "max_intensity_threshold": 1.0,
    "exclude_visitors": False,
    "slave_name": ("ch2",),
    "channel_names": ("a",),
    "line_style": "--",
    "hues": (0.5,),
    "colormap": (1.0, 0.0, 0.0),
    "colormap_fill": (1.0, 0.0, 0.0),
    "display_mode": "print",
    "display_all": True,
    "track_index": (1.0,),
    "cutoff_f": 3,
    "alpha": 0.01,
    "y_tick": (0.0, 1.0),
    "y_lim": (0.0, 1.0),
    "parent": 1.0,
    "remove_outliers": True,
    "plot_x_label": False,
    "show_pct": False,
    "show_stats": False,
    "avg_fun": len,
    "lft_data_name": "other.mat",
    "amplitude_correction": (1.0,),
    "align": "right",
}


@pytest.mark.parametrize("name", list(NOT_DEFAULT))
def test_plot_intensity_cohorts_takes_its_other_options_at_their_defaults_only(name):
    from petakit5d.detection.lifetime import plot_intensity_cohorts

    with pytest.raises(NotImplementedError, match="only the default is ported"):
        plot_intensity_cohorts({"source": "movie_missing/"}, **{name: NOT_DEFAULT[name]})


def test_every_option_of_plot_intensity_cohorts_but_its_ported_ones_is_guarded():
    from petakit5d.params import find

    ported = {"data", "analysis_path"}
    assert {p.py for p in find("plotIntensityCohorts").params} - ported == set(NOT_DEFAULT)


# ----------------------------------------------------------------------------------------------
# deskewData and loadTrackSettings (third-party), XR_correctXZoffsetData3D (the lab's)


@pytest.mark.parametrize(
    ("name", "python", "third"),
    [
        ("deskewData", "petakit5d.detection.wrapper.deskew_data", True),
        ("loadTrackSettings", "petakit5d.detection.tracking.settings.track_settings", True),
        ("XR_correctXZoffsetData3D", "petakit5d.detection.offsets.correct_xz_offset_data", False),
    ],
)
def test_the_detection_helpers_are_registered(extract, name, python, third):
    """MATLAB's names, kinds, order and defaults (check_registry's rules); third-party code without
    validator texts and with a check of the port's own for each withheld validator; the lab's code
    with its validators"""
    from petakit5d.params import find
    from petakit5d.params.check import check_signature_defaults

    spec = find(name)
    fn = extract[name]
    assert spec.python_name == python
    assert check_spec(spec, fn) == [] and check_signature_defaults(spec) == []
    assert (fn.third_party is not None) is third
    if third:
        assert all(p.validator is None for p in spec.params)
        assert {p.matlab for p in spec.params if p.accepts} == {
            m.name for m in fn.params if m.withheld
        }
    else:
        ours = [p.validator for p in spec.params if not p.python_only]
        assert ours == [m.validator for m in fn.params]


DESKEW_NOT_DEFAULT = {
    "movie_selector": "Ex",
    "crop": True,
    "rotate": True,
    "check_frame_mismatch": True,
    "load_settings": True,
    "llff_correction": True,
    "lower_limit": 0.5,
    **{f"ls_image_ch_{c}": "ls.tif" for c in range(1, 5)},
    **{f"background_ch_{c}": "bg.tif" for c in range(1, 5)},
    "save_16bit": True,
    "zoffset_correction": True,
}


@pytest.mark.parametrize("name", list(DESKEW_NOT_DEFAULT))
def test_deskew_data_takes_its_other_options_at_their_defaults_only(name):
    from petakit5d.detection.condition import Condition
    from petakit5d.detection.wrapper import deskew_data

    with pytest.raises(NotImplementedError, match="only the default is ported"):
        deskew_data(Condition([]), **{name: DESKEW_NOT_DEFAULT[name]})


def test_every_option_of_deskew_data_but_its_ported_ones_is_guarded():
    from petakit5d.params import find

    ported = {"data", "reverse", "overwrite", "skew_angle", "s_cmos_camera_flip", "aname"}
    ported.add("num_threads")
    assert {p.py for p in find("deskewData").params} - ported == set(DESKEW_NOT_DEFAULT)


def test_deskew_data_without_data_is_matlabs_interactive_load():
    from petakit5d.detection.condition import Condition
    from petakit5d.detection.wrapper import WrapperError, deskew_data

    with pytest.raises(WrapperError) as err:
        deskew_data()
    assert err.value.identifier == "MATLAB:services:MissingRequiredCapability"
    assert tuple(deskew_data(Condition([]), True)) == ()  # data, then Reverse by position


@pytest.mark.parametrize(
    ("name", "value"),
    [("gap_radius", (3.0, 6.0)), ("link_radius", (2.0, 4.0)), ("cost_type", "other")],
)
def test_track_settings_takes_its_other_options_at_their_defaults_only(name, value):
    from petakit5d.detection.tracking import track_settings

    with pytest.raises(NotImplementedError, match="only the default is ported"):
        track_settings(**{name: value})


def test_track_settings_takes_matlabs_name_value_options():
    """loadTrackSettings has name-value options only; Radius [] (None) is the default radius"""
    import numpy as np

    from petakit5d.detection.tracking import track_settings
    from petakit5d.params import find

    np.testing.assert_equal(track_settings(), track_settings(radius=(5, 10)))
    with pytest.raises(TypeError):
        track_settings((3, 6), 2)
    got = find("loadTrackSettings").call_matlab("Radius", [3.0, 6.0], "MaxGapLength", 1.0)
    np.testing.assert_equal(got, track_settings(radius=(3, 6), max_gap_length=1))
