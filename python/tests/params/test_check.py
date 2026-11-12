"""check_registry against the stored MATLAB extract; the mapping table is up to date."""

from __future__ import annotations

import copy
import importlib.util
import os
from pathlib import Path

import pytest

from petakit5d.params.check import check_all, check_spec, read_extract, same_python
from petakit5d.params.matlab_source import extract_tree, mcc_mode

PYROOT = Path(__file__).resolve().parents[2]


def _tool(name: str):
    spec = importlib.util.spec_from_file_location(name, PYROOT / "tools" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _mutated(spec, name: str, **changes):
    s = copy.copy(spec)
    s._params = [p.replace(**changes) if p.matlab == name else p for p in spec.params]
    return s


def test_worked_example_is_in_step(sample_crop, extract):
    spec = sample_crop.REG.find("XR_crop_dataset")
    assert check_spec(spec, extract["XR_crop_dataset"]) == []


@pytest.mark.parametrize(
    ("name", "changes", "expect"),
    [
        ("batchSize", {"default": (512, 512, 512)}, "default (512, 512, 512)"),
        ("pad", {"default": 0}, "default 0"),
        ("cpusPerTask", {"default": 2.5}, "default 2.5"),
        ("channelPatterns", {"default": ("CamA_ch0",)}, "default"),
        ("uuid", {"default": None}, "default None"),
        ("pad", {"validator": "@(x) islogical(x)"}, "validator"),
        ("cropType", {"matlab": "croptype"}, "spells it"),
        ("pad", {"kind": "optional"}, "optional here"),
        ("batchSize", {"mcc": "none"}, "mcc parser converts"),
        ("jobLogDir", {"matlab_default": "'../logs'"}, "default text"),
    ],
)
def test_drift_is_found(sample_crop, extract, name, changes, expect):
    spec = _mutated(sample_crop.REG.find("XR_crop_dataset"), name, **changes)
    problems = check_spec(spec, extract["XR_crop_dataset"])
    assert any(expect in p for p in problems), problems


def test_missing_and_extra_parameters_are_found(sample_crop, extract):
    spec = sample_crop.REG.find("XR_crop_dataset")
    s = copy.copy(spec)
    s._params = spec.params[:-1]
    assert "missing ['configFile']" in check_spec(s, extract["XR_crop_dataset"])[0]
    s._params = list(reversed(spec.params))
    assert "order differs" in check_spec(s, extract["XR_crop_dataset"])[0]
    assert check_spec(spec, None) != []


def test_differs_skips_a_site_specific_default(sample_crop, extract):
    spec = _mutated(
        sample_crop.REG.find("XR_crop_dataset"), "jobLogDir", default="/x", differs="site"
    )
    assert check_spec(spec, extract["XR_crop_dataset"]) == []


def test_generator_rule_matches_every_mcc_parser(extract):
    """The command-line conversion rule (mcc_mode) gives what every parser file does."""
    n = 0
    for name, fn in extract.items():
        if not fn.parser:
            continue
        for p in fn.params:
            if isinstance(p.validator, str):
                assert mcc_mode(p.validator) == fn.parser["conversions"].get(p.name, "none"), (
                    name,
                    p.name,
                )
                n += 1
    assert n > 1000


def test_parsers_have_the_functions_defaults(extract):
    """No MATLAB drift between the mcc parsers and their functions (the CLI relies on it)."""
    from petakit5d.params.check import _norm

    for name, fn in extract.items():
        if not fn.parser:
            continue
        for p in fn.params:
            pd = fn.parser["defaults"].get(p.name)
            if p.kind != "required" and isinstance(pd, str):
                assert _norm(pd) == _norm(p.default), (name, p.name)


def test_same_python():
    assert same_python((1, 2), (1.0, 2.0))
    assert not same_python(True, 1) and not same_python(0, False)
    assert same_python(float("nan"), float("nan"))
    assert not same_python(("a",), "a")


def test_check_all_and_the_tools(sample_crop, extract):
    spec = sample_crop.REG.find("XR_crop_dataset")
    assert check_all([spec], extract) == {"XR_crop_dataset": []}
    assert _tool("check_registry").main([]) == 0
    assert _tool("gen_mapping").main(["--check"]) == 0, "run python tools/gen_mapping.py"


def test_skeleton_registers_and_checks(extract, tmp_path):
    from petakit5d.params import Param, Registry
    from petakit5d.params.skeleton import skeleton

    reg = Registry()
    fn = extract["XR_crop_dataset"]
    code = skeleton(fn).replace("@register(", "@reg.register(")
    ns = {"reg": reg, "Param": Param}
    exec(compile(code, "skeleton", "exec"), ns)  # noqa: S102 - code we generated
    spec = reg.find("XR_crop_dataset")
    assert check_spec(spec, fn) == []


def test_stored_extract_is_the_oracle_trees(tmp_path):
    """On the cluster, the stored extract equals a fresh one from the oracle tree."""
    root = os.environ.get("PETAKIT5D_ORACLE_DIR")
    if not root or not Path(root, "mcc").is_dir():
        pytest.skip("no oracle tree (PETAKIT5D_ORACLE_DIR)")
    from petakit5d.params.check import site_pattern, write_extract

    stored = read_extract(PYROOT / "tests" / "params" / "data" / "inputparser_extract.json")
    out = tmp_path / "extract.json"
    write_extract(root, set(stored), out, "", site_pattern())
    fresh = read_extract(out)
    assert {k: v.to_json() for k, v in fresh.items()} == {k: v.to_json() for k, v in stored.items()}
    tree = extract_tree(root, {"XR_crop_dataset"})
    assert tree["XR_crop_dataset"].params == stored["XR_crop_dataset"].params


@pytest.mark.parametrize(
    ("text", "site"),
    [
        ("'/mnt/projects/x.sh'", True),
        ("'~/Downloads/'", True),
        ("'C:\\MATLAB'", True),
        ("'-p cpu_queue --qos normal -n1'", True),
        ("'-A lab_acct'", True),
        ("'--partition=cpu'", True),
        ("'../job_logs'", False),
        ("'DSR/'", False),
        ("[dataPath, '/FF/']", False),
        ("pi/2", False),
        ("1/2/3", False),
        ("[1 -2]", False),
        ("'-x,y,z'", False),
        ("{'CamA_ch0', 'CamB_ch0'}", False),
    ],
)
def test_site_pattern(text, site, monkeypatch):
    from petakit5d.params.check import SITE_PATTERN_ENV, site_pattern

    monkeypatch.delenv(SITE_PATTERN_ENV, raising=False)
    assert bool(site_pattern().search(text)) is site
    monkeypatch.setenv(SITE_PATTERN_ENV, "CamB")
    assert site_pattern().search("{'CamA_ch0', 'CamB_ch0'}")


SITE_FUNC = """\
function [] = XR_site_fake(dataPaths, varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @ischar);
ip.addParameter('SlurmParam', '-p cpu_queue --qos normal -n1', @ischar);
ip.addParameter('MCCMasterStr', '/mnt/projects/x.sh', @ischar);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('MCRParam', 'C:\\MATLAB\\R2024b', @ischar);
ip.addParameter('outDir', 'DSR/', @ischar);
ip.parse(dataPaths, varargin{:});
end
"""


def test_site_values_are_masked(tmp_path):
    """Site-specific MATLAB defaults are stored as hashes and need differs= in a spec."""
    from petakit5d.params import Param, Registry
    from petakit5d.params.check import is_site_value, site_pattern, write_extract
    from petakit5d.params.mapping import _planned_section

    (tmp_path / "XR_site_fake.m").write_text(SITE_FUNC)
    out = tmp_path / "extract.json"
    write_extract(tmp_path, {"XR_site_fake"}, out, "", site_pattern())
    text = out.read_text()
    for value in ("cpu_queue", "/mnt/projects", "R2024b"):
        assert value not in text
    fn = read_extract(out)["XR_site_fake"]
    assert [is_site_value(p.default) for p in fn.params] == [False, True, True, False, True, False]
    table = "\n".join(_planned_section("XR_site_fake", "petakit5d.x.site_fake", fn))
    assert "<site-specific>" in table and "/mnt/projects" not in table

    def make(differs):
        reg = Registry()

        @reg.register(
            "XR_site_fake",
            Param("dataPaths", "@ischar"),
            Param("SlurmParam", "@ischar", differs=differs),
            Param("MCCMasterStr", "@ischar", differs=differs),
            Param("jobLogDir", "@ischar"),
            Param("MCRParam", "@ischar", differs=differs),
            Param("outDir", "@ischar"),
        )
        def site_fake(
            data_paths, *, slurm_param="", mcc_master_str="", job_log_dir="../job_logs",
            mcr_param="", out_dir="DSR/",
        ):  # fmt: skip
            return locals()

        return reg.find("XR_site_fake")

    assert check_spec(make("site-specific MATLAB default"), fn) == []
    problems = check_spec(make(None), fn)
    assert len(problems) == 3 and all("site-specific" in p for p in problems)
    # a masked text that changes on the MATLAB side is still found through its hash
    fn.params[0].validator = "<site value sha256:0000000000000000>"
    assert any("validator" in p for p in check_spec(make("site"), fn))
