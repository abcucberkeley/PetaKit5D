"""The mutation campaign's tool (tools/mutation): its sites, the mutants it makes and applies,
its sampling, the test selection from coverage contexts and a run end to end on a small tree."""

import importlib.util
import json
import os
import sys
import textwrap
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[2] / "tools" / "mutation"


def _load(name):
    # (a module of its own name in sys.modules, as dataclasses need; the tool imports its
    # mutants module from its folder)
    spec = importlib.util.spec_from_file_location(f"pk5d_mutation_{name}", TOOL / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    sys.path.insert(0, str(TOOL))
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.path.remove(str(TOOL))
    return mod


C = _load("campaign")
M = C.M
CPP = _load("cpp")


def test_cpp_sites_skip_comments_strings_templates_and_the_preprocessor():
    code = textwrap.dedent(
        """\
        #include <vector>  // a < b
        #define N 4
        /* x + y
           a < b */
        std::vector<int> v(3);
        if (a < b && c != 0) s += "x < y" [0];
        bool f = true;  // 7
        """
    )
    got = {(s.line, s.op, s.old, s.new) for s in CPP.find_cpp_sites(code)}
    assert got == {
        (5, "num", "3", "4"),
        (6, "cmp", "<", "<="),
        (6, "bool", "&&", "||"),
        (6, "cmp", "!=", "=="),
        (6, "num", "0", "1"),
        (6, "arith", "+=", "-="),
        (7, "bool", "true", "false"),
    }
    site = next(s for s in CPP.find_cpp_sites(code) if s.op == "cmp")
    assert M.apply_site(code, site).split("\n")[5].startswith("if (a <= b")


def test_cpp_groups_are_found_from_their_loaders(tmp_path):
    root = tmp_path / "py"
    (root / "src" / "petakit5d" / "k").mkdir(parents=True)
    (root / "src" / "petakit5d" / "k" / "m.py").write_text(
        'def f():\n    from petakit5d._ext import load\n\n    return load(\n        "mip"\n    )\n'
    )
    assert CPP.loaders(root)["mip"] == ["src/petakit5d/k/m.py"]
    assert CPP.group_of("cpp/mip/module.cpp") == "mip"


CODE = textwrap.dedent(
    '''\
    """a module"""
    from typing import TYPE_CHECKING

    if TYPE_CHECKING:
        LIMIT = 3

    __all__ = ["f", "g"]


    def f(a: int, b: int = 2) -> tuple[int, int]:
        """doc 1 < 2"""
        assert a > 0
        if a < b and b != 0:
            a += 1
        print(f"{a + 1} {b}", a - 1)
        out = (a * b, max(a, b))
        return out


    def g(x, y):
        if isinstance(x, (int, float)):
            return not x
        y.v = round(x // 2)
        return y is None
    '''
)


def _sites(code=CODE, ops=M.OPERATORS):
    return {(s.line, s.op, s.old, s.new) for s in M.find_sites(code, ops)}


def test_sites_of_each_operator():
    s = _sites()
    assert {
        (13, "cmp", "<", "<="),
        (13, "cmp", "!=", "=="),
        (13, "bool", "and", "or"),
        (13, "cond", "a < b and b != 0", "not (a < b and b != 0)"),
        (13, "num", "0", "1"),
        (14, "arith", "+=", "-="),
        (14, "num", "1", "2"),
        (14, "stmt", "a += 1", "pass"),
        (16, "arith", "*", "/"),
        (16, "call", "max", "min"),
        (17, "stmt", "return out", "pass"),
        (22, "bool", "not ", ""),
        (23, "stmt", "y.v = round(x // 2)", "pass"),
        (23, "call", "round(x // 2)", "(x // 2)"),
        (24, "cmp", "is", "is not"),
    } <= s
    # (a tuple of items that are not numbers or names is not swapped)
    assert not [x for x in s if x[0] == 16 and x[1] == "swap"]


def test_no_sites_in_text_annotations_checks_and_prints():
    s = _sites()
    # docstrings, TYPE_CHECKING, __all__, assert, print
    assert not {x[0] for x in s} & {1, 4, 5, 7, 11, 12, 15}
    # annotations: only the default value; the types of isinstance are not swapped
    assert {x for x in s if x[0] == 10} == {(10, "num", "2", "3")}
    assert {x[1] for x in s if x[0] == 21} == {"cond"}


def test_swaps_of_simple_items():
    s = _sites("p = (2, 0, 1)\nq, r = r, q\n", ["swap"])
    assert s == {
        (1, "swap", "2, 0", "0, 2"),
        (1, "swap", "0, 1", "1, 0"),
        (2, "swap", "q, r", "r, q"),
        (2, "swap", "r, q", "q, r"),
    }


def test_byte_columns_and_application():
    code = 'x = "µ" + str(1 < 2)\n'
    (site,) = M.find_sites(code, ["cmp"])
    assert M.apply_site(code, site) == 'x = "µ" + str(1 <= 2)\n'
    with pytest.raises(ValueError):
        M.apply_site(code.replace("<", ">"), site)


def test_valid_mutants_compile_and_change_the_tree():
    code = "def f(x):\n    return x\n"
    (site,) = M.find_sites(code, ["stmt"])
    assert M.valid(code, site)
    assert not M.valid(code, M.Site(2, 4, 10, "stmt", "return", "retur"))
    assert not M.valid(code, M.Site(2, 11, 12, "num", "x", "x"))


def test_sampling_is_deterministic_and_takes_distinct_lines():
    texts = {"a.py": CODE, "b.py": "def h(z):\n    return z + 1 if z > 2 else z - 3\n"}
    sites = {f: M.find_sites(t) for f, t in texts.items()}
    one = list(M.sample_sites(sites, texts, 6, "seed"))
    assert one == list(M.sample_sites(sites, texts, 6, "seed"))
    assert len(one) == 6 and len({(f, s.line) for f, s in one}) == 6
    assert list(M.sample_sites(sites, texts, 6, "other")) != one


def test_a_mutant_is_applied_only_to_its_version_of_the_file():
    site = next(s for s in M.find_sites(CODE) if s.op == "cmp")
    m = M.make_mutant("x-001", "x", "src/petakit5d/x/m.py", CODE, site)
    assert m.mut_line.count("<=") == 1
    assert M.apply_mutant(CODE, m).splitlines()[site.line - 1] == m.mut_line
    assert M.Mutant.from_record(json.loads(json.dumps(m.record()))) == m
    with pytest.raises(ValueError, match="not the version"):
        M.apply_mutant(CODE + "\n", m)


def test_areas_and_junit_node_ids(tmp_path):
    assert C.area_of("src/petakit5d/stitch/volume.py") == "stitch"
    assert C.area_of("src/petakit5d/__init__.py") == "package"
    assert C.area_of("src/petakit5d/_ext/__init__.py") == "package"
    (tmp_path / "tests" / "io").mkdir(parents=True)
    (tmp_path / "tests" / "io" / "test_tiff.py").write_text("")
    mod = "tests/io/test_tiff.py"
    assert C.junit_nodeid("tests.io.test_tiff", "a[x.y]", tmp_path) == f"{mod}::a[x.y]"
    assert C.junit_nodeid("tests.io.test_tiff.TestB", "c", tmp_path) == f"{mod}::TestB::c"


def _tree(tmp_path):
    root = tmp_path / "py"
    (root / "src" / "petakit5d" / "foo").mkdir(parents=True)
    (root / "tests" / "foo").mkdir(parents=True)
    (root / "src" / "petakit5d" / "__init__.py").write_text("")
    (root / "src" / "petakit5d" / "foo" / "__init__.py").write_text("")
    (root / "src" / "petakit5d" / "foo" / "bar.py").write_text(
        "def clip(x, hi):\n    if x > hi:\n        return hi\n    return x\n\n\n"
        "def unused(x):\n    return x * 2\n"
    )
    (root / "tests" / "foo" / "test_bar.py").write_text(
        "from petakit5d.foo.bar import clip\n\n\ndef test_clip():\n    assert clip(5, 3) == 3\n"
        "\n\ndef test_clip_below():\n    assert clip(1, 3) == 1\n"
    )
    (root / "pyproject.toml").write_text("[tool.pytest.ini_options]\naddopts = []\n")
    return root


T1, T2 = "tests/foo/test_bar.py::test_clip", "tests/foo/test_bar.py::test_clip_below"


def test_generate_puts_the_prefix_before_the_ids(tmp_path):
    root = _tree(tmp_path)
    out = tmp_path / "m.jsonl"
    args = ["--root", str(root), "--out", str(out), "--per-area", "2", "--prefix", "p2/"]
    assert C.main(["generate", *args]) == 0
    assert [m["id"] for m in C.read_jsonl(out)] == ["p2/foo-001", "p2/foo-002"]


@pytest.mark.skipif(os.name == "nt", reason="the campaign runs on Linux (process groups)")
def test_a_campaign_end_to_end(tmp_path):
    root = _tree(tmp_path)
    out = tmp_path / "m.jsonl"
    assert C.main(["generate", "--root", str(root), "--out", str(out), "--per-area", "9"]) == 0
    mutants = C.read_jsonl(out)
    assert {m["line"] for m in mutants} == {2, 3, 4, 8}
    # the plan as the coverage run would give it: the tests that run each line
    tests = {2: [T1, T2], 3: [T1], 4: [T2], 8: []}
    plan = tmp_path / "plan.jsonl"
    with open(plan, "w") as fh:
        for m in mutants:
            t = tests[m["line"]]
            p = dict(m, covered=bool(t), area_tests=t, area_seconds=1.0, suite_tests=[])
            fh.write(json.dumps(dict(p, suite_seconds=0.0)) + "\n")
    res = tmp_path / "res.jsonl"
    run = ["run", "--root", str(root), "--plan", str(plan), "--out", str(res)]
    assert C.main([*run, "--work", str(tmp_path / "w"), "--parallel", "2"]) == 0
    st = {r["line"]: r for r in C.read_jsonl(res)}
    assert st[8]["status"] == "not_covered"
    assert st[3]["status"] == "killed" and st[3]["phases"]["area"]["first"] == T1
    assert st[4]["status"] == "killed"
    # (a rerun skips what is done)
    assert C.main([*run, "--work", str(tmp_path / "w")]) == 0
    assert len(C.read_jsonl(res)) == 4
    md = tmp_path / "report.md"
    disp = tmp_path / "d.json"
    unused = next(m["id"] for m in mutants if m["line"] == 8)
    disp.write_text(json.dumps({unused: {"kind": "unreachable", "why": "never called"}}))
    rep = ["report", "--plan", str(plan), "--results", str(res), "--dispositions", str(disp)]
    assert C.main([*rep, "--out", str(md)]) == 0
    text = md.read_text()
    assert "| foo | 4 |" in text and "unreachable: never called" in text


def test_report_leaves_out_what_does_not_compile_and_names_the_new_tests(tmp_path):
    def jsonl(name, rows):
        (tmp_path / name).write_text("".join(json.dumps(r) + "\n" for r in rows))
        return str(tmp_path / name)

    base = {"area": "cpp/x", "file": "cpp/x/a.cpp", "line": 1, "op": "num", "old": "1", "new": "2"}
    plan = jsonl("plan.jsonl", [dict(base, id=f"cpp/x-00{i}") for i in (1, 2, 3)])
    res = jsonl(
        "res.jsonl",
        [
            {"id": "cpp/x-001", "status": "uncompiled"},
            {"id": "cpp/x-002", "status": "survived"},
            {"id": "cpp/x-003", "status": "killed", "by": "area"},
        ],
    )
    new = {"id": "cpp/x-002", "status": "killed", "phases": {"area": {"first": "t.py::t"}}}
    md = tmp_path / "r.md"
    rep = ["report", "--plan", plan, "--results", res, "--rechecks", jsonl("re.jsonl", [new])]
    assert C.main([*rep, "--out", str(md)]) == 0
    text = md.read_text()
    assert "| cpp/x | 2 | 1 | 0 | 50% | 0 | 0 | 1 | 0 | 0 | 100% | 100% |" in text
    assert "don't compile: cpp/x-001." in text and "test: `t.py::t`" in text


def test_plan_selects_the_tests_that_run_each_line(tmp_path):
    coverage = pytest.importorskip("coverage")
    root = _tree(tmp_path)
    bar = str(root / "src" / "petakit5d" / "foo" / "bar.py")
    data = coverage.CoverageData(str(tmp_path / "cov.db"))
    for ctx, lines in ((f"{T1}|run", [1, 2, 3]), (f"{T2}|run", [2, 4]), ("", [1, 7])):
        data.set_context(ctx)
        data.add_lines({bar: lines})
    data.write()
    junit = tmp_path / "junit.xml"
    cases = "".join(
        f'<testcase classname="tests.foo.test_bar" name="{n}" time="0.01"/>'
        for n in ("test_clip", "test_clip_below")
    )
    junit.write_text(f"<testsuites><testsuite>{cases}</testsuite></testsuites>")
    out = tmp_path / "m.jsonl"
    assert C.main(["generate", "--root", str(root), "--out", str(out), "--per-area", "9"]) == 0
    plan = tmp_path / "plan.jsonl"
    args = ["--mutants", str(out), "--coverage", data.data_filename(), "--junit", str(junit)]
    assert C.main(["plan", "--root", str(root), *args, "--out", str(plan)]) == 0
    by_line = {p["line"]: p for p in C.read_jsonl(plan)}
    assert by_line[3]["area_tests"] == [T1]
    assert sorted(by_line[2]["area_tests"]) == [T1, T2]
    assert by_line[4]["area_tests"] == [T2]
    assert not by_line[8]["covered"] and by_line[8]["area_tests"] == []


def test_plan_gives_signature_lines_the_header_tests(tmp_path):
    coverage = pytest.importorskip("coverage")
    root = tmp_path / "py"
    for d in ("src/petakit5d/foo", "tests/foo", "tests/params"):
        (root / d).mkdir(parents=True)
    for f in ("src/petakit5d/__init__.py", "src/petakit5d/foo/__init__.py"):
        (root / f).write_text("")
    for f in ("tests/foo/test_f.py", "tests/params/test_check.py"):
        (root / f).write_text("")
    mod = root / "src" / "petakit5d" / "foo" / "f.py"
    mod.write_text("def f(x, k=2):\n    return x * k\n")
    data = coverage.CoverageData(str(tmp_path / "cov.db"))
    for ctx, lines in (("tests/foo/test_f.py::test_f|run", [2]), ("", [1])):
        data.set_context(ctx)
        data.add_lines({str(mod): lines})
    data.write()
    junit = tmp_path / "junit.xml"
    cases = '<testcase classname="tests.foo.test_f" name="test_f" time="0.1"/>'
    cases += '<testcase classname="tests.params.test_check" name="test_registry" time="0.2"/>'
    junit.write_text(f"<testsuites><testsuite>{cases}</testsuite></testsuites>")
    out, plan = tmp_path / "m.jsonl", tmp_path / "plan.jsonl"
    assert C.main(["generate", "--root", str(root), "--out", str(out), "--per-area", "9"]) == 0
    args = ["--mutants", str(out), "--coverage", data.data_filename(), "--junit", str(junit)]
    args += ["--out", str(plan), "--header-tests", "tests/params/test_check.py"]
    assert C.main(["plan", "--root", str(root), *args]) == 0
    by_line = {p["line"]: p for p in C.read_jsonl(plan)}
    assert by_line[1]["scope"] == "header"
    assert by_line[1]["area_tests"] == ["tests/foo/test_f.py::test_f"]
    assert by_line[1]["suite_tests"] == ["tests/params/test_check.py::test_registry"]
    assert by_line[2]["suite_tests"] == []
