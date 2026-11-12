#!/usr/bin/env python3
"""The map of the MATLAB test suite: every MATLAB test with its fate in the port.

The map lives in tests/matlab_map/: suite.json (the oracle commit and the counts of the MATLAB
suite) and one <folder>.json per folder of the MATLAB tests/ (classes, their methods and the
parameterized cases of each, as matlab.unittest names them). Every class has a category, an
owner (the port task that accounts for it) and a flag `finished`; every method (or single
case) has a fate, a list of

    pytest:<path under python/>::<test>   ported, e.g. pytest:tests/crop/test_crop.py::test_crop3d
    oracle:<path under python/>::<case>   covered by an oracle case (a case id, '*' globs allowed)
    na:<reason>                           not applicable to the port, and why

An empty fate is pending. A class's `fate` applies to methods without one of their own; a
method's `case_fates` ({parameter string: fate}) to single cases.

    check_test_map.py check [--strict]   validate (CI): every MATLAB test is in the map once, the
                                         counts match suite.json, references exist; tests of a
                                         finished class (all with --strict) have a fate
    check_test_map.py report             counts by folder, category, owner and fate
    check_test_map.py docs               write docs/test_map/README.md from the map
    check_test_map.py update --suite DUMP --tests-dir DIR --commit C
                                         merge a new MATLAB suite (pkpy_dump_suite.m on a new
                                         oracle tree): new tests come in pending (new classes
                                         unassigned), removed ones are dropped and listed,
                                         fates of the others are kept
"""

from __future__ import annotations

import argparse
import ast
import collections
import fnmatch
import json
import re
import sys
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
MAP_DIR = PYTHON_ROOT / "tests" / "matlab_map"
DOCS = PYTHON_ROOT / "docs" / "test_map" / "README.md"

CATEGORIES = {
    "data": "data functions and wrappers (numeric and file level): ported in intent, numeric "
    "checks as oracle cases",
    "framework_fake_slurm": "job logic against the fake Slurm: ported nearly one to one",
    "framework_design": "framework tests tied to MATLAB's design (runners, pools, real cluster, "
    "isolation): ported to the Python design",
    "matlab_only": "MATLAB-only mechanisms (mex fall-backs and packaging, path shadowing, mcc): "
    "not applicable",
    "unclassified": "new class, not classified yet",
}
FATE_KINDS = ("pytest", "oracle", "na")
KNOWN_ISSUE = "KnownIssue"


# ------------------------------------------------------------------------------------------
# loading
# ------------------------------------------------------------------------------------------


def load_map(map_dir: Path = MAP_DIR) -> tuple[dict, dict[str, dict]]:
    suite = json.loads((map_dir / "suite.json").read_text())
    folders = {}
    for f in sorted(map_dir.glob("*.json")):
        if f.name == "suite.json":
            continue
        folders[f.stem] = json.loads(f.read_text())
    return suite, folders


def write_folder(map_dir: Path, folder: str, data: dict) -> None:
    (map_dir / f"{folder}.json").write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n")


def case_name(cls: str, method: str, params: str | None) -> str:
    return f"{cls}/{method}" + (f"({params})" if params else "")


def iter_cases(folders: dict[str, dict]):
    """(folder, class, method, params or None, class record, method record) of every case"""
    for folder, fd in folders.items():
        for cls, cr in fd["classes"].items():
            for meth, mr in cr["methods"].items():
                for p in mr["cases"] if mr.get("cases") else [None]:
                    yield folder, cls, meth, p, cr, mr


def case_fate(cr: dict, mr: dict, params: str | None) -> list[str]:
    if params is not None and params in mr.get("case_fates", {}):
        return mr["case_fates"][params]
    return mr.get("fate") or cr.get("fate") or []


def method_value(cr: dict, mr: dict, key: str):
    return mr.get(key, cr.get(key))


# ------------------------------------------------------------------------------------------
# check
# ------------------------------------------------------------------------------------------


def _test_names(path: Path, cache: dict) -> set[str]:
    """test functions of a pytest file: 'test_x' and 'TestY::test_x'"""
    if path not in cache:
        names: set[str] = set()
        try:
            tree = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            cache[path] = None
            return None
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                names.add(node.name)
            elif isinstance(node, ast.ClassDef):
                for sub in node.body:
                    if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        names.add(f"{node.name}::{sub.name}")
        cache[path] = names
    return cache[path]


def check_fate(fate: list, root: Path, cache: dict) -> list[str]:
    """problems with one fate list"""
    errs = []
    if not isinstance(fate, list):
        return [f"fate must be a list, not {type(fate).__name__}"]
    for f in fate:
        if not isinstance(f, str) or ":" not in f:
            errs.append(f"bad fate {f!r}")
            continue
        kind, ref = f.split(":", 1)
        if kind not in FATE_KINDS:
            errs.append(f"bad fate kind {kind!r} in {f!r}")
            continue
        if kind == "na":
            if len(ref.strip()) < 5:
                errs.append(f"{f!r}: give the reason")
            continue
        if "::" not in ref:
            errs.append(f"{f!r}: needs <path>::<name>")
            continue
        rel, name = ref.split("::", 1)
        path = root / rel
        if not path.is_file():
            errs.append(f"{f!r}: no file {rel}")
            continue
        if kind == "pytest":
            names = _test_names(path, cache)
            base = re.sub(r"\[.*\]$", "", name)
            if names is None:
                errs.append(f"{f!r}: {rel} does not parse")
            elif base not in names:
                errs.append(f"{f!r}: no test {base} in {rel}")
        else:
            text = cache.setdefault(("text", path), path.read_text())
            # the id up to the first '*' is in the file, or each of its dot-separated parts is
            # (ids made from a table, e.g. f"golden.{scenario}.{what}")
            literal = name.split("*", 1)[0]
            parts = [p for p in literal.split(".") if p]
            if not parts or (literal not in text and not all(p in text for p in parts)):
                errs.append(f"{f!r}: no case {name} in {rel}")
    return errs


def check(
    map_dir: Path = MAP_DIR, root: Path = PYTHON_ROOT, strict: bool = False
) -> tuple[list[str], dict]:
    """(problems, counts) of the map"""
    errs: list[str] = []
    suite, folders = load_map(map_dir)
    if sorted(folders) != sorted(suite.get("folders", {})):
        errs.append(
            f"folders {sorted(folders)} differ from suite.json {sorted(suite.get('folders', {}))}"
        )
    seen: set[str] = set()
    n_all = n_default = 0
    per_folder: collections.Counter = collections.Counter()
    cache: dict = {}
    for folder, fd in folders.items():
        if fd.get("folder") != folder:
            errs.append(f"{folder}.json: folder field {fd.get('folder')!r}")
        for cls, cr in fd.get("classes", {}).items():
            where = f"{folder}/{cls}"
            cat = cr.get("category")
            if cat not in CATEGORIES:
                errs.append(f"{where}: category {cat!r} not one of {sorted(CATEGORIES)}")
            elif cat == "unclassified":
                errs.append(f"{where}: unclassified")
            if not cr.get("owner") or cr.get("owner") == "unassigned":
                errs.append(f"{where}: no owner")
            if not isinstance(cr.get("finished"), bool):
                errs.append(f"{where}: finished must be true or false")
            errs += [f"{where}: {e}" for e in check_fate(cr.get("fate", []), root, cache)]
            if not cr.get("methods"):
                errs.append(f"{where}: no methods")
            for meth, mr in cr.get("methods", {}).items():
                mw = f"{where}/{meth}"
                if "category" in mr and mr["category"] not in CATEGORIES:
                    errs.append(f"{mw}: category {mr['category']!r}")
                errs += [f"{mw}: {e}" for e in check_fate(mr.get("fate", []), root, cache)]
                cases = mr.get("cases")
                if cases is not None and (not isinstance(cases, list) or not cases):
                    errs.append(f"{mw}: cases must be null or a non-empty list")
                    continue
                for p, fate in mr.get("case_fates", {}).items():
                    if p not in (cases or []):
                        errs.append(f"{mw}: case_fates names unknown case ({p})")
                    errs += [f"{mw}({p}): {e}" for e in check_fate(fate, root, cache)]
                ki = mr.get("known_issue", False)
                for p in cases or [None]:
                    name = case_name(cls, meth, p)
                    if name in seen:
                        errs.append(f"{name}: listed twice")
                    seen.add(name)
                    n_all += 1
                    per_folder[folder] += 1
                    is_ki = ki is True or (isinstance(ki, list) and p in ki)
                    n_default += not is_ki
                    if (strict or cr.get("finished")) and not case_fate(cr, mr, p):
                        errs.append(
                            f"{name}: no fate"
                            + ("" if strict else f" (class finished, owner {cr.get('owner')})")
                        )
    if n_all != suite.get("n_all"):
        errs.append(f"the map has {n_all} tests, suite.json {suite.get('n_all')}")
    if n_default != suite.get("n_default"):
        errs.append(
            f"the map has {n_default} tests without KnownIssue, suite.json {suite.get('n_default')}"
        )
    for folder, n in suite.get("folders", {}).items():
        if per_folder[folder] != n:
            errs.append(f"{folder}: {per_folder[folder]} tests in the map, suite.json {n}")
    return errs, {"n_all": n_all, "n_default": n_default}


# ------------------------------------------------------------------------------------------
# report and docs
# ------------------------------------------------------------------------------------------


def fate_kind(fate: list[str]) -> str:
    kinds = sorted({f.split(":", 1)[0] for f in fate})
    return "+".join(kinds) if kinds else "pending"


def tallies(folders: dict[str, dict]) -> dict[str, collections.Counter]:
    t: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    classes: dict[str, set] = collections.defaultdict(set)
    for folder, cls, _meth, p, cr, mr in iter_cases(folders):
        cat = method_value(cr, mr, "category")
        owner = method_value(cr, mr, "owner")
        ki = mr.get("known_issue", False)
        is_ki = ki is True or (isinstance(ki, list) and p in ki)
        kind = fate_kind(case_fate(cr, mr, p))
        for key, val in (("folder", folder), ("category", cat), ("owner", owner), ("fate", kind)):
            t[key][val] += 1
            if not is_ki:
                t[key + "_default"][val] += 1
            classes[f"{key}:{val}"].add(cls)
        t["folder_category"][(folder, cat)] += 1
        if not is_ki:
            t["folder_category_default"][(folder, cat)] += 1
    t["classes"] = collections.Counter({k: len(v) for k, v in classes.items()})
    return t


def _table(rows: list[list], head: list[str]) -> list[str]:
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def report_lines(suite: dict, folders: dict[str, dict]) -> list[str]:
    t = tallies(folders)
    cats = [c for c in CATEGORIES if t["category"][c]]
    n_classes = sum(len(f["classes"]) for f in folders.values())
    n_methods = sum(len(c["methods"]) for f in folders.values() for c in f["classes"].values())
    lines = [
        f"MATLAB suite of oracle {suite['oracle_commit']} ({suite['matlab']}): "
        f"{suite['n_all']} tests, {suite['n_default']} in a default run (without {KNOWN_ISSUE}), "
        f"{n_classes} classes, {n_methods} methods.",
        "",
    ]
    class_cat = collections.Counter(
        cr["category"] for fd in folders.values() for cr in fd["classes"].values()
    )
    lines += [
        "By category (tests: default run, all; classes of the category, with such tests):",
        "",
    ]
    lines += _table(
        [
            [
                c,
                t["category_default"][c],
                t["category"][c],
                class_cat[c],
                t["classes"][f"category:{c}"],
            ]
            for c in cats
        ],
        ["category", "default run", "all", "classes", "classes with such tests"],
    )
    lines += ["", "By MATLAB folder and category (default run):", ""]
    rows = []
    for folder in sorted(folders, key=lambda f: -t["folder"][f]):
        rows.append(
            [folder, t["folder_default"][folder]]
            + [t["folder_category_default"][(folder, c)] or "" for c in cats]
            + [t["folder"][folder]]
        )
    rows.append(
        ["total", suite["n_default"]] + [t["category_default"][c] for c in cats] + [suite["n_all"]]
    )
    lines += _table(rows, ["folder", "default run"] + cats + ["all"])
    lines += ["", "By owner (default run / all):", ""]
    lines += _table(
        [
            [o, t["owner_default"][o], t["owner"][o], t["classes"][f"owner:{o}"]]
            for o in sorted(t["owner"])
        ],
        ["owner", "default run", "all", "classes"],
    )
    lines += ["", "By fate (all):", ""]
    lines += _table([[k, v] for k, v in sorted(t["fate"].items())], ["fate", "tests"])
    return lines


def docs_text(suite: dict, folders: dict[str, dict]) -> str:
    lines = [
        "# The MATLAB tests in the port",
        "",
        "Generated by `python tools/check_test_map.py docs` from `tests/matlab_map/`; edit the "
        "map, not this file.",
        "",
        "Every MATLAB test (class, method and parameterized case) has a fate: ported as a pytest "
        "test, covered by an oracle case, or not applicable with a reason. The area task named as "
        "owner fills its classes and marks them finished; CI "
        "(`tests/matlab_map/test_matlab_map.py`) fails when a test of a finished class has no "
        "fate.",
        "",
        "Categories:",
        "",
    ]
    lines += [f"- `{k}`: {v}" for k, v in CATEGORIES.items() if k != "unclassified"]
    lines += ["", "## Counts", ""] + report_lines(suite, folders)
    lines += ["", "## Classes", ""]
    rows = []
    for folder in sorted(folders):
        for cls, cr in sorted(folders[folder]["classes"].items()):
            n = sum(len(m["cases"]) if m.get("cases") else 1 for m in cr["methods"].values())
            kinds = collections.Counter(
                fate_kind(case_fate(cr, mr, p))
                for _, c, _, p, cr2, mr in iter_cases({folder: {"classes": {cls: cr}}})
            )
            rows.append(
                [
                    f"{folder}/{cls}",
                    cr["category"],
                    cr["owner"],
                    len(cr["methods"]),
                    n,
                    "yes" if cr["finished"] else "no",
                    ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())),
                ]
            )
    lines += _table(rows, ["class", "category", "owner", "methods", "tests", "finished", "fates"])
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------------------------------
# update from a MATLAB suite dump (pkpy_dump_suite.m)
# ------------------------------------------------------------------------------------------

_FUNC = re.compile(r"^\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?(\w+)\s*\(")


def _comment_block(lines: list[str], start: int) -> str:
    out = []
    for line in lines[start:]:
        s = line.strip()
        if not s.startswith("%"):
            break
        s = s.lstrip("%").strip()
        if not s and out:
            break
        if s:
            out.append(s)
    return " ".join(out)


def _humanize(name: str) -> str:
    return re.sub(r"(?<=[a-z0-9])(?=[A-Z])", " ", name).lower()


def parse_test_file(path: Path) -> tuple[str, dict[str, dict]]:
    """the class comment and, per function, its leading comment and body text"""
    lines = path.read_text(errors="replace").splitlines()
    cls_pins = _comment_block(lines, 1)
    funcs: dict[str, dict] = {}
    starts = [(i, m.group(1)) for i, line in enumerate(lines) if (m := _FUNC.match(line))]
    for j, (i, name) in enumerate(starts):
        end = starts[j + 1][0] if j + 1 < len(starts) else len(lines)
        funcs.setdefault(
            name,
            {"pins": _comment_block(lines, i + 1), "body": "\n".join(lines[i:end]), "line": i + 1},
        )
    return cls_pins, funcs


def _pins(text: str, limit: int = 400) -> str:
    text = re.sub(r"\s+", " ", text).strip()
    return text if len(text) <= limit else text[: limit - 3].rstrip() + "..."


def _as_list(v) -> list:
    return v if isinstance(v, list) else [] if v is None else [v]


def normalize_dump(dump: dict) -> dict[tuple[str, str], dict]:
    """(class, method) -> {folder, tags, cases: [params or None]} from pkpy_dump_suite's JSON"""
    out: dict[tuple[str, str], dict] = {}
    for t in _as_list(dump["tests"]):
        name, cls, meth = t["name"], t["class"], t["procedure"]
        m = re.fullmatch(re.escape(cls) + r"/" + re.escape(meth) + r"(?:\((.*)\))?", name)
        if not m:
            raise ValueError(f"unexpected test name {name}")
        rec = out.setdefault(
            (cls, meth), {"folder": t["folder"], "tags": set(), "cases": [], "ki": []}
        )
        tags = set(_as_list(t["tags"]))
        rec["tags"] |= tags
        rec["cases"].append(m.group(1))
        if KNOWN_ISSUE in tags:
            rec["ki"].append(m.group(1))
    return out


def update(
    dump_path: Path,
    tests_dir: Path,
    commit: str,
    map_dir: Path = MAP_DIR,
    rules: dict | None = None,
) -> list[str]:
    """merge a MATLAB suite dump into the map; returns what changed"""
    dump = json.loads(dump_path.read_text())
    tests = normalize_dump(dump)
    try:
        _, folders = load_map(map_dir)
    except FileNotFoundError:
        folders = {}
    rules = rules or {}
    notes: list[str] = []
    old_index = {
        (c, m): (f, cr, mr)
        for f, fd in folders.items()
        for c, cr in fd["classes"].items()
        for m, mr in cr["methods"].items()
    }
    new_folders: dict[str, dict] = {}
    parsed: dict[str, tuple] = {}
    for (cls, meth), rec in sorted(tests.items(), key=lambda kv: (kv[1]["folder"], kv[0])):
        folder = rec["folder"]
        fd = new_folders.setdefault(folder, {"folder": folder, "classes": {}})
        if cls not in parsed:
            p = tests_dir / folder / f"{cls}.m"
            parsed[cls] = parse_test_file(p) if p.is_file() else ("", {})
        cls_pins, funcs = parsed[cls]
        old_cls = next(
            (folders[f]["classes"][cls] for f in folders if cls in folders[f]["classes"]), None
        )
        if cls not in fd["classes"]:
            rule = rules.get(cls, {})
            base = {k: v for k, v in (old_cls or {}).items() if k != "methods"}
            cr = {
                "file": f"tests/{folder}/{cls}.m",
                "pins": _pins(cls_pins) or base.get("pins", ""),
                "category": base.get("category", rule.get("category", "unclassified")),
                "owner": base.get("owner", rule.get("owner", "unassigned")),
                "phase": base.get("phase", rule.get("phase", 1)),
                "finished": base.get("finished", rule.get("finished", False)),
                "fate": base.get("fate", rule.get("fate", [])),
                "methods": {},
            }
            for k, v in base.items():
                cr.setdefault(k, v)
            if old_cls is None:
                notes.append(
                    f"new class {folder}/{cls}" + ("" if rule else " (unclassified, no owner)")
                )
            fd["classes"][cls] = cr
        cr = fd["classes"][cls]
        cases = rec["cases"]
        old = old_index.get((cls, meth))
        fn = funcs.get(meth)
        # the method's comment; without the source the map's text; else the method name
        pins = fn["pins"] if fn else (old[2].get("pins", "") if old else "")
        mr = {
            "pins": _pins(pins or _humanize(meth)),
            "tags": sorted(rec["tags"]),
            "cases": None if cases == [None] else cases,
            "fate": [],
        }
        if rec["ki"]:
            mr["known_issue"] = True if len(rec["ki"]) == len(cases) else rec["ki"]
        mrule = _method_rule(rules.get(cls, {}), meth, fn["body"] if fn else "")
        mr.update(mrule)
        if old is not None:
            _, _, omr = old
            for k, v in omr.items():
                if k not in ("pins", "tags", "cases", "known_issue"):
                    mr[k] = v
            if mr["cases"] != omr.get("cases"):
                gone = set(omr.get("cases") or []) - set(mr["cases"] or [])
                n_old = len(omr.get("cases") or [None])
                notes.append(f"{cls}/{meth}: cases changed ({n_old} -> {len(cases)})")
                mr["case_fates"] = {
                    p: f for p, f in mr.get("case_fates", {}).items() if p not in gone
                }
                if not mr["case_fates"]:
                    mr.pop("case_fates")
        else:
            if old_cls is not None:
                notes.append(f"new method {cls}/{meth}")
                if cr.get("finished"):
                    cr["finished"] = False
                    notes.append(f"{cls}: no longer finished (new method {meth})")
        cr["methods"][meth] = mr
    for (cls, meth), (f, _, _) in sorted(old_index.items()):
        if (cls, meth) not in tests:
            notes.append(f"removed {f}/{cls}/{meth}")
    counts = collections.Counter(rec["folder"] for rec in tests.values() for _ in rec["cases"])
    n_all = sum(counts.values())
    n_ki = sum(len(rec["ki"]) for rec in tests.values())
    if n_all != dump.get("n_all", n_all) or n_all - n_ki != dump.get("n_default", n_all - n_ki):
        raise ValueError("the dump's counts disagree with its tests")
    suite = {
        "oracle_commit": commit,
        "matlab": dump.get("matlab", ""),
        "computer": dump.get("computer", ""),
        "n_all": n_all,
        "n_default": n_all - n_ki,
        "folders": dict(sorted(counts.items())),
    }
    map_dir.mkdir(parents=True, exist_ok=True)
    for f in map_dir.glob("*.json"):
        if f.stem not in new_folders and f.name != "suite.json":
            f.unlink()
    (map_dir / "suite.json").write_text(json.dumps(suite, indent=1) + "\n")
    for folder, fd in new_folders.items():
        write_folder(map_dir, folder, fd)
    return notes


def _method_rule(rule: dict, meth: str, body: str) -> dict:
    """method fields from a class rule: 'methods' [(glob, {fields})], first match wins;
    'fake_slurm_category': the category of methods that use the fake Slurm"""
    out: dict = {}
    for pat, fields in rule.get("methods", []):
        if fnmatch.fnmatchcase(meth, pat):
            out.update(fields)
            break
    if rule.get("fake_slurm_category") and "FakeSlurmFixture" in body and "category" not in out:
        out["category"] = rule["fake_slurm_category"]
    return out


# ------------------------------------------------------------------------------------------
# command line
# ------------------------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--map", type=Path, default=MAP_DIR, help="map folder (default tests/matlab_map)"
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("--strict", action="store_true", help="every test needs a fate")
    sub.add_parser("report")
    d = sub.add_parser("docs")
    d.add_argument("--out", type=Path, default=DOCS)
    u = sub.add_parser("update")
    u.add_argument("--suite", type=Path, required=True, help="JSON written by pkpy_dump_suite.m")
    u.add_argument("--tests-dir", type=Path, required=True, help="tests/ of the oracle tree")
    u.add_argument("--commit", required=True, help="the oracle tree's commit")
    u.add_argument(
        "--rules", type=Path, help="JSON of class rules (category, owner, phase, fate, methods)"
    )
    ns = ap.parse_args(argv)
    if ns.cmd == "check":
        errs, counts = check(ns.map, strict=ns.strict)
        for e in errs:
            print(e)
        print(
            f"{counts['n_all']} tests ({counts['n_default']} without {KNOWN_ISSUE}): "
            + ("ok" if not errs else f"{len(errs)} problems")
        )
        return 1 if errs else 0
    suite_folders = None
    if ns.cmd == "update":
        rules = json.loads(ns.rules.read_text()) if ns.rules else None
        for line in update(ns.suite, ns.tests_dir, ns.commit, ns.map, rules):
            print(line)
        return 0
    suite_folders = load_map(ns.map)
    if ns.cmd == "report":
        print("\n".join(report_lines(*suite_folders)))
    elif ns.cmd == "docs":
        ns.out.parent.mkdir(parents=True, exist_ok=True)
        ns.out.write_text(docs_text(*suite_folders))
        print(f"wrote {ns.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
