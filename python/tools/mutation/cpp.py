#!/usr/bin/env python3
"""Single-line mutants of the C++ extension modules (``python/cpp/<group>/``), built and run
against the Python tests of the modules that load each group (docs/hardening.md).

    cpp.py generate --root PY --out mutants.jsonl [--seed S] [--per-group auto|N]
                    [--group G]... [--exclude GLOB]...
    cpp.py plan     --root PY --mutants mutants.jsonl --coverage DATA --junit XML...
                    --out plan.jsonl [--budget S]
    cpp.py run      --root PY --plan plan.jsonl --out results.jsonl --work DIR [--jobs J]
                    [--shard I/N] [--cmake-arg -DCMAKE_CXX_COMPILER=g++]...

The results and plans have the fields campaign.py's ``report`` reads (area ``cpp/<group>``).

generate  sites are found by tokens, not a parser: comments, strings, characters and
          preprocessor lines are left alone; binary operators only with blanks on both sides
          (``a < b``, not ``std::vector<int>``). Operators: cmp (< <=, > >=, == !=), arith
          (+ -, * /, and the compound forms), num (a decimal integer n -> n + 1), bool (&& ||,
          true false). ``auto`` takes ceil(lines with sites / 120) per group, 2 to 6. A mutant
          that doesn't compile is reported as such and left out of the scores.
plan      the tests of a group are those that run any line of the Python modules that load
          it (``load("<group>")``, from the coverage contexts), fastest first, up to --budget
          seconds; tests that failed or were skipped there are left out.
run       one copy of python/ and one CMake build per job (--work, node-local), then per
          mutant: the file mutated, an incremental build, pytest -x on its tests with the
          job's build as PETAKIT5D_EXT_DIR, the file restored.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import campaign as C  # noqa: E402
import mutants as M  # noqa: E402

EXTENSIONS = (".cpp", ".hpp", ".h", ".inc", ".cc")
_CMP = {"<": "<=", "<=": "<", ">": ">=", ">=": ">", "==": "!=", "!=": "=="}
_ARITH = {"+": "-", "-": "+", "*": "/", "/": "*", "+=": "-=", "-=": "+=", "*=": "/=", "/=": "*="}
_BOOL = {"&&": "||", "||": "&&", "true": "false", "false": "true"}
_OP_RE = re.compile(
    r"(?<=\s)(<=|>=|==|!=|\+=|-=|\*=|/=|&&|\|\||<|>|\+|-|\*|/)(?=\s)"
    r"|(?<![\w.])(\d+)(?![\w.'])"
    r"|\b(true|false)\b"
)


def _code_lines(text: str) -> list[str]:
    """each line with comments, string and char literals blanked and preprocessor lines
    emptied (same length, so columns stay)"""
    out = []
    in_block = False
    for line in text.split("\n"):
        chars = list(line)
        i, n = 0, len(chars)
        if not in_block and line.lstrip().startswith("#"):
            out.append(" " * n)
            continue
        while i < n:
            if in_block:
                if line.startswith("*/", i):
                    chars[i] = chars[i + 1] = " "
                    in_block = False
                    i += 2
                    continue
                chars[i] = " "
                i += 1
                continue
            if line.startswith("//", i):
                for k in range(i, n):
                    chars[k] = " "
                break
            if line.startswith("/*", i):
                in_block = True
                chars[i] = chars[i + 1] = " "
                i += 2
                continue
            if chars[i] in "\"'":
                q = chars[i]
                k = i + 1
                while k < n and chars[k] != q:
                    k += 2 if chars[k] == "\\" else 1
                for j in range(i, min(k + 1, n)):
                    chars[j] = " "
                i = k + 1
                continue
            i += 1
        out.append("".join(chars))
    return out


def find_cpp_sites(text: str) -> list[M.Site]:
    """the candidate sites of a C++ source (byte columns of the line's UTF-8 text)"""
    sites = []
    raw = text.split("\n")
    for ln, (code, line) in enumerate(zip(_code_lines(text), raw, strict=True), 1):
        for m in _OP_RE.finditer(code):
            tok = m.group(0)
            if m.group(2) is not None:
                op, new = "num", str(int(tok) + 1)
            elif tok in _CMP:
                op, new = "cmp", _CMP[tok]
            elif tok in _ARITH:
                op, new = "arith", _ARITH[tok]
            else:
                op, new = "bool", _BOOL[tok]
            col = len(line[: m.start()].encode("utf-8"))
            sites.append(M.Site(ln, col, col + len(tok.encode("utf-8")), op, tok, new))
    return sites


def _sources(root: Path) -> list[str]:
    out = []
    for p in sorted((root / "cpp").rglob("*")):
        rel = p.relative_to(root).as_posix()
        # (not the vendored libraries, the GPL add-on or the headers every group shares)
        skip = "third_party" in p.parts or "/exact/" in rel or p.parts[-2:-1] == ("common",)
        if p.suffix in EXTENSIONS and not skip and "/common/" not in rel:
            out.append(rel)
    return out


def group_of(rel: str) -> str:
    return Path(rel).parts[1]


def cmd_generate(a: argparse.Namespace) -> int:
    import random

    root = Path(a.root)
    by_group: dict[str, dict[str, list[M.Site]]] = defaultdict(dict)
    texts = {}
    for rel in _sources(root):
        g = group_of(rel)
        if (a.group and g not in a.group) or (a.exclude and C._matches(rel, a.exclude)):
            continue
        text = (root / rel).read_text(encoding="utf-8")
        sites = find_cpp_sites(text)
        if sites:
            texts[rel] = text
            by_group[g][rel] = sites
    out = []
    for g in sorted(by_group):
        sites = by_group[g]
        lines = sorted({(f, s.line) for f, ss in sites.items() for s in ss})
        n = (
            int(a.per_group)
            if a.per_group != "auto"
            else min(6, max(2, math.ceil(len(lines) / 120)))
        )
        rng = random.Random(f"{a.seed}:{g}")
        for k, (f, ln) in enumerate(rng.sample(lines, min(n, len(lines))), 1):
            site = rng.choice([s for s in sites[f] if s.line == ln])
            out.append(M.make_mutant(f"cpp/{g}-{k:03d}", f"cpp/{g}", f, texts[f], site))
        print(f"cpp/{g:10s} {len(lines):5d} lines with sites, {min(n, len(lines))} mutants")
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        for m in out:
            fh.write(json.dumps(m.record()) + "\n")
    print(f"{len(out)} mutants -> {a.out}")
    return 0


_LOAD_RE = re.compile(r"""load\(\s*["'](\w+)["']""")
# loaders the pattern doesn't see (a module name through a variable)
_LOADERS = {"deskew": ["src/petakit5d/kernels/deskew.py"]}


def loaders(root: Path) -> dict[str, list[str]]:
    """group -> the package modules (python/ relative) that load it"""
    out: dict[str, list[str]] = defaultdict(list)
    for rel in C.modules(root):
        text = (root / rel).read_text(encoding="utf-8")
        for g in sorted(set(_LOAD_RE.findall(text))):
            out[g].append(rel)
    for g, rels in _LOADERS.items():
        out[g] = sorted(set(out[g]) | {r for r in rels if (root / r).is_file()})
    return out


def cmd_plan(a: argparse.Namespace) -> int:
    from coverage import CoverageData

    root = Path(a.root)
    data = CoverageData(a.coverage)
    data.read()
    by_rel = {}
    for f in data.measured_files():
        i = f.replace(os.sep, "/").find("/" + C.PACKAGE + "/")
        if i >= 0:
            by_rel[f.replace(os.sep, "/")[i + 1 :]] = f
    seconds, bad = C.junit_tests(a.junit, root)
    load = loaders(root)
    tests_of: dict[str, set[str]] = {}
    out = []
    for rec in C.read_jsonl(a.mutants):
        g = group_of(rec["file"])
        if g not in tests_of:
            ts = set()
            for rel in load.get(g, []):
                f = by_rel.get(rel)
                if f:
                    for cs in data.contexts_by_lineno(f).values():
                        ts |= {C._context_test(c) for c in cs if c}
            tests_of[g] = ts - bad
        ts = sorted(tests_of[g], key=lambda t: (seconds.get(t, 1.0), t))
        chosen, total = [], 0.0
        for t in ts:
            s = seconds.get(t, 1.0)
            if chosen and total + s > a.budget:
                break
            chosen.append(t)
            total += s
        p = dict(rec, covered=bool(chosen), area_tests=chosen, area_seconds=total)
        p.update(area_all=len(ts), suite_tests=[], suite_seconds=0.0, loaders=load.get(g, []))
        out.append(p)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        for p in out:
            fh.write(json.dumps(p) + "\n")
    print(f"{len(out)} mutants -> {a.out}")
    return 0


def _build(build: Path, jobs: int, log: Path) -> bool:
    with open(log, "w", encoding="utf-8") as fh:
        r = subprocess.run(["cmake", "--build", str(build), "-j", str(jobs)], stdout=fh,
                           stderr=subprocess.STDOUT)  # fmt: skip
    return r.returncode == 0


def cmd_run(a: argparse.Namespace) -> int:
    plan = C.read_jsonl(a.plan)
    if a.shard:
        i, n = (int(x) for x in a.shard.split("/"))
        plan = [p for k, p in enumerate(plan) if k % n == i - 1]
    done = {r["id"] for r in C.read_jsonl(a.out)} if os.path.exists(a.out) else set()
    plan = [p for p in plan if p["id"] not in done]
    work = Path(a.work).resolve()
    tree, build = work / "tree", work / "build"
    if not (build / "CMakeCache.txt").exists():
        skip = shutil.ignore_patterns("__pycache__", "build", ".pytest_cache")
        shutil.copytree(Path(a.root), tree, dirs_exist_ok=True, ignore=skip)
        cfg = ["cmake", "-S", str(tree), "-B", str(build), "-G", "Ninja",
               "-DCMAKE_BUILD_TYPE=Release", f"-DPython_EXECUTABLE={sys.executable}",
               *a.cmake_arg]  # fmt: skip
        if subprocess.run(cfg, capture_output=True).returncode != 0:
            raise SystemExit("cpp.py: configure failed")
    if not _build(build, a.jobs, work / "build_base.log"):
        raise SystemExit("cpp.py: the base build failed")
    a.root, a.threads, a.mem_gb, a.tests, a.given_budget = str(tree), a.jobs, 0, [], 0
    runner = C.Runner(a)
    os.environ["PETAKIT5D_EXT_DIR"] = str(build / "petakit5d" / "_ext")
    print(f"{len(plan)} mutants to run, base built", flush=True)
    for p in plan:
        t0 = time.monotonic()
        m = M.Mutant.from_record(p)
        path = tree / m.file
        text = path.read_text(encoding="utf-8")
        res = {"id": m.id, "area": m.area, "file": m.file, "line": m.line, "op": m.op}
        try:
            path.write_text(M.apply_mutant(text, m), encoding="utf-8", newline="")
            blog = work / f"build_{m.id.replace('/', '_')}.log"
            if not _build(build, a.jobs, blog):
                res["status"] = "uncompiled"
                lines = blog.read_text(encoding="utf-8", errors="replace").splitlines()
                res["error"] = "\n".join(x for x in lines if "error" in x.lower())[-2000:]
            elif not p["area_tests"]:
                res["status"] = "no_tests"
            else:
                limit = min(a.max_seconds, 120 + 4 * float(p["area_seconds"]))
                tmp = work / "m" / m.id.replace("/", "_")
                r = runner.pytest(tree / "src", p["area_tests"], tmp, limit)
                r["n"] = len(p["area_tests"])
                res["phases"] = {"area": r}
                rc = r["rc"]
                if r["timeout"]:
                    res["status"], res["by"] = "timeout", "area"
                elif rc == 0:
                    res["status"] = "survived"
                elif rc in (1, 2, 3) or (rc is not None and rc < 0):
                    res["status"], res["by"] = "killed", "area"
                else:
                    res["status"], res["error"] = "error", f"pytest exit {rc}"
        finally:
            path.write_text(text, encoding="utf-8", newline="")
        res["wall"] = round(time.monotonic() - t0, 1)
        with open(a.out, "a", encoding="utf-8", newline="\n") as fh:
            fh.write(json.dumps(res) + "\n")
        print(f"{m.id:18s} {res['status']:12s} {res['wall']:7.1f} s", flush=True)
    _build(build, a.jobs, work / "build_restored.log")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--root", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("--seed", default="0")
    g.add_argument("--per-group", default="auto")
    g.add_argument("--group", action="append")
    g.add_argument("--exclude", action="append")
    p = sub.add_parser("plan")
    p.add_argument("--root", required=True)
    p.add_argument("--mutants", required=True)
    p.add_argument("--coverage", required=True)
    p.add_argument("--junit", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--budget", type=float, default=900)
    r = sub.add_parser("run")
    r.add_argument("--root", required=True)
    r.add_argument("--plan", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--work", required=True)
    r.add_argument("--jobs", type=int, default=4)
    r.add_argument("--shard", default="")
    r.add_argument("--max-seconds", type=float, default=3600)
    r.add_argument("--cmake-arg", action="append", default=[])
    a = ap.parse_args(argv)
    return {"generate": cmd_generate, "plan": cmd_plan, "run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
