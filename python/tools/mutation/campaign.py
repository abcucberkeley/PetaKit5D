#!/usr/bin/env python3
"""A mutation campaign over the port: single-line mutants, each in a copy of the source tree,
each run against the tests that execute its line (docs/hardening.md).

    campaign.py generate --root PY --out mutants.jsonl [--seed S] [--per-area auto|N]
                         [--area A]... [--only GLOB]... [--exclude GLOB]... [--ops OP,...]
    campaign.py plan     --root PY --mutants mutants.jsonl --coverage DATA --junit XML...
                         --out plan.jsonl [--area-budget S] [--suite-budget S]
    campaign.py run      --root PY --plan plan.jsonl --out results.jsonl [--work DIR]
                         [--parallel K] [--threads T] [--shard I/N] [--ids ID,...]
                         [--tests NODEID]... [--recheck MAP.json] [--no-suite] [--no-verify]
                         [--mem-gb G]
    campaign.py report   --plan plan.jsonl --results RESULTS... [--rechecks RESULTS...]
                         [--dispositions D.json] [--out report.md] [--json summary.json]
    campaign.py show     --root PY --mutants mutants.jsonl ID...
    campaign.py tree     --root PY --mutants mutants.jsonl ID --dest DIR

PY is a python/ folder of the repo (src/petakit5d, tests/). Areas are the subpackages of
petakit5d (plus "package" for the top-level modules); an area's tests are tests/<area>/.

generate  samples mutants per area (mutants.py has the operators and the sampling): ``auto``
          is ceil(lines with sites / 100), at least --min (12) and at most --max (50). The
          same --seed gives the same mutants for the same sources.
plan      selects each mutant's tests from a coverage run with per-test contexts (pytest-cov's
          --cov-context=test; subprocesses' data labelled ``<nodeid>|sub``) and the run's JUnit
          files (durations; tests that failed or were skipped there are left out): the area
          phase runs the tests of tests/<area>/ that execute the line, fastest first, up to
          --area-budget seconds; the suite phase the other tests that execute it, up to
          --suite-budget. A line of a function's signature or decorators (it runs at import
          time and takes effect in calls) gets the tests that run the function's body; a line
          at module or class level gets the whole tests/<area>/ and, in the suite phase, the
          tests that run any line of its module; a line no test runs is "not covered" and is
          not run.
run       makes each mutant's tree (hard links of the base src/ with the one module rewritten,
          on --work, a node-local folder), runs pytest -x on the area tests and, if the
          mutant survives, on the suite tests, from PY with PYTHONPATH pointing at the
          mutant's src/; a test that kills a mutant is run once on the base tree too (a test
          that fails there doesn't count). Results are appended as they come (a rerun skips
          the ids already in --out). --tests runs the given tests instead of the plan's;
          --recheck FILE ({id: [node ids]}) runs each listed mutant against its own tests
          (the new tests written for survivors).
report    per area: mutants, killed by the area's tests, killed only by other tests, the score
          of the suite, survived, not covered, killed by new tests (--rechecks: the results of
          run --recheck), given a reason (--dispositions: {id: {"kind": "equivalent"|
          "unreachable"|"diagnostic"|"gpu"|"owner"|..., "why": text}}), still open, and the
          score with the new tests; then every mutant the suite did not kill, with its
          disposition.

Statuses: killed (a test failed, or collection failed), timeout (counted as killed),
survived, not_covered (no test runs the line), no_tests (the line runs only in tests left out:
failed or skipped in the coverage run), base_fails (the killing test fails on the base tree
too: rerun or look), error (the tool could not run it).
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import mutants as M  # noqa: E402

PACKAGE = "src/petakit5d"
KILLED = ("killed", "timeout")
# dispositions that are no bug a test could show (left out of the last score)
EQUIVALENT_KINDS = ("equivalent", "unreachable", "diagnostic")


# --- areas -----------------------------------------------------------------------------------


def area_of(rel: str) -> str:
    """the area of a module path relative to python/ (src/petakit5d/<area>/...)"""
    parts = Path(rel).parts
    if len(parts) > 3 and parts[2] != "_ext":
        return parts[2]
    return "package"


def area_tests_prefix(area: str) -> tuple[str, ...]:
    if area == "package":
        return ("tests/test_",)
    return (f"tests/{area}/",)


def modules(root: Path) -> list[str]:
    return sorted(p.relative_to(root).as_posix() for p in (root / PACKAGE).rglob("*.py"))


def _matches(rel: str, globs: list[str]) -> bool:
    return any(fnmatch.fnmatch(rel, g) for g in globs)


# --- generate --------------------------------------------------------------------------------


def cmd_generate(a: argparse.Namespace) -> int:
    root = Path(a.root)
    ops = a.ops.split(",") if a.ops else M.OPERATORS
    by_area: dict[str, dict[str, list[M.Site]]] = defaultdict(dict)
    texts: dict[str, str] = {}
    for rel in modules(root):
        if a.only and not _matches(rel, a.only):
            continue
        if a.exclude and _matches(rel, a.exclude):
            continue
        area = area_of(rel)
        if a.area and area not in a.area:
            continue
        text = (root / rel).read_text(encoding="utf-8")
        sites = M.find_sites(text, ops)
        if sites:
            texts[rel] = text
            by_area[area][rel] = sites
    out = []
    for area in sorted(by_area):
        sites = by_area[area]
        n_lines = len({(f, s.line) for f, ss in sites.items() for s in ss})
        if a.per_area == "auto":
            n = min(a.max, max(a.min, math.ceil(n_lines / 100)))
        else:
            n = int(a.per_area)
        picked = list(M.sample_sites(sites, texts, n, f"{a.seed}:{area}"))
        for k, (rel, site) in enumerate(picked, 1):
            out.append(M.make_mutant(f"{a.prefix}{area}-{k:03d}", area, rel, texts[rel], site))
        print(f"{area:14s} {n_lines:6d} lines with sites, {len(picked):3d} mutants")
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        for m in out:
            fh.write(json.dumps(m.record()) + "\n")
    print(f"{len(out)} mutants -> {a.out}")
    return 0


def read_jsonl(path: str | Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(x) for x in fh if x.strip()]


# --- plan ------------------------------------------------------------------------------------


def junit_tests(paths: list[str], root: Path) -> tuple[dict[str, float], set[str]]:
    """(seconds per node id, node ids that failed, errored or were skipped) of JUnit files"""
    seconds: dict[str, float] = {}
    bad: set[str] = set()
    for p in paths:
        for tc in ET.parse(p).getroot().iter("testcase"):
            cls, name = tc.get("classname", ""), tc.get("name", "")
            nodeid = junit_nodeid(cls, name, root)
            seconds[nodeid] = seconds.get(nodeid, 0.0) + float(tc.get("time") or 0)
            if any(ch.tag in ("failure", "error", "skipped") for ch in tc):
                bad.add(nodeid)
    return seconds, bad


def junit_nodeid(classname: str, name: str, root: Path) -> str:
    """pytest's node id from JUnit's classname and name (``tests.a.test_b.TestC`` and
    ``test_d[x]``): the module is the longest dotted prefix that is a file under root, the
    classes follow it"""
    parts = classname.split(".")
    for k in range(len(parts), 0, -1):
        mod = "/".join(parts[:k]) + ".py"
        if (root / mod).is_file():
            return "::".join([mod, *parts[k:], name])
    return f"{classname}::{name}"


def _context_test(ctx: str) -> str:
    return ctx.rsplit("|", 1)[0] if "|" in ctx else ctx


def line_scopes(text: str) -> dict[int, tuple[str, tuple[int, int]]]:
    """where each line of a module runs: ("body", ...) inside a function's body; ("header",
    (first, last)) in a function's decorators or signature, which run at import time and take
    effect in calls: (first, last) are the lines of its body; lines of neither are ("module",
    ...): module and class level, run at import time"""
    import ast

    out: dict[int, tuple[str, tuple[int, int]]] = {}
    tree = ast.parse(text)
    funcs = [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    # (outer functions first, so that a nested function's lines are its own)
    for f in sorted(funcs, key=lambda n: (n.lineno, -n.end_lineno)):
        start = min([f.lineno, *(d.lineno for d in f.decorator_list)])
        b0, b1 = f.body[0].lineno, f.end_lineno
        for ln in range(start, b0):
            out[ln] = ("header", (b0, b1))
        for ln in range(b0, b1 + 1):
            out[ln] = ("body", (b0, b1))
    return out


def cmd_plan(a: argparse.Namespace) -> int:
    from coverage import CoverageData
    from coverage.parser import PythonParser

    root = Path(a.root)
    data = CoverageData(a.coverage)
    data.read()
    measured = {f.replace(os.sep, "/"): f for f in data.measured_files()}
    by_rel = {}
    for f_norm, f in measured.items():
        i = f_norm.find("/" + PACKAGE + "/")
        if i >= 0:
            by_rel[f_norm[i + 1 :]] = f
    seconds, bad = junit_tests(a.junit, root)
    header_files = set(a.header_tests or [])
    folder_seconds: dict[str, float] = defaultdict(float)
    for t, s in seconds.items():
        folder_seconds["/".join(t.split("/")[:2]) + "/"] += s
    contexts_cache: dict[str, dict[int, list[str]]] = {}
    parsers: dict[str, PythonParser] = {}
    scopes: dict[str, dict[int, tuple[str, tuple[int, int]]]] = {}
    out = []
    for rec in read_jsonl(a.mutants):
        m = M.Mutant.from_record(rec)
        rel = m.file
        prefix = area_tests_prefix(m.area)
        p = dict(rec)
        f = by_rel.get(rel)
        if f is None:
            ctxs: list[str] = []
        else:
            if f not in contexts_cache:
                contexts_cache[f] = data.contexts_by_lineno(f)
                text = (root / rel).read_text(encoding="utf-8")
                parsers[f] = PythonParser(text=text, filename=f)
                parsers[f].parse_source()
                scopes[f] = line_scopes(text)
            first = parsers[f].first_line(m.line)
            ctxs = contexts_cache[f].get(first, []) or contexts_cache[f].get(m.line, [])
        scope, (b0, b1) = scopes.get(f, {}).get(m.line, ("module", (0, 0))) if f else ("", (0, 0))
        p["scope"] = scope
        if scope == "header" and ctxs:
            # (defaults, decorators: they take effect where the function is called)
            body = contexts_cache[f]
            ctxs = sorted({c for ln in range(b0, b1 + 1) for c in body.get(ln, [])} | {""})
        elif scope == "module" and ctxs:
            # (constants, tables, class bodies: they take effect wherever the module is used)
            ctxs = sorted({c for cs in contexts_cache[f].values() for c in cs} | {""})
        tests = {_context_test(c) for c in ctxs if c}
        if scope in ("header", "module") and ctxs:
            # (tests that read signatures and tables without calling them: a registry check)
            tests |= {t for t in seconds if t.split("::")[0] in header_files}
        tests = sorted(t for t in tests if t not in bad)
        import_time = "" in ctxs
        area_t = [t for t in tests if t.startswith(prefix)]
        suite_t = [t for t in tests if not t.startswith(prefix)]

        def budget(ts: list[str], limit: float) -> tuple[list[str], float, int]:
            ts = sorted(ts, key=lambda t: (seconds.get(t, 1.0), t))
            chosen, total = [], 0.0
            for t in ts:
                s = seconds.get(t, 1.0)
                if chosen and total + s > limit:
                    break
                chosen.append(t)
                total += s
            return chosen, total, len(ts)

        p["covered"] = bool(ctxs)
        p["import_time"] = import_time
        p["area_tests"], p["area_seconds"], p["area_all"] = budget(area_t, a.area_budget)
        p["suite_tests"], p["suite_seconds"], p["suite_all"] = budget(suite_t, a.suite_budget)
        p["area_folder"] = False
        if (not area_t or scope == "module") and import_time and m.area != "package":
            folder = prefix[0]
            if (root / folder).is_dir():
                p["area_tests"] = [folder.rstrip("/")]
                p["area_seconds"] = folder_seconds.get(folder, 600.0)
                p["area_folder"] = True
        out.append(p)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        for p in out:
            fh.write(json.dumps(p) + "\n")
    n_cov = sum(p["covered"] for p in out)
    print(f"{len(out)} mutants, {n_cov} on covered lines -> {a.out}")
    return 0


# --- run -------------------------------------------------------------------------------------


class Runner:
    def __init__(self, a: argparse.Namespace) -> None:
        self.a = a
        self.root = Path(a.root).resolve()
        self.work = Path(a.work).resolve()
        self.lock = threading.Lock()
        self.base_src = self.work / "base_src"
        self.base_results: dict[str, bool] = {}
        self.base_lock = threading.Lock()
        self.recheck: dict[str, list[str]] = {}
        if getattr(a, "recheck", ""):
            self.recheck = json.loads(Path(a.recheck).read_text(encoding="utf-8"))

    def prepare(self) -> None:
        self.work.mkdir(parents=True, exist_ok=True)
        if self.base_src.exists():
            shutil.rmtree(self.base_src)
        shutil.copytree(
            self.root / "src", self.base_src, ignore=shutil.ignore_patterns("__pycache__")
        )

    def env(self, src: Path) -> dict[str, str]:
        env = dict(os.environ)
        env["PYTHONPATH"] = str(src)
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        t = str(self.a.threads)
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            env[v] = t
        env.setdefault("CUDA_VISIBLE_DEVICES", "")
        return env

    def pytest(self, src: Path, tests: list[str], tmp: Path, limit: float) -> dict:
        """runs pytest -x on tests with PYTHONPATH=src; exit status, seconds, the first
        failing test (its node id from the run's JUnit file; a collection error: the file),
        whether collection failed, the tail of the output"""
        tmp.mkdir(parents=True, exist_ok=True)
        args_file = tmp / "args.txt"
        args_file.write_text("\n".join(tests) + "\n", encoding="utf-8")
        junit = tmp / "junit.xml"
        per_test = max(120, int(limit))
        cmd = [
            sys.executable, "-m", "pytest", "-x", "-q", "-p", "no:cacheprovider",
            f"--basetemp={tmp / 'pt'}", f"--timeout={per_test}", "-rfE",
            f"--junitxml={junit}", f"@{args_file}",
        ]  # fmt: skip
        t0 = time.monotonic()
        with open(tmp / "out.txt", "w", encoding="utf-8") as out:
            proc = subprocess.Popen(
                cmd, cwd=self.root, env=self.env(src), stdout=out, stderr=subprocess.STDOUT,
                start_new_session=True,
            )  # fmt: skip
            if self.a.mem_gb and hasattr(os, "sched_getaffinity"):
                # (a memory cap for the mutant's pytest and the children it starts later;
                # set after the start: preexec_fn is not safe in a threaded process)
                import resource

                n = int(self.a.mem_gb * 2**30)
                try:
                    resource.prlimit(proc.pid, resource.RLIMIT_AS, (n, n))
                except (OSError, ValueError):
                    pass
            timed_out = False
            try:
                proc.wait(timeout=limit)
            except subprocess.TimeoutExpired:
                timed_out = True
                _kill_group(proc)
        text = (tmp / "out.txt").read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        first = self._first_failure(junit)
        collection = any(re.match(r"^ERROR \S+\.py( - |$)", x) for x in lines)
        if not first and collection:
            first = next(x.split()[1] for x in lines if re.match(r"^ERROR \S", x))
        return {
            "rc": None if timed_out else proc.returncode,
            "timeout": timed_out,
            "seconds": round(time.monotonic() - t0, 1),
            "first": first,
            "collection": collection,
            "tail": "\n".join(lines[-12:])[-3000:],
        }

    def _first_failure(self, junit: Path) -> str:
        try:
            cases = ET.parse(junit).getroot().iter("testcase")
        except (OSError, ET.ParseError):
            return ""
        for tc in cases:
            if any(ch.tag in ("failure", "error") for ch in tc):
                return junit_nodeid(tc.get("classname", ""), tc.get("name", ""), self.root)
        return ""

    def base_passes(self, test: str) -> bool:
        """whether a test passes on the base tree (cached)"""
        with self.base_lock:
            if test in self.base_results:
                return self.base_results[test]
        r = self.pytest(self.base_src, [test], self.work / "base_runs" / _safe(test), 1800)
        ok = r["rc"] == 0
        with self.base_lock:
            self.base_results[test] = ok
        return ok

    def one(self, p: dict) -> dict:
        m = M.Mutant.from_record(p)
        res = {"id": m.id, "area": m.area, "file": m.file, "line": m.line, "op": m.op}
        tests_given = self.recheck.get(m.id) if self.recheck else self.a.tests
        if not tests_given and not p.get("covered") and not p.get("area_tests"):
            res["status"] = "not_covered"
            return res
        mdir = self.work / "m" / m.id.replace("/", "_")
        if mdir.exists():
            shutil.rmtree(mdir)
        src = mdir / "src"
        try:
            shutil.copytree(self.base_src, src, copy_function=os.link)
            target = src / Path(m.file).relative_to("src")
            text = (self.base_src / Path(m.file).relative_to("src")).read_text(encoding="utf-8")
            mutated = M.apply_mutant(text, m)
            target.unlink()
            target.write_text(mutated, encoding="utf-8", newline="")
        except (OSError, ValueError) as e:
            res["status"] = "error"
            res["error"] = f"{type(e).__name__}: {e}"
            return res
        phases = (
            [("given", list(tests_given), self.a.given_budget)]
            if tests_given
            else [("area", p["area_tests"], p["area_seconds"])]
            + ([] if self.a.no_suite else [("suite", p["suite_tests"], p["suite_seconds"])])
        )
        res["status"] = "no_tests"
        for phase, tests, seconds in phases:
            if not tests:
                continue
            res["status"] = "survived"
            limit = min(self.a.max_seconds, 120 + 4 * float(seconds))
            r = self.pytest(src, tests, mdir / phase, limit)
            r["n"] = len(tests)
            res.setdefault("phases", {})[phase] = r
            if r["timeout"]:
                res["status"], res["by"] = "timeout", phase
                break
            if r["rc"] == 0:
                continue
            rc = r["rc"]
            if rc in (1, 2, 3) or (rc == 4 and r["collection"]) or (rc is not None and rc < 0):
                # (a failed test, a module that no longer imports, a crash)
                res["status"], res["by"] = "killed", phase
                first = r["first"]
                if first and not self.a.no_verify and "::" in first and rc == 1:
                    if not self.base_passes(first):
                        res["status"] = "base_fails"
                break
            res["status"] = "error"
            res["error"] = f"pytest exit {r['rc']} in the {phase} phase"
            break
        shutil.rmtree(mdir, ignore_errors=True)
        return res


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)[-150:]


def _kill_group(proc: subprocess.Popen) -> None:
    """ends the process group the tool started for a mutant (its pytest and the children)"""
    for sig, wait in ((signal.SIGTERM, 10), (signal.SIGKILL, 30)):
        try:
            os.killpg(proc.pid, sig)
        except (ProcessLookupError, PermissionError, AttributeError):
            return
        try:
            proc.wait(timeout=wait)
            return
        except subprocess.TimeoutExpired:
            continue


def cmd_run(a: argparse.Namespace) -> int:
    plan = read_jsonl(a.plan)
    if a.recheck:
        want = set(json.loads(Path(a.recheck).read_text(encoding="utf-8")))
        plan = [p for p in plan if p["id"] in want]
    if a.ids:
        want = set(a.ids.split(","))
        plan = [p for p in plan if p["id"] in want]
    if a.shard:
        i, n = (int(x) for x in a.shard.split("/"))
        plan = [p for k, p in enumerate(plan) if k % n == i - 1]
    done = set()
    if os.path.exists(a.out) and not a.tests and not a.recheck:
        done = {r["id"] for r in read_jsonl(a.out)}
    todo = [p for p in plan if p["id"] not in done]
    run = Runner(a)
    run.prepare()
    print(f"{len(todo)} mutants to run ({len(done)} done), {a.parallel} at once", flush=True)

    def task(p: dict) -> None:
        t0 = time.monotonic()
        try:
            r = run.one(p)
        except Exception as e:  # the tool's own failure: record it, go on
            r = {"id": p["id"], "area": p["area"], "status": "error", "error": repr(e)}
        r["wall"] = round(time.monotonic() - t0, 1)
        with run.lock:
            with open(a.out, "a", encoding="utf-8", newline="\n") as fh:
                fh.write(json.dumps(r) + "\n")
            print(f"{r['id']:18s} {r['status']:12s} {r['wall']:7.1f} s", flush=True)

    with ThreadPoolExecutor(max(1, a.parallel)) as ex:
        list(ex.map(task, todo))
    return 0


# --- report ----------------------------------------------------------------------------------


def cmd_report(a: argparse.Namespace) -> int:
    plan = {p["id"]: p for p in read_jsonl(a.plan)}
    results: dict[str, dict] = {}
    for f in a.results:
        for r in read_jsonl(f):
            results[r["id"]] = r  # (a later file wins: reruns)
    new_kills: dict[str, str] = {}  # id: the new test that killed it
    for f in a.rechecks or []:
        for r in read_jsonl(f):
            if r["status"] in KILLED:
                firsts = [ph.get("first", "") for ph in r.get("phases", {}).values()]
                new_kills[r["id"]] = next((t for t in firsts if t), "")
    disp = json.loads(Path(a.dispositions).read_text()) if a.dispositions else {}
    rows = defaultdict(Counter)
    open_rows = []
    uncompiled = []  # (C++ mutants that don't build: not mutants, left out of the counts)
    for mid, p in plan.items():
        r = results.get(mid)
        st = r["status"] if r else "not_run"
        if st == "uncompiled":
            uncompiled.append(mid)
            continue
        if r and st in KILLED:
            st = "killed_area" if r.get("by") in ("area", "given") else "killed_suite"
        elif mid in new_kills:
            st = "killed_new"
        elif mid in disp:
            rows[p["area"]]["reason"] += 1
            if disp[mid].get("kind") in EQUIVALENT_KINDS:
                rows[p["area"]]["reason_equivalent"] += 1
        rows[p["area"]][st] += 1
        if st not in ("killed_area", "killed_suite"):
            open_rows.append((p, st, disp.get(mid, {})))
    lines = [
        "| area | mutants | killed by area tests | by other tests | score | survived "
        "| not covered | killed by new tests | with a reason | open | score now "
        "| without equivalents |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    tot = Counter()
    summary = {}
    for area in [*sorted(rows), "all"]:
        c = tot if area == "all" else rows[area]
        if area != "all":
            tot.update(c)
        n = sum(v for k, v in c.items() if not k.startswith("reason"))
        k_old = c["killed_area"] + c["killed_suite"]
        rest = n - k_old - c["killed_new"]
        left = rest - c["reason"]
        lines.append(
            f"| {area} | {n} | {c['killed_area']} | {c['killed_suite']} | {_pct(k_old, n)} "
            f"| {c['survived'] + c['no_tests'] + c['base_fails'] + c['error'] + c['not_run']} "
            f"| {c['not_covered']} | {c['killed_new']} | {c['reason']} | {left} "
            f"| {_pct(k_old + c['killed_new'], n)} "
            f"| {_pct(k_old + c['killed_new'], n - c['reason_equivalent'])} |"
        )
        summary[area] = dict(c, mutants=n)
    if uncompiled:
        lines += ["", f"Left out, as they don't compile: {', '.join(sorted(uncompiled))}."]
    lines += ["", "| id | status | where | mutation | disposition |", "|---|---|---|---|---|"]
    for p, st, d in sorted(open_rows, key=lambda x: x[0]["id"]):
        where = f"`{p['file'].removeprefix('src/petakit5d/')}:{p['line']}`"
        mut = f"`{p['old']}` -> `{p['new']}` ({p['op']})"
        if st == "killed_new":
            kind, why = "test", f"`{new_kills[p['id']]}`" if new_kills[p["id"]] else ""
        else:
            kind, why = d.get("kind", "**open**"), d.get("why", "")
        test = f" `{d['test']}`" if d.get("test") else ""
        lines.append(f"| {p['id']} | {st} | {where} | {_md(mut)} | {kind}: {_md(why)}{test} |")
    text = "\n".join(lines) + "\n"
    if a.out:
        Path(a.out).write_text(text, encoding="utf-8", newline="\n")
    if a.json:
        Path(a.json).write_text(json.dumps(summary, indent=1), encoding="utf-8", newline="\n")
    print(text)
    return 0


def _pct(k: int, n: int) -> str:
    return f"{100 * k / n:.0f}%" if n else "-"


def _md(s: str) -> str:
    return s.replace("|", "\\|").replace("\n", " ")


# --- show, tree ------------------------------------------------------------------------------


def _find(mutants_file: str, mid: str) -> M.Mutant:
    for rec in read_jsonl(mutants_file):
        if rec["id"] == mid:
            return M.Mutant.from_record(rec)
    raise SystemExit(f"no mutant {mid} in {mutants_file}")


def cmd_show(a: argparse.Namespace) -> int:
    for mid in a.id:
        m = _find(a.mutants, mid)
        print(f"{m.id} ({m.op}) {m.file}:{m.line}\n- {m.src_line}\n+ {m.mut_line}")
    return 0


def cmd_tree(a: argparse.Namespace) -> int:
    root = Path(a.root)
    m = _find(a.mutants, a.id)
    dest = Path(a.dest)
    shutil.copytree(root / "src", dest / "src", ignore=shutil.ignore_patterns("__pycache__"))
    text = (root / m.file).read_text(encoding="utf-8")
    (dest / m.file).write_text(M.apply_mutant(text, m), encoding="utf-8", newline="")
    print(f"{dest / 'src'}: {m.file}:{m.line} mutated (PYTHONPATH={dest / 'src'})")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--root", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("--seed", default="0")
    g.add_argument("--per-area", default="auto")
    g.add_argument("--min", type=int, default=12)
    g.add_argument("--max", type=int, default=50)
    g.add_argument("--area", action="append")
    g.add_argument("--only", action="append")
    g.add_argument("--exclude", action="append")
    g.add_argument("--ops", default="")
    g.add_argument("--prefix", default="", help="put before the ids (as p2/ for a second pass)")
    p = sub.add_parser("plan")
    p.add_argument("--root", required=True)
    p.add_argument("--mutants", required=True)
    p.add_argument("--coverage", required=True)
    p.add_argument("--junit", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--area-budget", type=float, default=1200)
    p.add_argument("--suite-budget", type=float, default=900)
    p.add_argument(
        "--header-tests",
        action="append",
        help="a test file whose tests join every signature, decorator and module-level line's "
        "(one that checks signatures against a reference, as tests/params/test_check.py)",
    )
    r = sub.add_parser("run")
    r.add_argument("--root", required=True)
    r.add_argument("--plan", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--work", default=os.path.join(os.environ.get("TMPDIR", "/tmp"), "pk5d_mut"))
    r.add_argument("--parallel", type=int, default=2)
    r.add_argument("--threads", type=int, default=2)
    r.add_argument("--shard", default="")
    r.add_argument("--ids", default="")
    r.add_argument("--tests", action="append", default=[])
    r.add_argument("--recheck", default="")
    r.add_argument("--given-budget", type=float, default=600)
    r.add_argument("--max-seconds", type=float, default=3600)
    r.add_argument("--no-suite", action="store_true")
    r.add_argument("--no-verify", action="store_true")
    r.add_argument("--mem-gb", type=float, default=0)
    rep = sub.add_parser("report")
    rep.add_argument("--plan", required=True)
    rep.add_argument("--results", nargs="+", required=True)
    rep.add_argument("--rechecks", nargs="*", default=[])
    rep.add_argument("--dispositions", default="")
    rep.add_argument("--out", default="")
    rep.add_argument("--json", default="")
    s = sub.add_parser("show")
    s.add_argument("--root", required=True)
    s.add_argument("--mutants", required=True)
    s.add_argument("id", nargs="+")
    t = sub.add_parser("tree")
    t.add_argument("--root", required=True)
    t.add_argument("--mutants", required=True)
    t.add_argument("id")
    t.add_argument("--dest", required=True)
    a = ap.parse_args(argv)
    cmds = {
        "generate": cmd_generate,
        "plan": cmd_plan,
        "run": cmd_run,
        "report": cmd_report,
        "show": cmd_show,
        "tree": cmd_tree,
    }
    return cmds[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
