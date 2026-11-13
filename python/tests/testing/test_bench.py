"""Tests of the benchmark harness (petakit5d.testing.bench): marks and phases, statistics,
comparisons, gates parsed from task texts, a real interleaved run of two small Python
variants, and the history table.

    python -m pytest tests/testing/test_bench.py
"""

import json
import os
import sys

import pytest

from petakit5d.testing import bench


def _run(variant, launch, first, ready, last, exit_, rnd, records=(), rc=0, warmup=False):
    marks = {"first": first, "last": last}
    if ready is not None:
        marks["ready"] = ready
    return {
        "variant": variant,
        "impl": "python",
        "round": rnd,
        "warmup": warmup,
        "rc": rc,
        "launch": launch,
        "exit": exit_,
        "marks": marks,
        "records": list(records),
        "maxrss_kb": 1024 * 100,
        "user_s": 0.1,
        "sys_s": 0.05,
    }


def _result(runs, variants=("a", "b"), compare=None, gate=None):
    spec = {"name": "t", "variants": [{"id": v, "impl": "python"} for v in variants]}
    if compare is not None:
        spec["compare"] = compare
    if gate is not None:
        spec["gate"] = gate
    return {
        "schema": bench.SCHEMA,
        "name": "t",
        "task": "X",
        "started": "20260101_000000",
        "spec": spec,
        "runs": runs,
        "node": {
            "host": "node12",
            "node_type": "c24",
            "cpu_model": "cpu",
            "cpus_affinity": 4,
            "slurm": {},
            "loadavg": ["0", "0", "0"],
        },
        "source": {"commit": "abc", "dirty": "clean"},
    }


def test_phases_from_marks():
    p = bench.phases(_run("a", 10.0, 10.5, 11.0, 13.0, 13.25, 0))
    assert p["start"] == pytest.approx(0.5)
    assert p["setup"] == pytest.approx(0.5)
    assert p["startup"] == pytest.approx(1.0)
    assert p["task"] == pytest.approx(2.0)
    assert p["teardown"] == pytest.approx(0.25)
    assert p["overhead"] == pytest.approx(1.25)
    assert p["total"] == pytest.approx(3.25)


def test_phases_without_ready_start_at_first_and_own_marks_are_steps():
    r = _run("a", 0.0, 1.0, None, 4.0, 5.0, 0)
    r["marks"]["context"] = 2.5
    p = bench.phases(r)
    assert p["setup"] is None
    assert p["startup"] == pytest.approx(1.0)
    assert p["task"] == pytest.approx(3.0)
    assert p["seg:first->context"] == pytest.approx(1.5)
    assert p["seg:context->last"] == pytest.approx(1.5)


def test_missing_mark_gives_none():
    r = _run("a", 0.0, 1.0, None, 4.0, 5.0, 0)
    del r["marks"]["last"]
    p = bench.phases(r)
    assert p["task"] is None and p["teardown"] is None and p["overhead"] is None
    assert p["total"] == pytest.approx(5.0)


def test_stats_median_quartiles_and_range():
    s = bench.stats([5.0, 1.0, 3.0, 2.0, 4.0])
    assert s["median"] == 3.0 and s["min"] == 1.0 and s["max"] == 5.0 and s["n"] == 5
    assert s["q1"] == pytest.approx(1.75) and s["q3"] == pytest.approx(4.25)
    assert bench.stats([]) is None
    assert bench.stats([2.0, 4.0])["median"] == 3.0


def _ab_runs(fa=1.0, fb=1.0, n=5):
    runs = []
    for r in range(n):
        for v, f in (("a", fa), ("b", fb)):
            t = 100.0 * r
            runs.append(
                _run(
                    v,
                    t,
                    t + 0.1 * f,
                    t + 0.2 * f,
                    t + 1.0 * f,
                    t + 1.1 * f,
                    r,
                    records=[
                        {"op": "x", "threads": 4, "times": [f * (1 + 0.01 * r)] * 3},
                        {"op": "y", "threads": 1, "times": [2 * f]},
                    ],
                )
            )
    runs.append(_run("a", 0, 5, 6, 7, 50, -1, warmup=True))
    runs.append(_run("b", 0, 1, 2, 3, 4, 0, rc=1))
    return runs


def test_summary_pools_ops_and_keeps_warmup_and_failures_apart():
    s = bench.summarize(_result(_ab_runs()))
    assert s["a"]["n"] == 5 and s["b"]["n"] == 5 and s["b"]["failed"] == 1
    assert s["a"]["first"][0]["teardown"] == pytest.approx(43.0)
    assert s["a"]["phases"]["teardown"]["median"] == pytest.approx(0.1)
    op = s["a"]["ops"]["x@4"]
    assert op["n"] == 15 and op["median"] == pytest.approx(1.02)
    assert op["round_min"] == pytest.approx(1.0) and op["round_max"] == pytest.approx(1.04)
    assert s["a"]["maxrss_mb"]["median"] == pytest.approx(100.0)


def test_compare_ratios_and_op_map():
    res = _result(_ab_runs(fa=1.2, fb=1.0), compare=[{"a": "a", "b": "b", "op_map": {"y": "y"}}])
    rows = {r["metric"]: r for r in bench.compare(res)[0]["rows"]}
    assert rows["op:x@4"]["ratio"] == pytest.approx(1.2)
    assert rows["op:x@4"]["paired"] == pytest.approx(1.2)
    assert rows["op:y@1"]["ratio"] == pytest.approx(1.2)
    assert rows["task"]["ratio"] == pytest.approx(1.2)
    assert rows["maxrss"]["ratio"] == pytest.approx(1.0)
    # an op whose name differs on the reference side
    res = _result(_ab_runs(), compare=[{"a": "a", "b": "b", "op_map": {"x": "y"}}])
    rows = {r["metric"]: r for r in bench.compare(res)[0]["rows"]}
    assert "op:x@4" not in rows  # y was measured at 1 thread only


def test_default_pairs_are_python_against_the_one_matlab_variant():
    res = _result([], variants=("p1", "p2", "m"))
    res["spec"]["variants"][2]["impl"] = "matlab"
    assert bench._pairs(res) == [{"a": "p1", "b": "m"}, {"a": "p2", "b": "m"}]
    res["spec"]["variants"][1]["impl"] = "matlab"
    assert bench._pairs(res) == []


@pytest.mark.parametrize(
    ("text", "want"),
    [
        (
            "no slower than the mex (1.05x) at 1, 4, 24 threads",
            {"ratio_max": 1.05, "threads": [1, 4, 24]},
        ),
        (
            "Python no slower than MATLAB: median ratio at most 1.05 at 1, 4 and 24 threads",
            {"ratio_max": 1.05, "threads": [1, 4, 24]},
        ),
        ("iterations no slower than the MATLAB GPU mex (1.10x) on A100", {"ratio_max": 1.10}),
        ("no slower than MATLAB on a 512 x 1800 x 500 tile", {"ratio_max": 1.05}),
        ("compute as above, and start-up + teardown at most 1 s (CPU)", {"overhead_max": 1.0}),
        ("peak memory at most 1.1 x MATLAB's", {"mem_ratio_max": 1.1}),
    ],
)
def test_parse_gate(text, want):
    g = bench.parse_gate(text)
    for k, v in want.items():
        assert g[k] == pytest.approx(v), (k, g)
    assert g["manual"]


def test_parse_gate_none():
    assert bench.parse_gate("none (no compute path of its own)")["none"]
    assert bench.evaluate_gate(_result([]), {"none": True})["met"] is None


def test_task_gate_reads_tasks_json(tmp_path):
    p = tmp_path / "tasks.json"
    p.write_text(
        json.dumps({"tasks": [{"id": "K9", "benchmark_gate": "no slower (1.02x)"}]}), newline="\n"
    )
    assert bench.task_gate("K9", str(p))["ratio_max"] == pytest.approx(1.02)
    assert bench.task_gate("K8", str(p)) is None


def test_gate_met_and_missed():
    res = _result(_ab_runs(fa=1.04))
    res["spec"]["compare"] = [["a", "b"]]
    ok = bench.evaluate_gate(res, {"ratio_max": 1.05})
    assert ok["met"] is True and len(ok["checks"]) == 2
    bad = bench.evaluate_gate(res, {"ratio_max": 1.03})
    assert bad["met"] is False
    assert {c["metric"] for c in bad["checks"] if not c["ok"]} == {"op:x@4", "op:y@1"}
    missing = bench.evaluate_gate(res, {"ratio_max": 1.05, "threads": [4]})
    assert missing["met"] is False  # y was measured at 1 thread only
    assert [c["metric"] for c in missing["checks"] if not c["ok"]] == ["threads:y"]
    over = bench.evaluate_gate(res, {"overhead_max": 0.2})
    assert over["met"] is False and over["checks"][0]["pair"] == "a"
    assert bench.evaluate_gate(res, {"overhead_max": 0.5})["met"] is True


def test_gate_metric_that_matches_no_row_fails():
    # N06's case: rows carry @threads, so "op:x" matches nothing and only the thread counts were
    # checked: the gate said met for an op 1.3x slower
    only_x = [
        {**r, "records": [x for x in r["records"] if x["op"] == "x"]} for r in _ab_runs(fa=1.3)
    ]
    res = _result(only_x, compare=[["a", "b"]])
    g = bench.evaluate_gate(res, {"ratio_max": 1.05, "threads": [4], "metrics": ["op:x"]})
    assert g["met"] is False
    res = _result(_ab_runs(fa=1.3), compare=[["a", "b"]])
    g = bench.evaluate_gate(res, {"ratio_max": 1.05, "threads": [4], "metrics": ["op:x"]})
    assert g["met"] is False
    bad = [c for c in g["checks"] if not c["ok"]]
    assert [c["metric"] for c in bad] == ["unmatched:op:x"]
    assert "op:x@4" in bad[0]["limit"]
    # the right name checks the ratio (and only the gated op's threads)
    g = bench.evaluate_gate(res, {"ratio_max": 1.05, "threads": [4], "metrics": ["op:x@*"]})
    assert g["met"] is False
    assert {c["metric"] for c in g["checks"]} == {"op:x@4", "threads:x"}
    assert not next(c for c in g["checks"] if c["metric"] == "op:x@4")["ok"]
    res = _result(_ab_runs(fa=1.0), compare=[["a", "b"]])
    assert bench.evaluate_gate(res, {"ratio_max": 1.05, "threads": [4], "metrics": ["op:x@*"]})[
        "met"
    ]
    # one name matching, another not: still fails
    g = bench.evaluate_gate(res, {"ratio_max": 1.05, "metrics": ["op:x@*", "op:z@*"]})
    assert g["met"] is False and [c["metric"] for c in g["checks"] if not c["ok"]] == [
        "unmatched:op:z@*"
    ]


def test_gate_on_ratios_without_comparisons_or_ops_fails():
    res = _result(_ab_runs(), compare=[])
    g = bench.evaluate_gate(res, {"ratio_max": 1.05})
    assert g["met"] is False and g["checks"][0]["metric"] == "unmatched:comparison"
    # compared, but no operation timed on both sides
    runs = [r for r in _ab_runs() if r["variant"] == "a"] + [
        {**r, "records": []} for r in _ab_runs() if r["variant"] == "b"
    ]
    g = bench.evaluate_gate(_result(runs, compare=[["a", "b"]]), {"ratio_max": 1.05})
    assert g["met"] is False and g["checks"][0]["metric"] == "unmatched:op:*"
    # a gate on process costs alone needs no comparison, and its metrics select nothing
    assert bench.evaluate_gate(res, {"overhead_max": 0.5})["met"] is True
    assert bench.evaluate_gate(res, {"overhead_max": 0.5, "metrics": ["startup"]})["met"] is True
    compared = _result(_ab_runs(), compare=[["a", "b"]])
    g = bench.evaluate_gate(compared, {"overhead_max": 0.5, "metrics": ["op:none"]})
    assert g["met"] is True  # metrics only select rows for ratio checks


def test_gate_within_for_a_build_against_itself():
    res = _result(_ab_runs(fa=1.02), compare=[["a", "b"]])
    assert bench.evaluate_gate(res, {"within": 0.03, "metrics": ["op:*", "task"]})["met"]
    res = _result(_ab_runs(fa=0.96), compare=[["a", "b"]])
    g = bench.evaluate_gate(res, {"within": 0.03, "metrics": ["op:*"]})
    assert g["met"] is False and len(g["checks"]) == 2


def test_placeholders_leave_other_braces_alone(monkeypatch):
    monkeypatch.setenv("PK5D_X", "/opt/x")
    ctx = {"work": "/w", "src": "/s"}
    assert bench._fmt("a{work}b {x{1}} {env[PK5D_X]}", ctx) == "a/wb {x{1}} /opt/x"
    assert bench._fmt(["{src}/f", {"k": "{work}"}], ctx) == ["/s/f", {"k": "/w"}]
    monkeypatch.delenv("PK5D_X")
    with pytest.raises(ValueError):
        bench._fmt("{env[PK5D_X]}", ctx)


def test_matlab_command_wraps_the_code_in_marks():
    argv = bench._matlab_argv(
        {"code": "y = x{1}", "paths": ["benchmarks/io"]}, {"src": "/s", "work": "/w"}
    )
    stmt = argv[-1]
    assert argv[-2] == "-batch" and "-nojvm" in argv
    assert stmt.index("'first'") < stmt.index("setup([])") < stmt.index("'ready'")
    assert stmt.index("'ready'") < stmt.index("y = x{1};") < stmt.index("'last'")
    assert "addpath('/s/benchmarks/common', '/s/benchmarks/io')" in stmt
    nosetup = bench._matlab_argv({"code": "", "setup": False, "flags": []}, {"src": "/s"})
    assert "setup([])" not in nosetup[-1] and nosetup[1] == "-batch"
    r = bench._matlab_argv({"code": "x = 1", "launch": "r"}, {"src": "/s"})
    assert r[-2] == "-r" and r[-1].startswith("try, ") and r[-1].endswith("exit(0)")
    assert r[-1].index("'last'") < r[-1].index("exit(0)")


def test_mark_and_record_write_lines(tmp_path, monkeypatch):
    f = tmp_path / "m.jsonl"
    bench.mark("ready")  # no file: nothing happens
    monkeypatch.setenv("PK5D_BENCH_FILE", str(f))
    bench.mark("ready")
    bench.record("op", [0.5, 0.25], threads=2, size="8x8")
    bench.info(x=1)
    marks, recs, inf = bench._read_marks(str(f))
    assert "ready" in marks and inf == {"x": 1}
    assert recs == [{"op": "op", "threads": 2, "times": [0.5, 0.25], "size": "8x8"}]


def test_timeit_skips_warmup():
    calls = []
    t = bench.timeit(lambda: calls.append(1), repeats=3, warmup=2, before=lambda: calls.append(0))
    assert len(t) == 3 and calls == [0, 1] * 5


class _SlowFree:
    def __del__(self):
        __import__("time").sleep(0.05)


def test_timeit_frees_the_result_after_the_timer():
    t = bench.timeit(_SlowFree, repeats=3, warmup=0)
    assert max(t) < 0.04, t
    t = bench.timeit(_SlowFree, repeats=3, warmup=0, include_free=True)
    assert min(t) >= 0.045, t


BUDDY = """\
Node 0, zone      DMA      0      0      0      0      0      0      0      0      1      1      2
Node 0, zone    DMA32      7      5      8      8     11      6      8      6      6      7    395
Node 0, zone   Normal  32599   8188   2556   1018    401    193     50     31  56121     12      3
Node 1, zone   Normal   3472   1722    897    697    388    149     69  56045  52419  28182  63065
"""


def _hp_files(tmp_path, buddy=BUDDY, mems="0-1"):
    (tmp_path / "buddyinfo").write_text(buddy, newline="\n")
    thp = tmp_path / "thp"
    thp.mkdir(exist_ok=True)
    (thp / "enabled").write_text("always [madvise] never\n", newline="\n")
    (thp / "defrag").write_text("always defer defer+madvise [madvise] never\n", newline="\n")
    (tmp_path / "status").write_text(f"Name:\tpython\nMems_allowed_list:\t{mems}\n", newline="\n")
    return str(tmp_path / "buddyinfo"), str(thp), str(tmp_path / "status")


def test_hugepage_info_counts_free_2mb_blocks_per_numa_node(tmp_path):
    hp = bench.hugepage_info(*_hp_files(tmp_path))
    # order 9 counts once, order 10 twice; all zones of a node add up
    assert hp["free_2mb_blocks"] == {
        "0": (1 + 2 * 2) + (7 + 2 * 395) + (12 + 2 * 3),
        "1": 28182 + 2 * 63065,
    }
    assert hp["thp_enabled"] == "madvise" and hp["thp_defrag"] == "madvise"
    assert hp["mems_allowed"] == [0, 1]
    assert (
        bench.hugepage_info(str(tmp_path / "none"), str(tmp_path / "none"))["free_2mb_blocks"] == {}
    )


def test_warning_when_a_numa_node_is_short_of_2mb_blocks(tmp_path):
    low = BUDDY.replace("395", "0").replace("    12      3", "     5      0")
    hp_low = bench.hugepage_info(*_hp_files(tmp_path, low))
    assert hp_low["free_2mb_blocks"]["0"] == 1 + 4 + 7 + 5
    res = _result(_ab_runs(), compare=[["a", "b"]])
    res["node"]["hugepages"] = hp_low
    w = bench.hugepage_warnings(res)
    assert len(w) == 1 and w[0].startswith(
        "NUMA node 0 had 17 free 2 MB blocks (34 MB) at the start"
    )
    assert "4 KB pages" in w[0] and "**Warning:** NUMA node 0" in bench.report(res)
    # enough blocks for a small output: no warning; a node the job may not use: no warning
    assert bench.hugepage_warnings(res, need_bytes=32 << 20) == []
    assert len(bench.hugepage_warnings(res, need_bytes=40 << 20)) == 1  # 34 MB < 40 MB
    res["node"]["hugepages"] = bench.hugepage_info(*_hp_files(tmp_path, low, mems="1"))
    assert bench.hugepage_warnings(res) == []
    # the end of the run counts too, and the spec can say how large an output is
    res["node"]["hugepages"] = bench.hugepage_info(*_hp_files(tmp_path))
    res["node_after"] = {"hugepages": hp_low}
    res["spec"]["output_bytes"] = 20 << 20
    assert bench.hugepage_warnings(res) == []
    res["spec"]["output_bytes"] = 100 << 20  # start: 1,640 MB free; end: 34 MB
    assert [x.split(" at the ")[1][:3] for x in bench.hugepage_warnings(res)] == ["end"]


def test_node_info_records_hugepages():
    d = bench.node_info()
    assert set(d["hugepages"]) == {"free_2mb_blocks", "thp_enabled", "thp_defrag", "mems_allowed"}


def test_load_spec_checks(tmp_path):
    with pytest.raises(ValueError):
        bench.load_spec(
            {
                "name": "x",
                "variants": [{"id": "a", "impl": "python"}, {"id": "a", "impl": "python"}],
            }
        )
    with pytest.raises(ValueError):
        bench.load_spec({"name": "x", "variants": [{"id": "a", "impl": "java"}]})
    p = tmp_path / "s.py"
    p.write_text(
        "SPEC = {'name': 'n', 'variants': [{'id': 'a', 'impl': 'python', 'code': ''}]}",
        newline="\n",
    )
    assert bench.load_spec(str(p))["name"] == "n"
    with pytest.raises(ValueError):
        bench.run(
            {"name": "x", "repeats": 3, "variants": [{"id": "a", "impl": "python", "code": ""}]}
        )


def test_run_interleaves_and_writes_a_result(tmp_path, monkeypatch):
    monkeypatch.setenv("PK5D_SRC", str(tmp_path))
    script = tmp_path / "s.py"
    script.write_text(
        "import sys, time\n"
        "from petakit5d.testing.bench import mark, record\n"
        "mark('ready')\n"
        "time.sleep(0.02)\n"
        "record('sleep', [0.01, 0.01], threads=int(sys.argv[1]))\n"
        "if sys.argv[2] == 'fail': sys.exit(3)\n",
        newline="\n",
    )
    pp = os.pathsep.join(p for p in sys.path if p)
    spec = {
        "name": "unit",
        "repeats": 5,
        "variants": [
            {
                "id": "a",
                "impl": "python",
                "script": "s.py",
                "args": ["2", "ok"],
                "env": {"PYTHONPATH": pp},
            },
            {
                "id": "b",
                "impl": "python",
                "code": "import sys; assert sys.argv[1:] == ['x']",
                "args": ["x"],
            },
            {
                "id": "c",
                "impl": "python",
                "script": "s.py",
                "args": ["1", "fail"],
                "env": {"PYTHONPATH": pp},
            },
            {"id": "d", "impl": "command", "argv": ["true"], "when": "PK5D_NOT_SET_ANYWHERE"},
        ],
        "compare": [["a", "b"]],
        "gate": {"ratio_max": 1e9, "metrics": ["task"]},
    }
    lines = []
    r = bench.run(spec, out_dir=str(tmp_path / "out"), log=lines.append)
    seq = [(x["variant"], x["round"]) for x in r["runs"]]
    assert seq[:3] == [("a", -1), ("b", -1), ("c", -1)]
    assert seq[3:6] == [("a", 0), ("b", 0), ("c", 0)] and seq[6:9] == [("c", 1), ("b", 1), ("a", 1)]
    assert r["skipped"] == [{"id": "d", "why": "PK5D_NOT_SET_ANYWHERE not set"}]
    s = r["summary"]
    assert s["a"]["n"] == 5 and s["b"]["n"] == 5 and s["c"]["n"] == 0 and s["c"]["failed"] == 6
    assert s["a"]["phases"]["task"]["median"] >= 0.02
    assert s["a"]["ops"]["sleep@2"]["n"] == 10
    assert all(x["rc"] == 3 for x in r["runs"] if x["variant"] == "c")
    assert "log_tail" in r["runs"][2]
    assert r["gate"]["met"] is True
    with open(r["path"]) as f:
        saved = json.load(f)
    assert saved["schema"] == bench.SCHEMA and len(saved["runs"]) == 18
    logs = os.listdir(r["path"][: -len(".json")])
    assert sum(n.startswith("marks_") for n in logs) == 18
    text = bench.report(saved)
    assert "| a | 5 |" in text and "a / b:" in text


def test_evict_counts_files(tmp_path):
    (tmp_path / "d").mkdir()
    for i in range(3):
        (tmp_path / "d" / f"f{i}").write_bytes(b"x" * 10)
    e = bench.evict([str(tmp_path)])
    if e["supported"]:
        assert e["files"] == 3


def test_warm_variants_are_rewarmed_after_a_cold_one(tmp_path, monkeypatch):
    monkeypatch.setenv("PK5D_SRC", str(tmp_path))
    (tmp_path / "d").mkdir()
    (tmp_path / "d" / "f").write_text("x", newline="\n")
    spec = {
        "name": "cold",
        "repeats": 5,
        "variants": [
            {"id": "warm", "impl": "python", "code": ""},
            {"id": "cold", "impl": "python", "code": "", "evict": [str(tmp_path / "d")]},
            {"id": "warm2", "impl": "python", "code": ""},
        ],
    }
    r = bench.run(spec, log=lambda s: None)
    seq = [(x["variant"], x["round"], x["rewarm"]) for x in r["runs"]]
    # a warm variant runs untimed first when the cold one ran since its last run, else not
    last = {}
    for i, (vid, rnd, again) in enumerate(seq):
        if vid.startswith("warm") and not again:
            j = last.get(vid, -1)
            evicted = any(v == "cold" for v, _, _ in seq[j + 1 : i])
            assert (seq[i - 1] == (vid, rnd, True)) == evicted, (i, seq)
            last[vid] = i
    # rounds alternate (warm cold warm2 / warm2 cold warm): after the warm-up, warm2 in rounds
    # 0, 2, 4 and warm in 0, 1, 3
    assert sum(again for _, _, again in seq) == 7
    assert not any(again for v, _, again in seq if v == "cold")
    s = r["summary"]
    assert s["warm"]["n"] == 5 and s["warm2"]["n"] == 5 and s["cold"]["n"] == 5
    assert len(s["warm2"]["first"]) == 1  # the re-warm runs are not first runs
    assert all(x["evicted"]["files"] == 1 for x in r["runs"] if x["variant"] == "cold")


@pytest.mark.skipif(os.name != "posix", reason="process groups and sh")
def test_timeout_stops_the_process_group(tmp_path, monkeypatch):
    monkeypatch.setenv("PK5D_SRC", str(tmp_path))
    spec = {
        "name": "slow",
        "repeats": 1,
        "warmup": 0,
        "allow_few": True,
        "variants": [
            {
                "id": "s",
                "impl": "command",
                "argv": ["sh", "-c", "sleep 60 & sleep 60"],
                "timeout": 1,
            },
        ],
    }
    t0 = __import__("time").time()
    r = bench.run(spec, log=lambda s: None)
    assert __import__("time").time() - t0 < 20
    run = r["runs"][0]
    assert run["rc"] != 0 and run["timed_out"] == pytest.approx(1, abs=0.5)


def test_history_flags_a_ratio_that_grew(tmp_path):
    for i, (fa, commit) in enumerate([(1.0, "c1"), (1.02, "c2"), (1.2, "c3")]):
        res = _result(_ab_runs(fa=fa), compare=[["a", "b"]])
        res["started"] = f"2026010{i + 1}_000000"
        res["source"]["commit"] = commit
        res["summary"] = bench.summarize(res)
        res["comparisons"] = bench.compare(res, res["summary"])
        (tmp_path / f"r{i}.json").write_text(json.dumps(res), newline="\n")
    (tmp_path / "other.json").write_text("{}", newline="\n")
    text = bench.history(str(tmp_path), write=True)
    rows = json.loads((tmp_path / "history.json").read_text())
    x = [r for r in rows if r["metric"] == "a:op:x@4"]
    assert [r["commit"] for r in x] == ["c1", "c2", "c3"]
    assert x[1]["flag"] == "" and x[2]["flag"] == "REGRESSION"
    assert x[2]["change"] == pytest.approx(1.2 / 1.02 - 1)
    assert "REGRESSION" in text and (tmp_path / "HISTORY.md").exists()
    # another node type starts its own line of results
    res = _result(_ab_runs(fa=2.0), compare=[["a", "b"]])
    res["started"], res["node"]["node_type"] = "20260109_000000", "c48"
    res["comparisons"] = bench.compare(res)
    (tmp_path / "r9.json").write_text(json.dumps(res), newline="\n")
    bench.history(str(tmp_path), write=True)
    rows = json.loads((tmp_path / "history.json").read_text())
    c48 = [r for r in rows if r["metric"] == "a:op:x@4" and r["ntype"] == "c48"]
    assert len(c48) == 1 and c48[0]["change"] is None and c48[0]["flag"] == ""
    # so does a run with interleaved memory on the same node type
    res = _result(_ab_runs(fa=2.0), compare=[["a", "b"]])
    res["started"], res["node"]["numa_policy"] = "20260110_000000", "interleave:0-1"
    res["comparisons"] = bench.compare(res)
    (tmp_path / "r10.json").write_text(json.dumps(res), newline="\n")
    bench.history(str(tmp_path), write=True)
    rows = json.loads((tmp_path / "history.json").read_text())
    ilv = [r for r in rows if r["metric"] == "a:op:x@4" and r["ntype"] == "c24 ilv"]
    assert len(ilv) == 1 and ilv[0]["change"] is None
    assert "memory interleave:0-1" in bench.report(res)


def test_history_ignores_changes_of_a_few_milliseconds(tmp_path):
    for i, fa in enumerate((0.001, 0.0015)):
        res = _result(_ab_runs(fa=fa), compare=[["a", "b"]])
        res["started"] = f"2026010{i + 1}_000000"
        res["comparisons"] = bench.compare(res)
        (tmp_path / f"r{i}.json").write_text(json.dumps(res), newline="\n")
    bench.history(str(tmp_path), write=True)
    rows = json.loads((tmp_path / "history.json").read_text())
    x = [r for r in rows if r["metric"] == "a:op:x@4"]
    assert x[1]["change"] == pytest.approx(0.5) and x[1]["flag"] == ""


def test_marks_first_and_ready_keep_the_first_last_keeps_the_last(tmp_path):
    f = tmp_path / "m.jsonl"
    f.write_text(
        "\n".join(
            json.dumps({"mark": k, "t": t})
            for k, t in [
                ("first", 1.0),
                ("ready", 2.0),
                ("first", 3.0),
                ("ready", 4.0),
                ("last", 5.0),
                ("last", 6.0),
            ]
        )
        + "\nnot json\n",
        newline="\n",
    )
    marks, _, _ = bench._read_marks(str(f))
    assert marks == {"first": 1.0, "ready": 2.0, "last": 6.0}


def test_gate_limit_is_inclusive():
    res = _result(_ab_runs(fa=1.25), compare=[["a", "b"]])
    ratio = bench.compare(res)[0]["rows"]
    r = next(x["ratio"] for x in ratio if x["metric"] == "op:x@4")
    assert bench.evaluate_gate(res, {"ratio_max": r, "metrics": ["op:x@4"]})["met"] is True


def test_round_medians_are_medians():
    runs = [
        _run("a", 0, 1, 2, 3, 4, r, records=[{"op": "x", "times": [1.0, 2.0, 9.0]}])
        for r in range(5)
    ]
    op = bench.summarize(_result(runs, variants=("a",)))["a"]["ops"]["x"]
    assert op["round_min"] == 2.0 and op["round_max"] == 2.0 and op["median"] == 2.0


def test_unknown_placeholders_stay():
    assert bench._fmt("{nothing} {work}", {"work": "/w"}) == "{nothing} /w"


def test_cold_variants_are_never_rewarmed(tmp_path, monkeypatch):
    monkeypatch.setenv("PK5D_SRC", str(tmp_path))
    (tmp_path / "f").write_text("x", newline="\n")
    spec = {
        "name": "cold2",
        "repeats": 5,
        "variants": [
            {"id": c, "impl": "python", "code": "", "evict": [str(tmp_path)]} for c in ("c1", "c2")
        ],
    }
    r = bench.run(spec, log=lambda s: None)
    assert not any(x["rewarm"] for x in r["runs"]) and len(r["runs"]) == 12
