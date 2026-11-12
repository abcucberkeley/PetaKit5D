"""Benchmark harness: MATLAB and Python processes interleaved on one node, with process costs.

A benchmark is a spec (JSON, or a Python file that defines ``SPEC``) with variants: processes
that run MATLAB code (``matlab -batch`` from the oracle tree), a Python script, module or code,
or any command. ``run`` starts every variant once to warm up, then in rounds (5 or more) in an
order that alternates, so drift and the node's state hit both sides alike. Each process gets
``PK5D_BENCH_FILE``, a JSON-lines file where it leaves time marks and timings:

- the harness: ``launch`` (just before the process starts) and ``exit`` (when it has been
  reaped, with its peak memory and CPU time);
- ``first``: the process's first line (Python: the harness's bootstrap, before any import;
  MATLAB: the first statement of ``-batch``);
- ``ready``: the first line of the work (Python: the script calls ``mark('ready')`` after its
  imports; MATLAB: after ``setup([])``);
- ``last``: the last line (Python: when the script or module returns; MATLAB: after the code);
- timings of calls inside the process (``record``; MATLAB ``pk5d_bench_record``).

From the marks: ``start`` (launch to first), ``setup`` (first to ready), ``startup`` (launch to
ready: process start to the first line of the work), ``task`` (ready to last), ``teardown``
(last to exit), ``overhead`` (startup + teardown) and ``total``. The result (JSON) holds every
run, the node (host, CPU model, CPUs, memory, load, huge pages, GPUs), the versions and the
source commit, medians and spreads, the comparisons and the gate. ``history`` makes a table
of all results so changes between tasks show. Details: ``python/docs/benchmarks.md``.

In a benchmark script (cheap to import; does nothing without ``PK5D_BENCH_FILE``)::

    from petakit5d.testing.bench import mark, record, timeit
    import numpy as np
    mark("ready")
    record("fftn", timeit(lambda: fftn(x, threads=4), repeats=5), threads=4)

Command line: ``python -m petakit5d.testing.bench run|summary|gate|history|node ...``.
"""

from __future__ import annotations

import os
import sys
import time

SCHEMA = "pk5d-bench/1"
MATLAB_FLAGS = ["-nodisplay", "-nosplash", "-nodesktop", "-nojvm"]  # the framework's launch
PHASES = ("start", "setup", "startup", "task", "teardown", "overhead", "total")
DEFAULT_RATIO = 1.05  # PLAN: "no slower" is a median ratio of at most 1.05
MIN_REPEATS = 5

# The first lines of every Python variant: a mark before any import, then the script, module
# or code, then the last mark (also when it raises or exits).
BOOT = """\
import os as _o, time as _t
_f = _o.environ.get("PK5D_BENCH_FILE")
def _m(label):
    if _f:
        with open(_f, "a") as h:
            h.write('{"mark": "%s", "t": %.6f}\\n' % (label, _t.time()))
_m("first")
import sys as _s
_k, _x = _s.argv[1], _s.argv[2]
try:
    if _k == "code":
        _s.argv = ["-c"] + _s.argv[3:]
        exec(compile(_x, "<bench>", "exec"), {"__name__": "__main__"})
    else:
        import runpy as _r
        _s.argv = [_x] + _s.argv[3:]
        if _k == "module":
            _r.run_module(_x, run_name="__main__", alter_sys=True)
        else:
            _s.path[0] = _o.path.dirname(_o.path.abspath(_x))
            _r.run_path(_x, run_name="__main__")
finally:
    _m("last")
"""


# ---------------------------------------------------------------------------------------------
# inside a benchmarked process


def _append(line: str) -> None:
    f = os.environ.get("PK5D_BENCH_FILE")
    if f:
        with open(f, "a") as h:
            h.write(line + "\n")


def mark(label: str) -> None:
    """Leave a time mark (``ready``, or a mark of your own between ``ready`` and the end)."""
    _append(f'{{"mark": "{label}", "t": {time.time():.6f}}}')


def record(op: str, times, threads: int | None = None, **extra) -> None:
    """Leave the timings (seconds per call) of an operation measured inside the process."""
    import json

    _append(json.dumps({"op": op, "threads": threads, "times": [float(t) for t in times], **extra}))


def info(**values) -> None:
    """Leave facts about the process (versions, settings) for the result."""
    import json

    _append(json.dumps({"info": values}))


def timeit(fn, repeats: int = 5, warmup: int = 1, before=None) -> list[float]:
    """Seconds per call of ``fn()`` for ``repeats`` calls after ``warmup`` calls; ``before()``
    runs untimed before each call."""
    out = []
    for i in range(warmup + repeats):
        if before is not None:
            before()
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        if i >= warmup:
            out.append(dt)
    return out


# ---------------------------------------------------------------------------------------------
# statistics


def _median(v):
    import statistics

    return statistics.median(v) if v else None


def stats(values) -> dict | None:
    """Median, quartiles, range and count."""
    v = sorted(x for x in values if x is not None)
    if not v:
        return None
    n = len(v)

    def q(p):  # linear interpolation between order statistics (as MATLAB's quantile, type 5)
        h = n * p - 0.5
        if h <= 0:
            return v[0]
        if h >= n - 1:
            return v[-1]
        i = int(h)
        return v[i] + (h - i) * (v[i + 1] - v[i])

    return {"median": _median(v), "q1": q(0.25), "q3": q(0.75), "min": v[0], "max": v[-1], "n": n}


def phases(run: dict) -> dict:
    """The phase durations of one run from its marks (None where a mark is missing)."""
    m = run.get("marks", {})
    launch, ex = run.get("launch"), run.get("exit")
    first, ready, last = m.get("first"), m.get("ready"), m.get("last")
    begin = ready if ready is not None else first

    def d(a, b):
        return None if a is None or b is None else b - a

    p = {
        "start": d(launch, first),
        "setup": d(first, ready),
        "startup": d(launch, begin),
        "task": d(begin, last),
        "teardown": d(last, ex),
        "total": d(launch, ex),
    }
    p["overhead"] = (
        None if p["startup"] is None or p["teardown"] is None else (p["startup"] + p["teardown"])
    )
    # the steps between marks of the process's own (in time order)
    seq = sorted((t, k) for k, t in m.items() if t is not None)
    for (t0, a), (t1, b) in zip(seq, seq[1:], strict=False):
        if not (a in ("first", "ready") and b in ("ready", "last")):
            p[f"seg:{a}->{b}"] = t1 - t0
    return p


def _ok_runs(result: dict, variant: str, warmup: bool = False) -> list[dict]:
    return [
        r
        for r in result["runs"]
        if r["variant"] == variant and r["rc"] == 0 and bool(r["warmup"]) == warmup
    ]


def summarize(result: dict) -> dict:
    """Per variant: phase statistics, peak memory, CPU time, the first run (warm-up) apart, and
    per operation and thread count the timings inside the processes (all calls pooled; the
    spread is the range of the per-round medians)."""
    out = {}
    for v in result["spec"]["variants"]:
        vid = v["id"]
        runs = _ok_runs(result, vid)
        warm = [r for r in _ok_runs(result, vid, warmup=True) if not r.get("rewarm")]
        ph = [phases(r) for r in runs]
        keys = sorted({k for p in ph for k in p}, key=lambda k: (k not in PHASES, k))
        s = {
            "n": len(runs),
            "failed": sum(1 for r in result["runs"] if r["variant"] == vid and r["rc"] != 0),
            "phases": {k: stats([p.get(k) for p in ph]) for k in keys},
            "first": [phases(r) for r in warm],
            "first_in_job": [
                phases(r)
                for r in result["runs"]
                if r.get("first_in_job") and r["variant"] == vid and r["rc"] == 0
            ],
            "maxrss_mb": stats([r.get("maxrss_kb", 0) / 1024 for r in runs if r.get("maxrss_kb")]),
            "cpu_s": stats([r["user_s"] + r["sys_s"] for r in runs if "user_s" in r]),
            "ops": {},
        }
        pooled, rounds = {}, {}
        for r in runs:
            for rec in r.get("records", []):
                key = _opkey(rec["op"], rec.get("threads"))
                pooled.setdefault(key, []).extend(rec["times"])
                if rec["times"]:
                    rounds.setdefault(key, []).append(_median(rec["times"]))
        for key in sorted(pooled):
            st = stats(pooled[key]) or {}
            rm = rounds.get(key, [])
            st["round_min"], st["round_max"] = (min(rm), max(rm)) if rm else (None, None)
            st["round_medians"] = rm
            s["ops"][key] = st
        out[vid] = s
    return out


def _opkey(op: str, threads) -> str:
    return f"{op}@{threads}" if threads is not None else op


def _split_opkey(key: str) -> tuple[str, int | None]:
    op, _, t = key.rpartition("@")
    if op and t.isdigit():
        return op, int(t)
    return key, None


def _pairs(result: dict) -> list[dict]:
    spec = result["spec"]
    if spec.get("compare") is not None:
        return [c if isinstance(c, dict) else {"a": c[0], "b": c[1]} for c in spec["compare"]]
    ml = [v["id"] for v in spec["variants"] if v["impl"] == "matlab"]
    if len(ml) != 1:
        return []
    return [{"a": v["id"], "b": ml[0]} for v in spec["variants"] if v["impl"] == "python"]


def compare(result: dict, summary: dict | None = None) -> list[dict]:
    """Candidate (a) against reference (b): ratio a/b of the medians of each phase and of each
    operation both have (an ``op_map`` maps a's op names to b's), the median of the round by
    round ratios (the runs are paired by round), and whether the round medians overlap."""
    summary = summary or summarize(result)
    out = []
    for c in _pairs(result):
        a, b = summary.get(c["a"]), summary.get(c["b"])
        if a is None or b is None:
            continue
        rows = []
        for k in PHASES:
            sa, sb = a["phases"].get(k), b["phases"].get(k)
            if sa and sb and sb["median"]:
                rows.append(
                    {
                        "metric": k,
                        "a": sa["median"],
                        "b": sb["median"],
                        "ratio": sa["median"] / sb["median"],
                        "diff": sa["median"] - sb["median"],
                        "overlap": sa["min"] <= sb["max"] and sb["min"] <= sa["max"],
                    }
                )
        if a["maxrss_mb"] and b["maxrss_mb"]:
            rows.append(
                {
                    "metric": "maxrss",
                    "a": a["maxrss_mb"]["median"],
                    "b": b["maxrss_mb"]["median"],
                    "ratio": a["maxrss_mb"]["median"] / b["maxrss_mb"]["median"],
                }
            )
        op_map = c.get("op_map", {})
        for key, sa in a["ops"].items():
            op, t = _split_opkey(key)
            names = op_map.get(op, op)
            for bop in [names] if isinstance(names, str) else names:
                sb = b["ops"].get(_opkey(bop, t))
                if not sb or not sb.get("median"):
                    continue
                pr = [
                    x / y
                    for x, y in zip(sa["round_medians"], sb["round_medians"], strict=False)
                    if y
                ]
                rows.append(
                    {
                        "metric": "op:" + key,
                        "op": op,
                        "ref_op": bop,
                        "threads": t,
                        "a": sa["median"],
                        "b": sb["median"],
                        "ratio": sa["median"] / sb["median"],
                        "paired": _median(pr),
                        "a_range": [sa["round_min"], sa["round_max"]],
                        "b_range": [sb["round_min"], sb["round_max"]],
                        "overlap": sa["round_min"] <= sb["round_max"]
                        and sb["round_min"] <= sa["round_max"],
                    }
                )
        out.append({"a": c["a"], "b": c["b"], "rows": rows})
    return out


# ---------------------------------------------------------------------------------------------
# gates


def parse_gate(text: str, gpu: bool = False) -> dict:
    """The machine-checkable parts of a ``benchmark_gate`` text of tasks.json: the median ratio
    (an explicit ``(1.05x)``; "no slower" is 1.05 as PLAN.md defines it), the thread counts
    ("at 1, 4 and 24 threads"), a start-up + teardown limit in seconds, a peak-memory ratio.
    ``manual`` keeps the text, which a reviewer reads as well."""
    import re

    t = (text or "").strip()
    g: dict = {"text": t}
    if not t or t.lower().startswith("none"):
        g["none"] = True
        return g
    m = re.search(r"\((\d+(?:\.\d+)?)\s*x\)", t) or re.search(r"\b(\d\.\d+)\s*x\b", t)
    if m:
        g["ratio_max"] = float(m.group(1))
    elif re.search(r"no slower|not slower", t, re.I):
        g["ratio_max"] = DEFAULT_RATIO
    m = re.search(r"\bat ((?:\d+(?:\s*,\s*|\s*,?\s*and\s*|\s+))+)threads", t)
    if m:
        g["threads"] = [int(x) for x in re.findall(r"\d+", m.group(1))]
    m = re.search(
        r"start-?up(?: \+ teardown| and teardown)?[^.;]*?at most (\d+(?:\.\d+)?)\s*s\b", t, re.I
    )
    if m:
        g["overhead_max"] = float(m.group(1))
    m = re.search(r"memory[^.;]*?(\d+(?:\.\d+)?)\s*x", t, re.I)
    if m:
        g["mem_ratio_max"] = float(m.group(1))
    g["manual"] = True  # the text may say more than these parts (sizes, data, modes)
    return g


def task_gate(task: str, tasks_json: str | None = None) -> dict | None:
    """The gate of a task from tasks.json (``PK5D_TASKS_JSON``), parsed."""
    import json

    path = tasks_json or os.environ.get("PK5D_TASKS_JSON")
    if not path or not os.path.exists(path):
        return None
    with open(path) as f:
        d = json.load(f)
    for t in d.get("tasks", d if isinstance(d, list) else []):
        if t.get("id") == task:
            return parse_gate(t.get("benchmark_gate", ""))
    return None


def evaluate_gate(result: dict, gate: dict | None = None) -> dict:
    """Check the comparisons against a gate: ``ratio_max`` on the compared operations (or on
    ``metrics``), ``threads`` that must all be measured, ``overhead_max`` (seconds) on the
    candidate's start-up + teardown, ``mem_ratio_max`` on peak memory, ``within`` (|ratio - 1|
    at most this, for a build against itself). ``met`` is True, False, or None (nothing
    machine-checkable)."""
    gate = dict(gate if gate is not None else result.get("spec", {}).get("gate") or {})
    if gate.get("none"):
        return {"gate": gate, "met": None, "checks": [], "note": "no gate"}
    comps = result.get("comparisons") or compare(result)
    summ = result.get("summary") or summarize(result)
    checks = []
    metrics = gate.get("metrics")
    for c in comps:
        ops = [r for r in c["rows"] if r["metric"].startswith("op:")]
        if metrics is None:
            sel = ops
        else:
            sel = [r for r in c["rows"] if _metric_selected(r["metric"], metrics)]
        if "ratio_max" in gate:
            for r in sel:
                checks.append(
                    {
                        "pair": f"{c['a']}/{c['b']}",
                        "metric": r["metric"],
                        "value": r["ratio"],
                        "limit": gate["ratio_max"],
                        "ok": r["ratio"] <= gate["ratio_max"],
                    }
                )
        if "within" in gate:
            for r in sel or [r for r in c["rows"] if r["metric"] != "maxrss"]:
                checks.append(
                    {
                        "pair": f"{c['a']}/{c['b']}",
                        "metric": r["metric"],
                        "value": r["ratio"],
                        "limit": [1 - gate["within"], 1 + gate["within"]],
                        "ok": abs(r["ratio"] - 1) <= gate["within"],
                    }
                )
        if "threads" in gate:
            for opn in sorted({r["op"] for r in ops}):
                have = {r["threads"] for r in ops if r["op"] == opn}
                missing = [t for t in gate["threads"] if t not in have]
                checks.append(
                    {
                        "pair": f"{c['a']}/{c['b']}",
                        "metric": f"threads:{opn}",
                        "value": sorted(t for t in have if t is not None),
                        "limit": gate["threads"],
                        "ok": not missing,
                    }
                )
        if "mem_ratio_max" in gate:
            for r in c["rows"]:
                if r["metric"] == "maxrss":
                    checks.append(
                        {
                            "pair": f"{c['a']}/{c['b']}",
                            "metric": "maxrss",
                            "value": r["ratio"],
                            "limit": gate["mem_ratio_max"],
                            "ok": r["ratio"] <= gate["mem_ratio_max"],
                        }
                    )
    if "overhead_max" in gate:
        cands = {c["a"] for c in comps} or {
            v["id"] for v in result["spec"]["variants"] if v["impl"] == "python"
        }
        for vid in sorted(cands):
            st = summ.get(vid, {}).get("phases", {}).get("overhead")
            if st:
                checks.append(
                    {
                        "pair": vid,
                        "metric": "overhead",
                        "value": st["median"],
                        "limit": gate["overhead_max"],
                        "ok": st["median"] <= gate["overhead_max"],
                    }
                )
    met = all(c["ok"] for c in checks) if checks else None
    return {"gate": gate, "met": met, "checks": checks}


def _metric_selected(metric: str, metrics) -> bool:
    import fnmatch

    return any(fnmatch.fnmatchcase(metric, m) for m in metrics)


# ---------------------------------------------------------------------------------------------
# node, versions, source


def _read(path: str) -> str | None:
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return None


def node_info() -> dict:
    """Host, CPU model, CPUs of the process, memory, load, huge pages and NUMA free memory,
    the Slurm job, GPUs (nvidia-smi, when the process may see one)."""
    import glob
    import platform
    import socket

    d: dict = {
        "host": socket.gethostname(),
        "system": platform.system(),
        "kernel": platform.release(),
        "machine": platform.machine(),
    }
    cpu = _read("/proc/cpuinfo") or ""
    d["cpu_model"] = next(
        (ln.split(":", 1)[1].strip() for ln in cpu.splitlines() if ln.startswith("model name")),
        platform.processor(),
    )
    d["cpus_online"] = os.cpu_count()
    try:
        d["cpus_affinity"] = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        d["cpus_affinity"] = None
    mem = _read("/proc/meminfo") or ""
    for ln in mem.splitlines():
        k, _, v = ln.partition(":")
        if k in ("MemTotal", "MemAvailable", "HugePages_Free", "HugePages_Total"):
            d[k] = int(v.split()[0])
    d["loadavg"] = (_read("/proc/loadavg") or "").split()[:3]
    thp = _read("/sys/kernel/mm/transparent_hugepage/enabled")
    d["thp"] = thp.strip() if thp else None
    numa = {}
    for f in sorted(glob.glob("/sys/devices/system/node/node*/meminfo")):
        txt = _read(f) or ""
        free = [ln.split()[-2] for ln in txt.splitlines() if "MemFree" in ln]
        if free:
            numa[os.path.basename(os.path.dirname(f))] = int(free[0])
    d["numa_free_kb"] = numa
    d["slurm"] = {
        k: os.environ[k]
        for k in (
            "SLURM_JOB_ID",
            "SLURM_JOB_NAME",
            "SLURM_JOB_PARTITION",
            "SLURM_CPUS_PER_TASK",
            "SLURM_MEM_PER_NODE",
            "SLURMD_NODENAME",
            "SLURM_JOB_GPUS",
        )
        if k in os.environ
    }
    d["node_type"] = os.environ.get("PK5D_BENCH_NODETYPE")
    d["cuda_visible_devices"] = os.environ.get("CUDA_VISIBLE_DEVICES")
    if os.environ.get("CUDA_VISIBLE_DEVICES", None) != "":
        d["gpus"] = _nvidia_smi()
    return d


def _nvidia_smi() -> list[dict] | None:
    import shutil
    import subprocess

    exe = shutil.which("nvidia-smi")
    if not exe:
        return None
    q = "index,name,memory.used,memory.total,utilization.gpu,driver_version,persistence_mode"
    try:
        out = subprocess.run(
            [exe, f"--query-gpu={q}", "--format=csv,noheader,nounits"],
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    keys = q.split(",")
    return [
        dict(zip(keys, (x.strip() for x in ln.split(",")), strict=False))
        for ln in out.splitlines()
        if ln.strip()
    ]


def _busy_processes(limit: int = 5) -> list[dict]:
    """The processes using the most CPU on the node, as own or other (no user names)."""
    import subprocess

    try:
        out = subprocess.run(
            ["ps", "-eo", "uid,pcpu,rss,comm", "--sort=-pcpu"],
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.splitlines()[1:]
    except (OSError, subprocess.SubprocessError):
        return []
    me = os.getuid() if hasattr(os, "getuid") else None
    rows = []
    for ln in out[:limit]:
        p = ln.split(None, 3)
        if len(p) == 4:
            rows.append(
                {
                    "own": me is not None and int(p[0]) == me,
                    "pcpu": float(p[1]),
                    "rss_mb": int(p[2]) // 1024,
                    "comm": p[3],
                }
            )
    return rows


def versions() -> dict:
    """Python and the installed packages (without importing them)."""
    from importlib import metadata

    d = {"python": sys.version.split()[0], "executable": sys.executable, "prefix": sys.prefix}
    for name in (
        "numpy",
        "scipy",
        "scikit-image",
        "cpp-zarr",
        "cpp-tiff",
        "h5py",
        "cupy-cuda12x",
        "nvidia-cufft-cu12",
    ):
        try:
            d[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pass
    try:
        from petakit5d import __version__

        d["petakit5d"] = __version__
    except ImportError:
        pass
    d["freeze"] = sorted(
        f"{x.metadata['Name']}=={x.version}"
        for x in metadata.distributions()
        if x.metadata and x.metadata["Name"]
    )
    d["PETAKIT5D_EXT_DIR"] = os.environ.get("PETAKIT5D_EXT_DIR")
    return d


def source_info() -> dict:
    """The commit and local changes of the benchmarked tree (set by the cluster helper)."""
    return {
        k[len("PK5D_BENCH_") :].lower(): os.environ[k]
        for k in ("PK5D_BENCH_COMMIT", "PK5D_BENCH_DIRTY", "PK5D_BENCH_WORKTREE", "PK5D_BENCH_TASK")
        if k in os.environ
    }


# ---------------------------------------------------------------------------------------------
# running a benchmark


def load_spec(path_or_spec) -> dict:
    """A spec from a dict, a JSON file or a Python file that defines ``SPEC``."""
    import json

    if isinstance(path_or_spec, dict):
        spec = dict(path_or_spec)
    elif str(path_or_spec).endswith(".py"):
        import runpy

        spec = dict(runpy.run_path(str(path_or_spec))["SPEC"])
    else:
        with open(path_or_spec) as f:
            spec = json.load(f)
    if not spec.get("name") or not spec.get("variants"):
        raise ValueError("a benchmark spec needs a name and variants")
    ids = [v.get("id") for v in spec["variants"]]
    if len(set(ids)) != len(ids) or not all(ids):
        raise ValueError(f"variant ids must be unique and given: {ids}")
    for v in spec["variants"]:
        if v.get("impl") not in ("python", "matlab", "command"):
            raise ValueError(f"variant {v['id']}: impl must be python, matlab or command")
    return spec


def _fmt(s, ctx: dict):
    """Placeholders ``{name}`` of ``ctx`` and ``{env[NAME]}`` of the environment; other braces
    (MATLAB cell indexing, Python dicts) stay as they are."""
    if isinstance(s, str):
        import re

        def sub(m):
            k = m.group(1)
            if k.startswith("env["):
                name = k[4:-1]
                if name not in os.environ:
                    raise ValueError(f"environment variable {name} not set for {s!r}")
                return os.environ[name]
            return ctx[k] if k in ctx else m.group(0)

        return re.sub(r"\{(env\[[A-Za-z_][A-Za-z0-9_]*\]|[A-Za-z_][A-Za-z0-9_]*)\}", sub, s)
    if isinstance(s, list):
        return [_fmt(x, ctx) for x in s]
    if isinstance(s, dict):
        return {k: _fmt(x, ctx) for k, x in s.items()}
    return s


def _matlab_argv(v: dict, ctx: dict) -> list[str]:
    common = os.path.join(ctx["src"], "benchmarks", "common")
    paths = [common] + [
        p if os.path.isabs(p) else os.path.join(ctx["src"], p)
        for p in _fmt(v.get("paths", []), ctx)
    ]
    q = "', '".join(p.replace("'", "''") for p in paths)
    code = _fmt(v.get("code", ""), ctx).strip()
    if code and not code.endswith((";", ",")):
        code += ";"
    setup = "setup([]); pk5d_bench_mark('setup'); " if v.get("setup", True) else ""
    stmt = (
        f"addpath('{common}'); pk5d_bench_mark('first'); {setup}addpath('{q}'); "
        f"pk5d_bench_info(); pk5d_bench_mark('ready'); {code} pk5d_bench_mark('last');"
    )
    exe = v.get("matlab") or os.environ.get("PETAKIT5D_MATLAB", "matlab")
    flags = v.get("flags", MATLAB_FLAGS)
    if v.get("launch", "batch") == "r":
        # the framework's launch before -batch: -r with an explicit exit
        stmt = f"try, {stmt} catch ME, disp(getReport(ME)), exit(1), end, exit(0)"
        return [exe, *flags, "-r", stmt]
    return [exe, *flags, "-batch", stmt]


def _python_argv(v: dict, ctx: dict, python: str) -> list[str]:
    args = [str(a) for a in _fmt(v.get("args", []), ctx)]
    pyflags = v.get("python_flags", [])
    if "code" in v:
        return [python, *pyflags, "-c", BOOT, "code", _fmt(v["code"], ctx), *args]
    if "module" in v:
        return [python, *pyflags, "-c", BOOT, "module", v["module"], *args]
    script = _fmt(v["script"], ctx)
    if not os.path.isabs(script):
        script = os.path.join(ctx["src"], script)
    return [python, *pyflags, "-c", BOOT, "script", script, *args]


def evict(paths) -> dict:
    """Drop the page-cache pages of every file under ``paths`` (posix_fadvise DONTNEED; the
    kernel's file metadata stays cached, and on network storage so does the server's cache).
    Returns the number of files and the seconds it took."""
    t0 = time.perf_counter()
    n = 0
    adv = getattr(os, "posix_fadvise", None)
    if adv is None:
        return {"files": 0, "seconds": 0.0, "supported": False}
    seen = set()
    for top in paths:
        for root, _dirs, files in os.walk(top):
            for name in files:
                p = os.path.join(root, name)
                try:
                    st = os.stat(p)
                except OSError:
                    continue
                if (st.st_dev, st.st_ino) in seen or not os.path.isfile(p):
                    continue
                seen.add((st.st_dev, st.st_ino))
                try:
                    fd = os.open(p, os.O_RDONLY)
                except OSError:
                    continue
                try:
                    adv(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                    n += 1
                finally:
                    os.close(fd)
    return {"files": n, "seconds": time.perf_counter() - t0, "supported": True}


def local_python_copy(dest: str) -> dict:
    """Copy this interpreter's installation, its virtual environment, the source tree on
    PYTHONPATH and the extension folder to ``dest`` (node-local disk), with the environment
    pointed at the copied interpreter. Returns the python executable and the environment
    variables for it."""
    import shutil

    os.makedirs(dest, exist_ok=True)
    base = os.path.realpath(sys.base_prefix)
    out = {"env": {}}
    nbase = os.path.join(dest, "base")
    shutil.copytree(base, nbase, symlinks=True, dirs_exist_ok=True)
    if sys.prefix != sys.base_prefix:
        venv = os.path.realpath(sys.prefix)
        nvenv = os.path.join(dest, "venv")
        shutil.copytree(venv, nvenv, symlinks=True, dirs_exist_ok=True)
        cfg = os.path.join(nvenv, "pyvenv.cfg")
        with open(cfg) as f:
            lines = f.read().splitlines()
        with open(cfg, "w") as f:
            for ln in lines:
                k = ln.split("=", 1)[0].strip()
                if k == "home":
                    ln = "home = " + os.path.join(nbase, "bin")
                elif k in ("executable", "command"):
                    continue
                f.write(ln.replace(base, nbase) + "\n")
        bindir = os.path.join(nvenv, "bin")
        for name in os.listdir(bindir):
            p = os.path.join(bindir, name)
            if os.path.islink(p) and os.path.realpath(p).startswith(base):
                t = os.path.join(nbase, os.path.relpath(os.path.realpath(p), base))
                os.unlink(p)
                os.symlink(t, p)
        out["python"] = os.path.join(bindir, os.path.basename(sys.executable))
    else:
        out["python"] = os.path.join(nbase, os.path.relpath(os.path.realpath(sys.executable), base))
    pp = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]
    npp = []
    for i, p in enumerate(pp):
        q = os.path.join(dest, f"path{i}")
        shutil.copytree(
            p, q, symlinks=True, dirs_exist_ok=True, ignore=shutil.ignore_patterns("__pycache__")
        )
        npp.append(q)
    if npp:
        out["env"]["PYTHONPATH"] = os.pathsep.join(npp)
    ext = os.environ.get("PETAKIT5D_EXT_DIR")
    if ext and os.path.isdir(ext):
        q = os.path.join(dest, "ext")
        shutil.copytree(ext, q, symlinks=True, dirs_exist_ok=True)
        out["env"]["PETAKIT5D_EXT_DIR"] = q
    return out


def _read_marks(path: str) -> tuple[dict, list, dict]:
    import json

    marks, records, inf = {}, [], {}
    try:
        with open(path) as f:
            lines = f.read().splitlines()
    except OSError:
        return marks, records, inf
    for ln in lines:
        try:
            d = json.loads(ln)
        except ValueError:
            continue
        if "mark" in d:
            # first and ready: the first one; last: the last one; own marks: the last one
            if d["mark"] in ("first", "ready") and d["mark"] in marks:
                continue
            marks[d["mark"]] = d["t"]
        elif "op" in d:
            records.append(d)
        elif "info" in d:
            inf.update(d["info"])
    return marks, records, inf


def _stop_group(p, timed_out: list, done) -> None:
    """Stop a benchmarked process that ran past its timeout, with the processes it started
    (its own process group): SIGTERM, then SIGKILL after 10 s if it has not been reaped."""
    import signal

    timed_out.append(time.time())
    for sig, wait in ((signal.SIGTERM, 10.0), (getattr(signal, "SIGKILL", signal.SIGTERM), 0)):
        if done.is_set():
            return
        try:
            if hasattr(os, "killpg"):
                os.killpg(p.pid, sig)
            else:
                p.kill()
        except OSError:
            return
        done.wait(wait)


def _run_process(argv, cwd, env, marks_file, log_file, timeout: float | None = None) -> dict:
    import subprocess
    import threading

    timer, timed_out, done = None, [], threading.Event()
    with open(log_file, "ab") as lf:
        t0 = time.time()
        p = subprocess.Popen(
            argv,
            cwd=cwd,
            env=env,
            stdout=lf,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=bool(timeout) and os.name == "posix",
        )
        if timeout:
            timer = threading.Timer(timeout, _stop_group, (p, timed_out, done))
            timer.daemon = True
            timer.start()
        if hasattr(os, "wait4"):
            _, status, ru = os.wait4(p.pid, 0)
            t1 = time.time()
            p.returncode = os.waitstatus_to_exitcode(status)
            usage = {
                "maxrss_kb": ru.ru_maxrss,
                "user_s": ru.ru_utime,
                "sys_s": ru.ru_stime,
                "majflt": ru.ru_majflt,
                "minflt": ru.ru_minflt,
            }
        else:
            p.wait()
            t1 = time.time()
            usage = {}
    done.set()
    if timer is not None:
        timer.cancel()
    marks, records, inf = _read_marks(marks_file)
    if timed_out:
        usage["timed_out"] = timed_out[0] - t0
    return {
        "launch": t0,
        "exit": t1,
        "rc": p.returncode if not timed_out else (p.returncode or -1),
        "marks": marks,
        "records": records,
        "info": inf,
        **usage,
    }


def run(
    spec,
    out_dir: str | None = None,
    repeats: int | None = None,
    keep: bool = False,
    tasks_json: str | None = None,
    log=print,
) -> dict:
    """Run a benchmark spec (see the module doc) and write its result JSON to ``out_dir``."""
    import shutil
    import subprocess
    import tempfile

    spec = load_spec(spec)
    repeats = int(repeats or spec.get("repeats", MIN_REPEATS))
    if repeats < MIN_REPEATS and not spec.get("allow_few"):
        raise ValueError(f"at least {MIN_REPEATS} repeats (allow_few in the spec for a test)")
    warmup = int(spec.get("warmup", 1))
    work = tempfile.mkdtemp(prefix="pk5d_bench_")
    src = os.environ.get("PK5D_SRC") or os.getcwd()
    ctx = {
        "src": src,
        "work": work,
        "python": sys.executable,
        "oracle": os.environ.get("PETAKIT5D_ORACLE_DIR", ""),
        "data": os.environ.get("PETAKIT5D_DEMO_DATA", ""),
        "common": os.path.join(src, "benchmarks", "common"),
    }
    ctx.update({k: str(x) for k, x in spec.get("params", {}).items()})
    started = time.strftime("%Y%m%d_%H%M%S")
    result: dict = {
        "schema": SCHEMA,
        "name": spec["name"],
        "task": spec.get("task") or os.environ.get("PK5D_BENCH_TASK"),
        "started": started,
        "spec": spec,
        "node": node_info(),
        "versions": versions(),
        "source": source_info(),
        "busy_before": _busy_processes(),
        "setup": [],
        "runs": [],
        "skipped": [],
    }
    log(
        f"bench {spec['name']}: {result['node']['host']} ({result['node']['cpu_model']}, "
        f"{result['node']['cpus_affinity']} CPUs), {repeats} rounds, work {work}"
    )
    try:
        for cmd in spec.get("setup", []):
            c = _fmt(cmd, ctx)
            t0 = time.perf_counter()
            p = subprocess.run(c, shell=True, cwd=work, capture_output=True, text=True)
            result["setup"].append(
                {
                    "cmd": c,
                    "rc": p.returncode,
                    "seconds": time.perf_counter() - t0,
                    "output": (p.stdout + p.stderr)[-4000:],
                }
            )
            if p.returncode:
                raise RuntimeError(f"setup failed ({p.returncode}): {c}\n{p.stdout}{p.stderr}")
        variants = []
        for v in spec["variants"]:
            need = v.get("when")
            if need and not os.environ.get(need):
                result["skipped"].append({"id": v["id"], "why": f"{need} not set"})
                continue
            variants.append(v)
        local = None
        prepared = {}
        for v in variants:
            env = os.environ.copy()
            env.update({k: str(x) for k, x in _fmt(v.get("env", {}), ctx).items()})
            cwd = _fmt(v.get("cwd", "{oracle}" if v["impl"] == "matlab" else "{work}"), ctx)
            if v["impl"] == "python":
                python = _fmt(v.get("python", sys.executable), ctx)
                if v.get("env_copy") == "local":
                    if local is None:
                        t0 = time.perf_counter()
                        local = local_python_copy(os.path.join(work, "localpy"))
                        result["setup"].append(
                            {
                                "cmd": "local python copy",
                                "rc": 0,
                                "seconds": time.perf_counter() - t0,
                                "output": local["python"],
                            }
                        )
                    python = local["python"]
                    env.update(local["env"])
                argv = _python_argv(v, ctx, python)
                ev = v.get("evict")
                if ev is None and v.get("cache") == "cold":
                    ev = [os.path.realpath(sys.base_prefix), os.path.realpath(sys.prefix)]
                    ev += [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p]
                    if env.get("PETAKIT5D_EXT_DIR"):
                        ev.append(env["PETAKIT5D_EXT_DIR"])
                    if v.get("env_copy") == "local":
                        ev = [os.path.join(work, "localpy")]
            elif v["impl"] == "matlab":
                argv = _matlab_argv(v, ctx)
                ev = v.get("evict")
            else:
                argv = [str(a) for a in _fmt(v["argv"], ctx)]
                ev = v.get("evict")
            prepared[v["id"]] = (argv, cwd, env, _fmt(ev, ctx) if ev else None)
        seq = [(v, -1 - w) for w in range(warmup) for v in variants]
        order = spec.get("order", "alternate")
        for r in range(repeats):
            vs = list(variants)
            if order == "rotate":
                vs = vs[r % len(vs) :] + vs[: r % len(vs)]
            elif r % 2:
                vs.reverse()
            seq += [(v, r) for v in vs]
        seen_matlab = False
        stale: set = set()  # warm variants whose files a cold variant evicted since their run
        k = 0
        for v, r in seq:
            argv, cwd, env, ev = prepared[v["id"]]
            rewarm = not ev and v["id"] in stale
            for again in ([True] if rewarm else []) + [False]:
                e = dict(env)
                mf = os.path.join(work, f"marks_{k:04d}_{v['id']}.jsonl")
                lf = os.path.join(work, f"log_{k:04d}_{v['id']}.txt")
                e["PK5D_BENCH_FILE"] = mf
                e["PK5D_BENCH_ROUND"] = str(r)
                ev_info = evict(ev) if ev else None
                if ev:
                    stale = {
                        x["id"] for x in variants if x["id"] != v["id"] and x["impl"] == "python"
                    }
                rec = _run_process(argv, cwd, e, mf, lf, timeout=v.get("timeout"))
                rec.update(
                    {
                        "variant": v["id"],
                        "impl": v["impl"],
                        "round": r,
                        "seq": k,
                        "warmup": r < 0 or again,
                        "rewarm": again,
                        "evicted": ev_info,
                    }
                )
                k += 1
                if v["impl"] == "matlab" and not seen_matlab:
                    rec["first_in_job"] = True
                    seen_matlab = True
                if rec["rc"]:
                    with open(lf, errors="replace") as f:
                        rec["log_tail"] = f.read()[-3000:]
                    log(f"  {v['id']} round {r}: exit {rec['rc']}\n{rec['log_tail'][-1500:]}")
                result["runs"].append(rec)
                ph = phases(rec)
                what = "re-warm" if again else ("warm-up" if r < 0 else f"round {r}")
                log(
                    f"  {v['id']:<18} {what:<9} "
                    + " ".join(
                        f"{p} {ph[p]:.3f}"
                        for p in ("startup", "task", "teardown")
                        if ph.get(p) is not None
                    )
                )
            stale.discard(v["id"])
        result["busy_after"] = _busy_processes()
        result["node_after"] = {"loadavg": node_info()["loadavg"]}
        if result["node"].get("gpus") is not None:
            result["node_after"]["gpus"] = _nvidia_smi()
        result["summary"] = summarize(result)
        result["comparisons"] = compare(result, result["summary"])
        gate = spec.get("gate")
        if gate is None and result["task"]:
            gate = task_gate(result["task"], tasks_json)
        result["gate"] = evaluate_gate(result, gate) if gate else None
        result["finished"] = time.strftime("%Y%m%d_%H%M%S")
        if out_dir:
            result["path"] = _write_result(result, out_dir, work)
            log(f"result {result['path']}")
    finally:
        if not keep:
            shutil.rmtree(work, ignore_errors=True)
    return result


def _write_result(result: dict, out_dir: str, work: str) -> str:
    import json
    import shutil

    os.makedirs(out_dir, exist_ok=True)
    job = result["node"]["slurm"].get("SLURM_JOB_ID", str(os.getpid()))
    stem = f"{result['started']}_{result['name']}_{job}"
    logs = os.path.join(out_dir, stem)
    os.makedirs(logs, exist_ok=True)
    for f in os.listdir(work):
        if f.startswith(("log_", "marks_")):
            shutil.copy2(os.path.join(work, f), logs)
    path = os.path.join(out_dir, stem + ".json")
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(result, f, indent=1, default=str)
    os.replace(tmp, path)
    return path


# ---------------------------------------------------------------------------------------------
# reports


def _s(x, digits: int = 3) -> str:
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.{digits}f}" if abs(x) < 100 else f"{x:.1f}"
    return str(x)


_VERDICT = {True: "met", False: "MISSED", None: "nothing to check"}


def report(result: dict) -> str:
    """Markdown tables: the node, each variant's phases (median [min-max]), the first run,
    the operations, the comparisons and the gate."""
    s = result.get("summary") or summarize(result)
    n = result["node"]
    dirty = result["source"].get("dirty")
    diff = f" + diff {dirty[:12]}" if dirty not in (None, "", "clean") else ""
    lines = [
        f"## {result['name']} ({result.get('task') or '-'}), {result['started']}",
        "",
        f"Node {n['host']} ({n.get('node_type') or '?'}), {n['cpu_model']}, "
        f"{n['cpus_affinity']} CPUs, job {n['slurm'].get('SLURM_JOB_ID', '-')}, load "
        f"{' '.join(n['loadavg'])}; source {result['source'].get('commit', '?')[:12]}"
        f"{diff}",
        "",
    ]
    if n.get("gpus"):
        mine = set(n["slurm"].get("SLURM_JOB_GPUS", "").split(","))
        lines += [
            "GPUs at the start (* = this job's, from SLURM_JOB_GPUS): "
            + "; ".join(
                f"{g.get('index')}{'*' if g.get('index') in mine else ''} {g.get('name')} "
                f"({g.get('memory.used')} MiB used, {g.get('utilization.gpu')}% busy)"
                for g in n["gpus"]
            ),
            "",
        ]
    for sk in result.get("skipped", []):
        lines.append(f"Skipped {sk['id']}: {sk['why']}")
    cols = ["start", "setup", "startup", "task", "teardown", "overhead", "total"]
    lines += [
        "| variant | n | " + " | ".join(cols) + " | first run: startup / teardown | "
        "peak MB | CPU s |",
        "|" + "---|" * (len(cols) + 5),
    ]
    for vid, v in s.items():
        cells = []
        for c in cols:
            st = v["phases"].get(c)
            cells.append(f"{_s(st['median'])} [{_s(st['min'])}-{_s(st['max'])}]" if st else "-")
        fr = v["first_in_job"] or v["first"]
        first = (
            (
                f"{_s(fr[0].get('startup'))} / {_s(fr[0].get('teardown'))}"
                + (" (first MATLAB of the job)" if v["first_in_job"] else "")
            )
            if fr
            else "-"
        )
        mem = _s(v["maxrss_mb"]["median"], 0) if v["maxrss_mb"] else "-"
        cpu = _s(v["cpu_s"]["median"], 2) if v["cpu_s"] else "-"
        fail = f" ({v['failed']} failed)" if v["failed"] else ""
        lines.append(
            f"| {vid} | {v['n']}{fail} | " + " | ".join(cells) + f" | {first} | {mem} | {cpu} |"
        )
    segs = sorted({k for v in s.values() for k in v["phases"] if k.startswith("seg:")})
    if segs:
        lines += [
            "",
            "Steps between the processes' own marks (median s):",
            "",
            "| variant | " + " | ".join(x[4:] for x in segs) + " |",
            "|" + "---|" * (len(segs) + 1),
        ]
        for vid, v in s.items():
            if any(k in v["phases"] for k in segs):
                lines.append(
                    f"| {vid} | "
                    + " | ".join(
                        _s(v["phases"][k]["median"]) if v["phases"].get(k) else "-" for k in segs
                    )
                    + " |"
                )
    if any(v["ops"] for v in s.values()):
        lines += [
            "",
            "| variant | op @ threads | median s | round medians | n |",
            "|---|---|---|---|---|",
        ]
        for vid, v in s.items():
            for key, st in v["ops"].items():
                lines.append(
                    f"| {vid} | {key} | {_s(st['median'], 4)} | "
                    f"{_s(st['round_min'], 4)}-{_s(st['round_max'], 4)} | {st['n']} |"
                )
    for c in result.get("comparisons") or []:
        lines += [
            "",
            f"{c['a']} / {c['b']}:",
            "",
            "| metric | a | b | ratio | paired | overlap |",
            "|---|---|---|---|---|---|",
        ]
        for r in c["rows"]:
            lines.append(
                f"| {r['metric']} | {_s(r['a'], 4)} | {_s(r['b'], 4)} | "
                f"**{r['ratio']:.3f}** | {_s(r.get('paired'))} | "
                f"{'' if r.get('overlap') is None else ('yes' if r['overlap'] else 'no')} |"
            )
    g = result.get("gate")
    if g:
        lines += [
            "",
            f"Gate: {g['gate'].get('text') or g['gate']}",
            f"Machine-checked: **{_VERDICT[g['met']]}**"
            + (" (the text may ask for more: read it)" if g["gate"].get("manual") else ""),
        ]
        for c in g["checks"]:
            if not c["ok"]:
                lines.append(
                    f"- missed: {c['pair']} {c['metric']} = {_s(c['value'])} (limit {c['limit']})"
                )
    return "\n".join(lines)


def _pct(x) -> str:
    return "" if x is None else f"{x:+.1%}"


def history(
    bench_dir: str, write: bool = False, threshold: float = 0.05, min_seconds: float = 0.005
) -> str:
    """A table of every result under ``bench_dir``: per benchmark and metric, each result's
    date, task, commit, node type and node, the candidate's median and its ratio to the
    reference, and the change from the previous result of the same benchmark, metric and node
    type; a ratio (or, without a reference, a median on the same node type) that grew by more
    than ``threshold``, with the median changed by more than ``min_seconds``, is flagged."""
    import glob
    import json

    rows = []
    for f in sorted(glob.glob(os.path.join(bench_dir, "**", "*.json"), recursive=True)):
        try:
            with open(f) as h:
                r = json.load(h)
        except (OSError, ValueError):
            continue
        if not isinstance(r, dict) or r.get("schema") != SCHEMA:
            continue
        base = {
            "bench": r["name"],
            "task": r.get("task") or "-",
            "date": r["started"],
            "commit": (r.get("source", {}).get("commit") or "?")[:9]
            + ("+" if r.get("source", {}).get("dirty") not in (None, "", "clean") else ""),
            "ntype": r["node"].get("node_type") or "?",
            "node": r["node"]["host"].split(".")[0],
            "file": os.path.relpath(f, bench_dir),
        }
        comps = r.get("comparisons") or []
        done = set()
        for c in comps:
            for x in c["rows"]:
                if x["metric"] == "maxrss":
                    continue
                rows.append(
                    {
                        **base,
                        "metric": f"{c['a']}:{x['metric']}",
                        "value": x["a"],
                        "ref": x["b"],
                        "ratio": x["ratio"],
                    }
                )
                done.add((c["a"], x["metric"]))
        for vid, v in (r.get("summary") or {}).items():
            for k in ("startup", "teardown", "task"):
                st = v["phases"].get(k)
                if st and (vid, k) not in done:
                    rows.append(
                        {
                            **base,
                            "metric": f"{vid}:{k}",
                            "value": st["median"],
                            "ref": None,
                            "ratio": None,
                        }
                    )
            for key, st in v["ops"].items():
                if (vid, "op:" + key) not in done:
                    rows.append(
                        {
                            **base,
                            "metric": f"{vid}:op:{key}",
                            "value": st["median"],
                            "ref": None,
                            "ratio": None,
                        }
                    )
    rows.sort(key=lambda x: (x["bench"], x["metric"], x["ntype"], x["date"]))
    prev: dict = {}
    for x in rows:
        k = (x["bench"], x["metric"], x["ntype"])
        p = prev.get(k)
        x["change"] = None
        if p is not None:
            if x["ratio"] is not None and p["ratio"]:
                x["change"] = x["ratio"] / p["ratio"] - 1
            elif p["value"]:
                x["change"] = x["value"] / p["value"] - 1
        # a change of a few milliseconds is noise, whatever its percentage
        big = p is not None and abs(x["value"] - p["value"]) > min_seconds
        x["flag"] = (
            "REGRESSION" if x["change"] is not None and x["change"] > threshold and big else ""
        )
        prev[k] = x
    out = [
        "# Benchmark history",
        "",
        f"Every result under this folder, per benchmark and metric (candidate variant:metric); "
        f"change = against the previous result of the same metric on the same node type "
        f"(the ratio to the reference when there is one); flagged above {threshold:.0%} when "
        f"the median moved by more than {min_seconds * 1000:g} ms.",
        "",
    ]
    flagged = [x for x in rows if x["flag"]]
    out += [f"Flagged: {len(flagged)}", ""]
    bench = None
    for x in rows:
        if x["bench"] != bench:
            bench = x["bench"]
            out += [
                "",
                f"## {bench}",
                "",
                "| metric | date | task | commit | node | median s | reference s | ratio | "
                "change | |",
                "|---|---|---|---|---|---|---|---|---|---|",
            ]
        out.append(
            f"| {x['metric']} | {x['date']} | {x['task']} | {x['commit']} | "
            f"{x['ntype']} {x['node']} | {_s(x['value'], 4)} | {_s(x['ref'], 4)} | "
            f"{_s(x['ratio'])} | "
            f"{_pct(x['change'])} | {x['flag']} |"
        )
    text = "\n".join(out) + "\n"
    if write:
        for name, content in (("HISTORY.md", text), ("history.json", json.dumps(rows, indent=1))):
            import socket

            p = os.path.join(bench_dir, name)
            tmp = f"{p}.tmp_{socket.gethostname()}_{os.getpid()}"
            with open(tmp, "w") as f:
                f.write(content)
            os.replace(tmp, p)
    return text


# ---------------------------------------------------------------------------------------------
# command line


def main(argv=None) -> int:
    import argparse
    import json

    ap = argparse.ArgumentParser(
        prog="python -m petakit5d.testing.bench",
        description="Interleaved MATLAB/Python benchmarks.",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="run a spec and write its result JSON")
    p.add_argument("spec")
    p.add_argument("--out", help="folder for the result (default: none, print only)")
    p.add_argument("--repeats", type=int)
    p.add_argument("--keep", action="store_true", help="keep the work folder")
    p.add_argument("--tasks", help="tasks.json for the gate (default PK5D_TASKS_JSON)")
    p.add_argument("--json", action="store_true", help="print the summary as JSON")
    p = sub.add_parser("summary", help="tables of a result")
    p.add_argument("result")
    p = sub.add_parser("gate", help="check a result against a gate")
    p.add_argument("result")
    p.add_argument("--task", help="the task whose gate applies (default: the result's)")
    p.add_argument("--tasks", help="tasks.json (default PK5D_TASKS_JSON)")
    p.add_argument("--text", help="a gate text instead of the task's")
    p = sub.add_parser("history", help="a table of all results in a folder")
    p.add_argument("bench_dir")
    p.add_argument("--write", action="store_true", help="write HISTORY.md and history.json")
    p.add_argument("--threshold", type=float, default=0.05)
    p.add_argument("--min-seconds", type=float, default=0.005)
    sub.add_parser("node", help="print the node, versions and source as JSON")
    a = ap.parse_args(argv)
    if a.cmd == "run":
        r = run(
            a.spec,
            out_dir=a.out,
            repeats=a.repeats,
            keep=a.keep,
            tasks_json=a.tasks,
            log=lambda s: print(s, flush=True),
        )
        print(json.dumps(r["summary"], indent=1, default=str) if a.json else report(r))
        bad = sum(1 for x in r["runs"] if x["rc"])
        return 1 if bad else 0
    if a.cmd in ("summary", "gate"):
        with open(a.result) as f:
            r = json.load(f)
        if a.cmd == "summary":
            print(report(r))
            return 0
        g = parse_gate(a.text) if a.text else task_gate(a.task or r.get("task"), a.tasks)
        if g is None:
            g = r.get("spec", {}).get("gate")
        ev = evaluate_gate(r, g)
        print(json.dumps(ev, indent=1, default=str))
        return {True: 0, False: 1, None: 0}[ev["met"]]
    if a.cmd == "history":
        print(history(a.bench_dir, write=a.write, threshold=a.threshold, min_seconds=a.min_seconds))
        return 0
    print(
        json.dumps({"node": node_info(), "versions": versions(), "source": source_info()}, indent=1)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
