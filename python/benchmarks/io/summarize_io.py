"""Summarize bench_io.py / bench_io_mex.m runs: one folder per storage, holding
python_r<k>.json, matlab_r<k>.json and process_times.jsonl (written by the driver).

    python summarize_io.py <run folder> [<run folder> ...] [--json out.json]

Prints per storage, op and thread count the median over all timed calls of all
rounds, the spread (min-max of the per-round medians), and Python / MATLAB.
"""

import argparse
import glob
import json
import os
import statistics
from collections import defaultdict

# the Python op that stands for each MATLAB op (copy to or from (y, x, z) included)
PAIRS = {
    "tiff_read": ["tiff_read", "tiff_read_tconv", "tiff_read_raw"],
    "tiff_write": ["tiff_write", "tiff_write_direct"],
    "zarr_read": ["zarr_read"],
    "zarr_write": ["zarr_write"],
}


def load(folder):
    runs = []
    for fn in sorted(glob.glob(os.path.join(folder, "*_r*.json"))):
        r = json.load(open(fn))
        r["round"] = int(os.path.basename(fn).split("_r")[1].split(".")[0])
        runs.append(r)
    procs = (
        [json.loads(line) for line in open(os.path.join(folder, "process_times.jsonl"))]
        if os.path.exists(os.path.join(folder, "process_times.jsonl"))
        else []
    )
    return runs, procs


def summarize(folder):
    runs, procs = load(folder)
    times = defaultdict(list)
    round_medians = defaultdict(list)
    for r in runs:
        rows = r["rows"] if isinstance(r["rows"], list) else [r["rows"]]
        for row in rows:
            ts = row["times"] if isinstance(row["times"], list) else [row["times"]]
            key = (row["impl"], row["op"], int(row["threads"]))
            times[key] += ts
            round_medians[key].append(statistics.median(ts))
    stats = {
        k: {
            "median": statistics.median(v),
            "n": len(v),
            "lo": min(round_medians[k]),
            "hi": max(round_medians[k]),
        }
        for k, v in times.items()
    }
    # start-up (launch to the first line; Python: imports included) and teardown (last line to exit)
    life = defaultdict(lambda: defaultdict(list))
    for r in runs:
        p = next((p for p in procs if p["impl"] == r["impl"] and p["round"] == r["round"]), None)
        if not p:
            continue
        first, end = (
            (r["t_imported"], r["t_end"]) if r["impl"] == "python" else (r["t_start"], r["t_end"])
        )
        # (MATLAB's posixtime of a datetime without a time zone is off by the UTC offset)
        tz = round((first - p["t_launch"]) / 900) * 900
        life[r["impl"]]["startup"].append(first - tz - p["t_launch"])
        life[r["impl"]]["teardown"].append(p["t_exit"] - (end - tz))
    node = (
        open(os.path.join(folder, "node.txt")).read().strip()
        if os.path.exists(os.path.join(folder, "node.txt"))
        else ""
    )
    return (
        stats,
        {k: {m: statistics.median(v) for m, v in d.items()} for k, d in life.items()},
        node,
        len(runs),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("folders", nargs="+")
    p.add_argument("--json")
    a = p.parse_args()
    out = {}
    for folder in a.folders:
        stats, life, node, nruns = summarize(folder)
        name = os.path.basename(folder.rstrip("/"))
        print(f"\n## {name}: {node}; {nruns} processes")
        head = f"{'op':10s} {'thr':>3s} {'MATLAB s':>9s} {'python op':18s} {'Python s':>9s}"
        print(f"{head} {'ratio':>6s}  spread M / P (round medians)")
        rows = []
        for mop, pops in PAIRS.items():
            for t in sorted({k[2] for k in stats}):
                m = stats.get(("matlab", mop, t))
                for pop in pops:
                    q = stats.get(("python", pop, t))
                    if not m or not q:
                        continue
                    ratio = q["median"] / m["median"]
                    rows.append(
                        {
                            "op": mop,
                            "threads": t,
                            "python_op": pop,
                            "matlab": m,
                            "python": q,
                            "ratio": ratio,
                        }
                    )
                    cells = f"{mop:10s} {t:3d} {m['median']:9.3f} {pop:18s} {q['median']:9.3f}"
                    spread = f"{m['lo']:.3f}-{m['hi']:.3f} / {q['lo']:.3f}-{q['hi']:.3f}"
                    print(f"{cells} {ratio:6.2f}  {spread}  (n {m['n']}/{q['n']})")
        for impl, d in life.items():
            up, down = d.get("startup", float("nan")), d.get("teardown", float("nan"))
            print(f"{impl}: start-up {up:.2f} s, teardown {down:.2f} s (medians)")
        out[name] = {"node": node, "rows": rows, "lifecycle": life}
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
