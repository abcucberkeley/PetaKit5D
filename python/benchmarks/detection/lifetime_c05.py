"""Lifetime analysis benchmark (C05), the Python side: ``run_lifetime_analysis`` (MATLAB's
runLifetimeAnalysis3D with 'overwrite' true: the tracks read again each run) and
``plot_intensity_cohorts`` on the demo condition prepared by lifetime_c05_data.py, each run
timed whole (the movies' data from the folder: source, frame rate from the frames' msec
tokens, movie length, channels).

usage: python lifetime_c05.py FOLDER ANALYSIS --repeats 1
"""

from __future__ import annotations

import argparse
import os
import re
import time

from petakit5d.testing.bench import mark, record

mark("ready")


def condition(folder: str) -> list[dict]:
    cond = os.path.join(folder, sorted(os.listdir(folder))[0])
    data = []
    for ex in sorted(e for e in os.listdir(cond) if e.startswith("Ex")):
        ch1 = os.path.join(cond, ex, "ch1")
        frames = sorted(f for f in os.listdir(ch1) if f.endswith(".tif"))
        ms = [int(re.search(r"_(\d+)msec_", f).group(1)) for f in frames[:2]]
        data.append({
            "source": ch1 + "/", "framerate": (ms[1] - ms[0]) / 1000, "movieLength": len(frames),
            "channels": [os.path.join(cond, ex, c) + "/" for c in ("ch1", "ch2")],
        })  # fmt: skip
    return data


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("analysis")
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.detection import lifetime

    data = condition(a.folder)
    t1, t2 = [], []
    for _ in range(a.repeats):
        t0 = time.perf_counter()
        lifetime.run_lifetime_analysis(data, analysis_path=a.analysis, overwrite=True)
        t1.append(time.perf_counter() - t0)
        t0 = time.perf_counter()
        lifetime.plot_intensity_cohorts(data, analysis_path=a.analysis)
        t2.append(time.perf_counter() - t0)
    record("lifetime_analysis", t1, threads=1)
    record("intensity_cohorts", t2, threads=1)


if __name__ == "__main__":
    main()
