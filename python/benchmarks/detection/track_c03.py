"""Tracking benchmark (C03), the Python side: run_tracking of FOLDER/Analysis/Detection3D.mat with
the demos' settings, the result removed before each run (untimed).

usage: python track_c03.py FOLDER --repeats 1
"""

from __future__ import annotations

import argparse
import os
import time

from petakit5d.testing.bench import mark, record

mark("ready")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.detection.tracking import run_tracking, track_settings

    rp = os.path.join(a.folder, "Analysis") + os.sep
    data = {"source": a.folder + os.sep}
    settings = track_settings((3, 6), 2)
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    times = []
    for _ in range(a.repeats):
        fn = os.path.join(rp, "trackedFeatures.mat")
        if os.path.exists(fn):
            os.remove(fn)
        t0 = time.perf_counter()
        run_tracking(
            data, settings, file_name="trackedFeatures.mat", overwrite=True, results_path=[rp]
        )
        times.append(time.perf_counter() - t0)
    record("track", times, threads=cpus)


if __name__ == "__main__":
    main()
