"""N5 benchmark (B09), the Python side, at each thread count (the process limited to that many
CPUs): zarr_to_n5 of FOLDER/t.zarr and n5_to_zarr of FOLDER/s.n5 (the defaults, the tasks in
this process), each result folder removed before each run (untimed).

usage: python n5_b09.py FOLDER --threads 24 --repeats 1
"""

from __future__ import annotations

import argparse
import os
import shutil
import time

from petakit5d.testing.bench import mark, record

mark("ready")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--threads", type=int, nargs="+", default=[24])
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.core.threads import refresh_cpus
    from petakit5d.io.n5 import n5_to_zarr, zarr_to_n5

    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    for t in a.threads:
        if cpus:
            os.sched_setaffinity(0, cpus[:t])
            refresh_cpus()
        tz, tn = [], []
        for _ in range(a.repeats):
            shutil.rmtree(os.path.join(a.folder, "n5"), ignore_errors=True)
            t0 = time.perf_counter()
            zarr_to_n5(os.path.join(a.folder, "t.zarr"), parse_cluster=False)
            tz.append(time.perf_counter() - t0)
            shutil.rmtree(os.path.join(a.folder, "zarr"), ignore_errors=True)
            t0 = time.perf_counter()
            n5_to_zarr(os.path.join(a.folder, "s.n5"), parse_cluster=False)
            tn.append(time.perf_counter() - t0)
        record("zarr_to_n5", tz, threads=t)
        record("n5_to_zarr", tn, threads=t)
    if cpus:
        os.sched_setaffinity(0, cpus)
        refresh_cpus()


if __name__ == "__main__":
    main()
