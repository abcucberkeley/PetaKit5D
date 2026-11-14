"""Puncta removal benchmark (B08), the Python side: puncta_removal_zarr of FOLDER/expanded.zarr
(the default method and parameters, the batches in this process) at each thread count, the
process limited to that many CPUs; the result folder removed before each run (untimed).

usage: python puncta_b08.py FOLDER --threads 1 4 24 --repeats 1
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
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.core.threads import refresh_cpus
    from petakit5d.tools.puncta_removal import puncta_removal_zarr

    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    fn = os.path.join(a.folder, "expanded.zarr")
    out = os.path.join(a.folder, "Puncta_Removed")
    for t in a.threads:
        if cpus:
            os.sched_setaffinity(0, cpus[:t])
            refresh_cpus()
        times = []
        for _ in range(a.repeats):
            shutil.rmtree(out, ignore_errors=True)
            t0 = time.perf_counter()
            puncta_removal_zarr(fn, parse_cluster=False)
            times.append(time.perf_counter() - t0)
        record("puncta_removal_zarr", times, threads=t)
    if cpus:
        os.sched_setaffinity(0, cpus)
        refresh_cpus()


if __name__ == "__main__":
    main()
