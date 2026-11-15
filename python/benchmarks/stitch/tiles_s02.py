"""Tile preparation benchmark (S02), the Python side: stitch_process_tiles (process_tiles) of the
16 tiles with the flat-field and erosion function file, tiff tiles to zarr (``tiff_tiles``) and
zarr tiles (``zarr_tiles``), the tasks in this process, at each thread count (the process
limited to that many CPUs); the result folders removed before each run (untimed). Times are
per tile.

usage: python tiles_s02.py FOLDER --threads 4 24 --repeats 1
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import time

from petakit5d.testing.bench import mark, record

mark("ready")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--threads", type=int, nargs="+", default=[4, 24])
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.core.threads import refresh_cpus
    from petakit5d.stitch.tiles import process_tiles

    d = os.path.abspath(a.folder).replace("\\", "/")
    tifs = sorted(glob.glob(f"{d}/tif/*.tif"))
    zarrs = sorted(glob.glob(f"{d}/zarr/python/*.zarr"))
    common = {
        "block_size": (125, 250, 250),
        "process_fun_path": [[f"{d}/fun.mat"]],
        "parse_cluster": False,
    }
    runs = [
        ("tiff_tiles", tifs, f"{d}/tif/zarr", {"zarr_file": False, "result_dir_name": "zarr"}),
        (
            "zarr_tiles",
            zarrs,
            f"{d}/zarr/python/proc",
            {"zarr_file": True, "result_dir_name": "proc"},
        ),
    ]
    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    for t in a.threads:
        if cpus:
            os.sched_setaffinity(0, cpus[:t])
            refresh_cpus()
        for op, files, out, extra in runs:
            times = []
            for _ in range(a.repeats):
                shutil.rmtree(out, ignore_errors=True)
                t0 = time.perf_counter()
                process_tiles(files, **common, **extra)
                times.append((time.perf_counter() - t0) / len(files))
            record(op, times, threads=t)
    if cpus:
        os.sched_setaffinity(0, cpus)
        refresh_cpus()


if __name__ == "__main__":
    main()
