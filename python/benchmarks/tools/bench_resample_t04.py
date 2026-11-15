"""The resampling benchmark (T04), Python side (MATLAB: bench_resample_t04.m): in a process
pinned to the CPUs of ``<work>/cpus_<t>.txt`` (before any OpenMP runtime loads), each
operation once untimed and then timed:

    frame   resample_file of the tile by (2, 2, 2) (linear, uint16 tiff to uint16 tiff)
    zarr    resample_zarr of the tile as native zarr by (2, 2, 2) (batches 512^3, chunks
            256^3, local)

usage: python bench_resample_t04.py <work folder> <threads>
"""

import os
import shutil
import sys
import time

WORK, T = sys.argv[1], int(sys.argv[2])
with open(os.path.join(WORK, f"cpus_{T}.txt")) as f:
    os.sched_setaffinity(0, [int(c) for c in f.read().split()])

from petakit5d.testing.bench import mark, record  # noqa: E402
from petakit5d.tools import resample  # noqa: E402

mark("ready")
TILE = os.path.join(WORK, "tile.tif")
ZARR = os.path.join(WORK, "zarr", "python", "tile.zarr")


def frame(out):
    resample.resample_file(TILE, out, (2, 2, 2))


def zarr(out):
    resample.resample_zarr(ZARR, out, (2, 2, 2), parse_cluster=False)


def remove(path):
    if os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)
    elif os.path.exists(path):
        os.remove(path)


def timed(op, fn, ext):
    warm = os.path.join(WORK, f"w_py_{op}{ext}")
    fn(warm)
    remove(warm)
    out = os.path.join(WORK, f"t_py_{op}{ext}")
    t0 = time.perf_counter()
    fn(out)
    dt = time.perf_counter() - t0
    remove(out)
    record(op, [dt], threads=T)


timed("frame", frame, ".tif")
timed("zarr", zarr, ".zarr")
