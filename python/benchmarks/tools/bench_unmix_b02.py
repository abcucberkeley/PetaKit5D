"""The unmixing benchmark (B02), Python side (MATLAB: bench_unmix_b02.m): in a process pinned to
the CPUs of ``<work>/cpus_<t>.txt`` (before any OpenMP runtime loads), each operation once
untimed and then timed:

    frame         unmix_channels_file of one demo tile pair (tiff, [1 -0.1], backgrounds
                  [100 100], offset 1, uint16)
    zarr          unmix_channels_zarr of the same pair as native zarr (batches 1024^3, chunks
                  256^3, local)
    dataset_demo  (with --dataset, timed once) unmix_channels_dataset of the demo's 8 tile
                  pairs with the demo's settings (local)

usage: python bench_unmix_b02.py <work folder> <threads> [--dataset]
"""

import os
import shutil
import sys
import time

WORK, T = sys.argv[1], int(sys.argv[2])
DATASET = "--dataset" in sys.argv[3:]
with open(os.path.join(WORK, f"cpus_{T}.txt")) as f:
    os.sched_setaffinity(0, [int(c) for c in f.read().split()])

from petakit5d.testing.bench import mark, record  # noqa: E402
from petakit5d.tools import unmix  # noqa: E402

mark("ready")
DEMO = os.path.join(WORK, "demo")
TILE = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs"
    "_000x_003y_000z_0000t"
)
FRAMES = [os.path.join(DEMO, TILE.format(c) + ".tif") for c in ("CamA", "CamB")]
ZARRS = [os.path.join(WORK, "zarr", "python", TILE.format(c) + ".zarr") for c in ("CamA", "CamB")]
SETTINGS = {"const_background": (100, 100), "const_offset": 1, "save_16bit": True}


def frame(name):
    unmix.unmix_channels_file(FRAMES, (1, -0.1), result_dir_name=name, **SETTINGS)


def zarr(name):
    unmix.unmix_channels_zarr(ZARRS, (1, -0.1), result_dir_name=name, parse_cluster=False)


def dataset(name):
    unmix.unmix_channels_dataset(
        DEMO,
        unmix_factors=(1, -0.1),
        mode="linear",
        result_dir_name=name,
        channel_patterns=("CamA", "CamB"),
        channel_ind=0,
        parse_cluster=False,
        master_compute=True,
        **SETTINGS,
    )


def timed(op, fn, folder, warm=True):
    if warm:
        fn("Uw_py")
        shutil.rmtree(os.path.join(folder, "Uw_py"), ignore_errors=True)
    t0 = time.perf_counter()
    fn("Ut_py")
    dt = time.perf_counter() - t0
    shutil.rmtree(os.path.join(folder, "Ut_py"), ignore_errors=True)
    record(op, [dt], threads=T)


timed("frame", frame, DEMO)
timed("zarr", zarr, os.path.join(WORK, "zarr", "python"))
if DATASET:
    timed("dataset_demo", dataset, DEMO, warm=False)
