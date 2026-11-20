"""The SIM benchmark (C06), Python side (MATLAB: sim_recon_c06_ml.m): in a process pinned to the
CPUs of ``<work>/cpus_<t>.txt``, sim_recon_file of <work>/data/Cell_5phase_small.tif with the
PSF file (the OTF made in the call), the 2020 workflow's values (the defaults), DS, on the CPU,
once untimed and then timed, the result folder removed after each (untimed).

usage: python sim_recon_c06.py <work folder> <threads>
"""

import os
import shutil
import sys
import time

WORK, T = sys.argv[1], int(sys.argv[2])
with open(os.path.join(WORK, f"cpus_{T}.txt")) as f:
    os.sched_setaffinity(0, [int(c) for c in f.read().split()])

from petakit5d.sim import sim_recon_file  # noqa: E402
from petakit5d.testing.bench import mark, record  # noqa: E402

mark("ready")
DATA = os.path.join(WORK, "data")


def run() -> None:
    sim_recon_file(
        os.path.join(DATA, "Cell_5phase_small.tif"),
        os.path.join(DATA, "PSF_5phase.tif"),
        use_gpu=False,
        ds=True,
        num_threads=T,
    )


for timed in (False, True):
    t0 = time.perf_counter()
    run()
    dt = time.perf_counter() - t0
    shutil.rmtree(os.path.join(DATA, "sim_recon"))
    if timed:
        record("frame", [dt], threads=T)
