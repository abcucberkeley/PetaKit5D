"""The chromatic shift benchmark (B03), Python side (MATLAB: chromatic_b03_ml.m): in a process
pinned to the CPUs of ``<work>/cpus_<t>.txt``, chromatic_shift_correction_dataset of <work>/ds/
as demo_useful_tools.m calls XR_chromatic_shift_correction_data_wrapper (offsets [0 0 0; 1 -5 -2]
in [y x z], 'valid' mode, tiff results, tasks in this process), then chromatic_shift_estimation
of the demo's two PSFs (maxOffset 20, the demo's), each once untimed and then timed, the results
removed after each (untimed).

usage: python chromatic_b03.py <work folder> <threads>
"""

import os
import shutil
import sys
import time

WORK, T = sys.argv[1], int(sys.argv[2])
with open(os.path.join(WORK, f"cpus_{T}.txt")) as f:
    os.sched_setaffinity(0, [int(c) for c in f.read().split()])

from petakit5d.testing.bench import mark, record  # noqa: E402
from petakit5d.tools import chromatic_shift_correction_dataset, chromatic_shift_estimation  # noqa: E402

mark("ready")
DS = os.path.join(WORK, "ds")
PSF = [os.path.join(WORK, "psf", n) for n in ("560_c.tif", "488_2_c.tif")]
RES = os.path.join(WORK, "py_measurement.mat")


def correct() -> None:
    chromatic_shift_correction_dataset(
        [DS + "/"],
        chromatic_offset=((0, 0, 0), (-2, 1, -5)),
        result_dir_name="Chromatic_Shift_Corrected",
        mode="valid",
        channel_patterns=("CamA", "CamB"),
        max_offset=(20, 20, 20),
        parse_cluster=False,
        master_compute=True,
        cpus_per_task=3,
    )


def estimate() -> None:
    chromatic_shift_estimation(PSF, result_fullpath=RES, max_offset=(20, 20, 20), num_threads=T)


for timed in (False, True):
    t0 = time.perf_counter()
    correct()
    dt = time.perf_counter() - t0
    shutil.rmtree(os.path.join(DS, "Chromatic_Shift_Corrected"))
    if timed:
        record("dataset", [dt], threads=T)
for timed in (False, True):
    t0 = time.perf_counter()
    estimate()
    dt = time.perf_counter() - t0
    os.remove(RES)
    if timed:
        record("estimation", [dt], threads=T)
