"""The FFT spectrum benchmark (B04), Python side (MATLAB: fft_spectrum_b04_ml.m): in a process
pinned to the CPUs of ``<work>/cpus_<t>.txt``, fft_spectrum_dataset of <work>/python/ds/ as
demo_useful_tools.m calls XR_fftSpectrumComputingWrapper (zarr input, outSize 1001^3, voxels of
0.1 * 0.488 / 1.33 um, no 3D stack, tasks in this process), once untimed and then timed, the
result folder removed after each (untimed).

usage: python fft_spectrum_b04.py <work folder> <threads>
"""

import os
import shutil
import sys
import time

WORK, T = sys.argv[1], int(sys.argv[2])
with open(os.path.join(WORK, f"cpus_{T}.txt")) as f:
    os.sched_setaffinity(0, [int(c) for c in f.read().split()])

from petakit5d.matlabcompat.elementary import sind  # noqa: E402
from petakit5d.testing.bench import mark, record  # noqa: E402
from petakit5d.tools import fft_spectrum_dataset  # noqa: E402

mark("ready")
DS = os.path.join(WORK, "python", "ds")


def run() -> None:
    fft_spectrum_dataset(
        [DS + "/"],
        result_dir_name="FFT",
        xy_pixel_size=0.108,
        dz=0.3 * float(sind(32.45)),
        zarr_file=True,
        out_size=(1001, 1001, 1001),
        out_pixel_size=0.1 * 0.488 / 1.33,
        save_3d_stack=False,
        parse_cluster=False,
        master_compute=True,
        cpus_per_task=24,
    )


for timed in (False, True):
    t0 = time.perf_counter()
    run()
    dt = time.perf_counter() - t0
    shutil.rmtree(os.path.join(DS, "FFT"))
    if timed:
        record("dataset", [dt], threads=T)
