"""Inputs of the resampling benchmark (T04): a uint16 512 x 1800 x 500 tile (MATLAB's
[y x z]; (500, 512, 1800) here) of noise around a background, as tiff and as zarr in each
side's layout (chunks 256^3), and the CPU lists the processes pin themselves to (the first 4
and 24 CPUs of the job's mask, or all it has).

usage: python bench_resample_t04_setup.py <work folder>
"""

import os
import sys

import numpy as np

from petakit5d.testing import synth


def main() -> None:
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    rng = np.random.default_rng(4)
    vol = rng.integers(90, 400, size=(500, 512, 1800), dtype=np.uint16)
    synth.write_paired_inputs(work, "tile", vol, "tif")
    synth.write_paired_inputs(os.path.join(work, "zarr"), "tile", vol, "zarr", (256, 256, 256))
    cpus = sorted(os.sched_getaffinity(0))
    for t in (4, 24):
        with open(os.path.join(work, f"cpus_{t}.txt"), "w") as f:
            f.write(" ".join(str(c) for c in cpus[: min(t, len(cpus))]) + "\n")


if __name__ == "__main__":
    main()
