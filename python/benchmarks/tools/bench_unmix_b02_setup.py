"""Inputs of the unmixing benchmark (B02): the demo's tiles copied to the work folder, one tile
pair as zarr in each side's layout (chunks 256^3), and the CPU lists the processes pin
themselves to (the first 1, 4 and 24 CPUs of the job's mask, or all it has).

usage: python bench_unmix_b02_setup.py <demo dataset folder> <work folder>
"""

import glob
import os
import shutil
import sys

from petakit5d.io import read_tiff
from petakit5d.testing import synth

TILE = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs"
    "_000x_003y_000z_0000t"
)


def main() -> None:
    src, work = sys.argv[1], sys.argv[2]
    demo = os.path.join(work, "demo")
    os.makedirs(demo, exist_ok=True)
    for f in sorted(glob.glob(os.path.join(src, "Scan_Iter_*_Cam*_ch0_*.tif"))):
        shutil.copyfile(f, os.path.join(demo, os.path.basename(f)))
    for cam in ("CamA", "CamB"):
        name = TILE.format(cam)
        vol = read_tiff(os.path.join(demo, f"{name}.tif"))
        synth.write_paired_inputs(os.path.join(work, "zarr"), name, vol, "zarr", (256, 256, 256))
    cpus = sorted(os.sched_getaffinity(0))
    for t in (1, 4, 24):
        with open(os.path.join(work, f"cpus_{t}.txt"), "w") as f:
            f.write(" ".join(str(c) for c in cpus[: min(t, len(cpus))]) + "\n")


if __name__ == "__main__":
    main()
