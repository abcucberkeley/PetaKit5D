"""Inputs of the FFT spectrum benchmark (B04): the first CamA ch0 tile of the demo cell image
dataset (PETAKIT5D_DEMO_DATA) as zarr in each side's layout (chunks 256^3), in
<work>/matlab/ds/ and <work>/python/ds/, and the CPU lists the processes pin themselves to (the
first 4 and 24 CPUs of the job's mask, or all it has).

usage: python fft_spectrum_b04_data.py <work folder>
"""

import os
import sys

from petakit5d.io import read_tiff
from petakit5d.testing import synth


def main() -> None:
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    demo = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "PetaKit5D_demo_cell_image_dataset")
    tile = sorted(n for n in os.listdir(demo) if "CamA_ch0" in n and n.endswith(".tif"))[0]
    vol = read_tiff(os.path.join(demo, tile))
    print(f"tile {tile}: (z, y, x) {vol.shape} {vol.dtype}", flush=True)
    for side in ("matlab", "python"):
        os.makedirs(os.path.join(work, side, "ds"), exist_ok=True)
    p = synth.write_paired_inputs(work, "ds/tile", vol, "zarr", (256, 256, 256))
    print(p, flush=True)
    cpus = sorted(os.sched_getaffinity(0))
    for t in (4, 24):
        with open(os.path.join(work, f"cpus_{t}.txt"), "w") as f:
            f.write(" ".join(str(c) for c in cpus[: min(t, len(cpus))]) + "\n")


if __name__ == "__main__":
    main()
