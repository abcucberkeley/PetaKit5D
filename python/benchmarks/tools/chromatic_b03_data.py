"""Inputs of the chromatic shift benchmark (B03): the eight tiles of the first time point of the
demo cell image dataset (PETAKIT5D_DEMO_DATA; CamA and CamB) in <work>/ds/, its two PSFs (CamA:
560_c.tif, CamB: 488_2_c.tif) in <work>/psf/, and the CPU lists the processes pin themselves to
(the first 4 and 24 CPUs of the job's mask, or all it has). Both sides read the same tiffs.

usage: python chromatic_b03_data.py <work folder>
"""

import os
import shutil
import sys


def main() -> None:
    work = sys.argv[1]
    demo = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "PetaKit5D_demo_cell_image_dataset")
    ds = os.path.join(work, "ds")
    psf = os.path.join(work, "psf")
    os.makedirs(ds, exist_ok=True)
    os.makedirs(psf, exist_ok=True)
    tiles = sorted(
        n for n in os.listdir(demo) if n.startswith("Scan_Iter_0000_0000_") and n.endswith(".tif")
    )
    for n in tiles:
        shutil.copyfile(os.path.join(demo, n), os.path.join(ds, n))
    for n in ("560_c.tif", "488_2_c.tif"):
        shutil.copyfile(os.path.join(demo, "PSF", n), os.path.join(psf, n))
    print(f"{len(tiles)} tiles, 2 PSFs", flush=True)
    cpus = sorted(os.sched_getaffinity(0))
    for t in (4, 24):
        with open(os.path.join(work, f"cpus_{t}.txt"), "w") as f:
            f.write(" ".join(str(c) for c in cpus[: min(t, len(cpus))]) + "\n")


if __name__ == "__main__":
    main()
