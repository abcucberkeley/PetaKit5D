"""Inputs of the SIM benchmark (C06): the 2020 lattice SIM data of the demo data
(PETAKIT5D_DEMO_DATA/sim3D_demo_data: Cell_5phase_small.tif and PSF_5phase.tif) copied to
<work>/data/ (tiff: one file for both sides), and the CPU lists the processes pin themselves to
(the first 4 and 24 CPUs of the job's mask, or all it has).

usage: python sim_recon_c06_data.py <work folder>
"""

import os
import shutil
import sys


def main() -> None:
    work = sys.argv[1]
    data = os.path.join(work, "data")
    os.makedirs(data, exist_ok=True)
    src = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "sim3D_demo_data")
    for name in ("Cell_5phase_small.tif", "PSF_5phase.tif"):
        shutil.copy(os.path.join(src, name), data)
    cpus = sorted(os.sched_getaffinity(0))
    for t in (4, 24):
        with open(os.path.join(work, f"cpus_{t}.txt"), "w") as f:
            f.write(" ".join(str(c) for c in cpus[: min(t, len(cpus))]) + "\n")


if __name__ == "__main__":
    main()
