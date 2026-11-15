"""MIP mask benchmark (B10), the Python side, at each thread count (the process limited to that
many CPUs):

- ``masks``: generate_mip_masks of FOLDER/mips (the defaults, the tasks in this process), the
  masks folder removed before each run (untimed); ``--no-figures`` leaves out the masks' figures,
  as MATLAB does under the job framework's ``-nojvm``;
- ``seg``: visualize_otf_mask_segmentation of FOLDER/psf/560_c.tif (0.88, skewed, the figure
  saved and not shown), as the demo calls it.

usage: python mip_masks_b10.py FOLDER {masks,seg} --threads 1 24 --repeats 1 [--no-figures]
"""

from __future__ import annotations

import argparse
import os
import shutil
import time

from petakit5d.testing.bench import mark, record

mark("ready")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("what", choices=("masks", "seg"))
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 24])
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args()
    from petakit5d.core.threads import refresh_cpus
    from petakit5d.tools import figures
    from petakit5d.tools.figures import visualize_otf_mask_segmentation
    from petakit5d.tools.mip_masks import generate_mip_masks

    if a.no_figures:
        figures.available = lambda: False
    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    mips = os.path.join(a.folder, "mips")
    psf = os.path.join(a.folder, "psf", "560_c.tif")
    for t in a.threads:
        if cpus:
            os.sched_setaffinity(0, cpus[:t])
            refresh_cpus()
        times = []
        for _ in range(a.repeats):
            if a.what == "masks":
                shutil.rmtree(os.path.join(mips, "masks"), ignore_errors=True)
                t0 = time.perf_counter()
                generate_mip_masks(mips, parse_cluster=False)
            else:
                shutil.rmtree(os.path.join(a.folder, "psf", "figures"), ignore_errors=True)
                t0 = time.perf_counter()
                visualize_otf_mask_segmentation(psf, 0.88, True, visible=False, save_fig=True)
            times.append(time.perf_counter() - t0)
        record(
            "generate_mip_masks" if a.what == "masks" else "otf_mask_segmentation", times, threads=t
        )
    if cpus:
        os.sched_setaffinity(0, cpus)
        refresh_cpus()


if __name__ == "__main__":
    main()
