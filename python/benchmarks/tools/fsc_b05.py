"""FSC benchmark (B05), the Python side: fsc_analysis_dataset of FOLDER/iters with the parameters
of demo_useful_tools.m's FSC section, at each thread count, the process limited to that many
CPUs; the results and figures removed before each run (untimed), the kept window and shells
forgotten (as MATLAB's ``clear XR_FSC``). ``--no-figures``: without matplotlib's figures.

usage: python fsc_b05.py FOLDER --threads 1 24 --repeats 1 [--no-figures]
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import time

from petakit5d.testing.bench import mark, record

mark("ready")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 24])
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args()
    from petakit5d.core.threads import refresh_cpus
    from petakit5d.matlabcompat.elementary import sind
    from petakit5d.tools import figures
    from petakit5d.tools import fsc as FS

    if a.no_figures:
        figures.available = lambda: False
    iters = os.path.join(a.folder, "iters") + "/"
    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    dz = 0.3 * float(sind(32.45))
    for t in a.threads:
        if cpus:
            os.sched_setaffinity(0, cpus[:t])
            refresh_cpus()
        times = []
        for _ in range(a.repeats):
            for d in ("FSCs", "figures"):
                shutil.rmtree(iters + d, ignore_errors=True)
            FS.clear_caches()
            t0 = time.perf_counter()
            FS.fsc_analysis_dataset(
                [iters],
                result_dir_name="FSCs",
                xy_pixel_size=0.108,
                dz=dz,
                dr=5,
                dtheta=math.pi / 6,
                res_thresh_method="one-bit",
                res_thresh=0.2,
                half_size=(201, 201, 201),
                input_bbox=None,
                channel_patterns=("Iter",),
                channels=(488,),
                suffix="decon",
                iter_interval=2,
                parse_cluster=False,
                master_compute=True,
                cpus_per_task=24,
            )
            times.append(time.perf_counter() - t0)
        record("fsc_analysis_dataset", times, threads=t)
    if cpus:
        os.sched_setaffinity(0, cpus)
        refresh_cpus()


if __name__ == "__main__":
    main()
