"""PSF analysis benchmark (B06), the Python side, at each thread count (the process limited to
that many CPUs):

- ``analysis``: psf_analysis_dataset of FOLDER/PSF with the demo's parameters (demo_useful_tools'
  "PSF analysis": dz 0.3, channel pattern '488', deskew, the 515 nm Richards & Wolf PSF), the
  figures not shown, the tasks in this process; the DS folder removed before each run
  (untimed);
- ``detection``: psf_detection_and_cropping of FOLDER/beads/beads.tif (crops of 31 x 64 x 64,
  isolation distance 40 x 64 x 64), its result folder removed before each run (untimed).

usage: python psf_analysis_b06.py FOLDER {analysis,detection} --threads 1 24 --repeats 1
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
    ap.add_argument("what", choices=("analysis", "detection"))
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 24])
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.core.threads import refresh_cpus
    from petakit5d.tools.psf_analysis import psf_analysis_dataset, psf_detection_and_cropping

    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    psf = os.path.join(a.folder, "PSF") + "/"
    beads = os.path.join(a.folder, "beads")
    for t in a.threads:
        if cpus:
            os.sched_setaffinity(0, cpus[:t])
            refresh_cpus()
        times = []
        for _ in range(a.repeats):
            if a.what == "analysis":
                shutil.rmtree(os.path.join(psf, "DS"), ignore_errors=True)
                t0 = time.perf_counter()
                psf_analysis_dataset(
                    [psf],
                    xy_pixel_size=0.108,
                    dz=0.3,
                    skew_angle=32.45,
                    channel_patterns=("488",),
                    channels=(488,),
                    deskew=True,
                    objective_scan=False,
                    source_str="Demo PSF",
                    rw_fn=(f"{psf}RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif",),
                    visible=False,
                    parse_cluster=False,
                    master_compute=True,
                    cpus_per_task=24,
                )
            else:
                shutil.rmtree(os.path.join(beads, "Cropped"), ignore_errors=True)
                t0 = time.perf_counter()
                psf_detection_and_cropping(
                    os.path.join(beads, "beads.tif"),
                    os.path.join(beads, "Cropped") + "/",
                    crop_size=(31, 64, 64),
                    dist_thresh=(40, 64, 64),
                    prefix="beads",
                )
            times.append(time.perf_counter() - t0)
        record("psf_analysis" if a.what == "analysis" else "psf_detection", times, threads=t)
    if cpus:
        os.sched_setaffinity(0, cpus)
        refresh_cpus()


if __name__ == "__main__":
    main()
