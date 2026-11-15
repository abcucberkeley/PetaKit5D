"""Flat-field benchmark (B07), the Python side, at each thread count (the process limited to
that many CPUs), the results removed before each run (untimed):

- ``flatfield_background_wrapper``: flatfield_background_dataset of FOLDER/dark (two cameras,
  in this process);
- ``phase_flat_field``: phase_image_flatfield_correction of each time point of FOLDER/phase
  (the demo's loop, offset 6000) into FOLDER/phase/ff_corrected_1/.

With B07_KEEP set, the last run's phase results are copied to $B07_KEEP/python.

usage: python flatfield_b07.py FOLDER --threads 1 24 --repeats 1
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
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 24])
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.core.threads import refresh_cpus
    from petakit5d.flatfield.basic import phase_image_flatfield_correction
    from petakit5d.flatfield.estimation import flatfield_background_dataset

    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    dark = os.path.join(a.folder, "dark")
    phase = os.path.join(a.folder, "phase") + "/"
    out = phase + "ff_corrected_1/"
    tiles = sorted(f for f in os.listdir(phase) if f.endswith(".tif"))
    groups = [[phase + f for f in tiles if f"Iter_{i:04d}" in f] for i in range(2)]
    for t in a.threads:
        if cpus:
            os.sched_setaffinity(0, cpus[:t])
            refresh_cpus()
        times = []
        for _ in range(a.repeats):
            shutil.rmtree(os.path.join(dark, "matlab_flat_field_estimation"), ignore_errors=True)
            t0 = time.perf_counter()
            flatfield_background_dataset(
                dark, channel_patterns=("CamA_ch0", "CamB_ch0"), parse_cluster=False
            )
            times.append(time.perf_counter() - t0)
        record("flatfield_background_wrapper", times, threads=t)
        times = []
        for _ in range(a.repeats):
            shutil.rmtree(out, ignore_errors=True)
            os.makedirs(out)
            t0 = time.perf_counter()
            for fns in groups:
                phase_image_flatfield_correction(fns, out, 6000)
            times.append(time.perf_counter() - t0)
        record("phase_flat_field", times, threads=t)
    if cpus:
        os.sched_setaffinity(0, cpus)
        refresh_cpus()
    keep = os.environ.get("B07_KEEP")
    if keep:
        shutil.copytree(out, os.path.join(keep, "python"), dirs_exist_ok=True)


if __name__ == "__main__":
    main()
