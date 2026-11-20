"""Setup of the C04a benchmark: the demo dataset's first experiment, its first frames of both
channels, copied as ``<work>/master/Ex1_demo_z0p5/ch1|ch2/`` (the folders are found as listed, so
that no dataset name is written here).

usage: python bench_c04a_setup.py <work folder> [frames]
"""

import glob
import os
import shutil
import sys

WORK = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 3
root = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "PetaKit5D_point_detection_demo_dataset")
cond = sorted(e.path for e in os.scandir(root) if e.is_dir() and not e.name.startswith("."))[0]
ex = sorted(glob.glob(os.path.join(cond, "Ex1*")))[0]
for ch in ("ch1", "ch2"):
    dst = os.path.join(WORK, "master", "Ex1_demo_z0p5", ch)
    os.makedirs(dst, exist_ok=True)
    for p in sorted(glob.glob(os.path.join(ex, ch, "*.tif")))[:N]:
        shutil.copyfile(p, os.path.join(dst, os.path.basename(p)))
