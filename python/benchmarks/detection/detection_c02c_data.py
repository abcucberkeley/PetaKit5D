"""Setup of the frame-parallel detection benchmark (C02c): the first frames of both channels of an
experiment of the point-detection demo (its raw and DS frames, its settings file), copied to the
work folder as an experiment (work/Ex1_bench_z0p5), and the PSF sigma of the experiment's
Detection3D.mat (work/sigma.txt). No folder or file name of the dataset is spelled out here.
usage: python detection_c02c_data.py <experiment folder> <work folder> <frames>
"""

import glob
import os
import shutil
import sys

import numpy as np

from petakit5d.io.matfile import read_mat

EX, WORK, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
dst = os.path.join(WORK, "Ex1_bench_z0p5")  # the loader reads the z step from the name
for c in ("ch1", "ch2"):
    os.makedirs(os.path.join(dst, c, "DS"), exist_ok=True)
    raw = sorted(glob.glob(os.path.join(EX, c, "*.tif")))[:N]
    ds = sorted(glob.glob(os.path.join(EX, c, "DS", "*.tif")))[:N]
    if len(raw) < N or len(ds) < N:
        sys.exit(f"{c}: fewer than {N} raw or DS frames")
    for p in raw:
        shutil.copyfile(p, os.path.join(dst, c, os.path.basename(p)))
    for p in ds:
        shutil.copyfile(p, os.path.join(dst, c, "DS", os.path.basename(p)))
for p in glob.glob(os.path.join(EX, "*.txt")):
    shutil.copyfile(p, os.path.join(dst, os.path.basename(p)))
fi = read_mat(os.path.join(EX, "ch1", "Analysis", "Detection3D.mat"), "frameInfo")["frameInfo"]
fi = [fi] if isinstance(fi, dict) else list(fi)
np.savetxt(os.path.join(WORK, "sigma.txt"), np.asarray(fi[0]["s"], float), fmt="%.17g")
