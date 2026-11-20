"""The C04a benchmark, Python side (MATLAB: bench_c04a.m): in one process, on a fresh copy of the
master condition (bench_c04a_setup.py), timed one after the other:

    load      load_condition_data (XR_loadConditionData3D)
    deskew    deskew_data with the wrapper's options (deskewData)
    offsets   correct_xz_offset_data (XR_correctXZoffsetData3D)

usage: python bench_c04a.py <work folder>
"""

import os
import shutil
import sys
import time

from petakit5d.detection import condition, offsets, wrapper
from petakit5d.testing.bench import mark, record

mark("ready")
WORK = sys.argv[1]
run = os.path.join(WORK, f"run_py_{os.getpid()}")
shutil.copytree(os.path.join(WORK, "master"), run)
movie = os.path.join(run, "Ex1_demo_z0p5")
t0 = time.perf_counter()
data = condition.load_condition_data(movie, ["ch1", "ch2"], ["gfp", "rfp"])
t1 = time.perf_counter()
wrapper.deskew_data(data, overwrite=False, skew_angle=31.5, aname="Analysis")
t2 = time.perf_counter()
offsets.correct_xz_offset_data(data, overwrite=False)
t3 = time.perf_counter()
record("load", [t1 - t0])
record("deskew", [t2 - t1])
record("offsets", [t3 - t2])
shutil.rmtree(run)
