"""The frame-parallel detection benchmark (C02c), Python side (MATLAB: detection_c02c_ml.m):
run_detection_frame_parallel on the setup's condition ('lowSNR' without the background
correction, XR_runDetTrack3d's detection), the results removed before each run (untimed); once
untimed and then timed.
usage: python detection_c02c.py <work folder> <threads>
"""

import os
import shutil
import sys
import time

WORK, T = sys.argv[1], int(sys.argv[2])

import numpy as np  # noqa: E402

from petakit5d.detection import detect, run  # noqa: E402
from petakit5d.detection.condition import load_condition_data  # noqa: E402
from petakit5d.testing.bench import mark, record  # noqa: E402

mark("ready")
sigma = np.loadtxt(os.path.join(WORK, "sigma.txt"), ndmin=2)
data = load_condition_data(
    os.path.join(WORK, "Ex1_bench_z0p5") + os.sep, ["ch1", "ch2"], ["gfp", "rfp"]
)
rp = data[0].source + "Analysis" + os.sep


def once() -> int:
    if os.path.isdir(rp):
        shutil.rmtree(rp)
    detect.clear_threshold_cache()
    (entries,) = run.run_detection_frame_parallel(
        data, sigma, [rp], detection_method="lowSNR", background_correction=False, num_threads=T
    )
    return sum(int(np.asarray(e["x"]).size) for e in entries)


once()
if os.path.isdir(rp):
    shutil.rmtree(rp)
t0 = time.perf_counter()
n = once()
record("detect_movie", [time.perf_counter() - t0], threads=T)
print("detections", n)
mark("done")
