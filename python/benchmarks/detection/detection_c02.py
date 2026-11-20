"""The detection benchmark (C02), Python side (MATLAB: detection_c02_ml.m): the background estimate
and the lowSNR detector on one demo frame as XR_runDetection3D runs them (sigma [1.3758 1.105],
'BackgroundInfo' [bg 3]), once untimed and then timed.
usage: python detection_c02.py <work folder> <threads>
"""

import os
import sys
import time

WORK, T = sys.argv[1], int(sys.argv[2])

import numpy as np  # noqa: E402

from petakit5d.detection import detect  # noqa: E402
from petakit5d.io import read_tiff  # noqa: E402
from petakit5d.testing.bench import mark, record  # noqa: E402

mark("ready")
vol = read_tiff(os.path.join(WORK, "frame.tif")).astype(np.float64)
SIGMA = (1.3758153893715306, 1.1049523723652759)


def run() -> int:
    detect.clear_threshold_cache()
    bi = detect.estimate_background(vol, 1)
    ps, _ = detect.point_source_detection(vol, SIGMA, background_info=np.r_[bi, 3.0], num_threads=T)
    return int(np.asarray(ps["x"]).size)


run()
t0 = time.perf_counter()
n = run()
record("detect_frame", [time.perf_counter() - t0], threads=T)
print("detections", n)
mark("done")
