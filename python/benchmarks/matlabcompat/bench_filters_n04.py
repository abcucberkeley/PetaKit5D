"""The N04 benchmark operations, Python side (the MATLAB side is bench_filters_n04.m): the same
data in (z, y, x) C order, the same op names and thread counts.

    python bench_filters_n04.py --threads 1 4 24 --repeats 3
"""

import argparse

import numpy as np

from petakit5d.matlabcompat import filters as F
from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
ap.add_argument("--repeats", type=int, default=3)
args = ap.parse_args()

# MATLAB's mod(y*7 + x*13 + z*17, 1000)/1000 in single, as (z, y, x)
y = np.arange(1, 513, dtype=np.float32)[None, :, None]
x = np.arange(1, 513, dtype=np.float32)[None, None, :]
z = np.arange(1, 401, dtype=np.float32)[:, None, None]
s = y * np.float32(7) + x * np.float32(13) + z * np.float32(17)
V = np.ascontiguousarray(np.mod(s, np.float32(1000)) / np.float32(1000))
U = np.ascontiguousarray(np.mod(s, np.float32(4000)).astype(np.uint16))
y2 = np.arange(1, 2049, dtype=np.float32)[:, None]
x2 = np.arange(1, 2049, dtype=np.float32)[None, :]
P2 = np.ascontiguousarray(
    np.mod(y2 * np.float32(7) + x2 * np.float32(13), np.float32(1000)) / np.float32(1000)
)
P1 = np.ascontiguousarray(P2[:1024, :1024])

mark("ready")
r = args.repeats
for t in args.threads:
    record(
        "ig3_unmix_single",
        timeit(lambda t=t: F.imgaussfilt3(V, 1.0, num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "ig3_puncta_uint16",
        timeit(lambda t=t: F.imgaussfilt3(U, 2.5, num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "ig2_phaseff_2048",
        timeit(lambda t=t: F.imgaussfilt(P2, 100.0, num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "ig2_phaseff_1024",
        timeit(lambda t=t: F.imgaussfilt(P1, 100.0, num_threads=t), repeats=r),
        threads=t,
    )
