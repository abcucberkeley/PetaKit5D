"""The N01 benchmark operations, Python side (the MATLAB side is bench_resize_n01.m): the same
uint16 volume in (z, y, x) C order, the same op names, sizes and thread counts.

    python bench_resize_n01.py --threads 4 24 --repeats 2
"""

import argparse

import numpy as np

from petakit5d.core import dtypes
from petakit5d.matlabcompat import resize as R
from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, nargs="+", default=[4, 24])
ap.add_argument("--repeats", type=int, default=2)
args = ap.parse_args()

# MATLAB's uint16(mod(y*7 + x*13 + z*17, 4000)) in single, as (z, y, x)
y = np.arange(1, 513, dtype=np.float32)[None, :, None]
x = np.arange(1, 1801, dtype=np.float32)[None, None, :]
z = np.arange(1, 501, dtype=np.float32)[:, None, None]
V = np.ascontiguousarray(
    np.mod(y * np.float32(7) + x * np.float32(13) + z * np.float32(17), np.float32(4000))
)
V = np.ascontiguousarray(np.asarray(dtypes.cast(V, "uint16")))
sz = V.shape


def rnd(v):
    return int(dtypes.mround(np.float64(v)))


half = tuple(rnd(s * 0.5) for s in sz)
z2 = (rnd(sz[0] * 2), sz[1], sz[2])
zh = (rnd(sz[0] / 2), sz[1], sz[2])

mark("ready")
r = args.repeats
for t in args.threads:
    record(
        "rs3_half_cubic",
        timeit(lambda t=t: R.imresize3(V, half, num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "rs3_half_linear",
        timeit(lambda t=t: R.imresize3(V, half, "linear", num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "rs3_z2_cubic", timeit(lambda t=t: R.imresize3(V, z2, num_threads=t), repeats=r), threads=t
    )
    record(
        "rs3_z2_linear",
        timeit(lambda t=t: R.imresize3(V, z2, "linear", num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "rs3_zhalf_linear",
        timeit(lambda t=t: R.imresize3(V, zh, "linear", num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "fast_half_cubic",
        timeit(lambda t=t: R.imresize3_fast(V, half, num_threads=t), repeats=r),
        threads=t,
    )
    record(
        "fast_z2_cubic",
        timeit(lambda t=t: R.imresize3_fast(V, z2, num_threads=t), repeats=r),
        threads=t,
    )
