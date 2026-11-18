"""matlabcompat.reduce's sums before and after the per-output NaN payloads (reduce_x05.json runs
this script on two source trees): single volumes without NaN and with a few, summed along each
axis and over all elements, their mean, var and cumsum.

    python bench_reduce_x05.py --threads 1 4 24 --repeats 3
"""

import argparse

import numpy as np

from petakit5d.matlabcompat import reduce as R
from petakit5d.matlabcompat import stats as S
from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
ap.add_argument("--repeats", type=int, default=3)
args = ap.parse_args()

# mod(y*7 + x*13 + z*17, 1000)/1000 in single, (z, y, x)
z = np.arange(1, 201, dtype=np.float32)[:, None, None]
y = np.arange(1, 513, dtype=np.float32)[None, :, None]
x = np.arange(1, 513, dtype=np.float32)[None, None, :]
V = np.ascontiguousarray(
    np.mod(y * np.float32(7) + x * np.float32(13) + z * np.float32(17), np.float32(1000))
    / np.float32(1000)
)
# a few NaNs of two payloads (one in 4,096 voxels), as a masked volume has
N = V.copy()
flat = N.reshape(-1)
flat[::4096] = np.uint32(0xFFC00000).view(np.float32)
flat[2048::4096] = np.uint32(0x7FC00001).view(np.float32)
D = np.ascontiguousarray(V[:, :256, :256], dtype=np.float64)

mark("ready")
r = args.repeats
for t in args.threads:
    for name, a in (("", V), ("_nan", N)):
        for ax in (0, 1, 2):
            record(
                f"sum_single_axis{ax}{name}",
                timeit(lambda a=a, ax=ax, t=t: R.sum(a, ax, num_threads=t), repeats=r),
                threads=t,
            )
        record(
            f"sum_single_all{name}",
            timeit(lambda a=a, t=t: R.sum(a, num_threads=t), repeats=r),
            threads=t,
        )
        record(
            f"sum_single_yx{name}",
            timeit(lambda a=a, t=t: R.sum(a, (1, 2), num_threads=t), repeats=r),
            threads=t,
        )
        record(
            f"mean_single_axis0{name}",
            timeit(lambda a=a, t=t: R.mean(a, 0, num_threads=t), repeats=r),
            threads=t,
        )
        record(
            f"var_single_axis0{name}",
            timeit(lambda a=a, t=t: S.var(a, 0, 0, num_threads=t), repeats=r),
            threads=t,
        )
    record(
        "sum_double_all",
        timeit(lambda t=t: R.sum(D, outtype="double", num_threads=t), repeats=r),
        threads=t,
    )
record("cumsum_single_axis0", timeit(lambda: R.cumsum(V, 0), repeats=r))
record("cumsum_single_axis0_nan", timeit(lambda: R.cumsum(N, 0), repeats=r))
