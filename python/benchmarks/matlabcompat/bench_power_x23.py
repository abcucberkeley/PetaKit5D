"""matlabcompat.elementary.power and the resize kernels' pow before and after their subnormal
results took fdlibm's scalbn instead of the C runtime's ldexp (power_x23.json runs this script
on two source trees and their builds): X06's million positive doubles with scalar and array
exponents (some results subnormal), a million results that are all subnormal, many calls with
one pair (the scalar path, with ordinary and with subnormal results), and the kernels' pow of a
million cubes (ordinary and subnormal), their cubic kernel and imresize3's cubic resampling.

    python bench_power_x23.py --repeats 3
"""

import argparse

import numpy as np

from petakit5d._ext import load
from petakit5d.matlabcompat import elementary as E
from petakit5d.matlabcompat import resize as R
from petakit5d.testing.bench import mark, record, timeit
from petakit5d.testing.synth import CounterRNG

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=3)
ap.add_argument("--n", type=int, default=1_000_000)
args = ap.parse_args()

n = args.n
k = np.arange(1, n + 1, dtype=np.int64)
# X06's inputs: integers of 21 bits times 2^(-30..30), exponents from -31.25 to 31.25
X = np.ldexp((k * 7919 % 2**21 + 1).astype(np.float64), (k * 31 % 61 - 30 - 21).astype(np.int32))
Y = (k * 13 % 4001 - 2000) / 64.0
# 2^y for y in [-1075, -1022): every result subnormal (or 0)
YS = -1075 + (k % 3392) / 64.0
pairs = [(float(a), float(b)) for a, b in zip(X[:10000], Y[:10000], strict=True)]
sub_pairs = [(2.0, float(b)) for b in YS[:10000]]
# the kernels' pow: cubes of distances (ordinary results) and of tiny values (subnormal results)
D = (k % 4001) / 1000.0 - 2.0
AD = np.abs(D)
TINY = np.ldexp(1 + (k % 8192) * 2.0**-13, (-342 - k % 17).astype(np.int32))
rs = load("resize")
out = np.empty(n)
V = CounterRNG(0, "x23/bench").uniform((48, 160, 144), np.float32)

mark("ready")
r = args.repeats
with np.errstate(all="ignore"):
    record("power_scalar_exp_2p5", timeit(lambda: E.power(X, 2.5), repeats=r))
    record("power_array_exp", timeit(lambda: E.power(X, Y), repeats=r))
    record("power_subnormal_results", timeit(lambda: E.power(2.0, YS), repeats=r))
    record("power_10k_scalar_calls", timeit(lambda: [E.power(a, b) for a, b in pairs], repeats=r))
    record(
        "power_10k_scalar_subnormal",
        timeit(lambda: [E.power(a, b) for a, b in sub_pairs], repeats=r),
    )
    record("kernel_pow_cubes", timeit(lambda: rs.fdlibm_pow(AD, 3.0, out), repeats=r))
    record("kernel_pow_subnormal_cubes", timeit(lambda: rs.fdlibm_pow(TINY, 3.0, out), repeats=r))
    record("kernel_cubic", timeit(lambda: rs.cubic(D, out), repeats=r))
    for t in (1, 24):
        record(
            "imresize3_cubic",
            timeit(lambda t=t: R.imresize3(V, (40, 200, 90), "cubic", num_threads=t), repeats=r),
            threads=t,
        )
