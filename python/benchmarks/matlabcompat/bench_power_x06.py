"""matlabcompat.elementary.power before and after its element-wise form (power_x06.json runs this
script on two source trees): a million positive doubles with scalar and array exponents, a
million singles, and many calls with one pair (the scalar path).

    python bench_power_x06.py --repeats 2
"""

import argparse

import numpy as np

from petakit5d.matlabcompat import elementary as E
from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("--repeats", type=int, default=2)
ap.add_argument("--n", type=int, default=1_000_000)
args = ap.parse_args()

n = args.n
k = np.arange(1, n + 1, dtype=np.int64)
# integers of 21 bits times 2^(-30..30): every magnitude of a window or a weight
X = np.ldexp((k * 7919 % 2**21 + 1).astype(np.float64), (k * 31 % 61 - 30 - 21).astype(np.int32))
Y = (k * 13 % 4001 - 2000) / 64.0
S = X.astype(np.float32)
pairs = [(float(a), float(b)) for a, b in zip(X[:10000], Y[:10000], strict=True)]

mark("ready")
r = args.repeats
record("power_scalar_exp_2p5", timeit(lambda: E.power(X, 2.5), repeats=r))
record("power_scalar_exp_third", timeit(lambda: E.power(X, 1 / 3), repeats=r))
record("power_array_exp", timeit(lambda: E.power(X, Y), repeats=r))
record("power_single_2p5", timeit(lambda: E.power(S, 2.5), repeats=r))
record("power_10k_scalar_calls", timeit(lambda: [E.power(a, b) for a, b in pairs], repeats=r))
