"""The harness against itself: the same build as two variants (``self_check.json``) must agree
within 3%. Each process times, after one warm-up call, an FFT pair of the port
(petakit5d.matlabcompat fftn and ifftn) at the given thread counts and a C-order copy of a
(z, y, x) volume of the same size.

    python self_check.py [--size 256 256 256] [--threads 1 4 24] [--repeats 5]
"""

import argparse

import numpy as np

from petakit5d.matlabcompat.fft import fftn, ifftn
from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("--size", type=int, nargs=3, default=[256, 256, 256])
ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
ap.add_argument("--repeats", type=int, default=5)
a = ap.parse_args()
mark("ready")

rng = np.random.default_rng(0)
x = np.asfortranarray(rng.random(a.size, dtype=np.float32))
for t in a.threads:
    record(
        "fftn_ifftn",
        timeit(lambda t=t: ifftn(fftn(x, threads=t), threads=t), repeats=a.repeats),
        threads=t,
        size="x".join(map(str, a.size)),
    )
c = np.ascontiguousarray(x)
record("copy_c_order", timeit(lambda: c.copy(order="C"), repeats=a.repeats), threads=1)
