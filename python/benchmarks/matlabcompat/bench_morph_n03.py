"""The N03 benchmark operations, Python side (the MATLAB side is bench_morph_n03.m): the same
data in (z, y, x) C order, the same op names and thread counts.

    python bench_morph_n03.py --threads 1 4 24 --repeats 3
"""

import argparse

import numpy as np

from petakit5d.matlabcompat import distance as D
from petakit5d.matlabcompat import label as L
from petakit5d.matlabcompat import morphology as M
from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
ap.add_argument("--repeats", type=int, default=3)
args = ap.parse_args()

# MATLAB's 512 x 1800 x 500 band mask as (z, y, x)
x = np.arange(1, 1801, dtype=np.float64)[None, None, :]
z = np.arange(1, 501, dtype=np.float64)[:, None, None]
band = (x > 0.6 * z + 90) & (x < 0.6 * z + 900)
A = np.ascontiguousarray(np.broadcast_to(band, (500, 512, 1800)))
se = M.strel("disk", 11)

yy, xx = np.meshgrid(np.arange(1, 2049), np.arange(1, 2049), indexing="ij")
S = (yy * 7 + xx * 13) % 997 == 0

y2 = np.arange(1, 129)[None, :, None]
x2 = np.arange(1, 129)[None, None, :]
z2 = np.arange(1, 102)[:, None, None]
P = np.ascontiguousarray(((y2 * 7 + x2 * 13 + z2 * 17) % 4096).astype(np.float32))

y3 = np.arange(1, 513)[None, :, None]
x3 = np.arange(1, 513)[None, None, :]
z3 = np.arange(1, 201)[:, None, None]
V = np.ascontiguousarray((y3 * 7 + x3 * 13 + z3 * 17) % 31 < 12)

mark("ready")
r = args.repeats
for t in args.threads:
    record("imerode_disk11", timeit(lambda t=t: M.imerode(A, se, threads=t), repeats=r), threads=t)
    record("bwdist_2048", timeit(lambda t=t: D.bwdist(S, threads=t), repeats=r), threads=t)
    record("medfilt3_psf", timeit(lambda t=t: M.medfilt3(P, threads=t), repeats=r), threads=t)
    record("bwlabeln_vol", timeit(lambda t=t: L.bwlabeln(V, 26, threads=t), repeats=r), threads=t)
    record(
        "bwconncomp_vol", timeit(lambda t=t: L.bwconncomp(V, 26, threads=t), repeats=r), threads=t
    )
