"""W03 trial: the detection fits timed in one process (W02's detection family: 3000 windows of 9^3,
the fitGaussian3D mex in modes xyzAsrc, xyzAc and Ac, and the CD_Newton fit), at the given thread
counts, through petakit5d.testing.bench (record, timeit), so that the harness interleaves the
three builds' processes (MinGW-w64 with the UCRT's erf, erfc and log, MinGW-w64 with libmingwex's,
MSVC).

    python bench_fits.py --threads 1 4 --repeats 3
"""

import argparse
import faulthandler

import numpy as np

faulthandler.dump_traceback_later(1800, exit=True)

import petakit5d.core.threads as th  # noqa: E402
from petakit5d._ext import load  # noqa: E402
from petakit5d.testing.bench import info, mark, record, timeit  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--threads", type=int, nargs="+", default=[1, 4])
ap.add_argument("--repeats", type=int, default=3)
args = ap.parse_args()
r = args.repeats
rng = np.random.default_rng(20261121)


def set_threads(t: int) -> None:
    th._affinity_now = lambda: t
    th.refresh_cpus()


from petakit5d.detection import fit  # noqa: E402

k = 3000
w = 9
z, y, x = np.meshgrid(*(np.arange(w, dtype=np.float64),) * 3, indexing="ij")
windows, p0 = [], []
for _ in range(k):
    c = 4 + rng.uniform(-0.6, 0.6, 3)
    a, s, sz = rng.uniform(200, 900), rng.uniform(1.1, 1.6), rng.uniform(1.4, 2.2)
    g = a * np.exp(
        -((x - c[2]) ** 2 + (y - c[1]) ** 2) / (2 * s * s) - (z - c[0]) ** 2 / (2 * sz * sz)
    )
    windows.append(np.ascontiguousarray(g + 100 + rng.normal(0, 10, g.shape)))
    p0.append([4.0, 4.0, 4.0, a * 0.8, 1.3, 1.8, 100.0])
p0 = np.array(p0)
info(windows=k, size=w, compiler=load("core").build_info()["compiler"])
mark("ready")
for t in args.threads:
    set_threads(t)
    for mode in ("xyzAsrc", "xyzAc", "Ac"):
        record(
            f"fit_mex_{mode}_{k}",
            timeit(lambda t=t, m=mode: fit.fit_gaussian3d_mex(windows, p0, m, num_threads=t), r),
            threads=t,
        )
    record(f"fit_cd_{k}", timeit(lambda t=t: fit.fit_cd(windows, p0, num_threads=t), r), threads=t)
