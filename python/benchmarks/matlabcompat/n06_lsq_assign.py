"""N06 benchmark: the 900-tile assignment problems through the caller's systems (built as the
MATLAB code builds them, the test module's replica) and petakit5d.matlabcompat.lsq; the same
code at every thread count (it is serial).

usage: n06_lsq_assign.py <python/ tree> <work folder> [--threads N ...] [--repeats R]
"""

import argparse
import importlib.util

from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("src")
ap.add_argument("work")
ap.add_argument("--threads", type=int, nargs="+", default=[1])
ap.add_argument("--repeats", type=int, default=5)
a = ap.parse_args()

import scipy.io as sio  # noqa: E402

spec = importlib.util.spec_from_file_location("test_lsq", f"{a.src}/tests/matlabcompat/test_lsq.py")
tl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tl)
mark("ready")
for name in ("grid", "free", "tight"):
    m = sio.loadmat(f"{a.work}/n06_p900_{name}.mat", squeeze_me=False)
    p = {
        "nF": float(m["nF"][0, 0]),
        "mx": m["mx"],
        "ash": m["ash"],
        "pairs": m["pairs"],
        "mas": m["mas"],
        "thr": float(m["thr"][0, 0]),
        "tileIdx": m["tileIdx"],
        "aw": m["aw"].ravel(),
    }
    for t in a.threads:
        times = timeit(lambda p=p: tl.assign(p, "grid"), repeats=a.repeats)
        record(f"assign_900_{name}", times, threads=t)
