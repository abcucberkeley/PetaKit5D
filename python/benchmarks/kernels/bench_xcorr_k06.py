"""normxcorr3_max_shift (petakit5d.kernels.xcorr) on the pairs of xcorr_k06_setup.py, at each
thread count, on one FFT backend; one warm-up call per pair and thread count.

usage: bench_xcorr_k06.py WORK --backend permissive|exact --threads 1 4 24 --repeats 3
"""

import argparse
from pathlib import Path

import numpy as np

from petakit5d.io import read_tiff
from petakit5d.kernels import xcorr
from petakit5d.testing.bench import mark, record, timeit

# MATLAB's maxShifts ([y x z] rows of lower and upper bounds) as (z, y, x)
PAIRS = {
    "ovl300": np.array([[0, 0, 0], [20, 40, 40]], dtype=np.float64),
    "demo2p": np.array([[0, -98, -98], [20, 302, 302]], dtype=np.float64),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--backend", default="permissive")
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
    ap.add_argument("--repeats", type=int, default=3)
    a = ap.parse_args()
    mark("ready")
    for name, ms in PAIRS.items():
        T = read_tiff(str(Path(a.work) / f"{name}_T.tif"))
        A = read_tiff(str(Path(a.work) / f"{name}_A.tif"))
        for t in a.threads:
            out = {}

            def run(t=t, ms=ms, T=T, A=A, out=out):
                out["r"] = xcorr.normxcorr3_max_shift(T, A, ms, threads=t, backend=a.backend)

            times = timeit(run, repeats=a.repeats, warmup=1)
            off, mx, _ = out["r"]
            record(f"max_shift_{name}", times, threads=t, offset=str(off), peak=float(mx))


if __name__ == "__main__":
    main()
