"""The stitching registration (petakit5d.stitch.registration) on the inputs of
registration_s04_setup.py, at each thread count, on one FFT backend; one warm-up run each.

- pair300: cross_correlation_registration_3d of the two tiles (native zarr);
- demo: the demo's 12 registration tasks (cross_correlation_registration, one per leading tile of
  each of its 4 volumes), as the stitching frame runs them, each into a new folder.

usage: registration_s04.py WORK --case pair300|demo --backend permissive|exact --threads 2 4 24
       --repeats 2
"""

import argparse
import json
import os
import tempfile
from pathlib import Path

from petakit5d import fft as pfft
from petakit5d.stitch import registration as sr
from petakit5d.testing.bench import mark, record, timeit


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--case", default="pair300")
    ap.add_argument("--backend", default="permissive")
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
    ap.add_argument("--repeats", type=int, default=2)
    a = ap.parse_args()
    work = Path(a.work)
    pfft.set_backend(a.backend)
    mark("ready")
    if a.case == "pair300":
        s = json.loads((work / "pair300.json").read_text())
        f1, f2 = str(work / "python" / "t1.zarr"), str(work / "python" / "t2.zarr")
        dorder = (2, 1, 0)
        for t in a.threads:
            out = {}

            def run(t=t, out=out):
                out["r"] = sr.cross_correlation_registration_3d(
                    f1, f2, "", s["c1"], s["c2"], s["co"], s["vs"], dorder,
                    down_sample=(1, 1, 1), max_offset=(10, 20, 20), num_threads=t,
                )  # fmt: skip

            times = timeit(run, repeats=a.repeats, warmup=1)
            rel, mx = out["r"]
            record("pair300", times, threads=t, shift=str(rel.tolist()), peak=float(mx))
        return
    calls = json.loads((work / "demo_python.json").read_text())
    for t in a.threads:

        def before():
            os.chdir(tempfile.mkdtemp(prefix="py_out_", dir=work))

        def run(t=t):
            for c in calls:
                sr.cross_correlation_registration(*c["args"], **c["kwargs"], num_threads=t)

        times = timeit(run, repeats=a.repeats, warmup=1, before=before)
        os.chdir(work)
        record("demo16", times, threads=t, pairs=len(calls))


if __name__ == "__main__":
    main()
