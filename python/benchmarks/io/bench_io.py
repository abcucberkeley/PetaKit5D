"""Read and write times of petakit5d.io on the port's (z, y, x) arrays (tasks IO1, IO3). The
MATLAB side is bench_io_mex.m (readtiff, writetiff, readzarr, writezarr of MATLAB's own
layout); a driver alternates the two on one node, both under ``numactl --interleave=all`` at
24 threads (the benchmark standard).

    python bench_io.py --input tile.tif --workdir /scratch/x --threads 1 4 24 --repeats 3 \
        --out py.json

Each op runs once untimed and then ``--repeats`` times; ``second`` is the first timed call,
the second call in the process, which both sides compare (a first read into new memory pays
its page faults). The result of a call is freed after the timer stops, as MATLAB's ``clear``
after ``toc``. Every call records its minor page faults and the process's AnonHugePages.

Ops:
  tiff_read          read_tiff (cpp-tiff's (z, y, x) array, no copy)
  tiff_read_reuse    read_tiff into one ReadBuffer over the calls
  tiff_write         write_tiff (LZW) of the volume
  zarr_read          read_zarr of a native array (256^3 chunks)
  zarr_read_reuse    the same into one ReadBuffer over the calls
  zarr_read_mlayout  read_zarr of MATLAB's (y, x, z) order-F array of the same data
  zarr_write         write_zarr of a new native array (blockSize 256^3)
"""

import time

T_START = time.time()

import argparse  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import platform  # noqa: E402
import resource  # noqa: E402
import shutil  # noqa: E402
import socket  # noqa: E402
import statistics  # noqa: E402

import numpy as np  # noqa: E402

import petakit5d.io as pio  # noqa: E402

T_IMPORTED = time.time()

OPS = [
    "tiff_read",
    "tiff_read_reuse",
    "tiff_write",
    "zarr_read",
    "zarr_read_reuse",
    "zarr_read_mlayout",
    "zarr_write",
]


def cpu_model():
    try:
        with open("/proc/cpuinfo") as f:
            return next(
                line.split(":", 1)[1].strip() for line in f if line.startswith("model name")
            )
    except (OSError, StopIteration):
        return platform.processor()


def huge_pages_kb() -> int:
    """the process's AnonHugePages (kB), -1 where /proc has no smaps_rollup"""
    try:
        with open("/proc/self/smaps_rollup") as f:
            for line in f:
                if line.startswith("AnonHugePages:"):
                    return int(line.split()[1])
    except OSError:
        pass
    return -1


def timed(f, repeats, before=None):
    """``repeats`` timed calls after one untimed: (times, records of every call). The result
    is freed after the timer stops."""
    times, calls = [], []
    for i in range(repeats + 1):
        if before:
            before()
        f0 = resource.getrusage(resource.RUSAGE_SELF).ru_minflt
        t0 = time.perf_counter()
        res = f()
        dt = time.perf_counter() - t0
        f1 = resource.getrusage(resource.RUSAGE_SELF).ru_minflt
        calls.append({"t": dt, "minflt": f1 - f0, "huge_kb": huge_pages_kb()})
        del res
        if i:
            times.append(dt)
    return times, calls


def _ops(t, tif_in, vol, tif_out, zarr_in, zarr_m, zarr_out, bufs):
    """The timed calls for t threads."""
    return {
        "tiff_read": lambda: pio.read_tiff(tif_in, num_threads=t),
        "tiff_read_reuse": lambda: pio.read_tiff(tif_in, num_threads=t, out=bufs["tiff"]),
        "tiff_write": lambda: pio.write_tiff(vol, tif_out, num_threads=t),
        "zarr_read": lambda: pio.read_zarr(zarr_in, num_threads=t),
        "zarr_read_reuse": lambda: pio.read_zarr(zarr_in, num_threads=t, out=bufs["zarr"]),
        "zarr_read_mlayout": lambda: pio.read_zarr(zarr_m, num_threads=t),
        "zarr_write": lambda: pio.write_zarr(
            vol, zarr_out, block_size=(256, 256, 256), num_threads=t
        ),
    }


def main():
    import cppzarr

    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--workdir", required=True)
    p.add_argument("--threads", type=int, nargs="+", default=[1, 4, 24])
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--ops", nargs="+", default=OPS)
    p.add_argument("--out", required=True)
    p.add_argument("--label", default="")
    a = p.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    tif_out = os.path.join(a.workdir, "py_out.tif")
    zarr_out = os.path.join(a.workdir, "py_out.zarr")
    zarr_in = os.path.join(a.workdir, "py_in.zarr")
    zarr_m = os.path.join(a.workdir, "py_in_matlab.zarr")
    vol = pio.read_tiff(a.input)
    pio.write_zarr(vol, zarr_in, block_size=(256, 256, 256))
    # MATLAB's array of the same data, as writezarr writes it (cpp-zarr's order-F write)
    m = np.transpose(vol, (1, 2, 0))
    cppzarr.write_zarr(zarr_m, m, order="F", chunks=[256, 256, 256], dimension_separator=".")
    bufs = {"tiff": pio.ReadBuffer(), "zarr": pio.ReadBuffer()}

    def rm_zarr():
        shutil.rmtree(zarr_out, ignore_errors=True)

    rows = []
    for t in a.threads:
        ops = _ops(t, a.input, vol, tif_out, zarr_in, zarr_m, zarr_out, bufs)
        for op in a.ops:
            ts, calls = timed(ops[op], a.repeats, rm_zarr if op == "zarr_write" else None)
            med = statistics.median(ts)
            rows.append(
                {
                    "impl": "python",
                    "op": op,
                    "threads": t,
                    "times": ts,
                    "median": med,
                    "second": ts[0],
                    "calls": calls,
                }
            )
            each = [f"{x:.3f}" for x in ts]
            flt = [c["minflt"] for c in calls]
            print(
                f"python {op:18s} t={t:2d} median {med:.3f} s  second {ts[0]:.3f}  {each}  "
                f"minflt {flt}",
                flush=True,
            )
    rm_zarr()
    for z in (zarr_in, zarr_m):
        shutil.rmtree(z, ignore_errors=True)
    if os.path.exists(tif_out):
        os.remove(tif_out)
    res = {
        "label": a.label,
        "impl": "python",
        "host": socket.gethostname(),
        "cpu": cpu_model(),
        "affinity": len(os.sched_getaffinity(0)),
        "shape": list(vol.shape),
        "dtype": str(vol.dtype),
        "input": os.path.basename(a.input),
        "t_start": T_START,
        "t_imported": T_IMPORTED,
        "t_end": time.time(),
        "rows": rows,
    }
    with open(a.out, "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
