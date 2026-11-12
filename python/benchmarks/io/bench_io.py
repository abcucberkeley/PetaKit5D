"""Read and write times of the cpp-tiff and cpp-zarr Python packages, the way PetaKit5D
uses them: (y, x, z) Fortran-ordered volumes in and out (task IO1). The MATLAB side
is bench_io_mex.m; a driver alternates the two on one node.

    python bench_io.py --input tile.tif --workdir /scratch/x --threads 1 4 24 --repeats 3 --out py.json

Ops (seconds per call, median of --repeats after one warm-up call):
  tiff_read            read_tiff + copy to (y, x, z) Fortran order, copy on one thread
  tiff_read_tconv      the same, the copy split over the threads
  tiff_read_raw        read_tiff alone ((z, y, x) C order, no copy)
  tiff_write           write_tiff of a (y, x, z) Fortran volume (the package copies it to C order)
  tiff_write_direct    the extension's transpose=True entry on the same memory (no copy)
  zarr_read            read_zarr ((y, x, z) Fortran order, no copy)
  zarr_write           write_zarr of a new array, 256^3 chunks, order F, '.' separator
"""

import time

T_START = time.time()

import argparse  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import platform  # noqa: E402
import shutil  # noqa: E402
import socket  # noqa: E402
import statistics  # noqa: E402
import sys  # noqa: E402

import numpy as np  # noqa: E402
import cpptiff  # noqa: E402
import cppzarr  # noqa: E402

T_IMPORTED = time.time()

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tests', 'io', 'bindings'))
import iobind_cases as ic  # noqa: E402


def cpu_model():
    try:
        with open('/proc/cpuinfo') as f:
            return next(line.split(':', 1)[1].strip() for line in f if line.startswith('model name'))
    except (OSError, StopIteration):
        return platform.processor()


def timed(f, repeats, before=None):
    out = []
    for i in range(repeats + 1):
        if before:
            before()
        t0 = time.perf_counter()
        f()
        dt = time.perf_counter() - t0
        if i:
            out.append(dt)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--workdir', required=True)
    p.add_argument('--threads', type=int, nargs='+', default=[1, 4, 24])
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--ops', nargs='+', default=['tiff_read', 'tiff_read_tconv', 'tiff_read_raw', 'tiff_write',
                                                'tiff_write_direct', 'zarr_read', 'zarr_write'])
    p.add_argument('--out', required=True)
    p.add_argument('--label', default='')
    a = p.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    tif_out = os.path.join(a.workdir, 'py_out.tif')
    zarr_out = os.path.join(a.workdir, 'py_out.zarr')
    zarr_in = os.path.join(a.workdir, 'py_in.zarr')
    im = ic.tiff_read_yxz(a.input, num_threads=max(a.threads), convert_threads=max(a.threads))
    cppzarr.write_zarr(zarr_in, im, order='F', chunks=[256, 256, 256], dimension_separator='.')

    def rm_zarr():
        shutil.rmtree(zarr_out, ignore_errors=True)

    rows = []
    for t in a.threads:
        ops = {
            'tiff_read': lambda: ic.tiff_read_yxz(a.input, num_threads=t, convert_threads=1),
            'tiff_read_tconv': lambda: ic.tiff_read_yxz(a.input, num_threads=t, convert_threads=t),
            'tiff_read_raw': lambda: cpptiff.read_tiff(a.input, num_threads=t),
            'tiff_write': lambda: ic.tiff_write_yxz(tif_out, im, 'lzw', t),
            'tiff_write_direct': lambda: ic.tiff_write_yxz_direct(tif_out, im, 'lzw', t),
            'zarr_read': lambda: cppzarr.read_zarr(zarr_in, num_threads=t),
            'zarr_write': lambda: cppzarr.write_zarr(zarr_out, im, order='F', chunks=[256, 256, 256],
                                                     dimension_separator='.', num_threads=t),
        }
        for op in a.ops:
            ts = timed(ops[op], a.repeats, rm_zarr if op == 'zarr_write' else None)
            rows.append(dict(impl='python', op=op, threads=t, times=ts, median=statistics.median(ts)))
            print(f'python {op:18s} t={t:2d} median {statistics.median(ts):.3f} s  {["%.3f" % x for x in ts]}', flush=True)
    rm_zarr()
    shutil.rmtree(zarr_in, ignore_errors=True)
    if os.path.exists(tif_out):
        os.remove(tif_out)
    res = dict(label=a.label, impl='python', host=socket.gethostname(), cpu=cpu_model(), affinity=len(os.sched_getaffinity(0)),
               versions=ic.versions(), shape=list(im.shape), dtype=str(im.dtype), input=os.path.basename(a.input),
               t_start=T_START, t_imported=T_IMPORTED, t_end=time.time(), rows=rows)
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
