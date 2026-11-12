"""Check the GPU installation: python -m petakit5d.gpu.selftest [--json FILE] [--size NY NX NZ].

Reports the CuPy, CUDA, NVRTC and cuFFT versions and libraries, the GPUs and the one chosen
(with its free memory), runs a kernel compiled at run time and an FFT round trip against
NumPy, and times each step. Exit status 0 when all checks pass, 1 when one fails, 2 without a
usable GPU.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any

_SOURCE = r"""
template <typename T>
__global__ void add_scaled(T* x, const T* y, const T a, const long long n) {
    const long long step = (long long)blockDim.x * gridDim.x;
    for (long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x; i < n; i += step)
        x[i] = x[i] + a * y[i];
}
"""


def run(size: tuple[int, ...] = (96, 80, 72)) -> dict[str, Any]:
    """The checks as a dict ('ok': all passed); see the module."""
    t = {}
    r: dict[str, Any] = {
        "ok": False,
        "times": t,
        "env": {
            k: os.environ.get(k)
            for k in (
                "CUDA_VISIBLE_DEVICES",
                "PETAKIT5D_CACHE_DIR",
                "PETAKIT5D_GPU_TARGET",
                "PETAKIT5D_CUFFT_LIBRARY",
            )
        },
    }
    t0 = time.perf_counter()
    from petakit5d.gpu import device, fft, kernels

    t["import_petakit5d_gpu"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    count = device.gpu_count()
    t["import_cupy_and_count"] = time.perf_counter() - t0
    r["count"] = {
        "usable": count.usable,
        "found": count.found,
        "reason": count.reason,
        "devices": [d._asdict() for d in count.devices],
    }
    if count.usable == 0:
        r["error"] = count.reason or "No GPU is found."
        return r
    import cupy as cp
    import numpy as np

    r["versions"] = {
        "cupy": cp.__version__,
        "numpy": np.__version__,
        "cuda_driver": cp.cuda.runtime.driverGetVersion(),
        "cuda_runtime": cp.cuda.runtime.runtimeGetVersion(),
        "nvrtc": ".".join(map(str, kernels.nvrtc_version())),
        "nvrtc_archs": kernels.nvrtc_archs(),
    }
    t0 = time.perf_counter()
    dev = device.select_device()
    t["context"] = time.perf_counter() - t0
    r["device"] = dev._asdict()
    # a kernel compiled at run time: x + a * y in float32 and float64, exact for these values
    t0 = time.perf_counter()
    mod = kernels.load_source(
        _SOURCE, name="selftest.cu", name_expressions=("add_scaled<float>", "add_scaled<double>")
    )
    t["kernel_load"] = time.perf_counter() - t0
    r["kernel"] = mod.info._asdict()
    n = 1 << 20
    ok = True
    for dt, name in ((np.float32, "add_scaled<float>"), (np.float64, "add_scaled<double>")):
        x = cp.arange(n, dtype=dt)
        y = cp.full(n, 2, dtype=dt)
        mod[name]((n + 255) // 256, 256, (x, y, dt(3), np.int64(n)))
        ok = ok and bool((x.get() == np.arange(n, dtype=dt) + dt(6)).all())
    cp.cuda.Device().synchronize()
    t["kernel_run"] = time.perf_counter() - t0 - t["kernel_load"]
    r["kernel_ok"] = ok
    # FFT round trip against NumPy (the half spectrum along y, as numpy's rfftn over z, x, y)
    t0 = time.perf_counter()
    info = fft.cufft_info()
    t["cufft_load"] = time.perf_counter() - t0
    r["cufft"] = info
    rng = np.random.default_rng(0)
    a = np.asfortranarray(rng.random(size, dtype=np.float32))
    t0 = time.perf_counter()
    plan = fft.get_plan(size)
    t["fft_plan"] = time.perf_counter() - t0
    g = cp.asarray(a, order="F")
    t0 = time.perf_counter()
    s = plan.forward(g)
    back = plan.inverse(s.copy(order="F"))
    cp.cuda.Device().synchronize()
    t["fft_run"] = time.perf_counter() - t0
    ref = np.fft.rfftn(a.astype(np.float64), axes=tuple(range(a.ndim - 1, -1, -1)))
    sg = s.get()
    fwd_err = float(np.abs(sg - ref).max() / np.abs(ref).max())
    inv_err = float(np.abs(back.get() / a.size - a).max() / np.abs(a).max())
    r["fft"] = {
        "size": list(size),
        "half_shape": list(plan.half_shape),
        "plan_dims": list(plan.dims),
        "work_bytes": plan.work_bytes,
        "forward_rel_err": fwd_err,
        "roundtrip_rel_err": inv_err,
    }
    r["fft_ok"] = s.shape == ref.shape and fwd_err < 1e-5 and inv_err < 1e-5
    free, total = cp.cuda.runtime.memGetInfo()
    r["memory_after"] = {"free": int(free), "total": int(total)}
    r["ok"] = r["kernel_ok"] and r["fft_ok"]
    return r


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m petakit5d.gpu.selftest", description=__doc__.splitlines()[0]
    )
    ap.add_argument("--json", help="also write the results to this file")
    ap.add_argument("--size", type=int, nargs="+", default=[96, 80, 72], help="FFT size (y x z)")
    args = ap.parse_args(argv)
    t0 = time.perf_counter()
    try:
        r = run(tuple(args.size))
    except Exception as e:  # (report any failure)
        r = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    r["total_seconds"] = time.perf_counter() - t0
    text = json.dumps(r, indent=1, default=str)
    print(text)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            f.write(text)
    if r.get("ok"):
        return 0
    return 2 if "count" in r and r["count"]["usable"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
