"""Which GPU this process may use, and why none.

Ports utils/usable_gpu_num.m and decon/decon_gpu_num.m with the "required" and "auto" use of
the parked MATLAB change (useGPU true: a usable GPU is required, else decon:noUsableGPU; []:
a GPU if one is usable, else the CPU with a warning that gives the reason).

The GPUs are those CUDA lists (CUDA_VISIBLE_DEVICES, which Slurm sets for the GPUs it gave the
job). Unlike MATLAB, a process in a Slurm job that has no GPUs of its own and no
CUDA_VISIBLE_DEVICES uses none: on nodes that don't hide their GPUs from such jobs, the GPUs
CUDA lists there belong to other jobs. Setting CUDA_VISIBLE_DEVICES uses them anyway.

A GPU counts as usable when NVRTC compiles for it (GPUs newer than NVRTC knows run PTX, see
gpu.kernels; MATLAB needs forward compatibility for them) and it isn't in prohibited compute
mode. Selecting one creates its CUDA context and checks its free memory: GPUs are not
necessarily isolated between jobs, so memory another process holds is seen there.
"""

from __future__ import annotations

import importlib
import os
import threading
import warnings
from collections.abc import Mapping
from typing import NamedTuple

from petakit5d.gpu._cupy import GPUUnavailableError, cupy_or_none, import_cupy

__all__ = [
    "GPUCount",
    "GPUFallbackWarning",
    "GPUInfo",
    "GPUUnavailableError",
    "NoUsableGPUError",
    "gpu_count",
    "gpu_for",
    "job_gpu_reason",
    "select_device",
    "usable_gpu_count",
    "use_gpu_mode",
]

_SLURM_GPU_VARS = ("SLURM_JOB_GPUS", "SLURM_STEP_GPUS", "SLURM_GPUS_ON_NODE")
_CUDA_ERROR_NO_DEVICE = 100
_COMPUTE_MODE_PROHIBITED = 2

_lock = threading.Lock()
_count: GPUCount | None = None


class NoUsableGPUError(GPUUnavailableError):
    """A GPU is required and none is usable (MATLAB: decon:noUsableGPU)."""


class GPUFallbackWarning(UserWarning):
    """A GPU would be used if one were usable: the work runs on the CPU."""


class GPUInfo(NamedTuple):
    """A GPU as CUDA lists it in this process."""

    index: int  # CUDA's device number in this process (after CUDA_VISIBLE_DEVICES)
    name: str
    compute_capability: int  # major * 10 + minor (100 for 10.0)
    total_bytes: int
    pci_bus_id: str
    free_bytes: int | None = None  # free memory when selected, our pool's cached blocks included

    def describe(self) -> str:
        return (
            f"GPU {self.index} ({self.name}, compute capability {self.compute_capability // 10}."
            f"{self.compute_capability % 10})"
        )


class GPUCount(NamedTuple):
    """GPUs as MATLAB's usable_gpu_num counts them (usable, found), and why none is usable."""

    usable: int
    found: int
    reason: str
    devices: tuple[GPUInfo, ...] = ()


def job_gpu_reason(env: Mapping[str, str], found: int) -> str:
    """Why the GPUs CUDA lists are not this process's to use, '' if they are (see the module)."""
    if found == 0 or "CUDA_VISIBLE_DEVICES" in env or not env.get("SLURM_JOB_ID"):
        return ""
    # GPU indices ('0' is GPU 0) or a count ('0' is none)
    if (
        env.get("SLURM_JOB_GPUS", "").strip()
        or env.get("SLURM_STEP_GPUS", "").strip()
        or env.get("SLURM_GPUS_ON_NODE", "").strip() not in ("", "0")
    ):
        return ""
    return (
        f"The Slurm job {env['SLURM_JOB_ID']} has no GPU (none of {', '.join(_SLURM_GPU_VARS)} "
        f"is set) and CUDA_VISIBLE_DEVICES is not set, so the {found} GPUs CUDA lists belong to "
        "other jobs: ask for a GPU (--gres=gpu:1), or set CUDA_VISIBLE_DEVICES to use them anyway."
    )


def gpu_count(*, refresh: bool = False) -> GPUCount:
    """The GPUs this process may use (counted once per process unless refresh)."""
    global _count
    with _lock:
        if _count is None or refresh:
            _count = _count_gpus(os.environ)
        return _count


def usable_gpu_count() -> tuple[int, str, int]:
    """(usable, reason, found) as MATLAB's [gpuNum, reason, gpuAll] = usable_gpu_num()."""
    c = gpu_count()
    return c.usable, c.reason, c.found


def _count_gpus(env: Mapping[str, str]) -> GPUCount:
    cp, reason = cupy_or_none()
    if cp is None:
        return GPUCount(0, 0, reason)
    rt = cp.cuda.runtime
    try:
        found = rt.getDeviceCount()
    except rt.CUDARuntimeError as e:
        if getattr(e, "status", None) == _CUDA_ERROR_NO_DEVICE:
            return GPUCount(0, 0, "")
        return GPUCount(0, 0, f"Counting the GPUs failed: {e}")
    except Exception as e:  # (CUDA libraries that fail to load)
        return GPUCount(0, 0, f"Counting the GPUs failed: {type(e).__name__}: {e}")
    if found == 0:
        return GPUCount(0, 0, "")
    why = job_gpu_reason(env, found)
    if why:
        return GPUCount(0, found, why)
    try:
        from petakit5d.gpu.kernels import nvrtc_archs

        oldest = nvrtc_archs()[0]
    except Exception as e:
        return GPUCount(
            0, found, f"The CUDA compiler (NVRTC) failed to load: {type(e).__name__}: {e}"
        )
    devices, reasons = [], []
    for i in range(found):
        try:
            p = rt.getDeviceProperties(i)
        except Exception as e:
            reasons.append(f"Reading the properties of GPU {i} failed: {e}")
            continue
        name = p["name"].decode() if isinstance(p["name"], bytes) else str(p["name"])
        d = GPUInfo(
            i,
            name,
            int(p["major"]) * 10 + int(p["minor"]),
            int(p["totalGlobalMem"]),
            f"{int(p.get('pciDomainID', 0)):04x}:{int(p['pciBusID']):02x}:"
            f"{int(p.get('pciDeviceID', 0)):02x}.0",
        )
        if d.compute_capability < oldest:
            reasons.append(
                f"{d.describe()} is older than the oldest GPU this CUDA compiler supports "
                f"({oldest // 10}.{oldest % 10})."
            )
        elif int(p.get("computeMode", 0)) == _COMPUTE_MODE_PROHIBITED:
            reasons.append(f"{d.describe()} is in prohibited compute mode.")
        else:
            devices.append(d)
    return GPUCount(len(devices), found, "" if devices else " ".join(reasons), tuple(devices))


def _gib(n: int) -> str:
    return f"{n / 2**30:.2f} GiB"


def select_device(index: int | None = None, need_bytes: int = 0) -> GPUInfo:
    """Make a usable GPU current (in this thread) and return it with its free memory.

    index: CUDA's device number; None: the first usable GPU in CUDA's order with need_bytes
    free (memory our pool caches counts as free). Raises GPUUnavailableError with the reason.
    """
    c = gpu_count()
    if c.usable == 0:
        raise GPUUnavailableError(c.reason or "No GPU is found.")
    cands = [d for d in c.devices if index is None or d.index == index]
    if not cands:
        raise GPUUnavailableError(
            f"GPU {index} is not a usable GPU of this process "
            f"(usable: {', '.join(str(d.index) for d in c.devices)})."
        )
    cp = import_cupy()
    reasons = []
    for d in cands:
        try:
            cp.cuda.Device(d.index).use()
            free, total = cp.cuda.runtime.memGetInfo()  # (creates the context)
            free += cp.get_default_memory_pool().free_bytes()
        except Exception as e:  # (busy, broken, out of memory for a context)
            reasons.append(f"Selecting GPU {d.index} failed: {' '.join(str(e).split())}")
            continue
        if free < need_bytes:
            reasons.append(
                f"{d.describe()} has {_gib(free)} free of {_gib(total)}, and {_gib(need_bytes)} "
                "are needed."
            )
            continue
        _cupy_ptx_for_new_gpus(d.compute_capability)
        return d._replace(free_bytes=int(free))
    raise GPUUnavailableError(" ".join(reasons))


def _cupy_ptx_for_new_gpus(cc: int) -> None:
    # CuPy (14.2) compiles its own kernels as a cubin for the newest architecture NVRTC knows,
    # which a newer GPU can't run; PTX lets the driver compile them (CUPY_COMPILE_WITH_PTX)
    from petakit5d.gpu.kernels import nvrtc_archs

    if cc > nvrtc_archs()[-1]:
        compiler = importlib.import_module("cupy.cuda.compiler")
        if getattr(compiler, "_use_ptx", True) is False:
            compiler._use_ptx = True


def use_gpu_mode(use_gpu: bool | None) -> str:
    """'required' for True, 'auto' for None (MATLAB's []), 'off' for False (the parked change)."""
    if use_gpu is None:
        return "auto"
    if not isinstance(use_gpu, bool):
        raise TypeError(f"use_gpu must be True, False or None, not {use_gpu!r}")
    return "required" if use_gpu else "off"


def gpu_for(
    use_gpu: bool | None,
    *,
    purpose: str = "the computation",
    need_bytes: int = 0,
    index: int | None = None,
) -> GPUInfo | None:
    """The GPU to run on, made current, or None for the CPU.

    use_gpu False: None. None ("auto"): a usable GPU, else None with a GPUFallbackWarning that
    gives the reason. True ("required"): a usable GPU, else NoUsableGPUError. purpose names the
    work in the messages ('the deconvolution'); need_bytes and index as for select_device.
    """
    mode = use_gpu_mode(use_gpu)
    if mode == "off":
        return None
    try:
        return select_device(index, need_bytes)
    except GPUUnavailableError as e:
        reason = str(e)
    if reason == "No GPU is found.":
        reason = ""
    msg = (
        f"No GPU is available ({gpu_count().found} found, "
        f"CUDA_VISIBLE_DEVICES='{os.environ.get('CUDA_VISIBLE_DEVICES', '')}')"
    )
    tail = f" {reason}" if reason else ""
    if mode == "required":
        raise NoUsableGPUError(f"{msg}, but {purpose} is to run on a GPU (use_gpu=True).{tail}")
    warnings.warn(f"{msg}, so {purpose} runs on the CPU.{tail}", GPUFallbackWarning, stacklevel=2)
    return None
