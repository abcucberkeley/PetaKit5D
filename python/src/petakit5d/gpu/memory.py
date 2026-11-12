"""GPU memory between tasks: CuPy's pool, free memory, peak use.

Policy: inside a task, arrays come from CuPy's memory pool (freed arrays are reused at once).
Between tasks of one process, task_scope returns the pool's cached blocks to the GPU, so the
next task's free-memory check sees what is really free and other processes on a shared GPU
get it back; the CUDA context (the costly part to make) stays, and so do the cached FFT plans
unless asked (consecutive tasks usually have the same size). All of it acts on the current GPU.
"""

from __future__ import annotations

import sys
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any

from petakit5d.gpu._cupy import import_cupy

__all__ = ["PeakUse", "available_bytes", "pool_bytes", "release_cached", "task_scope", "track_peak"]


def pool_bytes() -> dict[str, int]:
    """Our pool on the current GPU: 'used' by arrays, 'cached' (free blocks it keeps), 'fft_plans'
    (the work areas of the cached FFT plans)."""
    cp = import_cupy()
    from petakit5d.gpu import fft

    pool = cp.get_default_memory_pool()
    return {
        "used": int(pool.used_bytes()),
        "cached": int(pool.free_bytes()),
        "fft_plans": fft.cached_plan_bytes(),
    }


def available_bytes() -> int:
    """Free memory of the current GPU plus the blocks our pool caches (what new arrays can take)."""
    cp = import_cupy()
    free, _ = cp.cuda.runtime.memGetInfo()
    return int(free) + int(cp.get_default_memory_pool().free_bytes())


def release_cached(*, plans: bool = False) -> int:
    """Return the pool's cached blocks (and pinned host blocks) to the system; with plans, also
    drop the cached FFT plans first (ours and cupy.fft's). Returns the bytes given back."""
    cp = import_cupy()
    from petakit5d.gpu import fft

    pool = cp.get_default_memory_pool()
    before = pool.total_bytes()
    if plans:
        fft.clear_plan_cache()
        cp.fft.config.get_plan_cache().clear()
    pool.free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    return int(before - pool.total_bytes())


@contextmanager
def task_scope(*, keep_plans: bool = True) -> Iterator[None]:
    """Run one task; afterwards give the pool's cached blocks back (see the module)."""
    try:
        yield
    finally:
        if "cupy" in sys.modules:  # (a task that never used the GPU doesn't start CUDA here)
            release_cached(plans=not keep_plans)


class PeakUse:
    """Peak bytes of our pool's arrays while tracked (FFT work areas included; the CUDA context and
    cuFFT's internal buffers are not)."""

    def __init__(self, start: int) -> None:
        self.start = start
        self.used = start
        self.peak = start

    @property
    def peak_above_start(self) -> int:
        return self.peak - self.start


@contextmanager
def track_peak() -> Iterator[PeakUse]:
    """with track_peak() as p: ...; p.peak is the largest pool use on the current GPU meanwhile."""
    cp = import_cupy()
    from cupy.cuda import memory_hook

    use = PeakUse(int(cp.get_default_memory_pool().used_bytes()))

    class _Hook(memory_hook.MemoryHook):
        name = "petakit5d_peak"

        def malloc_postprocess(self, **kw: Any) -> None:
            use.used += int(kw["mem_size"])
            use.peak = max(use.peak, use.used)

        def free_postprocess(self, **kw: Any) -> None:
            use.used -= int(kw["mem_size"])

    with _Hook():
        yield use
