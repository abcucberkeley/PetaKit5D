"""CuPy, imported when a GPU path first needs it, with the reason when it can't be."""

from __future__ import annotations

import importlib
import threading
from types import ModuleType

_lock = threading.Lock()
_state: tuple[ModuleType | None, str] | None = None


class GPUUnavailableError(RuntimeError):
    """No GPU can be used: CuPy is missing or broken, or no device is usable."""


def cupy_or_none() -> tuple[ModuleType | None, str]:
    """CuPy and '', or None and why it can't be imported (tried once per process)."""
    global _state
    with _lock:
        if _state is None:
            try:
                cp = importlib.import_module("cupy")
            except ModuleNotFoundError as e:
                if e.name in ("cupy", "cupy_backends"):
                    _state = (
                        None,
                        "CuPy is not installed (the GPU code needs it: "
                        'pip install "petakit5d[gpu]" or see the GPU section of the README).',
                    )
                else:
                    _state = (None, f"CuPy failed to import: {e}")
            except Exception as e:  # (a broken install raises anything)
                _state = (None, f"CuPy failed to import: {type(e).__name__}: {e}")
            else:
                _state = (cp, "")
        return _state


def import_cupy() -> ModuleType:
    """CuPy, or GPUUnavailableError with the reason."""
    cp, reason = cupy_or_none()
    if cp is None:
        raise GPUUnavailableError(reason)
    return cp
