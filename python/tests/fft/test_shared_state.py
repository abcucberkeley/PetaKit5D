"""The process's one CPU FFT state: petakit5d._ext.fft exports its backends, planner lock and plan
cache as a capsule, and every other extension module that links the FFT driver attaches to it at
import (cpp/fft/shared.h), so that FFTW is loaded once and planned under one lock (its planner is
not thread-safe)."""

import ctypes
import importlib
import threading

import pytest

from petakit5d import _ext
from petakit5d import fft as P

CAPSULE = b"petakit5d._ext.fft._shared_state"
R2C, C2CF, C2CB, C2R = 0, 1, 2, 3


def capsule_pointer(cap):
    get = ctypes.pythonapi.PyCapsule_GetPointer
    get.restype = ctypes.c_void_p
    get.argtypes = [ctypes.py_object, ctypes.c_char_p]
    return get(cap, CAPSULE)


def fft_users():
    """The extension modules other than fft that use the FFT state (they have the hooks)."""
    mods = []
    for g in _ext.groups():
        if g == "fft":
            continue
        try:
            m = importlib.import_module(f"petakit5d._ext.{g}")
        except ImportError:
            continue
        if hasattr(m, "_fft_state_id"):
            mods.append(m)
    return mods


def test_the_fft_module_exports_its_state():
    ext = P._ext()
    assert ext._fft_state_id() == capsule_pointer(ext._shared_state)


def test_every_module_that_links_the_driver_uses_the_fft_modules_state():
    ext = P._ext()
    users = fft_users()
    for m in users:
        assert m._fft_state_id() == ext._fft_state_id(), m.__name__


def run_in_threads(jobs):
    errors = []

    def run(f, args):
        try:
            errors.extend(f(*args))
        except Exception as err:  # reported below
            errors.append(repr(err))

    threads = [threading.Thread(target=run, args=job) for job in jobs]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    return errors


@pytest.mark.parametrize("backend", P.available())
def test_planning_from_two_modules_without_the_interpreter_lock(backend):
    # the fft module and another one plan new lengths at the same time, each without the
    # interpreter lock: one planner lock and one plan cache serialize FFTW's planner, and a plan
    # made in one module is the other's
    users = fft_users()
    if not users:
        pytest.skip("no other extension module links the FFT driver")
    ext, other = P._ext(), users[0]
    be = P._backend_id(backend)
    base = 4001 + 1000 * be
    a = list(range(base, base + 300, 2))
    b = list(range(base + 1, base + 301, 2))
    errors = run_in_threads(
        [
            (ext._fft_plan_lengths, (be, a)),
            (other._fft_plan_lengths, (be, b)),
            (ext._fft_plan_lengths, (be, b)),
        ]
    )
    assert not errors, errors
    for n in (a[0], b[0], b[-1]):
        for kind in (R2C, C2CF, C2CB, C2R):
            assert other._fft_cached_plan(be, 1, kind, n, 4) == ext._fft_cached_plan(
                be, 1, kind, n, 4
            )


def test_planning_from_threads_of_the_fft_module():
    be = P._backend_id(P.available()[-1])
    lengths = [list(range(6001 + k, 6201, 4)) for k in range(4)]
    errors = run_in_threads([(P._ext()._fft_plan_lengths, (be, ls)) for ls in lengths])
    assert not errors, errors
