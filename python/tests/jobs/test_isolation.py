"""Tasks that share a process (jobs.isolation): what a task changes in its process is given back
after it, and the harness shows that tasks run one after another in one process give what fresh
processes give, bit for bit (the idea of MATLAB's tests/isolation)."""

import ctypes
import gc
import io
import json
import logging
import os
import platform
import signal
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d.jobs import isolation
from petakit5d.jobs.isolation import ProcessState, run_harness, run_in_process
from petakit5d.jobs.tasks import Task, task_record

LINUX = sys.platform.startswith("linux")


@pytest.fixture
def state():
    """the state of this process, given back after the test"""
    st = ProcessState.capture()
    yield st
    st.restore()


def test_folder_environment_and_path_are_given_back(state, tmp_path):
    os.environ["PK_ISO_NEW"] = "1"
    os.environ["PATH"] = os.environ["PATH"] + os.pathsep + "/nowhere"
    removed = next(iter(state.environ))
    del os.environ[removed]
    os.chdir(tmp_path)
    sys.path.insert(0, "/nowhere")
    changed = state.restore()
    assert {"folder", "environment", "path"} <= set(changed)
    assert os.getcwd() == state.cwd and dict(os.environ) == state.environ
    assert sys.path == state.path


def test_warning_filters_are_given_back(state):
    warnings.simplefilter("error")
    assert "warnings" in state.restore()
    assert warnings.filters == state.filters


def _warn():
    warnings.warn("pk iso once", UserWarning, stacklevel=1)


def test_warnings_shown_once_show_again_after_the_restore():
    # (as in a new process)
    shown = []
    with warnings.catch_warnings():
        warnings.simplefilter("default")
        warnings.showwarning = lambda message, *a, **k: shown.append(str(message))
        st = ProcessState.capture()
        _warn()
        _warn()
        assert st.restore() == []
        _warn()
    assert shown.count("pk iso once") == 2


@pytest.mark.skipif(
    not hasattr(signal, "SIGUSR1") or not hasattr(signal, "setitimer"), reason="POSIX signals"
)
def test_umask_signals_alarm_and_stdio_are_given_back():
    # (pytest-timeout's timer is paused here, and set again after the test)
    outer = signal.setitimer(signal.ITIMER_REAL, 0)
    try:
        state = ProcessState.capture()
        os.umask(0o077)
        old = signal.getsignal(signal.SIGUSR1)
        signal.signal(signal.SIGUSR1, lambda *a: None)
        signal.setitimer(signal.ITIMER_REAL, 100)
        sys.stdout = io.StringIO()
        changed = state.restore()
        assert {"umask", "alarm", "stdio"} <= set(changed)
        assert any("signal" in c for c in changed)
        assert signal.getsignal(signal.SIGUSR1) is old
        assert signal.getitimer(signal.ITIMER_REAL)[0] == 0
        assert sys.stdout is state.stdio[1]
        u = os.umask(0o022)
        os.umask(u)
        assert u == state.umask
        # a timer pending before the task is the process's own: it stays
        signal.setitimer(signal.ITIMER_REAL, 100)
        state = ProcessState.capture()
        assert "alarm" not in state.restore()
        assert signal.getitimer(signal.ITIMER_REAL)[0] > 0
    finally:
        signal.setitimer(signal.ITIMER_REAL, *outer)


def test_interpreter_settings_are_given_back(state):
    sys.setrecursionlimit(state.recursion + 100)
    sys.setswitchinterval(state.switch * 2)
    gc.set_threshold(1, 2, 3)
    sys.excepthook = lambda *a: None
    tempfile.tempdir = "/nowhere"
    logging.getLogger().setLevel(logging.DEBUG)
    logging.getLogger().addHandler(logging.NullHandler())
    changed = set(state.restore())
    assert {"recursion limit", "switch interval", "gc", "hooks", "tempfile", "logging"} <= changed
    assert sys.getrecursionlimit() == state.recursion and gc.get_threshold() == state.gc_state[1]
    assert tempfile.tempdir == state.tempdir and sys.excepthook is state.hooks[0]


@pytest.mark.skipif(not hasattr(os, "sched_setaffinity"), reason="no CPU affinity here")
def test_cpus_are_given_back(state):
    cpus = os.sched_getaffinity(0)
    if len(cpus) < 2:
        pytest.skip("one CPU")
    os.sched_setaffinity(0, {min(cpus)})
    assert "threads" in state.restore()
    assert os.sched_getaffinity(0) == cpus


def _rounding():
    # libm and its FE_TOWARDZERO
    lib = ctypes.CDLL("libm.so.6")
    return lib, {"x86_64": 0xC00, "aarch64": 0xC00000}[platform.machine()]


@pytest.mark.skipif(
    not LINUX or platform.machine() not in ("x86_64", "aarch64"), reason="glibc fenv"
)
def test_rounding_mode_is_given_back(state):
    lib, up = _rounding()
    one_third = (np.float32(1) / np.float32(3)).tobytes()
    lib.fesetround(up)
    assert (np.float32(1) / np.float32(3)).tobytes() != one_third
    assert "fp" in state.restore()
    assert lib.fegetround() == 0
    assert (np.float32(1) / np.float32(3)).tobytes() == one_third


def test_numpy_and_decimal_state_stay_in_the_tasks_context(state, tmp_path, monkeypatch):
    import decimal

    (tmp_path / "ctxtask.py").write_text(
        "import contextvars, decimal, numpy as np\n"
        "SETTING = contextvars.ContextVar('pk_setting', default='as it starts')\n"
        "def run():\n"
        "    np.set_printoptions(precision=2)\n"
        "    np.seterr(all='raise')\n"
        "    decimal.getcontext().prec = 5\n"
        "    SETTING.set('the task')\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    import ctxtask

    rec = task_record(Task("ctxtask:run", path=[str(tmp_path)]), "python")
    before = (np.get_printoptions(), np.geterr(), decimal.getcontext().prec)
    r = run_in_process(rec, state)
    assert r.error is None
    assert (np.get_printoptions(), np.geterr(), decimal.getcontext().prec) == before
    # (any module's context variables)
    assert ctxtask.SETTING.get() == "as it starts"


def test_random_numbers_are_seeded_anew(state):
    import random

    # (NumPy's legacy global random numbers: what a task could seed)
    random.seed(1)
    np.random.seed(1)  # noqa: NPY002
    state.restore()
    a = (random.random(), np.random.random())  # noqa: NPY002
    random.seed(1)
    np.random.seed(1)  # noqa: NPY002
    state.restore()
    assert (random.random(), np.random.random()) != a  # noqa: NPY002


TASKS = """
import os


def ok(out):
    with open(out, "w") as fh:
        fh.write("ok")


def fail():
    raise RuntimeError("on purpose")


def leave(code):
    raise SystemExit(code)


def interrupt():
    raise KeyboardInterrupt
"""


@pytest.fixture
def tasks(tmp_path, monkeypatch):
    (tmp_path / "isotiny.py").write_text(TASKS)
    monkeypatch.syspath_prepend(str(tmp_path))

    def rec(name, *args):
        return task_record(Task(f"isotiny:{name}", list(args), path=[str(tmp_path)]), "python")

    return rec


def test_run_in_process_returns_a_tasks_error(state, tasks, tmp_path, capfd):
    calls = []
    r = run_in_process(
        tasks("fail"), state, lambda: calls.append("before"), lambda t, e: calls.append(e)
    )
    assert isinstance(r.error, RuntimeError) and "on purpose" in r.text
    assert calls[0] == "before" and calls[1] is r.error
    assert "RuntimeError: on purpose" in capfd.readouterr().err
    assert run_in_process(tasks("ok", str(tmp_path / "o")), state).error is None
    assert run_in_process(tasks("leave", 0), state).error is None
    assert isinstance(run_in_process(tasks("leave", 3), state).error, SystemExit)


def test_interrupts_go_through(state, tasks):
    seen = []
    with pytest.raises(KeyboardInterrupt):
        run_in_process(tasks("interrupt"), state, None, lambda t, e: seen.append(type(e)))
    assert seen == [KeyboardInterrupt]


# --- the harness -------------------------------------------------------------------------------

HARNESS_TASKS = '''
"""task functions of the port, file in, file out, for the isolation harness"""

import json
import os

import numpy as np


def prepare(work):
    from petakit5d.io.tiff import write_tiff

    rng = np.random.default_rng(7)
    os.makedirs(f"{work}/in", exist_ok=True)
    vol = np.asfortranarray(rng.integers(0, 4000, size=(48, 40, 12), dtype=np.uint16))
    write_tiff(vol, f"{work}/in/vol.tif")
    np.save(f"{work}/in/vol.npy", np.asfortranarray(rng.standard_normal((32, 24, 10))))
    np.save(f"{work}/in/vol32.npy", np.asfortranarray(rng.standard_normal((20, 18, 6)), np.float32))


def tiff_to_zarr(work, name, block):
    from petakit5d.io.tiff import read_tiff
    from petakit5d.io.zarr import write_zarr

    write_zarr(read_tiff(f"{work}/in/vol.tif"), f"{work}/out/{name}.zarr", block_size=block)


def tiff_mip(work, name, axis):
    from petakit5d.io.tiff import read_tiff, write_tiff

    vol = read_tiff(f"{work}/in/vol.tif")
    write_tiff(np.asfortranarray(vol.max(axis=axis)), f"{work}/out/{name}.tif")


def fftn_of(work, name, src):
    from petakit5d.matlabcompat.fft import fftn

    np.save(f"{work}/out/{name}.npy", fftn(np.load(f"{work}/in/{src}")))


def ifftn_of(work, name, src):
    from petakit5d.matlabcompat.fft import fftn, ifftn

    np.save(f"{work}/out/{name}.npy", ifftn(fftn(np.load(f"{work}/in/{src}"))))


def cast_tiff(work, name, cls):
    from petakit5d.core.dtypes import cast
    from petakit5d.io.tiff import write_tiff

    x = np.load(f"{work}/in/vol32.npy") * 1000
    write_tiff(np.asfortranarray(cast(x, cls)), f"{work}/out/{name}.tif")


def info(work, name):
    from petakit5d.io.imagesize import get_image_bounding_box, get_image_size
    from petakit5d.io.tiff import read_tiff

    vol = read_tiff(f"{work}/in/vol.tif")
    vol[vol < 3000] = 0
    with open(f"{work}/out/{name}.json", "w") as fh:
        json.dump({"size": list(get_image_size(f"{work}/in/vol.tif")),
                   "bbox": list(get_image_bounding_box(vol))}, fh)


def published(work, name):
    from petakit5d.jobs.publish import publish_result, temp_path

    final = f"{work}/out/{name}"
    tmp = temp_path(final)
    os.makedirs(tmp)
    with open(f"{tmp}/part.txt", "w") as fh:
        fh.write(name)
    publish_result(tmp, final)


def leak(work, name):
    """what a careless task leaves: folder, environment, umask, print options, warnings, the
    rounding mode"""
    import ctypes
    import platform
    import warnings

    os.chdir("/")
    os.environ["PK_ISO_LEAK"] = "leaked"
    os.umask(0o077)
    np.set_printoptions(precision=3)
    warnings.simplefilter("ignore")
    import sys

    if sys.platform.startswith("linux") and platform.machine() in ("x86_64", "aarch64"):
        zero = {"x86_64": 0xC00, "aarch64": 0xC00000}[platform.machine()]
        ctypes.CDLL("libm.so.6").fesetround(zero)
    with open(f"{work}/out/{name}.txt", "w") as fh:
        fh.write("left things behind")


def victim(work, name):
    """what the next task sees of them"""
    import stat
    import warnings

    made = f"{work}/out/{name}_mode.txt"
    with open(made, "w") as fh:
        fh.write("here")
    seen_warnings = []
    old = warnings.showwarning
    warnings.showwarning = lambda *a, **k: seen_warnings.append(a[0])
    try:
        warnings.warn("pk iso", UserWarning, stacklevel=1)
    finally:
        warnings.showwarning = old
    seen = {
        "in the work folder": os.getcwd() == work,
        "env": os.environ.get("PK_ISO_LEAK"),
        "mode": oct(stat.S_IMODE(os.stat(made).st_mode)),
        "print": str(np.array([1 / 3, 2 / 3])),
        "third": (np.float32(1) / np.float32(3)).tobytes().hex(),
        "warned": len(seen_warnings),
    }
    with open(f"{work}/out/{name}.json", "w") as fh:
        json.dump(seen, fh)
'''


def harness_spec(tmp_path, with_leaks=True):
    (tmp_path / "isotasks.py").write_text(HARNESS_TASKS)
    work = str(tmp_path / "work")

    def rec(func, *args):
        t = Task(f"isotasks:{func}", [work, *args], path=[str(tmp_path)])
        return task_record(t, "python")

    tasks = [
        rec("tiff_to_zarr", "zarr_a", [16, 16, 4]),
        rec("tiff_to_zarr", "zarr_b", [48, 40, 12]),
        rec("tiff_mip", "mip_z", 2),
        rec("tiff_mip", "mip_y", 0),
        rec("fftn_of", "fft64", "vol.npy"),
        rec("fftn_of", "fft32", "vol32.npy"),
        rec("ifftn_of", "ifft64", "vol.npy"),
        rec("cast_tiff", "cast16", "uint16"),
        rec("cast_tiff", "cast8", "uint8"),
        rec("info", "info"),
        rec("published", "result"),
        rec("victim", "victim_1"),
    ]
    if with_leaks:
        tasks += [rec("leak", "leak_1"), rec("victim", "victim_2"), rec("leak", "leak_2")]
    return {"prepare": rec("prepare"), "tasks": tasks}, work


# (the harness runs where the runners run: Linux; POSIX file modes and signals)
posix = pytest.mark.skipif(os.name == "nt", reason="POSIX")


@posix
def test_harness_five_random_orders_equal_fresh_processes(tmp_path):
    spec, work = harness_spec(tmp_path)
    res = run_harness(spec, work, str(tmp_path / "results"), orders=5, seed=11, jobs=4)
    assert len(res.runs) == 5
    assert res.ok, res.text()
    assert all(r["files"] >= 15 for r in res.runs)
    # (the leaks were seen and given back)
    changed = {c for r in res.runs for t in r["tasks"] for c in t["changed"]}
    assert {"folder", "environment", "umask", "warnings"} <= changed
    if LINUX and platform.machine() in ("x86_64", "aarch64"):
        assert "fp" in changed
    summary = json.loads(Path(res.results, "summary.json").read_text())
    assert [r["name"] for r in summary["runs"]] == [f"order_{k}" for k in range(1, 6)]


@posix
def test_harness_sees_what_tasks_leave_without_the_restore(tmp_path):
    spec, work = harness_spec(tmp_path)
    n = len(spec["tasks"])
    # the leaking task first, then everything else
    order = [n - 3] + [i for i in range(n) if i != n - 3]
    res = run_harness(spec, work, str(tmp_path / "results"), orders=[order], reset=False)
    assert not res.ok
    text = res.text()
    assert "victim_1" in text and "DIFFERENT" in text


@posix
def test_harness_command_line(tmp_path):
    spec, work = harness_spec(tmp_path, with_leaks=False)
    spec_file = tmp_path / "spec.json"
    spec_file.write_text(json.dumps(spec))
    assert isolation.main(["run", str(spec_file), work, str(tmp_path / "res"), "1", "3"]) == 0
    assert isolation.main([]) == 2
