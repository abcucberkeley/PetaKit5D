"""The package skeleton: a light import, the build of the extension modules and its flags."""

import json
import os
import re
import subprocess
import sys

import pytest

import petakit5d
from petakit5d import _ext


def _run_python(code: str, tmp_path, env=None, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, *args, "-c", code],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )


def test_version():
    from petakit5d._version import __version__

    assert petakit5d.__version__ == __version__
    assert re.fullmatch(r"\d+\.\d+\.\d+((a|b|rc)\d+)?(\.post\d+)?(\.dev\d+)?", __version__)


def test_import_is_light(tmp_path):
    code = (
        "import json, sys, petakit5d; "
        "print(json.dumps(sorted(m for m in sys.modules "
        "if m.split('.')[0] in ('numpy', 'scipy', 'skimage', 'cupy', 'cppzarr', 'cpptiff') "
        "or m.startswith('petakit5d.'))))"
    )
    result = _run_python(code, tmp_path, None, "-X", "importtime")
    assert json.loads(result.stdout) == ["petakit5d._version"]
    # -X importtime: "import time: self [us] | cumulative | imported package"
    cumulative = [
        int(line.split("|")[1])
        for line in result.stderr.splitlines()
        if line.startswith("import time:") and line.split("|")[-1].strip() == "petakit5d"
    ]
    assert len(cumulative) == 1
    assert cumulative[0] < 50_000, f"import petakit5d took {cumulative[0]} us"


def test_lazy_subpackage_access():
    assert petakit5d.kernels.__name__ == "petakit5d.kernels"
    with pytest.raises(AttributeError):
        petakit5d.no_such_area  # noqa: B018


def test_core_extension_builds():
    assert _ext.available("core")
    assert "core" in _ext.groups()
    core = _ext.load("core")
    info = core.build_info()
    assert info["compiler"] != "unknown"
    # OpenMP 3.1: collapse, min/max reductions; MSVC's /openmp:llvm has them but keeps _OPENMP
    # at 2.0's 200203
    msvc_llvm = info["openmp"] == "MSVC /openmp:llvm"
    assert info["openmp_version"] >= (200203 if msvc_llvm else 201107), info


def test_modules_load_after_the_other_libraries(tmp_path):
    # NumPy, SciPy and the I/O packages load their own C++ and OpenMP runtimes first, as in a
    # real task; every module must still load (Linux: libstdc++ linked statically)
    code = (
        "import importlib\n"
        "for name in ('numpy', 'scipy.ndimage', 'scipy.fft', 'skimage.transform', 'cppzarr',"
        " 'cpptiff'):\n"
        "    try:\n"
        "        importlib.import_module(name)\n"
        "    except ImportError:\n"
        "        pass\n"
        "from petakit5d import _ext\n"
        "for group in _ext.groups():\n"
        "    _ext.load(group)\n"
        "print(_ext.load('core').parallel_threads(2))\n"
    )
    result = _run_python(code, tmp_path)
    assert result.stdout.split()[-1] in ("1", "2")


def test_missing_group_has_hint():
    assert not _ext.available("no_such_group")
    with pytest.raises(ModuleNotFoundError, match="not built"):
        _ext.load("no_such_group")


def test_kernel_flags():
    info = _ext.load("core").build_info()
    flags = info["flags"].split()
    if info["compiler"].startswith("msvc"):
        assert "/fp:precise" in flags and "/fp:fast" not in flags
    else:
        assert "-ffp-contract=off" in flags
        assert "-O3" in flags
    assert not any(f.startswith("-march=native") or f in ("-ffast-math", "-Ofast") for f in flags)
    assert info["avx2"] is False  # only groups declared AVX2 get it


def test_no_fused_multiply_add():
    assert _ext.load("core").fma_contracted() is False


def test_subnormals_kept():
    core = _ext.load("core")
    assert core.subnormals_flushed() is False
    import numpy  # noqa: F401  libraries loaded later must not switch on flush-to-zero

    assert core.subnormals_flushed() is False


@pytest.mark.parametrize("num_threads", [1, 2, 3])
def test_parallel_region_uses_the_requested_threads(num_threads):
    core = _ext.load("core")
    expected = num_threads if core.build_info()["openmp_version"] else 1
    assert core.parallel_threads(num_threads) == expected


def test_parallel_threads_rejects_zero():
    with pytest.raises(ValueError):
        _ext.load("core").parallel_threads(0)


def test_thread_count_follows_the_cpus_of_the_process(tmp_path):
    # OpenMP's default team equals the CPUs the process may use (Slurm's -c, taskset)
    env = {k: v for k, v in os.environ.items() if k not in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT")}
    code = (
        "import json; from petakit5d._ext import load; c = load('core'); "
        "print(json.dumps([c.affinity_cpu_count(), c.omp_num_procs(), c.omp_max_threads()]))"
    )
    cpus, procs, max_threads = json.loads(_run_python(code, tmp_path, env).stdout)
    assert cpus >= 1
    if hasattr(os, "sched_getaffinity"):
        assert cpus == len(os.sched_getaffinity(0))
    if sys.platform != "darwin":  # no affinity masks on macOS
        assert procs == cpus
        assert max_threads == cpus


def test_cuda_source_dirs_are_folders():
    assert all(d.is_dir() for d in _ext.cuda_source_dirs())
