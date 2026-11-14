"""One OpenMP runtime in a port process on Linux: the extension modules share the libgomp of the
cpp-zarr and cpp-tiff wheels (cmake/Pk5dOpenMP.cmake, docs/decisions/0008-openmp-runtime.md)."""

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from petakit5d import _ext

pytestmark = pytest.mark.skipif(
    not os.path.exists("/proc/self/maps"), reason="the mapped libraries are read from /proc (Linux)"
)

CMAKE_DIR = Path(__file__).resolve().parents[1] / "cmake"

# the OpenMP runtimes mapped after each step, in a fresh process: ours first or the I/O first,
# each calling a parallel region of 2 threads (our parallel_threads, a zarr and a tiff of 8 chunks
# or pages written and read with 2 threads)
_PROBE = r"""
import json, os, sys
order, work = sys.argv[1], sys.argv[2]
def runtimes():
    found = []
    with open("/proc/self/maps") as f:
        for line in f:
            p = line.split(None, 5)
            name = os.path.basename(p[5].strip()) if len(p) == 6 else ""
            if name.startswith(("libgomp", "libomp", "libiomp5")):
                if p[5].strip() not in found:
                    found.append(p[5].strip())
    return found
def ours():
    from petakit5d import _ext
    assert _ext.load("core").parallel_threads(2) == 2
def io():
    import numpy as np, cppzarr, cpptiff
    vol = np.arange(8 * 64 * 64, dtype=np.uint16).reshape(8, 64, 64)
    z = os.path.join(work, order + ".zarr")
    cppzarr.write_zarr(z, vol, chunks=[2, 32, 32], num_threads=2)
    assert (cppzarr.read_zarr(z, num_threads=2) == vol).all()
    t = os.path.join(work, order + ".tif")
    cpptiff.write_tiff(t, vol, num_threads=2)
    assert (cpptiff.read_tiff(t, num_threads=2) == vol).all()
steps = {}
parts = (("ours", ours), ("io", io))
for name, step in parts if order == "ours_first" else parts[::-1]:
    step()
    steps[name] = runtimes()
print(json.dumps(steps))
"""


def _build_info():
    if not _ext.available("core"):
        pytest.skip("the extension modules aren't built")
    return _ext.load("core").build_info()


def _need_wheels():
    for name in ("cppzarr", "cpptiff"):
        if importlib.util.find_spec(name) is None:
            pytest.skip(f"{name} isn't installed")


@pytest.mark.parametrize("order", ["ours_first", "io_first"])
def test_one_openmp_runtime_with_the_io_packages(order, tmp_path):
    info = _build_info()
    _need_wheels()
    if not info["openmp_version"]:
        pytest.skip("the modules are built without OpenMP")
    run = subprocess.run(
        [sys.executable, "-c", _PROBE, order, str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert run.returncode == 0, run.stderr
    steps = json.loads(run.stdout.splitlines()[-1])
    last = steps["io" if order == "ours_first" else "ours"]
    if len(last) > 1 and info["openmp_shared"] == "off":
        pytest.skip("built with PK5D_OPENMP_SHARE=OFF: the toolchain's runtime is a second one")
    assert len(last) == 1, (
        f"{len(last)} OpenMP runtimes in one process: {last}; the modules were built with "
        f"{info['openmp']!r}. Build them with cpp-zarr installed for the build's Python "
        "(PK5D_OPENMP_SHARE=AUTO), so that they share its libgomp"
    )
    if info["openmp_shared"] not in ("", "off"):
        assert os.path.basename(last[0]) == info["openmp_shared"]
        # before the I/O packages, our modules load the copy of it next to them
        if order == "ours_first":
            first = steps["ours"]
            assert [os.path.basename(p) for p in first] == [info["openmp_shared"]]
            assert os.path.dirname(first[0]) in [os.path.realpath(p) for p in _ext.__path__]


def test_shared_runtime_sits_next_to_the_modules():
    info = _build_info()
    if info["openmp_shared"] in ("", "off"):
        pytest.skip("the modules link the toolchain's OpenMP runtime")
    folders = [Path(p) for p in _ext.__path__ if (Path(p) / info["openmp_shared"]).is_file()]
    assert folders, f"{info['openmp_shared']} is not next to the modules ({_ext.__path__})"
    assert info["openmp_shared"] in info["openmp"]


def test_build_info_names_the_runtime_the_modules_need():
    info = _build_info()
    if not info["openmp_version"]:
        pytest.skip("the modules are built without OpenMP")
    dyn = subprocess.run(
        [_readelf(), "-d", "-W", _ext.load("core").__file__],
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    ).stdout
    needed = re.findall(r"\(NEEDED\).*\[(lib(?:gomp|omp|iomp5)[^\]]*)\]", dyn)
    if info["openmp_shared"] in ("", "off"):
        # the toolchain's runtime, under its own name, not a wheel's renamed copy
        assert needed and not any(re.match(r"lib\w+-[0-9a-f]{8}\.", n) for n in needed), needed
    else:
        assert needed == [info["openmp_shared"]]


def test_finder_names_the_io_packages_libgomp():
    _need_wheels()
    run = subprocess.run(
        [sys.executable, str(CMAKE_DIR / "pk5d_io_openmp.py")],
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    )
    lines = [ln.split("\t") for ln in run.stdout.splitlines()]
    assert lines, "no libgomp found for the installed cppzarr and cpptiff"
    assert lines[0][0] == "cppzarr"
    assert all(len(ln) == 2 and os.path.basename(ln[1]).startswith("libgomp") for ln in lines)
    assert all(os.path.isabs(ln[1]) and os.path.isfile(ln[1]) for ln in lines)
    # a file once: cpptiff's line only for a libgomp that cppzarr's modules didn't map
    assert len({ln[1] for ln in lines}) == len(lines)
    # the packages aren't imported (their import needs NumPy): only their compiled modules load
    probe = (
        "import runpy, sys; runpy.run_path(sys.argv[1], run_name='__main__'); "
        "print(sorted(m for m in sys.modules if m.split('.')[0] in "
        "('cppzarr', 'cpptiff', 'numpy')))"
    )
    run = subprocess.run(
        [sys.executable, "-c", probe, str(CMAKE_DIR / "pk5d_io_openmp.py")],
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    )
    assert run.stdout.splitlines()[-1] == "[]"


def test_finder_prints_nothing_without_the_packages():
    # an interpreter that can't see them: isolated, without site-packages
    run = subprocess.run(
        [sys.executable, "-I", "-S", str(CMAKE_DIR / "pk5d_io_openmp.py")],
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    )
    assert run.stdout == ""


# --- the build's check of each module (cmake/Pk5dCheckOpenMP.cmake) ---


def _cmake():
    exe = shutil.which("cmake")
    if exe is None:
        pytest.skip("no cmake")
    return exe


def _readelf():
    exe = shutil.which("readelf")
    if exe is None:
        pytest.skip("no readelf")
    return exe


def _check(module, soname):
    """The check's exit status and its output with the line breaks of CMake's messages undone."""
    run = subprocess.run(
        [
            _cmake(),
            f"-DPK5D_MODULE={module}",
            f"-DPK5D_SONAME={soname}",
            f"-DPK5D_READELF={_readelf()}",
            "-P",
            str(CMAKE_DIR / "Pk5dCheckOpenMP.cmake"),
        ],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    return run.returncode, " ".join((run.stdout + run.stderr).split())


def _core_module():
    info = _build_info()
    if info["openmp_shared"] in ("", "off"):
        pytest.skip("the modules link the toolchain's OpenMP runtime")
    return _ext.load("core").__file__, info["openmp_shared"]


def test_check_passes_a_module_of_the_build():
    module, soname = _core_module()
    rc, out = _check(module, soname)
    assert rc == 0, out


def test_check_fails_a_module_of_another_runtime():
    module, soname = _core_module()
    rc, out = _check(module, "libgomp.so.1")
    assert rc != 0
    assert f"needs the OpenMP runtimes '{soname}' instead of the shared libgomp.so.1 alone" in out


def _compile(tmp_path, name, source, *args):
    cc = shutil.which("cc") or shutil.which("gcc")
    if cc is None:
        pytest.skip("no C compiler")
    (tmp_path / f"{name}.c").write_text(source)
    out = tmp_path / f"lib{name}.so"
    subprocess.run(
        [cc, "-shared", "-fPIC", "-o", str(out), str(tmp_path / f"{name}.c"), *args],
        check=True,
        capture_output=True,
        timeout=60,
    )
    return out


# --- the build's choice of a runtime (cmake/Pk5dOpenMP.cmake), in a project of its own ---

_PROJECT = """cmake_minimum_required(VERSION 3.26)
project(pk5d_openmp_choice C CXX)
find_package(Python COMPONENTS Interpreter REQUIRED)
set(PK5D_EXT_OUTPUT_DIR "${CMAKE_BINARY_DIR}/ext")
list(APPEND CMAKE_MODULE_PATH "@DIR@")
include(Pk5dOpenMP)
message(STATUS "pk5d: found=[${PK5D_OPENMP_FOUND}] shared=[${PK5D_OPENMP_SHARED}]")
"""


def _configure(tmp_path, *defs):
    """Configures a project that includes Pk5dOpenMP.cmake: exit status, output (one line)."""
    _readelf()
    project = tmp_path / "project"
    project.mkdir(parents=True)
    (project / "CMakeLists.txt").write_text(_PROJECT.replace("@DIR@", CMAKE_DIR.as_posix()))
    gen = ["-G", "Ninja"] if shutil.which("ninja") else []
    run = subprocess.run(
        [_cmake(), "-S", str(project), "-B", str(tmp_path / "build"), *gen,
         f"-DPython_EXECUTABLE={sys.executable}", *defs],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )  # fmt: skip
    out = " ".join((run.stdout + run.stderr).split())
    if run.returncode == 0 and "found=[TRUE]" not in out:
        pytest.skip("no OpenMP for the C++ compiler")
    return run.returncode, out


def test_build_shares_the_io_packages_runtime(tmp_path):
    _need_wheels()
    first = subprocess.run(
        [sys.executable, str(CMAKE_DIR / "pk5d_io_openmp.py")],
        capture_output=True,
        text=True,
        timeout=60,
        check=True,
    ).stdout.splitlines()[0]
    soname = os.path.basename(first.split("\t")[1])
    rc, out = _configure(tmp_path)
    assert rc == 0, out
    assert f"shared=[{soname}]" in out
    copy = tmp_path / "build" / "ext" / soname
    assert copy.read_bytes() == Path(first.split("\t")[1]).read_bytes()
    rc, out = _configure(tmp_path / "off", "-DPK5D_OPENMP_SHARE=OFF")
    assert rc == 0, out
    assert "shared=[off]" in out


def test_build_refuses_a_runtime_without_the_kernels_entry_points(tmp_path):
    # a libgomp with only GOMP_parallel and the omp_ calls: the probe can't link the dynamic
    # schedules and the single construct
    stub = _compile(
        tmp_path,
        "gomp-stub",
        "void GOMP_parallel(void (*f)(void *), void *d, unsigned n, unsigned fl)"
        " { (void)n; (void)fl; f(d); }\n"
        "int omp_get_thread_num(void) { return 0; }\n"
        "int omp_get_num_threads(void) { return 1; }\n"
        "int omp_get_max_threads(void) { return 1; }\n"
        "int omp_get_num_procs(void) { return 1; }\n"
        "double omp_get_wtime(void) { return 0.0; }\n",
        "-Wl,-soname,libgomp-stub.so.1",
    )
    rc, out = _configure(tmp_path / "c", f"-DPK5D_OPENMP_SHARE={stub}")
    assert rc != 0
    assert "can't link the kernels' OpenMP constructs: undefined" in out
    assert "GOMP_loop" in out or "GOMP_single_start" in out


def test_build_leaves_a_runtime_named_as_the_toolchains(tmp_path):
    # a file named libgomp.so.1 is shared by the loader already (and AUTO leaves it alone)
    stub = _compile(tmp_path, "gomp-named", "int omp_get_num_procs(void) { return 1; }\n",
                    "-Wl,-soname,libgomp.so.1")  # fmt: skip
    rc, out = _configure(tmp_path / "c", f"-DPK5D_OPENMP_SHARE={stub}")
    assert rc != 0
    assert "libgomp is named libgomp.so.1" in out
    rc, out = _configure(tmp_path / "d", f"-DPK5D_OPENMP_SHARE={tmp_path / 'no_such.so'}")
    assert rc != 0
    assert f"no file {tmp_path / 'no_such.so'}" in out


def test_check_fails_an_entry_point_the_runtime_lacks(tmp_path):
    # linked against the shared runtime, but calling a function it doesn't define: the symbol
    # stays without a version
    module, soname = _core_module()
    lib = Path(module).parent / soname
    source = (
        "extern void GOMP_parallel(void (*)(void *), void *, unsigned, unsigned);\n"
        "extern void GOMP_no_such_entry(void);\n"
        "static void body(void *p) { (void)p; }\n"
        "void run(void) { GOMP_parallel(body, 0, 2, 0); GOMP_no_such_entry(); }\n"
    )
    rc, out = _check(_compile(tmp_path, "probe", source, str(lib)), soname)
    assert rc != 0
    assert f"the shared runtime {soname} lacks: GOMP_no_such_entry;" in out
    # without the unknown call it passes
    probe2 = _compile(tmp_path, "probe2", source.replace(" GOMP_no_such_entry();", ""), str(lib))
    rc, out = _check(probe2, soname)
    assert rc == 0, out
