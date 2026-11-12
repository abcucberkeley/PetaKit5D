"""The GPU foundation (petakit5d.gpu): device choice, compiled kernels, cuFFT plans, memory.

Ports the intent of tests/tools/UsableGPUNumTest.m (counting, reasons, no GPU found) and the GPU
fallback parts of tests/decon/DeconFunctionTest.m (gpuRequestWithoutGPUFallsBackToCPU,
noUsableGPUWarningGivesMATLABsReason), with the "required" use of the parked change. The
forward compatibility tests of UsableGPUNumTest have no counterpart: GPUs newer than NVRTC
knows run PTX compiled by the driver (test_arch_cap_gives_ptx_for_the_older_arch and the B200
runs of the decision record).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
import warnings
from pathlib import Path

import pytest

import petakit5d.gpu
from petakit5d.gpu import device, fft, kernels, memory

SRC_ROOT = str(Path(petakit5d.gpu.__file__).resolve().parents[2])


def _python(code: str, **env: str | None) -> subprocess.CompletedProcess[str]:
    e = dict(os.environ)
    e["PYTHONPATH"] = SRC_ROOT + os.pathsep + e.get("PYTHONPATH", "")
    for k, v in env.items():
        if v is None:
            e.pop(k, None)
        else:
            e[k] = v
    return subprocess.run(
        [sys.executable, "-c", code], env=e, capture_output=True, text=True, timeout=300
    )


@pytest.fixture
def no_count_cache(monkeypatch):
    monkeypatch.setattr(device, "_count", None)
    yield
    device._count = None


# ---- without CuPy and without GPUs -------------------------------------------------------------


def test_import_imports_neither_cupy_nor_numpy():
    r = _python(
        "import sys\nimport petakit5d.gpu, petakit5d.gpu.device, petakit5d.gpu.kernels, "
        "petakit5d.gpu.fft, petakit5d.gpu.memory, petakit5d.gpu.selftest\n"
        'print(sorted(m for m in ("cupy", "numpy") if m in sys.modules))'
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "[]"


_BLOCK_CUPY = """
import importlib.abc, sys, warnings
class _NoCuPy(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name.split('.')[0] in ('cupy', 'cupy_backends'):
            raise ModuleNotFoundError(f"No module named {name!r}", name=name)
sys.meta_path.insert(0, _NoCuPy())
"""


def test_without_cupy_nothing_is_usable_and_says_why():
    code = (
        _BLOCK_CUPY
        + """
import json
from petakit5d.gpu import device, kernels, memory
out = {'count': list(device.usable_gpu_count())}
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    out['auto'] = device.gpu_for(None, purpose='the deconvolution')
    out['warning'] = [str(x.message) for x in w
                      if issubclass(x.category, device.GPUFallbackWarning)]
try:
    device.gpu_for(True, purpose='the deconvolution')
except device.NoUsableGPUError as e:
    out['required'] = str(e)
out['off'] = device.gpu_for(False)
for f in (device.select_device, lambda: kernels.load_source('')):
    try:
        f()
    except device.GPUUnavailableError as e:
        out.setdefault('errors', []).append(str(e))
with memory.task_scope():
    pass
print(json.dumps(out))
"""
    )
    r = _python(code, CUDA_VISIBLE_DEVICES="0")
    assert r.returncode == 0, r.stderr
    out = json.loads(r.stdout)
    assert out["count"][0] == 0 and out["count"][2] == 0
    assert "CuPy is not installed" in out["count"][1]
    assert out["auto"] is None and out["off"] is None
    assert len(out["warning"]) == 1
    assert out["warning"][0].startswith(
        "No GPU is available (0 found, CUDA_VISIBLE_DEVICES='0'), so the "
        "deconvolution runs on the CPU. CuPy is not installed"
    )
    assert out["required"].startswith(
        "No GPU is available (0 found, CUDA_VISIBLE_DEVICES='0'), but the "
        "deconvolution is to run on a GPU (use_gpu=True). CuPy is not installed"
    )
    assert all("CuPy is not installed" in e for e in out["errors"]) and len(out["errors"]) == 2


def test_no_gpu_found_with_cuda_visible_devices_empty():
    # UsableGPUNumTest/noGPUFound: none usable, none found, no reason (CuPy installed or not)
    r = _python(
        "import json\nfrom petakit5d.gpu import device, _cupy\n"
        "print(json.dumps([list(device.usable_gpu_count()), _cupy.cupy_or_none()[0] is not None]))",
        CUDA_VISIBLE_DEVICES="",
    )
    assert r.returncode == 0, r.stderr
    (n, reason, found), has_cupy = json.loads(r.stdout)
    assert (n, found) == (0, 0)
    assert reason == "" if has_cupy else "CuPy is not installed" in reason


@pytest.mark.parametrize(
    ("env", "blocked"),
    [
        ({}, False),  # not in Slurm: a workstation's GPUs
        ({"SLURM_JOB_ID": "7"}, True),  # a Slurm job without GPUs on a node that doesn't hide them
        ({"SLURM_JOB_ID": "7", "SLURM_GPUS_ON_NODE": "0"}, True),
        ({"SLURM_JOB_ID": "7", "SLURM_GPUS_ON_NODE": "1"}, False),
        ({"SLURM_JOB_ID": "7", "SLURM_JOB_GPUS": "0"}, False),  # GPU 0, not "none"
        ({"SLURM_JOB_ID": "7", "SLURM_STEP_GPUS": "3"}, False),
        ({"SLURM_JOB_ID": "7", "SLURM_JOB_GPUS": " "}, True),
        ({"SLURM_JOB_ID": "7", "CUDA_VISIBLE_DEVICES": "2"}, False),  # set by Slurm or by hand
        ({"SLURM_JOB_ID": "7", "CUDA_VISIBLE_DEVICES": ""}, False),  # (CUDA then lists none)
        ({"SLURM_JOB_ID": ""}, False),
    ],
)
def test_gpus_of_a_slurm_job(env, blocked):
    reason = device.job_gpu_reason(env, 4)
    assert bool(reason) == blocked
    if blocked:
        assert (
            "job 7 has no GPU" in reason and "4 GPUs" in reason and "CUDA_VISIBLE_DEVICES" in reason
        )
    assert device.job_gpu_reason(env, 0) == ""


def test_use_gpu_modes():
    assert [device.use_gpu_mode(v) for v in (True, None, False)] == ["required", "auto", "off"]
    with pytest.raises(TypeError):
        device.use_gpu_mode(1)


def test_no_usable_gpu_warning_and_error_give_the_reason(monkeypatch, no_count_cache):
    # DeconFunctionTest/noUsableGPUWarningGivesMATLABsReason; decon:noUsableGPU of the parked change
    monkeypatch.setattr(
        device, "_count", device.GPUCount(0, 1, "Selecting GPU 0 failed: stub reason.")
    )
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "5")
    with pytest.warns(device.GPUFallbackWarning) as w:
        assert device.gpu_for(None, purpose="the deconvolution") is None
    msg = str(w[0].message)
    assert msg.startswith("No GPU is available (1 found, CUDA_VISIBLE_DEVICES='5')")
    assert msg.endswith(
        "so the deconvolution runs on the CPU. Selecting GPU 0 failed: stub reason."
    )
    with pytest.raises(
        device.NoUsableGPUError,
        match=r"but the deconvolution is to run on a GPU \(use_gpu=True\)\. "
        r"Selecting GPU 0 failed: stub reason\.$",
    ):
        device.gpu_for(True, purpose="the deconvolution")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert device.gpu_for(False) is None
    assert device.usable_gpu_count() == (0, "Selecting GPU 0 failed: stub reason.", 1)


def test_no_gpu_found_has_no_reason_in_the_message(monkeypatch, no_count_cache):
    monkeypatch.setattr(device, "_count", device.GPUCount(0, 0, ""))
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    with pytest.warns(
        device.GPUFallbackWarning,
        match=r"^No GPU is available \(0 found, CUDA_VISIBLE_DEVICES=''\), "
        r"so the computation runs on the CPU\.$",
    ):
        device.gpu_for(None)


@pytest.mark.parametrize(
    ("cc", "archs", "target", "expected"),
    [
        (80, (70, 75, 80, 90), None, (80, "cubin")),
        (86, (70, 75, 80, 86, 90), None, (86, "cubin")),
        (100, (70, 80, 90), None, (90, "ptx")),  # B200 and NVRTC 12.2: PTX for the driver
        (103, (80, 90, 100, 120), None, (100, "ptx")),  # newer than one known, older than another
        (120, (80, 90, 100, 120), None, (120, "cubin")),
        (70, (70, 80), "ptx", (70, "ptx")),
        (100, (70, 80, 90, 100), 90, (90, "ptx")),  # as if NVRTC knew nothing above 9.0
        (80, (70, 80, 90, 100), 90, (80, "cubin")),
        (75, (70, 80), None, (70, "ptx")),
    ],
)
def test_compile_target(cc, archs, target, expected):
    assert kernels.compile_target(cc, archs, target) == expected


def test_compile_target_rejects_gpus_older_than_known():
    with pytest.raises(
        kernels.KernelCompileError,
        match=r"compute capability 6\.1 is older than the oldest NVRTC "
        r"compiles for \(7\.5\)",
    ):
        kernels.compile_target(61, (75, 80))
    with pytest.raises(ValueError):
        kernels.compile_target(80, (80,), "sass")


def test_target_from_the_environment(monkeypatch):
    monkeypatch.setenv("PETAKIT5D_GPU_TARGET", "PTX")
    assert kernels._target_from_env() == "ptx"
    monkeypatch.setenv("PETAKIT5D_GPU_TARGET", "90")
    assert kernels._target_from_env() == 90
    monkeypatch.setenv("PETAKIT5D_GPU_TARGET", "x")
    with pytest.raises(ValueError):
        kernels._target_from_env()
    monkeypatch.delenv("PETAKIT5D_GPU_TARGET")
    assert kernels._target_from_env() is None


@pytest.mark.parametrize(
    ("shape", "dims"),
    [
        ((20, 40, 30), (20, 40, 30)),  # a C array is cuFFT's row-major order
        ((1, 40, 30), (40, 30)),  # size-1 axes other than the last left out, as in the mex
        ((20, 1, 30), (20, 30)),
        ((1, 1, 30), (30,)),
        ((20, 40, 1), (20, 40, 1)),
        ((40, 30), (40, 30)),
        ((7,), (7,)),
    ],
)
def test_plan_dims_as_the_mex(shape, dims):
    assert fft.plan_dims(shape) == dims


@pytest.mark.parametrize(
    ("native", "matlab", "mex_dims"),
    [
        # MATLAB's (y, x, z) Fortran volume is the C array (z, x, y); the mex's plan is cuFFT's
        # row-major (nz, nx, ny) with size-1 dimensions other than y left out
        ((20, 40, 30), (20, 30, 40), (20, 30, 40)),
        ((1, 40, 30), (1, 30, 40), (30, 40)),
        ((20, 40, 1), (20, 1, 40), (20, 40)),
        ((20, 1, 30), (20, 30, 1), (20, 30, 1)),
        ((40, 30), (30, 40), (30, 40)),
    ],
)
def test_matlab_shape_gives_the_mex_plan(native, matlab, mex_dims):
    assert fft.matlab_shape(native) == matlab
    assert fft.plan_dims(fft.matlab_shape(native)) == mex_dims
    with pytest.raises(ValueError):
        fft.matlab_shape((7,))


@pytest.mark.parametrize("shape", [(), (4, 0, 2), (1, 2, 3, 4), (-1,)])
def test_plan_dims_rejects_bad_sizes(shape):
    with pytest.raises(ValueError):
        fft.plan_dims(shape)


def test_cache_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("PETAKIT5D_CACHE_DIR", str(tmp_path / "c"))
    assert kernels.cache_dir() == tmp_path / "c" / "cuda-kernels"
    monkeypatch.delenv("PETAKIT5D_CACHE_DIR")
    if sys.platform not in ("win32", "darwin"):
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "x"))
        assert kernels.cache_dir() == tmp_path / "x" / "petakit5d" / "cuda-kernels"


def test_cache_round_trip_and_damage(monkeypatch, tmp_path):
    monkeypatch.setenv("PETAKIT5D_CACHE_DIR", str(tmp_path))
    key = "a" * 64
    assert kernels._cache_read(key) is None
    kernels._cache_write(
        key, {"kind": "cubin", "names": {"k<float>": "_Z1kIfEvv"}}, b"\x7fELF payload"
    )
    header, payload = kernels._cache_read(key)
    assert (
        payload == b"\x7fELF payload"
        and header["names"] == {"k<float>": "_Z1kIfEvv"}
        and header["kind"] == "cubin"
    )
    f = kernels.cache_dir() / f"{key}.bin"
    assert [p.name for p in kernels.cache_dir().iterdir()] == [f.name]  # no temporary files left
    blob = f.read_bytes()
    for damaged in (blob[:-1], blob[:-1] + b"X", b"XXXXXXXX" + blob[8:], blob[:20]):
        f.write_bytes(damaged)
        assert kernels._cache_read(key) is None


def test_cache_shared_by_concurrent_writers(monkeypatch, tmp_path):
    monkeypatch.setenv("PETAKIT5D_CACHE_DIR", str(tmp_path))
    key = "b" * 64
    payload = os.urandom(1 << 20)
    errors = []

    def write():
        try:
            for _ in range(5):
                kernels._cache_write(key, {"kind": "ptx"}, payload)
                r = kernels._cache_read(key)
                assert r is None or r[1] == payload
        except Exception as e:
            errors.append(e)

    threads = [threading.Thread(target=write) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors
    assert kernels._cache_read(key)[1] == payload
    assert [p.name for p in kernels.cache_dir().iterdir()] == [f"{key}.bin"]


@pytest.mark.skipif(
    sys.platform == "win32" or (hasattr(os, "geteuid") and os.geteuid() == 0),
    reason="needs a folder the process cannot write",
)
def test_unwritable_cache_warns_once(monkeypatch, tmp_path):
    d = tmp_path / "ro"
    d.mkdir()
    d.chmod(0o500)
    monkeypatch.setenv("PETAKIT5D_CACHE_DIR", str(d))
    monkeypatch.setattr(kernels, "_warned", set())
    try:
        with pytest.warns(RuntimeWarning, match="not writable"):
            kernels._cache_write("c" * 64, {}, b"x")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            kernels._cache_write("d" * 64, {}, b"x")
    finally:
        d.chmod(0o700)


def test_header_digests_follow_header_content(tmp_path):
    (tmp_path / "k.cu").write_text('#include "h.cuh"\n')
    (tmp_path / "h.cuh").write_text("#define A 1\n")
    (tmp_path / "notes.txt").write_text("x")
    a = kernels._header_digests([str(tmp_path)])
    assert list(a) == [str(tmp_path / "h.cuh")]
    (tmp_path / "h.cuh").write_text("#define A 2\n")
    assert kernels._header_digests([str(tmp_path)]) != a
    assert kernels._digest({"o": ("-arch=sm_80",)}) != kernels._digest({"o": ("-arch=sm_90",)})


def test_cuda_source_lookup(monkeypatch, tmp_path):
    import petakit5d._ext

    built, tree = tmp_path / "build" / "cuda", tmp_path / "cpp" / "cuda"
    for d in (built, tree):
        (d / "grp").mkdir(parents=True)
    (tree / "grp" / "k.cu").write_text("")
    monkeypatch.setattr(petakit5d._ext, "cuda_source_dirs", lambda: [built, tree])
    assert kernels.cuda_source("grp", "k.cu") == tree / "grp" / "k.cu"
    (built / "grp" / "k.cu").write_text("")
    assert kernels.cuda_source("grp", "k.cu") == built / "grp" / "k.cu"  # the build's copy first
    with pytest.raises(FileNotFoundError, match=r"grp/none\.cu is not in .*build"):
        kernels.cuda_source("grp", "none.cu")
    monkeypatch.setattr(petakit5d._ext, "cuda_source_dirs", list)
    with pytest.raises(FileNotFoundError, match="none of"):
        kernels.cuda_source("grp", "k.cu")


# ---- NVRTC without a GPU (CuPy needed) --------------------------------------------------------

_FMA_SRC = r"""
template <typename T>
__global__ void axpb(const T* a, const T* x, const T* b, T* y, const long long n) {
    const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i < n) y[i] = a[i] * x[i] + b[i];
}
extern "C" __global__ void scale(float* x, const float s, const long long n) {
    const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i < n) x[i] = x[i] * s;
}
"""


def test_nvrtc_compiles_without_fma_and_lowers_names():
    pytest.importorskip("cupy")
    archs = kernels.nvrtc_archs()
    arch = max(a for a in archs if a <= 90)
    ptx, names = kernels._nvrtc_compile(
        _FMA_SRC,
        "fma.cu",
        (f"-arch=compute_{arch}", *kernels.DEFAULT_OPTIONS),
        ("axpb<float>", "axpb<double>"),
        "ptx",
    )
    text = ptx.decode("utf-8", "replace")
    assert f".target sm_{arch}" in text
    assert "fma.rn" not in text and "mul.rn.f32" in text  # --fmad=false: no contraction
    assert ".ftz" not in text  # no flush of subnormals
    assert set(names) == {"axpb<float>", "axpb<double>"} and names["axpb<float>"].startswith("_Z")
    cubin, _ = kernels._nvrtc_compile(
        _FMA_SRC, "fma.cu", (f"-arch=sm_{arch}", *kernels.DEFAULT_OPTIONS), (), "cubin"
    )
    assert cubin[:4] == b"\x7fELF"
    with pytest.raises(kernels.KernelCompileError, match=r"bad\.cu failed to compile"):
        kernels._nvrtc_compile(
            "__global__ void k() { undefined_name(); }",
            "bad.cu",
            (f"-arch=compute_{arch}",),
            (),
            "ptx",
        )


# ---- on a GPU ---------------------------------------------------------------------------------


@pytest.fixture
def gpu(monkeypatch, tmp_path):
    cp = pytest.importorskip("cupy")
    n, reason, _ = device.usable_gpu_count()
    if n == 0:
        pytest.skip(f"no usable GPU: {reason}")
    monkeypatch.setenv("PETAKIT5D_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.delenv("PETAKIT5D_GPU_TARGET", raising=False)
    monkeypatch.setattr(kernels, "_modules", {})  # (no modules of other tests)
    d = device.select_device()
    yield cp, d
    fft.clear_plan_cache()
    cp.get_default_memory_pool().free_all_blocks()


def _axpb(cp, mod, dtype, n=100_003):
    import numpy as np

    rng = np.random.default_rng(3)
    a, x, b = (cp.asarray(rng.standard_normal(n).astype(dtype)) for _ in range(3))
    y = cp.empty_like(a)
    name = f"axpb<{'float' if dtype == np.float32 else 'double'}>"
    mod[name]((n + 127) // 128, 128, (a, x, b, y, np.int64(n)))
    expected = (a.get() * x.get()).astype(dtype) + b.get()  # rounded product, then the sum (no FMA)
    return y.get(), expected


@pytest.mark.gpu
def test_kernel_runs_with_numpy_results(gpu):
    import numpy as np

    cp, d = gpu
    mod = kernels.load_source(
        _FMA_SRC, name="fma.cu", name_expressions=("axpb<float>", "axpb<double>")
    )
    arch, kind = kernels.compile_target(d.compute_capability, kernels.nvrtc_archs())
    assert (
        mod.info.target == f"{'sm' if kind == 'cubin' else 'compute'}_{arch}"
        and mod.info.jit == (kind == "ptx")
    )
    assert not mod.info.from_cache and Path(mod.info.cache_file).is_file()
    for dt in (np.float32, np.float64):
        y, expected = _axpb(cp, mod, dt)
        assert np.array_equal(y, expected)
    assert (
        kernels.load_source(
            _FMA_SRC, name="fma.cu", name_expressions=("axpb<double>", "axpb<float>")
        )
        is mod
    )


@pytest.mark.gpu
def test_ptx_compiled_by_the_driver_gives_the_same_bits(gpu):
    import numpy as np

    cp, d = gpu
    exprs = ("axpb<float>", "axpb<double>")
    native = kernels.load_source(_FMA_SRC, name="fma.cu", name_expressions=exprs)
    ptx = kernels.load_source(_FMA_SRC, name="fma.cu", name_expressions=exprs, target="ptx")
    arch = kernels.compile_target(d.compute_capability, kernels.nvrtc_archs(), "ptx")[0]
    assert ptx.info.jit and ptx.info.target == f"compute_{arch}"
    for dt in (np.float32, np.float64):
        assert np.array_equal(_axpb(cp, native, dt)[0], _axpb(cp, ptx, dt)[0])
    # the driver's cubin is cached next to the PTX (a GPU newer than NVRTC knows gets PTX anyway)
    assert (ptx is native) == native.info.jit
    assert len(list(kernels.cache_dir().glob("*.bin"))) == (2 if native.info.jit else 3)


@pytest.mark.gpu
def test_arch_cap_gives_ptx_for_the_older_arch(gpu):
    # what a GPU newer than NVRTC knows gets (a B200 with NVRTC 12.2: compute_90)
    import numpy as np

    cp, d = gpu
    older = [a for a in kernels.nvrtc_archs() if a < d.compute_capability]
    if not older:
        pytest.skip("NVRTC knows no architecture older than this GPU")
    mod = kernels.load_source(
        _FMA_SRC, name="fma.cu", name_expressions=("axpb<float>",), target=older[-1]
    )
    assert mod.info.target == f"compute_{older[-1]}" and mod.info.jit
    y, expected = _axpb(cp, mod, np.float32)
    assert np.array_equal(y, expected)


@pytest.mark.gpu
def test_a_new_process_loads_from_the_cache(gpu, tmp_path):
    code = (
        "import json\nfrom petakit5d.gpu import device, kernels\ndevice.select_device()\n"
        f'm = kernels.load_source({_FMA_SRC!r}, name="fma.cu", name_expressions=("axpb<float>",), '
        "target={target})\nprint(json.dumps([m.info.from_cache, m.info.jit]))"
    )
    for i, target in enumerate(("None", "'ptx'")):
        env = {"PETAKIT5D_CACHE_DIR": str(tmp_path / f"shared{i}"), "PETAKIT5D_GPU_TARGET": None}
        runs = [_python(code.replace("{target}", target), **env) for _ in range(2)]
        assert all(r.returncode == 0 for r in runs), runs[0].stderr + runs[1].stderr
        assert [json.loads(r.stdout)[0] for r in runs] == [False, True]


@pytest.mark.gpu
def test_our_kernels_keep_subnormals_unlike_cupy_elementwise(gpu):
    # CuPy compiles its element-wise kernels with -ftz=true: glue code that must keep subnormals
    # (as MATLAB's CPU code does) can't use them for such values
    import numpy as np

    cp, _ = gpu
    mod = kernels.load_source(_FMA_SRC, name="fma.cu")
    tiny = np.float32(1e-40)
    assert tiny != 0 and tiny < np.finfo(np.float32).tiny
    x = cp.full(4, tiny, dtype=cp.float32)
    mod["scale"](1, 4, (x, np.float32(1), np.int64(4)))
    assert (x.get() == tiny).all()
    assert (cp.full(4, tiny, dtype=cp.float32) * cp.float32(1)).get().max() == 0


SIZES = [  # (nz, ny, nx), never cubic
    (40, 48, 64),
    (17, 50, 63),
    (30, 40, 1),
    (30, 1, 40),
    (1, 30, 40),
    (30, 40),
    (33,),
    (3, 5, 7),
    (2, 96, 128),
]


@pytest.mark.gpu
@pytest.mark.parametrize("shape", SIZES)
def test_fft_round_trip_against_numpy(gpu, shape):
    import numpy as np

    cp, _ = gpu
    rng = np.random.default_rng(1)
    a = rng.random(shape, dtype=np.float32) - np.float32(0.5)
    plan = fft.RealFFTPlan(shape)
    assert plan.half_shape == (*shape[:-1], shape[-1] // 2 + 1)  # halved along x
    assert plan.dims == fft.plan_dims(shape)
    g = cp.asarray(a)
    s = plan.forward(g)
    assert s.dtype == cp.complex64 and s.shape == plan.half_shape and s.flags.c_contiguous
    assert np.array_equal(g.get(), a)  # the input stays as it is
    ref = np.fft.rfftn(a.astype(np.float64))
    assert np.abs(s.get() - ref).max() <= 1e-5 * max(np.abs(ref).max(), 1)
    back = plan.inverse(s)
    assert back.shape == shape and back.dtype == cp.float32 and back.flags.c_contiguous
    assert np.abs(back.get() / a.size - a).max() <= 1e-5
    plan.close()
    assert plan.closed
    with pytest.raises(ValueError, match="closed"):
        plan.forward(g)


@pytest.mark.gpu
@pytest.mark.parametrize(
    "shape", [(40, 48, 64), (17, 50, 63), (1, 30, 40), (20, 30, 1), (30, 40), (33,)]
)
def test_complex_fft_round_trip_against_numpy(gpu, shape):
    import numpy as np

    cp, _ = gpu
    rng = np.random.default_rng(4)
    a = (rng.random(shape, dtype=np.float32) + 1j * rng.random(shape, dtype=np.float32)).astype(
        np.complex64
    )
    plan = fft.ComplexFFTPlan(shape)
    z = cp.asarray(a)
    s = plan.forward(z)
    assert np.array_equal(z.get(), a) and s.flags.c_contiguous and s.dtype == cp.complex64
    ref = np.fft.fftn(a.astype(np.complex128))
    assert np.abs(s.get() - ref).max() <= 1e-5 * np.abs(ref).max()
    back = plan.inverse(s, out=s)  # in place
    assert back is s
    assert np.abs(back.get() / a.size - a).max() <= 1e-5
    # fftn of real data: the forward transform of complex(x), its first half the R2C half spectrum
    x = cp.asarray(np.ascontiguousarray(a.real))
    full = plan.forward(x.astype(cp.complex64)).get()
    half = fft.RealFFTPlan(shape).forward(x).get()
    assert np.abs(full[..., : shape[-1] // 2 + 1] - half).max() <= 1e-5 * np.abs(half).max()
    assert fft.get_complex_plan(shape) is fft.get_complex_plan(shape) is not fft.get_plan(shape)


MATLAB_SIZES = [
    (40, 48, 64),
    (17, 50, 63),
    (9, 101, 181),
    (1, 30, 40),
    (20, 1, 30),
    (20, 30, 1),
    (30, 40),
]


@pytest.mark.gpu
@pytest.mark.parametrize("shape", MATLAB_SIZES)
def test_matlab_plan_gives_the_bits_of_matlabs_plan_on_matlabs_layout(gpu, shape):
    # MatlabFFTPlan on (z, y, x) volumes against MATLAB's plan (ComplexFFTPlan of matlab_shape) on
    # MATLAB's layout (the C array (z, x, y)), with ifftn's division by N of each part (in double,
    # rounded to single): bit for bit (the plan whose results equal MATLAB's gpuArray fftn, ifftn)
    import numpy as np

    from petakit5d.core.arrays import matlab_order

    cp, _ = gpu
    rng = np.random.default_rng(7)
    a = rng.random(shape, dtype=np.float32)
    m = fft.MatlabFFTPlan(shape)
    assert m.spectrum_shape == fft.matlab_shape(shape)
    ref_plan = fft.ComplexFFTPlan(fft.matlab_shape(shape))
    assert m.dims == ref_plan.dims and m.work_bytes == ref_plan.work_bytes
    g = cp.asarray(a)
    s = m.forward(g)
    assert np.array_equal(g.get(), a)
    ma = cp.asarray(np.ascontiguousarray(matlab_order(a)).astype(np.complex64))
    ref = ref_plan.forward(ma)
    assert s.shape == ref.shape and s.flags.c_contiguous
    assert np.array_equal(s.get().view(np.uint32), ref.get().view(np.uint32))
    # complex input: the same bits as real input with a zero imaginary part
    assert np.array_equal(
        m.forward(g.astype(cp.complex64)).get().view(np.uint32), s.get().view(np.uint32)
    )
    # the spectrum is fftn's, in MATLAB's layout
    full = m.to_native(s).get()
    np_ref = np.fft.fftn(a.astype(np.float64))
    assert np.abs(full - np_ref).max() <= 1e-5 * np.abs(np_ref).max()
    assert np.array_equal(
        m.from_native(cp.asarray(full)).get().view(np.uint32), s.get().view(np.uint32)
    )
    # ifftn: C2C inverse, then each part divided by N
    inv_ref = ref_plan.inverse(ref.copy()).get()
    parts = (inv_ref.view(np.float32).astype(np.float64) / a.size).astype(np.float32)
    parts = parts.view(np.complex64)
    want = np.ascontiguousarray(np.swapaxes(parts, -1, -2))  # back to (z, y, x)
    got = m.inverse(s.copy())
    assert got.shape == shape and got.dtype == cp.complex64 and got.flags.c_contiguous
    assert np.array_equal(got.get().view(np.uint32), want.view(np.uint32))
    got_real = m.inverse(s.copy(), real=True).get()
    assert got_real.dtype == np.float32 and np.array_equal(
        got_real.view(np.uint32), want.real.view(np.uint32)
    )
    assert np.abs(got_real - a).max() <= 1e-5
    unscaled = m.inverse(s.copy(), scale=False).get()
    assert np.array_equal(
        unscaled.view(np.uint32), np.ascontiguousarray(np.swapaxes(inv_ref, -1, -2)).view(np.uint32)
    )
    assert fft.get_matlab_plan(shape) is fft.get_matlab_plan(shape)


@pytest.mark.gpu
def test_matlab_ifftn_divides_in_double(gpu):
    # N = 101 * 181 * 1031 is no single: MATLAB's ifftn divides in double and rounds to single
    # (dividing by single(N) changes about two thirds of the values, MATLAB R2024b on 3 GPUs)
    import numpy as np

    cp, _ = gpu
    shape = (101, 181, 1031)
    n = 101 * 181 * 1031
    assert float(np.float32(n)) != n
    m = fft.MatlabFFTPlan(shape)
    s = m.forward(cp.random.default_rng(3).random(shape, dtype=cp.float32))
    unscaled = m.inverse(s.copy(), scale=False).view(cp.float32)
    got = m.inverse(s).view(cp.float32)
    assert cp.array_equal(got, (unscaled.astype(cp.float64) / n).astype(cp.float32))
    assert not cp.array_equal(got, unscaled / cp.float32(n))
    m.close()


@pytest.mark.gpu
def test_matlab_plan_checks_its_arrays(gpu):
    cp, _ = gpu
    m = fft.get_matlab_plan((6, 8, 10))
    with pytest.raises(ValueError, match="float32 or complex64"):
        m.forward(cp.zeros((6, 8, 10), dtype=cp.float64))
    with pytest.raises(ValueError, match="C-contiguous"):
        m.forward(cp.zeros((6, 10, 8), dtype=cp.float32).transpose(0, 2, 1))
    with pytest.raises(ValueError, match="shape"):
        m.inverse(cp.zeros((6, 8, 10), dtype=cp.complex64))  # a spectrum is in MATLAB's layout
    with pytest.raises(ValueError, match="out must be"):
        m.inverse(cp.zeros((6, 10, 8), dtype=cp.complex64), out=cp.zeros((6, 8, 10), cp.float32))
    with pytest.raises(ValueError):
        fft.MatlabFFTPlan((7,))
    with fft.MatlabFFTPlan((6, 8, 10)) as closed:
        pass
    for call in (closed.forward, closed.to_native, closed.from_native):
        with pytest.raises(ValueError, match="closed"):
            call(cp.zeros((6, 8, 10), dtype=cp.complex64))


@pytest.mark.gpu
def test_fft_plans_share_one_work_area_from_the_pool(gpu):
    cp, _ = gpu
    pool = cp.get_default_memory_pool()
    shape = (100, 256, 900)  # (a work area in cuFFT 11.0 and 11.4)
    fwd = fft.RealFFTPlan(shape, inverse=False)
    used = pool.used_bytes()
    both = fft.RealFFTPlan(shape)
    assert (
        len(fwd.work_sizes) == 1
        and len(both.work_sizes) == 2
        and both.work_sizes[0] == fwd.work_bytes
    )
    assert both.work_bytes == max(both.work_sizes)  # one area for the two plans
    assert pool.used_bytes() - used >= both.work_bytes
    assert pool.used_bytes() - used < both.work_bytes + 1024
    used = pool.used_bytes()
    both.close()
    assert used - pool.used_bytes() >= both.work_bytes
    with pytest.raises(ValueError, match="without the inverse"):
        fwd.inverse(cp.zeros(fwd.half_shape, dtype=cp.complex64))
    # MATLAB's plan: one C2C work area, as F09's plan of the same (MATLAB) shape
    used = pool.used_bytes()
    m = fft.MatlabFFTPlan(shape)
    assert m.work_sizes == fft.ComplexFFTPlan(fft.matlab_shape(shape)).work_sizes
    assert pool.used_bytes() - used >= m.work_bytes


@pytest.mark.gpu
def test_fft_checks_its_arrays(gpu):
    cp, _ = gpu
    plan = fft.get_plan((4, 6, 8))
    for bad in (
        cp.zeros((4, 6, 8), dtype=cp.float64),
        cp.zeros((8, 6, 4), dtype=cp.float32).transpose(2, 1, 0),
        cp.zeros((4, 6, 9), dtype=cp.float32),
    ):
        with pytest.raises(ValueError, match="C-contiguous float32 array of shape"):
            plan.forward(bad)
    with pytest.raises(ValueError):
        plan.forward(
            cp.zeros((4, 6, 8), dtype=cp.float32), out=cp.zeros((4, 6, 8), dtype=cp.complex64)
        )


@pytest.mark.gpu
def test_fft_runs_on_the_current_stream(gpu):
    import numpy as np

    cp, _ = gpu
    a = np.random.default_rng(2).random((20, 40, 60), dtype=np.float32)
    plan = fft.get_plan(a.shape)
    m = fft.get_matlab_plan(a.shape)
    expected = plan.forward(cp.asarray(a)).get()
    mexpected = m.forward(cp.asarray(a)).get()
    s = cp.cuda.Stream(non_blocking=True)
    with s:
        out = plan.forward(cp.asarray(a))
        mout = m.forward(cp.asarray(a))
    s.synchronize()
    assert np.array_equal(out.get(), expected)  # also: the same bits from run to run
    assert np.array_equal(mout.get(), mexpected)


@pytest.mark.gpu
def test_plan_cache(gpu):
    fft.set_plan_cache_size(2)
    p1 = fft.get_plan((16, 8, 4))
    assert fft.get_plan((16, 8, 4)) is p1
    p2 = fft.get_plan((16, 8, 5))
    p3 = fft.get_plan((16, 8, 6))  # evicts p1, the oldest (it stays usable while held)
    p4 = fft.get_plan((16, 8, 4))  # a new plan; evicts p2
    assert p4 is not p1 and not p1.closed
    assert fft.get_plan((16, 8, 5)) is not p2 and not p2.closed
    assert fft.get_plan((16, 8, 6)) is not p3  # evicted by the last two
    assert len(fft._plans) == 2
    assert fft.cached_plan_bytes() == sum(
        p.work_bytes for p in fft._plans.values()
    )  # (0: small sizes, cuFFT 11.4)
    fft.clear_plan_cache()
    assert fft.cached_plan_bytes() == 0


@pytest.mark.gpu
def test_cufft_info(gpu):
    info = fft.cufft_info()
    major = int(info["version"].split(".")[0])
    assert major >= 11 and info["version_code"] >= major * 1000
    assert info["library"]


@pytest.mark.gpu
def test_selected_device_and_free_memory(gpu):
    _, d = gpu
    assert d.free_bytes > 0 and d.total_bytes >= d.free_bytes - (1 << 30)
    huge = 4 * d.total_bytes
    with pytest.raises(device.GPUUnavailableError, match=r"free of .* GiB, and .* GiB are needed"):
        device.select_device(need_bytes=huge)
    with pytest.warns(
        device.GPUFallbackWarning, match=r"so the deconvolution runs on the CPU\. .* are needed"
    ):
        assert device.gpu_for(None, purpose="the deconvolution", need_bytes=huge) is None
    with pytest.raises(
        device.NoUsableGPUError, match=r"is to run on a GPU \(use_gpu=True\)\. .* are needed"
    ):
        device.gpu_for(True, purpose="the deconvolution", need_bytes=huge)
    assert device.gpu_for(True).index == d.index
    with pytest.raises(device.GPUUnavailableError, match="not a usable GPU"):
        device.select_device(index=999)


@pytest.mark.gpu
def test_memory_between_tasks(gpu):
    cp, _ = gpu
    pool = cp.get_default_memory_pool()
    with memory.task_scope(), memory.track_peak() as peak:
        x = cp.zeros(16 << 20, dtype=cp.float32)  # 64 MiB, freed, then again
        del x
        x = cp.zeros(16 << 20, dtype=cp.float32)
        del x
        cached = pool.free_bytes()
        assert cached >= 64 << 20
    assert 64 << 20 <= peak.peak_above_start < 96 << 20
    assert pool.free_bytes() <= cached - (
        64 << 20
    )  # given back after the task (parts of chunks in use stay)
    plan = fft.get_plan((64, 32, 16))
    work = plan.work_bytes
    del plan
    assert memory.pool_bytes()["fft_plans"] == work
    with memory.task_scope(keep_plans=True):
        pass
    assert len(fft._plans) == 1 and memory.pool_bytes()["fft_plans"] == work
    with memory.task_scope(keep_plans=False):
        pass
    after = memory.pool_bytes()
    assert not fft._plans and after["fft_plans"] == 0 and after["cached"] < 1 << 20
    assert memory.available_bytes() > 0


@pytest.mark.gpu
def test_cached_pool_blocks_count_as_free(gpu):
    cp, _ = gpu
    x = cp.empty(256 << 20, dtype=cp.uint8)  # 256 MiB, then cached by the pool
    del x
    free = cp.cuda.runtime.memGetInfo()[0]
    assert cp.get_default_memory_pool().free_bytes() >= 256 << 20
    assert device.select_device(need_bytes=free + (128 << 20)).free_bytes >= free + (256 << 20)


@pytest.mark.gpu
def test_selftest(gpu, tmp_path):
    r = _python(
        "import sys\nfrom petakit5d.gpu import selftest\n"
        f"sys.exit(selftest.main(['--json', {str(tmp_path / 's.json')!r}, "
        "'--size', '10', '20', '30']))",  # (z, y, x)
        PETAKIT5D_CACHE_DIR=str(tmp_path / "c"),
    )
    assert r.returncode == 0, r.stdout + r.stderr
    out = json.loads((tmp_path / "s.json").read_text())
    assert (
        out["ok"]
        and out["kernel_ok"]
        and out["fft_ok"]
        and out["fft"]["half_shape"] == [10, 20, 16]
        and out["fft"]["matlab_plan_dims"] == [10, 30, 20]
    )


def test_selftest_without_a_gpu():
    r = _python(
        "import sys\nfrom petakit5d.gpu import selftest\nsys.exit(selftest.main([]))",
        CUDA_VISIBLE_DEVICES="",
    )
    assert r.returncode == 2, r.stdout + r.stderr
