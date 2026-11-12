"""CUDA kernels compiled at run time by NVRTC, cached on disk.

The .cu sources ship with the package (cpp/cuda/<group>/ in the source tree) and are compiled
for the GPU in use the first time a process on a GPU of that kind loads them, so one install
serves every NVIDIA GPU and needs no CUDA toolkit to build.

Target: the GPU's own architecture when NVRTC knows it (a cubin). A GPU newer than NVRTC
knows (B200, compute capability 10.0, with NVRTC 12.2) gets PTX for the newest architecture
NVRTC knows below it, which the driver compiles for the GPU; that result is cached as a cubin
too, so the driver's own JIT cache isn't needed. Cubins are preferred: the driver can compile
PTX only from an NVRTC no newer than itself, while it runs cubins of newer ones.

Compile options: --std=c++17 --fmad=false. No fused multiply-adds (as -ffp-contract=off in
the C++ kernels) and no flush of subnormals, unlike CuPy's RawModule, which adds -ftz=true and
leaves NVRTC's --fmad=true.

Cache: one file per compiled module, named by a SHA-256 of everything that decides it (source,
headers next to it, options, NVRTC version; for a driver-compiled cubin also the driver version
and the GPU's compute capability), written to a temporary name and renamed, so parallel jobs
can share the folder: PETAKIT5D_CACHE_DIR/cuda-kernels, else the user's cache folder
(~/.cache/petakit5d on Linux, ~/Library/Caches/petakit5d on macOS, %LOCALAPPDATA%/petakit5d
on Windows). A damaged file is compiled again.

PETAKIT5D_GPU_TARGET (testing): 'ptx' compiles PTX even when NVRTC knows the GPU, and a number
(e.g. 90) acts as if NVRTC knew no architecture above it.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import struct
import sys
import threading
import time
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any, NamedTuple

from petakit5d.gpu._cupy import import_cupy

DEFAULT_OPTIONS = ("--std=c++17", "--fmad=false")
_FORMAT = 1
_MAGIC = b"PK5DCU01"
_HEADER_SUFFIXES = (".cuh", ".h", ".hpp")
_CUDA_ERROR_UNSUPPORTED_PTX_VERSION = 222

_lock = threading.RLock()
_modules: dict[tuple[int, str], Module] = {}
_warned: set[str] = set()


class KernelCompileError(RuntimeError):
    """A CUDA source failed to compile or load, with the compiler's log."""


class BuildInfo(NamedTuple):
    """How a module was made for this GPU."""

    target: str  # 'sm_80': a cubin from NVRTC; 'compute_90': PTX, compiled by the driver
    jit: bool  # the driver compiled PTX for this GPU
    from_cache: bool  # loaded from the disk cache (no compile in this process)
    seconds: float  # time of load_module (compile, driver compile, cache read, module load)
    cache_file: str | None


class Kernel:
    """One kernel of a module: kernel(grid, block, args, shared_mem=0, stream=None).

    args are CuPy arrays and NumPy scalars of the kernel's parameter types (np.int64(n),
    np.float32(x)); grid and block are ints or tuples of up to 3 ints.
    """

    __slots__ = ("_function", "name")

    def __init__(self, function: Any, name: str) -> None:
        self._function = function
        self.name = name

    def __call__(
        self,
        grid: int | Sequence[int],
        block: int | Sequence[int],
        args: Sequence[Any],
        *,
        shared_mem: int = 0,
        stream: Any = None,
    ) -> None:
        self._function(_dim3(grid), _dim3(block), tuple(args), shared_mem, stream)


class Module:
    """A compiled CUDA module on one GPU; its kernels by name or by name expression."""

    def __init__(self, module: Any, names: dict[str, str], info: BuildInfo, device: int) -> None:
        self._module = module
        self._names = names
        self._kernels: dict[str, Kernel] = {}
        self.info = info
        self.device = device

    def get_function(self, name: str) -> Kernel:
        k = self._kernels.get(name)
        if k is None:
            k = self._kernels[name] = Kernel(
                self._module.get_function(self._names.get(name, name)), name
            )
        return k

    __getitem__ = get_function


def _dim3(v: int | Sequence[int]) -> tuple[int, int, int]:
    t = (int(v),) if isinstance(v, int) else tuple(int(x) for x in v)
    if not 1 <= len(t) <= 3:
        raise ValueError(f"grid and block take 1 to 3 dimensions, not {len(t)}")
    return t + (1,) * (3 - len(t))


def nvrtc_version() -> tuple[int, int]:
    """NVRTC's version, e.g. (12, 9)."""
    import_cupy()
    from cupy_backends.cuda.libs import nvrtc

    return tuple(nvrtc.getVersion())


def nvrtc_archs() -> tuple[int, ...]:
    """The architectures NVRTC compiles for (compute capability as major * 10 + minor)."""
    import_cupy()
    from cupy_backends.cuda.libs import nvrtc

    return tuple(sorted(int(a) for a in nvrtc.getSupportedArchs()))


def compile_target(
    cc: int, archs: Sequence[int], target: str | int | None = None
) -> tuple[int, str]:
    """(architecture, 'cubin' or 'ptx') to compile for a GPU of compute capability cc.

    cc and archs as major * 10 + minor. The GPU's own architecture as a cubin if NVRTC knows it,
    else PTX of the newest known architecture below it, for the driver to compile. target:
    'ptx' for PTX always, a number to ignore known architectures above it (see the module).
    """
    known = sorted(int(a) for a in archs)
    if isinstance(target, int) and not isinstance(target, bool):
        known = [a for a in known if a <= target]
    elif target not in (None, "ptx"):
        raise ValueError(f"target must be None, 'ptx' or an architecture number, not {target!r}")
    if cc in known and target != "ptx":
        return cc, "cubin"
    older = [a for a in known if a <= cc]
    if not older:
        oldest = f"{known[0] // 10}.{known[0] % 10}" if known else "none"
        raise KernelCompileError(
            f"A GPU of compute capability {cc // 10}.{cc % 10} is older than the oldest "
            f"NVRTC compiles for ({oldest})."
        )
    return older[-1], "ptx"


def _target_from_env() -> str | int | None:
    v = os.environ.get("PETAKIT5D_GPU_TARGET", "").strip().lower()
    if not v:
        return None
    if v == "ptx":
        return "ptx"
    try:
        return int(v)
    except ValueError:
        raise ValueError(
            f"PETAKIT5D_GPU_TARGET must be 'ptx' or an architecture number, not {v!r}"
        ) from None


def cache_dir() -> Path:
    """The folder of the compiled kernels (see the module)."""
    d = os.environ.get("PETAKIT5D_CACHE_DIR")
    if d:
        base = Path(d)
    elif sys.platform == "win32":
        base = (
            Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local") / "petakit5d"
        )
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches" / "petakit5d"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "petakit5d"
    return base / "cuda-kernels"


def cuda_source(group: str, filename: str) -> Path:
    """The path of a shipped CUDA source: cpp/cuda/<group>/<filename>.

    Looked up where the build puts the sources (petakit5d/_ext/cuda, a build tree's first), then
    in the source tree (petakit5d._ext.cuda_source_dirs).
    """
    from petakit5d._ext import cuda_source_dirs

    dirs = cuda_source_dirs()
    for root in dirs:
        p = root / group / filename
        if p.is_file():
            return p
    where = (
        ", ".join(str(d) for d in dirs)
        or "any folder (none of petakit5d._ext.cuda_source_dirs() exists)"
    )
    raise FileNotFoundError(f"The CUDA source {group}/{filename} is not in {where}.")


def load_module(
    path: str | os.PathLike[str],
    *,
    name_expressions: Sequence[str] = (),
    options: Sequence[str] = (),
    include_dirs: Sequence[str | os.PathLike[str]] = (),
    target: str | int | None = None,
) -> Module:
    """Compile (or load from the cache) a .cu file for the current GPU.

    Headers are found next to the file and in include_dirs. name_expressions: template
    instantiations to compile ('kernel<float>'), then found by that name. options: NVRTC
    options added to DEFAULT_OPTIONS. target: see compile_target (default: PETAKIT5D_GPU_TARGET).
    """
    p = Path(path)
    return load_source(
        p.read_text(encoding="utf-8"),
        name=p.name,
        name_expressions=name_expressions,
        options=options,
        include_dirs=(p.parent, *include_dirs),
        target=target,
    )


def load_source(
    code: str,
    *,
    name: str = "kernels.cu",
    name_expressions: Sequence[str] = (),
    options: Sequence[str] = (),
    include_dirs: Sequence[str | os.PathLike[str]] = (),
    target: str | int | None = None,
) -> Module:
    """As load_module, for source code in a string (name: its file name in compiler messages)."""
    cp = import_cupy()
    dev = cp.cuda.Device()
    cc = int(dev.compute_capability)
    nv = nvrtc_version()
    arch, kind = compile_target(
        cc, nvrtc_archs(), target if target is not None else _target_from_env()
    )
    tgt = f"{'sm' if kind == 'cubin' else 'compute'}_{arch}"
    dirs = [str(Path(d).resolve()) for d in include_dirs]
    opts = (f"-arch={tgt}", *DEFAULT_OPTIONS, *(f"-I{d}" for d in dirs), *options)
    exprs = tuple(sorted(set(name_expressions)))
    key = _digest(
        {
            "format": _FORMAT,
            "code": code,
            "headers": _header_digests(dirs),
            "options": opts,
            "nvrtc": nv,
            "names": exprs,
        }
    )
    driver = cp.cuda.runtime.driverGetVersion()
    load_key = (
        key
        if kind == "cubin"
        else _digest({"format": _FORMAT, "ptx": key, "driver": driver, "cc": cc})
    )
    with _lock:
        m = _modules.get((dev.id, load_key))
        if m is not None:
            return m
        t0 = time.perf_counter()
        entry = _cache_read(load_key)
        from_cache = entry is not None
        if entry is None:
            src = _cache_read(key) if kind == "ptx" else None
            if src is None:
                payload, names = _nvrtc_compile(code, name, opts, exprs, kind)
                src = (
                    {"kind": kind, "target": tgt, "jit": False, "names": names, "nvrtc": list(nv)},
                    payload,
                )
                _cache_write(key, *src)
            if kind == "ptx":
                cubin = _driver_compile(src[1], name, driver)
                entry = (dict(src[0], kind="cubin", jit=True, driver=driver, cc=cc), cubin)
                _cache_write(load_key, *entry)
            else:
                entry = src
        header, cubin = entry
        module = cp.cuda.Module()
        try:
            module.load(cubin)
        except Exception as e:  # (CUDA driver errors) re-raised
            raise KernelCompileError(
                f"{name}: loading the module ({tgt}) on GPU {dev.id} failed: {e}"
            ) from e
        f = cache_dir() / f"{load_key}.bin"
        info = BuildInfo(
            target=tgt,
            jit=kind == "ptx",
            from_cache=from_cache,
            seconds=time.perf_counter() - t0,
            cache_file=str(f) if f.is_file() else None,
        )
        m = _modules[(dev.id, load_key)] = Module(
            module, dict(header.get("names", {})), info, dev.id
        )
        return m


def _nvrtc_compile(
    code: str, name: str, opts: Sequence[str], exprs: Sequence[str], kind: str
) -> tuple[bytes, dict[str, str]]:
    from cupy_backends.cuda.libs import nvrtc

    prog = nvrtc.createProgram(code, name, [], [])
    try:
        for e in exprs:
            nvrtc.addNameExpression(prog, e)
        try:
            nvrtc.compileProgram(prog, list(opts))
        except nvrtc.NVRTCError as e:
            raise KernelCompileError(
                f"{name} failed to compile ({' '.join(opts)}): {e}\n{nvrtc.getProgramLog(prog)}"
            ) from None
        payload = nvrtc.getCUBIN(prog) if kind == "cubin" else nvrtc.getPTX(prog)
        names = {e: nvrtc.getLoweredName(prog, e) for e in exprs}
    finally:
        nvrtc.destroyProgram(prog)
    return bytes(payload), names


def _driver_compile(ptx: bytes, name: str, driver: int) -> bytes:
    # the driver compiles the PTX for the current GPU (cuLink), giving a cubin to cache
    import cupy as cp

    try:
        ls = cp.cuda.function.LinkState()
        ls.add_ptr_data(ptx, name.rsplit(".", 1)[0] + ".ptx")
        return bytes(ls.complete())
    except cp.cuda.driver.CUDADriverError as e:
        if getattr(e, "status", None) == _CUDA_ERROR_UNSUPPORTED_PTX_VERSION:
            raise KernelCompileError(
                f"{name}: the driver (CUDA {driver // 1000}.{driver % 1000 // 10}) cannot compile "
                f"PTX from NVRTC {'.'.join(map(str, nvrtc_version()))}, which is newer; install an "
                "NVRTC no newer than the driver (or a newer driver)."
            ) from e
        raise KernelCompileError(
            f"{name}: the driver failed to compile the PTX for this GPU: {e}"
        ) from e


def _header_digests(dirs: Sequence[str]) -> dict[str, str]:
    out = {}
    for d in dirs:
        with contextlib.suppress(OSError):
            for p in sorted(Path(d).iterdir()):
                if p.suffix in _HEADER_SUFFIXES and p.is_file():
                    out[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def _digest(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=list).encode("utf-8")).hexdigest()


def _cache_read(key: str) -> tuple[dict[str, Any], bytes] | None:
    try:
        blob = (cache_dir() / f"{key}.bin").read_bytes()
    except OSError:
        return None
    try:
        if blob[:8] != _MAGIC:
            return None
        (n,) = struct.unpack("<I", blob[8:12])
        header = json.loads(blob[12 : 12 + n].decode("utf-8"))
        payload = blob[12 + n :]
        if (
            header.get("format") != _FORMAT
            or hashlib.sha256(payload).hexdigest() != header["sha256"]
        ):
            return None
    except (ValueError, KeyError, struct.error, UnicodeDecodeError):
        return None
    return header, payload


def _cache_write(key: str, header: dict[str, Any], payload: bytes) -> None:
    # temporary name and rename: readers see a whole file or none, and writers of the same key
    # write the same content
    d = cache_dir()
    final = d / f"{key}.bin"
    h = json.dumps(
        dict(header, format=_FORMAT, size=len(payload), sha256=hashlib.sha256(payload).hexdigest()),
        sort_keys=True,
    ).encode("utf-8")
    tmp = None
    try:
        import tempfile

        d.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=d, prefix=f".{key[:16]}.", suffix=".tmp")
        with os.fdopen(fd, "wb") as f:
            f.write(_MAGIC + struct.pack("<I", len(h)) + h + payload)
        os.replace(tmp, final)
        tmp = None
    except OSError as e:
        if not final.is_file() and str(d) not in _warned:
            _warned.add(str(d))
            warnings.warn(
                f"The CUDA kernel cache {d} is not writable ({e}); kernels are compiled in every "
                f"process (PETAKIT5D_CACHE_DIR sets the folder).",
                RuntimeWarning,
                stacklevel=3,
            )
    finally:
        if tmp is not None:
            with contextlib.suppress(OSError):
                os.unlink(tmp)
