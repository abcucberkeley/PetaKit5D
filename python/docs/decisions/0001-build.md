# 0001: Build system, compilers and flags

Status: accepted for Linux; macOS and Windows chosen, to be confirmed by CI (see the end).

## Context

The kernels move from MATLAB mex files to Python extension modules with their loops unchanged,
and must give the mex results bit for bit on Linux x86-64 and the same results on macOS ARM and
Windows wherever the C library allows. Several tasks add kernels at the same time, so no task
may need to edit a shared build file. The package must import and run on machines with or
without a GPU, and its import must stay cheap (start-up time is a goal of the port).

## Decision

**Build:** scikit-build-core, CMake (3.26 or newer) and nanobind, C++17 and C11.
`pip install .` builds a wheel; developers build a CMake tree outside the repository and set
`PETAKIT5D_EXT_DIR` to its `petakit5d/_ext` folder, which `petakit5d._ext` searches before its
own folder (`PYTHONPATH=python/src` for the Python sources).

**One module per kernel group:** each folder `cpp/<group>/` with a `CMakeLists.txt` is found by
`cpp/CMakeLists.txt` (also a folder added after configuring) and calls
`pk5d_add_extension(<group> <sources>...)`, giving `petakit5d._ext.<group>` (target
`pk5d_ext_<group>`). Options: `AVX2`, `AVX2_SOURCES`, `NO_OPENMP`, `LIBRARIES`, `INCLUDE_DIRS`,
`DEFINITIONS`. `cpp/common/pk5d/common.hpp` has the array types (`fvolume<T>`: 3D Fortran
order), `affinity_cpu_count`, `check_threads`, and `PK5D_BINDINGS` for groups whose bindings
are spread over several files (the primitives of `cpp/core/ops/` register into `core` this way).
`cpp/cuda/` holds GPU kernel sources, which are not compiled at build time: they are copied
next to the modules (`petakit5d/_ext/cuda/`) and compiled at run time
(`petakit5d._ext.cuda_source_dirs()`).

**Flags of every kernel** (`pk5d_flags`):

| | GCC, Clang | MSVC |
|---|---|---|
| optimization | `-O3` (the mex: `-O3`) | `/O2` |
| floating point | `-ffp-contract=off` | `/fp:precise` (VS 2022 does not contract under it) |
| OpenMP | `-fopenmp` (Apple Clang: `-Xpreprocessor -fopenmp` + libomp) | `/openmp:llvm` |
| AVX2 | `-mavx2` for groups or sources declared `AVX2`, x86-64 only | `/arch:AVX2` |

- No fused multiply-adds anywhere: on x86-64 that is what the mex do (built without FMA); on
  ARM, where GCC and Clang contract by default, it gives the x86 results.
- Never `-march=native`, fast-math, `-Ofast`, `-mfma`, `-ffp-contract=fast|on`, `/fp:fast` or
  `/fp:contract`: configuring fails if `CFLAGS`/`CXXFLAGS` or a toolchain file brings one in.
- nanobind's default `-Os` for module code is turned off (`NOMINSIZE`).
- Linux: `-static-libstdc++ -static-libgcc` and `--exclude-libs,ALL`. NumPy, SciPy and the I/O
  packages load the system's `libstdc++.so.6` first, and the loader reuses a loaded library of
  that name whatever a module's RPATH says, so a module built with a newer GCC than the system's
  failed to import after NumPy (`CXXABI_1.3.13 not found`). The OpenMP runtime stays shared:
  the cpp-zarr wheel's libgomp, copied next to the modules, so that a process has one runtime
  (`docs/decisions/0008-openmp-runtime.md`), else the toolchain's, whose folder goes on the
  modules' RPATH when it isn't a system folder.

**Compilers per system:**

- **Linux:** GCC with libgomp (reference: GCC 13.2 on Rocky 8 x86-64; CI: Ubuntu 24.04 GCC 13).
- **macOS ARM:** Apple Clang with libomp (`brew install libomp`, or `OpenMP_ROOT`). The mex use
  Homebrew GCC with static libgomp; for Python modules Apple Clang is the mainstream choice: it
  matches the SDK and libc++ that CPython and the other packages use. GCC remains possible
  (`CC=gcc-14 CXX=g++-14`). With `-ffp-contract=off` both give IEEE results without
  contraction; elementary functions (libm) differ by platform in any case.
- **Windows:** MSVC (Visual Studio 2022) with `/openmp:llvm`. The kernels need OpenMP 3.1
  (`collapse`, `min`/`max` reductions, unsigned loop indices), which MSVC's default `/openmp`
  (2.0) lacks. Its runtime `libomp140.<arch>.dll` is not part of Windows: the build finds it in
  the toolset's redistributable folders and puts it next to the modules (the DLL folder is
  searched when Python loads an extension). MinGW-w64, which the mex use, was not chosen:
  CPython and the binary packages are built with MSVC, nanobind and scikit-build-core don't
  support MinGW as a first-class target, and MinGW modules would need its C++ and OpenMP runtimes
  linked statically or shipped.

**Run-time checks** (`petakit5d._ext.core`, tested in `tests/test_smoke.py` on every system):
`build_info()` (compiler, flags, OpenMP and the shared runtime), `fma_contracted()` (a multiply-add whose result
differs when fused), `subnormals_flushed()`, `parallel_threads(n)`, and the CPU counts
(`affinity_cpu_count`, `omp_num_procs`, `omp_max_threads`: inside a Slurm job with `-c 2` all
three are 2).

**Environments:** ranges in `pyproject.toml` (Python 3.11 or newer; NumPy 2); a hash-pinned lock
for Python 3.12 on Linux x86-64 (`requirements/linux-x86_64-py312.txt`, from `requirements/*.in`
with the `uv pip compile` command in its header) for reproducible runs and the Linux 3.12 CI job.
Building needs a Python with its development headers; where the system Python has none, a
standalone CPython build (for example one installed by `uv python install`) works.

**CI** (`.github/workflows/python.yml`, only for changes under `python/`): Ubuntu 24.04,
macOS 15 (Apple silicon) and Windows Server 2025, Python 3.12 (Linux: the lock file) and 3.14
(newest packages); build with `pip install -v`, print `build_info()`, run the tests; ruff.

## Consequences

- Modules with `AVX2` don't load on x86-64 CPUs without AVX2, as the mex.
- OpenMP's static schedule may split loops differently between runtimes (libgomp, libomp);
  kernel results must not depend on the split (fixed-order sums), as the port's rules require.
- On macOS a second libomp from another package (loaded from another path) can conflict with
  ours; wheels (E07) bundle and rename the runtimes (delocate, delvewheel, auditwheel).
- nanobind is linked statically into each module; kernels never pass C++ objects between
  modules.

## What CI showed

Pending: the first push of this skeleton. To fill in from the CI run: compilers and OpenMP
versions per system (`build_info()`), the libomp140 DLL found on Windows, test results.
