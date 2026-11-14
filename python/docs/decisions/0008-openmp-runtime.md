# 0008: One OpenMP runtime with the I/O packages (Linux)

Status: proposed (task X05, 2026-11-13). Linux x86-64; measured on two-socket Xeon Gold 6146
(Skylake-SP, 24 cores) and Xeon Gold 6342 (Ice Lake-SP, 24 of 48 cores) nodes, Python 3.12,
cpp-zarr and cpp-tiff 1.7.3, the modules built with GCC 13.2.

## Context

- **Two runtimes in a port process** (record 0007). The cpp-zarr and cpp-tiff wheels bundle one
  libgomp under the name auditwheel gave it, `libgomp-e985bcbb.so.1.0.0` (the first 8 hex digits
  of the original file's SHA-256; it is also the copy's SONAME). The file is AlmaLinux 8's
  libgomp 8.5.0-28, from the manylinux_2_28 image that cibuildwheel 4.1.0 pins (2026.06.04-1).
  The modules linked the toolchain's `libgomp.so.1`, so at 24 CPUs a process had two pools of 23
  threads (70 threads in all, against 47), each spinning after its regions on the cores that the
  other's team needed.
- **The loader maps a library once per SONAME.** A library that needs `libgomp-e985bcbb.so.1.0.0`
  gets the one already loaded under that name, wherever its file is, and loads it from its own
  RUNPATH otherwise.
- **The symbol versions fit.** The modules call 16 OpenMP entry points, of the versions GOMP_1.0,
  GOMP_4.0, GOMP_4.5, OMP_1.0 and OMP_2.0 (built with GCC 13.2, 13.3 and 14.2). The wheels'
  libgomp defines every one of them with its version; its version nodes go up to GOMP_5.0 and
  OMP_5.0, those of GCC 13's libgomp to GOMP_5.1.1 and OMP_5.2.

## Decision

On Linux, with GCC's libgomp, the build links the modules against the libgomp that the cpp-zarr
wheel loads (else cpp-tiff's), puts a copy of that file next to the modules under its SONAME, and
gives the modules the RUNPATH `$ORIGIN`. Whichever loads first, a cpp-zarr or cpp-tiff library or
one of the modules, the others get the same runtime.

- **Finding the file.** `cmake/pk5d_io_openmp.py`, run with the build's Python, loads the
  packages' compiled modules with ctypes (not an import: the packages' import needs NumPy) and
  reads the libgomp files that the loader mapped.
- **Where the build sees cpp-zarr.** `pip install` builds in an isolated environment, which
  gets cpp-zarr as a build requirement on Linux x86-64 (the platform with cpp-zarr wheels).
  `--no-build-isolation`, conda builds and CMake build trees use the environment's own.
- **Checks.**
  - The file must have a SONAME other than `libgomp.so.1`: under that name the loader shares the
    runtime already.
  - It must link `cmake/pk5d_openmp_probe.cpp` as an executable: the constructs and calls the
    kernels use (parallel regions, static and dynamic loops with `int` and `size_t` indices,
    collapse, single, barrier, max/min/or reductions, `omp_get_wtime` and the thread queries).
  - After the link of each module, `cmake/Pk5dCheckOpenMP.cmake` checks that the module needs no
    other OpenMP runtime and imports no OpenMP symbol that the runtime lacks: the linker leaves
    such a symbol without a version, and the module would not load.
- **`PK5D_OPENMP_SHARE`.**
  - `AUTO`, the default: falls back to the toolchain's runtime, with a status line, when no
    package's libgomp is found or the checks fail.
  - `OFF`: the toolchain's runtime.
  - A file: that libgomp, and an error if it can't be used.
- **`build_info()["openmp_shared"]`** names the shared runtime (its SONAME), or is `"off"`, or
  empty (the toolchain's runtime, or not Linux).
- **Wheels.** Built in any manylinux image with cpp-zarr in the build environment, the wheel
  carries `petakit5d/_ext/libgomp-e985bcbb.so.1.0.0`, and auditwheel leaves it there: it is
  inside the wheel. Tested with auditwheel 6.7.0 in the image of the cpp-zarr wheels, installed
  next to the 1.7.3 wheels: one runtime in either import order.
- **Tests:** `tests/test_openmp_runtime.py`. It maps the runtimes in a fresh process, with our
  modules first and with the I/O packages first, each calling a region of 2 threads. It fails
  when two runtimes are mapped, and skips where the wheels aren't installed.

## Consequences

- **One pool of threads per calling thread** serves the three libraries. The `GOMP_*` and
  `OMP_*` environment acts on that one runtime, and libgomp's throttled spin counts all of its
  threads. The thread counts are no one else's: the kernels pass `num_threads` to every region,
  and a cpp-zarr call with 3 threads left the default team size (`omp_get_max_threads()`) at 24.
- **The runtime is the wheels' libgomp**, AlmaLinux 8's 8.5.0 package, not the toolchain's.
  Kernels that need OpenMP entry points newer than GOMP_5.0 and OMP_5.0 fail the build's check;
  they need `OFF`, or a newer runtime in the I/O wheels.
- **Drift.** A release of cpp-zarr or cpp-tiff with another libgomp file gives another name: the
  current manylinux_2_28 image has libgomp 8.5.0-31, which auditwheel names
  `libgomp-cb35b7ff.so.1.0.0`. A port built against the old file still loads, with its own copy,
  but has two runtimes again until it is rebuilt against the new one; the test reports it. The
  same split would come between cpp-zarr and cpp-tiff if they were released from different
  images.
- **License.** libgomp is under the GPL 3 with the GCC Runtime Library Exception, which allows
  shipping it with the modules; auditwheel ships it so in every Linux wheel that uses OpenMP,
  among them cpp-zarr's and cpp-tiff's with this same file.
- **macOS and Windows: unchanged.** The modules use LLVM's libomp there, the wheels GCC's libgomp
  (macOS: a copy in each package, so two; Windows: MinGW's `libgomp-1.dll`), which doesn't spin.
  libomp spins up to `KMP_BLOCKTIME` (200 ms) after each region.

## Measurements

Task X05's report has the tables and jobs. Four settings, interleaved in 5 rounds on one node,
every process under `numactl --interleave=all`: **one** runtime (this build) or **two** (the
toolchain's), each with libgomp's own spin or with record 0007's `GOMP_SPINCOUNT=100000`.

S07's feather blending of 4 volumes, as a ratio to MATLAB's time on the same node:

| threads | Skylake-SP: two, own spin / two, 100000 / one, own spin / one, 100000 | Ice Lake-SP: the same |
|---|---|---|
| 24 | 2.22 / 1.36 / 1.02 / 0.96 | 1.53 / 1.06 / 0.96 / 0.90 |
| 4 | 1.16 / 0.97 / 0.89 / 0.89 | 1.00 / 0.94 / 0.89 / 0.88 |
| 1 | 0.83-0.84 in every setting | 0.84-0.85 in every setting |

X02's microbenchmark at 24 threads (4 threads in brackets), as a ratio to two runtimes with
100000, the port after record 0007:

| | Skylake-SP: one, own spin / one, 100000 | Ice Lake-SP: the same |
|---|---|---|
| a read-copy-MIP-write loop of 512^3 blocks | 0.82 / 0.82 (0.98 / 0.98) | 0.91 / 0.92 (0.99 / 0.99) |
| a y MIP right after a read | 0.34 / 0.34 | 0.94 / 0.97 |
| a read right after a kernel of ours | 0.82 / 0.82 | 0.95 / 0.95 |
| a 16 MB MIP after an idle gap of 2 ms, us | 92 / 92 | 81 / 337 |
| the same after a gap of 10 ms, us | 92 / 236 | 338 / 336 |
| CPU time the idle pool spins in the 60 ms after a read, ms | 240 / 61 | 83 / 17 |

- One runtime closes S07's gap at 24 threads: within 1.05 of MATLAB with libgomp's own spin.
- `GOMP_SPINCOUNT=100000` still takes 6-8% off S07 at 24 threads with one runtime (paired
  ratios 1.06 and 1.08), changes nothing at 4 and 1 threads or in the loop, and spins a quarter
  of the CPU time. Its cost is the one record 0007 measured: a wake-up for a short kernel after
  an idle gap longer than its spin and shorter than libgomp's. So record 0007 stays.
- At 1 thread no team exists; differences of up to 6% between processes there don't follow the
  setting, as X02 found.

## Alternatives

- **The same manylinux image, without the copy** (`PK5D_OPENMP_SHARE=OFF` in the image of the
  cpp-zarr wheels, gcc-toolset-14, auditwheel 6.7.0): gcc-toolset-14 links the image's
  `/usr/lib64/libgomp.so.1`, which auditwheel bundles under the same name, byte for byte the
  wheels' file. One runtime, but only for that image: the next one gives another name, silently.
  Source builds aren't covered.
- **Loading the wheels' libgomp `RTLD_GLOBAL` before the modules** (X02's emulation): binds every
  OpenMP library loaded later in the process to it, while the modules still load the toolchain's
  `libgomp.so.1` that they need.
- **Rewriting the modules' needed library at install time** (patchelf): the same result as the
  build's link, with one more tool at install time.
