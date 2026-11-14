# 0007: OpenMP's spin after a parallel region (`GOMP_SPINCOUNT`)

Status: proposed (task X02, 2026-11-13). Measured on Linux x86-64, on two-socket Xeon Gold 6146
(Skylake-SP, 24 cores) and Xeon Gold 6342 (Ice Lake-SP, 24 of 48 cores) nodes, Python 3.12,
cpp-zarr and cpp-tiff 1.7.3, the extension modules built with GCC 13.2.

## Context

- **Two OpenMP runtimes in one process.** On Linux the cpp-zarr and cpp-tiff wheels bundle the
  same libgomp file under one hashed name, which the dynamic linker loads once for both. The
  port's extension modules use the toolchain's `libgomp.so.1`. At 24 CPUs a process holds two
  pools of 23 worker threads.
- **The spin.** After each parallel region a libgomp thread spins `GOMP_SPINCOUNT` times on the
  CPU's pause instruction, then sleeps: 300,000 by default, about 11 ms on Skylake-SP (36 ns a
  pause) and 3 ms on Ice Lake-SP (10 ns). Each runtime counts only its own threads, so neither
  shortens its spin while the other's team works.
- **The cost.** A cpp-zarr or cpp-tiff call next to one of the port's kernels runs both teams on
  the same cores while the idle one spins: feather blending took 2.3 times MATLAB's time at 24
  threads (S07), a MIP right after a read 4.6 times as long as with a short spin. MATLAB's mex
  files share one runtime (and MATLAB runs its own with `KMP_BLOCKTIME=1`, 1 ms).
- **The cost of too short a spin.** A sleeping team takes time to start again: about 0.13 ms at
  24 threads on Skylake-SP, up to 0.5 ms on Ice Lake-SP (futex wake-ups, the cores' C6 exit).
  When the spin is shorter than that, threads that wait for a late teammate fall asleep too and
  every region starts slowly. On Ice Lake, with 10,000 spins, the port's 2-6 ms MIP and pooling
  kernels took 1.04-1.36 times as long, a 2x2x2 one 37 times, a small cpp-zarr read up to twice;
  with 1,000 spins or `OMP_WAIT_POLICY=PASSIVE` even back-to-back calls took 7 times as long.
  On Skylake the same counts are harmless.

## Decision

`import petakit5d` sets `GOMP_SPINCOUNT=100000` in the process's environment unless
`GOMP_SPINCOUNT` or `OMP_WAIT_POLICY` is already set: about 3.6 ms of spin on Skylake-SP and
1 ms on Ice Lake-SP, a third of libgomp's default and still longer than a team's wake-up. The
package imports nothing at that point; the extension modules and the I/O packages are imported
lazily, so every runtime that loads afterwards reads it.

- **The user's setting wins.** `GOMP_SPINCOUNT` (300000 is libgomp's own default) or
  `OMP_WAIT_POLICY` set before Python starts, or in `os.environ` before the first `import
  petakit5d`, is left alone (libgomp lets `GOMP_SPINCOUNT` override the policy's spin, so with
  `OMP_WAIT_POLICY` set nothing is added).
- **Only runtimes loaded later.** libgomp reads its environment once, when it is loaded. In a
  process that imported cppzarr or cpptiff before petakit5d, their runtime keeps libgomp's
  default; import petakit5d first, or set the variable in the shell.
- **Process-wide.** Other libgomp-based libraries loaded after petakit5d (scikit-learn's,
  PyTorch's CPU threads) spin 100,000 times too, still 1-3.6 ms.
- **Child processes** inherit it: the job framework's local workers and isolation harness
  (`dict(os.environ)`), subprocesses, and Slurm jobs (sbatch exports the submitting
  environment); each of them also imports petakit5d itself. `PETAKIT5D_GOMP_SPINCOUNT` notes
  the value the port set, so that a program that isn't the port can be started without it
  (`petakit5d._without_openmp_spin`): the benchmark harness starts MATLAB so, since
  PetaKit5D's Linux mex files use libgomp too and MATLAB's times must stay those of MATLAB as
  the user's shell starts it. A benchmark of the port with libgomp's own spin sets
  `GOMP_SPINCOUNT=300000` in its variant.
- **macOS and Windows.** The wheels' libgomp there (GCC's macOS and MinGW builds) blocks at
  once and never spins, and the port's modules use LLVM's libomp, which ignores
  `GOMP_SPINCOUNT`: nothing changes there. libomp's own spin is `KMP_BLOCKTIME` (200 ms by
  default); whether a shorter one helps there is to be measured on those systems.

## Measurements

X02's report has every table, with its jobs. At 24 threads, as ratios to libgomp's default in
the same job (4 threads in brackets, where measured):

| | Skylake-SP | Ice Lake-SP |
|---|---|---|
| S07's feather blending of 4 volumes | 0.61 (0.83-0.85) | 0.70 (0.94) |
| T02's MIPs of a tiff tile, saved (`save_mip_tiff`) | 0.71 (0.96) | |
| a read-copy-MIP-write loop of 512^3 blocks | 0.79 (0.96) | 0.81-0.87 (0.95-1.00) |
| T03's crop of a large zarr in 1024^3 batches (`crop_zarr`) | | 0.85 (0.98) |
| D04, S06: deskew/rotate of datasets and large files, distance maps | | 0.89-0.99 (1.00) |
| R03, R05, R07: decon in memory, of a file, of a large zarr | | 0.94-0.99 (0.98-1.00) |
| K02's 2-6 ms MIP and pooling kernels | 0.99-1.01 | 0.95-1.09 |
| small cpp-zarr and cpp-tiff calls (0.2-3 ms) | 0.94-1.00 | 1.01-1.05 |
| a 16 MB MIP (80 us) back to back, or after an idle gap of 0.5 ms | 1.00-1.01 | 0.99-1.02 |
| the same after an idle gap of 2 ms / 10 ms | 1.01 / 2.9 | 2.6-3.0 / 0.98 |

At 1 thread nothing changes (no team exists). The last row is what any shorter spin costs: a
kernel of well under a millisecond, called again after an idle gap longer than the new spin but
shorter than libgomp's default, pays the team's wake-up (0.13-0.16 ms here).

## Alternatives

- **A shorter count (`GOMP_SPINCOUNT=10000` or 1000, or `OMP_WAIT_POLICY=PASSIVE`)**: best on
  Skylake (S07 0.43 of the default with 10000: 0.97 of MATLAB's time), but on Ice Lake the
  port's short kernels and the wheels' small calls fall into the slow regime above.
- **Two counts** (10,000 for the wheels' runtime, set around their first import by a meta path
  finder, 100,000 for the rest): S07 0.47 of the default on Skylake and 0.63 on Ice Lake, but on
  Ice Lake the wheels' small calls (under a few ms at 24 threads) run up to twice as long.
- **One shared runtime.** A port wheel repaired in the same manylinux image as the cpp-zarr and
  cpp-tiff wheels would bundle the identical libgomp under the identical name, and the dynamic
  linker would load one runtime for all three, with no contention to begin with (emulated: S07
  0.45 of the default on Skylake with libgomp's own spin, and the port's kernels run as fast on
  the wheels' older libgomp). That holds only while all three use the same toolchain file, only
  for those Linux wheels (source builds use the toolchain's `libgomp.so.1`, delocate on macOS
  keeps a copy per package, Windows mixes MinGW's libgomp with MSVC's libomp), and any drift
  splits it again silently. Loading the wheels' libgomp `RTLD_GLOBAL` before the port's modules
  does the same at run time, but binds every later OpenMP library of the process to it.
- **A spin by time**: the count's duration and the wake-up it must outlast both differ by CPU,
  so a count calibrated per machine (in an OpenMP-free extension, before the first OpenMP module
  loads) could give Skylake's short spin and Ice Lake's long one.
- **Releasing a pool after each call** (`omp_pause_resource_all`): its threads are created again
  by the next call, which costs about as much as the spin.
