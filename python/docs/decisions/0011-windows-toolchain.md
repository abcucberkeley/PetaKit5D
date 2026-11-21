# 0011: Windows: MinGW-w64's GCC beside MSVC

Status: proposed (task W02, 2026-11-21): both toolchains build and are tested; the user chooses the
one the release wheels and the Windows CI use, on the measurements below. Decision 0001 chose MSVC
for Windows; this record adds MinGW-w64 and does not remove MSVC.

## Context

- **The user, 2026-11-20:** "I think I would prefer MinGW instead especially if it can be faster".
  The lab's cpp-zarr and cpp-tiff build their Windows wheels with MinGW-w64's GCC: the MinGW-builds
  project's GCC 12.2.0 (UCRT, posix threads, SEH, binutils 2.39), with GCC's libgomp.
- **MSVC builds the port** with `/openmp:llvm` (LLVM's libomp140) and without GCC's vector
  extensions and target attributes: the MIP, flat-field, resize, unmix and filter kernels run their
  baseline loops (no AVX2 dispatch), morphology and the median filter work lane by lane, and
  pocketfft transforms one line at a time ([platforms](../platforms.md), "Toolchains").
- **Two OpenMP runtimes in an MSVC process:** the port's libomp140 and the I/O packages'
  libgomp-1.dll, each with a pool of threads; libomp spins 200 ms after each region.

## Decision

The build takes MinGW-w64's GCC as a second Windows toolchain: `CC=gcc CXX=g++` with the Ninja
generator (`CMAKE_GENERATOR=Ninja`), for `pip install`, cibuildwheel and CMake build trees alike.
`tools/get_mingw.sh` fetches the toolchain of the libraries' wheels by URL and SHA-256 (CI uses it).

**What the build does with it** (`cmake/Pk5dExtension.cmake`, `cmake/Pk5dOpenMP.cmake`):

- GCC's flags, as on Linux (`-O3 -ffp-contract=off`, `-mavx2` for the AVX2 groups and the run-time
  dispatch), with Linux's instruction set and tuning (`-march=x86-64 -mtune=generic`: MinGW-builds
  default to nocona and core2) unless `CFLAGS`/`CXXFLAGS` name their own.
- **Unaligned vector moves** (`-Wa,-muse-unaligned-vector-move`, binutils 2.38 and newer, required):
  GCC can't align the stack beyond 16 bytes with Windows' SEH but spills 32-byte AVX values with
  aligned moves (GCC bug 54412), which fault on a stack slot that isn't 32-byte aligned.
- libstdc++ and libgcc linked statically, as in the libraries; winpthreads stays a DLL, which
  libgomp loads too.
- **One OpenMP runtime with the I/O packages:** the toolchain's `libgomp-1.dll`, with the DLLs it
  loads (`libwinpthread-1.dll`, `libgcc_s_seh-1.dll`), is copied next to the modules under its own
  name (`petakit5d/_ext/`, also in the wheels: delvewheel leaves them). The cpp-zarr and cpp-tiff
  libraries load `libgomp-1.dll` of the same toolchain by that name, and Windows loads a DLL of a
  name once per process, so the three share one runtime whichever is imported first
  (`tests/test_openmp_runtime.py`). The three files are byte for byte those of the cpp-tiff 1.8.0
  wheels (cpp-zarr's `libgomp-1.dll` differs only in the import delvewheel renamed in it).
- **The exact add-on** (FFTW) builds with it too: DLLs named as MSVC names them, FFTW's API
  exported alone, FFTW's own aligned allocator (`WITH_OUR_MALLOC`, FFTW's advice for MinGW), libgcc
  linked in; the DLLs load only the system's.

**What MinGW-w64 needed in the code** (measured on Windows 10 with GCC 12.2.0; all of it compiles to
nothing elsewhere, and the Linux, macOS and MSVC builds compile the code they compiled before):

- **The x87 unit's precision differs by thread.** A Python built with MSVC runs its threads with
  the x87 precision of 53 bits (Windows' default, control word 0x027F); winpthreads gives every
  thread it starts, OpenMP's among them, MinGW's 64 bits (0x037F). MinGW-w64's own C math functions
  (libmingwex) compute on the x87 unit, so the same call gave other bits on other threads: exp in
  17% of 20,000 random arguments, pow in 48%, erfc in 4%, and the detection fits' Anderson-Darling
  statistic changed with the thread count. `cpp/common/mingw_libm.c` wraps them (`--wrap`): each
  runs at 64 bits, MinGW's and Linux's precision, where it comes closest to glibc (exp: 18 of the
  20,000 values differ from glibc's, the UCRT's 93; sin 45 and 485; cos 10 and 867; erfc 905 and
  8,768), and an invalid operation gives x86's default NaN, as glibc and the UCRT do (MinGW-w64's
  `log(-1)` and `sin(Inf)` were positive NaNs). pocketfft's twiddle factors start from an angle
  computed in long double: its plans are made at 64 bits too (`pk5d/x87.hpp`), so they are Linux's.
- **MinGW-w64's `fma` is not fused:** it differed from the instruction in 11% (double) and 22%
  (float) of random near-cancelling inputs; the wrapper calls the UCRT's, which equals the
  instruction. The filters' fused taps use it where the CPU has no FMA.
- **MinGW-w64 10's `round`** gives 1 for 0.49999999999999994 (the libraries met it, cpp-zarr
  22049b9): the wrappers round halves away from zero exactly.
- **thread_local objects with destructors** are destroyed after emulated TLS has freed them when a
  thread ends (OpenMP's threads end as their teams change): winpthreads runs emutls' key destructor
  before the module's TLS callback runs theirs. pocketfft's backend kept a line buffer per thread
  so, and freed through freed memory (heap corruption after about 20 transforms). On MinGW-w64 the
  buffers live under pthread keys instead.
- **pocketfft's thread pool** (the native half spectra) is a static object. At the process's exit
  Windows ends every other thread, then runs the modules' static destructors, and winpthreads'
  `pthread_cond_destroy` waited forever for the ended threads: the process never exited (and can't
  be ended: TerminateProcess refuses a process that is exiting). The pool's threads are ended at
  the interpreter's exit (`Py_AtExit`) instead.

## The same bits

[platforms](../platforms.md), "Measured: the kernels' bits on each system". With MinGW-w64 on
Windows 10 (Xeon Gold 6442Y, AVX-512), Python 3.12: every pinned case of
`tests/platform/test_kernel_hashes.py` has Linux's digest, and the cases recorded, not pinned, gave
Linux's digests for the FFTs; the platform and smoke tests: 194 passed, 47 skipped (the skips of
Linux x86-64's builds). On GitHub's windows-2025 runners: see "What CI showed".

## Speed

To be filled in from the interleaved runs of the two wheels on one Windows 10 machine (8 vCPUs).

## What CI showed

To be filled in from the wheels-trial runs.

## Consequences

- MinGW-w64 builds use the filters' FMA dispatch, so the fused chains' NaN rules hold on Windows
  as on Linux (`tests/platform/test_nan_rules.py` runs them).
- The Windows wheels carry GCC's and MinGW-w64's runtime libraries, and the toolchain's notices
  of them as it ships them: the build copies them from the toolchain's `licenses/` folder into
  `petakit5d/_ext/licenses/` (the add-on's: `petakit5d_fftw/licenses/`), checking their SHA-256
  where the toolchain is MinGW-builds' GCC 12.2.0, and the release checks require them in a wheel
  whose binaries GNU ld linked, at those exact paths, which they leave out of their search for
  folders and site names (THIRD_PARTY_NOTICES.md).
- A process can mix a MinGW-built port with other packages' `libgomp-1.dll` of another GCC: the
  first loaded serves all, and libgomp is backward compatible (the modules call the entry points of
  GCC 4.9 and later).
- New C math functions in the kernels: `cmake/Pk5dExtension.cmake`'s `PK5D_MINGW_LIBM` lists the
  wrapped ones; a function outside it, other than the exact ones (floor, ceil, trunc, nearbyint,
  rint, ldexp, frexp, modf, fmod, fabs, copysign, sqrt), needs adding there.
