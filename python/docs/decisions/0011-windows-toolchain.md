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
  (`tests/test_openmp_runtime.py`). The three files are byte for byte those of the cpp-tiff wheels
  (1.7.3 and 1.8.0; cpp-zarr's `libgomp-1.dll` differs only in the import delvewheel renamed in it).
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
- **`nearbyint` and `ldexp`** are exact at any precision, but MinGW-w64's run the x87 unit
  (`nearbyint` sets and restores its control word around `frndint`): 27 and 5.8 ns a call on
  Windows 10, against 1.4 for their SSE2 forms, which the wrappers use (the same bits as
  MinGW-w64's for a million arguments, subnormal results among them). The detection fits'
  correctly rounded exp calls them four times a value. (The UCRT's take 17 and 19 ns.)
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

[platforms](../platforms.md), "Measured: the kernels' bits on each system". On Windows 10 (Xeon
Gold 6442Y, AVX-512), with the two wheels of one commit built by the wheel workflow, Python 3.12:

- **MinGW-w64:** all 65 pinned cases of `tests/platform/test_kernel_hashes.py` give Linux's
  digests (the platform, smoke and OpenMP-runtime tests: 196 passed, 47 skipped, as Linux
  x86-64's builds skip; the same with the cp311 wheel on Python 3.11). Of the 9 cases recorded,
  not pinned, 6 give Linux's bits too: the FFTs (`fft.matlab_fft`, `fft.half_spectra`,
  `specials.fft_and_decon`), the normalized cross-correlation, the deconvolution's OTF and
  Richardson-Lucy, and the least squares. The other 3 differ where the C library or NumPy
  computes: the Anderson-Darling statistics (`erf`, `erfc`, `log`: 3 and 3 values) and edgetaper
  (NumPy's `np.fft` and complex `np.abs`).
- **MSVC:** 64 of the 65 (`detection.cr_exp` fails: Windows 10's `ldexp`, [platforms](../platforms.md),
  "What differs by system"), 187 passed, 55 skipped, 1 failed; 2 of the 9 recorded cases give
  Linux's bits.

With the SSE2 `nearbyint` and `ldexp` (above) the MinGW-w64 wheel gives the bits of the one
without them in all 74 cases, and its smoke, OpenMP-runtime, platform and detection tests pass
there (657 passed, 151 skipped). On GitHub's windows-2025 runners: see "What CI showed".

## Speed

The two wheels of one commit, built by the wheel workflow, on one Windows 10 machine (8 vCPUs of a
Xeon Gold 6442Y), measured 2026-11-21 through `petakit5d.testing.bench`: one process per round,
the wheels' processes in alternating order, a warm-up round and 5 rounds, each operation timed 3
to 5 times a round at 1 thread and at 8. The detection fits were measured again with the SSE2
`nearbyint` and `ldexp` above (a wheel of this tree built there with the same toolchain), against
both wheels in the same rounds. The table gives the
medians in seconds and MinGW-w64's time over MSVC's (the median of the rounds' paired ratios; "~"
where the two wheels' round medians overlap).

| family | operation | 1 thread: MSVC s | MinGW-w64 s | MinGW/MSVC | 8 threads: MSVC s | MinGW-w64 s | MinGW/MSVC |
|---|---|---|---|---|---|---|---|
| decon | Richardson-Lucy, 10 iterations, 81 x 320 x 320, PSF 41 x 31 x 31 | 4.41 | 1.89 | 0.43 | 0.795 | 0.404 | 0.52 |
| decon | OMW, 2 iterations, same data | 0.827 | 0.324 | 0.39 | 0.153 | 0.0699 | 0.49 |
| deskew/rotate | deskew, uint16, 301 x 256 x 512 | 0.786 | 0.53 | 0.67 | 0.118 | 0.0906 | 0.77 |
| deskew/rotate | deskew and rotate, same data | 0.252 | 0.232 | ~0.92 | 0.0713 | 0.0707 | ~1.01 |
| filters | imgaussfilt3, sigma 1, single, 400 x 512 x 512 | 6.31 | 0.938 | 0.15 | 0.907 | 0.267 | 0.28 |
| filters | imgaussfilt3, sigma 2.5, uint16, same size | 13.6 | 1.94 | 0.14 | 1.97 | 0.535 | 0.26 |
| filters | imgaussfilt, sigma 100, 1024^2 (spatial) | 3.15 | 0.157 | 0.05 | 0.424 | 0.0319 | 0.08 |
| filters | imgaussfilt, sigma 100, 2048^2 (frequency domain) | 0.324 | 0.165 | 0.49 | 0.107 | 0.0838 | 0.80 |
| FFT | fft along one axis, single, 96 x 512 x 512 | 0.212 | 0.157 | 0.77 | 0.0452 | 0.0368 | 0.83 |
| FFT | fftn, single, 96 x 512 x 512 | 0.559 | 0.373 | 0.65 | 0.0985 | 0.0934 | ~0.97 |
| FFT | ifftn, complex single, 96 x 512 x 512 | 0.699 | 0.455 | 0.67 | 0.117 | 0.108 | ~0.95 |
| FFT | half spectrum, forward, single, 96 x 512 x 512 | 0.338 | 0.146 | 0.44 | 0.0524 | 0.0389 | 0.73 |
| FFT | half spectrum, inverse | 0.362 | 0.17 | 0.47 | 0.0844 | 0.0653 | 0.77 |
| resize | imresize3 to half, cubic, uint16 500 x 512 x 1800 | 2.9 | 2.12 | 0.73 | 0.419 | 0.331 | 0.80 |
| resize | imresize3 to half, linear | 2.32 | 1.43 | 0.61 | 0.358 | 0.253 | 0.70 |
| resize | imresize3, z twice, cubic | 5.22 | 3.19 | 0.62 | 0.717 | 0.464 | 0.64 |
| resize | imresize3, z twice, linear | 4.88 | 2.33 | 0.48 | 0.67 | 0.359 | 0.54 |
| resize | imresize3, z half, linear | 1.31 | 0.779 | 0.59 | 0.194 | 0.124 | 0.65 |
| resize | imresize3_fast to half, cubic | 2.87 | 2.14 | 0.75 | 0.429 | 0.331 | 0.80 |
| resize | imresize3_fast, z twice, cubic | 5.28 | 3.15 | 0.59 | 0.747 | 0.471 | 0.61 |
| detection fits | fitGaussian3D mex, xyzAsrc, 3000 windows of 9^3 | 1.36 | 1.36 | ~1.00 | 0.198 | 0.196 | ~1.00 |
| detection fits | fitGaussian3D mex, xyzAc | 1.35 | 1.35 | ~1.01 | 0.191 | 0.19 | ~1.00 |
| detection fits | fitGaussian3D mex, Ac | 0.354 | 0.435 | 1.23 | 0.0695 | 0.0718 | ~1.04 |
| detection fits | CD_Newton fit | 1.07 | 1.05 | ~1.00 | 0.151 | 0.156 | ~1.03 |
| morphology | imerode, disk 11, 500 x 512 x 1800 logical | 0.462 | 0.525 | 1.17 | 0.0736 | 0.0921 | ~1.21 |
| morphology | bwdist, 2048^2 | 0.0511 | 0.0425 | 0.85 | 0.0079 | 0.0076 | ~0.94 |
| morphology | medfilt3 3 x 3 x 3, single, 101 x 128 x 128 | 0.0441 | 0.226 | 5.29 | 0.0068 | 0.0332 | 4.70 |
| morphology | bwlabeln 26, 200 x 512 x 512 | 0.779 | 0.745 | ~0.95 | 0.631 | 0.65 | ~1.03 |
| morphology | bwconncomp 26, same | 0.869 | 0.789 | ~0.94 | 0.332 | 0.326 | ~1.01 |

- **Faster with MinGW-w64:** the kernels where GCC's vector extensions and the run-time
  instruction-set dispatch apply, which MSVC builds without: the Gaussian filters' fused taps
  (the FMA instruction, where MSVC's build calls the C library's `fma` for each tap: up to 20
  times as fast), the deconvolutions (2 to 2.6 times at 1 thread, 2 at 8), pocketfft's vector lines
  (MSVC transforms one line at a time), the resize kernels and the deskew.
- **Slower with MinGW-w64:** the 3 x 3 x 3 median filter (5 times) and binary erosion (1.2),
  whose GCC code is Linux's: GCC compiles the median network's 16-lane vectors into scalar
  `minss` and `maxss`, while MSVC vectorizes the lane-by-lane loops it gets instead. A change to
  those kernels would make Linux faster too (outside this record).
- **The detection fits:** even, but for the fits of amplitude and background alone at 1 thread
  (1.23), whose time is mostly the Anderson-Darling test: MinGW-w64's `erf` and `log` take 62 and
  29 ns a call at 64 bits, the UCRT's 18 and 4. The UCRT's would make those fits faster than
  MSVC's build, but give the UCRT's bits there instead of results closer to glibc's (in the two
  recorded cases of the statistic, 59 values differ from Linux's with MSVC's build, 6 with
  MinGW-w64's), which this build does not do. Without the SSE2 `nearbyint` and `ldexp` the fits
  took 1.11 to 1.35 times MSVC's time at 1 thread.

## What CI showed

The trial of 2026-11-21 (the wheel and Python workflows on GitHub's runners, both toolchains):

- **Wheels:** every job passed. On windows-2025, MinGW-w64 built, repaired, audited, tested and
  checked the abi3 and cp311 wheels in 6.5 minutes (MSVC: 7.9) and the add-on in 3.7 (MSVC: 7.5).
  The MinGW-w64 wheels carry `libgomp-1.dll`, `libwinpthread-1.dll` and `libgcc_s_seh-1.dll` with
  the SHA-256 of the cpp-tiff wheels' copies; check_wheel found nothing in any wheel.
- **Python** (the source build with pip and the whole CPU suite, windows-2025): every Windows job
  passed, Python 3.12 and 3.14: MinGW-w64 27,432 passed, 6,019 skipped, 15 xfailed (31 to 32
  minutes of tests); MSVC 27,421 passed, 6,030 skipped, 15 xfailed (32 to 34). MinGW-w64 runs 11
  tests that MSVC's build skips: the instruction-set dispatch's (filters, resize, unmix), the NaN
  rules of the filters' fused chains, and the one `libgomp-1.dll` with the I/O packages.
- **A second trial** with the SSE2 `nearbyint` and `ldexp`: the wheels again all green, and the
  same suites but for one oracle case in MinGW-w64's Python 3.12 job whose run drew a Windows uuid
  (4 characters) equal to a value of the reference, which the comparison then masked on one side
  only: the test harness's, not the toolchain's (the same case passed in the other three Windows
  jobs, and in that job's re-run: 27,433 passed).

## Consequences

- MinGW-w64 builds use the filters' FMA dispatch, so the fused chains' NaN rules hold on Windows
  as on Linux (`tests/platform/test_nan_rules.py` runs them).
- The Windows wheels carry GCC's and MinGW-w64's runtime libraries, and the toolchain's notices
  of them as it ships them: the build copies them from the toolchain's `licenses/` folder into
  `petakit5d/_ext/licenses/` (the add-on's: `petakit5d_fftw/licenses/`), checking their SHA-256
  where the toolchain is MinGW-builds' GCC 12.2.0, and the release checks require them in a wheel
  whose binaries GNU ld linked, at those exact paths, which they leave out of their search for
  folders and site names (THIRD_PARTY_NOTICES.md).
- The MinGW-w64 wheels are larger: 5.9 MB against 4.2 for the abi3 wheel (the modules 10.3 MB
  against 5.6 unpacked: each links libstdc++ and its vector and dispatch code).
- A process can mix a MinGW-built port with other packages' `libgomp-1.dll` of another GCC: the
  first loaded serves all, and libgomp is backward compatible (the modules call the entry points of
  GCC 4.9 and later).
- New C math functions in the kernels: `cmake/Pk5dExtension.cmake`'s `PK5D_MINGW_LIBM` lists the
  wrapped ones; a function outside it, other than the exact ones (floor, ceil, trunc, rint, frexp,
  modf, fmod, fabs, copysign, sqrt), needs adding there.
