# Platforms: Linux, macOS ARM and Windows

The port gives the same bits on every supported system wherever its own code computes them.
This page lists the supported systems, the toolchains, what each system does differently and by
how much, and how the port keeps the results the same. The build flags are in
`docs/decisions/0001-build.md`.

## Supported systems and Python versions

| | Linux x86-64 | macOS ARM (Apple silicon) | Windows x86-64 |
|---|---|---|---|
| system | glibc 2.28 or newer (the manylinux_2_28 wheels of cpp-tiff and cpp-zarr) | macOS 13 or newer (their macosx_13_0_arm64 wheels) | 64-bit Windows 10 or newer (their win_amd64 wheels) |
| CI (GitHub Actions) | ubuntu-24.04: gates every push | macos-15: reports, does not gate (see "macOS" below) | windows-2025: gates every push (since 2026-11-14) |
| also tested on | the reference cluster (glibc 2.28, AVX-512); a runner-like container (glibc 2.39, AVX2) | an M2 Pro Mac, macOS 13.7 | |

- **Python:** 3.11 or newer (`requires-python`). CI tests 3.12, with the pinned environment on
  Linux (`requirements/linux-x86_64-py312.txt`), and 3.14 with the newest packages. cpp-tiff and
  cpp-zarr have wheels for 3.9 to 3.14 and for the free-threaded 3.14t; the port is not tested
  on free-threaded Python.
- **Not supported:** macOS on Intel processors, Linux on ARM, Windows on ARM and 32-bit systems:
  there are no cpp-tiff and cpp-zarr wheels for them.
- **The GPU paths** (CuPy, CUDA) are tested on Linux only (`docs/decisions/0003-gpu.md`); without
  a usable GPU the CPU paths run.

## Toolchains

| | Linux x86-64 | macOS ARM (Apple silicon) | Windows x86-64 |
|---|---|---|---|
| compiler | GCC (reference: 13) | Apple Clang | MSVC (Visual Studio 2022 or 2026) |
| OpenMP | libgomp: the cpp-zarr wheel's, shared with cpp-zarr and cpp-tiff (`docs/decisions/0008-openmp-runtime.md`) | libomp: Homebrew's (`brew install libomp`), or another found through `OpenMP_ROOT` (a folder with `include/omp.h` and `lib/libomp.dylib`, such as conda-forge's llvm-openmp) | `/openmp:llvm`, runtime `libomp140.x86_64.dll` next to the modules |
| no contraction | `-ffp-contract=off` | `-ffp-contract=off` | `/fp:precise` |
| vector code | GCC vector extensions, AVX2 dispatch | Clang vector extensions (NEON) | the same operations lane by lane |

`petakit5d._ext.core.build_info()` names the compiler, the flags and the OpenMP runtime, and
`fma_contracted()` checks at run time that no multiply-add was fused (tests/test_smoke.py).

**Windows: MSVC, not clang-cl.** CPython, NumPy and SciPy are built with MSVC, and `pip install`
from source uses it with a stock Visual Studio. Every kernel has a path for compilers without
GCC's vector extensions and target attributes: the run-time instruction-set dispatch of the MIP
and flat-field kernels runs the baseline loops, binary morphology and the 3 x 3 x 3 median
filter do their compare-exchanges lane by lane, and pocketfft runs one line at a time. These
are the same operations in the same order (selections, or element-wise IEEE operations), so
the bits are the same; only the speed differs. clang-cl would take the vector paths too (it
defines `__clang__`); it needs Visual Studio's Clang component, an LLVM OpenMP runtime to ship,
and compiler-rt's builtins for `__builtin_cpu_supports`, and is left for the wheels (E07) if
Windows speed matters.

## What differs by system, and what the port does about it

- **Fused multiply-adds.** Clang on ARM contracts `a * b + c` by default; the build turns it
  off everywhere, and refuses flags that turn it back on (`-ffast-math`, `-ffp-contract=fast`,
  `/fp:fast`, ...).
- **The NaN of an invalid operation** (Inf - Inf, 0 * Inf, 0 / 0, Inf / Inf). x86's default
  NaN has the sign bit set (0xFFF8000000000000, single 0xFFC00000); it is MATLAB's NaN. ARM's
  has the sign bit clear (0x7FF8000000000000). A NaN operand keeps its payload, quieted, on
  both. Where such a NaN can reach an output, the port writes x86's: `core.dtypes`' arithmetic,
  and `core.dtypes.x86_nan(r, *operands)` after other NumPy operations (nothing to do on x86).
  Tests: tests/platform/test_nan_rules.py, and the special-value cases of
  tests/platform/test_kernel_hashes.py (`specials.*`). The places, and what is still open, are
  in the table below.

| Where | Invalid operation | MATLAB's result has it when | Status |
|---|---|---|---|
| `core.dtypes` plus, minus, times, rdivide, ldivide | all | the operands hold Inf or 0 | x86's NaN |
| sums (`kernels.ops`, C++ `ops_sum`; reduce's own NaN rule) | Inf + -Inf | the data hold both infinities | x86's NaN |
| decon's acceleration sums (C++ `decon/accel.cpp`) | Inf + -Inf | the iterates hold both infinities | x86's NaN |
| `params.literal` (MATLAB's text values) | `Inf-Inf`, `Inf*0` | the text says so | x86's NaN |
| `matlabcompat.stats.prctile` | Inf - Inf in the interpolation | -Inf and Inf are neighbours | x86's NaN |
| `matlabcompat.neighbors` sph2cart | Inf * 0 | an Inf radius at a zero sine or cosine | x86's NaN |
| `stitch.geometry` batch grid | 0 / 0 | a stitched size with a zero | x86's NaN |
| `kernels.stitch` feather blending, MATLAB's code | 0 / 0, Inf * 0, Inf - Inf | a voxel without weight; Inf data or weights | x86's NaN |
| tile blending (C++ `stitch/tile_blending.cpp`: 'none', and the last step `nv .* (1 - z) + nvf .* z`) | Inf * 0 | the tiles hold Inf | x86's NaN |
| feather blending (C++ `stitch/feather_blending.cpp`) | Inf * 0, Inf / Inf | Inf weights | x86's NaN |
| flat-field correction (C++ `flatfield/kernels.inc`, and its NumPy reference) | Inf - Inf | Inf in the data and the background | x86's NaN |
| scale and offset (C++ `core/ops/ops_convert.cpp`) | Inf * 0, Inf - Inf | Inf data with a zero scale or an infinite offset | x86's NaN |
| the deskew/rotate warp and the skewed space interpolation (C++ `deskew`) | Inf * 0, Inf - Inf | float data with Inf | x86's NaN (measured on ARM: `specials.deskew`) |
| `matlabcompat.warp.imwarp`, linear; `matlabcompat.resize` (imresize, imresize3); `matlabcompat.filters` (convn, imfilter, imgaussfilt) | Inf * 0, Inf - Inf | float data with Inf | x86's NaN (measured on ARM: `specials.warp_and_resize`, `specials.filters`) |
| `matlabcompat.warp.imwarp`, cubic | Inf * 0, Inf - Inf | float data with Inf | **ARM's NaN** (`specials.warp_cubic`, expected to fail on ARM) |
| integral image (C++ `stitch/integral_image.cpp`) | Inf - Inf | single tiles with both infinities | **ARM's NaN** (`specials.integral_image`) |
| the resize kernels' fdlibm pow (C++ `resize`) | 0 / 0 | a negative base with a fractional exponent: never (the cubic weights take \|x\|^3) | **ARM's NaN** (`specials.resize_pow`) |
| the stitching distance maps (`stitch.distance`, the products with the box weights) | Inf * 0 | a negative weight degree (the weights' floor is Inf in single) | x86's NaN (tests/platform/test_nan_rules.py; found on ARM by tests/stitch/test_distance.py `slab.wd-1`) |
| the stitching registration's shift and peak (`stitch.registration`) | 0 / 0 | the slabs that see an axis weigh 0 | x86's NaN (tests/platform/test_nan_rules.py; found on ARM by tests/stitch/test_registration.py `weights.nan.3`, `.6`) |
| decon's ratio `x ./ max(v, eps)`, FFTs of data with Inf | Inf / Inf, Inf * 0, Inf - Inf | the data hold Inf | FFT outputs: recorded (`specials.fft_and_decon`), see the FFT section |

  Real microscope data (uint16) hold no Inf or NaN, so outside the feather blending of
  MATLAB's code these NaNs come from special values in float inputs. In C++ the step is
  `PK5D_X86_NAN(r, a, b)` (`cpp/common/pk5d/fpexact.hpp`), which expands to `r` except on ARM,
  so the x86-64 builds compile the same code (checked by disassembly). No hardware mode gives a
  negative default NaN: AArch64's FPCR.DN gives the positive one and drops payloads. The rows
  marked **ARM's NaN** differ from MATLAB only in the sign bit of those NaNs, and only on ARM.
- **64-bit integers to single.** AArch64 compilers convert vectors of int64 to float through
  double, which rounds twice above 2^53 (2^53 + 2^29 + 1 became 2^53, not 2^53 + 2^30). The
  port's cast (`kernels.ops.convert(x, 'single', mode='cast')`) rounds once on every processor
  (`pk5d::int64_to_float`); MATLAB's `single()` rounds through double on purpose and still
  does. NumPy's `astype(np.float32)` of int64 rounds twice on ARM too: don't use it where MATLAB
  rounds once. Tests: tests/platform/test_conversions.py.
- **MSVC's selects.** MSVC (19.5x, also with `/fp:strict`) compiles `a >= b ? a : b` and
  `a <= b ? a : b` to `maxss`/`minss`, which return b on a tie: `max(-0, 0)` gave +0 where
  MATLAB keeps -0 (decon's acceleration step, the flat-field correction, the clamps of
  `kernels.ops`). `PK5D_GE_SELECT` and `PK5D_LE_SELECT` make the select on the bits under MSVC
  and are the plain expressions elsewhere.
- **Line ends.** git on Windows checks text files out with CRLF by default (core.autocrlf), which
  changed the bytes of the MATLAB helper and harness files that oracle keys hash.
  `python/.gitattributes` (`* -text`) checks every file out as committed. Text files that tests
  or the port write must be written with `newline="\n"` (Python's text mode writes CRLF on
  Windows); the oracle cases' input generators do. cpp-zarr writes `.zarray` and the other zarr
  metadata in text mode, so with CRLF on Windows: MATLAB reads them the same, and the oracle
  keys and comparisons take zarr metadata with LF line ends (`testing.cases.tree_hash`,
  `testing.compare`).
- **MATLAB's parameter files on Windows.** `writeJsonFile` writes the JSON with
  `fprintf(fid, s)`, which turns the `\\` of a Windows path into `\`, so `parameters.json`
  is not valid JSON where a parameter holds a Windows path; the port writes the same bytes.
- **Tests that need POSIX tools** (the fake Slurm, the fake rsync, bash tasks) skip on Windows.
- **Task runners** (`jobs.single`'s pools) run on Linux only; elsewhere every single job runs
  its own task.
- **OpenMP version.** MSVC's `/openmp:llvm` has OpenMP 3.1's features but reports `_OPENMP`
  200203 (2.0).
- **The C library's elementary functions** (exp, log, sin, pow, hypot) differ in the last bit
  between glibc versions, Apple's libm and the UCRT. Where MATLAB's bits matter the port uses
  its own fdlibm-based functions (`petakit5d.matlabcompat.elementary`), not libm or NumPy's
  transcendental ufuncs. The FFTs' twiddle factors are the exception (next section).
- **NumPy's SIMD dispatch and OpenBLAS's kernels** depend on the CPU (AVX-512, AVX2, NEON). The
  port's sums that feed outputs run in a fixed order in its own kernels. Dense LAPACK is used
  in one place, the null spaces of `matlabcompat.lsq` (next section).
- **File systems.** Windows and macOS are case-insensitive by default; the few tests that need
  names differing only in case, or symbolic links, skip there.

## Measured: the kernels' bits on each system

tests/platform/test_kernel_hashes.py runs every kernel on fixed inputs and compares the digest
of its outputs with the one recorded on Linux x86-64 (`kernel_hashes.json`). The inputs come
from exact arithmetic (the counter-based generator of `testing.synth`, literals) and hold the
special values (MATLAB's NaN, a NaN with a payload, both infinities, -0, subnormals, the largest
values) where the kernels take float data; volumes are never cubic and long enough for the
vector loops, their tails and the fixed-order sums' blocks. The cases reach every binding of
the extension modules and every function of `matlabcompat` (two tests check that by recording
the calls), so a new kernel needs a case. A digest that differs fails its case by name;
`python tests/platform/_kernel_cases.py --save out.npz` on both systems and `--compare a.npz
b.npz` measure the difference.

| System | Pinned cases (non-FFT kernels) | Result |
|---|---|---|
| Linux x86-64: the cluster (Xeon Gold 6146, AVX-512, glibc 2.28, GCC 13.2), Python 3.12 | 41 | the pins (made there) |
| Linux x86-64: runner-like container (Ubuntu 24.04, glibc 2.39, GCC 13.3, AVX2), Python 3.12 and 3.14 | 41 | the same bits |
| macOS ARM (M2 Pro, macOS 13.7, Apple Clang 14, libomp 23), Python 3.12 and 3.14 | 41 | 38 the same bits; the 3 `specials.*` cases of the **ARM's NaN** rows differ in NaN signs only, as expected (strict xfails) |
| Windows x86-64 (GitHub's windows-2025, MSVC), Python 3.12 and 3.14 | 41 | from the CI run of this test: see "How Windows is tested" |

The cases whose digests are recorded, not pinned, and how much they differ (largest difference
in ulps of the largest |value| of the output; FFT outputs in single precision included):

| Case | What it holds | Linux glibc 2.28 against 2.39 | Linux against macOS ARM |
|---|---|---|---|
| `fft.matlab_fft` | fft, ifft, fftn, ifftn | the same bits | single: the same bits; double: 35 to 85% of the values, by at most 4 ulp (4.7e-16 of the largest) |
| `fft.half_spectra` | half spectra of FFT convolutions, OTFs | the same bits | single: the same bits; double: 73 to 87%, at most 5 ulp (6.0e-16) |
| `filters.edgetaper` | edgetaper | the same bits | single: the same bits; double: 43%, at most 2 ulp (2.5e-16) |
| `decon.otf_and_rl` | decon's OTF, RL and OMW steps (single) | the same bits | the same bits |
| `xcorr.ncc` | the registration's NCC | the same bits | single: the same bits; double: 81%, at most 1.5 ulp (2.8e-16) |
| `specials.fft_and_decon` | FFTs and an RL step of data with Inf and NaN | the same bits | the single FFT and the RL step: the same bits; the double fftn: its NaNs from Inf - Inf have ARM's sign, its finite values differ by at most 0.5 ulp |
| `lsq.lsqlin` | lsqlin and the assignment's exact optimum (dense SVD) | 4 of 6 values, by at most 1.5 ulp (1.8e-16); the bounded assignment 2 of 6, 0.5 ulp | 3 of 6 values, 1 ulp (1.2e-16); the bounded assignment: the same bits |

- **Why FFTs differ:** pocketfft (the default backend) computes its twiddle factors with the C
  library's `cos` and `sin` (in double, for single and double data), from an angle computed in
  `long double` (80 bits on x86-64, 64 on ARM and with MSVC). Single-precision transforms round
  those factors to float, so their last-bit differences almost never show; double-precision ones
  keep them. glibc 2.28 and 2.39 gave the same factors for these lengths; Apple's libm did not.
  FFTW (the exact add-on) takes its factors from the C library as well
  (`testing.fft_reference`: glibc 2.28 and 2.35 differ from 2.39 at 944 arguments up to 2048).
  The decon outputs, single precision throughout, have the same bits on the systems measured.
- **Why lsqlin differs:** the exact optimum of `matlabcompat.lsq` (with equality rows, or when
  the least squares break a bound) takes null spaces by NumPy's dense SVD, so its last bits
  follow the LAPACK and BLAS of NumPy's build and their kernels for the CPU: OpenBLAS's Haswell
  and SkylakeX kernels already differ, and macOS 13 takes NumPy's OpenBLAS wheel (NEON kernels;
  macOS 14 and newer take its Accelerate wheel). The stitching shifts are rounded
  from these values, and the ones within 1e-12 of a rounding boundary are flagged
  (`near_round_boundary`). The assignment whose bounds hold (sparse LU and sparse products only)
  is pinned (`lsq.assignment`).
- **The tolerance** for FFT-based outputs on macOS ARM and Windows against Linux is that of the
  permissive backend against MATLAB (`docs/decisions/0006-fft-backends.md`); the measured
  differences above are far inside it.

## macOS: cpp-tiff and cpp-zarr 1.7.3 and 1.8.0

- **With the 1.7.3 wheels** (the published ones), whichever of cpp-tiff and cpp-zarr is imported
  second crashes (segmentation fault) on its first write: both wheels export the same weak C++
  template code, and macOS merges such symbols across loaded modules (Linux loads modules with
  RTLD_LOCAL; Windows does not merge). CI's macOS jobs therefore run the tests in pytest-xdist
  workers (`-n 3 --max-worker-restart=2000`), so that a crash ends only its worker: about 1,100
  tests crash, each restart costs about 8.5 s, and the step takes from 30 minutes to more than
  its 110-minute limit, so its results are incomplete and vary from run to run. Those jobs
  report without gating.
- **With the 1.8.0 builds** (built 2026-11-18 and 19, not yet published: the libraries export
  only their module init and hide their C++ runtime), the crash is gone: both import orders
  write and read back in fresh processes, and the whole CPU suite runs in one process, without
  workers, in about 13 minutes on an M2 Pro (Python 3.12 and 3.14). Its results on test's head
  of 2026-11-19 (25,363 passed, 1,225 failed, 4,988 skipped):
  - 1,221 failures are those of the same suite on Linux with the 1.8.0 wheels, so they come with
    1.8.0, not with macOS: cpp-tiff 1.8.0 writes other bytes than 1.7.3 (the tests' input files,
    so 1,016 oracle cases get new keys; the port's tiff outputs, so about 180 cases see MATLAB's
    voxels in other bytes), and 12 binding tests expect 1.7.3's messages and gaps (10 strict
    xfails of gaps that 1.8.0 fixes now pass). Moving the port to 1.8.0 updates them.
  - 4 failed on macOS only, and pass or skip now: the NaN signs of the stitching distance maps
    and the registration (3 cases; their table rows above), and tests/jobs/test_single.py's
    `test_call_ends_its_pools_and_folders`, which expects a task runner's folder (Linux only)
    and skips elsewhere.
- **Once 1.8.0 is published and the port requires it,** the macOS jobs run the tests like the
  others (no workers) and can gate.

## How Windows is tested

- **On GitHub's runners only:** windows-2025, Python 3.12 and 3.14, the port built from source
  with MSVC as `pip install` builds it (`.[test,figures,n5]`), and the whole CPU suite
  (`pytest -v --timeout=900`), on every push. The job gates since 2026-11-14 (after two green
  runs in a row). No local Windows machine runs the port's tests.
- **What the CI shows** besides the tests: the build information (compiler, flags, OpenMP
  runtime), the CPU and NumPy's SIMD dispatch and OpenBLAS kernels, the OpenMP runtimes found,
  and the "Oracle key probe (Windows)" step (tests/platform/key_probe.py), which prints what
  goes into the keys of oracle cases whose keys differed on Windows (each input file's size and
  hash), to compare with Linux.
- **The kernels' bits on Windows** are those of tests/platform/test_kernel_hashes.py there:
  every pinned case must hold. The exact oracle cases, which compare the port's outputs with
  MATLAB's digests, hold on Windows too (the CI is green); on macOS ARM with the 1.8.0
  libraries they fail only where they fail on Linux with them, apart from the NaN signs above.
  The FFT-based digests can differ by the UCRT's cos and sin; `--save` and `--compare` measure
  them where a run's files can be fetched.

## Writing tests that pass on every system

- Make inputs from exact arithmetic: `testing.synth`'s `CounterRNG` (Philox, the same bits
  everywhere; its `normal` uses correctly rounded operations) and literals. No NumPy
  transcendental ufuncs or libm for inputs, and no `re + 1j * im` for complex data with
  infinities (`0 * Inf` is an invalid operation): set the real and imaginary parts.
- Write text files with `newline="\n"`, build paths with `pathlib`, compare normalized paths,
  take temporary folders from pytest (`tmp_path`), and run no shell commands; a test that needs
  POSIX tools skips on Windows with that reason.
- Where a result can only match on the reference platform, skip with a precise reason where a
  fingerprint shows the platform differs (`testing.fft_reference.digests_hold_here()` for
  MATLAB's FFT digests).
- A difference that a system still has is a strict xfail with its reason on that system
  (`ARM_NAN` in tests/platform/_kernel_cases.py), so that the fix shows.

## Known problems

- **macOS: cpp-zarr and cpp-tiff 1.7.3 in one process** crash on the first write of the second
  one imported (above). The fix is in the packages (1.8.0), not in the port.
- **ARM: three kernels still give ARM's NaN** for an invalid operation where MATLAB gives x86's:
  the **ARM's NaN** rows of the table above. They need special values in float data (the resize
  pow none that the port gives it), and differ from MATLAB only in the NaN's sign.
