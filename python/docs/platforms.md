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
- **Two NaN operands** (NaNs of other bits meeting in a sum: MATLAB's NaN and one with a
  payload, such as NumPy's `np.nan`). Both processors give the first operand's NaN (ARM a
  signalling second one's), but which operand of `a + b` or `a * b` is first is the compiler's
  choice: GCC on x86-64 and Clang on ARM swap them at some sums, so the NaN follows the build,
  not the code. MATLAB gives the first NaN operand as written (measured: the integral image's
  mex and imwarp, 99 inputs where NaNs meet at each sum, single and double). Where the port
  computes MATLAB's NaN there, a value that came out NaN is computed again with each
  operation's NaN made so, `pk5d::x86_nans(r, a, b)` (`cpp/common/pk5d/fpexact.hpp`), as the
  resize kernels' `nan_sum` does; values without NaN keep their arithmetic:
  - the integral image on every system (an input with NaN is computed again): 21 of the 21
    inputs as MATLAB (x86-64's build matched 18 before; ARM's matched them, but swapped sums of
    other inputs);
  - imwarp on ARM (each voxel that comes out NaN): 78 of 78. On x86-64 imwarp computes no voxel
    again, 29 of the 78 (32 before the cubic code of 2026-11-19: the compiler's choices follow
    the code); the recomputation there cost 5 to 10% of cubic's speed when it was measured
    (Known problems);
  - the linear filters on ARM (each output value that comes out NaN), with
    `pk5d::x86_fma_nans(r, a, b, c)` for the fused taps. MATLAB's filters take a tap's NaN
    before the running sum's: a fused tap's factors, a plain tap's product, a chain of a padded
    column before the chains above it (measured: convn, imfilter along each axis with each
    padding and on its sparse path, imgaussfilt3, single and double; NaNs of three kinds and
    infinities meeting). x86's fused multiply-adds give that too, in every form (132, 213, 231),
    and keep the addend's NaN at 0 * Inf + NaN; ARM's take the addend's first and give their own
    default NaN at 0 * Inf + NaN (measured: the cluster's Xeon Gold 6146, the M2 Pro). ARM's
    filters now give MATLAB's bits in these outputs (before: 12 to 591 NaNs differed per output)
    but where x86-64's differ too: 12 NaNs of the zero-padded double column pass, whose last
    rows' tap order is not MATLAB's (parity row N5). On x86-64 the filters keep their code: the
    fused chains and the vector loops' plain sums give MATLAB's bits, the one-value path's plain
    sums the compiler's choice (GCC: 13 of 159 NaNs of imfilter's sparse path along x differ
    from MATLAB);
  - the deskew kernels in every build without the mex's AVX2 code (ARM, and `deskew_scalar` on
    x86-64), each float value that comes out NaN: the first NaN operand of the mex's formula as
    written. Against the three mex files (NaNs of three kinds and infinities meeting): the same
    bits in the interpolation to a step size and in every row of the warp's and the fused
    kernel's lines but the mex's scalar tail (the rows past floor(len / 8) * 8), where the mex
    binary takes other operands (76 and 54 of about 12,000 NaNs of the warp, 535 and 554 of
    about 41,500 of the fused kernel, forward and reverse). The AVX2 build keeps the mex's code:
    the compiler's choice, in every row (341 and 424 NaNs of the warp, 552 and 544 of the
    interpolation, 3,531 and 3,408 of the fused kernel differ from the mex).
  - convn and conv2 of a kernel with Inf or NaN (and of an image with them where a kernel
    larger than the image swaps the roles) on every system: every operation's NaN written out
    (`taps_nd_inside`), a product's the image value's before the tap's, a sum's the product's
    (a fused tap's factors') before the running sum's. MATLAB gives these in every probe output
    (NaNs of three kinds and infinities meeting, plain and fused sums, single and double), but
    where its small arrays' sums differ from the port's (Known problems).

  Tests: tests/platform/test_nan_rules.py (a reference of the integral image's passes in
  `core.dtypes`' arithmetic; imwarp's and the deskew kernels' NaNs side by side, on the builds
  that compute the rule; references of the filters' tap programs and tap correlations, the
  fused ones on x86-64 too), `specials.meeting_nans` and `specials.conv_overlap`.
  **The pinned cases keep NaNs of other bits apart** (no output sees two, except where the
  kernel fixes their order itself: the ordered sums, reduce, the integral image, convn and conv2
  of kernels with Inf or NaN), since a sum's NaN is otherwise the compiler's choice, which
  another build can change (the stable-ABI build did for imwarp cubic): kernels with a footprint
  take the payload NaN on an input of its own (`payload_only`) or out of the other NaNs' reach.

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
| the deskew/rotate warp and the skewed space interpolation (C++ `deskew`, both layouts: without AVX2 a float value that comes out NaN computed again with `pk5d::x86_nans`, `deskew/nans_x86.hpp`) | Inf * 0, Inf - Inf, and NaNs that meet | float data with Inf (Inf * 0 at whole coordinates and zero interpolation weights) or NaN | x86's NaN (`specials.deskew_infinities`, tests/platform/test_nan_rules.py); where NaNs of other bits meet, without AVX2 (ARM, `deskew_scalar`) the mex's formula as written, the mex's bits but in its scalar tails (measured), and in the AVX2 build the compiler's choice |
| `matlabcompat.resize` (imresize, imresize3: C++ `resize`, a NaN sum computed again by `nan_sum`) | Inf * 0, Inf - Inf | float data with Inf | x86's NaN (measured on ARM with many infinities, and `specials.warp_and_resize`) |
| `matlabcompat.filters` (convn, conv2, imfilter, imgaussfilt: the multiply-add chains of C++ `filters`; on ARM an output that comes out NaN computed again with `pk5d::x86_fma_nans` and `pk5d::x86_nans`) | Inf * 0, Inf - Inf, and NaNs that meet | float data with Inf and a zero tap, or both infinities in a window, or NaN | x86's NaN (`specials.filters_infinities`, tests/platform/test_nan_rules.py); where NaNs of other bits meet, MATLAB's (a tap's NaN before the running sum's, measured) on ARM, and on x86-64 too but in the one-value plain sums (the compiler's choice); with MSVC's C runtime's `fma` where NaNs meet not measured |
| convn and conv2 of a kernel with Inf or NaN (C++ `filters`' `taps_nd_inside`: MATLAB's sums cover the overlap of the two arrays only, so the taps that would read the zero padding are left out; each operation's NaN written out on every system) | Inf * 0 with the padding (none: left out), Inf * 0 and Inf - Inf in the overlap, and NaNs that meet | a kernel (or, for a kernel larger than the image, an image) with Inf or NaN | MATLAB's values and NaN (measured): a product's NaN the image value's before the tap's, then the product's before the running sum's, x86's default NaN for an invalid operation (`specials.conv_overlap`; oracle cases `conv2.nonfinite.*`, `convn.nonfinite.*`); but where MATLAB's small-array sums differ (Known problems) |
| `matlabcompat.warp.imwarp`, linear and cubic (C++ `warp`: on ARM a voxel that comes out NaN computed again with `pk5d::x86_nans`) | Inf * 0, Inf - Inf, and NaNs that meet | float data with Inf (linear: at a whole coordinate, the second sample's weight 0) or NaN; cubic: double data whose sums pass the largest double (MATLAB's sums take twice Keys' weights and the 1/8 only at the end, so a sum overflows where Keys' would pass a half, a quarter or an eighth of it, along x, y and z: Inf, and NaN where infinities of both signs meet) | x86's NaN (`specials.warp_cubic`, `specials.warp_and_resize`, `specials.infinities`; oracle cases `cubic_extreme_*`); where NaNs of other bits meet, MATLAB's on ARM and the compiler's on x86-64 (tests/platform/test_nan_rules.py, on ARM) |
| integral image (C++ `stitch/integral_image.cpp`: `PK5D_X86_NAN`; an input with NaN computed again with `pk5d::x86_nans`) | Inf - Inf, and NaNs that meet | single tiles with both infinities, or NaN | x86's NaN, the first NaN operand as written (`specials.integral_image`, `specials.meeting_nans`, `specials.infinities`) |
| the resize kernels' fdlibm pow and cubic kernel (C++ `resize`: `PK5D_X86_NAN`) | 0 / 0, Inf - Inf, Inf * 0 | pow: a negative base with a fractional exponent, (+-1)^(+-Inf); the kernel: \|x\|^3 overflows. Never: the cubic weights take \|x\|^3 of distances of at most 2 | x86's NaN (`specials.resize_pow`, `specials.infinities`) |
| the stitching distance maps (`stitch.distance`, the products with the box weights) | Inf * 0 | a negative weight degree (the weights' floor is Inf in single) | x86's NaN (tests/platform/test_nan_rules.py; found on ARM by tests/stitch/test_distance.py `slab.wd-1`) |
| the stitching registration's shift and peak (`stitch.registration`) | 0 / 0 | the slabs that see an axis weigh 0 | x86's NaN (tests/platform/test_nan_rules.py; found on ARM by tests/stitch/test_registration.py `weights.nan.3`, `.6`) |
| decon's ratio `x ./ max(v, eps)`, FFTs of data with Inf | Inf / Inf, Inf * 0, Inf - Inf | the data hold Inf | FFT outputs: recorded (`specials.fft_and_decon`), see the FFT section |

  Real microscope data (uint16) hold no Inf or NaN, so outside the feather blending of
  MATLAB's code these NaNs come from special values in float inputs. In C++ the step is
  `PK5D_X86_NAN(r, a, b)` (`cpp/common/pk5d/fpexact.hpp`), which expands to `r` except on ARM,
  so the x86-64 builds compile the same code (checked by disassembly). The recomputations of the
  filters (ARM only) and of the deskew kernels (builds without AVX2) are left out of the x86-64
  builds the same way: their modules disassemble as before. No hardware mode gives a negative
  default NaN: AArch64's FPCR.DN gives the positive one and drops payloads. The
  `specials.*_infinities` cases hold many infinities (`infinities` in
  tests/platform/_kernel_cases.py, a tenth of the values): with one infinity of each sign, as in
  the other `specials.*` cases, an interpolation or filter makes no invalid operation unless a
  zero weight or the other infinity happens to be in reach.
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
| Linux x86-64: the cluster (Xeon Gold 6146, AVX-512, glibc 2.28, GCC 13.2), Python 3.12 | 55 | the pins (made there) |
| Linux x86-64: runner-like container (Ubuntu 24.04, glibc 2.39, GCC 13.3, AVX2), Python 3.12 and 3.14 | 55 | the same bits |
| macOS ARM (M2 Pro, macOS 13.7, Apple Clang 14, libomp 23), Python 3.12 and 3.14 | 45 | the same bits (before the filters' and the deskew kernels' NaN fix the 2 `specials.*_infinities` cases differed in NaN signs only); `specials.conv_overlap`, `specials.warp_cubic` since its sums are MATLAB's, and point detection's cases (`detection.*`, `neighbors.rangesearch`, `signal.*`), not measured yet |
| Windows x86-64 (GitHub's windows-2025, MSVC), Python 3.12 and 3.14 | 45 | the same bits (the CI, measured 2026-11-19: see "How Windows is tested"); `specials.conv_overlap`, `specials.warp_cubic` since its sums are MATLAB's, and point detection's cases (`detection.*`, `neighbors.rangesearch`, `signal.*`), not measured yet |

The cases whose digests are recorded, not pinned, and how much they differ (largest difference
in ulps of the largest |value| of the output; FFT outputs in single precision included):

| Case | What it holds | Linux glibc 2.28 against 2.39 | Linux against macOS ARM |
|---|---|---|---|
| `fft.matlab_fft` | fft, ifft, fftn, ifftn | the same bits | single: the same bits; double: 35 to 85% of the values, by at most 4 ulp (4.7e-16 of the largest) |
| `fft.half_spectra` | half spectra of FFT convolutions, OTFs | the same bits | single: the same bits; double: 73 to 87%, at most 5 ulp (6.0e-16) |
| `filters.edgetaper` | edgetaper | the same bits | single: the same bits; double: 43%, at most 2 ulp (2.5e-16) |
| `decon.otf_and_rl` | decon's OTF, RL and OMW steps (single) | the same bits | the same bits |
| `xcorr.ncc` | the registration's NCC | the same bits | single: the same bits; double: 81%, at most 1.5 ulp (2.8e-16) |
| `specials.fft_and_decon` | FFTs and an RL step of data with Inf and NaN | the same bits | the single FFT and the RL step: the same bits; the double fftn: NaN signs (Inf - Inf's are ARM's; the FFT negates the payload NaN and adds it to itself, so its sign follows the build), finite values by at most 1 ulp |
| `detection.ad_stat` | the Anderson-Darling statistic of point detection's fits (`adtest_mex`'s, with the C library's `erfc` and `log`, as MATLAB's mex) | the same bits | not measured yet |
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
    xfails of gaps that 1.8.0 fixes now pass). Since the oracle keys take tiff and zarr inputs by
    their decoded content, output tiffs are compared by voxels and tags and the binding tests
    follow the installed version, all but 8 of them pass with the 1.8.0 wheels on Linux. The 8
    are cases whose MATLAB functions report the input files' sizes in bytes, which 1.8.0 changes
    by a few bytes; their inputs need a writer independent of the library before the port moves
    to 1.8.0.
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
  every pinned case must hold, and does: with MSVC, on Python 3.12 and 3.14, every pinned
  digest equals Linux's (the CI, measured 2026-11-19; before the collision-free cases only
  `specials.warp_cubic` differed there). The exact oracle cases, which compare the port's
  outputs with MATLAB's digests, hold on Windows too (the CI is green); on macOS ARM with the
  1.8.0 libraries they fail only where they fail on Linux with them, apart from the NaN signs
  above. The FFT-based digests can differ by the UCRT's cos and sin; `--save` and `--compare`
  measure them where a run's files can be fetched.

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
- **x86-64: where NaNs of other bits meet, the filters' one-value plain sums and the deskew
  kernels' AVX2 build** give the compiler's choice of operand, not MATLAB's (the counts above):
  the filters along x on imfilter's sparse path (13 of 159 NaNs measured), and every deskew
  kernel (3 to 9% of the NaNs where many meet). ARM computes MATLAB's; x86-64 keeps the code it
  compiled before (its modules disassemble as before), as for imwarp below. Data from MATLAB hold
  MATLAB's NaN only, for which the choice changes no bits.
- **The deskew mex's scalar tails:** where NaNs of other bits meet in the rows past
  floor(len / 8) * 8 of the warp's and the fused kernel's lines, the mex binary takes other
  operands than its source as written (the interpolation mex's tails do not), so the builds
  without AVX2 differ from it there (0.5 to 1.3% of the NaNs measured, all in those rows); their
  other rows give the mex's bits.
- **MATLAB's conv2 and convn of small arrays** take other sums than the port's, on every system
  (measured, black box): below about 24 rows and columns (2-D) or 16 planes (3-D) MATLAB adds
  plain products (each rounded) over the overlap, with the array of fewer elements as the taps,
  in column-major order (3-D: z, then x descending, then y). The port has fused chains over the
  kernel there, and for a kernel larger than the image in every dimension plain sums in
  column-major order in 3-D too. Random data differ in the last bits in up to about half of the
  values, and where NaNs of other bits meet the NaN can differ (a plain sum keeps a product's
  default NaN at 0 * Inf + NaN, a fused tap the addend's). A kernel larger in every dimension
  than an image of 10 or more rows (2-D) or 8 or more planes (3-D), and 9 x 9 kernels, take forms
  not identified yet. The port's own code calls neither function. With Inf or NaN in a kernel,
  which values are Inf, NaN or finite was MATLAB's in every case probed (small and large, 2-D
  and 3-D, swapped; no 9 x 9 kernel): the sums cover the overlap only, as MATLAB's.
- **Windows: the filters' fused taps call the C runtime's `fma`** (MSVC builds have no dispatch).
  Its invalid operations give x86's NaN (`specials.filters_infinities` holds there); its NaN where
  NaNs of other bits meet, and at 0 * Inf + NaN, is not measured (the fused chains' tests skip
  there); glibc's (with the FMA instruction) is the instruction's.
- **x86-64: imwarp's NaN where NaNs of other bits meet** is the compiler's choice of operand,
  not MATLAB's first as written, in 49 of 78 measured inputs (above), and changes with the build
  (the stable-ABI build flipped 2 voxels of a case that let them meet). ARM computes MATLAB's.
  The same recomputation on x86-64 slowed cubic down by 5 to 10% when it was measured (uint16
  too, whose code it does not touch), through GCC's inlining of `Extended`'s end values into the
  voxel loop, which followed any change to `cpp/warp/warp.cpp`. The rows that need end values
  are now fetched outside the loop (`extended_row`), so that cost may be gone; it has not been
  measured again. Data from MATLAB hold MATLAB's NaN only, for which the choice changes no bits.
