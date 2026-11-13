# Platforms: Linux, macOS ARM and Windows

The port gives the same bits on every supported system wherever its own code computes them.
This page lists the toolchains, what each system does differently, and how the port keeps the
results the same. The build flags are in `docs/decisions/0001-build.md`.

## Toolchains

| | Linux x86-64 | macOS ARM (Apple silicon) | Windows x86-64 |
|---|---|---|---|
| compiler | GCC (reference: 13) | Apple Clang | MSVC (Visual Studio 2022 or 2026) |
| OpenMP | libgomp | libomp (Homebrew's, or `OpenMP_ROOT`) | `/openmp:llvm`, runtime `libomp140.x86_64.dll` next to the modules |
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
  Tests: tests/platform/test_nan_rules.py. The places, and what is still open, are in the
  table below.

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
| not covered by tests yet: decon's ratio `x ./ max(v, eps)`, the interpolations (warp, deskew), FFTs of data with Inf | Inf / Inf, Inf * 0 | the data hold Inf | to check (E04) |

  Real microscope data (uint16) hold no Inf or NaN, so outside the feather blending of
  MATLAB's code these NaNs come from special values in float inputs. In C++ the step is
  `PK5D_X86_NAN(r, a, b)` (`cpp/common/pk5d/fpexact.hpp`), which expands to `r` except on ARM,
  so the x86-64 builds compile the same code (checked by disassembly). No hardware mode gives a
  negative default NaN: AArch64's FPCR.DN gives the positive one and drops payloads.
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
- **OpenMP version.** MSVC's `/openmp:llvm` has OpenMP 3.1's features but reports `_OPENMP`
  200203 (2.0).
- **The C library's elementary functions** (exp, log, sin, pow, hypot) differ in the last bit
  between glibc versions, Apple's libm and the UCRT. Where MATLAB's bits matter the port uses
  its own fdlibm-based functions (`petakit5d.matlabcompat.elementary`), not libm or NumPy's
  transcendental ufuncs.
- **NumPy's SIMD dispatch and OpenBLAS's kernels** depend on the CPU (AVX-512, AVX2, NEON). The
  port's sums that feed outputs run in a fixed order in its own kernels.
- **File systems.** Windows and macOS are case-insensitive by default; the few tests that need
  names differing only in case, or symbolic links, skip there.

## Known problems

- **macOS: cpp-zarr and cpp-tiff in one process.** The 1.7.3 wheels both export the same weak
  C++ template code; macOS merges such symbols across loaded modules, so whichever of the two
  is imported second crashes (segmentation fault) on its first write. Linux and Windows are
  unaffected. The fix is in the packages (export only `PyInit_*`, hide the bundled C++
  runtime), not in the port.
