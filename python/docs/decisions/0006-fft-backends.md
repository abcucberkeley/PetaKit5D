# 0006: Two CPU FFT backends: pocketfft by default, FFTW as an optional add-on

Status: accepted (2026-11-12): pocketfft is the default CPU FFT backend, FFTW an optional GPL
add-on; the tolerances below are proposed. Measured on Linux x86-64 (24-core nodes, Xeon Gold
6146, glibc 2.28, GCC 13.2) against MATLAB R2024b Update 5.

## Context

The CPU FFTs ran only on FFTW 3.3.8 (0002), which gives MATLAB's bits but is GPL-2.0-or-later,
so every binary of the port was covered by the GPL. The user asked whether a permissive library
is "fast and very close in result". A study measured permissive libraries behind MATLAB's driver
and on their own against MATLAB: speed at 1, 4 and 24 threads, closeness per transform, and
deconvolution and registration end to end. The same change moved the FFT API to the port's
(z, y, x) C arrays (`conventions.md`).

| Library | License | Result |
|---|---|---|
| pocketfft (C++ header, the library NumPy bundles) | BSD-3-Clause | chosen: behind MATLAB's driver, 0.27-0.43 of MATLAB's time at 24 threads (study, six sizes) |
| `numpy.fft` (pocketfft) | BSD-3-Clause | one thread only: `fftn` + `ifftn` of 512 x 512 x 500 single 13.4 s, 2.6x pocketfft behind the driver on one thread |
| `scipy.fft` (SciPy 1.18: ducc0's FFT) | BSD-3-Clause OR GPL-2.0-or-later | as fast; not chosen: its own stages, no MATLAB driver, and a dependency on SciPy's choice of library |
| ducc0 itself | GPL-2.0-or-later | excluded |
| Intel MKL (FFTW interface) | Intel Simplified Software License | not taken: x86 only, a large binary dependency |

## Decision

### Two backends behind one driver

Both backends run MATLAB's driver (`cpp/fft/driver.cpp`, 0002): its stages, its line order, its
groups of four lines and lines alone, its real-to-complex and complex-to-real stages, the
Hermitian fill of `fftn`, the exact symmetry test of `ifft` and `ifftn`, the division after every
inverse stage and the real result when every imaginary part is zero. The backend only computes
the 1D transforms of one line or a group of four, through a table of the FFTW API functions the
driver calls (`cpp/fft/backend.h`).

- **`permissive`** (the default): pocketfft C++ (`cpp/fft/third_party/pocketfft`, commit
  4d865cba), built in, behind an adapter with FFTW's conventions (`cpp/fft/pocketfft_api.cpp`).
  pocketfft computes the lanes of a group as single lines, so the results do not depend on
  MATLAB's FFT mode, the groups, the vector width (SSE2, AVX and scalar builds give the same
  bits) or our thread count.
- **`exact`**: FFTW 3.3.8 with MATLAB's codelet sets (0002), built by `cpp/fft/exact/fftw.cmake`
  as two shared libraries and loaded at run time. With it the results are MATLAB's bit for bit.
  It is a package of its own, `petakit5d-fftw` (`pip install python/cpp/fft/exact`; in a build
  tree `-DPK5D_FFT_EXACT=ON`). No file outside `cpp/fft/exact` includes `fftw3.h`.
- **Choosing:** `petakit5d.fft.backend()`, `set_backend(name)`, `using(name)`; `PETAKIT5D_FFT`
  (`permissive`, `exact`, or `auto`: exact when installed) at import; and `backend=` on every
  call, so both backends can run in one process. The exact libraries come from
  `PETAKIT5D_FFTW` (a folder, or the two paths), the installed `petakit5d_fftw` package, or
  `PETAKIT5D_TEST_FFTW` (MATLAB's own `libmwfftw3f` / `libmwfftw3`).

### One FFT state per process

The backends, the default backend, FFTW's planner lock and the plan cache are one per process.
They live in `petakit5d._ext.fft`, which exports them as the capsule
`petakit5d._ext.fft._shared_state`: a versioned table of the functions of `cpp/fft/backend.h`.
Another extension module that links the driver library attaches to it when it is imported
(`cpp/fft/shared.h`: `attach_to_fft_module()`), so that FFTW is loaded once, every plan is made
under one lock (FFTW's planner is not thread-safe, and the transforms run without the interpreter
lock) and plans are cached once and destroyed when the exact backend is unloaded. Kernels get
MATLAB's plans with `cached_plan(backend, precision, kind, n, lines)`. A module built from other
sources than `_ext.fft` refuses to import. `tests/fft/test_shared_state.py` plans from two modules
at once without the interpreter lock.

### The license boundary

- The package (wheel and source distribution) is BSD-3-Clause and contains no GPL code. It
  ships pocketfft's BSD notice (`license-files`). pocketfft's header holds one piece derived from
  FFTW (odd-sized DCT-IV), which FFTW's authors released under BSD-3-Clause for it (its notice);
  the port uses no DCTs, so that code is not even compiled. The source distribution leaves
  `cpp/fft/exact` out; FFTW's sources are fetched (URL and SHA-256) only when the add-on is
  built. Checked on the built wheel (no FFTW code, names or license text; pocketfft's license
  included) and the source distribution.
- The add-on's wheel bundles FFTW and is GPL-2.0-or-later, with FFTW's license file. A program
  distributed together with it is covered by the GPL.
- Tests in CI: the extension contains no FFTW code, codelet names or version strings, and on
  Linux no FFTW library is mapped unless the add-on is loaded
  (`tests/fft/test_backends.py::test_the_default_build_contains_no_fftw`); no file outside
  `cpp/fft/exact` includes `fftw3.h`.

### The API on (z, y, x) arrays

`petakit5d.matlabcompat.fft.fft(x, n, axis)`, `ifft(x, n, axis, symmetric)`, `fftn(x, shape)` and
`ifftn(x, shape, symmetric)` take and return the port's C arrays. MATLAB's dimension 1 (y) is
axis -2, dimension 2 (x) axis -1, dimension k >= 3 axis -k; the default axis is MATLAB's first
dimension longer than 1, in MATLAB's order; `shape` is the Python shape of the padded or truncated
array (MATLAB's `siz`, mapped; one length is a row, as a 1-D array is). The result is
`from_matlab` of MATLAB's result of `to_matlab(x)`.

Inside, the driver keeps MATLAB's stages and numbers the lines as MATLAB does (its linear order
of the other dimensions); its threads visit them in the C array's memory order. Lines along y
are strided: a group of four is four neighbouring columns, gathered as four contiguous lanes.
Lines along x are contiguous, lines along z strided by a plane.

### Half spectra for the kernels (`petakit5d.fft.HalfSpectrum`)

For FFT convolutions (deconvolution, K05/K06): `forward(x)`, `multiply(s, h, conjugate)`,
`inverse(s)`, `otf(psf)` (as `decon_psf2otf`), `convolve(x, h, conjugate)`.

- exact: the elements of MATLAB's `fftn(x)` with index <= n/2 along MATLAB's last transformed
  dimension (z of a volume, x of an image) and MATLAB's `ifftn` of the conjugate-symmetric
  product, with MATLAB's complex product (no fused multiply-add): `convolve` equals MATLAB's
  `real(ifftn(H .* fftn(x)))` bit for bit. The transforms cost about what `fftn` and `ifftn` do
  (MATLAB's stages before the last run on the whole array); the spectra, OTFs and products are
  half the size.
- permissive: pocketfft's real-to-complex and complex-to-real transforms on the native layout,
  halved along x (a third of the exact backend's time; speed below).
- The product is a pass of its own over the half spectrum (it is not fused into the transforms'
  loops).

A spectrum's layout is its backend's (`spectrum_shape`); OTFs and spectra of one object belong
together.

### Windows

MSVC has no GCC or Clang vector extensions, which pocketfft's vector code needs: it compiles
pocketfft's scalar code (the CMake configure step says so). The results are the same bits. The
FFTs take 1.18 to 1.33 times as long: pocketfft built without vectors against the vector build,
on Linux, `fftn` + `ifftn` of 512 x 512 x 500, 2304 x 2304 and 3082 x 11530 x 4 single at 1, 4
and 24 threads (3 interleaved repeats). clang-cl (Visual Studio's "C++ Clang tools for Windows")
compiles the vector code; the build has not been tried with it (no Windows machine in this task):
from a Developer prompt, `pip install ./python -C cmake.define.CMAKE_C_COMPILER=clang-cl
-C cmake.define.CMAKE_CXX_COMPILER=clang-cl`; the other kernels then need clang-cl's OpenMP.

## Tolerances (proposed)

For results of the permissive backend against MATLAB's:

- single precision: max |difference| <= 1e-5 of max |MATLAB's|, and relative L2 <= 1e-6;
  double precision: 1e-13 and 1e-14;
- a result MATLAB returns as real (every imaginary part exactly zero) may come back complex with
  tiny imaginary parts, and the reverse: compare the values;
- deconvolution (up to about 30 iterations): relative L2 <= 1e-4 of the result, uint16 outputs
  off by at most 1 in at most 0.1% of the voxels;
- registration: the same integer peak where the normalized cross-correlation is well behaved,
  peak values within 5e-3.

The exact backend has no tolerance. Both backends' bits follow the C library: FFTW and pocketfft
make their twiddle factors with its cosine and sine, in double (GCC merges the two calls into one
call of `sincos`, in MATLAB's FFTW as in ours). glibc's `sincos` is not always its `sin` and
`cos`: on glibc 2.28 and 2.35 they differ in the last bit at 944 of FFTW's arguments for lengths
up to 2048, while glibc 2.39 (Ubuntu 24.04) gives `sin`'s and `cos`'s values there. MATLAB's own
FFTW gives other bits on such a system too. Tests compare with MATLAB's FFT bits only where the
function FFTW calls gives the values of the system they were recorded on
(`petakit5d.testing.fft_reference.digests_hold_here`, `twiddle_fingerprint`).

## Results

### Exact backend: MATLAB's bits through the (z, y, x) API

F08's oracle cases, given in MATLAB's form and mapped as the oracle's roles map them
(`from_matlab`, the dimension's axis, the size's shape, `to_matlab`):

| Cases | MATLAB's FFT modes | Our threads | Identical |
|---|---|---|---|
| 1D lengths 2..1024, every kind, both precisions | parallel, single-thread | 1 and 24 (tests: 1, 3, 8) | 8,768 of 8,768 |
| MATLAB's rules (real results, symmetry test, padding, `siz`, `'symmetric'`, signed zeros, ...) | parallel, single-thread | 1 and 24 | 779 of 779 |
| 3D sizes of the demos, the decon tests and a real decon batch | parallel, single-thread | 24; 1 and 3 | 130 of 130 |
| Argument rules (errors, empty and scalar inputs, `n` and `siz`) | parallel | any | 40 of 40 |

The half spectra (`HalfSpectrum`, exact) equal MATLAB's `fftn` and `real(ifftn(H .* fftn(x)))`
bit for bit in the tests, and 20 Richardson-Lucy iterations on the demo data with half spectra
give the same bits as with MATLAB's full spectra (confocal 231 x 551 x 551, widefield
501 x 571 x 501).

### Permissive backend against MATLAB

The 9,677 pack cases (MATLAB's parallel mode; the 1D and rules packs give the same numbers against
MATLAB's single-thread mode): the 10 where MATLAB raises an error raise one here too; the others
per output (MATLAB's class and realness):

| Output | Cases | Identical | Largest max abs difference / max abs | Largest relative L2 |
|---|---|---|---|---|
| `fft`, single, complex | 2,289 | 17 | 5.1e-7 | 3.7e-7 |
| `fft`, single, real | 35 | 32 | 1.7e-7 | 1.2e-7 |
| `ifft`, single, complex | 1,214 | 16 | 5.4e-7 | 3.7e-7 |
| `ifft`, single, real | 1,085 | 36 | 4.8e-7 | 3.4e-7 |
| `fftn`, single, complex | 72 | 5 | 5.4e-7 | 5.3e-7 |
| `fftn`, single, real | 15 | 15 | 0 | 0 |
| `ifftn`, single, complex | 35 | 2 | 5.4e-7 | 5.3e-7 |
| `ifftn`, single, real | 88 | 28 | 5.5e-7 | 4.8e-7 |
| double, every function | 4,834 | 147 | 1.2e-15 | 1.1e-15 |

Every case is within the tolerances. In 3 cases (`fft` of conjugate-even input) MATLAB's result
and ours differ in realness: one has imaginary parts of about 1e-16 of the largest value where
the other has exact zeros.

Deconvolution, 20 iterations of `decon_lucy_function`'s CPU path in single precision, against the
exact backend (MATLAB's bits), result as MATLAB's `uint16`:

| Data | FFTs | Relative L2 | Max abs difference / max abs | uint16 voxels off by 1 | off by more |
|---|---|---|---|---|---|
| confocal, 231 x 551 x 551 | full spectra | 2.8e-6 | 5.3e-6 | 1,500 of 70.1 M (0.002%) | 0 |
| | half spectra | 3.1e-6 | 4.5e-6 | 1,589 (0.002%) | 0 |
| widefield, 501 x 571 x 501 | full spectra | 1.3e-5 | 1.6e-5 | 86,646 of 143.3 M (0.060%) | 0 |
| | half spectra | 1.4e-5 | 1.7e-5 | 91,239 (0.064%) | 0 |

Registration (the study; pocketfft behind F08's driver, whose results equal this backend's: the
1D and rules packs give the same digests): MATLAB's normalized cross-correlation of the demo
data's phase, confocal and 2P tiles and of synthetic pairs. On the 284 well-conditioned pairs
(124 of real phase tiles, 160 synthetic) every integer peak is MATLAB's and the peak values are
within 1.4e-3. The 104 2P pairs are degenerate in single precision: MATLAB's own peak values
reach 14.5 (a correlation coefficient is at most 1), its double-precision peaks differ from its
single-precision ones in 103 of 104, and every backend finds other peaks there.

### Speed

`fftn` and `ifftn` (of the conjugate-symmetric result) of single-precision arrays, per call and
the pair `ifftn(fftn(x))` (end to end), on (z, y, x) arrays; MATLAB R2024b in its parallel mode
and each backend in processes of their own, interleaved, warm (one warm-up per size), on a 24-core
node (two Xeon Gold 6146). "exact" is FFTW 3.3.8 behind the driver on (z, y, x) arrays, "before"
F08's driver on MATLAB's layout; the half spectra are `HalfSpectrum.forward` and `inverse`.

The gate holds: at 512 x 512 x 500 the permissive backend takes 0.30, 0.63 and 0.61 of MATLAB's
time for the pair at 24, 4 and 1 threads (per call: `fftn` 0.47 to 0.75, `ifftn` 0.19 to 0.51).
At the other sizes its pair takes 0.29 to 0.81 of MATLAB's time, its `fftn` alone up to 1.03
(1800 x 512 x 500 at 4 threads). The layout change costs the exact backend up to 1.18 of F08's
time where its first stage, real-to-complex along y, runs on strided lines of long rows
(3082 x 11530 x 4; 1800 x 512 x 500 at 24 threads 1.11) and nothing at 512 x 512 x 500 and
2304 x 2304 (0.81 to 0.97). With pocketfft the half spectra, the kernels' path, take 0.12 to 0.37
of MATLAB's full transforms.

(medians of 5 repeats)

Gate, 512 x 512 x 500 single, seconds (fftn / ifftn / pair) and the permissive backend's time over MATLAB's:

| Threads | MATLAB | permissive | permissive / MATLAB |
|---|---|---|---|
| 24 | 0.411 / 0.619 / 1.042 | 0.193 / 0.119 / 0.312 | 0.47 / 0.19 / 0.30 |
| 4 | 1.187 / 1.180 / 2.383 | 0.887 / 0.604 / 1.496 | 0.75 / 0.51 / 0.63 |
| 1 | 3.845 / 3.705 / 7.550 | 2.719 / 1.847 / 4.572 | 0.71 / 0.50 / 0.61 |

Every size, time over MATLAB's (fftn / ifftn / pair):

| Size (MATLAB's) | Threads | permissive | exact | half spectra, permissive | half spectra, exact |
|---|---|---|---|---|---|
| 512 x 512 x 500 | 24 | 0.47 / 0.19 / 0.30 | 0.46 / 0.19 / 0.31 | 0.15 / 0.10 / 0.12 | 0.38 / 0.17 / 0.28 |
| 512 x 512 x 500 | 4 | 0.75 / 0.51 / 0.63 | 0.66 / 0.45 / 0.56 | 0.25 / 0.26 / 0.25 | 0.67 / 0.39 / 0.53 |
| 512 x 512 x 500 | 1 | 0.71 / 0.50 / 0.61 | 0.61 / 0.44 / 0.52 | 0.27 / 0.28 / 0.27 | 0.60 / 0.37 / 0.49 |
| 1800 x 512 x 500 | 24 | 0.54 / 0.26 / 0.38 | 0.56 / 0.26 / 0.38 | 0.14 / 0.10 / 0.12 | 0.56 / 0.21 / 0.36 |
| 1800 x 512 x 500 | 4 | 1.03 / 0.61 / 0.81 | 0.91 / 0.56 / 0.73 | 0.29 / 0.25 / 0.27 | 0.93 / 0.51 / 0.72 |
| 1800 x 512 x 500 | 1 | 0.90 / 0.64 / 0.76 | 0.81 / 0.58 / 0.69 | 0.35 / 0.32 / 0.33 | 0.83 / 0.48 / 0.65 |
| 2304 x 2304 | 24 | 0.46 / 0.20 / 0.29 | 0.42 / 0.18 / 0.27 | 0.18 / 0.12 / 0.16 | 0.45 / 0.15 / 0.25 |
| 2304 x 2304 | 4 | 0.97 / 0.62 / 0.78 | 0.89 / 0.53 / 0.69 | 0.26 / 0.24 / 0.26 | 0.94 / 0.45 / 0.68 |
| 2304 x 2304 | 1 | 1.01 / 0.66 / 0.81 | 0.88 / 0.57 / 0.71 | 0.30 / 0.28 / 0.29 | 0.88 / 0.44 / 0.64 |
| 3082 x 11530 x 4 | 24 | 0.32 / 0.26 / 0.29 | 0.39 / 0.33 / 0.36 | 0.13 / 0.11 / 0.12 | 0.40 / 0.32 / 0.36 |
| 3082 x 11530 x 4 | 4 | 0.80 / 0.61 / 0.70 | 1.01 / 0.79 / 0.88 | 0.31 / 0.25 / 0.28 | 1.01 / 0.77 / 0.88 |
| 3082 x 11530 x 4 | 1 | 0.94 / 0.70 / 0.80 | 1.23 / 0.93 / 1.06 | 0.42 / 0.32 / 0.37 | 1.24 / 0.91 / 1.05 |

The layout change: the exact backend on (z, y, x) arrays over F08's driver on MATLAB's layout (fftn / ifftn / pair):

| Size (MATLAB's) | 24 threads | 4 threads | 1 thread |
|---|---|---|---|
| 512 x 512 x 500 | 0.87 / 1.05 / 0.97 | 0.91 / 1.00 / 0.97 | 0.94 / 1.01 / 0.97 |
| 1800 x 512 x 500 | 1.12 / 1.09 / 1.11 | 0.99 / 1.03 / 1.02 | 0.94 / 1.04 / 0.99 |
| 2304 x 2304 | 0.93 / 0.78 / 0.88 | 0.82 / 0.85 / 0.81 | 1.00 / 0.97 / 0.97 |
| 3082 x 11530 x 4 | 1.07 / 1.11 / 1.07 | 1.07 / 1.14 / 1.10 | 1.16 / 1.20 / 1.18 |

## Not checked, and open

- macOS ARM and Windows: the permissive backend runs there in CI (tolerance tests); the exact
  digests are skipped there, and the add-on was not built there.
- clang-cl on Windows.
- The exact backend's layout cost where its first stage runs on strided lines of long rows (up
  to 1.18 of F08's time, 3082 x 11530 x 4): that stage on a scratch in MATLAB's layout (one
  transpose of the real input) would remove it, at the price of a second complex array.
