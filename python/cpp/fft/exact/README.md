# petakit5d-fftw: the exact CPU FFT backend (optional, GPL)

FFTW 3.3.8, built as MATLAB R2024b builds its `libmwfftw3f` / `libmwfftw3` (the codelet sets
`fftw-3.3.8-sse2-avx`: SSE2 and AVX, no FMA), as two shared libraries. With them, petakit5d's FFTs
give MATLAB's bits on Linux x86-64 (`petakit5d.fft.set_backend("exact")` or `PETAKIT5D_FFT=exact`):
MATLAB's own driver (`python/cpp/fft/driver.cpp`, BSD-3) runs on FFTW's plans through the backend
table. Without them petakit5d uses pocketfft (BSD-3), whose results differ from MATLAB's in the last
bits (`python/docs/decisions/0006-fft-backends.md`).

**License: GPL-2.0-or-later** (FFTW's). This package is built and distributed apart from petakit5d,
which is BSD-3 and contains no GPL code. A program that loads these libraries is covered by FFTW's
license when it is distributed with them.

## Build

```
pip install python/cpp/fft/exact            # the add-on's wheel (petakit5d_fftw/ with the libraries)
```

or in the port's own build tree, for development and the oracle tests:

```
cmake -S python -B <build> -DPK5D_FFT_EXACT=ON && cmake --build <build>
export PETAKIT5D_FFTW=<build>/fftw_exact
```

`fftw.cmake` fetches the release tarball by URL and SHA-256 (`PK5D_FFTW_VERSION`, 3.3.8 by default;
3.3.10 gives the same plans and bits and is not faster), compiles FFTW's own source lists with the
SIMD flags only on the SIMD codelets, `-ffp-contract=off` and no fast-math, and leaves out FFTW's
threads and Fortran wrappers. The F08 report (`python/docs/decisions/0002-fft.md`) has the
measurements: 603 codelets per precision, those of MATLAB's libraries; 65,728 plan-level cases and
9,717 oracle cases bit-identical with MATLAB.

FFTW makes its twiddle factors with the C library's `cos` and `sin`, as MATLAB's FFTW does, so the
bits are MATLAB's where these functions give the values of the platform the references were made on
(glibc 2.28, x86-64); the tests check that before comparing digests.
