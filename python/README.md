# petakit5d

The Python port of PetaKit5D: petabyte-scale 5D image processing for light-sheet microscopy
(deskew and rotation, deconvolution, stitching, conversions and the cluster job framework),
with the same results and the same files on disk as the MATLAB version. The hot paths are
C++ and CUDA kernels behind Python bindings.

Status: in development. The MATLAB code in the parent folder is the reference until the port
reaches parity.

## Install from source

Python 3.11 or newer and a C/C++17 compiler with OpenMP:

| System | Compiler | Before `pip install .` |
|---|---|---|
| Linux | GCC 9 or newer (or Clang with libomp) | nothing |
| macOS (Apple silicon) | Apple Clang (Xcode command line tools) | `brew install libomp` |
| Windows | MinGW-builds' GCC 12.2.0 with the UCRT (`bash tools/get_mingw.sh` fetches it), else MSVC from Visual Studio 2022 ("Desktop development with C++"): `docs/platforms.md` | nothing |

```
pip install ./python            # or: pip install "./python[test]" for the test tools
```

The build uses scikit-build-core, CMake and nanobind; pip fetches them. A GPU is optional:
the package imports and runs every CPU path without CuPy and without a GPU.

On Linux the modules share the OpenMP runtime of the cpp-zarr and cpp-tiff wheels, so that a
process has one pool of threads: the build links the libgomp that the cpp-zarr wheel bundles
(pip fetches cpp-zarr for the build; with `--no-build-isolation` or CMake directly, install it
first) and puts a copy next to the modules. `-DPK5D_OPENMP_SHARE=OFF` (with pip:
`-Ccmake.define.PK5D_OPENMP_SHARE=OFF`) links the toolchain's instead
(`docs/decisions/0008-openmp-runtime.md`).

### CPU FFTs

The CPU FFTs (`petakit5d.matlabcompat.fft`, `petakit5d.fft`) run MATLAB's FFT driver on one of
two backends:

- `permissive`, the default: pocketfft (BSD-3), built in. Results differ from MATLAB's in the
  last bits (single precision: relative differences up to about 5e-7); the measurements and
  tolerances are in `docs/decisions/0006-fft-backends.md`.
- `exact`: FFTW 3.3.8 built as MATLAB builds it, an optional add-on of two shared libraries
  (GPL-2.0-or-later, so not part of this package). With it the FFTs give MATLAB's bits on Linux
  x86-64. Build and install it with `pip install ./python/cpp/fft/exact`, then choose it with
  `PETAKIT5D_FFT=exact` (or `auto`: exact when installed) or `petakit5d.fft.set_backend("exact")`;
  `PETAKIT5D_FFTW` can point at its libraries instead (a folder, or the two paths). Every call
  may also name its backend (`backend="exact"`), so both can be used in one process.

On Windows, MSVC compiles pocketfft without its vector code: the results are the same, the FFTs
take 1.2 to 1.3 times as long (`docs/decisions/0006-fft-backends.md`). clang-cl (Visual Studio's
"C++ Clang tools for Windows") compiles the vector code; the build has not been tried with it yet.

## Development

Build the extension modules in a build tree outside the repository (never inside it: MATLAB's
`setup.m` puts every folder of the repository on its path) and point the package at it:

```
pip install -r python/requirements/linux-x86_64-py312.txt   # or: pip install "./python[dev]"
cmake -S python -B ../pk5d-build -G Ninja
cmake --build ../pk5d-build
export PYTHONPATH=$PWD/python/src PETAKIT5D_EXT_DIR=$PWD/../pk5d-build/petakit5d/_ext
python -m pytest python/tests
```

`cmake --build` is incremental, and picks up new kernel folders by itself. Lint and format
with `ruff check python` and `ruff format python`.

### Layout

```
src/petakit5d/   the package: core, params, matlabcompat, io, kernels, gpu, jobs, tools,
                 flatfield, deskew_rotate, decon, stitch, pipeline, cli, testing
cpp/             C++ kernels, one folder per kernel group; common/ shared headers;
                 cuda/ GPU kernel sources (compiled at run time)
tests/           pytest, one folder per area
docs/            conventions, decisions, parity exceptions, mapping to the MATLAB names
requirements/    pinned environments (*.in are the inputs of the lock files)
```

### Kernel groups

Each folder `cpp/<group>/` with a `CMakeLists.txt` becomes the module `petakit5d._ext.<group>`:

```cmake
pk5d_add_extension(<group> <sources>...
                   [AVX2] [AVX2_SOURCES <sources>...] [NO_OPENMP]
                   [LIBRARIES <libs>...] [INCLUDE_DIRS <dirs>...] [DEFINITIONS <defs>...])
```

Every kernel is built with `-O3`, OpenMP and `-ffp-contract=off` (MSVC: `/O2 /fp:precise`),
never with `-march=native` or fast-math, so that results don't depend on the machine;
`AVX2` adds `-mavx2` on x86-64 only. Kernels take `(y, x, z)` Fortran-ordered arrays, release
the interpreter lock and take their thread count as an argument. See
`docs/decisions/0001-build.md`.

### Tests

Markers: `gpu` (skipped without a usable GPU), `cluster` (skipped unless
`PETAKIT5D_TEST_SLURM_CONFIG` is set), `demodata` (skipped unless `PETAKIT5D_DEMO_DATA` is set),
`slow`, `oracle`. Site-specific values come only from environment variables:
`PETAKIT5D_ORACLE_DIR`, `PETAKIT5D_GOLDEN_DIR`, `PETAKIT5D_DEMO_DATA`,
`PETAKIT5D_TEST_SLURM_CONFIG`, `PETAKIT5D_TEST_SLURM_SCRATCH`, `PETAKIT5D_TEST_FFTW`.

## Citation

petakit5d is the Python port of PetaKit5D. Please cite the PetaKit5D paper if you find the
software useful for your research:

`Ruan, X., Mueller, M., Liu, G., Görlitz, F., Fu, T., Milkie, D.E., Lillvis, J.L., Kuhn, A., Chong, J.G., Hong, J.L., Herr, C.Y.A., Hercule, W., Nienhaus, M., Killilea, A.N., Betzig, E. and Upadhyayula, S. Image processing tools for petabyte-scale light sheet microscopy data. Nature Methods 21, 2342–2352 (2024). https://doi.org/10.1038/s41592-024-02475-4`

## License

BSD 3-Clause from the first release on; see `LICENSE`. (The history's first commits, of
2026-11-11, carried a copy of the MATLAB repository's license file until it was replaced the same
day; no release was made under it.) The notices of the third-party work the package derives from or
bundles are in `THIRD_PARTY_NOTICES.md` (included with the package).

The CPU FFTs use pocketfft (BSD 3-Clause, `cpp/fft/third_party/pocketfft/LICENSE.md`, included
with the package). The package contains no GPL code. FFTW, which gives MATLAB's exact FFT bits,
is GPL-2.0-or-later and comes only as the separate, optional add-on `petakit5d-fftw`
(`cpp/fft/exact`, a package of its own, left out of this package's source distribution). A
program that is distributed together with the add-on is covered by the GPL. The GPU path uses
NVIDIA's cuFFT, which is not GPL.
