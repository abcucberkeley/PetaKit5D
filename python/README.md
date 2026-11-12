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
| Windows | MSVC from Visual Studio 2022 ("Desktop development with C++") | nothing |

```
pip install ./python            # or: pip install "./python[test]" for the test tools
```

The build uses scikit-build-core, CMake and nanobind; pip fetches them. A GPU is optional:
the package imports and runs every CPU path without CuPy and without a GPU.

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

## License

GPL-3.0-or-later; see `LICENSE`.
