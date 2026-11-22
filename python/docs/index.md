# petakit5d documentation

petakit5d is the Python port of PetaKit5D: petabyte-scale 5D image processing for light-sheet
microscopy, with MATLAB PetaKit5D's results and files.

## For users

- [guide/index.md](guide/index.md): the user guide, and the three ways to call a function
  (Python API, MATLAB names, command line);
  - [guide/install.md](guide/install.md): the base package, the extras (`figures`, `n5`,
    `formats`, `basic`), the GPU packages, the exact FFT add-on and its GPL note;
  - [guide/examples.md](guide/examples.md): the MATLAB demos in Python, on the demo datasets;
  - [guide/cluster.md](guide/cluster.md): Slurm jobs and cluster configuration files.
- [migration.md](migration.md): moving from MATLAB PetaKit5D, PetaKit5D-GUI or PyPetaKit5D:
  unchanged configs plus `PythonLaunchStr`, the GUI and PyPetaKit5D switch, the files that differ.
- [MATLAB_MAPPING.md](MATLAB_MAPPING.md): every MATLAB function and parameter with its Python
  name and default, and the functions not ported with their reasons (generated from the
  registry).
- [cli.md](cli.md): the `pk5d` and `pk5d-master` command lines and the `run_mccMaster.sh` shim.
- [parity/index.md](parity/index.md): every known difference from MATLAB's results and files,
  with its status.

## How the port works

- [conventions.md](conventions.md): array order, boxes, MATLAB's classes and arithmetic, numbers
  as text, paths, threads.
- [jobs.md](jobs.md): the job framework (local, pool and cluster modes, task records, runners).
- [platforms.md](platforms.md): Linux, macOS and Windows: toolchains and the same bits on each.
- [detection_scope.md](detection_scope.md): Phase 3's scope: point detection, tracking and lifetime
  analysis as MATLAB runs them, traced from the demos (stages, inputs and outputs, what is reached).
- [release.md](release.md): the release plan and checklist: versions, the wheels, the licenses of
  what they bundle or require, and the move into a repository of its own.
- [testing.md](testing.md): the oracle harness that compares the port with MATLAB, and the test
  suite; [test_map/README.md](test_map/README.md): what became of each MATLAB test.
- [hardening.md](hardening.md): what the suite misses: the mutation campaign, coverage per
  area, the isolation of every task function and the soak of the job framework.
- [benchmarks.md](benchmarks.md): the benchmark harness and the speed gates;
  [benchmarks_e2e.md](benchmarks_e2e.md): whole workflows against MATLAB.

## Decision records

| Record | Subject |
|---|---|
| [0001](decisions/0001-build.md) | build system, compilers and flags |
| [0002](decisions/0002-fft.md) | FFTs identical to MATLAB's, without MATLAB |
| [0003](decisions/0003-gpu.md) | GPU support: CuPy, run-time compiled kernels, cuFFT plans, device choice |
| [0004](decisions/0004-io-bindings.md) | tiff and zarr through cpp-tiff and cpp-zarr; native zarr |
| [0005](decisions/0005-sums-and-math-library.md) | MATLAB's sums, statistics and elementary functions, exactly |
| [0006](decisions/0006-fft-backends.md) | two CPU FFT backends: pocketfft by default, FFTW as an optional add-on; their tolerances |
| [0007](decisions/0007-openmp-spin.md) | OpenMP's spin after a parallel region |
| [0008](decisions/0008-openmp-runtime.md) | one OpenMP runtime with cpp-zarr and cpp-tiff |
| [0009](decisions/0009-n5-layout.md) | standard N5 |
| [0010](decisions/0010-module-names.md) | a function's module names its MATLAB folder or kind (the renames before the first release) |
| [0011](decisions/0011-windows-toolchain.md) | Windows: MinGW-w64's GCC for the release wheels and source builds, MSVC as the fallback (the one runtime with cpp-zarr and cpp-tiff, the x87 precision, the toolchain's libraries) |

## Per-area differences

[decon.md](parity/decon.md), [demos.md](parity/demos.md), [flatfield.md](parity/flatfield.md),
[matlabcompat.md](parity/matlabcompat.md), [phase2.md](parity/phase2.md) (the demos' Phase 2
steps), [realdata.md](parity/realdata.md), [stitch.md](parity/stitch.md),
[tools.md](parity/tools.md); the index of all of them is [parity/index.md](parity/index.md).
