# Installing petakit5d

The package is installed from source until wheels are published. It needs Python 3.11 or newer,
a C/C++17 compiler with OpenMP and about a minute to build its C++ modules.

## The base package

| System | Compiler | Before `pip install` |
|---|---|---|
| Linux x86-64 | GCC 9 or newer (or Clang with libomp) | nothing |
| macOS (Apple silicon) | Apple Clang (Xcode command line tools) | `brew install libomp` |
| Windows x86-64 | MSVC from Visual Studio 2022 ("Desktop development with C++") | nothing |

From a clone of the repository, in a fresh virtual environment:

```
python -m venv /opt/envs/petakit5d
source /opt/envs/petakit5d/bin/activate        # Windows: \opt\envs\petakit5d\Scripts\activate
pip install ./python
```

This installs NumPy, SciPy, scikit-image, h5py, cpp-tiff and cpp-zarr, and the two commands
`pk5d` and `pk5d-master` ([cli.md](../cli.md)). Check the install with

```
pk5d list                     # the ported functions, by Python and MATLAB name
python -c "import petakit5d; print(petakit5d.__version__)"
```

On Linux the C++ modules share the OpenMP runtime of the cpp-zarr wheel, so a process has one
pool of threads ([decisions/0008](../decisions/0008-openmp-runtime.md)). The build flags and the
toolchains of each system are in [platforms.md](../platforms.md) and
[decisions/0001](../decisions/0001-build.md).

**Known problem on macOS:** the 1.7.3 wheels of cpp-zarr and cpp-tiff crash when both are used
in one process (the second one imported, on its first write). Linux and Windows are not
affected. The fix belongs to those packages ([platforms.md](../platforms.md), "Known problems").

## Optional extras

Each extra adds the dependencies of a group of functions. Without it the rest of the package
works, and a function that needs the extra says which one to install.

| Extra | Install | For | Licenses |
|---|---|---|---|
| `figures` | `pip install "./python[figures]"` | the figures of the decon wrappers (back projector and OTF masks), the MIP masks, the PSF and FSC plots, tile positions (`petakit5d.tools.figures`): matplotlib | matplotlib (PSF-based, permissive) |
| `n5` | `pip install "./python[n5]"` | N5 with blosc, MATLAB's N5 compression, and zstd before Python 3.14 (`petakit5d.io.n5`): numcodecs | MIT |
| `formats` | `pip install "./python[formats]"` | ND2 and CZI files (`petakit5d.io.bioformats`): nd2, czifile, imagecodecs. Bio-Formats itself is not ported | BSD-3 (their dependencies BSD-3 or MIT) |
| `basic` | `pip install "./python[basic]"`, then `pip install --no-deps basicpy==2.0.0` | the phase-image flat field (`petakit5d.flatfield.basic`): BaSiCPy and PyTorch | MIT, BSD |
| `gpu` | `pip install "./python[gpu]"` | the GPU paths (deconvolution, its OTFs and back projectors): CuPy, CUDA 12.9's runtime and NVRTC, cuFFT 11.0.8 ([GPU](#gpu) below) | CuPy MIT, cuda-pathfinder Apache-2.0; NVIDIA's CUDA license for the CUDA libraries |
| `test` | `pip install "./python[test]"` | pytest, to run the test suite | MIT |
| `dev` | `pip install "./python[dev]"` | pytest and ruff, for development | MIT |

Extras combine: `pip install "./python[figures,n5,formats]"`.

**BaSiCPy:** BaSiCPy 2.0.0 declares `scipy<1.13`, which this package's `scipy>=1.13` rules out,
although it uses only `scipy.fft`. The `basic` extra therefore installs its dependencies, and
BaSiCPy goes in after them with `--no-deps`. The phase flat field matches MATLAB's BaSiC to about
a percent, not bit for bit (`parity/flatfield.md`).

## GPU

A GPU is optional: every CPU path imports and runs without CuPy and without a GPU. The GPU paths
(deconvolution, its OTFs and back projectors) need an NVIDIA GPU, a driver for CUDA 12 and the
`gpu` extra's packages, all from PyPI (no CUDA toolkit is needed):

```
pip install "./python[gpu]"
```

The extra pins cuFFT to 11.0.8.103 and keeps CuPy, CUDA's runtime and NVRTC and cuda-pathfinder
within the releases the GPU results were measured with (CuPy 14, CUDA 12.9). The measured
versions themselves:

```
pip install cupy-cuda12x==14.2.0 nvidia-cuda-runtime-cu12==12.9.79 \
    nvidia-cuda-nvrtc-cu12==12.9.86 nvidia-cufft-cu12==11.0.8.103 cuda-pathfinder==1.8.6
```

cuFFT 11.0.8 is the version MATLAB R2024b bundles: with it the GPU FFTs give MATLAB's bits.
CUDA 12.9's NVRTC compiles for every GPU from Volta to Blackwell. The reasons and measurements
are in [decisions/0003](../decisions/0003-gpu.md); `python -m petakit5d.gpu.selftest` checks a
GPU environment.

## The exact FFT add-on (FFTW, GPL)

The CPU FFTs run MATLAB's FFT driver on one of two backends:

- **`permissive`, the default:** pocketfft (BSD-3), built in. Results differ from MATLAB's in the
  last bits: single-precision transforms by up to about 5e-7 relative, deconvolved uint16 volumes
  by 1 count in a few thousandths of a percent of the voxels. The tolerances are record 0006
  ([decisions/0006](../decisions/0006-fft-backends.md)); the exceptions that follow from them are
  listed in [parity/index.md](../parity/index.md). One of them: many OMW iterations on data where
  they don't settle, such as the widefield demo's 20, amplify the last bits in a few bright spots,
  as they amplify the difference between MATLAB's own CPU and GPU results
  ([parity/decon.md](../parity/decon.md)).
- **`exact`:** FFTW 3.3.8 built as MATLAB builds it, the separate package `petakit5d-fftw`. With
  it the FFTs, and so deconvolution, registration, spectra and FSC, give MATLAB's bits on Linux
  x86-64.

Install and choose the add-on:

```
pip install ./python/cpp/fft/exact          # fetches FFTW's release by URL and checksum and builds it
export PETAKIT5D_FFT=exact                  # or auto: exact when installed, else permissive
```

`petakit5d.fft.set_backend("exact")` chooses it in a running process, and functions that take
FFTs accept `backend="exact"` per call. `PETAKIT5D_FFTW` can point at the add-on's two libraries
instead (a folder, or the two paths).

**License note.** FFTW is GPL-2.0-or-later, and so is the add-on. petakit5d itself is
BSD-3-Clause and contains no GPL code; the add-on is built and distributed as a package of its
own and is left out of petakit5d's source distribution. Using it on your own machines puts no
obligation on you; a program or environment that you distribute together with the add-on is
covered by the GPL. The GPU path uses NVIDIA's cuFFT, which is not GPL.

## Development installs

For work on the package, build the C++ modules in a build tree outside the repository and point
the package at it; the steps are in the package's `README.md` ("Development"). Never build
inside the repository: MATLAB's `setup.m` puts every folder of the repository on its path.
