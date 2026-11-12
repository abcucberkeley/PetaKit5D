# CUDA kernel sources

The GPU kernels are plain `.cu` files, one folder per kernel group (`decon/`, ...), compiled
at run time by NVRTC through `petakit5d.gpu.kernels`. Nothing here is compiled when the package
is built: no CUDA toolkit is needed to build or install, and one install serves every NVIDIA
GPU, also those newer than the NVRTC in use (see `docs/decisions/0003-gpu.md`).

## Layout

```
cpp/cuda/<group>/<name>.cu      kernels
cpp/cuda/<group>/<name>.cuh     headers of the group (#include "name.cuh")
```

The build copies the `.cu`, `.cuh`, `.h` and `.hpp` files to `petakit5d/_ext/cuda/<group>/`
(next to the extension modules); `kernels.cuda_source(group, filename)` looks there first, then
in this source tree (`petakit5d._ext.cuda_source_dirs()` gives the order).

## Using a kernel

```python
from petakit5d.gpu import device, kernels

device.select_device()  # the job's GPU, made current
mod = kernels.load_module(
    kernels.cuda_source("decon", "rl_rfft.cu"), name_expressions=("init_kernel<InputU16>",)
)
mod["multiply_kernel"](grid, 256, (s, h, np.int64(nh), np.bool_(False)))
mod["init_kernel<InputU16>"](grid, 256, (j, y, inp, np.int64(n)))
```

Arguments are CuPy arrays and NumPy scalars of the parameter's exact type (`np.int64`,
`np.float32`, `np.uint16`, ...). Kernels run on the current CuPy stream unless one is given
(`stream=`).

## Rules for the sources

- **Device code only.** NVRTC has no host headers: no `cufft.h`, `cuda_runtime.h`, `mex.h`,
  `stdint.h`, `float.h` or `math.h`. Use `float2` for complex values (cuFFT's `cufftComplex`
  has that layout), `unsigned short` / `long long` for fixed-size integers, and literal
  constants for `DBL_EPSILON` and the like. CUDA's math functions (`floorf`, `__fmul_rn`, ...)
  are built in.
- **Names.** Kernels that Python calls by name are `extern "C"`; templates are compiled per
  instantiation listed in `name_expressions` and called by that expression.
- **Floating point.** Sources are compiled with `--std=c++17 --fmad=false`: no fused
  multiply-adds and no flush of subnormals, as the C++ kernels (`-ffp-contract=off`). Where a
  MATLAB mex used the rounding intrinsics (`__fmul_rn`, `__fadd_rn`, `__fdiv_rn`), keep them:
  they pin the result whatever the options.
- **Ported from a mex**, a kernel keeps the mex's loop and arithmetic so the two compare bit
  for bit; say in a comment which mex function it is.
- **Indices** are `long long` (volumes over 2^31 voxels), grid-stride loops as in the mex.

## Testing

GPU tests are marked `gpu` and run in short GPU jobs only. To check the path a GPU newer than
the installed NVRTC takes (PTX compiled by the driver) on any GPU, set
`PETAKIT5D_GPU_TARGET=ptx`, or `PETAKIT5D_GPU_TARGET=<arch>` (e.g. `90`) to act as if NVRTC
knew no architecture above it. `python -m petakit5d.gpu.selftest` checks an installation.
