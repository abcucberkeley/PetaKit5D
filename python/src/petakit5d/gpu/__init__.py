"""GPU support: device choice, run-time compiled CUDA kernels, cuFFT plans, the memory pool.

CuPy is optional. Importing this package or any of its modules imports neither CuPy nor
NumPy; the functions that need a GPU import CuPy, and the ones that only ask whether a GPU
can be used (and why not) also work without it.

Modules:
    device    which GPU this process may use (usable_gpu_count, select_device, gpu_for)
    kernels   CUDA sources compiled at run time by NVRTC and cached on disk (load_module)
    fft       single precision real-to-complex cuFFT plans with one shared work area, a plan cache
    memory    the memory pool between tasks, free memory, peak use
    selftest  a check of the installation: python -m petakit5d.gpu.selftest

The decisions behind this are in docs/decisions/0003-gpu.md.
"""
