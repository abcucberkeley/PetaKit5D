"""What a Python task process costs before and after its work: the imports of a level, then
``mark('ready')`` (the first line of the work), nothing, and the end. Run by the harness
(``process_costs.json``), which marks the first line and times launch and exit.

    python startup_probe.py LEVEL

Levels:
  skeleton  the petakit5d package (loads nothing else)
  numpy     NumPy
  scipy     NumPy, SciPy, scipy.ndimage and scipy.fft
  io        NumPy, cpp-tiff, cpp-zarr and petakit5d.io (what a file task needs)
  full      io plus scipy.ndimage and scikit-image's filters
  cupy      NumPy and CuPy, then the GPU: device count, context, one elementwise kernel (from
            CuPy's kernel cache after the first run) and one cuFFT plan; each step is marked
"""

import sys

level = sys.argv[1] if len(sys.argv) > 1 else "skeleton"

if level == "skeleton":
    import petakit5d  # noqa: F401
elif level == "numpy":
    import numpy  # noqa: F401
elif level == "scipy":
    import numpy  # noqa: F401
    import scipy.fft  # noqa: F401
    import scipy.ndimage  # noqa: F401
elif level == "io":
    import cpptiff  # noqa: F401
    import cppzarr  # noqa: F401
    import numpy  # noqa: F401

    import petakit5d.io  # noqa: F401
elif level == "full":
    import cpptiff  # noqa: F401
    import cppzarr  # noqa: F401
    import numpy  # noqa: F401
    import scipy.ndimage  # noqa: F401
    import skimage.filters  # noqa: F401

    import petakit5d.io  # noqa: F401
elif level == "cupy":
    import cupy
    import numpy  # noqa: F401

    from petakit5d.testing.bench import info, mark

    mark("imported")
    n = cupy.cuda.runtime.getDeviceCount()
    mark("devices")
    cupy.cuda.Device(0).use()
    cupy.cuda.runtime.free(0)  # makes the context
    mark("context")
    x = cupy.arange(1024, dtype=cupy.float32)
    float((x * 2 + 1).sum())
    mark("kernel")
    cupy.fft.fftn(cupy.ones((64, 64, 64), dtype=cupy.complex64))
    cupy.cuda.Device(0).synchronize()
    mark("fft")
    props = cupy.cuda.runtime.getDeviceProperties(0)
    info(
        gpu=props["name"].decode(),
        gpus=n,
        cupy=cupy.__version__,
        cuda_runtime=cupy.cuda.runtime.runtimeGetVersion(),
        cuda_driver=cupy.cuda.runtime.driverGetVersion(),
    )
else:
    raise SystemExit(f"unknown level {level}")

from petakit5d.testing.bench import mark  # noqa: E402

mark("ready")
