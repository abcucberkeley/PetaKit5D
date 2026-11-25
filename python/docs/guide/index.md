# User guide

petakit5d is the Python port of PetaKit5D: deskew and rotation, deconvolution, stitching, the
conversions and tools, the microscope pipeline and the cluster job framework, with MATLAB's
results and MATLAB's files on disk. This guide covers:

- [install.md](install.md): the base package, the optional extras (`figures`, `n5`, `formats`,
  `basic`), the GPU packages, and the exact FFT add-on with its GPL note;
- [examples.md](examples.md): the Python versions of the MATLAB demos, and how to run them on
  the demo datasets;
- [cluster.md](cluster.md): running on a Slurm cluster with a cluster configuration file;
- coming from MATLAB, the GUI or PyPetaKit5D: [../migration.md](../migration.md).

## Three ways to call a function

Every user-facing function is declared once, in the parameter registry, and can be called in
three forms. [MATLAB_MAPPING.md](../MATLAB_MAPPING.md) lists each function by its MATLAB name with
its Python name and parameters.

**1. The Python API** (snake_case names, Python values):

```python
from petakit5d.deskew_rotate import deskew_rotate_dataset
from petakit5d.tools import crop_dataset

deskew_rotate_dataset(
    ["/mnt/projects/embryo_2024/"],
    xy_pixel_size=0.108,
    dz=0.3,
    skew_angle=32.45,
    deskew=True,
    rotate=True,
    dsr_combined=True,
    channel_patterns=["CamA_ch0"],
    save_16bit=True,
    parse_cluster=False,
)
# the box is (z0, y0, x0, z1, y1, x1), 0-based and end-exclusive
crop_dataset(
    ["/mnt/projects/embryo_2024/"],
    (100, 200, 100, 400, 1200, 412),
    result_dir_name="Cropped",
    parse_cluster=False,
)
```

**2. MATLAB's names and values** (`petakit5d.matlab`), as a MATLAB script or PyPetaKit5D calls
them: name-value pairs with MATLAB's parameter names in any case and unique prefixes, or keyword
arguments, and inputParser's errors:

```python
from petakit5d import matlab

# the box is MATLAB's [ymin xmin zmin ymax xmax zmax], 1-based and inclusive
matlab.XR_crop_dataset(
    "/mnt/projects/embryo_2024/",
    [201, 101, 101, 1200, 412, 400],
    "resultDirName",
    "Cropped",
    parseCluster=False,
)
```

**3. The command line:** `pk5d <function> ...` takes the Python names and values (JSON);
`pk5d-master <MATLAB name> <strings>` takes the strings the compiled MATLAB master takes, for
the GUI and scripts ([cli.md](../cli.md)). `pk5d doc <function>` prints a function's parameters,
defaults and MATLAB names.

## Values: Python's or MATLAB's, depending on the layer

| | Python API, `pk5d` | `petakit5d.matlab`, `pk5d-master`, config and parameter files |
|---|---|---|
| volumes | NumPy `(z, y, x)`, C order | (passed as files) |
| boxes | `(z0, y0, x0, z1, y1, x1)`, 0-based, end-exclusive | `[ymin xmin zmin ymax xmax zmax]`, 1-based, inclusive |
| sizes and per-axis vectors (block, batch, border sizes, resample factors) | `(z, y, x)` | `[y x z]` (or what MATLAB takes for that parameter) |
| parameter names | `xy_pixel_size`, `dsr_dir_name`, `save_16bit` | `xyPixelSize`, `DSRDirName`, `save16bit` |
| defaults | MATLAB's | MATLAB's |

The registry converts between the two forms at the layer's edge, so a parameter file written by
the port holds MATLAB's values, as MATLAB writes them. Details: [conventions.md](../conventions.md).

## Results

The port's outputs are MATLAB's bit for bit on Linux x86-64 with the exact FFT add-on, apart from
the exceptions listed in [parity/index.md](../parity/index.md). With the default (permissive)
FFT backend, results of FFT-based steps differ in the last bits, within the tolerances of
record 0006. The same code gives the same bits on Linux, macOS and Windows wherever the port
computes them itself ([platforms.md](../platforms.md)).

Some functions take two parameters MATLAB doesn't have (marked Python-only in `pk5d doc`):
`num_threads`, the threads of the call, and `out`, an array to reuse for the result where a
function is called frame after frame.

## Errors and warnings

Where MATLAB raises an error or gives a warning, the port raises a Python exception or gives a
Python warning whose `identifier` attribute is MATLAB's identifier, unchanged: MathWorks'
(`MATLAB:InputParser:UnmatchedParameter`, `images:imresize3:badScaleOrSize`) and those of
PetaKit5D's MATLAB code (`crop:inputError`, `PetaKit5D:mexFallback`), so code that checks them
carries over. The port's own errors and warnings, where MATLAB has none or the port refuses what
MATLAB takes, have identifiers in one style, `petakit5d:<area>:<name>`: the area is the package
or module that raises it, in lowercase, and the name is lowerCamelCase
(`petakit5d:params:conversionRefused`, `petakit5d:detection:amiraNotWritten`). `pk5d-master`
prints an error's identifier on its last line, as the compiled MATLAB master does
([cli.md](../cli.md)).
