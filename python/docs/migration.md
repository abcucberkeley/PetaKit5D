# Moving from MATLAB PetaKit5D to petakit5d

For users of MATLAB PetaKit5D, of PetaKit5D-GUI and of PyPetaKit5D. The port runs the same
functions with the same parameters and defaults, and writes the same folders and files, so most
work moves over without changes. This page lists what to change and what is different.

## In short

| | What to do |
|---|---|
| MATLAB scripts | call the same functions from Python, under their MATLAB names (`petakit5d.matlab`) or their Python names ([MATLAB_MAPPING.md](MATLAB_MAPPING.md)) |
| cluster configuration files | keep them; add `PythonLaunchStr` |
| PetaKit5D-GUI | replace the GUI's `run_mccMaster.sh` with the port's shim (`pk5d mcc-shim`) |
| PyPetaKit5D | the same shim; or call `petakit5d.matlab` with the same keyword arguments |
| data processed by MATLAB | read as before; zarr and N5 written by the port have another axis order (below) |
| results | MATLAB's bits with the exact FFT add-on, apart from the listed exceptions ([parity/index.md](parity/index.md)) |

## Functions and parameters

Each MATLAB function of the port keeps its parameters, their names and MATLAB's defaults.
[MATLAB_MAPPING.md](MATLAB_MAPPING.md), generated from the parameter registry, lists every
function: its Python name, the MATLAB name, every parameter in both names with its default, and
the functions that are not ported, each with its reason. Names: the `XR_` prefix is dropped,
`Frame` functions are `*_file`, wrappers `*_dataset`, and parameters are MATLAB's names in
snake_case (`xyPixelSize` is `xy_pixel_size`, `DSRDirName` is `dsr_dir_name`).

A MATLAB call carries over as it is through the MATLAB-name layer:

```python
from petakit5d import matlab

matlab.XR_deskew_rotate_data_wrapper(
    "/mnt/projects/embryo_2024/", "xyPixelSize", 0.108, "dz", 0.3, "skewAngle", 32.45
)
```

Cell arrays are lists, structs are dicts, logicals are `True`/`False`; boxes, sizes and indices
keep MATLAB's form in this layer. The snake_case API takes the Python forms instead
(`(z, y, x)` sizes, 0-based end-exclusive boxes; see [guide/index.md](guide/index.md)).

Functions that MATLAB takes as code text (a tile processing function in a `.mat` or `.txt`
file, `processFunPath`) are read, not evaluated: the standard forms (the edge erosion
`erodeVolumeBy2DProjection`, the flat-field correction `processFFCorrectionFrame`, their
composition and constant offsets) are recognized; anything else needs a Python function, given
as `"package.module:function"`. The generic framework's tasks are Python functions or bash
commands, not MATLAB code strings.

Not ported in Phase 1 (also in the mapping table's "Not ported" list where they are MATLAB
functions):

| | Why | Instead |
|---|---|---|
| point detection, tracking, lifetime analysis, SIM reconstruction | later phases | MATLAB |
| the Imaris conversion (`XR_imaris_conversion_data_wrapper`) | planned (Phase 2, a wrapper of the Parallel_Imaris_Writer binary) | MATLAB, or the binary |
| the Wiener-Butterworth back projector of demo_RL_deconvolution (`BackProjector`) | WBDeconvolution is GPL | OMW (`RLMethod 'omw'`) |
| the Java and Spark stitching wrappers, the cudaDecon and cppDecon wrappers | external tools | the tools themselves |
| legacy code without callers (`XR_stitching_frame_v1`, `XR_deskewRotateFrame_old`, ...) | not used | the current functions |
| `setup.m`, the mcc and mex machinery, MATLAB-only readers (`readtiff_matlab`, blockedImage, the zarr-python adapters) | MATLAB-only | `pip install`; `petakit5d.io` |
| `get_uuid` | a utility | `petakit5d.core.ids.get_uuid` |

## Renamed before the first release

Three modules were renamed on 2026-11-18, so that a function's module names its MATLAB folder or
kind, as the names above do (`Frame` functions in `*_file` modules;
[decisions/0010](decisions/0010-module-names.md)). The functions, their parameters and their
MATLAB names stay as they were; only these module paths changed:

| MATLAB | Was | Now |
|---|---|---|
| `XR_RotateFrame3D` | `petakit5d.deskew_rotate.volume.rotate_file` | `petakit5d.deskew_rotate.file.rotate_file` |
| `XR_RLdeconFrame3D` | `petakit5d.decon.frame.decon_file` | `petakit5d.decon.file.decon_file` |
| `RLdecon` | `petakit5d.decon.file.rl_decon_file` (and `RLDeconResult`) | `petakit5d.decon.rl_file.rl_decon_file` |
| `XR_generate_image_list_wrapper` and its generators | `petakit5d.stitch.imagelist_generate` | `petakit5d.stitch.imagelist_generation` |

The names the packages export (`petakit5d.deskew_rotate.rotate_file`,
`petakit5d.decon.decon_file`, `petakit5d.decon.rl_decon_file`,
`petakit5d.stitch.generate_image_list`, ...) and the MATLAB-name layer (`petakit5d.matlab`)
did not change. The port had not been released, so the old module paths are not kept: there are
no compatibility shims. Code or a task record that names a module by its old path
(`"petakit5d.decon.frame:decon_file"`) takes the new one; note that `petakit5d.decon.file` now
holds `decon_file`, not `rl_decon_file`.

## Cluster configuration files

The files stay as they are: JSON (or `.mat`) with the framework's parameters in MATLAB's names.
Add one key, `PythonLaunchStr`, what a job runs to start Python (module loads, the environment,
then the interpreter; empty: the submitting interpreter):

```json
{
  "PythonLaunchStr": "source /opt/envs/petakit5d/bin/activate; python",
  "SlurmParam": "-p compute --qos normal -n1 --mem-per-cpu=21418M",
  "jobLogDir": "/mnt/projects/job_logs",
  "maxCPUNum": 24, "memPerCPU": 20.9, "wholeNodeJob": true, "GNUparallel": true
}
```

`MatlabLaunchStr`, `MCRParam`, `MCCMasterStr`, `MCRCacheRoot` and `mccMode` are ignored by the
port, and MATLAB ignores `PythonLaunchStr`, so one file serves both while they run side by side.
Job logs keep their names and layout and say `Python command:` where they said
`Matlab command:`. More: [guide/cluster.md](guide/cluster.md).

## PetaKit5D-GUI and PyPetaKit5D: the switch

Both start the compiled MATLAB master as
`<install>/PetaKit5D/mcc/linux/run_mccMaster.sh <MATLAB Runtime folder> <function> <strings>`
(macOS: `mac/run_mccMaster.sh`). `pk5d-master` takes the same function names and argument
strings, and the port writes a launcher with that calling convention:

```
pk5d mcc-shim /opt/PetaKit5D-GUI/PetaKit5D/mcc/linux/run_mccMaster.sh
```

1. Keep a copy of the original file (to switch back).
2. Run `pk5d mcc-shim <path>` from the environment that has petakit5d (or pass
   `--command`). The launcher ignores the Runtime folder, handles the callers' quoting as
   MATLAB's launcher does, and runs `python -m petakit5d.cli.master`.
3. For the GUI, do the same for `linux_with_jvm/run_mccMaster.sh`, which it uses for the OTF mask
   figures, the FSC analysis and the PSF analysis.

Nothing in the GUI changes; its jobs now run in Python. Functions that are not ported (the
Imaris conversion, SIM reconstruction) give MATLAB's "The parser for function ... does not
exist!" error with the reason: keep the original launcher for them. On Windows the GUI starts
`mcc/windows/mccMaster.exe` without a Runtime folder: a copy of the `pk5d-master.exe` console
script under that name has the same convention (not tested yet). The full grammar and the
tests of the GUI's and PyPetaKit5D's calls are in [cli.md](cli.md).

**PyPetaKit5D without the MATLAB Runtime:** its wrappers take MATLAB's parameter names as keyword
arguments, and so do the functions of `petakit5d.matlab`, so a script changes only its import:

```python
from petakit5d.matlab import XR_decon_data_wrapper

XR_decon_data_wrapper(
    ["/mnt/projects/embryo_2024/"],
    xyPixelSize=0.108,
    dz=0.3,
    psfFullpaths=["/mnt/projects/psfs/488.tif"],
    deconIter=2,
    channelPatterns=["CamA_ch0"],
)
```

A task starts in about 0.1 s instead of the compiled master's 3 to 4 s, and the port runs on
Linux, macOS and Windows.

## Files: what stays and what differs

Everything a user, the GUI or a later stage reads keeps its name, place and format: the result
folders (`DS/`, `DSR/`, `matlab_decon/`, `matlab_stitch/`, `MIPs/`, ...), tiff files (BigTIFF,
LZW, byte for byte the mex writer's), MIP names and orientation, `parameters.json` as MATLAB's
`jsonencode` prints it, flag folders, temporary names and job logs. The differences are
deliberate choices or follow from the language:

| File | MATLAB | petakit5d | Reading across |
|---|---|---|---|
| zarr arrays written | `(y, x, z)`: shape `[ny, nx, nz]`, `"order": "F"`, chunk files `iy.ix.iz` | native `(z, y, x)`: shape `[nz, ny, nx]`, `"order": "C"`, chunk files `iz.iy.ix` (what zarr-python, napari and OME-Zarr readers expect) | the port reads MATLAB's arrays (detected by `"order": "F"`); MATLAB reads the port's with its axes permuted ([decisions/0004](decisions/0004-io-bindings.md)) |
| sharded zarr (`shardSize`) | Zarr v2 `.zarray` with a `sharding_indexed` codec | standard Zarr v3 with sharding (`zarr.json`) | the port reads MATLAB's through cpp-zarr |
| N5 | `dimensions` `[nz, nx, ny]`, chunk files `iz/ix/iy` | standard N5 `[nx, ny, nz]`, x fastest (what Fiji, BigStitcher and BigDataViewer expect) | read MATLAB's N5 with `axes='yxz'` (it can't be detected); the BigStitcher XML's `data_order` defaults to `'zyx'` ([decisions/0009](decisions/0009-n5-layout.md)) |
| `.mat` files (`parameters.mat`, `stitchInfo/`, `xcorr/`, block information, flat-field and process-function files) | MATLAB's | still written, as v7.3 files MATLAB loads, for compatibility; the port reads MATLAB's | a saved inputParser object (`ip`) becomes a struct of its results; a function text (`usrFun`) is a char array, not a string object |
| figures | `.png` and `.fig` | `.png` only, with the `figures` extra | |
| tiff in the `'libtiff'` and `'imwrite'` modes, and uncompressed | classic TIFF through `saveastiff` | cpp-tiff's BigTIFF; the voxels and tags are the same | both read either |
| job logs | `Matlab command:` | `Python command:` | |

**Mixing MATLAB and Python stages** on the same data: tiff passes between them unchanged; zarr
and N5 written by one are read by the other with their axes in another order. A dataset that
MATLAB processed partly may hold arrays of both layouts; the port reads each by its own
metadata. Resuming a half-finished MATLAB run in Python works where a documented input carries
it (a stitch information file), not in general.

## Results

The port reproduces MATLAB's arithmetic, including its rounding, its sums' order and its quirks,
so outputs are MATLAB's bit for bit on Linux x86-64 with the exact FFT add-on. Without the add-on
(the default, so that the package stays BSD-licensed) FFT-based results (deconvolution,
registration peaks, spectra, FSC curves) differ in the last bits within the approved tolerances.
Every known difference, with its size and status, is in [parity/index.md](parity/index.md).
