# The examples: the MATLAB demos in Python

`python/examples/phase1/` has one Python file per MATLAB demo (`demos/demo_*.m`). Each follows
its demo step by step with the same parameters, through the snake_case API; comments give
MATLAB's values where the Python form differs (boxes, sizes). They need no MATLAB.
`python/examples/phase2/` has the demos' Phase 2 steps, to run after them on the same data.

| Example | MATLAB demo | What it shows |
|---|---|---|
| `demo_generic_computing_framework.py` (with `crop_deskew_rotate_demo_function.py`) | demo_generic_computing_framework.m | a user function run on many files by the job framework, locally or as Slurm jobs; an example cluster configuration file |
| `demo_fast_tiff_zarr_readers_writers.py` | demo_fast_tiff_zarr_readers_writers.m | cpp-tiff and cpp-zarr readers and writers, reads into a reused array (`out=`) |
| `demo_geometric_transformation.py` | demo_geometric_transformation.m | combined and separate deskew and rotation of datasets and of large single files |
| `demo_rl_deconvolution.py` | demo_RL_deconvolution.m | PSF preparation, the OMW back projector, RL (Biggs) and OMW deconvolution, in memory and on datasets |
| `demo_skewed_space_stitching.py` | demo_skewed_space_stitching.m | stitching in skewed space with a tile processing function |
| `demo_zarr_stitching.py` | demo_zarr_stitching.m | zarr stitching in skewed and DSR space (the demo's BigStitcher-Spark and stitching-spark setup, after its `return`, is left out) |
| `demo_large_scale_processing.py` | demo_large_scale_processing.m | large-file deconvolution in place and deskew/rotation of stitched volumes |
| `demo_useful_tools.py` | demo_useful_tools.m | resampling, cropping, MIPs, tiff/zarr conversion and resaving, image lists, unmixing, and the decon and DSR steps they use |
| `demo_phase_and_2photon_stitching.py` | demo_phase_and_2photon_stitching.m | image lists from tile positions, 2D phase and 3D 2-photon stitching; the phase flat field with the `basic` extra |
| `demo_widefield_and_confocal_deconvolution.py` | demo_widefield_and_confocal_deconvolution.m | OMW deconvolution of widefield and confocal data |
| `demo_microscope_pipeline_single_dataset.py`, `..._nostitch.py` | the same | the microscope pipeline (`automatic_processing`): flat field, deskew/rotation and stitching of a dataset |
| `demo_microscope_pipeline_multiple_dataset.py` | the same | the pipeline on several datasets, given as arguments |
| `run_all_demos_for_the_paper.py` | demo_run_all_demos_for_the_paper.m | the paper's nine demos in order, in one process |
| `compare_demo_trees.py` | | compares a Python run's output tree with a MATLAB run's |

The Phase 2 examples, in `python/examples/phase2/`:

| Example | MATLAB demo | What it shows |
|---|---|---|
| `demo_useful_tools.py` | demo_useful_tools.m | its Phase 2 steps, on the results of the phase 1 example: the FFT power spectra, the PSF and FSC analyses and the chromatic shift correction; the Bio-Formats conversion (`bioformats_conversion`), not run, as in MATLAB |
| `run_all_demos_for_the_paper.py` | demo_run_all_demos_for_the_paper.m | the paper's nine demos in order with their Phase 2 steps, in one process |
| `compare_phase2_trees.py` | | compares the outputs of a Python run's Phase 2 steps with a MATLAB run's, step by step |

The phase 1 examples of the phase and 2-photon stitching and of the widefield and confocal
deconvolution have their demos' Phase 2 steps already: the phase flat field (with the `basic`
extra) and the OTF mask figures (with the `figures` extra).

## Running them

1. Download the demo datasets with the MATLAB demos' downloaders (`demos/demo_*_downloader.m`)
   or from the links in them, and copy
   them to a folder you can write to: the demos write their results into the dataset folders.
2. Point `PETAKIT5D_DEMO_DATA` at the folder that holds `PetaKit5D_demo_cell_image_dataset` and
   `PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets` (by default `~/Downloads/`, as in the
   MATLAB demos).
3. Run an example from its folder, or all the paper's demos:

```
export PETAKIT5D_DEMO_DATA=/mnt/projects/demo_copy
cd python/examples/phase1
python demo_rl_deconvolution.py
python run_all_demos_for_the_paper.py                         # all nine, in the paper's order
python run_all_demos_for_the_paper.py demo_useful_tools       # some of them
python demo_microscope_pipeline_multiple_dataset.py /mnt/projects/embryo_2024 /mnt/projects/embryo_2025
```

Later demos use the results of earlier ones (the tools use the deconvolved tiles), as in MATLAB,
so run them in the paper's order on one copy of the data. They need a large-memory machine: the
reference runs took about 20 minutes on 24 cores with 450 GB.

The Phase 2 examples run the same way from `python/examples/phase2/`, on the same copy of the
data and after the phase 1 examples whose results they use. Its `run_all_demos_for_the_paper.py`
runs each demo's phase 1 example and then its Phase 2 steps, in place of phase 1's:

```
cd python/examples/phase2
python run_all_demos_for_the_paper.py                         # all nine with their Phase 2 steps
python demo_useful_tools.py                                   # after phase 1's demo_useful_tools.py
```

Other variables:

| Variable | Effect |
|---|---|
| `PETAKIT5D_DEMO_TIMINGS` | a file the time of every step is appended to (one JSON object per line) |
| `PETAKIT5D_DEMO_SAVE_DIR` | a folder the in-memory results are written to as tiff files |
| `PETAKIT5D_FFT` | `exact` for MATLAB's bits with the FFTW add-on ([install.md](install.md)) |

The generic-framework example writes an example cluster configuration with placeholder values
to the current folder; it runs locally unless you edit it for your cluster
([cluster.md](cluster.md)).

## What the examples leave out

Steps whose functions are not in the port (or not yet) are skipped, and the example says so:
the conventional MATLAB tiff reader and the blockedImage and zarr-python readers, the legacy
separate deskew/rotation (`XR_deskewRotateFrame_old`), the Wiener-Butterworth back projector
(WBDeconvolution, GPL), the Spark stitching tools, and the Imaris conversion (its port waits for
the next release of Parallel_Imaris_Writer). The Bio-Formats conversion, which the data tools
demo shows after its `return`, is the phase 2 example's `bioformats_conversion`, not run (the
demo data has no such files). Figures are drawn only with the `figures` extra (the PSF analysis,
whose results are figures, needs it), the phase flat field only with the `basic` extra.
[parity/demos.md](../parity/demos.md) lists the steps left out per demo,
[parity/phase2.md](../parity/phase2.md) the Phase 2 steps and what they leave out.

## How close they are to MATLAB

The examples and the MATLAB demos were run in local mode on the same data and their output
trees compared file by file (`compare_demo_trees.py`): with the exact FFT add-on every output of
every ported step is MATLAB's bit for bit, apart from a few listed file differences. With the
default FFT backend the deconvolved outputs differ by 1 count in a few thousandths of a percent
of their voxels. The numbers and every difference are in [parity/demos.md](../parity/demos.md).
The Phase 2 steps were compared step by step (`compare_phase2_trees.py`): with the exact add-on
every output is MATLAB's bit for bit apart from the phase flat field (BaSiCPy in place of
MATLAB's BaSiC, within its parity row N8). On the default backend the FFT spectra and FSC
curves, made from that backend's decon results, are within a proposed row (P14), the other
outputs within their rows. The numbers are in [parity/phase2.md](../parity/phase2.md).
