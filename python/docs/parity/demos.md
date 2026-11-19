# Differences from MATLAB: the demos

The examples in `examples/phase1/` follow the MATLAB demos (`demos/demo_*.m`) step by step with
the same parameters, through the snake_case API. This page lists how their output trees differ
from MATLAB's, with numbers. Rows marked "proposed" wait for the user's approval.

## How the trees are compared

- **The runs.** MATLAB R2024b ran the demos in local mode (Slurm hidden, as in
  `tests/+pktest/NoSlurmFixture`): the 9 demos of `demo_run_all_demos_for_the_paper.m` in one
  session on one copy of both demo datasets, in the paper's order (later demos use earlier
  results), then `demo_microscope_pipeline_single_dataset` and its `_nostitch` variant, each on a
  fresh copy of the cell dataset. The examples ran the same way
  (`run_all_demos_for_the_paper.py`, then the two pipeline examples), on 24 CPUs.
- **The comparison** (`examples/phase1/compare_demo_trees.py`, on
  `petakit5d.testing.compare`): every file and zarr array MATLAB wrote, by name; tiff files by
  voxels and tags, zarr arrays by voxels and metadata (MATLAB's `(y, x, z)` Fortran arrays against
  the port's native ones, as decided on 2026-11-12); `.mat` files by value (`t`, `t1` and the
  inputParser object `ip` left out); JSON and text by content. Each run's root folders are
  masked, and MATLAB's `~` is the home folder (the demos pass `~/Downloads/...` on). Empty
  folders are not compared (flag folders; MATLAB's harness keeps no empty folders).
- **Registration.** These runs predate the port's tile registration (2026-11-18), so the
  stitching steps that register their tiles were given MATLAB's merged results
  (`xcorr/<name>.mat`); everything else in those steps is the port's. The examples now register
  their tiles themselves, as the MATLAB demos do (row 4).

## Results

With the exact FFT add-on, every output of every Phase 1 step is MATLAB's bit for bit, apart
from rows 2 to 4 below:

| Demo | Outputs compared | Identical (exact add-on) |
|---|---|---|
| demo_generic_computing_framework | 16 | 16 |
| demo_fast_tiff_zarr_readers_writers | 3 | 2 (row 2) |
| demo_geometric_transformation | 32 | 32 |
| demo_RL_deconvolution | 119 (+ 6 in-memory results) | 119 (+ 6, row 6) |
| demo_skewed_space_stitching | 41 | 39 (row 3) |
| demo_zarr_stitching | 31 | 31 |
| demo_large_scale_processing | 70 | 70 |
| demo_useful_tools | 191 | 191 |
| demo_phase_and_2photon_stitching | 84 | 83 (row 3) |
| demo_widefield_and_confocal_deconvolution | 16 | 16 |
| demo_microscope_pipeline_single_dataset | 70 | 70 |
| demo_microscope_pipeline_single_dataset_nostitch | 32 | 32 |

With the default (permissive) FFT backend the deconvolved outputs differ as rows 1 and 8 say;
everything else is identical there too.

## Differences

| # | Output | Difference | Status | Why |
|---|---|---|---|---|
| 1 | deconvolved tiles and volumes and what is made from them (their DSR, MIPs and conversions), on the default FFT backend: demo_RL_deconvolution, demo_large_scale_processing, demo_useful_tools, the confocal decon | uint16 voxels off by 1 and never more: 0.0083% (RL demo, 64 outputs), 0.0023% (large-scale, 56), 0.0029% (tools, 52), 0.0047% (confocal); per decon result at most 0.0102% (RL, 30 iterations); z MIPs up to 0.175%. The single-precision OTFs and back projectors in `psfgen/` differ in most values by up to 3.0e-7 | record 0006 | the FFT backend (pocketfft): within 0006's rule for decon results (off by 1 in at most 0.1% of voxels) |
| 2 | `replicated/*_conventional.tif` (the libtiff-mode writer, demo_fast_tiff_zarr_readers_writers) | voxels and tags identical; the LZW bytes of the last, partial strip of each page differ (MATLAB 119,946 bytes, the port 119,057 on the first page): the files are 8,092,569,376 and 8,081,857,536 bytes | approved | the two writers' LZW streams of the partial strip (264 of 512 rows) differ, by about 1 KB a page; the full strips have the same offsets and sizes. Readers get the same voxels. Follow-up for the io owner |
| 3 | `tmp/processFunc/processFunction_*.mat` (demo_skewed_space_stitching, demo_phase_and_2photon_stitching) | `usrFun` is a char array; MATLAB's demo saves a string object (`sprintf` with a double-quoted format) | approved | `petakit5d.io.write_mat` writes no MATLAB string objects. The port and MATLAB read both forms as the same function text |
| 4 | `xcorr/<volume>/xcorr_tile_*.mat` (the stitching steps with registration) | not compared: in these runs the port was given MATLAB's merged results and wrote no per-pair results | until the demos run again | the registration is ported since (`petakit5d.stitch.registration`) and writes them: on the cell dataset in skewed space they are MATLAB's bit for bit with the exact FFT add-on, with the same shifts on the default backend ([realdata.md](realdata.md)); the 2D registration of the phase tiles uses `normxcorr2` ([stitch.md](stitch.md)) |
| 5 | `*.png`, `*.fig` (back-projector figures of the decon wrappers, FSC figures) | not drawn | the figures extra | the port draws the decon figures with matplotlib (the optional `figures` extra), which was not installed |
| 6 | demo_RL_deconvolution's in-memory results (`decon_lucy`, `omw_backprojector`, `decon_lucy_omw` on one tile) | the back projector differs in 99.7% of its single values (largest difference 6.0e-8), and the uint16 results by 1 in 0.0036% (OMW) and 0.0082% (RL) of voxels, with the exact add-on | accepted (N05, 2026-11-12, item 3) | the demo normalizes the PSF by `sum(psf, 'all')`, which MATLAB computes at 24 threads as 4747003 and at 1 thread as 4747002.5; the port gives the 1-thread value at any thread count. With the 1-thread sum MATLAB's results are the port's bit for bit (exact add-on) |
| 7 | `Phase/ff_corrected_1/` and the `_flatfield.mat` files (demo_phase_and_2photon_stitching) | not made | the basic extra | the phase flat field needs BaSiCPy (the optional `basic` extra, record of B07), not installed here |
| 8 | the widefield OMW decon (demo_widefield_and_confocal_deconvolution, 20 iterations, `hann_win_bounds` (0.8, 1.2)), on the default FFT backend | 0.353% of uint16 voxels off by 1, 4,586 off by more, the largest by 1,782 (MATLAB 22,714); relative L2 1.7e-3; its z MIP 3.3% | approved 2026-11-19 (default backend, parity row P3); outside record 0006 | last-bit differences grow in a few bright spots from about 8 OMW iterations on this data (the exact add-on is bit for bit). MATLAB's own CPU result is as far from the same OMW in exact arithmetic (relative L2 1.1e-3 at 20 iterations) and from MATLAB's GPU result; MATLAB's GPU result (1.3e-4) and the default backend's iterations (9.6e-5) are close to exact arithmetic. Nothing short of the exact add-on reduces the difference ([decon.md](decon.md), OMW with many iterations). RL on this dataset stays off by 1 |

## Steps left out

These steps' functions are not in Phase 1 of the port (or not ported at all); the examples skip
them and say so, and their MATLAB outputs are not compared.

| Demo | Step | Why |
|---|---|---|
| demo_fast_tiff_zarr_readers_writers | `readtiff_matlab` (the conventional reader), zarr through blockedImage and zarr-python | MATLAB only |
| demo_geometric_transformation | the separate deskew and rotation `XR_deskewRotateFrame_old` (6 files) | legacy code, not ported |
| demo_RL_deconvolution | the Wiener-Butterworth back projector (`BackProjector`) and its decon | WBDeconvolution is not ported |
| demo_zarr_stitching | BigStitcher-Spark and stitching-spark (after the demo's `return`) | Spark tools |
| demo_useful_tools | FFT spectra (57 files), PSF analysis (8), FSC analysis (25), the Imaris converter (4), the chromatic shift correction (16), Bio-Formats (after the `return`) | Phase 2 steps: `examples/phase2/demo_useful_tools.py`, compared in [phase2.md](phase2.md); the Imaris converter is not ported yet (T5) |
| demo_microscope_pipeline_multiple_dataset | the whole demo | it points at the lab's own data; the example takes the dataset folders as arguments |

## Notes

- Parameters MATLAB's demos pass that the wrappers ignore (`Reverse` and `fixIter` for
  `XR_decon_data_wrapper`) are left out of the examples.
- Paths in parameter files are the paths the caller gave, as in MATLAB. The MATLAB demos pass
  `~/Downloads/...` on unexpanded; the examples expand `~` in `PETAKIT5D_DEMO_DATA` (and in their
  default, `~/Downloads/`), so their files hold absolute paths. The comparison masks both.
