# Differences from MATLAB: the demos' Phase 2 steps

The steps of the MATLAB demos that need Phase 2 of the port are examples too.
`examples/phase2/demo_useful_tools.py` has those of `demos/demo_useful_tools.m`, to run after
`examples/phase1/demo_useful_tools.py` on the same data (it uses that example's decon results),
and `examples/phase2/run_all_demos_for_the_paper.py` runs the paper's demos with them. The
phase 1 examples of `demo_phase_and_2photon_stitching.m` and
`demo_widefield_and_confocal_deconvolution.m` have their Phase 2 steps already, with the
`basic` and `figures` extras. This page compares these steps' outputs with MATLAB's, with
numbers. Rows marked "proposed" wait for the user's approval.

## The steps

| Demo | Step (MATLAB) | Python | Example | Status |
|---|---|---|---|---|
| demo_useful_tools | FFT analysis (`XR_fftSpectrumComputingWrapper`) | `tools.fft_spectrum_dataset` | phase2: `fft_analysis` | ported |
| | PSF analysis (`XR_psf_analysis_wrapper`) | `tools.psf_analysis_dataset` | phase2: `psf_analysis` | ported; needs the `figures` extra (its results are figures) |
| | FSC analysis (`XR_FSC_analysis_wrapper`) | `tools.fsc_analysis_dataset` | phase2: `fsc_analysis` | ported; figures with the `figures` extra |
| | Imaris converter (`XR_imaris_conversion_data_wrapper`) | (planned: `tools.imaris_conversion_dataset`) | left out | **pending**: its port waits for the next release of Parallel_Imaris_Writer (row 4) |
| | Chromatic Shift Correction (`XR_chromatic_shift_correction_data_wrapper`) | `tools.chromatic_shift_correction_dataset` | phase2: `chromatic_shift_correction` | ported |
| | the Bio-Formats conversion, after the demo's `return` (`XR_bioformats_to_tiff_or_zarr_wrapper`) | `io.bioformats.bioformats_to_tiff_or_zarr_dataset` | phase2: `bioformats_conversion` | ported for ND2 and CZI (the `formats` extra); not run, as in MATLAB (no demo images) |
| demo_phase_and_2photon_stitching | flat field correction for phase images (`XR_phase_image_flat_field_correction`) | `flatfield.phase_image_flatfield_correction` | phase1 | ported; needs the `basic` extra |
| demo_widefield_and_confocal_deconvolution | Steps A2 and B2, the OTF masks (`XR_visualize_OTF_mask_segmentation`) | `tools.visualize_otf_mask_segmentation` | phase1 | ported; a figure, no files; needs the `figures` extra |

Every function these three demos call is in the mapping table
([MATLAB_MAPPING.md](../MATLAB_MAPPING.md)): the Phase 2 ones registered,
`XR_imaris_conversion_data_wrapper` planned.

## How the outputs are compared

- **The runs.** MATLAB R2024b (the oracle tree of 2026-11-18) ran the three demos in local mode
  (Slurm hidden, as in `tests/+pktest/NoSlurmFixture`) in one session, on a fresh copy of both
  demo datasets with the `DSR/` zarr tiles that `demo_geometric_transformation` writes in the
  paper's order (demo_useful_tools' zarr unmixing and resaving read them; the same copy for
  both sides, counted as inputs). The examples ran the same way
  (`examples/phase2/run_all_demos_for_the_paper.py` with the three demos, on 24 CPUs), with the
  `figures` and `basic` extras, once with the exact FFT add-on and on the default backend.
- **The comparison** (`examples/phase2/compare_phase2_trees.py`, on
  `petakit5d.testing.compare`): every output of a Phase 2 step, by name; tiff files by voxels
  and tags, `.mat` files by value (`t`, `t1` and the inputParser object `ip` left out), figures
  (`.png`, `.pdf`) by name only (their pixels are matplotlib's). Each run's root folders and
  MATLAB's `~` are masked. The phase 1 steps' outputs of the same runs were compared with
  `examples/phase1/compare_demo_trees.py` (below).
- **The steps alone.** On the default backend the FFT spectra and the FSC curves are made from
  that run's decon results, which differ from MATLAB's as [demos.md](demos.md) row 1 says. So the
  two steps were also run on MATLAB's own decon results (the deskewed tiles and the decon
  iterations of MATLAB's run) and compared with `--matlab-inputs`, which measures the steps
  themselves.

## Results

With the exact FFT add-on every output of every ported Phase 2 step is MATLAB's bit for bit,
apart from the phase flat field (row 1, BaSiCPy in place of MATLAB's BaSiC):

| Demo: step | MATLAB outputs | Compared | Identical (exact add-on) | Figures, by name | Not compared |
|---|---|---|---|---|---|
| useful_tools: FFT analysis | 57 | 57 | 57 | | |
| useful_tools: PSF analysis | 8 | 8 | 4 | 4 (2 `.png`, 2 `.pdf`) | |
| useful_tools: FSC analysis | 25 | 23 | 11 (the 10 iterations' results, their resolution table) | 12 `.png` | 2 `.fig` (parity row F6) |
| useful_tools: Imaris converter | 4 | | | | 4 (row 4) |
| useful_tools: chromatic shift correction | 16 | 16 | 16 | | |
| phase_and_2photon_stitching: phase flat field | 42 | 42 | 0 (row 1) | | |
| widefield_and_confocal_deconvolution: Steps A2, B2 | (no files) | | | | |
| the phase 1 steps' figures (the decon wrappers' back projectors) | 4 | 4 | | 4 `.png` | |

On the default FFT backend: the PSF analysis within P10, the chromatic shift correction
identical, the phase flat field as with the add-on (row 1); the FFT spectra and the FSC curves
within their own rows on MATLAB's inputs and within row 2 end to end; the FSC resolutions
within N3 (row 3).

The phase 1 steps of these runs (`examples/phase1/compare_demo_trees.py`, as in
[demos.md](demos.md)): with the exact add-on all of their outputs are MATLAB's except the process-function file of the 2-photon stitching
(`usrFun` as char, parity row F5); on the default backend the decon outputs differ as
[demos.md](demos.md) rows 1 and 8 say. The per-pair registration results
(`xcorr/<volume>/xcorr_tile_*.mat`, demos.md row 4), compared now that the port registers the
tiles: the 7 of the 2-photon stitching (3D) identical with the add-on and differing only in
their correlation peaks on the default backend (at most 1.8e-7, parity row P7); the 38 of the
phase stitching (2D, `normxcorr2`) differing only in their peaks, by at most 2.5e-13 on both
backends (parity row P13); every shift identical.

## Differences

| # | Output | Difference | Status | Why |
|---|---|---|---|---|
| 1 | `Phase/ff_corrected_1/` (40 uint16 tiles) and its 2 `BaSiC_flatfield/*_flatfield.mat` (demo_phase_and_2photon_stitching) | the flat fields differ by at most 5.6% (99th percentile 2.8%, median 0.33%), relative L2 0.64%; `minval` by at most 0.002; the tiles in 98.6% of their values, by at most 736 counts (5.3%; 99th percentile 2.7%), relative L2 0.65% | approved (N8, 2026-11-18) | BaSiCPy (MIT) in place of MATLAB's BaSiC (not ported: CC BY-NC-ND): the same method, other solvers ([flatfield.md](flatfield.md)); the same numbers as the measurement there |
| 2 | the FFT spectra (48 MIPs and central slices, gamma 0.5) and the FSC curves (`fsc_mu`, `fsc_cell` of 10 files) of demo_useful_tools, on the default backend, end to end | the spectra in 99.9% of their values, by at most 3.5e-5 (relative L2 at most 5.4e-5; as `\|F\|/max` the largest difference within P9's 1e-5, the relative L2 up to 1.0e-5, beyond P9's 1e-6); the FSC curves in all values, by at most 6.9e-6 (relative L2 at most 2.3e-6, beyond P11). On MATLAB's own decon results the steps are within their rows: the spectra by at most 5.8e-6 (relative L2 1.2e-5; within P9), the curves by at most 8.9e-16 of the largest value (within P11) | proposed (P14) | the inputs: the default backend's decon results (deskewed for the spectra) differ from MATLAB's as [demos.md](demos.md) row 1 says, and the steps carry that on. Proposed: P2's relative L2 (at most 1e-4) for what is made from decon results covers them; the comparison applies it end to end, and P9 or P11 with `--matlab-inputs` |
| 3 | the FSC resolutions of demo_useful_tools, on the default backend, end to end | 2 of the 240 radii of the four sub-lattice pairs (`res_cell`: 6 directions, 10 iterations) one step of 1e-4 away, their resolutions by at most 1.1e-4 relative; the means (`res_mu`) and the resolution table identical; on MATLAB's own decon results all identical | approved (N3) | row 2's inputs: the curves meet the threshold one grid step apart |
| 4 | `matlab_decon_omw/DSR_for_useful_tools/imaris/` (the `.ims` file, `parameters.json`, `parameters.mat`, `timepoints.txt`) | not made | temporary (parity row T5) | the port of `XR_imaris_conversion_data_wrapper` waits for the next release of Parallel_Imaris_Writer |
| 5 | MATLAB's `.fig` files (2, FSC) | not made | approved (F6) | the port writes `.png` figures only |

## Not compared

- the Imaris converter's outputs (row 4) and MATLAB's `.fig` files (row 5);
- the pixels of the figures (their names are);
- the Bio-Formats conversion: MATLAB's demo has it after its `return`, for files the demo data
  does not have; `bioformats_conversion` in the example is the same call, not run.
