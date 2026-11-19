# Differences from MATLAB: the production workflow on real data

`tests/e2e/test_realdata.py` runs the shape of a production run on the demo cell dataset (16
skewed light-sheet tiles: 4 tiles x 2 cameras x 2 time points, about 7 GB) through the MATLAB
wrappers' ports, by their MATLAB names and with the values a MATLAB script passes. MATLAB R2024b
runs the same calls from the same description. This page lists how the runs' output trees differ,
with numbers. Rows marked "proposed" wait for the user's approval.

## The workflow

1. the flat-field function files (`processFFCorrectionFrame` with the dataset's flat fields and
   camera backgrounds, a constant offset of 1 and a lower limit of 0.4), as a production script
   writes them;
2. `XR_tiffToZarr_wrapper` with the flat-field correction (`zarr_ff/`, 16 tiles);
3. `XR_matlab_stitching_wrapper` in skewed space from the raw tiffs: the flat field applied while
   the tiles are converted (`zarr/`), the CamB tiles of the first time point registered by
   cross-correlation (its 3 neighbouring pairs as jobs, `xcorr/`) and every volume placed by that
   result, feather blending (`matlab_stitch/`, four volumes of 6957 x 590 x 520 in MATLAB's
   order);
4. `XR_decon_data_wrapper`: large-file OMW decon (2 iterations) of the stitched volumes in place,
   its blocks as CPU jobs (`matlab_decon_omw/`) or on a GPU (`matlab_decon_omw_<gpu>/`);
5. `XR_deskew_rotate_data_wrapper`: large-file deskew/rotation of the deconvolved volumes;
6. `XR_MIP_wrapper`: the MIPs of the stitched volumes (zarr) and of the deskewed volumes (large
   zarr) along all three axes (`MIPs_wrapper/`).

The CPU steps ran in cluster mode, each wrapper's jobs submitted to Slurm from a cluster
configuration with the same Slurm parameters for both (MATLAB's or the port's launch string),
and for the port also in local mode. The GPU decon ran in local mode on one GPU.

## How the trees are compared

- Every file and zarr array either side wrote under the dataset folder, by name: tiff files by
  voxels, tags and bytes; zarr arrays by voxels and metadata (MATLAB's `(y, x, z)` Fortran
  arrays against the port's native ones); `.mat` files by value (`t`, `t1` and the inputParser
  object `ip` left out); JSON and text by content. Each run's root folder is masked.
- Not compared: job logs (job ids, hosts, times) and figures (`.png`, `.fig`: demos.md row 5).
- **Registration.** Both sides register the tiles themselves (the port with
  `petakit5d.stitch.registration`), and their registration results are compared: the per-pair
  files `xcorr/<volume>/xcorr_tile_*.mat` and the merged `xcorr/<volume>.mat`.

## Results

The whole workflow in cluster mode on CPUs, both sides registering the tiles (2026-11-18,
MATLAB from the oracle tree of 45b2bb1):

| Comparison | Outputs | Result |
|---|---|---|
| MATLAB (cluster) against the port (cluster), exact FFT add-on | 172 / 172 | identical apart from row 2 |
| MATLAB (cluster) against the port (cluster), default FFT | 172 / 172 | identical apart from rows 1, 2 and 3 |

The comparisons below are of the earlier runs of the same day (oracle 4a3bef3), made before the
registration was ported: the port's stitching was given MATLAB's merged registration results,
and the per-pair results were not compared.

| Comparison | Outputs | Result |
|---|---|---|
| the port in local mode against cluster mode (exact add-on) | 161 / 161 | identical apart from row 4, which MATLAB's modes have too |
| MATLAB (cluster) against the port (local), exact add-on | 161 / 161 | identical apart from rows 2 and 4 |
| GPU decon on a B200: MATLAB against the port, exact add-on | 40 / 40 | identical (row 5) |
| GPU decon on a B200: MATLAB against the port, default FFT | 40 / 40 | row 5 |

With the exact add-on every stitched volume, flat-field-corrected tile, registration result,
deconvolved volume (CPU and GPU), deskewed volume and MIP is MATLAB's, voxel for voxel and value
for value.

## Differences

| # | Output | Difference | Status | Why |
|---|---|---|---|---|
| 1 | the deconvolved volumes and what is made from them (their MIPs, DSR, DSR MIPs), and the OTFs and back projectors in `psfgen/`, on the default FFT backend | uint16 voxels off by 1 and never more: decon CamA 49-51 voxels of 2.1 G (0.0000024%), CamB 21,623-21,945 (0.0010%); DSR CamA 76-81 of 4.2 G, CamB 32,236-32,594 (0.00077%); MIPs at most 0.082% (a CamB y MIP). OTFs and back projectors (single) differ in 92-99.8% of their values, by at most 1.8e-7 (relative L2 1.4e-7 to 3.2e-7) | record 0006 | the FFT backend (pocketfft): within 0006 (decon: off by 1 in at most 0.1% of voxels; single transforms: max 1e-5 of the largest value, relative L2 1e-6). OMW with 2 iterations doesn't amplify it (compare row 6) |
| 2 | 12 zarr arrays (the stitched, deconvolved and deskewed volumes) | voxels and metadata identical; MATLAB writes 252 or 448 chunk files, the port 99-213: it leaves out the chunks that are all zero (MATLAB writes them; reading gives the same zeros) | decided (cpp-zarr's Z1) | until cpp-zarr's release with sparse writes, the port's zarr writes skip all-zero chunks (`petakit5d.io.convert`) |
| 3 | the registration's per-pair results `matlab_stitch/xcorr/<volume>/xcorr_tile_*.mat`, on the default FFT backend | the correlation peaks (`max_xcorr_mat`) of all 3 pairs off by 1.2e-7 to 1.8e-7 (2 to 3 single-precision ulp); the shifts, the merged `xcorr/<volume>.mat`, the stitch information and every stitched volume identical. All identical with the exact add-on | approved 2026-11-12 (P7) | the cross-correlation's FFTs on the default backend (record 0006); P7 allows the same shifts with peaks within 5e-3 |
| 4 | local against cluster mode: `stitchInfo/*_task_size_<n>*` and `PerBlockInfoFullpaths` in `stitchInfo/<volume>.mat`; `parseCluster` and `configFile` in the parameter files | the stitching's block information is made per blending task: one task of all 56 batches locally (`task_size_56`), tasks of 20 batches in cluster mode (`task_size_20`, three per-block files); the parameter files hold each call's own mode and configuration | as MATLAB | MATLAB's local and cluster modes differ the same way; every volume, MIP and tile is the same in both modes |
| 5 | the GPU decon (`matlab_decon_omw_<gpu>/`), B200 | exact add-on: bit for bit MATLAB's GPU results (4 volumes, their MIPs, `psfgen/`). Default FFT: uint16 off by 1 in at most 0.00055% of voxels (CamB 11,526-11,652 of 2.1 G, CamA 24-28), MIPs at most 0.043%, `psfgen/` as row 1 | record 0006 / R04 | cuFFT 11.0.8 gives MATLAB's bits; on the default backend only the back projector and OTF, made on the CPU, differ (row 1). For scale: MATLAB's own CPU and GPU results differ by 1 in up to 0.0012% of voxels (CamB 25,158-25,363), and so do the port's |
| 6 | OMW decon with many iterations on the default FFT backend (the widefield demo's 20; demos.md row 8) | at 20 iterations 0.35% of voxels off by 1 and 4,586 by more, the largest by 1,782; at 1-5 iterations only off by 1 (at most 0.009%), at 10 82 voxels by more (largest 68). Bit for bit with the exact add-on, on the CPU and on the GPU | approved 2026-11-19 (default backend, parity row P3) | inherent to the method: 20 multiplicative OMW iterations amplify rounding-level differences in a few bright spots they build. MATLAB's own CPU and GPU results differ the same way (20 iterations: 3,061 voxels by more, largest 817), and so does the exact add-on (MATLAB's bits) with the PSF moved by one ulp in a single value, 2.4e-12 of its sum (2,369 voxels, largest 807). The differing voxels sit in the same few bright spots in all three cases. MATLAB's CPU result is as far from the same OMW in exact arithmetic, and no option short of the exact add-on reduces the difference (X14: [decon.md](decon.md), OMW with many iterations) |
