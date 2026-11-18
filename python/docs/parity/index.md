# Differences from MATLAB: the index

The port aims at MATLAB's results bit for bit. On the reference platform (MATLAB R2024b, Linux
x86-64) with the exact FFT add-on it gets there for every output not listed below. This page
lists every known exception in one place: deliberate file-format choices, the tolerances of the
default (permissive) FFT backend, the few tolerances that remain with the exact backend, and
outputs that are not produced yet. Each row names the page or record with the measurements.

**Status:**
- **approved**: accepted by the user, either for this row or under the standing rule of
  2026-11-11 ("bit-identical if possible, a small tolerance OK" where exactness is impossible),
  on the date given;
- **proposed**: measured and documented, waiting for the user's decision;
- **temporary**: an output not produced, or a behaviour not matched, until a planned piece of
  work lands.

Tests name their row in the `reason` of their tolerance (`petakit5d.testing.compare.tol`,
[testing.md](../testing.md)). The area pages: [decon.md](decon.md), [demos.md](demos.md),
[flatfield.md](flatfield.md), [matlabcompat.md](matlabcompat.md), [tools.md](tools.md).

## 1. Files: deliberate format choices

These change what is on disk, not the values in it. [migration.md](../migration.md) says how to
read files across the two implementations.

| # | Output | Difference | Status | Source |
|---|---|---|---|---|
| F1 | zarr arrays the port writes | native `(z, y, x)` C order (shape `[nz, ny, nx]`, chunk files `iz.iy.ix`) instead of MATLAB's `(y, x, z)` Fortran order; MATLAB's arrays are still read | approved 2026-11-12 | [decisions/0004](../decisions/0004-io-bindings.md) |
| F2 | sharded zarr (`shardSize`) | standard Zarr v3 with sharding instead of MATLAB's v2 `.zarray` with a `sharding_indexed` codec | approved 2026-11-12 | [decisions/0004](../decisions/0004-io-bindings.md) |
| F3 | N5 the port writes | standard N5 (`dimensions` `[nx, ny, nz]`) instead of MATLAB's `[nz, nx, ny]`; MATLAB's N5 read with `axes='yxz'`; the BigStitcher XML's `data_order` defaults to `'zyx'` | approved 2026-11-17 | [decisions/0009](../decisions/0009-n5-layout.md) |
| F4 | `.mat` files with a saved inputParser (`parameters.mat`, `stitchInfo/<volume>.mat`) | the object `ip` is a struct of its results; the files are still written as v7.3 `.mat` for MATLAB | approved 2026-11-11 (plan, section 4.4) | [migration.md](../migration.md) |
| F5 | process-function files `tmp/processFunc/processFunction_*.mat` | `usrFun` is a char array, not a MATLAB string object | approved 2026-11-18 | [demos.md](demos.md), row 3 |
| F6 | figures | `.png` only, no `.fig` | approved 2026-11-11 (plan) | [migration.md](../migration.md) |
| F7 | job logs | `Python command:` where MATLAB's say `Matlab command:` | approved 2026-11-11 (plan) | [jobs.md](../jobs.md) |
| F8 | tiff in the `'libtiff'` and `'imwrite'` modes and uncompressed tiff | cpp-tiff's BigTIFF instead of classic TIFF through `saveastiff`; voxels and tags the same | approved 2026-11-11 (plan) | [decisions/0004](../decisions/0004-io-bindings.md) |
| F9 | the libtiff-mode writer's files (`replicated/*_conventional.tif` of the readers/writers demo) | the LZW bytes of each page's last, partial strip differ (about 1 KB a page); voxels and tags identical | approved 2026-11-18 | [demos.md](demos.md), row 2 |

## 2. The default FFT backend (pocketfft) against MATLAB's FFTW

With the exact add-on none of these apply: the results are MATLAB's bits.

| # | Output | Tolerance | Measured | Status | Source |
|---|---|---|---|---|---|
| P1 | FFT-based outputs in general | single: max <= 1e-5 of the largest value, relative L2 <= 1e-6; double: 1e-13 and 1e-14 | single <= 5.5e-7, double <= 1.2e-15 | approved 2026-11-12 | [decisions/0006](../decisions/0006-fft-backends.md) |
| P2 | deconvolution (up to about 30 iterations) and what is made from it (DSR, MIPs, conversions) | relative L2 <= 1e-4; uint16 off by at most 1 in at most 0.1% of voxels | demos: 0.0023-0.0083% of voxels off by 1; z MIPs up to 0.175% | approved 2026-11-12 | [decisions/0006](../decisions/0006-fft-backends.md); [demos.md](demos.md), row 1 |
| P3 | the widefield OMW decon of the widefield/confocal demo (20 iterations, `hann_win_bounds` (0.8, 1.2)) | outside P2 | 0.35% of voxels off by 1, 4,586 by more (largest 1,782), relative L2 1.7e-3 | **proposed** | [demos.md](demos.md), row 8 |
| P4 | the OMW back projector and its masked Wiener filter | max and relative L2 <= 1e-5; masks identical | | approved 2026-11-13 | tests of `decon` (`PERMISSIVE_BP`) |
| P5 | all-NaN back projectors (degenerate `hann_win_bounds`, `l == u`) | NaN positions compared, not NaN signs | | approved 2026-11-13 | `decon/backprojector.py` |
| P6 | a PSF resampled to one plane (`new_nz = 1`) | difference <= 1e-6 of the raw PSF's largest value | MATLAB's own result there is rounding noise (it takes the spectrum's last plane, not DC) | approved 2026-11-12 | `decon/psf.py`, tests of `decon/test_psf.py` |
| P7 | 3D registration (normalized cross-correlation) | the same integer shift where the correlation is well behaved; peak values within 5e-3 | shifts identical | approved 2026-11-12 | [decisions/0006](../decisions/0006-fft-backends.md) |
| P8 | chromatic shift estimation: correlation peak values | absolute 1e-4 | <= 1.8e-7; shifts and corrected files exact | approved 2026-11-18 | [tools.md](tools.md), chromatic shift |
| P9 | FFT power spectra: MIPs, central slices, 3D stack | the normalized spectrum `\|F\| / max` within P1's single tolerance | largest difference 4.8e-8 | approved 2026-11-17 | [tools.md](tools.md), FFT power spectra |
| P10 | PSF analysis: OTF slices, line cuts, info files | P1's double tolerance | <= 7.6e-16 of the largest value | approved 2026-11-12 | [tools.md](tools.md), PSF analysis |
| P11 | FSC curves | P1's double tolerance | <= 5.9e-15 of the largest value; radii identical | approved 2026-11-12 | [tools.md](tools.md), Fourier shell correlation |
| P12 | `deconvlucy` (RLdecon's `'original'` method) | P2 | | approved 2026-11-12 | [decon.md](decon.md) |
| P13 | 2D tile registration (`normxcorr2`): correlation peak values | absolute 1e-7 | shifts identical; peaks <= 5.3e-8 | **proposed** | stitch.md, with the tile registration (T1) |

## 3. Tolerances on both backends

MATLAB functions whose last bits could not be reproduced from their documentation and black-box
probes, and MATLAB results that depend on MATLAB's thread count.

| # | Output | Tolerance | Measured | Status | Source |
|---|---|---|---|---|---|
| N1 | single-precision `log10` (the FFT spectra's 3D stack, the figures' values) | within 2 ulp; the port's `log10` is correctly rounded, the same on every platform | 3,692 of 226,981 stack voxels differ, all within 2 ulp | approved 2026-11-17 | [tools.md](tools.md), FFT power spectra |
| N2 | `deconvlucy`'s double results (exact backend) | max <= 1e-13 of the largest value, relative L2 <= 1e-13 | <= 18 ulp (MATLAB's acceleration factor comes from MKL dot products) | **proposed** | [decon.md](decon.md) |
| N3 | FSC resolutions | each radius within one 1e-4 step of MATLAB's; resolutions within 1% | all 626 radii per backend equal MATLAB's | approved 2026-11-18 | [tools.md](tools.md), Fourier shell correlation |
| N4 | PSF detection (`radon`, `tand`, `atand` inside) | the decisions (peaks, crops, MIPs, info file) exact | bit for bit in all cases | approved 2026-11-17 | [tools.md](tools.md), PSF detection |
| N5 | filters: `imfilter` with zero padding or `'full'`, small 2-D `imgaussfilt` kernels, general 2-D and 3-D kernels, long-kernel `conv2`, the frequency domain and `edgetaper`, sparse arrays, uint16 along z | per row, 2 to 16 ulp, relative 1e-6 to 2e-6, integers 1 count | 352 of 399 cases bit-identical, the call sites of unmixing, puncta removal and PSF analysis included | approved 2026-11-13 | [matlabcompat.md](matlabcompat.md), filters |
| N6 | puncta removal: pixels filled by `fillmissing2` `'linear'` | 1 ulp, relative 1e-12, in at most 2% of values (uint16: 1 count in at most 1%) | 1 ulp in 26 of 10,647 probe fills | **proposed** | [matlabcompat.md](matlabcompat.md), puncta removal |
| N7 | puncta removal: background thresholds of batches whose y size is a multiple of 16 | relative 1e-6 | relative 6.1e-8 | **proposed** | [matlabcompat.md](matlabcompat.md), puncta removal |
| N8 | the phase-image flat field (BaSiCPy in place of MATLAB's BaSiC) | flat field: relative 0.06 (0.12 under strong vignetting), relative L2 0.01 (0.03); uint16 tiles 1 count more | about a percent; up to 10.7% at single pixels | approved 2026-11-18 | [flatfield.md](flatfield.md) |
| N9 | `median`, `mode`, `medfilt2`, `medfilt3` at ties of -0 and +0 | the sign of a zero may differ; values agree | 3 of 269 probe medians | approved 2026-11-12 | [decisions/0005](../decisions/0005-sums-and-math-library.md) |
| N10 | MATLAB's sums, means and standard deviations of 89,000 or more elements, which depend on MATLAB's thread count | the port gives MATLAB's one-thread result at any thread count | up to 1.2e-4 relative between MATLAB's own thread counts (single) | approved 2026-11-12 | [decisions/0005](../decisions/0005-sums-and-math-library.md) |
| N11 | demo_RL_deconvolution's in-memory results (the demo's `sum(psf, 'all')` at 24 threads) | N10: the one-thread sum | the back projector differs by <= 6.0e-8; uint16 results by 1 in <= 0.0082% of voxels; equal to MATLAB's with its one-thread sum | approved 2026-11-12 | [demos.md](demos.md), row 6 |
| N12 | the stitching's global shift assignment (`stitch_assignment_lsqlin`): its linear solves | absolute 1e-9 | sparse solvers (MATLAB's Cholesky and LU, SuperLU here) differ in the last bits | **proposed** | tests of `matlabcompat/test_lsq.py` (`ASSIGN_TOL`) |
| N13 | `sph2cart` of NaN inputs | NaN where MATLAB has NaN, NaN payloads not compared | | **proposed** | tests of `matlabcompat/test_neighbors.py` |

## 4. GPU

The GPU paths use cuFFT 11.0.8, the version MATLAB R2024b bundles, and are MATLAB's bits by
default ([decisions/0003](../decisions/0003-gpu.md)).

| # | Output | Tolerance | Measured | Status | Source |
|---|---|---|---|---|---|
| G1 | GPU deconvolution in the native `(z, y, x)` layout (an option; the default is MATLAB's layout, which is exact) | single: max and relative L2 <= 1e-5 (OTFs 1e-6); uint16 at most 1 count in at most 0.01% of voxels | relative L2 1.0e-6; 0.0016% of voxels | approved 2026-11-12 | tests of `gpu/test_decon_gpu_kernels.py` |

## 5. Not produced yet, and pending fixes

| # | Output | Until | Source |
|---|---|---|---|
| T1 | per-pair registration results `xcorr/<volume>/xcorr_tile_*.mat`, and registration itself | the tile registration is ported (planned: `cross_correlation_registration_wrapper`) | [demos.md](demos.md), row 4 |
| T2 | figures (back projector and OTF masks, FSC, PSF, MIP masks) | the `figures` extra is installed | [demos.md](demos.md), row 5 |
| T3 | the phase flat field (`ff_corrected_*/`, `_flatfield.mat`) | the `basic` extra is installed | [demos.md](demos.md), row 7 |
| T4 | a box with a NaN end (`[... NaN]`) is refused where MATLAB passes it on | the fix of 2026-11-18: accepted where MATLAB accepts it | [cli.md](../cli.md) |
