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
- **documented**: a property of the default backend, kept as it is and described with
  its measurements instead of a tolerance (the exact add-on is MATLAB's bits), on the
  date given;
- **temporary**: an output not produced, or a behaviour not matched, until a planned piece of
  work lands.

Tests name their row in the `reason` of their tolerance (`petakit5d.testing.compare.tol`,
[testing.md](../testing.md)). The area pages: [decon.md](decon.md), [demos.md](demos.md),
[detection.md](detection.md) (Phase 3), [sim.md](sim.md) (structured illumination, Phase 3), [flatfield.md](flatfield.md), [io.md](io.md) (the MATLAB-name layer of the io functions), [matlabcompat.md](matlabcompat.md), [phase2.md](phase2.md) (the
demos' Phase 2 steps), [realdata.md](realdata.md), [stitch.md](stitch.md),
[tools.md](tools.md).

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
| F8 | tiff in the `'libtiff'` and `'imwrite'` modes and uncompressed tiff | cpp-tiff's BigTIFF instead of classic TIFF through `saveastiff`; voxels and tags the same, but for the tags that name the writer (MATLAB's blockedImage writes `Software: MATLAB blockedImage`), which the comparison leaves out | approved 2026-11-11 (plan) | [decisions/0004](../decisions/0004-io-bindings.md) |
| F9 | the libtiff-mode writer's files (`replicated/*_conventional.tif` of the readers/writers demo) | the LZW bytes of each page's last, partial strip differ (about 1 KB a page); voxels and tags identical | approved 2026-11-18 | [demos.md](demos.md), row 2 |

## 2. The default FFT backend (pocketfft) against MATLAB's FFTW

With the exact add-on none of these apply: the results are MATLAB's bits.

| # | Output | Tolerance | Measured | Status | Source |
|---|---|---|---|---|---|
| P1 | FFT-based outputs in general | single: max <= 1e-5 of the largest value, relative L2 <= 1e-6; double: 1e-13 and 1e-14 | single <= 5.5e-7, double <= 1.2e-15 | approved 2026-11-12 | [decisions/0006](../decisions/0006-fft-backends.md) |
| P2 | deconvolution (up to about 30 iterations) and what is made from it (DSR, MIPs, conversions) | relative L2 <= 1e-4; uint16 off by at most 1 in at most 0.1% of voxels | demos: 0.0023-0.0083% of voxels off by 1; z MIPs up to 0.175% | approved 2026-11-12 | [decisions/0006](../decisions/0006-fft-backends.md); [demos.md](demos.md), row 1 |
| P3 | the widefield OMW decon of the widefield/confocal demo (20 iterations, `hann_win_bounds` (0.8, 1.2)), and OMW with many iterations on data where they don't settle, and RL with many iterations (from 26 on the widefield data) | outside P2 | 0.35% of voxels off by 1, 4,586 by more (largest 1,782), relative L2 1.7e-3; MATLAB's own CPU result is as far from the same OMW in exact arithmetic (1.1e-3), MATLAB's GPU result 1.3e-4, the default backend's iterations 9.6e-5; RL at 26 to 30 iterations on the widefield data: off by 1 in 0.103 to 0.134% of voxels (P2's share is 0.1%), never more | approved 2026-11-19 | [demos.md](demos.md), row 8; [decon.md](decon.md), OMW with many iterations; [realdata.md](realdata.md), row 6 |
| P4 | the OMW back projector and its masked Wiener filter | max and relative L2 <= 1e-5; masks identical | | approved 2026-11-13 | tests of `decon` (`PERMISSIVE_BP`) |
| P5 | all-NaN back projectors (degenerate `hann_win_bounds`, `l == u`) | NaN positions compared, not NaN signs | | approved 2026-11-13 | `decon/backprojector.py` |
| P6 | a PSF resampled to one plane (`new_nz = 1`) | difference <= 1e-6 of the raw PSF's largest value | MATLAB's own result there is rounding noise (it takes the spectrum's last plane, not DC) | approved 2026-11-12 | `decon/psf.py`, tests of `decon/test_psf.py` |
| P7 | 3D registration (normalized cross-correlation) | the same integer shift where the correlation is well behaved; peak values within 5e-3 | shifts identical; peaks <= 4.2e-7 (1.8e-7 on the real-data workflow) | approved 2026-11-12 | [decisions/0006](../decisions/0006-fft-backends.md); [stitch.md](stitch.md); [realdata.md](realdata.md), row 3 |
| P8 | chromatic shift estimation: correlation peak values | absolute 1e-4 | <= 1.8e-7; shifts and corrected files exact | approved 2026-11-18 | [tools.md](tools.md), chromatic shift |
| P9 | FFT power spectra: MIPs, central slices, 3D stack | the normalized spectrum `\|F\| / max` within P1's single tolerance | largest difference 4.8e-8 | approved 2026-11-17 | [tools.md](tools.md), FFT power spectra |
| P10 | PSF analysis: OTF slices, line cuts, info files | P1's double tolerance | <= 7.6e-16 of the largest value | approved 2026-11-12 | [tools.md](tools.md), PSF analysis |
| P11 | FSC curves | P1's double tolerance | <= 5.9e-15 of the largest value; radii identical | approved 2026-11-12 | [tools.md](tools.md), Fourier shell correlation |
| P12 | `deconvlucy` (RLdecon's `'original'` method) | P2 | | approved 2026-11-12 | [decon.md](decon.md) |
| P13 | 2D tile registration (`normxcorr2`): correlation peak values | absolute 1e-7 | shifts identical; peaks <= 5.3e-8 | approved 2026-11-18 | [stitch.md](stitch.md) |
| P14 | the FFT spectra and FSC curves of decon results (demo_useful_tools' FFT and FSC steps), end to end | P2's relative L2 (<= 1e-4) of what is made from decon results | spectra relative L2 <= 5.4e-5 (as `\|F\| / max` 1.0e-5, beyond P9's 1e-6); FSC curves <= 2.3e-6 (beyond P11); on MATLAB's own decon results within P9 and P11 | approved 2026-11-18 | [phase2.md](phase2.md), row 2 |

## 3. Tolerances on both backends

MATLAB functions whose last bits could not be reproduced from their documentation and black-box
probes, and MATLAB results that depend on MATLAB's thread count.

| # | Output | Tolerance | Measured | Status | Source |
|---|---|---|---|---|---|
| N1 | single-precision `log10` (the FFT spectra's 3D stack, the figures' values) | within 2 ulp; the port's `log10` is correctly rounded, the same on every platform | 3,692 of 226,981 stack voxels differ, all within 2 ulp | approved 2026-11-17 | [tools.md](tools.md), FFT power spectra |
| N2 | `deconvlucy`'s double results (exact backend) | max <= 1e-13 of the largest value, relative L2 <= 1e-13 | <= 18 ulp (MATLAB's acceleration factor comes from MKL dot products) | approved 2026-11-18 | [decon.md](decon.md) |
| N3 | FSC resolutions | each radius within one 1e-4 step of MATLAB's; resolutions within 1% | all 626 radii per backend equal MATLAB's | approved 2026-11-18 | [tools.md](tools.md), Fourier shell correlation |
| N4 | PSF detection (`radon`, `tand`, `atand` inside) | the decisions (peaks, crops, MIPs, info file) exact | bit for bit in all cases | approved 2026-11-17 | [tools.md](tools.md), PSF detection |
| N5 | filters: `imfilter` with zero padding or `'full'`, small 2-D `imgaussfilt` kernels, general 2-D and 3-D kernels, long-kernel `conv2`, the frequency domain and `edgetaper`, sparse arrays, uint16 along z | per row, 2 to 16 ulp, relative 1e-6 to 2e-6, integers 1 count | 352 of 399 cases bit-identical, the call sites of unmixing, puncta removal and PSF analysis included | approved 2026-11-13 | [matlabcompat.md](matlabcompat.md), filters |
| N6 | puncta removal: pixels filled by `fillmissing2` `'linear'` | 1 ulp, relative 1e-12, in at most 2% of values (uint16: 1 count in at most 1%) | 1 ulp in 26 of 10,647 probe fills | approved 2026-11-18 | [matlabcompat.md](matlabcompat.md), puncta removal |
| N8 | the phase-image flat field (BaSiCPy in place of MATLAB's BaSiC) | flat field: relative 0.06 (0.12 under strong vignetting), relative L2 0.01 (0.03); uint16 tiles 1 count more | about a percent; up to 10.7% at single pixels | approved 2026-11-18 | [flatfield.md](flatfield.md) |
| S1-S10 | point detection's statistics functions (`detection.stats`: `norminv`, `tcdf`, `tinv`, `chi2inv`, `fcdf`, `ksdensity`; MATLAB's closed special functions) | per row, from 4 ulp to `rel` 2e-4 (`ksdensity` from 500 values: 2%) | `ecdf`, `mad`, `kstest2`, the edge values and one degree of freedom bit for bit; where the differences pass a few ulp MATLAB is the less exact side (50-digit checks) | approved 2026-11-20 | [detection.md](detection.md) |
| D1-D13 | point detection (`detection.detect`, `detection.fit`, `detection.sigma`, `detection.run`): `fspecial3('log')`'s sums, `imfilter` of the detector's 3-D kernels, the detector's fits, the `fitGaussian3D` mex's linear and nonlinear modes (MINPACK's Levenberg-Marquardt, where the mex takes GSL's), the PSF sigma estimate, the 'speedup' detector, `frameInfo` with `XR_runDetection3D`'s defaults, the camera search on Windows, the bad-point filter's second valley search, 'lowSNR' without the background correction, the frame-parallel runner's three fixed bugs (D13: without a camera threshold it stops or finds nothing, where XR_runDetection3D takes each frame's noise) | per row, from `rel` 2e-15 (the kernel) to the fits' measured differences | the step-by-step demo's detections MATLAB's in all 60 frames, the multiple demo's in all but 4 of its 180 (4 detections: 3 fits on another path, D3; 1 that MATLAB's detector leaves out later); the master channel within 1e-9 | approved 2026-11-20 (D10-D12 on 2026-11-21; D13 on 2026-11-22 as the port's fix, measured per case with MATLAB's intent) | [detection.md](detection.md) |
| TK1-TK4 | tracking (`detection.tracking`): ties in the assignment solver, `KDTreeBallQuery`'s order at equal distances, the gap-closing costs' last bits, detections at x = 0 | none on the outputs | the demo's `trackedFeatures.mat` and every case bit for bit | approved 2026-11-20 | [detection.md](detection.md) |
| W1 | the condition loading, deskew and channel offsets (`detection.condition`, `detection.wrapper`, `detection.offsets`): `corr`'s Pearson correlations, whose BLAS order MATLAB's CPU and threads decide | `abs` 4e-15 (single `abs` 2e-7) | the other 104 of 110 cases bit for bit; the offsets the correlations decide equal in every case | approved 2026-11-20 | [detection.md](detection.md) |
| W2-W4 | MATLAB's slips the port fixes (`detection.wrapper`, `detection.condition`): reversed scans deskewed as forward ones, the tracking stages' folders shifted after a one-frame movie, the loader failing on a one-frame movie with a frame rate given | the corrected result, not MATLAB's (`matlab_quirks=True` gives MATLAB's for W2 and W3) | bit for bit the corrected MATLAB call (`deskewData` with `'Reverse'`, `XR_runDetTrack3d` on the movies of more than one frame, the load without the rate) | approved 2026-11-20 | [detection.md](detection.md) |
| W11-W20 | track processing (`detection.processing`, `detection.amira`): the fits at gap frames and in the buffers, their `pval_Ar`, the fits that take another path, `catIdx`, the Amira files, the MAT-files as files, options beyond the demos', one-frame branches, the Amira line ends on Windows; and a slip the port fixes, reversed scans' tracks rotated by the forward angle (W20: the mirrored angle, MATLAB's under `matlab_quirks`) | per row, from `rel` 1e-11 and `abs` 1e-10 (inside D3's bound) | every field at the detected frames bit for bit, every `catIdx` (W14: none since 2026-11-21); W13: 2 of the step-by-step demo's 42,412 fits on another path, and 2, 12 and 5 of the multiple demo's 42,412, 71,360 and 58,154 (the second's 13 before the keep rule's bound of 2026-11-22: MATLAB refits a free fit at exactly 2 voxels, and the port since), two of them on the other side of the keep rule (W13's multiple-demo addition approved 2026-11-22) | approved 2026-11-20 (W14 none since 2026-11-21) | [detection.md](detection.md) |
| W5-W6 | MATLAB's failures the port does not copy (`detection.wrapper`): `deskewData` with `aname` left out (its default, the handle `@isstr`) and on a condition without movies | the port's results, not MATLAB's failures (`matlab_quirks=True` gives MATLAB's) | the port's default bit for bit MATLAB's `'aname', []`; MATLAB's identifiers under `matlab_quirks` | approved 2026-11-22 | [detection.md](detection.md) |
| A1 | MATLAB's arguments-block conversions the port does not copy (the MATLAB-name layer of `readtiff`, `writetiff`, `readzarr`, `writezarr` and `createzarr`): numbers, a cell of char arrays or a function handle taken as text, characters taken as logical values | refused with the port's error (`petakit5d:params:conversionRefused`); every other value converted or refused as MATLAB's arguments block does, with MATLAB's identifiers and messages | MATLAB's outcome recorded for each refused value (each one MATLAB converts); the copied values' outcomes MATLAB's, in every form | proposed 2026-11-22 | [io.md](io.md) |
| L1-L7 | the lifetime analysis (`detection.lifetime`: `runLifetimeAnalysis3D`, `plotIntensityCohorts`; `matlabcompat.optim.lsqnonlin`) | per row: the fits `rel` 1e-6 where the port's path is MATLAB's, up to `rel` 0.03 (mode fits) where it is not; the threshold `rel` 1e-4 (0.1 where its fit ends in another local minimum); the cohorts `rel` 5e-3 | on the demo data every `lftRes` field bit for bit but the threshold (9.9e-10); histograms and counts bit for bit; the CCP and CS classification the same in every run but one synthetic condition (2 of about 130 tracks) | approved 2026-11-20 | [detection.md](detection.md) |
| E1 | the point-detection demos end to end (`examples/phase3`): the fits from each side's PSF sigmas (D6) and what is made from them | `rel` 1e-8 with `abs` 1e-8 on the cut-down oracle cases of `tests/e2e/test_demos_phase3.py` | largest there 1.2e-9; on the demo dataset a tail lies beyond it, listed (the step-by-step demo: 10 to 13 positions, 977 amplitudes and 203 p-values of 122,706 values per field; the multiple demo's three experiments: 10 to 152 positions and 127 to 977 amplitudes of 122,706 to 255,100) | approved 2026-11-21 | [phase3.md](phase3.md), row 1 |
| R1-R6 | structured illumination (`petakit5d.sim`): FFT and LAPACK rounding in the double and single steps, the OTF's bead fit, the masked registration, `sim_PSFtoOTF_gen`, and `simRecon`'s files on the CPU and the GPU | per row, from max 1e-13 (double steps) to rel L2 1e-4 (the reconstructions) | about a tenth of each bound; the separation matrix, masks, OTF resampling, `fminsearch`'s traces, the phase deskew and the phase split and merge bit for bit | approved 2026-11-20 | [sim.md](sim.md) |
| N9 | `median`, `mode`, `medfilt2`, `medfilt3` at ties of -0 and +0 | the sign of a zero may differ; values agree | 3 of 269 probe medians | approved 2026-11-12 | [decisions/0005](../decisions/0005-sums-and-math-library.md) |
| N10 | MATLAB's sums, means and standard deviations of 89,000 or more elements, which depend on MATLAB's thread count | the port gives MATLAB's one-thread result at any thread count | up to 1.2e-4 relative between MATLAB's own thread counts (single) | approved 2026-11-12 | [decisions/0005](../decisions/0005-sums-and-math-library.md) |
| N11 | demo_RL_deconvolution's in-memory results (the demo's `sum(psf, 'all')` at 24 threads) | N10: the one-thread sum | the back projector differs by <= 6.0e-8; uint16 results by 1 in <= 0.0082% of voxels; equal to MATLAB's with its one-thread sum | approved 2026-11-12 | [demos.md](demos.md), row 6 |
| N12 | the stitching's global shift assignment (`stitch_assignment_lsqlin`): its linear solves | absolute 1e-9 | sparse solvers (MATLAB's Cholesky and LU, SuperLU here) differ in the last bits | approved 2026-11-18 | tests of `matlabcompat/test_lsq.py` (`ASSIGN_TOL`) |
| N13 | `sph2cart` of NaN inputs | NaN where MATLAB has NaN, NaN payloads not compared | | approved 2026-11-18 | tests of `matlabcompat/test_neighbors.py` |
| N14 | `conv2` and `convn` of small arrays (below about 24 rows and columns, or 16 planes; a 2-D kernel larger than a small image in every dimension excepted), of 3-D kernels larger than the image, of kernels larger than images of 10 rows or 8 planes and more, and of 9 x 9 kernels: the last bits (no call site of the port; which values are Inf, NaN or finite was MATLAB's in every case probed) | not compared (no oracle case; the port's own code calls neither) | MATLAB picks its sums by size: below about 24 rows and columns (or 16 planes) it adds plain products over the overlap with the smaller array's taps in column-major order (3-D: x descending), where the port uses fused chains; the forms of the larger swaps and of 9 x 9 kernels were not identified (measured on 2026-11-19) | approved 2026-11-20 | [platforms.md](../platforms.md), Known problems |

N7 (the puncta removal's background thresholds of batches whose y size is a multiple of 16) was
removed on 2026-11-18: `imgaussfilt3` is bit-identical there since 2026-11-13, and those cases
compare exactly ([matlabcompat.md](matlabcompat.md), puncta removal).

## 4. GPU

The GPU paths use cuFFT 11.0.8, the version MATLAB R2024b bundles, and are MATLAB's bits by
default ([decisions/0003](../decisions/0003-gpu.md)).

| # | Output | Tolerance | Measured | Status | Source |
|---|---|---|---|---|---|
| G1 | GPU deconvolution in the native `(z, y, x)` layout (an option; the default is MATLAB's layout, which is exact) | single: max and relative L2 <= 1e-5 (OTFs 1e-6); uint16 at most 1 count in at most 0.01% of voxels | relative L2 1.0e-6; 0.0016% of voxels | approved 2026-11-12 | tests of `gpu/test_decon_gpu_kernels.py` |

## 5. Not produced yet, and pending fixes

| # | Output | Until | Source |
|---|---|---|---|
| T2 | figures (back projector and OTF masks, FSC, PSF, MIP masks) | the `figures` extra is installed | [demos.md](demos.md), row 5; with the extra, the Phase 2 demos' figures are all drawn ([phase2.md](phase2.md)) |
| T3 | the phase flat field (`ff_corrected_*/`, `_flatfield.mat`) | the `basic` extra is installed | [demos.md](demos.md), row 7; with the extra within N8 ([phase2.md](phase2.md)) |
| T4 | a box with a NaN end (`[... NaN]`) is refused where MATLAB passes it on | the fix of 2026-11-18: accepted where MATLAB accepts it | [cli.md](../cli.md) |
| T5 | the Imaris file of demo_useful_tools (`matlab_decon_omw/DSR_for_useful_tools/imaris/`: the `.ims` file, `parameters.json`, `parameters.mat`, `timepoints.txt`) | the port of `XR_imaris_conversion_data_wrapper`, which waits for the next release of Parallel_Imaris_Writer | [phase2.md](phase2.md) |

T1 (the tile registration and its per-pair results `xcorr/<volume>/xcorr_tile_*.mat`) was
removed on 2026-11-18: the registration is ported (`petakit5d.stitch.registration`), and the
port writes those files as MATLAB does ([stitch.md](stitch.md), [realdata.md](realdata.md)). The
demo runs, made before, were given MATLAB's merged results, so their per-pair files are compared
when the demos run again ([demos.md](demos.md), row 4); those of the phase and 2-photon
stitching are compared in [phase2.md](phase2.md) (P7, P13; every shift identical).

T6 (the point-detection demos' detections and processed tracks beyond D3 and W11-W15) was
removed on 2026-11-22: rows 2, 8 and 9 of [phase3.md](phase3.md) and row 11's categories are
fixed, and what is left of rows 10 and 11 is D3's and W13's fits that take another path (the
multiple demo's 4 detections and two processed tracks on the other side of the keep rule; a third
was the keep rule's bound, the port's since 2026-11-22).
