# Point detection, tracking and lifetime analysis: the scope, traced (task C01)

This page is the specification of Phase 3's detection workstream (tasks C02 to C05): what
MATLAB PetaKit5D runs for 3D point detection, tracking, track processing and lifetime analysis,
stage by stage, with the inputs and outputs of each stage, the black-box references recorded for
it, and what each task implements.

**Provenance.** The point-detection and u-track code is GPL code of other authors, or the lab's
code derived from it. It was not read. Everything here comes from running MATLAB as a black
box:

- MATLAB's profiler: the functions that ran, their files, call counts, times, which callers
  called them, and the line *numbers* executed (no line text);
- logging wrappers put in front of the functions on MATLAB's path: the class and size of each
  argument and result, the values of small ones (option names, scalars), and the callers;
- the stages' input and output files, and oracle cases (below).

The published methods are the specification: Aguet et al. 2013, Dev. Cell 26:279 (model-based
detection with statistical tests), and Jaqaman et al. 2008, Nat. Methods 5:695 (u-track: linking
and gap closing as linear assignment problems with Kalman motion models). MathWorks functions are
reproduced from their documentation and black-box probes (`petakit5d.detection.stats`,
`petakit5d.matlabcompat`).

## 1. How it was traced

MATLAB R2024b from the oracle tree (commit 45b2bb1), on a copy of the demos' dataset (one
condition of three experiments, each two channels of 60 frames of 512 x 320 x 35 raw voxels,
0.104 um pixels, 0.5 um sample scan, 31.5 degree skew, 2.8 s frame interval; and the two
channels' bead PSFs), one fresh MATLAB session per run, no parallel pool (parfor bodies run in
the session):

| Run | What | Frames | MATLAB time |
|---|---|---|---|
| by step | `demo_point_detection_3d_by_step` (each stage called in turn) | 1 experiment, 2 x 60 | 531 s |
| integrated | `demo_point_detection_3d_integrated` (`XR_runDetTrack3d` with the demo's options) | 1 experiment, 2 x 60 | 722 s (with wrappers) |
| multiple experiments | `demo_point_detection_3d_multiple_experiments` (`XR_runDetTrack3d`, lifetime analysis, intensity cohorts) | 3 experiments, 2 x 60 each | 2,577 s (with wrappers; lifetime analysis 34 s, cohorts 1 s) |
| defaults | `XR_runDetTrack3d(data, 'Sigma', sigma)`, every other option its default | 1 experiment, 2 x 60 | 565 s (with wrappers) |
| test | `PointDetectionTest` (tests/point_detection; 6 of 6 pass) | synthetic | 11 s |

The demos load their data interactively (`XR_loadConditionData3D(rt)` asks for the channel
folders); the runs call `XR_loadConditionData3D(rt, {'ch1', 'ch2'}, {'gfp', 'rfp'})`, its
non-interactive form, which gives the same answers. Runs with the logging wrappers and without
gave identical output files (45 of 45 files of a 6-frame run, and 251 of 251 of the full by-step
run, deskewed frames aside; `.mat` contents decoded).

Time is spent in track processing (216 s of the by-step run), detection (183 s, 1.5 s per
frame and channel), the channel offset estimate (59 s), tracking (37 s), track rotation (19 s)
and deskew (15 s).

## 2. The pipeline

`XR_runDetTrack3d(data, ...)` runs the stages in this order; the by-step demo calls them one by
one. `<ch>` is a channel folder, `<results>` the analysis folder under the first channel
(`'aname'`, default `Analysis`).

| # | Stage | MATLAB entry | Reads | Writes | Task |
|---|---|---|---|---|---|
| 1 | Load the condition | `XR_loadConditionData3D(root, chNames, markers)` | folder and frame names, the first frame's size, `*Settings.txt` | the `data` struct (no file) | C04 |
| 2 | PSF sigma estimate | `GU_estimateSigma3D(psfFile, [])` | a bead PSF tiff | `[sigmaXY, sigmaZ, XYZ]` (no file) | C02 |
| 3 | Deskew | `deskewData(data, ...)` | raw frames | `<ch>/DS/<frame>.tif` | C04 (D02's deskew kernel) |
| 4 | Channel offsets | `XR_correctXZoffsetData3D(data, 'Overwrite', ...)` | DS frames | `<ch>/DS/xz_offset.mat`, DS frames of the other channels shifted | C04 |
| 5 | Detection | `XR_runDetection3D` (by step) or `XR_runDetection3D_frame_parallel` (`XR_runDetTrack3d`) | DS frames | `<results>/Detection3D.mat`, `<results>/Masks/dmask_NN.tif`, `<results>/dynRange.mat`, `<results>/settings.mat` | C02 |
| 6 | Tracking | `XR_runTracking3D(data, loadTrackSettings(...), ...)` | `Detection3D.mat` | `<results>/trackedFeatures.mat` | C03 |
| 7 | Track processing | `XR_runTrackProcessing3D(data, ...)` | `trackedFeatures.mat`, `Detection3D.mat`, DS frames | `<results>/ProcessedTracks.mat`, `<results>/Amira/Tracking/tracking_t_NNNN.am` | C04 |
| 8 | Track rotation | `rotateTracks3D(data, ...)` | `ProcessedTracks.mat` | `<results>/RotatedTracks.mat` | C04 |
| 9 | Lifetime analysis | `runLifetimeAnalysis3D(data, 'AnalysisPath', aname, 'overwrite', true)`, `plotIntensityCohorts(data, 'AnalysisPath', aname)` | every movie's `ProcessedTracks.mat` | `<results>/lifetimeData.mat` per movie; the `lftRes` and cohort structs (returned) | C05 |

### Options as the demos pass them (recorded by the wrappers)

- `XR_runDetTrack3d` with the demos: `'Sigma'` (2 x 2, `[sigmaXY, sigmaZ/zAniso]` per channel),
  `'DetectionMethod', 'lowSNR'`, `'aname'`, `'FitMixtures', false`, `'BackgroundCorrection',
  true`, and `'XZoffsetCorrection', true` (integrated demo). The multiple-experiments demo passes
  `'zoffsetCorrection', true`, which `XR_runDetTrack3d` does not know and drops (it keeps
  unmatched options): that demo runs without the offset correction. Its defaults (the
  `defaults` run): `DetectionMethod 'lowSNR'`, `BackgroundCorrection false`, `Track true`,
  `RotateTracks true`, `TrackingRadius [3 6]`, `TrackingGapLength 2`, `Buffer [3 3]`,
  `ParallelMethod 'lock_distributed'`, `Mode 'xyzAc'`, `MaxMixtures 2`.
- Deskew: `'Overwrite', 'SkewAngle', 31.5, 'aname', 'Rotate', false, 'sCMOSCameraFlip', false,
  'Crop', false`, the light-sheet flat-field options off, `'zoffsetCorrection', false`.
- Detection, per frame and channel (from `XR_runDetection3D_frame_parallel`):
  `XR_estimate_overall_background_information(vol, 'randseed', 1)` (3 values: two background
  levels and a noise level), `XR_get_camera_standard_deviation_threshold(data)` (3 for the
  camera type `'SCMOS'` that `get_camera_type_from_setting_file` finds), then
  `XR_pointSourceDetection3D_3(vol, sigma, 'Alpha', 0.05, 'Mask', [], 'RemoveRedundant', true,
  'RedundancyRadius', 0.5, 'RefineMaskLoG', false, 'WindowSize', [], 'Mode', 'xyzAc',
  'FitMixtures', false, 'MaxMixtures', 2, 'BackgroundCorrection', true,
  'FilterByResidualSigma', true, 'ResidualSigmaThresh', 3, 'BackgroundInfo', [b1 b2 b3 3])`.
- Sigma estimate: `XR_pointSourceDetection3D(psf, s, 'Mode', 'xyzAsrc', 'Alpha', 0.05, 'Mask',
  [], 'RemoveRedundant', true, 'RefineMaskLoG', false, 'WindowSize', [], 'FitGaussianMethod',
  'original', 'BackgroundCorrection', false)`, three times, starting from `s = [1.2 2.5]` and
  each time from the previous result.
- Tracking: `loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2)` gives two cost matrices
  (`costMatLinearMotionLink2`, `costMatLinearMotionCloseGaps2`), `gapCloseParam` with
  `timeWindow 3, mergeSplit 1, minTrackLen 1`, link parameters `linearMotion 0,
  minSearchRadius 3, maxSearchRadius 6, brownStdMult 3`, and the Kalman functions
  (`kalmanResMemLM`, `kalmanInitLinearMotion`, `kalmanGainLinearMotion`,
  `kalmanReverseLinearMotion`). Merging and splitting are on.
- Track processing: `'TrackerOutput', 'trackedFeatures.mat', 'FileName',
  'ProcessedTracks.mat', 'Buffer', [3 3], 'BufferAll', false, 'FitMixtures', false,
  'WindowSize', []`, and `'ParallelMethod', 'lock_distributed'` from `XR_runDetTrack3d`.

## 3. The stages

Sizes are those of the demo data: DS frames are 512 x 460 x 35 (y, x, z) doubles in detection;
about a thousand detections per frame and channel; a few thousand tracks per experiment, a
sixth of them valid.

### 3.1 Load the condition (C04)

`XR_loadConditionData3D` (lab-edited, help text read up to its author line) lists the movie
folders (`Ex*`, `*_z0p5` gives `dz`), the channel folders, the frames (sorted by their frame
number; the `msec` token gives the frame rate) and calls `XR_loadConditionData`. Output per
movie: `source` (the first channel's folder), `channels`, `date`, `framerate` (s),
`imagesize` (`[ny nx nz]` of the raw frames), `movieLength`, `markers`, `framePaths`,
`maskPaths`, `pixelSize`, `dz`, `zAniso` (`dz * sind(angle) / pixelSize`),
`angle` (31.5), `objectiveScan`, `reversed`, `framePathsDS`, `framePathsDSR`.
Quirk: `date` is taken from digits in the path (a 6-digit token), so it changes with the folder
the data is in; compare it masked.

### 3.2 PSF sigma estimate (C02)

`GU_estimateSigma3D(file, [])` reads the PSF (128 x 128 x 101 in the demo), fits it three times
with `XR_pointSourceDetection3D` in mode `'xyzAsrc'` (sigma free), and returns `sigmaXY` and
`sigmaZ` and the bead's position. The demos divide `sigmaZ` by `zAniso`.

### 3.3 Deskew and channel offsets (C04)

`deskewData` (GPL, not read) runs, per frame in a parfor loop, `readtiff`,
`deskewFrame3D` (the copy in `llsmtools/transforms`; D02 has the port of its behaviour: an
affine `imwarp`) and `writetiff` into `<ch>/DS/`. 512 x 320 x 35 raw frames become
512 x 460 x 35. `XR_runDetTrack3d` passes `'Overwrite', false` (the traced runs started without
DS frames, so what it does with existing ones was not traced).

`XR_correctXZoffsetData3D` (lab) calls `XZoffsetEstimation3D(ds1, ds2)` on three frames: Pearson
correlations (`corr`, 165 per frame) of the two channels' DS volumes over shifts in x and z,
then `imregionalmax`; the result (integer offsets in x and z) goes to
`<ch>/DS/xz_offset.mat`, and the other channels' DS frames are rewritten shifted.

### 3.4 Detection (C02)

Per frame and channel, `XR_pointSourceDetection3D_3` (method `'lowSNR'`;
`XR_pointSourceDetection3D` for `'speedup'` and the sigma estimate) runs:

1. `estGaussianAmplitude3D` (llsmtools copy): two separable filterings with `conv3fast`
   (imageProcessing's mex; per axis a kernel given by 4 or 5 coefficients from 1 down): the
   Gaussian amplitude and background estimates;
2. filtering with `fspecial3` kernels and `imfilter` (N04's `matlabcompat.filters` covers the
   form), `nanmedian` of the volume (`'all'` and per z plane: `prctile(x, 50, ...)`);
3. local maxima: `XR_locmax3d` (B08's `tools.puncta_removal.locmax3d`) with `imdilate`, then
   connected components (`bwconncomp`, `labelmatrix`, `sub2ind`/`ind2sub`);
4. thresholds: `chi2inv(0.9, n)` with `n` the window's voxel count (343 for 7 x 7 x 7),
   `norminv(0.975)`, and `tinv(alpha, df)` cached per alpha (below);
5. `XR_fitGaussians3D` fits every candidate: `XR_fitGaussian3D` (an M-file fit, about 11
   iterations of `compute_gradient_hessian_mex` per candidate) in mode `'xyzAc'`
   (`XR_fitGaussian3D_Ac`), the `fitGaussian3D` mex (no source) in the sigma estimate's
   `'xyzAsrc'` fits; residual tests: `adtest_mex(residuals, 0, sigma, 0.05, 2)`
   (Anderson-Darling), and a Welch t-test of the amplitude against the background
   (`tcdf(t, df)` on vectors, non-integer df of 300 to 800); output fields `pval_Ar`,
   `hval_Ar`, `hval_AD`, `isPSF`;
6. `filter_bad_points_in_detection(pos, A, 3, bgInfo, vol, ...)`: neighbours by `knnsearch`
   (k 3) and `rangesearch`, `histcounts` (60 bins), `findpeaks`, `medfilt1`, `prctile`,
   `median`, `std`, `conv`;
7. `KDTreeBallQuery` (CMEAnalysis' mex, no source) removes detections within 0.5 voxel of
   another (`RemoveRedundant`).

Result per frame: a struct with `x, y, z, A, s (2 x n), c, x_pstd, y_pstd, z_pstd, A_pstd,
s_pstd, c_pstd, x_init, y_init, z_init, sigma_r, SE_sigma_r, RSS, pval_Ar, hval_Ar, hval_AD,
isPSF` and the mask (logical, the frame's size). `XR_runDetection3D*` writes them, with both
channels in each field's rows and the tracker's coordinates (`zCoord` = z times `zAniso`), to
`Detection3D.mat` (`frameInfo`), and the masks to `Masks/dmask_NN.tif`; the Amira writer of the
detections fails with a warning without identifier (third-party code: its text is not kept
here) and leaves an empty `Amira/Detection` folder (MATLAB's behaviour at 45b2bb1).

**Quirks for C02:**

- the t threshold is computed once per MATLAB session and alpha and reused (`tinv` ran twice in
  a whole run: once in the sigma estimate, once in detection), so the second channel's
  threshold is the first channel's, whatever its window and degrees of freedom
  (PointDetectionTest clears the function for this reason). On the synthetic references a second
  call gave the same detections as a fresh one (case `psd_threshold_cache`);
- `XR_estimate_overall_background_information` samples the volume with MATLAB's random number
  generator seeded with 1 (`'randseed', 1`), so its values depend on MATLAB's Mersenne Twister
  stream and its sampling functions;
- `tcdf` of the large negative t statistics gives values far below `realmin` (MATLAB keeps
  subnormal results; SciPy rounds them to 0): compare `pval_Ar` with an absolute floor.

### 3.5 Tracking (C03)

`XR_runTracking3D` loads `Detection3D.mat` and calls `trackCloseGapsKalmanSparse`:

- frame-to-frame linking: `linkFeaturesKalmanSparse` three times (forward, backward, forward),
  `costMatLinearMotionLink2` per frame pair (`pdist2`, `KDTreeBallQuery`), the Kalman filter
  (`kalmanInitLinearMotion` per new feature, `kalmanGainLinearMotion` per frame,
  `kalmanResMemLM`, `kalmanReverseLinearMotion`), and `lap` with `mexLap` (a sparse linear
  assignment solver, no source; `lap(cost, -5, 0)`: cost matrices of up to 2,000 x 2,000 with
  birth and death blocks);
- gap closing with merging and splitting: `costMatLinearMotionCloseGaps2` (450 lines run;
  `createDistanceMatrix` mex, `nanmean`, `mean`, `prctile`, `mpower`), `getAveDispEllipseAll2`
  (search ellipses), `estimTrackTypeParamLM2` with `asymDeterm2D3D` (track type from the
  asymmetry of its positions), `getTrackSEL`, then `lap` once on the gap-closing problem (tens
  of thousands of track starts and ends per movie);
- output `trackedFeatures.mat`: `tracksFinal` (per track `tracksFeatIndxCG`,
  `tracksCoordAmpCG` with 8 values per frame, `seqOfEvents`), as PointDetectionTest reads it.

### 3.6 Track processing and rotation (C04)

`XR_runTrackProcessing3D` (lab, 593 lines run and its local functions `interpTrack`,
`mergeStructs`, `roundConstr`) converts the tracker's output, classifies tracks (valid,
gaps, merges and splits, touching the movie's ends), interpolates gaps (`interp1`), reads out
3 frames before and after each track ("buffer") with Gaussian fits at fixed positions
(`estGaussianAmplitude3D`, `XR_locmax3d`, the `fitGaussian3D` mex 34,000 times,
`compute_gradient_hessian_mex`, `adtest_mex`, `tcdf`), computes per-track statistics
(`nanmedian`, `nansum`, `std`, `mean`, `prctile`, `norminv`), and writes `ProcessedTracks.mat`
(struct array `tracks` with `t, f, x, y, z, A, c, *_pstd, sigma_r, SE_sigma_r, pval_Ar, isPSF,
tracksFeatIndxCG, gapVect, gapStatus, ...`) and one Amira file per frame
(`XR_amiraWriteTracks`, `dlmwrite`, in a parfor loop). `rotateTracks3D` writes
`RotatedTracks.mat`. On a 6-frame cut of the data every track touches the movie's ends, and
none is valid.

## 4. Which copy of a duplicated name runs

51 function names have files in more than one folder of the path (112 files: `.m`, `.mexa64`
and their MATLAB fall-backs, listed by name only). 25 of them are reached; the profiler's file
for each (`CME/` is `pointDetection/CMEAnalysis_104_XR/software/`, `llsm/` is
`pointDetection/llsmtools/`):

| Name | The copy that ran | Not run |
|---|---|---|
| `KDTreeBallQuery` | `CME/mex/KDTreeBallQuery.mexa64` | its `.m`; `llsm/mathfunc` mex |
| `binterp` | `imageProcessing/interpolation/binterp.mexa64` | its `.m`; `CME/mex` m and mex |
| `conv3fast` | `imageProcessing/conv3fast.mexa64` | its `.m`; `llsm/psdetect3d` m and mex |
| `createDistanceMatrix` | `pointDetection/CMEAnalysis_104_XR/createDistanceMatrix.mexa64` | `CME/mex` m and mex; `pointDetection/mex/mexCreateDistanceMatrix` m and mex |
| `deskewData` | `pointDetection/deskewData.m` | `llsm/cmeAnalysis3D` |
| `deskewFrame3D` | `llsm/transforms/deskewFrame3D.m` | `pointDetection/src` |
| `estGaussianAmplitude3D` | `llsm/psdetect3d/estGaussianAmplitude3D.m` | `pointDetection/src` |
| `fitGaussian3D` | `llsm/psdetect3d/fitGaussian3D.mexa64` | its `.m`; `pointDetection/mex` mex |
| `getCellDir`, `getDirFromPath`, `getMovieName`, `getShortPath`, `recursiveDir`, `sortStringsByToken` | `CME/` | `llsm/pathfunc` |
| `getFluorPropStruct`, `getLifetimeData`, `getVisitorIndex`, `loadTrackSettings`, `loadTracks`, `plotIntensityCohorts`, `plotLifetimes` | `CME/` | `llsm/cmeAnalysis3D` |
| `loadFigureSettings`, `setupFigure` | `CME/` | `llsm/graphics` |
| `scaleEDFs` | `CME/` | `llsm/mathfunc` |
| `interpln` | `llsm/mathfunc/interpln.m` | `pointDetection/src` |

The other 26 duplicated names (`pointSourceDetection3D`, `locmax3d`, `runDetection3D`,
`fitGaussians3D`, `scaleContrast`, `filterGauss*`, ...) are not reached by any run: the `XR_`
versions in `pointDetection/src` run instead.

## 5. What C02 to C05 implement

Lines run are the profiler's executed line numbers per file, over all runs (no text read):
6,230 lines in 75 files, plus 11 mex files.

| Task | MATLAB reached (files, lines run) | Mex without source | MathWorks functions it needs (port status) |
|---|---|---|---|
| C02 detection | `XR_runDetection3D` (198), `XR_runDetection3D_frame_parallel` (246), `XR_pointSourceDetection3D_3` (195), `XR_pointSourceDetection3D` (120), `XR_fitGaussians3D` (143), `XR_fitGaussian3D` (257), `XR_fitGaussian3D_Ac` (49), `filter_bad_points_in_detection` (219), `XR_estimate_overall_background_information` (51), `estGaussianAmplitude3D` (23), `XR_locmax3d` (25; B08 has it: `tools.puncta_removal.locmax3d`), `GU_estimateSigma3D` (27), `XR_get_camera_standard_deviation_threshold` (11), `get_camera_type_from_setting_file` (40): about 1,600 lines | `fitGaussian3D` (34,000 calls in track processing, 31 in the sigma estimate), `KDTreeBallQuery`, `conv3fast` (imageProcessing); `adtest_mex` and `compute_gradient_hessian_mex` have C sources in `mex_functions` (their headers decide), `replace_nan_with_value` is the lab's (K04) | `imfilter` (N04, `matlabcompat.filters`), `fspecial3` (not yet: N04 has `fspecial`), `imdilate`/`strel` sphere (`matlabcompat.morphology`), `bwconncomp`, `labelmatrix`, `bwlabeln` (`matlabcompat.label`), `prctile`, `median`, `mean`, `std` (`matlabcompat`), `chi2inv`, `norminv`, `tinv`, `tcdf` (`detection.stats`), `knnsearch` (`matlabcompat.neighbors`), `rangesearch` (not yet), `histcounts`, `conv`, `ismember`, `findpeaks`, `medfilt1` (Signal Processing Toolbox, in `filter_bad_points_in_detection`), `meshgrid`, MATLAB's `rng` and random sampling (open) |
| C03 tracking | `XR_runTracking3D` (63), `loadTrackSettings` (58), `trackCloseGapsKalmanSparse` (366), `linkFeaturesKalmanSparse` (186), `costMatLinearMotionLink2` (104), `costMatLinearMotionCloseGaps2` (450), `kalmanGainLinearMotion` (65), `kalmanInitLinearMotion` (46), `kalmanResMemLM` (8), `kalmanReverseLinearMotion` (5), `lap` (46), `getAveDispEllipseAll2` (202), `estimTrackTypeParamLM2` (53), `asymDeterm2D3D` (26), `getTrackSEL` (30), `coordAmpMatFromIndicesSparse` (15), `progressText`, `sec2struct`: about 1,800 lines | `mexLap` (the assignment solver), `KDTreeBallQuery`, `createDistanceMatrix` | `pdist2`, `knnsearch`, `nanmean`, `mean`, `var`, `nancov`, `prctile`, `sortrows`, `unique`, `setdiff`, `intersect`, `spdiags`, sparse matrices |
| C04 wrapper, loading, processing | `XR_runDetTrack3d` (94), `XR_loadConditionData3D` (101), `XR_loadConditionData` (128), `recursiveDir`, `sortStringsByToken`, `getDirFromPath`, `getShortPath`, `deskewData` (119), `deskewFrame3D` (26, D02), `XR_correctXZoffsetData3D` (82), `XZoffsetEstimation3D` (45), `XR_runTrackProcessing3D` (593), `binarySegmentLengths`, `locmax1d`, `interpln`, `getFluorPropStruct`, `rotateTracks3D` (70), `XR_amiraWriteTracks` (26): about 1,500 lines | `fitGaussian3D` (track processing's fits) | `corr`, `imregionalmax` (B08: `tools.puncta_removal`), `imwarp` (D02), `interp1`, `conv`, `intersect`, `ismember`, `nanmedian`, `nansum`, `std`, `prctile`, `tcdf`, `norminv`, `dlmwrite` |
| C05 lifetime analysis | `runLifetimeAnalysis3D` (250), `getLifetimeData` (108), `loadTracks` (54), `scaleEDFs` (136), `fitGaussianModeToCDF` (23), `adtest1` (40), `getVisitorIndex` (15), `getCellDir`, `getMovieName`, `getParentDir`, `isposint`, and `plotIntensityCohorts` (286); the figures: `plotLifetimes` (91), `setupFigure` (109), `loadFigureSettings` (36), `formatTickLabels` (49), `rotateXTickLabels` (44), `getFluorophoreHues`, `name2wavelength`, `wavelength2rgb`: about 700 lines of numbers and 400 of figures | `binterp` (imageProcessing, 1,948 calls) | `lsqnonlin` (trust-region-reflective with bounds, 69 calls), `ecdf` (73 calls; `detection.stats`, exact), `hist` with given centers, `interp1` (`'linear'`), `normcdf`, `prctile`, `nanmean`, `nanmedian`, `nanstd` |

## 6. Lifetime analysis (C05)

`runLifetimeAnalysis3D(data, 'AnalysisPath', aname, 'overwrite', true)` over the three
experiments: `loadTracks` reads every movie's `ProcessedTracks.mat` (24 s of the 34 s),
`getLifetimeData` (a parfor over movies) writes `<results>/lifetimeData.mat` per movie, then:

- intensity thresholds from the first frames (the run prints the maximum intensity threshold
  on the first 2 frames and the 95th percentile of the first frame's distribution), `prctile`,
  `normcdf`;
- per movie and category, empirical distributions with `ecdf` (column vectors of 5 to 910
  values) interpolated on 999 points (`interp1`, `'linear'`), aligned between movies by
  `scaleEDFs` (`lsqnonlin` on `[scale shift]` with bounds `[0 -1]` to `[Inf 1]`, `TolFun` and
  `TolX` 1e-6), and fitted by `fitGaussianModeToCDF` (`lsqnonlin` on three parameters, bounds
  `[0 0 0]` to `[1 Inf Inf]`, 1e-8);
- `adtest1` (CMEAnalysis' Anderson-Darling test, 10 calls) and `binterp`;
- lifetime histograms with `hist` on 48 centers (the frame times).

It returns `lftRes`: `t` (1 x 48), `cellArea`, `lftHist_Ia`, `lftHist_Ib`, `lftHist_IIa`,
`nSamples_Ia`, `MaxIntensityThreshold`, `pctCCP`, `pctCS`, `pctVisit`, `nCCP`, `lftHistCCP`,
`lftHistCS`, `lftHistCCP_counts`, `lftHistCS_counts`, `histWeight`, `lftHistVisit`,
`slaveCombs`, `lftHistSlaveAll`, `lftHistSlaveCCP`, `lftHistSlaveCS` and their `_counts`,
`pctSlaveCCP`, `pctSlaveCS`, `startsPerFrameAll`, `initDensityAll`, `initDensityIa`,
`initDensityCCP` (NaN without cell masks), `meanLftHistCCP`, `meanLftHistCS`,
`meanLftHistVisit`; and draws figures (`plotLifetimes`).
`plotIntensityCohorts(data, 'AnalysisPath', aname)` returns, per channel, `bounds` (1 x 7
cohort bounds in s), `t`, `Amin`, `Aplus` and `A` (2 x 6 cells), with `ecdf` and `lsqnonlin`,
and draws its figure. The numbers are C05's; the figures are optional.

On a small synthetic condition (three experiments of 30 frames, 22 tracks each), MATLAB's
lifetime analysis stops with `stats:ecdf:NotEnoughData` ("Sample data contains only missing
data"): a category without tracks (case `multiple.lifetime_empty_category`). The lifetime
reference is the demo data's run (section 8).

## 7. Statistics functions (`petakit5d.detection.stats`, this task)

`norminv`, `tcdf`, `tinv`, `chi2inv`, `fcdf`, `ksdensity`, `ecdf`, `mad` and `kstest2` are
reproduced from MATLAB's documentation and black-box probes (no MathWorks source read), with
SciPy's special functions where MATLAB uses its closed special functions; so their last bits
differ. Measured against MATLAB R2024b (`tests/detection/test_stats.py`, 155 oracle cases,
about 600,000 values: probabilities in 1e-4 steps, decades to 1e-6 and their complements, the
tails to `realmin` and the subnormals; t values in 0.01 steps over [-40, 40] and their tails to
1e300; degrees of freedom from 0.5 to 1e10 and Inf, non-integer ones of the Welch tests; edge
values 0, 1, outside [0, 1], NaN and Inf; samples of 37 to 1,000 values with ties and NaN):

- **Exact** (bit for bit on every value): `ecdf` (MATLAB computes the Kaplan-Meier form
  `1 - cumprod(1 - d./r)`, probed: not `cumsum(d)/n`), `mad` (both flags, vectors and
  columns, NaN), `kstest2` (all three tails, `'Alpha'`; its series summed in MATLAB's order),
  `ksdensity`'s evaluation points and bandwidth; the edge values (NaN, Inf, 0, 1, invalid
  parameters) of the distribution functions; `tcdf(x, 1)` (`atan2(1, -x)/pi` below 0,
  `1 - atan2(1, x)/pi` above, fdlibm `atan2`) and `tinv(p, 1)` (`tan(pi*(p - 0.5))`, fdlibm
  `tan`): found by probing.
- **Within a tolerance**, proposed (rows S1 to S10; the oracle cases name them):

| Row | Function and range | Measured (largest) | Proposed | Cause |
|---|---|---|---|---|
| S1 | `norminv`, `p` in [realmin, 1], `tinv(p, Inf)` | 4 ulp; `mu + sigma*z` near 0: 2e-15 of the largest value | 4 ulp or 2e-15 of the largest value | MATLAB's `erfcinv` against SciPy's `ndtri` |
| S2 | `norminv` of subnormal `p` | 1.3e-4 relative | 2e-4 relative | MATLAB loses digits below `realmin`; SciPy keeps them |
| S3 | `tcdf`, `nu` up to 1e4, above 1e7 (normal cdf) | 4.5e-12 relative (`nu` 555.75); values below 1e-300 | 1e-11 relative, 1e-300 absolute | MATLAB's `betainc`; MATLAB keeps subnormal tails that SciPy rounds to 0 |
| S4 | `tcdf`, `nu` 1e6 to 1.5e7 | 1.3e-7 relative | 2e-7 relative | MATLAB's form for large `nu` (SciPy is exact there, checked at 50 digits) |
| S5 | `tinv`, `p` in [1e-6, 1 - 1e-6] | 1.2e-10 relative (`nu` 10, `p` 1e-6); 1.2e-12 at `nu` 100 to 1,000 | 2e-10 relative, 1e-14 absolute | MATLAB's inverse is the less exact one (SciPy's inversion refined by Newton steps on the cdf: within 2 ulp of the 50-digit value) |
| S6 | `chi2inv`, `p` in [1e-6, 1] | 7.3e-12 relative (`nu` 1e6), 1.7e-15 at `nu` 343 | 1e-11 relative | MATLAB's `gammaincinv` |
| S7 | `fcdf` | 8.4e-12 relative | 1e-11 relative, 1e-300 absolute | MATLAB's `betainc` |
| S8 | `fcdf`, `nu1` 1e6 | 7.9e-9 relative | 1e-8 relative | MATLAB's form for a large `nu1` (SciPy exact) |
| S9 | `ksdensity`'s density, samples of 100 and 200 | 10 ulp | 16 ulp | the kernel sum's order |
| S10 | `ksdensity`'s density, 500 values | 1.5% (far tails), 0.01% in the middle | 2% relative | MATLAB evaluates larger samples otherwise (not probed further: not called by the traced code) |

- **Not compared** (documented, the port keeps the exact values): `tinv` for `p` below 1e-6
  (MATLAB returns -Inf below about 2^-54 for fewer than about 600 degrees of freedom, and is off
  by up to 23% at 1e-16 and 0.2% up to 1e-12), `chi2inv` below 1e-6 (NaN below 1e-150 for
  hundreds of degrees of freedom, 4% off at subnormal `p`). The detection code uses `p` 0.05,
  0.8, 0.9, 0.975 and the Welch tests' degrees of freedom, all inside the compared ranges.

The traced calls: `chi2inv(0.9, n)` (window voxel counts), `norminv(0.975, 0, 1)` (equal to
MATLAB's bits), `tinv(0.05, df)` and `tinv(0.8, df)` (once per session), `tcdf(t, df)` with
vectors of 1 to 1,300 values (detection) and scalars (track processing), `ecdf(y)` with column
vectors (lifetime analysis, cohorts), `normcdf` inside `tcdf` and in the lifetime analysis,
`gaminv` only inside `chi2inv`. `ksdensity`, `mad`, `kstest2` and `fcdf` are not called by any
traced run; they are implemented and checked in their documented default forms.

## 8. Behaviour references

**Oracle cases on synthetic data** (MATLAB-only cases, `test_references.py`; the module is
kept with the task's notes until C02 to C05 take it, see the task report): a synthetic
condition (one experiment, two channels, 16 frames of 16 x 56 x 64 voxels; 10 tracks of
Gaussian spots with sigma (z, yx) (1.6, 1.3) and (1.9, 1.5) on a background of 100 with
Gaussian noise of 5, amplitudes 210 to 380, tracks that start and end inside the movie and one
with a one-frame gap; the second channel shifted by 0.4 and -0.3 voxels at 0.7 times the
amplitude; each frame written as raw and as DS frame, so deskew is skipped; a settings file
naming the camera; two bead PSFs), and a three-experiment condition for the lifetime analysis.
Each case clears the detection functions first (persistent caches). All 17 pass in both modes
(store, and manifest only).

| Case | Stage(s) | Outputs | For |
|---|---|---|---|
| `load_condition` | 1 | the `data` struct | C04 |
| `estimate_sigma` | 2 | `sigmaXY`, `sigmaZ`, position, both PSFs | C02 |
| `camera_and_background` | 5 (set-up) | camera type (`'SCMOS'`), threshold (3), background information | C02 |
| `track_settings` | 6 (set-up) | `loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2)` and its defaults | C03 |
| `psd_lowsnr.ch{1,2}.frame{0,8}` | 5 | background information, detection struct, mask | C02 |
| `psd_threshold_cache` | 5 | ch1 then ch2 in one session, ch2 fresh, and whether they agree (they do here) | C02 |
| `psd_speedup_form.ch{1,2}` | 5 (`'speedup'`, PointDetectionTest's form) | detection struct, mask | C02 |
| `psd_sigma_estimate_form.ch{1,2}` | 2 (one fit) | detection struct, mask | C02 |
| `chain.by_step` | 1, 2, 5 to 8 (by-step demo's calls) | `data`, sigma; `Detection3D.mat`, masks, `trackedFeatures.mat`, `ProcessedTracks.mat`, `RotatedTracks.mat`, Amira files | C02 to C04 |
| `chain.integrated` | 1, 2, 5 to 8 (`XR_runDetTrack3d`'s calls) | the same files | C02 to C04 |
| `chain.speedup_no_background` | 1, 5, 6 (PointDetectionTest's form) | `Detection3D.mat`, masks, `trackedFeatures.mat` | C02, C03 |
| `multiple.lifetime_empty_category` | 1 to 9 (the multiple-experiments demo's calls) | MATLAB's error (section 6) | C04, C05 |

Every stage's input is in the store with its output: a later stage's input is the earlier
stage's output file of the same case, so C03 can feed MATLAB's `Detection3D.mat` to its tracker
and C04 MATLAB's `trackedFeatures.mat` to its track processing.

**Demo-data runs** (cluster only, outside the repository): the five runs of section 1 keep
every output file (DS frames, `Detection3D.mat`, masks, `trackedFeatures.mat`,
`ProcessedTracks.mat`, `RotatedTracks.mat`, `lifetimeData.mat`, the `lftRes` and cohort
structs) next to the traces, for C02 to C05's comparisons on the demo data ("the detected set
equal, positions, amplitudes and backgrounds within a proposed tolerance", "lifetime
distributions within a proposed tolerance").

## 9. Out of scope

- The GUI branches and plotting (the figures of `plotIntensityCohorts` and `plotLifetimes`,
  `setupFigure` and the label helpers) except the numbers they compute;
- the Amira writers (`XR_amiraWriteDetections` fails in MATLAB and writes nothing;
  `XR_amiraWriteTracks` writes visualization files): optional, after the numbers;
- detection methods and options the demos do not reach: `'original'`, `'adaptive'`,
  `FitMixtures` (`XR_fitGaussianMixtures3D`, the `fitGaussianMixture3D` mex), `RefineMaskLoG`,
  `objectiveScan` data, the light-sheet flat-field options of `deskewData`, `'Rotate'` and
  `'Crop'`;
- `XR_correctZoffsetData3D` and the commented-out `GU_runTrackProcessing3D`;
- the `'lock_distributed'` mechanics (lock files of `chunk_lock_clean`) beyond what the job
  framework gives.

## 10. MATLAB quirks and surprises

1. `XR_runDetTrack3d` drops unknown options: the multiple-experiments demo's
   `'zoffsetCorrection', true` (the PDF documentation's spelling too) does nothing.
2. The detection's Amira output fails (a warning without identifier, its third-party text
   not kept here) and leaves an empty folder.
3. The t threshold cache (section 3.4).
4. `data.date` comes from digits in the path.
5. MATLAB's `tcdf`, `tinv`, `chi2inv` and `fcdf` are less exact than SciPy's in places
   (section 7).
6. A parallel pool starts when the session's settings allow it (deskew's, the Amira writer's
   and the lifetime analysis' parfor loops); results do not depend on it.
7. The lifetime analysis fails on a condition with an empty track category
   (`stats:ecdf:NotEnoughData`, section 6).
8. `XR_estimate_overall_background_information` draws from MATLAB's random number generator
   (`'randseed', 1`; section 3.4).
