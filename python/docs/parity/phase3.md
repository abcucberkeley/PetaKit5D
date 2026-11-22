# Differences from MATLAB: the point-detection demos (Phase 3)

The point-detection demos (`demos/demo_point_detection_3d_*.m`) are examples in
`examples/phase3/`, one function per demo step, with the demos' parameters
([guide/examples.md](../guide/examples.md)). This page compares each step's outputs with
MATLAB's on the demo dataset, end to end and each step on MATLAB's own inputs, and gives each
step's speed against MATLAB's. Rows marked "proposed" wait for the user's approval; "open"
marks a difference not explained yet, reported to the step's owner.

## The steps

| Demo | Step (MATLAB) | Python | Example | Status |
|---|---|---|---|---|
| all three | Step 1: load data (`XR_loadConditionData3D`, interactive in the demos) | `detection.load_condition_data` with the demos' answers (channels `ch1`, `ch2`; markers `gfp`, `rfp`) | `_phase3.load_data` | ported |
| | Step 2: PSF sigmas (`GU_estimateSigma3D` of both bead images, `sigmaZ ./ zAniso`) | `detection.sigma.estimate_sigma3d` | `_phase3.psf_sigmas` | ported |
| demo_point_detection_3d_by_step | Step 3: deskew (`deskewData`) | `detection.deskew_data` | `deskew` | ported |
| | Step 3: XZ offsets (`XR_correctXZoffsetData3D`) | `detection.correct_xz_offset_data` | `xz_offsets` | ported |
| | Step 4: detection (`XR_runDetection3D`, 'lowSNR' with the background correction) | `detection.run_detection` | `detection` | ported |
| | Step 5.1: tracking (`XR_runTracking3D`, `loadTrackSettings('Radius', [3, 6], 'MaxGapLength', 2)`) | `detection.tracking.run_tracking`, `track_settings` | `tracking` | ported |
| | Step 5.2: track processing (`XR_runTrackProcessing3D`) | `detection.run_track_processing` | `track_processing` | ported |
| | Step 5.3: rotate tracks (`rotateTracks3D`) | `detection.rotate_tracks` | `rotate` | ported |
| demo_point_detection_3d_integrated | Step 3: `XR_runDetTrack3d` (deskew, XZ offsets, the frame-parallel detection, tracking, track processing, rotation) | `detection.run_detection_tracking` | `detection_tracking` | ported |
| demo_point_detection_3d_multiple_experiments | Step 3: `XR_runDetTrack3d` of the three experiments, without the offsets (row 6) | `detection.run_detection_tracking` | `detection_tracking` | ported |
| | Step 4: lifetime analysis (`runLifetimeAnalysis3D`) and intensity cohorts (`plotIntensityCohorts`, its numbers) | `detection.lifetime.run_lifetime_analysis`, `plot_intensity_cohorts` | `lifetime_analysis`, `intensity_cohorts` | ported |

## The functions in the mapping table

The demos call these functions; [MATLAB_MAPPING.md](../MATLAB_MAPPING.md) has the registered
ones by their MATLAB name. Those without an entry are ported and called through the Python API
only.

| MATLAB | Python | Mapping table |
|---|---|---|
| `XR_loadConditionData3D` | `detection.condition.load_condition_data` | registered (interface only: third-party code) |
| `GU_estimateSigma3D` | `detection.sigma.estimate_sigma3d` | no entry |
| `deskewData` | `detection.wrapper.deskew_data` | no entry (third-party code) |
| `XR_correctXZoffsetData3D` | `detection.offsets.correct_xz_offset_data` | no entry |
| `XR_runDetection3D`, `XR_runDetection3D_frame_parallel` | `detection.run.run_detection`, `run_detection_frame_parallel` | registered (interface only) |
| `loadTrackSettings` | `detection.tracking.settings.track_settings` | no entry |
| `XR_runTracking3D` | `detection.tracking.run.run_tracking` | registered (interface only) |
| `XR_runTrackProcessing3D` | `detection.processing.run_track_processing` | registered (interface only) |
| `rotateTracks3D` | `detection.processing.rotate_tracks` | no entry (third-party code) |
| `XR_runDetTrack3d` | `detection.wrapper.run_detection_tracking` | registered |
| `runLifetimeAnalysis3D` | `detection.lifetime.run_lifetime_analysis` | registered (interface only) |
| `plotIntensityCohorts` | `detection.lifetime.plot_intensity_cohorts` | no entry (third-party code) |

## How the outputs are compared

- **The runs.** MATLAB R2024b (the oracle tree of 2026-11-19) ran each demo's steps as the demo
  runs them, with the interactive loading replaced by its call with the demo's answers, on a
  fresh copy of the demo dataset on a node's local disk (the step-by-step and integrated demos:
  the first experiment, 60 frames of two channels; the multiple-experiments demo: all three).
  The examples ran the same way on a copy of their own, on the same node: two rounds of each
  side, alternating, for the step-by-step and integrated demos, whose first round's trees are
  compared; for the multiple-experiments demo a round of each side in each of two jobs on one
  node, whose second job's trees are compared (the first job's were not kept; its comparisons
  gave the same counts).
- **End to end** (`examples/phase3/compare_phase3_trees.py`, on `petakit5d.testing.compare`):
  every output of every step, by name; tiff files by voxels and tags, `.mat` files by value,
  the Amira files by their numbers, the in-memory results (the condition, the PSF sigmas) as
  MAT-files; each run's root folder masked. Each step starts from the step before it on its
  own side, so the PSF sigmas' differences (D6) reach every fit after them.
- **Each step alone** (`--matlab-inputs`): the detection on MATLAB's deskewed frames with
  MATLAB's sigmas, the tracking on MATLAB's `Detection3D.mat`, the track processing on MATLAB's
  detections, tracks and masks, the rotation on MATLAB's `ProcessedTracks.mat`, the lifetime
  analysis and the cohorts on MATLAB's `ProcessedTracks.mat`: what each step itself does
  differently. (The load, the sigmas, the deskew and the offsets take the same inputs on both
  sides end to end.)
- **The lifetime files** (`lifetimeData.mat`, `Parameters.mat`, `lftRes`, the cohorts) are
  compared field by field: MATLAB's bit for bit but the threshold `MaxIntensityThreshold` (L4)
  and the cohorts' `A`, `Amin` and `Aplus` (L7), as in `tests/detection/test_lifetime.py`;
  end to end every field also within E1.
- **Cut-down data in the tests**: `tests/e2e/test_demos_phase3.py` runs the three examples
  against MATLAB's demo steps on a small synthetic dataset of the demo's layout (oracle cases:
  one experiment, two for the multiple-experiments demo, of 16 frames of 24 x 56 x 48 raw
  voxels per channel, with 12 and 5 tracks), and the examples' own logic. The
  multiple-experiments case stops after step 3: MATLAB's lifetime analysis needs tracks in
  every category and stops with `stats:ecdf:NotEnoughData` on small conditions
  ([detection_scope.md](../detection_scope.md), section 6), so step 4 is compared on the demo
  dataset only.

## Results

### demo_point_detection_3d_by_step (the first experiment, 60 frames)

| Step | MATLAB's outputs | End to end | Each step alone | Rows |
|---|---|---|---|---|
| Step 1: load data | the condition (`data`) | identical | | |
| Step 2: PSF sigmas | `sigma_mat` | within 1.2e-9 relative | | D6 (approved) |
| Step 3: deskew | 120 DS frames, `settings.mat`, `dynRange.mat` | 122 of 122 identical | | |
| Step 3: XZ offsets | `xz_offset.mat` of both channels | 2 of 2 identical | | |
| Step 4: detection | `Detection3D.mat`, 60 masks | 58 masks identical; `Detection3D.mat` and the masks of frames 33 and 55 differ (rows 1 to 3) | the same files differ (rows 2, 3) | D3 (approved); row 2 (open); E1 (proposed) |
| Step 5.1: tracking | `trackedFeatures.mat` | 10,915 tracks against 10,917 (row 1) | identical | TK1-TK4 (approved); E1 (proposed) |
| Step 5.2: track processing | `ProcessedTracks.mat`, 60 Amira files | 1 Amira file identical, the others differ (row 1) | the same: W13's fits, W14's four categories (row 4) | W11-W15 (approved); E1 (proposed) |
| Step 5.3: rotate tracks | `RotatedTracks.mat` | differs (row 1) | identical | W16 (approved); E1 (proposed) |

MATLAB's track processing also leaves an empty folder `ch1/Analysis/` (row 5); empty folders
are not compared.

### demo_point_detection_3d_integrated

`XR_runDetTrack3d` gives the step-by-step demo's results in every file (on both sides: MATLAB's
`Detection3D.mat`, `trackedFeatures.mat`, `ProcessedTracks.mat` and `RotatedTracks.mat` of the
two demos are the same, and so are the port's), and adds the frame-parallel detection's 60
per-frame files (`Detections/Detection_F_<k>.mat`), which differ as `Detection3D.mat` does. So
the comparison is the step-by-step demo's: the load, the sigmas, the deskew and the offsets
identical or within D6; the detection, the tracking, the track processing and the rotation as
in rows 1 to 4, end to end and alone, with the same counts (frames 33 and 55; 10,915 against
10,917 tracks end to end; the same 55 values of W13 and the same four categories of W14 alone).

### demo_point_detection_3d_multiple_experiments (the three experiments, 60 frames each)

| Step | MATLAB's outputs | End to end | Each step alone | Rows |
|---|---|---|---|---|
| Step 1: load data | the condition (3 movies) | identical | | |
| Step 2: PSF sigmas | `sigma_mat` | within 1.2e-9 relative | | D6 (approved) |
| Step 3: deskew | 360 DS frames, 3 `settings.mat`, 3 `dynRange.mat` | 366 of 366 identical | | |
| Step 3: detection | 3 `Detection3D.mat`, 180 masks, 180 per-frame files | 148 identical (masks); the first experiment's counts as in the step-by-step demo (frames 33 and 55: row 2), the second's and third's other in 20 and 10 frames (row 10), and the fits after D6 (row 1) | the same frames and counts | D3 (approved); rows 2 and 10 (open); E1 (proposed) |
| Step 3: tracking | 3 `trackedFeatures.mat` | 10,915, 17,163 and 14,293 tracks against MATLAB's 10,917, 17,079 and 14,271 | identical | TK1-TK4 (approved); E1 (proposed) |
| Step 3: track processing | 3 `ProcessedTracks.mat`, 180 Amira files | the same 3,636 and 5,193 tracks; 6,558 against 6,553 in the second experiment; their values as row 1 | the same tracks; the first experiment's as the step-by-step demo's (W13's 55 values in 11 tracks, W14's four categories), the second's and third's beyond W13 in 90 and 15 tracks and with 11 and 19 other categories (row 11) | W11-W15 (approved); row 11 (open); E1 (proposed) |
| Step 3: rotation | 3 `RotatedTracks.mat` | as the processed tracks | identical | W16 (approved); E1 (proposed) |
| Step 4: lifetime analysis | 3 `lifetimeData.mat`, `Parameters.mat`, `lftRes` | follow the tracks (row 12) | the 3 `lifetimeData.mat` identical; `MaxIntensityThreshold` 9.9e-10 relative apart, every other value identical | L4 (approved) |
| Step 4: intensity cohorts | `cohorts` (`A`, `Amin`, `Aplus` per channel and cohort) | follow the tracks (row 12) | 428 of 656 values differ, by at most 2.3e-8 relative (5.1e-8 absolute) | L7 (approved) |

The demo corrects no offsets (row 6), so its first experiment's frames and detections are not
the step-by-step demo's; MATLAB's and the port's are compared all the same way. Each step alone
on MATLAB's results, the lifetime analysis is MATLAB's but for the threshold and the cohorts,
by the numbers of its port's demo case (`lifetime.demo` of `tests/detection/test_lifetime.py`,
on the tracks of an earlier MATLAB run of this demo). End to end, the tracks of the second and
third experiments differ in number (row 10), and what is made from them with them.

## Speed

Per step, the median of the rounds (and their range), on one node each (24 cores, Xeon
Gold 6146, 512 GB), both sides with 24 threads under `numactl --interleave=all`, MATLAB R2024b
with the parallel pool its settings start (24 workers, from deskewData's `parfor`; its start
is part of MATLAB's deskew). Times include each step's file reading and writing.

**demo_point_detection_3d_by_step** (2 rounds of each side, alternating):

| Step | MATLAB (s) | Python (s) | Python / MATLAB |
|---|---|---|---|
| Step 1: load data | 0.42 (0.34-0.49) | 0.19 (0.11-0.26) | 0.45 |
| Step 2: PSF sigmas | 1.79 (1.72-1.87) | 1.11 (1.06-1.16) | 0.62 |
| Step 3: deskew | 32.3 (31.3-33.4) | 4.08 (4.03-4.12) | 0.13 |
| Step 3: XZ offsets | 77.8 (65.2-90.4) | 64.0 (62.7-65.3) | 0.82 |
| Step 4: detection | 182.0 (177.1-186.9) | 103.3 (103.1-103.4) | 0.57 |
| Step 5.1: tracking | 33.8 (33.5-34.1) | 7.43 (7.36-7.49) | 0.22 |
| Step 5.2: track processing | 193.0 (192.5-193.5) | 77.7 (77.5-77.9) | 0.40 |
| Step 5.3: rotate tracks | 18.6 (18.2-18.9) | 6.75 (6.72-6.77) | 0.36 |
| all steps | 539.7 | 264.5 | 0.49 |

**demo_point_detection_3d_integrated** (2 rounds of each side, alternating):

| Step | MATLAB (s) | Python (s) | Python / MATLAB |
|---|---|---|---|
| Step 1: load data | 0.36 (0.36-0.36) | 0.17 (0.10-0.23) | 0.47 |
| Step 2: PSF sigmas | 1.86 (1.65-2.06) | 1.10 (1.06-1.13) | 0.59 |
| Step 3: `XR_runDetTrack3d` | 586.8 (576.1-597.6) | 264.8 (261.8-267.7) | 0.45 |

**demo_point_detection_3d_multiple_experiments** (2 rounds of each side, in two jobs on one
node; Python's step 4 also run 3 times more on its second round's processed tracks, on the same
node):

| Step | MATLAB (s) | Python (s) | Python / MATLAB |
|---|---|---|---|
| Step 1: load data | 0.43 (0.41-0.45) | 0.22 (0.22-0.22) | 0.52 |
| Step 2: PSF sigmas | 1.75 (1.69-1.81) | 1.06 (1.06-1.06) | 0.61 |
| Step 3: `XR_runDetTrack3d` of the three experiments | 1853.2 (1844.3-1862.2) | 772.7 (761.1-784.3) | 0.42 |
| Step 4: lifetime analysis | 20.5 (20.1-20.9) | 9.30 (9.25-9.65) | 0.45 |
| Step 4: intensity cohorts | 0.67 (0.65-0.69) | 0.17 (0.17-0.17) | 0.26 |

The stages' own benchmarks (`benchmarks/detection/`: five or more interleaved repeats per
stage on one node) gave the same picture: a detection frame 0.62 of MATLAB's time at 4
threads and 0.41 at 24, tracking 0.24, track processing 0.42, rotation 0.35, the offsets
0.82. The XZ offsets are the one step near MATLAB's time: their 165 correlations per frame
pair are summed in MATLAB's order.

## Differences

| # | Output | Difference | Status | Why |
|---|---|---|---|---|
| 1 | end to end, what the PSF sigmas feed: `Detection3D.mat` (and the per-frame files), `trackedFeatures.mat`, `ProcessedTracks.mat`, the Amira files, `RotatedTracks.mat` | each side's sigmas (row D6: 1.2e-9 relative apart) give every fit other last bits: of the 118,604 values per field of the 58 frames with the same detections, all but 46 to 51 positions, 1,017 amplitudes (the other channel's mostly; largest 27% of max(\|A\|, 1)) and 195 p-values within `rel` 1e-8 or `abs` 1e-8; 10,915 tracks against 10,917 (by their positions 14 of MATLAB's and 12 of the port's differ: the tracks through row 2's detections and through fits that moved); the processed tracks the same 3,636, with W14's categories, all but 41 to 729 values per field (in 619 tracks) within E1's bound | proposed (E1, 2026-11-21): `rel` 1e-8 with `abs` 1e-8 on the cut-down data (`tests/e2e/test_demos_phase3.py`: largest 1.2e-9); on the demo dataset these numbers | D6's sigmas set the windows and starts of every fit after them, and D3's fits (the ones that take another path) and row 2 carry on into the tracks. Each step alone is within its own rows (rows 2 to 4) |
| 2 | `Detection3D.mat` and the masks of frames 33 and 55, each step alone and end to end | one detection fewer in the port in each of the two frames: MATLAB's detector (called alone on the frame, as its runner calls it) has a fit at z 33.4 of the frame's 35 planes and one at z 1.0; the port's nearest detections are 6.0 and 13.2 voxels away. The other 58 frames have MATLAB's detections | open: for the detection's owner | not found yet: both fits end at a border plane of the deskewed frame |
| 3 | `Detection3D.mat`, each step alone: the matched detections | positions of the master channel bit for bit or within 1e-12 voxels in 59,284 of 59,302 detections, 16 beyond 1e-9 (largest 9.7e-5 voxels, A 1.6e-4 relative); the other channel's 46 beyond 1e-9 (largest 0.57 voxels, A 0.66 relative); one detection's Anderson-Darling test (`hval_AD`, `isPSF`) in frames 15 and 48; 6 standard deviations NaN on one side (frames 5, 12, 15) | approved (D3: "demo frames: the same detections, the largest differences listed"; this is that list for 60 frames) | the fits' last bits (D2, D3): the restarted fits and the other channel's `'xyzAc'` fits from the master's positions amplify them |
| 4 | `ProcessedTracks.mat` and the Amira files, each step alone | the same 3,636 tracks; 55 values of 11 tracks beyond W11's `rel` 1e-11 with `abs` 1e-10; `catIdx` of the 4 tracks of W14 (1849, 2181, 2322, 2968, 0-based); the Amira files: 59 of 60 differ, in 48,037 of their 1,484,766 numbers | approved (W11 to W15) | W13's fits that take another path, W14's undecoded rules, W15's printed values |
| 5 | an empty folder `<source>Analysis/` beside the results folder | MATLAB's `XR_runTrackProcessing3D` makes it with any `ResultsPath` (probe: it is not there after the tracking, and there after the track processing); the port does not | documented: empty folders are not compared ([demos.md](demos.md)) | a folder nothing is written to (the default results folder, made although `ResultsPath` is given) |
| 6 | the multiple-experiments demo's channel offsets | the demo passes `'zoffsetCorrection', true`, which is not an option of `XR_runDetTrack3d`; its parser drops it, so that demo corrects no offsets (unlike the integrated demo's `'XZoffsetCorrection', true`); the example does the same | documented: MATLAB's result kept | the demo's option name (`deskewData` has a `zoffsetCorrection`, which `XR_runDetTrack3d` does not pass on) |
| 7 | the step-by-step demo's step 5 after a one-frame movie | the demo drops the one-frame movies from `data` but not their folders from `apath`, so later movies would get the folders before theirs (W3's slip, in the demo script); the example drops both | documented: the example as W3 fixes the wrapper | the demo dataset has no one-frame movie: no difference on it |
| 8 | `Detection3D.mat` of a cut-down synthetic movie (found while making the tests' dataset) | MATLAB's `XR_runDetection3D` and `XR_runDetection3D_frame_parallel` leave out one weak detection (amplitude 5.3 over a background of 100, residual sigma 2.86, the only one of the movie's weak noise detections below the camera threshold of 3; z 4.6 of 24 planes) that MATLAB's `XR_pointSourceDetection3D_3`, called alone on the frame with the runner's options, returns; with any background seed, residual-sigma threshold or mask, with the other channel or without it; the port's runners keep it. On that frame MATLAB's detector also has 10 candidates before its residual-sigma filter, the port's 9 | open: for the detection's owner | not found: the runners' measured behaviour (14 demo frames, C02) never showed it. The tests' dataset keeps away from it (its weak noise detections all lie above the camera threshold) |
| 9 | `ProcessedTracks.mat` and `RotatedTracks.mat` without tracks (every track dropped by the border filter; a cut-down movie of 12 planes) | none since 2026-11-21: MATLAB writes `tracks` as a 0 x 0 struct array of the processed tracks' fields but `significantMaster`, `significantVsBackground` and `significantSlave`, with `processingInfo`, and rotates it into the same 0 x 0 array (no `processingInfo`); its Amira files are the header alone (no vertex, edge or point); the port writes the same (it wrote a 1 x 0 array with every field) | fixed (C04c): oracle case `run.empty` of `tests/detection/test_processing.py`, every file bit for bit | an empty result's layout, measured by the case |
| 10 | `Detection3D.mat` and the per-frame files of the multiple-experiments demo's second and third experiments, each step alone and end to end | the port's detector finds other numbers of detections in 20 of the second experiment's 60 frames (134,349 detections against MATLAB's 134,259; frame 16: 2,203 against 2,180) and 10 of the third's (93,974 against 93,953); the same frames and counts on MATLAB's deskewed frames with MATLAB's sigmas as end to end. Matched by position (0.5 voxels), each step alone: weak detections (A 4.9 to 8.0 over a background of about 100) that only the port has, 99 in 10 frames of the second experiment (23 in frame 16) and 28 in 5 of the third, whose residual sigmas lie between 2.76 and 4.26 (the third's all between 2.80 and 2.99, just under the camera threshold of 3), and 3 that only MATLAB has (the third's frame 10, residual sigmas 2.85 to 2.94); bright detections (A 10 to 96) that only MATLAB has, 12 in the second and 4 in the third experiment, 7 of them on the frame's first or last plane (as row 2), and 3 (A 55 to 1,380) that only the port has | open: for the detection's owner | not found yet. D3's measurements were the first experiment's frames, where the counts differ in frames 33 and 55 only (row 2); the second and third experiments are denser (about 2,200 and 1,570 detections per frame against 1,020). The weak ones are of row 8's kind (a weak detection near the camera threshold that MATLAB's runner leaves out), decided per frame the other way; the bright ones on a border plane of row 2's |
| 11 | `ProcessedTracks.mat` of the multiple-experiments demo, each step alone | the same tracks (3,636, 6,553, 5,193) from MATLAB's detections and tracks. Since 2026-11-21 (C04c) every `catIdx` is MATLAB's (before: 4, 11 and 19 tracks differed). Beyond W11, before and after alike (C07's comparison): the first experiment's 55 values in 11 tracks, the second's 456 in 90 and the third's 75 in 15, by up to 1.33 voxels in the tracks' positions (end buffers 1.77), 0.55 relative in their amplitudes (end buffers 0.63) and 0.064 relative in `MotionAnalysis.MSD` | the categories fixed (W14); the values W13's, its multiple-demo addition proposed 2026-11-21 | the categories: the buffer test compares a buffer frame's intensity `A + c` with the track's largest, not its amplitude (W14). The values: MATLAB's track processing, recorded on the 239 tracks whose fits differ, reads every window and start of their 4,320 fits as the port does, bit for bit, and makes the same keep-or-refit decision in 171,923 of the demo's 171,926 fits; the free fits (`fit_cd`, D3) end elsewhere than MATLAB's on 465 of those windows (weak spots, mostly), and three of them on the other side of the keep rule's 2 voxels, which gives the voxel-sized differences (W13; reported to the fits' owner) |
| 12 | the lifetime analysis' and the cohorts' files of the multiple-experiments demo, end to end | they follow the tracks: the per-movie arrays differ in size (the second experiment's 6,558 processed tracks against 6,553, row 10; the first's valid tracks 557 against 555, from W14's four categories; the third's 827 against 818), `MaxIntensityThreshold` by 0.2% (31.39 against 31.33), `nCCP` by 1.6%, the lifetime histograms in up to 47 of their 48 bins, the cohorts' amplitudes by up to 0.78 relative (per cohort a relative L2 of up to 0.11) | follows from rows 10 and 11 and W14 | each step alone the lifetime analysis is MATLAB's (L4, L7) |

MATLAB's behaviour found on the way, off the demos' path: `XR_pointSourceDetection3D_3` with
`'RemoveRedundant', false` fails (`MATLAB:UndefinedFunction`).

## The Phase 3 rows of the exception index

[index.md](index.md), section 3, has the rows of Phase 3's functions; their measurements are in
[detection.md](detection.md) and [sim.md](sim.md):

- S1-S10: the statistics functions of point detection;
- D1-D12: point detection, the `fitGaussian3D` mex, the PSF sigma estimate, the runners;
- TK1-TK4: tracking;
- W1, W2-W4, W11-W20: the condition loading, deskew and offsets, the wrapper's fixed slips, track
  processing and the tracks' rotation;
- L1-L7: the lifetime analysis;
- R1-R6: structured illumination (no demo).

This page adds the demo-level rows of the Differences table above: E1 (row 1, proposed) and the
open rows 2 and 8 to 11, reported to the stages' owners.

## Not compared

- the figures of `plotIntensityCohorts` (the port gives its numbers);
- the empty `Amira/Detection/` folders (MATLAB's Amira writer of the detections fails and
  leaves them; the port makes them too) and the empty `ch1/Analysis/` folder (row 5).
