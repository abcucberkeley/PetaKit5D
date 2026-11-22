# Differences from MATLAB: point detection

The port matches MATLAB bit for bit unless an oracle case declares a tolerance
(`petakit5d.testing.compare.tol`, see [testing.md](../testing.md)). This page lists the
tolerances of `petakit5d.detection`, with what was measured and why the bits differ.
Each case's `reason` names its row as `detection_scope.md statistics Sn` (the scope document
has the same rows). All rows were approved by the user on 2026-11-20
([index.md](index.md)), but D13, proposed (2026-11-21; measured per case, with the evidence of
MATLAB's intent, on 2026-11-22).

## Statistics functions (`detection.stats`, C01)

MATLAB's statistics functions are reproduced from MATLAB's documentation and black-box probes
(no MathWorks source read). MATLAB's special functions (`erfcinv`, `betainc`, `gammainc` and
their inverses) are closed; the port takes SciPy's. Where the differences are large, MATLAB is
the one off the exact value (checked with 50-digit arithmetic).

Exact, bit for bit on every probed value: `mad` (flag 0 and 1, vectors, columns, NaN), `ecdf`
(the Kaplan-Meier form `1 - cumprod(1 - d./r)`, ties, NaN), `kstest2` (statistic and p-value,
all three tails, `'Alpha'`), the edge values of the five distribution functions, `tcdf(x, 1)`
and `tinv(p, 1)`.

| Row | Function and range | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| S1 | `norminv`, `p` in [realmin, 1]; `tinv(p, Inf)` | 4 ulp; `mu + sigma*z` near 0: 2e-15 of the largest value | `ulp` 4 or `rel_max` 2e-15 | MATLAB's `erfcinv` against SciPy's `ndtri` |
| S2 | `norminv` of subnormal `p` | 1.3e-4 relative | `rel` 2e-4 | MATLAB loses digits below `realmin` |
| S3 | `tcdf`, `nu` to 1e4 and above 1e7 | 4.5e-12 relative; values under 1e-300 | `rel` 1e-11, `abs` 1e-300 | MATLAB's `betainc`; MATLAB keeps subnormal tails that SciPy rounds to 0 |
| S4 | `tcdf`, `nu` 1e6 to 1.5e7 | 1.3e-7 relative | `rel` 2e-7 | MATLAB's form for large `nu` |
| S5 | `tinv`, `p` in [1e-6, 1 - 1e-6] | 1.2e-10 relative (`nu` 10, `p` 1e-6) | `rel` 2e-10, `abs` 1e-14 | MATLAB's inverse; the port is within 2 ulp of the exact value |
| S6 | `chi2inv`, `p` in [1e-6, 1] | 7.3e-12 relative | `rel` 1e-11 | MATLAB's `gammaincinv` |
| S7 | `fcdf` | 8.4e-12 relative | `rel` 1e-11, `abs` 1e-300 | MATLAB's `betainc` |
| S8 | `fcdf`, `nu1` 1e6 | 7.9e-9 relative | `rel` 1e-8 | MATLAB's form for a large `nu1` |
| S9 | `ksdensity` (normal kernel, default bandwidth), samples of 100 to 200 | 16 ulp of the density | `ulp` 16 | the kernel sum's order (MATLAB's grid and bandwidth are exact) |
| S10 | `ksdensity` from 500 values on | 1.5% (2% bound) | `rel` 2e-2 | MATLAB evaluates the density another way for long inputs; no traced run calls `ksdensity` |

The S rows decide nothing on the demo: on its frame 1 (MATLAB's own `tcdf` inputs through
`detection.stats`), none of 3,274 `p < alpha` decisions flips (largest relative difference
2.3e-11, the closest p-value 21.5 % away from alpha); `tinv(0.8, df)` differs by 6.7e-16 and the
mask voxel nearest to the threshold is 1e-4 away; `chi2inv(0.9, 343)` is equal.

## Detection and the Gaussian fits (`detection.detect`, `detection.fit`, `detection.run`, C02)

Reproduced from black-box probes of MATLAB R2024b (logging wrappers with full values, profiler
line numbers, designed inputs; no GPL or MathWorks source read). Bit for bit: the overall
background estimate (both channels of the demo's frame 1, 3 oracle cases), the background
correction, the kernel, the mask (32,967 of 32,967 voxels), the fits' starts (1,264 of 1,264),
`XR_locmax3d`, the NaN rule of the fits, the redundancy filter and the amplitude test's
selection (1,040 of 1,040 in order), `adtest_mex` (3,274 traced calls, 18 oracle cases),
`histcounts`, `findpeaks`, `medfilt1` (2,665, 2,932 and 85 recorded calls), `rangesearch`, and
the frame assembly of `XR_runDetection3D` given MATLAB's detector and fit outputs (14 real
frames, every field of `frameInfo`, the dmask files).

| Row | What | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| D1 | `fspecial3('log')` | 1.5e-15 relative | `abs` 1e-16, `rel` 2e-15 | the kernel's sums in another order (the formula is MATLAB's) |
| D2 | `imfilter` of 3-D kernels (the detector's window sums of `v`, `v^2` and the LoG) | 5e-15 relative (the sums of ones exact); the fits' starting A and c 1.1e-13 | none of its own (carried by D3) | MATLAB's order of the 343 products is not measured (no sequential order matches) |
| D3 | the detector's fits and results (positions, A, c, their standard deviations, RSS, sigma_r, p-values) | demo frame 1: positions of 1,043 of 1,043 plain fits and 188 of 221 restarted ones within 1e-12; over all fits positions 2.5e-10, A 2.7e-7, standard deviations 1.3e-7 relative; the oracle cases within 1e-9. Re-measured with C08's stop test (2026-11-21): of 23,969 traced fits of 8 demo frames and C02's probe 14, 23,737 within 1e-12 of MATLAB's parameters (23,510 with C02d's box, 23,171 before it); the step-by-step demo's 60 frames each step alone (phase3.md row 3): the master channel's positions within 1e-12 in 61,345 of 61,353 detections (61,334 before), 6 and 16 fits of the two channels beyond 1e-9 (17 and 48). Decoded fit by fit (C09, "The fits still on another path" below): the 232 traced fits beyond 1e-12 all part from MATLAB's path in their last bits first (225 at MATLAB's update of A, 6 at damped steps, 1 at the box itself), and 229 of them part decisively at the box, where MATLAB's position lies an ulp inside a face and the port's on it; the step-by-step demo's 22: of the master channel's 6, 5 so and 1 from its start (D2); of the other channel's 16, 13 from their starts (MATLAB given the port's window and start ends at the port's end: their `'Ac'` starts differ, D4) and 3 so | `rel` 1e-9, `abs` 1e-9 on the oracle cases; demo frames: the same detections, the largest differences listed | D2's last bits, the restarts' direction solves (the damped steps' eigenvalues) and the A update's sums, whose order MATLAB keeps in ways not reproduced (measured, C09: none of 20 orders and forms gives MATLAB's A in more than 31% of the traced iterations); on paths whose last bits differ so, MATLAB's cut steps land an ulp inside the box at times where the port's land on it, and MATLAB takes one more step there (C02d, C08, C09) |

## The fitGaussian3D mex, the PSF sigma estimate and the 'speedup' detector (C02b)

`fitGaussian3D` is a mex without source. Measured as a black box (its inputs and outputs, and
its own calls of `exp`, `erf` and `log` observed through an interposed math library): the
model is evaluated per axis, `((A ex) ey) ez + c - d` with `ex(k) = exp(-(k - x)^2 / (2 s^2))`
on 0-based window coordinates, bit for bit given MATLAB's parameters; the parameters of the
mode's letters are fitted by Levenberg-Marquardt in Moré's trust-region form (a Jacobian at
every accepted point, at most ten rejected trials per iteration, 500 iterations); the outputs
are `inv(J'J)`, `prmStd = sqrt(diag(C) RSS / (n - p - 1))`, the residuals' RSS, mean, standard
deviation and an Anderson-Darling test of the sorted residuals over their standard deviation
(critical value 2.308: the same decision on 33,686 of 33,686 traced calls). The mex calls the C
library's `exp` (glibc's, correctly rounded for 99.85% of the arguments it saw); the port's is
correctly rounded. The solver's internal arithmetic (its QR factorization and norms) is not
reproduced, which the rows below carry.

Bit for bit (no row): `estGaussianAmplitude3D` (9 traced calls, 5 oracle cases) and its
`conv3fast` filters (18 traced, 2 oracle cases; mirrored borders, centre-out symmetric sums,
passes x, y, z; NaN with the input's payload), `fspecial3('gaussian')` (12 traced calls), the
camera threshold (`XR_get_camera_standard_deviation_threshold`, 12 oracle cases with MATLAB's
errors), the 'speedup' detector's background correction, mask (the demo frames' 20,094, 19,906
and 20,230 voxels; the PSFs' 306, 448, 452, 114, 198 and 196) and starts.

| Row | What | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| D4 | the mex's linear modes (`'Ac'`, `'A'`, `'c'`; track processing's `'Ac'` fits) | 33,686 traced calls: A within 8.5 ulp and c within 2 ulp of max(\|A\|, \|c\|) (A bit for bit in 8 %, c in 92 %); prmStd, C, RSS and std within 7e-15 relative; hAD the same in all | `rel` 1e-13, `rel_max` 1e-13, `abs` 1e-12 (the residuals' mean is 0 to rounding) | the solver's QR and norm arithmetic, from the first step on |
| D5 | the mex's nonlinear modes (`'xyzAsrc'` of the sigma estimate, `'xyzAc'`) | the sigma estimate's 31 traced fits: parameters within 9e-8 relative; the residual window 1.4e-7 | `rel` 1e-6, `abs` 1e-6 | D4's last bits steer the trust-region path; the fit stops where the residuals' decrease is at rounding level, about 1e-8 from the minimum, so the two stop at different points |
| D6 | `GU_estimateSigma3D`: `sigmaXY`, `sigmaZ`, the bead's position | the demo's PSFs: 1.5e-10 and 6.7e-11 (560 nm), 3.0e-10 and 1.2e-9 (642 nm) relative; the oracle case within 1e-7 | `rel` 1e-7 | D5 (three rounds of the 'speedup' detector with the mex's fits) |
| D7 | the 'speedup' detector `XR_pointSourceDetection3D` | `'CD_Newton'` (XR_runDetection3D's): 4 demo frames, the same 661, 649 and 668 detections, values within 1.5e-14 for 99 %, restarted fits 3.3e-9; with the mex (the sigma estimate's form) as D5 | `rel` 1e-9, `abs` 1e-9 (D3's) with `'CD_Newton'`; D5's with the mex | D2's window sums, D3's fits, D5's; the t threshold's `tinv` (S5: 5.7e-13) changes no mask voxel |
| D8 | `XR_runDetection3D`'s `frameInfo` with its defaults ('speedup', no background correction) | 2 demo frames, both channels: the same 661 and 649 detections; the master's values within 3.3e-9 relative (D7), the other channel's A within 1.6e-4, c 9.3e-7, positions 6e-15; every test (`hval_AD`, `isPSF`) the same | D7's, and `rel` 1e-3 for the other channel's A and c (no case yet) | D3 in the other channel's `'xyzAc'` fits from the master's positions, which move up to 4 voxels; one in 600 stops where the residuals barely change |
| D9 | the camera search on Windows (`get_camera_type_from_setting_file`: `data.source` split at `filesep`, each folder up listed with `dir`) | deduced from the lab code (no MATLAB on Windows here): MATLAB splits at `\` only, so it walks a source written with `/` (the oracle harness writes its paths so) no higher than its last `\` and misses a settings file above it (`MATLAB:UndefinedFunction`); `strsplit`'s merged separators would make a UNC path's `\\server` `\server` | none: the port splits at both separators, which MATLAB's `dir` and `fileparts` take on Windows, and keeps a UNC prefix; the oracle cases on Windows find the settings files MATLAB finds on Linux | MATLAB's bug, fixed in the port |

The camera search's walk ends at an absolute path's first part: on Linux `''`, which `dir`
lists as nothing (`MATLAB:subsdeldimmismatch` without a settings file, measured); on Windows the
drive `C:`, which lists the drive's root, without `.` and `..` there (so
`MATLAB:UndefinedFunction`), and for a UNC path `\\server`, no folder
(`MATLAB:subsdeldimmismatch`): deduced from the lab code, Windows' file systems and MATLAB's
single current folder, and tested with Windows' rules on Linux.

The runner (`detection.run.run_detection`, measured on probes rd1 and rd2 of two-movie and
one-movie conditions): the default `ResultsPath` `[source 'Analysis' filesep]`, the default
method 'speedup' (MATLAB's default call and its 'speedup' call write the same file; 'original'
another), `Overwrite` false by default and a movie with a `Detection3D.mat` skipped, the file at
`ResultsPath + 'Detection3D.mat'` (a folder without the separator gets the file beside it) and
the masks in `ResultsPath/Masks`, a char `ResultsPath` refused, a sigma row per channel with the
rows past the channels kept in `s` only. An `'xyzAc'` fit of `XR_fitGaussians3D` is NaN when it
ends `w` or more from its window's centre, the rounded start: 0 of 2,315 traced fits from the
other channel's positions differ (the port measured from the unrounded start before, 282
differed). Open: a given `WindowSize` (MATLAB finds 4,830 detections on a demo frame with
`[2 2]`, the port's window 3,128; refused); 'lowSNR' without the background correction: D10.

## The bad-point filter's second valley search, the frame-parallel runner (C02c)

The bad-point filter (`detection.badpoints`) now runs MATLAB's second valley search (lines 221-267
of `filter_bad_points_in_detection`; line numbers and keep masks only, no source read), measured
with about 620 designed histograms (MATLAB's histogram and peak list mocked, every point given k
neighbours so that the keep mask shows the removal thresholds), two synthetic inputs without
mocks and C02's recorded calls. On the A paths: when the three bins from `max(bp - 1, 1)` hold
more than 0.12 of all counts, or the first two of them more than 0.09, the first bin j after the
background peak and before the next listed peak below `cbp / d` (d = 4, 3.5, 3, 2.5, 2 in turn)
gives `T2`, its upper edge. Found at one of the first four levels, the path is B with `T2` when
`T2 < 2 cam` (cam, the camera threshold, not 5; above 5 no bound was found up to 1.56 cam).
Otherwise, for cam below 5 and the first bin's centre below 1.5 cam (the background peak was bin
1 in the designs; of two recorded calls with it at bin 9 or 10, the one with cam 3 passed this
test, the one with cam 6 did not): B with `T2` when bin j's count is below `cbp / 5`; else, with
j found at one of the first three levels, the path A262 (lines 256-266): A's removal (2 or 3
neighbours) below `min(noise, 1.5 cam)`, the noise rule below at most that, B's median bound
(lines 351-357).
Otherwise, cam 5 included, the first valley's decision stands; the first valley's B skips cam 5
too (two designed histograms). The demo's 120 frames with the background correction stay bit for
bit (channel 1 frame 28 takes the search: 16 points more removed), 249 of 285 replayed variants
(one more, none lost), 942 of 1,011 designed decisions, 567 of the 580 designed keep masks with
mocked neighbour counts, and the seven oracle cases (camera thresholds 3, 2.4, 2.1, 2 and 1.85;
A262 with noise 2.6 and 3.2, on real neighbour counts). Open: 55 mismatches of C02's
provisional designs (e8-e11), as before; the smoothing choice of line 104 and B's median bound
on one frame (D10).

| Row | What | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| D10 | 'lowSNR' without the background correction (`XR_runDetTrack3d`'s default detection): `run_detection_frame_parallel` against `XR_runDetection3D_frame_parallel` on the demo's 60 frames, both channels | 53 frames with MATLAB's number of detections, 41 of them with every position within 1e-6 (the others: one or two fits apart, D3); of MATLAB's 61,260 detections 61,218 matched, the port has 61,290; frames 1, 33 and 55: 1, 9 and 1 detections differ, from the fits (with MATLAB's own points the port's filter keeps MATLAB's set on 56 of the 60 frames, probe lt4); frames 17, 21, 50 and 51: 19, 3, 25 and 14 detections, from the filter's open rules (the smoothing choice of line 104 on 17, 50 and 51: with that line's median mocked from 5 to 60, MATLAB smooths frame 17 with medfilt1 and frame 50 with conv at every value, frame 51 with conv only at 5, so the choice does not follow that median alone, and no rule of one quantity against it fits the 1,386 recorded and designed calls, the best 1,259; B's median bound on 21: 103 candidates, median count 8, MATLAB takes it with U 9 or 10); matched values: the master's positions within 3.4e-10 relative, A 9.1e-5, the other channel's A 0.27 (a few refits). Re-measured after the filter's few-point rules and line 160's bound (C02c's follow-up, job 4636924): the same counts on every frame (camera threshold 3, about 1,000 points a frame: none of those rules changes these calls). Re-measured with C02d's rules (the runners' residual threshold, the fits' box and stop test, the standard deviations of a negative diagonal, line 104's dip; job 4637860): 59 frames with MATLAB's number of detections, 42 of them with every detection matched; of MATLAB's 61,260 detections 61,239 matched, the port has 61,263; frames 1, 17, 33, 50, 51 and 55 now as MATLAB's, frame 21 (B's median bound) 982 against 979. Re-measured with C08's stop test (job 4638877): the same 59 frames with MATLAB's number of detections, 54 of them with every detection matched (42 before); of MATLAB's 61,260 detections 61,255 matched (61,239 before), the port has 61,263; frame 21 as before | none on the oracle cases (bit for bit); demo frames: listed | D3 (a restarted fit that moved 2.9 voxels and stopped 2.9e-4 away on frame 1, replayed from MATLAB's own window and start), which the filter's histogram and neighbour counts pass on; the filter's smoothing choice and B's bound (measured on the frames with the background correction) |
| D11 | `run_detection_frame_parallel(..., overwrite=True)` when `Detection3D.mat` or per-frame files exist | MATLAB keeps them (its job framework counts existing outputs as done, so `Overwrite` true changes nothing: probes fp1-fp3) | none: the port removes them first and detects every frame again, as `XR_runDetection3D` and the lab's other wrappers do with `Overwrite` | MATLAB's bug (an option ignored), fixed in the port |
| D12 | a frame without detections in `run_detection_frame_parallel` | MATLAB's entry repeats the session's previous frame's detections (with the frame's own `dRange`; probe fp3) | none: the port gives the empty entry, as `XR_runDetection3D` does | MATLAB's bug (a wrong result), fixed in the port |

The frame-parallel runner (`detection.run.run_detection_frame_parallel`, measured on probes fp1 to
fp4 of the synthetic condition): MATLAB runs its frames through its job framework, one task per
frame whose output is `[ResultsPath 'Detections' filesep 'Detection_F_<k>.mat']` (k not padded;
one variable `frameInfo_k`, the frame's entry; a marker `Detection_F_<k>.tmp` while it runs),
with `Detection3D.mat` as the final output, then gathers the files into `Detection3D.mat`; the
masks, the Amira folder and the skipping (`Overwrite` false and `Detection3D.mat`) as
`XR_runDetection3D`'s; a frame whose file exists is taken as it is. The port runs the same tasks
through `petakit5d.jobs` (local by default, as MATLAB's session; a process pool or Slurm with the
framework's options). MATLAB's default `DetectionMethod` runs `XR_pointSourceDetection3D_1`,
'original' `pointSourceDetection3D`, 'adaptive' `XR_pointSourceDetection3D_4`: not ported.
Measured too: a foreign fresh marker made MATLAB's runner leave its frame and end after about 14
s without `Detection3D.mat`; the port's framework waits for the marker's age (3 minutes) and
then takes the frame over (the framework's rule).

## The bad-point filter on few points, line 160's bound, histcounts of one value (C02c)

Measured the same way (designed and mocked calls of MATLAB's filter, line numbers and keep masks
only, no source read; 118 calls on few points, 318 designed histograms for line 160):

- **Fewer than three points:** MATLAB returns one `true` whatever the points (one and two, NaN
  among them; no MathWorks call), and the detector keeps them all: the port does the same
  (`detect.postfit` keeps every fit for a filter result of one element).
- **The short path (lines 34-39):** when `prctile(sr, 95) <= 2 prctile(sr, 5)` (taken up to the
  ratio 2 itself, not from 2.0000000003, whatever the spread or the count), lines 40-279 are
  skipped and the points are removed as on the path A with thr (lines 273 and 318), but a point
  with three neighbours only below `min(thr, cam, 2 noise)` (kept below thr 4 with camera
  thresholds 3 to 5.5 and noise 1 or 2.6, removed below thr 3 with camera 3 and noises 1.5 to 4;
  the long path removes it below thr, camera 2 included); a NaN sigma counts as below every
  threshold there (removed with 1 to 8 neighbours, kept with more).
- **The long path's early return:** with fewer than 10 sigmas below their 99th percentile (5 to
  10 points measured; 11 went on) every point is kept (lines 58-60 to 115), with no neighbour
  count.
- **NaN on the long path:** left out of the percentiles, the histogram and the bin count
  (MATLAB's `max` and `min` leave out the NaN that `std` gives), and never removed for its sigma
  (five NaN in channel 1 frame 28's call, each point given two neighbours: kept).
- **Line 160's bound** (C, the first valley's early end, when the background peak's bin centre is
  at or above it): `2 max(cam, noise)` while the noise is at most 2 camera thresholds (to 0.02 at
  noises 2.4 to 3.5 and camera thresholds 1.5 to 6; thr 2, 3 and 4 alike), `2.5 noise` above
  that (between 8.02 and 9.02 at noise 3.5 and camera 1.5, between 13.02 and 15.02 at 6 and 2.5),
  and no C with a camera threshold of 1.25 or less and the noise above it (A up to the largest
  centres designed: 5.6 at 1.25, 33 at 1; where between 1.25 and 1.5 that ends is open). C02's
  rule (`2 cam`) gave C there from 2 cam: on the 5,793 recorded and designed calls that reach
  line 160 the bound gives MATLAB's label on 5,774 (C02's 5,530). Open: 11 recorded variants
  whose sigmas' 30th percentile lies in a far mass of counts (C at a centre below the bound;
  they take the branch of lines 154-155), 3 designed histograms whose bin 1 holds one count (A,
  line 146) and 5 designs with the sigmas scaled at camera threshold 1 (C from 8.8; line 153).
- **Checks:** the demo's 120 frames with the background correction stay bit for bit; 248 of the
  285 replayed variants (one lost: noise 6.9 at camera threshold 3, where line 160 now takes A as
  MATLAB does and the B path's removal with the noise above 2 camera thresholds is not ported, 4
  points); 955 of the 1,011 designed decisions (942 before) and all 580 designed keep masks with
  mocked neighbour counts (567 before); the 118 designed calls of few points.
- **histcounts of one value** (`matlabcompat.signal`; equal on 640 of the 712 calls of probes
  hc1 and hc2, the other 72 refused): the limits are `v - 0.5` and `v + 0.5`, each to the
  nearest multiple of 0.5 with a tie toward `v`, split into equal bins; non-finite values are
  neither counted nor set the edges; with no finite value the edges are 0 to nbins. Bins
  narrower than about 1000 ulp of the value (`1e15 + 0.25` with 13 bins or more) are refused:
  MATLAB widens its limits there, not measured.

## The runners' residual threshold, the fits at the box, line 104's smoothing (C02d)

Measured the same way (traced calls with full values, line numbers, designed and recorded
inputs; no GPL or MathWorks source read), on the multiple-experiments demo's frames and
synthetic ones. These were the causes of the demos' differences in detection (phase3.md's rows
2, 8 and 10):

- **The runners' residual threshold:** `XR_runDetection3D` and
  `XR_runDetection3D_frame_parallel` give the detector a `ResidualSigmaThresh` of their own per
  frame, from the frame's noise (`BackgroundInfo(3)`, the third value of
  `XR_estimate_overall_background_information`) and the camera threshold: 0.9 camera
  thresholds while the noise is below 1/1.5 of it (the bound measured between 0.662 and 0.682
  camera thresholds), the camera threshold up to it, the noise up to 1.5 camera thresholds
  (between 1.477 and 1.506), the camera threshold beyond (92 frames of the demo and of
  synthetic conditions, camera thresholds 3 and 65, both runners; the values at exactly the
  bounds are not measured). The port gave the camera threshold always; now
  `run.residual_threshold`. Without a camera threshold (no settings file above the source, none
  with a camera line: the camera search's errors, which `XR_runDetection3D` catches; or a camera
  it does not know) `XR_runDetection3D` takes the noise, with NaN as `BackgroundInfo(4)`; the
  frame-parallel runner: D13.
- **The fits' box** (`XR_fitGaussian3D`'s CD_Newton, the port's `fit.fit_cd`): the positions are
  held from the window's first voxel to its last valid one on each axis, so a window whose last
  planes are NaN (the deskewed frame's empty wedge) keeps its fit off them, while NaN first
  planes do not move the box's start (the port held the fits in the whole window). A fit held
  on the box (its step cut to nothing, no start left to switch to) goes on while its amplitude
  and background still move: the stop test's step is the direction it was held from (the
  port's was the position's move, 0 there). On 23,969 traced fits (8 frames of the first two
  experiments and C02's probe 14) the port's results within 1e-12 of MATLAB's: 23,171 before,
  23,510 after (beyond 1e-3: 79, 69). MATLAB's cut steps land an ulp inside the box at times and
  take one more step there before they switch; the port puts a position within 1e-12 of the box
  on it (with positions put on the box only from beyond it: 23,376; without either: 23,273).
- **The standard deviations of a fit with a negative entry on its covariance's diagonal** (a fit
  that ends on the box where the RSS is not convex): MATLAB gives the magnitude of the complex
  root, `sqrt(|C_kk| RSS / (n - 6))`; the port gave NaN, which made the amplitude test's
  p-value NaN and dropped the detection (row 2's two detections). With the magnitude the port's
  `hval_Ar` is MATLAB's on all 14,673 fits evaluated at MATLAB's final positions.
- **Line 104's smoothing choice** (`badpoints.spike`, `badpoints.dip`): conv when the largest of
  the first three bins is a spike above the median of the list's counts (both neighbours below
  0.75 of it), or when medfilt1 raises one of the first three bins by more than that median;
  else medfilt1. MATLAB's choice on 1,456 of the 1,479 recorded, designed and mocked calls that
  reach the line, all 64 of the demo frames' calls replayed (the spike rule alone 1,453, C02's
  rule 1,168); the misses: 21 designed histograms with a spike in bin 3 (20 of C02c's grid, where
  MATLAB takes conv only when bin 1 is above bin 2), and 2 of C02c's mocked medians. The rule is measured, not
  known. Its oracle cases (`demo.smoothing.*`, with the demo data) give both sides the
  filter's inputs the port's detector makes of the second experiment's frame 16, which the spike
  decides (23 points apart from C02's rule's keep mask), and of the third's frames 27 and 40,
  which the dip decides (5 and 3 points apart from the spike rule's): MATLAB's keep masks.
- **The results:** the multiple demo's 180 frames (the master channel on MATLAB's DS frames with
  MATLAB's sigmas): the frames whose detections differ, 2, 21 and 10 per experiment before, 0,
  3 and 1 after (the first experiment: 2 detections only MATLAB has before, none after; the
  second: 12 and 102 before, 2 and 1 after; the third: 7 and 28 before, 1 and 0 after); the 4
  left are fits that take another path (D3). D10's frames: 59 of 60 with MATLAB's detections
  (53 before).

| Row | What | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| D13 | `XR_runDetection3D_frame_parallel` (`XR_runDetTrack3d`'s default detection on Linux and macOS) without a camera threshold | per case (probe d13: a condition of three frames, both runners with the demos' options, the detector's calls traced; cases of `tests/detection/test_runner_camera.py`): (a) no settings file above the source: the camera search stops with `MATLAB:subsdeldimmismatch`; `XR_runDetection3D` catches it and gives each frame's detector the frame's noise (4, 3 and 3 detections); the frame-parallel runner stops with the error before any frame and writes nothing, and `XR_runDetTrack3d` stops with it; (b) a settings file without a camera line: `MATLAB:TooManyOutputsDueToMissingBraces`, the same; (c) a camera the search does not know, or an empty camera name: the threshold NaN; `XR_runDetection3D` takes the noise; the frame-parallel runner gives each frame's detector NaN and finds nothing (0, 0, 0) without a message, writing the per-frame files, `Detection3D.mat` and the masks; `XR_runDetTrack3d` then prints "Empty movie. Nothing to track." and leaves only `.tmp` markers of its tracks, and with `'ParallelMethod', 'parfor'` (its default on Windows) it detects and tracks; (d) a threshold of NaN (the threshold function stubbed): as (c); (e) a threshold of 0 (stubbed): both runners give the detector 0, which keeps every candidate (4, 4 and 4) | none: in (a) to (d) the port's frame-parallel runner does what `XR_runDetection3D` does (its files are MATLAB's `XR_runDetection3D`'s, cases `camera.*.frame_parallel`); in (e) both of the port's runners are MATLAB's (`residual_sigma_threshold=0`, cases `camera.stub_zero.*`) | MATLAB's bug in (a) to (d), fixed in the port. The lab's camera search means to give no camera there: its code returns no camera class when it finds no settings file or no camera line, but fails before those returns in (a) and (b); the lab's threshold function says that without camera information the threshold is to be estimated from the image. The user's fix of 2026-10-21 made `XR_runDetection3D` stop "dropping all detections without a camera noise threshold" (the commit's words); before it, `XR_runDetection3D` found 0, 0 and 0 detections in (a) to (d) too (measured on its previous version), as the frame-parallel runner still does; and `XR_runDetTrack3d` takes one runner or the other by its parallel method alone. (e): MATLAB's, kept |

## The fits' stop test (C08)

Measured the same way (traced calls with full values, designed and recorded inputs; no GPL or
MathWorks source read), on track processing's recorded calls of `XR_fitGaussian3D` in the
multiple-experiments demo (C04c's: the 4,320 fits of the 239 tracks whose fits differed, every
window and start as MATLAB had them; track processing calls it without the 7th argument that
`XR_fitGaussians3D` passes, 10, and MATLAB's results are the same with it on all 4,320) and on the
detector's traced calls (C02d's of 8 demo frames and C02's probe 14: 23,969 fits):

- **The rule:** the CD_Newton fit stops after an iteration in which the RSS dropped by less than
  1e-7 of the first RSS with a step below 1e-4 (the position's largest move; held on the box, the
  direction it was held from), or in which the RSS changed by less than 1e-6 in absolute value.
  C02 had the second bound as 1e-10 of the first RSS, which is 1e-6 at the detector's windows'
  RSS of about 1e4; track processing's windows of weak spots have RSS of 3,000 to 7,000, where
  the port's fits took one iteration more than MATLAB's (290 of the 4,320, after C02d's box and
  stop test), and brighter ones of up to 1.7e5, where they took one fewer.
- **The iteration limit:** a fit that has not stopped ends after 2,000 iterations (the port had
  500): traced, MATLAB stops 3 of 50 fits of the recorded windows at 1/64 of their values at
  exactly 2,000 iterations, and of 12 designs that the port now runs to 2,000, 9 end where
  MATLAB's do (probe d1). Only values far below the demo's creep that long: no fit of the demo
  data takes 500 iterations (at most 150 and 263 in the two sets), and the three demos give the
  same values with either limit, end to end and each step alone (phase3.md's runs).
- **How it was found:**
  - MATLAB's iterations traced (136 of the recorded calls; C02's tracer on the gradient/Hessian
    MEX, which MATLAB calls once per iteration, each state decoded exactly from its inputs): in
    the 60 fits that MATLAB ended an iteration before the port (56 with every state decoded), its
    last iteration changed the RSS by 1.0e-10 to 2.6e-10 of the first RSS, above C02's bound,
    with steps of 1.0e-4 to 6.0e-4, and by 3.4e-7 to 9.97e-7; in the 107 traced fits with every
    state decoded, none of the 1,326 iterations before the last changed it by less than 1e-6
    (the smallest 1.03e-6).
  - The bound scanned, both sets together: with 0.99e-6 or 1.01e-6 in its place 10 and 2 fewer
    fits end at MATLAB's, with 0.95e-6 or 1.05e-6 43 and 22, with the signed drop in place of its
    absolute value 6,136, with both bounds (this one and 1e-10 of the first RSS) 281. The first
    test's 1e-7 and 1e-4 checked the same way: 5e-8, 2e-7, 5e-5, 2e-4 and an absolute drop of
    1e-3 in their places all give fewer.
  - Windows scaled by powers of 2 (the window and the start's A and c times 1/64 to 64: every
    residual and RSS of the fit scales exactly): MATLAB's ends move with the scale (at twice and
    a quarter of the values, 2,388 and 1,517 of the 4,320 end at the unscaled positions), so the
    bound is absolute, not relative to the RSS; with the rule the port's ends are MATLAB's at
    scales 1/4 to 64 in 4,239 to 4,263 of the 4,320 (the base's rule: 3,281 to 4,063).
- **The results** (within 1e-12 of MATLAB's parameters, relative to max(\|v\|, 1)): track
  processing's 4,320 recorded fits 4,256 (3,890 before); C04c's 465 reproducers 408 (42); the
  detector's 23,969 fits 23,737 (23,510); no fit moved away from MATLAB's end. The demos
  (phase3.md): the step-by-step demo's detections each step alone, the fits beyond 1e-9 in the
  two channels, 6 and 16 (17 and 48: D3's list); W13's windows on MATLAB's inputs beyond W11, 2
  of the step-by-step demo's (24) and 2, 13 and 5 of the multiple demo's three experiments (24,
  186 and 36); D10's frames with every detection matched, 54 (42).
- **What is left:** 64 of track processing's fits and 232 of the detector's. Traced again (the 64
  of track processing), MATLAB's path and the port's part in 56 where a position lies an ulp
  inside the box on one side and on it on the other (MATLAB's x 1.1102230246251565e-16 where the
  port's is 0, for example), after paths whose last bits already differ by about 1e-14 from a
  damped or restarted step (D3), and in 8 where those last bits grow apart without the box:
  MATLAB updates A and c once more on the face before it switches to its next start, and ends
  an iteration behind with another best record. This is
  C02d's open point (the port puts a position within 1e-12 of the box on it); without that, or
  with positions put on the box only from beyond it, 83 and 53 more of track processing's fits
  end elsewhere. Decoded fit by fit, with the last bits' first cause (MATLAB's update of A) and
  the fits of the demos' rows: "The fits still on another path, fit by fit (C09)" below.
- **Measured, not reproduced:**
  - *A negative starting amplitude:* MATLAB switches to its next start from A of about -1e-3
    down at the demo's scale (5 of 6 recorded windows; -1e-4 in the sixth) and takes a step for
    A of -1e-4 and -1e-6; at 1/64 of the values it takes steps down to -0.003 to -0.03 and
    switches from -0.01 to -0.1 (probe a2: 6 windows, 11 amplitudes, both scales). The bound
    grows as the values shrink, so it is neither a bound on A nor one relative to the values; the
    port switches at any A below 0. Of the recorded fits one starts with A in (-1e-3, 0) and 7
    in (-1e-2, 0), all with MATLAB's end.
  - *Values far below the demo's:* at 1/16 and 1/64 of the recorded windows' values the port's
    ends are MATLAB's in 3,961 and 3,238 of the 4,320 (3,866 and 2,753 with a limit of 500
    iterations, 2,201 and 991 before both rules); the paths part where MATLAB keeps a slightly
    negative A and at damped steps.
- **The cases** (`tests/detection/test_fit.py`, `xr_fit.stop*` and `xr_fit.demo_*`): three
  designed weak spots, each at a quarter, once and four times its values (the 1e-6 bound stops
  them at their scale an iteration before the bound of 1e-10 of the first RSS would, and later
  at 4 times), and 6 of track processing's recorded windows, cut from the demo's frames as the
  port deskews them (`demodata`; windows of gap frames and buffers, both channels, three that
  the bound stops an iteration earlier than C02's and three an iteration later): MATLAB's results
  within D3's 1e-9. Three designs at 1/64 and 1/256 of their values (`xr_fit.iters*`), whose
  fits take 1,104, 2,000 (the limit) and 1,232 iterations.

## The fits still on another path, fit by fit (C09)

Measured the same way (C02's tracer on the gradient/Hessian MEX, each iteration's state decoded
exactly, and MATLAB's executed line numbers only; no GPL or MathWorks source read), on every fit
that still ends beyond 1e-12 of MATLAB's: track processing's 64 of the 4,320 recorded fits (C08's
traces, 142 traced fits ending at MATLAB's beside them), with 19 of W13's 20 windows among them (the
20th cut from MATLAB's inputs as the port reads it, and traced); the detector's 232 of its 23,969
(traced again, with 30 that end at MATLAB's); D3's list and the windows of row 10's detections
([phase3.md](phase3.md)), cut from MATLAB's DS frames as the port cuts them and traced. The port's
iterates come from a mirror of `fit_cd` in Python, bit for bit the port's (x, y, z, A, c, the RSS
and the iteration count of all 28,289 recorded fits).

- **From MATLAB's state, MATLAB's step is the port's.** At each of 10,280 decoded MATLAB iterations
  the port's gradient and Hessian are MATLAB's bit for bit, and so is a Newton step's next position
  in 2,072 of the 2,152 off the box (the others' directions a few ulps apart; a Cholesky solve gives
  71%, the port's LU 96%). Where the step is cut to the box, MATLAB's next position is the cut `x +
  alpha d` as the port computes it, not snapped onto the face: in all 286 cut steps where the port's
  snap moves the position (MATLAB's lie 1.2e-32 to 2.2e-16 inside the face, or up to 1.1e-16 beyond
  it), the unsnapped cut is MATLAB's position; in 24 more MATLAB's is on the face where the port's
  cut from MATLAB's state is an ulp off it (damped steps: their directions differ in the last bits).
  From beyond the face with its direction outward (a cut length below 0) MATLAB switches to its next
  start, as from the face.
- **The paths' last bits part first, before the box.** MATLAB's update of A and c at the new
  position is the port's in value, not in its last bits: the port's A is MATLAB's in 17% of 3,000
  traced iterations and within 5 ulps in 75%, and none of 20 orders and forms tried (sums in index
  order, in 2 to 16 lanes, with fused multiply-adds, pairwise, a QR least-squares solve, the
  expanded sums) gives MATLAB's A in more than 31% of them (its c: 72%). Of the 296 fits that end
  elsewhere, 272 have A or c other while their positions are still equal, after which their steps
  part; 14 part at a damped step (the eigenvalues), 9 at a damped step that the box cuts, and 1 at
  the box itself, from states equal in every bit.
- **The decisive parting** (the positions beyond 1e-6) **is at the box**: in 56 of track
  processing's 64 fits, 229 of the detector's 232 and 8 of the other 9 (D3's list, row 10), the
  paths stay within 1e-6 until a cut step lands on a face where MATLAB's position (its own unsnapped
  cut, from its other last bits) lies an ulp inside and the port's on the face (its own cut landing
  there, or its snap). MATLAB then takes one more step on the face, which updates A and c once more,
  and switches to its next start an iteration later, or goes back into the window where that update
  has turned its direction (row 10's first detection below); the paths run an iteration apart from
  there, and their best records and stops differ. C08's reading is so, fit by fit. The others: 3 of
  track processing's fits and 1 of the detector's part at a switch at a slightly negative amplitude
  (A -4.0e-5, -3.3e-4 and -3.6e-4, where MATLAB takes a step and the port switches: the bound
  measured above, not reproduced), and 8 end apart within 1e-6.
- **MATLAB's rule at the box cannot be followed without moving fits away.** With the cut not snapped
  and a cut length of 0 or below holding the fit, the mirror gives MATLAB's end in 4,204 of track
  processing's 4,320 recorded fits (4,256 with the snap: none moves to MATLAB's end, 52 away) and in
  23,595 of the detector's 23,969 (23,737: 28 to, 170 away); unsnapped alone, 4,173 and 23,487. The
  port's paths reach the faces with other last bits than MATLAB's, so the side MATLAB's cut lands on
  is not the port's to follow there; the snap gives the most fits MATLAB's ends and stays. Following
  MATLAB at the box would need its update of A to the bit.
- **W13's windows**, one by one: the step-by-step demo's 2 (the multiple demo's first experiment's)
  part decisively at the box (a start buffer, 0.018 voxels, after a damped first step cut by the
  box; a gap frame, 1.6e-7, after A's bits); the second experiment's 13: 9 at the box, a channel 2
  end buffer (frame 58, 6.1e-6 voxels) at a switch at A -3.6e-4, and the three on either side of the
  keep rule: channel 2's end buffer at frame 33 (track 3414, 1.77 voxels) at a switch at A -3.3e-4,
  channel 2's gap frame 43 (track 1585, 1.33 voxels) at a switch at A -4.0e-5, and channel 1's end
  buffer at frame 30 (track 3526, 2.27 voxels) not on another path at all: MATLAB's free fit of the
  port's window and start ends at the port's (z 5.0 exactly, the window's last valid plane; x an ulp
  apart), exactly 2 voxels from the window's centre, where MATLAB refits at the given position and
  the port kept its free fit (`processing.fit_windows`: `<= 2`). That is the only one of the
  multiple demo's 171,926 gap and buffer fits at exactly 2 voxels: MATLAB's bound is `< 2` (measured
  on that one). The port's bound is `< 2` since 2026-11-22 (C10; the oracle case `run.keep_bound`
  of `tests/detection/test_processing.py` refits 12 designed windows at exactly 2 voxels, as
  MATLAB does): that window is MATLAB's, which leaves 12 of the second experiment's windows beyond
  W11. The third experiment's 5: 4 at the box, 1 apart within 1e-6.
- **D3's 22** (the step-by-step demo's 60 frames, each step alone): of the master channel's 6, 5
  part at the box (from the same window and start), and one gives the port's end from the port's
  start (1.6e-14), so its own start differs (D2's last bits of the start's A and c); of the other
  channel's 16, 13 give the port's end from the port's window and start (within 1.3e-15), so their
  `'xyzAc'` starts differ (the `'Ac'` fits' A and c, D4), and 3 part at the box.
- **Row 10's 4 detections:** the second experiment's frame 5 (only MATLAB, A 24.5 on the first
  plane): MATLAB's fit of a window the port cuts too ends at its detection ([3.070, 0.0027, 0.0] in
  the window), the port's at [2.930, 0.0, 0.343] with A 20.88: a cut lands on a face at iteration 3
  with A 2 ulps apart, MATLAB's 1.1e-16 inside, and its one more update of A (20.88 to 21.72) turns
  its direction into the window, to its detection, while the port, held on the face, switches. Frame
  8 (only MATLAB, A 29.4) and the third experiment's frame 30 (only MATLAB, A 10.0): the paths agree
  but in A's last bits until a cut at iteration 26 and 22 (MATLAB 3.5e-18 and 1.1e-19 inside; the
  port on the face, the second time by its snap); one iteration apart, the port's best record is
  another when both stop. Frame 6 (only the port, A 1,380): not on another path: MATLAB's fit of the
  port's window and start ends at the port's end (4.9e-16), and MATLAB leaves the detection out at a
  later stage of its detector (not measured further).

## Tracking (`detection.tracking`, C03)

The tracker is a new implementation of the method of Jaqaman et al. (Nature Methods 5, 695;
2008), with the details measured from MATLAB's tracker as a black box (no u-track source used).
MATLAB's mex files without source are replaced: `lap` (`mexLap`) by SciPy's
`min_weight_full_bipartite_matching`, `KDTreeBallQuery` by SciPy's `cKDTree`, and
`createDistanceMatrix` by NumPy.

Exact, bit for bit: the whole `trackedFeatures.mat` of the demo's detections (60 frames, 61,353
detections, 10,917 compound tracks: `tracksFinal`, `kalmanInfoLink` and the saved settings), and
every case of `tests/detection/test_tracking.py`: the tracks and filter states of seven
synthetic movies (spots in Brownian motion with births, deaths and missing frames; a dense
field; merges and splits; PointDetectionTest's movie; frames with no or one detection; the
search radius rules; the rules of gap closing and merging), `XR_runTracking3D`'s file,
`loadTrackSettings`, `createDistanceMatrix`, `KDTreeBallQuery`, and `lap` with ties and zero
costs. No case declares a tolerance. The rows below are where the two could still differ; where
the cases reach them (`lap`'s ties), the results are still MATLAB's.

| Row | Where | Measured | Tolerance | Cause |
|---|---|---|---|---|
| TK1 | `lap`: assignment problems with more than one optimal assignment | the same assignment as `mexLap` in the tie cases, in 60 of 60 tie probes and in 169 of the demo's 178 problems (3 linking passes of 59 frame pairs, and gap closing); in the other 9 only the pairing of birth and death entries differs (same links, same tracks) | none on the outputs. A case whose costs tie between links compares the total cost of the assignment (`rel` 1e-12) and the number of links instead of the tracks | SciPy's solver and `mexLap` may break ties differently |
| TK2 | `KDTreeBallQuery`: the order of the points found at exactly the same distance from a query point | not reproduced in general; the cases (random points, points on the radius, 2-D self query) bit for bit | none. A case with such ties compares the points at each distance as a set | MATLAB's order follows its tree, the port's the order of the points (a stable sort by distance). The tracker does not depend on it |
| TK3 | the gap-closing, merging and splitting costs (inside the tracker, neither returned nor saved) | the same 24,349 candidates as MATLAB's on the demo; the costs within 9 ulp (bit for bit: 82% of the gaps, 69% of the merges, 66% of the splits) | none on the outputs: all tracks equal on the demo and the synthetic movies | the order of MATLAB's arithmetic in these costs was not recovered from black-box probes. It can change a track only where two solutions' total costs are within a few ulp |
| TK4 | detections at exactly x = 0 | MATLAB never links them: a spot at x = 0 in every frame has no track, and one at x = 0 in one frame has a gap there (probes; y = 0, z = 0, x = 1e-300 and negative x are tracked) | none: the port tracks them like any other. Detections are pixel coordinates from 1 (the demo's 61,353 all have x from 45.6 on), and the cases keep away from x = 0 | MATLAB takes x = 0 for no detection |

`trackedFeatures.mat` holds the same values as MATLAB's. MATLAB's `save` stores integer-valued
doubles, such as the feature indices `tracksFeatIndxCG`, as 8- or 16-bit integers inside the
file; the port stores them as doubles. Both load as double in MATLAB and with
`petakit5d.io.matfile.read_mat`; `scipy.io.loadmat` shows the stored types unless it is given
`mat_dtype=True`.

## Condition loading, deskew and channel offsets (C04a)

Bit for bit, in every oracle case (synthetic conditions and the demo dataset): the condition
structs (`detection.condition`; `date` checked against MATLAB's probed rule, since it is the
digits of each side's own folder path), the deskewed frames with `settings.mat` and
`dynRange.mat` (`detection.wrapper.deskew_data`), the channel offsets, `xz_offset.mat` and the
padded frames (`detection.offsets`), and `interp1` (`matlabcompat.interp`). The cases' `reason`
names its row as `detection.md Wn`. Approved by the user on 2026-11-20.

| Row | Function and range | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| W1 | `corr` (Pearson, two columns): the 165 correlations of `XZoffsetEstimation3D` and the probed columns | 17 ulp where the sum cancels (small integer columns); 16 ulp at 65,536 values; 2 ulp at 5.7 million; single: 1 ulp of single | `abs` 4e-15 (single `abs` 2e-7) | MATLAB takes `x'*y / norm(x) / norm(y)` from its BLAS, whose order of additions depends on the CPU and, from 4,096 values, on its thread count (MATLAB's own results differ between 1 and 8 threads); the port sums in MATLAB's order, platform-independent. The offsets the correlations decide were equal in every case |

### MATLAB's slips the port fixes (W2 to W4)

Three places where MATLAB's result is wrong and the port gives the correct one instead, each
measured on synthetic data (probe job 4634531, oracle tree 45b2bb1).
- The oracle cases compare the port with the corrected MATLAB call.
- `run_detection_tracking(..., matlab_quirks=True)` gives MATLAB's output for W2 and W3;
  `deskew_data(..., reverse=False)`, deskewData's default, gives it for one deskew.
- Approved by the user on 2026-11-20.

| Row | Where | MATLAB | The port | Evidence |
|---|---|---|---|---|
| W2 | `run_detection_tracking`'s deskew of a reversed scan (a condition loaded with `reversed=True`) | `XR_runDetTrack3d` calls `deskewData` without its `'Reverse'` option, and `deskewData` leaves `data.reversed` alone, so reversed scans are deskewed as forward ones: sheared the wrong way | each movie is deskewed in its direction (`deskew_data(reverse=None)`), bit for bit `deskewData` with `'Reverse', true`, which is `XR_deskewRotateFrame` with `reverse` true | six beads at depths 3 to 19 of a 21-plane scan (31.5 degrees, dz 0.5, 0.104 um pixels). With `'Reverse'` true, every bead's x is the forward scan's (0.00 px apart; z mirrored, the acquisition order). Without it, the beads are displaced from -34 px to -132 px (13.7 um) relative to the first, by twice the shift per plane times their depth. Cases `deskew.reverse_option`, `deskew.wrapper_reversed_scan`; MATLAB's output in `deskew.reversed` |
| W3 | `run_detection_tracking` with a one-frame movie before others | it drops the one-frame movies from the tracking (tracking needs two frames) but passes `XR_runTracking3D`, `XR_runTrackProcessing3D` and `rotateTracks3D` the folders of all movies. Its i-th kept movie gets the i-th folder: an indexing slip | each movie of more than one frame is tracked, processed and rotated into its own folder | movies of 6, 1 and 6 frames. MATLAB wrote a one-frame tracking (from the one-frame movie's 9 detections) into the one-frame movie's folder, left the third movie without tracks, and failed in track processing (`MATLAB:numel:BadSubscriptingIndex`, a `ProcessedTracks.mat.tmp` left behind). Case `dettrack.one_frame_movie` (against `XR_runDetTrack3d` on the movies of more than one frame, plus the one-frame movie with `Track` false) |
| W4 | `load_condition_data(..., frame_rate=...)` with a one-frame movie | fails (`MATLAB:cellRefFromNonCell`, its dz lookup indexes the frame paths), even with the movie alone | loads it with the given rate: the struct of MATLAB's load without the rate, its `framerate` set | case `load.rate_given_one_frame_movie` (MATLAB's error kept in `load.rate_given_one_frame_movie.matlab`) |

Two more measured behaviours:
- `XR_loadConditionData3D` asks whether unequal frame rates are right ("Unequal frame rates
  detected: 2.000 NaN 2.000", then `input`). It asks for movies whose rates differ, and for a
  one-frame movie's NaN among others, two NaNs included; not for one movie alone. In batch mode
  that is `MATLAB:services:MissingRequiredCapability`.
  - The port raises the same error (`ConditionError`): cases `load.err_rates_*`.
  - A frame rate given (`frame_rate`) makes the rates equal, as in MATLAB.
  - A condition with a one-frame movie among others is then loaded one movie at a time and
    joined (`Condition(...)`).
- `deskewData` takes `'Reverse'` (and `'reverse'`) but refuses `'Reversed'`
  (`MATLAB:InputParser:UnmatchedParameter`).

### MATLAB's failures the port does not copy (W5, W6; proposed 2026-11-22, for the user)

Two calls that MATLAB fails on and the port carries out, measured on synthetic conditions
(probe jobs 4640495, 4640506 and 4640539, oracle tree 45b2bb1). `deskew_data(...,
matlab_quirks=True)` (`'matlab_quirks', true` in the MATLAB-name layer) gives MATLAB's failures.

| Row | Where | MATLAB | The port | Evidence |
|---|---|---|---|---|
| W5 | `deskew_data` with `aname` left out | `deskewData`'s default `aname` is the handle `@isstr` (a validator in the place of a default), which it cannot join to a folder name: every call that leaves `aname` out fails at its first movie, before it writes anything (`MATLAB:UnableToConvert`; one movie or two, and a movie it would skip). `XR_runDetTrack3d` and the step-by-step demo pass `aname` (`XR_runDetTrack3d` its own, `[]` by default), so they never meet it | the default `aname` is `None`, MATLAB's `[]`: `settings.mat` and `dynRange.mat` in each movie's source folder. With `matlab_quirks` an `aname` left out fails as MATLAB's (`''` is then the source folder) | case `deskew.default_aname` (the port's default against MATLAB's `'aname', []`); MATLAB's failure in `deskew.no_aname` (the MATLAB-name layer), `deskew.no_options` and `deskew.no_aname_skips`, compared with `matlab_quirks` |
| W6 | `deskew_data` of a condition without movies | fails whatever its options, `aname` given or not (`MATLAB:Java:InvalidInput`) | returns the condition unchanged; with `matlab_quirks` MATLAB's failure | case `deskew.no_movies` (with `matlab_quirks`) |

## Track processing (`detection.processing`, `detection.amira`, C04b)

`XR_runTrackProcessing3D` and `rotateTracks3D` are new implementations from MATLAB's behaviour as
a black box (both headers name an outside author): logging wrappers with full values in front of
the functions they call, the profiler's executed line numbers and counts, 717 tracks in designed
probes, and the files MATLAB writes; no GPL or MathWorks source read. `XR_amiraWriteTracks` is
the lab's and was read. The gap and buffer fits are `detection.fit`'s: `fit_cd` (row D3) for
the free fits, the `fitGaussian3D` mex's mode `'Ac'` (row D4) for the refits at the given
positions; the interpolations are `matlabcompat.interp.interp1`'s (bit for bit: C04a's section).
Approved by the user on 2026-11-20.

Bit for bit on the demo (10,917 tracker tracks, 3,636 processed tracks): the set and order of the
tracks (single-frame and border tracks dropped, compound tracks decoupled, the stable sort by
lifetime), every field at the detected frames, `tracksFeatIndxCG`, `gapVect`, `gapStatus`, `gapIdx`,
`seqOfEvents`, `nSeg`, `visibility`, the times, the read-out windows and the fits' starts (42,412 of
42,412), the keep-or-refit decision of each fit, `hasMultiplePeaks`, the three significance fields,
the lifetime cohorts' thresholds, `processingInfo`, every field's class and size (a buffer's
`hval_AD` double or logical, as MATLAB's) and the field order; `hval_AD` at the gap frames and in
the buffers (3,980 of 3,980 buffers); `totalDisplacement` (3,331 of 3,331). Given MATLAB's processed
tracks: the rotation (323,568 positions, 7 probes of 4 geometries, and `RotatedTracks.mat` made from
MATLAB's `ProcessedTracks.mat`: 3,636 of 3,636 tracks, every field) and the Amira files (60 of 60
frames). The oracle cases `run.designed` (gaps of one and two frames with and without signal,
buffers with signal, multiple peaks, cut, two-frame, single-frame and border tracks, a split and a
merge, large steps) and `run.compound` (34 compound structures: branches that split off and die, are
born and merge, end or start with the main segment, continue it, merge into a segment born 1 to 3
frames before, loops, two branches, gaps, a join whose overlap frame needs the interpolation, the
same join with a branch split off the merging segment, a structure over the whole movie, the movie's
first and last frames): every field but those of W11 and W12, in both files.

The multiple-experiments demo's three experiments, each processed on MATLAB's detections, tracks,
masks and frames (C04c, 2026-11-21): the same tracks (3,636, 6,553 and 5,193), every category
(15,382 of 15,382, since the buffer test's intensity rule of W14), the keep-or-refit decision of
171,924 of the 171,926 gap and buffer fits (the two others W13's; 171,923 before the keep rule's
bound of 2026-11-22, below), and, recorded from MATLAB on
the 239 tracks whose fits differ, every read-out window (values and the masked voxels) and every
start of those tracks' 4,320 fits bit for bit. The oracle cases `run.buffers` (spots in the
farthest buffer frames whose amplitude and intensity compare with the track's the other way round),
`run.empty` (every track dropped: MATLAB's files without tracks, row 9 of [phase3.md](phase3.md)),
`run.late` (the first frames without detections: MATLAB's tracks keep every channel's row, where
the port took the channels from the first frame's detections and failed on tracks with gaps) and
`run.demo_cut` (the demo's frames 41 to 60 of its first experiment around a track whose buffers
the amplitudes and the intensities decide the other way round, deskewed as the demos deskew them
and cut to 52 x 51 voxels, with the 55 tracks of a MATLAB run of the multiple-experiments demo
that lie in them): every field but those of W11 and W12, in both files.

The keep rule's bound (C10, 2026-11-22): a free fit is kept when its position is less than 2
voxels from the window's centre in x, y and z, as MATLAB's: it refits one at exactly 2 (`< 2`,
MATLAB's bound as C09 found it on the multiple demo's one such fit, W13; the port's had been
`<= 2`). The oracle
case `run.keep_bound`: 12 designed windows (the first frames of end buffers and a gap frame,
along z, x and y, both channels) whose last valid voxel along one axis lies 2 voxels beyond the
centre (a block of the detection mask covers the voxels after it), with a bright spot beyond
that face, so that the free fit ends on the face, exactly 2 voxels from the centre: MATLAB
refits all 12 at the given position, and so does the port (with `<= 2` it keeps them, and the
case fails); every field but those of W11 and W12, in both files.

| Row | What | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| W11 | the fits at the gap frames and in the buffers: x, y, z, A, c, A_pstd, c_pstd, sigma_r, SE_sigma_r; what follows from them: `MotionAnalysis.MSD` and `MSDstd`, the rotated x and z, the end buffer's `Y` | demo: 80% (A) to 99.8% (x) of the 83,976 gap-frame values bit for bit, all but W13's within 1e-11 of max(\|v\|, 1); `run.designed` within 2.4e-13 | `rel` 1e-11 with `abs` 1e-10 (`TRACKS_TOL` of `tests/detection/test_processing.py`, on both track files) | the free fits' last bits (D3: `fit_cd`) and the refits' (D4: the mex's `'Ac'`); the tolerance is inside D3's |
| W12 | `pval_Ar` of the gap frames and buffers | demo: 2.8e-13 (absolute), 82,139 of 83,976 gap-frame values bit for bit; `run.designed` bit for bit | `abs` 1e-10 (the same `TRACKS_TOL`) | the Welch test's degrees of freedom in another order (91% bit for bit even from MATLAB's own fits), `tcdf` (S3), and W11 |
| W13 | the demos' fits that take another path | the by-step demo on MATLAB's inputs: 2 of 42,412 windows beyond W11 since C08's stop test (2026-11-21; 24 on the merged port before it; C04b's measurement 25: 11 gap frames, 6 start-buffer and 8 end-buffer frames, positions to 2.8e-4, A to 1.6e-2 relative, MSD and MSDstd of their tracks, 50 of 18,866 values, to 1.2e-4): a start-buffer frame at 0.018 voxels and A 1.6% relative, a gap frame at 1.6e-7. The multiple-experiments demo on MATLAB's inputs (proposed addition, 2026-11-21; C08's measurement with its stop test and, since 2026-11-22, the keep rule's bound; before them in brackets): 2, 12 and 5 of its experiments' 42,412, 71,360 and 58,154 windows (the second's 13 before the bound; 24, 186 and 36 before the stop test; C04c measured 25, 185 and 37 before C02d's box), in 2, 12 and 5 tracks; the first experiment's as the by-step demo's (its frames and inputs are the same), the third's within 8.6e-8, the second's but two within 7.3e-4 voxels and A within 5.0e-3 relative; in those two (an end-buffer frame and a gap frame of the second experiment's channel 2) one side keeps its free fit where the other refits at the given position, the free fits ending on either side of the keep rule's 2 voxels: positions 1.77 and 1.33 voxels and A 63% and 40% relative apart. A third, an end-buffer frame of channel 1 (2.27 voxels and A 615% relative apart), was the keep rule's bound: the free fits agree, exactly at 2 voxels, where MATLAB refits and the port kept its free fit (C09); it is MATLAB's since the port's bound is MATLAB's `< 2` (C10, 2026-11-22: the oracle case `run.keep_bound`, 12 designed windows whose free fits end exactly 2 voxels from the centre, refitted on both sides). (The third experiment's gap frame of C04c's measurement, 1.04 voxels, is MATLAB's since C02d's box; the 2.27 one came with it.) | none on the oracle cases; the demos: listed | the free fits (D3: `fit_cd`) end elsewhere than MATLAB's on the same windows and starts. Of the 4,320 fits of the 239 tracks above, recorded from MATLAB (every window and start bit for bit), 465 free fits ended beyond 1e-9 of MATLAB's (C04c), weak spots mostly (median \|A\| 3.7, three quarters below 10; 221 of them refitted at the given position, which hides the difference); 422 after C02d's box, and 56 since C08's stop test: in 290 MATLAB's end had been `fit_cd`'s iterate one iteration before `fit_cd` stopped, and the stop test's bound on the RSS's change is absolute (1e-6), not 1e-10 of the first RSS ("The fits' stop test (C08)" above). The 56 left take another path; the 64 beyond 1e-12, decoded fit by fit (C09, "The fits still on another path" above): their last bits part first at MATLAB's update of A (or at a damped step), and 56 part decisively at the box (MATLAB's position an ulp inside a face, the port's on it), 3 at a switch at a slightly negative A, 5 stay within 1e-6. Of the three that were on either side of the keep rule, two part at such a switch (A -3.3e-4 and -4.0e-5); the third (channel 1's end buffer at frame 30, not among the tracks recorded) was not on another path: MATLAB's free fit of the port's window and start is the port's, exactly 2 voxels from the window's centre, where MATLAB refits and the port kept it (MATLAB's bound `< 2`, the port's `<= 2` in `processing.fit_windows` until 2026-11-22 and `< 2` since, C10; the only such fit of the demo's 171,926). Not the inputs' last bits: with the window or the start changed by 1 to 3 ulp, 5 of the 465 ends move, 1 to MATLAB's (C04c) |
| W14 | `catIdx` | none: every category equal (the by-step demo's 3,636 tracks, the multiple-experiments demo's 3,636, 6,553 and 5,193 on MATLAB's inputs, 717 probe tracks, the six oracle cases `run.*`) | none | fixed 2026-11-21 (C04c; was 4 of the by-step demo's tracks, and 11 and 19 of the multiple demo's second and third experiments): the buffer test compares a buffer frame's intensity, `A + c`, with the track's largest `A + c` over all its columns (the gap frames' fits too), not its amplitude with the track's largest amplitude. This rule alone gives MATLAB's 15,382 categories of the multiple demo; the amplitudes, the largest `A` plus the largest `c`, the detected frames only or the frame of the largest `A` leave 2 to 11 of its 6,986 single visibility-1 tracks of three frames or more in the other category. Cases `run.buffers` (spots in the farthest buffer frames, which the significance test does not read, dimmer than the track alone but brighter on a raised background, and brighter alone but dimmer on a lowered one) and `run.demo_cut` (two of its real tracks decided the other way round by the amplitudes) |
| W15 | the Amira files (`Amira/Tracking/tracking_t_<k>.am`) | the text differs where W11's values are printed (16 digits) | not compared on the oracle case; `test_amira_on_matlab_tracks` checks the writer byte for byte on MATLAB's tracks | W11 |
| W16 | `ProcessedTracks.mat`, `RotatedTracks.mat` as files | the same variables, values, classes, sizes, field order and NaN bits | none (compared by content) | the HDF5 layout (object order, times) is the writer's; `RotatedTracks.mat` is `ProcessedTracks.mat` copied with the rotated values written into it, which MATLAB loads to the same values as the file it writes itself (`isequaln`) |
| W17 | options other than the demos': `Buffer` other than `[3 3]`, `FitMixtures`, `WindowSize`, `BufferAll`, `Preprocess` and `Postprocess` false, `Frames` other than all; `rotateTracks3D` of objective scans and with `'Crop'` true | not measured | `NotImplementedError` | not measured |
| W18 | a single track whose span is longer than its columns (a one-frame branch that merges into a segment at that segment's first frame is dropped, the span kept) | MATLAB's `XR_amiraWriteTracks` stops with `MATLAB:badsubscript` in its `parfor`, after `ProcessedTracks.mat` is saved and some frames' files are written, and `XR_runTrackProcessing3D` raises it; the port writes every frame's file (the track's points up to its columns) | none: the oracle cases leave the structure out (`run.compound`) | MATLAB indexes the amplitudes by the span |
| W19 | the Amira files' line ends on Windows | not measured: MATLAB writes each file's header and section marks with `fprintf` to a file opened `'w'`/`'a'` (`\n`) and the data rows with `dlmwrite`, whose line terminator on Windows its documentation does not give (no Windows MATLAB was used); the port writes `\n` everywhere, as MATLAB does on Linux | none (the files are compared on Linux only) | not measured |
| W20 | `rotate_tracks` of a reversed scan (the movie's `reversed`) | MATLAB's `rotateTracks3D` does not read `reversed` (its output the same for `reversed` true and false): it rotates every scan by the angle. A reversed scan deskewed in its direction (`XR_deskewRotateFrame` with `reverse`, as the port deskews it: W2) needs the mirrored angle; MATLAB's tracks of that scan miss the beads of its DSR volume (`XR_deskewRotateFrame` with `reverse` and `rotate`) by up to 20.2 px in x and z | a reversed scan rotated by the mirrored angle about the same centres (`rotation_matrix(mirrored=True)`: z mirrored in the DS volume and in the rotated box); `rotate_tracks(matlab_quirks=True)`, and so `run_detection_tracking(matlab_quirks=True)`, give MATLAB's | probe jobs 4635205 and 4635214 (oracle tree 45b2bb1): six blobs in a synthetic 24-plane scan (31.5 degrees, dz 0.5, 0.104 um pixels), their centroids in the DS and DSR volumes (the separate and the combined deskew and rotation), tracks of the DS centroids; the frame of the rotated tracks to the DSR volume's voxels fitted on the forward scan (residuals 0.34 px). Reversed: the mirrored angle within 0.49 px of the beads, MATLAB's rotation up to 20.2 px off. The DSR volume of a reversed scan is itself the z mirror of the forward scan's of the same sample (within 0.4 px; the raw stack in the other order), so a reversed scan's volume and tracks both come out upside down against a forward scan's. `deskew_rotate_file` gives the probe's DS and DSR volumes bit for bit. Case `rotate.reversed` (MATLAB's output, `matlab_quirks`); `test_reversed_scan_rotation` (the port's deskew and rotation) |

## Lifetime analysis (`detection.lifetime`, C05)

`runLifetimeAnalysis3D` and `plotIntensityCohorts` (the numbers; the figures are not ported),
measured as black boxes: the GPL functions they reach (`getLifetimeData`, `loadTracks`,
`scaleEDFs`, `fitGaussianModeToCDF`, `getVisitorIndex`, `adtest1`, the `binterp` mex) through
their inputs and outputs, and the MathWorks functions (`lsqnonlin`, `hist`, `interp1`,
`normcdf`, `erf`) from documentation and probes. Cases: `tests/detection/test_lifetime.py`
(synthetic conditions of one, two and three movies, the functions alone, and `lifetime.demo`:
MATLAB's `ProcessedTracks.mat` of the multiple-experiments demo, three movies, 1,160 tracks of
category Ia). Approved by the user on 2026-11-20.

Exact, bit for bit: the files (`lifetimeData.mat`, `Parameters.mat` but its threshold), the
lifetime histograms of every category and their counts, the visitors, the CCP and CS
classification and every histogram that follows from it (but in one condition, L4), the starts
per frame, the cohorts' statistics of one movie (up to `binterp`, L7), `hist`
(centers; a value on an edge counts in the lower bin), `erf`, `erfc` (fdlibm's, 409,402 values)
and `normcdf`, `interp1` (`'linear'`, `matlabcompat.interp`), `getVisitorIndex`, and on the demo
data every field of `lftRes` but `MaxIntensityThreshold`.

The fits are `lsqnonlin`'s trust-region-reflective algorithm, rebuilt from Coleman and Li's
method and the rules measured on 469 problems with every evaluation recorded. Where the port's
path through a fit is MATLAB's, the results agree to about 1e-9; MATLAB sometimes takes another
step than the exact two-dimensional subproblem (the rule was not found), and on the non-smooth
objective of the mode fit (`fitGaussianModeToCDF` fits the CDF below the mode only) another step
can end in another local minimum.

| Row | Output | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| L1 | `matlabcompat.optim.lsqnonlin` on smooth problems (exponential decay, a Gaussian CDF, Rosenbrock with an active bound) | 8.8e-10 relative; resnorm 1e-9 | `rel` 1e-6, `abs` 1e-12 | the rounding of MATLAB's products and solves (BLAS and LAPACK) |
| L2 | `scaleEDFs` (the movies' scale factors `a` and offsets `c`) | where the path is MATLAB's (the demo's 8 fits, 85 of 100 random problems, the cases of lognormal samples): `a` 8.3e-10, `c` 7.9e-7; another path (3 cases: synthetic maxima, a grid past the reference's values): `a` 7.3e-4, `c` 3.6e-3 absolute | `rel` 1e-6, `abs` 1e-9; another path `rel` 1e-3, `abs` 5e-3 | the subproblem's steps (above); `c` is near 0 |
| L3 | `fitGaussianModeToCDF` (the Gaussian mode: mu, sigma) | the path is MATLAB's in 26 of the demo's 61 fits and 100 of 300 random problems (within 1e-6); otherwise up to 2.1% (cases: 3.3e-4, 6.4e-5, 2.1%), and up to 74% on the demo's two-frame cohorts | `rel` 0.03 (cases) | a non-smooth objective: another step, another local minimum |
| L4 | `MaxIntensityThreshold` (`lftRes`, `Parameters.mat`) | demo 9.9e-10; synthetic conditions 5e-6 to 7.9e-5; a synthetic condition of two movies 8.4% | `rel` 1e-4; where the fit ends in another local minimum `rel` 0.1, the fields of the CCP and CS classification not compared | L3 and L2 (the threshold is mu + sigma `norminv(0.99)` of the mode of the scaled maxima): in the two-movie condition the scale factors took another path (7.7e-4) and the threshold's fit, on inputs differing by that much, another minimum, so 2 of about 130 tracks were classified otherwise; elsewhere the classification was the same in every run, and so were `nFramesThreshold` and the ten Anderson-Darling decisions before it |
| L5 | `adtest1`'s statistic and p-value | 2e-14 relative | `rel` 1e-13 | MATLAB's `log` (138 of 144 probed statistics bit for bit); the decision `H` exact |
| L6 | `binterp` (cubic B-spline, mirror boundaries) | 2.2e-14 relative (cases), 7.6e-13 of the value on the demo's 1,948 calls | `rel_max` 2e-13, `abs` 1e-13 | the recursive filter's order (Thevenaz, Blu and Unser's algorithm, not the mex's code) |
| L7 | the cohorts' `A`, `Amin`, `Aplus` | several movies: demo 2.3e-8 relative (5e-8 absolute), synthetic 1.9e-3; one movie (no scaling): 4.7e-14 | several movies `rel` 5e-3, `abs` 1e-9; one movie `rel_max` 2e-13, `abs` 1e-12 | the scaled amplitudes (L2's factors); with one movie `binterp` (L6) |

Errors kept from MATLAB: a lifetime cohort without tracks raises `stats:ecdf:NotEnoughData`
(`ecdf` of nothing), constant samples `optimlib:commonMsgs:UsrObjUndefAtX0` (the mode fit's
start), fewer than five values `adtest1`'s error, and `binterp` NaN samples and positions outside
[0, n + 1]. The last three are MATLAB errors without an identifier whose messages are
third-party text: their cases compare the identifier only (the empty one), and their manifest
entries withhold the messages (`Case(error_text="withheld")`; decided by the user on
2026-11-20). Not reproduced, as MATLAB's result depends on memory beyond the input: `binterp` of
one sample, and of fewer than four samples at positions outside [1, n].

Options of `runLifetimeAnalysis3D` (registered with MATLAB's interface; the port's checks of its
options fitted to MATLAB's parser on 20 probe cases, `tests/params/test_third_party_checks.py`):
`AnalysisPath`, `ProcessedTracks` and `LifetimeData` (case `lifetime.file_names`), `Overwrite`,
and `Display` (`lifetime.display_off`: MATLAB's numbers and printed values are the default's)
are ported. The others are not measured: at other values than their defaults they raise
`NotImplementedError`, and so does a cell mask at `<source><MaskPath>` (the cell areas and
initiation densities stay NaN).

**MATLAB's quirks kept**: the tracks are sorted by lifetime, longest first, before anything
else (`index` numbers that order); the visitors are not scaled; the second channel's CS
histograms (`lftHistSlaveCS`) are normalized by their own count, not by the CS tracks of all
combinations; `startsPerFrameAll` counts the starts of frames 6 to `movieLength` - 2;
`plotIntensityCohorts`' standard error divides by the square root of the number of movies, also
where a movie has no track in the cohort, and with one movie its curves are another statistic:
the median of the tracks between their 25th and 75th percentiles; the mode fit starts from half
the standard deviation.
