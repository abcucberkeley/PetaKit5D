# Differences from MATLAB: point detection

The port matches MATLAB bit for bit unless an oracle case declares a tolerance
(`petakit5d.testing.compare.tol`, see [testing.md](../testing.md)). This page lists the
proposed tolerances of `petakit5d.detection`, with what was measured and why the bits differ.
Each case's `reason` names its row as `detection_scope.md statistics Sn` (the scope document
has the same rows). All rows are **proposed** (2026-11-19), waiting for the user's decision
([index.md](index.md)).

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
| D3 | the detector's fits and results (positions, A, c, their standard deviations, RSS, sigma_r, p-values) | demo frame 1: positions of 1,043 of 1,043 plain fits and 188 of 221 restarted ones within 1e-12; over all fits positions 2.5e-10, A 2.7e-7, standard deviations 1.3e-7 relative; the oracle cases within 1e-9 | `rel` 1e-9, `abs` 1e-9 on the oracle cases; demo frames: the same detections, the largest differences listed | D2's last bits, the restarts' direction solves (the damped steps' eigenvalues) and the A update's sums, whose order MATLAB keeps in ways not all measured |

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
| D8 | `XR_runDetection3D`'s `frameInfo` with its defaults ('speedup', no background correction) | 2 demo frames, both channels: the same 661 and 649 detections; the master's values within 3.3e-9 relative (D7), the other channel's A within 1.6e-4, c 9.3e-7, positions 6e-15; every test (`hval_AD`, `isPSF`) the same | D7's, and `rel` 1e-3 for the other channel's A and c (proposed; no case yet) | D3 in the other channel's `'xyzAc'` fits from the master's positions, which move up to 4 voxels; one in 600 stops where the residuals barely change |

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
`[2 2]`, the port's window 3,128; refused), and 'lowSNR' without the background correction (the
wrapper's default): frame 1 has MATLAB's 1,006 detections against 1,005, frame 2 the same 969;
given MATLAB's points the bad-point filter keeps MATLAB's 1,006 and 969, and all but one of the
1,040 points into it agree within 1e-6 (one fit, moved 2.9 voxels, 2.9e-4 off: D3), which the
filter passes on.

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
names its row as `detection.md Wn`. Proposed (2026-11-19), waiting for the user's decision.

| Row | Function and range | Measured (largest) | Tolerance | Cause |
|---|---|---|---|---|
| W1 | `corr` (Pearson, two columns): the 165 correlations of `XZoffsetEstimation3D` and the probed columns | 17 ulp where the sum cancels (small integer columns); 16 ulp at 65,536 values; 2 ulp at 5.7 million; single: 1 ulp of single | `abs` 4e-15 (single `abs` 2e-7) | MATLAB takes `x'*y / norm(x) / norm(y)` from its BLAS, whose order of additions depends on the CPU and, from 4,096 values, on its thread count (MATLAB's own results differ between 1 and 8 threads); the port sums in MATLAB's order, platform-independent. The offsets the correlations decide were equal in every case |
