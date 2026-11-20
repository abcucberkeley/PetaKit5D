# Differences from MATLAB: SIM reconstruction (`petakit5d.sim`)

MATLAB has no tests of `third_parties/sim3D`; the oracle cases are the port's own
(`tests/sim/test_sim_steps.py`, `tests/sim/test_sim_recon.py`): synthetic 5-phase lattice data
with a known pattern (`sim.synthetic`, also 3 phases and 3-beam), the 2020 lattice SIM data of
the demo data (5 and 3 phases) with the settings checked against the original 2020 workflow, and
GPU cases. Each case's `reason` names its row. The rows are **proposed** (C06, 2026-11-19); the
measured values are the largest of the cases of the row, the port's default FFT backend against
MATLAB R2024b on Linux x86-64 (GPU: an A100).

## Bit for bit

- `make_forward_separation_matrix` (fdlibm's sine and cosine, MATLAB's signed zeros: a value
  whose imaginary part is zero is real in MATLAB, so its exponential has a +0 imaginary part);
- `cylMask`, `tukwin` (MATLAB's colon), `resampleOTF` (linear: along x, then y, then z,
  `a w0 + b w1`), the apodization of `apodizeEllipse` on the lattice case;
- `imtranslate_function` for whole-voxel moves, `fminsearch` (every iterate, function value,
  iteration and evaluation count of MATLAB's traces through an output function, on five
  problems, one of them piecewise constant, where the points compared are equal);
- `deskewPhasesFrame` (deskew and rotation, single and uint16), `separatePhases`,
  `combinePhases`, and the files `simReconFrame` writes for skipped, existing, missing and
  blank data and for erosion before the reconstruction.

## Tolerances

| Row | Outputs | Tolerance | Measured | Why |
|---|---|---|---|---|
| R1 | the steps in double (`fourierShift3D` with `gpuPrecision` double, `perdecomp_3D`, `imtranslate_function` fractions, the pseudo-inverse, the 3-beam apodization window) | max <= 1e-13 of the largest value | 8.8e-16 | FFT rounding (the default backend; the exact add-on has FFTW's codelets but the same arithmetic around them), LAPACK's SVD (pinv), the order of the trilinear blend's operations |
| R1 | the steps in single (`fourierShift3D`, `perdecomp_3D` in single) | max <= 2e-6 of the largest value, relative L2 <= 1e-6 | 2.4e-7, 1.4e-7 | MATLAB's single-precision sine and cosine of the phase ramp are not the double functions rounded (1 ulp apart for 6.2 % of 200,001 probe angles); single FFT rounding |
| R2 | `Lsq_GaussianFit_3D` (and the bead centre of `sim_PSFtoOTF_gen`) | relative 1e-6 per parameter | 1.8e-9 | MATLAB sums the residuals (1.65 million for the 2020 PSF) on all its threads, the port in MATLAB's order at one thread, and exp is the C library's: the objective's last bits differ, and the simplex's path with them (fminsearch itself is bit for bit) |
| R3 | `MaskedTranslationRegistration2D_fit` and `normxcorr2_masked` in double: translation, peak, correlation image | max <= 1e-9 of the largest value | translation 3.2e-14, correlation 1.2e-14; `normxcorr2_masked` 6.4e-11 (its smallest overlaps, which nothing thresholds there) | FFT rounding in the correlation sums, divided by small overlap sums |
| R3 | the same in single (the default `gpuPrecision` of the registration) | max <= 5e-5 of the largest value, or 1e-5 | translation 4.6e-6 (1.5e-5 px), peak 7.1e-7, correlation 4.8e-6 | single precision in the sums |
| R4 | `sim_PSFtoOTF_gen` | max <= 1e-5 of the largest value, relative L2 <= 1e-5 | double 7.6e-11, 8.9e-11; single 1.0e-7, 1.6e-7; GPU 4.1e-7, 3.1e-7 | R1 and R2 together |
| R5 | `simRecon` and the files of `simReconFrame` on the CPU | max <= 2e-4 of the largest value, relative L2 <= 1e-4 | `simRecon` 2.2e-6, 3.0e-6; files 2.2e-6 (small chunks 2.2e-5, 1.1e-5; odd chunks 1.4e-5, 9.3e-6, with NaN in MATLAB's places, 21 % of the voxels; uint16 1.5e-5: +-1 in 0.2 % of the voxels); the 2020 data 1.5e-6, 8.4e-7 | R1 to R4 through the reconstruction: the refined shift vectors (R3) and the weights fitted on the overlaps move and scale the bands (more in small chunks, each refined on its own); rounding to uint16 at .5 |
| R6 | the same on a GPU (CuPy against MATLAB's gpuArray) | max <= 1e-4 of the largest value, relative L2 <= 5e-5 | 2.4e-6, 2.2e-6 (A100; the 2020 data 1.2e-6, 6.6e-7) | single precision in both, cuFFT plans (11.0.8, decision 0003) and element-wise kernels of CuPy against MATLAB's |

## Notes

- `simRecon` with `wVer` 2 fails in MATLAB ("Reference to a cleared variable denom.",
  `MATLAB:refClearedVar`); the port raises the same error before computing.
- The dataset wrappers and `simReconAutomaticProcessing` fail in today's MATLAB tree (they pass
  `cpuOnlyNodes`, which the job framework no longer takes); the port runs them as written. Their
  oracle cases record MATLAB's error only; the port's wrappers are tested on their own, and the
  frame functions they run are compared with MATLAB's.
- sim3D's reconstruction is not the original 2020 workflow's: by default each order is
  Wiener-filtered by its own OTF (`wVer` 1; the original: by the sum of the shifted OTFs, `wVer`
  0), the OTF is interpolated linearly (the original: cubic) and the masks are moved by a rounded
  copy. `simRecon` with `wVer` 0 and `gpuPrecision` double reproduces the original's result on
  the 2020 data to 2.6 % (relative L2). The port follows sim3D.
