# Differences from MATLAB: flatfield

The port matches MATLAB bit for bit unless an oracle case declares a tolerance
(`petakit5d.testing.compare.tol`, see [testing.md](../testing.md)). This page lists the
tolerances of `petakit5d.flatfield`, with what was measured and why the values differ. Each
case's `reason` names its row.

The light-sheet flat-field correction (`flatfield.correction`, D01) and the camera background
estimation (`flatfield.estimation`, B07: `XR_flatfield_background_estimation` and
`XR_flatfield_background_wrapper`) need no tolerance: all their oracle cases are bit-identical.

## The phase-image flat field (`flatfield.basic`, B07; approved 2026-11-18)

`XR_phase_image_flat_field_correction` estimates the flat field of a time point's phase tiles
with BaSiC (`third_parties/BaSiC-master`, CC BY-NC-ND 4.0, so not ported and not read). The
port uses BaSiCPy (MIT), another implementation of the same method (Peng et al., Nat. Commun.
8, 14836, 2017), with the options MATLAB passes (no dark field) and BaSiCPy's defaults
otherwise (`fitting_mode='approximate'`, the automatic smoothness, working size 128). BaSiCPy's
tighter tolerances (`optimization_tol` 1e-6, `reweighting_tol` 1e-3) and its `ladmap` mode
were tried: none brings it closer to MATLAB's result. The two implementations differ in their
solvers and resampling, so the flat fields agree to about a percent, not bit for bit.

The pre-filter before BaSiC (`imgaussfilt` sigma 100) is `matlabcompat`'s frequency-domain row
(rounding level), and so is `minval`, the smallest pre-filtered value (`abs` 0.01).

Measured with BaSiCPy 2.0.0 (PyTorch 2.15 on the CPU, one thread) against MATLAB R2024b on
synthetic tiles with vignetting (relative differences of the flat field; the corrected uint16
tiles differ by the same relative amounts plus their rounding):

| Case | Tiles | Flat field: max, 99th percentile, median | rel. L2 | Tolerance |
|---|---|---|---|---|
| `phase.tiles_128x8`: 128 x 128, strong vignetting (0.6 to 1) | 8 | 10.7%, 6.8%, 1.7% | 2.5% | `rel` 0.12, `rel_l2` 0.03 |
| `phase.tiles_96x160x12`, offset 6000 | 12 | 2.4%, 1.5%, 0.36% | 0.57% | `rel` 0.06, `rel_l2` 0.01 |
| `phase.tiles_200x256x20`, offset 6000 | 20 | 1.6%, 1.0%, 0.26% | 0.39% | `rel` 0.06, `rel_l2` 0.01 |
| the phase demo (`demo_phase_and_2photon_stitching.m`, not an oracle case), both time points: 2304 x 2304, offset 6000 | 20 | 5.6%, 2.8%, 0.33% | 0.64% | (the same) |

The uint16 tiles get `abs` 1 on top (rounding); on the demo they differ by 5.3% at most. BaSiCPy's own result changes in the last bits
with PyTorch's thread count (6e-7 relative between 1 and 8 threads), so the port fits with
one thread.
