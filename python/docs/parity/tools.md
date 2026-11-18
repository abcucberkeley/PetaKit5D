# Differences from MATLAB: tools

The port matches MATLAB bit for bit unless an oracle case declares a tolerance
(`petakit5d.testing.compare.tol`, see [testing.md](../testing.md)). This page lists the
tolerances of `petakit5d.tools` beyond record 0006's FFT backends, with what was measured and
why the values differ. Each case's `reason` names its row.

## FFT power spectra (`tools.fft_spectrum`, approved 2026-11-17)

| Output | Backend | Tolerance | Status | Why |
|---|---|---|---|---|
| MIPs and central slices (gamma 0.5) | exact | bit for bit | measured: every such file of the 27 exact cases identical | MATLAB's driver and FFTW, MATLAB's abs of complex single (msun's hypotf), the max and the normalization in single |
| the 3D stack (`save3DStack`: log10 of the normalized spectrum) | exact | within 2 ulp (`LOG10_ULP`) | approved | MATLAB's single `log10` is not correctly rounded: on 851,973 probe values in (0, 1] it is 1 ulp off the correctly rounded value at 13.9% and 2 ulp at 0.17%. The port's is correctly rounded (the C library's double log10, with ties decided in 40 digits), so it is the same on every platform. Measured on the stack cases: 3,692 of 226,981 voxels differ (frame.stack), all within 2 ulp. B10 found the same for the figures' log10 |
| MIPs, central slices and the 3D stack | permissive | the normalized spectrum `\|F\| / max` (the images squared, 10^ of the stack) within 0006's single tolerance: largest difference <= 1e-5, relative L2 <= 1e-6 (`GAMMA_CLOSE`, `STACK_CLOSE`) | approved | the FFT backend. Gamma 0.5 enlarges differences near 0 (the square root of a difference of 1e-9 is 3e-5), so the images are compared as the spectrum they show. Measured against the exact backend: largest difference of `\|F\| / max` 4.8e-8, relative L2 at most 2.9e-7 (cases up to 257^3); of the gamma images 1.5e-4 (the cosine's near-zero values) |

The probe case `probe.log10_sqrt_single` pins both functions: `x .^ 0.5` of single data is the
correctly rounded square root (bit for bit), `log10` is MATLAB's within 2 ulp. The closest
formula found for MATLAB's single log10 is fdlibm's single-precision `log10f`
(`k*log10_2hi + (k*log10_2lo + ivln10*logf(m))`) with fdlibm's `logf`, which reproduces
99.83% of the probe's values; correctly rounded log10 reproduces 86.1%, so the port keeps the
correctly rounded one.

## PSF analysis (`tools.psf_analysis`)

| Output | Backend | Tolerance | Status | Why |
|---|---|---|---|---|
| PSF slices (xy, xz, yz) of `Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update` | both | bit for bit | measured: 8 cases, 2 dataset cases | no FFT: MATLAB's median, `imresize3` linear, `imgaussfilt3` and divisions |
| OTF slices (xy, xz, yz) and line cuts | exact | bit for bit | measured: 7 cases, 2 dataset cases | `fftn` on FFTW with MATLAB's driver, `abs` as msun's `hypot`, `interpn` linear (measured: along y, then x, then z, `(1 - t) a + t b`) |
| OTF slices and line cuts, the info files and `RW_info.mat` | permissive | record 0006's double tolerance: max <= 1e-13 of the largest value, relative L2 <= 1e-14 | record 0006 | the FFT backend. Measured: at most 7.6e-16 of the largest value, relative L2 at most 1.1e-15, in 5-12% of the values of the 1001 x 1001 slices and 40-90% of the line cuts |

## PSF detection (`tools.psf_analysis`, approved 2026-11-17)

`XR_psf_detection_and_cropping` finds the PSFs with MATLAB's `radon` at the skewed angle
`atand(tand(a) * xy / (dz * sind(a)))`, `imregionalmax`, `median` and `std`. The port's
Radon transform is written from MathWorks' documented method (four subpixels per pixel, each
split linearly between the two nearest radial bins), without its exact arithmetic:

| Step | Measured | Tolerance |
|---|---|---|
| `radon` (internal, not an output) | single-pixel probes: the radial coordinates differ from MATLAB's by up to about 8 ulps of the terms (18 of 72 probes bit for bit); sums of random images within 4e-14 of their values | none: the decisions are compared |
| `tand` and `atand` (internal) | within 1 ulp; 81% of 20,016 and 20,009 probed values bit for bit | none |
| the decisions: the peaks (`pcoords`), the crops and their MIPs, the info file | bit for bit in the 4 cases (two beads, seven beads with a close pair, one at the border, a saturated and a dim one, the defaults at another angle, an image in two parts) | exact |

A peak can come out differently only where a last-bit difference of the Radon transform decides:
a regional maximum against an equal neighbour, a peak at the local threshold, or a voxel at the
10-pixel distance from its line. None of the cases has one; real bead images are noisy enough
that such ties are rare. Not compared bit for bit; approved on 2026-11-17 without MATLAB's `radon` to the last bit.
