# Approved differences: matlabcompat

The port matches MATLAB bit for bit unless an oracle case declares a tolerance
(`petakit5d.testing.compare.tol`, see [testing.md](../testing.md)). This page lists the approved
tolerances of `petakit5d.matlabcompat`, with what was measured and why the bits differ. Each
case's `reason` names its row.

## Filters (`matlabcompat.filters`)

The filters are reproduced from MATLAB's behaviour only, by probing the built-ins on chosen
inputs. 306 of the 361 oracle cases are bit-identical, including the unmixing (`imgaussfilt3` in
single, sigma `[1 1 1]` and anisotropic) and puncta removal (uint16, sigma 2.5) call sites. The
others stop at the best black-box result: MATLAB's vector code reorders a few taps in ways the
probes don't pin down, its long-kernel `conv2` isn't a direct sum, and its frequency-domain bits
match no FFT size tried.

| Output | Measured | Tolerance |
|---|---|---|
| `imfilter` column pass, rows after the last vector (double 9 and 13 taps, 17 taps and more; `'full'`; constant padding; the PSF analysis site in double, sigma 2) | at most 3 ulp (single 4) in 0.3–12% of values of 64-row columns | `ulp` 5, `max_frac` per case |
| `imgaussfilt` with small 2-D kernels (fewer than 49 values, sigma below about 1.25) | at most 8 ulp in 50–77% of values | `ulp` 12; integers ±1 in at most 0.6% |
| `imfilter` with general 2-D and 3-D kernels (point detection's form) | at most 10 ulp | `ulp` 16 |
| `conv2` with long kernels (`imgaussfilt`'s separable path at 21 taps or more: the phase flat field's spatial path, the puncta MIP at sigma 50) | relative 1.4e-6; uint16 ±1 in 2 of 102,000 values | `rel_max` 2e-6; integers ±1 |
| The frequency domain (`imgaussfilt`, `imgaussfilt3`) and `edgetaper` | rounding level | `rel_max` 1e-6; integers ±1 |
| Sparse arrays (less than a tenth nonzero) with replicate borders | at most 3 ulp in 5.5% of values | `ulp` 4 |
| `imfilter` of uint16 along z | ±1 in at most 0.05% of values (cause not found) | `abs` 1, `max_frac` 1e-3 |
