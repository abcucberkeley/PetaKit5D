# Differences from MATLAB: decon

The port matches MATLAB bit for bit unless an oracle case declares a tolerance
(`petakit5d.testing.compare.tol`, see [testing.md](../testing.md)). This page lists the
tolerances of `petakit5d.decon` beyond record 0006's FFT backends, with what was measured and
why the values differ. Each case's `reason` names its row.

## deconvlucy (the 'original' method, `decon.other_methods`, B11; approved 2026-11-18)

| Output | Backend | Tolerance | Status | Why |
|---|---|---|---|---|
| single, uint8, uint16, int16 results; error ids; OTF realness | exact | bit for bit | measured: 53 of 53 such outputs identical | the iterations are MATLAB's double arithmetic; the last bits of lambda vanish in the rounding to the class |
| double results (double data; the cell form's double estimates and changes) | exact | max <= 1e-13 of the largest value, relative L2 <= 1e-13 (`DOUBLE_TOL`) | approved 2026-11-18 (parity row N2) | MATLAB's acceleration factor comes from BLAS (MKL) dot products, whose order follows the CPU's SIMD and the thread count; the port sums in a fixed order. Measured: <= 18 ulp, max 1.1e-15 of the largest value on the estimates, 1.3e-14 on the changes (relative L2 <= 1.1e-14) |
| all results | permissive | 0006 (relative L2 <= 1e-4; integers at most 1 count in at most 0.1%) | record 0006 | the FFT backend |

Not compared on the permissive backend: a PSF with negative values (relative L2 3e-4: negative
PSF values amplify FFT round-off; exact backend bit-identical).

## OMW with many iterations on the default FFT backend (P3; approved 2026-11-19)

The widefield decon of demo_widefield_and_confocal_deconvolution (`WF_raw.tif`, 501 x 571 x 501;
OMW, wienerAlpha 0.02, OTFCumThresh 0.65, hannWinBounds (0.8, 1.2), background 100, 20
iterations) is the one decon of the demos and the real-data workflow whose default-backend result
is outside record 0006 ([demos.md](demos.md), row 8). Measured iteration by iteration against
MATLAB's CPU result (task X14), as uint16 voxels off by 1 / off by more / the largest difference.
"Exact arithmetic" is the same OMW in double precision on the same single-precision PSF and back
projector (pocketfft's and SciPy's double transforms agree on it within 3.5e-13).

| Iterations | default backend | its iterations alone (MATLAB's back projector) | MATLAB's GPU result | exact arithmetic | MATLAB's bits, one PSF value 1 ulp up |
|---|---|---|---|---|---|
| 2 | 6,801 / 0 / 1 | 6,435 / 0 / 1 | 12,327 / 0 / 1 | 4,869 / 0 / 1 | 3,387 / 0 / 1 |
| 5 | 12,652 / 0 / 1 | 11,717 / 0 / 1 | 16,165 / 0 / 1 | 8,765 / 0 / 1 | 5,594 / 0 / 1 |
| 10 | 18,185 / 82 / 68 | 17,650 / 64 / 48 | 24,973 / 69 / 52 | 12,952 / 64 / 55 | 7,215 / 46 / 35 |
| 20 | 505,369 / 4,586 / 1,782 | 337,279 / 3,130 / 866 | 332,901 / 3,061 / 817 | 345,464 / 3,182 / 864 | 245,604 / 2,369 / 807 |

- **MATLAB's CPU result is the one far from exact arithmetic.** At 20 iterations it is 1.1e-3
  (relative L2) from it, with 3,182 voxels off by more than 1; MATLAB's GPU result is 1.3e-4 from
  it (123 voxels) and the default backend's iterations 9.6e-5 (75 voxels, the largest by 68). At
  the brightest differing voxel, (z, y, x) = (297, 78, 123), MATLAB's CPU result is 48,379, exact
  arithmetic 47,515, MATLAB's GPU result 47,562 and the default backend's 47,599. The difference
  to MATLAB's CPU result is that result's own rounding, amplified. The back projector, whose last
  bits differ (P4), adds about what moving one PSF value by one ulp does.
- **Users' measures.** At 20 iterations the default backend's result has a PSNR of 108 dB against
  MATLAB's (MATLAB's GPU result 112 dB), a mean SSIM of 1 - 1.3e-8 with no voxel's local SSIM
  under 0.995, and 3.3% of its z MIP's pixels differ (MATLAB's GPU result: 2.2%). The voxels off
  by more than 1 lie in 157 spots that hold 0.9% of the intensity; each spot's integrated
  intensity agrees within 3.2% (median 0.3%) and its centroid within 0.13 voxels. Changing the
  data by -1, 0 or +1 count per voxel (their noise is about 7 counts) changes the result 35 times
  as much (relative L2 5.8e-2, PSNR 77 dB).
- **When.** Up to 7 iterations at most 6 voxels are off by more than 1, by 2; from 8 on the
  differences grow, fastest from 14 to 17 (at 16 the largest is 23,789 counts), and from 17 on
  4,600 to 10,700 voxels are off by more. On the same data, hannWinBounds (0.8, 1.0) with
  wienerAlpha 0.02 gives 460 voxels off by more at 20 iterations (the largest by 1,541), and
  (0.8, 1.2) with 0.05 gives 2,049 (12,858); (0.8, 1.0) with 0.005 diverges on both backends
  (overflow from 8 iterations). Light-sheet data stay within 0006 with 20 OMW iterations
  (demo_useful_tools' debug decon: off by 1 in 0.0031% of voxels), and RL on the widefield data
  stays off by 1 at 30 iterations (in 0.134% of the voxels, above P2's 0.1% from 26 iterations
  on). That RL case is listed with P3 too (approved 2026-11-19); P2 is unchanged.
- **Not reduced by** the FFTs and OTFs in double precision (they give exact arithmetic's result:
  3,182 voxels off by more at 20 iterations, at 2.1 to 3.8 times the time per iteration) or
  pocketfft behind MATLAB's FFT driver (4,073, the largest by 1,429; 2.2 times): an accurate
  computation lands near exact arithmetic, and so as far from MATLAB's CPU result as that result
  is from exact arithmetic. The exact add-on gives MATLAB's bits, at 3.0 times the default
  backend's time per iteration (24 threads).
