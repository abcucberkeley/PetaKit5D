# Differences from MATLAB: decon

The port matches MATLAB bit for bit unless an oracle case declares a tolerance
(`petakit5d.testing.compare.tol`, see [testing.md](../testing.md)). This page lists the
tolerances of `petakit5d.decon` beyond record 0006's FFT backends, with what was measured and
why the values differ. Each case's `reason` names its row.

## deconvlucy (the 'original' method, `decon.other_methods`, B11; proposed)

| Output | Backend | Tolerance | Status | Why |
|---|---|---|---|---|
| single, uint8, uint16, int16 results; error ids; OTF realness | exact | bit for bit | measured: 53 of 53 such outputs identical | the iterations are MATLAB's double arithmetic; the last bits of lambda vanish in the rounding to the class |
| double results (double data; the cell form's double estimates and changes) | exact | max <= 1e-13 of the largest value, relative L2 <= 1e-13 (`DOUBLE_TOL`) | proposed | MATLAB's acceleration factor comes from BLAS (MKL) dot products, whose order follows the CPU's SIMD and the thread count; the port sums in a fixed order. Measured: <= 18 ulp, max 1.1e-15 of the largest value on the estimates, 1.3e-14 on the changes (relative L2 <= 1.1e-14) |
| all results | permissive | 0006 (relative L2 <= 1e-4; integers at most 1 count in at most 0.1%) | record 0006 | the FFT backend |

Not compared on the permissive backend: a PSF with negative values (relative L2 3e-4: negative
PSF values amplify FFT round-off; exact backend bit-identical).
