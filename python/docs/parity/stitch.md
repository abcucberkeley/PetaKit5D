# Differences from MATLAB: stitching

All other stitching outputs are bit-identical to MATLAB. The 3D registration (the NCC of
`kernels.xcorr`) is bit-identical with the exact FFT backend and within decision record 0006 on the
permissive one (same integer shifts, peaks within 5e-3; measured 4.2e-7). The entries below are
for the 2D registration, which uses MathWorks' `normxcorr2`.

## normxcorr2 (proposed)

`petakit5d.stitch.xcorr.normxcorr2` is a new implementation of the Image Processing Toolbox
function, from its documentation and the method it cites (J. P. Lewis, "Fast Normalized
Cross-Correlation", 1995). MATLAB's function is MathWorks' code (black box): its FFTs and running
sums are its own, and its results change in the last bits with its thread count (up to 7e-16 between
1 and 4 threads for a 2304 x 54 template in a 2304 x 354 image). So the port's coefficients are not
MATLAB's bits.

| Output | Tolerance | Measured |
|---|---|---|
| `normxcorr2` coefficients | absolute 1e-7 | 2.6e-13 at most where the data have contrast (16 oracle cases and the probes; against the exact NCC of integer data MATLAB's error is up to 2.1e-13, the port's up to 2.6e-13); 9.2e-8 on 16-bit data at 2^15 with a contrast of 12, where it is MATLAB's own error (the port's is 1.6e-15 against the exact NCC) |
| 2D registration peaks (`max_xcorr`) | absolute 1e-7 | 5.3e-8 (offset data), 1e-13 elsewhere (87 pairs) |
| 2D registration shifts | 0 | 0 (87 of 87 pairs identical) |

Known differences, not covered by the tolerance (no case of the registration's oracle set has them):

1. **MATLAB loses a peak of 1 on offset data.** MATLAB sets coefficients to 0 where its rounding
   takes them beyond 1 (seen: a coefficient whose exact value is 1, at a template's own place in
   16-bit data at 2^15 with a contrast of 12; MATLAB kept 1 + 1.2e-8). The port's coefficients are
   exact to about 1e-15, so it keeps the peak. A 2D registration of such data would find the
   template's place in the port and another place in MATLAB.
2. **Flat and nearly flat windows.** Windows whose values are all equal give 0 in the port (found
   exactly). MATLAB gives rounding noise in them (up to 4e-7 in all-zero windows next to data with
   a contrast of 1e3, 2e-8 in windows of a constant 500) and, in windows whose variance is below the
   rounding of its sums (a variation of 1e-4 on 500), noise up to 0.45 or 0. A 2D registration
   whose whole window of allowed shifts is flat gives the first allowed shift in the port and a
   place of MATLAB's noise in MATLAB.
