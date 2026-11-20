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
