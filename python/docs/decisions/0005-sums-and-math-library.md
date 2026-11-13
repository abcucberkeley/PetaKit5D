# 0005: MATLAB's sums, statistics and elementary functions, exactly

Status: proposed (2026-11-12). Measured on Linux x86-64 (Xeon Gold 6146, glibc 2.28) against
MATLAB R2024b Update 5, with probes at 1, 2, 3, 4, 8, 12 and 24 computational threads and 653
oracle cases (`tests/matlabcompat/test_reduce.py`, `test_stats.py`, `test_elementary.py`).
Implemented in `petakit5d.matlabcompat.reduce`, `.stats` and `.elementary`.

## Context

Outputs of the pipeline rest on MATLAB's `sum`, `mean`, `std`, `var`, `median`, `prctile`,
`mode`, `cumsum`, `max`/`min` with their index, `sort`, `unique`, on its window functions and on
its elementary functions. Their order of operations, tie and NaN rules and math library are not
documented. NumPy differs in all of them: `np.sum` adds pairwise and in SIMD lanes that depend on
the CPU's dispatch, `np.median` and `np.percentile` interpolate otherwise, `np.argsort` is not
stable by default, and `np.sin`, `np.exp`, `np.power` are the C library's or NumPy's SIMD code.
The port's results must not depend on the C library, the CPU's SIMD dispatch, the thread count
or the Python version (decision of 2026-11-12).

## Decision

### Order of additions (single and double alike)

Found by testing order models against MATLAB's sums of inputs that are integers times powers
of two (each order gives its own result); the same structure holds for `sum`, `mean`, `std`,
`var` and outtype 'double', for rows and columns:

| Elements reduced | MATLAB's order at one thread | Threads |
|---|---|---|
| contiguous in MATLAB's memory, n <= 43 (`x(:)`, `'all'`, dimension 1, dimensions that only singletons precede) | index order | independent |
| the same, 44 <= n <= 88,999 | four lanes: element t in lane mod(t, 4), each lane in index order, then ((l1 + l2) + l3) + l4 | independent |
| the same, n >= 89,000 | blocks of 4,096 elements in four lanes, the block sums in index order | **depends on the thread count** (blocks split among threads) |
| a dimension that a non-singleton one precedes (`sum(x, 2)`, `sum(x, 3)`) | out(i) = x(i, 1) + x(i, 2) + ... in index order; in chunks of 1,024 (chunks added in order) from 2,048 along the dimension, and from 1,025 in arrays of more than 20,000 elements whose preceding dimensions hold 14 or more elements or that have more than one outer slice | independent |

- Adjacent reduced dimensions count as one (`[1 2]`: the vector rule on `x(:, :, k)`; `[2 3]`:
  the slice rule over `nx * nz` slices); runs of them are reduced one after the other in
  increasing order (`[1 3]`: the vector rule along y, then the slice rule along z). Singletons
  between runs do not merge them.
- Single data adds in single, double, integer and logical data in double (integers saturate at
  every addition with 'native'). `cumsum` is sequential in the class at any size.
- `mean` is the sum divided by the count in double, rounded to the class: `single(double(s)/n)`
  (single division differs above 2^24 elements). `var` is `sum(abs(x - mean(x)).^2)/(n - 1)` (or
  `/n`), both sums in the rules above; `abs` clears the sign of a NaN result. `std` is its
  square root. Integer data raise `MATLAB:var:integerClass`.
- Signed zeros: a sum of -0s is -0, an empty sum +0; `omitnan` skips NaN, and a NaN that starts a
  sum (a lane, a block, a chunk) counts as +0.
- MATLAB's multithreaded sums are deterministic for a given thread count (20 of 20 repeats equal
  at 24 threads) but the block split is not a simple static partition (no contiguous split model
  fits); the port reproduces only the one-thread order, at any thread count.

The port adds one IEEE operation at a time per lane, slice or block (NumPy element-wise `add`
and cumulative sums, never NumPy reductions), in MATLAB's linear order on the `(z, y, x)` arrays;
whole double sums use K08's `kernels.ops.sum_double`.

### MATLAB disagreeing with itself across thread counts

MATLAB's result at T threads against its result at one thread, n >= 89,000 (probe inputs with
strong cancellation; ordinary positive data shift less):

| Class | Function | 2 to 24 threads: cases differing | Largest difference |
|---|---|---|---|
| single | `sum`, `mean` | 23 to 31 of 25 to 31 | 1,929 ulp, 1.2e-4 relative |
| single | `std`, `var` (24 threads) | 3 of 6 | 3 ulp, 2.6e-7 relative |
| double | `sum`, `mean`, outtype 'double' | 20 to 30 of 25 to 31 | 256 ulp, 3.3e-14 relative |
| double | `std`, `var` (24 threads) | 2 to 4 of 6 | 2 ulp, 3.2e-16 relative |

On a 2 GiB single volume of positive values (`mod(7y + 13x + 31z, 1000) + 0.25`), MATLAB's
`sum(x(:))` is 268,314,083,000 at one thread and 268,303,352,000 at 4 threads (4.0e-5 relative),
its `mean` 499.773926 against 499.753937 at 4 and 499.751373 at 24 threads, its `std` 288.687408
against 288.672028 at 4 threads (5.3e-5 relative); each value is the same in 5 of 5 sessions at
a given thread count. The port gives the one-thread values at 1, 4 and 24 threads.

### Phase 1 sites whose MATLAB result depends on the thread count

At oracle commit d67e1b7. Every other reduction on Phase 1 paths is thread-independent (slice
rules, sums below 89,000 elements, selections, integer and logical counts).

| Site | Reduction | Effect | Proposal |
|---|---|---|---|
| `normxcorr3_fast.m:88-89` | `mean(T(:), 'double')`, `mean(A(:), 'double')` of single template and image | rounded to single afterwards (`cast(..., 'like', T)`): a 1e-14 change flips the single only next to a rounding boundary | a fixed order in MATLAB: `decon_sum_double(x)/numel(x)` (bug fix, as the parked decon change) |
| `normxcorr3_fast.m:136-137` | `sum(T0(:), 'double')`, `std(double(T0(:)))` | NCC values in the last bits (about 1e-14 relative); integer shifts only at exact ties | fixed order (`decon_sum_double`) in MATLAB |
| `normxcorr3_fast.m:146` | `var(A(:), 1)` of the single image | the flatness threshold `eps(single) * var(A(:), 1)`: windows at the threshold can flip (about 1e-5 relative in single) | fixed order in MATLAB, or the variance in double through `decon_sum_double` |
| `crop_subregion_by_intensity.m:29` | `std(single(region_2(:)))` (raw data, usually uint16) | the threshold `mu + 3*std`: voxels at it can flip and move the crop box (the mean of uint16 data, line 28, is exact in any order) | fixed order in MATLAB; else a tolerance on the crop box is not workable (it is an integer box), so the fix is preferred |
| `decon_lucy_function.m:213`, `decon_lucy_omw_function.m:189, 290` | `mean((double(J) - double(I)).^2, 'all')` | `err_mat` only, with `debug` | `decon_sum_double` (cosmetic) |
| decon PSF/OTF normalization, lambda | fixed already by `decon_sum_double` (the parked change) | none | done |

The port follows MATLAB's one-thread order at these sites, which is what MATLAB gives after the
proposed fixes.

### Ties, NaN and signed zeros

- `max`/`min`: NaN ignored unless all are NaN (then the first NaN); equal values keep the first
  (`max([0 -0])` is +0, `max([-0 0])` -0); over all elements, the first in MATLAB's linear order.
  Two arguments: a tie gives the first argument, NaN gives the other one; an integer and a
  double scalar compare in the integer class.
- `sort`: stable for equal values (-0 and +0 included), NaN last ascending and first
  descending. `unique`: sorted, each NaN a value of its own, -0 and +0 one value (the first in
  sorted order); 'first', 'last' and 'stable' as MATLAB.
- `median`: the middle element of the stable sort; for an even count `a + (b - a)/2`, or
  `(a + b)/2` when the signs differ or a value is infinite (`median([Inf Inf])` is Inf,
  `median([-realmax realmax])` 0, `median([-0 -0])` +0); integers `(a + b)/2` in double,
  rounded half away from zero. NaN gives NaN unless 'omitnan'.
- `prctile` ('exact', its default): NaN left out; `q = p/100*n`, `k = floor(q + 0.5)`,
  `r = q - (k - 0.5)`; `x(1)` below, `x(n)` above, `x(k)` for `r = 0` or equal neighbours, else
  `(1 - r)*x(k) + r*x(k+1)` (single: each product rounded to single, then added in single;
  integers in double, rounded). This operation order matched all 1,988 probe values; other
  algebraically equal forms missed 89 to 800 of them.
- `mode`: the smallest of the most frequent values, each NaN counted alone; with one or two
  outputs the last of the equal values in the stable sort, with three outputs the first.

### Elementary functions

MATLAB's functions of double data against the C library (glibc 2.28) and NumPy 2.5 on 318,207
dense and random values each (up to 116,010 for `exp`; 30,000 per exponent for powers):

| Function | MATLAB's | glibc 2.28 differs in | NumPy differs in | Port |
|---|---|---|---|---|
| `sin`, `cos` | fdlibm | 2.1 % (1 ulp) | 2.1 % | fdlibm in NumPy (E06's); exact up to 2^20 pi/2 |
| `tan`, `cot` (`1/tan`) | fdlibm | 2.4 % | 2.5 % | fdlibm `__kernel_tan` in NumPy; exact up to 2^20 pi/2 |
| `atan`, `atan2` | fdlibm | 1.1 %, 15 % | 1.1 %, 15 % | fdlibm (E06's) |
| `hypot` | msun | 0 | 0 | msun (E06's; glibc 2.35 and later differ) |
| `exp` | fdlibm, except `exp(1)` = e rounded | 4.2 % | 5.8 % | fdlibm with `exp(1)` |
| `sqrt` | IEEE | 0 | 0 | IEEE |
| `x.^p`, p = 2, 0.5, -1, 0, 1 | `x*x`, `sqrt`, `1/x` | 0 | 0 | the same |
| `x.^p`, other p (integer p too) | fdlibm `pow` | 9.5 to 10 % | 10.5 to 11 % | fdlibm `e_pow.c`, exact on 125 exponents |
| `sind`, `cosd` | the rule below (S03's, measured) | n/a | n/a | exact on every angle measured: a dense 0.1-degree grid of 7,201, 59,451 dense and special values to 1e16, S03's 1,979 angles, 27,632 random and boundary angles (to 2^27, multiples of 7.5 and their neighbours, +-45 + k ulp, -0, subnormals, 1e300, Inf, NaN) |

- `sind` and `cosd` (S03 found the rule; N05 completed it for `cosd`, `-0` and the rounding at
  45): exact values at multiples of 30 degrees (`sqrt(3)/2` rounded at 60); otherwise
  `a = fmod(|x|, 360)` (exact), `n = round(a/90)` (half away from zero, exact: `floor(a/90 +
  0.5)` rounds up just below 45 and leaves [-45, 45]), `d = a - 90n`, at `d = +-45` the cosine
  kernel (an odd quadrant for `sind`, an even one for `cosd`), and fdlibm's `__kernel_sin` or
  `__kernel_cos` of `d*(pi/180)` with the tail `d*lo(pi/180)` (pi/180 to 60 digits; the
  rounding error of the product is not part of the tail). `sind` is odd (`sind(-0)` is -0),
  `cosd` even (+0 at odd multiples of 90). Inf and NaN give MATLAB's NaN in `sind` and the
  positive NaN in `cosd`. The earlier `sin(z*(pi/180))` was 1 ulp off in 8 % / 4 % of a dense
  grid; S01's image-list code (`math.sin` of the old reduction) in 8 % of angles.
- One public module: `matlabcompat.elementary` holds all of these, E06's fdlibm and msun code
  moved there from `matlabcompat.neighbors`, with fdlibm's kernels (`kernel_sin`,
  `kernel_cos`, `kernel_tan`, `rem_pio2`) and its sine and cosine as they are (`fdlibm_sin`,
  `fdlibm_cos`) public; `neighbors`, the stitching geometry and the image-list generation take
  their trigonometry from it.
- Repeated multiplication does not give MATLAB's `x.^p` for integer `p >= 3` (up to 10 ulp off
  at `p = 15`); PetaKit5D's own `fastPower` does binary powering in the class, which the port
  reproduces as written.
- MATLAB returns the positive NaN for `cos`, `tan` and `cot` of Inf and NaN (fdlibm's `x - x`
  gives the negative default NaN, which MATLAB keeps for `sin`).
- Single data: MATLAB's single `sin`, `cos`, `tan`, `exp`, `log`, `atan2`, `hypot`, `.^` and
  `cosd` (2 of 44 single angles) are neither `single(f(double(x)))` (2 to 20 % of values differ, by 1 to 4 ulp) nor glibc's
  float functions (1 to 20 %). No Phase 1 code applies them to single data; the port computes them in double
  and rounds, a listed difference.
- Arguments above 2^20 pi/2: fdlibm's Payne-Hanek reduction is not written out; the port uses
  the C library's sine and cosine there (3 % of such values differ by 1 ulp). PetaKit5D's angles
  are below 2*pi.

### Windows and shifts

- `tukeywin(n, r)`, 0 < r < 1: with `t = (0:n-1)/(n-1)`, `(1 + cos(pi/per*(t - per)))/2` on the
  first `floor(per*(n-1)) + 1` points and `(1 + cos(pi/per*(t - 1 + per)))/2` on as many last
  points (`per = r/2`); exact on 737 of 737 probe windows. r >= 1 (the Hann window): 777 of 877
  windows of all r match with the mirrored half `0.5 - 0.5*cos(2*pi*t)`; for some n above 11 MATLAB's
  Hann differs in the last bit (its own cosine). Phase 1 calls `tukeywin(nw, 0.5)` only.
- `padarray` (constant, replicate, symmetric, circular; pre, post, both; pads larger than the
  array) and `circshift` move elements only and are exact.

## Proposed exceptions (for the user)

(The `sind`/`cosd` exception of the first version, 1 ulp in 8 % / 4 % of a dense grid, is gone:
both are exact with S03's rule.)

1. `median` and `mode` with both -0 and +0 among the equal values at the median or mode in
   arrays of more than a few elements: MATLAB's sign of zero follows its selection code (3 of 269
   probe medians, 504 of 1,616 column modes of rounded data with many zeros); values always
   agree. Same kind as N03's approved `medfilt` exception.
2. MATLAB's own sums of 89,000 elements or more at several threads (table above): the port gives
   the one-thread result; references of such sums are made at one thread.

## Consequences

- Port code takes sums, means, standard deviations, medians, percentiles and modes that feed
  outputs from `matlabcompat.reduce` and `.stats`, elementary functions from `.elementary`.
- Speed (2 GiB single volume, 24-core node, both processes under `numactl --interleave=all`, 5
  interleaved rounds): `median` takes 0.79 and `prctile` 0.18 to 0.56 of MATLAB's time at 1, 4
  and 24 threads. The NumPy sums of the first version took 15 to 22 times MATLAB's time for
  `sum` and `mean` and 8 to 30 times for `std`; K08's C++ kernels (single and double
  accumulators in MATLAB's order, the fused `abs(x - mean)^2` of `std` and `var`, the slice
  sums), wired into `reduce` and `stats` and bit-identical to them, give these ratios of Python
  to MATLAB at 1 / 4 / 24 threads (output frees not timed):

  | Reduction | `(z, y, x)` arrays | MATLAB's layout |
  |---|---|---|
  | `sum(x(:))` | 1.99 / 1.88 / 1.92 | 1.17 / 0.91 / 0.97 |
  | `mean(x(:))` | 2.06 / 1.86 / 1.86 | |
  | `std(x(:))` | 0.57 / 0.40 / 0.17 | 0.33 / 0.21 / 0.09 |
  | `sum(x, 2)` | 0.59 / 0.57 / 1.20 | |
  | `sum(x, 3)` | 0.56 / 0.51 / 1.18 | |
  | `mean(x, 3)` | 0.59 / 0.53 / 1.25 | |
  | `std` along dimensions 2 and 3 | 0.14 to 0.22 | |

- Open: the whole-array `sum` and `mean` of `(z, y, x)` arrays (1.9 to 2.1 times MATLAB's time)
  are the accepted cost of reading in MATLAB's linear order, a transposing read; hot paths keep
  MATLAB's layout, where the sum is at 0.91 to 1.17.
- To be judged end to end: the MATLAB-layout sum at one thread (1.17) and the slice sums at 24
  threads (1.18 to 1.25).
