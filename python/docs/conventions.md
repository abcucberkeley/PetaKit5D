# Conventions of the Python port

This is the one place that states the array order, the box conventions and MATLAB's class
rules. The code in `petakit5d.core` implements them; every area codes against it. Each rule
below was measured on MATLAB R2024b (Linux x86-64) at the oracle commit; the oracle cases in
`tests/core/f02_oracle_cases.py` (with references in `tests/core/oracle_f02.json`) keep them
checked.

## 1. Arrays

- **Volumes are `(y, x, z)`, Fortran-ordered** (y fastest in memory), as in MATLAB and in the
  zarr files PetaKit5D writes (`"order": "F"`, shape `[ny, nx, nz]`). A 2D image is `(ny, nx)`
  and counts as a volume with `nz = 1`.
- **Linear indices are MATLAB's minus one**, in column-major order: `a.ravel(order="F")[k - 1]`
  is MATLAB's `a(k)`. `arrays.sub2ind` and `arrays.ind2sub` take and return MATLAB's 1-based
  values. Ties such as "the first maximum" are first in this order (`matlabcompat` owns them).
- **A 1-D NumPy array stands for a row vector** (`1 x n`).
- **MATLAB's size rules** (`arrays.msize`, `mndims`, `msqueeze`): trailing singleton dimensions
  beyond the second don't exist (`zeros(3, 4, 1)` is `3 x 4`); `size(a, 3)` of an image is 1;
  `squeeze` leaves 2D arrays alone (`1 x 5` stays) and keeps 2 dimensions (`1 x 1 x 5` gives
  `5 x 1`).
- **At a public boundary** check with `arrays.check_volume` (class, dimensions, Fortran order;
  it never converts) or convert with `arrays.as_volume` (a view when the layout is right).
  Allocate with `arrays.fzeros`, `fempty`, `fones`, `ffull` (MATLAB class names accepted).
- **Tiff pages** come from cpp-tiff as C-ordered `(z, y, x)`: `arrays.zyx_to_yxz` and
  `yxz_to_zyx` are the reference conversions (IO1 owns the fast one).
- **Implicit expansion** aligns leading dimensions in MATLAB (`3x4x5 .* 3x4` expands along z);
  NumPy aligns trailing ones. `dtypes.plus`, `times`, ... expand as MATLAB does.
- `arrays.mpermute(a, order)` is `permute` with a 1-based order (a Fortran-ordered copy);
  `axis_order_mapping(order, 'yxz')` is the permutation that stores yxz data in `order`.

## 2. Bounding boxes and batch geometry

- **Boxes in parameters keep MATLAB's values**: `[ymin, xmin, zmin, ymax, xmax, zmax]`,
  1-based, inclusive. Configs, GUI settings and the oracle pass the same numbers to both
  implementations. Code converts with `bbox.bbox_to_slices` (0-based, half-open slices) and
  `bbox.crop` (a view), never by hand; `slices_to_bbox`, `bbox_size`, `full_bbox`,
  `clamp_bbox`, `bbox_intersect` and `pad_widths` (for boxes that reach outside, as the crop
  wrappers' pad option) cover the rest. A 2D image takes a box with `zmin = zmax = 1`.
- **Batch and block coordinates** of the large-file paths come from `blocks`, with MATLAB's
  values and order:
  - `zarr_chunk_coordinates` (`XR_zarrChunkCoordinatesExtraction`): batches of whole blocks
    (the batch size is first limited to the region and rounded up to whole blocks), numbered
    column-major (y fastest); borders limited to the volume; with `same_batch_size` the last
    batch on an axis starts earlier so that every batch has the same size before borders.
  - `sub_volume_coordinates` (`XR_subVolumeCoordinatesExtraction`): x fastest, then y, then
    z (its `meshgrid` order), sizes made good FFT sizes.
  - `bbox_to_blocks` (`bboxToBlocks`): the blocks a box touches, column-major.
  - `find_good_factor_number`: sizes with prime factors 2, 3, 5, 7 only, odd multiples of 7
    skipped by default.
  All return int64 values, 1-based; inputs must be integral.

## 3. Classes and arithmetic (`petakit5d.core.dtypes`)

NumPy differs from MATLAB in each of these: `np.round` rounds half to even, `astype` truncates
and wraps, a float64 scalar or array promotes a float32 array, and a Python float meets a
float32 array as a float32 (rounded before the operation). Never use them where MATLAB code
has `round`, a class conversion or mixed-class arithmetic; ruff bans `np.round`, `np.rint`,
`np.around`, `np.percentile` and `np.median` outside `core`, `matlabcompat` and the tests.

### 3.1 Conversions to integer classes

| MATLAB call | Python | Rounding |
|---|---|---|
| `cast(x, 'uint16')`, `feval('uint16', x)`, `z(:) = x` into a uint16 array, `cast(x, 'like', ...)` | `dtypes.cast(x, 'uint16')` | half away from zero, exactly |
| `uint16(x)`, `int16(x)`, `uint8(x)`, `int8(x)` on an array (2 or more elements), also through `@uint16` and `@(v) uint16(v)` | `dtypes.to_class(x, cls)` or `dtypes.uint16(x)`, ... | double data: exact; **single data: `floor(abs(x) + 0.5)` in single precision**, sign restored |
| the same class functions on a scalar the code knows as one (a local variable, `x(i)`, `x(end)`, a struct field, a cell element, a loop body, an anonymous function in a script) | `dtypes.uint16(x)` (exact for scalars) | exact |
| a scalar that is the argument of a function file (`function y = f(x)`, `uint16(x)` inside), or one in an anonymous function made by `eval`'d code | `dtypes.to_class(x, cls, scalar_fast=True)` | as arrays |
| `int32`, `uint32`, `int64`, `uint64` of anything | `dtypes.int32(x)`, ... | exact |

The scalar rows are how MATLAB's execution engine compiles the call, not a property of the
value: the same scalar rounds either way depending on the code around it. Only `0.49999997`
and `-0.49999997` are affected; a port follows the array rule for arrays and the row that
matches the MATLAB code for scalars.

- All of them saturate (`uint16(-3)` is 0, `uint16(70000)` 65535, `int8(Inf)` 127) and map NaN
  to 0. 64-bit targets saturate at the exact limits (`int64(2^63)` is `intmax('int64')`).
- **The single-array case differs from `cast` at one value**: `floor(0.49999997 + 0.5)` is 1
  in single precision, so MATLAB's `uint16(x)` gives 1 (and `int16(-x)` gives -1) where
  `cast(x, 'uint16')` gives 0. Measured: in a random 2,000,000-element single array, the two
  differ exactly at the 2,000 elements equal to `0.49999997`, and `uint16(x)` equals
  `uint16(floor(x + 0.5))` everywhere. No other single value and no double value differs (the
  sum `|x| + 0.5` is exact for every other value that does not saturate in these classes).
- `logical(x)`: nonzero is true; NaN raises (`MATLAB:nologicalnan`). `single(x)` and
  `double(x)` round to nearest (beyond the single range: Inf).
- `dtypes.mround` is `round`: half away from zero, class kept, and `-0` kept (`round(-0.3)` is
  `-0`).

### 3.2 Single with double

- **A single array with a double scalar or a double array is computed in double and rounded to
  single** (`single(double(a) op b)`), for `+ - .* ./ .\`, scalar or array on either side
  (measured on 10^6 values with 8 scalars and with an array: all equal to the double-then-round
  result; NumPy's float32 arithmetic with the scalar rounded first differs in most of them).
  `dtypes.plus`, `minus`, `times`, `rdivide`, `ldivide` do this. Single with single is single
  arithmetic (the same bits as computing in double and rounding).
- **`eps` is the double `eps`** (`dtypes.EPS`, 2^-52), also next to single data:
  `x + eps` with single `x` is computed in double and rounded. `dtypes.eps(x)` is `eps(x)`
  (`eps('single')` is 2^-23).
- Single with logical gives single; logical with logical gives double.

### 3.3 Integer classes

- **Integers combine only with the same class or a double scalar** (`MATLAB:mixedClasses`
  otherwise, also with a double array, any single value and logical); char counts as double.
- **8 to 32-bit results are computed in double, then rounded half away from zero and
  saturated**: `uint16(5)/2` is 3, `uint16(7)/2` is 4, `uint16(0)/0` is 0, `uint16(5)/0` is
  65535, `-int8(-128)` is 127, and `uint16(1) + 0.49999999999999994` is 2 (the double sum is
  1.5).
- **64-bit results are exact** (rational), then rounded and saturated:
  `int64(9007199254740993) * 1.5` is 13510798882111490. MATLAB has one measured exception
  (3.5), which `dtypes` reproduces.
- `dtypes.idivide(a, b, opt)` is `idivide` ('fix' by default, 'round', 'floor', 'ceil').

### 3.4 NaN

- **MATLAB's NaN has the sign bit set** (`0xFFF8000000000000`, the x86 default NaN); Python's
  `float('nan')` and `np.nan` do not. Where MATLAB code writes `NaN`, Python writes
  `dtypes.NAN`; `str2double` returns it for text that is not a number.
- A NaN operand keeps its payload: the result is the first NaN operand, quieted (MATLAB on
  x86; NumPy may swap the operands of `+` and `*`). The `dtypes` operations follow this.

### 3.5 MATLAB bugs kept

- **int64 with `-0.49999999999999994`**: MATLAB R2024b adds a positive int64 integer and
  `-0.49999999999999994` (the double below -0.5's magnitude) to one more than the exact
  rounding: `int64(5) - 0.49999999999999994`, `int64(5) + (-0.49999999999999994)`,
  `-0.49999999999999994 + int64(5)` and `-0.49999999999999994 - int64(-5)` are 6 (exact: 5).
  No other value of a 31-value near-half grid differs, nor uint64, `.*` or `./`. `dtypes`
  reproduces it.
- `writeJsonFile` passes the JSON to `fprintf` as a format (4); `int2str` of small values
  errors (4); `simplifyPath` of a missing path is `/` and `mkdir_recursive` of a relative path
  starts at the root (5); `readTextFile` of an empty file is `{-1}` (4).

## 4. Numbers as text (`petakit5d.core.text`, `jsonenc`)

The wrappers hand parameters to their tasks as text, and folder names and parameter files hold
numbers as text. Apply MATLAB's formatting at the same boundary:
`text.through_text(value, '%.10f')` is `str2double(sprintf('%.10f', value))`, and
`through_mat2str` is the trip through `mat2str` (15 significant digits).

- `sprintf` (`text.sprintf`): C formatting, plus MATLAB's rules:
  - `%d` and `%i` print integer values in `[-2^63, 2^63)`; other values switch to `%e` with the
    same flags, width and precision (`%.20d` of 0.108 is `1.07999999999999998779e-01`, of 32
    `00000000000000000032`). `%u`, `%x`, `%o` take `[0, 2^64)`; negative values switch to `%e`.
  - `%s` and `%c` of a number print the character for integers in `[0, 65535]`, else `%e`; a
    char argument is consumed whole by `%s` and one character at a time by the others.
  - `Inf`, `-Inf`, `NaN` print as such (width and `-` kept, `0` ignored, `+` only for Inf).
  - Arrays are consumed column-major; the format repeats while data remain and stops at the
    first conversion without data (with no data argument at all, at the first conversion;
    with only empty data, the format is printed once with empty conversions).
  - Escapes `\n \t \\ \" \a \b \f \r \v \xHH \NNN` are processed in the format; an unknown escape
    (`\q`, `\/`, `\u`) or an invalid conversion (`%"`, `%q`) ends the output there.
- `num2str`, `int2str`, `mat2str`, `mat2str_comma`: MATLAB's algorithms (`num2str(x)` uses
  `%<d+7>.<d>g` with `d = min(12, max(1, floor(log10(max|x|)) + 1)) + 4`, integers through
  `int2str`; `num2str(x, fmt)` formats each element alone and right-justifies the columns;
  `mat2str` uses `%.15g`, `true`/`false`, `zeros(r,c)` for empty matrices). MATLAB's own
  `int2str` errors on a vector whose values are all below 0.01; so does ours.
- `str2double`: whitespace and commas are dropped (`'1,000'` is 1000), hexadecimal reads
  (`'0x10'`), `Inf` and `NaN` in any case but not `Infinity`.
- **`jsonencode`** (`jsonenc.jsonencode`, compact or `pretty_print`): a double prints with
  `%.15g` if that reads back as the same double, else `%.17g` (single: `%.6g`, else `%.9g`, read
  back in single); decimal exponents -4 to 5 print plainly, others as `1E-7`, `1.0E+20`,
  `1.5E+6`; integer classes print in full; NaN and Inf are `null`. Vectors (either orientation)
  are flat arrays, a 1-element array is its value, matrices nest by rows. Strings escape `"`,
  `\` and control characters (`\b \t \n \f \r`, else `\u00xx`); `/` and non-ASCII are written
  as they are. PrettyPrint indents by 2 spaces with `"key": value`. The Python form of MATLAB
  values is in the module's docstring (dict for struct, list for cell, `None` for `[]`).
- **`writeJsonFile`** writes the JSON with `fprintf(fid, s)`, so the JSON is a format:
  `textfiles.write_json_file` does the same (a string with `%03d` truncates the file there;
  `\\` in the JSON becomes `\`). This is a MATLAB bug kept on purpose.
- `readTextFile` / `writeTextFile` (`textfiles`): UTF-8; lines end at `\n`, `\r\n` and a lone
  `\r`; an empty file reads as `[-1]`; a list is written with `\n` after every line.
- The field order of `inputParser`'s `Results` (and so of `parameters.json`): see the oracle
  case `inputParser.order` (MATLAB sorts the names; F03 builds parameter files in that order).

## 5. Paths (`petakit5d.core.paths`)

- `fileparts` is MATLAB's (`'.bashrc'` is name `''`, extension `'.bashrc'`; `'/a'` has folder
  `'/'`; on Windows MATLAB's own rules), not `os.path.splitext`.
- `matlab_dir` lists what `dir` lists, in MATLAB's order: names sorted by UTF-16 code units
  (uppercase before lowercase, `' sp'` and `'#h'` before `.` and `..`); a file's folder has
  symbolic links resolved, a folder's does not; broken links are listed.
- `simplify_path`: an absolute path, a folder without a trailing separator; a pattern that
  matches two or more entries gives their folder; a path that does not exist gives `/` (both
  kept from MATLAB).
- `dir_recursive`: every file (or folder) below a folder, as `<given>/<...>`, sorted as
  `unique` sorts; links to folders are followed. `mkdir_recursive` builds the path from the root
  even for a relative path (MATLAB's `strjoin({'', ...}, '/')`; kept).
- MATLAB code often joins `[folder, '/', name]` with folders that end with `/`;
  `with_trailing_sep` and `strip_trailing_sep` make those forms explicit.

## 6. Threads, memory, ids

- **Thread counts** come from the machine: `threads.available_cpus()` is the CPU affinity
  (Slurm's cpuset, taskset) limited by the cgroup CPU quota (v1 and v2), never an environment
  variable. Kernels take a thread count argument; `threads.resolve_threads(n)` turns a user's
  `n` (None or 0: all available) into it. Nothing sets global thread state.
- **`memory.get_system_memory()`** is `getSystemMemory`: GB (2^30 bytes) of the job's cgroup
  limit in a Slurm job (if the cgroup is the job's own), else `SLURM_MEM_PER_NODE` (unless
  `SLURM_MEM_PER_CPU` is also set), else `scontrol show job`'s `mem=`; `MemTotal` without
  Slurm; `hw.memsize` on macOS; NaN when unknown. Equal to MATLAB's value inside Slurm jobs on
  both node types of the lab cluster in every case measured.
- `ids.get_uuid()`: a random version-4 UUID (36 characters; 4 on Windows, as MATLAB keeps);
  never touches global random state. `ids.get_hostname()` prints `Hostname: <name> `.

## 7. How the rules were measured

Four probes (the forms of the calls, the positions of the rounded-up values, the formats) and
an oracle pack of about 700 cases (`tests/core/f02_oracle_cases.py`) ran in MATLAB R2024b from
the oracle tree with 4 computational threads. The pack covers: every class conversion of
200,000 mixed values plus the edge values of each class (halves and their neighbours, limits,
NaN, Inf), by direct call, `cast` and scalar forms; integer arithmetic tables of every class
with every class's edge values and with 20 double scalars, near-half grids and `idivide`;
single with double scalars and arrays on both sides; 27 `sprintf` formats, `num2str`,
`int2str`, `mat2str`, `mat2str_comma` and `jsonencode` of about 600 values per class, rows and
matrices; `jsonencode` of structures and `writeJsonFile`, `writeTextFile`, `readTextFile`
bytes; `dir` order, `dir_recursive`, `simplifyPath` and `fileparts` on a file tree with case,
punctuation, hidden files and links; and the three geometry functions on 800, 500 and 800
generated parameter sets.
