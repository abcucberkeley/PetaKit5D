# Conventions of the Python port

This is the one place that states the array order, the box conventions and MATLAB's class
rules (sections 1 and 2 follow the user's decisions of 2026-11-11/12: `(z, y, x)` C arrays and
0-based, end-exclusive boxes). The code in `petakit5d.core` implements them; every area codes against it. Each rule
below was measured on MATLAB R2024b (Linux x86-64) at the oracle commit; the oracle cases in
`tests/core/f02_oracle_cases.py` (with references in `tests/core/oracle_f02.json`) keep them
checked.

## 1. Arrays

- **A volume is a NumPy array `(z, y, x)`, C-contiguous**: x is fastest in memory, as in a tiff
  page. An image is `(y, x)`. Extra dimensions lead: a stack of tiles is `(tile, z, y, x)`.
  Allocate with plain `np.zeros` / `np.empty` (C order); check at a public boundary with
  `arrays.check_volume` (dimensions, class, C order; it never converts and names MATLAB's
  layout when it sees a Fortran array), convert with `arrays.as_volume` (a view when the layout
  is right; an image becomes `(1, ny, nx)`).
- **MATLAB's array of the same data is `(y, x, z)` in Fortran order, and it is not the same
  bytes**: every z plane is transposed. `arrays.from_matlab` maps MATLAB's `(y, x, z, d4, ...)`
  to ours, `(..., d4, z, y, x)`:
  `np.ascontiguousarray(np.transpose(matlab, (n-1, ..., 2, 0, 1)))` (3D: `(2, 0, 1)`; an image
  keeps its shape; MATLAB's 4th dimension becomes the first axis). `arrays.to_matlab` is the
  inverse, `np.transpose(python, (n-2, n-1, n-3, ..., 0))`, Fortran-contiguous. Both copy large
  arrays in several threads. As values, `from_matlab(a)` is MATLAB's `permute(a, [n:-1:3, 1, 2])`.
- **Only image arrays are mapped.** Matrices and vectors of parameters (coordinate tables,
  4 x 4 transforms, lists of sizes, boxes in MATLAB files) keep `M(i, j) = P[i - 1, j - 1]`; a
  1-D array stands for a row (`1 x n`).
- **MATLAB's linear order** (y fastest, then x, then z) is the C order of
  `arrays.matlab_order(vol)`, the view `np.swapaxes(vol, -1, -2)` with axes `(z, x, y)`:
  `matlab_order(vol).ravel()` is MATLAB's `vol(:)`, and "first in MATLAB's order" means first in
  that view. The voxel `(z, y, x)` has MATLAB's 0-based linear index `y + ny*(x + nx*z)` (C
  offset `x + nx*(y + ny*z)`); `arrays.zyx_to_matlab_index` and `matlab_index_to_zyx` convert.
- **Sizes**: our shape `(nz, ny, nx)` is MATLAB's `[ny nx nz]`; `arrays.size_to_matlab` applies
  MATLAB's size rules on the way (no trailing singletons: a volume with nz = 1 is `[ny nx]`),
  `size_from_matlab(msize, ndim)` goes back.
- **Implicit expansion is NumPy's broadcasting** under this mapping: `(nz, ny, nx)` with
  `(ny, nx)`, `(ny, 1)`, `(nx,)` or `(nz, 1, 1)` expands as MATLAB expands `ny x nx x nz` with
  `ny x nx`, `ny x 1`, `1 x nx` or `1 x 1 x nz` (checked on 12 shape pairs against MATLAB).
- **Axes**: MATLAB's dimension 1 (y) is axis -2, dimension 2 (x) axis -1, dimension 3 (z) axis
  -3, dimension k > 3 axis -k (`axes.dim_to_axis`, `axis_to_dim`): MATLAB's `max(a, [], 3)` is
  our `vol.max(axis=-3)`.
- **Order-dependent steps keep MATLAB's order inside**: each output voxel keeps MATLAB's
  arithmetic, and a kernel may visit voxels in any order, except for sums and other
  accumulating reductions (MATLAB's linear order), "first" rules and label numbering,
  recurrences along an axis (same axis, direction and order of passes: y, then x, then z),
  position-dependent arithmetic (keyed by the MATLAB index, not the memory offset), the exact
  FFT, and batch numbering. Read strided lines or use one scratch buffer in MATLAB's layout (a C
  array with axes `(z, x, y)`); the speed gates count these copies.
- **MATLAB's size rules on MATLAB-form values** stay as functions: `msize`, `mndims`, `mnumel`,
  `msqueeze`, `sub2ind`, `ind2sub` (1-based, column-major), `mpermute` (MATLAB's `permute`), and
  `axis_order_mapping` on MATLAB's axis strings. They are for parameters, `.mat` contents and
  the MATLAB-name layer, not for our arrays.
- **The layout rule is a test** (`tests/core/test_layout_rules.py`): outside `core`, `io` and
  `testing`, no `np.asfortranarray` and no `order="F"` in an array creation or copy
  (`reshape`/`ravel` with `order="F"` remain for MATLAB's element order of values); no
  `f_contig` arrays in C++ outside `cpp/common` and `cpp/io`. Kernels take
  `pk5d::cvolume<T>` / `cimage<T>` (C-contiguous `(nz, ny, nx)` / `(ny, nx)`, bound with
  `.noconvert()`); `fvolume` / `fimage` stay in `common.hpp` until the last module has moved.

## 2. Boxes, ranges, per-axis values and batch geometry

- **A box is `(z0, y0, x0, z1, y1, x1)`**, 0-based, end-exclusive: it selects
  `vol[z0:z1, y0:y1, x0:x1]`. A 2D box is `(y0, x0, y1, x1)`. A range is `(start, stop)`. Indices
  and offsets are 0-based. An open end is `None` (MATLAB's `Inf`).
- **MATLAB's values convert at the boundary**, with `bbox`: `box_from_matlab` /
  `box_to_matlab` (`[ymin xmin zmin ymax xmax zmax]`, 1-based, inclusive, is
  `(zmin-1, ymin-1, xmin-1, zmax, ymax, xmax)`; 2D `[ymin xmin ymax xmax]` is
  `(ymin-1, xmin-1, ymax, xmax)`), `boxes_from_matlab` / `boxes_to_matlab` for `(n, 6)` tables,
  `range_from_matlab` / `range_to_matlab` (`[first last]` is `(first-1, last)`). The rest work on
  our boxes: `box_slices`, `crop` (a view), `box_shape`, `full_box`, `check_box`, `clamp_box`,
  `intersect_boxes`, `pad_widths` (for boxes that reach outside, as the crop wrappers' pad option).
- **Per-axis vectors follow `(z, y, x)`**: sizes `(nz, ny, nx)`, block, batch and chunk sizes,
  pool sizes, borders and paddings, resample factors, voxel sizes given as a vector, shifts and
  offsets in voxels, per-axis flags (a z MIP is `(1, 0, 0)`). `axes.from_yxz` / `to_yxz` convert
  MATLAB's `[y x z]` parameters, `from_xyz` / `to_xyz` the `[x y z]` ones; a scalar applies to
  every axis and passes through; 2 values are 2D (`(y, x)`).
- **Exceptions, each by name**:
  1. scalars that name their axes keep their names: `xy_pixel_size`, `dz`, `skew_angle`,
     `xy_max_offset`, `z_max_offset`;
  2. stage and physical coordinates stay x, y, z by name, as in the image lists, and the
     `axis_order` strings (`'x,y,z'`, `'-x,y,z'`) are unchanged;
  3. axis strings that name an array's axes (`data_order`, `input_axis_order`,
     `output_axis_order`) name them in our order: the Python defaults are `'z,y,x'` and `'zyx'`
     where MATLAB's are `'y,x,z'` and `'yxz'` (`axes.axis_string_from_matlab` /
     `axis_string_to_matlab`; `axes.axis_permutation(src, dst)` gives the `np.transpose` axes);
  4. transform matrices built for `imwarp` and the deskew kernels (the 4 x 4 `tmat`) keep
     MATLAB's `(x, y, z)` convention inside the code: their entries are compared with MATLAB's;
  5. files in MATLAB's formats keep MATLAB's values and orders (`parameters.json`, `.mat` files,
     stitch and block information, flag files, image lists): convert at the file boundary.
- **The MATLAB-name layer** (`petakit5d.matlab`, the command line, cluster and dataset configs)
  keeps MATLAB's values and converts through the parameter registry; code inside the port only
  sees the Python forms.
- **Batch and block coordinates** of the large-file paths come from `blocks`, in our forms and
  in MATLAB's batch order and numbering:
  - `zarr_chunk_coordinates` (`XR_zarrChunkCoordinatesExtraction`): batches of whole blocks
    (the batch size is first limited to the region and rounded up to whole blocks), numbered
    with y fastest, then x, then z; borders limited to the volume; with `same_batch_size` the
    last batch on an axis starts earlier so that every batch has the same size before borders.
    Returns `(batch_boxes, region_boxes, local_boxes)` as `(n, 6)` arrays of our boxes.
  - `sub_volume_coordinates` (`XR_subVolumeCoordinatesExtraction`): an `(n, 6)` array of boxes,
    x fastest, then y, then z (its `meshgrid` order), sizes made good FFT sizes.
  - `bbox_to_blocks` (`bboxToBlocks`): the blocks a box touches, in MATLAB's order; block
    subscripts `(z, y, x)` 0-based and `block_ind` MATLAB's linear block number minus 1.
  - `find_good_factor_number`: sizes with prime factors 2, 3, 5, 7 only, odd multiples of 7
    skipped by default.
  The `*_matlab` forms (`zarr_chunk_coordinates_matlab`, ...) take and return MATLAB's values
  for files in MATLAB's format; they are the reference the Python forms convert around (the
  2,100 parameter sets of the oracle cases pass through the conversions).

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
  `double(x)` round to nearest (beyond the single range: Inf). **`single(x)` of int64 and uint64
  rounds twice**, to double and then to single (`dtypes.single`, `to_class`): `single(int64(2^53 +
  2^29 + 1))` is 2^53, where `cast(x, 'single')` rounds once and gives 2^53 + 2^30. Measured on
  arrays and on scalars passed to a function file or an `eval`'d anonymous function; `double(x)`
  rounds once.
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
    char argument is consumed whole by `%s` and one character at a time by the others. `%s` of
    a numeric argument takes the run of character codes that starts at its current element as
    one field (`sprintf('[%s]', [65 66 67])` is `[ABC]`, `[%05s]` of `[65 66]` `[000AB]`,
    `[%.1s]` of `[65 66 67]` `[A]`); an element that is no code (a fraction, a negative value,
    one above 65535, NaN, Inf) goes alone, with `%e` (`[%s]` of `[65 66.5 67]` is
    `[A][6.650000e+01][C]`). Logical values are codes. Text fills its width with zeros under
    the `0` flag, on the right with `-` (`%05s` of `'ab'` is `000ab`, `%-05s` `ab000`).
  - `Inf`, `-Inf`, `NaN` print as such (width and `-` kept, `0` ignored, `+` only for Inf).
  - Arrays are consumed column-major; the format repeats while data remain and stops at the
    first conversion without data (with no data argument at all, at the first conversion).
  - An empty argument (`[]`, `''`, also among other arguments) holds the place of one
    conversion, which prints the field of empty text: `sprintf('a[%d]b%d', [], 5)` is
    `a[]b5`, `sprintf('a%db%dc', [])` `ab`, `sprintf('[%5d]', [], 1)` `[     ][    1]`, and
    `[%05d]` of `[]` `[00000]`. A `*` width or precision taken from an empty argument ends the
    output there.
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
