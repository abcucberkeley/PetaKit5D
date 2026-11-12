# Testing against MATLAB: the oracle harness

MATLAB PetaKit5D is the reference of the port. The oracle harness runs MATLAB and Python on the
same inputs and compares their results, bit for bit unless a case declares a tolerance. Its
parts:

| Part | Where | What |
|---|---|---|
| cases | `petakit5d.testing.cases` | what MATLAB runs, what Python runs, their inputs, how to compare |
| oracle | `petakit5d.testing.oracle` | the golden store, MATLAB packs, the repo manifests, the pytest check, the command line |
| compare | `petakit5d.testing.compare` | MATLAB values in Python, value files, comparison with readable diffs, file trees |
| MATLAB side | `tests/oracle_matlab/+pkpy_oracle` | runs a pack of cases in one MATLAB session and publishes each result |
| oracle commit | `tests/oracle_matlab/ORACLE_COMMIT` | the MATLAB commit the repo's references were made with |

Three places hold MATLAB's results:

- **the golden store** (a shared folder on the cluster, `PETAKIT5D_GOLDEN_DIR`): every output
  in full, with the files a function wrote, its warnings and error. Written only by the harness.
- **manifests** in the repo, next to each test module (`test_x.py` → `test_x.oracle.json`): per
  case the SHA-256 of every output, its class and size, and small arrays and texts inline. CI
  has no MATLAB and no store; it checks the cases against the manifests (bit for bit by hash;
  the readable element-by-element diff of a large value needs the store).
- nothing else: a test never stores references of its own.

## Writing cases

```python
import numpy as np

from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, Work
from petakit5d.testing.compare import tol, tree


def crop(im, bbox):  # the Python side (normally the port's own function)
    b = [int(x) for x in bbox]
    return im[b[0] - 1 : b[3], b[1] - 1 : b[4], b[2] - 1 : b[5]]


CASES = [
    Case(
        "crop_uint16_odd",  # unique in the module; no spaces or slashes
        matlab="crop3d_mex",  # the MATLAB function (feval)
        args=[In("im"), [2, 3, 1, 9, 10, 7]],
        python=crop,  # or "petakit5d.tools.crop:crop_3d"
        inputs={"im": Value(lambda: make_volume((10, 12, 8), np.uint16))},
    ),
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)
```

`oracle.parametrize` makes one test per case (pytest id = case id), marks them `oracle` (and
`gpu`, `demodata` or the case's `marks`), and registers the cases for the command line.
`oracle.check` finds the reference, runs the Python side in a fresh work folder under
`tmp_path`, compares, and fails the test with the report.

### Case fields

| Field | Meaning |
|---|---|
| `id` | unique in its module |
| `matlab` | MATLAB function name, called with `feval` (`pkg.fn` for package functions) |
| `code` | instead of `matlab`: MATLAB code that reads `in` (cell of the arguments) and `work`, and sets `out` (cell of the outputs) |
| `args` | MATLAB arguments: Python values (rules below) and placeholders |
| `nargout` | number of MATLAB outputs (default 1); `outputs` names them (default `out1`, ...) |
| `python` | the Python callable, or `"package.module:function"`; None: a reference-only case |
| `py_args`, `py_kwargs` | Python arguments (default: `args`); placeholders work here too |
| `py_out` | maps the Python result to the list of `nargout` outputs (default: the result, or the tuple it returns) |
| `inputs` | name → `Value`, `File`, `Tree` or `Data` (below) |
| `compare` | a spec for every output, or `{output name or index: spec}` (default exact) |
| `files`, `tree` | compare the files written in the work folder (default: when the case uses `Work`, `InPath` or file inputs); `tree(...)` sets how |
| `expect_error` | MATLAB must raise (a string: with that identifier); otherwise a MATLAB error fails the test as a broken case |
| `py_error`, `error_message` | the exception classes Python may raise when MATLAB raises; compare messages too |
| `warnings` | `'record'` (default) or `'ids'`: Python must give warnings with the same identifiers (set `identifier` on the warning). MATLAB's warnings are read from its command window text (kept as the result's `stdout.txt` and parsed again when a reference is read): identifiers are exact; messages as displayed, long ones wrapped at the window width, the last one exact. A warning issued with verbose off (MATLAB prints no "Type warning off" line) has an empty identifier unless it is the last one, which gets `lastwarn`'s |
| `helpers` | folders of MATLAB helper files (relative to `python/`) put on MATLAB's path for the pack; part of the key |
| `device` | `'cpu'` (default), `'gpu'` or `'gpu:<model>'` (the MATLAB GPU's name must contain it) |
| `marks` | more pytest marks (`slow`, ...) |

### Inputs and placeholders

Inputs are made by deterministic Python functions and written once; MATLAB reads them from the
store, Python makes them again (in CI too) and the case key checks that the bytes are the same.
**A generator must give the same bytes every time**: no times, uuids, random seeds from the
clock or library versions that change the output. `oracle run` makes each input twice and
refuses to run MATLAB when they differ (`scipy.io.savemat` writes the time into the header, for
example: overwrite it).

| Input | Made by | MATLAB gets | Python gets |
|---|---|---|---|
| `Value(make)` | `make()` returns an array or value | the value (`In(name)`) | a copy of what `make` returned |
| `File(make, name)` | `make(path)` writes the file | a copy in the work folder named `name` (`InPath(key)`) | the same |
| `Tree(make, name)` | `make(folder)` fills the folder | a copy in the work folder named `name` (`InPath(key)`, with a trailing `/`) | the same |
| `Data(path, name, copy=True, content=False)` | existing data, `$VARIABLES` expanded | a copy in the work folder, or the path itself (`copy=False`) | the same |

`InPath(key)` names an input by its key in `inputs` and gives the path of its copy, which has
the input's `name` (`inputs={"m": File(make, "m.mat")}`: `InPath("m")` is `<work>/m.mat`), on
both sides.

`Work(sub)` is the case's work folder with a trailing `/` (and `sub`). Both sides run with the
work folder as the current folder; the files a call writes or changes there are its file
output. Paths in outputs, texts and files are compared with each side's work folder masked as
`<WORK>`.

### Python values as MATLAB arguments

| Python | MATLAB |
|---|---|
| `None` | `[]` (0x0 double) |
| `bool` | logical scalar |
| `int`, `float` | double scalar (an int beyond 2^53 is refused: use a NumPy scalar) |
| `complex` | complex double |
| `str` | char row (`''` is 0x0) |
| NumPy array or scalar | its class and shape; 1-D is a row, 0-D is 1x1; `<U1` 2-D arrays are char matrices |
| list or tuple of numbers / of bools | double / logical row |
| other lists and tuples | cell row |
| dict, dataclass, `SimpleNamespace` | 1x1 struct |
| `mcell(...)`, `mcol(values, cls)`, `mtyped(x, cls)`, `mstring(...)`, `mstruct(...)`, `mempty(cls, size)` | a cell, a column, a typed array, a string array, a struct, an empty array |

Arrays are `(y, x, z)` with MATLAB's sizes; give Fortran-ordered arrays where the port uses them
(values are compared by index, not by memory layout).

### Reading Python results

A Python result is read as the MATLAB reference where Python has no counterpart, and only
there: a 1-D array or list matches a row or a column; a Python int takes the reference's
integer class if it fits; `''`, `None` and empty 1-D arrays match MATLAB's empty sizes; a list
of dicts is a struct array; a dict is a struct or a `containers.Map`; a value compared with a
MATLAB object (a saved `inputParser`) is compared with the object's value, with a note. Classes
are never changed otherwise: a float64 result for a single reference is a class difference.

## Comparisons

- `exact()` (the default): class, size and every bit of every value, NaN payloads and the sign
  of zero included (MATLAB's NaN is `0xFFF8000000000000`; `np.nan` is `0x7FF8000000000000`).
- `tol(reason=..., abs=, rel=, rel_max=, ulp=, rel_l2=, max_frac=, max_count=)`: a declared
  tolerance; `reason` names the approved exception in `docs/parity/<area>.md`. A differing value
  is within it if its difference is at most the largest of `abs`, `rel` times the reference value
  and `rel_max` times the largest reference magnitude, or at most `ulp` units in the last place
  away; at most `max_count` values and a fraction `max_frac` may differ at all; the relative L2
  error is at most `rel_l2`. Class and size must still match; NaN must match NaN.
- `ignore(reason)`: not compared (a function handle, a time).
- `Custom(fn)`: `fn(ref, act)` returns a list of problems.

A report lists class, size and kind differences first, then for values: how many differ, the
first index (0-based, column-major, and as MATLAB subscripts), the next ones, the largest
absolute and relative difference, a histogram (counts for integers, ULP distance for floats),
NaN and signed-zero categories, and hints: values in C order, axes permuted, data shifted by
one along an axis, a constant offset, a missing last slice. Nested values are named by path
(`out2(3).field{2}`).

### Files

`tree(values=, per_file={glob: spec}, ignore=[regex], ignore_vars=('t', 't1'), ignore_fields=(),
tiff_bytes=True, zarr_chunks='note', text_format=True, masks=[(regex, replacement)])`:

- names are matched with uuids and 14-digit time stamps masked; files only one side wrote,
  and files removed from the work folder, are reported;
- **tiff**: equal bytes, else the decoded voxels are compared and a byte difference with equal
  voxels is a layout problem (`tiff_bytes=False`: a note);
- **zarr**: decoded voxels and the metadata texts (`.zarray`, `.zattrs`, ...); different
  chunk-file sets with equal data are a note (`zarr_chunks='strict'`: a problem), since the
  Python writer skips all-zero chunks;
- **.mat**: the variables, decoded (MATLAB's own reading of its file; the Python side through
  SciPy for v4-v7 and h5py for v7.3), without `ignore_vars`;
- **.json** and text files: the text with masks applied (uuids, times, the work folder, the
  host); for JSON, a text difference with equal content is a format problem
  (`text_format=False`: compared by content);
- **png** without its time chunks; anything else byte for byte.

## Making references: `orun.sh` and `python -m petakit5d.testing.oracle`

References are made on the cluster, from the oracle tree (an immutable copy of the reference
commit), never from a worktree. The cluster helper `orun.sh` (outside the repo, with the
cluster's own paths) snapshots your `python/`, submits one CPU job (a GPU job with `-g`) and runs

```
python -m petakit5d.testing.oracle run [--manifest] [--dry-run] <pytest selection>
```

there. `run` collects the cases with pytest (paths, node ids, `-k`, `-m` work as in pytest),
checks that their inputs are deterministic, puts the inputs in the store, and runs every case the
store lacks for the current oracle commit in **one MATLAB session** (a pack); each result is
published as soon as it is complete. With `--manifest` it writes the manifests of the selected
modules (`orun.sh -M` copies them into your worktree). Other commands: `manifest` (write
manifests from the store), `status` (store and manifest state of the selected cases),
`clean-staging --hours N` (remove leftovers of stopped sessions).

The session: `cd` to the oracle tree, `setup([])`, the pack's copy of `+pkpy_oracle` (and the
cases' `helpers`) on the path, then `pkpy_oracle.run_pack`. Each case starts as a MATLAB test
does (`rng(0, 'twister')`, the session's warning state, the decon OTF caches cleared, figures
closed afterwards) in a fresh work folder on local disk. MATLAB needs 4 or more CPUs (it fixes
its FFT mode at the first FFT). A session starts in about 11 to 16 s; then a small case costs
about 8 ms beyond its own compute (median; the first case of a session about 1 s more).

Typical loop:

1. write cases; `orun.sh -a <task> tests/<area>/test_x.py` makes the missing references;
2. run the tests on the cluster (`prun.sh ... 'python -m pytest tests/<area>'`): they compare
   with the store;
3. `orun.sh -a <task> -M tests/<area>/test_x.py` writes the manifest; commit it with the tests.

A stored reference is never replaced: a reference changes only with a new oracle tree (a new
`ORACLE_COMMIT`; the orchestrator reruns the packs and reviews what changed) or with a new case
key format (`FORMAT` in `cases.py`, bumped when the meaning of a key or a stored result changes).

## How a case is identified

The **key** is the SHA-256 of what determines MATLAB's result: the key format, the function or
code, `nargout`, the arguments (values by their hash, placeholders by name), the content hash
of every input and its name in the work folder, whether files are captured, the device, the
helper folders' contents and, for the harness's own `pkpy_oracle.*` functions, their file. The
Python side, the comparison and the case id are not part of it. Results are stored by oracle
commit and key; the manifest keeps the key, so CI notices a case or input that changed since its
reference was made, and every entry names its oracle commit, so a reference made with another
commit is flagged (`check_manifests`, and `test_repo_manifests_use_the_current_oracle`).

## The golden store and concurrency

```
<store>/inputs/<hh>/<hash>.pkv|.file|.tree     inputs, by content hash
<store>/<commit7>/<hh>/<key>/                  one case: meta.json, outputs.pkv, files/, decoded/, stdout.txt
<store>/<commit7>/.staging/<pack>/<key>/       a result being written (never read)
```

Each result is written into the session's own staging folder and published with **one
rename** onto its final name. A rename onto an existing result fails, so of two sessions that
run the same case the first wins and the second drops its copy (`lost_race`); a stopped session
leaves only its staging folder, which nothing reads. Inputs are written to a temporary name and
renamed the same way. Readers look only at published folders (`meta.json` is in every one).

## Value files (`.pkv`)

The exchange format between Python and MATLAB (`write_pkv`/`read_pkv`,
`pkpy_oracle.write_value`/`read_value`): 8 bytes `PKV1\0\0\0\0`; the raw data of every array,
each at an 8-byte aligned offset, column-major, little-endian (char as UTF-16 code units,
logical as one byte, complex as the real block then the imaginary block); a UTF-8 JSON
description; its length as uint64; 8 bytes `PKV1END\0`. JSON nodes:

| kind | node |
|---|---|
| array | `{"k":"a", "c": class, "s": size, "o": offset, "n": bytes, "x": complex}` |
| cell | `{"k":"cell", "s": size, "v": [nodes, column-major]}` |
| struct | `{"k":"struct", "s": size, "f": [fields], "v": [nodes: element by element, field by field]}` |
| string | `{"k":"string", "s": size, "v": [text or null]}` |
| sparse | `{"k":"sparse", "c": class, "s": size, "i": rows, "j": columns (1-based doubles), "v": values}` |
| function handle | `{"k":"func", "t": text}` |
| containers.Map | `{"k":"map", "kt": key type, "vt": value type, "keys": cell, "values": cell}` |
| object | `{"k":"object", "c": class, "v": its value (inputParser: Results) or null, "t": display text}` |
| placeholder | `{"k":"ph", "t": "work" or "path", "n": name, "sub": text}` (arguments only) |
| reference | `{"k":"ref", "p": path of another value file, "h": its hash}` (arguments only) |

Hashes (`vhash`): array `sha256("pk1:a:<class>:<size>:<r|c>:" + raw data)`, size written as
`10,12,8`; cell `sha256("pk1:cell:<size>:" + child hashes joined by ";")`; struct
`sha256("pk1:struct:<size>:<sorted fields>:" + per element the children's hashes in sorted
field order joined by ",", elements joined by ";")` (field order is not part of it; a report
notes it); other kinds similarly (see `vhash`). Manifests keep descriptors: kind, class, size,
hash, and the data of arrays up to 2 KB (base64; char rows as text) and texts up to 4 KB;
elements whose data is kept carry no hash of their own (it is computed again). A cell, struct,
string array or map whose descriptor would take more than 8 KB (compact JSON) keeps only its
kind, size, struct fields or map types, and hash: CI checks it bit for bit by hash and reports
kind, size and field differences; the element-by-element report needs the store. A manifest has
one line per case (a changed case is one changed line); older layouts are read too.

## Environment variables

| Variable | Meaning |
|---|---|
| `PETAKIT5D_GOLDEN_DIR` | the golden store (cluster); without it tests compare with the manifests |
| `PETAKIT5D_ORACLE_DIR` | the oracle tree; its commit must equal `ORACLE_COMMIT` |
| `PETAKIT5D_ORACLE_COMMIT` | overrides the oracle commit |
| `PETAKIT5D_MATLAB` | the MATLAB command (default `matlab`) |
| `PETAKIT5D_ORACLE_WORK` | where MATLAB's work folders go (default: the temp folder) |
| `PETAKIT5D_ORACLE_MISSING` | `skip`: cases without a reference are skipped (default: they fail) |
| `PETAKIT5D_ORACLE_REPORT_DIR` | write a JSON report of every compared case there |
| `PETAKIT5D_ORACLE_HARNESS` | the MATLAB side (default `tests/oracle_matlab`) |
| `PETAKIT5D_ORACLE_TEST_STORE` | the harness's own tests make their private stores there (to test on network storage) |

## Troubleshooting

- *no MATLAB reference*: run `orun.sh` for the module, then `orun.sh -M`.
- *the case or its inputs changed since its reference was made*: rerun `orun.sh -M`; if you
  changed nothing, an input generator is not deterministic on this platform.
- *reference made with oracle X*: the manifest predates the current oracle tree.
- *MATLAB raised ...; set expect_error=True*: the case calls MATLAB wrongly, or tests an error.
- a session that crashes on a case stops the pack; the cases before it are published and the
  log is the pack's `matlab.log`. The next run runs the rest, that case included: select it
  alone (`-k`) to debug it.
- MATLAB's launcher is a shell script: the harness passes nothing but the pack's path on its
  command line (case code goes through `pack.json`).
