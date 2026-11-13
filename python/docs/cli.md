# The command lines: `pk5d` and `pk5d-master`

Two console scripts come with the package (`petakit5d.cli`):

| Command | For | Names and values |
|---|---|---|
| `pk5d <function> ...` | people and scripts | snake_case functions and options, JSON values in the Python API's form |
| `pk5d-master <function> ...` | the GUI, PyPetaKit5D, existing scripts | MATLAB names and MATLAB-style strings, as the compiled MATLAB master takes them |

Both reach the same registered functions (`docs/MATLAB_MAPPING.md`), with the checks and
conversions of the parameter registry (`petakit5d.params`).

## `pk5d`

```
pk5d list [--all]                      the registered functions (--all: and the unported ones)
pk5d doc <function>                    arguments, options, defaults, MATLAB names
pk5d <function> [value ...] [--option value ...]
pk5d parse <function> <strings ...>    how pk5d-master reads MATLAB-style strings (JSON)
pk5d master <function> <strings ...>   the same as pk5d-master
pk5d mcc-shim <path> [--command CMD]   write a run_mccMaster.sh that runs pk5d-master
```

A function is named by its Python name (`crop_dataset` or `petakit5d.tools.crop.crop_dataset`)
or by its MATLAB name (`XR_crop_dataset`). Values are JSON; text that isn't JSON is taken as
a string, so paths need no quotes (`'"123"'` is the string 123). Options are the snake_case
parameter names with `-` or `_`; an option without a value is `true`; `--` ends the options.

```
pk5d crop_dataset /mnt/projects/embryo_2024 '[0, 0, 0, 50, 100, 100]' --batch-size '[512, 512, 512]' --pad
```

Values are the Python API's (docs/conventions.md, PLAN 3.2): boxes `(z0, y0, x0, z1, y1, x1)`,
0-based and end-exclusive; sizes and per-axis vectors `(z, y, x)`. The call above is
`XR_crop_dataset(..., [1 1 1 100 100 50], 'BatchSize', [512 512 512], 'pad', true)` in MATLAB.

Exit status: 0; 1 when the function raised (its traceback is printed); 2 for a usage or
argument error. A result is printed as JSON (or its repr).

## `pk5d-master`

```
pk5d-master <MATLAB function name> <argument strings ...>
```

The argument grammar of the compiled master (`mcc/parsers/mccMaster.m` and the generated
`mcc/parsers/*_parser.m`):

- **The function name** is a MATLAB name, exact case (mccMaster's `switch`). Any registered
  function is accepted, not only those in mccMaster's list; one that isn't ported gives
  MATLAB's error (`The parser for function X does not exist!`) and the registry's reason.
- **Positional strings**: the first n strings go to the function's n required and optional
  arguments, in order. All n are needed (mccMaster indexes `varargin{1..n}`: too few strings
  give `MATLAB:badsubscript`), an optional one too; give `""` or `[]` for its default.
- **Name-value pairs** follow, with names in any case and unique prefixes, as inputParser
  matches them (`bachSize` is `MATLAB:InputParser:UnmatchedParameter`, `save` with more than
  one `save...` parameter `MATLAB:InputParser:AmbiguousParameter`, a name without a value
  `MATLAB:InputParser:ParamMissingValue`). The last of two values wins.
- **Each string is converted by its parameter's rule** in the parser files (from the
  function's validator; `Param(mcc=...)` overrides): `str2num` for numbers and logicals
  (`[1,2,3]`, `[1 2 3]`, `1e-3`, `true`, `[true,false]`, `Inf`, `NaN`), `eval` when the text
  starts with `{` (cells: `{'a','b'}`, `{'a';'b'}`), `[` or `@` where the validator allows
  those, else the text itself (`Cropped`, a path). Then the function's own validator checks
  the value (`pad 1` and `pad TRUE` are `MATLAB:InputParser:ArgumentFailedValidation`: 1 is
  a double and `TRUE` is no MATLAB name, so `str2num` gives `[]`).
- **MATLAB's values**: boxes `[ymin xmin zmin ymax xmax zmax]` 1-based, sizes `[y x z]`, as in
  MATLAB; the registry's conversion kinds give the function the port's values.
- **Quirks kept from MATLAB**: a pair that names a positional argument (`inputBbox [...]` after
  the box) passes the parser and is ignored, since the parser hands the function its
  positional variable; a function whose parser takes no name-value pairs
  (`stitch_global_grid_assignment_wrapper`, ...) ignores strings beyond its arguments.
- **Function handles** in strings (`@(x) ...`) are refused (PLAN open question 6).

Output, as the compiled master prints it:

```
Hostname: node12 
Function name: XR_crop_dataset
Variables: {'/mnt/projects/embryo_2024'} [1,1,1,100,100,50] batchSize [512,512,512] 
...the function's own output...
Elapsed time is 12.345678 seconds.

```

(`Hostname` and `Variables` end with a blank, as in MATLAB.) The exit status is 0, or **249**,
the status of the Linux compiled master (MATLAB Runtime R2024b) after an uncaught error. The
error goes to stderr as MATLAB prints it:

```
Error using XR_crop_dataset
'bachSize' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.

MATLAB:InputParser:UnmatchedParameter
```

(MATLAB adds `(line N)` and `Error in ...` stack lines.) An error of the Python code prints
its traceback first; the last line is its MATLAB error identifier when it has one.

`PETAKIT5D_MASTER_DRY_RUN=1` binds the strings and prints the Python call instead of making
it. `pk5d parse <function> <strings>` prints the bound values in both forms as JSON.

Start-up: only the registry is imported before the function's module, so a task starts in
about 0.1 s (process start to the function's first line, environment on network storage;
0.08 s on local disk), where the compiled master takes 3.4 to 3.8 s.

## The run_mccMaster.sh shim: GUI and PyPetaKit5D

The GUI (PetaKit5D-GUI) and PyPetaKit5D start
`<install>/PetaKit5D/mcc/linux/run_mccMaster.sh <MATLAB Runtime folder> <function> <strings>`
(mac: `mac/run_mccMaster.sh`). The shim has that call convention and runs pk5d-master:

```
pk5d mcc-shim /opt/PetaKit5D-GUI/PetaKit5D/mcc/linux/run_mccMaster.sh
```

writes a launcher that

- ignores the Runtime folder (`run_mccMaster.sh` with no arguments prints MATLAB's usage);
- treats the other arguments as MATLAB's launcher does: each in double quotes, the line
  evaluated by `sh`. This undoes PyPetaKit5D's escapes of `\ " $ ``` and turns the GUI's
  `""` (its quoting of an empty value) into an empty argument, so both callers' strings
  arrive as the compiled master got them (and a GUI path with `$` is expanded, as before);
- runs `<this environment's python> -m petakit5d.cli.master`; `--command` or the environment
  variable `PETAKIT5D_MASTER` give another command.

Replacing the file in a GUI install (keep the original) makes the GUI run the port; nothing
in the GUI changes. The GUI's `linux_with_jvm/run_mccMaster.sh` is used for
XR_visualize_OTF_mask_segmentation, XR_FSC_analysis_wrapper and the PSF analysis, which are not
ported. On Windows the GUI starts `mcc/windows/mccMaster` (an .exe) without a Runtime folder:
a copy of the `pk5d-master.exe` console script under that name has the same convention
(untested). A one-line GUI change to read the master's path from its settings would make the
copy unnecessary (PLAN open question 5).

## What the corpus test covers

`tests/cli/test_corpus.py` runs, for every Phase 1 function the GUI and PyPetaKit5D start,
the calls they make (built as they build them: `tests/cli/callers.py`, from the GUI's
`mainwindow.cpp` submit functions and PyPetaKit5D's generated wrappers in
`tests/cli/data/pypetakit5d_wrappers.json`) plus the grammar's edge cases, through MATLAB's
mccMaster and parser in the oracle tree (the task function replaced by a recorder, its
inputParser applied) and through pk5d-master, and compares the values the function gets
(in MATLAB's form), the errors and the `Function name / Variables` lines.

Known differences (strict xfail cases): a hex literal (`0x100`, a uint16 in MATLAB) is not read
by the port's literal parser; a NaN box end is refused by the `bbox` kind where MATLAB passes
it on. A GUI call that MATLAB rejects too: the pipeline tab's bounding-box crop passes
`boundboxCrop`, which `XR_microscopeAutomaticProcessing` doesn't have.
