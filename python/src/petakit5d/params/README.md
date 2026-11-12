# The parameter registry: how to declare a user-facing function

Every user-facing function of the port is declared once, next to its code in its own area
folder. The declaration gives the snake_case API, `petakit5d.matlab.<MATLAB name>`, the
command line (`pk5d-master`), the saved parameters struct and the mapping table
(`docs/MATLAB_MAPPING.md`); `tools/check_registry.py` keeps it in step with MATLAB.

**Values follow the layer** (PLAN.md 3.2, 3.3; CONVENTION_BRIEF.md): the snake_case function
takes and returns Python values (boxes `(z0, y0, x0, z1, y1, x1)` 0-based and end-exclusive,
sizes and per-axis vectors `(z, y, x)`, NumPy axes, `'zyx'` axis strings, `(z, y, x)` C
arrays). The MATLAB-name layer, the command line, configs and parameters.json/.mat hold
MATLAB's values. Each parameter's conversion kind (`convert=`) turns one into the other at
that boundary, so the code of the port only ever sees the Python forms.

## The format

```python
from petakit5d.params import Param, register, matlab_results


@register(
    "XR_crop_dataset",  # MATLAB name, a plain string literal
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1).", convert="bbox"),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("pad", "@islogical", "Pad where the box leaves the image."),
    Param("channelPatterns", "@iscell", "Channel identifiers."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) per task.", convert="size_yxz"),
    task="T03",
)
def crop_dataset(
    data_paths,
    input_bbox,
    *,
    result_dir_name="Cropped",
    pad=False,
    channel_patterns=("CamA_ch0", "CamB_ch0"),
    batch_size=(1024, 1024, 1024),
):
    """Crop every file of the datasets with one bounding box."""
    pr = matlab_results(crop_dataset, locals())  # MATLAB's ip.Results, for parameters.mat/.json
    ...
```

Rules:

- **One `Param` per parameter, in signature order**, with the MATLAB name and the MATLAB
  validator text copied verbatim from the `ip.addRequired/addOptional/addParameter` line
  (for an `arguments` block: the declaration, e.g. `"(1, :) {mustBeNumeric}"`; for a
  function without either: `None`). The third field is the meaning (from the MATLAB header).
- **The signature gives argument kinds and defaults.** Required arguments: positional without
  a default. MATLAB's `addOptional`: positional with a default. Name-value parameters:
  keyword-only (after `*`). Python names are `snake_name(MATLAB name)`
  (`xyPixelSize` -> `xy_pixel_size`, `DSRDirName` -> `dsr_dir_name`, `save16bit` ->
  `save_16bit`); a different name needs `py="..."` (Python keywords: `lambda` -> `lambda_`).
- **A conversion kind for every value with axes or 1-based numbers** (`convert=`):

  | kind | MATLAB | Python |
  |---|---|---|
  | `bbox` | `[ymin xmin zmin ymax xmax zmax]` (2D `[ymin xmin ymax xmax]`) | `(z0, y0, x0, z1, y1, x1)` (2D `(y0, x0, y1, x1)`); Inf ends: None |
  | `bboxes` | rows of boxes (`n x 6`) | tuple of boxes |
  | `size_yxz` | `[y x z]`; 6 or 9 values in groups of three (pool sizes); 2D `[y x]` | `(z, y, x)` per group; 2D `(y, x)` |
  | `per_axis_xyz` | `[x y z]` | `(z, y, x)` |
  | `resample_yxz`, `resample_xyz` | `[y x z]` / `[x y z]`, or `[xy z]` | `(z, y, x)`, or `(z, xy)` |
  | `flags_yxz` | `[y x z]` flags (a z MIP is `[0 0 1]`) | `(z, y, x)` (`(1, 0, 0)`) |
  | `index1`, `point1_yxz` | 1-based indices, 1-based `[y x z]` points | 0-based, `(z, y, x)` |
  | `range1` | `[first last]` | `(start, stop)`; Inf: None |
  | `dim` | dimension numbers (1 y, 2 x, 3 z) | NumPy axes (-2, -1, -3) |
  | `axis_string`, `data_order` | `'yxz'`, `'y,x,z'` | `'zyx'`, `'z,y,x'` |
  | `volume` | `(y, x, z)` Fortran image array | `(z, y, x)` C array |
  | `same` | checked: no conversion (stage coordinates in x, y, z, counts, intensities) | |

  Decide each kind from the MATLAB code that uses the value, not from its name: most
  per-axis vectors are `[y x z]`, some are `[x y z]` (the resample factors of the deskew
  functions and of stitching). `params/planned.py` holds the kinds decided so far for the
  Phase 1 functions; the skeleton uses them. Anything else is a pair of functions
  `convert=(from_matlab, to_matlab)`. check_registry fails for a parameter whose name looks
  like a box, a size, an axis or an array and has no kind.
- **Function-level conversions**: `register(..., matlab_in=f, matlab_out=g)` (or module-level
  functions named `matlab_in` and `matlab_out`, as tools/memory_estimate.py declares them):
  `f(kwargs) -> kwargs` after the kinds, `g(result) -> result` for what the MATLAB-name layer
  and the command line return (a size the function returns, say).
- **Defaults equal MATLAB's after conversion**, written as plain Python values in the
  Python form: `''` -> `""`, `true` -> `True`, `[]` -> `None`, `{}` -> `()`, `{'a', 'b'}` ->
  `("a", "b")`, `[500 500 250]` with `size_yxz` -> `(250, 500, 500)`, `'yxz'` with
  `axis_string` -> `'zyx'`. Never a list or dict (mutable). A MATLAB default that is not a
  literal (`pwd`, `get_uuid()`) needs `matlab_default="<MATLAB text>"` and the Python default
  that matches it. A site-specific one (an absolute path, Slurm options) never goes into the
  repository: declare `differs="<reason>"` with a generic default (`SlurmParam`: `""`, site
  values come from the cluster config; MATLAB-only keys such as `MCCMasterStr` and
  `MCRParam`: accepted and ignored). The stored extract keeps such MATLAB texts only as
  `<site value sha256:...>`, compared by hash, the literal corpus leaves them out, and the
  mapping table shows `<site-specific>`. The rules are generic (`check.site_pattern`); a
  local scan can add its own expression through `PETAKIT5D_SITE_PATTERN`, set outside the
  repository.
- **The function receives canonical values**: lists become tuples, MATLAB values from the
  MATLAB layer or the command line become str, bool, int (integer-valued doubles), float,
  tuples and None, then the kinds convert them; NumPy arrays pass through untouched (a
  `volume` is converted). Values are checked against the validator's class tests on
  MATLAB's form (`ischar`, `islogical`, `isnumeric`, `iscell`, ... ; other terms such as
  `numel(x) == 3` are not checked), as inputParser checks them: `pad=1` is an error.
- Options of `register`: `doc=` (one line for the table; default: the docstring's first
  paragraph), `keep_unmatched=True` (MATLAB's `ip.KeepUnmatched = true`), `task=`,
  `matlab_in=`, `matlab_out=`.
- `Param` options: `convert=`, `py=`, `matlab_default=`, `mcc=` (the command-line conversion
  when the MATLAB parser file differs from the generator's rule; check_registry tells you),
  `python_only=True` (a parameter MATLAB does not have; needs the orchestrator's OK),
  `differs=`.
- A MATLAB function that is not ported: `not_ported("XR_name", "reason", "area")` in your
  area's module.

## Getting a skeleton

```
cd python
PYTHONPATH=src python tools/check_registry.py skeleton XR_crop_dataset --matlab-root <MATLAB tree>
```

prints the decorator and signature with MATLAB's names, validators, header meanings, the
planned kinds and the defaults in the Python form, and the planned function name
(params/planned.py). A parameter whose kind is open gets a comment. Edit, paste.

## Checking

```
PYTHONPATH=src python tools/check_registry.py                          # against the stored extract (CI)
PYTHONPATH=src python tools/check_registry.py --matlab-root <MATLAB tree>   # against a tree
PYTHONPATH=src python tools/gen_mapping.py                              # regenerate the table
```

The stored extract (`tests/params/data/inputparser_extract.json`) covers every planned
function, every mcc parser function and the registered ones. A function outside it needs
`tools/check_registry.py extract --matlab-root <oracle tree> --commit <commit>` (F03, or ask
the orchestrator). Regenerate `docs/MATLAB_MAPPING.md` after registering (CI checks it).

## The three ways in

```python
from petakit5d.tools import crop_dataset  # snake_case, Python values

crop_dataset("/mnt/projects/embryo_2024", (0, 0, 0, 50, 100, 100), batch_size=(512, 512, 512))

from petakit5d import matlab  # MATLAB names (any case, unique prefixes) and values

matlab.XR_crop_dataset(
    "/mnt/projects/embryo_2024", [1, 1, 1, 100, 100, 50], "BatchSize", [512, 512, 512]
)

from petakit5d.params import find  # command-line strings, as mccMaster passes them

find("XR_crop_dataset").call_strings(
    ["/mnt/projects/embryo_2024", "[1,1,1,100,100,50]", "batchSize", "[512,512,512]"]
)
```

All three make the same call. Unknown names are errors with a suggestion (`'bachSize' is not
a recognized parameter. ... Did you mean 'batchSize'?`), with inputParser's error
identifiers on the MATLAB and command-line paths. `spec.parameters_json(values)` is the text
MATLAB writes to parameters.json for the same call.
