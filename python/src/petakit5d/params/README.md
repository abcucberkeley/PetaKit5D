# The parameter registry: how to declare a user-facing function

Every user-facing function of the port is declared once, next to its code in its own area
folder. The declaration gives the snake_case API, `petakit5d.matlab.<MATLAB name>`, the
command line (`pk5d-master`), the saved parameters struct and the mapping table
(`docs/MATLAB_MAPPING.md`); `tools/check_registry.py` keeps it in step with MATLAB.

## The format

```python
from petakit5d.params import Param, register, matlab_results


@register(
    "XR_crop_dataset",  # MATLAB name, a plain string literal
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param("inputBbox", "@isnumeric", "Box [ymin, xmin, zmin, ymax, xmax, zmax]."),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("pad", "@islogical", "Pad where the box leaves the image."),
    Param("channelPatterns", "@iscell", "Channel identifiers."),
    Param("batchSize", "@isnumeric", "Batch size per task."),
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
- **The signature gives kinds and defaults.** Required arguments: positional without a
  default. MATLAB's `addOptional`: positional with a default. Name-value parameters:
  keyword-only (after `*`). Python names are `snake_name(MATLAB name)`
  (`xyPixelSize` -> `xy_pixel_size`, `DSRDirName` -> `dsr_dir_name`, `save16bit` ->
  `save_16bit`); a different name needs `py="..."` (Python keywords: `lambda` -> `lambda_`).
- **Defaults equal MATLAB's**, written as plain Python: `''` -> `""`, `true` -> `True`,
  `[]` -> `None`, `{}` -> `()`, `[1024, 1024, 1024]` -> `(1024, 1024, 1024)`,
  `{'a', 'b'}` -> `("a", "b")`. Never a list or dict (mutable). A MATLAB default that is
  not a literal (`pwd`, `get_uuid()`) needs `matlab_default="<MATLAB text>"` and the Python
  default that matches it. A site-specific one (an absolute path, Slurm options) never goes
  into the repository: declare `differs="<reason>"` with a generic default (`SlurmParam`:
  `""`, site values come from the cluster config; MATLAB-only keys such as `MCCMasterStr`
  and `MCRParam`: accepted and ignored). The stored extract keeps such MATLAB texts only as
  `<site value sha256:...>`, compared by hash, the literal corpus leaves them out, and the
  mapping table shows `<site-specific>`. The rules are generic (`check.site_pattern`); a
  local scan can add its own expression through `PETAKIT5D_SITE_PATTERN`, set outside the
  repository.
- **The function receives canonical values**: lists become tuples, MATLAB values from the
  MATLAB layer or the command line become str, bool, int (integer-valued doubles), float,
  tuples and None; NumPy arrays pass through untouched. Values are checked against the
  validator's class tests (`ischar`, `islogical`, `isnumeric`, `iscell`, ... ; other terms
  such as `numel(x) == 3` are not checked), as inputParser checks them: `pad=1` is an error.
- Options of `register`: `doc=` (one line for the table; default: the docstring's first
  paragraph), `keep_unmatched=True` (MATLAB's `ip.KeepUnmatched = true`), `task=`.
- `Param` options: `py=`, `matlab_default=`, `mcc=` (the command-line conversion when the
  MATLAB parser file differs from the generator's rule; check_registry tells you),
  `python_only=True` (a parameter MATLAB does not have; needs the orchestrator's OK),
  `differs=`.
- A MATLAB function that is not ported: `not_ported("XR_name", "reason", "area")` in your
  area's module.

## Getting a skeleton

```
cd python
PYTHONPATH=src python tools/check_registry.py skeleton XR_crop_dataset --matlab-root <MATLAB tree>
```

prints the decorator and signature with MATLAB's names, validators, defaults and header
meanings, and the planned function name (params/planned.py). Edit the meanings, paste it.

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
from petakit5d.tools import crop_dataset  # snake_case

crop_dataset("/mnt/projects/embryo_2024", (1, 1, 1, 100, 100, 50), batch_size=(512, 512, 512))

from petakit5d import matlab  # MATLAB names, any case, unique prefixes

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
identifiers on the MATLAB and command-line paths.
