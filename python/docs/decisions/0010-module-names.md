# 0010: A function's module names its MATLAB folder or kind

Status: accepted by the user, 2026-11-18 ("Yes, and note that"), on the mapping table review of
task E05. Three modules renamed before the first release.

## Context

The port's names follow one rule (PLAN.md 3.3; the registry's planned names): the `XR_` prefix
is dropped, data wrappers are `*_dataset`, `Frame` functions `*_file`, in-memory functions
`*_volume`, single large zarr files `*_zarr`, their block tasks `*_block`, and a function sits in
the module of its area and kind. The review of the mapping table found three modules that broke
it: `XR_RotateFrame3D`'s `rotate_file` in `deskew_rotate.volume` (the in-memory functions),
decon's two file-level functions in each other's modules (`XR_RLdeconFrame3D`'s `decon_file` in
`decon.frame`, `RLdecon`'s `rl_decon_file` in `decon.file`), and `stitch.imagelist_generate`,
whose name reads backwards. Its new name keeps the stem of `stitch.imagelist`, the parsing
of the lists, beside it: `stitch.image_list`, first chosen, would have sat an underscore away
from it.

## Decision

| MATLAB | Was | Now |
|---|---|---|
| `XR_RotateFrame3D` | `petakit5d.deskew_rotate.volume.rotate_file` | `petakit5d.deskew_rotate.file.rotate_file` |
| `XR_RLdeconFrame3D` | `petakit5d.decon.frame.decon_file` | `petakit5d.decon.file.decon_file` |
| `RLdecon` | `petakit5d.decon.file.rl_decon_file` | `petakit5d.decon.rl_file.rl_decon_file` |
| `XR_generate_image_list_wrapper` and its generators | `petakit5d.stitch.imagelist_generate` | `petakit5d.stitch.imagelist_generation` |

- The functions, their parameters, the names the packages export and the MATLAB-name layer
  stay as they were; the test modules follow their modules (`tests/decon/test_file.py`,
  `test_rl_file.py`, `tests/stitch/test_imagelist_generation.py`).
- **No compatibility shims:** the port had not been released, so nothing outside the repository
  imports the old paths. The migration guide lists the renames
  ([migration.md](../migration.md), "Renamed before the first release").
