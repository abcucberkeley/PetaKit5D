"""Toy functions in a registry of their own, for the command-line tests: they record their calls."""

from __future__ import annotations

from typing import Any

from petakit5d.params import Param, Registry

REG = Registry()
CALLS: list[tuple[str, dict[str, Any]]] = []


class ToyError(RuntimeError):
    identifier = "pkpy:toy:failed"


@REG.register(
    "XR_toy_dataset",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1).", convert="bbox"),
    Param("resultDirName", "@ischar", "Result folder."),
    Param("save16bit", "@islogical", "Save 16 bit."),
    Param("saveZarr", "@islogical", "Save zarr."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x).", convert="size_yxz"),
    Param("channelPatterns", "@iscell", "Channel identifiers."),
    Param("fail", "@islogical", "Raise in the function."),
    Param("inputType", "@ischar", "Input type (a name that shares a prefix with inputBbox)."),
    task="J12",
)
def toy_dataset(
    data_paths,
    input_bbox,
    *,
    result_dir_name="Toy",
    save_16bit=True,
    save_zarr=False,
    batch_size=(1024, 1024, 1024),
    channel_patterns=("CamA_ch0", "CamB_ch0"),
    fail=False,
    input_type="tiff",
):
    """Record a call of the toy dataset function."""
    CALLS.append(("toy_dataset", dict(locals())))
    if fail:
        raise ToyError("the toy failed")
    return {"batch": batch_size}


@REG.register(
    "toy_positional_task",
    Param("inputFullpath", "@ischar", "Input."),
    Param("threshold", "@isscalar", "Threshold."),
    Param("frame", "@isnumeric", "Optional frame."),
    task="J12",
)
def toy_positional(input_fullpath, threshold, frame=None):
    """Record a call of a function without name-value pairs."""
    CALLS.append(("toy_positional", dict(locals())))


def _percent_to_fraction(kwargs: dict[str, Any]) -> dict[str, Any]:
    if "level" in kwargs:
        kwargs = {**kwargs, "level": kwargs["level"] / 100}
    return kwargs


@REG.register(
    "toy_hooked_task",
    Param("level", "@isnumeric", "A level in percent (MATLAB) or as a fraction (Python)."),
    task="J12",
    matlab_in=_percent_to_fraction,
)
def toy_hooked(level=0.5):
    """Record a call of a function with a function-level conversion."""
    CALLS.append(("toy_hooked", dict(locals())))
