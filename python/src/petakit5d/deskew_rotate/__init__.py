"""Deskew and rotation of light-sheet volumes and files.

========================================  =====================================
``deskew_rotate_volume``                  deskewRotateFrame3D
``deskew_volume``                         deskewFrame3D
``rotate_volume``                         rotateFrame3D
``rotate_file``                           XR_RotateFrame3D
========================================  =====================================

``deskew_rotate_volume_matlab_layout`` is ``deskew_rotate_volume`` on MATLAB's layout, for
file-level code; ``skewed_space_interp_matlab`` the MATLAB function
``skewed_space_interp_defined_stepsize.m``. The geometry (output sizes, the skewed space
interpolation step, the transforms) is in ``deskew_rotate.geometry``. Nothing is imported until a
function is first used.
"""

_EXPORTS = {
    "deskew_rotate_volume": "volume",
    "deskew_rotate_volume_matlab_layout": "volume",
    "deskew_volume": "volume",
    "rotate_volume": "volume",
    "rotate_file": "volume",
    "skewed_space_interp_matlab": "volume",
}
# the planned file, dataset and large-file functions (PLAN.md 3.3), in the modules their tasks
# add: found here once those modules exist
_PLANNED = {
    "deskew_rotate_file": "file",
    "deskew_rotate_dataset": "dataset",
    "deskew_rotate_zarr": "zarr",
    "deskew_rotate_block": "zarr",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    module = _EXPORTS.get(name) or _PLANNED.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    try:
        mod = importlib.import_module(f"{__name__}.{module}")
    except ModuleNotFoundError as err:
        if err.name != f"{__name__}.{module}":
            raise
        raise AttributeError(f"module {__name__!r} has no attribute {name!r} (yet)") from None
    value = getattr(mod, name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_EXPORTS))
