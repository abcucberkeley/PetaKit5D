"""Deskew and rotation of light-sheet volumes, files and datasets.

========================================  ========================================================
``deskew_rotate_volume``                  deskewRotateFrame3D
``deskew_volume``                         deskewFrame3D
``rotate_volume``                         rotateFrame3D
``rotate_file``                           XR_RotateFrame3D
``deskew_rotate_file``                    XR_deskewRotateFrame
``deskew_rotate_dataset``                 XR_deskew_rotate_data_wrapper
``deskew_rotate_zarr``                    XR_deskewRotateZarr
``deskew_rotate_block``                   XR_deskewRotateBlock
``deskew_rotate_boxes_from_masks``        XR_getDeskeRotateBoxesFromMasks
``scale_contrast``                        scaleContrast (the rescale of ``rescale_rotate``)
``estimate_deskew_rotate_memory``         XR_estimateDeskewRotateMemory (measured for the port)
========================================  ========================================================

``deskew_rotate_volume_matlab_layout`` is ``deskew_rotate_volume`` on MATLAB's layout, for
file-level code; ``skewed_space_interp_matlab`` the MATLAB function
``skewed_space_interp_defined_stepsize.m``. The geometry (output sizes, the skewed space
interpolation step, the transforms) is in ``deskew_rotate.geometry``, a large file's batches in
``large_file_geometry``, the resources of its job in ``controller_resources``. Nothing is
imported until a name is first used.
"""

_EXPORTS = {
    "deskew_rotate_volume": "volume",
    "deskew_rotate_volume_matlab_layout": "volume",
    "deskew_volume": "volume",
    "rotate_volume": "volume",
    "rotate_file": "volume",
    "skewed_space_interp_matlab": "volume",
    "deskew_rotate_file": "file",
    "scale_contrast": "file",
    "estimate_deskew_rotate_memory": "memory",
    "DeskewRotateDatasetError": "dataset",
    "controller_resources": "dataset",
    "deskew_rotate_dataset": "dataset",
    "DeskewRotateZarrError": "zarr",
    "deskew_rotate_block": "zarr",
    "deskew_rotate_boxes_from_masks": "zarr",
    "deskew_rotate_zarr": "zarr",
    "large_file_geometry": "zarr",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    value = getattr(importlib.import_module(f"{__name__}.{module}"), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_EXPORTS))
