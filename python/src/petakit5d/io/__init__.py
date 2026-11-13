"""Tiff and zarr files (over cpp-tiff and cpp-zarr), image sizes and types.

Arrays are ``(z, y, x)`` in C order (images ``(y, x)``); boxes are ``(z0, y0, x0, z1, y1, x1)``
and ranges ``(start, stop)``, 0-based, end-exclusive; sizes and block sizes are ``(z, y, x)``.
Tiff files are the mex's; zarr arrays written by the port are native (``[nz, ny, nx]``, order
C), and MATLAB's (``[ny, nx, nz]``, order F) are still read (docs/decisions/0004-io-bindings.md).

========================  =========================
``read_tiff``             readtiff
``read_tiff_parallel``    readtiff_parallel
``write_tiff``            writetiff
``read_zarr``             readzarr
``write_zarr``            writezarr
``create_zarr``           createzarr
``write_zarr_block``      writeZarrBlock
``get_image_size``        getImageSize
``get_image_size_batch``  getImageSizeBatch
``get_image_dtype``       getImageDataType
``get_zarr_info``         getZarrInfo
``get_image_bounding_box``  getImageBoundingBox
``ReadBuffer``            (one buffer reused by the reads of a loop: ``out=``)
``read_mat``              load (.mat files)
``write_mat``             save (.mat files)
``save_ascii``            save -ascii
========================  =========================

Nothing is imported until a function is first used.
"""

_EXPORTS = {
    "read_tiff": "tiff",
    "read_tiff_parallel": "tiff",
    "write_tiff": "tiff",
    "read_zarr": "zarr",
    "write_zarr": "zarr",
    "create_zarr": "zarr",
    "write_zarr_block": "zarr",
    "get_image_size": "imagesize",
    "get_image_size_batch": "imagesize",
    "get_image_dtype": "imagesize",
    "get_zarr_info": "imagesize",
    "get_image_bounding_box": "imagesize",
    "ZarrInfo": "imagesize",
    "ReadBuffer": "imagesize",
    "read_mat": "matfile",
    "write_mat": "matfile",
    "save_ascii": "matfile",
    "mat_version": "matfile",
    "StructArray": "matfile",
    "MatlabObject": "matfile",
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
