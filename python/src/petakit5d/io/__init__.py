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
``n5_to_zarr``            N5ToZarr
``n5_to_zarr_block``      N5ToZarrBlock
``zarr_to_n5``            zarrToN5
``zarr_to_n5_block``      zarrToN5Block
``read_n5``               (N5Adapter's getIORegion: MATLAB's N5 layout)
``write_n5``              (N5Adapter's setRegion)
``create_n5``             (blockedImage with N5Adapter, mode 'w')
``construct_xml_file``    constructXmlFile (the BigStitcher XML)
``write_struct_xml``      writestruct (of its struct)
``remove_text_tags``      removeTextTags
``clean_xml_file``        cleanXmlFile
``bioformats_to_tiff_or_zarr_dataset``  XR_bioformats_to_tiff_or_zarr_wrapper (ND2, CZI)
========================  =========================

N5 blosc and zstd (before Python 3.14) need numcodecs, the ``n5`` extra; ND2 and CZI need
nd2, czifile and imagecodecs, the ``formats`` extra (Bio-Formats is not ported).

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
    "n5_to_zarr": "n5",
    "n5_to_zarr_block": "n5",
    "zarr_to_n5": "n5",
    "zarr_to_n5_block": "n5",
    "read_n5": "n5",
    "write_n5": "n5",
    "create_n5": "n5",
    "open_n5": "n5",
    "N5Array": "n5",
    "N5Error": "n5",
    "construct_xml_file": "n5",
    "write_struct_xml": "n5",
    "remove_text_tags": "n5",
    "clean_xml_file": "n5",
    "bioformats_to_tiff_or_zarr_dataset": "bioformats",
    "bioformats_to_tiff_or_zarr_file": "bioformats",
    "BioformatsError": "bioformats",
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
