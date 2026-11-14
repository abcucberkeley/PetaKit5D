"""Stitching (``microscopeDataProcessing/stitch``).

Image lists (``imagelist``, ``imagelist_generate``):

================================  ===============================================
``parse_image_list``              stitch_parse_image_list_information
``parse_multi_loc_image_list``    stitch_parse_multi_loc_image_list_information
``parse_tile_info``               parseTileInfo
``process_filenames``             stitch_process_filenames
``generate_image_list``           XR_generate_image_list_wrapper
``estimate_actual_step_size_from_encoder``  XR_estimate_actual_step_size_from_encoder
``read_table``, ``write_table``   readtable, writetable (image list csv files)
``Table``, ``ImageListInfo``      the tables and the parse results
================================  ===============================================

One volume (``volume``):

==============================================  ====================================
``stitch_volume``                               XR_stitching_frame_zarr_dev_v1
``read_stitch_info``, ``write_stitch_info``     the stitch information files (``stitchInfo/``)
``read_xcorr_result``, ``write_xcorr_result``   the registration result files (``xcorr/``)
``StitchVolumeError``, ``StitchVolumeWarning``  the frame's errors and warnings
==============================================  ====================================

Datasets (``dataset``):

==============================================  ====================================
``stitch_dataset``                              XR_matlab_stitching_wrapper
``StitchDatasetError``                          its errors
==============================================  ====================================

Nothing is imported until a name is first used.
"""

_EXPORTS = {
    "ImageListError": "imagelist",
    "ImageListInfo": "imagelist",
    "MatlabWarning": "imagelist",
    "Table": "imagelist",
    "contains_pattern": "imagelist",
    "make_valid_name": "imagelist",
    "parse_image_list": "imagelist",
    "parse_multi_loc_image_list": "imagelist",
    "parse_tile_info": "imagelist",
    "process_filenames": "imagelist",
    "read_table": "imagelist",
    "sortrows": "imagelist",
    "write_table": "imagelist",
    "estimate_actual_step_size_from_encoder": "imagelist_generate",
    "generate_image_list": "imagelist_generate",
    "generate_image_list_from_encoder": "imagelist_generate",
    "generate_image_list_from_encoder_ism": "imagelist_generate",
    "generate_image_list_from_sqlite": "imagelist_generate",
    "generate_image_list_from_sqlite_old_format": "imagelist_generate",
    "generate_image_list_from_tile_list": "imagelist_generate",
    "generate_image_list_from_tile_positions": "imagelist_generate",
    "StitchVolumeError": "volume",
    "StitchVolumeWarning": "volume",
    "read_stitch_info": "volume",
    "read_xcorr_result": "volume",
    "stitch_volume": "volume",
    "write_stitch_info": "volume",
    "write_xcorr_result": "volume",
    "StitchDatasetError": "dataset",
    "stitch_dataset": "dataset",
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
