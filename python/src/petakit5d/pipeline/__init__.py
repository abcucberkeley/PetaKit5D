"""The microscope pipeline: dataset file names, setting files and automatic processing.

==========================================  ===================================
``parse_image_filenames``                   parseImageFilenames
``parse_dataset_filenames``                 XR_parseImageFilenames
``group_image_filenames_by_channels``       groupImageFilenamesByChannels
``matches_channel_pattern``                 channel pattern test of the wrappers
``matlab_regex``                            MATLAB regular expressions
``parse_filename_fields``                   tile name grammar
``parse_setting_files``                     XR_parseSettingFiles
``parse_image_setting_files``               XR_parseSettingFiles_wrapper
``automatic_processing``                    XR_microscopeAutomaticProcessing
==========================================  ===================================

Nothing is imported until a function is first used.
"""

_EXPORTS = {
    "DatasetFilenames": "filenames",
    "ImageFilenames": "filenames",
    "MatlabError": "filenames",
    "group_image_filenames_by_channels": "filenames",
    "matches_channel_pattern": "filenames",
    "matlab_regex": "filenames",
    "parse_dataset_filenames": "filenames",
    "parse_filename_fields": "filenames",
    "parse_image_filenames": "filenames",
    "parse_image_setting_files": "settings",
    "parse_setting_files": "settings",
    "automatic_processing": "automatic",  # P02
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    try:
        mod = importlib.import_module(f"{__name__}.{module}")
    except ModuleNotFoundError as err:
        if err.name != f"{__name__}.{module}":
            raise
        raise AttributeError(f"{name} is not ported yet ({__name__}.{module})") from None
    value = getattr(mod, name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_EXPORTS))
