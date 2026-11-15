"""Flat-field correction of light-sheet frames.

========================================  ===================================================
``flatfield_correct_volume``              XR_LSFlatFieldCorrection
``flatfield_correct_file``                processFFCorrectionFrame
``prepare_flat_field``                    (XR_LSFlatFieldCorrection's flat field)
``flatfield_background_dataset``          XR_flatfield_background_wrapper
``flatfield_background_estimation``       XR_flatfield_background_estimation
``estimate_background``                   (XR_flatfield_background_estimation's estimate)
``phase_image_flatfield_correction``      XR_phase_image_flat_field_correction
``phase_prefilter``                       (its Gaussian pre-filter)
``basic_flatfield``                       (its flat field, by BaSiCPy in place of BaSiC)
========================================  ===================================================

Nothing is imported until a function is first used.
"""

_EXPORTS = {
    "flatfield_correct_volume": "correction",
    "flatfield_correct_file": "correction",
    "prepare_flat_field": "correction",
    "flatfield_background_dataset": "estimation",
    "flatfield_background_estimation": "estimation",
    "estimate_background": "estimation",
    "phase_image_flatfield_correction": "basic",
    "phase_prefilter": "basic",
    "basic_flatfield": "basic",
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
