"""Flat-field correction of light-sheet frames.

========================================  =====================================
``flatfield_correct_volume``              XR_LSFlatFieldCorrection
``flatfield_correct_file``                processFFCorrectionFrame
``prepare_flat_field``                    (XR_LSFlatFieldCorrection's flat field)
========================================  =====================================

Nothing is imported until a function is first used.
"""

_EXPORTS = {
    "flatfield_correct_volume": "correction",
    "flatfield_correct_file": "correction",
    "prepare_flat_field": "correction",
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
