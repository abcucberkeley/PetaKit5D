"""Deconvolution: PSFs, OTFs, the OMW back projector and Richardson-Lucy on the CPU.

==============================  ========================================================
``psf_gen``                     psf_gen_new (``generate_psf`` with backend and threads)
``crop_psf_support``            RLdecon's crop of a generated PSF to its support
``crop_psf_to_data``            RLdecon's crop of a PSF larger than the data
``read_tiff_cached``            decon_readtiff_cached (``clear_read_cache``)
``rotate_psf``                  XR_rotate_PSF
``psf2otf``                     decon_psf2otf
``omw_backprojector``           omw_backprojector_generation (``generate_backprojector``)
``omw_backprojector_files``     omw_backprojector_files (``backprojector_file_names``)
``decon_lucy``                  decon_lucy_function (``lucy`` with backend and threads)
``decon_lucy_omw``              decon_lucy_omw_function (``lucy_omw``)
``clear_otf_cache``             ``clear decon_lucy_function decon_lucy_omw_function``
``decon_mask_edge_erosion``     decon_mask_edge_erosion
``rl_decon_file``               RLdecon (``RLDeconResult``)
``decon_file``                  XR_RLdeconFrame3D
``decon_large_in_memory``       RLdecon_large_in_memory
``decon_large_in_place``        RLdecon_large_in_place
``decon_zarr_block``            RLdecon_for_zarr_block
==============================  ========================================================

Nothing is imported until a name is first used.
"""

_EXPORTS = {
    "PSFError": "psf",
    "clear_read_cache": "psf",
    "crop_psf_support": "psf",
    "crop_psf_to_data": "psf",
    "generate_psf": "psf",
    "psf_gen": "psf",
    "read_tiff_cached": "psf",
    "rotate_psf": "psf",
    "pad_shift_psf": "otf",
    "psf2otf": "otf",
    "BackProjector": "backprojector",
    "BackProjectorError": "backprojector",
    "backprojector_file_names": "backprojector",
    "generate_backprojector": "backprojector",
    "omw_backprojector": "backprojector",
    "omw_backprojector_files": "backprojector",
    "DeconError": "rl",
    "OMWResult": "rl",
    "RLResult": "rl",
    "clear_otf_cache": "rl",
    "decon_lucy": "rl",
    "decon_lucy_omw": "rl",
    "lucy": "rl",
    "lucy_omw": "rl",
    "decon_mask_edge_erosion": "erosion",
    "RLDeconResult": "file",
    "rl_decon_file": "file",
    "decon_file": "frame",
    "DeconLargeError": "large",
    "decon_large_in_memory": "large",
    "decon_large_in_place": "large",
    "in_place_controller_resources": "large",
    "decon_zarr_block": "block",
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
