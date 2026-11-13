"""Tools: maximum intensity projections, memory estimates, rsync, partial volumes.

========================================  =====================================
``mip_dataset``                           XR_MIP_wrapper
``mip_zarr``                              XR_MIP_zarr
``mip_block``                             MIP_block
``save_mip_tiff``                         saveMIP_tiff
``save_mip_zarr``                         saveMIP_zarr
``save_mip_frame``                        saveMIP_frame
``project_3d_image_to_2d``                project3DImageto2D
``mip_zarr_batch_size``                   MIP_zarr_batch_size
``block_mips``                            MIP_zarr_block_mips
``save_block_mips``                       MIP_zarr_save_block_mips
``merge_block_mips``                      MIP_zarr_merge_block_mips
``merge_memory``                          MIP_zarr_merge_memory
``estimate_computing_memory``             XR_estimateComputingMemory
``parallel_rsync``                        XR_parallel_rsync_wrapper
``group_partial_volume_files``            groupPartialVolumeFiles
========================================  =====================================

Nothing is imported until a function is first used.
"""

_EXPORTS = {
    "mip_dataset": "mip",
    "mip_zarr": "mip",
    "mip_block": "mip",
    "save_mip_tiff": "mip",
    "save_mip_zarr": "mip",
    "save_mip_frame": "mip",
    "project_3d_image_to_2d": "mip",
    "mip_zarr_batch_size": "mip",
    "block_mips": "mip",
    "save_block_mips": "mip",
    "merge_block_mips": "mip",
    "merge_memory": "mip",
    "axis_max": "mip",
    "estimate_computing_memory": "memory_estimate",
    "parallel_rsync": "rsync",
    "group_partial_volume_files": "partial_volumes",
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
