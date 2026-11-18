"""Tools: maximum intensity projections, cropping, memory estimates, rsync, partial volumes,
channel unmixing, resampling, MIP masks, figures and FFT power spectra.

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
``crop_dataset``                          XR_crop_dataset
``crop_file``                             XR_crop_frame
``crop_zarr``                             XR_crop_zarr
``crop_block``                            XR_crop_block
``trim_border``                           trimBorder
``estimate_computing_memory``             XR_estimateComputingMemory
``parallel_rsync``                        XR_parallel_rsync_wrapper
``group_partial_volume_files``            groupPartialVolumeFiles
``unmix_channels_dataset``                XR_unmix_channels_data_wrapper
``unmix_channels_file``                   XR_unmix_channels_frame
``unmix_channels_zarr``                   XR_unmix_channels_zarr
``unmix_channels_block``                  unmix_channels_block
``unmix_channels_gaussian_block``         unmix_channels_gaussian_block
``unmix_channels_sum``                    unmix_channels_sum_mex (the weighted sum)
``resample_dataset``                      XR_resample_dataset
``resample_file``                         XR_resampleFrame
``resample_zarr``                         XR_resampleSingleZarr
``resample_zarr_block``                   resampleZarrBlock
``check_resample_setting``                XR_checkResampleSetting
``generate_mip_masks``                    XR_generate_MIP_masks
``generate_single_mip_mask``              gererate_single_MIP_mask
``visualize_otf_mask_segmentation``       XR_visualize_OTF_mask_segmentation
``visualize_otf_and_mask_outline``        visualize_OTF_and_mask_outline
``visualize_tile_positions``              XR_visualize_tile_positions
``fft_spectrum_dataset``                  XR_fftSpectrumComputingWrapper
``fft_spectrum_file``                     XR_fftSpectrumComputingFrame
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
    "crop_dataset": "crop",
    "crop_file": "crop",
    "crop_zarr": "crop",
    "crop_block": "crop",
    "trim_border": "crop",
    "estimate_computing_memory": "memory_estimate",
    "parallel_rsync": "rsync",
    "group_partial_volume_files": "partial_volumes",
    "UnmixError": "unmix",
    "unmix_channels_block": "unmix",
    "unmix_channels_dataset": "unmix",
    "unmix_channels_file": "unmix",
    "unmix_channels_gaussian_block": "unmix",
    "unmix_channels_sum": "unmix",
    "unmix_channels_zarr": "unmix",
    "ResampleError": "resample",
    "check_resample_setting": "resample",
    "resample_dataset": "resample",
    "resample_file": "resample",
    "resample_zarr": "resample",
    "resample_zarr_block": "resample",
    "MipMaskError": "mip_masks",
    "generate_mip_masks": "mip_masks",
    "generate_single_mip_mask": "mip_masks",
    "FigureError": "figures",
    "visualize_otf_and_mask_outline": "figures",
    "visualize_otf_mask_segmentation": "figures",
    "visualize_tile_positions": "figures",
    "FFTSpectrumError": "fft_spectrum",
    "fft_spectrum_dataset": "fft_spectrum",
    "fft_spectrum_file": "fft_spectrum",
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
