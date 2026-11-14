"""The planned Python names of the user-facing MATLAB functions (PLAN.md section 3.3).

Until an area task registers a function, the mapping table shows its planned name from here,
so the names can be reviewed before the ports exist. A registration overrides its entry; a
registered name that differs from the plan is listed in the table (test_mapping fails on it
until this list or the registration changes). Only F03 edits this file; areas never need to.

Naming rule: drop ``XR_``; say what the function acts on (``_dataset`` for the data
wrappers, ``_file`` for the frame functions, ``_volume`` for in-memory arrays, ``_zarr`` for
single large zarr files, ``_block`` for their block tasks); parameters keep MATLAB's names
in snake_case (``spec.snake_name``).
"""

from __future__ import annotations

# (MATLAB name, Python name, task, phase)
PLANNED: tuple[tuple[str, str, str, int], ...] = (
    # job framework
    ("generic_computing_frameworks_wrapper", "petakit5d.jobs.run_tasks", "J02", 1),
    # I/O
    ("readtiff", "petakit5d.io.read_tiff", "IO2", 1),
    ("writetiff", "petakit5d.io.write_tiff", "IO2", 1),
    ("readzarr", "petakit5d.io.read_zarr", "IO2", 1),
    ("writezarr", "petakit5d.io.write_zarr", "IO2", 1),
    ("createzarr", "petakit5d.io.create_zarr", "IO2", 1),
    ("getImageSize", "petakit5d.io.image_size", "IO2", 1),
    ("XR_tiffToZarr_wrapper", "petakit5d.io.tiff_to_zarr_dataset", "T01", 1),
    ("tiffToZarr", "petakit5d.io.tiff_to_zarr", "T01", 1),
    ("XR_zarrToTiff_wrapper", "petakit5d.io.zarr_to_tiff_dataset", "T01", 1),
    ("zarrToTiff", "petakit5d.io.zarr_to_tiff", "T01", 1),
    ("XR_resave_zarr_wrapper", "petakit5d.io.resave_zarr_dataset", "T01", 1),
    ("XR_resaveSingleZarr", "petakit5d.io.resave_zarr", "T01", 1),
    ("resaveZarrBlock", "petakit5d.io.resave_zarr_block", "T01", 1),
    # tools
    ("XR_MIP_wrapper", "petakit5d.tools.mip_dataset", "T02", 1),
    ("XR_MIP_zarr", "petakit5d.tools.mip_zarr", "T02", 1),
    ("saveMIP_zarr", "petakit5d.tools.save_mip_zarr", "T02", 1),
    ("saveMIP_tiff", "petakit5d.tools.save_mip_tiff", "T02", 1),
    ("MIP_block", "petakit5d.tools.mip_block", "T02", 1),
    ("XR_crop_dataset", "petakit5d.tools.crop_dataset", "T03", 1),
    ("XR_crop_frame", "petakit5d.tools.crop_file", "T03", 1),
    ("XR_crop_zarr", "petakit5d.tools.crop_zarr", "T03", 1),
    ("XR_crop_block", "petakit5d.tools.crop_block", "T03", 1),
    ("XR_resample_dataset", "petakit5d.tools.resample_dataset", "T04", 1),
    ("XR_resampleFrame", "petakit5d.tools.resample_file", "T04", 1),
    ("XR_resampleSingleZarr", "petakit5d.tools.resample_zarr", "T04", 1),
    ("resampleZarrBlock", "petakit5d.tools.resample_zarr_block", "T04", 1),
    ("XR_parallel_rsync_wrapper", "petakit5d.tools.parallel_rsync", "T05", 1),
    ("XR_estimateComputingMemory", "petakit5d.tools.estimate_computing_memory", "T05", 1),
    ("groupPartialVolumeFiles", "petakit5d.tools.group_partial_volume_files", "T05", 1),
    # flat field, deskew/rotate
    ("processFFCorrectionFrame", "petakit5d.flatfield.flatfield_correct_file", "D01", 1),
    ("XR_LSFlatFieldCorrection", "petakit5d.flatfield.flatfield_correct_volume", "D01", 1),
    ("deskewRotateFrame3D", "petakit5d.deskew_rotate.deskew_rotate_volume", "D02", 1),
    ("deskewFrame3D", "petakit5d.deskew_rotate.deskew_volume", "D02", 1),
    ("rotateFrame3D", "petakit5d.deskew_rotate.rotate_volume", "D02", 1),
    ("XR_RotateFrame3D", "petakit5d.deskew_rotate.rotate_file", "D02", 1),
    ("XR_deskewRotateFrame", "petakit5d.deskew_rotate.deskew_rotate_file", "D03", 1),
    ("XR_deskew_rotate_data_wrapper", "petakit5d.deskew_rotate.deskew_rotate_dataset", "D04", 1),
    ("XR_deskewRotateZarr", "petakit5d.deskew_rotate.deskew_rotate_zarr", "D04", 1),
    ("XR_deskewRotateBlock", "petakit5d.deskew_rotate.deskew_rotate_block", "D04", 1),
    # pipeline
    ("XR_microscopeAutomaticProcessing", "petakit5d.pipeline.automatic_processing", "P02", 1),
    # deconvolution
    ("psf_gen_new", "petakit5d.decon.psf_gen", "R01", 1),
    ("XR_rotate_PSF", "petakit5d.decon.rotate_psf", "R01", 1),
    ("omw_backprojector_generation", "petakit5d.decon.omw_backprojector", "R02", 1),
    ("decon_lucy_function", "petakit5d.decon.decon_lucy", "R03", 1),
    ("decon_lucy_omw_function", "petakit5d.decon.decon_lucy_omw", "R03", 1),
    ("RLdecon", "petakit5d.decon.rl_decon_file", "R05", 1),
    ("XR_RLdeconFrame3D", "petakit5d.decon.decon_file", "R05", 1),
    ("XR_decon_data_wrapper", "petakit5d.decon.decon_dataset", "R06", 1),
    ("RLdecon_large_in_memory", "petakit5d.decon.decon_large_in_memory", "R07", 1),
    ("RLdecon_large_in_place", "petakit5d.decon.decon_large_in_place", "R07", 1),
    ("RLdecon_for_zarr_block", "petakit5d.decon.decon_zarr_block", "R07", 1),
    # stitching
    ("XR_generate_image_list_wrapper", "petakit5d.stitch.generate_image_list", "S01", 1),
    ("stitch_process_zarr_tile", "petakit5d.stitch.process_zarr_tile", "S02", 1),
    ("erodeVolumeBy2DProjection", "petakit5d.stitch.erode_volume_by_2d_projection", "S02", 1),
    (
        "cross_correlation_registration_wrapper",
        "petakit5d.stitch.cross_correlation_registration",
        "S04",
        1,
    ),
    (
        "stitch_global_grid_assignment_wrapper",
        "petakit5d.stitch.global_grid_assignment",
        "S05",
        1,
    ),
    ("compute_tile_bwdist", "petakit5d.stitch.compute_tile_bwdist", "S06", 1),
    ("compute_tile_bwdist_mip_slabs", "petakit5d.stitch.compute_tile_bwdist_mip_slabs", "S06", 1),
    ("stitch_organize_block_info", "petakit5d.stitch.organize_block_info", "S07", 1),
    ("processStitchBlock", "petakit5d.stitch.process_stitch_block", "S07", 1),
    ("XR_stitching_frame_zarr_dev_v1", "petakit5d.stitch.stitch_volume", "S08", 1),
    ("XR_matlab_stitching_wrapper", "petakit5d.stitch.stitch_dataset", "S09", 1),
    # Phase 2
    ("XR_imaris_conversion_data_wrapper", "petakit5d.tools.imaris_conversion_dataset", "B01", 2),
    ("XR_unmix_channels_data_wrapper", "petakit5d.tools.unmix_channels_dataset", "B02", 2),
    ("XR_unmix_channels_frame", "petakit5d.tools.unmix_channels_file", "B02", 2),
    ("XR_unmix_channels_zarr", "petakit5d.tools.unmix_channels_zarr", "B02", 2),
    ("unmix_channels_block", "petakit5d.tools.unmix_channels_block", "B02", 2),
    ("unmix_channels_gaussian_block", "petakit5d.tools.unmix_channels_gaussian_block", "B02", 2),
    (
        "XR_chromatic_shift_correction_data_wrapper",
        "petakit5d.tools.chromatic_shift_correction_dataset",
        "B03",
        2,
    ),
    ("XR_fftSpectrumComputingWrapper", "petakit5d.tools.fft_spectrum_dataset", "B04", 2),
    ("XR_fftSpectrumComputingFrame", "petakit5d.tools.fft_spectrum_file", "B04", 2),
    ("XR_FSC_analysis_wrapper", "petakit5d.tools.fsc_analysis_dataset", "B05", 2),
    ("XR_one_image_FSC_analysis_frame", "petakit5d.tools.fsc_analysis_file", "B05", 2),
    ("XR_psf_analysis_wrapper", "petakit5d.tools.psf_analysis_dataset", "B06", 2),
    (
        "XR_psf_detection_and_analysis_wrapper",
        "petakit5d.tools.psf_detection_and_analysis_dataset",
        "B06",
        2,
    ),
    ("XR_psf_detection_and_cropping", "petakit5d.tools.psf_detection_and_cropping", "B06", 2),
    ("XR_psf_analysis_plot", "petakit5d.tools.psf_analysis_plot", "B06", 2),
    (
        "XR_flatfield_background_wrapper",
        "petakit5d.flatfield.flatfield_background_dataset",
        "B07",
        2,
    ),
    (
        "XR_phase_image_flat_field_correction",
        "petakit5d.flatfield.phase_image_flatfield_correction",
        "B07",
        2,
    ),
    ("XR_ExM_PunctaRemoval_zarr", "petakit5d.tools.puncta_removal_zarr", "B08", 2),
    (
        "XR_bioformats_to_tiff_or_zarr_wrapper",
        "petakit5d.io.bioformats_to_tiff_or_zarr_dataset",
        "B09",
        2,
    ),
    ("N5ToZarr", "petakit5d.io.n5_to_zarr", "B09", 2),
    ("zarrToN5", "petakit5d.io.zarr_to_n5", "B09", 2),
    ("XR_generate_MIP_masks", "petakit5d.tools.generate_mip_masks", "B10", 2),
    (
        "XR_visualize_OTF_mask_segmentation",
        "petakit5d.tools.visualize_otf_mask_segmentation",
        "B10",
        2,
    ),
    # Phase 3
    ("XR_runDetTrack3d", "petakit5d.detection.run_detection_tracking", "C04", 3),
    ("XR_runDetection3D", "petakit5d.detection.run_detection", "C02", 3),
    ("XR_runTracking3D", "petakit5d.detection.run_tracking", "C03", 3),
    ("XR_runTrackProcessing3D", "petakit5d.detection.run_track_processing", "C04", 3),
    ("XR_loadConditionData3D", "petakit5d.detection.load_condition_data", "C04", 3),
    ("runLifetimeAnalysis3D", "petakit5d.detection.run_lifetime_analysis", "C05", 3),
    ("simReconAutomaticProcessing", "petakit5d.sim.sim_recon_automatic_processing", "C06", 3),
    ("simReconFrame", "petakit5d.sim.sim_recon_file", "C06", 3),
    ("deskewPhasesFrame", "petakit5d.sim.deskew_phases_file", "C06", 3),
    # demos
    (
        "crop_deskew_rotate_demo_function",
        "petakit5d.examples.crop_deskew_rotate_demo_function",
        "E01",
        1,
    ),
)


def planned() -> dict[str, tuple[str, str, int]]:
    """{MATLAB name: (Python name, task, phase)}."""
    return {m: (p, t, ph) for m, p, t, ph in PLANNED}


# -- conversion kinds of the Phase 1 parameters (params/kinds.py), decided from the MATLAB code
# A kind by parameter name; KINDS_BY_FUNCTION overrides it where a function differs. "?<task>"
# marks a parameter whose kind its area task decides (it reads the code that uses it). The
# skeletons and the mapping table use these; a registration's own convert= wins.

KINDS_BY_NAME: dict[str, str] = {
    # boxes [ymin xmin zmin ymax xmax zmax] (1-based, inclusive)
    **dict.fromkeys(
        ["inputBbox", "InputBbox", "tileOutBbox", "outBbox", "bbox", "deconBbox", "inBbox",
         "distBbox"],
        "bbox",
    ),
    # rows of boxes (n x 6)
    **dict.fromkeys(
        ["BatchBBoxes", "RegionBBoxes", "inBatchBBoxes", "outBatchBBoxes", "outRegionBBoxes",
         "outLocalBboxes", "distBboxes"],
        "bboxes",
    ),
    # sizes and per-axis vectors [y x z] (pool sizes: groups of three); the stitching's
    # xcorrDownsample is in the tiles' axis order (the registration max-pools the tiles' arrays
    # by it, [y x z] for 'y,x,z' tiles)
    **dict.fromkeys(
        ["batchSize", "blockSize", "shardSize", "zarrSubSize", "mipBatchSize", "imSize",
         "borderSize", "BorderSize", "stitchPadSize", "poolSize", "dataSize", "xcorrDownsample"],
        "size_yxz",
    ),
    # per-axis flags [y x z] (a z MIP is [0 0 1])
    **dict.fromkeys(["axis", "mipAxis", "GenMaxZproj", "stitchMIP"], "flags_yxz"),
    "inputAxisOrder": "axis_string",
    "outputAxisOrder": "axis_string",
    "dataOrder": "data_order",
    # 1-based indices and coordinates
    **dict.fromkeys(["batchInds", "blockInds", "tileInd", "imdistFileIdx"], "index1"),
    **dict.fromkeys(["startCoord", "lastStartCoords"], "point1_yxz"),
    "range": "range1",
    # image arrays (y, x, z)
    **dict.fromkeys(
        ["vol", "data", "img", "I", "PSF", "PSF_f", "PSF_b", "psf", "Rawdata", "LSImage",
         "rawdata", "frame"],
        "volume",
    ),
    # checked, no conversion: stage and physical values in x, y, z by name (exception 2),
    # counts, intensities, per-channel or per-step flags
    **dict.fromkeys(
        ["axisOrder", "tileInterval", "taskSize", "Overlap", "esize", "TileOffset", "tileOffset",
         "deconOffset", "save3Dstack", "save16bit", "flippedTile", "locIds", "Background",
         "xyz_voxelsizes", "constOffset", "resampleType", "usePrimaryCoords"],
        "same",
    ),
    # decided by the area task
    "downSample": "?S04",
    "MaxOffset": "?S04",
    "cuboid_1": "?S04",
    "cuboid_2": "?S04",
    "cuboid_overlap_12": "?S04",
    "pair_indices": "?S04",
    "data_order_mat": "?S04",
    "dataOrderMat": "?S06",
    "axisWeight": "?S05",
    "halfOrder": "?S03",
    "tileIdx": "?S08",
    "tileIndices": "?S01",
    "overlapSize": "?S01",
    "outSize": "?D02",
}  # fmt: skip

_DESKEW_RS = "resample_xyz"  # passed on to deskewRotateFrame3D's [x y z] resample
KINDS_BY_FUNCTION: dict[str, dict[str, str]] = {
    # resample factors: [y x z] (or [xy z]) where they divide an image size or are read in the
    # tiles' axis order (the stitching frame's xf, yf, zf are resampleFactor(data_order_mat),
    # and its wrapper passes them on), [x y z] where they reach the deskew/rotate transform
    **{
        f: {"resampleFactor": "resample_yxz"}
        for f in ("XR_resample_dataset", "XR_resampleFrame", "XR_resampleSingleZarr",
                  "resampleZarrBlock", "tiffToZarr", "XR_tiffToZarr_wrapper",
                  "XR_stitching_frame_zarr_dev_v1", "XR_matlab_stitching_wrapper")
    },
    "stitch_process_zarr_tile": {"resample": "resample_yxz"},
    **{
        f: {"resampleFactor": _DESKEW_RS}
        for f in ("XR_deskewRotateFrame", "XR_deskew_rotate_data_wrapper", "XR_deskewRotateZarr",
                  "XR_deskewRotateBlock", "XR_microscopeAutomaticProcessing")
    },
    "deskewRotateFrame3D": {"resample": _DESKEW_RS, "resampleFactor": _DESKEW_RS},
    # the memory order flag of a new zarr ('F' in MATLAB; the port writes 'C'): IO3 decides
    "createzarr": {"order": "?IO3"},
    # a crop box or size of the in-memory large-file decon ([] by default): R07 decides
    "RLdecon_large_in_memory": {"Crop": "?R07"},
    "XR_RotateFrame3D": {"resample": _DESKEW_RS},
    "RLdecon": {"Resample": _DESKEW_RS},
    # a size of voxels only in XR_estimateComputingMemory (a count there)
    "XR_estimateComputingMemory": {"dataSize": "same"},
    # an image (the flat field's background) there; a scalar elsewhere
    "XR_LSFlatFieldCorrection": {"background": "volume"},
    "XR_decon_data_wrapper": {"background": "same"},
    # resample of rotateFrame3D: GPL source, not read; D02 decides
    "rotateFrame3D": {"resample": "?D02", "zxRatio": "?D02"},
}  # fmt: skip


def planned_kind(function: str, param: str) -> str | None:
    """The planned conversion kind of a parameter ("?<task>": the area task decides; None:
    no decision needed by name)."""
    own = KINDS_BY_FUNCTION.get(function, {})
    if param in own:
        return own[param]
    return KINDS_BY_NAME.get(param)
