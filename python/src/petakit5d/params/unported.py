"""MATLAB functions that are not ported, with the reason (PLAN.md section 4.3).

These have no area code to sit next to, so they are declared here; an area task declares its
own with ``not_ported(...)`` in its own module.
"""

from petakit5d.params.registry import not_ported

_LEGACY = "legacy code without callers (PLAN.md 4.3)"
_V1 = "called only by the legacy stitcher XR_stitching_frame_v1 (PLAN.md 4.3)"

not_ported("XR_stitching_frame_v1", _LEGACY, "stitch")
not_ported("cross_correlation_stitching_matching", _V1, "stitch")
not_ported("RLdecon_for_ExMPipeline", _V1, "decon")
not_ported("GU_resampleStack3D", _V1, "tools")
not_ported("multires_cross_correlation_registration_imblock", _LEGACY, "stitch")
not_ported("normxcorr3_updated", _LEGACY, "stitch")
not_ported("normxcorr3", _LEGACY + "; third_parties copy", "stitch")
not_ported("writeZarrBlock", _LEGACY, "io")
not_ported("XR_normalize_z_stack", _LEGACY, "stitch")
not_ported("readtiff_parallel", _LEGACY, "io")
not_ported("decon_otf2psf", _LEGACY, "decon")
not_ported("XR_deskewRotateFrame_old", _LEGACY, "deskew_rotate")
not_ported("feather_distance_map_resize_3d_w_space_mex", _LEGACY, "stitch")
not_ported("rescale_zarr_file", "tools/rescale cannot run (PLAN.md 4.3)", "tools")
not_ported("rescale_zarr_block", "tools/rescale cannot run (PLAN.md 4.3)", "tools")
not_ported("XR_cudaDeconFrame3D", "wrapper of an external decon binary (PLAN.md 4.3)", "decon")
not_ported("XR_cppDeconFrame3D", "wrapper of an external decon binary (PLAN.md 4.3)", "decon")
not_ported("XR_java_stitching_wrapper", "Java stitching tools (PLAN.md 4.3)", "stitch")
not_ported("java_stitching_frame_wrapper", "Java stitching tools (PLAN.md 4.3)", "stitch")
_NOSRC = "listed for mcc, but its source is not in the repository (McCParserConsistencyTest)"
not_ported("image_intensity_correction_block", _NOSRC, "tools")
not_ported("replicateDataBlock", _NOSRC, "tools")
not_ported("multiple_region_spectra_encoding_block", _NOSRC + "; stale parser", "tools")
not_ported("mccMaster", "the compiled master; pk5d-master takes its arguments (J12)", "cli")
