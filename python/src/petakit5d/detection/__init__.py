"""3D point source detection, tracking and lifetime analysis (Phase 3).

New implementations from the published methods (Aguet et al. 2013 for the model-based detection,
Jaqaman et al. 2008 for the tracking) and from MATLAB's behaviour recorded as black-box
references; no code of the GPL sources is used. ``docs/detection_scope.md`` lists the stages,
their inputs and outputs, and the references.

Modules:

- ``stats``: MATLAB's statistics functions as the detection code calls them (norminv, tcdf, tinv,
  chi2inv, fcdf, ksdensity, ecdf, mad, kstest2).
- ``detect``: the lowSNR detector (``XR_pointSourceDetection3D_3``): the overall background
  estimate, the voxel-wise model fits, the mask, the starts, the fits and the filters after them.
- ``fit``: the Gaussian fits (``XR_fitGaussian3D``'s 'CD_Newton', ``XR_fitGaussian3D_Ac``,
  ``XR_fitGaussians3D``) and ``adtest_mex``'s test; the batched kernel ``fit_cd`` is the API
  for other stages.
- ``run``: ``XR_runDetection3D`` (``run_detection``): every frame of a condition, the other
  channels at the master's detections, ``Detection3D.mat`` and the masks.
- ``condition``: the loading of a condition's movies (``XR_loadConditionData3D``,
  ``XR_loadConditionData``: ``load_condition_data``).
- ``offsets``: the channels' x and z offsets (``XR_correctXZoffsetData3D``,
  ``XZoffsetEstimation3D``).
- ``wrapper``: the detection and tracking wrapper (``XR_runDetTrack3d``:
  ``run_detection_tracking``) and the deskew of a condition's frames (``deskew_data``).
- ``processing``: track processing (``XR_runTrackProcessing3D``: ``run_track_processing``,
  ``ProcessedTracks.mat``) and the tracks' rotation (``rotateTracks3D``: ``rotate_tracks``,
  ``RotatedTracks.mat``).
- ``amira``: the processed tracks as Amira files (``XR_amiraWriteTracks``: ``write_tracks``).
- ``lifetime``: the lifetime analysis of processed tracks (``runLifetimeAnalysis3D``:
  ``run_lifetime_analysis``; ``plotIntensityCohorts``' numbers: ``plot_intensity_cohorts``).
"""

from petakit5d.detection.condition import load_condition_data
from petakit5d.detection.lifetime import plot_intensity_cohorts, run_lifetime_analysis
from petakit5d.detection.offsets import correct_xz_offset_data
from petakit5d.detection.processing import rotate_tracks, run_track_processing
from petakit5d.detection.run import run_detection, run_detection_frame_parallel
from petakit5d.detection.wrapper import deskew_data, run_detection_tracking

__all__ = [
    "correct_xz_offset_data",
    "deskew_data",
    "load_condition_data",
    "plot_intensity_cohorts",
    "rotate_tracks",
    "run_detection",
    "run_detection_frame_parallel",
    "run_detection_tracking",
    "run_lifetime_analysis",
    "run_track_processing",
]
