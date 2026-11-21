"""The point-detection demo step by step (MATLAB: ``demos/demo_point_detection_3d_by_step.m``):
the frames of the dataset's first experiment deskewed, the channels' x and z offsets
corrected, the diffraction-limited spots detected and tracked, the tracks processed and rotated,
one function per step, with the demo's parameters.

    python demo_point_detection_3d_by_step.py

The data: ``PETAKIT5D_DEMO_DATA`` (``_phase3``). Results, as MATLAB writes them: the deskewed
frames in ``ch1/DS/`` and ``ch2/DS/`` with ``xz_offset.mat``, and in the master channel's
``Analysis_update_3/``: ``settings.mat`` and ``dynRange.mat`` (the deskew), ``Detection3D.mat``
and ``Masks/`` (the detection), ``trackedFeatures.mat`` (the tracking), ``ProcessedTracks.mat``
and the Amira files of the tracks (``Amira/Tracking/``), and ``RotatedTracks.mat``.

Step 5 tracks the movies of more than one frame. MATLAB's demo keeps the results folders of all
movies for them (``apath`` is not filtered with ``data``), so after a one-frame movie each movie
would get the folder of the movie before it, the slip of ``XR_runDetTrack3d`` that the port's
wrapper fixes (``docs/parity/detection.md`` W3); the example keeps each movie with its own
folder. The demo dataset has no one-frame movie, so the results are the same.
"""

from __future__ import annotations

import _phase3

# isort: split
import _common

DEMO = "demo_point_detection_3d_by_step"


def deskew(data) -> None:
    """Step 3a: the frames deskewed (``deskewData``) into each channel's ``DS/`` folder, with
    ``settings.mat`` and ``dynRange.mat`` in the results folder"""
    from petakit5d.detection import deskew_data

    deskew_data(
        data,
        overwrite=True,
        skew_angle=31.5,
        aname=_phase3.ANAME,
        rotate=False,
        s_cmos_camera_flip=False,
        crop=False,
    )


def xz_offsets(data) -> None:
    """Step 3b: the channels' offsets in x and z, from the maximum cross-correlation of their
    deskewed frames (``XR_correctXZoffsetData3D``); frames with an offset are rewritten
    shifted"""
    from petakit5d.detection import correct_xz_offset_data

    correct_xz_offset_data(data, overwrite=True)


def results_paths(data) -> list[str]:
    """each movie's results folder, as the demo builds it:
    ``arrayfun(@(i) [i.source filesep aname filesep], data, 'unif', 0)`` (``source`` ends with
    the separator, so the folder has two in a row, as in MATLAB)"""
    from petakit5d.core.paths import filesep

    return [m.source + filesep + _phase3.ANAME + filesep for m in data]


def detection(data, sigma, apath: list[str]) -> None:
    """Step 4: the spots of every frame (``XR_runDetection3D``: the 'lowSNR' detector with the
    background correction along z, the fits of x, y, z, A and c), the other channel read out at
    them, into ``Detection3D.mat`` and ``Masks/``"""
    from petakit5d.detection import run_detection

    run_detection(
        data,
        sigma=sigma,
        overwrite=True,
        window_size=None,
        results_path=apath,
        mode="xyzAc",
        fit_mixtures=False,
        max_mixtures=1,
        detection_method="lowSNR",
        background_correction=True,
    )


def tracked_movies(data, apath: list[str]):
    """the movies step 5 tracks (``data(~([data.movieLength]==1))``: tracking needs two
    frames), as MATLAB's struct array, with their own results folders"""
    from petakit5d.detection.condition import Condition, struct_array

    keep = [k for k, m in enumerate(data) if m.movie_length != 1]
    return struct_array(Condition(data[k] for k in keep)), [apath[k] for k in keep]


def tracking(movies, apath: list[str]) -> None:
    """Step 5.1: the detections linked into tracks (``XR_runTracking3D`` with
    ``loadTrackSettings('Radius', [3, 6], 'MaxGapLength', 2)``), into ``trackedFeatures.mat``"""
    from petakit5d.detection.tracking import run_tracking, track_settings

    run_tracking(
        movies,
        track_settings(radius=(3, 6), max_gap_length=2),
        file_name="trackedFeatures.mat",
        overwrite=True,
        results_path=apath,
    )


def track_processing(movies, apath: list[str]) -> None:
    """Step 5.2: the tracks processed (``XR_runTrackProcessing3D``): gaps filled, buffers of 3
    frames read out before and after each track, the tracks grouped by gaps, merges and
    splits, into ``ProcessedTracks.mat`` and the Amira files"""
    from petakit5d.detection import run_track_processing

    run_track_processing(
        movies,
        overwrite=True,
        tracker_output="trackedFeatures.mat",
        file_name="ProcessedTracks.mat",
        # frames before and after each track
        buffer=(3, 3),
        buffer_all=False,
        fit_mixtures=False,
        window_size=None,
        results_path=apath,
    )


def rotate(movies, apath: list[str]) -> None:
    """Step 5.3: the processed tracks rotated into the frame of the rotated stacks, coverslip
    horizontal (``rotateTracks3D``), into ``RotatedTracks.mat``"""
    from petakit5d.detection import rotate_tracks

    rotate_tracks(movies, overwrite=True, crop=False, results_path=apath)


def main() -> None:
    steps = _common.Steps(DEMO)
    rt = _phase3.experiment_path(1)
    with steps("Step 1: load data"):
        data = _phase3.load_data(rt)
    _phase3.save_values(DEMO, "load", {"data": data})
    with steps("Step 2: PSF sigmas"):
        sigma = _phase3.psf_sigmas(data, _phase3.psf_paths())
    _phase3.save_values(DEMO, "sigma", {"sigma_mat": sigma})
    with steps("Step 3: deskew"):
        deskew(data)
    with steps("Step 3: XZ offsets"):
        xz_offsets(data)
    apath = results_paths(data)
    with steps("Step 4: detection"):
        detection(data, sigma, apath)
    movies, apath = tracked_movies(data, apath)
    with steps("Step 5.1: tracking"):
        tracking(movies, apath)
    with steps("Step 5.2: track processing"):
        track_processing(movies, apath)
    with steps("Step 5.3: rotate tracks"):
        rotate(movies, apath)


if __name__ == "__main__":
    _common.run_demo(main)
