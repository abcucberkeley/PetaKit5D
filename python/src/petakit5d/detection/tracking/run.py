"""``XR_runTracking3D``: track the detections of each movie (``<results>/Detection3D.mat``) and
write ``<results>/trackedFeatures.mat`` with the variables MATLAB writes there:
``costMatrices``, ``gapCloseParam``, ``kalmanFunctions`` (the settings), ``tracksFinal``
(n x 1 struct array: ``tracksFeatIndxCG``, ``tracksCoordAmpCG``, ``seqOfEvents``) and
``kalmanInfoLink`` (1 x frames struct array of the last linking pass's filter state).

The wrapper's behaviour comes from its help text and its measured outputs (the file's header
names an outside author, so its code was not used). ``ParallelMethod`` is accepted: the movies
are tracked one after the other in this process.
"""

from __future__ import annotations

import os

from petakit5d.params import Param, register

__all__ = ["run_tracking", "track_detections", "write_tracked_features"]


def _movies(data) -> list[dict]:
    if isinstance(data, dict):
        return [data]
    return [dict(d) for d in data]


def track_detections(frame_info, settings: dict | None = None) -> tuple[list[dict], list[dict]]:
    """the tracks of one movie's detections (``frameInfo`` of ``Detection3D.mat``, a sequence
    of per-frame dicts): MATLAB's ``tracksFinal`` elements and ``kalmanInfoLink`` per frame"""
    from petakit5d.detection.tracking.settings import track_settings
    from petakit5d.detection.tracking.tracks import (
        frames_from_frame_info,
        kalman_to_matlab,
        track_frames,
    )

    if settings is None:
        settings = track_settings(radius=(3, 6), max_gap_length=2)
    frames = frames_from_frame_info(list(frame_info))
    tracks, kal = track_frames(frames, settings)
    return tracks, kalman_to_matlab(kal, frames)


def write_tracked_features(
    path: str, settings: dict, tracks: list[dict], kalman: list[dict]
) -> None:
    """``trackedFeatures.mat`` as MATLAB writes it (v7), the feature indices double as in MATLAB.
    MATLAB's save stores integer-valued doubles as small integers inside the file, the port as
    doubles; both load as double (``docs/parity/detection.md``, Tracking)."""
    from petakit5d.io.matfile import StructArray, write_mat

    cm = settings["costMatrices"]
    cm_rows = [{"funcName": c["funcName"], "parameters": dict(c["parameters"])} for c in cm]
    tf = StructArray(
        (len(tracks), 1), ("tracksFeatIndxCG", "tracksCoordAmpCG", "seqOfEvents"), tracks
    )
    if not tracks:
        tf = StructArray.empty(("tracksFeatIndxCG", "tracksCoordAmpCG", "seqOfEvents"), (0, 1))
    kl = StructArray(
        (1, len(kalman)), ("stateVec", "stateCov", "noiseVar", "stateNoise", "scheme"), kalman
    )
    out = {
        "costMatrices": StructArray((1, len(cm_rows)), ("funcName", "parameters"), cm_rows),
        "gapCloseParam": dict(settings["gapCloseParam"]),
        "kalmanFunctions": dict(settings["kalmanFunctions"]),
        "tracksFinal": tf,
        "kalmanInfoLink": kl,
    }
    write_mat(path, out, version="7")


@register(
    "XR_runTracking3D",
    # third-party code: the port's own checks (accepts=), made from MATLAB's behaviour
    # (tests/params/test_third_party_checks.py)
    Param(
        "data",
        None,
        "Movies (the structure XR_loadConditionData3D returns).",
        convert="same",
        accepts="struct",
    ),
    Param(
        "settings",
        None,
        "Tracker settings (default: loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2)).",
        accepts="struct",
    ),
    Param("Overwrite", None, "Track movies whose result exists again.", accepts="logical"),
    Param("FileName", None, "Result file name.", accepts="char"),
    Param(
        "Frames",
        None,
        "Frames to track (0-based; default all).",
        convert="index1",
        accepts="vector",
    ),
    Param("DetectionFile", None, "Detection file name in the results folder.", accepts="char"),
    Param(
        "ResultsPath",
        None,
        "Results folder of each movie (default <source>/Analysis/).",
        matlab_default="arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0)",
        accepts="cell",
    ),
    Param(
        "ParallelMethod",
        None,
        "Accepted; movies are tracked in this process.",
        accepts="char",
    ),
    task="C03",
)
def run_tracking(
    data,
    settings=None,
    *,
    overwrite: bool = False,
    file_name: str = "trackedFeatures",
    frames=None,
    detection_file: str = "Detection3D.mat",
    results_path=None,
    parallel_method: str = "lock_distributed",
) -> None:
    """``XR_runTracking3D``: track each movie's detections into ``<results>/<file_name>``"""
    from petakit5d.io.matfile import read_mat

    del parallel_method
    movies = _movies(data)
    if results_path is None:
        results_path = [os.path.join(str(m["source"]), "Analysis") + os.sep for m in movies]
    elif isinstance(results_path, str):
        results_path = [results_path]
    name = str(file_name)
    if not os.path.splitext(name)[1]:
        name += ".mat"
    if settings is None:
        from petakit5d.detection.tracking.settings import track_settings

        settings = track_settings(radius=(3, 6), max_gap_length=2)
    for k, _m in enumerate(movies):
        rp = str(results_path[k])
        out = os.path.join(rp, name)
        if os.path.exists(out) and not overwrite:
            continue
        fi = read_mat(os.path.join(rp, str(detection_file)), "frameInfo")["frameInfo"]
        fi = [fi] if isinstance(fi, dict) else list(fi)
        if frames is not None:
            fi = [fi[int(i)] for i in (frames if hasattr(frames, "__iter__") else [frames])]
        tracks, kal = track_detections(fi, settings)
        tmp = out + ".tmp_" + str(os.getpid())
        write_tracked_features(tmp, settings, tracks, kal)
        os.replace(tmp, out)
