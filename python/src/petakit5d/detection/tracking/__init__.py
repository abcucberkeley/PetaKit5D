"""Tracking of point detections: a new implementation of the u-track method (Jaqaman et al.
2008: frame-to-frame linking and gap closing with merging and splitting as linear assignment
problems, with a Kalman filter for Brownian motion), as MATLAB PetaKit5D runs it for the point
detection demos (``XR_runTracking3D`` with ``loadTrackSettings('Radius', [3 6],
'MaxGapLength', 2)``).

Written from the published method and from MATLAB's behaviour measured as a black box (the
inputs and outputs of each function, recorded on the demo data and on chosen inputs); no code of
the GPL sources is used. ``docs/detection_scope.md`` (section 3.5) lists the stages.

====================================  ===============================================
``run_tracking``                      XR_runTracking3D (Detection3D.mat -> trackedFeatures.mat)
``track_detections``                  trackCloseGapsKalmanSparse on a movie's frameInfo
``track_settings``                    loadTrackSettings
``kd_tree_ball_query``                KDTreeBallQuery (mex)
``create_distance_matrix``            createDistanceMatrix (mex)
``solve_assignment``                  lap / mexLap (SciPy's LAPJVsp)
====================================  ===============================================

Modules: ``linking`` (the Kalman filter and the linking passes), ``gapclosing`` (the gap,
merge and split costs), ``tracks`` (the passes, gap closing and the compound tracks), ``lap``,
``mex``, ``settings``, ``run``. Nothing is imported until a function is first used.
"""

_EXPORTS = {
    "run_tracking": "run",
    "track_detections": "run",
    "write_tracked_features": "run",
    "track_settings": "settings",
    "kd_tree_ball_query": "mex",
    "create_distance_matrix": "mex",
    "solve_assignment": "lap",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    if name in _EXPORTS:
        import importlib

        mod = importlib.import_module(f"{__name__}.{_EXPORTS[name]}")
        return getattr(mod, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(_EXPORTS))
