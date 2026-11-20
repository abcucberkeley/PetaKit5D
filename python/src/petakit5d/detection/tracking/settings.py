"""The tracker's settings as MATLAB's ``loadTrackSettings`` returns them (measured as a black box:
the defaults and the variants of ``'Radius'`` and ``'MaxGapLength'``; its source was not read).

``track_settings(radius=(5, 10), max_gap_length=2)``: the search radius ``[min max]`` of both
cost functions and the time window ``max_gap_length + 1``, which is also the window of the
nearest-neighbour distances and the time at which the Brownian search radius stops growing;
3 standard deviations per frame of the time window. The demos call it with ``radius=(3, 6)``.
"""

from __future__ import annotations

import numpy as np

__all__ = ["track_settings"]


def track_settings(radius=(5, 10), max_gap_length: int = 2) -> dict:
    """``loadTrackSettings('Radius', radius, 'MaxGapLength', max_gap_length)``: a dict with
    ``costMatrices`` (two dicts of ``funcName`` and ``parameters``), ``gapCloseParam`` and
    ``kalmanFunctions``, the fields and values (as MATLAB classes in NumPy) of MATLAB's struct"""
    r = np.asarray(radius, dtype=np.float64).ravel()
    if r.size != 2:
        raise ValueError("radius must be [min max]")
    tw = float(int(max_gap_length) + 1)
    n = int(tw)
    empty = np.zeros((0, 0))
    link = {
        "linearMotion": 0.0,
        "minSearchRadius": float(r[0]),
        "maxSearchRadius": float(r[1]),
        "brownStdMult": 3.0,
        "useLocalDensity": 1.0,
        "nnWindow": tw,
        "kalmanInitParam": empty,
        "diagnostics": empty,
    }
    gaps = {
        "linearMotion": 0.0,
        "minSearchRadius": float(r[0]),
        "maxSearchRadius": float(r[1]),
        "brownStdMult": np.full((n, 1), 3.0),
        "brownScaling": np.array([[0.5, 0.01]]),
        "timeReachConfB": tw,
        "ampRatioLimit": np.array([[0.0, np.inf]]),
        "lenForClassify": 5.0,
        "useLocalDensity": 1.0,
        "nnWindow": tw,
        "linStdMult": np.full((n, 1), 3.0),
        "linScaling": np.array([[1.0, 0.01]]),
        "timeReachConfL": tw,
        "maxAngleVV": 45.0,
        "gapPenalty": empty,
        "resLimit": empty,
    }
    return {
        "costMatrices": [
            {"funcName": "costMatLinearMotionLink2", "parameters": link},
            {"funcName": "costMatLinearMotionCloseGaps2", "parameters": gaps},
        ],
        "gapCloseParam": {
            "timeWindow": tw,
            "mergeSplit": 1.0,
            "minTrackLen": 1.0,
            "diagnostics": 0.0,
        },
        "kalmanFunctions": {
            "reserveMem": "kalmanResMemLM",
            "initialize": "kalmanInitLinearMotion",
            "calcGain": "kalmanGainLinearMotion",
            "timeReverse": "kalmanReverseLinearMotion",
        },
    }
