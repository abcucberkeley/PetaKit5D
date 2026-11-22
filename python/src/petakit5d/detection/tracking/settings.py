"""The tracker's settings as MATLAB's ``loadTrackSettings`` returns them (measured as a black box:
the defaults and the variants of ``'Radius'``, ``'GapRadius'``, ``'LinkRadius'``,
``'MaxGapLength'`` and ``'CostType'``; its source was not read).

``track_settings(radius=(5, 10), max_gap_length=2)``: the search radius ``[min max]`` of both
cost functions and the time window ``max_gap_length + 1``, which is also the window of the
nearest-neighbour distances and the time at which the Brownian search radius stops growing;
3 standard deviations per frame of the time window. The demos call it with ``radius=(3, 6)``.

Registered with MATLAB's interface (third-party code: its header names an outside author; MATLAB
validates none of its options). The search radii (measured): ``link_radius`` is the linking cost
function's ``[min max]`` (``costMatrices(1)``), ``gap_radius`` the gap closing one's
(``costMatrices(2)``), ``[5 10]`` each by default; a ``radius`` that is not empty (None is
MATLAB's ``[]``) is both and wins over them. Of each, the first two elements are read in MATLAB's
order; fewer fail (``TrackSettingsError``, ``MATLAB:badsubscript``). ``CostType`` is taken at
``'linear'`` only (in any case, as MATLAB): MATLAB's other values give the same settings without
the cost functions' names, which tracking cannot use.
"""

from __future__ import annotations

import numpy as np

from petakit5d.params import Param, register

__all__ = ["TrackSettingsError", "track_settings"]

_DEFAULT_RADIUS = (5.0, 10.0)


class TrackSettingsError(IndexError):
    """An error MATLAB's loadTrackSettings raises; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _pair(value, name: str) -> tuple[float, float]:
    """a search radius [min max]: its first two elements in MATLAB's order (more are not read,
    fewer fail as MATLAB's indexing does)"""
    a = np.asarray(value, dtype=np.float64).ravel(order="F")
    if a.size < 2:
        raise TrackSettingsError(
            f"track_settings: {name} needs two elements, [min max]; it has {a.size} "
            "(MATLAB's indexing fails).",
            "MATLAB:badsubscript",
        )
    return float(a[0]), float(a[1])


@register(
    "loadTrackSettings",
    Param(
        "Radius",
        None,
        "Search radius [min max] of both cost functions (None: GapRadius and LinkRadius).",
    ),
    Param("GapRadius", None, "Search radius [min max] of the gap closing (without Radius)."),
    Param("LinkRadius", None, "Search radius [min max] of the linking (without Radius)."),
    Param("MaxGapLength", None, "Longest gap closed, in frames (the time window less one)."),
    Param(
        "CostType",
        None,
        "Not ported: 'linear' only (in any case): MATLAB's other values give settings without "
        "the cost functions' names, which tracking cannot use.",
    ),
    task="C03",
)
def track_settings(
    *,
    radius=None,
    gap_radius=_DEFAULT_RADIUS,
    link_radius=_DEFAULT_RADIUS,
    max_gap_length: int = 2,
    cost_type: str = "linear",
) -> dict:
    """``loadTrackSettings('Radius', radius, 'GapRadius', gap_radius, 'LinkRadius',
    link_radius, 'MaxGapLength', max_gap_length)``: a dict with ``costMatrices`` (two dicts of
    ``funcName`` and ``parameters``), ``gapCloseParam`` and ``kalmanFunctions``, the fields and
    values (as MATLAB classes in NumPy) of MATLAB's struct (module docstring)"""
    if not isinstance(cost_type, str) or cost_type.lower() != "linear":
        raise NotImplementedError(
            f"track_settings: CostType {cost_type!r}: only the default is ported"
        )
    if radius is None or np.size(radius) == 0:
        link_r = _pair(link_radius, "LinkRadius")
        gap_r = _pair(gap_radius, "GapRadius")
    else:
        link_r = gap_r = _pair(radius, "Radius")
    tw = float(int(max_gap_length) + 1)
    n = int(tw)
    empty = np.zeros((0, 0))
    link = {
        "linearMotion": 0.0,
        "minSearchRadius": link_r[0],
        "maxSearchRadius": link_r[1],
        "brownStdMult": 3.0,
        "useLocalDensity": 1.0,
        "nnWindow": tw,
        "kalmanInitParam": empty,
        "diagnostics": empty,
    }
    gaps = {
        "linearMotion": 0.0,
        "minSearchRadius": gap_r[0],
        "maxSearchRadius": gap_r[1],
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
