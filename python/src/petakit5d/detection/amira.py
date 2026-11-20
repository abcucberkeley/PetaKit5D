"""``XR_amiraWriteTracks`` (the lab's): one AmiraMesh ASCII file per frame,
``<basename>_t_%04d.am``, with the tracks present in that frame as an ``HxSpatialGraph``: each
track one edge from its first position to its position in the frame, the positions up to that
frame as the edge's points, and its id, lifetime, amplitude and category.

As MATLAB writes them, byte for byte (the help text and the code of the lab's file):

- a track is on in frame ``k`` when ``start < k <= end``; its points are its first ``k - start
  + 1`` columns;
- the positions are ``(p - 1) * s`` with ``s = [1 1 -zAniso]`` and ``nz * zAniso`` added to z
  (``FlipZ``, ``nz`` the slices of the first DS frame), else ``s = [1 1 zAniso]``; MATLAB indexes
  the channel x frame position matrices linearly here (``tr.x(k)``, column-major), so with two
  channels the edge's end and its points take alternate channels' values (kept);
- the amplitude is channel ``Ch``'s, minus 1, times ``AmpMultiplier``;
- numbers as ``dlmwrite(..., 'precision', 16)`` writes them (``%.16g``, rows of values
  separated by a space, ``\\n`` line ends).
"""

from __future__ import annotations

import os

import numpy as np

__all__ = ["frame_text", "write_tracks"]

_HEADER = (
    "# AmiraMesh 3D ASCII 2.0\n\n"
    "define VERTEX %(nv)d\n"
    "define EDGE %(ne)d \n"
    "define POINT %(np)d\n"
    "\n"
    'Parameters { ContentType "HxSpatialGraph"}\n\n'
    "VERTEX { float[3] VertexCoordinates } @1\n"
    "EDGE { int[2] EdgeConnectivity } @2\n"
    "EDGE { int NumEdgePoints } @3\n"
    "POINT { float[3] EdgePointCoordinates } @4\n"
    "EDGE { int trackId } @5\n"
    "VERTEX {int startEnd } @6\n"
    "EDGE { int lifetime } @7\n"
    "POINT { int Amplitude } @8\n"
    "EDGE { int catIdx } @9\n"
)


def _g(v: float) -> str:
    """MATLAB's ``sprintf('%.16g', v)``"""
    v = float(v)
    if v != v:
        return "NaN"
    if v in (float("inf"), float("-inf")):
        return "Inf" if v > 0 else "-Inf"
    return f"{v:.16g}"


def _rows(a: np.ndarray) -> str:
    a = np.asarray(a, dtype=np.float64)
    if a.ndim == 1:
        a = a.reshape(-1, 1)
    return "".join(" ".join(_g(v) for v in row) + "\n" for row in a)


def _lin(m, idx):
    """``m(idx)`` of a MATLAB matrix: linear, column-major, 1-based"""
    a = np.asarray(m, dtype=np.float64)
    a = a.reshape(1, -1) if a.ndim < 2 else a
    return a.T.ravel()[np.atleast_1d(np.asarray(idx, dtype=np.int64)) - 1]


def frame_text(tracks, k: int, s, z_offset: float, ch: int = 1, amp_multiplier: float = 1.0) -> str:
    """the file of frame ``k`` (1-based) for ``tracks`` (dicts with x, y, z, A, start, end,
    lifetime_s, catIdx)"""
    s = [float(v) for v in s]
    starts = np.array([float(t["start"]) for t in tracks])
    ends = np.array([float(t["end"]) for t in tracks])
    on = np.nonzero((ends >= k) & (k > starts))[0]
    n = on.size
    text = _HEADER % {"nv": n * 2, "ne": n, "np": int(np.sum(k - starts[on] + 1))}
    if n == 0:
        return text
    ends_xyz = np.empty((2 * n, 3))
    pts, amp, npts = [], [], []
    for c, i in enumerate(on):
        t = tracks[i]
        e = int(k - starts[i] + 1)
        ends_xyz[2 * c] = [
            (_lin(t["x"], 1)[0] - 1) * s[0],
            (_lin(t["y"], 1)[0] - 1) * s[1],
            (_lin(t["z"], 1)[0] - 1) * s[2] + z_offset,
        ]
        ends_xyz[2 * c + 1] = [
            (_lin(t["x"], e)[0] - 1) * s[0],
            (_lin(t["y"], e)[0] - 1) * s[1],
            (_lin(t["z"], e)[0] - 1) * s[2] + z_offset,
        ]
        idx = np.arange(1, e + 1)
        pts.append(
            np.stack(
                [
                    (_lin(t["x"], idx) - 1) * s[0],
                    (_lin(t["y"], idx) - 1) * s[1],
                    (_lin(t["z"], idx) - 1) * s[2] + z_offset,
                ],
                axis=1,
            )
        )
        A = np.asarray(t["A"], dtype=np.float64)
        A = A.reshape(1, -1) if A.ndim == 1 else A
        amp.append((A[ch - 1, :e] - 1) * amp_multiplier)
        npts.append(e)
    start_end = np.ones(2 * n)
    start_end[0::2] = 0
    edges = np.arange(2 * n, dtype=np.float64).reshape(n, 2)
    parts = [
        text,
        "\n@1\n", _rows(ends_xyz),
        "\n@2\n", _rows(edges),
        "\n@3\n", _rows(np.array(npts, dtype=np.float64)),
        "\n@4\n", _rows(np.vstack(pts)),
        "\n@5\n", _rows(on + 1.0),
        "\n@6\n", _rows(start_end),
        "\n@7\n", _rows(np.array([float(tracks[i]["lifetime_s"]) for i in on])),
        "\n@8\n", _rows(np.concatenate(amp)),
        "\n@9\n", _rows(np.array([float(tracks[i]["catIdx"]) for i in on])),
    ]  # fmt: skip
    return "".join(parts)


def write_tracks(
    filename: str,
    tracks,
    movie_length: int,
    z_aniso: float,
    n_slices: int,
    *,
    flip_z: bool = True,
    scales=None,
    ch: int = 1,
    amp_multiplier: float = 1.0,
) -> list[str]:
    """``XR_amiraWriteTracks(filename, tracks, data, 'FlipZ', flip_z, ...)``: frames 1 to
    ``movie_length``; ``n_slices`` is the first DS frame's slice count (MATLAB reads the frame).
    Returns the files written."""
    if flip_z:
        s = list(scales) if scales is not None else [1.0, 1.0, -float(z_aniso)]
        z_offset = n_slices * float(z_aniso)
    else:
        s = list(scales) if scales is not None else [1.0, 1.0, float(z_aniso)]
        z_offset = 0.0
    folder, name = os.path.split(filename)
    name = os.path.splitext(name)[0]
    os.makedirs(folder, exist_ok=True)
    out = []
    for k in range(1, int(movie_length) + 1):
        path = os.path.join(folder, f"{name}_t_{k:04d}.am")
        with open(path, "w", newline="\n") as fh:
            fh.write(frame_text(tracks, k, s, z_offset, ch, amp_multiplier))
        out.append(path)
    return out
