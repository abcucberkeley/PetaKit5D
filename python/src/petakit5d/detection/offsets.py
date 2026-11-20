"""Channel offsets in x and z between the deskewed frames of a condition's channels (the lab's
``XR_correctXZoffsetData3D`` and ``XZoffsetEstimation3D``, ported from their code).

``xz_offset_estimation_3d(frame_1, frame_2)`` correlates two ``(z, y, x)`` volumes (NaN where
there is no data) over the shifts -7..7 in x and -5..5 in z: for each shift the Pearson
correlation (``matlabcompat.correlation.corr``) of the overlapping voxels that are not NaN in
either volume, taken in MATLAB's order. The offset is the regional maximum of that 15 x 11
matrix (8-connected, the matrix's border left out); of several, the one nearest to the
matrix's middle as MATLAB measures it (``(j - 11/2)^2 + (i - 15/2)^2`` on 1-based indices,
which leans towards negative shifts), the first in MATLAB's order on a tie; none gives
``(0, 0)`` (``pick_offset(c)``: this step on a given matrix). A NaN correlation (constant or
NaN-only overlaps) is MATLAB's ``imregionalmax`` error.

``correct_xz_offset_data(data)`` estimates the offsets of the channels against channel 1 on
three frames (the first, ``round((nF - 1)/2)`` and the last; all of them below three), takes
``fix(median)`` over them, and pads every channel's DS frames (single, written to
``<DS>_tmp_xz_correct/``, which then replaces ``<DS>``) so that the channels line up, the
padding on the side the offsets need; ``<DS>/xz_offset.mat`` (v7.3: ``xzoff_mat``, ``pad_mat``,
``channels``) records it in each channel. MATLAB's behaviour is kept in every detail, its quirks
included:

- only channel 1's ``xz_offset.mat`` decides whether a movie is done (without ``overwrite``);
- the lock ``<DS of channel 1>/xz_offset.tmp`` makes a second run skip the movie; ``*.tmp``
  files under that folder older than ``5 * nF * nCh`` minutes are removed first (MATLAB's
  ``chunk_lock_clean``; MATLAB misses them when the path holds the text ``tmp`` elsewhere);
- channels are sorted by name for the estimate (``data`` is returned as given);
- with one chosen frame and three or more channels the median runs along the row (MATLAB's
  ``median`` of a row vector), and the padding step fails as in MATLAB;
- ``overwrite`` estimates again from the frames as they are (padded ones included) and pads again;
- a rewritten frame is named ``<name>.tif`` whatever its extension was.
"""

from __future__ import annotations

import os
import shutil
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

__all__ = [
    "XZOFF_SHIFTS",
    "ZOFF_SHIFTS",
    "OffsetError",
    "correct_xz_offset_data",
    "pick_offset",
    "xz_correlations",
    "xz_offset_estimation_3d",
]

XZOFF_SHIFTS = tuple(range(-7, 8))
ZOFF_SHIFTS = tuple(range(-5, 6))


class OffsetError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _ranges(off: int, n: int) -> tuple[slice, slice]:
    """MATLAB's index ranges of the two volumes for a shift (``xa_1``, ``xa_2``), 0-based"""
    if off <= 0:
        return slice(0, max(n + off, 0)), slice(min(-off, n), n)
    return slice(min(off, n), n), slice(0, max(n - off, 0))


def _matlab_layout(f: np.ndarray) -> np.ndarray:
    """a (z, y, x) volume as a C array (z, x, y): its C order is MATLAB's linear order"""
    return np.ascontiguousarray(np.swapaxes(f, 1, 2))


def xz_correlations(frame_1: Any, frame_2: Any) -> np.ndarray:
    """The 15 x 11 matrix ``xcorr_mat`` of ``XZoffsetEstimation3D``: row i is the x shift
    ``XZOFF_SHIFTS[i]``, column j the z shift ``ZOFF_SHIFTS[j]``."""
    from petakit5d.matlabcompat.correlation import corr

    f1 = np.asarray(frame_1, dtype=np.float64)
    f2 = np.asarray(frame_2, dtype=np.float64)
    if f1.ndim == 2:
        f1, f2 = f1[None], f2[None]
    nz, _, nx = f1.shape
    m1, m2 = _matlab_layout(f1), _matlab_layout(f2)
    exact_nan = bool(np.isinf(m1).any() or np.isinf(m2).any())
    n1, n2 = np.isnan(m1), np.isnan(m2)
    out = np.zeros((len(XZOFF_SHIFTS), len(ZOFF_SHIFTS)))
    for i, xoff in enumerate(XZOFF_SHIFTS):
        xa1, xa2 = _ranges(xoff, nx)
        for j, zoff in enumerate(ZOFF_SHIFTS):
            za1, za2 = _ranges(zoff, nz)
            a = m1[za1, xa1]
            b = m2[za2, xa2]
            if exact_nan:  # MATLAB's ~isnan(a + b) also drops Inf + -Inf
                with np.errstate(invalid="ignore"):
                    keep = ~np.isnan(a + b)
            else:
                keep = ~(n1[za1, xa1] | n2[za2, xa2])
            out[i, j] = corr(a[keep], b[keep])
    return out


def xz_offset_estimation_3d(frame_1: Any, frame_2: Any, zoff_mat: Any = None) -> tuple[int, int]:
    """``XZoffsetEstimation3D(frame_1, frame_2)``: the x and z offsets of ``frame_2`` against
    ``frame_1`` (module docstring). Volumes are ``(z, y, x)``; NaN marks voxels without data.
    MATLAB's third argument leaves its x shifts undefined: giving one raises its error."""
    if zoff_mat is not None:
        raise OffsetError(
            "Unrecognized function or variable 'xoff_mat'.", "MATLAB:UndefinedFunction"
        )
    return pick_offset(xz_correlations(frame_1, frame_2))


def pick_offset(c: Any) -> tuple[int, int]:
    """The offsets of a correlation matrix (rows the x shifts, columns the z shifts): its regional
    maxima (``imregionalmax``, 8-connected) off the border, the one nearest to the middle as
    MATLAB measures it, the first in MATLAB's order on a tie; (0, 0) without one."""
    from petakit5d.tools.puncta_removal import imregionalmax

    c = np.asarray(c, dtype=np.float64)
    try:
        bw = imregionalmax(c)
    except RuntimeError as err:  # MATLAB:imregionalmax:expectedNonNaN
        raise OffsetError(str(err), getattr(err, "identifier", "")) from err
    bw[:, [0, -1]] = False
    bw[[0, -1], :] = False
    lz, lx = np.nonzero(bw.T)  # MATLAB's find: column-major
    if lx.size == 0:
        return 0, 0
    if lx.size > 1:
        d = (lz + 1 - c.shape[1] / 2) ** 2 + (lx + 1 - c.shape[0] / 2) ** 2
        k = int(np.argmin(d))
        lx, lz = lx[k : k + 1], lz[k : k + 1]
    return XZOFF_SHIFTS[int(lx[0])], ZOFF_SHIFTS[int(lz[0])]


def _read_frame(path: str, flip: bool) -> np.ndarray:
    from petakit5d.core import dtypes
    from petakit5d.io import read_tiff

    f = np.asarray(read_tiff(path), dtype=np.float64)
    if f.ndim == 2:
        f = f[None]
    f[f == 0] = dtypes.NAN
    return f[::-1] if flip else f


def _lock_clean(path: str, minutes: float) -> None:
    """MATLAB's ``chunk_lock_clean(path, duration)``: remove the ``*.tmp`` files under
    ``path`` older than ``minutes``."""
    now = time.time()
    for root, _, files in os.walk(path):
        for f in files:
            if f.endswith(".tmp"):
                p = os.path.join(root, f)
                try:
                    if (now - os.path.getmtime(p)) / 60 > minutes:
                        os.remove(p)
                except OSError:
                    pass


def _median_rows(m: np.ndarray) -> np.ndarray:
    """MATLAB's ``median(A)`` of a frames x channels matrix: along the rows, but along the
    only row of a 1 x n matrix"""
    from petakit5d.matlabcompat.stats import median

    if m.shape[0] == 1:
        return np.atleast_1d(np.asarray(median(m.reshape(-1)), dtype=np.float64)).reshape(-1)
    return np.asarray(median(m, axis=0), dtype=np.float64).reshape(-1)


def _short_path(source: str, n: int = 3) -> str:
    from petakit5d.detection.condition import short_path

    return short_path(source, n)


def _sorted_channels(movie) -> Any:
    import copy

    ch = list(movie.channels)
    if ch == sorted(ch):
        return movie
    order = sorted(range(len(ch)), key=lambda i: ch[i])
    m = copy.copy(movie)
    m.channels = [ch[i] for i in order]
    m.frame_paths = [movie.frame_paths[i] for i in order]
    m.frame_paths_ds = [movie.frame_paths_ds[i] for i in order]
    m.frame_paths_dsr = [movie.frame_paths_dsr[i] for i in order]
    return m


def _pad_frames(paths: Sequence[str], ds_path: str, pad: np.ndarray) -> None:
    """write each frame padded by ``pad`` (rows x, z; columns before, after) as single into
    ``<ds_path>_tmp_xz_correct/``, then replace ``ds_path`` by that folder"""
    from petakit5d.io import read_tiff, write_tiff

    out_dir = ds_path + "_tmp_xz_correct"
    os.makedirs(out_dir, exist_ok=True)
    (xb, xt), (zb, zt) = (int(pad[0, 0]), int(pad[0, 1])), (int(pad[1, 0]), int(pad[1, 1]))
    for p in paths:
        im = np.asarray(read_tiff(p))
        if im.ndim == 2:
            raise OffsetError(
                "Arrays have incompatible sizes for this operation.",
                "MATLAB:sizeDimensionsMustMatch",
            )
        nz, ny, nx = im.shape
        im2 = np.zeros((nz + zb + zt, ny, nx + xb + xt), dtype=np.float32)
        im2[zb : zb + nz, :, xb : xb + nx] = im.astype(np.float64)
        name = os.path.splitext(os.path.basename(p))[0]
        write_tiff(im2, os.path.join(out_dir, name + ".tif"))
    shutil.rmtree(ds_path)
    os.replace(out_dir, ds_path)


def _save_flag(path: str, xzoff: np.ndarray, pad_mat: np.ndarray, channels: Sequence[str]) -> None:
    from petakit5d.io import write_mat

    write_mat(
        path,
        {
            "xzoff_mat": xzoff,  # MATLAB's size; write_mat takes the values by index
            "pad_mat": pad_mat,
            "channels": list(channels),
        },
        version="7.3",
    )


def correct_xz_offset_data(
    data: Sequence[Any],
    *,
    overwrite: bool = False,
    ds: bool = True,
    s_cmos_camera_flip: bool = False,
) -> Sequence[Any]:
    """``XR_correctXZoffsetData3D(data, 'Overwrite', overwrite, 'DS', ds, 'sCMOSCameraFlip',
    flip)``: line the channels' DS frames up in x and z (module docstring). ``data`` is the
    condition (``condition.Movie`` items); it is returned unchanged."""
    for flag, name in (
        (overwrite, "Overwrite"),
        (ds, "DS"),
        (s_cmos_camera_flip, "sCMOSCameraFlip"),
    ):
        if not isinstance(flag, bool | np.bool_):
            raise OffsetError(
                f"The value of '{name}' is invalid. It must satisfy the function: islogical.",
                "MATLAB:InputParser:ArgumentFailedValidation",
            )
    for movie in data:
        m = _sorted_channels(movie)
        n_ch = len(m.channels)
        n_f = len(m.frame_paths[0])
        print(f"Correct xz-offset for '{_short_path(m.source)}'")
        chosen = list(range(1, n_f + 1)) if n_f < 3 else [1, int(_round((n_f - 1) / 2)), n_f]
        xzoff = np.zeros((len(chosen), n_ch - 1, 2))
        ds_path = os.path.dirname(m.frame_paths_ds[0][0])
        if os.path.isfile(os.path.join(ds_path, "xz_offset.mat")) and not overwrite:
            continue
        _lock_clean(ds_path, 5 * len(m.frame_paths_ds[0]) * n_ch)
        tmp = os.path.join(ds_path, "xz_offset.tmp")
        if os.path.exists(tmp):
            continue
        open(tmp, "w").close()
        src = m.frame_paths_ds if ds else m.frame_paths
        for k, f in enumerate(chosen):
            prim = _read_frame(src[0][f - 1], s_cmos_camera_flip)
            for j in range(1, n_ch):
                other = _read_frame(src[j][f - 1], s_cmos_camera_flip)
                xzoff[k, j - 1, :] = xz_offset_estimation_3d(prim, other)
        xmove = np.concatenate(([0.0], np.trunc(_median_rows(xzoff[:, :, 0]))))
        zmove = np.concatenate(([0.0], np.trunc(_median_rows(xzoff[:, :, 1]))))
        xpad = np.abs(xmove[:, None] - np.array([xmove.min(), xmove.max()])[None, :])
        zpad = np.abs(zmove[:, None] - np.array([zmove.min(), zmove.max()])[None, :])
        pad_mat = np.stack((xpad, zpad), axis=2)  # nCh x [before after] x [x z]
        for j in range(n_ch - 1, -1, -1):
            ds_j = os.path.dirname(m.frame_paths_ds[j][0])
            flag_j = os.path.join(ds_j, "xz_offset.mat")
            if os.path.isfile(flag_j) and not overwrite:
                continue
            if j >= pad_mat.shape[0]:
                raise OffsetError(
                    "Index in position 1 exceeds array bounds.", "MATLAB:badsubscript"
                )
            cur = pad_mat[j].T  # rows x, z; columns before, after
            if np.all(cur == 0):
                _save_flag(flag_j, xzoff, pad_mat, m.channels)
                continue
            _pad_frames(m.frame_paths_ds[j][:n_f], ds_j, cur)
            _save_flag(flag_j, xzoff, pad_mat, m.channels)
        if os.path.exists(tmp):
            os.remove(tmp)
        print("Done!", end="")
    return data


def _round(x: float) -> float:
    from petakit5d.core.dtypes import mround

    return float(mround(np.float64(x)))
