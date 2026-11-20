"""``XR_runDetection3D``: detection of every frame of a condition, the other channels read out at
the master channel's detections, and the results in MATLAB's files (``Detection3D.mat`` with the
struct array ``frameInfo``, ``Masks/dmask_<k>.tif``, an empty ``Amira/Detection`` folder).

Measured as a black box from traced MATLAB runs and designed-input runs (task C02; no GPL
source read): bit for bit on 14 real frames when given MATLAB's detector and fit outputs.
``frameInfo`` keeps MATLAB's values (1-based positions), as the file holds them.

Per frame k (1-based) of the master channel (the channel whose folder is ``data['source']``):
the DS frame as double with zeros as NaN; ``estimate_background(vol, k)``; the lowSNR detector
with ``'BackgroundInfo' [bg, thr]`` (``thr`` the camera's residual threshold, 3 for the demo's);
then, for each other channel, ``'Ac'`` fits at the master's positions and ``'xyzAc'`` fits from
them (default options: windows ``ceil(2 sigma)``, ``kLevel = norminv(0.975)``), and per entry
the choice between them (``other_channel``).
"""

from __future__ import annotations

import os

import numpy as np

__all__ = [
    "FIELD_ORDER",
    "assemble_frame",
    "dmask_name",
    "frame_volume",
    "run_detection",
    "runner_sigma",
]

# frameInfo's fields in MATLAB's order
FIELD_ORDER = (
    "x",
    "y",
    "z",
    "A",
    "s",
    "c",
    "x_pstd",
    "y_pstd",
    "z_pstd",
    "A_pstd",
    "c_pstd",
    "x_init",
    "y_init",
    "z_init",
    "sigma_r",
    "SE_sigma_r",
    "RSS",
    "pval_Ar",
    "hval_Ar",
    "hval_AD",
    "isPSF",
    "xCoord",
    "yCoord",
    "zCoord",
    "amp",
    "dRange",
)
# fields with one row per channel
ROW_FIELDS = (
    "x",
    "y",
    "z",
    "A",
    "c",
    "x_pstd",
    "y_pstd",
    "z_pstd",
    "A_pstd",
    "c_pstd",
    "sigma_r",
    "SE_sigma_r",
    "RSS",
    "pval_Ar",
    "hval_Ar",
    "hval_AD",
)
_LOGICAL = ("hval_Ar", "hval_AD", "isPSF")
_POS = ("x", "y", "z", "x_init", "y_init", "z_init")


def frame_volume(raw) -> np.ndarray:
    """a DS frame as the runner uses it: double, zeros as NaN"""
    v = np.asarray(raw).astype(np.float64)
    v[v == 0] = np.nan
    return v


def runner_sigma(sigma) -> np.ndarray:
    """the runner's sigma check ``if any(sigma < 1.1)`` with MATLAB's meaning: ``any`` along the
    first non-singleton dimension (per column of an nCh x 2 matrix), and ``if`` of a vector is
    true only when every element is; then every value below 1.1 becomes 1.1"""
    s = np.array(sigma, dtype=np.float64, ndmin=2)
    lt = s < 1.1
    dims = [d for d in range(lt.ndim) if lt.shape[d] != 1]
    v = np.atleast_1d(lt.any(axis=dims[0]) if dims else lt.ravel())
    if v.size and bool(np.all(v)):
        s = np.where(lt, 1.1, s)
    return s


def to_matlab_fields(ps: dict) -> dict:
    """detections (``detect.point_source_detection``'s form) as MATLAB's struct: 1 x n rows,
    1-based positions (exact: x + 1 of a position)"""
    out = {}
    for k, v in ps.items():
        a = np.asarray(v)
        if k in _POS:
            a = a + 1
        out[k] = a if a.ndim == 2 else a.reshape(1, -1)
    return out


def _mmin(v) -> int:
    """MATLAB's ``[~, i] = min(v)``: the first minimum, NaN skipped (the first entry if all NaN)"""
    v = np.asarray(v, dtype=np.float64)
    ok = ~np.isnan(v)
    if not ok.any():
        return 0
    return int(np.nonzero(v == np.min(v[ok]))[0][0])


def _ball_query(P: np.ndarray, radius: float) -> list[np.ndarray]:
    """per point, the points within ``radius`` (itself included), by increasing distance (ties
    by index)"""
    from scipy.spatial import cKDTree

    out = []
    if P.shape[0] == 0:
        return out
    for i, b in enumerate(cKDTree(P).query_ball_point(P, radius)):
        b = np.asarray(b, dtype=np.int64)
        d = np.sqrt(np.sum((P[b] - P[i]) ** 2, axis=1))
        b, d = b[d <= radius], d[d <= radius]
        out.append(b[np.lexsort((b, d))])
    return out


def other_channel(
    det: dict, ac: dict, xyz: dict, sigma_xy: float, radius: float = 0.5
) -> np.ndarray:
    """which fit each master detection takes in another channel (MATLAB forms of the fields):
    -1 its Ac fit, j >= 0 the xyzAc fit of entry j. Candidates: xyzAc fits whose xy shift from
    the master's position is below ``2 sigma_xy`` and whose A exceeds the Ac fit's; in every ball
    of ``radius`` among the candidates (in query order) with more than one member, the lowest-RSS
    fit goes to the member whose master position is nearest to it (later balls overwrite) and
    the other members take their Ac fit."""
    F = lambda d, k: np.asarray(d[k], dtype=np.float64).ravel()  # noqa: E731
    M = np.stack([F(det, "x"), F(det, "y"), F(det, "z")], axis=1)
    Q = np.stack([F(xyz, "x"), F(xyz, "y"), F(xyz, "z")], axis=1)
    n = M.shape[0]
    src = np.full(n, -1, dtype=np.int64)
    with np.errstate(invalid="ignore"):
        dxy = np.sqrt((Q[:, 0] - M[:, 0]) ** 2 + (Q[:, 1] - M[:, 1]) ** 2)
        cand = np.nonzero((dxy < 2.0 * float(sigma_xy)) & (F(xyz, "A") > F(ac, "A")))[0]
    rss = F(xyz, "RSS")
    src[cand] = cand
    marked = np.zeros(n, dtype=bool)
    for o in _ball_query(Q[cand], radius):
        if len(o) < 2:
            continue
        ent = cand[o]
        best = ent[_mmin(rss[ent])]
        rcv = ent[_mmin(np.sqrt(np.sum((M[ent] - Q[best]) ** 2, axis=1)))]
        src[rcv] = best
        marked[ent[ent != rcv]] = True
    src[marked] = -1
    return src


def _drange(vol: np.ndarray) -> np.ndarray:
    v = frame_volume(vol)
    return np.array([[np.nanmin(v), np.nanmax(v)]])


def dmask(mask: np.ndarray, det) -> np.ndarray:
    """uint8 (z, y, x): 255 on the mask's components (26-connected) that hold a detection's
    starting voxel (1-based ``x_init ...`` of the MATLAB fields); the whole mask without
    detections"""
    from scipy import ndimage

    m = np.asarray(mask, dtype=bool)
    if det is None:
        return m.astype(np.uint8) * 255
    lab, _ = ndimage.label(m, structure=np.ones((3, 3, 3), dtype=bool))
    idx = [
        np.asarray(det[k], dtype=np.float64).ravel().astype(np.int64) - 1
        for k in ("z_init", "y_init", "x_init")
    ]
    hit = np.unique(lab[tuple(idx)])
    hit = hit[hit > 0]
    return (np.isin(lab, hit) & m).astype(np.uint8) * 255


def dmask_name(k: int, nframes: int) -> str:
    """``dmask_<k>.tif``, k zero-padded to the digits of the frame count (dmask_01 of 11)"""
    return f"dmask_{k:0{len(str(int(nframes)))}d}.tif"


def assemble_frame(
    det, mask, vols, sigma, fits: dict, z_aniso: float, master: int = 0
) -> tuple[dict, np.ndarray]:
    """one frame's ``frameInfo`` entry and dmask. ``det``: the master's detections (MATLAB
    fields, ``to_matlab_fields``) or None; ``vols``: each channel's DS frame (z, y, x);
    ``sigma`` nCh x 2 after ``runner_sigma``; ``fits``: {channel: (Ac fields, xyzAc fields)}"""
    sigma = np.array(sigma, dtype=np.float64, ndmin=2)
    nch = sigma.shape[0]
    if det is None:
        out = {k: np.zeros((0, 0)) for k in FIELD_ORDER}
        out["s"] = sigma.copy()
        out["dRange"] = [
            _drange(v) if ch == master else np.array([[float(np.min(v[0])), float(np.max(v[0]))]])
            for ch, v in enumerate(vols)
        ]
        return out, dmask(mask, None)
    n = np.asarray(det["x"]).size
    out: dict = {}
    choice = {ch: other_channel(det, ac, xyz, sigma[master, 0]) for ch, (ac, xyz) in fits.items()}
    for k in ROW_FIELDS:
        a = np.zeros((nch, n), dtype=bool if k in _LOGICAL else np.float64)
        a[master] = np.asarray(det[k]).ravel()
        for ch, (ac, xyz) in fits.items():
            src = choice[ch]
            row = np.asarray(ac[k]).ravel().copy()
            take = src >= 0
            row[take] = np.asarray(xyz[k]).ravel()[src[take]]
            a[ch] = row
        out[k] = a
    out["s"] = sigma.copy()
    for k in ("x_init", "y_init", "z_init"):
        out[k] = np.asarray(det[k], dtype=np.float64).reshape(1, n)
    isp = np.zeros((nch, n), dtype=bool)
    isp[master] = np.asarray(det["isPSF"]).ravel()
    for ch in fits:
        isp[ch] = ~out["hval_AD"][ch]
    out["isPSF"] = isp
    out["xCoord"] = np.stack([out["x"][master], out["x_pstd"][master]], axis=1)
    out["yCoord"] = np.stack([out["y"][master], out["y_pstd"][master]], axis=1)
    out["zCoord"] = np.stack([out["z"][master], out["z_pstd"][master]], axis=1) * z_aniso
    out["amp"] = np.stack([out["A"][master], out["A_pstd"][master]], axis=1)
    out["dRange"] = [_drange(v) for v in vols]
    return {k: out[k] for k in FIELD_ORDER}, dmask(mask, det)


def detect_frame(
    vols: list,
    k: int,
    sigma: np.ndarray,
    z_aniso: float,
    master: int = 0,
    residual_sigma_threshold: float = 3.0,
    num_threads: int | None = None,
) -> tuple[dict, np.ndarray]:
    """``frameInfo`` entry and dmask of frame ``k`` (1-based; the background's seed) from the
    channels' DS frames ``vols`` (as read, (z, y, x))"""
    from petakit5d.detection import detect as _d
    from petakit5d.detection import fit as _fit

    V = [frame_volume(v) for v in vols]
    sig = runner_sigma(sigma)
    bi = _d.estimate_background(V[master], k)
    ps, mask = _d.point_source_detection(
        V[master],
        sig[master],
        background_info=np.r_[bi, residual_sigma_threshold],
        residual_sigma_threshold=residual_sigma_threshold,
        num_threads=num_threads,
    )
    if np.asarray(ps["x"]).size == 0:
        return assemble_frame(None, mask, vols, sig, {}, z_aniso, master)
    P = np.stack([np.asarray(ps[c], dtype=np.float64) for c in ("x", "y", "z")], axis=1)
    fits = {}
    for ch in range(sig.shape[0]):
        if ch == master:
            continue
        nan = np.full(P.shape[0], np.nan)
        ac = _fit.fit_gaussians3d(V[ch], P, nan, sig[ch], nan, "Ac", num_threads=num_threads)
        xyz = _fit.fit_gaussians3d(
            V[ch], P, ac["A"], sig[ch], ac["c"], "xyzAc", num_threads=num_threads
        )
        fits[ch] = (to_matlab_fields(ac), to_matlab_fields(xyz))
    return assemble_frame(to_matlab_fields(ps), mask, vols, sig, fits, z_aniso, master)


def run_detection(
    data: dict,
    sigma,
    results_path: str | os.PathLike,
    *,
    frames=None,
    residual_sigma_threshold: float = 3.0,
    num_threads: int | None = None,
) -> list[dict]:
    """``XR_runDetection3D(data, 'Sigma', sigma, 'ResultsPath', ..., 'DetectionMethod', 'lowSNR',
    'Mode', 'xyzAc', 'FitMixtures', false, 'BackgroundCorrection', true)`` for one condition.

    ``data``: MATLAB's condition struct fields ``source`` (the master channel's folder),
    ``channels`` (the channels' folders), ``framePathsDS`` (per channel the DS frames' paths in
    frame order) and ``zAniso``; ``sigma`` nCh x 2 ``[sigmaXY sigmaZ]``. Writes
    ``<results_path>/Detection3D.mat`` (``frameInfo``, a 1 x nFrames struct array, MATLAB's
    values), ``Masks/dmask_<k>.tif`` and an empty ``Amira/Detection``; returns the entries.
    ``residual_sigma_threshold``: the camera's (``XR_get_camera_standard_deviation_threshold``
    gives 3 for the demo's camera; not ported)."""
    from petakit5d.io import matfile, read_tiff, write_tiff

    chans = [os.path.normpath(str(c)) for c in data["channels"]]
    src = os.path.normpath(str(data["source"]))
    master = chans.index(src) if src in chans else 0
    paths = data["framePathsDS"]
    nframes = len(paths[master])
    sel = range(1, nframes + 1) if frames is None else [int(f) for f in frames]
    os.makedirs(os.path.join(results_path, "Masks"), exist_ok=True)
    os.makedirs(os.path.join(results_path, "Amira", "Detection"), exist_ok=True)
    entries = []
    for k in sel:
        vols = [read_tiff(paths[ch][k - 1]) for ch in range(len(paths))]
        entry, dm = detect_frame(
            vols, k, sigma, float(data["zAniso"]), master, residual_sigma_threshold, num_threads
        )
        write_tiff(dm, os.path.join(results_path, "Masks", dmask_name(k, nframes)))
        entries.append(entry)
    rows = [{f: (e[f] if f != "dRange" else list(e[f])) for f in FIELD_ORDER} for e in entries]
    fi = matfile.StructArray((1, len(rows)), FIELD_ORDER, rows)
    matfile.write_mat(
        os.path.join(results_path, "Detection3D.mat"), {"frameInfo": fi}, version="7.3"
    )
    return entries
