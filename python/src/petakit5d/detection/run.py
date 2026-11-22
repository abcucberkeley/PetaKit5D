"""``XR_runDetection3D``: detection of every frame of a condition, the other channels read out at
the master channel's detections, and the results in MATLAB's files (``Detection3D.mat`` with the
struct array ``frameInfo``, ``Masks/dmask_<k>.tif``, an empty ``Amira/Detection`` folder).

Measured as a black box from traced MATLAB runs and designed-input runs (task C02; no GPL
source read): bit for bit on 14 real frames when given MATLAB's detector and fit outputs.
``frameInfo`` keeps MATLAB's values (1-based positions), as the file holds them.

Per frame k (1-based) of the master channel (the channel whose folder is ``data['source']``):
the DS frame as double with zeros as NaN; ``estimate_background(vol, k)``; the lowSNR detector
with ``'BackgroundInfo' [bg, cam]`` (``cam`` the camera threshold, 3 for the demo's) and the
residual threshold of the frame's noise and ``cam`` (``residual_threshold``: the camera
threshold for the demo's frames, the noise for a noisier frame, task C02d);
then, for each other channel, ``'Ac'`` fits at the master's positions and ``'xyzAc'`` fits from
them (default options: windows ``ceil(2 sigma)``, ``kLevel = norminv(0.975)``), and per entry
the choice between them (``other_channel``).

``run_detection_frame_parallel``: ``XR_runDetection3D_frame_parallel`` (``XR_runDetTrack3d``'s
default detection; task C02c, black box): the same frames as tasks of ``petakit5d.jobs``, one
``Detections/Detection_F_<k>.mat`` (variable ``frameInfo_k``) per frame, gathered into
``Detection3D.mat``; ``Overwrite`` and frames without detections as ``XR_runDetection3D`` has
them (MATLAB's frame-parallel runner ignores the one and repeats the previous frame for the
other: docs/parity/detection.md D11, D12).
"""

from __future__ import annotations

import math
import os

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "FIELD_ORDER",
    "CameraSettingsError",
    "assemble_frame",
    "camera_threshold",
    "camera_type",
    "detections_file",
    "dmask_name",
    "frame_task",
    "frame_volume",
    "residual_threshold",
    "run_detection",
    "run_detection_frame_parallel",
    "runner_camera_threshold",
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


class CameraSettingsError(RuntimeError):
    """MATLAB's error while finding the camera type (``identifier`` holds its id)"""

    def __init__(self, identifier: str, message: str):
        super().__init__(message)
        self.identifier = identifier


def _os_listdir(folder: str) -> list[str] | None:
    """the system's names in a folder (never '.' or '..'); None for anything else"""
    if not os.path.isdir(folder):
        return None
    try:
        return os.listdir(folder)
    except OSError:
        return None


def _matlab_dir_names(folder: str, pc: bool, listdir=_os_listdir) -> list[str]:
    """the names MATLAB's ``dir(folder)`` lists, in its order (UTF-16 code units, as
    ``core.paths.matlab_dir``; measured on Linux: '#name' comes before '.'): a folder's names
    with '.' and '..'; nothing for '' or a name that is not a folder. On Windows (deduced: no
    MATLAB runs there) a drive's root lists no '.' and '..', which its file systems do not have,
    and a bare drive 'C:' lists the drive's root (MATLAB keeps one current folder, not one per
    drive)."""
    import re

    from petakit5d.core.paths import matlab_sort_key

    if folder == "":
        return []
    root = pc and re.fullmatch(r"[A-Za-z]:[\\/]*", folder) is not None
    names = listdir(folder[:2] + "\\" if root else folder)
    if names is None:
        return []
    if not root:
        names = [*names, ".", ".."]
    return sorted(names, key=matlab_sort_key)


def _settings_file(source, *, pc: bool | None = None, listdir=_os_listdir) -> str:
    """the file ``get_camera_type_from_setting_file`` (the lab's code) reads: ``data.source``
    split at ``filesep`` (its ``strsplit``) without a trailing empty part, then from the whole
    path to its first part (joined with ``filesep``) the first folder whose ``dir`` listing, its
    first two names dropped as '.' and '..', has a name ending in 'settings.txt' (any case):
    ``[folder filesep name]`` with the first such name. ``pc``: Windows' rules (default: this
    platform's; ``core.paths.filesep``), ``listdir``: the system's listing (tests give their own).

    MATLAB's errors: a folder on the way that lists fewer than two names raises
    ``MATLAB:subsdeldimmismatch`` (the deletion of the first two), a walk without a settings
    file ``MATLAB:UndefinedFunction`` (``setting_filename`` is never set). An absolute path's
    walk ends on Linux at its first part '', which lists nothing (``MATLAB:subsdeldimmismatch``,
    measured); on Windows at the drive 'C:', which lists the drive's root (so
    ``MATLAB:UndefinedFunction``), and for a UNC path at '\\\\server', no folder
    (``MATLAB:subsdeldimmismatch``; both deduced).

    Windows, where the port differs (docs/parity/detection.md D9): MATLAB splits at '\\' only,
    so it walks a source written with '/' (which its ``dir`` and ``fileparts`` take there) no
    higher than the last '\\'; the port splits at both separators. Repeated separators give
    empty parts here (MATLAB's ``strsplit`` merges them by default: the same folders, but a
    UNC path's '\\\\server' would become '\\server')."""
    import re

    from petakit5d.core import paths

    pc = paths.filesep == "\\" if pc is None else bool(pc)
    sep = "\\" if pc else "/"
    text = os.fspath(source) if isinstance(source, os.PathLike) else str(source)
    parts = re.split(r"[\\/]" if pc else "/", text)
    if parts[-1] == "":
        parts = parts[:-1]
    for i in range(len(parts)):
        cur = sep.join(parts[: len(parts) - i])
        listing = _matlab_dir_names(cur, pc, listdir)
        if len(listing) < 2:
            raise CameraSettingsError(
                "MATLAB:subsdeldimmismatch", "Matrix index is out of range for deletion."
            )
        hits = [nm for nm in listing[2:] if re.search(r"settings\.txt$", nm, flags=re.IGNORECASE)]
        if hits:
            return cur + sep + hits[0]
    raise CameraSettingsError(
        "MATLAB:UndefinedFunction", "Unrecognized function or variable 'setting_filename'."
    )


def camera_type(
    source: str | os.PathLike,
    *,
    names: tuple[str, ...] = ("Orca", "Andor"),
    classes: tuple[str, ...] = ("SCMOS", "EMCCD"),
) -> str:
    """``get_camera_type_from_setting_file(data)`` (the lab's code): the settings file of the
    channel folder ``data.source`` or of the first folder above it that has one
    (``_settings_file``: MATLAB's walk per platform, its errors kept); its first line starting
    with 'Camera Type' (any case, leading white space skipped); the class of the last of
    ``names`` the line contains (any case), else ''. A settings file without a 'Camera Type' line
    raises ``MATLAB:TooManyOutputsDueToMissingBraces``."""
    import re

    found = _settings_file(source)
    try:
        with open(found, "rb") as fh:
            text = fh.read().decode("latin-1")
    except OSError:
        # fopen gives -1 and textscan fails (a listed name that cannot be read; deduced)
        raise CameraSettingsError(
            "MATLAB:FileIO:InvalidFid",
            "Invalid file identifier. Use fopen to generate a valid file identifier.",
        ) from None
    lines = [ln.lstrip(" \t\b") for ln in text.split("\n")]
    match = [ln for ln in lines if re.match(r"Camera Type", ln, flags=re.IGNORECASE)]
    if not match:
        raise CameraSettingsError(
            "MATLAB:TooManyOutputsDueToMissingBraces",
            "Insufficient number of outputs from right hand side of equal sign to satisfy "
            "assignment.",
        )
    cls = ""
    for nm, c in zip(names, classes, strict=True):
        if nm.lower() in match[0].lower():
            cls = c
    return cls


def camera_threshold(
    source: str | os.PathLike, *, scmos: float = 3.0, emccd: float = 65.0
) -> float:
    """``XR_get_camera_standard_deviation_threshold(data)`` (the lab's code): the detector's
    residual sigma threshold from the camera class, ``SCMOS_DC_Thresh`` 3 for 'SCMOS',
    ``EMCCD_DC_Thres`` 65 for 'EMCCD', NaN when the class is unknown"""
    from petakit5d.core.dtypes import NAN

    cls = camera_type(source)
    if cls == "SCMOS":
        return float(scmos)
    if cls == "EMCCD":
        return float(emccd)
    return float(NAN)


def runner_camera_threshold(source: str | os.PathLike) -> float:
    """the camera threshold as the runners take it: ``camera_threshold(source)``, or NaN when
    the camera search stops with one of MATLAB's errors (no settings file above ``source``, or
    one without a camera line): ``XR_runDetection3D`` goes on without a camera threshold then
    (task C02d, probe r8d), with ``residual_threshold``'s noise. MATLAB's frame-parallel runner
    stops with the error there instead, its bug of D13, which the port does not reproduce: the
    lab's camera search means to give no camera class in both cases (its code returns '' for
    them, but fails before it), and the lab's threshold function to estimate the threshold from
    the image without camera information (task C09, probe d13)."""
    from petakit5d.core.dtypes import NAN

    try:
        return camera_threshold(source)
    except CameraSettingsError:
        return float(NAN)


def residual_threshold(noise: float, camera: float) -> float:
    """the residual sigma threshold the runners give the lowSNR detector for a frame
    (``'ResidualSigmaThresh'``, which sets its candidates' filter and the bad-point filter's
    ``thr``): from the frame's background noise (the third value of
    ``XR_estimate_overall_background_information``) and the camera threshold
    (``camera_threshold``). Measured as a black box (task C02d, probes r8a to r8c: the runners'
    detector calls traced on frames of noise 1.5 to 111 with camera thresholds 3, 65 and none):
    0.9 camera thresholds while the noise is below 1/1.5 of it (bounds measured between 0.662
    and 0.682 of it), the camera threshold up to it, the noise above it up to 1.5 camera
    thresholds (between 1.477 and 1.506), the camera threshold beyond. Without a camera
    threshold (NaN: an unknown camera, or none found, ``runner_camera_threshold``) the noise,
    as ``XR_runDetection3D`` has it since the fix of its dropped detections without one;
    MATLAB's frame-parallel runner passes NaN there and so finds nothing, or stops with the
    camera search's error, its bug, which the port does not reproduce (docs/parity/detection.md
    D13). A camera threshold of 0 is 0 for every frame, as both of MATLAB's runners pass it on
    (task C09, probe d13: a stubbed threshold function). The values at exactly the bounds are
    not measured."""
    noise, camera = float(noise), float(camera)
    if math.isnan(camera):
        return noise
    if noise < camera / 1.5:
        return 0.9 * camera
    if noise > 1.5 * camera:
        return camera
    # MATLAB's max, which leaves NaN out
    return camera if math.isnan(noise) else max(noise, camera)


def frame_volume(raw) -> np.ndarray:
    """a DS frame as the runner uses it: double, zeros as NaN (MATLAB's)"""
    from petakit5d.core.dtypes import NAN

    v = np.asarray(raw).astype(np.float64)
    v[v == 0] = NAN
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
    ``sigma`` after ``runner_sigma``, a row per channel (``s`` holds it whole, rows past the
    channels too, as MATLAB's does); ``fits``: {channel: (Ac fields, xyzAc fields)}"""
    sigma = np.array(sigma, dtype=np.float64, ndmin=2)
    nch = len(vols)
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
    camera_thr: float = 3.0,
    num_threads: int | None = None,
    *,
    detection_method: str = "lowSNR",
    background_correction: bool = True,
    window=None,
) -> tuple[dict, np.ndarray]:
    """``frameInfo`` entry and dmask of frame ``k`` (1-based; the background's seed) from the
    channels' DS frames ``vols`` (as read, (z, y, x)); ``detection_method`` 'lowSNR'
    (``XR_pointSourceDetection3D_3`` after the background estimate, with the residual threshold
    of the frame's noise and the camera threshold ``camera_thr`` (``residual_threshold``) and
    ``'BackgroundInfo'`` the estimate and ``camera_thr``) or 'speedup'
    (``XR_pointSourceDetection3D``, redundancy radius 0.5, no background estimate);
    ``window``: the detector's ``WindowSize`` (empty: from sigma)"""
    from petakit5d.detection import detect as _d
    from petakit5d.detection import fit as _fit

    V = [frame_volume(v) for v in vols]
    sig = runner_sigma(sigma)
    if detection_method == "speedup":
        ps, mask = _d.point_source_detection_speedup(
            V[master],
            sig[master],
            mode="xyzAc",
            alpha=0.05,
            remove_redundant_fits=True,
            redundancy_radius=0.5,
            window=window,
            background_correction=background_correction,
            num_threads=num_threads,
        )
    elif detection_method == "lowSNR":
        bi = _d.estimate_background(V[master], k)
        ps, mask = _d.point_source_detection(
            V[master],
            sig[master],
            window=window,
            background_correction=background_correction,
            background_info=np.r_[bi, camera_thr],
            residual_sigma_threshold=residual_threshold(bi[2], camera_thr),
            num_threads=num_threads,
        )
    else:
        raise NotImplementedError(
            f"DetectionMethod {detection_method!r}: 'lowSNR' and 'speedup' are ported"
        )
    if np.asarray(ps["x"]).size == 0:
        return assemble_frame(None, mask, vols, sig, {}, z_aniso, master)
    P = np.stack([np.asarray(ps[c], dtype=np.float64) for c in ("x", "y", "z")], axis=1)
    fits = {}
    for ch in range(len(V)):
        if ch == master:
            continue
        nan = np.full(P.shape[0], np.nan)
        ac = _fit.fit_gaussians3d(V[ch], P, nan, sig[ch], nan, "Ac", num_threads=num_threads)
        xyz = _fit.fit_gaussians3d(
            V[ch], P, ac["A"], sig[ch], ac["c"], "xyzAc", num_threads=num_threads
        )
        fits[ch] = (to_matlab_fields(ac), to_matlab_fields(xyz))
    return assemble_frame(to_matlab_fields(ps), mask, vols, sig, fits, z_aniso, master)


def _strings(v) -> list[str]:
    """a path or MATLAB's cell of paths (any shape, column-major) as a list of str"""
    if isinstance(v, str | os.PathLike):
        return [os.fspath(v)]
    if isinstance(v, np.ndarray):
        return [os.fspath(p) if isinstance(p, os.PathLike) else str(p) for p in v.ravel(order="F")]
    return [os.fspath(p) if isinstance(p, os.PathLike) else str(p) for p in v]


def _frame_paths(v) -> list[list[str]]:
    """``framePathsDS``: per channel its frames' paths (MATLAB's 1 x nCh cell of nFrames x 1
    cells, ``condition.Movie``'s list of columns, or one channel's paths)"""
    items = list(v.ravel(order="F")) if isinstance(v, np.ndarray) and v.dtype == object else list(v)
    if items and all(isinstance(p, str | os.PathLike) for p in items):
        return [_strings(items)]
    return [_strings(c) for c in items]


def _movies(data) -> list[dict]:
    """the condition's movies as MATLAB's structs: one movie's dict, a ``StructArray`` (its
    elements in MATLAB's linear order), or a sequence of dicts or ``condition.Movie``s"""
    from collections.abc import Mapping

    from petakit5d.io.matfile import StructArray

    if isinstance(data, StructArray):
        return data.ravel()
    if isinstance(data, list | tuple):  # a condition.Condition too (registered as a Mapping)
        return [m.to_matlab() if hasattr(m, "to_matlab") else dict(m) for m in data]
    if isinstance(data, Mapping):
        return [dict(data)]
    if hasattr(data, "to_matlab"):
        return [data.to_matlab()]
    raise TypeError(f"run_detection: data is a struct or struct array of movies, not {type(data)}")


def _results_paths(results_path, movies: list[dict]) -> list[str]:
    """each movie's results folder: by default ``[source 'Analysis' filesep]``, else a cell
    (list) of one folder per movie"""
    if results_path is None:
        return [str(m["source"]) + "Analysis" + os.sep for m in movies]
    if isinstance(results_path, str | os.PathLike):
        raise TypeError(
            "run_detection: results_path is a list of folders, one per movie (MATLAB refuses a "
            "char 'ResultsPath': it must be a cell)"
        )
    paths = _strings(results_path)
    if len(paths) < len(movies):
        raise IndexError(
            f"Index exceeds the number of array elements. Index must not exceed {len(paths)}."
        )
    return paths[: len(movies)]


def _is_default(value, default) -> bool:
    """Whether a MATLAB-side value equals the default ([] for None; a logical or number by
    value; text by text)."""
    if default is None:
        return value is None or np.size(value) == 0
    if isinstance(default, str):
        return isinstance(value, str) and value == default
    a = np.asarray(value).reshape(-1)
    return a.size == 1 and float(a[0]) == float(default)


# the two runners' parameters the port takes at MATLAB's default only
_DEFAULTS_ONLY = (
    ("RemoveRedundant", True),
    ("RedundancyRadius", 0.5),
    ("RemoveSecondaryRedundant", True),
    ("Master", None),
    ("Alpha", 0.05),
    ("CellMask", None),
    ("FileName", "Detection3D.mat"),
    ("FilterByResidualSigma", True),
)


def _defaults_only(runner: str, values) -> None:
    """NotImplementedError for a parameter of ``_DEFAULTS_ONLY`` (``values`` in its order) off
    its default"""
    for (name, default), value in zip(_DEFAULTS_ONLY, values, strict=True):
        if not _is_default(value, default):
            raise NotImplementedError(f"{runner}: {name} {value!r}: the default is ported")


@register(
    "XR_runDetection3D",
    Param("data", None, "Movies (the structure XR_loadConditionData3D returns).", convert="same"),
    Param("Sigma", None, "PSF sigmas, one row [sigma_xy, sigma_z] per channel.", convert="same"),
    Param("Mode", None, "Fitted parameters ('xyzAc' only)."),
    Param("FitMixtures", None, "Fit mixtures (false only)."),
    Param("MaxMixtures", None, "Most components of a mixture fit."),
    Param("RemoveRedundant", None, "Remove redundant detections (true only)."),
    Param("RedundancyRadius", None, "Radius of redundant detections (the default only)."),
    Param("RemoveSecondaryRedundant", None, "The same in the other channels (true only)."),
    Param("Overwrite", None, "Detect movies whose result exists again."),
    Param("Master", None, "Master channel ([] only: the channel of data.source)."),
    Param("Alpha", None, "Significance level of the detection tests (the default only)."),
    Param("CellMask", None, "Cell mask ([] only)."),
    Param("FileName", None, "Result file name (the default only)."),
    Param("WindowSize", None, "Detection window size ([] only).", convert="same"),
    Param(
        "ResultsPath",
        None,
        "Results folder of each movie (default <source>/Analysis/).",
        matlab_default="arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0)",
    ),
    Param("DetectionMethod", None, "'speedup' or 'lowSNR'."),
    Param("WriteAmira", None, "Amira files (MATLAB's writer fails: an empty folder stays)."),
    Param("BackgroundCorrection", None, "Correct the background before the detection."),
    Param("FilterByResidualSigma", None, "Filter detections by their residuals (true only)."),
    Param("Frames", None, "Only these frames (1-based).", python_only=True),
    Param(
        "ResidualSigmaThreshold",
        None,
        "The camera's residual threshold (None: from the settings file above source).",
        python_only=True,
    ),
    Param("numThreads", None, "Threads (None: the CPUs available).", python_only=True),
    task="C02",
)
def run_detection(
    data,
    *,
    sigma=None,
    mode: str = "xyzAc",
    fit_mixtures: bool = False,
    max_mixtures: int = 3,
    remove_redundant: bool = True,
    redundancy_radius: float = 0.5,
    remove_secondary_redundant: bool = True,
    overwrite: bool = False,
    master=None,
    alpha: float = 0.05,
    cell_mask=None,
    file_name: str = "Detection3D.mat",
    window_size=None,
    results_path=None,
    detection_method: str = "speedup",
    write_amira: bool = True,
    background_correction: bool = False,
    filter_by_residual_sigma: bool = True,
    frames=None,
    residual_sigma_threshold: float | None = None,
    num_threads: int | None = None,
) -> list[list[dict] | None]:
    """``XR_runDetection3D(data, 'Sigma', sigma, 'Overwrite', ow, 'WindowSize', ws,
    'ResultsPath', rp, 'Mode', mode, 'FitMixtures', fm, 'MaxMixtures', mm, 'DetectionMethod',
    method, 'BackgroundCorrection', bc, 'WriteAmira', wa)``: every frame of every movie of a
    condition. Measured as a black box: the defaults, the folders, the skipping.

    ``data``: the condition as ``condition.load_condition_data`` gives it, one movie's dict or a
    1 x n ``StructArray`` (or a list of movies); per movie MATLAB's fields ``source`` (the master
    channel's folder), ``channels``, ``framePathsDS`` and ``zAniso``. ``sigma``: nCh x 2
    ``[sigmaXY sigmaZ]``, a row per channel (``frameInfo.s`` holds it whole, rows past the
    channels too). ``results_path``: a list of one folder per movie (MATLAB refuses a char
    array), by default ``[source 'Analysis' filesep]``. As MATLAB does, ``Detection3D.mat``
    goes to ``results_path + 'Detection3D.mat'`` (a folder without the trailing separator gives
    a file beside it) and the masks to ``<results_path>/Masks/dmask_<k>.tif``; with
    ``overwrite`` false (the default) a movie whose ``Detection3D.mat`` exists is skipped (None
    in the result). ``detection_method``: 'speedup' (MATLAB's default) or 'lowSNR'; the ported
    settings are ``window_size`` empty, ``mode`` 'xyzAc' and ``fit_mixtures`` false
    (``max_mixtures`` matters only with mixtures). ``write_amira``: MATLAB's Amira writer fails
    and writes nothing; an empty ``Amira/Detection`` folder stays. ``frames``: only these frames
    (1-based; the port's); ``residual_sigma_threshold``: the camera threshold, by default as
    MATLAB finds it from the settings file above ``source`` (``camera_threshold``: 3 for sCMOS
    cameras; none, NaN, without such a file: ``runner_camera_threshold``); each frame's
    detector gets the residual threshold of it and the frame's noise (``residual_threshold``).
    Returns per movie its ``frameInfo`` entries."""
    if sigma is None or np.size(sigma) == 0:
        raise ValueError("run_detection: give sigma (MATLAB asks for it interactively)")
    if mode != "xyzAc":
        raise NotImplementedError(f"run_detection: Mode {mode!r}: 'xyzAc' is ported")
    if bool(np.asarray(fit_mixtures).reshape(-1)[0]):
        raise NotImplementedError("run_detection: FitMixtures true is not ported")
    if detection_method not in ("speedup", "lowSNR"):
        raise NotImplementedError(
            f"run_detection: DetectionMethod {detection_method!r}: 'speedup' and 'lowSNR' are "
            "ported"
        )
    _defaults_only(
        "run_detection",
        (remove_redundant, redundancy_radius, remove_secondary_redundant, master, alpha,
         cell_mask, file_name, filter_by_residual_sigma),
    )  # fmt: skip
    if window_size is not None and np.size(window_size) != 0:
        # MATLAB's detectors take 'WindowSize' otherwise than the port's window (measured: 4,830
        # detections against 3,128 for [2 2] on a demo frame)
        raise NotImplementedError("run_detection: WindowSize: the default (empty) is ported")
    movies = _movies(data)
    paths = _results_paths(results_path, movies)
    out: list[list[dict] | None] = []
    for mv, rp in zip(movies, paths, strict=True):
        if not overwrite and os.path.isfile(rp + "Detection3D.mat"):
            out.append(None)
            continue
        out.append(
            _run_movie(
                mv,
                sigma,
                rp,
                window_size=window_size,
                detection_method=detection_method,
                background_correction=bool(background_correction),
                write_amira=write_amira,
                frames=frames,
                residual_sigma_threshold=residual_sigma_threshold,
                num_threads=num_threads,
            )
        )
    return out


def _run_movie(
    data: dict,
    sigma,
    results_path: str,
    *,
    window_size,
    detection_method: str,
    background_correction: bool,
    write_amira: bool,
    frames,
    residual_sigma_threshold,
    num_threads,
) -> list[dict]:
    """one movie of ``run_detection``"""
    from petakit5d.io import read_tiff, write_tiff

    chans = [os.path.normpath(c) for c in _strings(data["channels"])]
    src = os.path.normpath(str(data["source"]))
    master = chans.index(src) if src in chans else 0
    paths = _frame_paths(data["framePathsDS"])
    nframes = len(paths[master])
    sel = range(1, nframes + 1) if frames is None else [int(f) for f in frames]
    thr = residual_sigma_threshold
    if thr is None:
        thr = runner_camera_threshold(data["source"])
    masks = results_path + os.sep + "Masks"
    os.makedirs(masks, exist_ok=True)
    if write_amira:
        os.makedirs(os.path.join(results_path, "Amira", "Detection"), exist_ok=True)
    entries = []
    for k in sel:
        vols = [read_tiff(paths[ch][k - 1]) for ch in range(len(paths))]
        entry, dm = detect_frame(
            vols,
            k,
            sigma,
            float(np.asarray(data["zAniso"]).reshape(-1)[0]),
            master,
            float(thr),
            num_threads,
            detection_method=detection_method,
            background_correction=background_correction,
            window=window_size,
        )
        write_tiff(dm, os.path.join(masks, dmask_name(k, nframes)))
        entries.append(entry)
    _write_detection3d(results_path + "Detection3D.mat", entries)
    return entries


def _entry_row(e: dict) -> dict:
    """a frameInfo entry as write_mat takes it (MATLAB's field order, dRange a cell row)"""
    return {f: (e[f] if f != "dRange" else list(e[f])) for f in FIELD_ORDER}


def _write_detection3d(path: str, entries: list[dict]) -> None:
    from petakit5d.io import matfile

    rows = [_entry_row(e) for e in entries]
    fi = matfile.StructArray((1, len(rows)), FIELD_ORDER, rows)
    matfile.write_mat(path, {"frameInfo": fi}, version="7.3")


# -- XR_runDetection3D_frame_parallel ------------------------------------------------------------


def detections_file(results_path: str, k: int) -> str:
    """the frame-parallel runner's per-frame file ``[ResultsPath 'Detections' filesep
    'Detection_F_<k>.mat']`` (k 1-based, not padded; a results folder without the trailing
    separator gives ``...AnalysisDetections``, as MATLAB's)"""
    return f"{results_path}Detections{os.sep}Detection_F_{int(k)}.mat"


def frame_task(
    frame_paths: list[str],
    k: int,
    nframes: int,
    sigma,
    z_aniso: float,
    master: int,
    threshold: float,
    results_path: str,
    *,
    detection_method: str = "lowSNR",
    background_correction: bool = False,
    window_size=None,
    num_threads: int | None = None,
) -> None:
    """One frame of ``run_detection_frame_parallel`` (a task of the job framework): frame ``k``
    (1-based) of every channel (``frame_paths``, the DS frames in channel order) detected as
    ``run_detection`` does it (``detect_frame``), its dmask written to
    ``<results_path>/Masks/dmask_<k>.tif``, then its entry published as
    ``detections_file(results_path, k)`` (MAT v7.3, the variable ``frameInfo_k``), last and by
    one rename, since the file's existence marks the task done."""
    from petakit5d.io import matfile, read_tiff, write_tiff
    from petakit5d.jobs.publish import publish_result, temp_path

    vols = [read_tiff(p) for p in frame_paths]
    entry, dm = detect_frame(
        vols,
        int(k),
        np.asarray(sigma, dtype=np.float64),
        float(z_aniso),
        int(master),
        float(threshold),
        num_threads,
        detection_method=detection_method,
        background_correction=bool(background_correction),
        window=window_size,
    )
    write_tiff(dm, os.path.join(results_path + os.sep + "Masks", dmask_name(int(k), int(nframes))))
    out = detections_file(results_path, int(k))
    tmp = temp_path(out)
    matfile.write_mat(tmp, {"frameInfo_k": _entry_row(entry)}, version="7.3")
    publish_result(tmp, out)


def _read_frame_entry(path: str) -> dict:
    from petakit5d.io import matfile

    e = matfile.read_mat(path, "frameInfo_k")["frameInfo_k"]
    return {f: (list(e[f]) if f == "dRange" else e[f]) for f in FIELD_ORDER}


@register(
    "XR_runDetection3D_frame_parallel",
    Param("data", None, "Movies (the structure XR_loadConditionData3D returns).", convert="same"),
    Param("Sigma", None, "PSF sigmas, one row [sigma_xy, sigma_z] per channel.", convert="same"),
    Param("Mode", None, "Fitted parameters ('xyzAc' only)."),
    Param("FitMixtures", None, "Fit mixtures (false only)."),
    Param("MaxMixtures", None, "Most components of a mixture fit."),
    Param("RemoveRedundant", None, "Remove redundant detections (true only)."),
    Param("RedundancyRadius", None, "Radius of redundant detections (the default only)."),
    Param("RemoveSecondaryRedundant", None, "The same in the other channels (true only)."),
    Param("Overwrite", None, "Detect movies whose result exists again, and their frames."),
    Param("Master", None, "Master channel ([] only: the channel of data.source)."),
    Param("Alpha", None, "Significance level of the detection tests (the default only)."),
    Param("CellMask", None, "Cell mask ([] only)."),
    Param("FileName", None, "Result file name (the default only)."),
    Param("WindowSize", None, "Detection window size ([] only).", convert="same"),
    Param(
        "ResultsPath",
        None,
        "Results folder of each movie (default <source>/Analysis/).",
        matlab_default="arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0)",
    ),
    Param(
        "DetectionMethod",
        None,
        "'speedup' or 'lowSNR' (the default 'update_1' runs XR_pointSourceDetection3D_1: not "
        "ported).",
    ),
    Param("WriteAmira", None, "Amira files (MATLAB's writer fails: an empty folder stays)."),
    Param("BackgroundCorrection", None, "Correct the background before the detection."),
    Param("FilterByResidualSigma", None, "Filter detections by their residuals (true only)."),
    Param("Frames", None, "Only these frames (1-based).", python_only=True),
    Param(
        "ResidualSigmaThreshold",
        None,
        "The camera's residual threshold (None: from the settings file above source).",
        python_only=True,
    ),
    Param("numThreads", None, "Threads of a frame (None: the CPUs available).", python_only=True),
    Param(
        "parseCluster",
        None,
        "Run the frames as Slurm jobs (false: here, as MATLAB's runner computes them).",
        python_only=True,
    ),
    Param("localWorkerNum", None, "Processes of a local pool (1: this one).", python_only=True),
    Param("jobLogDir", None, "Job log folder of the job framework.", python_only=True),
    Param("cpusPerTask", None, "CPUs per cluster task.", python_only=True),
    Param(
        "memAllocate", None, "Memory per cluster task (None: the framework's).", python_only=True
    ),
    Param("configFile", None, "The job framework's configuration file.", python_only=True),
    task="C02c",
)
def run_detection_frame_parallel(
    data,
    *,
    sigma=None,
    mode: str = "xyzAc",
    fit_mixtures: bool = False,
    max_mixtures: int = 3,
    remove_redundant: bool = True,
    redundancy_radius: float = 0.5,
    remove_secondary_redundant: bool = True,
    overwrite: bool = False,
    master=None,
    alpha: float = 0.05,
    cell_mask=None,
    file_name: str = "Detection3D.mat",
    window_size=None,
    results_path=None,
    detection_method: str = "update_1",
    write_amira: bool = True,
    background_correction: bool = False,
    filter_by_residual_sigma: bool = True,
    frames=None,
    residual_sigma_threshold: float | None = None,
    num_threads: int | None = None,
    parse_cluster: bool = False,
    local_worker_num: int = 1,
    job_log_dir: str = "",
    cpus_per_task: int = 1,
    mem_allocate: float | None = None,
    config_file: str = "",
) -> list[list[dict] | None]:
    """``XR_runDetection3D_frame_parallel(data, 'Sigma', sigma, 'Overwrite', ow, 'WindowSize',
    ws, 'ResultsPath', rp, 'Mode', mode, 'FitMixtures', fm, 'MaxMixtures', mm,
    'DetectionMethod', method, 'BackgroundCorrection', bc)``, ``XR_runDetTrack3d``'s detection
    stage: the frames of every movie as tasks of the job framework (``petakit5d.jobs``), then
    their entries gathered into ``Detection3D.mat``. Measured as a black box (probes of MATLAB
    R2024b); the same detections and files as ``run_detection`` with the same options.

    Per movie (results folder as ``run_detection``'s): with ``overwrite`` false a movie whose
    ``Detection3D.mat`` exists is skipped (None in the result). Each frame ``k`` is a task
    (``frame_task``) whose output is ``<results_path>Detections/Detection_F_<k>.mat`` (one
    variable ``frameInfo_k``, the frame's ``frameInfo`` entry); a frame whose file exists is not
    run again, and a task another caller holds (its ``Detection_F_<k>.tmp`` marker, fresh) is
    left to it, as MATLAB's job framework does, so several callers share a movie's frames.
    When every frame's file exists, ``Detection3D.mat`` is written from them (else the movie
    gives None); ``Detection3D.mat`` existing counts as all frames done.

    Three differences from MATLAB (its quirks that are bugs; docs/parity/detection.md D11 to
    D13): ``overwrite`` true removes the movie's ``Detection3D.mat`` and per-frame files first,
    so the frames are detected again (MATLAB keeps them and so ignores the flag); a frame
    without detections gets the empty entry, as ``run_detection`` gives it (MATLAB's
    frame-parallel runner repeats the previous frame's detections there); and without a camera
    threshold (an unknown camera, or no settings file with a camera line above ``source``) the
    frames are detected with their noise as the residual threshold, as ``run_detection`` does
    (MATLAB's frame-parallel runner passes NaN and finds nothing, or stops with the camera
    search's error).

    ``detection_method`` 'lowSNR' (what ``XR_runDetTrack3d`` passes) or 'speedup'; MATLAB's
    default 'update_1' runs ``XR_pointSourceDetection3D_1``, which is not ported. The other
    options as ``run_detection``'s (its parameters, in MATLAB's order; ``sigma`` and the rest
    by name). The job framework's (the port's): ``parse_cluster`` (MATLAB's runner computes in
    its session: False), ``local_worker_num`` (a pool of that many processes), ``job_log_dir``,
    ``cpus_per_task``, ``mem_allocate``, ``config_file`` (``run_tasks``'s). In one process the
    detector's t threshold is computed once (MATLAB's per-session quirk); a pool's or cluster's
    processes each compute their own, as MATLAB's workers do. Returns per movie its
    ``frameInfo`` entries."""
    if detection_method == "update_1":
        raise NotImplementedError(
            "run_detection_frame_parallel: MATLAB's default DetectionMethod 'update_1' runs "
            "XR_pointSourceDetection3D_1, which is not ported; give 'lowSNR' or 'speedup'"
        )
    if sigma is None or np.size(sigma) == 0:
        raise ValueError(
            "run_detection_frame_parallel: give sigma (MATLAB asks for it interactively)"
        )
    if mode != "xyzAc":
        raise NotImplementedError(f"run_detection_frame_parallel: Mode {mode!r}: 'xyzAc' is ported")
    if bool(np.asarray(fit_mixtures).reshape(-1)[0]):
        raise NotImplementedError("run_detection_frame_parallel: FitMixtures true is not ported")
    if detection_method not in ("speedup", "lowSNR"):
        raise NotImplementedError(
            f"run_detection_frame_parallel: DetectionMethod {detection_method!r}: 'speedup' and "
            "'lowSNR' are ported"
        )
    _defaults_only(
        "run_detection_frame_parallel",
        (remove_redundant, redundancy_radius, remove_secondary_redundant, master, alpha,
         cell_mask, file_name, filter_by_residual_sigma),
    )  # fmt: skip
    if window_size is not None and np.size(window_size) != 0:
        raise NotImplementedError(
            "run_detection_frame_parallel: WindowSize: the default (empty) is ported"
        )
    del max_mixtures  # (only with mixtures)
    movies = _movies(data)
    paths = _results_paths(results_path, movies)
    jobs = {
        "parse_cluster": bool(parse_cluster),
        "local_worker_num": int(local_worker_num),
        "cpus_per_task": cpus_per_task,
        "mem_allocate": mem_allocate,
        "config_file": config_file,
    }
    if job_log_dir:
        jobs["job_log_dir"] = job_log_dir
    out: list[list[dict] | None] = []
    for mv, rp in zip(movies, paths, strict=True):
        final = rp + "Detection3D.mat"
        if os.path.isfile(final) and not overwrite:
            out.append(None)
            continue
        out.append(
            _run_movie_frame_parallel(
                mv,
                sigma,
                rp,
                overwrite=bool(overwrite),
                window_size=window_size,
                detection_method=detection_method,
                background_correction=bool(background_correction),
                write_amira=write_amira,
                frames=frames,
                residual_sigma_threshold=residual_sigma_threshold,
                num_threads=num_threads,
                jobs=jobs,
            )
        )
    return out


def _run_movie_frame_parallel(
    data: dict,
    sigma,
    results_path: str,
    *,
    overwrite: bool,
    window_size,
    detection_method: str,
    background_correction: bool,
    write_amira: bool,
    frames,
    residual_sigma_threshold,
    num_threads,
    jobs: dict,
) -> list[dict] | None:
    """one movie of ``run_detection_frame_parallel``"""
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task

    chans = [os.path.normpath(c) for c in _strings(data["channels"])]
    src = os.path.normpath(str(data["source"]))
    master = chans.index(src) if src in chans else 0
    fpaths = _frame_paths(data["framePathsDS"])
    nframes = len(fpaths[master])
    sel = list(range(1, nframes + 1)) if frames is None else [int(f) for f in frames]
    thr = residual_sigma_threshold
    if thr is None:
        # (MATLAB's frame-parallel runner stops on the camera search's errors and drops every
        # detection without a camera threshold: D13; the port follows XR_runDetection3D)
        thr = runner_camera_threshold(data["source"])
    final = results_path + "Detection3D.mat"
    outputs = [detections_file(results_path, k) for k in sel]
    if overwrite:
        for p in [final, *outputs]:
            if os.path.isfile(p):
                os.remove(p)
    os.makedirs(results_path + os.sep + "Masks", exist_ok=True)
    os.makedirs(os.path.dirname(outputs[0]) if outputs else results_path, exist_ok=True)
    s = np.array(sigma, dtype=np.float64, ndmin=2).tolist()
    z_aniso = float(np.asarray(data["zAniso"]).reshape(-1)[0])
    tasks = [
        Task.of(
            frame_task,
            [fpaths[ch][k - 1] for ch in range(len(fpaths))],
            k,
            nframes,
            s,
            z_aniso,
            master,
            float(thr),
            results_path,
            detection_method=detection_method,
            background_correction=background_correction,
            window_size=None,
            num_threads=num_threads,
        )
        for k in sel
    ]
    framework.run_tasks(
        [fpaths[master][k - 1] for k in sel], outputs, tasks, final_out_fullpath=final, **jobs
    )
    if write_amira:
        os.makedirs(os.path.join(results_path, "Amira", "Detection"), exist_ok=True)
    if not all(os.path.isfile(p) for p in outputs):
        return None
    entries = [_read_frame_entry(p) for p in outputs]
    if not os.path.isfile(final):
        _write_detection3d(final, entries)
    return entries
