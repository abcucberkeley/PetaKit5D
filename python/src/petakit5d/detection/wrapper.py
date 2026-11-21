"""The detection and tracking wrapper (the lab's ``XR_runDetTrack3d``, ported from its code) and its
deskew stage (``deskewData``, another author's GPL function: from its behaviour, measured by
black-box probes).

``run_detection_tracking(data, ...)`` runs, as MATLAB does:

1. the results folder of each movie, ``<source><aname>/`` (``aname`` None: ``Analysis``;
   with ``apath``: ``<apath><the last three folders of the first channel>/<aname>/``), made;
2. ``deskew_data`` (with the ``aname`` as given, so that with the default its ``settings.mat``
   and ``dynRange.mat`` land in the channel folder itself, not in the results folder), each
   movie deskewed in its scan direction (its ``reversed``). MATLAB passes deskewData no
   ``'Reverse'``, so it deskews reversed scans as forward ones, which shears them the wrong way:
   a slip the port fixes (``docs/parity/detection.md`` W2);
3. ``correct_xz_offset_data`` when ``x_zoffset_correction``;
4. detection: ``run_detection`` for ``parallel_method`` ``'parfor'`` (always on Windows),
   ``run_detection_frame_parallel`` for ``'lock_distributed'`` (the default), none for any other
   text;
5. with ``track``, on the movies of more than one frame: ``run_tracking`` with
   ``track_settings(radius=tracking_radius, max_gap_length=tracking_gap_length)``,
   ``run_track_processing``, and ``rotate_tracks`` when ``rotate_tracks``, each movie with its
   own results folder. MATLAB drops the one-frame movies (tracking needs two frames) but passes
   the folders of all movies, so after a one-frame movie each movie is paired with the folder of
   the one before it: it tracks the one-frame movie's detections into that movie's folder,
   leaves its own folder without tracks, and fails in track processing. An indexing slip the
   port fixes (W3).

``matlab_quirks=True`` gives MATLAB's output for both slips (reversed scans deskewed as forward
ones, the tracking stages given the folders of all movies), and the rotation stage's for a
reversed scan (rotateTracks3D rotates it by the forward angle, W20).

``overwrite`` is one flag for all stages or four (deskew and offsets, detection, tracking, track
processing and rotation). Unknown MATLAB options are dropped, as MATLAB's parser keeps them
unmatched (the multiple-experiments demo's ``'zoffsetCorrection'`` does nothing). MATLAB's
interactive branches (no ``data``, an empty ``sigma``) are errors here.

The stages of steps 4 and 5 are other tasks' (detection C02, tracking C03, track processing and
rotation C04b). They are looked up in ``STAGES`` by name and called as MATLAB calls them,
once each, with the condition as the port holds MATLAB's struct array
(``condition.struct_array``: one movie's dict, more movies a 1 x n ``StructArray``, MATLAB's
field names), ``results_path`` the list of the movies' results folders, and keyword arguments
named as MATLAB's options in snake_case (the docstring of ``_call`` lists them). A stage that is
not there yet raises ``NotImplementedError``. Tests put stand-ins there.

``deskew_data(data, ...)`` (MATLAB's ``deskewData`` as the wrapper and the by-step demo call
it): for each movie, unless channel 1's DS folder already holds as many ``*.tif*`` files as the
movie has frames and ``overwrite`` is off (then the movie is left as it is, without writing
anything), it writes ``<source><aname>/settings.mat`` (``dz``, ``xyPixelSize``,
``objectiveScan``, ``reversed`` of the movie), each frame of each channel deskewed
(``deskew_rotate.deskew_volume`` of the frame as double with ``skew_angle``, the movie's ``dz``
and pixel size and the scan direction; an objective scan's frames as they are) as single into
``<channel>DS/<name>``, and ``<source><aname>/dynRange.mat`` (``dRange``: frames x [min max] x
channels of the deskewed doubles). The direction: ``reverse`` None takes each
movie's ``reversed`` (what the wrapper passes); True or False is deskewData's ``'Reverse'``
option, one direction for all movies, False by default as in MATLAB (which XR_runDetTrack3d
leaves: reversed scans then come out sheared the wrong way, W2). Existing DS frames of a movie
that is not skipped are written again.
``data`` is returned unchanged. Bit for bit MATLAB's in the probes and the oracle cases: DS
frames of uint16 and single raw frames, angles 31.5 and 45, the movie's dz and pixel size,
reversed scans with ``'Reverse'`` true (the lab's deskew of a reversed scan: measured equal to
``XR_deskewRotateFrame`` with ``reverse`` true); ``dRange`` from the doubles. Not ported:
``s_cmos_camera_flip`` (MATLAB calls a missing ``GU_sCMOSCameraFlip``: the same error),
``crop`` (MATLAB waits for an interactive region), ``rotate`` (the DSR frames),
``llff_correction`` and ``zoffset_correction``.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = ["STAGES", "WrapperError", "deskew_data", "run_detection_tracking"]


class WrapperError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# the stages of other tasks, by name: (module, function) until a stand-in or the port is set
_STAGE_HOMES = {
    "run_detection": ("petakit5d.detection.run", "run_detection", "C02"),
    "run_detection_frame_parallel": (
        "petakit5d.detection.run",
        "run_detection_frame_parallel",
        "C02",
    ),
    "track_settings": ("petakit5d.detection.tracking", "track_settings", "C03"),
    "run_tracking": ("petakit5d.detection.tracking", "run_tracking", "C03"),
    "run_track_processing": ("petakit5d.detection.processing", "run_track_processing", "C04b"),
    "rotate_tracks": ("petakit5d.detection.processing", "rotate_tracks", "C04b"),
}
STAGES: dict[str, Callable[..., Any] | None] = dict.fromkeys(_STAGE_HOMES)
# what a stage that is not there yet can be done with instead
_STAGE_HINTS = {
    "run_detection_frame_parallel": (
        "XR_runDetection3D_frame_parallel, the detection of parallel_method 'lock_distributed' "
        "(the default); parallel_method='parfor' runs run_detection (XR_runDetection3D) "
        "instead"
    ),
}


def _stage(name: str) -> Callable[..., Any]:
    fn = STAGES.get(name)
    if fn is not None:
        return fn
    module, attr, task = _STAGE_HOMES[name]
    try:
        import importlib

        fn = getattr(importlib.import_module(module), attr)
    except (ImportError, AttributeError):
        hint = f" ({_STAGE_HINTS[name]})" if name in _STAGE_HINTS else ""
        raise NotImplementedError(
            f"{name} ({module}.{attr}, task {task}) is not ported yet{hint}"
        ) from None
    return fn


def _call(name: str, *args: Any, **kwargs: Any) -> Any:
    """call a stage. The calls, with MATLAB's options in snake_case:

    - ``run_detection(data, sigma=, overwrite=, window_size=, results_path=, mode=,
      fit_mixtures=, max_mixtures=, detection_method=, background_correction=)`` (and
      ``run_detection_frame_parallel`` the same): ``XR_runDetection3D`` /
      ``XR_runDetection3D_frame_parallel``;
    - ``track_settings(radius=, max_gap_length=)``: ``loadTrackSettings``;
    - ``run_tracking(data, settings, file_name='trackedFeatures.mat', overwrite=,
      results_path=, parallel_method=)``: ``XR_runTracking3D``;
    - ``run_track_processing(data, overwrite=, tracker_output='trackedFeatures.mat',
      file_name='ProcessedTracks.mat', buffer=, buffer_all=, fit_mixtures=, window_size=,
      results_path=, parallel_method=)``: ``XR_runTrackProcessing3D``;
    - ``rotate_tracks(data, overwrite=, crop=, results_path=, matlab_quirks=)``:
      ``rotateTracks3D`` (``matlab_quirks``: its rotation of reversed scans, W20).

    ``data`` is the condition as MATLAB's struct array (``condition.struct_array``),
    ``results_path`` the list of each movie's results folder (ending with the separator),
    MATLAB's ``ResultsPath``."""
    return _stage(name)(*args, **kwargs)


# -- deskewData ---------------------------------------------------------------------------------


def _logical(v: Any, name: str) -> bool:
    if not isinstance(v, bool | np.bool_) and not (
        isinstance(v, np.ndarray) and v.dtype == np.bool_
    ):
        raise WrapperError(
            f"The value of '{name}' is invalid. It must satisfy the function: islogical.",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )
    return bool(np.asarray(v).reshape(-1)[0])


def _ds_done(movie: Any) -> bool:
    """MATLAB's skip test: channel 1's DS folder holds exactly as many ``*.tif*`` files as the
    movie has frames (measured: 2 or 4 for 3 frames do not skip, 3 ``.tiff`` files do)"""
    import fnmatch

    d = os.path.dirname(movie.frame_paths_ds[0][0])
    if not os.path.isdir(d):
        return False
    n = sum(1 for e in os.scandir(d) if e.is_file() and fnmatch.fnmatchcase(e.name, "*.tif*"))
    return n == len(movie.frame_paths[0])


def _save_settings(path: str, movie: Any) -> None:
    from petakit5d.io import write_mat

    write_mat(
        path,
        {
            "dz": movie.dz,
            "xyPixelSize": movie.pixel_size,
            "objectiveScan": movie.objective_scan,
            "reversed": movie.reversed,
        },
        version="7",
    )


def _deskew_frame(
    src: str, dst: str, movie: Any, skew_angle: float, reverse: bool
) -> tuple[float, float]:
    from petakit5d.deskew_rotate.volume import deskew_volume
    from petakit5d.io import read_tiff, write_tiff

    v = np.asarray(read_tiff(src), dtype=np.float64)
    if bool(np.asarray(movie.objective_scan).reshape(-1)[0]):
        ds = v
    else:
        ds = deskew_volume(
            v, float(skew_angle), float(movie.dz), float(movie.pixel_size), reverse, crop=False
        )
    write_tiff(ds.astype(np.float32), dst)
    return float(ds.min()), float(ds.max())


def deskew_data(
    data: Sequence[Any],
    *,
    overwrite: bool = False,
    skew_angle: float = 31.5,
    aname: str | None = None,
    rotate: bool = False,
    s_cmos_camera_flip: bool = False,
    crop: bool = False,
    llff_correction: bool = False,
    zoffset_correction: bool = False,
    reverse: bool | None = False,
    num_threads: int | None = None,
) -> Sequence[Any]:
    """``deskewData(data, 'Overwrite', overwrite, 'SkewAngle', skew_angle, 'aname', aname,
    ...)`` with the options the wrapper and the demos pass (module docstring): the DS frames of
    every movie, ``settings.mat`` and ``dynRange.mat``. Returns ``data`` unchanged.
    ``reverse``: deskewData's ``'Reverse'`` for all movies (False, as in MATLAB); None, each
    movie's scan direction (``reversed``: the wrapper's, W2). ``num_threads`` frames are
    deskewed at once (None: the CPUs available)."""
    ow = _logical(overwrite, "Overwrite")
    rev = None if reverse is None else _logical(reverse, "Reverse")
    for flag, name in ((rotate, "Rotate"), (crop, "Crop"), (llff_correction, "LLFFCorrection")):
        if _logical(flag, name):
            raise NotImplementedError(f"deskew_data: '{name}' is not ported")
    if _logical(zoffset_correction, "zoffsetCorrection"):
        raise NotImplementedError("deskew_data: 'zoffsetCorrection' is not ported")
    flip = _logical(s_cmos_camera_flip, "sCMOSCameraFlip")
    from concurrent.futures import ThreadPoolExecutor

    from petakit5d.core.threads import resolve_threads

    for movie in data:
        out = movie.source + ("" if aname is None or aname == "" else str(aname)) + os.sep
        os.makedirs(out, exist_ok=True)
        if not ow and _ds_done(movie):
            continue
        r = bool(np.asarray(movie.reversed).reshape(-1)[0]) if rev is None else rev
        _save_settings(os.path.join(out, "settings.mat"), movie)
        n_f = len(movie.frame_paths[0])
        n_ch = len(movie.channels)
        rng = np.zeros((n_f, 2, n_ch))
        for c in range(n_ch):
            os.makedirs(os.path.dirname(movie.frame_paths_ds[c][0]), exist_ok=True)
            if flip:  # MATLAB fails in its first frame (a parallel pool reports it otherwise)
                raise WrapperError(
                    "Undefined function 'GU_sCMOSCameraFlip' for input arguments of type 'double'.",
                    "MATLAB:UndefinedFunction",
                )
            jobs = list(zip(movie.frame_paths[c][:n_f], movie.frame_paths_ds[c][:n_f], strict=True))
            nt = max(1, min(resolve_threads(num_threads), len(jobs)))
            with ThreadPoolExecutor(nt) as ex:
                res = list(
                    ex.map(
                        lambda j, mv=movie, r=r: _deskew_frame(j[0], j[1], mv, skew_angle, r), jobs
                    )
                )
            for f, (lo, hi) in enumerate(res):
                rng[f, 0, c], rng[f, 1, c] = lo, hi
        from petakit5d.io import write_mat

        write_mat(os.path.join(out, "dynRange.mat"), {"dRange": rng}, version="7")
    return data


# -- XR_runDetTrack3d ---------------------------------------------------------------------------


def _movies_from_matlab(x: Any) -> Any:
    from petakit5d.detection.condition import Movie

    if x is None or isinstance(x, list | tuple) and x and isinstance(x[0], Movie):
        return x
    if isinstance(x, Movie):
        return [x]
    if isinstance(x, dict):
        return [Movie.from_matlab(x)]
    return [Movie.from_matlab(dict(e)) for e in x]


def _movies_to_matlab(x: Any) -> Any:
    from petakit5d.detection.condition import movies_to_matlab

    return None if x is None else movies_to_matlab(x)


def _is_pc() -> bool:
    """MATLAB's ``ispc``: XR_runDetTrack3d takes 'parfor' there whatever ``ParallelMethod``"""
    return os.name == "nt"


def _four(overwrite: Any) -> tuple[bool, bool, bool, bool]:
    a = np.asarray(overwrite)
    if a.dtype != np.bool_:
        raise WrapperError(
            "The value of 'Overwrite' is invalid. It must satisfy the function: islogical.",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )
    flags = a.reshape(-1)
    if flags.size == 1:
        flags = np.repeat(flags, 4)
    if flags.size < 4:
        raise WrapperError(
            f"Index exceeds the number of array elements. Index must not exceed {flags.size}.",
            "MATLAB:badsubscript",
        )
    return tuple(bool(v) for v in flags[:4])  # type: ignore[return-value]


@register(
    "XR_runDetTrack3d",
    Param(
        "data",
        "@isstruct",
        "The condition (load_condition_data's movies).",
        convert=(_movies_from_matlab, _movies_to_matlab),
    ),
    Param("apath", None, "Analysis root folder (None: each movie's first channel folder)."),
    Param("aname", None, "Results folder name under the first channel (None: 'Analysis')."),
    Param(
        "Overwrite",
        "@islogical",
        "Overwrite: one flag, or four (deskew, detection, tracking, processing).",
    ),
    Param("ParallelMethod", "@ischar", "'parfor' or 'lock_distributed' (detection's form)."),
    Param("SkewAngle", "@isscalar", "Skew angle in degrees."),
    Param("Sigma", None, "PSF sigmas, one row [sigma_xy, sigma_z] per channel.", convert="same"),
    Param("WindowSize", None, "Detection window size ([]: from sigma).", convert="same"),
    Param("Mode", "@ischar", "Fitted parameters of the detection."),
    Param("FitMixtures", "@islogical", "Fit mixtures of Gaussians."),
    Param("MaxMixtures", "@(x) numel(x)==1 && x>0 && round(x)==x", "Mixtures fitted at most."),
    Param("DetectionMethod", "@ischar", "'original', 'speedup', 'lowSNR' or 'adaptive'."),
    Param("BackgroundCorrection", "@islogical", "Correct the illumination along z."),
    Param("Track", "@islogical", "Track after detection."),
    Param("RotateTracks", "@islogical", "Rotate the tracks."),
    Param("TrackingRadius", "@(x) numel(x)==2", "Linking radius [min max]."),
    Param("TrackingGapLength", "@(x) numel(x)==1", "Frames a gap may span."),
    Param("Buffer", None, "Frames read out before and after each track.", convert="same"),
    Param("BufferAll", "@islogical", "Buffer all tracks."),
    Param("sCMOSCameraFlip", "@islogical", "sCMOS camera flip (MATLAB's deskew fails on it)."),
    Param("Crop", "@islogical", "Crop for the deskew (interactive in MATLAB; not ported)."),
    Param("Rotate", "@islogical", "Rotate the deskewed frames (not ported)."),
    Param("LLFFCorrection", "@islogical", "Light-sheet flat-field correction (not ported)."),
    Param("LowerLimit", "@isnumeric", "Flat-field mask's lower bound."),
    Param("LSImageCh1", "@isstr", "Light-sheet image of channel 1.", convert="same"),
    Param("LSImageCh2", "@isstr", "Light-sheet image of channel 2.", convert="same"),
    Param("LSImageCh3", "@isstr", "Light-sheet image of channel 3.", convert="same"),
    Param("LSImageCh4", "@isstr", "Light-sheet image of channel 4.", convert="same"),
    Param("BackgroundCh1", "@isstr", "Background image of channel 1."),
    Param("BackgroundCh2", "@isstr", "Background image of channel 2."),
    Param("BackgroundCh3", "@isstr", "Background image of channel 3."),
    Param("BackgroundCh4", "@isstr", "Background image of channel 4."),
    Param("XZoffsetCorrection", "@islogical", "Correct the channels' x and z offsets."),
    Param(
        "matlab_quirks",
        None,
        "MATLAB's two slips: reversed scans deskewed as forward ones (W2), the tracking stages "
        "given the folders of all movies (W3); and its rotation of reversed scans' tracks by the "
        "forward angle (W20).",
        python_only=True,
    ),
    keep_unmatched=True,
    task="C04",
)
def run_detection_tracking(
    data: Any = None,
    apath: str | None = None,
    aname: str | None = None,
    *,
    overwrite: Any = False,
    parallel_method: str = "lock_distributed",
    skew_angle: float = 31.5,
    sigma: Any = ((1.26, 1.32), (1.41, 1.38), (1.58, 1.608)),
    window_size: Any = None,
    mode: str = "xyzAc",
    fit_mixtures: bool = False,
    max_mixtures: int = 2,
    detection_method: str = "lowSNR",
    background_correction: bool = False,
    track: bool = True,
    rotate_tracks: bool = True,
    tracking_radius: Any = (3, 6),
    tracking_gap_length: Any = 2,
    buffer: Any = (3, 3),
    buffer_all: bool = False,
    s_cmos_camera_flip: bool = False,
    crop: bool = False,
    rotate: bool = False,
    llff_correction: bool = False,
    lower_limit: Any = 0.4,
    ls_image_ch_1: str = "",
    ls_image_ch_2: str = "",
    ls_image_ch_3: str = "",
    ls_image_ch_4: str = "",
    background_ch_1: str = "",
    background_ch_2: str = "",
    background_ch_3: str = "",
    background_ch_4: str = "",
    x_zoffset_correction: bool = False,
    matlab_quirks: bool = False,
) -> None:
    """Deskew, channel offsets, detection, tracking, track processing and track rotation of a
    condition (MATLAB's ``XR_runDetTrack3d``; module docstring). ``matlab_quirks``: MATLAB's
    output for its two slips (W2, W3) and for the rotation of reversed scans (W20) instead of
    the fixed ones."""
    from petakit5d.detection.condition import Condition, short_path
    from petakit5d.detection.offsets import correct_xz_offset_data

    method = "parfor" if _is_pc() else parallel_method
    ow = _four(overwrite)
    if sigma is None or np.size(sigma) == 0:
        raise WrapperError(
            "MATLAB asks for the Gaussian s.d. of the detection interactively; give sigma.",
            "MATLAB:services:MissingRequiredCapability",
        )
    if data is None or len(data) == 0:
        raise WrapperError(
            "MATLAB loads the condition interactively; give data (load_condition_data).",
            "MATLAB:services:MissingRequiredCapability",
        )
    data = Condition(data)
    if apath is None or np.size(apath) == 0:
        paths = [m.source for m in data]
    else:
        paths = [f"{apath}{short_path(m.channels[0], 3)}{os.sep}" for m in data]
    name = "Analysis" if aname is None or np.size(aname) == 0 else str(aname)
    paths = [f"{p}{name}{os.sep}" for p in paths]
    for p in paths:
        os.makedirs(p, exist_ok=True)
    data = Condition(
        deskew_data(
            data,
            overwrite=ow[0],
            skew_angle=skew_angle,
            aname=aname,
            rotate=rotate,
            s_cmos_camera_flip=s_cmos_camera_flip,
            crop=crop,
            llff_correction=llff_correction,
            zoffset_correction=False,
            reverse=False if matlab_quirks else None,
        )
    )
    if x_zoffset_correction:
        correct_xz_offset_data(data, overwrite=ow[0])
    det = {
        "sigma": sigma,
        "overwrite": ow[1],
        "window_size": window_size,
        "mode": mode,
        "fit_mixtures": fit_mixtures,
        "max_mixtures": max_mixtures,
        "detection_method": detection_method,
        "background_correction": background_correction,
    }
    from petakit5d.detection.condition import struct_array

    fields = list(data[0].to_matlab())
    stage = {"parfor": "run_detection", "lock_distributed": "run_detection_frame_parallel"}
    if method in stage:
        _call(stage[method], struct_array(data, fields), results_path=list(paths), **det)
    if not track:
        return
    # the movies of more than one frame with their own results folders. MATLAB passes the
    # folders of all movies, so that its i-th movie gets the i-th folder (W3)
    keep = [m for m in data if m.movie_length != 1]
    own = [p for m, p in zip(data, paths, strict=True) if m.movie_length != 1]
    structs = struct_array(keep, fields)
    folders = list(paths) if matlab_quirks else own
    settings = _call("track_settings", radius=tracking_radius, max_gap_length=tracking_gap_length)
    _call(
        "run_tracking",
        structs,
        settings,
        file_name="trackedFeatures.mat",
        overwrite=ow[2],
        results_path=folders,
        parallel_method=method,
    )
    _call(
        "run_track_processing",
        structs,
        overwrite=ow[3],
        tracker_output="trackedFeatures.mat",
        file_name="ProcessedTracks.mat",
        buffer=buffer,
        buffer_all=buffer_all,
        fit_mixtures=fit_mixtures,
        window_size=window_size,
        results_path=folders,
        parallel_method=method,
    )
    if rotate_tracks:
        _call(
            "rotate_tracks",
            structs,
            overwrite=ow[3],
            crop=crop,
            results_path=folders,
            matlab_quirks=matlab_quirks,
        )
