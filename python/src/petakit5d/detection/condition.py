"""Loading a condition: the movies of a condition folder, their channels and frames, as MATLAB's
``XR_loadConditionData3D`` and ``XR_loadConditionData`` give them (lab-edited files with parts of
another author: written from their help text and from MATLAB's behaviour, measured by black-box
probes; ``load_condition_data`` is registered with the interface of F03's stored inputParser
extract of ``XR_loadConditionData3D``: its parameters' names, order, validators and defaults).

A condition is a list of ``Movie`` (``Movie2D`` for ``load_condition_data_2d``), one per movie,
MATLAB's ``data`` struct array; ``to_matlab()`` gives the struct with MATLAB's field names and
values (``imagesize`` ``[ny nx nz]``, columns of paths) and ``Movie.from_matlab`` reads one back.
``struct_array(movies)`` is the condition as the port holds MATLAB's struct array (``io.matfile``'s
form: one movie's dict, more movies a 1 x n ``StructArray``); MATLAB's layer returns that.

What the loaders do (measured, oracle tree 45b2bb1):

- **movies**: the condition folder itself and its subfolders down to ``maxlevel`` (3) levels whose
  names contain ``movie_selector`` (any case), level by level; within a level ordered by the
  number after the selector (``ignore_exp_number``: by name, ASCII). Without
  ``ignore_exp_number`` a condition folder without a number is no movie, and an unnumbered
  movie folder below it is an error here (MATLAB errs or returns another movie twice). No
  movie: ``No movies found in: <folder>``. A relative folder is taken from the
  current folder (``..`` kept as it is); folders end with the separator.
- **channels**: ``<movie>/<name>/`` for each of ``ch_names``; ``source`` is the first. A missing
  channel folder is MATLAB's ``MATLAB:Java:InvalidInput``. The deskewed and rotated frames are
  ``<channel>DS/<name>.tif`` and ``<channel>DSR/<name>.tif`` whatever a frame's extension.
- **frames**: the ``.tif`` and ``.tiff`` files (any case) of each channel folder that do not start
  with ``.``, ordered by name (ASCII; not by the stack or time tokens); every channel keeps as
  many as channel 1 has (MATLAB's logical index: more are cut, fewer are its
  ``MATLAB:badsubscript``); an empty folder is ``MATLAB:cellRefFromNonCell``.
- **date**: every run of six digits of the movie's folder path, in order, joined
  (``c_111111/Ex2_222222`` gives ``'111111222222'``, ``20150407`` gives ``'201504'``); None
  (MATLAB's ``[]``) without one. The text depends on where the data sits; compare it masked.
- **frame rate** (3D): the median of the differences of the first ``<digits>msec`` token of
  channel 1's frame paths, in the frames' order, over 1000 (negative when the names are not in
  time order; NaN for one frame); a frame without one is ``MATLAB:badsubscript``. ``frame_rate``
  replaces it. The 2D loader takes the first ``<digits>ms`` of the movie folder's path over 1000,
  2 without one.
- **dz**: ``dz`` when given, else ``0.<digits>`` of the first ``_z0p<digits>`` (any case) in
  channel 1's folder path (the condition's folders, the movie's, the channel's); NaN without one.
  ``zAniso = dz * sin(angle*pi/180) / pixelSize`` (fdlibm's sine: not ``sind``);
  ``dz / pixelSize`` for objective scans.
- **image size**: ``(nz, ny, nx)`` of channel 1's first frame; one-page frames are MATLAB's
  ``MATLAB:badsubscript`` in the 3D loader, and so is a channel 1 without a ``.tif`` frame (all
  ``.tiff`` or ``.TIF``).
- **masks**: ``<source>Analysis/Masks/dmask_%03d.tif`` (3D), ``<source>Detection/Masks/
  dmask_%0<n>d.tif`` with ``n`` the digits of the frame count (2D).
- **markers** are kept as given and must name one of MATLAB's fluorophores (``FLUOROPHORES``, any
  case); an unknown one, missing markers or channels make MATLAB ask interactively, which the
  port does not do (``MATLAB:services:MissingRequiredCapability``, MATLAB's error in batch mode).
- options are kept as given (``angle``, ``pixel_size``, ``objective_scan``, ``reversed``);
  ``channel_patterns`` changes nothing (measured on flat and channel-folder layouts).
"""

from __future__ import annotations

import math
import os
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, fields
from typing import Any, ClassVar

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "FLUOROPHORES",
    "Condition",
    "ConditionError",
    "Movie",
    "Movie2D",
    "find_movies",
    "load_condition_data",
    "load_condition_data_2d",
    "movies_to_matlab",
    "short_path",
    "struct_array",
]

# the names getFluorPropStruct() returns in MATLAB (a black-box call), which the loaders accept
FLUOROPHORES = (
    "bfp", "ebfp", "dapi", "cfp", "egfp", "gfp", "mneongreen", "alexa488", "yfp", "meos",
    "alexa555", "cy3", "turborfp", "tmr", "dtomato", "tdtomato", "dsred", "tagrfp", "cyofp",
    "mruby2", "alexa568", "fluosphere605", "rfp", "mrfp", "mcherry", "texasred", "alexa594",
    "alexa647", "cy5",
)  # fmt: skip

_TIFF = re.compile(r"\.tiff?$", re.IGNORECASE)
_SIX = re.compile(r"\d{6}")
_MSEC = re.compile(r"(\d+)msec")
_MS = re.compile(r"(\d+)ms")
_DZ = re.compile(r"_z0p(\d+)", re.IGNORECASE)
_INTERACTIVE = "MATLAB:services:MissingRequiredCapability"


class ConditionError(ValueError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _sep(p: str) -> str:
    return p if p.endswith(("/", os.sep)) else p + os.sep


def _cell_column(paths: Sequence[str]) -> np.ndarray:
    a = np.empty((len(paths), 1), dtype=object)
    for i, p in enumerate(paths):
        a[i, 0] = p
    return a


def _empty_double() -> np.ndarray:
    return np.zeros((0, 0))


@dataclass
class Movie2D:
    """One movie as ``XR_loadConditionData`` gives it (an element of its ``data``)."""

    source: str
    channels: list[str]
    date: str | None
    framerate: float
    image_size: tuple[int, ...]  # (nz, ny, nx) of channel 1's first frame
    movie_length: int
    markers: list[str]
    frame_paths: list[list[str]]
    mask_paths: list[str]
    na: Any
    m: Any
    pixel_size: Any

    _MATLAB: ClassVar[dict[str, str]] = {
        "image_size": "imagesize",
        "movie_length": "movieLength",
        "frame_paths": "framePaths",
        "mask_paths": "maskPaths",
        "na": "NA",
        "m": "M",
        "pixel_size": "pixelSize",
    }

    def to_matlab(self) -> dict[str, Any]:
        """MATLAB's struct: its field names, ``imagesize`` ``[ny nx nz]``, path columns"""
        out: dict[str, Any] = {}
        for f in fields(self):
            v = getattr(self, f.name)
            key = self._MATLAB.get(f.name, f.name)
            if f.name == "image_size":
                s = tuple(int(n) for n in v)
                v = np.array([s[1], s[2], s[0]] if len(s) == 3 else [s[0], s[1]], dtype=np.float64)
            elif f.name in ("movie_length",):
                v = float(v)
            elif f.name == "date":
                v = _empty_double() if v is None else v
            elif f.name in ("frame_paths", "frame_paths_ds", "frame_paths_dsr"):
                v = [_cell_column(c) for c in v]
            elif f.name == "mask_paths":
                v = _cell_column(v)
            out[key] = v
        return out

    @classmethod
    def from_matlab(cls, d: Mapping[str, Any]):
        """a movie from MATLAB's struct (``to_matlab``'s form)"""
        inv = {v: k for k, v in cls._MATLAB.items()}
        kw: dict[str, Any] = {}
        for k, v in d.items():
            name = inv.get(k, k)
            if name == "image_size":
                s = [int(x) for x in np.asarray(v).reshape(-1)]
                v = (s[2], s[0], s[1]) if len(s) == 3 else tuple(s)
            elif name == "movie_length":
                v = int(np.asarray(v).reshape(-1)[0])
            elif name == "date":
                v = None if not isinstance(v, str) else v
            elif name in ("frame_paths", "frame_paths_ds", "frame_paths_dsr"):
                v = [[str(p) for p in np.asarray(c, dtype=object).reshape(-1)] for c in v]
            elif name in ("mask_paths", "channels", "markers"):
                v = [str(p) for p in np.asarray(v, dtype=object).reshape(-1)]
            kw[name] = v
        return cls(**kw)


@dataclass
class Movie(Movie2D):
    """One movie of a condition as ``XR_loadConditionData3D`` gives it (an element of MATLAB's
    ``data``): ``Movie2D``'s fields without ``NA`` and ``M``, and the deskew's geometry and the
    deskewed and rotated frames' paths."""

    na: Any = field(default=None, repr=False)
    m: Any = field(default=None, repr=False)
    pixel_size: Any = 0.104
    dz: Any = None
    z_aniso: Any = None
    angle: Any = 31.5
    objective_scan: Any = False
    reversed: Any = False
    frame_paths_ds: list[list[str]] = field(default_factory=list)
    frame_paths_dsr: list[list[str]] = field(default_factory=list)

    _MATLAB: ClassVar[dict[str, str]] = {
        **Movie2D._MATLAB,
        "z_aniso": "zAniso",
        "objective_scan": "objectiveScan",
        "frame_paths_ds": "framePathsDS",
        "frame_paths_dsr": "framePathsDSR",
    }

    def to_matlab(self) -> dict[str, Any]:
        d = super().to_matlab()
        del d["NA"], d["M"]
        return d


def movies_to_matlab(movies: Sequence[Movie2D]) -> list[dict[str, Any]]:
    """MATLAB's struct array of a condition (a list of structs)"""
    return [m.to_matlab() for m in movies]


def struct_array(movies: Sequence[Movie2D], fields: Sequence[str] = ()) -> Any:
    """the movies as the port holds MATLAB's struct array (``io.matfile``'s form): one movie's
    dict (``Movie.to_matlab()``), else a 1 x n ``StructArray`` of them (empty: of ``fields``)"""
    from petakit5d.io.matfile import StructArray

    rows = [m.to_matlab() for m in movies]
    if len(rows) == 1:
        return rows[0]
    return StructArray((1, len(rows)), list(rows[0]) if rows else list(fields), rows)


class Condition(list):
    """A condition: its movies (``Movie`` or ``Movie2D``), MATLAB's ``data`` struct array. It
    counts as a MATLAB struct for the parameter registry's ``isstruct`` check (registered as a
    ``Mapping``; it is a list otherwise)."""

    def to_matlab(self) -> list[dict[str, Any]]:
        return movies_to_matlab(self)


Mapping.register(Condition)


def short_path(source: str, n: int = 3) -> str:
    """MATLAB's ``getShortPath(data, n)`` (another author's GPL function; from its outputs): the
    last ``n`` folder names of ``source``'s folder, joined by ``/``, without the text before the
    first separator (a relative path's first folder)."""
    parts = re.split(r"[\\/]", source)
    parts = parts[1:-1]
    return "/".join(parts[-n:]) if n > 0 else ""


# -- the 2D loader ------------------------------------------------------------------------------


def _absolute(folder: str) -> str:
    if not os.path.isabs(folder):
        folder = os.path.join(os.getcwd(), folder)
    return _sep(folder)


def _movie_number(name: str, selector: str) -> float | None:
    m = re.search(re.escape(selector) + r"(\d+)", name, re.IGNORECASE)
    return None if m is None else float(m.group(1))


def _name(p: str) -> str:
    return os.path.basename(p.rstrip("/\\"))


def find_movies(
    cond_dir: str, movie_selector: str = "Ex", maxlevel: int = 3, ignore_exp_number: bool = False
) -> list[str]:
    """The movie folders of a condition folder (module docstring), each ending with the
    separator, in MATLAB's order."""
    root = _absolute(cond_dir)
    sel = movie_selector.lower()
    levels: list[list[str]] = [[root]]
    for _ in range(int(maxlevel)):
        nxt = []
        for d in levels[-1]:
            try:
                entries = sorted(os.scandir(d), key=lambda e: e.name)
            except OSError:
                continue
            nxt += [_sep(e.path) for e in entries if e.is_dir()]
        levels.append(nxt)
    found = [[p for p in level if sel in _name(p).lower()] for level in levels]
    if not any(found):
        raise ConditionError(f"No movies found in: {root}")
    if ignore_exp_number:
        return [p for level in found for p in level]
    out: list[str] = []
    unnumbered = False
    for depth, level in enumerate(found):
        nums = [_movie_number(_name(p), movie_selector) for p in level]
        # a condition folder without a number is no movie; an unnumbered movie folder below it
        # breaks MATLAB's ordering (an error, or another movie twice: not reproduced)
        unnumbered |= depth > 0 and any(n is None for n in nums)
        ranked = sorted(((n, p) for n, p in zip(nums, level, strict=True) if n is not None),
                        key=lambda t: t[0])  # fmt: skip
        out += [p for _, p in ranked]
    if not out:
        raise ConditionError(
            "Dot indexing is not supported for variables of this type.",
            "MATLAB:structRefFromNonStruct",
        )
    if unnumbered:
        raise ConditionError(
            "Unable to perform assignment because the left and right sides have a different "
            "number of elements.",
            "MATLAB:matrix:singleSubscriptNumelMismatch",
        )
    return out


def _frames(channel: str) -> list[str]:
    try:
        names = sorted(
            e.name
            for e in os.scandir(channel)
            if e.is_file() and not e.name.startswith(".") and _TIFF.search(e.name)
        )
    except FileNotFoundError:
        raise ConditionError(
            "Scalar input must be a java object", "MATLAB:Java:InvalidInput"
        ) from None
    return [channel + n for n in names]


def _check_inputs(ch_names: Any, markers: Any) -> tuple[list[str], list[str]]:
    if not isinstance(ch_names, list | tuple) or len(ch_names) == 0:
        raise ConditionError(
            "MATLAB asks for the channel folders interactively; give ch_names.", _INTERACTIVE
        )
    if not isinstance(markers, list | tuple) or len(markers) == 0:
        raise ConditionError(
            "MATLAB asks for the fluorescent markers interactively; give markers.", _INTERACTIVE
        )
    if len(markers) < len(ch_names):
        raise ConditionError(
            f"Index exceeds the number of array elements. Index must not exceed {len(markers)}.",
            "MATLAB:badsubscript",
        )
    bad = [mk for mk in markers if str(mk).lower() not in FLUOROPHORES]
    if bad:
        raise ConditionError(
            f"Unknown fluorescent marker(s) {bad}: MATLAB asks for them interactively.",
            _INTERACTIVE,
        )
    return [str(c) for c in ch_names], [str(m) for m in markers]


def _image_size(path: str) -> tuple[int, ...]:
    """MATLAB's imfinfo size: (pages, height, width); one page: (height, width)"""
    import cpptiff

    s = tuple(int(n) for n in cpptiff.get_image_shape(path))
    return s[1:3] if s[0] == 1 else s[:3]


def _load_movies(
    cond_dir: Any,
    ch_names: Any,
    markers: Any,
    *,
    movie_selector: str,
    ignore_exp_number: bool,
    maxlevel: Any,
) -> list[tuple[str, list[str], list[list[str]], str | None]]:
    if not isinstance(cond_dir, str) or not cond_dir:
        raise ConditionError(
            "MATLAB asks for the condition folder interactively; give cond_dir.", _INTERACTIVE
        )
    chs, _ = _check_inputs(ch_names, markers)
    level = 3 if maxlevel is None or np.size(maxlevel) == 0 else maxlevel
    out = []
    for movie in find_movies(cond_dir, movie_selector, level, ignore_exp_number):
        channels = [_sep(os.path.join(movie, c)) for c in chs]
        lists = []
        for c in channels:
            if not os.path.isdir(c):
                raise ConditionError(
                    "Scalar input must be a java object", "MATLAB:Java:InvalidInput"
                )
            lists.append(_frames(c))
        n = len(lists[0])
        if n == 0 or any(len(fr) == 0 for fr in lists):
            raise ConditionError(
                "Brace indexing is not supported for variables of this type.",
                "MATLAB:cellRefFromNonCell",
            )
        for k, fr in enumerate(lists[1:], start=1):
            if len(fr) < n:
                raise ConditionError(
                    "The logical indices contain a true value outside of the array bounds.",
                    "MATLAB:badsubscript",
                )
            lists[k] = fr[:n]
        digits = _SIX.findall(movie)
        out.append((movie, channels, lists, "".join(digits) if digits else None))
    return out


def load_condition_data_2d(
    cond_dir: Any = None,
    ch_names: Any = (),
    markers: Any = (),
    maxlevel: Any = 3,
    *,
    parameters: Sequence[float] = (1.4, 150, 1.6e-05),
    movie_selector: str = "Ex",
    ignore_exp_number: bool = False,
    frame_rate: Any = None,
) -> Condition:
    """MATLAB's ``XR_loadConditionData(condDir, chNames, markers, ...)`` (module docstring):
    the movies of a condition with ``NA``, ``M`` and the camera pixel size of ``parameters``
    (``[NA M pixelSize]``) and the frame rate from the movie folder's ``<digits>ms`` token."""
    if np.size(parameters) != 3:
        raise ConditionError(
            "The value of 'Parameters' is invalid. It must satisfy the function: @(x)numel(x)==3.",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )
    na, mag, ps = (float(v) for v in np.asarray(parameters, dtype=np.float64).reshape(-1))
    out = []
    for movie, channels, lists, date in _load_movies(
        cond_dir,
        ch_names,
        markers,
        movie_selector=movie_selector,
        ignore_exp_number=ignore_exp_number,
        maxlevel=maxlevel,
    ):
        if frame_rate is not None:
            rate = frame_rate
        else:
            mt = _MS.search(movie)
            rate = float(mt.group(1)) / 1000 if mt else 2.0
        n = len(lists[0])
        digits = len(str(n))
        out.append(
            Movie2D(
                source=channels[0],
                channels=channels,
                date=date,
                framerate=rate,
                image_size=_image_size(lists[0][0]),
                movie_length=n,
                markers=[str(m) for m in markers],
                frame_paths=lists,
                mask_paths=[
                    f"{channels[0]}Detection{os.sep}Masks{os.sep}dmask_{k:0{digits}d}.tif"
                    for k in range(1, n + 1)
                ],
                na=na,
                m=mag,
                pixel_size=ps,
            )
        )
    return Condition(out)


# -- the 3D loader ------------------------------------------------------------------------------


def _stem(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def _frame_rate(paths: Sequence[str]) -> float:
    from petakit5d.matlabcompat.stats import median

    t = []
    for p in paths:
        mt = _MSEC.search(p)
        if mt is None:
            raise ConditionError("Index exceeds array bounds.", "MATLAB:badsubscript")
        t.append(float(mt.group(1)))
    d = np.diff(np.array(t, dtype=np.float64))
    return float(np.asarray(median(d)).reshape(-1)[0]) / 1000


def _dz_token(source: str) -> float:
    from petakit5d.core import dtypes

    mt = _DZ.search(source)
    return float("0." + mt.group(1)) if mt else dtypes.NAN


def _z_aniso(dz: Any, angle: Any, pixel_size: Any, objective_scan: Any) -> Any:
    from petakit5d.matlabcompat.elementary import _sin

    dzf = float(np.asarray(dz, dtype=np.float64).reshape(-1)[0])
    ps = _as_number(pixel_size)
    if bool(np.asarray(objective_scan).reshape(-1)[0]):
        return dzf / ps
    a = _as_number(angle)
    s = float(_sin(np.array([a * math.pi / 180]))[0])
    return dzf * s / ps


def _as_number(v: Any) -> float:
    """a MATLAB scalar as a double (a char counts as its code, as in MATLAB's arithmetic)"""
    if isinstance(v, str):
        return float(ord(v)) if len(v) == 1 else float("nan")
    return float(np.asarray(v, dtype=np.float64).reshape(-1)[0])


@register(
    "XR_loadConditionData3D",
    Param(
        "condDir",
        "@(x) ischar(x) && ~any(strcmpi(x,     {'Parameters', 'MovieSelector', "
        "'IgnoreEmptyFolders', 'FrameRate'}))",
        "The condition folder (its movie folders), or one movie's folder.",
    ),
    Param("chNames", None, "The channel folders' names, the master channel first."),
    Param("markers", None, "The channels' fluorescent markers (getFluorPropStruct's names)."),
    Param(
        "maxlevel",
        "@(x) isempty(x) || (isnumeric(x) && abs(round(x))==x)",
        "Folder levels below condDir searched for movies.",
    ),
    Param("MovieSelector", "@ischar", "Text in a movie folder's name."),
    Param("IgnoreExpNumber", "@islogical", "Order the movies by name, not by their number."),
    Param("ChannelPatterns", "@iscell", "Channel patterns (no effect: measured)."),
    Param("Angle", "@isscalar", "Skew angle in degrees."),
    Param("PixelSize", "@isscalar", "Pixel size."),
    Param("dz", "@isnumeric", "Scan step ([]: from the folders' _z0p<digits> name)."),
    Param("ObjectiveScan", "@isscalar", "Objective scan (zAniso = dz / PixelSize)."),
    Param("Reversed", "@isscalar", "Reversed scan."),
    Param("FrameRate", "@isscalar", "Frame rate ([]: from the frames' <digits>msec names)."),
    keep_unmatched=True,
    task="C04",
    matlab_out=struct_array,
)
def load_condition_data(
    cond_dir: Any = None,
    ch_names: Any = (),
    markers: Any = (),
    maxlevel: Any = 3,
    *,
    movie_selector: str = "Ex",
    ignore_exp_number: bool = False,
    channel_patterns: Sequence[str] = (),
    angle: Any = 31.5,
    pixel_size: Any = 0.104,
    dz: Any = None,
    objective_scan: Any = False,
    reversed: Any = False,  # noqa: A002
    frame_rate: Any = None,
) -> Condition:
    """Load a condition's movies (MATLAB's ``XR_loadConditionData3D(condDir, chNames,
    markers, ...)``; module docstring): one ``Movie`` per movie folder, with the frames of each
    channel, the frame rate, the scan step and the z anisotropy of the deskewed frames, and the
    paths of the deskewed (``DS``) and rotated (``DSR``) frames. MATLAB's interactive branches
    (no folder, channels or markers) are not ported: they are errors."""
    base = _load_movies(
        cond_dir,
        ch_names,
        markers,
        movie_selector=movie_selector,
        ignore_exp_number=ignore_exp_number,
        maxlevel=maxlevel,
    )
    out = []
    for _movie, channels, lists, date in base:
        n = len(lists[0])
        if not any(f.endswith(".tif") for f in lists[0]):
            # (measured: frames that are all .tiff or .TIF fail in MATLAB's 3D loader)
            raise ConditionError("Index exceeds array bounds.", "MATLAB:badsubscript")
        size = _image_size(lists[0][0])
        if len(size) < 3:
            raise ConditionError(
                "Index exceeds the number of array elements. Index must not exceed 2.",
                "MATLAB:badsubscript",
            )
        rate = frame_rate if frame_rate is not None else _frame_rate(lists[0])
        d = _dz_token(channels[0]) if dz is None or np.size(dz) == 0 else dz
        ds = [
            [f"{c}DS{os.sep}{_stem(p)}.tif" for p in fr]
            for c, fr in zip(channels, lists, strict=True)
        ]
        dsr = [
            [f"{c}DSR{os.sep}{_stem(p)}.tif" for p in fr]
            for c, fr in zip(channels, lists, strict=True)
        ]
        out.append(
            Movie(
                source=channels[0],
                channels=channels,
                date=date,
                framerate=rate,
                image_size=size,
                movie_length=n,
                markers=[str(m) for m in markers],
                frame_paths=lists,
                mask_paths=[
                    f"{channels[0]}Analysis{os.sep}Masks{os.sep}dmask_{k:03d}.tif"
                    for k in range(1, n + 1)
                ],
                pixel_size=pixel_size,
                dz=d,
                z_aniso=_z_aniso(d, angle, pixel_size, objective_scan),
                angle=angle,
                objective_scan=objective_scan,
                reversed=reversed,
                frame_paths_ds=ds,
                frame_paths_dsr=dsr,
            )
        )
    return Condition(out)
