"""The port's own checks of the parameters of third-party MATLAB functions, against MATLAB.

The repository keeps no validator text of third-party MATLAB code (decided on 2026-11-20). A
parameter that MATLAB validates gets a check of the port's own instead (``Param(check=...)``),
made from MATLAB's behaviour. Each case here calls the MATLAB function of one such parameter
with a battery of designed values (``matlab/pkpy_x27_probe.m``) and records, for each value,
whether the call gets past MATLAB's inputParser or the identifier of inputParser's error, never
a message (inputParser's messages print the validator). The Python side binds the same values
through the spec's MATLAB-name layer: its names, positions and checks, without the conversions
to Python's form. Values outside the batteries may still be taken otherwise.

``data`` of the tracking and track-processing wrappers: MATLAB computes the defaults of
``ResultsPath`` and ``Frames`` from ``data`` before it checks it, so a value that is not a
struct fails there, outside inputParser, unless it is empty. Its battery holds structs and the
empty values of each class, which reach the checks.

``runLifetimeAnalysis3D``'s options are probed with one movie whose tracks are missing: its
``data`` passes MATLAB's check (which takes the movies' frame rates, so the port's function
checks ``data`` itself and its spec has no check of the port's), and each call that the parser
takes fails reading the tracks. The helper counts any error other than inputParser's as taken,
an error that a check raises itself too: the port's checks take what MATLAB's parser does not
refuse with its identifiers (``Display``'s cells of other than one element, ``Binning``'s cells,
structs, strings and handles).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from petakit5d.params import MatlabArgumentError
from petakit5d.params import values as mv
from petakit5d.params.registry import REGISTRY
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, mcell, mcol, mempty, mstring, mtyped
from petakit5d.testing.compare import (
    MArray,
    MCell,
    MFunc,
    MPlaceholder,
    MString,
    MStruct,
    MValue,
    to_mvalue,
)

HELPERS = ["tests/params/matlab"]


def _cell(items: list, size: tuple[int, ...] | None = None) -> MCell:
    vals = [x if isinstance(x, MValue) else to_mvalue(x) for x in items]
    return MCell(size or ((1, len(vals)) if vals else (0, 0)), vals)


def _structs(size: tuple[int, int], fields: list[str], rows: list[dict]) -> MStruct:
    return MStruct(size, fields, [{k: to_mvalue(v) for k, v in r.items()} for r in rows])


EMPTY_CELL = _cell([])

# the battery every probed parameter gets: classes, sizes and values
GENERIC: dict[str, object] = {
    "true": True,
    "false": False,
    "logical_row": [True, False],
    "logical_empty": mempty("logical"),
    "zero": 0.0,
    "one": 1.0,
    "two": 2.0,
    "half": 2.5,
    "negative": -3.0,
    "nan": math.nan,
    "inf": math.inf,
    "row2": [1.0, 2.0],
    "row3": [1.0, 2.0, 3.0],
    "row3_uneven": [1.0, 2.0, 4.0],
    "row3_const": [5.0, 5.0, 5.0],
    "col2": mcol([1.0, 2.0]),
    "matrix": np.array([[1.0, 2.0], [3.0, 4.0]]),
    "empty": None,
    "empty_1x0": mempty("double", (1, 0)),
    "uint8": mtyped(3, "uint8"),
    "int16_row": mtyped([1, 2], "int16"),
    "single": mtyped(1.5, "single"),
    "char": "abc",
    "char_empty": "",
    "char_one": "a",
    "string": mstring("abc"),
    "cell_empty": EMPTY_CELL,
    "cell_char": mcell("a"),
    "cell_chars": mcell("a", "b"),
    "cell_nums": mcell(1.0, 2.0),
    "struct": {"a": 1.0},
    "struct_row": _structs((1, 2), ["a"], [{"a": 1.0}, {"a": 2.0}]),
    "struct_empty": MStruct((0, 0), [], []),
    "handle": MFunc("@sin"),
}

# designed values by what a parameter is
INTERP = {"linear": "linear", "cubic": "cubic", "nearest": "nearest", "upper": "LINEAR",
          "mixed_case": "Cubic", "string_linear": mstring("linear"), "cell_linear": mcell("linear"),
          "cell_cubic_case": mcell("Cubic"), "cell_linear_other": mcell("linear", "abc"),
          "cell_others": mcell("abc", "def"), "cell_row_empty": MCell((1, 0), [])}  # fmt: skip
FRAMES = {"stride2": [1.0, 3.0, 5.0], "descending": [3.0, 2.0, 1.0], "five": 5.0,
          "uint8_row": mtyped([1, 2, 3], "uint8"), "cell_of_row": mcell([1.0, 2.0, 3.0]),
          "nan_row": [1.0, math.nan, 3.0], "inf_row": [1.0, math.inf],
          "inf_pair": [math.inf, math.inf],
          "matrix_uneven": np.array([[1.0, 2.0], [3.0, 5.0]])}  # fmt: skip
PARALLEL = {"default": "lock_distributed", "upper": "LOCK_DISTRIBUTED"}
LEVELS = {"default": 3.0, "int8_negative": mtyped(-2, "int8"), "row_half": [1.0, 2.5],
          "minus_inf": -math.inf, "big": 1e6, "row_negative": [-1.0, -2.0],
          "row_nan": [1.0, math.nan], "col_half": mcol([1.0, 2.5]),
          "matrix_negative": np.array([[-1.0, 2.0], [3.0, 4.0]]),
          "int8_row_negative": mtyped([-1, -2], "int8"),
          "single_row": mtyped([1.5, 2.5], "single")}  # fmt: skip
PAIRS = {"default": [3.0, 3.0], "char2": "ab"}
# text that names a parameter in another case or in part (the loader's first argument)
COND_TEXT = {"folder": "cond_missing", "lower_MovieSelector": "movieselector",
             "upper_FrameRate": "FRAMERATE", "prefix_MovieSelector": "MovieSel"}  # fmt: skip
# the lifetime analysis' text options: the default, other words, their case, string and cell
ON_OFF = {"on": "on", "off": "off", "upper_off": "OFF", "string_off": mstring("off"),
          "cell_off": mcell("off")}  # fmt: skip
SCREEN = {"screen": "screen", "print": "print", "upper_print": "PRINT",
          "string_print": mstring("print"), "cell_screen": mcell("screen")}  # fmt: skip
TYPES = {"all": "all", "Ia": "Ia", "upper_all": "ALL", "string_all": mstring("all"),
         "cell_all": mcell("all")}  # fmt: skip
DENSITY = {"mean": "mean", "median": "median", "upper_mean": "MEAN",
           "string_mean": mstring("mean"), "cell_mean": mcell("mean")}  # fmt: skip

# the leading arguments of each function's calls
VOL = np.arange(120, dtype=np.uint16).reshape(6, 5, 4)
DESKEW = [VOL, 32.45, 0.5, 0.108]
ROTATE = [VOL, 32.45, 2.5]
NO_MOVIES = MStruct((0, 0), ["source", "movieLength"], [])  # a condition without movies
MOVIE = {"source": "movie_missing/", "movieLength": 2.0}
COND = ["cond_missing"]
# a movie whose tracks are missing: the lifetime analysis checks its frame rate, then fails
# reading the tracks (after its inputParser)
LIFETIME_MOVIE = [
    _structs(
        (1, 1),
        ["source", "channels", "framerate", "movieLength"],
        [{"source": "movie_missing/ch1/", "channels": mcell("movie_missing/ch1/",
          "movie_missing/ch2/"), "framerate": 2.0, "movieLength": 60.0}],
    )
]  # fmt: skip
DATA = {
    "no_movies": NO_MOVIES,
    "one_movie": _structs((1, 1), ["source", "movieLength"], [MOVIE]),
    "two_movies": _structs((1, 2), ["source", "movieLength"], [MOVIE, MOVIE]),
    "column_of_movies": _structs((2, 1), ["source", "movieLength"], [MOVIE, MOVIE]),
    "struct_empty": MStruct((0, 0), [], []),
    "empty": None,
    "empty_1x0": mempty("double", (1, 0)),
    "logical_empty": mempty("logical"),
    "uint8_empty": mempty("uint8"),
    "char_empty": "",
    "cell_empty": EMPTY_CELL,
}

# the name-value parameters of each function with a value each takes (its default): text at an
# optional argument may be a name, alone or followed by a value
NAMES = {
    "deskewFrame3D": {"Crop": False, "Interp": "linear"},
    "rotateFrame3D": {"Crop": True, "resample": None, "ObjectiveScan": False, "outSize": None,
                      "Interp": "linear"},
    "XR_runTracking3D": {"Overwrite": False, "FileName": "trackedFeatures", "Frames": None,
                         "DetectionFile": "Detection3D.mat", "ParallelMethod": "lock_distributed"},
    "XR_loadConditionData3D": {"MovieSelector": "Ex", "IgnoreExpNumber": False,
                               "ChannelPatterns": EMPTY_CELL, "Angle": 31.5, "PixelSize": 0.104,
                               "dz": None, "ObjectiveScan": False, "Reversed": False,
                               "FrameRate": None},
}  # fmt: skip

# function, parameter, leading arguments, how the value is given ("pos": the next argument,
# "pair": name and value, "data": the first argument), designed values
PROBED = [
    ("deskewFrame3D", "reverse", DESKEW, "pos", {}),
    ("deskewFrame3D", "Crop", DESKEW, "pair", {}),
    ("deskewFrame3D", "Interp", DESKEW, "pair", INTERP),
    ("rotateFrame3D", "reverse", ROTATE, "pos", {}),
    ("rotateFrame3D", "Crop", ROTATE, "pair", {}),
    ("rotateFrame3D", "resample", ROTATE, "pair", {"ones3": [1.0, 1.0, 1.0]}),
    ("rotateFrame3D", "ObjectiveScan", ROTATE, "pair", {}),
    ("rotateFrame3D", "outSize", ROTATE, "pair", {"size3": [10.0, 12.0, 8.0]}),
    ("rotateFrame3D", "Interp", ROTATE, "pair", INTERP),
    ("XR_runTracking3D", "data", [], "data", {}),
    ("XR_runTracking3D", "settings", [NO_MOVIES], "pos", {}),
    ("XR_runTracking3D", "Overwrite", [NO_MOVIES], "pair", {}),
    ("XR_runTracking3D", "FileName", [NO_MOVIES], "pair", {"default": "trackedFeatures"}),
    ("XR_runTracking3D", "Frames", [NO_MOVIES], "pair", FRAMES),
    ("XR_runTracking3D", "DetectionFile", [NO_MOVIES], "pair", {"default": "Detection3D.mat"}),
    ("XR_runTracking3D", "ResultsPath", [NO_MOVIES], "pair", {"cell_path": mcell("results/")}),
    ("XR_runTracking3D", "ParallelMethod", [NO_MOVIES], "pair", PARALLEL),
    ("XR_runTrackProcessing3D", "data", [], "data", {}),
    ("XR_runTrackProcessing3D", "Buffer", [NO_MOVIES], "pair", PAIRS),
    ("XR_runTrackProcessing3D", "Overwrite", [NO_MOVIES], "pair", {}),
    ("XR_runTrackProcessing3D", "FitMixtures", [NO_MOVIES], "pair", {}),
    ("XR_runTrackProcessing3D", "TrackerOutput", [NO_MOVIES], "pair",
     {"default": "trackedFeatures.mat"}),
    ("XR_runTrackProcessing3D", "FileName", [NO_MOVIES], "pair",
     {"default": "ProcessedTracks.mat"}),
    ("XR_runTrackProcessing3D", "DetectionFile", [NO_MOVIES], "pair",
     {"default": "Detection3D.mat"}),
    ("XR_runTrackProcessing3D", "Frames", [NO_MOVIES], "pair", FRAMES),
    ("XR_runTrackProcessing3D", "Preprocess", [NO_MOVIES], "pair", {}),
    ("XR_runTrackProcessing3D", "Postprocess", [NO_MOVIES], "pair", {}),
    ("XR_runTrackProcessing3D", "BufferAll", [NO_MOVIES], "pair", {}),
    ("XR_runTrackProcessing3D", "ResultsPath", [NO_MOVIES], "pair",
     {"cell_path": mcell("results/")}),
    ("XR_runTrackProcessing3D", "ParallelMethod", [NO_MOVIES], "pair", PARALLEL),
    ("XR_runTrackProcessing3D", "WriteAmira", [NO_MOVIES], "pair", {}),
    ("XR_loadConditionData3D", "condDir", [], "pos", COND_TEXT),
    ("XR_loadConditionData3D", "maxlevel", [*COND, EMPTY_CELL, EMPTY_CELL], "pos", LEVELS),
    ("XR_loadConditionData3D", "MovieSelector", COND, "pair",
     {"default": "Ex", "cell_default": mcell("Ex")}),
    ("XR_loadConditionData3D", "IgnoreExpNumber", COND, "pair", {}),
    ("XR_loadConditionData3D", "ChannelPatterns", COND, "pair", {"cell_pattern": mcell("CamA")}),
    ("XR_loadConditionData3D", "Angle", COND, "pair", {"default": 31.5}),
    ("XR_loadConditionData3D", "PixelSize", COND, "pair", {"default": 0.104}),
    ("XR_loadConditionData3D", "dz", COND, "pair", {}),
    ("XR_loadConditionData3D", "ObjectiveScan", COND, "pair", {}),
    ("XR_loadConditionData3D", "Reversed", COND, "pair", {}),
    ("XR_loadConditionData3D", "FrameRate", COND, "pair", {}),
    ("runLifetimeAnalysis3D", "Display", LIFETIME_MOVIE, "pair", ON_OFF),
    ("runLifetimeAnalysis3D", "DisplayMode", LIFETIME_MOVIE, "pair", SCREEN),
    ("runLifetimeAnalysis3D", "AnalysisPath", LIFETIME_MOVIE, "pair", {"default": "Analysis"}),
    ("runLifetimeAnalysis3D", "ProcessedTracks", LIFETIME_MOVIE, "pair",
     {"default": "ProcessedTracks.mat"}),
    ("runLifetimeAnalysis3D", "LifetimeData", LIFETIME_MOVIE, "pair",
     {"default": "lifetimeData.mat"}),
    ("runLifetimeAnalysis3D", "Type", LIFETIME_MOVIE, "pair", TYPES),
    ("runLifetimeAnalysis3D", "Cutoff_f", LIFETIME_MOVIE, "pair", {"default": 3.0, "five": 5.0}),
    ("runLifetimeAnalysis3D", "Print", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "Overwrite", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "ShowThresholdRange", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "Rescale", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "RemoveOutliers", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "ExcludeVisitors", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "PoolDatasets", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "ShowStatistics", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "SelectIndex", LIFETIME_MOVIE, "pair",
     {"cell_rows": mcell([1.0, 2.0])}),
    ("runLifetimeAnalysis3D", "CorrectObservationBias", LIFETIME_MOVIE, "pair", {}),
    ("runLifetimeAnalysis3D", "InitDensity", LIFETIME_MOVIE, "pair", DENSITY),
    ("runLifetimeAnalysis3D", "MaskPath", LIFETIME_MOVIE, "pair",
     {"default": "Detection/cellmask.tif"}),
    ("runLifetimeAnalysis3D", "Binning", LIFETIME_MOVIE, "pair", {"default": 1.0}),
]  # fmt: skip


def probes(fname: str, param: str, how: str, extra: dict) -> dict[str, list]:
    """{probe name: the arguments after the leading ones} of one parameter's battery"""
    values = {**(DATA if how == "data" else GENERIC), **extra}
    if how == "pair":
        return {k: [param, v] for k, v in values.items()}
    out = {k: [v] for k, v in values.items()}
    if how == "pos":
        names = NAMES.get(fname, {})
        for name, v in names.items():
            out[f"name_{name}"] = [name]
            out[f"pair_{name}"] = [name, v]
        if names:
            first = next(iter(names))
            out["prefix_pair"] = [first[:2], names[first]]
        out["unknown_pair"] = ["foo", 1.0]
    return out


def _arg(x: object) -> MValue:
    return x if isinstance(x, MValue) else to_mvalue(x)


def _case(fname: str, param: str, lead: list, how: str, extra: dict) -> Case:
    battery = probes(fname, param, how, extra)
    names = [f"p_{k}" for k in battery]
    return Case(
        f"{fname}.{param}",
        matlab="pkpy_x27_probe",
        args=[fname, _cell(lead), mcell(*names), _cell([_cell(a) for a in battery.values()])],
        python=port_outcomes,
        helpers=HELPERS,
        files=False,
    )


# ----------------------------------------------------------------------------------------------
# the Python side


def _py(v: object) -> object:
    """A MATLAB value of a case as the MATLAB-name layer takes it: char rows as str, other
    arrays as the registry's MArray (MATLAB's class and size), cells as Cell, string scalars
    as MString, a 1x1 struct as a dict, other structs as StructArray, a handle as a callable;
    paths given by placeholder relative to the work folder."""
    from petakit5d.io.matfile import StructArray

    if isinstance(v, MPlaceholder):
        return (Path.cwd() / (v.name or v.sub)).as_posix()
    if isinstance(v, MArray):
        if v.mclass == "char":
            return v.text()
        flat = v.data.reshape(-1, order="F")
        items = [bool(x) for x in flat] if v.mclass == "logical" else [x.item() for x in flat]
        return mv.MArray(v.mclass, tuple(v.size), items)
    if isinstance(v, MCell):
        return mv.Cell([_py(x) for x in v.items], tuple(v.size))
    if isinstance(v, MString):
        return mv.MString(v.items[0])
    if isinstance(v, MStruct):
        els = [{k: _py(x) for k, x in e.items()} for e in v.items]
        if tuple(v.size) == (1, 1):
            return els[0]
        if not els:
            return StructArray.empty(tuple(v.fields), tuple(v.size))
        return StructArray(tuple(v.size), tuple(v.fields), els)
    if isinstance(v, MFunc):
        return math.sin
    raise TypeError(f"no Python form for {v!r}")


def port_outcomes(fname: str, lead: MCell, names: MCell, battery: MCell) -> dict[str, str]:
    """How the spec's MATLAB-name layer takes the same calls: 'accepted', or the identifier"""
    spec = REGISTRY.find(fname)
    head = [_py(x) for x in lead.items]
    out = {}
    for name, call in zip(names.items, battery.items, strict=True):
        try:
            spec._bind([*head, *(_py(x) for x in call.items)], {})
        except MatlabArgumentError as err:
            out[name.text()] = err.identifier
        else:
            out[name.text()] = "accepted"
    return out


CASES = [_case(*p) for p in PROBED]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_every_check_of_the_port_is_probed(extract):
    """each parameter with a check of the port's own (accepts=) has a case, each case's
    parameter has one, and each is a parameter whose MATLAB validator the extract withholds"""
    checked = {
        (s.matlab_name, p.matlab)
        for s in REGISTRY.specs()
        for p in s.params
        if p.accepts is not None
    }
    assert {(f, p) for f, p, *_ in PROBED} == checked
    for f, p in checked:
        assert any(m.name == p and m.withheld for m in extract[f].params), (f, p)


# ----------------------------------------------------------------------------------------------
# the words of the port's checks


def _spec(accepts: str, *, keep_unmatched: bool = False, optional: bool = False):
    from petakit5d.params import Param, Registry

    reg = Registry()
    params = (Param("x", None, accepts=accepts), Param("Flag", None))
    if optional:

        @reg.register("XR_words_fake", *params, keep_unmatched=keep_unmatched)
        def words_opt(x=None, *, flag=False):
            return x

    else:

        @reg.register("XR_words_fake", *params, keep_unmatched=keep_unmatched)
        def words(x, *, flag=False):
            return x

    return reg.find("XR_words_fake")


def _takes(accepts: str, value: object) -> bool:
    try:
        _spec(accepts).param("x").check(value, "XR_words_fake")
    except MatlabArgumentError:
        return False
    return True


@pytest.mark.parametrize(
    ("accepts", "value", "takes"),
    [
        ("logical", True, True),
        ("logical", mv.MArray("logical", (0, 0), []), True),
        ("logical", 1.0, False),
        ("numeric", mv.MArray("uint8", (1, 1), [3]), True),
        ("numeric", True, False),
        ("char", "abc", True),
        ("char", mv.MString("abc"), False),
        ("string", mv.MString("abc"), True),
        ("string", "abc", False),  # a Python str is MATLAB's char
        ("cell", mv.Cell([], (0, 0)), True),
        ("struct", {"a": 1.0}, True),
        ("handle", math.sin, True),
        ("scalar", "a", True),
        ("scalar", "ab", False),
        ("vector", mv.MArray("double", (1, 0), []), True),  # MATLAB's isvector: n >= 0
        ("vector", mv.MArray("double", (0, 1), []), True),
        ("vector", None, False),  # 0x0
        ("vector", mv.MArray("double", (2, 2), [1.0, 2.0, 3.0, 4.0]), False),
        ("numel(2)", "ab", True),
        ("numel(2)", (1.0, 2.0, 3.0), False),
        ("evenly-spaced", (1.0, 3.0, 5.0), True),
        ("evenly-spaced", (1.0, 2.0, 4.0), False),
        ("evenly-spaced", (5.0, 5.0), True),
        ("evenly-spaced", 5.0, False),
        ("evenly-spaced", "abc", True),  # char codes
        ("evenly-spaced", (True, False), True),
        ("evenly-spaced", (1.0, math.nan, 3.0), False),  # each NaN difference is a value
        ("evenly-spaced", (1.0, math.nan), True),
        ("evenly-spaced", (math.inf, math.inf), True),  # one NaN difference
        # along the first dimension that is not 1: [1 2; 3 4] steps 2 and 2, [1 2; 3 5] 2, 3
        ("evenly-spaced", mv.MArray("double", (2, 2), [1.0, 3.0, 2.0, 4.0]), True),
        ("evenly-spaced", mv.MArray("double", (2, 2), [1.0, 3.0, 2.0, 5.0]), False),
        ("evenly-spaced", mv.MString("abc"), False),
        ("nonnegative-whole", (0.0, 2.0, math.inf), True),
        ("nonnegative-whole", (1.0, 2.5), False),
        ("nonnegative-whole", -1.0, False),
        ("nonnegative-whole", math.nan, False),
        ("nonnegative-whole", True, False),
        ("positive-whole", (1.0, 2.0, math.inf), True),
        ("positive-whole", 0.0, False),
        ("positive-whole", (1.0, 2.5), False),
        ("positive-whole", -2.0, False),
        ("positive-whole", True, False),
        ("positive-whole", "a", False),
        ("in(cubic,linear)", "LINEAR", True),
        ("in(cubic,linear)", mv.MArray("char", (1, 6), list("linear")), True),
        ("evenly-spaced", mv.MArray("char", (1, 3), list("abc")), True),
        ("in(cubic,linear)", mv.MString("cubic"), True),
        ("in(cubic,linear)", "nearest", False),
        ("in(cubic,linear)", mv.Cell(["linear", "Cubic"]), True),  # every element
        ("in(cubic,linear)", mv.Cell(["linear", "x"]), False),
        ("in(cubic,linear)", mv.Cell([], (0, 0)), True),
        ("in(cubic,linear)", 1.0, False),
        ("char & !in(FrameRate,MovieSelector)", "movieselector", False),
        ("char & !in(FrameRate,MovieSelector)", "MovieSel", True),
        ("char & !in(FrameRate,MovieSelector)", 1.0, False),
        ("!scalar", (1.0, 2.0), True),
        ("!scalar", 1.0, False),
        ("empty | numeric & !scalar | numeric & nonnegative-whole", "", True),
        ("empty | numeric & !scalar | numeric & nonnegative-whole", "a", False),
        ("empty | numeric & !scalar | numeric & nonnegative-whole", 3.0, True),
        ("empty | numeric & !scalar | numeric & nonnegative-whole", 2.5, False),
        ("empty | numeric & !scalar | numeric & nonnegative-whole", (-1.0, 2.5), True),
    ],
)
def test_the_words_of_the_ports_checks(accepts, value, takes):
    assert _takes(accepts, value) is takes


def test_a_check_of_the_port_and_its_use():
    from petakit5d.params import Param, Registry
    from petakit5d.params.spec import RegistryError, accepts_alternatives, accepts_mcc

    assert accepts_alternatives("char | numeric & scalar") == [
        frozenset({"ischar"}),
        frozenset({"isnumeric", "isscalar"}),
    ]
    with pytest.raises(RegistryError, match="unknown word 'nonsense'"):
        accepts_alternatives("char | nonsense")
    # the command-line conversion follows the generator's rule on the classes it takes
    assert accepts_mcc(accepts_alternatives("logical")) == "str2num"
    assert accepts_mcc(accepts_alternatives("char")) == "none"
    assert accepts_mcc(accepts_alternatives("in(a,b)")) == "none"
    assert accepts_mcc(accepts_alternatives("cell")) == "eval:{"
    assert accepts_mcc(accepts_alternatives("numel(2)")) == "str2num"
    assert accepts_mcc(accepts_alternatives("char & !in(a,b)")) == "none"
    assert accepts_mcc(accepts_alternatives("empty | numeric & nonnegative-whole")) == "str2num"
    assert accepts_mcc(accepts_alternatives("numeric & positive-whole")) == "str2num"
    assert accepts_alternatives("char & !in(a, b)") == [
        frozenset({"ischar", "port:not:port:in:a,b"})
    ]
    spec = _spec("logical")
    assert spec.param("x").mcc == "str2num" and spec.param("x").cell_hint() is None
    assert _spec("cell").param("x").cell_hint() == "cell"
    with pytest.raises(MatlabArgumentError, match="must pass the port's check: logical"):
        spec.bind_matlab([1.0])
    reg = Registry()
    with pytest.raises(RegistryError, match="a validator text and the port's own check"):

        @reg.register("XR_both_fake", Param("x", "@islogical", accepts="logical"))
        def both(x):
            return x

        reg.find("XR_both_fake").params  # noqa: B018 - resolves the spec


def test_with_keep_unmatched_an_optional_gives_way_to_any_name_or_struct():
    """MATLAB's inputParser with KeepUnmatched: text that the optional argument refuses and
    that names no parameter, followed by a value, or a struct of any fields, starts the
    name-value part (what matches no parameter is kept apart); without KeepUnmatched the
    optional refuses them"""
    keep = _spec("numeric", keep_unmatched=True, optional=True)
    assert keep._bind(["foo", 1.0], {}) == {}
    assert keep._bind([{"zeta": 1.0}], {}) == {}
    assert keep._bind([{"flag": True}], {}) == {"flag": True}
    assert keep._bind(["Flag", True], {}) == {"flag": True}
    strict = _spec("numeric", optional=True)
    for args in (["foo", 1.0], [{"zeta": 1.0}]):
        with pytest.raises(MatlabArgumentError) as err:
            strict._bind(args, {})
        assert err.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
    assert strict._bind(["Flag", True], {}) == {"flag": True}


def test_third_party_code_takes_the_ports_checks_only(extract):
    """check_spec: a validator text where MATLAB's is withheld is a problem, and so is a check
    of the port's own where MATLAB's validator is compared"""
    from petakit5d.params import Param, Registry
    from petakit5d.params.check import check_spec
    from petakit5d.params.matlab_source import MatlabFunction, MatlabParam

    def fn(withheld):
        v = None if withheld else "@islogical"
        return MatlabFunction(
            "XR_fake", "pointDetection/llsmtools/XR_fake.m" if withheld else "x/XR_fake.m",
            ["x", "varargin"],
            [MatlabParam("required", "x", None, v, withheld)],
        )  # fmt: skip

    def make(validator, accepts):
        reg = Registry()

        @reg.register("XR_fake", Param("x", validator, accepts=accepts))
        def fake(x):
            return x

        return reg.find("XR_fake")

    assert check_spec(make(None, "logical"), fn(True)) == []
    assert check_spec(make(None, None), fn(True)) == []  # interface only
    problems = check_spec(make("@islogical", None), fn(True))
    assert any("instead of a validator" in p for p in problems), problems
    problems = check_spec(make(None, "logical"), fn(False))
    assert any("accepts= is the port's own check" in p for p in problems), problems
