"""Setting files of a microscope acquisition (``<name>_Settings.txt``).

=================================  ==============================
``parse_setting_files``            XR_parseSettingFiles
``parse_image_setting_files``      XR_parseSettingFiles_wrapper
=================================  ==============================

A setting file gives the X stage (``X Stage Offset, ...``, stage scans) or XZ stage
(``XZ Stage Offset, ...``, sample scans with the objective) offset, interval in um and the
number of excitation pixels on one tab-separated line; the sign of the interval decides
whether the z stack is flipped. Results are dicts with MATLAB's field names (``stageOff``,
``StageInterval``, ``numExiciation``), values as ``str2double`` reads them.
"""

from __future__ import annotations

import os
import re
from collections.abc import Sequence

from petakit5d.core.dtypes import NAN
from petakit5d.core.paths import fileparts, matlab_sort_key
from petakit5d.core.text import sprintf, str2double
from petakit5d.pipeline.filenames import (
    MatlabError,
    _warn,
    matlab_regex,
    parse_filename_fields,
    regexp_names,
)

__all__ = ["parse_image_setting_files", "parse_setting_files"]

FIELDS = ("stageOff", "StageInterval", "numExiciation")

# XR_parseSettingFiles's expressions, as MATLAB writes them (regexpi; translated by
# matlab_regex: '.' also matches a newline or a tab)
STAGE_LINE = (
    r"\n(XZ? +Stage Offset, Interval \(um\), # of Pixels for Excitation \(0\) :"
    r"\t\d\t-?\d*?.?\d+\t\d+)\n"
)
STAGE_FIELDS = r"(?<stageOff>\d+)\t(?<interval>-?\d*?.?\d+)\t(?<numExiciation>\d+)"
MAX_LINES = 100
_EOL = re.compile(r"\r\n|\r|\n")  # fgetl's line ends


def _fgetl_lines(fn: str) -> list[str]:
    """The lines parseSingleSettingFile keeps: ``fgetl`` until the end of the file (-1) or the
    100th line; fgetl ends a line at ``\\n``, ``\\r\\n`` or ``\\r`` and returns it without its end.
    """
    try:
        if os.path.isdir(fn):
            raise OSError(fn)
        with open(fn, "rb") as fh:
            data = fh.read()
    except OSError:
        raise MatlabError(
            "Invalid file identifier. Use fopen to generate a valid file identifier.",
            "MATLAB:FileIO:InvalidFid",
        ) from None
    text = data.decode("utf-8", errors="replace")
    lines: list[str] = []
    pos, n = 0, len(text)
    while len(lines) < MAX_LINES and pos < n:
        m = _EOL.search(text, pos)
        if m is None:
            lines.append(text[pos:])
            break
        lines.append(text[pos : m.start()])
        pos = m.end()
    return lines


def _parse_single_setting_file(fn: str) -> dict:
    # each line ends with a newline (the stage line can also be the last one)
    allline = "\n".join(_fgetl_lines(fn)) + "\n"
    m = matlab_regex(STAGE_LINE, True).search(allline)
    if m is None:
        # str{1}{1} of no match
        raise MatlabError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    tmp = regexp_names(m.group(1), matlab_regex(STAGE_FIELDS, True))
    if not tmp:
        raise MatlabError("Not enough input arguments.", "MATLAB:minrhs")
    if len(tmp) > 1:
        # str2double of a comma-separated list of several matches
        raise MatlabError("Too many input arguments.", "MATLAB:TooManyInputs")
    t = tmp[0]
    return {
        "stageOff": str2double(t["stageOff"]),
        "StageInterval": str2double(t["interval"]),
        "numExiciation": str2double(t["numExiciation"]),
    }


def parse_setting_files(setting_filenames: Sequence[str]) -> list[dict]:
    """The stage offset, interval and excitation pixels of each setting file.

    One dict per file. As in MATLAB: a missing file is skipped with a warning that names
    it; its entry gets None values when a later file is read (and no entry when none is).
    A file without the stage line fails.
    """
    from petakit5d.jobs.files import batch_file_exist

    setting_filenames = list(setting_filenames)
    exists = batch_file_exist(setting_filenames)
    info: list[dict] = []
    for i, fn in enumerate(setting_filenames):
        if not exists[i]:
            _warn(f"Setting file {fn} does not exist, skip it.")
            continue
        s = _parse_single_setting_file(fn)
        while len(info) <= i:
            info.append(dict.fromkeys(FIELDS))
        info[i].update(s)
    return info


NAN_INFO = {"stageOff": NAN, "StageInterval": NAN, "numExiciation": NAN}


def parse_image_setting_files(image_filenames: Sequence[str]) -> list[dict]:
    """The setting file values of each image file (XR_parseSettingFiles_wrapper).

    Partial volume files (``_part<digits>.tif``) are left out first. An image's setting file
    is ``<name>_Settings.txt``, else ``<prefix>Scan_Iter_<iter>_<x>x_<y>y_<z>z_<t>t_Settings.txt``,
    else ``<prefix>Scan_Iter_<iter>_Settings.txt``, from the fields of its tile name (no such
    files for a name that matches no grammar). Each setting file is read once. Images
    without a setting file get NaN values, with one warning that names them.
    """
    from petakit5d.jobs.files import batch_file_exist, exist

    rx = matlab_regex(r"^(?!.*_part[0-9]*.tif).*$")
    image_filenames = [f for f in image_filenames if rx.search(f)]
    parts = [fileparts(f) for f in image_filenames]
    pathstrs = [p[0] for p in parts]
    fsnames = [p[1] for p in parts]
    sfns = [f"{p}{os.sep}{s}_Settings.txt" for p, s in zip(pathstrs, fsnames, strict=True)]
    sfns_exist = batch_file_exist(sfns)

    tmp: list[dict[str, str] | None] = []
    sfns_1: list[str] = []
    sfns_1_exist: list[bool] = []
    if not all(sfns_exist):
        tmp = parse_filename_fields(fsnames)
        sfns_1 = [
            ""
            if t is None
            else sprintf(
                "%s/%sScan_Iter_%s_%sx_%sy_%sz_%st_Settings.txt",
                pathstrs[x],
                *(t[k] for k in ("prefix", "Iter", "x", "y", "z", "t")),
            )
            for x, t in enumerate(tmp)
        ]
        matched = [x for x, t in enumerate(tmp) if t is not None]
        sfns_1_exist = [False] * len(fsnames)
        for x, e in zip(matched, batch_file_exist([sfns_1[x] for x in matched]), strict=True):
            sfns_1_exist[x] = e

    sfn_cell: list[str | None] = [None] * len(image_filenames)
    for i in range(len(image_filenames)):
        if sfns_exist[i]:
            sfn_cell[i] = sfns[i]
            continue
        if sfns_1_exist[i]:
            sfn_cell[i] = sfns_1[i]
            continue
        t = tmp[i]
        if t is None:
            continue
        sfn_2 = sprintf("%s/%sScan_Iter_%s_Settings.txt", pathstrs[i], t["prefix"], t["Iter"])
        if exist(sfn_2):
            sfn_cell[i] = sfn_2

    missing = [f for f, s in zip(image_filenames, sfn_cell, strict=True) if s is None]
    if missing:
        _warn(f"No setting file found for {', '.join(missing)}, skip it.")
    uniq = sorted({s for s in sfn_cell if s is not None}, key=matlab_sort_key)
    index = {s: k for k, s in enumerate(uniq)}
    uniq_info = parse_setting_files(uniq)
    out = []
    for s in sfn_cell:
        if s is None:
            out.append(dict(NAN_INFO))
            continue
        k = index[s]
        if k >= len(uniq_info):
            raise MatlabError(
                "Index exceeds the number of array elements. Index must not exceed "
                f"{len(uniq_info)}.",
                "MATLAB:badsubscript",
            )
        out.append(dict(uniq_info[k]))
    return out
