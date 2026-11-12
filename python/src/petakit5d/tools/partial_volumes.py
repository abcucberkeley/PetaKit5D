"""Grouping of partial-volume files (``groupPartialVolumeFiles``).

A volume saved in parts is ``<name>.tif``, ``<name>_part0001.tif``, ``<name>_part0002.tif``,
...; each group is the main file and its parts in sorted order.

The main names are matched as text (MATLAB since b24c975), and a channel pattern that is not a
valid regular expression is matched as text only (``channel_regexp_pattern``). MATLAB's
quirks are kept: the extension is used as a regular expression, the part pattern is not
anchored at the start (the parts of ``ab`` also go to ``b``), parts without a main file are
left out, the dates and sizes of the parts are in listing order while their names are sorted,
the data size of tiff files is the largest of the first and last file's image sizes and the
other files' bytes, and without parts the dates of a folder listing are a row (a column for a
list of files).
"""

from __future__ import annotations

import datetime
import os
import re
import time
import warnings
from collections.abc import Sequence
from typing import Any

from petakit5d.params import Param, register

__all__ = ["PartialVolumeWarning", "dir_datenum", "group_partial_volume_files"]


class PartialVolumeWarning(UserWarning):
    """MATLAB's warning without an identifier."""

    identifier = ""


def dir_datenum(mtime: float) -> float:
    """``dir``'s datenum of a modification time: local time, whole seconds."""
    lt = time.localtime(int(mtime // 1))
    days = datetime.date(lt.tm_year, lt.tm_mon, lt.tm_mday).toordinal() + 366
    return days + (lt.tm_hour * 3600 + lt.tm_min * 60 + lt.tm_sec) / 86400


def _stat(path: str) -> tuple[float, float]:
    # dir's datenum and bytes (0 for a folder)
    st = os.stat(path)
    size = 0.0 if os.path.isdir(path) else float(st.st_size)
    return dir_datenum(st.st_mtime), size


def _has_match(pattern: str, names: Sequence[str]) -> list[bool]:
    # ~cellfun(@isempty, regexpi(names, pattern))
    rx = re.compile(pattern, re.IGNORECASE)
    return [rx.search(n) is not None for n in names]


def _dbytes(dtype: str) -> int:
    return {"uint8": 1, "uint16": 2, "single": 4}.get(dtype, 8)


def _if_text(x: Any) -> bool:
    # MATLAB's if of a char (only char passes the validator) or of the default false
    if isinstance(x, str):
        return len(x) > 0 and "\0" not in x
    return bool(x)


@register(
    "groupPartialVolumeFiles",
    Param("dataPath", "@(x) isempty(x) || ischar(x)", "Folder of the files (or empty)."),
    Param(
        "fileFullpathList",
        "@(x) isempty(x) || iscell(x)",
        "Files of one folder, used when dataPath is empty.",
    ),
    Param("ext", "@(x) ischar(x)", "File extension."),
    Param(
        "onlyFirstTP",
        "@(x) ischar(x)",
        "Only the first time point (Iter_0000_): any nonempty text.",
        differs="MATLAB's default false fails its own validator @ischar; '' acts the same",
    ),
    Param("channelPatterns", "@(x) iscell(x)", "Keep files that match one of these."),
    task="T05",
)
def group_partial_volume_files(
    data_path: str | None,
    file_fullpath_list: Sequence[str] = (),
    *,
    ext: str = ".tif",
    only_first_tp: str = "",
    channel_patterns: Sequence[str] = (),
) -> tuple[Any, list, Any, Any]:
    """Check and group file names that are parts of the same volume.

    Returns ``(contain_part_volume, grouped_fnames, grouped_datenum, grouped_datasize)``:
    a flag per file (False for no files); without parts one list per file and the dates and
    sizes as arrays; with parts one list per main file and an array of dates and of sizes
    per group.
    """
    import numpy as np

    from petakit5d.core.paths import fileparts, filesep, matlab_dir, matlab_sort_key

    print("Check partial volume files... ", end="", flush=True)
    files = list(file_fullpath_list or ())
    if not data_path:
        if _if_text(only_first_tp):
            files = [f for f in files if "Iter_0000_" in f]
        data_path, _, ext = fileparts(files[0])
        n = len(files)
        fnames, datenum, datasize = [], np.zeros(n), np.zeros(n)
        for i, f in enumerate(files):
            entries = matlab_dir(f)
            fnames.append(entries[0].name)
            datenum[i], datasize[i] = _stat(f)
        full = files
    else:
        expression = "*_Iter_0000_*" if _if_text(only_first_tp) else "*"
        entries = matlab_dir(f"{data_path}{filesep}{expression}{ext}")
        fnames = [e.name for e in entries]
        stats = [_stat(os.path.join(e.folder, e.name)) for e in entries]
        datenum = np.array([s[0] for s in stats], dtype=np.float64)
        datasize = np.array([s[1] for s in stats], dtype=np.float64)

        if channel_patterns:
            fns = [f"{data_path}{filesep}{x}" for x in fnames]
            keep = [False] * len(fns)
            for c in channel_patterns:
                try:
                    rx = re.compile(c)
                except re.error:  # matched as text (channel_regexp_pattern)
                    rx = re.compile(re.escape(c))
                keep = [
                    k or c in f or rx.search(f) is not None for k, f in zip(keep, fns, strict=True)
                ]
            fnames = [f for f, k in zip(fnames, keep, strict=True) if k]
            datenum = datenum[np.array(keep, dtype=bool)]
            datasize = datasize[np.array(keep, dtype=bool)]

        if not fnames:
            warnings.warn(
                "The input image list is empty or none of the images matches channelPatterns.",
                PartialVolumeWarning,
                stacklevel=2,
            )
            return False, [], np.zeros(0), np.zeros(0)

        from petakit5d.io.imagesize import get_image_dtype, get_image_size

        dbytes = _dbytes(get_image_dtype(f"{data_path}{filesep}{fnames[0]}"))
        # tiff: only the first and last file, and the largest size for all
        if ext in (".tif", ".tiff"):
            for f in (0, len(fnames) - 1):
                fn = f"{data_path}{filesep}{fnames[f]}"
                datasize[f] = _prod(get_image_size(fn)) * dbytes
            datasize[:] = datasize.max()
        elif ext == ".zarr":
            for f in range(len(fnames)):
                fn = f"{data_path}{filesep}{fnames[f]}"
                try:
                    datasize[f] = _prod(get_image_size(fn)) * dbytes
                except Exception as err:  # MATLAB displays the error and goes on
                    print(err)
        full = [f"{data_path}{filesep}{x}" for x in fnames]

    contain = np.array(_has_match(f"_part\\d+{ext}$", full), dtype=bool)
    # without partial volume files the lists come back in their original form
    if not contain.any():
        print("\nNo partial volume files found.")
        return contain, [[x] for x in fnames], datenum, datasize

    fsnames = [fileparts(f)[1] for f in full]
    main = [not m for m in _has_match(r"_part\d+$", fsnames)]
    grouped_fnames, grouped_datenum, grouped_datasize = [], [], []
    for i in [j for j, m in enumerate(main) if m]:
        mstr = fsnames[i]
        parts = [j for j, m in enumerate(_has_match(f"{re.escape(mstr)}_part\\d+$", fsnames)) if m]
        if parts:
            names = sorted((fsnames[j] for j in parts), key=matlab_sort_key)
            grouped_fnames.append([f"{x}{ext}" for x in [mstr, *names]])
            grouped_datenum.append(np.array([datenum[i], *datenum[parts]]))
            grouped_datasize.append(np.array([datasize[i], *datasize[parts]]))
        else:
            grouped_fnames.append([f"{mstr}{ext}"])
            grouped_datenum.append(np.array([datenum[i]]))
            grouped_datasize.append(np.array([datasize[i]]))

    print("Done!")
    return contain, grouped_fnames, grouped_datenum, grouped_datasize


def _prod(sz: Sequence[Any]) -> float:
    p = 1.0
    for s in sz:
        p = p * float(s)
    return p
