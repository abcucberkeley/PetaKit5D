"""Paths as MATLAB handles them: ``fileparts``, ``dir`` order, ``simplifyPath``, ``dir_recursive``.

- ``fileparts``: MATLAB's split (``'.bashrc'`` is name ``''`` and extension ``'.bashrc'``,
  ``'/a'`` has folder ``'/'``), not ``os.path``'s.
- ``matlab_dir``: the entries of ``dir(path)`` in MATLAB's order (names sorted by UTF-16 code
  units, so uppercase before lowercase), with ``.`` and ``..`` for a folder.
- ``simplify_path``: ``simplifyPath``: an absolute path, folders without a trailing separator.
- ``dir_recursive``: every file (or folder) below a folder, sorted as MATLAB's ``unique`` sorts.
- ``mkdir_recursive``: ``mkdir_recursive`` with the group-write option.

MATLAB writes folder paths with a trailing ``/`` in many parameters (``[dataPath, '/']``);
``with_trailing_sep`` and ``strip_trailing_sep`` make those forms explicit.
"""

from __future__ import annotations

import dataclasses
import fnmatch
import os
import stat
import sys

__all__ = [
    "DirEntry",
    "dir_recursive",
    "fileparts",
    "filesep",
    "matlab_dir",
    "matlab_sort_key",
    "mkdir_recursive",
    "simplify_path",
    "strip_trailing_sep",
    "with_trailing_sep",
]

_IS_PC = sys.platform.startswith("win")
filesep = "\\" if _IS_PC else "/"


def matlab_sort_key(s: str) -> tuple[int, ...]:
    """MATLAB orders text by UTF-16 code units (``sort`` and ``unique`` of cellstr)."""
    b = s.encode("utf-16-be", errors="surrogatepass")
    return tuple(int.from_bytes(b[i : i + 2], "big") for i in range(0, len(b), 2))


def strip_trailing_sep(path: str) -> str:
    """The path without trailing separators (the root stays ``/``)."""
    seps = "/\\" if _IS_PC else "/"
    out = path.rstrip(seps)
    return out if out else path[:1]


def with_trailing_sep(path: str) -> str:
    """The path with exactly one trailing ``/`` (MATLAB's ``[folder, '/']`` form)."""
    return strip_trailing_sep(path) + "/"


# -- fileparts ---------------------------------------------------------------------------------


def _split_ext(name: str) -> tuple[str, str]:
    i = name.rfind(".")
    if i < 0:
        return name, ""
    return name[:i], name[i:]


def _fileparts_unix(file: str) -> tuple[str, str, str]:
    i = file.rfind("/")
    if i < 0:
        path, name = "", file
    else:
        path, name = file[:i], file[i + 1 :]
        if path == "":
            path = "/"
    name, ext = _split_ext(name)
    return path, name, ext


def _fileparts_pc(file: str) -> tuple[str, str, str]:
    # MATLAB's legacy scalar path on Windows (filesep is '\\' there)
    filesep = "\\"
    ext = ""
    path = ""
    ind = max(file.rfind("/"), file.rfind("\\"))
    if ind < 0:
        c = file.rfind(":")
        if c >= 0:
            path = file[: c + 1]
            ind = c
        else:
            ind = None
    else:
        if ind == 1 and file[0] in "/\\":
            path = file
            ind = len(file) - 1
        else:
            path = file[:ind]
    if ind is None:
        name = file
    else:
        if path and path[-1] == ":":
            if len(path) > 2:
                path = path + filesep
            elif len(file) >= 3 and file[2] in "/\\":
                path = path + file[2]
        elif path.rstrip(" \0") == "":
            path = filesep
        name = file[ind + 1 :]
    if name:
        name, ext = _split_ext(name)
    return path, name, ext


def fileparts(file: str | os.PathLike) -> tuple[str, str, str]:
    """MATLAB's ``[pathstr, name, ext] = fileparts(file)``."""
    file = os.fspath(file)
    if file == "":
        return "", "", ""
    return _fileparts_pc(file) if _IS_PC else _fileparts_unix(file)


# -- dir ---------------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class DirEntry:
    """One entry of MATLAB's ``dir``: ``name``, ``folder`` (absolute) and ``isdir``."""

    name: str
    folder: str
    isdir: bool


def _abs_folder(path: str) -> str:
    # MATLAB's dir gives the absolute folder, with '.' and '..' resolved and no trailing
    # separator (symbolic links kept)
    return os.path.normpath(os.path.abspath(path))


def _isdir(full: str) -> bool:
    try:
        return stat.S_ISDIR(os.stat(full).st_mode)
    except OSError:
        return False


def _exists(full: str) -> bool:
    return os.path.lexists(full)


def matlab_dir(path: str | os.PathLike = ".") -> list[DirEntry]:
    """The entries ``dir(path)`` lists, in MATLAB's order.

    A folder lists ``.``, ``..`` and its entries; a file lists itself; a name with ``*``
    lists the matching entries of its folder; anything else lists nothing.
    """
    path = os.fspath(path)
    if path == "":
        path = "."
    if any(c in path for c in "*?") and not os.path.lexists(path):
        folder, pattern = os.path.split(path)
        folder = folder or "."
        try:
            names = os.listdir(folder)
        except OSError:
            return []
        # MATLAB's wildcard is '*' only
        names = [n for n in names if fnmatch.fnmatchcase(n, pattern.replace("?", "[?]"))]
        af = _abs_folder(folder)
        return [
            DirEntry(n, af, _isdir(os.path.join(folder, n)))
            for n in sorted(names, key=matlab_sort_key)
        ]
    if _isdir(path):
        af = _abs_folder(path)
        try:
            names = os.listdir(path) + [".", ".."]
        except OSError:
            return []
        # '.' and '..' sort with the other names (' sp' and '#h' come first)
        return [
            DirEntry(n, af, n in (".", "..") or _isdir(os.path.join(path, n)))
            for n in sorted(names, key=matlab_sort_key)
        ]
    if _exists(path):
        # a file (also a broken link): its folder with links resolved, its name as given
        parent, name = os.path.split(_abs_folder(path))
        return [DirEntry(name, os.path.realpath(parent), False)]
    return []


def simplify_path(path: str | os.PathLike) -> str:
    """``simplifyPath``: the absolute path; a folder has no trailing separator.

    Kept from MATLAB: a pattern that matches two or more entries gives their folder, and a path
    that does not exist gives ``filesep``.
    """
    entries = matlab_dir(path)
    if len(entries) >= 2:
        return entries[0].folder
    if not entries:
        return filesep
    return entries[0].folder + filesep + entries[0].name


def dir_recursive(
    given_dir: str | os.PathLike, fmode: str = "file", use_parpool: bool = False
) -> list[str]:
    """``dir_recursive``: the files (``fmode='file'``) or folders (``'folder'``) below a folder.

    Paths are ``<given_dir>/<...>`` (``given_dir`` as given, without a trailing separator),
    sorted as MATLAB's ``unique`` sorts them. Links to folders are followed, as ``dir`` does.
    ``use_parpool`` is accepted and ignored (the result is the same).
    """
    given_dir = os.fspath(given_dir)
    check_dir = fmode == "folder"
    if given_dir.endswith(filesep):
        given_dir = given_dir[:-1]
    out: list[str] = []
    stack: list[str] = [given_dir] if _isdir(given_dir) else []
    while stack:
        current = stack.pop()
        try:
            names = sorted(os.listdir(current), key=matlab_sort_key)
        except OSError:
            continue
        isdir = [_isdir(os.path.join(current, n)) for n in names]
        subdirs = [current + filesep + n for n, d in zip(names, isdir, strict=True) if d]
        stack.extend(subdirs)
        if check_dir:
            out.extend(subdirs)
        else:
            out.extend(current + filesep + n for n, d in zip(names, isdir, strict=True) if not d)
    return sorted(set(out), key=matlab_sort_key)


def mkdir_recursive(full_path: str | os.PathLike, group_write: bool = False) -> int:
    """``mkdir_recursive``: make every missing folder of the path; with ``group_write`` each
    folder it makes gets group write permission. Returns 1, as MATLAB does."""
    full_path = os.fspath(full_path)
    cur = ""
    # (a relative path from the current folder, not from the root)
    if not _IS_PC and full_path and full_path[0] not in "/~":
        cur = "."
    for part in full_path.split("/"):
        if part == "~":
            cur = "~"
            continue
        # strjoin({cur, part}, '/')
        cur = part if (_IS_PC and cur == "") else cur + "/" + part
        target = os.path.expanduser(cur)
        if os.path.isdir(target):
            continue
        try:
            # makedirs: on Windows a path with backslashes is one part here, and MATLAB's
            # mkdir makes the missing folders above it too
            os.makedirs(target)
        except FileExistsError:
            # another process made it after the check (parallel jobs make the same folders);
            # MATLAB's mkdir doesn't fail on a folder that exists
            if not os.path.isdir(target):
                raise
            continue
        if group_write and not _IS_PC:
            os.chmod(target, os.stat(target).st_mode | stat.S_IWGRP)
    return 1
