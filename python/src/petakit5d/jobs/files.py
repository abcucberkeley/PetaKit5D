"""File checks of the job framework: batch_file_exist and wait_for_files.

batch_file_exist gives the result of MATLAB's ``exist(path, 'file') > 0`` (files and folders;
trailing separators don't matter) for many paths at once. On Linux each parent folder is
opened (or, for more than 100 paths in it, listed) before its paths are checked by name: on
NFS a check by name can miss a file another node made for the attribute cache time (30 to
60 s), while a file in a folder opened just before is seen at once. Elsewhere the paths are
checked one by one, as in MATLAB.

Differences from MATLAB that come from the language: a relative path is relative to the
current folder (MATLAB also finds names on its path), and ``~`` is the home folder as for
MATLAB's exist.
"""

from __future__ import annotations

import os
import re
import sys
import time
from collections.abc import Sequence

# more paths than this in a folder: list it instead of checking each path (MATLAB's value)
MAX_NAME_CHECK_NUM = 100

_IS_LINUX = sys.platform.startswith("linux")
_SEPS = "/\\" if os.name == "nt" else "/"

# wait_for_files' clock and pause (tests replace them)
_pause = time.sleep
_clock = time.monotonic


def fileparts(path: str) -> tuple[str, str, str]:
    """MATLAB's fileparts: (folder, name, extension). The folder of ``'/a'`` is ``'/'``; the
    extension starts at the last dot of the name (``'.hidden'`` is all extension)."""
    i = max(path.rfind(s) for s in _SEPS)
    if i < 0:
        folder, rest = "", path
        if os.name == "nt" and len(path) >= 2 and path[1] == ":":
            folder, rest = path[:2], path[2:]
    else:
        folder, rest = path[:i], path[i + 1 :]
        if i == 0 or (os.name == "nt" and i == 2 and path[1] == ":"):
            folder = path[: i + 1]
    j = rest.rfind(".")
    if j < 0:
        return folder, rest, ""
    return folder, rest[:j], rest[j:]


def exist(path: str) -> bool:
    """``exist(path, 'file') > 0`` of MATLAB for one path: a file or a folder (symbolic links
    followed), trailing separators ignored, ``~`` expanded."""
    if not path:
        return False
    if path.startswith("~"):
        path = os.path.expanduser(path)
    stripped = path.rstrip(_SEPS)
    return os.path.exists(stripped or path[0])


def batch_file_exist(
    file_fullpaths: Sequence[str] | Sequence[Sequence[str]],
    out_fullpath: str | None = None,
    use_parpool: bool = False,
    use_java: bool = True,
) -> list[bool]:
    """Whether each path exists as a file or folder (MATLAB's ``exist(..., 'file') > 0``).

    file_fullpaths : paths; for several tasks with several inputs each (lists), the first
                     path of each is checked, as in MATLAB.
    out_fullpath   : with more than one path, a .mat file the result is saved to
                     (``file_exist_mat`` and ``fileFullpaths``). MATLAB's known issue is
                     kept: if that file exists already, the call fails.
    use_parpool, use_java : accepted for MATLAB's signature; they change how MATLAB checks,
                     not the result.
    """
    del use_parpool, use_java
    paths = list(file_fullpaths)
    if not paths:
        return []
    if len(paths) == 1:
        if not isinstance(paths[0], (str, os.PathLike)):
            # MATLAB's unique() fails on a cell in a cell here
            raise TypeError("batch_file_exist: one entry must be a path, not a list of paths")
        paths = [os.fspath(paths[0])]
        return _folder_file_exist(paths) if _IS_LINUX else [exist(paths[0])]

    if out_fullpath and exist(out_fullpath):
        # KNOWN ISSUE of MATLAB (BatchFileExistTest/existingOutputFileReturnsCachedResult):
        # it returns before it assigns its output, so the call errors
        raise RuntimeError(
            f'Output argument "file_exist_mat" (and possibly others) not assigned: the '
            f"result file {out_fullpath} exists already"
        )

    if not isinstance(paths[0], (str, os.PathLike)):
        paths = [p[0] for p in paths]
    paths = [os.fspath(p) for p in paths]

    flags = _folder_file_exist(paths) if _IS_LINUX else [exist(p) for p in paths]

    if out_fullpath:
        _save_exist_mat(out_fullpath, flags, list(file_fullpaths))
    return flags


def _folder_file_exist(paths: list[str]) -> list[bool]:
    # each unique path once. Absolute paths are grouped by folder: a folder with few paths
    # is opened and the paths checked by name; a larger one is listed and the paths whose
    # names the listing has (ignoring case, for case-insensitive file systems) are checked
    # by name (a listing also shows broken symlinks), as are names a listing may spell
    # differently (non-ASCII, short names with ~). Other paths (relative, ~, '/', '.',
    # '..') are checked one by one.
    uniq = sorted(set(paths))
    found = dict.fromkeys(uniq, False)
    by_folder: dict[str, list[tuple[str, str, str]]] = {}
    for p in uniq:
        full = re.sub("/+$", "", p)
        folder, name, ext = fileparts(full)
        fn = name + ext
        if full.startswith("/") and fn not in ("", ".", ".."):
            by_folder.setdefault(folder, []).append((p, full, fn))
        else:
            found[p] = exist(p)

    for folder, items in by_folder.items():
        if not os.path.isdir(folder):
            _open_parent_folders(folder)
            if not os.path.isdir(folder):
                continue
        if len(items) <= MAX_NAME_CHECK_NUM:
            _open_folder(folder)
        else:
            listed = _list_folder(folder)
            if listed is not None:
                lower = {x.lower() for x in listed}
                items = [it for it in items if it[2].lower() in lower or _odd_name(it[2])]
        for p, full, _ in items:
            found[p] = os.path.exists(full)
    return [found[p] for p in paths]


def _odd_name(name: str) -> bool:
    return "~" in name or any(ord(c) > 127 for c in name)


def _open_folder(folder: str) -> None:
    # opening a folder revalidates its cached attributes on NFS, so a check by name right
    # after it sees the files other nodes made in it
    try:
        with os.scandir(folder):
            pass
    except OSError:
        pass


def _open_parent_folders(folder: str) -> None:
    # a folder another node made is hidden to a check by name until its parent is
    # revalidated: open the parents, from the first one that exists
    parent = fileparts(folder)[0]
    if not parent or parent == folder:
        return
    if not os.path.isdir(parent):
        _open_parent_folders(parent)
    _open_folder(parent)


def _list_folder(folder: str) -> list[str] | None:
    # names in a folder, None if it can't be listed (e.g. no read permission)
    try:
        return os.listdir(folder)
    except OSError:
        return None


def _save_exist_mat(out_fullpath: str, flags: list[bool], paths: list) -> None:
    # as MATLAB: -v7.3 to <out without its 4-character extension>_uuid_<uuid>.mat, then a
    # rename
    import uuid

    from petakit5d.io import matfile  # F07

    tmp = f"{out_fullpath[:-4]}_uuid_{uuid.uuid4()}.mat"
    matfile.write_mat(tmp, {"file_exist_mat": flags, "fileFullpaths": paths}, version="7.3")
    os.replace(tmp, out_fullpath)


def wait_for_files(file_groups: Sequence[Sequence[str]] | None, wait_time: float) -> bool:
    """Waits ``wait_time`` seconds, or until all files of one of the groups exist (then
    returns True).

    file_groups : groups, each a list of paths of files or folders; a path whose name has
                  ``*`` is a pattern (as for MATLAB's dir). Empty groups are left out;
                  without groups it is one pause of wait_time.

    It looks for the files it hasn't found yet after each second of the wait, and less often
    when a look takes longer than a tenth of a second, so that looking takes at most a tenth
    of the wait; there is no look at the end (the caller's next round looks then).
    """
    groups = [[g] if isinstance(g, str) else list(g) for g in (file_groups or [])]
    groups = [g for g in groups if g]
    if not groups:
        _pause(wait_time)
        return False

    paths = sorted({os.fspath(p) for g in groups for p in g})
    index = {p: i for i, p in enumerate(paths)}
    members = [[index[os.fspath(p)] for p in g] for g in groups]
    found = [False] * len(paths)
    parts = [fileparts(p) for p in paths]
    names = [n + e for _, n, e in parts]
    is_pattern = ["*" in n for n in names]

    interval = 1.0
    waited = 0.0
    while wait_time - waited > 1e-6:
        t = min(interval, wait_time - waited)
        _pause(t)
        waited += t
        if wait_time - waited <= 1e-6:
            break

        t0 = _clock()
        # batch_file_exist sees files another node just made; so does a listing
        inds = [i for i in range(len(paths)) if not found[i] and not is_pattern[i]]
        if inds:
            for i, ok in zip(inds, batch_file_exist([paths[i] for i in inds]), strict=True):
                found[i] = ok
        folders: dict[str, list[int]] = {}
        for i in range(len(paths)):
            if not found[i] and is_pattern[i]:
                folders.setdefault(parts[i][0], []).append(i)
        for folder, inds in folders.items():
            listed = _list_folder(folder) or []
            for i in inds:
                found[i] = any(_name_matches(x, names[i]) for x in listed) and bool(
                    _dir_matches(folder, names[i], listed)
                )
        if any(all(found[i] for i in g) for g in members):
            return True
        interval = max(1.0, 10 * (_clock() - t0))
    return False


def _name_matches(name: str, pattern: str) -> bool:
    # names that match pattern (* for any characters; with more than one *, MATLAB's
    # regexptranslate('wildcard') also takes ? for one character)
    parts = pattern.split("*")
    if len(parts) == 2:
        return (
            name.startswith(parts[0])
            and name.endswith(parts[1])
            and len(name) >= len(parts[0]) + len(parts[1])
        )
    rx = "".join(".*" if c == "*" else "." if c == "?" else re.escape(c) for c in pattern)
    return re.fullmatch(rx, name, re.DOTALL) is not None


def _dir_matches(folder: str, pattern: str, listed: list[str]) -> list[str]:
    # what MATLAB's dir(<folder>/<pattern>) finds: names matching the * wildcards that exist
    rx = "".join(".*" if c == "*" else re.escape(c) for c in pattern)
    return [
        x
        for x in listed
        if re.fullmatch(rx, x, re.DOTALL) and os.path.exists(os.path.join(folder, x))
    ]
