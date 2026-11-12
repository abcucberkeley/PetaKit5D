"""Publishing a task's result: MATLAB's publish_result and keep_other_result.

A task writes its result under a temporary name in the result's folder (``temp_path``:
``<final>_<uuid>``), finishes everything else it writes (MIPs, flags, info files) and then
moves the result into place with one rename, because the job framework counts the task as
done, and may cancel its job, once the result exists.

Several writers can make the same result (a call started again on the same data while the
first call's jobs run, retries, extra copies of a task). They make the same result, so the
first one published is kept: it is never deleted or replaced, as a reader (the next step of a
pipeline) may already use it. The later writers remove their own temporary result.

The rename never replaces an existing result (``rename_no_replace``):

- folders: ``rename`` onto an existing folder fails on Linux and macOS unless that folder is
  empty, and always on Windows;
- files on Linux and macOS: a hard link to the final name (it fails if the name exists), then
  the temporary name is removed; where hard links are not supported, a check and a rename;
- files on Windows: ``rename`` fails if the name exists.

A failed rename onto an existing result means another writer published it. On network file
systems a request can be answered twice (a lost reply, a retransmission), so a failure is
checked against the file identity of the temporary result: if the final name is it, this
writer published.

Differences from MATLAB: MATLAB's movefile moved a temporary folder into an existing result
(and publish_result removed it from there); here the rename fails and nothing enters the
result. Removing the temporary result and the flag folder is best effort (a folder that
can't be removed, e.g. for a file another process holds open on NFS, is reported and left),
so that a task that published doesn't fail afterwards.
"""

from __future__ import annotations

import errno
import os
import shutil
import stat
import sys
import time
import uuid as _uuid

from petakit5d.jobs import files

__all__ = [
    "keep_other_result",
    "publish_result",
    "remove_path",
    "rename_no_replace",
    "temp_path",
]

_WINDOWS = os.name == "nt"
_SEPS = "/\\" if _WINDOWS else "/"

# errors of a rename onto an existing name
_EXISTS = {errno.EEXIST, errno.ENOTEMPTY}
# errors of os.link where the file system has no hard links
_NO_LINKS = {errno.EPERM, errno.EOPNOTSUPP, errno.ENOTSUP, errno.ENOSYS, errno.EMLINK}
# a rename on Windows that a scanner or indexer holds up: tries and the first pause (s)
RENAME_TRIES = 5
RENAME_PAUSE = 0.1
# removal of a temporary result or a flag folder: tries and the pause between them (s)
REMOVE_TRIES = 3
REMOVE_PAUSE = 0.5


def temp_path(final: str | os.PathLike[str], uuid: str | None = None) -> str:
    """the temporary name of a result: ``<final>_<uuid>`` (a new uuid if none is given),
    in the result's folder, so that publishing it is one rename"""
    final = os.fspath(final).rstrip(_SEPS) or os.fspath(final)
    return f"{final}_{uuid if uuid else _uuid.uuid4()}"


def keep_other_result(
    result_fullpath: str | os.PathLike[str],
    tmp_fullpath: str | os.PathLike[str] = "",
    flag_path: str | os.PathLike[str] = "",
) -> bool:
    """Whether another writer of the same result already put its result in place.

    Before a task moves its temporary result ``tmp_fullpath`` into place as
    ``result_fullpath``: if the result exists (``batch_file_exist``, which sees a result
    another node just made), it was made by another writer, as a task that overwrites removes
    its result when it starts. Then this task keeps that result and removes its own temporary
    result and its flag folder ``flag_path`` (both optional). Returns True if it kept one.
    """
    result_fullpath = os.fspath(result_fullpath)
    tmp_fullpath = os.fspath(tmp_fullpath)
    flag_path = os.fspath(flag_path)
    kept = files.batch_file_exist([result_fullpath])[0]
    if not kept:
        return False
    print(f"The result {result_fullpath} was made by another job, keep it.", flush=True)
    if tmp_fullpath:
        _check_not_result(tmp_fullpath, result_fullpath)
        remove_path(tmp_fullpath)
    if flag_path:
        remove_path(flag_path, folder_only=True)
    return True


def publish_result(
    tmp_fullpath: str | os.PathLike[str],
    result_fullpath: str | os.PathLike[str],
    flag_path: str | os.PathLike[str] = "",
) -> bool:
    """Moves the temporary result into place as the result, unless another writer did.

    tmp_fullpath    : the temporary result (a folder, e.g. a zarr, or a file), in the
                      result's folder (``temp_path``)
    result_fullpath : the result
    flag_path       : the task's flag folder, removed at the end (optional)

    Returns True if this call published the result. If another writer's result is there
    (``keep_other_result``), or comes between that check and the rename (the rename fails),
    it is kept and the temporary result removed. Without a temporary result nothing is
    published (False); the flag folder is removed in every case.
    """
    tmp_fullpath = os.fspath(tmp_fullpath)
    result_fullpath = os.fspath(result_fullpath)
    flag_path = os.fspath(flag_path)
    _check_not_result(tmp_fullpath, result_fullpath)
    if keep_other_result(result_fullpath, tmp_fullpath, flag_path):
        return False
    published = False
    if os.path.lexists(_strip(tmp_fullpath)):
        published = rename_no_replace(tmp_fullpath, result_fullpath)
        if not published:
            print(f"The result {result_fullpath} was made by another job, keep it.", flush=True)
            remove_path(tmp_fullpath)
    if flag_path:
        remove_path(flag_path, folder_only=True)
    return published


def rename_no_replace(src: str | os.PathLike[str], dst: str | os.PathLike[str]) -> bool:
    """Renames src (a file or a folder) to dst unless dst exists: True if src is now dst,
    False if dst existed (another writer's result; src is left as it is).

    Other errors are raised (src missing, a dst folder that doesn't exist, another file
    system). Linux and macOS replace an existing *empty* folder dst, as their rename does;
    a published result is never empty.
    """
    src = _strip(os.fspath(src))
    dst = _strip(os.fspath(dst))
    st = os.lstat(src)
    if _WINDOWS:
        return _rename_windows(src, dst, st)
    if stat.S_ISDIR(st.st_mode):
        try:
            os.rename(src, dst)
            return True
        except OSError as err:
            return _after_failed_rename(src, dst, st, err)
    try:
        os.link(src, dst)
    except FileExistsError as err:
        if not _is_same(dst, st):
            return _after_failed_rename(src, dst, st, err)
    except OSError as err:
        if err.errno not in _NO_LINKS:
            return _after_failed_rename(src, dst, st, err)
        # no hard links here (some network and FUSE file systems): a check, then a rename
        if os.path.lexists(dst):
            return False
        try:
            os.rename(src, dst)
            return True
        except OSError as err2:
            return _after_failed_rename(src, dst, st, err2)
    # (published: the temporary name goes, best effort)
    remove_path(src)
    return True


def _rename_windows(src: str, dst: str, st: os.stat_result) -> bool:
    # MoveFileEx without MOVEFILE_REPLACE_EXISTING: fails if dst exists, for files and
    # folders. A file that a scanner or indexer has open makes it fail for a moment.
    pause = RENAME_PAUSE
    for k in range(RENAME_TRIES):
        try:
            os.rename(src, dst)
            return True
        except FileExistsError as err:
            return _after_failed_rename(src, dst, st, err)
        except PermissionError as err:
            if os.path.lexists(dst) or k == RENAME_TRIES - 1:
                return _after_failed_rename(src, dst, st, err)
            time.sleep(pause)
            pause *= 2
    return False  # not reached


def _after_failed_rename(src: str, dst: str, st: os.stat_result, err: OSError) -> bool:
    # a failed rename or link: dst is src (the answer to a repeated request), another
    # writer's result, or the error is real
    if _is_same(dst, st):
        if not stat.S_ISDIR(st.st_mode) and _is_same(src, st):
            # (a file published by a link: the temporary name goes, best effort)
            remove_path(src)
        return True
    exists_error = isinstance(err, FileExistsError) or err.errno in _EXISTS
    if os.path.lexists(dst) and (exists_error or _WINDOWS or not os.path.lexists(src)):
        return False
    raise err


def _is_same(path: str, st: os.stat_result) -> bool:
    # whether path is the file or folder st was taken of
    try:
        s = os.lstat(path)
    except OSError:
        return False
    if _WINDOWS and not (s.st_ino and st.st_ino):
        return False
    return (s.st_dev, s.st_ino) == (st.st_dev, st.st_ino)


def remove_path(path: str | os.PathLike[str], folder_only: bool = False) -> bool:
    """Removes a file or a folder with everything in it, best effort: a path that is gone (or
    goes while it is removed, by another writer) is fine; an error is tried again
    (``REMOVE_TRIES``) and then reported on stderr. With folder_only, a file is left alone
    (MATLAB removed a flag path only if it was a folder). Returns True if the path is gone."""
    path = _strip(os.fspath(path))
    last: OSError | None = None
    for k in range(REMOVE_TRIES):
        try:
            st = os.lstat(path)
        except FileNotFoundError:
            return True
        except OSError as err:
            last = err
        else:
            if folder_only and not stat.S_ISDIR(st.st_mode):
                return False
            try:
                if stat.S_ISDIR(st.st_mode):
                    _rmtree(path)
                else:
                    os.unlink(path)
                return True
            except FileNotFoundError:
                if not os.path.lexists(path):
                    return True
            except OSError as err:
                last = err
        if k < REMOVE_TRIES - 1:
            time.sleep(REMOVE_PAUSE)
    print(f"Could not remove {path}: {last}", file=sys.stderr, flush=True)
    return False


def _rmtree(path: str) -> None:
    # rmtree that ignores entries another process removed first, and makes read-only
    # entries writable (Windows)
    def on_error(func, p, exc):
        if isinstance(exc, FileNotFoundError):
            return
        if _WINDOWS and isinstance(exc, PermissionError):
            try:
                os.chmod(p, stat.S_IWRITE)
                func(p)
                return
            except OSError:
                pass
        raise exc

    if sys.version_info >= (3, 12):
        shutil.rmtree(path, onexc=on_error)
    else:
        shutil.rmtree(path, onerror=lambda f, p, ei: on_error(f, p, ei[1]))


def _strip(path: str) -> str:
    # trailing separators off (a root stays)
    return path.rstrip(_SEPS) or path


def _check_not_result(tmp_fullpath: str, result_fullpath: str) -> None:
    # the temporary result must be another path than the result, and not hold it
    if not tmp_fullpath:
        return
    tmp = os.path.normcase(os.path.abspath(_strip(tmp_fullpath)))
    res = os.path.normcase(os.path.abspath(_strip(result_fullpath)))
    if tmp == res or res.startswith(tmp.rstrip(os.sep) + os.sep):
        raise ValueError(
            f"The temporary result {tmp_fullpath} is the result {result_fullpath} or holds it"
        )
