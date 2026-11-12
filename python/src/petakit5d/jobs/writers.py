"""Records of the writers of a result: MATLAB's result_writers.

A task that makes a result in parts with jobs of its own (large-file decon, deskew/rotate and
crop) records its process in its flag folder, so that another task of the same result (a call
on the same data started while the first call's jobs still run) waits for the one at work
instead of making the result a second time, and so that the block jobs of a writer that has
ended skip their work.

The record is ``<result folder>/zarr_flag/<name>_<uuid>/writer.txt``, as MATLAB writes it::

    result=<result path>
    uuid=<uuid>
    owner=<job>_<array task>:<node>:<pid>:<process start>
    ended=1                      (once the writer has ended without publishing)

It goes with the flag folder once the result is published. Tasks of the same uuid (one call
and its retries) share the record. MATLAB and Python writers of the same data read each
other's records: the format, the file names and the order of records (their modification
time in whole seconds, then the uuid) are MATLAB's.

A writer "may still work" if its record isn't marked ended, and its owner is a Slurm job
other than this process's own that Slurm lists (pending or running), or (outside Slurm, or in
this process's own job) a process other than this one that runs on this node and started at
the recorded time (Linux; elsewhere such writers are not waited for, as in MATLAB).

On a network file system a listing can miss a record another node made within the same
second as an earlier listing (NFS servers with whole-second times keep the folder's time,
so the client keeps its listing up to its attribute cache time): two writers that register
within a second of each other may then both compute. Their publish stays safe
(``jobs.publish``). Before the wait lists the flag folder, it changes the folder from this
client (``refresh_listing``, Linux), which makes the client list it again.

One change from MATLAB: a squeue that fails (no answer, a busy controller) counts as "the job
may still work" for the wait, as the framework treats a failed query as unknown; MATLAB
counted the job as gone. Without a squeue command at all, Slurm owners count as gone.

The squeue answers of ``orphaned`` are kept for ``max_age`` seconds per job in a session
object (``ResultWriters``); ``reset()`` forgets them (tests).
"""

from __future__ import annotations

import os
import platform
import re
import shutil
import sys
import time
import uuid as _uuid
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass

from petakit5d.jobs import files

__all__ = [
    "FLAG_FOLDER",
    "RECORD_NAME",
    "REFRESH_PREFIX",
    "SQUEUE_MAX_AGE",
    "WAIT_INTERVAL",
    "ResultWriters",
    "live",
    "my_owner",
    "normalize_result_path",
    "orphaned",
    "recorded_writer",
    "refresh_listing",
    "register",
    "release",
    "reset",
    "wait",
]

FLAG_FOLDER = "zarr_flag"
RECORD_NAME = "writer.txt"
SQUEUE_MAX_AGE = 60.0  # s an answer of squeue about a job is kept (orphaned)
WAIT_INTERVAL = 60.0  # s between looks at the other writers while the result doesn't come
# the file refresh_listing makes and removes in a flag folder (no record's name: it starts
# with a dot)
REFRESH_PREFIX = ".pk5d_refresh_"

_WINDOWS = os.name == "nt"
_LINUX = sys.platform.startswith("linux")
_SEPS = "/\\" if _WINDOWS else "/"
_RECORD = re.compile(r"result=(?P<result>[^\n]*)\nuuid=(?P<uuid>[^\n]*)\nowner=(?P<owner>[^\n]*)")
_UUID_OWNER = re.compile(r"uuid=(?P<uuid>[^\n]*)\nowner=(?P<owner>[^\n]*)")
_OWNER_JOB = re.compile(r"owner=(\d+)_(\d*):")
_OWNER_LINE = re.compile(r"owner=([^\n]*)")
_OWNER = re.compile(r"^(\d+)_(\d*):([^:]*):(\d+):?(\d*)$")
_ENDED = re.compile(r"^ended=", re.MULTILINE)
# a record replaced on Windows while a reader has it open: tries and the pause between (s)
_REPLACE_TRIES = 5
_REPLACE_PAUSE = 0.1


@dataclass
class _Record:
    owner: str
    uuid: str
    time: int  # modification time of the record, whole seconds


def normalize_result_path(result_fullpath: str | os.PathLike[str]) -> str:
    """the result path as the records hold it: repeated '/' as one (a leading '//' stays),
    no trailing '/' (the tasks spell paths with repeated slashes)"""
    p = re.sub(r"(?<!^)/+", "/", os.fspath(result_fullpath))
    return p.rstrip("/")


class ResultWriters:
    """result_writers with its state (squeue's answers for ``orphaned``) in this object"""

    def __init__(self) -> None:
        # job key -> (time of the answer, orphaned)
        self._answers: dict[str, tuple[float, bool]] = {}

    def reset(self) -> None:
        """forget squeue's answers"""
        self._answers.clear()

    # the record of this process

    def register(
        self, result_fullpath: str | os.PathLike[str], uuid: str, flag_path: str | os.PathLike[str]
    ) -> None:
        """records this process as the writer of uuid of the result, in flag_path (a folder
        that exists: the task's flag folder)"""
        result = normalize_result_path(result_fullpath)
        _write_record(flag_path, f"result={result}\nuuid={uuid}\nowner={my_owner()}\n")

    def release(
        self, result_fullpath: str | os.PathLike[str], uuid: str, flag_path: str | os.PathLike[str]
    ) -> None:
        """marks the record in flag_path ended if it is this process's (of uuid) and not ended
        yet; nothing if it is gone (the result was published). Never raises: it runs when the
        task ends, also by an error."""
        del result_fullpath
        try:
            fn = _record_file(flag_path)
            txt = _read(fn)
            if txt is None:
                return
            rec = _UUID_OWNER.search(txt)
            if (
                rec is None
                or rec["uuid"] != uuid
                or rec["owner"] != my_owner()
                or _ENDED.search(txt)
            ):
                return
            _write_record(flag_path, f"{txt.strip()}\nended=1\n")
        except Exception:
            pass

    # the other writers

    def live(self, result_fullpath: str | os.PathLike[str], uuid: str) -> list[str]:
        """the owners of the records of other uuids of the result that may still work"""
        return self._live(normalize_result_path(result_fullpath), uuid, earlier_only=False)

    def wait(self, result_fullpath: str | os.PathLike[str], uuid: str) -> bool:
        """After ``register``: waits while a writer recorded before this one may still work on
        the result. True if the result exists then (another writer made it), False if none
        works on it any more (this writer goes on)."""
        result = normalize_result_path(result_fullpath)
        announced = False
        while True:
            if files.batch_file_exist([result])[0]:
                print(f"The result {result} was made by another job.", flush=True)
                return True
            writers = self._live(result, uuid, earlier_only=True)
            if not writers:
                return False
            if not announced:
                print(
                    f"Another job ({writers[0]}) makes the result {result}, wait for it...",
                    flush=True,
                )
                announced = True
            # (back at once when the result comes)
            files.wait_for_files([[result]], WAIT_INTERVAL)

    def orphaned(self, flag_fullpath: str | os.PathLike[str], max_age: float | None = None) -> bool:
        """For a block task of a recorded writer (flag_fullpath: the block's flag, in the
        writer's flag folder): True if the writer has ended (its record says so, or its Slurm
        job surely ended: squeue answered without it), so that nothing would use the block.
        A failed squeue counts as running. squeue's answer about a job is kept for max_age
        seconds (default 60); the record is read each time."""
        max_age = SQUEUE_MAX_AGE if max_age is None else float(max_age)
        fn = _join(files.fileparts(os.fspath(flag_fullpath))[0], RECORD_NAME)
        if not os.path.isfile(fn):
            return False
        txt = _read(fn)
        if txt is None:
            return False
        if _ENDED.search(txt):
            m = _OWNER_LINE.search(txt)
            owner = m[1] if m else ""
            print(f"The task that runs this block ({owner}) has ended, skip the block.", flush=True)
            return True
        rec = _OWNER_JOB.search(txt)
        if rec is None or int(rec[1]) <= 0:
            return False
        my_job, my_arr = own_job()
        if int(rec[1]) == my_job and rec[2] == my_arr:
            return False
        job = f"{rec[1]}_{rec[2]}" if rec[2] else rec[1]
        t = time.monotonic()
        answer = self._answers.get(job)
        if answer is not None and t - answer[0] < max_age:
            tf = answer[1]
        else:
            from petakit5d.jobs.slurm import Slurm

            state = Slurm().job_state(job)
            tf = state == ""
            self._answers[job] = (t, tf)
        if tf:
            print(
                f"The task that runs this block (job {job}) has ended, skip the block.", flush=True
            )
        return tf

    def _live(self, result: str, uuid: str, earlier_only: bool) -> list[str]:
        # owners of the other uuids' records of the result that may still work (with
        # earlier_only, of those recorded before this uuid's record)
        result_path, fsname, _ = files.fileparts(result)
        flag_root = _join(result_path, FLAG_FOLDER)
        if not os.path.isdir(flag_root):
            return []
        if earlier_only:
            # (the wait decides whether this writer computes: list afresh)
            refresh_listing(flag_root)
        try:
            names = sorted(
                e.name
                for e in os.scandir(flag_root)
                if e.name.startswith(f"{fsname}_") and e.is_dir()
            )
        except OSError:
            return []
        records: list[_Record] = []
        mine: int | None = None
        for name in names:
            fn = _join(_join(flag_root, name), RECORD_NAME)
            try:
                mtime = os.stat(fn).st_mtime_ns // 1_000_000_000
            except OSError:
                continue
            txt = _read(fn)
            if txt is None:
                continue
            rec = _RECORD.search(txt)
            if rec is None or rec["result"] != result:
                continue
            if rec["uuid"] == uuid:
                mine = mtime
                continue
            # (a writer that has ended works on it no more)
            if _ENDED.search(txt):
                continue
            records.append(_Record(rec["owner"], rec["uuid"], mtime))
        if earlier_only and mine is not None:
            # (recorded at the same time: the smaller uuid first)
            records = [r for r in records if r.time < mine or (r.time == mine and r.uuid < uuid)]
        live = _owners_live([r.owner for r in records])
        return [r.owner for r, ok in zip(records, live, strict=True) if ok]


def _owners_live(owners: Sequence[str]) -> list[bool]:
    # whether the owners <job>_<array task>:<node>:<pid>[:<start>] may still work: a Slurm
    # job other than this process's own that is pending or running (one squeue for all array
    # jobs), or else a process other than this one that runs on this node
    tf = [False] * len(owners)
    my_job, my_arr = own_job()
    arrays: list[tuple[int, int, int]] = []
    singles: list[tuple[int, int]] = []
    for i, owner in enumerate(owners):
        m = _OWNER.match(owner)
        if m is None:
            continue
        job = int(m[1])
        if job <= 0 or (job == my_job and m[2] == my_arr):
            tf[i] = _process_live(m[3], m[4], m[5])
        elif m[2]:
            arrays.append((i, job, int(m[2])))
        else:
            singles.append((i, job))
    if not arrays and not singles:
        return tf
    from petakit5d.jobs.slurm import Slurm

    slurm = Slurm()
    # (unknown: a failed squeue may still have the job, unless there is no squeue at all)
    unknown_live = shutil.which("squeue") is not None
    if arrays:
        status = slurm.batch_status([a[1] for a in arrays], [a[2] for a in arrays])
        for (i, _, _), s in zip(arrays, status, strict=True):
            tf[i] = unknown_live if s is None else s >= 0
    for i, job in singles:
        s = slurm.job_status(job)
        tf[i] = unknown_live if s is None else s >= 0
    return tf


def _process_live(node: str, pid: str, start: str) -> bool:
    # whether the process pid runs on this node and is not this one: on Linux (/proc), and
    # the one that started at start if it is known (a process id is used again once its
    # process ends)
    return (
        node == node_name()
        and int(pid) != os.getpid()
        and os.path.isdir(f"/proc/{pid}")
        and (not start or start == process_start(pid))
    )


def my_owner() -> str:
    """this process as the owner of a record: <job>_<array task>:<node>:<pid>:<start>"""
    job, arr = own_job()
    pid = os.getpid()
    return f"{job}_{arr}:{node_name()}:{pid}:{process_start(pid)}"


def own_job() -> tuple[int, str]:
    """this process's Slurm job and array task ('' if not an array job; (0, '') without a
    job)"""
    job = _job_number(os.environ.get("SLURM_ARRAY_JOB_ID", ""))
    arr = os.environ.get("SLURM_ARRAY_TASK_ID", "")
    if job is None:
        job = _job_number(os.environ.get("SLURM_JOB_ID", ""))
        arr = ""
    if job is None:
        return 0, ""
    return job, arr


def _job_number(text: str) -> int | None:
    # MATLAB's str2double of a job id (NaN: None)
    try:
        return int(text.strip())
    except ValueError:
        try:
            x = float(text)
        except ValueError:
            return None
        return int(x) if x == x and abs(x) != float("inf") else None


def process_start(pid: int | str) -> str:
    """the start time of a process (clock ticks after boot, the 22nd field of
    /proc/<pid>/stat on Linux; '' elsewhere or if it can't be read)"""
    if not _LINUX:
        return ""
    try:
        with open(f"/proc/{pid}/stat", encoding="utf-8", errors="replace") as f:
            txt = f.read()
    except OSError:
        return ""
    # (the fields after the name in parentheses, which may have spaces)
    fields = txt[txt.rfind(")") + 1 :].split()
    return fields[19] if len(fields) >= 20 else ""


def node_name() -> str:
    """the node name, as the task runners record it"""
    node = os.environ.get("SLURMD_NODENAME", "")
    if not node and _WINDOWS:
        node = os.environ.get("COMPUTERNAME", "")
    elif not node:
        node = platform.node()
    return node


def _record_file(flag_path: str | os.PathLike[str]) -> str:
    return _join(os.fspath(flag_path).rstrip(_SEPS), RECORD_NAME)


def _join(folder: str, name: str) -> str:
    # folder/name with '/' as MATLAB's code spells it
    return f"{folder.rstrip(_SEPS)}/{name}" if folder.rstrip(_SEPS) else f"{folder}{name}"


def refresh_listing(folder: str) -> None:
    """Makes the next listing of folder fresh on Linux. An NFS client keeps a folder's listing
    for its attribute cache time (30 to 60 s), so that it misses the records other nodes just
    made there (seen: a record 6 s old missing on another node); a change of the folder from
    this client (an empty file made and removed) makes the client list it again."""
    if not _LINUX:
        return
    fn = _join(folder, f"{REFRESH_PREFIX}{_uuid.uuid4().hex}")
    try:
        os.close(os.open(fn, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666))
    except OSError:
        return
    try:
        os.unlink(fn)
    except OSError:
        pass


def _read(fn: str) -> str | None:
    try:
        with open(fn, encoding="utf-8", errors="replace", newline="") as f:
            return f.read()
    except OSError:
        return None


def _write_record(flag_path: str | os.PathLike[str], txt: str) -> None:
    # writer.txt in flag_path written at once (a temporary file renamed into place); nothing
    # if the folder can't be written (MATLAB's fopen failing)
    fn = _record_file(flag_path)
    tmp = f"{fn[:-4]}_{_uuid.uuid4()}.txt"
    try:
        f = open(tmp, "w", encoding="utf-8", newline="")
    except OSError:
        return
    with f:
        f.write(txt)
    for k in range(_REPLACE_TRIES):
        try:
            os.replace(tmp, fn)
            return
        except PermissionError:
            if not _WINDOWS or k == _REPLACE_TRIES - 1:
                raise
            time.sleep(_REPLACE_PAUSE)


# the session of this process

_session = ResultWriters()


def reset() -> None:
    """forget squeue's answers kept for orphaned (tests)"""
    _session.reset()


def register(
    result_fullpath: str | os.PathLike[str], uuid: str, flag_path: str | os.PathLike[str]
) -> None:
    """result_writers('register', ...): see ResultWriters.register"""
    _session.register(result_fullpath, uuid, flag_path)


def release(
    result_fullpath: str | os.PathLike[str], uuid: str, flag_path: str | os.PathLike[str]
) -> None:
    """result_writers('release', ...): see ResultWriters.release"""
    _session.release(result_fullpath, uuid, flag_path)


def live(result_fullpath: str | os.PathLike[str], uuid: str) -> list[str]:
    """result_writers('live', ...): see ResultWriters.live"""
    return _session.live(result_fullpath, uuid)


def wait(result_fullpath: str | os.PathLike[str], uuid: str) -> bool:
    """result_writers('wait', ...): see ResultWriters.wait"""
    return _session.wait(result_fullpath, uuid)


def orphaned(flag_fullpath: str | os.PathLike[str], max_age: float | None = None) -> bool:
    """result_writers('orphaned', ...): see ResultWriters.orphaned"""
    return _session.orphaned(flag_fullpath, max_age)


@contextmanager
def recorded_writer(
    result_fullpath: str | os.PathLike[str], uuid: str, flag_path: str | os.PathLike[str]
) -> Iterator[None]:
    """``register`` on entry, ``release`` on exit (MATLAB's onCleanup of the task): a task
    that ends other than by publishing marks its record ended, so that nobody waits for it"""
    register(result_fullpath, uuid, flag_path)
    try:
        yield
    finally:
        release(result_fullpath, uuid, flag_path)
