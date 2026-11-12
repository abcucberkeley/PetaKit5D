"""Claims of the task runners: which process runs which task of a call (MATLAB's claims of
task_runner.sh and task_runner_run).

The claims of a call are in ``<call folder>/claims/``, shared by all its jobs and its master:

    t<task>      a symbolic link whose target is the claimer (its owner); made atomically, also
                 on NFS (it fails if the name exists), so that one process runs the task
    f<task>      a symbolic link to the owner: the task failed in that process
    r<batch>     a file: the master submitted the batch again; other jobs leave its tasks to
                 the batch's own job
    s<job>_<a>   a symbolic link to the owner: the job's runner took tasks of other batches

Tasks and batches are 1-based in the names (MATLAB's) and 0-based (tasks) in the functions
here. An owner is ``<job>_<array task>:<node>:<pid>:<slot>``, with job and array task 0
outside Slurm. A claimed task whose output exists is skipped by its claimer. Only the master
removes the claims of processes that ended (``release``); a runner removes only its own claim
of a task it gives back.

An owner may still run its task (its claim is "live") while its Slurm job is in the queue
(running or pending), or, outside a job or in the caller's own job, while its process runs on
this node and is not the caller. A claim marked failed is not live, and an owner this module
didn't write is taken as live (MATLAB's rules). When squeue gives no answer, liveness is
unknown (None) and nothing may be released on it (MATLAB counted such claims as ended).

This module imports only the standard library and the framework's writers module (the job
entry point imports it before any task).
"""

from __future__ import annotations

import os
import re
import sys
from collections.abc import Callable, Iterable, Sequence

from petakit5d.jobs.writers import node_name, own_job

__all__ = [
    "CLAIMS_FOLDER",
    "Claims",
    "claimed_set",
    "claims_folder",
    "job_key_of",
    "owner",
    "owner_job_keys",
    "owners_live",
    "parse_owner",
    "resubmitted_set",
]

CLAIMS_FOLDER = "claims"

_OWNER = re.compile(r"^(\d+)_(\d+):([^:]*):(\d+):(\d*)")
_LINUX = sys.platform.startswith("linux")


def claims_folder(call_folder: str | os.PathLike[str]) -> str:
    return os.path.join(os.fspath(call_folder), CLAIMS_FOLDER)


def owner(slot: int = 1) -> str:
    """this process as a claimer: ``<job>_<array task>:<node>:<pid>:<slot>`` (job 0 outside
    Slurm, array task 0 when the job is not an array job)"""
    job, arr = own_job()
    return f"{job}_{arr or 0}:{node_name()}:{os.getpid()}:{int(slot)}"


def parse_owner(text: str) -> tuple[int, int, str, int, str] | None:
    """(job, array task, node, pid, slot) of an owner, None if this module didn't write it"""
    m = _OWNER.match(text or "")
    if m is None:
        return None
    return int(m[1]), int(m[2]), m[3], int(m[4]), m[5]


def job_key_of(text: str) -> str | None:
    """'<job>_<array task>' of an owner in a Slurm job, else None"""
    p = parse_owner(text)
    if p is None or p[0] <= 0:
        return None
    return f"{p[0]}_{p[1]}"


def owner_job_keys(owners: Sequence[str]) -> set[str]:
    """the keys '<job>_<array task>' of the owners whose liveness needs Slurm (in a job other
    than this process's own)"""
    my_job, my_arr = own_job()
    my_arr_n = int(my_arr) if my_arr.isdigit() else 0
    out = set()
    for text in owners:
        p = parse_owner(text)
        if p is not None and p[0] > 0 and not (p[0] == my_job and p[1] == my_arr_n):
            out.add(f"{p[0]}_{p[1]}")
    return out


def owners_live(owners: Sequence[str], listed: set[str] | None) -> list[bool | None]:
    """Whether each owner may still run its task: its Slurm job (other than this process's
    own) is in ``listed`` (the job keys Slurm has, running or pending; None: Slurm didn't
    answer, and then such owners are unknown, None), or (outside a job, or in this process's
    own job) its process runs on this node and is not this process."""
    out: list[bool | None] = [False] * len(owners)
    my_job, my_arr = own_job()
    my_arr_n = int(my_arr) if my_arr.isdigit() else 0
    for i, text in enumerate(owners):
        p = parse_owner(text)
        if p is None:
            out[i] = True  # (an owner the runners didn't write: taken as live)
            continue
        job, arr, node, pid, _ = p
        if job > 0 and not (job == my_job and arr == my_arr_n):
            out[i] = None if listed is None else f"{job}_{arr}" in listed
        else:
            out[i] = _process_live(node, pid)
    return out


def _process_live(node: str, pid: int) -> bool:
    return node == node_name() and pid != os.getpid() and _LINUX and os.path.isdir(f"/proc/{pid}")


class Claims:
    """The claims folder of a call (made on first use)."""

    def __init__(self, call_folder: str | os.PathLike[str]) -> None:
        self.folder = claims_folder(call_folder)
        self._made = False

    def _path(self, name: str) -> str:
        return os.path.join(self.folder, name)

    def _make(self) -> None:
        if not self._made:
            os.makedirs(self.folder, exist_ok=True)
            self._made = True

    # -- claims ------------------------------------------------------------------------------

    def usable(self) -> bool:
        """whether claims can be made here: the folder exists or can be made, and takes
        symbolic links (a probe, made and removed)"""
        try:
            self._make()
            probe = self._path(f".probe_{os.getpid()}_{id(self)}")
            os.symlink("probe", probe)
            os.unlink(probe)
        except OSError:
            return False
        return True

    def claim(self, index: int, who: str) -> bool:
        """claims task ``index`` (0-based) for ``who``; False if another process has it.
        Other errors (a file system without symbolic links) raise OSError."""
        self._make()
        try:
            os.symlink(who, self._path(f"t{index + 1}"))
        except FileExistsError:
            return False
        return True

    def claim_many(self, indices: Iterable[int], who: str) -> list[bool]:
        return [self.claim(i, who) for i in indices]

    def give_back(self, index: int, who: str) -> None:
        """removes this process's own claim of a task it doesn't run (only if it is ``who``'s)"""
        p = self._path(f"t{index + 1}")
        try:
            if os.readlink(p) == who:
                os.unlink(p)
        except OSError:
            pass

    def owner_of(self, index: int) -> str | None:
        """the owner of task ``index``'s claim, None without a claim"""
        try:
            return os.readlink(self._path(f"t{index + 1}"))
        except OSError:
            return None

    def mark_failed(self, index: int, who: str) -> None:
        self._make()
        try:
            os.symlink(who, self._path(f"f{index + 1}"))
        except OSError:
            pass

    def mark_stealer(self, job_key: str, who: str) -> None:
        """records that the job's runner took tasks of other batches (for the master)"""
        self._make()
        try:
            os.symlink(who, self._path(f"s{job_key}"))
        except OSError:
            pass

    def mark_resubmitted(self, batch: int) -> None:
        """batch (1-based) was submitted again: other jobs' runners leave its tasks to its job"""
        self._make()
        try:
            with open(self._path(f"r{int(batch)}"), "w"):
                pass
        except OSError:
            pass

    def resubmitted(self, batch: int) -> bool:
        return os.path.lexists(self._path(f"r{int(batch)}"))

    def release(self, indices: Iterable[int]) -> None:
        """removes the claims and failed marks of tasks ``indices`` (the master, for claims
        whose owners ended)"""
        for i in indices:
            for kind in ("t", "f"):
                try:
                    os.unlink(self._path(f"{kind}{i + 1}"))
                except FileNotFoundError:
                    pass

    # -- reading -----------------------------------------------------------------------------

    def names(self) -> set[str]:
        """the names in the claims folder (listing it shows names other nodes just made)"""
        try:
            return set(os.listdir(self.folder))
        except FileNotFoundError:
            return set()

    def owners(self, indices: Sequence[int], names: set[str]) -> dict[int, str]:
        """{task: owner} of the claimed tasks among ``indices`` that are not marked failed"""
        out = {}
        for i in indices:
            if f"t{i + 1}" in names and f"f{i + 1}" not in names:
                who = self.owner_of(i)
                if who is not None:
                    out[i] = who
        return out

    def states(
        self,
        indices: Sequence[int],
        names: set[str],
        listed: Callable[[set[str]], set[str] | None],
    ) -> tuple[list[bool], list[bool | None]]:
        """for tasks ``indices`` (0-based): whether each is claimed, and whether its claimer
        may still run it (False for claims marked failed; None: unknown). ``listed(keys)``
        gives the job keys Slurm has (None: no answer); it is called only if an owner is in
        another Slurm job."""
        claimed = [f"t{i + 1}" in names for i in indices]
        owned = self.owners(indices, names)
        live: list[bool | None] = [False] * len(indices)
        if owned:
            keys = owner_job_keys(list(owned.values()))
            seen = listed(keys) if keys else set()
            states = dict(zip(owned, owners_live(list(owned.values()), seen), strict=True))
            for k, i in enumerate(indices):
                if i in states:
                    live[k] = states[i]
        return claimed, live

    def stealers(self, names: set[str] | None = None) -> dict[str, str]:
        """{job key: owner} of the jobs whose runners took tasks of other batches"""
        names = self.names() if names is None else names
        out = {}
        for n in names:
            if n.startswith("s"):
                try:
                    out[n[1:]] = os.readlink(self._path(n))
                except OSError:
                    continue
        return out


def claimed_set(names: Iterable[str]) -> set[int]:
    """0-based tasks with a claim among folder names"""
    out = set()
    for n in names:
        if n.startswith("t") and n[1:].isdigit():
            out.add(int(n[1:]) - 1)
    return out


def resubmitted_set(names: Iterable[str]) -> set[int]:
    """1-based batches marked resubmitted among folder names"""
    out = set()
    for n in names:
        if n.startswith("r") and n[1:].isdigit():
            out.add(int(n[1:]))
    return out
