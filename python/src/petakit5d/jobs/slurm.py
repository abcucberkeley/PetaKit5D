"""The Slurm commands of the job framework, with timeouts.

Each method issues the command line MATLAB's framework issues (same options in the same order;
SlurmParam, SlurmConstraint and the time option are split into words as a shell would), without
a shell, and parses the answer as MATLAB does. One rule differs: a command that fails without a
Slurm answer (a timeout, a busy or unreachable controller, a missing command) gives ``None``,
"unknown", never "job gone", so that nothing is resubmitted on a failed query. An answer that a
job id is invalid (``Invalid job id``) means the job is gone, as in MATLAB.

Job states as MATLAB's check functions give them: RUNNING 1, PENDING 0, GONE -1 (any other
state, or not in the queue). Job keys are strings ``'<job>_<array task>'`` (or ``'<job>'``).
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
import threading
import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime
from typing import Any

__all__ = [
    "DEFAULT_JOB_LOG_DIR",
    "GONE",
    "PENDING",
    "QUERY_BATCH",
    "QUERY_TIMEOUT",
    "RUNNING",
    "SUBMIT_TIMEOUT",
    "CommandResult",
    "Slurm",
    "SlurmCheck",
    "SubmitResult",
    "TestOnlyResult",
    "check_slurm_cluster",
    "job_key",
    "job_log_dir_or_default",
    "job_log_names",
    "parse_key_values",
    "split_options",
]

RUNNING, PENDING, GONE = 1, 0, -1

QUERY_TIMEOUT = 60.0  # s, sinfo, squeue, scontrol show, sacct
SUBMIT_TIMEOUT = 120.0  # s, sbatch and scontrol update (a busy controller answers slowly)
CANCEL_TIMEOUT = 60.0  # s
QUERY_BATCH = 5000  # job keys per squeue -j call (check_batch_slurm_jobs_status)
TERMINATE_WAIT = 5.0  # s a timed-out command has to exit after SIGTERM before it is left alone

# the job log folder of the wrappers when none is given (relative to the current folder)
DEFAULT_JOB_LOG_DIR = "../job_logs"

_INVALID_ID = "Invalid job id"
_SUBMITTED = re.compile(r"Submitted batch job (\d+)\n")
_TEST_ONLY_START = re.compile(r"to start at (\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)")
# a line of squeue's default output (check_batch_slurm_jobs_status)
_SQUEUE_LINE = re.compile(r"^ +(\d+_?\d+) +.* (R|PD|CG|CD|F|PR|S|ST) ")
_PENDING_FIELDS = "JobArrayID:|,Partition:|,QOS:|,NumCPUs:|,tres-alloc:|,Reason:|"
_ACCOUNTING_FIELDS = ("JobIDRaw", "Submit", "Start", "End", "State")


def job_key(job_id: int | str, array_id: int | str | None = None) -> str:
    """'<job>_<array task>', or '<job>' without an array task"""
    if array_id is None or array_id == "":
        return str(int(job_id))
    return f"{int(job_id)}_{int(array_id)}"


def split_options(text: str | None) -> list[str]:
    """the words of SlurmParam/SlurmConstraint/a time option, as a POSIX shell splits them

    (no variable or glob expansion; MATLAB passed these through a shell)
    """
    return shlex.split(text) if text else []


def parse_key_values(line: str) -> dict[str, str]:
    """the Key=Value tokens of one line of scontrol's one-line output ((null) as '')"""
    out: dict[str, str] = {}
    for m in re.finditer(r"(?:^|\s)([A-Za-z][A-Za-z0-9_/:.]*)=(\S*)", line):
        out.setdefault(m.group(1), "" if m.group(2) == "(null)" else m.group(2))
    return out


@dataclass(frozen=True)
class CommandResult:
    """one Slurm command: its words, exit code (None: not run to its end), output and errors"""

    argv: tuple[str, ...]
    returncode: int | None
    stdout: str = ""
    stderr: str = ""
    timed_out: bool = False
    error: str = ""  # why it could not run (e.g. command not found)

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    @property
    def output(self) -> str:
        """stdout and stderr together, as MATLAB's system() returns them"""
        return self.stdout + self.stderr

    @property
    def answered(self) -> bool:
        """whether Slurm answered (the command ran to its end)"""
        return self.returncode is not None and not self.error


@dataclass(frozen=True)
class SubmitResult:
    """sbatch: state 'submitted' (job_id set), 'failed' (Slurm refused or the command failed)
    or 'unknown' (no answer: a job may have been submitted)"""

    state: str
    job_id: int | None
    result: CommandResult

    @property
    def submitted(self) -> bool:
        return self.state == "submitted"


@dataclass(frozen=True)
class TestOnlyResult:
    """sbatch --test-only: the start Slurm expects (None if it gave none) and the answer"""

    start: datetime | None
    result: CommandResult

    @property
    def refused(self) -> bool:
        """Slurm answered and gave no start (bad partition, account, QOS, ...)"""
        return self.result.answered and self.start is None

    __test__ = False  # (not a pytest test class)


@dataclass(frozen=True)
class SlurmCheck:
    """checkSlurmCluster: whether to use Slurm, the job log folder and the -o/-e names"""

    parse_cluster: bool
    job_log_fname: str
    job_log_error_fname: str
    job_log_dir: str | None


def _reap(proc: subprocess.Popen) -> None:
    try:
        proc.wait()
    except Exception:
        pass


class Slurm:
    """the Slurm commands, run from PATH (env: the environment of the commands, default this
    process's) with a timeout each (None: the defaults above)"""

    def __init__(
        self,
        *,
        query_timeout: float | None = None,
        submit_timeout: float | None = None,
        env: Mapping[str, str] | None = None,
    ) -> None:
        self.query_timeout = QUERY_TIMEOUT if query_timeout is None else query_timeout
        self.submit_timeout = SUBMIT_TIMEOUT if submit_timeout is None else submit_timeout
        self.env = dict(env) if env is not None else None

    # running a command

    def run(self, argv: Sequence[str], timeout: float | None = None) -> CommandResult:
        """run one command and wait for it, at most timeout s (default: query_timeout)

        On a timeout the command gets SIGTERM (as coreutils timeout sends) and is left to exit;
        the result is unknown (returncode None).
        """
        argv = tuple(str(a) for a in argv)
        timeout = self.query_timeout if timeout is None else timeout
        try:
            proc = subprocess.Popen(
                argv,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=self.env,
                text=True,
                errors="replace",
            )
        except OSError as err:
            return CommandResult(argv, None, error=f"{type(err).__name__}: {err}")
        try:
            out, err = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.terminate()
            try:
                out, err = proc.communicate(timeout=TERMINATE_WAIT)
            except subprocess.TimeoutExpired:
                out, err = "", ""
                threading.Thread(target=_reap, args=(proc,), daemon=True).start()
            return CommandResult(argv, None, out or "", err or "", timed_out=True)
        return CommandResult(argv, proc.returncode, out, err)

    def start(self, argv: Sequence[str]) -> CommandResult:
        """start a command without waiting for it (MATLAB's 'cmd &'); returncode None"""
        argv = tuple(str(a) for a in argv)
        try:
            proc = subprocess.Popen(
                argv,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                env=self.env,
            )
        except OSError as err:
            return CommandResult(argv, None, error=f"{type(err).__name__}: {err}")
        threading.Thread(target=_reap, args=(proc,), daemon=True).start()
        return CommandResult(argv, None)

    # sinfo

    def available(self) -> bool | None:
        """sinfo: True if it succeeds, False if it fails or is missing, None on a timeout"""
        r = self.run(["sinfo"])
        if r.timed_out:
            return None
        return r.ok

    # sbatch

    def sbatch(self, argv: Sequence[str]) -> SubmitResult:
        """an sbatch command line (argv[0] 'sbatch'), as the cluster loop and the single-job
        wrapper build it; the job id from 'Submitted batch job N'"""
        r = self.run(argv, self.submit_timeout)
        m = _SUBMITTED.search(r.output)
        if m:
            return SubmitResult("submitted", int(m.group(1)), r)
        if r.timed_out:
            return SubmitResult("unknown", None, r)
        return SubmitResult("failed", None, r)

    def test_only(
        self,
        *,
        cpus_per_task: int,
        slurm_param: str = "",
        slurm_constraint: str = "",
        time_str: str = "",
    ) -> TestOnlyResult:
        """sbatch --test-only --cpus-per-task=<n> <SlurmParam> <SlurmConstraint> <time>
        --wrap=true: when Slurm would start such a job (it submits nothing)"""
        argv = [
            "sbatch",
            "--test-only",
            f"--cpus-per-task={int(cpus_per_task)}",
            *split_options(slurm_param),
            *split_options(slurm_constraint),
            *split_options(time_str),
            "--wrap=true",
        ]
        r = self.run(argv)
        start = None
        m = _TEST_ONLY_START.search(r.output)
        if m and r.ok:
            start = datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%S")
        return TestOnlyResult(start, r)

    # squeue

    def batch_status(
        self, job_ids: Sequence[int], array_ids: Sequence[int] | None = None
    ) -> list[int | None]:
        """check_batch_slurm_jobs_status: per job 1 running, 0 pending, -1 other or gone, None
        unknown (its squeue call failed); ids <= 0 (no job) are -1

        One 'squeue -j <keys>' per 5000 distinct jobs. (MATLAB returned the states in the
        order of the sorted job keys when no key repeated; here each job gets its own.)
        """
        n = len(job_ids)
        if array_ids is not None and len(array_ids) != n:
            raise ValueError("job_ids and array_ids differ in length")
        status: list[int | None] = [GONE] * n
        keys: list[str | None] = []
        for i in range(n):
            jid = int(job_ids[i])
            if jid <= 0:
                keys.append(None)
            elif array_ids is None:
                keys.append(job_key(jid))
            else:
                keys.append(job_key(jid, array_ids[i]))
        order = sorted({k for k in keys if k is not None}, key=_key_order)
        found: dict[str, int | None] = {}
        for b in range(0, len(order), QUERY_BATCH):
            batch = order[b : b + QUERY_BATCH]
            states = self._squeue_default(batch)
            for k in batch:
                found[k] = None if states is None else states.get(k, GONE)
        for i, k in enumerate(keys):
            if k is not None:
                status[i] = found[k]
        return status

    def _squeue_default(self, keys: Sequence[str]) -> dict[str, int] | None:
        """squeue -j <keys> (default format): {key: 1 or 0 or -1} for listed jobs, None unknown"""
        r = self.run(["squeue", "-j", ",".join(keys)])
        if _INVALID_ID in r.output:
            return {}
        if not r.ok:
            return None
        lines = r.output.strip().split("\n")[1:]  # (MATLAB drops the first line: the header)
        states: dict[str, int] = {}
        wanted = set(keys)
        for line in lines:
            m = _SQUEUE_LINE.match(line)
            if m is None or m.group(1) not in wanted:
                continue
            # (a job listed more than once takes its last running or pending line)
            if m.group(2) in ("R", "PD"):
                states[m.group(1)] = RUNNING if m.group(2) == "R" else PENDING
        return states

    def job_status(self, job_id: int, array_id: int | None = None) -> int | None:
        """check_slurm_job_status: 1 running, 0 pending or another state, -1 gone, None unknown"""
        key = job_key(job_id, array_id)
        r = self.run(["squeue", "-j", key])
        out = r.output
        if key in out and _INVALID_ID not in out:
            return RUNNING if " R " in out else PENDING
        if _INVALID_ID in out or r.ok:
            return GONE
        return None

    def job_state(self, key: str) -> str | None:
        """squeue -h -j <key> -o %T (result_writers): the job's state, '' if Slurm doesn't
        list it (or says its id is invalid), None unknown"""
        r = self.run(["squeue", "-h", "-j", key, "-o", "%T"])
        if _INVALID_ID in r.output:
            return ""
        if not r.ok:
            return None
        return r.stdout.strip()

    def live_jobs(self, keys: Sequence[str]) -> set[str] | None:
        """squeue -h -r -o %i -j <keys> (task runners, 1000 keys per call): the keys Slurm still
        has in any state, or None unknown"""
        left: set[str] = set()
        for b in range(0, len(keys), 1000):
            batch = list(keys[b : b + 1000])
            r = self.run(["squeue", "-h", "-r", "-o", "%i", "-j", ",".join(batch)])
            if not r.ok and _INVALID_ID not in r.output:
                return None
            listed = {x.strip() for x in r.stdout.splitlines()}
            left.update(k for k in batch if k in listed)
        return left

    def pending_reasons(self, user: str | None = None) -> dict[str, str] | None:
        """squeue -r -h -u <user> -t PENDING -o '%i %r' (task_runner.sh): {key: reason} of the
        user's pending jobs, or None unknown"""
        if user is None:
            user = os.environ.get("USER") or _login_name()
        r = self.run(["squeue", "-r", "-h", "-u", user, "-t", "PENDING", "-o", "%i %r"])
        if not r.ok:
            return None
        out: dict[str, str] = {}
        for line in r.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
        return out

    def pending_jobs(self) -> list[list[str]] | None:
        """squeue --me -t PD -h -O 'JobArrayID:|,Partition:|,QOS:|,NumCPUs:|,tres-alloc:|,
        Reason:|' (SlurmAlternatives' view): the fields of each line (trailing '|' and spaces
        removed, fields not stripped), or None unknown"""
        r = self.run(["squeue", "--me", "-t", "PD", "-h", "-O", _PENDING_FIELDS])
        if not r.ok:
            return None
        rows = []
        for line in r.stdout.strip().splitlines():
            rows.append(re.sub(r"\|\s*$", "", line).split("|"))
        return rows

    # scancel

    def cancel(self, keys: Iterable[str | int], background: bool = False) -> CommandResult:
        """scancel <key> ... (background: started and not waited for, as 'scancel ... &')"""
        argv = ["scancel", *(str(k) for k in keys)]
        if len(argv) == 1:
            return CommandResult(tuple(argv), 0)
        if background:
            return self.start(argv)
        return self.run(argv, CANCEL_TIMEOUT)

    # scontrol

    def show_nodes(self) -> list[str] | None:
        """scontrol show node -o: the lines with NodeName=, or None unknown"""
        r = self.run(["scontrol", "show", "node", "-o"])
        if not r.ok:
            return None
        return [x for x in r.stdout.strip().splitlines() if "NodeName=" in x]

    def show_config(self) -> str | None:
        """scontrol show config, or None unknown"""
        r = self.run(["scontrol", "show", "config"])
        return r.stdout if r.ok else None

    def show_job(self, key: str) -> str | None:
        """scontrol show job <key>, or None (unknown, or no such job)"""
        r = self.run(["scontrol", "show", "job", str(key)])
        return r.stdout if r.ok else None

    def update_job(self, key: str, fields: str | Mapping[str, Any]) -> CommandResult:
        """scontrol update JobId=<key> <fields> (a string as MATLAB builds it, or a mapping)"""
        words = (
            split_options(fields)
            if isinstance(fields, str)
            else [f"{k}={v}" for k, v in fields.items()]
        )
        return self.run(["scontrol", "update", f"JobId={key}", *words], self.submit_timeout)

    # sacct

    def accounting(
        self, job_ids: Sequence[str | int], fields: Sequence[str] = _ACCOUNTING_FIELDS
    ) -> list[dict[str, str]] | None:
        """sacct -n -P -X --format=<fields> -j <ids> (500 ids per call): one dict per line, or
        None if a call failed"""
        rows: list[dict[str, str]] = []
        ids = [str(x) for x in job_ids]
        for b in range(0, len(ids), 500):
            r = self.run(
                [
                    "sacct",
                    "-n",
                    "-P",
                    "-X",
                    f"--format={','.join(fields)}",
                    "-j",
                    ",".join(ids[b : b + 500]),
                ]
            )
            if not r.ok:
                return None
            for line in r.stdout.strip().splitlines():
                cells = line.split("|")
                if len(cells) >= len(fields):
                    rows.append(dict(zip(fields, cells, strict=False)))
        return rows

    def accounting_state(self, job_id: str | int) -> str | None:
        """sacct -j <id> --format=State (chunk_lock_clean): the output, or None unknown"""
        r = self.run(["sacct", "-j", str(job_id), "--format=State"])
        return r.output if r.ok else None


def _key_order(key: str) -> tuple[int, ...]:
    return tuple(int(x) for x in key.split("_"))


def _login_name() -> str:
    try:
        import getpass

        return getpass.getuser()
    except Exception:
        return ""


def job_log_dir_or_default(job_log_dir: str | None) -> str:
    """an empty jobLogDir (as "jobLogDir": "" in a config file) is the default '../job_logs'"""
    return DEFAULT_JOB_LOG_DIR if not job_log_dir else job_log_dir


def job_log_names(job_log_dir: str) -> tuple[str, str]:
    """the -o and -e names of a job log folder: <dir>/job_%A_%a.out and .err"""
    return f"{job_log_dir}/job_%A_%a.out", f"{job_log_dir}/job_%A_%a.err"


def check_slurm_cluster(
    data_path: str, job_log_dir: str | None = None, slurm: Slurm | None = None
) -> SlurmCheck:
    """checkSlurmCluster.m: whether Slurm is there (sinfo), and the job log folder and names

    Without Slurm parse_cluster is False and there are no names. Without job_log_dir only
    Slurm is checked. A job log folder that doesn't exist is replaced by <data_path>/job_logs
    (made, writable for all); the folder is returned as an absolute path. A timed-out sinfo
    counts as Slurm being there (its commands exist; later queries that fail are unknown).
    """
    slurm = Slurm() if slurm is None else slurm
    there = slurm.available()
    if there is False:
        warnings.warn(
            "A slurm-based computing cluster is not exist. Set parseCluster as false.", stacklevel=2
        )
        return SlurmCheck(False, "", "", job_log_dir)
    if not job_log_dir:
        return SlurmCheck(True, "", "", job_log_dir)
    if not os.path.isdir(job_log_dir):
        warnings.warn(
            f"The job log directory does not exist, use {data_path}/job_logs as job log directory",
            stacklevel=2,
        )
        job_log_dir = f"{data_path}/job_logs"
        if not os.path.isdir(job_log_dir):
            os.makedirs(job_log_dir, exist_ok=True)
            try:
                os.chmod(job_log_dir, os.stat(job_log_dir).st_mode | 0o222)
            except OSError:
                pass
    folder = _dir_folder(job_log_dir)
    out, err = job_log_names(folder)
    return SlurmCheck(True, out, err, folder)


def _dir_folder(path: str) -> str:
    """dir(path).folder: the absolute path of a folder, with '/' separators"""
    return os.path.abspath(path).replace(os.sep, "/")
