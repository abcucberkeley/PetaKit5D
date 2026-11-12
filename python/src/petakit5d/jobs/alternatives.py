"""SlurmAlternatives: the jobs of a call on more than one place of a Slurm cluster.

The port of MATLAB's slurm_alternatives. A cluster config lists other places for its jobs
after its main settings, in order of preference::

    {"SlurmParam": "-p cpu --qos normal -n1 --mem-per-cpu=21418M", "memPerCPU": 20.9,
     "SlurmAlternatives": [{"SlurmParam": "-p gpu --qos normal -n1 --mem-per-cpu=7800M",
                            "memPerCPU": 7.6, "maxCPUNum": 32, "cpuJobs": true}]}

An entry has any of SlurmParam, SlurmConstraint, memPerCPU, maxCPUNum, minCPUNum, gpuHeadroom
and cpuJobs (names matched without regard to case; others are ignored with a warning). Keys it
leaves out are the main settings', except memPerCPU (the entry's --mem-per-cpu, else the
smallest RealMemory / CPUs of its nodes) and cpuJobs (false). The main settings are the first
place; a config with one place behaves exactly as one without the key.

- **Choice at submission.** A job goes to the first place (the main one first) where Slurm can
  start it now, else to the first place a CPU-only job may wait in (one that isn't of GPU nodes
  only), else to the main place. A job is shaped for the main place; at a place with as much
  memory per CPU as its shape (1% less counts) it gets the shape a config of only that place
  gives it (CPUs from memAllocate and the place's memPerCPU, minCPUNum and maxCPUNum, memory
  from its SlurmParam), at one with less it keeps its CPUs and memory.
- **Moves.** A pending job that waits for resources or a QOS or association limit, has waited 2
  minutes in its place (since its submission or last move) and wasn't looked at in the last
  minute moves to the first other place that can start it now, by one ``scontrol update`` of
  partition, QOS, account, CPUs, memory, GRES, features, node lists and job name (places that
  differ in other options don't trade jobs). Only jobs recorded with ``submitted`` (the
  caller's own) are looked at, and between looks no Slurm command runs.
- **Room.** Whether a place can start a job now comes from a view of the cluster, read at most
  every 20 s: the nodes (``scontrol show node -o``: free CPUs, memory and GPUs; not down,
  drained, reserved or PLANNED for a pending job), after the user's pending jobs that start
  first (``squeue --me -t PD -O ...``, each on the first node of its partition where it fits),
  minus the jobs placed since; a QOS or association limit that a pending job of the user waits
  for, or a job of the user that fits no node, closes its partition; then ``sbatch
  --test-only`` must say "now" (60 s), once per view for a place and shape. A place Slurm
  refuses (partition, account, QOS) is not used again in the session, the main place never.
- **GPU nodes are for jobs that need a GPU.** A CPU-only job (no GPUs in its request at that
  place) never goes to, or moves to, an alternative of GPU nodes only unless its entry says
  ``"cpuJobs": true``; it leaves each idle GPU of a node its share (gpuHeadroom: [CPUs, GB] per
  idle GPU, by default the node's CPUs and memory divided by its GPUs; 0 for none); it doesn't
  go to a partition where GPU jobs of the user are pending (held ones too); and one waiting at
  GPU nodes without room there, or where GPU jobs of the user are pending, moves to the first
  place that isn't of GPU nodes only.

The state (the view, the recorded jobs, the places Slurm refused and the timing settings) is a
``Session``; the module has one for all calls of the process (MATLAB's persistent variables),
``reset()`` for tests. The cluster loop gets its places from ``cluster_places`` (the protocol
of ``jobs.cluster._OnePlace``); single jobs from ``single_job_places`` and their status checks
call ``move`` and ``forget``.

What differs from MATLAB, on purpose:

- A failed or timed-out ``scontrol show node`` or ``squeue`` makes the view unknown until the
  next one (20 s): no place has room (jobs go to the main place, as without alternatives) and
  no job moves. MATLAB read it as no nodes, or as no pending jobs, which could send CPU-only
  jobs to a partition where GPU jobs of the user wait.
- An ``sbatch --test-only`` without an answer (a timeout, a missing command) is an error of
  the moment: the place is not left out.
- Commands run without a shell, so options with ``|``, ``&``, ``*`` or brackets reach Slurm as
  written (MATLAB's shell made ``-C a|b`` a pipe); an alternative whose options a shell can't
  split (an unbalanced quote) is left out with a warning.
- Entries with values of the wrong type (SlurmConstraint, memPerCPU, maxCPUNum, minCPUNum or
  gpuHeadroom that aren't text or numbers) raise an error (MATLAB used them as they came).
- An ``scontrol update`` without an answer leaves the job for 2 minutes (it may have moved);
  the next look finds where it is.
"""

from __future__ import annotations

import math
import numbers
import re
import struct
import sys
import threading
import time
import warnings
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from datetime import datetime
from typing import Any

from petakit5d.jobs.shapes import (
    TIME_FORMAT_BATCH,
    TIME_FORMAT_MCC,
    TIME_FORMAT_SINGLE,
    job_cpus,
    time_str,
)
from petakit5d.jobs.slurm import Slurm, parse_key_values

__all__ = [
    "ENTRY_KEYS",
    "LOOK_INTERVAL",
    "MIN_WAIT",
    "START_TOLERANCE",
    "VIEW_AGE",
    "AlternativesError",
    "AlternativesWarning",
    "ClusterPlaces",
    "Options",
    "Place",
    "Session",
    "Shape",
    "SingleJobPlaces",
    "cluster_places",
    "default_session",
    "expand_hostlist",
    "forget",
    "memory_mb",
    "move",
    "move_fields",
    "parse_options",
    "reset",
    "single_job_places",
    "time_strs",
    "with_memory",
]

# --- timing (s), MATLAB's; Session.set_timing changes them --------------------------------------

VIEW_AGE = 20.0  # a view of the cluster older than this is read again
MIN_WAIT = 120.0  # a job moves only after this long in its place (since submission or last move)
LOOK_INTERVAL = 60.0  # a waiting job is looked at at most this often
START_TOLERANCE = 60.0  # an sbatch --test-only start within this many s is "now"

# --- rules ------------------------------------------------------------------------------------

# the keys of an entry
ENTRY_KEYS = (
    "SlurmParam",
    "SlurmConstraint",
    "memPerCPU",
    "maxCPUNum",
    "minCPUNum",
    "gpuHeadroom",
    "cpuJobs",
)
# a place with this fraction of a shape's memory per CPU has as much (memPerCPU values are
# rounded)
SHAPE_RATIO_TOLERANCE = 0.99
# tasks of a pending array that count (more never start at once)
PENDING_TASKS_MAX = 1000
# node states without room (PLANNED: the backfill scheduler keeps the node for a pending job)
BAD_STATES = (
    "DOWN",
    "DRAIN",
    "FAIL",
    "MAINT",
    "RESERVED",
    "NOT_RESPONDING",
    "POWER",
    "REBOOT",
    "FUTURE",
    "INVALID",
    "PLANNED",
    "UNKNOWN",
    "NO_RESPOND",
    "BLOCKED",
)
# reasons of pending jobs that start first (they wait for resources)
STARTING_REASONS = ("Resources", "Priority", "None")
# reasons of pending jobs that may move (and QOS... / Assoc... limits)
MOVABLE_REASONS = ("Resources", "Priority")
# words of an sbatch --test-only error of the moment (a busy controller): the place stays
TRANSIENT_WORDS = ("timed out", "unable to contact", "try again", "temporarily", "busy")

# MATLAB's NaN (sign bit set: the bits of petakit5d.core.dtypes.NAN, without importing NumPy)
_NAN: float = struct.unpack("<d", struct.pack("<Q", 0xFFF8000000000000))[0]
_LIMIT_REASON = re.compile(r"(QOS|Assoc)")
# MATLAB's whitespace (strtrim, strsplit)
_WS = " \f\n\r\t\v"
_WS_RUN = re.compile(r"[ \f\n\r\t\v]+")
_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?")
_MEMORY = re.compile(r"([0-9.]+)([KMGT]?)B?", re.IGNORECASE)
_MEMORY_UNITS = {"K": 1 / 1024, "M": 1.0, "G": 1024.0, "T": 1024.0**2}
_TRAILING_DIGITS = re.compile(r"([0-9]+)\Z")
_PENDING_ID = re.compile(r"([0-9]+)_\[?([^\]]*)\]?")
_HOST_PART = re.compile(r"[^,\[]+(\[[^\]]*\])?[^,]*")
_HOST_RANGE = re.compile(r"([^\[]*)\[([^\]]*)\](.*)", re.S)
_TRES_MEM = re.compile(r"(?:^|,)mem=([^,]+)")

# the options a SlurmParam may say where a job runs and what it takes with (MATLAB's names:
# the short ones with a separate or an attached value, the long ones with = or a separate one)
_SHORT = {
    "-p": "partition",
    "-q": "qos",
    "-A": "account",
    "-C": "constraint",
    "-x": "exclude",
    "-w": "nodelist",
    "-J": "job_name",
    "-n": "ntasks",
    "-t": "time",
    "-G": "gpus",
}
_LONG = {
    "--partition": "partition",
    "--qos": "qos",
    "--account": "account",
    "--constraint": "constraint",
    "--exclude": "exclude",
    "--nodelist": "nodelist",
    "--job-name": "job_name",
    "--ntasks": "ntasks",
    "--time": "time",
    "--gpus": "gpus",
    "--mem": "mem",
    "--mem-per-cpu": "mem_per_cpu",
    "--gres": "gres",
    "--gpus-per-node": "gpus",
}


class AlternativesError(ValueError):
    """A SlurmAlternatives entry MATLAB refuses; ``identifier`` is MATLAB's error id."""

    def __init__(self, message: str, identifier: str = "slurm_alternatives:entry") -> None:
        super().__init__(message)
        self.identifier = identifier


class AlternativesWarning(UserWarning):
    """A warning of MATLAB's slurm_alternatives; ``identifier`` is MATLAB's warning id."""

    def __init__(self, message: str, identifier: str) -> None:
        super().__init__(message)
        self.identifier = identifier


def _warn(message: str, identifier: str) -> None:
    warnings.warn(AlternativesWarning(message, identifier), stacklevel=3)


# --- MATLAB's arithmetic and text functions on scalars -------------------------------------------


def _isnan(x: float) -> bool:
    return isinstance(x, float) and math.isnan(x)


def _mmin(a: float, b: float) -> float:
    """min of two scalars, NaN ignored (MATLAB)"""
    if _isnan(a):
        return b
    if _isnan(b):
        return a
    return b if b < a else a


def _mmax(a: float, b: float) -> float:
    if _isnan(a):
        return b
    if _isnan(b):
        return a
    return b if b > a else a


def _mdiv(a: float, b: float) -> float:
    """a / b of doubles (x / 0 is Inf, -Inf or NaN)"""
    if b == 0:
        if a == 0 or _isnan(a):
            return _NAN
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


def _mfloor(x: float) -> float:
    return x if math.isinf(x) or _isnan(x) else float(math.floor(x))


def _mceil(x: float) -> float:
    return x if math.isinf(x) or _isnan(x) else float(math.ceil(x))


def _round_half_away(x: float) -> int:
    return int(math.copysign(math.floor(abs(x) + 0.5), x))


def _sprintf_d(x: float) -> str:
    """sprintf('%d') of a double (integral values as integers, others as %e, Inf/NaN by name)"""
    if _isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if float(x).is_integer():
        return str(int(x))
    return f"{x:e}"


def _int(x: float) -> int | float:
    """an integral double as int (CPU counts), else as it is"""
    return int(x) if isinstance(x, float) and x.is_integer() else x


def _strtrim(text: str) -> str:
    return text.strip(_WS)


def _strsplit(text: str, delimiter: str | None = None) -> list[str]:
    """MATLAB's strsplit (repeated delimiters count as one; whitespace without delimiter)"""
    if delimiter is None:
        return _WS_RUN.split(text)
    return re.split(f"(?:{re.escape(delimiter)})+", text)


def _str2double(text: str) -> float:
    """MATLAB's str2double for the numbers Slurm prints (NaN if text is none)"""
    s = _strtrim(text)
    if _NUMBER.fullmatch(s):
        return float(s.replace("d", "e").replace("D", "e"))
    low = s.lower()
    if low in ("inf", "+inf"):
        return math.inf
    if low == "-inf":
        return -math.inf
    return _NAN


def memory_mb(value: str) -> float | None:
    """a Slurm memory size in MB (K, M, G, T; M by default), None if it is none"""
    m = _MEMORY.fullmatch(value)
    if m is None:
        return None
    return _str2double(m.group(1)) * _MEMORY_UNITS[m.group(2).upper() or "M"]


def _gpu_count(value: str) -> float:
    """the count of a --gpus value ([type:]count), 1 without one"""
    m = _TRAILING_DIGITS.search(value)
    return _str2double(m.group(1)) if m else 1.0


def _tres_count(tres: str, name: str) -> float:
    """the count of name in a TRES list (cpu=4,mem=10G,gres/gpu=1); 0 for a GPU or node count
    that isn't there, NaN for others"""
    m = re.search(f"(?:^|,){re.escape(name)}=([^,]+)", tres)
    if m is None:
        return 0.0 if name in ("gres/gpu", "node") else _NAN
    return _str2double(m.group(1))


def _colon(a: float, b: float) -> list[float]:
    """MATLAB's a : b (NaN : x and x : NaN are NaN; a range without end gives nothing)"""
    if _isnan(a) or _isnan(b):
        return [_NAN]
    if math.isinf(a) or math.isinf(b) or b < a:
        return []
    return [a + i for i in range(int(math.floor(b - a)) + 1)]


def _sprintf_0d(x: float, width: int) -> str:
    """sprintf('%0*d', width, x) of a double: integral values zero-padded, NaN and Inf by name
    (padded with spaces), others as %e"""
    if _isnan(x) or math.isinf(x):
        return _sprintf_d(x).rjust(width)
    if float(x).is_integer():
        return f"{int(x):0{width}d}"
    return f"{x:0{width}e}"


def _is_empty(x: Any) -> bool:
    if x is None:
        return True
    if isinstance(x, (str, list, tuple, dict)):
        return len(x) == 0
    size = getattr(x, "size", None)
    return isinstance(size, int) and size == 0


def _unwrap(x: Any) -> Any:
    # a JSON list of one value is the value (MATLAB's jsondecode gives a scalar)
    while isinstance(x, (list, tuple)) and len(x) == 1:
        x = x[0]
    return x


def _is_number(x: Any) -> bool:
    return isinstance(x, numbers.Real) or type(x).__name__ == "bool_"


# --- the options of a SlurmParam ---------------------------------------------------------------


@dataclass(frozen=True)
class Options:
    """what a SlurmParam (and SlurmConstraint) says about where a job runs and what it takes:
    partition, QOS, account, memory (MB per CPU or per node, and the option as given), GPUs,
    constraint, node lists, job name, tasks and time limit; other tokens are kept in order
    (places that differ in them don't trade jobs)"""

    partition: str = ""
    qos: str = ""
    account: str = ""
    mem_per_cpu: float | None = None
    mem: float | None = None
    mem_token: str = ""
    gres: str = ""
    gpus: float = 0.0
    constraint: str = ""
    exclude: str = ""
    nodelist: str = ""
    job_name: str = ""
    ntasks: float = 1.0
    time: str = ""
    other: tuple[str, ...] = ()


def parse_options(text: str) -> Options:
    """parse_options of slurm_alternatives.m: the options of SlurmParam + ' ' + SlurmConstraint
    (tokens split at whitespace as MATLAB's strsplit splits them; the last of an option wins)"""
    tokens = [t for t in _strsplit(_strtrim(text)) if t]
    values: dict[str, Any] = {}
    other: list[str] = []
    i, n = 0, len(tokens)
    while i < n:
        t = tokens[i]
        name, value = "", ""
        if t.startswith("--"):
            eq = t.find("=")
            name = _LONG.get(t if eq < 0 else t[:eq], "")
            if name:
                if eq >= 0:
                    value = t[eq + 1 :]
                elif i < n - 1:
                    i += 1
                    value = tokens[i]
        elif len(t) >= 2 and t[0] == "-":
            name = _SHORT.get(t[:2], "")
            if name:
                if len(t) > 2:
                    value = t[2:]
                elif i < n - 1:
                    i += 1
                    value = tokens[i]
        if not name:
            other.append(t)
        elif name == "mem_per_cpu":
            values["mem_per_cpu"] = memory_mb(value)
            values["mem_token"] = f"--mem-per-cpu={value}"
        elif name == "mem":
            values["mem"] = memory_mb(value)
            values["mem_token"] = f"--mem={value}"
        elif name == "ntasks":
            values["ntasks"] = _str2double(value)
        elif name == "gpus":
            values["gpus"] = _gpu_count(value)
        elif name == "gres":
            values["gres"] = value
            # gpu, gpu:2, gpu:a100:1
            for g in _strsplit(value, ","):
                parts = _strsplit(g, ":")
                if parts[0] == "gpu":
                    gpus = _str2double(parts[-1])
                    values["gpus"] = 1.0 if len(parts) == 1 or _isnan(gpus) else gpus
        else:
            values[name] = value
        i += 1
    return Options(**values, other=tuple(other))


def with_memory(param: str, mem_token: str) -> str:
    """a SlurmParam with its memory option (--mem-per-cpu or --mem) replaced by mem_token (added
    if it has none); tokens split and joined at whitespace as MATLAB does"""
    tokens = _strsplit(_strtrim(param))
    keep = [True] * len(tokens)
    at = -1
    for i, t in enumerate(tokens):
        if not keep[i]:
            continue
        if t in ("--mem", "--mem-per-cpu"):
            keep[i] = False
            if i + 1 < len(tokens):
                keep[i + 1] = False
        elif t.startswith(("--mem=", "--mem-per-cpu=")):
            keep[i] = False
        else:
            continue
        if at < 0:
            at = i
    if at < 0:
        tokens = [*tokens, mem_token]
    else:
        tokens = [
            *tokens[:at],
            mem_token,
            *(t for t, k in zip(tokens[at:], keep[at:], strict=True) if k),
        ]
    return " ".join(tokens)


_SPLITTABLE: dict[str, bool] = {}


def _shell_splittable(text: str) -> bool:
    """whether a shell can split text into words (no unbalanced quote)"""
    ok = _SPLITTABLE.get(text)
    if ok is None:
        import shlex

        try:
            shlex.split(text)
            ok = True
        except ValueError:
            ok = False
        if len(_SPLITTABLE) < 10000:
            _SPLITTABLE[text] = ok
    return ok


def expand_hostlist(text: str) -> list[str]:
    """the names of a Slurm host list (a,b and nd[01-03,07].x)"""
    hosts: list[str] = []
    for m in _HOST_PART.finditer(text):
        part = m.group(0)
        r = _HOST_RANGE.fullmatch(part)
        if r is None:
            hosts.append(part)
            continue
        prefix, ranges, suffix = r.groups()
        for span in _strsplit(ranges, ","):
            ends = _strsplit(span, "-")
            width = len(ends[0])
            for x in _colon(_str2double(ends[0]), _str2double(ends[-1])):
                hosts.append(f"{prefix}{_sprintf_0d(x, width)}{suffix}")
    return hosts


def _short_name(name: str) -> str:
    # a node name without its domain
    dot = name.find(".")
    return name if dot < 0 else name[:dot]


# --- places and shapes ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Place:
    """one place for jobs (the main settings, or an entry with the main settings for the keys
    it leaves out); cpus_per_task is the CPUs a config of only this place gives a job (set by
    the cluster framework)"""

    slurm_param: str
    slurm_constraint: str
    mem_per_cpu: float | None
    max_cpu_num: float
    min_cpu_num: float
    gpu_headroom: tuple[float, ...] | None
    cpu_jobs: bool
    opts: Options
    cpus_per_task: int | float | None = None

    @property
    def name(self) -> str:
        """its partition, or its SlurmParam"""
        return self.opts.partition or self.slurm_param


@dataclass(frozen=True)
class Shape:
    """a job's CPUs and memory at a place: its memory option, memory per CPU (GB) and the
    SlurmParam that asks for it"""

    cpus: int | float
    mem_token: str
    ratio: float | None
    param: str


def _main_shape(places: Sequence[Place], own_cpus: Sequence[float]) -> Shape:
    """a job at the main place: its own CPUs and memory there (the main memPerCPU per CPU if
    its SlurmParam has no memory option)"""
    p = places[0]
    token = p.opts.mem_token
    if not token:
        ratio = _NAN if p.mem_per_cpu is None else p.mem_per_cpu
        token = f"--mem-per-cpu={_sprintf_d(_mfloor(ratio * 1024))}M"
    return Shape(own_cpus[0], token, p.mem_per_cpu, p.slurm_param)


def _ratio_at_least(a: float | None, b: float | None) -> bool:
    # a >= b * 0.99 of doubles ([] and NaN: false)
    if a is None or b is None:
        return False
    return a >= b * SHAPE_RATIO_TOLERANCE


def _place_shape(
    places: Sequence[Place], own_cpus: Sequence[float], k: int, current: Shape
) -> Shape:
    """a job's shape at place k, for a job of shape current: the place's own shape if it has as
    much memory per CPU, else the current CPUs and memory (more CPUs can mean more workers and
    more memory)"""
    p = places[k]
    if _ratio_at_least(p.mem_per_cpu, current.ratio):
        return Shape(own_cpus[k], p.opts.mem_token, p.mem_per_cpu, p.slurm_param)
    return Shape(
        current.cpus,
        current.mem_token,
        current.ratio,
        with_memory(p.slurm_param, current.mem_token),
    )


def _shape_memory(place: Place, shape: Shape) -> float:
    """the memory (MB) of a job of a shape at a place"""
    o = parse_options("x " + shape.mem_token)
    if o.mem_per_cpu is not None:
        return o.mem_per_cpu * shape.cpus * place.opts.ntasks
    if o.mem is not None:
        return o.mem
    ratio = _NAN if place.mem_per_cpu is None else place.mem_per_cpu
    return ratio * 1024 * shape.cpus * place.opts.ntasks


def move_fields(a_place: Place, b_place: Place, current: Shape, shape: Shape) -> dict[str, str]:
    """the scontrol update fields that move a job of shape current from place a to place b with
    a shape, in order ({} if the places differ in options a pending job can't change)

    The QOS is always set when the partition, account or QOS changes (Slurm keeps a QOS the new
    partition's association doesn't allow, and can't clear one), so a place without a QOS
    can't take jobs of another partition, account or QOS; nor can a place without a memory
    option take jobs of another partition (its default memory can't be set back)."""
    a, b = a_place.opts, b_place.opts
    if a.other != b.other or a.ntasks != b.ntasks or a.time != b.time or not b.partition:
        return {}
    f = {"Partition": b.partition}
    if a.partition != b.partition or a.account != b.account or a.qos != b.qos:
        if not b.qos:
            return {}
        f["QOS"] = b.qos
    if a.account != b.account:
        if not b.account:
            return {}
        f["Account"] = b.account
    f["NumCPUs"] = _sprintf_d(shape.cpus * b.ntasks)
    f["CPUsPerTask"] = _sprintf_d(shape.cpus)
    m = parse_options("x " + shape.mem_token)
    if m.mem_per_cpu is not None:
        f["MinMemoryCPU"] = _sprintf_d(_mceil(m.mem_per_cpu))
    elif m.mem is not None:
        f["MinMemoryNode"] = _sprintf_d(_mceil(m.mem))
    elif current.mem_token or a.partition != b.partition:
        return {}
    if a.gres != b.gres:
        if not b.gres or not a.gres:
            return {}
        f["Gres"] = b.gres
    elif a.gpus != b.gpus:
        return {}
    for attr, key in (
        ("constraint", "Features"),
        ("exclude", "ExcNodeList"),
        ("nodelist", "ReqNodeList"),
        ("job_name", "JobName"),
    ):
        va, vb = getattr(a, attr), getattr(b, attr)
        if va != vb:
            f[key] = "wrap" if attr == "job_name" and not vb else vb
    return f


def time_strs(
    places: Sequence[Place], job_time_limit: float | None, fmt: str, main_str: str
) -> list[str]:
    """the time limit option of each place as the framework writes it (fmt: one of
    shapes.TIME_FORMAT_*; none for a place whose SlurmParam has a time limit; the main place's
    is main_str)"""
    return [main_str] + [time_str(job_time_limit, p.slurm_param, fmt) for p in places[1:]]


# --- the view of the cluster ---------------------------------------------------------------------


@dataclass
class _Pending:
    key: str
    partition: str
    qos: str
    cpus: float
    mem: float
    gpus: float
    reason: str
    node: int = -1  # the node it takes in the view (-1: none)


@dataclass
class _View:
    """the cluster's nodes and the user's pending jobs at a time (known: Slurm answered both)"""

    time: float
    known: bool = True
    names: list[str] = field(default_factory=list)
    parts: list[list[str]] = field(default_factory=list)
    feats: list[list[str]] = field(default_factory=list)
    cpu_tot: list[float] = field(default_factory=list)
    mem_tot: list[float] = field(default_factory=list)
    gpu_tot: list[float] = field(default_factory=list)
    free_cpu: list[float] = field(default_factory=list)
    free_mem: list[float] = field(default_factory=list)
    free_gpu: list[float] = field(default_factory=list)
    usable: list[bool] = field(default_factory=list)
    pend: list[_Pending] = field(default_factory=list)
    pend_index: dict[str, int] = field(default_factory=dict)
    limit_parts: list[str] = field(default_factory=list)
    limit_qos: list[str] = field(default_factory=list)
    wait_parts: set[str] = field(default_factory=set)
    gpu_parts: set[str] = field(default_factory=set)
    test_only: dict[str, bool] = field(default_factory=dict)
    masks: dict[str, list[bool]] = field(default_factory=dict)


def _pending_keys(job_array_id: str) -> tuple[list[str], int]:
    """the job keys (job_task) of squeue's JobArrayID of pending jobs: 12_[3], 12_3, 12 (no
    array: never one of ours), 12_[1-4,7%2]; and the number of tasks"""
    m = _PENDING_ID.fullmatch(job_array_id)
    if m is None:
        return [], 1
    spec = re.sub(r"%[0-9]+$", "", m.group(2))
    tasks: list[float] = []
    for part in _strsplit(spec, ","):
        ends = [_str2double(x) for x in _strsplit(part, "-")]
        if len(ends) == 2:
            tasks.extend(_colon(ends[0], ends[1]))
        elif len(ends) == 1 and not _isnan(ends[0]):
            tasks.extend(_colon(ends[0], ends[0]))
    tasks = tasks[:PENDING_TASKS_MAX]
    return [f"{m.group(1)}_{_sprintf_d(x)}" for x in tasks], len(tasks)


def _read_nodes(v: _View, lines: Iterable[str]) -> None:
    for line in lines:
        kv = parse_key_values(line)
        state = kv.get("State", "").upper()
        cpu_tot = _str2double(kv.get("CPUEfctv", ""))
        if _isnan(cpu_tot):
            cpu_tot = _str2double(kv.get("CPUTot", ""))
        mem_spec = _str2double(kv.get("MemSpecLimit", ""))
        mem_tot = _str2double(kv.get("RealMemory", "")) - (0.0 if _isnan(mem_spec) else mem_spec)
        gpu_tot = _tres_count(kv.get("CfgTRES", ""), "gres/gpu")
        free_cpu = cpu_tot - _str2double(kv.get("CPUAlloc", ""))
        free_mem = mem_tot - _str2double(kv.get("AllocMem", ""))
        v.names.append(kv.get("NodeName", ""))
        v.parts.append(_strsplit(kv.get("Partitions", ""), ","))
        v.feats.append(_strsplit(kv.get("AvailableFeatures", ""), ","))
        v.cpu_tot.append(cpu_tot)
        v.mem_tot.append(mem_tot)
        v.gpu_tot.append(gpu_tot)
        v.free_cpu.append(free_cpu)
        v.free_mem.append(free_mem)
        v.free_gpu.append(gpu_tot - _tres_count(kv.get("AllocTRES", ""), "gres/gpu"))
        v.usable.append(
            not any(b in state for b in BAD_STATES)
            and not _isnan(free_cpu)
            and not _isnan(free_mem)
        )


def _read_pending(v: _View, rows: Iterable[Sequence[str]]) -> None:
    for fields in rows:
        # (MATLAB's strsplit takes repeated delimiters as one)
        f = _strsplit("|".join(fields), "|")
        if len(f) < 6:
            continue
        keys, n_tasks = _pending_keys(_strtrim(f[0]))
        tres = _strtrim(f[4])
        cpus = _tres_count(tres, "cpu")
        if _isnan(cpus):
            cpus = _str2double(f[3])
        m = _TRES_MEM.search(tres)
        mem = memory_mb(m.group(1)) if m else None
        if mem is None:
            mem = 0.0
        nodes = _mmax(1.0, _tres_count(tres, "node"))
        keys += [""] * (max(1, n_tasks) - len(keys))
        row = (
            _strtrim(f[1]),
            _strtrim(f[2]),
            cpus / nodes,
            mem / nodes,
            _tres_count(tres, "gres/gpu") / nodes,
            _strtrim("|".join(f[5:])),
        )
        v.pend.extend(_Pending(key, *row) for key in keys)
    for i, p in enumerate(v.pend):
        if p.key:
            v.pend_index[p.key] = i


def _start_pending(v: _View) -> None:
    """the user's pending jobs start first: each takes its demand from the first usable node
    of its partition where it fits (one that fits nowhere takes nothing)"""
    part_masks: dict[str, list[bool]] = {}
    # (free resources only shrink here: a demand that fit no node fits none later)
    fits_nowhere: set[tuple[str, float, float, float]] = set()
    n = len(v.names)
    for p in v.pend:
        if p.reason not in STARTING_REASONS:
            continue
        demand = (p.partition, p.cpus, p.mem, p.gpus)
        if demand in fits_nowhere:
            continue
        mask = part_masks.get(p.partition)
        if mask is None:
            mask = [v.usable[i] and p.partition in v.parts[i] for i in range(n)]
            part_masks[p.partition] = mask
        for j in range(n):
            if (
                mask[j]
                and v.free_cpu[j] >= p.cpus
                and v.free_mem[j] >= p.mem
                and v.free_gpu[j] >= p.gpus
            ):
                _take(v, j, p.cpus, p.mem, p.gpus)
                p.node = j
                break
        else:
            fits_nowhere.add(demand)
    limited = [p for p in v.pend if _LIMIT_REASON.match(p.reason)]
    v.limit_parts = [p.partition for p in limited]
    v.limit_qos = [p.qos for p in limited]
    # partitions where jobs of the user wait for resources (they fit no node): a new job there
    # waits behind them (Slurm starts jobs by priority, and its backfill tests only so many)
    v.wait_parts = {p.partition for p in v.pend if p.reason in STARTING_REASONS and p.node < 0}
    # partitions where GPU jobs of the user are pending, for any reason: CPU-only jobs keep out
    v.gpu_parts = {p.partition for p in v.pend if p.gpus > 0}


def _take(v: _View, j: int, cpus: float, mem: float, gpus: float) -> None:
    v.free_cpu[j] = v.free_cpu[j] - cpus
    v.free_mem[j] = v.free_mem[j] - mem
    v.free_gpu[j] = v.free_gpu[j] - gpus


def _place_nodes(v: _View, place: Place) -> list[bool]:
    """the nodes of a place: in its partition, with its features, within its node lists (a
    place without a partition: none, its room comes from Slurm's check alone)"""
    o = place.opts
    key = f"{o.partition}|{o.constraint}|{o.exclude}|{o.nodelist}"
    mask = v.masks.get(key)
    if mask is not None:
        return mask
    n = len(v.names)
    mask = [o.partition in parts for parts in v.parts] if o.partition else [False] * n
    c = o.constraint
    if c:
        # a&b, a|b and single features (others are left to Slurm's check)
        if "&" in c and "|" not in c:
            need = _strsplit(c, "&")
            mask = [m and all(x in f for x in need) for m, f in zip(mask, v.feats, strict=True)]
        elif not any(ch in "[]()*&" for ch in c):
            either = _strsplit(c, "|")
            mask = [m and any(x in f for x in either) for m, f in zip(mask, v.feats, strict=True)]
    short = [_short_name(x) for x in v.names]
    if o.exclude:
        out = {_short_name(x) for x in expand_hostlist(o.exclude)}
        mask = [m and s not in out for m, s in zip(mask, short, strict=True)]
    if o.nodelist:
        inside = {_short_name(x) for x in expand_hostlist(o.nodelist)}
        mask = [m and s in inside for m, s in zip(mask, short, strict=True)]
    v.masks[key] = mask
    return mask


def _headroom(value: Any) -> tuple[float, ...] | None:
    """gpuHeadroom as (CPUs, GB) per idle GPU, None for the node's share ([] or None)"""
    if _is_empty(value):
        return None
    flat: list[Any] = []

    def add(x: Any) -> None:
        if isinstance(x, (list, tuple)):
            for y in x:
                add(y)
        else:
            flat.append(x)

    add(value)
    if not flat or not all(_is_number(x) for x in flat):
        raise AlternativesError(
            f"gpuHeadroom must be a number or [CPUs, GB], not {value!r}.",
            "slurm_alternatives:entry",
        )
    h = [float(x) for x in flat]
    return (h[0], h[0]) if len(h) == 1 else (h[0], h[1])


def _cpu_job_room(v: _View, place: Place, gpus: float) -> tuple[list[float], list[float]]:
    """each node's CPUs and memory a job may take: CPU-only jobs leave each idle GPU its share
    (gpuHeadroom, or the node's CPUs and memory per GPU)"""
    if gpus > 0:
        return v.free_cpu, v.free_mem
    h = place.gpu_headroom
    cpu, mem = [], []
    for i in range(len(v.names)):
        idle = _mmax(0.0, v.free_gpu[i])
        if h is None:
            share = _mmax(1.0, v.gpu_tot[i])
            h_cpu, h_mem = v.cpu_tot[i] / share, v.mem_tot[i] / share
        else:
            h_cpu, h_mem = h[0], h[1] * 1024
        cpu.append(v.free_cpu[i] - idle * h_cpu)
        mem.append(v.free_mem[i] - idle * h_mem)
    return cpu, mem


def _place_room(v: _View, place: Place, shape: Shape) -> tuple[float, int]:
    """the jobs of a shape that a place can start now by the view, and the first node with
    room for one (-1)"""
    o = place.opts
    mem = _mmax(_shape_memory(place, shape), 1.0)
    need = shape.cpus * o.ntasks
    mask = _place_nodes(v, place)
    avail_cpu, avail_mem = _cpu_job_room(v, place, o.gpus)
    slots, first = 0.0, -1
    for i in range(len(v.names)):
        per = _mmin(_mfloor(_mdiv(avail_cpu[i], need)), _mfloor(_mdiv(avail_mem[i], mem)))
        if o.gpus > 0:
            per = _mmin(per, _mfloor(_mdiv(v.free_gpu[i], o.gpus)))
        if not (mask[i] and v.usable[i]):
            per = 0.0
        per = _mmax(0.0, per)
        slots += per
        if first < 0 and per > 0:
            first = i
    if slots > 0:
        # a QOS or association limit that a pending job of the user waits for here holds new
        # jobs too, and so do jobs of the user that wait here for resources
        limited = any(
            lp == o.partition and (not o.qos or lq == o.qos)
            for lp, lq in zip(v.limit_parts, v.limit_qos, strict=True)
        )
        if limited or o.partition in v.wait_parts:
            return 0.0, -1
        # CPU-only jobs keep out of partitions where GPU jobs of the user are pending
        if o.gpus == 0 and o.partition in v.gpu_parts:
            return 0.0, -1
    return slots, first


def _take_room(v: _View, place: Place, shape: Shape, j: int) -> None:
    if j >= 0:
        _take(v, j, shape.cpus * place.opts.ntasks, _shape_memory(place, shape), place.opts.gpus)


def _gpu_nodes_only(v: _View, place: Place) -> bool:
    mask = _place_nodes(v, place)
    nodes = [v.gpu_tot[i] for i, m in enumerate(mask) if m]
    return bool(nodes) and all(g > 0 for g in nodes)


def _takes_job(v: _View, place: Place) -> bool:
    """alternatives of GPU nodes only take CPU-only jobs only if they say so (cpuJobs)"""
    return place.cpu_jobs or place.opts.gpus > 0 or not _gpu_nodes_only(v, place)


def _job_key(job_id: int | float | str, array_id: int | float | str) -> str:
    return f"{int(job_id)}_{int(array_id)}"


# --- the session ---------------------------------------------------------------------------------


@dataclass
class _Entry:
    """a job of this session at a place: the places, its own CPUs at each, their time
    options, its place now (k) and shape, when it got there and when it was looked at"""

    places: tuple[Place, ...]
    cpus: tuple[int | float, ...]
    time_strs: tuple[str, ...]
    k: int
    shape: Shape
    placed: float
    looked: float


class Session:
    """The state of SlurmAlternatives in a process: the view of the cluster (read at most every
    view_age s), the jobs submitted to one of several places, the places Slurm refused, and the
    timing settings. ``slurm`` runs the commands; ``clock`` gives the time (s)."""

    def __init__(
        self, slurm: Slurm | None = None, clock: Callable[[], float] | None = None
    ) -> None:
        self._lock = threading.RLock()
        self.slurm = Slurm() if slurm is None else slurm
        self.clock = time.time if clock is None else clock
        self.view_age = VIEW_AGE
        self.min_wait = MIN_WAIT
        self.look_interval = LOOK_INTERVAL
        self.start_tolerance = START_TOLERANCE
        self.view: _View | None = None
        self.registry: dict[str, _Entry] = {}
        self.disabled: set[str] = set()

    def reset(self) -> None:
        """as a new session (MATLAB: clear slurm_alternatives)"""
        with self._lock:
            self.__init__()

    def set_timing(
        self,
        *,
        view_age: float | None = None,
        min_wait: float | None = None,
        look_interval: float | None = None,
        start_tolerance: float | None = None,
    ) -> None:
        """timing settings (s); the view is read again (MATLAB's 'set')"""
        with self._lock:
            if view_age is not None:
                self.view_age = float(view_age)
            if min_wait is not None:
                self.min_wait = float(min_wait)
            if look_interval is not None:
                self.look_interval = float(look_interval)
            if start_tolerance is not None:
                self.start_tolerance = float(start_tolerance)
            self.view = None

    # -- places -----------------------------------------------------------------------------

    def places(self, main: Mapping[str, Any], alternatives: Any = None) -> list[Place]:
        """the main place and the alternatives (MATLAB's 'places')

        main: SlurmParam, SlurmConstraint, memPerCPU, maxCPUNum, minCPUNum and gpuHeadroom of
        the call; alternatives: a list of entries (dicts), one entry, or empty. An entry
        without memPerCPU and --mem-per-cpu reads the view for its nodes' memory per CPU.
        """
        with self._lock:
            return self._places(main, alternatives)

    def _places(self, main: Mapping[str, Any], alternatives: Any) -> list[Place]:
        base = {name: main.get(name) for name in ENTRY_KEYS}
        if _is_empty(base["SlurmConstraint"]):
            base["SlurmConstraint"] = ""
        param = "" if base["SlurmParam"] is None else str(base["SlurmParam"])
        constraint = str(base["SlurmConstraint"])
        first = Place(
            param,
            constraint,
            _number(base["memPerCPU"], "memPerCPU", 0),
            _number(base["maxCPUNum"], "maxCPUNum", 0, 24.0),
            _number(base["minCPUNum"], "minCPUNum", 0, 1.0),
            _headroom(base["gpuHeadroom"]),
            True,  # (the main place takes any job, as without alternatives)
            parse_options(f"{param} {constraint}"),
        )
        out = [first]
        if _is_empty(alternatives):
            return out
        if isinstance(alternatives, Mapping):
            alternatives = [alternatives]
        elif not isinstance(alternatives, (list, tuple)):
            raise AlternativesError(
                "SlurmAlternatives must be a list of entries (JSON objects).",
                "slurm_alternatives:entry",
            )
        for a, alt in enumerate(alternatives, 1):
            if not isinstance(alt, Mapping):
                raise AlternativesError(
                    "Each entry of SlurmAlternatives must be a struct (a JSON object).",
                    "slurm_alternatives:entry",
                )
            out.append(self._entry_place(base, first, alt, a))
        return out

    def _entry_place(
        self, base: Mapping[str, Any], first: Place, alt: Mapping[str, Any], a: int
    ) -> Place:
        p = dict(base)
        p["memPerCPU"] = None
        p["cpuJobs"] = False
        for key, value in alt.items():
            name = next((n for n in ENTRY_KEYS if n.lower() == str(key).lower()), None)
            if name is None:
                _warn(
                    f"SlurmAlternatives entry {a}: the key {key} is ignored (an entry has "
                    f"{', '.join(ENTRY_KEYS)}).",
                    "slurm_alternatives:key",
                )
                continue
            p[name] = value
        if not isinstance(p["SlurmParam"], str):
            raise AlternativesError(
                f"SlurmAlternatives entry {a}: SlurmParam must be text.", "slurm_alternatives:entry"
            )
        constraint = p["SlurmConstraint"]
        if _is_empty(constraint):
            constraint = ""
        if not isinstance(constraint, str):
            raise AlternativesError(
                f"SlurmAlternatives entry {a}: SlurmConstraint must be text.",
                "slurm_alternatives:entry",
            )
        cpu_jobs = _unwrap(p["cpuJobs"])
        cpu_jobs = _is_number(cpu_jobs) and not isinstance(cpu_jobs, str) and bool(cpu_jobs == 1)
        opts = parse_options(f"{p['SlurmParam']} {constraint}")
        place = Place(
            p["SlurmParam"],
            constraint,
            _number(p["memPerCPU"], "memPerCPU", a),
            _number(p["maxCPUNum"], "maxCPUNum", a, first.max_cpu_num),
            _number(p["minCPUNum"], "minCPUNum", a, first.min_cpu_num),
            _headroom(p["gpuHeadroom"]),
            cpu_jobs,
            opts,
        )
        if place.mem_per_cpu is None and opts.mem_per_cpu is not None:
            place = replace(place, mem_per_cpu=opts.mem_per_cpu / 1024)
        if place.mem_per_cpu is None:
            # the memory per CPU of the place's nodes (the smallest of them)
            v = self._current_view()
            mask = _place_nodes(v, place)
            ratio = _NAN
            ratios = [_mdiv(v.mem_tot[i], v.cpu_tot[i]) for i, m in enumerate(mask) if m]
            for r in ratios:
                ratio = _mmin(ratio, r)
            ratio = ratio / 1024
            mem_per_cpu = first.mem_per_cpu
            if ratios and ratio > 0:
                mem_per_cpu = ratio
            place = replace(place, mem_per_cpu=mem_per_cpu)
        if not _shell_splittable(f"{place.slurm_param} {place.slurm_constraint}"):
            _warn(
                f'Slurm refuses jobs with "{place.slurm_param}" (its options can\'t be split as a '
                "shell splits them); this place is not used.",
                "slurm_alternatives:testOnly",
            )
            self.disabled.add(place.slurm_param)
        return place

    # -- the view -------------------------------------------------------------------------------

    def current_view(self) -> _View:
        """the view of the cluster, read again if older than view_age s"""
        with self._lock:
            return self._current_view()

    def _current_view(self) -> _View:
        t = self.clock()
        if self.view is not None and t - self.view.time < self.view_age:
            return self.view
        v = _View(t)
        lines = self.slurm.show_nodes()
        rows = None if lines is None else self.slurm.pending_jobs()
        if lines is None or rows is None:
            v.known = False
            print(
                "Warning: Slurm did not answer (scontrol show node, squeue); the jobs go to the "
                "main place and none moves until the cluster is read again.",
                file=sys.stderr,
            )
        else:
            _read_nodes(v, lines)
            _read_pending(v, rows)
            _start_pending(v)
        self.view = v
        return v

    # -- Slurm's check ---------------------------------------------------------------------------

    def _test_only_key(self, place: Place, shape: Shape, time_option: str) -> str:
        return f"{shape.param}|{place.slurm_constraint}|{_sprintf_d(shape.cpus)}|{time_option}"

    def _slurm_says_now(
        self, v: _View, place: Place, shape: Shape, time_option: str, is_main: bool
    ) -> bool:
        """sbatch --test-only (once per view for a place and shape): whether the job would
        start now; a place Slurm refuses is not used again in this session (the main place,
        and a place at an error of the moment, are)"""
        key = self._test_only_key(place, shape, time_option)
        if key in v.test_only:
            return v.test_only[key]
        now_ok = False
        try:
            r = self.slurm.test_only(
                cpus_per_task=int(shape.cpus),
                slurm_param=shape.param,
                slurm_constraint=place.slurm_constraint,
                time_str=time_option,
            )
        except ValueError as err:  # (options a shell can't split)
            out, transient, start = str(err), False, None
        else:
            out, start = r.result.output, r.start
            transient = not r.result.answered or any(w in out.lower() for w in TRANSIENT_WORDS)
        if start is None:
            if not transient and not is_main:
                self.disabled.add(place.slurm_param)
                _warn(
                    f'Slurm refuses jobs with "{place.slurm_param}" ({_strtrim(out)}); this '
                    "place is not used.",
                    "slurm_alternatives:testOnly",
                )
        else:
            now = datetime.fromtimestamp(self.clock())
            now_ok = (start - now).total_seconds() <= self.start_tolerance
        v.test_only[key] = now_ok
        return now_ok

    def _has_room(
        self, v: _View, place: Place, shape: Shape, time_option: str, is_main: bool
    ) -> tuple[bool, int]:
        """whether a place can start a job of a shape now, and the node it takes (-1)"""
        if not self._usable(place):
            return False, -1
        slots, j = _place_room(v, place, shape)
        if slots < 1 and place.opts.partition:
            return False, j
        ok = self._slurm_says_now(v, place, shape, time_option, is_main)
        # (a place without a partition: Slurm's check alone, one job per view)
        if ok and not place.opts.partition:
            v.test_only[self._test_only_key(place, shape, time_option)] = False
        return ok, j

    def _wait_place(self, v: _View, places: Sequence[Place]) -> int:
        """where a job waits when no place has room: the first place that isn't of GPU nodes
        only (for CPU-only jobs), else the main place"""
        for i, p in enumerate(places):
            if not self._usable(p):
                continue
            if p.opts.gpus > 0 or not _gpu_nodes_only(v, p):
                return i
        return 0

    def _usable(self, place: Place) -> bool:
        # not refused by Slurm, and its options split as a shell splits them
        return place.slurm_param not in self.disabled and _shell_splittable(
            f"{place.slurm_param} {place.slurm_constraint}"
        )

    # -- choice, records, moves ------------------------------------------------------------------

    def choose(
        self, places: Sequence[Place], own_cpus: Sequence[int | float], time_options: Sequence[str]
    ) -> tuple[int, int | float, str]:
        """the place (index, 0 the main one) for a job whose own CPUs at place k (as a config
        of only place k gives them) are own_cpus[k], its CPUs there and the SlurmParam of its
        shape (MATLAB's 'choose'; MATLAB's place index is this one + 1)"""
        with self._lock:
            v = self._current_view()
            current = _main_shape(places, own_cpus)
            k, shape = 0, _place_shape(places, own_cpus, 0, current)
            if v.known:
                k, shape = self._choose_place(v, places, own_cpus, time_options, current)
            # (a shape a shell can't split goes to the main place, as without alternatives)
            if k > 0 and not _shell_splittable(f"{shape.param} {places[k].slurm_constraint}"):
                k, shape = 0, _place_shape(places, own_cpus, 0, current)
            if k > 0:
                print(
                    f"The job goes to {places[k].name} ({_sprintf_d(shape.cpus)} CPUs): the main "
                    "place has no room for it now."
                )
            return k, _int(shape.cpus), shape.param

    def _choose_place(
        self,
        v: _View,
        places: Sequence[Place],
        own_cpus: Sequence[int | float],
        time_options: Sequence[str],
        current: Shape,
    ) -> tuple[int, Shape]:
        for k, place in enumerate(places):
            if not _takes_job(v, place):
                continue
            shape = _place_shape(places, own_cpus, k, current)
            ok, j = self._has_room(v, place, shape, time_options[k], k == 0)
            if ok:
                _take_room(v, place, shape, j)
                return k, shape
        k = self._wait_place(v, places)
        return k, _place_shape(places, own_cpus, k, current)

    def submitted(
        self,
        job_id: int,
        array_id: int,
        places: Sequence[Place],
        own_cpus: Sequence[int | float],
        k: int,
        time_options: Sequence[str],
    ) -> None:
        """records a job submitted to place k (index), so that it can move later"""
        with self._lock:
            t = self.clock()
            places, cpus = tuple(places), tuple(own_cpus)
            shape = _place_shape(places, cpus, k, _main_shape(places, cpus))
            self.registry[_job_key(job_id, array_id)] = _Entry(
                places, cpus, tuple(time_options), k, shape, t, t
            )

    def forget(self, job_ids: Sequence[int], array_ids: Sequence[int]) -> None:
        """drops the records of jobs that ended"""
        with self._lock:
            if not self.registry:
                return
            for job_id, array_id in zip(job_ids, array_ids, strict=True):
                self.registry.pop(_job_key(job_id, array_id), None)

    def move(self, job_ids: Sequence[int], array_ids: Sequence[int]) -> None:
        """recorded jobs that the caller found pending, that have waited in their place for
        min_wait s and weren't looked at in the last look_interval s, move to the first other
        place that has room for them now (one scontrol update each), as many as there is room
        for; other jobs are never looked at, and between looks no Slurm command runs"""
        with self._lock:
            if not self.registry:
                return
            t = self.clock()
            cands: list[str] = []
            for job_id, array_id in zip(job_ids, array_ids, strict=True):
                key = _job_key(job_id, array_id)
                e = self.registry.get(key)
                if (
                    e is not None
                    and key not in cands
                    and t - e.placed >= self.min_wait
                    and t - e.looked >= self.look_interval
                ):
                    cands.append(key)
            if not cands:
                return
            v = self._current_view()
            if not v.known:
                return
            counts = [0, 0, 0]  # room where they wait, moved, no room elsewhere
            for key in cands:
                self._look(v, key, t, counts)
            if sum(counts) > 0:
                n_fit, n_moved, n_stay = counts
                print(
                    f"{datetime.fromtimestamp(t):%H:%M:%S}: of {sum(counts)} waiting jobs, "
                    f"{n_fit} have room where they wait, {n_moved} moved, {n_stay} have no room "
                    "elsewhere."
                )

    def _look(self, v: _View, key: str, t: float, counts: list[int]) -> None:
        e = self.registry.get(key)
        if e is None:
            return
        e.looked = t
        p = v.pend_index.get(key)
        if p is None:
            return
        row = v.pend[p]
        places, j = e.places, e.k
        here = places[j]
        # only jobs that wait for resources or a QOS or association limit (not held, waiting
        # for a dependency or start time, or not yet seen by the scheduler), and that don't fit
        # where they are (on a node of their place; a CPU-only job at GPU nodes where GPU jobs
        # of the user are pending doesn't stay there)
        mask = _place_nodes(v, here)
        gpu_work = (
            here.opts.gpus == 0 and here.opts.partition in v.gpu_parts and _gpu_nodes_only(v, here)
        )
        if not gpu_work and row.node >= 0 and (mask[row.node] or not here.opts.partition):
            counts[0] += 1
            return
        if row.reason not in MOVABLE_REASONS and not _LIMIT_REASON.match(row.reason):
            return
        # a job that something else moved (e.g. a user's own script) is left to it
        if here.opts.partition and row.partition != here.opts.partition:
            del self.registry[key]
            return
        target, node, shape = -1, -1, e.shape
        for k, place in enumerate(places):
            if k == j or not _takes_job(v, place):
                continue
            shape = _place_shape(places, e.cpus, k, e.shape)
            if not move_fields(here, place, e.shape, shape):
                continue
            ok, node = self._has_room(v, place, shape, e.time_strs[k], k == 0)
            if ok:
                target = k
                break
        # a CPU-only job doesn't wait in the queue of GPU nodes
        if target < 0 and here.opts.gpus == 0 and _gpu_nodes_only(v, here):
            k = self._wait_place(v, places)
            shape = _place_shape(places, e.cpus, k, e.shape)
            if k != j and move_fields(here, places[k], e.shape, shape):
                target, node = k, -1
        if target < 0:
            counts[2] += 1
            return
        to = places[target]
        fields = move_fields(here, to, e.shape, shape)
        try:
            # (the fields as MATLAB's shell splits them: quotes are removed)
            r = self.slurm.update_job(key, " ".join(f"{k}={v}" for k, v in fields.items()))
        except ValueError as err:
            print(f"Job {key} was not moved to {to.name}: {err}")
            counts[2] += 1
            return
        if not r.ok:
            print(f"Job {key} was not moved to {to.name}: {_strtrim(r.output)}")
            if not r.answered:
                # (it may have moved: left alone for min_wait, the next look finds it)
                e.placed = t
            counts[2] += 1
            return
        print(
            f"Job {key} moved from {here.name} ({_sprintf_d(e.shape.cpus)} CPUs) to {to.name} "
            f"({_sprintf_d(shape.cpus)} CPUs)."
        )
        _take_room(v, to, shape, node)
        e.k, e.shape, e.placed = target, shape, t
        counts[1] += 1


def _number(value: Any, what: str, a: int, default: float | None = None) -> float | None:
    """a number of an entry (a list of one is its value; empty: default)"""
    value = _unwrap(value)
    if _is_empty(value):
        return default
    if not _is_number(value) or isinstance(value, str):
        where = "the main settings" if a == 0 else f"SlurmAlternatives entry {a}"
        raise AlternativesError(
            f"{where}: {what} must be a number, not {value!r}.", "slurm_alternatives:entry"
        )
    return float(value)


# --- the session of the process ------------------------------------------------------------------

_SESSION = Session()


def default_session() -> Session:
    """the session of this process (MATLAB's persistent variables)"""
    return _SESSION


def reset() -> None:
    """forget the view, the recorded jobs, the refused places and the timing settings"""
    _SESSION.reset()


def move(job_ids: Sequence[int], array_ids: Sequence[int]) -> None:
    """Session.move of the process's session (for single jobs' status checks)"""
    _SESSION.move(job_ids, array_ids)


def forget(job_ids: Sequence[int], array_ids: Sequence[int]) -> None:
    """Session.forget of the process's session"""
    _SESSION.forget(job_ids, array_ids)


# --- the cluster loop's places ---------------------------------------------------------------


def _split_keys(job_keys: Sequence[str]) -> tuple[list[int], list[int]]:
    jobs, arrays = [], []
    for key in job_keys:
        job, _, array = str(key).partition("_")
        if job.isdigit() and array.isdigit():
            jobs.append(int(job))
            arrays.append(int(array))
    return jobs, arrays


class ClusterPlaces:
    """the places of a cluster call's jobs (the protocol of jobs.cluster._OnePlace): each
    submission's place, CPUs, SlurmParam, SlurmConstraint and time option from ``choose``; the
    job recorded by ``submitted``; pending jobs looked at by ``move``; ``forget`` at the end"""

    def __init__(
        self,
        session: Session,
        places: Sequence[Place],
        time_options: Sequence[str],
        slot_mode: bool,
        gpu_job: bool,
    ) -> None:
        self.session = session
        self.places = list(places)
        self.time_options = list(time_options)
        self.slot_mode = slot_mode
        self.gpu_job = gpu_job
        self._chosen: dict[int, tuple[int, list[int | float]]] = {}
        self._last: tuple[int, list[int | float]] | None = None

    def own_cpus(self, trial: int) -> list[int | float]:
        """each place's CPUs for a job's try after trial earlier ones: CPU jobs x (trial + 1)
        and GPU jobs as they are (in-job slots: the other way round), within maxCPUNum"""
        if self.slot_mode:
            factor = trial + 1 if self.gpu_job else 1
        else:
            factor = 1 if self.gpu_job else trial + 1
        return [_int(_mmin(p.max_cpu_num, float(p.cpus_per_task) * factor)) for p in self.places]

    def choose(self, cpus: int, trial: int, batch: int | None = None) -> tuple[int, str, str, str]:
        """(CPUs, SlurmParam, SlurmConstraint, time option) of a job's try; cpus: its CPUs at
        the main place"""
        own = self.own_cpus(trial)
        own[0] = cpus
        k, c, param = self.session.choose(self.places, own, self.time_options)
        if batch is None:
            self._last = (k, own)
        else:
            self._chosen[batch] = (k, own)
        return int(c), param, self.places[k].slurm_constraint, self.time_options[k]

    def submitted(self, job_id: int, array_id: int, trial: int, batch: int | None = None) -> None:
        """records the job of the last choice (of batch); a job whose choice isn't known is
        not recorded, and so never moves"""
        if batch is None:
            choice, self._last = self._last, None
        else:
            choice = self._chosen.pop(batch, None)
        if choice is not None:
            k, own = choice
            self.session.submitted(job_id, array_id, self.places, own, k, self.time_options)

    def move(self, job_keys: Sequence[str]) -> None:
        self.session.move(*_split_keys(job_keys))

    def forget(self, job_keys: Sequence[str]) -> None:
        self.session.forget(*_split_keys(job_keys))


def cluster_places(
    call: Any, shape: Any, time_option: str, requested_cpus: float | None = None
) -> ClusterPlaces | None:
    """The places of a cluster call (generic_computing_frameworks_wrapper's slurmPlaces), or
    None for one place. shape: the call's jobs.cluster.JobShape; time_option: the main place's
    time limit option; requested_cpus: the CPUs the call asked for (cpus_per_task before the
    main place's shaping; without it, the main place's CPUs).

    Each place gets the CPUs a config of only that place gives a job (job_cpus from
    memAllocate and its memPerCPU, minCPUNum, maxCPUNum, rounded), the --nice of the call (in
    its SlurmParam, not in its options: with a nice, the main place and the alternatives don't
    trade jobs, as in MATLAB), and the time option of its SlurmParam."""
    s = call.settings
    alternatives = s.get("slurm_alternatives")
    if _is_empty(alternatives) or s.get("whole_node_job"):
        return None
    session = default_session()
    main = {
        "SlurmParam": s.get("slurm_param", ""),
        "SlurmConstraint": s.get("slurm_constraint", ""),
        "memPerCPU": s.get("mem_per_cpu", 20.9),
        "maxCPUNum": s.get("max_cpu_num", 24),
        "minCPUNum": s.get("min_cpu_num", 1),
        "gpuHeadroom": s.get("gpu_headroom"),
    }
    places = session.places(main, alternatives)
    if len(places) < 2:
        return None
    request = s.get("cpus_per_task", 1) if requested_cpus is None else requested_cpus
    nice = s.get("nice", 0) or 0
    shaped = []
    for k, p in enumerate(places):
        if k == 0:
            cpus: float = shape.cpus_per_task
        else:
            cpus = job_cpus(
                request, s.get("mem_allocate"), p.mem_per_cpu, p.min_cpu_num, p.max_cpu_num, True
            )
        param = p.slurm_param
        if nice > 0 and "--nice=" not in param:
            param = f"{param} --nice={_round_half_away(nice)}"
        shaped.append(replace(p, slurm_param=param, cpus_per_task=_int(float(cpus))))
    fmt = TIME_FORMAT_MCC if shape.slot_mode else TIME_FORMAT_BATCH
    options = time_strs(shaped, s.get("job_time_limit"), fmt, time_option)
    return ClusterPlaces(session, shaped, options, shape.slot_mode, bool(s.get("gpu_job")))


# --- single jobs ---------------------------------------------------------------------------------


class SingleJobPlaces:
    """the places of a single job (generic_single_job_submit_wrapper): ``choose`` gives the
    place index, CPUs, SlurmParam (with the call's --nice), SlurmConstraint and time option of
    a submission, ``submitted(job_id, array_id, k)`` records it; the status checks call the
    module's ``move`` and ``forget``"""

    def __init__(
        self,
        session: Session,
        places: Sequence[Place],
        own_cpus: Sequence[int | float],
        time_options: Sequence[str],
        nice: float,
    ) -> None:
        self.session = session
        self.places = list(places)
        self.own_cpus = list(own_cpus)
        self.time_options = list(time_options)
        self.nice = nice

    def choose(self) -> tuple[int, int | float, str, str, str]:
        k, cpus, param = self.session.choose(self.places, self.own_cpus, self.time_options)
        if self.nice > 0 and "--nice=" not in param:
            param = f"{param} --nice={_round_half_away(self.nice)}"
        return k, cpus, param, self.places[k].slurm_constraint, self.time_options[k]

    def submitted(self, job_id: int, array_id: int, k: int) -> None:
        self.session.submitted(job_id, array_id, self.places, self.own_cpus, k, self.time_options)


def single_job_places(
    *,
    alternatives: Any,
    slurm_param: str,
    slurm_constraint: str,
    mem_per_cpu: float,
    max_cpu_num: float,
    min_cpu_num: float,
    gpu_headroom: Any,
    requested_cpus: float,
    mem_allocate: float | None,
    main_cpus: int,
    job_time_limit: float | None,
    main_time_option: str,
    nice: float = 0,
    session: Session | None = None,
) -> SingleJobPlaces | None:
    """The places of a single job, or None for one place. slurm_param: the call's SlurmParam
    without its --nice; requested_cpus: cpusPerTask as asked; main_cpus: the job's CPUs at the
    main place; main_time_option: its time limit option (shapes.TIME_FORMAT_SINGLE). A place's
    own CPUs are min(maxCPUNum, max(minCPUNum, max(cpusPerTask, ceil(memAllocate / memPerCPU))))
    of that place (single jobs aren't rounded)."""
    if _is_empty(alternatives):
        return None
    session = default_session() if session is None else session
    main = {
        "SlurmParam": slurm_param,
        "SlurmConstraint": slurm_constraint,
        "memPerCPU": mem_per_cpu,
        "maxCPUNum": max_cpu_num,
        "minCPUNum": min_cpu_num,
        "gpuHeadroom": gpu_headroom,
    }
    places = session.places(main, alternatives)
    if len(places) < 2:
        return None
    mem = 0.0 if mem_allocate is None else _mmax(float(_unwrap(mem_allocate)), 0.0)
    own: list[int | float] = [main_cpus]
    for p in places[1:]:
        per_mem = _mceil(_mdiv(mem, _NAN if p.mem_per_cpu is None else p.mem_per_cpu))
        cpus = _mmin(p.max_cpu_num, _mmax(p.min_cpu_num, _mmax(float(requested_cpus), per_mem)))
        own.append(_int(cpus))
    options = time_strs(places, job_time_limit, TIME_FORMAT_SINGLE, main_time_option)
    return SingleJobPlaces(session, places, own, options, nice)
