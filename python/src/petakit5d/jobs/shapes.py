"""The shape of a Slurm job as MATLAB's job framework makes it: CPUs, retries, array ids and the
time limit option.

Values follow MATLAB's double arithmetic (IEEE double, as Python's float); integral results are
returned as int.
"""

from __future__ import annotations

import math
from typing import Any

__all__ = [
    "ARRAY_ID_MODULUS",
    "TIME_FORMAT_BATCH",
    "TIME_FORMAT_MCC",
    "TIME_FORMAT_SINGLE",
    "array_id",
    "has_time_option",
    "job_cpus",
    "retry_cpus",
    "single_job_cpus",
    "time_str",
]

# --array=<batch mod 5000>: Slurm's MaxArraySize is 5001 by default
ARRAY_ID_MODULUS = 5000

# the -t option of the three MATLAB submitters (hours, minutes): the batch framework,
# the compiled-mode framework and the single-job wrapper
TIME_FORMAT_BATCH = " -t {h}:{m:02}:00 "
TIME_FORMAT_MCC = " -t {h:02}:{m:02}:00 "
TIME_FORMAT_SINGLE = " -t {h}:{m}:00 "


def _num(x: float) -> int | float:
    return int(x) if isinstance(x, float) and x.is_integer() else x


def _rem(x: float, y: float) -> float:
    """MATLAB rem (sign of x); rem(x, 0) is NaN"""
    if y == 0 or math.isinf(x) or math.isnan(x) or math.isnan(y):
        return math.nan
    if math.isinf(y):
        return x
    return math.fmod(x, y)


def _prime_factor_count(n: float) -> int:
    """numel(factor(n)) for a positive integer n (factor(1) is 1, one element)"""
    if n != n or n < 0 or not float(n).is_integer():
        raise ValueError(f"factor: {n} is not a nonnegative integer")
    n = int(n)
    if n < 4:
        return 1
    count, p = 0, 2
    while p * p <= n:
        while n % p == 0:
            n //= p
            count += 1
        p += 1
    return count + (1 if n > 1 else 0)


def _div(a: float, b: float) -> float:
    """MATLAB's a / b for doubles (x / 0 is Inf, -Inf or NaN)"""
    if b == 0:
        if a == 0 or a != a:
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


def _ceil(x: float) -> float:
    return x if math.isinf(x) or math.isnan(x) else float(math.ceil(x))


def _mmin(a: float, b: float) -> float:
    """MATLAB min of two scalars (NaN is ignored)"""
    if a != a:
        return b
    if b != b:
        return a
    return b if b < a else a


def _mmax(a: float, b: float) -> float:
    if a != a:
        return b
    if b != b:
        return a
    return b if b > a else a


def _empty(x: Any) -> bool:
    return x is None or (hasattr(x, "__len__") and len(x) == 0)


def _scalar(x: Any) -> float:
    if hasattr(x, "__len__"):
        if len(x) != 1:
            raise ValueError(f"expected a scalar, got {x!r}")
        x = x[0]
    return float(x)


def job_cpus(
    cpus_per_task: float,
    mem_allocate: float | None,
    mem_per_cpu: float,
    min_cpu_num: float,
    max_cpu_num: float,
    round_cpus: bool,
) -> int | float:
    """job_cpus of generic_computing_frameworks_wrapper.m

    Enough CPUs for mem_allocate GB at mem_per_cpu GB per CPU (None: no memory request),
    within [min_cpu_num, max_cpu_num]; with round_cpus (cluster jobs, not of whole nodes),
    rounded up toward a factor of max_cpu_num (at most doubled), or to max_cpu_num above 75%
    of it.
    """
    cpus = float(cpus_per_task)
    mem_per_cpu = float(mem_per_cpu)
    if not _empty(mem_allocate):
        mem = _scalar(mem_allocate)
        if cpus * mem_per_cpu < mem:
            cpus = _ceil(_div(mem, mem_per_cpu))
    max_cpu = float(max_cpu_num)
    cpus = _mmin(max_cpu, _mmax(float(min_cpu_num), cpus))
    if round_cpus and _rem(max_cpu, cpus) != 0 and _prime_factor_count(max_cpu) > 1:
        cpus_orig = cpus
        if _div(max_cpu, cpus) > 2:
            while _rem(max_cpu, cpus) > 1:
                if _div(cpus, cpus_orig) >= 2:
                    break
                cpus += 1
        elif _div(cpus, max_cpu) > 0.75:
            cpus = max_cpu
    return _num(cpus)


def single_job_cpus(
    cpus_per_task: float,
    mem_allocate: float | None,
    mem_per_cpu: float,
    min_cpu_num: float,
    max_cpu_num: float,
) -> int | float:
    """the CPUs of a single job (generic_single_job_submit_wrapper): job_cpus without rounding"""
    return job_cpus(cpus_per_task, mem_allocate, mem_per_cpu, min_cpu_num, max_cpu_num, False)


def retry_cpus(cpus_per_task: float, max_cpu_num: float, trial: int, gpu_job: bool) -> int | float:
    """the CPUs of a batch job's try after trial earlier ones (slurm_cluster_..._wrapper)

    CPU jobs get cpus x (trial + 1), so the memory of --mem-per-cpu grows with each try; GPU
    jobs keep their CPUs. Never above max_cpu_num.
    """
    factor = 1 if gpu_job else trial + 1
    return _num(_mmin(float(max_cpu_num), float(cpus_per_task) * factor))


def array_id(batch: int) -> int:
    """the --array index of batch (1-based): rem(batch, 5000)"""
    return int(batch) % ARRAY_ID_MODULUS


def has_time_option(slurm_param: str) -> bool:
    """whether SlurmParam sets the time limit itself (' -t ' or '--time', as MATLAB checks)"""
    return " -t " in slurm_param or "--time" in slurm_param


def _mformat_int(x: float) -> str:
    """sprintf('%d') of a double: integral values as integers, Inf/NaN by name"""
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return str(int(x))


def time_str(
    job_time_limit: float | None, slurm_param: str = "", fmt: str = TIME_FORMAT_BATCH
) -> str:
    """the time limit option for jobTimeLimit hours ('' for none, or if SlurmParam has one)

    Rounded to minutes, at least one minute; as MATLAB, 1.999 h gives ' -t 1:60:00 '. fmt is
    one of the TIME_FORMAT_* strings (fields h and m).
    """
    if _empty(job_time_limit) or has_time_option(slurm_param):
        return ""
    t = _mmax(_scalar(job_time_limit), 1 / 60)
    h = math.floor(t) if math.isfinite(t) else t
    m = _matlab_round((t - h) * 60) if math.isfinite(t) else math.nan
    if math.isfinite(h) and math.isfinite(m):
        return fmt.format(h=int(h), m=int(m))
    # (sprintf with Inf/NaN fields; the padded forms print them unpadded, as %02d does)
    return (
        fmt.replace("{h:02}", "{h}")
        .replace("{m:02}", "{m}")
        .format(h=_mformat_int(h), m=_mformat_int(m))
    )


def _matlab_round(x: float) -> float:
    """MATLAB round of a double: half away from zero (stand-in for petakit5d.matlabcompat)"""
    a = abs(x)
    r = math.floor(a)
    if a - r >= 0.5:
        r += 1
    return math.copysign(r, x)
