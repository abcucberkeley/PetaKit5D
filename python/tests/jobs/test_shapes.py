"""Job shapes: CPUs from memory (job_cpus), retries, array ids and the time limit option.

Ported from the MATLAB suite: SingleJobSubmitWrapperTest/cpusPerTaskIsClampedToLimits (its four
cases) and the CPU counts of the config tests (test_config.py). The grids below are oracle cases:
MATLAB R2024b computed them with job_cpus of generic_computing_frameworks_wrapper.m and
time_strings of slurm_alternatives.m (both local functions, run verbatim; 12 framework calls with
the fake Slurm gave the same CPUs as the extracted job_cpus), and the hashes are of its results.
"""

import hashlib
import math
import struct

import pytest

from petakit5d.jobs.shapes import (
    TIME_FORMAT_BATCH,
    TIME_FORMAT_MCC,
    TIME_FORMAT_SINGLE,
    array_id,
    has_time_option,
    job_cpus,
    retry_cpus,
    single_job_cpus,
    time_str,
)

# --- oracle grids (stand-in for the F04 harness: inputs here, MATLAB's results as hashes) ---

JOB_CPUS_GRID = {
    "cpus": list(range(1, 49)),
    "mems": [
        None,
        0,
        0.5,
        1,
        5,
        10,
        20,
        20.9,
        21,
        41.8,
        41.81,
        50,
        62.7,
        100,
        150,
        200,
        300,
        500,
        1000,
        2000,
        5000,
    ],
    "mem_per_cpu": [20.9, 7.6],
    "max_cpu": [1, 4, 16, 20, 24, 32, 48, 128],
    "min_cpu": [1, 4],
}
# MATLAB's 64512 values in the order below, as int32 little-endian
JOB_CPUS_SHA256 = "c37a8e1abd0a5387e5551b83a099a74e39e020d5ef4621c22ad71f4309238c57"

TIME_LIMITS = [
    None, 0, 1e-3, 1 / 120, 1 / 60, 0.0083, 0.00834, 0.5, 1, 1.5, 1.9999, 1.99, 2.0083, 2.0084,
    12, 23.99, 24, 100.5, 0.25, 1 / 3, 2 / 3, 7.0125, math.nan, math.inf,
]  # fmt: skip
TIME_PARAMS = [
    "-p q -n1",
    "-p q -t 5",
    "-p q --time=1:00",
    " -t 1:00 -p q",
    "-t 1:00 -p q",
    "-p q -t5",
]
TIME_FORMATS = [TIME_FORMAT_BATCH, TIME_FORMAT_MCC, TIME_FORMAT_SINGLE]
# MATLAB's 432 strings (format fastest, then SlurmParam, then limit) joined by newlines
TIME_SHA256 = "ec36107cdb43966d7c5e181faa22f495ea183893fec2ad8b6d335a29f0eedea6"


def job_cpus_grid() -> list[int]:
    g = JOB_CPUS_GRID
    return [
        job_cpus(c, m, mpc, mn, mx, r)
        for r in (False, True)
        for mx in g["max_cpu"]
        for mn in g["min_cpu"]
        for mpc in g["mem_per_cpu"]
        for m in g["mems"]
        for c in g["cpus"]
    ]


def time_grid() -> list[str]:
    return [time_str(t, p, f) for t in TIME_LIMITS for p in TIME_PARAMS for f in TIME_FORMATS]


@pytest.mark.oracle
def test_job_cpus_equals_matlab_on_the_grid():
    values = job_cpus_grid()
    assert len(values) == 64512
    assert all(isinstance(v, int) for v in values)
    blob = b"".join(struct.pack("<i", v) for v in values)
    assert hashlib.sha256(blob).hexdigest() == JOB_CPUS_SHA256


@pytest.mark.oracle
def test_time_option_equals_matlab_on_the_grid():
    values = time_grid()
    assert len(values) == 432
    assert hashlib.sha256("\n".join(values).encode()).hexdigest() == TIME_SHA256


# --- SingleJobSubmitWrapperTest/cpusPerTaskIsClampedToLimits (defaults: 1 CPU, 20.9 GB/CPU, 1..24)


@pytest.mark.parametrize(
    "kwargs, cpus",
    [
        ({"cpus_per_task": 40}, 24),  # cpusAboveMax
        ({"min_cpu_num": 4}, 4),  # minCPUNum
        ({"mem_allocate": 100}, 5),  # memAllocate
        ({"mem_allocate": 1000}, 24),  # memAllocateAboveMax
    ],
    ids=["cpusAboveMax", "minCPUNum", "memAllocate", "memAllocateAboveMax"],
)
def test_cpus_per_task_is_clamped_to_limits(kwargs, cpus):
    args = {
        "cpus_per_task": 1,
        "mem_allocate": None,
        "mem_per_cpu": 20.9,
        "min_cpu_num": 1,
        "max_cpu_num": 24,
    }
    assert single_job_cpus(**{**args, **kwargs}) == cpus


# --- the rounding of batch jobs (values checked by the MATLAB grid as well)


@pytest.mark.parametrize(
    "cpus, mem, mem_per_cpu, max_cpu, expected",
    [
        (1, 50, 5, 24, 12),  # 10 -> 12, a factor of 24
        (1, 100, 20.9, 24, 6),  # 5 -> 6
        (7, None, 20.9, 24, 8),  # 7 -> 8
        (9, None, 20.9, 24, 12),  # 9: 24 / 9 > 2, rounded to 12
        (17, None, 20.9, 24, 17),  # 17: 24 / 17 < 2 and 17 / 24 < 0.75: kept
        (19, None, 20.9, 24, 24),  # above 75% of 24
        (13, None, 20.9, 16, 16),
        (5, None, 20.9, 23, 5),  # a prime maxCPUNum: no rounding
        (1, 0.5, 20.9, 1, 1),
        (3, None, 20.9, 128, 4),
        (3, None, 20.9, 14, 6),  # 3..6 all leave a remainder > 1: stops at twice the CPUs
    ],
)
def test_batch_job_cpus_round_to_a_factor_of_max_cpu_num(cpus, mem, mem_per_cpu, max_cpu, expected):
    assert job_cpus(cpus, mem, mem_per_cpu, 1, max_cpu, True) == expected


def test_empty_memory_requests_and_edges():
    assert job_cpus(2, [], 20.9, 1, 24, True) == 2
    assert job_cpus(2, [50], 20.9, 1, 24, False) == 3
    assert job_cpus(1, 0, 20.9, 1, 24, True) == 1
    # MATLAB's x / 0: Inf CPUs, clamped to maxCPUNum
    assert job_cpus(1, 50, 0, 1, 24, False) == 24
    assert job_cpus(0, None, 20.9, 0, 24, True) == 0


def test_retry_cpus():
    assert [retry_cpus(4, 24, t, False) for t in range(4)] == [4, 8, 12, 16]
    assert retry_cpus(12, 24, 2, False) == 24
    assert [retry_cpus(4, 24, t, True) for t in range(3)] == [4, 4, 4]
    assert retry_cpus(40, 24, 0, True) == 24


def test_array_id():
    assert [array_id(b) for b in (1, 4999, 5000, 5001, 12345)] == [1, 4999, 0, 1, 2345]


def test_time_option():
    assert time_str(24) == " -t 24:00:00 "
    assert time_str(12, fmt=TIME_FORMAT_MCC) == " -t 12:00:00 "
    assert time_str(0.5, fmt=TIME_FORMAT_MCC) == " -t 00:30:00 "
    assert time_str(2.0833, fmt=TIME_FORMAT_SINGLE) == " -t 2:5:00 "
    assert time_str(0) == " -t 0:01:00 "  # at least one minute
    assert time_str(1.9999) == " -t 1:60:00 "  # as MATLAB rounds the minutes
    assert time_str(None) == "" and time_str([]) == ""
    assert time_str(12, "-p q --time=2:00:00") == ""
    assert time_str(12, "-p q -t 2:00:00") == ""
    # (MATLAB looks for ' -t ' only: a SlurmParam that starts with -t gets a second one)
    assert time_str(12, "-t 2:00:00 -p q") == " -t 12:00:00 "
    assert has_time_option(" -t 1:00") and not has_time_option("-p q -t5")

    assert time_str(math.inf) == " -t Inf:NaN:00 "  # (what MATLAB's sprintf makes of it)
    # half minutes round up (MATLAB's round; 7.125 h is 7 h 7.5 min exactly)
    assert time_str(7.125) == " -t 7:08:00 "
    assert time_str(0.125, fmt=TIME_FORMAT_SINGLE) == " -t 0:8:00 "
