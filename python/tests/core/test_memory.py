import math
import os
import re
import subprocess
import sys

import pytest

from petakit5d.core import memory

GB = 2**30


def fake_fs(files: dict):
    def read(path):
        if path not in files:
            raise OSError(path)
        return files[path]

    return read


def no_run(cmd):
    raise AssertionError(f"unexpected command {cmd}")


V2_CGROUP = "0::/system.slice/slurmstepd.scope/job_77/step_batch/user/task_0\n"
V1_CGROUP = (
    "12:pids:/slurm/uid_1/job_77/step_batch\n"
    "11:memory:/slurm/uid_1/job_77/step_batch/task_0\n"
    "4:cpu,cpuacct:/slurm\n"
)


def test_slurm_job_memory_is_in_gb():
    # SystemMemoryTest/slurmJobMemoryIsInGB: --mem gives SLURM_MEM_PER_NODE in MB
    env = {"SLURM_JOB_ID": "123", "SLURM_MEM_PER_NODE": "32768"}
    assert memory.linux_memory(env, fake_fs({"/proc/self/cgroup": V2_CGROUP}), no_run) == 32


def test_cgroup_v2_limit_of_the_job_wins():
    # SystemMemoryTest/memoryIsNotTheValueOfAParentJob: an inherited SLURM_MEM_PER_NODE is not taken
    base = "/sys/fs/cgroup/system.slice/slurmstepd.scope/job_77"
    files = {
        "/proc/self/cgroup": V2_CGROUP,
        base + "/memory.max": f"{20 * GB}\n",
        base + "/step_batch/memory.max": "max\n",
        base + "/step_batch/user/memory.max": f"{24 * GB}\n",
        base + "/step_batch/user/task_0/memory.max": "max\n",
    }
    env = {"SLURM_JOB_ID": "77", "SLURM_MEM_PER_NODE": "1"}
    assert memory.linux_memory(env, fake_fs(files), no_run) == 20


def test_cgroup_v1_limit_and_no_limit():
    base = "/sys/fs/cgroup/memory/slurm/uid_1/job_77"
    files = {
        "/proc/self/cgroup": V1_CGROUP,
        base + "/memory.limit_in_bytes": f"{8 * GB}",
        base + "/step_batch/memory.limit_in_bytes": "9223372036854771712",
        base + "/step_batch/task_0/memory.limit_in_bytes": "9223372036854771712",
    }
    assert memory.job_cgroup_memory("77", fake_fs(files)) == 8
    files[base + "/memory.limit_in_bytes"] = "9223372036854771712"
    assert math.isnan(memory.job_cgroup_memory("77", fake_fs(files)))
    # the path must be this job's (job_7 is not job_77)
    assert math.isnan(memory.job_cgroup_memory("7", fake_fs(files)))


def test_memory_per_node_next_to_memory_per_cpu_is_not_taken(capsys):
    # SystemMemoryTest/memoryPerNodeNextToMemoryPerCPUIsNotTaken: scontrol decides
    env = {
        "SLURM_JOB_ID": "123",
        "SLURM_MEM_PER_NODE": "32768",
        "SLURM_MEM_PER_CPU": "1024",
        "SLURM_ARRAY_TASK_ID": "4",
    }
    seen = []

    def run(cmd):
        seen.append(cmd)
        return "JobId=123 ReqTRES=cpu=4,mem=16G,node=1 AllocTRES=cpu=4,mem=24G"

    assert memory.linux_memory(env, fake_fs({"/proc/self/cgroup": V2_CGROUP}), run) == 16
    assert seen == [["scontrol", "show", "job", "123_4"]]
    for out, value in (("mem=1.5T", 1536), ("mem=512M", 0.5), ("mem=1048576K", 1)):
        assert memory.linux_memory(env, fake_fs({}), lambda c, o=out: o) == value
    assert math.isnan(
        memory.linux_memory(
            env, fake_fs({}), lambda c: "slurm_load_jobs error: Invalid job id specified"
        )
    )
    assert "No match found" in capsys.readouterr().out


def test_without_slurm_memtotal():
    files = {"/proc/meminfo": "MemTotal:       527033388 kB\nMemFree: 1 kB\n"}
    assert memory.linux_memory({}, fake_fs(files), no_run) == 527033388 / 1024**2


def test_memory_is_a_number():
    # SystemMemoryTest/memoryIsANumber, on every system CI runs
    m = memory.get_system_memory()
    assert isinstance(m, float) and m > 0


@pytest.mark.skipif(sys.platform != "darwin", reason="macOS: hw.memsize")
def test_macos_memory_is_hw_memsize():
    out = subprocess.run(
        ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=True
    ).stdout
    assert memory.get_system_memory() == int(out) / 1024**3


@pytest.mark.skipif(
    not sys.platform.startswith("linux") or not os.environ.get("SLURM_JOB_ID"),
    reason="in a Slurm job on Linux",
)
def test_in_a_slurm_job_the_cgroup_limit_is_used():
    txt = open("/proc/self/cgroup").read()
    job = re.search(r"/job_(\d+)(/|\s)", txt)
    if not job:
        pytest.skip("no cgroup of the job")
    env = dict(os.environ, SLURM_JOB_ID=job.group(1), SLURM_MEM_PER_NODE="1")
    assert memory.linux_memory(env) > 0.5
