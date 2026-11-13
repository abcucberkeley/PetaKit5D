"""Set-up of the real-Slurm tests of the job framework (SlurmClusterTest and the soak test).

They submit real jobs, so they run only when PETAKIT5D_TEST_SLURM_CONFIG names a cluster
config (marker ``cluster``; Linux and macOS):

  PETAKIT5D_TEST_SLURM_CONFIG      JSON cluster config passed to the framework as config_file
                                   (SlurmParam, maxCPUNum, memPerCPU, PythonLaunchStr, ...), e.g.
                                   tests/ci/github_actions_slurm.json. A config's values win over
                                   the caller's, so it must not set what the tests set: jobLogDir,
                                   SlurmConstraint, jobTimeLimit, maxTrialNum, maxJobNum,
                                   queryInterval, cpusPerTask, taskBatchNum, masterCompute.
  PETAKIT5D_TEST_SLURM_SCRATCH     folder on storage that the compute nodes share, for the test
                                   data and job logs (default: pytest's temporary folder, which
                                   only works where jobs run on this machine, as on the CI's
                                   single-node cluster)
  PETAKIT5D_TEST_SLURM_JOB_PREFIX  prefix of the tests' job names (default pk5dtest_)
  PETAKIT5D_TEST_SLURM_SOAK        the soak test's size, <tasks>x<jobs> (e.g. 2000x200); without
                                   it the soak test is skipped

Each test names its jobs ``<prefix><8 hex digits>`` (SlurmConstraint ``-J <name>``) and cancels
the jobs of that name and user when it ends, also when it fails. A watchdog cancels them every
20 s once a run has taken longer than its limit, so that the framework runs out of tries and
returns.
"""

from __future__ import annotations

import getpass
import json
import os
import re
import shutil
import threading
import time
import uuid
from collections.abc import Iterator, Sequence
from pathlib import Path

import pytest

from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.slurm import Slurm
from petakit5d.jobs.tasks import Task
from petakit5d.testing import slurm_fake

real_slurm = slurm_fake.real_slurm  # (the fixture)

JOB_PREFIX_VARIABLE = "PETAKIT5D_TEST_SLURM_JOB_PREFIX"
DEFAULT_JOB_PREFIX = "pk5dtest_"
WATCHDOG_MINUTES = 10.0
WATCHDOG_PERIOD = 20.0  # s between the watchdog's cancels
LEAVE_WAIT = 120.0  # s cancelled jobs get to leave the queue

# a task as bash: sleeps, then writes the Slurm environment of its job to its output (by a
# temporary name, so that the output appears only when complete); or records the attempt and
# fails without output
BASH_TASK = """sleep {sleep}
limit=$(scontrol show job "${{SLURM_JOB_ID:-0}}" 2>/dev/null | grep -o "TimeLimit=[^ ]*" \\
    | cut -d= -f2)
printf "array_job=%s\\ncpus=%s\\ntime_limit=%s\\n" "${{SLURM_ARRAY_JOB_ID:-}}" \\
    "${{SLURM_CPUS_PER_TASK:-}}" "$limit" > "$1.part"
mv "$1.part" "$1"
"""
BASH_FAIL = """echo "cpus=${{SLURM_CPUS_PER_TASK:-}}" >> {attempts}
exit 1
"""

# the same tasks in Python, a module of the test's folder
PYTHON_TASKS = '''"""tasks of the real-Slurm tests"""

import os
import re
import subprocess
import time


def _time_limit():
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        return ""
    try:
        r = subprocess.run(["scontrol", "show", "job", job], capture_output=True, text=True,
                           timeout=60)
    except (OSError, subprocess.SubprocessError):
        return ""
    m = re.search(r"TimeLimit=(\\S+)", r.stdout)
    return m.group(1) if m else ""


def _publish(text, path):
    with open(path + ".part", "w") as f:
        f.write(text)
    os.replace(path + ".part", path)


def record(output, sleep_seconds=0, final=""):
    """writes the Slurm environment of the job to output (and to final, if given)"""
    time.sleep(sleep_seconds)
    text = (f"array_job={os.environ.get('SLURM_ARRAY_JOB_ID', '')}\\n"
            f"cpus={os.environ.get('SLURM_CPUS_PER_TASK', '')}\\n"
            f"time_limit={_time_limit()}\\n")
    _publish(text, output)
    if final:
        _publish(text, final)


def fail(attempts):
    """records the attempt and fails"""
    with open(attempts, "a") as f:
        f.write(f"cpus={os.environ.get('SLURM_CPUS_PER_TASK', '')}\\n")
    raise RuntimeError("this task fails on purpose")


def copy(source, output, log):
    """copies source to output and logs the call (one line per call)"""
    with open(source) as f:
        text = f.read()
    with open(log, "a") as f:
        f.write(f"{output} {os.environ.get('SLURM_ARRAY_JOB_ID', '')}\\n")
    _publish(text, output)
'''


def _job_prefix() -> str:
    prefix = os.environ.get(JOB_PREFIX_VARIABLE) or DEFAULT_JOB_PREFIX
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", prefix):
        raise ValueError(f"{JOB_PREFIX_VARIABLE} must be letters, digits, '_', '.' or '-'")
    return prefix


def job_ids_of(text: str) -> list[int]:
    """the job ids of the sbatch answers in text"""
    return [int(x) for x in re.findall(r"Submitted batch job (\d+)", text)]


def read_field(paths: Sequence[str | os.PathLike[str]], key: str) -> list[str]:
    """the value of the key=value line in each file ('' if missing)"""
    values = []
    for p in paths:
        value = ""
        if os.path.isfile(p):
            m = re.search(rf"(?m)^{re.escape(key)}=([^\r\n]*)$", Path(p).read_text())
            if m:
                value = m.group(1).strip()
        values.append(value)
    return values


class Tee:
    """sys.stdout that also keeps what is written (the framework's sbatch answers), so that
    a long run still prints as it goes"""

    def __init__(self, stream) -> None:
        self.stream = stream
        self.parts: list[str] = []

    def write(self, text: str) -> int:
        self.parts.append(text)
        return self.stream.write(text)

    def flush(self) -> None:
        self.stream.flush()

    def __getattr__(self, name: str):
        return getattr(self.stream, name)

    def text(self) -> str:
        return "".join(self.parts)


def remove_folder(path: Path) -> None:
    """removes a test folder; network file systems can take a moment to release files that
    jobs just closed, so a few tries"""
    for _ in range(5):
        if not path.exists():
            return
        shutil.rmtree(path, ignore_errors=True)
        if path.exists():
            time.sleep(2)


class Watchdog:
    """cancels the test's jobs every WATCHDOG_PERIOD s once a run takes longer than its
    limit (jobs that hang, e.g. pending forever), so that the framework runs out of tries"""

    def __init__(self, cluster: SlurmCluster, minutes: float) -> None:
        self.cluster = cluster
        self.fired = False
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, args=(minutes * 60,), daemon=True)
        self._thread.start()

    def _run(self, delay: float) -> None:
        if self._stop.wait(delay):
            return
        while not self._stop.is_set():
            self.fired = True
            self.cluster.cancel_jobs()
            self._stop.wait(WATCHDOG_PERIOD)

    def stop(self) -> None:
        self._stop.set()
        self._thread.join(timeout=WATCHDOG_PERIOD + 70)


class SlurmCluster:
    """the real Slurm of one test: its config, job name, folders and helpers"""

    def __init__(self, config_file: str, scratch: Path) -> None:
        self.config_file = config_file
        self.conf = json.loads(Path(config_file).read_text())
        self.max_cpu_num = float(
            next((v for k, v in self.conf.items() if k.lower() == "maxcpunum"), 24)
        )
        self.scratch = scratch
        self.user = getpass.getuser()
        self.job_name = f"{_job_prefix()}{uuid.uuid4().hex[:8]}"
        self.slurm = Slurm()
        self.folders: list[Path] = []
        self.watchdogs: list[Watchdog] = []

    # jobs of the test

    def jobs(self, states: str = "") -> list[str] | None:
        """the ids (<job>_<array task>) of the test's jobs in the queue, None: no answer"""
        argv = ["squeue", "-h", "-r", "-u", self.user, f"--name={self.job_name}", "-o", "%i"]
        if states:
            argv += ["-t", states]
        r = self.slurm.run(argv)
        if not r.ok:
            return None
        return [x.strip() for x in r.stdout.splitlines() if x.strip()]

    def cancel_jobs(self) -> None:
        self.slurm.run(["scancel", "-u", self.user, f"--name={self.job_name}"], 60)

    def wait_until_gone(self, seconds: float = LEAVE_WAIT) -> list[str]:
        """waits up to seconds until no job of the test is in the queue; the jobs left"""
        deadline = time.monotonic() + seconds
        while True:
            left = self.jobs()
            if left == []:
                return []
            if time.monotonic() > deadline:
                return left if left is not None else ["(squeue gave no answer)"]
            time.sleep(2)

    def finish(self) -> None:
        for w in self.watchdogs:
            w.stop()
        self.cancel_jobs()
        self.wait_until_gone()
        for f in self.folders:
            remove_folder(f)

    # folders and tasks

    def shared_folder(self) -> Path:
        """a folder on storage the compute nodes reach, with job_logs/ and the Python tasks;
        removed after the test once its jobs are gone (they may still write job logs)"""
        folder = self.scratch / f"pk5d_slurm_test_{uuid.uuid4().hex[:8]}"
        (folder / "job_logs").mkdir(parents=True)
        self.folders.append(folder)
        self.task_module = f"j11_tasks_{folder.name[-8:]}"
        (folder / f"{self.task_module}.py").write_text(PYTHON_TASKS)
        return folder

    def tasks(
        self,
        folder: Path,
        n: int,
        language: str,
        sleep_seconds: Sequence[float] | None = None,
        fail: Sequence[bool] | None = None,
        final: str = "",
    ) -> tuple[list[str], list[str], list]:
        """n tasks that record the Slurm environment of their jobs in their outputs (bash
        scripts or Python tasks); task i fails without output if fail[i] (each attempt appends
        its CPUs to attempts_<i>.txt); task 1 also writes final"""
        sleep_seconds = list(sleep_seconds or [0] * n)
        fail = list(fail or [False] * n)
        inputs, outputs, tasks = [], [], []
        for i in range(1, n + 1):
            inp, out = folder / f"in_{i}.txt", folder / f"out_{i}.txt"
            inp.touch()
            attempts = folder / f"attempts_{i}.txt"
            if language == "bash":
                script = folder / f"task_{i}.sh"
                if fail[i - 1]:
                    body = BASH_FAIL.format(attempts=attempts)
                else:
                    body = BASH_TASK.format(sleep=int(sleep_seconds[i - 1]))
                script.write_text(f"#!/bin/bash\n{body}")
                task = f"bash {script} {out}"
                if final and i == 1:
                    task = f"{task} && cp {out} {final}"
            else:
                mod = self.task_module
                if fail[i - 1]:
                    task = Task(f"{mod}:fail", [str(attempts)], path=[str(folder)])
                else:
                    args = [str(out), sleep_seconds[i - 1], final if i == 1 else ""]
                    task = Task(f"{mod}:record", args, path=[str(folder)])
            inputs.append(str(inp))
            outputs.append(str(out))
            tasks.append(task)
        return inputs, outputs, tasks

    def run_framework(
        self, folder: Path, inputs, outputs, tasks, *, language: str, **params
    ) -> tuple[list[bool], bool]:
        """run_tasks in cluster mode with the test's settings; (done flags, whether the
        watchdog had to cancel jobs)"""
        minutes = params.pop("watchdog_minutes", WATCHDOG_MINUTES)
        settings = {
            "parse_cluster": True,
            "master_compute": False,
            "language": language,
            "job_log_dir": str(folder / "job_logs"),
            "config_file": self.config_file,
            "max_trial_num": 2,
            "job_time_limit": 0.25,
            "query_interval": 2,
            "slurm_constraint": f"-J {self.job_name}",
            **params,
        }
        dog = Watchdog(self, minutes)
        self.watchdogs.append(dog)
        try:
            done = run_tasks(inputs, outputs, tasks, **settings)
        finally:
            dog.stop()
        return done, dog.fired

    job_ids_of = staticmethod(job_ids_of)
    read_field = staticmethod(read_field)
    Tee = Tee

    @staticmethod
    def count_job_logs(log_dir: Path, text: str) -> int:
        """the job logs (job_*.out) in log_dir that mention text (each job echoes its
        command, which names the call folder, to its log)"""
        return sum(text in p.read_text(errors="replace") for p in log_dir.glob("job_*.out"))


def _config_and_scratch(tmp_path_factory: pytest.TempPathFactory) -> tuple[str, Path]:
    if os.name == "nt":
        pytest.skip("Slurm cluster mode is not supported on Windows.")
    config = os.environ.get("PETAKIT5D_TEST_SLURM_CONFIG", "")
    if not config:
        pytest.skip("PETAKIT5D_TEST_SLURM_CONFIG names no cluster config (the tests submit jobs).")
    assert os.path.isfile(config), f"Cluster configuration file {config} does not exist."
    r = Slurm().run(["sinfo"])
    assert r.ok, f"PETAKIT5D_TEST_SLURM_CONFIG is set but sinfo fails:\n{r.output}"
    scratch = os.environ.get("PETAKIT5D_TEST_SLURM_SCRATCH")
    if scratch:
        assert os.path.isdir(scratch), f"PETAKIT5D_TEST_SLURM_SCRATCH {scratch} is no folder"
        return config, Path(scratch)
    return config, tmp_path_factory.mktemp("slurm_scratch")


@pytest.fixture
def cluster(real_slurm, tmp_path_factory) -> Iterator[SlurmCluster]:
    """the real Slurm for one test; its jobs are cancelled and its folders removed after it"""
    config, scratch = _config_and_scratch(tmp_path_factory)
    c = SlurmCluster(config, scratch)
    try:
        yield c
    finally:
        c.finish()
