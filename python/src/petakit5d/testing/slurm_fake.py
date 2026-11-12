"""A stand-in Slurm for tests: the shell stubs of the MATLAB suite's FakeSlurmFixture, and the
no-Slurm and real-Slurm fixtures.

``FakeSlurm`` writes executables first on the PATH: sinfo succeeds, sbatch prints a new job id
(1000, 1001, ...), squeue lists the jobs whose states a test set with ``set_job_states``, and
scancel does nothing. Every call is logged (``<folder>/<command>.log``, one line of arguments
per call, and ``history.log``, all calls in order), so tests can count them. Not available on
Windows (the stubs need a POSIX shell).

A job's states are used one per squeue call that names it (-j), and the last one stays:
``set_job_states('123_7', ['PD', 'PD', 'R'])`` reports the job pending twice and then running.
The state 'done' means the job is gone. A pending job's reason (squeue -o %r) is Resources, or
the one after a colon ('PD:None': not looked at by the scheduler yet). squeue without -j lists
all jobs in their current state (it uses no state), with -t and -o (%i, %T, %t, %r) and -h. The
jobs run nothing, but a job can write (empty) output files when it is first reported running,
and ``reject_submissions`` makes sbatch fail.

For jobs on several partitions: sbatch records each job's partition, QOS, CPUs, memory and GPUs,
which squeue --me -O (JobArrayID, Partition, QOS, NumCPUs, tres-alloc, Reason; with -t and -p)
lists, and scontrol update changes while the job is pending (else it fails, as for a job that
started); scontrol show node gives the nodes a test set with ``set_nodes`` (none by default),
sbatch --test-only answers per partition as ``set_test_only`` says (start now by default; it
takes no job id), and other scontrol commands run the real scontrol if there is one.

The stubs are those of the MATLAB fixture, with its default partition and node name as
parameters (made-up names by default), and two additions for failure tests: ``hang`` (a command
that doesn't answer) and ``fail`` (a command that exits 1 with a message), both logged first.
"""

# ruff: noqa: E501  (the shell lines of the stubs are the MATLAB fixture's)
from __future__ import annotations

import os
import re
import shutil
from collections.abc import Iterator, Sequence
from pathlib import Path

import pytest

__all__ = [
    "DEFAULT_NODE",
    "DEFAULT_PARTITION",
    "SLURM_JOB_VARIABLES",
    "FakeSlurm",
    "fake_slurm",
    "hide_slurm",
    "no_slurm",
    "no_slurm_session",
    "real_slurm",
    "replay",
    "replay_history",
    "stub_text",
]

DEFAULT_PARTITION = "cpu_queue"
DEFAULT_NODE = "node12"

# the variables of a job that PetaKit5D and its scripts read (NoSlurmFixture)
SLURM_JOB_VARIABLES = (
    "SLURM_JOB_ID",
    "SLURM_JOBID",
    "SLURM_ARRAY_JOB_ID",
    "SLURM_ARRAY_TASK_ID",
    "SLURMD_NODENAME",
    "SLURM_MEM_PER_NODE",
    "SLURM_MEM_PER_CPU",
    "SLURM_CPUS_ON_NODE",
    "SLURM_JOB_START_TIME",
    "SLURM_JOB_END_TIME",
)
STUBS_VARIABLE = "PETAKIT5D_TEST_SLURM_STUBS"

_STUBS: dict[str, tuple[str, ...]] = {
    "sinfo": ('echo "PARTITION AVAIL TIMELIMIT NODES STATE NODELIST"',),
    "sbatch": (
        'tomb() { v="$1"; case "$v" in *G) echo $(( ${v%G} * 1024 )) ;; *M) echo "${v%M}" ;; *T) echo $(( ${v%T} * 1048576 )) ;; *) echo "$v" ;; esac; }',
        'part=@PARTITION@; qos=normal; cpus=1; mpc=""; mem=""; g=""; arr=""; test_only=0; prev=""',
        'for a in "$@"; do',
        '  case "$prev" in -p|--partition) part="$a" ;; -q|--qos) qos="$a" ;; --mem) mem="$a" ;; --mem-per-cpu) mpc="$a" ;; --gres) g="$a" ;; esac',
        '  case "$a" in',
        "    --test-only) test_only=1 ;;",
        '    --partition=*) part="${a#--partition=}" ;;',
        '    --qos=*) qos="${a#--qos=}" ;;',
        '    -p?*) part="${a#-p}" ;;',
        '    --cpus-per-task=*) cpus="${a#--cpus-per-task=}" ;;',
        '    --mem-per-cpu=*) mpc="${a#--mem-per-cpu=}" ;;',
        '    --mem=*) mem="${a#--mem=}" ;;',
        '    --gres=*) g="${a#--gres=}" ;;',
        '    --array=*) arr="${a#--array=}" ;;',
        "  esac",
        '  prev="$a"',
        "done",
        'if [ -n "$mpc" ]; then mem=$(( $(tomb "$mpc") * cpus )); elif [ -n "$mem" ]; then mem=$(tomb "$mem"); else mem=$(( 1000 * cpus )); fi',
        "gpus=0",
        'case "$g" in gpu*) gpus="${g##*:}"; case "$gpus" in ""|*[!0-9]*) gpus=1 ;; esac ;; esac',
        'if [ "$test_only" = 1 ]; then',
        '  ans=$(cat "$D/testonly/$part" 2>/dev/null || echo now)',
        '  case "$ans" in',
        '    error) echo "allocation failure: Invalid qos specification" >&2; exit 1 ;;',
        "    later) t=2099-01-01T00:00:00 ;;",
        "    *) t=$(date +%Y-%m-%dT%H:%M:%S) ;;",
        "  esac",
        '  echo "sbatch: Job 99999 to start at $t a using $cpus processors on nodes fake1 in partition $part" >&2',
        "  exit 0",
        "fi",
        'r=$(cat "$D/reject" 2>/dev/null || echo 0)',
        'if [ "$r" -gt 0 ]; then',
        '  echo $((r - 1)) > "$D/reject"',
        '  echo "sbatch: error: Batch job submission failed: Socket timed out on send/recv operation" >&2',
        "  exit 1",
        "fi",
        'n=$(cat "$D/next_id" 2>/dev/null || echo 1000)',
        'echo $((n + 1)) > "$D/next_id"',
        'echo "partition=$part qos=$qos cpus=$cpus mem=$mem gpus=$gpus" > "$D/jobs/${n}_${arr}"',
        'if [ -f "$D/make_logs" ]; then',
        '  out=""; err=""; arr=""',
        "  while [ $# -gt 0 ]; do",
        '    case "$1" in',
        '      -o) out="$2"; shift ;;',
        '      -e) err="$2"; shift ;;',
        '      --array=*) arr="${1#--array=}" ;;',
        "    esac",
        "    shift",
        "  done",
        '  T="$(dirname "$D")"',
        '  P="$(cd "$T" && pwd -P)"',
        '  for f in "$out" "$err"; do',
        '    g=$(printf \'%s\' "$f" | sed "s/%A/$n/g; s/%a/$arr/g; s/%j/$n/g")',
        '    case "$g" in "$T"/*|"$P"/*) : > "$g" 2> /dev/null ;; esac',
        "  done",
        "fi",
        'echo "Submitted batch job $n"',
    ),
    "scancel": (),
    "scontrol": (
        'if [ "$1" = show ] && [ "$2" = node ]; then',
        '  i=$(cat "$D/nodes_next" 2>/dev/null || echo 1)',
        '  [ -f "$D/nodes_$i" ] || exit 0',
        '  cat "$D/nodes_$i"',
        '  [ -f "$D/nodes_$((i + 1))" ] && echo $((i + 1)) > "$D/nodes_next"',
        "  exit 0",
        "fi",
        'if [ "$1" = update ]; then',
        '  id=""',
        '  for a in "$@"; do case "$a" in JobId=*) id="${a#JobId=}" ;; esac; done',
        '  f="$D/state/$id"',
        '  if [ ! -f "$f" ]; then echo "Invalid job id specified" >&2; exit 1; fi',
        '  case "$(head -n 1 "$f")" in PD*) ;; *) echo "$id: Job is no longer pending execution" >&2; exit 1 ;; esac',
        '  if [ -f "$D/update_fail" ]; then echo "$id: Invalid qos specification" >&2; exit 1; fi',
        '  part=@PARTITION@; qos=normal; cpus=1; mem=1000; gpus=0; mpc=""',
        '  for kv in $(cat "$D/jobs/$id" 2>/dev/null); do',
        '    case "$kv" in partition=*) part="${kv#*=}" ;; qos=*) qos="${kv#*=}" ;; cpus=*) cpus="${kv#*=}" ;; mem=*) mem="${kv#*=}" ;; gpus=*) gpus="${kv#*=}" ;; esac',
        "  done",
        '  for a in "$@"; do',
        '    case "$a" in Partition=*) part="${a#*=}" ;; QOS=*) qos="${a#*=}" ;; NumCPUs=*) cpus="${a#*=}" ;; MinMemoryCPU=*) mpc="${a#*=}" ;; MinMemoryNode=*) mem="${a#*=}" ;; esac',
        "  done",
        '  [ -n "$mpc" ] && mem=$(( mpc * cpus ))',
        '  echo "partition=$part qos=$qos cpus=$cpus mem=$mem gpus=$gpus" > "$D/jobs/$id"',
        "  exit 0",
        "fi",
        'for d in $(echo "$PATH" | tr ":" " "); do',
        '  [ "$d" = "$D" ] && continue',
        '  [ -x "$d/scontrol" ] && exec "$d/scontrol" "$@"',
        "done",
        'echo "scontrol: command not found" >&2',
        "exit 127",
    ),
    "squeue": (
        'ids=""; states=""; fmt=""; hdr=1; peek=0; long=""; parts=""',
        "while [ $# -gt 0 ]; do",
        '  case "$1" in',
        '    -j) ids="$2"; shift ;;',
        '    -t) states="$2"; shift ;;',
        '    -o) fmt="$2"; shift ;;',
        '    -O) long="$2"; shift ;;',
        '    -p) parts="$2"; shift ;;',
        "    -u) shift ;;",
        "    -h) hdr=0 ;;",
        "  esac",
        "  shift",
        "done",
        'if [ -n "$long" ]; then',
        '  for id in $(ls "$D/state" | grep -v "\\\\."); do',
        '    st=$(head -n 1 "$D/state/$id")',
        "    reason=Resources",
        '    case "$st" in *:*) reason=${st#*:}; st=${st%%:*} ;; esac',
        '    [ "$st" = "done" ] && continue',
        '    [ "$st" = "R" ] && reason=None',
        '    [ -n "$states" ] && [ "$states" != "$st" ] && continue',
        "    part=@PARTITION@; qos=normal; cpus=1; mem=1000; gpus=0",
        '    for kv in $(cat "$D/jobs/$id" 2>/dev/null); do',
        '      case "$kv" in partition=*) part="${kv#*=}" ;; qos=*) qos="${kv#*=}" ;; cpus=*) cpus="${kv#*=}" ;; mem=*) mem="${kv#*=}" ;; gpus=*) gpus="${kv#*=}" ;; esac',
        "    done",
        '    case ",$parts," in ,,|*",$part,"*) ;; *) continue ;; esac',
        '    jid="$id"',
        '    case "$id" in *_*) [ "$st" = PD ] && jid="${id%%_*}_[${id#*_}]" ;; esac',
        '    tres="cpu=$cpus,mem=${mem}M,node=1,billing=$cpus"',
        '    [ "$gpus" -gt 0 ] && tres="$tres,gres/gpu=$gpus"',
        '    line=""',
        '    for fld in $(echo "$long" | tr "," " "); do',
        '      case "${fld%%:*}" in',
        '        JobArrayID) v="$jid" ;; Partition) v="$part" ;; QOS) v="$qos" ;; NumCPUs) v="$cpus" ;;',
        '        tres-alloc) v="$tres" ;; Reason) v="$reason" ;; *) v="" ;;',
        "      esac",
        '      line="$line$v|"',
        "    done",
        '    echo "$line"',
        "  done",
        "  exit 0",
        "fi",
        'if [ -z "$ids" ]; then',
        "  peek=1",
        '  ids=$(ls "$D/state" | grep -v "\\." | tr "\\n" ",")',
        "fi",
        '[ "$hdr" = 1 ] && [ -z "$fmt" ] && echo "             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)"',
        'for id in $(echo "$ids" | tr "," " "); do',
        '  f="$D/state/$id"',
        '  [ -f "$f" ] || continue',
        '  st=$(head -n 1 "$f")',
        '  if [ "$peek" = 0 ]; then',
        '    if [ $(wc -l < "$f") -gt 1 ]; then tail -n +2 "$f" > "$f.next" && mv "$f.next" "$f"; fi',
        '    if [ "${st%%:*}" = "R" ] && [ -f "$f.outputs" ]; then',
        '      while IFS= read -r o; do : > "$o"; done < "$f.outputs"',
        '      rm "$f.outputs"',
        "    fi",
        "  fi",
        "  reason=Resources",
        '  case "$st" in *:*) reason=${st#*:}; st=${st%%:*} ;; esac',
        '  [ "$st" = "done" ] && continue',
        "  long=PENDING",
        '  [ "$st" = "R" ] && long=RUNNING && reason=None',
        '  [ -n "$states" ] && [ "$states" != "$long" ] && [ "$states" != "$st" ] && continue',
        '  if [ -n "$fmt" ]; then',
        '    echo "$fmt" | sed -e "s/%i/$id/g" -e "s/%T/$long/g" -e "s/%t/$st/g" -e "s/%r/$reason/g"',
        "  else",
        '    echo "             $id       @PARTITION@     wrap     user  $st       0:01      1 @NODE@"',
        "  fi",
        "done",
    ),
}

_HEADER = (
    "#!/bin/sh",
    'D="$(dirname "$0")"',
    'echo "$*" >> "$D/@CMD@.log"',
    'echo "@CMD@ $*" >> "$D/history.log"',
)
# (Python additions: a command that doesn't answer, or fails)
_FAILURES = (
    '[ -f "$D/@CMD@.hang" ] && exec sleep "$(cat "$D/@CMD@.hang")"',
    'if [ -f "$D/@CMD@.fail" ]; then cat "$D/@CMD@.fail" >&2; exit 1; fi',
)


def stub_text(cmd: str, partition: str = DEFAULT_PARTITION, node: str = DEFAULT_NODE) -> str:
    """the script of one fake command"""
    lines = [*_HEADER, *_FAILURES, *_STUBS[cmd]]
    text = "\n".join(lines) + "\n"
    return text.replace("@CMD@", cmd).replace("@PARTITION@", partition).replace("@NODE@", node)


def _write_executable(path: Path, text: str) -> None:
    path.write_text(text, newline="\n")
    path.chmod(path.stat().st_mode | 0o111)


def _number(text: str) -> int | float | str:
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


class FakeSlurm:
    """the fake Slurm commands in folder (see the module docstring)"""

    def __init__(
        self,
        folder: str | os.PathLike[str],
        partition: str = DEFAULT_PARTITION,
        node: str = DEFAULT_NODE,
    ):
        self.folder = Path(folder)
        self.partition = partition
        self.node = node

    def install(self) -> FakeSlurm:
        if os.name == "nt":
            raise RuntimeError("The fake Slurm needs a POSIX shell.")
        for sub in ("state", "jobs", "testonly"):
            (self.folder / sub).mkdir(parents=True, exist_ok=True)
        for cmd in _STUBS:
            _write_executable(self.folder / cmd, stub_text(cmd, self.partition, self.node))
        return self

    def path_with(self, path: str | None = None) -> str:
        """PATH with the fake commands first"""
        path = os.environ.get("PATH", "") if path is None else path
        return f"{self.folder}{os.pathsep}{path}" if path else str(self.folder)

    # the job states and what the commands answer

    def set_job_states(
        self, job_id: str, states: Sequence[str], outputs: Sequence[str] | None = None
    ) -> None:
        """job_id as squeue shows it ('123_7'); outputs: files the job writes (empty) when
        squeue first reports it running"""
        (self.folder / "state" / job_id).write_text("".join(f"{s}\n" for s in states))
        if outputs is not None:
            (self.folder / "state" / f"{job_id}.outputs").write_text(
                "".join(f"{os.fspath(o)}\n" for o in outputs)
            )

    def clear_calls(self) -> None:
        """forget the calls so far (their logs) and number the next jobs from 1000 again"""
        for log in self.folder.glob("*.log"):
            log.unlink()
        (self.folder / "next_id").unlink(missing_ok=True)

    def reject_submissions(self, n: int) -> None:
        """the next n sbatch calls fail without submitting a job (a busy controller)"""
        (self.folder / "reject").write_text(f"{int(n)}\n")

    def set_nodes(self, *node_sets: Sequence[str]) -> None:
        """the lines of scontrol show node -o (see node_line); with several sets each call
        reports the next one, and the last one stays"""
        for old in self.folder.glob("nodes_*"):
            old.unlink()
        for i, lines in enumerate(node_sets, 1):
            (self.folder / f"nodes_{i}").write_text("".join(f"{x}\n" for x in lines))

    def set_test_only(self, partition: str, answer: str) -> None:
        """what sbatch --test-only answers for partition: 'now', 'later' (2099) or 'error'"""
        (self.folder / "testonly" / partition).write_text(f"{answer}\n")

    def set_job_attributes(
        self,
        job_id: str,
        partition: str,
        cpus: int,
        mem_mb: int,
        gpus: int = 0,
        qos: str = "normal",
    ) -> None:
        """a job's resources as squeue -O lists them (sbatch records its own jobs')"""
        (self.folder / "jobs" / job_id).write_text(
            f"partition={partition} qos={qos} cpus={int(cpus)} mem={int(mem_mb)} gpus={int(gpus)}\n"
        )

    def job_attributes(self, job_id: str) -> dict[str, int | float | str]:
        """a job's partition, qos, cpus, mem (MB) and gpus as recorded ({} if none)"""
        fn = self.folder / "jobs" / job_id
        if not fn.is_file():
            return {}
        out: dict[str, int | float | str] = {}
        for kv in fn.read_text().split():
            k, _, v = kv.partition("=")
            out[k] = _number(v)
        return out

    def fail_updates(self) -> None:
        """from now on scontrol update fails (as for a QOS the new partition doesn't allow)"""
        (self.folder / "update_fail").touch()

    def make_log_files(self) -> None:
        """from now on sbatch makes the files of its -o and -e options (%A, %a, %j filled in),
        as Slurm does when a job starts; only files under the fake's parent folder"""
        (self.folder / "make_logs").touch()

    def hang(self, cmd: str, seconds: float = 3600) -> None:
        """from now on cmd doesn't answer for seconds (Python addition)"""
        (self.folder / f"{cmd}.hang").write_text(f"{seconds:g}\n")

    def fail(
        self, cmd: str, message: str = "slurm_load_jobs error: Unable to contact slurm controller"
    ) -> None:
        """from now on cmd exits 1 after printing message to stderr (Python addition)"""
        (self.folder / f"{cmd}.fail").write_text(f"{message}\n")

    def recover(self, cmd: str) -> None:
        """undo hang and fail for cmd"""
        (self.folder / f"{cmd}.hang").unlink(missing_ok=True)
        (self.folder / f"{cmd}.fail").unlink(missing_ok=True)

    # the calls

    def _lines(self, name: str) -> list[str]:
        fn = self.folder / name
        if not fn.is_file():
            return []
        text = fn.read_text()
        return text.split("\n")[:-1] if text.endswith("\n") else text.split("\n")

    def calls(self, cmd: str) -> list[str]:
        """the arguments of every call of cmd, one string per call"""
        return self._lines(f"{cmd}.log")

    def call_count(self, cmd: str) -> int:
        return len(self.calls(cmd))

    def history(self) -> list[str]:
        """every call of all commands in call order, as 'command arguments'"""
        return self._lines("history.log")

    @staticmethod
    def node_line(
        name: str,
        partitions: str,
        cpus: int,
        cpu_alloc: int,
        mem_mb: int,
        alloc_mb: int,
        gpus: int = 0,
        gpu_alloc: int = 0,
        state: str = "MIXED",
        features: str = "(null)",
    ) -> str:
        """a node as scontrol show node -o shows it"""
        cfg = f"cpu={cpus},mem={mem_mb}M,billing={cpus}"
        alloc = f"cpu={cpu_alloc},mem={alloc_mb}M"
        gres = "(null)"
        if gpus > 0:
            cfg = f"{cfg},gres/gpu={gpus}"
            gres = f"gpu:{gpus}"
        if gpu_alloc > 0:
            alloc = f"{alloc},gres/gpu={gpu_alloc}"
        return (
            f"NodeName={name} Arch=x86_64 CoresPerSocket={cpus} CPUAlloc={cpu_alloc} CPUEfctv={cpus} "
            f"CPUTot={cpus} AvailableFeatures={features} ActiveFeatures={features} Gres={gres} OS=Linux "
            f"4.18.0 #1 SMP RealMemory={mem_mb} AllocMem={alloc_mb} State={state} Partitions={partitions} "
            f"CfgTRES={cfg} AllocTRES={alloc}"
        )


# --- replaying recordings of the MATLAB fixture -------------------------------------------------
#
# A recording is a list of events in call order: {"op": <MATLAB setter>, "args": [...]} (a call
# of the fixture's setters) or {"cmd": <command>, "argv": [...], "stdout", "stderr", "rc"} (one
# command and its answer). Paths in them are placeholders: <F> the fake's folder, <B> its parent,
# <T1>, <T2>, ... other folders in <B> (the tests' folders); times of day are <NOW>.

_OPS = {
    "setJobStates": lambda f, a: f.set_job_states(a[0], _as_list(a[1]), _as_list(a[2]) or None),
    "clearCalls": lambda f, a: f.clear_calls(),
    "rejectSubmissions": lambda f, a: f.reject_submissions(a[0]),
    "setNodes": lambda f, a: f.set_nodes(*[_as_list(x) for x in a]),
    "setTestOnly": lambda f, a: f.set_test_only(a[0], a[1]),
    "setJobAttributes": lambda f, a: f.set_job_attributes(a[0], a[1], a[2], a[3], a[4], a[5]),
    "failUpdates": lambda f, a: f.fail_updates(),
    "makeLogFiles": lambda f, a: f.make_log_files(),
}
_NOW = re.compile(r"\b(?!2099-)\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\b")
_PATH_ARG = re.compile(r"^<(?:T\d+|B)>/[^\s<>]*$")


def _as_list(x: object) -> list:
    if x is None:
        return []
    return list(x) if isinstance(x, list) else [x]


def replay(
    events: Sequence[dict], root: str | os.PathLike[str], env: dict[str, str] | None = None
) -> list[str]:
    """run a recording against a new FakeSlurm in root (empty); the events whose answers
    differ, as text (empty: identical answers)"""
    import subprocess

    root = Path(root)
    fake = FakeSlurm(root / "fake").install()
    places = {"<F>": str(fake.folder), "<B>": str(root)}

    def path_of(tag: str) -> str:
        if tag not in places:
            d = root / tag.strip("<>").lower()
            d.mkdir(exist_ok=True)
            places[tag] = str(d)
        return places[tag]

    def fill(x):
        if isinstance(x, str):
            out = re.sub(r"<(?:F|B|T\d+)>", lambda m: path_of(m.group(0)), x)
            # (the folders of the files the stubs write were the test's: made here)
            if _PATH_ARG.match(x) and "/" in x:
                os.makedirs(os.path.dirname(out), exist_ok=True)
            return out
        if isinstance(x, list):
            return [fill(v) for v in x]
        return x

    def unfill(text: str) -> str:
        for tag, p in sorted(places.items(), key=lambda kv: -len(kv[1])):
            text = text.replace(p, tag)
        return _NOW.sub("<NOW>", text)

    run_env = dict(os.environ if env is None else env)
    run_env["PATH"] = fake.path_with(run_env.get("PATH", ""))
    differ = []
    for i, e in enumerate(events):
        if "op" in e:
            _OPS[e["op"]](fake, fill(e["args"]))
            continue
        r = subprocess.run(
            [e["cmd"], *fill(e["argv"])], capture_output=True, text=True, env=run_env
        )
        got = (unfill(r.stdout), unfill(r.stderr), r.returncode)
        want = (e["stdout"], e["stderr"], e["rc"])
        if got != want:
            call = f"{e['cmd']} {' '.join(e['argv'])}"[:200]
            differ.append(f"event {i} ({call}): got {got!r}, recorded {want!r}")
    return differ


def replay_history(root: str | os.PathLike[str]) -> list[str]:
    """history.log of a replay, with its paths as placeholders (compare with the recording's)"""
    root = Path(root)
    fake = FakeSlurm(root / "fake")
    places = {"<F>": str(fake.folder), "<B>": str(root)}
    places.update(
        {f"<{d.name.upper()}>": str(d) for d in root.iterdir() if re.fullmatch(r"t\d+", d.name)}
    )
    out = []
    for line in fake.history():
        for tag, p in sorted(places.items(), key=lambda kv: -len(kv[1])):
            line = line.replace(p, tag)
        out.append(_NOW.sub("<NOW>", line))
    return out


def _install_failing_stubs(folder: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for cmd in ("sinfo", "sbatch", "squeue", "scancel"):
        _write_executable(folder / cmd, "#!/bin/sh\nexit 1\n")


def hide_slurm(mp: pytest.MonkeyPatch, folder: Path) -> None:
    """NoSlurmFixture: failing sinfo, sbatch, squeue and scancel first on the PATH (in folder,
    named by PETAKIT5D_TEST_SLURM_STUBS), and the variables of the tests' own Slurm job unset,
    so that code under test runs locally and as outside a job. Nothing changes on Windows."""
    if os.name == "nt":
        return
    _install_failing_stubs(folder)
    mp.setenv("PATH", f"{folder}{os.pathsep}{os.environ.get('PATH', '')}")
    mp.setenv(STUBS_VARIABLE, str(folder))
    for name in SLURM_JOB_VARIABLES:
        mp.delenv(name, raising=False)


@pytest.fixture
def fake_slurm(
    tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
) -> Iterator[FakeSlurm]:
    """a FakeSlurm first on the PATH for one test (its folder in pytest's base temp folder, so
    that log files of make_log_files may go to any test's tmp_path)"""
    if os.name == "nt":
        pytest.skip("The fake Slurm needs a POSIX shell.")
    fake = FakeSlurm(tmp_path_factory.mktemp("fake_slurm")).install()
    monkeypatch.setenv("PATH", fake.path_with())
    yield fake
    shutil.rmtree(fake.folder, ignore_errors=True)


@pytest.fixture
def no_slurm(
    tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
) -> Iterator[None]:
    """Slurm hidden from one test (see hide_slurm)"""
    hide_slurm(monkeypatch, tmp_path_factory.mktemp("no_slurm"))
    yield


@pytest.fixture(scope="session")
def no_slurm_session(tmp_path_factory: pytest.TempPathFactory) -> Iterator[None]:
    """Slurm hidden from all tests of a session (as the MATLAB suite's shared NoSlurmFixture);
    meant for an autouse fixture in the suite's conftest"""
    with pytest.MonkeyPatch.context() as mp:
        hide_slurm(mp, tmp_path_factory.mktemp("no_slurm"))
        yield


@pytest.fixture
def real_slurm(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """RealSlurmFixture: the real Slurm commands visible again (the no-Slurm stubs off PATH)"""
    stubs = os.environ.get(STUBS_VARIABLE)
    if stubs:
        entries = os.environ.get("PATH", "").split(os.pathsep)
        monkeypatch.setenv("PATH", os.pathsep.join(e for e in entries if e != stubs))
    yield
