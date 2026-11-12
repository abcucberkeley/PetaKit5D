"""In-job slots: the tasks of a batch as one process each, several at once.

This replaces GNU parallel in MATLAB's compiled mode (``parallel --ungroup --retries 2 --jobs P
--delay 0.2 < input_<b>.txt``): up to P task processes run at once, starts are at least
SLOT_DELAY s apart, a task whose process fails runs once more (SLOT_TRIES runs in all), and the
processes write to the job's log directly. A task whose output (or the call's final output)
exists when its turn comes is skipped, as the task lines of the compiled mode do. In GPU jobs
slot k gets one GPU: the k-th device of the job's CUDA_VISIBLE_DEVICES if it lists enough,
else device k - 1 (MATLAB's ``CUDA_VISIBLE_DEVICES=$(({%} - 1))``). Each process also gets
GNU parallel's PARALLEL_JOBSLOT (the slot) and PARALLEL_SEQ (the task's place in the batch).

With runExtraTasks (call.json ``"slots": {"extra_tasks": true}``) a batch's tasks are followed
by tasks of later batches, in MATLAB's order, so that a job whose own tasks are done helps
with the others (each is skipped once its output exists).

Jobs run ``python -m petakit5d.jobs.slots <call folder> <batch> [--slots P] [--tries T]
[--gpu]`` (exit status: the number of tasks that failed, at most 101, as GNU parallel's); the
master runs ``run_batch`` in its process. A task process is ``python -m petakit5d.jobs.slots
--task <call folder> <index>``. This module imports only the standard library and the
framework's task records.
"""

from __future__ import annotations

import os
import queue
import shutil
import subprocess
import sys
import threading
import time
from collections import deque
from collections.abc import Sequence
from typing import Any

__all__ = [
    "MAX_EXIT_STATUS",
    "SLOT_DELAY",
    "SLOT_TRIES",
    "batch_tasks",
    "gpu_for_slot",
    "main",
    "run_batch",
    "run_slots",
    "run_task",
]

# s between the starts of two task processes (parallel --delay 0.2)
SLOT_DELAY = 0.2
# runs of a task whose process fails (parallel --retries 2: two runs in all)
SLOT_TRIES = 2
# the exit status of a batch is the number of failed tasks, at most this (GNU parallel's)
MAX_EXIT_STATUS = 101
# extra tasks of a batch after the first come from at most this many tasks' worth of batches
EXTRA_TASK_LIMIT = 500

USAGE = (
    "usage: python -m petakit5d.jobs.slots <call folder> <batch> [--slots P] [--tries T] [--gpu]\n"
    "       python -m petakit5d.jobs.slots --task <call folder> <index>"
)

# the pause and clock of the slot loop (tests replace them)
_sleep = time.sleep
_clock = time.monotonic


def gpu_for_slot(slot: int, n_slots: int, visible: str | None) -> str:
    """CUDA_VISIBLE_DEVICES of slot (1-based) of n_slots: the slot-th entry of the job's
    CUDA_VISIBLE_DEVICES if it lists n_slots or more, else slot - 1"""
    if visible:
        devices = [d.strip() for d in visible.split(",") if d.strip()]
        if len(devices) >= n_slots:
            return devices[slot - 1]
    return str(slot - 1)


def extra_batches(b: int, nb: int, outputs_existed: bool, batch_size: int) -> list[int]:
    """the batches whose tasks follow batch b's own with runExtraTasks (MATLAB's order): the
    later batches from the last down, interleaved b by b, and for b > 1 only the first
    max(2, ceil(500 / batch_size)) of them"""
    sets = list(range(nb - b + 1, b, -1))
    if outputs_existed:
        more = list(range(nb, max(b + 1, nb - b + 2) - 1, -1))
        seen: set[int] = set()
        sets = [x for x in sets + more if not (x in seen or seen.add(x))]
    if not sets:
        return []
    k = len(sets)
    rows = -(-k // b)
    order = [r * b + c for c in range(b) for r in range(rows)]
    sets = [sets[i] for i in order if i < k]
    if b > 1:
        sets = sets[: min(max(2, -(-EXTRA_TASK_LIMIT // batch_size)), len(sets))]
    return sets


def batch_tasks(meta: dict[str, Any], batch: int) -> list[int]:
    """0-based task indices a batch runs: its own, then (runExtraTasks) the extra ones, each
    extra batch's tasks last first"""
    n = int(meta["n_tasks"])
    bs = int(meta["task_batch_num"])
    nb = -(-n // bs)

    def own(x: int) -> list[int]:
        return list(range((x - 1) * bs, min(x * bs, n)))

    tasks = own(batch)
    opts = meta.get("slots") or {}
    if opts.get("extra_tasks"):
        for x in extra_batches(batch, nb, bool(opts.get("outputs_existed")), bs):
            tasks += own(x)[::-1]
    return tasks


def _output_exists(paths: Sequence[str]) -> bool:
    from petakit5d.jobs.files import batch_file_exist

    return any(batch_file_exist([p])[0] for p in paths if p)


def run_slots(
    call_folder: str,
    indices: Sequence[int],
    slots: int,
    tries: int = SLOT_TRIES,
    gpu: bool = False,
    delay: float = SLOT_DELAY,
    python: str | None = None,
) -> int:
    """Runs the tasks at indices (0-based) of a call folder, up to slots at once, each in its
    own process; returns the number of tasks that failed in every try."""
    from petakit5d.jobs.tasks import CallFolder

    cf = CallFolder(call_folder)
    lines = {x["index"] - 1: x for x in cf.lines(sorted(set(indices)))}
    final_out = cf.meta.get("final_out_fullpath") or ""
    bash = cf.language == "bash"
    launch = cf.meta.get("bash_launch_str") or ""
    python = python or sys.executable
    slots = max(1, int(slots))
    tries = max(1, int(tries))
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")

    todo: deque[tuple[int, int, int]] = deque((i, 1, seq) for seq, i in enumerate(indices, 1))
    free = list(range(1, slots + 1))
    running: dict[int, tuple[subprocess.Popen, int, int, int]] = {}
    ended: queue.Queue[int] = queue.Queue()
    failed = 0
    last_start = -float("inf")

    def watch(slot: int, proc: subprocess.Popen) -> None:
        proc.wait()
        ended.put(slot)

    sys.stdout.flush()
    sys.stderr.flush()
    while todo or running:
        while todo and free:
            i, attempt, seq = todo[0]
            out = lines[i]["output"]
            if _output_exists([out, final_out]):
                todo.popleft()
                print(f'output "{out}" already exists', flush=True)
                continue
            wait = last_start + delay - _clock()
            if wait > 0:
                _sleep(wait)
            todo.popleft()
            slot = free.pop(0)
            env = dict(os.environ, PARALLEL_JOBSLOT=str(slot), PARALLEL_SEQ=str(seq))
            if gpu:
                env["CUDA_VISIBLE_DEVICES"] = gpu_for_slot(slot, slots, visible)
            if bash:
                command = lines[i]["task"]["bash"]
                if launch:
                    command = f"{launch}; {command}"
                argv = [shutil.which("bash") or "bash", "-c", command]
            else:
                argv = [python, "-u", "-m", "petakit5d.jobs.slots", "--task", call_folder, str(i)]
            last_start = _clock()
            try:
                proc = subprocess.Popen(argv, env=env, stdin=subprocess.DEVNULL)
            except OSError as err:
                print(f"Task {i + 1} could not start: {err}", file=sys.stderr, flush=True)
                free.append(slot)
                free.sort()
                if attempt < tries:
                    todo.appendleft((i, attempt + 1, seq))
                else:
                    failed += 1
                continue
            running[slot] = (proc, i, attempt, seq)
            threading.Thread(target=watch, args=(slot, proc), daemon=True).start()
        if not running:
            continue
        slot = ended.get()
        proc, i, attempt, seq = running.pop(slot)
        free.append(slot)
        free.sort()
        if proc.returncode != 0:
            if attempt < tries:
                todo.appendleft((i, attempt + 1, seq))
            else:
                failed += 1
    return failed


def run_batch(
    call_folder: str,
    batch: int,
    slots: int,
    tries: int = SLOT_TRIES,
    gpu: bool = False,
    delay: float = SLOT_DELAY,
) -> int:
    """Runs batch (1-based) of a call folder in slots (with its extra tasks); returns the
    number of tasks that failed."""
    from petakit5d.jobs.tasks import CallFolder

    cf = CallFolder(call_folder)
    return run_slots(call_folder, batch_tasks(cf.meta, batch), slots, tries, gpu, delay)


def run_task(call_folder: str, index: int) -> int:
    """Runs one Python task (0-based index) of a call folder in this process; 0 if it
    returns, 1 if it raises (traceback on stderr)."""
    import traceback

    from petakit5d.jobs.tasks import CallFolder, run_record

    line = CallFolder(call_folder).lines([index])[0]
    try:
        run_record(line["task"])
    except Exception:
        sys.stdout.flush()
        print(f"Error in task {line['index']} ({line['task']['func']}):", file=sys.stderr)
        traceback.print_exc()
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else list(argv)
    try:
        if args and args[0] == "--task":
            if len(args) != 3:
                raise ValueError("--task takes a call folder and a task index")
            return run_task(args[1], int(args[2]))
        if len(args) < 2:
            raise ValueError("a call folder and a batch are needed")
        folder, batch = args[0], int(args[1])
        slots, tries, gpu = 1, SLOT_TRIES, False
        rest = args[2:]
        while rest:
            opt = rest.pop(0)
            if opt == "--gpu":
                gpu = True
            elif opt in ("--slots", "--tries") and rest:
                value = int(rest.pop(0))
                if opt == "--slots":
                    slots = value
                else:
                    tries = value
            else:
                raise ValueError(f"unknown option {opt!r}")
    except ValueError as err:
        print(f"{USAGE}\n{err}", file=sys.stderr)
        return 2
    try:
        failed = run_batch(folder, batch, slots, tries, gpu)
    except (OSError, KeyError) as err:
        print(f"petakit5d.jobs.slots: {folder}: {err}", file=sys.stderr)
        return 2
    return min(failed, MAX_EXIT_STATUS)


if __name__ == "__main__":
    sys.exit(main())
