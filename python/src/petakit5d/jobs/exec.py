"""The job entry point: ``python -m petakit5d.jobs.exec <call folder> <batch>``.

Runs the tasks of one batch (1-based) of a call folder (see ``jobs.tasks``) in this process,
as a MATLAB job runs its batch's function strings: Python tasks in order, stopping at the
first exception (its traceback goes to stderr; exit status 1); bash tasks as one shell
command, ``<BashLaunchStr>; <task 1>;<task 2>...``, with the shell's exit status. The last
line is MATLAB's ``toc``: ``Elapsed time is <s> seconds.`` The job script around it (J03)
writes the ``Start time:`` and ``End time:`` lines of the job log.

Before any task, the thread pools of the numeric libraries (OpenBLAS, OpenMP, MKL, numexpr)
are limited to the CPUs this process may use, unless their variables are set. With task
metrics (call settings ``task_metrics`` naming the run folder) each Python task, or the bash
batch, gets a record of mode ``slurm_python``. SIGTERM (scancel, the time limit) during a
Python task becomes the exception ``Terminated``, so that the task's ``finally`` blocks run;
the exit status is then 143. A bash batch ends on SIGTERM as before.

If call.json names a ``"runner"`` (``"package.module:function"``), the batch goes to it
instead, as ``runner(call_folder, batch)``, and its return value is the exit status (the
task runners of J06 plug in here).

Exit status 2: wrong arguments, or not a call folder.
"""

from __future__ import annotations

import os
import sys
import time
import traceback

USAGE = "usage: python -m petakit5d.jobs.exec <call folder> <batch>"

# thread pools of the numeric libraries, sized at their import (OpenBLAS starts one thread
# per CPU the process may use when NumPy is imported)
THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


class Terminated(BaseException):
    """SIGTERM (scancel, the job's time limit) as an exception, so that the running task's
    ``finally`` blocks and context managers run (a writer record is then marked ended at
    once). A BaseException: a task's ``except Exception`` doesn't catch it."""


def _on_sigterm(signum, frame) -> None:
    raise Terminated(f"signal {signum} (SIGTERM)")


def _default_sigterm() -> None:
    import signal

    try:
        if signal.getsignal(signal.SIGTERM) is _on_sigterm:
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
    except ValueError:  # (not the main thread)
        pass


def limit_library_threads() -> None:
    """Sets THREAD_VARIABLES to the CPUs this process may use (affinity and cgroup quota),
    except those already set; called before any task imports NumPy. Without it, importing
    NumPy with 24 CPUs costs 1.4 to 2.0 CPU-s of OpenBLAS threads, which packed task
    processes take from each other."""
    from petakit5d.core.threads import available_cpus

    cpus = str(available_cpus())
    for name in THREAD_VARIABLES:
        if not os.environ.get(name):
            os.environ[name] = cpus


def run_batch(call_folder: str, batch: int) -> int:
    """Runs batch ``batch`` of the call in ``call_folder``; returns the exit status."""
    limit_library_threads()
    from petakit5d.jobs.tasks import CallFolder, resolve, run_record

    cf = CallFolder(call_folder)
    runner = cf.meta.get("runner")
    if runner:
        return int(resolve({"func": runner})(call_folder, batch) or 0)
    indices = cf.batch(batch)
    if batch < 1 or not indices:
        print(f"petakit5d.jobs.exec: batch {batch} is not in 1..{_batch_num(cf)}", file=sys.stderr)
        return 2
    lines = cf.lines(indices)
    metrics = _metrics(cf)
    t0 = time.perf_counter()
    status = 0
    if cf.language == "bash":
        from petakit5d.jobs.framework import run_shell

        command = ";".join(x["task"]["bash"] for x in lines)
        launch = cf.meta.get("bash_launch_str") or ""
        if launch:
            command = f"{launch}; {command}"
        tm = metrics.task_start(cf, list(indices), "slurm_python") if metrics else None
        _default_sigterm()  # (the shell and this process end on SIGTERM, as before)
        status = run_shell(command)
        if metrics:
            metrics.task_stop(tm, status)
    else:
        for i, x in zip(indices, lines, strict=True):
            tm = metrics.task_start(cf, [i], "slurm_python") if metrics else None
            try:
                run_record(x["task"])
            except Terminated as err:
                if metrics:
                    metrics.task_stop(tm, err)
                raise
            except Exception as err:
                if metrics:
                    metrics.task_stop(tm, err)
                sys.stdout.flush()
                print(f"Error in task {x['index']} ({x['task']['func']}):", file=sys.stderr)
                traceback.print_exc()
                status = 1
                break
            if metrics:
                metrics.task_stop(tm, None)
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)
    return status


def _metrics(cf):
    # the metrics module if the call records task metrics (its settings name the run folder)
    folder = cf.meta.get("settings", {}).get("task_metrics")
    if not isinstance(folder, str) or not folder:
        return None
    try:
        from petakit5d.jobs import metrics
    except ImportError:
        print("Warning: task metrics are not available (petakit5d.jobs.metrics)", file=sys.stderr)
        return None
    return metrics


def _batch_num(cf) -> int:
    bs = int(cf.meta["task_batch_num"])
    return -(-cf.n // bs)


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print(USAGE, file=sys.stderr)
        return 2
    try:
        batch = int(args[1])
    except ValueError:
        print(f"{USAGE}\nbatch must be an integer, not {args[1]!r}", file=sys.stderr)
        return 2
    import signal

    # (only in the main thread; a bash batch's shell gets Slurm's signal itself)
    previous = signal.signal(signal.SIGTERM, _on_sigterm)
    try:
        return run_batch(args[0], batch)
    except (OSError, ValueError, KeyError) as err:
        print(f"petakit5d.jobs.exec: {args[0]}: {err}", file=sys.stderr)
        return 2
    except Terminated as err:
        sys.stdout.flush()
        print(f"petakit5d.jobs.exec: terminated by {err}", file=sys.stderr, flush=True)
        return 128 + signal.SIGTERM
    finally:
        signal.signal(signal.SIGTERM, previous)


if __name__ == "__main__":
    sys.exit(main())
