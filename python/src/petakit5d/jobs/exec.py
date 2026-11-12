"""The job entry point: ``python -m petakit5d.jobs.exec <call folder> <batch>``.

Runs the tasks of one batch (1-based) of a call folder (see ``jobs.tasks``) in this process,
as a MATLAB job runs its batch's function strings: Python tasks in order, stopping at the
first exception (its traceback goes to stderr; exit status 1); bash tasks as one shell
command, ``<BashLaunchStr>; <task 1>;<task 2>...``, with the shell's exit status. The last
line is MATLAB's ``toc``: ``Elapsed time is <s> seconds.`` The job script around it (J03)
writes the ``Start time:`` and ``End time:`` lines of the job log.

If call.json names a ``"runner"`` (``"package.module:function"``), the batch goes to it
instead, as ``runner(call_folder, batch)``, and its return value is the exit status (the
task runners of J06 plug in here).

Exit status 2: wrong arguments, or not a call folder.
"""

from __future__ import annotations

import sys
import time
import traceback

USAGE = "usage: python -m petakit5d.jobs.exec <call folder> <batch>"


def run_batch(call_folder: str, batch: int) -> int:
    """Runs batch ``batch`` of the call in ``call_folder``; returns the exit status."""
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
    t0 = time.perf_counter()
    status = 0
    if cf.language == "bash":
        from petakit5d.jobs.framework import run_shell

        command = ";".join(x["task"]["bash"] for x in lines)
        launch = cf.meta.get("bash_launch_str") or ""
        if launch:
            command = f"{launch}; {command}"
        status = run_shell(command)
    else:
        for x in lines:
            try:
                run_record(x["task"])
            except Exception:
                sys.stdout.flush()
                print(f"Error in task {x['index']} ({x['task']['func']}):", file=sys.stderr)
                traceback.print_exc()
                status = 1
                break
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)
    return status


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
    try:
        return run_batch(args[0], batch)
    except (OSError, ValueError, KeyError) as err:
        print(f"petakit5d.jobs.exec: {args[0]}: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
