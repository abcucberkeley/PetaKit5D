"""One framework task as a process, for the start-up and teardown table of the three
implementations (``task_costs.json``): the work of a task (a tile's deskew/rotation, or a tiny
crop), run the way a job runs it.

    python task_probe.py dsr TILE        deskew_rotate_file of TILE (its DSR/ result removed
                                         first), with the demo's parameters
    python task_probe.py tiny IN OUT     crop_file of IN to a 32^3 box in OUT
    python task_probe.py call FOLDER KIND ARGS...
                                         writes a call folder of one such task, for
                                         ``python -m petakit5d.jobs.exec FOLDER 1``

Run by the harness, the script marks ``ready`` after its imports (the first line of the work);
the harness marks the first line of the process and times launch and exit.
"""

from __future__ import annotations

import os
import shutil
import sys

from petakit5d.testing.bench import mark

# the deskew/rotation of the demos' tiles (demo_geometric_transformation)
DSR_ARGS = (0.108, 0.3)
DSR_KWARGS = {"skew_angle": 32.45, "reverse": True, "rotate": True, "save_16bit": True}
TINY_BOX = (0, 0, 0, 32, 32, 32)


def _clear_dsr(tile: str) -> None:
    shutil.rmtree(os.path.join(os.path.dirname(os.path.abspath(tile)), "DSR"), ignore_errors=True)


def main(argv: list[str]) -> None:
    kind = argv[0]
    if kind == "call":
        from petakit5d.jobs.tasks import Call, Task, task_record, write_call

        folder, what, rest = argv[1], argv[2], argv[3:]
        if what == "dsr":
            task = Task(
                "petakit5d.deskew_rotate.file:deskew_rotate_file",
                [rest[0], *DSR_ARGS],
                DSR_KWARGS,
            )
            out = os.path.join(os.path.dirname(os.path.abspath(rest[0])), "DSR", "x.tif")
            inp = rest[0]
        else:
            task = Task(
                "petakit5d.tools.crop:crop_file",
                [rest[0], rest[1], list(TINY_BOX)],
                {"overwrite": True},
            )
            inp, out = rest[0], rest[1]
        write_call(Call([inp], [out], [task_record(task, "python")], {}), folder)
        return
    if kind == "dsr":
        _clear_dsr(argv[1])
        from petakit5d.deskew_rotate.file import deskew_rotate_file

        mark("ready")
        deskew_rotate_file(argv[1], *DSR_ARGS, **DSR_KWARGS)
    elif kind == "tiny":
        from petakit5d.tools.crop import crop_file

        mark("ready")
        crop_file(argv[1], argv[2], TINY_BOX, overwrite=True)
    else:
        raise SystemExit(f"unknown kind {kind!r}")


if __name__ == "__main__":
    main(sys.argv[1:])
