"""One end-to-end workflow in cluster mode (or locally without a config file), as the port runs
it: the same steps and parameters as ``e2e_cluster_workflow.m`` (MATLAB jobs or compiled jobs),
for the speed comparison of the three implementations.

    python e2e_cluster_workflow.py WORKFLOW DATA_PATH CONFIG_FILE TIMING_FILE [EXAMPLES]

WORKFLOW: ``generic`` (demo_generic_computing_framework's 16 tasks: crop and deskew/rotate each
tile, the user function from ``EXAMPLES``, the folder of ``crop_deskew_rotate_demo_function.py``),
``dsr`` (demo_geometric_transformation's dataset deskew/rotation, zarr results) or ``dsr_large``
(the large-file deskew/rotation of the zarr volumes of two channels CamA_ch0 and CamB_ch0,
xy pixel size 0.098, dz 0.5). ``CONFIG_FILE`` ``''``: local mode. Appends ``<workflow>
<seconds>`` to ``TIMING_FILE``.
"""

from __future__ import annotations

import os
import sys
import time


def main(argv: list[str]) -> None:
    workflow, data_path, config_file, timing_file = argv[:4]
    parse_cluster = bool(config_file)
    if not data_path.endswith("/"):
        data_path += "/"
    t0 = time.perf_counter()
    if workflow == "generic":
        from petakit5d.core.paths import matlab_dir
        from petakit5d.core.text import sprintf
        from petakit5d.io import get_image_size
        from petakit5d.jobs import run_tasks
        from petakit5d.jobs.tasks import Task

        examples = os.path.abspath(argv[4])
        fsns = [e.name for e in matlab_dir(data_path + "*.tif")]
        inputs = [data_path + f for f in fsns]
        out_path = data_path + "Cropped_DSR/"
        os.makedirs(out_path, exist_ok=True)
        outputs = [out_path + f for f in fsns]
        # MATLAB's task text has the parameters through %d
        angle, px, dz = (float(sprintf("%d", v)) for v in (32.45, 0.108, 0.3))
        nz, ny, nx = get_image_size(inputs[0])
        box = [0, 100, 0, min(500, nz), min(1700, ny), min(512, nx)]
        tasks = [
            Task(
                "crop_deskew_rotate_demo_function:crop_deskew_rotate_demo_function",
                [fi, fo, angle, px, dz, box],
                path=[examples],
            )
            for fi, fo in zip(inputs, outputs, strict=True)
        ]
        run_tasks(
            inputs,
            outputs,
            tasks,
            parse_cluster=parse_cluster,
            master_compute=True,
            mem_allocate=20,
            config_file=config_file,
            mcc_mode=False,
        )
    elif workflow in ("dsr", "dsr_large"):
        from petakit5d.deskew_rotate import deskew_rotate_dataset

        if workflow == "dsr":
            kw = {
                "deskew": True,
                "rotate": True,
                "dsr_combined": True,
                "xy_pixel_size": 0.108,
                "dz": 0.3,
                "skew_angle": 32.45,
                "objective_scan": False,
                "reverse": True,
                "channel_patterns": ["CamA", "CamB"],
                "large_file": False,
                "zarr_file": False,
                "save_zarr": True,
                "block_size": (256, 256, 256),
                "save_16bit": True,
            }
        else:
            kw = {
                "xy_pixel_size": 0.098,
                "dz": 0.5,
                "reverse": True,
                "channel_patterns": ["CamA_ch0", "CamB_ch0"],
                "large_file": True,
                "zarr_file": True,
                "save_zarr": True,
                "save_16bit": True,
            }
        deskew_rotate_dataset(
            [data_path],
            parse_cluster=parse_cluster,
            master_compute=True,
            config_file=config_file,
            mcc_mode=False,
            **kw,
        )
    else:
        raise SystemExit(f"unknown workflow {workflow!r}")
    with open(timing_file, "a", encoding="utf-8", newline="\n") as fh:
        fh.write(f"{workflow} {time.perf_counter() - t0:.3f}\n")


if __name__ == "__main__":
    main(sys.argv[1:])
