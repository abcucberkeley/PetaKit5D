"""Run a user-defined function on many files with the generic computing framework
(MATLAB: ``demos/demo_generic_computing_framework.m``).

The framework runs one task per file, here on this machine (``parse_cluster=False``), or as
Slurm jobs with a cluster configuration file. The user function is
``crop_deskew_rotate_demo_function.py`` next to this file. Results:
``<data>/PetaKit5D_demo_cell_image_dataset/Cropped_DSR/``.
"""

from __future__ import annotations

import json
import os

import _common

DEMO = "demo_generic_computing_framework"
HERE = os.path.dirname(os.path.abspath(__file__))


def write_example_config(path: str) -> None:
    """An example configuration file for a Slurm cluster (skip this when running locally).

    The keys are the framework's parameters (``petakit5d.jobs.run_tasks``), in MATLAB's names;
    ``PythonLaunchStr`` is what a job runs before Python (module loads, the environment) and
    the interpreter. The values below are placeholders for your cluster.
    """
    config = {
        # what a worker job runs: activate the environment, then the interpreter
        "PythonLaunchStr": "source /opt/envs/petakit5d/bin/activate; python",
        # bash commands run before the tasks of a job (here: GNU parallel for whole-node jobs)
        "BashLaunchStr": "module load parallel",
        # Slurm parameters of the jobs: partition, qos, memory per CPU
        "SlurmParam": "-p compute --qos normal -n1 --mem-per-cpu=21418M",
        # the folder of the job logs, outside the code
        "jobLogDir": "/mnt/projects/job_logs",
        # largest and smallest number of CPU cores per job (at most a node's)
        "maxCPUNum": 24,
        "minCPUNum": 1,
        # allocate whole nodes (and run several tasks in one with GNU parallel), or cores
        "wholeNodeJob": True,
        # memory per CPU core in GB (502 GB on 24 cores: 20.9)
        "memPerCPU": 20.9,
        "GNUparallel": True,
        # job time limit in hours, and the most jobs submitted at once
        "jobTimeLimit": 12,
        "maxJobNum": 5000,
    }
    with open(path, "w", newline="\n") as f:
        json.dump(config, f, indent=4)


def main() -> None:
    from petakit5d.core.paths import matlab_dir
    from petakit5d.core.text import sprintf
    from petakit5d.io import get_image_size
    from petakit5d.jobs import run_tasks
    from petakit5d.jobs.tasks import Task

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()

    with steps("example configuration file for a Slurm cluster"):
        # MATLAB writes it to demos/ of the code; this example to the current folder
        write_example_config(os.path.join(os.getcwd(), "test_slurm_cpu_config_file.json"))

    with steps("crop and deskew/rotate the tiles with the generic computing framework"):
        fsns = [e.name for e in matlab_dir(data_path + "*.tif")]
        print(fsns)
        input_fullnames = [data_path + f for f in fsns]

        out_path = data_path + "Cropped_DSR/"
        os.makedirs(out_path, exist_ok=True)
        output_fullnames = [out_path + f for f in fsns]

        # skew angle, xy pixel size and z scan step; MATLAB passes them in the task's text with
        # %d (so 32.45 becomes 3.245000e+01), which is kept here
        skew_angle = float(sprintf("%d", 32.45))
        xy_pixel_size = float(sprintf("%d", 0.108))
        dz = float(sprintf("%d", 0.3))
        # the crop box (z0, y0, x0, z1, y1, x1); MATLAB's [ymin xmin zmin ymax xmax zmax] is
        # [101, 1, 1, 1700, 512, 500], its end limited to the image size
        nz, ny, nx = get_image_size(input_fullnames[0])
        in_bbox = [0, 100, 0, min(500, nz), min(1700, ny), min(512, nx)]

        # one task per file: the function by its import path, and its arguments
        tasks = [
            Task(
                "crop_deskew_rotate_demo_function:crop_deskew_rotate_demo_function",
                [fi, fo, skew_angle, xy_pixel_size, dz, in_bbox],
                path=[HERE],
            )
            for fi, fo in zip(input_fullnames, output_fullnames, strict=True)
        ]

        # parse_cluster=False: run on this machine; True: Slurm jobs set up by config_file.
        # master_compute: this process also runs tasks; mem_allocate: the estimated memory of a
        # task in GB
        run_tasks(
            input_fullnames,
            output_fullnames,
            tasks,
            parse_cluster=False,
            master_compute=True,
            mem_allocate=20,
            config_file="",
            mcc_mode=False,
        )


if __name__ == "__main__":
    _common.run_demo(main)
