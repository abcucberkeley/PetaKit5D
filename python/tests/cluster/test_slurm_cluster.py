"""The job framework in Slurm cluster mode on a real Slurm (SlurmClusterTest, task J11).

Tasks are submitted with sbatch, polled with squeue, retried with more CPUs and cancelled at
the end. The tests submit real jobs: they run only with PETAKIT5D_TEST_SLURM_CONFIG (see
conftest.py for the variables, the job names and the clean-up).

The tasks record the Slurm environment of the job that ran them (array job id, CPUs per task,
time limit) in their outputs. Each test of MATLAB's bash tasks runs with bash tasks and with
Python tasks (task records of a module in the test's folder, run by the job entry point); the
two tests whose MATLAB jobs started MATLAB run Python tasks here (MATLAB jobs could not run on
GitHub's runners; Python ones can). Python tasks run with the task runners (the default), except
where a test checks which job runs which task: runners take tasks of jobs that waited 30 s or
more, which a busy cluster makes happen.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest

from petakit5d.jobs.slurm import GONE, PENDING, RUNNING, split_options

pytestmark = pytest.mark.cluster

LANGUAGES = ("bash", "python")


def wait_for_status(cluster, job_id: int, array_id: int, target: int, seconds: float) -> int:
    deadline = time.monotonic() + seconds
    while True:
        status = cluster.slurm.batch_status([job_id], [array_id])[0]
        if status == target or time.monotonic() > deadline:
            return status
        time.sleep(2)


def ran_in_jobs(jobs: list[str]) -> bool:
    return all(j != "" for j in jobs)


def runner_off(language: str) -> dict:
    """for tests that check which job runs which task (see the module's description)"""
    return {"task_runner": False} if language == "python" else {}


# --- Unit ------------------------------------------------------------------------------------


def test_job_status_follows_squeue(cluster):
    # check_batch_slurm_jobs_status (Slurm.batch_status) parses squeue's default output
    folder = cluster.shared_folder()
    argv = [
        "sbatch", "--hold", "--array=3", "-J", cluster.job_name,
        "-o", str(folder / "job_%A_%a.out"), *split_options(cluster.conf.get("SlurmParam", "")),
        "-t", "0:05:00", "--wrap=sleep 120",
    ]  # fmt: skip
    r = cluster.slurm.sbatch(argv)
    assert r.submitted, r.result.output
    job = r.job_id

    assert cluster.slurm.batch_status([job], [3]) == [PENDING], "A held job should be pending."
    assert cluster.slurm.run(["scontrol", "release", f"{job}_3"]).ok
    assert wait_for_status(cluster, job, 3, RUNNING, 300) == RUNNING, (
        "The released job should be running."
    )
    cluster.slurm.cancel([job])
    assert wait_for_status(cluster, job, 3, GONE, 60) == GONE, "A cancelled job should be gone."


# --- Integration -----------------------------------------------------------------------------


@pytest.mark.parametrize("language", LANGUAGES)
def test_tasks_run_as_slurm_jobs(cluster, language):
    folder = cluster.shared_folder()
    inputs, outputs, tasks = cluster.tasks(folder, 3, language)
    done, fired = cluster.run_framework(
        folder, inputs, outputs, tasks, language=language, **runner_off(language)
    )

    assert not fired, "Jobs did not finish within the watchdog's limit and were cancelled."
    assert done == [True] * 3
    jobs = cluster.read_field(outputs, "array_job")
    assert ran_in_jobs(jobs), "Some tasks did not run as Slurm jobs."
    assert len(set(jobs)) == 3, "Each task should run as its own job."
    assert list((folder / "job_logs").glob("job_*.out")), "No job logs in jobLogDir."


@pytest.mark.parametrize("language", LANGUAGES)
def test_tasks_in_a_batch_share_a_job(cluster, language):
    folder = cluster.shared_folder()
    inputs, outputs, tasks = cluster.tasks(folder, 4, language)
    done, fired = cluster.run_framework(
        folder, inputs, outputs, tasks, language=language, task_batch_num=2,
        **runner_off(language),
    )  # fmt: skip

    assert not fired
    assert done == [True] * 4
    jobs = cluster.read_field(outputs, "array_job")
    assert ran_in_jobs(jobs), "Some tasks did not run as Slurm jobs."
    assert jobs[0] == jobs[1], "Tasks 1 and 2 should share a job."
    assert jobs[2] == jobs[3], "Tasks 3 and 4 should share a job."
    assert jobs[0] != jobs[2], "Each batch should run as its own job."


@pytest.mark.parametrize("language", LANGUAGES)
def test_existing_outputs_are_not_submitted(cluster, language):
    folder = cluster.shared_folder()
    inputs, outputs, tasks = cluster.tasks(folder, 3, language)
    Path(outputs[1]).write_text("existing\n")
    done, fired = cluster.run_framework(folder, inputs, outputs, tasks, language=language)

    assert not fired
    assert done == [True] * 3
    assert Path(outputs[1]).read_text().strip() == "existing", "A task with an output was run."
    jobs = cluster.read_field([outputs[0], outputs[2]], "array_job")
    assert ran_in_jobs(jobs), "The other tasks did not run as Slurm jobs."


@pytest.mark.parametrize("language", LANGUAGES)
def test_master_computes_part_of_the_tasks(cluster, language):
    # with master_compute the caller runs a batch itself instead of waiting for its job
    folder = cluster.shared_folder()
    inputs, outputs, tasks = cluster.tasks(folder, 3, language)
    done, fired = cluster.run_framework(
        folder, inputs, outputs, tasks, language=language, master_compute=True
    )

    assert not fired
    assert done == [True] * 3
    jobs = cluster.read_field(outputs, "array_job")
    assert "" in jobs, "No task was run by the caller itself."


@pytest.mark.parametrize("language", LANGUAGES)
def test_requested_resources_reach_the_job(cluster, language):
    if cluster.max_cpu_num < 2:
        pytest.skip("The cluster configuration allows only 1 CPU per task.")
    folder = cluster.shared_folder()
    inputs, outputs, tasks = cluster.tasks(folder, 2, language)
    done, fired = cluster.run_framework(
        folder, inputs, outputs, tasks, language=language, cpus_per_task=2, job_time_limit=0.25
    )

    assert not fired
    assert done == [True] * 2
    assert cluster.read_field(outputs, "cpus") == ["2", "2"]
    assert cluster.read_field(outputs, "time_limit") == ["00:15:00", "00:15:00"]


@pytest.mark.parametrize("language", LANGUAGES)
def test_final_output_ends_the_run(cluster, language):
    # the first task writes the final output; the others would run for minutes, but the
    # framework stops once the final output exists
    folder = cluster.shared_folder()
    final = str(folder / "final.txt")
    inputs, outputs, tasks = cluster.tasks(folder, 3, language, [0, 600, 600], final=final)
    t0 = time.monotonic()
    done, fired = cluster.run_framework(
        folder, inputs, outputs, tasks, language=language, final_out_fullpath=final
    )

    assert not fired
    assert time.monotonic() - t0 < 300, "The framework kept waiting after the final output."
    assert done == [True] * 3


@pytest.mark.parametrize("language", LANGUAGES)
def test_final_output_cancels_remaining_jobs(cluster, language, capfd):
    # once the final output exists, the framework returns and cancels the jobs still running
    folder = cluster.shared_folder()
    final = str(folder / "final.txt")
    inputs, outputs, tasks = cluster.tasks(folder, 3, language, [0, 600, 600], final=final)
    try:
        cluster.run_framework(
            folder, inputs, outputs, tasks, language=language, final_out_fullpath=final
        )
    finally:
        printed = capfd.readouterr().out
    assert len(cluster.job_ids_of(printed)) == 3, printed[-3000:]

    # cancelled jobs can take a few seconds to leave the queue
    left = cluster.wait_until_gone(60)
    assert left == [], f"Jobs still queued after the final output existed: {left}"


# --- Integration, slow -----------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("language", LANGUAGES)
def test_failed_task_is_retried_with_more_cpus(cluster, language):
    # a job that ends without its output is submitted again (after a grace period) with
    # cpus_per_task * (try + 1) CPUs, up to max_trial_num tries
    if cluster.max_cpu_num < 2:
        pytest.skip("The cluster configuration allows only 1 CPU per task.")
    folder = cluster.shared_folder()
    inputs, outputs, tasks = cluster.tasks(folder, 2, language, fail=[False, True])
    done, fired = cluster.run_framework(
        folder, inputs, outputs, tasks, language=language, max_trial_num=2
    )

    assert not fired
    assert done == [True, False]
    attempts = (folder / "attempts_2.txt").read_text().split()
    assert attempts == ["cpus=1", "cpus=2"], (
        "The failed task should be tried twice, with 1 and then 2 CPUs."
    )


# --- Python tasks where MATLAB's jobs started MATLAB -------------------------------------------


def test_python_tasks_run_as_jobs(cluster):
    # (matlabTasksRunAsJobs) jobs start Python with PythonLaunchStr and run the task records:
    # a function of a module outside the package, found through the record's path
    from petakit5d.jobs.tasks import Task

    folder = cluster.shared_folder()
    log = folder / "calls.log"
    inputs, outputs, tasks = [], [], []
    for i in (1, 2):
        inp, out = folder / f"in_{i}.txt", folder / f"out_{i}.txt"
        inp.write_text(f"task {i}\n")
        inputs.append(str(inp))
        outputs.append(str(out))
        tasks.append(
            Task(f"{cluster.task_module}:copy", [str(inp), str(out), str(log)], path=[str(folder)])
        )
    done, fired = cluster.run_framework(folder, inputs, outputs, tasks, language="python")

    assert not fired
    assert done == [True, True]
    for inp, out in zip(inputs, outputs, strict=True):
        assert Path(out).read_text() == Path(inp).read_text()
    calls = log.read_text().splitlines()
    assert len(calls) == 2, "Each task should run once."
    # the outputs alone would also appear if the framework ran the tasks in this process
    assert all(c.split(" ")[1:] != [""] for c in calls), calls
    assert cluster.count_job_logs(folder / "job_logs", str(folder)) >= 2, (
        "The tasks did not run as Slurm jobs."
    )


@pytest.mark.slow
def test_deskew_rotate_dataset_runs_as_jobs(cluster):
    # (deskewRotateWrapperRunsAsJobs) a data wrapper end to end on the cluster gives the
    # in-memory function's result
    import json

    from petakit5d.core.dtypes import cast
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.deskew_rotate.dataset import deskew_rotate_dataset
    from petakit5d.io import read_tiff, write_tiff
    from petakit5d.testing import synth

    folder = cluster.shared_folder()
    # the wrapper takes no SlurmConstraint: its jobs get the test's name through the config
    conf = dict(cluster.conf)
    conf["SlurmParam"] = f"{conf.get('SlurmParam', '')} -J {cluster.job_name}"
    config = folder / "config.json"
    config.write_text(json.dumps(conf))
    fsn = (
        "Scan_Iter_0000_0000_Cam{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs"
        "_000x_000y_000z_0000t"
    )
    frames = {}
    for i, cam in enumerate(("A", "B"), start=1):
        vol = synth.beads((40, 32, 48), [(11, 15, 19), (27, 9, 29)], 2, 1000.0, 100.0)
        frames[cam] = cast(vol + np.float32(50 * i), "uint16")
        write_tiff(frames[cam], str(folder / f"{fsn.format(cam)}.tif"))
    deskew_rotate_dataset(
        f"{folder}/", deskew=True, rotate=True, dsr_combined=True, xy_pixel_size=0.108, dz=0.2,
        skew_angle=32.45, reverse=True, channel_patterns=("CamA", "CamB"), save_16bit=True,
        parse_cluster=True, master_compute=False, config_file=str(config),
        job_log_dir=str(folder / "job_logs"),
    )  # fmt: skip

    for cam, frame in frames.items():
        expected = deskew_rotate_volume(frame, 32.45, 0.2, 0.108, reverse=True, save_16bit=True)
        got = read_tiff(str(folder / "DSR" / f"{fsn.format(cam)}.tif"))
        assert got.dtype == expected.dtype and got.shape == expected.shape
        np.testing.assert_array_equal(got, expected)
    # the frames were made by jobs whose logs are in job_log_dir (the master computes none; a
    # job's task runner may take the other frame's task, so one job can make both)
    logs = [p.read_text(errors="replace") for p in (folder / "job_logs").glob("job_*.out")]
    made = sum(t.count(" done in ") for t in logs if str(folder) in t)
    assert made >= 2, (
        "The frames were not processed as Slurm jobs, or their logs are not in jobLogDir."
    )
