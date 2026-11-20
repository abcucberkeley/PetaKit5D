"""Dataset wrappers of ``third_parties/sim3D``: SIM reconstruction (``simRecon_data_wrapper``),
deskewing of phases (``deskewPhases_data_wrapper``) and rotation (``Rotate3D_data_wrapper``) of
every tiff file of the datasets that matches a channel pattern, one task per file through the
job framework (``petakit5d.jobs.run_tasks``), retried once when not every task finished.

The tasks get their values as MATLAB's function strings give them (``%.10f`` text, booleans as
``true``/``false``) and only the parameters the strings name: ``simRecon_data_wrapper`` passes
``lattice_period``, ``phase_step``, ``norders``, ``nphases``, ``overlap``, ``chunk_size``,
``edge_taper``, ``edge_taper_val``, ``perdecomp``, ``use_gpu``, ``ds`` and ``background`` (the
other values of the wrapper are not used: MATLAB's quirk). Results go to ``<path>/sim_recon/``,
``<path>/DS/`` (``DSR`` with ``rotate``) and ``<path>/Rotated/``.

MATLAB's three wrappers fail in today's MATLAB tree, because they still pass ``cpuOnlyNodes``,
which the job framework no longer takes (``MATLAB:InputParser:UnmatchedParameter``); the port
runs them as they were written, GPU tasks as GPU jobs (MATLAB's site-specific ``SlurmParam``
comes from the cluster configuration).
"""

from __future__ import annotations

import math
import os
import warnings
from collections.abc import Sequence
from typing import Any

from petakit5d.params import Param, register

__all__ = ["deskew_phases_dataset", "rotate_3d_dataset", "sim_recon_dataset"]


def _paths(data_paths) -> list[str]:
    if isinstance(data_paths, (str, os.PathLike)):
        return [os.fspath(data_paths)]
    return [os.fspath(p) for p in data_paths]


def _patterns(p) -> list[str]:
    return [p] if isinstance(p, str) else list(p)


def _files(folder: str, pattern: str) -> list[str]:
    """``dir([folder filesep '*' pattern '*.tif'])``: names in MATLAB's order"""
    from petakit5d.core.paths import matlab_dir

    return [e.name for e in matlab_dir(f"{folder}{os.sep}*{pattern}*.tif")]


def _job_log_dir(job_log_dir: str, data_path: str, *, group: bool) -> str:
    """the job log folder, or ``<data path>/job_logs`` (with MATLAB's warning) without one"""
    if os.path.isdir(job_log_dir):
        return job_log_dir
    from petakit5d.core.text import sprintf
    from petakit5d.sim.recon import SimWarning

    warnings.warn(
        SimWarning(
            sprintf(
                "The job log directory does not exist, use %s/job_logs as job log directory.",
                data_path,
            )
        ),
        stacklevel=3,
    )
    d = f"{data_path}/job_logs"
    if not os.path.isdir(d):
        os.makedirs(d, exist_ok=True)
        if group:
            from petakit5d.io.tiff import group_writable

            group_writable(d)
    return d


def _f10(v: float) -> float:
    """a value through ``%.10f`` text"""
    from petakit5d.core.text import through_text

    return float(through_text(float(v), "%.10f"))


def _run(
    inputs,
    outputs,
    tasks,
    *,
    use_gpu: bool,
    cpus_per_task: float,
    task_batch_num: float,
    master_compute: bool,
    parse_cluster: bool,
    job_log_dir: str | None = None,
) -> list[bool]:
    from petakit5d.jobs import framework

    kw: dict[str, Any] = {
        "cpus_per_task": cpus_per_task,
        "max_job_num": math.inf,
        "task_batch_num": task_batch_num,
        "master_compute": master_compute,
        "parse_cluster": parse_cluster,
        "gpu_job": bool(use_gpu),
        "gpu_num": 1 if use_gpu else 0,
    }
    if job_log_dir is not None:
        kw["job_log_dir"] = job_log_dir
    if not inputs:
        return []
    done = framework.run_tasks(inputs, outputs, tasks, **kw)
    if not all(done):
        done = framework.run_tasks(inputs, outputs, tasks, **kw)
    return list(done)


@register(
    "simRecon_data_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param("PSFs", None, "One PSF tiff file per channel pattern.", py="psfs"),
    Param("ChannelPatterns", "@iscell", "Channel identifiers of the files."),
    Param("islattice", "@islogical", "Not passed to the tasks (MATLAB's)."),
    Param("NA_det", "@isnumeric", "Not passed to the tasks."),
    Param("NA_ext", "@isnumeric", "Not passed to the tasks."),
    Param("nimm", "@isnumeric", "Not passed to the tasks."),
    Param("wvl_em", "@isnumeric", "Not passed to the tasks."),
    Param("wvl_ext", "@isnumeric", "Not passed to the tasks."),
    Param("w", "@isnumeric", "Not passed to the tasks."),
    Param("apodize", "@islogical", "Not passed to the tasks."),
    Param("DS", "@islogical", "The data are deskewed."),
    Param("nphases", "@isnumeric", "Phases."),
    Param("norders", "@isnumeric", "Information orders."),
    Param("norientations", "@isnumeric", "Not passed to the tasks."),
    Param("lattice_period", "@isnumeric", "Lattice period in um."),
    Param("lattice_angle", "@isnumeric", "Not passed to the tasks."),
    Param("phase_step", "@isnumeric", "Phase step in um."),
    Param("pxl_dim_data", "@isnumeric", "Not passed to the tasks (z, y, x).", convert="size_yxz"),
    Param("pxl_dim_PSF", "@isnumeric", "Not passed to the tasks (z, y, x).", convert="size_yxz"),
    Param("Background", "@isnumeric", "Background.", convert="same"),
    Param("normalize_orientations", "@islogical", "Not passed to the tasks."),
    Param("perdecomp", "@islogical", "Periodic-plus-smooth decomposition."),
    Param("edgeTaper", "@islogical", "Tukey edge taper."),
    Param("edgeTaperVal", "@isnumeric", "Roll-off of the taper."),
    Param("useGPU", "@islogical", "Run the tasks on GPUs."),
    Param("Overwrite", "@islogical", "Not used."),
    Param("Save16bit", "@islogical", "Not used."),
    Param("flipZstack", "@islogical", "Not used."),
    Param("EdgeSoften", "@isnumeric", "Not used."),
    Param("zEdgeSoften", "@isnumeric", "Not used."),
    Param("Crop", "@isnumeric", "Not used."),
    Param("zFlip", "@islogical", "Not used."),
    Param("GenMaxZproj", "@isnumeric", "Not used.", convert="flags_yxz"),
    Param("ResizeImages", "@isnumeric", "Not used."),
    Param("EdgeErosion", "@isnumeric", "Not used."),
    Param("ErodeMaskfile", "@ischar", "Not used."),
    Param("SaveMaskfile", "@islogical", "Not used."),
    Param("ChunkSize", "@isvector", "Chunk size (z, y, x) of the tasks.", convert="size_yxz"),
    Param("Overlap", "@isnumeric", "Chunk overlap.", convert="same"),
    Param("maxSubVolume", "@isnumeric", "Not used."),
    Param("CPUMaxMem", "@isnumeric", "Not used."),
    Param("splitJobsByChannelPattern", "@islogical", "Not used."),
    Param("largeFile", "@islogical", "Not used."),
    Param("parseCluster", "@islogical", "Run the tasks as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("masterCPU", "@islogical", "Not used."),
    Param("jobLogDir", "@ischar", "Job log folder (<path>/job_logs if it does not exist)."),
    Param("cpusPerTask", "@isnumeric", "Not used (5 for GPU tasks, 24 otherwise: MATLAB's)."),
    Param("uuid", "@ischar", "Not used."),
    Param("maxTrialNum", "@isnumeric", "Not used."),
    Param("unitWaitTime", "@isnumeric", "Not used."),
    task="C06",
)
def sim_recon_dataset(
    data_paths,
    psfs,
    *,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    islattice: bool = True,
    na_det: float = 1.0,
    na_ext: float = 0.55,
    nimm: float = 1.33,
    wvl_em: float = 0.605,
    wvl_ext: float = 0.56,
    w: float = 0.005,
    apodize: bool = True,
    ds: bool = False,
    nphases: float = 5,
    norders: float = 5,
    norientations: float = 1,
    lattice_period: float = 1.2021,
    lattice_angle: Any = 1.5707963267948966,
    phase_step: float = 0.232,
    pxl_dim_data: Sequence[float] = (0.16074803849369898, 0.11, 0.11),
    pxl_dim_psf: Sequence[float] = (0.10716535899579933, 0.11, 0.11),
    background: float = 105,
    normalize_orientations: bool = False,
    perdecomp: bool = False,
    edge_taper: bool = False,
    edge_taper_val: float = 0.1,
    use_gpu: bool = True,
    overwrite: bool = False,
    save_16bit: bool = False,
    flip_zstack: bool = False,
    edge_soften: float = 5,
    z_edge_soften: float = 2,
    crop: Any = None,
    z_flip: bool = False,
    gen_max_zproj: Sequence[float] = (1, 0, 0),
    resize_images: Any = None,
    edge_erosion: float = 0,
    erode_maskfile: str = "",
    save_maskfile: bool = False,
    chunk_size: Sequence[float] = (250, 250, 250),
    overlap: float = 10,
    max_sub_volume: float = 300000000,
    cpu_max_mem: float = 500,
    split_jobs_by_channel_pattern: bool = False,
    large_file: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = False,
    master_cpu: bool = False,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 5,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 2,
) -> list[bool]:
    """SIM reconstructions (``sim_recon_file``) of every ``*<pattern>*.tif`` of the datasets
    with the pattern's PSF, into ``<path>/sim_recon/`` (``simRecon_data_wrapper``); returns the
    done flags of the last pattern's tasks."""
    from petakit5d.core.paths import fileparts
    from petakit5d.jobs.tasks import Task
    from petakit5d.sim.file import sim_recon_file

    psf_list = _paths(psfs)
    done: list[bool] = []
    cs = [_f10(v) for v in chunk_size]
    kw = {
        "lattice_period": _f10(lattice_period),
        "phase_step": _f10(phase_step),
        "norders": _f10(norders),
        "nphases": _f10(nphases),
        "overlap": _f10(overlap),
        "chunk_size": tuple(cs),
        "edge_taper": bool(edge_taper),
        "edge_taper_val": _f10(edge_taper_val),
        "perdecomp": bool(perdecomp),
        "use_gpu": bool(use_gpu),
        "ds": bool(ds),
        "background": _f10(background),
    }
    for d in _paths(data_paths):
        for c, patt in enumerate(_patterns(channel_patterns)):
            names = _files(d, patt)
            log_dir = _job_log_dir(job_log_dir, d, group=False) if parse_cluster else None
            if log_dir is not None:
                job_log_dir = log_dir
            inputs, outputs, tasks = [], [], []
            for name in names:
                _, fsname, ext = fileparts(name)
                full = f"{d}{os.sep}{name}"
                inputs.append(full)
                outputs.append(f"{d}{os.sep}sim_recon{os.sep}{fsname}_recon{ext}")
                tasks.append(Task.of(sim_recon_file, full, psf_list[c], **kw))
            done = _run(
                inputs, outputs, tasks, use_gpu=use_gpu, cpus_per_task=5 if use_gpu else 24,
                task_batch_num=5 if use_gpu else 1, master_compute=master_compute,
                parse_cluster=parse_cluster,
            )  # fmt: skip
    return done


@register(
    "deskewPhases_data_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param("xyPixelSize", None, "Pixel size in um."),
    Param("dz", None, "Scan step in um."),
    Param("SkewAngle", "@isscalar", "Skew angle in degrees."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("nphases", "@isnumeric", "Phases."),
    Param("ChannelPatterns", "@iscell", "Channel identifiers of the files."),
    Param("Rotate", "@islogical", "Rotate after the deskew (DSR)."),
    Param("Overwrite", "@islogical", "Not used."),
    Param("Save16bit", "@islogical", "Not passed to the tasks (MATLAB's)."),
    Param("flipZstack", "@islogical", "Not used."),
    Param("EdgeSoften", "@isnumeric", "Not used."),
    Param("zEdgeSoften", "@isnumeric", "Not used."),
    Param("Crop", "@isnumeric", "Not used."),
    Param("zFlip", "@islogical", "Not used."),
    Param("GenMaxZproj", "@isnumeric", "Not used.", convert="flags_yxz"),
    Param("ResizeImages", "@isnumeric", "Not used."),
    Param("EdgeErosion", "@isnumeric", "Not used."),
    Param("ErodeMaskfile", "@ischar", "Not used."),
    Param("SaveMaskfile", "@islogical", "Not used."),
    Param("ChunkSize", "@isvector", "Not used.", convert="size_yxz"),
    Param("Overlap", "@isnumeric", "Not used.", convert="same"),
    Param("maxSubVolume", "@isnumeric", "Not used."),
    Param("CPUMaxMem", "@isnumeric", "Not used."),
    Param("largeFile", "@islogical", "Not used."),
    Param("parseCluster", "@islogical", "Run the tasks as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("masterCPU", "@islogical", "Not used."),
    Param("jobLogDir", "@ischar", "Job log folder (<path>/job_logs if it does not exist)."),
    Param("cpusPerTask", "@isnumeric", "Not used (24: MATLAB's)."),
    Param("uuid", "@ischar", "Not used."),
    Param("maxTrialNum", "@isnumeric", "Not used."),
    Param("unitWaitTime", "@isnumeric", "Not used."),
    task="C06",
)
def deskew_phases_dataset(
    data_paths,
    xy_pixel_size,
    dz,
    skew_angle: float = 32.45,
    reverse: bool = True,
    *,
    nphases: float = 5,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    rotate: bool = False,
    overwrite: bool = False,
    save_16bit: bool = False,
    flip_zstack: bool = False,
    edge_soften: float = 5,
    z_edge_soften: float = 2,
    crop: Any = None,
    z_flip: bool = False,
    gen_max_zproj: Sequence[float] = (1, 0, 0),
    resize_images: Any = None,
    edge_erosion: float = 0,
    erode_maskfile: str = "",
    save_maskfile: bool = False,
    chunk_size: Sequence[float] = (250, 250, 250),
    overlap: float = 10,
    max_sub_volume: float = 300000000,
    cpu_max_mem: float = 500,
    large_file: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = False,
    master_cpu: bool = False,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 5,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 2,
) -> list[bool]:
    """Deskew every ``*<pattern>*.tif`` of the datasets phase by phase (``deskew_phases_file``)
    into ``<path>/DS/`` (``DSR`` with ``rotate``; ``deskewPhases_data_wrapper``); returns the
    done flags of the last pattern's tasks."""
    from petakit5d.core.paths import fileparts
    from petakit5d.jobs.tasks import Task
    from petakit5d.sim.file import deskew_phases_file

    fol = "DSR" if rotate else "DS"
    done: list[bool] = []
    kw = {
        "skew_angle": _f10(skew_angle),
        "reverse": bool(reverse),
        "nphases": _f10(nphases),
        "rotate": bool(rotate),
    }
    for d in _paths(data_paths):
        for patt in _patterns(channel_patterns):
            names = _files(d, patt)
            if parse_cluster:
                job_log_dir = _job_log_dir(job_log_dir, d, group=True)
            inputs, outputs, tasks = [], [], []
            for name in names:
                _, fsname, ext = fileparts(name)
                full = f"{d}{os.sep}{name}"
                inputs.append(full)
                outputs.append(f"{d}{os.sep}{fol}{os.sep}{fsname}{ext}")
                tasks.append(Task.of(deskew_phases_file, full, _f10(xy_pixel_size), _f10(dz), **kw))
            done = _run(
                inputs, outputs, tasks, use_gpu=False, cpus_per_task=24, task_batch_num=1,
                master_compute=master_compute, parse_cluster=parse_cluster,
            )  # fmt: skip
    return done


@register(
    "Rotate3D_data_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param("xyPixelSize", None, "Pixel size in um."),
    Param("dz", None, "Scan step in um."),
    Param("SkewAngle", "@isscalar", "Skew angle in degrees."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("ChannelPatterns", "@iscell", "Channel identifiers of the files."),
    Param("ObjectiveScan", "@islogical", "Objective scan."),
    Param("sCMOSCameraFlip", "@islogical", "Not used."),
    Param("Save16bit", "@islogical", "Save as uint16."),
    Param("Overwrite", "@islogical", "Not used."),
    Param("flipZstack", "@islogical", "Not used."),
    Param("EdgeSoften", "@isnumeric", "Not used."),
    Param("zEdgeSoften", "@isnumeric", "Not used."),
    Param("Crop", "@isnumeric", "Not used."),
    Param("maxSubVolume", "@isnumeric", "Not used."),
    Param("CPUMaxMem", "@isnumeric", "Not used."),
    Param("largeFile", "@islogical", "Not used."),
    Param("parseCluster", "@islogical", "Run the tasks as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("masterCPU", "@islogical", "Not used."),
    Param("jobLogDir", "@ischar", "Job log folder (<path>/job_logs if it does not exist)."),
    Param("cpusPerTask", "@isnumeric", "Not used (24: MATLAB's)."),
    Param("uuid", "@ischar", "Not used."),
    Param("maxTrialNum", "@isnumeric", "Not used."),
    Param("unitWaitTime", "@isnumeric", "Not used."),
    task="C06",
)
def rotate_3d_dataset(
    data_paths,
    xy_pixel_size,
    dz,
    skew_angle: float = 32.45,
    reverse: bool = True,
    *,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    objective_scan: bool = False,
    s_cmos_camera_flip: bool = False,
    save_16bit: bool = False,
    overwrite: bool = False,
    flip_zstack: bool = False,
    edge_soften: float = 5,
    z_edge_soften: float = 2,
    crop: Any = None,
    max_sub_volume: float = 300000000,
    cpu_max_mem: float = 500,
    large_file: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = False,
    master_cpu: bool = False,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 5,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 2,
) -> list[bool]:
    """Rotate every ``*<pattern>*.tif`` of the datasets (``rotate_file``, MATLAB's
    ``XR_RotateFrame3D``) into ``<path>/Rotated/`` (``Rotate3D_data_wrapper``); returns the
    done flags of the last pattern's tasks."""
    from petakit5d.core.paths import fileparts
    from petakit5d.deskew_rotate.file import rotate_file
    from petakit5d.jobs.tasks import Task

    done: list[bool] = []
    kw = {
        "skew_angle": _f10(skew_angle),
        "reverse": bool(reverse),
        "objective_scan": bool(objective_scan),
        "save_16bit": bool(save_16bit),
    }
    for d in _paths(data_paths):
        for patt in _patterns(channel_patterns):
            names = _files(d, patt)
            if parse_cluster:
                job_log_dir = _job_log_dir(job_log_dir, d, group=False)
            inputs, outputs, tasks = [], [], []
            for name in names:
                _, fsname, ext = fileparts(name)
                full = f"{d}{os.sep}{name}"
                inputs.append(full)
                outputs.append(f"{d}{os.sep}Rotated{os.sep}{fsname}{ext}")
                tasks.append(Task.of(rotate_file, full, _f10(xy_pixel_size), _f10(dz), **kw))
            done = _run(
                inputs, outputs, tasks, use_gpu=False, cpus_per_task=24, task_batch_num=1,
                master_compute=master_compute, parse_cluster=parse_cluster,
            )  # fmt: skip
    return done
