"""The SIM pipeline of ``third_parties/sim3D`` (``simReconAutomaticProcessing``): deskew the
phases of every matching file, then reconstruct the deskewed files, each stage as tasks of the
job framework run in the background, optionally streaming (new files as they arrive).

Each loop (every 5 s after the first) lists, per data path and channel pattern, the
``*<pattern>*.tif`` files not handled yet (while streaming without the newest one, which may
still be written, and without files younger than ``unit_wait_time`` minutes for the
reconstruction), and hands them as one job to a background worker (MATLAB's ``parfeval`` on a
pool of ``par_pool_size`` workers; here threads): ``deskew_phases_file`` into ``<path>/DS``
(``DSR`` with ``rotate``), ``sim_recon_file`` from there into ``<path>/DS/<results_dir_name>``
(the data path itself without ``deskew``). A stage whose results all exist is skipped. The loop
ends when every worker has finished and, while streaming, when no file changed for
``max_modify_time`` minutes. A worker that finished reports "Finished ..." (also when it failed:
MATLAB's, whose futures' errors are not read).

The values reach the tasks as MATLAB's function strings give them (``%.10f`` text). MATLAB's
own runs of this function fail in today's tree: its framework call passes ``cpuOnlyNodes``,
which the framework no longer takes (an error each worker keeps to itself, so MATLAB reports
the stages finished without results); the port runs the stages as written.
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Sequence
from typing import Any

from petakit5d.params import Param, register

__all__ = ["sim_recon_automatic_processing"]

_PAUSE = 5.0  # seconds between loops (MATLAB's pause(5))


def _minutes_old(paths: Sequence[str]) -> list[float]:
    now = time.time()
    return [(now - os.path.getmtime(p)) / 60 for p in paths]


@register(
    "simReconAutomaticProcessing",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders."),
    Param(
        "PSFs", "@(x) ischar(x) || iscell(x)", "One PSF tiff file per channel pattern.", py="psfs"
    ),
    Param("Deskew", "@islogical", "Deskew the phases first."),
    Param("Recon", "@islogical", "Reconstruct."),
    Param("Streaming", "@islogical", "Keep processing new files as they arrive."),
    Param("resultsDirName", "@ischar", "Result folder of the reconstructions."),
    Param("reconBatchNum", "@isnumeric", "Reconstructions per GPU job."),
    Param("parPoolSize", "@isnumeric", "Background workers.", convert="same"),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Scan step in um."),
    Param("SkewAngle", "@isscalar", "Skew angle in degrees."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("Rotate", "@islogical", "Rotate after the deskew (DSR)."),
    Param("ChannelPatterns", "@iscell", "Channel identifiers of the files."),
    Param("islattice", "@islogical", "Lattice light sheet data (else 3-beam SIM)."),
    Param("NA_det", "@isnumeric", "Detection NA."),
    Param("NA_ext", "@isnumeric", "Excitation NA."),
    Param("nimm", "@isnumeric", "Refractive index of the immersion medium."),
    Param("wvl_em", "@isnumeric", "Emission wavelength in um."),
    Param("wvl_ext", "@isnumeric", "Excitation wavelength in um."),
    Param("w", "@isnumeric", "Wiener coefficient for regularization."),
    Param("apodize", "@islogical", "Apodize the final spectrum."),
    Param("DS", "@islogical", "The reconstructed data are deskewed."),
    Param("nphases", "@isnumeric", "Phases."),
    Param("norders", "@isnumeric", "Information orders."),
    Param("norientations", "@isnumeric", "Pattern orientations."),
    Param("lattice_period", "@isnumeric", "Lattice period in um (the coarsest)."),
    Param("lattice_angle", "@isnumeric", "Angle of the pattern modulation (rad)."),
    Param("phase_step", "@isnumeric", "Phase step in um."),
    Param(
        "pxl_dim_data", "@isnumeric", "Voxel size of the data (z, y, x) in um.", convert="size_yxz"
    ),
    Param(
        "pxl_dim_PSF", "@isnumeric", "Voxel size of the PSF (z, y, x) in um.", convert="size_yxz"
    ),
    Param("Background", "@isnumeric", "Background subtracted.", convert="same"),
    Param("normalize_orientations", "@islogical", "Normalize the total intensity per orientation."),
    Param("perdecomp", "@islogical", "Periodic-plus-smooth decomposition."),
    Param("edgeTaper", "@islogical", "Tukey edge taper."),
    Param("edgeTaperVal", "@isnumeric", "Roll-off of the taper."),
    Param("useGPU", "@islogical", "Reconstruct on GPUs."),
    Param("Overwrite", "@islogical", "Not used."),
    Param("Save16bit", "@islogical", "Save as uint16."),
    Param("gpuPrecision", "@ischar", "Precision of GPU arrays and phase ramps."),
    Param("flipZstack", "@islogical", "Not used."),
    Param("EdgeSoften", "@isnumeric", "Not used."),
    Param("zEdgeSoften", "@isnumeric", "Not used."),
    Param("Crop", "@isnumeric", "Not used."),
    Param("zFlip", "@islogical", "Not used."),
    Param("GenMaxZproj", "@isnumeric", "Not used.", convert="flags_yxz"),
    Param("ResizeImages", "@isnumeric", "Not used."),
    Param("EdgeErosion", "@isnumeric", "Radius of the erosion of the data's mask (0: none)."),
    Param("ErodeBefore", "@islogical", "Apply the eroded mask to the raw data."),
    Param("ErodeAfter", "@islogical", "Apply the eroded mask to the result."),
    Param("ErodeMaskfile", "@ischar", "An existing mask file turns the erosion off."),
    Param("SaveMaskfile", "@islogical", "Save the eroded masks."),
    Param("ChunkSize", "@isvector", "Chunk size (z, y, x) of one phase.", convert="size_yxz"),
    Param("Overlap", "@isnumeric", "Chunk overlap.", convert="same"),
    Param("maxSubVolume", "@isnumeric", "Not used."),
    Param("CPUMaxMem", "@isnumeric", "Not used."),
    Param("splitJobsByChannelPattern", "@islogical", "Not used."),
    Param("largeFile", "@islogical", "Not used."),
    Param("parseCluster", "@islogical", "Run the tasks as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("masterCPU", "@islogical", "Not used."),
    Param("jobLogDir", "@ischar", "Job log folder (<path>/job_logs if it does not exist)."),
    Param("cpusPerTask", "@isnumeric", "CPUs per task (5 for GPU reconstructions)."),
    Param("uuid", "@ischar", "Not used."),
    Param("maxModifyTime", "@isnumeric", "Streaming ends after this many minutes without changes."),
    Param("maxTrialNum", "@isnumeric", "Passed to the reconstructions (not used there)."),
    Param("unitWaitTime", "@isnumeric", "Minutes a new file must be old to be processed."),
    Param("MatlabLaunchStr", "@ischar", "MATLAB-only: accepted and not used."),
    Param(
        "SlurmParam",
        "@ischar",
        "Not used (the stages' jobs take the cluster configuration's).",
        differs="site-specific MATLAB default (from the cluster config)",
    ),
    Param("intThresh", "@isnumeric", "A chunk is reconstructed if its sum is above this."),
    Param("occThresh", "@isnumeric", "... and its nonzero fraction at least this."),
    task="C06",
)
def sim_recon_automatic_processing(
    data_paths,
    *,
    psfs: Any = ("",),
    deskew: bool = True,
    recon: bool = True,
    streaming: bool = False,
    results_dir_name: str = "sim_recon",
    recon_batch_num: float = 5,
    par_pool_size: float = 24,
    xy_pixel_size: float = 0.108,
    dz: float = 0.5,
    skew_angle: float = 32.45,
    reverse: bool = True,
    rotate: bool = False,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    islattice: bool = True,
    na_det: float = 1.0,
    na_ext: float = 0.55,
    nimm: float = 1.33,
    wvl_em: float = 0.605,
    wvl_ext: float = 0.56,
    w: float = 0.005,
    apodize: bool = True,
    ds: bool = True,
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
    gpu_precision: str = "single",
    flip_zstack: bool = False,
    edge_soften: float = 5,
    z_edge_soften: float = 2,
    crop: Any = None,
    z_flip: bool = False,
    gen_max_zproj: Sequence[float] = (1, 0, 0),
    resize_images: Any = None,
    edge_erosion: float = 0,
    erode_before: bool = False,
    erode_after: bool = True,
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
    max_modify_time: float = 10,
    max_trial_num: float = 3,
    unit_wait_time: float = 2,
    matlab_launch_str: str = (
        "module load matlab/r2022a; matlab -nodisplay -nosplash -nodesktop -nojvm -r"
    ),
    slurm_param: str = "",
    int_thresh: float = 1,
    occ_thresh: float = 0.8,
) -> None:
    """Deskew and reconstruct the SIM datasets (``simReconAutomaticProcessing``): the stages'
    tasks in background jobs, optionally streaming."""
    from concurrent.futures import ThreadPoolExecutor

    from petakit5d.core.paths import fileparts, matlab_dir
    from petakit5d.core.text import sprintf
    from petakit5d.core.threads import available_cpus
    from petakit5d.jobs.tasks import Task
    from petakit5d.sim.dataset import _f10, _job_log_dir, _paths, _patterns, _run
    from petakit5d.sim.file import deskew_phases_file, sim_recon_file

    del overwrite, flip_zstack, edge_soften, z_edge_soften, crop, z_flip, gen_max_zproj
    del resize_images, max_sub_volume, cpu_max_mem, split_jobs_by_channel_pattern, large_file
    del master_cpu, uuid, matlab_launch_str, slurm_param
    paths = [p.replace("\\", "/") for p in _paths(data_paths)]
    psf_list = [p.replace("\\", "/") for p in _paths(psfs)]
    patterns = _patterns(channel_patterns)
    fol_str = "DSR" if rotate else "DS"
    paths_ds = [f"{p}/{fol_str}" for p in paths] if deskew and recon else list(paths)

    deskew_kw = {
        "skew_angle": _f10(skew_angle),
        "reverse": bool(reverse),
        "nphases": _f10(nphases),
        "rotate": bool(rotate),
        "save_16bit": bool(save_16bit),
    }
    recon_kw = {
        "islattice": bool(islattice),
        "na_det": _f10(na_det),
        "na_ext": _f10(na_ext),
        "nimm": _f10(nimm),
        "wvl_em": _f10(wvl_em),
        "wvl_ext": _f10(wvl_ext),
        "w": _f10(w),
        "apodize": bool(apodize),
        "norientations": _f10(norientations),
        "lattice_period": _f10(lattice_period),
        "lattice_angle": _f10(lattice_angle),
        "phase_step": _f10(phase_step),
        "pxl_dim_data": tuple(_f10(v) for v in pxl_dim_data),
        "pxl_dim_psf": tuple(_f10(v) for v in pxl_dim_psf),
        "normalize_orientations": bool(normalize_orientations),
        "norders": _f10(norders),
        "nphases": _f10(nphases),
        "overlap": _f10(overlap),
        "chunk_size": tuple(_f10(v) for v in chunk_size),
        "max_trial_num": _f10(max_trial_num),
        "unit_wait_time": _f10(unit_wait_time),
        "int_thresh": _f10(int_thresh),
        "occ_thresh": _f10(occ_thresh),
        "edge_taper": bool(edge_taper),
        "edge_taper_val": _f10(edge_taper_val),
        "perdecomp": bool(perdecomp),
        "use_gpu": bool(use_gpu),
        "save_16bit": bool(save_16bit),
        "gpu_precision": gpu_precision,
        "ds": bool(ds),
        "background": _f10(background),
        "edge_erosion": _f10(edge_erosion),
        "erode_before": bool(erode_before),
        "erode_after": bool(erode_after),
        "erode_maskfile": erode_maskfile,
        "save_maskfile": bool(save_maskfile),
        "results_dir_name": results_dir_name,
    }

    workers = int(par_pool_size)
    cpus = available_cpus()
    if workers > cpus:
        print(
            sprintf(
                "Could not allocate %d workers. Allocating %d workers instead (max allowed by "
                "current system).\n",
                workers,
                cpus,
            ),
            end="",
            flush=True,
        )
        workers = cpus
    pool = ThreadPoolExecutor(max_workers=max(1, workers))
    jobs: list[tuple[Any, str]] = []
    submitted: set[str] = set()
    first = True
    latest = 0.0

    def stage_files(folder: str, pattern: str, newest_only: bool) -> list[str] | None:
        entries = matlab_dir(f"{folder}/*{pattern}*.tif")
        if not entries:
            return None
        full = [f"{e.folder}/{e.name}" for e in entries]
        if streaming:
            ages = _minutes_old(full)
            if min(ages) < unit_wait_time:
                if newest_only:
                    full = [f for f, a in zip(full, ages, strict=True) if a != min(ages)]
                else:
                    full = [f for f, a in zip(full, ages, strict=True) if a >= unit_wait_time]
        return [f for f in full if f not in submitted]

    def submit(label: str, inputs, outputs, tasks, *, gpu: bool, folder: str, pattern: str):
        print(
            sprintf(
                "Attempting %s on %d file(s) for pattern '%s' in folder '%s'\n",
                label,
                len(inputs),
                pattern,
                folder,
            ),  # fmt: skip
            end="",
            flush=True,
        )
        fut = pool.submit(
            _run, inputs, outputs, tasks, use_gpu=gpu, cpus_per_task=cpus_per_task,
            task_batch_num=recon_batch_num if gpu else 1, master_compute=master_compute,
            parse_cluster=parse_cluster, job_log_dir=job_log_dir,
        )  # fmt: skip
        msg = sprintf("Finished %s on %d file(s) for pattern '%s' in folder '%s'\n", label,
                      len(inputs), pattern, folder)  # fmt: skip
        jobs.append((fut, msg))

    try:
        while first or jobs or (streaming and latest < max_modify_time):
            if first:
                first = False
            else:
                time.sleep(_PAUSE)
            if deskew:
                for i, d in enumerate(paths):
                    for patt in patterns:
                        full = stage_files(d, patt, True)
                        if full is None:
                            print(
                                sprintf(
                                    "No Files found for channel pattern: '%s'. Skipping to "
                                    "next pattern.\n",
                                    patt,
                                ),
                                end="",
                                flush=True,
                            )
                            continue
                        if not full:
                            continue
                        submitted.update(full)
                        if parse_cluster:
                            job_log_dir = _job_log_dir(job_log_dir, paths[i], group=True)
                        inputs, outputs, tasks = [], [], []
                        for f in full:
                            os.makedirs(f"{d}/{fol_str}", exist_ok=True)
                            _, fsname, ext = fileparts(f)
                            inputs.append(f"{d}/{os.path.basename(f)}")
                            outputs.append(f"{d}/{fol_str}/{fsname}{ext}")
                            tasks.append(
                                Task.of(
                                    deskew_phases_file,
                                    inputs[-1],
                                    _f10(xy_pixel_size),
                                    _f10(dz),
                                    **deskew_kw,
                                )  # fmt: skip
                            )
                        if all(os.path.exists(o) for o in outputs):
                            print(
                                sprintf(
                                    "Skipping Deskew for pattern '%s' in folder '%s' because it "
                                    "already exists.\n",
                                    patt,
                                    d,
                                ),
                                end="",
                                flush=True,
                            )
                            continue
                        submit("Deskew", inputs, outputs, tasks, gpu=False, folder=d, pattern=patt)
            if recon:
                for i, d in enumerate(paths_ds):
                    for c, patt in enumerate(patterns):
                        full = stage_files(d, patt, False)
                        if full is None:
                            if not streaming:
                                print(
                                    sprintf(
                                        "No Files found for channel pattern: '%s'. Skipping to "
                                        "next pattern.\n",
                                        patt,
                                    ),
                                    end="",
                                    flush=True,
                                )
                            continue
                        if not full:
                            continue
                        submitted.update(full)
                        if parse_cluster:
                            job_log_dir = _job_log_dir(job_log_dir, paths_ds[i], group=True)
                        inputs, outputs, tasks = [], [], []
                        for f in full:
                            _, fsname, ext = fileparts(f)
                            inputs.append(f"{d}/{os.path.basename(f)}")
                            outputs.append(f"{d}/{results_dir_name}/{fsname}_recon{ext}")
                            os.makedirs(f"{d}/{results_dir_name}", exist_ok=True)
                            psf = psf_list[c] if c < len(psf_list) else ""
                            tasks.append(Task.of(sim_recon_file, inputs[-1], psf, **recon_kw))
                        if all(os.path.exists(o) for o in outputs):
                            print(
                                sprintf(
                                    "Skipping Recon for pattern '%s' in folder '%s' because it "
                                    "already exists.\n",
                                    patt,
                                    d,
                                ),
                                end="",
                                flush=True,
                            )
                            continue
                        submit("Recon", inputs, outputs, tasks, gpu=bool(use_gpu), folder=d,
                               pattern=patt)  # fmt: skip
            if streaming:
                latest = math.inf
                for folders in (paths, paths_ds if deskew and recon else []):
                    for d in folders:
                        if not os.path.isdir(d):
                            continue
                        for patt in patterns:
                            files = [
                                f"{e.folder}/{e.name}" for e in matlab_dir(f"{d}/*{patt}*.tif")
                            ]
                            if files:
                                latest = min(latest, min(_minutes_old(files)))
            still = []
            for fut, msg in jobs:
                if fut.done():
                    print(msg, end="", flush=True)
                else:
                    still.append((fut, msg))
            jobs[:] = still
    finally:
        pool.shutdown(wait=True)
