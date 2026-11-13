"""Zarr resaving with other chunks, a crop box or as the port's native layout:
``XR_resave_zarr_wrapper``, ``XR_resaveSingleZarr`` and ``resaveZarrBlock``.

Arrays are ``(z, y, x)``; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive; block
and batch sizes ``(z, y, x)``. The MATLAB-name layer converts MATLAB's values. Inputs of either
layout are read (MATLAB's ``[ny, nx, nz]`` arrays come back as ``(z, y, x)``), and results are
native zarr arrays (``[nz, ny, nx]``, order C), so resaving also turns MATLAB's arrays into the
port's.

``resave_zarr`` cuts the output into batches of whole chunks (MATLAB's
``XR_zarrChunkCoordinatesExtraction``, batches numbered y fastest, then x, then z), hands groups
of batches to ``resave_zarr_block`` tasks through the job framework (each task's flag marks it
done), and moves the temporary result ``<result>_<uuid>`` into place when all are done; a
result another writer published is kept. Without a cluster, neighboring batches are read and
written together (up to one chunk per thread, within 1/16 of this process's share of the
memory), as MATLAB does, and their flags name the group size. Flags keep MATLAB's names and
values (1-based batch numbers, ``[y x z]`` sizes).
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = ["ResaveError", "resave_zarr", "resave_zarr_block", "resave_zarr_dataset"]

# GB a batch task asks for: this many times its batch (MATLAB's)
BATCH_MEMORY_FACTOR = 1.5
# GB a file's task asks for in the dataset wrapper: this many times its batch or file
FILE_MEMORY_FACTOR = 1.9


class ResaveError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _ints(v: Any) -> list[int]:
    return [int(x) for x in np.asarray(v, dtype=np.float64).ravel()]


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    return np.size(v) == 0


def _print(msg: str, end: str = "\n") -> None:
    print(msg, end=end, flush=True)


# ---------------------------------------------------------------------------------------------
# resaveZarrBlock
# ---------------------------------------------------------------------------------------------


@register(
    "resaveZarrBlock",
    Param("batchInds", "@isnumeric", "Batch numbers of the task (0-based).", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr (either layout)."),
    Param("resultFullpath", "@(x) ischar(x)", "Output zarr (native, existing)."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param("BatchBBoxes", "@isnumeric", "Input box of each batch.", convert="bboxes"),
    Param("RegionBBoxes", "@isnumeric", "Output box of each batch.", convert="bboxes"),
    Param("overwrite", "@islogical", "Redo batches whose flag exists."),
    task="T01",
)
def resave_zarr_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    result_fullpath: str,
    flag_fullname: str,
    batch_b_boxes: Sequence[Sequence[int]],
    region_b_boxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
) -> None:
    """Copy the input boxes ``batch_b_boxes`` of a zarr array into the output boxes
    ``region_b_boxes`` of an existing native array (``resaveZarrBlock``), then write the
    task's flag. The flag's folder must exist (else MATLAB's error: the run is over)."""
    from petakit5d.io import read_zarr, write_zarr

    inds = _ints(batch_inds)
    flag_fullname = os.fspath(flag_fullname)
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise ResaveError(
            f"The batch flag directory {flag_path} does not exist, skip the processing!"
        )
    if os.path.lexists(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            _print(f"The batch files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!")
            return
    if not os.path.isdir(zarr_fullpath):
        raise ResaveError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    if not os.path.isdir(result_fullpath):
        raise ResaveError(f"The output zarr file {result_fullpath} doesnot exist!")

    inb = [_ints(b) for b in batch_b_boxes]
    outb = [_ints(b) for b in region_b_boxes]
    for i, bi in enumerate(inds):
        _print(f"Process Batch {bi + 1}... ", end="")
        t0 = time.perf_counter()
        data = read_zarr(zarr_fullpath, input_bbox=inb[i])
        write_zarr(data, result_fullpath, bbox=outb[i], create=False)
        del data
        _print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.")
    with open(flag_fullname, "w"):
        pass


# ---------------------------------------------------------------------------------------------
# XR_resaveSingleZarr
# ---------------------------------------------------------------------------------------------


def _comp_threads() -> int:
    """MATLAB's maxNumCompThreads: the CPUs this process may use"""
    from petakit5d.core.threads import available_cpus

    return available_cpus()


def _node_cores() -> float:
    """SLURM_CPUS_ON_NODE, else the machine's cores (MATLAB's feature('numcores'))"""
    from petakit5d.core.text import str2double

    n = str2double(os.environ.get("SLURM_CPUS_ON_NODE", ""))
    if math.isnan(n):
        n = float(os.cpu_count() or 1)
    return n


def _memory_gb() -> float:
    """getSystemMemory(), at most the memory available now (Linux: MemAvailable); NaN unknown"""
    import sys

    from petakit5d.core.memory import get_system_memory

    mem = get_system_memory()
    if sys.platform.startswith("linux"):
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        avail = float(line.split()[1]) / 2**20
                        mem = mem if math.isnan(mem) else min(mem, avail)
                        break
        except OSError:
            pass
    return mem


def local_batch_size(
    out_size: Sequence[int], batch_size: Sequence[int], block_size: Sequence[int], byte_num: int
) -> list[int]:
    """The batch size for resaving without a cluster (MATLAB's ``local_batch_size``):
    neighboring batches together (along y, then x, then z), up to one chunk per thread of this
    process, with their data within 1/16 of this process's share of the memory (the memory
    times its share of the node's CPUs). Sizes are ``(z, y, x)``."""
    mem = _memory_gb()
    threads = _comp_threads()
    share = min(1.0, threads / _node_cores())
    per_batch = math.prod(b / k for b, k in zip(batch_size, block_size, strict=True))
    a = math.floor(threads / per_batch)
    if math.isnan(mem):
        return list(batch_size)
    b = math.floor(mem * share * 2**30 / 16 / (math.prod(batch_size) * byte_num))
    max_batch_num = min(a, b)
    if not max_batch_num > 1:
        return list(batch_size)
    sub = [math.ceil(o / b) for o, b in zip(out_size, batch_size, strict=True)]
    n = [1, 1, 1]
    for d in (1, 2, 0):  # y, x, z
        n[d] = min(sub[d], math.floor(max_batch_num / math.prod(n)))
    return [b * k for b, k in zip(batch_size, n, strict=True)]


@register(
    "XR_resaveSingleZarr",
    Param("zarrFullpath", "@ischar", "Input zarr (either layout)."),
    Param("resultFullpath", "@ischar", "Output zarr (native)."),
    Param(
        "inputBbox",
        "@isnumeric",
        "Region (z0, y0, x0, z1, y1, x1) to resave (None: all).",
        convert="bbox",
    ),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the output.", convert="size_yxz"),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of a block task.", convert="size_yxz"),
    Param("shardSize", "@isnumeric", "Not used (reserved).", convert="size_yxz"),
    Param("parseCluster", "@islogical", "Run the batches as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes batches too."),
    Param("cpusPerTask", "@isscalar", "CPUs of a batch job."),
    Param("uuid", "@ischar", "Id of the temporary result and flags (empty: a new one)."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T01",
)
def resave_zarr(
    zarr_fullpath: str,
    result_fullpath: str,
    *,
    input_bbox: Sequence[int] | None = None,
    block_size: Sequence[int] = (256, 256, 256),
    batch_size: Sequence[int] = (512, 512, 512),
    shard_size: Sequence[int] | None = (512, 512, 512),
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 1,
    uuid: str = "",
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Resave a zarr array (or its box ``input_bbox``) as a native array with chunks
    ``block_size`` (``XR_resaveSingleZarr``), in batches of whole chunks through the job
    framework. An existing result is kept."""
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.axes import to_yxz
    from petakit5d.core.blocks import zarr_chunk_coordinates
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive
    from petakit5d.io import create_zarr, get_image_dtype, get_zarr_info
    from petakit5d.jobs import framework
    from petakit5d.jobs.publish import publish_result, remove_path
    from petakit5d.jobs.tasks import Task

    del shard_size  # (reserved in MATLAB)
    zarr_fullpath = os.fspath(zarr_fullpath)
    result_fullpath = os.fspath(result_fullpath)
    if not uuid:
        uuid = get_uuid()
    if os.path.isdir(result_fullpath):
        _print("The output result exists, skip it!")
        return
    result_path, result_fsname, _ = fileparts(result_fullpath)
    tmp_path = f"{result_fullpath.removesuffix('/')}_{uuid}"
    if not os.path.isdir(zarr_fullpath):
        raise ResaveError(f"The input file {zarr_fullpath} does not exist, pleast double check!")

    info = get_zarr_info(zarr_fullpath)
    dtype = info.dtype
    sz = [int(v) for v in info.size]
    out_size = list(sz)
    start = [0, 0, 0]
    if not _empty(input_bbox):
        b = _ints(input_bbox)
        out_size = [b[i + 3] - b[i] for i in range(3)]
        start = b[:3]

    block = [min(o, b) for o, b in zip(out_size, _ints(block_size), strict=True)]
    batch = [min(o, b) for o, b in zip(out_size, _ints(batch_size), strict=True)]
    batch = [math.ceil(b / k) * k for b, k in zip(batch, block, strict=True)]
    flag_name = "batches"
    if not parse_cluster:
        local = local_batch_size(out_size, batch, block, data_type_to_byte_number(dtype))
        if local != batch:
            flag_name = "batches_{}_{}_{}".format(*[int(v) for v in to_yxz(local)])
        batch = local
    num_batch = math.prod(math.ceil(o / b) for o, b in zip(out_size, batch, strict=True))

    batch_boxes, region_boxes, _ = zarr_chunk_coordinates(
        out_size, batch_size=batch, block_size=block, same_batch_size=False, border_size=0
    )
    batch_boxes = np.asarray(batch_boxes, dtype=np.int64) + np.array(start + start)

    _print("Process resaving zarr data...")
    flag_path = f"{result_path}/zarr_flag/{result_fsname}_{uuid}/"
    if not os.path.isdir(flag_path):
        mkdir_recursive(flag_path)

    if os.path.isdir(tmp_path):
        tmp_info = get_zarr_info(tmp_path)
        if list(tmp_info.block_size) != block or list(tmp_info.size) != out_size:
            remove_path(tmp_path)
            remove_path(flag_path)
            os.makedirs(flag_path, exist_ok=True)
    if not os.path.isdir(tmp_path):
        nchunks = math.prod(math.ceil(o / b) for o, b in zip(out_size, block, strict=True))
        create_zarr(
            tmp_path,
            data_size=out_size,
            block_size=block,
            dtype=dtype,
            dim_separator="/" if nchunks > 10000 else ".",
        )

    task_size = max(5, math.ceil(num_batch / 5000))
    num_tasks = math.ceil(num_batch / task_size)
    tasks, outputs = [], []
    for t in range(num_tasks):
        inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
        flag = f"{flag_path}/{flag_name}_{inds[0] + 1}_{inds[-1] + 1}.mat"
        outputs.append(flag)
        tasks.append(
            Task.of(
                resave_zarr_block,
                inds,
                zarr_fullpath,
                tmp_path,
                flag,
                batch_boxes[inds].tolist(),
                np.asarray(region_boxes)[inds].tolist(),
            )
        )
    inputs = [zarr_fullpath] * num_tasks

    byte = data_type_to_byte_number(get_image_dtype(zarr_fullpath, True))
    mem_allocate = math.prod(batch) * byte / 2**30 * BATCH_MEMORY_FACTOR
    common = {
        "final_out_fullpath": result_fullpath,
        "cpus_per_task": cpus_per_task,
        "mem_allocate": mem_allocate,
        "parse_cluster": parse_cluster,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = framework.run_tasks(inputs, outputs, tasks, **common)
    if not all(done):
        done = framework.run_tasks(inputs, outputs, tasks, master_compute=master_compute, **common)
    if not all(done):
        raise ResaveError(f"{sum(not d for d in done)} / {len(done)} tasks failed!")

    # another copy of this task may have published the result: it is kept (the job framework
    # then reports all batches as done, so this copy's temporary result can be incomplete)
    publish_result(tmp_path, result_fullpath, flag_path)


# ---------------------------------------------------------------------------------------------
# XR_resave_zarr_wrapper
# ---------------------------------------------------------------------------------------------


def _save_parameters(result_path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) and parameters.json of the call (MATLAB's ``pr``)"""
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(resave_zarr_dataset, values).items()}
    write_mat(result_path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, result_path + "/parameters.json")


@register(
    "XR_resave_zarr_wrapper",
    Param("dataPaths", "@(x) iscell(x) || ischar(x)", "Dataset folders (one or several)."),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of a block task.", convert="size_yxz"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the output.", convert="size_yxz"),
    Param("shardSize", "@isnumeric", "Not used (reserved).", convert="size_yxz"),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to resave."),
    Param(
        "inputBbox",
        "@isnumeric",
        "Region (z0, y0, x0, z1, y1, x1) to resave (None: all).",
        convert="bbox",
    ),
    Param("parseCluster", "@islogical", "Run the files and batches as Slurm jobs."),
    Param("largeFile", "@islogical", "Resave in batches (else whole files when they are few)."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param("uuid", "@ischar", "Id of the temporary results (empty: a new one)."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T01",
)
def resave_zarr_dataset(
    data_paths: str | Sequence[str],
    *,
    result_dir_name: str = "resaved_zarr",
    batch_size: Sequence[int] = (512, 512, 512),
    block_size: Sequence[int] = (256, 256, 256),
    shard_size: Sequence[int] | None = None,
    channel_patterns: Sequence[str] = (".zarr",),
    input_bbox: Sequence[int] | None = None,
    parse_cluster: bool = True,
    large_file: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 1,
    uuid: str = "",
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Resave the zarr files of datasets that match a channel pattern
    (``XR_resave_zarr_wrapper``) as ``<data path>/<result_dir_name>/<name>.zarr``, one
    ``resave_zarr`` task per file through the job framework. parameters.mat and .json of the
    call go into each result folder."""
    values = {k: v for k, v in locals().items() if k in _WRAPPER_ARGS}
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import simplify_path
    from petakit5d.core.text import through_text
    from petakit5d.io import get_image_dtype, get_image_size
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_image_filenames

    if isinstance(data_paths, (str, os.PathLike)):
        data_paths = [data_paths]
    paths = [os.fspath(p) for p in data_paths]
    result_paths = []
    data_path = ""
    for d, p in enumerate(paths):
        data_path = p
        if not p.endswith(os.sep):
            p = p + os.sep
        paths[d] = simplify_path(p) + "/"
        result_path = f"{paths[d]}/{result_dir_name}/"
        result_paths.append(result_path)
        os.makedirs(result_path, exist_ok=True)
        group_writable(result_path)
        _save_parameters(result_path, values)

    if not uuid:
        uuid = get_uuid()
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(data_path, job_log_dir).parse_cluster)

    _, fsns, fd_inds, filepaths, _ = parse_image_filenames(paths, True, channel_patterns)
    n_f = len(fsns)
    sz = [int(v) for v in get_image_size(filepaths[0])[:3]]
    byte = data_type_to_byte_number(get_image_dtype(filepaths[0]))

    batch = _ints(batch_size)
    block = _ints(block_size)
    if not large_file:
        # few batches: whole files (the batch is the output)
        out_size = list(sz)
        if not _empty(input_bbox):
            b = _ints(input_bbox)
            out_size = [b[i + 3] - b[i] for i in range(3)]
        lim = [max(512, max(a, k)) for a, k in zip(batch, block, strict=True)]
        if math.prod(math.ceil(o / m) for o, m in zip(out_size, lim, strict=True)) > 10:
            large_file = True
        else:
            batch = out_size

    results = [f"{result_paths[int(fd_inds[i])]}/{fsns[i]}.zarr" for i in range(n_f)]
    cpus = through_text(cpus_per_task, "%d")
    cpus = int(cpus) if float(cpus) == int(cpus) else float(cpus)
    tasks = [
        Task.of(
            resave_zarr,
            filepaths[i],
            results[i],
            input_bbox=None if _empty(input_bbox) else _ints(input_bbox),
            block_size=block,
            batch_size=batch,
            shard_size=None if _empty(shard_size) else _ints(shard_size),
            parse_cluster=bool(parse_cluster),
            master_compute=bool(master_compute),
            cpus_per_task=cpus,
            uuid=uuid,
            mcc_mode=bool(mcc_mode),
            config_file=str(config_file),
        )
        for i in range(n_f)
    ]
    voxels = math.prod(batch) if large_file else math.prod(sz)
    mem_allocate = voxels * byte / 2**30 * FILE_MEMORY_FACTOR
    framework.run_tasks(
        list(filepaths),
        results,
        tasks,
        parse_cluster=parse_cluster,
        job_log_dir=job_log_dir,
        master_compute=master_compute,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )


_WRAPPER_ARGS = (
    "data_paths",
    "result_dir_name",
    "batch_size",
    "block_size",
    "shard_size",
    "channel_patterns",
    "input_bbox",
    "parse_cluster",
    "large_file",
    "master_compute",
    "job_log_dir",
    "cpus_per_task",
    "uuid",
    "mcc_mode",
    "config_file",
)
