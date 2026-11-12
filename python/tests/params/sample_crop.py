"""The worked example of a spec: XR_crop_dataset as the tools area (T03) declares it.

It lives in a registry of its own, so that the package registry gets the real one from T03.
The function body only records its call, so tests can compare the three ways in.
"""

from __future__ import annotations

from typing import Any

from petakit5d.params import Param, Registry

REG = Registry()
CALLS: list[dict[str, Any]] = []


@REG.register(
    "XR_crop_dataset",
    Param(
        "dataPaths",
        "@(x) ischar(x) || iscell(x)",
        "Directory paths for the datasets: one path or a cell array of paths with the "
        "same experimental settings.",
    ),
    Param(
        "inputBbox",
        "@isnumeric",
        "Bounding box of the crop: (z0, y0, x0, z1, y1, x1) (MATLAB: [ymin xmin zmin ymax xmax "
        "zmax]).",
        convert="bbox",
    ),
    Param("resultDirName", "@ischar", "Crop result directory under each data path."),
    Param(
        "cropType",
        "@ischar",
        "'fixed', 'moving' (moves linearly to lastStartCoords at the last time point) or "
        "'center' (the centre of the image with the size of inputBbox).",
    ),
    Param("pad", "@islogical", "Pad an empty region where the box leaves the image."),
    Param(
        "lastStartCoords",
        "@isnumeric",
        "Start coordinates (z, y, x) of the box at the last time point, for cropType 'moving'.",
        convert="point1_yxz",
    ),
    Param("channelPatterns", "@iscell", "Channel identifiers of the included channels."),
    Param("zarrFile", "@islogical", "Zarr files as input."),
    Param(
        "largeFile",
        "@islogical",
        "Large-scale crop that reads and writes the cropped regions in place (zarr only).",
    ),
    Param("saveZarr", "@islogical", "Save the results as zarr files."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) per task.", convert="size_yxz"),
    Param(
        "blockSize",
        "@isnumeric",
        "Block (chunk) size (z, y, x) of zarr output.",
        convert="size_yxz",
    ),
    Param("save16bit", "@islogical", "Save 16-bit results."),
    Param("parseCluster", "@islogical", "Use the Slurm cluster."),
    Param("masterCompute", "@islogical", "The master job takes part in the processing."),
    Param("jobLogDir", "@ischar", "Folder of the Slurm job logs."),
    Param("cpusPerTask", "@isnumeric", "CPU cores per task for Slurm jobs."),
    Param("uuid", "@ischar", "uuid in the temporary result paths (empty: a new one)."),
    Param("mccMode", "@islogical", "Run tasks with the compiled master."),
    Param("configFile", "@ischar", "Config file for job submission."),
    task="T03",
)
def crop_dataset(
    data_paths,
    input_bbox,
    *,
    result_dir_name="Cropped",
    crop_type="fixed",
    pad=False,
    last_start_coords=None,
    channel_patterns=("CamA_ch0", "CamB_ch0"),
    zarr_file=False,
    large_file=False,
    save_zarr=False,
    batch_size=(1024, 1024, 1024),
    block_size=(256, 256, 256),
    save_16bit=True,
    parse_cluster=True,
    master_compute=True,
    job_log_dir="../job_logs",
    cpus_per_task=2,
    uuid="",
    mcc_mode=False,
    config_file="",
):
    """Crop every file of the datasets with one bounding box."""
    args = dict(locals())
    CALLS.append(args)
    return args
