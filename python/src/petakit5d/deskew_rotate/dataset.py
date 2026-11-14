"""Deskew and rotation of datasets (``XR_deskew_rotate_data_wrapper``).

``deskew_rotate_dataset`` deskews and/or rotates every file of one or several dataset folders
whose path matches a channel pattern, one task per file through the job framework (locally or
as Slurm jobs), into ``<data>/DS/`` and ``<data>/DSR/``:

- regular files (tiff or zarr): ``deskew_rotate_file`` tasks, with the file's partial volume
  files (``<name>_part0001.tif``, ...) read as one volume, the z flip from its setting file
  (``parse_setting_file``), the flat-field images of its channel (``ff_correction``: the
  first channel pattern the path contains) and ``crop``;
- large zarr files (``large_file``): ``deskew_rotate_zarr`` tasks, each running its blocks as
  tasks of their own (the combined path only, zarr results, sample scan only), with the file's
  z flip.

The files are listed as the pipeline lists them (``pipeline.filenames``: by listing each folder,
not one file at a time). Existing results are skipped by the tasks; ``overwrite`` removes the
result folders first. ``parameters.mat`` and ``parameters.json`` of the call (MATLAB's values)
go into the DSR folder (the DS folder for deskew only).

The jobs ask for the memory of the largest file (all its partial files, after the input box)
from ``estimate_deskew_rotate_memory``, which sets their CPUs; the job of a large file asks for
what it holds (the process and the merge of the MIP slabs, and a block's estimate if it
computes blocks itself) and goes after its block jobs (``coordinator_nice``). Each round of
tries doubles the CPUs and the memory, three rounds of two tries, as MATLAB does.

Values reach the tasks as MATLAB's task strings pass them: ``constOffset`` with 10 decimals,
the skew angle of large files and their CPUs with ``%d`` (7 significant digits), vectors with
``mat2str``'s 15 digits; the other numbers exactly.

Arrays and boxes are the port's: ``input_bbox`` ``(z0, y0, x0, z1, y1, x1)``, sizes ``(z, y,
x)``, ``resample_factor`` ``(z, y, x)``, axis orders ``'zyx'``.

MATLAB quirks kept (each has a test):

- with ``objective_scan`` the result folder becomes ``rotated/`` only when it is given as
  ``'DSR/'`` (the default ``'DSR'`` stays);
- the flat field of a file is the one of the first channel pattern its path contains; a file
  that contains none raises MATLAB's error;
- the parameters files hold the call's values as given (before the changes for objective scans
  and large files).
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = ["DeskewRotateDatasetError", "controller_resources", "deskew_rotate_dataset"]

_INTERP = "@(x) any(strcmpi(x, {'cubic', 'linear'}))"


class DeskewRotateDatasetError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _logical(x: Any) -> bool:
    return bool(np.asarray(x).reshape(-1)[0]) if np.size(x) else False


def _first(x: Any) -> float:
    return float(np.asarray(x, dtype=np.float64).reshape(-1)[0])


def _empty(x: Any) -> bool:
    return x is None or np.size(x) == 0


def _through(value: float, fmt: str) -> float:
    from petakit5d.core.text import through_text

    return float(through_text(float(value), fmt))


def _through_mat2str(value: Any) -> Any:
    from petakit5d.core.text import through_mat2str

    if _empty(value):
        return None
    v = np.asarray(value, dtype=np.float64).reshape(-1)
    out = tuple(float(x) for x in np.asarray(through_mat2str(v)).reshape(-1))
    return out[0] if np.ndim(value) == 0 else out


def _first_selected(cell: Sequence[str], hit: list[int]) -> str:
    """``x = cell{mask}`` with a logical mask whose true places are ``hit``: the first of them
    (a comma-separated list gives its first value), MATLAB's errors for none and for a place
    past the end"""
    if any(k >= len(cell) for k in hit):
        raise DeskewRotateDatasetError(
            f"Index exceeds the number of array elements. Index must not exceed {len(cell)}.",
            "MATLAB:badsubscript",
        )
    if not hit:
        raise DeskewRotateDatasetError(
            "Insufficient number of outputs from right hand side of equal sign to satisfy "
            "assignment.",
            "MATLAB:needMoreRhsOutputs",
        )
    return cell[hit[0]]


def _ints(v: Any) -> list[int]:
    return [int(x) for x in np.asarray(v).reshape(-1)]


def _box_size(box: Sequence[int]) -> tuple[int, int, int]:
    b = _ints(box)
    return (b[3] - b[0], b[4] - b[1], b[5] - b[2])


def _save_parameters(path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) and parameters.json of the call (MATLAB's ``pr``)"""
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(deskew_rotate_dataset, values).items()}
    write_mat(path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, path + "/parameters.json")


def _make_folder(path: str) -> None:
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
        if os.name != "nt":
            from petakit5d.io.tiff import group_writable

            group_writable(path)


def controller_resources(
    in_size: Sequence[int],
    xy_pixel_size: float,
    dz: float,
    skew_angle: float,
    save_16bit: bool,
    interp_method: str,
    resample_factor: Sequence[float] | None,
    input_bbox: Sequence[int] | None,
    batch_size: Sequence[int],
    block_size: Sequence[int],
    save_mip: bool,
    data_type: str,
    computes: bool,
    cpus_per_task: float,
) -> tuple[float, float]:
    """The request ``(cpus_per_task, mem_allocate)`` of the job of a large file
    (``deskew_rotate_zarr``, the wrapper's ``dsr_controller_resources``), which runs its y
    batches as jobs (``coordinator_resources``): the process and the merge of the MIP slabs,
    and if it computes batches itself (``computes``) the estimate of its largest batch on the
    CPUs of a batch's job. The geometry is ``deskew_rotate_zarr``'s, without MIP masks (which
    only make the result smaller), its resample factor completed to ``[x y z]`` as there.
    ``in_size`` is the input's ``(nz, ny, nx)``. (MATLAB also counts a batch's flipped copy
    for ``flip_zstack``; the port's blocks flip in place.)"""
    from petakit5d.deskew_rotate.geometry import _cos, _round, _sin
    from petakit5d.deskew_rotate.memory import estimate_deskew_rotate_memory
    from petakit5d.deskew_rotate.zarr import _resample_matlab, _rs_xyz
    from petakit5d.jobs.coordinator import coordinator_resources
    from petakit5d.tools.mip import merge_memory

    nz, ny, nx = (int(v) for v in in_size)
    if not _empty(input_bbox):
        nz, ny, nx = _box_size(input_bbox)
    # 1, 2 or 3 values in xyz order, as deskew_rotate_zarr takes them
    rx, ry, rz = _rs_xyz(_resample_matlab(resample_factor))[:3]
    rs = [ry, rx, rz]  # rs([2, 1, 3]): [y x z]
    theta = float(skew_angle) * math.pi / 180
    st = _sin(abs(theta))
    z_aniso = st * float(dz) / float(xy_pixel_size)
    o = [
        _round(ny),
        _round((nx - 1) * _cos(theta) + (nz - 1) * z_aniso / st),
        _round((nx - 1) * st - 4),
    ]
    out = [int(max(1.0, _round(o[k] / rs[k]))) for k in range(3)]  # [y x z]
    bs = _ints(batch_size)  # (z, y, x)
    blk = _ints(block_size)
    batch_out = [min(bs[1], out[0]), out[1], out[2]]
    blk_yxz = [blk[1], blk[2], blk[0]]
    in_yxz = [ny, nx, nz]
    out_block = [min(blk_yxz[k], in_yxz[k], batch_out[k]) for k in range(3)]
    batch_in = [int(min(_round((batch_out[0] + 2) * ry), ny)), nx, nz]
    zyx = lambda v: (int(v[2]), int(v[0]), int(v[1]))  # noqa: E731
    block_mem = estimate_deskew_rotate_memory(
        zyx(batch_in),
        xy_pixel_size,
        dz,
        skew_angle=skew_angle,
        save_16bit=save_16bit,
        interp_method=interp_method,
        resample_factor=resample_factor,
        out_size=zyx(batch_out),
        large_file=True,
        block_size=zyx(out_block),
        data_type=data_type,
    )[0]
    merge_mem = 0.0
    if save_mip:
        dtype = "uint16" if save_16bit else "single"
        merge_mem = merge_memory(zyx(out), zyx(batch_out), zyx(out_block), dtype)
    return coordinator_resources(merge_mem, computes, cpus_per_task, block_mem)


@register(
    "XR_deskew_rotate_data_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param("DSDirName", "@ischar", "Folder of the deskewed results under each data path."),
    Param("DSRDirName", "@ischar", "Folder of the deskewed and rotated results."),
    Param("deskew", "@islogical", "Deskew the data."),
    Param("rotate", "@islogical", "Rotate the deskewed data."),
    Param("overwrite", "@islogical", "Remove the result folders first."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to process."),
    Param("dz", "@isscalar", "Scan interval in um."),
    Param("xyPixelSize", "@(x) isvector(x) && numel(x) <= 2", "Pixel size in um."),
    Param("skewAngle", "@isscalar", "Skew angle in degrees."),
    Param("objectiveScan", "@islogical", "Objective scan."),
    Param("zStageScan", "@islogical", "Z stage scan (orthogonal to objective scan)."),
    Param("reverse", "@islogical", "Reverse scan direction."),
    Param("flipZstack", "@islogical", "Flip the z stacks."),
    Param("parseSettingFile", "@islogical", "Flip each file's z stack as its setting file says."),
    Param("crop", "@islogical", "Crop the empty triangles of the separate deskew."),
    Param("inputAxisOrder", "@ischar", "Axis order of the input.", convert="axis_string"),
    Param("outputAxisOrder", "@ischar", "Axis order of the output.", convert="axis_string"),
    Param("DSRCombined", "@(x) islogical(x)", "Combined deskew, rotation and resampling."),
    Param("FFCorrection", "@islogical", "Flat-field correction (regular files)."),
    Param("BKRemoval", "@islogical", "Camera background removal (known issue: MATLAB errors)."),
    Param("lowerLimit", "@isnumeric", "Lowest value of the normalized flat field."),
    Param(
        "constOffset",
        "@(x) isnumeric(x)",
        "If set, a constant background instead of the camera's.",
        convert="same",
    ),
    Param("FFImagePaths", "@iscell", "Flat-field images, one per channel pattern."),
    Param("backgroundPaths", "@iscell", "Background images, one per channel pattern."),
    Param("save16bit", "@islogical", "Save 16-bit results (else single).", convert="same"),
    Param("save3DStack", "@islogical", "Save the volumes (else empty result files)."),
    Param("saveMIP", "@islogical", "Save the MIPs of the results."),
    Param("largeFile", "@islogical", "Large zarr files: deskew and rotate in blocks."),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("saveZarr", "@islogical", "Save the results as zarr."),
    Param("batchSize", "@isvector", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("blockSize", "@isvector", "Chunk size (z, y, x) of zarr results.", convert="size_yxz"),
    Param(
        "inputBbox",
        "@(x) isempty(x) || isvector(x)",
        "Box (z0, y0, x0, z1, y1, x1) of the input.",
        convert="bbox",
    ),
    Param("taskSize", "@isnumeric", "Batches per task of large files.", convert="same"),
    Param(
        "resampleFactor",
        "@(x) isempty(x) || isnumeric(x)",
        "Resampling after rotation (z, y, x).",
        convert="resample_xyz",
    ),
    Param("interpMethod", _INTERP, "Interpolation: 'linear' or 'cubic'."),
    Param("maskFullpaths", "@iscell", "MIP masks (xy, xz, yz) that set the boxes of large files."),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("parseParfor", "@islogical", "Local process pool of large files' blocks."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    keep_unmatched=True,
    task="D04",
)
def deskew_rotate_dataset(
    data_paths: str | Sequence[str],
    *,
    ds_dir_name: str = "DS",
    dsr_dir_name: str = "DSR",
    deskew: bool = True,
    rotate: bool = True,
    overwrite: bool = False,
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamA_ch1", "CamB_ch0", "CamB_ch1"),
    dz: float = 0.5,
    xy_pixel_size: float = 0.108,
    skew_angle: float = 32.45,
    objective_scan: bool = False,
    z_stage_scan: bool = False,
    reverse: bool = False,
    flip_zstack: bool = False,
    parse_setting_file: bool = False,
    crop: bool = False,
    input_axis_order: str = "zyx",
    output_axis_order: str = "zyx",
    dsr_combined: bool = True,
    ff_correction: bool = False,
    bk_removal: bool = False,
    lower_limit: Any = 0.4,
    const_offset: Any = None,
    ff_image_paths: Sequence[str] = ("", "", ""),
    background_paths: Sequence[str] = ("", "", ""),
    save_16bit: bool = True,
    save_3d_stack: bool = True,
    save_mip: bool = True,
    large_file: bool = False,
    zarr_file: bool = False,
    save_zarr: bool = False,
    batch_size: Sequence[int] = (1024, 1024, 1024),
    block_size: Sequence[int] = (256, 256, 256),
    input_bbox: Sequence[int] | None = None,
    task_size: Any = None,
    resample_factor: Sequence[float] | None = None,
    interp_method: str = "linear",
    mask_fullpaths: Sequence[str] = (),
    parse_cluster: bool = True,
    parse_parfor: bool = False,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 1,
    uuid: str = "",
    debug: bool = False,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Deskew and/or rotate the files of the datasets ``data_paths`` that match a channel
    pattern (``XR_deskew_rotate_data_wrapper``), one task per file through the job framework;
    see the module's description. Raises if files are still not done after the retries."""
    values = {k: v for k, v in locals().items() if k in _ARGS}
    from petakit5d.core.ids import get_uuid
    from petakit5d.deskew_rotate.file import deskew_rotate_file
    from petakit5d.deskew_rotate.memory import estimate_deskew_rotate_memory
    from petakit5d.deskew_rotate.zarr import _OBJECTIVE, deskew_rotate_zarr
    from petakit5d.io import get_image_dtype, get_image_size
    from petakit5d.jobs import framework
    from petakit5d.jobs.config import cluster_config_value
    from petakit5d.jobs.coordinator import coordinator_nice
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_dataset_filenames

    paths = [data_paths] if isinstance(data_paths, str | os.PathLike) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    paths = [p if p.endswith("/") else p + "/" for p in paths]
    deskew, rotate = _logical(deskew), _logical(rotate)
    objective_scan, large_file = _logical(objective_scan), _logical(large_file)
    dsr_combined, save_zarr = _logical(dsr_combined), _logical(save_zarr)
    save_16bit, save_mip = _logical(save_16bit), _logical(save_mip)
    if not uuid:
        uuid = get_uuid()
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(paths[-1], job_log_dir).parse_cluster)

    if objective_scan:
        deskew = False
        if dsr_dir_name == "DSR/":
            dsr_dir_name = "rotated/"
    if large_file:
        if not rotate:
            raise DeskewRotateDatasetError(
                "Currently only combined DSR is supported for large scale deskew/rotation!\\n"
            )
        if objective_scan:
            raise DeskewRotateDatasetError(_OBJECTIVE)
        dsr_combined = True
        save_zarr = True
    # deskew only is done by the separate path (the combined path only writes DSR)
    if not rotate:
        dsr_combined = False
    ext = ".zarr" if save_zarr else ".tif"
    save_ds = deskew and (not rotate or not dsr_combined)
    if not save_ds and not rotate:
        print("Both deskew and rotation are disabled, skip the processing!", flush=True)
        return

    ds_paths: list[str] = []
    dsr_paths: list[str] = []
    for d in paths:
        ds_path = f"{d}{ds_dir_name}/"
        dsr_path = f"{d}{dsr_dir_name}/"
        ds_paths.append(ds_path)
        dsr_paths.append(dsr_path)
        if _logical(overwrite):
            import shutil

            if save_ds and os.path.isdir(ds_path):
                shutil.rmtree(ds_path)
            if rotate and os.path.isdir(dsr_path):
                shutil.rmtree(dsr_path)
        if save_ds:
            _make_folder(ds_path)
        if rotate:
            _make_folder(dsr_path)
        _save_parameters(ds_path if save_ds and not rotate else dsr_path, values)

    # the files: listed per folder, with their partial files, sizes and setting-file flips
    names = parse_dataset_filenames(
        paths,
        list(channel_patterns),
        _logical(parse_setting_file),
        _logical(flip_zstack),
        False,
        [],
        False,
        1,
        _logical(zarr_file),
    )
    n_files = len(names.fnames)

    px = _first(xy_pixel_size)
    box = None if _empty(input_bbox) else tuple(_ints(input_bbox))
    rs_task = _through_mat2str(resample_factor)
    block = tuple(_ints(block_size))
    inputs: list[str] = []
    outputs: list[str] = []
    tasks: list[Any] = []
    for f in range(n_files):
        fname = names.fnames[f]
        fsname = os.path.splitext(fname)[0]
        fd = int(names.fd_inds[f])
        data_path = paths[fd]
        frame = f"{data_path}{fname}"
        inputs.append(frame)
        if save_ds and not rotate:
            outputs.append(f"{ds_paths[fd]}{fsname}{ext}")
        if rotate:
            outputs.append(f"{dsr_paths[fd]}{fsname}{ext}")
        if large_file:
            masks = tuple(os.fspath(p) for p in mask_fullpaths) or ("",)
            tasks.append(
                Task.of(
                    deskew_rotate_zarr,
                    frame,
                    _through(px, "%.20d"),
                    _through(dz, "%.20d"),
                    result_dir_str=dsr_dir_name,
                    overwrite=_logical(overwrite),
                    objective_scan=objective_scan,
                    skew_angle=_through(skew_angle, "%d"),
                    reverse=_logical(reverse),
                    flip_zstack=bool(names.flip_zstack[f]),
                    crop=_logical(crop),
                    dsr_combined=dsr_combined,
                    save_16bit=save_16bit,
                    save_mip=save_mip,
                    save_zarr=save_zarr,
                    batch_size=list(_ints(batch_size)),
                    block_size=list(block),
                    input_bbox=None if box is None else list(box),
                    task_size=None if _empty(task_size) else _first(task_size),
                    resample_factor=rs_task,
                    interp_method=interp_method,
                    mask_fullpaths=masks,
                    parse_cluster=bool(parse_cluster),
                    parse_parfor=_logical(parse_parfor),
                    master_compute=_logical(master_compute),
                    job_log_dir=job_log_dir,
                    cpus_per_task=_through(_first(cpus_per_task), "%d"),
                    uuid=uuid,
                    debug=_logical(debug),
                    mcc_mode=_logical(mcc_mode),
                    config_file=config_file,
                )
            )
            continue
        ff_image = background_image = ""
        if _logical(ff_correction):
            # the first channel pattern the path contains (MATLAB's FFImagePaths{mask})
            hit = [k for k, p in enumerate(channel_patterns) if p in frame]
            ff_image = _first_selected(ff_image_paths, hit)
            background_image = _first_selected(background_paths, hit)
        # a single volume file, or the group of its partial files
        frame_paths = [frame]
        if bool(names.partialvols[f]):
            frame_paths = [f"{data_path}{x}" for x in names.gfnames[f]]
        tasks.append(
            Task.of(
                deskew_rotate_file,
                frame_paths,
                _through(px, "%.20d"),
                _through(dz, "%.20d"),
                ds_dir_name=ds_dir_name,
                dsr_dir_name=dsr_dir_name,
                skew_angle=_through(skew_angle, "%.20d"),
                objective_scan=objective_scan,
                z_stage_scan=_logical(z_stage_scan),
                reverse=_logical(reverse),
                crop=_logical(crop),
                input_axis_order=input_axis_order,
                output_axis_order=output_axis_order,
                ff_correction=_logical(ff_correction),
                bk_removal=_logical(bk_removal),
                lower_limit=_through(_first(lower_limit), "%.20d"),
                const_offset=None
                if _empty(const_offset)
                else _through(_first(const_offset), "%0.10f"),
                ff_image=ff_image,
                background_image=background_image,
                rotate=rotate,
                resample_factor=rs_task,
                dsr_combined=dsr_combined,
                input_bbox=None if box is None else list(box),
                flip_zstack=bool(names.flip_zstack[f]),
                save_16bit=save_16bit,
                save_3d_stack=_logical(save_3d_stack),
                save_zarr=save_zarr,
                block_size=list(block),
                save_mip=save_mip,
                interp_method=interp_method,
            )
        )

    if not inputs:
        raise DeskewRotateDatasetError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    coord_nice = 0
    if large_file:
        # (the files' jobs wait for their block jobs: these go first)
        if parse_cluster:
            coord_nice = coordinator_nice()
        # each file's job runs its blocks as jobs: it holds the process and the merge of the
        # MIP slabs, and a block's job's resources if it computes blocks itself
        cpus, mem_allocate = controller_resources(
            get_image_size(inputs[0]),
            px,
            dz,
            skew_angle,
            save_16bit,
            interp_method,
            resample_factor,
            box,
            batch_size,
            block,
            save_mip,
            get_image_dtype(inputs[0]),
            _logical(cluster_config_value(config_file, "masterCompute", master_compute)),
            _first(cpus_per_task),
        )
    else:
        # the arrays of the largest input (with all its partial files), plus the process
        fmax = int(np.argmax(np.asarray(names.data_sizes, dtype=np.float64)))
        fns = [inputs[fmax]]
        if bool(names.partialvols[fmax]):
            fns = [f"{paths[int(names.fd_inds[fmax])]}{x}" for x in names.gfnames[fmax]]
        sizes = [tuple(int(v) for v in get_image_size(p)) for p in fns]
        read_size = (sum(s[0] for s in sizes), sizes[0][1], sizes[0][2])
        dsz = read_size if box is None else _box_size(box)
        preprocess = (
            _logical(ff_correction)
            or _logical(bk_removal)
            or bool(np.any(names.flip_zstack))
            or input_axis_order != "zyx"
        )
        mem_allocate = estimate_deskew_rotate_memory(
            dsz,
            px,
            dz,
            skew_angle=skew_angle,
            objective_scan=objective_scan,
            z_stage_scan=_logical(z_stage_scan),
            dsr_combined=dsr_combined,
            rotate=rotate,
            save_16bit=save_16bit,
            interp_method=interp_method,
            resample_factor=resample_factor,
            read_size=read_size,
            partial_file=len(fns) > 1,
            preprocess=preprocess,
            output_axis_order=output_axis_order,
            save_zarr=save_zarr,
            block_size=block,
            data_type=get_image_dtype(fns[0]),
        )[0]
        cpus = _first(cpus_per_task)

    # three rounds of tries, each with twice the resources
    done: list[bool] = [False]
    for i in range(3):
        if not all(done):
            done = framework.run_tasks(
                inputs,
                outputs,
                tasks,
                cpus_per_task=cpus * 2**i,
                mem_allocate=mem_allocate * 2**i,
                max_trial_num=2,
                master_compute=_logical(master_compute),
                parse_cluster=parse_cluster,
                job_log_dir=job_log_dir,
                mcc_mode=mcc_mode,
                config_file=config_file,
                nice=coord_nice,
            )
    if not all(done):
        raise DeskewRotateDatasetError(f"{sum(not d for d in done)} / {len(done)} failed")


_ARGS = (
    "data_paths",
    "ds_dir_name",
    "dsr_dir_name",
    "deskew",
    "rotate",
    "overwrite",
    "channel_patterns",
    "dz",
    "xy_pixel_size",
    "skew_angle",
    "objective_scan",
    "z_stage_scan",
    "reverse",
    "flip_zstack",
    "parse_setting_file",
    "crop",
    "input_axis_order",
    "output_axis_order",
    "dsr_combined",
    "ff_correction",
    "bk_removal",
    "lower_limit",
    "const_offset",
    "ff_image_paths",
    "background_paths",
    "save_16bit",
    "save_3d_stack",
    "save_mip",
    "large_file",
    "zarr_file",
    "save_zarr",
    "batch_size",
    "block_size",
    "input_bbox",
    "task_size",
    "resample_factor",
    "interp_method",
    "mask_fullpaths",
    "parse_cluster",
    "parse_parfor",
    "master_compute",
    "job_log_dir",
    "cpus_per_task",
    "uuid",
    "debug",
    "mcc_mode",
    "config_file",
)
