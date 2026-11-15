"""MIP masks (MATLAB's ``XR_generate_MIP_masks`` and ``gererate_single_MIP_mask``): masks of the
major objects in the MIPs of a dataset, by intensity threshold and object size, with their
bounding boxes, for post-processing of big data.

The mask of a MIP image ``im`` (``generate_single_mip_mask``), as MATLAB computes it:
``im > int_thresh``; opened with ``strel('disk', 5)``; components below ``vol_thresh`` pixels
removed (``bwareaopen``, 8-connected); dilated by ``strel('disk', 5 + d)`` when both
dilate sizes are ``d``, else by ``strel('disk', 5)`` and lines of ``2 d + 1`` pixels along the
image's rows and columns; holes filled; small components removed again. MATLAB's strels are
periodic-line decompositions, not Euclidean disks (``petakit5d.matlabcompat.morphology``).

Files keep MATLAB's names and contents: the mask ``<result>/<name>.tif`` (``.zarr`` with
``save_zarr``; uint16, or single without ``save_16bit``), ``<name>_bbox.mat`` (v7.3: ``fn``,
``fnout``, ``intThresh``, ``volThresh``, ``sz`` and ``bbox`` with MATLAB's values: the size
``[ny nx]`` and the box ``[ymin xmin ymax xmax]``, 1-based, inclusive, ``[]`` for an empty
mask) and the figure ``<name>_plot.png`` (the MIP and its mask; only with matplotlib, the
optional ``figures`` dependency, as MATLAB draws it only with its JVM). The box file and the
figure are written before the mask, which marks the task as done; each file is written under
a temporary name and renamed.

A MIP image is the 2D image of the tiff as MATLAB reads it (``petakit5d.tools.mip``: along z
``ny x nx``, along y ``nx x nz``, along x ``ny x nz``), so ``dilate_size`` of
``generate_single_mip_mask`` is in the image's (rows, columns). The dataset function takes its
per-axis values in ``(z, y, x)``; the MIP of an axis gets the threshold of that axis and the
dilate sizes of the other two, in MATLAB's order (``y x z`` without the projected axis: the
image's rows and columns).

Kept from MATLAB: the axis of each file comes from the list of files before the axis filter
(``mip_axes`` is not filtered with the file names), so with an axis left out the files after
it get another axis's threshold and dilate sizes; the image size for the memory estimate is
that of the first ``.tif`` (``.zarr``) of the last folder, channel pattern or not; with a
single MIP file in all MATLAB fails (``fileparts`` of a 1x1 cell gives a cell, which it wraps
in another, and ``sprintf`` of the mask's name fails: ``MATLAB:printf:invalidInputType``),
after parameters.mat is written.
"""

from __future__ import annotations

import os
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = ["MipMaskError", "generate_mip_masks", "generate_single_mip_mask", "mip_mask"]

_AXIS_STRS = ("y", "x", "z")
_N_SE = 5


class MipMaskError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


class MipMaskWarning(UserWarning):
    """a warning MATLAB gives, with its identifier ('' for ``warning(message)``)"""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# the mask of an image
# ---------------------------------------------------------------------------------------------


def _greater(im: np.ndarray, thresh: float) -> np.ndarray:
    """MATLAB's ``im > thresh`` of a numeric image and a double, exactly: in double (measured:
    ``single(100.3) > 100.3`` is true; NumPy would round the threshold to single); 8- to 32-bit
    integers and single values are exact in double"""
    a = np.asarray(im)
    if a.dtype == np.bool_:
        a = a.astype(np.uint8)
    return a.astype(np.float64, copy=False) > float(thresh)


def mip_mask(
    im: np.ndarray,
    int_thresh: float = 100,
    vol_thresh: float = 1e5,
    dilate_size: Sequence[float] = (0, 0),
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """The mask (bool, the image's shape) of ``gererate_single_MIP_mask`` for an image ``im``
    (``(y, x)``; a ``(z, y, x)`` volume is processed as MATLAB processes one): see the module.
    ``dilate_size``: the dilation in the image's (rows, columns)."""
    from petakit5d.matlabcompat import label as L
    from petakit5d.matlabcompat import morphology as M

    t = num_threads
    d = [float(v) for v in np.asarray(dilate_size, dtype=np.float64).ravel()]
    if len(d) < 2:
        raise MipMaskError(
            f"Index exceeds the number of array elements. Index must not exceed {len(d)}.",
            "MATLAB:badsubscript",
        )
    bw = _greater(im, int_thresh)
    bw = M.imopen(bw, M.strel("disk", _N_SE), threads=t)
    bw = L.bwareaopen(bw, vol_thresh, threads=t)
    if d[0] == d[1]:
        bw = M.imdilate(bw, M.strel("disk", _N_SE + d[0]), threads=t)
    else:
        bw = M.imdilate(bw, M.strel("disk", _N_SE), threads=t)
        if d[0] > 0:
            bw = M.imdilate(bw, M.strel("line", d[0] * 2 + 1, 90), threads=t)
        if d[1] > 0:
            bw = M.imdilate(bw, M.strel("line", d[1] * 2 + 1, 0), threads=t)
    bw = L.imfill(bw, "hole", threads=t)
    return L.bwareaopen(bw, vol_thresh, threads=t)


def _matlab_bbox(bw: np.ndarray) -> np.ndarray:
    """getImageBoundingBox of the mask in MATLAB's values (double; ``[]`` when empty)"""
    from petakit5d.core.bbox import box_to_matlab
    from petakit5d.io.imagesize import get_image_bounding_box

    box = get_image_bounding_box(bw)
    if not box:
        return np.zeros((0, 0), dtype=np.float64)
    return np.asarray(box_to_matlab(box), dtype=np.float64)


def _read(fn: str) -> np.ndarray:
    from petakit5d.core.paths import fileparts

    ext = fileparts(fn)[2]
    if ext in (".tif", ".tiff"):
        from petakit5d.io import read_tiff

        return read_tiff(fn)
    if ext == ".zarr":
        from petakit5d.io import read_zarr

        return read_zarr(fn)
    # MATLAB: im is undefined
    raise MipMaskError("Unrecognized function or variable 'im'.", "MATLAB:UndefinedFunction")


def generate_single_mip_mask(
    fn: str,
    fnout: str,
    int_thresh: float = 100,
    vol_thresh: float = 1e5,
    dilate_size: Sequence[float] = (0, 0),
    block_size: Sequence[int] = (1, 256, 256),
    save_16bit: bool = True,
    uuid: str = "",
) -> None:
    """The mask of the MIP file ``fn`` as ``fnout`` (``.tif`` or ``.zarr``) with
    ``<name>_bbox.mat`` and, with matplotlib, ``<name>_plot.png`` next to it
    (``gererate_single_MIP_mask``; see the module). Nothing if ``fnout`` exists.
    ``dilate_size``: (rows, columns) of the image; ``block_size``: the zarr mask's blocks
    ``(z, y, x)``."""
    from petakit5d.core.arrays import size_to_matlab
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.io import write_mat
    from petakit5d.jobs.files import exist

    if not uuid:
        uuid = get_uuid()
    if exist(fnout):
        print(f"The output file {fnout} already exists, skip it!", flush=True)
        return
    im = _read(fn)
    sz = np.asarray(size_to_matlab(im.shape), dtype=np.float64)
    bw = mip_mask(im, int_thresh, vol_thresh, dilate_size)

    out_path, fsn, ext = fileparts(fnout)
    bbox = _matlab_bbox(bw)
    bbox_fnout = f"{out_path}/{fsn}_bbox.mat"
    bbox_tmp = f"{out_path}/{fsn}_bbox_{uuid}.mat"
    write_mat(
        bbox_tmp,
        {
            "fn": fn,
            "fnout": fnout,
            "intThresh": float(int_thresh),
            "volThresh": float(vol_thresh),
            "sz": sz,
            "bbox": bbox,
        },
    )
    os.replace(bbox_tmp, bbox_fnout)

    # the figure of the input and the mask (MATLAB: with its JVM)
    from petakit5d.tools import figures

    if figures.available():
        figures.save_mip_mask_figure(im, bw, f"{out_path}/{fsn}_plot.png")

    tmp = f"{out_path}/{fsn}_{uuid}{ext}"
    out = bw.astype(np.uint16) if save_16bit else bw.astype(np.float32)
    if ext == ".tif":
        from petakit5d.io import write_tiff

        write_tiff(out, tmp)
        os.replace(tmp, fnout)
    elif ext == ".zarr":
        from petakit5d.io import create_zarr, write_zarr
        from petakit5d.jobs.publish import publish_result

        vol = out.reshape((1, *out.shape)) if out.ndim == 2 else out
        create_zarr(
            tmp,
            data_size=vol.shape,
            block_size=tuple(int(v) for v in block_size),
            dtype="uint16" if save_16bit else "single",
        )
        write_zarr(vol, tmp)
        publish_result(tmp, fnout)
    else:
        # MATLAB writes nothing and then fails to move the missing file
        raise MipMaskError(
            f"No matching files or folders were found for '{tmp}'.", "MATLAB:MOVEFILE:FileNotFound"
        )


# ---------------------------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------------------------


def _yxz3(v: Any, name: str) -> list[float]:
    """a (z, y, x) value (or a scalar: repeated) as MATLAB's [y x z] list"""
    from petakit5d.core.axes import to_yxz

    a = np.asarray(v, dtype=np.float64).ravel()
    if a.size == 1:
        return [float(a[0])] * 3
    if a.size != 3:
        raise MipMaskError(f"{name}: one value or three (z, y, x)")
    return [float(x) for x in to_yxz(tuple(a))]


def _names(pattern: str) -> list[str]:
    from petakit5d.core.paths import matlab_dir

    return [e.name for e in matlab_dir(pattern)]


def _save_parameters(result_path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) of the call (MATLAB's ``pr``)"""
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(generate_mip_masks, values).items()}
    write_mat(result_path + "/parameters.mat", {"pr": pr})


@register(
    "XR_generate_MIP_masks",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders with the MIPs."),
    Param("ResultDirStr", "@ischar", "Result folder under each data path."),
    Param("axis", "@isnumeric", "MIP axes (z, y, x) to mask (nonzero).", convert="flags_yxz"),
    Param("channelPatterns", "@iscell", "Channel identifiers of the MIP files."),
    Param("zarrFile", "@islogical", "The MIPs are zarr."),
    Param("saveZarr", "@islogical", "Save the masks as zarr."),
    Param("blockSize", "@isnumeric", "Block size (z, y, x) of zarr masks.", convert="size_yxz"),
    Param(
        "intThresh",
        "@isnumeric",
        "Intensity threshold of the MIPs along z, y and x (or one for all).",
        convert="size_yxz",
    ),
    Param("volThresh", "@isnumeric", "Smallest object, in pixels."),
    Param(
        "dilateSize",
        "@isnumeric",
        "Dilation (z, y, x) for buffer room (or one for all).",
        convert="size_yxz",
    ),
    Param("save16bit", "@islogical", "uint16 masks (else single).", convert="same"),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="B10",
)
def generate_mip_masks(
    data_paths: str | Sequence[str],
    *,
    result_dir_str: str = "masks",
    axis: Sequence[float] = (1, 1, 1),
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamA_ch1", "CamB_ch0", "CamB_ch1"),
    zarr_file: bool = False,
    save_zarr: bool = False,
    block_size: Sequence[int] = (1, 256, 256),
    int_thresh: float | Sequence[float] = (100, 100, 100),
    vol_thresh: float = 100,
    dilate_size: float | Sequence[float] = 100,
    save_16bit: bool = True,
    uuid: str = "",
    debug: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 3,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Masks of the MIPs ``<name>_MIP_{y,x,z}.tif`` (``.zarr`` with ``zarr_file``) of datasets
    (``XR_generate_MIP_masks``): for each MIP of a flagged axis whose name matches a channel
    pattern, ``generate_single_mip_mask`` into ``<data path>/<result_dir_str>/`` with the
    threshold of its axis and the dilate sizes of the other two, one task per file through the
    job framework. parameters.mat of the call goes into each result folder."""
    values = {k: v for k, v in locals().items() if k in _WRAPPER_ARGS}
    from petakit5d.core.axes import to_yxz
    from petakit5d.core.dtypes import mround
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, simplify_path
    from petakit5d.core.text import through_mat2str, through_text
    from petakit5d.io import get_image_size
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import matches_channel_pattern

    if not uuid:
        uuid = get_uuid()
    paths = [data_paths] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    axis_yxz = [float(v) for v in to_yxz(tuple(axis))]

    result_paths = []
    for d, data_path in enumerate(paths):
        result_path = f"{data_path}/{result_dir_str}/"
        os.makedirs(result_path, exist_ok=True)
        paths[d] = simplify_path(data_path)
        result_paths.append(simplify_path(result_path))
        group_writable(result_path)
        _save_parameters(result_path, values)

    fnames: list[str] = []
    fd_inds: list[int] = []
    mip_axes: list[int] = []
    sz: tuple[int, ...] = ()
    ext = ".zarr" if zarr_file else ".tif"
    for d, data_path in enumerate(paths):
        names = _names(f"{data_path}/*{ext}")
        if not names:
            raise MipMaskError(
                "Index exceeds the number of array elements. Index must not exceed 0.",
                "MATLAB:badsubscript",
            )
        sz = tuple(get_image_size(f"{data_path}/{names[0]}"))
        flags = [matches_channel_pattern(names, p, True) for p in channel_patterns]
        names = [n for k, n in enumerate(names) if any(f[k] for f in flags)]
        # the axis of each file (0: none; the last matching pattern wins)
        axes_d = [0] * len(names)
        keep = [False] * len(names)
        for i in range(3):
            hit = [f"mip_{_AXIS_STRS[i]}" in n.lower() for n in names]
            for k in range(len(names)):
                if hit[k]:
                    axes_d[k] = i + 1
            if axis_yxz[i] > 0:
                if not any(hit):
                    warnings.warn(
                        MipMaskWarning(
                            f"MIPs for axis {_AXIS_STRS[i]} for data {data_path} do not exist!"
                        ),
                        stacklevel=2,
                    )
                keep = [a or b for a, b in zip(keep, hit, strict=True)]
        fnames.extend(n for n, k in zip(names, keep, strict=True) if k)
        fd_inds.extend(d for k in keep if k)
        # (MATLAB's quirk: not filtered by axis)
        mip_axes.extend(axes_d)

    n_files = len(fnames)
    if n_files == 1:
        # MATLAB: fileparts of a 1x1 cell gives a 1x1 cell, wrapped in another one, and the
        # sprintf of the result path fails
        raise MipMaskError(
            "Function is not defined for 'cell' inputs.", "MATLAB:printf:invalidInputType"
        )
    int_t = _yxz3(int_thresh, "int_thresh")
    dil = _yxz3(dilate_size, "dilate_size")
    blk = through_mat2str(np.asarray(to_yxz(tuple(block_size)), dtype=np.float64))
    block_py = tuple(int(v) for v in (blk[2], blk[0], blk[1]))
    vol_t = float(through_text(vol_thresh, "%d"))

    frame_paths, out_paths, tasks = [], [], []
    for f in range(n_files):
        d = fd_inds[f]
        fn = f"{paths[d]}/{fnames[f]}"
        fsn = fileparts(fnames[f])[1]
        fnout = f"{result_paths[d]}/{fsn}{'.zarr' if save_zarr else '.tif'}"
        ax = mip_axes[f]
        if ax == 0:
            raise MipMaskError(
                "Index in position 1 is invalid. Array indices must be positive integers or "
                "logical values.",
                "MATLAB:index:expected_one_positive_integer",
            )
        dil_f = [v for i, v in enumerate(dil) if i != ax - 1]
        dil_f = [float(v) for v in through_mat2str(np.asarray(dil_f, dtype=np.float64))]
        int_f = float(through_text(int_t[ax - 1], "%d"))
        frame_paths.append(fn)
        out_paths.append(fnout)
        tasks.append(
            Task.of(
                generate_single_mip_mask,
                fn,
                fnout,
                int_f,
                vol_t,
                dil_f,
                list(block_py),
                bool(save_16bit),
                uuid,
            )
        )

    task_batch_num = max(1, int(mround(n_files / 1000)))
    mem_allocate = float(np.prod(sz, dtype=np.float64)) * 4 / 1024**3 * 20
    framework.run_tasks(
        frame_paths,
        out_paths,
        tasks,
        parse_cluster=parse_cluster,
        master_compute=master_compute,
        task_batch_num=task_batch_num,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        mcc_mode=mcc_mode,
        config_file=config_file,
    )


_WRAPPER_ARGS = (
    "data_paths",
    "result_dir_str",
    "axis",
    "channel_patterns",
    "zarr_file",
    "save_zarr",
    "block_size",
    "int_thresh",
    "vol_thresh",
    "dilate_size",
    "save_16bit",
    "uuid",
    "debug",
    "parse_cluster",
    "master_compute",
    "cpus_per_task",
    "mcc_mode",
    "config_file",
)
