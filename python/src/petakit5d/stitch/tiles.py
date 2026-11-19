"""Tile preparation of the stitching: ``stitch_process_tiles`` and ``stitch_process_zarr_tile``.

- ``process_tiles`` (stitch_process_tiles, called by the frame ``stitch.volume.stitch_volume``):
  tiff tiles are converted to zarr by the tiff wrapper (``io.convert.tiff_to_zarr_dataset``,
  XR_tiffToZarr_wrapper) with the processing; zarr tiles that need processing (a result folder
  given) become one ``process_zarr_tile`` task each, run by the job framework in two rounds
  (the second with twice the CPUs and memory);
- ``process_zarr_tile`` (stitch_process_zarr_tile): one zarr tile, in MATLAB's order: read
  (the box ``input_bbox``), flip in z, crop to ``tile_out_bbox``, resample (MATLAB's
  ``imresize3`` to ``round(size ./ resample)``, ``matlabcompat.resize``), the user function,
  the cast back to the input's class, and the write of a native zarr under a temporary name,
  moved into place with one rename (an existing result is kept).

Tile functions (``process_fun_path``: function files, rows for channels, columns for
alternatives picked by their ``keywords``) are MATLAB code of the standard forms
(``stitch.tile_functions.parse_tile_function``) or Python functions by import path
(``io.convert.user_function``). A tile of one slice is an image ``(y, x)`` for the function,
as MATLAB has it. Locally, the tiles of both kinds are read into one buffer kept from tile to
tile (``io.convert.keep_read_buffer``), and a tile function written as MATLAB code runs on the
tile's own volume in place (the edge erosion): the same results.

Conventions: arrays ``(z, y, x)``; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive
(rows of boxes: one per tile location); ``block_size`` and ``shard_size`` ``(z, y, x)``;
``resample`` ``(z, y, x)`` or ``(z, xy)``. Values MATLAB passes to its tasks as text
(``mat2str``) keep its 15 significant digits.

The memory per task is MATLAB's: the last tile's voxels as single, times 2.5 (twice that
without ``big_data``), and twice that again in the second round.
"""

from __future__ import annotations

import contextlib
import math
import os
import shutil
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = ["TileError", "process_tiles", "process_zarr_tile"]

# memory of a zarr tile task: its voxels as single times this (twice without big data)
MEM_FACTOR = 2.5


class TileError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _print(msg: str = "") -> None:
    print(msg, flush=True)


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    return np.size(np.asarray(v, dtype=object)) == 0


# ---------------------------------------------------------------------------------------------
# stitch_process_zarr_tile
# ---------------------------------------------------------------------------------------------


@register(
    "stitch_process_zarr_tile",
    Param(
        "inputFilename",
        "@(x) ischar(x) || iscell(x)",
        "Zarr tile (or a list whose first is read).",
    ),
    Param("zarrFilename", "@ischar", "Zarr result."),
    Param(
        "frame",
        "@isnumeric",
        "Data in memory instead of the tile.",
        convert="volume",
    ),
    Param("Overwrite", "@islogical", "Replace an existing result."),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x).", convert="size_yxz"),
    Param(
        "shardSize", "@isnumeric", "Shard size (z, y, x) (Zarr v3 with shards).", convert="size_yxz"
    ),
    Param("expand2dDim", "@islogical", "Store an image as one z slice ([1, ny, nx])."),
    Param("flipZstack", "@islogical", "Flip the tile in z."),
    Param(
        "resample",
        "@(x) isempty(x) || isnumeric(x)",
        "Resampling factors (z, y, x) or (z, xy) (None: none).",
        convert="resample_yxz",
    ),
    Param(
        "inputBbox", "@(x) isnumeric(x)", "Region (z0, y0, x0, z1, y1, x1) to read.", convert="bbox"
    ),
    Param(
        "tileOutBbox",
        "@(x) isempty(x) || isnumeric(x)",
        "Region (z0, y0, x0, z1, y1, x1) to keep after the flip.",
        convert="bbox",
    ),
    Param("compressor", "@ischar", "Zarr compressor."),
    Param(
        "usrFcn",
        "@(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x)",
        "Tile function applied after the resampling (MATLAB code, callable or import path).",
    ),
    Param("uuid", "@ischar", "Id of the temporary result (empty: a new one)."),
    task="S02",
)
def process_zarr_tile(
    input_filename: str | Sequence[str],
    zarr_filename: str,
    frame: np.ndarray | None = None,
    *,
    overwrite: bool = False,
    block_size: Sequence[int] = (250, 500, 500),
    shard_size: Sequence[int] | None = None,
    expand_2d_dim: bool = True,
    flip_zstack: bool = False,
    resample: Sequence[float] | None = None,
    input_bbox: Sequence[int] | None = None,
    tile_out_bbox: Sequence[int] | None = None,
    compressor: str = "lz4",
    usr_fcn: Any = "",
    uuid: str = "",
) -> None:
    """Process a zarr tile (stitch_process_zarr_tile): read it (its box ``input_bbox``), flip
    it in z, crop it to ``tile_out_bbox``, resample it to ``round(size ./ resample)`` (MATLAB's
    ``imresize3``), apply ``usr_fcn``, cast it back to its class and write it as a native zarr
    with chunks of ``block_size`` (cut to the result's size). An existing result is kept
    (``overwrite``: replaced first). A ``frame`` (data in memory) is processed instead of the
    tile."""
    from petakit5d.core.dtypes import cast
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import convert as C
    from petakit5d.io import create_zarr, read_zarr, write_zarr

    t0 = time.perf_counter()
    zarr_filename = os.fspath(zarr_filename).rstrip("/")
    if os.path.isdir(zarr_filename):
        if not overwrite:
            _print(f"Zarr result {zarr_filename} already exists, skip it!")
            return
        shutil.rmtree(zarr_filename)
    if not uuid:
        uuid = get_uuid()
    if isinstance(input_filename, (str, os.PathLike)):
        first = os.fspath(input_filename)
    else:
        files = [os.fspath(f) for f in input_filename]
        first = files[0] if files else ""
    _print(f"Process zarr tile for {first if first else zarr_filename} ...")

    # MATLAB's block size [y x z] (or [y x])
    bs = [float(v) for v in C._msize(C._ints(block_size))]
    have_frame = frame is not None and np.size(frame) > 0
    if have_frame:
        # the data given: its size, the chunks cut to it only at the end (an image takes
        # [y x] chunks)
        im = np.asarray(frame)
        sz = C._msize(im.shape)
        if im.ndim == 2:
            bs = bs[:2]
    else:
        # (into a loop's kept buffer: memory the tiles before it touched)
        im = read_zarr(
            first,
            input_bbox=None if _empty(input_bbox) else _read_box(input_bbox),
            out=C._kept_buffer(),
        )
        sz = C._msize(im.shape)
        if im.ndim == 2:
            if expand_2d_dim:
                bs = (bs + [1.0])[:3]
                bs[2] = 1
                sz = [*sz, 1]
            else:
                bs = bs[:2]
        bs = C._mmin(bs, sz)
    dtype = C._class(im)
    if any(s == 0 for s in sz):
        _print("The input image is empty, skip it!")
        return

    view = im
    if flip_zstack and view.ndim == 3:
        view = view[::-1]
    if not _empty(tile_out_bbox):
        view = _crop(view, _box(tile_out_bbox))
    im = C._materialize(view)
    del view

    rs = C._resample_factors(resample)
    if rs is not None:
        if im.ndim != 3:
            raise TileError(
                "Arrays have incompatible sizes for this operation.",
                "MATLAB:sizeDimensionsMustMatch",
            )
        from petakit5d.matlabcompat.resize import imresize3

        im = imresize3(im, [C._mround(n / r) for n, r in zip(im.shape, rs, strict=True)])

    fcn = C.user_function(usr_fcn)
    if fcn is not None:
        _print(C._function_text(usr_fcn))
        # (the tile is this call's own unless it is the caller's frame)
        owned = not have_frame or not np.may_share_memory(im, frame)
        im = C._m2d(np.asarray(C._call(fcn, im, owned)))
    sz = C._msize(im.shape)
    if len(sz) == 2:
        if expand_2d_dim:
            sz = [*sz, 1]
            bs = (list(bs) + [1.0])[:3]
            bs[2] = 1
        else:
            prod = math.prod(bs)
            if sz[0] <= sz[1]:
                bs = C._mmin(sz, [sz[0], C._mround(prod / sz[0])])
            else:
                bs = C._mmin(sz, [C._mround(prod / sz[1]), sz[1]])
    else:
        bs = C._mmin(sz, bs)
    bs = [int(v) for v in bs]

    tmp = f"{zarr_filename}_{uuid}"
    if not os.path.isdir(tmp):
        nchunks = math.prod(math.ceil(s / b) for s, b in zip(sz, bs, strict=True))
        create_zarr(
            tmp,
            data_size=C._pysize(sz),
            block_size=C._pysize(bs),
            shard_size=None if _empty(shard_size) else shard_size,
            dtype=dtype,
            expand_2d_dim=expand_2d_dim,
            compressor=compressor,
            dim_separator="/" if nchunks > 10000 else ".",
        )
    data = C._materialize(im if C._class(im) == dtype else np.asarray(cast(im, dtype)))
    del im
    write_zarr(data, tmp, block_size=C._pysize(bs))
    del data
    # (a result another writer published meanwhile is kept, also with overwrite)
    C._publish(tmp, zarr_filename, "Zarr")
    _print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.")
    _print("Done!")


def _box(v: Any) -> list[int | None]:
    """a box's values (None: an open end)"""
    return [None if x is None else int(x) for x in v]


def _read_box(v: Any) -> list[Any]:
    """a box to read (None: an open end; NaN and Inf kept, for read_zarr's refusal)"""
    return [x if x is None or not math.isfinite(float(x)) else int(x) for x in v]


def _crop(a: np.ndarray, box: Sequence[int | None]) -> np.ndarray:
    """crop3d (a view; MATLAB's error outside the data)"""
    from petakit5d.io import convert as C

    try:
        return C._crop(a, box)
    except C.ConvertError as err:
        raise TileError(str(err), err.identifier) from None


# ---------------------------------------------------------------------------------------------
# stitch_process_tiles
# ---------------------------------------------------------------------------------------------


def _function_table(
    process_fun_path: Any, n_channels: int
) -> tuple[list[list[Any]], list[list[Any]]]:
    """MATLAB's usrFcn_strs and keywords cells: n_channels (or more) rows of the function
    files' columns; ``None`` stands for a cell MATLAB fills with [] when a row beyond the
    channels is assigned"""
    from petakit5d.io import convert as C

    if isinstance(process_fun_path, (str, os.PathLike)) and os.fspath(process_fun_path):
        # (MATLAB: processFunPath{i, j} of a char)
        raise TileError(
            "Brace indexing is not supported for variables of this type.",
            "MATLAB:cellRefFromNonCell",
        )
    rows = C._fun_rows(process_fun_path)
    n_cols = max((len(r) for r in rows), default=0)
    usr: list[list[Any]] = [[""] * n_cols for _ in range(n_channels)]
    kws: list[list[Any]] = [[[""] for _ in range(max(1, n_cols))] for _ in range(n_channels)]
    _print("Process function paths:")
    for r in rows:
        _print("    " + "    ".join(r))
    for i, r in enumerate(rows):
        for j, fn in enumerate(r):
            if not fn:
                continue
            fun, kw = C._read_function_file(fn)
            while len(usr) <= i:
                usr.append([None] * n_cols)
                kws.append([None] * max(1, n_cols))
            if fun:
                usr[i][j] = fun
            if kw is not None:
                kws[i][j] = C._keywords(kw)
    return usr, kws


def _pick_function(path: str, usr: list[list[Any]], kws: list[list[Any]], cind: list[bool]):
    """the tile's function: of its channels' rows (column-major), the first whose keywords the
    path contains, else the first"""
    rows = [k for k, c in enumerate(cind) if c]
    n_cols = len(usr[0]) if usr else 0
    funs = [usr[r][j] for j in range(n_cols) for r in rows]
    keys = [kws[r][j] for j in range(max(1, n_cols)) for r in rows]
    if not funs:
        # (MATLAB: usrFcn_strs_c{1} of an empty cell)
        raise TileError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    for j, f in enumerate(funs):
        if any(k in path for k in (keys[j] or [])):
            return f
    return funs[0]


def _check_functions(process_fun_path: Any) -> None:
    """MATLAB code of a form the port doesn't read fails before any tile is converted"""
    from petakit5d.io import convert as C

    if isinstance(process_fun_path, (str, os.PathLike)):
        return
    for row in C._fun_rows(process_fun_path):
        for fn in row:
            if fn:
                fun, _ = C._read_function_file(fn)
                if fun:
                    C.user_function(fun)


def _text_values(v: Any) -> list[float] | None:
    """a numeric setting after MATLAB's mat2str (15 significant digits)"""
    from petakit5d.core.text import through_mat2str

    if _empty(v):
        return None
    return [float(x) for x in np.ravel(through_mat2str(np.asarray(v, dtype=np.float64)))]


def process_tiles(
    input_fullpaths: Sequence[str],
    *,
    result_dir_name: str = "zarr",
    zarr_file: bool = False,
    loc_ids: Sequence[float] | None = None,
    block_size: Sequence[int] = (250, 500, 500),
    shard_size: Sequence[int] | None = None,
    flipped_tile: Sequence[bool] | None = None,
    resample: Sequence[float] | None = None,
    partial_file: bool = False,
    channel_patterns: Sequence[str] = (".",),
    input_bbox: Any = None,
    tile_out_bbox: Any = None,
    process_fun_path: Any = "",
    parse_cluster: bool = True,
    big_data: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 1,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Prepare the tiles of a stitching (stitch_process_tiles): tiff tiles to zarr under
    ``<data>/<result_dir_name>/`` with the processing (the tiff wrapper), or zarr tiles through
    one ``process_zarr_tile`` task each when there is a result folder and something to do (a
    box or a function file). ``input_bbox`` and ``tile_out_bbox``: a box, or one row per tile
    location (``loc_ids``) when there are several."""
    from petakit5d.core.paths import fileparts
    from petakit5d.io import convert as C

    paths = [os.fspath(p) for p in input_fullpaths]
    _check_functions(process_fun_path)
    if not zarr_file:
        _print("\nConvert tiles from Tiff to Zarr...")
        C.tiff_to_zarr_dataset(
            "",
            tiff_fullpaths=paths,
            result_dir_name=result_dir_name,
            loc_ids=None if loc_ids is None else list(loc_ids),
            block_size=tuple(block_size),
            shard_size=None if _empty(shard_size) else tuple(shard_size),
            flipped_tile=None if flipped_tile is None else list(flipped_tile),
            resample_factor=None if _empty(resample) else tuple(resample),
            partial_file=partial_file,
            channel_patterns=tuple(channel_patterns),
            input_bbox=input_bbox,
            tile_out_bbox=tile_out_bbox,
            process_fun_path=process_fun_path,
            job_log_dir=job_log_dir,
            parse_cluster=parse_cluster,
            master_compute=master_compute,
            big_data=big_data,
            cpus_per_task=cpus_per_task,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )
        return
    if not result_dir_name:
        return

    _print("\nPreprocess Zarr tiles...")
    data_path = fileparts(paths[0])[0]
    if parse_cluster:
        from petakit5d.jobs.slurm import check_slurm_cluster

        parse_cluster = bool(check_slurm_cluster(data_path, job_log_dir).parse_cluster)

    patterns = list(channel_patterns)
    usr, kws = _function_table(process_fun_path, len(patterns))
    in_rows, out_rows = C._box_rows(input_bbox), C._box_rows(tile_out_bbox)
    if not in_rows and not out_rows and not any(len(r) for r in usr):
        return

    loc = [float(v) for v in (loc_ids if loc_ids is not None else ())]
    uniq = sorted(set(loc))
    n_loc = len(uniq)
    loc_specific = n_loc > 1 and (len(in_rows) == n_loc or len(out_rows) == n_loc)
    if loc_specific:
        for what, r in (("InputBbox", in_rows), ("tileOutBbox", out_rows)):
            if len(r) > 1 and len(r) != n_loc:
                raise TileError(
                    f"The number of {what} ({len(r)}) does not match the number of tile "
                    f"locations ({n_loc})!"
                )
        in_rows = in_rows if len(in_rows) > 1 else [in_rows[0] if in_rows else None] * n_loc
        out_rows = out_rows if len(out_rows) > 1 else [out_rows[0] if out_rows else None] * n_loc

    from petakit5d.jobs.tasks import Task

    flips = [] if flipped_tile is None else [bool(v) for v in flipped_tile]
    rs = _text_values(resample)
    outs: list[str] = []
    tasks = []
    for i, p in enumerate(paths):
        dpath, fsn, _ = fileparts(p)
        zpath = f"{dpath}/{result_dir_name}"
        if not os.path.isdir(zpath):
            from petakit5d.io.tiff import group_writable

            os.makedirs(zpath, exist_ok=True)
            group_writable(zpath)
        outs.append(f"{dpath}/{result_dir_name}/{fsn}.zarr/")
        cind = [pat in p for pat in patterns]
        if not any(cind):
            raise TileError(
                f"The file {p} does not match any channel patterns {' '.join(patterns)}"
            )
        if loc_specific:
            k = uniq.index(loc[i])
            in_i, out_i = in_rows[k], out_rows[k]
        else:
            if len(in_rows) > 1 or len(out_rows) > 1:
                # (MATLAB passes the rows on as one box, which the tile's read rejects)
                raise TileError(
                    "The bounding box must have 6 values (several boxes need one per tile "
                    "location)."
                )
            in_i = in_rows[0] if in_rows else None
            out_i = out_rows[0] if out_rows else None
        fun_i = _pick_function(p, usr, kws, cind)
        tasks.append(
            Task.of(
                process_zarr_tile,
                [p],
                outs[-1],
                None,
                block_size=[int(v) for v in block_size],
                shard_size=None if _empty(shard_size) else [int(v) for v in shard_size],
                flip_zstack=flips[i] if flips else False,
                resample=rs,
                input_bbox=None if _empty(in_i) else list(in_i),
                tile_out_bbox=None if _empty(out_i) else list(out_i),
                compressor="zstd",
                usr_fcn=fun_i or "",
            )
        )

    # (MATLAB sizes the tasks by the last tile: its loop's variable)
    from petakit5d.io import get_image_size

    s = [float(v) for v in get_image_size(paths[-1])]
    s3 = s if len(s) == 3 else [1.0, *s]
    mem = s3[1] * s3[2] * s3[0] * 4 / 1024**3 * MEM_FACTOR
    if not big_data:
        mem *= 2

    from petakit5d.jobs import framework

    common = {
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "max_trial_num": 2,
        "job_log_dir": job_log_dir,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    # locally, the tiles' reads keep their buffer from tile to tile
    with C.keep_read_buffer() if not parse_cluster else contextlib.nullcontext():
        done = framework.run_tasks(
            paths, outs, tasks, cpus_per_task=cpus_per_task, mem_allocate=mem, **common
        )
        if not all(done):
            framework.run_tasks(
                paths, outs, tasks, cpus_per_task=cpus_per_task * 2, mem_allocate=mem * 2, **common
            )
