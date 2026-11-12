"""Image list generation: csv files of stage coordinates for stitching.

Ports ``XR_generate_image_list_wrapper`` and the generators it calls:

- ``tile_list``: tile names with their indices (t, c, x, y, z) and the tile interval in um;
- ``tile_position``: the tiles in the folders, their indices from the names (the LLSM grammar,
  else the ``tilePatterns``), positions from the image size and the overlap;
- ``encoder``: the first MCS2 position of each tile in the ``SliceList_*part<N>.csv`` files
  (``XR_estimate_actual_step_size_from_encoder`` writes ``encoder_info/`` first);
- ``sqlite``: the acquisition database (``*.sqlite3``, view ``ImageListView_NoNull``, through
  the standard ``sqlite3`` module), with names mapped to the standard grammar.

The csv files are written as R2024b's ``writetable`` writes them (``imagelist.write_table``),
through a temporary name and one rename. Also here: ``stitch_generate_imagelist_from_encoder_ISM``
and ``stitch_generate_imagelist_from_sqlite_old_format``, which nothing calls.
"""

from __future__ import annotations

import math
import os
import sqlite3
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.core.dtypes import NAN, mround
from petakit5d.core.paths import fileparts, matlab_dir, simplify_path
from petakit5d.params import Param, register
from petakit5d.stitch.imagelist import (
    TILE_EXPRESSION_CAM,
    TILE_EXPRESSION_NOCAM,
    ImageListError,
    Table,
    _parse_names,
    _sprintf,
    _warn,
    contains_pattern,
    matlab_regex,
    parse_tile_info,
    read_table,
    unique_num,
    unique_text,
    write_table,
)

__all__ = [
    "estimate_actual_step_size_from_encoder",
    "generate_image_list",
    "generate_image_list_from_encoder",
    "generate_image_list_from_encoder_ism",
    "generate_image_list_from_sqlite",
    "generate_image_list_from_sqlite_old_format",
    "generate_image_list_from_tile_list",
    "generate_image_list_from_tile_positions",
]

_STAGE = (
    "StageX_um_",
    "StageY_um_",
    "StageZ_um_",
    "ObjectiveX_um_",
    "ObjectiveY_um_",
    "ObjectiveZ_um_",
)


def _texts(values) -> np.ndarray:
    values = list(values)
    out = np.empty(len(values), dtype=object)
    out[:] = list(values)
    return out


def _group_writable(path: str) -> None:
    from petakit5d.io.tiff import group_writable

    group_writable(path)


def _publish(t: Table, folder: str, name: str, uuid: str) -> str:
    """writetable to ``<folder>/<name>_<uuid>.csv``, group-writable, renamed to ``<name>.csv``."""
    tmpout = f"{folder}/{name}_{uuid}.csv"
    fnout = f"{folder}/{name}.csv"
    write_table(t, tmpout)
    _group_writable(tmpout)
    os.replace(tmpout, fnout)
    return fnout


def _stage_table(filepaths, filenames, xyz, mapped=None) -> Table:
    cols: dict[str, Any] = {"Filepath": _texts(filepaths), "Filename": _texts(filenames)}
    if mapped is not None:
        cols["mappedFilename"] = _texts(mapped)
    xyz = np.asarray(xyz, dtype=np.float64).reshape(-1, 3)
    for j, name in enumerate(_STAGE):
        cols[name] = xyz[:, j % 3].copy()
    return Table(cols, height=len(filepaths))


def _uuid(uuid: str) -> str:
    if uuid:
        return uuid
    from petakit5d.core.ids import get_uuid

    return get_uuid()


def _sind(x: float) -> float:
    """MATLAB's ``sind``: exact at multiples of 30 and 45 degrees, else the sine of the angle
    reduced to the nearest multiple of 90 (in radians as ``pi/180 * x``)."""
    if not math.isfinite(x):
        return NAN
    sign = -1.0 if x < 0 else 1.0
    a = math.fmod(abs(x), 360.0)
    exact = {0: 0.0, 30: 0.5, 90: 1.0, 150: 0.5, 180: 0.0, 210: -0.5, 270: -1.0, 330: -0.5}
    if a == math.floor(a) and int(a) in exact:
        return sign * exact[int(a)]
    n = int(math.floor(a / 90.0 + 0.5))
    t = math.radians(a - 90.0 * n)
    return sign * (math.sin(t), math.cos(t), -math.sin(t), -math.cos(t))[n % 4]


def _mode_rows(a: np.ndarray) -> np.ndarray:
    """``mode(a)`` of a matrix: per column the most frequent value, the smallest at ties."""
    out = []
    for col in a.T:
        vals, counts = np.unique(col, return_counts=True)
        out.append(vals[np.argmax(counts)])
    return np.asarray(out, dtype=np.float64)


def _positions(tab: Table, ts: np.ndarray, step: np.ndarray) -> np.ndarray:
    """xyz of each tile, per time point, camera, stack and channel: (index - ts) .* step."""
    xyz = np.zeros((tab.height, 3))
    full_iter = unique_text(tab["fullIter"])
    ch = unique_num(tab["ch"])
    cam = sorted(set(tab["camera"].tolist()))
    stackn = unique_num(tab["stack"])
    idx = np.column_stack([tab["x"], tab["y"], tab["z"]])
    for n in full_iter:
        for c_cam in cam:
            for s in stackn:
                for c in ch:
                    inds = (
                        (tab["ch"] == c)
                        & (tab["camera"] == c_cam)
                        & (tab["fullIter"] == n).astype(bool)
                        & (tab["stack"] == s)
                    )
                    xyz[inds, :] = (idx[inds, :] - ts) * step
    return xyz


def _check_cellstr_column(col: np.ndarray) -> None:
    # unique of a cell holding [] (a row never assigned) is an error in MATLAB
    if any(v is None for v in col):
        raise ImageListError(
            "Cell array input must be a cell array of character vectors.",
            "MATLAB:UNIQUE:InputClass",
        )


# ------------------------------------------------------------------------------------------
# tile list
# ------------------------------------------------------------------------------------------


def _mapped_name(full_iter: str, cam: str | None, ch, stack, x, y, z, t) -> str:
    if cam is None:
        return _sprintf(
            "Scan_Iter_%s_ch%d_CAM1_stack%04d_488nm_000000msec_000000msecAbs_%03dx_%03dy_%03dz_%04dt.tif",
            full_iter,
            ch,
            stack,
            x,
            y,
            z,
            t,
        )
    return _sprintf(
        "Scan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_488nm_000000msec_000000msecAbs_%03dx_%03dy_%03dz_%04dt.tif",
        full_iter,
        cam,
        ch,
        stack,
        x,
        y,
        z,
        t,
    )


def generate_image_list_from_tile_list(
    data_paths: str | Sequence[str],
    tile_filenames: str | Sequence[str],
    tile_indices,
    tile_interval,
    *,
    uuid: str = "",
) -> str:
    """``stitch_generate_imagelist_from_tile_list``: ``ImageList_from_tile_list.csv`` in the
    first folder, positions ``(index - min index) .* tile_interval`` (xyz, um), the names
    mapped to the standard grammar (camera A). Returns the file's path."""
    uuid = _uuid(uuid)
    if isinstance(data_paths, str):
        data_paths = [data_paths]
    fnames = [tile_filenames] if isinstance(tile_filenames, str) else list(tile_filenames)
    nf = len(fnames)
    ti = np.asarray(tile_indices, dtype=np.float64).reshape(-1, 5)
    if ti.shape[0] != nf:
        raise ImageListError(
            _sprintf(
                "The number of rows of tileIndices (%d) does not match the number of tiles (%d)!",
                ti.shape[0],
                nf,
            )
        )
    interval = np.asarray(tile_interval, dtype=np.float64).reshape(-1)
    filepaths = [data_paths[0] + x for x in fnames]
    full_iter = [_sprintf("%04d", v) for v in ti[:, 0]]
    tab = Table(
        {
            "Filepath": _texts(filepaths),
            "Filename": _texts(fnames),
            "prefix": _texts(["Scan"] * nf),
            "Iter": ti[:, 0],
            "subIter": _texts([""] * nf),
            "fullIter": _texts(full_iter),
            "camera": np.full(nf, "A", dtype="<U1"),
            "ch": ti[:, 1],
            "stack": np.zeros(nf),
            "x": ti[:, 2],
            "y": ti[:, 3],
            "z": ti[:, 4],
            "t": ti[:, 0],
        },
        height=nf,
    )
    tab["mappedFilename"] = _texts(
        _mapped_name(full_iter[f], "A", ti[f, 1], 0.0, ti[f, 2], ti[f, 3], ti[f, 4], ti[f, 0])
        for f in range(nf)
    )
    idx = np.column_stack([tab["x"], tab["y"], tab["z"]])
    ts = idx.min(axis=0) if nf else np.zeros(3)
    xyz = _positions(tab, ts, interval)
    t = _stage_table(tab["Filepath"], tab["Filename"], xyz, tab["mappedFilename"])
    fnout = _publish(t, data_paths[0], "ImageList_from_tile_list", uuid)
    print(f"Image list file {fnout} is generated!")
    return fnout


# ------------------------------------------------------------------------------------------
# tile positions
# ------------------------------------------------------------------------------------------


def _parse_image_filenames(data_paths: Sequence[str], zarr_file: bool, channel_patterns):
    """``parseImageFilenames`` (stand-in until P01's port lands): the tiff (or zarr) names of
    the folders that match a channel pattern (any case), with their folder index and path."""
    print("Parse image filenames for data path(s):\n    " + "\n    ".join(data_paths))
    fnames: list[str] = []
    fd_inds: list[int] = []
    for d, dp in enumerate(data_paths):
        ext = ".zarr" if zarr_file else ".tif"
        names = [e.name for e in matlab_dir(f"{dp}/*{ext}")]
        flag = np.zeros(len(names), dtype=bool)
        for p in channel_patterns:
            flag |= contains_pattern(names, p, ignore_case=True)
        sel = [n for n, k in zip(names, flag, strict=True) if k]
        fnames.extend(sel)
        fd_inds.extend([d] * len(sel))
    fsns = [fileparts(f)[1] for f in fnames]
    if not fnames:
        _warn(
            "No Zarr format (*.zarr) image files found in the input data folder(s)! Please "
            "check if dataPaths are correct!"
            if zarr_file
            else "No Tiff format (*.tif) image files found in the input data folder(s)! Please "
            "check if dataPaths are correct!"
        )
    filepaths = [f"{data_paths[d]}/{f}" for f, d in zip(fnames, fd_inds, strict=True)]
    if not fnames:
        _warn(
            "The input data folder(s) is empty, or no file matches the channel pattern, or file "
            "type (tiff/zarr) is not correctly setup. Please check."
        )
    return fnames, fsns, fd_inds, filepaths


def generate_image_list_from_tile_positions(
    data_paths: str | Sequence[str],
    *,
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamA_ch1", "CamB_ch0"),
    tile_patterns: Sequence[str] = ("0000t", "ch0", "000x", "000y", "000z"),
    ds: bool = False,
    dsr: bool = False,
    xy_pixel_size=0.108,
    dz=0.2,
    skew_angle=32.45,
    axis_order: str = "x,y,z",
    data_order: str = "z,y,x",
    objective_scan: bool = False,
    io_scan: bool = False,
    zarr_file: bool = False,
    overlap_size=None,
    overlap_size_type: str = "pixel",
    uuid: str = "",
) -> str:
    """``stitch_generate_imagelist_from_tile_positions``: ``ImageList_from_tile_positions.csv``
    in the first folder; tile steps are the tile size (the most common one, in stage
    coordinates of the scan type) minus the overlap. Names outside the LLSM grammar get a
    ``mappedFilename`` column. ``overlap_size`` is ``(z, y, x)``, in pixels or in um
    (``overlap_size_type``). ``axis_order`` and ``data_order`` (``'z,y,x'``, the array's axes)
    are accepted and unused, as in MATLAB. Returns the file's path."""
    from petakit5d.io.imagesize import get_image_size_batch

    uuid = _uuid(uuid)
    if isinstance(data_paths, str):
        data_paths = [data_paths]
    px = np.asarray(xy_pixel_size, dtype=np.float64).reshape(-1)
    pixel_sizes = np.array([px[0], px[-1], float(np.asarray(dz).reshape(-1)[0])])

    fnames, fsns, _, filepaths = _parse_image_filenames(data_paths, zarr_file, channel_patterns)
    nf = len(fnames)
    tile_info, cam_patterns = parse_tile_info(fsns, tile_patterns)

    cam_rx = matlab_regex(r"_Cam\w_ch")
    ch_rx = matlab_regex(r"_ch[0-9]_")
    specify_cam = True
    if all(cam_rx.search(s) for s in fsns):
        expression = TILE_EXPRESSION_CAM
    elif all(ch_rx.search(s) for s in fsns):
        expression, specify_cam = TILE_EXPRESSION_NOCAM, False
    else:
        expression, specify_cam = "", False
    rx = matlab_regex(expression) if expression else None
    standard_naming = all(rx is not None and rx.search(s) for s in fsns)

    base = Table({"Filepath": _texts(filepaths), "Filename": _texts(fnames)}, height=nf)
    if standard_naming:
        tab = _parse_names(base, fsns, expression, specify_cam)
    else:
        tab = _nonstandard_tab(base, tile_info, cam_patterns, specify_cam)

    nf = tab.height
    # rows (nz, ny, nx), (ny, nx) for a 2D zarr array; mode per column (order-free)
    im_sizes = get_image_size_batch(list(tab["Filepath"])).astype(np.float64)
    idx = np.column_stack([tab["x"], tab["y"], tab["z"]])
    ts = idx.min(axis=0)
    sz = _mode_rows(im_sizes)
    if sz.size < 3:
        sz = np.concatenate((np.ones(3 - sz.size), sz))

    # coordinates conversion for skewed space and DS space data
    if objective_scan or io_scan:
        z_aniso = pixel_sizes[2] / pixel_sizes[0]
        ds = dsr = False
    else:
        z_aniso = _sind(float(skew_angle)) * pixel_sizes[2] / pixel_sizes[0]
    theta = float(skew_angle) * math.pi / 180

    # overlap (z, y, x) -> xyz, in um
    ov = np.asarray(overlap_size if overlap_size is not None else [], dtype=np.float64)
    ov = ov.reshape(-1)
    if overlap_size_type.lower() == "pixel":
        ov = ov[[2, 1, 0]] * pixel_sizes[0] * np.array([1.0, 1.0, z_aniso])
    elif overlap_size_type.lower() == "um":
        ov = ov[::-1]
    else:
        raise ImageListError("Unknown overlap size type, must be pixel or um!")

    # physical size in xyz order, from (nz, ny, nx)
    sz_phy = sz[[2, 1, 0]] * pixel_sizes[0] * np.array([1.0, 1.0, z_aniso])
    sz_dsr = sz_phy
    st, ct = math.sin(theta), math.cos(theta)
    if objective_scan or ds:
        # size in DS space to DSR space (stage coordinates)
        sz_dsr = np.array(
            [sz_phy[0] * ct + sz_phy[2] * st, sz_phy[1], sz_phy[0] * st - sz_phy[2] * ct]
        )
    if not io_scan and not objective_scan and not ds and not dsr:
        # size in skewed space to DSR space (stage coordinates)
        sz_dsr = np.array([sz_phy[0] * ct + sz_phy[2] / st, sz_phy[1], sz_phy[0] * st])
    sz_1 = sz_dsr - ov

    _check_cellstr_column(tab["fullIter"])
    xyz = _positions(tab, ts, sz_1)
    mapped = None if standard_naming else tab["mappedFilename"]
    t = _stage_table(tab["Filepath"], tab["Filename"], xyz, mapped)
    fnout = _publish(t, data_paths[0], "ImageList_from_tile_positions", uuid)
    print(f"Image list file {fnout} is generated!")
    return fnout


def _nonstandard_tab(base: Table, tile_info, cam_patterns, specify_cam) -> Table:
    """The table of names outside the grammar, from parseTileInfo's indices. Rows with a
    missing index keep MATLAB's fill values (``None`` for ``[]``), and are not removed."""
    n = base.height
    cols: dict[str, Any] | None = None
    for f in range(n):
        info = tile_info[f]
        if np.any(np.isnan(info)):
            continue
        if cols is None:
            cols = {
                "prefix": [None] * n,
                "Iter": np.zeros(n),
                "subIter": [None] * n,
                "fullIter": [None] * n,
                "camera": np.full(n, "\x00", dtype="<U1"),
                "ch": np.zeros(n),
                "stack": np.zeros(n),
                "x": np.zeros(n),
                "y": np.zeros(n),
                "z": np.zeros(n),
                "t": np.zeros(n),
                "mappedFilename": [None] * n,
            }
        cols["prefix"][f] = "Scan"
        cols["Iter"][f] = info[0]
        cols["subIter"][f] = ""
        cols["fullIter"][f] = _sprintf("%04d", info[0])
        if specify_cam:
            k = info[5]
            if not cam_patterns:
                # camPatterns is [] without a camera in the channel pattern
                raise ImageListError(
                    "Brace indexing is not supported for variables of this type.",
                    "MATLAB:cellRefFromNonCell",
                )
            if k < 1 or k != int(k):
                raise ImageListError(
                    "Index in position 1 is invalid. Array indices must be positive integers or "
                    "logical values.",
                    "MATLAB:index:outOfBounds",
                )
            cols["camera"][f] = cam_patterns[int(k) - 1]
        else:
            cols["camera"][f] = "A"  # A stands for the camera when names have none
        cols["ch"][f] = info[1]
        cols["stack"][f] = 0.0
        cols["x"][f] = info[2]
        cols["y"][f] = info[3]
        cols["z"][f] = info[4]
        cols["t"][f] = info[0]
        cam = str(cols["camera"][f]) if specify_cam else None
        cols["mappedFilename"][f] = _mapped_name(
            cols["fullIter"][f], cam, info[1], 0.0, info[2], info[3], info[4], info[0]
        )
    out = base.copy()
    if cols is None:
        raise ImageListError(
            "Unrecognized table variable name 'x'.", "MATLAB:table:UnrecognizedVarName"
        )
    for k, v in cols.items():
        out[k] = _texts(v) if isinstance(v, list) else v
    return out


# ------------------------------------------------------------------------------------------
# encoder
# ------------------------------------------------------------------------------------------


def _first_row(path: str, names: Sequence[str]) -> Table:
    """The first data row of a slice list, split at commas (MATLAB's fgetl/strsplit path)."""
    from petakit5d.core.text import str2double

    line = None
    for _ in range(3):
        with open(path, encoding="utf-8", errors="replace", newline="") as f:
            f.readline()
            line = f.readline()
        line = line[:-1] if line.endswith("\n") else (line or None)
        if line:
            break
        time.sleep(5)
    if not line:
        raise ImageListError("First argument must be text.", "MATLAB:strsplit:InvalidStringType")
    fields = line.split(",")
    if len(fields) != len(names):
        raise ImageListError(
            "The VariableNames property must contain one name for each variable in the table.",
            "MATLAB:table:IncorrectNumberOfVarNames",
        )
    cols: dict[str, Any] = {names[0]: _texts([fields[0]])}
    for n, v in zip(names[1:], fields[1:], strict=True):
        cols[n] = np.array([str2double(v)])
    return Table(cols, height=1)


def _fmin(a: np.ndarray, mx: bool = False) -> np.ndarray:
    a = a[~np.isnan(a)] if a.size else a
    if a.size == 0:
        return np.empty(0)
    return np.array([a.max() if mx else a.min()])


def _estimate(
    data_path: str,
    *,
    coord_type: str = "FPGA",
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0"),
    pixel_size=0.108,
    dz=0.2,
    angle=32.45,
    only_coordinates: bool = False,
    x_thresh=0.3,
    skip_missing_file: bool = True,
    only_first_row: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task=2,
    uuid: str = "",
):
    from petakit5d.core.ids import get_uuid
    from petakit5d.io.matfile import StructArray, save_ascii, write_mat

    # ip.Results: fields in alphabetical order (any case)
    pr = {
        "angle": angle,
        "channelPatterns": list(channel_patterns),
        "CoordType": coord_type,
        "cpusPerTask": cpus_per_task,
        "dataPath": data_path,
        "dz": dz,
        "masterCompute": bool(master_compute),
        "OnlyCoordinates": bool(only_coordinates),
        "onlyFirstRow": bool(only_first_row),
        "parseCluster": bool(parse_cluster),
        "pixelSize": pixel_size,
        "skipMissingFile": bool(skip_missing_file),
        "uuid": uuid,
        "xThresh": x_thresh,
    }
    dz = float(np.asarray(dz).reshape(-1)[0])
    x_thresh = float(np.asarray(x_thresh).reshape(-1)[0])
    data_path = simplify_path(data_path) + "/"
    result_path = data_path + "encoder_info/"
    if not os.path.isdir(result_path):
        os.makedirs(result_path, exist_ok=True)
        _group_writable(result_path)
    write_mat(result_path + "/parameters.mat", {"pr": pr})
    _group_writable(result_path + "/parameters.mat")
    result_path = simplify_path(result_path) + "/"

    print("Process slice coordinates...")
    loc_fsns = [e.name for e in matlab_dir(data_path + "/SliceList*csv")]
    tok = matlab_regex(r"(SliceList_.+part\d+).csv")
    prefix = []
    for s in loc_fsns:
        m = tok.search(s)
        if m is None:
            raise ImageListError(
                "Index exceeds the number of array elements. Index must not exceed 0.",
                "MATLAB:badsubscript",
            )
        prefix.append(m.group(1))
    names = None
    if only_first_row:
        names = read_table(f"{data_path}/{prefix[0]}.csv").names
    tabs = []
    for p in prefix:
        loc_fn = f"{data_path}/{p}.csv"
        tabs.append(_first_row(loc_fn, names) if only_first_row else read_table(loc_fn))
    t = Table.vcat(tabs)

    filenames = list(t["Filename"])
    t_fns = unique_text(filenames)
    include = np.zeros(len(t_fns), dtype=bool)
    for c in channel_patterns:
        include |= contains_pattern(t_fns, c)
    if not skip_missing_file:
        for f, fn in enumerate(t_fns):
            if include[f] and not os.path.isfile(data_path + fn):
                include[f] = False
    t_fns = [x for x, k in zip(t_fns, include, strict=True) if k]
    if not t_fns:
        raise ImageListError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )

    # time point and tile indices (MATLAB's first test is always true)
    rx = matlab_regex(r"_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t", ignore_case=True)
    fn_arr = np.array(filenames, dtype=object)
    nz = int(np.sum(fn_arr == t_fns[0]))
    skipped: list[str] = []
    coords: list[dict] = []
    ch_col, cam_col = t["Channel"], t["Is_CamB_"]
    cols = {
        "xfpga": t["XFPGAPosition_um_"],
        "yfpga": t["YFPGAPosition_um_"],
        "zfpga": t["ZFPGAPosition_um_"],
        "xmcs2": t["XMCS2Position_um_"],
        "ymcs2": t["YMCS2Position_um_"],
        "zmcs2": t["ZMCS2Position_um_"],
    }
    for fn in t_fns:
        inds = fn_arr == fn
        if int(inds.sum()) < nz:
            skipped.append(fn)
            continue
        first = int(np.flatnonzero(inds)[0])
        c: dict[str, Any] = {"Filename": fn, "ch": ch_col[first], "cam": cam_col[first]}
        for k, v in cols.items():
            c[k] = v[inds].reshape(-1, 1)
        src = "fpga" if coord_type == "FPGA" else "mcs2"
        for a in "xyz":
            c[a] = c[a + src]
        m = rx.search(fn)
        if m is None:
            raise ImageListError(
                "Dot indexing is not supported for variables of this type.",
                "MATLAB:structRefFromNonStruct",
            )
        c["tind"] = np.array([[float(m["x"]), float(m["y"]), float(m["z"])]])
        c["t"] = float(m["t"])
        coords.append(c)

    # remove duplicate slices in the beginning and end
    nf = len(t_fns) - len(skipped)
    duplicate_mat = np.ones((nf, nz), dtype=bool)
    x_st = np.zeros((nf, 2))
    q = int(mround(nz / 4))
    for f in range(nf):
        t_x = coords[f]["x"].reshape(-1)
        d_x = np.diff(t_x)
        d_inds = ((dz > 0) & (d_x < dz * x_thresh)) | ((dz < 0) & (d_x > dz * x_thresh))
        head = np.flatnonzero(d_inds[:q])
        s = int(head[-1]) + 2 if head.size else 1
        tail = np.flatnonzero(d_inds[nz - q :])
        e = int(tail[0]) + 1 + nz - q if tail.size else nz
        duplicate_mat[f, s - 1 : e] = False
        x_st[f, :] = coords[f]["x"].reshape(-1)[[s - 1, e - 1]]
        coords[f]["duplicate"] = duplicate_mat[f : f + 1, :].copy()

    # the common region to keep across time points
    if dz > 0:
        global_xrange = np.array([[x_st[:, 0].max(), x_st[:, 1].min()]])
    else:
        global_xrange = np.array([[x_st[:, 0].min(), x_st[:, 1].max()]])

    # the global minimum for each tile (same tile indices) after removing the duplicates
    tinds_all = np.concatenate([c["tind"] for c in coords], axis=0)
    encoded = tinds_all[:, 0] * 1e6 + tinds_all[:, 1] * 1e3 + tinds_all[:, 2]
    for u in unique_num(encoded):
        sel = np.flatnonzero(encoded == u)
        keep = ~duplicate_mat[sel, :].T
        mins, maxs = [], []
        for a in "xyz":
            v = np.concatenate([coords[i][a] for i in sel], axis=1)[keep]
            mins.append(_fmin(v))
            maxs.append(_fmin(v, mx=True))
        glbmin = np.concatenate(mins).reshape(1, -1)
        glbmax = np.concatenate(maxs).reshape(1, -1)
        for i in sel:
            coords[i]["glbmin"] = glbmin
            coords[i]["glbmax"] = glbmax

    # actual step size (dz) without the top and bottom 10% of the coordinates
    parts = []
    for f in range(nf):
        d_x = np.diff(coords[f]["x"].reshape(-1))
        s = int(mround(d_x.size * 0.10)) + 1
        e = int(mround(d_x.size * 0.90))
        parts.append(d_x[s - 1 : e])
    single = [coords[f]["x"].size == 1 for f in range(nf)]
    if all(single):
        # diff of one slice is 0x0, its part 1x0: cat(1, ...) is nF x 0 and its mean, along
        # the first dimension, 1x0 (NaN for one tile)
        dz_actual: Any = NAN if nf == 1 else np.zeros((1, 0))
    elif any(single):
        raise ImageListError(
            "Dimensions of arrays being concatenated are not consistent.",
            "MATLAB:catenate:dimensionMismatch",
        )
    else:
        d_all = np.abs(np.concatenate(parts))
        total = 0.0
        for v in d_all.tolist():
            total += v
        dz_actual = total / d_all.size if d_all.size else NAN

    global_duplicate = np.any(duplicate_mat, axis=0).reshape(1, -1)
    fields = (
        "Filename", "ch", "cam", "xfpga", "yfpga", "zfpga", "xmcs2", "ymcs2", "zmcs2",
        "x", "y", "z", "tind", "t", "duplicate", "glbmin", "glbmax",
    )  # fmt: skip
    coords_struct = StructArray.from_list(
        [{k: c[k] for k in fields} for c in coords], fields=fields
    )
    info_fn = f"{result_path}/coordinates_info.mat"
    tmp_fn = f"{result_path}/coordinates_info_{get_uuid()}.mat"
    write_mat(
        tmp_fn,
        {
            "coords_struct": coords_struct,
            "skipped_fns": np.array(skipped, dtype=object).reshape(-1, 1) if skipped else [],
            "CoordType": coord_type,
            "global_xrange": global_xrange,
            "x_st": x_st,
            "dz_actual": dz_actual,
            "duplicate_mat": duplicate_mat,
            "global_duplicate": global_duplicate,
        },
    )
    _group_writable(tmp_fn)
    # sprintf('%sdz_actual_%d.txt', resultPath, dz_actual): an empty step prints as nothing
    # (core's sprintf stops at an empty argument; reported to F02)
    step = _sprintf("%d", dz_actual) if np.size(dz_actual) else ""
    txt = f"{result_path}dz_actual_{step}.txt"
    save_ascii(txt, dz_actual)
    _group_writable(txt)
    os.replace(tmp_fn, info_fn)
    return dz_actual, coords


def estimate_actual_step_size_from_encoder(
    data_path: str,
    *,
    coord_type: str = "FPGA",
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0"),
    pixel_size=0.108,
    dz=0.2,
    angle=32.45,
    only_coordinates: bool = False,
    x_thresh=0.3,
    skip_missing_file: bool = True,
    only_first_row: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task=2,
    uuid: str = "",
) -> float:
    """``XR_estimate_actual_step_size_from_encoder``: the mean step between the slices of each
    tile (the middle 80%) from the ``SliceList*csv`` files; writes ``encoder_info/``
    (parameters.mat, coordinates_info.mat with the per-tile coordinates and duplicate slices,
    ``dz_actual_<dz>.txt``). With ``only_first_row`` only the first slice of each file is read
    (the step is then NaN)."""
    dz_actual, _ = _estimate(
        data_path,
        coord_type=coord_type,
        channel_patterns=channel_patterns,
        pixel_size=pixel_size,
        dz=dz,
        angle=angle,
        only_coordinates=only_coordinates,
        x_thresh=x_thresh,
        skip_missing_file=skip_missing_file,
        only_first_row=only_first_row,
        parse_cluster=parse_cluster,
        master_compute=master_compute,
        cpus_per_task=cpus_per_task,
        uuid=uuid,
    )
    return dz_actual


def _unique_rows(t: Table) -> Table:
    """``unique(t, 'rows')``: the distinct rows, sorted by the columns in order."""
    from petakit5d.core.paths import matlab_sort_key

    keys = []
    for n in t.names:
        c = t[n]
        if c.dtype == object:
            keys.append([(0, matlab_sort_key(s)) if s is not None else (1, ()) for s in c])
        else:
            keys.append([(1, 0.0) if v != v else (0, v) for v in c.tolist()])
    rows = [tuple(k[i] for k in keys) for i in range(t.height)]
    order = sorted(range(t.height), key=lambda i: rows[i])
    keep = []
    for i in order:
        r = rows[i]
        nan = any(isinstance(x, tuple) and x[0] == 1 and x[1] == 0.0 for x in r)
        if keep and not nan and rows[keep[-1]] == r:
            continue
        keep.append(i)
    return t.rows(np.asarray(keep, dtype=np.intp))


def generate_image_list_from_encoder(
    data_path: str, dz, channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0")
) -> str:
    """``stitch_generate_imagelist_from_encoder``: ``ImageList_from_encoder.csv`` with the
    first MCS2 position of each tile (stage and objective columns alike). Returns its path."""
    from petakit5d.core.ids import get_uuid

    _, coords = _estimate(
        data_path,
        dz=dz,
        channel_patterns=channel_patterns,
        only_first_row=True,
        skip_missing_file=True,
    )
    fpth = [f"{data_path}/{c['Filename']}" for c in coords]
    xyz = [[c["xmcs2"][0, 0], c["ymcs2"][0, 0], c["zmcs2"][0, 0]] for c in coords]
    t = _stage_table(fpth, [c["Filename"] for c in coords], np.array(xyz).reshape(-1, 3))
    t = _unique_rows(t)
    return _publish(t, data_path, "ImageList_from_encoder", get_uuid())


def generate_image_list_from_encoder_ism(data_path: str, position_type: str = "MCS2") -> str:
    """``stitch_generate_imagelist_from_encoder_ISM``: ``ImageList_ISM_from_encoder.csv`` from
    the first slice (``DataImageIndex`` 0) of every ``SliceList_*.csv`` (``data_path`` ends
    with a separator, as in MATLAB). Returns its path."""
    names = [e.name for e in matlab_dir(f"{data_path}SliceList_*.csv")]
    tabs = []
    for n in names:
        t = read_table(f"{data_path}{n}")
        tabs.append(t[t["DataImageIndex"] == 0])
    tin = Table.vcat(tabs)
    n = tin.height
    fns = [s.replace("_0000ISMslice", "") for s in tin["Filename"]] if n else []
    src = {"mcs2": "MCS2", "target": "Target"}.get(position_type.lower())
    xyz = np.zeros((n, 3))
    if src is not None:
        for j, a in enumerate("XYZ"):
            xyz[:, j] = tin[f"{a}{src}Position_um_"]
    t = _stage_table(fns, fns, xyz)
    fnout = f"{data_path}/ImageList_ISM_from_encoder.csv"
    write_table(t, fnout)
    return fnout


# ------------------------------------------------------------------------------------------
# sqlite
# ------------------------------------------------------------------------------------------


def _fetch(dbfile: str) -> tuple[list[str], list[tuple]]:
    conn = sqlite3.connect(f"file:{dbfile}?mode=ro", uri=True)
    try:
        cur = conn.execute("SELECT * FROM ImageListView_NoNull")
        cols = [d[0] for d in cur.description]
        rows = cur.fetchall()
    finally:
        conn.close()
    return cols, rows


def generate_image_list_from_sqlite(
    data_path: str,
    *,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    map_tilename: bool = True,
) -> str:
    """``stitch_generate_imagelist_from_sqlite``: ``ImageList_from_sqlite.csv`` from the
    acquisition databases of the folder, with the target positions; with ``map_tilename`` the
    names mapped to the standard grammar (channels numbered per camera by channel and laser).
    Returns its path."""
    from petakit5d.core.ids import get_uuid

    fsns = [e.name for e in matlab_dir(f"{data_path}/*.sqlite3")]
    fname: list[str] = []
    pos: list[list[float]] = []
    info: list[list[float]] = []
    keys = ("Camera", "Channel", "Laser_nm", "File_Stack_Index", "Elapsed_ms", "Absolute_ms")
    keys += ("Tile_index_X", "Tile_index_Y", "Tile_index_Z", "Tile_index_T")
    for f in fsns:
        dbfile = f"{data_path}/{f}"
        cols, rows = _fetch(dbfile)
        if not rows:
            _warn(f"{dbfile} is empty, skip it!")
            continue
        ci = {c: i for i, c in enumerate(cols)}
        for r in rows:
            fname.append(str(r[ci["Filename"]]))
            pos.append([float(r[ci[f"Mean_{a}_Target_Position_um"]]) for a in "XYZ"])
            if map_tilename:
                info.append([float(r[ci[k]]) for k in keys])
    fsn = list(fname)
    full = [f"{data_path}/{x}" for x in fsn]
    include = np.zeros(len(fsn), dtype=bool)
    for c in channel_patterns:
        include |= contains_pattern(full, c)
    sel = np.flatnonzero(include) if not include.all() else np.arange(len(fsn))
    fname = [fname[i] for i in sel]
    xyz = np.array([pos[i] for i in sel]).reshape(-1, 3)
    mapped = None
    if map_tilename:
        ti = np.array([info[i] for i in np.flatnonzero(include)]).reshape(-1, len(keys))
        ti[:, 1] = _map_channels(ti[:, 0:3])
        cam_strs = ("A", "B", "C", "D")
        mapped = [
            _sprintf(
                "Scan_Iter_%04d_Cam%s_ch%d_CAM1_stack%04d_%03dnm_%06dmsec_%08dmsecAbs_%03dx_%03dy_%03dz_%04dt.tif",
                r[9],
                cam_strs[int(r[0])],
                r[1],
                r[3],
                r[2],
                r[4],
                r[5],
                r[6],
                r[7],
                r[8],
                r[9],
            )
            for r in ti
        ]
    t = _stage_table(fname, fname, xyz, mapped)
    # deduplicate tile names (unique keeps the first of each, sorted)
    first: dict[str, int] = {}
    for i, s in enumerate(t["Filename"]):
        first.setdefault(s, i)
    order = [first[s] for s in unique_text(first)]
    t = t.rows(np.asarray(order, dtype=np.intp))
    return _publish(t, data_path, "ImageList_from_sqlite", get_uuid())


def _map_channels(ccl: np.ndarray) -> np.ndarray:
    """Channels numbered 0, 1, ... per camera over the sorted unique (camera, channel, laser)."""
    uniq = sorted({tuple(r) for r in ccl.tolist()})
    mapped = ccl[:, 1].copy()
    for cam in sorted({u[0] for u in uniq}):
        group = [u for u in uniq if u[0] == cam]
        for k, u in enumerate(group):
            mapped[np.all(ccl == np.array(u), axis=1)] = k
    return mapped


def generate_image_list_from_sqlite_old_format(data_path: str) -> str:
    """``stitch_generate_imagelist_from_sqlite_old_format``: columns 2 and 14 to 19 of the
    view (name, MCS2 and target positions); ``data_path`` ends with a separator. Returns the
    file's path."""
    from petakit5d.core.ids import get_uuid

    fsns = [e.name for e in matlab_dir(f"{data_path}*.sqlite3")]
    fname: list[str] = []
    vals: list[list[float]] = []
    for f in fsns:
        dbfile = f"{data_path}{f}"
        _, rows = _fetch(dbfile)
        if not rows:
            _warn(f"{dbfile} is empty, skip it!")
            continue
        for r in rows:
            fname.append(str(r[1]))
            vals.append([float(v) for v in r[13:19]])
    v = np.array(vals).reshape(-1, 6)
    cols: dict[str, Any] = {"Filepath": _texts(fname), "Filename": _texts(fname)}
    for j, name in enumerate(_STAGE):
        cols[name] = v[:, j].copy()
    t = Table(cols, height=len(fname))
    return _publish(t, data_path, "ImageList_from_sqlite", get_uuid())


# ------------------------------------------------------------------------------------------
# XR_generate_image_list_wrapper
# ------------------------------------------------------------------------------------------


def _matlab_in(kwargs: dict[str, Any]) -> dict[str, Any]:
    """MATLAB-name layer: ``overlapSize`` is ``[y x z]`` in pixels and ``[x y z]`` in um
    (``overlapSizeType``, default pixel); the Python value is ``(z, y, x)`` for both units."""
    from petakit5d.params.kinds import KINDS

    ov = kwargs.get("overlap_size")
    if ov is None or np.ndim(ov) == 0:
        return kwargs
    unit = str(kwargs.get("overlap_size_type", "pixel")).lower()
    kind = KINDS["per_axis_xyz" if unit == "um" else "size_yxz"]
    return {**kwargs, "overlap_size": kind.from_matlab(ov)}


@register(
    "XR_generate_image_list_wrapper",
    Param(
        "dataPaths",
        "@(x) ischar(x) || iscell(x)",
        "Dataset folders (one, or several with the same settings).",
    ),
    Param(
        "generationMethod",
        "@(x) ischar(x)",
        "'encoder', 'sqlite', 'tile_position' or 'tile_list': from encoder positions, the "
        "acquisition database, the overlap between neighboring tiles, or the given tile "
        "names, indices and interval.",
    ),
    Param("channelPatterns", "@iscell", "Channel identifiers of the included channels."),
    Param("mapTilename", "@islogical", "Map tile names to standard ones (sqlite only)."),
    Param(
        "tilePatterns",
        "@iscell",
        "Patterns for time, channel, x, y and z in the tile names ([a-zA-Z]*[0-9]* or "
        "[0-9]*[a-zA-Z]*).",
    ),
    Param("tileFilenames", "@iscell", "Tile file names (tile_list)."),
    Param(
        "tileIndices",
        "@(x) isempty(x) || (isnumeric(x) && size(x, 2) == 5)",
        "Tile indices of the tile files, #file x 5, order t, c, x, y, z (tile_list).",
        convert="same",
    ),
    Param(
        "tileInterval",
        "@(x) isempty(x) || (isnumeric(x) && size(x, 1) == 1 && size(x, 2) == 3)",
        "Interval between neighboring tiles in um, order x, y, z (tile_list).",
        convert="same",
    ),
    Param("DS", "@islogical", "Data is in deskewed space."),
    Param("DSR", "@islogical", "Data is in deskew/rotated space (with stage coordinates)."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Scan interval in um."),
    Param("skewAngle", "@isnumeric", "Skew angle (in degrees) of the stage."),
    Param(
        "axisOrder",
        "@ischar",
        "Axis order of the coordinates in the image list ('-yxz': -y maps to x, x to y, z to z).",
        convert="same",
    ),
    Param(
        "dataOrder",
        "@ischar",
        "Axis order of the data ('z,y,x' in Python, MATLAB's 'y,x,z').",
        convert="data_order",
    ),
    Param("objectiveScan", "@islogical", "Objective scan."),
    Param(
        "IOScan",
        "@islogical",
        "Inverted objective scan (the scan with the stage coordinates, DSR space).",
    ),
    Param("zarrFile", "@islogical", "Zarr files as input."),
    Param(
        "overlapSize",
        "@isnumeric",
        "Overlap between tiles in pixels or um ((z, y, x) in Python; MATLAB's [y x z] in "
        "pixels and [x y z] in um, converted by the unit).",
    ),
    Param(
        "overlapSizeType",
        "@(x) ischar(x) && ismember(lower(x), {'pixel', 'um'})",
        "Unit of the overlap size: 'pixel' or 'um'.",
    ),
    Param("uuid", "@ischar", "uuid in the temporary result paths."),
    task="S01",
    matlab_in=_matlab_in,
)
def generate_image_list(
    data_paths,
    generation_method,
    *,
    channel_patterns=("CamA_ch0", "CamA_ch1", "CamB_ch0"),
    map_tilename=True,
    tile_patterns=("0000t", "ch0", "000x", "000y", "000z"),
    tile_filenames=(),
    tile_indices=None,
    tile_interval=None,
    ds=False,
    dsr=False,
    xy_pixel_size=0.108,
    dz=0.2,
    skew_angle=32.45,
    axis_order="x,y,z",
    data_order="z,y,x",
    objective_scan=False,
    io_scan=False,
    zarr_file=False,
    overlap_size=(10, 10, 10),
    overlap_size_type="pixel",
    uuid="",
):
    """Generate the image list (stage coordinates of the tiles) of datasets as
    ``ImageList_from_<method>.csv`` in their folder (the first folder for tile_position and
    tile_list)."""
    if isinstance(data_paths, str):
        data_paths = [data_paths]
    data_paths = list(data_paths)
    method = generation_method.lower()
    if method == "encoder":
        for d in data_paths:
            generate_image_list_from_encoder(d, dz, channel_patterns)
    elif method == "sqlite":
        for d in data_paths:
            generate_image_list_from_sqlite(
                d, channel_patterns=channel_patterns, map_tilename=map_tilename
            )
    elif method == "tile_position":
        generate_image_list_from_tile_positions(
            data_paths,
            channel_patterns=channel_patterns,
            tile_patterns=tile_patterns,
            xy_pixel_size=xy_pixel_size,
            dz=dz,
            skew_angle=skew_angle,
            axis_order=axis_order,
            data_order=data_order,
            ds=ds,
            dsr=dsr,
            objective_scan=objective_scan,
            io_scan=io_scan,
            zarr_file=zarr_file,
            overlap_size=overlap_size,
            overlap_size_type=overlap_size_type,
            uuid=uuid,
        )
    elif method == "tile_list":
        if tile_indices is None or tile_interval is None:
            raise ImageListError(
                "The value of 'tileIndices' is invalid. It must satisfy the function: "
                "@(x)isnumeric(x)&&size(x,2)==5.",
                "MATLAB:InputParser:ArgumentFailedValidation",
            )
        generate_image_list_from_tile_list(
            data_paths, tile_filenames, tile_indices, tile_interval, uuid=uuid
        )
    else:
        raise ImageListError(
            "Unknown gerneration method, it must be choosen from: encoder, sqlite, or "
            "tile_position."
        )
