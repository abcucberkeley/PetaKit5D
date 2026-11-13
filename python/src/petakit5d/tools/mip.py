"""Maximum intensity projections (MATLAB's ``tools/MIP``): of arrays (``saveMIP_frame``), tiff
files (``saveMIP_tiff``), zarr files (``saveMIP_zarr``), datasets (``XR_MIP_wrapper``) and large
zarr files in batches with MIP slabs (``XR_MIP_zarr``, ``MIP_block``), and the MIPs that block
tasks of large-file processing save with their flags (``MIP_zarr_block_mips``,
``MIP_zarr_save_block_mips``, ``MIP_zarr_merge_block_mips``, ``MIP_zarr_merge_memory``,
``MIP_zarr_batch_size``).

Arrays are ``(z, y, x)`` C order; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive;
sizes, batch, block and pool sizes and the axis flags are in ``(z, y, x)`` (a z MIP is
``axis=(1, 0, 0)``); an axis is a NumPy axis. The MATLAB-name layer converts MATLAB's values.

Each MIP is MATLAB's ``max`` (``range_max_3d``): NaN are ignored (a line of NaN gives its
first NaN), and of equal values (0 and -0) the first along the axis is kept. MIPs of slabs and
of blocks are merged in increasing position along the axis, which gives the same bits.

Files keep MATLAB's names, contents and orientation (the user's decision of 2026-11-12):

- ``<name>_MIP_y.tif``, ``_MIP_x.tif``, ``_MIP_z.tif`` of ``save_mip_tiff`` and
  ``save_mip_zarr`` are MATLAB's squeezed 2D images: along z ``ny x nx``, along y MATLAB's
  ``nx x nz`` and along x ``ny x nz`` (the transposes of the natural ``(nz, nx)`` and
  ``(nz, ny)`` projections: one 2D transpose when a MIP is written). ``save_mip_frame`` writes
  MATLAB's unsqueezed projections, ``(nz, 1, nx)`` and ``(nz, ny, 1)`` as they are.
- the MIP slabs ``<MIPs>/<name>_{y,x,z}.zarr`` (one element per batch along their axis) are
  native zarr arrays, as every zarr the port writes;
- the flags of block tasks with the MIPs of their blocks are MATLAB's ``.mat`` files with
  MATLAB's values (boxes ``[ymin xmin zmin ymax xmax zmax]``, sizes ``[y x z]``, MIPs as
  MATLAB's ``(y, x, z)`` arrays), so MATLAB's and the port's flags read alike.

A MIP file is written under a temporary name and renamed into place, as the job framework
counts a MIP task as done once its MIP exists; ``mip_zarr`` publishes its slabs before its MIPs,
and the MIP of the last axis, the wrapper's flag, last.

Not ported: the resampling of ``MIP_block`` before pooling (``pool_size`` with 9 values and a
factor above 1, MATLAB's ``imresize3``), and the ``mean_*`` methods of ``project3DImageto2D``
(``nanmean``); both raise ``NotImplementedError``.
"""

from __future__ import annotations

import math
import os
import shutil
import time
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "axis_max",
    "block_mips",
    "merge_block_mips",
    "merge_memory",
    "mip_block",
    "mip_dataset",
    "mip_zarr",
    "mip_zarr_batch_size",
    "project_3d_image_to_2d",
    "save_block_mips",
    "save_mip_frame",
    "save_mip_tiff",
    "save_mip_zarr",
]

# MATLAB's order of the MIP files (y, x, z), their NumPy axes and their place in (z, y, x)
_LABELS = ("y", "x", "z")
_AXES = (-2, -1, -3)
_ZYX = (1, 2, 0)
# the classes range_max_3d takes (MATLAB's max for the others: integers, no NaN)
_RANGE_DTYPES = frozenset(np.dtype(t) for t in ("uint8", "uint16", "float32", "float64"))
# MATLAB's whos: bytes of a cell element's header (R2024b, 64-bit)
_CELL_HEADER_BYTES = 112
# MIP images from this size are turned to MATLAB's orientation in threads
_THREADED_BYTES = 1 << 22
# GB a MIP batch task asks for: this many times its batch. MATLAB asks for 2.5 times in its
# compiled mode and 4 times otherwise (2.5 + 1.5 * ~mccMode); the port's batch task holds one
# batch and its three MIPs: peak RSS 17.1 GB for two 2048^3 uint16 batches (16 GiB each) in one
# mip_zarr process, 1.07 times a batch (c24 node, 24 threads, 2026-11-12), plus a margin.
BATCH_MEMORY_FACTOR = 1.2


class MipError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _flags(axis: Any) -> list[bool]:
    """``axis(i) ~= 0`` for MATLAB's y, x, z (the flags are ``(z, y, x)``)."""
    a = np.asarray(axis, dtype=np.float64).ravel()
    if a.size != 3:
        raise MipError(
            f"axis needs 3 flags (z, y, x), got {a.size}: Index exceeds the number of array "
            f"elements. Index must not exceed {a.size}.",
            "MATLAB:badsubscript",
        )
    return [bool(a[_ZYX[i]] != 0) for i in range(3)]


def _np_axis(axis: int) -> int:
    """a NumPy axis of a (z, y, x) volume as 0, 1 or 2"""
    if isinstance(axis, (bool, np.bool_)) or not isinstance(axis, (int, np.integer)):
        raise ValueError(f"axis must be a NumPy axis of a (z, y, x) volume, not {axis!r}")
    a = int(axis)
    if not -3 <= a <= 2:
        raise ValueError(f"axis must be a NumPy axis of a (z, y, x) volume, not {axis!r}")
    return a % 3


def _ints(v: Any) -> list[int]:
    return [int(x) for x in np.asarray(v).ravel()]


def _num(x: float) -> int | float:
    """a MATLAB double as a Python int where it is integral"""
    return int(x) if float(x).is_integer() else float(x)


def _vol3(a: np.ndarray) -> np.ndarray:
    """an image (y, x) as the volume (1, y, x)"""
    return a[np.newaxis] if a.ndim == 2 else a


def _cast(a: np.ndarray, dtype: str) -> np.ndarray:
    """``cast(a, dtype)`` (round half away from zero, saturate, NaN to 0); none for its class"""
    from petakit5d.core.arrays import matlab_class
    from petakit5d.core.dtypes import cast

    if matlab_class(a.dtype) == dtype:
        return a
    return np.ascontiguousarray(cast(a, dtype))


def axis_max(vol: np.ndarray, axis: int, *, num_threads: int | None = None) -> np.ndarray:
    """MATLAB's ``max(vol, [], dim)`` along a NumPy axis (-3 z, -2 y, -1 x), the axis kept with
    size 1: ``range_max_3d(vol, axis, [0])`` for uint8, uint16, single and double (NaN ignored,
    the first of equal values), NumPy's max for the other classes (integers and logical, where
    no rule of MATLAB's matters). An image ``(y, x)`` along z is the image itself."""
    from petakit5d.kernels.mip import range_max_3d

    a = np.asarray(vol)
    if a.dtype in _RANGE_DTYPES and a.size:
        return range_max_3d(a, axis, [0], num_threads=num_threads)
    if a.ndim < 3 and axis == -3:
        return a.copy()
    return np.max(a, axis=axis, keepdims=True)


def _matlab_image(mip: np.ndarray) -> np.ndarray:
    """The C-ordered (rows, columns) image that MATLAB's ``squeeze`` of a projection is in
    MATLAB's orientation (``writetiff`` writes its rows as rows): a y MIP ``(nz, 1, nx)`` is
    MATLAB's ``nx x nz``, an x MIP ``(nz, ny, 1)`` ``ny x nz``, a z MIP ``(1, ny, nx)`` and an
    image as they are; ``1 x 1 x n`` gives ``n x 1``, as MATLAB squeezes it."""
    from petakit5d.core.arrays import size_to_matlab

    s = list(size_to_matlab(mip.shape))
    if len(s) > 2:
        s = [n for n in s if n != 1]
        s += [1] * (2 - len(s))
    rows, cols = s
    if mip.ndim < 2:
        lin = mip
    else:
        # MATLAB's linear order (y fastest, then x, then z) is the C order of this view
        lin = np.swapaxes(mip, -1, -2)
    view = lin.reshape(cols, rows).T
    if view.flags.c_contiguous or view.nbytes < _THREADED_BYTES:
        return np.ascontiguousarray(view)
    return _copy_rows(view)


def _copy_rows(view: np.ndarray) -> np.ndarray:
    """a C-ordered copy of a strided 2D view (the transpose of a y or x MIP), its row blocks
    copied in threads (NumPy's copy releases the GIL)"""
    from concurrent.futures import ThreadPoolExecutor

    from petakit5d.core.threads import resolve_threads

    out = np.empty(view.shape, dtype=view.dtype)
    n = min(resolve_threads(None), view.shape[0])
    if n < 2:
        out[...] = view
        return out
    edges = np.linspace(0, view.shape[0], 4 * n + 1).astype(int)

    def copy(k: int) -> None:
        out[edges[k] : edges[k + 1]] = view[edges[k] : edges[k + 1]]

    with ThreadPoolExecutor(max_workers=n) as ex:
        list(ex.map(copy, range(4 * n)))
    return out


def _publish_tiff(img: np.ndarray, final: str, uuid: str) -> None:
    """``writetiff`` to ``<final without .tif>_<uuid>.tif``, then ``movefile`` onto final"""
    from petakit5d.io import write_tiff

    tmp = f"{final[:-4]}_{uuid}.tif"
    write_tiff(img, tmp)
    os.replace(tmp, final)


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _toc(t0: float) -> None:
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)


# ---------------------------------------------------------------------------------------------
# projections of arrays and files
# ---------------------------------------------------------------------------------------------


def project_3d_image_to_2d(frame: np.ndarray, method: str) -> np.ndarray:
    """``project3DImageto2D``: a 2D view of a ``(z, y, x)`` volume (or a ``(y, x)`` image), with
    the projected axis removed: ``'mip_xy'`` (along z) ``(ny, nx)``, ``'mip_xz'`` (along y)
    ``(nz, nx)``, ``'mip_yz'`` (along x) ``(nz, ny)``; ``'central_xy'``, ``'central_xz'``,
    ``'central_yz'`` the central slice (index ``floor((n + 1) / 2)`` 1-based); any other method
    the frame itself. These are the natural NumPy projections (MATLAB's of xz and yz are their
    transposes, ``nx x nz`` and ``ny x nz``). The ``mean_*`` methods (``nanmean``) are not
    ported."""
    a = np.asarray(frame)
    m = method.lower()
    v = _vol3(a)
    if m.startswith("mean_") and m in ("mean_xy", "mean_yz", "mean_xz"):
        raise NotImplementedError(f"project3DImageto2D's {method!r} (nanmean) is not ported")
    if m == "central_xy":
        return v[(v.shape[0] + 1) // 2 - 1]
    if m == "central_yz":
        return v[:, :, (v.shape[2] + 1) // 2 - 1]
    if m == "central_xz":
        return v[:, (v.shape[1] + 1) // 2 - 1, :]
    if m == "mip_xy":
        return axis_max(v, -3)[0]
    if m == "mip_yz":
        return axis_max(v, -1)[:, :, 0]
    if m == "mip_xz":
        return axis_max(v, -2)[:, 0, :]
    return a


def save_mip_frame(
    im: np.ndarray, mip_fullname: str, *, dtype: str = "uint16", axis: Sequence[float] = (1, 0, 0)
) -> None:
    """``saveMIP_frame``: the MIPs of a ``(z, y, x)`` array along the axes flagged in ``axis``
    (``(z, y, x)``), as ``<mip_fullname without its last 10 characters>_MIP_{y,x,z}.tif``.

    As MATLAB: nothing if ``mip_fullname`` exists or no axis is flagged; the array is cast to
    ``dtype`` first; the projections are not squeezed (``(nz, 1, nx)`` along y, ``(nz, ny, 1)``
    along x, MATLAB's ``1 x nx x nz`` and ``ny x 1 x nz`` tiffs); each is written under a
    temporary name and moved into place. Its folder is made if it is missing.
    """
    from petakit5d.core.ids import get_uuid

    if _exist(mip_fullname):
        print(f"The mip file {mip_fullname} already exists, skip it!", flush=True)
        return
    flags = _flags(axis)
    if not any(flags):
        return
    mip_path = os.path.dirname(mip_fullname)
    if mip_path and not os.path.isdir(mip_path):
        os.makedirs(mip_path, exist_ok=True)
        if os.name != "nt":
            from petakit5d.io.tiff import group_writable

            group_writable(mip_path)
    a = _cast(np.asarray(im), dtype)
    uuid = get_uuid()
    for i in range(3):
        if flags[i]:
            final = f"{mip_fullname[:-10]}_MIP_{_LABELS[i]}.tif"
            _publish_tiff(axis_max(a, _AXES[i]), final, uuid)


def _crop_yx(im: np.ndarray, box: Sequence[int]) -> np.ndarray:
    """crop3d with the box's y and x (z is read already)"""
    from petakit5d.core.bbox import crop

    if im.ndim == 2:
        return crop(im, (box[1], box[2], box[4], box[5]))
    return crop(im, (0, box[1], box[2], im.shape[0], box[4], box[5]))


@register(
    "saveMIP_tiff",
    Param("frameFullname", "@ischar", "Tiff file."),
    Param("MIPFullname", "@ischar", "MIP file name (<name>_MIP_z.tif); skipped if it exists."),
    Param("dtype", "@ischar", "Class of the MIPs."),
    Param("axis", "@isnumeric", "Axes (z, y, x) to project along (nonzero).", convert="flags_yxz"),
    Param(
        "inputBbox",
        "@(x) isempty(x) || isvector(x)",
        "Region (z0, y0, x0, z1, y1, x1) to project (None: all).",
        convert="bbox",
    ),
    task="T02",
)
def save_mip_tiff(
    frame_fullname: str,
    mip_fullname: str,
    *,
    dtype: str = "uint16",
    axis: Sequence[float] = (1, 0, 0),
    input_bbox: Sequence[int] | None = None,
) -> None:
    """The MIPs of a tiff file (or of its box) along the flagged axes, cast to ``dtype``, as
    ``<name>_MIP_{y,x,z}.tif`` in MATLAB's orientation (``saveMIP_tiff``); nothing if
    ``mip_fullname`` exists. Each MIP is written under a temporary name and moved into place."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import read_tiff

    if _exist(mip_fullname):
        print(f"The mip file {mip_fullname} already exists, skip it!", flush=True)
        return
    flags = _flags(axis)
    uuid = get_uuid()
    if input_bbox is None or len(input_bbox) == 0:
        im = read_tiff(frame_fullname)
    else:
        box = [int(v) for v in input_bbox]
        im = _crop_yx(read_tiff(frame_fullname, (box[0], box[3])), box)
    for i in range(3):
        if not flags[i]:
            continue
        mip = _cast(axis_max(im, _AXES[i]), dtype)
        mip_fullname = f"{mip_fullname[:-10]}_MIP_{_LABELS[i]}.tif"
        _publish_tiff(_matlab_image(mip), mip_fullname, uuid)


# ---------------------------------------------------------------------------------------------
# MIPs of zarr files
# ---------------------------------------------------------------------------------------------


def _slab_starts(zarr_fullname: str, box: Sequence[int] | None) -> tuple[list[int], list[int]]:
    """the starts of saveMIP_zarr's z slabs (0-based) and the box they cover: whole chunk
    layers, at least 64 chunks per slab (``ceil(64 / chunks per layer)`` layers)"""
    from petakit5d.io import get_zarr_info

    info = get_zarr_info(zarr_fullname)
    b = list(box) if box is not None and len(box) else [0, 0, 0, *info.size]
    if len(info.size) != 3 or len(b) != 6:
        return [], b
    cz, cy, cx = (int(v) for v in info.block_size)
    layer = (-(-b[4] // cy) - b[1] // cy) * (-(-b[5] // cx) - b[2] // cx)
    step = math.ceil(64 / layer) * cz
    return [b[0], *range((b[0] // step + 1) * step, b[3], step)], b


def _slab_mips(zarr_fullname: str, flags: list[bool], box: Sequence[int] | None) -> list[Any]:
    """saveMIP_zarr_slabs: the MIPs (z, y, x) of the image or its box, read in z slabs with a
    running max along z and the y and x MIPs of the slabs placed along z (one read for files
    that aren't 3D or can't be described)"""
    from petakit5d.io import ReadBuffer, read_zarr

    zs: list[int] = []
    b: list[int] = []
    try:
        zs, b = _slab_starts(zarr_fullname, box)
    except Exception as err:  # (MATLAB: disp(ME))
        print(err, flush=True)
    mips: list[Any] = [None, None, None]
    if len(zs) <= 1:
        im = read_zarr(zarr_fullname, input_bbox=box)
        for i in range(3):
            if flags[i]:
                mips[i] = axis_max(im, _AXES[i])
        return mips
    ze = [*zs[1:], b[3]]
    nz = b[3] - b[0]
    # (one buffer for the slabs: its pages are touched once)
    buf = ReadBuffer()
    for s, (za, zb) in enumerate(zip(zs, ze, strict=True)):
        box_s = (za, b[1], b[2], zb, b[4], b[5])
        im = _vol3(read_zarr(zarr_fullname, input_bbox=box_s, out=buf))
        for i in range(3):
            if not flags[i]:
                continue
            mip = axis_max(im, _AXES[i])
            if i == 2:
                # max(MIPs{3}, mip): the earlier slab first (for integers in place, where no
                # rule of NaN or signed zeros applies)
                if s == 0:
                    mips[i] = mip
                elif mip.dtype.kind in "biu":
                    np.maximum(mips[i], mip, out=mips[i])
                else:
                    mips[i] = axis_max(np.concatenate((mips[i], mip)), -3)
            else:
                if s == 0:
                    mips[i] = np.zeros((nz, *mip.shape[1:]), dtype=mip.dtype)
                mips[i][za - b[0] : zb - b[0]] = mip
        del im
    return mips


def _axis_mip_blockwise(zarr_fullname: str, i: int, box: Sequence[int] | None) -> np.ndarray:
    """saveMIP_zarr_axis: the MIP along one axis in blocks of up to ten chunks (5000 voxels)
    across it, for images that can't be read whole"""
    from petakit5d.io import ReadBuffer, get_zarr_info, read_zarr

    if box is not None and len(box):
        raise MipError("Batch MIP is not supported with input BBox yet!")
    info = get_zarr_info(zarr_fullname)
    size = [int(v) for v in info.size]
    chunks = [int(v) for v in info.block_size]
    if len(size) == 2:
        size, chunks = [1, *size], [1, *chunks]
    p = _ZYX[i]
    blk = [min(5000, c * 10) for c in chunks]
    blk[p] = size[p]
    out_shape = list(size)
    out_shape[p] = 1
    out = None
    buf = ReadBuffer()
    for z0 in range(0, size[0], blk[0]):
        for y0 in range(0, size[1], blk[1]):
            for x0 in range(0, size[2], blk[2]):
                bb = (z0, y0, x0, min(z0 + blk[0], size[0]), min(y0 + blk[1], size[1]))
                bb = (*bb, min(x0 + blk[2], size[2]))
                mip = axis_max(_vol3(read_zarr(zarr_fullname, input_bbox=bb, out=buf)), _AXES[i])
                if out is None:
                    out = np.zeros(out_shape, dtype=mip.dtype)
                sl = [slice(z0, bb[3]), slice(y0, bb[4]), slice(x0, bb[5])]
                sl[p] = slice(0, 1)
                out[tuple(sl)] = mip
    return out


@register(
    "saveMIP_zarr",
    Param("zarrFullname", "@(x) ischar(x)", "Zarr file."),
    Param("MIPFullname", "@(x) ischar(x)", "MIP file name (<name>_MIP_z.tif)."),
    Param("dtype", "@(x) ischar(x) || isstring(x)", "Class of the MIPs (empty: the zarr's class)."),
    Param(
        "axis",
        "@(x) isvector(x) || numel(x) == 3",
        "Axes (z, y, x) to project along (nonzero).",
        convert="flags_yxz",
    ),
    Param(
        "inputBbox",
        "@(x) isempty(x) || isvector(x)",
        "Region (z0, y0, x0, z1, y1, x1) to project (None: all).",
        convert="bbox",
    ),
    task="T02",
)
def save_mip_zarr(
    zarr_fullname: str,
    mip_fullname: str,
    *,
    dtype: str = "",
    axis: Sequence[float] = (1, 0, 0),
    input_bbox: Sequence[int] | None = None,
) -> None:
    """The MIPs of a zarr file (or of its box) along the flagged axes, cast to ``dtype`` (the
    zarr's class if empty), as ``<name>_MIP_{y,x,z}.tif`` in MATLAB's orientation
    (``saveMIP_zarr``). The image is read in z slabs of whole chunk layers (at least 64 chunks
    a read) with a running max along z, so it is not held whole. Each MIP is written under a
    temporary name and moved into place."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import get_image_dtype

    if not dtype:
        dtype = get_image_dtype(zarr_fullname)
    flags = _flags(axis)
    box = None if input_bbox is None or len(input_bbox) == 0 else [int(v) for v in input_bbox]
    uuid = get_uuid()
    try:
        mips = _slab_mips(zarr_fullname, flags, box)
    except Exception as err:  # (MATLAB: disp(ME_1))
        print(err, flush=True)
        mips = [None, None, None]
    for i in range(3):
        if not flags[i]:
            continue
        print(f"Generate MIP {_LABELS[i]}... ", end="", flush=True)
        mip = mips[i]
        mips[i] = None
        if mip is None:
            mip = _axis_mip_blockwise(zarr_fullname, i, box)
        print("Done!", flush=True)
        mip_fullname = f"{mip_fullname[:-10]}_MIP_{_LABELS[i]}.tif"
        _publish_tiff(_matlab_image(_cast(mip, dtype)), mip_fullname, uuid)


# ---------------------------------------------------------------------------------------------
# batch sizes and the MIPs of blocks (large-file processing)
# ---------------------------------------------------------------------------------------------


def _yxz(v: Any) -> list[float]:
    from petakit5d.core.axes import to_yxz

    return [float(x) for x in np.asarray(to_yxz(tuple(_floats(v))), dtype=np.float64).ravel()]


def _floats(v: Any) -> list[float]:
    return [float(x) for x in np.asarray(v, dtype=np.float64).ravel()]


def _zyx(v: Sequence[float]) -> tuple[int | float, ...]:
    """MATLAB's [y x z] doubles as (z, y, x), ints where integral"""
    return (_num(v[2]), _num(v[0]), _num(v[1]))


def _slab_gb(in_size: list[float], sub: list[float], byte: int) -> list[float]:
    # prod([inSize(setdiff(1 : 3, i)), bSubSz(i)]) * byteNum / 1024^3, MATLAB's order
    out = []
    for i in range(3):
        a, b = (in_size[j] for j in range(3) if j != i)
        out.append(a * b * sub[i] * byte / 1024**3)
    return out


def mip_zarr_batch_size(
    in_size: Sequence[int],
    block_size: Sequence[int],
    dtype: str,
    *,
    batch_size: Sequence[int] = (2048, 2048, 2048),
    mip_slab: bool = False,
) -> tuple[int | float, ...]:
    """``MIP_zarr_batch_size``: the batch size ``(z, y, x)`` of the large-file MIPs of a region
    of ``in_size`` with blocks of ``block_size`` (zarr chunks or pooling blocks): the requested
    batch size rounded up to whole blocks and capped at ``in_size``; for the final MIPs (not
    ``mip_slab``) increased by a block along the axis of the largest MIP slab while a slab is
    over 100 GiB and the batch is under 100 GiB (MATLAB's first largest of y, x, z)."""
    from petakit5d.core.arrays import data_type_to_byte_number

    byte = data_type_to_byte_number(dtype)
    ins, blk, bs = _yxz(in_size), _yxz(block_size), _yxz(batch_size)
    bs = [min(n, max(b, math.ceil(b / k) * k)) for n, b, k in zip(ins, bs, blk, strict=True)]
    if mip_slab:
        return _zyx(bs)
    gb = _slab_gb(ins, [math.ceil(n / b) for n, b in zip(ins, bs, strict=True)], byte)
    while any(g > 100 for g in gb) and bs[0] * bs[1] * bs[2] * byte / 1024**3 < 100:
        ind = gb.index(max(gb))
        bs = [min(ins[j], bs[j] + (blk[j] if j == ind else 0)) for j in range(3)]
        gb = _slab_gb(ins, [math.ceil(n / b) for n, b in zip(ins, bs, strict=True)], byte)
    return _zyx(bs)


def merge_memory(
    im_size: Sequence[int], region_size: Sequence[int], block_size: Sequence[int], dtype: str
) -> float:
    """``MIP_zarr_merge_memory``: the memory (GB) ``mip_zarr`` holds to merge the MIP slabs of
    a zarr of ``im_size`` from the MIPs of its blocks (a grid of ``region_size``; chunks
    ``block_size``): the slab, its coverage count and the blocks' MIPs, for the axis that needs
    the most."""
    from petakit5d.core.arrays import data_type_to_byte_number

    ims = _yxz(im_size)
    region = [min(a, b) for a, b in zip(ims, _yxz(region_size), strict=True)]
    byte = data_type_to_byte_number(dtype)
    mbs = _yxz(mip_zarr_batch_size(im_size, block_size, dtype))
    mem = 0.0
    total = ims[0] * ims[1] * ims[2]
    for i in range(3):
        n_mip = math.ceil(ims[i] / mbs[i])
        n_region = math.ceil(ims[i] / region[i])
        plane = total / ims[i]
        slab_voxels = plane * n_mip
        mip_voxels = plane * (n_region + n_mip - 1)
        mem = max(mem, (slab_voxels * (byte + 1) + mip_voxels * byte) / 2**30)
    return mem


def _range_max_parts(a: np.ndarray, p: int, starts: list[int]) -> np.ndarray:
    """max along axis p over the ranges that begin at starts (MATLAB's max for other classes)"""
    from petakit5d.kernels.mip import range_max_3d

    if a.dtype in _RANGE_DTYPES:
        return range_max_3d(a, p - 3, starts)
    if len(starts) == 1:
        return np.max(a, axis=p, keepdims=True)
    ends = [*starts[1:], a.shape[p]]
    return np.concatenate(
        [
            np.max(np.take(a, range(s, e), axis=p), axis=p, keepdims=True)
            for s, e in zip(starts, ends, strict=True)
        ],
        axis=p,
    )


def block_mips(
    block: np.ndarray, box: Sequence[int], mip_batch_size: Sequence[int], dtype: str
) -> list[np.ndarray]:
    """``MIP_zarr_block_mips``: the MIPs ``[z, y, x]`` of a block written to the region ``box``
    of a zarr, for MIP slabs with batches of ``mip_batch_size`` (``(z, y, x)``): along each
    axis, the max over each part of the block within one batch (the block's shape with one
    element per batch it overlaps along the axis), cast to the zarr's class ``dtype`` as the
    writer casts the block (rounded and saturated, which keeps the order of values)."""
    a = _vol3(np.ascontiguousarray(block))
    b = _ints(box)
    bs = _ints(mip_batch_size)
    out: list[Any] = [None, None, None]
    for i in range(3):
        p = _ZYX[i]
        b0, b1, k = b[p], b[p + 3], bs[p]
        starts = [b0, *range((b0 // k) * k + k, b1, k)]
        mip = _range_max_parts(a, p, [s - b0 for s in starts])
        out[p] = _cast(mip, dtype)
    return out


def _matlab_cell(arrays: Sequence[np.ndarray]) -> np.ndarray:
    """a MATLAB n x 1 cell of the arrays in MATLAB's (y, x, z) form"""
    from petakit5d.core.arrays import to_matlab

    cell = np.empty((len(arrays), 1), dtype=object)
    for k, m in enumerate(arrays):
        cell[k, 0] = to_matlab(m)
    return cell


def save_block_mips(
    flag_fullname: str,
    block_boxes: Sequence[Sequence[int]],
    mips: Sequence[Sequence[np.ndarray]],
    mip_batch_size: Sequence[int],
    t: float | None = None,
) -> None:
    """``MIP_zarr_save_block_mips``: the flag of a block task with the MIPs of its blocks
    (``mips``: per block the ``[z, y, x]`` of ``block_mips``), as MATLAB saves it: ``t``,
    ``blockBboxes`` (MATLAB's boxes), ``mipBatchSize`` (``[y x z]``), ``yMIPs``, ``xMIPs`` and
    ``zMIPs`` (cells of MATLAB's arrays), uncompressed (MAT v6 up to 2 GB, else v7.3), under a
    temporary name and moved into place (the flag marks the task as done)."""
    from petakit5d.core.bbox import boxes_to_matlab
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import write_mat

    rows = np.asarray(block_boxes, dtype=np.int64).reshape(-1, 6)
    boxes = boxes_to_matlab(rows).astype(np.float64)
    mbs = np.asarray(_yxz(mip_batch_size), dtype=np.float64).reshape(1, 3)
    per_axis = [[m[_ZYX[i]] for m in mips] for i in range(3)]
    nbytes = sum(int(np.asarray(m).nbytes) + _CELL_HEADER_BYTES for m3 in mips for m in m3)
    data = {
        # (MATLAB's save writes them in the order of its workspace)
        "blockBboxes": boxes,
        "mipBatchSize": mbs,
        "t": None if t is None else float(t),
        "yMIPs": _matlab_cell(per_axis[0]),
        "xMIPs": _matlab_cell(per_axis[1]),
        "zMIPs": _matlab_cell(per_axis[2]),
    }
    tmp = f"{flag_fullname[:-4]}_{get_uuid()}.mat"
    write_mat(tmp, data, version="6" if nbytes < 2**31 else "7.3", compression=False)
    os.replace(tmp, flag_fullname)


def _cell_items(v: Any) -> list[Any]:
    """the elements of a cell read from a MAT-file, in MATLAB's linear order"""
    if isinstance(v, list):
        return v
    if isinstance(v, np.ndarray) and v.dtype == object:
        return list(v.ravel(order="F"))
    return [v]


def _vol_from_matlab(m: Any) -> np.ndarray | None:
    """a MATLAB array of up to 3 dimensions as (z, y, x)"""
    from petakit5d.core.arrays import from_matlab

    a = np.asarray(m)
    if a.dtype == object or a.ndim > 3:
        return None
    return from_matlab(a.reshape(a.shape + (1,) * (3 - a.ndim)))


def _load(path: str, names: Sequence[str]) -> dict[str, Any] | None:
    from petakit5d.io import read_mat

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return read_mat(path, list(names))
        except Exception:
            return None


def merge_block_mips(
    flag_fullpaths: Sequence[str],
    axis: int,
    im_size: Sequence[int],
    batch_size: Sequence[int],
    dtype: str,
    max_mem: float = math.inf,
) -> np.ndarray | None:
    """``MIP_zarr_merge_block_mips``: the MIP slab along ``axis`` (a NumPy axis of ``(z, y,
    x)``) of a zarr of ``im_size`` with MIP batches of ``batch_size``, merged from the MIPs the
    block tasks that wrote it saved in their flags (``save_block_mips``): element k along the
    axis is the max over batch k. Parts in the same batch are merged in their order along the
    axis (the first of equal values, the first NaN). None (the caller reads the zarr) if a flag
    has no MIPs or MIPs for other batches, if the blocks don't cover the zarr exactly once, or
    if the merge would hold more than ``max_mem`` GB (NaN: no limit)."""
    from petakit5d.core.arrays import data_type_to_byte_number, matlab_class, numpy_dtype

    p = _np_axis(axis)
    ims = _ints(im_size)
    bs = _ints(batch_size)
    want = _yxz(batch_size)
    slab_size = list(ims)
    slab_size[p] = -(-ims[p] // bs[p])
    var = f"{_LABELS[_ZYX.index(p)]}MIPs"

    box_parts = []
    for f in flag_fullpaths:
        a = _load(f, ("blockBboxes", "mipBatchSize"))
        if a is None or "blockBboxes" not in a or "mipBatchSize" not in a:
            return None
        got = np.asarray(a["mipBatchSize"], dtype=np.float64).ravel(order="F")
        if got.shape != (3,) or [float(x) for x in got] != want:
            return None
        rows = np.asarray(a["blockBboxes"], dtype=np.float64)
        box_parts.append(rows.reshape(-1, 6) if rows.size else np.zeros((0, 6)))
    rows = np.concatenate(box_parts) if box_parts else np.zeros((0, 6))
    if rows.shape[0] == 0 or not np.all(rows == np.floor(rows)):
        return None
    # MATLAB's rows [ymin xmin zmin ymax xmax zmax] as (z0, y0, x0) starts and (z1, y1, x1) stops
    starts = rows[:, [2, 0, 1]].astype(np.int64) - 1
    stops = rows[:, [5, 3, 4]].astype(np.int64)
    ext = stops - starts
    if (
        np.any(starts < 0)
        or np.any(stops > np.asarray(ims))
        or int(np.prod(ext, axis=1).sum()) != int(np.prod(ims))
    ):
        return None

    k1 = starts[:, p] // bs[p]
    k2 = (stops[:, p] - 1) // bs[p]
    mip_sizes = ext.copy()
    mip_sizes[:, p] = k2 - k1 + 1
    byte = data_type_to_byte_number(dtype)
    need = (
        float(np.prod(slab_size)) * (byte + 1) + float(np.prod(mip_sizes, axis=1).sum()) * byte
    ) / 2**30
    if need > max_mem:
        return None

    mips: list[Any] = []
    for f, part in zip(flag_fullpaths, box_parts, strict=True):
        a = _load(f, (var,))
        if a is None or var not in a:
            return None
        items = _cell_items(a[var])
        if len(items) != part.shape[0]:
            return None
        mips.extend(items)

    parts: list[np.ndarray] = []
    part_boxes = []
    part_starts = []
    for b in range(rows.shape[0]):
        m = _vol_from_matlab(mips[b])
        if m is None or matlab_class(m.dtype) != dtype or m.shape != tuple(mip_sizes[b]):
            return None
        mips[b] = None
        for k in range(int(k1[b]), int(k2[b]) + 1):
            sl = [slice(None)] * 3
            sl[p] = slice(k - int(k1[b]), k - int(k1[b]) + 1)
            parts.append(m[tuple(sl)])
            lo, hi = starts[b].copy(), stops[b].copy()
            lo[p], hi[p] = k, k + 1
            part_boxes.append((*lo, *hi))
            part_starts.append(int(starts[b, p]))

    # parts of blocks in the same batch have the same region; the regions cover the slab once
    groups, ic = np.unique(np.asarray(part_boxes, dtype=np.int64), axis=0, return_inverse=True)
    ic = ic.ravel()
    covered = np.zeros(slab_size, dtype=np.uint8)
    for g in groups:
        v = covered[g[0] : g[3], g[1] : g[4], g[2] : g[5]]
        np.add(v, 1, out=v, where=v < 255)
    if np.any(covered != 1):
        return None
    del covered

    # merge the parts of each region in their order along the axis
    slab = np.zeros(slab_size, dtype=numpy_dtype(dtype))
    order = np.lexsort((np.asarray(part_starts), ic))
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and ic[order[j + 1]] == ic[order[i]]:
            j += 1
        group = [parts[q] for q in order[i : j + 1]]
        mip = (
            group[0] if len(group) == 1 else _range_max_parts(np.concatenate(group, axis=p), p, [0])
        )
        g = groups[ic[order[i]]]
        slab[g[0] : g[3], g[1] : g[4], g[2] : g[5]] = mip
        for q in order[i : j + 1]:
            parts[q] = None  # type: ignore[call-overload]
        i = j + 1
    return slab


# ---------------------------------------------------------------------------------------------
# large zarr files: batches, MIP slabs, final MIPs
# ---------------------------------------------------------------------------------------------


def _paths_from_matlab(v: Any) -> tuple[str, ...]:
    """MATLAB's {y, x, z} paths as (z, y, x)"""
    p = list(v)
    return (p[2], p[0], p[1]) if len(p) == 3 else tuple(p)


def _paths_to_matlab(v: Any) -> tuple[str, ...]:
    p = list(v)
    return (p[1], p[2], p[0]) if len(p) == 3 else tuple(p)


def _slab_shape(path: str) -> list[float]:
    from petakit5d.io import get_zarr_info

    try:
        return [int(v) for v in get_zarr_info(path).size]
    except FileNotFoundError:
        return [math.inf, math.inf, math.inf]  # (write_zarr says the slab is missing)


def _fit(out: np.ndarray, shape: Sequence[int]) -> np.ndarray:
    """MIP_block's padarray(..., 'post', 'replicate') and crop3d to the box's size"""
    d = [s - n for s, n in zip(shape, out.shape, strict=True)]
    if any(x > 0 for x in d):
        out = np.pad(out, [(0, max(x, 0)) for x in d], mode="edge")
    if any(x < 0 for x in d):
        out = np.ascontiguousarray(out[: shape[0], : shape[1], : shape[2]])
    return out


@register(
    "MIP_block",
    Param("batchInds", "@isnumeric", "Batch numbers (0-based).", convert="index1"),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr."),
    Param(
        "MIPFullpaths",
        "@(x) iscell(x)",
        "MIP slab zarrs along z, y and x.",
        convert=(_paths_from_matlab, _paths_to_matlab),
    ),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param("BatchBBoxes", "@isnumeric", "Box of each batch.", convert="bboxes", py="batch_bboxes"),
    Param("startCoord", "@isnumeric", "Start (z, y, x) of the region.", convert="point1_yxz"),
    Param(
        "poolSize",
        "@isnumeric",
        "Pool sizes (z, y, x): along the MIP axis[, across it[, resampling]].",
        convert="size_yxz",
    ),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("uuid", "@ischar", "Id of the call."),
    Param("debug", "@islogical", "Not used."),
    task="T02",
)
def mip_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    mip_fullpaths: Sequence[str],
    flag_fullname: str,
    batch_bboxes: Sequence[Sequence[int]],
    start_coord: Sequence[int],
    pool_size: Sequence[float],
    *,
    overwrite: bool = False,
    uuid: str = "",
    debug: bool = False,
) -> list[bool] | None:
    """``MIP_block``: for each batch (``batch_bboxes``, in the zarr), its MIP (or max pooling)
    along each axis written into that axis's slab zarr (``mip_fullpaths``, ``(z, y, x)``) at the
    batch's place relative to ``start_coord``; then the empty flag file. Pool sizes: 3 values
    along the MIP axis, 6 with the pooling across it, 9 with resampling (not ported).
    Nothing is done if the flag exists (unless ``overwrite``)."""
    from petakit5d.io import ReadBuffer, read_zarr, write_zarr
    from petakit5d.kernels.mip import max_pooling_3d

    del uuid, debug
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise MipError(f"The block directory {flag_path} does not exist, skip the processing!")
    inds = _ints(batch_inds)
    if _exist(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            print(
                f"The batch files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!",
                flush=True,
            )
            return None
    if not os.path.isdir(zarr_fullpath):
        raise MipError(f"The input zarr file {zarr_fullpath} doesnot exist!")
    boxes = np.asarray(batch_bboxes, dtype=np.int64).reshape(-1, 6)
    start = _ints(start_coord)
    pool = _floats(pool_size)
    pool_orig = pool + [1.0, 1.0, 1.0] if len(pool) == 3 else list(pool)
    slab_shapes = [_slab_shape(f) for f in mip_fullpaths]

    # if the pool size is smaller than the batch size, pool; with 6 values only, pool across
    first = boxes[0, 3:] - boxes[0, :3]
    larger = any(float(first[q]) > pool[q] for q in range(3))
    max_pooling = larger or (len(pool) == 6 and any(v > 1 for v in pool[3:6]))
    if larger and len(pool) == 9 and any(v > 1 for v in pool[6:9]):
        raise NotImplementedError("MIP_block's resampling before pooling (imresize3) is not ported")
    pool_1 = [1.0, 1.0, 1.0]
    pre_mip = False
    if max_pooling and len(pool) >= 6:
        pool_1 = pool[3:6]
        pool = pool[0:3]
        if all(math.fmod(a, b) == 0 for a, b in zip(pool, pool_1, strict=True)) and any(
            v > 1 for v in pool_1
        ):
            pool = [a / b for a, b in zip(pool, pool_1, strict=True)]
            pre_mip = True

    done = [False] * len(inds)
    # (one buffer for the batches' reads; each batch is written before the next is read)
    buf = ReadBuffer()
    for i, bi in enumerate(inds):
        print(f"Process Batch {bi + 1} ... ", end="", flush=True)
        t0 = time.perf_counter()
        box = [int(v) for v in boxes[i]]
        in_batch = _vol3(read_zarr(zarr_fullpath, input_bbox=box, out=buf))
        if pre_mip:
            in_batch = max_pooling_3d(in_batch, pool_1)
        for j in range(3):
            p = _ZYX[j]
            pool_j = [1.0, 1.0, 1.0] if pre_mip else list(pool_1)
            pool_j[p] = pool[p]
            if max_pooling:
                out = max_pooling_3d(in_batch, pool_j)
            else:
                out = axis_max(in_batch, _AXES[j])
            pool_out = list(pool_orig[3:6])
            pool_out[p] = pool_orig[p]
            ob0 = [math.floor((box[q] - start[q]) / pool_out[q]) for q in range(3)]
            ob1 = [math.floor((box[q + 3] - 1 - start[q]) / pool_out[q]) + 1 for q in range(3)]
            shape = [b - a for a, b in zip(ob0, ob1, strict=True)]
            if list(out.shape) != shape:
                warnings.warn(
                    "The MIP out size does not match that in the bounding box!", stacklevel=2
                )
                out = _fit(out, shape)
            if any(b > n for b, n in zip(ob1, slab_shapes[p], strict=True)):
                # (MATLAB: a box beyond the slab fails in writezarr)
                raise MipError("Upper bound is invalid", "zarr:inputError")
            write_zarr(out, mip_fullpaths[p], bbox=(*ob0, *ob1), create=False)
        # (one batch in memory: MATLAB holds the last one while it reads the next)
        del in_batch, out
        done[i] = True
        _toc(t0)
    if all(done):
        open(flag_fullname, "w").close()
    return done


def _mround(x: float) -> int:
    from petakit5d.core.dtypes import mround

    return int(mround(x))


@register(
    "XR_MIP_zarr",
    Param("zarrFullpath", "@(x) ischar(x)", "Zarr file."),
    Param("resultDirName", "@ischar", "Result folder next to the zarr."),
    Param("axis", "@isnumeric", "Axes (z, y, x) to project along (nonzero).", convert="flags_yxz"),
    Param(
        "inputBbox",
        "@(x) isempty(x) || isvector(x)",
        "Region (z0, y0, x0, z1, y1, x1) to project (None: all).",
        convert="bbox",
    ),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x).", convert="size_yxz"),
    Param(
        "poolSize",
        "@isnumeric",
        "Pool sizes (z, y, x) for the slabs: along the axis[, across it[, resampling]].",
        convert="size_yxz",
    ),
    Param("mipSlab", "@islogical", "Only the MIP slabs (no MIP tiffs)."),
    Param("parseCluster", "@islogical", "Run the batches as Slurm jobs."),
    Param("parseParfor", "@islogical", "Run the batches on a local process pool."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("masterCompute", "@islogical", "This process computes batches too."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param(
        "blockFlagFullpaths",
        "@iscell",
        "Flags of the block tasks that wrote the zarr, with the MIPs of their blocks.",
    ),
    task="T02",
)
def mip_zarr(
    zarr_fullpath: str,
    *,
    result_dir_name: str = "MIPs",
    axis: Sequence[float] = (1, 1, 1),
    input_bbox: Sequence[int] | None = None,
    batch_size: Sequence[int] = (2048, 2048, 2048),
    pool_size: Sequence[float] | None = None,
    mip_slab: bool = False,
    parse_cluster: bool = True,
    parse_parfor: bool = False,
    job_log_dir: str = "../job_logs/",
    master_compute: bool = True,
    mcc_mode: bool = False,
    config_file: str = "",
    uuid: str = "",
    debug: bool = False,
    block_flag_fullpaths: Sequence[str] = (),
) -> None:
    """MIPs of a large zarr file (``XR_MIP_zarr``): the MIPs of its batches along all three axes
    (``mip_block`` tasks through the job framework) into the MIP slabs
    ``<result_dir_name>/<name>_{y,x,z}.zarr``, then the MIPs ``<name>_MIP_{y,x,z}.tif`` of the
    flagged axes from the slabs. With the flags of the block tasks that wrote the zarr (and no
    box or pool size) the slabs are merged from the MIPs of their blocks instead of reading the
    zarr. Existing MIPs (or slabs, with ``mip_slab``) are kept; existing slabs give the MIPs; a
    slab another writer published is kept. The slabs are published before the MIPs."""
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.memory import get_system_memory
    from petakit5d.core.paths import fileparts, mkdir_recursive, strip_trailing_sep
    from petakit5d.io import create_zarr, get_zarr_info, write_zarr
    from petakit5d.jobs import framework
    from petakit5d.jobs.publish import publish_result
    from petakit5d.jobs.tasks import Task

    flags = _flags(axis)
    if not uuid:
        uuid = get_uuid()
    zarr_fullpath = strip_trailing_sep(zarr_fullpath)
    data_path, fsname, _ext = fileparts(zarr_fullpath)
    mip_path = f"{data_path}/{result_dir_name}/"
    if not os.path.isdir(mip_path):
        os.makedirs(mip_path, exist_ok=True)
    mip_fullpaths = [f"{mip_path}/{fsname}_MIP_{a}.tif" for a in _LABELS]
    slab_paths = [f"{mip_path}/{fsname}_{a}.zarr" for a in _LABELS]

    done_flag = [
        (not flags[i])
        or ((not mip_slab and _exist(mip_fullpaths[i])) or (mip_slab and _exist(slab_paths[i])))
        for i in range(3)
    ]
    if all(done_flag):
        print("The output results exist, skip it!", flush=True)
        return

    print(f"Start Large-file MIP for {zarr_fullpath}...", flush=True)
    zarr_flag_path = f"{mip_path}/zarr_flag/{fsname}_{uuid}/"
    if not os.path.isdir(zarr_flag_path):
        mkdir_recursive(zarr_flag_path)

    info = get_zarr_info(zarr_fullpath)
    im_size = [int(v) for v in info.size]
    dtype = info.dtype
    byte = data_type_to_byte_number(dtype)
    in_size = list(im_size)
    start = [0, 0, 0]
    has_box = input_bbox is not None and len(input_bbox) > 0
    if has_box:
        bb = _ints(input_bbox)
        in_size = [bb[q + 3] - bb[q] for q in range(3)]
        start = bb[:3]

    # pool sizes along the MIP axis (pool), across it (pool_1) and resampling factors
    pool_1 = [1, 1, 1]
    ds = [1, 1, 1]
    pool: list[Any] | None = None
    if pool_size is None or len(pool_size) == 0:
        inblock = [int(v) for v in info.block_size]
    else:
        pv = [_num(v) for v in _floats(pool_size)]
        if len(pv) == 9:
            ds = pv[6:9]
        if len(pv) >= 6:
            pool_1 = pv[3:6]
        pool = pv[0:3]
        inblock = [math.lcm(int(a), int(b)) for a, b in zip(pool, pool_1, strict=True)]
        inblock_1 = [math.lcm(a, int(b)) for a, b in zip(inblock, info.block_size, strict=True)]
        # limit the block size to 10 GB
        if float(np.prod(inblock_1, dtype=np.float64)) * byte / 2**30 < 10:
            inblock = inblock_1

    batch = [
        _num(v)
        for v in mip_zarr_batch_size(
            in_size, inblock, dtype, batch_size=batch_size, mip_slab=mip_slab
        )
    ]
    b_sub = [math.ceil(n / b) for n, b in zip(in_size, batch, strict=True)]
    num_batch = b_sub[0] * b_sub[1] * b_sub[2]
    # batches in MATLAB's order (y fastest, then x, then z); the last ends at the end of the
    # input box, within the image (MATLAB since 11c7d0b; before, at the image's end)
    end = [min(start[q] + in_size[q], im_size[q]) for q in range(3)]
    batch_boxes = []
    for iz in range(b_sub[0]):
        for ix in range(b_sub[2]):
            for iy in range(b_sub[1]):
                s = [iz * batch[0] + start[0], iy * batch[1] + start[1], ix * batch[2] + start[2]]
                batch_boxes.append(
                    [int(v) for v in s] + [int(min(s[q] + batch[q], end[q])) for q in range(3)]
                )

    tmp_paths = [f"{mip_path}/{fsname}_{a}_{uuid}.zarr" for a in _LABELS]
    zarr_done = [(not flags[i]) or os.path.isdir(slab_paths[i]) for i in range(3)]

    # all slabs exist: the MIPs from them
    if all(zarr_done) and not mip_slab:
        for i in range(3):
            if flags[i]:
                save_mip_zarr(slab_paths[i], mip_fullpaths[i], dtype=dtype, axis=_one_axis(i))
        return

    # merge the slabs from the MIPs the block tasks saved in their flags, if they all have them
    fuse = False
    if pool is None:
        pool = list(batch)
        fuse = len(block_flag_fullpaths) > 0 and not has_box

    # size, chunks and dimension separator of the slabs (MATLAB's order y, x, z)
    out_sizes, block_sizes, seps = [], [], []
    for i in range(3):
        p = _ZYX[i]
        out = [math.ceil(in_size[q] / (pool[q] if q == p else pool_1[q])) for q in range(3)]
        blk = [math.ceil(batch[q] / (pool[q] if q == p else pool_1[q])) for q in range(3)]
        n_chunks = 1
        for o, b in zip(out, blk, strict=True):
            n_chunks *= math.ceil(o / b)
        out_sizes.append(out)
        block_sizes.append(blk)
        seps.append("/" if n_chunks > 10000 else ".")

    if fuse:
        # up to half of the memory (no limit if it is unknown)
        max_mem = get_system_memory() / 2
        for i in range(3):
            slab = merge_block_mips(block_flag_fullpaths, _ZYX[i], in_size, batch, dtype, max_mem)
            if slab is None:
                fuse = False
                break
            create_zarr(
                tmp_paths[i],
                data_size=out_sizes[i],
                block_size=block_sizes[i],
                dtype=dtype,
                dim_separator=seps[i],
            )
            write_zarr(slab, tmp_paths[i], bbox=(0, 0, 0, *out_sizes[i]), create=False)
            del slab
        if fuse:
            print("Merged the MIP slabs from the MIPs of the blocks.", flush=True)
        else:
            print(
                "The MIPs of the blocks can't be merged, compute the MIPs from the zarr file.",
                flush=True,
            )

    if not fuse:
        for i in range(3):
            if not os.path.isdir(tmp_paths[i]):
                create_zarr(
                    tmp_paths[i],
                    data_size=out_sizes[i],
                    block_size=block_sizes[i],
                    dtype=dtype,
                    dim_separator=seps[i],
                )
        # the batches of a task
        task_size = max(5, min(10, _mround(num_batch / 5000)))
        num_tasks = math.ceil(num_batch / task_size)
        tasks, outputs = [], []
        slabs_zyx = [tmp_paths[2], tmp_paths[0], tmp_paths[1]]
        pools = [*pool, *pool_1, *ds]
        for t in range(num_tasks):
            inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
            flag = f"{zarr_flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
            outputs.append(flag)
            tasks.append(
                Task.of(
                    mip_block,
                    inds,
                    zarr_fullpath,
                    slabs_zyx,
                    flag,
                    [batch_boxes[k] for k in inds],
                    start,
                    pools,
                    uuid=uuid,
                    debug=bool(debug),
                )
            )
        inputs = [zarr_fullpath] * num_tasks
        done: list[bool] = []
        if parse_cluster or not parse_parfor:
            mem = float(np.prod(batch, dtype=np.float64)) * byte / 2**30 * BATCH_MEMORY_FACTOR
            min_task_job_num = 1
            if not mcc_mode:
                min_task_job_num = max(min(num_tasks, 5), _mround(num_tasks / 50))
            for i in range(3):
                if done and all(done):
                    break
                done = framework.run_tasks(
                    inputs,
                    outputs,
                    tasks,
                    mem_allocate=mem * 2**i,
                    max_job_num=math.inf,
                    task_batch_num=1,
                    min_task_job_num=min_task_job_num,
                    master_compute=master_compute,
                    parse_cluster=parse_cluster,
                    job_log_dir=job_log_dir,
                    mcc_mode=mcc_mode,
                    config_file=config_file,
                )
        else:
            done = framework.run_tasks(
                inputs,
                outputs,
                tasks,
                cluster_type="parfor",
                task_batch_num=1,
                gpu_job=False,
                local_worker_num=12,
                uuid=uuid,
            )
        if not all(done):
            raise MipError("Some blocks are not finished!")

    # a slab another writer put in place (of the same result) is kept
    for i in range(3):
        if os.path.isdir(tmp_paths[i]):
            publish_result(tmp_paths[i], slab_paths[i])
    if os.path.isdir(zarr_flag_path):
        shutil.rmtree(zarr_flag_path)
    if mip_slab:
        return
    for i in range(3):
        if flags[i]:
            save_mip_zarr(slab_paths[i], mip_fullpaths[i], dtype=dtype, axis=_one_axis(i))


def _one_axis(i: int) -> tuple[int, int, int]:
    """the (z, y, x) flags of MATLAB's axis i alone (``(1 : 3) == i``)"""
    f = [0, 0, 0]
    f[_ZYX[i]] = 1
    return (f[0], f[1], f[2])


# ---------------------------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------------------------


def _save_parameters(result_path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) and parameters.json of the call (MATLAB's ``pr``)"""
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(mip_dataset, values).items()}
    write_mat(result_path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, result_path + "/parameters.json")


@register(
    "XR_MIP_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("axis", "@isnumeric", "Axes (z, y, x) to project along (nonzero).", convert="flags_yxz"),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to project."),
    Param(
        "inputBbox",
        "@(x) isempty(x) || isvector(x)",
        "Region (z0, y0, x0, z1, y1, x1) to project (None: all).",
        convert="bbox",
    ),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("largeFile", "@islogical", "Large zarr files: MIPs in batches (mip_zarr)."),
    Param("batchSize", "@isvector", "Batch size (z, y, x) of large files.", convert="size_yxz"),
    Param("save16bit", "@islogical", "uint16 MIPs (else single).", convert="same"),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("parseParfor", "@islogical", "Run the files on a local process pool first."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("cpusPerTask", "@isscalar", "CPUs of a task's job."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Not used."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T02",
)
def mip_dataset(
    data_paths: str | Sequence[str],
    *,
    result_dir_name: str = "MIPs",
    axis: Sequence[float] = (1, 0, 0),
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamA_ch1", "CamB_ch0", "CamB_ch1"),
    input_bbox: Sequence[int] | None = None,
    zarr_file: bool = False,
    large_file: bool = False,
    batch_size: Sequence[int] = (2048, 2048, 2048),
    save_16bit: bool = True,
    parse_cluster: bool = True,
    parse_parfor: bool = False,
    master_compute: bool = True,
    cpus_per_task: float = 3,
    job_log_dir: str = "../job_logs/",
    uuid: str = "",
    debug: bool = False,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """The MIPs of every file of the datasets that matches a channel pattern
    (``XR_MIP_wrapper``), as ``<data path>/<result_dir_name>/<name>_MIP_{y,x,z}.tif``: tiff
    files (``save_mip_tiff``), zarr files (``save_mip_zarr``) or large zarr files
    (``mip_zarr``, each file's task runs its batches as jobs), one task per file through the job
    framework, uint16 (or single without ``save_16bit``). parameters.mat and .json of the call
    go into each result folder."""
    values = {k: v for k, v in locals().items() if k in _WRAPPER_ARGS}
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import simplify_path
    from petakit5d.io import get_image_size, get_zarr_info
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.config import cluster_config_value
    from petakit5d.jobs.coordinator import coordinator_nice, coordinator_resources
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_image_filenames

    _flags(axis)
    if not uuid:
        uuid = get_uuid()
    paths = [data_paths] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    dtype = "uint16" if save_16bit else "single"

    result_paths = []
    for d, data_path in enumerate(paths):
        result_path = f"{data_path}/{result_dir_name}/"
        os.makedirs(result_path, exist_ok=True)
        paths[d] = simplify_path(data_path)
        result_paths.append(simplify_path(result_path))
        group_writable(result_path)
        _save_parameters(result_path, values)

    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(paths[0], job_log_dir).parse_cluster)

    _, fsns, fd_inds, filepaths, _ = parse_image_filenames(paths, zarr_file, channel_patterns)
    n_files = len(fsns)
    sz = get_image_size(filepaths[0])

    # the last flagged axis names the MIP that marks a file as done
    flags = _flags(axis)
    axis_str = _LABELS[max(i for i in range(3) if flags[i])] if any(flags) else None
    if axis_str is None:
        raise MipError(
            "Index exceeds the number of array elements. Index must not exceed 3.",
            "MATLAB:badsubscript",
        )
    mip_fullpaths, tasks = [], []
    for f in range(n_files):
        frame = filepaths[f]
        mip = f"{result_paths[int(fd_inds[f])]}/{fsns[f]}_MIP_{axis_str}.tif"
        mip_fullpaths.append(mip)
        if zarr_file and large_file:
            kw: dict[str, Any] = {
                "result_dir_name": result_dir_name,
                "axis": list(axis),
                "input_bbox": None if input_bbox is None else list(input_bbox),
                "batch_size": list(batch_size),
                "parse_cluster": bool(parse_cluster),
                "parse_parfor": bool(parse_parfor),
                "job_log_dir": job_log_dir,
                "mcc_mode": bool(mcc_mode),
                "config_file": config_file,
                "uuid": uuid,
            }
            # (only then, as the default of mip_zarr is true)
            if not master_compute:
                kw["master_compute"] = False
            tasks.append(Task.of(mip_zarr, frame, **kw))
        else:
            func = save_mip_zarr if zarr_file else save_mip_tiff
            tasks.append(
                Task.of(
                    func,
                    frame,
                    mip,
                    dtype=dtype,
                    axis=list(axis),
                    input_bbox=None if input_bbox is None else list(input_bbox),
                )
            )

    if parse_parfor and not large_file:
        framework.run_tasks(
            filepaths,
            mip_fullpaths,
            tasks,
            cluster_type="parfor",
            gpu_job=False,
            local_worker_num=12,
        )
    task_batch_num = max(1, _mround(n_files / 1000))
    coord_nice = 0
    if zarr_file and large_file:
        # each file's job (mip_zarr) runs its MIP batches as jobs: it holds the process, and if
        # it computes batches itself (master_compute, the config's if it sets it), cpus_per_task
        # and a batch's estimate
        info = get_zarr_info(filepaths[0])
        in_size = [int(v) for v in info.size]
        if input_bbox is not None and len(input_bbox):
            bb = _ints(input_bbox)
            in_size = [bb[q + 3] - bb[q] for q in range(3)]
        mbs = mip_zarr_batch_size(in_size, info.block_size, info.dtype, batch_size=batch_size)
        batch_mem = (
            float(np.prod(mbs, dtype=np.float64))
            * data_type_to_byte_number(info.dtype)
            / 2**30
            * BATCH_MEMORY_FACTOR
        )
        computes = cluster_config_value(config_file, "masterCompute", master_compute)
        cpus_per_task, mem_allocate = coordinator_resources(0, computes, cpus_per_task, batch_mem)
        # (and they wait for their batch jobs: these go first)
        if parse_cluster:
            coord_nice = coordinator_nice()
    else:
        mem_allocate = float(np.prod(sz, dtype=np.float64)) * 4 / 1024**3 * 4
    framework.run_tasks(
        filepaths,
        mip_fullpaths,
        tasks,
        parse_cluster=parse_cluster,
        job_log_dir=job_log_dir,
        master_compute=master_compute,
        task_batch_num=task_batch_num,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        mcc_mode=mcc_mode,
        config_file=config_file,
        nice=coord_nice,
    )


_WRAPPER_ARGS = (
    "data_paths",
    "result_dir_name",
    "axis",
    "channel_patterns",
    "input_bbox",
    "zarr_file",
    "large_file",
    "batch_size",
    "save_16bit",
    "parse_cluster",
    "parse_parfor",
    "master_compute",
    "cpus_per_task",
    "job_log_dir",
    "uuid",
    "debug",
    "mcc_mode",
    "config_file",
)
