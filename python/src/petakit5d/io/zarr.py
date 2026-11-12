"""Zarr reading, writing and creation through cpp-zarr (ports readzarr, writezarr,
createzarr and writeZarrBlock).

PetaKit5D's arrays: Zarr v2, order F, shape ``[ny, nx, nz]`` (2D data as ``[ny, nx, 1]``),
blosc zstd level 1, chunk files ``y.x.z``. Data is ``(y, x, z)`` Fortran order and goes to
and from cpp-zarr without a copy.

As in MATLAB, a new array is created first (``create_zarr``: the ``.zarray`` text the mex
writes, written to a temporary name and renamed) and the data is then written into it as a
region. Until cpp-zarr closes the gaps of decision 0004:

- region writes into an array of another type convert as MATLAB does (rounding half away
  from zero, saturating) before cpp-zarr sees the data;
- all-zero chunks are not written (cpp-zarr 1.7.3 has no ``sparse`` option; it is passed
  as soon as the installed cpp-zarr takes it): the same data, fewer chunk files than
  MATLAB's ``writezarr`` writes;
- ``uuid=False`` is passed only when cpp-zarr takes it (chunks are always written to a
  temporary name and renamed otherwise, the default).
"""

from __future__ import annotations

import inspect
import json
import os
import shutil
import uuid as _uuid
import warnings
from collections.abc import Sequence
from functools import cache

import numpy as np

from petakit5d.io.imagesize import (
    matlab_class,
    matlab_shape,
    matlab_squeeze,
    numpy_dtype,
    path_str,
    threads,
)

__all__ = ["create_zarr", "read_zarr", "write_zarr", "write_zarr_block"]

# createzarr's types: MATLAB class -> .zarray dtype
_ZARR_DTYPES = {"double": "<f8", "single": "<f4", "uint16": "<u2", "uint8": "|u1"}
_BLOSC = ("lz4", "blosclz", "lz4hc", "zlib", "zstd")


@cache
def _write_options() -> frozenset:
    import cppzarr

    try:
        return frozenset(inspect.signature(cppzarr.write_zarr).parameters)
    except (TypeError, ValueError):
        return frozenset()


def _is_array(path: str) -> bool:
    return os.path.isfile(os.path.join(path, ".zarray")) or os.path.isfile(
        os.path.join(path, "zarr.json")
    )


def _ints(v, what: str) -> list[int]:
    out = []
    for x in np.asarray(v, dtype=np.float64).ravel():
        if not np.isfinite(x) or x != int(x):
            raise ValueError(f"{what} must hold whole numbers, not {v!r}")
        out.append(int(x))
    return out


def _bbox_coords(bbox, shape: Sequence[int]) -> tuple[list[int], list[int]]:
    """MATLAB's [starts, ends] (1-based, inclusive) as cpp-zarr's start and end (0-based,
    end exclusive), fitted to the array's axes (missing trailing axes are singletons)."""
    b = _ints(bbox, "bbox")
    if len(b) % 2:
        raise ValueError(f"bbox must have an even number of values, not {len(b)}")
    n = len(b) // 2
    start, end = [s - 1 for s in b[:n]], b[n:]
    nd = len(shape)
    if n < nd:
        start += [0] * (nd - n)
        end += [1] * (nd - n)
    elif n > nd:
        if any(s != 0 for s in start[nd:]) or any(e != 1 for e in end[nd:]):
            raise ValueError(f"bbox {b} has more axes than the array {list(shape)}")
        start, end = start[:nd], end[:nd]
    if any(s < 0 or e <= s or e > d for s, e, d in zip(start, end, shape, strict=False)):
        raise ValueError(f"bbox {b} is outside the array of size {list(shape)}")
    return start, end


def _meta(path: str) -> dict:
    with open(os.path.join(path, ".zarray"), "rb") as f:
        return json.loads(f.read().decode())


# --- stand-in for petakit5d.core.dtypes (F02): MATLAB's class conversions ------------------


def _round_matlab(v: np.ndarray) -> np.ndarray:
    # round half away from zero, in the input's precision (v - trunc(v) is exact)
    t = np.trunc(v)
    with np.errstate(invalid="ignore"):
        return t + np.copysign((np.abs(v - t) >= 0.5).astype(v.dtype), v)


def cast_like_matlab(a: np.ndarray, dtype) -> np.ndarray:
    """``cast(a, dtype)`` as MATLAB converts: to integers rounding half away from zero,
    saturating, NaN to 0; between floats by IEEE rounding."""
    dt = np.dtype(dtype)
    a = np.asarray(a)
    if a.dtype == dt:
        return a
    if dt.kind == "f":
        return a.astype(dt)
    if dt.kind not in "iu":
        raise TypeError(f"no MATLAB conversion to {dt}")
    info = np.iinfo(dt)
    if a.dtype.kind == "b":
        return a.astype(dt)
    if a.dtype.kind in "iu":
        lo, hi = max(info.min, np.iinfo(a.dtype).min), min(info.max, np.iinfo(a.dtype).max)
        return np.clip(a, a.dtype.type(lo), a.dtype.type(hi)).astype(dt)
    if a.dtype.kind != "f":
        raise TypeError(f"no MATLAB conversion from {a.dtype}")
    r = _round_matlab(a)
    hi_f, lo_f = float(info.max), float(info.min)
    over = r >= hi_f
    under = r <= lo_f
    nan = np.isnan(r)
    # (NaN to 0; values at the limits as floats may lie past them: set after the cast)
    safe = np.where(over | under | nan, 0, r)
    out = safe.astype(dt)
    out[over] = info.max
    out[under] = info.min
    return out


# --- createzarr ----------------------------------------------------------------------------


def _class_name(dtype) -> str:
    if isinstance(dtype, str) and dtype in _ZARR_DTYPES:
        return dtype
    try:
        return matlab_class(numpy_dtype(dtype))
    except (TypeError, ValueError):
        return str(dtype)


def zarray_text(
    data_size: Sequence[int],
    block_size: Sequence[int],
    dtype: str = "uint16",
    order: str = "F",
    compressor: str = "zstd",
    clevel: int = 1,
    dim_separator: str = ".",
    shard_size: Sequence[int] | None = None,
    zarr_sub_size: Sequence[int] | None = None,
) -> str:
    """The ``.zarray`` the mex's ``createZarrFile`` writes, byte for byte (nlohmann's
    ``setw(4)`` dump: sorted keys, 4 spaces, a newline at the end)."""
    cls = _class_name(dtype)
    if cls not in _ZARR_DTYPES:
        raise ValueError(f"Unsupported data type {cls}")
    if order not in ("F", "C"):
        raise ValueError(f"order must be 'F' or 'C', not {order!r}")
    if dim_separator not in (".", "/"):
        raise ValueError("dimension_separator must be a . or /")
    level = int(clevel)
    meta: dict = {"chunks": [int(c) for c in block_size]}
    if compressor in _BLOSC:
        meta["compressor"] = {
            "blocksize": 0,
            "clevel": level,
            "cname": compressor,
            "id": "blosc",
            "shuffle": 1,
        }
    elif compressor == "gzip":
        meta["compressor"] = {"id": "gzip", "level": level}
    elif compressor == "none":
        meta["compressor"] = None
    else:
        raise ValueError(f'Compressor: "{compressor}" is not currently supported')
    if dim_separator == "/":
        meta["dimension_separator"] = "/"
    meta["dtype"] = _ZARR_DTYPES[cls]
    meta["fill_value"] = 0
    meta["filters"] = None
    meta["order"] = order
    meta["shape"] = [int(s) for s in data_size]
    meta["zarr_format"] = 2
    if zarr_sub_size is not None and len(zarr_sub_size) and any(int(s) for s in zarr_sub_size):
        meta["subfolders"] = [int(s) for s in zarr_sub_size]
    if shard_size is not None and len(shard_size):
        if compressor in _BLOSC:
            inner = {
                "configuration": {
                    "blocksize": 0,
                    "clevel": level,
                    "cname": compressor,
                    "shuffle": "shuffle",
                    "typesize": np.dtype(_ZARR_DTYPES[cls]).itemsize,
                },
                "name": "blosc",
            }
        elif compressor == "gzip":
            inner = {"configuration": {"level": level}, "name": "gzip"}
        else:
            # (the mex's codec of uncompressed shards: an empty configuration, no name)
            inner = {"configuration": None}
        meta["codecs"] = [
            {
                "configuration": {"chunk_shape": [int(s) for s in shard_size], "codecs": [inner]},
                "name": "sharding_indexed",
            }
        ]
    return json.dumps(meta, indent=4, sort_keys=True) + "\n"


def _fit_axes(v: Sequence[int], n: int) -> list[int]:
    # the mex's mexFitAxes: extra trailing values dropped, missing ones as 1
    v = list(v)
    return v[:n] if len(v) >= n else v + [1] * (n - len(v))


def create_zarr(
    filepath,
    *,
    overwrite: bool = False,
    data_size: Sequence[int] = (500, 500, 500),
    block_size: Sequence[int] = (500, 500, 500),
    shard_size: Sequence[int] | None = None,
    dtype="uint16",
    order: str = "F",
    expand_2d_dim: bool = True,
    group_write: bool = True,
    compressor: str = "zstd",
    clevel: int = 1,
    zarr_sub_size: Sequence[int] | None = None,
    dim_separator: str = ".",
) -> None:
    """Create an empty zarr array (MATLAB ``createzarr``): an existing folder of that name
    is removed, then the ``.zarray`` is written (to a temporary name, then renamed). Types:
    ``'uint8'``, ``'uint16'``, ``'single'``, ``'double'`` (or their NumPy types)."""
    del overwrite  # (MATLAB's createzarr always replaces the folder)
    path = path_str(filepath)
    size = _ints(data_size, "dataSize")
    block = _ints(block_size, "blockSize")
    if len(size) == 2:
        if expand_2d_dim:
            size.append(1)
            block = (block + [0, 0, 0])[:2] + [1] + block[3:]
        else:
            block = block[:2]
    if len(block) != len(size):
        raise ValueError(f"blockSize {block} does not fit dataSize {size}")
    block = [min(s, b) for s, b in zip(size, block, strict=False)]
    shard = None
    if shard_size is not None and len(shard_size):
        shard = [min(b, s) for b, s in zip(block, _ints(shard_size, "shardSize"), strict=False)]
        bad = [b % s != 0 for b, s in zip(block, shard, strict=False)]
        if any(bad):
            warnings.warn(
                "The shard size needs to be divisors of the block size, setthe size "
                "of the axis that is not divisor to be the block size.",
                stacklevel=2,
            )
            shard = [b if k else s for b, s, k in zip(block, shard, bad, strict=False)]
        shard = _fit_axes(shard, len(size))
    sub = None
    if zarr_sub_size is not None and len(zarr_sub_size):
        sub = _ints(zarr_sub_size, "zarrSubSize")
        nchunks = int(np.prod([-(-s // b) for s, b in zip(size, block, strict=False)]))
        sub = None if nchunks < int(np.prod(sub)) else _fit_axes(sub, len(size))
    text = zarray_text(size, block, dtype, order, compressor, clevel, dim_separator, shard, sub)
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)
    final = os.path.join(path, ".zarray")
    tmp = final + _uuid.uuid4().hex
    try:
        with open(tmp, "w", newline="\n") as f:
            f.write(text)
        os.replace(tmp, final)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    if group_write:
        from petakit5d.io.tiff import group_writable

        group_writable(path)


# --- readzarr ------------------------------------------------------------------------------


def read_zarr(
    filepath,
    *,
    input_bbox: Sequence[int] | None = None,
    sparse_data: bool = True,
    num_threads: int | None = None,
) -> np.ndarray:
    """Read a zarr array, or the region ``input_bbox`` of it (MATLAB's ``[ymin, xmin, zmin,
    ymax, xmax, zmax]``, 1-based, inclusive), as ``(y, x, z)`` Fortran order (MATLAB
    ``readzarr``). Trailing singleton axes are dropped, as MATLAB's ``size`` does.
    ``sparse_data`` only speeds up the reads of all-zero chunks in the mex; results are the
    same, so it is not used."""
    import cppzarr

    del sparse_data
    path = path_str(filepath)
    if not _is_array(path):
        raise FileNotFoundError(f"zarr file {path} does not exist!")
    n = threads(num_threads)
    if input_bbox is None or len(input_bbox) == 0:
        a = cppzarr.read_zarr(path, order="F", num_threads=n)
    else:
        shape = _meta(path)["shape"] if os.path.isfile(os.path.join(path, ".zarray")) else None
        if shape is None:
            with open(os.path.join(path, "zarr.json"), "rb") as f:
                shape = json.loads(f.read().decode())["shape"]
        start, end = _bbox_coords(input_bbox, shape)
        a = cppzarr.read_zarr(path, start, end, order="F", num_threads=n)
    return matlab_squeeze(a)


# --- writezarr -----------------------------------------------------------------------------


def _region_write(path: str, data: np.ndarray, bbox, sparse_data: bool, uuid: bool, n: int) -> None:
    import cppzarr

    meta = _meta(path)
    start, end = _bbox_coords(bbox, meta["shape"])
    want = np.dtype(meta["dtype"]).newbyteorder("=")
    if data.dtype.newbyteorder("=") != want:
        data = cast_like_matlab(data, want)
    if [e - s for s, e in zip(start, end, strict=False)] != list(matlab_shape(data.shape)) + [1] * (
        len(start) - len(matlab_shape(data.shape))
    ):
        raise ValueError(f"The region {list(bbox)} does not match the data size {list(data.shape)}")
    kw = {}
    opts = _write_options()
    if "sparse" in opts:
        kw["sparse"] = bool(sparse_data)
    if "uuid" in opts and not uuid:
        kw["uuid"] = False
    cppzarr.write_zarr(path, data, start, end, num_threads=n, **kw)


def write_zarr(
    data,
    filepath,
    *,
    overwrite: bool = False,
    block_size: Sequence[int] = (500, 500, 500),
    shard_size: Sequence[int] | None = None,
    zarr_sub_size: Sequence[int] | None = None,
    expand_2d_dim: bool = True,
    group_write: bool = True,
    bbox: Sequence[int] | None = None,
    sparse_data: bool = False,
    dim_separator: str = ".",
    create: bool = True,
    uuid: bool = True,
    num_threads: int | None = None,
) -> None:
    """Write a ``(y, x, z)`` array as a zarr array (MATLAB ``writezarr``).

    A missing array is created (``create_zarr`` with ``block_size`` capped at the size); the
    data is then written into the array (an existing one keeps its metadata and type, the
    data converted as MATLAB converts), or into its region ``bbox`` (1-based, inclusive).
    ``overwrite`` replaces an existing array (without ``bbox``).
    """
    a = np.asarray(data)
    if a.dtype.kind not in "iuf":
        raise TypeError(f"writezarr writes numeric data, not {a.dtype}")
    path = path_str(filepath)
    sz = list(matlab_shape(a.shape))
    block = _ints(block_size, "blockSize")
    if len(sz) == 2:
        if expand_2d_dim:
            sz.append(1)
            block = (block + [0, 0, 0])[:2] + [1] + block[3:]
        else:
            block = block[:2]
    if len(block) != len(sz):
        raise ValueError(f"blockSize {block} does not fit the data size {sz}")
    block = [min(s, b) for s, b in zip(sz, block, strict=False)]
    if overwrite and os.path.isdir(path) and (bbox is None or len(bbox) == 0):
        shutil.rmtree(path)
    if bbox is None or len(bbox) == 0:
        s3 = (list(a.shape) + [1, 1, 1])[:3]
        bbox = [1, 1, 1, *s3] if expand_2d_dim else [1, 1, *s3[:2]]
    new_file = False
    if not os.path.isdir(path):
        if not create:
            raise FileNotFoundError(f"The zarr file {path} does not exist!")
        create_zarr(
            path,
            data_size=sz,
            block_size=block,
            shard_size=shard_size,
            dtype=matlab_class(a.dtype),
            zarr_sub_size=zarr_sub_size,
            dim_separator=dim_separator,
            group_write=False,
        )
        new_file = True
    _region_write(path, a, bbox, sparse_data, uuid, threads(num_threads))
    if group_write and new_file and os.name != "nt":
        from petakit5d.io.tiff import group_writable

        group_writable(path)


def write_zarr_block(
    file_path,
    block_sub: Sequence[int],
    data,
    level: int = 1,
    mode: str = "w",
    *,
    num_threads: int | None = None,
) -> None:
    """Write one chunk of a zarr array by its 1-based chunk subscript (MATLAB
    ``writeZarrBlock``, on the array's path instead of a ``blockedImage``); data past the
    array's edge is cut off. ``level`` and ``mode`` are accepted for compatibility."""
    del level, mode
    path = path_str(file_path)
    meta = _meta(path)
    shape, chunks = meta["shape"], meta["chunks"]
    sub = _fit_axes(_ints(block_sub, "blockSub"), len(shape))
    a = np.asarray(data)
    dshape = (list(a.shape) + [1] * len(shape))[: len(shape)]
    start = [(s - 1) * c for s, c in zip(sub, chunks, strict=False)]
    end = [min(st + d, n) for st, d, n in zip(start, dshape, shape, strict=False)]
    if any(e <= st for st, e in zip(start, end, strict=False)):
        raise ValueError(f"block {sub} is outside the array of size {shape}")
    keep = tuple(slice(0, e - st) for st, e in zip(start, end, strict=False))
    a = a.reshape(dshape, order="A")[keep]
    bbox = [st + 1 for st in start] + end
    _region_write(path, a, bbox, False, True, threads(num_threads))
