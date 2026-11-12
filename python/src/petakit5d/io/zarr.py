"""Zarr reading, writing and creation through cpp-zarr (ports readzarr, writezarr,
createzarr and writeZarrBlock).

The port's arrays are native: Zarr v2, shape ``[nz, ny, nx]``, ``"order": "C"``, chunks
``(cz, cy, cx)``, chunk files named by ``(z, y, x)`` chunk index, blosc zstd level 1, ``.``
separator; an image is stored as ``[1, ny, nx]``. Data is ``(z, y, x)`` C order and goes to and
from cpp-zarr without a copy. Boxes are ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive.

Arrays MATLAB wrote (``"order": "F"``, shape ``[ny, nx, nz]``) are still read: they come back
as ``(z, y, x)`` through one threaded transpose of each z plane, done in place in the buffer
cpp-zarr reads into (one volume in memory). ``axes`` overrides the detection. Region writes
into such arrays are refused (they would mix the layouts).

As in MATLAB, a new array is created first (``create_zarr``: the ``.zarray`` text, written to
a temporary name and renamed, as the mex does) and the data is then written into it as a
region. Until cpp-zarr closes the gaps of decision 0004:

- region writes into an array of another type convert as MATLAB's mex does
  (``core.dtypes.cast``: rounding half away from zero, saturating, NaN to 0) before cpp-zarr
  sees the data;
- all-zero chunks are not written (cpp-zarr 1.7.3 has no ``sparse`` option; it is passed
  as soon as the installed cpp-zarr takes it): the same data, fewer chunk files than
  MATLAB's ``writezarr`` writes;
- ``uuid=False`` is passed only when cpp-zarr takes it (chunks are always written to a
  temporary name and renamed otherwise, the default).

Sharded arrays (``shard_size``, stitching's ``shardSize``) are written as standard Zarr v3 with
the ``sharding_indexed`` codec (the user's decision), not as MATLAB's v2 ``.zarray`` with that
codec; MATLAB's v2-sharded files are still read through cpp-zarr.
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

from petakit5d.core.axes import from_yxz
from petakit5d.io.imagesize import path_str, threads, zarr_layout, zarr_meta

__all__ = [
    "create_zarr",
    "read_zarr",
    "write_zarr",
    "write_zarr_block",
    "zarr_json_text",
    "zarray_text",
]

# createzarr's types: MATLAB class -> .zarray dtype
_ZARR_DTYPES = {"double": "<f8", "single": "<f4", "uint16": "<u2", "uint8": "|u1"}
_BLOSC = ("lz4", "blosclz", "lz4hc", "zlib", "zstd")
# MATLAB's defaults dataSize and blockSize [500, 500, 500] ([y x z]) in (z, y, x)
_SIZE = tuple(from_yxz((500, 500, 500)))
_BLOCK = tuple(from_yxz((500, 500, 500)))


@cache
def _write_options() -> frozenset:
    import cppzarr

    try:
        return frozenset(inspect.signature(cppzarr.write_zarr).parameters)
    except (TypeError, ValueError):
        return frozenset()


@cache
def _inplace_kernel():
    """The in-place plane transpose of petakit5d._ext.io, or None (core.arrays.from_matlab)."""
    from petakit5d import _ext

    if _ext.available("io"):
        return _ext.load("io").transpose_pages_inplace
    warnings.warn(
        "petakit5d._ext.io is not built: MATLAB-layout zarr reads use a slower NumPy copy",
        RuntimeWarning,
        stacklevel=3,
    )
    return None


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


def _region(box, shape: Sequence[int]) -> tuple[list[int], list[int]]:
    """A box ``(z0, y0, x0, z1, y1, x1)`` (or ``(y0, x0, y1, x1)``; ends may be None) as the
    start and stop in an array of ``shape`` (``(z, y, x)`` or ``(y, x)``). A 2D box on a volume
    of one slice takes that slice."""
    from petakit5d.core.bbox import check_box

    b = list(check_box(box))
    h = len(b) // 2
    shape = [int(s) for s in shape]
    if h == 2 and len(shape) == 3:
        if shape[0] != 1:
            raise ValueError(f"a 2D box {b} on an array of {shape[0]} slices; give z0 and z1")
        b = [0, *b[:2], 1, *b[2:]]
        h = 3
    elif h == 3 and len(shape) == 2 and b[0] == 0 and b[3] in (1, None):
        # a box of one slice on an image stored [ny, nx] (createzarr without expand2dDim)
        b = [b[1], b[2], b[4], b[5]]
        h = 2
    if h != len(shape):
        raise ValueError(f"box {b} does not fit an array of shape {shape}")
    start = b[:h]
    stop = [shape[i] if b[h + i] is None else b[h + i] for i in range(h)]
    if any(s < 0 or e <= s or e > n for s, e, n in zip(start, stop, shape, strict=True)):
        raise ValueError(f"box {b} is empty or outside the array of shape {shape}")
    return start, stop


def _squeeze_z(a: np.ndarray) -> np.ndarray:
    # an image (one slice) comes back as (y, x), as MATLAB's size has it (a view)
    return a[0] if a.ndim == 3 and a.shape[0] == 1 else a


def _convert(data: np.ndarray, want: np.dtype) -> np.ndarray:
    """``data`` in the array's type, as MATLAB's mex converts (``cast``: integers rounded half
    away from zero and saturated, NaN to 0); cpp-zarr itself would use NumPy's ``astype`` (Z3)."""
    from petakit5d.core.arrays import matlab_class
    from petakit5d.core.dtypes import cast

    return np.asarray(cast(data, matlab_class(want)))


# --- createzarr ----------------------------------------------------------------------------


def _class_name(dtype) -> str:
    from petakit5d.core.arrays import matlab_class

    if isinstance(dtype, str) and dtype in _ZARR_DTYPES:
        return dtype
    try:
        return matlab_class(np.dtype(dtype))
    except TypeError:
        return str(dtype)


def zarray_text(
    data_size: Sequence[int],
    block_size: Sequence[int],
    dtype: str = "uint16",
    order: str = "C",
    compressor: str = "zstd",
    clevel: int = 1,
    dim_separator: str = ".",
    zarr_sub_size: Sequence[int] | None = None,
) -> str:
    """The ``.zarray`` the mex's ``createZarrFile`` writes for these values, byte for byte
    (nlohmann's ``setw(4)`` dump: sorted keys, 4 spaces, a newline at the end). The sizes are
    the file's (``[nz, ny, nx]`` for the port's arrays)."""
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
    return json.dumps(meta, indent=4, sort_keys=True) + "\n"


_V3_TYPES = {"double": "float64", "single": "float32", "uint16": "uint16", "uint8": "uint8"}


def zarr_json_text(
    data_size: Sequence[int],
    block_size: Sequence[int],
    shard_size: Sequence[int],
    dtype: str = "uint16",
    compressor: str = "zstd",
    clevel: int = 1,
    dim_separator: str = ".",
) -> str:
    """The ``zarr.json`` cpp-zarr writes for a Zarr v3 array with sharding: shards of
    ``block_size`` (MATLAB's blockSize) holding chunks of ``shard_size`` (MATLAB's shardSize), C
    order, the ``default`` chunk keys (``c.iz.iy.ix``, or ``c/iz/iy/ix``)."""
    cls = _class_name(dtype)
    if cls not in _V3_TYPES:
        raise ValueError(f"Unsupported data type {cls}")
    if dim_separator not in (".", "/"):
        raise ValueError("dimension_separator must be a . or /")
    level = int(clevel)
    # (one-byte values have no byte order: cpp-zarr's bytes codec then has no configuration)
    one_byte = np.dtype(_V3_TYPES[cls]).itemsize == 1
    codecs: list = [
        {"name": "bytes"} if one_byte else {"configuration": {"endian": "little"}, "name": "bytes"}
    ]
    if compressor in _BLOSC:
        codecs.append(
            {
                "configuration": {
                    "blocksize": 0,
                    "clevel": level,
                    "cname": compressor,
                    "shuffle": "shuffle",
                    "typesize": np.dtype(_V3_TYPES[cls]).itemsize,
                },
                "name": "blosc",
            }
        )
    elif compressor == "gzip":
        codecs.append({"configuration": {"level": level}, "name": "gzip"})
    elif compressor != "none":
        raise ValueError(f'Compressor: "{compressor}" is not currently supported')
    meta = {
        "attributes": {},
        "chunk_grid": {
            "configuration": {"chunk_shape": [int(b) for b in block_size]},
            "name": "regular",
        },
        "chunk_key_encoding": {"configuration": {"separator": dim_separator}, "name": "default"},
        "codecs": [
            {
                "configuration": {
                    "chunk_shape": [int(c) for c in shard_size],
                    "codecs": codecs,
                    "index_codecs": [
                        {"configuration": {"endian": "little"}, "name": "bytes"},
                        {"name": "crc32c"},
                    ],
                    "index_location": "end",
                },
                "name": "sharding_indexed",
            }
        ],
        "data_type": _V3_TYPES[cls],
        "fill_value": 0.0 if cls in ("single", "double") else 0,
        "node_type": "array",
        "shape": [int(v) for v in data_size],
        "zarr_format": 3,
    }
    return json.dumps(meta, indent=4, sort_keys=True) + "\n"


def _fit_axes(v: Sequence[int], n: int) -> list[int]:
    # the mex's mexFitAxes: extra trailing values dropped, missing ones as 1
    v = list(v)
    return v[:n] if len(v) >= n else v + [1] * (n - len(v))


def _fit2d(block: list[int], size: list[int], expand: bool) -> tuple[list[int], list[int]]:
    # an image (y, x): stored [1, ny, nx] with chunks (1, cy, cx), or as it is without expand
    yx = block[-2:] if len(block) >= 2 else block * 2
    return ([1, *size], [1, *yx]) if expand else (size, yx)


def create_zarr(
    filepath,
    *,
    overwrite: bool = False,
    data_size: Sequence[int] = _SIZE,
    block_size: Sequence[int] = _BLOCK,
    shard_size: Sequence[int] | None = None,
    dtype="uint16",
    expand_2d_dim: bool = True,
    group_write: bool = True,
    compressor: str = "zstd",
    clevel: int = 1,
    zarr_sub_size: Sequence[int] | None = None,
    dim_separator: str = ".",
) -> None:
    """Create an empty native zarr array (MATLAB ``createzarr``): an existing folder of that
    name is removed, then the ``.zarray`` is written (to a temporary name, then renamed).

    ``data_size``, ``block_size``, ``shard_size`` and ``zarr_sub_size`` are in ``(z, y, x)``;
    an image size ``(ny, nx)`` is stored as ``[1, ny, nx]`` (as ``[ny, nx]`` without
    ``expand_2d_dim``). Types: ``'uint8'``, ``'uint16'``, ``'single'``, ``'double'`` (or their
    NumPy types). MATLAB's ``order`` option is not taken: native arrays are C order. With
    ``shard_size`` the array is Zarr v3: shards of ``block_size`` holding chunks of
    ``shard_size`` (``zarr.json``).
    """
    del overwrite  # (MATLAB's createzarr always replaces the folder)
    path = path_str(filepath)
    size = _ints(data_size, "dataSize")
    block = _ints(block_size, "blockSize")
    if len(size) == 2:
        size, block = _fit2d(block, size, expand_2d_dim)
    if len(block) != len(size):
        raise ValueError(f"blockSize {block} does not fit dataSize {size}")
    block = [min(s, b) for s, b in zip(size, block, strict=True)]
    shard = None
    if shard_size is not None and len(shard_size):
        sh = _ints(shard_size, "shardSize")
        if len(size) == 3 and len(sh) == 2:
            sh = [1, *sh]
        shard = [min(b, s) for b, s in zip(block, sh, strict=False)]
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
        nchunks = int(np.prod([-(-s // b) for s, b in zip(size, block, strict=True)]))
        sub = None if nchunks < int(np.prod(sub)) else _fit_axes(sub, len(size))
    if shard is not None:
        # sharded arrays are standard Zarr v3 (the user's decision; MATLAB writes a v2 .zarray
        # with a sharding_indexed codec, which is still read)
        if sub is not None:
            raise ValueError("zarrSubSize (subfolders) and shardSize cannot be combined in Zarr v3")
        name = "zarr.json"
        text = zarr_json_text(size, block, shard, dtype, compressor, clevel, dim_separator)
    else:
        name = ".zarray"
        text = zarray_text(size, block, dtype, "C", compressor, clevel, dim_separator, sub)
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)
    final = os.path.join(path, name)
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
    axes: str = "auto",
    num_threads: int | None = None,
) -> np.ndarray:
    """Read a zarr array, or its box ``input_bbox`` (``(z0, y0, x0, z1, y1, x1)``, 0-based,
    end-exclusive), as a ``(z, y, x)`` C-ordered array (MATLAB ``readzarr``); one slice comes
    back as an image ``(y, x)``.

    ``axes``: ``'auto'`` reads an array with ``"order": "F"`` as MATLAB's ``[ny, nx, nz]``
    (one threaded transposing pass) and others as native ``[nz, ny, nx]`` (no copy); ``'zyx'``
    or ``'yxz'`` say which it is. ``sparse_data`` only speeds up the mex's reads of all-zero
    chunks; results are the same, so it is not used.
    """
    import cppzarr

    del sparse_data
    path = path_str(filepath)
    if not _is_array(path):
        raise FileNotFoundError(f"zarr file {path} does not exist!")
    meta = zarr_meta(path)
    layout = zarr_layout(meta, axes)
    n = threads(num_threads)
    fshape = [int(v) for v in meta["shape"]]
    whole = input_bbox is None or len(input_bbox) == 0
    if layout == "zyx":
        if whole:
            a = cppzarr.read_zarr(path, order="C", num_threads=n)
        else:
            start, stop = _region(input_bbox, fshape)
            a = cppzarr.read_zarr(path, start, stop, order="C", num_threads=n)
        return _squeeze_z(a)
    # MATLAB's layout: [ny, nx] or [ny, nx, nz] (y fastest)
    if len(fshape) not in (2, 3):
        raise ValueError(f"MATLAB-layout zarr arrays have 2 or 3 dimensions, not {fshape}")
    if whole:
        a = cppzarr.read_zarr(path, order="F", num_threads=n)
    else:
        start, stop = _region(input_bbox, fshape[2:] + fshape[:2])
        if len(fshape) == 3:
            # (z, y, x) -> the file's (y, x, z)
            start, stop = start[1:] + start[:1], stop[1:] + stop[:1]
        a = cppzarr.read_zarr(path, start, stop, order="F", num_threads=n)
    return _squeeze_z(_planes_from_matlab(a, n))


def _planes_from_matlab(a: np.ndarray, n: int) -> np.ndarray:
    """A Fortran ``(ny, nx[, nz])`` array as ``(nz, ny, nx)`` C order: each plane transposed in
    the array's own buffer when it can be (no second volume), else ``core.arrays.from_matlab``."""
    if a.ndim == 2:
        a = a[:, :, np.newaxis]
    ny, nx, nz = a.shape
    kernel = _inplace_kernel()
    if kernel is not None and a.flags.f_contiguous and a.flags.writeable:
        # each plane is nx rows of ny values (y fastest); turned into ny rows of nx
        flat = a.reshape(-1, order="F")
        kernel(flat, ny, nx, nz, n)
        return flat.reshape((nz, ny, nx))
    from petakit5d.core.arrays import from_matlab

    return from_matlab(a, num_threads=n)


# --- writezarr -----------------------------------------------------------------------------


def _region_write(path: str, data: np.ndarray, box, sparse_data: bool, uuid: bool, n: int) -> None:
    import cppzarr

    meta = zarr_meta(path)
    if zarr_layout(meta) == "yxz":
        raise ValueError(
            f"{path} is in MATLAB's layout ([ny, nx, nz], order F): the port writes native "
            "(z, y, x) arrays only; write a new array instead"
        )
    shape = [int(v) for v in meta["shape"]]
    start, stop = _region(box, shape)
    region = [e - s for s, e in zip(start, stop, strict=True)]
    dshape = [1] * max(0, len(region) - data.ndim) + list(data.shape)
    if dshape != region:
        raise ValueError(f"The box {list(box)} does not match the data size {list(data.shape)}")
    data = data.reshape(region)
    want = np.dtype(meta["dtype"] if "dtype" in meta else meta["data_type"]).newbyteorder("=")
    if data.dtype.newbyteorder("=") != want:
        data = _convert(data, want)
    kw = {}
    opts = _write_options()
    if "sparse" in opts:
        kw["sparse"] = bool(sparse_data)
    if "uuid" in opts and not uuid:
        kw["uuid"] = False
    cppzarr.write_zarr(path, data, start, stop, num_threads=n, **kw)


def write_zarr(
    data,
    filepath,
    *,
    overwrite: bool = False,
    block_size: Sequence[int] = _BLOCK,
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
    """Write a volume ``(z, y, x)`` or an image ``(y, x)`` as a native zarr array (MATLAB
    ``writezarr``).

    A missing array is created (``create_zarr`` with ``block_size``, in ``(z, y, x)``, capped at
    the size; an image is stored as ``[1, ny, nx]``, as MATLAB's createzarr does also when
    ``expand_2d_dim`` is false); the data is then written into the array (an existing one keeps
    its metadata and type, the data converted as MATLAB converts), or into its box ``bbox``
    (``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive). ``overwrite`` replaces an existing
    array (without ``bbox``). Arrays in MATLAB's layout are not written into.
    """
    del expand_2d_dim  # (MATLAB: createzarr expands a 2D size either way)
    a = np.asarray(data)
    if a.dtype.kind not in "iuf":
        raise TypeError(f"writezarr writes numeric data, not {a.dtype}")
    if a.ndim not in (2, 3):
        raise ValueError(f"writezarr writes an image (y, x) or a volume (z, y, x), not {a.ndim}D")
    if a.ndim == 3 and a.flags.f_contiguous and not a.flags.c_contiguous:
        raise ValueError(
            "writezarr takes (z, y, x) data; this array is Fortran-ordered (MATLAB's layout? "
            "use core.arrays.from_matlab)"
        )
    path = path_str(filepath)
    sz = list(a.shape)
    block = _ints(block_size, "blockSize")
    if len(sz) == 2:
        sz, block = _fit2d(block, sz, True)
    if len(block) != len(sz):
        raise ValueError(f"blockSize {block} does not fit the data size {sz}")
    block = [min(s, b) for s, b in zip(sz, block, strict=True)]
    if overwrite and os.path.isdir(path) and (bbox is None or len(bbox) == 0):
        shutil.rmtree(path)
    if bbox is None or len(bbox) == 0:
        bbox = (0, 0, 0, *sz)
    new_file = False
    if not os.path.isdir(path):
        if not create:
            raise FileNotFoundError(f"The zarr file {path} does not exist!")
        from petakit5d.core.arrays import matlab_class

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
    """Write one chunk of a native zarr array by its chunk index ``(iz, iy, ix)``, 0-based
    (MATLAB ``writeZarrBlock``, on the array's path instead of a ``blockedImage``); data past
    the array's edge is cut off. ``level`` and ``mode`` are accepted for compatibility."""
    del level, mode
    path = path_str(file_path)
    meta = zarr_meta(path)
    shape, chunks = [int(v) for v in meta["shape"]], [int(v) for v in meta["chunks"]]
    sub = _ints(block_sub, "blockSub")
    sub = [0] * (len(shape) - len(sub)) + sub if len(sub) < len(shape) else sub[: len(shape)]
    a = np.asarray(data)
    dshape = [1] * max(0, len(shape) - a.ndim) + list(a.shape)
    start = [s * c for s, c in zip(sub, chunks, strict=True)]
    stop = [min(st + d, n) for st, d, n in zip(start, dshape, shape, strict=True)]
    if any(e <= st for st, e in zip(start, stop, strict=True)):
        raise ValueError(f"block {sub} is outside the array of shape {shape}")
    keep = tuple(slice(0, e - st) for st, e in zip(start, stop, strict=True))
    _region_write(path, a.reshape(dshape)[keep], (*start, *stop), False, True, threads(num_threads))
