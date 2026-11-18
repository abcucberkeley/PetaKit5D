"""N5 and zarr conversion and the BigStitcher XML: ``N5ToZarr``, ``N5ToZarrBlock``,
``zarrToN5``, ``zarrToN5Block``, ``constructXmlFile``, ``removeTextTags`` and ``cleanXmlFile``,
with ``writestruct`` for the XML they write.

Arrays are ``(z, y, x)`` C order; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive;
block, batch and subfolder sizes ``(z, y, x)``. The MATLAB-name layer converts MATLAB's values.

**N5 on disk is standard N5** (docs/decisions/0009-n5-layout.md, a deliberate difference from
MATLAB's files, as for zarr): the dataset of a ``(z, y, x)`` volume has ``dimensions``
``[nx, ny, nz]`` (x fastest) and ``blockSize`` ``[bx, by, bz]``, chunk files named by the
chunk's index in that order (``ix/iy/iz``), what Fiji, BigStitcher and zarr-python's or
tensorstore's N5 drivers read as the image. Otherwise the files are MATLAB's: every chunk
stored whole (edge chunks padded with zeros, their headers giving the full block), big-endian,
blosc zstd level 1 with byte shuffle, and ``attributes.json`` as zarr-python 2 writes it
(``"n5": "2.0.0"`` included: the dataset is the store's root). MATLAB's N5Adapter (over
zarr-python 2's ``N5Store``) writes MATLAB's ``(y, x, z)`` array in C order instead:
``dimensions`` ``[nz, nx, ny]``. Readers take ``axes='yxz'`` for such files; nothing in the
file tells the two apart. The BigStitcher XML of ``construct_xml_file`` describes the port's
files with ``data_order`` ``'zyx'`` (MATLAB's ``'yxz'`` describes MATLAB's). Reading takes N5
datasets of any of the standard codecs (raw, gzip, bzip2, xz, blosc, zstd) and chunks stored
whole or cut at the edges (as n5-java writes them).

**Codecs**: raw, gzip, bzip2 and xz use the standard library; blosc (MATLAB's default, so
``zarr_to_n5`` needs it) and zstd before Python 3.14 need numcodecs (MIT), the ``n5`` extra
(``pip install petakit5d[n5]``). Chunks are coded and decoded on threads (numcodecs and the
copies release the GIL).

``n5_to_zarr`` and ``zarr_to_n5`` cut the output into batches of whole chunks (MATLAB's
``XR_zarrChunkCoordinatesExtraction``, batches numbered y fastest, then x, then z), hand groups
of them to block tasks through the job framework (each task's flag marks it done), and move the
temporary result ``<result>_<uuid>`` into place; a result another writer published is kept.
Zarr results are the port's native arrays; zarr inputs of both layouts are read. Quirks kept:
MATLAB reads and writes only uint16 and single N5 data (other types raise MATLAB's errors), and
publishes the result even when tasks still failed after the second try (with a warning here).
"""

from __future__ import annotations

import json
import math
import os
import struct
import time
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "N5Array",
    "N5Error",
    "clean_xml_file",
    "construct_xml_file",
    "create_n5",
    "n5_to_zarr",
    "n5_to_zarr_block",
    "open_n5",
    "read_n5",
    "remove_text_tags",
    "write_n5",
    "write_struct_xml",
    "zarr_to_n5",
    "zarr_to_n5_block",
]

# zarr-python 2's N5 names of NumPy types
_N5_TYPES = {
    "uint8": "u1",
    "uint16": "u2",
    "uint32": "u4",
    "uint64": "u8",
    "int8": "i1",
    "int16": "i2",
    "int32": "i4",
    "int64": "i8",
    "float32": "f4",
    "float64": "f8",
}
# the classes MATLAB's N5Adapter maps (N5 type -> MATLAB class, MATLAB class -> N5 type)
_ADAPTER_CLASS = {"float32": "single", "uint16": "uint16"}
_ADAPTER_TYPE = {"single": "float32", "uint16": "uint16"}
# the compression MATLAB's N5Adapter writes (numcodecs.Blosc(cname='zstd', clevel=1) with
# MATLAB's double clevel, as zarr-python 2 writes it to attributes.json)
MATLAB_COMPRESSION = {"blocksize": 0, "clevel": 1.0, "cname": "zstd", "shuffle": 1, "type": "blosc"}
N5_VERSION = "2.0.0"
# the layouts of a volume in N5: the axes of the dataset's C-order array (its dimensions
# reversed). 'zyx': standard N5, dimensions [nx, ny, nz] (the port's); 'yxz': MATLAB's
# N5Adapter, MATLAB's (y, x, z) array in C order, dimensions [nz, nx, ny]
N5_AXES = ("zyx", "yxz")
_TO_VIEW = {"zyx": (0, 1, 2), "yxz": (1, 2, 0)}
_FROM_VIEW = {"zyx": (0, 1, 2), "yxz": (2, 0, 1)}
# GB a block task asks for: this many times its batch in single (MATLAB's)
BATCH_MEMORY_FACTOR = 2.5


class N5Error(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _print(msg: str, end: str = "\n") -> None:
    print(msg, end=end, flush=True)


def _ints(v: Any) -> list[int]:
    return [int(x) for x in np.asarray(v, dtype=np.float64).ravel()]


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    return np.size(v) == 0


# ---------------------------------------------------------------------------------------------
# codecs
# ---------------------------------------------------------------------------------------------


def _numcodecs(what: str):
    try:
        import numcodecs
    except ImportError as err:
        raise ImportError(
            f"N5 {what} needs numcodecs, the n5 extra: pip install petakit5d[n5]"
        ) from err
    return numcodecs


def _zstd_module():
    try:
        from compression import zstd  # Python 3.14

        return zstd
    except ImportError:
        return None


def _decompress(comp: dict, data: bytes) -> bytes:
    kind = comp.get("type", "raw")
    if kind == "raw":
        return data
    if kind == "gzip":
        import zlib

        # gzip or zlib streams (n5's useZlib)
        return zlib.decompress(data, 47)
    if kind == "bzip2":
        import bz2

        return bz2.decompress(data)
    if kind == "xz":
        import lzma

        return lzma.decompress(data)
    if kind == "blosc":
        return bytes(_numcodecs("blosc").blosc.decompress(data))
    if kind == "zstd":
        zstd = _zstd_module()
        if zstd is not None:
            return zstd.decompress(data)
        return bytes(_numcodecs("zstd").Zstd().decode(data))
    raise N5Error(f"N5 compression {kind!r} is not supported")


def _compress(comp: dict, block: np.ndarray) -> bytes:
    """the chunk's payload: ``block`` (big-endian, C order) compressed by ``comp``"""
    kind = comp.get("type", "raw")
    if kind == "raw":
        return block.tobytes()
    if kind == "gzip":
        level = int(comp.get("level", -1))
        level = 6 if level < 0 else level
        if comp.get("useZlib", False):
            import zlib

            return zlib.compress(block.tobytes(), level)
        import gzip

        return gzip.compress(block.tobytes(), compresslevel=level, mtime=0)
    if kind == "bzip2":
        import bz2

        return bz2.compress(block.tobytes(), int(comp.get("blockSize", 9)))
    if kind == "xz":
        import lzma

        return lzma.compress(block.tobytes(), preset=int(comp.get("preset", 6)))
    if kind == "blosc":
        nc = _numcodecs("blosc")
        codec = nc.Blosc(
            cname=str(comp.get("cname", "lz4")),
            clevel=int(comp.get("clevel", 5)),
            shuffle=int(comp.get("shuffle", 1)),
            blocksize=int(comp.get("blocksize", 0)),
        )
        # (the array, so blosc's typesize is the element size, as zarr-python passes it)
        return bytes(codec.encode(block))
    if kind == "zstd":
        level = int(comp.get("level", 3))
        zstd = _zstd_module()
        if zstd is not None:
            return zstd.compress(block.tobytes(), level=level)
        return bytes(_numcodecs("zstd").Zstd(level=level).encode(block.tobytes()))
    raise N5Error(f"N5 compression {kind!r} is not supported")


# ---------------------------------------------------------------------------------------------
# N5 datasets
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class N5Array:
    """An N5 dataset: ``dims`` and ``block`` in N5's order (fastest first), its N5
    ``data_type`` and ``compression``, and ``axes``, its layout (``N5_AXES``). ``shape`` and
    ``block_size`` are the ``(z, y, x)`` of the port's array (``dims`` ``[nx, ny, nz]`` for
    ``'zyx'``, ``[nz, nx, ny]`` for MATLAB's ``'yxz'``)."""

    path: str
    dims: tuple[int, ...]
    block: tuple[int, ...]
    data_type: str
    compression: dict
    axes: str = "zyx"

    @property
    def ndim(self) -> int:
        return len(self.dims)

    @property
    def view_shape(self) -> tuple[int, ...]:
        """the C-order array of the file: ``dims`` reversed"""
        return tuple(reversed(self.dims))

    @property
    def view_block(self) -> tuple[int, ...]:
        return tuple(reversed(self.block))

    @property
    def shape(self) -> tuple[int, ...]:
        return _from_view_shape(self.view_shape, self.axes)

    @property
    def block_size(self) -> tuple[int, ...]:
        return _from_view_shape(self.view_block, self.axes)

    @property
    def dtype(self) -> np.dtype:
        return np.dtype(_N5_TYPES[self.data_type])

    @property
    def matlab_class(self) -> str:
        """the class MATLAB's N5Adapter gives the data ('' for types it doesn't map)"""
        return _ADAPTER_CLASS.get(self.data_type, "")


def _check_axes(axes: str) -> str:
    if axes not in N5_AXES:
        raise ValueError(f"axes must be one of {N5_AXES}, got {axes!r}")
    return axes


def _from_view_shape(v: Sequence[int], axes: str) -> tuple[int, ...]:
    v = tuple(int(x) for x in v)
    return tuple(v[i] for i in _FROM_VIEW[axes]) if len(v) == 3 else v


def _to_view_shape(s: Sequence[int], axes: str) -> tuple[int, ...]:
    s = tuple(int(x) for x in s)
    return tuple(s[i] for i in _TO_VIEW[axes]) if len(s) == 3 else s


def _to_view(a: np.ndarray, axes: str) -> np.ndarray:
    return a.transpose(_TO_VIEW[axes]) if a.ndim == 3 and axes != "zyx" else a


def _attributes_text(meta: dict) -> str:
    """attributes.json as zarr-python 2 writes it"""
    return json.dumps(meta, indent=4, sort_keys=True, ensure_ascii=True, separators=(",", ": "))


def _dataset_path(path: str) -> str:
    """N5Adapter's openToRead: a folder with a .zgroup is read at its s0"""
    if os.path.isfile(os.path.join(path, ".zgroup")):
        if os.path.isdir(os.path.join(path, "s0")):
            return os.path.join(path, "s0")
        raise N5Error(f"Unable to open {path}")
    return path


def open_n5(path: str | os.PathLike, axes: str = "zyx") -> N5Array:
    """The N5 dataset at ``path`` (its attributes.json), read in the layout ``axes``: ``'zyx'``
    (standard N5, the port's) or ``'yxz'`` (MATLAB's N5Adapter). The layout can't be told
    from the file: N5 has no order flag, and MATLAB's attributes are those zarr-python 2
    writes for any C-order array."""
    _check_axes(axes)
    path = _dataset_path(os.fspath(path))
    try:
        with open(os.path.join(path, "attributes.json"), encoding="utf-8") as f:
            meta = json.load(f)
    except (OSError, ValueError) as err:
        raise N5Error(f"Unable to open {path}") from err
    if "dimensions" not in meta or "dataType" not in meta:
        raise N5Error(f"{path} is not an N5 dataset")
    comp = meta.get("compression")
    if comp is None:
        # N5 before 1.0
        comp = {"type": meta.get("compressionType", "raw")}
    dims = tuple(int(v) for v in meta["dimensions"])
    block = tuple(int(v) for v in meta.get("blockSize", dims))
    dt = str(meta["dataType"])
    if dt not in _N5_TYPES:
        raise N5Error(f"N5 data type {dt!r} is not supported")
    return N5Array(path, dims, block, dt, dict(comp), axes)


def create_n5(
    path: str | os.PathLike,
    size: Sequence[int],
    block_size: Sequence[int],
    dtype: Any = "uint16",
    *,
    compression: dict | None = None,
    axes: str = "zyx",
) -> N5Array:
    """Create an empty N5 dataset of the port's ``size`` and ``block_size`` (``(z, y, x)``):
    its attributes.json (as zarr-python 2 writes it), no chunks. Standard N5 (``axes``
    ``'zyx'``: dimensions ``[nx, ny, nz]``), or MATLAB's layout (``'yxz'``). ``dtype``: a
    MATLAB class or a NumPy type; ``compression``: MATLAB's blosc zstd by default. An existing
    folder is kept (an existing dataset's attributes are replaced, as zarr-python's mode 'w'
    does)."""
    _check_axes(axes)
    path = os.fspath(path)
    name = np.dtype(_matlab_to_numpy(dtype)).name
    if name not in _N5_TYPES:
        raise N5Error(f"N5 data type {name!r} is not supported")
    view = _to_view_shape(size, axes)
    vblock = _to_view_shape(
        [min(int(b), int(s)) for b, s in zip(block_size, size, strict=True)], axes
    )
    comp = dict(MATLAB_COMPRESSION if compression is None else compression)
    meta = {
        "blockSize": list(reversed(vblock)),
        "compression": comp,
        "dataType": name,
        "dimensions": list(reversed(view)),
        "n5": N5_VERSION,
    }
    os.makedirs(path, exist_ok=True)
    final = os.path.join(path, "attributes.json")
    tmp = f"{final}.{os.getpid()}.{time.monotonic_ns()}"
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        f.write(_attributes_text(meta))
    os.replace(tmp, final)
    return N5Array(path, tuple(reversed(view)), tuple(reversed(vblock)), name, comp, axes)


def _matlab_to_numpy(dtype: Any) -> Any:
    return {"single": "float32", "double": "float64", "logical": "uint8"}.get(dtype, dtype)


def _chunk_file(arr: N5Array, vidx: Sequence[int]) -> str:
    return os.path.join(arr.path, *[str(i) for i in reversed(vidx)])


def _read_chunk(arr: N5Array, vidx: Sequence[int]) -> np.ndarray | None:
    """the chunk at view index ``vidx`` as a big-endian C array of its stored view shape
    (whole or cut at the edge; callers copy it once, swapping the bytes as they go), None
    where it has no file"""
    try:
        with open(_chunk_file(arr, vidx), "rb") as f:
            raw = f.read()
    except FileNotFoundError:
        return None
    mode, nd = struct.unpack_from(">HH", raw, 0)
    if mode not in (0, 1):
        raise N5Error(f"N5 chunk mode {mode} is not supported")
    dims = struct.unpack_from(f">{nd}I", raw, 4)
    off = 4 + 4 * nd + (4 if mode == 1 else 0)
    n = math.prod(dims)
    if mode == 1:
        n = struct.unpack_from(">I", raw, 4 + 4 * nd)[0]
    data = _decompress(arr.compression, raw[off:])
    be = arr.dtype.newbyteorder(">")
    a = np.frombuffer(data, dtype=be, count=n)
    if mode == 1:
        # (a varlength chunk: the elements of the whole block, C order of the view)
        a = a[: math.prod(dims)]
    return a.reshape(tuple(reversed(dims)))


def _write_chunk(arr: N5Array, vidx: Sequence[int], block: np.ndarray) -> None:
    """write the whole chunk ``block`` (big-endian, C order, its full view block shape)"""
    be = np.ascontiguousarray(block, dtype=arr.dtype.newbyteorder(">"))
    head = struct.pack(f">HH{block.ndim}I", 0, block.ndim, *reversed(block.shape))
    payload = _compress(arr.compression, be)
    final = _chunk_file(arr, vidx)
    os.makedirs(os.path.dirname(final), exist_ok=True)
    tmp = f"{final}.{os.getpid()}.{time.monotonic_ns()}.partial"
    with open(tmp, "wb") as f:
        f.write(head)
        f.write(payload)
    os.replace(tmp, final)


def _view_box(
    box: Sequence[int] | None, shape: Sequence[int], axes: str
) -> tuple[list[int], list[int]]:
    """the view's (start, stop) of a box of the port's array (None: all)"""
    n = len(shape)
    if box is None or _empty(box):
        lo, hi = [0] * n, [int(s) for s in shape]
    else:
        b = _ints(box)
        if len(b) != 2 * n:
            raise ValueError(f"a box of {2 * n} values is needed, got {len(b)}")
        lo, hi = b[:n], b[n:]
        for i in range(n):
            if not 0 <= lo[i] <= hi[i] <= shape[i]:
                raise N5Error("The bounding box is out of the bound of the N5 dataset!")
    return list(_to_view_shape(lo, axes)), list(_to_view_shape(hi, axes))


def _chunks_of(lo: Sequence[int], hi: Sequence[int], vblock: Sequence[int]) -> list[tuple]:
    import itertools

    ranges = [
        range(lo[d] // vblock[d], (hi[d] - 1) // vblock[d] + 1) if hi[d] > lo[d] else range(0)
        for d in range(len(lo))
    ]
    return list(itertools.product(*ranges))


def _map(fn, items: list, num_threads: int | None) -> None:
    from petakit5d.core.threads import resolve_threads

    n = min(resolve_threads(num_threads), len(items))
    if n <= 1:
        for it in items:
            fn(it)
        return
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=n) as ex:
        for _ in ex.map(fn, items):
            pass


def read_n5(
    path: str | os.PathLike | N5Array,
    *,
    bbox: Sequence[int] | None = None,
    num_threads: int | None = None,
    out: np.ndarray | None = None,
    axes: str = "zyx",
) -> np.ndarray:
    """Read an N5 dataset, or its box ``bbox`` (``(z0, y0, x0, z1, y1, x1)``), as the port's
    ``(z, y, x)`` C array (MATLAB's N5Adapter ``getIORegion``); missing chunks read as 0.
    ``axes``: the layout of a path (``open_n5``; an ``N5Array`` has its own). ``out``: a
    C-contiguous array of the box's shape and the data's type to read into."""
    arr = path if isinstance(path, N5Array) else open_n5(path, axes)
    lo, hi = _view_box(bbox, arr.shape, arr.axes)
    shape = _from_view_shape([h - b for b, h in zip(lo, hi, strict=True)], arr.axes)
    if out is None:
        res = np.zeros(shape, dtype=arr.dtype)
    else:
        res = out
        if res.shape != shape or res.dtype != arr.dtype or not res.flags.c_contiguous:
            raise ValueError(f"out must be a C-contiguous {arr.dtype} array of shape {shape}")
    view = _to_view(res, arr.axes)
    vblock = arr.view_block

    def one(c: tuple) -> None:
        a = [c[d] * vblock[d] for d in range(arr.ndim)]
        s0 = [max(a[d], lo[d]) for d in range(arr.ndim)]
        s1 = [min(a[d] + vblock[d], hi[d]) for d in range(arr.ndim)]
        dst = view[tuple(slice(s0[d] - lo[d], s1[d] - lo[d]) for d in range(arr.ndim))]
        blk = _read_chunk(arr, c)
        if blk is None:
            if out is not None:
                dst[...] = 0
            return
        # (a chunk cut at the edge may hold less than the box needs: the rest is 0)
        e1 = [min(s1[d], a[d] + blk.shape[d]) for d in range(arr.ndim)]
        if out is not None and any(e1[d] < s1[d] for d in range(arr.ndim)):
            dst[...] = 0
        if all(e1[d] > s0[d] for d in range(arr.ndim)):
            src = blk[tuple(slice(s0[d] - a[d], e1[d] - a[d]) for d in range(arr.ndim))]
            dst[tuple(slice(0, e1[d] - s0[d]) for d in range(arr.ndim))] = src

    _map(one, _chunks_of(lo, hi, vblock), num_threads)
    return res


def write_n5(
    data: np.ndarray,
    path: str | os.PathLike | N5Array,
    *,
    start: Sequence[int] | None = None,
    num_threads: int | None = None,
    axes: str = "zyx",
) -> None:
    """Write the port's ``(z, y, x)`` array ``data`` into an existing N5 dataset at ``start``
    (``(z0, y0, x0)``; None: the origin), as MATLAB's N5Adapter ``setRegion`` does: every chunk
    the region touches is written whole (chunks only partly covered are read first). The data
    is cast to the dataset's type with MATLAB's rounding and saturation. ``axes``: the layout
    of a path (an ``N5Array`` has its own)."""
    from petakit5d.core.dtypes import cast

    arr = path if isinstance(path, N5Array) else open_n5(path, axes)
    if data.ndim == 2 and arr.ndim == 3:
        data = data[np.newaxis]
    if data.ndim != arr.ndim:
        raise ValueError(f"a {arr.ndim}D region is needed, got {data.ndim}D")
    if data.dtype != arr.dtype:
        data = cast(data, arr.dtype)
    s = [0] * arr.ndim if start is None else _ints(start)
    box = s + [s[d] + data.shape[d] for d in range(arr.ndim)]
    lo, hi = _view_box(box, arr.shape, arr.axes)
    view = _to_view(data, arr.axes)
    vshape, vblock = arr.view_shape, arr.view_block

    def one(c: tuple) -> None:
        a = [c[d] * vblock[d] for d in range(arr.ndim)]
        e = [min(a[d] + vblock[d], vshape[d]) for d in range(arr.ndim)]
        s0 = [max(a[d], lo[d]) for d in range(arr.ndim)]
        s1 = [min(e[d], hi[d]) for d in range(arr.ndim)]
        whole = all(s0[d] == a[d] and s1[d] == e[d] for d in range(arr.ndim))
        blk = np.zeros(vblock, dtype=arr.dtype.newbyteorder(">"))
        if not whole:
            old = _read_chunk(arr, c)
            if old is not None:
                k = tuple(slice(0, min(old.shape[d], vblock[d])) for d in range(arr.ndim))
                blk[k] = old[k]
        blk[tuple(slice(s0[d] - a[d], s1[d] - a[d]) for d in range(arr.ndim))] = view[
            tuple(slice(s0[d] - lo[d], s1[d] - lo[d]) for d in range(arr.ndim))
        ]
        _write_chunk(arr, c, blk)

    _map(one, _chunks_of(lo, hi, vblock), num_threads)


# ---------------------------------------------------------------------------------------------
# N5ToZarrBlock, zarrToN5Block
# ---------------------------------------------------------------------------------------------


def _check_block_task(inds: list[int], flag_fullname: str, overwrite: bool) -> bool:
    """the checks of the block functions; False: the task is done already"""
    flag_path = os.path.dirname(flag_fullname)
    if not os.path.isdir(flag_path):
        raise N5Error(f"The block directory {flag_path} does not exist, skip the processing!")
    if os.path.lexists(flag_fullname):
        if overwrite:
            os.remove(flag_fullname)
        else:
            _print(f"The batch files ({inds[0] + 1} - {inds[-1] + 1}) already exist, skip them!")
            return False
    return True


def _adapter_class(arr: N5Array) -> str:
    """N5Adapter's getInfo: the data's MATLAB class (MATLAB fails for the types it lacks)"""
    if not arr.matlab_class:
        raise N5Error('Unrecognized field name "Datatype".', "MATLAB:nonExistentField")
    return arr.matlab_class


@register(
    "N5ToZarrBlock",
    Param(
        "blockInds",
        "@isnumeric",
        "Batch numbers of the task (0-based).",
        py="batch_inds",
        convert="index1",
    ),
    Param("zarrFullpath", "@(x) ischar(x)", "Input N5 dataset.", py="n5_fullpath"),
    Param("outFullpath", "@(x) ischar(x)", "Output zarr (native, existing)."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "inBatchBBoxes",
        "@isnumeric",
        "Input box of each batch.",
        py="in_batch_boxes",
        convert="bboxes",
    ),
    Param(
        "outBatchBBoxes",
        "@isnumeric",
        "Output box of each batch.",
        py="out_batch_boxes",
        convert="bboxes",
    ),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    Param("flipEmptyValue", "@islogical", "Set 65535 (empty regions) to 0."),
    Param(
        "axes",
        None,
        "Layout of the N5 dataset: 'zyx' (standard N5) or 'yxz' (MATLAB's zarrToN5).",
        python_only=True,
    ),
    task="B09",
)
def n5_to_zarr_block(
    batch_inds: Sequence[int],
    n5_fullpath: str,
    out_fullpath: str,
    flag_fullname: str,
    in_batch_boxes: Sequence[Sequence[int]],
    out_batch_boxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
    flip_empty_value: bool = False,
    axes: str = "zyx",
) -> list[bool]:
    """Copy the input boxes of an N5 dataset (layout ``axes``) into the output boxes of an
    existing native zarr array (``N5ToZarrBlock``), with 65535 set to 0 (``flip_empty_value``),
    then write the task's flag. The flag's folder must exist (else MATLAB's error: the run is
    over)."""
    from petakit5d.io import write_zarr

    inds = _ints(batch_inds)
    flag_fullname = os.fspath(flag_fullname)
    if not _check_block_task(inds, flag_fullname, overwrite):
        return []
    if not os.path.isdir(n5_fullpath):
        raise N5Error(f"The input n5 file {n5_fullpath} doesnot exist!")
    if not os.path.isdir(out_fullpath):
        raise N5Error(f"The output zarr file {out_fullpath} doesnot exist!")
    arr = open_n5(n5_fullpath, axes)
    _adapter_class(arr)
    inb = [_ints(b) for b in in_batch_boxes]
    outb = [_ints(b) for b in out_batch_boxes]
    done = [False] * len(inds)
    for i, bi in enumerate(inds):
        _print(f"Process Batch {bi + 1}... ", end="")
        t0 = time.perf_counter()
        data = read_n5(arr, bbox=inb[i])
        if flip_empty_value:
            # in_batch .* cast(in_batch ~= 65535, dtype)
            data[data == 65535] = 0
        write_zarr(data, out_fullpath, bbox=outb[i], create=False)
        del data
        done[i] = True
        _print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.")
    if all(done):
        with open(flag_fullname, "w"):
            pass
    return done


@register(
    "zarrToN5Block",
    Param(
        "blockInds",
        "@isnumeric",
        "Batch numbers of the task (0-based).",
        py="batch_inds",
        convert="index1",
    ),
    Param("zarrFullpath", "@(x) ischar(x)", "Input zarr (either layout)."),
    Param("outFullpath", "@(x) ischar(x)", "Output N5 dataset (existing)."),
    Param("flagFullname", "@(x) ischar(x)", "Flag file written when all batches are done."),
    Param(
        "inBatchBBoxes",
        "@isnumeric",
        "Input box of each batch.",
        py="in_batch_boxes",
        convert="bboxes",
    ),
    Param(
        "outBatchBBoxes",
        "@isnumeric",
        "Output box of each batch.",
        py="out_batch_boxes",
        convert="bboxes",
    ),
    Param("Overwrite", "@islogical", "Redo batches whose flag exists."),
    task="B09",
)
def zarr_to_n5_block(
    batch_inds: Sequence[int],
    zarr_fullpath: str,
    out_fullpath: str,
    flag_fullname: str,
    in_batch_boxes: Sequence[Sequence[int]],
    out_batch_boxes: Sequence[Sequence[int]],
    *,
    overwrite: bool = False,
) -> list[bool]:
    """Copy the input boxes of a zarr array into the output boxes of an existing N5 dataset
    (``zarrToN5Block``), then write the task's flag. The flag's folder must exist (else
    MATLAB's error: the run is over)."""
    from petakit5d.io import ReadBuffer, read_zarr

    inds = _ints(batch_inds)
    flag_fullname = os.fspath(flag_fullname)
    if not _check_block_task(inds, flag_fullname, overwrite):
        return []
    if not os.path.isdir(zarr_fullpath):
        raise N5Error(f"The input zarr file {zarr_fullpath} doesnot exist!")
    if not os.path.isdir(out_fullpath):
        raise N5Error(f"The output n5 file {out_fullpath} doesnot exist!")
    arr = open_n5(out_fullpath)
    inb = [_ints(b) for b in in_batch_boxes]
    outb = [_ints(b) for b in out_batch_boxes]
    buf = ReadBuffer()
    done = [False] * len(inds)
    for i, bi in enumerate(inds):
        _print(f"Process Batch {bi + 1}... ", end="")
        t0 = time.perf_counter()
        data = read_zarr(zarr_fullpath, input_bbox=inb[i], out=buf)
        write_n5(data, arr, start=outb[i][:3])
        del data
        done[i] = True
        _print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.")
    if all(done):
        with open(flag_fullname, "w"):
            pass
    return done


# ---------------------------------------------------------------------------------------------
# N5ToZarr, zarrToN5
# ---------------------------------------------------------------------------------------------


def _batches(
    sz: Sequence[int], bbox: Any, block_size: Sequence[int], batch_size: Sequence[int]
) -> tuple[list[int], list[int], list[int], np.ndarray, np.ndarray]:
    """MATLAB's framework of both converters: (out size, block, batch, input boxes, output
    boxes), boxes one row per batch in MATLAB's order"""
    from petakit5d.core.blocks import zarr_chunk_coordinates

    if not _empty(bbox):
        b = _ints(bbox)
        start, out_size = b[:3], [b[i + 3] - b[i] for i in range(3)]
    else:
        start, out_size = [0, 0, 0], [int(v) for v in sz]
    block = [min(o, b) for o, b in zip(out_size, _ints(block_size), strict=True)]
    batch = [min(o, b) for o, b in zip(out_size, _ints(batch_size), strict=True)]
    batch = [math.ceil(b / k) * k for b, k in zip(batch, block, strict=True)]
    out_boxes, _, _ = zarr_chunk_coordinates(
        out_size, batch_size=batch, block_size=block, same_batch_size=False, border_size=0
    )
    out_boxes = np.asarray(out_boxes, dtype=np.int64)
    in_boxes = out_boxes + np.array(start + start, dtype=np.int64)
    return out_size, block, batch, in_boxes, out_boxes


def _paths(fullpath: str, result_dir_str: str, ext: str) -> tuple[str, str, str, str]:
    """(data path, file name, result folder, result) as MATLAB builds them"""
    from petakit5d.core.paths import fileparts

    data_path, fsname, _ = fileparts(fullpath)
    out_path = f"{data_path}/{result_dir_str}"
    return data_path, fsname, out_path, f"{out_path}/{fsname}{ext}"


def _result_folder(out_path: str) -> None:
    from petakit5d.io.tiff import group_writable

    if not os.path.isdir(out_path):
        os.makedirs(out_path, exist_ok=True)
        group_writable(out_path)


def _run(
    func: Any,
    fullpath: str,
    tmp_path: str,
    flag_path: str,
    in_boxes: np.ndarray,
    out_boxes: np.ndarray,
    task_min: int,
    batch: Sequence[int],
    common: dict,
    **kwargs: Any,
) -> list[bool]:
    """the tasks of both converters through the job framework, with MATLAB's second try"""
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task

    num_batch = len(out_boxes)
    task_size = max(task_min, math.ceil(num_batch / 5000))
    num_tasks = math.ceil(num_batch / task_size)
    tasks, outputs = [], []
    for t in range(num_tasks):
        inds = list(range(t * task_size, min((t + 1) * task_size, num_batch)))
        flag = f"{flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
        outputs.append(flag)
        tasks.append(
            Task.of(
                func,
                inds,
                fullpath,
                tmp_path,
                flag,
                in_boxes[inds].tolist(),
                out_boxes[inds].tolist(),
                overwrite=False,
                **kwargs,
            )
        )
    inputs = [fullpath] * num_tasks
    mem_allocate = math.prod(batch) * 4 / 2**30 * BATCH_MEMORY_FACTOR
    done = framework.run_tasks(inputs, outputs, tasks, mem_allocate=mem_allocate, **common)
    if not all(done):
        done = framework.run_tasks(inputs, outputs, tasks, mem_allocate=mem_allocate * 2, **common)
    return done


def _finish(tmp_path: str, result_fullpath: str, flag_path: str, done: list[bool]) -> None:
    import warnings

    from petakit5d.jobs.publish import publish_result

    if not all(done):
        # (MATLAB moves the result into place without a check)
        warnings.warn(
            f"{sum(not d for d in done)} / {len(done)} tasks failed; the result is incomplete",
            stacklevel=3,
        )
    publish_result(tmp_path, result_fullpath, flag_path)
    _print("Done!")


_CLUSTER_PARAMS = (
    Param("parseCluster", "@islogical", "Run the batches as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes batches too."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a batch job."),
    Param("uuid", "@ischar", "Id of the temporary result and flags (empty: a new one)."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
)


@register(
    "N5ToZarr",
    Param("n5Fullpath", None, "N5 dataset (or a folder with a .zgroup: its s0).", py="n5_fullpath"),
    Param("resultDirStr", "@ischar", "Result folder next to the N5 dataset."),
    Param(
        "bbox",
        "@isnumeric",
        "Region (z0, y0, x0, z1, y1, x1) to convert (None: all).",
        convert="bbox",
    ),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the zarr.", convert="size_yxz"),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of a block task.", convert="size_yxz"),
    Param(
        "zarrSubSize", "@isnumeric", "Subfolder sizes (z, y, x) of the zarr.", convert="size_yxz"
    ),
    Param("save16bit", "@islogical", "Not used."),
    Param("flipEmptyValue", "@islogical", "Set 65535 (BigStitcher's empty regions) to 0."),
    *_CLUSTER_PARAMS,
    Param(
        "axes",
        None,
        "Layout of the N5 dataset: 'zyx' (standard N5) or 'yxz' (MATLAB's zarrToN5).",
        python_only=True,
    ),
    task="B09",
)
def n5_to_zarr(
    n5_fullpath: str,
    *,
    result_dir_str: str = "zarr/",
    bbox: Sequence[int] | None = None,
    block_size: Sequence[int] = (256, 256, 256),
    batch_size: Sequence[int] = (512, 512, 512),
    zarr_sub_size: Sequence[int] = (20, 20, 20),
    save_16bit: bool = True,
    flip_empty_value: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 1,
    uuid: str = "",
    mcc_mode: bool = False,
    config_file: str = "",
    axes: str = "zyx",
) -> None:
    """Convert an N5 dataset (or its box ``bbox``) to a native zarr array
    ``<folder>/<result_dir_str>/<name>.zarr`` (``N5ToZarr``), in batches of whole chunks
    through the job framework. ``axes``: the dataset's layout, ``'zyx'`` (standard N5, as
    the port writes it) or ``'yxz'`` (MATLAB's zarrToN5); it can't be told from the file. An
    existing result is kept."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import mkdir_recursive
    from petakit5d.io import create_zarr, get_zarr_info
    from petakit5d.jobs.publish import remove_path

    del save_16bit  # (not used in MATLAB)
    n5_fullpath = os.fspath(n5_fullpath)
    if not uuid:
        uuid = get_uuid()
    _, fsname, out_path, result_fullpath = _paths(n5_fullpath, result_dir_str, ".zarr")
    _result_folder(out_path)
    if os.path.isdir(result_fullpath):
        _print("The output result exists, skip it!")
        return
    tmp_path = f"{result_fullpath.removesuffix('/')}_{uuid}"
    if not os.path.isdir(n5_fullpath):
        raise N5Error(f"The input file {n5_fullpath} does not exist, pleast double check!")
    arr = open_n5(n5_fullpath, axes)
    dtype = _adapter_class(arr)
    out_size, block, batch, in_boxes, out_boxes = _batches(arr.shape, bbox, block_size, batch_size)

    _print(f"Convert n5 file {n5_fullpath} to zarr ...")
    flag_path = f"{out_path}/zarr_flag/{fsname}_{uuid}/"
    if not os.path.isdir(flag_path):
        mkdir_recursive(flag_path)
    if os.path.isdir(tmp_path):
        info = get_zarr_info(tmp_path)
        if list(info.block_size) != block or list(info.size) != out_size:
            remove_path(tmp_path)
            remove_path(flag_path)
            os.makedirs(flag_path, exist_ok=True)
    if not os.path.isdir(tmp_path):
        create_zarr(
            tmp_path,
            data_size=out_size,
            block_size=block,
            dtype=dtype,
            zarr_sub_size=None if _empty(zarr_sub_size) else _ints(zarr_sub_size),
        )
    common = {
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "cpus_per_task": cpus_per_task,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = _run(
        n5_to_zarr_block,
        n5_fullpath,
        tmp_path,
        flag_path,
        in_boxes,
        out_boxes,
        10,
        batch,
        common,
        flip_empty_value=bool(flip_empty_value),
        axes=axes,
    )
    _finish(tmp_path, result_fullpath, flag_path, done)


@register(
    "zarrToN5",
    Param("zarrFullpath", None, "Zarr array (either layout)."),
    Param("resultDirStr", "@ischar", "Result folder next to the zarr."),
    Param(
        "bbox",
        "@isnumeric",
        "Region (z0, y0, x0, z1, y1, x1) to convert (None: all).",
        convert="bbox",
    ),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the N5 dataset.", convert="size_yxz"),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x) of a block task.", convert="size_yxz"),
    Param("zarrSubSize", "@isnumeric", "Not used.", convert="size_yxz"),
    Param("save16bit", "@islogical", "Not used."),
    *_CLUSTER_PARAMS,
    task="B09",
)
def zarr_to_n5(
    zarr_fullpath: str,
    *,
    result_dir_str: str = "n5/",
    bbox: Sequence[int] | None = None,
    block_size: Sequence[int] = (256, 256, 256),
    batch_size: Sequence[int] = (512, 512, 512),
    zarr_sub_size: Sequence[int] = (20, 20, 20),
    save_16bit: bool = True,
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 1,
    uuid: str = "",
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Convert a zarr array (or its box ``bbox``) to a standard N5 dataset (dimensions
    ``[nx, ny, nz]``; MATLAB writes ``[nz, nx, ny]``, docs/decisions/0009-n5-layout.md),
    ``<folder>/<result_dir_str>/<name>.n5`` (``zarrToN5``), in batches of whole chunks through
    the job framework. Only uint16 and single data (MATLAB's). An existing result is kept."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import mkdir_recursive
    from petakit5d.io import get_image_dtype, get_image_size
    from petakit5d.jobs.publish import remove_path

    del zarr_sub_size, save_16bit  # (not used in MATLAB)
    zarr_fullpath = os.fspath(zarr_fullpath)
    if not uuid:
        uuid = get_uuid()
    _, fsname, out_path, result_fullpath = _paths(zarr_fullpath, result_dir_str, ".n5")
    _result_folder(out_path)
    if os.path.isdir(result_fullpath):
        _print("The output result exists, skip it!")
        return
    tmp_path = f"{result_fullpath.removesuffix('/')}_{uuid}"
    if not os.path.isdir(zarr_fullpath):
        raise N5Error(f"The input file {zarr_fullpath} does not exist, pleast double check!")
    dtype = get_image_dtype(zarr_fullpath, True)
    sz = get_image_size(zarr_fullpath)
    if len(sz) == 2:
        sz = (1, *sz)
    out_size, block, batch, in_boxes, out_boxes = _batches(sz, bbox, block_size, batch_size)

    _print(f"Convert zarr file {zarr_fullpath} to n5 ...")
    flag_path = f"{out_path}/zarr_flag/{fsname}_{uuid}/"
    if not os.path.isdir(flag_path):
        mkdir_recursive(flag_path)
    if os.path.isdir(tmp_path):
        try:
            old = open_n5(tmp_path)
            same = list(old.block_size) == block and list(old.shape) == out_size
        except N5Error:
            same = False
        if not same:
            remove_path(tmp_path)
            remove_path(flag_path)
            os.makedirs(flag_path, exist_ok=True)
    if not os.path.isdir(tmp_path):
        if dtype not in _ADAPTER_TYPE:
            # (N5Adapter's openToWrite maps single and uint16 only)
            raise N5Error("Unrecognized function or variable 'dtype'.", "MATLAB:UndefinedFunction")
        create_n5(tmp_path, out_size, block, _ADAPTER_TYPE[dtype])
    common = {
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "cpus_per_task": cpus_per_task,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = _run(
        zarr_to_n5_block, zarr_fullpath, tmp_path, flag_path, in_boxes, out_boxes, 5, batch, common
    )
    _finish(tmp_path, result_fullpath, flag_path, done)


# ---------------------------------------------------------------------------------------------
# the BigStitcher XML: constructXmlFile, writestruct, removeTextTags, cleanXmlFile
# ---------------------------------------------------------------------------------------------


def _sort_index(s: str) -> list[int]:
    """``[~, i] = sort(s)`` of a char row: 0-based, stable"""
    return sorted(range(len(s)), key=lambda i: ord(s[i]))


def _num_text(v: Any) -> str:
    """a number as writestruct writes it"""
    f = float(v)
    if f.is_integer() and abs(f) < 1e15:
        return str(int(f))
    return repr(f)


@register(
    "constructXmlFile",
    Param("dataPath", None, "Folder of the tiles (ends with a separator)."),
    Param("fileName", None, "The N5 container the XML names (its name and extension are used)."),
    Param("csvFile", None, "Image list (csv) of the tiles."),
    Param("voxelSize", None, "Voxel size (z, y, x).", convert="per_axis_xyz"),
    Param("ChannelPatterns", None, "Channel identifiers of the tiles to include."),
    Param("zarrFile", None, "The tiles are zarr (else tiff)."),
    Param(
        "axisOrder",
        None,
        "Axis order of the stage coordinates ('xyz', '-x,y,z', ...).",
        convert="same",
    ),
    Param(
        "dataOrder",
        None,
        "Axes of the N5 datasets' C-order arrays: 'zyx' (standard N5, zarr_to_n5's), 'yxz' "
        "(MATLAB's zarrToN5).",
        convert="same",
        differs="the port writes standard N5 (dimensions [nx, ny, nz]), MATLAB [nz, nx, ny]",
    ),
    Param("xcorrInfoFn", None, "Stitching xcorr result (.mat with dxyz_shift; '': none)."),
    task="B09",
)
def construct_xml_file(
    data_path: str,
    file_name: str,
    csv_file: str,
    voxel_size: Sequence[float],
    channel_patterns: Sequence[str] = ("CamA",),
    zarr_file: bool = False,
    axis_order: str = "xyz",
    data_order: str = "zyx",
    xcorr_info_fn: str = "",
) -> tuple[dict, Any, tuple[float, ...]]:
    """The BigStitcher XML of the tiles of an image list (``constructXmlFile``): returns the
    struct ``s`` (nested dicts and lists in MATLAB's field order; None is MATLAB's
    ``missing``) for ``write_struct_xml``, the tiles' table and the voxel size along the N5
    dimensions. Tiles are sorted by camera, channel, z, y, x, t and iteration and kept if they
    match a channel pattern and exist in ``data_path``; their stage coordinates go through
    ``axis_order`` and start at 0, and ``data_order`` (``'zyx'`` for standard N5, as
    ``zarr_to_n5`` writes it; ``'yxz'`` for MATLAB's zarrToN5) maps them, the voxel size and
    the tile sizes to the N5 dimensions; ``xcorr_info_fn`` adds the stitching shifts."""
    from petakit5d.core.paths import fileparts, matlab_dir
    from petakit5d.core.text import sprintf
    from petakit5d.io import get_image_size
    from petakit5d.stitch import imagelist as il

    data_path = os.fspath(data_path)
    channel_patterns = (
        [channel_patterns] if isinstance(channel_patterns, str) else list(channel_patterns)
    )
    vxyz = [float(v) for v in np.asarray(voxel_size, dtype=np.float64).ravel()]
    if len(vxyz) == 3:
        vxyz = vxyz[::-1]  # (z, y, x) -> MATLAB's [x y z]

    tab = il.read_table(csv_file)
    fn = list(tab["Filename"])
    ext = ".zarr" if zarr_file else ".tif"
    fsn = [fileparts(s)[1] for s in fn]
    use_map_fname = "mappedFilename" in tab
    if use_map_fname:
        fsn = [fileparts(s)[1] for s in tab["mappedFilename"]]
    cam_rx = il.matlab_regex(r"_Cam\w_ch")
    ch_rx = il.matlab_regex(r"_ch[0-9]_")
    if all(cam_rx.search(s) for s in fsn):
        expression, specify_cam = il.TILE_EXPRESSION_CAM, True
    elif all(ch_rx.search(s) for s in fsn):
        expression, specify_cam = il.TILE_EXPRESSION_NOCAM, False
    else:
        raise N5Error("Unrecognized function or variable 'expression'.", "MATLAB:UndefinedFunction")
    tab = il._parse_names(tab, fsn, expression, specify_cam)
    if zarr_file:
        tab["Filename"] = [s.replace(".tif", ext) for s in tab["Filename"]]
    tab, tinds = il.sortrows(tab, ("camera", "ch", "z", "y", "x", "t", "fullIter"))
    # (MATLAB indexes the names of all rows with the sorted rows' indices)
    fsn = [fsn[i] for i in tinds]

    image_fnames = [e.name for e in matlab_dir(f"{data_path}/*{ext}")]
    if not image_fnames:
        raise N5Error("The tiles do not exist!")
    if use_map_fname:
        full = [f"{data_path}/{x}{ext}" for x in fsn]
        include = il._channel_filter(full, channel_patterns)
        if not include.all():
            tab = tab[include]
    else:
        full = [f"{data_path}/{x}" for x in image_fnames]
        include = il._channel_filter(full, channel_patterns)
        image_fnames = [x for x, k in zip(image_fnames, include, strict=True) if k]
    exist = il._exists_in(image_fnames, [fileparts(s)[1] for s in tab["Filename"]])
    if not exist.all():
        il._warn("Some files in the image list file do not exist! Ignore them in the stitching")
        for s in tab["Filename"][~exist]:
            print(f"    {{'{s}'}}")
        print()
        tab = tab[exist]

    # coordinates in the axis order, from 0
    coords = np.stack([tab["StageX_um_"], tab["StageY_um_"], tab["StageZ_um_"]], axis=1)
    axis_order = axis_order.replace(" ", "")
    pure = axis_order.replace(",", "").replace("-", "")
    order = _sort_index(pure)
    if len(order) != 3:
        raise N5Error(
            "Unable to perform assignment because the size of the left side is 1-by-3 and the "
            f"size of the right side is 1-by-{len(order)}.",
            "MATLAB:subsassigndimmismatch",
        )
    signs = np.array([-1.0 if f"-{c}" in axis_order else 1.0 for c in "xyz"])
    xyz = coords[:, order] * signs
    xyz = xyz - _col_min(xyz)

    dpure = data_order.replace(",", "").replace("-", "")
    dom = _sort_index(dpure)
    flip = dom[::-1]
    dxyz_order = {"yxz": [0, 1, 2], "zyx": [2, 0, 1]}.get(data_order)
    xyz_xml = xyz[:, flip]
    xyz_xml = xyz_xml - _col_min(xyz_xml)
    vs_xml = [vxyz[i] for i in flip]

    shift_xml = None
    if xcorr_info_fn:
        from petakit5d.io import read_mat

        a = read_mat(xcorr_info_fn)
        if dxyz_order is None:
            raise N5Error(
                "Unrecognized function or variable 'dxyz_order_mat'.", "MATLAB:UndefinedFunction"
            )
        shift = np.atleast_2d(np.asarray(a["dxyz_shift"], dtype=np.float64))[:, dxyz_order]
        shift_xml = shift[:, flip]

    n = tab.height
    names = list(tab["Filename"])
    s: dict[str, Any] = {"versionAttribute": 0.2}
    s["BasePath"] = {"typeAttribute": "relative", "Text": "."}
    _, fstem, fext = fileparts(os.fspath(file_name))
    setups = []
    for i in range(n):
        size = get_image_size(data_path + names[i])
        size = tuple(size) if len(size) == 3 else (1, *size[-2:])
        # the tile's size along the N5 dimensions: its C-order view in data_order, reversed
        axes = dict(zip("zyx", size, strict=True))
        dims = [axes[c] for c in reversed(dpure)] if sorted(dpure) == ["x", "y", "z"] else []
        if not dims:
            raise N5Error(f"dataOrder {data_order!r} must name x, y and z")
        setups.append(
            {
                "id": i,
                "name": i,
                "size": sprintf("%d %d %d", *dims),
                "voxelSize": {"unit": "um", "size": sprintf("%f %f %f", *vs_xml)},
            }
        )
    s["SequenceDescription"] = {
        "ImageLoader": {
            "formatAttribute": "bdv.n5",
            "versionAttribute": 1.0,
            "n5": {"typeAttribute": "relative", "Text": fstem + fext},
        },
        "ViewSetups": {
            "ViewSetup": setups,
            "Attributes": [
                {
                    "nameAttribute": "illumination",
                    "Illumination": {"id": 0, "name": 0},
                    "Channel": None,
                    "Tile": None,
                    "Angle": None,
                },
                {
                    "nameAttribute": "channel",
                    "Illumination": None,
                    "Channel": {"id": 0, "name": 0},
                    "Tile": None,
                    "Angle": None,
                },
                {
                    "nameAttribute": "tile",
                    "Illumination": None,
                    "Channel": None,
                    "Tile": [{"id": i, "name": i} for i in range(n)],
                    "Angle": None,
                },
                {
                    "nameAttribute": "angle",
                    "Illumination": None,
                    "Channel": None,
                    "Tile": None,
                    "Angle": {"id": 0, "name": 0},
                },
            ],
        },
        "Timepoints": {"typeAttribute": "pattern", "integerpattern": 0},
        "MissingViews": "",
    }
    regs = []
    ident = "1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0"
    for i in range(n):
        t = [xyz_xml[i, k] / vs_xml[k] for k in range(3)]
        vt = [
            {
                "typeAttribute": "affine",
                "Name": "Translation to Regular Grid",
                "affine": sprintf("1.0 0.0 0.0 %f 0.0 1.0 0.0 %f 0.0 0.0 1.0 %f", *t),
            },
            {"typeAttribute": "affine", "Name": "calibration", "affine": ident},
        ]
        if shift_xml is not None:
            if i >= shift_xml.shape[0]:
                raise N5Error(
                    "Index in position 1 exceeds array bounds. Index must not exceed "
                    f"{shift_xml.shape[0]}.",
                    "MATLAB:badsubscript",
                )
            d = [shift_xml[i, k] / vs_xml[k] for k in range(3)]
            vt.append(
                {
                    "typeAttribute": "affine",
                    "Name": "Stitching Transform",
                    "affine": sprintf("1.0 0.0 0.0 %f 0.0 1.0 0.0 %f 0.0 0.0 1.0 %f", *d),
                }
            )
        regs.append({"timepointAttribute": 0, "setupAttribute": i, "ViewTransform": vt})
    s["ViewRegistrations"] = {"ViewRegistration": regs}
    s["ViewInterestPoints"] = ""
    s["BoundingBoxes"] = ""
    s["PointSpreadFunctions"] = ""
    s["IntensityAdjustments"] = ""
    return s, tab, tuple(vs_xml)


def _col_min(a: np.ndarray) -> np.ndarray:
    """``min(a, [], 1)``: NaN left out (NaN where a column has only NaN)"""
    if a.shape[0] == 0:
        return np.zeros((0, a.shape[1]))
    return np.fmin.reduce(a, axis=0)


def _xml_escape(s: str, attr: bool = False) -> str:
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return s.replace('"', "&quot;") if attr else s


def _xml_text(v: Any) -> str:
    """a value's text, not escaped"""
    if isinstance(v, (str, np.str_)):
        return str(v)
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    return _num_text(v)


def _xml_lines(name: str, v: Any, depth: int, lines: list[str]) -> None:
    pad = "    " * depth
    if v is None:
        return
    if isinstance(v, (list, tuple)):
        for e in v:
            _xml_lines(name, e, depth, lines)
        return
    if not isinstance(v, dict):
        lines.append(f"{pad}<{name}>{_xml_escape(_xml_text(v))}</{name}>")
        return
    attrs = "".join(
        f' {k[: -len("Attribute")]}="{_xml_escape(_xml_text(x), True)}"'
        for k, x in v.items()
        if k.endswith("Attribute") and x is not None
    )
    children = [(k, x) for k, x in v.items() if not k.endswith("Attribute") and x is not None]
    if not children:
        lines.append(f"{pad}<{name}{attrs}/>")
        return
    lines.append(f"{pad}<{name}{attrs}>")
    for k, x in children:
        _xml_lines(k, x, depth + 1, lines)
    lines.append(f"{pad}</{name}>")


def write_struct_xml(
    s: dict, filename: str | os.PathLike, struct_node_name: str = "struct"
) -> None:
    """MATLAB's ``writestruct(s, filename, "StructNodeName", struct_node_name)`` for the
    structs of ``construct_xml_file``: an XML declaration, one element per field (fields ending
    in ``Attribute`` are the element's attributes), struct arrays as repeated elements, values
    as text (numbers as MATLAB prints them), ``missing`` (None) left out, 4-space indents."""
    lines = ['<?xml version="1.0" encoding="UTF-8"?>']
    _xml_lines(struct_node_name, s, 0, lines)
    with open(filename, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


_STRIP = " \t\n\r\x0b\x0c"


@register("removeTextTags", Param("fileName", None, "XML file, changed in place."), task="B09")
def remove_text_tags(file_name: str) -> None:
    """Remove the ``<Text>`` tags of the lines that hold only a ``<Text>`` element (the text
    stays on its line) (``removeTextTags``)."""
    with open(file_name, encoding="utf-8", newline="") as f:
        text = f.read()
    lines = text.split("\n")
    if lines and lines[-1] == "" and text.endswith("\n"):
        lines = lines[:-1]
    lines = [ln[:-1] if ln.endswith("\r") else ln for ln in lines]
    out = []
    for ln in lines:
        t = ln.strip(_STRIP)
        if t.startswith("<Text>") and t.endswith("</Text>"):
            ln = ln.replace("<Text>", "").replace("</Text>", "")
        out.append(ln)
    with open(file_name, "w", encoding="utf-8", newline="\n") as f:
        f.write("".join(ln + "\n" for ln in out))


@register("cleanXmlFile", Param("fileName", None, "XML file, changed in place."), task="B09")
def clean_xml_file(file_name: str) -> None:
    """Remove the white space around the texts of ``<BasePath type="relative">`` and
    ``<n5 type="relative">`` (``cleanXmlFile``)."""
    from petakit5d.core.text import matlab_strtrim

    with open(file_name, "rb") as f:
        text = f.read().decode("utf-8")
    tags = ('<BasePath type="relative">', "</BasePath>", '<n5 type="relative">', "</n5>")
    for i in range(0, len(tags), 2):
        open_tag, close_tag = tags[i], tags[i + 1]
        starts = _find_all(text, open_tag)
        ends = _find_all(text, close_tag)
        for j, st in enumerate(starts):
            if j >= len(ends):
                raise N5Error(
                    "Index exceeds the number of array elements. Index must not exceed "
                    f"{len(ends)}.",
                    "MATLAB:badsubscript",
                )
            # (indices found before the loop, as MATLAB's strfind gave them)
            between = text[st + len(open_tag) : ends[j]]
            # (MATLAB's '\n' and '\t' in single quotes are backslash pairs)
            cleaned = matlab_strtrim(between.replace("\\n", "").replace("\\t", ""))
            if between:
                text = text.replace(between, cleaned)
    with open(file_name, "wb") as f:
        f.write(text.encode("utf-8"))


def _find_all(text: str, sub: str) -> list[int]:
    """``strfind``: 0-based starts, overlapping matches included"""
    out, i = [], text.find(sub)
    while i >= 0:
        out.append(i)
        i = text.find(sub, i + 1)
    return out
