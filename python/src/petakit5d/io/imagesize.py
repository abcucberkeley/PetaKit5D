"""Image sizes, data types and zarr metadata without reading the data.

Ports getImageSize, getImageSizeBatch, getImageDataType, getZarrInfo and
getImageBoundingBox. Sizes are ``(nz, ny, nx)`` for tiff files and for zarr arrays of both
layouts (the port's native ``[nz, ny, nx]`` arrays and MATLAB's ``[ny, nx, nz]`` ones, order
F); data types are MATLAB's class names (``'uint16'``, ``'single'``, ...); boxes are the
Python API's ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive.

Also the helpers ``tiff.py`` and ``zarr.py`` share: thread counts, file names, the layout of a
zarr array.
"""

from __future__ import annotations

import json
import os
import struct
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "ReadBuffer",
    "ZarrInfo",
    "get_image_bounding_box",
    "get_image_dtype",
    "get_image_size",
    "get_image_size_batch",
    "get_zarr_info",
]

_AXES = ("auto", "zyx", "yxz")


def threads(num_threads: int | None) -> int:
    """A call's thread count: the caller's (at most the CPUs this process may use), or all of
    them for None."""
    from petakit5d.core.threads import resolve_threads

    if num_threads is not None and (
        isinstance(num_threads, bool) or int(num_threads) != num_threads or num_threads < 1
    ):
        raise ValueError(f"num_threads must be a whole number of at least 1, not {num_threads!r}")
    return resolve_threads(num_threads)


def path_str(path) -> str:
    """A file name as text (``os.PathLike`` too; the I/O packages take ``str`` only)."""
    p = os.fspath(path)
    return p.decode() if isinstance(p, bytes) else p


def is_tiff_name(path: str) -> bool:
    # MATLAB: strcmp(filePath(end - 2 : end), 'tif') || strcmp(filePath(end - 3 : end), 'tiff')
    return path.endswith(("tif", "tiff"))


def is_zarr_name(path: str) -> bool:
    return path.endswith("zarr")


def zarr_meta(path: str) -> dict:
    """The array metadata: ``.zarray`` (Zarr v2), else ``zarr.json`` (v3)."""
    for name in (".zarray", "zarr.json"):
        fn = os.path.join(path, name)
        if os.path.isfile(fn):
            with open(fn, "rb") as f:
                return json.loads(f.read().decode())
    raise FileNotFoundError(f"zarr file {path} does not exist (no .zarray or zarr.json)")


def zarr_layout(meta: dict, axes: str = "auto") -> str:
    """``'zyx'`` (the port's native arrays, shape ``[nz, ny, nx]``) or ``'yxz'`` (MATLAB's,
    shape ``[ny, nx, nz]``). ``'auto'``: MATLAB's when a Zarr v2 array has ``"order": "F"``."""
    if axes not in _AXES:
        raise ValueError(f"axes must be one of {_AXES}, not {axes!r}")
    if axes != "auto":
        return axes
    v2 = int(meta.get("zarr_format", 2)) == 2 and "dtype" in meta
    return "yxz" if v2 and meta.get("order") == "F" else "zyx"


def _zyx(v: Sequence[int], layout: str) -> tuple[int, ...]:
    # a shape or chunk shape of the file in (z, y, x)
    v = [int(s) for s in v]
    if layout == "yxz" and len(v) >= 2:
        from petakit5d.core.arrays import size_from_matlab

        return tuple(size_from_matlab(v, len(v)))
    return tuple(v)


# --- buffers of reads --------------------------------------------------------------------


class ReadBuffer:
    """One buffer for a loop of reads: ``read_zarr(..., out=buf)`` and ``read_tiff(...,
    out=buf)`` read into a C-contiguous view of its first bytes, so the blocks or frames of a
    loop reuse memory that is already touched (no page faults for a new output, no free after
    it). The buffer grows to the largest read; a returned array is valid until the next read
    into the same buffer."""

    def __init__(self) -> None:
        self._buf: np.ndarray | None = None

    def view(self, shape: Sequence[int], dtype) -> np.ndarray:
        """A C-contiguous array of ``shape`` and ``dtype`` on the buffer's first bytes."""
        dt = np.dtype(dtype)
        shape = tuple(int(v) for v in shape)
        n = 1
        for v in shape:
            n *= v
        if self._buf is None or self._buf.dtype != dt or self._buf.size < n:
            self._buf = None
            self._buf = np.empty(max(n, 1), dtype=dt)
        return self._buf[:n].reshape(shape)

    def release(self) -> None:
        """Free the buffer (the next read allocates a new one)."""
        self._buf = None


def out_array(out, shape: Sequence[int], dtype, what: str) -> tuple[np.ndarray, np.ndarray]:
    """The array a read fills (C-contiguous, of ``shape``) and the one it returns, for ``out``:
    a ``ReadBuffer`` (a view of it; an image of one slice returned 2D), or an array of the
    read's dtype, writeable and C-contiguous, of ``shape`` or, for one slice, ``shape[1:]``
    (returned as it is)."""
    shape = tuple(int(v) for v in shape)
    dt = np.dtype(dtype).newbyteorder("=")
    if isinstance(out, ReadBuffer):
        a = out.view(shape, dt)
        return a, (a[0] if len(shape) == 3 and shape[0] == 1 else a)
    if not isinstance(out, np.ndarray):
        raise TypeError(f"{what}: out must be a NumPy array or a ReadBuffer, not {type(out)}")
    if out.dtype != dt:
        raise TypeError(f"{what}: out has class {out.dtype}, the data {dt}")
    if not out.flags.writeable or not out.flags.c_contiguous:
        raise ValueError(f"{what}: out must be writeable and C-contiguous (z, y, x)")
    ok = [shape] + ([shape[1:]] if len(shape) == 3 and shape[0] == 1 else [])
    if out.shape not in ok:
        raise ValueError(f"{what}: out has shape {out.shape}, the read gives {shape}")
    return out.reshape(shape), out


# --- the first directory of a tiff file ----------------------------------------------------


@dataclass(frozen=True)
class _TiffHeader:
    width: int
    length: int
    bits: int
    samples: int
    sample_format: int
    compression: int
    description: bytes | None


# tag type -> (struct code, bytes)
_TIFF_TYPES = {
    1: ("B", 1),
    2: ("c", 1),
    3: ("H", 2),
    4: ("I", 4),
    6: ("b", 1),
    7: ("B", 1),
    8: ("h", 2),
    9: ("i", 4),
    16: ("Q", 8),
    17: ("q", 8),
}


def _tiff_header(path: str) -> _TiffHeader:
    """The tags of the first image directory that sizes and types depend on (classic TIFF
    and BigTIFF, either byte order), with libtiff's defaults for missing tags."""
    with open(path, "rb") as f:
        head = f.read(16)
        if len(head) < 8 or head[:2] not in (b"II", b"MM"):
            raise ValueError(f'File "{path}" cannot be opened (not a TIFF file)')
        e = "<" if head[:2] == b"II" else ">"
        magic = struct.unpack(e + "H", head[2:4])[0]
        if magic == 42:
            big, offset = False, struct.unpack(e + "I", head[4:8])[0]
        elif magic == 43 and len(head) >= 16:
            big, offset = True, struct.unpack(e + "Q", head[8:16])[0]
        else:
            raise ValueError(f'File "{path}" cannot be opened (not a TIFF file)')
        f.seek(offset)
        n = struct.unpack(e + ("Q" if big else "H"), f.read(8 if big else 2))[0]
        size = 20 if big else 12
        raw = f.read(n * size)
        tags: dict[int, tuple[int, int, bytes]] = {}
        for i in range(n):
            ent = raw[i * size : (i + 1) * size]
            if len(ent) < size:
                break
            tag, typ = struct.unpack(e + "HH", ent[:4])
            count = struct.unpack(e + ("Q" if big else "I"), ent[4:12] if big else ent[4:8])[0]
            tags[tag] = (typ, count, ent[12:20] if big else ent[8:12])

        def values(tag):
            if tag not in tags or tags[tag][0] not in _TIFF_TYPES:
                return None
            typ, count, field = tags[tag]
            code, nbytes = _TIFF_TYPES[typ]
            total = nbytes * count
            if total <= len(field):
                data = field[:total]
            else:
                f.seek(struct.unpack(e + ("Q" if big else "I"), field)[0])
                data = f.read(total)
            if typ == 2:
                return data
            return struct.unpack(e + code * count, data)

        def first(tag, default):
            v = values(tag)
            return int(v[0]) if v else default

        desc = values(270) if 270 in tags else None
        return _TiffHeader(
            width=first(256, 1),
            length=first(257, 1),
            bits=first(258, 1),
            samples=first(277, 1),
            sample_format=first(339, 1),
            compression=first(259, 1),
            description=desc,
        )


def _imagej_images(desc: bytes | None) -> int:
    # cpp-tiff: an 'ImageJ' description's 'images=' value (strtol: 0 if no digits follow)
    if not desc or b"ImageJ" not in desc:
        return 0
    k = desc.find(b"images=")
    if k < 0:
        return 0
    s = desc[k + 7 :].split(b"\0", 1)[0].lstrip(b" \t\n\r\f\v")
    sign = -1 if s[:1] == b"-" else 1
    s = s[1:] if s[:1] in (b"+", b"-") else s
    digits = b""
    for c in s:
        if not 48 <= c <= 57:
            break
        digits += bytes([c])
    return sign * int(digits) if digits else 0


# --- getImageSize, getImageSizeBatch -------------------------------------------------------


def _tiff_size(path: str) -> tuple[int, ...]:
    import cpptiff

    if not os.path.isfile(path):
        raise FileNotFoundError(f'File "{path}" cannot be opened')
    h = _tiff_header(path)
    z = int(cpptiff.get_image_shape(path)[0])
    # getImageSizeMex: an ImageJ stack of one uncompressed directory holds all its slices
    # (cpp-tiff 1.7.3 reports 1, T8)
    images = _imagej_images(h.description)
    if z == 1 and h.compression == 1 and images:
        z = images
    if h.samples > 1:
        return (z, h.length, h.width, h.samples)
    return (z, h.length, h.width)


def get_image_size(file_path) -> tuple[int, ...]:
    """Size of a tiff file or zarr array without reading the data (MATLAB ``getImageSize``).

    Tiff: ``(nz, ny, nx)`` (an image has ``nz = 1``), with the slices of an ImageJ stack saved
    in one directory, and the samples as a fourth value for color images:
    ``(nz, ny, nx, c)``. Zarr: the array's shape in ``(z, y, x)``, for native arrays and for
    MATLAB's (``[ny, nx, nz]``, order F).
    """
    path = path_str(file_path)
    if is_tiff_name(path):
        return _tiff_size(path)
    if is_zarr_name(path) or os.path.isfile(os.path.join(path, ".zarray")):
        # (MATLAB reads .zarray only; the port's sharded arrays have a zarr.json)
        meta = zarr_meta(path)
        return _zyx(meta["shape"], zarr_layout(meta))
    raise ValueError(f"{path} is not a tiff file or a zarr array")


def get_image_size_batch(file_paths: Sequence, use_parpool: bool = True) -> np.ndarray:
    """Sizes of many files, one row ``(nz, ny, nx)`` each (MATLAB ``getImageSizeBatch``); all
    have the same number of dimensions. ``use_parpool`` reads them on threads (same result)."""
    paths = [path_str(p) for p in file_paths]
    if not paths:
        return np.zeros((0, 0), dtype=np.int64)
    if use_parpool and len(paths) > 64:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(min(16, threads(None) * 2)) as ex:
            sizes = list(ex.map(get_image_size, paths))
    else:
        sizes = [get_image_size(p) for p in paths]
    n = len(sizes[0])
    for p, s in zip(paths, sizes, strict=False):
        if len(s) != n:
            raise ValueError(f"{p} has {len(s)} dimensions, the first file {n}")
    return np.array(sizes, dtype=np.int64).reshape(len(paths), n)


# --- getImageDataType, getZarrInfo ---------------------------------------------------------


def _tiff_class(h: _TiffHeader) -> str:
    # MATLAB's class of a tiff image: bit depth and SampleFormat (2 signed, 3 float)
    if h.bits == 1:
        return "logical"
    if h.bits not in (8, 16, 32, 64):
        raise ValueError(f"tiff images of {h.bits} bits are not supported")
    if h.sample_format == 3 and h.bits in (32, 64):
        return "single" if h.bits == 32 else "double"
    return ("int" if h.sample_format == 2 else "uint") + str(h.bits)


def get_image_dtype(file_path, zarr_file: bool = False) -> str:
    """MATLAB class name of a tiff or zarr image (MATLAB ``getImageDataType``).

    ``zarr_file`` treats a folder without the ``.zarr`` extension as a zarr array.
    """
    path = path_str(file_path)
    if is_tiff_name(path):
        if not os.path.isfile(path):
            raise FileNotFoundError(f'File "{path}" cannot be opened')
        return _tiff_class(_tiff_header(path))
    if is_zarr_name(path) or zarr_file:
        return get_zarr_info(path).dtype
    raise ValueError("Unknown data type, currently only tiff and zarr files are supported!")


@dataclass(frozen=True)
class ZarrInfo:
    """The fields of MATLAB's ``getZarrInfo`` struct (``blockSize`` is ``block_size``), with
    ``size`` and ``block_size`` in ``(z, y, x)`` and ``axes`` the layout of the file: ``'zyx'``
    (native) or ``'yxz'`` (MATLAB's)."""

    size: tuple[int, ...]
    block_size: tuple[int, ...]
    ddtype: str
    order: str
    compressor: Any
    fill_value: Any
    filters: Any
    zarr_format: int
    dtype: str
    axes: str


_ZARR_CLASS = {"f8": "double", "f4": "single", "u2": "uint16", "u1": "uint8"}


_V3_CLASS = {"float64": "double", "float32": "single", "uint16": "uint16", "uint8": "uint8"}


def get_zarr_info(file_path) -> ZarrInfo:
    """The ``.zarray`` of a zarr array (MATLAB ``getZarrInfo``); ``dtype`` is the MATLAB
    class of the four types PetaKit5D writes, ``''`` (with a warning) for others. A Zarr v3
    array (the port's sharded arrays) is described from its ``zarr.json``: ``block_size`` is
    its shard (chunk grid) shape, ``compressor`` its codecs."""
    path = path_str(file_path)
    v2 = os.path.isfile(os.path.join(path, ".zarray"))
    v3 = not v2 and os.path.isfile(os.path.join(path, "zarr.json"))
    if not (v2 or v3):
        raise FileNotFoundError(
            f"The input zarr file {path} does not exist or the .zarray file is missing!"
        )
    meta = zarr_meta(path)
    if v3:
        ddtype = str(meta["data_type"])
        dtype = _V3_CLASS.get(ddtype, "")
        chunks = meta["chunk_grid"]["configuration"]["chunk_shape"]
        order, compressor, filters = "C", meta.get("codecs"), None
    else:
        ddtype = meta["dtype"]
        dtype = _ZARR_CLASS.get(ddtype[-2:], "")
        chunks, order = meta["chunks"], meta["order"]
        compressor, filters = meta["compressor"], meta["filters"]
    if not dtype:
        warnings.warn(
            f"Unknown/Unsupported zarr data type {ddtype}, set dtype as empty!", stacklevel=2
        )
    layout = zarr_layout(meta)
    return ZarrInfo(
        size=_zyx(meta["shape"], layout),
        block_size=_zyx(chunks, layout),
        ddtype=ddtype,
        order=order,
        compressor=compressor,
        fill_value=meta["fill_value"],
        filters=filters,
        zarr_format=int(meta["zarr_format"]),
        dtype=dtype,
        axes=layout,
    )


# --- getImageBoundingBox -------------------------------------------------------------------


def _positive_along(a: np.ndarray, axis: int) -> np.ndarray:
    """``sum(I, other axes) > 0`` along one axis, as MATLAB computes it."""
    others = tuple(i for i in range(a.ndim) if i != axis)
    if a.dtype.kind in "bu":
        # sums of non-negative values are positive exactly when one value is
        return np.any(a != 0, axis=others)
    if a.dtype.kind == "i":
        # MATLAB sums integers in double: exact below 2^53, as int64 sums are (8 to 32 bits)
        if a.dtype.itemsize < 8:
            return np.sum(a, axis=others, dtype=np.int64) > 0
        return np.sum(a, axis=others, dtype=np.float64) > 0
    if a.dtype.kind == "f" and not (a < 0).any():
        # no negative values: positive unless a NaN is there, or no value is positive
        pos = np.any(a > 0, axis=others)
        nan = np.any(np.isnan(a), axis=others)
        return pos & ~nan
    # stand-in for MATLAB's order of single/double sums with negative values (no matlabcompat
    # sum yet): NumPy's pairwise sum
    return np.sum(a, axis=others) > 0


def get_image_bounding_box(img: np.ndarray) -> tuple[int, ...]:
    """The box of the positive part of a volume ``(z, y, x)`` or an image ``(y, x)`` (MATLAB
    ``getImageBoundingBox``): ``(z0, y0, x0, z1, y1, x1)`` (``(y0, x0, y1, x1)`` for an image),
    0-based, end-exclusive; ``()`` when no value is positive. A plane counts as positive when
    its sum is (MATLAB's rule)."""
    a = np.asarray(img)
    if a.ndim not in (2, 3):
        raise ValueError(f"an image (y, x) or a volume (z, y, x), not shape {a.shape}")
    lo, hi = [], []
    for axis in range(a.ndim):
        idx = np.flatnonzero(_positive_along(a, axis))
        if idx.size == 0:
            return ()
        lo.append(int(idx[0]))
        hi.append(int(idx[-1]) + 1)
    return tuple(lo + hi)
