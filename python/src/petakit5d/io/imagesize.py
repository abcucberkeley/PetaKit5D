"""Image sizes, data types and zarr metadata without reading the data.

Ports getImageSize, getImageSizeBatch, getImageDataType, getZarrInfo and
getImageBoundingBox. Sizes are MATLAB's: ``(ny, nx, nz)`` (``(ny, nx, nz, samples)`` for
color tiff), data types MATLAB's class names (``'uint16'``, ``'single'``, ...).

Also the helpers ``tiff.py`` and ``zarr.py`` share: the default thread count, MATLAB's class
names and size rules. They stand in for ``petakit5d.core`` (task F02) until it lands.
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
    "ZarrInfo",
    "get_image_bounding_box",
    "get_image_dtype",
    "get_image_size",
    "get_image_size_batch",
    "get_zarr_info",
]

# --- stand-ins for petakit5d.core (F02) ----------------------------------------------------

# MATLAB class <-> NumPy type
MATLAB_TO_NUMPY = {
    "double": "float64",
    "single": "float32",
    "uint8": "uint8",
    "int8": "int8",
    "uint16": "uint16",
    "int16": "int16",
    "uint32": "uint32",
    "int32": "int32",
    "uint64": "uint64",
    "int64": "int64",
    "logical": "bool",
}
NUMPY_TO_MATLAB = {v: k for k, v in MATLAB_TO_NUMPY.items()}


def matlab_class(dtype) -> str:
    """MATLAB's class name of a NumPy type (``float32`` -> ``'single'``)."""
    name = np.dtype(dtype).newbyteorder("=").name
    try:
        return NUMPY_TO_MATLAB[name]
    except KeyError:
        raise ValueError(f"{name} has no MATLAB class") from None


def numpy_dtype(cls) -> np.dtype:
    """The NumPy type of a MATLAB class name (or of anything ``np.dtype`` takes)."""
    if isinstance(cls, str) and cls in MATLAB_TO_NUMPY:
        return np.dtype(MATLAB_TO_NUMPY[cls])
    return np.dtype(cls)


def default_threads() -> int:
    """The CPUs this process may run on (the affinity mask; F02 adds cgroup quotas)."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def threads(num_threads: int | None) -> int:
    """A call's thread count: the caller's, or every CPU the process may use."""
    if num_threads is None:
        return default_threads()
    if isinstance(num_threads, bool) or int(num_threads) != num_threads or num_threads < 1:
        raise ValueError(f"num_threads must be a whole number of at least 1, not {num_threads!r}")
    return int(num_threads)


def matlab_shape(shape: Sequence[int]) -> tuple[int, ...]:
    """A shape as MATLAB's ``size`` gives it: at least 2 axes, trailing singletons dropped."""
    shape = [int(s) for s in shape] + [1] * max(0, 2 - len(shape))
    while len(shape) > 2 and shape[-1] == 1:
        shape.pop()
    return tuple(shape)


def matlab_squeeze(a: np.ndarray) -> np.ndarray:
    """A view of ``a`` with MATLAB's shape (trailing singleton axes dropped; no copy)."""
    shape = matlab_shape(a.shape)
    if shape == a.shape:
        return a
    return a.reshape(shape, order="F" if a.flags.f_contiguous else "C")


def path_str(path) -> str:
    """A file name as text (``os.PathLike`` too; the I/O packages take ``str`` only)."""
    p = os.fspath(path)
    return p.decode() if isinstance(p, bytes) else p


def is_tiff_name(path: str) -> bool:
    # MATLAB: strcmp(filePath(end - 2 : end), 'tif') || strcmp(filePath(end - 3 : end), 'tiff')
    return path.endswith(("tif", "tiff"))


def is_zarr_name(path: str) -> bool:
    return path.endswith("zarr")


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
    images = _imagej_images(h.description)
    if z == 1 and h.compression == 1 and images:
        z = images
    if h.samples > 1:
        return (h.length, h.width, z, h.samples)
    return (h.length, h.width, z)


def _zarr_meta(path: str) -> dict:
    with open(os.path.join(path, ".zarray"), "rb") as f:
        return json.loads(f.read().decode())


def get_image_size(file_path) -> tuple[int, ...]:
    """Size of a tiff file or zarr array without reading the data (MATLAB ``getImageSize``).

    Tiff: ``(ny, nx, nz)``, with the number of slices of an ImageJ stack saved in one
    directory, and the samples as a fourth value for color images. Zarr: the ``shape`` of
    its ``.zarray``.
    """
    path = path_str(file_path)
    if is_tiff_name(path):
        return _tiff_size(path)
    if is_zarr_name(path) or os.path.isfile(os.path.join(path, ".zarray")):
        return tuple(int(s) for s in _zarr_meta(path)["shape"])
    raise ValueError(f"{path} is not a tiff file or a zarr array")


def get_image_size_batch(file_paths: Sequence, use_parpool: bool = True) -> np.ndarray:
    """Sizes of many files, one row each (MATLAB ``getImageSizeBatch``); all have the same
    number of dimensions. ``use_parpool`` reads them on threads (same result)."""
    paths = [path_str(p) for p in file_paths]
    if not paths:
        return np.zeros((0, 0), dtype=np.int64)
    if use_parpool and len(paths) > 64:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(min(16, default_threads() * 2)) as ex:
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
    """The fields of MATLAB's ``getZarrInfo`` struct (``blockSize`` is ``block_size``)."""

    size: tuple[int, ...]
    block_size: tuple[int, ...]
    ddtype: str
    order: str
    compressor: Any
    fill_value: Any
    filters: Any
    zarr_format: int
    dtype: str


_ZARR_CLASS = {"f8": "double", "f4": "single", "u2": "uint16", "u1": "uint8"}


def get_zarr_info(file_path) -> ZarrInfo:
    """The ``.zarray`` of a zarr array (MATLAB ``getZarrInfo``); ``dtype`` is the MATLAB
    class of the four types PetaKit5D writes, ``''`` (with a warning) for others."""
    path = path_str(file_path)
    if not (is_zarr_name(path) or os.path.isfile(os.path.join(path, ".zarray"))):
        raise FileNotFoundError(
            f"The input zarr file {path} does not exist or the .zarray file is missing!"
        )
    meta = _zarr_meta(path)
    ddtype = meta["dtype"]
    dtype = _ZARR_CLASS.get(ddtype[-2:], "")
    if not dtype:
        warnings.warn(
            f"Unknown/Unsupported zarr data type {ddtype}, set dtype as empty!", stacklevel=2
        )
    return ZarrInfo(
        size=tuple(int(s) for s in meta["shape"]),
        block_size=tuple(int(c) for c in meta["chunks"]),
        ddtype=ddtype,
        order=meta["order"],
        compressor=meta["compressor"],
        fill_value=meta["fill_value"],
        filters=meta["filters"],
        zarr_format=int(meta["zarr_format"]),
        dtype=dtype,
    )


# --- getImageBoundingBox -------------------------------------------------------------------


def _positive_along(a: np.ndarray, axis: int) -> np.ndarray:
    """``sum(I, other axes) > 0`` along one axis, as MATLAB computes it."""
    others = tuple(i for i in range(a.ndim) if i != axis)
    if a.dtype.kind in "bu":
        # sums of non-negative values are positive exactly when one value is
        return np.any(a != 0, axis=others) if a.ndim > 1 else a != 0
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
    # stand-in for MATLAB's single/double sum order (matlabcompat sums, N05)
    return np.sum(a, axis=others) > 0


def get_image_bounding_box(img: np.ndarray) -> tuple[int, ...]:
    """1-based inclusive box of the positive part of a 2D or 3D image (MATLAB
    ``getImageBoundingBox``): ``(ymin, xmin, zmin, ymax, xmax, zmax)``, or
    ``(ymin, xmin, ymax, xmax)`` for 2D; ``()`` when no value is positive."""
    a = np.asarray(img)
    nd = len(matlab_shape(a.shape))
    a3 = a.reshape(a.shape + (1,) * (3 - a.ndim), order="A") if a.ndim < 3 else a
    lo, hi = [], []
    for axis in range(3 if nd == 3 else 2):
        idx = np.flatnonzero(_positive_along(a3, axis))
        if idx.size:
            lo.append(int(idx[0]) + 1)
            hi.append(int(idx[-1]) + 1)
    return tuple(lo + hi)
