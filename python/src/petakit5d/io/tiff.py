"""Tiff reading and writing through cpp-tiff (ports readtiff, readtiff_parallel, writetiff).

Volumes are ``(y, x, z)`` Fortran order; one page comes back 2D, as MATLAB's ``size`` says.
The files are the mex's: BigTIFF, ``RowsPerStrip`` 512, LZW (the default), zstd or none.

cpp-tiff 1.7.3 reads and writes ``(z, y, x)`` C order. Until it takes ``order='F'``
(requested; detected at run time, once per process):

- reads are followed by a threaded page transpose in cpp-tiff's own buffer
  (``petakit5d._ext.io``; a NumPy copy when the extension is not built);
- writes pass the Fortran memory to the extension's ``transpose=True`` entry, which writes
  the same bytes as ``write_tiff`` without a copy (decision 0004).
"""

from __future__ import annotations

import inspect
import os
import stat
import warnings
from collections.abc import Sequence
from functools import cache

import numpy as np

from petakit5d.io.imagesize import (
    _imagej_images,
    _tiff_header,
    matlab_squeeze,
    path_str,
    threads,
)

__all__ = ["read_tiff", "read_tiff_parallel", "write_tiff"]

# the types cpp-tiff and MATLAB both have
_TIFF_TYPES = frozenset(
    np.dtype(t)
    for t in (
        "uint8",
        "int8",
        "uint16",
        "int16",
        "uint32",
        "int32",
        "uint64",
        "int64",
        "float32",
        "float64",
    )
)
_COMPRESSIONS = ("none", "lzw", "zstd")
_MODES = ("parallel", "libtiff", "imwrite")


# --- which entries the installed cpp-tiff has (decided once) -------------------------------


def _has_param(func, name: str) -> bool:
    try:
        return name in inspect.signature(func).parameters
    except (TypeError, ValueError):
        return False


@cache
def _reader() -> str:
    """'order' when cpp-tiff reads (y, x, z) itself, else 'transpose'."""
    import cpptiff

    return "order" if _has_param(cpptiff.read_tiff, "order") else "transpose"


@cache
def _writer() -> str:
    """'order' (cpp-tiff takes (y, x, z)), 'direct' (the extension's transpose=True entry)
    or 'copy' (the public write_tiff of a (z, y, x) view, which copies)."""
    import cpptiff

    if _has_param(cpptiff.write_tiff, "order"):
        return "order"
    try:
        from cpptiff.cpptiff import pybind11_write_tiff
    except ImportError:
        return "copy"
    doc = pybind11_write_tiff.__doc__ or ""
    return "direct" if "transpose" in doc else "copy"


@cache
def _inplace_kernel():
    """The in-place page transpose of petakit5d._ext.io, or None (NumPy copy instead)."""
    from petakit5d import _ext

    if _ext.available("io"):
        return _ext.load("io").transpose_pages_inplace
    warnings.warn(
        "petakit5d._ext.io is not built: tiff reads use a slower NumPy copy",
        RuntimeWarning,
        stacklevel=3,
    )
    return None


def _transpose_numpy(src: np.ndarray, dst: np.ndarray, num_threads: int) -> None:
    # (x, y, z) F -> (y, x, z) F, pages split among threads
    from concurrent.futures import ThreadPoolExecutor

    nz = src.shape[2]
    n = max(1, min(num_threads, nz))
    edges = [nz * i // n for i in range(n + 1)]

    def part(i):
        dst[:, :, edges[i] : edges[i + 1]] = src[:, :, edges[i] : edges[i + 1]].transpose(1, 0, 2)

    if n == 1:
        part(0)
    else:
        with ThreadPoolExecutor(n) as ex:
            list(ex.map(part, range(n)))


# --- readtiff ------------------------------------------------------------------------------


def _num_pages(path: str) -> int:
    # the bound the mex checks ranges against: the directories, or an ImageJ stack's slices
    import cpptiff

    h = _tiff_header(path)
    z = int(cpptiff.get_image_shape(path)[0])
    images = _imagej_images(h.description)
    if z == 1 and h.compression == 1 and images:
        return images
    return z


def _check_range(path: str, rng) -> list[int]:
    r = list(np.asarray(rng, dtype=np.float64).ravel())
    if len(r) != 2:
        raise ValueError("Input range is not 2")
    if any(not np.isfinite(v) or v != int(v) for v in r):
        raise ValueError(f"range must be two whole numbers, not {rng!r}")
    a, b = int(r[0]), int(r[1])
    if a < 1 or b < a or b > _num_pages(path):
        raise ValueError("Range is out of bounds")
    return [a - 1, b]


def read_tiff(
    filepath, range: Sequence[int] | None = None, *, num_threads: int | None = None
) -> np.ndarray:
    """Read a tiff file as ``(y, x, z)`` Fortran order (MATLAB ``readtiff``).

    ``range``: ``[first, last]`` slices, 1-based and inclusive, as in MATLAB. One slice
    comes back 2D ``(y, x)``; color images ``(y, x, c, z)``.
    """
    import cpptiff

    path = path_str(filepath)
    if not os.path.isfile(path):
        raise FileNotFoundError(f'File "{path}" cannot be opened')
    z_range = None if range is None or len(range) == 0 else _check_range(path, range)
    n = threads(num_threads)
    if _reader() == "order":
        return matlab_squeeze(cpptiff.read_tiff(path, z_range, num_threads=n, order="F"))

    r = cpptiff.read_tiff(path, z_range, num_threads=n)
    color = r.ndim == 4 or (r.ndim == 3 and r.shape[2] in (3, 4) and _tiff_header(path).samples > 1)
    if color:
        # (z, y, x, c) or one image (y, x, c) -> MATLAB's (y, x, c, z)
        r = r if r.ndim == 4 else r[np.newaxis]
        return matlab_squeeze(np.asfortranarray(r.transpose(1, 2, 3, 0)))
    nz, ny, nx = (1, *r.shape) if r.ndim == 2 else r.shape
    kernel = _inplace_kernel()
    if kernel is not None and r.flags.c_contiguous and r.flags.writeable:
        # the pages turned around in cpp-tiff's buffer: no second volume
        flat = r.reshape(-1)
        kernel(flat, nx, ny, nz, n)
        return matlab_squeeze(flat.reshape((ny, nx, nz), order="F"))
    # (z, y, x) C order is (x, y, z) F order in memory
    src = r.reshape((nz, ny, nx)).transpose(2, 1, 0)
    dst = np.empty((ny, nx, nz), dtype=r.dtype, order="F")
    _transpose_numpy(src, dst, n)
    return matlab_squeeze(dst)


def read_tiff_parallel(
    filepath,
    range: Sequence[int] | None = None,
    info=None,
    show_waitbar: bool = False,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """MATLAB ``readtiff_parallel``: the slices listed in ``range`` (1-based, any order),
    stacked along z; all slices without it. ``info`` and ``show_waitbar`` are accepted and
    not used."""
    if range is None or len(range) == 0:
        return read_tiff(filepath, num_threads=num_threads)
    pages = [int(p) for p in np.asarray(range).ravel()]
    # contiguous ascending runs, each one range read
    runs = [[pages[0], pages[0]]]
    for p in pages[1:]:
        if p == runs[-1][1] + 1:
            runs[-1][1] = p
        else:
            runs.append([p, p])
    parts = [read_tiff(filepath, r, num_threads=num_threads) for r in runs]
    parts = [p if p.ndim >= 3 else p[:, :, np.newaxis] for p in parts]
    out = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=2)
    return matlab_squeeze(np.asfortranarray(out))


# --- writetiff -----------------------------------------------------------------------------


def _check_image(img) -> np.ndarray:
    a = np.asarray(img)
    if a.dtype.kind == "b" or a.dtype.newbyteorder("=") not in _TIFF_TYPES:
        raise TypeError(f"writetiff writes numeric data of MATLAB's classes, not {a.dtype}")
    if a.ndim not in (2, 3):
        raise ValueError(f"writetiff writes 2D or 3D data, not {a.ndim}D")
    if 0 in a.shape:
        raise ValueError(f"writetiff cannot write empty data of size {a.shape}")
    if not a.dtype.isnative:
        a = a.astype(a.dtype.newbyteorder("="))
    return a


def group_writable(path: str) -> None:
    """MATLAB's ``fileattrib(path, '+w', 'g')`` (not on Windows), with its warning."""
    if os.name == "nt":
        return
    try:
        os.chmod(path, os.stat(path).st_mode | stat.S_IWGRP)
    except OSError:
        warnings.warn("Unable to change file attribe for group write!", stacklevel=3)


def write_tiff(
    img,
    filepath,
    *,
    compression: str = "lzw",
    mode: str = "parallel",
    group_write: bool = True,
    num_threads: int | None = None,
) -> None:
    """Write a ``(y, x)`` or ``(y, x, z)`` array as a tiff file (MATLAB ``writetiff``).

    ``compression``: ``'lzw'`` (default), ``'zstd'`` or ``'none'``. Every mode writes
    cpp-tiff's BigTIFF; MATLAB's ``'libtiff'`` and ``'imwrite'`` modes (and ``'none'`` in
    parallel mode) write classic TIFF through ``saveastiff`` instead, with LZW for
    ``'zstd'`` in those modes, as here.
    """
    import cpptiff

    if compression not in _COMPRESSIONS:
        raise ValueError(f"Compression must be one of {_COMPRESSIONS}, not {compression!r}")
    if mode not in _MODES:
        raise ValueError(f"Mode must be one of {_MODES}, not {mode!r}")
    # saveastiff (the libtiff modes) cannot write zstd: it writes lzw
    comp = compression if mode == "parallel" or compression == "none" else "lzw"
    a = _check_image(img)
    path = path_str(filepath)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    n = threads(num_threads)
    how = _writer()
    if how == "order":
        cpptiff.write_tiff(path, np.asfortranarray(a), compression=comp, num_threads=n, order="F")
    elif how == "direct":
        from cpptiff.cpptiff import pybind11_write_tiff

        if a.ndim == 2 and a.flags.c_contiguous and not a.flags.f_contiguous:
            # a C-ordered image is (x, y) in Fortran order already
            pybind11_write_tiff(path, a.T, False, comp, 0, n)
        else:
            a = np.asfortranarray(a)
            view = a.transpose(1, 0, 2) if a.ndim == 3 else a.T
            pybind11_write_tiff(path, view, True, comp, 0, n)
    else:
        cpptiff.write_tiff(
            path, a.transpose(2, 0, 1) if a.ndim == 3 else a, compression=comp, num_threads=n
        )
    if group_write:
        group_writable(path)
