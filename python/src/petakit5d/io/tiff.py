"""Tiff reading and writing through cpp-tiff (ports readtiff, readtiff_parallel, writetiff).

Volumes are ``(z, y, x)`` C order, cpp-tiff's own arrays, used as they are: no copy on reading
or writing. One page comes back as an image ``(y, x)``. Ranges of slices are ``(start, stop)``,
0-based, end-exclusive (MATLAB's ``[first last]`` is ``(first - 1, last)``). The files are the
mex's: BigTIFF, ``RowsPerStrip`` 512, LZW (the default), zstd or none; MATLAB's ``readtiff``
reads them as its ``(y, x, z)`` arrays of the same data (decision 0004).
"""

from __future__ import annotations

import os
import stat
import warnings
from collections.abc import Sequence

import numpy as np

from petakit5d.io.imagesize import _imagej_images, _tiff_header, path_str, threads

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
    r = list(rng)
    if len(r) != 2:
        raise ValueError(f"range is (start, stop), got {rng!r}")
    n = None
    if r[1] is None:
        n = _num_pages(path)
        r[1] = n
    if any(not np.isfinite(float(v)) or float(v) != int(v) for v in r):
        raise ValueError(f"range must be two whole numbers, not {rng!r}")
    a, b = int(r[0]), int(r[1])
    # (cpp-tiff 1.7.3 reports a range past the end as a damaged file, T7)
    if a < 0 or b <= a or b > (_num_pages(path) if n is None else n):
        raise ValueError("Range is out of bounds")
    return [a, b]


def read_tiff(
    filepath, range: Sequence[int] | None = None, *, num_threads: int | None = None
) -> np.ndarray:
    """Read a tiff file as a ``(z, y, x)`` C-ordered array (MATLAB ``readtiff``).

    ``range``: the slices ``(start, stop)``, 0-based, end-exclusive (``stop`` ``None``: to the
    end). One slice comes back as an image ``(y, x)``; color images as ``(z, y, x, c)``.
    """
    import cpptiff

    path = path_str(filepath)
    if not os.path.isfile(path):
        raise FileNotFoundError(f'File "{path}" cannot be opened')
    z_range = None if range is None or len(range) == 0 else _check_range(path, range)
    return cpptiff.read_tiff(path, z_range, num_threads=threads(num_threads))


def read_tiff_parallel(
    filepath,
    range: Sequence[int] | None = None,
    info=None,
    show_waitbar: bool = False,
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """MATLAB ``readtiff_parallel``: the slices listed in ``range`` (0-based indices, any order)
    stacked along z; all slices without it. ``info`` and ``show_waitbar`` are accepted and not
    used."""
    if range is None or len(range) == 0:
        return read_tiff(filepath, num_threads=num_threads)
    pages = [int(p) for p in np.asarray(range).ravel()]
    # contiguous ascending runs, each one range read
    runs = [[pages[0], pages[0] + 1]]
    for p in pages[1:]:
        if p == runs[-1][1]:
            runs[-1][1] = p + 1
        else:
            runs.append([p, p + 1])
    parts = [read_tiff(filepath, r, num_threads=num_threads) for r in runs]
    if len(parts) == 1:
        return parts[0]
    return np.concatenate([p if p.ndim == 3 else p[np.newaxis] for p in parts], axis=0)


# --- writetiff -----------------------------------------------------------------------------


def _check_image(img) -> np.ndarray:
    a = np.asarray(img)
    if a.dtype.kind == "b" or a.dtype.newbyteorder("=") not in _TIFF_TYPES:
        raise TypeError(f"writetiff writes numeric data of MATLAB's classes, not {a.dtype}")
    if a.ndim not in (2, 3):
        raise ValueError(f"writetiff writes an image (y, x) or a volume (z, y, x), not {a.ndim}D")
    if 0 in a.shape:
        # (cpp-tiff 1.7.3 writes from a 0-byte buffer, T4)
        raise ValueError(f"writetiff cannot write empty data of size {a.shape}")
    if a.flags.f_contiguous and not a.flags.c_contiguous:
        raise ValueError(
            "writetiff takes C-ordered (z, y, x) data; this array is Fortran-ordered "
            "(MATLAB's layout? use core.arrays.from_matlab)"
        )
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
    """Write an image ``(y, x)`` or a volume ``(z, y, x)`` as a tiff file, one page per z
    (MATLAB ``writetiff``).

    ``compression``: ``'lzw'`` (default), ``'zstd'`` or ``'none'``. Every mode writes
    cpp-tiff's BigTIFF; MATLAB's ``'libtiff'`` and ``'imwrite'`` modes (and ``'none'`` in
    parallel mode) write classic TIFF through ``saveastiff`` instead, with LZW for
    ``'zstd'`` in those modes, as here. A C-ordered array is written without a copy.
    """
    import cpptiff

    # only cpp-tiff's exact names: it writes LZW for any other (T5)
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
    # (strided views are copied to C order by cpp-tiff's write_tiff)
    cpptiff.write_tiff(path, a, compression=comp, num_threads=threads(num_threads))
    if group_write:
        group_writable(path)
