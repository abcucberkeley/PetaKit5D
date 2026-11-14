"""Tile functions of the stitching written as MATLAB code (a stand-in until task S02; not
registered: S02 registers erodeVolumeBy2DProjection as planned).

The stitching's frame (``stitch.volume.stitch_volume``) erodes the tiles' edges with a function
it writes as MATLAB code, ``@(x)erodeVolumeBy2DProjection(x,<n>)``, or for skewed space with a
tile offset ``@(x)erodeVolumeBy2DProjection(x+<offset>,<n>)``, and the tile preparation reads
user functions through ``parse_tile_function`` (``io.convert.user_function``). Only these two
forms are read here; any other MATLAB code raises NotImplementedError (S02 ports the parser).

``erode_volume_by_2d_projection`` is erodeVolumeBy2DProjection.m: the edges are eroded by the
erosion of the volume's max projection along y (the xz plane) instead of the whole volume, and
the first and last ``esize`` rows along y are set to 0.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import Any

import numpy as np

__all__ = ["erode_volume_by_2d_projection", "parse_tile_function"]

_ERODE = re.compile(
    r"@\(\s*x\s*\)\s*erodeVolumeBy2DProjection\(\s*x\s*(?:\+\s*([-+]?[0-9.eE+-]+)\s*)?,"
    r"\s*([-+]?[0-9.eE+-]+)\s*\)\s*"
)


def erode_volume_by_2d_projection(
    vol: np.ndarray, esize: Any, *, num_threads: int | None = None
) -> np.ndarray:
    """Erode the edges of a volume (z, y, x) by ``esize`` voxels through its max projection
    along y (erodeVolumeBy2DProjection): the voxels whose (z, x) column is eroded in that
    projection (a square of ``2 esize + 1``, outside the projection as empty) and the first and
    last ``esize`` rows along y become 0. An image (y, x) keeps its values (MATLAB's
    ``vol .* (vol ~= 0)``). Returns a new array of the input's class. ``num_threads`` (None:
    the CPUs this process may use) never changes the result."""
    from petakit5d.matlabcompat.morphology import imerode, strel

    vol = np.asarray(vol)
    e = float(np.asarray(esize).reshape(-1)[0])
    if e == 0:
        return vol
    e = int(e)
    if vol.ndim == 2:
        return np.multiply(vol, (vol != 0).astype(vol.dtype))
    nz, ny, nx = vol.shape
    if ny == 0 or nz == 0:
        return vol.copy()
    chunks = _chunks(nz, vol.nbytes, num_threads)
    # MATLAB's squeeze(max(vol, [], 1) > 0): the projection along y, here (z, x) (MATLAB's max
    # leaves NaN out, as fmax does)
    mip = np.empty((nz, nx), bool)
    _each(chunks, lambda a, b: np.greater(np.fmax.reduce(vol[a:b], axis=1), 0, out=mip[a:b]))
    if mip.all():
        erode = np.zeros_like(mip)
        erode[e : nz - e, e : nx - e] = mip[e : nz - e, e : nx - e]
    else:
        pad = np.zeros((nz + 2, nx + 2), bool)
        pad[1:-1, 1:-1] = mip
        pad = imerode(pad, strel("square", e * 2 + 1))
        erode = pad[1:-1, 1:-1]
    mask = erode.astype(vol.dtype)[:, np.newaxis, :]
    out = np.empty_like(vol)
    lo, hi = min(e, ny), max(ny - e, 0)

    def part(a: int, b: int) -> None:
        # vol .* mask, then the first and last esize rows along y set to 0
        np.multiply(vol[a:b], mask[a:b], out=out[a:b])
        out[a:b, :lo, :] = 0
        out[a:b, hi:, :] = 0

    _each(chunks, part)
    return out


def _chunks(nz: int, nbytes: int, num_threads: int | None) -> list[tuple[int, int]]:
    """z ranges for the threads (one range for small volumes)"""
    from petakit5d.core.threads import available_cpus

    n = available_cpus() if num_threads is None else max(1, int(num_threads))
    n = max(1, min(n, nz, nbytes // (1 << 22)))
    edges = np.linspace(0, nz, n + 1).astype(int)
    return [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:], strict=True) if b > a]


def _each(chunks: list[tuple[int, int]], fn) -> None:
    """fn(z0, z1) for each range, in threads (NumPy releases the GIL in its loops)"""
    if len(chunks) == 1:
        fn(*chunks[0])
        return
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(len(chunks)) as ex:
        for f in [ex.submit(fn, a, b) for a, b in chunks]:
            f.result()


def _plus(vol: np.ndarray, offset: float) -> np.ndarray:
    """MATLAB's ``x + offset`` (double): integer classes round and saturate"""
    from petakit5d.core.dtypes import cast

    if vol.dtype.kind in "iu":
        return np.asarray(cast(vol.astype(np.float64) + offset, str(vol.dtype)))
    return (vol + vol.dtype.type(offset)).astype(vol.dtype)


def parse_tile_function(text: str) -> Callable[[np.ndarray], np.ndarray]:
    """The function of MATLAB code ``@(x)erodeVolumeBy2DProjection(x,<n>)`` or
    ``@(x)erodeVolumeBy2DProjection(x+<offset>,<n>)``, on (z, y, x) volumes."""
    m = _ERODE.fullmatch(text.strip())
    if m is None:
        err = NotImplementedError(
            f"MATLAB code as a tile function ({text!r}): only the stitching's edge erosion "
            "('@(x)erodeVolumeBy2DProjection(x,<n>)', with x+<offset>) is read until the "
            "tile-function parser (task S02) is ported"
        )
        err.identifier = "petakit5d:convert:matlabCodeFunction"  # type: ignore[attr-defined]
        raise err
    esize = float(m.group(2))
    offset = None if m.group(1) is None else float(m.group(1))

    def fun(vol: np.ndarray) -> np.ndarray:
        if offset is not None:
            vol = _plus(np.asarray(vol), offset)
        return erode_volume_by_2d_projection(vol, esize)

    return fun
