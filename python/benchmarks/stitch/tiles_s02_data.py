"""Tile preparation benchmark (S02), the data: 16 DSR-like uint16 tiles (z, y, x) = (201, 512,
640) with an empty wedge on one side, as tiffs (shared) and as zarr (each side's layout), a
flat-field and a background image, and the function file of the flat field and the erosion
(``fun.mat``, with the absolute paths of this folder).

usage: python tiles_s02_data.py FOLDER
"""

from __future__ import annotations

import os
import sys

import numpy as np

SHAPE = (201, 512, 640)
N = 16


def main(folder: str) -> None:
    from petakit5d.io import write_mat
    from petakit5d.testing.synth import CounterRNG, write_paired_inputs, write_tiff

    folder = os.path.abspath(folder).replace("\\", "/")
    os.makedirs(f"{folder}/tif", exist_ok=True)
    z = np.arange(SHAPE[0])[:, None]
    x = np.arange(SHAPE[2])[None, :]
    wedge = np.broadcast_to((x < z // 2 + 3)[:, None, :], SHAPE)
    for k in range(N):
        v = CounterRNG(7, f"s02bench/{k}").integers(100, 4000, SHAPE).astype(np.uint16)
        v[wedge] = 0  # the empty wedge of a DSR tile
        name = f"tile_{k:02d}"
        write_tiff(f"{folder}/tif/{name}.tif", v)
        write_paired_inputs(f"{folder}/zarr", name, v, "zarr", (125, 250, 250))
    rng = CounterRNG(8, "s02bench/ff")
    write_tiff(f"{folder}/ff.tif", rng.integers(1000, 2000, SHAPE[1:]).astype(np.uint16))
    write_tiff(f"{folder}/bg.tif", np.full(SHAPE[1:], 100, np.uint16))
    fun = (
        f"@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'{folder}/ff.tif',"
        f"'{folder}/bg.tif',1,0.4,false),1)"
    )
    write_mat(f"{folder}/fun.mat", {"usrFun": fun})


if __name__ == "__main__":
    main(sys.argv[1])
