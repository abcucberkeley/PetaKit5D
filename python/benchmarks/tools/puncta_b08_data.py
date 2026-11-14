"""The data of the puncta removal benchmark (B08): an expanded-sample-like uint16 volume
(z, y, x) with stained tubes and blocks on a noisy background and bright puncta, written as
paired zarr inputs (<folder>/matlab/expanded.zarr in MATLAB's layout, <folder>/python/
expanded.zarr native).

usage: python puncta_b08_data.py FOLDER [--shape NZ NY NX]
"""

from __future__ import annotations

import argparse

import numpy as np

from petakit5d.testing import synth
from petakit5d.testing.synth import CounterRNG


def volume(shape: tuple[int, int, int]) -> np.ndarray:
    rng = CounterRNG(0, "b08/bench")
    nz, ny, nx = shape
    v = 100.0 + 20.0 * rng.normal(shape)
    z = np.arange(nz, dtype=np.float64)[:, None, None]
    y = np.arange(ny, dtype=np.float64)[None, :, None]
    x = np.arange(nx, dtype=np.float64)[None, None, :]
    # tubes along x and z, and blocks: the stained structure
    for k in range(6):
        cy, cz = (k + 1) * ny / 7, (k % 3 + 1) * nz / 4
        r2 = (y - cy) ** 2 + (z - cz) ** 2
        v = v + 4000.0 * (r2 < 36.0)
        cy2, cx2 = (k + 0.5) * ny / 6, (5 - k + 0.5) * nx / 6
        v = v + 3000.0 * (((y - cy2) ** 2 + (x - cx2) ** 2) < 25.0)
    for k in range(4):
        z0, y0, x0 = int(nz * (k + 1) / 6), int(ny * (k + 1) / 5), int(nx * (4 - k) / 5)
        v[z0 : z0 + 20, y0 : y0 + 40, x0 : x0 + 30] += 4500.0
    # bright puncta
    n = max(20, int(nz * ny * nx / 150000))
    u = rng.spawn("p").uniform((n, 3))
    pos = u * (np.array(shape, dtype=np.float64) - 1)
    v = v + synth.beads(shape, pos, 1.5, 20000.0, 0.0).astype(np.float64)
    return np.clip(np.floor(v + 0.5), 0, 65535).astype(np.uint16)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--shape", type=int, nargs=3, default=(200, 384, 384))
    a = ap.parse_args()
    synth.write_paired_inputs(a.folder, "expanded", volume(tuple(a.shape)), "zarr", (64, 64, 64))


if __name__ == "__main__":
    main()
