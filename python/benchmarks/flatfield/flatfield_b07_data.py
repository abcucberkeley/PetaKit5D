"""The data of the flat-field benchmark (B07): for each side (<folder>/matlab, <folder>/python)
a dark dataset of two cameras (uint16 stacks whose noise grows with z, bright slices in a few)
and a copy of the phase tiles of the demo (PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/
Phase/*.tif under the demo data folder).

usage: python flatfield_b07_data.py FOLDER DEMO_DATA [--shape NZ NY NX]
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil

import numpy as np

from petakit5d.testing import synth
from petakit5d.testing.synth import CounterRNG

NAME = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_"
    "000z_0000t.tif"
)
PHASE = "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Phase"


def dark_stack(shape: tuple[int, int, int], cam: str) -> np.ndarray:
    nz, ny, nx = shape
    rng = CounterRNG(3, f"b07/bench/{cam}")
    sig = 2.0 + 4.0 * np.arange(nz, dtype=np.float64) / max(nz - 1, 1)
    y, x = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    pat = 100.0 + (3 * y + 5 * x) % 21
    out = np.empty(shape, np.uint16)
    for z in range(nz):
        v = pat + sig[z] * rng.normal((ny, nx))
        if z % 37 == 5:
            v[ny // 4 : ny // 2, nx // 4 : nx // 2] += 2000.0
        out[z] = np.clip(np.floor(v + 0.5), 0, 65535).astype(np.uint16)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("demo")
    ap.add_argument("--shape", type=int, nargs=3, default=[256, 1024, 1024])
    a = ap.parse_args()
    tiles = sorted(glob.glob(os.path.join(a.demo, PHASE, "*.tif")))
    if not tiles:
        raise SystemExit(f"no phase tiles under {a.demo}")
    stacks = {cam: dark_stack(tuple(a.shape), cam) for cam in ("CamA", "CamB")}
    for side in ("matlab", "python"):
        dark = os.path.join(a.folder, side, "dark")
        phase = os.path.join(a.folder, side, "phase")
        os.makedirs(dark, exist_ok=True)
        os.makedirs(phase, exist_ok=True)
        for cam, st in stacks.items():
            synth.write_tiff(os.path.join(dark, NAME.format(cam)), st)
        for t in tiles:
            shutil.copy(t, phase)
    print(f"dark {a.shape} x 2 cameras, {len(tiles)} phase tiles")


if __name__ == "__main__":
    main()
