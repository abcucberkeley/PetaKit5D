"""The data of the PSF analysis benchmark (B06): the demo PSF folder of the demo cell image
dataset (PETAKIT5D_DEMO_DATA: its PSFs and Richards & Wolf PSFs), and a calibration image of
beads for the detection (512 x 256 x 201, 24 beads, uint16), in <folder>/matlab/ and
<folder>/python/ (the same files).

usage: python psf_analysis_b06_data.py FOLDER
"""

from __future__ import annotations

import argparse
import os
import shutil

import numpy as np


def beads_image(path: str) -> None:
    from petakit5d.core.dtypes import uint16
    from petakit5d.testing import synth

    rng = synth.CounterRNG(0, "b06/bench")
    sz = (201, 512, 256)
    n = 24
    u = rng.uniform((n, 3))
    pos = 20 + u * (np.array(sz, dtype=np.float64) - 40)
    im = synth.beads(sz, pos, (3.0, 1.6), 2000.0, 100.0)
    im = im + np.float32(30) * rng.uniform(sz, np.float32)
    synth.write_tiff(path, uint16(im))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    a = ap.parse_args()
    demo = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "PetaKit5D_demo_cell_image_dataset")
    first = os.path.join(a.folder, "matlab")
    psf = os.path.join(first, "PSF")
    os.makedirs(os.path.join(psf, "RW_PSFs"), exist_ok=True)
    for n in ("488_2_c.tif", "560_c.tif"):
        shutil.copy(os.path.join(demo, "PSF", n), os.path.join(psf, n))
    for n in os.listdir(os.path.join(demo, "PSF", "RW_PSFs")):
        shutil.copy(os.path.join(demo, "PSF", "RW_PSFs", n), os.path.join(psf, "RW_PSFs", n))
    os.makedirs(os.path.join(first, "beads"), exist_ok=True)
    beads_image(os.path.join(first, "beads", "beads.tif"))
    shutil.copytree(first, os.path.join(a.folder, "python"), dirs_exist_ok=True)


if __name__ == "__main__":
    main()
