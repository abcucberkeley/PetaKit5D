"""The data of the MIP mask benchmark (B10): the MIPs along y, x and z of the first two CamA ch0
tiles of the demo cell image dataset (PETAKIT5D_DEMO_DATA), and its PSF 560_c.tif, in
<folder>/matlab/ and <folder>/python/ (the same tiff files).

usage: python mip_masks_b10_data.py FOLDER [--tiles N]
"""

from __future__ import annotations

import argparse
import os
import shutil


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--tiles", type=int, default=2)
    a = ap.parse_args()
    from petakit5d.tools.mip import save_mip_tiff

    demo = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "PetaKit5D_demo_cell_image_dataset")
    tiles = sorted(n for n in os.listdir(demo) if "CamA_ch0" in n and n.endswith(".tif"))
    first = os.path.join(a.folder, "matlab")
    os.makedirs(os.path.join(first, "mips"), exist_ok=True)
    for n in tiles[: a.tiles]:
        stem = n[:-4]
        save_mip_tiff(
            os.path.join(demo, n),
            os.path.join(first, "mips", f"{stem}_MIP_z.tif"),
            axis=(1, 1, 1),
        )
    os.makedirs(os.path.join(first, "psf"), exist_ok=True)
    shutil.copy(os.path.join(demo, "PSF", "560_c.tif"), os.path.join(first, "psf", "560_c.tif"))
    shutil.copytree(first, os.path.join(a.folder, "python"), dirs_exist_ok=True)


if __name__ == "__main__":
    main()
