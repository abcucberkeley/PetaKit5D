"""The data of the N5 benchmark (B09): the first CamB ch0 tile of the demo cell image dataset
(PETAKIT5D_DEMO_DATA) as a zarr (chunks 256^3) in MATLAB's layout in <folder>/matlab/t.zarr and
as the native array in <folder>/python/t.zarr, and the N5 dataset of it each side writes as
s.n5 (MATLAB's layout, dimensions [nz, nx, ny], for MATLAB; standard N5 for Python; blosc zstd
level 1, chunks 256^3 either way).

usage: python n5_b09_data.py FOLDER
"""

from __future__ import annotations

import argparse
import os
import shutil


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    a = ap.parse_args()
    import cppzarr

    from petakit5d.core.arrays import to_matlab
    from petakit5d.io import read_tiff
    from petakit5d.io.n5 import create_n5, write_n5, zarr_to_n5

    demo = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "PetaKit5D_demo_cell_image_dataset")
    tile = sorted(n for n in os.listdir(demo) if "CamB_ch0" in n and n.endswith(".tif"))[0]
    vol = read_tiff(os.path.join(demo, tile))
    ml, py = os.path.join(a.folder, "matlab"), os.path.join(a.folder, "python")
    os.makedirs(ml, exist_ok=True)
    os.makedirs(py, exist_ok=True)
    cppzarr.write_zarr(os.path.join(ml, "t.zarr"), to_matlab(vol), chunks=[256] * 3, order="F")
    cppzarr.write_zarr(os.path.join(py, "t.zarr"), vol, chunks=[256] * 3, order="C")
    zarr_to_n5(os.path.join(py, "t.zarr"), parse_cluster=False)
    os.replace(os.path.join(py, "n5", "t.n5"), os.path.join(py, "s.n5"))
    shutil.rmtree(os.path.join(py, "n5"))
    arr = create_n5(os.path.join(ml, "s.n5"), vol.shape, (256, 256, 256), vol.dtype, axes="yxz")
    write_n5(vol, arr)


if __name__ == "__main__":
    main()
