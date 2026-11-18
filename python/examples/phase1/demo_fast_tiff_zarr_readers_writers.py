"""The fast tiff and zarr readers and writers (MATLAB:
``demos/demo_fast_tiff_zarr_readers_writers.m``): read a tile, replicate it to 10,000 frames,
and write and read it as tiff (libtiff mode and cpp-tiff) and as zarr (cpp-zarr).

Left out: MATLAB's own tiff reader (``readtiff_matlab``, the conventional reader the demo
compares with) and the zarr reader and writer of MATLAB's blockedImage, which have no Python
counterpart here. Results: ``<data>/PetaKit5D_demo_cell_image_dataset/replicated/``.
"""

from __future__ import annotations

import os

import _common

DEMO = "demo_fast_tiff_zarr_readers_writers"


def main() -> None:
    from petakit5d.io import create_zarr, read_tiff, read_zarr, write_tiff, write_zarr

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()
    fsn = _common.TILE_A
    fn = data_path + fsn + ".tif"

    steps.skip("read one tile with the conventional reader", "readtiff_matlab is MATLAB's own")
    with steps("read one tile with cpp-tiff"):
        im = read_tiff(fn)

    nframe = 10000
    with steps(f"replicate the tile to {nframe} frames"):
        # MATLAB's repmat(im, 1, 1, ceil(nframe / 500)), cropped to nframe slices
        im_rep = _common.replicate_frames(im, nframe)
        print(im_rep.shape)
        out_path = data_path + "replicated/"
        os.makedirs(out_path, exist_ok=True)

    with steps(f"write {nframe} frames in libtiff mode (the conventional writer)"):
        write_tiff(
            im_rep, f"{out_path}{fsn}_frame_number_{nframe}_conventional.tif", mode="libtiff"
        )

    fnout = f"{out_path}{fsn}_frame_number_{nframe}.tif"
    with steps(f"write {nframe} frames with cpp-tiff"):
        write_tiff(im_rep, fnout)

    # MATLAB clears im_rep and reads it again, into memory its memory manager kept; here the
    # reads go into the same array (out=), which saves the page faults of a new 18 GB array
    steps.skip(f"read {nframe} frames with the conventional reader", "readtiff_matlab")
    with steps(f"read {nframe} frames with cpp-tiff"):
        im_rep = read_tiff(fnout, out=im_rep)

    # zarr: chunks of 256^3 (z, y, x), zstd compression (the defaults)
    zarr_fnout = f"{out_path}{fsn}_frame_number_{nframe}.zarr"
    create_zarr(zarr_fnout, data_size=im_rep.shape, block_size=(256, 256, 256))
    with steps(f"write {nframe} frames with cpp-zarr"):
        write_zarr(im_rep, zarr_fnout)

    with steps(f"read {nframe} frames with cpp-zarr"):
        im_rep = read_zarr(zarr_fnout, out=im_rep)
    del im_rep

    steps.skip("zarr with MATLAB's blockedImage and zarr-python", "MATLAB only")


if __name__ == "__main__":
    _common.run_demo(main)
