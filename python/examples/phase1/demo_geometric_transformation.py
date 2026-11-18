"""Deskew and rotation (MATLAB: ``demos/demo_geometric_transformation.m``): the combined
deskew/rotation of the whole dataset, then of one tile, and of the tile replicated to 2,000 and
5,000 frames.

Left out: the conventional separate deskew and rotation (``XR_deskewRotateFrame_old``) that the
demo times for comparison; it is legacy code and not ported. Results: ``DSR/`` (the dataset, as
zarr), ``DSR_combined/`` and ``replicated/DSR_combined/`` in the dataset folder.
"""

from __future__ import annotations

import os
import shutil

import _common

DEMO = "demo_geometric_transformation"


def move_into(src: str, folder: str) -> None:
    """MATLAB's movefile(src, folder)"""
    os.replace(src, os.path.join(folder, os.path.basename(src)))


def main() -> None:
    from petakit5d.deskew_rotate import deskew_rotate_dataset, deskew_rotate_file
    from petakit5d.io import read_tiff, write_tiff

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()

    with steps("deskew/rotate the dataset (combined, zarr output)"):
        deskew_rotate_dataset(
            [data_path],
            deskew=True,
            rotate=True,
            dsr_combined=True,
            xy_pixel_size=0.108,
            dz=0.3,
            skew_angle=32.45,
            objective_scan=False,
            reverse=True,
            # MATLAB's 'yxz' axis orders are the arrays' own order, 'zyx' here
            input_axis_order="zyx",
            output_axis_order="zyx",
            channel_patterns=["CamA", "CamB"],
            # flat-field correction (off here) and its images
            ff_correction=False,
            lower_limit=0.4,
            const_offset=1,
            ff_image_paths=[
                data_path + "FF/averaged/ff_CamA_ch0_CAM1_stack0000_488nm_0000000msec_"
                "0096911533msecAbs_000x_000y_000z_0017t.tif",
                data_path + "FF/averaged/ff_CamB_ch0_CAM1_stack0000_488nm_0000000msec_"
                "0096911533msecAbs_000x_000y_000z_0017t.tif",
            ],
            background_paths=[
                _common.background_path(data_path, "CamA"),
                _common.background_path(data_path, "CamB"),
            ],
            large_file=False,
            zarr_file=False,
            save_zarr=True,
            block_size=(256, 256, 256),
            save_16bit=True,
            parse_cluster=False,
            master_compute=True,
            config_file="",
            mcc_mode=False,
        )

    fsn = _common.TILE_A
    fn = data_path + fsn + ".tif"
    with steps("read one tile"):
        im = read_tiff(fn)

    xy_pixel_size = 0.108
    dz = 0.3
    for d in ("DS", "DSR"):
        if os.path.exists(f"{data_path}{d}/{fsn}.tif"):
            os.remove(f"{data_path}{d}/{fsn}.tif")
    steps.skip("conventional deskew/rotation of 500 frames", "XR_deskewRotateFrame_old")

    with steps("combined deskew/rotation of 500 frames"):
        deskew_rotate_file(fn, xy_pixel_size, dz, rotate=True, reverse=True, save_16bit=True)
    if os.path.isdir(data_path + "DSR_combined"):
        shutil.rmtree(data_path + "DSR_combined")
    os.makedirs(data_path + "DSR_combined")
    move_into(f"{data_path}DSR/{fsn}.tif", data_path + "DSR_combined")

    out_path = data_path + "replicated/"
    for nframe, tag in ((2000, "2k"), (5000, "5k")):
        with steps(f"replicate the tile to {nframe} frames"):
            # MATLAB's repmat(im, 1, 1, ceil(nframe / 500)), cropped to nframe slices
            im_rep = _common.replicate_frames(im, nframe)
            print(im_rep.shape)
            os.makedirs(out_path, exist_ok=True)
            fnrep = f"{out_path}{fsn}_{tag}.tif"
            write_tiff(im_rep, fnrep)
            del im_rep

        if nframe == 2000:
            steps.skip(f"conventional deskew/rotation of {nframe} frames", "legacy, not ported")
        with steps(f"combined deskew/rotation of {nframe} frames"):
            deskew_rotate_file(fnrep, xy_pixel_size, dz, rotate=True, reverse=True, save_16bit=True)
        os.makedirs(out_path + "DSR_combined/", exist_ok=True)
        move_into(f"{out_path}DSR/{fsn}_{tag}.tif", out_path + "DSR_combined")


if __name__ == "__main__":
    _common.run_demo(main)
