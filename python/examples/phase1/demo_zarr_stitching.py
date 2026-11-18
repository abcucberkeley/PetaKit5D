"""Stitching (MATLAB: ``demos/demo_zarr_stitching.m``): the skewed-space stitching
(demo_skewed_space_stitching), the deskew/rotation of its result, and the stitching of the tiles
in DSR space (deskewed and rotated tiles as zarr, feather blending, registration).

Left out: the parts after the demo's ``return`` (BigStitcher-Spark and stitching-spark, which
run on a Spark cluster). The registration of the tiles is not ported yet: see
``_common.registration``. Results: ``matlab_stitch/`` with its ``DSR/``, ``DSR/`` (the tiles)
and ``matlab_stitch_dsr/`` in the dataset folder.
"""

from __future__ import annotations

import _common
import demo_skewed_space_stitching

DEMO = "demo_zarr_stitching"


def main() -> None:
    from petakit5d.deskew_rotate import deskew_rotate_dataset
    from petakit5d.stitch import stitch_dataset

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()
    dsr_common = {
        "xy_pixel_size": 0.108,
        "dz": 0.3,
        "reverse": True,
        "channel_patterns": ["CamA", "CamB"],
        "large_file": False,
        "save_zarr": True,
        "block_size": (256, 256, 256),
        "save_16bit": True,
        "parse_cluster": False,
        "master_compute": True,
        "config_file": "",
        "mcc_mode": False,
    }

    with steps("Step 1: stitching in skewed space"):
        demo_skewed_space_stitching.main()

    with steps("Step 2: deskew/rotate the stitched volumes (zarr)"):
        deskew_rotate_dataset([data_path + "matlab_stitch/"], zarr_file=True, **dsr_common)

    with steps("Step 3: deskew/rotate the tiles (zarr output)"):
        deskew_rotate_dataset([data_path], zarr_file=False, **dsr_common)

    xcorr_shift = _common.registration(data_path, "matlab_stitch_dsr")
    with steps("Step 4: stitching in DSR space"):
        stitch_dataset(
            data_path,
            data_path + "ImageList_from_encoder.csv",
            ds=False,
            dsr=True,
            processed_dir_str="DSR",
            zarr_file=True,
            skew_angle=32.45,
            xy_pixel_size=0.108,
            dz=0.3,
            reverse=True,
            parse_setting_file=False,
            axis_order="-x,y,z",
            data_order="z,y,x",
            blend_method="feather",
            result_dir_name="matlab_stitch_dsr",
            xcorr_shift=xcorr_shift,
            primary_ch="CamB_ch0",
            channel_patterns=["CamA_ch0", "CamB_ch0"],
            save_16bit=True,
            resample_factor=None,
            resample_type="isotropic",
            xy_max_offset=150,
            z_max_offset=150,
            xcorr_downsample=(2, 2, 2),
            input_bbox=None,
            tile_out_bbox=None,
            block_size=(256, 256, 256),
            process_fun_path="",
            tile_offset=0,
            parse_cluster=False,
        )

    steps.skip("BigStitcher-Spark and stitching-spark", "after the demo's return; Spark tools")


if __name__ == "__main__":
    _common.run_demo(main)
