"""Stitching in skewed space (MATLAB: ``demos/demo_skewed_space_stitching.m``): the 16 tiles
(4 tiles x 2 cameras x 2 time points) stitched without deskewing them first, with feather
blending, the CamB channel registered and CamA placed by its result.

Results: ``matlab_stitch/`` in the dataset folder. The demo_zarr_stitching and
demo_large_scale_processing examples run this one first.
"""

from __future__ import annotations

import datetime
import os

import _common

DEMO = "demo_skewed_space_stitching"


def process_function_files(data_path: str, channel_patterns: list[str]) -> list[str]:
    """Files of a function applied to each tile before stitching, here the flat-field
    correction (MATLAB's processFFCorrectionFrame). The demo makes them as an example and then
    stitches without them (process_fun_path='')."""
    from petakit5d.core.text import sprintf
    from petakit5d.io import write_mat

    tmp_path = data_path + "/tmp/processFunc/"
    os.makedirs(tmp_path, exist_ok=True)
    ls_image_paths = [
        data_path + "FF/averaged/ff_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_"
        "000x_000y_000z_0017t.tif",
        data_path + "FF/averaged/ff_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_"
        "000x_000y_000z_0017t.tif",
    ]
    background_paths = [
        _common.background_path(data_path, "CamA"),
        _common.background_path(data_path, "CamB"),
    ]
    lower_limit = 0.4
    ls_background = False
    tile_offset = 1
    # MATLAB's datetime format 'yyyyMMddHHmmSS' (SS: hundredths of a second)
    now = datetime.datetime.now()
    dt = f"{now:%Y%m%d%H%M}{now.microsecond // 10000:02d}"
    paths = []
    for i in range(len(channel_patterns)):
        usr_fun = sprintf(
            "@(x)processFFCorrectionFrame(x,'%s','%s',%d,%d,%s)",
            ls_image_paths[i],
            background_paths[i],
            tile_offset,
            lower_limit,
            "true" if ls_background else "false",
        )
        fn = sprintf("%s/processFunction_c%04d_%s.mat", tmp_path, i + 1, dt)
        write_mat(fn, {"usrFun": usr_fun}, version="7.3")
        paths.append(fn)
    return paths


def main() -> None:
    from petakit5d.stitch import stitch_dataset

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()
    image_list_fullpath = data_path + "ImageList_from_encoder.csv"
    channel_patterns = ["CamA_ch0", "CamB_ch0"]

    with steps("process function files (flat-field correction; not used below)"):
        process_function_files(data_path, channel_patterns)

    with steps("stitch the tiles in skewed space"):
        stitch_dataset(
            data_path,
            image_list_fullpath,
            ds=False,
            dsr=False,
            skew_angle=32.45,
            reverse=True,
            xy_pixel_size=0.108,
            dz=0.3,
            parse_setting_file=False,
            # stage axes to tile axes; the tiles' arrays are (z, y, x) (MATLAB's 'y,x,z')
            axis_order="-x,y,z",
            data_order="z,y,x",
            blend_method="feather",
            result_dir_name="matlab_stitch",
            xcorr_shift=True,
            primary_ch="CamB_ch0",
            channel_patterns=channel_patterns,
            zarr_file=False,
            save_16bit=True,
            resample_factor=None,
            resample_type="isotropic",
            xy_max_offset=150,
            z_max_offset=150,
            xcorr_downsample=(2, 2, 2),
            input_bbox=None,
            tile_out_bbox=None,
            batch_size=(512, 512, 512),
            block_size=(256, 256, 256),
            tile_offset=0,
            process_fun_path="",
            parse_cluster=False,
            master_compute=True,
            config_file="",
            mcc_mode=False,
        )


if __name__ == "__main__":
    _common.run_demo(main)
