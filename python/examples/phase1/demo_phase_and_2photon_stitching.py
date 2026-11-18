"""Stitching of phase and 2-photon images (MATLAB:
``demos/demo_phase_and_2photon_stitching.m``), on the second demo dataset: the phase images'
flat-field correction (BaSiC), their image list from the tile positions and their 2D stitching;
then the 2-photon tiles' image list and their 3D stitching with a flat-field correction and an
edge erosion applied to each tile.

The phase flat-field correction needs BaSiCPy (the ``basic`` extra) and is left out without it;
the phase stitching reads the dataset's own corrected images (``Phase/ff_corrected/``), as in
MATLAB. The registration of the tiles is not ported yet: see ``_common.registration``. Results:
``Phase/ff_corrected_1/``, ``Phase/ff_corrected/matlab_stitch/`` and ``2P/matlab_stitch/``.
"""

from __future__ import annotations

import datetime
import os

import _common

DEMO = "demo_phase_and_2photon_stitching"


def main() -> None:
    from petakit5d.core.paths import matlab_dir
    from petakit5d.core.text import sprintf
    from petakit5d.flatfield import phase_image_flatfield_correction
    from petakit5d.flatfield.basic import basic_available
    from petakit5d.io import write_mat
    from petakit5d.stitch import generate_image_list, stitch_dataset

    steps = _common.Steps(DEMO)
    data_path = _common.data_path(_common.OTHER)

    # ---- phase images
    curr_data_path = data_path + "/Phase/"
    if basic_available():
        with steps("Phase: flat-field correction (BaSiC) per time point"):
            fns = [curr_data_path + e.name for e in matlab_dir(curr_data_path + "*.tif")]
            out_path = sprintf("%s/ff_corrected_1/", curr_data_path)
            os.makedirs(out_path, exist_ok=True)
            for i in range(2):
                fns_i = [f for f in fns if sprintf("Iter_%04d", i) in f]
                phase_image_flatfield_correction(fns_i, out_path, 6000)
    else:
        steps.skip("Phase: flat-field correction", "BaSiCPy is not installed (the basic extra)")

    with steps("Phase: image list from the tile positions"):
        generate_image_list(
            curr_data_path,
            "tile_position",
            xy_pixel_size=0.108,
            channel_patterns=["Scan"],
            tile_patterns=["0000t", "ch0", "000x", "000y", "000z"],
            # (z, y, x); MATLAB's [100, 100, 1] is [y x z]
            overlap_size=(1, 100, 100),
            overlap_size_type="pixel",
            objective_scan=False,
            io_scan=True,
        )

    phase_dir = data_path + "Phase/ff_corrected/"
    xcorr_shift = _common.registration(phase_dir, "matlab_stitch", _common.OTHER)
    with steps("Phase: 2D stitching"):
        stitch_dataset(
            phase_dir,
            data_path + "Phase/ImageList_all_timepoints.csv",
            ds=False,
            dsr=False,
            processed_dir_str="",
            reverse=True,
            objective_scan=False,
            io_scan=True,
            axis_order="x,y,z",
            xcorr_shift=xcorr_shift,
            xcorr_mode="primary",
            xy_max_offset=100,
            z_max_offset=10,
            xcorr_downsample=(1, 1, 1),
            channel_patterns=["CamB_ch0"],
            resample_type="isotropic",
            xy_pixel_size=0.108,
            blend_method="feather",
            # (z, y, x); MATLAB's [512, 512, 1] and [2048, 2048, 1]
            block_size=(1, 512, 512),
            batch_size=(1, 2048, 2048),
            result_dir_name="matlab_stitch",
            parse_setting_file=False,
            process_fun_path="",
            tile_offset=0,
            edge_artifacts=0,
            save_16bit=True,
            parse_cluster=False,
        )

    # ---- 2-photon images
    curr_data_path = data_path + "/2P/"
    with steps("2P: image list from the tile positions"):
        generate_image_list(
            curr_data_path,
            "tile_position",
            xy_pixel_size=0.108,
            dz=0.5,
            channel_patterns=["Tile"],
            tile_patterns=["0000t", "ch0", "000x", "000y", "000z"],
            # (z, y, x) in um; MATLAB's [10, 10, 5]
            overlap_size=(5, 10, 10),
            overlap_size_type="um",
            objective_scan=False,
            io_scan=True,
        )

    with steps("2P: the tile function file (flat field, then edge erosion)"):
        tmp_path = curr_data_path + "/tmp/processFunc/"
        os.makedirs(tmp_path, exist_ok=True)
        ls_image_path = curr_data_path + "/flatfield/flatfield.tif"
        background_path = curr_data_path + "/flatfield/background.tif"
        tile_offset, lower_limit, ls_background, edge_artifacts = 1, 0.4, False, 1
        usr_fun = sprintf(
            "@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'%s','%s',%d,%d,%s),%0.10d)",
            ls_image_path,
            background_path,
            tile_offset,
            lower_limit,
            "true" if ls_background else "false",
            edge_artifacts,
        )
        dt = f"{datetime.datetime.now():%Y%m%d%H%M%S}"
        fn = sprintf("%s/processFunction_c%04d_%s.mat", tmp_path, 1, dt)
        write_mat(fn, {"usrFun": usr_fun}, version="7.3")

    xcorr_shift = _common.registration(curr_data_path, "matlab_stitch", _common.OTHER)
    with steps("2P: 3D stitching"):
        stitch_dataset(
            curr_data_path,
            curr_data_path + "ImageList_from_tile_positions.csv",
            ds=False,
            dsr=False,
            processed_dir_str="",
            reverse=True,
            objective_scan=False,
            io_scan=True,
            axis_order="-x,y,-z",
            channel_patterns=["ch0"],
            primary_ch="ch0",
            xcorr_shift=xcorr_shift,
            xcorr_downsample=(1, 1, 1),
            shift_method="grid",
            xy_max_offset=200,
            z_max_offset=10,
            xy_pixel_size=0.15,
            dz=0.5,
            blend_method="feather",
            result_dir_name="matlab_stitch",
            process_fun_path=[fn],
            save_16bit=False,
            block_size=(256, 256, 256),
            batch_size=(512, 512, 512),
            parse_cluster=False,
            master_compute=True,
            config_file="",
            mcc_mode=False,
        )


if __name__ == "__main__":
    _common.run_demo(main)
