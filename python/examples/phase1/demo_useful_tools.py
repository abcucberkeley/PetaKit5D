"""The data tools (MATLAB: ``demos/demo_useful_tools.m``, its Phase 1 tools): resampling,
cropping, MIPs, the decon with intermediate results, tiff to zarr and back, image lists from a
tile list, channel unmixing (tiff and zarr) and the resaving of zarr files.

The demo's Phase 2 steps (the FFT spectra, the PSF and FSC analyses, the chromatic shift
correction and the Bio-Formats converter) are ``examples/phase2/demo_useful_tools.py``, to run
after this one; the Imaris converter is not ported yet. The steps that MATLAB runs with
``parseCluster=true`` run on this machine when there is no Slurm, as in MATLAB.
Results: in the dataset folder (``resampled/``, ``Cropped/``, ``MIPs/``, ``tiff_to_zarr/``,
``Unmixed/``, ...).
"""

from __future__ import annotations

import _common

DEMO = "demo_useful_tools"


def main() -> None:
    from petakit5d.decon import decon_dataset
    from petakit5d.deskew_rotate import deskew_rotate_dataset
    from petakit5d.io.convert import tiff_to_zarr_dataset, zarr_to_tiff_dataset
    from petakit5d.io.resave import resave_zarr_dataset
    from petakit5d.stitch import generate_image_list
    from petakit5d.tools import (
        crop_dataset,
        mip_dataset,
        resample_dataset,
        unmix_channels_dataset,
        unmix_channels_zarr,
    )

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()
    ff_image_paths = [
        data_path + f"FF/averaged/ff_{c}_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_"
        "000x_000y_000z_0017t.tif"
        for c in ("CamA", "CamB")
    ]
    background_paths = [_common.background_path(data_path, c) for c in ("CamA", "CamB")]
    decon_common = {
        "xy_pixel_size": 0.108,
        "dz": 0.3,
        "dz_psf": 0.3,
        "parse_setting_file": False,
        "psf_fullpaths": [data_path + "PSF/488_2_c.tif"],
        "rl_method": "omw",
        "skewed": True,
        "background": 100,
        "edge_erosion": 0,
        "save_16bit": True,
        "zarr_file": False,
        "save_zarr": False,
        "cpus_per_task": 24,
        "parse_cluster": False,
        "large_file": False,
        "gpu_job": False,
        "config_file": "",
        "gpu_config_file": "",
        "mcc_mode": False,
    }
    # (MATLAB's demo also passes 'fixIter', which the wrapper does not use)

    with steps("decon used by some tools"):
        decon_dataset(
            [data_path + "/"],
            result_dir_name="matlab_decon_omw",
            channel_patterns=["CamB_ch0"],
            wiener_alpha=0.005,
            otf_cum_thresh=0.9,
            decon_iter=2,
            debug=False,
            **decon_common,
        )

    with steps("deskew/rotation used by some tools (separate deskew and rotation)"):
        deskew_rotate_dataset(
            [data_path + "matlab_decon_omw/"],
            ds_dir_name="DS_for_useful_tools/",
            dsr_dir_name="DSR_for_useful_tools/",
            deskew=True,
            rotate=True,
            dsr_combined=False,
            xy_pixel_size=0.108,
            dz=0.3,
            skew_angle=32.45,
            objective_scan=False,
            reverse=True,
            channel_patterns=["CamA", "CamB"],
            ff_correction=False,
            lower_limit=0.4,
            const_offset=1,
            ff_image_paths=ff_image_paths,
            background_paths=background_paths,
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

    with steps("resampling tool"):
        resample_dataset(
            [data_path],
            (2, 2, 2),
            result_dir_name="resampled",
            channel_patterns=["CamA_ch0", "CamB_ch0"],
            interp_method="linear",
            save_16bit=True,
            zarr_file=False,
            large_file=False,
            save_zarr=False,
            block_size=(256, 256, 256),
            batch_size=(512, 512, 512),
            border_size=(5, 5, 5),
            parse_cluster=False,
            master_compute=True,
            cpus_per_task=4,
            mcc_mode=False,
            config_file="",
        )

    with steps("crop tool"):
        crop_dataset(
            [data_path],
            # (z0, y0, x0, z1, y1, x1); MATLAB's [201, 101, 101, 1200, 412, 400]
            (100, 200, 100, 400, 1200, 412),
            result_dir_name="Cropped",
            channel_patterns=["CamA_ch0", "CamB_ch0"],
            zarr_file=False,
            large_file=False,
            save_zarr=False,
            block_size=(256, 256, 256),
            save_16bit=True,
            parse_cluster=False,
            master_compute=True,
            cpus_per_task=4,
            mcc_mode=False,
            config_file="",
        )

    with steps("MIP tool"):
        mip_dataset(
            [data_path],
            result_dir_name="MIPs",
            axis=(1, 1, 1),
            channel_patterns=["CamA_ch0", "CamB_ch0"],
            zarr_file=False,
            large_file=False,
            batch_size=(1024, 1024, 1024),
            save_16bit=True,
            parse_cluster=False,
            master_compute=True,
            cpus_per_task=4,
            mcc_mode=False,
            config_file="",
        )

    steps.skip("FFT analysis", "a Phase 2 step: examples/phase2/demo_useful_tools.py")
    steps.skip("PSF analysis", "a Phase 2 step: examples/phase2/demo_useful_tools.py")

    with steps("decon with intermediate results (for the FSC analysis)"):
        decon_dataset(
            [data_path + "/"],
            result_dir_name="matlab_decon_omw_debug",
            channel_patterns=[_common.TILE_B],
            wiener_alpha=0.05,
            otf_cum_thresh=0.9,
            decon_iter=20,
            save_step=2,
            debug=True,
            **decon_common,
        )
    steps.skip("FSC analysis", "a Phase 2 step: examples/phase2/demo_useful_tools.py")

    with steps("tiff to zarr conversion"):
        tiff_to_zarr_dataset(
            [data_path + "/"],
            result_dir_name="tiff_to_zarr",
            channel_patterns=["CamA_ch0", "CamB_ch0"],
            input_bbox=None,
            block_size=(256, 256, 256),
            read_whole_tiff=True,
            parse_cluster=False,
            master_compute=True,
            cpus_per_task=2,
            mcc_mode=False,
            config_file="",
        )

    with steps("zarr to tiff conversion"):
        zarr_to_tiff_dataset(
            [data_path + "matlab_decon_omw/DSR_for_useful_tools/"],
            result_dir_name="tiffs",
            channel_patterns=["CamA_ch0", "CamB_ch0"],
            parse_cluster=False,
            master_compute=True,
            cpus_per_task=2,
            mcc_mode=False,
            config_file="",
        )

    steps.skip(
        "Imaris converter",
        "XR_imaris_conversion_data_wrapper is not ported yet: it waits for the next release of "
        "Parallel_Imaris_Writer",
    )

    with steps("image list from a tile list"):
        # the 16 tiles' file names: time points, cameras, then y tiles 3 to 0
        stamps = {
            "0000": ("0106060251", "0106082240", "0106104333", "0106126432"),
            "0001": ("0106148492", "0106170570", "0106192817", "0106215043"),
        }
        tile_filenames = [
            f"Scan_Iter_0000_{t}_{cam}_ch0_CAM1_stack0000_488nm_0000000msec_"
            f"{m}msecAbs_000x_{y:03d}y_000z_{t}t.tif"
            for t, ms in stamps.items()
            for cam in ("CamA", "CamB")
            for y, m in zip((3, 2, 1, 0), ms, strict=True)
        ]
        # one row per tile: [t, c, x, y, z] indices; the interval of the tiles in x, y, z (um)
        tile_indices = [[t, c, 0, y, 0] for t in (0, 1) for c in (0, 1) for y in (3, 2, 1, 0)]
        generate_image_list(
            [data_path],
            "tile_list",
            tile_filenames=tile_filenames,
            tile_indices=tile_indices,
            tile_interval=[0, 189.4, 0],
        )

    steps.skip("chromatic shift correction", "a Phase 2 step: examples/phase2/demo_useful_tools.py")

    with steps("channel unmixing"):
        unmix_channels_dataset(
            data_path,
            unmix_factors=[1, -0.1],
            mode="linear",
            unmix_sigmas=None,
            result_dir_name="Unmixed",
            channel_patterns=["CamA", "CamB"],
            # the channel to unmix (0-based; MATLAB's channelInd 1)
            channel_ind=0,
            ff_correction=False,
            lower_limit=0.4,
            ff_image_paths=ff_image_paths,
            background_paths=background_paths,
            const_background=[100, 100],
            const_offset=1,
            zarr_file=False,
            large_file=False,
            save_zarr=False,
            save_16bit=True,
            batch_size=(1024, 1024, 1024),
            block_size=(256, 256, 256),
            border_size=(0, 0, 0),
            parse_cluster=True,
            master_compute=True,
            job_log_dir="../job_logs",
            cpus_per_task=3,
            config_file="",
            mcc_mode=False,
            uuid="",
            debug=False,
        )

    with steps("channel unmixing of zarr files"):
        unmix_channels_zarr(
            [data_path + f"DSR/{_common.TILE_A}.zarr", data_path + f"DSR/{_common.TILE_B}.zarr"],
            [1, -0.1],
            result_dir_name="Unmixed",
            mode="linear",
            unmix_sigmas=None,
            channel_ind=0,
            batch_size=(1024, 1024, 1024),
            block_size=(256, 256, 256),
            border_size=(0, 0, 0),
            parse_cluster=True,
            master_compute=True,
            cpus_per_task=4,
            job_log_dir="../job_logs",
            config_file="",
            mcc_mode=False,
            uuid="",
            debug=False,
        )

    with steps("resaving zarr files"):
        resave_zarr_dataset(
            data_path + "DSR/",
            result_dir_name="resaved_zarr",
            batch_size=(512, 512, 512),
            block_size=(512, 512, 512),
            channel_patterns=[".zarr"],
            input_bbox=None,
            parse_cluster=True,
            large_file=True,
            master_compute=True,
            cpus_per_task=4,
            uuid="",
            mcc_mode=False,
            config_file="",
        )

    steps.skip(
        "Bio-Formats conversion",
        "after the demo's return: examples/phase2/demo_useful_tools.py has it, not run",
    )


if __name__ == "__main__":
    _common.run_demo(main)
