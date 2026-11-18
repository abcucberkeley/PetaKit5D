"""Large-scale processing (MATLAB: ``demos/demo_large_scale_processing.m``): the skewed-space
stitching (demo_skewed_space_stitching), the large-file OMW decon of the stitched volumes
(in place, block by block) and the large-file deskew/rotation of the results.

GPU: as in MATLAB, the decon uses a GPU if this process can use one. Results:
``matlab_stitch/matlab_decon_omw/`` and its ``DSR/`` in the dataset folder.
"""

from __future__ import annotations

import _common
import demo_skewed_space_stitching

DEMO = "demo_large_scale_processing"


def main() -> None:
    from petakit5d.decon import decon_dataset
    from petakit5d.deskew_rotate import deskew_rotate_dataset

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()

    with steps("Step 1: stitching in skewed space"):
        demo_skewed_space_stitching.main()

    with steps("Step 3.1: OTF masks of the PSFs"):
        _common.otf_mask_figure(data_path + "PSF/560_c.tif", 0.88, True)
        _common.otf_mask_figure(data_path + "PSF/488_2_c.tif", 0.9, True)

    with steps("Step 3.3: large-file OMW decon of the stitched volumes"):
        decon_dataset(
            [data_path + "matlab_stitch/"],
            result_dir_name="matlab_decon_omw",
            xy_pixel_size=0.108,
            dz=0.3,
            dz_psf=0.3,
            parse_setting_file=False,
            channel_patterns=["CamA_ch0", "CamB_ch0"],
            psf_fullpaths=[data_path + "PSF/560_c.tif", data_path + "PSF/488_2_c.tif"],
            rl_method="omw",
            # per channel
            wiener_alpha=[0.05, 0.005],
            otf_cum_thresh=[0.88, 0.9],
            skewed=True,
            background=100,
            decon_iter=2,
            edge_erosion=20,
            save_16bit=True,
            zarr_file=True,
            # (z, y, x); MATLAB's [512, 768, 768] is [y x z]
            batch_size=(768, 512, 768),
            block_size=(256, 256, 256),
            parse_cluster=False,
            large_file=True,
            large_method="inplace",
            gpu_job=_common.can_use_gpu(),
            debug=False,
            cpus_per_task=24,
            config_file="",
            gpu_config_file="",
            mcc_mode=False,
        )
    # MATLAB's demo also passes 'Reverse' and 'fixIter', which the wrapper does not use

    with steps("Step 4: large-file deskew/rotation of the deconvolved volumes"):
        # the demo sets ChannelPatterns but does not pass them: the wrapper's default patterns
        deskew_rotate_dataset(
            [data_path + "matlab_stitch/matlab_decon_omw/"],
            xy_pixel_size=0.108,
            dz=0.3,
            reverse=True,
            large_file=True,
            zarr_file=True,
            save_zarr=True,
            batch_size=(512, 512, 512),
            block_size=(256, 256, 256),
            save_16bit=True,
            parse_cluster=False,
            master_compute=True,
            config_file="",
            mcc_mode=False,
        )


if __name__ == "__main__":
    _common.run_demo(main)
