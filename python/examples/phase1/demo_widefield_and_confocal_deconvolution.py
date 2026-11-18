"""Deconvolution of widefield and confocal images (MATLAB:
``demos/demo_widefield_and_confocal_deconvolution.m``): OMW decon of the widefield (20
iterations) and the confocal (3 iterations) demo images, on the second demo dataset.

GPU: as in MATLAB, the widefield decon uses a GPU if this process can use one; the confocal
decon runs on the CPU. Results: ``Widefield/matlab_decon_omw/`` and
``Confocal/matlab_decon_omw/`` in the dataset folder.
"""

from __future__ import annotations

import _common

DEMO = "demo_widefield_and_confocal_deconvolution"


def main() -> None:
    from petakit5d.decon import decon_dataset

    steps = _common.Steps(DEMO)
    data_path = _common.data_path(_common.OTHER)
    common = {
        "result_dir_name": "matlab_decon_omw",
        "xy_pixel_size": 0.157,
        "dz": 0.1,
        "dz_psf": 0.1,
        "parse_setting_file": False,
        "rl_method": "omw",
        "skewed": True,
        "background": 100,
        "edge_erosion": 0,
        "save_16bit": True,
        "zarr_file": False,
        "save_zarr": False,
        "parse_cluster": False,
        "large_file": False,
        "debug": False,
        "cpus_per_task": 24,
        "config_file": "",
        "gpu_config_file": "",
        "mcc_mode": False,
    }
    # (MATLAB's demo also passes 'Reverse' and 'fixIter', which the wrapper does not use)

    with steps("Step A2: OTF mask of the widefield PSF"):
        _common.otf_mask_figure(data_path + "Widefield/WF_PSF.tif", 0.65, False, min_int_thrsh=1e-3)

    with steps("Step A3.2: widefield OMW decon, 20 iterations"):
        decon_dataset(
            [data_path + "Widefield/"],
            channel_patterns=["WF_raw"],
            psf_fullpaths=[data_path + "Widefield/WF_PSF.tif"],
            wiener_alpha=0.02,
            otf_cum_thresh=0.65,
            hann_win_bounds=(0.8, 1.2),
            decon_iter=20,
            gpu_job=_common.can_use_gpu(),
            **common,
        )

    with steps("Step B2: OTF mask of the confocal PSF"):
        _common.otf_mask_figure(
            data_path + "Confocal/Confocal_PSF.tif", 0.875, False, min_int_thrsh=1e-3
        )

    with steps("Step B3.2: confocal OMW decon, 3 iterations"):
        decon_dataset(
            [data_path + "Confocal/"],
            channel_patterns=["Confocal_raw"],
            psf_fullpaths=[data_path + "Confocal/Confocal_PSF.tif"],
            wiener_alpha=0.004,
            otf_cum_thresh=0.875,
            decon_iter=3,
            gpu_job=False,
            **common,
        )


if __name__ == "__main__":
    _common.run_demo(main)
