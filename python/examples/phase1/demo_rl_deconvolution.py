"""Richardson-Lucy deconvolution (MATLAB: ``demos/demo_RL_deconvolution.m``): OMW and
traditional (Biggs) RL decon of the CamB tiles with the dataset wrapper, the deskew/rotation of
the raw and deconvolved tiles, and the decon functions on one tile in memory.

GPU: as in MATLAB, the wrapper uses a GPU if this process can use one (CuPy and a GPU).

Left out: the Wiener-Butterworth back projector (``BackProjector``, from WBDeconvolution, which
is not ported); the OMW back projector stands in for it in the comparison. Results:
``matlab_decon_omw/``, ``matlab_decon_traditional/`` and their ``DSR/`` folders, and ``DSR/``,
in the dataset folder; the in-memory results with ``PETAKIT5D_DEMO_SAVE_DIR``.
"""

from __future__ import annotations

import _common

DEMO = "demo_RL_deconvolution"


def crop_figure(name: str, images: dict, raw_low: float | None) -> None:
    """the demo's side-by-side crops: clim from the 40th to the 99th percentile, for the raw
    image from ``raw_low`` (None: the 1st percentile)"""
    from petakit5d.matlabcompat.stats import prctile

    clims = {}
    for i, (title, im) in enumerate(images.items()):
        v = im.ravel()
        low = (prctile(v, 1) if raw_low is None else raw_low) if i == 0 else prctile(v, 40)
        clims[title] = (low, prctile(v, 99))
    _common.show_images(DEMO, name, images, clims)


def main() -> None:
    import numpy as np

    from petakit5d.decon import (
        decon_dataset,
        decon_lucy,
        decon_lucy_omw,
        omw_backprojector,
        psf_gen,
    )
    from petakit5d.deskew_rotate import deskew_rotate_dataset, deskew_rotate_volume
    from petakit5d.io import read_tiff
    from petakit5d.kernels.crop import crop_3d
    from petakit5d.matlabcompat.reduce import sum as msum

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()

    with steps("OTF masks of the PSFs"):
        _common.otf_mask_figure(data_path + "PSF/560_c.tif", 0.88, True)
        _common.otf_mask_figure(data_path + "PSF/488_2_c.tif", 0.88, True)

    # the decon settings both runs share
    common = {
        "xy_pixel_size": 0.108,
        "dz": 0.3,
        "dz_psf": 0.3,
        "parse_setting_file": False,
        "channel_patterns": ["CamB_ch0"],
        "psf_fullpaths": [data_path + "PSF/488_2_c.tif"],
        "background": 100,
        "edge_erosion": 0,
        "save_16bit": True,
        "zarr_file": False,
        "save_zarr": False,
        "cpus_per_task": 24,
        "parse_cluster": False,
        "large_file": False,
        "gpu_job": _common.can_use_gpu(),
        "debug": False,
        "config_file": "",
        "gpu_config_file": "",
        "mcc_mode": False,
    }
    # MATLAB's demo also passes 'Reverse', which the wrapper does not use

    with steps("Step 3.2: OMW decon, 2 iterations"):
        decon_dataset(
            [data_path + "/"],
            result_dir_name="matlab_decon_omw",
            rl_method="omw",
            wiener_alpha=0.005,
            otf_cum_thresh=0.87,
            skewed=True,
            decon_iter=2,
            **common,
        )

    with steps("Step 4.2: traditional (Biggs) RL decon, 30 iterations"):
        decon_dataset(
            [data_path + "/"],
            result_dir_name="matlab_decon_traditional",
            rl_method="simplified",
            decon_iter=30,
            **common,
        )

    with steps("Step 5: deskew/rotate the raw and deconvolved tiles"):
        deskew_rotate_dataset(
            [data_path, data_path + "matlab_decon_omw/", data_path + "matlab_decon_traditional/"],
            xy_pixel_size=0.108,
            dz=0.3,
            reverse=True,
            channel_patterns=["CamB_ch0"],
            large_file=False,
            zarr_file=False,
            save_zarr=False,
            save_16bit=True,
            parse_cluster=False,
            master_compute=True,
            config_file="",
            mcc_mode=False,
        )

    fsn = _common.TILE_B
    # MATLAB's crop box [1040, 250, 133, 1240, 450, 133], one z slice of the DSR volumes
    bbox = (132, 1039, 249, 133, 1240, 450)
    with steps("Step 6: crops of the deskewed/rotated results"):
        crops = {}
        for title, folder in (
            ("Raw", ""),
            ("Traditional (Biggs) 30 iters", "matlab_decon_traditional/"),
            ("OMW 2 iters", "matlab_decon_omw/"),
        ):
            crops[title] = crop_3d(read_tiff(f"{data_path}{folder}DSR/{fsn}.tif"), bbox)[0]
        crop_figure("dsr_crops", crops, 100)

    # ---- the decon functions on one tile, in memory
    fn = data_path + fsn + ".tif"
    psf_fn = data_path + "PSF/488_2_c.tif"
    with steps("remove the background and clean the PSF"):
        im = read_tiff(fn)
        psf = read_tiff(psf_fn)
        bg = 100
        # MATLAB's max(im - bg, 0) for uint16 (the subtraction saturates at 0)
        im = np.maximum(im, bg) - np.uint16(bg)
        psf = psf_gen(psf, 0.3, 0.3, 1.5, "masked")
        psf = psf / np.asarray(msum(psf), dtype=psf.dtype)

    with steps("traditional (Biggs) RL decon, 30 iterations"):
        rl_decon = decon_lucy(im, psf, 30, use_gpu=None).decon
    _common.save_result(DEMO, "rl_decon", rl_decon)

    steps.skip("WB RL decon, 2 iterations", "the Wiener-Butterworth back projector is not ported")

    with steps("OMW back projector and RL decon, 2 iterations"):
        _common.otf_mask_figure(psf_fn, 0.88, True)
        bp = omw_backprojector(psf, 0.005, True, otf_cum_thresh=0.88, hann_win_bounds=(0.8, 1.0))
        omw_decon = decon_lucy_omw(im, psf, bp.b_omw, 2, use_gpu=None).decon
    _common.save_result(DEMO, "b_omw", bp.b_omw)
    _common.save_result(DEMO, "omw_decon", omw_decon)

    with steps("deskew/rotate the in-memory results and crop them"):
        crops = {}
        for title, vol in (("Raw", im), ("Traditional (Biggs) 30 iters", rl_decon),
                           ("OMW 2 iters", omw_decon)):  # fmt: skip
            dsr = deskew_rotate_volume(
                vol,
                32.45,
                0.3,
                0.108,
                reverse=True,
                crop=True,
                objective_scan=False,
                resample_factor=None,
                interp_method="linear",
                save_16bit=True,
            )
            crops[title] = crop_3d(dsr, bbox)[0]
            del dsr
        for title, name in zip(crops, ("raw_dsr_crop", "rl_dsr_crop", "omw_dsr_crop"), strict=True):
            _common.save_result(DEMO, name, crops[title])
        crop_figure("inmem_crops", crops, None)


if __name__ == "__main__":
    _common.run_demo(main)
