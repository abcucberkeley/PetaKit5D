"""The data of the FSC benchmark (B05): the demo's decon iterations, as demo_useful_tools.m makes
them before its FSC section: the port's XR_decon_data_wrapper with the demo's parameters (OMW, 20
iterations, every 2nd saved with debug) of the CamB ch0 000x_003y tile of the demo cell image
dataset (PETAKIT5D_DEMO_DATA), its PSF 488_2_c.tif; the ten Iter_*.tif in <folder>/matlab/iters/
and <folder>/python/iters/ (the same files).

usage: python fsc_b05_data.py FOLDER
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil

TILE = (
    "Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_"
    "000z_0000t"
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    a = ap.parse_args()
    from petakit5d import matlab

    demo = os.path.join(os.environ["PETAKIT5D_DEMO_DATA"], "PetaKit5D_demo_cell_image_dataset")
    src = os.path.join(a.folder, "demo") + "/"
    os.makedirs(src + "PSF", exist_ok=True)
    shutil.copy(os.path.join(demo, TILE + ".tif"), src)
    shutil.copy(os.path.join(demo, "PSF", "488_2_c.tif"), src + "PSF/")
    matlab.XR_decon_data_wrapper(
        [src],
        "resultDirName", "matlab_decon_omw_debug",
        "xyPixelSize", 0.108,
        "dz", 0.3,
        "ChannelPatterns", [TILE],
        "PSFFullpaths", [src + "PSF/488_2_c.tif"],
        "dzPSF", 0.3,
        "parseSettingFile", False,
        "RLmethod", "omw",
        "wienerAlpha", 0.05,
        "OTFCumThresh", 0.9,
        "skewed", True,
        "Background", 100,
        "DeconIter", 20,
        "fixIter", True,
        "EdgeErosion", 0,
        "Save16bit", True,
        "zarrFile", False,
        "saveZarr", False,
        "parseCluster", False,
        "largeFile", False,
        "GPUJob", False,
        "saveStep", 2,
        "debug", True,
        "cpusPerTask", 24,
        "mccMode", False,
    )  # fmt: skip
    iters = sorted(glob.glob(src + f"matlab_decon_omw_debug/{TILE}_debug/Iter_*.tif"))
    if len(iters) != 10:
        raise SystemExit(f"expected 10 decon iterations, found {len(iters)}")
    for side in ("matlab", "python"):
        dst = os.path.join(a.folder, side, "iters")
        os.makedirs(dst, exist_ok=True)
        for f in iters:
            shutil.copy(f, dst)
    shutil.rmtree(src)


if __name__ == "__main__":
    main()
