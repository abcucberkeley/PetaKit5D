"""The automatic processing pipeline on several datasets at once (MATLAB:
``demos/demo_microscope_pipeline_multiple_dataset.m``): deskew/rotation and stitching of three
experiments, each with its own image list.

The MATLAB demo points at the lab's own data; here the experiments are the folders given on the
command line (each with its ``ImageList_<folder name>.csv``), and the flat-field correction is off
as in the demo::

    python demo_microscope_pipeline_multiple_dataset.py /mnt/data/Exp01 /mnt/data/Exp02
"""

from __future__ import annotations

import os
import sys

import _common

DEMO = "demo_microscope_pipeline_multiple_dataset"


def main() -> None:
    from petakit5d.pipeline import automatic_processing

    data_paths = [p.rstrip("/") + "/" for p in sys.argv[1:]]
    if not data_paths:
        sys.exit("usage: demo_microscope_pipeline_multiple_dataset.py <dataset folder>...")
    image_lists = [f"{p}ImageList_{os.path.basename(p.rstrip('/'))}.csv" for p in data_paths]

    steps = _common.Steps(DEMO)
    with steps("deskew/rotation and stitching of the datasets"):
        automatic_processing(
            data_paths,
            xy_pixel_size=0.108,
            dz=0.3,
            reverse=True,
            channel_patterns=["CamB_ch0", "CamA_ch0"],
            save_16bit=[True, True],
            overwrite=False,
            streaming=False,
            deskew=True,
            rotate=True,
            stitch=True,
            image_list_fullpaths=image_lists,
            axis_order="-x,y,z",
            dsr_combined=True,
            # flat-field images per channel (off here): the averaged flat fields and the
            # cameras' dark images
            ff_correction=False,
            lower_limit=0.4,
            const_offset=1,
            ff_image_paths=[
                "/mnt/calibration/FF/AVG_488nm.tif",
                "/mnt/calibration/FF/AVG_560nm.tif",
            ],
            background_paths=["/mnt/calibration/AVG_Cam_B.tif", "/mnt/calibration/AVG_Cam_A.tif"],
            min_modify_time=1,
            max_modify_time=10,
            max_wait_loop_num=10,
            parse_cluster=True,
            cpus_per_task=24,
            config_file="",
            mcc_mode=False,
        )


if __name__ == "__main__":
    _common.run_demo(main)
