"""The automatic processing pipeline on one dataset (MATLAB:
``demos/demo_microscope_pipeline_single_dataset.m``): flat-field correction, combined
deskew/rotation of every tile and stitching of each time point, as the microscope writes them
(here with streaming off, so the finished dataset is processed once).

``parse_cluster=True`` as in the demo: Slurm jobs where Slurm is available, else this machine.
Results: ``DSR/`` and the stitching result folder in the dataset folder. The variant without
stitching is demo_microscope_pipeline_single_dataset_nostitch.
"""

from __future__ import annotations

import _common

DEMO = "demo_microscope_pipeline_single_dataset"


def pipeline_settings(data_path: str) -> dict:
    """the settings both single-dataset demos share"""
    ff = "FF/averaged/ff_{}_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif"  # noqa: E501
    return {
        "xy_pixel_size": 0.108,
        "dz": 0.3,
        "reverse": True,
        "channel_patterns": ["CamB_ch0", "CamA_ch0"],
        "save_16bit": [True, True],
        "overwrite": False,
        "streaming": False,
        "parse_cluster": True,
        "cpus_per_task": 24,
        "config_file": "",
        "mcc_mode": False,
        "deskew": True,
        "rotate": True,
        "ff_correction": True,
        "lower_limit": 0.4,
        "const_offset": 1,
        "ff_image_paths": [
            data_path + "/" + ff.format("CamB"),
            data_path + "/" + ff.format("CamA"),
        ],
        "background_paths": [
            _common.background_path(data_path, "CamB", "/FF/"),
            _common.background_path(data_path, "CamA", "/FF/"),
        ],
        # streaming waits: files must be 1 to 10 min old, and at most 10 rounds without new ones
        "min_modify_time": 1,
        "max_modify_time": 10,
        "max_wait_loop_num": 10,
    }


def main() -> None:
    from petakit5d.pipeline import automatic_processing

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()

    with steps("flat-field correction, deskew/rotation and stitching"):
        automatic_processing(
            data_path,
            stitch=True,
            image_list_fullpaths=data_path + "ImageList_from_encoder.csv",
            axis_order="-x,y,z",
            **pipeline_settings(data_path),
        )


if __name__ == "__main__":
    _common.run_demo(main)
