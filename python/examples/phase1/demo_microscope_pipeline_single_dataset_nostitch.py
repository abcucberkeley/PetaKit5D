"""The automatic processing pipeline on one dataset without stitching (MATLAB:
``demos/demo_microscope_pipeline_single_dataset_nostitch.m``): flat-field correction and the
combined deskew/rotation of every tile. Results: ``DSR/`` in the dataset folder.
"""

from __future__ import annotations

import _common
from demo_microscope_pipeline_single_dataset import pipeline_settings

DEMO = "demo_microscope_pipeline_single_dataset_nostitch"


def main() -> None:
    from petakit5d.pipeline import automatic_processing

    steps = _common.Steps(DEMO)
    data_path = _common.data_path()

    with steps("flat-field correction and deskew/rotation"):
        automatic_processing(
            data_path, stitch=False, dsr_combined=True, **pipeline_settings(data_path)
        )


if __name__ == "__main__":
    _common.run_demo(main)
