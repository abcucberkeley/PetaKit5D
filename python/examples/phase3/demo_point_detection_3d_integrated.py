"""The point-detection demo with the integrated function (MATLAB:
``demos/demo_point_detection_3d_integrated.m``): the dataset's first experiment loaded, the PSF
sigmas estimated, then the deskew, the channels' x and z offsets, the detection, the tracking,
the track processing and the tracks' rotation in one call (``XR_runDetTrack3d``:
``run_detection_tracking``), with the demo's parameters.

    python demo_point_detection_3d_integrated.py

The data: ``PETAKIT5D_DEMO_DATA`` (``_phase3``). The results are those of the step-by-step
demo (``demo_point_detection_3d_by_step.py``), in the same folders, with the wrapper's
defaults: the detection runs frame by frame (``ParallelMethod`` 'lock_distributed', the
default: ``run_detection_frame_parallel``, which also keeps each frame's detections in
``Analysis_update_3/Detections/``), and files that exist are kept (``Overwrite`` false).
"""

from __future__ import annotations

import _phase3

# isort: split
import _common

DEMO = "demo_point_detection_3d_integrated"


def detection_tracking(data, sigma) -> None:
    """Step 3: deskew, XZ offsets, detection, tracking, track processing and track rotation
    (``XR_runDetTrack3d``)"""
    from petakit5d.detection import run_detection_tracking

    run_detection_tracking(
        data,
        sigma=sigma,
        detection_method="lowSNR",
        aname=_phase3.ANAME,
        fit_mixtures=False,
        background_correction=True,
        x_zoffset_correction=True,
    )


def main() -> None:
    steps = _common.Steps(DEMO)
    rt = _phase3.experiment_path(1)
    with steps("Step 1: load data"):
        data = _phase3.load_data(rt)
    _phase3.save_values(DEMO, "load", {"data": data})
    with steps("Step 2: PSF sigmas"):
        sigma = _phase3.psf_sigmas(data, _phase3.psf_paths())
    _phase3.save_values(DEMO, "sigma", {"sigma_mat": sigma})
    with steps("Step 3: deskew, detection, tracking and track processing"):
        detection_tracking(data, sigma)


if __name__ == "__main__":
    _common.run_demo(main)
