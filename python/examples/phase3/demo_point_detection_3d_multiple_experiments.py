"""The point-detection demo for multiple experiments (MATLAB:
``demos/demo_point_detection_3d_multiple_experiments.m``): the three experiments of the
dataset's condition loaded together, the PSF sigmas estimated, each movie deskewed, detected,
tracked and its tracks processed and rotated (``XR_runDetTrack3d``:
``run_detection_tracking``), then the lifetime analysis of all the experiments and the
intensity cohorts, with the demo's parameters.

    python demo_point_detection_3d_multiple_experiments.py

The data: ``PETAKIT5D_DEMO_DATA`` (``_phase3``). Each experiment gets the step-by-step demo's
results in its own folders; the lifetime analysis adds ``lifetimeData.mat`` to each results
folder and ``Parameters.mat`` to the first's.

MATLAB's demo passes ``'zoffsetCorrection', true``, which is not an option of
``XR_runDetTrack3d``: its parser drops it, so this demo runs without the channels' offset
correction (``'XZoffsetCorrection'`` false, the default), unlike the integrated demo. The
example leaves it out and gives MATLAB's results.
"""

from __future__ import annotations

import _phase3

# isort: split
import _common

DEMO = "demo_point_detection_3d_multiple_experiments"


def detection_tracking(data, sigma) -> None:
    """Step 3: deskew, detection, tracking, track processing and track rotation of every
    movie (``XR_runDetTrack3d``)"""
    from petakit5d.detection import run_detection_tracking

    run_detection_tracking(
        data,
        sigma=sigma,
        detection_method="lowSNR",
        aname=_phase3.ANAME,
        fit_mixtures=False,
        background_correction=True,
    )


def lifetime_analysis(data) -> dict:
    """Step 4: the lifetime analysis of all the experiments (``runLifetimeAnalysis3D``): the
    tracks' lifetimes, classified by their intensity and by the other channel's signal, in
    ``lftRes``; ``lifetimeData.mat`` in each results folder, ``Parameters.mat`` in the first"""
    from petakit5d.detection.lifetime import run_lifetime_analysis

    return run_lifetime_analysis(data, analysis_path=_phase3.ANAME, overwrite=True)


def intensity_cohorts(data) -> list[dict]:
    """Step 4: the intensity cohorts (``plotIntensityCohorts``: the mean intensity over time of
    the tracks of each lifetime cohort, the numbers of its figure), one per combination of the
    other channel's significance"""
    from petakit5d.detection.lifetime import plot_intensity_cohorts

    return plot_intensity_cohorts(data, analysis_path=_phase3.ANAME)


def main() -> None:
    steps = _common.Steps(DEMO)
    rt = _phase3.condition_path()
    with steps("Step 1: load data"):
        data = _phase3.load_data(rt)
    _phase3.save_values(DEMO, "load", {"data": data})
    with steps("Step 2: PSF sigmas"):
        sigma = _phase3.psf_sigmas(data, _phase3.psf_paths())
    _phase3.save_values(DEMO, "sigma", {"sigma_mat": sigma})
    with steps("Step 3: deskew, detection, tracking and track processing"):
        detection_tracking(data, sigma)
    with steps("Step 4: lifetime analysis"):
        lft_res = lifetime_analysis(data)
    _phase3.save_values(DEMO, "lifetime", {"lftRes": lft_res})
    with steps("Step 4: intensity cohorts"):
        cohorts = intensity_cohorts(data)
    _phase3.save_values(DEMO, "cohorts", {"cohorts": cohorts})


if __name__ == "__main__":
    _common.run_demo(main)
