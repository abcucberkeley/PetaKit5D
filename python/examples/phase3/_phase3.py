"""Shared helpers of the Phase 3 examples: the point-detection demos
(``demos/demo_point_detection_3d_*.m``).

The demos analyse the point-detection demo dataset: one condition folder with three experiments
(``Ex1``, ``Ex2``, ``Ex3``), each with two channels of 60 frames (``ch1``, ``ch2``) taken by
stage scanning, and the two channels' bead images (``560totalPSF.tif``, ``642totalPSF.tif``) in
the condition folder. The examples find it through ``PETAKIT5D_DEMO_DATA``: the folder that
holds ``PetaKit5D_point_detection_demo_dataset``, or that folder itself (by default
``~/Downloads/``, as in the other demos). They write their results into the experiment folders,
as MATLAB does, so point it at a copy you can write to.

The MATLAB demos load the data interactively: they ask for the number of channels, the channel
folders and the fluorescent markers. The examples pass what the demos' comments say to enter:
two channels, ``ch1`` and ``ch2``, with the markers ``gfp`` and ``rfp``.

``examples/phase1`` is put first on the import path for its ``_common`` helpers (the step times
and ``PETAKIT5D_DEMO_TIMINGS``, ``run_demo``). ``PETAKIT5D_DEMO_SAVE_DIR`` names a folder the
demos' in-memory results are written to, as MAT-files (``<demo>/<name>.mat``).
"""

from __future__ import annotations

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE1 = os.path.normpath(os.path.join(HERE, os.pardir, "phase1"))
if not sys.path or sys.path[0] != PHASE1:
    if PHASE1 in sys.path:
        sys.path.remove(PHASE1)
    sys.path.insert(0, PHASE1)

DATASET = "PetaKit5D_point_detection_demo_dataset"
# what the demos' comments say to enter when the loader asks
CHANNELS = ["ch1", "ch2"]
MARKERS = ["gfp", "rfp"]
# the bead images of channel 1 (560 nm) and channel 2 (642 nm), in the condition folder
PSFS = ("560totalPSF.tif", "642totalPSF.tif")
# the name of the results folder, under the master channel's folder
ANAME = "Analysis_update_3"


def dataset_path() -> str:
    """the demo dataset's folder (without a trailing separator)"""
    import _common

    path = _common.dest_path(DATASET) + DATASET
    if not os.path.isdir(path):
        sys.exit(
            f"The point-detection demo dataset {DATASET} is not in {_common.dest_path(DATASET)}: "
            "set PETAKIT5D_DEMO_DATA to the folder that holds it (or to the folder itself)."
        )
    return path


def _folders(parent: str, pattern: str = "") -> list[str]:
    return sorted(
        e.name
        for e in os.scandir(parent)
        if e.is_dir() and not e.name.startswith(".") and re.match(pattern, e.name)
    )


def condition_path() -> str:
    """the condition folder, the dataset's only folder (the multiple-experiments demo's root
    ``rt``), written as the demos write their roots: without a trailing separator"""
    root = dataset_path()
    names = _folders(root)
    if not names:
        sys.exit(f"{root} holds no condition folder.")
    return root + "/" + names[0]


def experiment_path(n: int = 1) -> str:
    """the folder of experiment ``n`` of the condition (``Ex<n>...``: the root ``rt`` of the
    step-by-step and integrated demos, which take the first)"""
    cond = condition_path()
    names = _folders(cond, rf"Ex{n}(?!\d)")
    if not names:
        sys.exit(f"{cond} has no experiment {n} (a folder Ex{n}...).")
    return cond + "/" + names[0]


def psf_paths() -> list[str]:
    """the bead images of the two channels (the demos' ``ch1_psf_filename``,
    ``ch2_psf_filename``)"""
    cond = condition_path()
    return [cond + "/" + name for name in PSFS]


def load_data(rt: str):
    """Step 1: the movies of ``rt`` (one experiment's folder, or the condition folder with all
    of them), with the demos' channels and markers (``XR_loadConditionData3D(rt)`` and the
    answers its prompts ask for)"""
    from petakit5d.detection import load_condition_data

    return load_condition_data(rt, CHANNELS, MARKERS)


def psf_sigmas(data, psf_files):
    """Step 2: the Gaussian sigmas of the channels' PSFs from their bead images
    (``GU_estimateSigma3D(file, [])``), as the rows ``[sigma_xy, sigma_z]`` of the demos'
    ``sigma_mat``: the bead images are about isotropic, the deskewed frames are not, so
    ``sigma_z`` is divided by their z anisotropy (``data(1).zAniso``)"""
    import numpy as np

    from petakit5d.detection.sigma import estimate_sigma3d

    z_aniso = float(data[0].z_aniso)
    rows = []
    for path in psf_files:
        sigma_xy, sigma_z = estimate_sigma3d(path)[:2]
        rows.append([sigma_xy, sigma_z / z_aniso])
    return np.array(rows, dtype=np.float64)


def _struct_array(v):
    """a condition, or a list of dicts with the same fields, as MATLAB's 1 x n struct array (one
    element: its dict); anything else as it is"""
    from petakit5d.detection.condition import Condition, struct_array
    from petakit5d.io.matfile import StructArray

    if isinstance(v, Condition):
        return struct_array(v)
    if isinstance(v, list) and v and all(isinstance(e, dict) for e in v):
        if len({tuple(e) for e in v}) == 1:
            return v[0] if len(v) == 1 else StructArray((1, len(v)), list(v[0]), v)
    return v


def save_values(demo: str, name: str, values: dict) -> None:
    """writes in-memory results to ``PETAKIT5D_DEMO_SAVE_DIR/<demo>/<name>.mat`` (if set), as
    MATLAB's ``save`` writes them (version 7); a condition, or a list of structs (the intensity
    cohorts), as MATLAB's struct array"""
    root = os.environ.get("PETAKIT5D_DEMO_SAVE_DIR")
    if not root:
        return
    from petakit5d.io import write_mat

    out = {k: _struct_array(v) for k, v in values.items()}
    os.makedirs(os.path.join(root, demo), exist_ok=True)
    write_mat(os.path.join(root, demo, name + ".mat"), out, version="7")
