"""Shared helpers of the Phase 1 demo examples.

Each example in this folder follows one MATLAB demo (``demos/demo_*.m``) step by step, with the
same parameters, through the snake_case API. Python takes arrays as ``(z, y, x)`` and boxes as
``(z0, y0, x0, z1, y1, x1)`` (0-based, end-exclusive), so the comments give MATLAB's values where
they differ. Steps whose functions are not in Phase 1 of the port are left out, and say so.

Environment:

``PETAKIT5D_DEMO_DATA``
    the folder that holds the demo datasets (``PetaKit5D_demo_cell_image_dataset``,
    ``PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets``), or one dataset folder; by default
    ``~/Downloads/``, as in the MATLAB demos. The demos write their results into the dataset
    folders, so point it at a copy you can write to.
``PETAKIT5D_DEMO_TIMINGS``
    a file the time of every step is appended to (one JSON object per line).
``PETAKIT5D_DEMO_SAVE_DIR``
    a folder the demos' in-memory results are written to as tiff files (``<demo>/<name>.tif``).
"""

from __future__ import annotations

import contextlib
import glob
import importlib.util
import json
import os
import sys
import time
from collections.abc import Iterator

CELL = "PetaKit5D_demo_cell_image_dataset"
OTHER = "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets"

# the tile the demos look at closely
TILE_A = (
    "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_"
    "000x_003y_000z_0000t"
)
TILE_B = TILE_A.replace("CamA", "CamB")


def dest_path(dataset: str = CELL) -> str:
    """the folder that holds ``dataset`` (MATLAB's destPath), with a trailing slash"""
    p = os.path.expanduser(os.environ.get("PETAKIT5D_DEMO_DATA") or "~/Downloads/")
    p = p.replace("\\", "/")
    if os.path.basename(p.rstrip("/")) == dataset:
        p = os.path.dirname(p.rstrip("/"))
    return p if p.endswith("/") else p + "/"


def data_path(dataset: str = CELL) -> str:
    """the dataset folder, built as the MATLAB demos build it (``[destPath, '/<dataset>/']``)"""
    path = dest_path(dataset) + "/" + dataset + "/"
    if not os.path.isdir(path):
        sys.exit(
            f"The demo dataset {dataset} is not in {dest_path(dataset)}: download it (see the "
            "MATLAB demos' data downloaders) and set PETAKIT5D_DEMO_DATA to its parent folder."
        )
    return path


def background_path(data_path: str, cam: str, ff: str = "FF/") -> str:
    """the dataset's averaged dark image of camera ``cam``, the demos' background image:
    ``data_path + ff + '<folder>/AVG_DF400_<cam>_10ms.tif'`` (the path as the MATLAB demos write
    it, so that parameter files hold the same text)"""
    name = f"AVG_DF400_{cam}_10ms.tif"
    hits = sorted(glob.glob(os.path.join(data_path + ff, "*", name)))
    if not hits:
        sys.exit(f"The demo dataset in {data_path} has no background image for {cam}.")
    return f"{data_path}{ff}{os.path.basename(os.path.dirname(hits[0]))}/{name}"


class Steps:
    """times the steps of a demo: ``with steps("name"): ...`` prints the step's time and
    appends it to ``PETAKIT5D_DEMO_TIMINGS``"""

    def __init__(self, demo: str):
        self.demo = demo
        self.total = 0.0

    @contextlib.contextmanager
    def __call__(self, name: str) -> Iterator[None]:
        print(f"\n[{self.demo}] {name}...", flush=True)
        t0 = time.perf_counter()
        yield
        dt = time.perf_counter() - t0
        self.total += dt
        print(f"[{self.demo}] {name}: {dt:.2f} s", flush=True)
        log = os.environ.get("PETAKIT5D_DEMO_TIMINGS")
        if log:
            with open(log, "a", newline="\n") as f:
                f.write(json.dumps({"demo": self.demo, "step": name, "seconds": dt}) + "\n")

    def skip(self, name: str, reason: str) -> None:
        print(f"\n[{self.demo}] {name}: skipped ({reason})", flush=True)


def save_result(demo: str, name: str, im) -> None:
    """writes an in-memory result to ``PETAKIT5D_DEMO_SAVE_DIR/<demo>/<name>.tif`` (if set)"""
    root = os.environ.get("PETAKIT5D_DEMO_SAVE_DIR")
    if not root:
        return
    from petakit5d.io import write_tiff

    os.makedirs(os.path.join(root, demo), exist_ok=True)
    write_tiff(im, os.path.join(root, demo, name + ".tif"))


def replicate_frames(im, nframe: int):
    """MATLAB's ``repmat(im, 1, 1, ceil(nframe / nz))`` cropped to ``nframe`` z slices: the
    ``(z, y, x)`` stack repeated along z, copied in threads (NumPy's copies release the GIL)"""
    from concurrent.futures import ThreadPoolExecutor

    import numpy as np

    from petakit5d.core.threads import available_cpus

    nz = im.shape[0]
    out = np.empty((nframe, *im.shape[1:]), dtype=im.dtype)

    def part(z0: int) -> None:
        z1 = min(z0 + nz, nframe)
        out[z0:z1] = im[: z1 - z0]

    with ThreadPoolExecutor(max_workers=available_cpus()) as ex:
        list(ex.map(part, range(0, nframe, nz)))
    return out


def figures_available() -> bool:
    """whether matplotlib is installed (the ``figures`` extra), for the demos' figures"""
    return importlib.util.find_spec("matplotlib") is not None


def can_use_gpu() -> bool:
    """MATLAB's canUseGPU: a GPU this process may use"""
    from petakit5d.gpu.device import usable_gpu_count

    return usable_gpu_count()[0] > 0


def show_images(demo: str, name: str, images: dict, clims: dict) -> None:
    """the demo's figure of 2D images side by side (gray, with the given limits), if matplotlib
    is installed: saved to ``PETAKIT5D_DEMO_SAVE_DIR/<demo>/<name>.png`` if that is set, else
    shown"""
    if not figures_available():
        print(f"(figure {name} left out: matplotlib is not installed)", flush=True)
        return
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(images), figsize=(4 * len(images), 4), squeeze=False)
    for ax, (title, im) in zip(axes[0], images.items(), strict=True):
        ax.imshow(im, cmap="gray", vmin=clims[title][0], vmax=clims[title][1])
        ax.set_title(title)
        ax.axis("off")
    root = os.environ.get("PETAKIT5D_DEMO_SAVE_DIR")
    if root:
        os.makedirs(os.path.join(root, demo), exist_ok=True)
        fig.savefig(os.path.join(root, demo, name + ".png"))
        plt.close(fig)
    else:
        plt.show()


def otf_mask_figure(psf_fn: str, otf_cum_thresh: float, skewed: bool, **kwargs) -> None:
    """XR_visualize_OTF_mask_segmentation's figure, if matplotlib is installed"""
    if not figures_available():
        print("(OTF mask figure left out: matplotlib is not installed)", flush=True)
        return
    from petakit5d.tools import visualize_otf_mask_segmentation

    visualize_otf_mask_segmentation(psf_fn, otf_cum_thresh, skewed, **kwargs)


def run_demo(main) -> None:
    """runs a demo's ``main`` and prints its total time"""
    t0 = time.perf_counter()
    main()
    print(f"\nDone in {time.perf_counter() - t0:.1f} s.", flush=True)
