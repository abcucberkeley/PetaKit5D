"""Runs the paper's demos in order (MATLAB: ``demos/demo_run_all_demos_for_the_paper.m``), in
one process and on one copy of the data: later demos use the results of earlier ones (the tools
use the deconvolved tiles; demo_large_scale_processing the stitched volumes).

    python run_all_demos_for_the_paper.py [demo ...]

With names, only those demos run, in the paper's order.
"""

from __future__ import annotations

import importlib
import sys
import time

import _common

PAPER_DEMOS = (
    ("Generic computing framework", "demo_generic_computing_framework"),
    ("Fast Tiff/Zarr readers and writers", "demo_fast_tiff_zarr_readers_writers"),
    ("Geometric transformation (deskew/rotation)", "demo_geometric_transformation"),
    ("Deconvolution", "demo_rl_deconvolution"),
    ("Stitching", "demo_zarr_stitching"),
    ("Large-scale processing", "demo_large_scale_processing"),
    ("Useful tools", "demo_useful_tools"),
    ("Phase and 2-Photon stitching", "demo_phase_and_2photon_stitching"),
    ("Widefield and Confocal deconvolution", "demo_widefield_and_confocal_deconvolution"),
)


def main(names: list[str]) -> None:
    unknown = set(names) - {m for _, m in PAPER_DEMOS}
    if unknown:
        sys.exit(f"unknown demos: {sorted(unknown)}")
    for title, module in PAPER_DEMOS:
        if names and module not in names:
            continue
        print(f"\n{title} demo...", flush=True)
        t0 = time.perf_counter()
        importlib.import_module(module).main()
        print(f"\n{title} demo: {time.perf_counter() - t0:.1f} s", flush=True)
    print("\nAll demos are finished!")


if __name__ == "__main__":
    _common.run_demo(lambda: main(sys.argv[1:]))
