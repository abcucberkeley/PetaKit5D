"""Runs the paper's demos in order (MATLAB: ``demos/demo_run_all_demos_for_the_paper.m``) with
their Phase 2 steps, in one process and on one copy of the data: each demo's example of
``examples/phase1`` and then, where this folder has one, the example of its Phase 2 steps (the
data tools demo's FFT spectra, PSF and FSC analyses and chromatic shift correction).

    python run_all_demos_for_the_paper.py [demo ...]

With names, only those demos run, in the paper's order (``demo_useful_tools`` runs both parts).
The phase 1 examples of the phase and 2-photon stitching and of the widefield and confocal
deconvolution have their Phase 2 steps already: the phase flat field (with the ``basic`` extra)
and the OTF mask figures (with the ``figures`` extra).
"""

from __future__ import annotations

import importlib
import sys
import time

import _phase1

# isort: split
import _common

# the phase 1 example -> this folder's example of its Phase 2 steps
PHASE2_STEPS = {"demo_useful_tools": "demo_useful_tools"}


def paper_demos() -> tuple[tuple[str, str], ...]:
    """the paper's demos (title, phase 1 example), as phase 1's run-all lists them"""
    return importlib.import_module("run_all_demos_for_the_paper").PAPER_DEMOS


def main(names: list[str]) -> None:
    demos = paper_demos()
    unknown = set(names) - {m for _, m in demos}
    if unknown:
        sys.exit(f"unknown demos: {sorted(unknown)}")
    for title, module in demos:
        if names and module not in names:
            continue
        print(f"\n{title} demo...", flush=True)
        t0 = time.perf_counter()
        importlib.import_module(module).main()
        if module in PHASE2_STEPS:
            _phase1.load(PHASE2_STEPS[module]).main()
        print(f"\n{title} demo: {time.perf_counter() - t0:.1f} s", flush=True)
    print("\nAll demos are finished!")


if __name__ == "__main__":
    _common.run_demo(lambda: main(sys.argv[1:]))
