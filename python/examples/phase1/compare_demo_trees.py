"""Compares the output trees of the demo examples with those of the MATLAB demos, per demo.

    python compare_demo_trees.py MATLAB_ROOT PYTHON_ROOT --inputs INPUTS [options]

``MATLAB_ROOT`` and ``PYTHON_ROOT`` are the folders that hold the dataset folders after the
runs (the demos' ``Downloads``); ``INPUTS`` holds the untouched datasets, whose files are not
outputs. Every output is assigned to the demo that writes it (``RULES``, first match) and
compared with ``petakit5d.testing.compare``: file and folder names, tiff voxels and tags, zarr
voxels and metadata (MATLAB's ``(y, x, z)`` arrays against the port's native ones), the values
in .mat and JSON files and text. Each side's own root paths are masked (``--matlab-home``,
``--python-home``, ``--matlab-code``, ``--python-code``: the folders the runs' paths start with;
a leading ``~`` is the home folder). Empty folders are not compared.
MATLAB outputs of steps the port leaves out are listed by ``NOT_PORTED`` and not compared (the
demos' Phase 2 steps: ``examples/phase2/compare_phase2_trees.py``). The per-pair registration
results (``xcorr/<volume>/xcorr_tile_*.mat``) are compared within ``docs/parity/index.md``'s
rows for them: their correlation peaks within P13 (2D) and P7 (3D), their shifts exact.

Options: ``--only REGEX`` (names to compare), ``--demo NAME`` (one demo's run: every output is
its), ``--out DIR`` (``compare_<label>.txt`` and ``.json``), ``--label``, ``--workers``. Exit
status 1 if anything differs.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor

CELL = "PetaKit5D_demo_cell_image_dataset"
OTHER = "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets"

_C = re.escape(CELL) + "/"
_O = re.escape(OTHER) + "/"
# (regex on the name relative to the root, demo); first match wins. In the paper's order the
# demos write into the same dataset folders (demo_RL_deconvolution adds the CamB tiles' tiffs
# to DSR/, whose zarr files demo_geometric_transformation wrote).
RULES = [
    (_C + r"Cropped_DSR/", "generic_computing_framework"),
    (_C + r"replicated/[^/]*frame_number", "fast_tiff_zarr_readers_writers"),
    (_C + r"(replicated/|DSR_combined/|DS_separated/|DSR_separated/)", "geometric_transformation"),
    (_C + r"matlab_decon_omw/(DS|DSR)_for_useful_tools", "useful_tools"),
    (
        _C + r"(resampled|Cropped|MIPs|matlab_decon_omw_debug|tiff_to_zarr|Unmixed|FFT|FSC"
        r"|Chromatic_Shift_Corrected|ImageList_from_tile_list|PSF/)",
        "useful_tools",
    ),
    (_C + r"DSR/(resaved_zarr|Unmixed)/", "useful_tools"),
    (_C + r"matlab_stitch/matlab_decon_omw/", "large_scale_processing"),
    (_C + r"(matlab_stitch/DSR/|matlab_stitch_dsr/)", "zarr_stitching"),
    (_C + r"(matlab_stitch/|tmp/|zarr/)", "skewed_space_stitching"),
    (_C + r"DSR/[^/]*\.zarr", "geometric_transformation"),
    (_C + r"DSR/MIPs/[^/]*CamA", "geometric_transformation"),
    (_C + r"(matlab_decon_omw/|matlab_decon_traditional/|DSR/)", "RL_deconvolution"),
    (_O + r"(Phase|2P)/", "phase_and_2photon_stitching"),
    (_O + r"(Widefield|Confocal)/", "widefield_and_confocal_deconvolution"),
]

# MATLAB outputs of steps that are not in Phase 1 of the port: (regex, step)
NOT_PORTED = [
    (_C + r"(replicated/)?(DS_separated?|DSR_separated?)/", "XR_deskewRotateFrame_old (legacy)"),
    (_C + r"(replicated/)?DS/", "XR_deskewRotateFrame_old (legacy)"),
    (_C + r"(.*/)?FFT/", "XR_fftSpectrumComputingWrapper (Phase 2)"),
    (_C + r"PSF/(DS/|.*[aA]nalysis)", "XR_psf_analysis_wrapper (Phase 2)"),
    (_C + r"(.*/)?FSC|.*_debug/figures/fsc_", "XR_FSC_analysis_wrapper (Phase 2)"),
    (_C + r"(.*/)?imaris/", "XR_imaris_conversion_data_wrapper (Phase 2)"),
    (_C + r"Chromatic_Shift_Corrected/", "XR_chromatic_shift_correction_data_wrapper (Phase 2)"),
    # MATLAB's figures (png, fig): the port draws them with the optional figures extra
    (r".*\.(png|fig)$", "figures (MATLAB graphics; the port's figures extra)"),
]

# MATLAB outputs of steps that need an optional extra: not compared when Python has no such file
OPTIONAL = [
    (_O + r"Phase/(ff_corrected_1/|.*_flatfield\.mat$)", "the phase flat field (basic extra)"),
]


def tree_spec():
    """how the outputs are compared: parameters.mat without its time stamps and MATLAB's
    inputParser object; the per-pair registration results with their correlation peaks within
    P13 (the phase stitching's 2D registration, normxcorr2) and P7 (3D), integer shifts exact"""
    from petakit5d.testing.compare import tol, tree

    return tree(
        per_file={
            "*/Phase/*/xcorr/*/xcorr_tile_*.mat": tol(
                reason="P13 (proposed): 2D registration peaks of normxcorr2", abs=1e-7
            ),
            "*/xcorr/*/xcorr_tile_*.mat": tol(
                reason="P7 (approved): 3D registration peaks of the default FFT backend", abs=5e-3
            ),
        },
        ignore_vars=("t", "t1", "ip"),
    )


def demo_of(rel: str, default: str = "other") -> str:
    for pat, demo in RULES:
        if re.match(pat, rel):
            return demo
    return default


def not_ported(rel: str, python_rels: frozenset[str] = frozenset()) -> str:
    for pat, step in NOT_PORTED:
        if re.match(pat, rel):
            return step
    for pat, step in OPTIONAL:
        if re.match(pat, rel) and rel not in python_rels:
            return step + ", not installed"
    return ""


def main(argv: list[str]) -> int:
    from petakit5d.testing.compare import Masks, compare_trees, list_tree, scan_tree

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("matlab_root")
    ap.add_argument("python_root")
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--matlab-home", default="")
    ap.add_argument("--python-home", default="")
    ap.add_argument("--matlab-code", default="")
    ap.add_argument("--python-code", default="")
    ap.add_argument("--only", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--label", default="demos")
    ap.add_argument("--default-demo", default="other")
    ap.add_argument("--demo", default="")
    ap.add_argument("--workers", type=int, default=8)
    ns = ap.parse_args(argv)

    t0 = time.time()
    datasets = [d for d in (CELL, OTHER) if os.path.isdir(os.path.join(ns.python_root, d))]
    inputs = {
        f"{d}/{r}" for d in datasets for r in list_tree(os.path.join(ns.inputs, d))
    }  # fmt: skip

    def outputs(root: str) -> list[str]:
        # files and zarr arrays; empty folders are left out (a run's flag folders, and the
        # MATLAB harness keeps no empty folders)
        rels = [r for r in list_tree(root) if r not in inputs]
        return [r for r in rels if r.endswith(".zarr") or not os.path.isdir(os.path.join(root, r))]

    mrels = outputs(ns.matlab_root)
    prels = outputs(ns.python_root)
    if ns.only:
        mrels = [r for r in mrels if re.search(ns.only, r)]
        prels = [r for r in prels if re.search(ns.only, r)]
    # the MATLAB demos pass '~/Downloads/...' on: '~' is the home folder
    tilde = (r"(?<![\w.~-])~(?=/)", "<WORK>")
    mm = Masks().with_paths(work=ns.matlab_home, code=ns.matlab_code).plus(tilde)
    pm = Masks().with_paths(work=ns.python_home, code=ns.python_code).plus(tilde)
    spec = tree_spec()

    pset = frozenset(prels)
    skipped = {r: s for r in mrels if (s := not_ported(r, pset))}
    by_demo: dict[str, tuple[list[str], list[str]]] = {}
    for r in mrels:
        if r not in skipped:
            by_demo.setdefault((ns.demo or demo_of(r, ns.default_demo)), ([], []))[0].append(r)
    for r in prels:
        by_demo.setdefault((ns.demo or demo_of(r, ns.default_demo)), ([], []))[1].append(r)

    def one(item):
        demo, (mr, pr) = item
        ref, act = scan_tree(ns.matlab_root, mr), scan_tree(ns.python_root, pr)
        return demo, len(mr), len(pr), compare_trees(ref, act, spec, ref_masks=mm, act_masks=pm)

    lines = [f"compare {ns.label}: MATLAB {ns.matlab_root}, Python {ns.python_root}"]
    result = {"label": ns.label, "demos": {}, "not_ported": {}}
    bad = 0
    with ThreadPoolExecutor(ns.workers) as ex:
        for demo, nm, npy, rep in ex.map(one, sorted(by_demo.items())):
            lines.append(rep.text(f"\n== {demo} ({nm} MATLAB, {npy} Python outputs)"))
            result["demos"][demo] = {"matlab": nm, "python": npy, **rep.to_json()}
            bad += len(rep.problems)
    steps: dict[str, list[str]] = {}
    for r, step in skipped.items():
        steps.setdefault(step, []).append(r)
    lines.append("\n== MATLAB outputs of steps not in Phase 1 (not compared)")
    for step, rs in sorted(steps.items()):
        lines.append(f"  {step}: {len(rs)}, e.g. {rs[0]}")
        result["not_ported"][step] = rs
    lines.append(f"\ncompare {ns.label}: {bad} problems ({time.time() - t0:.0f} s)")
    text = "\n".join(lines)
    print(text)
    if ns.out:
        os.makedirs(ns.out, exist_ok=True)
        with open(os.path.join(ns.out, f"compare_{ns.label}.txt"), "w", newline="\n") as f:
            f.write(text + "\n")
        with open(os.path.join(ns.out, f"compare_{ns.label}.json"), "w", newline="\n") as f:
            json.dump(result, f, indent=1)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
