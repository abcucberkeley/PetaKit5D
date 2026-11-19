"""Compares the outputs of the demos' Phase 2 steps with those of the MATLAB demos, step by step.

    python compare_phase2_trees.py MATLAB_ROOT PYTHON_ROOT --inputs INPUTS [options]

``MATLAB_ROOT`` and ``PYTHON_ROOT`` hold the dataset folders after the runs (the demos'
``Downloads``); ``INPUTS`` holds untouched datasets, whose files are not outputs (repeatable:
also a folder of files a run started from). Every output of a Phase 2 step (``STEPS``, first
match) is compared with ``petakit5d.testing.compare``: tiff voxels and tags, .mat values, JSON
and text; figures (.png, .pdf) by name only, their pixels being matplotlib's. The phase 1
steps' outputs are left to ``examples/phase1/compare_demo_trees.py``.

``--backend`` is the FFT backend of the Python run and picks the tolerances of
``docs/parity/index.md``: with ``exact`` (the FFTW add-on) every output must be MATLAB's except
where a row holds for both backends (N3: the FSC resolutions; N8: the phase flat field of
BaSiCPy); ``default`` adds P10 (PSF analysis) and, for the FFT spectra and the FSC curves, which
are made from the decon results of the same run, P14 (P2's relative L2 of what is made from decon
results). With ``--matlab-inputs`` (the steps ran on MATLAB's own decon results) those two take
P9 and P11 instead, the steps' own rows. Outputs of steps that are not ported yet (``PENDING``)
and MATLAB-only files are listed, not compared. For every compared array that is not MATLAB's
bit for bit, the numbers are reported: the values that differ, the largest difference
(absolute, relative to MATLAB's largest magnitude, in ulp), the relative differences' 99th
percentile and median, and the relative L2 difference (``docs/parity/phase2.md``).

Options: ``--matlab-home``, ``--python-home``, ``--matlab-code``, ``--python-code`` (the folders
the runs' paths start with, masked; a leading ``~`` is the home folder), ``--only REGEX``,
``--out DIR`` (``compare_<label>.txt`` and ``.json``), ``--label``, ``--workers``. Exit status
1 if an output is missing, extra, or differs beyond its tolerance.
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

# the outputs of the Phase 2 steps: (regex on the name relative to the root, demo, step); first
# match wins. Figures the phase 1 steps draw with the figures extra (the decon wrappers' back
# projectors) are checked by name, with the rest of the figures.
STEPS = [
    (_C + r"matlab_decon_omw/DS_for_useful_tools/FFT/", "useful_tools", "FFT analysis"),
    (_C + r"PSF/DS/", "useful_tools", "PSF analysis"),
    (_C + r"matlab_decon_omw_debug/[^/]+_debug/(FSCs|figures)/", "useful_tools", "FSC analysis"),
    (_C + r"Chromatic_Shift_Corrected/", "useful_tools", "chromatic shift correction"),
    (_O + r"Phase/ff_corrected_1/", "phase_and_2photon_stitching", "phase flat field"),
    (r".*\.png$", "", "figures of the phase 1 steps"),
]

# MATLAB outputs that are not compared: (regex, reason)
PENDING = [
    (
        _C + r"matlab_decon_omw/DSR_for_useful_tools/imaris/",
        "the Imaris converter: XR_imaris_conversion_data_wrapper is not ported yet (it waits for "
        "the next release of Parallel_Imaris_Writer)",
    ),
    (r".*\.fig$", "MATLAB figure files: the port writes .png only (parity row F6)"),
]


def step_of(rel: str) -> tuple[str, str] | None:
    """(demo, step) of an output, or None for the outputs of phase 1 steps"""
    for pat, demo, step in STEPS:
        if re.match(pat, rel):
            return demo, step
    return None


def pending(rel: str) -> str:
    """why an output is not compared, or '' if it is"""
    for pat, why in PENDING:
        if re.match(pat, rel):
            return why
    return ""


# ---------------------------------------------------------------------------------------------
# tolerances (docs/parity/index.md)
# ---------------------------------------------------------------------------------------------


def _floats(v):
    """an MArray's values as a double vector in MATLAB's order"""
    import numpy as np

    from petakit5d.testing.compare import MArray

    if not isinstance(v, MArray) or v.data is None:
        raise TypeError("not an array with data")
    return np.asarray(v.data, dtype=np.float64).ravel(order="F")


def _power_close(power):
    """record 0006's single-precision tolerance on the normalized spectrum |F| / max: largest
    difference <= 1e-5, relative L2 <= 1e-6; ``power`` maps an image's values to it"""

    def fn(ref, act) -> list[str]:
        import numpy as np

        try:
            r, a = power(_floats(ref)), power(_floats(act))
        except TypeError as e:
            return [str(e)]
        if r.shape != a.shape:
            return [f"{r.size} values against {a.size}"]
        d = float(np.max(np.abs(r - a))) if r.size else 0.0
        l2 = float(np.linalg.norm(r - a) / max(np.linalg.norm(r), 1e-300)) if r.size else 0.0
        out = []
        if d > 1e-5:
            out.append(f"largest difference of |F|/max {d:g} > 1e-5")
        if l2 > 1e-6:
            out.append(f"relative L2 of |F|/max {l2:g} > 1e-6")
        return out

    return fn


def _res_problems(ref, act, steps: float = 1, rel: float = 0.01) -> list[str]:
    """resolution structs (or cells of them): comp_angles and thetas exact, every radius within
    ``steps`` steps of 1e-4, the resolutions within ``rel`` (NaN and Inf at the same places)"""
    import numpy as np

    from petakit5d.core.dtypes import mround
    from petakit5d.testing.compare import MCell, MStruct

    def structs(v, path):
        if isinstance(v, MCell):
            return [
                s for k, x in enumerate(v.items or []) for s in structs(x, f"{path}{{{k + 1}}}")
            ]
        if isinstance(v, MStruct):
            return [(path, v.items[0])]
        raise TypeError(f"{path}: not a struct or cell")

    try:
        rs, acts = structs(ref, ""), structs(act, "")
    except TypeError as e:
        return [str(e)]
    if len(rs) != len(acts):
        return [f"{len(rs)} structs in MATLAB's, {len(acts)} in Python's"]
    out = []
    for (path, r), (_, a) in zip(rs, acts, strict=True):
        if sorted(r) != sorted(a):
            out.append(f"{path}: fields {sorted(r)} against {sorted(a)}")
            continue
        for f in ("comp_angles", "thetas"):
            if _floats(r[f]).tobytes() != _floats(a[f]).tobytes():
                out.append(f"{path}.{f} differs")
        pr, pa = _floats(r["points"]), _floats(a["points"])
        if pr.shape != pa.shape:
            out.append(f"{path}.points: {pr.size} values against {pa.size}")
            continue
        if not np.array_equal(np.isnan(pr), np.isnan(pa)):
            out.append(f"{path}.points: NaN at other places")
        d = mround(np.abs(pr - pa) / 1e-4)  # steps of the grid
        if np.nanmax(d, initial=0) > steps:
            out.append(f"{path}.points: up to {np.nanmax(d):.0f} steps of 1e-4 (limit {steps:g})")
        out += _rel_problems(_floats(r["res"]), _floats(a["res"]), rel, f"{path}.res")
    return out


def _rel_problems(r, a, rel: float, path: str) -> list[str]:
    import numpy as np

    if r.shape != a.shape:
        return [f"{path}: {r.size} values against {a.size}"]
    if not (np.array_equal(np.isnan(r), np.isnan(a)) and np.array_equal(np.isinf(r), np.isinf(a))):
        return [f"{path}: Inf or NaN at other places"]
    fin = np.isfinite(r)
    if not np.array_equal(r[np.isinf(r)], a[np.isinf(a)]):
        return [f"{path}: other infinities"]
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.abs(r[fin] - a[fin]) / np.abs(r[fin])
    if d.size and np.nanmax(d) > rel:
        return [f"{path}: relative difference up to {np.nanmax(d):.3g} (limit {rel:g})"]
    return []


def per_variable(specs: dict, name: str, default=None):
    """a .mat file compared variable by variable, each with its spec (``default``: exact)"""
    from petakit5d.testing.compare import EXACT, Custom, MStruct, compare

    def fn(ref, act) -> list[str]:
        if not (isinstance(ref, MStruct) and isinstance(act, MStruct)):
            return ["not the variables of a MAT-file"]
        rv, av = ref.items[0], act.items[0]
        if sorted(rv) != sorted(av):
            return [f"variables {sorted(rv)} against {sorted(av)}"]
        out = []
        for k in sorted(rv):
            rep = compare(rv[k], av[k], specs.get(k, default or EXACT), path=k)
            out += [f"{p.path}: {p.what}: {'; '.join(p.lines[:2])}" for p in rep.problems]
        return out

    return Custom(fn, name)


def tolerances(backend: str, matlab_inputs: bool = False) -> dict:
    """the specs of the parity rows (docs/parity/index.md), for ``backend``; ``matlab_inputs``:
    the FFT and FSC steps ran on MATLAB's decon results, so their own rows (P9, P11) apply"""
    import numpy as np

    from petakit5d.testing.compare import EXACT, Custom, tol

    default = backend == "default"
    res = Custom(_res_problems, "N3 (approved): the FSC radii within a step of 1e-4, res 1%")
    res_values = Custom(
        lambda r, a: _rel_problems(_floats(r), _floats(a), 0.01, "res"),
        "N3 (approved): the FSC resolutions within 1%",
    )
    double_fft = tol(
        reason="P10/P11 (approved): double FFT results of the default backend, record 0006",
        rel_max=1e-13,
        rel_l2=1e-14,
    )
    # P14 (approved): the spectra and FSC curves of the default backend's decon results carry
    # the decon's differences (P2); P2's relative L2 for what is made from decon results
    from_decon = tol(
        reason="P14 (approved): made from the default backend's decon results, P2's relative L2",
        rel_l2=1e-4,
    )
    gamma = Custom(
        _power_close(np.square),
        "P9 (approved): the spectra's MIPs and slices as |F|/max within 0006's single tolerance",
    )
    curves = (double_fft if matlab_inputs else from_decon) if default else EXACT
    basic = "N8 (approved): BaSiCPy's flat field against MATLAB's BaSiC"
    return {
        "gamma": (gamma if matlab_inputs else from_decon) if default else EXACT,
        "double_fft": double_fft if default else EXACT,
        "fsc_file": per_variable(
            {"fsc_mu": curves, "fsc_cell": curves, "res_mu": res, "res_cell": res},
            "FSC results: curves "
            + (("P11" if matlab_inputs else "P14") if default else "exact")
            + ", resolutions N3",
        ),
        "fsc_info": per_variable(
            {"res_mat": res_values, "res_mu_mat": res_values}, "FSC resolution table: N3"
        ),
        "tiles": tol(reason=basic, abs=1, rel=0.06, rel_l2=0.01),
        "flatfield": per_variable(
            {
                "flatfield": tol(reason=basic, rel=0.06, rel_l2=0.01),
                "minval": tol(reason="matlabcompat: imgaussfilt frequency domain", abs=0.01),
            },
            basic,
        ),
    }


def step_spec(step: str, backend: str, matlab_inputs: bool = False):
    """how the outputs of ``step`` are compared"""
    from petakit5d.testing.compare import ignore, tree

    t = tolerances(backend, matlab_inputs)
    per_file: dict = {
        "*.png": ignore("figures: matplotlib's pixels"),
        "*.pdf": ignore("figures: matplotlib's"),
    }
    tiff_bytes = True
    if step == "FFT analysis" and backend == "default":
        per_file.update({"*_MIP_*.tif": t["gamma"], "*_center_*.tif": t["gamma"]})
        tiff_bytes = False
    elif step == "PSF analysis":
        per_file["*/PSFAnalysis/*.mat"] = t["double_fft"]
    elif step == "FSC analysis":
        per_file["*/FSCs/Iter_*.mat"] = t["fsc_file"]
        per_file["*/fsc_res_info_c*.mat"] = t["fsc_info"]
    elif step == "phase flat field":
        per_file.update({"*.tif": t["tiles"], "*_flatfield.mat": t["flatfield"]})
        tiff_bytes = False
    return tree(per_file=per_file, ignore_vars=("t", "t1", "ip"), tiff_bytes=tiff_bytes)


def strict_spec():
    """everything exact (figures by name): which outputs are MATLAB's bit for bit"""
    from petakit5d.testing.compare import ignore, tree

    return tree(
        per_file={"*.png": ignore("figures"), "*.pdf": ignore("figures")},
        ignore_vars=("t", "t1", "ip"),
        tiff_bytes=False,
    )


# ---------------------------------------------------------------------------------------------
# the numbers
# ---------------------------------------------------------------------------------------------


def array_stats(r, a) -> dict:
    """how two arrays of the same size differ: count, largest absolute difference, the same
    relative to the reference's largest magnitude and in ulp (floats), relative L2"""
    import numpy as np

    r, a = np.asarray(r), np.asarray(a)
    out: dict = {"n": int(r.size), "class": str(r.dtype)}
    if r.shape != a.shape:
        out["size"] = f"{r.shape} against {a.shape}"
        return out
    rf, af = r.astype(np.float64).ravel(), a.astype(np.float64).ravel()
    both_nan = np.isnan(rf) & np.isnan(af)
    diff = ~(both_nan | (rf == af))
    out["differ"] = int(np.count_nonzero(diff))
    if not out["differ"]:
        return out
    fin = np.isfinite(rf) & np.isfinite(af)
    d = np.abs(rf - af)[fin & diff]
    big = float(np.max(np.abs(rf[np.isfinite(rf)]), initial=0.0))
    out["nonfinite_differ"] = int(np.count_nonzero(diff & ~fin))
    out["max_abs"] = float(d.max(initial=0.0))
    out["max_rel_to_max"] = out["max_abs"] / big if big else None
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = d / np.abs(rf[fin & diff])
    out["max_rel"] = float(np.nanmax(rel[np.isfinite(rel)], initial=0.0)) if rel.size else 0.0
    if r.size >= 100:
        # the relative differences of all values (0 where equal): median and 99th percentile
        nz = fin & (rf != 0)
        allrel = np.abs(rf[nz] - af[nz]) / np.abs(rf[nz])
        for q, key in ((0.5, "rel_p50"), (0.99, "rel_p99")):
            k = int(q * (allrel.size - 1))
            out[key] = float(np.partition(allrel, k)[k]) if allrel.size else None
    nr = float(np.linalg.norm(rf[fin]))
    out["rel_l2"] = float(np.linalg.norm((rf - af)[fin]) / nr) if nr else None
    if r.dtype.kind == "f" and a.dtype == r.dtype:
        it = np.int32 if r.dtype.itemsize == 4 else np.int64

        def ordered(x):
            i = np.ascontiguousarray(x).ravel().view(it).astype(np.int64)
            return np.where(i < 0, np.iinfo(it).min - i, i)

        u = np.abs(ordered(r)[fin & diff] - ordered(a)[fin & diff])
        out["max_ulp"] = int(u.max(initial=0))
    return out


def value_stats(r, a, path: str = "") -> list[tuple[str, dict]]:
    """``array_stats`` of the differing arrays of two MATLAB values (structs and cells walked)"""
    from petakit5d.testing.compare import MArray, MCell, MStruct

    if isinstance(r, MStruct) and isinstance(a, MStruct) and len(r.items) == len(a.items):
        out = []
        for i, (ri, ai) in enumerate(zip(r.items, a.items, strict=True)):
            sub = f"({i + 1})" if len(r.items) > 1 else ""
            for k in ri:
                if k in ai:
                    out += value_stats(ri[k], ai[k], f"{path}{sub}.{k}" if path or sub else k)
        return out
    if isinstance(r, MCell) and isinstance(a, MCell) and len(r.items or []) == len(a.items or []):
        out = []
        for i, (ri, ai) in enumerate(zip(r.items or [], a.items or [], strict=True)):
            out += value_stats(ri, ai, f"{path}{{{i + 1}}}")
        return out
    if isinstance(r, MArray) and isinstance(a, MArray) and r.data is not None:
        if a.data is None or r.mclass == "char":
            return []
        s = array_stats(r.data, a.data)
        return [(path, s)] if s.get("differ") or "size" in s else []
    return []


def file_stats(mroot: str, mrel: str, proot: str, prel: str) -> list[tuple[str, dict]]:
    """the numbers of one output pair (tiff voxels, .mat variables)"""
    from petakit5d.testing.compare import read_mat

    m, p = os.path.join(mroot, mrel), os.path.join(proot, prel)
    low = mrel.lower()
    if low.endswith((".tif", ".tiff")) and os.path.getsize(m) and os.path.getsize(p):
        from petakit5d.io import read_tiff

        s = array_stats(read_tiff(m), read_tiff(p))
        return [("voxels", s)] if s.get("differ") or "size" in s else []
    if low.endswith(".mat"):
        try:
            return value_stats(read_mat(m), read_mat(p))
        except Exception as e:  # an unreadable file is compared as bytes by the tree compare
            return [("file", {"error": f"{type(e).__name__}: {e}"})]
    return []


def summary_text(numbers: dict) -> str:
    """the numbers of a step's differing arrays together: counts and the largest of each"""
    stats = [s for v in numbers.values() for _, s in v if "differ" in s]
    if not stats:
        return ""
    n, d = sum(s["n"] for s in stats), sum(s["differ"] for s in stats)
    parts = [f"{d} of {n} values of {len(stats)} arrays differ"]
    for key, label in (
        ("max_abs", "largest difference"),
        ("max_rel_to_max", "largest relative to the largest value"),
        ("max_rel", "largest relative"),
        ("rel_p99", "99th percentile relative"),
        ("rel_l2", "relative L2"),
        ("max_ulp", "ulp"),
    ):
        xs = [s[key] for s in stats if s.get(key) is not None]
        if xs:
            parts.append(f"{label} {max(xs):.3g}")
    return "  over the step: " + ", ".join(parts)


def stats_text(s: dict) -> str:
    if "error" in s:
        return s["error"]
    if "size" in s:
        return f"size {s['size']}"
    parts = [f"{s['differ']} of {s['n']} {s['class']} values differ"]
    if s.get("max_abs") is not None:
        parts.append(f"largest by {s['max_abs']:.3g}")
    if s.get("max_rel_to_max") is not None:
        parts.append(f"{s['max_rel_to_max']:.3g} of the largest value")
    if s.get("max_ulp") is not None:
        parts.append(f"{s['max_ulp']} ulp")
    if s.get("max_rel") is not None:
        parts.append(f"relative {s['max_rel']:.3g}")
    if s.get("rel_p99") is not None:
        parts.append(f"(99th percentile {s['rel_p99']:.3g}, median {s['rel_p50']:.3g})")
    if s.get("rel_l2") is not None:
        parts.append(f"relative L2 {s['rel_l2']:.3g}")
    if s.get("nonfinite_differ"):
        parts.append(f"{s['nonfinite_differ']} at Inf or NaN")
    return ", ".join(parts)


# ---------------------------------------------------------------------------------------------
# the comparison
# ---------------------------------------------------------------------------------------------


def outputs(root: str, inputs: set[str]) -> list[str]:
    """the files and zarr arrays under ``root`` that are not inputs (empty folders left out)"""
    from petakit5d.testing.compare import list_tree

    rels = [r for r in list_tree(root) if r not in inputs]
    return [r for r in rels if r.endswith(".zarr") or not os.path.isdir(os.path.join(root, r))]


def name_of(path: str, names) -> str | None:
    """the output a problem's path is about (``files/<name>``, then ``.field``, ``(i)``,
    ``{i}``, `` (voxels)`` or ``/<zarr part>``)"""
    for n in sorted(names, key=len, reverse=True):
        head = "files/" + n
        if path == head or (path.startswith(head) and path[len(head)] in " .({/"):
            return n
    return None


def compare_step(ns, step, mrels, prels, mm, pm) -> dict:
    """one step's outputs: counts, the identical ones, problems beyond the tolerances, and the
    numbers of those that differ"""
    from petakit5d.testing.compare import compare_trees, normalize_name, scan_tree

    ref, act = scan_tree(ns.matlab_root, mrels), scan_tree(ns.python_root, prels)
    spec = step_spec(step, ns.backend, ns.matlab_inputs)
    rep = compare_trees(ref, act, spec, ref_masks=mm, act_masks=pm)
    strict = compare_trees(ref, act, strict_spec(), ref_masks=mm, act_masks=pm)
    mnames = {normalize_name(r, mm): r for r in mrels}
    pnames = {normalize_name(r, pm): r for r in prels}
    both = set(mnames) & set(pnames)
    figs = {n for n in both if n.lower().endswith((".png", ".pdf"))}
    differing = sorted(
        {n for p in strict.problems if (n := name_of(p.path, both)) is not None} - figs
    )
    numbers = {
        n: file_stats(ns.matlab_root, mnames[n], ns.python_root, pnames[n]) for n in differing
    }
    return {
        "matlab": len(mrels),
        "python": len(prels),
        "figures": len(figs),
        "identical": len(both - figs) - len(differing),
        "differing": differing,
        "numbers": numbers,
        "report": rep,
    }


def main(argv: list[str]) -> int:
    from petakit5d.testing.compare import Masks, list_tree

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("matlab_root")
    ap.add_argument("python_root")
    ap.add_argument("--inputs", action="append", required=True)
    ap.add_argument("--backend", choices=("exact", "default"), default="exact")
    ap.add_argument("--matlab-inputs", action="store_true")
    ap.add_argument("--matlab-home", default="")
    ap.add_argument("--python-home", default="")
    ap.add_argument("--matlab-code", default="")
    ap.add_argument("--python-code", default="")
    ap.add_argument("--only", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--label", default="phase2")
    ap.add_argument("--workers", type=int, default=6)
    ns = ap.parse_args(argv)

    t0 = time.time()
    inputs: set[str] = set()
    for root in ns.inputs:
        for d in (CELL, OTHER):
            if os.path.isdir(os.path.join(root, d)):
                inputs |= {f"{d}/{r}" for r in list_tree(os.path.join(root, d))}
    mrels, prels = outputs(ns.matlab_root, inputs), outputs(ns.python_root, inputs)
    if ns.only:
        mrels = [r for r in mrels if re.search(ns.only, r)]
        prels = [r for r in prels if re.search(ns.only, r)]
    tilde = (r"(?<![\w.~-])~(?=/)", "<WORK>")
    mm = Masks().with_paths(work=ns.matlab_home, code=ns.matlab_code).plus(tilde)
    pm = Masks().with_paths(work=ns.python_home, code=ns.python_code).plus(tilde)

    skipped: dict[str, list[str]] = {}
    unexpected = [r for r in prels if pending(r)]
    by_step: dict[tuple[str, str], tuple[list[str], list[str]]] = {}
    for side, rels in ((0, mrels), (1, prels)):
        for r in rels:
            why = pending(r)
            if why:
                if side == 0:
                    skipped.setdefault(why, []).append(r)
                continue
            key = step_of(r)
            if key is not None:
                by_step.setdefault(key, ([], []))[side].append(r)
    others = sum(1 for r in mrels if step_of(r) is None and not pending(r))

    with ThreadPoolExecutor(ns.workers) as ex:
        jobs = {k: ex.submit(compare_step, ns, k[1], *v, mm, pm) for k, v in by_step.items()}
        results = {k: jobs[k].result() for k in sorted(jobs)}

    on = ", the steps on MATLAB's decon results" if ns.matlab_inputs else ""
    lines = [f"compare {ns.label} ({ns.backend} FFT backend{on}): MATLAB {ns.matlab_root}, "
             f"Python {ns.python_root}"]  # fmt: skip
    out = {"label": ns.label, "backend": ns.backend, "matlab_inputs": ns.matlab_inputs}
    out |= {"steps": {}, "not_compared": skipped}
    bad = 0
    for (demo, step), res in results.items():
        rep = res.pop("report")
        title = f"{demo}: {step}" if demo else step
        lines.append(
            rep.text(
                f"\n== {title}: {res['matlab']} MATLAB, {res['python']} Python outputs; "
                f"{res['identical']} identical, {len(res['differing'])} differ"
                + (f", {res['figures']} figures by name" if res["figures"] else "")
            )
        )
        if res["numbers"]:
            lines.append(summary_text(res["numbers"]))
        for n, numbers in res["numbers"].items():
            for k, s in numbers:
                lines.append(f"    {n} [{k}]: {stats_text(s)}")
        out["steps"][title] = {**res, **rep.to_json()}
        bad += len(rep.problems)
    lines.append("\n== MATLAB outputs not compared")
    for why, rs in sorted(skipped.items()):
        lines.append(f"  {why}: {len(rs)}, e.g. {rs[0]}")
    lines.append(f"  outputs of phase 1 steps (examples/phase1/compare_demo_trees.py): {others}")
    if unexpected:
        lines.append(f"\n== Python wrote outputs of steps listed as not ported: {len(unexpected)}")
        lines.extend(f"  {r}" for r in unexpected[:20])
        out["unexpected"] = unexpected
        bad += len(unexpected)
    lines.append(f"\ncompare {ns.label}: {bad} problems ({time.time() - t0:.0f} s)")
    text = "\n".join(lines)
    print(text)
    if ns.out:
        os.makedirs(ns.out, exist_ok=True)
        with open(os.path.join(ns.out, f"compare_{ns.label}.txt"), "w", newline="\n") as f:
            f.write(text + "\n")
        with open(os.path.join(ns.out, f"compare_{ns.label}.json"), "w", newline="\n") as f:
            json.dump(out, f, indent=1, default=str)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
