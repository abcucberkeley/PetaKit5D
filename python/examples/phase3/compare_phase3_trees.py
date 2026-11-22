"""Compares the outputs of a point-detection demo's run with those of the MATLAB demo, step by
step.

    python compare_phase3_trees.py MATLAB_ROOT PYTHON_ROOT --inputs INPUTS [options]

``MATLAB_ROOT`` and ``PYTHON_ROOT`` hold the dataset folder (``PetaKit5D_point_detection_demo_
dataset``) after the runs (each run's ``PETAKIT5D_DEMO_DATA``); ``INPUTS`` holds the untouched
dataset, whose files are not outputs (repeatable). ``--matlab-saved`` and ``--python-saved``
are the folders of the runs' in-memory results (``PETAKIT5D_DEMO_SAVE_DIR``: the condition, the
PSF sigmas, the lifetime results), compared as their steps' outputs. Every output is assigned
to its step (``STEPS``, first match) and compared with ``petakit5d.testing.compare``: tiff
voxels and tags, .mat values, text; the Amira files (``*.am``) by their numbers.

The tolerances are the parity rows of ``docs/parity/detection.md``, for the step's own
differences, and end to end the row E1 (``docs/parity/phase3.md``) for those it
inherits from the steps before it: the PSF sigmas (D6) feed the detection, whose fits (D3) feed
the tracks' coordinates, the track processing (whose own fits are W11 and W12), the rotation
and the lifetime analysis (whose threshold is L4, the cohorts' amplitudes L7; its files are
compared field by field). With ``--matlab-inputs`` (each step of the Python run made on
MATLAB's results of the steps before it) only each step's own rows apply.

For every compared array that is not MATLAB's bit for bit the numbers are reported: the values
that differ, the largest difference (absolute, relative), its 99th percentile and median, and
the relative L2 difference, per field of the detection and track files (their frames and
tracks together); and the detections per frame and the tracks that differ in number
(``docs/parity/phase3.md``).

The runs' paths are masked: ``--matlab-home`` and ``--python-home`` (the runs'
``PETAKIT5D_DEMO_DATA``). The condition's ``date`` is the six-digit runs of its folders' path,
so compare runs whose roots hold the same digits (or none). Options: ``--only REGEX``, ``--out
DIR`` (``compare_<label>.txt`` and ``.json``), ``--label``, ``--workers``. Exit status 1 if an
output is missing, extra, or differs beyond its tolerance.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

DATASET = "PetaKit5D_point_detection_demo_dataset"
SAVED = "saved/"
_E = re.escape(DATASET) + r"/[^/]+/Ex\d+[^/]*/"  # an experiment of the condition
_A = _E + r"ch1/(?!DS/)[^/]+/"  # its results folder (aname), under the master channel

# the outputs of each step: (regex on the name relative to the root, step); first match wins
STEPS = [
    (SAVED + r"[^/]+/load\.mat$", "load"),
    (SAVED + r"[^/]+/sigma\.mat$", "PSF sigmas"),
    (_E + r"ch\d+/DS/xz_offset\.mat$", "XZ offsets"),
    (_E + r"ch\d+/DS/[^/]+\.tif$", "deskew"),
    (_A + r"(settings|dynRange)\.mat$", "deskew"),
    (_A + r"(Detection3D\.mat|Detections/|Masks/)", "detection"),
    (_A + r"trackedFeatures\.mat$", "tracking"),
    (_A + r"(ProcessedTracks\.mat|Amira/Tracking/)", "track processing"),
    (_A + r"RotatedTracks\.mat$", "track rotation"),
    (_A + r"(lifetimeData|Parameters)\.mat$", "lifetime analysis"),
    (SAVED + r"[^/]+/lifetime\.mat$", "lifetime analysis"),
    (SAVED + r"[^/]+/cohorts\.mat$", "intensity cohorts"),
]
ORDER = [
    "load",
    "PSF sigmas",
    "deskew",
    "XZ offsets",
    "detection",
    "tracking",
    "track processing",
    "track rotation",
    "lifetime analysis",
    "intensity cohorts",
]


def step_of(rel: str) -> str | None:
    for pat, step in STEPS:
        if re.match(pat, rel):
            return step
    return None


# ---------------------------------------------------------------------------------------------
# tolerances (docs/parity/detection.md)
# ---------------------------------------------------------------------------------------------

# D6: GU_estimateSigma3D's sigmas (the fitGaussian3D mex's nonlinear fits, D5)
SIGMA_REL = 1e-7
# D3: the detector's fits on the oracle cases; W11 and W12: the track processing's gap and
# buffer fits (TRACKS_TOL of tests/detection/test_processing.py)
FIT_REL, FIT_ABS = 1e-9, 1e-9
TRACK_FIT_REL, TRACK_FIT_ABS = 1e-11, 1e-10
# E1: end to end, every fit of the detection and of the track processing starts from the
# PSF sigmas of each side (D6), and carries their differences
E2E_REL, E2E_ABS = 1e-8, 1e-8
# L4: the lifetime analysis' intensity threshold (its Gaussian-mode fit); L7: the intensity
# cohorts' amplitudes (scaled by L2's factors)
THRESHOLD_REL = 1e-4
COHORT_REL, COHORT_ABS = 5e-3, 1e-9
# the lifetime analysis' and the cohorts' files, compared field by field (lifetime_problems)
LIFETIME_FILES = ("*/lifetimeData.mat", "*/Parameters.mat", "*/lifetime.mat", "*/cohorts.mat")


def tolerances(matlab_inputs: bool = False) -> dict:
    """{glob: spec} of the outputs that may differ from MATLAB's (every other output exact):
    the steps' own rows with ``matlab_inputs``, else E1 for what the sigmas feed"""
    from petakit5d.testing.compare import tol

    out = {"*/sigma.mat": tol(reason="D6 (approved): GU_estimateSigma3D's sigmas", rel=SIGMA_REL)}
    if matlab_inputs:
        fits = tol(reason="D3 (approved): the detector's fits", rel=FIT_REL, abs=FIT_ABS)
        # the tracker on MATLAB's detections: TK1-TK4 (approved), none on the outputs
        track_fits = tol(
            reason="W11, W12 (approved): the gap and buffer fits and their p-values",
            rel=TRACK_FIT_REL,
            abs=TRACK_FIT_ABS,
        )
        out |= {"*/Detection3D.mat": fits, "*/Detections/*.mat": fits}
    else:
        track_fits = tol(
            reason="E1 (approved): the demos end to end, the fits from each side's PSF sigmas "
            "(D6) and what is made from them",
            rel=E2E_REL,
            abs=E2E_ABS,
        )
        out |= {"*/Detection3D.mat": track_fits, "*/Detections/*.mat": track_fits}
        out["*/trackedFeatures.mat"] = track_fits
    out |= {"*/ProcessedTracks.mat": track_fits, "*/RotatedTracks.mat": track_fits}
    out["*.am"] = track_fits
    return out


# path texts compare separator-insensitively: the port writes Windows' filesep there, as MATLAB
# does, and the references are made on Linux
SLASHES = ((r"\\", "/"),)


def tree_spec(matlab_inputs: bool = False, ignore: tuple = (), amira=None):
    """how the output files are compared: ``tolerances``, path texts with slashes, the Amira
    files by ``amira_problems`` (``amira``: their spec instead, for a tree comparison of its
    own)"""
    from petakit5d.testing.compare import ignore as not_compared
    from petakit5d.testing.compare import tree

    per_file = dict(tolerances(matlab_inputs))
    per_file["*.am"] = amira or not_compared("Amira files: by their numbers (amira_problems)")
    for g in LIFETIME_FILES:
        per_file[g] = not_compared("the lifetime files: field by field (lifetime_problems)")
    return tree(per_file=per_file, ignore=ignore, ignore_vars=("t", "t1", "ip"), masks=SLASHES)


def lifetime_specs(matlab_inputs: bool = False) -> tuple[dict, object]:
    """the lifetime files' fields that may differ, ``{variable or variable.field: spec}``, and
    the spec of the others: on MATLAB's tracks the threshold (L4) and the cohorts' amplitudes
    (L7), the others exact; end to end each also within E1, which the tracks' fits carry"""
    from petakit5d.testing.compare import EXACT, tol

    threshold = tol(reason="L4 (approved): the threshold's fit", rel=THRESHOLD_REL)
    cohort = tol(reason="L7 (approved): the scaled cohorts", rel=COHORT_REL, abs=COHORT_ABS)
    other = EXACT
    if not matlab_inputs:
        e1 = "; E1 (approved): the demos end to end"
        threshold = tol(reason=threshold.reason + e1, rel=THRESHOLD_REL, abs=E2E_ABS)
        cohort = tol(reason=cohort.reason + e1, rel=COHORT_REL, abs=max(COHORT_ABS, E2E_ABS))
        other = tol(
            reason="E1 (approved): the demos end to end, the lifetimes of each side's tracks",
            rel=E2E_REL,
            abs=E2E_ABS,
        )
    specs = {"lftRes.MaxIntensityThreshold": threshold, "MaxIntensityThreshold": threshold}
    specs |= {f"cohorts.{k}": cohort for k in ("A", "Amin", "Aplus")}
    return specs, other


def lifetime_problems(mv, pv, name: str, specs: dict, other, masks=None, act_masks=None):
    """a lifetime file (``mv``, ``pv``: ``testing.compare.read_mat``'s values) against MATLAB's,
    each variable by its spec, and the fields of a struct variable with fields in ``specs``
    (``variable.field``) each by theirs; the variables and fields in MATLAB's order"""
    from petakit5d.testing.compare import MStruct, Report, compare

    rep = Report()
    head = f"files/{name}"

    def names(where, r, a):
        if list(r) != list(a):
            only_m = [k for k in r if k not in a]
            only_p = [k for k in a if k not in r]
            what = f"MATLAB {list(r)}, Python {list(a)}"
            if only_m or only_p:
                what = f"MATLAB's only {only_m}, Python's only {only_p}"
            rep.add(where, "names", what)
        return [k for k in r if k in a]

    r, a = mv.items[0], pv.items[0]
    for var in names(head, r, a):
        rv, av, path = r[var], a[var], f"{head}.{var}"
        by_field = any(k.startswith(var + ".") for k in specs)
        if not (by_field and isinstance(rv, MStruct) and isinstance(av, MStruct)):
            rep.extend(
                compare(rv, av, specs.get(var, other), path, masks=masks, act_masks=act_masks)
            )
            continue
        if list(rv.size) != list(av.size):
            rep.add(path, "size", f"MATLAB {list(rv.size)}, Python {list(av.size)}")
            continue
        for i, (ri, ai) in enumerate(zip(rv.items, av.items, strict=True)):
            at = f"{path}({i + 1})" if len(rv.items) > 1 else path
            for k in names(at, ri, ai):
                spec = specs.get(f"{var}.{k}", other)
                rep.extend(
                    compare(ri[k], ai[k], spec, f"{at}.{k}", masks=masks, act_masks=act_masks)
                )
    return rep


# the note the comparison gives a value that differs within its tolerance
WITHIN = "values differ, within tolerance"


# ---------------------------------------------------------------------------------------------
# Amira files: text with numbers
# ---------------------------------------------------------------------------------------------

_NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$|^[-+]?(Inf|NaN)$")


def amira_problems(mtext: str, ptext: str, rel: float, abs_: float) -> tuple[list[str], dict]:
    """an Amira file against MATLAB's: the same tokens, the numbers within ``rel`` and
    ``abs_`` (integers exact); the problems and the numbers"""
    mt, pt = mtext.split(), ptext.split()
    if len(mt) != len(pt):
        return [f"{len(mt)} tokens in MATLAB's file, {len(pt)} in Python's"], {}
    n = differ = 0
    worst = 0.0
    problems = []
    for i, (a, b) in enumerate(zip(mt, pt, strict=True)):
        if a == b:
            if _NUMBER.match(a):
                n += 1
            continue
        if not (_NUMBER.match(a) and _NUMBER.match(b)):
            problems.append(f"token {i}: MATLAB {a!r}, Python {b!r}")
            continue
        n += 1
        differ += 1
        x, y = float(a), float(b)
        is_int = "." not in a + b and "e" not in (a + b).lower()
        d = abs(x - y)
        worst = max(worst, d / max(abs(x), 1.0))
        if is_int or d > max(abs_, rel * abs(x)):
            problems.append(f"token {i}: MATLAB {a}, Python {b} (by {d:.3g})")
    out = [f"{len(problems)} tokens beyond the tolerance; first: {problems[0]}"] if problems else []
    return out, {"n": n, "differ": differ, "max_rel_to_1": worst}


# ---------------------------------------------------------------------------------------------
# the numbers
# ---------------------------------------------------------------------------------------------


def array_stats(r, a) -> dict:
    """how two arrays of the same size differ: count, largest absolute and relative
    difference, the relative differences' 99th percentile and median, relative L2"""
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
    out["nonfinite_differ"] = int(np.count_nonzero(diff & ~fin))
    d = np.abs(rf - af)[fin & diff]
    out["max_abs"] = float(d.max(initial=0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = d / np.abs(rf[fin & diff])
    rel = rel[np.isfinite(rel)]
    out["max_rel"] = float(rel.max(initial=0.0))
    nz = fin & (rf != 0)
    if np.count_nonzero(nz) >= 100:
        allrel = np.abs(rf[nz] - af[nz]) / np.abs(rf[nz])
        for q, key in ((0.5, "rel_p50"), (0.99, "rel_p99")):
            k = int(q * (allrel.size - 1))
            out[key] = float(np.partition(allrel, k)[k])
    nr = float(np.linalg.norm(rf[fin]))
    out["rel_l2"] = float(np.linalg.norm((rf - af)[fin]) / nr) if nr else None
    return out


def _merge(stats: list[dict]) -> dict:
    """the stats of several arrays of one field (frames, tracks) together"""
    out = {"n": sum(s["n"] for s in stats), "differ": sum(s.get("differ", 0) for s in stats)}
    out["arrays"] = len(stats)
    for key in ("max_abs", "max_rel", "rel_p99", "rel_p50", "rel_l2"):
        xs = [s[key] for s in stats if s.get(key) is not None]
        if xs:
            out[key] = max(xs)
    sizes = [s["size"] for s in stats if "size" in s]
    if sizes:
        out["size"] = f"{len(sizes)} arrays of other sizes, e.g. {sizes[0]}"
    nf = sum(s.get("nonfinite_differ", 0) for s in stats)
    if nf:
        out["nonfinite_differ"] = nf
    out["class"] = stats[0].get("class", "")
    return out


def value_stats(r, a, path: str = "", out: dict | None = None) -> dict:
    """the stats of the arrays of two MATLAB values, by path with the struct array elements'
    indices left out (``frameInfo.x``: every frame's x together)"""
    from petakit5d.testing.compare import MArray, MCell, MStruct

    out = {} if out is None else out
    if isinstance(r, MStruct) and isinstance(a, MStruct) and len(r.items) == len(a.items):
        for ri, ai in zip(r.items, a.items, strict=True):
            for k in ri:
                if k in ai:
                    value_stats(ri[k], ai[k], f"{path}.{k}" if path else k, out)
    elif isinstance(r, MCell) and isinstance(a, MCell) and len(r.items or []) == len(a.items or []):
        for ri, ai in zip(r.items or [], a.items or [], strict=True):
            value_stats(ri, ai, path + "{}", out)
    elif isinstance(r, MArray) and isinstance(a, MArray) and r.data is not None:
        if a.data is not None and r.mclass != "char":
            out.setdefault(path, []).append(array_stats(r.data, a.data))
    return out


def counts(v) -> dict:
    """what a detection or track file holds: detections per frame (``frameInfo``), the number
    of tracks (``tracksFinal``, ``tracks``)"""
    from petakit5d.testing.compare import MStruct

    out: dict = {}
    if not isinstance(v, MStruct):
        return out
    top = v.items[0]
    fi = top.get("frameInfo")
    if isinstance(fi, MStruct):
        out["detections"] = [
            int(el["x"].size[1]) if "x" in el and len(el["x"].size) > 1 else 0 for el in fi.items
        ]
    for name in ("tracksFinal", "tracks"):
        t = top.get(name)
        if isinstance(t, MStruct):
            out[name] = len(t.items)
    return out


def file_numbers(mroot: str, mrel: str, proot: str, prel: str) -> dict:
    """the numbers of one output pair: {field path: merged stats}, and what the files hold"""
    from petakit5d.testing.compare import read_mat

    m, p = os.path.join(mroot, mrel), os.path.join(proot, prel)
    low = mrel.lower()
    if low.endswith((".tif", ".tiff")) and os.path.getsize(m) and os.path.getsize(p):
        from petakit5d.io import read_tiff

        s = array_stats(read_tiff(m), read_tiff(p))
        return {"fields": {"voxels": s} if s.get("differ") or "size" in s else {}}
    if low.endswith(".mat"):
        try:
            mv, pv = read_mat(m), read_mat(p)
        except Exception as e:  # an unreadable file: the tree comparison reports it
            return {"fields": {"file": {"error": f"{type(e).__name__}: {e}"}}}
        fields = {
            k: _merge(v)
            for k, v in value_stats(mv, pv).items()
            if any(s.get("differ") or "size" in s for s in v)
        }
        out = {"fields": fields}
        mc, pc = counts(mv), counts(pv)
        if mc or pc:
            out["matlab"], out["python"] = mc, pc
        return out
    return {"fields": {}}


def stats_text(s: dict) -> str:
    if "error" in s:
        return s["error"]
    parts = [f"{s['differ']} of {s['n']} {s.get('class', '')} values differ"]
    if "size" in s:
        parts.append(s["size"])
    if s.get("max_abs") is not None:
        parts.append(f"largest by {s['max_abs']:.3g}")
    if s.get("max_rel") is not None:
        parts.append(f"relative {s['max_rel']:.3g}")
    if s.get("rel_p99") is not None:
        parts.append(f"(99th percentile {s['rel_p99']:.3g}, median {s['rel_p50']:.3g})")
    if s.get("rel_l2") is not None:
        parts.append(f"relative L2 {s['rel_l2']:.3g}")
    if s.get("nonfinite_differ"):
        parts.append(f"{s['nonfinite_differ']} at Inf or NaN")
    return ", ".join(parts)


def counts_text(mc: dict, pc: dict) -> list[str]:
    lines = []
    if "detections" in mc or "detections" in pc:
        md, pd = mc.get("detections", []), pc.get("detections", [])
        lines.append(f"detections: MATLAB {sum(md)}, Python {sum(pd)} in {len(md)} frames")
        bad = [k for k in range(max(len(md), len(pd))) if md[k : k + 1] != pd[k : k + 1]]
        if bad:
            shown = ", ".join(
                f"{k + 1} ({md[k] if k < len(md) else '-'}/{pd[k] if k < len(pd) else '-'})"
                for k in bad[:10]
            )
            lines.append(f"frames with other counts (MATLAB/Python): {len(bad)}: {shown}")
    for name in ("tracksFinal", "tracks"):
        if name in mc or name in pc:
            lines.append(f"{name}: MATLAB {mc.get(name)}, Python {pc.get(name)}")
    return lines


# ---------------------------------------------------------------------------------------------
# the comparison
# ---------------------------------------------------------------------------------------------


def outputs(root: str, inputs: set[str]) -> list[str]:
    """the files under ``root`` that are not inputs (empty folders left out)"""
    from petakit5d.testing.compare import list_tree

    rels = [r for r in list_tree(root) if r not in inputs]
    return [r for r in rels if not os.path.isdir(os.path.join(root, r))]


def name_of(path: str, names) -> str | None:
    """the output a problem's path is about (``files/<name>``, then ``.field``, ``(i)``,
    ``{i}``, `` (voxels)``)"""
    for n in sorted(names, key=len, reverse=True):
        head = "files/" + n
        if path == head or (path.startswith(head) and path[len(head)] in " .({/"):
            return n
    return None


def _root_of(rel: str, root: str, saved: str) -> tuple[str, str]:
    return (saved, rel[len(SAVED) :]) if rel.startswith(SAVED) else (root, rel)


def compare_step(ns, step, mrels, prels, mm, pm) -> dict:
    """one step's outputs: counts, the identical ones, problems beyond the tolerances, and the
    numbers of those that differ"""
    from petakit5d.testing.compare import (
        FileItem,
        FileTree,
        Report,
        compare_trees,
        file_kind,
        normalize_name,
        read_mat,
    )

    def scan(root, saved, rels):
        items: dict = {}
        for rel in sorted(rels):
            base, sub = _root_of(rel, root, saved)
            p = os.path.join(base, sub)
            items.setdefault(normalize_name(rel), []).append(
                FileItem(rel, file_kind(rel, os.path.isdir(p)), Path(p))
            )
        return FileTree(items)

    ref = scan(ns.matlab_root, ns.matlab_saved, mrels)
    act = scan(ns.python_root, ns.python_saved, prels)
    spec = tree_spec(ns.matlab_inputs)
    rep = compare_trees(ref, act, spec, ref_masks=mm, act_masks=pm)
    mnames = {normalize_name(r, mm): r for r in mrels}
    pnames = {normalize_name(r, pm): r for r in prels}
    both = set(mnames) & set(pnames)
    # the lifetime files, field by field
    specs, other = lifetime_specs(ns.matlab_inputs)
    for n in sorted(both):
        if any(fnmatch.fnmatch(n, g) for g in LIFETIME_FILES):
            m, p = _root_of(mnames[n], ns.matlab_root, ns.matlab_saved)
            q, s = _root_of(pnames[n], ns.python_root, ns.python_saved)
            try:
                mv, pv = read_mat(os.path.join(m, p)), read_mat(os.path.join(q, s))
            except Exception as e:  # an unreadable file is a problem
                rep.add(f"files/{n}", "unreadable", f"{type(e).__name__}: {e}")
                continue
            rep.extend(lifetime_problems(mv, pv, n, specs, other, masks=mm, act_masks=pm))
    # the outputs that are not MATLAB's bit for bit: beyond their tolerance, or within it
    paths = [p.path for p in rep.problems]
    paths += [n.split(": ", 1)[0] for n in rep.notes if WITHIN in n]
    differing = sorted({n for p in paths if (n := name_of(p, both)) is not None})
    # the Amira files, by their numbers
    amira = sorted(n for n in both if n.endswith(".am"))
    am_tol = tolerances(ns.matlab_inputs)["*.am"]
    am_stats = {"files": len(amira), "differ": 0, "values": 0, "values_differ": 0, "max": 0.0}
    sub = Report()
    for n in amira:
        mb, mr = _root_of(mnames[n], ns.matlab_root, ns.matlab_saved)
        pb, pr = _root_of(pnames[n], ns.python_root, ns.python_saved)
        with open(os.path.join(mb, mr), encoding="utf-8") as f:
            mt = f.read()
        with open(os.path.join(pb, pr), encoding="utf-8") as f:
            pt = f.read()
        if mt == pt:
            continue
        am_stats["differ"] += 1
        differing.append(n)
        probs, st = amira_problems(mt, pt, am_tol.rel or 0.0, am_tol.abs or 0.0)
        am_stats["values"] += st.get("n", 0)
        am_stats["values_differ"] += st.get("differ", 0)
        am_stats["max"] = max(am_stats["max"], st.get("max_rel_to_1", 0.0))
        for line in probs:
            sub.add(f"files/{n}", "values", line)
    rep.extend(sub)
    numbers = {}
    for n in sorted(set(differing) - set(amira)):
        mb, mr = _root_of(mnames[n], ns.matlab_root, ns.matlab_saved)
        pb, pr = _root_of(pnames[n], ns.python_root, ns.python_saved)
        numbers[n] = file_numbers(mb, mr, pb, pr)
    return {
        "matlab": len(mrels),
        "python": len(prels),
        "identical": len(both) - len(set(differing)),
        "differing": sorted(set(differing)),
        "amira": am_stats if amira else None,
        "numbers": numbers,
        "report": rep,
    }


def _summary(numbers: dict) -> list[str]:
    """the fields that differ over the step's files, merged"""
    by_field: dict[str, list[dict]] = {}
    for res in numbers.values():
        for k, s in res.get("fields", {}).items():
            if "error" not in s:
                by_field.setdefault(k, []).append(s)
    return [f"    {k}: {stats_text(_merge(v))}" for k, v in sorted(by_field.items())]


def main(argv: list[str]) -> int:
    from petakit5d.testing.compare import Masks, list_tree

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("matlab_root")
    ap.add_argument("python_root")
    ap.add_argument("--inputs", action="append", required=True)
    ap.add_argument("--matlab-saved", default="")
    ap.add_argument("--python-saved", default="")
    ap.add_argument("--matlab-inputs", action="store_true")
    ap.add_argument("--matlab-home", default="")
    ap.add_argument("--python-home", default="")
    ap.add_argument("--only", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--label", default="phase3")
    ap.add_argument("--workers", type=int, default=4)
    ns = ap.parse_args(argv)

    t0 = time.time()
    inputs: set[str] = set()
    for root in ns.inputs:
        if os.path.isdir(os.path.join(root, DATASET)):
            inputs |= {f"{DATASET}/{r}" for r in list_tree(os.path.join(root, DATASET))}
    mrels, prels = outputs(ns.matlab_root, inputs), outputs(ns.python_root, inputs)
    for side, saved in ((mrels, ns.matlab_saved), (prels, ns.python_saved)):
        if saved and os.path.isdir(saved):
            side += [SAVED + r for r in outputs(saved, set())]
    if ns.only:
        mrels = [r for r in mrels if re.search(ns.only, r)]
        prels = [r for r in prels if re.search(ns.only, r)]
    mm = Masks().with_paths(work=ns.matlab_home)
    pm = Masks().with_paths(work=ns.python_home)

    by_step: dict[str, tuple[list[str], list[str]]] = {}
    other: list[tuple[str, str]] = []
    for side, rels in ((0, mrels), (1, prels)):
        for r in rels:
            step = step_of(r)
            if step is None:
                other.append(("MATLAB" if side == 0 else "Python", r))
                continue
            by_step.setdefault(step, ([], []))[side].append(r)

    with ThreadPoolExecutor(ns.workers) as ex:
        jobs = {k: ex.submit(compare_step, ns, k, *v, mm, pm) for k, v in by_step.items()}
        results = {k: jobs[k].result() for k in sorted(jobs, key=ORDER.index)}

    on = ", each step on MATLAB's inputs" if ns.matlab_inputs else ", end to end"
    lines = [f"compare {ns.label}{on}: MATLAB {ns.matlab_root}, Python {ns.python_root}"]
    out = {"label": ns.label, "matlab_inputs": ns.matlab_inputs, "steps": {}}
    bad = 0
    for step, res in results.items():
        rep = res.pop("report")
        head = (
            f"\n== {step}: {res['matlab']} MATLAB, {res['python']} Python outputs; "
            f"{res['identical']} identical, {len(res['differing'])} differ"
        )
        lines.append(rep.text(head))
        if res["numbers"]:
            lines.append("  the fields that differ, over the step's files:")
            lines.extend(_summary(res["numbers"]))
        for n, num in res["numbers"].items():
            for cl in counts_text(num.get("matlab", {}), num.get("python", {})):
                lines.append(f"    {n.rsplit('/', 1)[-1]}: {cl}")
        if res["amira"] and res["amira"]["differ"]:
            a = res["amira"]
            lines.append(
                f"  Amira files: {a['differ']} of {a['files']} differ, {a['values_differ']} of "
                f"{a['values']} numbers, the largest by {a['max']:.3g} of max(|v|, 1)"
            )
        out["steps"][step] = {**res, **rep.to_json()}
        bad += len(rep.problems)
    if other:
        lines.append(f"\n== outputs of no step: {len(other)}")
        lines.extend(f"  {who}: {r}" for who, r in other[:20])
        out["other"] = other
        bad += len(other)
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
