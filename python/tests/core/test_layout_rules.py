"""The layout rule: arrays are (z, y, x) C order; MATLAB's Fortran layout only at the boundary.

Outside ``core``, ``io`` (with ``io/matfile.py``) and ``testing``, Python code must not make
Fortran-ordered arrays (``np.asfortranarray``, ``order="F"`` in an array creation or copy) and C++
code must not take ``nb::f_contig`` arrays (``fvolume``, ``fimage``) outside ``cpp/common`` and
``cpp/io`` (which reads MATLAB-layout files).
``reshape``/``ravel``/``flatten``/``tobytes``/``ravel_multi_index``/``unravel_index`` with
``order="F"`` stay allowed: they put values in MATLAB's linear order, which parameters and
order-dependent steps need.

``KNOWN`` lists the uses that predate the rule, with the task that removes them; an entry that
no longer matches fails, so the list only shrinks.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]  # python/
SRC = ROOT / "src" / "petakit5d"
ALLOWED_PY = ("core", "io", "testing")
ELEMENT_ORDER = re.compile(r"\.(reshape|ravel|flatten|tobytes)\(|ravel_multi_index|unravel_index")
FORTRAN_PY = re.compile(r"asfortranarray|order\s*=\s*['\"]F['\"]")
FORTRAN_CPP = re.compile(r"\bf_contig\b|\bfvolume\b|\bfimage\b")

# path relative to python/ -> number of lines, and who removes them
KNOWN = {
    "src/petakit5d/matlabcompat/fft.py": (
        1,
        "F18: the exact FFT backend works in MATLAB's layout behind the interface",
    ),
    "src/petakit5d/gpu/fft.py": (3, "F18/R04: GPU FFT plans in MATLAB's layout"),
    "src/petakit5d/gpu/selftest.py": (3, "F18/R04: the GPU FFT self-test"),
}
KNOWN_CPP = {"cpp/fft/module.cpp": (1, "F18: the exact FFT backend")}
ALLOWED_CPP = ("common", "io")


def _py_hits(path: Path) -> int:
    n = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        code = line.split("#", 1)[0]
        if FORTRAN_PY.search(code) and not ELEMENT_ORDER.search(code):
            n += 1
    return n


def test_no_fortran_layout_outside_the_boundary():
    found = {}
    for p in SRC.rglob("*.py"):
        rel = p.relative_to(SRC).parts
        if rel[0] in ALLOWED_PY:
            continue
        n = _py_hits(p)
        if n:
            found[str(p.relative_to(ROOT))] = n
    unexpected = {k: v for k, v in found.items() if k not in KNOWN or v > KNOWN[k][0]}
    assert not unexpected, (
        "Fortran-ordered arrays outside core/io/testing (use core.arrays.from_matlab / "
        f"to_matlab): {unexpected}"
    )
    stale = {k: v[0] for k, v in KNOWN.items() if found.get(k, 0) < v[0]}
    assert not stale, f"remove or lower these entries of KNOWN (fewer uses now): {stale}"


def test_no_fortran_arrays_in_kernels():
    found = {}
    for p in (ROOT / "cpp").rglob("*"):
        if (
            p.suffix not in (".cpp", ".hpp", ".h", ".cu", ".cuh")
            or p.relative_to(ROOT / "cpp").parts[0] in ALLOWED_CPP
        ):
            continue
        n = sum(
            1
            for line in p.read_text(encoding="utf-8", errors="replace").splitlines()
            if FORTRAN_CPP.search(line.split("//", 1)[0])
        )
        if n:
            found[str(p.relative_to(ROOT))] = n
    unexpected = {k: v for k, v in found.items() if k not in KNOWN_CPP or v > KNOWN_CPP[k][0]}
    assert not unexpected, (
        f"kernels take cvolume/cimage (C order), not Fortran arrays: {unexpected}"
    )
    stale = {k: v[0] for k, v in KNOWN_CPP.items() if found.get(k, 0) < v[0]}
    assert not stale, f"remove or lower these entries of KNOWN_CPP: {stale}"


def test_the_rule_sees_what_it_should():
    assert FORTRAN_PY.search("x = np.asfortranarray(a)")
    assert FORTRAN_PY.search('b = np.empty(s, order="F")')
    assert not (
        FORTRAN_PY.search('v = a.ravel(order="F")')
        and not ELEMENT_ORDER.search('v = a.ravel(order="F")')
    )
    assert FORTRAN_CPP.search("void f(pk5d::fvolume<float> v)")
