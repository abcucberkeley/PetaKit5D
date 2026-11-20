"""The kernels' outputs on fixed inputs, pinned by hash: the same bits on Linux, macOS and Windows
(docs/platforms.md). The cases are in _kernel_cases.py, the pins in kernel_hashes.json (made on
Linux x86-64 with ``python tests/platform/_kernel_cases.py --write``). A case that fails here
names the kernels whose bits differ on this system; ``--save`` on both systems and ``--compare``
measure the difference.

The cases in RECORDED (FFTs, dense LAPACK) are recorded, not pinned: their last bits follow the
system's C library and BLAS. The cases in ARM_NAN are expected to fail where the CPU's default
NaN is not x86's (ARM), until their kernels give the same NaN on both (_kernel_cases.py has
why)."""

import importlib
import json
import sys
from pathlib import Path

import pytest

from petakit5d.core import dtypes

sys.path.insert(0, str(Path(__file__).parent))
import _kernel_cases as kc  # noqa: E402

PINS = json.loads(kc.HASHES.read_text())
PINNED = sorted(n for n in kc.CASES if n not in kc.RECORDED)
MATLABCOMPAT = (
    "distance",
    "elementary",
    "fft",
    "filters",
    "graph",
    "label",
    "lsq",
    "morphology",
    "neighbors",
    "reduce",
    "resize",
    "stats",
    "warp",
)


@pytest.mark.parametrize("name", PINNED)
def test_outputs_have_the_pinned_bits(name, request):
    if name in kc.ARM_NAN and not dtypes._cpu_nan_is_x86():
        request.applymarker(pytest.mark.xfail(strict=True, reason=kc.ARM_NAN[name]))
    got = kc.digest(kc.run(name))
    want = PINS["pinned"].get(name)
    assert want is not None, f"{name}: no pin (python tests/platform/_kernel_cases.py --write)"
    assert got == want, (
        f"{name}: the outputs differ from Linux x86-64's (digest {got}, pinned {want}); "
        "measure with _kernel_cases.py --save here and on Linux, then --compare"
    )


@pytest.mark.parametrize("name", sorted(kc.RECORDED))
def test_recorded_outputs_repeat(name):
    # not pinned (their last bits follow the C library or the BLAS), but the same on a second run
    assert kc.digest(kc.run(name)) == kc.digest(kc.run(name))


def test_the_pins_are_those_of_the_cases():
    assert sorted(PINS["pinned"]) == PINNED
    assert sorted(PINS["recorded"]) == sorted(kc.RECORDED)
    assert set(kc.ARM_NAN) <= set(PINNED)


def _record(monkeypatch, module, names, prefix, called):
    for n in names:
        f = getattr(module, n)

        def wrap(*args, _f=f, _key=f"{prefix}{n}", **kw):
            called.add(_key)
            return _f(*args, **kw)

        monkeypatch.setattr(module, n, wrap)


def test_every_binding_is_covered(monkeypatch):
    from petakit5d._ext import available, groups, load

    called: set[str] = set()
    for g in groups():
        m = load(g)
        names = [
            n
            for n in dir(m)
            if not n.startswith("__") and not n.isupper() and callable(getattr(m, n))
        ]
        _record(monkeypatch, m, names, f"{g}.", called)
    for name in kc.CASES:
        kc.run(name)
    monkeypatch.undo()
    names = set(kc.binding_names())
    exempt = set(kc.EXEMPT_BINDINGS) | set(kc.exact_only_bindings())
    assert sorted(names - called - exempt) == [], "bindings no case calls"
    # exemptions of modules built here that name no binding (deskew_scalar: x86-64 only)
    stale = {e for e in kc.EXEMPT_BINDINGS if available(e.split(".")[0]) and e not in names}
    assert sorted(stale) == []


def test_every_matlabcompat_function_is_covered(monkeypatch):
    called: set[str] = set()
    every = set()
    for mod in MATLABCOMPAT:
        m = importlib.import_module(f"petakit5d.matlabcompat.{mod}")
        names = [
            n for n in m.__all__ if callable(getattr(m, n)) and not isinstance(getattr(m, n), type)
        ]
        every |= {f"{mod}.{n}" for n in names}
        _record(monkeypatch, m, names, f"{mod}.", called)
    for name in kc.CASES:
        kc.run(name)
    monkeypatch.undo()
    assert sorted(every - called - set(kc.EXEMPT_FUNCTIONS)) == [], "functions no case calls"
    assert sorted(set(kc.EXEMPT_FUNCTIONS) - every) == []
