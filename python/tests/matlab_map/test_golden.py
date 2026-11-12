"""The golden references of the MATLAB suite (tests/golden/data/*.mat) as oracle cases.

GoldenOutputTest/outputMatchesReference compares 8 scenarios with references stored in the
MATLAB repo. Each scenario has three MATLAB cases (pkpy_golden.m in matlab/):

  golden.<scenario>.stored   the stored reference; Python loads the same .mat (exact)
  golden.<scenario>.compute  the scenario computed by the oracle tree (GoldenOutputTest)
  golden.<scenario>.inputs   the arrays the scenario starts from (MATLAB's rng(1, 'twister'))

An area task tests its port on a scenario with the inputs case's outputs and compares with
the compute case (exact) and the stored reference (relative L2 error <= 1e-4, as the MATLAB
test):

    ref = oracle.reference(GOLDEN_CASES['golden.decon_omw.inputs'])   # outputs: im, psf

Until the oracle harness (F04) is in the tree, the cases are stand-ins and their tests skip.
"""

from __future__ import annotations

import dataclasses
import os
from pathlib import Path
from typing import Any

import numpy as np
import pytest

try:  # the oracle harness (F04)
    from petakit5d.testing import cases as _cases
    from petakit5d.testing import oracle as _oracle
except ImportError:
    _cases = _oracle = None

# scenario: (output fields of the stored and computed struct, MATLAB classes, sizes; input fields)
SCENARIOS: dict[str, dict[str, Any]] = {
    "deskew_rotate_sample_scan": {"outputs": {"dsr": ("uint16", (32, 112, 21))}, "inputs": ["vol"]},
    "deskew_rotate_skewed_interp": {
        "outputs": {"dsr": ("uint16", (32, 148, 21))},
        "inputs": ["vol"],
    },
    "deskew_rotate_resampled": {"outputs": {"dsr": ("single", (16, 56, 11))}, "inputs": ["vol"]},
    "decon_rl_simplified": {
        "outputs": {"decon": ("single", (32, 32, 32))},
        "inputs": ["im", "psf"],
    },
    "decon_omw": {
        "outputs": {"backProjector": ("single", (21, 21, 21)), "decon": ("single", (32, 32, 32))},
        "inputs": ["im", "psf"],
    },
    "psf_gen_masked": {"outputs": {"psf": ("single", (63, 63, 41))}, "inputs": ["psf"]},
    "stitch_feather_xcorr": {"outputs": {"stitched": ("uint16", (38, 114, 38))}, "inputs": ["gt"]},
    "point_detection": {
        "outputs": {
            "amplitude": ("double", (5, 1)),
            "background": ("double", (5, 1)),
            "position": ("double", (5, 3)),
        },
        "inputs": ["vol"],
    },
}
WHAT = ("stored", "compute", "inputs")
HELPERS = ["tests/matlab_map/matlab"]
REL_TOL = 1e-4  # GoldenOutputTest's tolerance between platforms
_DTYPES = {"uint16": np.uint16, "single": np.float32, "double": np.float64}


def stored_path(scenario: str) -> Path | None:
    """the .mat of a scenario in the oracle tree (PETAKIT5D_ORACLE_DIR), or None"""
    root = os.environ.get("PETAKIT5D_ORACLE_DIR")
    if not root:
        return None
    p = Path(root, "tests", "golden", "data", f"{scenario}.mat")
    return p if p.is_file() else None


def load_stored(scenario: str) -> list[np.ndarray]:
    """the stored reference of a scenario, its fields in SCENARIOS' order (Python side of the
    stored case)"""
    import scipy.io

    p = stored_path(scenario)
    if p is None:
        raise FileNotFoundError(f"no reference for {scenario}: PETAKIT5D_ORACLE_DIR")
    m = scipy.io.loadmat(p, squeeze_me=False, struct_as_record=False)
    e = m["expected"][0, 0]
    return [np.asfortranarray(getattr(e, f)) for f in SCENARIOS[scenario]["outputs"]]


# ------------------------------------------------------------------------------------------
# cases (F04's Case, or a stand-in with the same fields until F04 lands)
# ------------------------------------------------------------------------------------------


@dataclasses.dataclass(eq=False)
class _StandinCase:
    id: str
    matlab: str | None = None
    args: tuple = ()
    nargout: int = 1
    python: Any = None
    py_args: tuple | None = None
    helpers: tuple = ()
    files: bool | None = None
    outputs: tuple | None = None
    note: str = ""


class _StandinWork:
    def __repr__(self) -> str:
        return "Work()"


def _make_case(**kw: Any) -> Any:
    if _cases is not None:
        return _cases.Case(**kw)
    return _StandinCase(**kw)


def _mcell(*items: str) -> Any:
    return _cases.mcell(*items) if _cases is not None else tuple(items)


def _work() -> Any:
    return _cases.Work() if _cases is not None else _StandinWork()


def golden_cases() -> list[Any]:
    out = []
    for s, spec in SCENARIOS.items():
        for what in WHAT:
            fields = list(spec["inputs"] if what == "inputs" else spec["outputs"])
            out.append(
                _make_case(
                    id=f"golden.{s}.{what}",
                    matlab="pkpy_golden",
                    args=[what, s, _mcell(*fields), _work()],
                    nargout=len(fields),
                    outputs=fields,
                    python=load_stored if what == "stored" else None,
                    py_args=[s] if what == "stored" else None,
                    helpers=HELPERS,
                    files=False,
                    note=f"GoldenOutputTest/outputMatchesReference(scenario={s}): {what}",
                )
            )
    return out


CASES = golden_cases()
GOLDEN_CASES = {c.id: c for c in CASES}


def _need_harness() -> None:
    if _oracle is None:
        pytest.skip("the oracle harness (petakit5d.testing.oracle, F04) is not in this tree")


def _need_tree() -> None:
    if not os.environ.get("PETAKIT5D_ORACLE_DIR"):
        pytest.skip(
            "PETAKIT5D_ORACLE_DIR is not set (the stored references are in the oracle tree)"
        )


# ------------------------------------------------------------------------------------------
# tests
# ------------------------------------------------------------------------------------------


def test_every_scenario_has_its_cases():
    ids = [c.id for c in CASES]
    assert len(ids) == len(set(ids)) == 3 * len(SCENARIOS)
    for c in CASES:
        assert c.nargout == len(c.outputs)


def _oracle_cases():
    if _oracle is None:
        return pytest.mark.parametrize("case", [pytest.param(c, id=c.id) for c in CASES])
    return _oracle.parametrize(CASES)


@_oracle_cases()
def test_golden_case(case, tmp_path):
    """MATLAB's result of each case exists (store or manifest); for the stored references,
    Python's reading of the .mat equals MATLAB's load"""
    _need_harness()
    if case.python is not None:
        _need_tree()
    _oracle.check(case, tmp_path)


@pytest.mark.parametrize("scenario", list(SCENARIOS))
def test_stored_reference_has_the_expected_fields(scenario):
    _need_tree()
    if stored_path(scenario) is None:
        pytest.fail(f"no tests/golden/data/{scenario}.mat in the oracle tree")
    arrays = load_stored(scenario)
    for (field, (mclass, size)), a in zip(
        SCENARIOS[scenario]["outputs"].items(), arrays, strict=True
    ):
        assert a.dtype == _DTYPES[mclass], field
        assert a.shape == size, field


def rel_l2(actual: np.ndarray, expected: np.ndarray) -> float:
    """GoldenOutputTest's measure: norm(a - e) / max(norm(e), eps), in double"""
    import math

    a = np.asarray(actual, dtype=np.float64).reshape(-1, order="F")
    e = np.asarray(expected, dtype=np.float64).reshape(-1, order="F")
    d = math.sqrt(math.fsum((a - e) ** 2))
    n = math.sqrt(math.fsum(e**2))
    return d / max(n, np.finfo(np.float64).eps)


@pytest.mark.oracle
@pytest.mark.parametrize("scenario", list(SCENARIOS))
def test_oracle_tree_reproduces_the_stored_reference(scenario):
    """the oracle tree's scenario agrees with the stored reference as GoldenOutputTest
    requires (class, size, relative L2 error <= 1e-4)"""
    _need_harness()
    computed = _oracle.reference(GOLDEN_CASES[f"golden.{scenario}.compute"])
    stored = _oracle.reference(GOLDEN_CASES[f"golden.{scenario}.stored"])
    if computed.source != "store" or stored.source != "store":
        pytest.skip("needs the golden store (the manifests hold hashes only)")
    for name, c, s in zip(
        SCENARIOS[scenario]["outputs"], computed.outputs, stored.outputs, strict=True
    ):
        assert c.mclass == s.mclass, name
        assert tuple(c.size) == tuple(s.size), name
        assert rel_l2(c.data, s.data) <= REL_TOL, name
