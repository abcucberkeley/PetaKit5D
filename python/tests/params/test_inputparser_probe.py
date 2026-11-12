"""The MATLAB-name layer binds arguments as inputParser does (stored MATLAB R2024b results).

oracle_standin/pkpy_f03_inputparser_probe.m ran these calls through an inputParser with
CaseSensitive false; the same calls go through a registered spec here, and the Results
(names, order, values) or the error identifier and message must agree. Our messages add a
"Did you mean" hint, which is stripped before comparing.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from petakit5d.params import MatlabArgumentError, MString, Param, Registry
from petakit5d.params.values import describe

pytestmark = pytest.mark.oracle

PROBE = json.loads((Path(__file__).parent / "data" / "inputparser_probe.json").read_text())

_NAMED = [
    Param("zeta", "@isnumeric"),
    Param("Alpha", "@isnumeric"),
    Param("beta", "@(x) islogical(x) || ischar(x)"),
    Param("save16bit", "@islogical"),
    Param("saveZarr", "@islogical"),
    Param("deskew", "@islogical"),
    Param("deskewRotate", "@islogical"),
    Param("Gamma", "@iscell"),
]


def _specs():
    reg = Registry()

    @reg.register("probe_num", Param("req", "@isnumeric"), Param("opt", "@isnumeric"), *_NAMED)
    def probe_num(req, opt=0, *, zeta=0, alpha=0, beta=False, save_16bit=True, save_zarr=False,
                  deskew=False, deskew_rotate=False, gamma=()):  # fmt: skip
        return locals()

    @reg.register("probe_char", Param("req", "@isnumeric"), Param("opt", "@ischar"), *_NAMED)
    def probe_char(req, opt="dflt", *, zeta=0, alpha=0, beta=False, save_16bit=True,
                   save_zarr=False, deskew=False, deskew_rotate=False, gamma=()):  # fmt: skip
        return locals()

    @reg.register("probe_none", Param("req", "@isnumeric"), *_NAMED)
    def probe_none(req, *, zeta=0, alpha=0, beta=False, save_16bit=True, save_zarr=False,
                   deskew=False, deskew_rotate=False, gamma=()):  # fmt: skip
        return locals()

    return {k: reg.find(f"probe_{k}") for k in ("num", "char", "none")}


SPECS = _specs()

# the calls of the probe, in Python (struct -> dict, "x" -> MString)
CALLS = {
    "base": [1],
    "opt_given": [1, 7],
    "opt_then_pair": [1, 7, "zeta", 5],
    "pair": [1, "zeta", 5],
    "case": [1, "ZETA", 5],
    "prefix": [1, "zet", 5],
    "prefix_mixed_case": [1, "al", 7],
    "ambiguous": [1, "save", True],
    "exact_over_prefix": [1, "deskew", True],
    "prefix_of_two_with_exact": [1, "deskewR", True],
    "unknown": [1, "foo", 1],
    "unknown_close": [1, "zetta", 1],
    "odd": [1, "zeta"],
    "dup": [1, "zeta", 1, "zeta", 2],
    "dup_case": [1, "zeta", 1, "ZETA", 2],
    "struct": [1, {"zeta": 3, "Alpha": 4}],
    "struct_then_pair": [1, {"zeta": 3}, "zeta", 9],
    "pair_then_struct": [1, "zeta", 9, {"zeta": 3}],
    "struct_unknown": [1, {"foo": 3}],
    "struct_case": [1, {"ZETA": 3}],
    "struct_prefix": [1, {"zet": 3}],
    "struct_mixed": [1, {"zeta": 3, "foo": 4}],
    "invalid": [1, "zeta", "abc"],
    "invalid_anon": [1, "beta", 1],
    "valid_anon_char": [1, "beta", "x"],
    "name_is_required": [1, "req", 5],
    "name_is_optional": [1, "opt", 5],
    "prefix_of_required": [1, "re", 5],
    "opt_invalid_char": [1, "hello", 5],
    "opt_invalid_struct": [1, {"zeta": 3}],
    "nonchar_name": None,  # per parser below
    "string_name": [1, MString("zeta"), 5],
    "empty_name": None,
    "required_missing": [],
    "required_invalid": ["x"],
    "cell_param": [1, "Gamma", ["a", "b"]],
    "cell_param_invalid": [1, "Gamma", "a"],
    "logical_param_numeric": [1, "save16bit", 1],
    "numeric_param_logical": [1, "zeta", True],
    "numeric_param_empty": [1, "zeta", None],
    "numeric_param_empty_char": [1, "zeta", ""],
    "opt_is_param_name": [1, "zeta", 5],
    "opt_string": [1, MString("zeta"), 5],
    "opt_nonchar": [1, 5],
    "pair_unknown_after": [1, "zeta", 1, "foo", 2],
    "string_value": [1, "Gamma", [MString("a")]],
    "too_many_positional": [1, 2],
}
PER_PARSER = {
    ("num", "nonchar_name"): [1, 7, 5, 5],
    ("num", "empty_name"): [1, 7, "", 5],
    ("none", "nonchar_name"): [1, 5, 5],
    ("none", "empty_name"): [1, "", 5],
    ("char", "opt_given"): [1, "hello"],
    ("char", "opt_then_pair"): [1, "hello", "zeta", 5],
}
# where the port differs on purpose
DEVIATIONS = {
    # ip.parse with nothing: inputParser's own id; a function with a required argument
    # fails earlier in MATLAB, with MATLAB:minrhs, as the port does
    ("num", "required_missing"): "minrhs",
    # values are normalized for the Python functions: a string scalar becomes a str (char)
    ("none", "string_value"): "string->char",
}


def _strip_hint(msg: str) -> str:
    return re.sub(r" Did you mean '[^']*'\?$", "", msg)


@pytest.mark.parametrize("r", PROBE, ids=[f"{r['parser']}:{r['name']}" for r in PROBE])
def test_binding_matches_inputparser(r):
    spec = SPECS[r["parser"]]
    args = PER_PARSER.get((r["parser"], r["name"]), CALLS[r["name"]])
    dev = DEVIATIONS.get((r["parser"], r["name"]))
    try:
        bound = spec.bind_matlab(args)
    except MatlabArgumentError as err:
        assert "error" in r, f"MATLAB accepted it: {r.get('values')}"
        if dev == "minrhs":
            assert err.identifier == "MATLAB:minrhs"
            return
        assert err.identifier == r["error"]
        assert _strip_hint(str(err)) == r["message"]
        return
    assert "error" not in r, f"MATLAB: {r.get('error')} {r.get('message')}"
    res = spec.matlab_results(bound)
    assert list(res) == r["fields"]
    got = {k: describe(v) for k, v in res.items()}
    if dev == "string->char":
        assert got["Gamma"]["d"][0]["c"] == "char"
        got["Gamma"]["d"][0]["c"] = "string"
    assert got == r["values"]
    defaults = sorted(p.matlab for p in spec.params if p.py not in bound)
    assert defaults == sorted(r["using_defaults"])


def test_unknown_names_get_a_suggestion():
    spec = SPECS["none"]
    with pytest.raises(MatlabArgumentError, match=r"Did you mean 'zeta'\?"):
        spec.bind_matlab([1, "zetta", 1])
    with pytest.raises(MatlabArgumentError, match=r"Did you mean 'deskewRotate'\?"):
        spec.bind_matlab([1, "deskewRotat3", True])
