"""parameters.json keeps MATLAB's values: a MATLAB-style call reaches the snake_case function
with Python's values, and ``parameters_json`` of what the function got is the text MATLAB's
``jsonencode(ip.Results, PrettyPrint=true)`` gives for the same call.

MATLAB side: XR_crop_dataset's inputParser block (from the stored extract of the oracle tree)
parses the case's arguments; Python side: the worked example's spec (sample_crop.py).
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

from petakit5d.testing import oracle
from petakit5d.testing.cases import Case

HERE = Path(__file__).parent


def _sample():
    if "sample_crop" in sys.modules:
        return sys.modules["sample_crop"]
    spec = importlib.util.spec_from_file_location("sample_crop", HERE / "sample_crop.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["sample_crop"] = mod
    spec.loader.exec_module(mod)
    return mod


def _parser_code(name: str) -> str:
    fn = json.loads((HERE / "data" / "inputparser_extract.json").read_text())["functions"][name]
    lines = ["ip = inputParser;", "ip.CaseSensitive = false;"]
    for p in fn["params"]:
        if p["kind"] == "required":
            lines.append(f"ip.addRequired('{p['name']}', {p['validator']});")
        else:
            lines.append(f"ip.addParameter('{p['name']}', {p['default']}, {p['validator']});")
    lines += ["ip.parse(in{:});", "out = {jsonencode(ip.Results, PrettyPrint=true)};"]
    return "\n".join(lines)


CODE = _parser_code("XR_crop_dataset")
ROOT = "/mnt/projects/embryo_2024"


def py_parameters_json(*args):
    spec = _sample().REG.find("XR_crop_dataset")
    got = spec.call_matlab(*args)  # the function returns the (Python) values it got
    return spec.parameters_json(got)


CALLS = {
    "defaults": [ROOT, [1, 2, 3, 100, 200, 50]],
    "sizes_and_cells": [
        [ROOT, "/mnt/projects/b"], [5, 6, 7, 8, 9, 10], "batchSize", [512, 256, 64],
        "BlockSize", [32, 64, 16], "channelPatterns", ["CamA_ch0"], "pad", True,
    ],
    "open_ends_and_point": [
        ROOT, [1, 1, 1, math.inf, 300, math.inf], "cropType", "moving",
        "lastStartCoords", [10, 20, 30], "cpusPerTask", 8, "uuid", "abc",
    ],
    "prefixes_and_case": [ROOT, [2, 2, 2, 9, 9, 9], "BATCHS", [100, 200, 300], "resultdir", "C"],
}  # fmt: skip

CASES = [
    Case(
        f"params_json_{name}",
        code=CODE,
        args=args,
        python=py_parameters_json,
        files=False,
    )
    for name, args in CALLS.items()
]


@oracle.parametrize(CASES)
def test_parameters_json_matches_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_python_function_gets_python_values():
    spec = _sample().REG.find("XR_crop_dataset")
    got = spec.call_matlab(*CALLS["open_ends_and_point"])
    assert got["input_bbox"] == (0, 0, 0, None, None, 300)
    assert got["last_start_coords"] == (29, 9, 19)
    got = spec.call_matlab(*CALLS["sizes_and_cells"][:6])
    assert got["batch_size"] == (64, 512, 256)


@pytest.mark.parametrize("name", list(CALLS))
def test_parameters_json_round_trips(name):
    """The MATLAB values come back from the Python values unchanged (no MATLAB needed)."""
    spec = _sample().REG.find("XR_crop_dataset")
    from petakit5d.matlab import MatlabFunction
    from petakit5d.params.values import describe, to_python

    got = MatlabFunction(spec)(*CALLS[name])
    res = spec.matlab_results(got)
    kw = {k: to_python(v) for k, v in res.items()}
    args = [kw.pop("dataPaths"), kw.pop("inputBbox")]
    again = spec.matlab_results(spec.bind_matlab(args, kw))
    assert {k: describe(v) for k, v in res.items()} == {k: describe(v) for k, v in again.items()}
