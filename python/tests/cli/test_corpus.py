"""The argument corpus: every call the GUI and PyPetaKit5D make for the Phase 1 functions (and
the edge cases of the grammar) binds in pk5d-master to the values MATLAB's compiled master
gives the function.

MATLAB side (tests/cli/matlab/pkpy_j12_master.m): the oracle tree's mccMaster with the
strings, its mcc parser, and the task function replaced by a recorder; the recorded arguments
are then checked by the function's inputParser (from the stored extract). The result is a
struct of what the parser passes (pos1..posN, then each name-value pair by name, vectors as
rows) and mccMaster's "Function name / Variables" lines.

Python side: ``master.bind`` (the registry's ``bind_strings``) and the values in MATLAB's form
(``spec.matlab_results``: the port's values converted back by their kinds). Parameters whose
MATLAB default is site-specific (Slurm options, the compiled master's paths) are left out on
both sides unless the call gives them. Functions that are not ported yet are skipped (their
references are made already).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from petakit5d.cli import master
from petakit5d.params import values as mv
from petakit5d.params.registry import REGISTRY
from petakit5d.params.spec import MatlabArgumentError
from petakit5d.testing import compare as cmp
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case

HERE = Path(__file__).parent
EXTRACT = json.loads((HERE.parent / "params" / "data" / "inputparser_extract.json").read_text())[
    "functions"
]


def _callers():
    if "pk5d_cli_callers" in sys.modules:
        return sys.modules["pk5d_cli_callers"]
    spec = importlib.util.spec_from_file_location("pk5d_cli_callers", HERE / "callers.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pk5d_cli_callers"] = mod
    spec.loader.exec_module(mod)
    return mod


C = _callers()
ROOT = C.ROOT

# ----------------------------------------------------------------------------------------------
# the MATLAB side's inputs


def _is_site(text: str | None) -> bool:
    return isinstance(text, str) and text.startswith("<site value")


def _tf(b: bool) -> str:
    return "true" if b else "false"


def ip_code(function: str) -> str:
    """MATLAB code that makes the function's inputParser as ``ip`` (from the stored extract;
    a site-specific default becomes '')."""
    e = EXTRACT[function]
    lines = [
        "ip = inputParser;",
        f"ip.CaseSensitive = {_tf(e['case_sensitive'])};",
        f"ip.KeepUnmatched = {_tf(e['keep_unmatched'])};",
        f"ip.PartialMatching = {_tf(e['partial_matching'])};",
        f"ip.StructExpand = {_tf(e['struct_expand'])};",
    ]
    for p in e["params"]:
        v = f", {p['validator']}" if p["validator"] else ""
        d = "''" if _is_site(p["default"]) else p["default"]
        if p["kind"] == "required":
            lines.append(f"ip.addRequired('{p['name']}'{v});")
        elif p["kind"] == "optional":
            lines.append(f"ip.addOptional('{p['name']}', {d}{v});")
        else:
            lines.append(f"ip.addParameter('{p['name']}', {d}{v});")
    return "\n".join(lines)


def site_names(function: str, argv: list[str]) -> list[str]:
    """Parameters whose MATLAB default is site-specific and that the call does not give."""
    e = EXTRACT[function]
    given = {a.lower() for a in argv}
    defaults = (e.get("parser") or {}).get("defaults", {})
    return [
        p["name"]
        for p in e["params"]
        if (_is_site(p["default"]) or _is_site(defaults.get(p["name"])))
        and p["name"].lower() not in given
    ]


# ----------------------------------------------------------------------------------------------
# the Python side


_NP = {"double": "float64", "single": "float32", "logical": "bool", "char": "<u2"}


def mvalue(v):
    """A model value of petakit5d.params.values as the harness's MATLAB value (rows for
    vectors, as the MATLAB side gives them)."""
    if isinstance(v, mv.MString):
        return cmp.MString((1, 1), [str(v)])
    if isinstance(v, mv.Cell):
        shape = v.shape
        if len(shape) == 2 and shape[1] == 1 and shape[0] > 1:
            shape = (1, shape[0])
        return cmp.MCell(shape, [mvalue(e) for e in v])
    if isinstance(v, dict):
        return cmp.MStruct((1, 1), list(v), [{k: mvalue(e) for k, e in v.items()}])
    shape = v.shape
    if len(shape) == 2 and shape[1] == 1 and shape[0] > 1:
        shape = (1, shape[0])
    data = [ord(c) for c in v.data] if v.mclass == "char" else v.data
    dtype = _NP.get(v.mclass, v.mclass)
    return cmp.MArray(v.mclass, np.array(data, dtype=dtype).reshape(shape, order="F"))


def py_master(ipcode, drop, function, *argv):
    """pk5d-master's binding of the strings, as the MATLAB side reports it."""
    spec = master.find_function(function)
    kwargs = master.bind(spec, argv)
    res = spec.matlab_results(kwargs)
    low = {d.lower() for d in drop}
    fields: dict = {}
    k = 0
    for p in spec.params:
        if p.python_only:
            continue
        if p.kind != "parameter":
            k += 1
            fields[f"pos{k}"] = mvalue(res[p.matlab])
        elif p.matlab.lower() not in low:
            fields[p.matlab] = mvalue(res[p.matlab])
    return cmp.MStruct((1, 1), list(fields), [fields]), master.header(function, list(argv))


# ----------------------------------------------------------------------------------------------
# the corpus

P1 = f"{ROOT}/cell1"
BBOX = "[1,1,1,100,200,50]"

# (name, function, arguments, MATLAB's error identifier or None, or True for an error without one)
EDGE = [
    ("names_in_any_case", "XR_crop_dataset",
     ["{'" + P1 + "'}", BBOX, "BATCHSIZE", "[64,128,32]", "Pad", "true", "channelpatterns",
      "{'CamA_ch0'}"], None),
    ("unique_prefixes", "XR_crop_dataset",
     ["{'" + P1 + "'}", BBOX, "batchS", "[64,128,32]", "resultDir", "out", "cpus", "3"], None),
    ("last_value_wins", "XR_crop_dataset",
     ["{'" + P1 + "'}", BBOX, "batchSize", "[1,2,3]", "batchSize", "[4,5,6]"], None),
    ("empty_text_value", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "jobLogDir", ""], None),
    ("blank_separated_numbers", "XR_crop_dataset",
     ["{'" + P1 + "'}", "[1 1 1 100 200 50]", "blockSize", "[32 64  16]"], None),
    ("column_and_range", "XR_crop_dataset",
     ["{'" + P1 + "'}", "[1;1;1;100;200;50]", "lastStartCoords", "[1:3]"], None),
    ("number_forms", "XR_crop_dataset",
     ["{'" + P1 + "'}", BBOX, "cpusPerTask", "1e1", "batchSize", "[2^9, 256, 1.5e2]"], None),
    ("hex_literal", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "batchSize", "[512, 0x100]"],
     None),
    ("text_path_not_cell", "XR_crop_dataset", [P1, BBOX], None),
    ("required_by_name", "XR_crop_dataset",
     ["{'" + P1 + "'}", BBOX, "inputBbox", "[1,1,1,5,5,5]", "dataPaths", "{'/mnt/x'}"], None),
    ("blank_separated_cell", "XR_crop_dataset",
     ["{'" + P1 + "' '" + ROOT + "/y'}", BBOX, "channelPatterns", "{'CamA' 'CamB'}"], None),
    ("empty_values", "XR_crop_dataset",
     ["{'" + P1 + "'}", BBOX, "lastStartCoords", "[]", "channelPatterns", "{}",
      "resultDirName", "''"], None),
    ("quote_in_cell", "XR_crop_dataset",
     ["{'" + ROOT + "/it''s', '" + ROOT + "/b'}", BBOX], None),
    ("column_cell", "XR_crop_dataset", ["{'" + P1 + "';'" + ROOT + "/x'}", BBOX], None),
    ("logical_array_for_numeric", "XR_MIP_wrapper", ["{'" + P1 + "'}", "axis", "[true,false,true]"],
     "MATLAB:InputParser:ArgumentFailedValidation"),  # axis is @isnumeric: the GUI sends 0/1
    ("inf_box_ends", "XR_crop_dataset", ["{'" + P1 + "'}", "[1,1,1,Inf,inf,50]"], None),
    ("nan_box_end", "XR_crop_dataset", ["{'" + P1 + "'}", "[1,1,1,Inf,inf,NaN]"], None),
    ("optional_positional", "XR_RLdeconFrame3D",
     ["{'" + P1 + "/a.tif'}", "0.108", "0.5", "", "PSFfile", f"{ROOT}/psf.tif",
      "Background", "100", "RLMethod", "omw", "wienerAlpha", "0.005", "GPUJob", "false"], None),
    ("optional_positional_given", "tiffToZarr",
     ["{'" + P1 + "/a.tif'}", f"{ROOT}/a.zarr", "[]", "blockSize", "[64,64,32]"], None),
    ("frame_task", "XR_deskewRotateFrame",
     [P1 + "/a.tif", "0.108", "0.5", "SkewAngle", "32.45", "Reverse", "true", "Rotate", "true",
      "resampleFactor", "[1,1,1]", "inputBbox", "[]"], None),
    ("no_name_value_pairs_extra_ignored", "stitch_global_grid_assignment_wrapper",
     [f"{ROOT}/in.mat", f"{ROOT}/out.mat", "0.25", "[1,0.1,10]", "u1", "surplus", "x"], None),
    # errors
    ("too_few_arguments", "XR_crop_dataset", ["{'" + P1 + "'}"], "MATLAB:badsubscript"),
    ("no_arguments", "XR_crop_dataset", [], "MATLAB:badsubscript"),
    ("unknown_name", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "bachSize", "[1,2,3]"],
     "MATLAB:InputParser:UnmatchedParameter"),
    ("ambiguous_prefix", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "save", "true"],
     "MATLAB:InputParser:AmbiguousParameter"),
    ("missing_value", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "pad"],
     "MATLAB:InputParser:ParamMissingValue"),
    ("number_for_logical", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "pad", "1"],
     "MATLAB:InputParser:ArgumentFailedValidation"),
    ("capital_true", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "pad", "TRUE"],
     "MATLAB:InputParser:ArgumentFailedValidation"),
    ("text_for_cell", "XR_crop_dataset", ["{'" + P1 + "'}", BBOX, "channelPatterns", "CamA"],
     "MATLAB:InputParser:ArgumentFailedValidation"),
    ("no_such_function", "no_such_function", ["a", "b"], True),
]  # fmt: skip


# where pk5d-master differs from the compiled master (reported; none is a call of the callers)
KNOWN_DIFFERENCES = {
    "edge_hex_literal": "MATLAB's hex literal 0x100 is a uint16 and makes [512, 0x100] uint16; "
    "petakit5d.params.literal reads it, but the binding gives the function Python numbers, so "
    "the class is not kept (double in MATLAB's form)",
    "edge_nan_box_end": "a NaN box end: MATLAB passes it on, the bbox kind refuses it (a box "
    "end in the port is an integer or None)",
}

# calls of the callers that MATLAB rejects too
CALLER_ERRORS = {
    # the GUI's "bound box crop" of the pipeline passes boundboxCrop, a name the pipeline
    # function does not have (it has no crop option)
    "gui_XR_microscopeAutomaticProcessing_options": "MATLAB:InputParser:UnmatchedParameter",
}


def _cases():
    calls = []
    for name, function, argv in C.gui_calls():
        cid = f"gui_{function}_{name}"
        calls.append((cid, function, argv, CALLER_ERRORS.get(cid)))
    for name, function, argv in C.pypetakit5d_calls():
        cid = f"pyp_{function}_{name}"
        calls.append((cid, function, argv, CALLER_ERRORS.get(cid)))
    for name, function, argv, err in EDGE:
        calls.append((f"edge_{name}", function, argv, err))
    cases = []
    for cid, function, argv, err in calls:
        known = function in EXTRACT
        args = [ip_code(function) if known else "", site_names(function, argv) if known else []]
        cases.append(
            Case(
                cid,
                matlab="pkpy_j12_master",
                nargout=2,
                outputs=("values", "header"),
                args=[*args, function, *argv],
                python=py_master,
                expect_error=err or False,
                py_error=(MatlabArgumentError, master.MasterError),
                helpers=["tests/cli/matlab"],
                files=False,
            )
        )
    return cases


CASES = _cases()
FUNCTION = {c.id: c.args[2] for c in CASES}


@oracle.parametrize(CASES)
def test_master_binds_as_mccmaster(case, tmp_path, request):
    function = FUNCTION[case.id]
    if function in EXTRACT and function not in REGISTRY:
        pytest.skip(f"{function} is not ported yet")
    if case.id in KNOWN_DIFFERENCES:
        request.applymarker(pytest.mark.xfail(strict=True, reason=KNOWN_DIFFERENCES[case.id]))
    oracle.check(case, tmp_path)


def test_corpus_covers_the_phase1_callers():
    """Every Phase 1 function the GUI or PyPetaKit5D starts has calls in the corpus."""
    got = {FUNCTION[c.id] for c in CASES if c.id.startswith(("gui_", "pyp_"))}
    want = {
        "XR_crop_dataset", "XR_decon_data_wrapper", "XR_deskew_rotate_data_wrapper",
        "XR_generate_image_list_wrapper", "XR_matlab_stitching_wrapper",
        "XR_microscopeAutomaticProcessing", "XR_MIP_wrapper", "XR_resample_dataset",
        "XR_resave_zarr_wrapper", "XR_tiffToZarr_wrapper", "XR_zarrToTiff_wrapper",
    }  # fmt: skip
    assert want <= got


def test_ip_code_declares_every_parameter():
    code = ip_code("XR_crop_dataset")
    assert code.count("ip.addRequired(") == 2
    assert "ip.addParameter('batchSize', [1024, 1024, 1024], @isnumeric);" in code
