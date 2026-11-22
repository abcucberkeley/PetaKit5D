"""MATLAB's isvector in the parameter registry: 1-by-n or n-by-1 with n >= 0 (1x0 and 0x1 too),
never 0x0, a matrix or an n-D array.

The registry's shared predicate (``spec._pred("isvector")``) serves every validator text that
calls isvector (the lab's code: ``LAB`` lists the parameters) and the word ``vector`` of the
port's own checks of third-party code. Until 2026-11-21 it refused 1x0 and 0x1 arrays, which
MATLAB takes. MATLAB R2024b decides (oracle cases): ``isvector`` itself on 1x0, 0x1, 0x0, 1x1,
2x2 and 1x1x0 arrays, and inputParser with each lab validator that reaches the predicate, on
empty arrays of several classes and sizes, against the registry's check of each parameter that
has it. The checks follow inputParser only: a function may still fail further on with a value
they take, in MATLAB as in the port.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.params import MatlabArgumentError
from petakit5d.params import values as mv
from petakit5d.params.registry import REGISTRY
from petakit5d.params.spec import _pred
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, mcell, mempty
from petakit5d.testing.compare import MCell

# name -> (MATLAB's argument, the registry's form of it)
SHAPES = {
    "1x0": (mempty("double", (1, 0)), mv.MArray("double", (1, 0), [])),
    "0x1": (mempty("double", (0, 1)), mv.MArray("double", (0, 1), [])),
    "0x0": (None, mv.MArray("double", (0, 0), [])),
    "1x1": (1.0, mv.MArray("double", (1, 1), [1.0])),
    "2x2": (np.array([[1.0, 2.0], [3.0, 4.0]]), mv.MArray("double", (2, 2), [1.0, 3.0, 2.0, 4.0])),
    "1x1x0": (mempty("double", (1, 1, 0)), mv.MArray("double", (1, 1, 0), [])),
}
# the values the lab validators get
VALUES = {
    **{f"double_{k}": SHAPES[k] for k in ("1x0", "0x1", "0x0", "1x1", "2x2", "1x1x0")},
    "uint16_1x0": (mempty("uint16", (1, 0)), mv.MArray("uint16", (1, 0), [])),
    "logical_0x1": (mempty("logical", (0, 1)), mv.MArray("logical", (0, 1), [])),
    "cell_1x0": (MCell((1, 0), []), mv.Cell([], (1, 0))),
    "cell_0x0": (MCell((0, 0), []), mv.Cell([], (0, 0))),
}
VALUES_PORT = {k: p for k, (_, p) in VALUES.items()}

# the lab's validator texts that reach the predicate, and the parameters that have them
V_EMPTY_OR = "@(x) isempty(x) || isvector(x)"
V_ANON = "@(x) isvector(x)"
V_TWO = "@(x) isvector(x) && numel(x) <= 2"
V_NAMED = "@isvector"
LAB = [
    ("compute_tile_bwdist", "blockSize", "required", V_NAMED),
    ("compute_tile_bwdist", "shardSize", "required", V_EMPTY_OR),
    ("compute_tile_bwdist", "distBbox", "required", V_EMPTY_OR),
    ("compute_tile_bwdist", "dataOrderMat", "required", V_ANON),
    ("compute_tile_bwdist_mip_slabs", "blockSize", "required", V_NAMED),
    ("compute_tile_bwdist_mip_slabs", "shardSize", "required", V_EMPTY_OR),
    ("compute_tile_bwdist_mip_slabs", "poolSize", "required", V_NAMED),
    ("compute_tile_bwdist_mip_slabs", "dataOrderMat", "required", V_ANON),
    ("deskewPhases_data_wrapper", "ChunkSize", "parameter", V_NAMED),
    ("saveMIP_tiff", "inputBbox", "parameter", V_EMPTY_OR),
    ("saveMIP_zarr", "inputBbox", "parameter", V_EMPTY_OR),
    ("simRecon_data_wrapper", "ChunkSize", "parameter", V_NAMED),
    ("simReconAutomaticProcessing", "ChunkSize", "parameter", V_NAMED),
    ("simReconFrame", "ChunkSize", "parameter", V_NAMED),
    ("stitch_global_grid_assignment_wrapper", "axisWeight", "required", V_ANON),
    ("XR_chromatic_shift_correction_data_wrapper", "batchSize", "parameter", V_NAMED),
    ("XR_chromatic_shift_correction_data_wrapper", "blockSize", "parameter", V_NAMED),
    ("XR_crop_zarr", "batchSize", "parameter", V_NAMED),
    ("XR_decon_data_wrapper", "xyPixelSize", "parameter", V_TWO),
    ("XR_decon_data_wrapper", "batchSize", "parameter", V_NAMED),
    ("XR_deskew_rotate_data_wrapper", "xyPixelSize", "parameter", V_TWO),
    ("XR_deskew_rotate_data_wrapper", "batchSize", "parameter", V_NAMED),
    ("XR_deskew_rotate_data_wrapper", "blockSize", "parameter", V_NAMED),
    ("XR_deskew_rotate_data_wrapper", "inputBbox", "parameter", V_EMPTY_OR),
    ("XR_deskewRotateZarr", "batchSize", "parameter", V_NAMED),
    ("XR_deskewRotateZarr", "blockSize", "parameter", V_NAMED),
    ("XR_deskewRotateZarr", "inputBbox", "parameter", V_EMPTY_OR),
    ("XR_ExM_PunctaRemoval_block", "volThrsh", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemoval_interpolated_background_block", "volThrsh", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemoval_interpolated_background_block", "localWinSize", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemoval_updated_block", "volThrsh", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemoval_updated_block", "localWinSize", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemoval_zarr", "localWinSize", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemoval_zarr", "bbox", "parameter", V_EMPTY_OR),
    ("XR_ExM_PunctaRemoval_zarr", "batchSize", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemoval_zarr", "blockSize", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemovalPointDetection", "localWinSize", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemovalPointDetection", "SigmaThrsh", "parameter", V_NAMED),
    ("XR_ExM_PunctaRemovalPointDetection", "intThrsh", "parameter", V_NAMED),
    ("XR_matlab_stitching_wrapper", "xyMaxOffset", "parameter", V_TWO),
    ("XR_MIP_wrapper", "inputBbox", "parameter", V_EMPTY_OR),
    ("XR_MIP_wrapper", "batchSize", "parameter", V_NAMED),
    ("XR_MIP_zarr", "inputBbox", "parameter", V_EMPTY_OR),
    ("XR_psf_analysis_plot", "PSFsubpix", "required", V_NAMED),
    ("XR_unmix_channels_data_wrapper", "batchSize", "parameter", V_NAMED),
    ("XR_unmix_channels_data_wrapper", "blockSize", "parameter", V_NAMED),
    ("XR_unmix_channels_data_wrapper", "borderSize", "parameter", V_NAMED),
    ("XR_unmix_channels_frame", "blockSize", "parameter", V_NAMED),
    ("XR_unmix_channels_zarr", "batchSize", "parameter", V_NAMED),
    ("XR_unmix_channels_zarr", "blockSize", "parameter", V_NAMED),
    ("XR_unmix_channels_zarr", "borderSize", "parameter", V_NAMED),
]

# inputParser with one argument of the parameter's kind and the lab's validator, on each value
_ADD = {"required": "ip.addRequired('x', %s);", "parameter": "ip.addParameter('x', [], %s);"}
_GIVE = {"required": "v{i}", "parameter": "'x', v{i}"}


def _parse_code(kind: str, validator: str) -> str:
    return (
        "ip = inputParser; " + _ADD[kind] % validator + " v = in{1}; n = in{2}; r = struct(); "
        "for i = 1 : numel(v); try; ip.parse(" + _GIVE[kind] + "); r.(n{i}) = 'accepted'; "
        "catch ME; r.(n{i}) = ME.identifier; end; end; out = {r};"
    )


def isvector(x) -> bool:
    """the registry's isvector (the Python side of the cases)"""
    return _pred("isvector", x)


def outcomes(function: str, param: str, values: dict) -> dict[str, str]:
    """the registry's check of one parameter on each value: 'accepted' or the identifier"""
    p = REGISTRY.find(function).param(param)
    out = {}
    for name, v in values.items():
        try:
            p.check(v, function)
        except MatlabArgumentError as err:
            out[name] = err.identifier
        else:
            out[name] = "accepted"
    return out


def _lab_case(function: str, param: str, kind: str, validator: str) -> Case:
    return Case(
        f"lab.{function}.{param}",
        code=_parse_code(kind, validator),
        args=[mcell(*[m for m, _ in VALUES.values()]), mcell(*VALUES)],
        python=lambda *_: outcomes(function, param, VALUES_PORT),
    )


CASES = [
    Case(f"isvector.{name}", matlab="isvector", args=[m], python=isvector, py_args=[p])
    for name, (m, p) in SHAPES.items()
] + [_lab_case(*row) for row in LAB]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- without MATLAB -------------------------------------------------------------------------------


def test_lab_lists_every_parameter_whose_validator_reaches_isvector():
    """the table is the registry's: each parameter whose validator's class checks hold isvector,
    with its kind and validator text (a new one needs a row and its MATLAB reference)"""
    found = {
        (s.matlab_name, p.matlab): (p.kind, p.validator)
        for s in REGISTRY.specs()
        for p in s.params
        if any("isvector" in alt for alt in p.alternatives or [])
    }
    assert found == {(f, p): (k, v) for f, p, k, v in LAB}


@pytest.mark.parametrize(("function", "param", "kind", "validator"), LAB)
def test_what_each_lab_parameter_takes_now(function, param, kind, validator):
    """1x0 and 0x1 values of any class, as MATLAB's isvector; 0x0 only where the validator
    takes an empty value; never a matrix or an n-D array"""
    took = {k: v == "accepted" for k, v in outcomes(function, param, VALUES_PORT).items()}
    for name in ("double_1x0", "double_0x1", "uint16_1x0", "logical_0x1", "cell_1x0"):
        assert took[name], name
    assert took["double_1x1"]
    empty_ok = validator == V_EMPTY_OR
    assert took["double_0x0"] is empty_ok and took["cell_0x0"] is empty_ok
    assert took["double_1x1x0"] is empty_ok  # empty, and not a vector
    assert not took["double_2x2"]


# the port's forms of values: (value, MATLAB's isvector of it)
FORMS = {
    "numpy_1d_empty": (np.zeros(0), True),  # a 1-D array is a row
    "numpy_1x0": (np.zeros((1, 0)), True),
    "numpy_0x1": (np.zeros((0, 1)), True),
    "numpy_0x0": (np.zeros((0, 0)), False),
    "numpy_1x1x0": (np.zeros((1, 1, 0)), False),
    "numpy_0x1x1": (np.zeros((0, 1, 1)), True),  # MATLAB's size drops trailing singletons: 0x1
    "numpy_3x1x1": (np.zeros((3, 1, 1)), True),
    "numpy_1x3x2": (np.zeros((1, 3, 2)), False),
    "numpy_row": (np.arange(3.0), True),
    "numpy_column": (np.arange(3.0).reshape(3, 1), True),
    "numpy_scalar": (np.float64(1.0), True),
    "numpy_matrix": (np.zeros((2, 3)), False),
    "tuple": ((1.0, 2.0), True),
    "empty_tuple": ((), False),  # MATLAB's {} or []
    "none": (None, False),  # MATLAB's []
    "number": (2.0, True),
    "char_row": ("abc", True),
    "char_empty": ("", False),  # '' is 0x0
    "string_scalar": (mv.MString("abc"), True),
    "char_1x0": (mv.MArray("char", (1, 0), []), True),
    "cell_0x3": (mv.Cell([], (0, 3)), False),
    "cell_3x1": (mv.Cell([1.0, 2.0, 3.0], (3, 1)), True),
    "struct": ({"a": 1.0}, True),
}


@pytest.mark.parametrize("name", list(FORMS))
def test_isvector_of_the_ports_forms(name):
    value, want = FORMS[name]
    assert isvector(value) is want


def test_the_word_vector_is_the_shared_predicate():
    from petakit5d.params.spec import accepts_alternatives

    (alt,) = accepts_alternatives("vector")
    for name, (value, want) in FORMS.items():
        assert all(_pred(q, value) for q in alt) is want, name
    for name, (_, value) in SHAPES.items():
        assert all(_pred(q, value) for q in alt) is isvector(value), name
