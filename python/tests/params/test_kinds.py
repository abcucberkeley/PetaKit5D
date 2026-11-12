"""Conversion kinds: MATLAB value -> Python value -> MATLAB value is the identity, and back."""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.params import MatlabArgumentError, Param, Registry
from petakit5d.params.kinds import KINDS, get_kind, looks_converted

INF = math.inf

# (kind, MATLAB values, Python values they become)
CASES = {
    "same": [((1, 2, 3), (1, 2, 3)), ("x,y,z", "x,y,z"), (None, None)],
    "bbox": [
        ((1, 2, 3, 10, 20, 30), (2, 0, 1, 30, 10, 20)),
        ((1, 1, 1, INF, 20, INF), (0, 0, 0, None, None, 20)),
        ((5, 6, 7, 8), (4, 5, 7, 8)),  # 2D [ymin xmin ymax xmax]
        ((1, 1, INF, INF), (0, 0, None, None)),
        (None, None),
    ],
    "bboxes": [
        (((1, 2, 3, 10, 20, 30), (4, 5, 6, 7, 8, 9)), ((2, 0, 1, 30, 10, 20), (5, 3, 4, 9, 7, 8))),
        (((1, 1, 2, 2),), ((0, 0, 2, 2),)),
        (None, None),
    ],
    "size_yxz": [
        ((512, 256, 64), (64, 512, 256)),
        ((1, 2, 3, 4, 5, 6), (3, 1, 2, 6, 4, 5)),  # two pool levels
        ((1, 2, 3, 4, 5, 6, 7, 8, 9), (3, 1, 2, 6, 4, 5, 9, 7, 8)),  # three
        ((100, 200), (100, 200)),  # 2D [y x]
        (500, 500),
        (None, None),
    ],
    "per_axis_xyz": [((1, 2, 3), (3, 2, 1)), ((1, 2), (2, 1)), (4, 4), (None, None)],
    "resample_yxz": [((1, 2, 3), (3, 1, 2)), ((2, 5), (5, 2)), (2, 2), (None, None)],
    "resample_xyz": [((1, 2, 3), (3, 2, 1)), ((2, 5), (5, 2)), (0.5, 0.5), (None, None)],
    "flags_yxz": [((0, 0, 1), (1, 0, 0)), ((True, False, False), (False, True, False))],
    "index1": [(1, 0), ((1, 5, 9), (0, 4, 8)), (((1, 2), (3, 4)), ((0, 1), (2, 3))), (None, None)],
    "point1_yxz": [((5, 6, 7), (6, 4, 5)), ((5, 6), (4, 5)), (None, None)],
    "range1": [((1, 10), (0, 10)), ((3, INF), (2, None)), (None, None)],
    "dim": [(1, -2), (2, -1), (3, -3), (4, -4), ((3, 1), (-3, -2))],
    "axis_string": [("yxz", "zyx"), ("xzy", "yxz"), ("", "")],
    "data_order": [("y,x,z", "z,y,x"), ("x,y,z", "z,x,y"), ("", "")],
}


@pytest.mark.parametrize(
    ("kind", "m", "py"),
    [(k, m, p) for k, rows in CASES.items() for m, p in rows],
    ids=[f"{k}-{i}" for k, rows in CASES.items() for i in range(len(rows))],
)
def test_kind_round_trip(kind, m, py):
    k = KINDS[kind]
    got = k.from_matlab(m)
    assert got == py and type(got) is type(py)
    back = k.to_matlab(got)
    assert back == m
    assert k.from_matlab(k.to_matlab(py)) == py


def test_every_kind_is_tested():
    assert set(CASES) | {"volume"} == set(KINDS)


def test_volume_kind():
    m = np.arange(2 * 3 * 4, dtype=np.uint16).reshape((2, 3, 4), order="F")  # (y, x, z)
    py = KINDS["volume"].from_matlab(m)
    assert py.shape == (4, 2, 3) and py.flags.c_contiguous
    assert py[1, 0, 2] == m[0, 2, 1]
    assert np.array_equal(KINDS["volume"].to_matlab(py), m)
    assert KINDS["volume"].from_matlab(None) is None


def test_numpy_inputs_keep_their_type():
    v = KINDS["size_yxz"].from_matlab(np.array([512, 256, 64]))
    assert isinstance(v, np.ndarray) and v.tolist() == [64, 512, 256]
    v = KINDS["index1"].from_matlab(np.array([[1, 2], [3, 4]]))
    assert v.tolist() == [[0, 1], [2, 3]]


@pytest.mark.parametrize(
    ("kind", "bad"),
    [("bbox", (1, 2, 3)), ("size_yxz", (1, 2, 3, 4)), ("range1", (1, 2, 3)), ("index1", 1.5)],
)
def test_bad_values_are_errors(kind, bad):
    with pytest.raises(ValueError):
        KINDS[kind].from_matlab(bad)


def test_custom_pair_and_unknown_kind():
    k = get_kind((lambda v: v * 2, lambda v: v // 2))
    assert k.from_matlab(3) == 6 and k.to_matlab(6) == 3
    with pytest.raises(ValueError, match="unknown conversion kind"):
        get_kind("boxx")


@pytest.mark.parametrize(
    ("name", "looks"),
    [
        ("inputBbox", True), ("batchSize", True), ("BatchBBoxes", True), ("dataOrder", True),
        ("inputAxisOrder", True), ("resampleFactor", True), ("mipAxis", True), ("axis", True),
        ("batchInds", True), ("poolSize", True), ("xcorrDownsample", True), ("vol", True),
        ("xyPixelSize", False), ("dz", False), ("skewAngle", False), ("xyMaxOffset", False),
        ("parseCluster", False), ("cpusPerTask", False), ("jobLogDir", False),
    ],
)  # fmt: skip
def test_looks_converted(name, looks):
    assert looks_converted(name) is looks


def test_phase1_parameters_have_planned_kinds(extract):
    """Every Phase 1 parameter whose name looks like a box, size, axis or array has a planned
    kind (or names the area task that decides it), and the planned defaults convert."""
    from petakit5d.params.literal import LiteralError, parse
    from petakit5d.params.planned import PLANNED, planned_kind
    from petakit5d.params.values import to_python

    missing, bad = [], []
    for f, _, _, phase in PLANNED:
        if phase != 1 or f not in extract:
            continue
        for p in extract[f].params:
            k = planned_kind(f, p.name)
            if k is None:
                if looks_converted(p.name):
                    missing.append(f"{f}.{p.name}")
                continue
            assert k.startswith("?") or k in KINDS, (f, p.name, k)
            if k in KINDS and p.default and p.kind != "required":
                try:
                    m = to_python(parse(p.default))
                except LiteralError:
                    continue  # a computed or site-specific default
                try:
                    v = KINDS[k].from_matlab(m)
                    assert KINDS[k].to_matlab(v) is not None or v is None
                except Exception as err:  # noqa: BLE001
                    bad.append(f"{f}.{p.name}={p.default}: {err}")
    assert missing == [] and bad == []


def _registry():
    reg = Registry()

    @reg.register(
        "XR_kind_fake",
        Param("data", "@isnumeric", convert="volume"),
        Param("inputBbox", "@isnumeric", convert="bbox"),
        Param("batchSize", "@isnumeric", convert="size_yxz"),
        Param("mipAxis", "@isnumeric", convert="flags_yxz"),
        Param("inputAxisOrder", "@ischar", convert="axis_string"),
        Param("dim", "@isnumeric", convert="dim"),
    )
    def kind_fake(
        data, input_bbox, *, batch_size=(256, 1024, 512), mip_axis=(1, 0, 0),
        input_axis_order="zyx", dim=-3,
    ):  # fmt: skip
        return locals()

    return reg.find("XR_kind_fake")


def test_matlab_layer_converts_every_kind():
    spec = _registry()
    m = np.zeros((3, 4, 5), order="F")
    got = spec.call_matlab(m, [1, 2, 3, 3, 4, 5], "batchSize", [8, 16, 4], "mipAxis", [1, 0, 0],
                           "inputAxisOrder", "xzy", "dim", 1)  # fmt: skip
    assert got["data"].shape == (5, 3, 4)
    assert got["input_bbox"] == (2, 0, 1, 5, 3, 4)
    assert got["batch_size"] == (4, 8, 16)
    assert got["mip_axis"] == (0, 1, 0)
    assert got["input_axis_order"] == "yxz"
    assert got["dim"] == -2
    # the defaults are MATLAB's after conversion
    texts = {p.matlab: p.default_text() for p in spec.params if not p.required}
    assert texts == {
        "batchSize": "[1024,512,256]",
        "mipAxis": "[0,0,1]",
        "inputAxisOrder": "'yxz'",
        "dim": "3",
    }
    res = spec.matlab_results({**got, "data": None})
    assert res["inputBbox"].data == [1.0, 2.0, 3.0, 3.0, 4.0, 5.0]
    assert res["batchSize"].data == [8.0, 16.0, 4.0]
    assert res["inputAxisOrder"].text() == "xzy"
    assert res["dim"].data == [1.0]


def test_bad_values_in_the_matlab_layer():
    spec = _registry()
    with pytest.raises(MatlabArgumentError) as e:
        spec.call_matlab(np.zeros((2, 2, 2)), [1, 2, 3])
    assert e.value.identifier == "PETAKIT5D:params:badValue"
    assert "inputBbox" in str(e.value)


def test_function_level_hooks():
    reg = Registry()
    seen = {}

    def m_in(kw):
        seen["in"] = dict(kw)
        return {**kw, "n": kw["n"] * 10}

    def m_out(r):
        return ("matlab", r)

    @reg.register("XR_hook_fake", Param("n", "@isnumeric"), matlab_in=m_in, matlab_out=m_out)
    def hook_fake(n):
        return n + 1

    spec = reg.find("XR_hook_fake")
    assert spec.call_matlab(2) == ("matlab", 21)
    assert spec.call_strings(["3"]) == ("matlab", 31)
    assert hook_fake(2) == 3  # the snake_case function is not converted
