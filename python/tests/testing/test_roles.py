"""Roles: MATLAB's arguments in Python's form for the Python call, Python's results in MATLAB's.

A case keeps its arguments as MATLAB takes them; ``Vol``, ``Box``, ``PerAxis``, ``Index1``,
``Range1``, ``Dim`` and ``LinIndex`` say what an argument means, ``out_roles`` what each result
is. Neither changes the case key. The comparison names the mapping's mistakes (a volume
without its role, a box converted twice, y and x swapped, a first index taken in C order) and
gives indices in both forms; zarr arrays of the two layouts compare voxel by voxel.
"""

import json

import numpy as np
import pytest

from petakit5d.core import arrays, axes, bbox
from petakit5d.testing import oracle
from petakit5d.testing.cases import (
    Box,
    Case,
    Dim,
    In,
    Index1,
    LinIndex,
    PerAxis,
    Range1,
    Value,
    Vol,
    Work,
    out_to_matlab,
)
from petakit5d.testing.compare import (
    DOUBLE_MAP,
    MISSING_MAP,
    YX_SWAP,
    Masks,
    compare,
    compare_trees,
    describe,
    describe_tree,
    from_descriptor,
    layout,
    scan_tree,
    to_mvalue,
    tree,
    tree_from_descriptor,
    vhash,
)


def mvol(shape, dtype=np.uint16, seed=0):
    """a MATLAB-form test volume (NumPy shape = MATLAB size), distinct deterministic values"""
    a = (np.arange(int(np.prod(shape)), dtype=np.int64) * 7919 + seed * 104729) % 65521
    return np.asfortranarray(a.reshape(shape, order="F").astype(dtype))


def fake_reference(case, outputs, source="store"):
    """MATLAB's result as the harness reads it: from the store (full values) or a manifest
    (descriptors: arrays over 2 KB by hash only)"""
    vals = [to_mvalue(x) for x in outputs]
    if source == "manifest":
        vals = [from_descriptor(describe(v)) for v in vals]
    return oracle.Reference(
        key=case.key(),
        status="ok",
        error=None,
        warnings=[],
        outputs=vals,
        files=None,
        meta={},
        source=source,
        masks=Masks(),
        commit="0" * 40,
    )


def run(case, matlab_outputs, tmp_path, source="store"):
    res = case.run_python(tmp_path)
    assert res.error is None, res.error_text
    return oracle.compare_case(case, fake_reference(case, matlab_outputs, source), res)


# --- arguments ----------------------------------------------------------------------------------


def seen(**kw):
    """a case whose Python side records what it was called with"""
    got = {}

    def fn(*args, **kwargs):
        got["args"], got["kwargs"] = args, kwargs
        return 0

    c = Case("seen", matlab="disp", python=fn, **kw)
    return c, got


def test_the_python_call_gets_each_role_in_python_form(tmp_path):
    m = mvol((4, 5, 3))
    c, got = seen(
        args=[
            Vol(In("im")),
            Box([2, 3, 1, 4, 5, 3]),
            PerAxis([0.1, 0.2, 0.5]),
            PerAxis([7, 8, 9], "xyz"),
            Index1([1, 4]),
            Range1([2, np.inf]),
            Dim(3),
            Dim(1, 3),
            LinIndex(6, In("im")),
            LinIndex([1, 60], [4, 5, 3]),
        ],
        inputs={"im": Value(lambda: m)},
    )
    c.run_python(tmp_path)
    a = got["args"]
    np.testing.assert_array_equal(a[0], arrays.from_matlab(m))
    assert a[0].flags.c_contiguous and a[0].shape == (3, 4, 5)
    assert a[1] == (0, 1, 2, 3, 4, 5)
    assert list(a[2]) == [0.5, 0.1, 0.2]
    assert list(a[3]) == [9, 8, 7]
    assert list(a[4]) == [0, 3]
    assert a[5] == (1, None)
    assert a[6] == -3 and a[7] == 1
    assert a[8] == (0, 1, 1)  # MATLAB's 6th voxel of a 4x5x3 array: y 2, x 2, z 1
    assert [tuple(map(int, t)) for t in zip(*a[9], strict=True)] == [(0, 0, 0), (2, 3, 4)]


def test_roles_in_py_args_and_py_kwargs(tmp_path):
    c, got = seen(
        args=[In("im"), "bbox", [1, 1, 1, 2, 2, 2]],
        py_args=[Vol(In("im"))],
        py_kwargs={"bbox": Box([1, 1, 1, 2, 2, 2]), "axis": Dim(2)},
        inputs={"im": Value(lambda: mvol((2, 3, 4)))},
    )
    c.run_python(tmp_path)
    assert got["args"][0].shape == (4, 2, 3)
    assert got["kwargs"] == {"bbox": (0, 0, 0, 2, 2, 2), "axis": -1}


def test_roles_change_neither_the_key_nor_matlabs_arguments():
    inputs = {"im": Value(lambda: mvol((4, 5, 3)))}
    plain = Case("k", matlab="crop3d_mex", args=[In("im"), [2, 3, 1, 4, 5, 3], 3], inputs=inputs)
    roles = Case(
        "k",
        matlab="crop3d_mex",
        args=[Vol(In("im")), Box([2, 3, 1, 4, 5, 3]), Dim(3)],
        inputs=inputs,
        out_roles=("vol",),
        python=lambda a, b, d: a,
    )
    assert plain.key() == roles.key()
    path = lambda name, x: f"/in/{name}"  # noqa: E731
    assert vhash(plain.matlab_args(path)) == vhash(roles.matlab_args(path))


def test_out_roles_map_each_kind_back():
    p = arrays.from_matlab(mvol((4, 5, 3)))
    np.testing.assert_array_equal(out_to_matlab("vol", p), mvol((4, 5, 3)))
    assert out_to_matlab("box", (0, 1, 2, 3, 4, 5)) == (2, 3, 1, 4, 5, 3)
    assert out_to_matlab("yxz", (3, 1, 2)) == (1, 2, 3)
    assert out_to_matlab("xyz", [3, 1, 2]) == [2, 1, 3]
    assert out_to_matlab("index1", [0, 4]) == [1, 5]
    assert out_to_matlab("range1", (0, None)) == (1, np.inf)
    assert out_to_matlab(("dim", 3), -3) == 3
    assert out_to_matlab(("linindex", (4, 5, 3)), (0, 1, 1)) == 6
    assert out_to_matlab(["index1", "yxz"], (2, 0, 1)) == (1, 2, 3)  # in order
    two = out_to_matlab(("each", "vol"), [p, p])
    assert len(two) == 2 and two[1].shape == (4, 5, 3)
    assert out_to_matlab({"im": "vol", "n": None}, {"im": p, "n": 2})["im"].shape == (4, 5, 3)
    assert out_to_matlab(lambda v: v + 1, 1) == 2
    with pytest.raises(ValueError):
        out_to_matlab("volume", p)


def test_out_roles_must_match_nargout():
    with pytest.raises(ValueError, match="out roles"):
        Case("n", matlab="disp", nargout=2, out_roles=("vol",))


# --- acceptance: a case with roles equals the case converted by hand --------------------------


def py_crop_sum(vol, box, axis):
    """a crop of a (z, y, x) volume and its cumulative sum along an axis"""
    c = vol[bbox.box_slices(box)]
    return c, np.cumsum(c, axis=axis, dtype=np.float64)


def test_a_case_with_roles_equals_the_case_converted_by_hand(tmp_path):
    m = mvol((7, 6, 5), np.uint16, 3)
    inputs = {"im": Value(lambda: m)}
    box, dim = [2, 1, 2, 6, 5, 4], 2

    def hand(im, b, d):  # the conversions written out, as cases did before roles
        out = py_crop_sum(arrays.from_matlab(im), bbox.box_from_matlab(b), axes.dim_to_axis(d, 3))
        return tuple(arrays.to_matlab(x) for x in out)

    by_hand = Case(
        "hand",
        matlab="pkpy_oracle.selftest.identity",
        args=[In("im"), box, dim],
        nargout=2,
        python=hand,
        inputs=inputs,
    )
    with_roles = Case(
        "roles",
        matlab="pkpy_oracle.selftest.identity",
        args=[Vol(In("im")), Box(box), Dim(dim, 3)],
        nargout=2,
        python=py_crop_sum,
        out_roles=("vol", "vol"),
        inputs=inputs,
    )
    assert by_hand.key() == with_roles.key()
    a = by_hand.run_python(tmp_path / "a").outputs
    b = with_roles.run_python(tmp_path / "b").outputs
    assert [vhash(to_mvalue(x)) for x in a] == [vhash(to_mvalue(x)) for x in b]
    # and both agree with MATLAB's result (made here in MATLAB's form)
    crop = m[1:6, 0:5, 1:4]
    want = [crop, np.cumsum(crop, axis=1, dtype=np.float64)]
    for case in (by_hand, with_roles):
        for source in ("store", "manifest"):
            rep = run(case, want, tmp_path / f"c_{case.id}_{source}", source)
            assert rep.ok, rep.text()


# --- acceptance: seeded mistakes give named hints ------------------------------------------------


def times_two(v):
    return v * 2


def test_a_volume_without_its_role_is_named(tmp_path):
    m = mvol((4, 5, 3), np.float64, 1)
    want = [m * 2]
    no_in_role = Case(
        "no_in_role",
        matlab="times",
        args=[In("x"), 2],
        py_args=[In("x")],
        python=times_two,
        out_roles=("vol",),
        inputs={"x": Value(lambda: m)},
    )
    for source in ("store", "manifest"):
        text = run(no_in_role, want, tmp_path / f"a_{source}", source).text()
        assert "without its Vol role" in text, text
    no_out_role = Case(
        "no_out_role",
        matlab="times",
        args=[Vol(In("x")), 2],
        py_args=[Vol(In("x"))],
        python=times_two,
        inputs={"x": Value(lambda: m)},
    )
    text = run(no_out_role, want, tmp_path / "b").text()
    assert MISSING_MAP in text, text
    # a cube: the sizes agree, the values say it
    c = mvol((4, 4, 4), np.float64, 2)
    cube = Case(
        "cube",
        matlab="times",
        args=[Vol(In("x")), 2],
        py_args=[Vol(In("x"))],
        python=times_two,
        inputs={"x": Value(lambda: c)},
    )
    text = run(cube, [c * 2], tmp_path / "c").text()
    assert MISSING_MAP in text, text


def test_a_box_converted_twice_is_named(tmp_path):
    def grow(box):  # a Python box one voxel larger at the end
        z0, y0, x0, z1, y1, x1 = box
        return (z0, y0, x0, z1 + 1, y1 + 1, x1 + 1)

    twice = Case(
        "box_twice",
        matlab="plus",
        args=[Box([2, 3, 1, 9, 10, 7]), [0, 0, 0, 1, 1, 1]],
        py_args=[Box([2, 3, 1, 9, 10, 7])],
        python=grow,
        py_out=lambda b: [bbox.box_to_matlab(b)],  # and the out role converts again
        out_roles=("box",),
    )
    text = run(twice, [np.array([[2.0, 3, 1, 10, 11, 8]])], tmp_path).text()
    assert "a box converted twice" in text, text


def test_y_and_x_swapped_is_named(tmp_path):
    m = mvol((4, 6, 3), np.float64, 4)
    swapped = Case(
        "yx",
        matlab="times",
        args=[Vol(In("x")), 2],
        py_args=[Vol(In("x"))],
        python=lambda v: np.swapaxes(v, -1, -2) * 2,
        out_roles=("vol",),
        inputs={"x": Value(lambda: m)},
    )
    text = run(swapped, [m * 2], tmp_path / "a").text()
    assert YX_SWAP in text, text
    vec = Case(
        "yx_vec",
        matlab="size",
        args=[In("x")],
        py_args=[Vol(In("x"))],
        python=lambda v: v.shape,
        out_roles=("xyz",),  # MATLAB's size is [y x z]
        inputs={"x": Value(lambda: m)},
    )
    text = run(vec, [np.array([[4.0, 6, 3]])], tmp_path / "b").text()
    assert "y and x swapped" in text, text


def test_a_first_maximum_in_c_order_is_named(tmp_path):
    m = np.zeros((4, 5, 3), np.uint16, order="F")
    m[2, 1, 0] = m[1, 3, 0] = m[0, 0, 2] = 9  # ties: MATLAB's first is (3,2,1), linear 7
    first = Case(
        "argmax",
        code="[~, k] = max(in{1}(:)); out = {k};",
        args=[In("x")],
        py_args=[Vol(In("x"))],
        python=lambda v: np.unravel_index(int(np.argmax(v)), v.shape),
        out_roles=(("linindex", "x"),),
        inputs={"x": Value(lambda: m)},
    )
    text = run(first, [7.0], tmp_path).text()
    assert "a first-index rule (first maximum, first nonzero) taken in C order" in text, text


# --- diff reports give both index forms -----------------------------------------------------------


def test_diff_reports_give_matlab_subscripts_and_zyx_indices():
    r = mvol((4, 5, 3))
    a = r.copy(order="F")
    a[2, 3, 1] += 1  # MATLAB (3,4,2) = (z, y, x) (1, 2, 3)
    text = compare(to_mvalue(r), a).text()
    assert "MATLAB (3,4,2) = (z, y, x) (1, 2, 3)" in text, text
    img = mvol((4, 5))
    b = img.copy(order="F")
    b[1, 2] += 1
    assert "(y, x)" not in compare(to_mvalue(img), b).text()  # an image: the same indices
    with layout("vol"):
        assert "(2,3) = (y, x) (1, 2)" in compare(to_mvalue(img), b).text()


def test_size_hints_name_the_mapping_without_the_data():
    r = to_mvalue(mvol((40, 50, 30)))  # over 2 KB: a manifest holds its hash only
    h = from_descriptor(describe(r))
    assert h.data is None
    zyx = np.ascontiguousarray(arrays.from_matlab(r.data))
    text = compare(h, zyx).text()
    assert MISSING_MAP in text, text
    text = compare(h, arrays.to_matlab(r.data)).text()
    assert DOUBLE_MAP in text, text


# --- zarr arrays of both layouts ------------------------------------------------------------------


cppzarr = pytest.importorskip("cppzarr")


def write_matlab_zarr(path, m):
    """a zarr as MATLAB's writer makes it: order F, shape and chunks (y, x, z)"""
    cppzarr.write_zarr(str(path), m, order="F", chunks=[4, 4, 2])


def write_python_zarr(path, p):
    """a native zarr: order C, shape and chunks (z, y, x)"""
    cppzarr.write_zarr(str(path), p, order="C", chunks=[2, 4, 4])


def trees(tmp_path, m, p):
    (tmp_path / "m").mkdir()
    (tmp_path / "p").mkdir()
    write_matlab_zarr(tmp_path / "m" / "a.zarr", m)
    write_python_zarr(tmp_path / "p" / "a.zarr", p)
    return scan_tree(tmp_path / "m", ["a.zarr"]), scan_tree(tmp_path / "p", ["a.zarr"])


def test_matlab_and_python_zarr_of_the_same_data_compare_equal(tmp_path):
    m = mvol((8, 8, 5), np.uint16, 6)
    ref, act = trees(tmp_path, m, arrays.from_matlab(m))
    assert json.loads((tmp_path / "m/a.zarr/.zarray").read_text())["order"] == "F"
    assert json.loads((tmp_path / "p/a.zarr/.zarray").read_text())["shape"] == [5, 8, 8]
    spec = tree(zarr_chunks="strict")
    rep = compare_trees(ref, act, spec)
    assert rep.ok, rep.text()
    assert any("compared after the mapping" in n for n in rep.notes)
    # as CI reads it: MATLAB's tree from its descriptor (voxels by hash)
    d = json.loads(json.dumps(describe_tree(ref, spec, Masks())))
    rep = compare_trees(tree_from_descriptor(d), act, spec)
    assert rep.ok, rep.text()


def test_one_transposed_plane_of_a_python_zarr_is_found(tmp_path):
    m = mvol((8, 8, 5), np.uint16, 7)
    p = arrays.from_matlab(m)
    p[3] = p[3].T.copy()  # z plane 3 with y and x exchanged
    ref, act = trees(tmp_path, m, p)
    text = compare_trees(ref, act).text()
    assert "values" in text and "(z, y, x) (3," in text, text


def test_a_python_zarr_with_matlabs_shape_is_named(tmp_path):
    m = mvol((8, 6, 5), np.uint16, 8)
    (tmp_path / "m").mkdir()
    (tmp_path / "p").mkdir()
    write_matlab_zarr(tmp_path / "m" / "a.zarr", m)
    cppzarr.write_zarr(str(tmp_path / "p" / "a.zarr"), np.ascontiguousarray(m), order="C")
    rep = compare_trees(
        scan_tree(tmp_path / "m", ["a.zarr"]), scan_tree(tmp_path / "p", ["a.zarr"])
    )
    assert "written without its mapping" in rep.text(), rep.text()


def test_a_python_zarr_in_matlabs_format_compares_as_before(tmp_path):
    m = mvol((8, 6, 5), np.uint16, 9)
    (tmp_path / "m").mkdir()
    (tmp_path / "p").mkdir()
    write_matlab_zarr(tmp_path / "m" / "a.zarr", m)
    write_matlab_zarr(tmp_path / "p" / "a.zarr", m)
    rep = compare_trees(
        scan_tree(tmp_path / "m", ["a.zarr"]), scan_tree(tmp_path / "p", ["a.zarr"])
    )
    assert rep.ok and not any("mapping" in n for n in rep.notes), rep.text()


def test_a_case_writing_a_native_zarr_against_matlabs(tmp_path):
    """the whole path: MATLAB's reference tree with its zarr, Python's case writing its own"""
    m = mvol((8, 8, 4), np.uint16, 10)
    mt = tmp_path / "matlab"
    mt.mkdir()
    write_matlab_zarr(mt / "out.zarr", m)
    case = Case(
        "zarr_layouts",
        matlab="writezarr",
        args=[Vol(In("x")), Work("out.zarr")],
        python=lambda p, path: write_python_zarr(path.rstrip("/"), p),
        nargout=0,
        inputs={"x": Value(lambda: m)},
        tree=tree(zarr_chunks="strict"),
    )
    res = case.run_python(tmp_path / "py")
    assert res.error is None, res.error_text
    ref = fake_reference(case, [])
    ref.files = scan_tree(mt, ["out.zarr"])
    rep = oracle.compare_case(case, ref, res)
    assert rep.ok, rep.text()


def test_a_volume_output_gives_python_indices_also_for_an_image(tmp_path):
    img = mvol((5, 6), np.float64, 5)
    c = Case(
        "img",
        matlab="times",
        args=[Vol(In("x")), 2],
        py_args=[Vol(In("x"))],
        python=lambda v: v * 2 + (np.arange(v.size).reshape(v.shape) == 7),
        out_roles=("vol",),
        inputs={"x": Value(lambda: img)},
    )
    text = run(c, [img * 2], tmp_path).text()
    assert "MATLAB (2,2) = (y, x) (1, 1)" in text, text
