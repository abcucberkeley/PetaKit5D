"""Unit tests of petakit5d.testing.compare (no MATLAB needed)."""

import hashlib
import json
import struct

import numpy as np
import pytest

from petakit5d.testing import compare as C
from petakit5d.testing.compare import (
    MArray,
    Masks,
    MCell,
    MMap,
    MObject,
    MString,
    MStruct,
    compare,
    compare_trees,
    describe,
    describe_tree,
    from_descriptor,
    mchar,
    msize,
    read_mat,
    read_pkv,
    scan_tree,
    subscripts,
    to_mvalue,
    tol,
    tree,
    tree_from_descriptor,
    vhash,
    write_pkv,
)


def from_bits(b: int) -> np.float64:
    return np.array([b], np.uint64).view(np.float64)[0]


MATLAB_NAN = from_bits(0xFFF8000000000000)  # what MATLAB's NaN and 0/0 are


def vol(shape=(10, 12, 8), dtype=np.uint16, seed=0):
    a = (np.arange(int(np.prod(shape)), dtype=np.int64) * 7919 + seed * 104729) % 4001
    return np.asfortranarray(a.reshape(shape, order="F").astype(dtype))


def report(ref, act, spec=C.EXACT):
    rep = compare(to_mvalue(ref), act, spec, "out1")
    return rep, rep.text("case")


# --- sizes and conversion ---------------------------------------------------------------------


def test_msize_drops_trailing_singletons_after_the_second():
    assert msize((5, 7, 1)) == (5, 7)
    assert msize((1, 1, 1, 1)) == (1, 1)
    assert msize((0, 3, 1, 2)) == (0, 3, 1, 2)
    with pytest.raises(ValueError):
        msize((5,))


def test_subscripts_are_matlab_column_major():
    assert subscripts(0, (10, 12, 8)) == (1, 1, 1)
    assert subscripts(517, (10, 12, 8)) == (8, 4, 5)


@pytest.mark.parametrize(
    "x, cls, size",
    [
        (None, "double", (0, 0)),
        (True, "logical", (1, 1)),
        (3, "double", (1, 1)),
        (3.5, "double", (1, 1)),
        (1 + 2j, "double", (1, 1)),
        ("abc", "char", (1, 3)),
        ("", "char", (0, 0)),
        (b"ab", "uint8", (1, 2)),
        ([1, 2, 3], "double", (1, 3)),
        ([True, False], "logical", (1, 2)),
        ([], "double", (0, 0)),
        (np.float32(2), "single", (1, 1)),
        (np.zeros(4, np.int8), "int8", (1, 4)),
        (np.zeros((5, 7, 1), np.uint16), "uint16", (5, 7)),
        (np.zeros((2, 3, 4, 5)), "double", (2, 3, 4, 5)),
        (np.array(True), "logical", (1, 1)),
        (np.array([["a", "b"], ["c", "d"]]), "char", (2, 2)),
    ],
)
def test_python_values_as_matlab_arguments(x, cls, size):
    v = to_mvalue(x)
    assert isinstance(v, MArray)
    assert (v.mclass, v.size) == (cls, size)


def test_lists_dicts_and_text():
    c = to_mvalue(["CamA", "CamB"])
    assert isinstance(c, MCell) and c.size == (1, 2) and c.items[1].text() == "CamB"
    m = to_mvalue([1, "a"])
    assert isinstance(m, MCell)
    s = to_mvalue({"a": 1, "b": [1, 2]})
    assert isinstance(s, MStruct) and s.fields == ["a", "b"] and s.items[0]["b"].size == (1, 2)
    assert to_mvalue("\U0001f600").size == (1, 2)  # UTF-16: a surrogate pair, as MATLAB's char
    assert np.asarray(to_mvalue(np.arange(6).reshape(2, 3)).data)[1, 2] == 5


def test_reading_python_results_as_the_reference():
    col = to_mvalue(np.zeros((4, 1)))
    assert to_mvalue(np.zeros(4), like=col).size == (4, 1)
    assert to_mvalue(np.zeros(4)).size == (1, 4)
    u16 = to_mvalue(np.zeros((1, 1), np.uint16))
    assert to_mvalue(7, like=u16).mclass == "uint16"
    assert to_mvalue(7).mclass == "double"
    empty_row = MArray("char", np.zeros((1, 0), "<u2"))
    assert to_mvalue("", like=empty_row).size == (1, 0)
    ec = MCell((0, 3), [])
    assert to_mvalue([], like=ec).size == (0, 3)
    sa = MStruct((1, 2), ["a"], [{"a": to_mvalue(1)}, {"a": to_mvalue(2)}])
    got = to_mvalue([{"a": 1}, {"a": 2}], like=sa)
    assert isinstance(got, MStruct) and got.size == (1, 2)
    mp = MMap("char", "any", to_mvalue(["a", "b"]), MCell((1, 2), [to_mvalue(1), to_mvalue(2)]))
    assert isinstance(to_mvalue({"b": 2, "a": 1}, like=mp), MMap)
    assert isinstance(to_mvalue("x", like=MString((1, 1), ["y"])), MString)


def test_values_without_matlab_counterpart_are_rejected():
    with pytest.raises(TypeError):
        to_mvalue(np.zeros(3, np.float16))
    with pytest.raises(ValueError):
        to_mvalue(2**60)
    with pytest.raises(TypeError):
        to_mvalue({1, 2})


# --- value files and hashes -------------------------------------------------------------------


def special_values():
    nd = np.asfortranarray(np.arange(120, dtype=np.float32).reshape(2, 3, 4, 5))
    return [
        to_mvalue(
            np.array(
                [MATLAB_NAN, np.nan, -0.0, np.inf, -np.inf, 5e-324, from_bits(0x7FF0000000000123)]
            )
        ),
        to_mvalue(np.array([np.iinfo(np.int64).min, np.iinfo(np.int64).max, -1], np.int64)),
        to_mvalue(np.array([2**64 - 1, 0], np.uint64)),
        to_mvalue(nd),
        to_mvalue(np.zeros((0, 3), np.uint16)),
        to_mvalue(np.zeros((1, 0, 2), np.single)),
        to_mvalue(np.array([[complex(1.5, 0.0), complex(0, -0.0)]], np.complex64)),
        to_mvalue(np.array([[True, False], [False, True]])),
        mchar("héllo \U0001f600"),
        mchar(""),
        to_mvalue({"a": [1, 2], "b": {"c": "x"}, "e": None}),
        MStruct(
            (2, 1),
            ["p", "q"],
            [
                {"p": to_mvalue(1), "q": mchar("a")},
                {"p": to_mvalue("two"), "q": MCell((1, 1), [to_mvalue(3)])},
            ],
        ),
        MStruct((0, 0), ["p"], []),
        MStruct((1, 1), [], [{}]),
        MCell((2, 2), [to_mvalue(i) for i in range(4)]),
        MCell((0, 3), []),
        MString((1, 3), ["a", None, ""]),
        C.MSparse("double", (3, 4), [0, 2], [1, 1], np.array([1.5, -2.0])),
        C.MFunc("@(x)x+1"),
        MObject("inputParser", to_mvalue({"alpha": 5}), ""),
        C.MPlaceholder("work", "", "out/"),
    ]


@pytest.mark.parametrize("i", range(len(special_values())))
def test_value_files_round_trip_bit_for_bit(tmp_path, i):
    v = special_values()[i]
    p = tmp_path / "v.pkv"
    write_pkv(p, v)
    back = read_pkv(p)
    if isinstance(v, C.MPlaceholder):
        assert (back.what, back.sub) == ("work", "out/")
        return
    assert vhash(back) == vhash(v)
    assert compare(v, back).ok


def test_value_file_layout_is_aligned_with_a_json_footer(tmp_path):
    p = tmp_path / "v.pkv"
    write_pkv(p, MCell((1, 2), [mchar("ab"), to_mvalue(np.arange(3, dtype=np.uint16))]))
    b = p.read_bytes()
    assert b[:8] == C.MAGIC and b[-8:] == C.END_MAGIC
    n = struct.unpack("<Q", b[-16:-8])[0]
    hdr = json.loads(b[-16 - n : -16])
    second = hdr["root"]["v"][1]
    assert second["o"] % 8 == 0 and second["n"] == 6
    assert np.frombuffer(b[second["o"] : second["o"] + 6], "<u2").tolist() == [0, 1, 2]


def test_hash_follows_the_documented_definition():
    a = np.array([[1, 2, 3], [4, 5, 6]], np.uint16)
    raw = np.asfortranarray(a).tobytes(order="F")
    want = hashlib.sha256(b"pk1:a:uint16:2,3:r:" + raw).hexdigest()
    assert vhash(to_mvalue(a)) == want
    h1, h2 = vhash(mchar("x")), vhash(to_mvalue(2.0))
    assert (
        vhash(MCell((1, 2), [mchar("x"), to_mvalue(2.0)]))
        == hashlib.sha256(f"pk1:cell:1,2:{h1};{h2}".encode()).hexdigest()
    )
    s = MStruct((1, 1), ["b", "a"], [{"b": mchar("x"), "a": to_mvalue(2.0)}])
    assert vhash(s) == hashlib.sha256(f"pk1:struct:1,1:a,b:{h2},{h1}".encode()).hexdigest()


def test_hash_ignores_struct_field_order_but_not_class_or_size():
    a = MStruct((1, 1), ["a", "b"], [{"a": to_mvalue(1), "b": to_mvalue(2)}])
    b = MStruct((1, 1), ["b", "a"], [{"a": to_mvalue(1), "b": to_mvalue(2)}])
    assert vhash(a) == vhash(b)
    assert vhash(to_mvalue(np.zeros((2, 3)))) != vhash(to_mvalue(np.zeros((3, 2))))
    assert vhash(to_mvalue(np.zeros(3))) != vhash(to_mvalue(np.zeros(3, np.float32)))


def test_descriptors_keep_small_arrays_and_hashes_of_large_ones():
    small, big = to_mvalue(vol((4, 5, 3))), to_mvalue(vol((40, 50, 30)))
    ds, db = describe(small), describe(big)
    assert "b64" in ds and "b64" not in db
    assert compare(from_descriptor(json.loads(json.dumps(ds))), small.data).ok
    assert compare(from_descriptor(db), big.data).ok
    changed = big.data.copy()
    changed[3, 4, 5] += 1
    rep = compare(from_descriptor(db), changed, path="out1")
    assert not rep.ok and "SHA-256 mismatch" in rep.text()
    assert "readable element-by-element diff is in the golden store" in rep.text()


def roundtrip(d):
    return from_descriptor(json.loads(json.dumps(d)))


def test_large_containers_keep_only_their_hash_and_compare_bit_for_bit():
    rows = [f"tile_{i:04d}_CamA_ch0_stack{i:04d}.tif" for i in range(1600)]
    c = to_mvalue(rows)
    d = describe(c)
    assert "v" not in d and len(json.dumps(d)) < 200
    assert d["h"] == describe(c, budget=10**9)["h"] == vhash(c)
    ref = roundtrip(d)
    assert compare(ref, rows).ok
    bad = list(rows)
    bad[7] = "x"
    text = compare(ref, bad, path="out1").text()
    assert "out1: values: the cell array [1 1600] differs from MATLAB's" in text
    assert "readable element-by-element diff is in the golden store" in text
    assert "cell size: MATLAB [1 1600], Python [1 1599]" in compare(ref, rows[:-1]).text()
    # few elements, but long: the size counts, not the number of elements
    long_rows = [f"{i:03d}" + "x" * 300 for i in range(100)]
    assert "v" not in describe(to_mvalue(long_rows))
    assert "v" in describe(to_mvalue(long_rows[:20]))
    # an old manifest kept every element: still read, with the element-by-element report
    full = roundtrip(describe(c, budget=10**9))
    assert "value{8}: text" in compare(full, bad).text()
    # a small field next to a large one keeps its data
    s = MStruct((1, 1), ["n", "names"], [{"n": to_mvalue(3.0), "names": c}])
    ds = describe(s)
    assert "b64" in ds["v"][0]["n"] and "v" not in ds["v"][0]["names"]
    assert len(json.dumps(ds, separators=(",", ":"))) <= C.DESCRIPTOR_BYTES
    assert "value.n: values" in compare(roundtrip(ds), {"n": 4.0, "names": rows}).text()


def test_char_rows_are_kept_as_text_and_element_hashes_computed_again():
    texts = ["abc", "héllo \U0001f600", ""]
    c = to_mvalue(texts)
    d = json.loads(json.dumps(describe(c)))
    assert [x.get("t") for x in d["v"]] == texts and all("h" not in x for x in d["v"])
    assert vhash(roundtrip(d)) == vhash(c) == d["h"]
    lone = MArray("char", np.array([[0xD800, 0x41]], "<u2"))  # a lone surrogate (UTF-16)
    assert vhash(roundtrip(describe(lone))) == vhash(lone)
    empty = MArray("char", np.zeros((1, 0), "<u2"))
    assert roundtrip(describe(empty)).size == (1, 0)


def test_large_struct_arrays_and_maps_keep_fields_and_types():
    els = [{"a": to_mvalue(float(i)), "b": mchar(f"x{i}")} for i in range(500)]
    sa = MStruct((1, 500), ["a", "b"], els)
    d = describe(sa)
    assert "v" not in d and d["f"] == ["a", "b"]
    py = [{"a": float(i), "b": f"x{i}"} for i in range(500)]
    assert compare(roundtrip(d), py).ok
    text = compare(roundtrip(d), [{"a": float(i), "c": f"x{i}"} for i in range(500)]).text()
    assert "only in MATLAB ['b'], only in Python ['c']" in text
    keys = [f"k{i:04d}" for i in range(2000)]
    mp = MMap(
        "char", "any", to_mvalue(keys), MCell((1, 2000), [to_mvalue(float(i)) for i in range(2000)])
    )
    dm = describe(mp)
    assert "keys" not in dm and dm["kt"] == "char"
    assert compare(roundtrip(dm), {k: float(i) for i, k in enumerate(keys)}).ok
    assert not compare(roundtrip(dm), dict.fromkeys(keys, 0.0)).ok


# a Python result is read as the MATLAB reference (to_mvalue's rules); with a hash-only manifest
# it must be read the same way, or a CI pass would mean something else than a store pass
READING = [
    ("empty list as empty cell", MCell((0, 0), []), []),
    (
        "struct field: empty list as empty cell",
        MStruct((1, 1), ["name", "other"], [{"name": mchar("abc"), "other": MCell((0, 0), [])}]),
        {"name": "abc", "other": []},
    ),
    ("int as uint16", to_mvalue(np.array([[7]], np.uint16)), 7),
    ("1-D as column", to_mvalue(np.arange(3.0).reshape(3, 1)), np.arange(3.0)),
    ("empty str as 1x0 char", MArray("char", np.zeros((1, 0), "<u2")), ""),
    ("None as empty struct", MStruct((0, 0), ["a"], []), None),
    (
        "list of dicts as struct array",
        MStruct((1, 2), ["a"], [{"a": to_mvalue(1.0)}, {"a": to_mvalue(2.0)}]),
        [{"a": 1.0}, {"a": 2.0}],
    ),
    (
        "dict as containers.Map",
        MMap("char", "double", to_mvalue(["a"]), MCell((1, 1), [to_mvalue(1.0)])),
        {"a": 1.0},
    ),
    ("str as string", MString((1, 1), ["x"]), "x"),
    (
        "dict as inputParser value",
        MObject("inputParser", to_mvalue({"alpha": 5.0}), ""),
        {"alpha": 5},
    ),
]


@pytest.mark.parametrize(("ref", "py"), [r[1:] for r in READING], ids=[r[0] for r in READING])
def test_results_are_read_alike_with_the_store_and_with_a_hash_only_manifest(ref, py):
    big = MCell((1, 400), [ref] * 400)  # large enough that the manifest keeps only its hash
    d = json.loads(json.dumps(describe(big)))
    assert "v" not in d
    store, ci = compare(big, [py] * 400), compare(from_descriptor(d), [py] * 400)
    assert store.ok, store.text()
    assert ci.ok, ci.text()
    assert not store.undecided and not ci.undecided


def test_hash_only_manifests_still_find_differences_and_leave_tolerances_undecided():
    big = MCell((1, 400), [MArray("char", np.zeros((1, 0), "<u2"))] * 400)
    py = [""] * 399 + ["x"]
    assert not compare(big, py).ok
    ci = compare(roundtrip(describe(big)), py)
    assert not ci.ok and "differs from MATLAB's" in ci.text()
    r = to_mvalue(np.linspace(1, 2, 1000))  # 8 KB: kept by hash only
    a = np.nextafter(r.data, 3.0)
    spec = tol(ulp=1, reason="parity/example.md: test")
    assert compare(r, a, spec).ok
    ci = compare(roundtrip(describe(r)), a, spec)
    assert ci.ok and ci.undecided and "not decided without the golden store" in ci.text()
    assert not compare(roundtrip(describe(r)), a).ok  # exact: a difference


def test_a_manifest_without_reading_signatures_decides_only_equal_hashes():
    # written before the signature: a matching hash still passes, any other is undecided
    big = MCell((1, 400), [MCell((0, 0), [])] * 400)
    old = describe(big)
    del old["r"]
    assert compare(roundtrip(old), [MCell((0, 0), [])] * 400).ok
    ci = compare(roundtrip(old), [[]] * 400)
    assert ci.ok and "lacks the reading signature" in ci.undecided[0]


# --- seeded mismatches give readable reports ----------------------------------------------------


def test_one_voxel():
    r = vol()
    a = r.copy()
    a[7, 3, 4] += 1
    rep, text = report(r, a)
    assert not rep.ok
    assert "1 of 960 values differ" in text
    assert "index 517 (0-based, column-major) = MATLAB (8,4,5)" in text
    assert "largest |difference| 1" in text and "histogram" in text


def test_wrong_dtype():
    r = vol()
    rep, text = report(r, r.astype(np.float64))
    assert "class: MATLAB uint16, Python double" in text
    assert "values are equal after conversion to double" in text


def test_shifted_axis():
    r = vol()
    a = np.roll(r, 1, axis=1)
    rep, text = report(r, a)
    assert "shifted by +1 along dimension 2 (x)" in text


def test_permuted_axes():
    r = vol()
    rep, text = report(r, np.transpose(r, (2, 0, 1)))
    assert "size: MATLAB [10 12 8], Python [8 10 12]" in text
    assert "permutation" in text and "permute(MATLAB, [3 1 2])" in text


def test_c_order_values():
    r = vol((6, 6, 6))
    a = r.reshape(-1, order="F").reshape(r.shape, order="C")
    rep, text = report(r, a)
    assert "C (row-major) order" in text


def test_size_off_by_one():
    r = vol()
    rep, text = report(r, r[:, :, :-1])
    assert "dimension 3 (z) differs by -1" in text
    assert "Python lacks the last 1 element(s) along dimension 3" in text


def test_nan_bits_and_signed_zero():
    r = np.array([MATLAB_NAN, 0.0, 1.0])
    rep, text = report(r, np.array([np.nan, -0.0, 1.0]))
    assert "differ only in NaN bits" in text and "sign of zero" in text


def test_nested_paths_name_the_element():
    r = to_mvalue(
        [{"a": 1}, {"a": 2}],
        like=MStruct((1, 2), ["a"], [{"a": to_mvalue(1)}, {"a": to_mvalue(1)}]),
    )
    a = [{"a": 1}, {"a": 3}]
    rep = compare(r, a, path="out2")
    assert "out2(2).a: values" in rep.text()
    c = MCell((1, 2), [mchar("x"), to_mvalue(np.arange(3.0))])
    rep = compare(c, ["x", np.array([0.0, 1, 5])], path="out1")
    assert "out1{2}: values" in rep.text()
    rep = compare(
        MStruct((1, 1), ["a", "b"], [{"a": to_mvalue(1), "b": to_mvalue(2)}]), {"a": 1, "c": 2}
    )
    assert "only in MATLAB ['b'], only in Python ['c']" in rep.text()


def test_kind_and_text_mismatches():
    rep = compare(mchar("abc"), "abd")
    assert "first difference at character 3" in rep.text()
    rep = compare(to_mvalue(["a"]), "a")
    assert "MATLAB: cell array [1 1], Python: char array [1 1]" in rep.text()


# --- tolerances ---------------------------------------------------------------------------------


def test_tolerances_need_a_reason():
    with pytest.raises(ValueError):
        tol(abs=1, reason="")


def test_tolerance_counts_and_magnitudes():
    r = vol((20, 20, 20))
    a = r.copy()
    a[1, 2, 3] += 1
    spec = tol(abs=1, max_frac=2e-4, reason="parity/example.md: test")
    rep = compare(to_mvalue(r), a, spec)
    assert rep.ok and "within tolerance" in rep.notes[0]
    a[4, 5, 6] += 2
    rep = compare(to_mvalue(r), a, spec)
    assert not rep.ok and "outside tolerance" in rep.text()
    a = r.copy()
    a[1, 2, 3] += 1
    a[2, 2, 3] += 1  # 2 of 8000 values: 2.5e-4, more than max_frac
    rep = compare(to_mvalue(r), a, spec)
    assert not rep.ok and "of values differ, more than 0.0002" in rep.text()


def test_tolerance_in_ulp_and_relative_l2():
    r = np.linspace(1, 2, 50).astype(np.float32)
    a = np.nextafter(r, np.float32(3))
    assert compare(to_mvalue(r), a, tol(ulp=1, reason="x")).ok
    assert not compare(to_mvalue(r), a, tol(ulp=0, reason="x")).ok
    assert compare(to_mvalue(r), a, tol(rel_l2=1e-6, reason="x")).ok
    assert not compare(to_mvalue(r), a * 1.01, tol(rel_l2=1e-6, reason="x")).ok
    assert not compare(
        to_mvalue(np.array([1.0, MATLAB_NAN])), np.array([1.0, 2.0]), tol(abs=10, reason="x")
    ).ok


# --- masks --------------------------------------------------------------------------------------


def test_masks_hide_the_work_folder_host_uuids_and_times():
    m = Masks().with_paths(work="/tmp/w1", host="node12")
    s = m.text(
        "/tmp/w1/a_0a1b2c3d-1111-2222-3333-444455556666.tif on node12 at 2026-11-11 10:00:01"
    )
    assert s == "<WORK>/a_<uuid>.tif on <host> at <time>"
    m2 = Masks().with_paths(work="C:\\Temp\\w2")
    assert m2.text("C:/Temp/w2/x and C:\\Temp\\w2\\y") == "<WORK>/x and <WORK>\\y"
    assert compare(
        mchar("/tmp/w1/x"),
        "/tmp/w2/x",
        masks=Masks().with_paths(work="/tmp/w1"),
        act_masks=Masks().with_paths(work="/tmp/w2"),
    ).ok


# --- file trees ---------------------------------------------------------------------------------


def write_tree(root, im, *, tiff=True, zarr=True, text="line 1\n", js=None, mat=True, extra=None):
    root.mkdir(parents=True, exist_ok=True)
    if tiff:
        cpptiff = pytest.importorskip("cpptiff")
        (root / "sub").mkdir(exist_ok=True)
        cpptiff.write_tiff(
            str(root / "sub" / "im.tif"), np.ascontiguousarray(np.transpose(im, (2, 0, 1)))
        )
    if zarr:
        cppzarr = pytest.importorskip("cppzarr")
        cppzarr.write_zarr(
            str(root / "im.zarr"),
            np.asfortranarray(im),
            chunks=[5, 6, 4],
            order="F",
            dimension_separator=".",
        )
    (root / "log.txt").write_text(text)
    (root / "params.json").write_text(js if js is not None else '{"a":1,"b":"x"}')
    if mat:
        from scipy.io import savemat

        savemat(root / "info.mat", {"name": "demo", "sz": np.array([[10.0, 12, 8]]), "t": 1.0})
    (root / "empty").mkdir(exist_ok=True)
    (root / "flag_0a1b2c3d-1111-2222-3333-444455556666.txt").write_text("")
    if extra:
        (root / extra).write_bytes(b"x")


def test_equal_trees_and_their_descriptors(tmp_path):
    im = vol()
    write_tree(tmp_path / "a", im)
    write_tree(tmp_path / "b", im)
    (tmp_path / "b" / "flag_0a1b2c3d-1111-2222-3333-444455556666.txt").rename(
        tmp_path / "b" / "flag_99999999-1111-2222-3333-444455556666.txt"
    )
    ra, rb = scan_tree(tmp_path / "a"), scan_tree(tmp_path / "b")
    rep = compare_trees(ra, rb)
    assert rep.ok, rep.text()
    d = json.loads(json.dumps(describe_tree(ra, tree(), Masks())))
    assert compare_trees(tree_from_descriptor(d), rb).ok


def test_missing_extra_and_changed_files(tmp_path):
    im = vol()
    write_tree(tmp_path / "a", im)
    b = im.copy()
    b[0, 0, 0] += 3
    write_tree(tmp_path / "b", b, text="line 2\n", js='{"a":2,"b":"x"}', extra="new.bin")
    (tmp_path / "b" / "info.mat").unlink()
    rep = compare_trees(scan_tree(tmp_path / "a"), scan_tree(tmp_path / "b"))
    text = rep.text("case")
    assert "files/info.mat: missing: MATLAB wrote this mat, Python did not" in text
    assert "files/new.bin: extra: Python wrote this, MATLAB did not" in text
    assert "files/sub/im.tif (voxels): values: 1 of 960 values differ" in text
    assert "files/im.zarr (voxels): values" in text
    assert (
        "files/log.txt: text: first difference on line 1: MATLAB 'line 1', Python 'line 2'" in text
    )
    assert "files/params.json: values: JSON .a: MATLAB 1, Python 2" in text


def test_json_formatting_and_mat_contents(tmp_path):
    from scipy.io import savemat

    for side, js, v in (("a", '{"a":1}', 2.0), ("b", '{"a": 1}', 3.0)):
        d = tmp_path / side
        d.mkdir()
        (d / "p.json").write_text(js)
        savemat(d / "m.mat", {"x": np.array([[1.0, v]]), "t": v})
    rep = compare_trees(scan_tree(tmp_path / "a"), scan_tree(tmp_path / "b"))
    text = rep.text()
    assert "files/p.json: format: JSON content is equal but the text differs" in text
    assert "files/m.mat.x: values: 1 of 2 values differ" in text  # t is ignored by default
    rep = compare_trees(
        scan_tree(tmp_path / "a"),
        scan_tree(tmp_path / "b"),
        tree(text_format=False, ignore_vars=("t", "x")),
    )
    assert rep.ok, rep.text()


def test_read_mat_v5_classes(tmp_path):
    from scipy.io import savemat

    savemat(
        tmp_path / "v.mat",
        {
            "u": np.array([[1, 2]], np.uint16),
            "c": np.array(["ab", "cd"]),
            "l": np.array([[True, False]]),
            "cell": np.array([1.0, "x"], dtype=object),
            "s": {"f": 1.0},
        },
    )
    m = read_mat(tmp_path / "v.mat")
    v = m.items[0]
    assert (v["u"].mclass, v["u"].size) == ("uint16", (1, 2))
    assert v["c"].mclass == "char" and v["c"].size == (2, 2)
    assert v["l"].mclass in ("logical", "uint8")
    assert isinstance(v["cell"], MCell) and v["cell"].size == (1, 2)
    assert isinstance(v["s"], MStruct) and v["s"].fields == ["f"]


def test_zarr_chunk_sets_are_a_note_by_default():
    # MATLAB writes all-zero chunks, cpp-zarr's Python writes skip them: same data, other files
    data = describe(to_mvalue(vol((4, 4, 2))))
    z = {"k": "zarr", "meta": {".zarray": "{}"}, "data": data}
    ref = tree_from_descriptor({"files": {"x.zarr": [dict(z, nchunks=4, chunks_sha="a")]}})
    act = tree_from_descriptor({"files": {"x.zarr": [dict(z, nchunks=3, chunks_sha="b")]}})
    rep = compare_trees(ref, act)
    assert rep.ok and "chunk files differ: MATLAB 4, Python 3" in rep.notes[0]
    rep = compare_trees(ref, act, tree(zarr_chunks="strict"))
    assert not rep.ok and "layout" in rep.text()


def test_compare_finds_different_outputs(tmp_path):
    # the intent of IsolationHarnessTest/compareFindsDifferentOutputs (MATLAB)
    cpptiff = pytest.importorskip("cpptiff")
    cppzarr = pytest.importorskip("cppzarr")
    from scipy.io import savemat

    im = np.asfortranarray(np.arange(1, 4001, dtype=np.uint16).reshape((10, 20, 20), order="F"))
    im2 = im.copy(order="F")
    im2[4, 0, 0] = 0
    a, b = tmp_path / "a", tmp_path / "b"

    def tif(path, x):
        cpptiff.write_tiff(str(path), np.ascontiguousarray(np.transpose(x, (2, 0, 1))))

    def zarr(path, x, chunks):
        cppzarr.write_zarr(str(path), x, chunks=chunks, order="F", dimension_separator=".")

    for d in (a, b):
        (d / "sub").mkdir(parents=True)
        tif(d / "same.tif", im)
        zarr(d / "same.zarr", im, [8, 8, 8])
        savemat(d / "sub" / "same.mat", {"value": 1.0})
    zarr(a / "chunks.zarr", im, [8, 8, 8])
    zarr(b / "chunks.zarr", im, [16, 16, 16])
    tif(a / "diff.tif", im)
    tif(b / "diff.tif", im2)
    zarr(a / "diff.zarr", im, [8, 8, 8])
    zarr(b / "diff.zarr", im2, [8, 8, 8])
    savemat(a / "diff.mat", {"value": 1.0})
    savemat(b / "diff.mat", {"value": 2.0})
    (a / "tmp_0d8c8a76-1b6e-4c3e-9a51-0f5f2b6c1e11.txt").write_text("")
    (b / "tmp_ffa1b2c3-1111-2222-3333-444455556666.txt").write_text("")
    (a / "onlyA.txt").write_text("")

    rep = compare_trees(scan_tree(a), scan_tree(b))
    text = rep.text()
    for name in ("same.tif", "same.zarr", "sub/same.mat", "tmp_<uuid>.txt"):
        assert f"files/{name}" not in text, text
    assert "files/chunks.zarr/.zarray: meta" in text  # zarr metadata (chunks) is compared
    assert "files/diff.tif (voxels): values: 1 of 4000 values differ" in text
    assert "files/diff.zarr (voxels): values: 1 of 4000 values differ" in text
    assert "files/diff.mat.value: values" in text
    assert "files/onlyA.txt: missing" in text


def test_an_empty_flag_file_named_mat_is_compared_by_its_bytes(tmp_path):
    """a 0-byte .mat (a flag file) is not a MAT-file: two of them are equal, in the store and
    in a manifest; an empty one against a real MAT-file, or other bytes, is a difference"""
    from scipy.io import savemat

    for side in ("a", "b", "c", "d"):
        (tmp_path / side).mkdir()
    (tmp_path / "a" / "done.mat").write_bytes(b"")
    (tmp_path / "b" / "done.mat").write_bytes(b"")
    savemat(tmp_path / "c" / "done.mat", {"x": np.array([[1.0]])})
    (tmp_path / "d" / "done.mat").write_bytes(b"not a mat")
    ref = scan_tree(tmp_path / "a")
    rep = compare_trees(ref, scan_tree(tmp_path / "b"))
    assert rep.ok, rep.text()
    assert any("not a readable MAT-file" in n for n in rep.notes)
    spec = tree()
    d = json.loads(json.dumps(describe_tree(ref, spec, C.Masks())))
    entry = d["files"]["done.mat"][0]
    assert "vars" not in entry and entry["sha"] == hashlib.sha256(b"").hexdigest()
    rep = compare_trees(tree_from_descriptor(d), scan_tree(tmp_path / "b"), spec)
    assert rep.ok, rep.text()
    text = compare_trees(ref, scan_tree(tmp_path / "c")).text()
    assert "MATLAB's file is not a readable MAT-file" in text, text
    text = compare_trees(scan_tree(tmp_path / "c"), ref).text()
    assert "Python's file is not a readable MAT-file" in text, text
    text = compare_trees(ref, scan_tree(tmp_path / "d")).text()
    assert "not a readable MAT-file on either side" in text and "bytes differ" in text, text
