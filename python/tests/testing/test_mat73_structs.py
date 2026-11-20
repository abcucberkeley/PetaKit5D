"""The harness reads struct arrays of v7.3 MAT-files as MATLAB R2024b's load gives them.

MATLAB writes a struct array without elements (``struct('a', {})``, 0x3, 0x2x3) and one without
fields (``struct()``, ``repmat(struct(), 1, 3)``) as a dataset of its size with MATLAB_empty set,
and the field names of an empty one in MATLAB_fields. The harness keeps those names, the size
and the elements without fields, as MATLAB's own decoding does: the oracle case saves such
values in MATLAB and in the port (``io.matfile.write_mat``), and its stored MATLAB file must read
as MATLAB decoded it.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.io.matfile import StructArray, write_mat
from petakit5d.testing import compare as C
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case
from petakit5d.testing.compare import MCell, MStruct, compare, read_mat, read_pkv, to_mvalue

F34 = tuple(f"f{i:02d}" for i in range(1, 35))  # as many fields as tracks of ProcessedTracks.mat
LONG = "a" * 63  # the longest MATLAB name
# the variables: name -> (MATLAB expression, the port's value)
VARIABLES = {
    "s11": ("struct()", dict),
    "s00": ("struct([])", lambda: StructArray.empty()),
    "e00f": ("struct('b', {}, 'a', {})", lambda: StructArray.empty(("b", "a"))),
    "e03f": ("reshape(struct('a', {}), 0, 3)", lambda: StructArray.empty(("a",), (0, 3))),
    "e30f": ("reshape(struct('a', {}), 3, 0)", lambda: StructArray.empty(("a",), (3, 0))),
    "e10f": ("struct('a', cell(1, 0))", lambda: StructArray.empty(("a",), (1, 0))),
    "e023f": ("reshape(struct('a', {}), 0, 2, 3)", lambda: StructArray.empty(("a",), (0, 2, 3))),
    "e03n": ("reshape(struct([]), 0, 3)", lambda: StructArray.empty((), (0, 3))),
    "tracks": (
        "reshape(cell2struct(cell(34, 0), {" + ", ".join(f"'{f}'" for f in F34) + "}, 1), 0, 0)",
        lambda: StructArray.empty(F34),
    ),
    "elong": (f"struct('{LONG}', {{}})", lambda: StructArray.empty((LONG,))),
    "n13": ("repmat(struct(), 1, 3)", lambda: StructArray((1, 3), (), [{}] * 3)),
    "n32": ("repmat(struct(), 3, 2)", lambda: StructArray((3, 2), (), [{}] * 6)),
    "a13": (
        "struct('a', {1, 2, 3})",
        lambda: StructArray((1, 3), ("a",), [{"a": float(i)} for i in (1, 2, 3)]),
    ),
    "a32": (
        "repmat(struct('b', 1, 'a', 'x'), 3, 2)",
        lambda: StructArray((3, 2), ("b", "a"), [{"b": 1.0, "a": "x"}] * 6),
    ),
    "nest1": (
        "struct('e', struct('b', {}, 'a', {}), 'u', 1)",
        lambda: {"e": StructArray.empty(("b", "a")), "u": 1.0},
    ),
    "nest2": (
        "struct('x', {1, 2}, 'e', {struct('q', {}), reshape(struct('r', {}, 's', {}), 0, 2)})",
        lambda: StructArray(
            (1, 2),
            ("x", "e"),
            [
                {"x": 1.0, "e": StructArray.empty(("q",))},
                {"x": 2.0, "e": StructArray.empty(("r", "s"), (0, 2))},
            ],
        ),
    ),
    "nestc": (
        "{struct('a', {}), 1, reshape(struct([]), 2, 0)}",
        lambda: [StructArray.empty(("a",)), 1.0, StructArray.empty((), (2, 0))],
    ),
    "nestn": (
        "struct('e', repmat(struct(), 1, 2))",
        lambda: {"e": StructArray((1, 2), (), [{}] * 2)},
    ),
}
FILE = "structs73.mat"
SAVE = (
    "".join(f"{n} = {e}; " for n, (e, _) in VARIABLES.items())
    + f"save('{FILE}', "
    + ", ".join(f"'{n}'" for n in VARIABLES)
    + ", '-v7.3'); out = {};"
)


def write_structs() -> None:
    write_mat(FILE, {n: make() for n, (_, make) in VARIABLES.items()}, version="7.3")


CASES = [Case("mat73.struct_arrays", code=SAVE, python=write_structs, nargout=0, files=True)]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def field_orders(v, path="") -> list:
    """(path, fields) of every struct in a value, in order (comparisons only note the order)"""
    out = []
    if isinstance(v, MStruct):
        out.append((path, list(v.fields)))
        for i, el in enumerate(v.items or []):
            for k in v.fields:
                out += field_orders(el[k], f"{path}({i + 1}).{k}")
    elif isinstance(v, MCell):
        for i, x in enumerate(v.items or []):
            out += field_orders(x, f"{path}{{{i + 1}}}")
    return out


@pytest.mark.oracle
def test_matlab_files_read_as_matlab_decoded_them():
    store = oracle.Store.from_env()
    case = CASES[0]
    if store is None or not store.has(case.key()):
        pytest.skip("needs the golden store with this module's reference")
    it = store.load(case.key()).files.items[FILE][0]
    theirs, mine = read_pkv(it.decoded), read_mat(it.path)
    rep = compare(theirs, mine, path=FILE)
    assert rep.ok and not rep.notes, rep.text()
    assert field_orders(mine) == field_orders(theirs)


# --- files as MATLAB lays them out, without MATLAB ------------------------------------------------


def read_one(tmp_path, value, version="7.3"):
    path = tmp_path / "v.mat"
    write_mat(path, {"v": value}, version=version)
    return read_mat(path).items[0]["v"]


@pytest.mark.parametrize(
    "size, fields",
    [((0, 0), ("b", "a")), ((0, 3), ("a",)), ((3, 0), ("a",)), ((0, 2, 3), ("a",)), ((0, 0), F34)],
)
def test_struct_arrays_without_elements_keep_their_fields_and_size(tmp_path, size, fields):
    v = read_one(tmp_path, StructArray.empty(fields, size))
    assert isinstance(v, MStruct)
    assert (v.size, v.fields, v.items) == (size, list(fields), [])
    assert C.vhash(v) == C.vhash(MStruct(size, list(fields), []))


@pytest.mark.parametrize("size", [(1, 1), (1, 3), (3, 2), (0, 3), (0, 0)])
def test_struct_arrays_without_fields_keep_their_elements(tmp_path, size):
    n = int(np.prod(size))
    value = {} if size == (1, 1) else StructArray(size, (), [{}] * n)
    v = read_one(tmp_path, value)
    assert (v.size, v.fields, v.items) == (size, [], [{}] * n)


def test_nested_struct_arrays_keep_their_fields_and_elements(tmp_path):
    path = tmp_path / "n.mat"
    nested = ("nest1", "nest2", "nestc", "nestn")
    write_mat(path, {n: VARIABLES[n][1]() for n in nested}, version="7.3")
    got = read_mat(path).items[0]
    assert got["nest1"].items[0]["e"].fields == ["b", "a"]
    e2 = [el["e"] for el in got["nest2"].items]
    assert [(e.size, e.fields) for e in e2] == [((0, 0), ["q"]), ((0, 2), ["r", "s"])]
    c = got["nestc"].items
    assert (c[0].size, c[0].fields) == ((0, 0), ["a"]) and (c[2].size, c[2].fields) == ((2, 0), [])
    assert got["nestn"].items[0]["e"].items == [{}, {}]
    # each equals the value the harness makes of MATLAB's (to_mvalue of the same structs)
    empty = MStruct((0, 0), ["b", "a"], [])
    want = MStruct((1, 1), ["e", "u"], [{"e": empty, "u": to_mvalue(1.0)}])
    assert C.vhash(got["nest1"]) == C.vhash(want)


def test_field_order_is_matlabs(tmp_path):
    # MATLAB_fields gives the order of a struct's fields (HDF5 lists the members by name)
    v = read_one(tmp_path, StructArray((1, 2), ("x", "e"), [{"x": 1.0, "e": 2.0}] * 2))
    assert v.fields == ["x", "e"]
    assert read_one(tmp_path, {"x": 1.0, "e": 2.0}).fields == ["x", "e"]


@pytest.mark.parametrize("value", [{"a": 1.0}, StructArray((1, 3), ("a",), [{"a": 1.0}] * 3)])
def test_a_struct_of_one_field_needs_no_field_names(tmp_path, value):
    # MATLAB writes no MATLAB_fields for a struct (array) of one field
    import h5py

    path = tmp_path / "one.mat"
    write_mat(path, {"v": value}, version="7.3")
    with h5py.File(path, "r+") as f:
        del f["v"].attrs["MATLAB_fields"]
    v = read_mat(path).items[0]["v"]
    size = value.shape if isinstance(value, StructArray) else (1, 1)
    assert (v.size, v.fields) == (size, ["a"])
    assert [el["a"].data.item() for el in v.items] == [1.0] * len(v.items)


def matlab_layout(path, size, fields=()) -> None:
    """a struct written as MATLAB writes one without elements or fields (h5py, MATLAB's
    attributes: the class as a fixed-length string, MATLAB_empty uint8, the size as uint64, the
    names as variable-length lists of one-character strings)"""
    import h5py

    with h5py.File(path, "w", userblock_size=512) as f:
        ds = f.create_dataset("v", data=np.array(size, np.uint64))
        tid = h5py.h5t.C_S1.copy()
        tid.set_size(6)
        tid.set_strpad(h5py.h5t.STR_NULLTERM)
        aid = h5py.h5a.create(ds.id, b"MATLAB_class", tid, h5py.h5s.create(h5py.h5s.SCALAR))
        aid.write(np.array(b"struct", dtype="S6"), mtype=tid)
        ds.attrs.create("MATLAB_empty", np.uint8(1))
        if fields:
            vt = h5py.vlen_dtype(np.dtype("S1"))
            names = np.empty(len(fields), dtype=vt)
            for i, k in enumerate(fields):
                names[i] = np.frombuffer(k.encode("ascii"), dtype="S1")
            ds.attrs.create("MATLAB_fields", names, dtype=vt)
    with open(path, "r+b") as fh:
        fh.write(b"MATLAB 7.3 MAT-file".ljust(116) + bytes(8) + b"\x00\x02IM")


def test_matlabs_layout_of_an_empty_struct_array_is_read(tmp_path):
    # tracks of ProcessedTracks.mat: 0x0 with 34 fields
    matlab_layout(tmp_path / "t.mat", (0, 0), F34)
    v = read_mat(tmp_path / "t.mat").items[0]["v"]
    assert (v.size, v.fields, v.items) == ((0, 0), list(F34), [])
    matlab_layout(tmp_path / "n.mat", (1, 3))
    v = read_mat(tmp_path / "n.mat").items[0]["v"]
    assert (v.size, v.fields, v.items) == ((1, 3), [], [{}, {}, {}])


def test_a_struct_with_fields_but_no_values_is_not_a_matlab_file(tmp_path):
    # elements and field names, but no values: not what MATLAB writes; the file is unreadable,
    # and a comparison takes it byte for byte
    matlab_layout(tmp_path / "x.mat", (1, 2), ("a",))
    with pytest.raises(ValueError, match="without their values"):
        read_mat(tmp_path / "x.mat")
