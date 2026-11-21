"""The harness reads struct arrays of v6 and v7 MAT-files as MATLAB R2024b's load gives them, and
reads the port's ``io.matfile.StructArray`` results as MATLAB's struct arrays.

SciPy's loadmat reads a struct array without fields (``struct()``, ``repmat(struct(), 1, 3)``,
``struct([])``) as an object array of None, an empty one as an empty cell, and an element that
MATLAB writes without data (those of ``cell(1, 3)``, never assigned) as a 1x0 double. Where a
value holds such arrays, the harness takes the classes from the file (its MAT 5 elements), so
those arrays are structs of their size and those elements ``[]``, in cells and structs too, and
empty cells stay cells (MATLAB's v7.3 file stores those elements as class "canonical empty",
also ``[]``). ``to_mvalue`` reads a StructArray as the struct array of its size and fields. The
oracle cases save such values in MATLAB and in the port (``io.matfile.write_mat``: v7, v6 and
v7.3) and return them from MATLAB and from Python; MATLAB's own files must read as MATLAB decoded
them, and the port's values must convert to what the readers give for MATLAB's files.
"""

from __future__ import annotations

import struct

import numpy as np
import pytest

from petakit5d.io.matfile import StructArray, write_mat
from petakit5d.testing import compare as C
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case
from petakit5d.testing.compare import (
    MArray,
    MCell,
    MHashed,
    MStruct,
    compare,
    describe,
    from_descriptor,
    read_mat,
    read_pkv,
    to_mvalue,
    vhash,
)

# the values: name -> (MATLAB expression, the port's value)
VARIABLES = {
    # struct arrays without fields (SciPy: object arrays of None, an empty one an empty cell)
    "s11": ("struct()", dict),
    "s00": ("struct([])", lambda: StructArray.empty()),
    "n13": ("repmat(struct(), 1, 3)", lambda: StructArray((1, 3), (), [{}] * 3)),
    "n31": ("repmat(struct(), 3, 1)", lambda: StructArray((3, 1), (), [{}] * 3)),
    "n32": ("repmat(struct(), 3, 2)", lambda: StructArray((3, 2), (), [{}] * 6)),
    "n222": ("repmat(struct(), [2 2 2])", lambda: StructArray((2, 2, 2), (), [{}] * 8)),
    "n300": ("repmat(struct(), 1, 300)", lambda: StructArray((1, 300), (), [{}] * 300)),
    "e03n": ("reshape(struct([]), 0, 3)", lambda: StructArray.empty((), (0, 3))),
    "e30n": ("reshape(struct([]), 3, 0)", lambda: StructArray.empty((), (3, 0))),
    "e10n": ("reshape(struct([]), 1, 0)", lambda: StructArray.empty((), (1, 0))),
    "e203n": ("reshape(struct([]), 2, 0, 3)", lambda: StructArray.empty((), (2, 0, 3))),
    # struct arrays with fields
    "a11": ("struct('b', 1, 'a', 'x')", lambda: {"b": 1.0, "a": "x"}),
    "a13": (
        "struct('a', {1, 2, 3})",
        lambda: StructArray((1, 3), ("a",), [{"a": float(i)} for i in (1, 2, 3)]),
    ),
    "a32": (
        "repmat(struct('b', 1, 'a', 'x'), 3, 2)",
        lambda: StructArray((3, 2), ("b", "a"), [{"b": 1.0, "a": "x"}] * 6),
    ),
    "a300": (
        "struct('a', num2cell(1 : 300))",
        lambda: StructArray((1, 300), ("a",), [{"a": float(i)} for i in range(1, 301)]),
    ),
    "e00f": ("struct('b', {}, 'a', {})", lambda: StructArray.empty(("b", "a"))),
    "e023f": ("reshape(struct('a', {}), 0, 2, 3)", lambda: StructArray.empty(("a",), (0, 2, 3))),
    # cells (empty ones too), which stay cells
    "c00": ("{}", list),
    "c03": ("cell(0, 3)", lambda: np.empty((0, 3), dtype=object)),
    "c13": ("cell(1, 3)", lambda: [None] * 3),
    "cc": ("{{}}", lambda: [[]]),
    # nested in cells and structs
    "mixc": (
        "{struct(), {}, struct([]), cell(0, 2), repmat(struct(), 2, 1), {struct()}}",
        lambda: [
            {},
            [],
            StructArray.empty(),
            np.empty((0, 2), dtype=object),
            StructArray((2, 1), (), [{}] * 2),
            [{}],
        ],
    ),
    "nests": (
        "struct('e', struct(), 'c', {{}}, 'n', repmat(struct(), 1, 2), 'z', struct([]))",
        lambda: {
            "e": {},
            "c": [],
            "n": StructArray((1, 2), (), [{}] * 2),
            "z": StructArray.empty(),
        },
    ),
    "nesta": (
        "struct('x', {1, 2}, 'e', {struct(), reshape(struct([]), 0, 2)})",
        lambda: StructArray(
            (1, 2),
            ("x", "e"),
            [{"x": 1.0, "e": {}}, {"x": 2.0, "e": StructArray.empty((), (0, 2))}],
        ),
    ),
    "deep": ("{struct('c', {{1, {struct()}}})}", lambda: [{"c": [1.0, [{}]]}]),
    "nestf": (
        "struct('e', struct('b', {}, 'a', {}), 'u', 1)",
        lambda: {"e": StructArray.empty(("b", "a")), "u": 1.0},
    ),
}
# without a reference an empty list is [] (to_mvalue's rule), with one (a cell) an empty cell
NEEDS_REFERENCE = {"c00", "cc", "mixc", "nests"}
FILES = {"structs7.mat": "7", "structs6.mat": "6", "structs73.mat": "7.3"}
DEFINE = "".join(f"{n} = {e}; " for n, (e, _) in VARIABLES.items())
NAMES = ", ".join(f"'{n}'" for n in VARIABLES)
SAVE = DEFINE + "".join(f"save('{f}', {NAMES}, '-v{v}'); " for f, v in FILES.items()) + "out = {};"
RETURN = DEFINE + "out = {" + ", ".join(VARIABLES) + "};"


# MATLAB's objects and function handles beside such arrays (the port writes neither): files
# whose layout the harness reads past objects and function handles
OBJECTS = {
    "o1": "{inputParser, struct(), {}, @sin, containers.Map, struct([])}",
    "o2": "struct('p', inputParser, 'e', struct(), 'c', {{}}, 'f', @sin)",
    "o3": "inputParser",
    "o4": "repmat(struct(), 1, 2)",
}
OBJECT_FILES = {"objects7.mat": "7", "objects6.mat": "6"}
SAVE_OBJECTS = (
    "".join(f"{n} = {e}; " for n, e in OBJECTS.items())
    + "".join(
        f"save('{f}', {', '.join(repr(n) for n in OBJECTS)}, '-v{v}'); "
        for f, v in OBJECT_FILES.items()
    )
    + "out = {};"
)


def write_structs() -> None:
    values = {n: make() for n, (_, make) in VARIABLES.items()}
    for f, v in FILES.items():
        write_mat(f, values, version=v)


def port_values() -> tuple:
    return tuple(make() for _, make in VARIABLES.values())


CASES = [
    Case("mat.struct_arrays", code=SAVE, python=write_structs, nargout=0, files=True),
    Case(
        "values.struct_arrays",
        code=RETURN,
        python=port_values,
        nargout=len(VARIABLES),
        outputs=tuple(VARIABLES),
    ),
    Case("mat.objects_beside_struct_arrays", code=SAVE_OBJECTS, nargout=0, files=True),
]


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


def stored(case):
    store = oracle.Store.from_env()
    if store is None or not store.has(case.key()):
        pytest.skip("needs the golden store with this module's references")
    return store.load(case.key())


@pytest.mark.oracle
def test_matlab_files_read_as_matlab_decoded_them():
    files = stored(CASES[0]).files
    for name in FILES:
        it = files.items[name][0]
        theirs, mine = read_pkv(it.decoded), read_mat(it.path)
        rep = compare(theirs, mine, path=name)
        assert rep.ok and not rep.notes, rep.text()
        assert field_orders(mine) == field_orders(theirs)


@pytest.mark.oracle
def test_struct_arrays_beside_matlabs_objects_read_as_matlab_decoded_them():
    # objects are not decoded by the harness; the cells and struct arrays beside them are
    files = stored(CASES[2]).files
    for name in OBJECT_FILES:
        it = files.items[name][0]
        theirs, mine = read_pkv(it.decoded).items[0], read_mat(it.path).items[0]
        pairs = [(theirs["o1"].items[i], mine["o1"].items[i]) for i in (1, 2, 5)]
        pairs += [(theirs["o2"].items[0][k], mine["o2"].items[0][k]) for k in ("e", "c")]
        pairs.append((theirs["o4"], mine["o4"]))
        assert [vhash(m) for _, m in pairs] == [vhash(t) for t, _ in pairs], name
        assert [type(m).__name__ for _, m in pairs] == ["MStruct", "MCell", "MStruct"] + [
            "MStruct",
            "MCell",
            "MStruct",
        ]


@pytest.mark.oracle
def test_the_ports_values_convert_to_what_the_readers_give_for_matlabs_files():
    files = stored(CASES[0]).files
    outputs = stored(CASES[1]).outputs
    read = {f: read_mat(files.items[f][0].path).items[0] for f in FILES}
    for i, (name, (_, make)) in enumerate(VARIABLES.items()):
        want = vhash(outputs[i])  # MATLAB's value, as MATLAB returned it
        assert {f: vhash(read[f][name]) for f in FILES} == dict.fromkeys(FILES, want), name
        assert vhash(to_mvalue(make(), like=outputs[i])) == want, name
        if name not in NEEDS_REFERENCE:
            assert vhash(to_mvalue(make())) == want, name


# --- files of the port and of SciPy, without MATLAB -----------------------------------------------


def read_one(tmp_path, value, version):
    path = tmp_path / f"v{version}.mat"
    write_mat(path, {"v": value}, version=version)
    return read_mat(path).items[0]["v"]


SIZES = [(1, 1), (1, 3), (3, 1), (3, 2), (2, 2, 2), (0, 0), (0, 3), (3, 0), (2, 0, 3)]


@pytest.mark.parametrize("version", ["7", "6"])
@pytest.mark.parametrize("size", SIZES)
def test_struct_arrays_without_fields_keep_their_size_and_elements(tmp_path, version, size):
    n = int(np.prod(size))
    v = read_one(tmp_path, {} if size == (1, 1) else StructArray(size, (), [{}] * n), version)
    assert isinstance(v, MStruct)
    assert (v.size, v.fields, v.items) == (size, [], [{}] * n)
    assert vhash(v) == vhash(MStruct(size, [], [{}] * n))


@pytest.mark.parametrize("version", ["7", "6"])
def test_empty_cells_stay_cells(tmp_path, version):
    v = read_one(tmp_path, [[], np.empty((0, 3), dtype=object), [[]], [None] * 2], version)
    c00, c03, cc, c12 = v.items
    assert [(type(x), x.size) for x in (c00, c03, cc, c12)] == [
        (MCell, (0, 0)),
        (MCell, (0, 3)),
        (MCell, (1, 1)),
        (MCell, (1, 2)),
    ]
    assert isinstance(cc.items[0], MCell) and cc.items[0].size == (0, 0)
    assert [(x.mclass, x.size) for x in c12.items] == [("double", (0, 0))] * 2


@pytest.mark.parametrize("version", ["7", "6"])
def test_struct_arrays_without_fields_nested_in_cells_and_structs(tmp_path, version):
    path = tmp_path / "n.mat"
    nested = ("mixc", "nests", "nesta", "deep", "nestf")
    write_mat(path, {n: VARIABLES[n][1]() for n in nested}, version=version)
    got = read_mat(path).items[0]
    kinds = [(type(x).__name__, x.size, getattr(x, "fields", None)) for x in got["mixc"].items]
    assert kinds == [
        ("MStruct", (1, 1), []),
        ("MCell", (0, 0), None),
        ("MStruct", (0, 0), []),
        ("MCell", (0, 2), None),
        ("MStruct", (2, 1), []),
        ("MCell", (1, 1), None),
    ]
    assert got["mixc"].items[5].items[0].fields == [] and got["mixc"].items[4].items == [{}, {}]
    s = got["nests"].items[0]
    assert [(type(s[k]).__name__, s[k].size) for k in ("e", "c", "n", "z")] == [
        ("MStruct", (1, 1)),
        ("MCell", (0, 0)),
        ("MStruct", (1, 2)),
        ("MStruct", (0, 0)),
    ]
    e = [el["e"] for el in got["nesta"].items]
    assert [(type(x).__name__, x.size, x.fields) for x in e] == [
        ("MStruct", (1, 1), []),
        ("MStruct", (0, 2), []),
    ]
    inner = got["deep"].items[0].items[0]["c"].items[1]
    assert isinstance(inner, MCell) and isinstance(inner.items[0], MStruct)
    assert got["nestf"].items[0]["e"].fields == ["b", "a"]
    # each is what to_mvalue makes of the port's value, and what the v7.3 reader gives
    write_mat(tmp_path / "n73.mat", {n: VARIABLES[n][1]() for n in nested}, version="7.3")
    got73 = read_mat(tmp_path / "n73.mat").items[0]
    for n in nested:
        mine = to_mvalue(VARIABLES[n][1](), like=got[n])
        assert vhash(got[n]) == vhash(mine) == vhash(got73[n]), n


def test_a_value_without_empty_or_none_object_arrays_reads_without_the_layout(
    tmp_path, monkeypatch
):
    def layout(path, names):
        raise AssertionError(f"layout read for {names}")

    # numbers, a cell of values and an empty struct array with fields: SciPy reads them right
    p, q = tmp_path / "p.mat", tmp_path / "q.mat"
    write_mat(
        p,
        {"x": np.arange(3.0), "s": {"a": [1.0, "b"]}, "t": StructArray.empty(("a",))},
        version="7",
    )
    write_mat(q, {"x": 1.0, "c": [[]]}, version="7")
    monkeypatch.setattr(C, "_mat5_layout", layout)
    got = read_mat(p).items[0]
    assert (got["t"].size, got["t"].fields) == ((0, 0), ["a"])
    with pytest.raises(AssertionError, match=r"layout read for \['c'\]"):
        read_mat(q)


def test_scipys_files_are_read_with_their_classes(tmp_path):
    from scipy.io import savemat

    cell = np.empty((1, 2), dtype=object)
    cell[0, 0], cell[0, 1] = np.zeros((0, 0)), np.empty((0, 0), dtype=object)
    savemat(tmp_path / "s.mat", {"c": cell, "e": np.empty((0, 3), dtype=object), "s": {"a": 1.0}})
    got = read_mat(tmp_path / "s.mat").items[0]
    assert [type(x).__name__ for x in got["c"].items] == ["MArray", "MCell"]
    assert (type(got["e"]).__name__, got["e"].size) == ("MCell", (0, 3))
    assert got["s"].fields == ["a"]


class Mat5:
    """MAT 5 (v6) elements in a byte order (">": as a big-endian machine writes them)"""

    def __init__(self, o: str):
        self.o = o

    def el(self, t: int, data: bytes) -> bytes:
        if 0 < len(data) <= 4:  # small data element: one word, the byte count in its upper half
            return struct.pack(self.o + "I", len(data) << 16 | t) + data.ljust(4, b"\x00")
        return struct.pack(self.o + "II", t, len(data)) + data + bytes(-len(data) % 8)

    def matrix(self, cls: int, dims, name: bytes, rest: bytes = b"") -> bytes:
        o = self.o
        body = self.el(6, struct.pack(o + "II", cls, 0)) + self.el(5, struct.pack(o + "2i", *dims))
        return self.el(14, body + self.el(1, name) + rest)

    def no_fields(self, dims) -> bytes:
        """struct array without fields, as MATLAB writes it (field name length 1, no names)"""
        return self.matrix(2, dims, b"", self.el(5, struct.pack(self.o + "i", 1)) + self.el(1, b""))

    def write(self, path, *variables: bytes) -> None:
        head = b"MATLAB 5.0 MAT-file".ljust(116) + bytes(8)
        head += b"\x01\x00MI" if self.o == ">" else b"\x00\x01IM"
        with open(path, "wb") as f:
            f.write(head + b"".join(variables))


def test_a_big_endian_file_is_read_in_its_byte_order(tmp_path):
    m = Mat5(">")
    m.write(
        tmp_path / "b.mat",
        m.matrix(1, (1, 2), b"c", m.no_fields((1, 1)) + m.matrix(1, (0, 0), b"")),
    )
    got = read_mat(tmp_path / "b.mat").items[0]["c"]
    assert [(type(x).__name__, x.size) for x in got.items] == [
        ("MStruct", (1, 1)),
        ("MCell", (0, 0)),
    ]


def test_elements_written_without_data_read_as_empty_doubles(tmp_path):
    # MATLAB writes the elements of cell(1, 3) (never assigned) as miMATRIX elements of no
    # bytes, which SciPy reads as 1x0 doubles and MATLAB's load as []; in a struct's field too
    m = Mat5("<")
    nothing = struct.pack("<II", 14, 0)
    s = m.matrix(2, (1, 1), b"s", m.el(5, struct.pack("<i", 2)) + m.el(1, b"a\x00") + nothing)
    m.write(tmp_path / "e.mat", m.matrix(1, (1, 3), b"c", nothing * 3), s)
    got = read_mat(tmp_path / "e.mat").items[0]
    assert [(x.mclass, x.size) for x in got["c"].items] == [("double", (0, 0))] * 3
    assert (got["s"].items[0]["a"].mclass, got["s"].items[0]["a"].size) == ("double", (0, 0))
    # a 1x0 double that the file holds as one stays 1x0
    write_mat(tmp_path / "z.mat", {"c": [np.zeros((1, 0)), []]}, version="7")
    z = read_mat(tmp_path / "z.mat").items[0]["c"]
    assert [(type(x).__name__, x.size) for x in z.items] == [("MArray", (1, 0)), ("MCell", (0, 0))]


def test_a_layout_that_does_not_match_scipys_value_is_an_error():
    none = np.empty((1, 2), dtype=object)
    with pytest.raises(ValueError, match=r"not the file's struct \[1 3\]"):
        C._from_scipy(none, (2, (1, 3), 0, None))
    filled = np.empty((1, 2), dtype=object)
    filled[0, 0] = np.zeros((1, 1))
    with pytest.raises(ValueError, match=r"not the file's struct \[1 2\]"):
        C._from_scipy(filled, (2, (1, 2), 0, None))
    with pytest.raises(ValueError, match=r"not the file's cell \[1 2\]"):
        C._from_scipy(np.zeros((1, 2)), (1, (1, 2), 1, [None, (2, (1, 1), 0, None)]))
    rec = np.empty((1, 1), dtype=[("a", object)])
    with pytest.raises(ValueError, match=r"not the file's struct \[1 1\]"):
        C._from_scipy(rec, (2, (1, 1), 2, [None, (2, (1, 1), 0, None)]))
    with pytest.raises(ValueError, match="not the file's element of no bytes"):
        C._from_scipy(np.zeros((1, 2)), C._NO_DATA)
    with pytest.raises(ValueError, match="not the file's element of no bytes"):
        C._from_scipy(np.empty((0, 0), dtype=object), C._NO_DATA)


def test_cell_elements_never_assigned_read_as_empty_doubles(tmp_path):
    # MATLAB's v7.3 layout of cell(1, 3): three references to one dataset, uint64 [0 0] of
    # class "canonical empty" (MATLAB_empty 1), which MATLAB's load gives as []
    import h5py

    path = tmp_path / "c.mat"
    with h5py.File(path, "w") as f:
        a = f.create_group("#refs#").create_dataset("a", data=np.zeros(2, np.uint64))
        a.attrs["MATLAB_class"] = np.bytes_("canonical empty")
        a.attrs.create("MATLAB_empty", np.uint8(1))
        c = f.create_dataset("c", data=np.array([[a.ref]] * 3, dtype=h5py.ref_dtype))
        c.attrs["MATLAB_class"] = np.bytes_("cell")
    got = read_mat(path).items[0]["c"]
    assert got.size == (1, 3)
    assert [(x.mclass, x.size) for x in got.items] == [("double", (0, 0))] * 3


def test_a_cell_element_that_is_not_a_matrix_is_an_error():
    # a cell's elements are miMATRIX elements (14); anything else is not a MAT-file's layout
    body = (
        struct.pack("<II", 6, 8)
        + struct.pack("<II", 1, 0)
        + struct.pack("<II", 5, 8)
        + struct.pack("<2i", 1, 1)
        + struct.pack("<II", 1, 0)
        + struct.pack("<II", 9, 8)
        + struct.pack("<d", 1.0)
    )
    with pytest.raises(ValueError, match="element type 9 in a cell or struct"):
        C._mat5_node(memoryview(body), "<")


# --- to_mvalue of StructArray ---------------------------------------------------------------------


@pytest.mark.parametrize("size", SIZES)
@pytest.mark.parametrize("fields", [(), ("b", "a")])
def test_struct_arrays_convert_to_their_size_and_fields(size, fields):
    n = int(np.prod(size))
    rows = [{k: float(i) for k in fields} for i in range(n)]
    v = to_mvalue(StructArray(size, fields, rows))
    assert isinstance(v, MStruct)
    assert (v.size, v.fields, len(v.items)) == (size, list(fields), n)
    want = MStruct(size, list(fields), [{k: to_mvalue(r[k]) for k in fields} for r in rows])
    assert vhash(v) == vhash(want)


def test_a_struct_array_with_trailing_singleton_dimensions_has_matlabs_size():
    assert to_mvalue(StructArray((2, 3, 1), (), [{}] * 6)).size == (2, 3)


def test_a_1d_struct_array_is_a_row_or_the_references_column():
    s = StructArray((3,), ("a",), [{"a": 1.0}] * 3)
    assert to_mvalue(s).size == (1, 3)
    col = to_mvalue(StructArray((3, 1), ("a",), [{"a": 1.0}] * 3))
    assert to_mvalue(s, like=col).size == (3, 1)
    assert to_mvalue(StructArray((0,), ()), like=MStruct((0, 3), [], [])).size == (0, 3)


def test_a_struct_array_reads_its_elements_as_the_reference():
    s = StructArray((1, 2), ("a", "n"), [{"a": [1.0, 2.0], "n": 3}, {"a": [], "n": 4}])
    ref = MStruct(
        (1, 2),
        ["a", "n"],
        [
            {
                "a": MCell((1, 2), [to_mvalue(1.0), to_mvalue(2.0)]),
                "n": MArray("uint8", np.full((1, 1), 3, np.uint8)),
            },
            {"a": MCell((0, 0), []), "n": MArray("uint8", np.full((1, 1), 4, np.uint8))},
        ],
    )
    v = to_mvalue(s, like=ref)
    assert vhash(v) == vhash(ref)
    assert compare(ref, s).ok
    plain = to_mvalue(s)  # without the reference: a double row, [] and doubles
    assert [(el["a"].mclass, el["n"].mclass) for el in plain.items] == [("double", "double")] * 2


def test_a_struct_array_is_its_list_of_dicts():
    rows = [{"a": 1.0, "b": "x"}, {"a": 2.0, "b": "y"}]
    like = MStruct((1, 2), ["a", "b"], None)
    assert vhash(to_mvalue(StructArray.from_list(rows), like=like)) == vhash(to_mvalue(rows, like))


def test_struct_arrays_in_dicts_lists_cells_and_struct_arrays():
    inner = StructArray((1, 2), (), [{}, {}])
    cell = np.empty((2, 1), dtype=object)
    cell[0, 0], cell[1, 0] = inner, StructArray.empty(("q",))
    v = to_mvalue(
        StructArray((1, 2), ("s",), [{"s": {"t": inner}}, {"s": [inner, cell]}]),
    )
    assert v.items[0]["s"].items[0]["t"].size == (1, 2)
    lst = v.items[1]["s"]
    assert (type(lst).__name__, lst.items[0].size, lst.items[1].size) == ("MCell", (1, 2), (2, 1))
    assert (lst.items[1].items[1].size, lst.items[1].items[1].fields) == ((0, 0), ["q"])


def test_a_hash_only_reference_reads_a_struct_array_as_the_full_one():
    s = StructArray((1, 300), ("a",), [{"a": float(i)} for i in range(300)])
    full = to_mvalue(s)
    d = describe(full)
    assert "v" not in d and "r" in d  # the manifest keeps the hash and the reading signature
    ref = from_descriptor(d)
    assert isinstance(ref, MHashed)
    assert vhash(to_mvalue(s, like=ref)) == vhash(full)
    assert compare(ref, s).ok


def test_a_struct_array_of_another_size_or_fields_is_a_difference():
    ref = to_mvalue(StructArray((1, 3), ("a",), [{"a": 1.0}] * 3))
    rep = compare(ref, StructArray((3, 1), ("a",), [{"a": 1.0}] * 3))
    assert not rep.ok and "size" in rep.text()
    rep = compare(ref, StructArray((1, 3), ("b",), [{"b": 1.0}] * 3))
    assert not rep.ok and "field" in rep.text()
