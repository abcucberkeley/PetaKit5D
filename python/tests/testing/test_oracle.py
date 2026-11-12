"""Tests of the oracle harness: its own cases against MATLAB, the store, manifests and packs.

The cases here check the harness end to end: values through MATLAB and back, an array
function, a function that writes files, functions that read input files, expected errors and
warnings. Their references are in tests/oracle_matlab/selftest.oracle.json (CI) and in the
golden store (cluster). Tests that start MATLAB themselves need PETAKIT5D_ORACLE_DIR and MATLAB.
"""

import json
import os
import shutil
import signal
import subprocess
import threading
import time
import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d.testing import oracle
from petakit5d.testing.cases import (
    FORMAT,
    Case,
    File,
    In,
    InPath,
    Tree,
    Value,
    Work,
    changed_files,
    mcell,
    mcol,
    mempty,
    mstring,
    mstruct,
    snapshot,
)
from petakit5d.testing.compare import (
    MCell,
    MStruct,
    compare,
    ignore,
    mchar,
    read_mat,
    read_pkv,
    to_mvalue,
    vhash,
    write_pkv,
)

ORACLE_MANIFEST = Path(__file__).resolve().parents[1] / "oracle_matlab" / "selftest.oracle.json"
MATLAB_NAN = np.array([0xFFF8000000000000], np.uint64).view(np.float64)[0]


def vol(shape, dtype=np.uint16, seed=0):
    """deterministic test volume (no random generator: the same bytes on every platform)"""
    a = (np.arange(int(np.prod(shape)), dtype=np.int64) * 7919 + seed * 104729) % 4001
    return np.asfortranarray(a.reshape(shape, order="F").astype(dtype))


class SelfTestError(RuntimeError):
    def __init__(self, identifier, message):
        super().__init__(message)
        self.identifier = identifier


# --- values through MATLAB and back --------------------------------------------------------------

ROUND_TRIP = {
    "nan_inf_zero": lambda: np.array([MATLAB_NAN, np.nan, -0.0, np.inf, -np.inf, 5e-324]),
    "int64_extremes": lambda: np.array(
        [np.iinfo(np.int64).min, np.iinfo(np.int64).max, -1], np.int64
    ),
    "uint64_max": lambda: np.array([[2**64 - 1], [0]], np.uint64),
    "single_4d": lambda: vol((2, 3, 4, 5), np.float32, 3),
    "uint16_3d_odd": lambda: vol((7, 9, 5), np.uint16, 4),
    "int8_row": lambda: np.arange(-5, 5, dtype=np.int8),
    "empty_0x3_uint16": lambda: np.zeros((0, 3), np.uint16),
    "empty_1x0x2_single": lambda: np.zeros((1, 0, 2), np.float32),
    "complex64": lambda: np.array([[1.5 - 0j, complex(0, -0.0)], [2j, -1]], np.complex64),
    "logical": lambda: np.array([[True, False, True]]),
    "char_unicode": lambda: "héllo \U0001f600",
    "char_empty": lambda: "",
    "struct_nested": lambda: {"a": [1, 2], "b": {"c": "x", "d": None}, "e": True},
    "struct_array": lambda: MStruct(
        (2, 1),
        ["p", "q"],
        [
            {"p": to_mvalue(1), "q": mchar("a")},
            {"p": to_mvalue("two"), "q": MCell((1, 1), [to_mvalue(3)])},
        ],
    ),
    "cell_mixed": lambda: mcell(1, "two", [3, 4], mcell()),
    "string_array": lambda: mstring("a", None, ""),
    "column": lambda: mcol([1, 2, 3], "int32"),
    "empty_struct_fields": lambda: MStruct((0, 1), ["a", "b"], []),
    "empty_char_1x0": lambda: mempty("char", (1, 0)),
    "struct_literal": lambda: mstruct(x=1.5, y="z"),
}

CASES = [
    Case(
        f"identity_{name}",
        matlab="pkpy_oracle.selftest.identity",
        args=[In("x")],
        python=lambda x: x,
        inputs={"x": Value(make)},
    )
    for name, make in ROUND_TRIP.items()
]
CASES.append(
    Case(
        "identity_three_outputs",
        matlab="pkpy_oracle.selftest.identity",
        args=[1, "b", [True]],
        nargout=3,
        python=lambda *a: a,
    )
)


# --- MATLAB values Python has no literal for -----------------------------------------------------


def py_make_values():
    import scipy.sparse as sps

    return (
        "hello",
        ["a", None, "ccc"],
        {"b": 1, "a": [2, 3]},
        None,
        {"alpha": 5, "name": "x"},
        sps.csc_matrix((np.array([1.5, -2.0]), (np.array([0, 2]), np.array([1, 1]))), shape=(3, 4)),
        [{"a": 1, "b": None}, {"a": "two", "b": [3]}],
        np.array([np.iinfo(np.int64).min, np.iinfo(np.int64).max, -1], np.int64),
        np.array([MATLAB_NAN, -0.0, np.inf, -np.inf, 1e-310, MATLAB_NAN]),
        np.complex64(complex(1.5, -0.0)),
        np.array([["a", "b"], ["c", "d"]]),
        np.array([[True, False], [False, True]]),
        None,
        [],
    )


CASES.append(
    Case(
        "make_values",
        matlab="pkpy_oracle.selftest.make_values",
        nargout=14,
        python=py_make_values,
        outputs=(
            "str",
            "strs",
            "map",
            "fh",
            "ip",
            "sp",
            "sa",
            "i64",
            "specials",
            "cs",
            "chm",
            "lg",
            "es",
            "ec",
        ),
        compare={"fh": ignore("a MATLAB function handle has no Python counterpart")},
    )
)


# --- array functions (PetaKit5D mex and MATLAB built-ins) ----------------------------------------


def py_crop(a, bbox):
    b = [int(x) for x in bbox]
    return a[b[0] - 1 : b[3], b[1] - 1 : b[4], b[2] - 1 : b[5]]


def py_sort3(x):
    i = np.argsort(x, axis=2, kind="stable")
    return np.take_along_axis(x, i, axis=2), (i + 1).astype(np.float64)


def py_max3(x):
    i = np.argmax(x, axis=2)
    return x.max(axis=2), (i + 1).astype(np.float64)


CASES += [
    Case(
        "crop3d_mex_uint16",
        matlab="crop3d_mex",
        args=[In("im"), [2, 3, 1, 9, 10, 7]],
        python=py_crop,
        inputs={"im": Value(lambda: vol((10, 12, 8), np.uint16, 1))},
    ),
    Case(
        "crop3d_mex_single",
        matlab="crop3d_mex",
        args=[In("im"), [1, 1, 2, 5, 3, 3]],
        python=py_crop,
        inputs={"im": Value(lambda: vol((5, 3, 4), np.float32, 2))},
    ),
    Case(
        "cumsum_dim2_double",
        matlab="cumsum",
        args=[In("x"), 2],
        python=lambda x, d: np.cumsum(x, axis=d - 1),
        inputs={"x": Value(lambda: vol((6, 50, 3), np.float64, 5) / 7)},
    ),
    Case(
        "sort_dim3_single",
        matlab="sort",
        args=[In("x"), 3],
        nargout=2,
        python=lambda x, d: py_sort3(x),
        inputs={"x": Value(lambda: vol((4, 5, 9), np.float32, 6) % 13)},
    ),
    Case(
        "max_dim3_uint16_ties",
        code="[m, i] = max(in{1}, [], 3); out = {m, i};",
        nargout=2,
        args=[In("x")],
        python=py_max3,
        inputs={"x": Value(lambda: vol((4, 5, 9), np.uint16, 7) % 5)},
    ),
    Case(
        "code_numel_class",
        code="out = {numel(in{1}), class(in{1})};",
        nargout=2,
        args=[In("x")],
        python=lambda x: (x.size, x.dtype.name),
        inputs={"x": Value(lambda: vol((3, 4), np.uint16))},
    ),
]


# --- a function that writes files ----------------------------------------------------------------


def py_write_files(out_path, im):
    import uuid

    import cpptiff
    import cppzarr
    from scipy.io import savemat

    out = Path(out_path)
    (out / "sub").mkdir()
    cpptiff.write_tiff(
        str(out / "sub" / "im.tif"), np.ascontiguousarray(np.transpose(im, (2, 0, 1)))
    )
    cppzarr.write_zarr(
        str(out / "im.zarr"),
        np.asfortranarray(im),
        chunks=list(im.shape),
        order="F",
        dimension_separator=".",
    )
    s = {
        "name": "demo",
        "sz": np.array([im.shape], dtype=np.float64),
        "vals": np.array([1.0, "a"], dtype=object),
        "t": 2.5,
    }
    savemat(out / "info.mat", s)
    savemat(out / "info_v6.mat", s)
    (out / "params.json").write_text(
        json.dumps({"a": 1, "b": "x", "c": [1, 2]}, separators=(",", ":"))
    )
    (out / "log.txt").write_text(f"line 1\nwork {out_path}\n")
    (out / "empty").mkdir()
    flag = f"flag_{uuid.uuid4()}.txt"
    (out / flag).write_text("")
    return ["sub/im.tif", "im.zarr", "info.mat", "info_v6.mat", "params.json", "log.txt", flag]


CASES.append(
    Case(
        "write_files",
        matlab="pkpy_oracle.selftest.write_files",
        args=[Work(), In("im")],
        python=py_write_files,
        inputs={"im": Value(lambda: vol((10, 12, 8), np.uint16, 2))},
    )
)


# --- functions that read input files -------------------------------------------------------------


def make_data_tree(d):
    import cpptiff
    import cppzarr

    cpptiff.write_tiff(
        str(d / "a.tif"),
        np.ascontiguousarray(np.transpose(vol((6, 7, 3), np.uint16, 8), (2, 0, 1))),
    )
    cppzarr.write_zarr(
        str(d / "b.zarr"),
        vol((5, 4, 3), np.float32, 9),
        chunks=[3, 3, 3],
        order="F",
        dimension_separator=".",
    )
    (d / "sub").mkdir()
    (d / "sub" / "c.txt").write_text("c\n")


def make_mat(p):
    # inputs must be the same bytes every time: savemat writes the time into the header
    from scipy.io import savemat

    savemat(p, {"x": np.arange(6.0).reshape(2, 3), "s": "hi"})
    b = p.read_bytes()
    p.write_bytes(b"MATLAB 5.0 MAT-file, PetaKit5D test input".ljust(116) + b[116:])


def py_read_inputs(data_path, mat_fn):
    import cpptiff
    import cppzarr
    from scipy.io import loadmat

    root = Path(data_path)
    names = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())
    tif = np.transpose(cpptiff.read_tiff(str(root / "a.tif")), (1, 2, 0))
    zr = cppzarr.read_zarr(str(root / "b.zarr"))
    m = loadmat(mat_fn)
    return names, tif, zr, {"s": str(m["s"][0]), "x": m["x"]}


CASES.append(
    Case(
        "read_inputs",
        matlab="pkpy_oracle.selftest.read_inputs",
        args=[InPath("data"), InPath("m")],
        nargout=4,
        python=py_read_inputs,
        inputs={"data": Tree(make_data_tree, "data"), "m": File(make_mat, "m.mat")},
    )
)


# --- expected errors and warnings ----------------------------------------------------------------


def py_fail(n):
    warnings.warn("about to fail", stacklevel=2)
    raise SelfTestError("pkpy:selftest:fail", f"failed on purpose with {n}")


def py_readtiff(fn):
    raise FileNotFoundError(fn)


def py_warn(x):
    for ident, msg in (
        ("pkpy:selftest:first", "first warning"),
        ("pkpy:selftest:second", f"second warning: {x.size} values"),
    ):
        w = UserWarning(msg)
        w.identifier = ident
        warnings.warn(w, stacklevel=2)
    return x * 2


CASES += [
    Case(
        "expected_error",
        matlab="pkpy_oracle.selftest.fail",
        args=[3],
        python=py_fail,
        expect_error="pkpy:selftest:fail",
        py_error=SelfTestError,
        error_message=True,
    ),
    Case(
        "readtiff_missing_file",
        matlab="readtiff",
        args=[Work("missing.tif")],
        python=py_readtiff,
        expect_error=True,
        py_error=FileNotFoundError,
    ),
    Case(
        "warnings",
        matlab="pkpy_oracle.selftest.warn",
        args=[In("x")],
        python=py_warn,
        warnings="ids",
        inputs={"x": Value(lambda: np.arange(4.0))},
    ),
]

# --- an input whose key is not its name; warnings in all the forms MATLAB prints -----------------


def make_csv(p):
    p.write_text("a b,c-d,e\n1,2,3\n4,5,6\n")  # headers readtable must rename: it warns


def py_warn_forms(csv):
    for ident, msg in (
        ("pkpy:selftest:long", "long " + "words " * 20),
        ("MATLAB:table:ModifiedAndSavedVarnames", "Column headers from the file were modified"),
        ("pkpy:selftest:quiet", "quiet [x]"),
    ):
        w = UserWarning(msg)
        w.identifier = ident
        warnings.warn(w, stacklevel=2)
    return float(len(Path(csv).read_text().splitlines()) - 1)


CASES += [
    Case(
        "inpath_by_name",
        matlab="pkpy_oracle.selftest.read_inputs",
        args=[InPath("tiles"), InPath("m")],
        nargout=4,
        python=py_read_inputs,
        inputs={"tiles": Tree(make_data_tree, "raw"), "m": File(make_mat, "m.mat")},
    ),
    Case(
        "warnings_wrapped_readtable_quiet",
        code=(
            "warning('pkpy:selftest:long', ['long ', repmat('words ', 1, 20)]); "
            "T = readtable(in{1}); w = warning('off', 'verbose'); "
            "warning('pkpy:selftest:quiet', 'quiet [x]'); warning(w); out = {height(T)};"
        ),
        args=[InPath("csv")],
        python=py_warn_forms,
        warnings="ids",
        inputs={"csv": File(make_csv, "t.csv")},
    ),
]

# --- 100 small cases (harness overhead per case) -------------------------------------------------

SMALL = [
    Case(
        f"small_plus_{i:03d}",
        matlab="plus",
        args=[In("a"), In("b")],
        python=np.add,
        inputs={
            "a": Value(lambda i=i: vol((8, 9, 3), np.float64, i) / 3),
            "b": Value(lambda i=i: vol((8, 9, 3), np.float64, i + 1000) / 7),
        },
    )
    for i in range(100)
]
CASES += SMALL

# --- a GPU pack (orun.sh -g): MATLAB's gpuArray; the case is marked gpu ---------------------------

CASES.append(
    Case(
        "gpu_times_two",
        code="out = {gather(gpuArray(in{1}) * 2)};",
        args=[In("x")],
        python=lambda x: x * 2,
        device="gpu",
        inputs={"x": Value(lambda: vol((16, 16, 4), np.float32, 11))},
    )
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- harness logic (no MATLAB) -------------------------------------------------------------------

BS = "\x08"
# command window text as MATLAB R2024b's evalc gives it (from a probe): readtable's warning,
# wrapped at 80 columns also inside the "Type" line; warnings issued with verbose off
READTABLE_TEXT = (
    f"[{BS}Warning: Column headers from the file were modified to make them valid MATLAB\n"
    "identifiers before creating variable names for the table. The original column\n"
    "headers are saved in the VariableDescriptions property.\n"
    "Set 'VariableNamingRule' to 'preserve' to use the original column headers as\n"
    "table variable names.\n"
    '(Type "<a href="matlab:warning off MATLAB:table:ModifiedAndSavedVarnames">warning off '
    'MATLAB:table:ModifiedAndSavedVarnames</a>" to suppress this\nwarning.)\n'
    f"]{BS} \n"
)
QUIET_TEXT = f"[{BS}Warning: quiet one]{BS} \n[{BS}Warning: quiet two [x]]{BS} \n"


def test_matlab_warnings_are_parsed_in_every_form():
    w = oracle.parse_matlab_warnings(READTABLE_TEXT)
    assert [x["identifier"] for x in w] == ["MATLAB:table:ModifiedAndSavedVarnames"]
    assert w[0]["message"].startswith("Column headers") and w[0]["message"].endswith("names.")
    w = oracle.parse_matlab_warnings(QUIET_TEXT, "pk:b", "quiet two [x]")
    assert w == [
        {"identifier": "", "message": "quiet one"},  # verbose off and not the last: no id
        {"identifier": "pk:b", "message": "quiet two [x]"},  # the last: lastwarn's
    ]
    exact = "Column headers from the file were modified to make them valid MATLAB identifiers"
    w = oracle.parse_matlab_warnings(
        READTABLE_TEXT + QUIET_TEXT, "MATLAB:table:ModifiedAndSavedVarnames", exact
    )
    assert [x["identifier"] for x in w] == ["MATLAB:table:ModifiedAndSavedVarnames", "", ""]
    plain = '[Warning: plain\n(Type "warning off a:b" to suppress this warning.)\n]\n'
    w = oracle.parse_matlab_warnings(plain)
    assert w == [{"identifier": "a:b", "message": "plain"}]
    w = oracle.parse_matlab_warnings(f"[{BS}Warning: cut off")
    assert w[-1]["identifier"] == "pkpy:oracle:unparsedWarning"
    # a warning without an identifier, wrapped: the last one's exact text is lastwarn's
    noid = (
        f"[{BS}Warning: Setting file\n<WORK>/s/x.txt does not\nexist, skip it.\n"
        '(Type "<a href="matlab:warning off ">warning off </a>" to suppress this warning.)\n'
        f"]{BS} \n"
    )
    exact = "Setting file <WORK>/s/x.txt does not exist, skip it."
    assert oracle.parse_matlab_warnings(noid, "", exact) == [{"identifier": "", "message": exact}]
    assert oracle.parse_matlab_warnings(noid, "", "another")[0]["message"] != "another"


def test_inpath_resolves_the_input_name_on_both_sides(tmp_path):
    case = next(c for c in CASES if c.id == "inpath_by_name")
    args = case.matlab_args(lambda name, x: "")
    assert [(a.what, a.name) for a in args.items] == [("path", "raw"), ("path", "m.mat")]
    py_args, _ = case.python_args(tmp_path)
    assert py_args == [tmp_path.as_posix() + "/raw", tmp_path.as_posix() + "/m.mat"]


def test_case_validation():
    with pytest.raises(ValueError):
        Case("a b", matlab="f")
    with pytest.raises(ValueError):
        Case("x", matlab="f", code="out = {};")
    with pytest.raises(ValueError):
        Case("x", matlab="f", args=[In("nope")])
    with pytest.raises(ValueError):
        Case("x", matlab="f", args=[In("m")], inputs={"m": File(make_mat, "m.mat")})
    with pytest.raises(ValueError):
        Case("x", matlab="f", nargout=2, outputs=("a",))


def test_key_depends_on_the_matlab_side_only():
    def mk(**kw):
        base = {
            "matlab": "plus",
            "args": [In("a"), 2],
            "python": np.add,
            "inputs": {"a": Value(lambda: np.arange(3.0))},
        }
        base.update(kw)
        return Case("k", **base)

    k = mk().key()
    assert mk(python=np.subtract, compare=ignore("x")).key() == k
    assert mk(args=[In("a"), 3]).key() != k
    assert mk(inputs={"a": Value(lambda: np.arange(4.0))}).key() != k
    assert mk(nargout=2).key() != k
    assert mk(matlab="minus").key() != k
    assert mk(device="gpu").key() != k


def test_changed_files_rules(tmp_path):
    w = tmp_path
    (w / "in.tif").write_bytes(b"1")
    (w / "z.zarr").mkdir()
    (w / "z.zarr" / ".zarray").write_text("{}")
    (w / "gone.txt").write_text("x")
    before = snapshot(w)
    (w / "out").mkdir()
    (w / "out" / "a.tif").write_bytes(b"2")
    (w / "z.zarr" / "0.0.0").write_bytes(b"c")
    (w / "new.zarr").mkdir()
    (w / "new.zarr" / ".zarray").write_text("{}")
    (w / "emptydir").mkdir()
    (w / "gone.txt").unlink()
    (w / "in.tif").write_bytes(b"12")  # an input the call changed
    files, deleted = changed_files(before, snapshot(w))
    assert files == ["emptydir", "in.tif", "new.zarr", "out/a.tif", "z.zarr"]
    assert deleted == ["gone.txt"]


def test_input_store_takes_concurrent_writers(tmp_path):
    store = oracle.Store(tmp_path / "golden", "5df18b227e5ac00e584e8274a783c3df5a09c97b")
    v = Value(lambda: vol((30, 30, 30), np.float32))
    t = Tree(make_text_tree, "d")
    paths, errors = [], []

    def put():
        try:
            paths.append(store.put_input(v))
            paths.append(store.put_input(t))
        except Exception as e:  # collected for the assertion below
            errors.append(e)

    ts = [threading.Thread(target=put) for _ in range(8)]
    [x.start() for x in ts]
    [x.join() for x in ts]
    assert not errors
    assert len(set(paths)) == 2
    assert vhash(read_pkv(store.put_input(v))) == v.hash()
    assert not list((tmp_path / "golden" / "inputs" / ".tmp").iterdir())


def make_text_tree(d):
    for i in range(20):
        (d / f"f{i}.txt").write_text(str(i) * 1000)


def fake_publish(store, case, outputs, status="ok", error=None, work="/scratch/pkpy_1/w/abc"):
    """write a case result as the MATLAB runner does (staging folder, then one rename)"""
    key = case.key()
    stage = store.staging("test") / key
    stage.mkdir(parents=True)
    if status == "ok":
        write_pkv(stage / "outputs.pkv", MCell((1, len(outputs)), [to_mvalue(x) for x in outputs]))
    meta = {
        "format": 1,
        "key": key,
        "id": case.id,
        "status": status,
        "error": error,
        "warnings": [],
        "capture": False,
        "files": [],
        "deleted": [],
        "work": work,
        "host": "node12",
        "oracle_commit": store.commit,
    }
    (stage / "meta.json").write_text(json.dumps(meta))
    final = store.case_dir(key)
    final.parent.mkdir(parents=True, exist_ok=True)
    os.rename(stage, final)


def test_reference_from_store_then_manifest(tmp_path, monkeypatch):
    commit = oracle.expected_commit()
    golden = tmp_path / "golden"
    golden.mkdir()
    monkeypatch.setenv(oracle.ENV_GOLDEN, str(golden))
    monkeypatch.setenv(oracle.ENV_COMMIT, commit)
    monkeypatch.delenv(oracle.ENV_ORACLE, raising=False)
    store = oracle.Store(golden, commit)
    mfile = tmp_path / "fake.oracle.json"
    case = Case(
        "fake_path",
        matlab="fullfile",
        args=[Work("a.tif")],
        python=lambda p: p,
        inputs={"x": Value(lambda: np.arange(3.0))},
    )
    case._module_globals = {
        "__file__": str(tmp_path / "test_fake.py"),
        "ORACLE_MANIFEST": str(mfile),
    }
    fake_publish(store, case, ["/scratch/pkpy_1/w/abc/a.tif"])
    with pytest.warns(UserWarning, match="manifest"):
        oracle.check(case, tmp_path / "run1")  # store; the work folder is masked
    oracle._REGISTRY[str(tmp_path / "test_fake.py")] = [case]
    oracle._write_manifests([("n", case)], store, None)
    entry = json.loads(mfile.read_text())["cases"]["fake_path"]
    assert entry["outputs"][0]["k"] == "a" and "/scratch" not in mfile.read_text()
    lines = mfile.read_text().splitlines()  # one line per case: a changed case is one line
    assert sum(ln.startswith('  "fake_path": {') for ln in lines) == 1 and len(lines) == 8
    monkeypatch.delenv(oracle.ENV_GOLDEN)
    oracle.check(case, tmp_path / "run2")  # manifest only, as in CI
    bad = Case(
        "fake_path",
        matlab="fullfile",
        args=[Work("a.tif")],
        python=lambda p: p + "x",
        inputs={"x": Value(lambda: np.arange(3.0))},
    )
    bad._module_globals = case._module_globals
    with pytest.raises(pytest.fail.Exception, match="text: MATLAB '<WORK>/a.tif'"):
        oracle.check(bad, tmp_path / "run3")
    doc = json.loads(mfile.read_text())
    doc["cases"]["fake_path"]["oracle"] = "0" * 40
    mfile.write_text(json.dumps(doc))
    with pytest.raises(pytest.fail.Exception, match="reference made with oracle 0000000"):
        oracle.check(case, tmp_path / "run4")
    doc["cases"]["fake_path"]["oracle"] = commit
    doc["cases"]["fake_path"]["key"] = "f" * 64
    mfile.write_text(json.dumps(doc))
    with pytest.raises(pytest.fail.Exception, match="changed since its reference was made"):
        oracle.check(case, tmp_path / "run5")
    assert oracle.check_manifests(tmp_path) == []


def test_check_manifests_flags_other_commits(tmp_path):
    (tmp_path / "a.oracle.json").write_text(
        json.dumps({"cases": {"c1": {"oracle": "1234567" * 6, "format": FORMAT}}})
    )
    probs = oracle.check_manifests(tmp_path)
    assert len(probs) == 1 and "made with oracle 1234567" in probs[0]


def test_repo_manifests_use_the_current_oracle():
    assert oracle.check_manifests() == []


def test_reference_errors_must_be_expected(tmp_path, monkeypatch):
    commit = oracle.expected_commit()
    monkeypatch.setenv(oracle.ENV_GOLDEN, str(tmp_path))
    monkeypatch.setenv(oracle.ENV_COMMIT, commit)
    store = oracle.Store(tmp_path, commit)
    case = Case("err", matlab="error", args=["a:b", "boom"], python=lambda *a: None)
    case._module_globals = {
        "__file__": str(tmp_path / "t.py"),
        "ORACLE_MANIFEST": str(tmp_path / "m.json"),
    }
    fake_publish(store, case, [], status="error", error={"identifier": "a:b", "message": "boom"})
    with pytest.warns(UserWarning, match="manifest"):
        with pytest.raises(pytest.fail.Exception, match="set expect_error=True"):
            oracle.check(case, tmp_path / "r")
    expected = Case(
        "err", matlab="error", args=["a:b", "boom"], python=lambda *a: None, expect_error="a:b"
    )
    expected._module_globals = case._module_globals
    with pytest.warns(UserWarning, match="manifest"):
        with pytest.raises(pytest.fail.Exception, match="Python returned normally"):
            oracle.check(expected, tmp_path / "r2")


# --- tests that need the golden store with the self-test references ------------------------------


def stored_reference(case_id):
    store = oracle.Store.from_env()
    case = next(c for c in CASES if c.id == case_id)
    if store is None or not store.has(case.key()):
        pytest.skip("needs the golden store with the self-test references")
    return store.load(case.key())


@pytest.mark.oracle
def test_mat_files_matlab_wrote_decode_as_matlab_reads_them():
    ref = stored_reference("write_files")
    for name in ("info.mat", "info_v6.mat"):
        it = ref.files.items[name][0]
        mine = read_mat(it.path)
        theirs = read_pkv(it.decoded)
        rep = compare(theirs, mine, path=name)
        assert rep.ok, rep.text()


# --- tests that start MATLAB ---------------------------------------------------------------------


@pytest.fixture
def matlab_oracle():
    od = oracle.oracle_dir()
    if (
        od is None
        or not od.is_dir()
        or shutil.which(os.environ.get(oracle.ENV_MATLAB, "matlab").split()[0]) is None
    ):
        pytest.skip("needs MATLAB and PETAKIT5D_ORACLE_DIR")
    return od


def private_store(tmp_path, monkeypatch):
    """a new, empty store; on network storage if PETAKIT5D_ORACLE_TEST_STORE names a folder there
    (the golden store's file system), else in the test's temporary folder"""
    base = os.environ.get("PETAKIT5D_ORACLE_TEST_STORE")
    golden = (Path(base) / f"{tmp_path.name}_{os.getpid()}") if base else tmp_path / "golden"
    golden.mkdir(parents=True)
    monkeypatch.setenv(oracle.ENV_GOLDEN, str(golden))
    return oracle.Store(golden, oracle.current_commit())


@pytest.mark.oracle
@pytest.mark.slow
def test_two_sessions_running_the_same_cases_publish_one_result(
    tmp_path, monkeypatch, matlab_oracle
):
    store = private_store(tmp_path, monkeypatch)
    cases = [
        Case(f"race_sleep_{i}", matlab="pkpy_oracle.selftest.sleep", args=[float(8 + 2 * i)])
        for i in range(2)
    ]
    cases += CASES[:3]
    packs = [
        oracle.write_pack(cases, store, tmp_path / f"pack{j}", tmp_path / f"work{j}")
        for j in range(2)
    ]
    for j in range(2):
        (tmp_path / f"work{j}").mkdir()
    rcs = [None, None]

    def run(j):
        rcs[j] = oracle.run_matlab(packs[j], matlab_oracle, echo=False)[0]

    ts = [threading.Thread(target=run, args=(j,)) for j in range(2)]
    [t.start() for t in ts]
    [t.join() for t in ts]
    assert rcs == [0, 0]
    rows = [
        json.loads((tmp_path / f"pack{j}" / "summary.json").read_text())["cases"] for j in range(2)
    ]
    for c in cases:
        k = c.key()
        st = sorted(r["status"] for rs in rows for r in rs if r["key"] == k)
        assert st[0] == "ok" or st[1] == "ok", st
        assert sum(s == "ok" for s in st) == 1, st  # one published, the other lost or found it
        assert set(st) <= {"ok", "cached", "lost_race"}
        assert len(list(store.case_dir(k).parent.glob(k + "*"))) == 1
        store.load(k)
    left = [p for p in store.dir.joinpath(".staging").rglob("*") if p.is_file()]
    assert left == []
    # the sessions start within seconds of each other, so both ran the sleeping cases at once
    sleeps = [r["status"] for rs in rows for r in rs if r["id"].startswith("race_sleep")]
    assert "lost_race" in sleeps, f"the sessions did not overlap ({sleeps}): no race was tested"


@pytest.mark.oracle
@pytest.mark.slow
def test_a_killed_session_leaves_no_partial_case(tmp_path, monkeypatch, matlab_oracle):
    store = private_store(tmp_path, monkeypatch)
    quick = CASES[0]
    slow = Case("killed_sleep", matlab="pkpy_oracle.selftest.sleep", args=[600.0])
    pj = oracle.write_pack([quick, slow], store, tmp_path / "pack", tmp_path / "work")
    (tmp_path / "work").mkdir()
    staged = Path(json.loads(pj.read_text())["staging"]) / slow.key()
    env = dict(os.environ)
    env.pop("DISPLAY", None)
    with open(tmp_path / "matlab.log", "w") as log:
        p = subprocess.Popen(
            [
                *os.environ.get(oracle.ENV_MATLAB, "matlab").split(),
                "-batch",
                oracle.matlab_command(pj, matlab_oracle),
            ],
            cwd=matlab_oracle,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            t0 = time.time()
            while (
                not (store.has(quick.key()) and staged.is_dir())
                and p.poll() is None
                and time.time() - t0 < 300
            ):
                time.sleep(0.5)
            assert p.poll() is None, (tmp_path / "matlab.log").read_text()
        finally:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
    assert store.has(quick.key())
    assert staged.is_dir()  # the stopped case's partial result ...
    assert not store.has(slow.key())  # ... is never visible in the store
    assert not store.case_dir(slow.key()).exists()
    store.load(quick.key())
