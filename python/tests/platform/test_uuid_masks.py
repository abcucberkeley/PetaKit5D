"""get_uuid's Windows form (its first 4 characters, as MATLAB there) in oracle comparisons: a
run's short uuids are aligned with the reference, made <uuid> only where the reference has
<uuid> at that place (temporary names and paths), so a real value equal to one (a 1024 in
parameters.json against the uuid "1024") compares as itself, with the golden store and with a
manifest alike. Linux's 36-character uuids are masked by UUID_RE on both sides."""

import functools
import hashlib
import json
import os
import sys
import types
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core import ids
from petakit5d.testing import oracle
from petakit5d.testing.cases import UUID_LOG_ENV, Case, PyResult, Work
from petakit5d.testing.compare import (
    Masks,
    MCell,
    MMap,
    MObject,
    MString,
    MStruct,
    Report,
    _cmp_zarr_meta,
    _masked_choices,
    align_json,
    align_name,
    align_text,
    align_value,
    compare,
    compare_trees,
    describe,
    from_descriptor,
    list_tree,
    mchar,
    scan_tree,
    short_uuids,
    to_mvalue,
    tree,
    uuid_masks,
    vhash,
    write_pkv,
)

FULL = "0f8fad5b-d9cb-469f-a165-70867728950e"  # MATLAB's uuid in the references (Linux)
U = ("1024",)


def sha(t):
    return hashlib.sha256(t.encode("utf-8")).hexdigest()


def test_only_the_runs_own_uuids_are_masked():
    m = Masks().plus(*uuid_masks(["6cac"]))
    assert m.text("zarr_flag/crop_6cac/writer.txt") == "zarr_flag/crop_<uuid>/writer.txt"
    assert m.text("im_6cac") == "im_<uuid>"
    for name in ("/mnt/projects/embryo_2024/cell1", "tile_0001.tif", "x_07d5", "ab6cacd"):
        assert m.text(name) == name


def test_full_uuids_are_left_to_the_usual_rule():
    assert uuid_masks([FULL]) == ()
    assert Masks().text(f"crop_{FULL}") == "crop_<uuid>"
    assert short_uuids([FULL, "6cac", "", "6cac", "0000"]) == ("0000", "6cac")
    assert Masks().with_uuids([FULL]).uuids == ()


def test_get_uuid_logs_what_it_gives(tmp_path, monkeypatch):
    log = tmp_path / "uuids.log"
    monkeypatch.setenv(UUID_LOG_ENV, str(log))
    monkeypatch.setattr(sys, "platform", "win32")  # the Windows form, also here
    a, b = ids.get_uuid(), ids.get_uuid()
    monkeypatch.undo()
    assert len(a) == 4 and len(b) == 4
    assert log.read_text().split() == [a, b]


def test_the_masks_leave_a_value_equal_to_the_runs_short_uuid_alone(tmp_path):
    # a Windows run whose uuid was "1024", against a reference value 1024 (in parameters.json):
    # the run's uuids go with its masks, to be aligned, and neither side substitutes them
    ref = oracle.Reference(
        key="k", status="ok", error=None, warnings=[], outputs=[], files=None, meta={},
        source="manifest", masks=Masks(), commit="",
    )  # fmt: skip
    res = PyResult(
        outputs=[], error=None, error_text="", warnings=[], files=None, work=tmp_path,
        elapsed=0.0, uuids=("1024", FULL),
    )  # fmt: skip
    rm, am = oracle.case_masks(ref, res)
    assert rm.text("    1024,") == am.text("    1024,") == "    1024,"
    assert (rm.uuids, am.uuids) == ((), U)
    assert am.plus(("x", "y")).uuids == am.with_paths(work="/w").uuids == U
    assert am.with_uuids(["0000"]).uuids == ("0000", "1024")


def test_a_uuid_is_masked_only_where_the_reference_has_one():
    ref = json.dumps({"batchSize": [13, 1024, 1024], "tmp": "<WORK>/x_<uuid>.tif"}, indent=2)
    act = ref.replace("<uuid>", "1024")
    assert align_text(ref, act, U) == ref  # the real values stay, the temporary name is masked
    assert align_text(ref, act, ()) == act  # no short uuids (Linux): nothing to align
    # the forms the port and MATLAB write a uuid in, and one glued to a name
    for form in ("uuid=<uuid>\n", "tmp/<uuid>/a.mat", "data_<uuid>_recon.tif", "ab<uuid>"):
        assert align_text(form, form.replace("<uuid>", "1024"), U) == form
    # MATLAB's real value where Python wrote its uuid, and the reverse, are differences
    assert align_text("    1024,", "    5f3a,", ("5f3a",)) == "    5f3a,"
    assert align_text("x_<uuid>.tif", "x_2048.tif", U) == "x_2048.tif"
    assert align_text("x_<uuid>.tif <uuid>", "x_1024.tif", U) == "x_1024.tif"
    # a text longer than the reference, or one that ends where a uuid should be, doesn't fit
    assert align_text("x_<uuid>.tif", "x_1024.tif.bak", U) == "x_1024.tif.bak"
    assert align_text("1024 x<uuid>", "1024 x", U) == "1024 x"
    # a text that differs beyond its uuids: aligned line by line, the real difference stays
    out = align_text(ref, act.replace("13,", "14,"), U)
    assert out == ref.replace("13,", "14,")
    # a uuid of the run against MATLAB's long one, and Linux's form masked by UUID_RE
    assert align_text("a <uuid>", Masks().text(f"a {FULL}"), U) == "a <uuid>"
    # uuids of several lengths, one a prefix of another (never get_uuid's): backtracked
    assert align_text("<uuid>12", "ab1212", ("ab", "ab12")) == "<uuid>12"


def test_a_text_known_by_its_hash_alone_takes_the_masking_that_gives_it():
    # a manifest keeps only the hash of a text over 4096 bytes; real 1024s beside the uuid 1024
    ref = "".join(f"task {i} of 1024: x_<uuid>.tif, 1024 done\n" for i in range(200))
    act = ref.replace("<uuid>", "1024")
    assert len(ref) > 4096
    h = sha(ref)
    assert align_text(None, act, U, lambda c: sha(c) == h) == ref
    assert align_text(None, act, U) == act  # no test of a candidate: as it is
    bad = act.replace("task 7 ", "task 8 ")
    assert sha(align_text(None, bad, U, lambda c: sha(c) == h)) != h  # a real difference stays
    # every choice of few occurrences is tried (all, none, in names, then by distance)
    assert len(set(_masked_choices("a 1 b_1 c 1 d=1", [(2, 3), (6, 7), (10, 11), (14, 15)]))) == 16
    assert align_text(None, "_1 1 1 _1", ("1",), lambda c: c == "_1 <uuid> 1 _1") == (
        "_1 <uuid> 1 _1"
    )
    # two of the run's uuids, one a real name's token (Iter_0000) on every line: by uuid
    for line in ("Iter_0000: x_<uuid>.tif\n", "Iter_0000: <uuid>\n", "Iter_0000: x_<uuid> 1024\n"):
        ref = line * 300
        act = ref.replace("<uuid>", "1024")
        same = functools.partial(lambda c, h: sha(c) == h, h=sha(ref))
        assert align_text(None, act, ("0000", "1024"), same) == ref


def test_names_are_aligned_with_the_references(tmp_path):
    names = {"d/Scan_Iter_0000_a.txt", "d/zarr_flag/crop_<uuid>/writer.txt", "d/x_0000_<uuid>.tif"}
    names |= {"d/y_0000.tif", "d/y_<uuid>.tif"}  # a real name and a temporary one
    u = ("0000",)
    assert align_name("d/Scan_Iter_0000_a.txt", names, u) == "d/Scan_Iter_0000_a.txt"
    assert align_name("d/zarr_flag/crop_0000/writer.txt", names, u) in names
    assert align_name("d/x_0000_0000.tif", names, u) == "d/x_0000_<uuid>.tif"
    assert align_name("d/y_0000.tif", names, u) == "d/y_0000.tif"  # itself first
    assert align_name("d/only_0000.tif", names, u) == "d/only_0000.tif"  # Python's own: as written
    # in a tree: the run's names as written, aligned when compared; removed files as they are
    for side, uid in (("ref", FULL), ("act", "0000")):
        for rel in ("Scan_Iter_0000_a.txt", f"zarr_flag/crop_{uid}/writer.txt"):
            (tmp_path / side / rel).parent.mkdir(parents=True, exist_ok=True)
            (tmp_path / side / rel).write_text(f"uuid={uid}\nbatch=0000\n", newline="\n")
    deleted = ["Iter_0000_in.tif"]
    ref = scan_tree(tmp_path / "ref", deleted=deleted)
    act = scan_tree(tmp_path / "act", deleted=deleted)
    rep = compare_trees(ref, act, act_masks=Masks().with_uuids(u))
    assert rep.ok, rep.text()
    assert not compare_trees(ref, act).ok  # the run's uuids are needed to match


def test_values_are_aligned_with_the_references():
    r = MCell((1, 3), [mchar("x_<uuid>.tif"), mchar("1024"), MString((1, 2), ["1024", "<uuid>"])])
    a = MCell((1, 3), [mchar("x_1024.tif"), mchar("1024"), MString((1, 2), ["1024", "1024"])])
    out = align_value(r, a, U)
    assert vhash(out) == vhash(r) and align_value(r, a, ()) is a
    s = to_mvalue({"file": "x_<uuid>.tif", "size": "1024"})
    t = align_value(s, to_mvalue({"file": "x_1024.tif", "size": "1024"}), U)
    assert isinstance(t, MStruct) and vhash(t) == vhash(s)
    # an object's value (an inputParser's Results against a struct), a containers.Map
    o = align_value(MObject("inputParser", s), to_mvalue({"file": "x_1024.tif", "size": "1024"}), U)
    assert vhash(o) == vhash(s)
    m = MMap("char", "any", MCell((1, 1), [mchar("1024")]), MCell((1, 1), [mchar("x_<uuid>")]))
    am = MMap("char", "any", MCell((1, 1), [mchar("1024")]), MCell((1, 1), [mchar("x_1024")]))
    assert vhash(align_value(m, am, U)) == vhash(m)
    # compare(): the uuids of the run's masks
    assert compare(mchar("x_<uuid>.tif"), "x_1024.tif", act_masks=Masks().with_uuids(U)).ok
    assert not compare(mchar("x_<uuid>.tif"), "x_1024.tif", act_masks=Masks()).ok
    # values a manifest keeps by their hash alone: a large cell, a long char row
    big = MCell((1, 400), [mchar("x_<uuid>.tif")] * 200 + [mchar("1024")] * 200)
    hb = from_descriptor(describe(big))
    act = MCell((1, 400), [mchar("x_1024.tif")] * 200 + [mchar("1024")] * 200)
    assert hb.kind == "hashed" and vhash(align_value(hb, act, U)) == vhash(big)
    two = MCell((1, 400), [mchar("Iter_0000_x.tif")] * 200 + [mchar("x_<uuid>.tif")] * 200)
    act = MCell((1, 400), [mchar("Iter_0000_x.tif")] * 200 + [mchar("x_1024.tif")] * 200)
    ht = from_descriptor(describe(two))  # two of the run's uuids, one a real name's: by uuid
    assert vhash(align_value(ht, act, ("0000", "1024"))) == vhash(two)
    row = "x_<uuid>.tif of 1024, " * 300
    hr = from_descriptor(describe(mchar(row)))
    assert hr.data is None and vhash(align_value(hr, mchar(row.replace("<uuid>", "1024")), U)) == (
        vhash(mchar(row))
    )
    two = mchar(row.replace("<uuid>", "1024").replace(" of", " in", 1))
    assert vhash(align_value(hr, two, U)) != vhash(mchar(row))  # a real difference stays


def test_json_is_aligned_by_its_content(tmp_path):
    # Python's JSON in another format (spacing): compared by content, its uuid aligned there
    for side, uid, ind in (("ref", FULL, None), ("act", "1024", 2)):
        (tmp_path / side).mkdir()
        p = {"batchSize": [13, 1024], "tmp": f"x_{uid}.tif", f"k_{uid}": 1024}
        (tmp_path / side / "p.json").write_text(json.dumps(p, indent=ind), newline="\n")
    ref, act = scan_tree(tmp_path / "ref"), scan_tree(tmp_path / "act")
    m = Masks().with_uuids(U)
    assert compare_trees(ref, act, tree(text_format=False), act_masks=m).ok
    assert not compare_trees(ref, act, tree(text_format=False)).ok
    assert align_json({"a": ["x_<uuid>", 1024]}, {"a": ["x_1024", 1024]}, U) == {
        "a": ["x_<uuid>", 1024]
    }
    # zarr metadata in another format: by its JSON, the uuid aligned there
    rep = Report()
    _cmp_zarr_meta('{"a": "x_<uuid>", "n": 1024}', '{\n "a": "x_1024",\n "n": 1024\n}', "m", rep, U)
    assert rep.ok and rep.notes, rep.text()


# ---------------------------------------------------------------------------------------------
# a case whose run made the uuid "1024", against a fake MATLAB reference, in both suite modes
# ---------------------------------------------------------------------------------------------

MATLAB_WORK = "/scratch/pkpy_1/w/abc"


def write_run(root, work, u, batch=1024):
    """what the case's call writes under <work>/d/DSR and returns: real values 1024 (the
    parameters, a log over 4096 bytes, a MAT-file, a zarr's shape, a file name) beside
    temporary names and texts with the uuid ``u``"""
    from scipy.io import savemat

    cppzarr = pytest.importorskip("cppzarr")
    d = Path(root) / "d" / "DSR"
    (d / "zarr_flag" / f"v_{u}").mkdir(parents=True)
    par = {"batchSize": [13, batch, 1024], "tmp": f"{work}/d/DSR/v.zarr_{u}", "uuid": u}
    (d / "parameters.json").write_text(json.dumps(par, indent=2), newline="\n")
    flag = f"result={work}/d/DSR/v.zarr\nuuid={u}\nbatch=1024\n"
    (d / "zarr_flag" / f"v_{u}" / "writer.txt").write_text(flag, newline="\n")
    log = "".join(f"task {i} of 1024: x_{u}.tif, 1024 done\n" for i in range(200))
    (d / "log.txt").write_text(log, newline="\n")
    (d / "Iter_1024_a.txt").write_text("1024\n", newline="\n")
    savemat(d / "info.mat", {"tmp": f"{work}/d/DSR/x_{u}.tif", "n": "1024"})
    vol = np.arange(2 * 1024 * 3, dtype=np.uint16).reshape(2, 1024, 3)
    cppzarr.write_zarr(
        str(d / "v.zarr"), np.asfortranarray(vol), chunks=[2, 1024, 3], order="F",
        dimension_separator=".",
    )  # fmt: skip
    (d / "v.zarr" / ".zattrs").write_text(json.dumps({"from": f"x_{u}.tif"}), newline="\n")
    return f"{work}/d/DSR/x_{u}.tif", "1024", {"file": f"x_{u}.tif", "size": "1024"}


def publish(store, case, error=None):
    """MATLAB's result of the case, as the MATLAB runner stores it (MATLAB's uuid is FULL):
    write_run's files and outputs, or the error"""
    key = case.key()
    stage = store.staging("test") / key
    stage.mkdir(parents=True)
    if error is None:
        outs = write_run(stage / "files", MATLAB_WORK, FULL)
        write_pkv(stage / "outputs.pkv", MCell((1, 3), [to_mvalue(x) for x in outs]))
    meta = {
        "format": 1, "key": key, "id": case.id, "status": "error" if error else "ok",
        "error": error, "warnings": [], "capture": error is None,
        "files": [] if error else list_tree(stage / "files"), "deleted": [],
        "work": MATLAB_WORK, "host": "node12", "oracle_commit": store.commit,
    }  # fmt: skip
    (stage / "meta.json").write_text(json.dumps(meta), newline="\n")
    store.case_dir(key).parent.mkdir(parents=True, exist_ok=True)
    os.rename(stage, store.case_dir(key))


@pytest.fixture
def harness(tmp_path, tmp_path_factory, monkeypatch):
    """a golden store; ``add(make, error)`` stores MATLAB's result of the case ``make()`` makes
    and writes its manifest"""
    commit = oracle.expected_commit()
    golden = tmp_path_factory.mktemp("g")  # short: a case's files below a key (Windows' paths)
    monkeypatch.setenv(oracle.ENV_GOLDEN, str(golden))
    monkeypatch.setenv(oracle.ENV_COMMIT, commit)
    monkeypatch.delenv(oracle.ENV_ORACLE, raising=False)
    store = oracle.Store(golden, commit)
    mod = str(tmp_path / "test_uuid_case.py")
    globs = {"__file__": mod, "ORACLE_MANIFEST": str(tmp_path / "uuid_case.oracle.json")}

    def add(make, error=None):
        def made(*a, **k):
            c = make(*a, **k)
            c._module_globals = globs
            return c

        case = made()
        publish(store, case, error)
        monkeypatch.setitem(oracle._REGISTRY, mod, [case])
        oracle._write_manifests([("n", case)], store, None)
        return made

    return types.SimpleNamespace(add=add, golden=golden, monkeypatch=monkeypatch)


def writer_case(batch=1024):
    def run(w):
        return write_run(w, w.rstrip("/"), ids.get_uuid(), batch)

    return Case("uuid1024", matlab="pk_fake_writer", args=[Work("")], python=run, nargout=3)


def force_uuid(monkeypatch, form):
    """get_uuid gives 1024abcd-... (Linux's form) or its first 4 characters (Windows')"""
    plat = "win32" if form == "win" else "linux"
    monkeypatch.setattr(ids, "sys", types.SimpleNamespace(platform=plat))
    full = "1024abcd-0000-4000-8000-00000000cafe"
    monkeypatch.setattr(ids, "uuid", types.SimpleNamespace(uuid4=lambda: full))
    return full[:4] if form == "win" else full


def run_case(h, case, mode, tmp):
    if mode == "store":
        h.monkeypatch.setenv(oracle.ENV_GOLDEN, str(h.golden))
    else:
        h.monkeypatch.delenv(oracle.ENV_GOLDEN, raising=False)
    ref = oracle.reference(case)
    assert ref.source == mode
    res = case.run_python(tmp)
    return res, oracle.compare_case(case, ref, res)


@pytest.mark.parametrize("form", ["win", "linux"])
def test_a_run_whose_uuid_equals_a_value_of_its_output_passes(harness, form, tmp_path):
    """the uuid 1024 in temporary names, flag texts, a MAT-file, a zarr's attributes and the
    outputs, beside real 1024s in parameters.json, a long log (hashed in the manifest), the
    MAT-file, a file name, the outputs and a zarr's .zarray: equal with the store and with the
    manifest"""
    make = harness.add(writer_case)
    u = force_uuid(harness.monkeypatch, form)
    case = make()
    for mode in ("store", "manifest"):
        res, rep = run_case(harness, case, mode, tmp_path / f"run_{mode}")
        assert res.uuids == (u,)
        assert rep.ok and not rep.undecided, rep.text()
        assert not [n for n in rep.notes if "format" in n], rep.notes  # equal as texts too
    oracle.check(case, tmp_path / "check")  # the pytest path, in manifest mode here


@pytest.mark.parametrize("form", ["win", "linux"])
def test_a_real_difference_beside_the_uuid_still_fails(harness, form, tmp_path):
    make = harness.add(writer_case)
    force_uuid(harness.monkeypatch, form)
    case = make(batch=2048)  # the same case for MATLAB, another value from Python
    for mode in ("store", "manifest"):
        _, rep = run_case(harness, case, mode, tmp_path / f"run_{mode}")
        text = rep.text()
        assert not rep.ok and "files/d/DSR/parameters.json" in text and "2048" in text, text
        assert len(rep.problems) == 1, text


class _NotFound(ValueError):
    identifier = "pk:notFound"


def not_found(work, u):
    return f"File '{work}/d/x_{u}.tif' not found (batch 1024)."


@pytest.mark.parametrize("form", ["win", "linux"])
def test_an_error_message_with_the_runs_uuid_compares_as_matlabs(harness, form, tmp_path):
    def error_case(u=None):
        def run(w):
            raise _NotFound(not_found(w.rstrip("/"), u or ids.get_uuid()))

        return Case(
            "uuid_error", matlab="pk_fake_reader", args=[Work("")], python=run,
            expect_error="pk:notFound", py_error=_NotFound, error_message=True,
        )  # fmt: skip

    error = {"identifier": "pk:notFound", "message": not_found(MATLAB_WORK, FULL)}
    make = harness.add(error_case, error)
    force_uuid(harness.monkeypatch, form)
    for mode in ("store", "manifest"):
        _, rep = run_case(harness, make(), mode, tmp_path / f"run_{mode}")
        assert rep.ok, rep.text()
        _, rep = run_case(harness, make("2048"), mode, tmp_path / f"bad_{mode}")
        assert [p.what for p in rep.problems] == ["message"], rep.text()
