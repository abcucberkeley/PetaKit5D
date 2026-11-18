"""get_uuid's Windows form (its first 4 characters, as MATLAB there) is masked in oracle
comparisons only where the run made it: 4 hex characters elsewhere stay names."""

import sys

from petakit5d.core import ids
from petakit5d.testing.cases import UUID_LOG_ENV
from petakit5d.testing.compare import Masks, uuid_masks


def test_only_the_runs_own_uuids_are_masked():
    m = Masks().plus(*uuid_masks(["6cac"]))
    assert m.text("zarr_flag/crop_6cac/writer.txt") == "zarr_flag/crop_<uuid>/writer.txt"
    assert m.text("im_6cac") == "im_<uuid>"
    for name in ("/mnt/projects/embryo_2024/cell1", "tile_0001.tif", "x_07d5", "ab6cacd"):
        assert m.text(name) == name


def test_full_uuids_are_left_to_the_usual_rule():
    full = "0f8fad5b-d9cb-469f-a165-70867728950e"
    assert uuid_masks([full]) == ()
    assert Masks().text(f"crop_{full}") == "crop_<uuid>"


def test_get_uuid_logs_what_it_gives(tmp_path, monkeypatch):
    log = tmp_path / "uuids.log"
    monkeypatch.setenv(UUID_LOG_ENV, str(log))
    monkeypatch.setattr(sys, "platform", "win32")  # the Windows form, also here
    a, b = ids.get_uuid(), ids.get_uuid()
    monkeypatch.undo()
    assert len(a) == 4 and len(b) == 4
    assert log.read_text().split() == [a, b]


def test_a_value_equal_to_the_runs_short_uuid_is_masked_on_both_sides(tmp_path):
    # a Windows run whose uuid was "1024", against a reference value 1024 (in parameters.json)
    from petakit5d.testing import oracle
    from petakit5d.testing.cases import PyResult

    ref = oracle.Reference(
        key="k", status="ok", error=None, warnings=[], outputs=[], files=None, meta={},
        source="manifest", masks=Masks(), commit="",
    )  # fmt: skip
    res = PyResult(
        outputs=[], error=None, error_text="", warnings=[], files=None, work=tmp_path,
        elapsed=0.0, uuids=("1024",),
    )  # fmt: skip
    rm, am = oracle.case_masks(ref, res)
    assert rm.text("    1024,") == am.text("    1024,") == "    <uuid>,"


def test_a_name_with_a_token_equal_to_the_runs_short_uuid_matches(tmp_path):
    # a Windows run whose uuid was "0000": its names are masked when scanned
    # (Scan_Iter_<uuid>_a.txt), so the reference's Scan_Iter_0000_a.txt is masked the same way
    from petakit5d.testing.compare import compare_trees, scan_tree

    for side in ("ref", "act"):
        (tmp_path / side).mkdir()
        (tmp_path / side / "Scan_Iter_0000_a.txt").write_text("x\n", newline="\n")
    um = Masks().plus(*uuid_masks(["0000"]))
    ref = scan_tree(tmp_path / "ref")
    act = scan_tree(tmp_path / "act", masks=um)
    rep = compare_trees(ref, act, ref_masks=um, act_masks=um)
    assert rep.ok, rep.text()
