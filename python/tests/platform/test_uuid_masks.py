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
