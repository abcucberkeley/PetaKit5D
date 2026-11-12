import pytest

from petakit5d.core import textfiles


def test_read_text_file_lines(tmp_path):
    p = tmp_path / "a.txt"
    p.write_bytes(b"a\r\nb\n\nc\rd")
    assert textfiles.read_text_file(p) == ["a", "b", "", "c", "d"]
    p.write_bytes(b"x\n")
    assert textfiles.read_text_file(p) == ["x"]
    p.write_bytes(bytes([233, 10, 120]))
    assert textfiles.read_text_file(p) == ["\ufffd", "x"]


def test_read_text_file_empty_is_minus_one(tmp_path):
    # readTextFile of an empty file: fgetl's -1, kept as MATLAB returns it
    p = tmp_path / "e.txt"
    p.write_bytes(b"")
    assert textfiles.read_text_file(p) == [-1]


def test_read_text_file_missing(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        textfiles.read_text_file(tmp_path / "nope.txt")


def test_write_text_file(tmp_path):
    p = tmp_path / "w.txt"
    textfiles.write_text_file(["line1", "l2 \u00e9", "", "x\\ny %d"], p)
    assert p.read_bytes() == "line1\nl2 \u00e9\n\nx\\ny %d\n".encode()
    textfiles.write_text_file("raw", p)
    assert p.read_bytes() == b"raw"
    textfiles.write_text_file([], p)
    assert p.read_bytes() == b"\n"
    textfiles.write_text_file([f"n{i}" for i in range(25)], p, 10)
    assert p.read_bytes().count(b"\n") == 25
    with pytest.raises(TypeError):
        textfiles.write_text_file(5, p)


def test_write_json_file_goes_through_fprintf(tmp_path):
    # writeJsonFile: fprintf(fid, s) takes the JSON as a format
    p = tmp_path / "p.json"
    textfiles.write_json_file({"a": "plain", "b": [1, 2], "c": 0.108}, p)
    assert p.read_text() == '{\n  "a": "plain",\n  "b": [\n    1,\n    2\n  ],\n  "c": 0.108\n}'
    textfiles.write_json_file({"a": "C:\\data\\new", "b": "x%03d.tif"}, p)
    assert p.read_bytes() == b'{\n  "a": "C:\\data\\new",\n  "b": "x'  # truncated at %03d
    textfiles.write_json_file({"p": "100%% done"}, p)
    assert p.read_text() == '{\n  "p": "100% done"\n}'
    textfiles.write_json_file({"ctl": "\x01x"}, p)
    assert p.read_text() == '{\n  "ctl": "'  # \u0001 is not an escape sprintf knows
