"""Line ends that differ by system: cpp-zarr writes zarr metadata in text mode (CRLF on
Windows), which MATLAB reads the same, so oracle keys and comparisons take it with LF."""

from petakit5d.testing.cases import tree_hash


def _tree(root, zarray: bytes, other: bytes):
    (root / "im.zarr").mkdir(parents=True)
    (root / "im.zarr" / ".zarray").write_bytes(zarray)
    (root / "list.csv").write_bytes(other)
    return root


def test_zarr_metadata_hashes_alike_with_crlf_and_lf(tmp_path):
    lf = _tree(tmp_path / "lf", b'{\n    "chunks": [1]\n}\n', b"a,b\n1,2\n")
    crlf = _tree(tmp_path / "crlf", b'{\r\n    "chunks": [1]\r\n}\r\n', b"a,b\n1,2\n")
    for content in (True, False):
        assert tree_hash(lf, content) == tree_hash(crlf, content)


def test_other_files_keep_their_bytes(tmp_path):
    lf = _tree(tmp_path / "lf", b"{}\n", b"a,b\n1,2\n")
    crlf = _tree(tmp_path / "crlf", b"{}\n", b"a,b\r\n1,2\r\n")
    assert tree_hash(lf) != tree_hash(crlf)
