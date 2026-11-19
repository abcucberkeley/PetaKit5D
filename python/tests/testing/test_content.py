"""Decoded content of tiff and zarr files (petakit5d.testing.content): what the comparison of
output tiffs takes (and the case keys of file inputs will), so that the bytes cpp-tiff and cpp-zarr
versions write for the same data (strip and IFD offsets, page placement, compressed chunks)
don't count, while a voxel, a tag or a metadata field of their own does."""

import json
import math
import struct

import numpy as np
import pytest

from petakit5d.testing import cases, content
from petakit5d.testing.compare import (
    Masks,
    compare_trees,
    describe_tree,
    scan_tree,
    tree,
    tree_from_descriptor,
)

cpptiff = pytest.importorskip("cpptiff")
cppzarr = pytest.importorskip("cppzarr")


def stack(dtype=np.uint16, shape=(3, 37, 23), seed=1):
    """a (pages, y, x) stack of exact values"""
    n = math.prod(shape)
    a = (np.arange(n, dtype=np.int64) * 7919 + seed * 104729) % 60001
    return a.reshape(shape).astype(dtype)


# --- an uncompressed tiff in the byte layout asked for ------------------------------------------

_TYPES = {2: ("s", 1), 3: ("H", 2), 4: ("I", 4), 5: ("I", 8), 16: ("Q", 8)}


def _pack(e, typ, vals):
    if typ == 2:
        b = vals.encode("latin-1") + b"\0"
        return len(b), b
    code, _ = _TYPES[typ]
    flat = [x for v in vals for x in v] if typ == 5 else list(vals)
    return len(vals), struct.pack(f"{e}{len(flat)}{code}", *flat)


def _ifd(e, big, entries, at, nxt):
    """an IFD at offset ``at`` with its long values right after it, and the offset of the next"""
    fsize, esize = (8, 20) if big else (4, 12)
    head = struct.pack(e + ("Q" if big else "H"), len(entries))
    body = len(head) + len(entries) * esize + fsize
    recs, extra = b"", b""
    for tag, typ, vals in sorted(entries, key=lambda t: t[0]):
        count, b = _pack(e, typ, vals)
        if len(b) <= fsize:
            field = b.ljust(fsize, b"\0")
        else:
            field = struct.pack(e + ("Q" if big else "I"), at + body + len(extra))
            extra += b + b"\0" * (len(b) % 2)
        recs += struct.pack(e + ("HHQ" if big else "HHI"), tag, typ, count) + field
    return head + recs + struct.pack(e + ("Q" if big else "I"), nxt) + extra


def raw_tiff(
    path,
    vol,
    *,
    endian="<",
    big=True,
    strips_first=True,
    rows=None,
    short=False,
    extra=(),
    later=(),
):
    """write a (pages, y, x) stack uncompressed by hand: byte order, classic TIFF or BigTIFF,
    every strip before the IFDs (as cpp-tiff 1.7.3 places them) or each page's after its IFD
    (as 1.8.0 does), rows per strip, sizes as SHORT or LONG; ``extra`` (tag, type, values)
    entries go into page 0's IFD, ``later`` ones into page 1's"""
    vol = np.asarray(vol)
    n, h, w = vol.shape
    rows = rows or h
    starts = list(range(0, h, rows))
    counts = [(min(r + rows, h) - r) * w * vol.dtype.itemsize for r in starts]
    page = h * w * vol.dtype.itemsize
    fmt = {"u": 1, "i": 2, "f": 3}[vol.dtype.kind]
    otype = 16 if big else 4
    dim = 3 if short else 4

    def entries(k, offs):
        es = [
            (256, dim, [w]),
            (257, dim, [h]),
            (258, 3, [vol.dtype.itemsize * 8]),
            (259, 3, [1]),
            (262, 3, [1]),
            (273, otype, offs),
            (277, 3, [1]),
            (278, 4, [rows]),
            (279, otype, counts),
            (284, 3, [1]),
            (339, 3, [fmt]),
        ]
        return es + list(extra if k == 0 else later if k == 1 else ())

    e = endian
    hlen = 16 if big else 8
    sizes = [len(_ifd(e, big, entries(k, [0] * len(starts)), 0, 0)) for k in range(n)]
    pages = [vol[k].astype(vol.dtype.newbyteorder(e)).tobytes() for k in range(n)]
    if strips_first:
        data_at = [hlen + k * page for k in range(n)]
        ifd_at = [hlen + n * page]
        for k in range(n - 1):
            ifd_at.append(ifd_at[-1] + sizes[k])
    else:
        ifd_at, data_at, at = [], [], hlen
        for k in range(n):
            ifd_at.append(at)
            data_at.append(at + sizes[k])
            at += sizes[k] + page
    blocks = {}
    for k in range(n):
        offs = [data_at[k] + sum(counts[:i]) for i in range(len(starts))]
        nxt = ifd_at[k + 1] if k + 1 < n else 0
        blocks[ifd_at[k]] = _ifd(e, big, entries(k, offs), ifd_at[k], nxt)
        blocks[data_at[k]] = pages[k]
    order = b"II" if e == "<" else b"MM"
    if big:
        head = order + struct.pack(e + "HHHQ", 43, 8, 0, ifd_at[0])
    else:
        head = order + struct.pack(e + "HI", 42, ifd_at[0])
    out = bytearray(head)
    for at in sorted(blocks):
        assert len(out) == at, (len(out), at)
        out += blocks[at]
    path.write_bytes(bytes(out))
    return path


LAYOUTS = {
    "little_big_strips_first_8_rows": {"rows": 8},
    "little_big_strips_after_each_ifd": {"strips_first": False},
    "big_endian_classic": {"endian": ">", "big": False},
    "classic_short_sizes_one_strip": {"big": False, "short": True},
}


# --- tiff content -----------------------------------------------------------------------------


@pytest.mark.parametrize("layout", list(LAYOUTS))
@pytest.mark.parametrize("dtype", [np.uint16, np.float32])
def test_a_tiff_of_another_byte_layout_has_the_same_content(tmp_path, layout, dtype):
    vol = stack(dtype)
    a = raw_tiff(tmp_path / "a.tif", vol)
    b = raw_tiff(tmp_path / "b.tif", vol, **LAYOUTS[layout])
    assert a.read_bytes() != b.read_bytes()
    assert np.array_equal(cpptiff.read_tiff(str(b)), vol)
    ca, cb = content.tiff_content(a), content.tiff_content(b)
    assert ca == cb and ca["shape"] == [3, 37, 23] and ca["dtype"] == np.dtype(dtype).name
    assert content.tiff_digest(a) == content.tiff_digest(b)


def test_the_layout_tags_are_left_out(tmp_path):
    pages = content.tiff_pages(raw_tiff(tmp_path / "a.tif", stack(), rows=8))
    assert len(pages) == 3
    assert sorted(map(int, pages[0])) == [256, 257, 258, 259, 262, 277, 284, 339]
    assert pages[0]["256"] == ["int", [23]] and pages[0]["339"] == ["int", [1]]


def test_the_tags_that_name_the_writer_are_left_out(tmp_path):
    # MATLAB's blockedImage writes classic TIFF with Software 'MATLAB blockedImage', the port
    # cpp-tiff's BigTIFF without: the same content
    vol = stack()
    writer = [(305, 2, "MATLAB blockedImage"), (306, 2, "2026:11:19 09:30:00"), (316, 2, "node12")]
    a = raw_tiff(tmp_path / "a.tif", vol, big=False, extra=writer)
    b = raw_tiff(tmp_path / "b.tif", vol)
    assert content.tiff_content(a) == content.tiff_content(b)
    assert "305" not in content.tiff_pages(a)[0]


def test_cpp_tiffs_files_have_the_content_of_the_same_data_written_by_hand(tmp_path):
    vol = stack()
    a = tmp_path / "lzw.tif"
    cpptiff.write_tiff(str(a), vol)
    b = raw_tiff(tmp_path / "raw.tif", vol)
    ca, cb = content.tiff_content(a), content.tiff_content(b)
    assert ca["voxels"] == cb["voxels"] and ca["shape"] == cb["shape"]
    # (cpp-tiff writes LZW: the compression is a tag of the file's content)
    assert ca["tags"]["p0"]["259"] == ["int", [5]] and cb["tags"]["p0"]["259"] == ["int", [1]]
    assert ca["tags"]["sha"] != cb["tags"]["sha"]


@pytest.mark.parametrize(
    "change",
    ["voxel", "description", "resolution", "one_page_less", "dtype"],
)
def test_a_voxel_a_tag_or_a_page_changes_the_content(tmp_path, change):
    vol = stack()
    desc = [(270, 2, "ImageJ=1.54f\nimages=3\n")]
    a = raw_tiff(tmp_path / "a.tif", vol, extra=desc)
    other = vol.copy()
    extra = desc
    if change == "voxel":
        other[2, 36, 22] += 1
    elif change == "description":
        extra = [(270, 2, "ImageJ=1.54f\nimages=4\n")]
    elif change == "resolution":
        extra = [*desc, (282, 5, [(1, 3)])]
    elif change == "one_page_less":
        other = vol[:2]
    else:
        other = vol.astype(np.int16)
    b = raw_tiff(tmp_path / "b.tif", other, strips_first=False, extra=extra)
    assert content.tiff_digest(a) != content.tiff_digest(b)


@pytest.mark.parametrize("what", ["empty", "text", "truncated", "loop"])
def test_files_that_are_not_tiffs_have_no_content(tmp_path, what):
    p = tmp_path / "x.tif"
    good = raw_tiff(tmp_path / "good.tif", stack()).read_bytes()
    p.write_bytes(
        {
            "empty": b"",
            "text": b"not a tiff\n",
            "truncated": good[:200],
            # the first IFD's next-IFD offset pointing back at itself
            "loop": good,
        }[what]
    )
    if what == "loop":
        pages = content.tiff_pages(p)
        assert len(pages) == 3
        first = struct.unpack("<Q", good[8:16])[0]
        n = struct.unpack("<Q", good[first : first + 8])[0]
        at = first + 8 + n * 20
        p.write_bytes(good[:at] + struct.pack("<Q", first) + good[at + 8 :])
    assert content.tiff_content(p) is None and content.tiff_digest(p) is None
    with pytest.raises(content.NotTiff):
        content.tiff_pages(p)


# --- zarr content -----------------------------------------------------------------------------


def zarr(path, vol, **kw):
    opts = {"chunks": [2, 16, 16], "order": "C", "dimension_separator": ".", **kw}
    cppzarr.write_zarr(str(path), vol, **opts)
    return path


def test_a_zarr_by_its_voxels_and_metadata_not_its_chunk_bytes(tmp_path):
    vol = stack(shape=(5, 37, 23))
    ref = content.zarr_digest(zarr(tmp_path / "a.zarr", vol))
    assert ref is not None
    same = {
        "lz4": {"cname": "lz4", "clevel": 5},
        "slash": {"dimension_separator": "/"},
    }
    for name, kw in same.items():
        p = zarr(tmp_path / f"{name}.zarr", vol, **kw)
        assert cases.tree_hash(p) != cases.tree_hash(tmp_path / "a.zarr"), name
        assert content.zarr_digest(p) == ref, name
    other = vol.copy()
    other[4, 36, 22] += 1
    differ = {
        "voxel": (other, {}),
        "chunks": (vol, {"chunks": [3, 16, 16]}),
        "dtype": (vol.astype(np.int32), {}),
    }
    for name, (v, kw) in differ.items():
        assert content.zarr_digest(zarr(tmp_path / f"d_{name}.zarr", v, **kw)) != ref, name
    # the order flag is metadata: MATLAB's (y, x, z) order F against a native order C array
    f = np.asfortranarray(np.transpose(vol, (1, 2, 0)))
    p = zarr(tmp_path / "f.zarr", f, chunks=[16, 16, 2], order="F")
    assert content.zarr_array_meta(p)["order"] == "F" and content.zarr_digest(p) != ref


# --- output tiffs: compared by content --------------------------------------------------------


def _trees(tmp_path, vol_a, vol_b, kw_a=None, kw_b=None):
    for side, v, kw in (("a", vol_a, kw_a or {}), ("b", vol_b, kw_b or {})):
        (tmp_path / side).mkdir()
        raw_tiff(tmp_path / side / "im.tif", v, **kw)
    return scan_tree(tmp_path / "a"), scan_tree(tmp_path / "b")


def test_output_tiffs_of_another_layout_are_equal_and_in_a_manifest_too(tmp_path):
    vol = stack()
    ref, act = _trees(tmp_path, vol, vol, None, {"strips_first": False, "rows": 5})
    rep = compare_trees(ref, act)
    assert rep.ok, rep.text()
    assert any("layout only" in n for n in rep.notes), rep.notes
    rep = compare_trees(ref, act, tree(tiff_bytes=True))
    assert not rep.ok and "files/im.tif: layout: voxels and tags are equal" in rep.text()
    d = json.loads(json.dumps(describe_tree(ref, tree(), Masks())))
    desc = d["files"]["im.tif"][0]
    assert desc["tags"]["n"] == 3 and desc["tags"]["p0"]["259"] == ["int", [1]]
    rep = compare_trees(tree_from_descriptor(d), act)
    assert rep.ok and not rep.undecided, rep.text()
    # a manifest written before tags were kept cannot decide them: skipped, not passed
    del desc["tags"]
    rep = compare_trees(tree_from_descriptor(d), act)
    assert rep.ok and rep.undecided and "keeps no tags" in rep.undecided[0]


def test_a_voxel_or_a_tag_difference_is_still_found(tmp_path):
    vol = stack()
    other = vol.copy()
    other[1, 2, 3] = 7
    ref, act = _trees(tmp_path, vol, other, None, {"strips_first": False})
    text = compare_trees(ref, act).text()
    assert "files/im.tif (voxels): values: 1 of 2553 values differ" in text, text
    for kw, line in (
        ({"extra": [(270, 2, "x")]}, "page 0 ImageDescription (270): MATLAB absent"),
        ({"extra": [(282, 5, [(72, 1)])]}, "page 0 XResolution (282): MATLAB absent"),
    ):
        sub = tmp_path / line[:10].replace(" ", "_")
        sub.mkdir()
        ref, act = _trees(sub, vol, vol, None, kw)
        rep = compare_trees(ref, act)
        assert not rep.ok and f"files/im.tif: tags: {line}" in rep.text(), rep.text()
        # and with the reference from a manifest
        d = json.loads(json.dumps(describe_tree(ref, tree(), Masks())))
        assert line in compare_trees(tree_from_descriptor(d), act).text()
    # a later page's tag: named with both files, found by the hash of all pages with a manifest
    # (which keeps the first page's tags)
    sub = tmp_path / "later"
    sub.mkdir()
    ref, act = _trees(sub, vol, vol, None, {"later": [(297, 3, [1, 3])]})
    text = compare_trees(ref, act).text()
    assert "files/im.tif: tags: page 1 PageNumber (297): MATLAB absent, Python [1,3]" in text
    d = json.loads(json.dumps(describe_tree(ref, tree(), Masks())))
    text = compare_trees(tree_from_descriptor(d), act).text()
    assert "files/im.tif: tags: the tags of a later page differ" in text, text


def test_lzw_against_uncompressed_names_the_compression(tmp_path):
    vol = stack()
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    cpptiff.write_tiff(str(tmp_path / "a" / "im.tif"), vol)
    raw_tiff(tmp_path / "b" / "im.tif", vol)
    rep = compare_trees(scan_tree(tmp_path / "a"), scan_tree(tmp_path / "b"))
    assert "files/im.tif: tags: page 0 Compression (259): MATLAB [5], Python [1]" in rep.text()
    assert "voxels" not in rep.text()
