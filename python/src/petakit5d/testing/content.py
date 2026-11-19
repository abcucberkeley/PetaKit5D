"""Decoded content of tiff and zarr files, without the bytes a library happens to write.

A tiff's content is what a reader gets from it: its voxels (data type, the pages' shape and
values) and its tags, without the byte layout (strip and tile offsets and byte counts, rows per
strip, tile sizes, where the IFDs and the strips sit in the file, byte order, classic TIFF or
BigTIFF) and the tags that name the writer (software, computer, time). A zarr array's content
is its decoded voxels and the metadata that decides them (format, shape, chunks, data type,
order, fill value; and its attributes), not the compressed chunk bytes. Two libraries that
write the same data with other bytes give files of the same content; a voxel or a tag of their
own differs.

Used by the case keys of file inputs (``cases.File``, ``cases.Tree``) and by the comparison of
the tiff files a call writes (``compare.compare_trees``, the manifests' tiff descriptors).
"""

from __future__ import annotations

import hashlib
import json
import os
import struct
from pathlib import Path
from typing import Any

import numpy as np

# tags that say where the bytes are, not what they hold
LAYOUT_TAGS = {
    273: "StripOffsets",
    278: "RowsPerStrip",
    279: "StripByteCounts",
    288: "FreeOffsets",
    289: "FreeByteCounts",
    322: "TileWidth",
    323: "TileLength",
    324: "TileOffsets",
    325: "TileByteCounts",
    330: "SubIFDs",
    513: "JPEGInterchangeFormat",
    514: "JPEGInterchangeFormatLength",
    34665: "ExifIFD",
    34853: "GPSIFD",
    40965: "InteroperabilityIFD",
}
# tags that name the program, the computer and the time that wrote the file, not what it holds
# (MATLAB's blockedImage writes Software 'MATLAB blockedImage', cpp-tiff none)
WRITER_TAGS = {305: "Software", 306: "DateTime", 316: "HostComputer"}
TAG_NAMES = {
    254: "NewSubfileType",
    255: "SubfileType",
    256: "ImageWidth",
    257: "ImageLength",
    258: "BitsPerSample",
    259: "Compression",
    262: "PhotometricInterpretation",
    266: "FillOrder",
    269: "DocumentName",
    270: "ImageDescription",
    271: "Make",
    272: "Model",
    274: "Orientation",
    277: "SamplesPerPixel",
    280: "MinSampleValue",
    281: "MaxSampleValue",
    282: "XResolution",
    283: "YResolution",
    284: "PlanarConfiguration",
    285: "PageName",
    296: "ResolutionUnit",
    297: "PageNumber",
    315: "Artist",
    317: "Predictor",
    320: "ColorMap",
    338: "ExtraSamples",
    339: "SampleFormat",
    340: "SMinSampleValue",
    341: "SMaxSampleValue",
    700: "XMP",
    34675: "ICCProfile",
    50838: "IJMetadataByteCounts",
    50839: "IJMetadata",
    **LAYOUT_TAGS,
    **WRITER_TAGS,
}
TIFF_SUFFIXES = (".tif", ".tiff")
ZARR_META = (".zarray", ".zattrs", ".zgroup", "zarr.json")
MAX_PAGES = 1 << 22  # IFDs followed before a file counts as broken (a loop)
INLINE_TAG = 256  # tag values up to this JSON size are kept as they are in descriptors

# TIFF field types: struct code, size in bytes, kind
_TYPES = {
    1: ("B", 1, "bytes"),
    2: ("s", 1, "ascii"),
    3: ("H", 2, "int"),
    4: ("I", 4, "int"),
    5: ("I", 8, "rational"),
    6: ("b", 1, "int"),
    7: ("B", 1, "bytes"),
    8: ("h", 2, "int"),
    9: ("i", 4, "int"),
    10: ("i", 8, "rational"),
    11: ("f", 4, "float"),
    12: ("d", 8, "float"),
    13: ("I", 4, "int"),
    16: ("Q", 8, "int"),
    17: ("q", 8, "int"),
    18: ("Q", 8, "int"),
}


class NotTiff(ValueError):
    """the file is not a TIFF this module can walk (no header, a broken or looping IFD chain)"""


def is_tiff_name(name: str | os.PathLike) -> bool:
    return os.fspath(name).lower().endswith(TIFF_SUFFIXES)


def tag_name(code: int | str) -> str:
    code = int(code)
    return f"{TAG_NAMES.get(code, 'tag')} ({code})"


def _canon(x: Any) -> str:
    return json.dumps(x, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha(*parts: str | bytes) -> str:
    h = hashlib.sha256()
    for p in parts:
        h.update(p.encode() if isinstance(p, str) else p)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------
# tiff
# ---------------------------------------------------------------------------------------------


def _value(f, size: int, e: str, big: bool, typ: int, cnt: int, field: bytes) -> list:
    """an entry's value as [kind, value]: integers as a list, ASCII as text (trailing NULs
    dropped), BYTE and UNDEFINED as hex, rationals as [numerator, denominator] pairs, floats
    as their exact hex form; the field type itself is not kept (SHORT and LONG are alike)"""
    if typ not in _TYPES:
        return [f"type{typ}", field.hex()]
    code, n, kind = _TYPES[typ]
    nbytes = cnt * n
    if nbytes <= len(field):
        data = field[:nbytes]
    else:
        off = struct.unpack(e + ("Q" if big else "I"), field)[0]
        if off + nbytes > size:
            raise NotTiff(f"a tag value past the end of the file ({off} + {nbytes} > {size})")
        f.seek(off)
        data = f.read(nbytes)
    if kind == "bytes":
        return [kind, data.hex()]
    if kind == "ascii":
        return [kind, data.rstrip(b"\0").decode("latin-1")]
    vals = struct.unpack(f"{e}{cnt * (2 if kind == 'rational' else 1)}{code}", data)
    if kind == "rational":
        return [kind, [list(vals[i : i + 2]) for i in range(0, len(vals), 2)]]
    if kind == "float":
        return [kind, [float(v).hex() for v in vals]]
    return [kind, list(vals)]


def tiff_pages(path: str | os.PathLike) -> list[dict[str, list]]:
    """the tags of each page (IFD) along the file's IFD chain, without the layout and writer
    tags: ``{str(tag): [kind, value]}`` per page (see ``_value``)"""
    path = Path(path)
    size = path.stat().st_size
    with open(path, "rb") as f:
        head = f.read(16)
        if head[:2] == b"II":
            e = "<"
        elif head[:2] == b"MM":
            e = ">"
        else:
            raise NotTiff("no TIFF header")
        if len(head) < 8:
            raise NotTiff("no TIFF header")
        magic = struct.unpack(e + "H", head[2:4])[0]
        if magic == 42:
            big, off = False, struct.unpack(e + "I", head[4:8])[0]
        elif magic == 43 and len(head) == 16 and struct.unpack(e + "H", head[4:6])[0] == 8:
            big, off = True, struct.unpack(e + "Q", head[8:16])[0]
        else:
            raise NotTiff("no TIFF header")
        pages: list[dict[str, list]] = []
        seen: set[int] = set()
        nsize, esize, fsize = (8, 20, 8) if big else (2, 12, 4)
        while off:
            if off in seen or off + nsize > size or len(pages) >= MAX_PAGES:
                raise NotTiff(f"a broken IFD chain (IFD {len(pages)} at {off})")
            seen.add(off)
            f.seek(off)
            n = struct.unpack(e + ("Q" if big else "H"), f.read(nsize))[0]
            raw = f.read(n * esize + fsize)
            if len(raw) < n * esize + fsize:
                raise NotTiff(f"IFD {len(pages)} ends past the end of the file")
            tags: dict[str, list] = {}
            for i in range(n):
                p = i * esize
                tag, typ, cnt = struct.unpack_from(e + ("HHQ" if big else "HHI"), raw, p)
                if tag in LAYOUT_TAGS or tag in WRITER_TAGS:
                    continue
                field = raw[p + esize - fsize : p + esize]
                here = f.tell()
                tags[str(tag)] = _value(f, size, e, big, typ, cnt, field)
                f.seek(here)
            pages.append(tags)
            off = struct.unpack(e + ("Q" if big else "I"), raw[n * esize :])[0]
    if not pages:
        raise NotTiff("no IFD")
    return pages


def _short(value: list) -> list:
    """a tag value for a descriptor: as it is, or its hash when it is long"""
    t = _canon(value[1])
    return value if len(t) <= INLINE_TAG else [value[0], "sha256:" + _sha(t)]


def tags_summary(pages: list[dict[str, list]]) -> dict:
    """the tags of a file for a descriptor: the number of pages, the hash of every page's tags
    and the first page's tags (long values by their hash)"""
    return {
        "n": len(pages),
        "sha": _sha("pk2:tifftags:", _canon(pages)),
        "p0": {k: _short(v) for k, v in pages[0].items()},
    }


def tags_differences(r: dict, a: dict, rpages: list | None = None, apages: list | None = None):
    """lines naming how two files' tags differ (MATLAB's first); [] when they are the same.
    With both files' pages the first page that differs is named, else the first page only"""
    if r["sha"] == a["sha"]:
        return []
    out = []
    if r["n"] != a["n"]:
        out.append(f"pages (IFDs): MATLAB {r['n']}, Python {a['n']}")
    if rpages is not None and apages is not None:
        pairs = list(enumerate(zip(rpages, apages, strict=False)))
    else:
        pairs = [(0, (r["p0"], a["p0"]))]
    for i, (rp, ap) in pairs:
        rows = []
        for k in sorted(set(rp) | set(ap), key=int):
            if rp.get(k) != ap.get(k):
                rv = "absent" if k not in rp else _canon(_short(rp[k])[1])[:120]
                av = "absent" if k not in ap else _canon(_short(ap[k])[1])[:120]
                rows.append(f"page {i} {tag_name(k)}: MATLAB {rv}, Python {av}")
        if rows:
            out.extend(rows[:8])
            break
    if not out:
        out.append(
            "the tags of a later page differ (the manifest keeps the first page's; the store "
            "names the page)"
        )
    return out


def tiff_voxels(path: str | os.PathLike, pages: list[dict[str, list]] | None = None) -> np.ndarray:
    """the decoded voxels as (pages, y, x) (with samples last when there are several), whatever
    shape the reader gives them (cpp-tiff versions shape an ImageJ stack differently)"""
    import cpptiff

    a = np.asarray(cpptiff.read_tiff(os.fspath(path), num_threads=1))
    p0 = (pages or tiff_pages(path))[0]
    h = p0.get("257", ["int", [0]])[1][0]
    w = p0.get("256", ["int", [0]])[1][0]
    s = p0.get("277", ["int", [1]])[1][0]
    per = h * w * s
    if per and a.size % per == 0 and a.size:
        a = a.reshape((a.size // per, h, w) + ((s,) if s > 1 else ()))
    return a


def tiff_content(path: str | os.PathLike) -> dict | None:
    """a tiff's content: data type, shape (pages, y, x[, samples]), the hash of its voxels and
    its tags (``tags_summary``); None when the file is not a tiff that decodes (an empty flag
    file, a broken file)"""
    try:
        pages = tiff_pages(path)
        a = tiff_voxels(path, pages)
    except ImportError:
        raise
    except Exception:  # not a tiff the readers take: the caller uses the bytes
        return None
    le = np.ascontiguousarray(a, dtype=a.dtype.newbyteorder("<"))
    h = hashlib.sha256(f"pk2:tiffvox:{a.dtype.name}:{list(a.shape)}:".encode())
    h.update(memoryview(le).cast("B"))
    return {
        "dtype": a.dtype.name,
        "shape": list(a.shape),
        "voxels": h.hexdigest(),
        "tags": tags_summary(pages),
    }


def tiff_digest(path: str | os.PathLike) -> str | None:
    """SHA-256 of a tiff's content (``tiff_content``), or None"""
    c = tiff_content(path)
    return None if c is None else _sha("pk2:tiff:", _canon(c))


# ---------------------------------------------------------------------------------------------
# zarr
# ---------------------------------------------------------------------------------------------


def is_zarr_folder(path: str | os.PathLike) -> bool:
    """a folder named *.zarr with zarr metadata in it"""
    p = Path(path)
    return p.name.lower().endswith(".zarr") and any((p / m).is_file() for m in ZARR_META)


def _json(p: Path) -> Any:
    return json.loads(p.read_bytes().decode("utf-8"))


def _attrs(p: Path) -> Any:
    """user attributes, with an empty set the same as none"""
    a = _json(p) if p.is_file() else None
    return a or None


def zarr_array_meta(path: str | os.PathLike) -> dict | None:
    """the metadata that decides a zarr array's voxels, or None when the folder is no array:
    v2's zarr_format, shape, chunks, dtype, order and fill_value (and .zattrs); v3's shape,
    data type, chunk grid, fill value, the order of a transpose codec and the inner chunks of
    a sharding codec (and the attributes). Compressors, filters, chunk-key encodings and
    dimension separators change only the bytes."""
    p = Path(path)
    if (p / ".zarray").is_file():
        m = _json(p / ".zarray")
        out = {k: m.get(k) for k in ("zarr_format", "shape", "chunks", "dtype", "order")}
        out["fill_value"] = m.get("fill_value")
        out["attributes"] = _attrs(p / ".zattrs")
        return out
    if (p / "zarr.json").is_file():
        m = _json(p / "zarr.json")
        if m.get("node_type") != "array":
            return None
        out = {k: m.get(k) for k in ("zarr_format", "shape", "data_type", "chunk_grid")}
        out["fill_value"] = m.get("fill_value")
        out["attributes"] = m.get("attributes") or None

        def walk(codecs: Any) -> None:
            for c in codecs or []:
                if not isinstance(c, dict):
                    continue
                cfg = c.get("configuration") or {}
                if c.get("name") == "transpose":
                    out.setdefault("transpose", []).append(cfg.get("order"))
                elif c.get("name") == "sharding_indexed":
                    out.setdefault("inner_chunks", []).append(cfg.get("chunk_shape"))
                    walk(cfg.get("codecs"))

        walk(m.get("codecs"))
        return out
    return None


def _zarr_array_digest(p: Path, meta: dict) -> str:
    import cppzarr

    a = np.asarray(cppzarr.read_zarr(os.fspath(p), num_threads=1))
    le = np.ascontiguousarray(a, dtype=a.dtype.newbyteorder("<"))
    h = hashlib.sha256(f"pk2:zarrvox:{a.dtype.name}:{list(a.shape)}:".encode())
    h.update(memoryview(le).cast("B"))
    return _sha("pk2:zarr:", _canon(meta), h.hexdigest())


def zarr_digest(path: str | os.PathLike) -> str | None:
    """SHA-256 of a zarr folder's content: an array by its metadata (``zarr_array_meta``) and
    decoded voxels; a group by its metadata and its members in order (arrays by their content,
    metadata files by their JSON, other files by their bytes). None when an array does not
    decode (the caller uses the bytes)."""
    p = Path(path)
    try:
        meta = zarr_array_meta(p)
        if meta is not None:
            return _zarr_array_digest(p, meta)
        lines = []
        for dp, dns, fns in os.walk(p):
            d = Path(dp)
            rel = d.relative_to(p).as_posix()
            if d != p and zarr_array_meta(d) is not None:
                lines.append(f"{rel}/\0{_zarr_array_digest(d, zarr_array_meta(d))}")
                dns[:] = []
                continue
            for fn in fns:
                q = d / fn
                r = q.relative_to(p).as_posix()
                if fn in ZARR_META:
                    lines.append(f"{r}\0json:{_canon(_json(q))}")
                else:
                    lines.append(f"{r}\0{_sha(q.read_bytes())}")
        return _sha("pk2:zarrgroup:", "\n".join(sorted(lines)))
    except ImportError:
        raise
    except Exception:  # not a zarr the readers take: the caller uses the bytes
        return None
