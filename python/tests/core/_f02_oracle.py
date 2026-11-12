"""Digests and comparisons for the core oracle cases (stand-in until F04's harness lands).

A reference is one manifest entry per case: ``array`` outputs by class, MATLAB size and the
SHA-256 of their column-major bytes; ``list`` and ``text`` outputs inline when small, else by
the SHA-256 of their canonical JSON.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np

from petakit5d.core.arrays import matlab_class, msize

INLINE_LIMIT = 4000


def canonical_json(v: Any) -> str:
    return json.dumps(v, ensure_ascii=False, separators=(",", ":"), allow_nan=True)


def array_digest(cls: str, data: bytes) -> str:
    return hashlib.sha256(cls.encode() + b"|" + data).hexdigest()


def _plain(v):
    if isinstance(v, np.generic):
        return v.item()
    return v


def entry_of_python(kind: str, out: Any) -> dict:
    """The manifest entry the Python output would have."""
    if kind == "array":
        a = np.asarray(out)
        cls = matlab_class(a)
        data = (
            a.astype(np.uint8).tobytes(order="F")
            if cls == "logical"
            else np.ascontiguousarray(a.ravel(order="F")).tobytes()
        )
        return {
            "kind": "array",
            "class": cls,
            "size": list(msize(a)),
            "numel": int(a.size),
            "sha256": array_digest(cls, data),
        }
    if kind == "list":
        lst = [_plain(x) for x in out]
        s = canonical_json(lst)
        e = {"kind": "list", "n": len(lst), "sha256": hashlib.sha256(s.encode()).hexdigest()}
        if len(s) <= INLINE_LIMIT:
            e["value"] = lst
        return e
    if kind == "text":
        e = {"kind": "text", "sha256": hashlib.sha256(out.encode()).hexdigest()}
        if len(out) <= INLINE_LIMIT:
            e["value"] = out
        return e
    raise ValueError(kind)


def same(ref: dict, got: dict) -> bool:
    if ref["kind"] != got["kind"]:
        return False
    if ref["kind"] == "array":
        if ref["class"] != got["class"] or ref["numel"] != got["numel"]:
            return False
        # a 1-D NumPy result stands for either MATLAB vector orientation
        if (
            len(got["size"]) > 2
            or (got["size"][0] != 1 and got["size"][1] != 1)
            or 0 in got["size"]
        ):
            if ref["size"] != got["size"]:
                return False
        return ref["sha256"] == got["sha256"]
    return ref["sha256"] == got["sha256"]


def describe(ref: dict, got: dict) -> str:
    """A readable difference (first differing items for inline values)."""
    lines = [
        f"MATLAB: { {k: v for k, v in ref.items() if k != 'value'} }",
        f"Python: { {k: v for k, v in got.items() if k != 'value'} }",
    ]
    if "value" in ref and "value" in got:
        rv, gv = ref["value"], got["value"]
        if isinstance(rv, list):
            diffs = [(i, a, b) for i, (a, b) in enumerate(zip(rv, gv, strict=False)) if a != b]
            lines.append(
                f"{len(diffs)} differing items"
                + ("" if len(rv) == len(gv) else f", lengths {len(rv)} and {len(gv)}")
            )
            lines += [f"  [{i}] MATLAB {a!r} Python {b!r}" for i, a, b in diffs[:10]]
        else:
            lines.append(f"MATLAB {rv!r}")
            lines.append(f"Python {gv!r}")
    return "\n".join(lines)
