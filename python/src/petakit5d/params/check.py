"""check_registry: every registered spec against the inputParser blocks of the MATLAB tree.

For each registered function: the MATLAB function exists; parameter names, kinds and order
agree; validator texts agree; defaults agree (by value, or by text for defaults that are not
literals); the inputParser settings agree; and where the function has an mcc parser, its
string conversions and defaults agree with the spec (the command line follows them). The
MATLAB side comes from a tree (the oracle tree on the cluster) or from a stored extract (CI).
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from petakit5d.params.literal import LiteralError, parse
from petakit5d.params.matlab_source import MatlabFunction, MatlabParam, extract_tree
from petakit5d.params.spec import _MISSING, FunctionSpec, ParamSpec, RegistryError
from petakit5d.params.values import to_python

EXTRACT_VERSION = 1


SITE_PATTERN_ENV = "PETAKIT5D_SITE_PATTERN"

# MATLAB texts that hold site values: absolute paths (POSIX with two or more components,
# Windows drives, home-relative) and Slurm options
_GENERIC_SITE = (
    r"(?<![\w.~])/[\w.~+-]+/[\w.~+-]+"
    r"|\b[A-Za-z]:[\\/]"
    r"|(?<!\w)~[\\/]"
    r"|--(?:qos|partition|account|constraint)\b"
    r"|(?-i:(?<![\w-])-[pA]\s+\w)"
)


def site_pattern() -> re.Pattern[str]:
    """What marks a MATLAB text as site-specific (it is then stored only as a hash).

    Generic rules only; a local scan can add a regular expression of its own through the
    environment variable PETAKIT5D_SITE_PATTERN, set outside the repository.
    """
    import os

    extra = os.environ.get(SITE_PATTERN_ENV, "")
    return re.compile(_GENERIC_SITE + (f"|(?:{extra})" if extra else ""), re.I)


_SITE_RE = re.compile(r"^<site value sha256:([0-9a-f]{16})>$")


def site_placeholder(text: str) -> str:
    """How a site-specific MATLAB text is stored in the repository: its short SHA-256
    (blanks outside quotes ignored), so a change on the MATLAB side is still seen."""
    return f"<site value sha256:{hashlib.sha256(_norm(text).encode()).hexdigest()[:16]}>"


def is_site_value(text: Any) -> bool:
    """Whether a stored text is the placeholder of a site-specific MATLAB text."""
    return isinstance(text, str) and bool(_SITE_RE.match(text))


def _norm(text: str | None) -> str:
    """MATLAB text without blanks outside quotes (how two validator texts are compared)."""
    if text is None:
        return ""
    parts = re.split(r"('(?:[^']|'')*')", text.strip())
    return "".join(q if i % 2 else re.sub(r"\s+", "", q) for i, q in enumerate(parts))


def _same_text(ours: str | None, theirs: Any) -> bool:
    if is_site_value(theirs):  # masked in the stored extract: compared by hash
        return site_placeholder(ours or "") == theirs
    return _norm(ours) == _norm(theirs)


def same_python(a: Any, b: Any) -> bool:
    """Equality of plain Python values that tells bool from number and NaN equals NaN."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, int | float) and isinstance(b, int | float):
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
            return True
        return a == b
    if isinstance(a, tuple) and isinstance(b, tuple):
        return len(a) == len(b) and all(same_python(x, y) for x, y in zip(a, b, strict=True))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same_python(a[k], b[k]) for k in a)
    return type(a) is type(b) and a == b


def _check_default(spec: FunctionSpec, p: ParamSpec, m: MatlabParam) -> list[str]:
    where = f"{spec.matlab_name}.{p.matlab}"
    if p.differs:
        return []
    if m.default is None:
        return [f"{where}: MATLAB has no default text"]
    if is_site_value(m.default) or site_pattern().search(m.default):
        return [f"{where}: MATLAB's default is site-specific; declare differs=<reason>"]
    if p.matlab_default is not None:
        if not _same_text(p.matlab_default, m.default):
            return [f"{where}: default text {p.matlab_default!r}, MATLAB has {m.default!r}"]
        return []
    try:
        mv = parse(m.default)
    except LiteralError:
        return [
            f"{where}: MATLAB's default {m.default!r} is not a literal; give "
            f"matlab_default={m.default!r} (and the Python default that matches it)"
        ]
    want = to_python(mv)
    got = p.default
    if not same_python(want, got):
        return [f"{where}: default {got!r}, MATLAB has {m.default} ({want!r})"]
    return []


def check_spec(spec: FunctionSpec, fn: MatlabFunction | None) -> list[str]:
    """Problems of one spec against MATLAB's extract of its function (empty: in step)."""
    name = spec.matlab_name
    try:
        params = spec.params
    except RegistryError as err:
        return [str(err)]
    if fn is None:
        return [f"{name}: no MATLAB function of that name in the tree or the extract"]
    out: list[str] = []
    ours = [p for p in params if not p.python_only]
    theirs = fn.params
    on = [p.matlab for p in ours]
    tn = [p.name for p in theirs]
    if [n.lower() for n in on] != [n.lower() for n in tn]:
        missing = [n for n in tn if n.lower() not in {x.lower() for x in on}]
        extra = [n for n in on if n.lower() not in {x.lower() for x in tn}]
        if missing or extra:
            out.append(f"{name}: parameters differ: missing {missing}, not in MATLAB {extra}")
        else:
            out.append(f"{name}: parameter order differs from MATLAB's: {tn}")
        return out
    for p, m in zip(ours, theirs, strict=True):
        where = f"{name}.{p.matlab}"
        if p.matlab != m.name:
            out.append(f"{where}: MATLAB spells it {m.name!r}")
        if m.kind == "positional":
            # no inputParser or arguments block: MATLAB fills defaults by nargin
            if p.kind == "parameter":
                out.append(f"{where}: positional in MATLAB, name-value here")
            if p.validator:
                out.append(f"{where}: MATLAB has no validator, the spec {p.validator!r}")
            continue
        if p.kind != m.kind:
            out.append(f"{where}: {p.kind} here, {m.kind} in MATLAB")
        if not _same_text(p.validator, m.validator):
            out.append(f"{where}: validator {p.validator!r}, MATLAB has {m.validator!r}")
        if p.kind != "required":
            out.extend(_check_default(spec, p, m))
    if fn.case_sensitive:
        out.append(f"{name}: MATLAB's inputParser is case sensitive; the registry is not")
    if fn.keep_unmatched != spec.keep_unmatched:
        out.append(f"{name}: KeepUnmatched is {fn.keep_unmatched} in MATLAB")
    if not fn.partial_matching or not fn.struct_expand:
        out.append(f"{name}: PartialMatching/StructExpand differ from inputParser's defaults")
    if fn.parser:
        out.extend(_check_parser(spec, ours, fn))
    return out


def _check_parser(spec: FunctionSpec, ours: list[ParamSpec], fn: MatlabFunction) -> list[str]:
    out = []
    parser = fn.parser or {}
    conv = parser.get("conversions", {})
    pdefaults = parser.get("defaults", {})
    for p, m in zip(ours, fn.params, strict=True):
        want = conv.get(m.name, "none")
        if p.mcc != want:
            out.append(
                f"{spec.matlab_name}.{p.matlab}: the mcc parser converts strings by {want!r}, "
                f"the spec by {p.mcc!r} (give mcc={want!r})"
            )
        pd = pdefaults.get(m.name)
        if m.kind != "required" and pd is not None and _norm(str(pd)) != _norm(str(m.default)):
            out.append(
                f"{spec.matlab_name}.{p.matlab}: MATLAB drift: the mcc parser's default is "
                f"{pd!r}, the function's {m.default!r}"
            )
    npos = parser.get("positional")
    nours = len([p for p in ours if p.kind != "parameter"])
    if npos is not None and npos != nours:
        out.append(
            f"{spec.matlab_name}: mccMaster passes {npos} positional arguments, the spec has "
            f"{nours} required and optional ones"
        )
    return out


def check_signature_defaults(spec: FunctionSpec) -> list[str]:
    """Problems on the Python side alone: the declaration fits the signature, defaults are
    canonical (tuples, not lists) and pass their own validators."""
    try:
        params = spec.params
    except RegistryError as err:
        return [str(err)]
    out = []
    for p in params:
        if p.default is _MISSING:
            continue
        if isinstance(p.default, list | dict | set):
            out.append(f"{spec.matlab_name}.{p.matlab}: mutable default {p.default!r}")
        try:
            p.check(p.default, spec.matlab_name)
        except Exception as err:  # noqa: BLE001 - reported, not raised
            out.append(f"{spec.matlab_name}.{p.matlab}: default fails its validator: {err}")
    return out


# ----------------------------------------------------------------------------------------------
# the stored extract (CI has no MATLAB tree)


def _redact(text: str | None, scan: re.Pattern[str] | None) -> Any:
    if text is not None and scan is not None and scan.search(text):
        return site_placeholder(text)
    return text


def write_extract(
    root: str | Path,
    names: Iterable[str],
    out: str | Path,
    commit: str = "",
    scan: re.Pattern[str] | None = None,
) -> dict[str, Any]:
    """Extract the inputParser blocks of ``names`` from a MATLAB tree into a JSON file.

    Texts that match ``scan`` (``site_pattern()``) are stored as
    ``<site value sha256:...>`` and compared by that hash.
    """
    tree = extract_tree(root, set(names))
    funcs = {}
    for name in sorted(tree, key=str.lower):
        d = tree[name].to_json()
        d["notes"] = [n for n in d["notes"] if not scan or not scan.search(n)]
        for p in d["params"]:
            p["default"] = _redact(p["default"], scan)
            p["validator"] = _redact(p["validator"], scan)
        if d.get("parser"):
            d["parser"]["defaults"] = {
                k: _redact(v, scan) for k, v in d["parser"]["defaults"].items()
            }
            d["parser"]["validators"] = {
                k: _redact(v, scan) for k, v in d["parser"]["validators"].items()
            }
        funcs[name] = d
    data = {"version": EXTRACT_VERSION, "oracle_commit": commit, "functions": funcs}
    Path(out).write_text(
        "{"
        + f'"version": {EXTRACT_VERSION}, "oracle_commit": {json.dumps(commit)}, "functions": {{\n'
        + ",\n".join(
            f"{json.dumps(k)}: {json.dumps(v, separators=(',', ':'))}" for k, v in funcs.items()
        )
        + "\n}}\n"
    )
    return data


def read_extract(path: str | Path) -> dict[str, MatlabFunction]:
    data = json.loads(Path(path).read_text())
    if data.get("version") != EXTRACT_VERSION:
        raise ValueError(f"{path}: extract version {data.get('version')}, want {EXTRACT_VERSION}")
    return {k: MatlabFunction.from_json(v) for k, v in data["functions"].items()}


def check_all(
    specs: Iterable[FunctionSpec], matlab: dict[str, MatlabFunction]
) -> dict[str, list[str]]:
    """{MATLAB name: problems} for every spec (empty lists when in step)."""
    out = {}
    for spec in specs:
        out[spec.matlab_name] = check_signature_defaults(spec) + check_spec(
            spec, matlab.get(spec.matlab_name)
        )
    return out
