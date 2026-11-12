"""A spec skeleton from a MATLAB function: the decorator and signature to paste and edit.

    python tools/check_registry.py skeleton XR_crop_dataset --matlab-root <MATLAB tree>

Names, validators and defaults come from the inputParser block, meanings from the header
comment; the function name from the plan (params/planned.py). Defaults that are not literals
get ``matlab_default=`` and a ``None`` to replace by hand.
"""

from __future__ import annotations

import keyword
import re
from pathlib import Path

from petakit5d.params.check import is_site_value, site_pattern
from petakit5d.params.kinds import KINDS, looks_converted
from petakit5d.params.literal import LiteralError, parse
from petakit5d.params.matlab_source import MatlabFunction
from petakit5d.params.planned import planned, planned_kind
from petakit5d.params.spec import snake_name
from petakit5d.params.values import to_python

_SCAN = site_pattern()
_DOC_RE = re.compile(r"^\s*%\s+(\w+)\s*:\s*(.+?)\s*$")


def header_docs(path: str | Path) -> dict[str, str]:
    """{name: meaning} from the comment lines ``%   name : type (default: x). meaning``."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    out: dict[str, str] = {}
    started = False
    for line in text.splitlines():
        s = line.strip()
        if not started:
            started = s.startswith("function")
            continue
        if not s:
            continue
        if not s.startswith("%"):
            break
        m = _DOC_RE.match(line)
        if m:
            doc = m.group(2)
            d = re.search(r"\(default:.*?\)\.?\s*(.*)$", doc)
            if d and d.group(1):
                doc = d.group(1)
            out.setdefault(m.group(1), doc.strip())
    return out


def python_name_of(matlab_name: str) -> tuple[str, str | None]:
    """(Python name, py override or None) of a MATLAB parameter name."""
    py = snake_name(matlab_name)
    if keyword.iskeyword(py):
        return py + "_", py + "_"
    return py, None


def skeleton(fn: MatlabFunction, path: str | Path | None = None, task: str = "") -> str:
    docs = header_docs(path) if path else {}
    plan = planned().get(fn.name)
    func = plan[0].rsplit(".", 1)[-1] if plan else snake_name(re.sub(r"^XR_", "", fn.name))
    task = task or (plan[1] if plan else "")
    lines = ["@register(", f'    "{fn.name}",']
    sig: list[str] = []
    star = False
    for p in fn.params:
        py, override = python_name_of(p.name)
        args = [repr(p.name), repr(p.validator) if p.validator else "None"]
        args.append(repr(docs.get(p.name, "")))
        if override:
            args.append(f"py={override!r}")
        kind = planned_kind(fn.name, p.name)
        conv = KINDS.get(kind) if kind else None
        if kind and not kind.startswith("?"):
            args.append(f"convert={kind!r}")
        elif kind is None and looks_converted(p.name):
            kind = "?"
        default = ""
        comment = f"  # convert: {kind[1:] or 'you'} decide(s)" if kind and kind[0] == "?" else ""
        required = p.kind in ("required", "positional")
        if not required and (is_site_value(p.default) or _SCAN.search(p.default or "")):
            # a site-specific MATLAB default never goes into the repository
            args.append('differs="site-specific MATLAB default (from the cluster config)"')
            default = '""'
        elif not required:
            try:
                v = to_python(parse(p.default or ""))
                default = repr(conv.from_matlab(v) if conv else v)
            except LiteralError:
                args.append(f"matlab_default={p.default!r}")
                default = "None"
                comment = f"  # MATLAB: {p.default}"
        lines.append(f"    Param({', '.join(args)}),")
        if p.kind == "parameter" and not star:
            sig.append("    *,")
            star = True
        sig.append(f"    {py},{comment}" if required else f"    {py}={default},{comment}")
    if task:
        lines.append(f'    task="{task}",')
    lines.append(")")
    lines.append(f"def {func}(")
    lines.extend(sig)
    lines.append("):")
    summary = (first_comment_line(path) if path else "") or "TODO: what it does, in one line."
    lines.append(f'    """{summary}"""')
    return "\n".join(lines) + "\n"


def first_comment_line(path: str | Path) -> str:
    """The first header comment line after the function line (the function's summary)."""
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines()[1:]:
        s = line.strip()
        if s.startswith("%"):
            return s.lstrip("% ").strip()
        if s:
            break
    return ""
