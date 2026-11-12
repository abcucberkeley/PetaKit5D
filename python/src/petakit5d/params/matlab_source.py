"""Reading MATLAB sources: inputParser blocks, mcc parser conversions, demo arguments.

Used by check_registry (names, defaults and validators of the registered functions against
the MATLAB tree), by the literal corpus (defaults of the mcc parsers, arguments of the demos)
and by the spec skeleton generator. Plain text processing: comments and ``...``
continuations are removed first, strings are kept intact.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from petakit5d.params.spec import mcc_mode  # noqa: F401 - also exported here

_ADD_METHODS = ("addRequired", "addOptional", "addParameter", "addParamValue")
_KIND = {
    "addRequired": "required",
    "addOptional": "optional",
    "addParameter": "parameter",
    "addParamValue": "parameter",
}
_VALUE_END = set("_)]}'\".")


def _is_value_end(c: str) -> bool:
    return c.isalnum() or c in _VALUE_END


def clean_source(text: str) -> str:
    """Remove comments (``%``, ``%{ ... %}``) and join ``...`` continuations.

    Strings stay as they are; a quote is a transpose after a value with no blank in between
    (or outside brackets after a blank, as MATLAB reads it).
    """
    out: list[str] = []
    i, n = 0, len(text)
    depth = 0  # [] and {} nesting, where a blank before ' starts a string
    lines_start = True
    block = 0
    while i < n:
        c = text[i]
        if lines_start:
            # block comments: a line with only %{ or %}
            j = text.find("\n", i)
            line = text[i : n if j < 0 else j].strip()
            if line == "%{":
                block += 1
                i = n if j < 0 else j + 1
                out.append("\n")
                continue
            if block and line == "%}":
                block -= 1
                i = n if j < 0 else j + 1
                out.append("\n")
                continue
            if block:
                i = n if j < 0 else j + 1
                out.append("\n")
                continue
            lines_start = False
        if c == "\n":
            out.append(c)
            lines_start = True
            i += 1
            continue
        if c == "%":
            j = text.find("\n", i)
            i = n if j < 0 else j
            continue
        if text.startswith("...", i):
            j = text.find("\n", i)
            out.append(" ")
            i = n if j < 0 else j + 1
            continue
        if c in "[{":
            depth += 1
        elif c in "]}" and depth:
            depth -= 1
        if c == '"' or (c == "'" and not _quote_is_transpose(out, depth)):
            j = i + 1
            while j < n and text[j] != "\n":
                if text[j] == c:
                    if j + 1 < n and text[j + 1] == c:
                        j += 2
                        continue
                    break
                j += 1
            out.append(text[i : j + 1])
            i = j + 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def _quote_is_transpose(out: list[str], depth: int) -> bool:
    k = len(out) - 1
    while k >= 0 and out[k] == "":
        k -= 1
    if k < 0:
        return False
    prev = out[k][-1]
    if prev in " \t":
        if depth:
            return False
        # outside brackets "a '" is still a transpose of a; look further back
        kk = k
        while kk >= 0 and out[kk][-1:] in (" ", "\t"):
            kk -= 1
        return kk >= 0 and _is_value_end(out[kk][-1]) and not _ends_with_keyword(out, kk)
    return _is_value_end(prev)


def _ends_with_keyword(out: list[str], k: int) -> bool:
    tail = "".join(out[max(0, k - 20) : k + 1])
    m = re.search(r"(\w+)$", tail)
    return bool(m) and m.group(1) in ("case", "otherwise", "return", "else", "end")


def split_args(text: str) -> list[str]:
    """Split an argument list at top-level commas (outside (), [], {} and strings)."""
    args: list[str] = []
    depth = 0
    start = 0
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == '"' or (c == "'" and not (i > 0 and _is_value_end(text[i - 1]))):
            j = i + 1
            while j < n:
                if text[j] == c:
                    if j + 1 < n and text[j + 1] == c:
                        j += 2
                        continue
                    break
                j += 1
            i = j
        elif c == "," and depth == 0:
            args.append(text[start:i].strip())
            start = i + 1
        i += 1
    tail = text[start:].strip()
    if tail or args:
        args.append(tail)
    return args


def balanced(text: str, open_pos: int) -> int:
    """Index of the bracket that closes the one at ``open_pos`` (strings skipped)."""
    pairs = {"(": ")", "[": "]", "{": "}"}
    stack = [pairs[text[open_pos]]]
    i, n = open_pos + 1, len(text)
    while i < n:
        c = text[i]
        if c in pairs:
            stack.append(pairs[c])
        elif c in ")]}":
            if c != stack[-1]:
                raise ValueError(f"unbalanced {c!r} at {i}")
            stack.pop()
            if not stack:
                return i
        elif c == '"' or (c == "'" and not _is_value_end(text[i - 1])):
            j = i + 1
            while j < n and text[j] != "\n":
                if text[j] == c:
                    if j + 1 < n and text[j + 1] == c:
                        j += 2
                        continue
                    break
                j += 1
            i = j
        i += 1
    raise ValueError("no closing bracket")


def unquote(text: str) -> str | None:
    """The value of a char or string literal, or None if ``text`` is not one."""
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        q = text[0]
        body = text[1:-1]
        if q + q in body.replace(q + q, ""):
            return None
        if q in body.replace(q + q, ""):
            return None
        return body.replace(q + q, q)
    return None


@dataclass
class MatlabParam:
    kind: str  # required, optional, parameter
    name: str
    default: str | None  # MATLAB text; None for required arguments
    validator: str | None  # MATLAB text; None if not given


@dataclass
class MatlabFunction:
    name: str
    file: str
    args: list[str]
    params: list[MatlabParam]
    case_sensitive: bool = False
    keep_unmatched: bool = False
    partial_matching: bool = True
    struct_expand: bool = True
    parser: dict[str, Any] | None = None  # the mcc parser file, if any
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> dict[str, Any]:
        d = asdict(self)
        return d

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> MatlabFunction:
        d = dict(d)
        d["params"] = [MatlabParam(**p) for p in d["params"]]
        return cls(**d)


_FUNC_RE = re.compile(
    r"^[ \t]*function\b\s*(?:(?:\[[^\]]*\]|\w+)\s*=\s*)?(\w+)\s*(?:\(([^)]*)\))?", re.M
)


def main_function_text(clean: str) -> tuple[str, str, list[str]] | None:
    """(name, text of the main function, argument names) of a cleaned function file."""
    m = _FUNC_RE.search(clean)
    if not m:
        return None
    name = m.group(1)
    args = [a.strip() for a in (m.group(2) or "").split(",") if a.strip()]
    nxt = _FUNC_RE.search(clean, m.end())
    body = clean[m.start() : nxt.start() if nxt else len(clean)]
    return name, body, args


def _setting(body: str, var: str, prop: str, default: bool) -> bool:
    m = re.search(rf"\b{re.escape(var)}\.{prop}\s*=\s*(true|false|1|0)\b", body)
    if not m:
        return default
    return m.group(1) in ("true", "1")


def _split_top_eq(text: str) -> tuple[str, str | None]:
    """Split ``decl = default`` at the first top-level ``=`` (not ``==``)."""
    depth = 0
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c in "'\"":
            j = text.find(c, i + 1)
            i = n if j < 0 else j
        elif c == "=" and depth == 0:
            if text[i + 1 : i + 2] == "=" or text[i - 1 : i] in "~<>=":
                i += 2
                continue
            return text[:i].strip(), text[i + 1 :].strip()
        i += 1
    return text.strip(), None


def extract_arguments_block(
    path: str | Path, root: str | Path | None = None, clean: str | None = None
) -> MatlabFunction | None:
    """The function's ``arguments`` block (positional and ``options.<name>`` name-value
    arguments), or, without one, its positional arguments; None for a script.

    The validator text of an argument is its declaration (``(1, :) {mustBeNumeric}``);
    arguments without a block have kind "positional" (MATLAB fills their defaults by nargin).
    """
    path = Path(path)
    if clean is None:
        clean = clean_source(path.read_text(encoding="utf-8", errors="replace"))
    found = main_function_text(clean)
    if not found:
        return None
    name, body, args = found
    rel = (str(path.relative_to(root)) if root else path.name).replace("\\", "/")
    params: list[MatlabParam] = []
    m = re.search(r"^[ \t]*arguments\b[^\n]*\n(.*?)^[ \t]*end\b", body, re.S | re.M)
    if m:
        for line in m.group(1).split("\n"):
            line = line.strip().rstrip(";")
            if not line:
                continue
            lm = re.match(r"^(\w+)(?:\.(\w+))?\s*(.*)$", line)
            if not lm:
                continue
            decl, default = _split_top_eq(lm.group(3))
            if lm.group(2):
                params.append(MatlabParam("parameter", lm.group(2), default, decl or None))
            else:
                kind = "required" if default is None else "optional"
                params.append(MatlabParam(kind, lm.group(1), default, decl or None))
        return MatlabFunction(name, rel, args, params, notes=["arguments block"])
    for a in args:
        if a != "varargin":
            params.append(MatlabParam("positional", a, None, None))
    return MatlabFunction(name, rel, args, params, notes=["positional arguments only"])


def extract_inputparser(path: str | Path, root: str | Path | None = None) -> MatlabFunction | None:
    """The inputParser block of the main function of a MATLAB file, or None."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    if "inputParser" not in text:
        return None
    clean = clean_source(text)
    found = main_function_text(clean)
    if not found:
        return None
    name, body, args = found
    m = re.search(r"\b(\w+)\s*=\s*inputParser\b", body)
    if not m:
        return None
    var = m.group(1)
    params: list[MatlabParam] = []
    notes: list[str] = []
    pat = re.compile(rf"\b{re.escape(var)}\s*\.\s*({'|'.join(_ADD_METHODS)})\s*\(")
    for cm in pat.finditer(body):
        method = cm.group(1)
        op = cm.end() - 1
        try:
            close = balanced(body, op)
        except ValueError:
            notes.append(f"unbalanced {method} call")
            continue
        parts = split_args(body[op + 1 : close])
        pname = unquote(parts[0]) if parts else None
        if pname is None:
            notes.append(f"{method} with a computed name: {parts[:1]}")
            continue
        kind = _KIND[method]
        if kind == "required":
            default = None
            validator = parts[1] if len(parts) > 1 else None
        else:
            default = parts[1] if len(parts) > 1 else None
            validator = parts[2] if len(parts) > 2 else None
            if len(parts) > 3:
                notes.append(f"{pname}: extra arguments {parts[3:]}")
        params.append(MatlabParam(kind, pname, default, validator))
    rel = str(path.relative_to(root)) if root else path.name
    return MatlabFunction(
        name=name,
        file=rel.replace("\\", "/"),
        args=args,
        params=params,
        case_sensitive=_setting(body, var, "CaseSensitive", False),
        keep_unmatched=_setting(body, var, "KeepUnmatched", False),
        partial_matching=_setting(body, var, "PartialMatching", True),
        struct_expand=_setting(body, var, "StructExpand", True),
        notes=notes,
    )


_CONV_RE = re.compile(
    r"\bif\s+ischar\((\w+)\)(?P<cond>[^\n]*)\n\s*(\w+)\s*=\s*(str2num|eval)\(\s*(\w+)\s*\)\s*;"
)


def extract_parser_conversions(path: str | Path) -> dict[str, Any]:
    """Conversions of a generated ``<name>_parser.m``: {arg: "str2num" | "eval:<chars>"}."""
    clean = clean_source(Path(path).read_text(encoding="utf-8", errors="replace"))
    conv: dict[str, str] = {}
    for m in _CONV_RE.finditer(clean):
        arg = m.group(1)
        if m.group(3) != arg or m.group(5) != arg:
            continue
        if m.group(4) == "str2num":
            conv[arg] = "str2num"
        else:
            chars = re.findall(r"strcmp\(\s*\w+\(1\)\s*,\s*'(.)'\s*\)", m.group("cond"))
            conv[arg] = "eval:" + "".join(sorted(set(chars)))
    ip = extract_inputparser(path)
    defaults = {p.name: p.default for p in ip.params} if ip else {}
    validators = {p.name: p.validator for p in ip.params} if ip else {}
    return {"conversions": conv, "defaults": defaults, "validators": validators}


def extract_mcc_master(path: str | Path) -> dict[str, int]:
    """Number of positional arguments mccMaster.m passes to each parser function."""
    clean = clean_source(Path(path).read_text(encoding="utf-8", errors="replace"))
    out: dict[str, int] = {}
    for m in re.finditer(r"case\s+'(\w+)'\s*\n\s*(\w+)_parser\(([^\n]*)\);", clean):
        singles = re.findall(r"varargin\{(\d+)\}", m.group(3))
        out[m.group(1)] = len(singles)
    return out


def extract_tree(root: str | Path, names: set[str] | None = None) -> dict[str, MatlabFunction]:
    """inputParser blocks of every function in a MATLAB tree (``names``: only these).

    The mcc parser of a function (mcc/parsers/<name>_parser.m) is attached when it exists.
    Files named ``*_parser.m`` under mcc/ and ``*_bak.m`` files are not functions of their own.
    """
    root = Path(root)
    out: dict[str, MatlabFunction] = {}
    parsers_dir = root / "mcc" / "parsers"
    master = parsers_dir / "mccMaster.m"
    positional = extract_mcc_master(master) if master.exists() else {}
    for path in sorted(root.rglob("*.m")):
        rel = path.relative_to(root).as_posix()
        if rel.startswith("mcc/") or path.stem.endswith("_bak"):
            continue
        if names is not None and path.stem not in names:
            continue
        try:
            fn = extract_inputparser(path, root)
            if fn is None and names is not None:
                fn = extract_arguments_block(path, root)
        except (OSError, ValueError) as err:  # unreadable or malformed file
            fn = MatlabFunction(path.stem, rel, [], [], notes=[f"not read: {err}"])
        if fn is None:
            continue
        if fn.name != path.stem:
            fn.notes.append(f"file {path.stem}.m defines {fn.name}")
            fn.name = path.stem
        if fn.name in out:
            # the first in path order is kept; MATLAB's path decides which one runs
            out[fn.name].notes.append(f"also in {fn.file}")
            continue
        pfile = parsers_dir / f"{fn.name}_parser.m"
        if pfile.exists():
            p = extract_parser_conversions(pfile)
            p["file"] = pfile.relative_to(root).as_posix()
            p["positional"] = positional.get(fn.name)
            fn.parser = p
        out[fn.name] = fn
    return out


# ----------------------------------------------------------------------------------------------
# demo scripts: literal arguments


_ASSIGN_RE = re.compile(r"^[ \t]*([A-Za-z]\w*)\s*=\s*(?!=)(.+?)\s*;?\s*$", re.M)


def demo_literals(path: str | Path) -> list[dict[str, str]]:
    """Candidate literal texts of a demo script.

    Right-hand sides of simple assignments (``dz = 0.3;``) and the arguments of function
    calls written as statements, name-value pairs included (``'Reverse', true`` and
    ``Reverse=true``). Each item: {"kind": "assign" | "arg", "where": ..., "text": ...}.
    """
    clean = clean_source(Path(path).read_text(encoding="utf-8", errors="replace"))
    items: list[dict[str, str]] = []
    for line in clean.split("\n"):
        stmts = _split_statements(line)
        for st in stmts:
            m = _ASSIGN_RE.match(st)
            if m and not re.match(r"\s*\w+\s*\(", m.group(1)):
                rhs = m.group(2).strip()
                items.append({"kind": "assign", "where": m.group(1), "text": rhs})
                call = re.match(r"(\w+)\s*\(", rhs)
                if call:
                    items.extend(_call_args(rhs, call))
                continue
            call = re.match(r"\s*(?:\[[^\]]*\]\s*=\s*)?(\w+)\s*\(", st)
            if call:
                items.extend(
                    _call_args(
                        st.strip(), re.match(r"(?:\[[^\]]*\]\s*=\s*)?(\w+)\s*\(", st.strip())
                    )
                )
    return items


def _call_args(st: str, call: re.Match[str] | None) -> list[dict[str, str]]:
    if call is None:
        return []
    op = call.end() - 1
    try:
        close = balanced(st, op)
    except ValueError:
        return []
    fname = call.group(1)
    out = []
    for a in split_args(st[op + 1 : close]):
        nv = re.match(r"^([A-Za-z]\w*)\s*=\s*(?!=)(.+)$", a)
        if nv:
            out.append({"kind": "arg", "where": f"{fname}:{nv.group(1)}", "text": nv.group(2)})
        else:
            out.append({"kind": "arg", "where": fname, "text": a})
    return out


def _split_statements(line: str) -> list[str]:
    """Split a cleaned line at top-level ``;`` and ``,`` that end statements."""
    out: list[str] = []
    depth = 0
    start = 0
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth = max(depth - 1, 0)
        elif c == '"' or (c == "'" and not (i > 0 and _is_value_end(line[i - 1]))):
            j = i + 1
            while j < n:
                if line[j] == c:
                    if j + 1 < n and line[j + 1] == c:
                        j += 2
                        continue
                    break
                j += 1
            i = j
        elif c == ";" and depth == 0:
            out.append(line[start:i])
            start = i + 1
        i += 1
    if line[start:].strip():
        out.append(line[start:])
    return [s for s in out if s.strip()]
