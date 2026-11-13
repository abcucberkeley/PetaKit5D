"""The literal corpus and its MATLAB references (a stand-in until the oracle harness, F04).

Corpus: every default of the mcc parser functions, every literal argument and assigned
literal of the demos (both read from the oracle tree), and hand-written edge cases. Texts
that name anything but the constants and constructors the parser knows are not literals
(variables, function calls) and are counted, not kept.

    python literal_corpus.py build <oracle tree> <corpus.json>
    (MATLAB) pkpy_f03_literals('<corpus.json>', '<matlab.json>')
    python literal_corpus.py manifest <corpus.json> <matlab.json> <literal_oracle.json>

The manifest (texts and MATLAB's descriptions) is what test_literal_oracle.py checks, in CI
without MATLAB. Site-specific texts (site_pattern: absolute paths, Slurm options) are left
out of the corpus and counted.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any

from petakit5d.params.check import site_pattern
from petakit5d.params.literal import LiteralError, format_literal, parse, str2num
from petakit5d.params.matlab_source import demo_literals, extract_inputparser
from petakit5d.params.values import describe

_SCAN = site_pattern()

_KNOWN = {
    "true", "false", "Inf", "inf", "NaN", "nan", "pi", "eps", "zeros", "ones", "cell", "char",
    "double", "single", "logical", "int8", "uint8", "int16", "uint16", "int32", "uint32",
    "int64", "uint64", "sind", "cosd", "sqrt", "abs", "floor", "ceil", "round", "fix",
}  # fmt: skip

EDGE_CASES = [
    "1", "-1", "+1", "-0", "0.1", ".5", "5.", "1e-3", "1E3", "3e8", "1e-310",
    "1.7976931348623157e308", "2.2250738585072014e-308", "9007199254740993", "0.1+0.2",
    "1/3", "2^10", "2^0.5", "10^-3",
    "-2^2", "2^-1", "1.1^3", "pi/2", "pi / 12", "[pi/2]", "Inf", "-Inf", "NaN", "eps",
    "[]", "[ ]", "{}", "''", '""', "[1,2,3]", "[1 2 3]", "[1 -2]", "[1 - 2]", "[1 , -2]",
    "[1 -  2]", "[1]", "[1;2;3]", "[1,2;3,4]", "[1 2\n3 4]", "[1, 2, ...\n 3]",
    "[1,2,]", "[1;;2]", "[[1,2],[3]]", "[[1;2],[3;4]]", "[true, false]", "[true, 1]",
    "[true; false]", "true", "false", "~true", "~[1 0 2]", "[1 2]'", "[1 2].'", "{1,2}'",
    "'it''s'", '"a""b"', "'a'", "['ab', 'cd']", "['ab'; 'cd']", "['a' 'b']", "['' 65]",
    "[65 'a']", "[[] 'a']", "[{1}, 2]", "[{1}, []]", "[{1}, {2}]", "[{1}; {2}]", "{1, 'a'}",
    "{'a'; 'b'}", "{1, 2; 3, 4}", "{{1}, {2, 3}}", "{[], 'a', {}}", "{1,}",
    "{'CamA_ch0', 'CamB_ch0'}",
    "{'', '', ''}", "0:3", "1:0", "0:0.1:1", "0:0.1:0.3", "1:2:10", "10:-3:1", "-1:0.25:1",
    "0:pi/4:pi", "1.5:4", "[1:3]", "[0:2:6]'", "zeros(1,0)", "zeros(0,3)", "zeros(2)", "ones(2,3)",
    "true(1,2)", "false(2,1)", "NaN(1,2)", "Inf(1,1)", "cell(1,2)", "cell(0,0)", "char(65)",
    "char([72 105])", "double('A')", "single(0.1)", "single([1.5 2.5])", "uint8(300)",
    "uint8(2.5)", "int8(-2.5)", "int16([1.5 -1.5])", "uint16(65535.6)", "logical([1 0 2])",
    "int32(7)", "single(1)", "[single(1), 2]", "[uint8(1), 2.7]", "[int8(5), single(3.6)]",
    "[1,2]*2", "[1 2] + [3 4]", "[1 2] .* [3 4]", "[2 4] ./ 2", "[1 2] .^ 2", "-[1 2]",
    "sind(32.4)", "0.3*sind(32.4)", "[0.11,0.11,0.2*sind(32.4)]", "sind(30)", "sind(45)",
    "cosd(45)", "cosd(60)", "sind(180)", "sind(-180)", "sind(-90)", "cosd(90)", "sind(270)",
    "sind(-32.45)", "cosd(32.45)", "sind(1e10)", "sqrt(2)", "abs(-3)", "floor(-0.5)",
    "ceil(-0.5)", "round(2.5)", "round(-0.4)", "fix(-2.5)", "floor([1.5 -1.5])",
    "0/0", "Inf-Inf", "-NaN", "[NaN 1]+1", "Inf*0", "1/0", "-1/0", "[1 2]/0", "1/-0",
    "[1,2;3]", "[1 2", "{1 2", "'abc", "1i", "x", "zeros(1,2)(1)", "a.b",
    # hexadecimal and binary literals (values up to 2^53: pkpy_f03_literals describes integers
    # as doubles), their classes, suffixes and errors, and integer classes in expressions
    "0x100", "0xFF", "0x0", "0x0001", "0x00FF", "0x000000FF", "0x10000", "0xFFFFFFFF",
    "0x100000000", "0b101", "0b11111111", "0b100000000", "0b0000000000000001", "0X1F", "0B11",
    "0xAbC", "0x1e3", "0xFFu16", "0x1u64", "0xFFs8", "0x80s8", "0x7Fs8", "0xFFs16", "0xFFFFs16",
    "0x80000000s32", "0b10000000s8", "0xFFU8", "0xFFS8", "0b00000001u8", "0x100u8", "0x00FFu8",
    "0b000000001u8", "0x0000000000000000001", "0xFFu7", "0xFFs", "0x", "0b", "0xG", "0b102",
    "0x1.5", "0x1F.", "0xFF_FF", "0b1e3", "00x10", "0xF'", "-0x10", "-0x10s8", "-0x80s8",
    "+0x10", "~0x10", "[0x10, 0x100]", "[0x100, 0x10]", "[1.5, 0x10]", "[0x10; 0x20]", "{0x10}",
    "[0x10 -1]", "[0x41 'a']", "[true, 0x2]", "[single(2.5), 0x10]", "[int8(-5), 0xFF]",
    "[0xFFs8, 300]", "[512, 0x100]", "[0x10 0x20]", "double(0x10)", "single(0x10)", "int8(0xFF)",
    "uint16(0xFFs8)", "logical(0x10)", "char(0x41)", "[int8(1), uint16(2)]",
    "[uint16(70000), int8(-3)]", "-int8(-128)", "~uint8([0 2])", "-[int8(1), int8(-128)]",
    "0x10 + 1", "0x1F:0x22",
]  # fmt: skip

STR2NUM_EXTRA = [
    "1 2 3", "1,2,3", "1;2", "true false", "true", "[1,2],[3]", "yes", "'abc'", "{1}",
    "{'a','b'}", "", " ", "1e-3", "[1024, 1024, 1024]", "-32.45", "Inf", "NaN", "false",
    "[true, false]", "1 -2", "1 - 2", "512 0x100", "0xFF, 0x100",
]  # fmt: skip


def _idents(text: str) -> set[str]:
    stripped = re.sub(r"'(?:[^'\n]|'')*'|\"(?:[^\"\n]|\"\")*\"", "''", text)
    stripped = re.sub(r"\d[\d.]*[eE][+-]?\d+", "0", stripped)
    return set(re.findall(r"\b[A-Za-z]\w*", stripped))


def build(root: Path) -> dict[str, Any]:
    items: list[dict[str, str]] = []
    seen: set[str] = set()
    counts: dict[str, Any] = {
        "parser_files": 0, "parser_defaults": 0, "parser_defaults_kept": 0,
        "parser_defaults_scan_dropped": 0, "parser_defaults_not_literal": [],
        "demo_files": 0, "demo_items": 0, "not_literal": 0, "scan_dropped": 0,
    }  # fmt: skip
    not_literal: dict[str, int] = {}

    def add(text: str, source: str, kind: str) -> str:
        text = text.strip()
        if not text:
            return "empty"
        if _SCAN.search(text):
            counts["scan_dropped"] += 1
            return "scan"
        extra = _idents(text) - _KNOWN
        if extra and kind != "edge":
            counts["not_literal"] += 1
            for name in extra:
                not_literal[name] = not_literal.get(name, 0) + 1
            return "not_literal"
        key = f"{kind}\0{text}"
        if key not in seen:
            seen.add(key)
            items.append({"text": text, "source": source, "kind": kind})
        return "kept"

    for pfile in sorted((root / "mcc" / "parsers").glob("*_parser.m")):
        fn = extract_inputparser(pfile, root)
        defaults = [p for p in (fn.params if fn else []) if p.default is not None]
        counts["parser_files"] += bool(defaults)
        for p in defaults:
            counts["parser_defaults"] += 1
            how = add(p.default, f"{pfile.name}:{p.name}", "eval")
            if how == "kept":
                counts["parser_defaults_kept"] += 1
            elif how == "scan":
                counts["parser_defaults_scan_dropped"] += 1
            else:
                counts["parser_defaults_not_literal"].append(f"{pfile.name}:{p.name}={p.default}")
    counts["demo_files"] = len(list((root / "demos").glob("*.m")))
    for demo in sorted((root / "demos").glob("*.m")):
        for it in demo_literals(demo):
            counts["demo_items"] += 1
            add(it["text"], f"{demo.name}:{it['where']}", "eval")
    for t in EDGE_CASES:
        add(t, "edge", "edge")
    for t in STR2NUM_EXTRA:
        key = f"str2num\0{t}"
        if key not in seen:
            seen.add(key)
            items.append({"text": t, "source": "str2num", "kind": "str2num"})
    for it in items:
        it.update(python_side(it["text"]))
    return {
        "counts": counts,
        "not_literal_names": dict(sorted(not_literal.items())),
        "items": items,
    }


def python_side(text: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    try:
        v = parse(text)
        out["py_eval"] = describe(v)
        out["formatted"] = format_literal(v)
    except LiteralError as err:
        out["py_eval"] = {"error": str(err)}
    x, ok, _ = str2num(text)
    out["py_str2num"] = {"value": describe(x), "ok": ok}
    return out


def manifest(corpus: dict[str, Any], matlab: list[dict[str, Any]]) -> dict[str, Any]:
    items = []
    for it, m in zip(corpus["items"], matlab, strict=True):
        if it["text"] != m["text"]:
            raise ValueError(f"corpus and MATLAB results out of step at {it['text']!r}")
        items.append(
            {
                "text": it["text"],
                "source": it["source"],
                "kind": it["kind"],
                "eval": m["eval"],
                "str2num": m["str2num"],
                "formatted": it.get("formatted"),
                "formatted_eval": m.get("formatted_eval"),
            }
        )
    return {"counts": corpus["counts"], "items": items}


def dump_lines(data: dict[str, Any]) -> str:
    """JSON with one item per line (small diffs when references change)."""
    head = {k: v for k, v in data.items() if k != "items"}
    items = ",\n".join(json.dumps(it, separators=(",", ":")) for it in data["items"])
    return json.dumps(head)[:-1] + ', "items": [\n' + items + "\n]}\n"


def main(argv: list[str]) -> int:
    if argv[0] == "build":
        data = build(Path(argv[1]))
        Path(argv[2]).write_text(json.dumps(data, indent=1) + "\n")
        print(json.dumps(data["counts"]), len(data["items"]), "items")
        print("not literal (names):", data["not_literal_names"])
        return 0
    if argv[0] == "manifest":
        corpus = json.loads(Path(argv[1]).read_text())
        matlab = json.loads(Path(argv[2]).read_text())
        Path(argv[3]).write_text(dump_lines(manifest(corpus, matlab)))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
