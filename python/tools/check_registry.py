"""Check the parameter registry against MATLAB's inputParser blocks.

    python tools/check_registry.py [--matlab-root TREE | --extract FILE]
        every registered spec against a MATLAB tree (on the cluster: the oracle tree) or the
        stored extract (default; CI); exit 1 on any difference
    python tools/check_registry.py extract --matlab-root TREE [--commit C] [--out FILE]
        write the stored extract (planned, registered and mcc functions) from a MATLAB tree
    python tools/check_registry.py skeleton NAME [NAME ...] --matlab-root TREE
        print a spec skeleton (decorator and signature) for a MATLAB function

Run from python/ with PYTHONPATH=src (or the package installed).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
DEFAULT_EXTRACT = HERE / "tests" / "params" / "data" / "inputparser_extract.json"


def _names(root: Path) -> set[str]:
    from petakit5d.params import REGISTRY
    from petakit5d.params.planned import PLANNED

    names = {m for m, *_ in PLANNED} | {s.matlab_name for s in REGISTRY.specs()}
    parsers = root / "mcc" / "parsers"
    names |= {p.name[: -len("_parser.m")] for p in parsers.glob("*_parser.m")}
    return names


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("command", nargs="?", default="check", choices=["check", "extract", "skeleton"])
    ap.add_argument("names", nargs="*")
    ap.add_argument("--matlab-root", type=Path)
    ap.add_argument("--extract", type=Path, default=DEFAULT_EXTRACT)
    ap.add_argument("--out", type=Path, default=DEFAULT_EXTRACT)
    ap.add_argument("--commit", default="")
    a = ap.parse_args(argv)

    from petakit5d.params import REGISTRY
    from petakit5d.params.check import (
        check_all,
        read_extract,
        site_pattern,
        write_extract,
    )
    from petakit5d.params.matlab_source import extract_tree

    if a.command == "extract":
        if not a.matlab_root:
            ap.error("extract needs --matlab-root")
        data = write_extract(a.matlab_root, _names(a.matlab_root), a.out, a.commit,
                             site_pattern())  # fmt: skip
        print(f"{len(data['functions'])} functions -> {a.out}")
        return 0
    if a.command == "skeleton":
        from petakit5d.params.skeleton import skeleton

        if not a.matlab_root or not a.names:
            ap.error("skeleton needs names and --matlab-root")
        tree = extract_tree(a.matlab_root, set(a.names))
        for n in a.names:
            if n not in tree:
                print(f"# {n}: no inputParser block found", file=sys.stderr)
                continue
            print(skeleton(tree[n], Path(a.matlab_root) / tree[n].file))
        return 0

    specs = REGISTRY.specs()
    if a.matlab_root:
        matlab = extract_tree(a.matlab_root, {s.matlab_name for s in specs})
        source = str(a.matlab_root)
    else:
        matlab = read_extract(a.extract)
        source = str(a.extract)
    problems = check_all(specs, matlab)
    bad = {k: v for k, v in problems.items() if v}
    for msgs in bad.values():
        for m in msgs:
            print(f"DRIFT {m}")
    print(f"{len(specs)} registered functions checked against {source}: "
          f"{len(specs) - len(bad)} in step, {len(bad)} with differences")  # fmt: skip
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
