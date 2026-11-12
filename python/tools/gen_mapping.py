"""Generate docs/MATLAB_MAPPING.md from the parameter registry (never edit it by hand).

    python tools/gen_mapping.py            write the table
    python tools/gen_mapping.py --check    exit 1 if the file is not what the registry gives

Run from python/ with PYTHONPATH=src (or the package installed). The parameters of planned
functions come from the stored inputParser extract (tools/check_registry.py extract).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
TABLE = HERE / "docs" / "MATLAB_MAPPING.md"
EXTRACT = HERE / "tests" / "params" / "data" / "inputparser_extract.json"


def render(extract: Path = EXTRACT) -> str:
    from petakit5d.params import REGISTRY
    from petakit5d.params.check import read_extract
    from petakit5d.params.mapping import generate

    return generate(REGISTRY, read_extract(extract) if extract.exists() else {})


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--out", type=Path, default=TABLE)
    ap.add_argument("--extract", type=Path, default=EXTRACT)
    a = ap.parse_args(argv)
    text = render(a.extract)
    if a.check:
        current = a.out.read_text() if a.out.exists() else ""
        if current != text:
            print(f"{a.out} is out of date: run python tools/gen_mapping.py")
            return 1
        print(f"{a.out} is up to date")
        return 0
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(text)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
