"""Print what goes into the keys of some oracle cases, to compare systems (a CI probe, not a test).

For each case: its key and, per input, the input's hash and every file it holds (relative path,
size, SHA-256), so that a key that differs between Linux, macOS and Windows names its file.

usage (from python/): python tests/platform/key_probe.py <test module>::<case id> ...
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
import sys
from pathlib import Path


def _module(path: str):
    spec = importlib.util.spec_from_file_location(Path(path).stem + "_probe", path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    sys.path.insert(0, str(Path(path).parent))
    spec.loader.exec_module(m)
    return m


def _cases(m) -> dict:
    out = {}
    for name, v in vars(m).items():
        if name.endswith("CASES") and isinstance(v, list | tuple):
            for c in v:
                if hasattr(c, "key"):
                    out.setdefault(c.id, c)
    return out


def _files(p: Path):
    if p.is_file():
        yield ".", p
        return
    for dp, dns, fns in sorted(os.walk(p)):
        dns.sort()
        for fn in sorted(fns):
            f = Path(dp) / fn
            yield f.relative_to(p).as_posix(), f


def main(args: list[str]) -> None:
    mods: dict[str, dict] = {}
    for a in args:
        mod, cid = a.split("::", 1)
        if mod not in mods:
            mods[mod] = _cases(_module(mod))
        c = mods[mod].get(cid)
        if c is None:
            print(f"{a}: no such case")
            continue
        print(f"== {a} key {c.key()[:12]}")
        for n, x in sorted(c.inputs.items()):
            line = f"   input {n} ({x.kind}) {x.hash()[:12]}"
            local = getattr(x, "local", None)
            if local is None:
                print(line)
                continue
            print(line)
            for rel, f in _files(Path(local())):
                b = f.read_bytes()
                print(f"      {rel} {len(b)} {hashlib.sha256(b).hexdigest()[:12]}")


if __name__ == "__main__":
    main(sys.argv[1:])
