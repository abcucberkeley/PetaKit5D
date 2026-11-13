"""pk5d: the registered functions on the command line, with snake_case options and JSON values.

    pk5d list [--all]                      the registered functions (--all: and the unported)
    pk5d doc <function>                    a function's arguments, options and defaults
    pk5d <function> [value ...] [--option value ...]
    pk5d parse <function> <strings ...>    how pk5d-master binds MATLAB-style strings (JSON)
    pk5d master <function> <strings ...>   the same as pk5d-master
    pk5d mcc-shim <path> [--command CMD]   write a run_mccMaster.sh that runs pk5d-master

A function is named by its Python name (``crop_dataset``, ``petakit5d.tools.crop_dataset``)
or its MATLAB name (``XR_crop_dataset``). Values are JSON (``[512, 512, 512]``, ``true``,
``null``, ``"text"``); text that is not JSON is taken as a string, so paths need no quotes.
Values are the Python API's: boxes ``(z0, y0, x0, z1, y1, x1)`` 0-based and end-exclusive,
sizes and per-axis vectors ``(z, y, x)`` (docs/conventions.md). Options are the snake_case
parameter names, with ``-`` or ``_``: ``--batch-size [256,512,512]``, ``--pad=true``; an
option without a value is ``true``. ``--`` ends the options.

Exit status: 0; 1 when the function raised (with its traceback); 2 for a usage or argument
error. The result, when there is one, is printed as JSON (or its repr).
"""

from __future__ import annotations

import json
import os
import sys
from collections.abc import Sequence
from typing import Any, TextIO

__all__ = ["main", "find_spec", "parse_value", "parse_call", "doc_text", "list_text"]

USAGE = """usage: pk5d <function> [value ...] [--option value ...]
       pk5d list [--all] | doc <function> | parse <function> <strings ...>
       pk5d master <function> <strings ...> | mcc-shim <path> [--command CMD]
       pk5d --help | --version

Runs a function of PetaKit5D in Python with snake_case options and JSON values; see
'pk5d doc <function>' and docs/cli.md.
"""

EXIT_OK, EXIT_FAILED, EXIT_USAGE = 0, 1, 2


class UsageError(Exception):
    """A command line pk5d can't run (exit status 2)."""


def _registry(registry: Any) -> Any:
    if registry is not None:
        return registry
    from petakit5d.params.registry import REGISTRY

    return REGISTRY


def find_spec(name: str, registry: Any = None) -> Any:
    """A spec by MATLAB name (exact, else ignoring case) or Python name (short or dotted)."""
    reg = _registry(registry)
    if name in reg:
        return reg.find(name)
    try:
        return reg.find_python(name)
    except LookupError:
        pass
    low = name.lower()
    matches = [n for n in reg.names() if n.lower() == low]
    if len(matches) == 1:
        return reg.find(matches[0])
    try:
        reg.find(name)  # the registry's message: not ported, planned, or a suggestion
    except LookupError as err:
        raise UsageError(str(err)) from None
    raise UsageError(f"no registered function {name}")  # pragma: no cover


def parse_value(text: str) -> Any:
    """A command-line value: JSON, else the text itself."""
    try:
        return json.loads(text)
    except ValueError:
        return text


def _option_name(spec: Any, raw: str) -> str:
    name = raw.replace("-", "_")
    py = {p.py for p in spec.params}
    if name in py:
        return name
    low = {p.matlab.lower(): p.py for p in spec.params}
    hint = f" (did you mean --{low[name.lower()]}?)" if name.lower() in low else ""
    if not hint:
        import difflib

        m = difflib.get_close_matches(name, sorted(py), n=1, cutoff=0.6)
        hint = f" (did you mean --{m[0]}?)" if m else ""
    raise UsageError(f"{spec.python_name.rsplit('.', 1)[-1]}: unknown option --{raw}{hint}")


def parse_call(spec: Any, argv: Sequence[str]) -> tuple[list[Any], dict[str, Any]]:
    """Positional values and options of ``pk5d <function> ...`` (JSON values)."""
    args: list[Any] = []
    kwargs: dict[str, Any] = {}
    i, argv, options = 0, list(argv), True
    while i < len(argv):
        a = argv[i]
        if options and a == "--":
            options = False
        elif options and a.startswith("--") and len(a) > 2:
            raw, eq, val = a[2:].partition("=")
            name = _option_name(spec, raw)
            if name in kwargs:
                raise UsageError(f"--{raw} is given twice")
            if eq:
                kwargs[name] = parse_value(val)
            elif i + 1 < len(argv) and not argv[i + 1].startswith("--"):
                i += 1
                kwargs[name] = parse_value(argv[i])
            else:
                kwargs[name] = True
        else:
            args.append(parse_value(a))
        i += 1
    return args, kwargs


def _jsonable(x: Any) -> Any:
    if isinstance(x, os.PathLike):
        return os.fspath(x)
    if type(x).__module__ == "numpy" and hasattr(x, "tolist"):
        if getattr(x, "size", 1) <= 4096:
            return x.tolist()
        return f"<array {x.dtype} {tuple(x.shape)}>"
    if isinstance(x, set | frozenset):
        return sorted(x, key=repr)
    raise TypeError(type(x).__name__)


def _json(x: Any) -> str:
    try:
        return json.dumps(x, default=_jsonable)
    except (TypeError, ValueError):
        return repr(x)


def _short(spec: Any) -> str:
    return spec.python_name.rsplit(".", 1)[-1]


def list_text(registry: Any = None, include_unported: bool = False) -> str:
    """The registered functions: Python name, MATLAB name, summary."""
    reg = _registry(registry)
    specs = reg.specs()
    rows = [(_short(s), s.matlab_name, s.python_name, s.doc) for s in specs]
    w0 = max((len(r[0]) for r in rows), default=0)
    w1 = max((len(r[1]) for r in rows), default=0)
    lines = [f"{a:<{w0}}  {b:<{w1}}  {d}" for a, b, _, d in sorted(rows)]
    if include_unported:
        lines += ["", "Not ported:"]
        for n in reg.not_ported_all().values():
            lines.append(f"  {n.matlab_name}: {n.reason}")
    return "\n".join(lines) + "\n"


def doc_text(name: str, registry: Any = None) -> str:
    """A function's arguments and options as ``pk5d`` takes them, with MATLAB's names."""
    spec = find_spec(name, registry)
    lines = [
        f"{_short(spec)}  ({spec.python_name}; MATLAB {spec.matlab_name})",
        "",
        spec.doc,
        "",
        f"usage: pk5d {_short(spec)} "
        + " ".join(
            p.py if p.required else f"[{p.py}]" for p in spec.params if p.kind != "parameter"
        )
        + " [--option value ...]",
        "",
        f"MATLAB: pk5d-master {spec.matlab_name} <strings>, with MATLAB's names and values",
    ]
    for title, kinds in (("Arguments", ("required", "optional")), ("Options", ("parameter",))):
        ps = [p for p in spec.params if p.kind in kinds]
        if not ps:
            continue
        lines += ["", f"{title}:"]
        for p in ps:
            name = p.py if p.kind != "parameter" else "--" + p.py.replace("_", "-")
            dflt = "required" if p.required else f"default {_json(p.default)}"
            notes = [f"MATLAB {p.matlab}"]
            if p.conv is not None and p.conv.name != "same":
                notes.append(f"converted: {p.conv.name}")
            if p.python_only:
                notes.append("Python only")
            if p.differs:
                notes.append(f"differs from MATLAB: {p.differs}")
            lines.append(f"  {name}  ({dflt}; {'; '.join(notes)})")
            if p.doc:
                lines.append(f"      {p.doc}")
    return "\n".join(lines) + "\n"


def _parse_text(name: str, strings: Sequence[str], registry: Any) -> str:
    from petakit5d.cli import master
    from petakit5d.params.values import to_plain

    spec = master.find_function(name, registry)
    kwargs = master.bind(spec, strings)
    res = spec.matlab_results(kwargs)
    given = {p.matlab for p in spec.params if p.py in kwargs}
    out = {
        "function": spec.matlab_name,
        "python": spec.python_name,
        "python_arguments": kwargs,
        "matlab_values": {k: to_plain(v) for k, v in res.items() if k in given},
    }
    return json.dumps(out, indent=1, default=_jsonable) + "\n"


def _run(spec: Any, argv: Sequence[str], out: TextIO, err: TextIO) -> int:
    args, kwargs = parse_call(spec, argv)
    try:
        bound = spec.bind_python(tuple(args), kwargs)
    except (TypeError, ValueError) as e:
        raise UsageError(str(e)) from None
    try:
        result = spec.func(**bound)
    except Exception:  # noqa: BLE001 - reported with the exit status
        import traceback

        traceback.print_exc(file=err)
        return EXIT_FAILED
    if result is not None:
        out.write(_json(result) + "\n")
    return EXIT_OK


def main(
    argv: Sequence[str] | None = None,
    *,
    registry: Any = None,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> int:
    """``pk5d``; returns the exit status."""
    argv = list(sys.argv[1:] if argv is None else argv)
    out = stdout if stdout is not None else sys.stdout
    err = stderr if stderr is not None else sys.stderr
    if not argv or argv[0] in ("-h", "--help", "help"):
        out.write(USAGE)
        return EXIT_OK if argv else EXIT_USAGE
    cmd, rest = argv[0], argv[1:]
    try:
        if cmd == "--version":
            from petakit5d import __version__

            out.write(f"pk5d (petakit5d {__version__})\n")
            return EXIT_OK
        if cmd == "list":
            if rest not in ([], ["--all"]):
                raise UsageError("usage: pk5d list [--all]")
            out.write(list_text(registry, include_unported=bool(rest)))
            return EXIT_OK
        if cmd == "doc":
            if len(rest) != 1:
                raise UsageError("usage: pk5d doc <function>")
            out.write(doc_text(rest[0], registry))
            return EXIT_OK
        if cmd == "master":
            from petakit5d.cli import master

            return master.main(rest, registry=registry, stdout=out, stderr=err)
        if cmd == "parse":
            if not rest:
                raise UsageError("usage: pk5d parse <function> <strings ...>")
            from petakit5d.cli.master import MasterError
            from petakit5d.params.spec import MatlabArgumentError

            try:
                out.write(_parse_text(rest[0], rest[1:], registry))
            except (MasterError, MatlabArgumentError) as e:
                ident = getattr(e, "identifier", "")
                raise UsageError(f"{e}" + (f" ({ident})" if ident else "")) from None
            return EXIT_OK
        if cmd == "mcc-shim":
            return _shim(rest, out)
        if cmd.startswith("-"):
            raise UsageError(f"unknown option {cmd}")
        return _run(find_spec(cmd, registry), rest, out, err)
    except UsageError as e:
        err.write(f"pk5d: error: {e}\n")
        return EXIT_USAGE


def _shim(rest: Sequence[str], out: TextIO) -> int:
    from petakit5d.cli import master

    rest = list(rest)
    command = None
    if "--command" in rest:
        i = rest.index("--command")
        if i + 1 >= len(rest):
            raise UsageError("--command needs a value")
        command = rest[i + 1]
        del rest[i : i + 2]
    if len(rest) != 1:
        raise UsageError("usage: pk5d mcc-shim <path> [--command CMD]")
    path = master.write_shim(rest[0], command)
    out.write(f"wrote {path}\n")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
