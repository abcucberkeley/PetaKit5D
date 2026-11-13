"""pk5d-master: the command line of the compiled MATLAB master (mcc/parsers/mccMaster.m).

    pk5d-master <function> <argument strings>

takes what the GUI and PyPetaKit5D pass to MATLAB's ``run_mccMaster.sh`` after its Runtime
argument: a MATLAB function name (as in mccMaster's list, exact case) and strings such as
``"{'/mnt/projects/a','/mnt/projects/b'}"``, ``"[1,2,3]"``, ``true`` or ``Cropped``. As
mccMaster does, the first strings go to the function's required and optional arguments by
position, the rest are name-value pairs (names in any case, unique prefixes), and each string
is converted by its parameter's rule in the mcc parser files (``str2num``, ``eval`` of a
leading ``{``, ``[`` or ``@``, or kept as text), then checked by the function's validator
(``FunctionSpec.bind_strings``). The function gets MATLAB's values converted to the port's
(``petakit5d.params.kinds``). A function whose parser takes no name-value pairs ignores
strings beyond its arguments, as mccMaster passes it only those.

Output and exit status are the compiled master's: ``Hostname: <host> ``, ``Function name:
<function>`` and ``Variables: <strings joined by spaces> `` first, ``Elapsed time is <s>
seconds.`` and an empty line after the task; exit status 0, or 249 (the Linux compiled
master's status for an uncaught error) with the error on stderr as MATLAB prints it: ``Error
using <where>``, the message, an empty line and the error identifier. Errors of Python code
also print their traceback first. Any registered function is accepted, not only mccMaster's.

``run_mccMaster.sh <Runtime> <function> <strings>`` is the launcher the GUI and PyPetaKit5D
start; ``shim_text`` / ``write_shim`` (``pk5d mcc-shim <path>``) make one that runs this
master instead of the compiled one, with the launcher's own argument handling (each string
quoted and evaluated by the shell, which undoes PyPetaKit5D's escapes and turns the GUI's
``""`` into an empty string). The Runtime argument is ignored.

``PETAKIT5D_MASTER_DRY_RUN=1`` binds and prints the call instead of making it.

Start-up: nothing but the registry is imported before dispatch; the function's own module is
imported when it is looked up.
"""

from __future__ import annotations

import os
import sys
import time
from collections.abc import Sequence
from typing import Any, TextIO

__all__ = [
    "DRY_RUN_ENV",
    "EXIT_ERROR",
    "MasterError",
    "bind",
    "find_function",
    "header",
    "hostname",
    "main",
    "shim_text",
    "write_shim",
]

EXIT_ERROR = 249  # the exit status of the Linux compiled master after an uncaught error
DRY_RUN_ENV = "PETAKIT5D_MASTER_DRY_RUN"
SHIM_ENV = "PETAKIT5D_MASTER"  # overrides the master command a shim runs

USAGE = """usage: pk5d-master <function> [argument strings ...]

The compiled MATLAB master's command line (mccMaster): a MATLAB function name and MATLAB-style
strings ("{'a','b'}", "[1,2,3]", true, text), required and optional arguments first, then
name-value pairs. See docs/cli.md; 'pk5d list' lists the functions, 'pk5d doc <function>'
their parameters.
"""


class MasterError(Exception):
    """An error raised by the master itself (as mccMaster raises it), with MATLAB's identifier."""

    def __init__(self, message: str, identifier: str = "", where: str | None = "mccMaster"):
        super().__init__(message)
        self.identifier = identifier
        self.where = where


def hostname() -> str:
    """The host name as mccMaster prints it (get_hostname: ``echo $HOSTNAME``, ``hostname``)."""
    if os.name == "nt":
        import platform

        return platform.node()
    return os.environ.get("HOSTNAME") or os.uname().nodename


def header(function_name: str, args: Sequence[str]) -> str:
    """mccMaster's lines before the task: ``Function name: %s\\nVariables: %s \\n``."""
    return f"Function name: {function_name}\nVariables: {' '.join(args)} \n"


def find_function(name: str, registry: Any = None) -> Any:
    """The spec of a function by its MATLAB name (exact, as mccMaster's switch)."""
    if registry is None:
        from petakit5d.params.registry import REGISTRY as registry
    try:
        return registry.find(name)
    except LookupError as err:
        raise MasterError(f"The parser for function {name} does not exist!\n{err}") from None


def bind(spec: Any, args: Sequence[str]) -> dict[str, Any]:
    """The Python keyword arguments of a call with mccMaster's argument strings.

    mccMaster passes ``varargin{1..n}`` to the n required and optional arguments of the
    function's parser and ``varargin{n+1:end}`` only where the parser takes name-value pairs;
    the strings are converted and checked by ``FunctionSpec.bind_strings``. A name-value pair
    that names one of the n arguments is ignored, as in the parser (MATLAB quirk).
    """
    args = [str(a) for a in args]
    npos = sum(1 for p in spec.params if p.kind != "parameter")
    if not any(p.kind == "parameter" for p in spec.params):
        return spec.bind_strings(args[:npos])
    # a pair that names a positional argument passes the parser's inputParser, but the parser
    # hands the function its positional variable: the pair is ignored
    rest = args[npos:]
    keep = args[:npos]
    for i in range(0, len(rest), 2):
        pair = rest[i : i + 2]
        p = _named(spec, pair[0])
        if len(pair) == 2 and p is not None and p.kind != "parameter":
            continue
        keep += pair
    return spec.bind_strings(keep)


def _named(spec: Any, name: str) -> Any:
    """The parameter a name-value name means (exact in any case, else a unique prefix)."""
    low = name.lower()
    exact = [p for p in spec.params if p.matlab.lower() == low]
    if exact:
        return exact[0]
    pre = [p for p in spec.params if p.matlab.lower().startswith(low)]
    return pre[0] if len(pre) == 1 else None


def _call(spec: Any, kwargs: dict[str, Any]) -> Any:
    hook = spec.matlab_in
    if hook is not None:
        kwargs = hook(kwargs)
    return spec.func(**kwargs)


def _error_text(where: str | None, message: str, identifier: str) -> str:
    # MATLAB's display of an uncaught error in a compiled application (stack lines left out)
    lines = [f"Error using {where}"] if where else []
    lines += [message, "", identifier]
    return "\n".join(lines) + "\n"


def _dry_run_text(spec: Any, kwargs: dict[str, Any]) -> str:
    args = "".join(f"    {k}={v!r},\n" for k, v in kwargs.items())
    return f"Dry run ({DRY_RUN_ENV}): not called\n{spec.python_name}(\n{args})\n"


def main(
    argv: Sequence[str] | None = None,
    *,
    registry: Any = None,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> int:
    """Run one task as the compiled master does; returns the exit status (0 or 249)."""
    argv = list(sys.argv[1:] if argv is None else argv)
    out = stdout if stdout is not None else sys.stdout
    err = stderr if stderr is not None else sys.stderr
    if argv and argv[0] in ("-h", "--help"):
        out.write(USAGE)
        return 0
    out.write(f"Hostname: {hostname()} \n")
    if not argv:
        out.flush()
        err.write(_error_text(None, "Not enough input arguments.", "MATLAB:minrhs"))
        return EXIT_ERROR
    name, args = argv[0], argv[1:]
    out.write(header(name, args))
    out.flush()
    t0 = time.perf_counter()
    from petakit5d.params.spec import MatlabArgumentError

    try:
        spec = find_function(name, registry)
        kwargs = bind(spec, args)
    except MasterError as e:
        err.write(_error_text(e.where, str(e), e.identifier))
        return EXIT_ERROR
    except MatlabArgumentError as e:
        err.write(_error_text(name, str(e), e.identifier))
        return EXIT_ERROR
    if os.environ.get(DRY_RUN_ENV, "") not in ("", "0"):
        out.write(_dry_run_text(spec, kwargs))
        out.flush()
    else:
        try:
            _call(spec, kwargs)
        except Exception as e:  # noqa: BLE001 - reported as the compiled master reports errors
            import traceback

            out.flush()
            traceback.print_exc(file=err)
            ident = str(getattr(e, "identifier", "") or "")
            err.write("\n" + _error_text(name, f"{type(e).__name__}: {e}", ident))
            return EXIT_ERROR
    out.write(f"Elapsed time is {time.perf_counter() - t0:.6f} seconds.\n\n")
    out.flush()
    return 0


# ----------------------------------------------------------------------------------------------
# the run_mccMaster.sh shim

_SHIM = """#!/bin/sh
# run_mccMaster.sh for PetaKit5D in Python: MATLAB's launcher of the compiled master, with the
# same arguments (<deployedMCRroot> <function> <arguments>) and the same handling of them, but
# it runs pk5d-master. The MATLAB Runtime folder is not used. Made by 'pk5d mcc-shim'.
exe_name=$0
exe_dir=`dirname "$0"`
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \\<deployedMCRroot\\> args
else
  pk5d_master={command}
  if [ -n "${env}" ]; then
    pk5d_master=${env}
  fi
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${{args}} \\"${{token}}\\""
      shift
  done
  eval "$pk5d_master" $args
fi
exit
"""


def _sh_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def shim_text(command: str | None = None) -> str:
    """A run_mccMaster.sh that runs ``command`` (shell words; default: this interpreter with
    ``-m petakit5d.cli.master``) on the evaluated arguments after the Runtime folder.
    ``$PETAKIT5D_MASTER``, when set, replaces the command."""
    if command is None:
        command = f"{_sh_quote(sys.executable)} -m petakit5d.cli.master"
    return _SHIM.format(command=_sh_quote(command), env=SHIM_ENV)


def write_shim(path: str | os.PathLike, command: str | None = None) -> str:
    """Write ``shim_text(command)`` to ``path`` (executable); returns the path."""
    path = os.fspath(path)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(shim_text(command))
    if os.name != "nt":
        os.chmod(path, 0o755)
    return path


if __name__ == "__main__":
    sys.exit(main())
