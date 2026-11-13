"""``get_uuid`` and ``get_hostname``.

``get_uuid`` is a random version-4 UUID in MATLAB's form (36 lowercase characters, as
``java.util.UUID``); on Windows only its first 4 characters, as MATLAB keeps (long paths). It
never touches NumPy's or Python's global random state.
"""

from __future__ import annotations

import os
import socket
import sys
import uuid

__all__ = ["get_hostname", "get_uuid"]


def get_uuid() -> str:
    u = str(uuid.uuid4())
    if sys.platform.startswith("win"):
        u = u[:4]
    log = os.environ.get("PETAKIT5D_UUID_LOG")
    if log:
        # the oracle harness masks the uuids a case's run made (Windows' 4 characters are
        # otherwise not told from other names); appended, also from the run's subprocesses
        with open(log, "a", encoding="utf-8", newline="\n") as fh:
            fh.write(u + "\n")
    return u


def get_hostname(verbose: bool = True) -> str:
    """``get_hostname``: the host name, printed as ``Hostname: <name> `` (MATLAB's line).

    MATLAB returns ``echo $HOSTNAME``'s output with its newline; this returns the name alone.
    """
    if sys.platform.startswith("win"):
        name = socket.gethostname()
    else:
        name = os.environ.get("HOSTNAME") or socket.gethostname()
    if verbose:
        print(f"Hostname: {name.strip()} ")
    return name.strip()
