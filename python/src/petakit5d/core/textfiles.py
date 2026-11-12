"""``readTextFile``, ``writeTextFile`` and ``writeJsonFile``, with MATLAB's results byte for byte.

Text is UTF-8 both ways (MATLAB's default encoding on Linux and macOS); bytes that are not
UTF-8 read as U+FFFD, as in MATLAB. Quirks kept on purpose (each has a test):

- ``read_text_file`` of an empty file returns ``[-1]`` (``fgetl``'s end-of-file value), as
  ``readTextFile`` does. ``\\r\\n`` ends a line like ``\\n``.
- ``write_text_file`` of a list writes every line followed by ``\\n`` (an empty list writes one
  ``\\n``); of a string, the string as it is.
- ``write_json_file`` passes the JSON text to ``fprintf`` as the format, as ``writeJsonFile``
  does: ``\\\\`` in the JSON becomes ``\\``, ``\\t`` a tab, and the file ends at the first ``%``
  that starts a conversion (``"x%03d.tif"`` truncates the file after ``"x``) or at an escape
  ``sprintf`` does not know (``\\u0001``).
"""

from __future__ import annotations

import os
import re

from petakit5d.core.jsonenc import jsonencode
from petakit5d.core.text import sprintf

__all__ = ["read_text_file", "write_json_file", "write_text_file"]


def read_text_file(filename: str | os.PathLike) -> list:
    """``readTextFile``: the lines of a text file, without line ends."""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{os.fspath(filename)} does not exist, please check the path!")
    with open(filename, "rb") as f:
        data = f.read()
    if not data:
        return [-1]
    text = data.decode("utf-8", errors="replace")
    # fgetl ends a line at \n, \r\n and a lone \r
    lines = re.split(r"\r\n|\r|\n", text)
    if lines[-1] == "":
        lines.pop()
    return lines


def write_text_file(text_lines, filename: str | os.PathLike, batch_size: int = 10000) -> None:
    """``writeTextFile``: a string as it is, or a list of lines each ended by ``\\n``.

    ``batch_size`` is accepted for compatibility (MATLAB writes long lists in batches; the
    bytes are the same).
    """
    if isinstance(text_lines, str):
        data = text_lines
    elif isinstance(text_lines, list | tuple):
        data = "\n".join(text_lines) + "\n"
    else:
        raise TypeError("text_lines must be a string or a list of strings")
    with open(filename, "wb") as f:
        f.write(data.encode("utf-8"))


def write_json_file(pr, filename: str | os.PathLike) -> None:
    """``writeJsonFile``: ``jsonencode(pr, PrettyPrint=true)`` written with ``fprintf(fid, s)``."""
    s = jsonencode(pr, pretty_print=True)
    with open(filename, "wb") as f:
        f.write(sprintf(s).encode("utf-8"))
