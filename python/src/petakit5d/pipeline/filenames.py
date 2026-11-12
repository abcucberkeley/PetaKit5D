"""File names of a dataset: listing, channel patterns, partial volumes and the name grammar.

=======================================  ===================================================
``parse_image_filenames``                parseImageFilenames (utils)
``parse_dataset_filenames``              XR_parseImageFilenames (the pipeline's file list)
``group_image_filenames_by_channels``    groupImageFilenamesByChannels
``matches_channel_pattern``              ``contains(s, p) | contains(s, regexpPattern(p))``
``matlab_regex``                         a MATLAB regular expression as a compiled ``re``
``parse_filename_fields``                the tile name grammar of XR_parseSettingFiles_wrapper
=======================================  ===================================================

Indices in results are 0-based (``fd_inds[i]`` indexes ``data_paths``); where MATLAB uses 0 for
"none" (``FTP_inds``), Python uses -1. Everything else keeps MATLAB's values and order.

MATLAB's regular expressions are translated, not passed to ``re`` as they are; the
differences (each pinned by an oracle case in tests/pipeline/test_filenames.py):

- ``.`` matches a newline; ``$`` is the end of the text only (``re``: also before a final
  newline); ``\\b`` is a backspace, ``\\<`` and ``\\>`` are the word boundaries;
- ``\\s`` is ``[ \\f\\n\\r\\t\\v]`` (``re`` adds the separators ``\\x1c-\\x1f`` and Unicode
  spaces); unknown escapes are the character itself (``\\q`` is ``q``; ``re`` fails);
- ``\\xN`` and ``\\oN`` take every hex or octal digit that follows (``re``: two hex digits);
- ``(?<name>...)`` and ``\\k<name>`` are named tokens (``re``: ``(?P<name>...)``,
  ``(?P=name)``); ``(?i)`` in the middle applies to the rest of its alternative (``re``
  refuses it);
- ``[]`` is an empty class and ``[^]`` any character (``re``: ``]`` is a member); a ``[``
  inside a class is a member that opens a level, closed by a ``]`` (a member too);
- a ``)`` without ``(``, a quantifier with nothing to repeat (``*a``, the second ``*`` of
  ``a**``) and ``{,n}`` are literal text (``re`` fails or reads ``{0,n}``);
- an unmatched ``(``, ``[`` or ``{`` is the character itself for ``regexp``, and an error
  with MATLAB's identifier for ``regexpPattern`` (``strict``).
"""

from __future__ import annotations

import fnmatch
import functools
import math
import os
import re
import statistics
import time
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from petakit5d.core.paths import fileparts, matlab_dir, matlab_sort_key

__all__ = [
    "DatasetFilenames",
    "ImageFilenames",
    "MatlabError",
    "contains_text",
    "group_image_filenames_by_channels",
    "matches_channel_pattern",
    "matlab_regex",
    "parse_dataset_filenames",
    "parse_filename_fields",
    "parse_image_filenames",
]


class MatlabError(Exception):
    """An error where MATLAB raises one; ``identifier`` is MATLAB's error identifier."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


class _MatlabWarning(UserWarning):
    pass


def _warn(message: str, identifier: str = "") -> None:
    w = _MatlabWarning(message)
    w.identifier = identifier  # type: ignore[attr-defined]
    warnings.warn(w, stacklevel=3)


# ---------------------------------------------------------------------------------------------
# MATLAB regular expressions
# ---------------------------------------------------------------------------------------------

_SPACE = r" \f\n\r\t\v"
# MATLAB's escapes of control characters ('\b' is the backspace, not a word boundary)
_CONTROL = {"a": "\x07", "b": "\x08", "f": "\f", "n": "\n", "r": "\r", "t": "\t", "v": "\v"}
_CLASS_ESC = {"w": r"\w", "W": r"\W", "d": r"\d", "D": r"\D"}
_HEX = "0123456789abcdefABCDEF"
_OCT = "01234567"
_FLAGS = set("imsx-")
_NAME = re.compile(r"[A-Za-z]\w*")
_COUNT = re.compile(r"\{(\d+)(,(\d*))?\}")


def _code_point(p: str, i: int, digits: str, base: int) -> tuple[str, int]:
    """``\\x{N}``, ``\\xN``, ``\\o{N}``, ``\\oN`` from position i (after the letter)."""
    if i < len(p) and p[i] == "{":
        j = p.find("}", i)
        if j > i + 1 and all(c in digits for c in p[i + 1 : j]):
            return chr(int(p[i + 1 : j], base)), j + 1
        return "", -1
    j = i
    while j < len(p) and p[j] in digits:
        j += 1
    if j == i:
        return "", -1
    return chr(int(p[i:j], base)), j


def _escape(p: str, i: int, in_class: bool) -> tuple[str, int]:
    """The Python text of the escape at p[i] == '\\'; returns it and the next position."""
    if i + 1 >= len(p):
        return r"\\", i + 1
    c = p[i + 1]
    if c in _CONTROL:
        return re.escape(_CONTROL[c]), i + 2
    if c in _CLASS_ESC:
        return _CLASS_ESC[c], i + 2
    if c == "s":
        return (_SPACE if in_class else f"[{_SPACE}]"), i + 2
    if c == "S":
        return (r"\S" if in_class else f"[^{_SPACE}]"), i + 2
    if c in "xo":
        ch, j = _code_point(p, i + 2, _HEX if c == "x" else _OCT, 16 if c == "x" else 8)
        if j > 0:
            return re.escape(ch), j
        return re.escape(c), i + 2
    if not in_class:
        if c == "<":
            return r"\b(?=\w)", i + 2
        if c == ">":
            return r"\b(?<=\w)", i + 2
        if c.isdigit() and c != "0":
            return "\\" + c, i + 2
        if c == "k" and p.startswith("<", i + 2):
            j = p.find(">", i + 3)
            if j > 0 and _NAME.fullmatch(p[i + 3 : j]):
                return f"(?P={p[i + 3 : j]})", j + 1
    return re.escape(c), i + 2


def _class_body(p: str, i: int) -> tuple[str, bool, int] | None:
    """The class at p[i] == '[': (its Python body, negated, next position), or None when it
    has no closing ']'. A ']' right after '[' or '[^' closes the class (MATLAB: ``[]`` is
    empty). A '[' inside is a member and opens a level that a ']' closes (also a member):
    the class ends at the ']' that closes the first level."""
    j = i + 1
    neg = j < len(p) and p[j] == "^"
    if neg:
        j += 1
    body: list[str] = []
    depth = 0
    while j < len(p):
        c = p[j]
        if c == "]":
            if depth == 0:
                return "".join(body), neg, j + 1
            depth -= 1
            body.append(r"\]")
            j += 1
            continue
        if c == "\\":
            t, j = _escape(p, j, True)
            body.append(t)
            continue
        if c == "[":
            depth += 1
            body.append(r"\[")
            j += 1
            continue
        # characters that are set operations for re (FutureWarning) are literal
        body.append("\\" + c if c in "&~|" else c)
        j += 1
    return None


def _char_class(p: str, i: int, strict: bool) -> tuple[str, int]:
    """The class starting at p[i] == '['; returns its Python text and the next position."""
    cls = _class_body(p, i)
    if cls is None:
        if strict:
            raise MatlabError(
                "Unmatched '['. Use '\\[' to specify literal '['.",
                "MATLAB:REGEXP:ErrorUnmatchedBracket",
            )
        return r"\[", i + 1  # regexp: an unclosed '[' is the character itself
    body, neg, j = cls
    if not body:
        return (r"[\s\S]" if neg else "(?!)"), j
    return f"[{'^' if neg else ''}{body}]", j


_UNMATCHED_PAREN = (
    "Unmatched '('. use '\\(' to specify literal '('.",
    "MATLAB:REGEXP:ErrorUnmatchedParen",
)


@functools.lru_cache(maxsize=4096)
def _translate(p: str, strict: bool) -> str:
    # regexp: an unclosed '(' is the character itself, and what follows it is read again
    literal: frozenset[int] = frozenset()
    while True:
        text, unclosed = _translate_once(p, strict, literal)
        if not unclosed:
            return text
        if strict:
            raise MatlabError(*_UNMATCHED_PAREN)
        literal = literal | frozenset(unclosed)


def _translate_once(p: str, strict: bool, literal: frozenset[int]) -> tuple[str, list[int]]:
    out: list[str] = []
    # open groups: [position of the '(' in p, the scoped flag groups ('(?i:') opened in the
    # group's current alternative]; the first is the pattern
    groups: list[list] = [[-1, []]]
    atom = False  # a quantifiable item precedes
    quant = 0  # 0: no quantifier just before; 1: a quantifier; 2: with its ? or + modifier
    i, n = 0, len(p)

    def open_group(py: str) -> None:
        out.append(py)
        groups.append([i, []])

    while i < n:
        c = p[i]
        if c == "(" and i in literal:
            out.append(r"\(")
            atom, quant = True, 0
            i += 1
            continue
        if c == "\\":
            t, i = _escape(p, i, False)
            out.append(t)
            atom, quant = t not in (r"\b(?=\w)", r"\b(?<=\w)"), 0
            continue
        if c == "[":
            t, i = _char_class(p, i, strict)
            out.append(t)
            atom, quant = True, 0
            continue
        if c == "(":
            atom, quant = False, 0
            if p.startswith("(?#", i):
                j = p.find(")", i)
                if j < 0:
                    if strict:
                        raise MatlabError(*_UNMATCHED_PAREN)
                    out.append(re.escape(p[i:]))
                    break
                i = j + 1
                continue
            if p.startswith("(?", i):
                j = i + 2
                while j < n and p[j] in _FLAGS:
                    j += 1
                if j > i + 2 and j < n and p[j] == ")":
                    # inline flags: up to the end of the alternative in the enclosing group
                    out.append(f"(?{p[i + 2 : j]}:")
                    groups[-1][1].append(out[-1])
                    i = j + 1
                    continue
                if p.startswith("(?<", i) and not p.startswith(("(?<=", "(?<!"), i):
                    m = _NAME.match(p, i + 3)
                    if m and m.end() < n and p[m.end()] == ">":
                        open_group(f"(?P<{m.group()}>")
                        i = m.end() + 1
                        continue
                if p.startswith(("(??", "(?@"), i):
                    raise MatlabError(
                        "MATLAB's dynamic expressions are not supported", "pk5d:regexp"
                    )
                for opener in ("(?:", "(?=", "(?!", "(?<=", "(?<!", "(?>", "(?("):
                    if p.startswith(opener, i):
                        open_group(opener)
                        i += len(opener)
                        break
                else:
                    open_group("(")
                    i += 1
                continue
            open_group("(")
            i += 1
            continue
        if c == ")":
            if len(groups) == 1:
                out.append(r"\)")  # a ')' without '(' is the character itself
            else:
                out.append(")" * (len(groups.pop()[1]) + 1))
            atom, quant = True, 0
            i += 1
            continue
        if c == "|":
            out.append(")" * len(groups[-1][1]) + "|")
            groups[-1][1] = []
            i += 1
            atom, quant = False, 0
            continue
        if c in "*+?":
            if quant == 1 and c in "?+":
                out.append(c)  # lazy or possessive
                quant = 2
            elif atom and quant == 0:
                out.append(c)
                quant = 1
            else:
                out.append("\\" + c)  # nothing to repeat: the character itself
                atom, quant = True, 0
            i += 1
            continue
        if c == "{":
            m = _COUNT.match(p, i)
            if m:
                if atom and quant == 0:
                    out.append(m.group())
                    quant = 1
                else:
                    out.append(re.escape(m.group()))
                    atom, quant = True, 0
                i = m.end()
                continue
            if strict and p.find("}", i) < 0:
                raise MatlabError(
                    "Unmatched '{'.  Use '\\{' to specify literal '{'.",
                    "MATLAB:REGEXP:ErrorUnmatchedCurly",
                )
            out.append(r"\{")
            i += 1
            atom, quant = True, 0
            continue
        if c == "$":
            out.append(r"\Z")
            i += 1
            atom, quant = False, 0
            continue
        if c == "^":
            out.append(r"\A")
            i += 1
            atom, quant = False, 0
            continue
        out.append(c)
        i += 1
        atom, quant = True, 0
    if len(groups) > 1:
        return "", [g[0] for g in groups[1:]]
    out.append(")" * len(groups[0][1]))
    return "".join(out), []


@functools.lru_cache(maxsize=4096)
def matlab_regex(pattern: str, ignore_case: bool = False, strict: bool = False) -> re.Pattern:
    """MATLAB's ``regexp(..., pattern)`` (``regexpi`` with ``ignore_case``) as a ``re`` pattern.

    MATLAB's defaults: ``.`` matches newlines, ``^`` and ``$`` are the start and end of the
    text. ``regexp`` reads an unmatched ``(``, ``[`` or ``{`` as the character itself;
    ``regexpPattern`` (strict) raises MatlabError with MATLAB's identifier for them.
    """
    flags = re.DOTALL | (re.IGNORECASE if ignore_case else 0)
    try:
        return re.compile(_translate(pattern, strict), flags)
    except re.error as err:
        raise MatlabError(f"invalid regular expression {pattern!r}: {err}", "pk5d:regexp") from err


def contains_text(texts: Sequence[str], pattern: str, ignore_case: bool = False) -> list[bool]:
    """``contains(texts, pattern)`` for a literal pattern."""
    if ignore_case:
        low = pattern.lower()
        return [low in t.lower() for t in texts]
    return [pattern in t for t in texts]


def matches_channel_pattern(
    texts: Sequence[str], pattern: str, ignore_case: bool = False
) -> list[bool]:
    """``contains(texts, p) | contains(texts, regexpPattern(p))``: the channel pattern test.

    The regular expression is always compiled (as regexpPattern does: strict), so a pattern
    that is not a valid MATLAB expression raises MATLAB's error even where the text contains
    it literally.
    """
    rx = matlab_regex(pattern, ignore_case, True)
    lit = contains_text(texts, pattern, ignore_case)
    return [a or rx.search(t) is not None for a, t in zip(lit, texts, strict=True)]


# ---------------------------------------------------------------------------------------------
# the tile name grammar
# ---------------------------------------------------------------------------------------------

# XR_parseSettingFiles_wrapper's expressions, as MATLAB writes them; translated by matlab_regex
FILENAME_EXPRESSION_CAM = (
    r"(?<prefix>\w*)Scan_Iter_(?<Iter>(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM1_stack"
    r"(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_"
    r"(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t"
)
FILENAME_EXPRESSION_NO_CAM = (
    r"(?<prefix>\w*)Scan_Iter_(?<Iter>(\d+_)*\d+?)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_"
    r"(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_"
    r"(?<z>-?\d+)z_(?<t>\d+)t"
)


def regexp_names(text: str, rx: re.Pattern) -> list[dict[str, str]]:
    """``regexp(text, expr, 'names')``: one dict per match, '' for a token that did not take
    part."""
    names = list(rx.groupindex)
    return [{k: m.group(k) or "" for k in names} for m in rx.finditer(text)]


def parse_filename_fields(fsnames: Sequence[str]) -> list[list[dict[str, str]]]:
    """The tile fields of file names without extension (XR_parseSettingFiles_wrapper).

    Names with ``_Cam<x>_ch`` (all of them) use the camera grammar, else names that all have
    ``_ch<digit>_`` the grammar without a camera; otherwise MATLAB fails (no expression).
    Each name gives its list of matches (``regexpi(..., 'names')``), each a dict of prefix,
    Iter, (Cam,) ch, stack, laser, abstime, fpgatime, x, y, z and t, as text.
    """
    if all(matlab_regex(r"_Cam\w_ch").search(s) for s in fsnames):
        expression = FILENAME_EXPRESSION_CAM
    elif all(matlab_regex(r"_ch[0-9]_").search(s) for s in fsnames):
        expression = FILENAME_EXPRESSION_NO_CAM
    else:
        raise MatlabError(
            "Unrecognized function or variable 'expression'.", "MATLAB:UndefinedFunction"
        )
    rx = matlab_regex(expression, True)
    return [regexp_names(s, rx) for s in fsnames]


# ---------------------------------------------------------------------------------------------
# parseImageFilenames
# ---------------------------------------------------------------------------------------------


@dataclass
class ImageFilenames:
    """Result of ``parse_image_filenames`` (parseImageFilenames's outputs)."""

    fnames: list[str]
    fsns: list[str]
    fd_inds: np.ndarray  # int64, 0-based index into data_paths
    filepaths: list[str]
    ch_inds: np.ndarray  # int64, 0-based index of the first channel pattern a file matches

    def __iter__(self):
        return iter((self.fnames, self.fsns, self.fd_inds, self.filepaths, self.ch_inds))


def _dir_names(pattern: str) -> list[str]:
    """The names ``dir(pattern)`` lists for ``<folder>/*<ext>``, in MATLAB's order (as
    ``matlab_dir``, without its per-entry checks of folders, which only the names need)."""
    if os.path.lexists(pattern):
        return [e.name for e in matlab_dir(pattern)]
    folder, name = os.path.split(pattern)
    try:
        names = os.listdir(folder or ".")
    except OSError:
        return []
    rx = re.compile(fnmatch.translate(name.replace("?", "[?]")), re.DOTALL)
    names = [n for n in names if rx.match(n)]
    return sorted(names, key=lambda s: s.encode("utf-16-be", errors="surrogatepass"))


def parse_image_filenames(
    data_paths: str | Sequence[str],
    zarr_file: bool,
    channel_patterns: Sequence[str],
    data_formats: Sequence[str] = (),
) -> ImageFilenames:
    """The image files of dataset folders that match channel patterns (parseImageFilenames).

    Lists ``*.tif`` (``*.zarr`` with zarr_file; every name ending with one of data_formats,
    any case, when they are given and zarr_file is false) in each folder, in MATLAB's ``dir``
    order, keeps names that contain a channel pattern as text or as a regular expression
    (any case), and gives each file its folder index and its first matching pattern.
    """
    if isinstance(data_paths, (str, os.PathLike)):
        data_paths = [os.fspath(data_paths)]
    data_paths = list(data_paths)
    if not isinstance(zarr_file, (bool, np.bool_)):
        raise MatlabError(
            "The value of 'zarrFile' is invalid. It must satisfy the function: islogical.",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )
    channel_patterns = list(channel_patterns)
    data_formats = list(data_formats)
    print("Parse image filenames for data path(s):\n    " + "\n    ".join(data_paths))

    nc = len(channel_patterns)
    fnames: list[str] = []
    fd_inds: list[int] = []
    ch_inds: list[int] = []
    for d, data_path in enumerate(data_paths):
        if zarr_file:
            names = _dir_names(data_path + "/*.zarr")
        elif not data_formats:
            names = _dir_names(data_path + "/*.tif")
        else:
            names = _dir_names(data_path + "/*")
            low = [n.lower() for n in names]
            names = [
                n
                for n, ln in zip(names, low, strict=True)
                if any(ln.endswith(f.lower()) for f in data_formats)
            ]
        flags = [matches_channel_pattern(names, p, True) for p in channel_patterns]
        for k, name in enumerate(names):
            hits = [c for c in range(nc) if flags[c][k]]
            if hits:
                fnames.append(name)
                fd_inds.append(d)
                ch_inds.append(hits[0])

    fsns = [fileparts(f)[1] for f in fnames]
    if not fnames:
        if zarr_file:
            _warn(
                "No Zarr format (*.zarr) image files found in the input data folder(s)! Please "
                "check if dataPaths are correct!"
            )
        else:
            _warn(
                "No Tiff format (*.tif) image files found in the input data folder(s)! Please "
                "check if dataPaths are correct!"
            )
    filepaths = [f"{data_paths[d]}/{f}" for d, f in zip(fd_inds, fnames, strict=True)]
    if not fnames:
        _warn(
            "The input data folder(s) is empty, or no file matches the channel pattern, or file "
            "type (tiff/zarr) is not correctly setup. Please check."
        )
    return ImageFilenames(
        fnames, fsns, np.array(fd_inds, np.int64), filepaths, np.array(ch_inds, np.int64)
    )


# ---------------------------------------------------------------------------------------------
# groupImageFilenamesByChannels
# ---------------------------------------------------------------------------------------------


def _code_units(names: list[str], reverse: bool, pad: int) -> np.ndarray:
    """names as rows of UTF-16 code units (MATLAB's chars), padded with ``pad`` up to one
    column past the longest name"""
    units = [np.frombuffer(s.encode("utf-16-le", errors="surrogatepass"), "<u2") for s in names]
    width = max((len(u) for u in units), default=0) + 1
    out = np.full((len(units), width), pad, np.int32)
    for k, u in enumerate(units):
        out[k, : len(u)] = u[::-1] if reverse else u
    return out


def group_image_filenames_by_channels(
    filenames: Sequence[str],
    channel_patterns: Sequence[str],
    ch_inds: Sequence[int] | np.ndarray | None = None,
) -> np.ndarray:
    """Group files of several channels into rows of one file per channel.

    ch_inds: the 0-based channel of each file (default: the first channel pattern each name
    matches, any case). Each file not yet grouped starts a row; for every other channel the
    row takes the free file of that channel with the longest common start plus end (the first
    such file). Returns an int64 array (groups x channels) of 0-based file indices.

    As MATLAB, it fails when a channel has no free file left for a row (channels with
    different numbers of files).
    """
    filenames = list(filenames)
    nf = len(filenames)
    nc = len(channel_patterns)
    if ch_inds is None:
        flags = [matches_channel_pattern(filenames, p, True) for p in channel_patterns]
        chs: list[int] = []
        for f in range(nf):
            hits = [c for c in range(nc) if flags[c][f]]
            if not hits:
                raise MatlabError(
                    f"Non-scalar in Uniform output, at index {f + 1}, output 1.\n"
                    "Set 'UniformOutput' to false.",
                    "MATLAB:arrayfun:NotAScalarOutput",
                )
            chs.append(hits[0] + 1)
    else:
        chs = [int(c) + 1 for c in np.asarray(ch_inds).ravel()]
    if len(set(chs)) != nc:
        raise MatlabError(
            "Number of channels of the images does not match the number of channel patterns!"
        )

    # MATLAB's 1-based logic; group_f_inds keeps its values from row to row. The common
    # start and end of two names (extractCommonStringFromStart of the names and of the names
    # flipped) are counted on rows of code units: padded with -1, and the name compared with
    # -2, so the first difference is never past the shorter name.
    included = np.zeros(nf, bool)
    ch_arr = np.array(chs, np.int64)
    fwd = _code_units(filenames, False, -1)
    rev = _code_units(filenames, True, -1)
    lens = np.array([len(s.encode("utf-16-le", errors="surrogatepass")) // 2 for s in filenames])
    group = [0] * nc
    rows: list[list[int]] = []
    for f in range(1, nf + 1):
        if included[f - 1]:
            continue
        rf, rr = fwd[f - 1].copy(), rev[f - 1].copy()
        rf[lens[f - 1] :] = -2
        rr[lens[f - 1] :] = -2
        for c in range(1, nc + 1):
            if chs[f - 1] == c:
                group[c - 1] = f
                continue
            cur = np.flatnonzero(~included & (ch_arr == c))
            if cur.size == 0:
                # max of an empty column gives an empty 0x1 index: the assignment fails
                raise MatlabError(
                    "Unable to perform assignment because the left and right sides have a "
                    "different number of elements.",
                    "MATLAB:matrix:singleSubscriptNumelMismatch",
                )
            counts = np.argmin(fwd[cur] == rf, axis=1) + np.argmin(rev[cur] == rr, axis=1)
            group[c - 1] = int(cur[np.argmax(counts)]) + 1  # the first largest, as max
        included[np.array(group) - 1] = True
        rows.append(list(group))
    if not rows:
        return np.zeros((0, 0), np.int64)
    return np.array(rows, np.int64) - 1


# ---------------------------------------------------------------------------------------------
# XR_parseImageFilenames
# ---------------------------------------------------------------------------------------------


@dataclass
class DatasetFilenames:
    """Result of ``parse_dataset_filenames`` (XR_parseImageFilenames's nine outputs)."""

    fnames: list[str]
    fd_inds: np.ndarray  # int64, 0-based index into data_paths
    gfnames: list[list[str]]  # each file with its partial volume files
    partialvols: np.ndarray  # bool
    data_sizes: np.ndarray  # float64, bytes of the decoded volume
    flip_zstack: np.ndarray  # bool
    latest_modify_times: np.ndarray  # float64, minutes, per data path
    ftp_inds: np.ndarray  # int64, 0-based index of each data path's first file; -1: none
    mask_fullpaths: list[str | None] = field(default_factory=list)

    def __iter__(self):
        return iter(
            (
                self.fnames,
                self.fd_inds,
                self.gfnames,
                self.partialvols,
                self.data_sizes,
                self.flip_zstack,
                self.latest_modify_times,
                self.ftp_inds,
                self.mask_fullpaths,
            )
        )


_DTYPE_BYTES = {"uint8": 1, "uint16": 2, "single": 4}


def _group_partial_volume_files(
    data_path: str, ext: str, channel_patterns: Sequence[str]
) -> tuple[bool, list[list[str]], list, list]:
    """groupPartialVolumeFiles(dataPath, 'ext', ext, 'channelPatterns', ...).

    A stand-in until T05's petakit5d.tools.partial_volumes lands (same results, the folder
    form only). Returns (contains partial volumes, groups of names, mtimes, sizes): without
    partial volumes mtimes and sizes are one value per file, with them one list per group.
    """
    from petakit5d.io.imagesize import get_image_dtype, get_image_size

    print("Check partial volume files... ", end="")
    entries = matlab_dir(data_path + os.sep + "*" + ext)
    fnames = [e.name for e in entries]
    mtimes = []
    sizes = []
    for e in entries:
        full = os.path.join(e.folder, e.name)
        try:
            st = os.stat(full)
            mtimes.append(st.st_mtime)
            sizes.append(0.0 if e.isdir else float(st.st_size))
        except OSError:
            mtimes.append(math.nan)
            sizes.append(0.0)
    if channel_patterns:
        fns = [data_path + os.sep + f for f in fnames]
        keep = [False] * len(fnames)
        for p in channel_patterns:
            keep = [a or b for a, b in zip(keep, matches_channel_pattern(fns, p), strict=True)]
        fnames = [x for x, k in zip(fnames, keep, strict=True) if k]
        mtimes = [x for x, k in zip(mtimes, keep, strict=True) if k]
        sizes = [x for x, k in zip(sizes, keep, strict=True) if k]
    if not fnames:
        _warn("The input image list is empty or none of the images matches channelPatterns.")
        return False, [], [], []

    dtype = get_image_dtype(data_path + os.sep + fnames[0], ext == ".zarr")
    dbytes = _DTYPE_BYTES.get(dtype, 8)
    if ext in (".tif", ".tiff"):
        for f in sorted({0, len(fnames) - 1}):
            sizes[f] = float(math.prod(get_image_size(data_path + os.sep + fnames[f]))) * dbytes
        sizes = [max(sizes)] * len(sizes)
    elif ext == ".zarr":
        for f, name in enumerate(fnames):
            try:
                sizes[f] = float(math.prod(get_image_size(data_path + os.sep + name))) * dbytes
            except Exception as err:  # MATLAB displays the error and goes on
                print(err)

    full = [data_path + os.sep + f for f in fnames]
    part_rx = matlab_regex("_part\\d+" + ext + "$", True)
    contain = [part_rx.search(f) is not None for f in full]
    if not any(contain):
        print("\nNo partial volume files found.")
        return False, [[f] for f in fnames], mtimes, sizes

    fsnames = [fileparts(f)[1] for f in full]
    main_rx = matlab_regex("_part\\d+$", True)
    mains = [i for i, s in enumerate(fsnames) if main_rx.search(s) is None]
    groups, gtimes, gsizes = [], [], []
    for i in mains:
        rx = matlab_regex(fsnames[i] + "_part\\d+$", True)
        parts = [k for k, s in enumerate(fsnames) if rx.search(s)]
        names = [fsnames[i]] + sorted((fsnames[k] for k in parts), key=matlab_sort_key)
        groups.append([s + ext for s in names])
        gtimes.append([mtimes[i]] + [mtimes[k] for k in parts])
        gsizes.append([sizes[i]] + [sizes[k] for k in parts])
    print("Done!")
    return True, groups, gtimes, gsizes


def parse_dataset_filenames(
    data_paths: Sequence[str],
    channel_patterns: Sequence[str],
    parse_setting_file: bool,
    flip_zstack: bool,
    decon: bool,
    decon_paths: Sequence[str],
    streaming: bool,
    min_modify_time: float,
    zarr_file: bool = False,
) -> DatasetFilenames:
    """The files the pipeline processes (XR_parseImageFilenames).

    For each folder: the ``.tif`` (``.zarr``) files whose path matches a channel pattern, with
    their partial volume files (``_part0001``, ...), the size of the decoded volume and the
    folder index. In streaming mode, files changed in the last ``min_modify_time`` minutes
    (times the median number of files per group) are left for later. data_paths end with a
    separator (names are appended as they are). With parse_setting_file, a file's z stack is
    flipped where its setting file gives a negative stage interval. Each folder's first file
    of the first matching channel pattern is its first time point (ftp_inds), with its eroded
    mask path in the decon folder when decon is set.
    """
    ext = ".zarr" if zarr_file else ".tif"
    data_paths = list(data_paths)
    channel_patterns = list(channel_patterns)
    nd = len(data_paths)
    t0 = time.perf_counter()
    print("Parse image filenames for data path(s):\n    " + "\n    ".join(data_paths))

    fnames_cell: list[list[str]] = [[] for _ in range(nd)]
    gf_cell: list[list[list[str]]] = [[] for _ in range(nd)]
    pv_cell: list[list[bool]] = [[] for _ in range(nd)]
    ds_cell: list[list[float]] = [[] for _ in range(nd)]
    latest = np.zeros(nd)
    # partialvols is logical if every folder had files: cat(1, ...) with the [] of a skipped
    # folder makes it double
    skipped = False
    for d, data_path in enumerate(data_paths):
        contain, groups, gtimes, gsizes = _group_partial_volume_files(
            data_path, ext, channel_patterns
        )
        fn_d = [g[0] for g in groups]
        times_d = [max(t) for t in gtimes] if contain else list(gtimes)
        sizes_d = list(gsizes)
        if not fn_d:
            skipped = True
            continue
        if streaming:
            now = time.time()
            ages = [(now - t) / 60 for t in times_d]
            lt = max(sorted(ages)[: len(channel_patterns)])
            latest[d] = lt
            ngf = statistics.median(len(g) for g in groups)
            if lt < min_modify_time * ngf:
                lt = max(lt, min_modify_time)
                keep = [a > lt for a in ages]
                fn_d = [x for x, k in zip(fn_d, keep, strict=True) if k]
                groups = [x for x, k in zip(groups, keep, strict=True) if k]
                sizes_d = [x for x, k in zip(sizes_d, keep, strict=True) if k]
        fnames_cell[d] = fn_d
        gf_cell[d] = groups
        if contain:
            pv_cell[d] = [len(g) > 1 for g in groups]
            ds_cell[d] = [float(sum(s)) for s in sizes_d]
        else:
            pv_cell[d] = [False] * len(fn_d)
            ds_cell[d] = [float(s) for s in sizes_d]

    fd_inds = [d for d in range(nd) for _ in fnames_cell[d]]
    fnames = [f for fl in fnames_cell for f in fl]
    gfnames = [g for gl in gf_cell for g in gl]
    partialvols = [p for pl in pv_cell for p in pl]
    data_sizes = [s for sl in ds_cell for s in sl]

    if not fnames:
        if zarr_file:
            _warn(
                "No Zarr format (*.zarr) image files found in the input data folder(s)! Please "
                "check if dataPaths are correct!"
            )
        else:
            _warn(
                "No Tiff format (*.tif) image files found in the input data folder(s)! Please "
                "check if dataPaths are correct!"
            )
        return DatasetFilenames(
            [],
            np.array(fd_inds, np.int64),
            [],
            np.array(partialvols, np.float64 if skipped else bool),
            np.array(data_sizes, np.float64),
            np.zeros(0, bool),
            latest,
            np.zeros(0, np.int64),
            [],
        )

    fullnames = [data_paths[d] + f for d, f in zip(fd_inds, fnames, strict=True)]
    include = [False] * len(fnames)
    for p in channel_patterns:
        include = [
            a or b for a, b in zip(include, matches_channel_pattern(fullnames, p), strict=True)
        ]
    pv_dtype = np.float64 if skipped else bool
    fnames = [x for x, k in zip(fnames, include, strict=True) if k]
    fd_inds = [x for x, k in zip(fd_inds, include, strict=True) if k]
    gfnames = [x for x, k in zip(gfnames, include, strict=True) if k]
    partialvols = [x for x, k in zip(partialvols, include, strict=True) if k]
    data_sizes = [x for x, k in zip(data_sizes, include, strict=True) if k]
    nf = len(fnames)

    flip = np.full(nf, bool(flip_zstack))
    if parse_setting_file:
        from petakit5d.pipeline.settings import parse_image_setting_files

        frame_fullpaths = [data_paths[fd_inds[i]] + fnames[i] for i in range(nf)]
        info = parse_image_setting_files(frame_fullpaths)
        if not info:
            # struct() without fields: no images, or only partial volume files
            raise MatlabError('Unrecognized field name "StageInterval".', "MATLAB:nonExistentField")
        # [settingInfo.StageInterval] skips the [] of setting files that were missing
        vals = [s["StageInterval"] for s in info if s.get("StageInterval") is not None]
        flip = np.array([v < 0 for v in vals], bool)

    ftp_inds = np.full(nd, -1, np.int64)
    mask_fullpaths: list[str | None] = [None] * nd
    for d in range(nd):
        if not fnames_cell[d]:
            continue
        ftp_name = ""
        full_d = [data_paths[d] + x for x in fnames_cell[d]]
        for p in channel_patterns:
            hits = matches_channel_pattern(full_d, p)
            if any(hits):
                ftp_name = fnames_cell[d][hits.index(True)]
                break
        if ftp_name:
            ind = [i for i in range(nf) if fnames[i] == ftp_name and fd_inds[i] == d]
            ftp_inds[d] = ind[0]
            if decon:
                mask_fullpaths[d] = f"{decon_paths[d]}/Masks/{fileparts(ftp_name)[1]}_eroded.zarr"
    print(f"Done! Elapsed time is {time.perf_counter() - t0:.6f} seconds.\n")

    if not streaming and nf == 0:
        _warn(
            "The input data folder(s) is empty, or no file matches the channel pattern, or file "
            "type (tiff/zarr) is not correctly setup. Please check."
        )
    return DatasetFilenames(
        fnames,
        np.array(fd_inds, np.int64),
        gfnames,
        np.array(partialvols, pv_dtype),
        np.array(data_sizes, np.float64),
        flip,
        latest,
        ftp_inds,
        mask_fullpaths,
    )
