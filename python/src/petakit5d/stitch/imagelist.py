"""Image lists of stitching: tables, reading them as MATLAB does, and parsing the tile names.

Ports ``stitch_parse_image_list_information``, ``stitch_parse_multi_loc_image_list_information``,
``parseTileInfo`` and ``stitch_process_filenames``, with the pieces of MATLAB they rely on:

- ``Table``: the subset of a MATLAB table the image lists use: named columns of equal height,
  in order. A column is a float64 array (double), an object array of ``str`` (a cell of char
  or a string array; ``None`` is a missing string) or a ``'<U1'`` array (a char column, such
  as ``camera``). Rows are selected with a mask, indices or a slice (``tab[mask]``), columns
  by name (``tab["ch"]``).
- ``read_table``: ``readtable(file, 'Delimiter', 'comma')``: header names made valid as MATLAB
  does (``StageX (um)`` is ``StageX_um_``, the original in ``descriptions``), the type of each
  column by its majority (numbers at least as many as texts: double, text in it NaN), empty
  cells (NaN, ``''``), quoted fields, leading lines that hold no number skipped, a first line
  with a number read as data (names ``Var1``, ...).
- ``write_table``: ``writetable(t, file, 'Delimiter', ',')`` as R2024b writes it: the texts of
  a column quoted when one needs it, numbers with 15 significant digits (``%.15g``), missing
  strings empty.
- ``sortrows``, ``unique`` of text, and ``contains`` of a channel pattern (plain text or a
  ``regexpPattern``).

``parse_image_list`` returns MATLAB's twelve outputs as an ``ImageListInfo`` (a named tuple, so
``tab, primary_tab, *rest = parse_image_list(...)`` works as in MATLAB). Quirks kept from MATLAB:
the tile name grammar puts the last digit of an iteration without sub-iteration into
``subIter`` (``Iter_0012_Cam`` gives ``Iter`` 1, ``subIter`` ``'2'``, ``fullIter`` ``'0012'``);
the file names used for the primary-channel checks and the mapped-name filter are not filtered
by the matched rows (``fsn(tinds)`` of the unfiltered names); a tile exists if its name without
extension is part of a file name of the folder; ``Cam`` is a char column (a ``str`` here, one
letter per camera). With several folders, ``prefix`` is the unique characters of the folders'
prefixes, as ``unique`` of a char matrix gives.
"""

from __future__ import annotations

import os
import re
import warnings
from collections.abc import Iterable, Sequence
from typing import Any, NamedTuple

import numpy as np

from petakit5d.core.dtypes import NAN
from petakit5d.core.paths import fileparts, matlab_dir, matlab_sort_key
from petakit5d.core.text import sprintf as _core_sprintf

__all__ = [
    "ImageListError",
    "ImageListInfo",
    "MatlabWarning",
    "Table",
    "contains_pattern",
    "make_valid_name",
    "matlab_regex",
    "parse_image_list",
    "parse_multi_loc_image_list",
    "parse_tile_info",
    "process_filenames",
    "read_table",
    "sortrows",
    "unique_text",
    "write_table",
]

# the tile name grammar of stitch_parse_image_list_information (MATLAB syntax)
TILE_EXPRESSION_CAM = (
    r"(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)"
    r"_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_"
    r"(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)"
)
TILE_EXPRESSION_NOCAM = (
    r"(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_ch(?<ch>\d+)"
    r"_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_"
    r"(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)"
)
_SORT_KEYS = ("camera", "ch", "z", "y", "x", "t", "fullIter")
XCORR_MODES = ("primaryFirst", "primary", "all", "stitchInfo")


class ImageListError(ValueError):
    """An error MATLAB raises; ``identifier`` is MATLAB's error id (``''`` for none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


class MatlabWarning(UserWarning):
    """A warning MATLAB gives; ``identifier`` is MATLAB's warning id (``''`` for none)."""

    def __init__(self, message: str, identifier: str = "") -> None:
        super().__init__(message)
        self.identifier = identifier


def _warn(message: str, identifier: str = "") -> None:
    warnings.warn(MatlabWarning(message, identifier), stacklevel=3)


def _sprintf(fmt: str, *args) -> str:
    try:
        return _core_sprintf(fmt, *args)
    except IndexError:
        # core.text.sprintf fails when a vector argument makes the format recycle over an
        # empty char argument (reported to F02); MATLAB prints the empty argument as nothing
        raise NotImplementedError(
            f"sprintf({fmt!r}, ...) recycles its format over an empty argument"
        ) from None


# ------------------------------------------------------------------------------------------
# tables
# ------------------------------------------------------------------------------------------


def _column(v: Any) -> np.ndarray:
    """A table column: float64, object (text) or '<U1' (char), always 1-D."""
    if isinstance(v, np.ndarray):
        if v.dtype == object:
            return v.reshape(-1)
        if v.dtype.kind == "U":
            if v.dtype.itemsize // 4 > 1:
                return np.array([str(s) for s in v.reshape(-1)], dtype=object)
            return v.reshape(-1).astype("<U1")
        if v.dtype.kind in "biuf":
            return v.reshape(-1).astype(np.float64)
        raise TypeError(f"unsupported column type {v.dtype}")
    if isinstance(v, str):
        return np.array(list(v), dtype="<U1")
    items = list(v)
    if items and all(isinstance(s, (str, type(None))) for s in items):
        out = np.empty(len(items), dtype=object)
        out[:] = items
        return out
    if not items:
        return np.empty(0, dtype=np.float64)
    return np.asarray(items, dtype=np.float64).reshape(-1)


def column_kind(col: np.ndarray) -> str:
    """``'double'``, ``'text'`` (a cell of char or strings) or ``'char'``."""
    if col.dtype == object:
        return "text"
    if col.dtype.kind == "U":
        return "char"
    return "double"


class Table:
    """A MATLAB table with named columns (see the module docstring)."""

    __slots__ = ("_cols", "_height", "descriptions")

    def __init__(
        self,
        columns: dict[str, Any] | Iterable[tuple[str, Any]] | None = None,
        descriptions: Sequence[str] = (),
        height: int | None = None,
    ) -> None:
        self._cols: dict[str, np.ndarray] = {}
        self._height = height
        self.descriptions: list[str] = []
        items = columns.items() if isinstance(columns, dict) else (columns or ())
        for k, v in items:
            self[k] = v
        if self._height is None:
            self._height = 0
        # VariableDescriptions: empty, or one per variable (the original headers)
        self.descriptions = list(descriptions)

    @property
    def names(self) -> list[str]:
        """The variable names, in order."""
        return list(self._cols)

    @property
    def height(self) -> int:
        return int(self._height or 0)

    def __len__(self) -> int:
        return self.height

    def __contains__(self, name: object) -> bool:
        return name in self._cols

    def __getitem__(self, key):
        if isinstance(key, str):
            try:
                return self._cols[key]
            except KeyError:
                raise ImageListError(
                    f"Unrecognized table variable name '{key}'.",
                    "MATLAB:table:UnrecognizedVarName",
                ) from None
        return self.rows(key)

    def __setitem__(self, name: str, value: Any) -> None:
        col = _column(value)
        if self._height is None:
            self._height = len(col)
        elif len(col) != self._height:
            raise ImageListError(
                f"column {name!r} has {len(col)} rows, the table {self._height}",
                "MATLAB:table:RowDimensionMismatch",
            )
        if name not in self._cols and getattr(self, "descriptions", None):
            self.descriptions.append("")  # a new variable has an empty description
        self._cols[name] = col

    def __delitem__(self, name: str) -> None:
        del self._cols[name]

    def rows(self, index) -> Table:
        """The rows of a mask, an index array (0-based) or a slice, as a new table."""
        if isinstance(index, slice):
            sel = np.arange(self.height)[index]
        else:
            index = np.asarray(index)
            sel = np.flatnonzero(index) if index.dtype == bool else index.astype(np.intp)
            if index.dtype == bool and len(index) != self.height:
                raise IndexError(f"mask of {len(index)} rows for a table of {self.height}")
        out = Table(height=len(sel), descriptions=self.descriptions)
        for k, c in self._cols.items():
            out._cols[k] = c[sel]
        return out

    def copy(self) -> Table:
        out = Table(height=self.height, descriptions=self.descriptions)
        for k, c in self._cols.items():
            out._cols[k] = c.copy()
        return out

    def kind(self, name: str) -> str:
        return column_kind(self[name])

    @staticmethod
    def vcat(tables: Sequence[Table]) -> Table:
        """``cat(1, tables{:})``: the same variables (matched by name), in the first's order."""
        tables = list(tables)
        if not tables:
            return Table()
        names = tables[0].names
        for t in tables[1:]:
            if sorted(t.names) != sorted(names):
                raise ImageListError(
                    "All tables being vertically concatenated must have the same variable names.",
                    "MATLAB:table:vertcat:UnequalVarNames",
                )
        out = Table(height=sum(t.height for t in tables), descriptions=tables[0].descriptions)
        for k in names:
            cols = [t[k] for t in tables]
            kinds = {column_kind(c) for c in cols}
            if len(kinds) > 1:
                raise ImageListError(
                    f"Cannot concatenate the table variable '{k}' because it is of "
                    "different types in the tables.",
                    "MATLAB:table:vertcat:VertcatMethodFailed",
                )
            out._cols[k] = np.concatenate(cols)
        return out

    def __repr__(self) -> str:
        return f"Table({self.height} rows: {', '.join(self.names)})"


# ------------------------------------------------------------------------------------------
# readtable
# ------------------------------------------------------------------------------------------

_WS = " \t\b"
_NUM_RE = re.compile(r"[+-]?(?:(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?|[iI][nN][fF]|[nN][aA][nN])\Z")
_DATE_RE = re.compile(
    r"(?:\d{4}-\d{1,2}-\d{1,2}|\d{1,2}/\d{1,2}/\d{2,4}|\d{1,2}-[A-Za-z]{3}-\d{2,4})"
    r"(?:[ T]\d{1,2}:\d{2}(?::\d{2}(?:\.\d+)?)?)?\Z"
)
_DURATION_RE = re.compile(r"-?\d+:\d{2}(?::\d{2}(?:\.\d+)?)?\Z")
_KEYWORDS = frozenset(
    "break case catch classdef continue else elseif end for function global if otherwise "
    "parfor persistent return spmd switch try while".split()
)
NAMELENGTHMAX = 63
_MODIFIED_MSG = (
    "Column headers from the file were modified to make them valid MATLAB identifiers before "
    "creating variable names for the table. The original column headers are saved in the "
    "VariableDescriptions property.\nSet 'VariableNamingRule' to 'preserve' to use the "
    "original column headers as table variable names."
)


def _parse_number(s: str) -> float | None:
    """A field as readtable reads a number (``None``: not a number)."""
    if not _NUM_RE.match(s):
        return None
    low = s.lower()
    if low.lstrip("+-") == "nan":
        return NAN
    return float(low.replace("d", "e"))


def make_valid_name(s: str) -> str:
    """``matlab.lang.makeValidName``: whitespace removed (a lowercase letter after it made
    uppercase), other characters that are not ``[A-Za-z0-9_]`` replaced by ``_``, ``x`` put
    before a name that does not start with a letter and before a keyword (whose first letter
    is made uppercase), at most ``namelengthmax`` characters."""
    out = []
    up = False
    for ch in s:
        if ch.isspace():
            up = bool(out)
            continue
        if up and "a" <= ch <= "z":
            ch = ch.upper()
        up = False
        out.append(ch if (ch.isascii() and (ch.isalnum() or ch == "_")) else "_")
    name = "".join(out)
    if name in _KEYWORDS:
        name = "x" + name[0].upper() + name[1:]
    elif not name or not ("a" <= name[0].lower() <= "z"):
        name = "x" + name
    return name[:NAMELENGTHMAX]


def _make_unique(names: list[str]) -> list[str]:
    """``matlab.lang.makeUniqueStrings``: later duplicates get ``_1``, ``_2``, ... (skipping
    names that exist)."""
    taken = set(names)
    seen: set[str] = set()
    out = []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
            continue
        k = 1
        while True:
            suffix = f"_{k}"
            cand = n[: NAMELENGTHMAX - len(suffix)] + suffix
            if cand not in taken:
                break
            k += 1
        taken.add(cand)
        seen.add(cand)
        out.append(cand)
    return out


def _decode(raw: bytes) -> str:
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]  # the UTF-8 byte order mark
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


_TOKEN = re.compile(r'"((?:[^"]|"")*)"|([^,"\r\n]+)|(,)|(\r\n|\n|\r)|(")')


def _field(parts: list[tuple[str, bool]]) -> str:
    """A field from its quoted and unquoted parts: whitespace trimmed outside the quotes."""
    if not parts:
        return ""
    parts = list(parts)
    text, quoted = parts[0]
    if not quoted:
        parts[0] = (text.lstrip(_WS), False)
    text, quoted = parts[-1]
    if not quoted:
        parts[-1] = (text.rstrip(_WS), False)
    return "".join(t for t, _ in parts)


def _records(text: str) -> list[list[str]]:
    """The records of a csv text: fields split at commas outside double quotes (``""`` is a
    quote inside them, line breaks inside them are kept); empty lines are skipped."""
    rows: list[list[str]] = []
    fields: list[str] = []
    parts: list[tuple[str, bool]] = []
    started = False
    for m in _TOKEN.finditer(text):
        q, plain, comma, nl, lone = m.groups()
        if q is not None:
            parts.append((q.replace('""', '"'), True))
            started = True
        elif plain is not None:
            parts.append((plain, False))
            started = True
        elif lone is not None:
            parts.append(('"', False))
            started = True
        elif comma is not None:
            fields.append(_field(parts))
            parts = []
            started = True
        else:
            if started:
                fields.append(_field(parts))
                if not (len(fields) == 1 and fields[0] == ""):
                    rows.append(fields)
            fields, parts, started = [], [], False
    if started:
        fields.append(_field(parts))
        if not (len(fields) == 1 and fields[0] == ""):
            rows.append(fields)
    return rows


def read_table(filename: str | os.PathLike, delimiter: str = ",") -> Table:
    """``readtable(filename, 'Delimiter', 'comma')`` (see the module docstring).

    Not supported (``NotImplementedError``): columns MATLAB would read as datetime or
    duration; rows with different numbers of fields. Complex numbers read as text here.
    """
    if delimiter not in (",", "comma"):
        raise NotImplementedError("only comma-delimited files are supported")
    with open(filename, "rb") as f:
        rows = _records(_decode(f.read()))
    if not rows:
        return Table()
    ncol = len(rows[0])
    if any(len(r) != ncol for r in rows):
        raise NotImplementedError(
            f"{os.fspath(filename)}: rows with different numbers of fields are not supported"
        )
    isnum = [[_parse_number(f) is not None for f in r] for r in rows]
    # data start at the first line with a number; a first line with one is data (no header)
    d0 = next((i for i, r in enumerate(isnum) if any(r)), None)
    if d0 is None:
        d0 = 1
    header = rows[0] if d0 >= 1 else None
    data = rows[d0:]
    dnum = isnum[d0:]
    # a trailing column without header and data is not a variable
    while ncol > 0 and header is not None and header[ncol - 1] == "":
        if any(r[ncol - 1] != "" for r in data):
            break
        ncol -= 1
    names: list[str]
    descriptions: list[str] = []
    if header is None:
        names = [f"Var{j + 1}" for j in range(ncol)]
    else:
        orig = [h if h != "" else f"Var{j + 1}" for j, h in enumerate(header[:ncol])]
        valid = [make_valid_name(h) if header[j] != "" else orig[j] for j, h in enumerate(orig)]
        names = _make_unique(valid)
        if names != orig:
            descriptions = orig
            _warn(_MODIFIED_MSG, "MATLAB:table:ModifiedAndSavedVarnames")
    tab = Table(height=len(data))
    for j in range(ncol):
        vals = [r[j] for r in data]
        n_num = sum(1 for r in dnum if r[j])
        n_text = sum(1 for v, r in zip(vals, dnum, strict=True) if v != "" and not r[j])
        if n_num >= n_text:
            col = np.array(
                [(_parse_number(v) if r[j] else NAN) for v, r in zip(vals, dnum, strict=True)],
                dtype=np.float64,
            )
        else:
            texts = [v for v in vals if v != ""]
            if texts and (
                all(_DATE_RE.match(v) for v in texts) or all(_DURATION_RE.match(v) for v in texts)
            ):
                raise NotImplementedError(
                    f"readtable reads column {names[j]!r} as datetime or duration; not supported"
                )
            col = np.empty(len(vals), dtype=object)
            col[:] = vals
        tab[names[j]] = col
    tab.descriptions = descriptions
    return tab


# ------------------------------------------------------------------------------------------
# writetable
# ------------------------------------------------------------------------------------------


def format_number(v: float) -> str:
    """A double as writetable writes it: ``%.15g``, ``NaN``, ``Inf``, ``-Inf``."""
    v = float(v)
    if v != v:
        return "NaN"
    if v in (float("inf"), float("-inf")):
        return "Inf" if v > 0 else "-Inf"
    return f"{v:.15g}"


def _needs_quotes(col: np.ndarray) -> bool:
    # writetable quotes every text of a column once one holds a comma, a quote or a newline
    return any(s is not None and any(c in s for c in ',"\n\r') for s in col)


def _text_field(s: str | None, quote: bool) -> str:
    if s is None:
        return ""
    return '"' + s.replace('"', '""') + '"' if quote else s


def write_table(tab: Table, filename: str | os.PathLike) -> None:
    """``writetable(t, filename, 'Delimiter', ',')`` as MATLAB R2024b writes a csv file: the
    texts of a column are quoted (quotes doubled) when one of them holds a comma, a quote or
    a line break, missing strings are empty."""
    kinds = [tab.kind(n) for n in tab.names]
    cols = [tab[n] for n in tab.names]
    quote = [k == "text" and _needs_quotes(c) for k, c in zip(kinds, cols, strict=True)]
    lines = [",".join(tab.names)]
    for i in range(tab.height):
        out = []
        for k, c, q in zip(kinds, cols, quote, strict=True):
            if k == "double":
                out.append(format_number(c[i]))
            elif k == "text":
                out.append(_text_field(c[i], q))
            else:
                out.append(str(c[i]))
        lines.append(",".join(out))
    with open(filename, "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(lines) + "\n")


# ------------------------------------------------------------------------------------------
# MATLAB's sort, unique, regexp and contains
# ------------------------------------------------------------------------------------------


def _num_key(v: float) -> tuple[int, float]:
    return (1, 0.0) if v != v else (0, v)


def _col_keys(col: np.ndarray) -> list:
    k = column_kind(col)
    if k == "double":
        return [_num_key(v) for v in col.tolist()]
    if k == "char":
        return [ord(c) if c else 0 for c in col.tolist()]
    return [matlab_sort_key(s or "") for s in col.tolist()]


def sortrows(tab: Table, names: Sequence[str]) -> tuple[Table, np.ndarray]:
    """``[tab, idx] = sortrows(tab, names)``: stable, ascending, NaN last; ``idx`` 0-based."""
    keys = [_col_keys(tab[n]) for n in names]
    order = sorted(range(tab.height), key=lambda i: tuple(k[i] for k in keys))
    idx = np.asarray(order, dtype=np.intp)
    return tab.rows(idx), idx


def unique_text(values: Iterable[str]) -> list[str]:
    """``unique`` of a cell of char: sorted (UTF-16 code units), each once."""
    return sorted(set(values), key=matlab_sort_key)


def unique_num(values) -> np.ndarray:
    """``unique`` of doubles: sorted, each once, every NaN kept (at the end)."""
    a = np.asarray(values, dtype=np.float64).reshape(-1)
    nan = a[np.isnan(a)]
    b = np.sort(a[~np.isnan(a)], kind="stable")
    if b.size:
        b = b[np.concatenate(([True], b[1:] != b[:-1]))]
    return np.concatenate((b, nan))


_NAMED_GROUP = re.compile(r"\(\?<(?=[A-Za-z_])")


def matlab_regex(pattern: str, ignore_case: bool = False) -> re.Pattern:
    """A MATLAB regular expression as a compiled Python one (``(?<name>`` groups, ASCII
    classes as MATLAB's ``\\w`` and ``\\d``)."""
    flags = re.ASCII | (re.IGNORECASE if ignore_case else 0)
    return re.compile(_NAMED_GROUP.sub("(?P<", pattern), flags)


def contains_pattern(strings: Sequence[str], pattern: str, ignore_case: bool = False) -> np.ndarray:
    """``contains(s, p) | contains(s, regexpPattern(p))`` (MATLAB's channel-pattern test)."""
    rx = matlab_regex(pattern, ignore_case)
    p = pattern.lower() if ignore_case else pattern
    out = np.zeros(len(strings), dtype=bool)
    for i, s in enumerate(strings):
        t = s.lower() if ignore_case else s
        out[i] = p in t or rx.search(s) is not None
    return out


def _channel_filter(strings: Sequence[str], patterns: Sequence[str]) -> np.ndarray:
    flag = np.zeros(len(strings), dtype=bool)
    for p in patterns:
        flag |= contains_pattern(strings, p)
    return flag


def _stem(path: str) -> str:
    return fileparts(path)[1]


# ------------------------------------------------------------------------------------------
# stitch_parse_image_list_information
# ------------------------------------------------------------------------------------------


class ImageListInfo(NamedTuple):
    """The outputs of ``stitch_parse_image_list_information``, in MATLAB's order.

    ``tab``: the tiles, sorted by camera, ch, z, y, x, t, fullIter, with the columns parsed
    from the names (prefix, Iter, subIter, fullIter, camera, ch, stack, laser, abstime,
    fpgatime, x, y, z, t, did) after the file's own. ``primary_tab``: the tiles of the primary
    channel at the first time point (``None`` unless ``xcorr_mode`` is ``'primaryFirst'``).
    ``full_iter``: the unique fullIter texts. ``ch``, ``stackn``, ``nz``: float64 arrays
    (``ch`` with the primary channel first). ``cam``: one letter per camera, the primary first
    (``'N'`` without cameras in the names). ``stitch_info_fullpath``: the primary channel's
    stitch information file (or the given one).
    """

    tab: Table
    primary_tab: Table | None
    full_iter: list[str]
    ch: np.ndarray
    cam: str
    stackn: np.ndarray
    nz: np.ndarray
    specify_cam: Any
    prefix: str
    zlayer_stitch: Any
    xcorr_mode: str
    stitch_info_fullpath: str


_NEW_COLUMNS = (
    ("prefix", "text"),
    ("Iter", "double"),
    ("subIter", "text"),
    ("fullIter", "text"),
    ("camera", "char"),
    ("ch", "double"),
    ("stack", "double"),
    ("laser", "double"),
    ("abstime", "double"),
    ("fpgatime", "double"),
    ("x", "double"),
    ("y", "double"),
    ("z", "double"),
    ("t", "double"),
    ("did", "double"),
)


def _blank(kind: str, n: int) -> np.ndarray:
    if kind == "double":
        return np.zeros(n, dtype=np.float64)
    if kind == "char":
        return np.full(n, "\x00", dtype="<U1")
    out = np.empty(n, dtype=object)
    out[:] = [""] * n
    return out


def _parse_names(tab: Table, fsn: Sequence[str], expression: str, specify_cam: bool) -> Table:
    """The loop over regexpi(fsn, expression, 'names'): new columns, unmatched rows removed."""
    rx = matlab_regex(expression, ignore_case=True)
    n = tab.height
    cols: dict[str, np.ndarray] | None = None
    matched = np.ones(n, dtype=bool)
    for f, s in enumerate(fsn):
        m = rx.search(s)
        if m is None:
            matched[f] = False
            continue
        if cols is None:
            cols = {}
            for name, kind in _NEW_COLUMNS:
                if name in tab and tab.kind(name) == kind:
                    cols[name] = tab[name].copy()
                else:
                    cols[name] = _blank(kind, n)
        g = m.groupdict()
        cols["prefix"][f] = g["prefix"]
        cols["Iter"][f] = float(g["Iter"])
        cols["subIter"][f] = g["subIter"]
        cols["fullIter"][f] = g["Iter"] + g["subIter"]
        if specify_cam:
            if len(g["Cam"]) != 1:
                raise ImageListError(
                    "Unable to perform assignment because the left and right sides have a "
                    "different number of elements.",
                    "MATLAB:index:assignmentDimensionMismatch",
                )
            cols["camera"][f] = g["Cam"]
        else:
            cols["camera"][f] = "N"  # N stands for the camera when names have none
        for k in ("ch", "stack", "laser", "abstime", "fpgatime", "x", "y", "z", "t"):
            cols[k][f] = float(g[k])
        cols["did"][f] = 1.0
    out = tab.copy()
    if cols is not None:
        for name, _ in _NEW_COLUMNS:
            out[name] = cols[name]
    return out.rows(matched)


def _exists_in(names: Sequence[str], stems: Sequence[str]) -> np.ndarray:
    """``any(contains(names, stem))`` for each stem."""
    if not names:
        return np.zeros(len(stems), dtype=bool)
    joined = "\0".join(names)
    return np.array([s in joined for s in stems], dtype=bool)


def _cellfun_str2double(texts: Sequence[str], lo: int, hi: int) -> np.ndarray:
    from petakit5d.core.text import str2double

    out = []
    for x in texts:
        if len(x) < hi:
            raise ImageListError(
                f"Index exceeds the number of array elements. Index must not exceed {len(x)}.",
                "MATLAB:badsubscript",
            )
        out.append(str2double(x[lo:hi]))
    return np.asarray(out, dtype=np.float64)


def _ismember_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.array([any(np.all(r == q) for q in b) for r in a], dtype=bool)


def _check_options(xcorr_mode: str) -> None:
    if xcorr_mode not in XCORR_MODES:
        raise ImageListError(
            "Invalid value for 'xcorrMode' argument. Value must be a member of this set:\n"
            + "\n".join(f"    '{m}'" for m in XCORR_MODES),
            "MATLAB:validators:mustBeMember",
        )


def parse_image_list(
    data_path: str,
    image_list_filename: str,
    *,
    streaming: bool = False,
    stitch_info_fullpath: str = "",
    stitch_info_path: str = "",
    zarr_file: bool = False,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    use_processed_data: bool = False,
    processed_dir_str: str = "",
    timepoints: Sequence[float] | None = None,
    subtimepoints: Sequence[float] | None = None,
    xcorr_mode: str = "primaryFirst",
    primary_ch: str = "",
    online_stitch: bool = False,
) -> ImageListInfo:
    """``stitch_parse_image_list_information``: the tiles of an image list csv file.

    The rows of partial files (``_part<N>.tif``) are dropped, the tile names parsed (the
    mapped names when the file has ``mappedFilename``), rows with names the grammar does not
    match dropped, the rest sorted; tiles not in ``data_path`` (or in
    ``data_path + processed_dir_str`` when ``data_path`` has none) are dropped unless
    ``streaming``; the primary channel is ``primary_ch`` or the first channel found.
    """
    _check_options(xcorr_mode)
    tab = read_table(image_list_filename)
    if tab.kind("Filename") != "text":
        raise ImageListError("All cells must be char row vectors.", "MATLAB:REGEXP:invalidInputs")
    fn = list(tab["Filename"])
    if zarr_file:
        ext = ".zarr"
    else:
        part = matlab_regex(r"_part\d+.tif", ignore_case=True)
        partialvol = np.array([part.search(s) is not None for s in fn], dtype=bool)
        ext = ".tif"
        # remove records for partial volume files
        if partialvol.any():
            tab = tab[~partialvol]
            fn = list(tab["Filename"])
    fsn = [_stem(s) for s in fn]

    use_map_fname = "mappedFilename" in tab
    if use_map_fname:
        fsn = [_stem(s) for s in tab["mappedFilename"]]

    cam_rx = matlab_regex(r"_Cam\w_ch")
    ch_rx = matlab_regex(r"_ch[0-9]_")
    if all(cam_rx.search(s) for s in fsn):
        expression, specify_cam = TILE_EXPRESSION_CAM, True
    elif all(ch_rx.search(s) for s in fsn):
        expression, specify_cam = TILE_EXPRESSION_NOCAM, False
    else:
        raise ImageListError(
            "Unrecognized function or variable 'expression'.", "MATLAB:UndefinedFunction"
        )

    tab = _parse_names(tab, fsn, expression, specify_cam)
    if zarr_file:
        tab["Filename"] = [s.replace(".tif", ext) for s in tab["Filename"]]
    # sort by camera, ch, tile indices, time, Iter
    tab, tinds = sortrows(tab, _SORT_KEYS)
    # fsn keeps the names of all rows (MATLAB indexes the unfiltered names with tinds)
    fsn = [fsn[i] for i in tinds]

    prefixes = unique_text(tab["prefix"])
    prefix = prefixes[0] if prefixes else ""
    full_iter = unique_text(tab["fullIter"])
    ch = unique_num(tab["ch"])
    cam = "".join(sorted(set(tab["camera"].tolist())))
    stackn = unique_num(tab["stack"])
    if streaming and online_stitch:
        nz = unique_num(tab["z"])
        zlayer_stitch = True
        xcorr_mode = "all"
    else:
        nz = np.array([1.0])
        zlayer_stitch = False

    # check whether the image files in the image list file exist
    image_fnames = [e.name for e in matlab_dir(f"{data_path}/*{ext}")]
    if not image_fnames:
        if use_processed_data:
            image_fnames = [e.name for e in matlab_dir(f"{data_path}{processed_dir_str}/*{ext}")]
        else:
            raise ImageListError("The tiles do not exist!")

    # filter filenames by channel patterns
    if use_map_fname:
        full = [f"{data_path}/{x}{ext}" for x in fsn]
        include = _channel_filter(full, channel_patterns)
        if not include.all():
            tab = tab[include]
    else:
        full = [f"{data_path}/{x}" for x in image_fnames]
        include = _channel_filter(full, channel_patterns)
        image_fnames = [x for x, k in zip(image_fnames, include, strict=True) if k]

    exist = _exists_in(image_fnames, [_stem(s) for s in tab["Filename"]])
    if not exist.all():
        _warn("Some files in the image list file do not exist! Ignore them in the stitching")
        for s in tab["Filename"][~exist]:
            print(f"    {{'{s}'}}")
        print()
        if not streaming:
            tab = tab[exist]
            full_iter = unique_text(tab["fullIter"])
            ch = unique_num(tab["ch"])
            cam = "".join(sorted(set(tab["camera"].tolist())))
            stackn = unique_num(tab["stack"])

    if stitch_info_fullpath:
        if os.path.exists(stitch_info_fullpath):
            xcorr_mode = "stitchInfo"
        else:
            raise ImageListError(
                f"The user defined stitch info file {stitch_info_fullpath} does not exist!"
            )

    if xcorr_mode.lower() in ("primary", "primaryfirst"):
        ch, cam, primary_ch = _primary_channel(fsn, ch, cam, primary_ch, specify_cam)

    # separate timepoints, subtimepoints
    if timepoints is not None and len(timepoints) > 0:
        tp = np.asarray(timepoints, dtype=np.float64).reshape(-1)
        it = _cellfun_str2double(full_iter, 0, 4)
        full_iter = [x for x, i in zip(full_iter, it, strict=True) if np.any(tp == i)]
        if subtimepoints is not None and len(subtimepoints) > 0:
            stp = np.asarray(subtimepoints, dtype=np.float64).reshape(-1)
            sub = _cellfun_str2double(full_iter, 5, 9)
            if len(it) != len(sub):
                raise ImageListError(
                    "Dimensions of arrays being concatenated are not consistent.",
                    "MATLAB:catenate:dimensionMismatch",
                )
            if len(tp) != len(stp):
                raise ImageListError(
                    "Dimensions of arrays being concatenated are not consistent.",
                    "MATLAB:catenate:dimensionMismatch",
                )
            keep = _ismember_rows(np.column_stack([it, sub]), np.column_stack([tp, stp]))
            full_iter = [x for x, k in zip(full_iter, keep, strict=True) if k]

    # predefine stitchInfo when xcorrMode is 'primaryFirst'
    primary_tab = None
    if xcorr_mode == "primaryFirst":
        primary_tab = _primary_rows(tab, ch, cam, full_iter, stackn, nz, zlayer_stitch)
        p_laser = unique_num(primary_tab["laser"])
        p_abstime = unique_num(primary_tab["abstime"])
        p_fpgatime = primary_tab["fpgatime"][0]
        p_z = unique_num(primary_tab["z"])[0]
        if len(p_laser) > 1:
            p_laser = p_laser[:1]
        pz_str = _sprintf("_%03dz", p_z) if zlayer_stitch else ""
        if specify_cam:
            stitch_info_fullpath = _sprintf(
                "%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat",
                stitch_info_path,
                prefix,
                full_iter[0],
                cam[0],
                ch[0],
                stackn[0],
                p_laser,
                p_abstime,
                p_fpgatime,
                pz_str,
            )
        else:
            stitch_info_fullpath = _sprintf(
                "%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat",
                stitch_info_path,
                prefix,
                full_iter[0],
                ch[0],
                stackn[0],
                p_laser,
                p_abstime,
                p_fpgatime,
                pz_str,
            )

    return ImageListInfo(
        tab,
        primary_tab,
        full_iter,
        ch,
        cam,
        stackn,
        nz,
        specify_cam,
        prefix,
        zlayer_stitch,
        xcorr_mode,
        stitch_info_fullpath,
    )


def _index1(a, what: str):
    if len(a) == 0:
        raise ImageListError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    return a[0]


def _primary_rows(tab, ch, cam, full_iter, stackn, nz, zlayer_stitch) -> Table:
    mask = (
        (tab["ch"] == _index1(ch, "Ch"))
        & (tab["camera"] == _index1(cam, "Cam"))
        & (tab["fullIter"] == _index1(full_iter, "fullIter"))
        & (tab["stack"] == _index1(stackn, "stackn"))
    )
    if zlayer_stitch:
        mask &= tab["z"] == _index1(nz, "nz")
    primary_tab = tab[np.asarray(mask, dtype=bool)]
    if primary_tab.height == 0:
        raise ImageListError(
            "The Image List Info for the primary channel for the first time point does not exist!"
        )
    return primary_tab


def _primary_channel(fsn, ch, cam, primary_ch, specify_cam):
    """The primary channel (given or the first found) moved to the front of ch and cam."""
    from petakit5d.core.text import str2double

    def found(pc: str) -> bool:
        key = f"_{pc}_".lower()
        return any(key in s.lower() for s in fsn)

    p_ch = None
    if primary_ch:
        print(f"The primary channel is {primary_ch}.")
        if specify_cam and not matlab_regex("Cam[a-z]_ch[0-9]", True).search(primary_ch):
            raise ImageListError("primaryCh must be empty or with format 'Cam[a-z]_ch[0-9]'")
        if not specify_cam and not matlab_regex("ch[0-9]", True).search(primary_ch):
            raise ImageListError("primaryCh must be empty or with format 'ch[0-9]'")
        p_cam = primary_ch[3] if specify_cam else None
        p_ch = str2double(primary_ch[-1])
        if not found(primary_ch):
            raise ImageListError(f"The given primary channel {primary_ch} does not exist!")
    else:
        _warn(
            "The primary channel is not set, use the first available channel as primary channel..."
        )
        p_cam = None
        for c_cam in cam:
            for c in ch:
                if specify_cam:
                    primary_ch = _sprintf("Cam%s_ch%d", c_cam, c)
                else:
                    primary_ch = _sprintf("ch%d", c)
                if found(primary_ch):
                    p_cam, p_ch = c_cam, c
                    break
            if p_cam is not None:
                break
        print(f"Set primary channel as {primary_ch}.")

    # reorder Ch and Cam to make primary channel the first
    first = _index1(ch, "Ch")
    if p_ch is None:
        raise ImageListError("Unrecognized function or variable 'pCh'.", "MATLAB:UndefinedFunction")
    if first != p_ch:
        ch = np.concatenate(([p_ch], ch[ch != p_ch]))
    if specify_cam and _index1(cam, "Cam").lower() != p_cam.lower():
        rx = matlab_regex(p_cam, ignore_case=True)
        inds = [m.start() + 1 for m in rx.finditer(cam)]
        if inds:
            i = inds[0]
            order = inds + list(range(1, i)) + list(range(i + 1, len(cam) + 1))
            cam = "".join(cam[k - 1] for k in order)
        else:
            cam = ""
    return ch, cam, primary_ch


# ------------------------------------------------------------------------------------------
# stitch_parse_multi_loc_image_list_information
# ------------------------------------------------------------------------------------------


def _unique_logical(values: Sequence[bool]):
    u = sorted({bool(v) for v in values})
    return u[0] if len(u) == 1 else np.array(u, dtype=bool)


def parse_multi_loc_image_list(
    data_paths: Sequence[str],
    image_list_filenames: Sequence[str],
    *,
    streaming: bool = False,
    stitch_info_fullpath: str = "",
    stitch_info_path: str = "",
    zarr_file: bool = False,
    channel_patterns: Sequence[str] = ("CamA", "CamB"),
    use_processed_data: bool = False,
    processed_dir_str: str = "",
    timepoints: Sequence[float] | None = None,
    subtimepoints: Sequence[float] | None = None,
    xcorr_mode: str = "primaryFirst",
    primary_ch: str = "",
    online_stitch: bool = False,
) -> ImageListInfo:
    """``stitch_parse_multi_loc_image_list_information``: one image list per folder, joined.

    ``did`` is the folder number (1-based). ``ch`` and ``cam`` keep the order of first
    appearance over the folders; ``specify_cam`` and ``zlayer_stitch`` are a bool when all
    folders agree (else MATLAB's ``[false; true]`` as an array); ``prefix`` is the sorted
    unique characters of the folders' prefixes (all must have the same length).
    """
    _check_options(xcorr_mode)
    infos = []
    for d, (dp, fn) in enumerate(zip(data_paths, image_list_filenames, strict=False)):
        info = parse_image_list(
            dp,
            fn,
            streaming=streaming,
            online_stitch=online_stitch,
            stitch_info_fullpath=stitch_info_fullpath,
            stitch_info_path=stitch_info_path,
            zarr_file=zarr_file,
            channel_patterns=channel_patterns,
            use_processed_data=use_processed_data,
            processed_dir_str=processed_dir_str,
            timepoints=timepoints,
            subtimepoints=subtimepoints,
            xcorr_mode=xcorr_mode,
            primary_ch=primary_ch,
        )
        xcorr_mode = info.xcorr_mode
        info.tab["did"] = np.full(info.tab.height, float(d + 1))
        infos.append(info)
    if len(image_list_filenames) < len(data_paths):
        raise ImageListError(
            "Index exceeds the number of array elements. "
            f"Index must not exceed {len(image_list_filenames)}.",
            "MATLAB:badsubscript",
        )

    tab = Table.vcat([i.tab for i in infos])
    full_iter = unique_text(x for i in infos for x in i.full_iter)
    ch = np.asarray(_stable_num(np.concatenate([i.ch for i in infos])), dtype=np.float64)
    cam = "".join(dict.fromkeys("".join(i.cam for i in infos)))
    stackn = unique_num(np.concatenate([i.stackn for i in infos]))
    nz = unique_num(np.concatenate([i.nz for i in infos]))
    specify_cam = _unique_logical([i.specify_cam for i in infos])
    prefixes = [i.prefix for i in infos]
    if len({len(p) for p in prefixes}) > 1:
        raise ImageListError(
            "Dimensions of arrays being concatenated are not consistent.",
            "MATLAB:catenate:dimensionMismatch",
        )
    prefix = "".join(sorted(set("".join(prefixes)), key=ord))
    zlayer_stitch = _unique_logical([i.zlayer_stitch for i in infos])
    stitch_info_fullpath = infos[0].stitch_info_fullpath
    primary_tab = infos[-1].primary_tab
    if xcorr_mode == "primaryFirst":
        if not isinstance(zlayer_stitch, bool):
            zl = bool(np.all(zlayer_stitch))
        else:
            zl = zlayer_stitch
        primary_tab = _primary_rows(tab, ch, cam, full_iter, stackn, nz, zl)
    return ImageListInfo(
        tab,
        primary_tab,
        full_iter,
        ch,
        cam,
        stackn,
        nz,
        specify_cam,
        prefix,
        zlayer_stitch,
        xcorr_mode,
        stitch_info_fullpath,
    )


def _stable_num(values: np.ndarray) -> list[float]:
    out: list[float] = []
    for v in values.tolist():
        if v != v or v not in out:
            out.append(v)
    return out


# ------------------------------------------------------------------------------------------
# parseTileInfo
# ------------------------------------------------------------------------------------------


def _parse_tcxyz(filenames: Sequence[str], tile_pattern: str, specify_cam: bool = False):
    m = re.search(r"\D+", tile_pattern, re.ASCII)
    if m is None:
        raise ImageListError(
            f"The tile pattern {tile_pattern} does not contain a non-digit character to "
            "specify the pattern!"
        )
    pat = f"{m.group(0)}(\\d+)" if m.start() == 0 else f"(\\d+){m.group(0)}"
    rx = matlab_regex(pat, ignore_case=True)
    tile_inds = np.empty(len(filenames), dtype=np.float64)
    for i, s in enumerate(filenames):
        t = rx.search(s)
        tile_inds[i] = float(t.group(1)) if t is not None else NAN
    cam_inds = None
    cam_patterns: list[str] = []
    if specify_cam:
        crx = matlab_regex("Cam([A-Z])", ignore_case=True)
        letters = []
        for s in filenames:
            t = crx.search(s)
            letters.append(t.group(1).upper() if t is not None else None)
        if any(x is None for x in letters):
            # ismember of a cell holding [] (a name without a camera) is an error in MATLAB
            raise ImageListError(
                "Input A of class cell and input B of class cell must be cell arrays of "
                "character vectors, unless one is a character vector.",
                "MATLAB:ISMEMBER:InputClass",
            )
        cam_patterns = unique_text(x for x in letters if x is not None)
        cam_inds = np.array(
            [cam_patterns.index(x) + 1 if x is not None else NAN for x in letters],
            dtype=np.float64,
        )
    return tile_inds, cam_inds, cam_patterns


def parse_tile_info(
    filenames: Sequence[str], tile_patterns: Sequence[str]
) -> tuple[np.ndarray, list[str]]:
    """``parseTileInfo``: tile indices of the t, c, x, y, z patterns (``'0000t'``, ``'ch0'``,
    ``'000x'``, ...) for each name, NaN where a pattern is missing; column 6 is the camera's
    1-based index in ``cam_patterns`` when the channel pattern names one (``CamA_ch0``)."""
    tile_info = np.zeros((len(filenames), 6), dtype=np.float64)
    cam_patterns: list[str] = []
    for i in range(5):
        pattern = tile_patterns[i]
        if not pattern:
            continue
        if i == 1:
            specify_cam = matlab_regex("Cam[A-Z]", ignore_case=True).search(pattern) is not None
            inds, cam_inds, cam_patterns = _parse_tcxyz(filenames, pattern, specify_cam)
            if specify_cam:
                tile_info[:, 5] = cam_inds
        else:
            inds, _, _ = _parse_tcxyz(filenames, pattern)
        tile_info[:, i] = inds
    return tile_info, cam_patterns


# ------------------------------------------------------------------------------------------
# stitch_process_filenames
# ------------------------------------------------------------------------------------------


def process_filenames(
    tile_fullpaths: Sequence[str],
    processed_dirstr: str,
    stitch_mip,
    resample,
    zarr_file: bool,
    process_tiles: bool,
) -> tuple[list[str], list[str], list[str], str]:
    """``stitch_process_filenames``: input paths, zarr paths, names and the zarr folder name of
    the tiles (``zarr``, ``zarr_processed`` for zarr input, ``zarr_<ry>_<rx>_<rz>`` when
    resampled, ``''`` without tile processing)."""
    ext = ".zarr" if zarr_file else ".tif"
    suffix_str = "z"
    input_fullpaths: list[str] = []
    zarr_fullpaths: list[str] = []
    fsnames: list[str] = []
    zarr_pathstr = None
    rs_in = np.asarray([] if resample is None else resample, dtype=np.float64).reshape(-1)
    for tile in tile_fullpaths:
        data_path, fsname, _ = fileparts(tile)
        data_path = data_path.rstrip("/")
        tile_path = f"{data_path}/{processed_dirstr}" if processed_dirstr else data_path
        fsnames.append(fsname)
        if process_tiles:
            if rs_in.size and np.any(rs_in != 1):
                # complete resample to 3d
                rs = np.concatenate((np.full(4 - rs_in.size, rs_in[0]), rs_in[1:]))
                zarr_pathstr = _sprintf("zarr_%d_%d_%d", rs[0], rs[1], rs[2])
            else:
                zarr_pathstr = "zarr_processed" if zarr_file else "zarr"
        else:
            zarr_pathstr = ""
        if np.any(np.asarray(stitch_mip)):
            input_fullpaths.append(f"{tile_path}/MIPs/{fsname}_MIP_{suffix_str}{ext}")
            zarr_fullpaths.append(f"{tile_path}/MIPs/{zarr_pathstr}/{fsname}_MIP_{suffix_str}.zarr")
        else:
            input_fullpaths.append(f"{tile_path}/{fsname}{ext}")
            zarr_fullpaths.append(f"{tile_path}/{zarr_pathstr}/{fsname}.zarr")
    if zarr_pathstr is None:
        raise ImageListError(
            'Output argument "zarrPathstr" (and possibly others) not assigned a value in the '
            'execution with "stitch_process_filenames" function.',
            "MATLAB:unassignedOutputs",
        )
    return input_fullpaths, zarr_fullpaths, fsnames, zarr_pathstr
