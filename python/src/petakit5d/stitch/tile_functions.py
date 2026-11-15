"""Tile functions of the stitching and the tiff/zarr conversion: the edge erosion
(``erodeVolumeBy2DProjection``) and the reader of tile functions written as MATLAB code.

MATLAB keeps a tile's processing function as MATLAB code in a function file (``usrFun`` of a
``.mat`` file, or the first line of a ``.txt`` file): the stitching's frame writes
``@(x)erodeVolumeBy2DProjection(x,<n>)`` (skewed space with a tile offset:
``@(x)erodeVolumeBy2DProjection(x+<offset>,<n>)``) for its edge erosion, and users write the
flat-field correction and its composition with the erosion, as the demos do::

    @(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'ff.tif','bg.tif',1,0.4,false),1)

MATLAB runs such code with ``str2func``. The port never evaluates code: ``parse_tile_function``
reads these standard forms with a small parser and builds the function from the port's own
functions. The forms read (whitespace anywhere, any variable name)::

    @(x) EXPR
    EXPR  := x
           | erodeVolumeBy2DProjection(EXPR, esize)
           | processFFCorrectionFrame(EXPR, LSImage, BackgroundImage, constOffset
                                      [, LowerLimit [, removeFFImBackground [, LSRescale
                                      [, castDataType]]]])
           | EXPR + n | EXPR - n | n + EXPR | n - EXPR | (EXPR)

with constants ``n`` and arguments of MATLAB's literal forms: numbers (``1``, ``-2.5``,
``.4``, ``4.000000e-01``, ``0000000001``, sums and differences of them), texts in ``'...'``
or ``"..."`` (doubled quotes inside), ``true``, ``false`` and ``[]``. ``+`` and ``-`` follow
MATLAB's classes: an integer volume rounds (half away from zero) and saturates, a single
volume stays single. Any other code raises ``NotImplementedError`` (identifier
``petakit5d:convert:matlabCodeFunction``); a function of the user's own is given as a Python
function instead (``io.convert.user_function``: ``'package.module:function'`` or a record of
its arguments, also in the function files).

``erode_volume_by_2d_projection`` is erodeVolumeBy2DProjection.m: the edges are eroded by the
erosion of the volume's max projection along y (the xz plane) instead of the whole volume, and
the first and last ``esize`` rows along y are set to 0.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "TileFunctionError",
    "erode_volume_by_2d_projection",
    "parse_tile_function",
]

UNSUPPORTED_ID = "petakit5d:convert:matlabCodeFunction"


class TileFunctionError(ValueError):
    """An error MATLAB raises when it runs the function; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _to_matlab_result(result: Any) -> Any:
    from petakit5d.core.arrays import to_matlab

    return to_matlab(result)


# ---------------------------------------------------------------------------------------------
# erodeVolumeBy2DProjection
# ---------------------------------------------------------------------------------------------


@register(
    "erodeVolumeBy2DProjection",
    Param("vol", "@isnumeric", "Volume (z, y, x), or an image (y, x).", convert="volume"),
    Param("esize", "@isnumeric", "Erosion size (voxels).", convert="same"),
    Param("out", None, "Array to write the result into (may be vol itself).", python_only=True),
    Param(
        "numThreads",
        None,
        "Threads (None: the CPUs available); never changes the result.",
        py="num_threads",
        python_only=True,
    ),
    task="S02",
    matlab_out=_to_matlab_result,
)
def erode_volume_by_2d_projection(
    vol: np.ndarray,
    esize: Any,
    *,
    out: np.ndarray | None = None,
    num_threads: int | None = None,
) -> np.ndarray:
    """Erode the edges of a volume (z, y, x) by ``esize`` voxels through its max projection
    along y (erodeVolumeBy2DProjection): the voxels whose (z, x) column is eroded in that
    projection (a square of ``2 esize + 1``, outside the projection as empty) and the first and
    last ``esize`` rows along y become 0. An image (y, x) keeps its values (MATLAB's
    ``vol .* (vol ~= 0)``). An ``esize`` of 0 returns ``vol`` itself. ``out``: a C-contiguous
    array of the input's shape and class to write into, ``vol`` itself included (the
    projection is taken first)."""
    from petakit5d.matlabcompat.morphology import imerode, strel

    vol = np.asarray(vol)
    e_val = float(np.asarray(esize, dtype=np.float64).reshape(-1)[0])
    if e_val == 0:
        if out is None or out is vol:
            return vol
        out[...] = vol
        return out
    if e_val != int(e_val) or e_val < 0:
        # (MATLAB indexes with it)
        raise TileFunctionError(
            "Array indices must be positive integers or logical values.", "MATLAB:badsubscript"
        )
    e = int(e_val)
    if out is None:
        out = np.empty_like(vol)
    elif out.shape != vol.shape or out.dtype != vol.dtype or not out.flags.c_contiguous:
        raise ValueError("out must be C-contiguous, of the volume's shape and class")
    if vol.ndim == 2:
        return np.multiply(vol, (vol != 0).astype(vol.dtype), out=out)
    nz, ny, nx = vol.shape
    if ny == 0 or nz == 0:
        out[...] = vol
        return out
    chunks = _chunks(nz, vol.nbytes, num_threads)
    # MATLAB's squeeze(max(vol, [], 1) > 0): the projection along y, here (z, x) (MATLAB's max
    # leaves NaN out, as fmax does)
    mip = np.empty((nz, nx), bool)
    _each(chunks, lambda a, b: np.greater(np.fmax.reduce(vol[a:b], axis=1), 0, out=mip[a:b]))
    if mip.all():
        erode = np.zeros_like(mip)
        erode[e : nz - e, e : nx - e] = mip[e : nz - e, e : nx - e]
    else:
        pad = np.zeros((nz + 2, nx + 2), bool)
        pad[1:-1, 1:-1] = mip
        pad = imerode(pad, strel("square", e * 2 + 1))
        erode = pad[1:-1, 1:-1]
    mask = erode.astype(vol.dtype)[:, np.newaxis, :]
    lo, hi = min(e, ny), max(ny - e, 0)

    def part(a: int, b: int) -> None:
        # vol .* mask, then the first and last esize rows along y set to 0
        np.multiply(vol[a:b], mask[a:b], out=out[a:b])
        out[a:b, :lo, :] = 0
        out[a:b, hi:, :] = 0

    _each(chunks, part)
    return out


def _chunks(nz: int, nbytes: int, num_threads: int | None) -> list[tuple[int, int]]:
    """z ranges for the threads (one range for small volumes)"""
    from petakit5d.core.threads import available_cpus

    n = available_cpus() if num_threads is None else max(1, int(num_threads))
    n = max(1, min(n, nz, nbytes // (1 << 22)))
    edges = np.linspace(0, nz, n + 1).astype(int)
    return [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:], strict=True) if b > a]


def _each(chunks: list[tuple[int, int]], fn) -> None:
    """fn(z0, z1) for each range, in threads (NumPy releases the GIL in its loops)"""
    if len(chunks) == 1:
        fn(*chunks[0])
        return
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(len(chunks)) as ex:
        for f in [ex.submit(fn, a, b) for a, b in chunks]:
            f.result()


# ---------------------------------------------------------------------------------------------
# MATLAB's volume + double and volume - double
# ---------------------------------------------------------------------------------------------


def _arith(vol: np.ndarray, op: str, d: float, num_first: bool) -> np.ndarray:
    """MATLAB's ``vol + d``, ``vol - d``, ``d + vol`` or ``d - vol`` for a double scalar ``d``
    (``core.dtypes``' rules), as a new array, by z ranges in threads"""
    from petakit5d.core import dtypes

    vol = np.asarray(vol)
    if vol.dtype == np.float32 and float(np.float32(d)) == d:
        # single with a double that single holds: rounding the double sum to single is the
        # single sum (both operands single; double has more than 2 * 24 + 1 bits)
        ufunc = np.add if op == "+" else np.subtract
        s = np.float32(d)
        if vol.ndim < 3:
            return ufunc(s, vol) if num_first else ufunc(vol, s)
        res = np.empty_like(vol)

        def part32(a: int, b: int) -> None:
            if num_first:
                ufunc(s, vol[a:b], out=res[a:b])
            else:
                ufunc(vol[a:b], s, out=res[a:b])

        _each(_chunks(vol.shape[0], vol.nbytes, None), part32)
        return res
    f = dtypes.plus if op == "+" else dtypes.minus
    if vol.ndim < 3 or vol.dtype.kind not in "iuf" or vol.dtype.itemsize == 8:
        return np.asarray(f(d, vol) if num_first else f(vol, d))
    res = np.empty_like(vol)
    if vol.dtype.kind in "iu" and d == int(d) and abs(d) < 2**40:
        # integers up to 32 bits and a whole d: the exact sum in int64, saturated
        info = np.iinfo(vol.dtype)
        k = int(d)

        def part(a: int, b: int) -> None:
            w = vol[a:b].astype(np.int64)
            if op == "+":
                w += k
            elif num_first:
                np.subtract(k, w, out=w)
            else:
                w -= k
            np.clip(w, info.min, info.max, out=w)
            res[a:b] = w

    else:

        def part(a: int, b: int) -> None:
            res[a:b] = f(d, vol[a:b]) if num_first else f(vol[a:b], d)

    # (dtypes works in double: z ranges of at most about 32 MB of doubles)
    n = max(1, min(vol.shape[0], vol.size * 8 // (1 << 25)))
    edges = np.linspace(0, vol.shape[0], n + 1).astype(int)
    ranges = [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:], strict=True) if b > a]
    from concurrent.futures import ThreadPoolExecutor

    from petakit5d.core.threads import available_cpus

    with ThreadPoolExecutor(max(1, min(available_cpus(), len(ranges)))) as ex:
        for fut in [ex.submit(part, a, b) for a, b in ranges]:
            fut.result()
    return res


# ---------------------------------------------------------------------------------------------
# the parser
# ---------------------------------------------------------------------------------------------

_TOKEN = re.compile(
    r"\s*(?:"
    r"(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)"
    r"|(?P<name>[A-Za-z]\w*)"
    r"|(?P<char>'(?:[^']|'')*')"
    r"|(?P<str>\"(?:[^\"]|\"\")*\")"
    r"|(?P<op>[@(),+\-\[\]])"
    r")"
)

# the functions read, with the number of their arguments (first: the volume)
_FUNCTIONS = {
    "erodeVolumeBy2DProjection": (2, 2),
    "processFFCorrectionFrame": (4, 8),
}


class _Unsupported(Exception):
    pass


def _tokens(text: str) -> list[tuple[str, Any]]:
    out: list[tuple[str, Any]] = []
    pos, n = 0, len(text)
    while pos < n:
        if text[pos:].strip() == "":
            break
        m = _TOKEN.match(text, pos)
        if m is None or m.end() == pos:
            raise _Unsupported(f"unexpected {text[pos:].strip()[:20]!r}")
        kind = m.lastgroup
        val = m.group(kind)
        if kind in ("char", "str") and out and out[-1][0] in ("num", "name", ")", "]"):
            # (a quote after a value is MATLAB's transpose)
            raise _Unsupported("transpose")
        if kind == "num":
            out.append(("num", float(val)))
        elif kind == "name":
            out.append(("name", val))
        elif kind == "char":
            out.append(("text", val[1:-1].replace("''", "'")))
        elif kind == "str":
            out.append(("text", val[1:-1].replace('""', '"')))
        else:
            out.append((val, val))
        pos = m.end()
    return out


class _Parser:
    """recursive descent over the tokens; nodes are tuples:
    ('var',), ('const', value), ('call', name, [nodes]), ('arith', op, vol node, d, num_first)"""

    def __init__(self, tokens: list[tuple[str, Any]]):
        self.t = tokens
        self.i = 0
        self.var = ""

    def peek(self) -> str | None:
        return self.t[self.i][0] if self.i < len(self.t) else None

    def take(self, kind: str) -> Any:
        if self.peek() != kind:
            got = self.t[self.i][1] if self.i < len(self.t) else "the end"
            raise _Unsupported(f"expected {kind!r}, got {got!r}")
        v = self.t[self.i][1]
        self.i += 1
        return v

    def function(self) -> tuple[str, Any]:
        self.take("@")
        self.take("(")
        self.var = self.take("name")
        self.take(")")
        node = self.expr()
        if self.peek() is not None:
            raise _Unsupported(f"unexpected {self.t[self.i][1]!r}")
        if node[0] == "const":
            raise _Unsupported("a function of no volume")
        return self.var, node

    def expr(self) -> tuple:
        node = self.unary()
        while self.peek() in ("+", "-"):
            op = self.take(self.peek())
            node = self.combine(op, node, self.unary())
        return node

    def combine(self, op: str, a: tuple, b: tuple) -> tuple:
        if a[0] == "const" and b[0] == "const":
            x, y = self.number(a), self.number(b)
            return ("const", x + y if op == "+" else x - y)
        if a[0] != "const" and b[0] == "const":
            return ("arith", op, a, self.number(b), False)
        if a[0] == "const" and b[0] != "const":
            return ("arith", op, b, self.number(a), True)
        raise _Unsupported("arithmetic of two volumes")

    @staticmethod
    def number(node: tuple) -> float:
        v = node[1]
        if isinstance(v, bool) or not isinstance(v, float):
            raise _Unsupported("arithmetic with a value that is not a number")
        return v

    def unary(self) -> tuple:
        if self.peek() in ("+", "-"):
            op = self.take(self.peek())
            node = self.unary()
            if node[0] != "const":
                raise _Unsupported("a sign on a volume")
            v = self.number(node)
            return ("const", -v if op == "-" else v)
        return self.primary()

    def primary(self) -> tuple:
        kind = self.peek()
        if kind == "num":
            return ("const", self.take("num"))
        if kind == "text":
            return ("const", self.take("text"))
        if kind == "(":
            self.take("(")
            node = self.expr()
            self.take(")")
            return node
        if kind == "[":
            self.take("[")
            self.take("]")
            return ("const", None)
        if kind == "name":
            name = self.take("name")
            if name == self.var:
                if self.peek() == "(":
                    raise _Unsupported("indexing the volume")
                return ("var",)
            if name in ("true", "false") and self.peek() != "(":
                return ("const", name == "true")
            if name not in _FUNCTIONS:
                raise _Unsupported(f"the function {name}")
            self.take("(")
            args = [self.expr()]
            while self.peek() == ",":
                self.take(",")
                args.append(self.expr())
            self.take(")")
            if args[0][0] == "const" or any(a[0] != "const" for a in args[1:]):
                raise _Unsupported(f"{name} of other than a volume and constants")
            return ("call", name, args)
        raise _Unsupported(f"unexpected {self.t[self.i][1]!r}" if kind else "an empty body")


def _evaluate(node: tuple, x: np.ndarray, owned: bool) -> tuple[np.ndarray, bool]:
    """the node's volume, and whether it is a new array (which later steps may overwrite)"""
    kind = node[0]
    if kind == "var":
        return x, owned
    if kind == "arith":
        _, op, sub, d, num_first = node
        v, _ = _evaluate(sub, x, owned)
        return _arith(v, op, d, num_first), True
    _, name, args = node
    v, own = _evaluate(args[0], x, owned)
    consts = [a[1] for a in args[1:]]
    if name == "erodeVolumeBy2DProjection":
        esize = consts[0]
        if isinstance(esize, (bool, str)) or esize is None:
            raise TileFunctionError(
                "The value of 'esize' is invalid. It must satisfy the function: isnumeric.",
                "MATLAB:InputParser:ArgumentFailedValidation",
            )
        reuse = own and v.flags.c_contiguous and v.flags.writeable
        r = erode_volume_by_2d_projection(v, esize, out=v if reuse else None)
        return r, own or r is not v
    from petakit5d.flatfield import flatfield_correct_file

    return np.asarray(flatfield_correct_file(v, *consts)), True


def _unsupported(text: str, why: str) -> NotImplementedError:
    err = NotImplementedError(
        f"MATLAB code as a tile function ({text!r}) is not one of the forms the port reads "
        f"({why}): @(x) with x, erodeVolumeBy2DProjection(...), processFFCorrectionFrame(...), "
        "+ and - of constants; give a Python function ('package.module:function') instead"
    )
    err.identifier = UNSUPPORTED_ID  # type: ignore[attr-defined]
    return err


def parse_tile_function(text: str) -> Callable[[np.ndarray], np.ndarray]:
    """The function of a tile function written as MATLAB code (the forms in the module's
    docstring), on (z, y, x) volumes or (y, x) images. Raises ``NotImplementedError`` (with
    ``identifier``) for other code; errors MATLAB raises only when it runs the function are
    raised when it runs here too."""
    src = str(text)
    try:
        _, node = _Parser(_tokens(src.strip())).function()
    except _Unsupported as why:
        raise _unsupported(src, str(why)) from None
    for name, args in _calls(node):
        lo, hi = _FUNCTIONS[name]
        if not lo <= len(args) <= hi:
            raise TileFunctionError(
                f"{'Too many' if len(args) > hi else 'Not enough'} input arguments.",
                "MATLAB:TooManyInputs" if len(args) > hi else "MATLAB:minrhs",
            )

    def tile_function(vol: np.ndarray) -> np.ndarray:
        return _evaluate(node, np.asarray(vol), False)[0]

    tile_function.matlab_code = src  # type: ignore[attr-defined]
    return tile_function


def _calls(node: tuple) -> Sequence[tuple[str, list]]:
    if node[0] == "call":
        return [(node[1], node[2]), *_calls(node[2][0])]
    if node[0] == "arith":
        return _calls(node[2])
    return []
