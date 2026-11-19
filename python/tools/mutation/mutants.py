"""Single-line mutants of Python modules: where they can go (sites) and how one is made.

A site is a span of one physical line of a module (1-based line, byte columns of the line's
UTF-8 text, end exclusive), the text it holds and its replacement. A mutant is a site of a
module, with the SHA-256 of the module's text, so that a mutant is never applied to another
version of the file. Operators:

  cmp     a comparison to its neighbour: < <=, > >=, == !=, is / is not, in / not in
  arith   a binary or augmented operator: + -, * /, / *, // /, % //, ** *, & |, | &, ^ &,
          << >>, @ *
  num     an int literal n -> n + 1; a float x -> 2 x (0.0 -> 1.0)
  bool    True <-> False, and <-> or, `not x` -> `x`
  stmt    a statement alone on its line dropped (pass): a call, an augmented assignment, an
          assignment to an attribute or an item, return, raise, continue, break
  cond    the test of an if, elif or while on one line negated: c -> not (c)
  swap    two neighbouring items of a tuple or list display of 2 to 4 numbers or names
  call    min <-> max (and the NumPy forms), floor <-> ceil, any <-> all; a one-argument
          rounding or cast dropped: mround(x), round(x), int(x), float(x), abs(x), floor(x),
          ceil(x), trunc(x), fix(x) -> (x)

Never mutated (their mutants are diagnostic only or not code): strings, f-strings, docstrings,
annotations, assert statements, ``if TYPE_CHECKING:`` blocks, ``__all__``, ``__repr__`` and
``__str__``, the types of isinstance and of except clauses, the arguments of print and of
logger calls, and flush calls. A site is used only if the mutated module compiles and its
syntax tree differs from the original's (``valid``).

This module imports only the standard library.
"""

from __future__ import annotations

import ast
import hashlib
import random
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass

OPERATORS = ("cmp", "arith", "num", "bool", "stmt", "cond", "swap", "call")

_CMP = {
    ast.Lt: ("<", "<="),
    ast.LtE: ("<=", "<"),
    ast.Gt: (">", ">="),
    ast.GtE: (">=", ">"),
    ast.Eq: ("==", "!="),
    ast.NotEq: ("!=", "=="),
    ast.Is: ("is", "is not"),
    ast.IsNot: ("is not", "is"),
    ast.In: ("in", "not in"),
    ast.NotIn: ("not in", "in"),
}
_ARITH = {
    ast.Add: ("+", "-"),
    ast.Sub: ("-", "+"),
    ast.Mult: ("*", "/"),
    ast.Div: ("/", "*"),
    ast.FloorDiv: ("//", "/"),
    ast.Mod: ("%", "//"),
    ast.Pow: ("**", "*"),
    ast.BitAnd: ("&", "|"),
    ast.BitOr: ("|", "&"),
    ast.BitXor: ("^", "&"),
    ast.LShift: ("<<", ">>"),
    ast.RShift: (">>", "<<"),
    ast.MatMult: ("@", "*"),
}
_SWAP_NAMES = {
    "min": "max",
    "max": "min",
    "minimum": "maximum",
    "maximum": "minimum",
    "fmin": "fmax",
    "fmax": "fmin",
    "nanmin": "nanmax",
    "nanmax": "nanmin",
    "amin": "amax",
    "amax": "amin",
    "argmin": "argmax",
    "argmax": "argmin",
    "floor": "ceil",
    "ceil": "floor",
    "any": "all",
    "all": "any",
}
_DROP_CALLS = {"mround", "round", "int", "float", "abs", "floor", "ceil", "trunc", "fix"}
# calls whose arguments are diagnostic only (printing, logging)
_DIAGNOSTIC_CALLS = {"print", "debug", "info", "warning", "error", "exception", "critical"}
_LOGGER_NAMES = {"log", "logger", "LOG", "LOGGER", "_log", "_logger", "logging"}


@dataclass(frozen=True)
class Site:
    """a mutation of one line: bytes [col, end) of the line hold ``old``, which becomes
    ``new``"""

    line: int
    col: int
    end: int
    op: str
    old: str
    new: str


@dataclass
class Mutant:
    """a site of a module (``file`` relative to the tree's python/ folder)"""

    id: str
    area: str
    file: str
    sha256: str
    line: int
    col: int
    end: int
    op: str
    old: str
    new: str
    src_line: str
    mut_line: str

    def record(self) -> dict:
        return asdict(self)

    @classmethod
    def from_record(cls, rec: dict) -> Mutant:
        return cls(**{k: rec[k] for k in cls.__dataclass_fields__})


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


# --- finding sites ---------------------------------------------------------------------------


class _Finder(ast.NodeVisitor):
    def __init__(self, lines: Sequence[bytes]) -> None:
        self.lines = lines
        self.sites: list[Site] = []

    # helpers

    def _text(self, line: int, col: int, end: int) -> str:
        return self.lines[line - 1][col:end].decode("utf-8")

    def _one_line(self, node: ast.AST) -> bool:
        return getattr(node, "lineno", None) is not None and node.lineno == node.end_lineno

    def _add(self, line: int, col: int, end: int, op: str, new: str) -> None:
        self.sites.append(Site(line, col, end, op, self._text(line, col, end), new))

    def _between(self, a: ast.AST, b: ast.AST, token: str) -> tuple[int, int, int] | None:
        """the span of ``token`` between the end of a and the start of b (one line; only
        brackets and blanks besides it)"""
        if a.end_lineno != b.lineno or a.lineno != a.end_lineno or b.lineno != b.end_lineno:
            return None
        seg = self.lines[a.end_lineno - 1][a.end_col_offset : b.col_offset].decode("utf-8")
        core = seg.strip(" \t()[]")
        if " ".join(core.split()) != token:
            return None
        start = a.end_col_offset + len(seg[: seg.index(core)].encode("utf-8"))
        return a.end_lineno, start, start + len(core.encode("utf-8"))

    # skipped subtrees

    def visit_JoinedStr(self, node: ast.JoinedStr) -> None:
        return

    def visit_Assert(self, node: ast.Assert) -> None:
        return

    def visit_FunctionDef(self, node: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
        if node.name in ("__repr__", "__str__"):
            return
        for d in node.decorator_list:
            self.visit(d)
        self._visit_arguments(node.args)
        for s in node.body:
            self.visit(s)

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_Lambda(self, node: ast.Lambda) -> None:
        self._visit_arguments(node.args)
        self.visit(node.body)

    def _visit_arguments(self, args: ast.arguments) -> None:
        # defaults only (annotations are not code)
        for d in [*args.defaults, *[d for d in args.kw_defaults if d is not None]]:
            self.visit(d)

    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        if node.value is not None:
            self.visit(node.value)

    def visit_arg(self, node: ast.arg) -> None:
        return

    def visit_If(self, node: ast.If) -> None:
        t = node.test
        if (isinstance(t, ast.Name) and t.id == "TYPE_CHECKING") or (
            isinstance(t, ast.Attribute) and t.attr == "TYPE_CHECKING"
        ):
            for s in node.orelse:
                self.visit(s)
            return
        self._negate(node.test, node)
        self.generic_visit(node)

    def visit_While(self, node: ast.While) -> None:
        self._negate(node.test, node)
        self.generic_visit(node)

    def visit_Assign(self, node: ast.Assign) -> None:
        if any(isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets):
            return
        if all(isinstance(t, (ast.Attribute, ast.Subscript)) for t in node.targets):
            self._drop(node)
        self.generic_visit(node)

    def visit_Expr(self, node: ast.Expr) -> None:
        v = node.value
        if isinstance(v, ast.Constant):  # docstrings and bare literals
            return
        if isinstance(v, ast.Call):
            if _diagnostic(v):
                return
            self._drop(node)
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        if _diagnostic(node):
            return
        f = node.func
        name = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else None
        if name in _SWAP_NAMES and self._one_line(f):
            end = f.end_col_offset
            self._add(f.end_lineno, end - len(name), end, "call", _SWAP_NAMES[name])
        if (
            name in _DROP_CALLS
            and len(node.args) == 1
            and not node.keywords
            and not isinstance(node.args[0], ast.Starred)
            and self._one_line(node)
        ):
            a = node.args[0]
            arg = self._text(a.lineno, a.col_offset, a.end_col_offset)
            self._add(node.lineno, node.col_offset, node.end_col_offset, "call", f"({arg})")
        if name in ("isinstance", "issubclass"):
            # (the types are a set: swapping them changes nothing)
            self.visit(node.func)
            if node.args:
                self.visit(node.args[0])
            return
        self.generic_visit(node)

    # operators

    def visit_Compare(self, node: ast.Compare) -> None:
        left: ast.AST = node.left
        for op, right in zip(node.ops, node.comparators, strict=True):
            old, new = _CMP[type(op)]
            span = self._between(left, right, old)
            if span:
                self._add(*span, "cmp", new)
            left = right
        self.generic_visit(node)

    def visit_BinOp(self, node: ast.BinOp) -> None:
        old, new = _ARITH[type(node.op)]
        span = self._between(node.left, node.right, old)
        if span:
            self._add(*span, "arith", new)
        self.generic_visit(node)

    def visit_AugAssign(self, node: ast.AugAssign) -> None:
        old, new = _ARITH[type(node.op)]
        span = self._between(node.target, node.value, old + "=")
        if span:
            self._add(*span, "arith", new + "=")
        self._drop(node)
        self.generic_visit(node)

    def visit_BoolOp(self, node: ast.BoolOp) -> None:
        old, new = ("and", "or") if isinstance(node.op, ast.And) else ("or", "and")
        for a, b in zip(node.values, node.values[1:], strict=False):
            span = self._between(a, b, old)
            if span:
                self._add(*span, "bool", new)
                break
        self.generic_visit(node)

    def visit_UnaryOp(self, node: ast.UnaryOp) -> None:
        if isinstance(node.op, ast.Not) and node.lineno == node.operand.lineno:
            line = self.lines[node.lineno - 1]
            end = node.col_offset + 3
            if line[node.col_offset : end] == b"not":
                while end < len(line) and line[end : end + 1] == b" ":
                    end += 1
                self._add(node.lineno, node.col_offset, end, "bool", "")
        self.generic_visit(node)

    def visit_Constant(self, node: ast.Constant) -> None:
        v = node.value
        if not self._one_line(node):
            return
        if isinstance(v, bool):
            self._add(node.lineno, node.col_offset, node.end_col_offset, "bool", str(not v))
        elif isinstance(v, int):
            text = self._text(node.lineno, node.col_offset, node.end_col_offset)
            new = hex(v + 1) if text.lower().startswith("0x") else str(v + 1)
            self._add(node.lineno, node.col_offset, node.end_col_offset, "num", new)
        elif isinstance(v, float):
            new = repr(2 * v) if v != 0 else "1.0"
            self._add(node.lineno, node.col_offset, node.end_col_offset, "num", new)

    def visit_Tuple(self, node: ast.Tuple | ast.List) -> None:
        elts = node.elts
        if 2 <= len(elts) <= 4 and self._one_line(node) and all(map(_simple, elts)):
            texts = [self._text(e.lineno, e.col_offset, e.end_col_offset) for e in elts]
            for i in range(len(elts) - 1):
                if texts[i] != texts[i + 1]:
                    a, b = elts[i], elts[i + 1]
                    mid = self._text(a.lineno, a.end_col_offset, b.col_offset)
                    new = texts[i + 1] + mid + texts[i]
                    self._add(node.lineno, a.col_offset, b.end_col_offset, "swap", new)
        self.generic_visit(node)

    visit_List = visit_Tuple

    def visit_ExceptHandler(self, node: ast.ExceptHandler) -> None:
        # (the exception types are a set: swapping them changes nothing)
        for s in node.body:
            self.visit(s)

    # statements

    def visit_Return(self, node: ast.Return) -> None:
        self._drop(node)
        self.generic_visit(node)

    def visit_Raise(self, node: ast.Raise) -> None:
        self._drop(node)
        self.generic_visit(node)

    def visit_Continue(self, node: ast.Continue) -> None:
        self._drop(node)

    def visit_Break(self, node: ast.Break) -> None:
        self._drop(node)

    def _drop(self, node: ast.stmt) -> None:
        if not self._one_line(node):
            return
        line = self.lines[node.lineno - 1]
        before = line[: node.col_offset].strip()
        after = line[node.end_col_offset :].strip()
        if before or (after and not after.startswith(b"#")):
            return
        self._add(node.lineno, node.col_offset, node.end_col_offset, "stmt", "pass")

    def _negate(self, test: ast.expr, stmt: ast.stmt) -> None:
        if not self._one_line(test) or test.lineno != stmt.lineno:
            return
        text = self._text(test.lineno, test.col_offset, test.end_col_offset)
        self._add(test.lineno, test.col_offset, test.end_col_offset, "cond", f"not ({text})")


def _simple(e: ast.AST) -> bool:
    if isinstance(e, ast.Constant):
        return isinstance(e.value, (int, float)) and not isinstance(e.value, bool)
    if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub):
        return _simple(e.operand)
    return isinstance(e, ast.Name)


def _diagnostic(call: ast.Call) -> bool:
    f = call.func
    if isinstance(f, ast.Name):
        return f.id == "print"
    if isinstance(f, ast.Attribute):
        if f.attr == "flush":
            return True
        if f.attr == "write" and isinstance(f.value, ast.Attribute):
            # sys.stdout.write / sys.stderr.write
            return f.value.attr in ("stdout", "stderr")
        if f.attr in _DIAGNOSTIC_CALLS:
            base = f.value
            name = base.id if isinstance(base, ast.Name) else getattr(base, "attr", "")
            return name in _LOGGER_NAMES
    return False


def find_sites(text: str, ops: Iterable[str] = OPERATORS) -> list[Site]:
    """the candidate sites of a module's text, sorted by line and column (not yet checked:
    see ``valid``)"""
    tree = ast.parse(text)
    f = _Finder(text.encode("utf-8").split(b"\n"))
    f.visit(tree)
    wanted = set(ops)
    seen = set()
    out = []
    for s in sorted(f.sites, key=lambda s: (s.line, s.col, s.end, s.op, s.new)):
        key = (s.line, s.col, s.end, s.new)
        if s.op in wanted and key not in seen:
            seen.add(key)
            out.append(s)
    return out


def valid(text: str, site: Site, base_dump: str | None = None) -> bool:
    """the mutated module compiles and its syntax tree differs from the original's"""
    try:
        t = ast.parse(apply_site(text, site))
        compile(t, "<mutant>", "exec")
    except (SyntaxError, ValueError):
        return False
    return ast.dump(t) != (base_dump if base_dump is not None else ast.dump(ast.parse(text)))


def apply_site(text: str, site: Site) -> str:
    """the module's text with the site's replacement (checks that the span holds ``old``)"""
    lines = text.encode("utf-8").split(b"\n")
    line = lines[site.line - 1]
    if line[site.col : site.end].decode("utf-8") != site.old:
        raise ValueError(f"line {site.line}: {site.old!r} is not at columns {site.col}:{site.end}")
    lines[site.line - 1] = line[: site.col] + site.new.encode("utf-8") + line[site.end :]
    return b"\n".join(lines).decode("utf-8")


def apply_mutant(text: str, m: Mutant) -> str:
    if sha256_text(text) != m.sha256:
        raise ValueError(f"{m.id}: {m.file} is not the version the mutant was made from")
    return apply_site(text, Site(m.line, m.col, m.end, m.op, m.old, m.new))


# --- sampling --------------------------------------------------------------------------------


def sample_sites(
    sites: dict[str, list[Site]], texts: dict[str, str], n: int, seed: str
) -> Iterator[tuple[str, Site]]:
    """n valid sites of different lines: lines drawn uniformly from all lines that have sites
    (so dense code is not favoured), then a valid site of the line drawn uniformly (a line
    without one is passed over); deterministic for a seed"""
    rng = random.Random(seed)
    lines = sorted({(f, s.line) for f, ss in sites.items() for s in ss})
    rng.shuffle(lines)
    dumps: dict[str, str] = {}
    k = 0
    for f, ln in lines:
        if k >= n:
            return
        cands = [s for s in sites[f] if s.line == ln]
        rng.shuffle(cands)
        if f not in dumps:
            dumps[f] = ast.dump(ast.parse(texts[f]))
        for s in cands:
            if valid(texts[f], s, dumps[f]):
                k += 1
                yield f, s
                break


def make_mutant(mid: str, area: str, rel: str, text: str, site: Site) -> Mutant:
    mutated = apply_site(text, site).split("\n")[site.line - 1]
    return Mutant(
        id=mid,
        area=area,
        file=rel,
        sha256=sha256_text(text),
        line=site.line,
        col=site.col,
        end=site.end,
        op=site.op,
        old=site.old,
        new=site.new,
        src_line=text.split("\n")[site.line - 1],
        mut_line=mutated,
    )
