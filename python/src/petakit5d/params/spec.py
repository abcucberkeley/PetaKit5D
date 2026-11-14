"""Function specs: the declaration of a user-facing function and its parameters.

A spec is written next to its function, in the function's own area folder::

    from petakit5d.params import Param, register

    @register(
        "XR_crop_dataset",
        Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset directories."),
        Param("inputBbox", "@isnumeric", "Bounding box [ymin, xmin, zmin, ymax, xmax, zmax]."),
        Param("resultDirName", "@ischar", "Result directory under each data path."),
        Param("pad", "@islogical", "Pad where the box leaves the image."),
    )
    def crop_dataset(data_paths, input_bbox, *, result_dir_name="Cropped", pad=False):
        ...

The function's signature gives the Python names and defaults (required arguments are
positional, MATLAB's optional arguments are positional with a default, name-value
parameters are keyword-only); the ``Param`` entries give, in the same order, each MATLAB
name, the MATLAB validator text (copied from the inputParser block) and the meaning. From
these come the snake_case API (the decorated function, which checks names and values as
inputParser does), ``petakit5d.matlab.<MATLAB name>`` (MATLAB names in any case,
name-value pairs, structs), the command-line conversion of MATLAB-style strings (as the mcc
parser functions convert them), the parameters struct MATLAB saves (``matlab_results``) and
the mapping table. ``tools/check_registry.py`` compares every spec with the MATLAB tree.
"""

from __future__ import annotations

import functools
import re
from collections.abc import Callable, Mapping, Sequence
from typing import Any

from petakit5d.params import values as _v

__all__ = [
    "Param",
    "FunctionSpec",
    "ParamSpec",
    "RegistryError",
    "MatlabArgumentError",
    "snake_name",
]


class _Missing:
    """No default (a required argument)."""

    __slots__ = ()

    def __repr__(self) -> str:
        return "<required>"


_MISSING = _Missing()


class RegistryError(Exception):
    """A spec that does not fit its function (a programming error in the port)."""


class MatlabArgumentError(TypeError):
    """Bad arguments, as inputParser rejects them; ``identifier`` is MATLAB's error id."""

    def __init__(self, message: str, identifier: str):
        super().__init__(message)
        self.identifier = identifier


# ----------------------------------------------------------------------------------------------
# names


_TOKEN_RE = re.compile(r"[A-Z]+(?![a-z])|[A-Z]?[a-z]+|[0-9]+[A-Z]+(?![a-z])|[0-9]+[a-z]*")


def snake_name(matlab_name: str) -> str:
    """The snake_case name of a MATLAB name: xyPixelSize -> xy_pixel_size,
    DSRDirName -> dsr_dir_name, save16bit -> save_16bit, save3DStack -> save_3d_stack."""
    parts: list[str] = []
    for chunk in matlab_name.split("_"):
        parts.extend(t.lower() for t in _TOKEN_RE.findall(chunk))
    return "_".join(p for p in parts if p)


# ----------------------------------------------------------------------------------------------
# validators: MATLAB's validator text read as alternatives of class predicates


_PRED_RE = re.compile(r"^(is\w+)\s*\(\s*x\s*\)$")
_ISA_RE = re.compile(r"^isa\s*\(\s*x\s*,\s*['\"](\w+)['\"]\s*\)$")
_KNOWN_PREDS = {
    "ischar",
    "isstr",
    "isstring",
    "islogical",
    "isnumeric",
    "iscell",
    "iscellstr",
    "isstruct",
    "isempty",
    "isscalar",
    "isvector",
    "isfloat",
    "isinteger",
    "isreal",
    "function_handle",
}


def _split_top(text: str, ops: tuple[str, ...]) -> list[str]:
    out, depth, start, i = [], 0, 0, 0
    while i < len(text):
        c = text[i]
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == "'":
            j = text.find("'", i + 1)
            i = len(text) if j < 0 else j
        elif depth == 0:
            for op in ops:
                if text.startswith(op, i):
                    out.append(text[start:i])
                    start = i + len(op)
                    i += len(op) - 1
                    break
        i += 1
    out.append(text[start:])
    return [s.strip() for s in out]


def _strip_parens(t: str) -> str:
    t = t.strip()
    while t.startswith("(") and t.endswith(")"):
        inner = t[1:-1]
        depth = 0
        ok = True
        for c in inner:
            depth += c == "("
            depth -= c == ")"
            if depth < 0:
                ok = False
                break
        if not ok:
            break
        t = inner.strip()
    return t


def validator_alternatives(validator: str | None) -> list[frozenset[str]] | None:
    """The validator as alternatives, each a set of class predicates that must all hold.

    ``@(x) ischar(x) || iscell(x)`` -> [{ischar}, {iscell}]; ``@isnumeric`` -> [{isnumeric}];
    ``@(x) isnumeric(x) && numel(x) == 3`` -> [{isnumeric}] (other terms are not checked);
    None means no check (no validator, or one this reading does not understand).
    """
    if not validator:
        return None
    v = validator.strip()
    if not v.startswith("@"):
        return _declaration_alternatives(v)
    m = re.match(r"^@\s*(\w+)$", v)
    if m:
        name = m.group(1)
        return [frozenset({name})] if name in _KNOWN_PREDS else None
    m = re.match(r"^@\s*\(\s*x\s*\)\s*(.+)$", v, re.S)
    if not m:
        return None
    alts: list[frozenset[str]] = []
    for alt in _split_top(m.group(1), ("||", "|")):
        preds = set()
        for term in _split_top(_strip_parens(alt), ("&&", "&")):
            term = _strip_parens(term)
            pm = _PRED_RE.match(term)
            if pm and pm.group(1) in _KNOWN_PREDS:
                preds.add(pm.group(1))
                continue
            am = _ISA_RE.match(term)
            if am and am.group(1) == "function_handle":
                preds.add("function_handle")
            # negations and other terms (numel, ismember, ...) are not checked
        if not preds:
            return None  # an alternative we cannot check accepts anything
        alts.append(frozenset(preds))
    return alts or None


# arguments-block classes (MATLAB converts to them, so the classes it accepts are wider)
_DECL_CLASS = {
    "char": [{"ischar"}, {"isstring"}],
    "string": [{"ischar"}, {"isstring"}],
    "logical": [{"islogical"}, {"isnumeric"}],
    "cell": [{"iscell"}],
    "struct": [{"isstruct"}],
    "function_handle": [{"function_handle"}],
    **{c: [{"isnumeric"}, {"islogical"}] for c in ("double", "single", *_v.INT_RANGES)},
}
_DECL_MUST = {
    "mustBeNumeric": [{"isnumeric"}],
    "mustBeNumericOrLogical": [{"isnumeric"}, {"islogical"}],
    "mustBeText": [{"ischar"}, {"isstring"}, {"iscellstr"}],
    "mustBeTextScalar": [{"ischar"}, {"isstring"}],
}


def _declaration_alternatives(decl: str) -> list[frozenset[str]] | None:
    """Alternatives of an arguments-block declaration such as ``(1, :) {mustBeNumeric}``
    (size, then class, then validation functions; sizes and other functions not checked)."""
    rest = re.sub(r"^\([^)]*\)\s*", "", decl.strip())
    m = re.match(r"^(\w+)?\s*(?:\{(.*)\})?\s*$", rest, re.S)
    if not m:
        return None
    sets: list[list[set[str]]] = []
    if m.group(1):
        if m.group(1) not in _DECL_CLASS:
            return None
        sets.append(_DECL_CLASS[m.group(1)])
    for fn in re.findall(r"\b(mustBe\w+)", m.group(2) or ""):
        if fn in _DECL_MUST:
            sets.append(_DECL_MUST[fn])
    if not sets:
        return None
    # every set must hold: the alternatives are the combinations
    alts: list[set[str]] = [set()]
    for options in sets:
        alts = [a | o for a in alts for o in options]
    return [frozenset(a) for a in alts]


def mcc_mode(validator: str | None) -> str:
    """The conversion the mcc parser generator writes for a validator (``none`` if none).

    The rules of generate_function_mcc_parser_file.m (check_and_update_actual_argment_type),
    on the text of the function's own validator: "str2num", "eval:<leading chars that are
    evaluated>" or "none" (the string stays text).
    """
    if not validator:
        return "none"
    v = validator
    t = [
        "ischar" in v or "isstring" in v,
        "isempty" in v,
        any(s in v for s in ("isnumeric", "islogical", "isvector", "isscalar")),
        "iscell" in v,
        "isstruct" in v,
        "function_handle" in v,
    ]
    if sum(t) == 0 or (sum(t) == 1 and t[0]) or (t[0] and any(t[1:3]) and not any(t[3:6])):
        return "none"
    if sum(t[1:6]) == 1:
        mode = t[1:6].index(True) + 2
    else:
        mode = 7 if (t[1] and t[2]) else 8
    return {
        2: "eval:[{",
        3: "str2num",
        4: "eval:{",
        5: "eval:[",
        6: "eval:@",
        7: "str2num",
        8: "eval:@[{",
    }[mode]


def _is_num(x: Any) -> bool:
    return isinstance(x, int | float) and not isinstance(x, bool)


def _np_kind(x: Any) -> str | None:
    if type(x).__module__ == "numpy" and hasattr(x, "dtype"):
        return x.dtype.kind
    return None


def _seq_all(x: Any, test: Callable[[Any], bool]) -> bool:
    if not isinstance(x, list | tuple):
        return False
    flat = []
    for e in x:
        if isinstance(e, list | tuple):
            flat.extend(e)
            if any(isinstance(f, list | tuple) for f in e):
                return False
        else:
            flat.append(e)
    return all(test(e) for e in flat)


def _pred(name: str, x: Any) -> bool:
    k = _np_kind(x)
    mclass = _v.matlab_class(x) if isinstance(x, _v.MArray | _v.Cell | _v.MString) else None
    if name in ("ischar", "isstr"):
        if mclass is not None:
            return mclass == "char"
        return isinstance(x, str) or k == "U"
    if name == "isstring":
        return isinstance(x, str) and (mclass in (None, "string"))
    if name == "islogical":
        if mclass is not None:
            return mclass == "logical"
        if k is not None:
            return k == "b"
        return isinstance(x, bool) or (bool(x) and _seq_all(x, lambda e: isinstance(e, bool)))
    if name == "isnumeric":
        if mclass is not None:
            return mclass in _v.NUMERIC_CLASSES
        if k is not None:
            return k in "iuf"
        if x is None or _is_num(x):
            return True
        return _seq_all(x, lambda e: isinstance(e, bool) or _is_num(e)) and not (
            x and _seq_all(x, lambda e: isinstance(e, bool))
        )
    if name == "isfloat":
        if mclass is not None:
            return mclass in _v.FLOAT_CLASSES
        if k is not None:
            return k == "f"
        return x is None or _is_num(x) or (_pred("isnumeric", x) and isinstance(x, list | tuple))
    if name == "isinteger":
        if mclass is not None:
            return mclass in _v.INT_RANGES
        return k in ("i", "u")
    if name == "isreal":
        return not isinstance(x, complex) and k != "c"
    if name == "iscell":
        return mclass == "cell" if mclass is not None else isinstance(x, list | tuple)
    if name == "iscellstr":
        if isinstance(x, _v.Cell):
            return all(_pred("ischar", e) for e in x)
        return isinstance(x, list | tuple) and all(isinstance(e, str) for e in x)
    if name == "isstruct":
        return isinstance(x, Mapping)
    if name == "function_handle":
        return callable(x) and not isinstance(x, type)
    if name == "isempty":
        return _numel(x) == 0
    if name == "isscalar":
        return _numel(x) == 1
    if name == "isvector":
        sz = _size(x)
        return len(sz) == 2 and (sz[0] == 1 or sz[1] == 1) and sz[0] * sz[1] >= 1
    raise RegistryError(f"unknown predicate {name}")


def _size(x: Any) -> tuple[int, ...]:
    if x is None:
        return (0, 0)
    if hasattr(x, "shape") and not isinstance(x, _v.MArray | _v.Cell):
        shp = tuple(int(n) for n in x.shape)
        return shp if len(shp) >= 2 else (1, shp[0]) if shp else (1, 1)
    if isinstance(x, Mapping) or (callable(x) and not isinstance(x, str)):
        return (1, 1)
    if isinstance(x, list | tuple) and x and all(isinstance(e, list | tuple) for e in x):
        return (len(x), len(x[0]))
    try:
        return _v.matlab_size(x)
    except TypeError:
        return (1, 1)


def _numel(x: Any) -> int:
    n = 1
    for d in _size(x):
        n *= d
    return n


def normalize(x: Any) -> Any:
    """Values as the snake_case functions receive them: model values to plain Python
    (``values.to_python``), lists to tuples (deeply); NumPy arrays and others unchanged."""
    if isinstance(x, _v.MArray | _v.Cell | _v.MString):
        return _v.to_python(x)
    if isinstance(x, list | tuple):
        return tuple(normalize(e) for e in x)
    if isinstance(x, dict) and any(
        isinstance(e, _v.MArray | _v.Cell | _v.MString | list) for e in x.values()
    ):
        return {k: normalize(e) for k, e in x.items()}
    return x


# ----------------------------------------------------------------------------------------------
# the declaration


class _Record:
    """A small record class (dataclasses would cost the import of inspect at start-up)."""

    __slots__ = ()
    _fields: tuple[str, ...] = ()
    _defaults: dict[str, Any] = {}

    def __init__(self, *args: Any, **kwargs: Any):
        cls = type(self).__name__
        if len(args) > len(self._fields):
            raise TypeError(f"{cls} takes at most {len(self._fields)} values")
        for name, val in zip(self._fields, args, strict=False):
            if name in kwargs:
                raise TypeError(f"{cls} got two values for {name!r}")
            kwargs[name] = val
        for name in self._fields:
            if name in kwargs:
                val = kwargs.pop(name)
            elif name in self._defaults:
                val = self._defaults[name]
            else:
                raise TypeError(f"{cls} needs {name!r}")
            object.__setattr__(self, name, val)
        if kwargs:
            raise TypeError(f"{cls} got unknown {sorted(kwargs)}")

    def replace(self, **changes: Any) -> _Record:
        return type(self)(**{**{f: getattr(self, f) for f in self._fields}, **changes})

    def __eq__(self, other: object) -> bool:
        if type(other) is not type(self):
            return NotImplemented
        return all(getattr(self, f) == getattr(other, f) for f in self._fields)

    def __repr__(self) -> str:
        shown = ", ".join(f"{f}={getattr(self, f)!r}" for f in self._fields if f != "alternatives")
        return f"{type(self).__name__}({shown})"


class Param(_Record):
    """One parameter of a spec, in the order of the function's signature.

    matlab: the MATLAB name (as in the inputParser block).
    validator: the MATLAB validator text, verbatim (``"@islogical"``); None if MATLAB has none.
    doc: what it means (from the MATLAB function's header).
    py: the Python name if it is not ``snake_name(matlab)``.
    matlab_default: the MATLAB default text when the Python default cannot show it (a default
        computed in MATLAB, like ``pwd``); otherwise the defaults are compared by value.
    mcc: the conversion of command-line strings when MATLAB's parser file differs from the
        generator's rule ("none", "str2num", "eval:{", "eval:[{", ...).
    python_only: a parameter MATLAB does not have (check_registry skips it; the table says so).
    differs: why the Python default differs from MATLAB's on purpose (a site-specific MATLAB
        default such as Slurm options); check_registry then skips the default and the table
        shows the reason. Needs the orchestrator's OK.
    convert: how the MATLAB value becomes the Python value and back: a kind of
        ``petakit5d.params.kinds`` (``"bbox"``, ``"size_yxz"``, ``"per_axis_xyz"``, ``"index1"``,
        ``"range1"``, ``"dim"``, ``"flags_yxz"``, ``"axis_string"``, ``"data_order"``, ...;
        ``"same"``: checked, no conversion) or a ``(from_matlab, to_matlab)`` pair. None: the
        name does not suggest a conversion (check_registry asks for one where it does).
    """

    __slots__ = _PARAM_FIELDS = (
        "matlab",
        "validator",
        "doc",
        "py",
        "matlab_default",
        "mcc",
        "python_only",
        "differs",
        "convert",
    )
    _fields = _PARAM_FIELDS
    _defaults = {
        "doc": "",
        "py": None,
        "matlab_default": None,
        "mcc": None,
        "python_only": False,
        "differs": None,
        "convert": None,
    }
    matlab: str
    validator: str | None
    doc: str
    py: str | None
    matlab_default: str | None
    mcc: str | None
    python_only: bool
    differs: str | None
    convert: Any

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("Param is immutable (use .replace)")

    def __hash__(self) -> int:
        return hash(tuple(getattr(self, f) for f in self._fields))


class ParamSpec(_Record):
    """A resolved parameter: the declaration joined with the function's signature.

    kind: required | optional | parameter; default: the Python default (``_MISSING`` for
    required arguments, else in Python's form); alternatives: the validator's class checks
    (validator_alternatives); conv: the conversion kind (``kinds.Kind``) or None.
    """

    __slots__ = _SPEC_FIELDS = (
        "matlab",
        "py",
        "kind",
        "validator",
        "doc",
        "default",
        "matlab_default",
        "mcc",
        "python_only",
        "differs",
        "alternatives",
        "conv",
    )
    _fields = _SPEC_FIELDS
    _defaults = {"alternatives": None, "conv": None}
    matlab: str
    py: str
    kind: str
    validator: str | None
    doc: str
    default: Any
    matlab_default: str | None
    mcc: str
    python_only: bool
    differs: str | None
    alternatives: list[frozenset[str]] | None
    conv: Any

    __hash__ = None  # type: ignore[assignment]

    @property
    def required(self) -> bool:
        return self.kind == "required"

    def cell_hint(self) -> str | None:
        """'cell' when a sequence for this parameter is a MATLAB cell, not a numeric array."""
        alts = self.alternatives or []
        preds = set().union(*alts) if alts else set()
        if "iscell" in preds and not preds & {"isnumeric", "islogical", "isfloat", "isvector"}:
            return "cell"
        if isinstance(self.default, tuple) and any(isinstance(e, str) for e in self.default):
            return "cell"
        if (
            isinstance(self.default, tuple)
            and any(isinstance(e, tuple) for e in self.default)
            and (self.conv is None or self.conv.name != "bboxes")
        ):
            return "cell"
        return None

    def from_matlab(self, value: Any) -> Any:
        """A (normalized) MATLAB value of this parameter in Python's form."""
        return value if self.conv is None else self.conv.from_matlab(value)

    def to_matlab(self, value: Any) -> Any:
        """A Python value of this parameter in MATLAB's form."""
        return value if self.conv is None else self.conv.to_matlab(value)

    def default_text(self) -> str:
        """The MATLAB text of the default (the Python default converted to MATLAB's form)."""
        if self.matlab_default is not None:
            return self.matlab_default
        if self.default is _MISSING:
            return ""
        from petakit5d.params.literal import format_literal

        hint = self.cell_hint() or _cell_if_empty_default(self)
        return format_literal(self.to_matlab(self.default), hint)

    def check(self, value: Any, fname: str) -> None:
        """Raise MatlabArgumentError if ``value`` fails the MATLAB validator's class checks."""
        alts = self.alternatives
        if alts is None:
            return
        for preds in alts:
            if all(_pred(p, value) for p in preds):
                return
        raise MatlabArgumentError(
            f"The value of '{self.matlab}' is invalid. It must satisfy the function: "
            f"{_matlab_validator_display(self.validator)}.",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )


def _kind_of(d: Param, where: str) -> Any:
    if d.convert is None:
        return None
    from petakit5d.params.kinds import get_kind

    try:
        return get_kind(d.convert)
    except ValueError as err:
        raise RegistryError(f"{where}: {d.matlab}: {err}") from None


def _passes(p: ParamSpec, value: Any) -> bool:
    if p.alternatives is None:
        return True
    return any(all(_pred(q, value) for q in alt) for alt in p.alternatives)


def _as_matlab_sees(p: ParamSpec, value: Any) -> Any:
    """A snake_case value as MATLAB's validator sees it: in MATLAB's form by the parameter's
    kind (an open end None is MATLAB's Inf). Image arrays (only transposed) and a function's
    own conversions are checked as given, and so is a value the kind cannot convert."""
    k = p.conv
    if k is None or k.name in ("same", "volume") or k.name.startswith("custom:"):
        return value
    try:
        return k.to_matlab(normalize(value))
    except (TypeError, ValueError):
        return value


def _cell_if_empty_default(p: ParamSpec) -> str | None:
    return "cell" if p.default == () else None


def _matlab_validator_display(validator: str | None) -> str:
    """The validator as inputParser's message shows it (func2str: a named function without
    its @, an anonymous one without blanks outside quotes)."""
    v = (validator or "").strip()
    m = re.fullmatch(r"@\s*(\w+)", v)
    if m:
        return m.group(1)
    parts = re.split(r"('(?:[^']|'')*')", v)
    return "".join(q if i % 2 else re.sub(r"\s+", "", q) for i, q in enumerate(parts))


class FunctionSpec:
    """A registered function: MATLAB name, Python function, parameters."""

    def __init__(
        self,
        matlab_name: str,
        func: Callable[..., Any],
        params: Sequence[Param],
        doc: str | None,
        keep_unmatched: bool,
        task: str | None,
        matlab_in: Callable[[dict[str, Any]], dict[str, Any]] | None = None,
        matlab_out: Callable[[Any], Any] | None = None,
    ):
        self.matlab_name = matlab_name
        self.func = func
        self.declared = tuple(params)
        self._doc = doc
        self.keep_unmatched = keep_unmatched
        self.task = task
        self._matlab_in = matlab_in
        self._matlab_out = matlab_out
        self._params: list[ParamSpec] | None = None
        self.wrapper: Callable[..., Any] | None = None

    # -- resolution (lazy, so that importing an area module stays cheap) ---------------------
    @property
    def python_name(self) -> str:
        return f"{self.func.__module__}.{self.func.__qualname__}"

    @property
    def doc(self) -> str:
        if self._doc is not None:
            return self._doc
        import inspect

        d = inspect.getdoc(self.func) or ""
        return d.split("\n\n")[0].replace("\n", " ").strip()

    def _hook(self, name: str) -> Callable[..., Any] | None:
        own = getattr(self, f"_{name}")
        if own is not None:
            return own
        # function-level hooks a module declares by name (as tools/memory_estimate.py does)
        import sys

        f = getattr(sys.modules.get(self.func.__module__), name, None)
        return f if callable(f) else None

    @property
    def matlab_in(self) -> Callable[[dict[str, Any]], dict[str, Any]] | None:
        """The function-level conversion of the bound keyword arguments, MATLAB -> Python."""
        return self._hook("matlab_in")

    @property
    def matlab_out(self) -> Callable[[Any], Any] | None:
        """The function-level conversion of the result, Python -> MATLAB."""
        return self._hook("matlab_out")

    @property
    def params(self) -> list[ParamSpec]:
        if self._params is None:
            self._params = self._resolve()
        return self._params

    def _resolve(self) -> list[ParamSpec]:
        import inspect

        sig = inspect.signature(self.func)
        sparams = [
            p
            for p in sig.parameters.values()
            if p.kind not in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD)
        ]
        where = f"spec of {self.matlab_name} ({self.python_name})"
        if len(sparams) != len(self.declared):
            raise RegistryError(
                f"{where}: the signature has {len(sparams)} parameters, the spec "
                f"{len(self.declared)}: {[p.name for p in sparams]} vs "
                f"{[d.matlab for d in self.declared]}"
            )
        out: list[ParamSpec] = []
        seen: set[str] = set()
        for sp, d in zip(sparams, self.declared, strict=True):
            default = _MISSING if sp.default is inspect.Parameter.empty else sp.default
            py = d.py or snake_name(d.matlab)
            if sp.name != py:
                raise RegistryError(
                    f"{where}: parameter {sp.name!r} is declared as {d.matlab!r}, whose Python "
                    f"name is {py!r} (give py={sp.name!r} to use another name)"
                )
            if sp.kind == inspect.Parameter.KEYWORD_ONLY:
                kind = "parameter"
                if default is _MISSING:
                    raise RegistryError(f"{where}: keyword-only {sp.name!r} needs a default")
            elif sp.kind == inspect.Parameter.POSITIONAL_ONLY or default is _MISSING:
                kind = "required"
                if out and out[-1].kind != "required":
                    raise RegistryError(f"{where}: required {sp.name!r} after optional ones")
            else:
                kind = "optional"
            key = d.matlab.lower()
            if key in seen:
                raise RegistryError(f"{where}: {d.matlab!r} is declared twice")
            seen.add(key)
            out.append(
                ParamSpec(
                    matlab=d.matlab,
                    py=py,
                    kind=kind,
                    validator=d.validator,
                    doc=d.doc,
                    default=default,
                    matlab_default=d.matlab_default,
                    mcc=d.mcc or mcc_mode(d.validator),
                    python_only=d.python_only,
                    differs=d.differs,
                    alternatives=validator_alternatives(d.validator),
                    conv=_kind_of(d, where),
                )
            )
        return out

    def param(self, name: str) -> ParamSpec:
        """A parameter by Python or MATLAB name (exact)."""
        for p in self.params:
            if name in (p.py, p.matlab):
                return p
        raise KeyError(name)

    # -- name resolution as inputParser does it --------------------------------------------
    def _resolve_matlab_name(self, name: Any) -> ParamSpec | None:
        """The parameter a name-value name means: exact (any case), else a unique prefix,
        among all arguments (required and optional ones can be given by name too)."""
        if isinstance(name, _v.MString):
            name = str(name)
        if not isinstance(name, str):
            raise MatlabArgumentError(
                "Expected a string scalar or character vector for the parameter name.",
                "MATLAB:InputParser:ParamMustBeChar",
            )
        matches = self._matches(name)
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            names = ", ".join(f"'{p.matlab}'" for p in results_order(matches))
            raise MatlabArgumentError(
                f"'{name}' matches multiple parameter names: {names}. To avoid ambiguity, "
                f"specify the complete name of the parameter.",
                "MATLAB:InputParser:AmbiguousParameter",
            )
        return self._unmatched(name)

    def _resolve_field(self, name: str) -> ParamSpec | None:
        """The parameter a struct field means: exact names only (any case)."""
        low = name.lower()
        for p in self.params:
            if p.matlab.lower() == low:
                return p
        return self._unmatched(name)

    def _unmatched(self, name: str) -> None:
        if self.keep_unmatched:
            return None
        raise MatlabArgumentError(
            f"'{name}' is not a recognized parameter. For a list of valid name-value pair "
            f"arguments, see the documentation for this function."
            + _suggest(name, [p.matlab for p in self.params]),
            "MATLAB:InputParser:UnmatchedParameter",
        )

    def _matches(self, name: str) -> list[ParamSpec]:
        low = name.lower()
        exact = [p for p in self.params if p.matlab.lower() == low]
        if exact:
            return exact
        return [p for p in self.params if p.matlab.lower().startswith(low)]

    # -- the three ways in --------------------------------------------------------------------
    def bind_python(self, args: tuple[Any, ...], kwargs: dict[str, Any]) -> dict[str, Any]:
        """Check and normalize a snake_case call; returns {python name: value} of what was given.
        Each value is checked in MATLAB's form, as the validator sees it in MATLAB."""
        params = self.params
        by_py = {p.py: p for p in params}
        for k in kwargs:
            if k not in by_py:
                hint = _suggest(k, list(by_py))
                if not hint:
                    lower = {p.matlab.lower(): p.py for p in params}
                    if k.lower() in lower:
                        hint = f" Did you mean '{lower[k.lower()]}'?"
                raise TypeError(
                    f"{self.func.__qualname__}() got an unexpected keyword argument {k!r}.{hint}"
                )
        positional = [p for p in params if p.kind != "parameter"]
        if len(args) > len(positional):
            raise TypeError(
                f"{self.func.__qualname__}() takes {len(positional)} positional arguments "
                f"but {len(args)} were given"
            )
        given: dict[str, Any] = {}
        for p, a in zip(positional, args, strict=False):
            given[p.py] = a
        for k, a in kwargs.items():
            if k in given:
                raise TypeError(f"{self.func.__qualname__}() got multiple values for {k!r}")
            given[k] = a
        missing = [p.py for p in positional if p.required and p.py not in given]
        if missing:
            raise TypeError(
                f"{self.func.__qualname__}() missing required arguments: "
                + ", ".join(repr(m) for m in missing)
            )
        out = {}
        for k, a in given.items():
            p = by_py[k]
            p.check(_as_matlab_sees(p, a), self.matlab_name)
            out[k] = normalize(a)
        return out

    def bind_matlab(
        self, args: Sequence[Any], kwargs: Mapping[str, Any] | None = None
    ) -> dict[str, Any]:
        """Bind a MATLAB-style call: required, then optional positionals, then name-value
        pairs or structs, then keyword arguments (MATLAB names, any case, unique prefixes).
        Returns {python name: normalized value} of what was given; later values win."""
        given = self._bind(list(args), dict(kwargs or {}))
        return self._from_matlab({k: normalize(a) for k, a in given.items()})

    def _from_matlab(self, values: dict[str, Any]) -> dict[str, Any]:
        by_py = {p.py: p for p in self.params}
        out = {}
        for k, v in values.items():
            p = by_py[k]
            try:
                out[k] = p.from_matlab(v)
            except (ValueError, TypeError) as err:
                raise MatlabArgumentError(
                    f"The value of '{p.matlab}' is not a valid {p.conv.doc or p.conv.name}: {err}",
                    "PETAKIT5D:params:badValue",
                ) from None
        return out

    def bind_strings(self, argv: Sequence[str]) -> dict[str, Any]:
        """Bind command-line strings as the compiled master and the mcc parser functions do.

        mccMaster passes the first arguments to every required and optional argument by
        position, the rest as name-value pairs (names as inputParser matches them). Each
        value is then converted by its parameter's rule (``str2num``, ``eval`` of a leading
        ``{``, ``[`` or ``@``, or kept as text) and checked by the function's own validator.
        """
        argv = [str(a) for a in argv]
        positional = [p for p in self.params if p.kind != "parameter"]
        if len(argv) < len(positional):
            raise MatlabArgumentError(
                f"Index exceeds the number of array elements. Index must not exceed {len(argv)}.",
                "MATLAB:badsubscript",
            )
        raw = {p.py: s for p, s in zip(positional, argv, strict=False)}
        rest = argv[len(positional) :]
        for i in range(0, len(rest), 2):
            q = self._resolve_matlab_name(rest[i])
            if i + 1 >= len(rest):
                raise MatlabArgumentError(
                    f"No value was given for '{rest[i]}'. Name-value pair arguments require a "
                    f"name followed by a value.",
                    "MATLAB:InputParser:ParamMissingValue",
                )
            if q is not None:
                raw[q.py] = rest[i + 1]
        given = {}
        by_py = {p.py: p for p in self.params}
        for k, s in raw.items():
            p = by_py[k]
            val = convert_string(p, s)
            p.check(val, self.matlab_name)
            given[k] = val
        return self._from_matlab({k: normalize(a) for k, a in given.items()})

    def _bind(self, args: list[Any], kwargs: dict[str, Any]) -> dict[str, Any]:
        params = self.params
        required = [p for p in params if p.kind == "required"]
        optional = [p for p in params if p.kind == "optional"]
        if len(args) < len(required):
            raise MatlabArgumentError("Not enough input arguments.", "MATLAB:minrhs")
        given: dict[str, Any] = {}

        def put(p: ParamSpec, val: Any) -> None:
            p.check(val, self.matlab_name)
            given[p.py] = val

        for p, a in zip(required, args, strict=False):
            put(p, a)
        rest = args[len(required) :]
        for p in optional:
            if not rest:
                break
            a = rest[0]
            if p.alternatives is None and not p.validator and isinstance(a, str):
                # an optional without a validator never takes text: inputParser reads it as
                # a parameter name, and fails for text that names no parameter
                if self._matches(str(a)):
                    break
                raise MatlabArgumentError(
                    f"The argument '{a}' is a string scalar or character vector and does not "
                    f"match any parameter names. It failed validation for the argument "
                    f"'{p.matlab}'.",
                    "MATLAB:InputParser:ArgumentFailedDefaultValidation",
                )
            if not _passes(p, a):
                # not this optional: a name followed by a value, or a struct of names,
                # starts the name-value part; anything else fails the optional's validator
                if isinstance(a, str) and len(rest) > 1 and len(self._matches(str(a))) == 1:
                    break
                if isinstance(a, Mapping) and all(
                    any(q.matlab.lower() == str(k).lower() for q in params) for k in a
                ):
                    break
                p.check(a, self.matlab_name)
            given[p.py] = a
            rest = rest[1:]
        i = 0
        while i < len(rest):
            a = rest[i]
            if isinstance(a, Mapping):
                for k, val in a.items():
                    q = self._resolve_field(str(k))
                    if q is not None:
                        put(q, val)
                i += 1
                continue
            q = self._resolve_matlab_name(a)
            if i + 1 >= len(rest):
                raise MatlabArgumentError(
                    f"No value was given for '{a}'. Name-value pair arguments require a name "
                    f"followed by a value.",
                    "MATLAB:InputParser:ParamMissingValue",
                )
            if q is not None:
                put(q, rest[i + 1])
            i += 2
        for k, val in kwargs.items():
            q = self._resolve_matlab_name(k)
            if q is not None:
                put(q, val)
        return given

    # -- calls ------------------------------------------------------------------------------
    def call_python(self, *args: Any, **kwargs: Any) -> Any:
        return self.func(**self.bind_python(args, kwargs))

    def call_matlab(self, *args: Any, **kwargs: Any) -> Any:
        """Call with MATLAB's arguments and values; the result in MATLAB's form."""
        return self._call_converted(self.bind_matlab(args, kwargs))

    def call_strings(self, argv: Sequence[str]) -> Any:
        """Call with command-line strings; the result in MATLAB's form."""
        return self._call_converted(self.bind_strings(argv))

    def _call_converted(self, kwargs: dict[str, Any]) -> Any:
        hook_in, hook_out = self.matlab_in, self.matlab_out
        if hook_in is not None:
            kwargs = hook_in(kwargs)
        result = self.func(**kwargs)
        return result if hook_out is None else hook_out(result)

    # -- what MATLAB saves --------------------------------------------------------------------
    def matlab_results(self, values: Mapping[str, Any]) -> dict[str, Any]:
        """``ip.Results`` as MATLAB builds it: {MATLAB name: model value}, in MATLAB's order.

        ``values`` maps Python names to the values the function got (``locals()`` at its
        top works); parameters missing from it take their defaults. Values are converted to
        MATLAB's form by their kinds (boxes, sizes, axis strings, ...), so the files hold
        MATLAB's values. Python-only parameters are left out. Pass the result to the .mat and
        JSON writers (parameters.mat/.json).
        """
        out = {}
        for p in results_order(self.params):
            if p.python_only:
                continue
            if p.py in values:
                v = values[p.py]
            elif p.default is not _MISSING:
                v = p.default
            else:
                raise KeyError(f"no value for required argument {p.py!r}")
            v = p.to_matlab(normalize(v))
            hint = p.cell_hint() or (_cell_if_empty_default(p) if v == () else None)
            out[p.matlab] = _v.from_python(v, hint)
        return out

    def parameters_json(self, values: Mapping[str, Any]) -> str:
        """``jsonencode(pr, PrettyPrint=true)`` of ``matlab_results(values)``: the text MATLAB's
        writeJsonFile puts in parameters.json (MATLAB's values)."""
        from petakit5d.core.jsonenc import jsonencode

        res = self.matlab_results(values)
        return jsonencode({k: _v.to_plain(v) for k, v in res.items()}, pretty_print=True)

    def __repr__(self) -> str:
        return f"<FunctionSpec {self.matlab_name} -> {self.python_name}>"


def results_order(params: Sequence[ParamSpec]) -> list[ParamSpec]:
    """The field order of inputParser's Results: names sorted ignoring case."""
    return sorted(params, key=lambda p: p.matlab.lower())


def convert_string(p: ParamSpec, s: str) -> Any:
    """A command-line string converted by the parameter's mcc rule (a model value or str)."""
    from petakit5d.params.literal import matlab_eval, str2num

    mode = p.mcc
    if mode == "str2num":
        return str2num(s)[0]
    if mode.startswith("eval:"):
        if s and s[0] in mode[5:]:
            if s[0] == "@":
                raise MatlabArgumentError(
                    f"function handles in strings are not supported ({p.matlab}: {s!r})",
                    "PETAKIT5D:params:functionHandleString",
                )
            return matlab_eval(s)
        return s
    return s


def _suggest(name: str, choices: Sequence[str]) -> str:
    low = {c.lower(): c for c in choices}
    import difflib

    m = difflib.get_close_matches(name.lower(), list(low), n=1, cutoff=0.6)
    return f" Did you mean '{low[m[0]]}'?" if m else ""


def make_wrapper(spec: FunctionSpec) -> Callable[..., Any]:
    func = spec.func

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        return func(**spec.bind_python(args, kwargs))

    wrapper.__pk5d_spec__ = spec  # type: ignore[attr-defined]
    return wrapper
