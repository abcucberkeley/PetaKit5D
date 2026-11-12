"""The registry: every user-facing function by MATLAB name and by Python name.

Specs register themselves when their module is imported (``@register``). To find a
function by MATLAB name without importing the whole package, the registry reads the package's
source files once for ``register("<name>"`` and ``not_ported("<name>"`` and imports only the
module that declares the name. So a spec must name its MATLAB function as the first argument
of ``register``, as a plain string literal (test_registry checks that the scan and the
imports agree).
"""

from __future__ import annotations

import importlib
import os
import re
import sys
from collections.abc import Callable, Iterator
from typing import Any

from petakit5d.params.spec import FunctionSpec, Param, RegistryError, _Record, make_wrapper

__all__ = ["Registry", "NotPorted", "REGISTRY", "register", "not_ported"]

_DECL_RE = re.compile(r"""\b(register|not_ported)\(\s*["']([A-Za-z]\w*)["']""")
_SKIP_DIRS = {"__pycache__", "_ext", "testing", "matlab"}


class NotPorted(_Record):
    """A MATLAB function that is not ported, with the reason, its area and where declared."""

    __slots__ = _fields = ("matlab_name", "reason", "area", "module")
    matlab_name: str
    reason: str
    area: str
    module: str


class Registry:
    """Specs by MATLAB name; ``root`` is the package folder scanned for declarations."""

    def __init__(self, package: str | None = None, root: str | os.PathLike | None = None):
        self.package = package
        self.root = os.fspath(root) if root else None
        self._specs: dict[str, FunctionSpec] = {}
        self._by_python: dict[str, FunctionSpec] = {}
        self._not_ported: dict[str, NotPorted] = {}
        self._index: dict[str, tuple[str, str]] | None = None

    # -- declaring ---------------------------------------------------------------------------
    def register(
        self,
        matlab_name: str,
        *params: Param,
        doc: str | None = None,
        keep_unmatched: bool = False,
        task: str | None = None,
        matlab_in: Callable[[dict[str, Any]], dict[str, Any]] | None = None,
        matlab_out: Callable[[Any], Any] | None = None,
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        """Decorator: declare ``func`` as the port of MATLAB's ``matlab_name``.

        params: one ``Param`` per parameter of the function, in signature order.
        doc: one line for the mapping table (default: the docstring's first paragraph).
        keep_unmatched: MATLAB's ``ip.KeepUnmatched = true`` (unknown names are ignored).
        task: the port task that owns the function (shown in the mapping table).
        matlab_in, matlab_out: function-level conversions for the MATLAB-name layer and the
            command line, applied after the parameters' kinds: ``matlab_in(kwargs) -> kwargs``
            (MATLAB values to Python's) and ``matlab_out(result) -> result`` (Python's to
            MATLAB's). Without them, functions of those names in the function's module are used.
        """
        if not isinstance(matlab_name, str) or not re.fullmatch(r"[A-Za-z]\w*", matlab_name):
            raise RegistryError(f"bad MATLAB name {matlab_name!r}")
        for p in params:
            if not isinstance(p, Param):
                raise RegistryError(f"{matlab_name}: parameters are Param(...), got {p!r}")

        def deco(func: Callable[..., Any]) -> Callable[..., Any]:
            spec = FunctionSpec(
                matlab_name, func, params, doc, keep_unmatched, task, matlab_in, matlab_out
            )
            old = self._specs.get(matlab_name)
            if old is not None and old.python_name != spec.python_name:
                raise RegistryError(
                    f"{matlab_name} is registered twice: {old.python_name} and {spec.python_name}"
                )
            if matlab_name in self._not_ported:
                raise RegistryError(f"{matlab_name} is registered and declared not ported")
            wrapper = make_wrapper(spec)
            spec.wrapper = wrapper
            self._specs[matlab_name] = spec
            self._by_python[spec.python_name] = spec
            return wrapper

        return deco

    def not_ported(self, matlab_name: str, reason: str, area: str = "", module: str = "") -> None:
        """Declare a MATLAB function that is not ported, with the reason (for the table)."""
        if matlab_name in self._specs:
            raise RegistryError(f"{matlab_name} is registered and declared not ported")
        self._not_ported[matlab_name] = NotPorted(matlab_name, reason, area, module)

    # -- finding -----------------------------------------------------------------------------
    def index(self) -> dict[str, tuple[str, str]]:
        """{MATLAB name: (module, "register" | "not_ported")} of every declaration in the
        package's sources (nothing imported)."""
        if self._index is None:
            self._index = scan_declarations(self.root, self.package) if self.root else {}
        return self._index

    def find(self, matlab_name: str) -> FunctionSpec:
        """The spec of a MATLAB function (imports its module on first use)."""
        spec = self._specs.get(matlab_name)
        if spec is not None:
            return spec
        entry = self.index().get(matlab_name)
        if entry is not None:
            importlib.import_module(entry[0])
            spec = self._specs.get(matlab_name)
            if spec is not None:
                return spec
        if matlab_name in self._not_ported:
            raise LookupError(
                f"{matlab_name} is not ported: {self._not_ported[matlab_name].reason}"
            )
        if self.root is not None:
            from petakit5d.params.planned import planned

            plan = planned().get(matlab_name)
            if plan is not None:
                raise LookupError(
                    f"{matlab_name} is not ported yet (planned as {plan[0]}, task {plan[1]})"
                )
        import difflib

        names = sorted(set(self.index()) | set(self._specs))
        m = difflib.get_close_matches(matlab_name, names, n=1, cutoff=0.6)
        hint = f" Did you mean {m[0]}?" if m else ""
        raise LookupError(f"no registered function {matlab_name}.{hint}")

    def find_python(self, name: str | Callable[..., Any]) -> FunctionSpec:
        """The spec of a Python function: the function itself, 'module.function' or a name."""
        if callable(name):
            spec = getattr(name, "__pk5d_spec__", None)
            if spec is None:
                raise LookupError(f"{name!r} is not a registered function")
            return spec
        if name in self._by_python:
            return self._by_python[name]
        self.load_all()
        if name in self._by_python:
            return self._by_python[name]
        for spec in self._specs.values():
            if spec.python_name.rsplit(".", 1)[-1] == name:
                return spec
        raise LookupError(f"no registered Python function {name}")

    def load_all(self) -> None:
        """Import every module that declares something."""
        for mod in sorted({m for m, _ in self.index().values()}):
            importlib.import_module(mod)

    def specs(self) -> list[FunctionSpec]:
        """Every registered spec (all declaring modules imported), by MATLAB name."""
        self.load_all()
        return [self._specs[k] for k in sorted(self._specs, key=str.lower)]

    def not_ported_all(self) -> dict[str, NotPorted]:
        self.load_all()
        return dict(sorted(self._not_ported.items(), key=lambda kv: kv[0].lower()))

    def names(self) -> list[str]:
        """MATLAB names of the registered functions (from the scan; nothing imported)."""
        scanned = {k for k, (_, kind) in self.index().items() if kind == "register"}
        return sorted(scanned | set(self._specs), key=str.lower)

    def __contains__(self, matlab_name: str) -> bool:
        entry = self.index().get(matlab_name)
        return matlab_name in self._specs or (entry is not None and entry[1] == "register")

    def __iter__(self) -> Iterator[FunctionSpec]:
        return iter(self.specs())


def scan_declarations(root: str, package: str | None) -> dict[str, tuple[str, str]]:
    """Read the package's .py files for ``register("X"`` and ``not_ported("X"``.

    The registry's own package is skipped (its docstrings show examples), except the module
    that declares the functions without an area (params/unported.py).
    """
    out: dict[str, tuple[str, str]] = {}
    root = os.fspath(root)
    pkg = package or os.path.basename(root)
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d not in _SKIP_DIRS)
        rel = os.path.relpath(dirpath, root)
        prefix = [pkg] + ([] if rel == "." else rel.replace(os.sep, "/").split("/"))
        for fn in sorted(filenames):
            if not fn.endswith(".py"):
                continue
            parts = prefix if fn == "__init__.py" else [*prefix, fn[:-3]]
            mod = ".".join(parts)
            if mod.startswith(f"{pkg}.params") and mod != f"{pkg}.params.unported":
                continue
            try:
                with open(os.path.join(dirpath, fn), encoding="utf-8") as f:
                    text = f.read()
            except OSError:
                continue
            if "register(" not in text and "not_ported(" not in text:
                continue
            for m in _DECL_RE.finditer(text):
                out.setdefault(m.group(2), (mod, m.group(1)))
    return out


REGISTRY = Registry("petakit5d", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def register(
    matlab_name: str,
    *params: Param,
    doc: str | None = None,
    keep_unmatched: bool = False,
    task: str | None = None,
    matlab_in: Callable[[dict[str, Any]], dict[str, Any]] | None = None,
    matlab_out: Callable[[Any], Any] | None = None,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Declare a user-facing function in the package registry (see ``Registry.register``)."""
    return REGISTRY.register(
        matlab_name,
        *params,
        doc=doc,
        keep_unmatched=keep_unmatched,
        task=task,
        matlab_in=matlab_in,
        matlab_out=matlab_out,
    )


def not_ported(matlab_name: str, reason: str, area: str = "") -> None:
    """Declare a MATLAB function that is not ported (see ``Registry.not_ported``)."""
    module = sys._getframe(1).f_globals.get("__name__", "")
    REGISTRY.not_ported(matlab_name, reason, area, module)


def find(matlab_name: str) -> FunctionSpec:
    """The spec of a registered MATLAB function (``Registry.find``)."""
    return REGISTRY.find(matlab_name)


def spec_of(func: Callable[..., Any] | str) -> FunctionSpec:
    """The spec of a registered function (the function, 'module.name' or a MATLAB name)."""
    if isinstance(func, str) and func in REGISTRY:
        return REGISTRY.find(func)
    return REGISTRY.find_python(func)


def matlab_results(func: Callable[..., Any] | str, values: dict[str, Any]) -> dict[str, Any]:
    """``ip.Results`` of a call, as MATLAB saves it in parameters.mat/.json.

    Inside a registered function, ``matlab_results(crop_dataset, locals())`` at the top.
    """
    return spec_of(func).matlab_results(values)
