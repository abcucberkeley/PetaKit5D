"""PetaKit5D functions under their MATLAB names, with MATLAB's calling conventions.

    from petakit5d import matlab
    matlab.XR_crop_dataset("/mnt/projects/embryo_2024", [1, 1, 1, 100, 100, 50],
                           "resultDirName", "Cropped", BatchSize=[512, 512, 512])

Required and optional arguments come first, then name-value pairs (MATLAB names in any case,
unique prefixes, structs as dicts) or keyword arguments, with inputParser's errors for
unknown, ambiguous and invalid ones. Each function calls its snake_case port (see the
mapping table, docs/MATLAB_MAPPING.md).
"""

from __future__ import annotations

from typing import Any


class MatlabFunction:
    """A registered function called the MATLAB way."""

    __slots__ = ("spec",)

    def __init__(self, spec: Any):
        self.spec = spec

    @property
    def __name__(self) -> str:
        return self.spec.matlab_name

    @property
    def __doc__(self) -> str:  # type: ignore[override]
        s = self.spec
        lines = [f"{s.matlab_name}: MATLAB-style call of {s.python_name}.", "", s.doc, ""]
        for p in s.params:
            dflt = "" if p.required else f" (default: {p.default_text()})"
            lines.append(f"{p.matlab}{dflt}: {p.doc}")
        return "\n".join(lines)

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        return self.spec.call_matlab(*args, **kwargs)

    def __repr__(self) -> str:
        return f"<MATLAB-name function {self.spec.matlab_name} -> {self.spec.python_name}>"


def __getattr__(name: str) -> MatlabFunction:
    if name.startswith("__"):
        raise AttributeError(name)
    from petakit5d.params.registry import REGISTRY

    try:
        spec = REGISTRY.find(name)
    except LookupError as err:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}: {err}") from None
    return MatlabFunction(spec)


def __dir__() -> list[str]:
    from petakit5d.params.registry import REGISTRY

    return REGISTRY.names()
