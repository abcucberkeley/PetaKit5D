"""PetaKit5D: petabyte-scale 5D image processing for light-sheet microscopy.

Importing the package loads nothing else: NumPy, the I/O packages and the compiled
extension modules are imported by the functions that need them. Subpackages load on first
use (``petakit5d.io``, ``petakit5d.deskew_rotate``, ...).
"""

from petakit5d._version import __version__

__all__ = ["__version__"]


def __getattr__(name: str):
    # subpackages on first attribute access, so that "import petakit5d" stays cheap
    if name.startswith("_"):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    try:
        return importlib.import_module(f"{__name__}.{name}")
    except ModuleNotFoundError as err:
        if err.name == f"{__name__}.{name}":
            raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
        raise
