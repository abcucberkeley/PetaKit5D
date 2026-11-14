"""PetaKit5D: petabyte-scale 5D image processing for light-sheet microscopy.

Importing the package loads nothing else: NumPy, the I/O packages and the compiled
extension modules are imported by the functions that need them. Subpackages load on first
use (``petakit5d.io``, ``petakit5d.deskew_rotate``, ...).

It sets one variable of the process's environment, before anything loads an OpenMP runtime:
``GOMP_SPINCOUNT=100000``, unless ``GOMP_SPINCOUNT`` or ``OMP_WAIT_POLICY`` is set already, and
notes that it did in ``PETAKIT5D_GOMP_SPINCOUNT`` (``docs/decisions/0007-openmp-spin.md``).
"""

import os as _os

from petakit5d._version import __version__

__all__ = ["__version__"]

# After each parallel region a libgomp team spins GOMP_SPINCOUNT times before it sleeps
# (300,000 by default: about 11 ms on Skylake servers, 3 ms on Ice Lake ones). A port process
# holds two libgomp runtimes, the one the cpp-zarr and cpp-tiff wheels share and the one of the
# port's modules, and each counts only its own threads, so the idle team of one spins on the
# cores the other's work needs. 100,000 (about 3.6 and 1 ms) still bridges the gaps between
# calls, which a sleeping team needs up to 0.5 ms to cross. libgomp reads the variable once,
# when it is loaded, so it is set here, before any extension module or I/O package is
# imported; a value the user set wins (OMP_WAIT_POLICY implies its own spin). The mark tells
# child processes that the value is the port's, so that programs that aren't the port (MATLAB,
# whose mex files use libgomp too) can be started without it.
_OPENMP_SPINCOUNT = "100000"
_OPENMP_MARK = "PETAKIT5D_GOMP_SPINCOUNT"


def _set_openmp_spin(environ=_os.environ) -> bool:
    """Set GOMP_SPINCOUNT, and the mark, in ``environ`` unless GOMP_SPINCOUNT or
    OMP_WAIT_POLICY is in it; True when it was set here."""
    if "GOMP_SPINCOUNT" in environ or "OMP_WAIT_POLICY" in environ:
        return False
    environ["GOMP_SPINCOUNT"] = environ[_OPENMP_MARK] = _OPENMP_SPINCOUNT
    return True


def _without_openmp_spin(environ) -> dict:
    """A copy of ``environ`` without the GOMP_SPINCOUNT that this package set (in this process
    or a parent), as the user's shell would start a program; the user's own value stays."""
    env = dict(environ)
    mark = env.pop(_OPENMP_MARK, None)
    if mark is not None and env.get("GOMP_SPINCOUNT") == mark:
        del env["GOMP_SPINCOUNT"]
    return env


_set_openmp_spin()


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
