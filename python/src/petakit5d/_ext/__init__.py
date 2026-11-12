"""The compiled extension modules: ``petakit5d._ext.<group>``, one per folder of ``python/cpp``.

Installed modules sit in this folder. A build tree outside the package
(``cmake -S python -B <build>``) is used by setting ``PETAKIT5D_EXT_DIR=<build>/petakit5d/_ext``
(several folders separated by ``os.pathsep``); it is searched before this folder.
"""

import os
from importlib import import_module
from importlib.machinery import EXTENSION_SUFFIXES
from importlib.util import find_spec
from pathlib import Path
from types import ModuleType

__all__ = ["available", "cuda_source_dirs", "groups", "load"]

_env_dirs = os.environ.get("PETAKIT5D_EXT_DIR", "")
__path__ = [os.path.abspath(p) for p in _env_dirs.split(os.pathsep) if p] + list(__path__)


def load(group: str) -> ModuleType:
    """Import the extension module of a kernel group, with a hint when it is not built."""
    try:
        return import_module(f"{__name__}.{group}")
    except ModuleNotFoundError as err:
        if err.name != f"{__name__}.{group}":
            raise
        raise ModuleNotFoundError(
            f"petakit5d extension module {group!r} is not built (searched {__path__}); "
            "build with 'pip install .' or set PETAKIT5D_EXT_DIR to a build tree's "
            "petakit5d/_ext folder",
            name=err.name,
        ) from None


def available(group: str) -> bool:
    """Whether the extension module of a group can be imported (without importing it)."""
    return find_spec(f"{__name__}.{group}") is not None


def groups() -> list[str]:
    """Names of the extension modules found on the search path."""
    names = set()
    for folder in __path__:
        try:
            entries = os.listdir(folder)
        except OSError:
            continue
        for entry in entries:
            for suffix in EXTENSION_SUFFIXES:
                if entry.endswith(suffix) and entry[: -len(suffix)].isidentifier():
                    names.add(entry[: -len(suffix)])
    return sorted(names)


def cuda_source_dirs() -> list[Path]:
    """Folders with the CUDA kernel sources compiled at run time, in search order.

    The build copies ``python/cpp/cuda`` to ``<module folder>/cuda``; without a build, the
    source tree's folder is used.
    """
    dirs = [Path(folder) / "cuda" for folder in __path__]
    dirs.append(Path(__file__).resolve().parents[3] / "cpp" / "cuda")
    return [d for d in dirs if d.is_dir()]
