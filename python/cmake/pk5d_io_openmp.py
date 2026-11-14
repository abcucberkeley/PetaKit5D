"""The OpenMP runtime the cpp-zarr and cpp-tiff packages load (cmake/Pk5dOpenMP.cmake, Linux).

Prints one line per libgomp file that loading a package's compiled modules maps into this
process, as ``<package>\t<path>``. The modules are loaded with ctypes, without importing the
packages (their import needs NumPy, which a build environment may not have), so the dynamic
loader resolves their libraries as it does in a port process. Prints nothing where neither
package is installed or /proc/self/maps doesn't exist.
"""

import ctypes
import importlib.machinery
import importlib.util
import os


def mapped() -> set[str]:
    found = set()
    try:
        with open("/proc/self/maps") as f:
            for line in f:
                parts = line.split(None, 5)
                if len(parts) == 6 and os.path.basename(parts[5].strip()).startswith("libgomp"):
                    found.add(parts[5].strip())
    except OSError:
        pass
    return found


def main() -> None:
    seen = mapped()
    for pkg in ("cppzarr", "cpptiff"):
        try:
            spec = importlib.util.find_spec(pkg)
        except (ImportError, ValueError):
            continue
        if spec is None or not spec.submodule_search_locations:
            continue
        for folder in spec.submodule_search_locations:
            try:
                names = sorted(os.listdir(folder))
            except OSError:
                continue
            for name in names:
                if not name.endswith(tuple(importlib.machinery.EXTENSION_SUFFIXES)):
                    continue
                try:
                    ctypes.CDLL(os.path.join(folder, name))
                except OSError:
                    continue
                now = mapped()
                for path in sorted(now - seen):
                    print(f"{pkg}\t{os.path.realpath(path)}")
                seen = now


if __name__ == "__main__":
    main()
