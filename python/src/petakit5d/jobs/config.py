"""Cluster configuration files (JSON or .mat), read and merged as MATLAB's job framework does.

A config file's keys override the caller's parameters of the same name, matched without regard
to case (``MemPerCPU`` sets ``memPerCPU``); a key that names no parameter is kept but has no
effect. Keys MATLAB uses only for its own jobs (``MATLAB_ONLY_KEYS``) are accepted and ignored.
``PythonLaunchStr`` replaces ``MatlabLaunchStr``: what a job runs before Python (module loads,
environment activation, exports) followed by the interpreter; by default the interpreter that
submits.

JSON is decoded as MATLAB's ``jsondecode`` names the fields: keys that are not MATLAB names are
made valid (``matlab.lang.makeValidName``) and repeated keys made unique (``a``, ``a_1``), so the
first of two equal keys is the one that matches. Values are Python's (``null`` is ``None``,
arrays are lists); ``is_empty`` is MATLAB's ``isempty`` for them.
"""

from __future__ import annotations

import json
import os
import shlex
import sys
import threading
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    "MATLAB_ONLY_KEYS",
    "PYTHON_LAUNCH_KEY",
    "ClusterConfig",
    "ConfigCache",
    "cluster_config_value",
    "clear_cache",
    "decode_json",
    "default_python_launch_str",
    "is_empty",
    "load_config",
    "matlab_valid_name",
]

# keys of MATLAB's compiled and MATLAB-job modes: accepted in a config file and ignored
MATLAB_ONLY_KEYS = ("MatlabLaunchStr", "MCRParam", "MCCMasterStr", "MCRCacheRoot", "mccMode")
PYTHON_LAUNCH_KEY = "PythonLaunchStr"

# MATLAB's namelengthmax
_NAME_LENGTH_MAX = 63


def is_empty(value: Any) -> bool:
    """MATLAB's isempty for decoded config values: None, '', [] and {} are empty"""
    if value is None:
        return True
    if isinstance(value, (str, bytes, list, tuple, dict)):
        return len(value) == 0
    size = getattr(value, "size", None)
    return isinstance(size, int) and size == 0


def matlab_valid_name(key: str) -> str:
    """the field name jsondecode gives an object key (MATLAB R2024b, checked by the oracle)

    Whitespace is removed, and a lowercase letter after whitespace inside the key becomes
    uppercase ('mem per cpu' -> 'memPerCpu'); other characters than ASCII letters, digits and
    '_' become '_'; a name that doesn't start with a letter gets an 'x' ('1abc' -> 'x1abc',
    '' -> 'x'); at most 63 characters. Keywords stay ('end').
    """
    out: list[str] = []
    upper_next = False
    for ch in key:
        if ch.isspace():
            upper_next = bool(out)
            continue
        if upper_next and "a" <= ch <= "z":
            ch = ch.upper()
        upper_next = False
        out.append(ch)
    name = "".join(c if (c.isascii() and (c.isalnum() or c == "_")) else "_" for c in out)
    if not name or not ("a" <= name[0].lower() <= "z"):
        name = "x" + name
    return name[:_NAME_LENGTH_MAX]


def _unique_names(names: list[str]) -> list[str]:
    """matlab.lang.makeUniqueStrings: later repeats get _1, _2, ... (avoiding taken names)"""
    taken = set(names)
    seen: set[str] = set()
    out = []
    for name in names:
        if name not in seen:
            seen.add(name)
            out.append(name)
            continue
        k = 1
        while f"{name}_{k}"[:_NAME_LENGTH_MAX] in taken:
            k += 1
        new = f"{name}_{k}"[:_NAME_LENGTH_MAX]
        taken.add(new)
        seen.add(new)
        out.append(new)
    return out


def _pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    names = _unique_names([matlab_valid_name(k) for k, _ in pairs])
    return dict(zip(names, (v for _, v in pairs), strict=True))


def decode_json(text: str) -> dict[str, Any]:
    """the top-level object of a config file, with MATLAB's field names (nested objects too)"""
    data = json.loads(text, object_pairs_hook=_pairs_hook)
    if not isinstance(data, dict):
        raise ValueError("a cluster configuration file must hold one JSON object")
    return data


def default_python_launch_str() -> str:
    """the default PythonLaunchStr: the interpreter of this process"""
    return shlex.quote(sys.executable)


@dataclass(frozen=True)
class ClusterConfig:
    """the keys and values of one cluster configuration file, in file order

    ``get``/``in`` match names without regard to case; with several keys of one name (in any
    case) the last one counts, as in MATLAB.
    """

    path: str | None = None
    entries: tuple[tuple[str, Any], ...] = ()
    _index: dict[str, int] = field(default_factory=dict, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        for i, (key, _) in enumerate(self.entries):
            self._index[key.lower()] = i

    # construction

    @classmethod
    def empty(cls) -> ClusterConfig:
        return cls()

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any], path: str | None = None) -> ClusterConfig:
        return cls(path, tuple(values.items()))

    @classmethod
    def load(cls, path: str | os.PathLike[str]) -> ClusterConfig:
        """read a config file (.json, or .mat through petakit5d.io.matfile); no cache"""
        path = os.fspath(path)
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"The cluster configuration file {path} does not exist, please check and "
                "provide the right one!"
            )
        ext = os.path.splitext(path)[1]
        if ext == ".json":
            # (MATLAB reads the bytes as characters; files are ASCII in practice)
            with open(path, encoding="utf-8") as f:
                values = decode_json(f.read())
        elif ext == ".mat":
            values = _read_mat_config(path)
        else:
            raise ValueError(f"The cluster configuration file {path} is not a .json or .mat file")
        return cls.from_mapping(values, path)

    # lookups

    def __bool__(self) -> bool:
        return bool(self.entries)

    def __len__(self) -> int:
        return len(self.entries)

    def __contains__(self, name: object) -> bool:
        return isinstance(name, str) and name.lower() in self._index

    def keys(self) -> list[str]:
        return [k for k, _ in self.entries]

    def get(self, name: str, default: Any = None) -> Any:
        """the value of key name (any case; the last such key), else default"""
        i = self._index.get(name.lower())
        return default if i is None else self.entries[i][1]

    def ignored_keys(self) -> list[str]:
        """the MATLAB-only keys of the file (accepted, without effect)"""
        ignored = {k.lower() for k in MATLAB_ONLY_KEYS}
        return [k for k, _ in self.entries if k.lower() in ignored]

    def overrides(self, names: Iterable[str]) -> dict[str, Any]:
        """{name: value} for each of names (MATLAB parameter names) that a key sets

        Keys are applied in file order, so of two keys for one name the later wins. MATLAB-only
        keys are left out, also when names has them.
        """
        by_lower: dict[str, str] = {}
        for n in names:
            by_lower.setdefault(n.lower(), n)
        ignored = {k.lower() for k in MATLAB_ONLY_KEYS}
        out: dict[str, Any] = {}
        for key, value in self.entries:
            low = key.lower()
            if low in ignored:
                continue
            name = by_lower.get(low)
            if name is not None:
                out[name] = value
        return out

    def unknown_keys(self, names: Iterable[str]) -> list[str]:
        """keys that name none of names and are not MATLAB-only (they have no effect)"""
        known = {n.lower() for n in names} | {k.lower() for k in MATLAB_ONLY_KEYS}
        return [k for k, _ in self.entries if k.lower() not in known]

    def merge(
        self, params: Mapping[str, Any], names: Mapping[str, str] | None = None
    ) -> dict[str, Any]:
        """params with the config's values in place of the caller's

        names maps MATLAB parameter names to the keys of params (e.g. ``{'memPerCPU':
        'mem_per_cpu'}``); without it the keys of params are the MATLAB names. Keys of the file
        that name no parameter are not added (see unknown_keys).
        """
        if names is None:
            names = {k: k for k in params}
        merged = dict(params)
        for name, value in self.overrides(names).items():
            merged[names[name]] = value
        return merged

    @property
    def python_launch_str(self) -> str:
        """PythonLaunchStr of the file, else the interpreter of this process"""
        value = self.get(PYTHON_LAUNCH_KEY)
        return default_python_launch_str() if is_empty(value) else str(value)


def _read_mat_config(path: str) -> dict[str, Any]:
    try:
        from petakit5d.io.matfile import read_mat
    except ImportError as err:  # (F07's module; until it exists only JSON configs are read)
        raise ImportError(f"reading {path} needs petakit5d.io.matfile ({err})") from err
    values = read_mat(path)
    return dict(values)


class ConfigCache:
    """config files read once, again only when a file changes (MATLAB's persistent cache)"""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._files: dict[str, tuple[tuple[int, int], ClusterConfig]] = {}

    def load(self, path: str | os.PathLike[str]) -> tuple[ClusterConfig, bool]:
        """the config of a file, and whether it was read now (a new or changed file)"""
        path = os.fspath(path)
        try:
            st = os.stat(path)
        except OSError:
            st = None
        stamp = (st.st_mtime_ns, st.st_size) if st is not None else (-1, -1)
        key = os.path.abspath(path)
        with self._lock:
            hit = self._files.get(key)
            if hit is not None and st is not None and hit[0] == stamp:
                return hit[1], False
        config = ClusterConfig.load(path)
        with self._lock:
            self._files[key] = (stamp, config)
        return config, True

    def clear(self) -> None:
        with self._lock:
            self._files.clear()


_CACHE = ConfigCache()


def load_config(
    path: str | os.PathLike[str] | None, cache: ConfigCache | None = _CACHE
) -> ClusterConfig:
    """the config of path ('' or None: an empty config), cached unless cache is None

    A path that is given but missing raises FileNotFoundError, as in MATLAB's framework.
    """
    if path is None or os.fspath(path) == "":
        return ClusterConfig.empty()
    if cache is None:
        return ClusterConfig.load(path)
    return cache.load(path)[0]


def clear_cache() -> None:
    """forget the cached config files (MATLAB: clear generic_computing_frameworks_wrapper)"""
    _CACHE.clear()


def cluster_config_value(
    config_file: str | os.PathLike[str] | None, name: str, value: Any = None
) -> Any:
    """cluster_config_value.m: the value of name in the file, else value

    No file, a missing file, another extension than .json/.mat or no such key give value.
    """
    if config_file is None or os.fspath(config_file) == "" or not os.path.exists(config_file):
        return value
    if os.path.splitext(os.fspath(config_file))[1] not in (".json", ".mat"):
        return value
    return ClusterConfig.load(config_file).get(name, value)
