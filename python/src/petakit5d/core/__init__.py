"""Core conventions every area codes against (see ``docs/conventions.md``).

Submodules, each importable on its own (``from petakit5d.core import bbox``):

- ``arrays``: (y, x, z) Fortran-ordered volumes, MATLAB size rules, class names, permutations
- ``dtypes``: MATLAB's class conversions, rounding, integer-class and single/double arithmetic
- ``bbox``: MATLAB boxes ``[ymin, xmin, zmin, ymax, xmax, zmax]`` (1-based, inclusive) and slices
- ``blocks``: batch and block coordinates of the large-file paths
- ``text``: ``sprintf``, ``num2str``, ``mat2str``, ``str2double`` as MATLAB formats and parses
- ``jsonenc``: ``jsonencode`` as MATLAB writes it
- ``textfiles``: ``readTextFile``, ``writeTextFile``, ``writeJsonFile``
- ``paths``: ``simplifyPath``, ``fileparts``, ``dir`` order, ``dir_recursive``, ``mkdir_recursive``
- ``ids``: ``get_uuid``, ``get_hostname``
- ``threads``: CPUs this process may use, per-call thread counts
- ``memory``: ``getSystemMemory`` (cgroup v1 and v2, Slurm, macOS, Windows)

Importing ``petakit5d.core`` loads none of them; names below load their module on first use.
"""

_EXPORTS = {
    # arrays
    "as_volume": "arrays",
    "check_volume": "arrays",
    "fzeros": "arrays",
    "fempty": "arrays",
    "msize": "arrays",
    "mndims": "arrays",
    "mnumel": "arrays",
    "msqueeze": "arrays",
    "mpermute": "arrays",
    "sub2ind": "arrays",
    "ind2sub": "arrays",
    "matlab_class": "arrays",
    "numpy_dtype": "arrays",
    "axis_order_mapping": "arrays",
    "data_type_to_byte_number": "arrays",
    # dtypes
    "cast": "dtypes",
    "to_class": "dtypes",
    "mround": "dtypes",
    "plus": "dtypes",
    "minus": "dtypes",
    "times": "dtypes",
    "rdivide": "dtypes",
    "idivide": "dtypes",
    "EPS": "dtypes",
    "NAN": "dtypes",
    # boxes and blocks
    "bbox_to_slices": "bbox",
    "slices_to_bbox": "bbox",
    "bbox_size": "bbox",
    "full_bbox": "bbox",
    "check_bbox": "bbox",
    "clamp_bbox": "bbox",
    "crop": "bbox",
    "zarr_chunk_coordinates": "blocks",
    "sub_volume_coordinates": "blocks",
    "bbox_to_blocks": "blocks",
    "find_good_factor_number": "blocks",
    # text and files
    "sprintf": "text",
    "num2str": "text",
    "int2str": "text",
    "mat2str": "text",
    "mat2str_comma": "text",
    "str2double": "text",
    "through_text": "text",
    "jsonencode": "jsonenc",
    "read_text_file": "textfiles",
    "write_text_file": "textfiles",
    "write_json_file": "textfiles",
    # system
    "simplify_path": "paths",
    "fileparts": "paths",
    "matlab_dir": "paths",
    "dir_recursive": "paths",
    "mkdir_recursive": "paths",
    "get_uuid": "ids",
    "get_hostname": "ids",
    "available_cpus": "threads",
    "resolve_threads": "threads",
    "get_system_memory": "memory",
}

_SUBMODULES = (
    "arrays",
    "dtypes",
    "bbox",
    "blocks",
    "text",
    "jsonenc",
    "textfiles",
    "paths",
    "ids",
    "threads",
    "memory",
)

__all__ = sorted(_EXPORTS) + list(_SUBMODULES)


def __getattr__(name: str):
    import importlib

    if name in _SUBMODULES:
        return importlib.import_module(f"{__name__}.{name}")
    if name in _EXPORTS:
        return getattr(importlib.import_module(f"{__name__}.{_EXPORTS[name]}"), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return __all__
