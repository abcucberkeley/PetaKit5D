"""Core conventions every area codes against (see ``docs/conventions.md``).

Submodules, each importable on its own (``from petakit5d.core import bbox``):

- ``arrays``: ``(z, y, x)`` C-ordered volumes and the mapping to and from MATLAB's arrays,
  MATLAB's size rules on MATLAB-form values, class names
- ``dtypes``: MATLAB's class conversions, rounding, integer-class and single/double arithmetic
- ``bbox``: boxes ``(z0, y0, x0, z1, y1, x1)`` and ranges ``(start, stop)`` (0-based,
  end-exclusive) and their MATLAB forms
- ``axes``: per-axis vectors in ``(z, y, x)``, MATLAB dimensions and NumPy axes, axis strings
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
    "from_matlab": "arrays",
    "to_matlab": "arrays",
    "matlab_order": "arrays",
    "size_to_matlab": "arrays",
    "size_from_matlab": "arrays",
    "matlab_index_to_zyx": "arrays",
    "zyx_to_matlab_index": "arrays",
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
    "box_from_matlab": "bbox",
    "box_to_matlab": "bbox",
    "range_from_matlab": "bbox",
    "range_to_matlab": "bbox",
    "box_slices": "bbox",
    "box_shape": "bbox",
    "full_box": "bbox",
    "check_box": "bbox",
    "clamp_box": "bbox",
    "intersect_boxes": "bbox",
    "pad_widths": "bbox",
    "crop": "bbox",
    "from_yxz": "axes",
    "to_yxz": "axes",
    "from_xyz": "axes",
    "to_xyz": "axes",
    "dim_to_axis": "axes",
    "axis_to_dim": "axes",
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
    "axes",
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


# the CPU affinity is taken now, before a kernel's OpenMP region can bind the main thread to one
# CPU (OMP_PROC_BIND); core.threads is light (no NumPy)
from petakit5d.core import threads as _threads  # noqa: E402,F401


def __getattr__(name: str):
    import importlib

    if name in _SUBMODULES:
        return importlib.import_module(f"{__name__}.{name}")
    if name in _EXPORTS:
        return getattr(importlib.import_module(f"{__name__}.{_EXPORTS[name]}"), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return __all__
