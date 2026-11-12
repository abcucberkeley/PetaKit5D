"""Memory estimate of a processing step (``XR_estimateComputingMemory``), in GB.

Image sizes are ``(nz, ny, nx)`` in the Python API (``im_size`` and the returned size);
MATLAB's ``[ny nx nz]`` is converted by ``matlab_in`` and ``matlab_out`` for the MATLAB-name
layer.

Each step among ``steps`` is estimated (``any(contains(steps, name))``, MATLAB since 4b79f9b).
The GPU estimate is None when deconvolution is not among them (MATLAB leaves it unassigned).
MATLAB's quirk is kept: ``steps`` as text (``'deskew'``) gives as many estimates as it has
characters, the first one set.
"""

from __future__ import annotations

import os
import warnings
from collections.abc import Sequence
from typing import Any

from petakit5d.params import Param, register

__all__ = ["estimate_computing_memory", "matlab_in", "matlab_out"]

GB = 1024.0**3


class MatlabError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _zyx_to_matlab(sz: Sequence[Any]) -> list[Any]:
    # (z, y, x[, c]) -> [y x z [c]]
    sz = list(sz)
    if len(sz) < 3:
        return sz
    return [sz[1], sz[2], sz[0], *sz[3:]]


def _matlab_to_zyx(sz: Sequence[Any]) -> tuple[Any, ...]:
    # [y x z [c]] -> (z, y, x[, c])
    sz = list(sz)
    if len(sz) < 3:
        return tuple(sz)
    return (sz[2], sz[0], sz[1], *sz[3:])


def matlab_in(kwargs: dict[str, Any]) -> dict[str, Any]:
    """MATLAB-name layer: ``imSize`` ``[ny nx nz]`` -> ``im_size`` ``(nz, ny, nx)``."""
    out = dict(kwargs)
    if out.get("im_size") is not None and len(out["im_size"]) > 0:
        out["im_size"] = _matlab_to_zyx(out["im_size"])
    return out


def matlab_out(result: tuple[Any, ...]) -> tuple[Any, ...]:
    """MATLAB-name layer: the returned ``(nz, ny, nx)`` size as MATLAB's ``[ny nx nz]``."""
    mem, gpu, raw, sz = result
    if sz is not None and len(sz) > 0:
        sz = tuple(_zyx_to_matlab(sz))
    return mem, gpu, raw, sz


def _contains(steps: Any, name: str) -> Any:
    # contains(steps, name, 'IgnoreCase', true): a bool for text, a list for a cell
    low = name.lower()
    if isinstance(steps, str):
        return low in steps.lower()
    if isinstance(steps, Sequence):
        return [isinstance(s, str) and low in s.lower() for s in steps]
    raise MatlabError("First argument must be text.", "MATLAB:string:MustBeCharCellArrayOrString")


def _any(v: Any) -> bool:
    # any(contains(...)): a bool for text, any of a cell's
    if isinstance(v, list):
        return any(v)
    return bool(v)


def _strcmpi(steps: Any, name: str) -> list[int]:
    # positions where strcmpi(steps, name) is true
    low = name.lower()
    if isinstance(steps, str):
        return [0] if steps.lower() == low else []
    return [i for i, s in enumerate(steps) if isinstance(s, str) and s.lower() == low]


def _numel(steps: Any) -> int:
    return len(steps)


def _prod(sz: Sequence[Any]) -> float:
    p = 1.0
    for s in sz:
        p = p * float(s)
    return p


@register(
    "XR_estimateComputingMemory",
    Param("filePath", None, "Image file (tiff or zarr)."),
    Param("steps", None, "Steps to estimate: 'deskew', 'rotate' or 'deconvolution'."),
    Param(
        "imSize",
        "@(x) isnumeric(x) && (isempty(x) || numel(x) == 3)",
        "Image size, used when the file does not exist ((nz, ny, nx) in Python).",
    ),
    Param(
        "dataSize",
        "@(x) isnumeric(x) && (isempty(x) || numel(x) == 1)",
        "Number of voxels (uint16), used without deskew.",
    ),
    Param("memFactors", None, "Memory factors of deskew, rotate and deconvolution."),
    Param("cudaDecon", "@islogical", "Deconvolution on the GPU (gives the GPU estimate)."),
    Param("GPUMemFactor", None, "GPU memory factor of deconvolution."),
    Param("GPUMaxMem", "@isnumeric", "GPU memory in GB (not used)."),
    task="T05",
)
def estimate_computing_memory(
    file_path: str,
    steps: Any = ("deskew", "rotate", "deconvolution"),
    *,
    im_size: Sequence[float] | None = None,
    data_size: float | None = None,
    mem_factors: Sequence[float] = (15, 5, 15),
    cuda_decon: bool = True,
    gpu_mem_factor: float = 1.5,
    gpu_max_mem: float = 12,
) -> tuple[Any, ...]:
    """Estimate the memory a step needs, in GB.

    Returns ``(est_required_memory, est_required_gpu_mem, raw_image_size, im_size)``: one
    estimate per step (a float64 array), the GPU estimate of deconvolution (NaN without
    ``cuda_decon``; None when no deconvolution estimate is made), the size of the image as
    single in GB, and the image size ``(nz, ny, nx)`` (from the file when it exists).
    """
    import numpy as np

    from petakit5d.core.dtypes import NAN
    from petakit5d.core.paths import fileparts

    del gpu_max_mem
    file_path = os.fspath(file_path) if file_path is not None else ""
    # MATLAB order inside: [ny nx nz]
    msize = None if im_size is None else _zyx_to_matlab(im_size)
    empty_data = data_size is None or (np.size(data_size) == 0)

    with warnings.catch_warnings():
        # MATLAB hides its warnings during the estimate
        warnings.simplefilter("ignore")
        if empty_data or _any(_contains(steps, "deskew")):
            if not (file_path and os.path.exists(file_path)):
                if msize is None or len(msize) == 0:
                    raise MatlabError(f"File {file_path} does not exist!")
            else:
                from petakit5d.io.imagesize import get_image_size

                msize = list(get_image_size(file_path))
            if len(msize) < 3:
                raise MatlabError(
                    "Index exceeds the number of array elements. Index must not exceed "
                    f"{len(msize)}.",
                    "MATLAB:badsubscript",
                )
            # in some cases only the size of the first slice is known
            ext = fileparts(file_path)[2]
            if float(msize[2]) == 1 and ext in (".tif", ".tiff"):
                if not os.path.isfile(file_path):
                    # dir() of a missing file is empty: MATLAB's dir_info.bytes fails
                    raise MatlabError("Not enough input arguments.", "MATLAB:minrhs")
                raw = os.stat(file_path).st_size / GB
            else:
                # in single size
                raw = _prod(msize) * 4 / GB
        else:
            raw = float(np.asarray(data_size).ravel()[0]) * 2 / GB

        est = np.zeros(_numel(steps), dtype=np.float64)
        gpu = None
        mf = [float(x) for x in np.asarray(mem_factors, dtype=np.float64).ravel(order="F")]
        if _any(_contains(steps, "deskew")):
            # 300 slices as threshold
            q = float(msize[2]) / 300
            est[_strcmpi(steps, "deskew")] = raw * mf[0] * max(1.0, q * q)
        if _any(_contains(steps, "rotate")):
            est[_strcmpi(steps, "rotate")] = raw * mf[1]
        if _any(_contains(steps, "deconvolution")):
            est[_strcmpi(steps, "deconvolution")] = raw * mf[2]
            gpu = float(NAN)
            if cuda_decon:
                gpu = raw * float(gpu_mem_factor)

    out_size = None if msize is None else _matlab_to_zyx(msize)
    return est, gpu, raw, out_size
