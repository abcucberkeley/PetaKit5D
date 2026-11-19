"""Flat-field correction of 2-D phase tiles: ``XR_phase_image_flat_field_correction``, with
BaSiCPy for the flat field.

- ``phase_prefilter(image, const)``: ``single(im) - imgaussfilt(single(im), 100) + const``, the
  tile with its large-scale background removed (``matlabcompat.filters.imgaussfilt``);
- ``basic_flatfield(stack)``: the flat field of a ``(t, y, x)`` stack of tiles, estimated with
  BaSiCPy (``basicpy.BaSiC``) with the options MATLAB passes to BaSiC (no dark field) and
  BaSiCPy's defaults otherwise (``fitting_mode='approximate'``, the original method's
  automatic smoothness, working size 128), on the CPU;
- ``phase_image_flatfield_correction(fns, out_path, const)``: the whole function, the tiles of
  a time point pre-filtered, their flat field estimated, each tile divided by it and saved as
  uint16 ``<out_path><name>.tif``, and ``minval`` and ``flatfield`` saved in
  ``<out_path>BaSiC_flatfield/<first name>_flatfield.mat``.

MATLAB uses BaSiC (``third_parties/BaSiC-master``, CC BY-NC-ND: not ported, not read); BaSiCPy
(MIT) implements the same published method (Peng et al., Nat. Commun. 8, 14836, 2017). The
flat fields of the two agree within the tolerance in ``docs/parity/flatfield.md`` (approved, N8),
not bit for bit. BaSiCPy is an optional dependency (the ``basic`` extra and
``pip install --no-deps basicpy``: BaSiCPy 2.0.0 declares ``scipy<1.13``); without it the
estimation raises ``ImportError``. Its fit runs with one PyTorch thread, so the flat field
does not depend on the machine's thread count.

``out_path`` is joined to the names as text, as MATLAB joins them (give it a trailing
separator). Results go through a temporary name next to the output (MATLAB uses ``/dev/shm``
on Linux); the files written are the same. Arrays are the port's: tiles ``(y, x)``, stacks
``(t, y, x)``.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "BASICPY_VERSION",
    "MatlabError",
    "basic_available",
    "basic_flatfield",
    "phase_image_flatfield_correction",
    "phase_prefilter",
]

BASICPY_VERSION = "2.0.0"  # the version the tolerance was measured with
GAUSS_SIGMA = 100
_F32 = np.dtype(np.float32)


class MatlabError(RuntimeError):
    """An error MATLAB raises, with its identifier."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def basic_available() -> bool:
    """Whether BaSiCPy (and PyTorch) can be imported."""
    import importlib.util

    return all(importlib.util.find_spec(m) is not None for m in ("basicpy", "torch"))


def _basicpy():
    try:
        import basicpy
    except ImportError as err:
        raise ImportError(
            "The phase-image flat field needs BaSiCPy: pip install 'petakit5d[basic]' and "
            f"pip install --no-deps basicpy=={BASICPY_VERSION}"
        ) from err
    return basicpy


def phase_prefilter(
    image: np.ndarray, const: Any = 1100, *, num_threads: int | None = None
) -> np.ndarray:
    """``imf = single(im); imf - imgaussfilt(imf, 100) + const`` of a ``(y, x)`` tile
    (``num_threads``: the filter's threads, all available by default)."""
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat.filters import imgaussfilt

    imf = np.asarray(image)
    if imf.dtype != _F32:
        imf = np.asarray(dtypes.to_class(imf, "single"))
    g = imgaussfilt(imf, GAUSS_SIGMA, num_threads=num_threads)
    d = np.subtract(imf, g, dtype=_F32)
    return np.asarray(dtypes.plus(d, np.float64(const)), dtype=_F32)


def basic_flatfield(stack: np.ndarray) -> np.ndarray:
    """The flat field ``(y, x)`` of a ``(t, y, x)`` stack of tiles, in double as MATLAB's BaSiC
    gives it: BaSiCPy's ``BaSiC(get_darkfield=False, fitting_mode='approximate')`` fit on the
    CPU with one PyTorch thread (BaSiC's ``'darkfield', 'false'``; the other options BaSiCPy's
    defaults, which tolerances of 1e-6 or reweighting at 1e-3 don't bring closer to MATLAB's).
    BaSiCPy's result (float32) depends on PyTorch's thread count in the last bits, hence the
    one thread."""
    basicpy = _basicpy()
    import torch

    im = np.ascontiguousarray(np.asarray(stack, dtype=_F32))
    threads = torch.get_num_threads()
    torch.set_num_threads(1)
    try:
        b = basicpy.BaSiC(get_darkfield=False, fitting_mode="approximate", device="cpu")
        b.fit(im)
        ff = np.asarray(b.flatfield)
    finally:
        torch.set_num_threads(threads)
    return np.ascontiguousarray(ff.reshape(im.shape[1:]), dtype=np.float64)


def _names(fns: Any) -> list[str]:
    if isinstance(fns, str | os.PathLike):
        # MATLAB indexes the characters of a char array with braces
        raise MatlabError(
            "Brace indexing is not supported for variables of this type.",
            "MATLAB:cellRefFromNonCell",
        )
    return [os.fspath(f) for f in fns]


@register(
    "XR_phase_image_flat_field_correction",
    Param("fns", None, "Tiff files of the tiles of a time point."),
    Param("outPath", None, "Output folder, with a trailing separator (joined as text)."),
    Param("const", None, "Offset added after the background subtraction.", convert="same"),
    task="B07",
)
def phase_image_flatfield_correction(fns: Sequence[str], out_path: str, const: Any = 1100) -> None:
    """Flat-field correction of the phase tiles of a time point
    (``XR_phase_image_flat_field_correction``): each tile minus its Gaussian blur (sigma 100)
    plus ``const``, the flat field of all of them (``basic_flatfield``), each divided by it and
    saved as uint16 ``<out_path><name>.tif``; ``minval`` (the smallest pre-filtered value) and
    ``flatfield`` go to ``<out_path>BaSiC_flatfield/<first name>_flatfield.mat``. Nothing is
    done when every output exists."""
    from concurrent.futures import ThreadPoolExecutor

    from petakit5d.core import dtypes
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.core.threads import resolve_threads
    from petakit5d.io import read_tiff, write_mat, write_tiff
    from petakit5d.jobs.files import exist
    from petakit5d.matlabcompat import reduce

    files = _names(fns)
    out_path = os.fspath(out_path)
    fsns = [fileparts(f)[1] for f in files]
    fnouts = [f"{out_path}{x}.tif" for x in fsns]
    uuid = get_uuid()
    tmpouts = [f"{out_path}{x}_{uuid}.tif" for x in fsns]
    if all(exist(f) for f in fnouts):
        print("All output files exist, skip it!")
        return

    # the tiles are read and pre-filtered, and later divided and converted, in threads (a tile
    # per thread, each filter on one thread): the values of one tile after the other
    nt = resolve_threads(None)

    def prefiltered(fn: str) -> np.ndarray:
        return phase_prefilter(read_tiff(fn), const, num_threads=1)

    with ThreadPoolExecutor(max_workers=nt) as ex:
        ims = list(ex.map(prefiltered, files))
    im = np.stack(ims)
    minval, _ = reduce.min(im)
    flatfield = basic_flatfield(im)
    del im

    def corrected(f: int) -> np.ndarray:
        imf = np.asarray(dtypes.rdivide(ims[f], flatfield), dtype=_F32)
        return np.asarray(dtypes.cast(imf, "uint16"))

    with ThreadPoolExecutor(max_workers=nt) as ex:
        outs = list(ex.map(corrected, range(len(files))))
    for f in range(len(files)):
        write_tiff(outs[f], tmpouts[f])
        os.replace(tmpouts[f], fnouts[f])

    ff_path = f"{out_path}BaSiC_flatfield/"
    if not os.path.isdir(ff_path):
        os.makedirs(ff_path)
    write_mat(
        f"{ff_path}/{fsns[0]}_flatfield.mat", {"minval": np.float32(minval), "flatfield": flatfield}
    )
