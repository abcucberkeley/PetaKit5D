"""Deskew and rotation of files (``XR_deskewRotateFrame``) and the intensity rescale of its
``rescaleRotate`` option (``scale_contrast``, after ``scaleContrast``).

``deskew_rotate_file`` reads one file, or several partial files of one volume (concatenated
along z in the order given), and writes MATLAB's outputs next to the first one:

- the combined deskew and rotation (``dsr_combined``, the default): ``<dir>/DSR/<name>.tif``
  (or ``.zarr``), through ``deskew_rotate_volume``;
- the separate path: the deskewed volume ``<dir>/DS/<name>.tif`` (``deskew_volume``), then,
  with ``rotate``, its rotation ``<dir>/DSR/<name>.tif`` (``rotate_volume``);
- the z MIP of each, ``<DS or DSR>/MIPs/<name>_MIP_z.tif`` (or ``.zarr``), with ``save_mip``.

Before the deskew, in MATLAB's order: the input box (``input_bbox``, MATLAB's ``crop3d``), the
flip of the z stack, the flat-field correction (``XR_LSFlatFieldCorrection`` with the
flat-field and background images), the camera background removal (MATLAB calls a function the
repository does not have: the port raises MATLAB's error, known issue) and the input axis order.
After it: the output axis order, the MIP, the rescale and the 16-bit or single conversion.

Files and names are MATLAB's: the folders ``<dir>/<DSDirName><surffix>/`` (MATLAB's text, a
doubled separator included), the temporary results ``<name>_<uuid>.tif`` (``.zarr``) in the
result's folder and the DSR MIP's ``<name>_MIP_z.tiff_<uuid>`` (``.zarr_<uuid>``); the DS MIP is
written under its final name, as in MATLAB. Each task writes its MIP before its result and
publishes the result with one rename (``jobs.publish``): a result another writer published
first is kept, also with ``overwrite``, which removes an existing DS before it is made again.
Without ``save_3d_stack`` the result is an empty file. An existing DSR is never made again (as
MATLAB, ``overwrite`` is the DS's only). Zarr results without a given uuid are private to the
call and their chunks are written in place (cpp-zarr's ``uuid=False``), as MATLAB writes them.

The combined path reads tiff data into MATLAB's layout (one threaded in-place transpose of
each page) and runs the deskew kernels on that layout at the mex's speed; the result is turned
back the same way for tiff, and written from a ``(z, y, x)`` view of that layout for zarr.
The ``imwarp`` paths (objective scan, resampling, cubic) and the separate path work on
``(z, y, x)`` arrays. Every output voxel is MATLAB's: bit for bit.

Differences from MATLAB, each reported:

- above 1000 slices the separate path deskews in blocks of rows as MATLAB's ``blockedImage``
  does (MATLAB's block size and 5-row borders: the same voxels); its DS file is written as
  every port tiff (cpp-tiff), where MATLAB's ``MPageTiffAdapter`` writes another tiff layout;
- the combined path's resampling of inputs above 1000 slices before the deskew uses MATLAB's
  ``imresize``, which the port does not have yet (``NotImplementedError``);
- the flat field and background are prepared once per process for the same files (path, size
  and modification time), where MATLAB reads and prepares them for every frame: the same
  result (``_flat_field``).

Arrays are ``(z, y, x)``; ``input_bbox`` is ``(z0, y0, x0, z1, y1, x1)`` of the file's array
(before the input axis order), ``block_size`` ``(z, y, x)``, ``resample_factor`` ``(z, y, x)``
(or one value, or ``(z, xy)``), axis orders our strings (``'zyx'``). ``num_threads`` (None: the
CPUs this process may use) never changes a result.
"""

from __future__ import annotations

import math
import os
import time
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = ["deskew_rotate_file", "scale_contrast"]

_INTERP = "@(x) any(strcmpi(x, {'cubic', 'linear'}))"
# the flat field and background prepared for the last files seen (see _flat_field)
_FF_CACHE: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}
# arrays from this size are copied in threads
_THREADED_BYTES = 1 << 22
# elements per chunk of the element-wise conversions (their float64 temporaries: 8 MB a chunk)
_CHUNK = 1 << 20


class MatlabError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id."""

    def __init__(self, message: str, identifier: str):
        super().__init__(message)
        self.identifier = identifier


def _logical(x: Any) -> bool:
    return bool(np.asarray(x).reshape(-1)[0]) if np.size(x) else False


# -- scaleContrast -----------------------------------------------------------------------------


def scale_contrast(
    a: np.ndarray,
    range_in: Sequence[Any] | None = None,
    range_out: Sequence[Any] = (0, 255),
    *,
    num_threads: int | None = None,
) -> np.ndarray:
    """``scaleContrast(A, rangeIn, rangeOut)``: ``A`` mapped linearly from ``range_in`` to
    ``range_out``, ``(double(A) - r1) ./ (r2 - r1) .* (o2 - o1) + o1``, in MATLAB's classes.

    A new implementation from MATLAB's behaviour (``scaleContrast`` is another author's GPL
    function, not read; black-box probes, kept as oracle cases). Measured:

    - the data enter as double, so the result's class is the range's combined with double:
      single for a single ``range_in`` (or ``range_out``), double otherwise; an integer
      ``range_in`` with more than one value of data raises ``MATLAB:mixedClasses``;
    - each operation rounds in that class, in that order; values outside ``range_in`` are not
      clipped; NaN and Inf follow IEEE arithmetic (an Inf end gives 0 and NaN);
    - ``r1 == r2`` (also ``-0`` and ``0``) gives double zeros of ``A``'s shape;
    - ``range_in`` defaults to ``[min(A(:)) max(A(:))]`` (NaN ignored), ``range_out`` to
      ``[0 255]``.

    Large arrays are computed in chunks of their elements into one result, in
    ``num_threads`` threads (None: the CPUs available; the same result).
    """
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat import reduce

    a = np.asarray(a)
    if range_in is None:
        lo, _ = reduce.min(a)
        hi, _ = reduce.max(a)
        range_in = (np.asarray(lo, dtype=a.dtype), np.asarray(hi, dtype=a.dtype))
    ri = _range_values(range_in)
    ro = _range_values(range_out)
    r1, r2 = ri[0], ri[1]
    if bool(np.asarray(dtypes.minus(r2, r1) == 0)):
        return np.zeros(a.shape, dtype=np.float64)
    d_in = dtypes.minus(r2, r1)
    d_out = dtypes.minus(ro[1], ro[0])
    o1 = ro[0]

    def scale(x: np.ndarray) -> np.ndarray:
        x = np.asarray(dtypes.double(x))
        try:
            t = dtypes.minus(x, r1)
        except TypeError as err:
            raise MatlabError(
                "Integers can only be combined with integers of the same class, or scalar doubles.",
                "MATLAB:mixedClasses",
            ) from err
        return dtypes.plus(dtypes.times(dtypes.rdivide(t, d_in), d_out), o1)

    if a.nbytes < _THREADED_BYTES:
        return np.asarray(scale(a))
    cls = np.asarray(scale(a.reshape(-1)[:2])).dtype
    return _elementwise(a, scale, cls, _threads(num_threads))


def _range_values(r: Any) -> list[Any]:
    """the first two values of a range in their MATLAB class (Python numbers are doubles)"""
    if isinstance(r, list | tuple) and all(
        isinstance(x, int | float) and not isinstance(x, bool) for x in r
    ):
        v = np.asarray(r, dtype=np.float64)
    else:
        v = np.asarray(r)
    v = v.reshape(-1)
    if v.size < 2:
        raise MatlabError(
            f"Index exceeds the number of array elements. Index must not exceed {v.size}.",
            "MATLAB:badsubscript",
        )
    return [v[0], v[1]]


# -- arrays ------------------------------------------------------------------------------------


def _threads(num_threads: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(num_threads)


def _copy(view: np.ndarray, n: int, out: np.ndarray | None = None) -> np.ndarray:
    """a C-contiguous copy of a strided view, its first axis split over ``n`` threads (NumPy's
    copies release the GIL)"""
    if out is None:
        out = np.empty(view.shape, dtype=view.dtype)
    if n < 2 or view.nbytes < _THREADED_BYTES or view.ndim == 0 or view.shape[0] < 2:
        out[...] = view
        return out
    from concurrent.futures import ThreadPoolExecutor

    edges = np.linspace(0, view.shape[0], min(view.shape[0], 4 * n) + 1).astype(int)

    def part(k: int) -> None:
        out[edges[k] : edges[k + 1]] = view[edges[k] : edges[k + 1]]

    with ThreadPoolExecutor(max_workers=n) as ex:
        list(ex.map(part, range(len(edges) - 1)))
    return out


def _to_class(vol: np.ndarray, cls: str, n: int) -> np.ndarray:
    """MATLAB's class function ``cls(vol)`` (``core.dtypes.to_class``: for single arrays to
    uint16, MATLAB's ``floor(|x| + 0.5)``), chunk by chunk in ``n`` threads into one result: its
    temporaries (float64 copies) stay the size of a chunk. Every chunk has two or more elements, so
    each takes the array rule of the whole."""
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import numpy_dtype

    if vol.size < _CHUNK:
        return np.ascontiguousarray(dtypes.to_class(vol, cls))
    return _elementwise(vol, lambda x: dtypes.to_class(x, cls), numpy_dtype(cls), n)


def _elementwise(a: np.ndarray, fn, dtype, n: int) -> np.ndarray:
    """``fn`` (an element-wise MATLAB operation) of ``a``, chunk by chunk of the elements in
    memory order (``_CHUNK`` each, every chunk two or more elements) in ``n`` threads, into one
    result of ``dtype``: the temporaries of ``fn`` stay the size of a chunk per thread."""
    a = np.ascontiguousarray(a)
    flat = a.reshape(-1)
    out = np.empty(a.shape, dtype=dtype)
    oflat = out.reshape(-1)
    edges = list(range(0, flat.size, _CHUNK)) + [flat.size]
    if len(edges) > 2 and edges[-1] - edges[-2] < 2:
        del edges[-2]

    def part(k: int) -> None:
        oflat[edges[k] : edges[k + 1]] = fn(flat[edges[k] : edges[k + 1]])

    if n < 2:
        for k in range(len(edges) - 1):
            part(k)
        return out
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=n) as ex:
        list(ex.map(part, range(len(edges) - 1)))
    return out


def _turn_pages(vol: np.ndarray, n: int) -> np.ndarray:
    """``(nz, a, b)`` C data as ``(nz, b, a)`` C data in the same memory: each page transposed
    in place (``petakit5d._ext.io``), else a threaded copy. ``(z, y, x)`` gives MATLAB's layout
    ``(z, x, y)`` (y fastest), and back."""
    from petakit5d import _ext

    nz, na, nb = vol.shape
    if _ext.available("io") and vol.flags.c_contiguous and vol.flags.writeable:
        flat = vol.reshape(-1)
        _ext.load("io").transpose_pages_inplace(flat, nb, na, nz, n)
        return flat.reshape(nz, nb, na)
    return _copy(np.swapaxes(vol, 1, 2), n)


# -- reading -----------------------------------------------------------------------------------


def _read_frame(paths: list[str], n: int, read_buffer: Any) -> np.ndarray:
    """``readtiff`` of each file and ``cat(3, ...)``, or ``readtiff`` / ``readzarr`` of one
    file by its extension, as a ``(z, y, x)`` array (an image as one slice). Partial files are
    read into one array, without the copy of the concatenation."""
    from petakit5d.io import get_image_dtype, get_image_size, read_tiff, read_zarr
    from petakit5d.jobs.files import fileparts

    if len(paths) == 1:
        ext = fileparts(paths[0])[2]
        if ext in (".tif", ".tiff"):
            vol = read_tiff(paths[0], num_threads=n, out=read_buffer)
        elif ext == ".zarr":
            vol = read_zarr(paths[0], num_threads=n, out=read_buffer)
        else:
            # MATLAB's switch has no other case: `frame` is never set
            raise MatlabError(
                "Unrecognized function or variable 'frame'.", "MATLAB:UndefinedFunction"
            )
        return vol if vol.ndim == 3 else vol[np.newaxis]
    sizes = [tuple(int(v) for v in get_image_size(p)) for p in paths]
    classes = [get_image_dtype(p) for p in paths]
    if any(s[1:3] != sizes[0][1:3] or len(s) != 3 for s in sizes):
        raise MatlabError(
            "Dimensions of arrays being concatenated are not consistent.",
            "MATLAB:catenate:dimensionMismatch",
        )
    if len(set(classes)) != 1:
        parts = [read_tiff(p, num_threads=n) for p in paths]
        return np.concatenate([p if p.ndim == 3 else p[np.newaxis] for p in parts], axis=0)
    from petakit5d.core.arrays import numpy_dtype

    shape = (sum(s[0] for s in sizes), sizes[0][1], sizes[0][2])
    dtype = numpy_dtype(classes[0])
    if read_buffer is not None:
        vol = read_buffer.view(shape, dtype)
    else:
        vol = np.empty(shape, dtype=dtype)
    z = 0
    for p, s in zip(paths, sizes, strict=True):
        read_tiff(p, num_threads=n, out=vol[z : z + s[0]])
        z += s[0]
    return vol


def _flat_field(ff_image: str, background_image: str, plane: tuple[int, int], lower_limit):
    """The flat field and background of ``XR_LSFlatFieldCorrection`` for planes ``plane``,
    prepared from the two tiff files (``flatfield.prepare_flat_field``).

    MATLAB reads and prepares them for every frame. They are kept here for the last files seen,
    keyed on each file's path, size and modification time (and the plane and the lower limit),
    and prepared again when any of them changes: the same images. The prepared images give the
    correction's result again with ``remove_ff_im_background`` and ``ls_rescale`` off (the crop
    is the identity on them, and ``max(LS, LowerLimit)`` gives the same values twice)."""
    from petakit5d.flatfield import prepare_flat_field
    from petakit5d.io import read_tiff

    def stat_key(p: str) -> tuple:
        st = os.stat(p)
        return (os.path.abspath(p), st.st_size, st.st_mtime_ns)

    try:
        key = (stat_key(ff_image), stat_key(background_image), plane, repr(lower_limit))
    except OSError:
        key = None
    if key is not None and key in _FF_CACHE:
        return _FF_CACHE[key]
    ls = read_tiff(ff_image)
    bg = read_tiff(background_image)
    prepared = prepare_flat_field(ls, bg, plane, lower_limit=lower_limit)
    _FF_CACHE.clear()
    if key is not None:
        _FF_CACHE[key] = prepared
    return prepared


# -- writing -----------------------------------------------------------------------------------


def _make_folder(path: str) -> None:
    """``mkdir`` and ``fileattrib(path, '+w', 'g')`` (not on Windows) of a missing folder"""
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
        if os.name != "nt":
            from petakit5d.io.tiff import group_writable

            group_writable(path)


def _write(data: np.ndarray, path: str, zarr: bool, block_size, private: bool, n: int) -> None:
    from petakit5d.io import write_tiff, write_zarr

    if zarr:
        write_zarr(data, path, block_size=block_size, uuid=not private, num_threads=n)
    else:
        write_tiff(data, path, num_threads=n)


def _empty_file(path: str) -> None:
    """``fclose(fopen(path, 'w'))``"""
    with open(path, "w"):
        pass


def _mip_z(vol: np.ndarray, save_16bit: bool, n: int) -> np.ndarray:
    """``max(vol, [], 3)``, made uint16 with ``save16bit`` (MATLAB's ``uint16``)"""
    from petakit5d.core import dtypes
    from petakit5d.tools.mip import axis_max

    mip = axis_max(vol, -3, num_threads=n)[0]
    if save_16bit and mip.dtype != np.uint16:
        mip = np.ascontiguousarray(dtypes.uint16(mip))
    return mip


def _to_output_class(vol: np.ndarray, save_16bit: bool, n: int) -> np.ndarray:
    """``uint16(vol)`` or ``single(vol)``"""
    want = np.uint16 if save_16bit else np.float32
    if vol.dtype == want:
        return vol
    return _to_class(vol, "uint16" if save_16bit else "single", n)


def _publish(tmp: str, final: str) -> None:
    from petakit5d.jobs.publish import publish_result

    publish_result(tmp, final)


# -- XR_deskewRotateFrame ----------------------------------------------------------------------


def _resample_completed(resample_factor) -> tuple[float, float, float] | None:
    """``[ones(1, 4 - numel(rs)) * rs(1), rs(2:end)]`` of MATLAB's ``[x y z]`` factor, in
    (z, y, x): one value for all axes, ``(z, xy)`` as ``(z, xy, xy)``"""
    if resample_factor is None or np.size(resample_factor) == 0:
        return None
    v = [float(x) for x in np.asarray(resample_factor, dtype=np.float64).reshape(-1)]
    if len(v) == 1:
        return (v[0], v[0], v[0])
    if len(v) == 2:
        return (v[0], v[1], v[1])
    return (v[0], v[1], v[2])


def _z_aniso(skew_angle: float, dz: float, px: float, objective_scan: bool, z_stage_scan: bool):
    from petakit5d.deskew_rotate.geometry import _cos, _sin

    theta = skew_angle * math.pi / 180
    if objective_scan:
        return dz / px
    if z_stage_scan:
        return _cos(abs(theta)) * dz / px
    return _sin(abs(theta)) * dz / px


def _rotate_out_size(shape, angle: float, z_aniso: float, objective_scan: bool):
    """the rotation's ``outSize`` of the separate path, ``(nz, ny, nx)``"""
    from petakit5d.deskew_rotate.geometry import _cos, _round, _sin

    nz, ny, nx = (float(v) for v in shape)
    t = angle * math.pi / 180
    c, s = _cos(t), _sin(abs(t))
    if not objective_scan:
        o = [ny, (nx - 1) * c + (nz - 1) * z_aniso / s, (nx - 1) * s - 4]
    else:
        o = [ny, nx * c + nz * z_aniso * s, nz * z_aniso * c + nx * s]
    o = [_round(v) for v in o]
    return (o[2], o[0], o[1])


def _deskew_blocks(frame: np.ndarray, angle, dz, px, reverse, crop, interp, n) -> np.ndarray:
    """MATLAB's ``blockedImage`` deskew of more than 1000 slices: blocks of
    ``min(ceil(ny / 8), 150)`` rows, each deskewed with 5 more rows on each side within the
    volume and cut back (``trimBorder``). The deskew moves voxels in x and z only, so the blocks
    give the rows of the deskew of the whole volume."""
    from petakit5d.deskew_rotate.volume import deskew_volume

    ny = frame.shape[1]
    rows = min(math.ceil(ny / 8), 150)
    out = None
    for y0 in range(0, ny, rows):
        y1 = min(ny, y0 + rows)
        b0, b1 = max(0, y0 - 5), min(ny, y1 + 5)
        block = _copy(frame[:, b0:b1, :], n)
        ds = deskew_volume(block, angle, dz, px, reverse, crop=crop, interp=interp)
        if out is None:
            out = np.empty((ds.shape[0], ny, ds.shape[2]), dtype=ds.dtype)
        out[:, y0:y1, :] = ds[:, y0 - b0 : y0 - b0 + (y1 - y0), :]
    return out


def _combined(
    frame,
    angle,
    dz,
    px,
    *,
    reverse,
    objective_scan,
    resample,
    interp,
    thresh,
    s16,
    n,
    keep_matlab_layout=False,
):
    """deskewRotateFrame3D on the frame: the kernels on MATLAB's layout where they run (the
    frame turned around in place and the result back, unless ``keep_matlab_layout``), else on
    (z, y, x). Returns the result and whether it is in MATLAB's layout ((z, x, y) C)."""
    from petakit5d.deskew_rotate import volume as dv
    from petakit5d.deskew_rotate.geometry import deskew_rotate_geometry

    kw = {
        "reverse": reverse,
        "objective_scan": objective_scan,
        "x_step_thresh": thresh,
        "resample_factor": resample,
        "save_16bit": s16,
        "interp_method": interp,
    }
    fast = False
    if frame.dtype in (np.uint16, np.float32) and isinstance(interp, str):
        try:
            g = deskew_rotate_geometry(
                frame.shape,
                angle,
                dz,
                px,
                reverse=reverse,
                objective_scan=objective_scan,
                x_step_thresh=thresh,
                resample_factor=resample,
                interp_method=interp,
            )
            fast = g.use_fast_method
        except (ValueError, IndexError):
            fast = False
    if not fast:
        return dv.deskew_rotate_volume(frame, angle, dz, px, crop=True, **kw), False
    vol_zxy = _turn_pages(frame, n)
    out = dv.deskew_rotate_volume_matlab_layout(vol_zxy, angle, dz, px, num_threads=n, **kw)
    if keep_matlab_layout:
        return out, True
    return _turn_pages(out, n), False


def _pre_resample(frame, resample, px, dz, angle):
    """MATLAB's binning in z and ``imresize`` in xy of more than 1000 slices before a combined
    deskew/rotation with resampling (``imresize`` is not ported yet)"""
    from petakit5d.deskew_rotate.geometry import _round
    from petakit5d.matlabcompat.elementary import sind

    rs = resample[::-1]  # [x y z]
    out_px = [r * px for r in rs]
    m = min(out_px)
    sd = float(np.asarray(sind(np.array([angle], dtype=np.float64))).reshape(-1)[0])
    pre = [m / px, m / px, m / (dz * sd)]
    pre[2] = max(1.0, _round(pre[2]))
    if all(p == 1 for p in pre):
        return frame, px, dz, resample
    raise NotImplementedError(
        "the combined deskew/rotation of more than 1000 slices with resampling bins and "
        "resizes the input first with MATLAB's imresize, which the port does not have yet"
    )


@register(
    "XR_deskewRotateFrame",
    Param("framePath", "@(x) ischar(x) || iscell(x)", "A file, or partial files of one volume."),
    Param("xyPixelSize", "@isscalar", "Pixel size in um."),
    Param("dz", "@isscalar", "Scan interval in um."),
    Param("DSDirName", "@ischar", "Folder of the deskewed results."),
    Param("DSRDirName", "@ischar", "Folder of the deskewed and rotated results."),
    Param("objectiveScan", "@islogical", "Objective scan."),
    Param("zStageScan", "@islogical", "Z stage scan (the angle becomes 90 - skewAngle)."),
    Param("overwrite", "@islogical", "Deskew again where the DS result exists."),
    Param("crop", "@islogical", "Crop the empty triangles of the separate deskew."),
    Param("skewAngle", "@isscalar", "Skew angle in degrees."),
    Param("reverse", "@islogical", "Reverse scan direction."),
    Param("inputAxisOrder", "@ischar", "Axis order of the input.", convert="axis_string"),
    Param("outputAxisOrder", "@ischar", "Axis order of the output.", convert="axis_string"),
    Param("rotate", "@islogical", "Rotate after the separate deskew."),
    Param("inputBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) of the input.", convert="bbox"),
    Param("flipZstack", "@islogical", "Flip the input along z."),
    Param("sCMOSCameraFlip", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("FFCorrection", "@islogical", "Flat-field correction."),
    Param("lowerLimit", "@isnumeric", "Lowest value of the normalized flat field."),
    Param("FFImage", "@ischar", "Flat-field image (tiff)."),
    Param("backgroundImage", "@ischar", "Camera background image (tiff)."),
    Param(
        "constOffset",
        "@(x) isnumeric(x)",
        "If set, a constant background instead of the camera's.",
        convert="same",
    ),
    Param("BKRemoval", "@islogical", "Camera background removal (known issue: MATLAB errors)."),
    Param("save16bit", "@islogical", "Save 16-bit results (else single).", convert="same"),
    Param("rescaleRotate", "@islogical", "Rescale the rotated data to [0 65535]."),
    Param("save3DStack", "@islogical", "Save the volumes (else empty result files)."),
    Param("saveMIP", "@islogical", "Save the z MIPs of the DS and DSR results."),
    Param("saveZarr", "@islogical", "Save as zarr."),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of zarr results.", convert="size_yxz"),
    Param("xStepThresh", "@isnumeric", "x step above which to interpolate in skewed space."),
    Param("zOffsetCorrection", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("DSRCombined", "@islogical", "Combined deskew, rotation and resampling."),
    Param(
        "resampleFactor",
        "@(x) isnumeric(x)",
        "Resampling after rotation (z, y, x).",
        convert="resample_xyz",
    ),
    Param("interpMethod", _INTERP, "Interpolation: 'linear' or 'cubic'."),
    Param("surffix", "@ischar", "Suffix of the result folders."),
    Param("uuid", "@ischar", "Id of the temporary results (empty: a new one, private)."),
    Param(
        "numThreads",
        None,
        "Threads of the reads, kernels and writes (None: the CPUs available).",
        py="num_threads",
        python_only=True,
    ),
    Param(
        "readBuffer",
        None,
        "A petakit5d.io.ReadBuffer to read the input into (reused over files).",
        py="read_buffer",
        python_only=True,
    ),
    task="D03",
)
def deskew_rotate_file(
    frame_path,
    xy_pixel_size: float,
    dz: float,
    *,
    ds_dir_name: str = "DS/",
    dsr_dir_name: str = "DSR/",
    objective_scan: bool = False,
    z_stage_scan: bool = False,
    overwrite: bool = False,
    crop: bool = False,
    skew_angle: float = 32.45,
    reverse: bool = False,
    input_axis_order: str = "zyx",
    output_axis_order: str = "zyx",
    rotate: bool = False,
    input_bbox: Sequence[int] | None = None,
    flip_zstack: bool = False,
    s_cmos_camera_flip: bool = False,
    ff_correction: bool = False,
    lower_limit: Any = 0.4,
    ff_image: str = "",
    background_image: str = "",
    const_offset: Any = None,
    bk_removal: bool = False,
    save_16bit: bool = True,
    rescale_rotate: bool = False,
    save_3d_stack: bool = True,
    save_mip: bool = True,
    save_zarr: bool = False,
    block_size: Sequence[int] = (500, 500, 500),
    x_step_thresh: float = 2,
    z_offset_correction: bool = False,
    dsr_combined: bool = True,
    resample_factor: Sequence[float] | None = None,
    interp_method: str = "linear",
    surffix: str = "",
    uuid: str = "",
    num_threads: int | None = None,
    read_buffer: Any = None,
) -> None:
    """Deskew and/or rotate one file, or the partial files of one volume
    (``XR_deskewRotateFrame``); see the module's description for the outputs and their order.

    ``frame_path``: a tiff or zarr file, or a list of tiff files concatenated along z.
    ``read_buffer``: a ``petakit5d.io.ReadBuffer`` the input is read into, so that a loop over
    files reuses its memory. A missing input warns and returns.
    """
    from petakit5d.core.axes import axis_permutation
    from petakit5d.core.ids import get_uuid
    from petakit5d.jobs.files import exist, fileparts

    del s_cmos_camera_flip, z_offset_correction
    n = _threads(num_threads)
    paths = (
        [os.fspath(frame_path)]
        if isinstance(frame_path, str | os.PathLike)
        else [os.fspath(p) for p in frame_path]
    )
    objective_scan, z_stage_scan = _logical(objective_scan), _logical(z_stage_scan)
    reverse, save_16bit = _logical(reverse), _logical(save_16bit)
    dsr_combined, rotate = _logical(dsr_combined), _logical(rotate)
    save_zarr, save_mip = _logical(save_zarr), _logical(save_mip)
    save_3d_stack = _logical(save_3d_stack)
    private = not uuid
    if not uuid:
        uuid = get_uuid()
    # with a full random uuid made here (shorter on Windows) no other process writes the
    # temporary results: zarr chunks are written without temporary files
    private = private and len(uuid) >= 32
    px, dz, skew_angle = float(xy_pixel_size), float(dz), float(skew_angle)
    z_aniso = _z_aniso(skew_angle, dz, px, objective_scan, z_stage_scan)
    if z_stage_scan:
        angle = 90 - skew_angle
        reverse = not reverse
    else:
        angle = skew_angle
    box = None if input_bbox is None or np.size(input_bbox) == 0 else tuple(input_bbox)
    block = tuple(int(v) for v in np.asarray(block_size).reshape(-1))

    for p in paths:
        if not exist(p):
            warnings.warn(f"{p} does not exist!", stacklevel=2)
            return
    rt, fsname, _ = fileparts(paths[0])
    ext = ".zarr" if save_zarr else ".tif"
    ds_path = ds_fullname = ""
    if not dsr_combined:
        ds_path = f"{rt}/{ds_dir_name}{surffix}/"
        _make_folder(ds_path)
        ds_fullname = f"{ds_path}{fsname}{ext}"

    frame = ds = dsr = None
    in_shape = None
    if dsr_combined or not exist(ds_fullname) or _logical(overwrite):
        t0 = time.perf_counter()
        frame = _read_frame(paths, n, read_buffer)
        print(f"Image read time: {time.perf_counter() - t0:0.6f} s", flush=True)
        frame = _preprocess(
            frame,
            paths[0],
            box,
            flip_zstack,
            ff_correction,
            ff_image,
            background_image,
            lower_limit,
            const_offset,
            bk_removal,
            input_axis_order,
            n,
        )
        in_shape = frame.shape

        if not dsr_combined:
            print(f"Deskew frame {paths[0]}...", flush=True)
            if _logical(overwrite) and exist(ds_fullname):
                from petakit5d.jobs.publish import remove_path

                remove_path(ds_fullname)
            t0 = time.perf_counter()
            split = in_shape[0] > 1000
            fr = frame if frame.dtype == np.float32 else _to_class(frame, "single", n)
            frame = None
            if not split:
                from petakit5d.deskew_rotate.volume import deskew_volume

                ds = deskew_volume(
                    fr, angle, dz, px, reverse, crop=_logical(crop), interp=interp_method
                )
            else:
                print("Use image block method for deskew...", flush=True)
                ds = _deskew_blocks(fr, angle, dz, px, reverse, _logical(crop), interp_method, n)
            fr = None
            print(f"Deskew processing time: {time.perf_counter() - t0:0.6f} s", flush=True)
            if save_mip:
                mip_path = f"{ds_path}/MIPs/"
                _make_folder(mip_path)
                mip = _mip_z(ds, save_16bit, n)
                _write(mip, f"{mip_path}{fsname}_MIP_z{ext}", save_zarr, (500, 500, 500), False, n)
                del mip
            t0 = time.perf_counter()
            tmp = f"{ds_path}{fsname}_{uuid}{ext}"
            if save_3d_stack:
                _write(_to_output_class(ds, save_16bit, n), tmp, save_zarr, block, private, n)
            else:
                _empty_file(tmp)
            _publish(tmp, ds_fullname)
            if split:
                # MATLAB keeps only the blocks on disk: the rotation reads DS back
                ds = None
            print(f"Deskew result write time: {time.perf_counter() - t0:0.6f} s", flush=True)

    if not (rotate or dsr_combined):
        return
    dsr_path = f"{rt}/{dsr_dir_name}{surffix}/"
    _make_folder(dsr_path)
    dsr_fullname = f"{dsr_path}{fsname}{ext}"
    resample = _resample_completed(resample_factor)
    done = os.path.isdir(dsr_fullname) if save_zarr else exist(dsr_fullname)
    if done:
        return
    t0 = time.perf_counter()
    ml_layout = False
    if not dsr_combined:
        print(f"rotate frame {paths[0]}...", flush=True)
        if ds is None:
            ds, in_shape = _read_ds(ds_fullname, paths, box, input_axis_order, save_zarr, n)
        out_size = _rotate_out_size(in_shape, angle, z_aniso, objective_scan)
        from petakit5d.deskew_rotate.volume import rotate_volume

        dsr = rotate_volume(
            ds,
            angle,
            z_aniso,
            reverse,
            crop=True,
            resample=resample,
            objective_scan=objective_scan,
            out_size=out_size,
            interp=interp_method,
        )
        ds = None
        print(f"Rotation processing time: {time.perf_counter() - t0:0.6f} s", flush=True)
    else:
        if resample is not None and any(r != 1 for r in resample) and frame.shape[0] > 1000:
            frame, px, dz, resample = _pre_resample(frame, resample, px, dz, angle)
        if not save_16bit and frame.dtype != np.uint16 and frame.dtype != np.float32:
            frame = _to_class(frame, "single", n)
        print(f"Deskew, rotate and resample for frame {paths[0]}...", flush=True)
        # a zarr result is written from the kernels' layout through a (z, y, x) view (cpp-zarr
        # writes strided data at about the cost of contiguous data), saving the page turn
        keep_ml = (
            save_zarr
            and axis_permutation("zyx", output_axis_order) == (0, 1, 2)
            and not _logical(rescale_rotate)
        )
        dsr, ml_layout = _combined(
            frame,
            angle,
            dz,
            px,
            reverse=reverse,
            objective_scan=objective_scan,
            resample=resample,
            interp=interp_method,
            thresh=float(x_step_thresh),
            s16=save_16bit,
            n=n,
            keep_matlab_layout=keep_ml,
        )
        frame = None
        print(
            f"Combined deskew/rotation processing time: {time.perf_counter() - t0:0.6f} s",
            flush=True,
        )

    perm = axis_permutation("zyx", output_axis_order)
    if perm != (0, 1, 2):
        dsr = _copy(np.transpose(dsr, perm), n)

    if save_mip:
        mip_path = f"{dsr_path}/MIPs/"
        _make_folder(mip_path)
        mip = _mip_z(dsr, save_16bit, n)
        if ml_layout:
            mip = np.ascontiguousarray(mip.T)  # (x, y) -> (y, x)
        final = f"{mip_path}{fsname}_MIP_z{ext}"
        tmp = f"{final}_{uuid}" if save_zarr else f"{mip_path}{fsname}_MIP_z.tiff_{uuid}"
        _write(mip, tmp, save_zarr, (500, 500, 500), False, n)
        if save_zarr:
            _publish(tmp, final)
        else:
            os.replace(tmp, final)
        del mip

    t0 = time.perf_counter()
    tmp = f"{dsr_path}{fsname}_{uuid}{ext}"
    if _logical(rescale_rotate):
        from petakit5d.matlabcompat import reduce

        lo, _ = reduce.min(dsr)
        hi, _ = reduce.max(dsr)
        rng = (np.asarray(lo, dtype=dsr.dtype), np.asarray(hi, dtype=dsr.dtype))
        dsr = scale_contrast(dsr, rng, (0, 65535), num_threads=n)
    dsr = _to_output_class(dsr, save_16bit, n)
    if ml_layout:
        dsr = dsr.transpose(0, 2, 1)  # the (z, y, x) view of MATLAB's layout
    if save_3d_stack:
        _write(dsr, tmp, save_zarr, block, private, n)
    else:
        _empty_file(tmp)
    _publish(tmp, dsr_fullname)
    print(f"Deskew/rotation result write time: {time.perf_counter() - t0:0.6f} s", flush=True)


def _preprocess(
    frame,
    path,
    box,
    flip_zstack,
    ff_correction,
    ff_image,
    background_image,
    lower_limit,
    const_offset,
    bk_removal,
    input_axis_order,
    n,
):
    """the input box, the flip, the flat field, the background removal and the axis order, in
    MATLAB's order"""
    from petakit5d.core.axes import axis_permutation

    if box is not None:
        from petakit5d.core.bbox import box_to_matlab
        from petakit5d.kernels.crop import crop3d_matlab

        frame = crop3d_matlab(frame, box_to_matlab(box), num_threads=n)
        frame = frame if frame.ndim == 3 else frame[np.newaxis]
    if _logical(flip_zstack):
        frame = _copy(frame[::-1], n)
    if _logical(ff_correction):
        from petakit5d.flatfield import flatfield_correct_volume

        print(
            f"Flat-field correction for frame {path}...\n  Flat-field image: {ff_image}\n"
            f"  Background image: {background_image}",
            flush=True,
        )
        ls, bg = _flat_field(ff_image, background_image, frame.shape[1:], lower_limit)
        frame = flatfield_correct_volume(
            frame,
            ls,
            bg,
            lower_limit=lower_limit,
            const_offset=const_offset,
            remove_ff_im_background=False,
            ls_rescale=False,
            num_threads=n,
        )
        frame = frame if frame.ndim == 3 else frame[np.newaxis]
    if _logical(bk_removal):
        from petakit5d.io import read_tiff

        read_tiff(background_image)
        # KNOWN ISSUE (kept): XR_deskewRotateFrame calls XR_CameraBackgroundRemoval, which the
        # repository does not have
        raise MatlabError(
            "Unrecognized function or variable 'XR_CameraBackgroundRemoval'.",
            "MATLAB:UndefinedFunction",
        )
    perm = axis_permutation(input_axis_order, "zyx")
    if perm != (0, 1, 2):
        frame = _copy(np.transpose(frame, perm), n)
    return frame


def _read_ds(ds_fullname, paths, box, input_axis_order, save_zarr, n):
    """the deskewed volume read back as single, and the size the deskew had (after the box, in
    (z, y, x)), as MATLAB rotates a DS result made earlier"""
    from petakit5d.core.axes import axis_permutation
    from petakit5d.io import get_image_size, read_tiff, read_zarr

    if save_zarr:
        ds = read_zarr(ds_fullname, num_threads=n)
    else:
        ds = read_tiff(ds_fullname, num_threads=n)
    ds = ds if ds.ndim == 3 else ds[np.newaxis]
    if ds.dtype != np.float32:
        ds = _to_class(ds, "single", n)
    sz = list(get_image_size(paths[0]))[:3]
    for p in paths[1:]:
        sz[0] += int(get_image_size(p)[0])
    if box is not None:
        sz = [box[i + 3] - box[i] for i in range(3)]
    perm = axis_permutation(input_axis_order, "zyx")
    return ds, tuple(int(sz[i]) for i in perm)
