"""Tiff and zarr conversion: ``tiffToZarr``, ``XR_tiffToZarr_wrapper``, ``zarrToTiff`` and
``XR_zarrToTiff_wrapper``.

Arrays are ``(z, y, x)`` C order; boxes ``(z0, y0, x0, z1, y1, x1)``, 0-based, end-exclusive;
block, shard and resample sizes ``(z, y, x)`` (a resample factor of two values is
``(z, xy)``). The MATLAB-name layer converts MATLAB's values. The zarr arrays written are the
port's native ones (``[nz, ny, nx]``, order C); tiff and zarr inputs of both layouts are read.

``tiff_to_zarr`` processes a whole image as an array, as MATLAB does since 89b7790, in
MATLAB's order: read (an ``input_bbox`` reads only its slices, then crops y and x), flip in z,
crop to ``tile_out_bbox``, resample, the user function, then the cast back to the input's
class (MATLAB's rounding and saturation). With ``read_whole_tiff`` false (and none of flip,
resampling or a user function) the tiff is read by z blocks of the zarr's chunk size, as
MATLAB's blockedImage loop reads it, without blockedImage: one read of the block's slices, a
crop of y and x, one region write per block.

Hooks:

- **user functions** (``usr_fcn``; the dataset wrapper's ``process_fun_path`` files): a Python
  callable, ``'package.module:function'``, or a record as JSON text
  ``'{"func": "package.module:function", "args": [...], "kwargs": {...}}'`` that calls
  ``func(volume, *args, **kwargs)`` (text, so it goes into job tasks and function files as
  MATLAB's code strings did). MATLAB code (``'@(x)...'``) is handed to the
  stitching's tile-function parser (``petakit5d.stitch.tile_functions.parse_tile_function``,
  task S02), and fails with ``NotImplementedError`` until that exists. The flat-field
  correction is such a function (D01's ``petakit5d.flatfield.flatfield_correct_file``,
  MATLAB's ``processFFCorrectionFrame``).
- **resampling** goes to ``petakit5d.tools.resample.imresize3_fast(vol, out_size)`` (T04's
  port of ``imresize3_fast``, cubic), and fails with ``NotImplementedError`` until it exists.

Results are written under ``<name>_<uuid>`` and moved into place with one rename; a result
another writer published first is kept (``jobs.publish``), also with ``overwrite`` (MATLAB
then deleted it). Until cpp-zarr releases the dev branch's options, the port's zarr writes
skip all-zero chunks (MATLAB writes them; the data read back is the same) and always write
each chunk to a temporary file (MATLAB writes a temporary result its call named in place,
``uuid`` false), as ``petakit5d.io.zarr`` describes.
"""

from __future__ import annotations

import math
import os
import shutil
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "ConvertError",
    "tiff_to_zarr",
    "tiff_to_zarr_dataset",
    "user_function",
    "zarr_to_tiff",
    "zarr_to_tiff_dataset",
]

# memory factors of XR_tiffToZarr_wrapper: no function, the flat-field function, others
_MEM_FACTOR_PLAIN = 2.5
_MEM_FACTOR_FLATFIELD = 2.75
_MEM_FACTOR_FUNCTION = 4
_FLATFIELD_NAMES = ("processFFCorrectionFrame", "flatfield_correct_file")
# copies of strided views from this size are made in threads
_THREADED_BYTES = 1 << 22


class ConvertError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _m2d(a: np.ndarray) -> np.ndarray:
    """MATLAB has no trailing singleton: a volume of one slice is an image (a view)"""
    return a[0] if a.ndim == 3 and a.shape[0] == 1 else a


def _msize(shape: Sequence[int]) -> list[int]:
    """MATLAB's size of a (z, y, x) or (y, x) shape: [ny nx nz] or [ny nx]"""
    s = [int(v) for v in shape]
    return s[1:] + s[:1] if len(s) == 3 else s


def _pysize(msize: Sequence[int]) -> tuple[int, ...]:
    """a MATLAB size [y x z] (or [y x]) as (z, y, x) (or (y, x))"""
    s = [int(v) for v in msize]
    return tuple(s[2:] + s[:2]) if len(s) == 3 else tuple(s)


def _mmin(a: Sequence[float], b: Sequence[float]) -> list[float]:
    """MATLAB's min of two vectors (a scalar expands; other sizes must match)"""
    a, b = list(a), list(b)
    if len(a) == 1:
        a = a * len(b)
    if len(b) == 1:
        b = b * len(a)
    if len(a) != len(b):
        raise ConvertError(
            "Arrays have incompatible sizes for this operation.", "MATLAB:sizeDimensionsMustMatch"
        )
    return [min(x, y) for x, y in zip(a, b, strict=True)]


def _mround(x: float) -> int:
    from petakit5d.core.dtypes import mround

    return int(mround(x))


def _materialize(view: np.ndarray) -> np.ndarray:
    """a C-contiguous copy of a strided view (flips and crops), planes copied in threads"""
    if view.flags.c_contiguous:
        return view
    out = np.empty(view.shape, dtype=view.dtype)
    if view.nbytes < _THREADED_BYTES or view.shape[0] < 2:
        out[...] = view
        return out
    from concurrent.futures import ThreadPoolExecutor

    from petakit5d.core.threads import resolve_threads

    n = min(resolve_threads(None), view.shape[0])
    edges = np.linspace(0, view.shape[0], 4 * n + 1).astype(int)

    def copy(k: int) -> None:
        out[edges[k] : edges[k + 1]] = view[edges[k] : edges[k + 1]]

    with ThreadPoolExecutor(max_workers=n) as ex:
        list(ex.map(copy, range(4 * n)))
    return out


def _crop_yx(a: np.ndarray, box: Sequence[int]) -> np.ndarray:
    """``crop3d(I, [ymin xmin 1 ymax xmax size(I, 3)])``: the y and x range of a box, all z"""
    from petakit5d.kernels.crop import crop_3d

    b = [int(v) for v in box]
    if a.ndim == 2:
        return crop_3d(a, (b[1], b[2], b[4], b[5]))
    return crop_3d(a, (0, b[1], b[2], a.shape[0], b[4], b[5]))


def _crop(a: np.ndarray, box: Sequence[int]) -> np.ndarray:
    """MATLAB's crop3d of a volume or an image with a box (MATLAB's error outside the data)"""
    from petakit5d.core.bbox import box_slices, check_box

    b = list(check_box(box))
    if a.ndim == 2:
        if len(b) == 6:
            # a 3D box on an image: its z range must be the one slice
            if b[0] != 0 or b[3] not in (1, None):
                raise ConvertError(
                    "The bounding box is out of the bound of the input array!", "crop:inputError"
                )
            b = [b[1], b[2], b[4], b[5]]
    elif len(b) == 4:
        b = [0, b[0], b[1], a.shape[0], b[2], b[3]]
    try:
        check_box(b, a.shape)
    except ValueError as err:
        raise ConvertError(
            "The bounding box is out of the bound of the input array!", "crop:inputError"
        ) from err
    return a[box_slices(b)]


def _ints(v: Any) -> list[int]:
    return [int(x) for x in np.asarray(v, dtype=np.float64).ravel()]


def _empty(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, str):
        return v == ""
    return np.size(v) == 0


def _print(msg: str) -> None:
    print(msg, flush=True)


def _publish(tmp: str, final: str, kind: str) -> None:
    """MATLAB's last step: a result that exists now (another writer's) is kept and the
    temporary one removed, else the temporary result is moved into place (one rename)"""
    from petakit5d.jobs.publish import publish_result, remove_path

    if os.path.lexists(final):
        _print(f"{kind} result {final} already exists, skip it!")
        remove_path(tmp)
        return
    publish_result(tmp, final)


# ---------------------------------------------------------------------------------------------
# user functions and resampling (hooks)
# ---------------------------------------------------------------------------------------------


def _import(name: str) -> Callable[..., Any]:
    import importlib

    if name.count(":") != 1:
        raise ValueError(
            f"a user function is 'package.module:function' (or a callable), not {name!r}"
        )
    module, qual = name.split(":")
    obj: Any = importlib.import_module(module)
    for part in qual.split("."):
        obj = getattr(obj, part)
    if not callable(obj):
        raise TypeError(f"{name!r} is not callable")
    return obj


def _matlab_code(text: str) -> Callable[..., Any]:
    """MATLAB code as a user function: through the stitching's tile-function parser (S02)"""
    try:
        from petakit5d.stitch.tile_functions import parse_tile_function  # type: ignore
    except ImportError:
        err = NotImplementedError(
            f"MATLAB code as a user function ({text!r}) needs the tile-function parser "
            "(petakit5d.stitch.tile_functions), which is not there; give a Python function "
            "('package.module:function' or a callable)"
        )
        err.identifier = "petakit5d:convert:matlabCodeFunction"  # type: ignore[attr-defined]
        raise err from None
    return parse_tile_function(text)


def user_function(spec: Any) -> Callable[[np.ndarray], np.ndarray] | None:
    """The function a user-function setting names, or None for none (``''``, None).

    ``spec``: a callable; ``'package.module:function'``; a record ``{"func": ..., "args": [...],
    "kwargs": {...}}`` (``func(volume, *args, **kwargs)``), as a mapping or as its JSON text;
    MATLAB code (``'@(x)...'``, through the stitching's tile-function parser)."""
    if spec is None or callable(spec):
        return spec
    if isinstance(spec, str) and spec.strip().startswith("{"):
        import json

        spec = json.loads(spec)
        if not isinstance(spec, Mapping) or "func" not in spec:
            raise ValueError(f"a user-function record has a 'func', not {spec!r}")
    if isinstance(spec, Mapping):
        func = _import(str(spec["func"]))
        args = list(spec.get("args", ()))
        kwargs = dict(spec.get("kwargs", {}))
        return lambda vol: func(vol, *args, **kwargs)
    if isinstance(spec, str):
        text = spec.strip()
        if not text:
            return None
        if text.startswith("@"):
            return _matlab_code(text)
        return _import(text)
    if isinstance(spec, Sequence) and len(spec) == 1:
        return user_function(spec[0])
    raise TypeError(f"a user function is a callable, an import path or a record, not {spec!r}")


def _function_text(spec: Any) -> str:
    """the setting as MATLAB displays it (disp(usrFcn)) and as the memory rule reads it"""
    if isinstance(spec, str):
        return spec
    if isinstance(spec, Mapping):
        import json

        return json.dumps(dict(spec))
    name = getattr(spec, "__qualname__", None) or repr(spec)
    return f"{getattr(spec, '__module__', '')}:{name}"


def _resize(vol: np.ndarray, out_size: Sequence[int]) -> np.ndarray:
    """``imresize3_fast(I, outSize)`` (cubic), T04's"""
    try:
        from petakit5d.tools.resample import imresize3_fast  # type: ignore
    except ImportError:
        err = NotImplementedError(
            "resampling (resample_factor) needs petakit5d.tools.resample.imresize3_fast "
            "(task T04), which is not there yet"
        )
        err.identifier = "petakit5d:convert:resampleNotPorted"  # type: ignore[attr-defined]
        raise err from None
    return imresize3_fast(vol, tuple(int(v) for v in out_size))


def _resample_factors(resample_factor: Any) -> list[float] | None:
    """(z, y, x) factors (MATLAB: [r r r] for one value, [xy xy z] for [xy z]); None for none
    or all ones"""
    if _empty(resample_factor):
        return None
    rs = [float(v) for v in np.asarray(resample_factor, dtype=np.float64).ravel()]
    if all(v == 1 for v in rs):
        return None
    if len(rs) == 1:
        return rs * 3
    if len(rs) == 2:  # (z, xy)
        return [rs[0], rs[1], rs[1]]
    return rs[:3] if len(rs) == 3 else rs


# ---------------------------------------------------------------------------------------------
# tiffToZarr
# ---------------------------------------------------------------------------------------------


def _read_whole(paths: list[str], input_bbox: Any, partial: bool) -> np.ndarray:
    from petakit5d.io import read_tiff

    if not partial:
        if not _empty(input_bbox):
            b = _ints(input_bbox)
            a = read_tiff(paths[0], (b[0], b[3]))
            return _m2d(_crop_yx(a, b))
        return read_tiff(paths[0])
    parts = [read_tiff(p) for p in paths]
    parts = [p if p.ndim == 3 else p[np.newaxis] for p in parts]
    a = parts[0] if len(parts) == 1 else np.concatenate(parts, axis=0)
    if not _empty(input_bbox):
        a = _materialize(_crop(a, _ints(input_bbox)))
    return _m2d(a)


@register(
    "tiffToZarr",
    Param("tifFilename", "@(x) ischar(x) || iscell(x)", "Tiff file, or its partial files."),
    Param("zarrFilename", "@ischar", "Zarr result."),
    Param(
        "frame",
        "@isnumeric",
        "Data in memory to convert instead (None: read the tiff).",
        convert="volume",
    ),
    Param("Overwrite", "@islogical", "Replace an existing result."),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x).", convert="size_yxz"),
    Param(
        "shardSize", "@isnumeric", "Shard size (z, y, x) (Zarr v3 with shards).", convert="size_yxz"
    ),
    Param("expand2dDim", "@islogical", "Store an image as one z slice ([1, ny, nx])."),
    Param("flipZstack", "@islogical", "Flip the volume in z."),
    Param(
        "resampleFactor",
        "@(x) isempty(x) || isnumeric(x)",
        "Resampling factors (z, y, x) or (z, xy) (None: none).",
        convert="resample_yxz",
    ),
    Param(
        "inputBbox", "@(x) isnumeric(x)", "Region (z0, y0, x0, z1, y1, x1) to read.", convert="bbox"
    ),
    Param(
        "tileOutBbox",
        "@(x) isempty(x) || isnumeric(x)",
        "Region (z0, y0, x0, z1, y1, x1) to keep after the flip.",
        convert="bbox",
    ),
    Param("readWholeTiff", "@islogical", "Read the tiff at once (else by z blocks)."),
    Param("compressor", "@ischar", "Zarr compressor."),
    Param(
        "usrFcn",
        "@(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x)",
        "User function applied last (callable, import path or record).",
    ),
    Param("uuid", "@ischar", "Id of the temporary result (empty: a new one)."),
    task="T01",
)
def tiff_to_zarr(
    tif_filename: str | Sequence[str],
    zarr_filename: str,
    frame: np.ndarray | None = None,
    *,
    overwrite: bool = False,
    block_size: Sequence[int] = (250, 500, 500),
    shard_size: Sequence[int] | None = None,
    expand_2d_dim: bool = True,
    flip_zstack: bool = False,
    resample_factor: Sequence[float] | float | None = None,
    input_bbox: Sequence[int] | None = None,
    tile_out_bbox: Sequence[int] | None = None,
    read_whole_tiff: bool = True,
    compressor: str = "zstd",
    usr_fcn: Any = "",
    uuid: str = "",
) -> None:
    """Convert a tiff file (its partial files joined along z, or ``frame``, a ``(z, y, x)``
    array) to a native zarr array (``tiffToZarr``), with the flip in z, the crop boxes,
    resampling and the user function in MATLAB's order, and the result cast back to the input's
    class. An existing result is kept (``overwrite``: replaced first)."""
    from petakit5d.core.axes import to_yxz
    from petakit5d.core.dtypes import cast
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import create_zarr, get_image_dtype, get_image_size, read_tiff, write_zarr

    zarr_filename = os.fspath(zarr_filename).rstrip("/")
    if os.path.isdir(zarr_filename):
        if not overwrite:
            _print(f"Zarr result {zarr_filename} already exists, skip it!")
            return
        shutil.rmtree(zarr_filename)

    given_uuid = uuid
    if not uuid:
        uuid = get_uuid()
    # with a full random uuid made here no other process writes the temporary result
    private_tmp = not given_uuid and len(uuid) >= 32

    if isinstance(tif_filename, (str, os.PathLike)):
        files = [] if os.fspath(tif_filename) == "" else [os.fspath(tif_filename)]
    else:
        files = [os.fspath(f) for f in tif_filename]
    have_frame = frame is not None and np.size(frame) > 0
    if not files:
        _print(f"Convert frame to Zarr for {zarr_filename} ...")
    else:
        _print(f"Convert Tiff to Zarr for {files[0]} ...")

    fcn = user_function(usr_fcn)
    if flip_zstack or _resample_factors(resample_factor) is not None or fcn is not None:
        read_whole_tiff = True

    # MATLAB's sizes and block size, [y x z] ([y x] for an image)
    bs = list(to_yxz(block_size)) if len(_ints(block_size)) == 3 else _ints(block_size)
    im = None
    region = None  # by z blocks: the file, its box (z0, y0, x0, z1, y1, x1)
    if have_frame:
        im = _m2d(np.asarray(frame))
        if im.ndim == 2:
            bs = bs[:2]
        bs = _mmin(bs, _msize(im.shape))
    elif len(files) == 1:
        if read_whole_tiff:
            im = _read_whole(files, input_bbox, False)
            sz = _msize(im.shape)
            if im.ndim == 2:
                if expand_2d_dim:
                    bs = (bs + [1])[:3]
                    bs[2] = 1
                    sz = [*sz, 1]
                else:
                    bs = bs[:2]
            bs = _mmin(bs, sz)
        else:
            full = [int(v) for v in get_image_size(files[0])[:3]]
            box = [0, 0, 0, *full]
            if not _empty(input_bbox):
                box = _sub_box(box, _ints(input_bbox))
            region = {"path": files[0], "box": box, "dtype": get_image_dtype(files[0])}
    elif not files:
        # (MATLAB: readtiff('') fails)
        raise FileNotFoundError('File "" cannot be opened')
    else:
        im = _read_whole(files, input_bbox, True)
        if im.ndim == 2:
            bs = bs[:2]
        bs = _mmin(bs, _msize(im.shape))

    if region is None:
        sz = _msize(im.shape)
        dtype = _class(im)
    else:
        b = region["box"]
        sz = _msize([b[3] - b[0], b[4] - b[1], b[5] - b[2]])
        dtype = region["dtype"]
    if any(s == 0 for s in sz):
        _print("The input image is empty, skip it!")
        return

    view = im
    if flip_zstack and view.ndim == 3:
        view = view[::-1]
    if not _empty(tile_out_bbox):
        if region is None:
            view = _crop(view, _ints(tile_out_bbox))
        else:
            region["box"] = _sub_box(region["box"], _ints(tile_out_bbox))
    if view is not None:
        im = _m2d(_materialize(view))
    del view

    rs = _resample_factors(resample_factor)
    if rs is not None:
        if im.ndim != 3:
            raise ConvertError(
                "Arrays have incompatible sizes for this operation.", "MATLAB:dimagree"
            )
        out_size = [_mround(n / r) for n, r in zip(im.shape, rs, strict=True)]
        im = _resize(im, out_size)
    if region is None:
        sz = _msize(im.shape)
    else:
        b = region["box"]
        sz = _msize([b[3] - b[0], b[4] - b[1], b[5] - b[2]])

    if fcn is not None:
        _print(_function_text(usr_fcn))
        im = np.asarray(fcn(im))
    if len(sz) == 2:
        if expand_2d_dim:
            sz = [*sz, 1]
            bs = (list(bs) + [1])[:3]
            bs[2] = 1
        else:
            prod = math.prod(bs)
            if sz[0] <= sz[1]:
                bs = _mmin(sz, [sz[0], _mround(prod / sz[0])])
            else:
                bs = _mmin(sz, [_mround(prod / sz[1]), sz[1]])
    else:
        bs = _mmin(sz, bs)
    bs = [int(v) for v in bs]

    tmp = f"{zarr_filename}_{uuid}"
    if not os.path.isdir(tmp):
        nchunks = math.prod(math.ceil(s / b) for s, b in zip(sz, bs, strict=True))
        create_zarr(
            tmp,
            data_size=_pysize(sz),
            block_size=_pysize(bs),
            shard_size=None if _empty(shard_size) else shard_size,
            dtype=dtype,
            expand_2d_dim=expand_2d_dim,
            compressor=compressor,
            dim_separator="/" if nchunks > 10000 else ".",
        )

    if region is None:
        data = _materialize(im if _class(im) == dtype else np.asarray(cast(im, dtype)))
        del im
        write_zarr(data, tmp, block_size=_pysize(bs), uuid=not private_tmp)
        del data
    else:
        z0, y0, x0, z1, y1, x1 = region["box"]
        nz, bz = z1 - z0, bs[2]
        for i in range(math.ceil(nz / bz)):
            zs, zt = i * bz, min((i + 1) * bz, nz)
            block = read_tiff(region["path"], (z0 + zs, z0 + zt))
            block = _crop_yx(block, (0, y0, x0, 0, y1, x1))
            write_zarr(
                block,
                tmp,
                bbox=(zs, 0, 0, zt, y1 - y0, x1 - x0),
                uuid=not private_tmp,
            )
            del block

    # (a result another writer published meanwhile is kept, also with overwrite)
    _publish(tmp, zarr_filename, "Zarr")


def _sub_box(box: Sequence[int], sub: Sequence[int]) -> list[int]:
    """the box ``sub`` (relative to ``box``'s start) inside ``box``, as blockedImage's crop of a
    cropped image; outside it is MATLAB's error"""
    b, s = [int(v) for v in box], [int(v) for v in sub]
    if len(s) != 6:
        raise ConvertError("The bounding box must have 6 values.", "")
    ext = [b[3] - b[0], b[4] - b[1], b[5] - b[2]]
    if any(s[i] < 0 or s[i + 3] > ext[i] or s[i] >= s[i + 3] for i in range(3)):
        raise ConvertError(
            "Expected input number 2, PixelStart, to be within the image.",
            "images:blockedImage:pixelSubscriptOutOfBounds",
        )
    return [b[i] + s[i] for i in range(3)] + [b[i] + s[i + 3] for i in range(3)]


def _class(a: np.ndarray) -> str:
    from petakit5d.core.arrays import matlab_class

    return matlab_class(a.dtype)


# ---------------------------------------------------------------------------------------------
# XR_tiffToZarr_wrapper
# ---------------------------------------------------------------------------------------------


def _bbox_rows_in(v: Any) -> Any:
    """one box, or one box per tile location (rows)"""
    from petakit5d.core.bbox import box_from_matlab

    if _empty(v):
        return None
    a = np.asarray(v, dtype=np.float64)
    if a.ndim == 2 and a.shape[0] > 1:
        return tuple(box_from_matlab(r) for r in a)
    return box_from_matlab(a.ravel())


def _bbox_rows_out(v: Any) -> Any:
    from petakit5d.core.bbox import box_to_matlab

    if _empty(v):
        return None
    if isinstance(v[0], (Sequence, np.ndarray)) and not isinstance(v[0], str):
        return tuple(tuple(box_to_matlab(r)) for r in v)
    return box_to_matlab(v)


_BBOX_ROWS = (_bbox_rows_in, _bbox_rows_out)


def _box_rows(v: Any) -> list[Any]:
    """the rows of a box setting (``size(x, 1)``): [] for none, [box] for one box"""
    if _empty(v):
        return []
    if isinstance(v[0], (Sequence, np.ndarray)) and not isinstance(v[0], str):
        return [None if _empty(r) else list(r) for r in v]
    return [list(v)]


def _fun_rows(process_fun_path: Any) -> list[list[str]]:
    """process_fun_path as MATLAB's cell: rows (channels) of columns (keyword alternatives); a
    plain sequence is one function file per channel"""
    if _empty(process_fun_path):
        return []
    if isinstance(process_fun_path, (str, os.PathLike)):
        return [[os.fspath(process_fun_path)]]
    rows = []
    for r in process_fun_path:
        if isinstance(r, (str, os.PathLike)) or r is None:
            rows.append(["" if r is None else os.fspath(r)])
        else:
            rows.append(["" if e is None else os.fspath(e) for e in r])
    return rows


def _read_function_file(path: str) -> tuple[str, Any]:
    """(usrFun, keywords or None) of a .mat (variables usrFun and keywords) or .txt (first
    line) function file"""
    from petakit5d.core.paths import fileparts

    ext = fileparts(path)[2]
    if ext == ".mat":
        from petakit5d.io import read_mat

        a = read_mat(path)
        fun = a["usrFun"]
        kw = a.get("keywords")
    elif ext == ".txt":
        from petakit5d.core.textfiles import read_text_file

        lines = read_text_file(path)
        fun, kw = lines[0], None
        if fun == -1:
            fun = ""
    else:
        raise ConvertError("Unknown file type for the processed function file!")
    if isinstance(fun, np.ndarray):
        fun = "".join(fun.ravel().tolist()) if fun.size else ""
    return str(fun), kw


def _keywords(kw: Any) -> list[str]:
    if kw is None:
        return [""]
    if isinstance(kw, str):
        return [kw]
    if isinstance(kw, np.ndarray):
        return [str(e) for e in kw.ravel(order="F")]
    return [str(e) for e in kw]


def _save_parameters(func: Callable[..., Any], result_path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) and parameters.json of the call (MATLAB's ``pr``)"""
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(func, values).items()}
    write_mat(result_path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, result_path + "/parameters.json")


def _list(v: Any) -> list[Any]:
    if v is None or (isinstance(v, str) and v == ""):
        return []
    if isinstance(v, (str, os.PathLike)):
        return [v]
    return list(np.asarray(v).ravel()) if isinstance(v, np.ndarray) else list(v)


@register(
    "XR_tiffToZarr_wrapper",
    Param("dataPaths", "@(x) iscell(x) || ischar(x)", "Dataset folders (one or several)."),
    Param(
        "tiffFullpaths",
        "@(x) iscell(x) || ischar(x)",
        "Tiff files (or a folder of them) when data_paths is empty.",
    ),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("locIds", "@isnumeric", "Location ids of the tiles (stitching).", convert="same"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x).", convert="size_yxz"),
    Param(
        "shardSize", "@isnumeric", "Shard size (z, y, x) (Zarr v3 with shards).", convert="size_yxz"
    ),
    Param(
        "flippedTile", "@(x) isempty(x) || islogical(x)", "Per file: flip it in z.", convert="same"
    ),
    Param(
        "resampleFactor",
        "@(x) isempty(x) || isnumeric(x)",
        "Resampling factors (z, y, x) or (z, xy) (None: none).",
        convert="resample_yxz",
    ),
    Param("partialFile", "@islogical", "Join each file's partial files (*_part[0-9]{3}.tif)."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to convert."),
    Param(
        "inputBbox",
        "@isnumeric",
        "Region (z0, y0, x0, z1, y1, x1) to read (or one per tile location).",
        convert=_BBOX_ROWS,
    ),
    Param(
        "tileOutBbox",
        "@isnumeric",
        "Region (z0, y0, x0, z1, y1, x1) to keep (or one per tile location).",
        convert=_BBOX_ROWS,
    ),
    Param("readWholeTiff", "@islogical", "Read each tiff at once (else by z blocks)."),
    Param(
        "processFunPath",
        "@(x) isempty(x) || isa(x,'function_handle') || ischar(x) || isstring(x) || iscell(x)",
        "User function files (.mat or .txt) per channel.",
    ),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("bigData", "@islogical", "Big dataset (less memory per task)."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param("uuid", "@ischar", "Not used."),
    Param("maxTrialNum", "@isnumeric", "Not used (two tries)."),
    Param("unitWaitTime", "@isnumeric", "Not used."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T01",
)
def tiff_to_zarr_dataset(
    data_paths: str | Sequence[str],
    *,
    tiff_fullpaths: str | Sequence[str] = "",
    result_dir_name: str = "zarr",
    loc_ids: Sequence[float] | None = None,
    block_size: Sequence[int] = (256, 256, 256),
    shard_size: Sequence[int] | None = None,
    flipped_tile: Sequence[bool] | None = None,
    resample_factor: Sequence[float] | float | None = None,
    partial_file: bool = False,
    channel_patterns: Sequence[str] = ("tif",),
    input_bbox: Any = None,
    tile_out_bbox: Any = None,
    read_whole_tiff: bool = True,
    process_fun_path: Any = "",
    parse_cluster: bool = True,
    big_data: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 1,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 3,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Convert the tiff files of datasets (those that match a channel pattern) to zarr
    (``XR_tiffToZarr_wrapper``): ``<data path>/<result_dir_name>/<name>.zarr``, one
    ``tiff_to_zarr`` task per file through the job framework, with per-file flips, per-location
    boxes and per-channel user functions (picked by keywords in the file names).
    parameters.mat and .json of the call go into each result folder."""
    values = {k: v for k, v in locals().items() if k in _T2Z_ARGS}
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.paths import fileparts, matlab_dir, simplify_path
    from petakit5d.core.text import through_mat2str
    from petakit5d.io import get_image_dtype, get_image_size
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_image_filenames

    paths = _list(data_paths)
    paths = [os.fspath(p) for p in paths]
    channel_patterns = [str(c) for c in channel_patterns]
    if paths:
        for d, p in enumerate(paths):
            paths[d] = simplify_path(p)
            result_path = f"{paths[d]}/{result_dir_name}/"
            os.makedirs(result_path, exist_ok=True)
            group_writable(result_path)
            _save_parameters(tiff_to_zarr_dataset, result_path, values)
        _, _, _, tiffs, _ = parse_image_filenames(paths, False, channel_patterns)
        tiffs = list(tiffs)
    else:
        tiffs = tiff_fullpaths
    if isinstance(tiffs, (str, os.PathLike)):
        tiffs = os.fspath(tiffs)
        if os.path.isdir(tiffs):
            names = [e.name for e in matlab_dir(f"{tiffs}/*.tif")]
            tiffs = [f"{tiffs}/{n}" for n in names]
        else:
            tiffs = [tiffs]
    tiffs = [os.fspath(t) for t in tiffs]
    if not tiffs:
        raise ConvertError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )

    data_path = fileparts(tiffs[0])[0]
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(data_path, job_log_dir).parse_cluster)

    # the user functions: channels (rows) x keyword alternatives (columns)
    n_c = len(channel_patterns)
    fun_rows = _fun_rows(process_fun_path)
    n_cols = max(1, max((len(r) for r in fun_rows), default=0))
    usr_strs = [[""] * n_cols for _ in range(max(n_c, len(fun_rows)))]
    keywords = [[[""] for _ in range(n_cols)] for _ in range(max(n_c, len(fun_rows)))]
    if fun_rows:
        _print("Process function paths:")
        for r in fun_rows:
            _print("    " + "    ".join(r))
    for i, r in enumerate(fun_rows):
        for j, fn in enumerate(r):
            if not fn:
                continue
            fun, kw = _read_function_file(fn)
            if fun:
                usr_strs[i][j] = fun
            if kw is not None:
                keywords[i][j] = _keywords(kw)

    n_f = len(tiffs)
    fnames_cell: list[list[str]] = []
    udata_paths: list[str] = []
    if partial_file:
        udata_paths = sorted({fileparts(t)[0] for t in tiffs})
        for d in udata_paths:
            fnames_cell.append([f"{d}/{e.name}" for e in matlab_dir(f"{d}/*tif")])

    compressor = "zstd"

    in_rows, out_rows = _box_rows(input_bbox), _box_rows(tile_out_bbox)
    uniq_loc = sorted(set(_list(loc_ids)))
    n_loc = len(uniq_loc)
    loc_specific = n_loc > 1 and (len(in_rows) == n_loc or len(out_rows) == n_loc)
    if loc_specific:
        if len(in_rows) <= 1:
            in_rows = [in_rows[0] if in_rows else None] * n_loc
        elif len(in_rows) != n_loc:
            raise ConvertError(
                f"The number of InputBbox ({len(in_rows)}) does not match the number of tile "
                f"locations ({n_loc})!"
            )
        if len(out_rows) <= 1:
            out_rows = [out_rows[0] if out_rows else None] * n_loc
        elif len(out_rows) != n_loc:
            raise ConvertError(
                f"The number of tileOutBbox ({len(out_rows)}) does not match the number of "
                f"tile locations ({n_loc})!"
            )
    loc_list = _list(loc_ids)
    flips = _list(flipped_tile)
    rs = None
    if not _empty(resample_factor):
        rs = np.atleast_1d(through_mat2str(np.asarray(resample_factor, np.float64))).tolist()

    zarr_fullpaths: list[str] = []
    tasks = []
    group: list[str] = []
    for i in range(n_f):
        tiff_i = tiffs[i]
        dpath, fsname, _ = fileparts(tiff_i)
        if partial_file:
            did = udata_paths.index(dpath) if len(udata_paths) > 1 else 0
            fnames = fnames_cell[did]
            group = [f for f in fnames if dpath in f and fsname in f]
        else:
            group = [tiff_i]
        zarr_path = f"{dpath}/{result_dir_name}"
        if not os.path.isdir(zarr_path):
            os.makedirs(zarr_path, exist_ok=True)
            group_writable(zarr_path)
        zarr_fullpaths.append(f"{dpath}/{result_dir_name}/{fsname}.zarr/")
        flip = bool(flips[i]) if flips else False

        cind = [p in tiff_i for p in channel_patterns]
        if not any(cind):
            raise ConvertError(
                f"The file {tiff_i} does not match any channel patterns "
                + " ".join(channel_patterns)
            )
        if loc_specific:
            k = uniq_loc.index(loc_list[i])
            in_i, out_i = in_rows[k], out_rows[k]
        else:
            if len(in_rows) > 1 or len(out_rows) > 1:
                # (MATLAB passes the rows on, and tiffToZarr indexes them linearly)
                raise ConvertError(
                    "Boxes per tile location need loc_ids with as many locations as boxes"
                )
            in_i = in_rows[0] if in_rows else None
            out_i = out_rows[0] if out_rows else None

        # the channel's functions, column-major (MATLAB's linear indexing of usrFcn_strs(cind, :))
        rows = [c for c in range(n_c) if cind[c]]
        strs = [usr_strs[r][j] for j in range(n_cols) for r in rows]
        kws = [keywords[r][j] for j in range(n_cols) for r in rows]
        usr_i = strs[0]
        for j in range(len(strs)):
            if any(k in tiff_i for k in kws[j]):
                usr_i = strs[j]
                break

        tasks.append(
            Task.of(
                tiff_to_zarr,
                group,
                zarr_fullpaths[i],
                None,
                block_size=_ints(block_size),
                shard_size=None if _empty(shard_size) else _ints(shard_size),
                flip_zstack=flip,
                resample_factor=rs,
                input_bbox=in_i,
                tile_out_bbox=out_i,
                read_whole_tiff=bool(read_whole_tiff),
                compressor=compressor,
                usr_fcn=usr_i,
            )
        )

    # the memory of a task: of the last file's group (MATLAB's), times the function's factor
    sizes = [get_image_size(f) for f in group]
    byte = data_type_to_byte_number(get_image_dtype(group[0]))
    voxels = float(sizes[0][1]) * float(sizes[0][2]) * float(sum(s[0] for s in sizes))
    if not usr_strs:
        mem_factor = _MEM_FACTOR_PLAIN
    elif any(n in usr_strs[0][0] for n in _FLATFIELD_NAMES):
        mem_factor = _MEM_FACTOR_FLATFIELD
    else:
        # (MATLAB: all(isempty(usrFcn_strs)) is false for any cell of functions, empty or not)
        mem_factor = _MEM_FACTOR_FUNCTION
    mem_allocate = voxels * byte / 1024**3 * mem_factor
    if not big_data:
        mem_allocate *= 2

    common = {
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "max_trial_num": 2,
        "job_log_dir": job_log_dir,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = framework.run_tasks(
        tiffs,
        zarr_fullpaths,
        tasks,
        cpus_per_task=cpus_per_task,
        mem_allocate=mem_allocate,
        **common,
    )
    if not all(done):
        framework.run_tasks(
            tiffs,
            zarr_fullpaths,
            tasks,
            cpus_per_task=cpus_per_task * 2,
            mem_allocate=mem_allocate * 2,
            **common,
        )


_T2Z_ARGS = (
    "data_paths",
    "tiff_fullpaths",
    "result_dir_name",
    "loc_ids",
    "block_size",
    "shard_size",
    "flipped_tile",
    "resample_factor",
    "partial_file",
    "channel_patterns",
    "input_bbox",
    "tile_out_bbox",
    "read_whole_tiff",
    "process_fun_path",
    "parse_cluster",
    "big_data",
    "master_compute",
    "job_log_dir",
    "cpus_per_task",
    "uuid",
    "max_trial_num",
    "unit_wait_time",
    "mcc_mode",
    "config_file",
)


# ---------------------------------------------------------------------------------------------
# zarrToTiff
# ---------------------------------------------------------------------------------------------


@register(
    "zarrToTiff",
    Param("zarrFilename", "@ischar", "Zarr file (either layout)."),
    Param("tiffFilename", "@ischar", "Tiff result."),
    Param("Overwrite", "@islogical", "Replace an existing result."),
    Param(
        "usrFcn",
        "@(x) isempty(x) || isa(x,'function_handle') || ischar(x)",
        "User function applied to the volume (callable, import path or record).",
    ),
    Param("uuid", "@ischar", "Id of the temporary result (empty: a new one)."),
    task="T01",
)
def zarr_to_tiff(
    zarr_filename: str,
    tiff_filename: str,
    *,
    overwrite: bool = False,
    usr_fcn: Any = "",
    uuid: str = "",
) -> None:
    """Convert a zarr array (the port's or MATLAB's layout) to a tiff file (``zarrToTiff``),
    after the user function if one is given. An existing result is kept (``overwrite``:
    replaced first)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive
    from petakit5d.io import read_zarr, write_tiff

    zarr_filename = os.fspath(zarr_filename)
    tiff_filename = os.fspath(tiff_filename)
    if os.path.lexists(tiff_filename):
        if not overwrite:
            _print(f"Tiff result {tiff_filename} already exists, skip it!")
            return
        os.remove(tiff_filename)
    pthstr = fileparts(tiff_filename)[0]
    if pthstr and not os.path.isdir(pthstr):
        mkdir_recursive(pthstr)
    if not uuid:
        uuid = get_uuid()

    _print(f"Convert Zarr to Tiff for {zarr_filename} ...")
    im = read_zarr(zarr_filename)
    if 0 in im.shape:
        _print("The input image is empty, skip it!")
        return
    fcn = user_function(usr_fcn)
    if fcn is not None:
        im = np.asarray(fcn(im))

    tmp = f"{tiff_filename[:-4]}_{uuid}.tif"
    write_tiff(im, tmp)
    del im
    _publish(tmp, tiff_filename, "Tiff")


@register(
    "XR_zarrToTiff_wrapper",
    Param("dataPaths", "@(x) iscell(x) || ischar(x)", "Dataset folders (one or several)."),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the files to convert."),
    Param(
        "usrFcn",
        "@(x) isempty(x) || isa(x,'function_handle') || ischar(x)",
        "User function applied to each volume (import path or record).",
    ),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param("uuid", "@ischar", "Not used."),
    Param("maxTrialNum", "@isnumeric", "Not used (two tries, then one)."),
    Param("unitWaitTime", "@isnumeric", "Not used."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="T01",
)
def zarr_to_tiff_dataset(
    data_paths: str | Sequence[str],
    *,
    result_dir_name: str = "tiffs",
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0"),
    usr_fcn: Any = "",
    parse_cluster: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 1,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 1,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Convert the zarr files of datasets (``<data path>/*.zarr`` that contain a channel
    pattern) to ``<data path>/<result_dir_name>/<name>.tif`` (``XR_zarrToTiff_wrapper``), one
    ``zarr_to_tiff`` task per file through the job framework."""
    from petakit5d.core.paths import matlab_dir
    from petakit5d.io import get_image_size
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import matches_channel_pattern

    paths = [os.fspath(p) for p in _list(data_paths)]
    zarrs: list[str] = []
    tiffs: list[str] = []
    for data_path in paths:
        names = [e.name for e in matlab_dir(f"{data_path}/*.zarr")]
        os.makedirs(f"{data_path}/{result_dir_name}", exist_ok=True)
        zarrs += [f"{data_path}/{n}" for n in names]
        tiffs += [f"{data_path}/{result_dir_name}/{n[:-5]}.tif" for n in names]

    keep = [False] * len(zarrs)
    for p in channel_patterns:
        keep = [a or b for a, b in zip(keep, matches_channel_pattern(zarrs, str(p)), strict=True)]
    zarrs = [z for z, k in zip(zarrs, keep, strict=True) if k]
    tiffs = [t for t, k in zip(tiffs, keep, strict=True) if k]
    if not zarrs:
        return
    fcn = usr_fcn if isinstance(usr_fcn, str) else _function_text(usr_fcn)
    tasks = [Task.of(zarr_to_tiff, z, t, usr_fcn=fcn) for z, t in zip(zarrs, tiffs, strict=True)]

    # each task converts one file
    mem_allocate = max(math.prod(get_image_size(z)[:3]) for z in zarrs) * 4 / 1024**3 * 2.5
    common = {
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "cpus_per_task": cpus_per_task,
        "mem_allocate": mem_allocate,
        "job_log_dir": job_log_dir,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = framework.run_tasks(zarrs, tiffs, tasks, max_trial_num=2, **common)
    if not all(done):
        framework.run_tasks(zarrs, tiffs, tasks, max_trial_num=1, **common)
