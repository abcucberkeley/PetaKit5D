"""Deconvolution of one file or an array (``RLdecon``), with the deskew and rotation after the
decon.

``rl_decon_file`` is RLdecon.m at the oracle commit, step by step:

- **Results and the skip rule.** The result is ``output_fn`` (a tiff, or a zarr with
  ``save_zarr``); with ``deskew`` (and ``dsr_combined`` False) also ``<result folder>/DS/<name>``,
  with ``rotate`` ``DSR/<name>`` (``Rotated/<name>`` for objective scans). ``save_3d_stack``
  ``(decon, DS, DSR)`` says which to write; nothing is computed when every one to write exists
  (or, without the decon, the DS and DSR ones). ``rotate`` without ``deskew`` (the decon
  wrapper's deconRotate) is the combined deskew and rotation for sample scans.
- **The method** is one of 'original', 'simplified', 'omw' and 'cudagen' (case-sensitive), checked
  with the arguments (``MATLAB:InputParser:ArgumentFailedValidation``) before anything is read
  or written.
- **The PSF** (a tiff; a ``.mat`` file's variable ``psf``; another suffix gives an empty PSF, as
  MATLAB's dummy): with ``psf_gen`` the generated PSF ``<result folder>/psfgen/<psf>_<method>.tif``
  of an earlier call is read (``read_tiff_cached``), else ``psf_gen`` (``'masked'``, medFactor
  1.5, the data's z step ``dz sind(angle)``; ``dz`` for objective scans), for methods other than
  OMW cropped to its support (``crop_psf_support``); a PSF larger than the data is cropped to it
  (``crop_psf_to_data``, with MATLAB's warning) and, for file inputs, written as that psfgen
  file if it is missing. Errors in these steps are printed and the PSF kept as far as it got,
  as MATLAB's ``try``/``catch``. OMW's back projector comes from ``omw_backprojector_files``
  (R02's ``decon.backprojector``), cropped to the data likewise. Data in memory without an
  output file has no result folder: the PSF and the back projector are made without their
  files (``generate_backprojector``), and nothing is written.
- **The data**: ``rawdata`` or the file, flipped in z (``flip_zstack``); for methods other than
  'simplified' and 'omw' ``max(single(data) - background, 0)``.
- **The decon**, for at least one iteration and data with a nonzero voxel: the PSF normalized by
  its fixed-order double sum, then 'simplified' (``decon.rl.decon_lucy``) or 'omw'
  (``decon_lucy_omw``) with ``use_gpu`` (R04's GPU path where it runs), the options passed on
  (``edge_erosion`` 0 where an existing eroded mask file is used: the result is then multiplied
  by ``mask > 0`` of that zarr); otherwise the data itself, cropped to ``decon_bbox``. With the
  decon result on disk and DS or DSR to make, the result is read instead.
- **After the decon**: the separate deskew (``deskew_volume`` of the single result, the DS MIP,
  the DS file) and rotation (``rotate_volume``), the combined deskew and rotation
  (``deskew_rotate_volume`` with ``resample``, ``interp_method``, ``x_step_thresh``,
  ``save_16bit``) or the objective-scan rotation; DSR in uint16 (``save_16bit``) or single, its
  MIP, the DSR file.
- **Info files** ('simplified' after a decon): ``decon_info/<name>_info.mat`` (v7.3: err_mat,
  nIter, fixIter, errThresh, debug) and ``decon_info/<name>_actual_iterations_<k>.txt``. Then the
  result in uint16 (``save_16bit``), its MIP ``MIPs/<input name>_MIP_z.tif`` (``mip_axis``) and
  the result file. Every file is written under a temporary name (``<name>_<uuid>``) and renamed
  into place after its MIPs and info files (the publish order). A tiff result replaces one that
  is there, as MATLAB's movefile does (with the same voxels: the decon result read back, or the
  same computation); a zarr result that is there is kept (MATLAB's movefile would move the new
  folder into it).

Returns the decon result, DS (None after the rotation), DSR and the output of ``overlap_fcn``.

Data in memory (no ``input_fn``) with an output file names its results after the output file
(its MIP ``MIPs/<output name>_MIP_z.tif``, its DS and DSR files and their MIPs); without an
output file, the deskew and rotation and ``debug`` raise MATLAB's undefined ``outFsn`` and
``fsname`` before anything is written (there is nowhere to write).

MATLAB's behaviour is kept where it is odd (see the report of task R05): without an output file
the result is returned as the decon gave it; without any 3D stack to save nothing else is
written; ``deskew`` without ``rotate`` under the default ``dsr_combined`` writes nothing; the
separate path ignores ``resample``; ``n_tapering`` (MATLAB's ``edgetaper``) and the 'original'
method (``deconvlucy``, Phase 2) are not ported (``NotImplementedError``); 'cudagen' raises
MATLAB's error for the missing ``decon_lucy_cuda_function``.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from typing import Any, NamedTuple

import numpy as np

from petakit5d.params import Param, register

__all__ = ["RLDeconResult", "rl_decon_file"]

_UNDEFINED = "MATLAB:UndefinedFunction"
_METHODS = ("original", "simplified", "omw", "cudagen")
_METHOD_VALIDATOR = "@(x) any(strcmp(x, {'original', 'simplified', 'omw', 'cudagen'}))"


class RLDeconResult(NamedTuple):
    """RLdecon's outputs: the decon result, the deskewed result (None after the rotation), the
    deskewed and rotated result, and the output of ``overlap_fcn`` (None where it didn't run)."""

    decon: np.ndarray | None
    ds: np.ndarray | None
    dsr: np.ndarray | None
    overlap_out: Any = None


def _undefined(name: str) -> Exception:
    from petakit5d.decon.rl import DeconError

    return DeconError(_UNDEFINED, f"Unrecognized function or variable '{name}'.")


def _defined(value: Any, name: str) -> Any:
    """a variable MATLAB reads: its undefined-variable error where it was never set"""
    if value is None:
        raise _undefined(name)
    return value


def _threads() -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(None)


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _mkdir_group(path: str) -> None:
    """mkdir and, off Windows, ``fileattrib(path, '+w', 'g')``"""
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
        if os.name != "nt":
            from petakit5d.io.tiff import group_writable

            group_writable(path)


def _scalar(v: Any, default: float) -> float:
    a = np.asarray(v if v is not None else default, dtype=np.float64).reshape(-1)
    return float(a[0]) if a.size else default


def _flags3(v: Any) -> tuple[bool, bool, bool]:
    """save3Dstack: one flag or three"""
    a = [bool(x) for x in np.asarray(v).reshape(-1)]
    return (a[0],) * 3 if len(a) == 1 else (a[0], a[1], a[2])


def _read(path: str, threads: int) -> np.ndarray:
    """readtiff or readzarr by the extension; None for another one (MATLAB's data stays [])"""
    from petakit5d import io

    ext = os.path.splitext(path.rstrip("/\\"))[1]
    if ext in (".tif", ".tiff"):
        return io.read_tiff(path, num_threads=threads)
    if ext == ".zarr":
        return io.read_zarr(path, num_threads=threads)
    return None


def _flip_z(a: np.ndarray, own: bool) -> np.ndarray:
    """MATLAB's ``flip(a, 3)``: z planes reversed (an image as it is); in place for an array of
    our own"""
    if a.ndim < 3 or a.shape[0] < 2:
        return a
    if not own or not a.flags.writeable or not a.flags.c_contiguous:
        return np.ascontiguousarray(a[::-1])
    n = a.shape[0]
    tmp = np.empty(a.shape[1:], dtype=a.dtype)
    for k in range(n // 2):
        tmp[...] = a[k]
        a[k] = a[n - 1 - k]
        a[n - 1 - k] = tmp
    return a


def _any_nonzero(a: np.ndarray) -> bool:
    """MATLAB's ``any(a, 'all')`` (NaN counts): the first plane first"""
    if a.size == 0:
        return False
    flat = a.reshape(-1)
    head = flat[: min(flat.size, 1 << 16)]
    if np.any(head != 0):
        return True
    return bool(np.any(a != 0))


def _single(a: np.ndarray, threads: int) -> np.ndarray:
    from petakit5d.kernels import ops

    if a.dtype == np.float32:
        return a
    return ops.convert(a, "single", num_threads=threads)


def _uint16(a: np.ndarray, threads: int) -> np.ndarray:
    from petakit5d.kernels import ops

    if a.dtype == np.uint16:
        return a
    return ops.convert(a, "uint16", num_threads=threads)


def _psf_dims(psf: np.ndarray) -> tuple[int, int, int]:
    """MATLAB's ``size(psf, 1 : 3)`` in (z, y, x)"""
    s = psf.shape
    if psf.ndim == 3:
        return s
    if psf.ndim == 2:
        return (1, *s)
    return (1, 1, max(1, psf.size)) if psf.ndim == 1 else (1, 1, 1)


def _image_size(input_fn: str, rawdata: np.ndarray | None) -> tuple[int, int, int]:
    """getImageSize of the file, or ``size(rawdata, 1 : 3)``, in (z, y, x)"""
    if input_fn:
        from petakit5d.io import get_image_size

        s = tuple(int(v) for v in get_image_size(input_fn))
        return s[:3]
    return _psf_dims(np.asarray(rawdata))


def _normalized(psf: np.ndarray, threads: int) -> np.ndarray:
    """``psf ./ decon_sum_double(psf)`` in the PSF's class (single, or double for a resampled or
    loaded double PSF)"""
    from petakit5d.core import dtypes
    from petakit5d.kernels import decon as K

    s = K.sum_double(psf, threads=threads)
    return np.ascontiguousarray(dtypes.rdivide(psf, s))


def _tmp_tiff(final: str, tag: str) -> str:
    return f"{final[:-4]}_{tag}.tif"


def _move(tmp: str, final: str) -> bool:
    """``movefile(tmp, final)`` of a result: a tiff replaces the file there (MATLAB's movefile);
    a zarr folder is renamed only where none is (MATLAB's movefile would move it into the
    existing folder): False if another one was there (ours is then removed)"""
    from petakit5d.jobs.publish import remove_path, rename_no_replace

    if not os.path.isdir(tmp):
        os.replace(tmp, final)
        return True
    if rename_no_replace(tmp, final):
        return True
    print(f"The result {final} exists, keep it.", flush=True)
    remove_path(tmp)
    return False


def _publish(write: Callable[[str], None], tmp: str, final: str, zarr: bool) -> None:
    """write the result to tmp and move it to final (``_move``); an existing zarr result is kept
    without writing it again"""
    if zarr and _exist(final):
        print(f"The result {final} exists, keep it.", flush=True)
        return
    write(tmp)
    _move(tmp, final)


def _write(data: np.ndarray, path: str, zarr: bool, block_size, threads: int) -> None:
    from petakit5d import io

    if zarr:
        io.write_zarr(data, path, block_size=block_size, num_threads=threads)
    else:
        io.write_tiff(data, path, num_threads=threads)


def _save_mip(im: np.ndarray, mip_fn: str, mip_axis) -> None:
    from petakit5d.tools.mip import save_mip_frame

    save_mip_frame(im, mip_fn, axis=mip_axis)


def _prepare_psf(
    psf_fn: Any,
    *,
    psf_gen: bool,
    input_fn: str,
    psfgen_path: str | None,
    method: str,
    dz_psf: float,
    dz_data: float,
    im_size: tuple[int, int, int],
    loc_uuid: str,
    threads: int,
) -> tuple[Any, str | None]:
    """RLdecon's PSF block: (psf, psf file name without extension); psf None where MATLAB's
    stays undefined (not a file name)"""
    from petakit5d import io
    from petakit5d.core.paths import fileparts
    from petakit5d.decon import psf as P

    if not isinstance(psf_fn, (str, os.PathLike)):
        return None, None
    psf_fn = os.fspath(psf_fn)
    _, psf_fsn, suffix = fileparts(psf_fn)
    if suffix == ".mat":
        from petakit5d.core.arrays import from_matlab
        from petakit5d.io.matfile import read_mat

        v = read_mat(psf_fn, "psf")
        if "psf" not in v:
            raise _undefined("psf")
        return from_matlab(np.asarray(v["psf"])), psf_fsn
    if suffix != ".tif":
        return np.zeros((0, 0), np.float64), psf_fsn  # MATLAB's dummy PSF []

    psfgen_fn = None
    if (psf_gen or input_fn) and psfgen_path:
        psfgen_fn = f"{psfgen_path}/{psf_fsn}_{method}.tif"
    psf = None
    if psf_gen and psfgen_fn is not None and _exist(psfgen_fn):
        print(f"Load existing generated PSF {psfgen_fn} for {psf_fn} ...", flush=True)
        try:
            psf = _single(np.asarray(P.read_tiff_cached(psfgen_fn)), threads)
            psf_gen = False
        except Exception as e:  # noqa: BLE001 (MATLAB's catch: disp(ME), read the raw PSF)
            print(e, flush=True)
    if psf is None:
        print(f"Load PSF {psf_fn} ...", flush=True)
        psf = _single(io.read_tiff(psf_fn, num_threads=threads), threads)
    try:
        if psf_gen:
            print(f"PSF generation for {psf_fn} ...", flush=True)
            psf = P.generate_psf(psf, dz_psf, dz_data, 1.5, "masked", threads=threads)
            if method != "omw":
                psf = P.crop_psf_support(psf)
        if psf.ndim in (2, 3):
            psf = P.crop_psf_to_data(psf, im_size)
        if input_fn and not _exist(psfgen_fn):
            from petakit5d.jobs.publish import remove_path, rename_no_replace

            tmp = f"{psfgen_path}/{psf_fsn}_{loc_uuid}.tif"
            io.write_tiff(psf, tmp, num_threads=threads)
            if not rename_no_replace(tmp, psfgen_fn):
                remove_path(tmp)
    except Exception as e:  # noqa: BLE001 (MATLAB's catch: disp(ME), the PSF as far as it got)
        print(e, flush=True)
    return psf, psf_fsn


def _backprojector(psf, psfgen_path, psf_fsn, psf_fn, im_size, kw) -> np.ndarray:
    """the OMW back projector, from (or into) the psfgen folder; made without its files where
    there is no folder (data in memory without an output file)"""
    from petakit5d.decon import backprojector as B
    from petakit5d.decon import psf as P

    if psf is None:
        raise _undefined("psf")
    if not psfgen_path:
        b = B.generate_backprojector(
            psf,
            kw["wiener_alpha"],
            kw["skewed"],
            kw["otf_cum_thresh"],
            hann_win_bounds=kw["hann_win_bounds"],
        ).b_omw
    else:
        b = B.omw_backprojector_files(
            psf,
            psfgen_path,
            _defined(psf_fsn, "psfFsn"),
            kw["wiener_alpha"],
            kw["otf_cum_thresh"],
            kw["skewed"],
            kw["hann_win_bounds"],
            os.fspath(psf_fn),
        ).b_omw
    return P.crop_psf_to_data(np.asarray(b), im_size)


def _decon_call(rawdata, psf, psf_b, n_iter, method, kw, threads):
    """the switch of RLdecon's methods: (result, err_mat or None, iterations, overlap output)"""
    from petakit5d.decon import rl

    common = {
        "background": kw["background"],
        "use_gpu": kw["use_gpu"],
        "save_16bit": kw["save_16bit"],
        "damp_factor": kw["damp_factor"],
        "scale_factor": kw["scale_factor"],
        "decon_offset": kw["decon_offset"],
        "bbox": kw["decon_bbox"],
        "edge_erosion": kw["edge_erosion"],
        "debug": kw["debug"],
        "debug_folder": kw["debug_folder"],
        "save_step": kw["save_step"],
        "overlap_fcn": kw["overlap_fcn"],
    }
    if method == "simplified":
        r = rl.decon_lucy(rawdata, psf, n_iter, **common)
        return r.decon, r.err_mat, r.iterations, r.overlap_out
    if method == "omw":
        r = rl.decon_lucy_omw(rawdata, psf, psf_b, n_iter, **common)
        return r.decon, r.err_mat, None, r.overlap_out
    if method == "original":
        raise NotImplementedError(
            "RLMethod 'original' (MATLAB's deconvlucy) is not ported yet (task B11, Phase 2); "
            "use 'simplified' or 'omw'"
        )
    # 'cudagen', the only other method the arguments let through
    from petakit5d.decon.rl import DeconError

    raise DeconError(
        _UNDEFINED,
        "Undefined function 'decon_lucy_cuda_function' for input arguments of type 'single'.",
    )


def _scale(dec: np.ndarray, scale: float) -> np.ndarray:
    from petakit5d.core import dtypes

    return np.asarray(dtypes.times(dec, scale))


def _apply_mask_file(dec: np.ndarray, mask_fn: str, threads: int) -> np.ndarray:
    """``deconvolved .* cast(mask > 0, class(deconvolved))`` of the eroded-mask zarr"""
    from petakit5d import io

    m = io.read_zarr(mask_fn, num_threads=threads)
    if m.shape != dec.shape and m.squeeze().shape != dec.squeeze().shape:
        from petakit5d.decon.rl import DeconError

        raise DeconError(
            "MATLAB:sizeDimensionsMustMatch", "Arrays have incompatible sizes for this operation."
        )
    keep = (m > 0).reshape(dec.shape)
    out = dec if dec.flags.writeable else dec.copy()
    np.multiply(out, keep.astype(out.dtype), out=out)
    return out


def _subtract_background(raw: np.ndarray, background: Any, threads: int) -> np.ndarray:
    """``max(single(rawdata) - Background, 0)`` (Background [] errors as MATLAB's empty operand)"""
    from petakit5d.decon import rl

    b = np.asarray(background) if background is not None else np.zeros((0, 0))
    if b.size == 0:
        from petakit5d.decon.rl import DeconError

        raise DeconError(
            "MATLAB:sizeDimensionsMustMatch", "Arrays have incompatible sizes for this operation."
        )
    x = _single(raw, threads)
    scalar, bg = rl._background(background)
    return rl._subtract_background(x, scalar, bg, threads, inplace=x is not raw)


def _crop(a: np.ndarray, box) -> np.ndarray:
    from petakit5d.decon.rl import _crop as crop

    return crop(a, box)


def _info_files(info_path: str, out_fsn: str, err_mat, n_iter, fix_iter, err_thresh, debug, k):
    from petakit5d.io.matfile import save_ascii, write_mat

    _mkdir(info_path)
    info_fn = f"{info_path}/{out_fsn}_info.mat"
    write_mat(
        info_fn,
        {
            "err_mat": np.asarray(err_mat, np.float64),
            "nIter": float(n_iter),
            "fixIter": bool(fix_iter),
            "errThresh": np.asarray(
                np.zeros((0, 0)) if err_thresh is None else err_thresh, np.float64
            ),
            "debug": bool(debug),
        },
    )
    save_ascii(f"{info_fn[:-8]}actual_iterations_{int(k)}.txt", float(k))


def _volume_out(v):
    from petakit5d.core.arrays import to_matlab

    return np.zeros((0, 0)) if v is None else to_matlab(v)


def _matlab_out(r: RLDeconResult):
    return (_volume_out(r.decon), _volume_out(r.ds), _volume_out(r.dsr), r.overlap_out)


@register(
    "RLdecon",
    Param("inputFn", None, "Input file (tiff or zarr); '' for data in rawdata."),
    Param("outputFn", None, "Result file; '' to only return the result."),
    Param("PSFfn", None, "PSF file (tiff, or a .mat file with the variable psf).", py="psf_fn"),
    Param("xyPixelSize", "@isnumeric", "Pixel size (um)."),
    Param("dz", "@isnumeric", "Scan step (um)."),
    Param("dzPSF", "@isnumeric", "z step of the PSF (um)."),
    Param(
        "rawdata", "@isnumeric", "Data (z, y, x) used instead of reading inputFn.", convert="volume"
    ),
    Param("nTapering", "@isnumeric", "Edge tapering in x-y (edgetaper; not ported)."),
    Param("save16bit", "@islogical", "Results in uint16 (else single).", convert="same"),
    Param("objectiveScan", "@islogical", "Objective scan data (rotation without deskew)."),
    Param("Deskew", "@islogical", "Deskew the result."),
    Param("Rotate", "@islogical", "Rotate the result (with the deskew for sample scans)."),
    Param("DSRCombined", "@islogical", "Deskew and rotate in one step."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("SkewAngle", "@isnumeric", "Skew angle (degrees)."),
    Param("flipZstack", "@islogical", "Flip the data in z first."),
    Param("interpMethod", "@(x) any(strcmpi(x, {'cubic', 'linear'}))", "Interpolation."),
    Param("Crop", "@islogical", "Crop of the separate deskew."),
    Param("xStepThresh", "@isnumeric", "x step above which the combined path interpolates."),
    Param(
        "Background", "@isnumeric", "Background ([]: 0 for 'simplified' and 'omw').", convert="same"
    ),
    Param(
        "Resample",
        "@isnumeric",
        "Resample factor of the combined deskew and rotation (z, y, x).",
        convert="resample_xyz",
    ),
    Param("DeconIter", "@isnumeric", "Number of iterations."),
    Param("RLMethod", _METHOD_VALIDATOR, "'simplified', 'omw', 'original' or 'cudagen'."),
    Param("wienerAlpha", "@isnumeric", "Wiener alpha of the OMW back projector."),
    Param("OTFCumThresh", "@isnumeric", "OTF cumulative sum threshold (OMW)."),
    Param("hannWinBounds", "@isnumeric", "Apodization range of the distance matrix (OMW)."),
    Param("skewed", "@(x) isempty(x) || islogical(x)", "Decon in skewed space (None: automatic)."),
    Param("fixIter", "@islogical", "Fixed number of iterations (saved in the info file)."),
    Param("errThresh", "@isnumeric", "Error threshold (saved in the info file)."),
    Param("saveZarr", "@islogical", "Save results as zarr."),
    Param("blockSize", "@isnumeric", "Block size (z, y, x) of zarr results.", convert="size_yxz"),
    Param("dampFactor", "@isnumeric", "Damp factor of the result."),
    Param("scaleFactor", "@isnumeric", "Scale factor of the result ([]: 1)."),
    Param("deconOffset", "@isnumeric", "Offset of the result.", convert="same"),
    Param("EdgeErosion", "@isnumeric", "Erosion of the edges of the mask of nonzero data."),
    Param("ErodeMaskfile", "@ischar", "Eroded mask zarr used instead of the erosion."),
    Param("deconBbox", "@isnumeric", "Box (z0, y0, x0, z1, y1, x1) of the result.", convert="bbox"),
    Param("debug", "@islogical", "Save the iterations' estimates."),
    Param("saveStep", "@isnumeric", "Save an estimate every saveStep iterations (debug)."),
    Param(
        "useGPU",
        "@(x) isempty(x) || islogical(x)",
        "True: on a GPU (an error if none is usable); False: the CPU; None: a GPU if one is "
        "usable.",
    ),
    Param(
        "save3Dstack",
        "@(x) islogical(x) && (numel(x) == 1 || numel(x) == 3)",
        "Save the (decon, DS, DSR) stacks.",
        convert="same",
        py="save_3d_stack",
    ),
    Param(
        "mipAxis", "@(x) isnumeric(x) && numel(x) == 3", "MIP axes (z, y, x).", convert="flags_yxz"
    ),
    Param("psfGen", "@islogical", "Generate the PSF (psf_gen) or use it as it is."),
    Param("uuid", "@ischar", "uuid of the temporary names."),
    Param(
        "overlapFcn",
        "@(x) isempty(x) || isa(x, 'function_handle')",
        "Function run on the CPU while the GPU deconvolves; its output is returned.",
    ),
    task="R05",
    matlab_out=_matlab_out,
)
def rl_decon_file(
    input_fn,
    output_fn,
    psf_fn,
    xy_pixel_size,
    dz,
    dz_psf,
    *,
    rawdata=None,
    n_tapering=0,
    save_16bit=True,
    objective_scan=False,
    deskew=False,
    rotate=False,
    dsr_combined=True,
    reverse=False,
    skew_angle=32.45,
    flip_zstack=False,
    interp_method="linear",
    crop=False,
    x_step_thresh=2,
    background=None,
    resample=None,
    decon_iter=15,
    rl_method="omw",
    wiener_alpha=0.005,
    otf_cum_thresh=0.9,
    hann_win_bounds=(0.8, 1),
    skewed=None,
    fix_iter=True,
    err_thresh=None,
    save_zarr=False,
    block_size=(256, 256, 256),
    damp_factor=1,
    scale_factor=None,
    decon_offset=0,
    edge_erosion=0,
    erode_maskfile="",
    decon_bbox=None,
    debug=False,
    save_step=5,
    use_gpu=None,
    save_3d_stack=(True, True, True),
    mip_axis=(1, 0, 0),
    psf_gen=True,
    uuid="",
    overlap_fcn=None,
) -> RLDeconResult:
    """Deconvolve a file (or ``rawdata``) and deskew and rotate the result (RLdecon): see the
    module. Arrays are (z, y, x); boxes, sizes and MIP axes in (z, y, x)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.matlabcompat.elementary import sind

    t = _threads()
    input_fn = os.fspath(input_fn) if input_fn else ""
    output_fn = os.fspath(output_fn) if output_fn else ""
    # (inputParser's check of RLMethod: before anything is read or written)
    if not (isinstance(rl_method, str) and rl_method in _METHODS):
        from petakit5d.params.spec import MatlabArgumentError

        raise MatlabArgumentError(
            "The value of 'RLMethod' is invalid. It must satisfy the function: "
            "@(x)any(strcmp(x,{'original','simplified','omw','cudagen'})).",
            "MATLAB:InputParser:ArgumentFailedValidation",
        )
    method = str(rl_method)
    if n_tapering is not None and _scalar(n_tapering, 0) > 0:
        raise NotImplementedError(
            "nTapering (MATLAB's fspecial and edgetaper) is not ported; no caller of RLdecon "
            "sets it"
        )

    angle = abs(float(skew_angle))
    dz = float(dz)
    xy = float(xy_pixel_size)
    dz_ratio = float(sind(angle)) if not objective_scan else 1.0
    z_aniso = dz_ratio * dz / xy
    if not uuid:
        uuid = get_uuid()
    loc_uuid = get_uuid()

    fsname = None
    if input_fn:
        data_path, fsname, _ = fileparts(input_fn)
    out_fsn = None
    if output_fn:
        decon_path, out_fsn, _ = fileparts(output_fn)
        # in memory (rawdata), the results are named after the output file
        if not input_fn:
            fsname = out_fsn
    elif input_fn:
        decon_path = f"{data_path}/matlab_decon/"
    else:
        # in memory without an output file: no folder (nothing is saved)
        decon_path = ""
    ext = ".zarr" if save_zarr else ".tif"

    save3 = list(_flags3(save_3d_stack))
    exists = [False, False, False]
    deskew, rotate, dsr_combined = bool(deskew), bool(rotate), bool(dsr_combined)
    # rotation without deskew (the decon wrapper's deconRotate) is the combined deskew and
    # rotation for sample scans
    if not deskew and rotate and not objective_scan:
        deskew = dsr_combined = True
    ds_path = ds_fn = dsr_path = dsr_fn = None
    if deskew:
        if not dsr_combined:
            ds_path = f"{decon_path}/DS/"
            ds_fn = f"{ds_path}/{_defined(out_fsn, 'outFsn')}{ext}"
        if rotate:
            dsr_path = f"{decon_path}/DSR/"
            dsr_fn = f"{dsr_path}/{_defined(out_fsn, 'outFsn')}{ext}"
    if objective_scan:
        deskew = False
        if rotate:
            dsr_path = f"{decon_path}/Rotated/"
            dsr_fn = f"{dsr_path}/{_defined(out_fsn, 'outFsn')}{ext}"
    if not deskew or dsr_combined:
        save3[1] = False
        exists[1] = True
    if not rotate:
        save3[2] = False
        exists[2] = True
    if save3[0] and _exist(output_fn):
        exists[0] = True
    if save3[1] and _exist(ds_fn):
        exists[1] = True
    if save3[2] and _exist(dsr_fn):
        exists[2] = True
    if all(exists) or (not save3[0] and (deskew or rotate) and all(exists[1:])):
        print(
            f"The output file(s) {output_fn} and related files already exist, skip the processing!",
            flush=True,
        )
        return RLDeconResult(None, None, None, None)

    if input_fn:
        im_size = _image_size(input_fn, None)
    else:
        if rawdata is None:
            raise ValueError("RLdecon needs an input file or rawdata")
        im_size = _image_size("", rawdata)

    # (without a folder, the PSF and the back projector are made without their files)
    psfgen_path = ""
    if (psf_gen or input_fn or method == "omw") and decon_path:
        psfgen_path = f"{decon_path}/psfgen/"
        _mkdir(psfgen_path)

    psf, psf_fsn = _prepare_psf(
        psf_fn,
        psf_gen=bool(psf_gen),
        input_fn=input_fn,
        psfgen_path=psfgen_path,
        method=method,
        dz_psf=float(dz_psf),
        dz_data=dz * dz_ratio,
        im_size=im_size,
        loc_uuid=loc_uuid,
        threads=t,
    )
    opts = {
        "wiener_alpha": float(wiener_alpha),
        "otf_cum_thresh": float(otf_cum_thresh),
        "skewed": skewed,
        "hann_win_bounds": tuple(float(v) for v in np.asarray(hann_win_bounds).reshape(-1)),
    }
    psf_b = None
    if method == "omw":
        psf_b = _backprojector(psf, psfgen_path, psf_fsn, psf_fn, im_size, opts)

    own = False
    if rawdata is None:
        raw = _read(input_fn, t)
        if raw is None:
            raise ValueError(f"RLdecon reads tiff and zarr files, not {input_fn}")
        own = True
    else:
        raw = np.asarray(rawdata)
    if flip_zstack:
        raw = _flip_z(raw, own)
    if method not in ("omw", "simplified"):
        raw = _subtract_background(raw, background, t)

    n_iter = _scalar(decon_iter, 15)
    use_mask_file = (
        _scalar(edge_erosion, 0) > 0 and bool(erode_maskfile) and _exist(str(erode_maskfile))
    )
    decon = ds = dsr = None
    err_mat = None
    iterations = None
    overlap_out = None
    if not exists[0]:
        if n_iter > 0 and _any_nonzero(raw):
            print(f"Deconvolution for frame {input_fn} with {method} Method ...", flush=True)
            sf = 1.0 if scale_factor is None or np.size(scale_factor) == 0 else scale_factor
            if psf is None:
                raise _undefined("psf")
            psf = _normalized(psf, t)
            n_save = save_step
            if debug:
                if fsname is None:
                    raise _undefined("fsname")
                debug_folder = f"{decon_path}/{fsname}_debug/"
                _mkdir(debug_folder)
                n_save = min(_scalar(save_step, 5), n_iter)
            else:
                debug_folder = "/tmp/"
            kw = {
                "background": background,
                "use_gpu": use_gpu,
                "save_16bit": bool(save_16bit),
                "damp_factor": damp_factor,
                "scale_factor": sf,
                "decon_offset": decon_offset,
                "decon_bbox": decon_bbox,
                "edge_erosion": 0 if use_mask_file else edge_erosion,
                "debug": bool(debug),
                "debug_folder": debug_folder,
                "save_step": n_save,
                "overlap_fcn": overlap_fcn,
            }
            decon, err_mat, iterations, overlap_out = _decon_call(
                raw, psf, psf_b, decon_iter, method, kw, t
            )
            raw = None
            if method not in ("simplified", "omw") and _scalar(sf, 1) != 1:
                decon = _scale(decon, sf)
            if use_mask_file:
                decon = _apply_mask_file(decon, str(erode_maskfile), t)
        else:
            decon = raw
            if decon_bbox is not None and len(decon_bbox) > 0:
                decon = _crop(decon, decon_bbox)
    elif deskew or rotate:
        decon = _read(output_fn, t)
    raw = None

    if deskew and (not rotate or not dsr_combined):
        from petakit5d.deskew_rotate import deskew_volume, rotate_volume

        if not exists[1]:
            if fsname is None:
                raise _undefined("fsname")
            print(f"Deskew deconvolved frame {fsname}...", flush=True)
            _mkdir_group(ds_path)
            ds = deskew_volume(
                _single(decon, t),
                angle,
                dz,
                xy,
                bool(reverse),
                crop=bool(crop),
                interp=interp_method,
            )
            _save_mip(ds, f"{ds_path}/MIPs/{fsname}_MIP_z.tif", mip_axis)
            if save3[1]:
                tmp = f"{ds_path}{fsname}_{uuid}{ext}"
                data = _uint16(ds, t) if save_16bit else ds
                _publish(lambda p: _write(data, p, save_zarr, block_size, t), tmp, ds_fn, save_zarr)
        elif rotate and not exists[2]:
            ds = _single(_read(ds_fn, t), t)
        if rotate:
            if fsname is None:
                raise _undefined("fsname")
            print(f"Rotate deskewed deconvolved frame {fsname}...", flush=True)
            dsr = rotate_volume(
                ds,
                angle,
                z_aniso,
                bool(reverse),
                crop=True,
                objective_scan=bool(objective_scan),
                interp=interp_method,
            )
            ds = None
    elif deskew and rotate and dsr_combined:
        from petakit5d.deskew_rotate import deskew_rotate_volume

        if fsname is None:
            raise _undefined("fsname")
        print(f"Deskew, rotate and resample for deconvolved frame {fsname}...", flush=True)
        dsr = deskew_rotate_volume(
            decon,
            angle,
            dz,
            xy,
            reverse=bool(reverse),
            crop=True,
            objective_scan=bool(objective_scan),
            resample_factor=resample,
            interp_method=interp_method,
            x_step_thresh=float(x_step_thresh),
            save_16bit=bool(save_16bit),
        )
    if objective_scan and rotate:
        from petakit5d.deskew_rotate import rotate_volume

        if fsname is None:
            raise _undefined("fsname")
        print(f"Rotate deconvolved frame {fsname}...", flush=True)
        dsr = rotate_volume(
            _single(decon, t),
            angle,
            z_aniso,
            bool(reverse),
            crop=True,
            objective_scan=True,
            interp=interp_method,
        )

    if rotate:
        dsr = _uint16(dsr, t) if save_16bit else _single(dsr, t)
        _save_mip(dsr, f"{dsr_path}/MIPs/{fsname}_MIP_z.tif", mip_axis)
        if save3[2]:
            tmp = f"{dsr_path}{fsname}_{uuid}{ext}"
            _publish(lambda p: _write(dsr, p, save_zarr, block_size, t), tmp, dsr_fn, save_zarr)

    if not output_fn or not any(save3):
        return RLDeconResult(decon, ds, dsr, overlap_out)

    if n_iter > 0 and method == "simplified" and err_mat is not None:
        _info_files(
            f"{decon_path}/decon_info/",
            out_fsn,
            err_mat,
            n_iter,
            fix_iter,
            err_thresh,
            debug,
            iterations,
        )
    if save_16bit:
        decon = _uint16(decon, t)
    if any(bool(v) for v in np.asarray(mip_axis).reshape(-1)):
        if fsname is None:
            raise _undefined("fsname")
        _save_mip(decon, f"{decon_path}/MIPs/{fsname}_MIP_z.tif", mip_axis)
    if save3[0]:
        tmp = f"{output_fn[:-4]}_{uuid}{ext}"
        _publish(lambda p: _write(decon, p, save_zarr, block_size, t), tmp, output_fn, save_zarr)
    return RLDeconResult(decon, ds, dsr, overlap_out)
