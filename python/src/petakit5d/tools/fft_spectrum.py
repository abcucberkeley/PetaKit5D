"""FFT power spectra (MATLAB's ``microscopeDataProcessing/fft_spectrum``): one file
(``XR_fftSpectrumComputingFrame``) and datasets (``XR_fftSpectrumComputingWrapper``).

A file's spectrum: the image (its central region, read from the tiff or zarr file) resampled to
isotropic voxels of ``out_pixel_size`` (``matlabcompat.resize.imresize3_fast``, only the region
of the resized image that is kept), padded with zeros or cropped to ``out_size`` around its
center, and ``|fftn|`` of it (``matlabcompat.fft``: MATLAB's driver on the ``petakit5d.fft``
backend; ``|.|`` is MATLAB's abs of complex single, msun's ``hypotf``). Saved next to the
result name:

- ``MIPs/<name>_MIP_{y,x,z}.tif``: the maxima of the spectrum along y, x and z;
- ``central_slices/<name>_center_{xy,yz,xz}.tif``: its planes through DC;
- both normalized to the spectrum's maximum, with gamma 0.5 (``sqrt``), fftshifted, single;
- the result ``<name>.tif``: the log10 of the normalized, fftshifted spectrum (single) with
  ``save_3d_stack``, else an empty file (published last, by one rename).

The projections are 2-D images in MATLAB's orientation (``squeeze``): the y MIP and the xz
slice are ``(nx, nz)``, the x MIP and the yz slice ``(ny, nz)``, the z MIP and the xy slice
``(ny, nx)``. Arrays are ``(z, y, x)`` C order; ``out_size`` is ``(z, y, x)`` (one number: every
axis).

Arithmetic is MATLAB's: the region read and the resized size in double (``ceil``, ``round``
half away from zero, the read box's non-integer ends truncated as the readers' mex do), the
background subtracted in double and rounded to single, the normalization in single.

MATLAB quirks kept (each has a test or an oracle case):

- an even ``out_size`` (or one that is not three numbers) errors when the resized image must
  be cropped or padded to it, and works when the resized image has that size;
- a zarr file's spectrum is a tiff file named ``<name>.zarr`` (the dataset wrapper keeps the
  input's file name);
- an image with one z plane errors (MATLAB's size of a 2-D array has two values);
- ``overwrite`` and ``debug`` of the dataset wrapper are not used: an existing result is kept;
- the wrapper passes the pixel sizes to its tasks as ``%0.20d`` text, ``background`` as ``%d``
  (``%e``, 7 digits, for a non-integer) and ``out_size`` through ``mat2str``.
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "FFTSpectrumError",
    "fft_spectrum_dataset",
    "fft_spectrum_file",
    "spectrum_region",
]

_DIMS = "MATLAB:sizeDimensionsMustMatch"
_DIMS_MSG = "Arrays have incompatible sizes for this operation."
_CHUNK = 1 << 22


class FFTSpectrumError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _floats(v: Any) -> np.ndarray:
    return np.asarray(v, dtype=np.float64).ravel()


def _yxz(v: Any) -> np.ndarray:
    """a (z, y, x) vector (or one number) as MATLAB's [y x z] (one number stays one)"""
    a = _floats(v)
    return a[[1, 2, 0]] if a.size == 3 else a


def _mround(v: np.ndarray) -> np.ndarray:
    from petakit5d.core.dtypes import mround

    return np.asarray(mround(np.asarray(v, dtype=np.float64)), dtype=np.float64)


def _pixel_sizes(xy_pixel_size: float, dz: float, out_pixel_size: Any) -> tuple[np.ndarray, Any]:
    """MATLAB's [xy, xy, dz] and the output pixel size (min(xy, dz) without one)"""
    xyz = np.array([float(_floats(xy_pixel_size)[0])] * 2 + [float(_floats(dz)[0])])
    if out_pixel_size is None or np.asarray(out_pixel_size).size == 0:
        px: Any = min(float(xyz[0]), float(xyz[2]))
    else:
        px = _floats(out_pixel_size)
    return xyz, px


def spectrum_region(
    im_size: Sequence[int],
    xy_pixel_size: float,
    dz: float,
    out_pixel_size: float | None,
    out_size: Sequence[float] | float,
) -> tuple[int, int, int, int, int, int]:
    """The box ``(z0, y0, x0, z1, y1, x1)`` of a ``(z, y, x)`` image that
    ``fft_spectrum_file`` reads: the center that, resized, covers ``out_size`` voxels, with 10%
    more (MATLAB's ``ceil(N ./ (pixel sizes ./ px)) .* 1.1``, its non-integer ends truncated)."""
    xyz, px = _pixel_sizes(xy_pixel_size, dz, out_pixel_size)
    sz = _yxz(im_size)
    n = _yxz(out_size)
    keep = np.minimum(sz, np.ceil(n / (xyz / px)) * 1.1)
    s = np.maximum(1.0, np.ceil((sz - keep) / 2))
    t = np.minimum(sz, s + keep - 1)
    # (the readers' mex take the box as uint64: the ends truncated)
    y0, x0, z0 = (int(v) - 1 for v in s)
    y1, x1, z1 = (int(math.floor(v)) for v in t)
    return z0, y0, x0, z1, y1, x1


def _read_region(fn: str, zarr_file: bool, box: tuple[int, ...]) -> np.ndarray:
    from petakit5d.io import read_tiff, read_zarr

    z0, y0, x0, z1, y1, x1 = box
    if zarr_file:
        return read_zarr(fn, input_bbox=box)
    im = read_tiff(fn, range=(z0, z1))
    if im.ndim == 2:
        return im[y0:y1, x0:x1]
    return im[:, y0:y1, x0:x1]


def _abs(c: np.ndarray, nt: int) -> np.ndarray:
    """MATLAB's abs of complex single (msun's hypotf; the decon kernels'), of real data the
    magnitude"""
    if not np.iscomplexobj(c):
        return np.abs(c)
    from petakit5d._ext import load

    flat = np.ascontiguousarray(c, dtype=np.complex64).reshape(-1).view(np.float32)
    return load("decon").omw_abs(flat, nt).reshape(c.shape)


def _max_all(a: np.ndarray) -> np.float32:
    """``max(a(:))``: NaN ignored, NaN if all are"""
    with np.errstate(invalid="ignore"):
        m = np.max(a)
    if not np.isnan(m):
        return m
    v = a[~np.isnan(a)]
    return v.max() if v.size else m


def _par(fn, n: int, nt: int) -> None:
    """``fn(start, stop)`` on pieces of ``range(n)``, on up to ``nt`` threads"""
    pieces = max(1, -(-n // _CHUNK))
    edges = [n * i // pieces for i in range(pieces + 1)]
    if nt <= 1 or pieces == 1:
        for i in range(pieces):
            fn(edges[i], edges[i + 1])
        return
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=min(nt, pieces)) as ex:
        for f in [ex.submit(fn, edges[i], edges[i + 1]) for i in range(pieces)]:
            f.result()


def _log10_single(q: np.ndarray, nt: int, divisor: np.float32 | None = None) -> np.ndarray:
    """``log10(q)`` (``log10(q ./ divisor)``, the quotient in single) of single data,
    correctly rounded: the double log10 rounded to single, and where that is within a few
    double ulps of a midpoint of singles, the log10 decided in 40 digits, so the bits don't
    depend on the C library. (MATLAB's single log10 is 1 or 2 ulp off it for some values.)"""
    out = np.empty(q.shape, np.float32)
    src = np.ascontiguousarray(q).reshape(-1)
    dst = out.reshape(-1)

    def run(i0: int, i1: int) -> None:
        x = src[i0:i1]
        with np.errstate(divide="ignore", invalid="ignore"):
            if divisor is not None:
                x = x / divisor
            d = np.log10(x.astype(np.float64))
            r = d.astype(np.float32)
            fin = np.isfinite(d) & (d != 0)
            if fin.any():
                # the midpoint of singles nearest to d, exact in double
                r64 = r.astype(np.float64)
                up = np.nextafter(r, np.float32(np.inf))
                nb = np.where(d > r64, up, np.nextafter(r, np.float32(-np.inf)))
                mid = (r64 + nb.astype(np.float64)) / 2
                near = fin & (np.abs(d - mid) <= 8 * np.spacing(np.abs(d)))
                for i in np.flatnonzero(near):
                    r[i] = _log10_decided(float(x[i]), float(mid[i]), r[i], nb[i])
        dst[i0:i1] = r

    _par(run, src.size, nt)
    return out


def _log10_decided(x: float, mid: float, r: np.float32, nb: np.float32) -> np.float32:
    """the single nearest to log10(x) of the two around the midpoint ``mid``"""
    from decimal import Context, Decimal

    v = Context(prec=40).log10(Decimal(x))
    m = Decimal(mid)
    if v == m:
        # (a tie: the even significand)
        return r if (int(np.float32(r).view(np.uint32)) & 1) == 0 else nb
    above = v > m
    return max(r, nb) if above else min(r, nb)


def _image(a: np.ndarray, axis: int) -> np.ndarray:
    """MATLAB's ``squeeze`` of a 3-D array with a singleton along ``axis``, as the C-ordered
    (rows, columns) image writetiff writes: z: ``(ny, nx)``, y: ``(nx, nz)``, x: ``(ny, nz)``"""
    if axis == -3:
        return np.ascontiguousarray(a[0])
    if axis == -2:
        return np.ascontiguousarray(a[:, 0, :].T)
    return np.ascontiguousarray(a[:, :, 0].T)


def _normalize_gamma(m: np.ndarray, im_max: np.float32) -> np.ndarray:
    """``fftshift(single((m ./ imMax) .^ 0.5))`` (single arithmetic)"""
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.sqrt(m / im_max)
    return np.fft.fftshift(g)


def _resized(
    im: np.ndarray,
    xyz: np.ndarray,
    px: Any,
    n: np.ndarray,
    interp_method: str,
    nt: int,
) -> np.ndarray:
    """the image resized to ``px`` voxels, padded and cropped to ``n`` (MATLAB's [y x z])"""
    from petakit5d.core.text import mat2str
    from petakit5d.matlabcompat.resize import imresize3_fast

    msz = np.array([im.shape[1], im.shape[2], im.shape[0]], dtype=np.float64)
    rsz = _mround(msz * xyz / px)
    if not np.any(rsz != n):
        return imresize3_fast(im, rsz[[2, 0, 1]], interp_method, num_threads=nt)
    if n.size != 3 or np.any(np.mod(n, 2) == 0):
        raise FFTSpectrumError(
            f"The output size {mat2str(n)} must be three odd numbers to crop or pad the "
            f"resized image (size {mat2str(rsz)}) to it."
        )
    hfn = (n - 1) / 2
    c = _mround((rsz + 1) / 2)
    s = np.maximum(1.0, c - hfn)
    t = np.minimum(rsz, c + hfn)
    box = [int(s[2]) - 1, int(s[0]) - 1, int(s[1]) - 1, int(t[2]), int(t[0]), int(t[1])]
    r = imresize3_fast(im, rsz[[2, 0, 1]], interp_method, box, num_threads=nt)
    if np.any(rsz < n):
        pre = np.maximum(0, np.floor((n - rsz) / 2)).astype(np.int64)
        post = np.maximum(0, np.ceil((n - rsz) / 2)).astype(np.int64)
        r = np.pad(r, [(pre[2], post[2]), (pre[0], post[0]), (pre[1], post[1])])
    return r


# ---------------------------------------------------------------------------------------------
# XR_fftSpectrumComputingFrame
# ---------------------------------------------------------------------------------------------


@register(
    "XR_fftSpectrumComputingFrame",
    Param("frameFullpath", "@ischar", "Input tiff or zarr file."),
    Param("spectrumFullname", "@ischar", "Result file (.tif) in the spectrum folder."),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Actual pixel size in z (um)."),
    Param("zarrFile", "@islogical", "The input is zarr."),
    Param(
        "outPixelSize",
        "@(x) isnumeric(x) || isempty(x)",
        "Output voxel size (um, every axis; None: the smaller of xyPixelSize and dz).",
    ),
    Param(
        "outSize",
        "@isnumeric",
        "Output size (z, y, x): odd numbers to crop or pad to.",
        convert="size_yxz",
    ),
    Param("save3DStack", "@islogical", "Save the log10 power spectrum (else an empty file)."),
    Param("background", "@isnumeric", "Background subtracted before the FFT.", convert="same"),
    Param("interpMethod", "@ischar", "Interpolation of the resampling ('linear', 'cubic', ...)."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    task="B04",
)
def fft_spectrum_file(
    frame_fullpath: str,
    spectrum_fullname: str,
    *,
    xy_pixel_size: float = 0.108,
    dz: float = 0.1,
    zarr_file: bool = False,
    out_pixel_size: float | None = None,
    out_size: Sequence[float] = (1001, 1001, 1001),
    save_3d_stack: bool = False,
    background: float = 0,
    interp_method: str = "linear",
    num_threads: int | None = None,
) -> None:
    """The power spectrum of one file (``XR_fftSpectrumComputingFrame``): the MIPs and central
    slices of ``|fftn|`` of the image resampled to ``out_pixel_size`` and cropped or padded to
    ``out_size``, normalized with gamma 0.5, in ``MIPs/`` and ``central_slices/`` next to
    ``spectrum_fullname``, and the result (the log10 spectrum with ``save_3d_stack``, else an
    empty file). Nothing is done if the result exists."""
    from petakit5d.core import dtypes
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.core.threads import resolve_threads
    from petakit5d.io import get_image_size, write_tiff
    from petakit5d.jobs.files import exist
    from petakit5d.jobs.publish import publish_result
    from petakit5d.kernels.ops import convert
    from petakit5d.matlabcompat.fft import fftn
    from petakit5d.tools.mip import axis_max

    uuid = get_uuid()
    frame_fullpath, spectrum_fullname = os.fspath(frame_fullpath), os.fspath(spectrum_fullname)
    if exist(spectrum_fullname):
        return
    nt = resolve_threads(num_threads)
    xyz, px = _pixel_sizes(xy_pixel_size, dz, out_pixel_size)
    n = _yxz(out_size)

    sz = get_image_size(frame_fullpath)[:3]
    box = spectrum_region(sz, xy_pixel_size, dz, out_pixel_size, out_size)
    im = _read_region(frame_fullpath, zarr_file, box)
    if im.ndim != 3:
        # (MATLAB's size of the 2-D image has two values)
        raise FFTSpectrumError(_DIMS_MSG, _DIMS)
    if background != 0:
        im = dtypes.minus(convert(im, "single", num_threads=nt), float(background))
        im = im * (im > 0)
    im = convert(im, "single", num_threads=nt)
    im = _resized(im, xyz, px, n, interp_method, nt)

    # |FFT| without fftshift: the MIPs and central slices are taken from it, then normalized
    # and shifted (the same values as from the normalized, shifted spectrum)
    a = _abs(fftn(im, threads=nt), nt)
    del im
    mips = {axis: axis_max(a, axis, num_threads=nt) for axis in (-2, -1, -3)}
    # (the maximum of the z MIP is the spectrum's)
    im_max = _max_all(mips[-3])

    fft_path, fsname, _ = fileparts(spectrum_fullname)
    tmp = os.path.join(fft_path, f"{fsname}_{uuid}.tif")
    if save_3d_stack:
        # log10(single(fftshift(im) ./ imMax)), shifted after the element-wise steps
        write_tiff(np.fft.fftshift(_log10_single(a, nt, im_max)), tmp, num_threads=nt)
    else:
        with open(tmp, "wb"):
            pass

    mip_path = os.path.join(fft_path, "MIPs")
    for axis, tag in ((-2, "y"), (-1, "x"), (-3, "z")):
        m = mips[axis]
        write_tiff(
            _image(_normalize_gamma(m, im_max), axis),
            os.path.join(mip_path, f"{fsname}_MIP_{tag}.tif"),
            num_threads=nt,
        )
    cs_path = os.path.join(fft_path, "central_slices")
    os.makedirs(cs_path, exist_ok=True)
    for sl, axis, tag in ((a[0:1], -3, "xy"), (a[:, :, 0:1], -1, "yz"), (a[:, 0:1, :], -2, "xz")):
        write_tiff(
            _image(_normalize_gamma(sl, im_max), axis),
            os.path.join(cs_path, f"{fsname}_center_{tag}.tif"),
            num_threads=nt,
        )
    publish_result(tmp, spectrum_fullname)


# ---------------------------------------------------------------------------------------------
# XR_fftSpectrumComputingWrapper
# ---------------------------------------------------------------------------------------------


@register(
    "XR_fftSpectrumComputingWrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param(
        "overwrite",
        "@(x) (numel(x) == 1) && islogical(x)",
        "Overwrite existing results (not used: existing results are kept).",
    ),
    Param("xyPixelSize", "@isnumeric", "Pixel size in um."),
    Param("dz", "@isnumeric", "Scan interval in um."),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param(
        "outPixelSize",
        "@(x) isnumeric(x) || isempty(x)",
        "Output voxel size (um, every axis; None: the smaller of xyPixelSize and dz).",
    ),
    Param(
        "outSize",
        "@isnumeric",
        "Output size (z, y, x): odd numbers to crop or pad to.",
        convert="size_yxz",
    ),
    Param(
        "channelPatterns",
        "@(x) ischar(x) || iscell(x)",
        "Channel identifiers of the files (empty: every file).",
    ),
    Param("save3DStack", "@islogical", "Save the log10 power spectra (in log10 scale)."),
    Param("background", "@isnumeric", "Background subtracted before the FFT.", convert="same"),
    Param(
        "interpMethod",
        "@ischar",
        "Interpolation of the resampling: 'linear', 'cubic', 'lanczos2' or 'lanczos3'.",
    ),
    Param("parseCluster", "@islogical", "Run the files as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("cpusPerTask", "@isscalar", "CPUs of a task's job."),
    Param("debug", "@islogical", "Debug mode (not used)."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    task="B04",
)
def fft_spectrum_dataset(
    data_paths: str | Sequence[str],
    *,
    result_dir_name: str = "FFT",
    overwrite: bool = False,
    xy_pixel_size: float = 0.108,
    dz: float = 0.1,
    zarr_file: bool = False,
    out_pixel_size: float | None = None,
    out_size: Sequence[float] = (1001, 1001, 1001),
    channel_patterns: str | Sequence[str] = (),
    save_3d_stack: bool = False,
    background: float = 0,
    interp_method: str = "linear",
    parse_cluster: bool = True,
    master_compute: bool = True,
    cpus_per_task: float = 3,
    debug: bool = False,
    mcc_mode: bool = False,
    config_file: str = "",
) -> list[bool]:
    """The power spectra of every tiff (or zarr) file of the datasets that matches a channel
    pattern (``XR_fftSpectrumComputingWrapper``), one ``fft_spectrum_file`` task per file
    through the job framework, into ``<data path>/<result_dir_name>/`` (parameters.mat of the
    call there too). Returns the tasks' done flags."""
    values = {k: v for k, v in locals().items() if k in _DATASET_ARGS}
    del overwrite, debug
    from petakit5d.core.paths import matlab_dir
    from petakit5d.core.text import sprintf, through_mat2str, through_text
    from petakit5d.io import get_image_size, write_mat
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.tasks import Task
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain
    from petakit5d.pipeline.filenames import matches_channel_pattern

    patterns = [channel_patterns] if isinstance(channel_patterns, str) else list(channel_patterns)
    if out_pixel_size is None or np.asarray(out_pixel_size).size == 0:
        out_pixel_size = min(float(xy_pixel_size), float(dz))
        print(sprintf("Out pixel size: %d um\n", out_pixel_size), end="", flush=True)

    paths = [data_paths] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    if os.name == "nt":
        paths = [p.replace("\\", "/") for p in paths]
    paths = [p if p.endswith("/") else p + "/" for p in paths]

    pr = {k: to_plain(v) for k, v in matlab_results(fft_spectrum_dataset, values).items()}
    inputs: list[str] = []
    results: list[str] = []
    for d_path in paths:
        fft_path = f"{d_path}{result_dir_name}/"
        os.makedirs(fft_path, exist_ok=True)
        group_writable(fft_path)
        write_mat(fft_path + "/parameters.mat", {"pr": pr})
        names = [e.name for e in matlab_dir(d_path + ("*zarr" if zarr_file else "*tif"))]
        full = [d_path + x for x in names]
        out = [fft_path + x for x in names]
        if patterns:
            keep = [False] * len(full)
            for p in patterns:
                keep = [a or b for a, b in zip(keep, matches_channel_pattern(full, p), strict=True)]
            full = [f for f, k in zip(full, keep, strict=True) if k]
            out = [f for f, k in zip(out, keep, strict=True) if k]
        inputs += full
        results += out

    # the values as the function strings give them
    kw = {
        "xy_pixel_size": float(through_text(float(xy_pixel_size), "%0.20d")),
        "dz": float(through_text(float(dz), "%0.20d")),
        "zarr_file": bool(zarr_file),
        "out_pixel_size": float(through_text(float(_floats(out_pixel_size)[0]), "%0.20d")),
        "out_size": tuple(float(v) for v in _zyx(through_mat2str(_yxz(out_size)))),
        "save_3d_stack": bool(save_3d_stack),
        "background": float(through_text(float(background), "%d")),
        "interp_method": interp_method,
    }
    tasks = [Task.of(fft_spectrum_file, inputs[i], results[i], **kw) for i in range(len(inputs))]
    if not inputs:
        raise FFTSpectrumError(
            "Index exceeds the number of array elements. Index must not exceed 0.",
            "MATLAB:badsubscript",
        )
    sz = np.asarray(get_image_size(inputs[0])[:3], dtype=np.float64)
    mem_allocate = float(np.prod(np.maximum(sz, _floats(out_size)))) * 4 / 1024**3 * 15
    common = {
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "cpus_per_task": cpus_per_task,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = framework.run_tasks(inputs, results, tasks, mem_allocate=mem_allocate, **common)
    if not all(done):
        done = framework.run_tasks(inputs, results, tasks, mem_allocate=mem_allocate * 2, **common)
    return list(done)


def _zyx(v: Any) -> np.ndarray:
    """MATLAB's [y x z] (or one number) as (z, y, x)"""
    a = _floats(v)
    return a[[2, 0, 1]] if a.size == 3 else a


_DATASET_ARGS = (
    "data_paths",
    "result_dir_name",
    "overwrite",
    "xy_pixel_size",
    "dz",
    "zarr_file",
    "out_pixel_size",
    "out_size",
    "channel_patterns",
    "save_3d_stack",
    "background",
    "interp_method",
    "parse_cluster",
    "master_compute",
    "cpus_per_task",
    "debug",
    "mcc_mode",
    "config_file",
)
