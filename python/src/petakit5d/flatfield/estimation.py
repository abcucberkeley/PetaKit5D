"""Camera background estimation from dark stacks: ``XR_flatfield_background_estimation`` and
``XR_flatfield_background_wrapper``.

- ``estimate_background(stack, ...)``: the estimate of a ``(z, y, x)`` single stack, as
  ``XR_flatfield_background_estimation`` computes it (no files);
- ``flatfield_background_estimation(frame_fullnames, result_fullname, ...)``: the per-file task,
  the frames read and concatenated along z, the estimate saved as a v7.3 MAT-file (``imp``,
  ``mu``, ``sigma``, ``med``, ``numSlice``, ``pmat``) and optionally the z MIP of the frames;
- ``flatfield_background_dataset(data_paths, ...)``: the wrapper, one task per file (or per
  camera and channel) through the job framework, then each dataset's camera and channel
  background as ``<prefix>Scan_Cam<c>_ch<n>_CAM1_background.tif``: the mean of the tasks' ``mu``
  weighted by their ``numSlice``.

The estimate, in MATLAB's order and classes (single data; ``petakit5d.matlabcompat`` for the
means, standard deviations, medians and percentiles)::

    sigma_z = std(im, 0, [1 2]);  im = im(:, :, sigma_z < min(sigmaThresh * 5, median(sigma_z)))
    mu = mean(im, 3); sigma = std(im, 0, 3); sigma_z = std(im - mu, 0, [1 2])
    repeat (at most 100 times): drop the slices with more than probThresh of their pixels above
        mu + sigmaFactor * sigma that are noisy (sigma_z > min(sigmaThresh, median of the first
        sigma_z)) or have a range of im - mu above 20 sigma_z; stop when under 0.1% would go
    med = median(im, 3); imp = prctile(im, 0 : 5 : 100, 3) for more than 200 slices, else []

Kept from MATLAB: a result file name without a folder is in the current folder, its MIP in the
current folder's ``MIPs`` (the caller makes that folder, as the wrapper does); the wrapper's
``cam_ch`` and ``came_ch_iter`` modes index their task lists by channel and camera (not by
dataset: the last dataset's task of a channel and camera is the one run, written under the
last dataset's result folder, and ``came_ch_iter`` takes the files of every dataset with the
current dataset's folder); a background averages the tasks of its camera and channel of every
dataset (MATLAB's ``contains`` on the task output paths); ``cam_ch_iter_stack`` and unknown modes
raise MATLAB's errors; a retry of failed tasks (twice the memory) is not counted, so the call
raises ``Some output files are missing!`` after it even when it succeeds (known issue).

Arrays are the port's: stacks ``(z, y, x)``, images ``(y, x)``. Sums over the 4,096-element
blocks of planes of 89,000 pixels or more are MATLAB's at one thread (``matlabcompat.reduce``);
MATLAB's own multithreaded sums there can differ in the last bit, which can move a slice across
a threshold only when it sits on one. Speed: the planes' standard deviations go plane by plane
through the ordered sums of ``kernels.ops``, the median is taken by selection (``stats.median``
where a column has NaN or the median is a zero), and the passes over slices run in threads;
no result depends on the thread count.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "MatlabError",
    "estimate_background",
    "flatfield_background_dataset",
    "flatfield_background_estimation",
]

_F32 = np.dtype(np.float32)
PMAT = tuple(range(0, 101, 5))  # 0 : 5 : 100
MAX_ITERATIONS = 100
PERCENTILE_SLICES = 200  # prctile for more slices than this


class MatlabError(RuntimeError):
    """An error MATLAB raises, with its identifier ('' for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# the estimate
# ---------------------------------------------------------------------------------------------


def _threads(num_threads: int | None) -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(num_threads)


def _parts(n: int, nt: int) -> list[range]:
    """n slices (or pixels) in up to 4 * nt runs"""
    k = max(1, min(n, 4 * nt))
    return [range(n * i // k, n * (i + 1) // k) for i in range(k)]


def _pmap(fn, parts: list[range], nt: int) -> list:
    if nt <= 1 or len(parts) <= 1:
        return [fn(p) for p in parts]
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(nt) as ex:
        return list(ex.map(fn, parts))


def _single(im: np.ndarray, nt: int) -> np.ndarray:
    """``single(im)``: 8- and 16-bit integers are exact in single (a copy, slice runs in
    threads), other classes by ``core.dtypes``"""
    from petakit5d.core import dtypes

    if im.dtype == _F32:
        return im
    if im.dtype.kind not in "iu" or im.dtype.itemsize > 2 or im.ndim != 3:
        return np.asarray(dtypes.to_class(im, "single"))
    out = np.empty(im.shape, _F32)

    def run(r: range) -> None:
        out[r.start : r.stop] = im[r.start : r.stop]

    _pmap(run, _parts(im.shape[0], nt), nt)
    return out


def _std_planes(d: np.ndarray, nt: int) -> np.ndarray:
    """``squeeze(std(d, 0, [1, 2]))``: one value per z (the vector rule on each ``d(:, :, k)``,
    as ``std`` of the plane)"""
    from petakit5d.matlabcompat import stats

    out = np.empty(d.shape[0], _F32)
    for k in range(d.shape[0]):
        out[k] = np.asarray(stats.std(d[k], 0, None, num_threads=nt)).reshape(-1)[0]
    return out


def _min_single(a: Any, s: np.float32) -> np.float32:
    """``min(a, s)`` of a double scalar and a single one: single (NaN ignored)"""
    from petakit5d.matlabcompat import reduce

    return np.float32(np.asarray(reduce.minimum(np.float64(a), np.float32(s))).reshape(-1)[0])


def _moments(im: np.ndarray, nt: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``mu = mean(im, 3); sigma = std(im, 0, 3)`` and ``im - mu``"""
    from petakit5d.matlabcompat import reduce, stats

    mu = np.asarray(reduce.mean(im, axis=0, num_threads=nt), dtype=_F32)
    sigma = np.asarray(stats.std(im, 0, axis=0, num_threads=nt), dtype=_F32)
    d = np.empty_like(im)

    def sub(r: range) -> None:
        np.subtract(im[r.start : r.stop], mu, out=d[r.start : r.stop])

    _pmap(sub, _parts(im.shape[0], nt), nt)
    return mu, sigma, d


def _plane_stats(im, d, cut, nt) -> tuple[np.ndarray, np.ndarray]:
    """per z: the pixels of ``im`` above ``cut``, and ``max(d, [], [1 2]) - min(d, [], [1 2])``
    (NaN ignored unless a plane is all NaN; a signed zero can't change the comparison it
    feeds)"""
    nz, npix = im.shape[0], im.shape[1] * im.shape[2]
    high = np.zeros(nz, np.int64)
    rng = np.zeros(nz, _F32)

    def run(r: range) -> None:
        a, b = r.start, r.stop
        if b <= a or npix == 0:
            return
        with np.errstate(invalid="ignore"):
            high[a:b] = np.count_nonzero((im[a:b] > cut).reshape(b - a, npix), axis=1)
            w = d[a:b].reshape(b - a, npix)
            rng[a:b] = np.subtract(np.fmax.reduce(w, axis=1), np.fmin.reduce(w, axis=1))

    _pmap(run, _parts(nz, nt), nt)
    return high, rng


def _median_z(im: np.ndarray, nt: int) -> np.ndarray:
    """``median(im, 3)``: the middle values by selection, pixel by pixel; where a column holds
    NaN or the median is a zero (whose sign MATLAB's stable sort decides), ``stats.median``"""
    from petakit5d.matlabcompat import stats

    n, ny, nx = im.shape
    npix = ny * nx
    if n == 0 or npix == 0:
        return np.asarray(stats.median(im, axis=0), dtype=_F32)
    flat = im.reshape(n, npix)
    out = np.empty(npix, _F32)
    kth = [n // 2] if n % 2 else [n // 2 - 1, n // 2]
    special = np.zeros(npix, bool)

    def run(r: range) -> None:
        a, b = r.start, r.stop
        w = np.ascontiguousarray(flat[:, a:b].T)
        nan = np.isnan(w).any(axis=1)
        w = np.partition(w, kth, axis=1)
        if n % 2:
            v = w[:, n // 2]
        else:
            v = np.asarray(stats._meanof(w[:, n // 2 - 1], w[:, n // 2], _F32))
        out[a:b] = v
        special[a:b] = nan | (v == 0)

    _pmap(run, _parts(npix, nt), nt)
    idx = np.flatnonzero(special)
    if idx.size:
        out[idx] = np.asarray(stats.median(flat[:, idx], axis=0), dtype=_F32)
    return out.reshape(ny, nx)


def estimate_background(
    stack: np.ndarray,
    *,
    sigma_thresh: Any = 8,
    sigma_factor: Any = 4,
    prob_thresh: Any = 0.001,
    num_threads: int | None = None,
) -> dict[str, Any]:
    """The background estimate of a dark ``(z, y, x)`` stack (made single), as
    ``XR_flatfield_background_estimation`` computes it: a dict of ``imp`` (``(21, y, x)``
    percentiles, or None for 200 slices or fewer), ``mu``, ``sigma``, ``med`` (``(y, x)``
    single images), ``numSlice`` (the slices kept, a float) and ``pmat`` (0, 5, ..., 100).
    ``num_threads`` (None: the CPUs available) never changes a result."""
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat import stats

    nt = _threads(num_threads)
    im = np.asarray(stack)
    if im.ndim == 2:
        im = im[None]
    if im.dtype != _F32:
        im = _single(im, nt)
    im = np.ascontiguousarray(im)
    st = float(np.asarray(sigma_thresh, dtype=np.float64).reshape(-1)[0])
    sf = np.asarray(sigma_factor, dtype=np.float64)
    pt = float(np.asarray(prob_thresh, dtype=np.float64).reshape(-1)[0])

    sigma_z = _std_planes(im, nt)
    im = im[sigma_z < _min_single(st * 5, stats.median(sigma_z))]

    mu, sigma, d = _moments(im, nt)
    sigma_z = _std_planes(d, nt)
    lim = _min_single(st, stats.median(sigma_z))
    for _ in range(MAX_ITERATIONS):
        npix = im.shape[1] * im.shape[2]
        cut = np.asarray(dtypes.plus(mu, dtypes.times(sf, sigma)), dtype=_F32)
        high, rng = _plane_stats(im, d, cut, nt)
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = high / np.float64(npix) if npix else np.full(im.shape[0], np.nan)
            bright = np.asarray(dtypes.times(sigma_z, 20), dtype=_F32)
            remove = (frac > pt) & ((sigma_z > lim) | (rng > bright))
        n = remove.size
        if n and np.count_nonzero(remove) / n < 0.001:
            break
        im = im[~remove]
        mu, sigma, d = _moments(im, nt)
        sigma_z = _std_planes(d, nt)
    del d

    med = _median_z(im, nt)
    num_slice = im.shape[0]
    imp = None
    if num_slice > PERCENTILE_SLICES:
        imp = np.asarray(stats.prctile(im, np.asarray(PMAT, dtype=np.float64), axis=0))
    return {
        "imp": imp,
        "mu": mu,
        "sigma": sigma,
        "med": med,
        "numSlice": float(num_slice),
        "pmat": np.asarray(PMAT, dtype=np.float64),
    }


def _read_frames(frames: Sequence[str], zarr_file: bool) -> np.ndarray:
    """``im = single(cat(3, im{:}))`` of the frames read (z, y, x)"""
    from petakit5d.io import read_tiff, read_zarr

    vols = []
    for fn in frames:
        v = np.asarray(read_zarr(fn) if zarr_file else read_tiff(fn))
        vols.append(v[None] if v.ndim == 2 else v)
    im = vols[0] if len(vols) == 1 else np.concatenate(vols, axis=0)
    return np.ascontiguousarray(_single(im, _threads(None)))


def _publish_tiff(img: np.ndarray, tmp: str, final: str) -> None:
    from petakit5d.io import write_tiff

    write_tiff(img, tmp)
    os.replace(tmp, final)


def _frames_list(frame_fullnames: Any) -> list[str]:
    if isinstance(frame_fullnames, str | os.PathLike):
        # MATLAB indexes a char array with braces
        raise MatlabError(
            "Brace indexing is not supported for variables of this type.",
            "MATLAB:cellRefFromNonCell",
        )
    return [os.fspath(f) for f in frame_fullnames]


@register(
    "XR_flatfield_background_estimation",
    Param("frameFullnames", None, "Files of the dark stack (concatenated along z)."),
    Param("resultFullname", None, "Result MAT-file."),
    Param("objectiveScan", "@islogical", "Not used."),
    Param("Overwrite", "@islogical", "Replace an existing result."),
    Param("zarrFile", "@islogical", "The frames are zarr."),
    Param("MovieSelector", "@ischar", "Not used."),
    Param(
        "sigmaThresh", "@isnumeric", "Camera noise threshold of a background slice.", convert="same"
    ),
    Param(
        "sigmaFactor",
        "@isnumeric",
        "Factor of the pixel sigma above which a pixel is bright.",
        convert="same",
    ),
    Param(
        "probThresh",
        "@isnumeric",
        "Fraction of bright pixels above which a slice is dropped.",
        convert="same",
    ),
    Param("save16bit", "@islogical", "Not used."),
    Param("saveMIP", "@islogical", "Save the z MIP of the frames (in the result folder's MIPs)."),
    Param("uuid", "@ischar", "Id of the temporary files (empty: a new one)."),
    task="B07",
)
def flatfield_background_estimation(
    frame_fullnames: Sequence[str],
    result_fullname: str,
    *,
    objective_scan: bool = False,
    overwrite: bool = False,
    zarr_file: bool = False,
    movie_selector: str = "cell",
    sigma_thresh: Any = 8,
    sigma_factor: Any = 4,
    prob_thresh: Any = 0.001,
    save_16bit: bool = True,
    save_mip: bool = False,
    uuid: str = "",
) -> None:
    """The camera background estimate of a dark stack (``XR_flatfield_background_estimation``):
    the frames read and concatenated along z (single), the estimate (``estimate_background``)
    saved as ``result_fullname`` (v7.3: ``imp``, ``mu``, ``sigma``, ``med``, ``numSlice``,
    ``pmat``), and with ``save_mip`` the z MIP of the frames as
    ``<folder>/MIPs/<name>_MIP_z.tif``. An existing result is kept unless ``overwrite``. A name
    without a folder is in the current folder."""
    from petakit5d.core.arrays import to_matlab
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.io import write_mat
    from petakit5d.jobs.files import exist
    from petakit5d.matlabcompat import reduce

    result_fullname = os.fspath(result_fullname)
    if not uuid:
        uuid = get_uuid()
    if exist(result_fullname) and not overwrite:
        print(f"The result {result_fullname} already exists, skip it!")
        return
    frames = _frames_list(frame_fullnames)
    im = _read_frames(frames, bool(zarr_file))

    if save_mip:
        mip, _ = reduce.max(im, axis=0)
        out_path, fsn, _ = fileparts(result_fullname)
        if not out_path:
            out_path = "."
        _publish_tiff(
            np.ascontiguousarray(mip),
            f"{out_path}/MIPs/{fsn}_MIP_z_{uuid}.tif",
            f"{out_path}/MIPs/{fsn}_MIP_z.tif",
        )

    r = estimate_background(
        im, sigma_thresh=sigma_thresh, sigma_factor=sigma_factor, prob_thresh=prob_thresh
    )
    del im
    pathstr, fname, ext = fileparts(result_fullname)
    if not pathstr:
        pathstr = "."
    tmp = f"{pathstr}/{fname}_{uuid}{ext}"
    imp = r["imp"]
    write_mat(
        tmp,
        {
            "imp": None if imp is None else to_matlab(imp),
            "mu": r["mu"],
            "sigma": r["sigma"],
            "med": r["med"],
            "numSlice": r["numSlice"],
            "pmat": r["pmat"],
        },
    )
    os.replace(tmp, result_fullname)


# ---------------------------------------------------------------------------------------------
# the wrapper
# ---------------------------------------------------------------------------------------------

_EXPR_CAM = (
    r"(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?\d+?_?\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)"
    r"_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_"
    r"(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t"
)
_EXPR_NO_CAM = (
    r"(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?\d+?_?\d+?)_ch(?<ch>\d+)_CAM1_stack"
    r"(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_"
    r"(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t"
)


def _table(fnames: Sequence[str], fd_inds: Sequence[int], zarr_file: bool):
    """the wrapper's table of the names that match its expression: (rows, specifyCam); a row is
    (name, prefix, Iter, camera, ch, folder index)"""
    from petakit5d.pipeline.filenames import matlab_regex, regexp_names

    ext = ".zarr" if zarr_file else ".tif"
    if all(matlab_regex(r"_Cam\w_ch").search(f) for f in fnames):
        expr, specify_cam = _EXPR_CAM + ext, True
    elif all(matlab_regex(r"_ch[0-9]_").search(f) for f in fnames):
        expr, specify_cam = _EXPR_NO_CAM + ext, False
    else:
        raise MatlabError(
            "Unrecognized function or variable 'expression'.", "MATLAB:UndefinedFunction"
        )
    rx = matlab_regex(expr, True)
    rows = []
    for f, name in enumerate(fnames):
        m = regexp_names(name, rx)
        if not m:
            continue
        m = m[0]
        cam = m["Cam"] if specify_cam else "N"
        if len(cam) != 1:
            raise MatlabError(
                "Unable to perform assignment because the left and right sides have a different "
                "number of elements.",
                "MATLAB:subsassignnumelmismatch",
            )
        rows.append((name, m["prefix"], int(m["Iter"]), cam, int(m["ch"]), int(fd_inds[f])))
    if not rows:
        raise MatlabError(
            "Unrecognized table variable name 'prefix'.", "MATLAB:table:UnrecognizedVarName"
        )
    return rows, specify_cam


def _cell_order(entries: dict[tuple[int, ...], Any]) -> list[Any]:
    """the non-empty cells of a MATLAB cell array assigned at these 1-based subscripts, in
    MATLAB's linear order (the first subscript fastest; the cell grows along its first
    dimension only)"""
    return [entries[k] for k in sorted(entries, key=lambda k: tuple(reversed(k)))]


def _task_lists(
    mode: str,
    paths: list[str],
    result_dir: str,
    bg_rt: str,
    fnames: list[str],
    fd_inds: list[int],
    gfnames: list[list[str]],
    rows: list,
    chs: list[int],
    cams: list[str],
    iters: list[int],
    zarr_file: bool,
) -> tuple[list[str], list[list[str]], list[str]]:
    """inputFullpaths, input_str_filenames and outputFullpaths of the estimation mode"""
    if mode == "per_image":
        cut = 4 + int(bool(zarr_file))
        ins = [f"{paths[fd_inds[x]]}/{fnames[x]}" for x in range(len(fnames))]
        strs = [[f"{paths[fd_inds[x]]}/{y}" for y in gfnames[x]] for x in range(len(fnames))]
        outs = [
            f"{paths[fd_inds[x]]}/{result_dir}/{fnames[x][: len(fnames[x]) - cut]}.mat"
            for x in range(len(fnames))
        ]
        return ins, strs, outs
    if mode == "cam_ch_iter_stack":
        raise MatlabError("Unrecognized function or variable 'nume'.", "MATLAB:UndefinedFunction")
    if mode not in ("cam_ch", "came_ch_iter"):
        raise MatlabError(
            "Unrecognized function or variable 'inputFullpaths'.", "MATLAB:UndefinedFunction"
        )
    # gfnames of the table's rows (rows keep the order of fnames)
    by_name = {}
    for x, name in enumerate(fnames):
        by_name.setdefault(name, gfnames[x])
    cells: dict[tuple[int, ...], tuple[str, list[str], str]] = {}
    for d, data_path in enumerate(paths):
        for c, ch in enumerate(chs, 1):
            for ncam, cam in enumerate(cams, 1):
                cur = [r for r in rows if r[5] == d and r[4] == ch and r[3] == cam]
                if not cur:
                    continue
                if mode == "cam_ch":
                    name = cur[0][0]
                    cells[(c, ncam)] = (
                        f"{data_path}/{name}",
                        [f"{data_path}/{g}" for r in cur for g in by_name[r[0]]],
                        f"{bg_rt}/{name[: len(name) - 4]}.mat",
                    )
                    continue
                for n, it in enumerate(iters, 1):
                    cur = [r for r in rows if r[4] == ch and r[3] == cam and r[2] == it]
                    if not cur:
                        continue
                    name = cur[0][0]
                    cells[(c, ncam, n)] = (
                        f"{data_path}/{name}",
                        [f"{data_path}/{g}" for r in cur for g in by_name[r[0]]],
                        f"{data_path}/{result_dir}/{name[: len(name) - 4]}.mat",
                    )
    items = _cell_order(cells)
    return [i[0] for i in items], [i[1] for i in items], [i[2] for i in items]


def _save_parameters(result_path: str, values: dict[str, Any]) -> None:
    """parameters.mat (v7.3) and parameters.json of the call (MATLAB's ``pr``)"""
    from petakit5d.core.textfiles import write_json_file
    from petakit5d.io import write_mat
    from petakit5d.params import matlab_results
    from petakit5d.params.values import to_plain

    pr = {k: to_plain(v) for k, v in matlab_results(flatfield_background_dataset, values).items()}
    write_mat(result_path + "/parameters.mat", {"pr": pr})
    write_json_file(pr, result_path + "/parameters.json")


def _background_images(
    paths: list[str],
    result_dir: str,
    prefix: str,
    chs: list[int],
    cams: list[str],
    specify_cam: bool,
    outs: list[str],
    uuid: str,
) -> None:
    """each dataset's background of a camera and a channel: the tasks' ``mu`` weighted by their
    ``numSlice`` (tasks without slices left out), ``sum(z .* w, 3)`` in single"""
    from petakit5d.core import dtypes
    from petakit5d.io import read_mat
    from petakit5d.jobs.files import exist
    from petakit5d.matlabcompat import reduce

    for data_path in paths:
        bg_rt = f"{data_path}/{result_dir}/"
        for ch in chs:
            for cam in cams:
                bg_fname = f"{bg_rt}/{prefix}Scan_Cam{cam}_ch{ch}_CAM1_background.tif"
                if exist(bg_fname):
                    continue
                key = f"Cam{cam}_ch{ch}" if specify_cam else f"ch{ch}"
                cur = [o for o in outs if key in o]
                if not cur:
                    continue
                mus, sizes = [], []
                for o in cur:
                    a = read_mat(o, ["mu", "numSlice"])
                    n = float(np.asarray(a["numSlice"]).reshape(-1)[0])
                    if n != 0:
                        mus.append(np.ascontiguousarray(a["mu"]))
                        sizes.append(n)
                z = np.stack(mus)
                sz = np.asarray(sizes, dtype=np.float64)
                tot = np.asarray(reduce.sum(sz))
                w = (sz / tot).reshape(-1, 1, 1)
                zw = np.asarray(dtypes.times(z, w))
                z_im = np.asarray(reduce.sum(zw, axis=0))
                _publish_tiff(
                    np.ascontiguousarray(z_im),
                    f"{bg_rt}/{prefix}Scan_Cam{cam}_ch{ch}_CAM1_background_{uuid}.tif",
                    bg_fname,
                )


@register(
    "XR_flatfield_background_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param("resultDirStr", "@ischar", "Result folder under each data path."),
    Param("channelPatterns", "@iscell", "Channel identifiers of the dark stacks."),
    Param(
        "estimationMode",
        "@ischar",
        "'per_image' (a task per file), 'cam_ch' or 'came_ch_iter' (per camera and channel).",
    ),
    Param("onlyFirstTP", "@islogical", "Only the first iteration ('came_ch_iter')."),
    Param("zarrFile", "@islogical", "The files are zarr."),
    Param("save16bit", "@islogical", "Not used."),
    Param("saveMIP", "@islogical", "Save the z MIP of each task's stack."),
    Param("parseCluster", "@islogical", "Run the tasks as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process computes tasks too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a task's job."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    task="B07",
)
def flatfield_background_dataset(
    data_paths: str | Sequence[str],
    *,
    result_dir_str: str = "matlab_flat_field_estimation",
    channel_patterns: Sequence[str] = ("CamA_ch0", "CamB_ch0"),
    estimation_mode: str = "per_image",
    only_first_tp: bool = False,
    zarr_file: bool = False,
    save_16bit: bool = True,
    save_mip: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = True,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 24,
    mcc_mode: bool = False,
    config_file: str = "",
    uuid: str = "",
) -> None:
    """Camera backgrounds of dark datasets (``XR_flatfield_background_wrapper``): an estimate
    per file (``flatfield_background_estimation``; or per camera and channel) in
    ``<data path>/<result_dir_str>/<name>.mat`` through the job framework, then per dataset,
    camera and channel ``<prefix>Scan_Cam<c>_ch<n>_CAM1_background.tif``, the estimates' means
    weighted by their slice counts. parameters.mat and .json of the call go into each result
    folder."""
    values = {k: v for k, v in locals().items() if k in _WRAPPER_ARGS}
    from petakit5d.core.arrays import data_type_to_byte_number
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import get_image_dtype, get_image_size
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_dataset_filenames

    paths = [data_paths] if isinstance(data_paths, str | os.PathLike) else list(data_paths)
    paths = [os.fspath(p) for p in paths]
    bg_rt = ""
    for data_path in paths:
        if not os.path.isdir(data_path):
            raise MatlabError(
                f"{data_path} does not exist, please double check the path to make sure it is "
                "correct!"
            )
        bg_rt = data_path + os.sep + result_dir_str
        if not os.path.isdir(bg_rt):
            os.makedirs(bg_rt)
            group_writable(bg_rt)
        _save_parameters(bg_rt, values)
        tmp = bg_rt + os.sep + "tmp"
        if not os.path.isdir(tmp):
            os.makedirs(tmp)
            group_writable(tmp)
        if save_mip:
            mip_path = bg_rt + os.sep + "MIPs"
            if not os.path.isdir(mip_path):
                os.makedirs(mip_path)
                group_writable(mip_path)

    if parse_cluster:
        chk = check_slurm_cluster(paths[0], job_log_dir)
        parse_cluster = bool(chk.parse_cluster)
        if chk.job_log_dir is not None:
            job_log_dir = chk.job_log_dir
    if not uuid:
        uuid = get_uuid()

    fnames, fd_inds, gfnames, *_ = parse_dataset_filenames(
        paths, list(channel_patterns), False, False, False, "", False, 0, bool(zarr_file)
    )
    fd = [int(i) for i in fd_inds]
    rows, specify_cam = _table(fnames, fd, bool(zarr_file))
    prefix = sorted({r[1] for r in rows})[0]
    iters = sorted({r[2] for r in rows})
    chs = sorted({r[4] for r in rows})
    cams = sorted({r[3] for r in rows})
    if only_first_tp:
        iters = iters[:1]

    ins, strs, outs = _task_lists(
        estimation_mode,
        paths,
        result_dir_str,
        bg_rt,
        list(fnames),
        fd,
        gfnames,
        rows,
        chs,
        cams,
        iters,
        bool(zarr_file),
    )
    tasks = [
        Task.of(
            flatfield_background_estimation,
            strs[x],
            outs[x],
            zarr_file=bool(zarr_file),
            save_mip=bool(save_mip),
        )
        for x in range(len(ins))
    ]

    first = paths[fd[0]]
    sizes = [get_image_size(f"{first}/{y}") for y in gfnames[0]]
    nz = sum(int(s[0]) for s in sizes)
    voxels = float(sizes[0][1]) * float(sizes[0][2]) * nz
    byte = data_type_to_byte_number(get_image_dtype(f"{first}/{gfnames[0][0]}"))
    mem_allocate = voxels * byte / 2**30 * 5

    kw = {
        "cpus_per_task": cpus_per_task,
        "parse_cluster": parse_cluster,
        "master_compute": master_compute,
        "job_log_dir": job_log_dir,
        "mcc_mode": mcc_mode,
        "config_file": config_file,
    }
    done = framework.run_tasks(ins, outs, tasks, mem_allocate=mem_allocate, **kw)
    if not all(done):
        # (its result is not kept: MATLAB's call does not assign it)
        framework.run_tasks(ins, outs, tasks, mem_allocate=mem_allocate * 2, **kw)
    if not all(done):
        raise MatlabError("Some output files are missing!")

    _background_images(paths, result_dir_str, prefix, chs, cams, specify_cam, outs, uuid)


_WRAPPER_ARGS = (
    "data_paths",
    "result_dir_str",
    "channel_patterns",
    "estimation_mode",
    "only_first_tp",
    "zarr_file",
    "save_16bit",
    "save_mip",
    "parse_cluster",
    "master_compute",
    "job_log_dir",
    "cpus_per_task",
    "mcc_mode",
    "config_file",
    "uuid",
)
