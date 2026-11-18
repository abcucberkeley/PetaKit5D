"""ND2 and CZI files to tiff or zarr (``XR_bioformats_to_tiff_or_zarr_wrapper``), without
Bio-Formats.

MATLAB's wrapper reads its files through Bio-Formats' ``bfopen`` (GPL Java), which the
BSD-3-Clause port neither copies nor depends on. The port reads the wrapper's two default
formats with permissive readers, the optional ``formats`` extra (``pip install
petakit5d[formats]``):

- ``.nd2``: ``nd2`` (BSD-3-Clause; with dask, ome-types, pydantic and xsdata: BSD-3 and MIT);
- ``.czi``: ``czifile`` (BSD-3-Clause) with ``imagecodecs`` (BSD-3-Clause; its bundled codec
  libraries under permissive licences).

Other formats raise ``BioformatsError`` naming the supported ones. Each file is one task of the
job framework (in this process by default, as MATLAB converts file after file).

What is converted, as MATLAB does it: each series (an ND2 position ``P``, a CZI scene) is a
stack of planes in Bio-Formats' usual order, channel fastest, then z, then t; the planes are
dealt to ``n_channels`` channels in turn (MATLAB's ``reshape(vol, ny, nx, nChannels, [])``,
whatever channels the file has), and each channel's planes, in order, are its z slices:
``<name>_ch<c>`` (``<name>_fov<series, 4 digits>_ch<c>`` for files of several series), as tiff or
as zarr (chunks ``block_size`` capped at the volume, the ``/`` separator above 10,000 chunks),
written under a temporary name and renamed; existing results are kept unless ``over_write``.
Differences from bfopen: planes come from the readers' axis labels (``T``, ``Z``, ``C``), so a
file whose Bio-Formats dimension order is not XYCZT, an RGB file, or a file with other axes of
more than one element (CZI mosaics, phases, rotations, illuminations) raises instead of being
read in another order; a CZI mosaic is not split into tiles.
"""

from __future__ import annotations

import math
import os
from collections.abc import Iterator, Sequence

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "SUPPORTED_FORMATS",
    "BioformatsError",
    "bioformats_to_tiff_or_zarr_dataset",
    "bioformats_to_tiff_or_zarr_file",
    "read_series",
]

SUPPORTED_FORMATS = (".nd2", ".czi")
_PLANE_AXES = ("T", "Z", "C")


class BioformatsError(RuntimeError):
    """An error of the conversion; ``identifier`` is MATLAB's error id where MATLAB has one."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def _unsupported(what: str) -> BioformatsError:
    return BioformatsError(
        f"{what}: only {', '.join(SUPPORTED_FORMATS)} files are converted (Bio-Formats is not "
        "part of the port; these are read by nd2 and czifile, the formats extra)",
        "petakit5d:bioformats:unsupportedFormat",
    )


def _import(name: str):
    import importlib

    try:
        return importlib.import_module(name)
    except ImportError as err:
        raise ImportError(
            f"{name} is needed to read this file: pip install petakit5d[formats]"
        ) from err


def _planes(a: np.ndarray, dims: Sequence[str], what: str) -> np.ndarray:
    """the (plane, y, x) stack of one series: other axes of one element dropped, planes in
    Bio-Formats' XYCZT order (c fastest, then z, then t)"""
    dims = [str(d) for d in dims]
    if len(dims) != a.ndim or "Y" not in dims or "X" not in dims:
        raise BioformatsError(f"{what}: unexpected axes {''.join(dims)}")
    extra = [
        d for d, n in zip(dims, a.shape, strict=True) if d not in (*_PLANE_AXES, "Y", "X") and n > 1
    ]
    if extra:
        raise BioformatsError(
            f"{what}: axes {''.join(extra)} with more than one element (RGB samples, mosaic "
            "tiles, ...) are not converted",
            "petakit5d:bioformats:unsupportedAxes",
        )
    keep = [i for i, d in enumerate(dims) if d in (*_PLANE_AXES, "Y", "X")]
    drop = tuple(i for i in range(a.ndim) if i not in keep)
    a = a.reshape([n for i, n in enumerate(a.shape) if i not in drop]) if drop else a
    kd = [dims[i] for i in keep]
    for d in _PLANE_AXES:
        if d not in kd:
            a = a[np.newaxis]
            kd = [d, *kd]
    a = np.transpose(a, [kd.index(d) for d in (*_PLANE_AXES, "Y", "X")])
    return np.ascontiguousarray(a).reshape(-1, a.shape[-2], a.shape[-1])


def _nd2_series(path: str) -> Iterator[np.ndarray]:
    nd2 = _import("nd2")
    with nd2.ND2File(path) as f:
        if f.is_rgb or f.components_per_channel > 1:
            raise BioformatsError(
                f"{path}: RGB files are not converted", "petakit5d:bioformats:unsupportedAxes"
            )
        a = np.asarray(f.asarray())
        dims = list(f.sizes)
    if "P" in dims:
        p = dims.index("P")
        rest = dims[:p] + dims[p + 1 :]
        for k in range(a.shape[p]):
            yield _planes(np.take(a, k, axis=p), rest, path)
    else:
        yield _planes(a, dims, path)


def _czi_series(path: str) -> Iterator[np.ndarray]:
    czifile = _import("czifile")
    with czifile.CziFile(path, squeeze=False) as f:
        for key in list(f.scenes.keys()):
            img = f.scenes[key]
            yield _planes(np.asarray(img.asarray()), list(img.dims), path)


def read_series(path: str) -> Iterator[np.ndarray]:
    """The series of an ND2 or CZI file, each as its ``(plane, y, x)`` stack in Bio-Formats'
    XYCZT plane order."""
    ext = os.path.splitext(os.fspath(path))[1].lower()
    if ext == ".nd2":
        return _nd2_series(os.fspath(path))
    if ext == ".czi":
        return _czi_series(os.fspath(path))
    raise _unsupported(os.fspath(path))


def _channels(planes: np.ndarray, n_channels: int) -> list[np.ndarray]:
    """MATLAB's ``reshape(vol, ny, nx, nChannels, [])`` and ``squeeze(vol(:, :, c, :))``"""
    n = planes.shape[0]
    if n_channels < 1 or n % n_channels:
        raise BioformatsError(
            f"Product of known dimensions, {planes.shape[1] * planes.shape[2] * n_channels}, not "
            f"divisible into total number of elements, {planes.size}.",
            "MATLAB:getReshapeDims:notDivisible",
        )
    out = []
    for c in range(n_channels):
        v = planes[c::n_channels]
        out.append(np.ascontiguousarray(v[0] if v.shape[0] == 1 else v))
    return out


def bioformats_to_tiff_or_zarr_file(
    filename: str,
    result_path: str,
    fsn: str,
    flag_fullname: str,
    *,
    n_channels: int = 2,
    save_zarr: bool = False,
    block_size: Sequence[int] = (256, 256, 256),
    over_write: bool = False,
    uuid: str = "",
) -> None:
    """One file of ``bioformats_to_tiff_or_zarr_dataset``: every series and channel written as
    tiff or zarr into ``result_path``, then the task's flag."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.io import write_tiff, write_zarr

    uuid = uuid or get_uuid()
    ext = ".zarr" if save_zarr else ".tif"
    series = list(read_series(filename))
    m_num = len(series)
    for m, planes in enumerate(series, start=1):
        for c, vc in enumerate(_channels(planes, int(n_channels))):
            if m_num == 1:
                fnout = f"{result_path}/{fsn}_ch{c}{ext}"
                tmpout = f"{result_path}/{fsn}_ch{c}_{uuid}{ext}"
            else:
                fnout = f"{result_path}/{fsn}_fov{m:04d}_ch{c}{ext}"
                tmpout = f"{result_path}/{fsn}_fov{m:04d}_ch{c}_{uuid}{ext}"
            if not over_write and os.path.lexists(fnout):
                continue
            if save_zarr:
                b = [int(v) for v in block_size]
                if vc.ndim != len(b):
                    raise BioformatsError(
                        "Arrays have incompatible sizes for this operation.",
                        "MATLAB:sizeDimensionsMustMatch",
                    )
                block = [min(s, k) for s, k in zip(vc.shape, b, strict=True)]
                nchunks = math.prod(s / k for s, k in zip(vc.shape, block, strict=True))
                sep = "/" if nchunks > 10000 else "."
                write_zarr(vc, tmpout, block_size=block, dim_separator=sep)
            else:
                write_tiff(vc, tmpout)
            if os.path.isdir(fnout) and over_write:
                import shutil

                shutil.rmtree(fnout)
            os.replace(tmpout, fnout)
    with open(flag_fullname, "w"):
        pass


@register(
    "XR_bioformats_to_tiff_or_zarr_wrapper",
    Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Dataset folders (one or several)."),
    Param("resultDirName", "@ischar", "Result folder under each data path."),
    Param("channelPatterns", "@iscell", "File name patterns of the files to convert."),
    Param("nChannels", "@isnumeric", "Channels the planes of each series are dealt to."),
    Param("dataFormats", "@iscell", "File formats to convert (.nd2 and .czi are supported)."),
    Param("saveZarr", "@islogical", "Write zarr (else tiff)."),
    Param("blockSize", "@isnumeric", "Zarr chunk size (z, y, x).", convert="size_yxz"),
    Param("overWrite", "@islogical", "Replace existing results.", py="over_write"),
    Param("uuid", "@ischar", "Id of the temporary results (empty: a new one)."),
    Param("parseCluster", None, "Run the files as Slurm jobs.", python_only=True),
    Param("masterCompute", None, "This process converts files too.", python_only=True),
    Param("cpusPerTask", None, "CPUs of a file's job.", python_only=True),
    Param("jobLogDir", None, "Folder of the job logs.", python_only=True),
    Param("configFile", None, "Cluster configuration file.", python_only=True),
    task="B09",
)
def bioformats_to_tiff_or_zarr_dataset(
    data_paths: str | Sequence[str],
    *,
    result_dir_name: str = "tiffs",
    channel_patterns: Sequence[str] = (".nd2", ".czi"),
    n_channels: float = 2,
    data_formats: Sequence[str] = (".nd2", ".czi"),
    save_zarr: bool = False,
    block_size: Sequence[int] = (256, 256, 256),
    over_write: bool = False,
    uuid: str = "",
    parse_cluster: bool = False,
    master_compute: bool = True,
    cpus_per_task: float = 1,
    job_log_dir: str = "../job_logs",
    config_file: str = "",
) -> None:
    """Convert the ND2 and CZI files of datasets (those that match a channel pattern and a data
    format) to one tiff or zarr per series and channel in ``<data path>/<result_dir_name>/``
    (``XR_bioformats_to_tiff_or_zarr_wrapper``, without Bio-Formats: the readers of the
    ``formats`` extra), one task per file through the job framework (in this process unless
    ``parse_cluster``)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import simplify_path
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs import framework
    from petakit5d.jobs.publish import remove_path
    from petakit5d.jobs.tasks import Task
    from petakit5d.pipeline.filenames import parse_image_filenames

    data_formats = [str(f) for f in data_formats]
    bad = [f for f in data_formats if f.lower() not in SUPPORTED_FORMATS]
    if bad:
        raise _unsupported("data formats " + ", ".join(bad))
    paths = (
        [os.fspath(data_paths)] if isinstance(data_paths, (str, os.PathLike)) else list(data_paths)
    )
    result_paths = []
    for d, p in enumerate(paths):
        p = os.fspath(p)
        if not p.endswith(("/", os.sep)):
            p = p + os.sep
        paths[d] = simplify_path(p) + "/"
        rp = f"{paths[d]}/{result_dir_name}/"
        os.makedirs(rp, exist_ok=True)
        group_writable(rp)
        result_paths.append(rp)
    if not uuid:
        uuid = get_uuid()

    info = parse_image_filenames(paths, False, list(channel_patterns), data_formats)
    fnames, fsns, fd_inds = info.fnames, info.fsns, info.fd_inds
    tasks, inputs, outputs = [], [], []
    flag_paths = sorted({f"{result_paths[int(d)]}bioformats_flag_{uuid}" for d in fd_inds})
    for fp in flag_paths:
        os.makedirs(fp, exist_ok=True)
    for f, fname in enumerate(fnames):
        d = int(fd_inds[f])
        fn = f"{paths[d]}/{fname}"
        flag = f"{result_paths[d]}bioformats_flag_{uuid}/{fsns[f]}.flag"
        inputs.append(fn)
        outputs.append(flag)
        tasks.append(
            Task.of(
                bioformats_to_tiff_or_zarr_file,
                fn,
                result_paths[d],
                fsns[f],
                flag,
                n_channels=int(n_channels),
                save_zarr=bool(save_zarr),
                block_size=[int(v) for v in block_size],
                over_write=bool(over_write),
                uuid=uuid,
            )
        )
    try:
        if tasks:
            mem = max(os.path.getsize(p) for p in inputs) / 2**30 * 3
            common = {
                "parse_cluster": parse_cluster,
                "master_compute": master_compute,
                "cpus_per_task": cpus_per_task,
                "job_log_dir": job_log_dir,
                "config_file": config_file,
                "mem_allocate": mem,
            }
            done = framework.run_tasks(inputs, outputs, tasks, **common)
            if not all(done):
                common["mem_allocate"] = mem * 2
                done = framework.run_tasks(inputs, outputs, tasks, **common)
            if not all(done):
                raise BioformatsError(f"{sum(not d for d in done)} / {len(done)} files failed!")
    finally:
        for fp in flag_paths:
            remove_path(fp)
